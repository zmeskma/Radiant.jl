"""
    feed(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,
    is_elastic::Bool,is_subshells::Bool)

Calculate the feed function 𝓕 (normalized probability of scattering from Ei into each
group gf) for each Legendre moment up to order L. Also calculate the energy weighted
feed function 𝓕ₑ for energy-deposition cross section.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `atz::Vector{Float64}` : atomic percent of the element(s) composing the material.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction informations.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` :  number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of the particles imply in the interaction.
- `type::String` : type of interaction (scattering or production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic::Bool` : boolean indicating if the scattering is elastic.
- `is_subshells::Bool` : boolean indicating if the cross-sections are subshells dependant.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function (per unit nuclei density).
- `𝓕ₑ::Vector{Float64}` : energy weighted feed function (per unit nuclei density).

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool)

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)
ΔQ = get_mass_energy_variation(interaction,type,true)

# Outgoing particle energy spectrum
is_dirac, Np, q_type = out_distribution_dispatch(interaction,type)
if is_dirac Np = 1; u = [0]; w = [2] else u,w = quadrature(Np,q_type) end

#----
# Feed function over all groups and under the cutoff energy
#----

# Heavy inelastic S (same scattered particle)
is_heavy_inelastic_S = (is_proton(incoming_particle) || is_alpha(incoming_particle)) && (incoming_particle == scattered_particle)

# Loop over the compound elements
Nz = length(Z)
for i in range(1,Nz)

    # Loop over subshells and outgoing groups
    Nshells,Zi,Ui,Ti,ri,_ = electron_subshells(Z[i],~is_subshells)
    for gf in range(1,Ng), δi in range(1,Nshells)

        # Final energy group
        Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
        Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui[δi],Ec,incoming_particle)
        if isSkip continue end
        ΔEf = Ef⁻ - Ef⁺

        # Integration over the energy group
        𝓕i = zeros(L+1)
        𝓕iₑ = 0

        # For heavy particles compute cache once and reuse for all analytic integrals
        analytic_A = 0.0
        M₁ = 0.0
        cache = nothing
        if is_heavy_inelastic_S
            cache = HeavyInelasticCache(Zi[δi], Ei, incoming_particle)
            analytic_A = integrate_A_over_W2_per_subshell(cache, Ef⁻, Ef⁺) * atz[i]
            M₁ = feed_first_moment_heavy_particle(cache, Ef⁻, Ef⁺) * atz[i]
        end

        # Use quadrature integration for all particles
        for n in range(1,Np)
            # Outgoing particle energy group
            if (is_elastic) Ef = Ei else Ef = (u[n]*ΔEf + (Ef⁻+Ef⁺))/2 end

            # Compute Legendre angular flux moments
            Σsᵢ = ΔEf .* w[n]/2 .* dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui[δi],Zi[δi],Ti[δi],ri[δi],Ec,incoming_particle) * atz[i]
            if is_dirac Σsᵢ /= ΔEf end
            𝓕i .+= Σsᵢ
            𝓕iₑ += Σsᵢ[1] * (Ef + ΔQ)
        end

        # Add analytic singular contribution
        if is_heavy_inelastic_S
            𝓕i .+= analytic_A .* ones(L+1)
            for l in range(0,L)
                lead_log = integrate_leading_1overW_per_subshell(cache, Ef⁻, Ef⁺, l) * atz[i]
                𝓕i[l+1] += lead_log
            end
            σ_analytic = feed_analytical_heavy_particle(cache, Ef⁻, Ef⁺) * atz[i]
            𝓕iₑ = Ei * σ_analytic - M₁

            if ~isapprox(𝓕i[1], σ_analytic; rtol=1e-3, atol=1e-12)
                rel = abs(𝓕i[1] - σ_analytic) / max(abs(σ_analytic), 1e-300)
                print("FEED_WARN: Z=$(Z[i]), δi=$(δi), gf=$(gf), σ_analytic=$(σ_analytic), numeric_l0=$(𝓕i[1]), rel_diff=$(rel)\n")
                𝓕i[1] = σ_analytic
            end
        end
        𝓕[gf,:] .+= 𝓕i
        𝓕ₑ[gf] += 𝓕iₑ
    end
end
return 𝓕, 𝓕ₑ
end

"""
    feed_elastic_scattering(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,
    Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,
    particles::Vector{Particle},type::String,incoming_particle::Particle,
    scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,
    is_subshells::Bool,A::Vector{Vector{Int64}},
    atpercentA::Vector{Vector{Float64}})

Calculate the elastic-scattering feed function 𝓕 from incident energy `Ei` into each
outgoing energy group and Legendre moment up to order `L`. Also calculate the
energy-weighted feed function 𝓕ₑ for energy-deposition cross sections, including isotope
fractions when isotope-resolved data are provided.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic number of the element(s) composing the material.
- `atz::Vector{Float64}` : atomic percent of the element(s) composing the material.
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : energy of the incoming particle [in mₑc²].
- `Eout::Vector{Float64}` : energy group boundaries [in mₑc²].
- `Ng::Int64` : number of groups.
- `interaction::Interaction` : interaction information.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` : number of groups for the incoming particle.
- `particles::Vector{Particle}` : list of particles involved in the interaction.
- `type::String` : type of interaction (`"S"` for scattering or `"P"` for production).
- `incoming_particle::Particle` : incoming particle.
- `scattered_particle::Particle` : scattered particle.
- `Ein::Vector{Float64}` : energy group boundaries corresponding to the incoming
  particle [in mₑc²].
- `Ec::Float64` : cutoff energy between soft and catastrophic interaction.
- `is_elastic::Bool` : boolean indicating if the outgoing particle energy is equal to `Ei`.
- `is_subshells::Bool` : boolean indicating if subshell-dependent cross sections are used.
- `A::Vector{Vector{Int64}}` : isotope mass numbers per element.
- `atpercentA::Vector{Vector{Float64}}` : isotope atomic fractions per element.

# Output Argument(s)
- `𝓕::Array{Float64}` : feed function.
- `𝓕ₑ::Vector{Float64}` : energy-weighted feed function.

# Reference(s)
- MacFarlane et al. (2021) : The NJOY Nuclear Data Processing System, Version 2012.

"""
function feed_elastic_scattering(Z::Vector{Int64},atz::Vector{Float64},L::Int64,Ei::Float64,Eout::Vector{Float64},Ng::Int64,interaction::Interaction,gi::Int64,Ngi::Int64,particles::Vector{Particle},type::String,incoming_particle::Particle,scattered_particle::Particle,Ein::Vector{Float64},Ec::Float64,is_elastic::Bool,is_subshells::Bool,A::Vector{Vector{Int64}},atpercentA::Vector{Vector{Float64}})

#----
# Initialization
#----
𝓕 = zeros(Ng+1,L+1)
𝓕ₑ = zeros(Ng+1)

#----
# Feed function over all groups and under the cutoff energy
#----

# Loop over the compound elements
Nz = length(Z)
for i in range(1,Nz)
    # Loop over isotopes
    for (Ai, atai) in zip(A[i], atpercentA[i])
        if type == "P" && !(Z[i] == 1 && Ai == 1)
            continue
        end
        δi = 0
        Ui = 0.0
        Zi = Z[i]
        Ti = 0.0
        ri = 0.0
        for gf in range(1,Ng)

            # Final energy group
            Ef⁻ = Eout[gf]; Ef⁺ = Eout[gf+1]
            M_target = get_mass(Z[i], Ai)
            Ef⁻,Ef⁺,isSkip = bounds_dispatch(interaction,Ef⁻,Ef⁺,Ei,gi,gf,type,Ui,Ec,incoming_particle,M_target)
            if isSkip continue end
            ΔEf = Ef⁻ - Ef⁺

            # Integration over the energy group
            𝓕i = zeros(L+1)
            𝓕iₑ = 0
            Ef = is_elastic ? Ei : (Ef⁻ + Ef⁺) / 2
            Σsᵢ = dcs_dispatch(interaction,L,Ei,Ef,Z[i],scattered_particle,type,i,particles,Ein,Ef⁻,Ef⁺,δi,Ui,Zi,Ti,ri,Ec,incoming_particle,Ai) * atz[i] * atai
            𝓕i .+= Σsᵢ
            𝓕iₑ += Σsᵢ[1] * Ef
            𝓕[gf,:] .+= 𝓕i
            𝓕ₑ[gf] += 𝓕iₑ
        end
    end
end
return 𝓕, 𝓕ₑ
end

"""
    spectrum_weights(E_out::Vector{Float64}, LEP::Int)

Quadrature weights integrating a tabulated ENDF MF=6 outgoing-energy spectrum, following
its interpolation law: for `LEP = 1` the spectrum is a histogram, constant over each bin,
and the weight of a point is the width of the bin it opens, the last point closing the
table with a zero weight; for any other law it is interpolated linear-linear and the
trapezoidal weights are used. Both reproduce ∫f dE_out = 1 for a normalized spectrum.

# Input Argument(s)
- `E_out::Vector{Float64}` : tabulated outgoing energies [eV], increasing.
- `LEP::Int` : ENDF interpolation law of the outgoing energy.

# Output Argument(s)
- `ΔE::Vector{Float64}` : quadrature weight of each tabulated point [eV].

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=1 continuum distributions.
"""
function spectrum_weights(E_out::Vector{Float64}, LEP::Int)
    N = length(E_out)
    ΔE = zeros(N)
    N < 2 && return ΔE
    if LEP == 1
        for j in 1:N-1
            ΔE[j] = E_out[j+1] - E_out[j]
        end
    else
        ΔE[1] = (E_out[2] - E_out[1]) / 2
        for j in 2:N-1
            ΔE[j] = (E_out[j+1] - E_out[j-1]) / 2
        end
        ΔE[N] = (E_out[N] - E_out[N-1]) / 2
    end
    return ΔE
end

"""
    feed_nuclear_reaction(Z::Vector{Int64}, atz::Vector{Float64},
    ωz::Vector{Float64}, ρ::Float64, state_of_matter::String,
    L::Int64, Ei::Float64, Eout::Vector{Float64}, Ng::Int64,
    interaction::Nuclear_Reaction, gi::Int64, Ngi::Int64,
    particles::Vector{Particle}, type::String,
    incoming_particle::Particle, Ein::Vector{Float64},
    I_eff::Float64,
    A::Union{Nothing,Vector{Vector{Int64}}},
    atpercentA::Union{Nothing,Vector{Vector{Float64}}})

Calculate the nuclear-reaction feed function 𝓕 for type `"P"` (equivalent-proton
production following Salvat & Quesada 2020 §4.2). For each charged-particle production
channel with a Kalbach-Mann energy distribution, each histogram bin in E_b (product
kinetic energy) is mapped to an equivalent proton energy E_eq via CSDA range-matching
with weight w = E_b/E_eq, then accumulated into the appropriate outgoing energy group.

The angular distribution is isotropic (Legendre moment L=0 only; higher moments are zero).

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements composing the material.
- `atz::Vector{Float64}` : atomic fraction of each element.
- `ωz::Vector{Float64}` : weight fraction of each element.
- `ρ::Float64` : material density [g/cm³].
- `state_of_matter::String` : material state (`"solid"`, `"liquid"`, or `"gas"`).
- `L::Int64` : Legendre truncation order.
- `Ei::Float64` : incoming particle energy [mₑc²].
- `Eout::Vector{Float64}` : outgoing energy group boundaries [mₑc²], decreasing.
- `Ng::Int64` : number of outgoing energy groups.
- `interaction::Nuclear_Reaction` : nuclear-reaction interaction structure.
- `gi::Int64` : incoming particle group index.
- `Ngi::Int64` : number of incoming particle groups.
- `particles::Vector{Particle}` : particles in the simulation.
- `type::String` : interaction type string (`"P"`).
- `incoming_particle::Particle` : incoming particle.
- `Ein::Vector{Float64}` : incoming energy group boundaries [mₑc²].
- `I_eff::Float64` : effective mean excitation energy override [mₑc²]; NaN ⟹ tables.
- `A::Union{Nothing,Vector{Vector{Int64}}}` : isotope mass numbers per element.
- `atpercentA::Union{Nothing,Vector{Vector{Float64}}}` : isotope atomic fractions per element.

# Output Argument(s)
- `𝓕::Array{Float64,2}` : feed function [cm²], size (Ng+1, L+1).
- `𝓕ₑ::Vector{Float64}` : energy-weighted feed function [cm² × mₑc²], size (Ng+1).

# Reference(s)
- Salvat & Quesada (2020), NIMB 475, 49–62, §4.2 equivalent-proton method.
- MacFarlane et al. (2021), The NJOY Nuclear Data Processing System, Version 2012.
"""
function feed_nuclear_reaction(Z::Vector{Int64}, atz::Vector{Float64},
        ωz::Vector{Float64}, ρ::Float64, state_of_matter::String,
        L::Int64, Ei::Float64, Eout::Vector{Float64}, Ng::Int64,
        interaction::Nuclear_Reaction, _gi::Int64, _Ngi::Int64,
        _particles::Vector{Particle}, _type::String,
        incoming_particle::Particle, _Ein::Vector{Float64},
        I_eff::Float64,
        A::Union{Nothing,Vector{Vector{Int64}}},
        atpercentA::Union{Nothing,Vector{Vector{Float64}}})

    mₑc² = M_E_C2_MEV   # 0.510999 MeV
    eV_per_mc2 = mₑc² * 1e6  # eV per mₑc²

    𝓕  = zeros(Ng+1, L+1)
    𝓕ₑ = zeros(Ng+1)

    ptype = get_type(incoming_particle)
    haskey(interaction.production_db, ptype) || return 𝓕, 𝓕ₑ

    E_in_eV = Ei * eV_per_mc2

    Nz = length(Z)
    for i in 1:Nz
        Ai_vec         = isnothing(A)           ? nothing : A[i]
        atpercentAi_vec = isnothing(atpercentA) ? nothing : atpercentA[i]

        # Collect (mass_number, isotope_atomic_fraction) pairs for this element
        iso_pairs = if isnothing(Ai_vec)
            isotopic_composition(Z[i])
        else
            collect(zip(Ai_vec, atpercentAi_vec))
        end

        for (Aii, ataiii) in iso_pairs
            key = (Z[i], Aii)
            haskey(interaction.production_db[ptype], key) || continue
            channels = interaction.production_db[ptype][key]
            isempty(channels) && continue

            # Build equivalent-proton table lazily (cached per ZAP × material)
            # Loop over channels, group by ZAP to share the cache per ZAP
            for channel in channels
                channel.kalbach === nothing && continue

                # Total channel cross-section at Ei [barn]
                σ_ch = interp_TAB1(E_in_eV, channel.E_xs, channel.S_xs,
                                   channel.NBT_xs, channel.INT_xs)
                σ_ch <= 0.0 && continue

                # Yield (products per reaction) at Ei
                y_ch = if isempty(channel.yield_E)
                    1.0
                else
                    interp_TAB1(E_in_eV, channel.yield_E, channel.yield_y,
                                Int[], Int[])
                end
                y_ch <= 0.0 && continue

                # Get or build the equivalent-proton lookup table for this (ZAP, material)
                (E_b_tbl, E_eq_tbl, w_tbl) = get_or_build_eq_cache!(
                    interaction, channel.ZAP, Z, ωz, ρ, state_of_matter, I_eff)

                # Find the nearest lower incident-energy index in the Kalbach table
                # (histogram interpolation in E_in: use spectrum from E_in[k])
                kalbach = channel.kalbach
                k_ein = searchsortedfirst(kalbach.E_in, E_in_eV) - 1
                k_ein = clamp(k_ein, 1, length(kalbach.E_in))

                E_out_k = kalbach.E_out[k_ein]
                f_k     = kalbach.f[k_ein]
                isempty(E_out_k) && continue

                # Quadrature weights of the tabulated spectrum, following its ENDF
                # interpolation law: histogram bins for LEP=1, trapezoids for LEP=2. Both
                # integrate f over the table exactly, so that Σⱼ f[j] × ΔE_b[j] = 1.
                ΔE_b = spectrum_weights(E_out_k, kalbach.LEP)

                # Iterate over the tabulated product kinetic energies
                for j in eachindex(f_k)
                    ΔE_b_eV = ΔE_b[j]
                    ΔE_b_eV <= 0.0 && continue

                    E_b_eV = E_out_k[j]
                    f_j    = f_k[j]  # [eV⁻¹], normalized spectrum
                    f_j <= 0.0 && continue

                    E_b_mc2 = E_b_eV / eV_per_mc2

                    # Equivalent proton energy and weight from precomputed table
                    E_eq_mc2 = interp_linear(E_b_tbl, E_eq_tbl, E_b_mc2)
                    w_j      = interp_linear(E_b_tbl, w_tbl, E_b_mc2)

                    # Equivalent proton energy must be positive and within output grid
                    E_eq_mc2 <= 0.0 && continue

                    # Find which outgoing group contains E_eq_mc2
                    # Eout is decreasing: Eout[gf] ≥ E_eq > Eout[gf+1]
                    gf = searchsortedlast(Eout, E_eq_mc2; rev=true)
                    (gf < 1 || gf > Ng) && continue

                    # Contribution to feed [cm²]:
                    # σ_ch [barn] × y_ch × f_j [eV⁻¹] × ΔE_b [eV] × w_j × atz × atai
                    contrib = σ_ch * BARN_TO_CM2 * y_ch * f_j * ΔE_b_eV * w_j *
                              atz[i] * ataiii

                    # Isotropic: only L=0 moment is non-zero
                    𝓕[gf, 1] += contrib
                    # Energy-weighted feed [cm² × mₑc²]:
                    # w_j × E_eq = E_b (energy balance), so contrib × E_b_mc2
                    𝓕ₑ[gf] += contrib * E_b_mc2 / w_j  # = contrib × E_eq_mc2
                    # equivalently: σ × y × f × ΔE_b × w × E_eq = σ × y × f × ΔE_b × E_b
                end
            end
        end
    end

    return 𝓕, 𝓕ₑ
end
