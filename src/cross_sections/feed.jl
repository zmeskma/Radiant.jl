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
    eval_spectrum(x::Vector{Float64}, g::Vector{Float64}, xq::Float64, LEP::Int)

Evaluates a tabulated spectrum at one abscissa, following its ENDF interpolation law:
constant over the bin it opens for `LEP = 1`, linear between its nodes otherwise. Outside
the table the spectrum is zero, the outgoing energies beyond the tabulated range being
kinematically inaccessible.

# Input Argument(s)
- `x::Vector{Float64}` : tabulated abscissa, increasing.
- `g::Vector{Float64}` : tabulated spectrum.
- `xq::Float64` : abscissa at which the spectrum is evaluated.
- `LEP::Int` : ENDF interpolation law of the outgoing energy.

# Output Argument(s)
- `gq::Float64` : value of the spectrum at `xq`.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=1 continuum distributions.
"""
function eval_spectrum(x::Vector{Float64}, g::Vector{Float64}, xq::Float64, LEP::Int)
    n = length(x)
    n == 0 && return 0.0
    (xq < x[1] || xq > x[n]) && return 0.0
    n == 1 && return g[1]
    j = clamp(searchsortedlast(x, xq), 1, n - 1)
    LEP == 1 && return g[j]
    Δ = x[j+1] - x[j]
    Δ <= 0.0 && return g[j]
    t = (xq - x[j]) / Δ
    return g[j] + t * (g[j+1] - g[j])
end

"""
    interpolate_spectrum(spectrum::KalbachMannTable, E_in::Float64)

Gives the outgoing-energy spectrum of a product at an arbitrary incident energy [eV], by
interpolating between the two tabulated spectra that bracket it.

The interpolation is done on a unit base, as ENDF prescribes for LAW=1 continuum
distributions: the outgoing-energy range of a spectrum grows with the incident energy — a
proton of 80 MeV on Fe-56 emits up to 77 MeV, one of 100 MeV up to 96 MeV — so
interpolating the two spectra at a fixed outgoing energy would let an intermediate incident
energy emit above its own kinematic limit. Each spectrum is therefore mapped onto
x = (E_out − E_min)/(E_max − E_min) ∈ [0,1] with the density rescaled accordingly,
the two are combined at constant x, and the result is mapped back onto the interpolated
range. The transform preserves ∫f dE_out = 1 exactly.

The weight given to each bracketing spectrum follows the interpolation law of the TAB2
record on the incident-energy axis: the lower spectrum is kept as it stands for a histogram
law, and the two are blended linearly, or logarithmically in energy, otherwise.

# Input Argument(s)
- `spectrum::KalbachMannTable` : tabulated spectra of one production channel.
- `E_in::Float64` : incident energy [eV].

# Output Argument(s)
- `E_out::Vector{Float64}` : outgoing energy grid at `E_in` [eV].
- `f::Vector{Float64}` : probability density on that grid [eV⁻¹], of unit integral.

# Reference(s)
- ENDF-6 Formats Manual, MF=6, LAW=1, interpolation between incident energies.
- MacFarlane et al. (2021), The NJOY Nuclear Data Processing System, unit-base
  interpolation of continuum distributions.
"""
function interpolate_spectrum(spectrum::KalbachMannTable, E_in::Float64)
    E_grid = spectrum.E_in
    n = length(E_grid)
    n == 0 && return Float64[], Float64[]

    # Bracketing spectra. searchsortedlast gives the last tabulated energy at or below
    # E_in, so an incident energy falling exactly on a node uses that node's spectrum.
    k = clamp(searchsortedlast(E_grid, E_in), 1, n)
    (k == n || E_in <= E_grid[1]) && return spectrum.E_out[k], spectrum.f[k]

    E_lo, E_hi = spectrum.E_out[k],   spectrum.E_out[k+1]
    f_lo, f_hi = spectrum.f[k],       spectrum.f[k+1]
    isempty(E_lo) && return E_hi, f_hi
    isempty(E_hi) && return E_lo, f_lo

    # Weight of the upper spectrum, from the law of the incident-energy axis.
    law = select_interpolation(spectrum.NBT_Ei, spectrum.INT_Ei, k)
    law == 1 && return E_lo, f_lo
    t = if law == 3 || law == 5
        (log(E_in) - log(E_grid[k])) / (log(E_grid[k+1]) - log(E_grid[k]))
    else
        (E_in - E_grid[k]) / (E_grid[k+1] - E_grid[k])
    end
    t = clamp(t, 0.0, 1.0)
    t == 0.0 && return E_lo, f_lo   # incident energy on a tabulated node
    t == 1.0 && return E_hi, f_hi

    # Unit-base transform of both spectra.
    span_lo = E_lo[end] - E_lo[1]
    span_hi = E_hi[end] - E_hi[1]
    (span_lo <= 0.0 || span_hi <= 0.0) && return E_lo, f_lo
    x_lo = (E_lo .- E_lo[1]) ./ span_lo
    x_hi = (E_hi .- E_hi[1]) ./ span_hi
    g_lo = f_lo .* span_lo
    g_hi = f_hi .* span_hi

    # Common abscissa, so that neither spectrum loses a node.
    x = sort!(unique!(vcat(x_lo, x_hi)))
    LEP = spectrum.LEP
    g = similar(x)
    for (j, xj) in enumerate(x)
        g[j] = (1 - t) * eval_spectrum(x_lo, g_lo, xj, LEP) +
                    t  * eval_spectrum(x_hi, g_hi, xj, LEP)
    end

    # Back to the interpolated outgoing-energy range.
    E_min = (1 - t) * E_lo[1]   + t * E_hi[1]
    E_max = (1 - t) * E_lo[end] + t * E_hi[end]
    span  = E_max - E_min
    span <= 0.0 && return E_lo, f_lo

    return E_min .+ x .* span, g ./ span
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
channel with a Kalbach-Mann energy distribution, the outgoing spectrum is interpolated at
the incident energy, then each of its intervals in E_b (product kinetic energy) is mapped
to an interval in equivalent-proton energy E_eq by CSDA range-matching, and integrated
over the outgoing groups that interval overlaps, with the weight w = E_b/E_eq.

Mapping intervals rather than points is what keeps the result a spectrum: the content of a
source interval is shared between the groups it covers, in proportion to the part of the
interval each group receives, so the transfer matrix stays smooth however fine the group
structure is. Both the number of equivalent protons and the energy they carry are
conserved exactly, whatever the outgoing grid.

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
                (E_b_tbl, E_eq_tbl, _) = get_or_build_eq_cache!(
                    interaction, channel.ZAP, Z, ωz, ρ, state_of_matter, I_eff)

                # Spectrum at this very incident energy, interpolated on a unit base
                # between the two tabulated spectra bracketing it.
                kalbach = channel.kalbach
                E_out_k, f_k = interpolate_spectrum(kalbach, E_in_eV)
                isempty(E_out_k) && continue

                # Each interval of the tabulated spectrum, not each of its points, is
                # transported: the interval [E_b(j), E_b(j+1)] maps to the equivalent-proton
                # interval [E_eq(j), E_eq(j+1)], and its content is shared between the
                # outgoing groups that interval overlaps. Sending it instead to the single
                # group holding one mapped point leaves a comb of empty and overfull groups
                # as soon as the group structure is finer than the ENDF spectrum.
                σy = σ_ch * BARN_TO_CM2 * y_ch * atz[i] * ataiii
                for j in 1:length(E_out_k)-1
                    E_b_lo_eV, E_b_hi_eV = E_out_k[j], E_out_k[j+1]
                    E_b_hi_eV <= E_b_lo_eV && continue
                    f_j = f_k[j]
                    # A histogram bin is empty when its own value is zero, but a
                    # linear-linear one rising from zero carries mass and must be kept.
                    f_hi = kalbach.LEP == 1 ? f_j : f_k[j+1]
                    (f_j <= 0.0 && f_hi <= 0.0) && continue

                    E_eq_lo = interp_linear(E_b_tbl, E_eq_tbl, E_b_lo_eV / eV_per_mc2)
                    E_eq_hi = interp_linear(E_b_tbl, E_eq_tbl, E_b_hi_eV / eV_per_mc2)
                    (E_eq_lo <= 0.0 || E_eq_hi <= E_eq_lo) && continue

                    # Groups overlapped by the mapped interval. Eout decreases, so the
                    # group g spans [Eout[g+1], Eout[g]].
                    g_hi = clamp(searchsortedlast(Eout, E_eq_lo; rev=true), 1, Ng)
                    g_lo = clamp(searchsortedlast(Eout, E_eq_hi; rev=true), 1, Ng)
                    for g in g_lo:g_hi
                        a = max(E_eq_lo, Eout[g+1])
                        b = min(E_eq_hi, Eout[g])
                        b <= a && continue

                        # Back to the product energies feeding this group, so that the
                        # spectrum is integrated over the source interval it comes from.
                        # The two ends of the interval are known exactly and must not be
                        # taken through the inverse map: rounding there would shave a
                        # sliver off every interval and lose part of the spectrum.
                        E_b_a = a <= E_eq_lo ? E_b_lo_eV :
                            clamp(interp_linear(E_eq_tbl, E_b_tbl, a) * eV_per_mc2,
                                  E_b_lo_eV, E_b_hi_eV)
                        E_b_b = b >= E_eq_hi ? E_b_hi_eV :
                            clamp(interp_linear(E_eq_tbl, E_b_tbl, b) * eV_per_mc2,
                                  E_b_lo_eV, E_b_hi_eV)
                        ΔE_b_eV = E_b_b - E_b_a
                        ΔE_b_eV <= 0.0 && continue

                        # The spectrum is linear over the sub-interval, so ∫f dE_b is the
                        # trapezoidal rule and ∫f E_b dE_b, whose integrand is quadratic,
                        # has the closed form below. For a histogram spectrum (LEP=1)
                        # f_hi equals f_j and both reduce to the exact rectangle rules.
                        u_a = (E_b_a - E_b_lo_eV) / (E_b_hi_eV - E_b_lo_eV)
                        u_b = (E_b_b - E_b_lo_eV) / (E_b_hi_eV - E_b_lo_eV)
                        f_a = f_j + u_a * (f_hi - f_j)
                        f_b = f_j + u_b * (f_hi - f_j)
                        mass = 0.5 * (f_a + f_b) * ΔE_b_eV

                        df = f_b - f_a
                        mass_E = ΔE_b_eV * (f_a * E_b_a + (f_a * ΔE_b_eV + df * E_b_a) / 2 +
                                            df * ΔE_b_eV / 3) / eV_per_mc2

                        E_b_mid = 0.5 * (E_b_a + E_b_b) / eV_per_mc2
                        E_eq_mid = 0.5 * (a + b)
                        w_mid = E_b_mid / E_eq_mid

                        # Isotropic: only the L=0 moment is non-zero.
                        𝓕[g, 1] += σy * mass * w_mid
                        # Energy-weighted feed [cm² × mₑc²]: w × E_eq = E_b, so the
                        # energy the equivalent protons carry is the one of the product.
                        𝓕ₑ[g] += σy * mass_E
                    end
                end
            end
        end
    end

    return 𝓕, 𝓕ₑ
end
