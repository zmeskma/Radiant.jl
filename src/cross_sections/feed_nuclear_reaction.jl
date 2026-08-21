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

                # Iterate over histogram bins in product kinetic energy
                n_bins = length(f_k)
                for j in 1:n_bins
                    # Bin width: use trapezoidal intervals between E_out points
                    if j < length(E_out_k)
                        ΔE_b_eV = E_out_k[j+1] - E_out_k[j]
                    else
                        ΔE_b_eV = j > 1 ? (E_out_k[j] - E_out_k[j-1]) : 0.0
                    end
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
