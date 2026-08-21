"""
    nuclear_reaction_endf(db_name::String, isotopes::Vector{Tuple{Int,Int}},
    particle::Particle, mt::Int; data_root::Union{Nothing,String}=nothing)

Builds the nuclear-reaction ENDF database for a set of target isotopes and one incident
particle, reading the MF=3 total cross-section curve for the requested reaction (MT)
from each isotope's ENDF file.

# Input Argument(s)
- `db_name::String` : ENDF database directory name.
- `isotopes::Vector{Tuple{Int,Int}}` : list of `(Z, A)` target isotopes.
- `particle::Particle` : incident particle type.
- `mt::Int` : ENDF reaction number (MT) to read from MF=3.
- `data_root::Union{Nothing,String}` : optional root directory containing ENDF databases.

# Output Argument(s)
- `db::NuclearReactionENDFDB` : nuclear-reaction ENDF database for the requested particle
  and MT.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 reaction cross-sections.
"""
function nuclear_reaction_endf(db_name::String, isotopes::Vector{Tuple{Int,Int}}, particle::Particle, mt::Int; data_root::Union{Nothing,String}=nothing)
    root = isnothing(data_root) ? normpath(joinpath(@__DIR__, "..", "..", "data")) : data_root
    db = Dict{Tuple{Int,Int},IsotopeNuclearReactionDB}()
    for (Z, A) in isotopes
        endf_path = endf_path_for_isotope(db_name, Z, A, particle; data_root=root)
        if !isfile(endf_path)
            error("ENDF file not found: $(endf_path)")
        end
        lines = readlines(endf_path)
        (ZA, AWR, E, S, NBT, INT) = read_mf3_mt(lines, mt)
        db[(Z, A)] = IsotopeNuclearReactionDB(E, S, NBT, INT)
    end
    return NuclearReactionENDFDB(mt, db)
end

# ZAP → (charge Z_b, rest mass in MeV/c²)
const _PRODUCT_ZM = Dict{Int,Tuple{Int,Float64}}(
    1001 => (1, 938.272),    # proton
    1002 => (1, 1875.613),   # deuteron
    1003 => (1, 2808.921),   # triton
    2003 => (2, 2808.391),   # He-3
    2004 => (2, 3727.379),   # alpha
)
const _CHARGED_ZAPS = Set{Int}([1001, 1002, 1003, 2003, 2004])

"""
    nuclear_production_endf(db_name::String, isotopes::Vector{Tuple{Int,Int}},
    particle::Particle, production_mts::Vector{Int};
    data_root::Union{Nothing,String}=nothing)

Reads MF=3 cross-sections and MF=6 Kalbach-Mann energy distributions for all requested
charged-particle production channels from each isotope's ENDF file. Only LAW=1 LANG=2
(Kalbach-Mann) product subsections are retained; other distributions are silently skipped.

# Input Argument(s)
- `db_name::String` : ENDF database directory name.
- `isotopes::Vector{Tuple{Int,Int}}` : list of `(Z, A)` target isotopes.
- `particle::Particle` : incident particle type.
- `production_mts::Vector{Int}` : ENDF MT numbers to search for charged-particle products.
- `data_root::Union{Nothing,String}` : optional root directory containing ENDF databases.

# Output Argument(s)
- `db::Dict{Tuple{Int,Int},Vector{IsotopeProductionChannelDB}}` : per-isotope list of
  charged-particle production channels.

# Reference(s)
- ENDF-6 Formats Manual, MF=6 product energy-angle distributions, LAW=1 LANG=2.
- Salvat & Quesada (2020), NIMB 475, 49–62, §4.2 equivalent-proton method.
"""
function nuclear_production_endf(db_name::String, isotopes::Vector{Tuple{Int,Int}},
        particle::Particle, production_mts::Vector{Int};
        data_root::Union{Nothing,String}=nothing)

    root = isnothing(data_root) ? normpath(joinpath(@__DIR__, "..", "..", "data")) : data_root
    db = Dict{Tuple{Int,Int}, Vector{IsotopeProductionChannelDB}}()

    for (Z, A) in isotopes
        endf_path = endf_path_for_isotope(db_name, Z, A, particle; data_root=root)
        if !isfile(endf_path)
            db[(Z, A)] = IsotopeProductionChannelDB[]
            continue
        end
        lines = readlines(endf_path)

        channels = IsotopeProductionChannelDB[]

        for mt in production_mts
            # Check MF=3 MT exists
            i_mf3 = find_first_mf(lines, 3)
            i_mf3 === nothing && continue
            find_section_start(lines, 3, mt; start_i=i_mf3) === nothing && continue

            # Read MF=3 total XS for this channel
            (_, _, E_xs, S_xs, NBT_xs, INT_xs) = read_mf3_mt(lines, mt)

            # Check MF=6 MT exists
            i_mf6 = find_first_mf(lines, 6)
            i_mf6 === nothing && continue
            find_section_start(lines, 6, mt; start_i=i_mf6) === nothing && continue

            # Read product subsections
            products = try
                read_mf6_mt(lines, mt)
            catch e
                @warn "nuclear_production_endf: failed MF6 MT=$mt for (Z=$Z, A=$A): $e"
                continue
            end

            for prod in products
                ZAP_i = round(Int, prod.ZAP)
                ZAP_i ∈ _CHARGED_ZAPS || continue
                prod.law1_LANG === nothing && continue
                prod.law1_LANG == 2       || continue
                isempty(prod.law1_E_in)   && continue

                kalbach = KalbachMannTable(
                    prod.law1_E_in, prod.law1_NBT_Ei, prod.law1_INT_Ei,
                    prod.law1_E_out, prod.law1_f, prod.law1_r)

                push!(channels, IsotopeProductionChannelDB(
                    mt, ZAP_i, prod.AWP,
                    E_xs, S_xs, NBT_xs, INT_xs,
                    prod.Eint, prod.yi,
                    kalbach))
            end
        end

        db[(Z, A)] = channels
    end

    return db
end

"""
    _build_proton_range_table(E_min_mc2::Float64, E_max_mc2::Float64, N_pts::Int,
    Z::Vector{Int64}, ωz::Vector{Float64}, ρ::Float64,
    state_of_matter::String, I_eff::Float64, proton::Particle)

Builds a cumulative proton CSDA-range table on a log-spaced energy grid using the
corrected Bethe stopping-power formula. The density does not need to be physically exact
for range-matching purposes; it is passed to `bethe` only to compute the density-effect
correction consistently for both product and proton.

Returns `(E_grid_mc2, R_grid_cm)` — incident energy [mₑc²] and cumulative CSDA range
[cm], with R_grid_cm[1] = 0.
"""
function _build_proton_range_table(E_min_mc2::Float64, E_max_mc2::Float64, N_pts::Int,
        Z::Vector{Int64}, ωz::Vector{Float64}, ρ::Float64,
        state_of_matter::String, I_eff::Float64, proton::Particle)

    E_grid = exp.(range(log(E_min_mc2), log(E_max_mc2), length=N_pts))
    R_grid = zeros(N_pts)

    # 8-point Gauss-Legendre for each sub-interval
    u_gl, w_gl = quadrature(8, "gauss-legendre")

    for k in 2:N_pts
        E_lo = E_grid[k-1]
        E_hi = E_grid[k]
        half = (E_hi - E_lo) / 2.0
        mid  = (E_hi + E_lo) / 2.0
        dR = 0.0
        for n in eachindex(u_gl)
            E_n = mid + half * u_gl[n]
            S_n = bethe(Z, ωz, ρ, E_n, proton, "fano", state_of_matter, I_eff)
            S_n > 0.0 && (dR += w_gl[n] * half / S_n)
        end
        R_grid[k] = R_grid[k-1] + dR
    end

    return E_grid, R_grid
end

"""
    _interp_inv_monotone(x_grid::Vector{Float64}, y_grid::Vector{Float64},
    y_query::Float64)

Inverse interpolation on a monotone increasing (x, y) table: returns x such that
y(x) ≈ y_query, clamped to the grid endpoints.
"""
function _interp_inv_monotone(x_grid::Vector{Float64}, y_grid::Vector{Float64},
        y_query::Float64)
    y_query <= y_grid[1]   && return x_grid[1]
    y_query >= y_grid[end] && return x_grid[end]
    k = searchsortedfirst(y_grid, y_query)
    k = clamp(k, 2, length(y_grid))
    t = (y_query - y_grid[k-1]) / (y_grid[k] - y_grid[k-1])
    return x_grid[k-1] + t * (x_grid[k] - x_grid[k-1])
end

"""
    get_or_build_eq_cache!(interaction::Nuclear_Reaction, ZAP::Int,
    Z::Vector{Int64}, ωz::Vector{Float64}, ρ::Float64,
    state_of_matter::String, I_eff::Float64)

Returns the equivalent-proton lookup table `(E_b_grid, E_eq_grid, weight_grid)` [all in
mₑc²] for the given product ZAP in the specified material, building and caching it on
first call.

The range-scaling relation used is the standard Bethe-based formula:
    R_product(E_b) = (M_b/M_p / Z_b²) × R_proton(E_b × M_p/M_b)
where M_p is the proton rest mass and the equality of γ at the same kinetic energy per
nucleon is used (exact in the non-relativistic limit).
"""
function get_or_build_eq_cache!(interaction, ZAP::Int,
        Z::Vector{Int64}, ωz::Vector{Float64}, ρ::Float64,
        state_of_matter::String, I_eff::Float64)

    cache_key = hash((ZAP, Z, round.(ωz, digits=6),
                      round(ρ, digits=4),
                      isnan(I_eff) ? 0.0 : round(I_eff, digits=4)))

    haskey(interaction.production_eq_cache, cache_key) &&
        return interaction.production_eq_cache[cache_key]

    mₑc²   = M_E_C2_MEV          # 0.510999 MeV
    M_p_MeV = 938.272
    proton  = Proton()

    # ---- Proton CSDA range table ----
    E_min_mc2 = 1e-4 / mₑc²      # 100 eV
    E_max_mc2 = 350.0 / mₑc²     # 350 MeV (covers all product energies)
    N_range   = 300
    E_p_grid, R_p_grid = _build_proton_range_table(
        E_min_mc2, E_max_mc2, N_range, Z, ωz, ρ, state_of_matter, I_eff, proton)

    # ---- Equivalent-proton table for this ZAP ----
    N_b       = 200
    E_b_min   = 1e-3 / mₑc²      # 1 keV
    E_b_max   = 300.0 / mₑc²     # 300 MeV
    E_b_grid  = exp.(range(log(E_b_min), log(E_b_max), length=N_b))
    E_eq_grid = similar(E_b_grid)
    w_grid    = similar(E_b_grid)

    if ZAP == 1001
        # Product is already a proton — identity mapping
        E_eq_grid .= E_b_grid
        w_grid    .= 1.0
    else
        Z_b, M_b_MeV = _PRODUCT_ZM[ZAP]
        M_b_over_M_p = M_b_MeV / M_p_MeV  # ≈ A_b (mass number)
        scale_factor  = M_b_over_M_p / Z_b^2

        for k in eachindex(E_b_grid)
            E_b = E_b_grid[k]
            # Proton energy at same γ (same velocity) as product at E_b
            E_p_same = E_b / M_b_over_M_p
            # Product CSDA range via scaling
            R_prod = scale_factor * interp_linear(E_p_grid, R_p_grid, E_p_same)
            # Equivalent proton energy by inverse interpolation
            E_eq = _interp_inv_monotone(E_p_grid, R_p_grid, R_prod)
            E_eq_grid[k] = E_eq
            w_grid[k]    = E_b / max(E_eq, E_b * 1e-9)
        end
    end

    result = (E_b_grid, E_eq_grid, w_grid)
    interaction.production_eq_cache[cache_key] = result
    return result
end
