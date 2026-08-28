"""
    MF3_SUM_RULES

ENDF-6 redundancy rules for MF=3 cross sections, mapping each redundant (summation)
reaction to the reactions it already contains. Used when the non-elastic total has to be
rebuilt from partial channels, so that a channel covered by a coarser one present in the
same evaluation is not counted twice.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 sum rules for redundant reactions.
"""
const MF3_SUM_RULES = Dict{Int,Vector{Int}}(
    27  => [18, 101],           # absorption = fission + disappearance
    101 => collect(102:117),    # disappearance = capture and charged-particle emission
    4   => collect(50:91),      # (z,n) total = discrete levels + continuum
    16  => collect(875:891),    # (z,2n) total = discrete levels
    18  => [19, 20, 21, 38],    # fission = first- to third-chance + (n,3nf)
    103 => collect(600:649),    # (z,p) total = discrete levels + continuum
    104 => collect(650:699),    # (z,d) total
    105 => collect(700:749),    # (z,t) total
    106 => collect(750:799),    # (z,3He) total
    107 => collect(800:849),    # (z,alpha) total
)

"""
    is_nonelastic_reaction_mt(mt::Int)

Tells whether an MT number designates a non-elastic reaction cross-section that may enter
the non-elastic total. The total (MT=1), the elastic channel (MT=2), the non-elastic total
itself (MT=3), the resonance-parameter section (MT=151) and the auxiliary quantities above
MT=200 (particle-production yields, mu-bar, xi, fission-energy data) are excluded.

# Input Argument(s)
- `mt::Int` : ENDF reaction number.

# Output Argument(s)
- `is_nonelastic::Bool` : true when the reaction is a non-elastic channel cross-section.

# Reference(s)
- ENDF-6 Formats Manual, MT reaction designations.
"""
function is_nonelastic_reaction_mt(mt::Int)
    mt in (1, 2, 3, 151) && return false
    (4 <= mt <= 200) && return true
    (600 <= mt <= 849) && return true
    (875 <= mt <= 891) && return true
    return false
end

"""
    nonelastic_partial_mts(mts_present::Vector{Int})

Selects, among the MF=3 sections of an evaluation, the non-elastic channels whose sum
reproduces the non-elastic total without double counting. Whenever a redundant reaction is
present, the channels it already contains, and their own descendants, are discarded in
favour of the coarser one.

# Input Argument(s)
- `mts_present::Vector{Int}` : MT numbers found in the MF=3 file of the evaluation.

# Output Argument(s)
- `mts::Vector{Int}` : sorted MT numbers to be summed into the non-elastic total.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 sum rules for redundant reactions.
"""
function nonelastic_partial_mts(mts_present::Vector{Int})
    candidates = Set{Int}(mt for mt in mts_present if is_nonelastic_reaction_mt(mt))

    # Mark every channel already contained in a redundant reaction present in the file.
    covered = Set{Int}()
    queue = Int[mt for mt in candidates if haskey(MF3_SUM_RULES, mt)]
    while !isempty(queue)
        parent = pop!(queue)
        for child in MF3_SUM_RULES[parent]
            child in covered && continue
            push!(covered, child)
            haskey(MF3_SUM_RULES, child) && push!(queue, child)
        end
    end

    return sort!(collect(setdiff(candidates, covered)))
end

"""
    sum_mf3_partials(lines::Vector{String}, mts::Vector{Int}, source::AbstractString)

Sums MF=3 cross-section tables over a set of reactions. Each channel is interpolated with
its own TAB1 law on the union of the channel energy grids, and contributes nothing below
its own threshold. The sum itself is tabulated linear-linear on that union grid, which is
exact at the grid points of every channel.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `mts::Vector{Int}` : MT numbers to be summed.
- `source::AbstractString` : source label used in error messages.

# Output Argument(s)
- `E::Vector{Float64}` : union energy grid [eV].
- `S::Vector{Float64}` : summed cross-section [barn].

# Reference(s)
- ENDF-6 Formats Manual, MF=3 reaction cross sections.
"""
function sum_mf3_partials(lines::Vector{String}, mts::Vector{Int}, source::AbstractString)
    channels = [read_mf3_mt(lines, mt; source=source) for mt in mts]

    E = Float64[]
    for (_ZA, _AWR, Ec, _Sc, _NBTc, _INTc) in channels
        append!(E, Ec)
    end
    sort!(E)
    unique!(E)

    S = zeros(length(E))
    for (_ZA, _AWR, Ec, Sc, NBTc, INTc) in channels
        isempty(Ec) && continue
        for (k, Ek) in enumerate(E)
            Ek < Ec[1] && continue # below the channel threshold
            S[k] += interp_TAB1(Ek, Ec, Sc, NBTc, INTc)
        end
    end

    return E, S
end

"""
    zero_mf3_table()

Builds a MF=3-like cross-section table that is identically zero over the whole energy
range, used when an evaluation carries no usable non-elastic data.

# Output Argument(s)
- `table::Tuple` : energy grid [eV], cross-section values [barn], and TAB1 interpolation
  metadata (linear-linear).
"""
function zero_mf3_table()
    return ([1.0e-5, 1.0e11], [0.0, 0.0], [2], [2])
end

"""
    mf3_nonelastic_table(lines::Vector{String}, mt::Int, source::AbstractString)

Provides the MF=3 cross-section table of the requested reaction, falling back on the
partial channels when the requested reaction is the non-elastic total (MT=3) and the
evaluation does not tabulate it:
- MF=3 MT=`mt` present: the tabulated cross-section is returned;
- MT=3 requested but absent, with non-elastic partial channels available: the non-elastic
  total is rebuilt as their sum, and the channels used are reported;
- no MF=3 data at all, or no non-elastic channel besides elastic scattering: the
  cross-section is set to zero and a warning is issued.

# Input Argument(s)
- `lines::Vector{String}` : ENDF file lines.
- `mt::Int` : requested MT reaction number.
- `source::AbstractString` : isotope and library label used in the messages.

# Output Argument(s)
- `table::Tuple` : energy grid [eV], cross-section values [barn], and TAB1 interpolation
  metadata.

# Reference(s)
- ENDF-6 Formats Manual, MF=3 reaction cross sections.
"""
function mf3_nonelastic_table(lines::Vector{String}, mt::Int, source::AbstractString)
    i_mf3 = find_first_mf(lines, 3)
    if isnothing(i_mf3)
        @warn "No MF=3 data found for $(source): the nuclear-reaction cross-section is set to zero."
        return zero_mf3_table()
    end

    if !isnothing(find_section_start(lines, 3, mt; start_i=i_mf3))
        (_ZA, _AWR, E, S, NBT, INT) = read_mf3_mt(lines, mt; source=source)
        return (E, S, NBT, INT)
    end

    if mt != 3
        error("MF=3 MT=$(mt) not found in $(source); reconstruction from partial channels is only defined for the non-elastic total (MT=3).")
    end

    mts = nonelastic_partial_mts(list_mf_mts(lines, 3))
    if isempty(mts)
        @warn "No MF=3 MT=3 and no non-elastic partial channel for $(source): the nuclear-reaction cross-section is set to zero."
        return zero_mf3_table()
    end

    @info "No MF=3 MT=3 for $(source): the non-elastic total is rebuilt as the sum of the MF=3 partial channels MT=$(mts)."
    if length(mts) == 1
        # A single channel is kept with its own grid and interpolation law.
        (_ZA, _AWR, E, S, NBT, INT) = read_mf3_mt(lines, mts[1]; source=source)
        return (E, S, NBT, INT)
    end
    (E, S) = sum_mf3_partials(lines, mts, source)
    return (E, S, [length(E)], [2])
end

"""
    warn_energy_range(db::NuclearReactionENDFDB, E_max::Float64, db_name::String,
    particle::Particle)

Warns for every isotope whose tabulated cross-section stops below the highest energy of the
simulation. Above the last tabulated energy the cross-section is held constant, as the TAB1
interpolation clamps to the last value, so the transport is then carried out on an
extrapolated cross-section.

# Input Argument(s)
- `db::NuclearReactionENDFDB` : nuclear-reaction ENDF database of one incident particle.
- `E_max::Float64` : highest energy of the group structure of that particle [MeV], as given
  by the user and stored in the group structure, before the conversion to mₑc² performed
  in multigroup().
- `db_name::String` : ENDF database directory name.
- `particle::Particle` : incident particle.

# Output Argument(s)
N/A
"""
function warn_energy_range(db::NuclearReactionENDFDB, E_max::Float64, db_name::String, particle::Particle)
    E_max_eV = E_max * 1.0e6 # MeV -> eV (ENDF native unit)
    for (Z, A) in sort!(collect(keys(db.isotopes)))
        iso = db.isotopes[(Z, A)]
        isempty(iso.E) && continue
        if iso.E[end] < E_max_eV
            @warn "The nuclear-reaction cross-section of $(atomic_symbol(Z))-$(A) (Z=$(Z), A=$(A)) of library $(db_name) is tabulated up to $(round(iso.E[end]/1.0e6, digits=3)) MeV only, below the highest energy of the $(get_type(particle)) group structure ($(round(E_max_eV/1.0e6, digits=3)) MeV): it is held constant above the last tabulated energy."
        end
    end
    return nothing
end

"""
    nuclear_reaction_endf(db_name::String, isotopes::Vector{Tuple{Int,Int}},
    particle::Particle, mt::Int; data_root::Union{Nothing,String}=nothing)

Builds the nuclear-reaction ENDF database for a set of target isotopes and one incident
particle, reading the MF=3 total cross-section curve for the requested reaction (MT)
from each isotope's ENDF file. When the non-elastic total (MT=3) is requested but not
tabulated, it is rebuilt from the partial channels, or set to zero when the evaluation
carries no non-elastic data; both cases are reported for the isotope concerned.

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
        source = "$(atomic_symbol(Z))-$(A) (Z=$(Z), A=$(A)) of library $(db_name)"
        lines = readlines(endf_path)
        (E, S, NBT, INT) = mf3_nonelastic_table(lines, mt, source)
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

Reads MF=3 cross-sections and MF=6 Kalbach-Mann energy distributions for the
charged-particle production channels of each isotope's ENDF file. Only LAW=1 LANG=2
(Kalbach-Mann) product subsections are retained; other distributions are silently skipped.

When `production_mts` is empty, the channels are chosen per evaluation as the non-elastic
reactions present in both MF=3 and MF=6, passed through the ENDF sum rules of
`nonelastic_partial_mts` so that a reaction already contained in a coarser one present in
the same file is not counted twice. This picks up MT=5, which carries the whole production
of the light elements and everything above the closing of the exclusive channels, around
30 MeV, for the heavier ones. MT=5 is by definition disjoint from the reactions given
explicitly, so summing it with them cannot double count.

# Input Argument(s)
- `db_name::String` : ENDF database directory name.
- `isotopes::Vector{Tuple{Int,Int}}` : list of `(Z, A)` target isotopes.
- `particle::Particle` : incident particle type.
- `production_mts::Vector{Int}` : ENDF MT numbers to search for charged-particle products;
  empty to let each evaluation decide, as described above.
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
        i_mf3 = find_first_mf(lines, 3)
        i_mf6 = find_first_mf(lines, 6)
        if isnothing(i_mf3) || isnothing(i_mf6)
            db[(Z, A)] = channels
            continue
        end

        # Channels of this evaluation: the requested ones, or the non-redundant non-elastic
        # reactions it tabulates in both MF=3 and MF=6.
        mts = if isempty(production_mts)
            nonelastic_partial_mts(intersect(list_mf_mts(lines, 6), list_mf_mts(lines, 3)))
        else
            production_mts
        end

        for mt in mts
            find_section_start(lines, 3, mt; start_i=i_mf3) === nothing && continue
            find_section_start(lines, 6, mt; start_i=i_mf6) === nothing && continue

            # Read MF=3 total XS for this channel
            (_, _, E_xs, S_xs, NBT_xs, INT_xs) = read_mf3_mt(lines, mt)

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
                    prod.law1_E_out, prod.law1_f, prod.law1_r,
                    isnothing(prod.law1_LEP) ? 2 : prod.law1_LEP)

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
