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
    _NEUTRAL_ZAPS

ENDF product identifiers of the neutral particles a nuclear reaction emits, the photon and
the neutron. No interaction type transports them; their channels are read so that the
energy they carry can be removed from the local deposition when it is taken to leave the
control volume.
"""
const _NEUTRAL_ZAPS = Set{Int}([0, 1])

"""
    nuclear_production_endf(db_name::String, isotopes::Vector{Tuple{Int,Int}},
    particle::Particle, production_mts::Vector{Int};
    data_root::Union{Nothing,String}=nothing)

Reads MF=3 cross-sections and MF=6 energy distributions for the production channels of
each isotope's ENDF file, for the charged particles the equivalent-proton method
transports and for the neutrons and photons whose energy may have to be removed from the
local deposition. LAW=1 subsections are retained whatever their representation, the
outgoing spectrum occupying the same slots for every LANG; other laws are skipped.

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
                (ZAP_i ∈ _CHARGED_ZAPS || ZAP_i ∈ _NEUTRAL_ZAPS) || continue
                prod.law1_LANG === nothing && continue
                isempty(prod.law1_E_in)    && continue
                isempty(prod.law1_E_out)   && continue

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
    equivalent_proton_spectra(interaction::Nuclear_Reaction, channels, E_in::Float64,
    Z::Vector{Int64}, ωz::Vector{Float64}, ρ::Float64, state_of_matter::String,
    I_eff::Float64)

Prepares, for one incident energy, everything the feed needs from the production channels
of an isotope: the channel cross-section and yield at that energy, the outgoing spectrum
interpolated there, and the energies of its points carried over to the equivalent-proton
scale by CSDA range matching.

The transformation of the energy grid is the same for every incident energy, so it rests
on the range table cached per product and material by `get_or_build_eq_cache!`. Doing it
here, once per channel, leaves the feed with nothing to do but integrate over its groups.

# Input Argument(s)
- `interaction::Nuclear_Reaction` : nuclear reaction structure, holding the range cache.
- `channels` : production channels of the isotope.
- `E_in::Float64` : incident energy [eV].
- `Z::Vector{Int64}` : atomic numbers of the elements composing the material.
- `ωz::Vector{Float64}` : weight fraction of each element.
- `ρ::Float64` : material density [g/cm³].
- `state_of_matter::String` : material state.
- `I_eff::Float64` : effective mean excitation energy [mₑc²]; NaN to use the tables.

# Output Argument(s)
- `spectra::Vector{NamedTuple}` : one entry per channel, with the product energies `E_b`
  [eV] and their spectrum `f` [eV⁻¹], the interpolation law `LEP`, the same energies on the
  equivalent-proton scale `E_eq` [mₑc²], the product of cross-section and yield `σy`
  [cm²], and the range table `E_b_tbl`, `E_eq_tbl` [mₑc²] used to invert the mapping.

# Reference(s)
- Salvat & Quesada (2020), NIMB 475, 49–62, §4.2 equivalent-proton method.
"""
function equivalent_proton_spectra(interaction::Nuclear_Reaction, channels, E_in::Float64,
        Z::Vector{Int64}, ωz::Vector{Float64}, ρ::Float64, state_of_matter::String,
        I_eff::Float64)

    eV_per_mc2 = M_E_C2_MEV * 1e6
    spectra = NamedTuple[]

    for channel in channels
        # Only the charged products are transported; the neutrons and the photons are
        # handled apart, by neutral_energy_release.
        channel.ZAP ∈ _CHARGED_ZAPS || continue
        channel.kalbach === nothing && continue

        σ_ch = interp_TAB1(E_in, channel.E_xs, channel.S_xs, channel.NBT_xs, channel.INT_xs)
        σ_ch <= 0.0 && continue
        y_ch = isempty(channel.yield_E) ? 1.0 :
               interp_TAB1(E_in, channel.yield_E, channel.yield_y, Int[], Int[])
        y_ch <= 0.0 && continue

        E_b, f = interpolate_spectrum(channel.kalbach, E_in)
        length(E_b) < 2 && continue

        (E_b_tbl, E_eq_tbl, _) = get_or_build_eq_cache!(
            interaction, channel.ZAP, Z, ωz, ρ, state_of_matter, I_eff)

        # The mapping is monotone, so the points keep their order on the new scale.
        E_eq = [interp_linear(E_b_tbl, E_eq_tbl, E / eV_per_mc2) for E in E_b]

        # The feed is built in cm², as the other feed functions are.
        push!(spectra, (E_b=E_b, f=f, LEP=channel.kalbach.LEP, E_eq=E_eq,
                        σy=σ_ch * BARN_TO_CM2 * y_ch,
                        E_b_tbl=E_b_tbl, E_eq_tbl=E_eq_tbl))
    end

    return spectra
end

"""
    integrate_equivalent_proton_group(spectra, Ef⁺::Float64, Ef⁻::Float64)

Integrates the equivalent-proton production of one outgoing energy group, that is over
`Ef⁺ ≤ E_eq ≤ Ef⁻` [mₑc²], over the prepared spectra of an isotope.

Only the part of a spectrum interval that falls inside the group is taken, and it is
integrated between the product energies bounding that part, so that a group narrower than
the tabulation still receives its share and one wider than it receives the whole. The
spectrum being linear over an interval, the number of products follows the trapezoidal
rule and the energy they carry the closed form of the quadratic integrand; for a histogram
spectrum both reduce to the exact rectangle rules.

# Input Argument(s)
- `spectra` : prepared spectra, from `equivalent_proton_spectra`.
- `Ef⁺::Float64` : lower boundary of the outgoing group [mₑc²].
- `Ef⁻::Float64` : upper boundary of the outgoing group [mₑc²].

# Output Argument(s)
- `𝓕i::Float64` : equivalent protons produced in the group [cm²].
- `𝓕iₑ::Float64` : energy they carry [cm² × mₑc²].

# Reference(s)
- Salvat & Quesada (2020), NIMB 475, 49–62, §4.2 equivalent-proton method.
"""
function integrate_equivalent_proton_group(spectra, Ef⁺::Float64, Ef⁻::Float64)

    eV_per_mc2 = M_E_C2_MEV * 1e6
    𝓕i = 0.0
    𝓕iₑ = 0.0
    Ef⁻ <= Ef⁺ && return 𝓕i, 𝓕iₑ

    for s in spectra
        # Intervals of this spectrum reaching into the group.
        n = length(s.E_eq)
        (s.E_eq[1] >= Ef⁻ || s.E_eq[n] <= Ef⁺) && continue
        j_first = clamp(searchsortedlast(s.E_eq, Ef⁺), 1, n - 1)
        j_last  = clamp(searchsortedlast(s.E_eq, Ef⁻), 1, n - 1)

        for j in j_first:j_last
            E_eq_lo, E_eq_hi = s.E_eq[j], s.E_eq[j+1]
            E_eq_hi <= E_eq_lo && continue
            E_b_lo, E_b_hi = s.E_b[j], s.E_b[j+1]
            E_b_hi <= E_b_lo && continue

            f_lo = s.f[j]
            # A histogram bin holds its own value across the interval; a linear-linear one
            # rises to the next point, and must be kept even when it starts from zero.
            f_up = s.LEP == 1 ? f_lo : s.f[j+1]
            (f_lo <= 0.0 && f_up <= 0.0) && continue

            # Part of the interval inside the group, and the product energies bounding it.
            a = max(E_eq_lo, Ef⁺)
            b = min(E_eq_hi, Ef⁻)
            b <= a && continue
            E_b_a = a <= E_eq_lo ? E_b_lo :
                clamp(interp_linear(s.E_eq_tbl, s.E_b_tbl, a) * eV_per_mc2, E_b_lo, E_b_hi)
            E_b_b = b >= E_eq_hi ? E_b_hi :
                clamp(interp_linear(s.E_eq_tbl, s.E_b_tbl, b) * eV_per_mc2, E_b_lo, E_b_hi)
            ΔE_b = E_b_b - E_b_a
            ΔE_b <= 0.0 && continue

            u_a = (E_b_a - E_b_lo) / (E_b_hi - E_b_lo)
            u_b = (E_b_b - E_b_lo) / (E_b_hi - E_b_lo)
            f_a = f_lo + u_a * (f_up - f_lo)
            f_b = f_lo + u_b * (f_up - f_lo)
            df = f_b - f_a

            mass   = 0.5 * (f_a + f_b) * ΔE_b
            mass_E = ΔE_b * (f_a * E_b_a + (f_a * ΔE_b + df * E_b_a) / 2 +
                             df * ΔE_b / 3) / eV_per_mc2

            # w × E_eq = E_b, so the equivalent protons carry the energy of the product.
            w = (0.5 * (E_b_a + E_b_b) / eV_per_mc2) / (0.5 * (a + b))

            𝓕i  += s.σy * mass * w
            𝓕iₑ += s.σy * mass_E
        end
    end

    return 𝓕i, 𝓕iₑ
end

"""
    neutral_energy_release(channels, E_in::Float64)

Gives the energy the neutrons and photons of a reaction carry away from one isotope, at
one incident energy, as the sum over their production channels of the channel
cross-section times the product yield times the mean outgoing energy of its spectrum.

Nothing transports these particles. The quantity is meant to be removed from the local
energy deposition when they are taken to leave the control volume, which is what the
`is_neutral_escape` option of `Nuclear_Reaction` asks for; left as it is, their energy
stays where the reaction happened.

Only the channels whose distribution the reader could keep contribute, so a product given
under a law other than LAW=1 is missing from the balance.

# Input Argument(s)
- `channels` : production channels of the isotope.
- `E_in::Float64` : incident energy [eV].

# Output Argument(s)
- `energy::Float64` : energy carried away by the neutrons and photons [cm² × mₑc²].

# Reference(s)
- ENDF-6 Formats Manual, MF=6 product energy distributions.
"""
function neutral_energy_release(channels, E_in::Float64)

    eV_per_mc2 = M_E_C2_MEV * 1e6
    energy = 0.0

    for channel in channels
        channel.ZAP ∈ _NEUTRAL_ZAPS || continue
        channel.kalbach === nothing && continue

        σ_ch = interp_TAB1(E_in, channel.E_xs, channel.S_xs, channel.NBT_xs, channel.INT_xs)
        σ_ch <= 0.0 && continue
        y_ch = isempty(channel.yield_E) ? 1.0 :
               interp_TAB1(E_in, channel.yield_E, channel.yield_y, Int[], Int[])
        y_ch <= 0.0 && continue

        E_out, f = interpolate_spectrum(channel.kalbach, E_in)
        length(E_out) < 2 && continue

        # Mean outgoing energy, ∫f E dE over the tabulated spectrum and its own law.
        mean_E = 0.0
        for j in 1:length(E_out)-1
            ΔE = E_out[j+1] - E_out[j]
            ΔE <= 0.0 && continue
            f_a = f[j]
            f_b = channel.kalbach.LEP == 1 ? f[j] : f[j+1]
            df  = f_b - f_a
            mean_E += ΔE * (f_a * E_out[j] + (f_a * ΔE + df * E_out[j]) / 2 +
                            df * ΔE / 3)
        end

        energy += σ_ch * BARN_TO_CM2 * y_ch * mean_E / eV_per_mc2
    end

    return energy
end
