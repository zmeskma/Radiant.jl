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
- `E_max::Float64` : highest energy of the group structure of that particle [mₑc²].
- `db_name::String` : ENDF database directory name.
- `particle::Particle` : incident particle.

# Output Argument(s)
N/A
"""
function warn_energy_range(db::NuclearReactionENDFDB, E_max::Float64, db_name::String, particle::Particle)
    E_max_eV = E_max * M_E_C2_MEV * 1.0e6 # mₑc² -> MeV -> eV (ENDF native unit)
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
