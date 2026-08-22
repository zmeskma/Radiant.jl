"""
    IsotopeNuclearReactionDB

Stores the ENDF MF=3 total cross-section curve for one target isotope and one reaction
(MT number).

# Field(s)
- `E::Vector{Float64}` : incident energy grid [eV].
- `S::Vector{Float64}` : tabulated cross-section [barn].
- `NBT::Vector{Int}` : interpolation-region breakpoints.
- `INT::Vector{Int}` : interpolation laws.
"""
struct IsotopeNuclearReactionDB
    E::Vector{Float64}
    S::Vector{Float64}
    NBT::Vector{Int}
    INT::Vector{Int}
end

"""
    NuclearReactionENDFDB

Stores the nuclear-reaction ENDF database for one incoming particle type.

# Field(s)
- `mt::Int` : ENDF reaction number (MT) this database was built for.
- `isotopes::Dict{Tuple{Int,Int},IsotopeNuclearReactionDB}` : isotope data indexed by
  `(Z, A)`.
"""
struct NuclearReactionENDFDB
    mt::Int
    isotopes::Dict{Tuple{Int,Int},IsotopeNuclearReactionDB}
end

"""
    Nuclear_Reaction

Structure used to define parameters for production of multigroup nuclear-reaction total
cross-sections, read directly from ENDF MF=3 data. Only the total cross-section is
produced — no scattering matrix, secondary-particle production, stopping power or
momentum transfer.

# Optional field(s) - with default values
- `interaction_types::Dict{Tuple{Type,Type},Vector{String}} = Dict((Proton,Proton) => ["A"])` : absorption-only registration.
- `library::String = "TENDL2023"` : nuclear-reaction ENDF library name.
- `data_path::String = "../../data"` : root directory for the ENDF data.
- `mt::Int64 = 3` : ENDF reaction number (MT) to read from MF=3 (default: non-elastic total).
"""
mutable struct Nuclear_Reaction <: Interaction

    # Variable(s)
    name::String
    incoming_particle::Vector{Type}
    interaction_particles::Vector{Type}
    interaction_types::Dict{Tuple{Type,Type},Vector{String}}
    is_CSD::Bool
    is_AFP::Bool
    is_AFP_decomposition::Bool
    is_elastic::Bool
    is_subshells_dependant::Bool
    scattering_model::String
    library::String
    data_path::String
    mt::Int64
    endf_db::Dict{Type,NuclearReactionENDFDB}

    # Constructor(s)
    function Nuclear_Reaction()
        this = new()
        this.name = "Nuclear_Reaction"
        this.interaction_types = Dict((Proton,Proton) => ["A"])
        this.incoming_particle = unique([t[1] for t in collect(keys(this.interaction_types))])
        this.interaction_particles = unique([t[2] for t in collect(keys(this.interaction_types))])
        this.is_CSD = false
        this.is_AFP = false
        this.is_AFP_decomposition = false
        this.is_elastic = false
        this.is_subshells_dependant = false
        this.scattering_model = "BTE"
        this.set_library("TENDL2023")
        this.set_data_path("../../data")
        this.set_mt(3)
        this.set_endf_db(Dict{Type,NuclearReactionENDFDB}())
        return this
    end
end

# Method(s)
"""
    set_interaction_types(this::Nuclear_Reaction,interaction_types::Dict{Tuple{DataType,DataType},Vector{String}})

To define the interaction types for nuclear-reaction processes.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `interaction_types::Dict{Tuple{DataType,DataType},Vector{String}}` : Dictionary of the interaction processes types, of the form (incident particle,outgoing particle) => associated list of interaction type, which can be:
    - `(Proton,Proton) => ["A"]` : absorption of incoming protons via nuclear reactions.

# Output Argument(s)
N/A

# Examples
```jldoctest
julia> nuclear_reaction = Nuclear_Reaction()
julia> nuclear_reaction.set_interaction_types( Dict((Proton,Proton) => ["A"]) )
```
"""
function set_interaction_types(this::Nuclear_Reaction,interaction_types)
    this.interaction_types = interaction_types
end

"""
    set_library(this::Nuclear_Reaction,library::String)

Define the nuclear-reaction database library name (e.g., TENDL2023).

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `library::String` : library name.

# Output Argument(s)
N/A
"""
function set_library(this::Nuclear_Reaction,library::String)
    this.library = library
end

"""
    get_library(this::Nuclear_Reaction)

Get the nuclear-reaction database library name.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `library::String` : library name.
"""
function get_library(this::Nuclear_Reaction)
    return this.library
end

"""
    set_data_path(this::Nuclear_Reaction,data_path::String)

Define the root directory containing the nuclear-reaction ENDF data files. Relative
paths are resolved from the location of `Nuclear_Reaction.jl`.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `data_path::String` : data root path. The default `../../data` resolves to the
  package data directory.

# Output Argument(s)
N/A
"""
function set_data_path(this::Nuclear_Reaction,data_path::String)
    path = isabspath(data_path) ? normpath(data_path) : normpath(joinpath(@__DIR__, data_path))
    this.data_path = path
    this.endf_db = Dict{Type,NuclearReactionENDFDB}()
end

"""
    get_data_path(this::Nuclear_Reaction)

Get the root directory containing the nuclear-reaction ENDF data files.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `data_path::String` : data root path.
"""
function get_data_path(this::Nuclear_Reaction)
    return this.data_path
end

"""
    set_mt(this::Nuclear_Reaction,mt::Int64)

Define the ENDF reaction number (MT) to read from MF=3.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `mt::Int64` : ENDF MT number (e.g. 3 for non-elastic total).

# Output Argument(s)
N/A
"""
function set_mt(this::Nuclear_Reaction,mt::Int64)
    this.mt = mt
    this.endf_db = Dict{Type,NuclearReactionENDFDB}()
end

"""
    get_mt(this::Nuclear_Reaction)

Get the ENDF reaction number (MT) read from MF=3.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `mt::Int64` : ENDF MT number.
"""
function get_mt(this::Nuclear_Reaction)
    return this.mt
end

"""
    set_endf_db(this::Nuclear_Reaction,endf_db::Dict{Type,NuclearReactionENDFDB})

Set the nuclear-reaction ENDF database.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `endf_db::Dict{Type,NuclearReactionENDFDB}` : per-particle ENDF database.

# Output Argument(s)
N/A
"""
function set_endf_db(this::Nuclear_Reaction,endf_db::Dict{Type,NuclearReactionENDFDB})
    this.endf_db = endf_db
end

"""
    get_endf_db(this::Nuclear_Reaction)

Get the nuclear-reaction ENDF database.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `endf_db::Dict{Type,NuclearReactionENDFDB}` : per-particle ENDF database.
"""
function get_endf_db(this::Nuclear_Reaction)
    return this.endf_db
end

"""
    initialize(this::Nuclear_Reaction, particles, isotopes::Vector{Tuple{Int,Int}};
    energy_boundaries::Union{Nothing,Vector{Vector{Float64}}}=nothing)

Initializes the nuclear-reaction ENDF databases required by the incoming particles.
Already initialized particle databases are reused. When the energy group structure is
provided, the tabulated energy range of each isotope is compared with the highest energy
of the corresponding particle and a warning is issued for the evaluations that stop below
it.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `particles` : particle objects used in the simulation.
- `isotopes::Vector{Tuple{Int,Int}}` : target isotope pairs `(Z, A)` required by the materials.
- `energy_boundaries::Union{Nothing,Vector{Vector{Float64}}}` : energy group boundaries of
  each particle [MeV], in the same order as `particles`.

# Output Argument(s)
N/A
"""
function initialize(this::Nuclear_Reaction, particles, isotopes::Vector{Tuple{Int,Int}};
    energy_boundaries::Union{Nothing,Vector{Vector{Float64}}}=nothing)
    endf_db = this.get_endf_db()
    for (n, particle) in enumerate(particles)
        ptype = get_type(particle)
        if ptype ∈ this.get_in_particles()
            if !haskey(endf_db, ptype)
                endf_db[ptype] = nuclear_reaction_endf(
                    this.get_library(), isotopes, particle, this.get_mt();
                    data_root=this.get_data_path())
            end
            if !isnothing(energy_boundaries)
                warn_energy_range(endf_db[ptype], maximum(energy_boundaries[n]),
                    this.get_library(), particle)
            end
        end
    end
    this.set_endf_db(endf_db)
    return nothing
end

"""
    in_distribution(this::Nuclear_Reaction)

Describe the energy discretization method for the incoming particle in the nuclear
reaction interaction.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.

# Output Argument(s)
- `is_dirac::Bool` : boolean describing if a Dirac distribution is used.
- `N::Int64` : number of quadrature points.
- `quadrature::String` : type of quadrature.
"""
function in_distribution(this::Nuclear_Reaction)
    is_dirac = false
    N = 8
    quadrature = "gauss-legendre"
    return is_dirac, N, quadrature
end

"""
    tcs(this::Nuclear_Reaction,Ei::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

Gives the total cross-section for the nuclear-reaction interaction.

# Input Argument(s)
- `this::Nuclear_Reaction` : nuclear reaction structure.
- `Ei::Float64` : incoming particle energy.
- `particle::Particle` : incoming particle.
- `Z::Int64` : atomic number.
- `A::Union{Nothing,Vector{Int64}}` : isotope mass numbers for this element.
- `atpercentA::Union{Nothing,Vector{Float64}}` : isotope atomic fractions.

# Output Argument(s)
- `σt::Float64` : total cross-section.
"""
function tcs(this::Nuclear_Reaction,Ei::Float64,particle::Particle,Z::Int64;
    A::Union{Nothing,Vector{Int64}}=nothing,
    atpercentA::Union{Nothing,Vector{Float64}}=nothing)

    if isnothing(A) || isnothing(atpercentA)
        pairs = isotopic_composition(Z)
        A = [p[1] for p in pairs]
        atpercentA = [p[2] for p in pairs]
    end

    ptype = get_type(particle)
    db = this.get_endf_db()[ptype]
    E_eV = Ei * M_E_C2_MEV * 1.0e6 # mₑc² -> MeV -> eV (ENDF native unit)

    σt = 0.0
    for (Ai, atai) in zip(A, atpercentA)
        iso = db.isotopes[(Z,Ai)]
        σt += atai * interp_TAB1(E_eV, iso.E, iso.S, iso.NBT, iso.INT)
    end

    return σt * BARN_TO_CM2
end
