#=
    plot_nuclear_reaction_cross_sections.jl

Test script for the "A" (absorption) interaction type of Nuclear_Reaction.

For one isotope of each element, it reads the absolute non-elastic (p,x) microscopic
cross-section from an ENDF proton library and plots it. This is the very cross-section the
"A" type feeds into Sigma_t and Sigma_a: the script calls mf3_nonelastic_table(), the same
Radiant function the interaction itself calls, so the curves show exactly what a transport
calculation would see, fallbacks included:

- MF=3 MT=3 tabulated          -> the evaluated non-elastic total,
- MT=3 absent, partials found  -> the sum of the non-elastic partial channels,
- no usable non-elastic data   -> zero.

Above the last tabulated energy the TAB1 interpolation clamps to the last value, so each
plot continues as a dashed line up to E_PLOT_MAX to show the constant extrapolation the
transport would actually use.

Outputs, written into a subdirectory of the library directory:
- per_element/Z###_Sym-A.png : one plot per element,
- overview_map.png           : log10(sigma) over all elements and energies,
- summary.csv                : isotope, provenance, energy range and a few sample values.

Usage:
    julia --project=. scripts/plot_nuclear_reaction_cross_sections.jl [data_root]

`data_root` is the directory holding ENDF/<library>/. It defaults to the RADIANT_ENDF_DATA
environment variable, then to the data/ directory of this repository. Plotting uses PyPlot
from the default Julia environment; Radiant itself needs no plotting dependency.
=#

using Radiant
using PyPlot
using Printf

# ---- User parameters ---------------------------------------------------------------

const LIBRARY     = "TENDL2023"
const OUT_SUBDIR  = "cross_section_plots"  # subdirectory created inside the library folder
const Z_RANGE     = 1:115                  # elements to process
const MT_TOTAL    = 3                      # ENDF MT of the non-elastic total
const E_PLOT_MIN  = 1.0e-1                 # MeV, lower end of the plots
const E_PLOT_MAX  = 3.0e2                  # MeV, upper end of the plots
const E_SAMPLES   = [10.0, 50.0, 100.0, 200.0]   # MeV, sample values written to summary.csv
const N_MAP       = 240                    # energy points of the overview map
const MAP_LOG_MIN = -4.0                   # log10(barn) floor of the overview colour scale
const MAP_LOG_MAX = 0.5                    # log10(barn) ceiling of the overview colour scale

# ---- Library inspection ------------------------------------------------------------

"""
    available_isotopes(library_dir::String, Z::Int)

Lists the mass numbers of the isotopes of element `Z` present in a proton ENDF library
directory, from the `p_ZZZ-Sym-AAA.endf` file names.
"""
function available_isotopes(library_dir::String, Z::Int)
    dir = joinpath(library_dir, "proton")
    isdir(dir) || return Int[]
    pattern = Regex("^p_" * @sprintf("%03d", Z) * raw"-[A-Za-z]+-(\d{3})\.endf$")
    As = Int[]
    for f in readdir(dir)
        m = match(pattern, f)
        isnothing(m) || push!(As, parse(Int, m.captures[1]))
    end
    return sort!(As)
end

"""
    select_isotope(Z::Int, As::Vector{Int})

Chooses the isotope of element `Z` to be plotted: the most abundant natural isotope when
the element has a natural composition and that isotope is present in the library, and
otherwise the median of the mass numbers available, as a representative evaluation.
"""
function select_isotope(Z::Int, As::Vector{Int})
    isempty(As) && return nothing
    pairs = try
        Radiant.isotopic_composition(Z)
    catch
        nothing # no natural composition: synthetic element
    end
    if !isnothing(pairs) && !isempty(pairs)
        A_nat = first(pairs[argmax([p[2] for p in pairs])])
        A_nat in As && return A_nat
    end
    return As[cld(length(As), 2)]
end

# ---- Cross-section retrieval -------------------------------------------------------

"""
    nonelastic_cross_section(data_root::String, Z::Int, A::Int, particle::Particle)

Reads the non-elastic (p,x) cross-section of one isotope through the same Radiant call the
"A" interaction type uses, and reports how the table was obtained.

Returns a named tuple with the energy grid `E` [MeV], the cross-section `S` [barn], the
TAB1 metadata `NBT` and `INT`, the provenance `origin` (`:tabulated`, `:rebuilt` or
`:zero`) and the channels `mts` used when the total had to be rebuilt.
"""
function nonelastic_cross_section(data_root::String, Z::Int, A::Int, particle::Particle)
    path = Radiant.endf_path_for_isotope(LIBRARY, Z, A, particle; data_root=data_root)
    isfile(path) || return nothing
    lines = readlines(path)

    # Provenance, established with the same selection rules as the reader itself.
    i_mf3 = Radiant.find_first_mf(lines, 3)
    origin, mts = if isnothing(i_mf3)
        (:zero, Int[])
    elseif !isnothing(Radiant.find_section_start(lines, 3, MT_TOTAL; start_i=i_mf3))
        (:tabulated, Int[])
    else
        partials = Radiant.nonelastic_partial_mts(Radiant.list_mf_mts(lines, 3))
        (isempty(partials) ? :zero : :rebuilt, partials)
    end

    label = "$(Radiant.atomic_symbol(Z))-$(A) (Z=$(Z), A=$(A)) of library $(LIBRARY)"
    (E_eV, S, NBT, INT) = Radiant.mf3_nonelastic_table(lines, MT_TOTAL, label)

    return (E=E_eV ./ 1.0e6, S=S, NBT=NBT, INT=INT, origin=origin, mts=mts)
end

"""
    sigma_at(xs, E_MeV::Float64)

Interpolates a cross-section table at one energy [MeV] exactly as tcs() does, clamping to
the last tabulated value above the end of the grid.
"""
function sigma_at(xs, E_MeV::Float64)
    return Radiant.interp_TAB1(E_MeV * 1.0e6, xs.E .* 1.0e6, xs.S, xs.NBT, xs.INT)
end

# ---- Plotting ----------------------------------------------------------------------

function origin_caption(xs)
    xs.origin === :tabulated && return "MF=3 MT=3 tabulated"
    xs.origin === :rebuilt   && return "rebuilt from MF=3 MT=" * string(xs.mts)
    return "no non-elastic data: cross-section set to zero"
end

"""
    plot_isotope(dir::String, Z::Int, A::Int, xs)

Draws the non-elastic cross-section of one isotope and saves it as a PNG. The tabulated
range is drawn solid; beyond the last tabulated energy the constant extrapolation used by
the transport is drawn dashed.
"""
function plot_isotope(dir::String, Z::Int, A::Int, xs)
    symbol = Radiant.atomic_symbol(Z)
    fig = figure(figsize=(6.4, 4.4))
    ax = fig.add_subplot(1, 1, 1)

    if all(iszero, xs.S)
        ax.set_xscale("log")
        ax.set_xlim(E_PLOT_MIN, E_PLOT_MAX)
        ax.set_ylim(0.0, 1.0)
        ax.plot([E_PLOT_MIN, E_PLOT_MAX], [0.0, 0.0], color="tab:red", linewidth=1.8)
        ax.text(0.5, 0.5, "no non-elastic data\nsigma = 0", transform=ax.transAxes,
                ha="center", va="center", fontsize=11, color="tab:red")
    else
        keep = xs.S .> 0.0
        ax.loglog(xs.E[keep], xs.S[keep], color="tab:blue", linewidth=1.8,
                  label="tabulated")
        E_end = xs.E[end]
        if E_end < E_PLOT_MAX
            ax.loglog([E_end, E_PLOT_MAX], [xs.S[end], xs.S[end]], color="tab:orange",
                      linewidth=1.6, linestyle="--", label="clamped extrapolation")
            ax.axvline(E_end, color="tab:orange", linewidth=0.9, linestyle=":")
        end
        ax.set_xlim(E_PLOT_MIN, E_PLOT_MAX)
        ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    end

    ax.set_xlabel("proton energy [MeV]")
    ax.set_ylabel(L"$\sigma_{\mathrm{non-elastic}}$ [barn]")
    ax.set_title("$(symbol)-$(A)  (Z=$(Z))  non-elastic (p,x) cross-section\n" *
                 origin_caption(xs), fontsize=10)
    ax.grid(true, which="both", linewidth=0.4, alpha=0.5)
    fig.tight_layout()

    name = @sprintf("Z%03d_%s-%03d.png", Z, symbol, A)
    fig.savefig(joinpath(dir, name), dpi=140)
    close(fig)
    return name
end

"""
    plot_overview(path::String, results)

Draws a map of log10(sigma) over all the elements processed and the plotted energy range,
to expose at a glance the evaluations that are zero or that stop early.
"""
function plot_overview(path::String, results)
    E = 10 .^ collect(range(log10(E_PLOT_MIN), log10(E_PLOT_MAX), length=N_MAP))
    Zs = [r.Z for r in results]
    M = fill(NaN, length(results), N_MAP)
    for (i, r) in enumerate(results), (j, Ej) in enumerate(E)
        s = sigma_at(r.xs, Ej)
        M[i, j] = s > 0.0 ? log10(s) : NaN
    end

    fig = figure(figsize=(9.0, 10.0))
    ax = fig.add_subplot(1, 1, 1)
    mesh = ax.pcolormesh(E, Zs, M, shading="auto", cmap="viridis",
                         vmin=MAP_LOG_MIN, vmax=MAP_LOG_MAX)
    ax.set_xscale("log")
    ax.set_xlabel("proton energy [MeV]")
    ax.set_ylabel("atomic number Z")
    ax.set_title("Non-elastic (p,x) cross-section used by Nuclear_Reaction type \"A\"\n" *
                 "$(LIBRARY), one isotope per element (white: zero or no data)", fontsize=11)
    cb = fig.colorbar(mesh, ax=ax, extend="min")
    cb.set_label(L"$\log_{10}(\sigma\ /\ \mathrm{barn})$")

    # Mark the end of each tabulated range.
    E_ends = [min(r.xs.E[end], E_PLOT_MAX) for r in results]
    ax.plot(E_ends, Zs, color="white", linewidth=0.8, linestyle="--", alpha=0.8,
            label="last tabulated energy")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    close(fig)
end

# ---- Main --------------------------------------------------------------------------

function main()
    data_root = if !isempty(ARGS)
        ARGS[1]
    else
        get(ENV, "RADIANT_ENDF_DATA", normpath(joinpath(@__DIR__, "..", "data")))
    end
    library_dir = joinpath(data_root, "ENDF", LIBRARY)
    isdir(library_dir) || error("Library directory not found: $(library_dir). Pass the " *
                                "data root as the first argument or set RADIANT_ENDF_DATA.")

    out_dir = joinpath(library_dir, OUT_SUBDIR)
    per_element_dir = joinpath(out_dir, "per_element")
    mkpath(per_element_dir)
    println("Library    : ", library_dir)
    println("Output     : ", out_dir)

    pygui(false)
    particle = Proton()
    results = NamedTuple[]

    for Z in Z_RANGE
        As = available_isotopes(library_dir, Z)
        A = select_isotope(Z, As)
        if isnothing(A)
            @warn "No evaluation available for Z=$(Z) in $(LIBRARY): element skipped."
            continue
        end
        xs = nonelastic_cross_section(data_root, Z, A, particle)
        isnothing(xs) && continue
        plot_isotope(per_element_dir, Z, A, xs)
        push!(results, (Z=Z, A=A, xs=xs))
    end

    plot_overview(joinpath(out_dir, "overview_map.png"), results)

    open(joinpath(out_dir, "summary.csv"), "w") do io
        println(io, "Z,symbol,A,origin,channels,E_min_MeV,E_max_MeV," *
                    join(["sigma_" * string(E) * "MeV_barn" for E in E_SAMPLES], ","))
        for r in results
            channels = isempty(r.xs.mts) ? "" : join(r.xs.mts, " ")
            samples = join([@sprintf("%.6g", sigma_at(r.xs, E)) for E in E_SAMPLES], ",")
            @printf(io, "%d,%s,%d,%s,%s,%.6g,%.6g,%s\n", r.Z, Radiant.atomic_symbol(r.Z),
                    r.A, string(r.xs.origin), channels, r.xs.E[1], r.xs.E[end], samples)
        end
    end

    n_tab = count(r -> r.xs.origin === :tabulated, results)
    n_reb = count(r -> r.xs.origin === :rebuilt, results)
    n_zero = count(r -> r.xs.origin === :zero, results)
    n_short = count(r -> r.xs.E[end] < E_PLOT_MAX, results)
    println()
    println("Elements plotted        : ", length(results))
    println("  MF=3 MT=3 tabulated   : ", n_tab)
    println("  rebuilt from partials : ", n_reb)
    println("  set to zero           : ", n_zero)
    println("  ending below ", E_PLOT_MAX, " MeV : ", n_short)

    # End-to-end check: the plotted table is the one the interaction serves to transport.
    r = results[argmax([r.Z for r in results])]
    db = Radiant.nuclear_reaction_endf(LIBRARY, [(r.Z, r.A)], particle, MT_TOTAL;
                                       data_root=data_root)
    iso = db.isotopes[(r.Z, r.A)]
    ok = iso.E ≈ r.xs.E .* 1.0e6 && iso.S ≈ r.xs.S
    println("Check on $(Radiant.atomic_symbol(r.Z))-$(r.A): plotted table identical to " *
            "the one nuclear_reaction_endf() builds: ", ok)

    return nothing
end

main()
