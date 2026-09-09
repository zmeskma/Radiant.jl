# Reference(s)
# - Fritsch, F. N. and Carlson, R. E. (1980), Monotone piecewise cubic
#   interpolation, SIAM Journal on Numerical Analysis 17(2), 238-246.
#   doi:10.1137/0717021
#   Already implemented as cubic_hermite_spline in src/tools/spline.jl; applied
#   here to log(E), log(S), as srim_sp.jl does.
#
# The stopping-power values themselves have no reference here, by construction:
# they are supplied by the user at runtime. Provenance for a given table belongs
# in that table's own comment lines, not in this file.

"""
    resolve_user_stopping_power_file(file::String)

Resolves the path of a user-supplied stopping-power table. An absolute path is used as
given; anything else is taken relative to `data/` in the package root.

# Input Argument(s)
- `file::String` : file name or path.

# Output Argument(s)
- `path::String` : resolved absolute path.

"""
function resolve_user_stopping_power_file(file::String)
    dir  = joinpath(find_package_root(),"data")
    path = isabspath(file) ? file : joinpath(dir,file)
    if !isfile(path)
        error("The stopping-power file '$file' was not found. Tried '$path'. A relative name " *
              "is taken relative to '$dir', which is where such a table is expected; give an " *
              "absolute path to read it from anywhere else.")
    end
    return abspath(path)
end

"""
    parse_user_stopping_power(path::String)

Parses and validates a two-column stopping-power table.

Column 1 is the projectile's **total kinetic energy in MeV** — not energy per nucleon, and
not the `mₑc²` that `bethe` works in. Column 2 is the **mass** stopping power in
**MeV cm²/g**. The separator is whitespace or a comma, blank lines are skipped, and a line
whose first non-blank character is `#` or `%` is a comment.

Nothing else may be in the file. A third column is an error rather than something to
ignore: a user who pasted in SRIM's output needs to be told, because the extra column there
is *nuclear* stopping and adding it in silently would be wrong by percent at the low end.

Energies must be strictly increasing and are never sorted silently — an out-of-order row is
far more often a corrupt paste than a descending table, and sorting turns it into a smooth
wrong curve.

Only the *file* is validated, never the numbers. A user reaches for this because their
numbers are not Radiant's, so there is no comparison against `bethe`, no warning when the
two disagree, and no check that the curve has a Bragg peak or falls like 1/β² at the top.

# Input Argument(s)
- `path::String` : resolved absolute path.

# Output Argument(s)
- `E_MeV::Vector{Float64}` : energies [in MeV], strictly increasing.
- `S_mass::Vector{Float64}` : mass stopping powers [in MeV cm²/g].

"""
function parse_user_stopping_power(path::String)

    # readlines rather than a streaming eachline: every check below can throw, and a stream
    # left open by an error keeps a lock on the file on Windows.
    E_MeV = Float64[]; S_mass = Float64[]; lines_of = Int64[]
    for (n,line) in enumerate(readlines(path))
        str = strip(line)
        isempty(str) && continue
        (str[1] == '#' || str[1] == '%') && continue
        fields = split(replace(str,","=>" "))
        if length(fields) != 2
            error("Line $n of '$path' has $(length(fields)) fields; a stopping-power table " *
                  "has exactly two, energy [MeV] and mass stopping power [MeV cm²/g]. If this " *
                  "came from SRIM, the extra column is nuclear stopping and must be removed " *
                  "rather than added in.")
        end
        e = tryparse(Float64,fields[1]); s = tryparse(Float64,fields[2])
        if isnothing(e) || isnothing(s)
            error("Line $n of '$path' does not parse as two numbers: '$str'.")
        end
        if !isfinite(e) || !isfinite(s) || e ≤ 0 || s ≤ 0
            error("Line $n of '$path' has a non-positive or non-finite value ($e, $s). The " *
                  "interpolation is in log-log and neither logarithm exists there.")
        end
        if !isempty(E_MeV) && e ≤ E_MeV[end]
            error("Line $n of '$path' has energy $e MeV, which does not exceed the $(E_MeV[end]) " *
                  "MeV on line $(lines_of[end]). Energies must be strictly increasing; they are " *
                  "not sorted for you, because an out-of-order row is more often a corrupt paste " *
                  "than a descending table.")
        end
        push!(E_MeV,e); push!(S_mass,s); push!(lines_of,n)
    end
    if length(E_MeV) < 3
        error("'$path' holds $(length(E_MeV)) data rows; at least 3 are needed, because a " *
              "monotone cubic needs a slope on both sides of a node.")
    end
    return E_MeV, S_mass
end

"""
    user_stopping_power_table(file::String)

Gives the parsed columns of a user-supplied stopping-power table, reading and validating the
file on first use and caching it by resolved path afterwards, so a file registered on ten
materials is read once.

# Input Argument(s)
- `file::String` : file name or path.

# Output Argument(s)
- `E_MeV::Vector{Float64}` : energies [in MeV].
- `S_mass::Vector{Float64}` : mass stopping powers [in MeV cm²/g].

"""
function user_stopping_power_table(file::String)
    path = resolve_user_stopping_power_file(file)
    id = string("user_sp_",path)
    if ~haskey(cache_radiant[],id)
        cache_radiant[][id] = parse_user_stopping_power(path)
    end
    return cache_radiant[][id]
end

"""
    user_stopping_power(file::String,Ei::Float64,ρ::Float64)

Gives the electronic stopping power read from a user-supplied table. Returns a *linear*
stopping power [in mₑc² per cm], as `bethe` does.

**This is the one path with no validation behind it.** Everything else in Radiant's
stopping-power chain was measured against NIST or SRIM; a user table is measured against
nothing, because its correctness is the user's to establish. Radiant validates the *file*;
it cannot validate the *numbers*.

The table must be a **total electronic** stopping power. `sp` computes `Sr = Stot - Sc`
where `Sc` is the catastrophic part it derives from its own differential cross-section, so
feeding in a restricted table subtracts the catastrophic part twice, and a nuclear-stopping
column — or the sum of the two — is the most likely thing to be pasted in by mistake, small
enough at first glance to survive a sanity check while being percent-level and rising below
about 100 keV/u. The soft/catastrophic split is therefore approximate on this path in
exactly the way it is for the `"icru49"` low-energy branch: `Stot` comes from the user's
source and `Sc` from Radiant's analytic cross-section, and the two were never reconciled.

**The density that scales the table is the material's, not the one the table came from.**
For everything except the density effect that is correct, because mass stopping power is
density-independent. But the density effect *is* baked into whatever numbers were supplied,
at whatever density their source assumed, and rescaling cannot touch it. A user who takes a
water-vapour table and hands it to a liquid-water material gets vapour physics scaled to
liquid density, silently. Radiant cannot detect this and does not pretend to.

Interpolation is monotone cubic Hermite (Fritsch-Carlson) on `log(E)` against `log(S)`,
built once and cached by file path. It passes through every tabulated point to round-off,
which is the whole reason a user hands in a table; it cannot overshoot, where a natural
cubic spline rings near the Bragg peak and produces a stopping power that rises where the
data falls; it is local, so one mistyped row corrupts two intervals rather than the whole
curve; and it is two passes over the array rather than a dense 4n x 4n solve. Log-log
because `log S` is close to linear in `log E` above the Bragg peak, so the curve the
interpolant has to follow is nearly straight, and because it guarantees `S > 0` between
nodes.

Outside the tabulated range this errors. It does not extrapolate, does not fall back to
`bethe` and does not clamp to the endpoint: the fix is always to extend the table or to
narrow the group structure, and a silent fallback would put a discontinuity in the middle of
a stopping-power curve where the user cannot see it.

# Input Argument(s)
- `file::String` : file name or path.
- `Ei::Float64` : incoming particle energy [in mₑc²].
- `ρ::Float64` : material density.

# Output Argument(s)
- `S::Float64` : linear stopping power [in mₑc² per cm].

# Reference(s)
- Fritsch and Carlson (1980), Monotone piecewise cubic interpolation.

"""
function user_stopping_power(file::String,Ei::Float64,ρ::Float64)

    mₑ = 0.51099895069
    E_MeV, S_mass = user_stopping_power_table(file)
    E = Ei*mₑ

    if E < E_MeV[1] || E > E_MeV[end]
        error("Energy $E MeV is outside the stopping-power table '$file', which covers " *
              "$(E_MeV[1]) to $(E_MeV[end]) MeV. The table registered with set_stopping_power " *
              "must span the whole energy range of the group structure for that particle: " *
              "extend the table, or narrow the group structure.")
    end

    spline = fast_spline(log.(E_MeV),log.(S_mass),1,string("user_sp_",resolve_user_stopping_power_file(file)))
    return exp(spline(log(E))) * ρ/mₑ
end
