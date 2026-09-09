# Reference(s)
# - Ziegler, J. F., Ziegler, M. D. and Biersack, J. P. (2010), SRIM - The
#   Stopping and Range of Ions in Matter, Nuclear Instruments and Methods B
#   268(11-12), 1818-1823. doi:10.1016/j.nimb.2010.02.091
#   Tables here were produced with SRIM-2013 (SR Module reporting
#   SRIM-2012.01), output units 6 = MeV/(g/cm2).
# - Ziegler, J. F. and Manoyan, J. M. (1988), The stopping of ions in compounds,
#   Nuclear Instruments and Methods B 35(3-4), 215-228.
#   doi:10.1016/0168-583X(88)90273-X
#   The "cores and bonds" model behind SRIM's compound correction. SRIM stores
#   only the per-molecule bond inventory (Data/Compound.dat); the bond stopping
#   coefficients are not in any shipped data file.
# - Fritsch, F. N. and Carlson, R. E. (1980), Monotone piecewise cubic
#   interpolation, SIAM Journal on Numerical Analysis 17(2), 238-246.
#   doi:10.1137/0717021
#   Already implemented as cubic_hermite_spline in src/tools/spline.jl; applied
#   here to log(E), log(S).

"""
    srim_phase(state_of_matter::String)

Maps Radiant's state of matter onto SRIM's binary phase flag. `"gaz"` selects `"gas"`;
`"solid"` and `"liquid"` both select `"solid"`, because SRIM has no liquid.

The phase follows the material the user declared and is never decided from the element.
All 92 elements have a different gas table, not only the ones that are gases at room
temperature: 87 of 92 differ by more than 5% for protons, and the largest differences sit
on lead, titanium, gadolinium and barium at about 8.8%. SRIM applies a near-uniform
gas-phase enhancement to whatever is declared a gas, so a solid containing oxygen must not
pick up oxygen's gas table.

# Input Argument(s)
- `state_of_matter::String` : material state.

# Output Argument(s)
- `phase::String` : `"solid"` or `"gas"`.

"""
function srim_phase(state_of_matter::String)
    if state_of_matter == "gaz"
        return "gas"
    elseif state_of_matter ∈ ["solid","liquid"]
        return "solid"
    else
        error("The state of matter is either solid, liquid or gaz.")
    end
end

"""
    srim_stopping_power(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,
    particle::Particle,state_of_matter::String)

Gives the electronic stopping power from SRIM-2013's own tabulated values, mixed from the
elements the way SRIM does it internally. Returns a *linear* stopping power [in mₑc² per
cm], as `bethe` does.

**The point is fidelity to SRIM, not accuracy.** A user turns this on to reproduce SRIM
numbers — to compare against a SRIM-based workflow, or to feed a code that was tuned on
them. Everything in the corrected Bethe path is about getting the physics right; this is
about getting *SRIM* right, and the two are not the same thing.

SRIM's compound calculation is, in full,

    S_elec = ( Σⱼ ωⱼ S_elecⱼ(E,phase) ) * ( 1 + (CC - 1)*g(E) )

and this implements it with `CC` pinned to 1, so the second factor is identically 1 and
never appears. Every part of the first factor was measured rather than assumed: Bragg
mixing is exact at CC = 1 (H2O from separate elemental runs against SRIM's own two-element
compound run, 114 points, max 0.058%), and the phase flag is nothing but a choice of
elemental table (C2H2 and CH4 in both phases, max 0.076%). `ωⱼ` are mass fractions, the
same rule Bragg additivity uses and for the same reason; this is not the Eq. (53) electron
weighting that applies inside the Bethe bracket.

**The compound correction is deliberately not implemented, and there is no way to set it.**
CC is a property of the molecule, not of the interaction, so one setting on
`Inelastic_Collision` would apply one molecule's bond correction to every material in the
geometry. The numbers are not available either: SRIM computes CC inside SR.exe from the
bond inventory in Data/Compound.dat by the Ziegler-Manoyan "cores and bonds" model, and the
bond stopping coefficients are in no shipped SRIM data file — the only way to obtain one is
to read `Bragg Correction = X%` off SRIM's Compound Dictionary by hand. And what SRIM would
supply is questionable anyway: its own help states the correction is computed for He ions
and that *"It is not known what is the compound correction for heavy ions (Z > 3)"*, yet it
applies that He-derived number to protons unchanged.

**What that costs, plainly.** For a molecular target the correction SRIM's dictionary would
have applied reaches about 6% below 1 MeV — liquid water is -6.00% flat from 10 keV to
400 keV — fading to zero above about 2.5 MeV/u. So below a few MeV this reproduces *SRIM
run on a composition*, which is not the same thing as *SRIM's compound library entry*. A
user comparing against a number read off the SRIM GUI for a compound in the dictionary will
see that difference.

**SRIM has no Fermi density effect.** Not a crude one — none. Above roughly 500 MeV/u these
tables are therefore wrong in a way that grows with energy: for liquid water the density
effect removes about 1.8% of the stopping number at 2 GeV and about 8.5% at 10 GeV, and
none of it is here. It is also why SRIM cannot tell liquid water from water vapour except
through the gas flag. The tables still run to 1 GeV/u because a user comparing against SRIM
wants SRIM's numbers over SRIM's whole range, including where they are wrong — but this is
the wrong switch for high-energy work, and the corrected Bethe path is the right one.

Interpolation is monotone cubic Hermite (Fritsch-Carlson) on `log(E)` against `log(S)`,
built once per (element, particle, phase) and cached. It is chosen over a natural cubic
spline for robustness rather than accuracy: it is local, so one bad point cannot ring
through the table, it is non-oscillatory by construction, and `natural_cubic_spline` would
assemble and solve a dense 4n x 4n system — 3120 x 3120 for 780 points — per element per
phase for no measurable gain. The shipped grid is never subsampled or resampled: it carries
exact nodes at SRIM's internal handovers, 1 MeV for protons and 8 MeV for alphas, with
helium having a third of its own at 2 keV.

Outside the tabulated range this errors rather than extrapolating or falling back to
`bethe`, because a quiet fallback puts a discontinuity in the middle of a stopping-power
curve where the user cannot see it.

# Input Argument(s)
- `Z::Vector{Int64}` : atomic numbers of the elements in the material.
- `ωz::Vector{Float64}` : weight fraction of the elements composing the material.
- `ρ::Float64` : material density.
- `Ei::Float64` : incoming particle energy [in mₑc²].
- `particle::Particle` : incoming particle.
- `state_of_matter::String` : material state.

# Output Argument(s)
- `S::Float64` : linear stopping power [in mₑc² per cm].

# Reference(s)
- Ziegler, Ziegler and Biersack (2010), SRIM - The Stopping and Range of Ions in Matter.
- Ziegler and Manoyan (1988), The stopping of ions in compounds.
- Fritsch and Carlson (1980), Monotone piecewise cubic interpolation.

"""
function srim_stopping_power(Z::Vector{Int64},ωz::Vector{Float64},ρ::Float64,Ei::Float64,particle::Particle,state_of_matter::String)

    mₑ = 0.51099895069

    #----
    # Resolve the ion and the phase
    #----
    if is_proton(particle)
        ion = "proton"
    elseif is_alpha(particle)
        ion = "alpha"
    else
        error("There is no SRIM stopping-power table for the '$(get_type(particle))' particle. " *
              "set_use_srim_sp is only supported for Proton and Alpha; turn it off, or restrict " *
              "the calculation to those particles.")
    end
    phase = srim_phase(state_of_matter)

    data = fast_load("srim_stopping_power.jld2")
    E_grid = data[ion]["E"]
    S_elec = data[ion][phase]["S_elec"]

    #----
    # Range check, on the energy and on the elements
    #----
    E_MeV = Ei*mₑ
    if E_MeV < E_grid[1] || E_MeV > E_grid[end]
        error("Energy $E_MeV MeV is outside the SRIM table for $ion, which covers " *
              "$(E_grid[1]) to $(E_grid[end]) MeV. Narrow the group structure, or turn " *
              "set_use_srim_sp off.")
    end
    for z in Z
        if z < 1 || z > 92
            error("SRIM stopping-power tables cover Z = 1 to 92 and the material contains " *
                  "Z = $z. Turn set_use_srim_sp off for this material.")
        end
    end

    #----
    # Bragg mix over the elements, in mass fractions
    #----
    # The interpolant is built once per (element, particle, phase) and cached: building it
    # inside sp would rebuild a 780-point derivative array at every quadrature point of
    # every energy group.
    logE = log(E_MeV)
    S = 0.0
    for i in eachindex(Z)
        spline = fast_spline(log.(E_grid),log.(S_elec[Z[i],:]),Z[i],string("srim_sp_",ion,"_",phase))
        S += ωz[i] * exp(spline(logE))
    end

    # The tables hold a mass stopping power and sp returns a linear one. The density here is
    # the material's, which is the correct one; the 1.0 g/cm³ of the generating runs never
    # enters, because SRIM's mass stopping power is exactly density-independent.
    return S * ρ/mₑ
end
