# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Interactive development
julia --project          # then: using Radiant

# Build documentation
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

To activate the environment in the REPL: `] activate .` then `instantiate`.

## Source layout

`src/Radiant.jl` loads four directories in this fixed order (types must exist before they are used):

| Directory | Role |
|---|---|
| `structures/` | Mutable structs with `set_*` / `get_*` methods and `build()` |
| `tools/` | Math utilities: quadratures, polynomials, interpolation, integrals, caching |
| `cross_sections/` | Physics models that compute multigroup cross-section data |
| `particle_transport/` | Transport sweeps, source construction, flux extraction |

Tabulated nuclear/atomic data lives in `data/*.jld2` (Seltzer-Berger, EPDL97, Mott, etc.) and is loaded lazily via `JLD2` when cross-sections are built.

## User-facing workflow

Every struct follows the same builder pattern: configure via `set_*()` methods, then call `build()`. `Computation_Unit.run()` auto-builds any component that hasn't been built yet.

```
Material
  └─ add_element / set_density
Cross_Sections
  └─ set_source / set_materials / set_particles / set_interactions / set_group_structure / set_legendre_order
Geometry
  └─ set_type / set_dimension / set_number_of_regions / set_voxels_per_region / set_region_boundaries / set_material_per_region / set_boundary_conditions
Discrete_Ordinates (one per particle)
  └─ set_particle / set_solver_type / set_quadrature / set_scheme / set_legendre_order
Solvers
  └─ add_solver (one per particle) / set_maximum_number_of_generations
Fixed_Sources
  └─ add_source (Volume_Source or Surface_Source)
Computation_Unit
  └─ set_cross_sections / set_geometry / set_solvers / set_sources / run()
```

## Python-style method notation

`src/tools/python_method_notation.jl` overloads `Base.getproperty` for all exported structs via the `RadiantObject` union type, making `obj.set_density(1.0)` and `set_density(obj, 1.0)` equivalent. Adding a new exported struct to `RadiantObject` in that file gives it this notation automatically.

## Interaction system

`src/structures/Interaction.jl` defines the abstract `Interaction` type. Each concrete subtype (`Compton`, `Elastic_Collision`, `Inelastic_Collision`, `Bremsstrahlung`, `Photoelectric`, `Pair_Production`, `Rayleigh`, `Annihilation`, `Relaxation`) holds model-selection flags and an `interaction_types` dict mapping `(incoming_particle_type, outgoing_particle_type) → [process_strings]`.

Adding a new interaction requires:
1. New struct in `src/structures/`
2. Physics implementation in `src/cross_sections/`
3. Entry in `src/cross_sections/interaction_interdependances.jl`

## Solver types

`Discrete_Ordinates` and `Spherical_Harmonics` (in `structures/`) configure per-particle numerical methods. Solver type strings:

- `"BTE"` – full Boltzmann scattering integral
- `"BFP"` – Boltzmann-Fokker-Planck (BTE + FP approximation for forward-peaked scattering)
- `"FP"` – pure Fokker-Planck
- `"CSD"` – continuous slowing-down only

Transport sweep kernels live in `particle_transport/sn_*D_bte.jl`, `sn_*D_bfp.jl`, `pn_*D_*.jl`. The `scheme_weights.jl` file computes closure-relation weights for the spatial and energy polynomial bases (Diamond Difference `"DD"`, Discontinuous Galerkin `"DG"`, Adaptive Weighted `"AWD"`).

## Coupled transport

`src/particle_transport/transport.jl` iterates over *generations*: each pass computes flux for every particle and feeds cross-particle scattering into the source for the next particle. Convergence is checked on flux, energy-deposition, or charge-deposition norms. Single-particle problems skip the generation loop entirely.
