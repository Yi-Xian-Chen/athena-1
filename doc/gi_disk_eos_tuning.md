# GI Disk EOS Atmosphere Tuning Notes

This tracks initialization hyperparameters for the 3D Helmholtz GI disk problem,
especially the low-density atmosphere/floor treatment.

## Test Objective

The immediate goal is not to produce a final GI disk model. The current test is
to check the cooling and electron-capture behavior of a slowly cooling disk with
the modified Helmholtz alpha EOS.

The intended development sequence is:

1. Start with a quiet rotating disk and verify that cooling, heating guards,
   `Ye` evolution, and alpha recombination/dissociation behave sensibly.
2. If the disk becomes dynamically violent before the source behavior can be
   interpreted, improve the hydrodynamic initialization first. Candidate fixes
   include pressure-corrected rotation and smoother/geometrically consistent
   disk tapering.
3. After the non-self-gravitating cooling test is under control, add
   self-gravity and let gravito-turbulence develop.

## Current Baseline

- Problem generator: `src/pgen/gi_disk_3d_eos_v2.cpp`
- Smoke input: `inputs/hydro/athinput.gi_disk_3d_eos_v2_smoke`
- Semi-production input: `inputs/hydro/athinput.gi_disk_3d_eos_v2_semiprod`
- Units: cgs
- Central mass: `M = 10 Msun`
- Disk region: `r_in = 50 rg`, `r_out = 151 rg`
- Domain outer radius in smoke input: `x1max = 180 rg`
- Initial disk electron fraction: `Ye0 = 0.1`
- Initial pressure ratio: `(P/rho)(100 rg) = 1.0e17 erg/g`
- Pressure-ratio slope: `(P/rho) proportional to r^-1`
- Density profile: fixed `rho_mid proportional to r^-3`, Gaussian in `z/H`
- Disk truncation: spherical shell `r_in < r < r_out`; Gaussian vertical
  structure is used only inside this shell
- Self-gravity: off
- Mesh generation: off
- Source terms: QW-style energy and Ye source terms available

## Atmosphere Hyperparameters

| Parameter | Meaning | Current/default | Notes |
| --- | --- | --- | --- |
| `rho_floor` | Minimum density used by initialization and floors | `1.0e6 g/cm^3` | Main knob for EOS conditioning and atmosphere inertia. |
| `T_vacuum` | Optional prescribed atmosphere temperature | `-1.0` disabled | If positive, atmosphere pressure/energy are initialized from `rho_floor, T_vacuum, Ye_vacuum`. |
| `Ye_vacuum` | Optional atmosphere electron fraction | `Ye0` | Used only in cells classified as atmosphere; disk cells keep `Ye0`. |
| `qw_rho_min` | Density below which net cooling is clipped off | `1.0e7 g/cm^3` | Allows heating but suppresses cooling in very low-density cells. Does not fix EOS inversion problems. |

## Composition Target

For the alpha Saha model in `helmholtz.cpp`, high temperature dissociates alpha
particles. Thus a hot collapsar-like atmosphere/funnel should be alpha-poor, not
alpha-rich.

At `rho = 1.0e6 g/cm^3` using the current Saha implementation:

| `T` | `Ye` | `Xalpha` behavior |
| --- | --- | --- |
| `6.0e9 K` | `0.1` | `Xalpha ~ 0.2`, alpha-saturated |
| `6.0e9 K` | `0.4` | `Xalpha ~ 0.8`, alpha-saturated |
| `8.0e9 K` | `0.1` | `Xalpha ~ 0.124`, transition regime |
| `8.0e9 K` | `0.4` | `Xalpha ~ 0.463`, transition regime |
| `1.0e10 K` | `0.1` | `Xalpha ~ 1.0e-4`, alpha-poor |
| `1.0e10 K` | `0.4` | `Xalpha ~ 7.0e-4`, alpha-poor |

## Known Run Behavior

With `rho_floor = 1.0e6 g/cm^3` and the `10 Msun`, `Ye0 = 0.1`,
`(P/rho)(100 rg) = 1.0e17` disk:

| Atmosphere setup | Result |
| --- | --- |
| `T_vacuum = -1.0` | Runs for at least 5 tiny cycles with mild bulk cooling. |
| `T_vacuum = 1.0e9 K` | Runs for a few tiny cycles, but this is a cool alpha-rich atmosphere. |
| `T_vacuum = 3.0e9 K` | Fails dynamically after cycle 1 in pressure/EOS recovery. |
| `T_vacuum = 5.0e9 K` | Fails immediately after cycle 0 with near-floor density pathology. |
| `T_vacuum = 6.0e9 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |
| `T_vacuum = 8.0e9 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |
| `T_vacuum = 1.0e10 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |

Changing the atmosphere electron fraction at `rho_floor = 1.0e6 g/cm^3`
does not fix the hot-atmosphere failures:

| `Ye_vacuum` | `T_vacuum` | Result |
| --- | --- | --- |
| `0.4` | `8.0e9 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |
| `0.4` | `1.0e10 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |
| `0.4` | `1.2e10 K` | Fails immediately in pressure inversion at `rho = 1.0e6`. |
| `0.5` | `8.0e9 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |
| `0.5` | `1.0e10 K` | Fails immediately in Helmholtz energy inversion at `rho = 1.0e6`. |
| `0.5` | `1.2e10 K` | Fails immediately in pressure inversion at `rho = 1.0e6`. |

At `rho_floor = 1.0e8 g/cm^3`, the hot-atmosphere branch becomes much more
EOS-invertible, but the first hydro update can still fail if the atmosphere
pressure jump is too large.

With `Ye_vacuum = Ye0 = 0.1`:

| `T_vacuum` | Initial atmosphere `Xalpha` | Result |
| --- | --- | --- |
| disabled | `0.2` | Runs for 5 tiny cycles. |
| `8.0e9 K` | `0.1999` | Runs for 5 tiny cycles. |
| `1.0e10 K` | `0.1901` | Runs for 5 tiny cycles. |
| `1.05e10 K` | not recorded in output | Fails after cycle 0 in pressure/sound-speed inversion. |
| `1.1e10 K` | not recorded in output | Fails after cycle 0 after producing negative density. |
| `1.2e10 K` | `0.0700` | Fails after cycle 0 after producing negative density. |

Thus for `Ye_vacuum = 0.1`, the small-cycle stability edge is between
`T_vacuum = 1.0e10 K` and `1.05e10 K`. However, even at `1.0e10 K`, the
atmosphere is still close to alpha-saturated because the higher density favors
recombination.

With `Ye_vacuum = 0.5`:

| `T_vacuum` | Initial atmosphere `Xalpha` | Result |
| --- | --- | --- |
| `8.0e9 K` | `0.9810` | Fails after cycle 0 in pressure/sound-speed inversion. |
| `1.0e10 K` | `0.8182` | Fails after cycle 0 after producing negative density. |
| `1.05e10 K` | `0.7246` | Fails after cycle 0 after producing negative density. |
| `1.1e10 K` | `0.6043` | Fails after cycle 0 after producing negative density. |

At this density, raising `Ye_vacuum` to `0.5` makes the atmosphere more
alpha-rich, not alpha-poor, and is dynamically less stable in the current setup.

## Pressure Scale Warning

At `T_vacuum = 1.0e10 K`, radiation pressure alone is

```text
P_rad = a T^4 / 3 ~ 2.5e25 erg/cm^3.
```

This is only about 10-100 times below the midplane pressure in the present thin
disk setup, so a hot alpha-poor atmosphere is not dynamically vacuum-like.

Approximate cold electron degeneracy pressure becomes comparable to this
radiation pressure around:

| `Ye` | density where `P_e,deg ~ P_rad(T=1e10 K)` |
| --- | --- |
| `0.5` | `rho ~ 1.0e8 g/cm^3` |
| `0.1` | `rho ~ 1.0e9 g/cm^3` |

Raising `rho_floor` may help EOS conditioning, but it also makes the atmosphere
more massive and dynamically important.

## Next Scan

The next useful tuning scan is:

```text
rho_floor = 1e6, 1e7, 1e8, 1e9
T_vacuum  = disabled, 8e9, 1e10, 1.2e10
Ye_vacuum = Ye0, 0.4, 0.5
```

Diagnostics to record:

- cycle-0 inversion success
- first 5-cycle stability
- atmosphere `Xalpha`
- atmosphere `P / P_mid`
- minimum/maximum `rho`, `P`, `T`, `Ye`, `Xalpha`
- whether low-density cells generate large velocities or energy changes

## Semi-Production Baseline

For now, the preferred semi-production atmosphere is deliberately simple:

```text
rho_floor = 1.0e7 g/cm^3
T_vacuum  = disabled
Ye_vacuum = Ye0
```

This avoids imposing a hot pressure-supported funnel by hand. The low-density
material starts as the continuation/floor of the disk initialization, and later
turbulent or irradiation heating can move it toward the alpha-poor branch if the
run supplies enough energy.

The matching runtime density floor should also be set:

```text
hydro/dfloor = 1.0e7
```
