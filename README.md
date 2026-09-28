athena
======
<!-- Jenkins Status Badge in Markdown (with view), unprotected, flat style -->
<!-- In general, need to be on Princeton VPN, logged into Princeton CAS, with ViewStatus access to Jenkins instance to click on unprotected Build Status Badge, but server is configured to whitelist GitHub -->
<!-- [![Jenkins Build Status](https://jenkins.princeton.edu/buildStatus/icon?job=athena/PrincetonUniversity_athena_jenkins_master)](https://jenkins.princeton.edu/job/athena/job/PrincetonUniversity_athena_jenkins_master/) -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4455880.svg)](https://doi.org/10.5281/zenodo.4455880) <!-- v21.0, not Concept DOI that tracks the "latest" version (erroneously sorted by DOI creation date on Zenodo). 10.5281/zenodo.4455879 -->
[![Travis CI Build Status](https://travis-ci.com/PrincetonUniversity/athena.svg?token=Ejzw3yndG1Fqub679gCB&branch=master)](https://travis-ci.com/PrincetonUniversity/athena)
[![codecov](https://codecov.io/gh/PrincetonUniversity/athena/branch/master/graph/badge.svg?token=ZzniY084kP)](https://codecov.io/gh/PrincetonUniversity/athena)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.0-4baaaa.svg)](code_of_conduct.md)

<!--[![Public GitHub  issues](https://img.shields.io/github/issues/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/issues)
[![Public GitHub pull requests](https://img.shields.io/github/issues-pr/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/pulls) -->

<p align="center">
	  <img width="345" height="345" src="https://user-images.githubusercontent.com/1410981/115276281-759d8580-a108-11eb-9fc9-833480b97f95.png">
</p>

Athena++ GRMHD code and adaptive mesh refinement (AMR) framework

Please read [our contributing guidelines](./CONTRIBUTING.md) for details on how to participate.

## Helmholtz-alpha GI disk v4 working notes

- The current problem generator is `gi_disk_3d_eos_v4`. It combines the v3 disk
  initialization and Helmholtz-alpha EOS, the v3.5 nonuniform polar mesh generator,
  the v3.6 masked polar cones, and the spherical-harmonic self-gravity port.
- The v4.1 variant `gi_disk_3d_eos_v41` adds
  `hydro/helm_Tfloor_temperature`, which defaults to the Helmholtz table minimum
  ($10^3\,\mathrm{K}$) and may be raised without changing the table. Its atmosphere
  and masked cones inherit this temperature unless `problem/T_vacuum` is set
  explicitly. Keep `efloor\le e(\rho_{\rm floor},T_{\rm floor},Y_{e,\rm vacuum})`
  and `pfloor\le P(\rho_{\rm floor},T_{\rm floor},Y_{e,\rm vacuum})`; otherwise
  the hydro floors produce a hotter atmosphere. The tested setup used
  $T_{\rm floor}=10^6\,\mathrm{K}$ and `efloor=pfloor=10^23` in cgs units.
  Its four-rank smoke run completed $1\,\Omega_{100}^{-1}$ at cycle 125; the
  atmosphere retained a $10^6\,\mathrm{K}$ minimum while dynamical atmosphere
  cells heated as high as $1.07\times10^{10}\,\mathrm{K}$. The high atmosphere
  still loaded the domain strongly, increasing $M/M_0$ to $1.18479$ in one orbit.
- The v4.2 variant `gi_disk_3d_eos_v42` optionally adds the disk's enclosed
  monopole mass to the initial centrifugal balance with
  `problem/rotation_include_disk_monopole=true`. The correction is evaluated as
  a function of cylindrical radius from the configured surface-density profile,
  rather than applying a constant velocity multiplier. In the coarse Model 5
  one-orbit test ($Q=0.5$, $Y_e=0.5$, $\hat c_s=0.05$), the retained mass improved
  from $M/M_0=0.53929$ without the correction to $0.78864$ with it. The supplied
  v4.2 input is the corrected, radial-diode case with $r_{\rm out}=200\,r_g$.
- Configure with spherical-polar coordinates, the general Helmholtz EOS, one passive
  scalar for $Y_e$, spherical-harmonic gravity, MPI, and FFT support. The Helmholtz
  table is read from `data/eos/helm_table.dat` when compiling and running from the
  repository root.
- The reproducible low-resolution diagnostic input is
  `inputs/hydro/athinput.gi_disk_v4_model1_r240_lowrho_inneroutflow_outerdiode_4omega_4rank`.
  It covers $40$--$240\,r_g$ with a constant $\Delta r/r$ grid of
  $20\times48\times16$ cells and four $10\times48\times8$ meshblocks/ranks. Each
  block spans the complete nonuniform theta domain, including both polar boundaries.
- Model 1 uses $Q_{\rm ref}=2$, $Y_e=0.1$, and
  $(P/\rho)_{100r_g}=8.987551787368176\times10^{16}\,\mathrm{erg\,g^{-1}}$.
  The active disk is initially truncated at $50$ and $151\,r_g$; the larger numerical
  radial domain supplies an outer buffer.
- The tested dilute-atmosphere settings are `rho_floor=1e6`, `dfloor=1e6`,
  `pfloor=1e23`, `efloor=1e23`, and `helm_Tfloor=true` in cgs units. Keep
  `qw_rho_min=2e7` so the dilute numerical atmosphere cannot cool or deleptonize
  through the QW source.
- The atmosphere `rho_floor` can also be treated as a tunable numerical mass-loading
  parameter. A moderately higher floor supplies gradual infall that may drive the disk
  toward $Q\sim1$ when physical cooling is too slow on the simulated timescale. This
  should be calibrated with $M(t)/M(0)$ and resolution/floor comparisons: it represents
  atmosphere- or boundary-fed loading, not a physical accretion prescription.
- The current radial-boundary experiment uses Athena++ built-in outflow at the inner
  radial face and the pgen user diode at the outer face. The theta mask uses outflow
  copying at `theta_cut=pi/2-1` while the physical theta boundaries remain polar.
- Commit `02cc33cc` makes `helm_Tfloor` tolerant of roundoff at the Helmholtz table
  minimum and prevents HLLC's temporary PVRS wave-speed estimates from querying the
  EOS with nonpositive density or pressure. The HLLC guard affects intermediate
  wave-speed estimates, not the evolved cell state.
- The four-rank diagnostic completed $4\,\Omega_{100}^{-1}$ with smooth timesteps.
  Its mass changed from $2.47382\times10^{33}$ to $2.57926\times10^{33}\,\mathrm{g}$,
  or $M/M_0=1.04262$. For comparison, the older $2\times10^7\,\mathrm{g\,cm^{-3}}$
  atmosphere reached $M/M_0=1.36530$ over the same interval, showing that atmosphere
  density dominated the artificial mass growth.
- Restarting the dilute-atmosphere diagnostic from $4$ to
  $10\,\Omega_{100}^{-1}$ also completed cleanly (cycle 1919). The mass peaked at
  $M/M_0=1.04620$ near $1.34\,\Omega_{100}^{-1}$ and declined to
  $2.55497\times10^{33}\,\mathrm{g}$, or $M/M_0=1.03280$, at the final time; there
  is no continuing secular mass growth in this low-resolution test.
- At $4\,\Omega_{100}^{-1}$ and $r\simeq103\,r_g$ on the midplane, the dilute run has
  $\rho=5.36\times10^8\,\mathrm{g\,cm^{-3}}$, $P=7.19\times10^{25}\,\mathrm{erg\,cm^{-3}}$,
  $T=1.26\times10^9\,\mathrm{K}$, $Y_e=0.1$, and $X_\alpha=0.2$.

## Citation
To cite Athena++ in your publication, please use the following BibTeX to refer to the code's [method paper](https://ui.adsabs.harvard.edu/abs/2020ApJS..249....4S/abstract):
```
@article{Stone2020,
	doi = {10.3847/1538-4365/ab929b},
	url = {https://doi.org/10.3847%2F1538-4365%2Fab929b},
	year = 2020,
	month = jun,
	publisher = {American Astronomical Society},
	volume = {249},
	number = {1},
	pages = {4},
	author = {James M. Stone and Kengo Tomida and Christopher J. White and Kyle G. Felker},
	title = {The Athena$\mathplus$$\mathplus$ Adaptive Mesh Refinement Framework: Design and Magnetohydrodynamic Solvers},
	journal = {The Astrophysical Journal Supplement Series},
}
```
Additionally, you can add a reference to `https://github.com/PrincetonUniversity/athena` in a footnote.

Finally, we have minted DOIs for each released version of Athena++ on Zenodo. This practice encourages computational reproducibility, since you can specify exactly which version of the code was used to produce the results in your publication. `10.5281/zenodo.4455879` is the DOI which cites _all_ versions of the code; it will always resolve to the latest release. Click on the Zenodo badge above to get access to BibTeX, etc. info related to these DOIs, e.g.:

```
@software{athena,
  author       = {Athena++ development team},
  title        = {{PrincetonUniversity/athena-public-version: 
                   Athena++ v21.0}},
  month        = jan,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {21.0},
  doi          = {10.5281/zenodo.4455880},
  url          = {https://doi.org/10.5281/zenodo.4455880}
}
```
