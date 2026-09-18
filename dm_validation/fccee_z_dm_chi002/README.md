# DM validation run: FCC-ee Z chi_e=0.002, replacing LCFA mean-kick with DM stochastic recoil

## Objective

This reruns the FCC-ee Z-pole long-term-tracking case published in
`LaserPol_v7.pdf` (Fig. 6/7, chi_e=0.002 operating point) -- but using the
density-matrix (DM) stochastic-recoil laser module (branch `density-matrix`
of this repo) instead of the original LCFA mean-kick module, to test the
paper's own stated open question:

> "The present BMAD implementation applies the mean LCFA recoil from the
> precomputed lookup table rather than sampling individual photon emissions.
> Therefore the results should be interpreted as a mean-kick beam-survival
> estimate; stochastic photon-number and photon-energy fluctuations, which
> may enhance energy-spread growth and rare-particle loss, will be treated
> in future work using the full density-matrix framework."

i.e.: does genuine photon-number/photon-energy stochasticity change the
particle-loss threshold or polarization buildup rate found with the
mean-kick LCFA table at this operating point?

**Do not silently substitute a stock BMAD build or the LCFA module.** This
run requires `bmad/laser/dm_interp_mod.f90` and the DM-only
`bmad/laser/laser_tracking_mod.f90` from THIS branch (`density-matrix`,
commit `a996febfc` at time of writing -- check `git log -1` after cloning),
plus the three DM lookup tables in `bmad/laser/` (`dm_energy_table.dat`,
`dm_spin_table.dat`, `dm_cdf_table.dat`).

## Provenance: what changed vs. the original LCFA reference, and why

The original validated LCFA run (`reference/xip05_2_out_cons_LCFA_reference.ave`,
handed off previously as `SDCC_xip05_2_1700` in the paper archive) used:
- `fccee_z.bmad`'s `laser_ip` element with `LASER_XI` ramped to **0.1**,
  `LASER_WL_UM := 10.6` (a 10.6 micron CO2 laser), theta implicit in the
  LCFA module's head-on-collision formula (no explicit angle input).
- Final result: `<Sy> = 0.006893103` (0.6893%) at turn 100000, with all
  particles retained (920 live at turn 0 in that specific archived run;
  the published paper text reports 2214 macro-particles retained with no
  loss for its chi_e=0.002 case).

**The LCFA module only ever consumed the resulting chi_e** -- computed
internally from `LASER_XI`, `LASER_WL_UM`, and the electron's gamma via a
head-on-collision formula -- not the wavelength or geometry themselves.
Per that formula, `LASER_XI=0.1` at 10.6 micron and E=45.6 GeV gives
chi_e ~ 0.0041 (a factor ~2 above the paper's stated 0.002 benchmark; the
paper's Table III geometry-aware chi_e differs from the module's simplified
head-on formula, and the exact 0.6893% figure in the archived reference is
close enough to the paper's stated "0.69%" that this is evidently the
intended benchmark case regardless of the small formula-convention gap).

**The DM tables are built at a fixed 1064 nm** (Nd:YAG) -- wavelength is not
a table axis, so this run cannot reuse the CO2 configuration directly.
Since only chi_e mattered physically for the original LCFA result, this run
targets the **same chi_e=0.002** at 1064 nm instead, using a round,
already-well-characterized geometry (theta=10 deg, matching the "FCC-ee
nominal" point used throughout the DensityMatrix/ Python validation work)
and solving for the LASER_XI (=a0) that gives chi_e=0.002 there:

```
$ python3 -c "
import sys; sys.path.insert(0,'DensityMatrix/src')
from constants_natural import gamma_from_energy, compute_chi_e
from scipy.optimize import brentq
gamma = gamma_from_energy(45.6)   # GeV in, dimensionless gamma out
f = lambda a0: compute_chi_e(gamma, a0, 1064.0, 10.0) - 0.002
print(brentq(f, 1e-4, 0.85))
"
0.3556232088098463
```

This value (`0.3556232088098463`) is well inside the DM tables' validated
range (a0 <= 0.854 usable, theta >= 6 deg usable -- see
`bmad/laser/dm_interp_mod.f90`'s module header), so no new table build was
needed for this comparison.

**Everything else is unchanged** from the validated LCFA reference: same
lattice (`fccee_z.bmad`, positron beam, E_tot=45.6 GeV), same RF (90 MV),
same radiation damping/fluctuations, spin tracking, and aperture settings,
same emittances/bunch length/energy spread, same 100-turns-on/1-turn-off
laser gating pattern and z-sweep, same random seed policy. Only `n_particle`
was set to 2214 (matching the published paper's text exactly, rather than
the 920/1700 used in some earlier exploratory archives of this same case),
and the laser element gained a `LASER_THETA_DEG := 10.0` attribute (new,
required by the DM module; see `fccee_z.bmad`'s comments for the full diff).

## Files

- `fccee_z.bmad`: patched FCC-ee Z lattice (custom_attribute11=LASER_THETA_DEG
  added; laser_ip element gets `LASER_THETA_DEG := 10.0`). Diff against the
  original is confined to the header comment block and the `laser_ip` line.
- `lat.bmad`: unchanged top-level wrapper.
- `laser_ramper_expression.bmad`: patched only in the `LASER_XI` peak value
  (0.1 -> 0.3556232088098463); timing/z-sweep logic unchanged.
- `long_term_tracking.init`: n_particle=2214, DM-run output filenames,
  otherwise identical to the validated LCFA reference's tracking init.
- `reference/xip05_2_out_cons_LCFA_reference.ave`: the original LCFA
  reference output (920-particle archived run) for side-by-side comparison.
- `run.sh`, `submit.slurm`, `validate_output.py`, `SHA256SUMS`: as before,
  adapted for this run's file/output names.

## Build (SDCC)

1. `git clone git@github.com:vranjbarbnl/bmad-ecosystem.git && cd bmad-ecosystem && git checkout density-matrix`
2. Build following the normal Bmad-ecosystem production build
   (`util/dist_build_production -r` or equivalent, per SDCC's Bmad module
   environment) -- confirm the resulting `long_term_tracking` links against
   a `bmad` library containing `laser_tracking_mod` (DM version, calls
   `dm_sample_energy_kick` / `dm_spin_transfer_matrix`, NOT `lcfa_map`) and
   `dm_interp_mod`.
3. Sanity check before the long run: `bmad/laser/test_dm_interp` (build with
   `gfortran -O2 -std=f2008 dm_interp_mod.f90 test_dm_interp.f90 -o
   test_dm_interp` from `bmad/laser/`) should print `16 passed, 0 failed`.

## Run

From this directory:
```
./run.sh
```
`run.sh` checks for the DM tables in `../../bmad/laser/` and symlinks them
into the run directory if not already present (the module looks for
`dm_energy_table.dat`, `dm_spin_table.dat`, `dm_cdf_table.dat` in the
working directory by default), then invokes `long_term_tracking`.

For batch: `sbatch submit.slurm` (fill in SDCC account/partition/time first).

## Acceptance criteria before treating results as meaningful

- `beam_init%n_particle = 2214`, `ltt%n_turns = 100000`, `ltt%rfcavity_on = T`.
- Radiation damping, radiation fluctuations, spin tracking, and aperture
  limiting are all on.
- Reference energy is 45.6 GeV; RF remains 90 MV (unchanged in `fccee_z.bmad`).
- `LASER_THETA_DEG = 10.0` (fixed) and `LASER_XI` peaks at
  `0.3556232088098463` during its ON window (100 turns on, 1 turn off,
  z-sweep unchanged).
- The three DM table files are present and match the checksums recorded
  for them in `bmad/laser/` on this branch (not re-verified here --
  `SHA256SUMS` in this directory covers only the files unique to this run).
- Fixed random seed 3972547 unless a deliberate statistical-replica run is
  requested.

## What to compare against the LCFA reference

From `reference/xip05_2_out_cons_LCFA_reference.ave`: final `<Sy> =
0.006893103` (0.6893%), 920/920 particles retained (turn 0 -> turn 100000,
no loss). The DM run is NOT expected to reproduce this number exactly --
different wavelength-equivalent operating point, different (stochastic vs.
mean) recoil model, different n_particle -- but should be examined for:

1. **Particle loss**: does the DM run's live-particle count drop below its
   own turn-0 value anywhere before turn 100000, where the LCFA run (at a
   comparable chi_e) showed none? This is the central question this
   validation run exists to answer.
2. **Polarization buildup rate**: is `<Sy>` vs. turn roughly consistent in
   scale/trend with the LCFA run's ~6.9e-6 %/turn rate, or does stochastic
   recoil measurably change it?
3. **Energy-spread growth**: compare `sigma_pz` vs. turn between the two
   runs' `.ave` output -- the paper's own text flags this as the mechanism
   stochastic recoil could affect that a mean kick cannot.

Record live particle counts at turn 0 and turn 100000 explicitly (as
`validate_output.py` does), not just final polarization, per the same
requested-vs-live distinction the original LCFA handoff emphasized.
