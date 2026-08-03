# xTB MD pipeline v2 — E2 screening

This revision is aligned with the current Zn(His)2 project directory tree.

## Main changes

- Two independent replicas per system by default.
- Each replica receives its own Packmol spherical water droplet.
- Packmol uses classic `inside sphere` syntax, so no `pbc` keyword is required.
- Independent Packmol seeds are stored in each `manifest.json`.
- The full droplet is centered at its center of mass before applying the xTB wall.
- The xTB wall radius is based on the actual outermost atom after centering.
- A short solvent-only GFN2-xTB pre-relaxation precedes MD by default.
- MD logs are checked for explicit instability/emergency-exit patterns.
- `xtb-trj.pdb` is archived when generated.
- Output root defaults to `md_screening/`.

## Expected source tree

```text
.
├── xtb_md_pipeline.py
├── O6_I_box/
│   ├── water.pdb
│   └── ZnHis2_I.pdb
├── O6_II_box/
│   └── ZnHis2_II.pdb
├── O6_III_box/
│   └── ZnHis2_III.pdb
├── T4_box_cristal/
│   └── cristal.pdb
└── T4_IV_box/
    └── ZnHis2_IV.pdb
```

The program finds `water.pdb` automatically from the current project folders.
You may instead pass `--water-pdb /path/to/water.pdb`.

## Prepare the four main systems

```bash
python3 xtb_md_pipeline.py --main
```

Default E2 droplet:

- spherical Packmol droplet;
- radius = 12 A;
- water count estimated to give an initial total mass density near 1.0 g/cm3;
- 2 independent replicas per system;
- Packmol tolerance = 2 A.

The current Zn(His)2 composition should give roughly ~220 waters for a 12 A droplet.
The exact value is computed from the actual PDB composition.

## Inspect before MD

For every replica inspect:

```text
md_screening/<SYSTEM>/replica_XX/packing/packed_sphere.pdb
md_screening/<SYSTEM>/replica_XX/system_centered.pdb
md_screening/<SYSTEM>/replica_XX/manifest.json
md_screening/<SYSTEM>/replica_XX/00_relax.inp
```

## Run after inspection

```bash
python3 xtb_md_pipeline.py --main --threads 8 --run
```

Completed stages are skipped on later invocations unless `--force` is used.
The `stage.done` marker means that a stage finished and passed all output/log
validations. The default `--thermostat-warning-policy allow` treats an isolated
`thermostating problem` as a warning, at any MD stage, only after xTB returns
zero, the log reports `normal exit of md()`, no fatal pattern is present, and
`xtb.trj`, a new valid `mdrestart`, and `xtbmdok` all pass the integrity checks.
The warning remains visible in the terminal and in
`stage_manifest.json` (`thermal_result.thermostating_problem=true` and
`warning_accepted=true`). It never bypasses output, restart, hash, or provenance
validation. The revised protocol always canonicalizes the policy to `allow`, so
an isolated warning cannot block 03 -> 04 after those checks pass. Legacy CLI
values `strict` and `ramp` are accepted with a deprecation warning and mapped to
`allow`. With `--force`, the old status marker and archive are explicitly
invalidated before the new attempt begins.
For MD, a completed archive must contain `xtb.trj`, `mdrestart`, and `xtbmdok`.
Restarted stages must replace the input `mdrestart` with a byte-different output.
Each MD archive records its inputs and restart hashes in `stage_manifest.json`;
therefore, rerunning an earlier stage can require later stages to be rerun when
their recorded input-restart hash no longer matches.

To resume a calculation prepared with the current protocol, without Packmol,
recentering, input regeneration, or a repeated `00_relax`, use:

```bash
python3 xtb_md_pipeline.py --main --replicas 2 --threads 8 --run --resume
```

`--resume` validates `stage.done`/`stage.failed`, `stage_manifest.json`, the
configuration and input hashes, the archived log, `xtb.trj`, `mdrestart`,
`mdrestart.input`, `xtbmdok`, and restart continuity. Compatible stages print
`REUSE`; the first absent stage prints `RUN`. A `stage.failed` whose only reason
is `thermostating problem` is promoted without executing xTB only if all the
new checks pass. For a legacy thermostat-only failure created before a stage
manifest was written, configuration provenance is reconstructed from the
untouched replica `manifest.json`, the deterministic historical stage input,
and the archived restart/output chain; this inference is recorded explicitly
under `recovery` in the new stage manifest. A crash, fatal log, missing
output/restart, unchanged or malformed restart, hash inconsistency, or
incompatible provenance is never promoted.

The new stage durations are part of each stage signature, so normal `--resume`
does not silently reuse an old 0.5/0.5/1.0 ps thermalization as the new
1.0/2.0/3.0 ps protocol. To explicitly run only the 5 ps screening from an
existing, valid historical `03_298K_equil` restart, use:

```bash
python3 xtb_md_pipeline.py \
    --main \
    --project md_screening \
    --replicas 2 \
    --threads 8 \
    --gfn 2 \
    --charge 0 \
    --uhf 0 \
    --run \
    --start-stage 04_298K_screen
```

The command above is exact for a default `md_screening` project prepared with
the four main systems, two replicas, eight threads, GFN2, charge 0, UHF 0, and
ALPB off. For another existing project, repeat its exact system selector,
`--project`, replica count, thread count, GFN, charge, UHF, and ALPB setting;
the provenance checks abort instead of guessing. This existing-only mode does
not call Packmol, rebuild the droplet, rerun `00_relax`, or rerun 01/02/03. It
validates the historical chain, safely promotes a thermostat-warning-only 03
when possible, preserves its recorded historical duration, and prints an
explicit warning that the restart came from an earlier duration configuration.
Do not add `--force` or `--repack`.

## Optional 20 ps continuation at 298.15 K

`05_298K_extended` is an opt-in continuation of the validated final
`mdrestart` from `04_298K_screen`. Normal preparation, `--run`, and `--resume`
still stop at stage 04; they do not launch the extra 20 ps automatically. The
continuation keeps the stage-04 Hamiltonian, charge/UHF, NVT settings, physical
H masses, SHAKE setting, 0.5 fs timestep, 10 fs dump interval, and the same
spherical wall. It uses `restart=true`, so coordinates and velocities come
directly from `stages/04_298K_screen/mdrestart`; no packing, optimization,
temperature ramp, velocity randomization, or solvation rebuild is performed.

For the currently completed `replica_01` data, whose historical manifests
record 25 threads, the new stage may use a different allocation. For example,
validate a 32-thread continuation without writing files or calling xTB:

```bash
python3 xtb_md_pipeline.py \
    --system O6_I T4_cristal T4_IV_Hudecova \
    --project md_screening \
    --replicas 1 \
    --threads 32 \
    --start-stage 05_298K_extended \
    --dry-run
```

Execute the extension only after inspecting the dry-run:

```bash
python3 xtb_md_pipeline.py \
    --system O6_I T4_cristal T4_IV_Hudecova \
    --project md_screening \
    --replicas 1 \
    --threads 32 \
    --start-stage 05_298K_extended \
    --run
```

Individual commands use the same interface:

```bash
python3 xtb_md_pipeline.py --system O6_I --project md_screening \
    --replicas 1 --threads 32 --start-stage 05_298K_extended --run
python3 xtb_md_pipeline.py --system T4_cristal --project md_screening \
    --replicas 1 --threads 32 --start-stage 05_298K_extended --run
python3 xtb_md_pipeline.py --system T4_IV_Hudecova --project md_screening \
    --replicas 1 --threads 32 --start-stage 05_298K_extended --run
```

`--threads` applies to the new xTB execution and may differ from stages 01--04.
Both the historical and newly requested counts remain recorded as execution
provenance, but thread count is not a scientific compatibility condition. All
physical settings must still match; those mismatches abort. A valid stage-05
`stage.done` is reused by default. `stage.failed`, `stage.running`, or an
unmarked nonempty archive aborts without overwriting it. To deliberately retry
only stage 05 while retaining stage 04, add `--force` to the stage-05 command.
The preflight also validates the complete 01--04 chain, restart coordinates and
velocities, composition, hashes, and a conservative free-disk estimate.

The workflow is:

```text
Packmol spherical droplet
  -> solvent-only GFN2-xTB pre-relaxation
  -> 01_100K          1 ps  initial thermal preparation (not production)
  -> 02_200K          2 ps  intermediate heating (not production)
  -> 03_298K_equil    3 ps  final-temperature equilibration
  -> 04_298K_screen   5 ps  structural screening at 298.15 K
  -> 05_298K_extended 20 ps optional direct continuation at 298.15 K
```

`04_298K_screen` has `restart=true` and receives positions and velocities from
the final `mdrestart` of `03_298K_equil`. It does not start from
`system_relaxed.pdb` and does not generate new random velocities.

During `00_relax`, all atoms belonging to Zn(His)2 are fixed and only the
explicit water molecules can move. The atom range is derived from the
unsolvated solute PDB and validated against the beginning of Packmol's output.
This removes artificial contacts and strain from initial packing without
allowing a zero-temperature optimization to change the coordination state that
the MD is intended to test. The same spherical log-Fermi wall used by MD remains
active. Defaults are `--opt loose`, 30 cycles, and no ALPB.

The native xTB optimization engine for `00_relax` can be selected with
`--relax-engine auto|rf|lbfgs|inertial`. The default `auto` writes no `$opt`
engine override. `inertial` selects xTB's native FIRE engine while preserving
the fixed solute and spherical wall. Compare engines in a separate project:

```bash
python3 xtb_md_pipeline.py --system O6_I --replicas 1 \
    --project md_fire_test --relax-engine inertial --relax-cycles 150
```

To reproduce the earlier workflow and start directly from
`system_centered.pdb`, use:

```bash
python3 xtb_md_pipeline.py --system O6_I --skip-relax --run
```

## Run only one system

```bash
python3 xtb_md_pipeline.py --system O6_I
```

Then:

```bash
python3 xtb_md_pipeline.py --system O6_I --threads 8 --run
```

## Control structure

```bash
python3 xtb_md_pipeline.py --controls
python3 xtb_md_pipeline.py --controls --threads 8 --run
```

## Override droplet size

Example: 13 A:

```bash
python3 xtb_md_pipeline.py --main --sphere-radius 13.0
```

## Override water count

```bash
python3 xtb_md_pipeline.py --main --waters 220
```

If `--waters` is supplied, the density-based estimate is disabled.

## Force new Packmol configurations

```bash
python3 xtb_md_pipeline.py --main --repack
```

This replaces the existing `packed_sphere.pdb` for the selected replicas.
Do not do this after starting production runs unless intentionally rebuilding them.
An existing packing is reused only when `packing/packing_manifest.json` proves
that its solute, water template, radius, water count, tolerance, and seed match
the current request. Otherwise the pipeline stops and requires `--repack`.

## Output structure

```text
md_screening/
└── O6_I/
    ├── replica_01/
    │   ├── packing/
    │   │   ├── solute.pdb
    │   │   ├── water.pdb
    │   │   ├── packmol_sphere.inp
    │   │   ├── packmol.out
    │   │   ├── packing_manifest.json
    │   │   └── packed_sphere.pdb
    │   ├── system_centered.pdb
    │   ├── system_relaxed.pdb
    │   ├── manifest.json
    │   ├── 00_relax.inp
    │   ├── 01_100K.inp
    │   ├── 02_200K.inp
    │   ├── 03_298K_equil.inp
    │   ├── 04_298K_screen.inp
    │   ├── 05_298K_extended.inp
    │   └── stages/
    │       ├── 00_relax/
    │       │   ├── 00_relax.out
    │       │   ├── xtbopt.*
    │       │   └── stage.done
    │       └── 05_298K_extended/       # only when explicitly requested
    │           ├── 05_298K_extended.inp
    │           ├── 05_298K_extended.out
    │           ├── xtb.trj
    │           ├── mdrestart.input
    │           ├── mdrestart
    │           ├── stage_manifest.json
    │           └── stage.done
    └── replica_02/
```

## E2 protocol

Per replica:

- solvent-only pre-relaxation, loose, at most 30 cycles
- initial thermal preparation at 100 K, 1.0 ps (2000 steps; not production)
- intermediate heating at 200 K, 2.0 ps (4000 steps; not production)
- 298.15 K equilibration, 3.0 ps (6000 steps)
- 298.15 K structural screening, 5.0 ps (10000 steps)
- optional 298.15 K extension, 20.0 ps (40000 steps; approximately 2000
  frames, allowing the xTB initial/final-frame convention)

The times in ps are the source of truth; step counts above are derived from the
0.5 fs timestep. GFN2-xTB remains the default, with physical H masses and SHAKE
off.

## Trajectory analysis and diagnostic plots

`xtb_analysis.py` reads archived `xtb.trj` files, calculates structural
descriptors, writes reusable CSV/JSON tables, and generates headless Matplotlib
plots. Plotting is enabled by default.

The trajectory parser accepts both coordinate-only extended XYZ frames and the
xTB layout in which each coordinate block is followed by one three-component
velocity record per atom. Velocities are preserved in `Frame.velocities`, while
`frames.csv` and `analysis_metadata.json` record their presence. They are not
used in the current structural analyses, no physical unit is assigned, and
all TRAVIS/VMD exports contain coordinates only.

Analyze the default screening stage of one replica:

```bash
python3 xtb_analysis.py \
    --replica md_screening/O6_I/replica_01 \
    --stage 04_298K_screen \
    --discard-first-ps 1.0 \
    --output-dir analysis/O6_I
```

Analyze stages 04 and 05 on one continuous output time axis:

```bash
python3 xtb_analysis.py \
    --replica md_screening/O6_I/replica_01 \
    --stage 04_298K_screen \
    --stage 05_298K_extended \
    --discard-first-ps 1.0 \
    --n-solute 39 \
    --output-dir analysis/O6_I_25ps
```

The current analyzer applies `--discard-first-ps` independently to every input
trajectory. Thus the command above discards 1 ps from stage 04 **and** 1 ps
from stage 05, even though the retained frames receive a continuous time axis.
There is currently no option to discard only from the first stage; this
behavior was documented here and was not changed as part of the MD extension.

Run numerical analysis without loading Matplotlib or creating graphics:

```bash
python3 xtb_analysis.py \
    --replica md_screening/O6_I/replica_01 \
    --stage 04_298K_screen \
    --discard-first-ps 1.0 \
    --output-dir analysis/O6_I \
    --no-plots
```

Named coordination sites use one-based atom indices. The values below are
placeholders: replace every `<INDEX>` after checking the replica PDB or a first
run's `atoms.csv`.

```bash
python3 xtb_analysis.py \
    --replica md_screening/O6_I/replica_01 \
    --stage 04_298K_screen \
    --discard-first-ps 1.0 \
    --site Nam1:<INDEX>:N \
    --site Ndelta1:<INDEX>:N \
    --site Ocarb1:<INDEX>:Ocarb \
    --site Nam2:<INDEX>:N \
    --site Ndelta2:<INDEX>:N \
    --site Ocarb2:<INDEX>:Ocarb \
    --output-dir analysis/O6_I
```

The plotting options are:

```text
--no-plots
--plot-format {png,pdf,svg}   default: png
--plot-dpi 300                PNG resolution
--plot-dir PATH               default: <output-dir>/plots
```

Coordinate export is enabled unless `--no-travis-export` is supplied. The
default alignment uses all non-hydrogen solute atoms. Export selection can be
controlled independently of the analysis tables:

```text
--alignment-selection {solute-heavy,solute-all}
--alignment-indices 1-5,10,12-19
--export-stride N
--export-start-ps START
--export-end-ps END
--no-travis-export
```

`--alignment-indices` overrides the named alignment selection and uses one-based
indices. At least three non-collinear atoms are required. The alignment
reference is the first analyzed frame; one Kabsch transform fitted on the
selected atoms is applied unchanged to every atom, including water. The
start/end window is inclusive, and `--export-stride` is applied after the main
analysis stride to frames inside that window. It does not remove rows from
`frames.csv`.

Outputs are conditional on the available data and requested analyses:

```text
analysis/O6_I/
├── atoms.csv
├── frames.csv
├── distance_summary.csv
├── Zn_Owater_radial_number.csv
├── coordination_states.csv
├── water_contact_events.csv
├── solute_RMSF.csv
├── analysis_metadata.json
├── trajectory_for_travis.xyz
├── travis/
│   ├── trajectory_full_coordinates.xyz
│   ├── trajectory_full_aligned_on_solute.xyz
│   ├── trajectory_solute_coordinates.xyz
│   ├── trajectory_solute_aligned.xyz
│   ├── reference_topology.pdb
│   ├── atom_map.csv
│   ├── travis_export_metadata.json
│   └── README_TRAVIS.md
├── vmd/
│   ├── convert_xyz_to_dcd.tcl
│   └── README_VMD.md
└── plots/
    ├── coordination_distances_vs_time.png
    ├── Zn_<group>_distances_vs_time.png
    ├── nearest_water_vs_time.png
    ├── coordination_number_vs_time.png
    ├── smooth_coordination_number_vs_time.png
    ├── coordination_state_vs_time.png
    ├── coordination_state_fraction.png
    ├── Zn_Owater_shell_count.png
    ├── Zn_Owater_cumulative_N.png
    ├── solute_RMSD_vs_time.png
    ├── solute_RMSF.png
    ├── tetrahedrality_vs_time.png
    ├── energy_vs_time.png
    ├── gnorm_vs_time.png
    └── distributions/
        ├── Zn_<site>_distance_histogram.png
        ├── nearest_water_distance_histogram.png
        ├── CN_smooth_histogram.png
        └── tetrahedrality_histogram.png
```

The four XYZ variants have distinct roles:

- `full_coordinates`: all atoms without transformation;
- `full_aligned_on_solute`: all atoms transformed using the solute-derived fit,
  suitable for solvent-relative inspection, SDF/TDO preparation, and videos;
- `solute_coordinates`: only the fixed first `n_solute` atoms, unaligned;
- `solute_aligned`: only the solute after the same configured fit, suitable for
  TDO and clean complex visualization.

The root `trajectory_for_travis.xyz` remains a byte-identical compatibility copy
of `travis/trajectory_full_coordinates.xyz`. Solute-dependent files are omitted
when `n_solute` is unknown, and aligned files are omitted when a valid alignment
selection cannot be constructed. `reference_topology.pdb` is copied from the
analysis topology when available; otherwise a minimal PDB is generated without
`CONECT` records. The generated VMD Tcl helper is never executed automatically.

Hard coordination numbers and states are calculated only when the user supplies
`--contact-cutoff-A`; the program does not assume a universal chemical cutoff.
The radial outputs are finite-droplet shell counts and cumulative `N(r)`, not a
bulk-normalized RDF or `g(r)`. State fractions and histograms from these short
screening trajectories are descriptive, not converged equilibrium populations.
See [`xtb_analysis_USAGE.md`](xtb_analysis_USAGE.md) for the remaining analysis
options and multi-stage examples.

## Scope boundary

`xtb_md_pipeline.py` is deliberately limited to E2 preparation/screening.

`xtb_analysis.py` consumes the stable directory layout without mixing simulation
execution with analysis logic. More advanced analysis, clustering, and E2->E3/E4
selection remain future work intended for integration with OOCCuPy.
