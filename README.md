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

The workflow is:

```text
Packmol spherical droplet
  -> solvent-only GFN2-xTB pre-relaxation
  -> 01_100K          1 ps  initial thermal preparation (not production)
  -> 02_200K          2 ps  intermediate heating (not production)
  -> 03_298K_equil    3 ps  final-temperature equilibration
  -> 04_298K_screen   5 ps  structural screening at 298.15 K
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
    │   └── stages/
    │       └── 00_relax/
    │           ├── 00_relax.out
    │           ├── xtbopt.*
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

The times in ps are the source of truth; step counts above are derived from the
0.5 fs timestep. GFN2-xTB remains the default, with physical H masses and SHAKE
off.

## Scope boundary

This script is deliberately limited to E2 preparation/screening.

A later analysis-oriented program can consume the stable directory layout and
control TRAVIS, coordination analysis, clustering and E2->E3/E4 selection without
mixing simulation execution with analysis logic.
