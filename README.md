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
validations. The default `--thermostat-warning-policy ramp` accepts an isolated
`thermostating problem` warning in `01_100K` and `02_200K`, after all output
integrity checks pass, but remains strict in `03_298K_equil` and
`04_298K_screen`. Use `strict` to reject it everywhere or `allow` to accept it
at any stage after the same checks. Accepted warnings remain explicit in
`stage_manifest.json`. With `--force`, the old status marker and archive are
explicitly invalidated before the new attempt begins.
For MD, a completed archive must contain `xtb.trj`, `mdrestart`, and `xtbmdok`.
Restarted stages must replace the input `mdrestart` with a byte-different output.
Each MD archive records its inputs and restart hashes in `stage_manifest.json`;
therefore, rerunning an earlier stage can require later stages to be rerun when
their recorded input-restart hash no longer matches.

To promote an existing `stage.failed` whose sole reason is
`thermostating problem`, without rerunning xTB, use
`--resume-thermostat-warning`. The archive, log, required outputs, restart
chain, and current warning policy are revalidated before `stage.done` is
created. Other failure reasons are never promoted.

The workflow is:

```text
Packmol spherical droplet
  -> solvent-only GFN2-xTB pre-relaxation
  -> 100 K
  -> 200 K
  -> 298 K equilibration
  -> 298 K screening
```

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
- 100 K, 0.5 ps
- 200 K, 0.5 ps
- 298.15 K equilibration, 1.0 ps
- 298.15 K screening, 5.0 ps

GFN2-xTB by default, 0.5 fs timestep, physical H masses, SHAKE off.

## Scope boundary

This script is deliberately limited to E2 preparation/screening.

A later analysis-oriented program can consume the stable directory layout and
control TRAVIS, coordination analysis, clustering and E2->E3/E4 selection without
mixing simulation execution with analysis logic.
