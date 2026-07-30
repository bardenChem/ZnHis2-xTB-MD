# AGENTS.md

## Scope

These instructions apply repository-wide unless a more specific `AGENTS.md` exists in a subdirectory.

This repository supports computational modeling of the biomimetic Zn(His)2 system in explicit water, with emphasis on xTB-based finite-droplet molecular dynamics, structural screening, coordination-state analysis, and preparation for later higher-level calculations (e.g. CP2K/DFT). The current pipeline is an E2 screening workflow, not a converged bulk-water production protocol.

The agent should optimize for **scientific traceability, reproducibility, conservative automation, and protection against silent reuse of incompatible results**.

---

## 1. General working rules

Before editing code:

1. Read the relevant file(s) completely enough to understand the existing control flow.
2. Prefer surgical changes over rewrites.
3. Preserve the existing architecture unless a redesign is explicitly requested.
4. Do not make unrelated cleanup/refactoring in the same change.
5. Do not silently change scientific parameters, defaults, file formats, stage ordering, or interpretation rules.
6. Distinguish clearly between:
   - a software bug;
   - a robustness/provenance improvement;
   - a scientific-method change.
7. If a requested fix would require changing scientific assumptions, stop and explain the consequence before implementing it.

Do not invent behavior of xTB, Packmol, CP2K, or other external programs. If a change depends on undocumented behavior, state the uncertainty and prefer a conservative implementation.

---

## 2. Git and repository safety

Unless explicitly requested by the user:

- do **not** commit;
- do **not** push;
- do **not** create or modify remotes;
- do **not** rewrite git history;
- do **not** delete user data or previous simulation results;
- do **not** run destructive commands on real production directories.

When substantial code changes are made, report:

```bash
git diff --stat
git diff -- <modified files>
```

If tests need temporary data, use `/tmp`, a dedicated test directory, or fixtures. Never test destructive behavior against `md_screening/` production data.

---

## 3. Computational cost policy

Do not launch long molecular dynamics, geometry optimizations, metadynamics, production trajectories, or other expensive calculations unless the user explicitly asks for execution.

By default, use only lightweight validation such as:

```bash
python3 -m py_compile xtb_md_pipeline.py
python3 xtb_md_pipeline.py --help
```

and unit/mock/preparation-only tests that do not invoke long xTB runs.

Short external-program smoke tests are acceptable only when explicitly useful and should be identified as smoke tests, not scientific production calculations.

---

## 4. Current scientific workflow

The current E2 workflow is conceptually:

```text
unsolvated Zn(His)2 structure
        ↓
Packmol spherical explicit-water droplet
        ↓
system_centered.pdb
        ↓
00_relax
  xTB optimization
  Zn(His)2 fixed
  solvent mobile
  same spherical log-Fermi wall
        ↓
system_relaxed.pdb
        ↓
01_100K
        ↓ mdrestart
02_200K
        ↓ mdrestart
03_298K_equil
        ↓ mdrestart
04_298K_screen
```

The primary structural systems are:

- `T4_cristal`
- `O6_I`
- `O6_II`
- `O6_III`

The literature tetrahedral structure `T4_IV_Hudecova` is a control.

Do not rename these identifiers unless explicitly asked.

---

## 5. Scientific defaults that must not be changed implicitly

Unless the user explicitly requests a methodological change, preserve the current defaults and stage definitions, including:

- GFN2-xTB as the default Hamiltonian;
- total charge `0`;
- UHF/unpaired-electron setting `0`;
- explicit water;
- no ALPB by default;
- finite spherical droplet rather than PBC;
- spherical log-Fermi wall;
- Packmol spherical packing;
- default Packmol sphere radius `12 Å`;
- wall margin `0.75 Å`;
- two replicas by default;
- `dt = 0.5 fs`;
- `dump = 10 fs`;
- `hmass = 1`;
- SHAKE off;
- `sccacc = 1.0`;
- 100 K / 0.5 ps;
- 200 K / 0.5 ps;
- 298.15 K / 1.0 ps equilibration;
- 298.15 K / 5.0 ps screening;
- solvent-only pre-relaxation with Zn(His)2 fixed;
- default relaxation level `loose`;
- current default relaxation-cycle count unless the user explicitly changes it.

A software-maintenance task must not silently become a protocol redesign.

---

## 6. Interpretation limits of E2

The E2 workflow is a **finite-droplet structural screening protocol**.

Do not describe its outputs as converged bulk equilibrium properties unless later evidence supports that claim.

In particular, E2 alone must not be used to claim quantitative:

- equilibrium populations;
- residence times;
- kinetic rate constants;
- free energies;
- bulk concentration effects;
- thermodynamic stability rankings.

Suitable E2 outputs include qualitative/semiquantitative assessment of:

- Zn coordination changes;
- Zn–N and Zn–O distances;
- carboxylate coordination/decoordination;
- water entry into the Zn first coordination shell;
- persistence of structural basins;
- differences between replicas;
- candidate states for later periodic or higher-level calculations.

Do not describe constrained optimizations or relaxed scans as free-energy calculations.

---

## 7. Provenance is a first-class requirement

Never silently reuse a computational result when its generating configuration cannot be verified.

### Packing provenance

A reusable `packed_sphere.pdb` should be tied to a packing signature containing at least:

- solute SHA256;
- water-template SHA256;
- sphere radius;
- water count;
- Packmol tolerance;
- Packmol seed.

If an existing packed structure has no provenance metadata, do not assume compatibility. Require explicit regeneration (`--repack`) or user intervention.

Changing a packing parameter must not cause a new manifest to describe an old packing.

### Relaxation provenance

A reusable `00_relax` result should retain a configuration signature including at least:

- input geometry SHA256;
- GFN method;
- charge;
- UHF;
- implicit-solvent setting;
- optimization level;
- requested maximum cycles;
- wall radius;
- fixed-atom selection.

Requested optimization cycles and cycles reported by xTB are distinct provenance values and should not be assumed equal.

### MD provenance

Each MD stage should retain a stage manifest/signature sufficient to determine whether a completed stage is compatible with the current configuration and restart chain. Relevant fields include:

- stage name;
- GFN;
- charge;
- UHF;
- ALPB setting;
- thread count as execution provenance;
- target temperature;
- simulation time;
- timestep;
- dump interval;
- hydrogen mass;
- SHAKE setting;
- SCC accuracy;
- restart setting;
- wall radius;
- input geometry hash;
- xTB input hash;
- input restart hash where applicable;
- output restart hash.

Do not reuse an old downstream stage if its recorded input restart does not correspond to the currently valid upstream stage.

---

## 8. Stage-status semantics

Status markers must have unambiguous meaning.

### `stage.done`

`stage.done` means that the stage:

1. completed numerically;
2. produced required outputs;
3. passed workflow validation;
4. is compatible with its recorded input configuration.

It must be written **only after all validation succeeds**.

### `stage.failed`

`stage.failed` may be used when outputs exist but the stage failed workflow validation.

A failed stage must never be treated as completed.

### Archives

Archiving outputs and marking a stage as valid are separate operations.

An archive function must not implicitly certify success.

If `--force` re-runs a stage, invalidate/remove old completion/failure markers before starting the new attempt. A failed forced rerun must not leave a stale `stage.done` that certifies the old result.

Do not automatically delete later stages merely because an upstream stage changed, but downstream reuse must fail if provenance/restart hashes show incompatibility.

---

## 9. MD success criteria

A zero xTB return code alone is not sufficient to approve a molecular-dynamics stage.

Before writing `stage.done`, validate at minimum:

- return code is zero;
- no configured fatal pattern occurs in the log;
- no `thermostating problem` occurs;
- `xtb.trj` exists;
- `mdrestart` exists;
- for restart stages, output `mdrestart` is demonstrably new relative to the input restart (prefer SHA256, not timestamp);
- `xtbmdok` exists;
- stage provenance has been written successfully.

If a restart-stage output `mdrestart` is byte-identical to its input restart, do not approve the stage.

Preserve diagnostic outputs when a stage fails validation.

---

## 10. `thermostating problem`

For this project, xTB reporting:

```text
thermostating problem
```

is treated as a scientific validation failure of the MD stage, even when xTB exits numerically without a crash.

Required behavior:

- archive available outputs;
- mark failure;
- do not create `stage.done`;
- do not continue automatically to the next temperature stage.

The user may inspect the trajectory and later choose whether to modify preparation or thermalization.

---

## 11. Restart integrity

The restart chain is scientifically important.

For stages with `restart=true`:

1. retrieve the validated `mdrestart` from the previous stage;
2. record/hash it as the input restart;
3. execute the new stage;
4. ensure a new output restart was produced;
5. record/hash the output restart;
6. archive both input provenance and output provenance where practical.

Existence alone is insufficient because the input `mdrestart` already exists before xTB starts.

Do not use file modification time as the primary integrity test.

---

## 12. `00_relax` behavior

The pre-relaxation is intended to reduce pathological Packmol contacts and solvent strain before MD, not necessarily to locate a rigorously converged 0 K minimum of the full droplet.

Maintain:

```text
$fix
 atoms: 1-N_SOLUTE
$end
```

with solvent atoms mobile.

The code should validate that the fixed solute did not move beyond the defined numerical tolerance.

Formal optimization non-convergence may be a warning rather than a fatal failure if:

- xTB returned normally;
- a valid final geometry exists;
- there are no fatal log patterns;
- fixed solute coordinates remain fixed within tolerance.

Do not report `converged=None` as confirmed convergence. Distinguish:

- `True`: convergence detected;
- `False`: non-convergence detected;
- `None`: convergence status could not be determined robustly.

Parsers should prefer `None` to guessing.

---

## 13. xTB-specific conventions

The reference production version for current reproducibility is **xTB 6.7.1** unless the user explicitly changes versions.

Do not describe it as the globally latest version without checking current official sources.

Supported optimization-level names should follow xTB syntax:

- `crude`
- `sloppy`
- `loose`
- `lax`
- `normal`
- `tight`
- `vtight`
- `extreme`

Do not use the invalid alias `verytight`.

Use `mdrestart` consistently; avoid the typo `mdrrestart` in messages and documentation.

When parsing xTB output:

- tolerate minor formatting variation;
- use case-insensitive regexes where appropriate;
- do not infer energies or cycle counts from ambiguous lines;
- retain `None` when a value cannot be extracted reliably.

---

## 14. Shell/environment concerns

On systems where xTB requires a larger process stack, the execution environment may need:

```bash
ulimit -s unlimited
```

and OpenMP settings such as:

```bash
export OMP_STACKSIZE=4G
```

Do not add Python `resource.setrlimit()` logic or other automatic OS-level stack modifications unless explicitly requested.

Do not hardcode host-specific installation paths such as `/home/.../xtb-dist/bin/xtb` into the pipeline. Use CLI configuration/PATH and preserve portability.

Thread count should remain configurable rather than hardcoded globally.

---

## 15. Output-file preservation

Production outputs are scientific data.

Before deleting, overwriting, or replacing simulation outputs:

- know which stage they belong to;
- know whether they are archived;
- preserve provenance;
- prefer explicit `--force`/`--repack` behavior over silent replacement.

A failed run may still contain diagnostically valuable trajectories, logs, restart files, or snapshots. Archive them where feasible before raising an error.

---

## 16. Coding style

Use Python standard library unless a new dependency is clearly necessary and explicitly justified.

Prefer:

- small functions with one responsibility;
- `pathlib.Path` for filesystem paths;
- explicit validation;
- deterministic metadata;
- SHA256 for content identity;
- clear exceptions with actionable messages;
- conservative parsing;
- idempotent preparation when inputs/configuration are unchanged.

Avoid:

- broad exception swallowing;
- silent fallback to incompatible data;
- hidden mutation of scientific parameters;
- duplicated stage-validation logic when a helper function is clearer;
- rewriting working code merely for style.

---

## 17. Required lightweight checks after edits

At minimum run:

```bash
python3 -m py_compile xtb_md_pipeline.py
```

For changes affecting the CLI, also inspect/test:

```bash
python3 xtb_md_pipeline.py --help
```

For changes affecting packing, provenance, status markers, restart logic, or parsers, add lightweight mock/unit tests where practical. Tests should cover failure paths as well as success paths.

Examples of expected test cases:

- identical input/output restart rejected;
- changed restart accepted;
- missing `xtbmdok` rejected;
- `stage.done` with incomplete archive rejected;
- `stage.failed` not reused;
- compatible MD stage manifest reused;
- incompatible temperature/time/input-restart hash rejected;
- Packmol signature mismatch requires `--repack`;
- non-converged xTB optimization log parsed conservatively;
- unknown convergence state remains `None`.

Do not claim tests passed unless they were actually executed.

---

## 18. Response/reporting after code changes

When reporting a completed coding task, summarize:

1. what problem was addressed;
2. what files changed;
3. what behavior changed;
4. what scientific behavior was intentionally left unchanged;
5. tests run and their actual results;
6. unresolved uncertainties or assumptions;
7. whether any external calculation was run;
8. whether any commit/push occurred.

For substantial edits, include or summarize:

```bash
git diff --stat
```

Do not claim that the workflow is scientifically validated solely because software tests pass. Software validation and scientific validation are separate.

---

## 19. Scientific review mindset

When code behavior affects interpretation of molecular simulations, reason conservatively.

Keep separate:

- reported xTB output;
- program-level validation;
- scientific interpretation;
- hypothesis for later testing.

Examples:

- `normal termination` means the program returned normally, not necessarily that the MD is physically acceptable;
- `stage.done` means workflow validation passed, not that equilibrium was reached;
- a stable coordination state over a short E2 trajectory is a screening observation, not an equilibrium population;
- semiempirical MD describes the chosen semiempirical Hamiltonian and should not be presented as DFT-level evidence without benchmarking.

When uncertain, preserve the data and surface the uncertainty instead of automating a strong conclusion.

---

## 20. Default agent behavior for future tasks

If the user asks to “fix”, “improve”, “review”, or “add” something to this repository without more detailed instructions:

1. inspect the current implementation;
2. identify the smallest correct change;
3. preserve the scientific protocol;
4. preserve provenance and stage integrity;
5. avoid expensive calculations;
6. add lightweight tests appropriate to the change;
7. show the relevant diff/result;
8. do not commit or push unless explicitly asked.

These defaults should be assumed so they do not need to be repeated in every prompt.
