#!/usr/bin/env python3
"""
xtb_md_pipeline.py

Pipeline for E2 screening of Zn(His)2 + explicit water with xTB.

Current scope
-------------
1. Selects the structural systems defined in SYSTEMS.
2. Creates TWO independent replicas by default.
3. For each replica, builds a spherical explicit-water droplet with Packmol.
4. Centers the packed droplet at its total center of mass.
5. Creates xTB MD inputs for:
      01_100K        0.5 ps
      02_200K        0.5 ps
      03_298K_equil  1.0 ps
      04_298K_screen 5.0 ps
6. Optionally runs xTB sequentially, chaining stages through mdrrestart.
7. Archives trajectories, logs, restart files and snapshots per stage.

Default Hamiltonian
-------------------
GFN2-xTB, charge 0, UHF 0, explicit H2O, no implicit solvent.

Important
---------
This is an E2 finite-droplet screening workflow, NOT a periodic bulk-water MD.
A spherical log-Fermi wall is used to confine the droplet.

The Packmol droplet is created with `inside sphere`, which is compatible with
older Packmol versions that do not support the newer `pbc` keyword.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Iterable

BOHR_PER_ANGSTROM = 1.8897261254578281
AVOGADRO = 6.02214076e23
WATER_MOLAR_MASS = 18.01528  # g/mol

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# Project systems: UNSOLVATED PDB structures.
# These paths match the directory tree supplied for this project.
# ---------------------------------------------------------------------------

SYSTEMS = {
    "T4_cristal": {
        "solute": ROOT / "T4_box_cristal" / "cristal.pdb",
    },
    "T4_IV_Hudecova": {
        "solute": ROOT / "T4_IV_box" / "ZnHis2_IV.pdb",
    },
    "O6_I": {
        "solute": ROOT / "O6_I_box" / "ZnHis2_I.pdb",
    },
    "O6_II": {
        "solute": ROOT / "O6_II_box" / "ZnHis2_II.pdb",
    },
    "O6_III": {
        "solute": ROOT / "O6_III_box" / "ZnHis2_III.pdb",
    },
}

MAIN_SYSTEMS = [
    "T4_cristal",
    "O6_I",
    "O6_II",
    "O6_III",
]

CONTROL_SYSTEMS = [
    "T4_IV_Hudecova",
]

# Existing water templates in the user's current directory layout.
WATER_CANDIDATES = [
    ROOT / "water.pdb",
    ROOT / "O6_I_box" / "water.pdb",
    ROOT / "O6_II_box" / "water.pdb",
    ROOT / "O6_III_box" / "water.pdb",
    ROOT / "T4_box_cristal" / "water.pdb",
    ROOT / "T4_IV_box" / "water.pdb",
]

ATOMIC_MASSES = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "Zn": 65.38,
}

STAGES = [
    {
        "name": "01_100K",
        "temp": 100.0,
        "time": 0.5,
        "restart": False,
    },
    {
        "name": "02_200K",
        "temp": 200.0,
        "time": 0.5,
        "restart": True,
    },
    {
        "name": "03_298K_equil",
        "temp": 298.15,
        "time": 1.0,
        "restart": True,
    },
    {
        "name": "04_298K_screen",
        "temp": 298.15,
        "time": 5.0,
        "restart": True,
    },
]

FATAL_MD_PATTERNS = [
    "MD is unstable",
    "emergency exit",
    "Runtime exception",
    "segmentation fault",
    "floating point exception",
]


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def infer_element(line: str) -> str:
    """Infer element from a PDB ATOM/HETATM line."""
    if len(line) >= 78:
        element = line[76:78].strip()
        if element:
            element = element[0].upper() + element[1:].lower()
            if element in ATOMIC_MASSES:
                return element

    atom_name = line[12:16].strip()
    letters = re.sub(r"[^A-Za-z]", "", atom_name)
    if not letters:
        raise ValueError(f"Cannot infer element from PDB line:\n{line}")

    if letters[:2].lower() == "zn":
        return "Zn"

    element = letters[0].upper()
    if element not in ATOMIC_MASSES:
        raise ValueError(
            f"Unsupported/inferred element '{element}' from atom name "
            f"'{atom_name}'. Add it to ATOMIC_MASSES if needed."
        )
    return element


def read_pdb_atoms(lines: Iterable[str]):
    atoms = []

    for idx, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError as exc:
            raise ValueError(
                f"Invalid PDB coordinates at line {idx + 1}:\n{line}"
            ) from exc

        element = infer_element(line)

        atoms.append(
            {
                "line_index": idx,
                "element": element,
                "mass": ATOMIC_MASSES[element],
                "xyz": (x, y, z),
            }
        )

    if not atoms:
        raise ValueError("No ATOM/HETATM records found in PDB.")

    return atoms


def pdb_molar_mass(pdb: Path) -> float:
    """Return molecular molar mass in g/mol from PDB atom identities."""
    lines = pdb.read_text().splitlines(keepends=True)
    atoms = read_pdb_atoms(lines)
    return sum(a["mass"] for a in atoms)


def center_pdb(input_pdb: Path, output_pdb: Path):
    """
    Center full droplet at total center of mass.

    This is the reference used for the spherical xTB wall.
    Analysis can later recenter frames on Zn independently.
    """
    lines = input_pdb.read_text().splitlines(keepends=True)
    atoms = read_pdb_atoms(lines)

    total_mass = sum(a["mass"] for a in atoms)

    com = [
        sum(a["mass"] * a["xyz"][k] for a in atoms) / total_mass
        for k in range(3)
    ]

    atom_by_line = {a["line_index"]: a for a in atoms}
    centered_xyz = []
    out_lines = []

    for idx, line in enumerate(lines):
        if idx not in atom_by_line:
            out_lines.append(line)
            continue

        x, y, z = atom_by_line[idx]["xyz"]
        xc = x - com[0]
        yc = y - com[1]
        zc = z - com[2]

        centered_xyz.append((xc, yc, zc))

        core = line.rstrip("\r\n")
        newline = "\n"

        if len(core) < 54:
            core = core.ljust(54)

        new_line = (
            core[:30]
            + f"{xc:8.3f}"
            + f"{yc:8.3f}"
            + f"{zc:8.3f}"
            + core[54:]
            + newline
        )
        out_lines.append(new_line)

    output_pdb.write_text("".join(out_lines))

    rmax = max(
        math.sqrt(x * x + y * y + z * z)
        for x, y, z in centered_xyz
    )

    return {
        "n_atoms": len(atoms),
        "total_mass_u": total_mass,
        "original_center_of_mass_A": com,
        "max_radius_from_COM_A": rmax,
    }


# ---------------------------------------------------------------------------
# Packmol spherical droplet
# ---------------------------------------------------------------------------

def find_water_pdb(explicit: Path | None = None) -> Path:
    if explicit is not None:
        if not explicit.exists():
            raise SystemExit(f"Water PDB not found: {explicit}")
        return explicit.resolve()

    for candidate in WATER_CANDIDATES:
        if candidate.exists():
            return candidate.resolve()

    raise SystemExit(
        "No water.pdb found. Put water.pdb in the project root or use "
        "--water-pdb /path/to/water.pdb."
    )


def sphere_volume_A3(radius_A: float) -> float:
    return (4.0 / 3.0) * math.pi * radius_A ** 3


def estimate_water_count(
    solute_pdb: Path,
    radius_A: float,
    density_g_cm3: float,
) -> int:
    """
    Estimate number of waters from target TOTAL mass density of the initial
    finite droplet.

    This is only an initialization rule for E2. A finite droplet with a wall
    is not a thermodynamic bulk-density calculation.
    """
    volume_cm3 = sphere_volume_A3(radius_A) * 1.0e-24
    target_mass_equiv_g_mol = density_g_cm3 * volume_cm3 * AVOGADRO

    solute_mm = pdb_molar_mass(solute_pdb)
    remaining_mm = target_mass_equiv_g_mol - solute_mm

    if remaining_mm <= 0:
        raise ValueError(
            f"Sphere radius {radius_A:.2f} A is too small for target density "
            f"{density_g_cm3:.3f} g/cm3 and solute molar mass "
            f"{solute_mm:.2f} g/mol."
        )

    nwater = int(round(remaining_mm / WATER_MOLAR_MASS))

    if nwater < 1:
        raise ValueError("Estimated water count is < 1.")

    return nwater


def packmol_input(
    solute_name: str,
    water_name: str,
    output_name: str,
    radius_A: float,
    nwater: int,
    tolerance_A: float,
    seed: int,
) -> str:
    """
    Packmol input using classic `inside sphere` syntax.

    `center` + `fixed 0 0 0 0 0 0` places the solute geometric barycenter
    at the origin without rotation.
    """
    return f"""# Spherical explicit-water droplet for xTB E2 screening
# Independent Packmol seed: {seed}

tolerance {tolerance_A:.3f}
filetype pdb
output {output_name}
seed {seed}

structure {solute_name}
  number 1
  center
  fixed 0.0 0.0 0.0 0.0 0.0 0.0
end structure

structure {water_name}
  number {nwater}
  inside sphere 0.0 0.0 0.0 {radius_A:.3f}
end structure
"""


def run_packmol(
    replica_dir: Path,
    solute_pdb: Path,
    water_pdb: Path,
    args,
    seed: int,
):
    packing_dir = replica_dir / "packing"
    packing_dir.mkdir(parents=True, exist_ok=True)

    solute_local = packing_dir / "solute.pdb"
    water_local = packing_dir / "water.pdb"
    packed_pdb = packing_dir / "packed_sphere.pdb"
    packmol_inp = packing_dir / "packmol_sphere.inp"
    packmol_log = packing_dir / "packmol.out"

    shutil.copy2(solute_pdb, solute_local)
    shutil.copy2(water_pdb, water_local)

    if args.waters is None:
        nwater = estimate_water_count(
            solute_pdb=solute_pdb,
            radius_A=args.sphere_radius,
            density_g_cm3=args.density,
        )
    else:
        nwater = args.waters

    inp_text = packmol_input(
        solute_name=solute_local.name,
        water_name=water_local.name,
        output_name=packed_pdb.name,
        radius_A=args.sphere_radius,
        nwater=nwater,
        tolerance_A=args.packmol_tolerance,
        seed=seed,
    )
    packmol_inp.write_text(inp_text)

    if packed_pdb.exists() and not args.repack:
        return packed_pdb, nwater, packmol_inp, packmol_log

    packmol_exe = shutil.which(args.packmol)
    if packmol_exe is None:
        raise RuntimeError(
            f"Packmol executable '{args.packmol}' not found in PATH. "
            "Use --packmol /path/to/packmol if needed."
        )

    with packmol_inp.open("r") as fin, packmol_log.open("w") as fout:
        result = subprocess.run(
            [packmol_exe],
            cwd=packing_dir,
            stdin=fin,
            stdout=fout,
            stderr=subprocess.STDOUT,
            text=True,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"Packmol failed in {replica_dir}. See {packmol_log}"
        )

    if not packed_pdb.exists():
        raise RuntimeError(
            f"Packmol ended without creating {packed_pdb}. "
            f"See {packmol_log}"
        )

    log_text = packmol_log.read_text(errors="replace")

    # Packmol normally prints "Success!" on successful convergence.
    # We warn instead of failing because older builds may format output differently.
    if "Success!" not in log_text:
        print(
            f"  WARNING: Packmol output for {replica_dir.name} does not contain "
            f"'Success!'. Inspect {packmol_log}."
        )

    return packed_pdb, nwater, packmol_inp, packmol_log


# ---------------------------------------------------------------------------
# xTB input and run
# ---------------------------------------------------------------------------

def md_input(stage, wall_radius_bohr: float) -> str:
    restart = "true" if stage["restart"] else "false"

    return f"""$md
   temp={stage['temp']:.2f}
   time={stage['time']:.3f}
   dump=10.0
   step=0.5
   velo=true
   nvt=true
   hmass=1
   shake=0
   sccacc=1.0
   restart={restart}
$end

$wall
   potential=logfermi
   sphere: {wall_radius_bohr:.3f}, all
$end
"""


def write_manifest(
    replica_dir: Path,
    system_name: str,
    replica_index: int,
    source: Path,
    packed_pdb: Path,
    water_pdb: Path,
    nwater: int,
    packmol_seed: int,
    geom_info: dict,
    wall_radius_A: float,
    args,
):
    data = {
        "workflow_stage": "E2_screening",
        "system": system_name,
        "replica": replica_index,
        "source_solute_pdb": str(source.resolve()),
        "water_template_pdb": str(water_pdb.resolve()),
        "packed_pdb": str(packed_pdb.relative_to(replica_dir)),
        "centered_pdb": "system_centered.pdb",
        "packing": {
            "shape": "sphere",
            "sphere_radius_A": args.sphere_radius,
            "water_count": nwater,
            "target_initial_total_density_g_cm3": (
                None if args.waters is not None else args.density
            ),
            "packmol_tolerance_A": args.packmol_tolerance,
            "packmol_seed": packmol_seed,
        },
        "geometry_after_centering": {
            "n_atoms": geom_info["n_atoms"],
            "total_mass_u": geom_info["total_mass_u"],
            "original_center_of_mass_A": geom_info["original_center_of_mass_A"],
            "new_center_of_mass_A": [0.0, 0.0, 0.0],
            "max_radius_from_COM_A": geom_info["max_radius_from_COM_A"],
        },
        "wall": {
            "type": "logfermi_sphere",
            "margin_A": args.wall_margin,
            "radius_A": wall_radius_A,
            "radius_bohr": wall_radius_A * BOHR_PER_ANGSTROM,
        },
        "xtb": {
            "gfn": args.gfn,
            "charge": args.charge,
            "uhf": args.uhf,
            "alpb": args.alpb,
            "threads": args.threads,
        },
        "md": {
            "step_fs": 0.5,
            "dump_fs": 10.0,
            "hmass": 1,
            "shake": 0,
            "sccacc": 1.0,
            "stages": STAGES,
        },
        "interpretation_note": (
            "Finite explicit-water droplet screening. No PBC. "
            "Do not interpret E2 state frequencies as converged equilibrium populations."
        ),
    }

    (replica_dir / "manifest.json").write_text(
        json.dumps(data, indent=2)
    )


def prepare_replica(
    system_name: str,
    solute_pdb: Path,
    water_pdb: Path,
    replica_index: int,
    project_dir: Path,
    args,
) -> Path:
    replica_dir = (
        project_dir
        / system_name
        / f"replica_{replica_index:02d}"
    )
    replica_dir.mkdir(parents=True, exist_ok=True)

    packmol_seed = args.seed_base + replica_index

    packed_pdb, nwater, _, _ = run_packmol(
        replica_dir=replica_dir,
        solute_pdb=solute_pdb,
        water_pdb=water_pdb,
        args=args,
        seed=packmol_seed,
    )

    centered_pdb = replica_dir / "system_centered.pdb"
    geom_info = center_pdb(packed_pdb, centered_pdb)

    # The wall follows the ACTUAL centered droplet, not just the requested
    # Packmol radius. This avoids putting initial atoms inside the wall.
    wall_radius_A = (
        geom_info["max_radius_from_COM_A"] + args.wall_margin
    )
    wall_radius_bohr = wall_radius_A * BOHR_PER_ANGSTROM

    for stage in STAGES:
        (replica_dir / f"{stage['name']}.inp").write_text(
            md_input(stage, wall_radius_bohr)
        )

    write_manifest(
        replica_dir=replica_dir,
        system_name=system_name,
        replica_index=replica_index,
        source=solute_pdb,
        packed_pdb=packed_pdb,
        water_pdb=water_pdb,
        nwater=nwater,
        packmol_seed=packmol_seed,
        geom_info=geom_info,
        wall_radius_A=wall_radius_A,
        args=args,
    )

    print(f"\n[{system_name} / replica_{replica_index:02d}]")
    print(f"  solute             : {solute_pdb}")
    print(f"  droplet radius (A) : {args.sphere_radius:.3f}")
    print(f"  waters             : {nwater}")
    print(f"  Packmol seed       : {packmol_seed}")
    print(f"  atoms total        : {geom_info['n_atoms']}")
    print(
        "  original COM (A)   : "
        + " ".join(
            f"{v:8.3f}"
            for v in geom_info["original_center_of_mass_A"]
        )
    )
    print(
        f"  r_max centered (A) : "
        f"{geom_info['max_radius_from_COM_A']:.3f}"
    )
    print(f"  wall radius (A)    : {wall_radius_A:.3f}")
    print(f"  prepared in        : {replica_dir}")

    return replica_dir


def xtb_command(args, input_name: str):
    cmd = [
        args.xtb,
        "system_centered.pdb",
        "--gfn", str(args.gfn),
        "--chrg", str(args.charge),
        "--uhf", str(args.uhf),
        "--md",
        "--input", input_name,
    ]

    if args.alpb:
        cmd += ["--alpb", args.alpb]

    return cmd


def stage_archive_dir(replica_dir: Path, stage_name: str) -> Path:
    return replica_dir / "stages" / stage_name


def restore_restart_from_previous(
    replica_dir: Path,
    stage_index: int,
):
    if stage_index == 0:
        restart_file = replica_dir / "mdrestart"
        if restart_file.exists():
            restart_file.unlink()
        return

    prev_name = STAGES[stage_index - 1]["name"]

    archived = (
        stage_archive_dir(replica_dir, prev_name)
        / "mdrestart"
    )

    if not archived.exists():
        raise RuntimeError(
            f"Previous restart not found: {archived}. "
            "Run/complete the previous stage first."
        )

    shutil.copy2(
        archived,
        replica_dir / "mdrestart",
    )


def check_md_log(log_path: Path):
    text = log_path.read_text(errors="replace")
    lowered = text.lower()

    for pattern in FATAL_MD_PATTERNS:
        if pattern.lower() in lowered:
            raise RuntimeError(
                f"Fatal MD pattern detected: '{pattern}'. "
                f"See {log_path}"
            )

    warnings = []

    if "thermostating problem" in lowered:
        warnings.append("thermostating problem")

    return warnings


def archive_stage_outputs(
    replica_dir: Path,
    stage_name: str,
    log_path: Path,
):
    archive = stage_archive_dir(replica_dir, stage_name)
    archive.mkdir(parents=True, exist_ok=True)

    core_outputs = [
        "xtb.trj",
        "xtb-trj.pdb",
        "mdrestart",
        "xtbmdok",
    ]

    for filename in core_outputs:
        src = replica_dir / filename
        if src.exists():
            shutil.copy2(src, archive / filename)

    if log_path.exists():
        shutil.copy2(
            log_path,
            archive / log_path.name,
        )

    scoord_dir = archive / "scoord"
    scoord_files = sorted(replica_dir.glob("scoord.*"))

    if scoord_files:
        scoord_dir.mkdir(exist_ok=True)

        for src in scoord_files:
            shutil.move(
                str(src),
                scoord_dir / src.name,
            )

    (archive / "stage.done").write_text("ok\n")


def run_replica(replica_dir: Path, args):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    env.setdefault("MKL_NUM_THREADS", str(args.threads))
    env.setdefault("OMP_STACKSIZE", "4G")

    for i, stage in enumerate(STAGES):
        name = stage["name"]
        archive = stage_archive_dir(replica_dir, name)

        if (
            (archive / "stage.done").exists()
            and not args.force
        ):
            print(
                f"  SKIP {replica_dir.parent.name}/"
                f"{replica_dir.name} {name}: already completed"
            )
            continue

        restore_restart_from_previous(
            replica_dir,
            i,
        )

        for transient in [
            "xtb.trj",
            "xtb-trj.pdb",
            "xtbmdok",
        ]:
            p = replica_dir / transient
            if p.exists():
                p.unlink()

        input_name = f"{name}.inp"
        log_path = replica_dir / f"{name}.out"
        cmd = xtb_command(args, input_name)

        print(
            f"  RUN  {replica_dir.parent.name}/"
            f"{replica_dir.name} {name}"
        )
        print(f"       {' '.join(cmd)}")

        with log_path.open("w") as log:
            result = subprocess.run(
                cmd,
                cwd=replica_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                text=True,
            )

        if result.returncode != 0:
            raise RuntimeError(
                f"xTB failed at {name} for "
                f"{replica_dir.parent.name}/{replica_dir.name}. "
                f"See {log_path}"
            )

        warnings = check_md_log(log_path)

        if not (replica_dir / "mdrestart").exists():
            raise RuntimeError(
                f"{name} ended without mdrrestart for "
                f"{replica_dir.parent.name}/{replica_dir.name}. "
                f"See {log_path}"
            )

        if not (replica_dir / "xtb.trj").exists():
            raise RuntimeError(
                f"{name} ended without xtb.trj for "
                f"{replica_dir.parent.name}/{replica_dir.name}. "
                f"See {log_path}"
            )

        archive_stage_outputs(
            replica_dir,
            name,
            log_path,
        )

        if warnings:
            print(
                f"  WARNING {name}: "
                + ", ".join(warnings)
            )

        print(f"  OK   {name}")

    final_archived = (
        stage_archive_dir(
            replica_dir,
            STAGES[-1]["name"],
        )
        / "mdrestart"
    )

    if final_archived.exists():
        shutil.copy2(
            final_archived,
            replica_dir / "mdrestart_final",
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Prepare spherical Packmol droplets and optionally run "
            "xTB MD E2 thermalization + screening."
        )
    )

    select = p.add_mutually_exclusive_group()

    select.add_argument(
        "--system",
        nargs="+",
        choices=list(SYSTEMS),
        help="Specific system(s) to prepare/run.",
    )

    select.add_argument(
        "--main",
        action="store_true",
        help="Use the four main systems.",
    )

    select.add_argument(
        "--controls",
        action="store_true",
        help="Use control systems.",
    )

    select.add_argument(
        "--all",
        action="store_true",
        help="Use all systems.",
    )

    p.add_argument(
        "--project",
        type=Path,
        default=ROOT / "md_screening",
        help="Output directory (default: ROOT/md_screening).",
    )

    p.add_argument(
        "--replicas",
        type=int,
        default=2,
        help="Independent replicas per system (default: 2).",
    )

    # Packmol
    p.add_argument(
        "--packmol",
        default="packmol",
        help="Packmol executable (default: packmol).",
    )

    p.add_argument(
        "--water-pdb",
        type=Path,
        default=None,
        help="Optional explicit water template PDB.",
    )

    p.add_argument(
        "--sphere-radius",
        type=float,
        default=12.0,
        help="Packmol droplet radius in A (default: 12.0).",
    )

    p.add_argument(
        "--density",
        type=float,
        default=1.0,
        help=(
            "Target initial total mass density in g/cm3 used to "
            "estimate water count (default: 1.0). Ignored with --waters."
        ),
    )

    p.add_argument(
        "--waters",
        type=int,
        default=None,
        help="Override automatic water count with a fixed number.",
    )

    p.add_argument(
        "--packmol-tolerance",
        type=float,
        default=2.0,
        help="Packmol tolerance in A (default: 2.0).",
    )

    p.add_argument(
        "--seed-base",
        type=int,
        default=271828,
        help=(
            "Base Packmol seed. Replica r uses seed_base + r "
            "(default: 271828)."
        ),
    )

    p.add_argument(
        "--repack",
        action="store_true",
        help="Re-run Packmol even if packed_sphere.pdb already exists.",
    )

    # xTB
    p.add_argument(
        "--run",
        action="store_true",
        help=(
            "Run xTB after preparation. Without --run, only Packmol "
            "packing + centering + input generation are performed."
        ),
    )

    p.add_argument(
        "--xtb",
        default="xtb",
        help="xTB executable (default: xtb).",
    )

    p.add_argument(
        "--threads",
        type=int,
        default=8,
        help="OMP threads per xTB run (default: 8).",
    )

    p.add_argument(
        "--gfn",
        type=int,
        default=2,
        choices=[0, 1, 2],
        help="GFN parametrization (default: 2).",
    )

    p.add_argument(
        "--charge",
        type=int,
        default=0,
        help="Total charge (default: 0).",
    )

    p.add_argument(
        "--uhf",
        type=int,
        default=0,
        help="Unpaired-electron/UHF setting (default: 0).",
    )

    p.add_argument(
        "--alpb",
        default=None,
        help=(
            "Optional ALPB solvent, e.g. --alpb water. "
            "Default: off."
        ),
    )

    p.add_argument(
        "--wall-margin",
        type=float,
        default=0.75,
        help=(
            "Wall radius beyond the outermost centered atom in A "
            "(default: 0.75)."
        ),
    )

    p.add_argument(
        "--force",
        action="store_true",
        help="Re-run xTB stages even if stage.done exists.",
    )

    return p.parse_args()


def select_systems(args):
    if args.all:
        selected = list(SYSTEMS)

    elif args.main:
        selected = MAIN_SYSTEMS.copy()

    elif args.controls:
        selected = CONTROL_SYSTEMS.copy()

    elif args.system:
        selected = args.system

    else:
        raise SystemExit(
            "No system selected.\n"
            "Use --system, --main, --controls or --all."
        )

    return list(dict.fromkeys(selected))


def validate_args(args):
    if args.replicas < 1:
        raise SystemExit("--replicas must be >= 1.")

    if args.sphere_radius <= 0:
        raise SystemExit("--sphere-radius must be > 0.")

    if args.density <= 0:
        raise SystemExit("--density must be > 0.")

    if args.waters is not None and args.waters < 1:
        raise SystemExit("--waters must be >= 1.")

    if args.packmol_tolerance <= 0:
        raise SystemExit("--packmol-tolerance must be > 0.")

    if args.wall_margin <= 0:
        raise SystemExit("--wall-margin must be > 0.")


def main():
    args = parse_args()
    validate_args(args)

    selected = select_systems(args)
    water_pdb = find_water_pdb(args.water_pdb)

    for name in selected:
        solute = SYSTEMS[name]["solute"]

        if not solute.exists():
            raise SystemExit(
                f"Unsolvated PDB not found for {name}:\n{solute}"
            )

    if args.run and shutil.which(args.xtb) is None:
        raise SystemExit(
            f"xTB executable '{args.xtb}' not found in PATH. "
            "Use --xtb /path/to/xtb if necessary."
        )

    if shutil.which(args.packmol) is None:
        raise SystemExit(
            f"Packmol executable '{args.packmol}' not found in PATH. "
            "Use --packmol /path/to/packmol if necessary."
        )

    args.project.mkdir(
        parents=True,
        exist_ok=True,
    )

    replica_dirs = []

    print(f"Water template: {water_pdb}")
    print(f"Output project: {args.project}")
    print(f"Replicas/system: {args.replicas}")

    for system_name in selected:
        solute = SYSTEMS[system_name]["solute"]

        for replica_index in range(
            1,
            args.replicas + 1,
        ):
            replica_dir = prepare_replica(
                system_name=system_name,
                solute_pdb=solute,
                water_pdb=water_pdb,
                replica_index=replica_index,
                project_dir=args.project,
                args=args,
            )

            replica_dirs.append(replica_dir)

    if args.run:
        print("\nStarting xTB MD E2 pipeline...")

        for replica_dir in replica_dirs:
            run_replica(
                replica_dir,
                args,
            )

        print("\nAll requested replicas completed.")

    else:
        print(
            "\nPreparation only completed.\n"
            "Inspect each packing/packed_sphere.pdb and "
            "system_centered.pdb before running xTB.\n"
            "Then repeat the same command with --run."
        )


if __name__ == "__main__":
    main()
