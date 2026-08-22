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
5. Creates a solvent-only xTB pre-relaxation input and xTB MD inputs for:
      00_relax       short optimization with Zn(His)2 fixed
      01_100K        1.0 ps initial thermal preparation (not production)
      02_200K        2.0 ps intermediate heating (not production)
      03_298K_equil  3.0 ps final-temperature equilibration
      04_298K_screen 5.0 ps
      05_298K_extended 20.0 ps opt-in continuation at 298.15 K
6. Optionally runs xTB sequentially, chaining stages through mdrestart.
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
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess
import sys
from typing import Iterable

BOHR_PER_ANGSTROM = 1.8897261254578281
AVOGADRO = 6.02214076e23
WATER_MOLAR_MASS = 18.01528  # g/mol
MD_STEP_FS = 0.5
MD_DUMP_FS = 10.0

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
        "time": 1.0,
        "restart": False,
        "purpose": (
            "Initial thermal preparation of the droplet; not 100 K production."
        ),
    },
    {
        "name": "02_200K",
        "temp": 200.0,
        "time": 2.0,
        "restart": True,
        "purpose": "Intermediate heating; not 200 K production.",
    },
    {
        "name": "03_298K_equil",
        "temp": 298.15,
        "time": 3.0,
        "restart": True,
        "purpose": "Equilibration at the final temperature.",
    },
    {
        "name": "04_298K_screen",
        "temp": 298.15,
        "time": 5.0,
        "restart": True,
        "purpose": "Structural screening at 298.15 K.",
    },
    {
        "name": "05_298K_extended",
        "temp": 298.15,
        "time": 20.0,
        "steps": 40000,
        "restart": True,
        "restart_from": "04_298K_screen",
        "preserve_velocities": True,
        "continuation": True,
        "continuation_of": "04_298K_screen",
        "coordinates_reinitialized": False,
        "velocities_reinitialized": False,
        "solvation_rebuilt": False,
        "temperature_ramp_applied": False,
        "cumulative_nominal_time_ps": 25.0,
        "purpose": (
            "Opt-in continuation of the final 04_298K_screen restart at "
            "298.15 K; not run by the default pipeline."
        ),
    },
]

CORE_STAGE_COUNT = 4
EXTENDED_STAGE_INDEX = CORE_STAGE_COUNT
EXTENDED_STAGE = STAGES[EXTENDED_STAGE_INDEX]
EXTENDED_STAGE_NAME = EXTENDED_STAGE["name"]
EXECUTION_STAGES = ["00_relax", *(stage["name"] for stage in STAGES)]
EXECUTION_PROVENANCE_ONLY_FIELDS = frozenset({"threads"})
CO2_PACK_STAGE = "06_CO2_shell_pack"
CO2_ACCOMMODATION_STAGE = "07_CO2_accommodation"
CO2_EQUIL_STAGE = "08_CO2_298K_equil"
CO2_SCREEN_STAGE = "09_CO2_298K_screen"
CO2_EXTENDED_STAGE = "10_CO2_298K_extended"
CO2_MD_STAGE_NAMES = (
    CO2_EQUIL_STAGE,
    CO2_SCREEN_STAGE,
    CO2_EXTENDED_STAGE,
)
CO2_START_STAGE_CHOICES = (
    CO2_ACCOMMODATION_STAGE,
    *CO2_MD_STAGE_NAMES,
)
CO2_WORKFLOW_NOTE = (
    "This CO2 workflow is a new condition derived from a previously "
    "equilibrated aqueous configuration. It is not a dynamical continuation "
    "of stage 05."
)
CO2_SITE_DIRECTED_NOTE = (
    "This condition uses a site-directed initial placement of one CO2 "
    "molecule. The directional restriction is applied only during Packmol "
    "generation and does not persist during xTB accommodation or MD. "
    "Therefore this trajectory tests the evolution of a deliberately "
    "site-proximal initial condition and must not be interpreted as "
    "spontaneous CO2 access from bulk."
)
PDB_COORDINATE_TOLERANCE_A = 0.002
SPEED_OF_LIGHT_CM_S = 2.99792458e10

FATAL_MD_PATTERNS = [
    "MD is unstable",
    "emergency exit",
    "Runtime exception",
    "segmentation fault",
    "floating point exception",
]

LEGACY_MDRESTART_COUNT_ERROR = re.compile(
    r"^invalid mdrestart atom count:\s*"
    r"found\s+(\d+),\s*expected\s+(\d+)\s*$",
    flags=re.IGNORECASE,
)

FATAL_XTB_PATTERNS = [
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


def pdb_atoms(path: Path):
    """Read PDB atoms while retaining template line information."""
    return read_pdb_atoms(path.read_text().splitlines(keepends=True))


def validate_packed_solute(solute_pdb: Path, packed_pdb: Path):
    """Validate Packmol's documented solute-first atom ordering."""
    solute = pdb_atoms(solute_pdb)
    packed = pdb_atoms(packed_pdb)
    n_solute = len(solute)

    if len(packed) <= n_solute:
        raise ValueError(
            f"Packed system must contain solvent: found {len(packed)} total "
            f"atoms and {n_solute} solute atoms in {packed_pdb}."
        )

    expected = [atom["element"] for atom in solute]
    observed = [atom["element"] for atom in packed[:n_solute]]
    if observed != expected:
        mismatch = next(
            i for i, (a, b) in enumerate(zip(expected, observed), start=1)
            if a != b
        )
        raise ValueError(
            "Packed PDB solute-first validation failed at atom "
            f"{mismatch}: expected {expected[mismatch - 1]}, found "
            f"{observed[mismatch - 1]}. Refusing to freeze possibly wrong atoms."
        )

    return n_solute, len(packed)


def replace_pdb_coordinates(
    template_pdb: Path,
    coordinates,
    output_pdb: Path,
    elements=None,
):
    """Write coordinates into a PDB template, preserving atom metadata/order."""
    lines = template_pdb.read_text().splitlines(keepends=True)
    atoms = read_pdb_atoms(lines)
    coordinates = list(coordinates)

    if len(coordinates) != len(atoms):
        raise ValueError(
            f"Geometry has {len(coordinates)} atoms; expected {len(atoms)}."
        )
    if elements is not None:
        normalized = [
            value[0].upper() + value[1:].lower() for value in elements
        ]
        expected = [atom["element"] for atom in atoms]
        if normalized != expected:
            raise ValueError(
                "Optimized geometry element sequence differs from PDB template."
            )

    atom_by_line = {
        atom["line_index"]: xyz for atom, xyz in zip(atoms, coordinates)
    }
    output = []
    for index, line in enumerate(lines):
        if index not in atom_by_line:
            output.append(line)
            continue
        x, y, z = atom_by_line[index]
        core = line.rstrip("\r\n").ljust(54)
        output.append(
            core[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + core[54:] + "\n"
        )
    output_pdb.write_text("".join(output))


def read_xyz_geometry(path: Path):
    lines = path.read_text(errors="replace").splitlines()
    try:
        n_atoms = int(lines[0].strip())
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid XYZ atom count in {path}.") from exc
    if len(lines) < n_atoms + 2:
        raise ValueError(f"Incomplete XYZ geometry in {path}.")

    elements, coordinates = [], []
    for line in lines[2:n_atoms + 2]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"Invalid XYZ atom record in {path}: {line}")
        elements.append(fields[0])
        coordinates.append(tuple(float(value) for value in fields[1:4]))
    return elements, coordinates


def read_coord_geometry(path: Path):
    """Read an xTB/Turbomole coord file (coordinates are in bohr)."""
    lines = path.read_text(errors="replace").splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == "$coord")
    except StopIteration as exc:
        raise ValueError(f"No $coord section in {path}.") from exc

    elements, coordinates = [], []
    for line in lines[start + 1:]:
        if line.lstrip().startswith("$"):
            break
        fields = line.split()
        if len(fields) < 4:
            continue
        coordinates.append(
            tuple(float(value) / BOHR_PER_ANGSTROM for value in fields[:3])
        )
        elements.append(fields[3])
    if not coordinates:
        raise ValueError(f"No atoms in $coord section of {path}.")
    return elements, coordinates


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


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def packing_signature(
    solute_pdb: Path,
    water_pdb: Path,
    radius_A: float,
    nwater: int,
    tolerance_A: float,
    seed: int,
) -> dict:
    return {
        "solute_sha256": file_sha256(solute_pdb),
        "water_sha256": file_sha256(water_pdb),
        "sphere_radius_A": radius_A,
        "water_count": nwater,
        "packmol_tolerance_A": tolerance_A,
        "packmol_seed": seed,
    }


def validate_packing_provenance(
    packed_pdb: Path,
    packing_manifest: Path,
    expected_signature: dict,
):
    if not packing_manifest.exists():
        raise RuntimeError(
            f"Existing {packed_pdb} has no provenance metadata. "
            "Use --repack once to establish a validated packing."
        )
    try:
        recorded = json.loads(packing_manifest.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Cannot read packing provenance from {packing_manifest}. "
            "Use --repack to regenerate it."
        ) from exc
    if recorded != expected_signature:
        raise RuntimeError(
            "Existing packed_sphere.pdb was generated with different packing "
            "parameters. Use --repack to regenerate it."
        )


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
    packing_manifest = packing_dir / "packing_manifest.json"

    if args.waters is None:
        nwater = estimate_water_count(
            solute_pdb=solute_pdb,
            radius_A=args.sphere_radius,
            density_g_cm3=args.density,
        )
    else:
        nwater = args.waters

    signature = packing_signature(
        solute_pdb=solute_pdb,
        water_pdb=water_pdb,
        radius_A=args.sphere_radius,
        nwater=nwater,
        tolerance_A=args.packmol_tolerance,
        seed=seed,
    )

    if packed_pdb.exists() and not args.repack:
        validate_packing_provenance(
            packed_pdb, packing_manifest, signature
        )
        return packed_pdb, nwater, packmol_inp, packmol_log

    shutil.copy2(solute_pdb, solute_local)
    shutil.copy2(water_pdb, water_local)

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

    packmol_exe = shutil.which(args.packmol)
    if packmol_exe is None:
        raise RuntimeError(
            f"Packmol executable '{args.packmol}' not found in PATH. "
            "Use --packmol /path/to/packmol if needed."
        )

    if args.repack and packed_pdb.exists():
        packed_pdb.unlink()

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

    packing_manifest.write_text(json.dumps(signature, indent=2))

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

def md_step_count(stage: dict) -> int:
    """Derive the integral MD step count from the time in ps."""
    exact_steps = stage["time"] * 1000.0 / MD_STEP_FS
    steps = round(exact_steps)
    if not math.isclose(exact_steps, steps, abs_tol=1.0e-9):
        raise ValueError(
            f"Stage {stage['name']} time ({stage['time']} ps) is not an "
            f"integral number of {MD_STEP_FS} fs steps."
        )
    if "steps" in stage and stage["steps"] != steps:
        raise ValueError(
            f"Stage {stage['name']} records {stage['steps']} steps, but "
            f"{stage['time']} ps at {MD_STEP_FS} fs requires {steps}."
        )
    return steps


def md_input(stage, wall_radius_bohr: float) -> str:
    restart = "true" if stage["restart"] else "false"
    dump_fs = stage.get("dump_fs", MD_DUMP_FS)

    return f"""$md
   temp={stage['temp']:.2f}
   time={stage['time']:.3f}
   dump={dump_fs:.1f}
   step={MD_STEP_FS:.1f}
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


def relax_input(
    n_solute_atoms: int,
    wall_radius_bohr: float,
    optimization_engine: str = "auto",
) -> str:
    opt_block = ""
    if optimization_engine != "auto":
        opt_block = f"""$opt
   engine={optimization_engine}
$end

"""

    return f"""{opt_block}$fix
   atoms: 1-{n_solute_atoms}
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
    n_solute_atoms: int,
    args,
):
    manifest_path = replica_dir / "manifest.json"
    previous_result = None
    if manifest_path.exists():
        try:
            previous_result = json.loads(
                manifest_path.read_text()
            ).get("relaxation_result")
        except (OSError, json.JSONDecodeError):
            pass

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
        "relaxation": {
            "enabled": not args.skip_relax,
            "mode": "solvent_only",
            "solute_fixed": not args.skip_relax,
            "fixed_atoms": (
                f"1-{n_solute_atoms}" if not args.skip_relax else None
            ),
            "n_solute_atoms": n_solute_atoms,
            "n_mobile_atoms": geom_info["n_atoms"] - n_solute_atoms,
            "gfn": args.gfn,
            "optimization_level": args.relax_level,
            "optimization_engine": args.relax_engine,
            "max_cycles": args.relax_cycles,
            "wall_type": "logfermi_sphere",
            "input_geometry": "system_centered.pdb",
            "output_geometry": (
                "system_relaxed.pdb" if not args.skip_relax else None
            ),
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
            "step_fs": MD_STEP_FS,
            "dump_fs": MD_DUMP_FS,
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
    if (
        previous_result is not None
        and (stage_archive_dir(replica_dir, "00_relax") / "stage.done").exists()
        and not args.force
    ):
        data["relaxation_result"] = previous_result

    manifest_path.write_text(json.dumps(data, indent=2))


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
    n_solute_atoms, n_total_atoms = validate_packed_solute(
        solute_pdb, centered_pdb
    )
    if n_total_atoms != geom_info["n_atoms"]:
        raise RuntimeError("Internal atom-count mismatch after centering.")

    # The wall follows the ACTUAL centered droplet, not just the requested
    # Packmol radius. This avoids putting initial atoms inside the wall.
    wall_radius_A = (
        geom_info["max_radius_from_COM_A"] + args.wall_margin
    )
    wall_radius_bohr = wall_radius_A * BOHR_PER_ANGSTROM

    (replica_dir / "00_relax.inp").write_text(
        relax_input(
            n_solute_atoms,
            wall_radius_bohr,
            args.relax_engine,
        )
    )

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
        n_solute_atoms=n_solute_atoms,
        args=args,
    )

    print(f"\n[{system_name} / replica_{replica_index:02d}]")
    print(f"  solute             : {solute_pdb}")
    print(f"  droplet radius (A) : {args.sphere_radius:.3f}")
    print(f"  waters             : {nwater}")
    print(f"  Packmol seed       : {packmol_seed}")
    print(f"  atoms total        : {geom_info['n_atoms']}")
    print(f"  solute atoms       : {n_solute_atoms}")
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


def xtb_command(
    args,
    geometry_name: str,
    input_name: str,
    *,
    optimize: bool = False,
):
    cmd = [
        args.xtb,
        geometry_name,
        "--gfn", str(args.gfn),
        "--chrg", str(args.charge),
        "--uhf", str(args.uhf),
    ]
    if optimize:
        cmd += [
            "--opt", args.relax_level,
            "--cycles", str(args.relax_cycles),
        ]
    else:
        cmd.append("--md")
    cmd += ["--input", input_name]

    if args.alpb:
        cmd += ["--alpb", args.alpb]

    return cmd


def update_manifest(replica_dir: Path, key: str, value):
    path = replica_dir / "manifest.json"
    data = json.loads(path.read_text())
    data[key] = value
    path.write_text(json.dumps(data, indent=2))


def valid_pdb(path: Path, expected_atoms: int) -> bool:
    try:
        return path.exists() and len(pdb_atoms(path)) == expected_atoms
    except (OSError, ValueError):
        return False


def materialize_relaxed_pdb(
    replica_dir: Path,
    expected_atoms: int,
    *,
    template_name: str = "system_centered.pdb",
    destination_name: str = "system_relaxed.pdb",
    preserve_template_metadata: bool = False,
) -> Path:
    """Find xTB's optimized geometry and create the operational relaxed PDB."""
    template = replica_dir / template_name
    destination = replica_dir / destination_name

    for name in ["xtbopt.pdb", "xtblast.pdb"]:
        candidate = replica_dir / name
        if valid_pdb(candidate, expected_atoms):
            candidate_atoms = pdb_atoms(candidate)
            template_elements = [a["element"] for a in pdb_atoms(template)]
            if [a["element"] for a in candidate_atoms] != template_elements:
                raise RuntimeError(
                    f"{name} element sequence differs from input."
                )
            if preserve_template_metadata:
                replace_pdb_coordinates(
                    template,
                    [atom["xyz"] for atom in candidate_atoms],
                    destination,
                    elements=[atom["element"] for atom in candidate_atoms],
                )
            else:
                shutil.copy2(candidate, destination)
            return destination

    errors = []
    for name, reader in [
        ("xtbopt.xyz", read_xyz_geometry),
        ("xtbopt.coord", read_coord_geometry),
        ("xtblast.xyz", read_xyz_geometry),
        ("xtblast.coord", read_coord_geometry),
    ]:
        candidate = replica_dir / name
        if not candidate.exists():
            continue
        try:
            elements, coordinates = reader(candidate)
            if len(coordinates) != expected_atoms:
                raise ValueError(
                    f"found {len(coordinates)} atoms, expected {expected_atoms}"
                )
            replace_pdb_coordinates(
                template, coordinates, destination, elements=elements
            )
            return destination
        except (OSError, ValueError) as exc:
            errors.append(f"{name}: {exc}")

    detail = (
        "; ".join(errors)
        if errors
        else "no xtbopt/xtblast PDB, XYZ, or coord geometry found"
    )
    raise RuntimeError(f"No valid final relaxation geometry: {detail}.")


def solute_displacement(
    initial_pdb: Path,
    relaxed_pdb: Path,
    n_solute_atoms: int,
):
    initial = pdb_atoms(initial_pdb)[:n_solute_atoms]
    relaxed = pdb_atoms(relaxed_pdb)[:n_solute_atoms]
    displacements = [
        math.dist(a["xyz"], b["xyz"]) for a, b in zip(initial, relaxed)
    ]
    return {
        "solute_rmsd_A": math.sqrt(
            sum(value * value for value in displacements) / len(displacements)
        ),
        "solute_max_displacement_A": max(displacements),
    }


def _last_float(pattern: str, text: str):
    matches = re.findall(pattern, text, flags=re.IGNORECASE | re.MULTILINE)
    if not matches:
        return None
    value = matches[-1]
    if isinstance(value, tuple):
        value = next((part for part in value if part), "")
    try:
        return float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        return None


def parse_relaxation_output(text: str, not_converged_file: bool = False):
    """Parse only unambiguous final optimization diagnostics."""
    lowered = text.lower()
    for pattern in FATAL_XTB_PATTERNS:
        if pattern.lower() in lowered:
            raise RuntimeError(
                f"Fatal xTB pattern detected: '{pattern}'."
            )

    failed = re.search(
        r"failed\s+to\s+converge\s+geometry\s+optimization"
        r"(?:\s+in\s+(\d+)\s+cycles?)?",
        text,
        flags=re.IGNORECASE,
    )
    converged_match = re.search(
        r"geometry\s+optimization\s+converged"
        r"(?:\s+(?:in|after)\s+(\d+)\s+(?:cycles?|iterations?))?",
        text,
        flags=re.IGNORECASE,
    )

    converged = None
    cycles = None
    if failed or not_converged_file:
        converged = False
        if failed and failed.group(1):
            cycles = int(failed.group(1))
    elif converged_match:
        converged = True
        if converged_match.group(1):
            cycles = int(converged_match.group(1))

    if cycles is None:
        cycle_numbers = [
            int(value) for value in re.findall(
                r"^\s*cycle\s+(\d+)\b",
                text,
                flags=re.IGNORECASE | re.MULTILINE,
            )
        ]
        cycles = max(cycle_numbers) if cycle_numbers else None

    number = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)"
    final_energy = _last_float(
        rf"^(?![^\n]*\bINITIAL\s+(?:TOTAL\s+)?ENERGY\b)"
        rf"[^\n]*\bTOTAL\s+ENERGY\s*:?\s*{number}\s*Eh\b",
        text,
    )
    final_gradient = _last_float(
        rf"\bGRADIENT\s+NORM\s*:?\s*{number}\s*Eh\s*(?:/|per)\s*(?:a0|α)\b",
        text,
    )
    initial_energy = _last_float(
        rf"\bINITIAL\s+(?:TOTAL\s+)?ENERGY\s*:?\s*{number}\s*Eh\b",
        text,
    )

    return {
        "converged": converged,
        "cycles": cycles,
        "initial_energy_Eh": initial_energy,
        "final_energy_Eh": final_energy,
        "delta_energy_Eh": (
            final_energy - initial_energy
            if initial_energy is not None and final_energy is not None
            else None
        ),
        "final_gradient_Eh_a0": final_gradient,
    }


def relaxation_diagnostics(log_path: Path, not_converged_file: bool):
    try:
        return parse_relaxation_output(
            log_path.read_text(errors="replace"),
            not_converged_file=not_converged_file,
        )
    except RuntimeError as exc:
        raise RuntimeError(f"{exc} See {log_path}") from exc


def relaxation_signature(replica_dir: Path, manifest: dict, args) -> dict:
    """Inputs that must match before a completed relaxation can be reused."""
    centered = replica_dir / "system_centered.pdb"
    return {
        "input_sha256": hashlib.sha256(centered.read_bytes()).hexdigest(),
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "optimization_level": args.relax_level,
        "optimization_engine": args.relax_engine,
        "max_cycles": args.relax_cycles,
        "wall_radius_bohr": manifest["wall"]["radius_bohr"],
        "fixed_atoms": manifest["relaxation"]["fixed_atoms"],
    }


def stage_archive_dir(replica_dir: Path, stage_name: str) -> Path:
    return replica_dir / "stages" / stage_name


def validate_historical_relaxation(replica_dir: Path):
    """Validate a completed relaxation against its own historical records."""
    archive = stage_archive_dir(replica_dir, "00_relax")
    done = archive / "stage.done"
    if not done.is_file():
        raise RuntimeError(
            f"Cannot reuse historical 00_relax: missing {done}."
        )
    failed = archive / "stage.failed"
    if failed.exists():
        raise RuntimeError(
            "Cannot reuse historical 00_relax: conflicting stage.done and "
            f"stage.failed markers exist in {archive}."
        )

    manifest_path = replica_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Cannot reuse historical 00_relax: missing or invalid "
            f"{manifest_path}."
        ) from exc
    if not isinstance(manifest, dict):
        raise RuntimeError(
            f"Cannot reuse historical 00_relax: invalid {manifest_path}."
        )

    relaxation = manifest.get("relaxation")
    if not isinstance(relaxation, dict):
        raise RuntimeError(
            "Cannot reuse historical 00_relax: manifest relaxation record "
            "is missing or invalid."
        )
    if relaxation.get("enabled") is not True:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: manifest records "
            "relaxation.enabled != true."
        )

    relaxation_result = manifest.get("relaxation_result")
    configuration = (
        relaxation_result.get("configuration")
        if isinstance(relaxation_result, dict)
        else None
    )
    if not isinstance(configuration, dict) or not configuration:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: "
            "relaxation_result.configuration is missing or invalid."
        )

    centered_pdb = replica_dir / "system_centered.pdb"
    relaxed_pdb = replica_dir / "system_relaxed.pdb"
    for geometry in (centered_pdb, relaxed_pdb):
        if not geometry.is_file():
            raise RuntimeError(
                f"Cannot reuse historical 00_relax: missing {geometry}."
            )

    try:
        geometry_record = manifest["geometry_after_centering"]
        wall_record = manifest["wall"]
        xtb_record = manifest["xtb"]
        n_total = geometry_record["n_atoms"]
        n_solute = relaxation["n_solute_atoms"]
        fixed_atoms = relaxation["fixed_atoms"]
        relaxation_gfn = relaxation["gfn"]
        optimization_level = relaxation["optimization_level"]
        optimization_engine = relaxation.get(
            "optimization_engine", "auto"
        )
        max_cycles = relaxation["max_cycles"]
        wall_radius_bohr = wall_record["radius_bohr"]
        xtb_gfn = xtb_record["gfn"]
        charge = xtb_record["charge"]
        uhf = xtb_record["uhf"]
        alpb = xtb_record["alpb"]
    except (KeyError, TypeError) as exc:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: required historical manifest "
            f"field is missing or invalid ({exc})."
        ) from exc

    if (
        not isinstance(n_total, int)
        or isinstance(n_total, bool)
        or n_total < 1
        or not isinstance(n_solute, int)
        or isinstance(n_solute, bool)
        or not 1 <= n_solute <= n_total
    ):
        raise RuntimeError(
            "Cannot reuse historical 00_relax: invalid historical atom "
            f"counts (n_atoms={n_total!r}, n_solute_atoms={n_solute!r})."
        )

    expected_fixed_atoms = f"1-{n_solute}"
    if fixed_atoms != expected_fixed_atoms:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: relaxation.fixed_atoms "
            f"is {fixed_atoms!r}, expected {expected_fixed_atoms!r} from "
            "relaxation.n_solute_atoms."
        )
    if configuration.get("fixed_atoms") != fixed_atoms:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: configuration.fixed_atoms "
            "does not match relaxation.fixed_atoms."
        )

    try:
        centered_atoms = pdb_atoms(centered_pdb)
        relaxed_atoms = pdb_atoms(relaxed_pdb)
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: invalid historical geometry "
            f"({exc})."
        ) from exc
    if len(centered_atoms) != n_total:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: system_centered.pdb has "
            f"{len(centered_atoms)} atoms, expected {n_total}."
        )
    if len(relaxed_atoms) != n_total:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: system_relaxed.pdb has "
            f"{len(relaxed_atoms)} atoms, expected {n_total}."
        )
    for geometry_name, atoms in (
        ("system_centered.pdb", centered_atoms),
        ("system_relaxed.pdb", relaxed_atoms),
    ):
        if any(
            not all(math.isfinite(value) for value in atom["xyz"])
            for atom in atoms
        ):
            raise RuntimeError(
                "Cannot reuse historical 00_relax: "
                f"{geometry_name} contains non-finite coordinates."
            )
    if [
        atom["element"] for atom in centered_atoms
    ] != [
        atom["element"] for atom in relaxed_atoms
    ]:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: element sequence differs "
            "between system_centered.pdb and system_relaxed.pdb."
        )

    centered_sha256 = file_sha256(centered_pdb)
    if configuration.get("input_sha256") != centered_sha256:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: system_centered.pdb SHA-256 "
            "does not match configuration.input_sha256."
        )

    if relaxation_gfn != xtb_gfn:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: relaxation.gfn does not "
            "match xtb.gfn."
        )

    expected_configuration = {
        "gfn": xtb_gfn,
        "charge": charge,
        "uhf": uhf,
        "alpb": alpb,
        "optimization_level": optimization_level,
        "optimization_engine": optimization_engine,
        "max_cycles": max_cycles,
        "wall_radius_bohr": wall_radius_bohr,
    }
    historical_configuration = dict(configuration)
    # Results predating --relax-engine used xTB's default engine selection.
    historical_configuration.setdefault("optimization_engine", "auto")
    mismatches = [
        field
        for field, expected in expected_configuration.items()
        if historical_configuration.get(field) != expected
    ]
    if mismatches:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: historical configuration "
            "does not match the manifest (fields: "
            + ", ".join(mismatches)
            + ")."
        )

    displacement = solute_displacement(
        centered_pdb,
        relaxed_pdb,
        n_solute,
    )
    if displacement["solute_max_displacement_A"] > 1.0e-3:
        raise RuntimeError(
            "Cannot reuse historical 00_relax: fixed solute atoms moved by "
            f"{displacement['solute_max_displacement_A']:.6f} A, exceeding "
            "the 0.001 A tolerance."
        )

    print("  REUSE 00_relax")


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


def parse_md_thermal_output(text: str) -> dict:
    lowered = text.lower()
    fatal_patterns = []
    for pattern in FATAL_MD_PATTERNS:
        if pattern.lower() in lowered:
            fatal_patterns.append(pattern)

    average_temperature = None
    average_sections = list(re.finditer(
        r"average\s+properties",
        text,
        flags=re.IGNORECASE,
    ))
    if average_sections:
        section = text[average_sections[-1].end():]
        number = (
            r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)"
            r"(?:[EeDd][-+]?\d+)?)"
        )
        match = re.search(
            rf"^\s*(?:\|\s*)?T\s*:\s*{number}\b",
            section,
            flags=re.IGNORECASE | re.MULTILINE,
        )
        if match:
            try:
                average_temperature = float(
                    match.group(1).replace("D", "E").replace("d", "e")
                )
            except ValueError:
                pass

    return {
        "fatal_patterns": fatal_patterns,
        "thermostating_problem": "thermostating problem" in lowered,
        "normal_exit_of_md": "normal exit of md()" in lowered,
        "average_temperature_K": average_temperature,
    }


def inspect_md_log(log_path: Path):
    return parse_md_thermal_output(
        log_path.read_text(errors="replace")
    )


def thermostat_warning_allowed(policy: str, stage_name: str) -> bool:
    return policy == "allow"


def md_thermal_result(
    stage: dict,
    validation: dict,
    policy: str,
    warning_accepted: bool,
) -> dict:
    return {
        "target_temperature_K": stage["temp"],
        "average_temperature_K": validation["average_temperature_K"],
        "thermostating_problem": validation["thermostating_problem"],
        "thermostat_warning_policy": policy,
        "warning_accepted": warning_accepted,
        "normal_exit_of_md": validation["normal_exit_of_md"],
    }


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
        "mdrestart.input",
        "xtbmdok",
        f"{stage_name}.inp",
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


def mark_stage_done(archive: Path):
    running = archive / "stage.running"
    if running.exists():
        running.unlink()
    failed = archive / "stage.failed"
    if failed.exists():
        failed.unlink()
    (archive / "stage.done").write_text("ok\n")


def mark_stage_failed(archive: Path, reason: str):
    running = archive / "stage.running"
    if running.exists():
        running.unlink()
    done = archive / "stage.done"
    if done.exists():
        done.unlink()
    archive.mkdir(parents=True, exist_ok=True)
    (archive / "stage.failed").write_text(reason.rstrip() + "\n")


def mark_stage_running(archive: Path, metadata: dict):
    archive.mkdir(parents=True, exist_ok=True)
    (archive / "stage.running").write_text(json.dumps(metadata, indent=2))


def md_stage_configuration(
    replica_dir: Path,
    stage: dict,
    geometry_name: str,
    input_name: str,
    input_restart_sha256: str | None,
    args,
) -> dict:
    manifest = json.loads((replica_dir / "manifest.json").read_text())
    configuration = {
        "stage": stage["name"],
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "threads": args.threads,
        "temp_K": stage["temp"],
        "time_ps": stage["time"],
        "step_fs": MD_STEP_FS,
        "dump_fs": MD_DUMP_FS,
        "hmass": 1,
        "shake": 0,
        "sccacc": 1.0,
        "restart": stage["restart"],
        "wall_radius_bohr": manifest["wall"]["radius_bohr"],
        "input_geometry": geometry_name,
        "input_geometry_sha256": file_sha256(
            replica_dir / geometry_name
        ),
        "xtb_input_sha256": file_sha256(replica_dir / input_name),
        "input_restart_sha256": input_restart_sha256,
    }
    if stage["name"] == EXTENDED_STAGE_NAME:
        configuration.update({
            "steps": md_step_count(stage),
            "restart_from": stage["restart_from"],
            "preserve_velocities": stage["preserve_velocities"],
            "continuation": stage["continuation"],
            "continuation_of": stage["continuation_of"],
            "coordinates_reinitialized": (
                stage["coordinates_reinitialized"]
            ),
            "velocities_reinitialized": stage["velocities_reinitialized"],
            "solvation_rebuilt": stage["solvation_rebuilt"],
            "temperature_ramp_applied": (
                stage["temperature_ramp_applied"]
            ),
            "cumulative_nominal_time_ps": (
                stage["cumulative_nominal_time_ps"]
            ),
            "nvt": True,
            "velocity_output": True,
            "randomize_velocities": False,
            "input_restart_path": str(
                Path("stages") / stage["restart_from"] / "mdrestart"
            ),
        })
    return configuration


def ensure_extended_stage_registered(replica_dir: Path):
    """Append the opt-in stage definition without rewriting prior stages."""
    manifest_path = replica_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
        md = manifest["md"]
        recorded_stages = md["stages"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise RuntimeError(
            f"Cannot register {EXTENDED_STAGE_NAME} in {manifest_path}."
        ) from exc

    if not isinstance(recorded_stages, list):
        raise RuntimeError(
            f"Cannot register {EXTENDED_STAGE_NAME}: md.stages is not a list "
            f"in {manifest_path}."
        )
    if md.get("step_fs") != MD_STEP_FS or md.get("dump_fs") != MD_DUMP_FS:
        raise RuntimeError(
            f"Cannot register {EXTENDED_STAGE_NAME}: manifest timestep/dump "
            "does not match the continuation protocol."
        )

    matches = [
        item for item in recorded_stages
        if isinstance(item, dict) and item.get("name") == EXTENDED_STAGE_NAME
    ]
    if len(matches) > 1:
        raise RuntimeError(
            f"Cannot register {EXTENDED_STAGE_NAME}: duplicate definitions "
            f"exist in {manifest_path}."
        )
    if matches:
        mismatched = configuration_mismatches(EXTENDED_STAGE, matches[0])
        if mismatched:
            raise RuntimeError(
                f"Existing {EXTENDED_STAGE_NAME} definition is incompatible "
                f"(fields: {', '.join(mismatched)})."
            )
        return False

    recorded_stages.append(dict(EXTENDED_STAGE))
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"  REGISTER {EXTENDED_STAGE_NAME} in manifest.json")
    return True


def configuration_mismatches(expected: dict, recorded: dict) -> list[str]:
    """Return compatibility differences, excluding execution-only metadata."""
    return [
        key for key, expected_value in expected.items()
        if key not in EXECUTION_PROVENANCE_ONLY_FIELDS
        if recorded.get(key) != expected_value
    ]


def historical_md_stage_configuration(
    replica_dir: Path,
    stage_index: int,
    input_restart_sha256: str | None,
) -> tuple[dict, dict]:
    """Reconstruct a stage signature from the untouched historical manifest."""
    manifest_path = replica_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
        recorded_stages = manifest["md"]["stages"]
        stage_name = STAGES[stage_index]["name"]
        stage = next(
            item for item in recorded_stages
            if item.get("name") == stage_name
        )
        xtb = manifest["xtb"]
        md = manifest["md"]
        wall_radius_bohr = manifest["wall"]["radius_bohr"]
        relaxation_enabled = manifest["relaxation"]["enabled"]
    except (OSError, json.JSONDecodeError, KeyError, StopIteration) as exc:
        raise RuntimeError(
            f"Cannot reconstruct historical provenance for "
            f"{STAGES[stage_index]['name']} from {manifest_path}."
        ) from exc

    try:
        geometry_name = (
            "system_relaxed.pdb"
            if stage_index == 0 and relaxation_enabled
            else "system_centered.pdb"
        )
        input_name = f"{stage['name']}.inp"
        input_path = replica_dir / input_name
        input_text = input_path.read_text()
        expected_input_text = md_input(stage, wall_radius_bohr)
        if input_text != expected_input_text:
            raise RuntimeError(
                f"Historical input {input_path} does not match the stage "
                "configuration recorded in manifest.json."
            )
        configuration = {
            "stage": stage["name"],
            "gfn": xtb["gfn"],
            "charge": xtb["charge"],
            "uhf": xtb["uhf"],
            "alpb": xtb["alpb"],
            "threads": xtb["threads"],
            "temp_K": stage["temp"],
            "time_ps": stage["time"],
            "step_fs": md["step_fs"],
            "dump_fs": md["dump_fs"],
            "hmass": md["hmass"],
            "shake": md["shake"],
            "sccacc": md["sccacc"],
            "restart": stage["restart"],
            "wall_radius_bohr": wall_radius_bohr,
            "input_geometry": geometry_name,
            "input_geometry_sha256": file_sha256(
                replica_dir / geometry_name
            ),
            "xtb_input_sha256": file_sha256(input_path),
            "input_restart_sha256": input_restart_sha256,
        }
        if stage["name"] == EXTENDED_STAGE_NAME:
            configuration.update({
                "steps": md_step_count(stage),
                "restart_from": stage["restart_from"],
                "preserve_velocities": stage["preserve_velocities"],
                "continuation": stage["continuation"],
                "continuation_of": stage["continuation_of"],
                "coordinates_reinitialized": (
                    stage["coordinates_reinitialized"]
                ),
                "velocities_reinitialized": (
                    stage["velocities_reinitialized"]
                ),
                "solvation_rebuilt": stage["solvation_rebuilt"],
                "temperature_ramp_applied": (
                    stage["temperature_ramp_applied"]
                ),
                "cumulative_nominal_time_ps": (
                    stage["cumulative_nominal_time_ps"]
                ),
                "nvt": True,
                "velocity_output": True,
                "randomize_velocities": False,
                "input_restart_path": str(
                    Path("stages") / stage["restart_from"] / "mdrestart"
                ),
            })
    except (KeyError, OSError) as exc:
        raise RuntimeError(
            f"Historical provenance inputs for {stage['name']} are "
            f"incomplete or unreadable in {replica_dir}."
        ) from exc
    return stage, configuration


def validate_current_md_input(
    replica_dir: Path,
    stage: dict,
    *,
    create_missing: bool = False,
):
    """Ensure an existing input is exactly the requested current protocol."""
    manifest_path = replica_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
        wall_radius_bohr = manifest["wall"]["radius_bohr"]
    except (OSError, json.JSONDecodeError, KeyError) as exc:
        raise RuntimeError(
            f"Cannot validate MD input from {manifest_path}."
        ) from exc

    input_path = replica_dir / f"{stage['name']}.inp"
    expected = md_input(stage, wall_radius_bohr)
    if create_missing and not input_path.exists():
        input_path.write_text(expected)
        print(f"  PREPARE {stage['name']} input")
    try:
        actual = input_path.read_text()
    except OSError as exc:
        raise RuntimeError(f"Missing MD input: {input_path}.") from exc
    if actual != expected:
        raise RuntimeError(
            f"{input_path} is incompatible with the current "
            f"{stage['name']} protocol. Resume/start-stage will not overwrite "
            "an existing input."
        )


def write_md_stage_manifest(
    archive: Path,
    configuration: dict,
    output_restart_sha256: str,
    thermal_result: dict,
    returncode: int = 0,
    recovery: dict | None = None,
    runtime_metadata: dict | None = None,
):
    archive.mkdir(parents=True, exist_ok=True)
    data = dict(configuration)
    data["returncode"] = returncode
    data["output_restart_sha256"] = output_restart_sha256
    data["thermal_result"] = thermal_result
    if runtime_metadata is not None:
        data.update(runtime_metadata)
    if recovery is not None:
        data["recovery"] = recovery
    (archive / "stage_manifest.json").write_text(
        json.dumps(data, indent=2)
    )


def validate_completed_md_stage(
    archive: Path,
    expected_configuration: dict,
    warning_policy: str,
) -> dict:
    stage_name = expected_configuration["stage"]
    required = [
        "stage.done",
        "mdrestart",
        "xtb.trj",
        "xtbmdok",
        f"{stage_name}.out",
        "stage_manifest.json",
    ]
    if expected_configuration["restart"]:
        required.append("mdrestart.input")
    conflicting_markers = [
        marker for marker in ("stage.failed", "stage.running")
        if (archive / marker).exists()
    ]
    if conflicting_markers:
        raise RuntimeError(
            f"Stage {stage_name} has stage.done plus conflicting marker(s) "
            f"{', '.join(conflicting_markers)}. "
            "Inspect the archive and rerun with --force."
        )
    missing = [name for name in required if not (archive / name).exists()]
    if missing:
        raise RuntimeError(
            f"Stage {stage_name} is marked done but archived outputs are "
            f"incomplete (missing: {', '.join(missing)}). Inspect the "
            "archive and rerun with --force."
        )

    path = archive / "stage_manifest.json"
    try:
        recorded = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Stage {stage_name} has an unreadable stage_manifest.json. "
            "Inspect the archive and rerun with --force."
        ) from exc

    mismatched = configuration_mismatches(
        expected_configuration,
        recorded,
    )
    if mismatched:
        raise RuntimeError(
            f"Completed stage {stage_name} is incompatible with the current "
            f"MD settings (fields: {', '.join(mismatched)}). "
            "Use --force to rerun it."
        )

    # Legacy successful manifests predate the returncode field. Their
    # stage.done marker was written only after the xTB return-code check;
    # retain that control-flow evidence so existing 01--03 archives remain
    # resumable. Failed-stage promotion is stricter and requires an explicit
    # recorded zero whenever a stage manifest already exists.
    if recorded.get("returncode", 0) != 0:
        raise RuntimeError(
            f"Stage {stage_name} records a nonzero xTB return code."
        )

    if stage_name == EXTENDED_STAGE_NAME:
        required_metadata = [
            "started_at", "completed_at", "command", "xtb_version",
            "hostname", "status", "output_paths", "trajectory_integrity",
            "preflight",
        ]
        missing_metadata = [
            key for key in required_metadata if key not in recorded
        ]
        if missing_metadata:
            raise RuntimeError(
                f"Completed stage {stage_name} lacks runtime provenance "
                f"fields: {', '.join(missing_metadata)}."
            )
        if recorded.get("status") != "completed":
            raise RuntimeError(
                f"Completed stage {stage_name} records status "
                f"{recorded.get('status')!r}, expected 'completed'."
            )
        if not isinstance(recorded.get("command"), list):
            raise RuntimeError(
                f"Completed stage {stage_name} has no valid command record."
            )
        if not isinstance(recorded.get("xtb_version"), str):
            raise RuntimeError(
                f"Completed stage {stage_name} has no xTB version record."
            )
        trajectory_integrity = validate_xtb_trajectory(
            archive / "xtb.trj",
            expected_atoms=restart_atom_count(archive.parent.parent),
            expected_frames=md_expected_frame_count(EXTENDED_STAGE),
            require_velocities=True,
        )
        recorded_integrity = recorded["trajectory_integrity"]
        integrity_fields = [
            "sha256", "size_bytes", "n_atoms", "frames",
            "expected_nominal_frames", "constant_atom_count",
            "constant_element_order", "velocities_present",
            "finite_coordinates_and_velocities",
        ]
        integrity_mismatches = [
            key for key in integrity_fields
            if recorded_integrity.get(key) != trajectory_integrity.get(key)
        ]
        if integrity_mismatches:
            raise RuntimeError(
                f"Completed stage {stage_name} trajectory provenance is "
                "inconsistent (fields: "
                f"{', '.join(integrity_mismatches)})."
            )

    archived_restart_hash = validate_output_restart(
        archive / "mdrestart",
        expected_configuration["input_restart_sha256"],
        expected_atoms=restart_atom_count(archive.parent.parent),
    )
    if recorded.get("output_restart_sha256") != archived_restart_hash:
        raise RuntimeError(
            f"Stage {stage_name} has an archived mdrestart inconsistent with "
            "stage_manifest.json. Inspect the archive and rerun with --force."
        )
    if expected_configuration["restart"]:
        archived_input_hash = file_sha256(archive / "mdrestart.input")
        if (
            archived_input_hash
            != expected_configuration["input_restart_sha256"]
        ):
            raise RuntimeError(
                f"Stage {stage_name} archived input restart does not match "
                "the validated output restart of the previous stage."
            )

    archived_input = archive / f"{stage_name}.inp"
    if (
        archived_input.exists()
        and file_sha256(archived_input)
        != expected_configuration["xtb_input_sha256"]
    ):
        raise RuntimeError(
            f"Stage {stage_name} archived xTB input is inconsistent with "
            "stage_manifest.json."
        )

    thermal_result = recorded.get("thermal_result")
    if not isinstance(thermal_result, dict):
        raise RuntimeError(
            f"Stage {stage_name} has no valid thermal_result in "
            "stage_manifest.json. Inspect the archive and rerun with --force."
        )

    validation = inspect_md_log(archive / f"{stage_name}.out")
    if validation["fatal_patterns"]:
        raise RuntimeError(
            f"Stage {stage_name} archived log contains fatal pattern(s): "
            f"{', '.join(validation['fatal_patterns'])}."
        )
    for key in ["thermostating_problem", "normal_exit_of_md"]:
        if thermal_result.get(key) != validation[key]:
            raise RuntimeError(
                f"Stage {stage_name} thermal_result field '{key}' is "
                "inconsistent with the archived log."
            )
    if bool(thermal_result.get("warning_accepted")) != bool(
        validation["thermostating_problem"]
    ):
        raise RuntimeError(
            f"Stage {stage_name} thermal_result has an inconsistent "
            "warning_accepted value."
        )

    if validation["thermostating_problem"]:
        if not validation["normal_exit_of_md"]:
            raise RuntimeError(
                f"Stage {stage_name} has 'thermostating problem' without "
                "'normal exit of md()'."
            )
        if not thermal_result.get("warning_accepted"):
            raise RuntimeError(
                f"Stage {stage_name} is marked done but its thermostat "
                "warning was not recorded as accepted."
            )
        if not thermostat_warning_allowed(warning_policy, stage_name):
            raise RuntimeError(
                f"Stage {stage_name} thermostat warning is not accepted by "
                f"policy '{warning_policy}'."
            )

    return recorded


def restart_atom_count(replica_dir: Path) -> int:
    try:
        manifest = json.loads((replica_dir / "manifest.json").read_text())
        count = manifest["geometry_after_centering"]["n_atoms"]
    except (OSError, json.JSONDecodeError, KeyError) as exc:
        raise RuntimeError(
            f"Cannot determine restart atom count from "
            f"{replica_dir}/manifest.json."
        ) from exc
    if not isinstance(count, int) or count < 1:
        raise RuntimeError(
            f"Invalid restart atom count in {replica_dir}/manifest.json."
        )
    return count


def parse_mdrestart(path: Path, expected_atoms: int) -> dict:
    """Validate xTB mdrestart records without changing the original file."""
    if not path.is_file():
        raise RuntimeError("missing mdrestart")
    if path.stat().st_size == 0:
        raise RuntimeError("empty mdrestart")

    records = [
        (line_number, line.split())
        for line_number, line in enumerate(
            path.read_text(errors="replace").splitlines(),
            start=1,
        )
        if line.strip()
    ]
    if not records:
        raise RuntimeError("empty mdrestart")

    control_records = 0
    first_line_number, first_fields = records[0]
    if len(first_fields) == 1:
        try:
            control_value = float(
                first_fields[0].replace("D", "E").replace("d", "e")
            )
        except ValueError as exc:
            raise RuntimeError(
                f"invalid mdrestart control record at line "
                f"{first_line_number}: non-numeric value"
            ) from exc
        if not math.isfinite(control_value):
            raise RuntimeError(
                f"invalid mdrestart control record at line "
                f"{first_line_number}: non-finite value"
            )
        control_records = 1

    atom_rows = records[control_records:]
    atom_records = len(atom_rows)
    if atom_records != expected_atoms:
        raise RuntimeError(
            f"invalid mdrestart atom count: found {atom_records} atom "
            f"records, expected {expected_atoms}"
        )

    for atom_index, (line_number, fields) in enumerate(
        atom_rows,
        start=1,
    ):
        if len(fields) != 6:
            raise RuntimeError(
                f"invalid mdrestart atom record {atom_index} at line "
                f"{line_number}: expected 6 numeric columns, "
                f"found {len(fields)}"
            )
        try:
            values = [
                float(value.replace("D", "E").replace("d", "e"))
                for value in fields
            ]
        except ValueError as exc:
            raise RuntimeError(
                f"invalid mdrestart atom record {atom_index} at line "
                f"{line_number}: non-numeric value"
            ) from exc
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(
                f"invalid mdrestart atom record {atom_index} at line "
                f"{line_number}: non-finite value"
            )

    return {
        "control_records": control_records,
        "atom_records": atom_records,
        "expected_atoms": expected_atoms,
        "coordinates_present": True,
        "velocities_present": True,
        "valid": True,
    }


def expected_md_frames(time_ps: float, dump_fs: float, stage_name: str) -> int:
    exact_frames = time_ps * 1000.0 / dump_fs
    frames = round(exact_frames)
    if not math.isclose(exact_frames, frames, abs_tol=1.0e-9):
        raise ValueError(
            f"Stage {stage_name} duration/dump interval does not yield "
            "an integral nominal frame count."
        )
    return frames


def md_expected_frame_count(stage: dict) -> int:
    return expected_md_frames(
        stage["time"],
        stage.get("dump_fs", MD_DUMP_FS),
        stage["name"],
    )


def _finite_values(fields, *, context: str):
    try:
        values = [
            float(value.replace("D", "E").replace("d", "e"))
            for value in fields
        ]
    except ValueError as exc:
        raise RuntimeError(f"{context}: non-numeric value") from exc
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"{context}: non-finite value")
    return values


def validate_xtb_trajectory(
    path: Path,
    *,
    expected_atoms: int,
    expected_frames: int,
    require_velocities: bool,
    expected_elements=None,
) -> dict:
    """Stream-validate xTB extended XYZ coordinates and velocity blocks."""
    if not path.is_file():
        raise RuntimeError(f"missing trajectory: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"empty trajectory: {path}")

    frame_count = 0
    element_order = None
    with path.open("r", errors="replace") as handle:
        while True:
            atom_count_line = handle.readline()
            if atom_count_line == "":
                break
            if not atom_count_line.strip():
                raise RuntimeError(
                    f"Invalid trajectory {path}: blank record before frame "
                    f"{frame_count + 1}."
                )
            try:
                atom_count = int(atom_count_line.strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"Invalid trajectory {path}: frame {frame_count + 1} "
                    "does not start with an atom count."
                ) from exc
            if atom_count != expected_atoms:
                raise RuntimeError(
                    f"Invalid trajectory {path}: frame {frame_count + 1} has "
                    f"{atom_count} atoms, expected {expected_atoms}."
                )

            if handle.readline() == "":
                raise RuntimeError(
                    f"Invalid trajectory {path}: missing comment for frame "
                    f"{frame_count + 1}."
                )

            frame_elements = []
            for atom_index in range(1, expected_atoms + 1):
                line = handle.readline()
                if line == "":
                    raise RuntimeError(
                        f"Invalid trajectory {path}: truncated coordinates "
                        f"in frame {frame_count + 1}."
                    )
                fields = line.split()
                if len(fields) < 4:
                    raise RuntimeError(
                        f"Invalid trajectory {path}: coordinate record "
                        f"{atom_index} in frame {frame_count + 1} has fewer "
                        "than four columns."
                    )
                _finite_values(
                    fields[1:4],
                    context=(
                        f"Invalid trajectory {path}, frame "
                        f"{frame_count + 1}, coordinate {atom_index}"
                    ),
                )
                frame_elements.append(fields[0])

            if element_order is None:
                element_order = tuple(frame_elements)
                if (
                    expected_elements is not None
                    and element_order != tuple(expected_elements)
                ):
                    raise RuntimeError(
                        f"Invalid trajectory {path}: first-frame element "
                        "order differs from the input geometry."
                    )
            elif tuple(frame_elements) != element_order:
                raise RuntimeError(
                    f"Invalid trajectory {path}: element order changes at "
                    f"frame {frame_count + 1}."
                )

            if require_velocities:
                for atom_index in range(1, expected_atoms + 1):
                    line = handle.readline()
                    if line == "":
                        raise RuntimeError(
                            f"Invalid trajectory {path}: missing velocity "
                            f"block in frame {frame_count + 1}."
                        )
                    fields = line.split()
                    if len(fields) != 3:
                        raise RuntimeError(
                            f"Invalid trajectory {path}: velocity record "
                            f"{atom_index} in frame {frame_count + 1} must "
                            "contain exactly three columns."
                        )
                    _finite_values(
                        fields,
                        context=(
                            f"Invalid trajectory {path}, frame "
                            f"{frame_count + 1}, velocity {atom_index}"
                        ),
                    )
            frame_count += 1

    if frame_count == 0:
        raise RuntimeError(f"Trajectory {path} contains no frames.")
    if abs(frame_count - expected_frames) > 1:
        raise RuntimeError(
            f"Trajectory {path} contains {frame_count} frames; expected "
            f"approximately {expected_frames} (allowing +/- 1 for the xTB "
            "initial/final-frame convention)."
        )

    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "size_bytes": path.stat().st_size,
        "n_atoms": expected_atoms,
        "frames": frame_count,
        "expected_nominal_frames": expected_frames,
        "frame_tolerance": 1,
        "constant_atom_count": True,
        "constant_element_order": True,
        "velocities_present": require_velocities,
        "finite_coordinates_and_velocities": True,
    }


def extract_xtb_version(log_path: Path, trajectory_path: Path) -> str | None:
    for path, patterns in (
        (
            log_path,
            [r"\bxtb\s+version\s+([^\n]+)"],
        ),
        (
            trajectory_path,
            [r"\bxtb\s*:\s*([^\n]+)"],
        ),
    ):
        try:
            with path.open("r", errors="replace") as handle:
                text = handle.read(128 * 1024)
        except OSError:
            continue
        for pattern in patterns:
            match = re.search(pattern, text, flags=re.IGNORECASE)
            if match:
                return match.group(1).strip()
    return None


def validate_output_restart(
    output_restart: Path,
    input_restart_sha256: str | None,
    *,
    expected_atoms: int,
) -> str:
    parse_mdrestart(output_restart, expected_atoms)
    output_sha256 = file_sha256(output_restart)
    if (
        input_restart_sha256 is not None
        and output_sha256 == input_restart_sha256
    ):
        raise RuntimeError(
            "output mdrestart is byte-identical to input mdrestart"
        )
    return output_sha256


def recover_legacy_restart_validation_failure(
    stage: dict,
    archive: Path,
    configuration: dict,
    args,
    recovery_requested: bool,
) -> bool:
    """Promote only the known control-record atom-count false positive."""
    failed = archive / "stage.failed"
    if not recovery_requested or not failed.exists():
        return False

    stage_name = stage["name"]
    if (archive / "stage.done").exists():
        raise RuntimeError(
            f"Cannot recover {stage_name}: archive contains both stage.done "
            "and stage.failed."
        )

    reason = failed.read_text(errors="replace").strip()
    reason_match = LEGACY_MDRESTART_COUNT_ERROR.fullmatch(reason)
    if reason_match is None:
        return False

    old_found = int(reason_match.group(1))
    old_expected = int(reason_match.group(2))
    expected_atoms = restart_atom_count(archive.parent.parent)
    if old_expected != expected_atoms or old_found != expected_atoms + 1:
        raise RuntimeError(
            f"Cannot recover {stage_name}: legacy atom-count reason is "
            "inconsistent with the current replica manifest."
        )

    required = ["xtb.trj", "mdrestart", "xtbmdok", f"{stage_name}.out"]
    if stage["restart"]:
        required.append("mdrestart.input")
    missing = [name for name in required if not (archive / name).exists()]
    if missing:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived outputs are incomplete "
            f"(missing: {', '.join(missing)})."
        )

    log_path = archive / f"{stage_name}.out"
    validation = inspect_md_log(log_path)
    if validation["fatal_patterns"]:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived log contains fatal "
            f"pattern(s): {', '.join(validation['fatal_patterns'])}."
        )
    warning_accepted = (
        validation["thermostating_problem"]
        and validation["normal_exit_of_md"]
        and thermostat_warning_allowed(
            args.thermostat_warning_policy,
            stage_name,
        )
    )
    if validation["thermostating_problem"] and not warning_accepted:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived log has an unaccepted "
            "thermostating problem."
        )

    manifest_path = archive / "stage_manifest.json"
    had_stage_manifest = manifest_path.exists()
    if had_stage_manifest:
        try:
            recorded = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(
                f"Cannot recover {stage_name}: unreadable "
                "stage_manifest.json."
            ) from exc
        mismatched = configuration_mismatches(configuration, recorded)
        if mismatched:
            raise RuntimeError(
                f"Cannot recover {stage_name}: archived configuration is "
                f"incompatible (fields: {', '.join(mismatched)})."
            )
        if recorded.get("returncode") != 0:
            raise RuntimeError(
                f"Cannot recover {stage_name}: archived return code is "
                f"{recorded.get('returncode')}."
            )

    input_restart_sha256 = configuration["input_restart_sha256"]
    if stage["restart"]:
        archived_input_hash = file_sha256(archive / "mdrestart.input")
        if archived_input_hash != input_restart_sha256:
            raise RuntimeError(
                f"Cannot recover {stage_name}: archived input restart does "
                "not match the currently valid previous stage."
            )

    restart_metadata = parse_mdrestart(
        archive / "mdrestart",
        expected_atoms,
    )
    if (
        restart_metadata["control_records"] != 1
        or old_found
        != (
            restart_metadata["control_records"]
            + restart_metadata["atom_records"]
        )
    ):
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived mdrestart does not match "
            "the legacy control-record false-positive signature."
        )
    output_restart_sha256 = validate_output_restart(
        archive / "mdrestart",
        input_restart_sha256,
        expected_atoms=expected_atoms,
    )
    if (
        had_stage_manifest
        and recorded.get("output_restart_sha256")
        != output_restart_sha256
    ):
        raise RuntimeError(
            f"Cannot recover {stage_name}: output restart hash is "
            "inconsistent with stage_manifest.json."
        )

    archived_input = archive / f"{stage_name}.inp"
    if (
        archived_input.exists()
        and file_sha256(archived_input)
        != configuration["xtb_input_sha256"]
    ):
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived xTB input hash is "
            "inconsistent with historical provenance."
        )

    thermal_result = md_thermal_result(
        stage,
        validation,
        args.thermostat_warning_policy,
        warning_accepted,
    )
    write_md_stage_manifest(
        archive,
        configuration,
        output_restart_sha256,
        thermal_result,
        returncode=0,
        recovery={
            "promoted_from_stage_failed": True,
            "original_failure_reason": reason,
            "reason": (
                "legacy parser counted mdrestart control record as an atom"
            ),
            "mdrestart_control_records": (
                restart_metadata["control_records"]
            ),
            "mdrestart_atom_records": restart_metadata["atom_records"],
            "xtb_recalculated": False,
            "configuration_source": (
                "existing stage_manifest.json"
                if had_stage_manifest
                else (
                    "legacy manifest.json plus deterministic stage input "
                    "and archived restart/output chain"
                )
            ),
            "returncode_evidence": (
                "recorded returncode == 0"
                if had_stage_manifest
                else (
                    "legacy restart-validation stage.failed marker was "
                    "created only after returncode == 0"
                )
            ),
        },
    )
    mark_stage_done(archive)
    print(
        f"  RECOVER {stage_name}: legacy mdrestart validation false "
        f"positive; {restart_metadata['control_records']} control record + "
        f"{restart_metadata['atom_records']} atom records; archived MD "
        "outputs passed integrity checks; no xTB recalculation performed."
    )
    return True


def recover_thermostat_warning_stage(
    stage: dict,
    archive: Path,
    configuration: dict,
    args,
    recovery_requested: bool,
) -> bool:
    failed = archive / "stage.failed"
    if not recovery_requested or not failed.exists():
        return False

    stage_name = stage["name"]
    if (archive / "stage.done").exists():
        raise RuntimeError(
            f"Cannot recover {stage_name}: archive contains both stage.done "
            "and stage.failed."
        )
    reason = failed.read_text(errors="replace").strip()
    if reason.lower() != "thermostating problem":
        raise RuntimeError(
            f"Cannot recover {stage_name}: stage.failed reason is "
            f"'{reason}', not 'thermostating problem'."
        )
    if not thermostat_warning_allowed(
        args.thermostat_warning_policy,
        stage_name,
    ):
        raise RuntimeError(
            f"Cannot recover {stage_name}: thermostat warning policy "
            f"'{args.thermostat_warning_policy}' does not accept this stage."
        )

    required = ["xtb.trj", "mdrestart", "xtbmdok", f"{stage_name}.out"]
    if stage["restart"]:
        required.append("mdrestart.input")
    missing = [name for name in required if not (archive / name).exists()]
    if missing:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived outputs are incomplete "
            f"(missing: {', '.join(missing)})."
        )

    log_path = archive / f"{stage_name}.out"
    validation = inspect_md_log(log_path)
    if validation["fatal_patterns"]:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived log contains fatal "
            f"pattern(s): {', '.join(validation['fatal_patterns'])}."
        )
    if not validation["thermostating_problem"]:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived log does not contain "
            "'thermostating problem'."
        )
    if not validation["normal_exit_of_md"]:
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived warning is not followed "
            "by 'normal exit of md()'."
        )

    manifest_path = archive / "stage_manifest.json"
    had_stage_manifest = manifest_path.exists()
    if had_stage_manifest:
        try:
            recorded = json.loads(manifest_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise RuntimeError(
                f"Cannot recover {stage_name}: unreadable "
                "stage_manifest.json."
            ) from exc
        mismatched = configuration_mismatches(configuration, recorded)
        if mismatched:
            raise RuntimeError(
                f"Cannot recover {stage_name}: archived configuration is "
                f"incompatible (fields: {', '.join(mismatched)})."
            )
        if recorded.get("returncode") != 0:
            raise RuntimeError(
                f"Cannot recover {stage_name}: archived return code is "
                f"{recorded.get('returncode')}."
            )

    input_restart_sha256 = configuration["input_restart_sha256"]
    if stage["restart"]:
        archived_input_hash = file_sha256(archive / "mdrestart.input")
        if archived_input_hash != input_restart_sha256:
            raise RuntimeError(
                f"Cannot recover {stage_name}: archived input restart does "
                "not match the currently valid previous stage."
            )

    output_restart_sha256 = validate_output_restart(
        archive / "mdrestart",
        input_restart_sha256,
        expected_atoms=restart_atom_count(archive.parent.parent),
    )
    if (
        manifest_path.exists()
        and recorded.get("output_restart_sha256")
        != output_restart_sha256
    ):
        raise RuntimeError(
            f"Cannot recover {stage_name}: output restart hash is "
            "inconsistent with stage_manifest.json."
        )
    archived_input = archive / f"{stage_name}.inp"
    if (
        archived_input.exists()
        and file_sha256(archived_input)
        != configuration["xtb_input_sha256"]
    ):
        raise RuntimeError(
            f"Cannot recover {stage_name}: archived xTB input hash is "
            "inconsistent with historical provenance."
        )
    thermal_result = md_thermal_result(
        stage,
        validation,
        args.thermostat_warning_policy,
        warning_accepted=True,
    )
    write_md_stage_manifest(
        archive,
        configuration,
        output_restart_sha256,
        thermal_result,
        recovery={
            "promoted_from_stage_failed": True,
            "original_failure_reason": reason,
            "configuration_source": (
                "existing stage_manifest.json"
                if had_stage_manifest
                else (
                    "legacy manifest.json plus deterministic stage input "
                    "and archived restart/output chain"
                )
            ),
            "returncode_evidence": (
                "recorded returncode == 0"
                if had_stage_manifest
                else (
                    "legacy pipeline created a thermostat-only stage.failed "
                    "marker only after returncode == 0"
                )
            ),
        },
    )
    mark_stage_done(archive)
    print(
        f"  RECOVER {stage_name}: archived stage had thermostat warning "
        "only; outputs passed integrity checks and current policy accepts it."
    )
    return True


def recover_failed_md_stage(
    stage: dict,
    archive: Path,
    configuration: dict,
    args,
    recovery_requested: bool,
) -> bool:
    if recover_legacy_restart_validation_failure(
        stage,
        archive,
        configuration,
        args,
        recovery_requested,
    ):
        return True
    return recover_thermostat_warning_stage(
        stage,
        archive,
        configuration,
        args,
        recovery_requested,
    )


def prepare_md_stage_attempt(
    replica_dir: Path,
    stage_name: str,
    force: bool,
    expected_configuration: dict,
    warning_policy: str,
):
    archive = stage_archive_dir(replica_dir, stage_name)
    done = archive / "stage.done"
    failed = archive / "stage.failed"

    if force:
        if archive.exists() and any(archive.iterdir()):
            print(
                f"  WARNING replacing archived {stage_name} outputs (--force)"
            )
            shutil.rmtree(archive)
        return archive, False

    if done.exists():
        validate_completed_md_stage(
            archive,
            expected_configuration,
            warning_policy,
        )
        return archive, True
    if failed.exists():
        raise RuntimeError(
            f"Previous failed stage exists for {stage_name}; inspect "
            "outputs or rerun with --force."
        )
    if archive.exists() and any(archive.iterdir()):
        raise RuntimeError(
            f"Incomplete archived outputs exist for {stage_name} without a "
            "valid status marker; inspect them or rerun with --force."
        )
    return archive, False


def archive_failed_md_stage(
    replica_dir: Path,
    stage_name: str,
    log_path: Path,
    reason: str,
):
    archive_stage_outputs(replica_dir, stage_name, log_path)
    archive = stage_archive_dir(replica_dir, stage_name)
    mark_stage_failed(archive, reason)


def run_relaxation(replica_dir: Path, args, env):
    manifest = json.loads((replica_dir / "manifest.json").read_text())
    config = manifest["relaxation"]
    n_solute = config["n_solute_atoms"]
    n_total = manifest["geometry_after_centering"]["n_atoms"]
    relaxed_pdb = replica_dir / "system_relaxed.pdb"
    archive = stage_archive_dir(replica_dir, "00_relax")
    done = archive / "stage.done"
    signature = relaxation_signature(replica_dir, manifest, args)

    if done.exists() and not args.force:
        previous = manifest.get("relaxation_result", {})
        previous_configuration = dict(
            previous.get("configuration", {})
        )
        # Results produced before --relax-engine existed necessarily used
        # xTB's default engine selection, i.e. the current "auto" mode.
        if "optimization_engine" not in previous_configuration:
            previous_configuration["optimization_engine"] = "auto"
            previous["configuration"] = previous_configuration
            update_manifest(replica_dir, "relaxation_result", previous)
        if previous_configuration != signature:
            raise RuntimeError(
                "Completed 00_relax is incompatible with the current geometry "
                "or settings. Use --force to rebuild it."
            )
        if not valid_pdb(relaxed_pdb, n_total):
            raise RuntimeError(
                f"{done} exists, but system_relaxed.pdb is absent or invalid. "
                "Use --force to rebuild 00_relax."
            )
        displacement = solute_displacement(
            replica_dir / "system_centered.pdb", relaxed_pdb, n_solute
        )
        if displacement["solute_max_displacement_A"] > 1.0e-3:
            raise RuntimeError(
                "Completed relaxation has moved fixed solute atoms by "
                f"{displacement['solute_max_displacement_A']:.6f} A."
            )
        if args.resume or args.start_stage is not None:
            print("  REUSE 00_relax")
        else:
            print(
                f"  SKIP {replica_dir.parent.name}/{replica_dir.name} "
                "00_relax: already completed"
            )
        return

    if relaxed_pdb.exists() and not args.force:
        raise RuntimeError(
            f"{relaxed_pdb} exists without a reusable completed stage. "
            "Refusing to overwrite it; inspect it or use --force."
        )

    archive.mkdir(parents=True, exist_ok=True)
    if args.force and done.exists():
        print(
            f"  WARNING replacing archived 00_relax outputs for "
            f"{replica_dir.parent.name}/{replica_dir.name} (--force)"
        )
    if args.force and relaxed_pdb.exists():
        relaxed_pdb.unlink()
    for name in [
        "00_relax.out",
        "xtbopt.pdb", "xtbopt.xyz", "xtbopt.coord", "xtbopt.log",
        "xtblast.pdb", "xtblast.xyz", "xtblast.coord",
        "NOT_CONVERGED",
    ]:
        path = replica_dir / name
        if path.exists():
            path.unlink()
        archived = archive / name
        if archived.exists():
            archived.unlink()
    if done.exists():
        done.unlink()

    cmd = xtb_command(
        args, "system_centered.pdb", "00_relax.inp", optimize=True
    )
    log_path = replica_dir / "00_relax.out"
    print(f"  RELAX {replica_dir.parent.name}/{replica_dir.name}")
    print(f"       fixed solute atoms  : 1-{n_solute}")
    print(f"       mobile solvent atoms: {n_solute + 1}-{n_total}")
    print(f"       optimization level : {args.relax_level}")
    print(f"       optimization engine: {args.relax_engine}")
    print(f"       max cycles          : {args.relax_cycles}")
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
            f"xTB failed at 00_relax for {replica_dir.parent.name}/"
            f"{replica_dir.name}. See {log_path}"
        )

    not_converged_file = (replica_dir / "NOT_CONVERGED").exists()
    diagnostics = relaxation_diagnostics(log_path, not_converged_file)
    diagnostics["configuration"] = signature
    materialize_relaxed_pdb(replica_dir, n_total)
    displacement = solute_displacement(
        replica_dir / "system_centered.pdb", relaxed_pdb, n_solute
    )
    diagnostics.update(displacement)
    if displacement["solute_max_displacement_A"] > 1.0e-3:
        raise RuntimeError(
            "00_relax changed fixed solute coordinates: maximum displacement "
            f"{displacement['solute_max_displacement_A']:.6f} A exceeds "
            "the 0.001 A tolerance. MD was not started."
        )

    shutil.copy2(log_path, archive / "00_relax.out")
    for name in [
        "xtbopt.pdb", "xtbopt.xyz", "xtbopt.coord", "xtbopt.log",
        "xtblast.pdb", "xtblast.xyz", "xtblast.coord",
        "NOT_CONVERGED",
    ]:
        path = replica_dir / name
        if path.exists():
            shutil.copy2(path, archive / name)
    update_manifest(replica_dir, "relaxation_result", diagnostics)
    done.write_text("ok\n")

    if diagnostics["converged"] is False:
        print(
            "  WARNING 00_relax did not converge fully; using the last valid "
            "geometry for MD."
        )
    else:
        print("  OK   00_relax")


def validate_historical_md_prefix(
    replica_dir: Path,
    stop_index: int,
    args,
    *,
    allow_recovery: bool = True,
) -> dict | None:
    """Validate/reuse historical MD stages preceding an explicit start."""
    previous_output_sha256 = None
    previous_configuration = None

    for stage_index in range(stop_index):
        stage, configuration = historical_md_stage_configuration(
            replica_dir,
            stage_index,
            previous_output_sha256,
        )
        validate_historical_stage_compatibility(
            configuration,
            STAGES[stage_index],
            args,
        )
        archive = stage_archive_dir(replica_dir, stage["name"])
        recover_failed_md_stage(
            stage,
            archive,
            configuration,
            args,
            recovery_requested=allow_recovery,
        )
        recorded = validate_completed_md_stage(
            archive,
            configuration,
            args.thermostat_warning_policy,
        )
        previous_output_sha256 = recorded["output_restart_sha256"]
        previous_configuration = configuration
        print(f"  REUSE {stage['name']}")

    return previous_configuration


def validate_historical_stage_compatibility(
    historical_configuration: dict,
    current_stage: dict,
    args,
) -> dict | None:
    """Require current E2 settings except for an explicit old duration."""
    expected = {
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "temp_K": current_stage["temp"],
        "step_fs": MD_STEP_FS,
        "dump_fs": MD_DUMP_FS,
        "hmass": 1,
        "shake": 0,
        "sccacc": 1.0,
        "restart": current_stage["restart"],
    }
    mismatched = configuration_mismatches(
        expected,
        historical_configuration,
    )
    if mismatched:
        raise RuntimeError(
            f"Cannot reuse historical {current_stage['name']}: "
            "provenance is incompatible with the requested E2 settings "
            f"(fields: {', '.join(mismatched)})."
        )

    historical_time = historical_configuration["time_ps"]
    if historical_time == current_stage["time"]:
        return None

    return {
        "previous_stage": current_stage["name"],
        "historical_time_ps": historical_time,
        "current_default_time_ps": current_stage["time"],
        "explicit_start_stage_override": True,
    }


def validate_start_predecessor_compatibility(
    previous_configuration: dict,
    previous_stage: dict,
    start_stage_name: str,
    args,
) -> dict | None:
    """Report an explicit historical-duration exception for the predecessor."""
    override = validate_historical_stage_compatibility(
        previous_configuration,
        previous_stage,
        args,
    )
    if override is None:
        return None

    previous_time = override["historical_time_ps"]
    print(
        f"  WARNING: starting {start_stage_name} from an existing "
        f"{previous_stage['name']} generated with a previous stage-duration "
        f"configuration ({previous_time:g} ps; current default "
        f"{previous_stage['time']:g} ps)."
    )
    return override


def validate_extended_stage_preconditions(
    replica_dir: Path,
    previous_configuration: dict,
    args,
) -> dict:
    """Validate the exact 04 -> 05 continuation before any stage mutation."""
    previous_name = EXTENDED_STAGE["restart_from"]
    previous_archive = stage_archive_dir(replica_dir, previous_name)
    done = previous_archive / "stage.done"
    failed = previous_archive / "stage.failed"
    if not done.is_file():
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: missing validated "
            f"{previous_name}/stage.done at {done}."
        )
    if failed.exists():
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: unresolved "
            f"{previous_name}/stage.failed exists at {failed}."
        )

    try:
        manifest = json.loads((replica_dir / "manifest.json").read_text())
        n_atoms = manifest["geometry_after_centering"]["n_atoms"]
        relaxation = manifest["relaxation"]
        n_solute = relaxation["n_solute_atoms"]
        n_waters = manifest["packing"]["water_count"]
        wall_radius_bohr = manifest["wall"]["radius_bohr"]
    except (OSError, json.JSONDecodeError, KeyError, TypeError) as exc:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: replica composition/wall "
            "metadata is missing or invalid."
        ) from exc

    for label, value in (
        ("total atom count", n_atoms),
        ("solute atom count", n_solute),
        ("water count", n_waters),
    ):
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise RuntimeError(
                f"Cannot start {EXTENDED_STAGE_NAME}: invalid {label} "
                f"{value!r} in manifest.json."
            )
    if n_solute + 3 * n_waters != n_atoms:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: manifest composition is "
            f"inconsistent ({n_solute} solute atoms + {n_waters} waters x 3 "
            f"!= {n_atoms} total atoms)."
        )

    expected_previous = {
        "stage": previous_name,
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "temp_K": EXTENDED_STAGE["temp"],
        "time_ps": STAGES[CORE_STAGE_COUNT - 1]["time"],
        "step_fs": MD_STEP_FS,
        "dump_fs": MD_DUMP_FS,
        "hmass": 1,
        "shake": 0,
        "sccacc": 1.0,
        "restart": True,
        "wall_radius_bohr": wall_radius_bohr,
    }
    mismatched = configuration_mismatches(
        expected_previous,
        previous_configuration,
    )
    if mismatched:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: {previous_name} does not "
            "match the inherited continuation protocol (fields: "
            f"{', '.join(mismatched)})."
        )

    restart_path = previous_archive / "mdrestart"
    restart_metadata = parse_mdrestart(restart_path, n_atoms)
    restart_sha256 = file_sha256(restart_path)
    stage_manifest_path = previous_archive / "stage_manifest.json"
    try:
        previous_stage_manifest = json.loads(stage_manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: missing or invalid "
            f"{stage_manifest_path}."
        ) from exc
    if previous_stage_manifest.get("output_restart_sha256") != restart_sha256:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: {previous_name} restart "
            "SHA-256 does not match its stage manifest."
        )

    previous_trajectory = previous_archive / "xtb.trj"
    if not previous_trajectory.is_file() or previous_trajectory.stat().st_size == 0:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: missing or empty "
            f"{previous_name} trajectory at {previous_trajectory}."
        )
    duration_ratio = (
        EXTENDED_STAGE["time"] / STAGES[CORE_STAGE_COUNT - 1]["time"]
    )
    estimated_trajectory_bytes = math.ceil(
        previous_trajectory.stat().st_size * duration_ratio
    )
    required_free_bytes = math.ceil(
        estimated_trajectory_bytes * 1.25
        + restart_path.stat().st_size * 3
    )
    free_bytes = shutil.disk_usage(replica_dir).free
    if free_bytes < required_free_bytes:
        raise RuntimeError(
            f"Cannot start {EXTENDED_STAGE_NAME}: only {free_bytes} bytes "
            f"free, below the conservative estimate of {required_free_bytes} "
            "bytes required for trajectory, restart copies, and margin."
        )

    return {
        "previous_stage": previous_name,
        "restart_path": restart_path,
        "restart_sha256": restart_sha256,
        "restart_size_bytes": restart_path.stat().st_size,
        "restart_validation": restart_metadata,
        "n_atoms": n_atoms,
        "n_solute_atoms": n_solute,
        "n_waters": n_waters,
        "wall_radius_bohr": wall_radius_bohr,
        "estimated_trajectory_bytes": estimated_trajectory_bytes,
        "required_free_bytes": required_free_bytes,
        "free_bytes": free_bytes,
        "previous_stage_threads": previous_configuration.get("threads"),
        "requested_stage_threads": args.threads,
        "thread_count_changed": (
            previous_configuration.get("threads") != args.threads
        ),
    }


def print_extended_dry_run(replica_dir: Path, args, preflight: dict):
    stage = EXTENDED_STAGE
    input_name = f"{EXTENDED_STAGE_NAME}.inp"
    geometry_name = "system_centered.pdb"
    cmd = xtb_command(args, geometry_name, input_name)
    print(f"\nDRY-RUN {replica_dir.parent.name}/{replica_dir.name}")
    print(f"  System: {replica_dir.parent.name}")
    print(f"  Replica: {replica_dir.name}")
    print(f"  Stage: {EXTENDED_STAGE_NAME}")
    print(f"  Restart: {preflight['restart_path']}")
    print(f"  Restart SHA256: {preflight['restart_sha256']}")
    print(f"  Atoms: {preflight['n_atoms']}")
    print(f"  Solute atoms: {preflight['n_solute_atoms']}")
    print(f"  Waters: {preflight['n_waters']}")
    print(f"  Temperature: {stage['temp']:.2f} K")
    print(f"  Duration: {stage['time']:.1f} ps")
    print(f"  Time step: {MD_STEP_FS:.1f} fs")
    print(f"  Steps: {md_step_count(stage)}")
    print(f"  Dump interval: {MD_DUMP_FS:.1f} fs")
    print(
        f"  Threads: {args.threads} "
        f"(previous stage: {preflight['previous_stage_threads']})"
    )
    print("  Randomize velocities: no")
    print("  Solvation rebuilt: no")
    print(f"  Command: {' '.join(cmd)}")
    print(f"  Planned input {input_name}:")
    print(md_input(stage, preflight["wall_radius_bohr"]).rstrip())


def md_stop_index_for_request(args) -> int:
    """Keep stage 05 opt-in while all existing execution modes stop at 04."""
    if args.start_stage == EXTENDED_STAGE_NAME:
        return EXTENDED_STAGE_INDEX
    return CORE_STAGE_COUNT - 1


def run_replica(replica_dir: Path, args):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    env.setdefault("MKL_NUM_THREADS", str(args.threads))
    env.setdefault("OMP_STACKSIZE", "4G")

    existing_only = args.resume or args.start_stage is not None
    if existing_only:
        try:
            replica_manifest = json.loads(
                (replica_dir / "manifest.json").read_text()
            )
            relaxation_enabled = replica_manifest["relaxation"]["enabled"]
        except (OSError, json.JSONDecodeError, KeyError) as exc:
            raise RuntimeError(
                f"Cannot resume {replica_dir}: missing or invalid manifest.json."
            ) from exc
    else:
        relaxation_enabled = not args.skip_relax

    md_start_index = 0
    md_stop_index = md_stop_index_for_request(args)
    start_stage_override = None
    extended_preflight = None

    if args.resume:
        if relaxation_enabled:
            validate_historical_relaxation(replica_dir)
    elif args.start_stage is not None:
        if args.start_stage == "00_relax":
            if not relaxation_enabled:
                raise RuntimeError(
                    "Cannot start at 00_relax: this replica manifest records "
                    "relaxation as disabled."
                )
            run_relaxation(replica_dir, args, env)
        else:
            md_start_index = next(
                index for index, item in enumerate(STAGES)
                if item["name"] == args.start_stage
            )
            if relaxation_enabled:
                validate_historical_relaxation(replica_dir)

            previous_configuration = validate_historical_md_prefix(
                replica_dir,
                md_start_index,
                args,
                allow_recovery=not args.dry_run,
            )
            if previous_configuration is not None:
                previous_stage = STAGES[md_start_index - 1]
                start_stage_override = (
                    validate_start_predecessor_compatibility(
                        previous_configuration,
                        previous_stage,
                        args.start_stage,
                        args,
                    )
                )
            if args.start_stage == EXTENDED_STAGE_NAME:
                extended_preflight = validate_extended_stage_preconditions(
                    replica_dir,
                    previous_configuration,
                    args,
                )
                if args.dry_run:
                    print_extended_dry_run(
                        replica_dir,
                        args,
                        extended_preflight,
                    )
                    return
                ensure_extended_stage_registered(replica_dir)
    elif relaxation_enabled:
        run_relaxation(replica_dir, args, env)

    for i, stage in enumerate(STAGES):
        if i < md_start_index or i > md_stop_index:
            continue

        name = stage["name"]
        input_name = f"{name}.inp"
        geometry_name = (
            "system_relaxed.pdb"
            if i == 0 and relaxation_enabled
            else "system_centered.pdb"
        )
        validate_current_md_input(
            replica_dir,
            stage,
            create_missing=(args.start_stage == name),
        )
        input_restart_sha256 = None
        if stage["restart"]:
            previous_restart = (
                stage_archive_dir(replica_dir, STAGES[i - 1]["name"])
                / "mdrestart"
            )
            if not previous_restart.exists():
                raise RuntimeError(
                    f"Previous restart not found: {previous_restart}. "
                    "Run/complete the previous stage first."
                )
            input_restart_sha256 = file_sha256(previous_restart)

        stage_configuration = md_stage_configuration(
            replica_dir=replica_dir,
            stage=stage,
            geometry_name=geometry_name,
            input_name=input_name,
            input_restart_sha256=input_restart_sha256,
            args=args,
        )
        archive = stage_archive_dir(replica_dir, name)
        recovery_requested = (
            (args.resume or args.start_stage is not None)
            and not args.force
        )
        if recovery_requested and (archive / "stage.failed").exists():
            historical_stage, historical_configuration = (
                historical_md_stage_configuration(
                    replica_dir,
                    i,
                    input_restart_sha256,
                )
            )
            mismatched = configuration_mismatches(
                stage_configuration,
                historical_configuration,
            )
            if mismatched:
                raise RuntimeError(
                    f"Cannot recover {name}: historical calculation is "
                    "incompatible with the requested protocol "
                    f"(fields: {', '.join(mismatched)})."
                )
        else:
            historical_stage = stage

        if recover_failed_md_stage(
            historical_stage,
            archive,
            stage_configuration,
            args,
            recovery_requested=recovery_requested,
        ):
            continue

        archive, already_done = prepare_md_stage_attempt(
            replica_dir,
            name,
            args.force,
            stage_configuration,
            args.thermostat_warning_policy,
        )
        if already_done:
            if recovery_requested:
                print(f"  REUSE {name}")
            else:
                print(
                    f"  SKIP {replica_dir.parent.name}/"
                    f"{replica_dir.name} {name}: already completed"
                )
            continue

        if i == md_start_index and start_stage_override is not None:
            stage_configuration["start_stage_override"] = start_stage_override

        for transient in [
            "xtb.trj",
            "xtb-trj.pdb",
            "xtbmdok",
            "mdrestart.input",
        ]:
            p = replica_dir / transient
            if p.exists():
                p.unlink()

        restore_restart_from_previous(
            replica_dir,
            i,
        )
        if stage["restart"]:
            restart_path = replica_dir / "mdrestart"
            copied_input_hash = file_sha256(restart_path)
            if copied_input_hash != input_restart_sha256:
                raise RuntimeError(
                    f"Input mdrestart changed while preparing {name}; "
                    "refusing to run an inconsistent restart chain."
                )
            shutil.copy2(
                restart_path,
                replica_dir / "mdrestart.input",
            )
            if name == "04_298K_screen":
                print(
                    "  START 04_298K_screen from restart of "
                    "03_298K_equil"
                )
                print("       previous stage: 03_298K_equil")
                print(
                    f"       input restart SHA256: {copied_input_hash}"
                )
            elif name == EXTENDED_STAGE_NAME:
                print(
                    f"  START {EXTENDED_STAGE_NAME} from final restart of "
                    f"{EXTENDED_STAGE['restart_from']}"
                )
                print(
                    f"       input restart SHA256: {copied_input_hash}"
                )

        log_path = replica_dir / f"{name}.out"
        cmd = xtb_command(args, geometry_name, input_name)
        started_at = datetime.now(timezone.utc).isoformat()
        running_metadata = {
            "stage": name,
            "status": "running",
            "started_at": started_at,
            "command": cmd,
            "hostname": socket.gethostname(),
        }
        mark_stage_running(archive, running_metadata)

        print(
            f"  RUN  {replica_dir.parent.name}/"
            f"{replica_dir.name} {name}"
        )
        print(f"       {' '.join(cmd)}")

        with log_path.open("w") as log:
            if name == EXTENDED_STAGE_NAME:
                continuation_message = (
                    "Continuing 04_298K_screen from its final restart: "
                    f"{stage['time']:.1f} ps at {stage['temp']:.2f} K, "
                    f"{md_step_count(stage)} steps, dt = {MD_STEP_FS:.1f} fs."
                )
                print(f"       {continuation_message}")
                log.write(continuation_message + "\n")
                log.flush()
            result = subprocess.run(
                cmd,
                cwd=replica_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                text=True,
            )

        validation = inspect_md_log(log_path)
        if result.returncode != 0:
            archive_failed_md_stage(
                replica_dir,
                name,
                log_path,
                f"xTB return code {result.returncode}",
            )
            raise RuntimeError(
                f"xTB failed at {name} for "
                f"{replica_dir.parent.name}/{replica_dir.name}. "
                f"Outputs were archived for inspection. See {log_path}"
            )

        if validation["fatal_patterns"]:
            patterns = ", ".join(validation["fatal_patterns"])
            archive_failed_md_stage(
                replica_dir,
                name,
                log_path,
                f"fatal log pattern(s): {patterns}",
            )
            raise RuntimeError(
                f"Fatal MD pattern(s) detected at {name}: {patterns}. "
                "Outputs were archived for inspection."
            )

        if not (replica_dir / "xtb.trj").exists():
            archive_failed_md_stage(
                replica_dir, name, log_path, "missing xtb.trj"
            )
            raise RuntimeError(
                f"{name} ended without xtb.trj for "
                f"{replica_dir.parent.name}/{replica_dir.name}. "
                f"Outputs were archived for inspection. See {log_path}"
            )

        trajectory_integrity = None
        if name == EXTENDED_STAGE_NAME:
            try:
                trajectory_integrity = validate_xtb_trajectory(
                    replica_dir / "xtb.trj",
                    expected_atoms=restart_atom_count(replica_dir),
                    expected_frames=md_expected_frame_count(stage),
                    require_velocities=True,
                )
            except RuntimeError as exc:
                archive_failed_md_stage(
                    replica_dir,
                    name,
                    log_path,
                    str(exc),
                )
                raise RuntimeError(
                    f"{name} produced an invalid trajectory ({exc}). "
                    "Outputs were archived for inspection."
                ) from exc

        restart_path = replica_dir / "mdrestart"
        try:
            output_restart_sha256 = validate_output_restart(
                restart_path,
                input_restart_sha256,
                expected_atoms=restart_atom_count(replica_dir),
            )
        except RuntimeError as exc:
            archive_failed_md_stage(
                replica_dir, name, log_path, str(exc)
            )
            if "byte-identical" in str(exc):
                raise RuntimeError(
                    f"{name} ended without producing a new mdrestart. "
                    "The output restart is byte-identical to the input "
                    "restart. Outputs were archived for inspection."
                ) from exc
            raise RuntimeError(
                f"{name} ended with invalid mdrestart ({exc}) for "
                f"{replica_dir.parent.name}/{replica_dir.name}. "
                f"Outputs were archived for inspection. See {log_path}"
            ) from exc

        if not (replica_dir / "xtbmdok").exists():
            archive_failed_md_stage(
                replica_dir, name, log_path, "missing xtbmdok"
            )
            raise RuntimeError(
                f"{name} ended without xtbmdok; refusing to mark the stage "
                "as complete. Outputs were archived for inspection."
            )

        warning_accepted = (
            validation["thermostating_problem"]
            and validation["normal_exit_of_md"]
            and thermostat_warning_allowed(
                args.thermostat_warning_policy,
                name,
            )
        )
        thermal_result = md_thermal_result(
            stage,
            validation,
            args.thermostat_warning_policy,
            warning_accepted,
        )
        runtime_metadata = None
        if name == EXTENDED_STAGE_NAME:
            if file_sha256(extended_preflight["restart_path"]) != (
                extended_preflight["restart_sha256"]
            ):
                archive_failed_md_stage(
                    replica_dir,
                    name,
                    log_path,
                    f"archived {EXTENDED_STAGE['restart_from']} restart "
                    "changed during continuation",
                )
                raise RuntimeError(
                    f"Archived {EXTENDED_STAGE['restart_from']} restart "
                    "changed during the continuation; refusing to certify "
                    f"{name}."
                )
            xtb_version = extract_xtb_version(
                log_path,
                replica_dir / "xtb.trj",
            )
            if xtb_version is None:
                archive_failed_md_stage(
                    replica_dir,
                    name,
                    log_path,
                    "xTB version could not be determined from log/trajectory",
                )
                raise RuntimeError(
                    f"{name} outputs do not report the xTB version; refusing "
                    "to mark the continuation as provenance-complete."
                )
            archived_trajectory_path = str(
                Path("stages") / name / "xtb.trj"
            )
            trajectory_integrity = dict(trajectory_integrity)
            trajectory_integrity["path"] = archived_trajectory_path
            preflight_record = {
                key: str(value) if isinstance(value, Path) else value
                for key, value in extended_preflight.items()
            }
            runtime_metadata = {
                "started_at": started_at,
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "command": cmd,
                "xtb_version": xtb_version,
                "hostname": socket.gethostname(),
                "status": "completed",
                "preflight": preflight_record,
                "trajectory_integrity": trajectory_integrity,
                "output_trajectory_sha256": trajectory_integrity["sha256"],
                "output_paths": {
                    "input": str(Path("stages") / name / input_name),
                    "log": str(Path("stages") / name / log_path.name),
                    "trajectory": archived_trajectory_path,
                    "input_restart_copy": str(
                        Path("stages") / name / "mdrestart.input"
                    ),
                    "output_restart": str(
                        Path("stages") / name / "mdrestart"
                    ),
                },
            }
        if (
            validation["thermostating_problem"]
            and not warning_accepted
        ):
            reason = (
                "thermostating problem"
                if validation["normal_exit_of_md"]
                else "thermostating problem without normal exit of md()"
            )
            write_md_stage_manifest(
                archive,
                stage_configuration,
                output_restart_sha256,
                thermal_result,
                runtime_metadata=runtime_metadata,
            )
            archive_stage_outputs(
                replica_dir,
                name,
                log_path,
            )
            mark_stage_failed(archive, reason)
            next_stage = (
                STAGES[i + 1]["name"]
                if i + 1 <= md_stop_index
                else "pipeline completion"
            )
            if not validation["normal_exit_of_md"]:
                raise RuntimeError(
                    f"{name} reported 'thermostating problem' without "
                    "'normal exit of md()'. Outputs were archived; refusing "
                    f"to continue to {next_stage}."
                )
            raise RuntimeError(
                f"{name} completed numerically but xTB reported "
                "'thermostating problem'. Outputs passed integrity checks, "
                f"but policy '{args.thermostat_warning_policy}' refuses to "
                f"continue to {next_stage}."
            )

        write_md_stage_manifest(
            archive,
            stage_configuration,
            output_restart_sha256,
            thermal_result,
            runtime_metadata=runtime_metadata,
        )
        archive_stage_outputs(
            replica_dir,
            name,
            log_path,
        )

        mark_stage_done(archive)
        if warning_accepted:
            next_stage = (
                STAGES[i + 1]["name"]
                if i + 1 <= md_stop_index
                else "pipeline completion"
            )
            print(
                f"  WARNING {name}: xTB reported 'thermostating problem', "
                "but MD exited normally and all restart/output integrity "
                f"checks passed. Continuing to {next_stage}."
            )
        else:
            print(f"  OK   {name}")

    final_archived = (
        stage_archive_dir(
            replica_dir,
            STAGES[md_stop_index]["name"],
        )
        / "mdrestart"
    )

    if final_archived.exists():
        final_name = (
            "mdrestart_final"
            if md_stop_index < EXTENDED_STAGE_INDEX
            else f"mdrestart_final_{EXTENDED_STAGE_NAME}"
        )
        shutil.copy2(
            final_archived,
            replica_dir / final_name,
        )


# ---------------------------------------------------------------------------
# Independent CO2 shell-screening workflow
# ---------------------------------------------------------------------------

def co2_md_stages(args) -> list[dict]:
    """Build the independent CO2 MD protocol from explicit CLI values."""
    equilibration = {
        "name": CO2_EQUIL_STAGE,
        "temp": 298.15,
        "time": args.co2_equil_time_ps,
        "dump_fs": args.co2_equil_dump_fs,
        "restart": False,
        "velocities_reinitialized": True,
        "continuation": False,
        "source_stage": CO2_ACCOMMODATION_STAGE,
        "equilibration_time_ps": args.co2_equil_time_ps,
        "production_time_ps": 0.0,
        "cumulative_production_time_ps": 0.0,
        "purpose": (
            "Short NVT equilibration after CO2 insertion and structural "
            "accommodation. New velocities are initialized at 298.15 K."
        ),
    }
    screen = {
        "name": CO2_SCREEN_STAGE,
        "temp": 298.15,
        "time": args.co2_screen_time_ps,
        "dump_fs": args.co2_production_dump_fs,
        "restart": True,
        "restart_from": CO2_EQUIL_STAGE,
        "velocities_reinitialized": False,
        "continuation": True,
        "continuation_of": CO2_EQUIL_STAGE,
        "equilibration_time_ps": args.co2_equil_time_ps,
        "production_time_ps": args.co2_screen_time_ps,
        "cumulative_production_time_ps": args.co2_screen_time_ps,
        "purpose": "Initial unbiased production/screening trajectory with CO2.",
    }
    extended = {
        "name": CO2_EXTENDED_STAGE,
        "temp": 298.15,
        "time": args.co2_extended_time_ps,
        "dump_fs": args.co2_production_dump_fs,
        "restart": True,
        "restart_from": CO2_SCREEN_STAGE,
        "velocities_reinitialized": False,
        "continuation": True,
        "continuation_of": CO2_SCREEN_STAGE,
        "equilibration_time_ps": args.co2_equil_time_ps,
        "production_time_ps": args.co2_extended_time_ps,
        "cumulative_production_time_ps": (
            args.co2_screen_time_ps + args.co2_extended_time_ps
        ),
        "purpose": "Opt-in extended unbiased production trajectory with CO2.",
    }
    for stage in (equilibration, screen, extended):
        stage["steps"] = md_step_count(stage)
        stage["expected_frames"] = md_expected_frame_count(stage)
    return [equilibration, screen, extended]


def co2_md_input(stage: dict, wall_radius_bohr: float) -> str:
    """Create an unconstrained CO2 MD input with the stage-specific dump."""
    return md_input(stage, wall_radius_bohr)


def co2_md_command(args, stage: dict) -> list[str]:
    command = [
        args.xtb,
        "system_CO2_accommodated.pdb",
        "--gfn", str(args.gfn),
        "--chrg", str(args.charge),
        "--uhf", str(args.uhf),
        "--md",
        "--input", f"{stage['name']}.inp",
    ]
    if args.alpb:
        command += ["--alpb", args.alpb]
    return command


def co2_sampling_metadata(dump_fs: float, velocities_present: bool) -> dict:
    sampling_frequency_hz = 1.0 / (dump_fs * 1.0e-15)
    return {
        "trajectory_sampling_interval_fs": dump_fs,
        "sampling_frequency_Hz": sampling_frequency_hz,
        "nyquist_wavenumber_cm-1": (
            sampling_frequency_hz / (2.0 * SPEED_OF_LIGHT_CM_S)
        ),
        "spectroscopy_sampling_ready": (
            dump_fs <= 2.0 and velocities_present
        ),
        "spectroscopy_sampling_note": (
            "Temporal sampling suitable for subsequent vibrational "
            "power-spectrum / VDOS analysis; this is not an IR spectrum."
        ),
    }


def co2_execution_resources(args) -> dict:
    parallel_jobs = getattr(args, "co2_parallel_jobs", 1)
    detected_cpus = os.cpu_count()
    maximum_requested = parallel_jobs * args.threads
    return {
        "co2_parallel_jobs": parallel_jobs,
        "xtb_threads_per_job": args.threads,
        "maximum_requested_cpus": maximum_requested,
        "detected_cpu_count": detected_cpus,
        "oversubscription_warning": (
            detected_cpus is not None and maximum_requested > detected_cpus
        ),
    }


def _validated_pdb_atoms(path: Path, description: str):
    try:
        atoms = pdb_atoms(path)
    except (OSError, ValueError) as exc:
        raise RuntimeError(
            f"Cannot read {description} PDB {path}: {exc}"
        ) from exc
    for index, atom in enumerate(atoms, start=1):
        if not all(math.isfinite(value) for value in atom["xyz"]):
            raise RuntimeError(
                f"{description} PDB {path} has non-finite coordinates at "
                f"atom {index}."
            )
    return atoms


def _element_sequence(atoms) -> list[str]:
    return [atom["element"] for atom in atoms]


def validate_co2_source_pdb(
    source_pdb: Path,
    n_solute_atoms: int,
    expected_solute_pdb: Path | None = None,
) -> dict:
    """Validate a full-droplet CO2 source without changing its atom order."""
    atoms = _validated_pdb_atoms(source_pdb, "CO2 source")
    if n_solute_atoms <= 0:
        raise RuntimeError("--co2-solute-atoms must be positive.")
    if len(atoms) <= n_solute_atoms:
        raise RuntimeError(
            f"CO2 source {source_pdb} must contain the {n_solute_atoms} "
            "solute atoms followed by explicit water."
        )

    elements = _element_sequence(atoms)
    zinc_indices = [i for i, element in enumerate(elements) if element == "Zn"]
    if len(zinc_indices) != 1:
        raise RuntimeError(
            f"CO2 source {source_pdb} must contain exactly one Zn atom; "
            f"found {len(zinc_indices)}."
        )
    zinc_index = zinc_indices[0]
    if zinc_index >= n_solute_atoms:
        raise RuntimeError(
            f"The only Zn atom in {source_pdb} is outside the first "
            f"{n_solute_atoms} solute atoms."
        )

    solvent_elements = elements[n_solute_atoms:]
    if len(solvent_elements) % 3:
        raise RuntimeError(
            f"The solvent portion of {source_pdb} has "
            f"{len(solvent_elements)} atoms, which is not divisible by 3."
        )
    water_oxygen_indices = []
    for offset in range(0, len(solvent_elements), 3):
        triplet = solvent_elements[offset:offset + 3]
        if sorted(triplet) != ["H", "H", "O"]:
            first = n_solute_atoms + offset + 1
            raise RuntimeError(
                f"Solvent atoms {first}-{first + 2} in {source_pdb} do not "
                f"form one H2O triplet (elements: {', '.join(triplet)})."
            )
        oxygen_local = triplet.index("O")
        water_oxygen_indices.append(n_solute_atoms + offset + oxygen_local)

    if expected_solute_pdb is not None:
        expected_atoms = _validated_pdb_atoms(
            expected_solute_pdb, "registered system solute"
        )
        if len(expected_atoms) != n_solute_atoms:
            raise RuntimeError(
                f"Registered solute {expected_solute_pdb} has "
                f"{len(expected_atoms)} atoms, but --co2-solute-atoms is "
                f"{n_solute_atoms}."
            )
        expected_elements = _element_sequence(expected_atoms)
        if elements[:n_solute_atoms] != expected_elements:
            mismatch = next(
                index for index, (observed, expected) in enumerate(
                    zip(elements[:n_solute_atoms], expected_elements), start=1
                )
                if observed != expected
            )
            raise RuntimeError(
                f"CO2 source solute element sequence differs from "
                f"{expected_solute_pdb} at atom {mismatch}: found "
                f"{elements[mismatch - 1]}, expected "
                f"{expected_elements[mismatch - 1]}."
            )

    return {
        "source_pdb": str(source_pdb.resolve()),
        "source_sha256": file_sha256(source_pdb),
        "n_source_atoms": len(atoms),
        "n_solute_atoms": n_solute_atoms,
        "n_water_atoms": len(solvent_elements),
        "n_waters": len(solvent_elements) // 3,
        "zinc_index": zinc_index + 1,
        "water_oxygen_indices": [index + 1 for index in water_oxygen_indices],
        "element_sequence": elements,
    }


def _angle_degrees(first, vertex, third) -> float:
    vector_a = tuple(a - b for a, b in zip(first, vertex))
    vector_b = tuple(c - b for c, b in zip(third, vertex))
    norm_a = math.sqrt(sum(value * value for value in vector_a))
    norm_b = math.sqrt(sum(value * value for value in vector_b))
    if norm_a == 0.0 or norm_b == 0.0:
        raise RuntimeError("Cannot compute a CO2 angle with a zero-length bond.")
    cosine = sum(a * b for a, b in zip(vector_a, vector_b)) / (
        norm_a * norm_b
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def validate_co2_template(co2_pdb: Path) -> dict:
    atoms = _validated_pdb_atoms(co2_pdb, "CO2 template")
    elements = _element_sequence(atoms)
    if len(atoms) != 3 or elements.count("C") != 1 or elements.count("O") != 2:
        raise RuntimeError(
            f"CO2 template {co2_pdb} must contain exactly three atoms: "
            "one C and two O atoms."
        )
    carbon_index = elements.index("C")
    oxygen_indices = [i for i, element in enumerate(elements) if element == "O"]
    carbon_xyz = atoms[carbon_index]["xyz"]
    distances = [
        math.dist(carbon_xyz, atoms[index]["xyz"])
        for index in oxygen_indices
    ]
    if any(distance <= 0.0 for distance in distances):
        raise RuntimeError(f"CO2 template {co2_pdb} has a zero C-O distance.")
    angle = _angle_degrees(
        atoms[oxygen_indices[0]]["xyz"],
        carbon_xyz,
        atoms[oxygen_indices[1]]["xyz"],
    )
    if angle < 170.0:
        raise RuntimeError(
            f"CO2 template {co2_pdb} is not sufficiently linear: "
            f"O-C-O = {angle:.3f} degrees (minimum 170 degrees)."
        )
    return {
        "co2_pdb": str(co2_pdb.resolve()),
        "co2_sha256": file_sha256(co2_pdb),
        "n_atoms": 3,
        "carbon_local_index": carbon_index + 1,
        "oxygen_local_indices": [index + 1 for index in oxygen_indices],
        "co_bond_lengths_A": distances,
        "oco_angle_degrees": angle,
        "element_sequence": elements,
    }


def co2_pack_seed(seed_base: int, co2_count: int, pack_index: int) -> int:
    """Return a stable, positive Packmol seed unique to count/replica pairs."""
    if seed_base <= 0 or co2_count <= 0 or pack_index <= 0:
        raise ValueError("CO2 seed base, count, and pack index must be positive.")
    pair_sum = co2_count + pack_index - 2
    pairing = pair_sum * (pair_sum + 1) // 2 + pack_index
    return ((seed_base - 1 + pairing) % 2_147_483_646) + 1


def co2_effective_target_distance(
    shell_inner_A: float,
    shell_outer_A: float,
    requested_distance_A: float | None,
) -> float:
    if requested_distance_A is None:
        return (shell_inner_A + shell_outer_A) / 2.0
    return requested_distance_A


def co2_site_direction_metadata(
    centered_pdb: Path,
    source_metadata: dict,
    direction_atom_index: int,
    target_distance_A: float,
    target_radius_A: float,
) -> dict:
    """Derive a site direction from 1-based atom order in centered geometry."""
    atoms = _validated_pdb_atoms(centered_pdb, "site-direction geometry")
    n_source_atoms = source_metadata["n_source_atoms"]
    if len(atoms) < n_source_atoms:
        raise RuntimeError(
            f"Site-direction geometry {centered_pdb} has {len(atoms)} atoms; "
            f"the source requires at least {n_source_atoms}."
        )
    if _element_sequence(atoms[:n_source_atoms]) != (
        source_metadata["element_sequence"]
    ):
        raise RuntimeError(
            "Site-direction geometry does not preserve the source atom order."
        )
    if direction_atom_index < 1:
        raise RuntimeError("CO2 direction atom index must be >= 1.")
    if direction_atom_index > n_source_atoms:
        raise RuntimeError(
            f"CO2 direction atom index {direction_atom_index} is outside "
            f"the source atom range 1-{n_source_atoms}."
        )
    zinc_index = source_metadata["zinc_index"]
    if direction_atom_index == zinc_index:
        raise RuntimeError(
            f"CO2 direction atom {direction_atom_index} is the Zn atom."
        )
    if not math.isfinite(target_distance_A) or target_distance_A <= 0.0:
        raise RuntimeError("CO2 target distance must be finite and > 0 A.")
    if not math.isfinite(target_radius_A) or target_radius_A <= 0.0:
        raise RuntimeError("CO2 target radius must be finite and > 0 A.")

    zinc_xyz = atoms[zinc_index - 1]["xyz"]
    direction_atom = atoms[direction_atom_index - 1]
    reference_xyz = direction_atom["xyz"]
    vector = tuple(
        reference - zinc
        for reference, zinc in zip(reference_xyz, zinc_xyz)
    )
    distance = math.sqrt(sum(value * value for value in vector))
    if distance <= 1.0e-6:
        raise RuntimeError(
            "CO2 direction atom is coincident with Zn within 1e-6 A."
        )
    unit_vector = tuple(value / distance for value in vector)
    target_xyz = tuple(
        zinc + target_distance_A * direction
        for zinc, direction in zip(zinc_xyz, unit_vector)
    )
    return {
        "zinc_index": zinc_index,
        "direction_atom_index": direction_atom_index,
        "direction_atom_element": direction_atom["element"],
        "zinc_xyz_A": list(zinc_xyz),
        "direction_atom_xyz_A": list(reference_xyz),
        "zn_direction_atom_distance_A": distance,
        "unit_vector": list(unit_vector),
        "target_xyz_A": list(target_xyz),
        "target_distance_A": target_distance_A,
        "target_radius_A": target_radius_A,
    }


def co2_targeted_carbon_global_atom_index(
    source_metadata: dict,
    template_metadata: dict,
) -> int:
    """Return the 1-based carbon index of targeted CO2 molecule 1."""
    return (
        source_metadata["n_source_atoms"]
        + template_metadata["carbon_local_index"]
    )


def co2_packmol_input(
    tolerance_A: float,
    seed: int,
    co2_count: int,
    carbon_local_index: int,
    zinc_xyz,
    shell_inner_A: float,
    shell_outer_A: float,
    *,
    placement_mode: str = "random-shell",
    site_direction: dict | None = None,
) -> str:
    zinc_x, zinc_y, zinc_z = zinc_xyz
    header = f"""tolerance {tolerance_A:.6f}
filetype pdb
output packed_CO2_shell.pdb
seed {seed}

structure source_centered.pdb
  number 1
  fixed 0. 0. 0. 0. 0. 0.
end structure
"""
    shell_constraints = f"""    outside sphere {zinc_x:.6f} {zinc_y:.6f} {zinc_z:.6f} {shell_inner_A:.6f}
    inside sphere {zinc_x:.6f} {zinc_y:.6f} {zinc_z:.6f} {shell_outer_A:.6f}"""
    if placement_mode == "random-shell":
        return header + f"""
structure co2.pdb
  number {co2_count}
  atoms {carbon_local_index}
{shell_constraints}
  end atoms
end structure
"""
    if placement_mode != "site-directed":
        raise ValueError(f"Unsupported CO2 placement mode: {placement_mode}.")
    if site_direction is None:
        raise ValueError("site-directed Packmol input requires site metadata.")
    target_x, target_y, target_z = site_direction["target_xyz_A"]
    target_radius_A = site_direction["target_radius_A"]
    targeted_block = f"""
structure co2.pdb
  number 1
  atoms {carbon_local_index}
{shell_constraints}
    inside sphere {target_x:.6f} {target_y:.6f} {target_z:.6f} {target_radius_A:.6f}
  end atoms
end structure
"""
    if co2_count == 1:
        return header + targeted_block
    background_block = f"""
structure co2.pdb
  number {co2_count - 1}
  atoms {carbon_local_index}
{shell_constraints}
  end atoms
end structure
"""
    return header + targeted_block + background_block


def _co2_carbon_indices(
    n_source_atoms: int,
    co2_count: int,
    carbon_local_index: int,
) -> list[int]:
    carbon_offset = carbon_local_index - 1
    return [
        n_source_atoms + molecule * 3 + carbon_offset
        for molecule in range(co2_count)
    ]


def co2_geometry_metrics(
    pdb_path: Path,
    source_metadata: dict,
    template_metadata: dict,
    co2_count: int,
) -> dict:
    atoms = _validated_pdb_atoms(pdb_path, "CO2-screen geometry")
    expected_atoms = source_metadata["n_source_atoms"] + 3 * co2_count
    if len(atoms) != expected_atoms:
        raise RuntimeError(
            f"{pdb_path} has {len(atoms)} atoms; expected {expected_atoms}."
        )
    zinc_index = source_metadata["zinc_index"] - 1
    carbon_indices = _co2_carbon_indices(
        source_metadata["n_source_atoms"],
        co2_count,
        template_metadata["carbon_local_index"],
    )
    zinc_xyz = atoms[zinc_index]["xyz"]
    zinc_carbon_distances = [
        math.dist(zinc_xyz, atoms[index]["xyz"])
        for index in carbon_indices
    ]
    water_oxygen_indices = [
        index - 1 for index in source_metadata["water_oxygen_indices"]
    ]
    nearest_water = None
    for oxygen_index in water_oxygen_indices:
        for molecule_index, carbon_index in enumerate(carbon_indices, start=1):
            distance = math.dist(
                atoms[oxygen_index]["xyz"], atoms[carbon_index]["xyz"]
            )
            candidate = (distance, oxygen_index + 1, molecule_index, carbon_index + 1)
            if nearest_water is None or candidate < nearest_water:
                nearest_water = candidate

    return {
        "zn_c_distances_A": zinc_carbon_distances,
        "zn_c_min_A": min(zinc_carbon_distances),
        "zn_c_mean_A": sum(zinc_carbon_distances) / len(zinc_carbon_distances),
        "zn_c_max_A": max(zinc_carbon_distances),
        "nearest_zn_c_co2_index": (
            zinc_carbon_distances.index(min(zinc_carbon_distances)) + 1
        ),
        "water_o_c_min_A": nearest_water[0] if nearest_water else None,
        "nearest_water_oxygen_atom_index": (
            nearest_water[1] if nearest_water else None
        ),
        "nearest_water_o_c_co2_index": (
            nearest_water[2] if nearest_water else None
        ),
        "nearest_water_o_c_carbon_atom_index": (
            nearest_water[3] if nearest_water else None
        ),
    }


def validate_co2_site_direction_geometry(
    pdb_path: Path,
    source_metadata: dict,
    template_metadata: dict,
    site_direction: dict,
    coordinate_tolerance_A: float,
) -> dict:
    """Recompute and validate targeted-CO2 geometry in the current frame."""
    current_direction = co2_site_direction_metadata(
        pdb_path,
        source_metadata,
        site_direction["direction_atom_index"],
        site_direction["target_distance_A"],
        site_direction["target_radius_A"],
    )
    atoms = _validated_pdb_atoms(pdb_path, "site-directed CO2 geometry")
    carbon_atom_index = co2_targeted_carbon_global_atom_index(
        source_metadata, template_metadata
    )
    carbon_xyz = atoms[carbon_atom_index - 1]["xyz"]
    zinc_xyz = atoms[current_direction["zinc_index"] - 1]["xyz"]
    direction_xyz = atoms[
        current_direction["direction_atom_index"] - 1
    ]["xyz"]
    zinc_carbon_distance = math.dist(zinc_xyz, carbon_xyz)
    target_point_distance = math.dist(
        carbon_xyz, current_direction["target_xyz_A"]
    )
    if (
        target_point_distance
        > current_direction["target_radius_A"] + coordinate_tolerance_A
    ):
        raise RuntimeError(
            "Targeted CO2 carbon is outside the site-directed target "
            f"sphere: C-target = {target_point_distance:.6f} A, radius = "
            f"{current_direction['target_radius_A']:.6f} A, tolerance = "
            f"{coordinate_tolerance_A:.6f} A."
        )
    angle = _angle_degrees(direction_xyz, zinc_xyz, carbon_xyz)
    return {
        "targeted_CO2_index": 1,
        "targeted_carbon_atom_index": carbon_atom_index,
        "targeted_Zn_C_distance_A": zinc_carbon_distance,
        "targeted_target_point_distance_A": target_point_distance,
        "target_radius_A": current_direction["target_radius_A"],
        "targeted_direction_angle_degrees": angle,
        "site_direction_valid": True,
        "recomputed_site_direction": current_direction,
    }


def validate_co2_packed_geometry(
    source_centered_pdb: Path,
    packed_pdb: Path,
    source_metadata: dict,
    template_metadata: dict,
    co2_count: int,
    shell_inner_A: float,
    shell_outer_A: float,
    coordinate_tolerance_A: float = PDB_COORDINATE_TOLERANCE_A,
    site_direction: dict | None = None,
) -> dict:
    source_atoms = _validated_pdb_atoms(source_centered_pdb, "centered source")
    packed_atoms = _validated_pdb_atoms(packed_pdb, "Packmol CO2 output")
    expected_count = len(source_atoms) + 3 * co2_count
    if len(packed_atoms) != expected_count:
        raise RuntimeError(
            f"Packmol CO2 output has {len(packed_atoms)} atoms; expected "
            f"{expected_count}."
        )
    expected_elements = _element_sequence(source_atoms) + (
        template_metadata["element_sequence"] * co2_count
    )
    observed_elements = _element_sequence(packed_atoms)
    if observed_elements != expected_elements:
        mismatch = next(
            index for index, (observed, expected) in enumerate(
                zip(observed_elements, expected_elements), start=1
            )
            if observed != expected
        )
        raise RuntimeError(
            f"Packmol CO2 atom ordering/composition differs at atom "
            f"{mismatch}: found {observed_elements[mismatch - 1]}, expected "
            f"{expected_elements[mismatch - 1]}."
        )

    source_displacements = [
        math.dist(source["xyz"], packed["xyz"])
        for source, packed in zip(source_atoms, packed_atoms)
    ]
    max_source_displacement = max(source_displacements)
    if max_source_displacement > coordinate_tolerance_A:
        raise RuntimeError(
            "Packmol changed the fixed source geometry by "
            f"{max_source_displacement:.6f} A (allowed "
            f"{coordinate_tolerance_A:.6f} A)."
        )

    metrics = co2_geometry_metrics(
        packed_pdb, source_metadata, template_metadata, co2_count
    )
    outside = [
        (index, distance) for index, distance in enumerate(
            metrics["zn_c_distances_A"], start=1
        )
        if distance < shell_inner_A - coordinate_tolerance_A
        or distance > shell_outer_A + coordinate_tolerance_A
    ]
    if outside:
        detail = ", ".join(
            f"CO2 {index}: {distance:.6f} A" for index, distance in outside
        )
        raise RuntimeError(
            f"Packmol placed CO2 carbon(s) outside the requested Zn shell "
            f"[{shell_inner_A:.6f}, {shell_outer_A:.6f}] A: {detail}."
        )

    n_source = source_metadata["n_source_atoms"]
    co2_atoms = packed_atoms[n_source:]
    min_co2_source = min(
        math.dist(co2_atom["xyz"], source_atom["xyz"])
        for co2_atom in co2_atoms
        for source_atom in source_atoms
    )
    min_co2_co2 = None
    if co2_count > 1:
        min_co2_co2 = min(
            math.dist(
                co2_atoms[first_molecule * 3 + first_atom]["xyz"],
                co2_atoms[second_molecule * 3 + second_atom]["xyz"],
            )
            for first_molecule in range(co2_count)
            for second_molecule in range(first_molecule + 1, co2_count)
            for first_atom in range(3)
            for second_atom in range(3)
        )

    validation = {
        "atom_count_valid": True,
        "element_order_valid": True,
        "source_coordinates_preserved": True,
        "source_max_displacement_A": max_source_displacement,
        "coordinate_tolerance_A": coordinate_tolerance_A,
        "shell_valid": True,
        "minimum_CO2_source_distance_A": min_co2_source,
        "minimum_CO2_CO2_distance_A": min_co2_co2,
        "metrics": metrics,
    }
    if site_direction is not None:
        validation.update(validate_co2_site_direction_geometry(
            packed_pdb,
            source_metadata,
            template_metadata,
            site_direction,
            coordinate_tolerance_A,
        ))
    return validation


def validate_co2_centered_geometry(
    packed_pdb: Path,
    centered_pdb: Path,
    source_metadata: dict,
    template_metadata: dict,
    co2_count: int,
    shell_inner_A: float | None = None,
    shell_outer_A: float | None = None,
    site_direction: dict | None = None,
) -> dict:
    packed_atoms = _validated_pdb_atoms(packed_pdb, "packed CO2 geometry")
    centered_atoms = _validated_pdb_atoms(centered_pdb, "centered CO2 geometry")
    if _element_sequence(packed_atoms) != _element_sequence(centered_atoms):
        raise RuntimeError("Final centering changed the CO2-system element order.")
    before = co2_geometry_metrics(
        packed_pdb, source_metadata, template_metadata, co2_count
    )
    after = co2_geometry_metrics(
        centered_pdb, source_metadata, template_metadata, co2_count
    )
    maximum_change = max(
        abs(first - second) for first, second in zip(
            before["zn_c_distances_A"], after["zn_c_distances_A"]
        )
    )
    if maximum_change > 2 * PDB_COORDINATE_TOLERANCE_A:
        raise RuntimeError(
            "Centering changed a Zn-C distance by "
            f"{maximum_change:.6f} A, beyond PDB precision tolerance."
        )
    if (shell_inner_A is None) != (shell_outer_A is None):
        raise ValueError("Both centered-shell bounds must be provided together.")
    if shell_inner_A is not None:
        shell_tolerance = 2 * PDB_COORDINATE_TOLERANCE_A
        outside = [
            (index, distance) for index, distance in enumerate(
                after["zn_c_distances_A"], start=1
            )
            if distance < shell_inner_A - shell_tolerance
            or distance > shell_outer_A + shell_tolerance
        ]
        if outside:
            raise RuntimeError(
                "Final centering did not preserve the validated Zn-C shell: "
                + ", ".join(
                    f"CO2 {index}: {distance:.6f} A"
                    for index, distance in outside
                )
            )
    validation = {
        "element_order_valid": True,
        "maximum_zn_c_distance_change_A": maximum_change,
        "shell_revalidated": shell_inner_A is not None,
        "metrics": after,
    }
    if site_direction is not None:
        validation.update(validate_co2_site_direction_geometry(
            centered_pdb,
            source_metadata,
            template_metadata,
            site_direction,
            2 * PDB_COORDINATE_TOLERANCE_A,
        ))
    return validation


def _read_json_dict(path: Path, description: str) -> dict:
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Cannot read {description} {path}.") from exc
    if not isinstance(data, dict):
        raise RuntimeError(f"Invalid {description} {path}: expected an object.")
    return data


def _write_json(path: Path, data: dict):
    path.write_text(json.dumps(data, indent=2) + "\n")


def co2_condition_dir(
    co2_project: Path,
    system_name: str,
    co2_count: int,
    pack_index: int,
) -> Path:
    return (
        co2_project
        / system_name
        / f"NCO2_{co2_count:02d}"
        / f"pack_{pack_index:02d}"
    )


def _co2_stage_active_entries(stage_dir: Path) -> list[Path]:
    if not stage_dir.exists():
        return []
    return sorted(
        (path for path in stage_dir.iterdir() if path.name != "attempts"),
        key=lambda path: path.name,
    )


def archive_co2_stage_attempt(stage_dir: Path, reason: str) -> Path | None:
    """Move a prior attempt aside without deleting scientific outputs."""
    entries = _co2_stage_active_entries(stage_dir)
    if not entries:
        return None
    attempts_dir = stage_dir / "attempts"
    attempts_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    archive = attempts_dir / stamp
    suffix = 1
    while archive.exists():
        archive = attempts_dir / f"{stamp}_{suffix:02d}"
        suffix += 1
    archive.mkdir()
    for path in entries:
        shutil.move(str(path), archive / path.name)
    (archive / "archive_reason.txt").write_text(reason.rstrip() + "\n")
    return archive


def _validate_co2_output_hashes(stage_dir: Path, manifest: dict):
    output_hashes = manifest.get("output_sha256")
    if not isinstance(output_hashes, dict) or not output_hashes:
        raise RuntimeError(
            f"CO2 stage manifest in {stage_dir} has no output SHA-256 map."
        )
    for filename, recorded_hash in output_hashes.items():
        path = stage_dir / filename
        if not path.is_file():
            raise RuntimeError(
                f"Completed CO2 stage {stage_dir} is missing {filename}."
            )
        if file_sha256(path) != recorded_hash:
            raise RuntimeError(
                f"Completed CO2 stage output changed after validation: {path}."
            )


def co2_pack_configuration(
    system_name: str,
    source_metadata: dict,
    template_metadata: dict,
    co2_count: int,
    pack_index: int,
    seed: int,
    args,
) -> dict:
    placement_mode = getattr(args, "co2_placement_mode", "random-shell")
    configuration = {
        "stage": CO2_PACK_STAGE,
        "system": system_name,
        "source_sha256": source_metadata["source_sha256"],
        "source_atom_count": source_metadata["n_source_atoms"],
        "solute_atom_count": source_metadata["n_solute_atoms"],
        "water_count": source_metadata["n_waters"],
        "co2_template_sha256": template_metadata["co2_sha256"],
        "co2_carbon_local_index": template_metadata["carbon_local_index"],
        "co2_count": co2_count,
        "pack_index": pack_index,
        "packmol_seed": seed,
        "packmol_tolerance_A": args.packmol_tolerance,
        "shell_inner_A": args.co2_shell_inner,
        "shell_outer_A": args.co2_shell_outer,
        "wall_margin_A": args.wall_margin,
        "placement_mode": placement_mode,
    }
    if placement_mode == "site-directed":
        direction_atom_index = getattr(args, "co2_direction_atom", None)
        if direction_atom_index is None:
            raise RuntimeError(
                "site-directed CO2 placement requires --co2-direction-atom."
            )
        if direction_atom_index < 1:
            raise RuntimeError("CO2 direction atom index must be >= 1.")
        if direction_atom_index > source_metadata["n_source_atoms"]:
            raise RuntimeError(
                f"CO2 direction atom index {direction_atom_index} is "
                "outside the source atom range 1-"
                f"{source_metadata['n_source_atoms']}."
            )
        if direction_atom_index == source_metadata["zinc_index"]:
            raise RuntimeError(
                f"CO2 direction atom {direction_atom_index} is the Zn atom."
            )
        target_distance_A = co2_effective_target_distance(
            args.co2_shell_inner,
            args.co2_shell_outer,
            getattr(args, "co2_target_distance", None),
        )
        target_radius_A = getattr(args, "co2_target_radius", 1.5)
        if (
            not math.isfinite(target_distance_A)
            or not (
                args.co2_shell_inner
                <= target_distance_A
                <= args.co2_shell_outer
            )
        ):
            raise RuntimeError(
                "CO2 target distance must be finite and inside the Zn shell."
            )
        if not math.isfinite(target_radius_A) or target_radius_A <= 0.0:
            raise RuntimeError("CO2 target radius must be finite and > 0 A.")
        configuration.update({
            "direction_atom_index": direction_atom_index,
            "direction_atom_element": source_metadata[
                "element_sequence"
            ][direction_atom_index - 1],
            "zinc_index": source_metadata["zinc_index"],
            "target_distance_A": target_distance_A,
            "target_radius_A": target_radius_A,
        })
    return configuration


def co2_site_direction_from_configuration(
    centered_pdb: Path,
    source_metadata: dict,
    template_metadata: dict,
    configuration: dict,
) -> dict | None:
    if configuration.get("placement_mode") != "site-directed":
        return None
    metadata = co2_site_direction_metadata(
        centered_pdb,
        source_metadata,
        configuration["direction_atom_index"],
        configuration["target_distance_A"],
        configuration["target_radius_A"],
    )
    metadata.update({
        "targeted_CO2_index": 1,
        "targeted_carbon_global_atom_index": (
            co2_targeted_carbon_global_atom_index(
                source_metadata, template_metadata
            )
        ),
    })
    return metadata


def validate_completed_co2_pack(
    stage_dir: Path,
    expected_configuration: dict,
    source_metadata: dict,
    template_metadata: dict,
) -> dict:
    done = stage_dir / "stage.done"
    if not done.is_file():
        raise RuntimeError(f"Missing completed-stage marker {done}.")
    conflicts = [
        marker for marker in ("stage.running", "stage.failed")
        if (stage_dir / marker).exists()
    ]
    if conflicts:
        raise RuntimeError(
            f"Completed {CO2_PACK_STAGE} has conflicting markers: "
            f"{', '.join(conflicts)}. Use --co2-repack after inspection."
        )
    manifest = _read_json_dict(
        stage_dir / "stage_manifest.json", "CO2 packing manifest"
    )
    if manifest.get("status") != "completed":
        raise RuntimeError(
            f"{CO2_PACK_STAGE} stage.done conflicts with manifest status "
            f"{manifest.get('status')!r}."
        )
    recorded_configuration = manifest.get("configuration")
    if not isinstance(recorded_configuration, dict):
        raise RuntimeError(f"{CO2_PACK_STAGE} manifest has no configuration.")
    comparison_configuration = dict(recorded_configuration)
    if (
        expected_configuration["placement_mode"] == "random-shell"
        and "placement_mode" not in comparison_configuration
    ):
        # Packing manifests predating the placement option can only have
        # been generated by the historical random-shell implementation.
        comparison_configuration["placement_mode"] = "random-shell"
    mismatches = configuration_mismatches(
        expected_configuration, comparison_configuration
    )
    if mismatches:
        raise RuntimeError(
            f"Cannot reuse {CO2_PACK_STAGE}: provenance is incompatible "
            f"(fields: {', '.join(mismatches)}). Use --co2-repack to "
            "archive and regenerate this packing."
        )
    _validate_co2_output_hashes(stage_dir, manifest)
    if file_sha256(stage_dir / "source_medoid.pdb") != (
        source_metadata["source_sha256"]
    ):
        raise RuntimeError(
            f"Cannot reuse {CO2_PACK_STAGE}: source_medoid.pdb does not "
            "match the requested source."
        )
    if file_sha256(stage_dir / "co2.pdb") != template_metadata["co2_sha256"]:
        raise RuntimeError(
            f"Cannot reuse {CO2_PACK_STAGE}: co2.pdb does not match the "
            "requested template."
        )
    site_direction = co2_site_direction_from_configuration(
        stage_dir / "source_centered.pdb",
        source_metadata,
        template_metadata,
        expected_configuration,
    )
    if site_direction is not None:
        if manifest.get("site_direction") != site_direction:
            raise RuntimeError(
                f"Cannot reuse {CO2_PACK_STAGE}: recorded site_direction "
                "metadata differs from the centered source geometry."
            )
        if manifest.get("site_directed_scientific_note") != (
            CO2_SITE_DIRECTED_NOTE
        ):
            raise RuntimeError(
                f"Cannot reuse {CO2_PACK_STAGE}: missing or incompatible "
                "site-directed scientific note."
            )
    centered_source_atoms = _validated_pdb_atoms(
        stage_dir / "source_centered.pdb", "centered CO2 source"
    )
    zinc_xyz = centered_source_atoms[
        source_metadata["zinc_index"] - 1
    ]["xyz"]
    expected_packmol_input = co2_packmol_input(
        expected_configuration["packmol_tolerance_A"],
        expected_configuration["packmol_seed"],
        expected_configuration["co2_count"],
        expected_configuration["co2_carbon_local_index"],
        zinc_xyz,
        expected_configuration["shell_inner_A"],
        expected_configuration["shell_outer_A"],
        placement_mode=expected_configuration["placement_mode"],
        site_direction=site_direction,
    )
    if (stage_dir / "06_CO2_shell_pack.inp").read_text() != (
        expected_packmol_input
    ):
        raise RuntimeError(
            f"Cannot reuse {CO2_PACK_STAGE}: archived Packmol input is "
            "incompatible with the requested placement."
        )
    validation = validate_co2_packed_geometry(
        stage_dir / "source_centered.pdb",
        stage_dir / "packed_CO2_shell.pdb",
        source_metadata,
        template_metadata,
        expected_configuration["co2_count"],
        expected_configuration["shell_inner_A"],
        expected_configuration["shell_outer_A"],
        site_direction=site_direction,
    )
    centered_validation = validate_co2_centered_geometry(
        stage_dir / "packed_CO2_shell.pdb",
        stage_dir / "system_CO2_centered.pdb",
        source_metadata,
        template_metadata,
        expected_configuration["co2_count"],
        expected_configuration["shell_inner_A"],
        expected_configuration["shell_outer_A"],
        site_direction=site_direction,
    )
    if site_direction is not None:
        recorded_validation = manifest.get("packing_validation", {})
        for field in (
            "targeted_CO2_index",
            "targeted_carbon_atom_index",
            "targeted_Zn_C_distance_A",
            "targeted_target_point_distance_A",
            "target_radius_A",
            "targeted_direction_angle_degrees",
            "site_direction_valid",
        ):
            if recorded_validation.get(field) != validation.get(field):
                raise RuntimeError(
                    f"Cannot reuse {CO2_PACK_STAGE}: recorded site-directed "
                    f"validation field {field} is inconsistent."
                )
    recorded_wall = manifest.get("wall", {})
    centered_atoms = _validated_pdb_atoms(
        stage_dir / "system_CO2_centered.pdb", "centered CO2 system"
    )
    maximum_radius = max(
        math.sqrt(sum(value * value for value in atom["xyz"]))
        for atom in centered_atoms
    )
    wall_radius_A = maximum_radius + expected_configuration["wall_margin_A"]
    if not math.isclose(
        recorded_wall.get("radius_A", math.nan),
        wall_radius_A,
        abs_tol=PDB_COORDINATE_TOLERANCE_A,
    ):
        raise RuntimeError(
            f"Cannot reuse {CO2_PACK_STAGE}: recorded wall radius is "
            "incompatible with the centered packed geometry."
        )
    if not math.isclose(
        recorded_wall.get("radius_bohr", math.nan),
        recorded_wall.get("radius_A", math.nan) * BOHR_PER_ANGSTROM,
        abs_tol=PDB_COORDINATE_TOLERANCE_A * BOHR_PER_ANGSTROM,
    ):
        raise RuntimeError(
            f"Cannot reuse {CO2_PACK_STAGE}: recorded wall radii in A and "
            "bohr are inconsistent."
        )
    manifest["reuse_validation"] = {
        "packing": validation,
        "centering": centered_validation,
    }
    return manifest


def run_co2_pack_stage(
    condition_dir: Path,
    system_name: str,
    source_pdb: Path,
    co2_pdb: Path,
    source_metadata: dict,
    template_metadata: dict,
    co2_count: int,
    pack_index: int,
    args,
) -> dict:
    stage_dir = condition_dir / CO2_PACK_STAGE
    seed = co2_pack_seed(args.co2_seed_base, co2_count, pack_index)
    configuration = co2_pack_configuration(
        system_name,
        source_metadata,
        template_metadata,
        co2_count,
        pack_index,
        seed,
        args,
    )

    if (stage_dir / "stage.done").is_file() and not args.co2_repack:
        manifest = validate_completed_co2_pack(
            stage_dir, configuration, source_metadata, template_metadata
        )
        print(
            f"  REUSE {system_name}/NCO2_{co2_count:02d}/"
            f"pack_{pack_index:02d}/{CO2_PACK_STAGE}"
        )
        return manifest
    active_entries = _co2_stage_active_entries(stage_dir)
    if active_entries and not args.co2_repack:
        names = ", ".join(path.name for path in active_entries[:6])
        raise RuntimeError(
            f"Existing incomplete or failed {CO2_PACK_STAGE} attempt in "
            f"{stage_dir} ({names}). Inspect it and use --co2-repack to "
            "archive it before a new packing."
        )

    packmol_exe = shutil.which(args.packmol)
    if packmol_exe is None:
        raise RuntimeError(
            f"Packmol executable '{args.packmol}' not found in PATH. "
            "Use --packmol /path/to/packmol if needed."
        )

    if active_entries:
        archive = archive_co2_stage_attempt(
            stage_dir, "Replaced explicitly with --co2-repack."
        )
        print(f"  ARCHIVE prior {CO2_PACK_STAGE} attempt at {archive}")
    stage_dir.mkdir(parents=True, exist_ok=True)
    source_local = stage_dir / "source_medoid.pdb"
    source_centered = stage_dir / "source_centered.pdb"
    co2_local = stage_dir / "co2.pdb"
    packed_pdb = stage_dir / "packed_CO2_shell.pdb"
    final_centered = stage_dir / "system_CO2_centered.pdb"
    input_path = stage_dir / "06_CO2_shell_pack.inp"
    log_path = stage_dir / "06_CO2_shell_pack.out"

    shutil.copy2(source_pdb, source_local)
    shutil.copy2(co2_pdb, co2_local)
    source_centering = center_pdb(source_local, source_centered)
    centered_source_atoms = _validated_pdb_atoms(
        source_centered, "centered CO2 source"
    )
    zinc_xyz = centered_source_atoms[source_metadata["zinc_index"] - 1]["xyz"]
    site_direction = co2_site_direction_from_configuration(
        source_centered,
        source_metadata,
        template_metadata,
        configuration,
    )
    input_path.write_text(co2_packmol_input(
        args.packmol_tolerance,
        seed,
        co2_count,
        template_metadata["carbon_local_index"],
        zinc_xyz,
        args.co2_shell_inner,
        args.co2_shell_outer,
        placement_mode=configuration["placement_mode"],
        site_direction=site_direction,
    ))
    command = [args.packmol]
    started_at = datetime.now(timezone.utc).isoformat()
    runtime = {
        "workflow": "CO2_shell_screening",
        "stage": CO2_PACK_STAGE,
        "status": "running",
        "started_at": started_at,
        "hostname": socket.gethostname(),
        "command": command,
        "configuration": configuration,
        "execution_resources": co2_execution_resources(args),
    }
    if site_direction is not None:
        runtime["site_direction"] = site_direction
    mark_stage_running(stage_dir, runtime)
    print(
        f"  PACK {system_name}/NCO2_{co2_count:02d}/pack_{pack_index:02d} "
        f"(seed {seed})"
    )
    print(f"       CO2 placement mode       : {configuration['placement_mode']}")
    if site_direction is not None:
        unit_vector = " ".join(
            f"{value:.6f}" for value in site_direction["unit_vector"]
        )
        target_center = " ".join(
            f"{value:.6f}" for value in site_direction["target_xyz_A"]
        )
        zinc_center = " ".join(
            f"{value:.6f}" for value in site_direction["zinc_xyz_A"]
        )
        reference_center = " ".join(
            f"{value:.6f}"
            for value in site_direction["direction_atom_xyz_A"]
        )
        print(f"       Zn atom index            : {site_direction['zinc_index']}")
        print(f"       Zn xyz                   : {zinc_center}")
        print(
            "       Direction atom index     : "
            f"{site_direction['direction_atom_index']}"
        )
        print(
            "       Direction atom element   : "
            f"{site_direction['direction_atom_element']}"
        )
        print(f"       Direction atom xyz       : {reference_center}")
        print(
            "       Zn-direction distance    : "
            f"{site_direction['zn_direction_atom_distance_A']:.6f} A"
        )
        print(f"       Direction unit vector    : {unit_vector}")
        print(
            "       Target distance from Zn  : "
            f"{site_direction['target_distance_A']:.6f} A"
        )
        print(
            f"       Target radius            : "
            f"{site_direction['target_radius_A']:.6f} A"
        )
        print(f"       Target center            : {target_center}")
        print("       Targeted CO2             : molecule 1")
        print("       Targeted CO2 molecules   : 1")
        print(
            f"       Background CO2 molecules : {co2_count - 1}"
        )

    returncode = None
    try:
        with input_path.open("r") as packmol_input_handle, log_path.open("w") as log:
            result = subprocess.run(
                [packmol_exe],
                cwd=stage_dir,
                stdin=packmol_input_handle,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )
        returncode = result.returncode
        if returncode != 0:
            raise RuntimeError(f"Packmol returned code {returncode}.")
        if not packed_pdb.is_file():
            raise RuntimeError("Packmol did not create packed_CO2_shell.pdb.")
        if "Success!" not in log_path.read_text(errors="replace"):
            print(
                f"  WARNING: {log_path} does not contain Packmol's usual "
                "'Success!' marker; geometric validation will decide reuse."
            )

        packing_validation = validate_co2_packed_geometry(
            source_centered,
            packed_pdb,
            source_metadata,
            template_metadata,
            co2_count,
            args.co2_shell_inner,
            args.co2_shell_outer,
            site_direction=site_direction,
        )
        final_centering = center_pdb(packed_pdb, final_centered)
        centered_validation = validate_co2_centered_geometry(
            packed_pdb,
            final_centered,
            source_metadata,
            template_metadata,
            co2_count,
            args.co2_shell_inner,
            args.co2_shell_outer,
            site_direction=site_direction,
        )
        wall_radius_A = final_centering["max_radius_from_COM_A"] + args.wall_margin
        wall_radius_bohr = wall_radius_A * BOHR_PER_ANGSTROM
        final_atoms = _validated_pdb_atoms(final_centered, "centered CO2 system")
        zinc_after_centering = final_atoms[
            source_metadata["zinc_index"] - 1
        ]["xyz"]
        finished_at = datetime.now(timezone.utc).isoformat()
        output_names = [
            source_local.name,
            source_centered.name,
            co2_local.name,
            input_path.name,
            log_path.name,
            packed_pdb.name,
            final_centered.name,
        ]
        manifest = {
            **runtime,
            "status": "completed",
            "finished_at": finished_at,
            "returncode": returncode,
            "source": source_metadata,
            "co2_template": template_metadata,
            "source_centering": source_centering,
            "zinc_center_in_packmol_frame_A": list(zinc_xyz),
            "zinc_coordinate_after_final_centering_A": list(
                zinc_after_centering
            ),
            "packing_validation": packing_validation,
            "final_centering": final_centering,
            "centered_validation": centered_validation,
            "wall": {
                "margin_A": args.wall_margin,
                "radius_A": wall_radius_A,
                "radius_bohr": wall_radius_bohr,
                "center_A": [0.0, 0.0, 0.0],
            },
            "composition": {
                "n_solute_atoms": source_metadata["n_solute_atoms"],
                "n_waters": source_metadata["n_waters"],
                "n_co2": co2_count,
                "n_total_atoms": source_metadata["n_source_atoms"] + 3 * co2_count,
            },
            "paths": {
                "source_pdb": source_local.name,
                "source_centered_pdb": source_centered.name,
                "co2_template_pdb": co2_local.name,
                "packmol_input": input_path.name,
                "packmol_log": log_path.name,
                "packed_pdb": packed_pdb.name,
                "centered_pdb": final_centered.name,
            },
            "composition_changed": True,
            "velocities_preserved": False,
            "restart_chain_continued": False,
            "mdrestart_used": False,
            "source_structure_type": "full-droplet representative frame",
            "scientific_note": CO2_WORKFLOW_NOTE,
            "output_sha256": {
                name: file_sha256(stage_dir / name) for name in output_names
            },
        }
        if site_direction is not None:
            manifest["site_direction"] = site_direction
            manifest["site_directed_scientific_note"] = (
                CO2_SITE_DIRECTED_NOTE
            )
        _write_json(stage_dir / "stage_manifest.json", manifest)
        mark_stage_done(stage_dir)
        if site_direction is not None:
            print(
                "       Targeted Zn-C            : "
                f"{packing_validation['targeted_Zn_C_distance_A']:.6f} A"
            )
            print(
                "       C-target distance        : "
                f"{packing_validation['targeted_target_point_distance_A']:.6f} A"
            )
            print(
                "       Zn-direction/C angle     : "
                f"{packing_validation['targeted_direction_angle_degrees']:.6f} deg"
            )
        print(f"  OK   {stage_dir}")
        return manifest
    except (OSError, RuntimeError, ValueError) as exc:
        failure = {
            **runtime,
            "status": "failed",
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "returncode": returncode,
            "failure_reason": str(exc),
            "source": source_metadata,
            "co2_template": template_metadata,
            "composition_changed": True,
            "velocities_preserved": False,
            "restart_chain_continued": False,
            "mdrestart_used": False,
            "source_structure_type": "full-droplet representative frame",
            "scientific_note": CO2_WORKFLOW_NOTE,
        }
        if site_direction is not None:
            failure["site_direction"] = site_direction
            failure["site_directed_scientific_note"] = (
                CO2_SITE_DIRECTED_NOTE
            )
        _write_json(stage_dir / "stage_manifest.json", failure)
        mark_stage_failed(stage_dir, str(exc))
        raise RuntimeError(
            f"{CO2_PACK_STAGE} failed in {stage_dir}: {exc} See "
            f"{log_path}. Consider reducing --co2-counts or widening "
            "--co2-shell-inner/--co2-shell-outer; the workflow will not "
            "remove water or relax the packing constraints silently."
        ) from exc


def validate_co2_accommodated_geometry(
    initial_pdb: Path,
    accommodated_pdb: Path,
    source_metadata: dict,
    template_metadata: dict,
    co2_count: int,
    fixed_tolerance_A: float = 1.0e-3,
) -> dict:
    initial_atoms = _validated_pdb_atoms(initial_pdb, "CO2 accommodation input")
    final_atoms = _validated_pdb_atoms(
        accommodated_pdb, "CO2 accommodation output"
    )
    expected_atoms = source_metadata["n_source_atoms"] + 3 * co2_count
    if len(initial_atoms) != expected_atoms or len(final_atoms) != expected_atoms:
        raise RuntimeError(
            "CO2 accommodation geometry has an invalid atom count: "
            f"input {len(initial_atoms)}, output {len(final_atoms)}, "
            f"expected {expected_atoms}."
        )
    if _element_sequence(initial_atoms) != _element_sequence(final_atoms):
        raise RuntimeError(
            "CO2 accommodation changed the input element sequence."
        )
    displacement = solute_displacement(
        initial_pdb, accommodated_pdb, source_metadata["n_solute_atoms"]
    )
    if displacement["solute_max_displacement_A"] > fixed_tolerance_A:
        raise RuntimeError(
            "CO2 accommodation changed fixed solute coordinates by "
            f"{displacement['solute_max_displacement_A']:.6f} A "
            f"(allowed {fixed_tolerance_A:.6f} A)."
        )
    return {
        "atom_count_valid": True,
        "element_order_valid": True,
        "fixed_solute_valid": True,
        "fixed_solute_tolerance_A": fixed_tolerance_A,
        **displacement,
        "metrics_before": co2_geometry_metrics(
            initial_pdb, source_metadata, template_metadata, co2_count
        ),
        "metrics_after": co2_geometry_metrics(
            accommodated_pdb, source_metadata, template_metadata, co2_count
        ),
    }


def co2_accommodation_configuration(
    packing_manifest_path: Path,
    packing_manifest: dict,
    args,
) -> dict:
    packing_configuration = packing_manifest["configuration"]
    return {
        "stage": CO2_ACCOMMODATION_STAGE,
        "system": packing_configuration["system"],
        "n_CO2": packing_configuration["co2_count"],
        "pack_index": packing_configuration["pack_index"],
        "input_geometry_sha256": file_sha256(
            packing_manifest_path.parent / "system_CO2_centered.pdb"
        ),
        "packing_manifest_sha256": file_sha256(packing_manifest_path),
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "threads": args.threads,
        "optimization_level": args.co2_accommodation_level,
        "optimization_engine": args.co2_accommodation_engine,
        "max_cycles": args.co2_accommodation_cycles,
        "wall_radius_A": packing_manifest["wall"]["radius_A"],
        "wall_radius_bohr": packing_manifest["wall"]["radius_bohr"],
        "fixed_atoms": f"1-{packing_manifest['source']['n_solute_atoms']}",
    }


def co2_accommodation_command(args) -> list[str]:
    command = [
        args.xtb,
        "system_CO2_centered.pdb",
        "--gfn", str(args.gfn),
        "--chrg", str(args.charge),
        "--uhf", str(args.uhf),
        "--opt", args.co2_accommodation_level,
        "--cycles", str(args.co2_accommodation_cycles),
        "--input", "07_CO2_accommodation.inp",
    ]
    if args.alpb:
        command += ["--alpb", args.alpb]
    return command


def _co2_stage_output_hashes(stage_dir: Path) -> dict:
    excluded = {
        "stage.running", "stage.done", "stage.failed", "stage_manifest.json"
    }
    return {
        path.name: file_sha256(path)
        for path in sorted(stage_dir.iterdir(), key=lambda item: item.name)
        if path.is_file() and path.name not in excluded
    }


def validate_completed_co2_accommodation(
    stage_dir: Path,
    expected_configuration: dict,
    source_metadata: dict,
    template_metadata: dict,
) -> dict:
    if not (stage_dir / "stage.done").is_file():
        raise RuntimeError(
            f"Missing completed-stage marker {stage_dir / 'stage.done'}."
        )
    conflicts = [
        marker for marker in ("stage.running", "stage.failed")
        if (stage_dir / marker).exists()
    ]
    if conflicts:
        raise RuntimeError(
            f"Completed {CO2_ACCOMMODATION_STAGE} has conflicting markers: "
            f"{', '.join(conflicts)}. Use --force after inspection."
        )
    manifest = _read_json_dict(
        stage_dir / "stage_manifest.json", "CO2 accommodation manifest"
    )
    if manifest.get("status") != "completed":
        raise RuntimeError(
            f"{CO2_ACCOMMODATION_STAGE} stage.done conflicts with manifest "
            f"status {manifest.get('status')!r}."
        )
    recorded_configuration = manifest.get("configuration")
    if not isinstance(recorded_configuration, dict):
        raise RuntimeError(
            f"{CO2_ACCOMMODATION_STAGE} manifest has no configuration."
        )
    mismatches = configuration_mismatches(
        expected_configuration, recorded_configuration
    )
    if mismatches:
        raise RuntimeError(
            f"Cannot reuse {CO2_ACCOMMODATION_STAGE}: provenance is "
            f"incompatible (fields: {', '.join(mismatches)}). Use --force "
            "to archive and rerun this accommodation."
        )
    _validate_co2_output_hashes(stage_dir, manifest)
    if file_sha256(stage_dir / "system_CO2_centered.pdb") != (
        expected_configuration["input_geometry_sha256"]
    ):
        raise RuntimeError(
            f"Cannot reuse {CO2_ACCOMMODATION_STAGE}: its input copy no "
            "longer matches the validated packing."
        )
    validate_co2_accommodated_geometry(
        stage_dir / "system_CO2_centered.pdb",
        stage_dir / "system_CO2_accommodated.pdb",
        source_metadata,
        template_metadata,
        expected_configuration["n_CO2"],
    )
    return manifest


def run_co2_accommodation_stage(
    condition_dir: Path,
    packing_manifest: dict,
    source_metadata: dict,
    template_metadata: dict,
    args,
    *,
    force: bool | None = None,
) -> dict:
    if force is None:
        force = args.force
    pack_stage_dir = condition_dir / CO2_PACK_STAGE
    packing_manifest_path = pack_stage_dir / "stage_manifest.json"
    stage_dir = condition_dir / CO2_ACCOMMODATION_STAGE
    configuration = co2_accommodation_configuration(
        packing_manifest_path, packing_manifest, args
    )
    co2_count = configuration["n_CO2"]
    pack_index = configuration["pack_index"]
    system_name = configuration["system"]

    if (stage_dir / "stage.done").is_file() and not force:
        manifest = validate_completed_co2_accommodation(
            stage_dir, configuration, source_metadata, template_metadata
        )
        print(
            f"  REUSE {system_name}/NCO2_{co2_count:02d}/"
            f"pack_{pack_index:02d}/{CO2_ACCOMMODATION_STAGE}"
        )
        return manifest
    active_entries = _co2_stage_active_entries(stage_dir)
    if active_entries and not force:
        names = ", ".join(path.name for path in active_entries[:6])
        raise RuntimeError(
            f"Existing incomplete or failed {CO2_ACCOMMODATION_STAGE} "
            f"attempt in {stage_dir} ({names}). Inspect it and use --force "
            "to archive it before a new accommodation."
        )

    xtb_exe = shutil.which(args.xtb)
    if xtb_exe is None:
        raise RuntimeError(
            f"xTB executable '{args.xtb}' not found. Use --xtb "
            "/path/to/xtb if needed."
        )
    if active_entries:
        archive = archive_co2_stage_attempt(
            stage_dir, "Replaced explicitly with --force."
        )
        print(f"  ARCHIVE prior {CO2_ACCOMMODATION_STAGE} attempt at {archive}")
    stage_dir.mkdir(parents=True, exist_ok=True)
    initial_pdb = stage_dir / "system_CO2_centered.pdb"
    input_path = stage_dir / "07_CO2_accommodation.inp"
    log_path = stage_dir / "07_CO2_accommodation.out"
    accommodated_pdb = stage_dir / "system_CO2_accommodated.pdb"
    shutil.copy2(pack_stage_dir / "system_CO2_centered.pdb", initial_pdb)
    input_path.write_text(relax_input(
        source_metadata["n_solute_atoms"],
        configuration["wall_radius_bohr"],
        args.co2_accommodation_engine,
    ))
    recorded_command = co2_accommodation_command(args)
    execution_command = [xtb_exe, *recorded_command[1:]]
    started_at = datetime.now(timezone.utc).isoformat()
    runtime = {
        "workflow": "CO2_shell_screening",
        "stage": CO2_ACCOMMODATION_STAGE,
        "status": "running",
        "started_at": started_at,
        "hostname": socket.gethostname(),
        "command": recorded_command,
        "configuration": configuration,
        "execution_resources": co2_execution_resources(args),
    }
    mark_stage_running(stage_dir, runtime)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    env["MKL_NUM_THREADS"] = str(args.threads)
    env["OMP_STACKSIZE"] = "4G"
    print(
        f"  ACCOMMODATE {system_name}/NCO2_{co2_count:02d}/"
        f"pack_{pack_index:02d}"
    )
    print(f"       {' '.join(recorded_command)}")

    returncode = None
    try:
        with log_path.open("w") as log:
            result = subprocess.run(
                execution_command,
                cwd=stage_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                text=True,
            )
        returncode = result.returncode
        if returncode != 0:
            raise RuntimeError(f"xTB returned code {returncode}.")
        diagnostics = relaxation_diagnostics(
            log_path, (stage_dir / "NOT_CONVERGED").exists()
        )
        materialize_relaxed_pdb(
            stage_dir,
            source_metadata["n_source_atoms"] + 3 * co2_count,
            template_name=initial_pdb.name,
            destination_name=accommodated_pdb.name,
            preserve_template_metadata=True,
        )
        geometry_validation = validate_co2_accommodated_geometry(
            initial_pdb,
            accommodated_pdb,
            source_metadata,
            template_metadata,
            co2_count,
        )
        finished_at = datetime.now(timezone.utc).isoformat()
        manifest = {
            **runtime,
            "status": "completed",
            "finished_at": finished_at,
            "returncode": returncode,
            "xtb_version": extract_xtb_version(log_path, accommodated_pdb),
            "source": source_metadata,
            "co2_template": template_metadata,
            "composition": packing_manifest["composition"],
            "wall": packing_manifest["wall"],
            "optimization_diagnostics": diagnostics,
            "geometry_validation": geometry_validation,
            "composition_changed": True,
            "velocities_preserved": False,
            "restart_chain_continued": False,
            "mdrestart_used": False,
            "source_structure_type": "full-droplet representative frame",
            "scientific_note": CO2_WORKFLOW_NOTE,
            "output_sha256": _co2_stage_output_hashes(stage_dir),
        }
        _write_json(stage_dir / "stage_manifest.json", manifest)
        mark_stage_done(stage_dir)
        if diagnostics["converged"] is False:
            print(
                f"  WARNING {CO2_ACCOMMODATION_STAGE}: formal optimization "
                "non-convergence reported; final geometry passed workflow "
                "validation."
            )
        elif diagnostics["converged"] is None:
            print(
                f"  WARNING {CO2_ACCOMMODATION_STAGE}: convergence status "
                "could not be determined robustly; final geometry passed "
                "workflow validation."
            )
        print(f"  OK   {stage_dir}")
        return manifest
    except (OSError, RuntimeError, ValueError) as exc:
        failure = {
            **runtime,
            "status": "failed",
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "returncode": returncode,
            "failure_reason": str(exc),
            "source": source_metadata,
            "co2_template": template_metadata,
            "composition_changed": True,
            "velocities_preserved": False,
            "restart_chain_continued": False,
            "mdrestart_used": False,
            "source_structure_type": "full-droplet representative frame",
            "scientific_note": CO2_WORKFLOW_NOTE,
        }
        _write_json(stage_dir / "stage_manifest.json", failure)
        mark_stage_failed(stage_dir, str(exc))
        raise RuntimeError(
            f"{CO2_ACCOMMODATION_STAGE} failed in {stage_dir}: {exc} "
            f"Outputs were preserved for diagnosis; see {log_path}."
        ) from exc


def co2_md_configuration(
    condition_dir: Path,
    stage: dict,
    packing_manifest: dict,
    accommodation_manifest: dict,
    previous_manifest: dict | None,
    args,
) -> dict:
    packing_configuration = packing_manifest["configuration"]
    composition = packing_manifest["composition"]
    geometry = condition_dir / CO2_ACCOMMODATION_STAGE / (
        "system_CO2_accommodated.pdb"
    )
    accommodation_manifest_path = (
        condition_dir / CO2_ACCOMMODATION_STAGE / "stage_manifest.json"
    )
    input_restart_sha256 = None
    predecessor_manifest_sha256 = None
    if stage["restart"]:
        if previous_manifest is None:
            raise RuntimeError(
                f"{stage['name']} requires a validated predecessor manifest."
            )
        input_restart_sha256 = previous_manifest.get("output_restart_sha256")
        predecessor_manifest_path = (
            condition_dir / stage["restart_from"] / "stage_manifest.json"
        )
        predecessor_manifest_sha256 = file_sha256(predecessor_manifest_path)

    configuration = {
        "workflow": "CO2_shell_screening_MD",
        "stage": stage["name"],
        "system": packing_configuration["system"],
        "n_CO2": packing_configuration["co2_count"],
        "pack_index": packing_configuration["pack_index"],
        "input_geometry": "system_CO2_accommodated.pdb",
        "input_geometry_sha256": file_sha256(geometry),
        "accommodation_manifest_sha256": file_sha256(
            accommodation_manifest_path
        ),
        "predecessor_manifest_sha256": predecessor_manifest_sha256,
        "n_atoms": composition["n_total_atoms"],
        "n_solute_atoms": composition["n_solute_atoms"],
        "n_waters": composition["n_waters"],
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "threads": args.threads,
        "temp_K": stage["temp"],
        "time_ps": stage["time"],
        "step_fs": MD_STEP_FS,
        "dump_fs": stage["dump_fs"],
        "steps": stage["steps"],
        "expected_frames": stage["expected_frames"],
        "nvt": True,
        "velo": True,
        "hmass": 1,
        "shake": 0,
        "sccacc": 1.0,
        "restart": stage["restart"],
        "restart_from": stage.get("restart_from"),
        "continuation": stage["continuation"],
        "continuation_of": stage.get("continuation_of"),
        "velocities_reinitialized": stage["velocities_reinitialized"],
        "mdrestart_input_used": stage["restart"],
        "input_restart_sha256": input_restart_sha256,
        "wall_radius_A": packing_manifest["wall"]["radius_A"],
        "wall_radius_bohr": packing_manifest["wall"]["radius_bohr"],
        "equilibration_time_ps": stage["equilibration_time_ps"],
        "production_time_ps": stage["production_time_ps"],
        "cumulative_production_time_ps": (
            stage["cumulative_production_time_ps"]
        ),
        "source_stage": (
            CO2_ACCOMMODATION_STAGE
            if not stage["restart"]
            else stage["restart_from"]
        ),
    }
    return configuration


def co2_extended_disk_preflight(
    condition_dir: Path,
    stage: dict,
    previous_manifest: dict,
) -> dict:
    previous_trajectory = condition_dir / CO2_SCREEN_STAGE / "xtb.trj"
    if not previous_trajectory.is_file():
        raise RuntimeError(
            f"Cannot estimate {stage['name']} disk space: missing "
            f"{previous_trajectory}."
        )
    previous_size = previous_trajectory.stat().st_size
    previous_integrity = previous_manifest.get("trajectory_integrity", {})
    if previous_integrity.get("size_bytes") != previous_size:
        raise RuntimeError(
            f"Cannot estimate {stage['name']} disk space: stage-09 "
            "trajectory size differs from its manifest."
        )
    if previous_integrity.get("sha256") != file_sha256(previous_trajectory):
        raise RuntimeError(
            f"Cannot estimate {stage['name']} disk space: stage-09 "
            "trajectory hash differs from its manifest."
        )
    previous_time = previous_manifest["configuration"]["time_ps"]
    previous_dump = previous_manifest["configuration"]["dump_fs"]
    scale = (stage["time"] / previous_time) * (
        previous_dump / stage["dump_fs"]
    )
    estimated_trajectory_bytes = math.ceil(previous_size * scale)
    safety_factor = 1.25
    required_free_bytes = math.ceil(
        estimated_trajectory_bytes * safety_factor
    )
    disk = shutil.disk_usage(condition_dir)
    record = {
        "source_stage": CO2_SCREEN_STAGE,
        "source_trajectory": str(previous_trajectory),
        "source_trajectory_bytes": previous_size,
        "duration_and_dump_scale": scale,
        "estimated_trajectory_bytes": estimated_trajectory_bytes,
        "safety_factor": safety_factor,
        "required_free_bytes": required_free_bytes,
        "available_free_bytes": disk.free,
        "sufficient": disk.free >= required_free_bytes,
    }
    return record


def _co2_md_required_outputs(stage: dict) -> list[str]:
    required = [
        "stage.done",
        "stage_manifest.json",
        f"{stage['name']}.inp",
        f"{stage['name']}.out",
        "system_CO2_accommodated.pdb",
        "xtb.trj",
        "mdrestart",
        "xtbmdok",
    ]
    if stage["restart"]:
        required.append("mdrestart.input")
    return required


def validate_completed_co2_md_stage(
    condition_dir: Path,
    stage: dict,
    expected_configuration: dict,
    expected_elements: list[str],
    args,
) -> dict:
    stage_dir = condition_dir / stage["name"]
    conflicts = [
        marker for marker in ("stage.failed", "stage.running")
        if (stage_dir / marker).exists()
    ]
    if conflicts:
        raise RuntimeError(
            f"Completed {stage['name']} has conflicting markers: "
            f"{', '.join(conflicts)}. Inspect it and use --force."
        )
    missing = [
        name for name in _co2_md_required_outputs(stage)
        if not (stage_dir / name).is_file()
    ]
    if missing:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: missing archived outputs "
            f"({', '.join(missing)})."
        )
    manifest = _read_json_dict(
        stage_dir / "stage_manifest.json", "CO2 MD stage manifest"
    )
    if manifest.get("status") != "completed":
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: manifest status is "
            f"{manifest.get('status')!r}."
        )
    recorded_configuration = manifest.get("configuration")
    if not isinstance(recorded_configuration, dict):
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: missing configuration."
        )
    mismatches = configuration_mismatches(
        expected_configuration, recorded_configuration
    )
    if mismatches:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: provenance is incompatible "
            f"(fields: {', '.join(mismatches)}). Use --force to archive "
            "and rerun this stage."
        )
    _validate_co2_output_hashes(stage_dir, manifest)
    geometry = stage_dir / "system_CO2_accommodated.pdb"
    if file_sha256(geometry) != expected_configuration["input_geometry_sha256"]:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: input geometry hash changed."
        )
    geometry_atoms = _validated_pdb_atoms(geometry, "CO2 MD input geometry")
    if len(geometry_atoms) != expected_configuration["n_atoms"]:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: input geometry atom count changed."
        )
    input_path = stage_dir / f"{stage['name']}.inp"
    expected_input = co2_md_input(
        stage, expected_configuration["wall_radius_bohr"]
    )
    if input_path.read_text() != expected_input:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: xTB input content changed."
        )

    input_restart_sha256 = expected_configuration["input_restart_sha256"]
    if stage["restart"]:
        input_restart = stage_dir / "mdrestart.input"
        parse_mdrestart(input_restart, expected_configuration["n_atoms"])
        if file_sha256(input_restart) != input_restart_sha256:
            raise RuntimeError(
                f"Cannot reuse {stage['name']}: mdrestart.input does not "
                "match the validated predecessor restart."
            )
    elif (stage_dir / "mdrestart.input").exists():
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: restart=false stage contains "
            "mdrestart.input."
        )

    output_restart_sha256 = validate_output_restart(
        stage_dir / "mdrestart",
        input_restart_sha256,
        expected_atoms=expected_configuration["n_atoms"],
    )
    if manifest.get("output_restart_sha256") != output_restart_sha256:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: output restart hash differs "
            "from the manifest."
        )
    trajectory = validate_xtb_trajectory(
        stage_dir / "xtb.trj",
        expected_atoms=expected_configuration["n_atoms"],
        expected_frames=stage["expected_frames"],
        require_velocities=True,
        expected_elements=expected_elements,
    )
    recorded_trajectory = manifest.get("trajectory_integrity", {})
    for field in (
        "sha256", "size_bytes", "n_atoms", "frames",
        "expected_nominal_frames", "velocities_present",
    ):
        if recorded_trajectory.get(field) != trajectory.get(field):
            raise RuntimeError(
                f"Cannot reuse {stage['name']}: trajectory integrity field "
                f"{field} differs from the manifest."
            )
    log_validation = inspect_md_log(stage_dir / f"{stage['name']}.out")
    if log_validation["fatal_patterns"]:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: log contains fatal patterns."
        )
    warning_accepted = (
        log_validation["thermostating_problem"]
        and log_validation["normal_exit_of_md"]
        and thermostat_warning_allowed(
            args.thermostat_warning_policy, stage["name"]
        )
    )
    if log_validation["thermostating_problem"] and not warning_accepted:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: unaccepted thermostating problem."
        )
    expected_thermal_result = md_thermal_result(
        stage,
        log_validation,
        args.thermostat_warning_policy,
        warning_accepted,
    )
    recorded_thermal_result = manifest.get("thermal_result")
    if recorded_thermal_result != expected_thermal_result:
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: thermal_result is inconsistent "
            "with the archived log and current warning policy."
        )
    if not isinstance(manifest.get("xtb_version"), str):
        raise RuntimeError(
            f"Cannot reuse {stage['name']}: xTB version is not recorded."
        )
    if stage["name"] in (CO2_SCREEN_STAGE, CO2_EXTENDED_STAGE):
        sampling = co2_sampling_metadata(
            stage["dump_fs"], trajectory["velocities_present"]
        )
        for field, expected_value in sampling.items():
            if manifest.get(field) != expected_value:
                raise RuntimeError(
                    f"Cannot reuse {stage['name']}: spectroscopy metadata "
                    f"field {field} is incompatible."
                )
    return manifest


def run_co2_md_stage(
    condition_dir: Path,
    stage: dict,
    packing_manifest: dict,
    accommodation_manifest: dict,
    previous_manifest: dict | None,
    source_metadata: dict,
    template_metadata: dict,
    args,
    *,
    force: bool,
    allow_run: bool,
) -> dict:
    stage_dir = condition_dir / stage["name"]
    configuration = co2_md_configuration(
        condition_dir,
        stage,
        packing_manifest,
        accommodation_manifest,
        previous_manifest,
        args,
    )
    accommodated_source = (
        condition_dir
        / CO2_ACCOMMODATION_STAGE
        / "system_CO2_accommodated.pdb"
    )
    expected_atoms = _validated_pdb_atoms(
        accommodated_source, "validated CO2 accommodation"
    )
    expected_elements = _element_sequence(expected_atoms)
    if len(expected_atoms) != configuration["n_atoms"]:
        raise RuntimeError(
            f"{stage['name']} composition mismatch: accommodation has "
            f"{len(expected_atoms)} atoms, expected {configuration['n_atoms']}."
        )

    if (stage_dir / "stage.done").is_file() and not force:
        manifest = validate_completed_co2_md_stage(
            condition_dir,
            stage,
            configuration,
            expected_elements,
            args,
        )
        print(
            f"  REUSE {configuration['system']}/NCO2_"
            f"{configuration['n_CO2']:02d}/pack_"
            f"{configuration['pack_index']:02d}/{stage['name']}"
        )
        return manifest

    active_entries = _co2_stage_active_entries(stage_dir)
    if not allow_run:
        raise RuntimeError(
            f"{stage['name']} must already be completed for "
            f"--co2-start-stage; no reusable stage was found in {stage_dir}."
        )
    if active_entries and not force:
        raise RuntimeError(
            f"Existing incomplete or failed {stage['name']} attempt in "
            f"{stage_dir}. Inspect it and use --force to archive and rerun."
        )
    xtb_exe = shutil.which(args.xtb)
    if xtb_exe is None:
        raise RuntimeError(
            f"xTB executable '{args.xtb}' not found. Use --xtb /path/to/xtb."
        )
    if active_entries:
        archive = archive_co2_stage_attempt(
            stage_dir, "Replaced explicitly with --force."
        )
        print(f"  ARCHIVE prior {stage['name']} attempt at {archive}")
    stage_dir.mkdir(parents=True, exist_ok=True)
    geometry = stage_dir / "system_CO2_accommodated.pdb"
    input_path = stage_dir / f"{stage['name']}.inp"
    log_path = stage_dir / f"{stage['name']}.out"
    shutil.copy2(accommodated_source, geometry)
    input_path.write_text(
        co2_md_input(stage, configuration["wall_radius_bohr"])
    )

    if stage["restart"]:
        predecessor_restart = (
            condition_dir / stage["restart_from"] / "mdrestart"
        )
        if not predecessor_restart.is_file():
            raise RuntimeError(
                f"Validated predecessor restart is missing: "
                f"{predecessor_restart}."
            )
        predecessor_hash = file_sha256(predecessor_restart)
        if predecessor_hash != configuration["input_restart_sha256"]:
            raise RuntimeError(
                f"{stage['name']} predecessor restart changed after "
                "validation."
            )
        shutil.copy2(predecessor_restart, stage_dir / "mdrestart")
        shutil.copy2(predecessor_restart, stage_dir / "mdrestart.input")
        if file_sha256(stage_dir / "mdrestart.input") != predecessor_hash:
            raise RuntimeError(
                f"{stage['name']} failed to preserve its input restart hash."
            )
    elif (stage_dir / "mdrestart").exists():
        raise RuntimeError(
            f"{stage['name']} restart=false preflight found stale mdrestart."
        )

    preflight_disk_space = None
    requested_command = co2_md_command(args, stage)
    execution_command = [xtb_exe, *requested_command[1:]]
    started_at = datetime.now(timezone.utc).isoformat()
    runtime = {
        "workflow": "CO2_shell_screening_MD",
        "stage": stage["name"],
        "status": "running",
        "started_at": started_at,
        "hostname": socket.gethostname(),
        "command": execution_command,
        "requested_xtb_executable": args.xtb,
        "configuration": configuration,
        "execution_resources": co2_execution_resources(args),
    }
    mark_stage_running(stage_dir, runtime)
    try:
        if stage["name"] == CO2_EXTENDED_STAGE:
            preflight_disk_space = co2_extended_disk_preflight(
                condition_dir, stage, previous_manifest
            )
            if not preflight_disk_space["sufficient"]:
                raise RuntimeError(
                    f"Insufficient disk space for {stage['name']}: "
                    "estimated requirement with margin is "
                    f"{preflight_disk_space['required_free_bytes']} bytes, "
                    "but only "
                    f"{preflight_disk_space['available_free_bytes']} bytes "
                    "are free. No files were deleted."
                )
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(args.threads)
        env["MKL_NUM_THREADS"] = str(args.threads)
        env["OMP_STACKSIZE"] = "4G"
        print(
            f"  RUN  {configuration['system']}/NCO2_"
            f"{configuration['n_CO2']:02d}/pack_"
            f"{configuration['pack_index']:02d} {stage['name']}"
        )
        print(f"       {' '.join(execution_command)}")
        with log_path.open("w") as log:
            result = subprocess.run(
                execution_command,
                cwd=stage_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                env=env,
                text=True,
            )
        if result.returncode != 0:
            raise RuntimeError(f"xTB returned code {result.returncode}.")
        log_validation = inspect_md_log(log_path)
        if log_validation["fatal_patterns"]:
            raise RuntimeError(
                "Fatal MD pattern(s): "
                + ", ".join(log_validation["fatal_patterns"])
            )
        trajectory = validate_xtb_trajectory(
            stage_dir / "xtb.trj",
            expected_atoms=configuration["n_atoms"],
            expected_frames=stage["expected_frames"],
            require_velocities=True,
            expected_elements=expected_elements,
        )
        input_restart_sha256 = configuration["input_restart_sha256"]
        output_restart_sha256 = validate_output_restart(
            stage_dir / "mdrestart",
            input_restart_sha256,
            expected_atoms=configuration["n_atoms"],
        )
        if not (stage_dir / "xtbmdok").is_file():
            raise RuntimeError("xTB ended without xtbmdok.")
        warning_accepted = (
            log_validation["thermostating_problem"]
            and log_validation["normal_exit_of_md"]
            and thermostat_warning_allowed(
                args.thermostat_warning_policy, stage["name"]
            )
        )
        if log_validation["thermostating_problem"] and not warning_accepted:
            raise RuntimeError("Unaccepted thermostating problem.")
        thermal_result = md_thermal_result(
            stage,
            log_validation,
            args.thermostat_warning_policy,
            warning_accepted,
        )
        xtb_version = extract_xtb_version(log_path, stage_dir / "xtb.trj")
        if xtb_version is None:
            raise RuntimeError(
                "xTB version could not be determined from log/trajectory."
            )
        trajectory = dict(trajectory)
        trajectory["path"] = "xtb.trj"
        manifest = {
            **configuration,
            **runtime,
            "status": "completed",
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "returncode": result.returncode,
            "xtb_version": xtb_version,
            "input_restart_path": (
                "mdrestart.input" if stage["restart"] else None
            ),
            "output_restart_sha256": output_restart_sha256,
            "trajectory_integrity": trajectory,
            "thermal_result": thermal_result,
            "preflight_disk_space": preflight_disk_space,
            "composition_changed_relative_to_aqueous": True,
            "solvation_rebuilt": False,
            "geometry_reoptimized": False,
            "coordinates_recentered": False,
            "output_sha256": _co2_stage_output_hashes(stage_dir),
        }
        if stage["name"] in (CO2_SCREEN_STAGE, CO2_EXTENDED_STAGE):
            manifest.update(co2_sampling_metadata(
                stage["dump_fs"], trajectory["velocities_present"]
            ))
        _write_json(stage_dir / "stage_manifest.json", manifest)
        mark_stage_done(stage_dir)
        if warning_accepted:
            print(
                f"  WARNING {stage['name']}: accepted thermostating problem "
                "after complete output/restart validation."
            )
        else:
            print(f"  OK   {stage['name']}")
        return manifest
    except (OSError, RuntimeError, ValueError) as exc:
        failure = {
            **configuration,
            **runtime,
            "status": "failed",
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "failure_reason": str(exc),
            "preflight_disk_space": preflight_disk_space,
            "composition_changed_relative_to_aqueous": True,
            "output_sha256": _co2_stage_output_hashes(stage_dir),
        }
        _write_json(stage_dir / "stage_manifest.json", failure)
        mark_stage_failed(stage_dir, str(exc))
        raise RuntimeError(
            f"{stage['name']} failed in {stage_dir}: {exc} Outputs were "
            "preserved for diagnosis."
        ) from exc


PACKING_SUMMARY_FIELDS = [
    "system", "n_CO2", "pack_index", "seed", "placement_mode",
    "direction_atom_index", "direction_atom_element", "target_distance_A",
    "target_radius_A", "targeted_CO2_index",
    "targeted_carbon_atom_index", "targeted_Zn_C_A",
    "targeted_target_distance_A", "targeted_direction_angle_deg",
    "site_direction_valid", "shell_inner_A", "shell_outer_A", "n_atoms",
    "min_Zn_C_A", "mean_Zn_C_A", "max_Zn_C_A",
    "min_CO2_source_A", "min_CO2_CO2_A", "packing_status",
]

ACCOMMODATION_SUMMARY_FIELDS = [
    "system", "n_CO2", "pack_index", "accommodation_status", "converged",
    "optimization_cycles", "min_Zn_C_before_A", "min_Zn_C_after_A",
    "mean_Zn_C_before_A", "mean_Zn_C_after_A",
    "min_waterO_C_before_A", "min_waterO_C_after_A", "initial_energy",
    "final_energy",
]


def _summary_value(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return f"{value:.10g}"
    return value


def _co2_effective_stage_status(stage_dir: Path, manifest: dict) -> str:
    markers = [
        name for name in ("stage.done", "stage.failed", "stage.running")
        if (stage_dir / name).exists()
    ]
    if len(markers) > 1:
        return "invalid_conflicting_markers"
    if markers == ["stage.done"]:
        return (
            "completed"
            if manifest.get("status") == "completed"
            else "invalid_done_manifest"
        )
    if markers == ["stage.failed"]:
        return "failed"
    if markers == ["stage.running"]:
        return "running"
    return str(manifest.get("status") or "incomplete")


def _write_tsv(path: Path, fields: list[str], rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _summary_value(row.get(field)) for field in fields})


def rebuild_co2_summaries(system_dir: Path):
    packing_rows = []
    accommodation_rows = []
    condition_dirs = sorted(
        system_dir.glob("NCO2_*/pack_*"),
        key=lambda path: (
            int(path.parent.name.removeprefix("NCO2_")),
            int(path.name.removeprefix("pack_")),
        ),
    )
    for condition_dir in condition_dirs:
        pack_stage_dir = condition_dir / CO2_PACK_STAGE
        pack_manifest_path = pack_stage_dir / "stage_manifest.json"
        if pack_manifest_path.is_file():
            manifest = _read_json_dict(
                pack_manifest_path, "CO2 packing manifest"
            )
            config = manifest.get("configuration", {})
            validation = manifest.get("packing_validation", {})
            metrics = validation.get("metrics", {})
            composition = manifest.get("composition", {})
            site_direction = manifest.get("site_direction", {})
            placement_mode = config.get("placement_mode", "random-shell")
            site_directed = placement_mode == "site-directed"
            packing_rows.append({
                "system": config.get("system"),
                "n_CO2": config.get("co2_count"),
                "pack_index": config.get("pack_index"),
                "seed": config.get("packmol_seed"),
                "placement_mode": placement_mode,
                "direction_atom_index": (
                    site_direction.get("direction_atom_index")
                    if site_directed else None
                ),
                "direction_atom_element": (
                    site_direction.get("direction_atom_element")
                    if site_directed else None
                ),
                "target_distance_A": (
                    site_direction.get("target_distance_A")
                    if site_directed else None
                ),
                "target_radius_A": (
                    site_direction.get("target_radius_A")
                    if site_directed else None
                ),
                "targeted_CO2_index": (
                    validation.get("targeted_CO2_index")
                    if site_directed else None
                ),
                "targeted_carbon_atom_index": (
                    validation.get("targeted_carbon_atom_index")
                    if site_directed else None
                ),
                "targeted_Zn_C_A": (
                    validation.get("targeted_Zn_C_distance_A")
                    if site_directed else None
                ),
                "targeted_target_distance_A": (
                    validation.get("targeted_target_point_distance_A")
                    if site_directed else None
                ),
                "targeted_direction_angle_deg": (
                    validation.get("targeted_direction_angle_degrees")
                    if site_directed else None
                ),
                "site_direction_valid": (
                    validation.get("site_direction_valid")
                    if site_directed else None
                ),
                "shell_inner_A": config.get("shell_inner_A"),
                "shell_outer_A": config.get("shell_outer_A"),
                "n_atoms": composition.get("n_total_atoms"),
                "min_Zn_C_A": metrics.get("zn_c_min_A"),
                "mean_Zn_C_A": metrics.get("zn_c_mean_A"),
                "max_Zn_C_A": metrics.get("zn_c_max_A"),
                "min_CO2_source_A": validation.get(
                    "minimum_CO2_source_distance_A"
                ),
                "min_CO2_CO2_A": validation.get(
                    "minimum_CO2_CO2_distance_A"
                ),
                "packing_status": _co2_effective_stage_status(
                    pack_stage_dir, manifest
                ),
            })

        accommodation_stage_dir = condition_dir / CO2_ACCOMMODATION_STAGE
        accommodation_manifest_path = accommodation_stage_dir / "stage_manifest.json"
        if accommodation_manifest_path.is_file():
            manifest = _read_json_dict(
                accommodation_manifest_path, "CO2 accommodation manifest"
            )
            config = manifest.get("configuration", {})
            diagnostics = manifest.get("optimization_diagnostics", {})
            geometry = manifest.get("geometry_validation", {})
            before = geometry.get("metrics_before", {})
            after = geometry.get("metrics_after", {})
            accommodation_rows.append({
                "system": config.get("system"),
                "n_CO2": config.get("n_CO2"),
                "pack_index": config.get("pack_index"),
                "accommodation_status": _co2_effective_stage_status(
                    accommodation_stage_dir, manifest
                ),
                "converged": diagnostics.get("converged"),
                "optimization_cycles": diagnostics.get("cycles"),
                "min_Zn_C_before_A": before.get("zn_c_min_A"),
                "min_Zn_C_after_A": after.get("zn_c_min_A"),
                "mean_Zn_C_before_A": before.get("zn_c_mean_A"),
                "mean_Zn_C_after_A": after.get("zn_c_mean_A"),
                "min_waterO_C_before_A": before.get("water_o_c_min_A"),
                "min_waterO_C_after_A": after.get("water_o_c_min_A"),
                "initial_energy": diagnostics.get("initial_energy_Eh"),
                "final_energy": diagnostics.get("final_energy_Eh"),
            })

    _write_tsv(system_dir / "packing_summary.tsv", PACKING_SUMMARY_FIELDS, packing_rows)
    accommodation_summary = system_dir / "accommodation_summary.tsv"
    if accommodation_rows or accommodation_summary.exists():
        _write_tsv(
            accommodation_summary,
            ACCOMMODATION_SUMMARY_FIELDS,
            accommodation_rows,
        )


CO2_MD_SUMMARY_FIELDS = [
    "system", "n_CO2", "pack_index",
    "equil_status", "equil_time_ps", "equil_dump_fs",
    "equil_average_temperature_K", "equil_frames",
    "equil_restart_sha256",
    "screen_status", "screen_time_ps", "screen_dump_fs",
    "screen_average_temperature_K", "screen_frames",
    "screen_trajectory_path", "screen_trajectory_sha256",
    "screen_restart_sha256",
    "extended_status", "extended_time_ps", "extended_dump_fs",
    "extended_average_temperature_K", "extended_frames",
    "extended_trajectory_path", "extended_trajectory_sha256",
    "extended_restart_sha256", "cumulative_production_time_ps",
    "spectroscopy_sampling_ready",
]


def _co2_md_summary_stage_values(
    condition_dir: Path,
    stage_name: str,
) -> tuple[dict | None, str]:
    stage_dir = condition_dir / stage_name
    manifest_path = stage_dir / "stage_manifest.json"
    if not manifest_path.is_file():
        return None, "not_requested"
    manifest = _read_json_dict(manifest_path, "CO2 MD stage manifest")
    return manifest, _co2_effective_stage_status(stage_dir, manifest)


def rebuild_co2_md_summary(system_dir: Path):
    rows = []
    condition_dirs = sorted(
        system_dir.glob("NCO2_*/pack_*"),
        key=lambda path: (
            int(path.parent.name.removeprefix("NCO2_")),
            int(path.name.removeprefix("pack_")),
        ),
    )
    for condition_dir in condition_dirs:
        equil, equil_status = _co2_md_summary_stage_values(
            condition_dir, CO2_EQUIL_STAGE
        )
        screen, screen_status = _co2_md_summary_stage_values(
            condition_dir, CO2_SCREEN_STAGE
        )
        extended, extended_status = _co2_md_summary_stage_values(
            condition_dir, CO2_EXTENDED_STAGE
        )
        if equil is None and screen is None and extended is None:
            continue
        identity = next(
            manifest for manifest in (equil, screen, extended)
            if manifest is not None
        )

        def config(manifest):
            return manifest.get("configuration", {}) if manifest else {}

        def thermal(manifest):
            return manifest.get("thermal_result", {}) if manifest else {}

        def trajectory(manifest):
            return manifest.get("trajectory_integrity", {}) if manifest else {}

        equil_config = config(equil)
        screen_config = config(screen)
        extended_config = config(extended)
        screen_trajectory = trajectory(screen)
        extended_trajectory = trajectory(extended)
        production_manifest = (
            extended
            if extended and extended.get("trajectory_integrity")
            else screen
        )
        rows.append({
            "system": identity.get("system"),
            "n_CO2": identity.get("n_CO2"),
            "pack_index": identity.get("pack_index"),
            "equil_status": equil_status,
            "equil_time_ps": equil_config.get("time_ps"),
            "equil_dump_fs": equil_config.get("dump_fs"),
            "equil_average_temperature_K": thermal(equil).get(
                "average_temperature_K"
            ),
            "equil_frames": trajectory(equil).get("frames"),
            "equil_restart_sha256": (
                equil.get("output_restart_sha256") if equil else None
            ),
            "screen_status": screen_status,
            "screen_time_ps": screen_config.get("time_ps"),
            "screen_dump_fs": screen_config.get("dump_fs"),
            "screen_average_temperature_K": thermal(screen).get(
                "average_temperature_K"
            ),
            "screen_frames": screen_trajectory.get("frames"),
            "screen_trajectory_path": (
                str(stage_relative_path(condition_dir, CO2_SCREEN_STAGE, "xtb.trj"))
                if screen_trajectory else None
            ),
            "screen_trajectory_sha256": screen_trajectory.get("sha256"),
            "screen_restart_sha256": (
                screen.get("output_restart_sha256") if screen else None
            ),
            "extended_status": extended_status,
            "extended_time_ps": extended_config.get("time_ps"),
            "extended_dump_fs": extended_config.get("dump_fs"),
            "extended_average_temperature_K": thermal(extended).get(
                "average_temperature_K"
            ),
            "extended_frames": extended_trajectory.get("frames"),
            "extended_trajectory_path": (
                str(stage_relative_path(
                    condition_dir, CO2_EXTENDED_STAGE, "xtb.trj"
                ))
                if extended_trajectory else None
            ),
            "extended_trajectory_sha256": extended_trajectory.get("sha256"),
            "extended_restart_sha256": (
                extended.get("output_restart_sha256") if extended else None
            ),
            "cumulative_production_time_ps": config(
                production_manifest
            ).get("cumulative_production_time_ps"),
            "spectroscopy_sampling_ready": (
                production_manifest.get("spectroscopy_sampling_ready")
                if production_manifest else None
            ),
        })
    summary_path = system_dir / "co2_md_summary.tsv"
    if rows or summary_path.exists():
        _write_tsv(summary_path, CO2_MD_SUMMARY_FIELDS, rows)


def stage_relative_path(
    condition_dir: Path,
    stage_name: str,
    filename: str,
) -> Path:
    return Path(condition_dir.parent.name) / condition_dir.name / stage_name / filename


def _co2_marker_status(stage_dir: Path) -> str:
    if (stage_dir / "stage.failed").exists():
        return "FAILED"
    if (stage_dir / "stage.running").exists():
        return "RUNNING"
    if (stage_dir / "stage.done").exists():
        return "OK"
    return "--"


def print_co2_workflow_status(
    system_dir: Path,
    condition_specs: list[dict],
    args,
):
    print("\nCO2 workflow summary")
    print("NCO2  pack  06      07      08      09      10")
    for spec in condition_specs:
        condition_dir = spec["condition_dir"]
        statuses = [
            _co2_marker_status(condition_dir / stage_name)
            for stage_name in (
                CO2_PACK_STAGE,
                CO2_ACCOMMODATION_STAGE,
                *CO2_MD_STAGE_NAMES,
            )
        ]
        print(
            f"{spec['co2_count']:<5d}  {spec['pack_index']:<4d}  "
            + "  ".join(f"{status:<6s}" for status in statuses)
        )
    resources = co2_execution_resources(args)
    print(f"\nConcurrent jobs: {resources['co2_parallel_jobs']}")
    print(f"Threads/job: {resources['xtb_threads_per_job']}")
    print(f"Maximum requested threads: {resources['maximum_requested_cpus']}")


def find_co2_pdb(explicit: Path | None = None) -> Path:
    if explicit is not None:
        if not explicit.is_file():
            raise RuntimeError(f"CO2 template PDB not found: {explicit}")
        return explicit.resolve()
    candidate = ROOT / "co2.pdb"
    if candidate.is_file():
        return candidate.resolve()
    raise RuntimeError(
        "No CO2 template PDB was supplied. Use --co2-pdb /path/to/co2.pdb "
        "or create co2.pdb in the repository root."
    )


def run_co2_condition(
    condition_spec: dict,
    system_name: str,
    source_pdb: Path,
    co2_pdb: Path,
    source_metadata: dict,
    template_metadata: dict,
    args,
) -> dict:
    """Run one independent CO2 branch; stages inside it stay sequential."""
    condition_dir = condition_spec["condition_dir"]
    co2_count = condition_spec["co2_count"]
    pack_index = condition_spec["pack_index"]
    if (
        args.co2_start_stage is not None
        and not (condition_dir / CO2_PACK_STAGE / "stage.done").is_file()
    ):
        raise RuntimeError(
            f"--co2-start-stage requires an existing validated "
            f"{CO2_PACK_STAGE} in {condition_dir}."
        )
    packing_manifest = run_co2_pack_stage(
        condition_dir,
        system_name,
        source_pdb,
        co2_pdb,
        source_metadata,
        template_metadata,
        co2_count,
        pack_index,
        args,
    )
    result = {
        "condition": f"NCO2_{co2_count:02d}/pack_{pack_index:02d}",
        CO2_PACK_STAGE: "completed",
    }
    if not args.run:
        return result

    accommodation_force = (
        args.force
        and args.co2_start_stage in (None, CO2_ACCOMMODATION_STAGE)
    )
    if args.co2_start_stage in CO2_MD_STAGE_NAMES:
        accommodation_configuration = co2_accommodation_configuration(
            condition_dir / CO2_PACK_STAGE / "stage_manifest.json",
            packing_manifest,
            args,
        )
        accommodation_manifest = validate_completed_co2_accommodation(
            condition_dir / CO2_ACCOMMODATION_STAGE,
            accommodation_configuration,
            source_metadata,
            template_metadata,
        )
        print(
            f"  REUSE {system_name}/NCO2_{co2_count:02d}/"
            f"pack_{pack_index:02d}/{CO2_ACCOMMODATION_STAGE}"
        )
    else:
        accommodation_manifest = run_co2_accommodation_stage(
            condition_dir,
            packing_manifest,
            source_metadata,
            template_metadata,
            args,
            force=accommodation_force,
        )
    result[CO2_ACCOMMODATION_STAGE] = "completed"

    if not args.co2_md:
        return result
    stages = co2_md_stages(args)
    stop_index = 2 if args.co2_extended else 1
    if args.co2_start_stage in CO2_MD_STAGE_NAMES:
        start_index = CO2_MD_STAGE_NAMES.index(args.co2_start_stage)
    else:
        start_index = 0
    previous_manifest = None
    for index, stage in enumerate(stages[:stop_index + 1]):
        allow_run = index >= start_index
        previous_manifest = run_co2_md_stage(
            condition_dir,
            stage,
            packing_manifest,
            accommodation_manifest,
            previous_manifest,
            source_metadata,
            template_metadata,
            args,
            force=(args.force and allow_run),
            allow_run=allow_run,
        )
        result[stage["name"]] = "completed"
    return result


def run_co2_workflow(args, system_name: str):
    source_pdb = args.co2_source_pdb.resolve()
    co2_pdb = find_co2_pdb(args.co2_pdb)
    expected_solute_pdb = SYSTEMS[system_name]["solute"].resolve()
    source_metadata = validate_co2_source_pdb(
        source_pdb, args.co2_solute_atoms, expected_solute_pdb
    )
    template_metadata = validate_co2_template(co2_pdb)
    if args.co2_start_stage is None and shutil.which(args.packmol) is None:
        raise RuntimeError(
            f"Packmol executable '{args.packmol}' not found in PATH. "
            "Use --packmol /path/to/packmol if needed."
        )
    if args.run and shutil.which(args.xtb) is None:
        raise RuntimeError(
            f"xTB executable '{args.xtb}' not found. Use --xtb "
            "/path/to/xtb if needed."
        )
    system_dir = args.co2_project.resolve() / system_name
    system_dir.mkdir(parents=True, exist_ok=True)

    print(f"CO2 shell-screen project: {args.co2_project.resolve()}")
    print(f"System: {system_name}")
    print(f"Source full droplet: {source_pdb}")
    print(f"CO2 template: {co2_pdb}")
    print(
        f"Shell around Zn: {args.co2_shell_inner:.3f}-"
        f"{args.co2_shell_outer:.3f} A"
    )
    print(f"CO2 placement mode: {args.co2_placement_mode}")
    print(CO2_WORKFLOW_NOTE)
    resources = co2_execution_resources(args)
    print(f"CO2 concurrent jobs     : {resources['co2_parallel_jobs']}")
    print(f"xTB threads per job     : {resources['xtb_threads_per_job']}")
    print(f"maximum requested CPUs  : {resources['maximum_requested_cpus']}")
    if resources["oversubscription_warning"]:
        print(
            "WARNING: CO2 parallel jobs x threads exceeds os.cpu_count() "
            f"({resources['detected_cpu_count']}). Execution will continue "
            "because HPC affinity/allocation may differ."
        )

    condition_specs = [
        {
            "co2_count": co2_count,
            "pack_index": pack_index,
            "condition_dir": co2_condition_dir(
                args.co2_project.resolve(),
                system_name,
                co2_count,
                pack_index,
            ),
        }
        for co2_count in sorted(args.co2_counts)
        for pack_index in range(1, args.co2_pack_replicas + 1)
    ]
    failures = []
    try:
        if args.co2_parallel_jobs == 1:
            for spec in condition_specs:
                run_co2_condition(
                    spec,
                    system_name,
                    source_pdb,
                    co2_pdb,
                    source_metadata,
                    template_metadata,
                    args,
                )
        else:
            with ThreadPoolExecutor(
                max_workers=args.co2_parallel_jobs,
                thread_name_prefix="co2-condition",
            ) as executor:
                future_to_spec = {
                    executor.submit(
                        run_co2_condition,
                        spec,
                        system_name,
                        source_pdb,
                        co2_pdb,
                        source_metadata,
                        template_metadata,
                        args,
                    ): spec
                    for spec in condition_specs
                }
                for future in as_completed(future_to_spec):
                    spec = future_to_spec[future]
                    try:
                        future.result()
                    except Exception as exc:
                        label = (
                            f"NCO2_{spec['co2_count']:02d}/"
                            f"pack_{spec['pack_index']:02d}"
                        )
                        failures.append((label, exc))
                        print(f"  FAILED {label}: {exc}")
    finally:
        rebuild_co2_summaries(system_dir)
        rebuild_co2_md_summary(system_dir)
        print_co2_workflow_status(system_dir, condition_specs, args)

    if failures:
        detail = "; ".join(
            f"{label}: {error}" for label, error in failures
        )
        raise RuntimeError(
            f"{len(failures)} independent CO2 condition(s) failed after "
            f"all submitted workers finished: {detail}"
        )

    if not args.run:
        print(
            "CO2 shell packing finished and was validated. Inspect "
            "system_CO2_centered.pdb and packing_summary.tsv before running "
            "the xTB accommodation with the same command plus --run."
        )
    elif args.co2_extended:
        print("CO2 workflow finished through 10_CO2_298K_extended.")
    elif args.co2_md:
        print("CO2 workflow finished through 09_CO2_298K_screen.")
    else:
        print("CO2 shell packing and structural accommodation finished.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Prepare spherical Packmol droplets and optionally run "
            "xTB MD E2 thermalization + screening, with an opt-in "
            "20 ps continuation, or run the independent CO2 shell screen."
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
        "--dry-run",
        action="store_true",
        help=(
            "Validate and display the opt-in 05_298K_extended plan without "
            "writing files or calling xTB; requires "
            "--start-stage 05_298K_extended."
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
        help=(
            "OMP threads for each new xTB run; recorded as execution "
            "provenance but not used to invalidate historical stages "
            "(default: 8)."
        ),
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
        "--skip-relax",
        action="store_true",
        help="Skip solvent-only 00_relax and start MD from system_centered.pdb.",
    )

    p.add_argument(
        "--relax-level",
        default="loose",
        choices=[
            "crude", "sloppy", "loose", "lax", "normal", "tight",
            "vtight", "extreme",
        ],
        help="xTB optimization level for 00_relax (default: loose).",
    )

    p.add_argument(
        "--relax-cycles",
        type=int,
        default=30,
        help="Maximum optimization cycles for 00_relax (default: 30).",
    )

    p.add_argument(
        "--relax-engine",
        default="auto",
        choices=["auto", "rf", "lbfgs", "inertial"],
        help=(
            "Optimization engine for 00_relax; inertial selects native xTB "
            "FIRE (default: auto)."
        ),
    )

    p.add_argument(
        "--thermostat-warning-policy",
        default="allow",
        choices=["strict", "ramp", "allow"],
        help=(
            "Record an isolated xTB 'thermostating problem' as a warning "
            "after all integrity checks pass. Legacy strict/ramp values are "
            "accepted but mapped to allow (default: allow)."
        ),
    )

    p.add_argument(
        "--resume",
        "--resume-thermostat-warning",
        dest="resume",
        action="store_true",
        help=(
            "Resume an existing project without Packmol/preparation; validate "
            "and reuse compatible stages, including safe promotion of a "
            "thermostat-warning-only failed stage."
        ),
    )

    p.add_argument(
        "--start-stage",
        choices=EXECUTION_STAGES,
        help=(
            "Start execution from an existing prepared stage. In particular, "
            "04_298K_screen may explicitly reuse a valid historical 03 "
            "restart, while 05_298K_extended is the opt-in 20 ps continuation "
            "from a validated 04 restart."
        ),
    )

    p.add_argument(
        "--force",
        action="store_true",
        help="Re-run xTB stages even if stage.done exists.",
    )

    co2 = p.add_argument_group(
        "independent CO2 shell screening",
        "Build full-droplet + CO2 conditions, optimize mobile water/CO2, "
        "and optionally run their independent 298 K MD branches.",
    )
    co2.add_argument(
        "--co2-shell-screen",
        action="store_true",
        help="Activate the independent CO2 workflow (06 onward).",
    )
    co2.add_argument(
        "--co2-source-pdb",
        type=Path,
        default=None,
        metavar="FILE",
        help="Full representative droplet PDB: Zn(His)2 plus explicit waters.",
    )
    co2.add_argument(
        "--co2-pdb",
        type=Path,
        default=None,
        metavar="FILE",
        help="Three-atom molecular CO2 template PDB (default: ROOT/co2.pdb).",
    )
    co2.add_argument(
        "--co2-counts",
        type=int,
        nargs="+",
        default=None,
        metavar="N",
        help="One or more CO2 molecule counts, e.g. 1 2 4 8.",
    )
    co2.add_argument(
        "--co2-placement-mode",
        choices=["random-shell", "site-directed"],
        default="random-shell",
        help=(
            "Initial CO2 Packmol placement: historical random shell or one "
            "site-directed molecule plus random-shell background "
            "(default: random-shell)."
        ),
    )
    co2.add_argument(
        "--co2-direction-atom",
        type=int,
        default=None,
        metavar="INT",
        help=(
            "1-based ordinal ATOM/HETATM index in --co2-source-pdb used "
            "to define site direction. Required for site-directed mode."
        ),
    )
    co2.add_argument(
        "--co2-target-distance",
        type=float,
        default=None,
        metavar="ANGSTROM",
        help=(
            "Zn-to-target distance in A; default is the midpoint of the "
            "requested Zn shell."
        ),
    )
    co2.add_argument(
        "--co2-target-radius",
        type=float,
        default=1.5,
        metavar="ANGSTROM",
        help="Site-directed target-sphere radius in A (default: 1.5).",
    )
    co2.add_argument(
        "--co2-shell-inner",
        type=float,
        default=4.0,
        help="Inner Zn-to-CO2-carbon packing radius in A (default: 4.0).",
    )
    co2.add_argument(
        "--co2-shell-outer",
        type=float,
        default=6.0,
        help="Outer Zn-to-CO2-carbon packing radius in A (default: 6.0).",
    )
    co2.add_argument(
        "--co2-pack-replicas",
        type=int,
        default=1,
        help="Independent Packmol placements for each CO2 count (default: 1).",
    )
    co2.add_argument(
        "--co2-project",
        type=Path,
        default=ROOT / "co2_screening",
        help="CO2 workflow output directory (default: ROOT/co2_screening).",
    )
    co2.add_argument(
        "--co2-solute-atoms",
        type=int,
        default=39,
        help="Leading biomimetic atoms fixed during accommodation (default: 39).",
    )
    co2.add_argument(
        "--co2-seed-base",
        type=int,
        default=314159,
        help="Deterministic independent base for CO2 Packmol seeds (default: 314159).",
    )
    co2.add_argument(
        "--co2-repack",
        action="store_true",
        help="Archive and regenerate only the CO2 packing stage.",
    )
    co2.add_argument(
        "--co2-accommodation-level",
        default="loose",
        choices=[
            "crude", "sloppy", "loose", "lax", "normal", "tight",
            "vtight", "extreme",
        ],
        help="xTB optimization level for CO2 accommodation (default: loose).",
    )
    co2.add_argument(
        "--co2-accommodation-cycles",
        type=int,
        default=30,
        help="Maximum CO2 accommodation optimization cycles (default: 30).",
    )
    co2.add_argument(
        "--co2-accommodation-engine",
        default="auto",
        choices=["auto", "rf", "lbfgs", "inertial"],
        help="CO2 accommodation optimizer; inertial selects xTB FIRE (default: auto).",
    )
    co2.add_argument(
        "--co2-md",
        action="store_true",
        help="After 06/07, run 08 equilibration and 09 CO2 screening MD.",
    )
    co2.add_argument(
        "--co2-extended",
        action="store_true",
        help="Also run the direct 10_CO2_298K_extended continuation.",
    )
    co2.add_argument(
        "--co2-equil-time-ps",
        type=float,
        default=1.0,
        help="Stage-08 equilibration duration in ps (default: 1.0).",
    )
    co2.add_argument(
        "--co2-screen-time-ps",
        type=float,
        default=5.0,
        help="Stage-09 production/screening duration in ps (default: 5.0).",
    )
    co2.add_argument(
        "--co2-extended-time-ps",
        type=float,
        default=20.0,
        help="Stage-10 additional production duration in ps (default: 20.0).",
    )
    co2.add_argument(
        "--co2-equil-dump-fs",
        type=float,
        default=10.0,
        help="Stage-08 trajectory dump interval in fs (default: 10.0).",
    )
    co2.add_argument(
        "--co2-production-dump-fs",
        type=float,
        default=2.0,
        help="Stage-09/10 trajectory dump interval in fs (default: 2.0).",
    )
    co2.add_argument(
        "--co2-parallel-jobs",
        type=int,
        default=1,
        help="Maximum independent CO2 condition pipelines in parallel (default: 1).",
    )
    co2.add_argument(
        "--co2-start-stage",
        choices=CO2_START_STAGE_CHOICES,
        default=None,
        help="Validate predecessors and begin/reuse the selected CO2 stage.",
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
    if args.thermostat_warning_policy != "allow":
        print(
            "WARNING: --thermostat-warning-policy "
            f"{args.thermostat_warning_policy} is deprecated by the revised "
            "E2 protocol and will be treated as 'allow'."
        )
        args.thermostat_warning_policy = "allow"

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

    if args.relax_cycles < 1:
        raise SystemExit("--relax-cycles must be >= 1.")

    if args.co2_shell_screen:
        if args.system is None or len(args.system) != 1:
            raise SystemExit(
                "--co2-shell-screen requires exactly one system selected "
                "with --system NAME."
            )
        if args.co2_source_pdb is None:
            raise SystemExit(
                "--co2-shell-screen requires --co2-source-pdb FILE."
            )
        if not args.co2_source_pdb.is_file():
            raise SystemExit(
                f"CO2 source PDB not found: {args.co2_source_pdb}"
            )
        if args.co2_counts is None:
            raise SystemExit(
                "--co2-shell-screen requires --co2-counts N [N ...]."
            )
        if any(count < 1 for count in args.co2_counts):
            raise SystemExit("Every --co2-counts value must be >= 1.")
        if len(set(args.co2_counts)) != len(args.co2_counts):
            raise SystemExit("--co2-counts must not contain duplicates.")
        if (
            not math.isfinite(args.co2_shell_inner)
            or args.co2_shell_inner <= 0
        ):
            raise SystemExit("--co2-shell-inner must be > 0.")
        if (
            not math.isfinite(args.co2_shell_outer)
            or args.co2_shell_outer <= args.co2_shell_inner
        ):
            raise SystemExit(
                "--co2-shell-outer must be greater than --co2-shell-inner."
            )
        if (
            args.co2_placement_mode == "site-directed"
            and args.co2_direction_atom is None
        ):
            raise SystemExit(
                "--co2-placement-mode site-directed requires "
                "--co2-direction-atom INT."
            )
        if (
            args.co2_direction_atom is not None
            and args.co2_direction_atom < 1
        ):
            raise SystemExit("--co2-direction-atom must be >= 1.")
        if (
            not math.isfinite(args.co2_target_radius)
            or args.co2_target_radius <= 0
        ):
            raise SystemExit("--co2-target-radius must be > 0.")
        target_distance_A = co2_effective_target_distance(
            args.co2_shell_inner,
            args.co2_shell_outer,
            args.co2_target_distance,
        )
        if not math.isfinite(target_distance_A) or target_distance_A <= 0:
            raise SystemExit("--co2-target-distance must be > 0.")
        if not (
            args.co2_shell_inner
            <= target_distance_A
            <= args.co2_shell_outer
        ):
            raise SystemExit(
                "--co2-target-distance must lie within the requested "
                "Zn shell."
            )
        if args.co2_pack_replicas < 1:
            raise SystemExit("--co2-pack-replicas must be >= 1.")
        if args.co2_solute_atoms < 1:
            raise SystemExit("--co2-solute-atoms must be >= 1.")
        if args.co2_seed_base < 1:
            raise SystemExit("--co2-seed-base must be >= 1.")
        seeds = [
            co2_pack_seed(args.co2_seed_base, count, pack_index)
            for count in args.co2_counts
            for pack_index in range(1, args.co2_pack_replicas + 1)
        ]
        if len(seeds) != len(set(seeds)):
            raise SystemExit(
                "The requested CO2 count/packing combinations produce a "
                "Packmol seed collision; reduce the extreme count/replica "
                "range."
            )
        if args.co2_accommodation_cycles < 1:
            raise SystemExit("--co2-accommodation-cycles must be >= 1.")
        for option, value in (
            ("--co2-equil-time-ps", args.co2_equil_time_ps),
            ("--co2-screen-time-ps", args.co2_screen_time_ps),
            ("--co2-extended-time-ps", args.co2_extended_time_ps),
            ("--co2-equil-dump-fs", args.co2_equil_dump_fs),
            ("--co2-production-dump-fs", args.co2_production_dump_fs),
        ):
            if not math.isfinite(value) or value <= 0:
                raise SystemExit(f"{option} must be > 0.")
        if args.co2_parallel_jobs < 1:
            raise SystemExit("--co2-parallel-jobs must be >= 1.")
        if args.co2_md and not args.run:
            raise SystemExit("--co2-md requires --run.")
        if args.co2_extended and not args.co2_md:
            raise SystemExit("--co2-extended requires --co2-md.")
        if args.co2_extended and not args.run:
            raise SystemExit("--co2-extended requires --run.")
        if args.co2_start_stage is not None and not args.run:
            raise SystemExit("--co2-start-stage requires --run.")
        if (
            args.co2_start_stage in CO2_MD_STAGE_NAMES
            and not args.co2_md
        ):
            raise SystemExit(
                "CO2 start stages 08-10 require --co2-md."
            )
        if (
            args.co2_start_stage == CO2_EXTENDED_STAGE
            and not args.co2_extended
        ):
            raise SystemExit(
                f"--co2-start-stage {CO2_EXTENDED_STAGE} requires "
                "--co2-extended."
            )
        if args.co2_start_stage is not None and args.co2_repack:
            raise SystemExit(
                "--co2-start-stage cannot be combined with --co2-repack; "
                "the existing packing is a required predecessor."
            )
        try:
            stages = co2_md_stages(args)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        for stage in stages:
            exact_dump_steps = stage["dump_fs"] / MD_STEP_FS
            if not math.isclose(
                exact_dump_steps,
                round(exact_dump_steps),
                abs_tol=1.0e-9,
            ):
                raise SystemExit(
                    f"{stage['name']} dump interval ({stage['dump_fs']} fs) "
                    f"is not an integral number of {MD_STEP_FS} fs steps."
                )
        if args.force and not args.run:
            raise SystemExit(
                "In CO2 mode, --force reruns 07_CO2_accommodation and "
                "therefore requires --run. Use --co2-repack for stage 06."
            )
        incompatible = []
        if args.resume:
            incompatible.append("--resume")
        if args.start_stage is not None:
            incompatible.append("--start-stage")
        if args.dry_run:
            incompatible.append("--dry-run")
        if args.repack:
            incompatible.append("--repack")
        if incompatible:
            raise SystemExit(
                "--co2-shell-screen cannot be combined with "
                + ", ".join(incompatible)
                + "."
            )
        return

    if (
        args.co2_repack
        or args.co2_source_pdb is not None
        or args.co2_pdb is not None
        or args.co2_counts
        or args.co2_placement_mode != "random-shell"
        or args.co2_direction_atom is not None
        or args.co2_target_distance is not None
        or args.co2_target_radius != 1.5
        or args.co2_shell_inner != 4.0
        or args.co2_shell_outer != 6.0
        or args.co2_pack_replicas != 1
        or args.co2_project != ROOT / "co2_screening"
        or args.co2_solute_atoms != 39
        or args.co2_seed_base != 314159
        or args.co2_accommodation_level != "loose"
        or args.co2_accommodation_cycles != 30
        or args.co2_accommodation_engine != "auto"
        or args.co2_md
        or args.co2_extended
        or args.co2_equil_time_ps != 1.0
        or args.co2_screen_time_ps != 5.0
        or args.co2_extended_time_ps != 20.0
        or args.co2_equil_dump_fs != 10.0
        or args.co2_production_dump_fs != 2.0
        or args.co2_parallel_jobs != 1
        or args.co2_start_stage is not None
    ):
        raise SystemExit(
            "CO2-specific source/count/repack options require "
            "--co2-shell-screen."
        )

    if args.run and args.dry_run:
        raise SystemExit("--run cannot be combined with --dry-run.")

    if args.dry_run and args.start_stage != EXTENDED_STAGE_NAME:
        raise SystemExit(
            f"--dry-run requires --start-stage {EXTENDED_STAGE_NAME}."
        )

    if args.resume and args.start_stage is not None:
        raise SystemExit("--resume cannot be combined with --start-stage.")

    if args.resume and args.dry_run:
        raise SystemExit("--resume cannot be combined with --dry-run.")

    if (
        (args.resume or args.start_stage is not None)
        and not (args.run or args.dry_run)
    ):
        raise SystemExit("--resume/--start-stage require --run or --dry-run.")

    if args.resume and args.force:
        raise SystemExit("--resume cannot be combined with --force.")

    if (
        args.start_stage is not None
        and args.force
        and args.start_stage != EXTENDED_STAGE_NAME
    ):
        raise SystemExit(
            "--start-stage can be combined with --force only for the isolated "
            f"{EXTENDED_STAGE_NAME} continuation."
        )

    if args.dry_run and args.force:
        raise SystemExit("--dry-run cannot be combined with --force.")

    if (args.resume or args.start_stage is not None) and args.repack:
        raise SystemExit(
            "--resume/--start-stage cannot be combined with --repack."
        )


def main():
    args = parse_args()
    validate_args(args)

    selected = select_systems(args)
    if args.co2_shell_screen:
        try:
            run_co2_workflow(args, selected[0])
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc
        return

    existing_only = args.resume or args.start_stage is not None

    if args.run and shutil.which(args.xtb) is None:
        raise SystemExit(
            f"xTB executable '{args.xtb}' not found in PATH. "
            "Use --xtb /path/to/xtb if necessary."
        )

    replica_dirs = []

    print(f"Output project: {args.project}")
    print(f"Replicas/system: {args.replicas}")

    if existing_only:
        if not args.project.is_dir():
            raise SystemExit(
                f"Existing project directory not found: {args.project}"
            )
        for system_name in selected:
            for replica_index in range(1, args.replicas + 1):
                replica_dir = (
                    args.project
                    / system_name
                    / f"replica_{replica_index:02d}"
                )
                if not replica_dir.is_dir():
                    raise SystemExit(
                        f"Existing replica directory not found: {replica_dir}"
                    )
                if not (replica_dir / "manifest.json").is_file():
                    raise SystemExit(
                        f"Existing replica has no manifest.json: {replica_dir}"
                    )
                replica_dirs.append(replica_dir)
        print("Existing-only mode: Packmol and preparation are not run.")
    else:
        water_pdb = find_water_pdb(args.water_pdb)
        print(f"Water template: {water_pdb}")

        for name in selected:
            solute = SYSTEMS[name]["solute"]
            if not solute.exists():
                raise SystemExit(
                    f"Unsolvated PDB not found for {name}:\n{solute}"
                )

        if shutil.which(args.packmol) is None:
            raise SystemExit(
                f"Packmol executable '{args.packmol}' not found in PATH. "
                "Use --packmol /path/to/packmol if necessary."
            )

        args.project.mkdir(parents=True, exist_ok=True)

        for system_name in selected:
            solute = SYSTEMS[system_name]["solute"]
            for replica_index in range(1, args.replicas + 1):
                replica_dir = prepare_replica(
                    system_name=system_name,
                    solute_pdb=solute,
                    water_pdb=water_pdb,
                    replica_index=replica_index,
                    project_dir=args.project,
                    args=args,
                )
                replica_dirs.append(replica_dir)

    if args.run or args.dry_run:
        print("\nStarting xTB relaxation + MD E2 pipeline...")

        for replica_dir in replica_dirs:
            run_replica(
                replica_dir,
                args,
            )

        if args.dry_run:
            print(
                "\nDry-run validation completed; no files were written and "
                "xTB was not called."
            )
        else:
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
