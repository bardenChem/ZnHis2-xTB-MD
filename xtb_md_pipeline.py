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
import hashlib
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


def relax_input(n_solute_atoms: int, wall_radius_bohr: float) -> str:
    return f"""$fix
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
        relax_input(n_solute_atoms, wall_radius_bohr)
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


def materialize_relaxed_pdb(replica_dir: Path, expected_atoms: int) -> Path:
    """Find xTB's optimized geometry and create the operational relaxed PDB."""
    template = replica_dir / "system_centered.pdb"
    destination = replica_dir / "system_relaxed.pdb"

    for name in ["xtbopt.pdb", "xtblast.pdb"]:
        candidate = replica_dir / name
        if valid_pdb(candidate, expected_atoms):
            candidate_atoms = pdb_atoms(candidate)
            template_elements = [a["element"] for a in pdb_atoms(template)]
            if [a["element"] for a in candidate_atoms] != template_elements:
                raise RuntimeError(
                    f"{name} element sequence differs from input."
                )
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
        "max_cycles": args.relax_cycles,
        "wall_radius_bohr": manifest["wall"]["radius_bohr"],
        "fixed_atoms": manifest["relaxation"]["fixed_atoms"],
    }


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


def inspect_md_log(log_path: Path):
    text = log_path.read_text(errors="replace")
    lowered = text.lower()

    fatal_patterns = []
    for pattern in FATAL_MD_PATTERNS:
        if pattern.lower() in lowered:
            fatal_patterns.append(pattern)
    return {
        "fatal_patterns": fatal_patterns,
        "thermostating_problem": "thermostating problem" in lowered,
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
    failed = archive / "stage.failed"
    if failed.exists():
        failed.unlink()
    (archive / "stage.done").write_text("ok\n")


def mark_stage_failed(archive: Path, reason: str):
    done = archive / "stage.done"
    if done.exists():
        done.unlink()
    archive.mkdir(parents=True, exist_ok=True)
    (archive / "stage.failed").write_text(reason.rstrip() + "\n")


def md_stage_configuration(
    replica_dir: Path,
    stage: dict,
    geometry_name: str,
    input_name: str,
    input_restart_sha256: str | None,
    args,
) -> dict:
    manifest = json.loads((replica_dir / "manifest.json").read_text())
    return {
        "stage": stage["name"],
        "gfn": args.gfn,
        "charge": args.charge,
        "uhf": args.uhf,
        "alpb": args.alpb,
        "threads": args.threads,
        "temp_K": stage["temp"],
        "time_ps": stage["time"],
        "step_fs": 0.5,
        "dump_fs": 10.0,
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


def write_md_stage_manifest(
    archive: Path,
    configuration: dict,
    output_restart_sha256: str,
):
    archive.mkdir(parents=True, exist_ok=True)
    data = dict(configuration)
    data["output_restart_sha256"] = output_restart_sha256
    (archive / "stage_manifest.json").write_text(
        json.dumps(data, indent=2)
    )


def validate_completed_md_stage(
    archive: Path,
    expected_configuration: dict,
):
    stage_name = expected_configuration["stage"]
    required = [
        "stage.done",
        "mdrestart",
        "xtb.trj",
        "xtbmdok",
        "stage_manifest.json",
    ]
    if (archive / "stage.failed").exists():
        raise RuntimeError(
            f"Stage {stage_name} has both success and failure markers. "
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

    mismatched = [
        key for key, expected in expected_configuration.items()
        if recorded.get(key) != expected
    ]
    if mismatched:
        raise RuntimeError(
            f"Completed stage {stage_name} is incompatible with the current "
            f"MD settings (fields: {', '.join(mismatched)}). "
            "Use --force to rerun it."
        )

    archived_restart_hash = file_sha256(archive / "mdrestart")
    if recorded.get("output_restart_sha256") != archived_restart_hash:
        raise RuntimeError(
            f"Stage {stage_name} has an archived mdrestart inconsistent with "
            "stage_manifest.json. Inspect the archive and rerun with --force."
        )


def validate_output_restart(
    output_restart: Path,
    input_restart_sha256: str | None,
) -> str:
    if not output_restart.exists():
        raise RuntimeError("missing mdrestart")
    output_sha256 = file_sha256(output_restart)
    if (
        input_restart_sha256 is not None
        and output_sha256 == input_restart_sha256
    ):
        raise RuntimeError(
            "output mdrestart is byte-identical to input mdrestart"
        )
    return output_sha256


def prepare_md_stage_attempt(
    replica_dir: Path,
    stage_name: str,
    force: bool,
    expected_configuration: dict,
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
        validate_completed_md_stage(archive, expected_configuration)
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
        if previous.get("configuration") != signature:
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


def run_replica(replica_dir: Path, args):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)
    env.setdefault("MKL_NUM_THREADS", str(args.threads))
    env.setdefault("OMP_STACKSIZE", "4G")

    if not args.skip_relax:
        run_relaxation(replica_dir, args, env)

    for i, stage in enumerate(STAGES):
        name = stage["name"]
        input_name = f"{name}.inp"
        geometry_name = (
            "system_relaxed.pdb"
            if i == 0 and not args.skip_relax
            else "system_centered.pdb"
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
        archive, already_done = prepare_md_stage_attempt(
            replica_dir,
            name,
            args.force,
            stage_configuration,
        )
        if already_done:
            print(
                f"  SKIP {replica_dir.parent.name}/"
                f"{replica_dir.name} {name}: already completed"
            )
            continue

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

        log_path = replica_dir / f"{name}.out"
        cmd = xtb_command(args, geometry_name, input_name)

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

        if validation["thermostating_problem"]:
            archive_failed_md_stage(
                replica_dir,
                name,
                log_path,
                "thermostating problem",
            )
            next_stage = (
                STAGES[i + 1]["name"]
                if i + 1 < len(STAGES)
                else "pipeline completion"
            )
            raise RuntimeError(
                f"{name} completed numerically but xTB reported "
                "'thermostating problem'. Outputs were archived for "
                f"inspection. Refusing to continue to {next_stage}."
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

        restart_path = replica_dir / "mdrestart"
        try:
            output_restart_sha256 = validate_output_restart(
                restart_path,
                input_restart_sha256,
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
                f"{name} ended without mdrestart for "
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

        write_md_stage_manifest(
            archive,
            stage_configuration,
            output_restart_sha256,
        )
        archive_stage_outputs(
            replica_dir,
            name,
            log_path,
        )

        mark_stage_done(archive)
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

    if args.relax_cycles < 1:
        raise SystemExit("--relax-cycles must be >= 1.")


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
        print("\nStarting xTB relaxation + MD E2 pipeline...")

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
