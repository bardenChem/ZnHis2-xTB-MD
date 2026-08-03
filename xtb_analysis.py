#!/usr/bin/env python3
"""xtb_analysis.py

Post-processing, plotting, and coordinate-export toolkit for xTB
molecular-dynamics trajectories.

Focus of v0.3:
- xTB extended-XYZ trajectories (xtb.trj), with coordinates only or with an
  additional per-frame velocity block
- stage archives from xtb_md_pipeline.py
- finite-droplet structural screening
- reusable CSV/JSON outputs for future OOCCuPy integration

Implemented analyses:
1. xTB trajectory parsing and provenance.
2. Zn--named-site distance time series.
3. Nearest Zn--water-O distance.
4. Zn-centered water-O shell counts and cumulative N(r).
5. Optional hard and smooth coordination numbers.
6. Optional tetrahedral q parameter for four nearest donor candidates.
7. Optional Zn--water contact episodes with hysteresis.
8. Kabsch-aligned solute RMSD and RMSF.
9. Standard XYZ export for TRAVIS.
10. Optional headless diagnostic plots from the calculated analysis tables.
11. Validated coordinate-only exports for TRAVIS and VMD workflows.

No bulk-normalized RDF is computed for finite droplets because a spherical,
inhomogeneous droplet with a wall is not equivalent to periodic bulk liquid.

Velocity records are preserved in ``Frame.velocities`` but are not analyzed in
v0.3, no physical unit is assigned to them, and the TRAVIS/VMD exports continue
to contain coordinates only.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from typing import Iterator

import numpy as np

DEFAULT_STAGE = "04_298K_screen"
DEFAULT_WATER_RESNAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "H2O"}
FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"


@dataclass(frozen=True)
class AtomInfo:
    index1: int
    element: str
    name: str | None = None
    resname: str | None = None
    resid: str | None = None


@dataclass
class Frame:
    index: int
    comment: str
    elements: list[str]
    xyz: np.ndarray
    metadata: dict
    velocities: np.ndarray | None = None


@dataclass(frozen=True)
class Site:
    label: str
    index0: int
    group: str


@dataclass(frozen=True)
class TrajectorySpec:
    path: Path
    stage: str
    dt_fs: float | None


@dataclass(frozen=True)
class KabschTransform:
    rotation_matrix: np.ndarray
    translation_before_rotation: np.ndarray
    translation_after_rotation: np.ndarray


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def element_name(value: str) -> str:
    s = re.sub(r"[^A-Za-z]", "", value)
    if s[:2].lower() == "zn":
        return "Zn"
    return s[:1].upper() + s[1:].lower() if s else value


def parse_comment(text: str) -> dict:
    out = {"comment": text.strip()}
    for key in ("energy", "gnorm"):
        m = re.search(rf"\b{key}\s*:\s*({FLOAT})", text, re.I)
        if m:
            out[key] = float(m.group(1).replace("D", "E").replace("d", "e"))
    m = re.search(r"\bxtb\s*:\s*(.+?)\s*$", text, re.I)
    if m:
        out["xtb"] = m.group(1).strip()
    return out


def parse_float_token(token: str) -> float:
    """Parse one finite float, including Fortran D/d exponents."""
    value = float(token.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError(f"non-finite numeric value: {token!r}")
    return value


def parse_numeric_triplet(line: str) -> np.ndarray | None:
    """Return exactly three finite numbers, or None for any other record."""
    fields = line.split()
    if len(fields) != 3:
        return None
    try:
        values = [parse_float_token(token) for token in fields]
    except ValueError:
        return None
    return np.asarray(values, dtype=float)


def read_next_nonblank_line(handle) -> str | None:
    """Read the next nonblank line, returning None only at end of file."""
    while True:
        line = handle.readline()
        if not line:
            return None
        if line.strip():
            return line


def velocity_layout(frames_total: int, frames_with_velocities: int) -> str:
    """Classify velocity presence without assigning units or interpretation."""
    if frames_with_velocities == 0:
        return "coordinates_only"
    if frames_with_velocities == frames_total:
        return "coordinates_and_velocities"
    return "mixed"


def iter_xtb_trj(path: Path) -> Iterator[Frame]:
    """Stream xTB extended XYZ with optional per-frame velocity blocks."""
    with path.open("r", errors="replace") as f:
        iframe = 0
        pending_line = None
        while True:
            line = pending_line
            pending_line = None
            if line is None:
                line = read_next_nonblank_line(f)
            if line is None:
                return
            try:
                nat = int(line.strip())
            except ValueError as exc:
                content = line.rstrip("\r\n")
                raise RuntimeError(
                    f"{path}: expected atom count at frame {iframe}, "
                    f"found {content!r}"
                ) from exc
            if nat < 1:
                raise RuntimeError(
                    f"{path}: atom count must be positive at frame {iframe}, found {nat}"
                )
            comment = f.readline()
            if not comment:
                raise RuntimeError(f"{path}: missing comment at frame {iframe}")
            xyz = np.empty((nat, 3), float)
            elements = []
            for i in range(nat):
                atom_line = f.readline()
                fields = atom_line.split()
                if len(fields) < 4:
                    content = atom_line.rstrip("\r\n")
                    raise RuntimeError(
                        f"{path}: malformed atom record at frame {iframe}, "
                        f"atom {i+1}/{nat}: {content!r}"
                    )
                elements.append(element_name(fields[0]))
                try:
                    xyz[i] = [parse_float_token(token) for token in fields[1:4]]
                except ValueError as exc:
                    content = atom_line.rstrip("\r\n")
                    raise RuntimeError(
                        f"{path}: malformed atom record at frame {iframe}, "
                        f"atom {i+1}/{nat}: {content!r}"
                    ) from exc

            velocities = None
            next_line = read_next_nonblank_line(f)
            if next_line is not None:
                first_velocity = parse_numeric_triplet(next_line)
                if first_velocity is not None:
                    velocities = np.empty((nat, 3), dtype=float)
                    velocities[0] = first_velocity
                    for i in range(1, nat):
                        velocity_line = f.readline()
                        if not velocity_line:
                            raise RuntimeError(
                                f"{path}: truncated velocity block at frame {iframe}, "
                                f"found {i} of {nat} records"
                            )
                        velocity = parse_numeric_triplet(velocity_line)
                        if velocity is None:
                            content = velocity_line.rstrip("\r\n")
                            raise RuntimeError(
                                f"{path}: malformed velocity record at frame {iframe}, "
                                f"atom {i+1}/{nat}: {content!r}"
                            )
                        velocities[i] = velocity
                else:
                    if len(next_line.split()) == 3:
                        content = next_line.rstrip("\r\n")
                        raise RuntimeError(
                            f"{path}: malformed velocity record at frame {iframe}, "
                            f"atom 1/{nat}: {content!r}"
                        )
                    pending_line = next_line

            yield Frame(
                iframe,
                comment.rstrip(),
                elements,
                xyz,
                parse_comment(comment),
                velocities,
            )
            iframe += 1


def parse_pdb(path: Path) -> list[AtomInfo]:
    atoms = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[12:16].strip() or None
        resname = line[17:20].strip() or None
        resid = line[22:26].strip() or None
        elem = line[76:78].strip() if len(line) >= 78 else ""
        if not elem:
            elem = name or ""
        atoms.append(AtomInfo(len(atoms) + 1, element_name(elem), name, resname, resid))
    if not atoms:
        raise RuntimeError(f"No PDB atoms found in {path}")
    return atoms


def parse_site(text: str) -> Site:
    fields = text.split(":")
    if len(fields) not in (2, 3):
        raise argparse.ArgumentTypeError("--site must be LABEL:INDEX[:GROUP]")
    try:
        index1 = int(fields[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("site index must be integer") from exc
    if index1 < 1:
        raise argparse.ArgumentTypeError("site index is one-based and must be >= 1")
    return Site(fields[0], index1 - 1, fields[2] if len(fields) == 3 else "site")


def parse_alignment_indices(text: str) -> list[int]:
    """Parse one-based comma/range syntax and return zero-based indices."""
    indices = []
    seen = set()
    for raw_field in text.split(","):
        field = raw_field.strip()
        if not field:
            raise argparse.ArgumentTypeError(
                "--alignment-indices contains an empty field"
            )
        match = re.fullmatch(r"(\d+)(?:-(\d+))?", field)
        if not match:
            raise argparse.ArgumentTypeError(
                "--alignment-indices must use syntax such as 1-5,10,12-19"
            )
        start = int(match.group(1))
        end = int(match.group(2) or start)
        if start < 1 or end < start:
            raise argparse.ArgumentTypeError(
                "alignment indices are one-based and ranges must be ascending"
            )
        for index1 in range(start, end + 1):
            index0 = index1 - 1
            if index0 in seen:
                raise argparse.ArgumentTypeError(
                    f"duplicate alignment index: {index1}"
                )
            seen.add(index0)
            indices.append(index0)
    return indices


def manifest(replica: Path) -> dict:
    try:
        return json.loads((replica / "manifest.json").read_text())
    except Exception as exc:
        raise RuntimeError(f"Cannot read {replica/'manifest.json'}") from exc


def topology_from_replica(replica: Path) -> Path | None:
    for p in (replica / "system_centered.pdb", replica / "system_relaxed.pdb"):
        if p.is_file():
            return p
    return None


def trajectory_specs(args) -> tuple[list[TrajectorySpec], Path | None, int | None]:
    if args.replica:
        rep = args.replica.resolve()
        m = manifest(rep)
        dt = args.frame_dt_fs
        if dt is None:
            try:
                dt = float(m["md"]["dump_fs"])
            except Exception:
                pass
        stages = args.stage or [DEFAULT_STAGE]
        specs = []
        for stage in stages:
            p = rep / "stages" / stage / "xtb.trj"
            if not p.is_file():
                if len(stages) == 1 and (rep / "xtb.trj").is_file():
                    p = rep / "xtb.trj"
                else:
                    raise RuntimeError(f"Missing trajectory: {p}")
            specs.append(TrajectorySpec(p, stage, dt))
        top = args.topology.resolve() if args.topology else topology_from_replica(rep)
        nsolute = args.n_solute
        if nsolute is None:
            try:
                nsolute = int(m["relaxation"]["n_solute_atoms"])
            except Exception:
                pass
        return specs, top, nsolute

    specs = [TrajectorySpec(p.resolve(), p.stem, args.frame_dt_fs) for p in args.trajectory]
    return specs, args.topology.resolve() if args.topology else None, args.n_solute


def unique_zn(elements: list[str], user_index: int | None) -> int:
    if user_index is not None:
        return user_index - 1
    idx = [i for i, e in enumerate(elements) if e.lower() == "zn"]
    if len(idx) != 1:
        raise RuntimeError(f"Expected exactly one Zn, found {len(idx)}; use --zn-index")
    return idx[0]


def water_oxygens(topology, elements, nsolute, water_resnames):
    warnings = []
    if topology and len(topology) == len(elements):
        idx = [a.index1 - 1 for a in topology if (a.resname or "").upper() in water_resnames and a.element == "O"]
        if idx:
            return idx, warnings
    if nsolute is not None:
        idx = []
        i = nsolute
        while i + 2 < len(elements):
            if [e.upper() for e in elements[i:i+3]] == ["O", "H", "H"]:
                idx.append(i)
                i += 3
            else:
                i += 1
        if idx:
            warnings.append("Water O atoms inferred from O-H-H triplets after the solute; verify topology/order.")
            return idx, warnings
    warnings.append("Water O atoms were not identified; water-specific analyses skipped.")
    return [], warnings


def kabsch_transform(mobile, reference) -> KabschTransform:
    """Return the rigid transform that aligns mobile row vectors to reference."""
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Kabsch inputs must have matching shape (n, 3)")
    mc = mobile.mean(0)
    rc = reference.mean(0)
    mob = mobile - mc
    ref = reference - rc
    u, _, vt = np.linalg.svd(mob.T @ ref)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        vt[-1] *= -1
        rot = u @ vt
    return KabschTransform(
        rotation_matrix=rot,
        translation_before_rotation=-mc,
        translation_after_rotation=rc,
    )


def apply_kabsch_transform(xyz, transform: KabschTransform):
    xyz = np.asarray(xyz, dtype=float)
    return (
        (xyz + transform.translation_before_rotation)
        @ transform.rotation_matrix
        + transform.translation_after_rotation
    )


def kabsch(mobile, reference):
    """Align coordinates while preserving the original public helper."""
    transform = kabsch_transform(mobile, reference)
    return apply_kabsch_transform(mobile, transform)


def validate_alignment_geometry(xyz, indices):
    if len(indices) < 3:
        raise RuntimeError("Alignment requires at least three atoms")
    selected = np.asarray(xyz, dtype=float)[indices]
    centered = selected - selected.mean(axis=0)
    if np.linalg.matrix_rank(centered) < 2:
        raise RuntimeError(
            "Alignment selection must contain at least three non-collinear atoms"
        )


def resolve_alignment_selection(
    elements,
    xyz,
    nsolute,
    selection,
    explicit_indices,
):
    """Return alignment indices and reproducible selection metadata."""
    nat = len(elements)
    if explicit_indices is not None:
        indices = list(explicit_indices)
        mode = "explicit-indices"
    elif nsolute is None:
        return None, {
            "mode": selection,
            "indices_one_based": [],
            "available": False,
            "reason": "n_solute is unavailable and no explicit indices were supplied",
        }
    elif selection == "solute-heavy":
        indices = [i for i in range(nsolute) if elements[i].upper() != "H"]
        mode = selection
    elif selection == "solute-all":
        indices = list(range(nsolute))
        mode = selection
    else:
        raise RuntimeError(f"Unsupported alignment selection: {selection}")

    if len(indices) != len(set(indices)):
        raise RuntimeError("Alignment selection contains duplicate atom indices")
    if any(index < 0 or index >= nat for index in indices):
        bad = [index + 1 for index in indices if index < 0 or index >= nat]
        raise RuntimeError(f"Alignment indices out of range for {nat} atoms: {bad}")
    validate_alignment_geometry(xyz, indices)
    return indices, {
        "mode": mode,
        "indices_one_based": [index + 1 for index in indices],
        "elements": [elements[index] for index in indices],
        "available": True,
        "reason": None,
    }


def tetra_q(vectors):
    if vectors.shape != (4, 3):
        return math.nan
    norm = np.linalg.norm(vectors, axis=1)
    if np.any(norm == 0):
        return math.nan
    u = vectors / norm[:, None]
    s = 0.0
    for j in range(3):
        for k in range(j + 1, 4):
            c = float(np.clip(u[j] @ u[k], -1.0, 1.0))
            s += (c + 1.0 / 3.0) ** 2
    return 1.0 - 3.0 * s / 8.0


def smooth_cn(distances, r0, n):
    distances = np.asarray(distances)
    return float(np.sum(1.0 / (1.0 + (distances / r0) ** n)))


def write_csv(path: Path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if fields:
            w.writeheader()
            w.writerows(rows)


def write_xyz_frame(handle, elements, xyz, comment):
    handle.write(f"{len(elements)}\n{comment}\n")
    for e, r in zip(elements, xyz):
        handle.write(f"{e:<3s} {r[0]: .12f} {r[1]: .12f} {r[2]: .12f}\n")


def export_comment(
    system,
    stage,
    source_frame,
    global_frame,
    time_ps,
    alignment,
):
    """Build a deterministic one-line XYZ comment."""
    clean = lambda value: re.sub(r"\s+", "_", str(value))
    fields = [
        f"system={clean(system)}",
        f"stage={clean(stage)}",
        f"source_frame={source_frame}",
        f"global_frame={global_frame}",
    ]
    if time_ps is not None and math.isfinite(time_ps):
        fields.append(f"time_ps={time_ps:.6f}")
    fields.append(f"alignment={clean(alignment)}")
    return " ".join(fields)


class CoordinateExporter:
    """Stream coordinate-only XYZ variants from already-retained frames."""

    def __init__(
        self,
        output_dir,
        nsolute,
        alignment_indices,
        alignment_mode,
        system,
    ):
        self.output_dir = output_dir
        self.nsolute = nsolute
        self.alignment_indices = alignment_indices
        self.alignment_mode = alignment_mode
        self.system = system
        self.reference_xyz = None
        self.reference_frame = None
        self.exported_frames = 0
        self.first_exported_time_ps = None
        self.last_exported_time_ps = None
        self.paths = {
            "full_coordinates": output_dir / "trajectory_full_coordinates.xyz",
        }
        if nsolute is not None:
            self.paths["solute_coordinates"] = (
                output_dir / "trajectory_solute_coordinates.xyz"
            )
        if alignment_indices is not None:
            self.paths["full_aligned_on_solute"] = (
                output_dir / "trajectory_full_aligned_on_solute.xyz"
            )
            if nsolute is not None:
                self.paths["solute_aligned"] = (
                    output_dir / "trajectory_solute_aligned.xyz"
                )

        managed_variants = {
            "trajectory_full_coordinates.xyz",
            "trajectory_full_aligned_on_solute.xyz",
            "trajectory_solute_coordinates.xyz",
            "trajectory_solute_aligned.xyz",
        }
        expected_names = {path.name for path in self.paths.values()}
        stale = sorted(
            str(output_dir / name)
            for name in managed_variants - expected_names
            if (output_dir / name).exists()
        )
        vmd_dir = output_dir.parent / "vmd"
        if alignment_indices is None and vmd_dir.exists():
            stale.extend(str(path) for path in sorted(vmd_dir.iterdir()))
        if stale:
            raise RuntimeError(
                "Existing export files are incompatible with the current "
                f"selection; use a clean --output-dir or remove them explicitly: {stale}"
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        self.handles = {}
        try:
            for name, path in self.paths.items():
                self.handles[name] = path.open("w")
        except Exception:
            self.close()
            raise

    def set_reference(self, xyz, context):
        if self.alignment_indices is None or self.reference_xyz is not None:
            return
        self.reference_xyz = np.asarray(xyz, dtype=float).copy()
        self.reference_frame = dict(context)

    def write(
        self,
        elements,
        xyz,
        stage,
        source_frame,
        global_frame,
        time_ps,
    ):
        unaligned_comment = export_comment(
            self.system,
            stage,
            source_frame,
            global_frame,
            time_ps,
            "none",
        )
        write_xyz_frame(
            self.handles["full_coordinates"],
            elements,
            xyz,
            unaligned_comment,
        )
        if "solute_coordinates" in self.handles:
            write_xyz_frame(
                self.handles["solute_coordinates"],
                elements[:self.nsolute],
                xyz[:self.nsolute],
                unaligned_comment,
            )

        if self.alignment_indices is not None:
            if self.reference_xyz is None:
                raise RuntimeError("Alignment reference frame was not initialized")
            transform = kabsch_transform(
                xyz[self.alignment_indices],
                self.reference_xyz[self.alignment_indices],
            )
            aligned = apply_kabsch_transform(xyz, transform)
            aligned_comment = export_comment(
                self.system,
                stage,
                source_frame,
                global_frame,
                time_ps,
                self.alignment_mode,
            )
            write_xyz_frame(
                self.handles["full_aligned_on_solute"],
                elements,
                aligned,
                aligned_comment,
            )
            if "solute_aligned" in self.handles:
                write_xyz_frame(
                    self.handles["solute_aligned"],
                    elements[:self.nsolute],
                    aligned[:self.nsolute],
                    aligned_comment,
                )

        self.exported_frames += 1
        if time_ps is not None and math.isfinite(time_ps):
            if self.first_exported_time_ps is None:
                self.first_exported_time_ps = time_ps
            self.last_exported_time_ps = time_ps

    def close(self):
        for handle in getattr(self, "handles", {}).values():
            if not handle.closed:
                handle.close()


def validate_coordinate_xyz(path, expected_frames, expected_atoms):
    """Validate a coordinate-only, constant-order XYZ trajectory."""
    frame_count = 0
    reference_elements = None
    pending_line = None
    with path.open("r", errors="replace") as handle:
        while True:
            atom_count_line = pending_line
            pending_line = None
            if atom_count_line is None:
                atom_count_line = read_next_nonblank_line(handle)
            if atom_count_line is None:
                break
            try:
                nat = int(atom_count_line.strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}: expected atom count while validating frame "
                    f"{frame_count}, found {atom_count_line.rstrip()!r}"
                ) from exc
            if nat != expected_atoms:
                raise RuntimeError(
                    f"{path}: frame {frame_count} has {nat} atoms; "
                    f"expected {expected_atoms}"
                )
            comment = handle.readline()
            if not comment:
                raise RuntimeError(
                    f"{path}: missing comment at exported frame {frame_count}"
                )
            elements = []
            for iatom in range(nat):
                coordinate_line = handle.readline()
                fields = coordinate_line.split()
                if len(fields) != 4:
                    raise RuntimeError(
                        f"{path}: exported coordinate record must contain exactly "
                        f"element x y z at frame {frame_count}, atom "
                        f"{iatom+1}/{nat}: {coordinate_line.rstrip()!r}"
                    )
                try:
                    [parse_float_token(token) for token in fields[1:4]]
                except ValueError as exc:
                    raise RuntimeError(
                        f"{path}: non-finite or malformed coordinate at frame "
                        f"{frame_count}, atom {iatom+1}/{nat}: "
                        f"{coordinate_line.rstrip()!r}"
                    ) from exc
                elements.append(element_name(fields[0]))
            if reference_elements is None:
                reference_elements = elements
            elif elements != reference_elements:
                raise RuntimeError(
                    f"{path}: element order changed at frame {frame_count}"
                )
            frame_count += 1

            next_line = read_next_nonblank_line(handle)
            if next_line is not None:
                if parse_numeric_triplet(next_line) is not None or len(
                    next_line.split()
                ) == 3:
                    raise RuntimeError(
                        f"{path}: exported XYZ unexpectedly contains velocity records"
                    )
                pending_line = next_line

    if frame_count != expected_frames:
        raise RuntimeError(
            f"{path}: found {frame_count} frames; expected {expected_frames}"
        )
    return {
        "frames": frame_count,
        "atoms_per_frame": expected_atoms,
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
        "coordinates_only": True,
    }


def write_minimal_pdb(path, elements, xyz):
    """Write order/element-preserving PDB records without inferred bonds."""
    with path.open("w") as handle:
        for i, (element, coord) in enumerate(zip(elements, xyz), 1):
            if not np.all(np.isfinite(coord)):
                raise RuntimeError(f"Cannot write non-finite PDB coordinate at atom {i}")
            atom_name = element.upper()[:2]
            handle.write(
                f"HETATM{i:5d} {atom_name:^4s} SYS A{1:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                f"  1.00  0.00          {element:>2s}\n"
            )
        handle.write("END\n")


def prepare_reference_topology(output_path, topology_path, first_frame):
    if topology_path is not None:
        source = topology_path.resolve()
        if source != output_path.resolve():
            shutil.copy2(source, output_path)
        provenance = {
            "status": "copied",
            "source_path": str(source),
            "source_sha256": sha256(source),
        }
    else:
        write_minimal_pdb(output_path, first_frame.elements, first_frame.xyz)
        provenance = {
            "status": "generated",
            "source_path": None,
            "source_sha256": None,
            "note": "Minimal PDB generated without inferred connectivity.",
        }
    exported_atoms = parse_pdb(output_path)
    if len(exported_atoms) != len(first_frame.elements):
        raise RuntimeError(
            f"{output_path}: topology has {len(exported_atoms)} atoms; "
            f"expected {len(first_frame.elements)}"
        )
    if [atom.element for atom in exported_atoms] != first_frame.elements:
        raise RuntimeError(
            f"{output_path}: topology element/order differs from trajectory"
        )
    provenance.update({
        "output_path": str(output_path),
        "output_sha256": sha256(output_path),
    })
    return provenance


def write_atom_map(path, elements, topology, nsolute, water_O, zn, sites):
    site_map = {}
    for site in sites:
        entry = site_map.setdefault(site.index0, {"labels": [], "groups": []})
        entry["labels"].append(site.label)
        entry["groups"].append(site.group)
    water_set = set(water_O)
    rows = []
    for i, element in enumerate(elements):
        atom = topology[i] if topology else None
        site = site_map.get(i, {"labels": [], "groups": []})
        rows.append({
            "index": i + 1,
            "element": element,
            "atom_name": atom.name if atom else "",
            "resname": atom.resname if atom else "",
            "resid": atom.resid if atom else "",
            "is_solute": bool(nsolute is not None and i < nsolute),
            "is_water_O": i in water_set,
            "is_Zn": i == zn,
            "site_label": ";".join(site["labels"]),
            "site_group": ";".join(site["groups"]),
        })
    fields = [
        "index", "element", "atom_name", "resname", "resid", "is_solute",
        "is_water_O", "is_Zn", "site_label", "site_group",
    ]
    write_csv(path, rows, fields)


def write_travis_readme(path, system, available_files):
    system_dir = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(system)) or "system"
    available = set(available_files)
    lines = [
        "# TRAVIS coordinate exports",
        "",
        "All XYZ files contain coordinates only. Velocities from `xtb.trj` were",
        "intentionally omitted. Run each TRAVIS analysis in a separate directory:",
        "",
        "```text",
        f"travis_runs/{system_dir}/",
        "├── tdo/",
        "├── adf/",
        "├── sdf/",
        "├── cdf/",
        "└── radial/",
        "```",
        "",
        "Do not combine different analyses in the same output directory, and do",
        "not automate menu numbers because they can change between TRAVIS versions.",
        "",
    ]
    if "solute_aligned" in available:
        lines += [
            "For TDO and solute-geometry inspection:",
            "",
            "```bash",
            f"mkdir -p travis_runs/{system_dir}/tdo",
            f"cd travis_runs/{system_dir}/tdo",
            f"travis -p ../../../analysis/{system_dir}/travis/trajectory_solute_aligned.xyz",
            "```",
            "",
        ]
    if "full_aligned_on_solute" in available:
        lines += [
            "For water-relative SDF exploration and aligned videos:",
            "",
            "```bash",
            f"mkdir -p travis_runs/{system_dir}/sdf",
            f"cd travis_runs/{system_dir}/sdf",
            f"travis -p ../../../analysis/{system_dir}/travis/trajectory_full_aligned_on_solute.xyz",
            "```",
            "",
        ]
    lines += [
        "## Files",
        "",
        "- `trajectory_solute_aligned.xyz`: solute-only coordinates aligned on the",
        "  configured atom selection; intended for TDO and complex geometry.",
        "- `trajectory_full_aligned_on_solute.xyz`: the same solute-derived rigid",
        "  transform applied to every atom; intended for solvent-relative views.",
        "- `trajectory_full_coordinates.xyz`: retained frames without alignment.",
        "- `trajectory_solute_coordinates.xyz`: unaligned solute-only frames.",
        "- `reference_topology.pdb`: atom names/residue identity when available; no",
        "  Zn-ligand connectivity is inferred.",
        "- `atom_map.csv`: one-based mapping for TRAVIS and VMD selections.",
        "",
        "Files dependent on `n_solute` or a valid alignment selection are omitted",
        "when that information is unavailable. Adjust the example analysis path if",
        "a non-default `--output-dir` was used.",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_vmd_support(vmd_dir):
    vmd_dir.mkdir(parents=True, exist_ok=True)
    script_path = vmd_dir / "convert_xyz_to_dcd.tcl"
    script_path.write_text(
        "set script_dir [file dirname [file normalize [info script]]]\n"
        "set topology [file join $script_dir .. travis reference_topology.pdb]\n"
        "set trajectory [file join $script_dir .. travis trajectory_full_aligned_on_solute.xyz]\n"
        "set output [file join $script_dir trajectory_full_aligned_on_solute.dcd]\n"
        "set molid [mol new $topology type pdb waitfor all]\n"
        "mol addfile $trajectory type xyz waitfor all molid $molid\n"
        "animate delete beg 0 end 0 $molid\n"
        "animate write dcd $output beg 0 end -1 skip 1 waitfor all $molid\n"
        "quit\n"
    )
    readme_path = vmd_dir / "README_VMD.md"
    readme_path.write_text(
        "# VMD conversion helper\n\n"
        "Run manually from the repository or analysis directory:\n\n"
        "```bash\n"
        "vmd -dispdev text -e vmd/convert_xyz_to_dcd.tcl\n"
        "```\n\n"
        "The helper loads the PDB for topology, adds the aligned coordinate-only\n"
        "XYZ, removes the initial PDB coordinate frame, and writes\n"
        "`trajectory_full_aligned_on_solute.dcd`. The DCD contains coordinates\n"
        "only and must be loaded together with `reference_topology.pdb`. It is\n"
        "intended primarily for VMD, not TRAVIS. The helper is never executed\n"
        "automatically.\n"
    )
    return [script_path, readme_path]


def load_pyplot():
    """Load matplotlib lazily so --no-plots does not require it."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib. Install it with: "
            "python -m pip install matplotlib"
        ) from exc
    return plt


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return cleaned.strip("._") or "unnamed"


def numeric_pairs(rows, xkey, ykey):
    pairs = []
    for row in rows:
        try:
            x = float(row[xkey])
            y = float(row[ykey])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(x) and math.isfinite(y):
            pairs.append((x, y))
    if not pairs:
        return np.array([], float), np.array([], float)
    values = np.asarray(pairs, float)
    return values[:, 0], values[:, 1]


def numeric_values(rows, key):
    values = []
    for row in rows:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return np.asarray(values, float)


def save_figure(plt, fig, path: Path, dpi: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.tight_layout()
        kwargs = {"dpi": dpi} if path.suffix.lower() == ".png" else {}
        fig.savefig(path, **kwargs)
    finally:
        plt.close(fig)
    return path


def style_axis(ax):
    ax.grid(True, alpha=0.25, linewidth=0.7)


def add_cutoff(ax, cutoff_A):
    if cutoff_A is not None:
        ax.axhline(
            cutoff_A,
            color="0.35",
            linestyle="--",
            linewidth=1.0,
            label="coordination cutoff",
        )


def plot_distance_timeseries(
    plt, rows, sites, output_dir, suffix, dpi, cutoff_A,
    time_key, time_label,
):
    available = []
    for site in sites:
        _, y = numeric_pairs(rows, time_key, f"d_Zn_{site.label}_A")
        if y.size:
            available.append(site)
    if not available:
        return []

    generated = []

    def plot_sites(selected, filename, title):
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        for site in selected:
            x, y = numeric_pairs(rows, time_key, f"d_Zn_{site.label}_A")
            ax.plot(x, y, linewidth=1.1, label=site.label)
        add_cutoff(ax, cutoff_A)
        ax.set(xlabel=time_label, ylabel="Zn–site distance / Å", title=title)
        style_axis(ax)
        ax.legend()
        generated.append(save_figure(
            plt, fig, output_dir / f"{filename}.{suffix}", dpi
        ))

    plot_sites(
        available,
        "coordination_distances_vs_time",
        "Zn–site coordination distances",
    )
    groups = {}
    for site in available:
        groups.setdefault(site.group, []).append(site)
    if len(groups) > 1:
        for group in sorted(groups):
            plot_sites(
                groups[group],
                f"Zn_{safe_filename(group)}_distances_vs_time",
                f"Zn–site distances: {group}",
            )
    return generated


def plot_distance_histograms(plt, rows, sites, output_dir, suffix, dpi):
    generated = []
    hist_dir = output_dir / "distributions"
    for site in sites:
        values = numeric_values(rows, f"d_Zn_{site.label}_A")
        if not values.size:
            continue
        fig, ax = plt.subplots(figsize=(6.4, 4.5))
        ax.hist(values, bins="auto", edgecolor="black", linewidth=0.5)
        ax.set(
            xlabel=f"Zn–{site.label} distance / Å",
            ylabel="Frame count",
            title=f"Zn–{site.label} distance distribution (descriptive)",
        )
        style_axis(ax)
        generated.append(save_figure(
            plt,
            fig,
            hist_dir / f"Zn_{safe_filename(site.label)}_distance_histogram.{suffix}",
            dpi,
        ))
    return generated


def plot_nearest_water(
    plt, rows, output_dir, suffix, dpi, cutoff_A, time_key, time_label,
):
    x, y = numeric_pairs(rows, time_key, "nearest_water_O_distance_A")
    if not y.size:
        return []
    generated = []
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(x, y, linewidth=1.1)
    add_cutoff(ax, cutoff_A)
    ax.set(
        xlabel=time_label,
        ylabel="Nearest Zn–Ow distance / Å",
        title="Nearest water oxygen to Zn",
    )
    style_axis(ax)
    if cutoff_A is not None:
        ax.legend()
    generated.append(save_figure(
        plt, fig, output_dir / f"nearest_water_vs_time.{suffix}", dpi
    ))

    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.hist(y, bins="auto", edgecolor="black", linewidth=0.5)
    ax.set(
        xlabel="Nearest Zn–Ow distance / Å",
        ylabel="Frame count",
        title="Nearest-water distance distribution (descriptive)",
    )
    style_axis(ax)
    generated.append(save_figure(
        plt,
        fig,
        output_dir / "distributions" / f"nearest_water_distance_histogram.{suffix}",
        dpi,
    ))
    return generated


def plot_coordination_number(
    plt, rows, output_dir, suffix, dpi, time_key, time_label,
):
    generated = []
    hard_keys = [
        ("CN_total_hard", "total"),
        ("CN_named_hard", "named sites"),
        ("CN_water_hard", "water O"),
    ]
    present = []
    for key, label in hard_keys:
        x, y = numeric_pairs(rows, time_key, key)
        if y.size:
            present.append((key, label, x, y))
    if present:
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        for _, label, x, y in present:
            ax.step(x, y, where="post", linewidth=1.1, label=label)
        ax.set(
            xlabel=time_label,
            ylabel="Coordination number",
            title="Hard-cutoff coordination numbers",
        )
        style_axis(ax)
        ax.legend()
        generated.append(save_figure(
            plt, fig, output_dir / f"coordination_number_vs_time.{suffix}", dpi
        ))

    x, smooth = numeric_pairs(rows, time_key, "CN_smooth")
    if smooth.size:
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot(x, smooth, linewidth=1.1)
        ax.set(
            xlabel=time_label,
            ylabel="Smooth coordination number",
            title="Smooth coordination number",
        )
        style_axis(ax)
        generated.append(save_figure(
            plt,
            fig,
            output_dir / f"smooth_coordination_number_vs_time.{suffix}",
            dpi,
        ))

        fig, ax = plt.subplots(figsize=(6.4, 4.5))
        ax.hist(smooth, bins="auto", edgecolor="black", linewidth=0.5)
        ax.set(
            xlabel="Smooth coordination number",
            ylabel="Frame count",
            title="Smooth-CN distribution (descriptive)",
        )
        style_axis(ax)
        generated.append(save_figure(
            plt,
            fig,
            output_dir / "distributions" / f"CN_smooth_histogram.{suffix}",
            dpi,
        ))
    return generated


def plot_coordination_states(
    plt, rows, state_rows, output_dir, suffix, dpi, time_key, time_label,
):
    timeline = []
    for row in rows:
        state = row.get("coordination_state")
        try:
            time = float(row[time_key])
        except (KeyError, TypeError, ValueError):
            continue
        if state and math.isfinite(time):
            timeline.append((time, str(state)))
    if not timeline:
        return []

    generated = []
    states = sorted({state for _, state in timeline})
    state_to_y = {state: i for i, state in enumerate(states)}
    x = np.asarray([item[0] for item in timeline])
    y = np.asarray([state_to_y[item[1]] for item in timeline])
    fig, ax = plt.subplots(figsize=(9.0, max(4.5, 0.38 * len(states) + 2.0)))
    ax.step(x, y, where="post", linewidth=1.0)
    ax.scatter(x, y, s=8)
    ax.set_yticks(range(len(states)), labels=states)
    ax.set(
        xlabel=time_label,
        ylabel="Coordination state (categorical)",
        title="Coordination state vs time",
    )
    style_axis(ax)
    generated.append(save_figure(
        plt, fig, output_dir / f"coordination_state_vs_time.{suffix}", dpi
    ))

    fractions = []
    for row in state_rows:
        try:
            value = float(row["fraction_of_analyzed_frames"])
        except (KeyError, TypeError, ValueError):
            continue
        fractions.append((str(row.get("state", "")), value))
    if fractions:
        fractions.sort(key=lambda item: item[0])
        labels = [item[0] for item in fractions]
        values = [item[1] for item in fractions]
        fig, ax = plt.subplots(figsize=(8.5, max(4.5, 0.38 * len(labels) + 2.0)))
        ax.barh(range(len(labels)), values)
        ax.set_yticks(range(len(labels)), labels=labels)
        ax.set(
            xlabel="Fraction of analyzed frames",
            ylabel="Coordination state",
            title="Descriptive state fractions (not equilibrium populations)",
        )
        style_axis(ax)
        generated.append(save_figure(
            plt,
            fig,
            output_dir / f"coordination_state_fraction.{suffix}",
            dpi,
        ))
    return generated


def plot_radial_number(plt, radial_rows, output_dir, suffix, dpi):
    r = numeric_values(radial_rows, "r_mid_A")
    shell = numeric_values(radial_rows, "mean_water_O_count_in_shell")
    cumulative = numeric_values(radial_rows, "cumulative_N_water_O")
    if not r.size or len(shell) != len(r) or len(cumulative) != len(r):
        return []
    generated = []
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.plot(r, shell, linewidth=1.2)
    ax.set(
        xlabel="r / Å",
        ylabel="Mean water O count in shell",
        title="Zn–Ow radial number distribution",
    )
    style_axis(ax)
    generated.append(save_figure(
        plt, fig, output_dir / f"Zn_Owater_shell_count.{suffix}", dpi
    ))

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.plot(r, cumulative, linewidth=1.2)
    ax.set(
        xlabel="r / Å",
        ylabel="Cumulative N_water_O",
        title="Zn–Ow cumulative N(r)",
    )
    style_axis(ax)
    generated.append(save_figure(
        plt, fig, output_dir / f"Zn_Owater_cumulative_N.{suffix}", dpi
    ))
    return generated


def plot_rmsd(plt, rows, output_dir, suffix, dpi, time_key, time_label):
    x, y = numeric_pairs(rows, time_key, "solute_RMSD_A")
    if not y.size:
        return []
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(x, y, linewidth=1.1)
    ax.set(
        xlabel=time_label,
        ylabel="Solute RMSD / Å",
        title="Solute RMSD after Kabsch alignment",
    )
    style_axis(ax)
    return [save_figure(
        plt, fig, output_dir / f"solute_RMSD_vs_time.{suffix}", dpi
    )]


def plot_rmsf(plt, rmsf_rows, output_dir, suffix, dpi):
    x = numeric_values(rmsf_rows, "index")
    y = numeric_values(rmsf_rows, "RMSF_A")
    if not x.size or len(y) != len(x):
        return []
    fig, ax = plt.subplots(figsize=(max(7.0, len(x) * 0.18), 4.8))
    ax.plot(x, y, marker="o", markersize=3, linewidth=1.0)
    if len(x) <= 40:
        labels = []
        for row in rmsf_rows:
            atom_label = row.get("name") or row.get("element") or "atom"
            labels.append(f"{row.get('index')}:{atom_label}")
        ax.set_xticks(x, labels=labels, rotation=60, ha="right")
    ax.set(
        xlabel="Solute atom index",
        ylabel="RMSF / Å",
        title="Kabsch-aligned solute RMSF",
    )
    style_axis(ax)
    return [save_figure(
        plt, fig, output_dir / f"solute_RMSF.{suffix}", dpi
    )]


def plot_tetrahedrality(
    plt, rows, output_dir, suffix, dpi, time_key, time_label,
):
    x, y = numeric_pairs(rows, time_key, "q_tetra4_nearest")
    if not y.size:
        return []
    generated = []
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.plot(x, y, linewidth=1.1)
    ax.axhline(1.0, color="0.35", linestyle="--", linewidth=1.0,
               label="ideal tetrahedron, q = 1")
    ax.set(
        xlabel=time_label,
        ylabel="q_tetra",
        title="Four-nearest-donor tetrahedrality descriptor",
    )
    style_axis(ax)
    ax.legend()
    generated.append(save_figure(
        plt, fig, output_dir / f"tetrahedrality_vs_time.{suffix}", dpi
    ))

    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.hist(y, bins="auto", edgecolor="black", linewidth=0.5)
    ax.set(
        xlabel="q_tetra",
        ylabel="Frame count",
        title="Tetrahedrality distribution (descriptive)",
    )
    style_axis(ax)
    generated.append(save_figure(
        plt,
        fig,
        output_dir / "distributions" / f"tetrahedrality_histogram.{suffix}",
        dpi,
    ))
    return generated


def plot_energy(plt, rows, output_dir, suffix, dpi, time_key, time_label):
    generated = []
    for key, filename, ylabel, title in [
        ("energy_Eh_from_trj", "energy_vs_time", "Energy / Eh",
         "xTB trajectory energy"),
        ("gnorm_from_trj", "gnorm_vs_time", "Gradient norm / Eh a0⁻¹",
         "xTB trajectory gradient norm"),
    ]:
        x, y = numeric_pairs(rows, time_key, key)
        if y.size < 2:
            continue
        fig, ax = plt.subplots(figsize=(8.0, 4.8))
        ax.plot(x, y, linewidth=1.1)
        ax.set(xlabel=time_label, ylabel=ylabel, title=title)
        style_axis(ax)
        generated.append(save_figure(
            plt, fig, output_dir / f"{filename}.{suffix}", dpi
        ))
    return generated


def generate_plots(
    rows,
    sites,
    radial_rows,
    state_rows,
    rmsf_rows,
    output_dir: Path,
    plot_format: str,
    dpi: int,
    cutoff_A=None,
    time_in_ps=True,
):
    """Generate plots only from already-calculated analysis data."""
    plt = load_pyplot()
    output_dir.mkdir(parents=True, exist_ok=True)
    time_key = "time_ps" if time_in_ps else "global_frame"
    time_label = "Time / ps" if time_in_ps else "Analyzed frame"
    generated = []
    generated += plot_distance_timeseries(
        plt, rows, sites, output_dir, plot_format, dpi, cutoff_A,
        time_key, time_label,
    )
    generated += plot_distance_histograms(
        plt, rows, sites, output_dir, plot_format, dpi
    )
    generated += plot_nearest_water(
        plt, rows, output_dir, plot_format, dpi, cutoff_A,
        time_key, time_label,
    )
    generated += plot_coordination_number(
        plt, rows, output_dir, plot_format, dpi, time_key, time_label
    )
    generated += plot_coordination_states(
        plt, rows, state_rows, output_dir, plot_format, dpi,
        time_key, time_label,
    )
    generated += plot_radial_number(
        plt, radial_rows, output_dir, plot_format, dpi
    )
    generated += plot_rmsd(
        plt, rows, output_dir, plot_format, dpi, time_key, time_label
    )
    generated += plot_rmsf(
        plt, rmsf_rows, output_dir, plot_format, dpi
    )
    generated += plot_tetrahedrality(
        plt, rows, output_dir, plot_format, dpi, time_key, time_label
    )
    generated += plot_energy(
        plt, rows, output_dir, plot_format, dpi, time_key, time_label
    )
    return generated


class ResidenceTracker:
    def __init__(self, oxygen_indices, enter_A, exit_A):
        self.oxygen_indices = oxygen_indices
        self.enter = enter_A
        self.exit = exit_A
        self.active = {}
        self.events = []

    def update(self, time_ps, distances):
        for atom, d in zip(self.oxygen_indices, distances):
            if atom not in self.active and d <= self.enter:
                self.active[atom] = [time_ps, float(d)]
            elif atom in self.active and d >= self.exit:
                start, dmin = self.active.pop(atom)
                self.events.append({"water_O_index": atom+1, "start_ps": start, "end_ps": time_ps,
                                    "duration_ps": time_ps-start, "min_distance_A": dmin,
                                    "censored_at_end": False})
            elif atom in self.active:
                self.active[atom][1] = min(self.active[atom][1], float(d))

    def finish(self, time_ps):
        for atom, (start, dmin) in sorted(self.active.items()):
            self.events.append({"water_O_index": atom+1, "start_ps": start, "end_ps": time_ps,
                                "duration_ps": time_ps-start, "min_distance_A": dmin,
                                "censored_at_end": True})
        self.active.clear()


def parser():
    p = argparse.ArgumentParser(description="Analyze xTB MD trajectories")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--replica", type=Path)
    src.add_argument("--trajectory", type=Path, nargs="+")
    p.add_argument("--stage", action="append", help=f"Repeatable; default {DEFAULT_STAGE}")
    p.add_argument("--topology", type=Path)
    p.add_argument("--frame-dt-fs", type=float)
    p.add_argument("--discard-first-ps", type=float, default=0.0)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--output-dir", type=Path, default=Path("xtb_analysis"))
    p.add_argument("--zn-index", type=int)
    p.add_argument("--n-solute", type=int)
    p.add_argument("--site", action="append", type=parse_site, default=[], metavar="LABEL:INDEX[:GROUP]")
    p.add_argument("--water-resname", action="append", default=[])
    p.add_argument("--contact-cutoff-A", type=float,
                   help="Optional descriptive Zn--donor cutoff; choose from observed first-shell distribution")
    p.add_argument("--smooth-r0-A", type=float)
    p.add_argument("--smooth-exponent", type=int, default=6)
    p.add_argument("--water-enter-A", type=float)
    p.add_argument("--water-exit-A", type=float)
    p.add_argument("--radial-dr-A", type=float, default=0.05)
    p.add_argument("--radial-rmax-A", type=float, default=6.0)
    p.add_argument("--no-rmsd", action="store_true")
    p.add_argument("--no-travis-export", action="store_true")
    p.add_argument(
        "--alignment-selection",
        choices=("solute-heavy", "solute-all"),
        default="solute-heavy",
        help="Atoms used for aligned exports (default: solute-heavy)",
    )
    p.add_argument(
        "--alignment-indices",
        type=parse_alignment_indices,
        help="Explicit one-based alignment selection, e.g. 1-5,10,12-19",
    )
    p.add_argument(
        "--export-stride",
        type=int,
        default=1,
        help="Export every Nth analyzed frame (default: 1)",
    )
    p.add_argument("--export-start-ps", type=float)
    p.add_argument("--export-end-ps", type=float)
    p.add_argument("--no-plots", action="store_true",
                   help="Do not generate diagnostic plots")
    p.add_argument("--plot-format", choices=("png", "pdf", "svg"), default="png",
                   help="Plot file format (default: png)")
    p.add_argument("--plot-dpi", type=int, default=300,
                   help="PNG resolution in dots per inch (default: 300)")
    p.add_argument("--plot-dir", type=Path,
                   help="Plot directory (default: OUTPUT_DIR/plots)")
    return p


def main():
    args = parser().parse_args()
    if args.stride < 1 or args.discard_first_ps < 0:
        raise SystemExit("Invalid stride/discard")
    if args.export_stride < 1:
        raise SystemExit("--export-stride must be >= 1")
    if (args.export_start_ps is not None and args.export_end_ps is not None
            and args.export_end_ps < args.export_start_ps):
        raise SystemExit("--export-end-ps must be >= --export-start-ps")
    if args.plot_dpi < 1:
        raise SystemExit("--plot-dpi must be >= 1")
    if (args.water_enter_A is None) != (args.water_exit_A is None):
        raise SystemExit("Supply --water-enter-A and --water-exit-A together")
    if args.water_enter_A is not None and args.water_exit_A <= args.water_enter_A:
        raise SystemExit("--water-exit-A must be larger than --water-enter-A")

    specs, topology_path, nsolute = trajectory_specs(args)
    if (not args.no_travis_export
            and (args.export_start_ps is not None or args.export_end_ps is not None)
            and not all(spec.dt_fs is not None for spec in specs)):
        raise SystemExit(
            "--export-start-ps/--export-end-ps require known frame timing"
        )
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    topology = parse_pdb(topology_path) if topology_path else None

    first = next(iter_xtb_trj(specs[0].path), None)
    if first is None:
        raise RuntimeError("Empty trajectory")
    nat = len(first.elements)
    if topology and len(topology) != nat:
        raise RuntimeError("Topology and trajectory atom counts differ")
    if topology:
        mismatches = [
            i + 1 for i, atom in enumerate(topology)
            if atom.element != first.elements[i]
        ]
        if mismatches:
            raise RuntimeError(
                "Topology and trajectory element/order differ at atom indices: "
                f"{mismatches[:10]}"
            )
    if nsolute is not None and not 1 <= nsolute <= nat:
        raise RuntimeError(f"n_solute must be between 1 and {nat}")
    zn = unique_zn(first.elements, args.zn_index)
    if not 0 <= zn < nat:
        raise RuntimeError("Zn index out of range")
    for s in args.site:
        if s.index0 >= nat:
            raise RuntimeError(f"Site {s.label} index out of range")

    if args.no_travis_export:
        alignment_indices = None
        alignment_metadata = {
            "mode": "explicit-indices" if args.alignment_indices is not None
                    else args.alignment_selection,
            "indices_one_based": [],
            "available": False,
            "reason": "coordinate export disabled by --no-travis-export",
        }
    else:
        alignment_indices, alignment_metadata = resolve_alignment_selection(
            first.elements,
            first.xyz,
            nsolute,
            args.alignment_selection,
            args.alignment_indices,
        )

    water_names = DEFAULT_WATER_RESNAMES | {x.upper() for x in args.water_resname}
    water_O, warnings = water_oxygens(topology, first.elements, nsolute, water_names)
    candidate = list(dict.fromkeys([s.index0 for s in args.site] + water_O))

    atoms = []
    water_set = set(water_O)
    for i, e in enumerate(first.elements):
        a = topology[i] if topology else None
        atoms.append({"index": i+1, "element": e, "name": a.name if a else "",
                      "resname": a.resname if a else "", "resid": a.resid if a else "",
                      "is_solute": bool(nsolute is not None and i < nsolute),
                      "is_water_O": i in water_set, "is_Zn": i == zn})
    write_csv(out / "atoms.csv", atoms)

    rows = []
    radial = []
    state_rows = []
    rmsf_rows = []
    site_values = {s.label: [] for s in args.site}
    water_all = []
    state_counts = {}
    ref = None
    rmsf_sum = rmsf_sum2 = None
    rmsf_n = 0
    tracker = None
    if water_O and args.water_enter_A is not None:
        tracker = ResidenceTracker(water_O, args.water_enter_A, args.water_exit_A)

    if args.replica:
        system_label = args.replica.resolve().parent.name
    elif len(specs) == 1:
        system_label = specs[0].path.stem
    else:
        system_label = "multiple_sources"

    exporter = None
    if not args.no_travis_export:
        exporter = CoordinateExporter(
            output_dir=out / "travis",
            nsolute=nsolute,
            alignment_indices=alignment_indices,
            alignment_mode=alignment_metadata["mode"],
            system=system_label,
        )
        if nsolute is None:
            warnings.append(
                "n_solute is unavailable; solute-only coordinate exports were skipped."
            )
        if alignment_indices is None:
            warnings.append(
                "Alignment selection is unavailable; aligned coordinate exports were skipped."
            )
    elif (out / "travis").exists() or (out / "vmd").exists() or (
        out / "trajectory_for_travis.xyz"
    ).exists():
        warnings.append(
            "Existing TRAVIS/VMD export files were left untouched because "
            "--no-travis-export was supplied."
        )

    time_offset = 0.0
    final_time = 0.0
    analyzed = 0
    export_window_frames = 0
    source_velocity_stats = []

    try:
        for spec in specs:
            seen = 0
            frames_with_velocities = 0
            for fr in iter_xtb_trj(spec.path):
                seen += 1
                if fr.velocities is not None:
                    frames_with_velocities += 1
                if fr.elements != first.elements:
                    raise RuntimeError(f"Atom order/elements changed in {spec.path}")
                local_ps = fr.index * spec.dt_fs / 1000.0 if spec.dt_fs is not None else math.nan
                if spec.dt_fs is not None and local_ps < args.discard_first_ps:
                    continue
                if fr.index % args.stride:
                    continue
                analyzed += 1
                time_ps = time_offset + local_ps if spec.dt_fs is not None else float(analyzed-1)
                export_time_ps = time_ps if spec.dt_fs is not None else None
                final_time = time_ps
                z = fr.xyz[zn]
                row = {"global_frame": analyzed-1, "source_frame": fr.index, "stage": spec.stage,
                       "time_ps": time_ps, "stage_time_ps": local_ps,
                       "energy_Eh_from_trj": fr.metadata.get("energy", ""),
                       "gnorm_from_trj": fr.metadata.get("gnorm", ""),
                       "has_velocities": fr.velocities is not None}

                if (exporter is not None
                        and alignment_indices is not None
                        and exporter.reference_xyz is None):
                    exporter.set_reference(fr.xyz, {
                        "stage": spec.stage,
                        "source_frame": fr.index,
                        "global_frame": analyzed - 1,
                        "time_ps": export_time_ps,
                    })

                named_dist = []
                contacts = []
                for s in args.site:
                    d = float(np.linalg.norm(fr.xyz[s.index0] - z))
                    row[f"d_Zn_{s.label}_A"] = d
                    site_values[s.label].append(d)
                    named_dist.append(d)
                    if args.contact_cutoff_A is not None and d <= args.contact_cutoff_A:
                        contacts.append(s.label)

                wd = np.array([], float)
                if water_O:
                    wd = np.linalg.norm(fr.xyz[water_O] - z[None, :], axis=1)
                    water_all.append(wd.copy())
                    j = int(np.argmin(wd))
                    row["nearest_water_O_index"] = water_O[j] + 1
                    row["nearest_water_O_distance_A"] = float(wd[j])
                    if tracker is not None and spec.dt_fs is not None:
                        tracker.update(time_ps, wd)

                if args.contact_cutoff_A is not None:
                    nw = int(np.count_nonzero(wd <= args.contact_cutoff_A)) if wd.size else 0
                    row["CN_named_hard"] = len(contacts)
                    row["CN_water_hard"] = nw
                    row["CN_total_hard"] = len(contacts) + nw
                    state = "+".join(contacts) if contacts else "no_named_sites"
                    state += f"|Ow={nw}|CN={len(contacts)+nw}"
                    row["coordination_state"] = state
                    state_counts[state] = state_counts.get(state, 0) + 1

                if args.smooth_r0_A is not None:
                    row["CN_smooth"] = smooth_cn(np.r_[named_dist, wd], args.smooth_r0_A, args.smooth_exponent)

                if len(candidate) >= 4:
                    cxyz = fr.xyz[candidate]
                    order = np.argsort(np.linalg.norm(cxyz - z[None, :], axis=1))[:4]
                    row["q_tetra4_nearest"] = tetra_q(cxyz[order] - z[None, :])

                if not args.no_rmsd and nsolute is not None and nsolute >= 2:
                    sol = fr.xyz[:nsolute]
                    if ref is None:
                        ref = sol.copy()
                    aligned = kabsch(sol, ref)
                    delta = aligned - ref
                    row["solute_RMSD_A"] = float(np.sqrt(np.mean(np.sum(delta*delta, axis=1))))
                    if rmsf_sum is None:
                        rmsf_sum = np.zeros_like(aligned)
                        rmsf_sum2 = np.zeros_like(aligned)
                    rmsf_sum += aligned
                    rmsf_sum2 += aligned*aligned
                    rmsf_n += 1

                rows.append(row)
                if exporter is not None:
                    in_export_window = (
                        (args.export_start_ps is None
                         or export_time_ps >= args.export_start_ps)
                        and (args.export_end_ps is None
                             or export_time_ps <= args.export_end_ps)
                    )
                    if in_export_window:
                        if export_window_frames % args.export_stride == 0:
                            exporter.write(
                                fr.elements,
                                fr.xyz,
                                spec.stage,
                                fr.index,
                                analyzed - 1,
                                export_time_ps,
                            )
                        export_window_frames += 1
            if spec.dt_fs is not None:
                time_offset += seen * spec.dt_fs / 1000.0
            layout = velocity_layout(seen, frames_with_velocities)
            source_velocity_stats.append({
                "frames_total": seen,
                "frames_with_velocities": frames_with_velocities,
                "frames_without_velocities": seen - frames_with_velocities,
                "velocity_layout": layout,
            })
            if layout == "mixed":
                warnings.append(
                    f"Mixed velocity layout in {spec.path}: "
                    f"{frames_with_velocities} of {seen} frames contain velocity blocks; "
                    "verify trajectory concatenation and completeness."
                )
    finally:
        if exporter is not None:
            exporter.close()

    if not rows:
        raise RuntimeError("No frames retained")

    velocity_frames_total = sum(
        item["frames_total"] for item in source_velocity_stats
    )
    velocity_frames_with = sum(
        item["frames_with_velocities"] for item in source_velocity_stats
    )
    if velocity_frames_with == 0:
        velocity_status = "absent"
    elif velocity_frames_with == velocity_frames_total:
        velocity_status = "present in all frames"
    else:
        velocity_status = "mixed"
        if not any(
            item["velocity_layout"] == "mixed" for item in source_velocity_stats
        ):
            warnings.append(
                "Velocity blocks are present in only some trajectory sources; "
                "verify source compatibility."
            )

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    write_csv(out / "frames.csv", rows, fields)

    summaries = []
    for label, values in site_values.items():
        a = np.asarray(values)
        if a.size:
            summaries.append({"site": label, "n": len(a), "mean_A": a.mean(), "std_A": a.std(ddof=1) if len(a)>1 else 0,
                              "min_A": a.min(), "q05_A": np.quantile(a,.05), "median_A": np.median(a),
                              "q95_A": np.quantile(a,.95), "max_A": a.max()})
    if water_O:
        a = np.asarray([r["nearest_water_O_distance_A"] for r in rows])
        summaries.append({"site": "nearest_water_O", "n": len(a), "mean_A": a.mean(), "std_A": a.std(ddof=1) if len(a)>1 else 0,
                          "min_A": a.min(), "q05_A": np.quantile(a,.05), "median_A": np.median(a),
                          "q95_A": np.quantile(a,.95), "max_A": a.max()})
    write_csv(out / "distance_summary.csv", summaries)

    if water_all:
        edges = np.arange(0, args.radial_rmax_A + args.radial_dr_A, args.radial_dr_A)
        counts = np.zeros(len(edges)-1)
        for d in water_all:
            counts += np.histogram(d, bins=edges)[0]
        counts /= len(water_all)
        cumulative = np.cumsum(counts)
        radial = [{"r_left_A": edges[i], "r_right_A": edges[i+1], "r_mid_A": (edges[i]+edges[i+1])/2,
                   "mean_water_O_count_in_shell": counts[i], "cumulative_N_water_O": cumulative[i]}
                  for i in range(len(counts))]
        write_csv(out / "Zn_Owater_radial_number.csv", radial)

    if state_counts:
        state_rows = [{"state": s, "frames": n, "fraction_of_analyzed_frames": n/len(rows)}
                      for s, n in sorted(state_counts.items(), key=lambda x:(-x[1],x[0]))]
        write_csv(out / "coordination_states.csv", state_rows)

    if tracker:
        tracker.finish(final_time)
        write_csv(out / "water_contact_events.csv", tracker.events)

    if rmsf_n and rmsf_sum is not None:
        mean = rmsf_sum/rmsf_n
        var = np.maximum(rmsf_sum2/rmsf_n - mean*mean, 0)
        rmsf = np.sqrt(np.sum(var, axis=1))
        for i, v in enumerate(rmsf):
            a = topology[i] if topology else None
            rmsf_rows.append({"index": i+1, "element": first.elements[i], "name": a.name if a else "",
                              "resname": a.resname if a else "", "resid": a.resid if a else "", "RMSF_A": v})
        write_csv(out / "solute_RMSF.csv", rmsf_rows)

    plot_dir = args.plot_dir.resolve() if args.plot_dir else out / "plots"
    generated_plots = []
    if not args.no_plots:
        try:
            generated_plots = generate_plots(
                rows=rows,
                sites=args.site,
                radial_rows=radial,
                state_rows=state_rows,
                rmsf_rows=rmsf_rows,
                output_dir=plot_dir,
                plot_format=args.plot_format,
                dpi=args.plot_dpi,
                cutoff_A=args.contact_cutoff_A,
                time_in_ps=all(spec.dt_fs is not None for spec in specs),
            )
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc

    sources_metadata = [
        {
            "path": str(spec.path),
            "sha256": sha256(spec.path),
            "stage": spec.stage,
            "frame_dt_fs": spec.dt_fs,
            **velocity_stats,
        }
        for spec, velocity_stats in zip(specs, source_velocity_stats)
    ]

    travis_export_summary = {
        "enabled": exporter is not None,
        "directory": str(out / "travis") if exporter is not None else None,
        "exported_frames": 0,
        "metadata_path": None,
        "generated_files": [],
    }
    if exporter is not None:
        if exporter.exported_frames == 0:
            raise RuntimeError(
                "No analyzed frames matched the requested export time window"
            )

        validations = {}
        for name, path in exporter.paths.items():
            expected_atoms = nsolute if name.startswith("solute_") else nat
            validations[str(path.relative_to(out))] = validate_coordinate_xyz(
                path,
                exporter.exported_frames,
                expected_atoms,
            )

        legacy_path = out / "trajectory_for_travis.xyz"
        shutil.copy2(exporter.paths["full_coordinates"], legacy_path)
        validations[str(legacy_path.relative_to(out))] = validate_coordinate_xyz(
            legacy_path,
            exporter.exported_frames,
            nat,
        )

        topology_output = out / "travis" / "reference_topology.pdb"
        topology_provenance = prepare_reference_topology(
            topology_output,
            topology_path,
            first,
        )
        atom_map_path = out / "travis" / "atom_map.csv"
        write_atom_map(
            atom_map_path,
            first.elements,
            topology,
            nsolute,
            water_O,
            zn,
            args.site,
        )
        travis_readme = out / "travis" / "README_TRAVIS.md"
        write_travis_readme(travis_readme, system_label, exporter.paths)

        generated_files = list(exporter.paths.values()) + [
            legacy_path,
            topology_output,
            atom_map_path,
            travis_readme,
        ]
        if "full_aligned_on_solute" in exporter.paths:
            generated_files.extend(write_vmd_support(out / "vmd"))

        export_metadata_path = out / "travis" / "travis_export_metadata.json"
        generated_files.append(export_metadata_path)
        export_metadata = {
            "program": "xtb_analysis.py",
            "version": "0.3",
            "sources": sources_metadata,
            "stages": [spec.stage for spec in specs],
            "dt_fs": [spec.dt_fs for spec in specs],
            "frames_read": velocity_frames_total,
            "analyzed_frame_count": len(rows),
            "exported_frame_count": exporter.exported_frames,
            "export_stride": args.export_stride,
            "requested_start_ps": args.export_start_ps,
            "requested_end_ps": args.export_end_ps,
            "first_exported_time_ps": exporter.first_exported_time_ps,
            "last_exported_time_ps": exporter.last_exported_time_ps,
            "n_atoms": nat,
            "n_solute": nsolute,
            "coordinate_units": "angstrom",
            "velocity_blocks": {
                "status": velocity_status,
                "frames_with_velocities": velocity_frames_with,
                "frames_total": velocity_frames_total,
            },
            "velocities_exported": False,
            "alignment": {
                "performed": alignment_indices is not None,
                "selection": alignment_metadata,
                "reference_frame": exporter.reference_frame,
                "transform_application": (
                    "One rigid transform fitted on the selected atoms is "
                    "applied unchanged to every atom in each full frame."
                    if alignment_indices is not None else None
                ),
            },
            "topology_provenance": topology_provenance,
            "generated_files": [
                str(path.relative_to(out)) for path in generated_files
            ],
            "xyz_validation": validations,
            "interpretation": (
                "The exported XYZ/PDB files contain coordinates only. "
                "Velocity records present in xtb.trj were intentionally omitted."
            ),
        }
        export_metadata_path.write_text(json.dumps(export_metadata, indent=2))
        travis_export_summary = {
            "enabled": True,
            "directory": str(out / "travis"),
            "exported_frames": exporter.exported_frames,
            "metadata_path": str(export_metadata_path),
            "generated_files": [
                str(path.relative_to(out)) for path in generated_files
            ],
        }

    meta = {
        "program": "xtb_analysis.py",
        "version": "0.3",
        "sources": sources_metadata,
        "topology": {"path": str(topology_path), "sha256": sha256(topology_path)} if topology_path else None,
        "n_atoms": nat, "n_solute_atoms": nsolute, "Zn_index": zn+1,
        "water_O_indices": [i+1 for i in water_O],
        "named_sites": [{"label": s.label, "index": s.index0+1, "group": s.group} for s in args.site],
        "frames_analyzed": len(rows), "discard_first_ps_per_trajectory": args.discard_first_ps, "stride": args.stride,
        "contact_cutoff_A": args.contact_cutoff_A,
        "smooth_CN": {"r0_A": args.smooth_r0_A, "exponent": args.smooth_exponent,
                      "formula": "sum 1/(1+(r/r0)^n)"} if args.smooth_r0_A is not None else None,
        "radial_analysis": {"type": "Zn-centered shell counts and cumulative N(r)", "bulk_RDF_normalization": False,
                            "reason": "finite droplet / inhomogeneous wall-confined system"},
        "trajectory_scalar_provenance": {
            "energy_Eh_from_trj": "energy field parsed from each xtb.trj XYZ comment line",
            "gnorm_from_trj": "gnorm field parsed from each xtb.trj XYZ comment line"
        },
        "trajectory_velocity_data": {
            "present_in_any_frame": velocity_frames_with > 0,
            "present_in_all_frames": velocity_frames_with == velocity_frames_total,
            "analysis_performed": False,
            "units_assigned": False,
            "note": (
                "Velocity records were parsed and preserved in each Frame but "
                "were not used in the current structural analyses."
            )
        },
        "travis_export": travis_export_summary,
        "plotting": {
            "enabled": not args.no_plots,
            "format": args.plot_format,
            "dpi": args.plot_dpi,
            "directory": str(plot_dir),
            "generated_files": [str(path) for path in generated_plots],
            "histogram_bins": "numpy auto (descriptive only)"
        },
        "warnings": warnings,
        "interpretation_limits": [
            "State fractions from short screening trajectories are descriptive, not converged equilibrium populations.",
            "Water contact durations are descriptive, not converged kinetic residence times.",
            "Hard coordination states depend on the user-selected cutoff.",
            "q_tetra4 is a geometric descriptor, not a chemical-state assignment."
        ]
    }
    (out / "analysis_metadata.json").write_text(json.dumps(meta, indent=2))

    print(f"Analyzed frames : {len(rows)}")
    print(f"Atoms           : {nat}")
    print(f"Solute atoms    : {nsolute if nsolute is not None else 'unknown'}")
    print(f"Zn index        : {zn+1}")
    print(f"Water oxygens   : {len(water_O)}")
    print(f"Velocity blocks : {velocity_status}")
    print(f"Output          : {out}")
    if exporter is not None:
        print(f"Exported frames : {exporter.exported_frames}")
    if not args.no_plots:
        print(f"Plots           : {len(generated_plots)} in {plot_dir}")
    for w in warnings:
        print(f"WARNING: {w}")


if __name__ == "__main__":
    main()
