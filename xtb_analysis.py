#!/usr/bin/env python3
"""xtb_analysis.py

Initial post-processing toolkit for xTB molecular-dynamics trajectories.

Focus of v0.1:
- xTB extended-XYZ trajectories (xtb.trj)
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

No bulk-normalized RDF is computed for finite droplets because a spherical,
inhomogeneous droplet with a wall is not equivalent to periodic bulk liquid.
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


def iter_xtb_trj(path: Path) -> Iterator[Frame]:
    """Stream an xTB extended-XYZ trajectory."""
    with path.open("r", errors="replace") as f:
        iframe = 0
        while True:
            line = f.readline()
            while line and not line.strip():
                line = f.readline()
            if not line:
                return
            try:
                nat = int(line.strip())
            except ValueError as exc:
                raise RuntimeError(f"{path}: expected atom count at frame {iframe}") from exc
            comment = f.readline()
            if not comment:
                raise RuntimeError(f"{path}: missing comment at frame {iframe}")
            xyz = np.empty((nat, 3), float)
            elements = []
            for i in range(nat):
                fields = f.readline().split()
                if len(fields) < 4:
                    raise RuntimeError(f"{path}: malformed atom record in frame {iframe}")
                elements.append(element_name(fields[0]))
                xyz[i] = [float(v.replace("D", "E").replace("d", "e")) for v in fields[1:4]]
            yield Frame(iframe, comment.rstrip(), elements, xyz, parse_comment(comment))
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


def kabsch(mobile, reference):
    mc = mobile.mean(0)
    rc = reference.mean(0)
    mob = mobile - mc
    ref = reference - rc
    u, _, vt = np.linalg.svd(mob.T @ ref)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        vt[-1] *= -1
        rot = u @ vt
    return mob @ rot + rc


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
    return p


def main():
    args = parser().parse_args()
    if args.stride < 1 or args.discard_first_ps < 0:
        raise SystemExit("Invalid stride/discard")
    if (args.water_enter_A is None) != (args.water_exit_A is None):
        raise SystemExit("Supply --water-enter-A and --water-exit-A together")
    if args.water_enter_A is not None and args.water_exit_A <= args.water_enter_A:
        raise SystemExit("--water-exit-A must be larger than --water-enter-A")

    specs, topology_path, nsolute = trajectory_specs(args)
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    topology = parse_pdb(topology_path) if topology_path else None

    first = next(iter_xtb_trj(specs[0].path), None)
    if first is None:
        raise RuntimeError("Empty trajectory")
    nat = len(first.elements)
    if topology and len(topology) != nat:
        raise RuntimeError("Topology and trajectory atom counts differ")
    zn = unique_zn(first.elements, args.zn_index)
    if not 0 <= zn < nat:
        raise RuntimeError("Zn index out of range")
    for s in args.site:
        if s.index0 >= nat:
            raise RuntimeError(f"Site {s.label} index out of range")

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
    site_values = {s.label: [] for s in args.site}
    water_all = []
    state_counts = {}
    ref = None
    rmsf_sum = rmsf_sum2 = None
    rmsf_n = 0
    tracker = None
    if water_O and args.water_enter_A is not None:
        tracker = ResidenceTracker(water_O, args.water_enter_A, args.water_exit_A)

    travis = None if args.no_travis_export else (out / "trajectory_for_travis.xyz").open("w")
    time_offset = 0.0
    final_time = 0.0
    analyzed = 0

    try:
        for spec in specs:
            seen = 0
            for fr in iter_xtb_trj(spec.path):
                seen += 1
                if fr.elements != first.elements:
                    raise RuntimeError(f"Atom order/elements changed in {spec.path}")
                local_ps = fr.index * spec.dt_fs / 1000.0 if spec.dt_fs is not None else math.nan
                if spec.dt_fs is not None and local_ps < args.discard_first_ps:
                    continue
                if fr.index % args.stride:
                    continue
                analyzed += 1
                time_ps = time_offset + local_ps if spec.dt_fs is not None else float(analyzed-1)
                final_time = time_ps
                z = fr.xyz[zn]
                row = {"global_frame": analyzed-1, "source_frame": fr.index, "stage": spec.stage,
                       "time_ps": time_ps, "stage_time_ps": local_ps,
                       "energy_Eh_from_trj": fr.metadata.get("energy", ""),
                       "gnorm_from_trj": fr.metadata.get("gnorm", "")}

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
                if travis:
                    write_xyz_frame(travis, fr.elements, fr.xyz,
                                    f"stage={spec.stage} frame={fr.index} time_ps={time_ps:.6f}")
            if spec.dt_fs is not None:
                time_offset += seen * spec.dt_fs / 1000.0
    finally:
        if travis:
            travis.close()

    if not rows:
        raise RuntimeError("No frames retained")

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
        rmsf_rows = []
        for i, v in enumerate(rmsf):
            a = topology[i] if topology else None
            rmsf_rows.append({"index": i+1, "element": first.elements[i], "name": a.name if a else "",
                              "resname": a.resname if a else "", "resid": a.resid if a else "", "RMSF_A": v})
        write_csv(out / "solute_RMSF.csv", rmsf_rows)

    meta = {
        "program": "xtb_analysis.py",
        "version": "0.1",
        "sources": [{"path": str(s.path), "sha256": sha256(s.path), "stage": s.stage, "frame_dt_fs": s.dt_fs} for s in specs],
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
    print(f"Zn index        : {zn+1}")
    print(f"Water oxygens   : {len(water_O)}")
    print(f"Output          : {out}")
    for w in warnings:
        print(f"WARNING: {w}")


if __name__ == "__main__":
    main()
