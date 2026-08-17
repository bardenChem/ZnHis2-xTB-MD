#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_xyz_frames(path):
    with Path(path).open() as f:
        frame = 0

        while True:
            line = f.readline()

            if not line:
                break

            line = line.strip()

            if not line:
                continue

            n_atoms = int(line)

            comment = f.readline().rstrip("\n")

            atoms = []

            for _ in range(n_atoms):
                atom_line = f.readline()

                if not atom_line:
                    raise RuntimeError(
                        f"Trajectory ended inside frame {frame}"
                    )

                atoms.append(atom_line.rstrip("\n"))

            yield frame, n_atoms, comment, atoms

            frame += 1


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("trajectory")
    parser.add_argument("frame", type=int)

    parser.add_argument(
        "-o",
        "--output",
        required=True,
    )

    args = parser.parse_args()

    for frame, n_atoms, comment, atoms in read_xyz_frames(
        args.trajectory
    ):
        if frame == args.frame:

            with open(args.output, "w") as out:
                out.write(f"{n_atoms}\n")
                out.write(
                    f"extracted_frame={frame} "
                    f"original_comment={comment}\n"
                )

                for line in atoms:
                    out.write(line + "\n")

            print(
                f"Extracted frame {frame}: "
                f"{n_atoms} atoms -> {args.output}"
            )

            return

    raise SystemExit(
        f"Frame {args.frame} not found."
    )


if __name__ == "__main__":
    main()
