#!/usr/bin/env python
from pathlib import Path
import argparse

import dace

SDFG_NAME = "velocity_tendencies"


def optimize(sdfg: dace.SDFG) -> dace.SDFG:
    # PLACEHOLDER: here goes your optimization pipeline: ...

    return sdfg


def main():
    parser = argparse.ArgumentParser(description=(f"Optimize '{SDFG_NAME}'"))

    parser.add_argument(
        "simplified_sdfg_path",
        help="Path to the simplified input SDFG",
        type=Path,
    )
    parser.add_argument(
        "optimized_sdfg_path",
        help="Output path to the optimized SDFG",
        type=Path,
    )

    args = parser.parse_args()

    sdfg = dace.SDFG.from_file(str(args.simplified_sdfg_path))

    sdfg = optimize(sdfg)

    # Makefile does this automatically no?
    # output_path = args.optimized_sdfg_path.resolve()
    # if not output_path.parent.exists():
    #    output_path.parent.mkdir(parents=True)

    sdfg.save(str(args.optimized_sdfg_path))


if __name__ == "__main__":
    main()
