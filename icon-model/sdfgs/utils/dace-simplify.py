#!/usr/bin/env python
from pathlib import Path
import argparse

import dace


def main():
    parser = argparse.ArgumentParser(description=("Simplify an SDFG"))

    parser.add_argument(
        "unsimplified_sdfg_path",
        help="Path to the unsimplified input SDFG",
        type=Path,
    )
    parser.add_argument(
        "simplified_sdfg_path",
        help="Path to the simplified ouput SDFG",
        type=Path,
    )

    args = parser.parse_args()

    sdfg = dace.SDFG.from_file(str(args.unsimplified_sdfg_path))
    sdfg.simplify()

    output_path = args.simplified_sdfg_path.resolve()
    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)

    sdfg.save(str(output_path))


if __name__ == "__main__":
    main()
