#!/usr/bin/env python
from pathlib import Path
import argparse

from yaml import load as load_yaml

try:
    from yaml import CLoader as YAML_Loader
except ImportError:
    from yaml import Loader as YAML_Loader


SIMPLE_FORMATERS = {
    "sdfg_names": lambda sdfg_name: sdfg_name,
    "f90_interfaces": lambda sdfg_name: f"sdfgs/build/{sdfg_name}_bindings.f90",
}


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extracts the integrations specified in an `integrations.yaml`, and prints "
            "them in a suitable format."
        )
    )

    parser.add_argument(
        "integrations_yaml_path",
        help="Path to `integrations.yaml` file",
        type=Path,
    )
    subparsers = parser.add_subparsers(
        title="format",
        description="The format to print the integrations",
        required=True,
        dest="format",
    )

    for name in SIMPLE_FORMATERS.keys():
        subparsers.add_parser(name)
    subparsers.add_parser("icon_f90_files")

    association_files_parser = subparsers.add_parser("association_files")
    association_files_parser.add_argument(
        "fortran_file_path",
        help="The Fortran file for which to print the association files",
        type=Path,
    )

    args = parser.parse_args()

    with open(args.integrations_yaml_path) as yaml_file:
        integrations_yaml = load_yaml(yaml_file, Loader=YAML_Loader)

    if args.format in SIMPLE_FORMATERS:

        sdfg_integrations = set(
            SIMPLE_FORMATERS[args.format](sdfg_name)
            for fortran_source_integrations in integrations_yaml.values()
            for sdfg_name in fortran_source_integrations.keys()
        )

        print("\n".join(sdfg_integrations))

    elif args.format == "icon_f90_files":
        for icon_f90_file in integrations_yaml.keys():
            # Normalizes the path (e.g., removes "./" in "./relative_path/src_file.f90")
            print(str(Path(icon_f90_file)))

    elif args.format == "association_files":

        association_files = set(
            f"sdfgs/{sdfg_name}_associations.yaml"
            for fotran_source_file, fortran_source_integrations in integrations_yaml.items()
            for sdfg_name in fortran_source_integrations.keys()
            if Path(fotran_source_file) == args.fortran_file_path
        )
        print("\n".join(association_files))

    else:
        assert False, f"Requested invalid format '{args.format}'"


if __name__ == "__main__":
    main()
