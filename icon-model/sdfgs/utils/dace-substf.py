#!/usr/bin/env python
from __future__ import annotations
from typing import (
    Dict,
    List,
    Tuple,
    Set,
)
from collections import defaultdict
import re
import sys
import argparse
from pathlib import Path
from shutil import copy as copy_file

from yaml import load as load_yaml
try:
    from yaml import CLoader as YAML_Loader
except ImportError:
    from yaml import Loader as YAML_Loader


assert (
    (3, 7) <= sys.version_info
), f"Unsupported Python version ({sys.version_info}), 3.7 or higher is required"


COMPILER_DEFINE_ENABLE = "DACE_SUBST_ENABLE"
COMPILER_DEFINE_VERIFICATION_MODE = "DACE_SUBST_VERIFY"


def join_wrapped(elems: List[str], seperator: str, before: str, after: str) -> str:
    source = seperator.join(elems)
    if 0 < len(elems):
        source = before + source + after
    return source


_IDENTIFIER_REGEX = "\\w+"
_WHITE_SPACE_OPTIONAL_REGEX = "\\s*"
_WHITE_SPACE_MANDATORY_REGEX = "\\s+"
_MODULE_LOWER_CASE_REGEX = (
    _WHITE_SPACE_OPTIONAL_REGEX +
    "module" +
    _WHITE_SPACE_MANDATORY_REGEX +
    _IDENTIFIER_REGEX +
    _WHITE_SPACE_OPTIONAL_REGEX
)


def process(
    fortran_file_path: Path,
    processed_output_path: Path,
    integrations: Dict[str, List[Tuple[int, int]]],
    associations: Dict[str, Dict[str, Dict[int, Dict[int, Dict[str, str]]]]],
) -> None:

    with open(fortran_file_path) as fortran_file:
        lines = fortran_file.readlines()

    imports = defaultdict(lambda: set())

    insert_before = defaultdict(lambda: [])
    insert_after = defaultdict(lambda: [])
    for sdfg_name, insertions in integrations.items():

        if 0 < len(insertions):
            imports[f"mo_{sdfg_name}_bindings"].update((
                f"run_{sdfg_name}",
                f"run_{sdfg_name}_verification",
                f"verify_{sdfg_name}",
            ))

        for start, end in insertions:
            start, end = start, end

            print(associations[sdfg_name])
            before_src = generate_start_substitution_src(sdfg_name, associations[sdfg_name][fortran_file_path.name][start][end])
            after_src = generate_end_substitution_src(sdfg_name, associations[sdfg_name][fortran_file_path.name][start][end])

            # `line_nr` will be `int`, so convert back
            # (we really should use YAML instead, because it supports int pair keys...)
            insert_before[start].append(before_src)
            insert_after[end].append(after_src)

    processed_lines = []
    for line_nr, line in enumerate(lines):
        line_nr += 1 # fix 1-based line numbering

        if line_nr in insert_before:
            processed_lines.extend(insert_before[line_nr])

        processed_lines.append(line)

        if line_nr in insert_after:
            processed_lines.extend(insert_after[line_nr])

        if re.fullmatch(_MODULE_LOWER_CASE_REGEX, line, flags=re.IGNORECASE):
            processed_lines.append(generate_imports_src(imports))

    processed_output_content = "".join(processed_lines)
    with open(processed_output_path, "w+") as processed_output_file:
        processed_output_file.write(processed_output_content)


def generate_imports_src(imports: Dict[str, Set[str]]) -> str:
    return "\n".join(
        f"  USE {module}, ONLY: &\n    " +
        ', &\n    '.join(functions)
        for module, functions in imports.items()
    )


def generate_start_substitution_src(name: str, arguments: Dict[str, str]) -> str:
    arguments_str = join_wrapped(
        [f"{name} = {arg}" for name, arg in arguments.items()],
        seperator=", &\n    ",
        before=" &\n    ",
        after=" &\n  ",
    )
    return f"""\
#if defined({COMPILER_DEFINE_ENABLE})
#if defined({COMPILER_DEFINE_VERIFICATION_MODE})
  PRINT *, "Enter velocity tendencies"
  CALL run_{name}_verification({arguments_str})
  PRINT *, "Exit velocity tendencies"
#else
  PRINT *, "Enter velocity tendencies"
  CALL run_{name}({arguments_str})
  PRINT *, "Exit velocity tendencies"
#endif
#endif

#if !defined({COMPILER_DEFINE_ENABLE}) || defined({COMPILER_DEFINE_VERIFICATION_MODE})
"""


def generate_end_substitution_src(name: str, arguments: Dict[str, str]) -> str:
    arguments_str = join_wrapped(
        [f"{name} = {arg}" for name, arg in arguments.items()],
        seperator=", &\n    ",
        before=" &\n    ",
        after=" &\n  ",
    )
    return f"""\
#endif

#if defined({COMPILER_DEFINE_ENABLE}) && defined({COMPILER_DEFINE_VERIFICATION_MODE})
  CALL verify_{name}({arguments_str})
#endif
"""


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Preprocesses a Fortran file introducing substitutions of Fortran code "
            "with calls to SDFGs."
        )
    )

    parser.add_argument(
        "integrations_yaml_path",
        help="Path to `integrations.yaml` file",
        type=Path,
    )
    parser.add_argument(
        "associations_yaml_paths",
        help="Path to the associations yaml file",
        type=Path,
        nargs="*",
    )
    parser.add_argument(
        "f90_file_path",
        help="Path to the Fortran file to be processed",
        type=Path,
    )
    parser.add_argument(
        "processed_file_path",
        help="Path to the output Fortran file after processing",
        type=Path,
    )

    args = parser.parse_args()

    with open(args.integrations_yaml_path) as integrattions_yaml_file:
        integrations_yaml = load_yaml(integrattions_yaml_file, Loader=YAML_Loader)
    associations = {}
    for associations_yaml_path in args.associations_yaml_paths:
        associations_yaml_path = Path(associations_yaml_path)
        associations_yaml_file_name = associations_yaml_path.name
        assert associations_yaml_file_name.endswith("_associations.yaml")
        sdfg_name = associations_yaml_file_name[:-len("_associations.yaml")]
        with open(associations_yaml_path) as associations_yaml_file:
            associations[sdfg_name] = load_yaml(associations_yaml_file, Loader=YAML_Loader)

    # We check only by name, collisions seem unlikely
    path_normalized_integrations_yaml = {
        Path(path).name: integrations
        for path, integrations in integrations_yaml.items()
    }

    if args.f90_file_path.name not in path_normalized_integrations_yaml.keys():
        # no integrations for this file -> no processing needed -> copy as is
        copy_file(args.f90_file_path, args.processed_file_path)
        return

    process(
        args.f90_file_path,
        args.processed_file_path,
        path_normalized_integrations_yaml[args.f90_file_path.name],
        associations,
    )


if __name__ == "__main__":
    main()
