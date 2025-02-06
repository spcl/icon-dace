#!/bin/sh

# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

""":"
for cmd in python3 python; do
  if command -v > /dev/null "${cmd}"; then
    exec "${cmd}" "$0" "$@"
  fi
done
echo "Error: could not find a python interpreter!" >&2
exit 1
":"""

import argparse
import os
import sys

from depgen import exhaust, map23, open23, zip_longest23


class ArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        try:
            # Drop everything after the first occurrence of #:
            arg_line = arg_line[: arg_line.index("#")]
        except ValueError:
            pass

        result = []
        # Do not regard consecutive whitespaces as a single separator:
        for arg in arg_line.split(" "):
            if arg:
                result.append(arg)
            elif result:
                # The previous argument has a significant whitespace:
                result[-1] += " "
        return result


def parse_args():
    parser = ArgumentParser(
        fromfile_prefix_chars="@",
        description="Generates OUTPUT makefile containing dependency rules for "
        "the INPUT source file. Recognizes preprocessor `#include`, `#if` and "
        "associated directives as well as Fortran `INCLUDE`, `USE` and "
        "`MODULE` statements.",
    )

    def comma_splitter(s):
        return list(filter(None, s.lower().split(",")))

    def path_splitter(s):
        return list(filter(None, s.split(":")))

    def default_obj_name(src_name):
        if src_name:
            src_no_ext_basename = os.path.splitext(os.path.basename(src_name))[
                0
            ]
            if src_no_ext_basename:
                return src_no_ext_basename + ".o"
        return None

    parser.add_argument(
        "--input",
        "-i",
        metavar="INPUT",
        nargs="+",
        help="input source file; if not specified, the program reads from the "
        "standard input stream",
    )
    parser.add_argument(
        "--output",
        "-o",
        metavar="OUTPUT",
        nargs="+",
        help="output makefile with generated dependency rules; if not "
        "specified, the program writes to the standard output stream",
    )
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        help="dump debug information to OUTPUT (commented with `#`)",
    )
    parser.add_argument(
        "--src-name",
        metavar="SRC_NAME",
        nargs="+",
        help="name of the source file, the prerequisite of the corresponding "
        "compilation rule as it will appear in the OUTPUT; normally (and by "
        "default) equals to the INPUT or to an empty string when the latter is "
        "set to the standard input stream",
    )
    parser.add_argument(
        "--obj-name",
        metavar="OBJ_NAME",
        nargs="+",
        help="name of the object file, the target of the corresponding "
        "compilation rule as it will appear in the OUTPUT; normally equals to "
        "the path to the object file that is supposed to be generated as a "
        "result of compilation (default: SRC_NAME without the directory-part "
        "and the file extension replaced with `.o`)",
    )
    parser.add_argument(
        "--dep-name",
        metavar="DEP_NAME",
        nargs="+",
        help="name of the generated makefile, the additional target of the "
        "corresponding compilation rule (for automatic dependency generation) "
        "as it will appear in the OUTPUT; normally (and by default) equals to "
        "the OUTPUT or to an empty string when the latter is set to the "
        "standard output stream)",
    )
    parser.add_argument(
        "--src-roots",
        metavar="SRC_ROOTS",
        type=path_splitter,
        help="colon-separated list of paths to directories; if specified and "
        "not empty, dependencies on files that do not reside in one of the "
        "specified directories will be ignored; applies only to files included "
        "using the preprocessor `#include` directive or the Fortran `INCLUDE` "
        "statement",
    )
    parser.add_argument(
        "--lc-enable",
        action="store_true",
        help="enable recognition of the preprocessor line control directives "
        "and generation of additional dependencies based on the detected "
        "filenames",
    )
    parser.add_argument(
        "flags",
        metavar="-- [$CPPFLAGS | $FCFLAGS]",
        nargs="?",
        default=argparse.SUPPRESS,
        help="actual flags to be used in compilation, i.e. $(CPPFLAGS) or "
        "$(FCFLAGS), must be given at the end of the command line following "
        "the double dash separator (--); the program searches these flags for "
        "(possibly multiple instances of) PP_INC_FLAG, PP_MACRO_FLAG, "
        "FC_INC_FLAG and FC_MOD_DIR_FLAG; any values found are used in the "
        "dependency generation (in the case of FC_MOD_DIR_FLAG, only the last "
        "value found is used)",
    )

    pp_arg_group = parser.add_argument_group("preprocessor arguments")
    pp_arg_group.add_argument(
        "--pp-enable",
        action="store_true",
        help="enable the preprocessing stage; if disabled (default), all "
        "arguments of this argument group are ignored",
    )
    pp_arg_group.add_argument(
        "--pp-eval-expr",
        action="store_true",
        help="enable evaluation of expressions that appear in preprocessor "
        "directives `#if` and `#elif` (does not apply to `#ifdef` and "
        "`#ifndef`, which are always evaluated); if disabled (default) or "
        "evaluation fails, both branches of the directives are included by the "
        "preprocessing stage",
    )
    pp_arg_group.add_argument(
        "--pp-inc-sys",
        action="store_true",
        help="enable recognition of dependencies specified with the "
        "angle-bracket form of the preprocessor `#include` directive (i.e. "
        "`#include <filename>`); the constraint set by SRC_ROOTS applies",
    )
    pp_arg_group.add_argument(
        "--pp-inc-order",
        default="inc,flg",
        metavar="ORDER_LIST",
        type=comma_splitter,
        help="directory search order of files included using the quoted form "
        "of the preprocessor `#include` directive (i.e. "
        '`#include "filename"`) ; ORDER_LIST is an ordered comma-separated '
        "list of keywords, the corresponding search paths of which are to be "
        "searched in the given order. The recognized keywords are: `cwd` (for "
        "the current working directory), `flg` (for the directories specified "
        "with PP_INC_FLAG compiler flag), `src` (for the directory containing "
        "the INPUT source file), and `inc` (for the directory containing the "
        "file with the `#include` directive). Default: `%(default)s`.",
    )
    pp_arg_group.add_argument(
        "--pp-inc-sys-order",
        default="flg",
        metavar="ORDER_LIST",
        type=comma_splitter,
        help="equivalent to the `--pp-inc-order` argument, only for the "
        "angle-bracket form of the preprocessor `#include` directive (i.e. "
        "`#include <filename>`, default: `%(default)s`)",
    )
    pp_arg_group.add_argument(
        "--pp-macro-flag",
        metavar="PP_MACRO_FLAG",
        default="-D",
        help="preprocessor flag used for macro definition; only flags that "
        "start with a single dash (-) and have no more than one trailing "
        "whitespace are supported (default: `%(default)s`)",
    )
    pp_arg_group.add_argument(
        "--pp-inc-flag",
        metavar="PP_INC_FLAG",
        default="-I",
        help="preprocessor flag used for setting search paths for the "
        "`#include` directive; only flags that start with a single dash (-) "
        "and have no more than one trailing whitespace are supported (default: "
        "`%(default)s`)",
    )

    fc_arg_group = parser.add_argument_group("Fortran arguments")
    fc_arg_group.add_argument(
        "--fc-enable",
        action="store_true",
        help="enable recognition of Fortran dependencies specified with "
        "`INCLUDE`, `USE` and `MODULE` statements; if disabled (default), all "
        "arguments of this argument group are ignored",
    )
    fc_arg_group.add_argument(
        "--fc-mod-stamp-name",
        metavar="FC_MOD_STAMP_NAME",
        nargs="+",
        help="name of the Fortran module stamp file (a.k.a witness or anchor), "
        "the prerequisite of the Fortran module and submodule files that are "
        "generated as a result of, and an extra target of the Fortran module "
        "and submodule files that are required for the compilation of SRC_NAME "
        "as it will appear in the OUTPUT; normally (and by default) equals to "
        "the OBJ_NAME",
    )
    fc_arg_group.add_argument(
        "--fc-mod-ext",
        default="mod",
        help="filename extension (without leading dot) of compiler-generated "
        "Fortran module files (default: `%(default)s`)",
    )
    fc_arg_group.add_argument(
        "--fc-mod-upper",
        choices=["yes", "no"],
        default="no",
        help="whether Fortran compiler-generated module files have uppercase "
        "names (default: `%(default)s`)",
    )
    fc_arg_group.add_argument(
        "--fc-smod-ext",
        default="smod",
        help="filename extension (without leading dot) of compiler-generated "
        "Fortran submodule files (default: `%(default)s`)",
    )
    fc_arg_group.add_argument(
        "--fc-smod-infix",
        default="@",
        help="filename infix of compiler-generated Fortran submodule files, "
        "i.e. a string in the basename of the submodule filename that appears "
        "between the name of the ancestor module and the name of the "
        "submodule; an empty value of the argument means that the compiler "
        "does not prefix submodule filenames with the names of their ancestor "
        "modules (default: `%(default)s`)",
    )
    fc_arg_group.add_argument(
        "--fc-root-smod",
        choices=["yes", "no"],
        default="yes",
        help="whether Fortran compiler generates submodule files for the root "
        "module ancestors (default: `%(default)s`)",
    )
    fc_arg_group.add_argument(
        "--fc-inc-order",
        default="src,flg",
        metavar="ORDER_LIST",
        type=comma_splitter,
        help="equivalent to the `--pp-inc-order` argument, only for the "
        "Fortran `INCLUDE` statement and FC_INC_FLAG (default: `%(default)s`)",
    )
    fc_intrisic_mods_default = (
        "iso_c_binding,iso_fortran_env,ieee_exceptions,"
        "ieee_arithmetic,ieee_features,omp_lib,"
        "omp_lib_kinds,openacc"
    )
    fc_arg_group.add_argument(
        "--fc-intrinsic-mods",
        metavar="FC_INTRINSIC_MODS_LIST",
        type=comma_splitter,
        action="append",
        help="comma-separated list of Fortran intrinsic modules. Fortran "
        "modules that are explicitly specified as intrinsic in the source file "
        "(i.e. `USE, INTRINSIC :: MODULENAME`) are ignored regardless of "
        "whether they are mentioned on the FC_INTRINSIC_MODS_LIST. Fortran "
        "modules that are mentioned on the FC_INTRINSIC_MODS_LIST are ignored "
        "only when their nature is not specified in the source file at all "
        "(i.e. `USE :: MODULENAME`). Fortran modules that need to be ignored "
        "unconditionally must be put on the FC_EXTERNAL_MODS_LIST (see "
        "`--fc-external-mods`). Default: `{0}`.".format(
            fc_intrisic_mods_default
        ),
    )
    fc_arg_group.add_argument(
        "--fc-external-mods",
        metavar="FC_EXTERNAL_MODS_LIST",
        type=comma_splitter,
        action="append",
        help="comma-separated list of external (to the project) Fortran "
        "modules that need to be unconditionally ignored when generating "
        "dependency rules (see also `--fc-intrinsic-mods`)",
    )
    fc_arg_group.add_argument(
        "--fc-mod-dir-flag",
        metavar="FC_MOD_DIR_FLAG",
        default="-J",
        help="Fortran compiler flag used to specify the directory where module "
        "files are saved; only flags that start with a single dash (-) and "
        "have no more than one trailing whitespace are supported (default: "
        "`%(default)s`)",
    )
    fc_arg_group.add_argument(
        "--fc-inc-flag",
        metavar="FC_INC_FLAG",
        default="-I",
        help="preprocessor flag used for setting search paths for the Fortran "
        "`INCLUDE` statement; only flags that start with a single dash (-) and "
        "have no more than one trailing whitespace are supported (default: "
        "`%(default)s`)",
    )

    unknown = []
    try:
        sep_idx = sys.argv.index("--")
        args = parser.parse_args(sys.argv[1:sep_idx])
        unknown = sys.argv[sep_idx + 1 :]
    except ValueError:
        args = parser.parse_args()

    if not args.input:
        args.input = (None,)

    if not args.src_name:
        args.src_name = args.input
    elif len(args.src_name) != len(args.input):
        parser.error(
            "number of SRC_NAME values is not equal to the number of "
            "INPUT values"
        )

    if not args.obj_name:
        args.obj_name = map23(default_obj_name, args.src_name)
    elif len(args.obj_name) != len(args.input):
        parser.error(
            "number of OBJ_NAME values is not equal to the number of "
            "INPUT values"
        )

    if not args.output:
        args.output = ()
    elif len(args.output) != len(args.input):
        parser.error(
            "number of OUTPUT values is not equal to the number of "
            "INPUT values"
        )

    if not args.dep_name:
        args.dep_name = args.output
    elif len(args.dep_name) != len(args.input):
        parser.error(
            "number of DEP_NAME values is not equal to the number of "
            "INPUT values"
        )

    compiler_arg_dests = dict()

    if args.pp_enable:
        compiler_arg_dests.update(
            pp_inc_dirs=args.pp_inc_flag, pp_macros=args.pp_macro_flag
        )

    if args.fc_enable:
        compiler_arg_dests.update(
            fc_inc_dirs=args.fc_inc_flag, fc_mod_dir=args.fc_mod_dir_flag
        )

    if compiler_arg_dests:
        compiler_args = dict()
        for dest, flag in compiler_arg_dests.items():
            if not flag.startswith("-") or flag.endswith("  "):
                parser.error("unsupported compiler/preprocessor flag " + flag)
            # Several dests might share the same flag and we want them to share
            # the same list of values in this the case:
            val_list = compiler_args.get(flag, None)
            if val_list is None:
                val_list = []
                compiler_args[flag] = val_list
            setattr(args, dest, val_list)

        appended_val_lists = []
        for arg in unknown:
            if arg.startswith("-"):
                appended_val_lists *= 0
                arg_ws = arg + " "
                for flag, val_list in compiler_args.items():
                    if flag == arg or flag == arg_ws:
                        # If the current argument equals to a flag, which might
                        # have a significant trailing whitespace, the next
                        # argument on the command line is the flag's value:
                        appended_val_lists.append(val_list)
                    elif arg.startswith(flag):
                        # If the current argument starts with a flag that does
                        # not have a trailing whitespace, the suffix of the
                        # argument is the flag's value:
                        val_list.append(arg[len(flag) :])
            elif appended_val_lists:
                for val_list in appended_val_lists:
                    val_list.append(arg)
                appended_val_lists *= 0

    if args.pp_enable and args.pp_macros:
        import re

        predefined_macros = dict()
        for m in args.pp_macros:
            match = re.match(r"^=*([a-zA-Z_]\w*)(\(.*\))?(?:=(.+))?$", m)
            if match:
                name = match.group(1)
                if name != "defined":
                    body = match.group(3) if match.group(3) else "1"
                    predefined_macros[name] = (match.group(2), body)
        args.pp_macros = predefined_macros

    if args.fc_enable:
        if not args.fc_mod_stamp_name:
            args.fc_mod_stamp_name = args.obj_name
        elif len(args.fc_mod_stamp_name) != len(args.input):
            parser.error(
                "number of FC_MOD_STAMP_NAME values is not equal to the number "
                "of INPUT values"
            )

        args.fc_mod_upper = args.fc_mod_upper == "yes"
        args.fc_root_smod = args.fc_root_smod == "yes"
        args.fc_mod_dir = args.fc_mod_dir[-1] if args.fc_mod_dir else None

        if args.fc_intrinsic_mods:
            args.fc_intrinsic_mods = [
                m for sublist in args.fc_intrinsic_mods for m in sublist
            ]
        else:
            args.fc_intrinsic_mods = comma_splitter(fc_intrisic_mods_default)

        if args.fc_external_mods:
            args.fc_external_mods = [
                m for sublist in args.fc_external_mods for m in sublist
            ]
    else:
        args.fc_mod_stamp_name = ()

    return args


def main():
    args = parse_args()

    included_files = set()

    def include_callback(filename):
        included_files.add(filename)

    lc_files = set()

    provided_modules = set()
    required_modules = set()
    provided_submodules = set()
    required_submodules = set()

    pp_debug_info = None
    lc_debug_info = None
    ftn_debug_info = None

    def init_debug_info(section):
        return ["#\n# {0}:\n".format(section)]

    def format_debug_line(line, msg):
        return "#  `{0}`:\t{1}\n".format(line.rstrip("\n"), msg)

    parser = None
    if args.pp_enable:
        from depgen.preprocessor import Parser

        parser = Parser(
            include_order=args.pp_inc_order,
            include_sys_order=args.pp_inc_sys_order,
            include_dirs=args.pp_inc_dirs,
            include_roots=args.src_roots,
            try_eval_expr=args.pp_eval_expr,
            inc_sys=args.pp_inc_sys,
            predefined_macros=args.pp_macros,
            subparser=parser,
        )

        parser.include_callback = include_callback

        if args.debug:
            parser.debug_callback = lambda line, msg: pp_debug_info.append(
                format_debug_line(line, msg)
            )

    if args.lc_enable:
        from depgen.line_control import Parser

        parser = Parser(include_roots=args.src_roots, subparser=parser)
        parser.lc_callback = lambda filename: lc_files.add(filename)

        if args.debug:
            parser.debug_callback = lambda line, msg: lc_debug_info.append(
                format_debug_line(line, msg)
            )

    if args.fc_enable:
        from depgen.fortran import Parser

        parser = Parser(
            include_order=args.fc_inc_order,
            include_dirs=args.fc_inc_dirs,
            include_roots=args.src_roots,
            intrinsic_mods=args.fc_intrinsic_mods,
            external_mods=args.fc_external_mods,
            subparser=parser,
        )

        parser.include_callback = include_callback
        parser.module_start_callback = lambda module: provided_modules.add(
            module
        )
        parser.submodule_start_callback = (
            lambda submodule, parent, module: provided_submodules.add(
                (module, submodule)
            )
            or (
                required_submodules.add((module, parent))
                if parent or args.fc_root_smod
                else required_modules.add(module)
            )
        )
        parser.module_use_callback = lambda module: required_modules.add(module)

        if args.fc_root_smod:
            parser.extendable_module_callback = (
                lambda module: provided_submodules.add((module, None))
            )

        if args.debug:
            parser.debug_callback = lambda line, msg: ftn_debug_info.append(
                format_debug_line(line, msg)
            )

    for inp, out, src_name, obj_name, dep_name, mod_stamp_name in zip_longest23(
        args.input,
        args.output,
        args.src_name,
        args.obj_name,
        args.dep_name,
        args.fc_mod_stamp_name,
    ):
        in_stream, in_stream_close = (
            (sys.stdin, False) if inp is None else (open23(inp), True)
        )

        if args.debug:
            if args.pp_enable:
                pp_debug_info = init_debug_info("Preprocessor")
            if args.lc_enable:
                lc_debug_info = init_debug_info("Line control")
            if args.fc_enable:
                ftn_debug_info = init_debug_info("Fortran")

        if parser:
            exhaust(parser.parse(in_stream, in_stream.name))

        not in_stream_close or in_stream.close()

        out_lines = gen_lc_deps(src_name, lc_files)

        include_targets = [obj_name, dep_name]
        if obj_name != mod_stamp_name:
            include_targets.append(mod_stamp_name)

        out_lines.extend(
            gen_include_deps(include_targets, src_name, included_files)
        )

        if (
            provided_modules
            or required_modules
            or provided_submodules
            or required_submodules
        ):
            out_lines.extend(
                gen_module_deps(
                    obj_name,
                    mod_stamp_name,
                    provided_modules,
                    required_modules,
                    provided_submodules,
                    required_submodules,
                    args.fc_mod_dir,
                    args.fc_mod_upper,
                    args.fc_mod_ext,
                    args.fc_smod_infix,
                    args.fc_smod_ext,
                )
            )

        if args.debug:
            out_lines.extend(
                [
                    "\n# Python version: ",
                    sys.version.replace("\n", " "),
                    "\n#\n",
                    "# Command:\n",
                    "#  ",
                    " ".join(sys.argv),
                    "\n#\n",
                    "# Parsed arguments:\n#  ",
                    "\n#  ".join(
                        [k + "=" + str(v) for k, v in vars(args).items()]
                    ),
                    "\n",
                ]
            )
            if pp_debug_info is not None:
                out_lines.extend(pp_debug_info)
            if lc_debug_info is not None:
                out_lines.extend(lc_debug_info)
            if ftn_debug_info is not None:
                out_lines.extend(ftn_debug_info)
            out_lines.append("\n")

        out_stream, out_stream_close = (
            (sys.stdout, False) if out is None else (open23(out, "w"), True)
        )
        out_stream.writelines(out_lines)
        not out_stream_close or out_stream.close()

        included_files.clear()
        lc_files.clear()
        provided_modules.clear()
        required_modules.clear()
        provided_submodules.clear()
        required_submodules.clear()


def gen_lc_deps(src_name, lc_files):
    result = []
    if src_name and lc_files:
        result.append("{0}: {1}\n".format(src_name, " ".join(lc_files)))
    return result


def gen_include_deps(include_targets, src_name, included_files):
    result = []
    targets = " ".join(filter(None, include_targets))
    if targets:
        prereqs = " ".join(filter(None, [src_name] + list(included_files)))
        if prereqs:
            result.append("{0}: {1}\n".format(targets, prereqs))
    return result


def gen_module_deps(
    obj_name,
    mod_stamp_name,
    provided_modules,
    required_modules,
    provided_submodules,
    required_submodules,
    mod_dir,
    mod_upper,
    mod_ext,
    smod_infix,
    smod_ext,
):
    result = []

    mod_stamp_name = mod_stamp_name or obj_name
    if mod_stamp_name:
        targets = list(
            modules_to_filenames(provided_modules, mod_dir, mod_upper, mod_ext)
        )
        targets.extend(
            submodules_to_filenames(
                provided_submodules, mod_dir, mod_upper, smod_infix, smod_ext
            )
        )

        if targets:
            result.append(
                "{0}: {1}\n".format(" ".join(targets), mod_stamp_name)
            )

        prereqs = list(
            modules_to_filenames(
                # Do not depend on the provided modules:
                [m for m in required_modules if m not in provided_modules],
                mod_dir,
                mod_upper,
                mod_ext,
            )
        )

        prereqs.extend(
            submodules_to_filenames(
                # Do not depend on the provided submodules:
                [
                    m
                    for m in required_submodules
                    if m not in provided_submodules
                ],
                mod_dir,
                mod_upper,
                smod_infix,
                smod_ext,
            )
        )

        if prereqs:
            result.append(
                "{0}: {1}\n".format(
                    " ".join(filter(None, (obj_name, mod_stamp_name))),
                    " ".join(prereqs),
                )
            )

    if obj_name:
        targets = list(
            modules_to_filenames(
                set(
                    module
                    for module, _ in provided_submodules
                    # Do not depend on the provided modules:
                    if module not in provided_modules
                ),
                mod_dir,
                mod_upper,
                mod_ext,
            )
        )

        if mod_stamp_name != obj_name:
            targets.append(mod_stamp_name)

        if targets:
            result.append(
                "{0}: #-hint {1}\n".format(" ".join(targets), obj_name)
            )

    return result


def modules_to_filenames(modules, directory, upprecase, extension):
    result = modules
    if upprecase:
        result = map(lambda s: s.upper(), result)
    if directory:
        result = map(lambda s: os.path.join(directory, s), result)
    if extension:
        result = map(lambda s: "{0}.{1}".format(s, extension), result)
    return result


def submodules_to_filenames(submodules, directory, uppercase, infix, extension):
    result = modules_to_filenames(
        map(
            lambda module_submodule: (
                "{1}{0}{2}".format(infix, *module_submodule)
                if infix and module_submodule[1]
                else module_submodule[1] or module_submodule[0]
            ),
            submodules,
        ),
        directory,
        uppercase,
        extension,
    )
    return result


if __name__ == "__main__":
    main()
