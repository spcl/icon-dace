#!/usr/bin/env python
import os
from typing import Any, Collection, Dict, List, Optional, Set, Tuple, Union

from functools import singledispatch
from pathlib import Path
import re
import argparse
import logging

import numpy as np

from yaml import load as load_yaml

try:
    from yaml import CLoader as YAML_Loader
except ImportError:
    from yaml import Loader as YAML_Loader

import dace

# FIXME(small): maybe use `contiguous` attribute  to better ensure C interoperatiblity
# FIXME: Analysis pass to find out which struct members are actually read & written by SDFG
# FIXME: Need free function (for verification deep copies & uncached dace shallow copies)!
# FIXME: add caching of struct & struct array arguments
# FIXME: have to see about copy back from dace structs/arrays to fortran
# FIXME: check if arrays are allocated/associated etc when copying in

_PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE = {
    dace.float32: "real(kind=c_float)",
    dace.float64: "real(kind=c_double)",
    dace.int32: "integer(kind=c_int)",
    dace.int64: "integer(kind=c_long)",
    dace.int8: "character(kind=c_char)",
}


def join_wrapped(
    elems: Collection[str], seperator: str, before: str, after: str
) -> str:
    source = seperator.join(elems)
    if 0 < len(elems):
        source = before + source + after
    return source


@singledispatch
def dace_type_to_fortran_c_var_type_decl(dace_type) -> str:
    raise NotImplementedError(
        f"Unable to convert dace type '{dace_type}' to type declaration for a "
        "fortran variable (with C interoperatiblity)"
    )


@dace_type_to_fortran_c_var_type_decl.register
def _(typeclass: dace.typeclass) -> str:
    return _PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE[typeclass]


@dace_type_to_fortran_c_var_type_decl.register
def _(scalar: dace.data.Scalar) -> str:
    return _PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE[scalar.dtype]


@dace_type_to_fortran_c_var_type_decl.register(dace.data.Array)
@dace_type_to_fortran_c_var_type_decl.register(dace.data.Structure)
def _(_) -> str:
    return "type(c_ptr)"


script_dir = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(script_dir, "../unused_names.txt")
unused_names = set()

# Read the file line by line and add each line to the set
with open(input_file, "r") as file:
    for line in file:
        # Strip any leading/trailing whitespace and add the line to the set
        unused_names.add(line.strip())

print(unused_names)


class ImportsCollector:

    required_symbols: Dict[str, None]
    initializaion_checks_ignore_list: Set[str]

    def __init__(self, initializaion_checks_ignore_list: Set[str]):
        self.initializaion_checks_ignore_list = initializaion_checks_ignore_list
        self.required_symbols = {}

    def require_symbol(self, name: str) -> None:
        self.required_symbols[name] = None

    def generate_imports_str(self, module_definitions: Dict[str, str]) -> str:
        imports_dict: Dict[str, List[str]] = {}

        missing_module_definitions = (
            set(self.required_symbols.keys())
            - set(module_definitions.keys())
            - self.initializaion_checks_ignore_list
        )
        if 0 < len(missing_module_definitions):
            logging.warning(
                "Failed to find modules for symbols:"
                + "".join(
                    f"\n- {symbol}" for symbol in sorted(missing_module_definitions)
                )
            )

        for symbol in self.required_symbols:
            if symbol in missing_module_definitions:
                continue
            imports_dict.setdefault(module_definitions[symbol], []).append(symbol)

        return "\n".join(
            f"  use {module_name}, only: &\n"
            + ", &\n".join(f"    {symbol}" for symbol in symbols)
            for module_name, symbols in imports_dict.items()
        )


@singledispatch
def dace_type_to_fortran_rich_var_type_decl(
    dace_type, is_local: bool = False, enable_inout_hack: bool = False
) -> str:
    raise NotImplementedError(
        f"Unable to convert dace type '{dace_type}' to type declaration for a fortran variable"
    )


@dace_type_to_fortran_rich_var_type_decl.register
def _(
    typeclass: dace.typeclass, is_local: bool = False, enable_inout_hack: bool = False
) -> str:
    return _PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE[typeclass]


@dace_type_to_fortran_rich_var_type_decl.register
def _(
    scalar: dace.data.Scalar, is_local: bool = False, enable_inout_hack: bool = False
) -> str:
    return _PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE[scalar.dtype]


@dace_type_to_fortran_rich_var_type_decl.register
def _(
    desc: dace.data.ContainerArray,
    is_local: bool = False,
    enable_inout_hack: bool = False,
) -> str:
    pointer_attr = "target" if not is_local else "pointer"

    shape_str = ",".join(":" * len(desc.shape))
    return f"type({desc.stype.name}), dimension({shape_str}), {pointer_attr}"


@dace_type_to_fortran_rich_var_type_decl.register
def _(
    desc: dace.data.Array, is_local: bool = False, enable_inout_hack: bool = False
) -> str:
    pointer_attr = "target" if not is_local else "pointer"

    if enable_inout_hack and desc.shape == (1,):
        # FIXME: hack for ```INOUT```/```OUT``` parameters
        return f"{_PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE[desc.dtype]}, {pointer_attr}"

    shape_str = ",".join(":" * len(desc.shape))
    return f"{_PRIMITIVE_DACE_TYPE_TO_FORRTRAN_C_TYPE[desc.dtype]}, dimension({shape_str}), {pointer_attr}"


@dace_type_to_fortran_rich_var_type_decl.register
def _(
    desc: dace.data.Structure, is_local: bool = False, enable_inout_hack: bool = False
) -> str:
    pointer_attr = "target" if not is_local else "pointer"
    return f"type({desc.name}), {pointer_attr}"


def collect_structs_and_arrays(
    desc: dace.data.Data,
    structs: Dict[dace.data.Structure, None],
    arrays: Dict[dace.data.Array, None],
    imports_collector: ImportsCollector,
) -> None:

    if isinstance(desc, dace.data.ContainerArray):
        if desc in arrays:
            return
        arrays[desc] = None
        collect_structs_and_arrays(desc.stype, structs, arrays, imports_collector)

    elif isinstance(desc, dace.data.Array):
        arrays[desc] = None

    elif isinstance(desc, dace.data.Structure):
        if desc in structs:
            return
        structs[desc] = None
        imports_collector.require_symbol(desc.name)
        for member in desc.members.values():
            collect_structs_and_arrays(member, structs, arrays, imports_collector)


def generate_dace_struct_definition(struct: dace.data.Structure) -> str:

    struct_members_src = ""

    for name, desc in sorted(struct.members.items()):
        name = fix_identifier(name)
        struct_members_src += (
            f"    {dace_type_to_fortran_c_var_type_decl(desc)} :: {name}\n"
        )

    return f"""\
  type, bind(c) :: dace_{struct.name}
{struct_members_src}
  end type dace_{struct.name}
"""


def fix_identifier(identifier: str) -> str:
    while identifier.startswith("_"):
        identifier = identifier[1:]
    return identifier


_F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX = fix_identifier("__f2dace_A_")
_F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX = fix_identifier("__f2dace_OA_")
_F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX = fix_identifier("__f2dace_SA_")
_F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX = fix_identifier("__f2dace_SOA_")

_IDENTIFIER_REGEX = "\\w+"
_NUMBER_REGEX = "\\d+"
_F2DACE_PARAM_ARRAY_SIZE_HELPER_PATTERN_STR = (
    f"{_F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX}"
    f"(?P<array_name>{_IDENTIFIER_REGEX})"
    "_d_"
    f"(?P<dim_num>{_NUMBER_REGEX})"
    f"_s_{_NUMBER_REGEX}"
)
_F2DACE_PARAM_ARRAY_OFFSET_HELPER_PATTERN_STR = (
    f"{_F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX}"
    f"(?P<array_name>{_IDENTIFIER_REGEX})"
    "_d_"
    f"(?P<dim_num>{_NUMBER_REGEX})"
    f"_s_{_NUMBER_REGEX}"
)
_F2DACE_STRUCT_ARRAY_SIZE_HELPER_PATTERN_STR = (
    f"{_F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX}"
    f"(?P<array_name>{_IDENTIFIER_REGEX})"
    "_d_"
    f"(?P<dim_num>{_NUMBER_REGEX})"
    f"_s_{_NUMBER_REGEX}"
)
_F2DACE_STRUCT_ARRAY_OFFSET_HELPER_PATTERN_STR = (
    f"{_F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX}"
    f"(?P<array_name>{_IDENTIFIER_REGEX})"
    "_d_"
    f"(?P<dim_num>{_NUMBER_REGEX})"
    f"_s_{_NUMBER_REGEX}"
)


def extract_array_helper_information(name: str, pattern_str: str) -> Tuple[str, int]:
    match = re.fullmatch(pattern_str, name)
    assert (
        match
    ), f"Failed to extract array name and dim number from helper array ('{name}')"

    return match.group("array_name"), int(match.group("dim_num"))


# Error: 'cells_plwa_verts' at (1) is not a member of the 't_int_state' structure; did you mean 'cells_aw_verts'?
# Error: 'geofac_qdiv' at (1) is not a member of the 't_int_state' structure; did you mean 'geofac_div'?
# Error: 'cn_e' at (1) is not a member of the 't_grid_edges' structure; did you mean 'fn_e'?
# Error: 'cz_c' at (1) is not a member of the 't_grid_cells' structure; did you mean 'f_c'?
# Error: 'ddt_ua_cen' at (1) is not a member of the 't_nh_diag' structure; did you mean 'ddt_ua_dyn'?
# Error: 'ddt_ua_cen_is_associated' at (1) is not a member of the 't_nh_diag' structure; did you mean 'ddt_ua_dyn_is_associated'?
# Error: 'ddt_va_cen' at (1) is not a member of the 't_nh_diag' structure; did you mean 'ddt_va_dyn'?
# Error: 'ddt_va_cen_is_associated' at (1) is not a member of the 't_nh_diag' structure; did you mean 'ddt_va_dyn_is_associated'?
# Error: 'ddt_vn_cen' at (1) is not a member of the 't_nh_diag' structure; did you mean 'ddt_vn_dyn'?
# Error: 'ddt_vn_cen_is_associated' at (1) is not a member of the 't_nh_diag' structure; did you mean 'ddt_vn_dyn_is_associated'?
# Error: 'vor_u' at (1) is not a member of the 't_nh_diag' structure; did you mean 'vor'?
# Error: 'vor_v' at (1) is not a member of the 't_nh_diag' structure; did you mean 'vor'?
# Error: 'deepatmo_t1ifc' at (1) is not a member of the 't_nh_metrics' structure; did you mean 'deepatmo_vol_mc'?
# Error: 'deepatmo_t1mc' at (1) is not a member of the 't_nh_metrics' structure; did you mean 'deepatmo_vol_mc'?
# Error: 'deepatmo_t2mc' at (1) is not a member of the 't_nh_metrics' structure; did you mean 'deepatmo_vol_mc'?
# Error: 'dzgpot_mc' at (1) is not a member of the 't_nh_metrics' structure; did you mean 'dgeopot_mc'?
# Error: 'fbk_dom_volume' at (1) is not a member of the 't_nh_metrics' structure
# Error: 'zgpot_ifc' at (1) is not a member of the 't_nh_metrics' structure; did you mean 'z_ifc'?
# Error: 'zgpot_mc' at (1) is not a member of the 't_nh_metrics' structure; did you mean 'dgeopot_mc'?


_STRUCT_MEMBER_TYPES_IGNORE_LIST = {
    dace.data.Scalar(dace.int8),  # this likely was a string, so we ignore it
}


# TODO: better checks that ```malloc``` memory is suitable for DaCe data descriptors
# e.g.,: alignment, padding, strides, etc
# This is relevant for the generation of copy in functions, but also the generation of copy in expressions
def generate_copy_in_function_struct(
    struct: dace.data.Structure,
    struct_ignore_list: Set[str],
    struct_members_use_null: Dict[str, Set[str]],
) -> str:
    global unused_names
    dace_type_name = f"dace_{struct.name}"

    members_use_null = struct_members_use_null.get(struct.name, set())
    copy_fields_src = ""
    for member_name, member_type in struct.members.items():

        if member_type in _STRUCT_MEMBER_TYPES_IGNORE_LIST:
            logging.warning(
                f"Ignoring member initialization '{member_name}' in struct '{struct.name}' "
                f"because type '{member_type}' is on the ignore list."
            )
            continue

        member_name = fix_identifier(member_name)

        if member_name not in unused_names:
            if member_name in members_use_null:
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = c_null_ptr
    """
            elif member_name.startswith(_F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX):
                array_name, dim_num = extract_array_helper_information(
                    member_name, _F2DACE_STRUCT_ARRAY_SIZE_HELPER_PATTERN_STR
                )
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = size(fortran_obj%{array_name}, dim={dim_num + 1})
    """

            elif member_name.startswith(
                _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX
            ):
                array_name, dim_num = extract_array_helper_information(
                    member_name, _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_PATTERN_STR
                )
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = lbound(fortran_obj%{array_name}, dim={dim_num + 1})
    """

            else:
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = {
            generate_copy_in_fortran_expr(
                member_type,
                expr=f"fortran_obj%{member_name}",
                steal_arrays_expr="steal_arrays",
                minimal_structs_expr="minimal_structs",
            )
    }
    """
        else:
            if member_name in members_use_null:
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = c_null_ptr
    """
            elif member_name.startswith(_F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX):
                array_name, dim_num = extract_array_helper_information(
                    member_name, _F2DACE_STRUCT_ARRAY_SIZE_HELPER_PATTERN_STR
                )
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = 0
    """

            elif member_name.startswith(
                _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX
            ):
                array_name, dim_num = extract_array_helper_information(
                    member_name, _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_PATTERN_STR
                )
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = 0
    """

            else:
                copy_fields_src += f"""\
        dace_rich_obj%{member_name} = c_null_ptr
    """

    if struct.name in struct_ignore_list:
        logging.warning(
            f"Disabling copy in of struct '{struct.name}' because it's on the ignore list."
        )
        copy_fields_src = ""

    return f"""\
  function copy_in_{struct.name}(fortran_obj, steal_arrays, minimal_structs) result(dace_obj_ptr)
    {dace_type_to_fortran_rich_var_type_decl(struct)} :: fortran_obj
    logical :: steal_arrays, minimal_structs
    type(c_ptr) :: dace_obj_ptr
    type({dace_type_name}), pointer :: dace_rich_obj

    dace_obj_ptr = malloc(c_sizeof(dace_rich_obj))
    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

{copy_fields_src}
  end function copy_in_{struct.name}

"""


class ArrayLoopHelper:

    @staticmethod
    def dimensions(rank: int) -> str:
        return ",".join(rank * ":")

    @staticmethod
    def indices_decl(rank: int, prefix: str = "i") -> str:
        return "integer :: " + ", ".join(f"{prefix}{i}" for i in range(rank))

    @staticmethod
    def indices_expr(rank: int, prefix: str = "i") -> str:
        return ", ".join(f"{prefix}{i}" for i in range(rank))

    @staticmethod
    def indices_copy_stmt(
        rank: int, lhs_prefix: str, rhs_prefix: str, rhs_suffix: str = ""
    ) -> str:
        return "; ".join(
            f"{lhs_prefix}{i} = {rhs_prefix}{i}{rhs_suffix}" for i in range(rank)
        )

    @staticmethod
    def loop_begins(rank: int, array_name: str, prefix: str = "i") -> str:
        loop_begins = ""

        for i in range(rank):
            loop_begins += f"""
    {i * '  '}do {prefix}{i} = 1, size({array_name}, dim={i + 1})"""

        return loop_begins

    @staticmethod
    def loop_ends(rank: int) -> str:
        loop_ends = ""

        for i in range(rank):
            loop_ends += f"""\
    {(rank - i - 1) * '  '}end do
"""

        return loop_ends


def generate_copy_in_function_array(
    array: Union[dace.data.ContainerArray, dace.data.Array]
) -> Tuple[str, str]:
    rank = len(array.shape)

    # FIXME(important!): take offsets into consideration!

    copy_in_str = ""
    copy_in_interface_str = ""
    if isinstance(array, dace.data.ContainerArray):
        base_type = array.stype
        base_name = base_type.name
        # cannot take pointer of struct arrays
        array_shallow_copy_str = ""
    else:
        assert isinstance(array, dace.data.Array)
        base_type = array.dtype
        base_name = base_type.to_string()
        array_shallow_copy_str = """
    if (steal_arrays .eqv. .true.) then
      dace_array_ptr = c_loc(fortran_array)
      return
    end if
"""
        # TODO: quite hacky to implement this here
        if array.dtype == dace.int32:
            # need to add fix because F2DaCe uses `int32` for Fortran `logical`
            copy_in_interface_str += f"""
interface logical_fix_{rank}d
  module procedure logical_to_int_{rank}d
  module procedure int_to_int_{rank}d
end interface logical_fix_{rank}d
"""

            copy_in_str += f"""\
  function logical_to_int_{rank}d(inp) result(out)
    logical(4), dimension({ArrayLoopHelper.dimensions(rank)}), target :: inp
    integer(kind=c_int), dimension({ArrayLoopHelper.dimensions(rank)}), pointer :: out

    call c_f_pointer(c_loc(inp), out, shape=shape(inp))
  end function logical_to_int_{rank}d

  function int_to_int_{rank}d(inp) result(out)
    integer(kind=c_int), dimension({ArrayLoopHelper.dimensions(rank)}), target :: inp
    integer(kind=c_int), dimension({ArrayLoopHelper.dimensions(rank)}), pointer :: out

    out => inp
  end function int_to_int_{rank}d

"""

    # FIXME: type declarations here are super messy!
    copy_in_str += f"""\
  function copy_in_{base_name}_{rank}d_array(fortran_array, steal_arrays, minimal_structs) result(dace_array_ptr)
    {dace_type_to_fortran_rich_var_type_decl(array)} :: fortran_array
    logical :: steal_arrays, minimal_structs
    {dace_type_to_fortran_c_var_type_decl(array)} :: dace_array_ptr
    {dace_type_to_fortran_c_var_type_decl(base_type)}, dimension({ArrayLoopHelper.dimensions(rank)}), pointer :: dace_rich_array

    {ArrayLoopHelper.indices_decl(rank)}

    if (.not. c_associated(c_loc(fortran_array))) then
      dace_array_ptr = c_null_ptr
      return
    end if

{array_shallow_copy_str}
    dace_array_ptr = malloc(c_sizeof(dace_array_ptr) * size(fortran_array))
    call c_f_pointer(dace_array_ptr, dace_rich_array, shape=shape(fortran_array))
{ArrayLoopHelper.loop_begins(rank, "fortran_array")}
    {rank * "  "}dace_rich_array({ArrayLoopHelper.indices_expr(rank)}) = {
        generate_copy_in_fortran_expr(
            base_type,
            expr=f"fortran_array({ArrayLoopHelper.indices_expr(rank)})",
            steal_arrays_expr="steal_arrays",
            minimal_structs_expr="minimal_structs",
        )
}
{ArrayLoopHelper.loop_ends(rank)}
  end function copy_in_{base_name}_{rank}d_array

"""

    return copy_in_str, copy_in_interface_str


@singledispatch
def generate_copy_in_fortran_expr(
    dace_type: dace.data.Data,
    expr: str,
    steal_arrays_expr: str,
    minimal_structs_expr: str,
    enable_inout_hack: bool = False,
) -> str:
    raise NotImplementedError(f"Unable to copy in fortran expression to '{dace_type}'")


_PRIMITIVE_FORTRAN_TO_DACE_COPY_IN_FUNCTIONS: Dict[
    dace.dtypes.typeclass, Optional[str]
] = {
    dace.float32: None,
    dace.float64: None,
    dace.int32: None,
    dace.int64: None,
    # dace.int8: None,
    # dace.bool: None,
    # dace.bool_: None,
}


@generate_copy_in_fortran_expr.register
def _(
    dtype: dace.dtypes.typeclass,
    expr: str,
    steal_arrays_expr: str,
    minimal_structs_expr: str,
    enable_inout_hack: bool = False,
) -> str:
    copy_in_func = _PRIMITIVE_FORTRAN_TO_DACE_COPY_IN_FUNCTIONS[dtype]
    if copy_in_func is not None:
        return f"{copy_in_func}({expr})"
    return expr


@generate_copy_in_fortran_expr.register
def _(
    scalar: dace.data.Scalar,
    expr: str,
    steal_arrays_expr: str,
    minimal_structs_expr: str,
    enable_inout_hack: bool = False,
) -> str:
    return generate_copy_in_fortran_expr(
        scalar.dtype,
        expr,
        steal_arrays_expr=steal_arrays_expr,
        minimal_structs_expr=minimal_structs_expr,
        enable_inout_hack=enable_inout_hack,
    )


@generate_copy_in_fortran_expr.register
def _(
    array: dace.data.Array,
    expr: str,
    steal_arrays_expr: str,
    minimal_structs_expr: str,
    enable_inout_hack: bool = False,
) -> str:
    if enable_inout_hack and array.shape == (1,):
        # FIXME: this really doesn't work for verification!
        return f"c_loc({expr})"

    rank = len(array.shape)
    if array.dtype == dace.int32:
        # could have been `logical` in fortran, might need cast
        expr = f"logical_fix_{rank}d({expr})"

    return f"""copy_in_{array.dtype.to_string()}_{rank}d_array( &
    fortran_array={expr}, &
    steal_arrays={steal_arrays_expr}, &
    minimal_structs={minimal_structs_expr} &
  )"""


@generate_copy_in_fortran_expr.register
def _(
    struct: dace.data.Structure,
    expr: str,
    steal_arrays_expr: str,
    minimal_structs_expr: str,
    enable_inout_hack: bool = False,
) -> str:
    return f"""copy_in_{struct.name}( &
    fortran_obj={expr}, &
    steal_arrays={steal_arrays_expr}, &
    minimal_structs={minimal_structs_expr} &
  )"""


@generate_copy_in_fortran_expr.register
def _(
    struct_array: dace.data.ContainerArray,
    expr: str,
    steal_arrays_expr: str,
    minimal_structs_expr: str,
    enable_inout_hack: bool = False,
) -> str:
    return f"""copy_in_{struct_array.stype.name}_{len(struct_array.shape)}d_array( &
    fortran_array={expr}, &
    steal_arrays={steal_arrays_expr}, &
    minimal_structs={minimal_structs_expr} &
  )"""


def generate_copy_back_subroutine_struct(
    struct: dace.data.Structure,
    struct_ignore_list: Set[str],
    struct_members_use_null: Dict[str, Set[str]],
) -> str:
    global unused_names
    copy_back_fields_src = ""

    if struct.name not in struct_ignore_list:

        members_use_null = struct_members_use_null.get(struct.name, set())
        for member_name, member_type in struct.members.items():
            member_name = fix_identifier(member_name)

            if (
                (
                    # only copy back scalar, struct & struct arays
                    isinstance(member_type, dace.data.Array)
                    and not isinstance(member_type, dace.data.ContainerArray)
                )
                or member_type in _STRUCT_MEMBER_TYPES_IGNORE_LIST
                or member_name in members_use_null
                or member_name.startswith(_F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX)
                or member_name.startswith(
                    _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX
                )
            ):
                continue

            assert (
                isinstance(member_type, dace.data.Scalar)
                or isinstance(member_type, dace.data.Structure)
                or (
                    isinstance(member_type, dace.data.ContainerArray)
                    and isinstance(member_type.stype, dace.data.Structure)
                )
            )

            if member_name not in unused_names:
                copy_back_fields_src += generate_copy_back_stmts(
                    member_type,
                    f"fortran_obj%{member_name}",
                    f"dace_rich_obj%{member_name}",
                )

    return f"""\
  subroutine copy_back_{struct.name}(fortran_obj, dace_obj_ptr)
    {dace_type_to_fortran_rich_var_type_decl(struct)} :: fortran_obj
    type(c_ptr) :: dace_obj_ptr

    type(dace_{struct.name}), pointer :: dace_rich_obj

    if (.not. c_associated(c_loc(fortran_obj))) then
      if (c_associated(dace_obj_ptr)) then
        print *, "copy_back_{struct.name}: Invalid allocation of {struct.name} by DaCe!"
      end if
      return
    end if

    call c_f_pointer(dace_obj_ptr, dace_rich_obj)

{copy_back_fields_src}

  end subroutine copy_back_{struct.name}
"""


def generate_copy_back_subroutine_struct_array(
    struct_array: dace.data.ContainerArray,
) -> str:
    rank = len(struct_array.shape)
    stype = struct_array.stype
    assert isinstance(stype, dace.data.Structure)

    return f"""
  subroutine copy_back_{stype.name}_{rank}d_array(fortran_struct_array, dace_struct_array_ptr)
    {dace_type_to_fortran_rich_var_type_decl(struct_array)}, intent(in) :: fortran_struct_array
    {dace_type_to_fortran_c_var_type_decl(struct_array)}, intent(in) :: dace_struct_array_ptr

    {ArrayLoopHelper.indices_decl(rank)}
    {dace_type_to_fortran_c_var_type_decl(stype)}, dimension({ArrayLoopHelper.dimensions(rank)}), pointer :: dace_struct_array_rich

    if (.not. c_associated(c_loc(fortran_struct_array))) then
      if (c_associated(dace_struct_array_ptr)) then
        print *, "copy_back_{stype.name}_{rank}d_array: Invalid allocation of {stype.name} array by DaCe!" 
      end if
      return
    end if

    call c_f_pointer(dace_struct_array_ptr, dace_struct_array_rich, shape=shape(fortran_struct_array))

{ArrayLoopHelper.loop_begins(rank, "fortran_struct_array")}
{generate_copy_back_stmts(
    stype,
    f"fortran_struct_array({ArrayLoopHelper.indices_expr(rank)})",
    f"dace_struct_array_rich({ArrayLoopHelper.indices_expr(rank)})",
)}
{ArrayLoopHelper.loop_ends(rank)}
  end subroutine copy_back_{stype.name}_{rank}d_array
"""


@singledispatch
def generate_copy_back_stmts(dace_type, fortran_expr: str, dace_expr: str) -> str:
    raise NotImplementedError(
        f"Unable to generate copy back statements for dace type {dace_type}"
    )


@generate_copy_back_stmts.register
def generate_copy_back_stmts_scalar(
    scalar: dace.data.Scalar, fortran_expr: str, dace_expr: str
) -> str:
    if scalar.dtype == dace.int32:
        # could have been `logical` in fortran, might need cast
        return f"""\
    {fortran_expr} = transfer({dace_expr}, mold={fortran_expr})
"""
    return f"""\
    {fortran_expr} = {dace_expr}
"""


@generate_copy_back_stmts.register
def generate_copy_back_stmts_struct(
    struct: dace.data.Structure, fortran_expr: str, dace_expr: str
) -> str:
    return f"""\
    call copy_back_{struct.name}({fortran_expr}, {dace_expr})
"""


@generate_copy_back_stmts.register
def generate_copy_back_stmts_struct_array(
    struct_array: dace.data.ContainerArray, fortran_expr: str, dace_expr: str
) -> str:
    rank = len(struct_array.shape)
    stype = struct_array.stype
    assert isinstance(stype, dace.data.Structure)
    return f"""\
    call copy_back_{stype.name}_{rank}d_array({fortran_expr}, {dace_expr})
"""


_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX = fix_identifier("__f2dace_OPTIONAL_")
_OPTIONAL_PROXY_PREFIX = fix_identifier("__DACE_OPT_PROXY_")


# FIXME(THRESHOLDS!): choose appropriate thresholds!
_COMPARISON_DEFAULT_THRESHOLDS_STR = """
  real(8), parameter :: float32_default_rel_threshold = 1.0e-8
  real(8), parameter :: float32_default_abs_threshold = 0.0

  real(8), parameter :: float64_default_rel_threshold = 1.0e-12
  real(8), parameter :: float64_default_abs_threshold = 0.0

  real(8), parameter :: int32_default_rel_threshold = 1.0e-8
  real(8), parameter :: int32_default_abs_threshold = 0.0

  real(8), parameter :: int64_default_rel_threshold = 1.0e-12
  real(8), parameter :: int64_default_abs_threshold = 0.0
"""
_COMPARISON_PRIMITIVE_FUNCTIONS_STR = "".join(
    f"""

  subroutine compare_{dtype.to_string()}_scalar( &
    actual, &
    ref, &
    result, &
    rel_threshold, &
    abs_threshold, &
    scalar_expr &
  )
    {fortran_type_str}, intent(in) :: actual, ref
    logical, intent(out) :: result
    real(8), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in), optional :: scalar_expr

    real(8) :: rel_error, abs_error, threshold_ratio
    real(8) :: actual_rel_threshold, actual_abs_threshold
    CHARACTER(len=5000) :: message_text = ''

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = {dtype.to_string()}_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = {dtype.to_string()}_default_abs_threshold
    end if

    result = abs(ref - actual) <= max(actual_rel_threshold * abs(ref), actual_abs_threshold)

    if (present(scalar_expr) .and. .not. result) then

      threshold_ratio = real(abs(ref - actual), kind=8) / max(actual_rel_threshold * abs(ref), actual_abs_threshold)
      rel_error = abs(real(ref - actual, kind=8)/ref)
      abs_error = abs(ref - actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,e28.20,a,e28.20)') &
        "Verification failed for scalar '", &
          trim(scalar_expr), &
        "'"//char(10)//"    - rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
          ", threshold_ratio = ", &
          threshold_ratio, &
        char(10)//"    - ref = ", &
          ref, &
          ", actual = ", &
          actual
        print *, "compare_{dtype.to_string()}_scalar"
        print *, message_text

    end if

  end subroutine compare_{dtype.to_string()}_scalar
"""
    for dtype, fortran_type_str in (
        (dace.float32, "real(kind=c_float)"),
        (dace.float64, "real(kind=c_double)"),
        (dace.int32, "integer(kind=c_int)"),
        (dace.int64, "integer(kind=c_long)"),
    )
)


@singledispatch
def generate_array_comparison_subroutine(dace_type) -> str:
    raise NotImplementedError(
        f"Unable to generate comparison subroutine for dace type {dace_type}"
    )


@generate_array_comparison_subroutine.register
def _(array: dace.data.Array) -> str:
    # FIXME: can we check sizes of Fortran & SDFG arrays?!
    rank = len(array.shape)
    dtype = array.dtype

    comparison_stmts = generate_comparison_check_stmts(
        dtype,
        actual_expr=f"actual_rich({ArrayLoopHelper.indices_expr(rank)})",
        ref_expr=f"ref({ArrayLoopHelper.indices_expr(rank)})",
        result_expr="local_result",
        rel_threshold_expr="actual_rel_threshold",
        abs_threshold_expr="actual_abs_threshold",
    )

    ll = len(ArrayLoopHelper.indices_expr(rank, "max_threshold_ratio_i").split(","))
    assert ll == rank

    return f"""
  subroutine compare_{dtype.to_string()}_{rank}d_array( &
    actual, &
    ref, &
    result, &
    array_expr, &
    rel_threshold, &
    abs_threshold &
  )
    {dace_type_to_fortran_c_var_type_decl(array)} :: actual
    {dace_type_to_fortran_rich_var_type_decl(array)}, intent(in) :: ref
    logical, intent(out) :: result
    real(kind=c_double), intent(in), optional :: rel_threshold, abs_threshold
    character(*), intent(in) :: array_expr

    real(kind=c_double) :: actual_rel_threshold, actual_abs_threshold
    logical :: local_result
    {ArrayLoopHelper.indices_decl(rank)}
    {dace_type_to_fortran_c_var_type_decl(dtype)}, dimension({ArrayLoopHelper.dimensions(rank)}), pointer :: actual_rich
    CHARACTER(len=5000) :: message_text = ''

    {dace_type_to_fortran_c_var_type_decl(dtype)} :: error_ref, error_actual
    integer, dimension(0:{rank - 1}) :: max_threshold_ratio_loc
    real(8) :: rel_error, abs_error, threshold_ratio

    {ArrayLoopHelper.indices_decl(rank, "max_threshold_ratio_i")}
    {ArrayLoopHelper.indices_decl(rank, "first_fail_i")}
    {ArrayLoopHelper.indices_decl(rank, "last_fail_i")}
    {ArrayLoopHelper.indices_decl(rank, "dim_i")}
    integer :: total_fails
    integer :: total_indices
    integer :: first_fail

    {"\n    ".join([f"{ex} = -1" for ex in ArrayLoopHelper.indices_expr(rank, "first_fail_i").split(", ")])}
    {"\n    ".join([f"{ex} = -1" for ex in ArrayLoopHelper.indices_expr(rank, "last_fail_i").split(", ")])}
    {"\n    ".join([f"{ex} = size(ref, dim={kk+1})" for kk, ex in enumerate(ArrayLoopHelper.indices_expr(rank, "dim_i").split(", "))])}
    total_fails = 0
    total_indices = 0
    first_fail = -1
    
    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_{dtype.to_string()}_{rank}d_array"
        print *, message_text
      end if

      return
    end if

    result = .true.

    if (present(rel_threshold)) then
      actual_rel_threshold = rel_threshold
    else
      actual_rel_threshold = {dtype.to_string()}_default_rel_threshold
    end if

    if (present(abs_threshold)) then
      actual_abs_threshold = abs_threshold
    else
      actual_abs_threshold = {dtype.to_string()}_default_abs_threshold
    end if

    call c_f_pointer(actual, actual_rich, shape=shape(ref))


{ArrayLoopHelper.loop_begins(rank, "ref")}
{comparison_stmts}
    result = result .and. local_result
{ArrayLoopHelper.loop_ends(rank)}

{ArrayLoopHelper.loop_begins(rank, "ref")}
    {comparison_stmts}
    if (.not. local_result) then
        if (first_fail == -1) then
            {"\n            ".join([f"first_fail_{i} = {i}" for i in ArrayLoopHelper.indices_expr(rank).split(", ")])}
            first_fail = 1
        endif
        {"\n        ".join([f"last_fail_{i} = {i}" for i in ArrayLoopHelper.indices_expr(rank).split(", ")])}
        total_fails = total_fails + 1
    endif
{ArrayLoopHelper.loop_ends(rank)}

    total_indices =  {" * ".join([f"size(ref, dim={i+1})" for i in range(rank)])}

    if (.not. result) then
      max_threshold_ratio_loc = maxloc(abs(ref - actual_rich) / max(actual_rel_threshold * abs(ref), actual_abs_threshold))
      {ArrayLoopHelper.indices_copy_stmt(rank, "max_threshold_ratio_i", "max_threshold_ratio_loc(", ")")}

      error_ref = ref({ArrayLoopHelper.indices_expr(rank, "max_threshold_ratio_i")})
      error_actual = actual_rich({ArrayLoopHelper.indices_expr(rank, "max_threshold_ratio_i")})

      threshold_ratio = real(abs(error_ref - error_actual), kind=8) / max(actual_rel_threshold * abs(error_ref), actual_abs_threshold)
      rel_error = abs(real(error_ref - error_actual, kind=8)/error_ref)
      abs_error = abs(error_ref - error_actual)

      write (message_text, '(a,a,a,e28.20,a,e28.20,a,e28.20,a,"(",{',", ",'.join(["i0"]*rank)},")",a,e28.20,a,e28.20,a,"(",{',", ",'.join(["i0"]*rank)},")",a,"(",{',", ",'.join(["i0"]*rank)},")",a,i0,a,i0,a,i0,a,"(",{',", ",'.join(["i0"]*rank)},")")') &
        "Verification failed for array '", &
          trim(array_expr), &
        "':"//char(10)//"    - max_threshold_ratio = ", &
          threshold_ratio, &
          ", rel_error = ", &
          rel_error, &
          ", abs_error = ", &
          abs_error, &
        char(10)//"    - at (", &
          {ArrayLoopHelper.indices_expr(rank, "max_threshold_ratio_i")}, &
          "), ref = ", &
          error_ref, &
          ", actual = ", &
          error_actual, &
          " first_fail_index: ", {ArrayLoopHelper.indices_expr(rank, "first_fail_i")}, &
          " last_fail_index: ", {ArrayLoopHelper.indices_expr(rank, "last_fail_i")}, &
          " total_fails: ", total_fails, &
          " total_indices: ", total_indices, &
          " call_to_size: ", size(ref), &
          " shape: ", {ArrayLoopHelper.indices_expr(rank, "dim_i")}
      print *, "compare_{dtype.to_string()}_{rank}d_array"
      print *, message_text

    end if

    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_{dtype.to_string()}_{rank}d_array
"""


def generate_comparison_routine_dace_struct(
    struct: dace.data.Structure,
    struct_ignore_list: Set[str],
    struct_members_use_null: Dict[str, Set[str]],
) -> str:

    compare_fields_src = ""

    members_use_null = struct_members_use_null.get(struct.name, set())
    for member_name, member_type in struct.members.items():

        # TODO: refactor `_COPY_IN_IGNORES_STRUCT_MEMBER_TYPE`
        if member_type in _STRUCT_MEMBER_TYPES_IGNORE_LIST:
            continue
        if member_name in members_use_null:
            continue

        member_name = fix_identifier(member_name)

        if member_name.startswith(
            _F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX
        ) or member_name.startswith(_F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX):
            # we don't verify the helper fields
            continue

        compare_fields_src += f"""
    write (member_expr, '(a,a)') &
      trim(struct_expr), &
      "%{member_name}"
{generate_comparison_check_stmts(
    member_type,
    actual_expr=f"actual_rich%{member_name}",
    ref_expr=f"ref%{member_name}",
    result_expr="local_result",
    var_expr="member_expr",
)}
    result = result .and. local_result
"""

    if struct.name in struct_ignore_list:
        compare_fields_src = ""

    # FIXME: maybe structs need a `thresholds` argument
    return f"""
  subroutine compare_{struct.name}_struct( &
    actual, &
    ref, &
    result, &
    struct_expr &
  )
    {dace_type_to_fortran_c_var_type_decl(struct)}, intent(in) :: actual
    {dace_type_to_fortran_rich_var_type_decl(struct)}, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_expr

    CHARACTER(len=5000) :: member_expr = ''
    type(dace_{struct.name}), pointer :: actual_rich
    logical :: local_result
    call c_f_pointer(actual, actual_rich)

    result = .true.
{compare_fields_src}

    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_{struct.name}_struct
"""


@generate_array_comparison_subroutine.register
def _(struct_array: dace.data.ContainerArray) -> str:
    # FIXME: can we check sizes of Fortran & SDFG arrays?!
    rank = len(struct_array.shape)
    stype = struct_array.stype
    assert isinstance(stype, dace.data.Structure)

    comparison_stmts = generate_comparison_check_stmts(
        stype,
        actual_expr=f"actual_rich({ArrayLoopHelper.indices_expr(rank)})",
        ref_expr=f"ref({ArrayLoopHelper.indices_expr(rank)})",
        result_expr=f"local_result",
        var_expr="member_expr",
    )

    return f"""
  subroutine compare_{stype.name}_{rank}d_array( &
    actual, &
    ref, &
    result, &
    struct_array_expr &
  )
    {dace_type_to_fortran_c_var_type_decl(struct_array)}, intent(in) :: actual
    {dace_type_to_fortran_rich_var_type_decl(struct_array)}, intent(in) :: ref
    logical, intent(out) :: result
    character(*), intent(in) :: struct_array_expr

    CHARACTER(len=5000) :: member_expr = ''
    CHARACTER(len=5000) :: message_text = ''
    logical :: local_result
    {ArrayLoopHelper.indices_decl(rank)}
    {dace_type_to_fortran_c_var_type_decl(stype)}, dimension({ArrayLoopHelper.dimensions(rank)}), pointer :: actual_rich

    if (.not. c_associated(c_loc(ref))) then
      result = .not. c_associated(actual)

      if (.not. result) then
        write (message_text, '(a,a,a)') &
          "Verification failed for array '", &
            trim(struct_array_expr), &
          "':"//char(10)//"    - ref was NULL, but actual was not!"
        print *, "compare_{stype.name}_{rank}d_array"
        print *, message_text
      end if

      return
    end if

    result = .true.

    call c_f_pointer(actual, actual_rich, shape=shape(ref))

{ArrayLoopHelper.loop_begins(rank, "ref")}
    write (member_expr, '(a,a,{rank}(i0:,", "),a)') &
      trim(struct_array_expr), &
      "(", &
      [{ArrayLoopHelper.indices_expr(rank)}], &
      ")"
{comparison_stmts}
    result = result .and. local_result
{ArrayLoopHelper.loop_ends(rank)}

    ! FIXME: should be separate
    call free(actual)

  end subroutine compare_{stype.name}_{rank}d_array
"""


@singledispatch
def generate_comparison_check_stmts(
    dace_type,
    actual_expr: str,
    ref_expr: str,
    result_expr: str,
    var_expr: Optional[str] = None,
) -> str:
    raise NotImplementedError(
        f"Unable to generate comparison statements for dace type {dace_type}"
    )


@generate_comparison_check_stmts.register
def _(
    dtype: dace.dtypes.typeclass,
    actual_expr: str,
    ref_expr: str,
    result_expr: str,
    var_expr: Optional[str] = None,
    rel_threshold_expr: Optional[str] = None,
    abs_threshold_expr: Optional[str] = None,
) -> str:
    optional_arguments = ""

    if var_expr is not None:
        optional_arguments += f""", &
      scalar_expr={var_expr}"""

    if rel_threshold_expr is not None:
        optional_arguments += f""", &
      rel_threshold={rel_threshold_expr}"""

    if abs_threshold_expr is not None:
        optional_arguments += f""", &
      abs_threshold={abs_threshold_expr}"""

    if dtype == dace.int32:
        # could have been `logical` in fortran, might need cast
        ref_expr = f"transfer({ref_expr}, mold=int(1, kind=4))"

    return f"""
    call compare_{dtype.to_string()}_scalar( &
      actual={actual_expr}, &
      ref={ref_expr}, &
      result={result_expr}{optional_arguments} &
    )
"""


@generate_comparison_check_stmts.register
def _(scalar: dace.data.Scalar, *args, **kwargs) -> str:
    return generate_comparison_check_stmts(scalar.dtype, *args, **kwargs)


@generate_comparison_check_stmts.register
def _(
    array: dace.data.Array,
    actual_expr: str,
    ref_expr: str,
    result_expr: str,
    var_expr: Optional[str] = None,
    rel_threshold_expr: Optional[str] = None,
    abs_threshold_expr: Optional[str] = None,
) -> str:
    assert var_expr is not None

    dtype = array.dtype
    rank = len(array.shape)

    if array.dtype == dace.int32:
        # could have been `logical` in fortran, might need cast
        ref_expr = f"logical_fix_{rank}d({ref_expr})"

    optional_arguments = ""
    if rel_threshold_expr is not None:
        optional_arguments += f""", &
      rel_threshold={rel_threshold_expr}"""

    if abs_threshold_expr is not None:
        optional_arguments += f""", &
      abs_threshold={abs_threshold_expr}"""

    return f"""
    call compare_{dtype.to_string()}_{rank}d_array( &
        actual={actual_expr}, &
        ref={ref_expr}, &
        result={result_expr}, &
        array_expr={var_expr}{optional_arguments} &
    )
"""


@generate_comparison_check_stmts.register
def _(
    struct: dace.data.Structure,
    actual_expr: str,
    ref_expr: str,
    result_expr: str,
    var_expr: Optional[str] = None,
) -> str:
    assert var_expr is not None
    return f"""
    call compare_{struct.name}_struct( &
        actual={actual_expr}, &
        ref={ref_expr}, &
        result={result_expr}, &
        struct_expr={var_expr} &
    )
"""


@generate_comparison_check_stmts.register
def _(
    struct_array: dace.data.ContainerArray,
    actual_expr: str,
    ref_expr: str,
    result_expr: str,
    var_expr: Optional[str] = None,
):
    assert var_expr is not None

    struct = struct_array.stype
    rank = len(struct_array.shape)

    return f"""
    call compare_{struct.name}_{rank}d_array( &
        actual={actual_expr}, &
        ref={ref_expr}, &
        result={result_expr}, &
        struct_array_expr={var_expr} &
    )
"""


_INITIALIZATION_CHECKS_TYPES_IGNORE_LIST = {
    dace.int8,  # this likely was a string, so we ignore it
}

# This is exclusive, i.e., the marker state is not considered for initializations
_END_INITIALIZATIONS_MARKER_STATE_LABELS = {
    "GlobalDefEnd",
}

# all of these types should be immutable
_INITIALIZATION_CHECKS_SUPPORTED_TYPES = (
    bool,
    np.bool_,
    int,
    np.int32,
    np.int64,
    np.float32,
    np.float64,
    float,
)


def _initializations_check_bool_fix(val, dtype: dace.dtypes.typeclass) -> Any:
    if isinstance(val, (bool, np.bool_)) and not issubclass(
        dtype.type, (bool, np.bool_)
    ):
        # some Fortran `LOGICAL` are encoded as integers in SDFG
        return val
    return dtype(val)


def extract_initializations(
    sdfg: dace.SDFG, initializaion_checks_ignore_list: Set[str]
) -> Dict[str, Any]:
    initializations: Dict[str, Any] = {}

    for state in sdfg.states():
        if state.label in _END_INITIALIZATIONS_MARKER_STATE_LABELS:
            break
    else:
        # if no marker state was found, don't do initialization checks
        # (SDFG might be a cut-out without initializations)
        logging.warning(
            f"Could not extract initializations of globals ('{sdfg.name}')."
        )
        return initializations

    state = sdfg.start_state
    assert isinstance(
        state, dace.sdfg.SDFGState
    ), f"Start must be an SDFGState ({state})"
    initialization_sequence: List[
        Union[dace.sdfg.SDFGState, dace.sdfg.InterstateEdge]
    ] = []

    # get initialization_sequence
    while state.label not in _END_INITIALIZATIONS_MARKER_STATE_LABELS:
        initialization_sequence.append(state)
        edges = sdfg.out_edges(state)
        assert (
            len(edges) == 1
        ), f"initialization states must have only one out edge ({len(edges)})"
        edge = edges[0]
        initialization_sequence.append(edge.data)
        state = edge.dst

    # interpret initialization_sequence
    for state_or_edge in initialization_sequence:
        if isinstance(state_or_edge, dace.sdfg.SDFGState):

            nodes = state_or_edge.nodes()
            if len(nodes) == 0:
                continue

            # check valid tasklet
            tasklets = [
                tasklet for tasklet in nodes if isinstance(tasklet, dace.nodes.Tasklet)
            ]
            assert 1 == len(
                tasklets
            ), f"Expected only one tasklet (found {len(tasklets)})"
            tasklet = tasklets[0]
            assert isinstance(
                tasklet, dace.nodes.Tasklet
            ), f"Expected tasklet ({tasklet})"
            assert (
                tasklet.code.language == dace.dtypes.Language.Python
            ), f"Expected Python tasklet ({tasklet})"
            assert (
                len(tasklet.out_connectors) == 1
            ), f"Expected single out connector ({tasklet.out_connectors})"
            out_connector = next(iter(tasklet.out_connectors.keys()))

            # check valid write from tasklet
            out_memlets = state_or_edge.out_edges(tasklet)
            assert (
                len(out_memlets) == 1
            ), f"Expected only one out memlet (found {len(out_memlets)})"
            out_memlet = out_memlets[0]
            access_node = out_memlet.dst
            assert (
                out_memlet.src == tasklet and out_memlet.dst == access_node
            ), f"Expected memlet to go from Tasklet to AccessNode ({out_memlet.src}, {out_memlet.dst})"
            assert isinstance(
                access_node, dace.nodes.AccessNode
            ), f"Expected AccessNode ({access_node})"

            # better to use new dict in case the tasklet has other assignments
            # technically, we should only copy the symbols, but not data
            local_tasklet_scope = dict(initializations)
            for in_memlet in state_or_edge.in_edges(tasklet):
                src_node = in_memlet.src
                assert isinstance(
                    src_node, dace.nodes.AccessNode
                ), f"Expected AccessCode ({src_node})"
                local_tasklet_scope[in_memlet.dst_conn] = initializations[src_node.data]

            name = access_node.data
            assert name not in initializations, f"Double initializations of {name}"
            if name in initializaion_checks_ignore_list:
                logging.warning(
                    f"Skipping initialization check for '{name}' because it's on the ignore list."
                )
                continue

            desc = sdfg.arrays[name]
            assert isinstance(
                desc, dace.data.Scalar
            ), f"Expected scalar for initialization ({desc})"
            if desc.dtype in _INITIALIZATION_CHECKS_TYPES_IGNORE_LIST:
                logging.warning(
                    f"Skipping initialization check for '{name}' because its dtype "
                    f"'{desc.dtype}' is on the ignore list."
                )
                continue

            exec(tasklet.code.as_string, local_tasklet_scope)
            value = local_tasklet_scope[out_connector]
            assert isinstance(
                value, _INITIALIZATION_CHECKS_SUPPORTED_TYPES
            ), f"Unexpected type for initialization value ({value})"
            value = _initializations_check_bool_fix(value, desc.dtype)
            initializations[name] = value

        elif isinstance(state_or_edge, dace.sdfg.InterstateEdge):
            for name, expr in state_or_edge.assignments.items():
                assert name not in initializations, f"Double initializations of {name}"

                if name in initializaion_checks_ignore_list:
                    logging.warning(
                        f"Skipping initialization check for '{name}' because it's on the ignore list."
                    )
                    continue

                dtype = sdfg.symbols[name]
                if dtype in _INITIALIZATION_CHECKS_TYPES_IGNORE_LIST:
                    logging.warning(
                        f"Skipping initialization check for '{name}' because its dtype "
                        f"'{dtype}' is on the ignore list."
                    )
                    continue

                value = eval(expr, initializations.copy())
                assert isinstance(
                    value, _INITIALIZATION_CHECKS_SUPPORTED_TYPES
                ), f"Unexpected type for initialization value ({value})"
                value = _initializations_check_bool_fix(value, dtype)
                initializations[name] = dtype(value)
        else:
            assert False

    return initializations


@singledispatch
def generate_initialization_check_stmts(value, var_name: str, routine_name: str) -> str:
    raise NotImplementedError(
        f"Unable to generate check expression for variable '{var_name}' with value {value}"
    )


@generate_initialization_check_stmts.register
def _(value: bool, var_name: str, routine_name: str) -> str:
    return f"""
    if ({var_name} .neqv. {'.true.' if value else '.false.'}) then
      write (message_text,'(a,l,a)') &
        "variable incorrectly initialized: {var_name} = '", &
        {var_name}, &
        "' (in SDFG = '{value}')"
      print *, "{routine_name}"
      print *, message_text
    end if
"""


@generate_initialization_check_stmts.register(int)
@generate_initialization_check_stmts.register(np.int32)
@generate_initialization_check_stmts.register(np.int64)
def _(value: Union[int, np.int32, np.int64], var_name: str, routine_name: str) -> str:
    return f"""
    if ({var_name} /= {value}) then
      write (message_text,'(a,i0,a)') &
        "variable incorrectly initialized: {var_name} = '", &
        {var_name}, &
        "' (in SDFG = '{value}')"
      print *, "{routine_name}"
      print *, message_text
    end if
"""


@generate_initialization_check_stmts.register(float)
@generate_initialization_check_stmts.register(np.float32)
@generate_initialization_check_stmts.register(np.float64)
def _(
    value: Union[float, np.float32, np.float64], var_name: str, routine_name: str
) -> str:
    fortran_float_suffix = "_sp" if isinstance(value, np.float32) else "_wp"
    return f"""
    if ({var_name} /= {value:.20e}{fortran_float_suffix}) then
      write (message_text,'(a,e28.20,a)') &
        "variable incorrectly initialized: {var_name} = '", &
        {var_name}, &
        "' (in SDFG = '{value:.20e}')"
      print *, "{routine_name}"
      print *, message_text
    end if
"""


def generate_initializations_check(
    sdfg: dace.SDFG,
    imports_collector: ImportsCollector,
    initializaion_checks_ignore_list: Set[str],
) -> str:

    checks_str = ""
    unsimplified = None
    build_folder = Path(sdfg.build_folder)

    if (build_folder / "unsimplified.sdfg").is_file():
        unsimplified = sdfg.from_file(str(build_folder / "unsimplified.sdfg"))
    elif (build_folder / "unsimplified.sdfgz").is_file():
        unsimplified = sdfg.from_file(str(build_folder / "unsimplified.sdfgz"))

    if unsimplified is None:
        logging.warning(
            f"Could not find unsimplified sdfg for initialization ('{sdfg.name}')."
        )
    else:

        initializations = extract_initializations(
            unsimplified, initializaion_checks_ignore_list
        )

        routine_name = f"dace_fortran_interface_{sdfg.name}"

        for name, value in initializations.items():
            assert (
                name not in initializaion_checks_ignore_list
            ), f"initialization for '{name}' should be ignored"

            imports_collector.require_symbol(name)
            checks_str += generate_initialization_check_stmts(value, name, routine_name)

    return f"""\
  subroutine check_initializations()
    CHARACTER(len=5000) :: message_text = ''
{checks_str}
  end subroutine check_initializations
"""


def generate_fortran_interface_source(
    sdfg: dace.SDFG,
    struct_ignore_list: Set[str],
    struct_members_use_null: Dict[str, Set[str]],
    initializaion_checks_ignore_list: Set[str],
    module_definitions: Dict[str, str],
) -> str:
    source = ""

    sdfg_name = sdfg.name
    module_name = f"mo_{sdfg_name}_bindings"

    imports_collector = ImportsCollector(initializaion_checks_ignore_list)
    imports_collector.require_symbol("warning")
    #imports_collector.require_symbol("MAX_CHAR_LENGTH") # Requires module from ICON cant be used in ECRAD

    # they are sorted according to the C interfaces of the relevant functions
    sdfg_parameters = {
        fix_identifier(name): desc
        for name, desc in sdfg.arglist().items()
        if (
            not fix_identifier(name).startswith(
                _F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX
            )
            and not fix_identifier(name).startswith(
                _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX
            )
            and not fix_identifier(name).startswith("tmp_struct_symbol")
        )
    }
    optionals = {
        name.removeprefix(_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX)
        for name in sdfg_parameters
        if name.startswith(_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX)
    }

    # Using `dict` instead of `set` so insertion order is preserved
    structs: Dict[dace.data.Structure, None] = {}
    arrays: Dict[dace.data.Array, None] = {}
    for desc in sdfg_parameters.values():
        collect_structs_and_arrays(desc, structs, arrays, imports_collector)

    # FIXME(medium): check that alignment & data layout are compatible with Fortran for arrays & structs

    # the arrays that need translation functions
    array_translations: Dict[
        Tuple[Union[dace.dtypes.typeclass, dace.data.Structure], int],
        Union[dace.data.Array, dace.data.ContainerArray],
    ] = {}
    for array in arrays.keys():
        assert (
            array.alignment == 0
        ), f"Unsupported alignment for array ({array.alignment})"

        rank = len(array.shape)
        dtype = (
            array.stype if isinstance(array, dace.data.ContainerArray) else array.dtype
        )

        if (dtype, rank) not in array_translations:
            array_translations[(dtype, rank)] = array

    ###################################################
    # Fortran types to DaCe types copy in functions
    ###################################################
    copy_in_functions_str = ""
    copy_in_functions_interface_str = ""
    copy_back_subroutines_src = ""
    for struct in structs:
        copy_in_functions_str += generate_copy_in_function_struct(
            struct,
            struct_ignore_list,
            struct_members_use_null,
        )
        copy_back_subroutines_src += generate_copy_back_subroutine_struct(
            struct,
            struct_ignore_list,
            struct_members_use_null,
        )
    # we want to generate only one copy in procedure per base type & rank
    for array in array_translations.values():
        copy_in_str, copy_in_interface_str = generate_copy_in_function_array(array)
        copy_in_functions_str += copy_in_str
        copy_in_functions_interface_str += copy_in_interface_str
        if isinstance(array, dace.data.ContainerArray):
            copy_back_subroutines_src += generate_copy_back_subroutine_struct_array(
                array
            )

    initializations_check_subroutine_str = generate_initializations_check(
        sdfg,
        imports_collector,
        initializaion_checks_ignore_list,
    )

    ###################################################
    # Comparison functionality
    ###################################################

    comparison_subroutines_str = ""
    comparison_subroutines_str += _COMPARISON_PRIMITIVE_FUNCTIONS_STR
    for struct in structs:
        comparison_subroutines_str += generate_comparison_routine_dace_struct(
            struct,
            struct_ignore_list,
            struct_members_use_null,
        )
    for array in array_translations.values():
        comparison_subroutines_str += generate_array_comparison_subroutine(array)

    ###################################################
    # Header
    ###################################################

    imports_str = imports_collector.generate_imports_str(module_definitions)

    # TODO: better breadcrumbs (how the fortran interface file was generated)
    source += f"""\
! Auto-generated file by "{__file__}"
module {module_name}

  use iso_c_binding

{imports_str}

  implicit none

  private
  public :: run_{sdfg_name}
  public :: run_{sdfg_name}_verification
  public :: verify_{sdfg_name}
  public :: dace_init_{sdfg_name}
  public :: dace_exit_{sdfg_name}
  public :: dace_program_{sdfg_name}

"""

    ###################################################
    # DaCe struct definitions
    ###################################################
    dace_struct_definitions_str = "\n".join(
        generate_dace_struct_definition(struct) for struct in structs
    )

    verification_deep_copies_declarations_src = ""
    cached_shallow_copyies_declarations_src = ""
    for param_name, desc in sdfg_parameters.items():
        if isinstance(desc, dace.data.Scalar):
            # no copies needed for scalars
            continue

        assert isinstance(
            desc, (dace.data.Array, dace.data.Structure, dace.data.ContainerArray)
        )
        assert not param_name.startswith(_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX)
        assert not param_name.startswith(_F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX)
        assert not param_name.startswith(_F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX)
        assert not param_name.startswith(_F2DACE_STRUCT_ARRAY_SIZE_HELPER_FIELD_PREFIX)
        assert not param_name.startswith(
            _F2DACE_STRUCT_ARRAY_OFFSET_HELPER_FIELD_PREFIX
        )

        if isinstance(desc, (dace.data.Structure, dace.data.ContainerArray)):
            cached_shallow_copyies_declarations_src += f"""\
  type(c_ptr) :: cached_shallow_copy_{param_name} = C_NULL_PTR
"""

        verification_deep_copies_declarations_src += f"""\
  type(c_ptr) :: verification_deep_copy_{param_name} = C_NULL_PTR
"""

    source += f"""\

{dace_struct_definitions_str}

  logical :: is_initialized = .false.
  type(c_ptr) :: dace_state = C_NULL_PTR

{cached_shallow_copyies_declarations_src}
{verification_deep_copies_declarations_src}
"""

    ###################################################
    # Direct bindings
    ###################################################

    direct_program_parameters_without_state_str = join_wrapped(
        list(sdfg_parameters.keys()),
        seperator=", &\n    ",
        before=" &\n    ",
        after=" &\n  ",
    )
    direct_program_parameters_str = join_wrapped(
        ["state"] + list(sdfg_parameters.keys()),
        seperator=", &\n    ",
        before=" &\n    ",
        after=" &\n  ",
    )

    direct_program_parameter_decls_str = "\n".join(
        f"    {dace_type_to_fortran_c_var_type_decl(desc)}, value :: {name}"
        for name, desc in sdfg_parameters.items()
    )

    source += f"""\
interface

  type(c_ptr) function malloc(size) &
    bind(c, name="malloc")
    use iso_c_binding
    integer(kind=c_size_t), value :: size
  end function malloc

  subroutine free(ptr) &
    bind(c, name="free")
    use iso_c_binding
    type(c_ptr), value :: ptr
  end subroutine free


  type(c_ptr) function dace_init_{sdfg_name}({direct_program_parameters_without_state_str}) &
    bind(c, name="__dace_init_{sdfg_name}")
    use iso_c_binding

{direct_program_parameter_decls_str}
  end function dace_init_{sdfg_name}

  integer(c_int) function dace_exit_{sdfg_name}(state) &
    bind(c, name="__dace_exit_{sdfg_name}")
    use iso_c_binding

    type(c_ptr), value :: state
  end function dace_exit_{sdfg_name}

  subroutine dace_program_{sdfg_name}({direct_program_parameters_str}) &
    bind(c, name="__program_{sdfg_name}")
    use iso_c_binding

    type(c_ptr), value :: state
{direct_program_parameter_decls_str}
  end subroutine dace_program_{sdfg_name}

end interface
{copy_in_functions_interface_str}
{_COMPARISON_DEFAULT_THRESHOLDS_STR}

contains

{copy_in_functions_str}
{copy_back_subroutines_src}
{initializations_check_subroutine_str}
{comparison_subroutines_str}
"""

    ###################################################
    # Convenience bindings
    ###################################################
    convenience_parameters = {
        name: desc
        for name, desc in sdfg_parameters.items()
        if (
            not name.startswith(_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX)
            and not name.startswith(_F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX)
            and not name.startswith(_F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX)
        )
    }
    convinience_parameters_str = join_wrapped(
        list(convenience_parameters.keys()),
        seperator=", &\n    ",
        before=" &\n    ",
        after=" &\n  ",
    )

    convenience_parameter_decls_str = ""
    for param_name, desc in convenience_parameters.items():
        optionals_attribute = ", optional" if param_name in optionals else ""
        convenience_parameter_decls_str += f"""\
    {dace_type_to_fortran_rich_var_type_decl(desc, enable_inout_hack=True)}{optionals_attribute} :: {param_name}
"""

    convenience_locals_decls_str = ""
    # initialize optionals helper arguments
    initialize_optionals_src = ""
    for param_name in optionals:

        desc = sdfg_parameters[param_name]

        convenience_parameter_decls_str += f"""\
    !Optional helper parameter
    !WARNING: HACKFIX, POSSIBLE EXPLOSION
    integer(kind=c_int) :: {_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX}{param_name}
"""
        #{dace_type_to_fortran_rich_var_type_decl(desc, enable_inout_hack=True)} :: {_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX}{param_name}


        # FIXME(small): try to use _obvious wrong value_ if not present
        assign_op = "="
        default_initialization = ""
        if not isinstance(desc, dace.data.Scalar):
            assign_op = "=>"
            default_initialization = " => NULL()"

        convenience_locals_decls_str += f"""\
    {dace_type_to_fortran_rich_var_type_decl(desc, is_local=True, enable_inout_hack=True)} :: {_OPTIONAL_PROXY_PREFIX + param_name}{default_initialization}
"""
        initialize_optionals_src += f"""
    if (present({param_name})) then
      {_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX + param_name} = 1
      {_OPTIONAL_PROXY_PREFIX + param_name} {assign_op} {param_name}
    else
      {_F2DACE_PARAM_OPTIONAL_HELPER_PREFIX + param_name} = 0
    end if
"""

    kw_args: List[str] = []
    kw_args_verification: List[str] = []
    verification_shallow_copies_copy_ins_src = ""
    verification_deep_copies_copy_ins_src = ""
    shallow_copies_copy_back_src = ""
    verification_copies_comparisons_src = ""
    for param_name, desc in sdfg_parameters.items():
        var_name = param_name
        if param_name in optionals:
            var_name = _OPTIONAL_PROXY_PREFIX + param_name

        if param_name.startswith(
            _F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX
        ) or param_name.startswith(_F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX):
            continue

        kw_arg_verification = (
            kw_arg
        ) = f"""{param_name} = {generate_copy_in_fortran_expr(
                desc,
                expr=var_name,
                steal_arrays_expr=".true.",
                minimal_structs_expr=".false.",
                enable_inout_hack=True,
            )}"""

        if not isinstance(desc, dace.data.Scalar):
            kw_arg_verification = f"{param_name} = verification_deep_copy_{param_name}"
            verification_shallow_copies_copy_ins_src += f"""\
    verification_deep_copy_{param_name} = {generate_copy_in_fortran_expr(
        desc,
        expr=var_name,
        steal_arrays_expr=".true.",
        minimal_structs_expr=".false.",
        enable_inout_hack=True,
    )}
"""
            verification_deep_copies_copy_ins_src += f"""\
    verification_deep_copy_{param_name} = {generate_copy_in_fortran_expr(
        desc,
        expr=var_name,
        steal_arrays_expr=".false.",
        minimal_structs_expr=".false.",
        enable_inout_hack=True,
    )}
"""
            verification_copies_comparisons_src += f"""\
{generate_comparison_check_stmts(
    desc,
    actual_expr=f"verification_deep_copy_{param_name}",
    ref_expr=var_name,
    result_expr=f"local_result",
    var_expr=f'"{param_name}"',
)}
    result = result .and. local_result
"""
            if isinstance(desc, (dace.data.Structure, dace.data.ContainerArray)):
                shallow_copies_copy_back_src += generate_copy_back_stmts(
                    desc, var_name, f"verification_deep_copy_{param_name}"
                )
        kw_args.append(kw_arg)
        kw_args_verification.append(kw_arg_verification)

    array_helper_parameters = {
        name: desc
        for name, desc in sdfg_parameters.items()
        if (
            name.startswith(_F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX)
            or name.startswith(_F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX)
        )
    }
    for array_helper_name, desc in array_helper_parameters.items():
        assert isinstance(desc, dace.data.Scalar)

        if array_helper_name.startswith(_F2DACE_PARAM_ARRAY_SIZE_HELPER_FIELD_PREFIX):
            array_name, dim_num = extract_array_helper_information(
                array_helper_name, _F2DACE_PARAM_ARRAY_SIZE_HELPER_PATTERN_STR
            )
            helper_init_expr = f"size({array_name}, dim={dim_num + 1})"
        elif array_helper_name.startswith(
            _F2DACE_PARAM_ARRAY_OFFSET_HELPER_FIELD_PREFIX
        ):
            array_name, dim_num = extract_array_helper_information(
                array_helper_name, _F2DACE_PARAM_ARRAY_OFFSET_HELPER_PATTERN_STR
            )
            helper_init_expr = f"lbound({array_name}, dim={dim_num + 1})"
        else:
            assert False

        if array_name not in sdfg_parameters:
            logging.warning(
                f"Initializing array helper '{array_helper_name}' with '-1' (it shouldn't be a parameter)"
            )
            kw_arg_verification = kw_arg = f"{array_helper_name} = -1"
        else:
            kw_arg_verification = kw_arg = f"{array_helper_name} = {helper_init_expr}"
        kw_args.append(kw_arg)
        kw_args_verification.append(kw_arg_verification)

    kw_args_str = join_wrapped(
        ["state = dace_state"] + kw_args,
        seperator=", &\n      ",
        before=" &\n      ",
        after=" &\n    ",
    )
    kw_args_verification_src = join_wrapped(
        ["state = dace_state"] + kw_args_verification,
        seperator=", &\n      ",
        before=" &\n      ",
        after=" &\n    ",
    )
    kw_args_verification_without_state_src = join_wrapped(
        kw_args_verification,
        seperator=", &\n      ",
        before=" &\n      ",
        after=" &\n    ",
    )

    # FIXME(small): decided whether initializations should be checked in production run -> no?...
    source += f"""
    
  subroutine run_{sdfg_name}({convinience_parameters_str})
{convenience_parameter_decls_str}
{convenience_locals_decls_str}
{initialize_optionals_src}

    call check_initializations()

{verification_shallow_copies_copy_ins_src}

    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_{sdfg_name}({kw_args_verification_without_state_src})
    end if

    call dace_program_{sdfg_name}({kw_args_verification_src})

{shallow_copies_copy_back_src}

  end subroutine run_{sdfg_name}

  subroutine run_{sdfg_name}_verification({convinience_parameters_str})
{convenience_parameter_decls_str}
{convenience_locals_decls_str}
{initialize_optionals_src}

    call check_initializations()

{verification_deep_copies_copy_ins_src}

    if (is_initialized .eqv. .false.) then
      is_initialized = .true.
      dace_state = dace_init_{sdfg_name}({kw_args_verification_without_state_src})
    end if

    call dace_program_{sdfg_name}({kw_args_verification_src})
  end subroutine run_{sdfg_name}_verification

  subroutine verify_{sdfg_name}({convinience_parameters_str})
{convenience_parameter_decls_str}
{convenience_locals_decls_str}
    logical :: local_result, result = .true.
{initialize_optionals_src}

    if (is_initialized .eqv. .false.) then
      print *, "verify_{sdfg_name}: dace state is not initialized"
    end if

    call check_initializations()

{verification_copies_comparisons_src}
    if (.not. result) then
      print *, "verify_{sdfg_name}: Failed verification"
    else
      print *, "verify_{sdfg_name}: Verification successful :)"
    end if

  end subroutine verify_{sdfg_name}

end module {module_name}
"""

    return source


def main():

    parser = argparse.ArgumentParser(
        description=("Generate F90 Fortran binginds for SDFG file.")
    )

    parser.add_argument(
        "sdfg_file_path",
        help="Path to the SDFG file",
        type=Path,
    )
    parser.add_argument(
        "meta_data_yaml_path",
        help="Path to the meta_data yaml file",
        type=Path,
    )
    parser.add_argument(
        "module_definitions_file_path",
        help="Path to the module definitions yaml file",
        type=Path,
    )
    parser.add_argument(
        "bindings_file_path",
        help="Path to the bindings file which should be generated",
        type=Path,
    )

    args = parser.parse_args()

    sdfg = dace.SDFG.from_file(str(args.sdfg_file_path))

    with open(args.meta_data_yaml_path) as meta_data_yaml_file:
        meta_data = load_yaml(meta_data_yaml_file, Loader=YAML_Loader)
    struct_ignore_list = set(meta_data["struct_ignore_list"])
    struct_members_use_null = {
        struct_name: set(members_use_null)
        for struct_name, members_use_null in meta_data[
            "struct_members_use_null"
        ].items()
    }
    initializaion_checks_ignore_list = set(
        meta_data["initialization_checks_ignore_list"]
    )

    with open(args.module_definitions_file_path) as module_definitions_file:
        module_definitions = load_yaml(module_definitions_file, Loader=YAML_Loader)

    source = generate_fortran_interface_source(
        sdfg,
        struct_ignore_list,
        struct_members_use_null,
        initializaion_checks_ignore_list,
        module_definitions,
    )

    with open(args.bindings_file_path, "w") as fortran_interface_file:
        fortran_interface_file.write(source)


if __name__ == "__main__":
    main()
