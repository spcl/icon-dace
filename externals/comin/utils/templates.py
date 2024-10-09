#  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Please see the file LICENSE in the root of the source tree for this code.
#  Where software is supplied by third parties, it is indicated in the
#  headers of the routines.

def f_to_c_type(typ, kind):
    if typ == "INTEGER":
        if kind == "(c_signed_char)":
            return "int8_t"
        return "int"
    if typ == "CHARACTER":
        return "char"
    if typ == "REAL" and kind == "(wp)":
        return "double"
    if typ == "LOGICAL":
        return "bool"
    raise Exception("Unknown type")


def check_exist(name, attr):
    if "POINTER" in attr:
        return f"ASSOCIATED({name})"
    elif "ALLOCATABLE" in attr:
        return f"ALLOCATED({name})"
    else:
        return ".TRUE."


def f90_get_ptr(namespace, name, has_jg):
    print(f"""
  FUNCTION comin_descrdata_get_{namespace}_{name}({"jg" if has_jg else ""})
    USE iso_c_binding, ONLY: C_INT
    {"INTEGER(C_INT), INTENT(IN), VALUE  :: jg" if has_jg else ""}
    TYPE(t_comin_descrdata_{namespace}), POINTER :: p => NULL()
    TYPE(t_comin_descrdata_{namespace}_{name}), POINTER :: comin_descrdata_get_{namespace}_{name}

    p => comin_descrdata_get_{namespace}({"jg" if has_jg else ""})
    comin_descrdata_get_{namespace}_{name} => p%{name}
  END FUNCTION comin_descrdata_get_{namespace}_{name}""")


def f90_get_array(namespace, name, typ, ndim, attr, has_jg):
    print(f"""
  SUBROUTINE comin_descrdata_get_{namespace}_{name}({"jg," if has_jg else ""} {name}, arr_size) &
      &  BIND(C, NAME="comin_descrdata_get_{namespace}_{name}")
{"    INTEGER(C_INT), INTENT(IN), VALUE  :: jg" if has_jg else ""}
    TYPE(C_PTR),    INTENT(OUT) :: {name}
    INTEGER(C_INT), INTENT(INOUT) :: arr_size({ndim})
    !
    TYPE(t_comin_descrdata_{namespace}), POINTER :: p => NULL()
    p => comin_descrdata_get_{namespace}({"jg" if has_jg else ""})
    IF (.NOT. {check_exist(f"p%{name}", attr)}) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_{namespace}", &
      &  "ERROR: Pointer of {name} not associated.")
    END IF
    {f'arr_size(1) = LEN_TRIM(p%{name})' if typ == 'CHARACTER' and ndim == 1 else
     f'arr_size = SHAPE(p%{name})'}
    {name} = C_LOC(p%{name})
  END SUBROUTINE comin_descrdata_get_{namespace}_{name}""")


def c_get_array(namespace, name, typ, kind, ndim, has_jg):
    c_typ = f_to_c_type(typ, kind)
    print(f"""  void comin_descrdata_get_{namespace}_{name}({"int jg, " if has_jg else ""}{c_typ}** {name}, int* arr_size);""")


def f90_get_scalar(namespace, name, typ, kind, attr, has_jg):
    print(f"""
  FUNCTION comin_descrdata_get_{namespace}_{name}({"jg" if has_jg else ""}) &
      &  BIND(C, NAME="comin_descrdata_get_{namespace}_{name}") &
      &  RESULT({name})
{"    INTEGER(C_INT), INTENT(IN), VALUE  :: jg" if has_jg else ""}
    {typ}{kind if kind is not None else ""}                      :: {name}
    !
    TYPE(t_comin_descrdata_{namespace}), POINTER :: p => NULL()

    p => comin_descrdata_get_{namespace}({"jg" if has_jg else ""})
    IF (.NOT. {check_exist(f"p%{name}", attr)}) THEN
      CALL comin_plugin_finish("Message of comin_descrdata_query_{namespace}", &
      & "ERROR: Pointer of {name} not associated.")
    END IF
    {name}  = p%{name}
  END FUNCTION comin_descrdata_get_{namespace}_{name}""")


def c_get_scalar(namespace, name, typ, kind, has_jg):
    c_typ = f_to_c_type(typ, kind)
    print(f"""  {c_typ} comin_descrdata_get_{namespace}_{name}({"int jg" if has_jg else ""});""")


def gen_f90(namespace, struct, has_jg):
    for name, substruct in struct.items():
        if type(substruct) is dict:
            f90_get_ptr(namespace, name, has_jg)
            gen_f90(namespace+"_"+name, substruct, has_jg)
        else:
            typ, kind, attr, ndims = substruct
            if typ == "CHARACTER" and ndims == 0 and kind is not None:
                ndims = 1
            if ndims == 0:
                f90_get_scalar(namespace, name, typ, kind, attr, has_jg)
            else:
                f90_get_array(namespace, name, typ, ndims, attr, has_jg)


def gen_c(namespace, struct, has_jg):
    for name, substruct in struct.items():
        if type(substruct) is dict:
            gen_c(namespace+"_"+name, substruct, has_jg)
        else:
            typ, kind, attr, ndims = substruct
            if typ == "CHARACTER" and ndims == 0 and kind is not None:
                ndims = 1
            if ndims == 0:
                c_get_scalar(namespace, name, typ, kind, has_jg)
            else:
                c_get_array(namespace, name, typ, kind, ndims, has_jg)


def gen_c_properties(namespace, struct, has_jg):
    for name, substruct in struct.items():
        if type(substruct) is dict:
            gen_c_properties(namespace+"_"+name, substruct, has_jg)
    print(f"""
  const struct comin_descrdata_property_t comin_descrdata_{namespace}_properties[] = {{""")
    for name, substruct in struct.items():
        if type(substruct) is dict:
            print(f'    {{"{name}", 0, "void", 0, {"true" if has_jg else "false"}, comin_descrdata_{namespace}_{name}_properties }},')
        else:
            typ, kind, attr, ndims = substruct
            if typ == "CHARACTER" and ndims == 0 and kind is not None:
                ndims = 1
            c_typ = f_to_c_type(typ, kind)
            print(f'    {{"{name}", (void*)&comin_descrdata_get_{namespace}_{name}, "{c_typ}", {ndims}, {"true" if has_jg else "false"}, 0 }},')
    print("    {0,0,0,0,0}};")
