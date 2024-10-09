#!/usr/bin/env python3
#  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Please see the file LICENSE in the root of the source tree for this code.
#  Where software is supplied by third parties, it is indicated in the
#  headers of the routines.

from fparser.two.parser import ParserFactory
from fparser.common.readfortran import FortranFileReader
from fparser.two.utils import walk, get_child
from fparser.two.Fortran2003 import *
from collections import namedtuple

component = namedtuple("component", ["type", "kind", "attributes", "ndims"],
                       defaults=[None, None, [], 0])


def find_type(filename, typename):
    reader = FortranFileReader(filename)
    f2008_parser = ParserFactory().create(std="f2008")
    parse_tree = f2008_parser(reader)
    #print(parse_tree.__repr__())
    types = walk(parse_tree, Derived_Type_Def)
    for t in types:
        name_smt = walk(t, Type_Name)[0]
        if name_smt is not None and name_smt.string == typename:
            components = {}
            for c in walk(t, Data_Component_Def_Stmt):
                t = get_child(c, Intrinsic_Type_Spec)
                name = get_child(walk(c, Component_Decl)[0], Name).string
                if t is None: # this is a derived type
                    dt_name = walk(c, Type_Name)[0].string
                    components[name] = find_type(filename, dt_name)
                else: # this is a intrinsic type
                    ndims = len(walk(c, Deferred_Shape_Spec))
                    attr = [e.string for e in walk(c, Component_Attr_Spec)]
                    kind = t.items[1].string if t.items[1] is not None else None
                    components[name] = component(type=t.items[0],
                                                 kind=kind,
                                                 attributes=attr,
                                                 ndims=ndims)
            return components


if __name__ == "__main__":
    import pprint
    t = find_type("../src/comin_descrdata_types.F90", "t_comin_descrdata_domain")
    pprint.pprint(t)
