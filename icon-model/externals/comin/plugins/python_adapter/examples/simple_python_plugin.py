"""
Test plugin for the ICON Community Interface (ComIn)

This simple test plugin shows how to use the basic features of
ComIn analogous to simple_c_plugin and simple_fortran_plugin.

@authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

SPDX-License-Identifier: BSD-3-Clause

Please see the file LICENSE in the root of the source tree for this code.
Where software is supplied by third parties, it is indicated in the
headers of the routines.
"""

import numpy as np
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--hello", help="just a dummy argument")
parser.add_argument("arg", help="just a dummy argument")

# allows the script to be called from the command line e.g. with "--help"
try:
    import comin
    args = parser.parse_args(comin.current_get_plugin_info().args)
except ImportError:
    parser.parse_args()

print(args)
print(comin.current_get_plugin_info())

# request to add a variable
comin.var_request_add(("simple_python_var", 1), False)

# request to add a tracer
comin.var_request_add(("simple_python_tracer", -1), False)
comin.metadata_set(("simple_python_tracer", -1), tracer=True)


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def simple_python_constructor():
    global pres, simple_python_var, simple_python_tracer
    print("simple_python_constructor called!", file=sys.stderr)
    pres = comin.var_get([comin.EP_ATM_WRITE_OUTPUT_BEFORE], ("pres", 1),
                         comin.COMIN_FLAG_READ)
    simple_python_var = comin.var_get([comin.EP_ATM_WRITE_OUTPUT_BEFORE],
                                      ("simple_python_var", 1),
                                      comin.COMIN_FLAG_READ | comin.COMIN_FLAG_WRITE)
    simple_python_tracer = comin.var_get([comin.EP_ATM_WRITE_OUTPUT_BEFORE],
                                         ("simple_python_tracer", 1),
                                      comin.COMIN_FLAG_READ | comin.COMIN_FLAG_WRITE)
    print("tracer pos:", simple_python_tracer.pos)
    print("tracer ncontained:", simple_python_tracer.ncontained)


@comin.register_callback(comin.EP_ATM_WRITE_OUTPUT_BEFORE)
def simple_python_diagfct():
    print("simple_python_diagfct called!", file=sys.stderr)
    np.asarray(simple_python_var)[:] = np.asarray(pres)+42.
    np.asarray(simple_python_tracer)[:] = np.asarray(simple_python_var) / 1337.


@comin.register_callback(comin.EP_DESTRUCTOR)
def simple_python_destructor():
    print("simple_python_destructor called!", file=sys.stderr)
