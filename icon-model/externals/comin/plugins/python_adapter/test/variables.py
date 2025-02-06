# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

import comin

print(comin.current_get_plugin_info())

comin.var_request_add(("test", 1), False)
comin.metadata_set(("test", 1), tracer=False, tracer_turb=False, units="%",
                   hgrid_id=comin.COMIN_HGRID_UNSTRUCTURED_CELL)
comin.metadata_set(("test", 1), custom_metadata=42)


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def secondary_constructor():
    global test
    test = comin.var_get([comin.EP_ATM_WRITE_OUTPUT_BEFORE], ("test", 1),
                         comin.COMIN_FLAG_READ)

    # print all variables
    for var_name, var_id in comin.var_descr_list():
        print(f"{var_name=}, {var_id=}, tracer={comin.metadata_get((var_name, var_id), 'tracer')} units={comin.metadata_get((var_name, var_id), 'units')}")

    # use the metadata as a dictionary:
    print(comin.metadata(("test", 1))["custom_metadata"])

    for key, data in comin.metadata(("test", 1)).items():
        print(f"{key}: {data}")



@comin.register_callback(comin.EP_ATM_WRITE_OUTPUT_BEFORE)
def output_before():
    print(f"{test.descriptor=}")
    print(f"{test.to_3d.shape=}")
