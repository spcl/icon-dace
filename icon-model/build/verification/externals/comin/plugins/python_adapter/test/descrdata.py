# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

import comin
import re
import sys


# As the CI checks the output for changes we remove the memory addresses
def dump(s):
    print(re.sub("at 0x[0-9a-f]*>", "at BEEEEEEP>", str(s)))


print(comin.current_get_plugin_info())

domain = comin.descrdata_get_domain(1)
print(f"{dir(domain)=}")
dump(f"{dir(domain.cells)=}")
dump(f"{domain.cells.clon=}")
dump(f"{ {name: getattr(domain, name) for name in dir(domain)}=}")

glob = comin.descrdata_get_global()
dump(f"{ {name: getattr(glob, name) for name in dir(glob)}=}")

print(f"{comin.descrdata_get_simulation_interval()=}")
print(f"{comin.descrdata_get_cell_indices(42, 1, 3, 2, 3, 4)=}")
print(f"{comin.descrdata_get_cell_npromz(1)=}")
print(f"{comin.descrdata_get_edge_npromz(1)=}")
print(f"{comin.descrdata_get_vert_npromz(1)=}")
print(f"{comin.descrdata_index_lookup_glb2loc_cell(1,42)=}")
print(f"{comin.setup_get_version()=}")
print(f"{comin.setup_get_verbosity_level()=}")

print(f"{comin.parallel_get_host_mpi_rank()=}")


@comin.register_callback(comin.EP_ATM_PHYSICS_BEFORE)
def phy():
    print(f"{comin.current_get_domain_id()=}")
    print(f"{comin.current_get_datetime()=}")
