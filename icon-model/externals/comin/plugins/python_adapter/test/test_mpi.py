# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

from mpi4py import MPI
import comin
import sys

rank = comin.parallel_get_host_mpi_rank()

host_comm = MPI.Comm.f2py(comin.parallel_get_host_mpi_comm())
plugin_comm = MPI.Comm.f2py(comin.parallel_get_plugin_mpi_comm())
print(f"host_comm: {host_comm.rank}/{host_comm.size}", file=sys.stderr)
print(f"plugin_comm: {plugin_comm.rank}/{plugin_comm.size}", file=sys.stderr)

assert rank == host_comm.rank


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def sec_constructor():
    host_comm = MPI.Comm.f2py(comin.parallel_get_host_mpi_comm())
    plugin_comm = MPI.Comm.f2py(comin.parallel_get_plugin_mpi_comm())
    print(f"host_comm: {host_comm.rank}/{host_comm.size}",
          file=sys.stderr)
    print(f"plugin_comm: {plugin_comm.rank}/{plugin_comm.size}",
          file=sys.stderr)
