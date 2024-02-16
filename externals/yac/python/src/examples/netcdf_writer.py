#!/usr/bin/env python3

##
# @file netcdf_writer.py
# @brief gets data from yac and write it to a netCDF file
#
# @copyright Copyright  (C)  2023 DKRZ, MPI-M
#
# @author Nils-Arne Dreier <dreier@dkrz.de>
#
#
# Keywords:
# Maintainer: Nils-Arne Dreier <dreier@dkrz.de>
# URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
#
# This file is part of YAC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are  permitted provided that the following conditions are
# met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# Neither the name of the DKRZ GmbH nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

from yac import *
from netCDF4 import Dataset

import datetime
import cftime
import numpy as np
from mpi4py import MPI

from yac.utils import read_grid


## todo interpolation (currently only NN supported)
## todo support collection size

class NetCDF_Writer:
    def __init__(self, filename, timestep, gridfile, variables,
                 compname = "netcdf_writer", gridname = "netcdf_writer_grid",
                 compression_kwargs = {}, yac = None):
        self.yac = yac or YAC.default_instance
        self.comp = self.yac.predef_comp(compname)
        self.compname = compname
        self.filename = filename
        self.timestep = timestep
        self.gridfile = gridfile
        self.variables = variables
        self.gridname = gridname
        self.compression_kwargs = compression_kwargs

    def setup(self):
        comm = self.comp.comp_comm
        if comm.size > 1:
            self.dataset = Dataset(self.filename, "w", parallel=True, comm=comm)
        else:
            self.dataset = Dataset(self.filename, "w")

        start = datetime.datetime.fromisoformat(self.yac.start_datetime)
        end = datetime.datetime.fromisoformat(self.yac.end_datetime)
        no_timesteps = int((end-start)/self.timestep)+1
        time_dim = self.dataset.createDimension("time", no_timesteps)
        time_var = self.dataset.createVariable("time", "f8", ("time",))
        time_var.units = "seconds since "+self.yac.start_datetime
        time_range = [start+i*self.timestep for i in range(no_timesteps)]
        time_var[:] = cftime.date2num(time_range, units=time_var.units)

        grid, self.idx, _ = read_grid(self.gridfile, self.gridname, self.comp.comp_comm)
        self.points = grid.cell_points

        self.dataset.createDimension("cell", self.points.size)

    def def_couples(self):
        nnn = InterpolationStack()
        nnn.add_nnn(NNNReductionType.AVG, 1, 1.)

        self.fields = []
        for var in self.variables:
            collection_size = self.yac.get_field_collection_size(*var)
            dt_ms = str(int(self.timestep.total_seconds())*1000)
            self.fields.append(Field.create(var[2], self.comp, self.points, 1,
                                            dt_ms, TimeUnit.MILLISECOND))
            self.yac.def_couple(*var,
                                self.compname, self.gridname, var[2],
                                dt_ms, TimeUnit.MILLISECOND,
                                0, nnn)

            self.dataset.createVariable(var[2], "f4", ("time", "cell" ),
                                        **self.compression_kwargs)
            self.time_counter = 0

    def step(self):
        for field in self.fields:
            print(f"writing {field.name} at {field.datetime}")
            buf, info = field.get()
            self.dataset[field.name][self.time_counter,self.idx] = buf
        self.time_counter += 1
        return field.datetime
