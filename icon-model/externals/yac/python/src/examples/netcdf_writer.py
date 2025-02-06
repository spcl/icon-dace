#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

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
        nnn.add_nnn(NNNReductionType.AVG, 1, 0., 1.)

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
