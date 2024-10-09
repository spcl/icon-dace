#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

from yac import *
from netCDF4 import Dataset

import cftime
from yac.utils import read_grid

class NetCDF_Reader:
    def __init__(self, filename, gridfile = None,
                 compname = "netcdf_reader", yac = None):
        self.yac = yac or YAC.default_instance
        self.comp = self.yac.predef_comp(compname)

        self.dataset = Dataset(filename, "r")
        time_var = self.dataset["time"]
        start = cftime.num2date(time_var[0], units=time_var.units, calendar="proleptic_gregorian")
        end = cftime.num2date(time_var[-1], units=time_var.units, calendar="proleptic_gregorian")
        self.timestep = cftime.num2date(time_var[1], units=time_var.units, calendar="proleptic_gregorian") - start

        def_calendar(Calendar.PROLEPTIC_GREGORIAN)
        self.yac.def_datetime(start.isoformat(), end.isoformat())

        # TODO: read it from dataset metadata
        self.gridfile = gridfile or filename
        self.compname = compname
        self.time_counter = 0

    def setup(self):
        grid, _, _ = read_grid(self.gridfile, f"{self.compname}_grid")

        self.fields = []
        for name, v in self.dataset.variables.items():
            if len(v.dimensions) < 2:
                continue
            if "time" != v.dimensions[0]:
                continue
            if "cell" == v.dimensions[1]:
                assert grid.cell_points, "cells (clon, clat) not defined in the grid"
                points = grid.cell_points
            elif "vertex" == v.dimensions[1]:
                points = grid.corner_points
            else:
                continue
            # todo check for level o.Ä.
            self.fields.append(Field.create(name, self.comp, points, 1,
                                            str(self.timestep.seconds), TimeUnit.SECOND))

    def def_couples(self):
        pass

    def step(self):
        for field in self.fields:
            print(f"reading {field.name} at {field.datetime}")
            field.put(self.dataset[field.name][self.time_counter, :])
        self.time_counter += 1
        return self.fields[0].datetime
