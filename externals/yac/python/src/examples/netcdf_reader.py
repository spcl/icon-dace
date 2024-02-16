#!/usr/bin/env python3

##
# @file netcdf_reader.py
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
