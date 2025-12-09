#!/usr/bin/env python3

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

from yac import YAC, Reg2dGrid, Location, Field, TimeUnit, InterpolationStack, NNNReductionType, Action
from netCDF4 import Dataset
import numpy as np
import argparse
import isodate
import cftime
import sys
import yaml

parser = argparse.ArgumentParser(description="Simple output component. Files can be inspected with ncview.")
parser.add_argument("filename", type=str, help="output filename")
parser.add_argument("src", type=lambda s: tuple(s.split(",")),
                    nargs="*", help="(component, grid, field) tuple of the source field (komma-separated)")
parser.add_argument("--bounds", type=float, nargs=4, help="bounds of the regular grid (min_lon, max_lon, min_lat, max_lat)",
                    default=(-np.pi, np.pi, -0.5*np.pi, 0.5*np.pi))
parser.add_argument("--res", type=int, nargs=2, help="resolution for lon-lat grid",
                    default=(360, 180))
parser.add_argument("--output-interval", type=str, required=False,
                    help="The interval for the time dimension in the output file")
args = parser.parse_args()

yac = YAC()

comp_name = f"simple_output_{hash(tuple(args.src))}"
comp = yac.def_comp(comp_name)

vlon = np.linspace(np.deg2rad(args.bounds[0]), np.deg2rad(args.bounds[1]), args.res[0], endpoint=True)
vlat = np.linspace(np.deg2rad(args.bounds[2]), np.deg2rad(args.bounds[3]), args.res[1], endpoint=True)
grid_name = f"{comp_name}_grid"
grid = Reg2dGrid(grid_name, vlon, vlat)
clon = 0.5*(vlon[1:]+vlon[:-1])
clat = 0.5*(vlat[1:]+vlat[:-1])
points = grid.def_points(Location.CELL, clon, clat)

yac.sync_def()

dt = args.output_interval or yac.get_field_timestep(*args.src[0]) # use timestep of first field by default
print(f"Timestep: {dt}", file=sys.stderr)

nnn = InterpolationStack()
nnn.add_nnn(NNNReductionType.AVG, 1, 0.0, 1.)

levels = None

fields = {}
for field_desc in args.src:
    collection_size = yac.get_field_collection_size(*field_desc)
    if collection_size > 1:
        assert levels is None or levels == collection_size
        levels = collection_size
    fields[field_desc] = Field.create(field_desc[2], comp, points, collection_size, dt, TimeUnit.ISO_FORMAT)
    yac.def_couple(*field_desc,
                   comp_name, grid_name, field_desc[2],
                   dt, TimeUnit.ISO_FORMAT, 0, nnn)

yac.enddef()

dataset = Dataset(args.filename, "w")

start = isodate.parse_datetime(yac.start_datetime)
end = isodate.parse_datetime(yac.end_datetime)
timestep = isodate.parse_duration(dt)
no_timesteps = int((end-start)/timestep)
print(f"{no_timesteps=}", file=sys.stderr)
time_dim = dataset.createDimension("time", no_timesteps)
time_var = dataset.createVariable("time", "i8", ("time",))
time_var.units = "seconds since "+yac.start_datetime
time_range = [start+i*timestep for i in range(no_timesteps)]
time_var[:] = cftime.date2num(time_range, units=time_var.units)

dataset.createDimension("clat", len(clat))
dataset.createDimension("clon", len(clon))
if levels is not None:
    dataset.createDimension("level", levels)
    level = ("level", )

var = {}
for field_desc in args.src:
    if fields[field_desc].collection_size > 1:
        var[field_desc] = dataset.createVariable(field_desc[2], "f4", ("time", "level", "clat", "clon"))
    else:
        var[field_desc] = dataset.createVariable(field_desc[2], "f4", ("time", "clat", "clon"))
    metadata = yaml.safe_load(yac.get_field_metadata(*field_desc)).items()
    for k, v in metadata:
        setattr(var[field_desc], k, v)

data = None
for t in range(no_timesteps):
    for field, v in zip(fields.values(), var.values()):
        print(f"Writing {field.name} at {field.datetime} ({t})", file=sys.stderr)
        data, info = field.get(data)
        v[t, :] = data.reshape((-1, len(clat), len(clon)))
print("Finish!", file=sys.stderr)
