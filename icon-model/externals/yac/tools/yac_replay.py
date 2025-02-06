#!/usr/bin/env python3
# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

from math import prod
import numpy as np
import xarray as xr
import grid_utils
import yac
import argparse
import logging


def parse_duration(s):
    import isodate
    print(s)
    dt = isodate.parse_duration(s)
    return np.timedelta64(int(dt.total_seconds()*10e6), "us")


parser = argparse.ArgumentParser("yac_replay",
                                 description="""Replay simulation from a dataset.
The given files are loaded with `xarray.open_mfdataset`.
Currently only cell based fields are supported that have 'cell' as dimension.
It is currently not supported to run this component in parallel.""")
parser.add_argument("files", type=str, nargs='+',
                    default="xarray URL for the data file")
parser.add_argument("--engine", type=str, required=False, default=None,
                    help="xarray backend engine (see xarray.open_mfdataset)")
parser.add_argument("--fallback_grid", type=str, required=False,
                    help="Description of grid (default: try to deduce from file")
parser.add_argument("--compname", type=str, default="replay", required=False,
                    help="Name for the yac component (default: replay)")
parser.add_argument("--gridprefix", type=str, default="replay_", required=False,
                    help="Prefix for the grid name (default: replay_)")
parser.add_argument("--start_date", type=np.datetime64, required=False,
                    help="start time for date")
parser.add_argument("--end_date", type=np.datetime64, required=False,
                    help="end time for date")
parser.add_argument("--dt", type=parse_duration, required=False,
                    help="timestep")
parser.add_argument('--debug', help="Print lots of debugging statements",
                    action="store_const", dest="loglevel", const=logging.DEBUG,
                    default=logging.WARNING)
parser.add_argument('-v', '--verbose', help="Be verbose",
                    action="store_const", dest="loglevel", const=logging.INFO)

args = parser.parse_args()
log_handler = logging.StreamHandler()
log_handler.setFormatter(logging.Formatter("%(asctime)-10s %(name)-10s %(levelname)-8s %(message)s"))
logging.basicConfig(level=args.loglevel, handlers=[log_handler])

logging.debug("set calendar to PROLEPTIC_GREGORIAN")
yac.def_calendar(yac.Calendar.PROLEPTIC_GREGORIAN)
logging.debug("Instantiate yac")
y = yac.YAC()

logging.debug(f"open dataset: {args.files}")
ds = xr.open_mfdataset(args.files, engine=args.engine)

start_date = args.start_date or np.datetime64(ds.time[0].data)
end_date = args.end_date or np.datetime64(ds.time[-1].data)

if args.dt:
    dt = args.dt
else:
    dt = np.diff(ds.coords["time"])
    assert np.all(dt[0] == dt), "Time coordinates are not equidistant"
    dt = dt[0]

logging.info(f"define datetime: {start_date} - {end_date}")
y.def_datetime(np.datetime_as_string(start_date)[:23],
               np.datetime_as_string(end_date)[:23])

logging.debug(f"define component: {args.compname}")
comp = y.def_comp(args.compname)

logging.info(f"Time interval: {y.start_datetime} - {y.end_datetime}")
logging.info(f"dt: {dt} timesteps: {(end_date - start_date)/dt}")


grids = {}
def grid_cache(grid_spec):
    global grids
    if grid_spec not in grids:
        logging.info(f"Defining grid {args.gridprefix+grid_spec}")
        grids[grid_spec] = grid_utils.get_grid(grid_spec).def_yac(args.gridprefix+grid_spec, [yac.Location.CELL, yac.Location.CORNER])
    return grids[grid_spec]


def get_grid(variable):
    global ds
    if {"lat", "lon"} <= set(variable.dims):
        lat = ds.lat
        lon = ds.lon
        return grid_cache(f"g{len(lon)},{len(lat)},{float(lon.min())},{float(lon.max())},{float(lat.min())},{float(lat.max())}")
    if "grid_mapping" in variable.attrs:
        crs = variable.attrs["grid_mapping"]
        if ds[crs].attrs["grid_mapping_name"] == "healpix":
            nside = ds[crs].attrs['healpix_nside']
            nest = "n" if ds[crs].attrs['healpix_order'] == "nest" else "r"
            return grid_cache(f"h{nside}{nest}")
    if args.fallback_grid:
        return grid_cache(args.fallback_grid)
    raise Exception(f"Cannot not determine grid for variable {variable.name}. Please specify a fallback grid.")


fields = []

for varname in ds.variables:
    var = ds[varname]
    logging.debug(f"checking variable {varname}")
    if "time" not in var.dims:
        continue
    if "cell" in var.dims or {"lat", "lon"} <= set(var.dims):
        cells, points = get_grid(var)
        point_id = cells
    else:
        continue
    collection_size = prod(var.sizes[d] for d in var.dims if d != "time" and d != "cell")
    logging.debug(f"define field with collection size: {collection_size}")
    fields.append((yac.Field.create(varname, comp, point_id, collection_size,
                                    str(int(dt / np.timedelta64(1, 's'))), yac.TimeUnit.SECOND),
                   var))

logging.info("calling enddef")
y.enddef()

fields = [(field, var) for field, var in fields if field.role is yac.ExchangeType.SOURCE]
logging.info(f"have {len(fields)} with role source")

for t in np.arange(start_date, end_date+dt, dt):
    logging.info(f"Time: {t}")
    for field, var in fields:
        logging.info(f"calling put for field {var.name}")
        field.put(var.sel(time=t).data)

logging.info("done")
