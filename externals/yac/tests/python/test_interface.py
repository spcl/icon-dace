#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

import yac
from mpi4py import MPI
import numpy as np
import os

rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size
assert size > 1

yac.def_calendar(yac.Calendar.PROLEPTIC_GREGORIAN)
assert yac.get_calendar() == yac.Calendar.PROLEPTIC_GREGORIAN
print(yac.version())

# class YAC
yac_def_instance = yac.YAC(default_instance=True)
yac_def_instance2 = yac.YAC.default_instance
yac_from_id = yac.YAC.from_id(yac_def_instance2.instance_id)
assert yac_def_instance.instance_id == yac_def_instance2.instance_id == yac_from_id.instance_id

yac_def_instance.set_config_output_file("foo.yaml",
                                        yac.ConfigOutputFormat.YAML,
                                        yac.ConfigOutputSyncLoc.ENDDEF)

del yac_def_instance, yac_def_instance2
yac_comm = yac.YAC(comm=MPI.COMM_WORLD)
del yac_comm
yac_instance = yac.YAC()

yac_instance.def_datetime("2020-01-01T00:00:00", "2020-01-02T00:00:00")

# components
if rank % 2 == 0:
    comp1 = yac_instance.predef_comp("comp1")
    comp2 = yac_instance.def_comp("comp2")
else:
    comp1, comp3 = yac_instance.def_comps(["comp1", "comp3"])

assert yac_instance.start_datetime == "2020-01-01T00:00:00.000"
assert yac_instance.end_datetime == "2020-01-02T00:00:00.000"

assert yac_instance.get_comps_comm(["comp2", "comp3"]).size == size
assert comp1.comp_comm.size == size
assert comp1.size == comp1.comp_comm.size
assert comp1.rank == comp1.comp_comm.rank


yac_instance.def_component_metadata("comp1", "COMP_METADATA".encode())

assert set(yac_instance.component_names) == {"comp1", "comp2", "comp3"}

# class Grid
x = np.array([-1., 0., 1.])
x.setflags(write=0)  # require a const memoryview in Reg2d
grid1 = yac.Reg2dGrid("grid1", x, x)
grid1.set_global_index(range(9), yac.Location.CORNER)
grid1.set_global_index(np.arange(4, dtype=np.intc), yac.Location.CELL)
print(f"{grid1.compute_grid_cell_areas()=}")
assert grid1.nbr_cells == 4
assert grid1.nbr_corners == 9
assert grid1.nbr_edges == 12
grid1.set_core_mask([True, False] * 2, yac.Location.CELL)
points1 = grid1.def_points(yac.Location.CELL, [-0.5, 0.5], [-0.5, 0.5])
assert points1.size == grid1.nbr_cells
points1.set_mask([True] * 4)

yac_instance.def_grid_metadata("grid1", "GRID_METADATA".encode())

grid2 = yac.UnstructuredGrid(
    "grid2",
    [3],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0, 1, 2],
    use_ll_edges=False,
)
mask1 = grid2.def_mask(yac.Location.CELL, [True], name="mask1")
points2 = grid2.def_points(yac.Location.CELL, [0.5], [0.5])

grid3 = yac.UnstructuredGrid(
    "grid3",
    [4],
    [0.0, 1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, 1.0],
    [0, 1, 2, 3],
    use_ll_edges=True,
)

grid4 = yac.CloudGrid("grid4",
                      [0.0, 1.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0, 1.0])

points4 = grid4.def_points([0.0, 1.0, 1.0, 0.0],
                           [0.0, 0.0, 1.0, 1.0])

# class Field
field1 = yac.Field.create("field1", comp1, points1, 1, "1", yac.TimeUnit.MINUTE)
yac_instance.def_field_metadata(
    "comp1", "grid1", "field1", "FIELD_METADATA".encode()
)
yac_instance.enable_field_frac_mask("comp1", "grid1", "field1", 42.0)

field2 = yac.Field.create(
    "field2", comp1, points2, 1, "1", yac.TimeUnit.HOUR, mask1
)
field3 = yac.Field.create(
    "field3", comp1, points1, 1, "1", yac.TimeUnit.HOUR
)
field4 = yac.Field.create(
    "field4", comp1, points2, 1, "1", yac.TimeUnit.HOUR, mask1
)

# class InterpolationStack
interp = yac.InterpolationStack()
interp.add_nnn(yac.NNNReductionType.AVG, 1, 0.0, 1.0)
yac.InterpolationStack().add_average(yac.AverageReductionType.AVG_DIST, 1)
yac.InterpolationStack().add_conservative(
    1, 1, 1, yac.ConservNormalizationType.FRACAREA
)
yac.InterpolationStack().add_spmap(
    1.0, 1.0, yac.SPMAPWeightType.AVG, yac.SPMAPScaleType.NONE, 1.0, 1.0,
)
yac.InterpolationStack().add_hcsbb()
#  interp.add_user_file()
yac.InterpolationStack().add_fixed(999.0)
yac.InterpolationStack().add_check("ctor", "search key")
yac.InterpolationStack().add_creep(42)

yac_instance.def_couple(
    "comp1",
    "grid1",
    "field1",
    "comp1",
    "grid2",
    "field2",
    "60",
    yac.TimeUnit.MINUTE,
    yac.Reduction.TIME_NONE,
    interp,
    src_lag=2,
)
del interp

interp_yaml = yac.InterpolationStack.from_string_yaml(
    "- nnn:\n    n: 1"
)
del interp_yaml

interp_json = yac.InterpolationStack.from_string_json(
    "[{\"nnn\": {\"n\": 1}}]"
)
del interp_json

if rank == 0:
    with open("config.yaml", "w") as f:
        f.write(
            """timestep_unit: minute
coupling:
  - src_component: comp1
    src_grid: grid1
    src_lag: 0
    tgt_component: comp1
    tgt_grid: grid2
    coupling_period: 60
    time_reduction: none
    interpolation:
      - nnn:
          n: 1
    field:
      src: field3
      tgt: field4
"""
        )
    yac_instance.read_config_yaml("config.yaml")
    os.remove("config.yaml")

yac_instance.sync_def()

assert yac_instance.get_component_metadata("comp1") == "COMP_METADATA"
assert yac_instance.get_grid_metadata("grid1") == "GRID_METADATA"
assert (
    yac_instance.get_field_metadata("comp1", "grid1", "field1")
    == "FIELD_METADATA"
)

assert "grid1" in yac_instance.grid_names
assert "grid1" in yac_instance.get_comp_grid_names("comp1")
assert "field1" in yac_instance.get_field_names("comp1", "grid1")
field1_id = yac_instance.get_field_id("comp1", "grid1", "field1")
assert (field1.component_name == "comp1")
assert (field1.grid_name == "grid1")
assert (
    yac_instance.get_field_timestep("comp1", "grid1", "field1")
    == field1.timestep
    == "PT01M"
)
assert (
    yac_instance.get_field_role("comp1", "grid1", "field1")
    == yac.ExchangeType.SOURCE
)
assert (yac_instance.get_field_collection_size("comp1", "grid1", "field1")
        == field1.collection_size
        == 1)
assert (
    yac_instance.get_field_frac_mask_fallback_value("comp1", "grid1", "field1")
    == 42.0
)

yac_instance.enddef()

assert (field1.role == yac.ExchangeType.SOURCE ) # it is set only after enddef

data, send_info, recv_info = yac.Field.exchange(field3, field4, np.random.rand(2, 2))
field1.update()
