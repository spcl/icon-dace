#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

import yac
import numpy as np
import asyncio

yac.def_calendar(yac.Calendar.PROLEPTIC_GREGORIAN)

yac_instance = yac.YAC()

yac_instance.def_datetime("2020-01-01T00:00:00", "2020-01-02T00:00:00")

comp1, comp2 = yac_instance.def_comps(["comp1", "comp2"])

grid1 = yac.Reg2dGrid("grid1", [-1, 0, 1], [-1, 0, 1])
points1 = grid1.def_points(yac.Location.CELL, [-0.5, 0.5], [-0.5, 0.5])

grid2 = yac.Reg2dGrid("grid2", [-1, 0, 1], [-1, 0, 1])
points2 = grid2.def_points(yac.Location.CELL, [-0.5, 0.5], [-0.5, 0.5])

field1 = yac.Field.create("field1", comp1, points1, 1, "1", yac.TimeUnit.HOUR)
field2 = yac.Field.create("field2", comp2, points2, 1, "1", yac.TimeUnit.HOUR)
field3 = yac.Field.create("field3", comp2, points2, 1, "2", yac.TimeUnit.HOUR)

interp = yac.InterpolationStack()
interp.add_nnn(yac.NNNReductionType.AVG, 1, 0.0, 1.0)

yac_instance.def_couple(
    "comp1",
    "grid1",
    "field1",
    "comp2",
    "grid2",
    "field2",
    "60",
    yac.TimeUnit.MINUTE,
    yac.Reduction.TIME_NONE,
    interp,
)

yac_instance.def_couple(
    "comp1",
    "grid1",
    "field1",
    "comp2",
    "grid2",
    "field3",
    "120",
    yac.TimeUnit.MINUTE,
    yac.Reduction.TIME_NONE,
    interp,
)

yac_instance.enddef()


async def recv_loop(field, delay=False):
    while True:
        buf, info = await field.get_coro()
        print(f"recv (field_id={field.field_id}):", buf, info)
        if delay:
            await asyncio.sleep(0.01)
        if info == yac.Action.GET_FOR_RESTART:
            break


async def send_loop(field):
    while True:
        info = await field.put_coro(np.random.rand(4))
        print(f"send (field_id={field.field_id}):", info)
        if info == yac.Action.PUT_FOR_RESTART:
            break


# asyncio entry point
async def main():
    # create task group
    await asyncio.gather(
        send_loop(field1), recv_loop(field2), recv_loop(field3, True)
    )


# entry point
asyncio.run(main())
print("done!")
