#!/usr/bin/env python3

from yac import YAC, Curve2dGrid, Location, Field, TimeUnit
from itertools import product
import math
import numpy as np
from mpi4py import MPI


def lonlat2xyz(lon, lat):
    clat = np.cos(lat)
    slat = np.sin(lat)
    return clat * np.cos(lon), clat * np.sin(lon), slat


def xyz2lonlat(x, y, z):
    lat = np.arcsin(z)
    lon = np.arctan2(y, x)
    return lon, lat


def rotmat(theta, plane_dims, d=2):
    e = np.eye(d)
    e[np.ix_(plane_dims, plane_dims)] = np.array([[np.cos(theta), -np.sin(theta)],
                                                  [np.sin(theta), np.cos(theta)]])
    return e


MPI.COMM_WORLD.Barrier()

yac = YAC()
my_comp_name = "PISM"
comp = yac.def_comp(my_comp_name)

size = comp.comp_comm.size
rank = comp.comp_comm.rank


def msg(m):
    if rank == 0:
        print("toy_multi_curve2d: ", m)


# factorize size
fac = int(np.sqrt(size))
while size % fac != 0:
    fac = fac-1
msg(f"decomposition: {fac} x {size//fac}")

dim2d = [100, 50]
size2d = [fac, size//fac]
rank2d = [rank % fac, rank//fac]
chunk2d = [math.ceil(dim2d[0]/size2d[0]), math.ceil(dim2d[1]/size2d[1])]
slice2d = [slice(rank2d[0]*chunk2d[0], (rank2d[0]+1)*chunk2d[0]+1),
           slice(rank2d[1]*chunk2d[1], (rank2d[1]+1)*chunk2d[1]+1)]

x = np.linspace(-.5, .5, dim2d[0])
y = np.linspace(-.5, .5, dim2d[1])
x = x[slice2d[0]]
y = y[slice2d[1]]
vx, vy = np.meshgrid(x, y)

r = rotmat(-.45*np.pi, [0, 2], 3)@rotmat(0.1*np.pi, [0, 1], 3)
v_xyz = np.stack(lonlat2xyz(vx, vy), axis=-1)
v_xyz = v_xyz@r
c_xyz = 0.25*(v_xyz[:-1, :-1, :]+v_xyz[1:, :-1, :]+v_xyz[:-1, 1:, :]+v_xyz[1:, 1:, :])
vx, vy = xyz2lonlat(v_xyz[..., 0], v_xyz[..., 1], v_xyz[..., 2])
cx, cy = xyz2lonlat(c_xyz[..., 0], c_xyz[..., 1], c_xyz[..., 2])

grid = Curve2dGrid(f"{my_comp_name}_grid", vx, vy)

corner_idx = np.arange(math.prod(dim2d)).reshape(dim2d)[slice2d[0], slice2d[1]]
grid.set_global_index(corner_idx.T.flatten(), Location.CORNER)

cell_idx = np.arange(math.prod([d-1 for d in dim2d])).reshape([d-1 for d in dim2d])[
    slice(rank2d[0]*chunk2d[0], (rank2d[0]+1)*chunk2d[0]),
    slice(rank2d[1]*chunk2d[1], (rank2d[1]+1)*chunk2d[1])]
grid.set_global_index(cell_idx.T.flatten(), Location.CELL)

points_cell = grid.def_points(Location.CELL, cx, cy)
points_vertex = grid.def_points(Location.CORNER, vx, vy)

interpolations = [("conserv", Location.CELL),
                  ("2nd_conserv", Location.CELL),
                  ("avg", Location.CORNER),
                  ("hcsbb", Location.CELL),
                  ("rbf", Location.CELL)]

MPI.COMM_WORLD.Barrier()
MPI.COMM_WORLD.Barrier()

field = {f"{interp[0]}_{comp_name}_field_out":
         Field.create(f"{interp[0]}_{comp_name}_field_out", comp,
                      points_vertex if interp[1] == Location.CORNER else points_cell,
                      1, "2", TimeUnit.SECOND)
         for comp_name, interp in product(yac.component_names, interpolations)}

MPI.COMM_WORLD.Barrier()
MPI.COMM_WORLD.Barrier()
yac.enddef()
MPI.COMM_WORLD.Barrier()


def data(lon, lat):
    return (np.sin(8*lon)*np.cos(12*lat)).flatten()


data_vertex = data(vx, vy)
data_cell = data(cx, cy)

MPI.COMM_WORLD.Barrier()
for interp in interpolations:
    field[f"{interp[0]}_{my_comp_name}_field_out"].put(data_vertex if interp[1] == Location.CORNER else data_cell)

for comp_name, interp in product(yac.component_names, interpolations):
    if comp_name == my_comp_name:
        continue
    data = field[f"{interp[0]}_{comp_name}_field_out"].get()
MPI.COMM_WORLD.Barrier()
