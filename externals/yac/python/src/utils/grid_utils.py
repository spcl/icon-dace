from yac import *
from netCDF4 import Dataset
import numpy as np

def read_grid(filename, grid_name, comm=None):
    dataset = Dataset(filename, "r")
    nv = dataset.dimensions["nv"].size
    index_base = np.min(dataset["vertex_of_cell"])

    if comm is not None and comm.size > 1:
        domain_idx = np.flatnonzero(dataset["cell_no_of_domains"]==comm.size)[0]
        cell_idx = np.flatnonzero(dataset["cell_domain_id"][:,domain_idx] == comm.rank)
        vertex_idx = np.flatnonzero(dataset["vertex_of_cell"][:,cell_idx].flatten())
        no_cells = len(cell_idx)
    else:
        cell_idx = slice(None)
        vertex_idx = slice(None)
        no_cells = dataset.dimensions["cell"].size

    grid = UnstructuredGrid(grid_name, np.ones(no_cells)*nv,
                            dataset["vlon"][vertex_idx], dataset["vlat"][vertex_idx],
                            dataset["vertex_of_cell"][:,cell_idx].T.flatten() - index_base)

    corner_points = grid.def_points(Location.CORNER,
                                    dataset["vlon"][vertex_idx],
                                    dataset["vlat"][vertex_idx])

    if "clon" in dataset.variables:
        cell_points = grid.def_points(Location.CELL,
                                      dataset["clon"][cell_idx],
                                      dataset["clat"][cell_idx])
    else:
        cell_points = None

    grid.corner_points = corner_points
    grid.cell_points = cell_points
    return grid, cell_idx, vertex_idx
