#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import re
import cartopy.crs as ccrs
from pathlib import Path


def to_numpy_deg(var):
    if "rad" in var.units:
        return np.rad2deg(var)
    else:
        return var   # assuming rad else


def lonlat2xyz(lon, lat):
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)
    clat = np.cos(lat)
    slat = np.sin(lat)
    return np.stack([clat * np.cos(lon), clat * np.sin(lon), slat], axis=-1)


def xyz2lonlat(xyz):
    xyz = np.array(xyz)
    lat = np.arcsin(xyz[..., 2])
    lon = np.arctan2(xyz[..., 1], xyz[..., 0])
    return np.rad2deg(lon), np.rad2deg(lat)


def spherical_mean(lon, lat):
    m = np.mean(lonlat2xyz(lon, lat), axis=1)
    norm_inv = 1./np.linalg.norm(m, axis=-1)
    return xyz2lonlat(norm_inv[:, None] * m)


class GridBase:
    """
    Base class for Grids
    """
    def __init__(self):
        self.transform = ccrs.Geodetic()

    def plot(self, ax, label=None, plot_kwargs={}):
        edge_vertices = self.edge_vertices
        x = self.vlon[edge_vertices]
        y = self.vlat[edge_vertices]

        # exclude edges outside the extent
        extent = ax.get_extent(crs=ccrs.PlateCarree())
        mask = ((x[0, :] >= extent[0]) * (x[0, :] <= extent[1]) *
                (y[0, :] >= extent[2]) * (y[0, :] <= extent[3])) + \
               ((x[1, :] >= extent[0]) * (x[1, :] <= extent[1]) *
                (y[1, :] >= extent[2]) * (y[1, :] <= extent[3]))

        ax.plot(x[:, mask], y[:, mask], transform=self.transform,
                **plot_kwargs)

        text_kwargs = {}
        if "color" in plot_kwargs:
            text_kwargs["color"] = plot_kwargs["color"]
        if "zorder" in plot_kwargs:
            text_kwargs["zorder"] = plot_kwargs["zorder"]

        def add_label_in_extent(x, y, txt):
            mask = ((x >= extent[0]) * (x <= extent[1]) *
                    (y >= extent[2]) * (y <= extent[3]))
            for xx, yy, tt in zip(x[mask], y[mask], txt[mask]):
                ax.text(xx, yy, tt, ha="center", va="center", transform=ccrs.PlateCarree(),
                        **text_kwargs)

        if label == "cell":
            add_label_in_extent(self.clon, self.clat,
                                np.arange(len(self.clon))+self.idx_offset)
        elif label == "edge":
            add_label_in_extent(self.elon, self.elat,
                                np.arange(len(self.elon))+self.idx_offset)
        elif label == "vertex":
            add_label_in_extent(self.vlon, self.vlat,
                                np.arange(len(self.vlon))+self.idx_offset)

    def save_netcdf(self, filename):
        with Dataset(filename, "w") as dataset:
            vdim = dataset.createDimension("vertex", self.vlon.shape[0])
            vlon = dataset.createVariable("vlon", "f8", dimensions=(vdim,))
            vlon[:] = np.deg2rad(self.vlon)
            vlon.units = "radian"
            vlon.long_name = "vertex longitude"
            vlon.standard_name = "grid_longitude"
            vlat = dataset.createVariable("vlat", "f8", dimensions=(vdim,))
            vlat.units = "radian"
            vlat.long_name = "vertex latitude"
            vlat.standard_name = "grid_latitude"
            vlat[:] = np.deg2rad(self.vlon)
            cdim = dataset.createDimension("cell", self.vertex_of_cell.shape[1])
            nv = dataset.createDimension("nv", self.vertex_of_cell.shape[0])
            vertex_of_cell = dataset.createVariable("vertex_of_cell", "i", dimensions=(nv,cdim))
            vertex_of_cell[:] = self.vertex_of_cell
            vertex_of_cell.long_name = "vertices of each cell"
            if hasattr(self, "clon") and hasattr(self, "clat"):
                clon = dataset.createVariable("clon", "f8", dimensions=(cdim,))
                clon.units = "radian"
                clon.long_name = "center longitude"
                clon.standard_name = "grid_longitude"
                clon[:] = np.deg2rad(self.clon)
                clat = dataset.createVariable("clat", "f8", dimensions=(cdim,))
                clat.units = "radian"
                clat.long_name = "center latitude"
                clat.standard_name = "grid_latitude"
                clat[:] = np.deg2rad(self.clon)
            if hasattr(self, "elon") and hasattr(self, "elat"):
                edim = dataset.createDimension("edge", self.elon.shape[0])
                elon = dataset.createVariable("elon", "f8", dimensions=(edim,))
                elon[:] = self.elon
                elat = dataset.createVariable("elat", "f8", dimensions=(edim,))
                elat[:] = self.elon


class ICONGrid(GridBase):
    """
    Read grid from netCDF file like it is stored in ICON with fields
    clon, clat, vlon, vlat, vertex_of_cell etc.
    """
    def __init__(self, string):
        super().__init__()
        filename = Path(string[5:]).expanduser()  # truncate "icon:"
        with Dataset(filename, "r") as dataset:
            self.no_cells = dataset.dimensions["cell"].size
            self.clon = to_numpy_deg(dataset["clon"])
            self.clat = to_numpy_deg(dataset["clat"])
            self.vlon = to_numpy_deg(dataset["vlon"])
            self.vlat = to_numpy_deg(dataset["vlat"])
            self.elon = to_numpy_deg(dataset["elon"])
            self.elat = to_numpy_deg(dataset["elat"])
            index_base = np.min(dataset["vertex_of_cell"])
            self.vertex_of_cell = dataset["vertex_of_cell"] - index_base
            self.edge_vertices = dataset["edge_vertices"] - index_base
            self.adjacent_cell_of_edge = dataset["adjacent_cell_of_edge"] - index_base
            self.idx_offset = index_base


class GaussianGrid(GridBase):
    """
    Generate Gaussian grid with an given resolution and extent
    (currently only on an closed interval in lon)
    """
    def __init__(self, string):
        super().__init__()
        m = re.match(r"([gG])([0-9]+),([0-9]+)(,(-?[0-9]+(\.[0-9]+)?),(-?[0-9]+(\.[0-9]+)?),(-?[0-9]+(\.[0-9]+)?),(-?[0-9]+(\.[0-9]+)?))?", string)
        if not m:
            raise Exception(f"Defintion of Gaussian grid {string} does not match regular expression.")
        self.idx_offset = 0 if m.groups()[0] == "g" else 1
        if m.groups()[3] is not None:
            extent = [float(m.groups()[4]), float(m.groups()[6]), float(m.groups()[8]), float(m.groups()[10])]
        else:
            extent = [-180, 180, -90, 90]
        res = [int(m.groups()[1]), int(m.groups()[2])]

        self.lon = np.linspace(extent[0], extent[1], res[0]+1, endpoint=True)
        self.lat = np.linspace(extent[2], extent[3], res[1]+1, endpoint=True)[::-1]
        self.vlon, self.vlat = [x.flatten() for x in np.meshgrid(self.lon, self.lat)]
        self.clon, self.clat = [x.flatten() for x in
                                np.meshgrid(self.lon[:-1] + 0.5*np.abs(self.lon[1]-self.lon[0]),
                                            self.lat[1:] + 0.5*np.abs(self.lat[1]-self.lat[0]))]
        self.no_cells = len(self.clat)
        yidx, xidx = np.mgrid[0:res[1], 0:res[0]]

        self.vertex_of_cell = np.stack([(xidx     + yidx*(res[0]+1)).flatten(),
                                        ((xidx+1) + yidx*(res[0]+1)).flatten(),
                                        (xidx     + (yidx+1)*(res[0]+1)).flatten(),
                                        ((xidx+1) + (yidx+1)*(res[0]+1)).flatten()])

        yidx, xidx = np.mgrid[0:res[1]+1, 0:res[0]]
        hedges = np.stack([(xidx     + yidx*(res[0]+1)).flatten(),
                           ((xidx+1) + yidx*(res[0]+1)).flatten()])
        yidx, xidx = np.mgrid[0:res[1], 0:res[0]+1]
        vedges = np.stack([(xidx + yidx    *(res[0]+1)).flatten(),
                           (xidx + (yidx+1)*(res[0]+1)).flatten()])
        self.edge_vertices = np.hstack([hedges, vedges])
        self.elon = 0.5*np.sum(self.vlon[self.edge_vertices], axis=0)
        self.elat = 0.5*np.sum(self.vlat[self.edge_vertices], axis=0)
        self.transform = ccrs.PlateCarree()


class HealPixGrid(GridBase):
    """
    HealPix Grid
    """
    def __init__(self, string):
        super().__init__()
        import healpy
        m = re.match(r"([hH]([0-9]+))([rn])?", string)
        if not m:
            raise Exception(f"Defintion of HealPix grid {string} does not match regular expression.")
        self.idx_offset = 0 if m.groups()[0] == "h" else 1
        self.nside = int(m.groups()[1])
        self.nest = m.groups()[2] != "r"
        ncells = healpy.pixelfunc.nside2npix(self.nside)
        centers_xyz = np.stack(
            healpy.pixelfunc.pix2vec(self.nside, range(ncells), nest=self.nest),
            axis=-1,
        )
        self.clon, self.clat = xyz2lonlat(centers_xyz)
        boundaries_xyz = (
            healpy.boundaries(self.nside, range(ncells), nest=self.nest)
            .transpose(0, 2, 1)
            .reshape(-1, 3)
        )
        verts_xyz, quads = np.unique(boundaries_xyz, return_inverse=True, axis=0)
        self.vlon, self.vlat = xyz2lonlat(verts_xyz)
        self.vertex_of_cell = quads.reshape(-1, 4).T
        edges = np.hstack([self.vertex_of_cell[(i, (i+1) % 4), :] for i in range(4)])
        edges = np.sort(edges, axis=0)
        self.edge_vertices = np.unique(edges, axis=1)


class DualGrid(GridBase):
    def __init__(self, string):
        super().__init__()
        self.simple_dual = string[0] == "D"
        prim_grid = get_grid(string[5:])  # truncate "dual:"
        self.clon = prim_grid.vlon
        self.clat = prim_grid.vlat
        # compute corners
        if hasattr(prim_grid, "clon") and hasattr(prim_grid, "clat"):
            elem_centers_lon, elem_centers_lat = prim_grid.clon, prim_grid.clat
        else:
            elem_centers_lon, elem_centers_lat = spherical_mean(
                prim_grid.vlon[prim_grid.vertex_of_cell], prim_grid.vlat[prim_grid.vertex_of_cell])

        boundary_edges = (prim_grid.adjacent_cell_of_edge < 0)
        is_boundary = sum(boundary_edges) != 0
        boundary_nodes, bnodes_inverse = np.unique(
            prim_grid.edge_vertices[:, boundary_edges[0, :] + boundary_edges[1, :]], return_inverse=True)
        bnodes_inverse = bnodes_inverse.reshape(2, -1)
        self.idx_offset = 1

        if self.simple_dual:
            if hasattr(prim_grid, "elon") and hasattr(prim_grid, "elat"):
                edge_centers_lon, edge_centers_lat = prim_grid.elon[is_boundary], prim_grid.elat[is_boundary]
            else:
                edge_centers_lon, edge_centers_lat = spherical_mean(
                    prim_grid.vlon[prim_grid.edge_vertices[:, is_boundary].T], prim_grid.vlat[prim_grid.edge_vertices[:, is_boundary].T])
        else:
            if hasattr(prim_grid, "elon") and hasattr(prim_grid, "elat"):
                edge_centers_lon, edge_centers_lat = prim_grid.elon, prim_grid.elat
            else:
                edge_centers_lon, edge_centers_lat = spherical_mean(
                    prim_grid.vlon[prim_grid.edge_vertices.T], prim_grid.vlat[prim_grid.edge_vertices.T])

        nelem = len(elem_centers_lon)
        nedges = len(edge_centers_lon)
        self.vlon = np.concatenate([edge_centers_lon, elem_centers_lon, prim_grid.vlon[boundary_nodes]])
        self.vlat = np.concatenate([edge_centers_lat, elem_centers_lat, prim_grid.vlat[boundary_nodes]])

        if self.simple_dual:
            self.edge_vertices = np.hstack([
                nedges + prim_grid.adjacent_cell_of_edge[:, ~is_boundary],  # connect adjacent cell centers of non-boundary edges
                nedges + nelem + bnodes_inverse,
                np.stack([range(sum(is_boundary)), nedges + prim_grid.adjacent_cell_of_edge[boundary_edges[0, is_boundary].astype(int), is_boundary]])
            ])
            self.adjacent_cell_of_edge = np.hstack([
                prim_grid.edge_vertices[:, ~is_boundary],
                np.stack([boundary_nodes[bnodes_inverse[0, :]], -1*np.ones(bnodes_inverse[0, :].shape, dtype=int)]),
                prim_grid.edge_vertices[:, is_boundary],
            ])
        else:
            self.edge_vertices = np.hstack(
                [np.stack([np.flatnonzero(~boundary_edges[0, :]), nedges+prim_grid.adjacent_cell_of_edge[0, ~boundary_edges[0, :]]]),
                 np.stack([np.flatnonzero(~boundary_edges[1, :]), nedges+prim_grid.adjacent_cell_of_edge[1, ~boundary_edges[1, :]]]),
                 np.stack([np.flatnonzero(is_boundary), nedges + nelem + bnodes_inverse[0, :]]),
                 np.stack([np.flatnonzero(is_boundary), nedges + nelem + bnodes_inverse[1, :]])
                 ])
            self.adjacent_cell_of_edge = np.hstack([
                    prim_grid.edge_vertices[:, ~boundary_edges[0, :]],
                    prim_grid.edge_vertices[:, ~boundary_edges[1, :]],
                    np.stack([boundary_nodes[bnodes_inverse[0, :]], -1*np.ones(bnodes_inverse[0, :].shape, dtype=int)]),
                    np.stack([boundary_nodes[bnodes_inverse[1, :]], -1*np.ones(bnodes_inverse[1, :].shape, dtype=int)]),
                ])


class FesomGrid(GridBase):
    def __init__(self, string):
        super().__init__()
        path = Path(string[6:]).expanduser()
        nodes = np.fromfile(path / "nod2d.out", sep=" ")[1:].reshape(-1, 4)
        self.vertex_of_cell = np.fromfile(path / "elem2d.out", sep=" ", dtype=int)[1:].reshape(-1, 3).T - 1
        self.edge_vertices = np.fromfile(path / "edges.out", sep=" ", dtype=int).reshape(-1, 2).T - 1
        self.vlon = (((180 + nodes[:, 1]) % 360) - 180)
        self.vlat = nodes[:, 2]
        self.clon, self.clat = spherical_mean(self.vlon[self.vertex_of_cell].T, self.vlat[self.vertex_of_cell].T)
        self.adjacent_cell_of_edge = np.fromfile(path / "edge_tri.out", sep=" ", dtype=int).reshape(-1, 2).T - 1


class ScripGrid(GridBase):
    def __init__(self, string):
        super().__init__()
        path = string[6:].split(':')[0]
        grid = string[6:].split(':')[1]
        with Dataset(path, "r") as dataset:
            x_size = dataset.dimensions["x_"+grid].size
            y_size = dataset.dimensions["y_"+grid].size
            self.no_cells = x_size * y_size
            self.no_crn = dataset.dimensions["crn_"+grid].size
            self.clon = np.asfortranarray(dataset[grid+".lon"]).flatten()
            self.clat = np.asfortranarray(dataset[grid+".lat"]).flatten()
            self.vlon = np.asfortranarray(dataset[grid+".clo"]).flatten()
            self.vlat = np.asfortranarray(dataset[grid+".cla"]).flatten()

        self.vertex_of_cell = np.arange(self.no_cells*self.no_crn).reshape(self.no_crn,-1)
        edges = np.hstack([self.vertex_of_cell[(i, (i+1) % self.no_crn), :] for i in range(self.no_crn)])
        edges = np.sort(edges, axis=0)
        self.edge_vertices = np.unique(edges, axis=1)
        self.elon = 0.5*np.sum(self.vlon[self.edge_vertices], axis=0)
        self.elat = 0.5*np.sum(self.vlat[self.edge_vertices], axis=0)
        self.idx_offset = 0


def get_grid(grid):
    if grid.lower().startswith("icon:"):
        return ICONGrid(grid)
    if grid.lower().startswith("g"):
        return GaussianGrid(grid)
    if grid.lower().startswith("h"):
        return HealPixGrid(grid)
    if grid.lower().startswith("fesom:"):
        return FesomGrid(grid)
    if grid.lower().startswith("dual:"):
        return DualGrid(grid)
    if grid.lower().startswith("scrip:"):
        return ScripGrid(grid)
    else:
        raise Exception(f"cannot parse grid description: {grid}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="grid_utils.py")
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')
    parser_save_grid = subparsers.add_parser('save_grid', help='save a grid to netCDF file')
    parser_save_grid.add_argument('grid', type=str, help='grid')
    parser_save_grid.add_argument('filename', type=Path, help='filename')
    args = parser.parse_args()
    if args.command == "save_grid":
        get_grid(args.grid).save_netcdf(args.filename)
