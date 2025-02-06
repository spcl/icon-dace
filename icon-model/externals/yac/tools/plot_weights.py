#!/usr/bin/env python3
# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from grid_utils import get_grid
import logging
from logging import info


def plot_weights(ax, weightsfile, src_grid, tgt_grid,
                 src_idx=None, tgt_idx=None, zoom=1, quiver_kwargs={}):
    weights = Dataset(weightsfile, "r")
    if 'version' in weights.__dict__:
        info(f"{weights.version=}")
        if weights.version != "yac weight file 1.0":
            print("WARNING: You are using an incompatible weight file version\n" +
                  weights.version + " != yac weight file 1.0")
    num_src_fields = weights.num_src_fields if "num_src_fields" in weights.variables else 1
    info(f"{num_src_fields=}")
    for src_field in range(num_src_fields):
        if "src_locations" in weights.variables:
            src_locstr = bytes(np.array(weights["src_locations"][src_field, :])).decode().rstrip("\0")
        else:
            src_locstr = "CELL"
        info(f"source loc: {src_locstr}")
        if src_locstr == "CELL":
            src_points = np.stack([src_grid.clon, src_grid.clat])
        elif src_locstr == "CORNER":
            src_points = np.stack([src_grid.vlon, src_grid.vlat])
        elif src_locstr == "EDGE":
            src_points = np.stack([src_grid.vlon, src_grid.vlat])
        else:
            raise f"Unknown location string {src_locstr}"
    if "dst_locations" in weights.variables:
        tgt_locstr = bytes(np.array(weights["dst_location"])).decode().rstrip("\0")
    else:
        tgt_locstr = "CELL"
    info(f"target loc: {tgt_locstr}")
    if tgt_locstr == "CELL":
        tgt_points = np.stack([tgt_grid.clon, tgt_grid.clat])
    elif tgt_locstr == "CORNER":
        tgt_points = np.stack([tgt_grid.vlon, tgt_grid.vlat])
    elif tgt_locstr == "EDGE":
        tgt_points = np.stack([tgt_grid.elon, tgt_grid.elat])
    else:
        raise f"Unknown location string {tgt_locstr}"

    yac_weights_format_address_offset = 1
    src_adr = np.asarray(weights["src_address"])-yac_weights_format_address_offset-src_grid.idx_offset
    tgt_adr = np.asarray(weights["dst_address"])-yac_weights_format_address_offset-tgt_grid.idx_offset

    info(f"Source indices in weight file range from {np.min(src_adr)} to {np.max(src_adr)}")
    info(f"Number of source points is {src_points.shape[1]}")
    info(f"Target indices in weight file range from {np.min(tgt_adr)} to {np.max(tgt_adr)}")
    info(f"Number of target points is {tgt_points.shape[1]}")

    if np.max(src_adr) > src_points.shape[1]:
        raise Exception(f"Source grid too small. Max index in weight file is {np.max(src_adr)}, number of source points are {src_points.shape[1]}")
    if np.max(tgt_adr) > tgt_points.shape[1]:
        raise Exception(f"Target grid too small. Max index in weight file is {np.max(tgt_adr)}, number of source points are {tgt_points.shape[1]}")

    # Remove redundant links as in SCRIP bilinear and bicubic weights
    mask = (src_adr >= 0)

    # Restrain plot for targeted cells
    if src_idx is not None:
        mask *= (src_adr == (src_idx-src_grid.idx_offset))
    if tgt_idx is not None:
        mask *= (tgt_adr == (tgt_idx-tgt_grid.idx_offset))

    if not any(mask):
        raise Exception("All points are mask.")

    if src_idx is not None or tgt_idx is not None:
        info(f"Mask: {sum(mask)}/{len(mask)}")
        src_adr = src_adr[mask, ...]
        tgt_adr = tgt_adr[mask, ...]
        # transform to the projection that is used
        src_t = ax.projection.transform_points(ccrs.PlateCarree(),
                                               src_points[0, src_adr], src_points[1, src_adr])
        info(f"src: {src_points[0, src_adr], src_points[1, src_adr]} -> {src_t[:, 0], src_t[:, 1]}")
        tgt_t = ax.projection.transform_points(ccrs.PlateCarree(),
                                               tgt_points[0, tgt_adr], tgt_points[1, tgt_adr])
        info(f"tgt: {tgt_points[0, tgt_adr], tgt_points[1, tgt_adr]} -> {tgt_t[:, 0], tgt_t[:, 1]}")
        e = np.array([min(np.min(src_t[:, 0]), np.min(tgt_t[:, 0])),
                      max(np.max(src_t[:, 0]), np.max(tgt_t[:, 0])),
                      min(np.min(src_t[:, 1]), np.min(tgt_t[:, 1])),
                      max(np.max(src_t[:, 1]), np.max(tgt_t[:, 1]))])
        c = np.array([0.5*(e[0]+e[1]), 0.5*(e[0]+e[1]),
                      0.5*(e[2]+e[3]), 0.5*(e[2]+e[3])])
        extent = (e-c)*1.15*zoom + c
        info(f"{c=}")
        info(f"{e=}")
        info(f"{extent=}")
        ax.set_extent(extent, crs=ax.projection)
    else:
        # transform to the projection that is used
        src_t = ax.projection.transform_points(ccrs.PlateCarree(),
                                               src_points[0, src_adr], src_points[1, src_adr])
        tgt_t = ax.projection.transform_points(ccrs.PlateCarree(),
                                               tgt_points[0, tgt_adr], tgt_points[1, tgt_adr])

        extent = ax.get_extent()
        info(f"Extent: {extent}")
        # Finalize restriction for arrays
        mask *= (((src_t[:, 0] >= extent[0]) * (src_t[:, 0] <= extent[1]) *
                  (src_t[:, 1] >= extent[2]) * (src_t[:, 1] <= extent[3])) +
                 ((tgt_t[:, 0] >= extent[0]) * (tgt_t[:, 0] <= extent[1]) *
                  (tgt_t[:, 1] >= extent[2]) * (tgt_t[:, 1] <= extent[3])))
        info(f"Mask: {sum(mask)}/{len(mask)}")

        src_adr = src_adr[mask, ...]
        tgt_adr = tgt_adr[mask, ...]
        src_t = src_t[mask, ...]
        tgt_t = tgt_t[mask, ...]

    if sum(mask) > 3000:
        logging.warning(f"Trying to display a lot of points ({sum(mask)}). "
                        "This may take some time. "
                        "Consider a smaller --zoom parameter to reduce the number of points")

    uv = tgt_t - src_t

    c = weights["remap_matrix"][mask, 0]

    norm = matplotlib.colors.Normalize()
    cm = matplotlib.cm.Oranges  # decide for colormap
    ax.quiver(src_t[:, 0], src_t[:, 1],
              uv[:, 0], uv[:, 1], angles='xy', scale_units='xy', scale=1,
              color=cm(norm(c)),
              **quiver_kwargs)
    sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    clb = plt.colorbar(sm, ax=ax)
    clb.ax.zorder = -1

    if src_idx is not None:
        ax.set_title(f"Interpolation for source index {src_idx}\nSum of weights: {sum(c):.8f}")
    if tgt_idx is not None:
        ax.set_title(f"Interpolation for target index {tgt_idx}\nSum of weights: {sum(c):.8f}")
    # add pop ups if restricted to one cell:
    if src_idx is not None or tgt_idx is not None:
        # add label
        for x, y, t in zip(src_t[:, 0] + 0.5*uv[:, 0],
                           src_t[:, 1] + 0.5*uv[:, 1], c):
            ax.text(x, y, f"{t:.3}",
                    horizontalalignment='center', verticalalignment='center')

        wlines = ax.plot([src_t[:,0], tgt_t[:,0]],
                         [src_t[:,1], tgt_t[:,1]],
                         color='white', alpha=0.)
        wlabel = np.vstack([src_adr, tgt_adr, c.data])
        annotation = ax.annotate(text='', xy=(0, 0), xytext=(15, 15),
                                 textcoords='offset points',
                                 bbox={'boxstyle': 'round', 'fc': 'w'},
                                 arrowprops={'arrowstyle': '->'},
                                 zorder=9999)
        annotation.set_visible('False')

        def motion_hover(event):
            annotation_visible = annotation.get_visible()
            if event.inaxes == ax:
                if idx := next((idx+1 for idx, wl in enumerate(wlines) if wl.contains(event)[0]), None):
                    annotation.xy = (event.xdata, event.ydata)
                    text_label = f"src: {wlabel[0, idx-1] + src_grid.idx_offset:.0f} tgt: {wlabel[1, idx-1] + tgt_grid.idx_offset:.0f}\nweight: {wlabel[2, idx-1]:.3f}"
                    annotation.set_text(text_label)
                    annotation.set_visible(True)
                    ax.figure.canvas.draw_idle()
                else:
                    if annotation_visible:
                        annotation.set_visible(False)
                        ax.figure.canvas.draw_idle()

        ax.figure.canvas.mpl_connect('motion_notify_event', motion_hover)

    src_mask = {f"{src_locstr.lower()}_idx": src_adr}
    tgt_mask = {f"{tgt_locstr.lower()}_idx": tgt_adr}

    return src_mask, tgt_mask


def cell_extent(grid, idx, zoom=4):
    print(f"cell {idx}: ", grid.clon[idx-grid.idx_offset],
          grid.clat[idx-grid.idx_offset])
    vidx = grid.vertex_of_cell[:, idx-grid.idx_offset]
    e = np.array([np.min(grid.vlon[vidx]), np.max(grid.vlon[vidx]),
                  np.min(grid.vlat[vidx]), np.max(grid.vlat[vidx])])
    c = np.array([0.5*(e[0]+e[1]), 0.5*(e[0]+e[1]),
                  0.5*(e[2]+e[3]), 0.5*(e[2]+e[3])])
    return (e-c)*zoom + c


def main(source_grid, target_grid, weights_file, center=None,
         source_idx=None, target_idx=None, zoom=1,
         label_src_grid=None, label_tgt_grid=None,
         coast_res="50m", projection="orthografic",
         stencil_only=False,
         save_as=None, log_level=logging.INFO):
    logging.basicConfig(level=log_level)
    src_grid = get_grid(source_grid)
    if target_grid is not None:
        tgt_grid = get_grid(target_grid)
    else:
        tgt_grid = None
    fig = plt.figure(figsize=[10, 10])
    if center is None:
        center = [0, 0]
    if source_idx is not None:
        center = [src_grid.clon[source_idx-src_grid.idx_offset],
                  src_grid.clat[source_idx-src_grid.idx_offset]]
    elif target_idx is not None:
        center = [tgt_grid.clon[target_idx-tgt_grid.idx_offset],
                  tgt_grid.clat[target_idx-tgt_grid.idx_offset]]

    if projection == "orthographic":
        proj = ccrs.Orthographic(*center)
    elif projection == "stereographic":
        proj = ccrs.Stereographic(*center[::-1])
    elif projection == "platecarree":
        proj = ccrs.PlateCarree(center[0])
    info(f"{center=}")
    ax = fig.add_subplot(1, 1, 1, projection=proj)

    if weights_file is None:
        if source_idx is not None:
            extent = cell_extent(src_grid, source_idx, zoom)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
        elif target_idx is not None:
            extent = cell_extent(tgt_grid, target_idx, zoom)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
    if source_idx is None and target_idx is None:
        ax.set_extent([center[0]-1000000*zoom, center[0]+1000000*zoom,
                       center[1]-1000000*zoom, center[1]+1000000*zoom], crs=proj)

    # Put a background image on for nice sea rendering.
    ax.set_facecolor("#9ddbff")
    if coast_res:
        feature = cfeature.NaturalEarthFeature(name='land',
                                               category='physical',
                                               scale=coast_res,
                                               edgecolor='#000000',
                                               facecolor='#cbe7be')
        ax.add_feature(feature, zorder=-999)
    glin = ax.gridlines(draw_labels=True, alpha=0)
    if source_idx is not None or target_idx is not None:
        glin.top_labels = False

    if weights_file is not None:
        src_mask, tgt_mask = plot_weights(ax, weights_file, src_grid, tgt_grid,
                                          source_idx, target_idx, zoom, quiver_kwargs={"zorder": 2})

    # Plot grids
    if src_grid is not None:
        gridname = src_grid.gridname if hasattr(src_grid, 'gridname') else type(src_grid).__name__
        src_grid.plot(ax, label=label_src_grid,
                      plot_kwargs={"color": "green",
                                   "zorder": 1,
                                   "linewidth": 3,
                                   "label": f"Source grid: {gridname}"},
                      key="S",
                      **(src_mask if stencil_only else {}))
    if tgt_grid is not None:
        gridname = tgt_grid.gridname if hasattr(tgt_grid, 'gridname') else type(tgt_grid).__name__
        tgt_grid.plot(ax, label=label_tgt_grid,
                      plot_kwargs={"color": "blue",
                                   "zorder": 2,
                                   "linewidth": 2,
                                   "label": f"Target grid: {gridname}"},
                      key="T",
                      **(tgt_mask if stencil_only else {}))

    # remove duplicate legend entry
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper center',
               fancybox=True, shadow=True, ncol=2)

    if save_as:
        plt.savefig(save_as)
    else:
        ax.text(0.0, 0.0, "Press 'S' or 'T' to show source or target grid indices",
                verticalalignment='bottom', horizontalalignment='left',
                transform=fig.transFigure)
        plt.show()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="plot_weights.py",
                                     description="""
                                     Plot grids and yac weights file.
                                     """)
    parser.add_argument("source_grid", type=str,
                        help="source grid (for an overview of grids and how they are specified here see grid_utils.py)")
    parser.add_argument("target_grid", type=str, nargs='?',
                        default=None,
                        help="target grid (for an overview of grids and how they are specified here see grid_utils.py)")
    parser.add_argument("weights_file", type=str, help="YAC weights file", nargs='?',
                        default=None)
    parser.add_argument("--center", "-c", type=float, nargs=2, help="center of the orthografic projection",
                        default=(0, 0), metavar=("LON", "LAT"))
    parser.add_argument("--source_idx", "-s", type=int,
                        help="index of source cell to focus")
    parser.add_argument("--target_idx", "-t", type=int,
                        help="index of target cell to focus")
    parser.add_argument("--zoom", "-z", type=float, default=1,
                        help="zoom around the cell")
    parser.add_argument("--label_src_grid", type=str, default=None,
                        choices=("vertex", "edge", "cell"),
                        help="Add labels at the source grid")
    parser.add_argument("--label_tgt_grid", type=str, default=None,
                        choices=("vertex", "edge", "cell"),
                        help="Add labels at the source grid")
    parser.add_argument("--coast_res", type=str, default="50m",
                        nargs='?',
                        choices=("10m", "50m", "110m"),
                        help="Resolution of coastlines (def 50m).\nOmit argument to disable coastlines.")
    parser.add_argument("--projection", type=str, default="orthographic", choices=("orthographic", "stereographic", "platecarree"),
                        nargs="?", help="Type of projection")
    parser.add_argument("--stencil_only", action="store_true")
    parser.add_argument("--save_as", type=argparse.FileType("wb"), help="Save to file instead of showing the figure")
    parser.add_argument("--log-level", default=logging.WARNING, type=lambda x: getattr(logging, x.upper()),
                        help="Configure the logging level.")
    args = parser.parse_args()
    main(**args.__dict__)
