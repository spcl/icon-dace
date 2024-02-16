#!/usr/bin/env python3

from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from grid_utils import get_grid


def plot_weights(fig, ax, weightsfile, src_grid, tgt_grid,
                 src_idx=None, tgt_idx=None, zoom=1, quiver_kwargs={}):
    weights = Dataset(weightsfile, "r")
    if 'version' in weights.__dict__:
        if weights.version != "yac weight file 1.0":
            print("WARNING: You are using an incompatible weight file version\n" +
                  weights.version + " != yac weight file 1.0")
    src_points = np.empty([2, 0])
    num_src_fields = weights.num_src_fields if "num_src_fields" in weights.variables else 1
    for src_field in range(num_src_fields):
        if "src_locations" in weights.variables:
            locstr = bytes(np.array(weights["src_locations"][src_field, :])).decode().rstrip("\0")
        else:
            locstr = "CELL"
        if locstr == "CELL":
            src_points = np.hstack([src_points, np.stack([src_grid.clon,
                                                          src_grid.clat])])
        elif locstr == "CORNER":
            src_points = np.hstack([src_points, np.stack([src_grid.vlon,
                                                          src_grid.vlat])])
        elif locstr == "EDGE":
            src_points = np.hstack([src_points, np.stack([src_grid.vlon,
                                                          src_grid.vlat])])
        else:
            raise f"Unknown location string {locstr}"
    if "dst_locations" in weights.variables:
        locstr = bytes(np.array(weights["dst_location"])).decode().rstrip("\0")
    else:
        locstr = "CELL"
    if locstr == "CELL":
        tgt_points = np.stack([tgt_grid.clon, tgt_grid.clat])
    elif locstr == "CORNER":
        tgt_points = np.stack([tgt_grid.vlon, tgt_grid.vlat])
    elif locstr == "EDGE":
        tgt_points = np.stack([tgt_grid.elon, tgt_grid.elat])
    else:
        raise f"Unknown location string {locstr}"

    yac_weights_format_address_offset = 1
    src_adr = np.asarray(weights["src_address"])-yac_weights_format_address_offset-src_grid.idx_offset
    tgt_adr = np.asarray(weights["dst_address"])-yac_weights_format_address_offset-tgt_grid.idx_offset

    # Remove redundant links as in SCRIP bilinear and bicubic weights
    mask = (src_adr >= 0)

    # Restrain plot for targeted cells
    if src_idx is not None:
        mask *= (src_adr == (src_idx-src_grid.idx_offset))
    if tgt_idx is not None:
        mask *= (tgt_adr == (tgt_idx-tgt_grid.idx_offset))

    if src_idx is not None or tgt_idx is not None:
        src_res = src_adr[mask]
        tgt_res = tgt_adr[mask]
        e = np.array([min(np.min(src_points[0, src_res]), np.min(tgt_points[0, tgt_res])),
                      max(np.max(src_points[0, src_res]), np.max(tgt_points[0, tgt_res])),
                      min(np.min(src_points[1, src_res]), np.min(tgt_points[1, tgt_res])),
                      max(np.max(src_points[1, src_res]), np.max(tgt_points[1, tgt_res]))])
        c = np.array([0.5*(e[0]+e[1]), 0.5*(e[0]+e[1]),
                      0.5*(e[2]+e[3]), 0.5*(e[2]+e[3])])
        extent = (e-c)*1.15*zoom + c
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    else:
        extent = ax.get_extent(crs=ccrs.PlateCarree())

    # Finalize restriction for arrays
    mask *= (((src_points[0, src_adr] >= extent[0]) * (src_points[0, src_adr] <= extent[1]) *
              (src_points[1, src_adr] >= extent[2]) * (src_points[1, src_adr] <= extent[3])) +
             ((tgt_points[0, tgt_adr] >= extent[0]) * (tgt_points[0, tgt_adr] <= extent[1]) *
              (tgt_points[1, tgt_adr] >= extent[2]) * (tgt_points[1, tgt_adr] <= extent[3])))

    src_adr = src_adr[mask]
    tgt_adr = tgt_adr[mask]

    ax_proj = ax.projection
    src_t = ax_proj.transform_points(ccrs.PlateCarree(),
                                     src_points[0, src_adr], src_points[1, src_adr])
    tgt_t = ax_proj.transform_points(ccrs.PlateCarree(),
                                     tgt_points[0, tgt_adr], tgt_points[1, tgt_adr])
    uv = tgt_t-src_t

    c = weights["remap_matrix"][mask, 0]

    norm = matplotlib.colors.Normalize()
    cm = matplotlib.cm.Oranges  # decide for colormap
    ax.quiver(src_t[:, 0], src_t[:, 1],
              uv[:, 0], uv[:, 1], angles='xy', scale_units='xy', scale=1,
              color=cm(norm(c)),
              **quiver_kwargs)
    sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax)

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

        wlines = ax.plot([src_points[0, src_adr], tgt_points[0, tgt_adr]],
                         [src_points[1, src_adr], tgt_points[1, tgt_adr]],
                         color='white', alpha=0.,
                         transform=ccrs.PlateCarree())
        wcenter = np.vstack([0.5*(src_points[0, src_adr]+tgt_points[0, tgt_adr]),
                             0.5*(src_points[1, src_adr]+tgt_points[1, tgt_adr])])
        wlabel = np.vstack([src_adr, tgt_adr, c.data])
        annotation = ax.annotate(text='', xy=(0, 0), xytext=(15, 15),
                                 textcoords='offset points',
                                 bbox={'boxstyle': 'round', 'fc': 'w'},
                                 arrowprops={'arrowstyle': '->'},
                                 transform=ccrs.PlateCarree(),
                                 zorder=9999)
        annotation.set_visible('False')

        def motion_hover(event):
            annotation_visible = annotation.get_visible()
            if event.inaxes == ax:
                is_on = False
                for idx, wl in enumerate(wlines):
                    if wl.contains(event)[0]:
                        is_on = True
                        break
                if is_on:
                    annotation.xy = (wcenter[0, idx], wcenter[1, idx])
                    text_label = 'src: {} tgt: {}\nweight: {:.3f}'.format(int(wlabel[0, idx]),
                                                                          int(wlabel[1, idx]),
                                                                          wlabel[2, idx])
                    annotation.set_text(text_label)
                    annotation.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if annotation_visible:
                        annotation.set_visible(False)
                        fig.canvas.draw_idle()

        fig.canvas.mpl_connect('motion_notify_event', motion_hover)


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
         coast_res="50m", save_as=None):
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

    proj = ccrs.Orthographic(*center)
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

    if weights_file is not None:
        plot_weights(fig, ax, weights_file, src_grid, tgt_grid,
                     source_idx, target_idx, zoom, quiver_kwargs={"zorder": 2})

    # Plot grids
    if src_grid is not None:
        src_grid.plot(ax, label=label_src_grid, plot_kwargs={"color": "green", "zorder": 1})
    if tgt_grid is not None:
        tgt_grid.plot(ax, label=label_tgt_grid, plot_kwargs={"color": "blue", "zorder": 1})

    if save_as:
        plt.savefig(save_as)
    else:
        plt.show()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog="plot_weights.py",
                                     description="""
                                     plot grids and yac weights file.

                                     Grids can be specified either a filename.
                                     This is iterpreted like an ICON grid file or
                                     a string like "g360,180" from which a structured
                                     grid is generated. Use a capital G for 1-based
                                     indexing and a small g for 0-basid indexing
                                     followed by the resolution (lon,lat). Optionally
                                     you can add an extent by adding further 4 numbers
                                     (min_lon, max_lon, min_lat, max_lat). E.g.
                                     g100,100,-50,-45,-5,5
                                     """)
    parser.add_argument("source_grid", type=str,
                        help="source grid (netCDF file or [gG]lon,lat[,min_lon,max_lon,min_lat,max_lat])")
    parser.add_argument("target_grid", type=str, nargs='?',
                        default=None,
                        help="target grid (netCDF file or [gG]lon,lat[,min_lon,max_lon,min_lat,max_lat])")
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
    parser.add_argument("--save_as", type=str, help="Save to file instead of showing the figure")
    args = parser.parse_args()
    main(**args.__dict__)
