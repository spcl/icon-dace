"""
CLI for plotting and checking model results against tolerance levels

This module provides functionality for:
- Deleting uninitialized fields from DataFrames.
- Generating plots to visualize the relative differences between model outputs
  and reference data, compared against specified tolerance levels.
"""

from pathlib import Path

import click
import numpy as np
from matplotlib import lines as mlines
from matplotlib import pyplot as plt

from util.click_util import cli_help
from util.constants import compute_statistics
from util.dataframe_ops import compute_rel_diff_dataframe, parse_check
from util.log_handler import logger
from util.visualize_util import create_figure


def delete_uninitialized_fields(df):
    mask = (df == 0) | df.isnull()
    uninitialized = df.index[mask.all(axis=1)]
    df_filtered = df.drop(uninitialized)
    return df_filtered


prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


@click.command()
@click.option(
    "--tolerance-file-name",
    help=cli_help["tolerance_file_name"],
)
@click.option(
    "--input-file-ref",
    help=cli_help["input_file_ref"],
)
@click.option(
    "--input-file-cur",
    help=cli_help["input_file_cur"],
)
@click.option(
    "--factor",
    type=float,
    help=cli_help["factor"],
)
@click.option(
    "--savedir",
    help=cli_help["savedir"],
)
def check_plot(tolerance_file_name, input_file_ref, input_file_cur, factor, savedir):

    df_tol, df_ref, df_cur = parse_check(
        tolerance_file_name, input_file_ref, input_file_cur, factor
    )

    # compute relative difference
    diff_df = compute_rel_diff_dataframe(df_ref, df_cur)
    # take maximum over height
    diff_df = diff_df.groupby(["file_ID", "variable"]).max()
    # Delete uninitialized fields as they are not interesting for plotting
    diff_df = delete_uninitialized_fields(diff_df)
    error_relative_to_tolerance = diff_df - df_tol

    # Select a subset of variables if there are too many variables
    max_nsubplots = 30
    if len(diff_df.index.values) > max_nsubplots:
        selected_keys = (
            error_relative_to_tolerance.max(axis=1, skipna=True)
            .nlargest(max_nsubplots)
            .keys()
        )
        diff_df = diff_df.loc[selected_keys]
        print(
            f"Limit the subplots to the graphs with the {max_nsubplots}"
            " largest error relative to tolerance"
        )

    times = np.array(diff_df.columns.levels[0])

    max_ncol = 3
    ncols = min(len(diff_df.index.values), max_ncol)
    nrows = int(np.ceil(len(diff_df.index.values) / ncols))
    legend_height = 2  # height in inches
    fig_height = 3.5 * nrows + legend_height

    fig, ax = create_figure(ncols, nrows, fig_width=5 * ncols, fig_height=fig_height)

    for i, e in enumerate(diff_df.index.values):
        for j, c in enumerate(compute_statistics):
            fid = e[0]
            vn = e[1]
            # TODO: why does matplotlib not understand panda nans?
            y_tol = df_tol.loc[(fid, vn), (slice(None), c)].values
            nn_tol = ~np.isnan(y_tol)
            y_diff = diff_df.loc[(fid, vn), (slice(None), c)].values
            nn_diff = ~np.isnan(y_diff)
            ax[i // ncols, i % ncols].semilogy(
                times[nn_tol], y_tol[nn_tol], c=colors[j]
            )
            ax[i // ncols, i % ncols].semilogy(
                times[nn_diff], y_diff[nn_diff], ls="--", c=colors[j]
            )
            ax[i // ncols, i % ncols].set_ylim(bottom=1e-15, top=1)
            ax[i // ncols, i % ncols].set_xlabel("timestep")
            ax[i // ncols, i % ncols].set_title(f"{fid}: {vn}")
        if i % ncols == 0:
            ax[i // ncols, 0].set_ylabel("relative error")

    lines = []
    for i, c in enumerate(compute_statistics):
        lines.append(mlines.Line2D([], [], color=colors[i], label=c))
    lines.append(mlines.Line2D([], [], color="k", label="tolerance"))
    lines.append(mlines.Line2D([], [], color="k", ls="--", label="model"))

    legcol = max_ncol if ncols == 1 else 5
    ax[-1, ncols // 2].legend(
        handles=lines, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=legcol
    )
    fig.tight_layout()
    fig.subplots_adjust(
        bottom=legend_height / fig_height
    )  # use an absolute margin for the legend
    if savedir:
        Path(savedir).mkdir(exist_ok=True, parents=True)
        path = f"{savedir}/check_plot.png"
        logger.info("saving figure to %s", path)
        fig.savefig(path)
    else:
        plt.show()
