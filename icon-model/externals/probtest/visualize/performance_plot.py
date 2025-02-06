"""
CLI for generating and saving performance plots

This script creates performance plots for various timing names from a timing
database. The performance data is visualized as time series graphs, where each
timing name is plotted alongside its associated timing data. The script also
highlights revision periods on the plot for easier comparison.

Usage:
    python script_name.py --timing-database <timing_database> --savedir
    <savedir> --timing-names <timing_names> --experiment-name <experiment_name>
    [--i-table <i_table>]
"""

from datetime import datetime, timedelta

import click
import matplotlib
import numpy as np
from matplotlib import pyplot as plt

from util.click_util import CommaSeperatedStrings, cli_help
from util.constants import DATETIME_FORMAT
from util.log_handler import logger
from util.tree import TimingTree
from util.utils import first_idx_of, last_idx_of, unique_elements


def colour_revs(times, revs, ax):
    urevs = unique_elements(revs)
    srevs = [r[:8] for r in revs]
    first_last = [(first_idx_of(revs, r), last_idx_of(revs, r)) for r in urevs]
    ntot = len(times)

    i = 0
    for f, l in first_last:
        c = "gray" if i % 2 == 0 else "white"
        t0 = matplotlib.dates.date2num(datetime.strptime(times[f], DATETIME_FORMAT))
        t1 = matplotlib.dates.date2num(
            datetime.strptime(times[min(ntot - 1, l + 1)], DATETIME_FORMAT)
        )

        r = srevs[f]
        assert revs[f] == revs[l]

        bottom, top = ax.get_ylim()
        loc = (t0, bottom) if i % 2 == 0 else (t0, top)
        va = "bottom" if i % 2 == 0 else "top"

        ax.axvspan(t0, t1, alpha=0.3, color=c)
        ax.annotate(r, loc, rotation=90, alpha=0.3, va=va)
        i += 1


@click.command()
@click.option(
    "--timing-database",
    help=cli_help["timing_database"],
)
@click.option(
    "--savedir",
    help=cli_help["savedir"],
)
@click.option(
    "--timing-names",
    type=CommaSeperatedStrings(),
    help=cli_help["timing_names"],
)
@click.option(
    "--experiment-name",
    help=cli_help["experiment_name"],
)
@click.option("--i-table", type=int, help=cli_help["i_table"], default=-1)
def performance_plot(timing_database, savedir, timing_names, experiment_name, i_table):
    tt = TimingTree.from_json(timing_database)

    dates = tt.get_sorted_finish_times()

    x = matplotlib.dates.date2num(dates)

    last_measurement = dates[-1]
    week = timedelta(days=7)
    day = timedelta(days=1)
    xticks_dt = [last_measurement + day - week * n for n in range(20)]
    xticks = matplotlib.dates.date2num(xticks_dt)
    xticklabels = [d.strftime("%d-%m-%Y") for d in xticks_dt]

    for timer in timing_names:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(10, 6))

        ax.set_xlim(left=xticks[-1], right=xticks[0])
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, rotation=90)
        ax.set_ylabel("runtime / s")

        t = tt.find(timer, i_table)
        y = tt.data[i_table].loc[(timer, slice(None)), "total max (s)"].values
        ymax, ymin = np.max(y), np.min(y)
        ax.plot(x, y, "k", label=timer, lw=2, linestyle="solid", marker=None)
        for c in t.children:
            y = (
                tt.data[i_table]
                .loc[(c.get_name(), slice(None)), "total max (s)"]
                .values
            )
            ymax, ymin = max(ymax, np.max(y)), min(ymin, np.min(y))

            ax.plot(x, y, label=c.get_name(), linestyle="solid", marker=None)

        fig.tight_layout()
        colour_revs(tt.meta_data["finish_time"], tt.meta_data["revision"], ax)

        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8]
        )

        # Put a legend below current axis
        ax.legend(
            loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=5, frameon=False
        )

        if savedir:
            path = f"{savedir}/perf_{experiment_name}_{timer}.png"
            logger.info("saving figure to %s", path)
            fig.savefig(path)
