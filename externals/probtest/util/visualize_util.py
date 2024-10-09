"""
This module provides a functions for visualization.
"""

import numpy as np
from matplotlib import pyplot as plt


def create_figure(ncols, nrows, fig_width, fig_height):
    fig, ax = plt.subplots(
        ncols=ncols,
        nrows=nrows,
        figsize=(fig_width, fig_height),
        sharey=True,
        sharex=True,
    )
    if nrows == 1:
        ax = np.expand_dims(ax, axis=0)
    if ncols == 1:
        ax = np.expand_dims(ax, axis=1)
    return fig, ax
