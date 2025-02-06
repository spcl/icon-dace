"""
This module contains unit tests for the `visualize_util.py` module.
"""

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from util.visualize_util import create_figure


def test_create_figure_single_row():
    ncols = 3
    nrows = 1
    fig_width = 10
    fig_height = 5

    fig, ax = create_figure(ncols, nrows, fig_width, fig_height)

    # Check that the output is a figure and a numpy array
    assert isinstance(fig, Figure)
    assert isinstance(ax, np.ndarray)

    # Check the shape of the axes array
    assert ax.shape == (1, ncols)

    # Check that each element in the array is an Axes object
    for a in ax.flatten():
        assert isinstance(a, Axes)


def test_create_figure_single_column():
    ncols = 1
    nrows = 3
    fig_width = 5
    fig_height = 10

    fig, ax = create_figure(ncols, nrows, fig_width, fig_height)

    # Check that the output is a figure and a numpy array
    assert isinstance(fig, Figure)
    assert isinstance(ax, np.ndarray)

    # Check the shape of the axes array
    assert ax.shape == (nrows, 1)

    # Check that each element in the array is an Axes object
    for a in ax.flatten():
        assert isinstance(a, Axes)


def test_create_figure_multiple_rows_and_columns():
    ncols = 3
    nrows = 2
    fig_width = 15
    fig_height = 10

    fig, ax = create_figure(ncols, nrows, fig_width, fig_height)

    # Check that the output is a figure and a numpy array
    assert isinstance(fig, Figure)
    assert isinstance(ax, np.ndarray)

    # Check the shape of the axes array
    assert ax.shape == (nrows, ncols)

    # Check that each element in the array is an Axes object
    for a in ax.flatten():
        assert isinstance(a, Axes)


def test_create_figure_single_cell():
    ncols = 1
    nrows = 1
    fig_width = 5
    fig_height = 5

    fig, ax = create_figure(ncols, nrows, fig_width, fig_height)

    # Check that the output is a figure and a numpy array
    assert isinstance(fig, Figure)
    assert isinstance(ax, np.ndarray)

    # Check the shape of the axes array
    assert ax.shape == (1, 1)

    # Check that the element is an Axes object
    assert isinstance(ax[0, 0], Axes)
