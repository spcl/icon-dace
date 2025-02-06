"""
This module contains unittests for verifying the behavior of the `perturb_array`
function and its application to numpy arrays and NetCDF4 files.
It ensures the correctness of array perturbations and validates the precision
and amplitude of the perturbations.
"""

import numpy as np
import pytest
from matplotlib import pyplot as plt
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from engine.perturb import perturb_array

atype = np.float32
AMPLITUDE = atype(1e-14)
ARRAY_DIM = 100


@pytest.fixture(name="create_nc_files")
def fixture_create_nc_files(tmp_dir):
    data = [None] * 2
    for i, dt in enumerate([np.float32, np.float64]):
        data[i] = Dataset(f"{tmp_dir}/dummy{i}.nc", "w")
        data[i].createDimension("x", size=ARRAY_DIM)
        data[i].createDimension("y", size=ARRAY_DIM)
        data[i].createVariable("z", dt, dimensions=("x", "y"))
        data[i].variables["z"][:] = np.ones((ARRAY_DIM, ARRAY_DIM))
    yield data
    for d in data:
        d.close()


def test_perturb_array():
    # create two arrays, perturb one.
    x1 = np.ones((ARRAY_DIM, ARRAY_DIM), dtype=atype)
    x2 = np.ones((ARRAY_DIM, ARRAY_DIM), dtype=atype)
    x_perturbed = perturb_array(x2, 10, AMPLITUDE)

    # compute some stats and do assertions
    diff1 = np.abs(x1 - x_perturbed) / x1
    diff2 = np.abs(x2 - x_perturbed) / x2
    mean_diff1 = np.mean(diff1)
    mean_diff2 = np.mean(diff2)
    assert np.max(diff2) < AMPLITUDE, "perturbation is larger than amplitude!"
    assert np.isclose(
        mean_diff1, AMPLITUDE * 0.5, atol=AMPLITUDE * 1e-2
    ), "perturbation is most likely too small!"
    assert np.isclose(
        mean_diff2, AMPLITUDE * 0.5, atol=AMPLITUDE * 1e-2
    ), "perturbation did not return a copy of input!"


def test_perturb_nc(tmp_dir, create_nc_files):

    data = create_nc_files
    zd = data[1].variables["z"][:]
    zd_perturb = perturb_array(zd, 10, AMPLITUDE)
    data[1].variables["z"][:] = zd_perturb

    # reopen the files to make sure we get the values from disk
    for i in range(2):
        data[i] = Dataset(f"{tmp_dir}/dummy{i}.nc", "r")
    zf = data[0].variables["z"][:]
    zd_perturb = data[1].variables["z"][:]

    # compute some stats and do assertions
    diff = zf - zd_perturb
    mean_diff = np.mean(np.abs(diff))
    assert zf.dtype == "float32", "zf is not float32!"
    assert zd.dtype == "float64", "zf is not float64!"
    assert np.max(diff) < AMPLITUDE, "perturbation is larger than amplitude!"
    assert np.isclose(
        mean_diff, AMPLITUDE * 0.5, atol=AMPLITUDE * 1e-2
    ), "perturbation is most likely too small!"

    # create a plot of the difference to be read out manually if you are scared
    fig, ax = plt.subplots(1, 1)
    xx, yy = np.meshgrid(np.linspace(0, 1, ARRAY_DIM), np.linspace(0, 1, ARRAY_DIM))
    cs = ax.contourf(xx, yy, diff)
    fig.colorbar(cs, ax=ax)
    fig.savefig(f"{tmp_dir}/diff_figure.pdf")
