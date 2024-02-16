import os
import unittest

import numpy as np
from netCDF4 import Dataset

from engine.stats import create_stats_dataframe

TIME_DIM_SIZE = 3
HOR_DIM_SIZE = 100
HEIGHT_DIM_SIZE = 5


def initialize_dummy_netcdf_file(name):
    data = Dataset(name, "w")

    data.createDimension("t", size=TIME_DIM_SIZE)
    data.createVariable("t", np.float64, dimensions="t")
    data.variables["t"][:] = np.arange(TIME_DIM_SIZE)

    data.createDimension("x", size=HOR_DIM_SIZE)
    data.createVariable("x", np.float64, dimensions="x")
    data.variables["x"][:] = np.arange(HOR_DIM_SIZE)

    data.createDimension("z", size=HEIGHT_DIM_SIZE)
    data.createVariable("z", np.float64, dimensions="z")
    data.variables["z"][:] = np.arange(HEIGHT_DIM_SIZE)

    return data


class TestStatsNetcdf(unittest.TestCase):
    nc_file_name = "test_stats.nc"
    nc_file_glob = "test_s*.nc"
    stats_file_names = "test_stats.csv"

    def setUp(self):
        data = initialize_dummy_netcdf_file(self.nc_file_name)

        data.createVariable("v1", np.float64, dimensions=("t", "z", "x"))
        data.variables["v1"][:] = np.ones(
            (TIME_DIM_SIZE, HEIGHT_DIM_SIZE, HOR_DIM_SIZE)
        )
        data.variables["v1"][:, :, 0] = 0
        data.variables["v1"][:, :, -1] = 2

        data.createVariable("v2", np.float64, dimensions=("t", "x"), fill_value=42)
        data.variables["v2"][:] = np.ones((TIME_DIM_SIZE, HOR_DIM_SIZE)) * 2
        data.variables["v2"][:, 0] = 1
        data.variables["v2"][:, 1] = 42  # shall be ignored in max-statistic
        data.variables["v2"][:, -1] = 3

        data.createVariable("v3", np.float64, dimensions=("t", "x"))
        data.variables["v3"][:] = np.ones((TIME_DIM_SIZE, HOR_DIM_SIZE)) * 3
        data.variables["v3"][:, 0] = 2
        data.variables["v3"][:, -1] = 4

        data.close()

    def TearDown(self):
        os.remove(self.nc_file_name)
        os.remove(self.stats_file_names)

    def test_stats(self):
        file_specification = {
            "Test data": dict(
                format="netcdf",
                time_dim="t",
                horizontal_dims=["x"],
                fill_value_key="_FillValue",  # This should be the name for fill_value.
            ),
        }

        df = create_stats_dataframe(
            input_dir=".",
            file_id=[["Test data", self.nc_file_glob]],
            stats_file_name=self.stats_file_names,
            file_specification=file_specification,
        )

        # check that the mean/max/min are correct
        expected = np.array(
            [
                [1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0],
                [2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0],
                [3.0, 4.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 2.0],
            ]
        )

        self.assertTrue(
            np.array_equal(df.values, expected),
            "stats dataframe incorrect. Difference:\n{}".format(df.values == expected),
        )


class TestStatsCsv(unittest.TestCase):
    dat_file_name = "test_stats_csv.dat"
    stats_file_name = "test_stats_csv.csv"

    def setUp(self):
        lines = (
            "time v1  v2 v3 v4 v5",
            "10   1.4 15 16 17 18",
            "20   2.4 25 26 27 28",
            "30   3.4 35 36 37 38",
        )
        with open(self.dat_file_name, "w") as f:
            f.write("\n".join(lines))

    def TearDown(self):
        os.remove(self.dat_file_name)
        os.remove(self.stats_file_name)

    def test_stats(self):
        file_specification = {
            "Test data": dict(
                format="csv",
                parser_args=dict(
                    delimiter="\\s+",
                    header=0,
                    index_col=0,
                ),
            ),
        }

        df = create_stats_dataframe(
            input_dir=".",
            file_id=[["Test data", self.dat_file_name]],
            stats_file_name=self.stats_file_name,
            file_specification=file_specification,
        )

        # check that the mean/max/min are correct (i.e. the same as in CSV)
        expected = np.array(
            [
                [1.4, 1.4, 1.4, 2.4, 2.4, 2.4, 3.4, 3.4, 3.4],
                [15, 15, 15, 25, 25, 25, 35, 35, 35],
                [16, 16, 16, 26, 26, 26, 36, 36, 36],
                [17, 17, 17, 27, 27, 27, 37, 37, 37],
                [18, 18, 18, 28, 28, 28, 38, 38, 38],
            ],
        )

        self.assertTrue(
            np.array_equal(df.values, expected),
            "stats dataframe incorrect. Difference:\n{}".format(df.values == expected),
        )


if __name__ == "__main__":
    unittest.main()
