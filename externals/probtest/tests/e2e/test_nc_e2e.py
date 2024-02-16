import os
import unittest
from pathlib import Path

import numpy as np
import xarray as xr


def load_netcdf(rel_path):
    nc_ref_path = Path(os.environ["PROBTEST_REF_DATA"]) / rel_path
    nc_cur_path = Path(os.environ["PROBTEST_CUR_DATA"]) / rel_path

    data_ref = xr.load_dataset(nc_ref_path)
    data_cur = xr.load_dataset(nc_cur_path)

    return data_ref, data_cur


def check_netcdf(data_ref, data_cur):
    diff_keys = set(data_ref.keys()) - set(data_cur.keys())
    same_keys = set(data_ref.keys()).intersection(set(data_cur.keys()))

    err = []
    for key in same_keys:
        diff = np.fabs(data_ref[key] - data_cur[key])
        if np.sum(diff) > 0:
            err.append(key)

    return list(diff_keys), err


class TestNcE2E(unittest.TestCase):
    def setUp(self) -> None:
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_perturb_e2e(self):
        exp = os.environ["PROBTEST_TEST_EXPERIMENT"]

        # Get a list of directories that include the experiment name
        perturb_path = Path(os.getcwd(), os.environ["PROBTEST_REF_DATA"], "perturb")
        exp_dirs = [d for d in os.listdir(perturb_path) if exp in d]
        self.assertNotEqual(
            len(exp_dirs),
            0,
            msg="No experiment folders found in {}".format(perturb_path),
        )

        for exp_dir in exp_dirs:
            data_ref, data_cur = load_netcdf(
                "perturb/{}/initial_condition.nc".format(exp_dir)
            )

            diff_keys, err = check_netcdf(data_ref, data_cur)
            self.assertEqual(
                err, [], msg="The following variables contain errors:\n{}".format(err)
            )
            self.assertEqual(
                diff_keys,
                [],
                msg="The following variables are not contained in both files:\n"
                + "{}".format(err),
            )
