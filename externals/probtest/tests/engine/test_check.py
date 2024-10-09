"""
This module contains unittests for verifying the behavior of data frame operations
related to the PROBTEST suite. It ensures the correctness of relative difference
calculations and the checking of variable values against specified tolerances.
"""

import unittest

import numpy as np
import pandas as pd

from util.constants import CHECK_THRESHOLD
from util.dataframe_ops import (
    check_intersection,
    check_variable,
    compute_rel_diff_dataframe,
)


class TestCheck(unittest.TestCase):
    """
    Unit tests for the functionality of tolerance checking in dataframes.

    This class uses the `unittest` framework to validate the correctness of
    functions involved in computing and checking tolerances between dataframes,
    especially in the context of relative differences between expected and
    observed data.
    """

    def setUp(self):
        index = pd.MultiIndex.from_arrays(
            [
                ["NetCDF:*atm_3d*.nc"] * 4,
                ["var_1"] * 3 + ["var_2"],
                list(range(3)) + [0],
            ],
            names=["file_ID", "variable", "height"],
        )
        columns = pd.MultiIndex.from_product(
            [range(4), ["max", "mean", "min"]], names=["time", "statistic"]
        )

        array1 = np.linspace(0.9, 1.9, 4 * 12).reshape(4, 12)
        array1[:2] *= -1  # make some test data negative
        self.df1 = pd.DataFrame(array1, index=index, columns=columns)
        array2 = np.linspace(1.1, 2.1, 4 * 12).transpose().reshape(4, 12)
        array2[:2] *= -1  # make some test data negative
        self.df2 = pd.DataFrame(array2, index=index, columns=columns)
        # Relative differences |df1-df2|/((1+|df1|) are between 0.069 and 0.105.

        self.tol_large = pd.DataFrame(
            np.ones((2, 12)) * 0.21, index=["var_1", "var_2"], columns=columns
        )
        self.tol_small = pd.DataFrame(
            np.ones((2, 12)) * 0.06, index=["var_1", "var_2"], columns=columns
        )

    def check(self, df1, df2):
        # compute relative difference
        diff_df = compute_rel_diff_dataframe(df1, df2)
        # take maximum over height
        diff_df = diff_df.groupby(["variable"]).max()

        out1, err1, _ = check_variable(diff_df, self.tol_large)
        out2, err2, _ = check_variable(diff_df, self.tol_small)
        self.assertTrue(
            out1,
            "Check with large tolerances did not validate! "
            + f"Here is the DataFrame:\n{err1}",
        )
        self.assertFalse(
            out2,
            "Check with small tolerances did validate! "
            + f"Here is the DataFrame:\n{err2}",
        )

    def test_check(self):
        self.check(self.df1, self.df2)

    def test_check_one_zero(self):
        """Probtest should not pass if any of the values is 0
        and the other is much larger"""
        df1 = self.df1.copy()
        df1.loc[("NetCDF:*atm_3d*.nc", "var_1", 2), (0, "max")] = 0
        df2 = self.df2.copy()
        diff_df = compute_rel_diff_dataframe(df1, df2)
        diff_df = diff_df.groupby(["variable"]).max()

        out, err, _ = check_variable(diff_df, self.tol_large)

        self.assertFalse(
            out,
            "Check with 0-value reference validated! "
            + f"Here is the DataFrame:\n{err}",
        )

        df2.loc[("NetCDF:*atm_3d*.nc", "var_1", 2), (0, "max")] = CHECK_THRESHOLD / 2
        # now, both data are comparable again
        self.check(df1, df2)

    def test_check_smalls(self):
        """both values are close to 0 and should be accepted even though
        their relative difference is large.

        Close to 0 means < util.constants.CHECK_THRESHOLD"""
        df1 = self.df1.copy()
        df1.loc["var_1", (2, "min")] = CHECK_THRESHOLD * 1e-5
        df2 = self.df2.copy()
        df2.loc["var_1", (2, "min")] = CHECK_THRESHOLD / -2
        self.check(df1, df2)

    def test_no_intersection(self):
        """Probtest should fail if the variables in the
        reference and test case have no intersection"""
        df1 = self.df1.copy()
        df1 = df1.rename(index={"var_1": "var_3", "var_2": "var_4"})
        skip_test, _, _ = check_intersection(df1, self.df2)

        self.assertNotEqual(
            skip_test,
            0,
            "No intersection of variables in reference "
            + "and test case but test didn't fail",
        )

    def test_missing_variables(self):
        """Probtest should through a warning if some variables
        are not in the reference and test case"""
        df1 = self.df1.copy()
        df1 = df1.drop("var_1", level="variable")

        expected_warning_msg = (
            "WARNING: The following variables are in the "
            "test case but not in the reference case and therefore not tested: var_1"
        )
        with self.assertWarnsRegex(UserWarning, expected_warning_msg):
            check_intersection(df1, self.df2)

        expected_warning_msg = (
            "WARNING: The following variables are in the "
            "reference case but not in the test case and therefore not tested: var_1"
        )
        with self.assertWarnsRegex(UserWarning, expected_warning_msg):
            check_intersection(self.df2, df1)


class TestCheckSwapped(TestCheck):
    """Test that all Checks are symmetrical"""

    def check(self, df1, df2):
        super().check(df2, df1)  # pylint: disable=arguments-out-of-order


if __name__ == "__main__":
    unittest.main()
