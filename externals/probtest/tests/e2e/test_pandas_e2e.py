import os
import unittest
from pathlib import Path

import pandas as pd

from util.tree import TimingTree

pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.float_format", "{:,.2e}".format)


def load_pandas(rel_path, index_col=[0, 1], header=[0, 1]):
    input_file_ref = Path(os.environ["PROBTEST_REF_DATA"]) / rel_path
    input_file_cur = Path(os.environ["PROBTEST_CUR_DATA"]) / rel_path
    df_ref = pd.read_csv(input_file_ref, index_col=index_col, header=header)
    df_cur = pd.read_csv(input_file_cur, index_col=index_col, header=header)

    return df_ref, df_cur


def pandas_error(df_ref, df_cur):
    diff = (df_ref - df_cur).abs()
    err_mask = (diff > 0).any(axis=1)

    return diff[err_mask]


class TestPandasE2E(unittest.TestCase):
    def setUp(self) -> None:
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_stats_e2e(self):
        exp_name = os.environ["PROBTEST_TEST_EXPERIMENT"]
        rel_path = "reference/{}_ref.csv".format(exp_name)
        df_ref, df_cur = load_pandas(rel_path, index_col=[0, 1, 2])

        err = pandas_error(df_ref, df_cur)

        self.assertEqual(
            len(err.values), 0, "Stats datasets are not equal!\n{}".format(err)
        )

        return

    def test_tolerance_e2e(self):
        exp_name = os.environ["PROBTEST_TEST_EXPERIMENT"]
        rel_path = "tolerance/{}.csv".format(exp_name)
        df_ref, df_cur = load_pandas(rel_path)

        err = pandas_error(df_ref, df_cur)

        self.assertEqual(
            len(err.values), 0, "Tolerance datasets are not equal!\n{}".format(err)
        )

        return

    def test_cdo_e2e(self):
        exp_name = os.environ["PROBTEST_TEST_EXPERIMENT"]
        rel_path = "cdo_table/{}.csv".format(exp_name)
        df_ref, df_cur = load_pandas(rel_path)

        err = pandas_error(df_ref, df_cur)

        self.assertEqual(
            len(err.values), 0, "CDO table datasets are not equal!\n{}".format(err)
        )

        return

    def test_performance_data_e2e(self):
        exp_name = os.environ["PROBTEST_TEST_EXPERIMENT"]

        df_ref = TimingTree.from_json(
            Path(os.environ["PROBTEST_REF_DATA"]) / "performance/{}".format(exp_name)
        )
        df_cur = TimingTree.from_json(
            Path(os.environ["PROBTEST_CUR_DATA"]) / "performance/{}".format(exp_name)
        )

        for i in range(len(df_ref.data)):
            err = pandas_error(df_ref.data[i], df_cur.data[i])

            self.assertEqual(
                len(err.values),
                0,
                "Performance datasets are not equal for table {}!\n{}".format(err, i),
            )
