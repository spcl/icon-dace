import os
import unittest

from util.tree import TimingTree


class TestTreeE2E(unittest.TestCase):
    def setUp(self) -> None:
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_performance_tree_e2e(self):
        exp_name = os.environ["PROBTEST_TEST_EXPERIMENT"]
        ref_tree = TimingTree.from_json(
            "{}/performance/{}".format(os.environ["PROBTEST_REF_DATA"], exp_name)
        )
        cur_tree = TimingTree.from_json(
            "{}/performance/{}".format(os.environ["PROBTEST_CUR_DATA"], exp_name)
        )

        for i in range(len(ref_tree.root)):
            ref_nodes = set(ref_tree.root[i].to_ancestry_name_list())
            cur_nodes = set(cur_tree.root[i].to_ancestry_name_list())

            diff_nodes = list(ref_nodes - cur_nodes)

            self.assertEqual(
                diff_nodes,
                [],
                msg="The test tree does not contain the following nodes:\n{}".format(
                    diff_nodes
                ),
            )
