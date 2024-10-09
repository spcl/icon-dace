"""
This module contains test cases to validate the functionality of member
selection and tolerance testing CLI commands.
"""

import os
import re

from tests.helpers import run_select_members_cli, run_tolerance_cli


def test_select_members(stats_file_set):
    run_select_members_cli(
        stats_file_set["stats"],
        stats_file_set["members"],
        stats_file_set["tol"],
    )

    assert os.path.isfile(
        stats_file_set["members"]
    ), "File 'selected_members.csv' was not created"

    with open(stats_file_set["members"], "r", encoding="utf-8") as file:
        content = file.read().strip()
    expected_content = "50,21,40,39,16\nexport FACTOR=5"
    assert content == expected_content, "The member selection failed"


def test_select_members_increase_factor(stats_file_set):
    run_select_members_cli(
        stats_file_set["stats"],
        stats_file_set["members"],
        stats_file_set["tol"],
        max_member_num=5,
        iterations=1,
    )

    assert os.path.isfile(
        stats_file_set["members"]
    ), "File 'selected_members.csv' was not created"

    with open(stats_file_set["members"], "r", encoding="utf-8") as file:
        content = file.read().strip()
    expected_content = "50,21,40,39,16\nexport FACTOR=10"
    assert (
        content == expected_content
    ), "Increasing the factor within the member selection failed"


def test_select_members_failure(stats_file_set, caplog):
    log = run_select_members_cli(
        stats_file_set["stats"],
        stats_file_set["members"],
        stats_file_set["tol"],
        max_member_num=5,
        iterations=1,
        max_factor=5,
        log=caplog,
    )
    assert "ERROR" in log, "The member selection did not fail, although it should have."


def test_test_tolerance(stats_file_set, caplog):
    run_tolerance_cli(stats_file_set["stats"], stats_file_set["tol"], member_num=5)

    log = run_select_members_cli(
        stats_file_set["stats"],
        stats_file_set["members"],
        stats_file_set["tol"],
        test_tolerance=True,
        log=caplog,
    )

    match = re.search(r"passed for (\d+) out of 50", log)
    passed_count = match.group(1) if match else "0"
    error_message = (
        f"The test-tolerance output is incorrect. It should pass for "
        f"49 out of 50 but it passed for {passed_count} out of 50."
    )

    assert "49 out of 50" in log, error_message
