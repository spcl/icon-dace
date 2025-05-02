"""
CLI for member selection

This script provides a command line interface for selecting members and testing
tolerance factors for statistical data.
It can find members and a corresponding tolerance factor that validate for all
stats files or test the tolerance of a given selection.
"""

import logging
import os
import sys
from collections import Counter
from datetime import datetime

import click
import numpy as np

from engine.tolerance import tolerance
from util.click_util import cli_help
from util.dataframe_ops import test_stats_file_with_tolerances
from util.log_handler import logger


# finds members and a corresponding tolerance factor validating for all stats files.
def find_members_and_factor_validating_for_all_stats_files(
    random_tolerance_file_name,
    stats_file_name,
    member_type,
    min_member_num,
    max_member_num,
    total_member_num,
    min_factor,
    max_factor,
    iterations,
):

    members = list(range(1, total_member_num + 1))

    # Iteratively change likelihood of members being selected
    # to come to a solution faster
    weights = np.ones(total_member_num) / total_member_num

    # Use fixed seed to make result reproducible
    np.random.seed(total_member_num)

    for f in range(
        int(min_factor), int(max_factor) + 1, 5
    ):  # Try with bigger factor if max_member_num is not enough
        logger.info("Set factor to %s", f)
        for mem_num in range(min_member_num, max_member_num + 1):
            logger.info("Try with %s members.", mem_num)
            max_passed = 1
            variables = []
            for iteration in range(iterations):
                random_members = np.random.choice(
                    members, size=mem_num, replace=False, p=weights
                )
                logger.info(
                    "Test %s with %s randomly selected members and factor %s.",
                    iteration + 1,
                    mem_num,
                    f,
                )
                # Create tolerances from random members
                context = click.Context(tolerance)
                context.invoke(
                    tolerance,
                    stats_file_name=stats_file_name,
                    tolerance_file_name=random_tolerance_file_name,
                    member_num=random_members,
                    member_type=member_type,
                )
                # Test selection (exclude random selection)
                validation_members = [
                    item for item in members if item not in random_members
                ]
                passed, new_vars = test_selection(
                    stats_file_name,
                    random_tolerance_file_name,
                    validation_members,
                    member_type,
                    f,
                )

                validation_members_np = np.array(validation_members)
                indices = [i for i, value in enumerate(passed) if value == 0]
                failed = validation_members_np[indices]
                # Increase weights for members which failed
                weights[failed - 1] += 1 / total_member_num
                # weights needs to sum up to 1 for np.random.choice
                weights = weights / sum(weights)

                variables.extend(new_vars)

                if (mem_num == max_member_num) and (f == max_factor):
                    duplicates = dict(Counter(variables).items())
                    sorted_duplicates = dict(
                        sorted(duplicates.items(), key=lambda x: x[1], reverse=True)
                    )
                # The following is to save computing time
                elif iteration < 32:
                    max_passed = max(max_passed, sum(passed))
                    # The more combs were tested
                    # the higher should the success rate be to continue
                    tested_stats = len(validation_members)
                    if max_passed < (iteration + 1) * 0.03 * tested_stats:
                        break

                if sum(passed) == len(validation_members):
                    return random_members, f
        # If factore needs to be increased, test only with max_members
        min_member_num = max_member_num

    max_count = max(sorted_duplicates.values())
    most_common_vars = [
        item for item, count in sorted_duplicates.items() if count == max_count
    ]
    logger.error(
        "ERROR: Could not find %s random members, which pass for all stat files. "
        + "The most sensitive variable(s) is/are %s, which failed for %s out of %s "
        + "random selections. Consider removing this/those variable(s) from the "
        + "experiment and run again.",
        max_member_num,
        most_common_vars,
        max_count,
        iteration + 1,
    )
    sys.exit(1)


# Tests how may stats files pass the tolerance test for the selected members
# Returns the number of passed stats files and the variables which failed
def test_selection(
    stats_file_name, tolerance_file_name, total_member_num, member_type, factor
):
    if isinstance(total_member_num, int):
        members = list(range(1, total_member_num + 1))
    else:
        members = total_member_num
        total_member_num = len(members)

    passed = [0] * total_member_num

    # Change level to not get whole output from test_stats_file_with_tolerances
    original_level = logging.getLogger().level
    logging.getLogger().setLevel(logging.ERROR)

    variables = []
    i = 0
    for m_num in members:
        m_id = str(m_num) if not member_type else member_type + "_" + str(m_num)

        out, err, _ = test_stats_file_with_tolerances(
            tolerance_file_name,
            stats_file_name.format(member_id="ref"),
            stats_file_name.format(member_id=m_id),
            factor,
        )

        if not out:
            var = set(index[1] for index in err[0].index)
            var = list(var)
            variables.extend(var)
        else:
            passed[i] = 1
        i = i + 1

    variables = list(set(variables))

    # Reset logger level
    logging.getLogger().setLevel(original_level)
    logger.info(
        "The tolerance test passed for %s out of %s references.",
        sum(passed),
        total_member_num,
    )
    return passed, variables


@click.command()
@click.option(
    "--experiment-name",
    help=cli_help["experiment_name"],
)
@click.option(
    "--test-tolerance/--no-test-tolerance",
    is_flag=True,
    help=cli_help["test_tolerance"],
)
@click.option(
    "--stats-file-name",
    help=cli_help["stats_file_name"],
)
@click.option(
    "--selected-members-file-name",
    help=cli_help["selected_members_file_name"],
)
@click.option(
    "--tolerance-file-name",
    help=cli_help["tolerance_file_name"],
)
@click.option(
    "--member-type",
    type=str,
    default="",
    help=cli_help["member_type"],
)
@click.option(
    "--min-member-num",
    type=int,
    default=5,
    help=cli_help["min_member_num"],
)
@click.option(
    "--max-member-num",
    type=int,
    default=15,
    help=cli_help["max_member_num"],
)
@click.option(
    "--total-member-num",
    type=int,
    default=50,
    help=cli_help["total_member_num"],
)
@click.option(
    "--factor",
    type=float,
    default=5.0,
    help=cli_help["factor"],
)
@click.option(
    "--min-factor",
    type=float,
    default=5.0,
    help=cli_help["min_factor"],
)
@click.option(
    "--max-factor",
    type=float,
    default=50.0,
    help=cli_help["max_factor"],
)
@click.option(
    "--iterations",
    type=int,
    default=50,
    help=cli_help["iterations"],
)
# Selects members and writes them to a file together with the tolerance factor
def select_members(
    experiment_name,
    test_tolerance,
    stats_file_name,
    selected_members_file_name,
    tolerance_file_name,
    member_type,
    min_member_num,
    max_member_num,
    total_member_num,
    factor,
    min_factor,
    max_factor,
    iterations,
):  # pylint: disable=unused-argument

    if min_member_num > max_member_num:
        logger.error(
            "ERROR: min_member_num must be equal or smaller than max_member_num"
        )
        sys.exit(1)

    if min_factor > max_factor:
        logger.error("ERROR: min_factor must be equal or smaller than max_factor")
        sys.exit(1)

    if max_member_num >= total_member_num:
        logger.error("ERROR: max_member_num must be smaller than total_member_num")
        sys.exit(1)

    if test_tolerance:
        # Test selection
        test_selection(
            stats_file_name, tolerance_file_name, total_member_num, member_type, factor
        )
    else:
        random_tolerance_file_name = f"random_tolerance_{experiment_name}.csv"
        start_time = datetime.now()
        selection, factor = find_members_and_factor_validating_for_all_stats_files(
            random_tolerance_file_name,
            stats_file_name,
            member_type,
            min_member_num,
            max_member_num,
            total_member_num,
            min_factor,
            max_factor,
            iterations,
        )
        end_time = datetime.now()
        elapsed_time = end_time - start_time
        logger.info("The member selection took %ss.", elapsed_time.total_seconds())

        # Write selection into a file
        selection = ",".join(map(str, selection))
        logger.info(
            "Writing selected members %s with tolerance factor %s to file %s",
            selection,
            int(factor),
            selected_members_file_name,
        )
        with open(selected_members_file_name, "w", encoding="utf-8") as file:
            file.write(selection + "\n")
            file.write("export FACTOR=" + str(int(factor)))

        # The last created file was successful
        logger.info("Writing tolerance file to %s", tolerance_file_name)
        os.rename(random_tolerance_file_name, tolerance_file_name)
