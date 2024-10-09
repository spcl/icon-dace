"""
This module contains unit tests for the `utils.py` module.
"""

from util.utils import process_member_num


def test_process_member_num_single_element():
    """
    Test case for a single element in the input list.
    It should generate a list from 1 to that number and convert each to a string.
    """
    input_data = [5]
    expected_output = [(1, "1"), (2, "2"), (3, "3"), (4, "4"), (5, "5")]
    assert process_member_num(input_data) == expected_output


def test_process_member_num_multiple_elements():
    """
    Test case for multiple elements in the input list.
    It should convert each number to a string.
    """
    input_data = [2, 3, 4]
    expected_output = [(2, "2"), (3, "3"), (4, "4")]
    assert process_member_num(input_data) == expected_output


def test_process_member_num_empty_list():
    """
    Test case for an empty input list.
    It should return an empty list.
    """
    input_data = []
    expected_output = []
    assert process_member_num(input_data) == expected_output


def test_process_member_num_single_element_zero():
    """
    Test case for a single element in the input list being zero.
    It should return an empty list.
    """
    input_data = [0]
    expected_output = []
    assert process_member_num(input_data) == expected_output
