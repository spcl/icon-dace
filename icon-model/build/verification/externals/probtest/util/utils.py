""""
This module provides utility functions for list and string operations, as well as
a function to generate seeds based on a member number for probtest.
"""

import re

import numpy as np


def unique_elements(inlist):
    unique = []
    for element in inlist:
        if element not in unique:
            unique.append(element)
    return unique


def first_idx_of(li, el):
    return li.index(el)


def last_idx_of(li, el):
    return len(li) - li[::-1].index(el) - 1


def numbers(s):
    "join all numbers from a string and return as an int (base 10)"
    return int("".join(re.findall(r"\d+", s)))


def get_seed_from_member_num(member_num, use_64_bits=True):
    # Ensure member_num is a valid positive integer
    if member_num <= 0 or member_num > 100:
        raise ValueError(
            "Invalid member number: only values between 1 and 100 are valid."
        )

    seeds = [
        4820203056152078343,
        -2050370216167377885,
        5288723146233497607,
        -4100974783045284818,
        2749886054240915477,
        -2698141363522547153,
        4211464910216594971,
        -3195019776896701904,
        947813232643208240,
        -3017506358830162380,
        6457971157868528684,
        -1953845546453850051,
        1921509193601203792,
        -5151669688969927610,
        7359535530781478017,
        -7391450849366071736,
        4462894312318012041,
        -1153799594491576246,
        6280275070608421016,
        -8484692618266580910,
        4631300305370527834,
        -3070534271375073186,
        9004873009304494755,
        -302121758898965408,
        2411332205156471480,
        -7450398516460191118,
        5399077892206697197,
        -200345814144568716,
        8109952707559379697,
        -8796238266403548550,
        3686067195567049466,
        -5224622705771655010,
        2050729923911099647,
        -6470169296301394303,
        3578726697261266567,
        -2955312497544182115,
        239410929803317998,
        -7685045599314604897,
        2592198181079681805,
        -5430321512478074716,
        2563722063726054062,
        -8106129836903808332,
        4356212043861455336,
        -4900825218918790988,
        2064819667132781234,
        -726244375936754502,
        7666486099514309817,
        -1999720451738726212,
        1511249452891823422,
        -1418548619670343491,
        7186232883256903869,
        -1352220800016635706,
        6539476031458837823,
        -5624549502336598326,
        7624451126352253709,
        -2042604993628170038,
        8112820910641906648,
        -1614233251291285295,
        7232445954078326991,
        -3846078357683954464,
        2354119075015407358,
        -2545232179950367008,
        3870498242096184092,
        -8315137910183253273,
        7682552012510876118,
        -2622863005352196890,
        2374653056844096748,
        -3589113757772089599,
        8407878844808879891,
        -3030008833914349273,
        6036827779073203513,
        -4080542692836830403,
        177902514776640316,
        -467067027994043070,
        7767918186497823696,
        -755117393590705851,
        6504560056028236748,
        -6737289754062327991,
        5360428069762684338,
        -7214175785528760502,
        444174484151156551,
        -4878173763282084013,
        7888957635420352936,
        -1192230971368046761,
        8389364988510469013,
        -9084871955485573255,
        7431885829203443597,
        -8113010284615350917,
        4711890985677493628,
        -5759250609871551630,
        1672287084338692485,
        -5484416819603667576,
        4766790831848565149,
        -5853819326041077352,
        189441956786792371,
        -6821719999876591705,
        7792553165469396381,
        -8011518733766621251,
        3037342685824194534,
        -4677366596323783742,
    ]

    seed = seeds[member_num - 1]
    if not use_64_bits:
        seed = np.uint32(seed & 0xFFFFFFFF)

    return seed


def process_member_num(member_num):
    """
    Processes the member numbers.
    If a single number is provided, it generates a list of numbers from 1 to
    that number.
    It then converts each number to a string.

    Args:
        member_num (list): A list containing either a single integer or multiple
                           integers.

    Returns:
        list: A list of strings representing the processed member numbers.
    """
    if len(member_num) == 1:
        member_num = list(range(1, member_num[0] + 1))
    return [(m_num, str(m_num)) for m_num in member_num]
