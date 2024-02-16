#!/bin/bash
set -e -x

# load python env
source /project/g110/icon/probtest/conda/miniconda/bin/activate
conda activate probtest

# generate testdata
hash=$(cat ${0%/*/*/*}/test_reference)
ICON_DATA=/project/g110/probtest_testdata/$hash/icon_data PROBTEST_DATA=./probtest_data ./scripts/generate_icon_testdata.sh
