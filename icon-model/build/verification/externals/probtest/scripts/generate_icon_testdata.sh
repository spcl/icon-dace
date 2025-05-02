#!/bin/bash

PROBTEST_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cd .. &> /dev/null && pwd )"

if [[ -z "$ICON_DATA" ]]; then
    echo "Path to the collected ICON run and output data must be set as 'ICON_DATA'"
    exit 1
fi

if [[ -z "$PROBTEST_DATA" ]]; then
    echo "Path where the output of probtest should be stored must be set as 'PROBTEST_DATA'"
    exit 1
fi

for exp in mch_opr_r04b07
do
    file_id=(
        --file-id NetCDF "*atm_3d_ml*.nc"
    )
    python $PROBTEST_DIR/probtest.py init \
        --codebase-install "$ICON_DATA" \
        --experiment-name "$exp" \
        "${file_id[@]}" \
        --reference "$PROBTEST_DATA" \
        --template-name "$PROBTEST_DIR/templates/testdata.jinja" \
        --member-num 2,5 \
        --timing-current "$ICON_DATA/probtest_testdata" \
        --timing-reference "$ICON_DATA/performance_reference/probtest_testdata" \
        || exit 1

    python "$PROBTEST_DIR/probtest.py" run-ensemble || exit 1
    python "$PROBTEST_DIR/probtest.py" stats || exit 1
    python "$PROBTEST_DIR/probtest.py" tolerance || exit 1
    python "$PROBTEST_DIR/probtest.py" cdo-table || exit 1
    python "$PROBTEST_DIR/probtest.py" performance || exit 1
    python "$PROBTEST_DIR/probtest.py" check-plot || exit 1
    python "$PROBTEST_DIR/probtest.py" perturb || exit 1

    # Use the `true` tool to test the non-dry code path. `true` does nothing and returns 0. See `man true`.
    python "$PROBTEST_DIR/probtest.py" run-ensemble --no-dry --submit-command=true || exit 1
done

# the following make use of long-running performance database for mch_ch_lowres
# python "$PROBTEST_DIR/probtest.py" performance-plot || exit 1
# python "$PROBTEST_DIR/probtest.py" performance-meta-data || exit 1
