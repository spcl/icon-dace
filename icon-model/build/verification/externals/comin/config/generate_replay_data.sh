#!/usr/bin/bash

set -ex

THIS_DIR=`pwd -P`
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

ICON_DIR=${ICON_DIR=`pwd -P`/icon}
cd $ICON_DIR
ICON_COMMIT_SHA=`git rev-parse HEAD`

# preparations for make_runscripts
cp $ICON_DIR/run/*.iconrun ${SCRIPTPATH}/replay_data_run/
# note: we restrict the job to 4 processes
echo "use_mpi_procs_pernode='4'" >> $ICON_DIR/run/set-up.info
./make_runscripts --all -r `realpath -s --relative-to=${ICON_DIR} ${SCRIPTPATH}/replay_data_run`

cd ${SCRIPTPATH}/replay_data_run
for runscript in *.run; do
    ./${runscript}
done

cd ${ICON_DIR}/experiments
rm -rf ${THIS_DIR}/predef_replay_data.cmake
for exp in *; do
    pushd ${exp}
    tar -cvzf ${exp}.tar.gz [0-9].nc vars_*.nc
    URL="${CI_API_V4_URL}/projects/${CI_PROJECT_ID}/packages/generic/replay_data/1/${ICON_COMMIT_SHA}_${exp}.tar.gz"
    curl --fail \
         --header "JOB-TOKEN: $CI_JOB_TOKEN" \
         --upload-file ${exp}.tar.gz \
         ${URL}
    cat >> ${THIS_DIR}/predef_replay_data.cmake <<EOF
comin_add_replay_data(NAME ${exp}
                      URL ${URL}
                      MD5HASH `md5sum ${exp}.tar.gz | awk '{ print $1 }'`)
EOF
    popd
done
