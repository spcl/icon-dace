# coupling_yaml_generate.sh - Generate YAML file for HD model in a coupled system with YAC
# 
# Copyright (C) 2022, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Ha Ho-Hagemann (Hereon, Germany)
# Contact: ha.hagemann@hereon.de
#_________________________________________
#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#
created_date="15-Jul-2022 09:50"
#
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
WRKDIR=$(pwd)

runstate=$1

if [ -e ./job_settings ] ; then 
  source ./job_settings
else
  echo 'Script job_settings does not exist in present directory scr!'
  exit
fi

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#
# (1) Fill model list
# -------------------
#
export EXPNAME=${EXP}

modelname[0]=${MODEL1[0]}
modelcomponent[0]=${MODEL1[1]}
modelgrid[0]=${MODEL1[2]}
model_id[0]=${MODEL1[3]}

modelname[1]=${MODEL2[0]}
modelcomponent[1]=${MODEL2[1]}
modelgrid[1]=${MODEL2[2]}
model_id[1]=${MODEL2[3]}

modelname[2]=${MODEL3[0]}
modelcomponent[2]=${MODEL3[1]}
modelgrid[2]=${MODEL3[2]}
model_id[2]=${MODEL3[3]}

tran_id[0]=${exchanged_var1[0]} # "1"
tran_id[1]=${exchanged_var2[0]} # "2"
tran_id[2]=${exchanged_var3[0]} # "3"

var[0]=${exchanged_var1[1]} 	# "RDC2NEMO"
var[1]=${exchanged_var2[1]} 	# "RUNOFF_S"
var[2]=${exchanged_var3[1]} 	# "RUNOFF_G"

collsize[0]=${exchanged_var1[2]} # "1"
collsize[1]=${exchanged_var2[2]} # "2"
collsize[2]=${exchanged_var3[2]} # "3"

# (2) start and end date+time of experiment
# -------------------------------------
if [ ${runstate} == "restart" ];then
  YYYY=`cat ${WRKDIR}/log/${EXP}.year`
fi
YYYY_NEXT=$((${YYYY}+1))

start_date=${start_date:="${YYYY}-01-01T00:00:00"}
#    end_date=${end_date:="${YYYY_NEXT}-01-01T00:00:00"}
    end_date=${end_date:="${YYYY}-12-31T00:00:00"}

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# (3) prepare coupling.yaml
#------ source_to_target_map.txt -------------------
interpolation_method_file="source_to_target_map.txt"
cat > ${interpolation_method_file} << EOF
      - source_to_target_map:
          max_search_distance: 3.6
      - fixed:
          user_value: 0.0
EOF

#------ bernstein_bezier.txt -------------------
interpolation_method_file="bernstein_bezier.txt"
cat > ${interpolation_method_file} << EOF
      - bernstein_bezier
      - nnn:
          n: 4
      - fixed:
          user_value: -999.0
EOF

#------ conservative.txt -------------------
interpolation_method_file="conservative.txt"
cat > ${interpolation_method_file} << EOF
      - conservative:
          order: 1
EOF

#------ user_file.txt -------------------
interpolation_method_file="user_file.txt"
cat > ${interpolation_method_file} << EOF
      - user_file:
EOF
#
# *** 3.2: coupling.yaml file
#
# component names in coupling.yaml must (!) match with modelname[*]
#
coupling_yaml_file="coupling_${EXPNAME}.yaml"

cat > ${coupling_yaml_file} << EOF
created date: ${created_date}
start_date: ${start_date}
end_date: ${end_date}
timestep_unit: second
calendar: proleptic-gregorian
coupling:
# ocean -> HD
  - src_component: ${modelname[2]}
    src_grid: ${modelgrid[2]}
    tgt_component: ${modelname[0]}
    tgt_grid: ${modelgrid[0]}
    coupling_period: ${couplingTimeStep}
    time_reduction: accumulate
    field: ${var[0]}
    src_lag: ${riv_lag}
    tgt_lag: ${oce_lag}
    interpolation:
EOF
if [ ${runstate} == "restart" ];then
  cat user_file.txt >> ${coupling_yaml_file}
cat >> ${coupling_yaml_file} << EOF
          filename: ${var[0]}.nc
EOF
else
  cat source_to_target_map.txt >> ${coupling_yaml_file}
cat >> ${coupling_yaml_file} << EOF
    weight_file_name: ${var[0]}.nc
EOF
fi
cat >> ${coupling_yaml_file} << EOF
# land -> HD
  - src_component: ${modelname[2]}
    src_grid: ${modelgrid[2]}
    tgt_component: ${modelname[0]}
    tgt_grid: ${modelgrid[0]}
    coupling_period: ${couplingTimeStep}
    time_reduction: accumulate
    field: ${var[0]}
    src_lag: ${lnd_lag}
    tgt_lag: ${riv_lag}
    interpolation:
EOF
if [ ${runstate} == "restart" ];then
  cat user_file.txt >> ${coupling_yaml_file}
  cat >> ${coupling_yaml_file} << EOF
          filename: ${var[1]}.nc
EOF
else
  cat conservative.txt >> ${coupling_yaml_file}
cat >> ${coupling_yaml_file} << EOF
    weight_file_name: ${var[1]}.nc
EOF
fi
cat >> ${coupling_yaml_file} << EOF
# land -> HD
  - src_component: ${modelname[2]}
    src_grid: ${modelgrid[2]}
    tgt_component: ${modelname[0]}
    tgt_grid: ${modelgrid[0]}
    coupling_period: ${couplingTimeStep}
    time_reduction: accumulate
    field: ${var[0]}
    src_lag: ${lnd_lag}
    tgt_lag: ${riv_lag}
    interpolation:
EOF
if [ ${runstate} == "restart" ];then
  cat user_file.txt >> ${coupling_yaml_file}
  cat >> ${coupling_yaml_file} << EOF
          filename: ${var[2]}.nc
EOF
else
  cat conservative.txt >> ${coupling_yaml_file}
cat >> ${coupling_yaml_file} << EOF
    weight_file_name: ${var[2]}.nc
EOF
fi
