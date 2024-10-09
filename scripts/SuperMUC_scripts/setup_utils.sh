#!/bin/ksh

# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

#==============================================================================
#
# Evaluates the user parameters
#
#==============================================================================


warning()
{
  if [ $1 != 0 ] ; then
    echo "   WARNING : $2"
  fi
}

finish()
{
  if [ $1 != 0 ] ; then 
    echo "-------------------------------------------------------"
    echo $1
    echo "    ERROR : $2. STOP"
    exit $1
  fi
}
#-------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------
include_targets()
{
  in_targets=$1
  #target_list="./config/set_target_*"
  include_targets="failure"
  for target in $in_targets
  do
    target_list=`ls $base_folder/config/set_target_${target}*`
    for target_filename in $target_list
    do
      if [[ $target_filename != *~ ]] ; then
        . $target_filename
        include_targets="OK"
      fi
    done
  done
}
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
check_target()
{
  target=$1
  target_list=`ls $base_folder/config/set_target_${target}*`
  check_target="not_found"
  for target_filename in $target_list
  do
    if [[ $target_filename != *~ ]] ; then
      check_target="found"
    fi
  done
  if [[ $check_target == "found" ]] ; then
    check_target="OK"
  fi
}
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
eval_arg ()
{
  #IFS=" ="
  case $1 in
  "with-mpi" )     use_mpi=$2
  ;;
  "with-openmp" )  use_openmp=$2
  ;;
  "with-fortran" ) use_fortran=$2
  ;;
  "target" )       use_target=$2
  ;;
  "optimize-level" ) optimize_level=$2
  ;;
  "debug-mode" )   debug_mode=$2
  ;;
  "with-define" )  user_definitions="$2 $3 $4 $5 $6 $7 $8 $9"
  ;;
  "with-flags" )   user_flags="$2 $3 $4 $5 $6 $7 $8 $9"
  ;;
  "cpu-time" )     use_cpu_time=$2
  ;;
  "mpi-procs-pernode" ) use_mpi_procs_pernode=$2
  ;;
  "mpi-procs" )    use_mpi_procs=$2
  ;;
  "nproma" ) use_nproma=$2
  ;;
  "openmp-threads" ) use_openmp_threads=$2
  ;;
  "omp-stacksize" ) use_omp_stacksize=$2
  ;;
  "resources" ) use_resources=$2
  ;;
  "memory-model" ) use_memory_model=$2
  ;;
  "memory" ) use_memory=$2
  ;;
  "node-usage" ) use_node_usage=$2
  ;;
  "no-of-nodes" ) use_nodes=$2
  ;;
  "queue" ) use_queue=$2
  ;;
  "run" ) run=$2
  ;;
  "out-name" ) output_name=$2
#     echo "output_name=$output_name"
  ;;
  "in-script" ) input_name="$input_name $2"
  ;;
  "job-name" ) job_name="$2"
  ;;
  "in-folder" ) input_folder=$2
  ;;
  "out-folder" ) output_folder=$2
  ;;
  "with-shell" ) use_shell=$2
  ;;
  *) let "no_of_free_variables=$no_of_free_variables+1"
     free_variable[$no_of_free_variables]=$1
     free_variable_value[$no_of_free_variables]=$2 
#      echo "free_variable=${free_variable[$no_of_free_variables]} with value=${free_variable_value[$no_of_free_variables]}"
#      exit -1
  ;;
  esac
}
#-------------------------------------------------------------------------------------

