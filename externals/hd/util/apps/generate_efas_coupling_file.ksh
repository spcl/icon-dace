#!/bin/ksh
#   
# generate_efas_coupling.ksh - Generate EFAS coupling/remapping file o a selected NEMO ocean grid 
# 
# Copyright (C) 2024, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
SCR=/home/g/g260122/hdmodel/util/generate_hd_coupling_file.ksh 
#
IDIS=2    # Discharge grid: 1=HD5.1-euro5min, 2=EFAS
INEMO=1   # Nemo grid: 1 = GCOAST, 2=North Sea, 3=Nordic 
#
SRAD=100.
clexb=''
#
# Discharge grid
case $IDIS in
   1) tagdis=hd5_1  
      SBOUND=0.0833333
      DNRIV=/work/gg0302/g260122/HD/input/euro5min/hdpara_vs5_1_euro5min.nc  ;;
   2) tagdis=efas  
      SBOUND=0.05
      DNRIV=/work/gg0302/g260122/HD/input/efas/mask_catchment_full_efas.nc  ;;
esac      
#
# NEMO grid
case $INEMO in 
   1 ) EXP=${tagdis}_to_nemo-gcoast
       SRAD=150.
       clexb='-x 1'
       DNO=/work/gg0302/g260122/HD/input/nemo/NEMO_grid_GCOAST.nc ;;
   2 ) EXP=${tagdis}_to_nemo-ns
       SRAD=67.    # Works for HD5 and EFAS
       DNO=/work/gg0302/g260122/HD/input/nemo/nemo_nordsee_bsh_domain_cfg.nc ;;
   3 ) EXP=${tagdis}_to_nemo-nordic
       SRAD=60.    # Works for HD5 and EFAS
       DNO=/work/gg0302/g260122/HD/input/nemo/slm_nemo-nordic.nc ;;
   * ) echo "Grid no. $INEMO does not exist --> TERMINATE" ; exit  ;;
esac
#
ulimit -s unlimited
cd /scratch/g/g260122/tmp
$SCR -s $DNRIV -o $DNO -r $SRAD -b $SBOUND $clexb -i $EXP

