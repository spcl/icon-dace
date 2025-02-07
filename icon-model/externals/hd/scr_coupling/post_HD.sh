# post_HD.sh - Script to store outputs and restarts of the HD model
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Author: Ha Ho-Hagemann (Hereon, Germany)
# Contact: ha.hagemann@hereon.de
#_________________________________________
#
#################################################
# store output and restarts
#################################################
WRKDIR=$(pwd)

runstate=$1

if [ -e ./job_settings ] ; then 
  source ./job_settings
else
  echo 'Script job_settings does not exist in present directory scr!'
  exit
fi

if [ ${runstate} == "restart" ];then
  YYYY=`cat ${WRKDIR}/log/${EXP}.year`
fi

echo "Post-processing for YYYY" ${YYYY}

HDOUT=${HDDIR}/${EXP}/out          # HD Output dir

cd $HDDIR/$EXP

##typeset -Z2 MM_NEXT
if (( ${IWORK} == 1 )) ; then
  let "YYYY_NEXT = ${YYYY} +1"
  MM_NEXT=01
  cyear=${YYYY}
elif (( ${IWORK} == 2 )) ; then
  if (( ${MM} == "12" )) ; then
    MM_NEXT=01
    let "YYYY_NEXT = ${YYYY} +1"
  else
    let "MM_NEXT = ${MM} +1"
    YYYY_NEXT=${YYYY}
  fi
  cyear=${YYYY}${MM}
elif (( ${IWORK} == 3 )) ; then
  let "YYYY_NEXT = ${YYYY} +1"
  MM_NEXT=01
  cyear=${YYYY}
  for dnout in "meanflowbin.srv" "outhd_${YYYY}-${MM}-02_hd_meanflow.nc" "outhd_${YYYY}-${MM}-01_hd_meanflow.nc"
  do
    if [ -s $dnout ] ; then
      cdo -f nc4 -z zip_2 settaxis,${YYYY}-01-01,12:00:00,1d -setcalendar,360_day $dnout tmp.nc
      cp tmp.nc $dnout
      if [ -e $dnout ];then
       rm -f tmp.nc
      else
       echo "$dnout is not yet there"
      fi
    fi
  done
fi
#
if [ -s hd_outflow_07.log ] ; then 
 cp hd_outflow_07.log ${HDOUT}/${EXP}_outflow_07_${cyear}.log
 if [ -e ${HDOUT}/${EXP}_outflow_07_${cyear}.log ];then
  rm -f hd_outflow_07.log
 else
  echo "${HDOUT}/${EXP}_outflow_07_${cyear}.log is not yet there"
 fi
fi

if [ -s hd_outflow_99.log ] ; then 
 cp hd_outflow_99.log ${HDOUT}/${EXP}_outflow_99_${cyear}.log
 if [ -e ${HDOUT}/${EXP}_outflow_99_${cyear}.log ];then
  rm -f hd_outflow_99.log
 else
  echo "${HDOUT}/${EXP}_outflow_99_${cyear}.log is not yet there"
 fi
fi
#
cp hdrestart.nc ${HDOUT}/${EXP}_hdrestart_${YYYY_NEXT}${MM_NEXT}.nc
if [ -e ${HDOUT}/${EXP}_hdrestart_${YYYY_NEXT}${MM_NEXT}.nc ];then
  rm -f hdrestart.nc
 else
  echo "${HDOUT}/${EXP}_hdrestart_${YYYY_NEXT}${MM_NEXT}.nc is not yet there"
fi

if [ -e hdrestart.nc ];then
 echo "hdrestart.nc for ${YYYY_NEXT}${MM_NEXT} is not yet copied. Please check!!!"
 exit
fi
#
if [ -s meanflowbin.srv ] ; then
  cdo -b 32 -f nc4 -z zip_2 setgrid,grid_hd.txt meanflowbin.srv ${HDOUT}/${EXP}_meanflow_${cyear}.nc
  rm meanflowbin.srv 
elif [ -s outhd_${YYYY}-${MM}-02_hd_meanflow.nc ] ; then   # *** ncks uses too much memory,cdo increases file size
##  ncks -O -h -v friv outhd_${YYYY}-${MM}-02_hd_meanflow.nc ${HDOUT}/${EXP}_meanflow_${cyear}.nc
  nccopy -V friv,lon,lat,time outhd_${YYYY}-${MM}-02_hd_meanflow.nc ${HDOUT}/${EXP}_meanflow_${cyear}.nc
  rm outhd_${YYYY}-${MM}-02_hd_meanflow.nc
elif [ -s outhd_${YYYY}-${MM}-01_hd_meanflow.nc ] ; then
  nccopy -V friv,lon,lat,time outhd_${YYYY}-${MM}-01_hd_meanflow.nc ${HDOUT}/${EXP}_meanflow_${cyear}.nc
  rm outhd_${YYYY}-${MM}-01_hd_meanflow.nc
fi
#
# Parameter-File-Version als Attribut 
ncatted -O -h -a parameter_file,global,o,c,"$DNPARA" ${HDOUT}/${EXP}_meanflow_${cyear}.nc
#
if (( ${IWORK} == 1 )) || (( ${IWORK} == 3 )) ; then
  cdo monmean ${HDOUT}/${EXP}_meanflow_${YYYY}.nc ${HDOUT}/mon_${EXP}_${YYYY}.nc
fi
if (( ${ICOUPLE} == 2 )) ; then
  cdo -b 32 -f nc4 -z zip_2 copy discharge_on_ocean.nc ${HDOUT}/${EXP}_discharge_on_ocean_${YYYY}.nc
  set +e 
  rm discharge_on_ocean.nc
  set -e 
fi
#
if [ -s ${WRKDIR}/log ] ; then echo "Log directory exists" 
else
  mkdir ${WRKDIR}/log
  echo "Log directory created"
fi

#
# *** Year (and month info)
cat > ${WRKDIR}/log/${EXP}.year << EOF
${YYYY_NEXT}
EOF
if (( ${IWORK} == 2 )) ; then
cat > ${WRKDIR}/log/${EXP}.month << EOF1
${MM_NEXT}
EOF1
fi
#
set +x
echo ${WRKDIR}
cd ${WRKDIR}
date
#
# *** End of simulation or submit next 
if [ ${YYYY_NEXT} -le ${YEND} ] ; then
  echo 'sbatch run_hdyac.ksh YEAR: ' ${YYYY_NEXT}
  sbatch run_hdyac.ksh ${YYYY_NEXT}
fi
if [ ${YYYY} -eq ${YEND} ] ; then
  unset -v verbose
  if ! command -v mailx &> /dev/null
  then
    echo "mailx could not be found"
    exit
  else
    mailx -s "run_hdmodel_${EXP} finished for year ${cyear}" ${HD_CONT} < ${WRKDIR}/log/${EXP}.year
  fi
fi

echo "--- Post-processing for HD finished ---"

exit
#
