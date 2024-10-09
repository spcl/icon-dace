#!/bin/ksh
#
# hd_post.ksh - Postprocessing script to store HD outputs and restarts 
# 
# Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
# SPDX-License-Identifier: Apache-2.0
# See ./LICENSES/ for license information
#
# Authors: Stefan Hagemann
# Contact: <stefan.hagemann@hereon.de>
#_________________________________________
#
#
#################################################
# store output and restarts
#################################################
#
# The experiment number must be passed as an argument
EXP=$1
#
if [ -e setuserdir.com ] ; then
  source ./setuserdir.com
elif [ -e ../scr/setuserdir.com ] ; then
  source ../scr/setuserdir.com
else
  echo 'Script setuserdir.com does not exist in HDMAIN/scr!'
  exit
fi
#
if [ -e ${HDMAIN}/${EXP}/hd_run_settings_${EXP}.ksh ] ; then    
  export ICALL=2  
  source ${HDMAIN}/${EXP}/hd_run_settings_${EXP}.ksh
  echo "Script hd_run_settings_${EXP}.ksh was read"
else
  echo "Script hd_run_settings_${EXP}.ksh does not exist in directory ${HDMAIN}/${EXP}"
  exit
fi

HDOUT=${HDDIR}/${EXP}/out          # HD Output dir

cd $HDDIR/$EXP

typeset -Z2 MM_NEXT
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
      mv tmp.nc $dnout
    fi
  done
elif (( ${IWORK} == 4 )) ; then
  if (( ${YYYY} < $YEND )) ; then
    let "YYYY_NEXT = ${YYYY} +1"
    MM_NEXT=01
    cyear=${YYYY}
  else
    YYYY_NEXT=${YYYY}
    if (( ${nday_final} <= 31 )) ; then MM_NEXT=02
    elif (( ${nday_final} <= 60 )) ; then MM_NEXT=03
    elif (( ${nday_final} <= 91 )) ; then MM_NEXT=04
    elif (( ${nday_final} <= 121 )) ; then MM_NEXT=05
    elif (( ${nday_final} <= 152 )) ; then MM_NEXT=06
    elif (( ${nday_final} <= 182 )) ; then MM_NEXT=07
    elif (( ${nday_final} <= 213 )) ; then MM_NEXT=08
    elif (( ${nday_final} <= 244 )) ; then MM_NEXT=09
    elif (( ${nday_final} <= 274 )) ; then MM_NEXT=10
    elif (( ${nday_final} <= 305 )) ; then MM_NEXT=11
    elif (( ${nday_final} <= 335 )) ; then MM_NEXT=12
    fi
    cyear=${YYYY}_${nday_final}days
  fi
fi
echo "Next year was calculated"

#
if [ -s hd_outflow_07.log ] ; then mv hd_outflow_07.log ${HDOUT}/${EXP}_outflow_07_${cyear}.log ; fi
if [ -s hd_outflow_99.log ] ; then mv hd_outflow_99.log ${HDOUT}/${EXP}_outflow_99_${cyear}.log ; fi
#
mv hdrestart.nc ${HDOUT}/${EXP}_hdrestart_${YYYY_NEXT}${MM_NEXT}.nc
echo "HD Log and Restart files were moved to ${HDOUT}"
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
echo "HD discharge file was moved (using nccopy) to ${HDOUT}"
#
# Bias corrected outflows?
if [ -s outhd_${YYYY}-${MM}-01_hd_bcflow.nc ] ; then
  # Calculate daily mean bc flow
  cdo -f nc4 -z zip_2 daymean -selvar,friv_bc outhd_${YYYY}-${MM}-01_hd_bcflow.nc ${HDOUT}/${EXP}_bcflow_${cyear}.nc
  ncrename -h -O -v friv_bc,friv ${HDOUT}/${EXP}_bcflow_${cyear}.nc
  ncatted -O -h -a parameter_file,global,o,c,"$DNPARA" ${HDOUT}/${EXP}_bcflow_${cyear}.nc
  rm outhd_${YYYY}-${MM}-01_hd_bcflow.nc
fi
#
# Parameter-File-Version als Attribut 
ncatted -O -h -a parameter_file,global,o,c,"$DNPARA" ${HDOUT}/${EXP}_meanflow_${cyear}.nc
# Start-File als Attribut in first year
if (( ${INEU} == 1 )) ; then
  ncatted -O -h -a start_file,global,o,c,"${HDFILE}/${HDSTART}" ${HDOUT}/${EXP}_meanflow_${cyear}.nc
fi
#
if (( ${IWORK} == 1 )) || (( ${IWORK} == 3 )) ; then
  cdo -f nc4 -z zip_2 monmean ${HDOUT}/${EXP}_meanflow_${YYYY}.nc ${HDOUT}/mon_${EXP}_${YYYY}.nc
  echo "Monthly mean HD discharge file was generated in ${HDOUT}"
elif (( ${IWORK} == 4 )) ; then
  cdo -f nc4 -z zip_2 monmean ${HDOUT}/${EXP}_meanflow_${YYYY}.nc ${HDOUT}/mon_${EXP}_${cyear}.nc
  echo "Monthly mean HD discharge file was generated in ${HDOUT}"
fi
if (( ${ICOUPLE} == 2 )) ; then
  mv discharge_on_ocean.nc ${HDOUT}/${EXP}_discharge_on_ocean_${YYYY}.nc
  echo "Discharge on ocean file was moved to ${HDOUT}"
fi
#
if [ -s ${HDMAIN}/log ] ; then echo "Log directory exists" 
else
  mkdir ${HDMAIN}/log
  echo "Log directory created"
fi

#
# *** Year (and month info)
cat > ${HDMAIN}/log/${EXP}.year << EOF
${YYYY_NEXT}
EOF
if (( ${IWORK} == 2 )) || (( ${IWORK} == 4 )); then
cat > ${HDMAIN}/log/${EXP}.month << EOF1
${MM_NEXT}
EOF1
fi
#
set +x
echo ${HDMAIN}/${EXP}
cd ${HDMAIN}/${EXP}
ls -al run_hdmodel_${EXP}.ksh
date
#
# *** End of simulation or submit next 
if [ ${YYYY_NEXT} -le ${YEND} ] ; then
  if (( ${IWORK} != 4 )) || (( ${YYYY} != ${YEND} )) ; then
    echo 'YEAR: ' ${YYYY_NEXT}
    sbatch run_hdmodel_${EXP}.ksh
  fi
fi
if [ ${YYYY} -eq ${YEND} ] ; then

  if ! command -v mailx  >/dev/null 2>&1
  then
    echo "mailx could not be found"
    exit
  else
    echo "Use mailx to email ${USER_EMAIL}"
    NODE=$(hostname)
    MAIL_SUBJECT="HD run_hdmodel_${EXP} finished for year ${cyear}"
    MAIL_ADDRESS="${USER_EMAIL}"
    MAIL_BODY="run_hdmodel_${EXP} finished for year ${cyear}"

    mailx -s "$MAIL_SUBJECT" "$MAIL_ADDRESS" <<< "$MAIL_BODY"
  fi
fi

echo "--- Post-processing for HD finshed ---"

exit
#
