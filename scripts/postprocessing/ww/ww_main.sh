#!/bin/ksh

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#
#
typeset          name=${0##*/} ; name=${name%.ne[uw]}

typeset modul_ww=$(dirname $0)/ww_main
modul_ww=${MODUL_WW:-${modul_ww}}

integer rc
typeset modul=$(whence ${modul_ww})
if (( rc != 0 )) ; then
    print -- "modul_ww=${modul_ww} not found! Set MODUL_WW" >&2
    exit 10
fi

#typeset dir='.'
typeset dir=${TMPDIR}
typeset pdatin=h00
integer hstop=12 hinc=1
integer hstart=9999
integer nlev=90
integer k1=0
integer khhl=1
integer iverb=1
typeset sky_cat
typeset -u model='ICON'
typeset rkgrenz1 rkgrenz2
typeset ww_output_gribs
typeset queue='ksh'

while getopts c:C:d:i:k:l:m:o:q:S:s:t:T:v:hD opts
do
    case ${opts} in
        c) dir=${OPTARG}    ;;
        C) sky_cat="${OPTARG}" ;;
        d) pdatin=${OPTARG} ;;
        i) hinc=${OPTARG}   ;;
        k) k1=${OPTARG}     ;;
        l) nlev=${OPTARG}   ;;
        m) model=${OPTARG}   ;;
        o) ww_output_gribs=${OPTARG}   ;;
        q) queue=${OPTARG}   ;;
        s) hstop=${OPTARG}  ;;
        S) hstart=${OPTARG} ;;
        t) rkgrenz1=${OPTARG}  ;;
        T) rkgrenz2=${OPTARG}  ;;
        v) iverb=${OPTARG}  ;;
        D) set -x
           if [[ ${iverb} -lt 1 ]]; then
               iverb=2
           fi
           ;;
        h) cat <<-HELP_EOF

	  ${name} [-c GRIB_dir] [-C sky_category]
	             [-d pdatin] [-i hinc] [-k k1] [-l nlev] [-m model]
	             [-S hstart] [-s hstop] [-v verbosity] [-D] [-h]

	Calculate siginificant weather, WW, with ww_main.

	The required input variables must be in files ww_input_STEP where STEP
	is the forecast step in hours with 4 digits.

	Options:
	   -c GRIB_dir    directory with input file ww_input_STEP. Default: ${dir}
	   -C category    SKY data base category. If category is set the data is retrieved from the data base.
	   -d pdatin      initital date and time. Default: ${pdatin}
	   -i hinc        step increment [h]. Default: ${hinc}
	   -k k1          Extract levels k1/to/nlev from data base. Default: 31 for ICON, 1 for ICON-EU.
	   -l nlev        number of levels. Default: ${nlev}
	   -m model       Use constant for this model. Default: ${model}
	                  Use TEST for new values of rkgrenz1 or rkgrenz2.
	   -o output_gribs  Name of output GRIB file. Default: ww_\$model.grb
	   -q queue       Submit to job to queue on lc. Default: ${queue} (run on current node)
	   -S hstart      start step [h]. Default: \${hinc}
	   -s hstop       last step [h].  Default: ${hstop}
	   -v verbosity   verbosity level. Default: ${iverb}
	   -t rkgrenz1    precip. limit for thunderstorms only if model=TEST
	   -T rkgrenz2    precip. limit for strong thunderstorms only if model=TEST
	   -D             Debug option. Sets verbosity to at least 2
	   -h             This help

	Examples:

	WW for operational ICON:
	${name} -c \$TMPDIR/WWicon -C icogl130l90_main_fc_rout -d h00 -i 1 -s 78 -m ICON -l 90 -q lc_big

	WW for operational ICON-EU nest:
	${name} -c \$TMPDIR/WWicon -C icreu_main_fc_rout -d h00 -i 1 -s 78 -m ICON -l 60 -q lc_big

	WW for ICON with new tuning parameters (data was already retrieved from data base)
	${name} -c \$TMPDIR/WWicon -d h00 -i 1 -s 78 -m TEST -l 90 -q lc_big -t 2.5 -T 6.0

	As above but for 3 hourly intervalls:
	${name} -c \$TMPDIR/WWicon -d h00 -i 3 -s 96 -m TEST -l 90 -q lc_big -t 2.5 -T 6.0

	  Helmut Frank, 31.07.2017
	HELP_EOF

#	WW for operational COSMO-EU:
#	${name} -c \$TMPDIR/wwLME -C c2_main_fc_rout -d h00 -i 1 -s 6 -m LME -l 40

#	WW for pre-operational COSMO-EU:
#	${name} -c \$TMPDIR/wwLMEp -C c2_main_fc_para -d h00 -i 1 -s 6 -m LME -l 40

           exit
           ;;
    esac
done

#  Set missing options
typeset cymdg
cymdg=$(datconv -Cy ${pdatin})
if [[ ${dir} = .* ]]; then
    dir=${PWD}/${dir}
fi
if [[ ! -d ${dir} ]]; then
    mkdir -p ${dir}
fi

if [[ ${hstart} -eq 9999 ]]; then
    hstart=${hinc}
fi

case ${model} in
    ICON|ICEU|TEST) rkgrenz1=${rkgrenz1:-0.75}
          rkgrenz2=${rkgrenz2:-2.0}
          ;;
    LME)  rkgrenz1=${rkgrenz1:-0.75}
          rkgrenz2=${rkgrenz2:-2.0}
          ;;
    LMK)  rkgrenz1=${rkgrenz1:-0.60}
          rkgrenz2=${rkgrenz2:-2.0}
          ;;
    GME)  rkgrenz1=${rkgrenz1:-0.25}
          rkgrenz2=${rkgrenz2:-1.0}
          ;;
esac

ww_output_gribs=${ww_output_gribs:-ww_${model}${cymdg}.grb}

typeset bank
case ${sky_cat} in
    *_rout)  bank='roma'  ;;
    *_para*) bank='parma' ;;
    *_vera*) bank='vera'  ;;
         *)  bank='numex' ;;
esac

typeset job_file out_file

if [[ ${queue} = 'ksh' ]]; then
    job_file="${dir}/ww_${model}.sh$$"
else
    job_file="${dir}/ww_${model}.job$$"
    out_file="ww_${model}_${sky_cat}${cymdg}.out"
fi
out_file=${WW_OUT_FILE:-${out_file}}

cat > ${job_file} <<EOFJOB
#PBS -S /bin/ksh
#PBS -N  WW_${model}
#PBS -q  ${queue}
##PBS -l cput=${cpu_time}
#PBS -l walltime=3600
#PBS -o ${out_file}
#PBS -j oe
#PBS -W umask=022
#PBS -m a

cd ${dir}

if [[ -n "${sky_cat}" ]]; then
#
#  Retrieve data from the data base
#
    typeset req_file=req4ww.\$$

    integer nlev1=$((${nlev}+1))
    integer nm1=$((${nlev}-1))
    integer nm10=$((${nlev}-10))
    integer nm30=$((${nlev}-30))
    integer k1=${k1}
    if (( \$k1 == 0 )); then
        if (( ${nlev} >= 90 )); then
            k1=31
            khhl=31
        elif [[ ${model} = 'LME' ]]; then
            k1=11
            khhl=11
        elif [[ ${model} = 'ICEU' ]]; then
            khhl=1
        fi
    fi

    integer vv
    if [[ ${hstart} -gt ${hinc} ]]; then
        vv=\$((${hstart}-${hinc}))
    else
        vv=${hstart}
    fi
    typeset -Z4 vvvv=\${vv}

    cat > \${req_file} <<REQ_EOF
reqColl proc=parallel
read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=0 bin f=${dir}/ww_input_hhl info=countPlus p=HHL lv1=\${khhl}/to/\${nlev1} lv=genv
REQ_EOF

    if [[ ${hstart} -gt ${hinc} ]]; then
        integer v0=\$((${hstart}-${hinc}))
        print -- "read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=\${v0} bin f=${dir}/ww_input_hhl info=countPlus p=rain_gsp,rain_con,snow_gsp,snow_con" >> \${req_file}
    fi

    integer vv=${hstart}
    while [[ \${vv} -le ${hstop} ]]
    do
        vvvv=\${vv}
        cat >> \${req_file} <<REQ_EOF
read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=\${vv} bin f=${dir}/ww_input_\${vvvv} info=countPlus p=T,P lv1=\${k1}/to/${nlev}   lv=genv
read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=\${vv} bin f=${dir}/ww_input_\${vvvv} info=countPlus p=U,V,QC lv1=\${nm10}/to/${nlev} lv=genv
read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=\${vv} bin f=${dir}/ww_input_\${vvvv} info=countPlus p=QV lv1=\${nm30}/to/${nlev} lv=genv
read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=\${vv} bin f=${dir}/ww_input_\${vvvv} info=countPlus p=CLC       lv1=\${nm1}/to/${nlev}  lv=genv
read db=${bank} cat=${sky_cat} d=${cymdg} s[h]=\${vv} bin f=${dir}/ww_input_\${vvvv} info=countPlus p=ps,t_2m,td_2m,t_g,clct,clcm,u_10m,v_10m,rain_gsp,rain_con,snow_gsp,snow_con,hbas_con,htop_con
REQ_EOF
        vv=\$((\${vv}+${hinc}))
    done

    cat \${req_file}

    sky -v \${req_file}
    rm \${req_file}

fi

#  Write the namelist
cat > NAMELIST_WW <<EOF_WW
&ww_namelist
 ymodel='${model}',
 ycat="${dir}",
 ydate='${cymdg}',
 hstart=${hstart},
 hstop=${hstop},
 hinc=${hinc},
 ke=${nlev}, 
 verbosity=${iverb},
 ww_output="${ww_output_gribs}",
/
&ww_tune_nml
 rkgrenz1  = ${rkgrenz1},
 rkgrenz2  = ${rkgrenz2},
!rf_fog    = 95.,   
!clc_fog   = 99.,  
!rain_l_m  = 2.5,   
!rain_m_s  = 10.0,   
!snow_l_m  = 1.00, 
!snow_m_s  = 5.0, 
!rash_lm_s = 2.5,
!rash_s_vs = 20.0,
!snsh_l_ms = 1.00, 
!driz_l_m  = 0.1, 
!driz_m_s  = 0.5,
!drif_l_ms = 0.1,
!raif_l_ms = 2.5,
!rgdiff_th1= 0.015, 
!rgdiff_th2= 0.050 
/
EOF_WW

# Run the WW program

#${modul_ww}
${modul}

exit
EOFJOB

# Run the job

if [[ ${queue} = 'ksh' ]]; then
    ksh ${job_file}
    rm ${job_file}
else
    cat ${job_file}
    qsub ${job_file}
fi
