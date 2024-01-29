#!/bin/ksh 

function set_cluster {
 ##
 echo " Running on " >&2
 HPC=`(hostname || uname -n) 2>/dev/null | sed 1q`
 case x"$HPC" in #(
   xuc1*) :
     echo "...UC1 at KIT"; CENTER="IMK" 
     input_folder="/lsdf/kit/imk/projects/icon/TESTSUITE"
     FILETYPE="4" 
     output_folder="${WORK}/TESTSUITE_OUTPUT"
     icon_data_poolFolder=/lsdf/kit/imk/projects/icon/INPUT/AMIP/amip_input
     aer_opt="${icon_data_poolFolder}"
     ;;
   xhk*) :
     echo "...HoreKa at KIT"; CENTER="IMK"
     input_folder="/lsdf/kit/imk/projects/icon/TESTSUITE"
     FILETYPE="5" 
     ws=$(ws_list -s)
     if [[ "${ws}" == "" ]]; then
         echo "No workspaces found! Help how to create one: https://www.nhr.kit.edu/userdocs/horeka/filesystems/#create-workspace"
         exit
     else
         ws_id=$(echo $ws | awk '{print $1}')
     fi
     WORK=/hkfs/work/workspace/scratch/$(whoami)-$ws_id
     output_folder="${WORK}/TESTSUITE_OUTPUT"
     icon_data_poolFolder=/lsdf/kit/imk/projects/icon/INPUT/AMIP/amip_input
     aer_opt="${icon_data_poolFolder}"
     ;;
   xxce*) :
     echo "...XCE at DWD"; CENTER="DWD"
         . /opt/modules/3.2.10.3/init/ksh
         input_folder="${SCRATCH}/TESTSUITE_INPUT/"
         output_folder="${WORK}/TESTSUITE_OUTPUT" 
	 FILETYPE="4" 
	 if [[ ":${PE_ENV}" != ':CRAY' ]] ; then
          module unload libdwd grib_api
          module unload cray-netcdf netcdf perftools stat
          module unload craype-hugepages2M
          case "{PE_ENV}" in
             GNU)   prgenv_="PrgEnv-gnu"   ;;
             INTEL) prgenv_="PrgEnv-intel" ;;
          esac
          module swap ${prgenv_} PrgEnv-cray
       fi
       module load cray-netcdf perftools stat
       module load libdwd
       module load grib_api
       module list
       ##       
       ;;
   xlce*) :
     echo "...LCE at DWD"; echo "  ERROR: Testsuite has to be run and compiled on XCE !!!", exit 201 ;; #(
   xjuwels*) :
     echo "...juwels at FZJ"; CENTER="FZJ" 
     icon_data_poolFolder=/gpfs/homea/hka21/hka211/INPUT/AMIP/atm_amip_base/
	 FILETYPE="4" 
     output_folder="${WORK}/TESTSUITE_OUTPUT"
	 ;; #(
   xmlogin*) :
     echo "...MISTRAL at DKRZ"; CENTER="DKRZ"
     icon_data_poolFolder=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/ 
     input_folder="${WORK}/TESTSUITE/"
	 FILETYPE="4" 
     output_folder="/scratch/b/${USER}/TESTSUITE_OUTPUT/"
     aer_opt="/pool/data/ICON/grids/public/mpim/independent"
	 ;;

   mlogin*)
     echo "...MISTRAL at DKRZ"; CENTER="DKRZ" 
     input_folder="/pfs/imk/ICON/TESTSUITE/"
	 FILETYPE="4" 
     output_folder="${SCRATCH}/TESTSUITE_OUTPUT"
	 icon_data_poolFolder=/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/
     aer_opt="/pool/data/ICON/grids/public/mpim/independent"
	 ;;
   xl40053*)
     echo "...LEVANTE at DKRZ"; CENTER="DKRZ"
     input_folder="${WORK}/../testcases/"
   FILETYPE="4"
     output_folder="${SCRATCH}/TESTSUITE_OUTPUT"
     aer_opt="/pool/data/ICON/grids/public/mpim/independent"
   ;;
   *) :
     echo "...unknown HPC" ; exit 202 ;; #(
 esac
 ##
 ICON_FOLDER=$(pwd)/../..
 ART_FOLDER=$ICON_FOLDER/externals/art


 }

function read_setup {
#declare -a experiment
#declare -a script
#awk -va="experiment=$experiment; script=$script" 'BEGIN{n=1;}
#     { FS=" "; $0=$0; script[n]=$1; experiment[n]=$2;
#		 print n $1 script[n];
#
#     }' 'setup_file.txt' 
#awk -F '{ print $1, $2}' 'setup_file.txt' | read var1 var2
#echo $var1 ' ' $var2
i=0

for line in `awk '{print $2}' $1`; do
  experiment[$i]=$line
  i=`expr $i + 1`
done
i=0

for line in `awk '{print $1}' $1`; do
  script[$i]=$line
  i=`expr $i + 1`
done

i=0
for line in `awk '{print $3}' $1`; do
  config[$i]=$line
  i=`expr $i + 1`
done

i=0
for line in `awk '{print $4}' $1`; do
  walltime[$i]=$line
  i=`expr $i + 1`
done


i=0
for line in `awk '{print $5}' $1`; do
  queue[$i]=$line
  i=`expr $i + 1`
done


i=0
for line in `awk '{print $6}' $1`; do
  nodes[$i]=$line
  i=`expr $i + 1`
done



number_of_experiments=$i

 }

 function usage
{
    echo "This is a script for creating buildbot / testsuite / experiment scripts"
    echo "for different clusters"
    echo "Usage:"
    echo "./create_testsuite_runs --input_file=YOUR_FILE"
    echo ""
}

function create_header
{
d=`date -d today +%Y%m%d`
complete_output_folder=${output_folder}/${d}/$1
OUTDIR=$complete_output_folder
output_script=$ICON_FOLDER/run/checksuite.icon-kit/runscripts/$1.run
icon_data_poolFolder=$icon_data_poolFolder
EXPERIMENT=$1
lart=$lart

cat > $output_script << EOF
#!/bin/bash
CENTER=$CENTER
basedir=$ICON_FOLDER
icon_data_poolFolder=$icon_data_poolFolder
aer_opt=$aer_opt
EXPNAME=atm_amip_test_kit
OUTDIR=$complete_output_folder
ICONFOLDER=$ICON_FOLDER
ARTFOLDER=$ART_FOLDER
INDIR=$input_folder
EXP=$EXPERIMENT
lart=$lart

FILETYPE=$FILETYPE
COMPILER=intel
restart=.False.
read_restart_namelists=.False.



# Remove folder ${EXP} from OUTDIR for postprocessing output
OUTDIR_PREFIX=`dirname \${OUTDIR}`

# Create output directory and go to this directory

if [ ! -d \$OUTDIR ]; then
    mkdir -p \$OUTDIR
fi

cd \$OUTDIR


EOF

}
function create_footer
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/runscripts/$1.run

cat >> $output_script << EOF
	
cp $ICON_FOLDER/bin/icon ./icon.exe
EOF
 case x"$HPC" in #(
      xjuwels*)
cat >> $output_script << EOF
	
cat > job_ICON << ENDFILE
#!/bin/bash -x
#SBATCH --nodes=4
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpi-out.%j
#SBATCH --error=mpi-err.%j
#SBATCH --time=00:15:00
#SBATCH --partition=batch
module load Intel
module load ParaStationMPI
module load Python/2.7.14
module load netCDF-Fortran/4.4.4
module load netCDF/4.6.1
srun ./icon.exe

ENDFILE


chmod +x job_ICON
sbatch job_ICON

EOF
;;

   xhk*)
cat >> $output_script << EOF
	   
cat > job_ICON << ENDFILE
#!/bin/bash -x
#SBATCH --$3
#SBATCH --time=$2
#SBATCH --ntasks-per-node=76
#SBATCH --partition=$4


$5 

mpirun --bind-to core --map-by core --report-bindings ./icon.exe

ENDFILE

chmod +x job_ICON
sbatch job_ICON

EOF
;;

  xxce*)
cat >> $output_script << EOF
PBS_O_WORKDIR='\${PBS_O_WORKDIR}'
cat > icon.job <<ENDFILEDWD
#!/bin/ksh
#-----------------------------------------------------------------------------
#PBS -q xc_norm_h
#PBS -l select=6:ompthreads=4
#PBS -l place=scatter
#PBS -l walltime=00:60:00
#PBS -j oe
# ----------------------------------------------------------------------
set -x

 export OMP_SCHEDULE="static"
 export OMP_DYNAMIC="false"
#export OMP_STACKSIZE=256M
 export ATP_ENABLED=1
 export MPICH_RMA_OVER_DMAPP=1

# for PBS change to directory where job was submitted
# (without this job is started in HOME)
 if [[ -n "\${PBS_O_WORKDIR}" ]] ; then
    cd "\${PBS_O_WORKDIR}"
 fi

# run ICON
 aprun  -n 72 -N 12 -j 2 -d 4 -m 3g icon.exe

ENDFILEDWD

      # submit the job
      qsub icon.job

      ##       
EOF
;;
   xl40053*)
account_id=$(id -gn)
cat >> $output_script << EOF

cat > job_ICON << ENDFILE
#!/bin/bash -x
#SBATCH --account=$account_id

#SBATCH --partition=$4
#SBATCH --$3
#SBATCH --exclusive
#SBATCH --ntasks-per-node=128
#SBATCH --time=$2

ulimit -s 102400
ulimit -c 0

$5
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl=self
export OMPI_MCA_osc="pt2pt"
export UCX_IB_ADDR_TYPE=ib_global

export OMPI_MCA_coll="^ml,hcoll"
export OMPI_MCA_coll_hcoll_enable="0"
export HCOLL_ENABLE_MCAST_ALL="0"
export HCOLL_MAIN_IB=mlx5_0:1
export UCX_NET_DEVICES=mlx5_0:1
export UCX_TLS=mm,knem,cma,dc_mlx5,dc_x,self
export UCX_UNIFIED_MODE=y
export HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_io="romio321"
export UCX_HANDLE_ERRORS=bt
module load eccodes/2.21.0-intel-2021.5.0

srun -l --cpu_bind=cores --distribution=plane=64 \
        --propagate=STACK,CORE ./icon.exe

ENDFILE
                 
chmod +x job_ICON
sbatch job_ICON  

EOF
;;
   xmlogin*)
account_id=$(id -gn)
cat >> $output_script << EOF
	   
cat > job_ICON << ENDFILE
#!/bin/bash -x
#SBATCH --account=$account_id

#SBATCH --partition=$4
#SBATCH --$3
#SBATCH --exclusive
#SBATCH --ntasks-per-node=24
#SBATCH --time=$2


ulimit -s 102400

$5 
export OMPI_MCA_pml=cm
export OMPI_MCA_mtl=mxm
export OMPI_MCA_mtl_mxm_np=0
export MXM_RDMA_PORTS=mlx5_0:1
export MXM_LOG_LEVEL=ERROR
# Disable GHC algorithm for collective communication
export OMPI_MCA_coll=^ghc

export MXM_HANDLE_ERRORS=bt
export UCX_HANDLE_ERRORS=bt

export OMPI_MCA_coll=^fca
export OMPI_MCA_coll_hcoll_enable=1
export OMPI_MCA_coll_hcoll_priority=95
export OMPI_MCA_coll_hcoll_np=8
export HCOLL_MAIN_IB=mlx5_0:1
export HCOLL_ENABLE_MCAST=1
export HCOLL_ENABLE_MCAST_ALL=1

export HCOLL_ML_DISABLE_BARRIER=1
export HCOLL_ML_DISABLE_IBARRIER=1
export HCOLL_ML_DISABLE_BCAST=1
export HCOLL_ML_DISABLE_REDUCE=1

srun -l --propagate=STACK --cpu_bind=cores \
  --distribution=block:cyclic ./icon.exe

ENDFILE

chmod +x job_ICON
sbatch job_ICON

EOF
;;

esac


chmod +x $output_script
}
function create_body

{
 cat ${ART_FOLDER}/$1 >> $output_script
}


function create_lart_loop_start
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/runscripts/$1.run
cat >> $output_script << EOF
	
for iter in 1 2 ; do


if [ \$iter -eq 1 ] ; then

lart=.False.
if [ ! -d $OUTDIR/wo_ART ]; then
    mkdir -p $OUTDIR/wo_ART
fi

cd $OUTDIR/wo_ART
OUTDIR=$OUTDIR/wo_ART
else

lart=.True.
if [ ! -d $OUTDIR/w_ART ]; then
    mkdir -p $OUTDIR/w_ART
fi

cd $OUTDIR/w_ART

OUTDIR=$OUTDIR/w_ART

fi

EOF
}

function create_lart_loop_end
{
output_script=$ICON_FOLDER/run/checksuite.icon-kit/runscripts/$1.run
cat >> $output_script << EOF

done

EOF
}

function check_action
{
output_script=runscripts/$2.run
	
	case x"$1" in 
		xrun*)
          ./$output_script
		  ;;
		*none*)
		  ;;
  esac
}

function read_configure
{
. $(pwd)/../set-up.info

output="module load ${use_load_modules}"
#i=0
#for element in $(echo ${use_load_modules}|sed -e " ");do
#   module_str[$i]="module load ${element}"
#  i=`expr $i + 1`
#done
#
#for item in ${!module_str[*]};do
#	 printf " %s \n" "${module_str[${item}]}"
#done

 }

