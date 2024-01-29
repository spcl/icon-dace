#! /bin/bash

# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

#____________________________________________________________________________________________________
#
# kornblueh1@jwlogin22% srun -N1 -n1 nvidia-smi topo -mp
# 	 GPU0	GPU1	GPU2	GPU3	mlx5_0	mlx5_1	mlx5_2	mlx5_3	CPU Affinity	NUMA Affinity
# GPU0	 X 	SYS	SYS	SYS	PIX	SYS	SYS	SYS	18-23,66-71	3
# GPU1	 SYS	 X 	SYS	SYS	SYS	PIX	SYS	SYS	6-11,54-59	1
# GPU2	 SYS	SYS	 X 	SYS	SYS	SYS	PIX	SYS	42-47,90-95	7
# GPU3	 SYS	SYS	SYS	 X 	SYS	SYS	SYS	PIX	30-35,78-83	5
# mlx5_0 PIX	SYS	SYS	SYS	 X 	SYS	SYS	SYS		
# mlx5_1 SYS	PIX	SYS	SYS	SYS	 X 	SYS	SYS		
# mlx5_2 SYS	SYS	PIX	SYS	SYS	SYS	 X 	SYS		
# mlx5_3 SYS	SYS	SYS	PIX	SYS	SYS	SYS	 X 		
#
# Legend:
#
# X    = Self
# SYS  = Connection traversing PCIe as well as the SMP interconnect between NUMA nodes (e.g., QPI/UPI)
# NODE = Connection traversing PCIe as well as the interconnect between PCIe Host Bridges within
#        a NUMA node
# PHB  = Connection traversing PCIe as well as a PCIe Host Bridge (typically the CPU)
# PXB  = Connection traversing multiple PCIe bridges (without traversing the PCIe Host Bridge)
# PIX  = Connection traversing at most a single PCIe bridge
#____________________________________________________________________________________________________
#
# nvidia-smi:
# NUMA Affinity is the relevant number --cpunodebind= ... 
#
# hwloc lstopo:
# NumaNode is the respective number and visual easier to detect
#____________________________________________________________________________________________________

while getopts n:o:e: argv
do
    case "${argv}" in
        n) tasks_per_node=${OPTARG};;
        o) io_tasks=${OPTARG};;
        e) executable=${OPTARG};;
    esac
done

set -eu

mpi_total_procs=$(wc -l $SLURM_HOSTFILE | awk '{print $1}')

lrank=$SLURM_LOCALID

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^vader,tcp,openib,smcuda"

# need to check in run script that the variables make sense and are
# exported!
(( compute_tasks = mpi_total_procs - io_tasks ))

# splitted the blocks for later use of juwels nodes as I/O server

if (( SLURM_PROCID < (compute_tasks - 1) ))
then

    echo Compute process $SLURM_LOCALID on $(hostname)

    numanode=(3 1 7 5)
    gpus=(0 1 2 3)
    nics=(mlx5_0:1 mlx5_1:1 mlx5_2:1 mlx5_3:1)
    reorder=(0 1 2 3)

    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]}) 

    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]}) 

  export UCX_NET_DEVICES=${nic_reorder[lrank]}
  export CUDA_VISIBLE_DEVICES=${gpus[${reorder[lrank]}]}

  export UCX_RNDV_SCHEME=put_zcopy
  export UCX_RNDV_THRESH=8192

  export UCX_TLS=rc_x,mm,cuda_ipc,cuda_copy,gdr_copy
  export UCX_MEMTYPE_CACHE=n

else

    echo IO process $SLURM_LOCALID on $(hostname)

    numanode=(3 1 7 5)
    nics=(mlx5_0:1 mlx5_1:1 mlx5_2:1 mlx5_3:1)    
    reorder=(0 1 2 3)
 
    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]}) 

    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]}) 
    
    export UCX_NET_DEVICES=${nic_reorder[lrank]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=8192

    export UCX_TLS=rc_x,mm,cuda_ipc,cuda_copy,gdr_copy    
    export UCX_MEMTYPE_CACHE=n
    
fi

numactl --cpunodebind=${numanode_reorder[$lrank]} --membind=${numanode_reorder[$lrank]} $executable

