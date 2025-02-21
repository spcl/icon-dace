! Module acc_device_management
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_acc_device_management

!------------------------------------------------------------------------------
!
! Description : 
!  This module provides routines to set up the GPU device as well wrapper
!  for cuda function calls. The module is empty if preprocessor variable 
!  _OPENACC is not set.
!  Note: routines using STOP in case of error instead of model_abort
!  as they may be called before mpi init
!  
! Routines contained:
!  - cudaGetErrorString     : get cuda error
!  - cudaGetDeviceCount     : get number of devices
!  - cudaSetDevice          : set device
!  - cudaDeviceReset        : reset device
!  - cudaMemGetInfo         : get memory usage
!  - cudaDeviceSynchronize  : halts CPU/Host thread until GPU has processed all 
!                             tasks
!  - my_cudaErrorCheck      : check and print cuda error if any
!  - initAccDevice          : init device
!  - checkAccDevice         : check that device is correctly initialized
!  - setAccDevice           : set device to current thread
!  - runSmallAccKernel      : run small test OpenACC kernel
!  - finalizeAccDevice      : finalize device
!  - printGPUMem            : print current GPU usage
!
!------------------------------------------------------------------------------


END MODULE mo_acc_device_management
