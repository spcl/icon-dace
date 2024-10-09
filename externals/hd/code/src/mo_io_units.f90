! mo_io_units.f90 - IO Units
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________

MODULE mo_io_units
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.

  IMPLICIT NONE
  
  PUBLIC
  
  ! This parameter is taken from /usr/include/stdio.h (ANSI C standard). If 
  ! problems with filename length appear, check the before mentioned file.
  
  INTEGER, PARAMETER :: filename_max = 1024
  
  ! Standard I/O-units
  
#ifdef hpux
  INTEGER, PARAMETER :: nerr  = 7     ! error output
#else
  INTEGER, PARAMETER :: nerr  = 0     ! error output
#endif
  INTEGER, PARAMETER :: nlog  = 1     ! standard log file unit
  INTEGER, PARAMETER :: nin   = 5     ! standard input
  INTEGER, PARAMETER :: nout  = 6     ! standard output  
  
  INTEGER, PARAMETER, PRIVATE :: none = -1  ! unit given back, when nothing 
                                            ! in the allowed range is available
  
END MODULE mo_io_units











