! mo_kind.f90 - Precision settings
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!_________________________________________
!>
!!  This routine comprises excerps of the corresponding routine (year 2017) from MPI-ESM, 
!!  the Earth System Model of the Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
!!  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
!!  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
!!  doi: 10.1029/2018MS001400.
!!
!!
MODULE mo_kind

  ! L. Kornblueh, MPI-M, August 2001, added working precision and comments 

  IMPLICIT NONE

  PUBLIC

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15  
  !                   exponent  = 37   exponent  =  307 
  !
  ! Most likely this are the only possible models.

  ! Floating point section: 

  INTEGER, PARAMETER :: ps = 6
  INTEGER, PARAMETER :: rs = 37

  INTEGER, PARAMETER :: pd = 12
  INTEGER, PARAMETER :: rd = 307

  INTEGER, PARAMETER :: pi4 = 9
  INTEGER, PARAMETER :: pi8 = 14

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(ps,rs) !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(pd,rd) !< double precision
  
  ! Floating point working precision

  INTEGER, PARAMETER :: wp = dp   
  
  ! Integer section
  ! ---------------
  
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(pi4)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
  !
  !--------------------------------------------------------------------

END MODULE mo_kind

