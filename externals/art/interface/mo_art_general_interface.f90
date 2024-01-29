!
! Provides the interface for ART-routines to ICON parts
!
! This module provides an interface to ICON parts.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
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

MODULE mo_art_general_interface
  
  USE mo_cdi,                           ONLY: DATATYPE_FLT32, DATATYPE_FLT64
  USE mo_io_config,                     ONLY: lnetcdf_flt64_output

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_general_interface'

  PUBLIC  :: getNetcdfPrecision 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
FUNCTION getNetcdfPrecision() RESULT(flt_prec)

  INTEGER :: flt_prec

  flt_prec = MERGE(DATATYPE_FLT64,DATATYPE_FLT32,lnetcdf_flt64_output)

END FUNCTION getNetcdfPrecision 
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_general_interface
