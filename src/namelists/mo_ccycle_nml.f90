! Read configuration parameters as Fortran namelist from an external file.
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

MODULE mo_ccycle_nml

  USE mo_ccycle_config    ,ONLY: ccycle_config, init_ccycle_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_ccycle_nml

  NAMELIST /ccycle_nml/ ccycle_config

CONTAINS

  SUBROUTINE nml_read(funit)
    INTEGER, INTENT(in) :: funit
    READ(funit, NML=ccycle_nml)
  END SUBROUTINE nml_read
  !
  SUBROUTINE nml_write(funit)
    INTEGER, INTENT(in) :: funit
    WRITE(funit, NML=ccycle_nml)
  END SUBROUTINE nml_write
  !
  SUBROUTINE process_ccycle_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_ccycle_config
    !
    CALL process_nml(filename, 'ccycle_nml', nml_read, nml_write)
    !
    !$ACC ENTER DATA COPYIN(ccycle_config)
  END SUBROUTINE process_ccycle_nml

END MODULE mo_ccycle_nml
