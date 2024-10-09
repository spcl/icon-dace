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

! Read configuration parameters as Fortran namelist from an external file.

MODULE mo_aes_rad_nml

  USE mo_aes_rad_config   ,ONLY: aes_rad_config, init_aes_rad_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_aes_rad_nml

  NAMELIST /aes_rad_nml/ aes_rad_config

CONTAINS

  SUBROUTINE nml_read(funit)
    INTEGER, INTENT(in) :: funit
    READ(funit, NML=aes_rad_nml)
  END SUBROUTINE nml_read
  !
  SUBROUTINE nml_write(funit)
    INTEGER, INTENT(in) :: funit
    WRITE(funit, NML=aes_rad_nml)
  END SUBROUTINE nml_write
  !
  SUBROUTINE process_aes_rad_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_aes_rad_config
    !
    CALL process_nml(filename, 'aes_rad_nml', nml_read, nml_write)
    !
  END SUBROUTINE process_aes_rad_nml

END MODULE mo_aes_rad_nml
