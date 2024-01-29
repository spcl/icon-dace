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

MODULE mo_aes_cov_nml

  USE mo_aes_cov_config   ,ONLY: aes_cov_config, init_aes_cov_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_aes_cov_nml

  NAMELIST /aes_cov_nml/ aes_cov_config

CONTAINS

  SUBROUTINE nml_read(funit)
    INTEGER, INTENT(in) :: funit
    READ(funit, NML=aes_cov_nml)
  END SUBROUTINE nml_read
  !
  SUBROUTINE nml_write(funit)
    INTEGER, INTENT(in) :: funit
    WRITE(funit, NML=aes_cov_nml)
  END SUBROUTINE nml_write
  !
  SUBROUTINE process_aes_cov_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_aes_cov_config
    !
    CALL process_nml(filename, 'aes_cov_nml', nml_read, nml_write)
    !
  END SUBROUTINE process_aes_cov_nml

END MODULE mo_aes_cov_nml
