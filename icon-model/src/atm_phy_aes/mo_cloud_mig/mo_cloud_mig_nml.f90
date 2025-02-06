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

MODULE mo_cloud_mig_nml

  USE mo_cloud_mig_config ,ONLY: cloud_mig_config, init_cloud_mig_config
  USE mo_process_nml      ,ONLY: process_nml
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: process_cloud_mig_nml

  NAMELIST /cloud_mig_nml/ cloud_mig_config

CONTAINS

  SUBROUTINE process_cloud_mig_nml(filename)
    !
    CHARACTER(LEN=*), INTENT(in) :: filename
    !
    CALL init_cloud_mig_config
    !
    CALL process_nml(filename, 'cloud_mig_nml', nml_read, nml_write)
    !
  CONTAINS
    !
    SUBROUTINE nml_read(funit)
      INTEGER, INTENT(in) :: funit
      READ(funit, NML=cloud_mig_nml)
    END SUBROUTINE nml_read
    !
    SUBROUTINE nml_write(funit)
      INTEGER, INTENT(in) :: funit
      WRITE(funit, NML=cloud_mig_nml)
    END SUBROUTINE nml_write
    !
  END SUBROUTINE process_cloud_mig_nml

END MODULE mo_cloud_mig_nml
