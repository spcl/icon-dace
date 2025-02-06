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

! This module contains debugging utilities, especially subroutines
! for writing REAL arrays to NetCDF files (for debugging purposes).

MODULE mo_util_debug

! enable the following directive for disabling
! NetCDF dumps at compile time:
!define DISABLE_DUMP 1

  !
  ! debugging utilities
  !
  USE mo_kind,           ONLY: wp
  USE mo_util_string,    ONLY: int2string
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH

  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: dump_array_to_netcdf
  PUBLIC :: debug_step
  PUBLIC :: ldebug_enable

  ! The global variables "debug_step", "ldebug_enable" are useful when debugging output is
  ! desired only for certain steps inside a loop. For example, one may set "debug_step" to 
  ! the current iteration and use this value at some other place (where the
  ! original counter is not available).

  INTEGER :: debug_step    = 0       !< global counter
  LOGICAL :: ldebug_enable = .TRUE.  !< enabling/disabling dumps (during runtime)

  INTERFACE dump_array_to_netcdf
    MODULE PROCEDURE dump_array_to_netcdf_1d
    MODULE PROCEDURE dump_array_to_netcdf_2d
    MODULE PROCEDURE dump_array_to_netcdf_3d
  END INTERFACE
  !
CONTAINS

  !> Dumps real array to NetCDF file
  !
  !  Note: When using this SR, take care of OpenMP regions!

  SUBROUTINE dump_array_to_netcdf_1d(zfilename, p_array)
    CHARACTER (len=*), INTENT(IN) :: zfilename
    REAL(wp)         , INTENT(IN) :: p_array(:)
    ! local variables
    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: idim, ncfile, ncid_var, &
         &     ncid_dim(ndims), icount(ndims)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_util_debug:dump_array_to_netcdf_1d'


    IF (.NOT. ldebug_enable) RETURN

#ifndef DISABLE_DUMP
    WRITE (*,*) "Dumping ", zfilename   
    ! create NetCDF file:
    CALL nf(nf90_create("00_"//TRIM(zfilename)//"_"//TRIM(int2string(debug_step))//".nc", &
      &               nf90_clobber, ncfile), routine)
    ! create dimensions:
    DO idim=1,ndims
      CALL nf(nf90_def_dim(ncfile, 'dim'//int2string(idim), SIZE(p_array,idim), &
        &     ncid_dim(idim)), routine)
      icount(idim) = SIZE(p_array,idim)
    END DO
    ! create variable:
    CALL nf(nf90_def_var(ncfile, "var", NF90_DOUBLE, ncid_dim, ncid_var), routine)
    ! End of definition mode
    CALL nf(nf90_enddef(ncfile), routine)
    ! put data:
    CALL nf(nf90_put_var(ncfile, ncid_var, p_array, &
      &                  (/ (1, idim=1,ndims) /), icount), &
      &     routine)
    ! close file
    CALL nf(nf90_close(ncfile), routine)
#endif

  END SUBROUTINE dump_array_to_netcdf_1d


  !> Dumps real array to NetCDF file
  !
  !  Note: When using this SR, take care of OpenMP regions!

  SUBROUTINE dump_array_to_netcdf_2d(zfilename, p_array)
    CHARACTER (len=*), INTENT(IN) :: zfilename
    REAL(wp)         , INTENT(IN) :: p_array(:,:)
    ! local variables
    INTEGER, PARAMETER :: ndims = 2
    INTEGER :: idim, ncfile, ncid_var, &
         &     ncid_dim(ndims), icount(ndims)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_util_debug:dump_array_to_netcdf_2d'


    IF (.NOT. ldebug_enable) RETURN

#ifndef DISABLE_DUMP
    WRITE (*,*) "Dumping ", zfilename
    ! create NetCDF file:
    CALL nf(nf90_create("00_"//TRIM(zfilename)//"_"//TRIM(int2string(debug_step))//".nc", &
      &               nf90_clobber, ncfile), routine)
    ! create dimensions:
    DO idim=1,ndims
      CALL nf(nf90_def_dim(ncfile, 'dim'//int2string(idim), SIZE(p_array,idim), &
        &     ncid_dim(idim)), routine)
      icount(idim) = SIZE(p_array,idim)
    END DO
    ! create variable:
    CALL nf(nf90_def_var(ncfile, "var", NF90_DOUBLE, ncid_dim, ncid_var), &
      &     routine)
    ! End of definition mode
    CALL nf(nf90_enddef(ncfile), routine)
    ! put data:
    CALL nf(nf90_put_var(ncfile, ncid_var, p_array, &
      &                  (/ (1, idim=1,ndims) /), icount), &
      &     routine)
    ! close file
    CALL nf(nf90_close(ncfile), routine)
#endif

  END SUBROUTINE dump_array_to_netcdf_2d


  !> Dumps real array to NetCDF file
  !
  !  Note: When using this SR, take care of OpenMP regions!

  SUBROUTINE dump_array_to_netcdf_3d(zfilename, p_array)
    CHARACTER (len=*), INTENT(IN) :: zfilename
    REAL(wp)         , INTENT(IN) :: p_array(:,:,:)
    ! local variables
    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: idim, ncfile, ncid_var, &
         &     ncid_dim(ndims), icount(ndims)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
      routine = 'mo_util_debug:dump_array_to_netcdf_3d'


    IF (.NOT. ldebug_enable) RETURN

#ifndef DISABLE_DUMP
    WRITE (*,*) "Dumping ", zfilename   
    ! create NetCDF file:
    CALL nf(nf90_create("00_"//TRIM(zfilename)//"_"//TRIM(int2string(debug_step))//".nc", &
      &               nf90_clobber, ncfile), routine)
    ! create dimensions:
    DO idim=1,ndims
      CALL nf(nf90_def_dim(ncfile, 'dim'//int2string(idim), SIZE(p_array,idim), &
        &     ncid_dim(idim)), routine)
      icount(idim) = SIZE(p_array,idim)
    END DO
    ! create variable:
    CALL nf(nf90_def_var(ncfile, "var", NF90_DOUBLE, ncid_dim, ncid_var), &
      &     routine)
    ! End of definition mode
    CALL nf(nf90_enddef(ncfile), routine)
    ! put data:
    CALL nf(nf90_put_var(ncfile, ncid_var, p_array, &
      &                  (/ (1, idim=1,ndims) /), icount), &
      &     routine)
    ! close file
    CALL nf(nf90_close(ncfile), routine)
#endif

  END SUBROUTINE dump_array_to_netcdf_3d

END MODULE mo_util_debug
