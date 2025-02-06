!> @file netcdf_utils.F90
!! @brief Utilities for using the netcdf library
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

MODULE netcdf_utils

  USE netcdf,                 ONLY: nf90_def_dim, nf90_def_var, nf90_inq_varid, &
       &                            nf90_inquire_variable, nf90_inquire_dimension, &
       &                            NF90_NOERR, nf90_strerror
  USE utils,                  ONLY: int2string
  USE comin_plugin_interface, ONLY: EP_FINISH, comin_current_get_ep, &
                                    comin_plugin_finish

  IMPLICIT NONE

  PUBLIC nf90_utils_def_var, nf90_utils_get_shape, nf90

CONTAINS

  FUNCTION nf90_utils_def_var(ncid, name, typ, dims, opt_prepend_dimids) &
    RESULT(var_id)
    INTEGER, INTENT(IN)           :: ncid
    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER, INTENT(IN)           :: typ
    INTEGER, INTENT(IN)           :: dims(:)
    INTEGER, INTENT(IN), OPTIONAL :: opt_prepend_dimids(:)

    INTEGER :: i, var_id
    INTEGER, ALLOCATABLE :: dim_ids(:)

    IF (PRESENT(opt_prepend_dimids)) THEN
      ALLOCATE(dim_ids(SIZE(dims) + SIZE(opt_prepend_dimids)))
      dim_ids(1:SIZE(opt_prepend_dimids)) = opt_prepend_dimids
    ELSE
      ALLOCATE(dim_ids(SIZE(dims)))
    ENDIF

    ! create dimensions
    DO i=1,SIZE(dims)
      CALL nf90(nf90_def_dim(ncid, name//"_dim_"//int2string(i), &
                             dims(i), dim_ids(i)))
    END DO
    IF (PRESENT(opt_prepend_dimids)) THEN
      dim_ids(SIZE(dims)+1:) = opt_prepend_dimids
    ENDIF

    CALL nf90(nf90_def_var(ncid, name, typ, dim_ids, var_id))
  END FUNCTION nf90_utils_def_var

  FUNCTION nf90_utils_get_shape(ncid, name, ndims, opt_varid) &
    RESULT(dims)
    INTEGER, INTENT(IN)            :: ncid
    CHARACTER(LEN=*),  INTENT(IN)  :: name
    INTEGER,           INTENT(IN)  :: ndims
    INTEGER, OPTIONAL, INTENT(OUT) :: opt_varid
    INTEGER :: varid, dims(ndims)

    INTEGER :: dimids(ndims), i

    CALL nf90(nf90_inq_varid(ncid, name, varid))
    IF (PRESENT(opt_varid)) opt_varid = varid
    CALL nf90(nf90_inquire_variable(ncid, varid, dimids=dimids))
    ! create dimensions
    DO i=1,ndims
      CALL nf90(nf90_inquire_dimension(ncid, dimids(i), len=dims(i)))
    END DO
  END FUNCTION nf90_utils_get_shape

  RECURSIVE SUBROUTINE nf90(status)
    INTEGER, INTENT(in) :: status
    IF (status /= NF90_NOERR) THEN
      IF (comin_current_get_ep() /= EP_FINISH) THEN
        CALL comin_plugin_finish('netCDF error', nf90_strerror(status))
      ELSE
        WRITE(0,*) 'netCDF error' // nf90_strerror(status)
      ENDIF
    ENDIF
  END SUBROUTINE nf90

END MODULE netcdf_utils
