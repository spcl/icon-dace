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

! This module provides wrappers for netcdf functions for reading a NetCDF file
! in a parallel run. These wrappers have the same interface as the corresponding
! NetCDF routines, but only one processor actually reads the file and broadcasts
! the data to the others.

MODULE mo_netcdf_parallel

USE mo_kind, ONLY: dp
USE mo_mpi, ONLY: p_pe_work, p_io, p_bcast, p_comm_work

USE mo_netcdf

IMPLICIT NONE

PRIVATE

PUBLIC :: p_nf90_open
PUBLIC :: p_nf90_close
PUBLIC :: p_nf90_inq_dimid
PUBLIC :: p_nf90_inquire_dimension
PUBLIC :: p_nf90_inquire_attribute
PUBLIC :: p_nf90_inquire_variable

PUBLIC :: p_nf90_get_att
INTERFACE p_nf90_get_att
  MODULE PROCEDURE p_nf90_get_att_real_dp
  MODULE PROCEDURE p_nf90_get_att_one_real_dp
  MODULE PROCEDURE p_nf90_get_att_one_int
  MODULE PROCEDURE p_nf90_get_att_text
END INTERFACE p_nf90_get_att

PUBLIC :: p_nf90_inq_varid

PUBLIC :: p_nf90_get_var
INTERFACE p_nf90_get_var
  MODULE PROCEDURE p_nf90_get_var_1D_real_dp
  MODULE PROCEDURE p_nf90_get_var_2D_real_dp
  MODULE PROCEDURE p_nf90_get_var_3D_real_dp
  MODULE PROCEDURE p_nf90_get_var_4D_real_dp
  MODULE PROCEDURE p_nf90_get_var_1D_int
  MODULE PROCEDURE p_nf90_get_var_2D_int
  MODULE PROCEDURE p_nf90_get_var_3D_int
  MODULE PROCEDURE p_nf90_get_var_text
  MODULE PROCEDURE p_nf90_get_var_1D_text
END INTERFACE p_nf90_get_var

PUBLIC :: p_nf90x_get_var_local
INTERFACE p_nf90x_get_var_local
  MODULE PROCEDURE p_nf90x_get_var_local_1D_real_dp
  MODULE PROCEDURE p_nf90x_get_var_local_1D_int
  MODULE PROCEDURE p_nf90x_get_var_local_2D_int
END INTERFACE p_nf90x_get_var_local

REAL(dp), ALLOCATABLE :: local_buffer_real_dp(:)
INTEGER, ALLOCATABLE :: local_buffer_int(:)

CONTAINS

INTEGER FUNCTION p_nf90_open(path, omode, ncid)

  CHARACTER(len=*), INTENT(in) :: path
  INTEGER, INTENT(in) :: omode
  INTEGER, INTENT(out) :: ncid

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
      res = nf90_open(path, omode, ncid)
  ELSE
      ncid = -1 ! set it to an invalid value
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_open = res

END FUNCTION p_nf90_open

INTEGER FUNCTION p_nf90_close(ncid)

  INTEGER, INTENT(in) :: ncid

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_close(ncid)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_close = res

  IF (ALLOCATED(local_buffer_real_dp)) DEALLOCATE(local_buffer_real_dp)
  IF (ALLOCATED(local_buffer_int)) DEALLOCATE(local_buffer_int)

END FUNCTION p_nf90_close

INTEGER FUNCTION p_nf90_inq_dimid(ncid, name, dimid)

  INTEGER, INTENT(in) :: ncid
  CHARACTER(len=*), INTENT(in) :: name
  INTEGER, INTENT(out) :: dimid

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_inq_dimid(ncid, name, dimid)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_inq_dimid = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(dimid, p_io, p_comm_work)

END FUNCTION p_nf90_inq_dimid

INTEGER FUNCTION p_nf90_inquire_dimension(ncid, dimid, name, len)

  INTEGER, INTENT(in) :: ncid, dimid
  CHARACTER(len=*), OPTIONAL, INTENT(out) :: name
  INTEGER, OPTIONAL, INTENT(out) :: len

  INTEGER :: res

  CHARACTER(len=NF90_MAX_NAME) :: local_name
  INTEGER :: local_len

  IF (p_pe_work == p_io) THEN
    res = nf90_inquire_dimension(ncid, dimid, &
                               & name = local_name, &
                               & len = local_len)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_inquire_dimension = res

  IF (res /= NF90_NOERR) RETURN

  IF (PRESENT(name)) THEN
    CALL p_bcast(local_name, p_io, p_comm_work)
    name = local_name
  ENDIF
  
  IF (PRESENT(len)) THEN
    CALL p_bcast(local_len, p_io, p_comm_work)
    len = local_len
  ENDIF

END FUNCTION p_nf90_inquire_dimension

INTEGER FUNCTION p_nf90_inquire_attribute(ncid, varid, name, xtype, len, attnum)

  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), INTENT(in) :: name
  INTEGER, OPTIONAL, INTENT(out) :: xtype, len, attnum

  INTEGER :: res

  ! Pack local values into a single array to broadcast them in one go:
  INTEGER :: local_values(3)  ! xtype, len, attnum

  IF (p_pe_work == p_io) THEN
    res = nf90_inquire_attribute(ncid, varid, name, &
                               & xtype = local_values(1), &
                               & len = local_values(2), &
                               & attnum = local_values(3))
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_inquire_attribute = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(local_values, p_io, p_comm_work)

  IF (PRESENT(xtype)) THEN
    xtype = local_values(1)
  ENDIF

  IF (PRESENT(len)) THEN
    len = local_values(2)
  ENDIF

  IF (PRESENT(attnum)) THEN
    attnum = local_values(3)
  ENDIF

END FUNCTION p_nf90_inquire_attribute

INTEGER FUNCTION p_nf90_inquire_variable(ncid, varid, name, xtype, ndims, dimids)

  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), OPTIONAL, INTENT(out) :: name
  INTEGER, OPTIONAL, INTENT(out) :: xtype, ndims, dimids(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  CHARACTER(len=NF90_MAX_NAME) :: local_name
  ! Pack local values into a single array to broadcast them in one go:
  INTEGER :: local_values(1 + 1 + NF90_MAX_VAR_DIMS)  ! xtype, ndims, dimids

  IF (p_pe_work == p_io) THEN
    res = nf90_inquire_variable(ncid, varid, local_name, &
                               & xtype = local_values(1), &
                               & ndims = local_values(2), &
                               & dimids = local_values(3:))
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_inquire_variable = res

  IF (res /= NF90_NOERR) RETURN

  IF (PRESENT(name)) THEN
    CALL p_bcast(local_name, p_io, p_comm_work)
    name = local_name
  ENDIF

  IF (PRESENT(xtype) .OR. PRESENT(ndims) .OR. PRESENT(dimids)) THEN
    CALL p_bcast(local_values, p_io, p_comm_work)
  ENDIF

  IF (PRESENT(xtype)) THEN
    xtype = local_values(1)
  ENDIF

  IF (PRESENT(ndims)) THEN
    ndims = local_values(2)
  ENDIF

  IF (PRESENT(dimids)) THEN
    IF (SIZE(dimids) .GE. local_values(2)) THEN
      dimids(:local_values(2)) = local_values(3:3 + local_values(2) - 1)
    ELSE
      p_nf90_inquire_variable = NF90_EINVAL
    ENDIF
  ENDIF

END FUNCTION p_nf90_inquire_variable

INTEGER FUNCTION p_nf90_get_att_real_dp(ncid, varid, name, values)

  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), INTENT(in) :: name
  REAL(dp), INTENT(out) :: values(:)

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_get_att(ncid, varid, name, values)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_att_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_att_real_dp

INTEGER FUNCTION p_nf90_get_att_one_real_dp(ncid, varid, name, values)

  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), INTENT(in) :: name
  REAL(dp), INTENT(out) :: values

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_get_att(ncid, varid, name, values)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_att_one_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_att_one_real_dp

INTEGER FUNCTION p_nf90_get_att_one_int(ncid, varid, name, values)

  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), INTENT(in) :: name
  INTEGER, INTENT(out) :: values

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_get_att(ncid, varid, name, values)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_att_one_int = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_att_one_int

INTEGER FUNCTION p_nf90_get_att_text(ncid, varid, name, values)

  INTEGER, INTENT(in) :: ncid, varid
  CHARACTER(len=*), INTENT(in) :: name
  CHARACTER(len=*), INTENT(out) :: values

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_get_att(ncid, varid, name, values)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_att_text = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_att_text

INTEGER FUNCTION p_nf90_inq_varid(ncid, name, varid)

  INTEGER, INTENT(in) :: ncid
  CHARACTER(len=*), INTENT(in) :: name
  INTEGER, INTENT(out) :: varid

  INTEGER :: res

  IF (p_pe_work == p_io) THEN
    res = nf90_inq_varid(ncid, name, varid)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_inq_varid = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(varid, p_io, p_comm_work)

END FUNCTION p_nf90_inq_varid

INTEGER FUNCTION p_nf90_get_var_1D_real_dp(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  REAL(dp), INTENT(out) :: values(:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_1D_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_1D_real_dp

INTEGER FUNCTION p_nf90_get_var_2D_real_dp(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  REAL(dp), INTENT(out) :: values(:,:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_2D_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_2D_real_dp

INTEGER FUNCTION p_nf90_get_var_3D_real_dp(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  REAL(dp), INTENT(out) :: values(:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_3D_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_3D_real_dp

INTEGER FUNCTION p_nf90_get_var_4D_real_dp(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  REAL(dp), INTENT(out) :: values(:,:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_4D_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_4D_real_dp

INTEGER FUNCTION p_nf90_get_var_1D_int(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  INTEGER, INTENT(out) :: values(:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_1D_int = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_1D_int

INTEGER FUNCTION p_nf90_get_var_2D_int(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  INTEGER, INTENT(out) :: values(:,:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_2D_int = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_2D_int

INTEGER FUNCTION p_nf90_get_var_3D_int(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  INTEGER, INTENT(out) :: values(:,:,:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER :: ndims

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    ndims = SIZE(SHAPE(values))
    local_start(:) = 1
    local_count(:ndims) = SHAPE(values)
    local_count(ndims + 1:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_3D_int = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_3D_int

INTEGER FUNCTION p_nf90_get_var_text(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  CHARACTER(len=*), INTENT(out) :: values
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    local_start(:) = 1
    local_count(1) = LEN(values)
    local_count(2:) = 1

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_text = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_text

INTEGER FUNCTION p_nf90_get_var_1D_text(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  CHARACTER(len=*), INTENT(out) :: values(:)
  INTEGER, OPTIONAL, INTENT(in) :: start(:), count(:)
  ! TODO: add the remaining OPTIONAL arguments of the original interface

  INTEGER :: res

  INTEGER, DIMENSION(NF90_MAX_VAR_DIMS) :: local_start, local_count
  INTEGER, PARAMETER :: ndims = 1

  IF (p_pe_work == p_io) THEN
    ! As in the original implementation, the user is responsible for all size
    ! mismatches (i.e. we do not query the file for the real number and the
    ! sizes of the dimenstions:
    local_start(:) = 1
    local_count(:ndims + 1) = [LEN(values(1)), SHAPE(values)]
    local_count(ndims + 2:) = 0

    IF (PRESENT(start)) local_start(:SIZE(start)) = start(:)
    IF (PRESENT(count)) local_count(:SIZE(count)) = count(:)

    res = nf90_get_var(ncid, varid, values, &
                     & start = local_start, count = local_count)
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90_get_var_1D_text = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(values, p_io, p_comm_work)

END FUNCTION p_nf90_get_var_1D_text

SUBROUTINE ensure_buffer_real_dp(buffer_size)

  INTEGER, INTENT(in) :: buffer_size

  IF (ALLOCATED(local_buffer_real_dp)) THEN
    IF (SIZE(local_buffer_real_dp) < buffer_size) &
    & DEALLOCATE(local_buffer_real_dp)
  ENDIF

  IF (.NOT. ALLOCATED(local_buffer_real_dp)) &
  & ALLOCATE(local_buffer_real_dp(buffer_size))

END SUBROUTINE ensure_buffer_real_dp

SUBROUTINE ensure_buffer_int(buffer_size)

  INTEGER, INTENT(in) :: buffer_size

  IF (ALLOCATED(local_buffer_int)) THEN
    IF (SIZE(local_buffer_int) < buffer_size) &
    & DEALLOCATE(local_buffer_int)
  ENDIF

  IF (.NOT. ALLOCATED(local_buffer_int)) &
  & ALLOCATE(local_buffer_int(buffer_size))

END SUBROUTINE ensure_buffer_int

SUBROUTINE get_slices_real_dp(tgt, src, dimlens, start, count)

  REAL(dp), INTENT(out) :: tgt(*)
  REAL(dp), INTENT(in) :: src(:)
  INTEGER, DIMENSION(:), INTENT(in) :: dimlens, start, count

  INTEGER, DIMENSION(SIZE(dimlens)) :: offsets, dimsizes
  INTEGER :: ii, jj

  offsets(:) = start(:SIZE(dimlens))
  dimsizes(1) = 1
  dimsizes(2:) = [(PRODUCT(dimlens(:jj)), jj = 1, SIZE(dimlens) - 1)]

  DO ii = 1, PRODUCT(count)
    DO jj = 1, SIZE(dimlens) - 1
      IF (offsets(jj) > start(jj) + count(jj) - 1) THEN
        offsets(jj) = start(jj)
        offsets(jj + 1) = offsets(jj + 1) + 1
      ELSE
        EXIT
      ENDIF
    ENDDO
    tgt(ii) = src(SUM(dimsizes * (offsets - 1)) + 1)
    offsets(1) = offsets(1) + 1
  ENDDO

END SUBROUTINE get_slices_real_dp

SUBROUTINE get_slices_int(tgt, src, dimlens, start, count)

  INTEGER, INTENT(out) :: tgt(*)
  INTEGER, INTENT(in) :: src(:)
  INTEGER, DIMENSION(:), INTENT(in) :: dimlens, start, count

  INTEGER, DIMENSION(SIZE(dimlens)) :: offsets, dimsizes
  INTEGER :: ii, jj

  offsets(:) = start(:SIZE(dimlens))
  dimsizes(1) = 1
  dimsizes(2:) = [(PRODUCT(dimlens(:jj)), jj = 1, SIZE(dimlens) - 1)]

  DO ii = 1, PRODUCT(count)
    DO jj = 1, SIZE(dimlens) - 1
      IF (offsets(jj) > start(jj) + count(jj) - 1) THEN
        offsets(jj) = start(jj)
        offsets(jj + 1) = offsets(jj + 1) + 1
      ELSE
        EXIT
      ENDIF
    ENDDO
    tgt(ii) = src(SUM(dimsizes * (offsets - 1)) + 1)
    offsets(1) = offsets(1) + 1
  ENDDO

END SUBROUTINE get_slices_int

INTEGER FUNCTION p_nf90x_get_var_local_assumed_real_dp(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  REAL(dp), INTENT(out) :: values(*)
  INTEGER, INTENT(in) :: start(:), count(:)

  INTEGER :: res

  INTEGER :: ii, dimdata(NF90_MAX_VAR_DIMS + 1)  ! ndims + dimids/dimlens

  IF (p_pe_work == p_io) THEN
    res = nf90_inquire_variable(ncid, varid, ndims = dimdata(1), dimids = dimdata(2:))
    IF (res /= NF90_NOERR) GOTO 999

    DO ii = 1, dimdata(1)
      res = nf90_inquire_dimension(ncid, dimdata(1 + ii), len = dimdata(1 + ii))
      IF (res /= NF90_NOERR) GOTO 999
    ENDDO
    dimdata(1 + dimdata(1) + 1:) = 1

    CALL ensure_buffer_real_dp(PRODUCT(dimdata(2:)))
    res = nf90_get_var(ncid, varid, local_buffer_real_dp, count = dimdata(2:))
  ENDIF

999 CONTINUE

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90x_get_var_local_assumed_real_dp = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(dimdata, p_io, p_comm_work)
  IF (p_pe_work /= p_io) THEN
    CALL ensure_buffer_real_dp(PRODUCT(dimdata(2:)))
  ENDIF

  CALL p_bcast(local_buffer_real_dp, p_io, p_comm_work)

  CALL get_slices_real_dp(values, local_buffer_real_dp, &
                        & dimdata(2:dimdata(1) + 1), start, count)

END FUNCTION p_nf90x_get_var_local_assumed_real_dp

INTEGER FUNCTION p_nf90x_get_var_local_assumed_int(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  INTEGER, INTENT(out) :: values(*)
  INTEGER, INTENT(in) :: start(:), count(:)

  INTEGER :: res

  INTEGER :: ii, dimdata(NF90_MAX_VAR_DIMS + 1)  ! ndims + dimids/dimlens

  IF (p_pe_work == p_io) THEN
    res = nf90_inquire_variable(ncid, varid, ndims = dimdata(1), dimids = dimdata(2:))
    IF (res /= NF90_NOERR) GOTO 999

    DO ii = 1, dimdata(1)
      res = nf90_inquire_dimension(ncid, dimdata(1 + ii), len = dimdata(1 + ii))
      IF (res /= NF90_NOERR) GOTO 999
    ENDDO
    dimdata(1 + dimdata(1) + 1:) = 1

    CALL ensure_buffer_int(PRODUCT(dimdata(2:)))
    res = nf90_get_var(ncid, varid, local_buffer_int, count = dimdata(2:))
  ENDIF

999 CONTINUE

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf90x_get_var_local_assumed_int = res

  IF (res /= NF90_NOERR) RETURN

  CALL p_bcast(dimdata, p_io, p_comm_work)
  IF (p_pe_work /= p_io) THEN
    CALL ensure_buffer_int(PRODUCT(dimdata(2:)))
  ENDIF

  CALL p_bcast(local_buffer_int, p_io, p_comm_work)

  CALL get_slices_int(values, local_buffer_int, &
                    & dimdata(2:dimdata(1) + 1), start, count)

END FUNCTION p_nf90x_get_var_local_assumed_int

INTEGER FUNCTION p_nf90x_get_var_local_1D_real_dp(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  REAL(dp), INTENT(out) :: values(:)
  INTEGER, INTENT(in) :: start(:), count(:)

  p_nf90x_get_var_local_1D_real_dp = &
  & p_nf90x_get_var_local_assumed_real_dp(ncid, varid, values, start, count)

END FUNCTION p_nf90x_get_var_local_1D_real_dp

INTEGER FUNCTION p_nf90x_get_var_local_1D_int(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  INTEGER, INTENT(out) :: values(:)
  INTEGER, INTENT(in) :: start(:), count(:)

  p_nf90x_get_var_local_1D_int = &
  & p_nf90x_get_var_local_assumed_int(ncid, varid, values, start, count)

END FUNCTION p_nf90x_get_var_local_1D_int

INTEGER FUNCTION p_nf90x_get_var_local_2D_int(ncid, varid, values, start, count)

  INTEGER, INTENT(in)  :: ncid, varid
  INTEGER, INTENT(out) :: values(:,:)
  INTEGER, INTENT(in) :: start(:), count(:)

  p_nf90x_get_var_local_2D_int = &
  & p_nf90x_get_var_local_assumed_int(ncid, varid, values, start, count)

END FUNCTION p_nf90x_get_var_local_2D_int

END MODULE mo_netcdf_parallel
