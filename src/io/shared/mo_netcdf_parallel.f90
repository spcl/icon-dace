! This module provides wrappers for netcdf functions.
!
! This module provides wrappers for netcdf functions
! for reading a NetCDF file in a parallel run.
! These wrappers have the same interface as the corresponding
! NetCDF routines, but only one processor actually reads the file
! and broadcasts the data to the others.
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

MODULE mo_netcdf_parallel
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind, ONLY: dp
USE mo_mpi,  ONLY: p_pe_work, p_io, p_bcast, p_comm_work

USE mo_netcdf

IMPLICIT NONE

PRIVATE

INTEGER, PARAMETER :: nf_read = NF_NOWRITE

!modules interface-------------------------------------------
!subroutines
PUBLIC :: p_nf_open
PUBLIC :: p_nf_close
PUBLIC :: p_nf_inq_dimid
PUBLIC :: p_nf_inq_dimlen
PUBLIC :: p_nf_inq_varid
PUBLIC :: p_nf_get_att_text
PUBLIC :: p_nf_get_att_int
PUBLIC :: p_nf_get_att_double
PUBLIC :: p_nf_get_var_int
PUBLIC :: p_nf_get_vara_int
PUBLIC :: p_nf_get_vara_text
PUBLIC :: p_nf_get_var_double
PUBLIC :: p_nf_get_vara_double
PUBLIC :: p_nf_get_vara_double_
PUBLIC :: p_nf_inq_attid

! constants
PUBLIC :: nf_read

! make some names from mo_netcdf also global
PUBLIC :: nf_nowrite, nf_global, nf_noerr, nf_strerror

INTERFACE p_nf_get_att_int
   MODULE PROCEDURE p_nf_get_att_int_0
   MODULE PROCEDURE p_nf_get_att_int_1
END INTERFACE

INTERFACE p_nf_get_att_double
   MODULE PROCEDURE p_nf_get_att_double_single
   MODULE PROCEDURE p_nf_get_att_double_array
END INTERFACE

!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!! Wrapper for nf_open.
!!
INTEGER FUNCTION p_nf_open(path, omode, ncid)

   CHARACTER(len=*), INTENT(in) :: path
   INTEGER, INTENT(in) :: omode
   INTEGER, INTENT(out) :: ncid

   INTEGER :: res

!-----------------------------------------------------------------------

   IF (p_pe_work == p_io) THEN
      res = nf_open(path, omode, ncid)
   ELSE
      ncid = -1 ! set it to an invalid value
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_open = res

END FUNCTION p_nf_open

!-------------------------------------------------------------------------
!! Wrapper for nf_close.
!!
INTEGER FUNCTION p_nf_close(ncid)

!
   INTEGER, INTENT(in) :: ncid

   INTEGER :: res

   IF (p_pe_work == p_io) THEN
      res = nf_close(ncid)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_close = res

END FUNCTION p_nf_close

!-------------------------------------------------------------------------
!
!

!! Wrapper for nf_inq_dimid.
!!
INTEGER FUNCTION p_nf_inq_dimid(ncid, name, dimid)

!
   INTEGER, INTENT(in) :: ncid
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: dimid

   INTEGER :: res

   IF (p_pe_work == p_io) THEN
      res = nf_inq_dimid(ncid, name, dimid)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_inq_dimid = res

   CALL p_bcast(dimid, p_io, p_comm_work)

END FUNCTION p_nf_inq_dimid

!-------------------------------------------------------------------------

!! Wrapper for nf_inq_dimlen.
!!
INTEGER FUNCTION p_nf_inq_dimlen(ncid, dimid, len)

!
   INTEGER, INTENT(in) :: ncid, dimid
   INTEGER, INTENT(out) :: len

   INTEGER :: res


!-----------------------------------------------------------------------

   IF (p_pe_work == p_io) THEN
      res = nf_inq_dimlen(ncid, dimid, len)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_inq_dimlen = res

   CALL p_bcast(len, p_io, p_comm_work)

END FUNCTION p_nf_inq_dimlen

!-------------------------------------------------------------------------
!! Wrapper for nf_inq_varid.
!!
INTEGER FUNCTION p_nf_inq_varid(ncid, name, varid)

!
   INTEGER, INTENT(in) :: ncid
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: varid

   INTEGER :: res


   IF (p_pe_work == p_io) THEN
      res = nf_inq_varid(ncid, name, varid)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_inq_varid = res

   CALL p_bcast(varid, p_io, p_comm_work)

END FUNCTION p_nf_inq_varid

!-------------------------------------------------------------------------
!! Wrapper for nf_get_att_text.
!!
INTEGER FUNCTION p_nf_get_att_text(ncid, varid, name, tval)
  INTEGER,          INTENT(in)  :: ncid, varid
  CHARACTER(len=*), INTENT(in)  :: name
  CHARACTER(len=*), INTENT(out) :: tval

  INTEGER :: res, tlen

  IF  (p_pe_work == p_io) THEN
    res = nf_inq_attlen(ncid, varid, name, tlen)
    IF (res == nf_noerr) THEN
      res = nf_get_att_text(ncid, varid, name, tval)
      tval = tval(1:tlen)
    END IF
  ENDIF

  CALL p_bcast(res, p_io, p_comm_work)
  p_nf_get_att_text = res

  CALL p_bcast(tval, p_io, p_comm_work)

END FUNCTION p_nf_get_att_text
!-----------------------------------------------------------

!-----------------------------------------------------------
!! Wrapper for nf_get_att_double
!!
INTEGER FUNCTION p_nf_get_att_double_single(ncid, varid, name, dvalue)

!
   INTEGER, INTENT(in) :: ncid, varid
   CHARACTER(len=*), INTENT(in) :: name
   REAL(dp), INTENT(out) :: dvalue

   INTEGER :: res

   IF (p_pe_work == p_io) THEN
      res = nfx_get_att(ncid, varid, name, dvalue)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_att_double_single = res

   CALL p_bcast(dvalue, p_io, p_comm_work)

END FUNCTION p_nf_get_att_double_single
!-----------------------------------------------------------------------

!-----------------------------------------------------------
!! Wrapper for nf_get_att_double
!!
INTEGER FUNCTION p_nf_get_att_double_array(ncid, varid, name, dvalue)

!
   INTEGER, INTENT(in) :: ncid, varid
   CHARACTER(len=*), INTENT(in) :: name
   REAL(dp), INTENT(out) :: dvalue(:)

   INTEGER :: res

   IF (p_pe_work == p_io) THEN
      res = nf_get_att_double(ncid, varid, name, dvalue)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_att_double_array = res

   CALL p_bcast(dvalue, p_io, p_comm_work)

END FUNCTION p_nf_get_att_double_array
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! Wrapper for nf_get_att_int.
!!
INTEGER FUNCTION p_nf_inq_attid(ncid, varid, name, ivals)

   INTEGER, INTENT(in) :: ncid, varid
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(inout) :: ivals

   INTEGER :: res


!-----------------------------------------------------------------------

   IF (p_pe_work == p_io) THEN
      res = nf_inq_attid(ncid, varid, name, ivals)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_inq_attid = res

END FUNCTION p_nf_inq_attid

!-----------------------------------------------------------------------
!! Wrapper for nf_get_att_int.
!!
INTEGER FUNCTION p_nf_get_att_int_0(ncid, varid, name, ivals)

   INTEGER, INTENT(in) :: ncid, varid
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: ivals

   INTEGER :: res


!-----------------------------------------------------------------------

   IF (p_pe_work == p_io) THEN
      res = nfx_get_att(ncid, varid, name, ivals)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_att_int_0 = res

   CALL p_bcast(ivals, p_io, p_comm_work)

END FUNCTION p_nf_get_att_int_0


!-------------------------------------------------------------------------
!! Wrapper for nf_get_att_int.
!!
INTEGER FUNCTION p_nf_get_att_int_1(ncid, varid, name, ivals)

!
   INTEGER, INTENT(in) :: ncid, varid
   CHARACTER(len=*), INTENT(in) :: name
   INTEGER, INTENT(out) :: ivals(:)

   INTEGER :: res, len


!-----------------------------------------------------------------------

   IF (p_pe_work == p_io) THEN
      ! First get the length of the attribute
      res = nf_inq_attlen (ncid, varid, name, len)
      IF(res == nf_noerr) &
         res = nf_get_att_int(ncid, varid, name, ivals)
   ENDIF

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_att_int_1 = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_work)
   CALL p_bcast(ivals(1:len), p_io, p_comm_work)

END FUNCTION p_nf_get_att_int_1

!-------------------------------------------------------------------------
!! Wrapper for nf_get_var_int.
!!
INTEGER FUNCTION p_nf_get_var_int(ncid, varid, ivals)

   INTEGER, INTENT(in) :: ncid, varid
   INTEGER, INTENT(out) :: ivals(*)

   INTEGER :: res, len, ndims, dimids(NF_MAX_VAR_DIMS), dimlen, i


   IF (p_pe_work == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      len = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen)
         IF(res /= nf_noerr) GOTO 9999
         len = len * dimlen
      ENDDO

      res = nf_get_var_int(ncid, varid, ivals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_var_int = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_work)
   CALL p_bcast(ivals(1:len), p_io, p_comm_work)

END FUNCTION p_nf_get_var_int

!-------------------------------------------------------------------------
!! Wrapper for nf_get_var_double.
!!
INTEGER FUNCTION p_nf_get_var_double(ncid, varid, dvals)

!
   INTEGER,  INTENT(in)  :: ncid, varid
   REAL(dp), INTENT(out) :: dvals(*)

   INTEGER :: res, len, ndims, dimids(NF_MAX_VAR_DIMS), dimlen, i

   IF (p_pe_work == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      len = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen)
         IF(res /= nf_noerr) GOTO 9999
         len = len * dimlen
      ENDDO

      res = nf_get_var_double(ncid, varid, dvals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_var_double = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_work)
   CALL p_bcast(dvals(1:len), p_io, p_comm_work)

END FUNCTION p_nf_get_var_double

!-------------------------------------------------------------------------
!! Wrapper for nf_get_vara_text.
!!
INTEGER FUNCTION p_nf_get_vara_text(ncid, varid, start, count, vals)

   INTEGER, INTENT(in)  :: ncid, varid, start(*), count(*)
   CHARACTER(len=*), INTENT(out) :: vals(:)

   INTEGER :: i, o, s, res, ndims, dimids(NF_MAX_VAR_DIMS), &
              start_(7), count_(7), dimlen(7)
   CHARACTER(len=size(vals,1) * len(vals(1))) :: hold

   IF (p_pe_work == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      dimlen = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen(i))
         IF(res /= nf_noerr) GOTO 9999
      ENDDO
   ENDIF

   IF (p_pe_work == p_io) THEN
      res = nf_get_var_text(ncid, varid, hold)
   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_vara_text = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) THEN
      RETURN
   END IF

   ! Broadcast number of values and values themselves

   CALL p_bcast(dimlen, p_io, p_comm_work)
   CALL p_bcast(ndims, p_io, p_comm_work)
   CALL p_bcast(hold, p_io, p_comm_work)

   o = 1
   s = len(vals(1))
   DO i = 1, size(vals,1)
     vals(i) = hold(o:o+s-1)
     o = o + s
   ENDDO

END FUNCTION p_nf_get_vara_text


!-------------------------------------------------------------------------
!! Wrapper for nf_get_vara_int.
!!
INTEGER FUNCTION p_nf_get_vara_int(ncid, varid, start, count, ivals)

!
   INTEGER, INTENT(in)  :: ncid, varid, start(*), count(*)
   INTEGER, INTENT(out) :: ivals(*)

   INTEGER :: i, res, ndims, dimids(NF_MAX_VAR_DIMS), &
              start_(7), count_(7), dimlen(7), len
   INTEGER, ALLOCATABLE :: t_ivals(:,:,:,:,:,:,:)

   IF (p_pe_work == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      dimlen = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen(i))
         IF(res /= nf_noerr) GOTO 9999
      ENDDO

      ALLOCATE(t_ivals(dimlen(1), dimlen(2), dimlen(3), dimlen(4), &
                       dimlen(5), dimlen(6), dimlen(7)))

      res = nf_get_var_int(ncid, varid, t_ivals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_vara_int = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) THEN
      IF ((p_pe_work == p_io) .AND. (ALLOCATED(t_ivals))) THEN
         DEALLOCATE(t_ivals)
      END IF
      return
   END IF

   ! Broadcast number of values and values themselves

   CALL p_bcast(dimlen, p_io, p_comm_work)
   CALL p_bcast(ndims, p_io, p_comm_work)

   IF (p_pe_work /= p_io) &
      ALLOCATE(t_ivals(dimlen(1), dimlen(2), dimlen(3), dimlen(4), &
                       dimlen(5), dimlen(6), dimlen(7)))

   CALL p_bcast(t_ivals, p_io, p_comm_work)

   start_ = 1
   count_ = 1
   start_(1:ndims) = start(1:ndims)
   count_(1:ndims) = count(1:ndims)

   len = 1
   DO i = 1, ndims
      len = len * count(i)
   END DO

   ivals(1:len) = RESHAPE(t_ivals(start_(1):start_(1)+count_(1) - 1, &
                                  start_(2):start_(2)+count_(2) - 1, &
                                  start_(3):start_(3)+count_(3) - 1, &
                                  start_(4):start_(4)+count_(4) - 1, &
                                  start_(5):start_(5)+count_(5) - 1, &
                                  start_(6):start_(6)+count_(6) - 1, &
                                  start_(7):start_(7)+count_(7) - 1), (/len/))

   DEALLOCATE(t_ivals)

END FUNCTION p_nf_get_vara_int


!-------------------------------------------------------------------------
!! Wrapper for nf_get_vara_double.
!!
INTEGER FUNCTION p_nf_get_vara_double_(ncid, varid, start, count, dvals)

!
   INTEGER, INTENT(in)  :: ncid, varid, start(*), count(*)
   REAL(dp), INTENT(out) :: dvals(*)

   INTEGER :: i, res, ndims, dimids(NF_MAX_VAR_DIMS), &
              start_(7), count_(7), dimlen(7), len
   REAL(dp), ALLOCATABLE :: t_dvals(:,:,:,:,:,:,:)

   IF (p_pe_work == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      dimlen = 1
      DO i = 1, ndims
         res = nf_inq_dimlen(ncid, dimids(i), dimlen(i))
         IF(res /= nf_noerr) GOTO 9999
      ENDDO

      ALLOCATE(t_dvals(dimlen(1), dimlen(2), dimlen(3), dimlen(4), &
                       dimlen(5), dimlen(6), dimlen(7)))

      res = nf_get_var_double(ncid, varid, t_dvals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_vara_double_ = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) THEN
      IF ((p_pe_work == p_io) .AND. (ALLOCATED(t_dvals))) THEN
         DEALLOCATE(t_dvals)
      END IF
      return
   END IF

   ! Broadcast number of values and values themselves

   CALL p_bcast(dimlen, p_io, p_comm_work)
   CALL p_bcast(ndims, p_io, p_comm_work)

   IF (p_pe_work /= p_io) &
      ALLOCATE(t_dvals(dimlen(1), dimlen(2), dimlen(3), dimlen(4), &
                       dimlen(5), dimlen(6), dimlen(7)))

   CALL p_bcast(t_dvals, p_io, p_comm_work)

   start_ = 1
   count_ = 1
   start_(1:ndims) = start(1:ndims)
   count_(1:ndims) = count(1:ndims)

   len = 1
   DO i = 1, ndims
      len = len * count(i)
   END DO

   dvals(1:len) = RESHAPE(t_dvals(start_(1):start_(1)+count_(1) - 1, &
                                  start_(2):start_(2)+count_(2) - 1, &
                                  start_(3):start_(3)+count_(3) - 1, &
                                  start_(4):start_(4)+count_(4) - 1, &
                                  start_(5):start_(5)+count_(5) - 1, &
                                  start_(6):start_(6)+count_(6) - 1, &
                                  start_(7):start_(7)+count_(7) - 1), (/len/))

   DEALLOCATE(t_dvals)

END FUNCTION p_nf_get_vara_double_

!-------------------------------------------------------------------------
!! Wrapper for nf_get_vara_double.
!!
INTEGER FUNCTION p_nf_get_vara_double(ncid, varid, start, count, dvals)

!
   INTEGER,  INTENT(in)  :: ncid, varid, start(*), count(*)
   REAL(dp), INTENT(out) :: dvals(*)

   INTEGER :: res, len, ndims, dimids(NF_MAX_VAR_DIMS), i


   IF (p_pe_work == p_io) THEN

      ! First get the length of the array

      res = nf_inq_varndims(ncid, varid, ndims)
      IF(res /= nf_noerr) GOTO 9999
      res = nf_inq_vardimid(ncid, varid, dimids)
      IF(res /= nf_noerr) GOTO 9999

      len = 1
      DO i = 1, ndims
         len = len * count(i)
      ENDDO

      res = nf_get_vara_double(ncid, varid, start, count, dvals)

   ENDIF

9999 CONTINUE

   CALL p_bcast(res, p_io, p_comm_work)
   p_nf_get_vara_double = res

   ! If there was an error, don't try to broadcast the values

   IF(res /= nf_noerr) return

   ! Broadcast number of values and values themselves

   CALL p_bcast(len, p_io, p_comm_work)
   CALL p_bcast(dvals(1:len), p_io, p_comm_work)

END FUNCTION p_nf_get_vara_double

!-------------------------------------------------------------------------

END MODULE mo_netcdf_parallel
