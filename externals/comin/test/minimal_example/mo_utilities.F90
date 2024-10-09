!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE mo_utilities

  USE mpi
  USE netcdf
#ifndef __NO_ICON_COMIN__
  USE comin_host_interface, ONLY: comin_callback_context_call, &
    &                             EP_FINISH, COMIN_DOMAIN_OUTSIDE_LOOP
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: wp
  PUBLIC :: int2string
  PUBLIC :: blk_no
  PUBLIC :: idx_no
  PUBLIC :: finish
  PUBLIC :: t_geographical_coordinates
  PUBLIC :: read_netcdf_scalar
  PUBLIC :: filename_max
  PUBLIC :: MAX_DATETIME_STR_LEN

  INTERFACE read_netcdf_scalar
    MODULE PROCEDURE read_netcdf_scalar_r
    MODULE PROCEDURE read_netcdf_scalar_i1
    MODULE PROCEDURE read_netcdf_scalar_i2
  END INTERFACE read_netcdf_scalar

  CHARACTER(LEN=*), PARAMETER :: modname = "mo_utilities"

  ! floating-point precision
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

  ! max. file name length
  INTEGER, PARAMETER :: filename_max = 1024

  !> string length for iso string of datatime
  INTEGER, PARAMETER           :: MAX_DATETIME_STR_LEN = 32

  ! geographical coordinate class
  TYPE t_geographical_coordinates
    REAL(wp) :: lon
    REAL(wp) :: lat
  END TYPE t_geographical_coordinates

  INTEGER, PUBLIC :: mpi_rank

CONTAINS
  ! --------------------------------------------------------------------------------
  ! returns integer n as a string (often needed in printing messages)
  ! --------------------------------------------------------------------------------
  FUNCTION int2string(n, opt_fmt)
    CHARACTER(:), ALLOCATABLE :: int2string ! result
    CHARACTER(len=128) :: res
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=128) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(i0)'
    END IF
    WRITE(res,fmt) n
    res = ADJUSTL(res)
    int2string = TRIM(res)
  END FUNCTION int2string

  ELEMENTAL INTEGER FUNCTION blk_no(j, nproma)
    INTEGER, INTENT(in) :: j, nproma
    blk_no = MAX((ABS(j)-1)/nproma + 1, 1) ! i.e. also 1 for j=0, nproma=1
  END FUNCTION blk_no

  ELEMENTAL INTEGER FUNCTION idx_no(j, nproma)
    INTEGER, INTENT(in) :: j, nproma
    IF(j==0) THEN
      idx_no = 0
    ELSE
      idx_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
    ENDIF
  END FUNCTION idx_no

  SUBROUTINE finish(routine, text)
    CHARACTER(LEN=*), INTENT(IN) :: routine
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text
    INTEGER :: ierr, iexit

#ifndef __NO_ICON_COMIN__
    CALL comin_callback_context_call(EP_FINISH, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)
#endif /* #ifndef __NO_ICON_COMIN__ */

    iexit = -1
    IF (mpi_rank == 0)  WRITE (0,*) routine, ": ", text, iexit
    CALL MPI_ABORT(MPI_COMM_WORLD, iexit, ierr)
  END SUBROUTINE finish

  !----------------------------------------------------------------
  !> Auxiliary routine. Check for NetCDF error.
  !----------------------------------------------------------------
  SUBROUTINE nf(status, routine)
    INTEGER, INTENT(in) :: status
    CHARACTER(len=*), INTENT(in) :: routine
    IF (status /= nf90_noerr) THEN
      CALL finish(TRIM(routine)//' netCDF error', nf90_strerror(status))
    ENDIF
  END SUBROUTINE nf

  !----------------------------------------------------------------
  !> Auxiliary routine. Read scalar REAL array. Blocked output.
  !----------------------------------------------------------------
  SUBROUTINE read_netcdf_scalar_r(name, nproma, out_buffer, ltype, opt_size, opt_nblks, opt_npromz, &
    &                             opt_ncfileID, opt_filename)
    CHARACTER(LEN=*),           INTENT(IN)    :: name            !< variable name
    INTEGER,                    INTENT(IN)    :: nproma          !< cache line length
    TYPE(t_geographical_coordinates),  ALLOCATABLE,     INTENT(INOUT) :: out_buffer(:,:) !< blocked result array (nproma, out_nblks).
    CHARACTER(LEN=*),           INTENT(IN)    :: ltype
    INTEGER,          OPTIONAL, INTENT(OUT)   :: opt_size        !< number of items
    INTEGER,          OPTIONAL, INTENT(OUT)   :: opt_nblks       !< block number of items
    INTEGER,          OPTIONAL, INTENT(OUT)   :: opt_npromz      !< size of last block
    INTEGER,          OPTIONAL, INTENT(IN)    :: opt_ncfileID    !< NetCDF file ID (already opened).
    CHARACTER(LEN=filename_max), OPTIONAL, INTENT(IN)    :: opt_filename    !< string identifier, user-defined
    !
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::read_netcdf_scalar_r"
    REAL(wp), ALLOCATABLE :: buffer1d(:)
    INTEGER :: varID, ndims, dimIDs(1), ilength, ncfileID, nblks

    IF (PRESENT(opt_ncfileID))  ncfileID = opt_ncfileID
    IF (PRESENT(opt_filename) .AND. .NOT. PRESENT(opt_ncfileID)) &
      CALL nf(nf90_open(TRIM(opt_filename), NF90_NOWRITE, ncfileID), routine)
    ! get variable info:
    CALL nf(nf90_inq_varid(ncfileID, TRIM(name), varID), routine//","//TRIM(name))
    CALL nf(nf90_inquire_variable(ncfileID, varID, ndims=ndims), routine//","//TRIM(name))
    ! get dimension info:
    IF (ndims > 1)  CALL finish(routine, "1D array expected!")
    CALL nf(nf90_inquire_variable(ncfileID, varID, dimids=dimIDs), routine//","//TRIM(name))
    CALL nf(nf90_inquire_dimension(ncfileID, dimids(1), len=ilength), routine//","//TRIM(name))
    ! read data:
    IF (.NOT. ALLOCATED(buffer1D))  ALLOCATE(buffer1D(ilength))
    IF (PRESENT(opt_size))  opt_size  = ilength
    nblks      = blk_no(ilength,nproma)
    IF (PRESENT(opt_nblks))  opt_nblks  = nblks
    IF (PRESENT(opt_npromz)) opt_npromz = idx_no(ilength,nproma)
    IF (.NOT. ALLOCATED(out_buffer))  ALLOCATE(out_buffer(nproma,nblks))
    CALL nf(nf90_get_var(ncfileID, varID, buffer1D), routine//","//TRIM(name))

    IF (TRIM(ltype)=='lon') out_buffer(:,:)%lon = RESHAPE(buffer1D, shape=[nproma,nblks], pad=[0._wp])
    IF (TRIM(ltype)=='lat') out_buffer(:,:)%lat = RESHAPE(buffer1D, shape=[nproma,nblks], pad=[0._wp])
    IF (PRESENT(opt_filename) .AND. .NOT. PRESENT(opt_ncfileID)) &
      CALL nf(nf90_close(ncfileID), routine)
  END SUBROUTINE read_netcdf_scalar_r

  !----------------------------------------------------------------
  !> Auxiliary routine. Read scalar INTEGER array. Blocked output.
  !----------------------------------------------------------------
  SUBROUTINE read_netcdf_scalar_i1(name, nproma, out_buffer, out_size, out_nblks, out_npromz, &
    &                              opt_ncfileID, opt_filename)
    CHARACTER(LEN=*),           INTENT(IN)    :: name              !< variable name
    INTEGER,                    INTENT(IN)    :: nproma            !< cache line length
    INTEGER,  ALLOCATABLE,      INTENT(INOUT) :: out_buffer(:,:) !< blocked result array.
    INTEGER,          OPTIONAL, INTENT(OUT)   :: out_size          !< number of items
    INTEGER,          OPTIONAL, INTENT(OUT)   :: out_nblks         !< block number of items
    INTEGER,          OPTIONAL, INTENT(OUT)   :: out_npromz        !< size of last block
    INTEGER,          OPTIONAL, INTENT(IN)    :: opt_ncfileID      !< NetCDF file ID (already opened).
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: opt_filename      !< string identifier, user-defined
    !
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::read_netcdf_scalar_i1"
    INTEGER, ALLOCATABLE :: buffer1d(:)
    INTEGER :: varID, ndims, dimIDs(1), ilength(1), ncfileID, i, jl, jb, nblks

    IF (PRESENT(opt_ncfileID))  ncfileID = opt_ncfileID
    IF (PRESENT(opt_filename) .AND. .NOT. PRESENT(opt_ncfileID)) &
      CALL nf(nf90_open(TRIM(opt_filename), NF90_NOWRITE, ncfileID), routine)
    ! get variable info:
    CALL nf(nf90_inq_varid(ncfileID, TRIM(name), varID), routine//","//TRIM(name))
    CALL nf(nf90_inquire_variable(ncfileID, varID, ndims=ndims), routine//","//TRIM(name))
    ! get dimension info:
    IF (ndims /= 1)  CALL finish(routine, "1D array expected!")
    CALL nf(nf90_inquire_variable(ncfileID, varID, dimids=dimIDs), routine//","//TRIM(name))
    CALL nf(nf90_inquire_dimension(ncfileID, dimids(1), len=ilength(1)), routine//","//TRIM(name))
    ! read data:
    IF (.NOT. ALLOCATED(buffer1D)) ALLOCATE(buffer1D(ilength(1)))
    IF (PRESENT(out_size))  out_size  = ilength(1)
    nblks     = blk_no(ilength(1),nproma)
    IF (PRESENT(out_nblks))  out_nblks  = nblks
    IF (PRESENT(out_npromz)) out_npromz = idx_no(ilength(1),nproma)
    IF (.NOT. ALLOCATED(out_buffer))  ALLOCATE(out_buffer(nproma,nblks))
    CALL nf(nf90_get_var(ncfileID, varID, buffer1D), routine//","//TRIM(name))
    DO i=1,ilength(1)
      jl=idx_no(i,nproma) ; jb=blk_no(i,nproma)
      out_buffer(jl,jb) = buffer1D(i)
    END DO
    IF (PRESENT(opt_filename) .AND. .NOT. PRESENT(opt_ncfileID)) &
      CALL nf(nf90_close(ncfileID), routine)
  END SUBROUTINE read_netcdf_scalar_i1

  !----------------------------------------------------------------
  !> Auxiliary routine. Read scalar INTEGER array. Blocked output.
  !----------------------------------------------------------------
  SUBROUTINE read_netcdf_scalar_i2(name, nproma, out_buffer, out_size, out_nblks, out_npromz, &
    &                              opt_ncfileID, opt_filename)
    CHARACTER(LEN=*),           INTENT(IN)    :: name              !< variable name
    INTEGER,                    INTENT(IN)    :: nproma            !< cache line length
    INTEGER,  ALLOCATABLE,      INTENT(INOUT) :: out_buffer(:,:,:) !< blocked result array.
    INTEGER,          OPTIONAL, INTENT(OUT)   :: out_size          !< number of items
    INTEGER,          OPTIONAL, INTENT(OUT)   :: out_nblks         !< block number of items
    INTEGER,          OPTIONAL, INTENT(OUT)   :: out_npromz        !< size of last block
    INTEGER,          OPTIONAL, INTENT(IN)    :: opt_ncfileID      !< NetCDF file ID (already opened).
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: opt_filename      !< string identifier, user-defined
    !
    CHARACTER(LEN=*), PARAMETER  :: routine = modname//"::read_netcdf_scalar_i2"
    INTEGER, ALLOCATABLE :: buffer2d(:,:)
    INTEGER :: varID, ndims, dimIDs(2), ilength(2), ncfileID, i, jl, jb, nblks

    IF (PRESENT(opt_ncfileID))  ncfileID = opt_ncfileID
    IF (PRESENT(opt_filename) .AND. .NOT. PRESENT(opt_ncfileID)) &
      CALL nf(nf90_open(TRIM(opt_filename), NF90_NOWRITE, ncfileID), routine)
    ! get variable info:
    CALL nf(nf90_inq_varid(ncfileID, TRIM(name), varID), routine//","//TRIM(name))
    CALL nf(nf90_inquire_variable(ncfileID, varID, ndims=ndims), routine//","//TRIM(name))
    ! get dimension info:
    IF (ndims /= 2)  CALL finish(routine, "2D array expected!")
    CALL nf(nf90_inquire_variable(ncfileID, varID, dimids=dimIDs), routine//","//TRIM(name))
    CALL nf(nf90_inquire_dimension(ncfileID, dimids(1), len=ilength(1)), routine//","//TRIM(name))
    CALL nf(nf90_inquire_dimension(ncfileID, dimids(2), len=ilength(2)), routine//","//TRIM(name))
    ! read data:
    IF (.NOT. ALLOCATED(buffer2D)) ALLOCATE(buffer2D(ilength(1), ilength(2)))
    IF (PRESENT(out_size))  out_size  = ilength(1)
    nblks     = blk_no(ilength(1),nproma)
    IF (PRESENT(out_nblks))  out_nblks  = nblks
    IF (PRESENT(out_npromz)) out_npromz = idx_no(ilength(1),nproma)
    IF (.NOT. ALLOCATED(out_buffer))  ALLOCATE(out_buffer(nproma,nblks,ilength(2)))
    CALL nf(nf90_get_var(ncfileID, varID, buffer2D), routine//","//TRIM(name))
    DO i=1,ilength(1)
      jl=idx_no(i,nproma) ; jb=blk_no(i,nproma)
      out_buffer(jl,jb,:) = buffer2D(i,:)
    END DO
    IF (PRESENT(opt_filename) .AND. .NOT. PRESENT(opt_ncfileID)) &
      CALL nf(nf90_close(ncfileID), routine)
  END SUBROUTINE read_netcdf_scalar_i2

END MODULE mo_utilities
