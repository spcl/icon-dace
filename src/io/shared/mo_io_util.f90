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

MODULE mo_io_util

  USE mo_exception,             ONLY: finish
  USE mo_cdi,                   ONLY: FILETYPE_NC, FILETYPE_NC2, FILETYPE_NC4,         &
    &                                 FILETYPE_GRB, FILETYPE_GRB2
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH
  USE mo_util_string,           ONLY: tolower
  USE mo_netcdf_errhandler,     ONLY: nf
  USE mo_netcdf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_filetype
  PUBLIC :: get_file_extension
  PUBLIC :: read_netcdf_int_1d
  PUBLIC :: t_netcdf_att_int

  ! module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_io_util'

  ! Derived data type for integer NetCDF attribute (variable-specific
  ! or global)
  !
  TYPE t_netcdf_att_int
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: var_name  ! variable name (empty string=global attribute)
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: att_name  ! attribute name string
    INTEGER, POINTER :: value ! (result) contains result value after read-in
  END TYPE t_netcdf_att_int

CONTAINS

  !-------------------------------------------------------------------------
  !> @return One of CDI's FILETYPE\_XXX constants. Possible values: 2
  !          (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)
  !
  !  The file type is determined by the setting of the "filetype"
  !  namelist parameter in "initicon_nml". If this parameter has not
  !  been set, we try to determine the file type by its extension
  !  "*.grb*" or ".nc".
  !
  FUNCTION get_filetype(filename)
    INTEGER :: get_filetype
    CHARACTER(LEN=*), INTENT(IN) :: filename
    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//'::get_filetype'
    INTEGER :: idx
    
    idx = INDEX(tolower(filename),'.nc')
    IF (idx==0) THEN
      idx = INDEX(tolower(filename),'.grb')
      IF (idx==0) THEN
        CALL finish(routine, "File type could not be determined!  File: " // trim (filename))
      ELSE
        get_filetype = FILETYPE_GRB2
      END IF
    ELSE
      get_filetype = FILETYPE_NC2
    END IF
  END FUNCTION get_filetype


  !> @return file extension string corresponding to file type
  !
  FUNCTION get_file_extension(filetype) RESULT(extn)
    CHARACTER(LEN=16)   :: extn
    INTEGER, INTENT(IN) :: filetype
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::get_file_extension"

    SELECT CASE (filetype)
    CASE (FILETYPE_NC)
      CALL finish(routine,'netCDF classic not supported')
    CASE (FILETYPE_NC2, FILETYPE_NC4)
      ! this is ok, both formats can write more than 2GB files
      extn = '.nc'
    CASE (FILETYPE_GRB)
      CALL finish(routine,'GRIB1 not supported')
    CASE (FILETYPE_GRB2)
      extn = '.grb'
    CASE default
      CALL finish(routine,'unknown output_type')
    END SELECT
  END FUNCTION get_file_extension


  !----------------------------------------------------------------
  !> Reads integer variables from NetCDF file
  !
  !  Note: opens and closes file and allocates output variables.
  !----------------------------------------------------------------
  SUBROUTINE read_netcdf_int_1d(filename, varname1, var1, opt_varname2, opt_var2, opt_att) 

    CHARACTER(len=*),                 INTENT(IN)    :: filename     ! NetCDF file name
    CHARACTER(len=*),                 INTENT(IN)    :: varname1     ! variable name string
    INTEGER, ALLOCATABLE,             INTENT(INOUT) :: var1(:)      ! output data
    CHARACTER(len=*),       OPTIONAL, INTENT(IN)    :: opt_varname2 ! variable name string
    INTEGER, ALLOCATABLE,   OPTIONAL, INTENT(INOUT) :: opt_var2(:)  ! output data
    TYPE(t_netcdf_att_int), OPTIONAL, INTENT(INOUT) :: opt_att(:)   ! optional attribute values
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine    = "read_netcdf_int_1d"
    INTEGER :: varid, ndims, dimids(1), dimlen, ncfileID, dim_aux, i, iret, optval
    LOGICAL :: l_exist

    ! consistency checks
    IF (PRESENT(opt_var2) .NEQV. PRESENT(opt_varname2)) THEN
      CALL finish(routine, "Internal error!")
    END IF

    ! --- open NetCDF file
    INQUIRE (FILE=filename, EXIST=l_exist)
    IF (.NOT. l_exist) THEN
      CALL finish(routine, 'file "'//TRIM(filename)//'" not found!')
    END IF
    CALL nf(nf_open(TRIM(FILENAME), NF_NOWRITE, ncfileID), routine)

    ! ----------------------
    ! --- variable "var1"
    ! ----------------------

    ! --- find out about variable dimensions:
    CALL nf(nf_inq_varid(ncfileID, TRIM(varname1), varid), routine)
    CALL nf(nf_inq_varndims(ncfileID, varID, ndims), routine)
    IF (ndims /= 1)  CALL finish(routine, "Variable '"//TRIM(varname1)//"' has more than one dimension!")
    CALL nf(nf_inq_vardimid(ncfileID, varID, dimids), routine)
    dim_aux = dimids(1)
    CALL nf(nf_inq_dimlen(ncfileID, dim_aux, dimlen), routine)

    ! --- allocate output variable, read data
    ALLOCATE(var1(dimlen))
    CALL nf(nf_get_var_int(ncfileID, varID, var1), routine)

    ! -------------------------
    ! --- variable "opt_var2"
    ! -------------------------

    IF (PRESENT(opt_var2)) THEN

      ! --- find out about variable dimensions:
      CALL nf(nf_inq_varid(ncfileID, TRIM(opt_varname2), varid), routine)
      CALL nf(nf_inq_varndims(ncfileID, varID, ndims), routine)
      IF (ndims /= 1)  CALL finish(routine, "Variable '"//TRIM(opt_varname2)//"' has more than one dimension!")
      CALL nf(nf_inq_vardimid(ncfileID, varID, dimids), routine)
      dim_aux = dimids(1)
      CALL nf(nf_inq_dimlen(ncfileID, dim_aux, dimlen), routine)
      
      ! --- allocate output variable, read data
      ALLOCATE(opt_var2(dimlen))
      CALL nf(nf_get_var_int(ncfileID, varID, opt_var2), routine)

    END IF

    ! -------------------------------------
    ! --- optional: read integer attributes
    ! -------------------------------------

    IF (PRESENT(opt_att)) THEN

      DO i=1,SIZE(opt_att)
        IF (LEN_TRIM(opt_att(i)%var_name) > 0) THEN
          CALL nf(nf_inq_varid(ncfileID, TRIM(opt_att(i)%var_name), varid), routine)
        ELSE
          varid = NF_GLOBAL
        END IF

        iret = nfx_get_att(ncfileID, varID, TRIM(opt_att(i)%att_name), optval)

        ! Only overwrite previous value if the attribute exists.
        IF (iret == nf_noerr) THEN
          opt_att(i)%value = optval
        END IF
      END DO

    END IF

    ! --- close NetCDF file
    CALL nf(nf_close(ncfileID), routine)

  END SUBROUTINE read_netcdf_int_1d


END MODULE mo_io_util

