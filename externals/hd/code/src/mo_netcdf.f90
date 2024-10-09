! mo_netcdf.f90 - Routines for dealing with IO to NETCDF files
!
! Copyright (C) 2014, MPI-M
! SPDX-License-Identifier: BSD-3-Clause
! See ./LICENSES/ for license information
!
! This file has been modified for the use in HD-model.
!_________________________________________

MODULE mo_netcdf
  !
  !  This routine originates (year 2014) from MPI-ESM, the Earth System Model of the 
  !  Max Planck Institute for Meteorology (Mauritsen et al. 2019). 
  !  Reference: Mauritsen, T., et al. (2019) Developments in the MPI-M Earth System Model 
  !  version 1.2 (MPI-ESM1.2) and its response to increasing CO2. J. Adv. Model. Earth Syst., 11, 
  !  doi: 10.1029/2018MS001400.
  !
  !  - 2020, Stefan Hagemann, HZG: Subroutine IO_GET_VAR_INT added
  !  - Aug. 2024, Stefan Hagemann, Hereon: JSBACH soil layer no.set to 5 (not used in HD)

  USE mo_kind,      ONLY: dp
  USE mo_exception, ONLY: finish, message, message_text
  USE mo_control,   ONLY: ldebugio
  USE mo_filename,  ONLY: NETCDF

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: io_get_varindx     ! mo_memory_base
  PUBLIC :: add_dim            !
  PUBLIC :: add_unknown_dim    !
  PUBLIC :: max_dim_name       ! mo_memory_gl
  PUBLIC :: FILE_INFO          ! mo_io
  PUBLIC :: nf_check           ! checks/print the status of netcdf operations
  PUBLIC :: nf_noerr           ! mo_io
  PUBLIC :: nf_enotnc          ! mo_io
  PUBLIC :: nf_nowrite         ! mo_io
  PUBLIC :: nf_write           ! mo_io
  PUBLIC :: nf_nofill          ! mo_io
  PUBLIC :: nf_max_var_dims    ! mo_read_netcdf
  PUBLIC :: io_ndim_ids        ! mo_io
  PUBLIC :: io_dim_ids         ! mo_io
  PUBLIC :: io_dim             ! io_dim_ids
  PUBLIC :: nf_max_name        ! mo_io
  PUBLIC :: nf_inq_varid       ! mo_io
  PUBLIC :: nf_inq_dimid       ! mo_io
  PUBLIC :: nf_inq_dimlen      !
  PUBLIC :: nf_get_var_double  !
  PUBLIC :: nf_get_vara_double !
  PUBLIC :: nf_get_var_int     !
  PUBLIC :: nf_get_att_int     !
  PUBLIC :: nf_global          ! mo_io
  PUBLIC :: nf_double          ! mo_io
  PUBLIC :: nf_create          ! mo_io
  PUBLIC :: nf__create         ! mo_io, version allows setting I/O buffer size
#ifdef PARALLEL_NC4
  PUBLIC :: nf_create_par      ! mo_io, version allows parallel I/O HDF5 based 
#endif
  PUBLIC :: nf_set_fill        ! mo_io
  PUBLIC :: nf_fill_real       ! mo_io
  PUBLIC :: nf_strerror        ! mo_io
  PUBLIC :: nf_open            ! mo_io
  PUBLIC :: nf__open           ! mo_io, version allows setting I/O buffer size
  PUBLIC :: ncopn              ! mo_io
  PUBLIC :: nf_clobber         ! mo_io
  PUBLIC :: nf_close           ! mo_io
  PUBLIC :: nf_64bit_offset    ! mo_io
!OTBHDF  PUBLIC :: nf_netcdf4         ! mo_io, version allows netcdf I/O HDF5 based 
  PUBLIC :: nf_netcdf4         ! mo_io, version allows netcdf I/O HDF5 based 
  PUBLIC :: nf_def_dim
  PUBLIC :: nf_def_var
  PUBLIC :: nf_def_var_deflate
  PUBLIC :: nf_enddef
  PUBLIC :: nf_inq_vardimid    ! mo_read_ncdf_f77
  PUBLIC :: nf_inq_varndims
  PUBLIC :: nf_inq_dimname
  PUBLIC :: nf_inq_varname
  PUBLIC :: nf_inq_nvars
  PUBLIC :: nf_inq_unlimdim
  PUBLIC :: nf_inq_attlen
  PUBLIC :: nf_get_att_text
  PUBLIC :: nf_get_var1_double
  PUBLIC :: nf_put_var_double
  PUBLIC :: nf_put_var1_double
  PUBLIC :: nf_put_vara_double
  PUBLIC :: nf_put_att_text
  PUBLIC :: nf_unlimited
  PUBLIC :: nf_noclobber
  PUBLIC :: nf_float
  PUBLIC :: nf_real
  PUBLIC :: io_inq_dimid       ! mo_o3clim
  PUBLIC :: io_inq_dimlen      ! mo_o3clim
  PUBLIC :: io_inq_varid       ! mo_o3clim
  PUBLIC :: io_get_var_double  ! mo_o3clim
  PUBLIC :: io_get_var_double1
  PUBLIC :: io_get_vara_double ! mo_o3clim
  PUBLIC :: io_get_var_int     ! mo_couple_to_ocean
  PUBLIC :: io_info_print      ! mo_sst
  PUBLIC :: message_text       ! mo_sst
  PUBLIC :: io_init_dims       ! inictl
  PUBLIC :: io_get_att_int     ! mo_io
  PUBLIC :: io_get_att_double  ! mo_io
  PUBLIC :: io_get_att_text    ! mo_io
  PUBLIC :: io_def_var         ! mo_io
  PUBLIC :: io_put_var_double  ! mo_io
  PUBLIC :: io_put_var_double1
  PUBLIC :: io_put_att_text    ! mo_io
  PUBLIC :: io_put_att_int     ! mo_io
  PUBLIC :: io_put_att_double  ! mo_io
  PUBLIC :: io_enddef          ! mo_io
  PUBLIC :: io_def_dim         ! mo_io
  PUBLIC :: io_info_construct  ! mo_io
  PUBLIC :: t_att_text         ! text attributes data type
  PUBLIC :: put_att            ! store attribute in attributes data type
  PUBLIC :: global_att         ! global text attributes

  PUBLIC :: initialsize        ! initial size of netCDF file
  PUBLIC :: chunksize          ! preferred netCDF I/O buffer size 

  PUBLIC :: cleanup_netcdf     ! deallocate module variables

  INCLUDE 'netcdf.inc'

  PUBLIC :: BELOWSUR, SURFACE, ABOVESUR2, ABOVESUR10, HYBRID, HYBRID_H, TILES, &
            SOILLEV, ROOTZONES, CANOPY

  !-----------------------------------------------------------------------------
  ! 
  ! Due to I/O performance problems, means insufficient return values by
  ! OS information in the stat() system call, we change the buffer size 
  ! to 16 MB fuer netcdf I/O buffer, cannot be PARAMETER, because the netCDF
  ! library is written in C and we get trouble passing parameters.

#if defined (__SX__) || defined (ES) || defined (__PGI)
  INTEGER :: initialsize = 33554432      ! that's 32 MByte   
  INTEGER :: chunksize   = 33554432      ! too
#else
  INTEGER :: initialsize =    32768      ! that's 32 kByte   
  INTEGER :: chunksize   =    32768      ! too
#endif

  !-----------------------------------------------------------------------------
  TYPE FILE_INFO
    
    LOGICAL :: opened                       ! open = .true. or closed = .false.
    LOGICAL :: parallel                     ! parallel .true. else .false.    
    INTEGER :: file_id                      ! netCDF file id 
    INTEGER :: access_mode                  ! access mode for that file
    INTEGER :: ncdims(NF_MAX_VAR_DIMS) 
    INTEGER :: FORMAT                       ! file format NETCDF
    
    CHARACTER(len=NF_MAX_NAME) :: creation_program ! name of this program
    CHARACTER(len=NF_MAX_NAME) :: creation_user    ! who has run this program
    CHARACTER(len=NF_MAX_NAME) :: creation_date    ! created netCDF file at
    CHARACTER(len=NF_MAX_NAME) :: binary_source    ! binary type (CRAY or IEEE)
    CHARACTER(len=NF_MAX_NAME) :: file_type        ! initital or restart file
    CHARACTER(len=NF_MAX_NAME) :: file_name        ! nc file name
    CHARACTER(len=NF_MAX_NAME) :: title
  END TYPE FILE_INFO
  !-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: max_dim_name = 32
  !-----------------------------------------------------------------------------
  TYPE IO_dim
    INTEGER                     :: dim_id   =   0      ! temporary NetCDF id
    INTEGER                     :: var_id   =   0      ! temporary NetCDF id
    INTEGER                     :: dim_len  =  -1
    CHARACTER(len=max_dim_name) :: dim_name =  ''
    CHARACTER(len=64)           :: longname =  ''
    CHARACTER(len=32)           :: units    =  ''
    INTEGER                     :: levtyp   =   0      ! GRIB level type
    LOGICAL                     :: single   =  .FALSE. ! single level
    REAL(dp)  ,POINTER          :: VALUE(:) => NULL()
  END TYPE IO_dim
  !-----------------------------------------------------------------------------
  INTEGER       ,PARAMETER    :: max_dim_ids = 50
  INTEGER                     :: IO_ndim_ids
  TYPE (IO_dim) ,TARGET ,SAVE :: IO_dim_ids (max_dim_ids)

  INTEGER :: BELOWSUR   != 111
  INTEGER :: SOILLEV    !=  71
  INTEGER :: SURFACE    !=   1
  INTEGER :: ABOVESUR2  != 105
  INTEGER :: ABOVESUR10 != 105
  INTEGER :: HYBRID     != 109 ?
  INTEGER :: HYBRID_H   != 110 ?
  INTEGER :: TILES      !=  70
  INTEGER :: ROOTZONES  !=  72
  INTEGER :: CANOPY     !=  73
  !-----------------------------------------------------------------------------
  ! data type to hold file attributes
  !----------------------------------
  TYPE t_att_text
    CHARACTER(len= 32) :: name = ''
    CHARACTER(len=128) :: text = ''
  END TYPE t_att_text

  TYPE (t_att_text) ,SAVE :: global_att (20)
  !-----------------------------------------------------------------------------
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE put_att (att, name, VALUE)
    TYPE (t_att_text), INTENT(inout) :: att(:)
    CHARACTER(len=*),  INTENT(in)    :: name
    CHARACTER(len=*),  INTENT(in)    :: VALUE
    !----------------------------------------------
    ! store an attribute in the attribute data type
    !----------------------------------------------
    INTEGER :: i
    DO i = 1, SIZE(att)
      IF (att(i)% name == '' .OR. att(i)% name == name) THEN
        att(i)% name = name
        att(i)% text = VALUE
        EXIT
      ENDIF
    END DO
  END SUBROUTINE put_att
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_info_construct(fileinfo)
    TYPE(FILE_INFO), INTENT(inout)  :: fileinfo

    fileinfo%opened           = .FALSE.
    fileinfo%parallel         = .FALSE.
    fileinfo%file_id          = 0
    fileinfo%access_mode      = 0
    fileinfo%ncdims(:)        = 0
    fileinfo%creation_program = '?'
    fileinfo%creation_user    = '?'
    fileinfo%creation_date    = '?'
    fileinfo%binary_source    = '?'
    fileinfo%file_type        = '?'
    fileinfo%file_name        = '?'
    fileinfo%title            = '?'

  END SUBROUTINE IO_info_construct
  !-----------------------------------------------------------------------------
  FUNCTION IO_get_varindx (dimname) RESULT(indx)
    CHARACTER(len=*), INTENT(in) :: dimname
    INTEGER :: indx

    DO indx = 1, IO_ndim_ids
      IF (IO_dim_ids(indx)%dim_name == dimname) RETURN
    END DO
    
    indx = -1
    
  END FUNCTION IO_get_varindx
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_init_dims
    !---------------------
    ! predifine dimensions
    !---------------------
    USE mo_control,    ONLY: ngl, nhgl, nlon, nlev, nlevp1, nsp, nvclev, nmp1

    INTEGER :: jpgrnd = 5   ! hard coded number of soil layers from JSBACH - not relevant in HD

    IO_ndim_ids = 0
print*,'OTBPR in IO_init_dims'
    !-------------------------------
    ! definitions used in rerun file
    !-------------------------------
    CALL add_dim ("lon",     nlon, "longitude", "degrees_east" )
    CALL add_dim ("lat",     ngl,  "latitude" , "degrees_north")
    CALL add_dim ("nhgl",    nhgl)
    CALL add_dim ("nlevp1",  nlevp1)
    CALL add_dim ("spc",     nsp,  "spectral coefficients")
    CALL add_dim ("nvclev" , nvclev)
    CALL add_dim ("complex", 2,    "real, imaginary part")
    CALL add_dim ("nmp1",    nmp1  )
    CALL add_dim ("belowsurface",  jpgrnd, "levels below the surface", &
                                        levtyp = 111,                  &
                                        indx   = BELOWSUR,             &
                                        units  = 'cm',                 &
                                        VALUE  = (/  3.0_dp, &
                                                    19.0_dp, &
                                                    78.0_dp, &
                                                   268.0_dp, &
                                                   698.0_dp/))
    !-------------------------------------------------------
    ! different/additional definitions used in output stream
    !-------------------------------------------------------
    CALL add_dim ("lev",    nlev,   "hybrid level at layer midpoints", &
                                                    units  = "level",  &
                                                    levtyp = 109,      &
                                                    indx   = HYBRID)
    CALL add_dim ("ilev",   nlev+1, "hybrid level at layer interfaces",&
                                                    units  = "level",  &
                                                    levtyp = 109,      &
                                                    indx   = HYBRID_H)
    !    call add_dim ("time",nf_unlimited)
    !--------------------
    ! single level fields
    !--------------------
    CALL add_dim ("surface",1,                      levtyp =    1,     &
                                                    single = .TRUE.,   &
                                                    indx   = SURFACE)
    CALL add_dim ("height2m",     1,                units  = "m",      &
                                                    levtyp = 105,      &
                                                    single = .TRUE.,   &
                                                    VALUE  = (/2.0_dp/),   &
                                                    indx   = ABOVESUR2)
    CALL add_dim ("height10m",    1,                units  = "m",      &
                                                    levtyp = 105,      &
                                                    single = .TRUE.,   &
                                                    VALUE  = (/10.0_dp/),  &
                                                    indx   = ABOVESUR10)
  END SUBROUTINE IO_init_dims
  !-----------------------------------------------------------------------------
  SUBROUTINE add_dim (name, len, longname, units, levtyp, single, VALUE, indx)
    !--------------------------------------------------
    ! define a new dimension for Netcdf and GRIB output
    !--------------------------------------------------
    CHARACTER(len=*) ,INTENT(in)            :: name      ! mnemonic
    INTEGER          ,INTENT(in)            :: len       ! size of dimension
    CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: longname  ! long name
    CHARACTER(len=*) ,INTENT(in)  ,OPTIONAL :: units     ! units
    INTEGER          ,INTENT(in)  ,OPTIONAL :: levtyp    ! GRIB level type
    LOGICAL          ,INTENT(in)  ,OPTIONAL :: single    ! single layer flag
    REAL(dp)         ,INTENT(in)  ,OPTIONAL :: VALUE (:) ! coordinates
    INTEGER          ,INTENT(out) ,OPTIONAL :: indx      ! index in IO_dim_ids 

    INTEGER                :: idx, i
    TYPE (IO_dim) ,POINTER :: dim
    !--------------------
    ! check for zero size
    !--------------------
    IF (len == 0) CALL finish ('add_dim',TRIM(name)//' len = 0')
    !-----------------------------------------------------------
    ! if dimension is already defined it must have the same size
    !-----------------------------------------------------------
    idx = 0
print*,'NETCDF ', len
!   IO_ndim_ids = SIZE(IO_dim_ids)
    DO i = 1, IO_ndim_ids
      IF (IO_dim_ids(i)% dim_name == name) THEN
        idx = i
        IF (IO_dim_ids(i)% dim_len /= len) &
          CALL finish ('add_dim',TRIM(name)//' redifined with different len')
        EXIT
      ENDIF
    END DO
    !--------------------------------------
    ! increase number of dimensions defined
    !--------------------------------------
    IF (idx == 0) THEN
      IO_ndim_ids = IO_ndim_ids + 1
      IF (IO_ndim_ids > SIZE(IO_dim_ids)) &
        CALL finish ('add_dim','increase size of IO_dim_ids')
      idx = IO_ndim_ids
    ENDIF
    dim => IO_dim_ids(idx)
    !---------------------------------------------------------
    ! define attributes, chose defaults or optional parameters
    !---------------------------------------------------------
    dim% dim_id   =   0   
    dim% var_id   =   0
    dim% dim_len  =  len
    dim% dim_name =  name
    dim% longname =  ''     ; IF (PRESENT (longname)) dim% longname = longname
    dim% units    =  ''     ; IF (PRESENT (units   )) dim% units    = units
    dim% levtyp   =   0     ; IF (PRESENT (levtyp  )) dim% levtyp   = levtyp
    dim% single   =  .FALSE.; IF (PRESENT (single  )) dim% single   = single
    IF (ASSOCIATED (dim% VALUE)) DEALLOCATE (dim% VALUE)
    IF (PRESENT (VALUE)) THEN
      IF (SIZE (VALUE) /= dim% dim_len) &
        CALL finish ('add_dim',TRIM(name)//' len /= size (value)')
      ALLOCATE (dim% VALUE (SIZE (VALUE)))
      dim% VALUE = VALUE
    ENDIF
    IF (PRESENT(indx)) indx = idx
  END SUBROUTINE add_dim
  !-----------------------------------------------------------------------------
  SUBROUTINE add_unknown_dim (len, indx)
    !-----------------------------------------
    ! add a dummy entry in the dimension table
    !-----------------------------------------
    INTEGER, INTENT(in)  :: len  ! size of dimension
    INTEGER, INTENT(out) :: indx ! index of entry in table

    CHARACTER(len=4) :: nam
    CHARACTER(len=4) :: form
    INTEGER          :: i

    SELECT CASE (len)
    CASE (1:9)
      form ='(i1)'
    CASE (10:99)
      form ='(i2)'
    CASE (100:999)
      form ='(i3)'
    CASE default
      CALL finish ('add_unknown_dim','len < 1 or len > 999')
    END SELECT
    nam = 'n' 
    WRITE (nam(2:4),form) len
    DO i=1, IO_ndim_ids
      IF (IO_dim_ids(i)% dim_name == nam) THEN
        indx = i
        RETURN
      ENDIF
    END DO
    CALL add_dim (nam, len, levtyp=109, indx=indx)    
  END SUBROUTINE add_unknown_dim
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_DEF_DIM (ncid, name, len, dimid)
    INTEGER :: ncid, len, dimid, status, indx
    CHARACTER(len=*) :: name

    status = NF_DEF_DIM (ncid, name, len, dimid)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, name, len
      CALL message ('IO_DEF_DIM', message_text)
      CALL message ('IO_DEF_DIM', NF_STRERROR(status))
      CALL finish  ('IO_DEF_DIM', 'Run terminated.')
    END IF

    DO indx = 1, IO_ndim_ids
      IF (name == IO_dim_ids(indx)%dim_name) THEN
        IO_dim_ids(indx)%dim_id = dimid
        RETURN
      END IF
    END DO

    WRITE (message_text,*) 'element ',name,' not available ...'
    CALL message ('IO_DEF_DIM', message_text)
    CALL finish('IO_DEF_DIM','IO_dim_ids error')

  END SUBROUTINE IO_DEF_DIM
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_INQ_DIMID (ncid, name, dimid)
    INTEGER :: ncid, dimid, status
    CHARACTER(len=*) :: name

    status = NF_INQ_DIMID (ncid, name, dimid)
    IF (status /= NF_NOERR) THEN
      CALL message ('IO_INQ_DIMID', name)
      CALL message ('IO_INQ_DIMID', NF_STRERROR(status))
      CALL finish  ('IO_INQ_DIMID', 'Run terminated.')
    END IF

  END SUBROUTINE IO_INQ_DIMID
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_INQ_DIMLEN (ncid, dimid, len)
    INTEGER :: ncid, dimid, len, status

    status = NF_INQ_DIMLEN (ncid, dimid, len)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, dimid
      CALL message ('IO_INQ_DIMLEN', NF_STRERROR(status))
      CALL finish  ('IO_INQ_DIMLEN', 'Run terminated.')
    END IF

  END SUBROUTINE IO_INQ_DIMLEN
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_INQ_VARID (ncid, name, varid, lstat)
    INTEGER :: ncid, varid, status
    CHARACTER(len=*) :: name
    LOGICAL, INTENT(out), OPTIONAL :: lstat

    status = NF_INQ_VARID (ncid, name, varid)

    IF (PRESENT(lstat)) THEN
      lstat = status == NF_NOERR
      RETURN
    ENDIF

    IF (ldebugio) THEN
      WRITE(message_text,*) '     Id=',ncid,' varid=',varid,'   ',name
      CALL message ('IO_INQ_VARID', message_text)
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, name
      CALL message ('IO_INQ_VARID', message_text)
      CALL message ('IO_INQ_VARID', NF_STRERROR(status))
      CALL finish  ('IO_INQ_VARID', 'Run terminated.')
    END IF

  END SUBROUTINE IO_INQ_VARID
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_DEF_VAR (ncid, name, xtype, nvdims, vdims, varid, lparallel)
    INTEGER :: ncid, varid, xtype, nvdims
    CHARACTER(len=*) :: name
    INTEGER :: vdims(*)
    INTEGER :: status
    LOGICAL, OPTIONAL, INTENT(in) :: lparallel

    status = NF_DEF_VAR (ncid, name, xtype, nvdims, vdims, varid)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, name, xtype, nvdims,vdims(1:nvdims)
      CALL message ('IO_DEF_VAR', message_text)
      CALL message ('IO_DEF_VAR', NF_STRERROR(status))
      CALL finish  ('IO_DEF_VAR', 'Run terminated.')
    ENDIF
    IF (PRESENT(lparallel)) THEN
      IF (lparallel) THEN
#ifdef PARALLEL_NC4
        status = NF_VAR_PAR_ACCESS(ncid, varid, NF_INDEPENDENT)
        IF (status /= NF_NOERR) THEN
          WRITE(message_text,*) ncid, name, xtype, nvdims,vdims(1:nvdims)
          CALL message ('IO_DEF_VAR collective', message_text)
          CALL message ('IO_DEF_VAR collective', NF_STRERROR(status))
          CALL finish  ('IO_DEF_VAR collective', 'Run terminated.')
        ENDIF
#else
        CALL finish('IO_DEV_VAR', 'Parallel write requires -DPARALLEL_NC4 enabled during compilation')
#endif
      ENDIF
    ENDIF

  END SUBROUTINE IO_DEF_VAR
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_ATT_TEXT (ncid, varid, name, text)
    INTEGER :: ncid, varid, TYPE, len
    CHARACTER(len=*) :: name, text
    INTEGER :: status

    status = NF_INQ_ATT(ncid, varid, name, TYPE, len)
    IF (status /= NF_NOERR) THEN  
      name = 'not available'
      RETURN
    ENDIF

    status = NF_GET_ATT_TEXT (ncid, varid, name, text)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid, name
      CALL message ('IO_GET_ATT_TEXT', message_text)
      CALL message ('IO_GET_ATT_TEXT', NF_STRERROR(status))
      CALL finish  ('IO_GET_ATT_TEXT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_ATT_TEXT
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_PUT_ATT_TEXT (ncid, varid, name, text)
    INTEGER :: ncid, varid
    CHARACTER(len=*) :: name, text
    INTEGER :: status
    INTEGER :: lentext

    lentext = LEN(TRIM(text))

    status = NF_PUT_ATT_TEXT (ncid, varid, name, lentext, text)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid, name, lentext, text
      CALL message ('IO_PUT_ATT_TEXT', message_text)
      CALL message ('IO_PUT_ATT_TEXT', NF_STRERROR(status))
      CALL finish  ('IO_PUT_ATT_TEXT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_ATT_TEXT
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_ATT_INT (ncid, varid, name, ival)
    INTEGER :: ncid, varid
    CHARACTER(len=*) :: name
    INTEGER :: ival, ival_buf(1)
    INTEGER :: status

    status = NF_GET_ATT_INT (ncid, varid, name, ival_buf)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid, name
      CALL message ('IO_GET_ATT_INT', message_text)
      CALL message ('IO_GET_ATT_INT', NF_STRERROR(status))
      CALL finish  ('IO_GET_ATT_INT', 'Run terminated.')
    END IF

    ival = ival_buf(1)

  END SUBROUTINE IO_GET_ATT_INT
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_PUT_ATT_INT (ncid, varid, name, ival)
    INTEGER :: ncid, varid
    CHARACTER(len=*) :: name
    INTEGER :: ival, ival_buf(1)
    INTEGER :: xtype, len
    INTEGER :: status

    len = 1
    xtype = NF_INT
    ival_buf(1) = ival


    status = NF_PUT_ATT_INT (ncid, varid, name, xtype, len, ival_buf)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid, name
      CALL message ('IO_PUT_ATT_INT', message_text)
      CALL message ('IO_PUT_ATT_INT', NF_STRERROR(status))
      CALL finish  ('IO_PUT_ATT_INT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_ATT_INT
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_ATT_DOUBLE (ncid, varid, name, rval)
    INTEGER :: ncid, varid
    CHARACTER(len=*) :: name
    REAL(dp) :: rval, rval_buf(1)
    INTEGER :: status

    status = NF_GET_ATT_DOUBLE (ncid, varid, name, rval_buf)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid,  name
      CALL message ('IO_GET_ATT_DOUBLE', message_text)
      CALL message ('IO_GET_ATT_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_GET_ATT_DOUBLE', 'Run terminated.')
    END IF

    rval = rval_buf(1)

  END SUBROUTINE IO_GET_ATT_DOUBLE
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_PUT_ATT_DOUBLE (ncid, varid, name, rval)
    INTEGER :: ncid, varid
    CHARACTER(len=*) :: name
    REAL(dp) :: rval, rval_buf(1)
    INTEGER :: xtype, len
    INTEGER :: status

    len = 1
    xtype = NF_DOUBLE
    rval_buf(1) = rval

    status = NF_PUT_ATT_DOUBLE (ncid, varid, name, xtype, len, rval_buf)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid, name
      CALL message ('IO_PUT_ATT_DOUBLE', message_text)
      CALL message ('IO_PUT_ATT_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_PUT_ATT_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_ATT_DOUBLE
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_VAR_DOUBLE (ncid, varid, dvals)
    INTEGER :: ncid, varid
    REAL(dp) :: dvals(*)
    INTEGER :: status

    status = NF_GET_VAR_DOUBLE (ncid, varid, dvals)

    IF (ldebugio) THEN
      WRITE(message_text,*) ' Id=',ncid,' varid=',varid
      CALL message ('IO_GET_VAR_DOUBLE', message_text)
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*)  ncid, varid
      CALL message ('IO_GET_VAR_DOUBLE', message_text)
      CALL message ('IO_GET_VAR_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_GET_VAR_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_VAR_DOUBLE
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_VAR_DOUBLE1 (ncid, varid, dval)
    INTEGER :: ncid, varid
    REAL(dp) :: dval
    INTEGER :: status

    REAL(dp) :: temp_dval(1)

    status = NF_GET_VAR_DOUBLE (ncid, varid, temp_dval)
    dval = temp_dval(1)

    IF (ldebugio) THEN
      WRITE(message_text,*) 'IO_GET_VAR_DOUBLE :',' Id=',ncid,' varid=',varid
      CALL message ('IO_GET_VAR_DOUBLE', message_text)
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid
      CALL message ('IO_GET_VAR_DOUBLE', message_text)
      CALL message ('IO_GET_VAR_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_GET_VAR_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_VAR_DOUBLE1
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_VARA_DOUBLE (ncid, varid, istart, icount, dvals)
    INTEGER :: ncid, varid
    INTEGER :: istart(*), icount(*)
    REAL(dp) :: dvals(*)
    INTEGER :: status

    status = NF_GET_VARA_DOUBLE (ncid, varid, istart, icount, dvals)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid, istart(1:4), icount(1:4)
      CALL message ('IO_GET_VARA_DOUBLE', message_text)
      CALL message ('IO_GET_VARA_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_GET_VARA_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_VARA_DOUBLE
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_GET_VAR_INT (ncid, varid, dvals)
    INTEGER :: ncid, varid
    INTEGER :: dvals(*)
    INTEGER :: status

    status = NF_GET_VAR_INT (ncid, varid, dvals)

    IF (ldebugio) THEN
      WRITE(message_text,*) ' Id=',ncid,' varid=',varid
      CALL message ('IO_GET_VAR_INT', message_text)
    END IF

    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*)  ncid, varid
      CALL message ('IO_GET_VAR_INT', message_text)
      CALL message ('IO_GET_VAR_INT', NF_STRERROR(status))
      CALL finish  ('IO_GET_VAR_INT', 'Run terminated.')
    END IF

  END SUBROUTINE IO_GET_VAR_INT
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_PUT_VARA_DOUBLE (ncid, varid, istart, icount, dvals)
    INTEGER :: ncid, varid
    INTEGER :: istart(*), icount(*)
    REAL(dp) :: dvals(*)
    INTEGER :: status
    CHARACTER(len=80) :: varname

    status = NF_PUT_VARA_DOUBLE (ncid, varid, istart, icount, dvals)
    IF (status /= NF_NOERR) THEN
      status = NF_INQ_VARNAME (ncid, varid, varname)
      WRITE(message_text,*) ncid, varid, varname
      CALL message ('IO_PUT_VARA_DOUBLE', message_text)
      CALL message ('IO_PUT_VARA_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_PUT_VARA_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_VARA_DOUBLE
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_PUT_VAR_DOUBLE (ncid, varid, dvals)
    INTEGER :: ncid, varid
    REAL(dp) :: dvals(*)
    INTEGER :: status
    CHARACTER(len=80) :: varname

    IF (ldebugio) THEN
      status = NF_INQ_VARNAME (ncid, varid, varname)
      WRITE(message_text,*) ncid, varid, varname
        CALL message ('IO_PUT_VAR_DOUBLE', message_text)
      IF (status /= NF_NOERR) THEN
        WRITE(message_text,*) ncid, varid
        CALL message ('IO_PUT_VAR_DOUBLE', message_text)
        CALL message ('IO_PUT_VAR_DOUBLE', NF_STRERROR(status))
        CALL finish  ('IO_PUT_VAR_DOUBLE', 'Run terminated.')
      END IF
    END IF

    status = NF_PUT_VAR_DOUBLE (ncid, varid, dvals)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) ncid, varid
      CALL message ('IO_PUT_VAR_DOUBLE', message_text)
      CALL message ('IO_PUT_VAR_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_PUT_VAR_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_VAR_DOUBLE
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_PUT_VAR_DOUBLE1 (ncid, varid, dval)
    INTEGER :: ncid, varid
    REAL(dp) :: dval
    INTEGER :: status
    CHARACTER (80) :: varname

    REAL(dp) :: temp_dval(1)

    IF (ldebugio) THEN
      status = NF_INQ_VARNAME (ncid, varid, varname)
      WRITE(message_text,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid, varname
      IF (status /= NF_NOERR) THEN
        WRITE(message_text,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
        CALL message ('IO_PUT_VAR_DOUBLE', NF_STRERROR(status))
        CALL finish  ('IO_PUT_VAR_DOUBLE', 'Run terminated.')
      END IF
    END IF

    temp_dval(1) = dval
    status = NF_PUT_VAR_DOUBLE (ncid, varid, temp_dval)
    IF (status /= NF_NOERR) THEN
      WRITE(message_text,*) 'IO_PUT_VAR_DOUBLE :', ncid, varid
      CALL message ('IO_PUT_VAR_DOUBLE', NF_STRERROR(status))
      CALL finish  ('IO_PUT_VAR_DOUBLE', 'Run terminated.')
    END IF

  END SUBROUTINE IO_PUT_VAR_DOUBLE1
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_ENDDEF (ncid)
    INTEGER :: ncid, status

    status = NF_ENDDEF (ncid)
    IF (status /= NF_NOERR) THEN
      CALL message ('IO_ENDDEF', NF_STRERROR(status))
      CALL finish  ('IO_ENDDEF', 'Run terminated.')
    END IF

  END SUBROUTINE IO_ENDDEF
  !-----------------------------------------------------------------------------
  SUBROUTINE IO_info_print(fileinfo)
    TYPE(FILE_INFO) :: fileinfo
    INTEGER :: i

    CALL message('',' ')
    WRITE(message_text,*) &
         'File name ',TRIM(fileinfo%file_name), &
         ' type ',     TRIM(fileinfo%file_type), &
         ' title ',    TRIM(fileinfo%title)
    CALL message('',message_text)
    SELECT CASE (fileinfo%format)
    CASE (NETCDF)
      CALL message('',' File type   : NETCDF')
    CASE default
      CALL message('',' File type   : UNKNOWN')
    END SELECT
    IF (fileinfo%opened) THEN
       CALL message('',' File open.')
    ELSE
       CALL message('',' File closed.')
    END IF
    WRITE(message_text,*) &
         'File id :',   fileinfo%file_id,&
         ' access mode :',fileinfo%access_mode
    CALL message('',message_text)
    DO i = 1, NF_MAX_VAR_DIMS
       IF (fileinfo%ncdims(i) /= 0) THEN
          WRITE(message_text,*) &
               'Variable id :',i,' field dimension :',fileinfo%ncdims(i)
          CALL message('',message_text)
       END IF
    END DO
    WRITE(message_text,*) &
         'Created by ',TRIM(fileinfo%creation_user), &
         ' with ',     TRIM(fileinfo%creation_program), &
         ' at ',       TRIM(fileinfo%creation_date), &
         ' source ',   TRIM(fileinfo%binary_source)
    CALL message('',message_text)
    CALL message('',' ')

  END SUBROUTINE IO_info_print
  !-----------------------------------------------------------------------------
  SUBROUTINE cleanup_netcdf
    !
    ! deallocate module variables
    !
    INTEGER :: i

    DO i=1, SIZE(IO_dim_ids)
      IF (ASSOCIATED (IO_dim_ids(i)% VALUE)) DEALLOCATE (IO_dim_ids(i)% VALUE)
    END DO

  END SUBROUTINE cleanup_netcdf
  !-----------------------------------------------------------------------------
  ! turns nf_* function into subroutine and checks status
  SUBROUTINE nf_check(status, fname) 
    INTEGER :: status
    CHARACTER(len=*), OPTIONAL :: fname

    IF (status /= nf_noerr) THEN
      IF (PRESENT(fname)) THEN
        CALL finish('netcdf error in '//TRIM(fname),nf_strerror(status))
      ELSE
        CALL finish('netcdf error',nf_strerror(status))
      END IF
    ENDIF

  END SUBROUTINE nf_check
  !-----------------------------------------------------------------------------
END MODULE mo_netCDF
