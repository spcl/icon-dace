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
!
! Background routines to read and write NetCDF output
! All calls to the netcdf library should be directed through here.
! The module opens (with open_nc), closes (with close_nc),
! and writes (with writevar_nc) anything between 0D (e.g. timeseries)
! and 4D (e.g. z,x,y,t) fields to file.
!
MODULE mo_write_netcdf

  use mo_kind
  use mo_exception, only: finish
  use mo_io_config,   only: lsync=>lkeep_in_sync
  USE mo_netcdf

  implicit none
  private

  real(wp), parameter    :: fillvalue_double = -32678. !< Fill value
  integer(i4), parameter :: fillvalue_int    = -32678  !< Fill value
  integer, parameter :: icomp_dbl   = 1 !< Write back doubles instead of floats
  integer, parameter :: icomp_int   = 2 !< Write back integers instead of floats

  public :: open_nc, close_nc
  public :: writevar_nc, addvar_nc

!> Interface to write into a netcdf file
  interface writevar_nc
    module procedure writevar0D_nc
    module procedure writevar1D_nc
  end interface writevar_nc

!> Interface to put attributes in a netcdf file
  interface putatt_nc
    module procedure putatt_single_nc
    module procedure putatt_double_nc
    module procedure putatt_int_nc
    module procedure putatt_short_nc
    module procedure putatt_str_nc
  end interface putatt_nc

!> Interface to put attributes in a netcdf file
  interface getatt_nc
    module procedure getatt_single_nc
    module procedure getatt_double_nc
    module procedure getatt_int_nc
    module procedure getatt_short_nc
    module procedure getatt_str_nc
    module procedure getatt_str_fname_nc
  end interface getatt_nc

contains

!-------------------------------------------------------------------------
!> Subroutine Open_NC: Opens a NetCDF File for writing
!! If the file already exists, the record number is being set to the current
!! time in the simulation
  subroutine open_nc(fname, ncid,nrec, rtimee, ldelete)
    integer, intent (out)          :: ncid   !< NetCDF file number
    integer, intent (out),optional :: nrec   !< The record number that corresponds to the current simulation time
    real(wp), intent(in), optional :: rtimee !< Simulation time
    logical, intent(in), optional  :: ldelete !< Whether to delete an old file when present
    character (len=*), intent (in) :: fname  !< File name
    integer :: ndims, nvars, n, nn, xtype, dimsize(7)
    integer, allocatable, dimension(:) :: start, dimids

    character (len=12) :: date, time
    integer            :: iret,  ncall, RecordDimID,timeid, totdimsize
    real(wp), allocatable  :: xtimes(:)
    logical                :: exists, ldef

    if (present(rtimee)) then ! Check whether file exists
      inquire(file=trim(fname),exist=exists)
      if (exists .and. present(ldelete)) then
        if (ldelete) then
          exists = .false.
          open(1, file=trim(fname), status='old')
          close(1,status='delete')
        end if
      end if
    else
      exists = .false.
    end if

    ncall = 0
    if (.not.exists) then ! Write the header of the file, and ensure we're in data mode

      call date_and_time(date,time)
      iret = nf90_create(trim(fname),nf90_noclobber,ncid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncid, 'title', trim(fname))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncid,'history','Created on '//date(1:8)//' at '//trim(time))
      ldef = ensuredata_nc(ncid)

    else ! If the file exists, only the record number needs to be set

      iret  = nf90_open (trim(fname), NF90_WRITE, ncid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret  = nf90_inquire_dimension(ncid, RecordDimID, len = nrec)

      ncall = nrec


    end if

    if (present(nrec)) nrec = ncall
    iret = nf90_sync(ncid)

  end subroutine open_nc

!-------------------------------------------------------------------------
!> Switch NetCDF file to define mode. Returns true if the dataset already
!  was in define mode, and false if it was in data mode
  logical function ensuredefine_nc(ncid)
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    ensuredefine_nc = .false.
    iret = nf90_redef(ncid)
    select case(iret)
    case(nf90_noerr)
      ensuredefine_nc = .false.
    case(nf90_eindefine)
      ensuredefine_nc = .true.
    case default
      call nchandle_error(ncid, iret)
    end select
 end function ensuredefine_nc

!-------------------------------------------------------------------------
!> Switch NetCDF file to data mode. Returns true if the dataset was in define
!  mode, and false if it was already in data mode
  logical function ensuredata_nc(ncid)
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    ensuredata_nc = .false.
    iret = nf90_enddef(ncid)
    select case(iret)
    case(nf90_noerr)
      ensuredata_nc = .true.
    case(nf90_enotindefine)
      ensuredata_nc = .false.
    case default
      call nchandle_error(ncid, iret)
    end select
 end function ensuredata_nc

!-------------------------------------------------------------------------
!> Switch NetCDF file to data mode
  subroutine enddefine_nc(ncid)
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret
    iret = nf90_enddef(ncid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
  end subroutine enddefine_nc

!-------------------------------------------------------------------------
!> Close a NetCDF file
  subroutine close_nc(ncid)
    integer, intent(in) :: ncid !< NetCDF file number
    integer :: iret, iformat
    iret = nf90_inquire(ncid, formatNum = iformat)
    if (iret /= nf90_noerr) return
    iret = nf90_close(ncid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
  end subroutine close_nc

!-------------------------------------------------------------------------
!> Add a variable to a NetCDF file. The procedure checks whether the variable(name) is already present,
!! and if not, adds it to the file. The same holds for the dimension(names) of the variable.
  subroutine addvar_nc(ncID, name, lname, unit, dimname, dimlongname, dimunit, dimsize, dimvalues, &
                       icompress)
    integer, intent (in)                              :: ncID        !< NetCDF file number
    character (*), intent (in)                        :: name        !< Netcdf name of the variable
    character (*), intent (in)                        :: lname       !< Longname of the variable
    character (*), intent (in)                        :: unit        !< Unit of the variable
    character (*), dimension(:), intent(in), optional :: dimname     !< NetCDF names of the dimensions of the variable (in array form)
    character (*), dimension(:), intent(in), optional :: dimlongname !< Longnames of the dimensions of the variable (in array form)
    character (*), dimension(:), intent(in), optional :: dimunit     !< Units of the dimensions of the variable (in array form)
    integer, dimension(:), intent(in), optional       :: dimsize     !< List of dimension sizes; 0 for unlimited
    real(wp), dimension(:,:), intent(in), optional    :: dimvalues   !< List of values of the dimension
    integer, optional                                 :: icompress   !< Choice of compression: 1 for double,
                                                                     !   2 for integer

    integer                                     :: iret, n, nrdim,VarID, icomp, datatype
    integer, allocatable, dimension(:)          :: ncdim
    logical                                     :: ldef

    ldef = ensuredefine_nc(ncid) ! Need to be in define mode

    if (present(dimname)) then
      nrdim = size(dimname)
    else
      nrdim = 0
    end if

    if (present(icompress)) then
      icomp = icompress
    else
      icomp = icomp_dbl
    end if

    allocate (ncdim(nrdim))

    do n=1,nrdim ! For every dimension, check whether it exists in the file. If not, inqdim_nc will write it.
      if (present(dimvalues)) then
        ncdim(n) = inqdimid_nc(ncid ,dimname(n) ,dimlongname(n) ,dimunit(n) ,dimvalues(1:dimsize(n),n))
      else
        ncdim(n) = inqdimid_nc(ncid ,dimname(n) ,dimlongname(n) ,dimunit(n), (/0._wp/))
      end if
    end do

    iret=nf90_inq_varid(ncid,trim(validate(name)),varid)

    if (iret /= nf90_noerr) then ! If the inquiry fails, the variable needs to be created

      select case (icomp)
      case(icomp_int)
        datatype = nf90_int
      case(icomp_dbl)
        datatype = nf90_double
      case default
        CALL finish('mo_write_netcdf:addvar_nc','wrong datatype!')
      end select

      iret=nf90_def_var(ncID, trim(validate(name)), datatype, ncdim, VarID)

      if (iret/=0) then
        write (*,*) 'Variable ',trim(name), ncdim
        call nchandle_error(ncid, iret)
      end if

      iret = putatt_nc(ncID, 'longname', lname, trim(name))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

      iret = putatt_nc(ncID, 'units', unit, trim(name))
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

      select case (icomp)
      case (icomp_int)
        iret = putatt_nc(ncID, '_FillValue', fillvalue_int, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      case (icomp_dbl)
        iret = putatt_nc(ncID, '_FillValue', fillvalue_double, trim(name))
        if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      end select

    end if

    if (ldef .eqv. .false.) then
      call enddefine_nc(ncid) ! Leave define mode again
    end if

    deallocate (ncdim)
  end subroutine addvar_nc

!-------------------------------------------------------------------------
!> Check whether in the give file a dimension with the given dimension name already exists.
!! If so, check whether the extent and values are correct.
!! If the dimension does not exist, create it and fill the related variable.
!! Inqdimid_nc sets the netcdf number of the resulting dimension as a return value.
  function inqdimid_nc(ncid,dimname,dimlongname,dimunit,dimvalues)
    integer                         :: inqdimid_nc !< The netcdf dimension number (return value)
    integer, intent(in)             :: ncid        !< Netcdf file number
    character(*), intent(in)        :: dimname     !< Name of the dimension
    character(*), intent(in)        :: dimlongname !< Netcdf long name of the dimension
    character(*), intent(in)        :: dimunit     !< Unit of the dimension
    real(wp), dimension(:), intent(in)  :: dimvalues   !< Values of the dimension. 0 for the unlimited dimension
    logical :: ltime, ldef
    integer :: iret,varid, dimsize
    real(wp), dimension(:), allocatable :: dimvar
    character(LEN=80) :: fname

    if (all(dimvalues == 0)) then ! If all values of this dimension are zero,
                                  ! it is the unlimited (time) dimension
      ltime = .true.
    else
      ltime = .false.
    end if

    iret = 0
    iret = nf90_inq_dimid(ncid, trim(validate(dimname)), inqdimid_nc)

    if (iret == 0) then !If the dimension already exists....

      if (.not. ltime) then
        iret = nf90_inquire_dimension(ncID, inqdimid_nc, len = dimsize)

        allocate(dimvar(dimsize))

        if (dimsize == size(dimvalues)) then ! Check whether the number of points along this dimension is correct
          ldef = ensuredata_nc(ncid)
          iret =  nf90_inq_varid(ncid, trim(validate(dimname)), varid)
          iret  = nf90_get_var(ncid,varid,dimvar)
          ldef = ensuredefine_nc(ncid)
          ! Check whether dimension in file matches with the desired values
          if (any((abs((dimvar - dimvalues)/(dimvar+epsilon(1.)))) > 1e-4)) then
            iret = getatt_nc(ncid,'title', fname)
            if (iret /= nf90_noerr) fname = ''
            print *, 'NetCDF error in file ' // trim(fname)
            print *, 'For dimension ', trim(dimname)
            print *, "Inqdimid_nc: Dimensions don't match with what's already in the file"
            call finish('mo_write_netcdf:inqdimid_nc','stoppping!')
          end if
        else
          iret = getatt_nc(ncid, 'title', fname)
          if (iret /= nf90_noerr)  fname = ''
          print *, 'NetCDF error in file ' // trim(fname)
          print *, 'For dimension ', trim(dimname)
          print *, "Number of points doesn't fit"
          call finish('mo_write_netcdf:inqdimid_nc','stoppping!')
        end if

        deallocate(dimvar)
      end if

    else !If the dimension does not exist yet, we need to create it.

      if (ltime) then
        iret = nf90_def_dim(ncID, trim(validate(dimname)), NF90_UNLIMITED, inqdimid_nc)
      else
        iret = nf90_def_dim(ncID, trim(validate(dimname)), size(dimvalues), inqdimid_nc)
      endif

      if (iret/=0) then
        print *,'For dimension ',trim(dimname)
        call nchandle_error(ncid, iret)
      end if

      iret = nf90_def_var(ncID, trim(validate(dimname)), nf90_double, inqdimid_nc, VarID)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, 'longname', dimlongname, dimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, 'units', dimunit, dimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)
      iret = putatt_nc(ncID, '_FillValue', fillvalue_double, dimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

     if (.not. ltime) then !Fill the dimension-values for any dimension but time
        ldef = ensuredata_nc(ncid)
        iret = nf90_put_var(ncid, varID, dimvalues)
        if (iret/=0) then
          print *,'For dimension ',trim(dimname)
          call nchandle_error(ncid, iret)
        end if
        ldef = ensuredefine_nc(ncid)
      end if
    end if

  end function inqdimid_nc

!-------------------------------------------------------------------------
!> Write down a number of variables that depend on (possibly) time and 0 other dimension.
  subroutine writevar0D_nc(ncid,ncname,var,nrec)
    integer, intent(in)              :: ncid   !< Netcdf file number
    character(len=*),intent(in)      :: ncname !< The variable names
    real(wp),intent(in)              :: var    !< The variables to be written to file
    integer, intent(inout), optional :: nrec   !< Netcdf Record number

    integer :: iret,varid, udimid, loc(1)
    character(len=20) :: udimname

    iret = nf90_inq_varid(ncid, trim(validate(ncname)), VarID)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

    if (present(nrec)) then
      iret = nf90_inquire(ncid, unlimitedDimId = udimid)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

      iret = nf90_inquire_dimension(ncid,udimid,name=udimname)
      if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

      if(trim(ncname)==trim(udimname)) nrec = nrec+1
      loc = nrec
    else
      loc = 1
    end if

    iret = nf90_put_var(ncid, VarID, var, loc)
    if (iret /= nf90_noerr .and. iret /= nf90_erange) call nchandle_error(ncid, iret)

    iret = sync_nc(ncid)

  end subroutine writevar0D_nc

!-------------------------------------------------------------------------
!> Write down a number of variables that depend on (possibly) time and 1 other dimension.
  subroutine writevar1D_nc(ncid,ncname,var,nrec)
    integer, intent(in)           :: ncid   !< Netcdf file number
    character(len=*),intent(in)   :: ncname !< The variable names
    real(wp),dimension(:),intent(in) :: var    !< The variables to be written to file
    integer, intent(in), optional    :: nrec   !< Netcdf Record number

    integer :: iret,varid
    integer :: dimids(2), dimsize(2)
    INTEGER :: loc(2), nd

    dimsize(2) = 1 !The time dimension has size 1

    iret = nf90_inq_varid(ncid, trim(validate(ncname)), varid)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

    iret = nf90_inquire_variable(ncid,varid,dimids=dimids)
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

    iret = nf90_inquire_dimension(ncid,dimids(1),len=dimsize(1))
    if (iret /= nf90_noerr) call nchandle_error(ncid, iret)

    loc(1) = 1
    IF (PRESENT(nrec)) THEN
      loc(2) = nrec
      nd = 2
    ELSE
      nd = 1
    END IF

    iret = nf90_put_var(ncid, VarID, var(1:dimsize(1)), loc(1:nd), dimsize)

    if (iret /= nf90_noerr .and. iret /= nf90_erange) call nchandle_error(ncid, iret)
    iret = sync_nc(ncid)

  end subroutine writevar1D_nc

!-------------------------------------------------------------------------
  integer function putatt_str_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    character(len=*), intent(in)           :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef

    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_str_nc = nf90_put_att(ncid, varid, validate(attrname), trim(attrval))

    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_str_nc

!-------------------------------------------------------------------------
  integer function putatt_single_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(sp), intent(in)                   :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef

    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_single_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_single_nc

!-------------------------------------------------------------------------
  integer function putatt_double_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(dp), intent(in)                   :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef

    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_double_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_double_nc

!-------------------------------------------------------------------------
  integer function putatt_int_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(i4), intent(in)                :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef

    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_int_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_int_nc

!-------------------------------------------------------------------------
  integer function putatt_short_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(i2), intent(in)                :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid
    logical :: ldef

    ldef = ensuredefine_nc(ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    putatt_short_nc = nf90_put_att(ncid, varid, validate(attrname), attrval)
    if (ldef .eqv. .false.) then
      ldef = ensuredata_nc(ncid)
    end if
  end function putatt_short_nc

!-------------------------------------------------------------------------
  integer function getatt_str_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    character(len=*), intent(out)          :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid

    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_str_nc = nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_str_nc

!-------------------------------------------------------------------------
  integer function getatt_single_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(sp), intent(out)                  :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid

    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_single_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_single_nc

!-------------------------------------------------------------------------
  integer function getatt_double_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    real(dp), intent(out)                  :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid

    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_double_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_double_nc

!-------------------------------------------------------------------------
  integer function getatt_int_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(i4), intent(out)               :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid

    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_int_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_int_nc

!-------------------------------------------------------------------------
  integer function getatt_short_nc(ncid, attrname, attrval, varname)
    integer, intent(in)                    :: ncid
    character(len=*), intent(in)           :: attrname
    integer(i2), intent(out)               :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid

    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_short_nc= nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_short_nc

!-------------------------------------------------------------------------
  integer function getatt_str_fname_nc(fname, attrname, attrval, varname)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: attrname
    character(len=*), intent(out)          :: attrval
    character(len=*), intent(in), optional :: varname
    integer :: iret, varid, ncid

    iret = nf90_open(trim(fname), NF90_NOWRITE, ncid)
    if (present(varname)) then
      iret  = nf90_inq_varid(ncid,validate(varname),varid)
    else
      varid = nf90_global
    end if
    getatt_str_fname_nc = nf90_get_att(ncid, varid, validate(attrname), attrval)
  end function getatt_str_fname_nc

!-------------------------------------------------------------------------
  function validate(input) !\todo Make this allocatable as soon as all common compilers allow it
    character(len=*), intent(in) :: input
    character(len=len_trim(input)):: validate
    character(len=2) :: cforbidden
    integer          :: n, nn
    cforbidden = '()'
    nn = 0
    do n=1,len_trim(input)
      validate(n:n) = ' '
      if (scan(input(n:n),cforbidden) == 0) then
        nn = nn + 1
        validate(nn:nn) = input(n:n)
      end if
    end do
  end function validate

!-------------------------------------------------------------------------
!> Synchronize the netcdf file on disk with the buffer.
  integer function sync_nc(ncid)
    integer, intent(in) :: ncid            !< Netcdf file number
    sync_nc = 1
    if (lsync) sync_nc = nf90_sync(ncid)
  end function sync_nc

!-------------------------------------------------------------------------
!> Write netcdf error to screen and stop the program
  subroutine nchandle_error(ncid, status)
    integer, intent(in) :: ncid
    integer, intent(in) :: status !< NetCDF error code
    integer       :: iret
    character(80) :: fname
    if(status /= nf90_noerr) then
      iret = getatt_nc(ncid, 'title', fname)
      if (iret /= nf90_noerr) fname = ''
      print *, 'NetCDF error in file ' // trim(fname)
      print *, status, trim(nf90_strerror(status))
      call finish('mo_write_netcdf:','stoppping!')
    end if
  end subroutine nchandle_error

!-------------------------------------------------------------------------
end module mo_write_netcdf
