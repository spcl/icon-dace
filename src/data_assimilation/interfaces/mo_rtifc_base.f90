! Basis for RTTOV interface modules
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

MODULE mo_rtifc_base

!---------------------------------------------------------------------------
!
! Description:
!   This module contains basic routines, parameters and variables.
!   for the interface modules. This might be stuff, that does not
!   depend on the RTTOV version, or DWD specific features, that will
!   not change (quite likely) in future RTTOV versions.
!---------------------------------------------------------------------------

!---------------------
! MACRO SETTINGS BEGIN
!---------------------

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
!+ Preprocessor macros for mo_rtifc_*.f90
!
!
! Description:
!   This file contains the Macro definitions
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email:  robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! @VERSION@    @DATE@     Robin Faulwetter
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Heritage: mo_rtifc_10.f90
!=======================================================================


!---------------
! MACRO SETTINGS
!---------------

! Take into account _DACE_ and __DACE__

! Icon without RTTOV

! Set macros for RTTOV coefficient distribution

! Select mpi routines to be used for coefficient distribution





!-------------------------
!--------------------------------

  !-------------
  ! Modules used
  !-------------

  use iso_fortran_env,    only: stdout => output_unit, &!
                                stderr => error_unit,  &!
                                iostat_end              !






  implicit none

  public

  !-----------
  ! Interfaces
  !-----------

  ! MPI routines



  !-----------
  ! Parameters
  !-----------

  ! data type precision
  integer,         parameter :: dp              = selected_real_kind(13)     ! Double Precision
  !integer,        parameter :: sp              = selected_real_kind(6)      ! Single Precision
  integer,         parameter :: wp              = dp                         ! Working precision

  ! Values for verbosity
  integer,         parameter :: silent          =  0
  integer,         parameter :: production      =  1
  integer,         parameter :: verbose         =  2
  integer,         parameter :: debug           =  3

  ! Flag variables for output arrays
  integer,         parameter :: OUT_ASB         = 0                          ! all sky brightness temp.
  integer,         parameter :: OUT_CSB         = 1                          ! clear sky brightness temp.
  integer,         parameter :: OUT_ASR         = 2                          ! all sky radiances
  integer,         parameter :: OUT_CSR         = 3                          ! clear sky radiances
  integer,         parameter :: OUT_VIS         = 4                          ! reflectances (all sky and clear sky)

  ! Physical parameters required for weighting function calculation
  real(kind=wp),   parameter :: rd              = 287.05_wp
  real(kind=wp),   parameter :: g               = 9.80665_wp

  ! Possible coefficient level numbers
  integer,         parameter :: levels_rttov(4) = (/ 44, 51, 54, 101 /)

  ! Max numbers
  integer,         parameter :: mx_chan         = 8700

  ! Define variables for no-rttov case, that are used from RTTOV modules otherwise

!  integer,         parameter :: jpim            = selected_int_kind(9)       ! standard integer type
!  integer,         parameter :: jprb            = selected_real_kind(13,300) ! standard real type
!  integer,         parameter :: jplm            = kind(.true.)               ! standard logical type
  integer,         parameter :: version         = 0                          ! RTTOV version
  integer,         parameter :: fastem_sp       = 5                          ! Number of FASTEM parameters
  integer,         parameter :: def_gas_unit    = 0
  real(kind=wp),   parameter :: qmin_rttov      = -1._wp
  real(kind=wp),   parameter :: qmax_rttov      = -1._wp
  real(kind=wp),   parameter :: tmin_rttov      = -1._wp
  real(kind=wp),   parameter :: tmax_rttov      = -1._wp
  integer,         parameter :: rts_land        = 0
  integer,         parameter :: rts_sea         = 1
  integer,         parameter :: rts_ice         = 2





  integer :: usd = stdout



  ! error codes and messages
  integer,         parameter :: nerr            = 24                        ! number of different error messages
  character(len=100)         :: err_msg(0:nerr)

  integer,parameter::NO_ERROR            =  0;data err_msg(  0)/ 'No error. Everything okay.'/
  integer,parameter::ERR_ALLOC           =  1;data err_msg(  1)/ 'Allocation error.'/
  integer,parameter::ERR_DIM             =  2;data err_msg(  2)/ 'Wrong dimension size of input array'/
  integer,parameter::ERR_RTTOV_SETUP     =  3;data err_msg(  3)/ 'Error in RTTOV setup'/
  integer,parameter::ERR_CLOUD_AERO_MISSM=  4;data err_msg(  4)/ 'Cloud/Aerosol class mismatch'/
  integer,parameter::ERR_RTTOV_CALL      =  5;data err_msg(  5)/ 'RTTOV call failed'/
  integer,parameter::WARN_RTTOV          =  6;data err_msg(  6)/ 'Warning: RTTOV error status /= 0'/
  integer,parameter::ERR_RTTOV_MPI       =  7;data err_msg(  7)/ 'MPI error'/
  integer,parameter::ERR_NO_RTTOV_LIB    =  8;data err_msg(  8)/ 'No RTTOV library available'/
  integer,parameter::ERR_RTTOV_PREC      =  9;data err_msg(  9)/ 'Mismatch of working precision and RTTOV precision'/
  integer,parameter::ERR_CLOUD_INCONS    = 10;data err_msg( 10)/ 'Cloud cover and cloud water content arrays inconsistent'/
  integer,parameter::ERR_GOD_FILE        = 11;data err_msg( 11)/ 'Failed to read god_par_file'/
  integer,parameter::ERR_WR_PROF         = 12;data err_msg( 12)/ 'Failed to write hdf5 profile file'/
  integer,parameter::ERR_INVALID_TSKIN   = 13;data err_msg( 13)/ 'Some invalid t_surf/T_G/ts_fg'/
  integer,parameter::ERR_INPUT           = 14;data err_msg( 14)/ 'Invalid/unsupported input'/
  integer,parameter::ERR_NO_ATLAS        = 15;data err_msg( 15)/ 'No atlas support in current configuration'/
  integer,parameter::ERR_ATLAS_INIT      = 16;data err_msg( 16)/ 'Atlas was not initialized'/
  integer,parameter::ERR_INVALID_INSTR   = 17;data err_msg( 17)/ 'Invalid instrument (e.g. not supported by atlas)'/
  integer,parameter::ERR_TRACEGAS_INCONS = 18;data err_msg( 18)/ 'Trace gase options inconsistent with input'/
  integer,parameter::ERR_INVALID_NLEVS   = 19;data err_msg( 19)/ 'Invalid number of levels'/
  integer,parameter::ERR_INVALID_VERSION = 20;data err_msg( 20)/ 'Invalid RTTOV version'/
  integer,parameter::ERR_PRECISION_INCONS= 21;data err_msg( 21)/ 'Real data precisions inconsistent (wp /= jprb).'/
  integer,parameter::ERR_NO_COEFS        = 22;data err_msg( 22)/ 'No matching coefs found.'/
  integer,parameter::ERR_MULTIPLE_COEFS  = 23;data err_msg( 23)/ 'Multiple matching coefs found.'/
  integer,parameter::ERR_NO_OPTS_TMPL    = 24;data err_msg( 24)/ 'Option templates not initialized so far'/
! integer,parameter::ERR_                = XX;data err_msg( XX)/ ''/




  !--------------------------
  ! Internal module variables
  !--------------------------
  ! Only for internal use
  integer                    :: pe_ifc                    = -1           ! mpi id of this processor
  ! For internal use, might be used by the user, but MUST not be modified by user!
  integer                    :: nlevs_top                 = 0            ! Number of coeff. levels above user levels (0 or 1)

  !------------------------------------
  ! External module variables
  ! Intended to be modified by the user
  !------------------------------------
  ! default profile values
  real(kind=wp)              :: default_wfetch            =  100000._wp ! wind fetch
  real(kind=wp)              :: default_fastem(fastem_sp) = (/3.0_wp,5.0_wp,15.0_wp,0.1_wp,0.3_wp/)
                                                                        ! fastem coefficients relevant for land/ice
  integer                    :: default_watertype         =       1     ! water type (fresh 0/ocean 1)
  real(kind=wp)              :: default_salinity          =       0._wp ! salinity
  real(kind=wp)              :: default_o3_surf           = 0.031438_wp ! o3 surface
  real(kind=wp)              :: default_satazim           =      0.0_wp ! satellite azimuth angle
  real(kind=wp)              :: default_sunzenangle       =      0.0_wp ! solar zenith angle
  real(kind=wp)              :: default_sunazangle        =      0.0_wp ! solar azimuth angle
  real(kind=wp)              :: default_ctp               =    500.0_wp ! cloud top pressure
  real(kind=wp)              :: default_cfraction         =      0.0_wp ! cloud fraction



  integer                    :: default_idg               =       4     ! Scheme for IWC to eff
                                                                        ! shape of the ice crystals RTTOVv12

  integer                    :: default_ice_scheme        =       1     ! ice particle scheme
  integer                    :: default_clw_scheme        =       2     ! cloud liquid water scheme
  integer                    :: default_gas_units         = def_gas_unit! default gas unit

  integer                    :: verbosity                 = production
  logical                    :: read1pe                   = .false.      ! Read coeffs.only on I/O PE
                                                  ! (Only effective with -D_RTIFC_DISTRIBCOEF)
  integer                    :: rtifc_alloc_mode          = 0

  ! T/q hard limits
  real(kind=wp)              :: qmin_ifc                  = qmin_rttov
  real(kind=wp)              :: qmax_ifc                  = qmax_rttov
  real(kind=wp)              :: tmin_ifc                  = tmin_rttov
  real(kind=wp)              :: tmax_ifc                  = tmax_rttov


  real(kind=wp)              :: min_od                    = -1._wp




  ! check on regularization limits
  real(kind=wp)              :: chk_plim_t                = -1._wp      ! do not check t for p < chk_plim_t
  real(kind=wp)              :: chk_plim_q                = -1._wp      ! do not check q for p < chk_plim_t
  integer                    :: chk_reg_lims              = 0           ! Check regularization limits in rtifc
                                                                        ! bit1 (1): print results (invalid profiles)
                                                                        ! bit2 (2): set flag for use in calling prog

  ! generalized optical depth
  character(len=300)         :: god_par_file              = ''
  logical                    :: wr_god                    = .false.
  character(len=300)         :: out_path                  = ''
  ! check influence of god smoothing
  real(kind=wp)              :: god_thresh                = 1._wp
  integer                    :: chk_god                   = 0           ! Check influence of god smoothing
                                                                        ! bit1 (1): print results (invalid profiles)
                                                                        ! bit2 (2): set flag for use in calling prog

  ! Atlas
  logical                    :: atlas_single_inst         = .false.

contains


  subroutine rtifc_check_config(vers, nlevs, status)
    integer, intent(in)  :: vers
    integer, intent(in)  :: nlevs
    integer, intent(out) :: status
    !--------------------------------------------
    ! Check RTTOV version and number of levels.
    ! Additionally set nlevs_top (in check_nlevs)
    !--------------------------------------------

    call check_version(vers, status)
    if (status /= NO_ERROR) return
    call check_nlevs(nlevs, nlevs_top, status)
    if (status /= NO_ERROR) return






  end subroutine rtifc_check_config


  subroutine check_version(vers, status)
    integer, intent(in)  :: vers
    integer, intent(out) :: status

    if (vers == version) then
      status = NO_ERROR
    else
      status = ERR_INVALID_VERSION
    end if

  end subroutine check_version


  subroutine check_nlevs(nlevs, nlevs_top, status)
    integer, intent(in)  :: nlevs
    integer, intent(out) :: nlevs_top
    integer, intent(out) :: status

    integer :: i

    status = ERR_INVALID_NLEVS

    do i = 1, size(levels_rttov)
      nlevs_top = (levels_rttov(i) - nlevs)
      if (any(nlevs_top == (/0, 1/))) then
        status = NO_ERROR
        return
      end if
    end do

  end subroutine check_nlevs


  function rttov_version() result(vers)
    character(len=11) :: vers

    vers = 'NO_RTTOV'



  end function rttov_version


  function errmsg(code) result(msg)
    character(len=120)            :: msg
    integer, intent(in)           :: code

    write(msg, '("ERROR (",I4,"):")') code
    if (abs(code) >= lbound(err_msg,1) .and. abs(code) <= ubound(err_msg,1)) then
      msg = trim(msg)//' '//trim(err_msg(abs(code)))
    else
      msg = trim(msg)//' '//'Unknown error'
    end if

  end function errmsg

  function rts_name(styp) result(name)
    character(len=6) :: name
    integer, intent(in) :: styp
    select case(styp)
    case(rts_land)
      name = 'land'
    case(rts_sea )
      name = 'sea'
    case(rts_ice)
      name = 'seaice'
    case default
      name = '??????'
    end select
  end function rts_name



  subroutine finish(proc, msg)
    character(len=*), intent(in) :: proc
    character(len=*), intent(in) :: msg

    write(stderr,'(/,80("*"),/)')
    write(stderr,'(2x,A,": ",A)') trim(proc),trim(msg)
    write(stderr,'(/,80("*"),/)')
    STOP 13

  end subroutine finish





!==============================
!==============================

end module mo_rtifc_base
