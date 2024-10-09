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

! Basis for RTTOV interface modules

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

#include "mo_rtifc_macros.incf"

  !-------------
  ! Modules used
  !-------------

  use iso_fortran_env,    only: stdout => output_unit, &!
                                stderr => error_unit,  &!
                                iostat_end              !

#if defined(__DACE__) || defined(__ICON__)
  use mo_exception,       only: finish
#endif

#if defined(__DACE__)
  use mo_rad,               only: usd
#endif

#if (_RTTOV_VERSION > 0)
  use rttov_const,        only: version,               &!
                                release,               &!
                                minor_version,         &!
                                fastem_sp,             &! max. number of fastem surface parameters
                                qmin_rttov => qmin,    &! minimum allowed water vapour
                                qmax_rttov => qmax,    &! maximum allowed water vapour
                                tmin_rttov => tmin,    &! minimum allowed temperature
                                tmax_rttov => tmax,    &! maximum allowed temperature
                                min_od,                &!
                                def_gas_unit => gas_unit_specconc,  &
                                rts_land     => surftype_land,      &
                                rts_sea      => surftype_sea,       &
                                rts_ice      => surftype_seaice
  use parkind1,           only: jprb, jpim, jplm, jpis
#endif

#if defined(_RTTOV_GOD)
  use rttov_types,        only: rttov_god_par
  use rttov_god,          only: rttov_o2god
#endif

#if defined(_RTIFC_DISTRIBCOEF)
  use mpi
#endif

  implicit none

  public

  !-----------
  ! Interfaces
  !-----------

  ! MPI routines
#if defined(_RTIFC_DISTRIBCOEF)
#if (_RTTOV_VERSION >= 12) || defined(_RTIFC_USE_MPIF)
  interface p_bcast
#if defined(_RTIFC_USE_MPIF)
    module procedure p_bcast_rttov_integer
    module procedure p_bcast_rttov_integer_1d
    module procedure p_bcast_rttov_integer_2d
    module procedure p_bcast_rttov_integer_3d
    module procedure p_bcast_rttov_integer_4d
    module procedure p_bcast_rttov_real_1d
    module procedure p_bcast_rttov_real_2d
    module procedure p_bcast_rttov_real_3d
    module procedure p_bcast_rttov_real_4d
    module procedure p_bcast_rttov_bool
    module procedure p_bcast_rttov_char_1d
#endif
#if (_RTTOV_VERSION >= 12)
    module procedure p_bcast_rttov_sinteger_4d
    module procedure p_bcast_rttov_sinteger_3d
    module procedure p_bcast_rttov_sinteger_2d
    module procedure p_bcast_rttov_complex_1d
#endif
  end interface
#endif

#endif /* _RTIFC_DISTRIBCOEF */

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
#if (_RTTOV_VERSION <= 0)
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
#else
  ! Used from rttov modules
#endif

#if !defined(__DACE__)
  integer :: usd = stdout
#endif


  ! error codes and messages
  integer,         parameter :: nerr            = 24                        ! number of different error messages
  character(len=100)         :: err_msg(0:nerr)
#define DEF_ERR(codename, code, msg) integer,parameter::codename=code;data err_msg(code)/msg/
  DEF_ERR(NO_ERROR            ,  0, 'No error. Everything okay.')
  DEF_ERR(ERR_ALLOC           ,  1, 'Allocation error.')
  DEF_ERR(ERR_DIM             ,  2, 'Wrong dimension size of input array')
  DEF_ERR(ERR_RTTOV_SETUP     ,  3, 'Error in RTTOV setup')
  DEF_ERR(ERR_CLOUD_AERO_MISSM,  4, 'Cloud/Aerosol class mismatch')
  DEF_ERR(ERR_RTTOV_CALL      ,  5, 'RTTOV call failed')
  DEF_ERR(WARN_RTTOV          ,  6, 'Warning: RTTOV error status /= 0')
  DEF_ERR(ERR_RTTOV_MPI       ,  7, 'MPI error')
  DEF_ERR(ERR_NO_RTTOV_LIB    ,  8, 'No RTTOV library available')
  DEF_ERR(ERR_RTTOV_PREC      ,  9, 'Mismatch of working precision and RTTOV precision')
  DEF_ERR(ERR_CLOUD_INCONS    , 10, 'Cloud cover and cloud water content arrays inconsistent')
  DEF_ERR(ERR_GOD_FILE        , 11, 'Failed to read god_par_file')
  DEF_ERR(ERR_WR_PROF         , 12, 'Failed to write hdf5 profile file')
  DEF_ERR(ERR_INVALID_TSKIN   , 13, 'Some invalid t_surf/T_G/ts_fg')
  DEF_ERR(ERR_INPUT           , 14, 'Invalid/unsupported input')
  DEF_ERR(ERR_NO_ATLAS        , 15, 'No atlas support in current configuration')
  DEF_ERR(ERR_ATLAS_INIT      , 16, 'Atlas was not initialized')
  DEF_ERR(ERR_INVALID_INSTR   , 17, 'Invalid instrument (e.g. not supported by atlas)')
  DEF_ERR(ERR_TRACEGAS_INCONS , 18, 'Trace gase options inconsistent with input')
  DEF_ERR(ERR_INVALID_NLEVS   , 19, 'Invalid number of levels')
  DEF_ERR(ERR_INVALID_VERSION , 20, 'Invalid RTTOV version')
  DEF_ERR(ERR_PRECISION_INCONS, 21, 'Real data precisions inconsistent (wp /= jprb).')
  DEF_ERR(ERR_NO_COEFS        , 22, 'No matching coefs found.')
  DEF_ERR(ERR_MULTIPLE_COEFS  , 23, 'Multiple matching coefs found.')
  DEF_ERR(ERR_NO_OPTS_TMPL    , 24, 'Option templates not initialized so far')
! DEF_ERR(ERR_                , XX, '')
#undef DEF_ERR



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
#if (_RTTOV_VERSION >= 13)
  integer                    :: default_icede_param       =       4     ! Same as default_idg, but for RTTOVv13 !CStrtifc_set_opt_template
#else
  integer                    :: default_idg               =       4     ! Scheme for IWC to eff
                                                                        ! shape of the ice crystals RTTOVv12
#endif
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

#if (_RTTOV_VERSION <= 0)
  real(kind=wp)              :: min_od                    = -1._wp
#else
  ! Used from rttov_const
#endif

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
#if (_RTTOV_VERSION > 0)
    if (wp /= jprb) then
      status = ERR_PRECISION_INCONS
      return
    end if
#endif
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
#if (_RTTOV_VERSION <= 0)
    vers = 'NO_RTTOV'
#else
    write(vers,'("RTTOV",I2.2,".",I1,".",I1)') version, release, minor_version
#endif
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


#if ! (defined(__DACE__) || defined(__ICON__))
  subroutine finish(proc, msg)
    character(len=*), intent(in) :: proc
    character(len=*), intent(in) :: msg

    write(stderr,'(/,80("*"),/)')
    write(stderr,'(2x,A,": ",A)') trim(proc),trim(msg)
    write(stderr,'(/,80("*"),/)')
    STOP 13

  end subroutine finish
#endif


#if (_RTTOV_VERSION > 0)

  subroutine get_weighting_function(p, t, transm, wf, height)
    real(kind=jprb)           :: p(:)
    real(kind=jprb)           :: t(:)
    real(kind=jprb)           :: transm(:)
    real(kind=jprb), optional :: wf(:)
    real(kind=jprb), optional :: height

    integer :: wf_method = 1  ! 0: Method proposed by Lucio Torrisi
                              ! 1: Leapfrog-like method
                              ! 2: spline interpolation
    integer :: hgt_method = 0 ! 0: Method proposed by Lucio Torrisi
                              ! 1: Maximum of weighting function
                              ! 2: spline based maximum of weighting function
    real(kind=jprb) :: wf_(size(p))
    real(kind=jprb) :: p_, t_, psum, wsum, max_wf
    integer         :: i, i1, i2, i3, nlev

    nlev = size(p)

    select case(wf_method)
    case(0)
      do i = 1, nlev
        if (i==1) then
          i1 = 1
          i2 = 2
          i3 = i1
        else
          i1 = i - 1
          i2 = i
          i3 = i2
        end if
        ! Hydrostatic assumption
        wf_(i) = - g * p(i3) / (rd * t(i3)) * &
             (transm(i1) - transm(i2)) / (p(i1) - p(i2))
      end do
    case(1)
      do i = 1, nlev
        if (i==1) then
          i1 = 1
          i2 = 2
          p_ = (p(i1) * p(i2)) * 0.5
          t_ = (t(i1) * t(i2)) * 0.5
        elseif (i==nlev) then
          i1 = i - 1
          i2 = i
          p_ = (p(i1) * p(i2)) * 0.5
          t_ = (t(i1) * t(i2)) * 0.5
        else
          i1 = i - 1
          i2 = i + 1
          p_ = p(i)
          t_ = t(i)
        end if
        ! Hydrostatic assumption
        wf_(i) = - g * p_ / (rd * t_) * &
             (transm(i1) - transm(i2)) / (p(i1) - p(i2))
      end do
    case(2)
    end select

    if (present(wf)) then
      wf(1:nlev) = wf_(1:nlev)
    end if

    if (present(height)) then
      select case(hgt_method)
      case(0)
        max_wf = maxval(wf_(:))
        psum = 0._jprb
        wsum = 0._jprb
        do i = 1, nlev
          if (wf_(i) > 0.8 * max_wf) then
            psum = psum + p(i) * wf_(i)
            wsum = wsum + wf_(i)
          end if
        end do
        if (wsum /= 0._jprb) then
          height = psum / wsum
        else
          height = -999._jprb
        end if
      case(1)
        i1 = maxloc(wf_(:), 1)
        if (i1 > 1 .and. i1 < nlev) then
          psum = 0._jprb
          wsum = 0._jprb
          do i = i1 - 1, i1 + 1
            psum = psum + p(i) * wf_(i)
            wsum = wsum + wf_(i)
          end do
          if (wsum /= 0._jprb) then
            height = psum / wsum
          else
            height = -999._jprb
          end if
        else
          height = p(i1)
        end if
      case(2)
        !> \todo spline-based weighting function
      end select
    end if

  end subroutine get_weighting_function
#endif

#if defined(_RTTOV_GOD)
  subroutine check_god_infl(god, istat, msg, debug, tau, od_ref)
    type(rttov_god_par), intent(in)           :: god(:)
    integer,             intent(out)          :: istat
    character(len=*),    intent(out)          :: msg
    real(kind=jprb),     intent(in), optional :: tau(:)
    real(kind=jprb),     intent(in), optional :: od_ref(:)
    logical,             intent(in), optional :: debug
    ! Calculate transmission and optical depth profile (without god smoothing) on the
    ! basis on god-smoothed profiles.
    character(len=14), parameter   :: proc = 'check_god_infl'
    real(kind=jprb),   allocatable :: tau0 (:)
    real(kind=jprb),   allocatable :: tau_ (:)
    real(kind=jprb),   allocatable :: infl (:)
    real(kind=jprb),   allocatable :: infl0(:)
    real(kind=jprb) :: s_infl, od, od0, od0_tot, od_tot, rdiff
    integer :: nl, i, l
    logical :: l_tau
    logical :: l_debug = .false.

    if (.not.present(tau) .and. .not.present(od_ref)) &
         call finish(proc, 'Either tau or od_ref must be present')
    l_tau = .not.present(od_ref)

    if (present(debug)) then
      l_debug = debug
    else
      l_debug = .false.
    end if

    istat = 0
    msg   = ''

    if (l_tau) then
      nl = size(tau)-1
    else
      nl = size(od_ref)
    end if
    allocate(tau_(nl+1), tau0(nl+1), infl(nl), infl0(nl))

    if (l_tau) then
      tau_(:) = tau(:)
      od0_tot = 0._jprb
      do i = 1,nl
        if (tau_(i) > tau_(i+1)) then
          od  = log(tau_(i+1)/tau_(i))
          od0 = min(rttov_o2god(od, p=god(i)), 0._jprb)
        else
          od  = 0._jprb
          od0 = 0._jprb
        end if
        tau0(i) = exp(od0_tot)
        od0_tot = od0_tot + od0
        if (l_debug) write(usd,*) 'check_god_infl tau -> tau0',i,od,od0, tau_(i), tau0(i)
      end do
      tau0(nl+1) = exp(od0_tot)
    else
      nl = size(od_ref)
      tau_(1) = 1._jprb
      tau0(1) = 1._jprb
      od0_tot = 0._jprb
      od_tot  = 0._jprb
      do i = 1, nl
        od  = -abs(od_ref(i))
        od_tot = od_tot + od
        tau_(i+1) = exp(od_tot)
        od0 = min(rttov_o2god(od, p=god(i)), 0._jprb)
        od0_tot = od0_tot + od0
        tau0(i+1) = exp(od0_tot)
        if (l_debug) write(usd,*) 'check_god_infl od_ref -> tau*',i,od,od0, tau_(i), tau0(i)
      end do
    end if

    s_infl = 0._jprb
    do i = 1,nl
      infl (i) = tau_(i) * (tau_(i) - tau_(i+1))
      infl0(i) = tau0(i) * (tau0(i) - tau0(i+1))
      s_infl = s_infl + infl0(i)
    end do

    if (l_debug) write(usd,*) i,s_infl
    if (s_infl > 0._jprb) then
      do i = 1, nl
        rdiff = (infl(i)-infl0(i))/s_infl
        if (l_debug) write(usd,*) i,infl(i),infl0(i),infl(i)-infl0(i),(infl(i)-infl0(i))/s_infl
        if (rdiff > god_thresh) then
          istat = 1
          l = len_trim(msg) + 1
          if (l+30 <= len(msg)) then
            if (l > 1) then
              msg = trim(msg)//','
              l = l + 1
            end if
            write(msg(l:),'("layer ",I3," d(infl)/infl=",f6.4)') i, rdiff
          end if
        end if
      end do
    else
      istat = 2
      write(msg,*) 'Unable to calculate influence, total influence==0 !!!'
    end if

  end subroutine check_god_infl

  subroutine l2c_god_tau(tau, god, l2c, l2c_corr, debug)
    real(kind=jprb),     intent(in)           :: tau(:)
    type(rttov_god_par), intent(in)           :: god(:)
    real,                intent(in)           :: l2c
    real,                intent(out)          :: l2c_corr
    logical,             intent(in), optional :: debug
    ! Calculate transmission and optical depth profile (without god smoothing) on the
    ! basis on god-smoothed profiles. Then determine the transmission (with god-smoothing)
    ! at the level "l2c". Finally, calculate a "corrected l2c" as the level where the transmission
    ! is obtained without god-smoothing.
    character(len=11), parameter :: proc = 'l2c_god_tau'
    real(kind=jprb),  parameter :: eps = 1.e-30
    real(kind=jprb)             :: tau0(size(tau))
    real(kind=jprb)             :: od, od0, od0_tot, tr_l2c
    real                        :: w0
    integer                     :: nl, i, l0, l1, igod
    logical                     :: l_debug = .false.

    if (present(debug)) then
      l_debug = debug
    else
      l_debug = .false.
    end if
    nl = size(tau)

    ! Calculate transmission profile without god-smoothing
    od0_tot = log(tau(1))
    do i = 1, nl-1
      igod = i + nlevs_top
      if (tau(i) > tau(i+1) .and. tau(i+1) > eps) then
        od  = log(tau(i+1)/tau(i))
        od0 = min(rttov_o2god(od, p=god(igod)), 0._jprb)
      else
        od  = 0._jprb
        od0 = 0._jprb
      end if
      tau0(i) = exp(od0_tot)
      od0_tot = od0_tot + od0
      if (l_debug) write(usd,*) proc//' tau',i, tau(i), tau0(i),od,od0,&
           god(igod)%tr(1:god(igod)%ntr)%opdep,rttov_o2god(od, p=god(igod))

    end do
    tau0(nl) = exp(od0_tot)

    ! Determine transmission at l2c
    l0 = int(l2c)
    l1 = l0 + 1
    if (l_debug) write(usd,*) proc//' l2c',l2c,l0,l1,nl
    if (l1 <= nl) then
      w0 = (l1 - l2c)
      tr_l2c = w0 * tau(l0) + (1._jprb - w0) * tau(l1)
    else
      l2c_corr = l2c
      if (l_debug) write(usd,*) proc//' l1 > nl',l2c_corr
      return
    end if
    if (l_debug) write(usd,*) proc//' tr_l2c',l2c,tr_l2c

    ! Calc. level where tau0 equals tr_l2c
    l2c_corr = -1.
    do i = 1, nl
      if (tau0(i) <= tr_l2c) then
        if (i >  1) then
          w0 = (tau0(i) - tr_l2c)/(tau0(i) - tau0(i-1))
          l2c_corr = (i-1) * w0 + (1.-w0) * i
          exit
        else
          l2c_corr = 1.
        end if
        exit
      end if
    end do
    if (l2c_corr <= 0.) l2c_corr = real(nl)

    if (l_debug) write(usd,*) proc//' l2c_corr',l2c_corr,l2c_corr-l2c,l2c_corr==l2c

  end subroutine l2c_god_tau

  subroutine l2c_god_opd(opd, god, lp_l2c, lp_rt, l2c, l2c_corr, debug)
    real(kind=jprb),     intent(in)           :: opd(:)
    type(rttov_god_par), intent(in)           :: god(:)
    real(kind=jprb),     intent(in)           :: lp_l2c(:)
    real(kind=jprb),     intent(in)           :: lp_rt(:)
    real,                intent(in)           :: l2c
    real,                intent(out)          :: l2c_corr
    logical,             intent(in), optional :: debug
    ! Here the optical depths are on RTTOV-levels - not on user levels as in l2c_god_tau.
    ! Determine the log(p) value "lp" corresponding to l2c. Calculate (total) optical depth
    ! profile od/od0 (with/without god smoothing) on the RTTOV levels. Find the od at lp,
    ! and the value "lp_corr" where the od == od0. Finally determine the level corresponding to
    ! lp_corr.
    character(len=11), parameter :: proc = 'l2c_god_opd'
    real(kind=jprb),   parameter :: eps = 1.e-10_jprb
    real(kind=jprb)              :: od(size(opd)+1),od0(size(opd)+1)
    real(kind=jprb)              :: lp, lp_corr, od_l2c
    real(kind=jprb)              :: od_, od0_
    real(kind=jprb)              :: w0
    integer                      :: nl, nlr, i, l0, l1, igod
    logical                      :: l_debug = .false.

    l2c_corr = l2c

    if (present(debug)) then
      l_debug = debug
    else
      l_debug = .false.
    end if
    nlr = size(opd)+1
    nl  = size(lp_l2c)
    if (l_debug) write(usd,*) proc//' nl',nl,nlr

    ! Calculate opdep (total( profile with/without god-smoothing
    od (1) = 0._jprb
    od0(1) = 0._jprb
    do i = 2, nlr
      igod = i-1 !+ nlevs_top
      od_    = -abs(opd(i-1))
      od(i)  = od(i-1) + od_
      od0_   = min(rttov_o2god(od_, p=god(igod)), 0._jprb)
      od0(i) = od0(i-1) + od0_
      if (l_debug) write(usd,*) proc//' od',i, od_, od(i), od0_, od0(i), &
           god(igod)%tr(1:god(igod)%ntr)%opdep,rttov_o2god(od_, p=god(igod))
    end do
    if (l_debug) write(usd,*) proc//' od_diff',abs(od(nlr) - od0(nlr)), eps
    if (abs(od(nlr) - od0(nlr)) < eps) return

    ! Determine ln(p) at l2c
    l0 = int(l2c)
    l1 = l0 + 1
    if (l_debug) write(usd,*) proc//' l2c',l2c,l0,l1,nl
    if (l1 <= nl) then
      w0 = (l1 - l2c)
      lp = w0 * lp_l2c(l0) + (1._jprb - w0) * lp_l2c(l1)
    else
      if (l_debug) write(usd,*) proc//' l1 > nl',l2c_corr
      return
    end if
    if (l_debug) write(usd,*) proc//' lp',l2c,lp,exp(lp)

    ! Find corresponding RT-level
    l1 = -1
    do i = 1, nlr
      if (lp_rt(i) >= lp) then
        l1 = i
        exit
      end if
    end do
    if (l1 <= 1) then
      if (l_debug) write(usd,*) proc//' l1 <= 0',l2c_corr
      return
    end if
    l0 = l1 - 1
    if (l_debug) write(usd,*) proc//' l0,l1',l0,l1,abs(od(l1) - od0(l1)),eps
    if (abs(od(l1) - od0(l1)) < eps) return

    ! Calc. od at l2c
    w0 = (lp_rt(l1) - lp)/(lp_rt(l1) - lp_rt(l0))
    if (l_debug) write(usd,*) proc//' calc_od_l2c',lp_rt(l0),lp_rt(l1),w0,od(l0),od(l1)
    od_l2c =  w0 * od(l0) + (1._jprb - w0) * od(l1)
    if (l_debug) write(usd,*) proc//' od_l2c',od_l2c

    ! Find RT level with od0 == od_l2c
    l1 = nlr
    do i = 1, nlr
      if (od0(i) <= od_l2c) then
        l1 = i
        exit
      end if
    end do
    ! In case of zero optical depths, i.e. constant od0 over several levels, we select
    ! the most "conservative" solution, i.e. the lowest level:
    i = l1
    do while (l1 < nlr)
      if (od0(i) == od0(l1+1)) then
        l1 = l1 + 1
      else
        exit
      end if
    end do
    l0 = l1 - 1
    if (l_debug) write(usd,*) proc//' l0,l1',l0,l1

    ! Calc. p at this level
    lp_corr = lp_rt(l1)
    if (l1 > 1 .and. od0(l1) <= od_l2c) then
      if (od0(l1) < od0(l0)) then
        w0 = (od0(l1) - od_l2c) / (od0(l1) - od0(l0))
        if (l_debug) write(usd,*) proc//' calc lp_corr',od0(l0),od0(l1),w0,lp_rt(l0),lp_rt(l1)
        lp_corr = w0 * lp_rt(l0) + (1._jprb - w0) * lp_rt(l1)
      end if
    end if
    if (l_debug) write(usd,*) proc//' lp_corr',lp_corr,exp(lp_corr)
    lp_corr = max(lp, lp_corr)

    ! Corrected l2c value
    l1 = nl
    do i = 1, nl
      if (l_debug) write(usd,*) proc//' lp_l2c',i,lp_l2c(i),exp(lp_l2c(i))
      if (lp_l2c(i) >= lp_corr) then
        l1 = i
        exit
      end if
    end do
    l0 = l1 - 1
    if (l_debug) write(usd,*) proc//' l0,l1',l0,l1
    l2c_corr = real(l1)
    if (l1 > 1 .and. lp_l2c(l1) >= lp_corr) then
      w0 = (lp_l2c(l1) - lp_corr) / (lp_l2c(l1) - lp_l2c(l0))
      l2c_corr = w0 * l0 + (1._jprb - w0) * l1
    end if
    if (l_debug) write(usd,*) proc//' l2c_corr',l2c_corr,l2c_corr-l2c,l2c_corr==l2c

  end subroutine l2c_god_opd
#endif

!==============================
#if defined(_RTIFC_DISTRIBCOEF)
!==============================

#if defined(_RTIFC_USE_MPIF)
  subroutine p_bcast_rttov_integer(buffer,source,comm)
    integer(jpim),    intent(inout)   :: buffer
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------
    ! Broadcast an integer vector across all available processors
    !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,1, MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_integer@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_integer


  subroutine p_bcast_rttov_integer_1d(buffer,source,comm)
    integer(jpim),    intent(inout)   :: buffer(:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------
    ! Broadcast an integer vector across all available processors
    !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_integer_1d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_integer_1d


  subroutine p_bcast_rttov_integer_2d(buffer,source,comm)
    integer(jpim),    intent(inout)   :: buffer(:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a real rank-2 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_integer_2d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_integer_2d


  subroutine p_bcast_rttov_integer_3d(buffer,source,comm)
    integer(jpim),    intent(inout)   :: buffer(:,:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a real rank-2 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_integer_3d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_integer_3d


  subroutine p_bcast_rttov_integer_4d(buffer,source,comm)
    integer(jpim),    intent(inout)   :: buffer(:,:,:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a real rank-2 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_INTEGER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_integer_4d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_integer_4d


  subroutine p_bcast_rttov_real_1d(buffer,source,comm)
    real(jprb),    intent(inout)   :: buffer(:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------
    ! Broadcast a real vector across all available processors
    !--------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_real_1d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_real_1d


  subroutine p_bcast_rttov_real_2d(buffer,source,comm)
    real(jprb),    intent(inout)   :: buffer(:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a real rank-2 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_real_2d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_real_2d


  subroutine p_bcast_rttov_real_3d(buffer,source,comm)
    real(jprb),    intent(inout)   :: buffer(:,:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a real rank-3 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_real_3d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_real_3d


  subroutine p_bcast_rttov_real_4d(buffer,source,comm)
    real(jprb),    intent(inout)   :: buffer(:,:,:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a real rank-3 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_PRECISION, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_real_4d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_real_4d


  subroutine p_bcast_rttov_bool(buffer,source,comm)
    logical,          intent(inout)   :: buffer
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------
    ! Broadcast an integer vector across all available processors
    !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,1, MPI_LOGICAL, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_bool@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_bool


  subroutine p_bcast_rttov_char_1d(buffer,source,comm)
    character(len=4), intent(inout)   :: buffer(:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------
    ! Broadcast an integer vector across all available processors
    !------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_CHARACTER, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_char_1d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_char_1d

#endif /* _RTIFC_USE_MPIF */

#if (_RTTOV_VERSION >= 12)
  subroutine p_bcast_rttov_sinteger_4d(buffer,source,comm)
    integer(jpis),    intent(inout)   :: buffer(:,:,:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a small integer rank-4 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_sinteger_4d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_sinteger_4d


  subroutine p_bcast_rttov_sinteger_3d(buffer,source,comm)
    integer(jpis),    intent(inout)   :: buffer(:,:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a small integer rank-4 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_sinteger_3d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_sinteger_3d


  subroutine p_bcast_rttov_sinteger_2d(buffer,source,comm)
    integer(jpis),    intent(inout)   :: buffer(:,:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !--------------------------------------------------------------
    ! Broadcast a small integer rank-4 array across all available processors
    !--------------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_sinteger_2d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_sinteger_2d


  subroutine p_bcast_rttov_complex_1d(buffer,source,comm)
    complex(jprb),    intent(inout)   :: buffer(:)
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !-----------------------------------------------------------
    ! Broadcast a complex vector across all available processors
    !-----------------------------------------------------------
    integer :: lcom, errorcode

    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(buffer), MPI_DOUBLE_COMPLEX, source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_complex_1d@mo_rtifc_base', 'MPI ERROR in MPI_Bcast')

  end subroutine p_bcast_rttov_complex_1d
#endif


!==============================
#endif /* _RTIFC_DISTRIBCOEF */
!==============================

end module mo_rtifc_base
