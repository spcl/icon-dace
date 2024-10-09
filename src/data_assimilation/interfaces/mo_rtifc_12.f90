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

! Interface for the RTTOV12 library

MODULE mo_rtifc_12

!----------------------------------------------------------------------------
!
! Description:
!   This module contains the RTTOV12 specific stuff for the mo_rtifc module
!
! Heritage: mo_rtifc_10.f90
!----------------------------------------------------------------------------

#include "mo_rtifc_macros.incf"


!=========================
#if (_RTTOV_VERSION == 12)
!=========================

! Define macro for distributed printout
!#if defined(_DACE_) && !defined(NOMPI)
!#define WR call add_line(trim(msg))
!#else
#define WR write(*,*) trim(msg)
!#endif

  !-------------
  ! modules used
  !-------------
  use mo_rtifc_base

#if defined(_DACE_)
  use mo_rad,               only: t_radv                  ! derived type to store radiance obs.
#endif

  use parkind1,             only: jpim,                  &! default integer;
                                  jpis,                  &! small integer
                                  jprb,                  &! default double precision (hopefully)
                                  jplm                    ! default logical type

  use rttov_const,          only: errorstatus_success,   &! RTTOV errorstatus (0)
                                  errorstatus_fatal,     &! RTTOV errorstatus (2)
                                  gas_id_mixed,          &!
                                  ngases_max,            &!
                                  ncldtyp,               &!
                                  baran_ngauss,          &!
                                  gas_unit_specconc,     &! specific concentration (kg/kg over wet air)
                                  gas_unit_ppmv,         &! ppmv over wet air
                                  gas_unit_ppmvdry,      &! ppmv over dry air
                                  gas_id_ozone,          &!
                                  gas_id_co2,            &!
                                  mair,                  &!
                                  mh2o,                  &!
                                  sensor_id_mw,          &!
                                  sensor_id_po,          &!
                                  sensor_id_ir,          &!
                                  sensor_id_hi

  use rttov_types,          only: rttov_options,         &! Structure for RTTOV options
                                  rttov_coef,            &! Structure for RTTOV basic coefficients
                                  rttov_coefs,           &! Structure for the RTTOV coefficients
                                  rttov_coef_pccomp,     &! Structure for RTTOV principal component coefficients
                                  rttov_coef_pccomp1,    &! Structure for RTTOV principal component coefficients
                                  rttov_chanprof,        &!
                                  rttov_coef_scatt_ir,   &! Structure for RTTOV IR scattering coefficients
                                  rttov_optpar_ir,       &! Structure for set of RTTOV IR scattering coefficients
                                  rttov_phasefn_int,     &! interpolation for phase functions
                                  rttov_coef_optpiclb,   &! interpolation factors for frequency
                                  rttov_fast_coef,       &! fast coefficients
                                  rttov_fast_coef_gas,   &! gas fast coefficients
                                  rttov_nlte_coef,       &! internal NLTE coeffs
                                  rttov_profile,         &! structure for atmospheric profiles
                                  rttov_transmission,    &! Transmissions and optical depths (unitless)
                                  rttov_radiance,        &! Radiance and corresponding brightness temperature
                                  rttov_radiance2,       &! upwelling and downwelling radiances
                                  rttov_emissivity,      &!
                                  rttov_coef_pccomp2      ! Structure for RTTOV principal component coefficients
#if defined(_DACE_) && !defined(__ICON__)
  use rttov_types,          only: ipr_deb,               &!
                                  pe_rt => pe
#endif

#if defined(_RTTOV_GOD)
  use rttov_types,          only: rttov_god_par           !

  use rttov_god,            only: rttov_god2o,         &!
                                  rttov_d_god2o,       &!
                                  rttov_god_init        !
#endif

#if defined(_RTTOV_ATLAS)
  use rttov_math_mod,       only: planck
  use mod_rttov_emis_atlas, only: rttov_emis_atlas_data,&! Data type to hold atlas info
                                  atlas_type_ir,       &!
                                  atlas_type_mw,       &!
                                  uwiremis_atlas_id,   &!
                                  camel_atlas_id,      &!
                                  telsem2_atlas_id,    &!
                                  cnrm_mw_atlas_id      !
  use mod_mwatlas_m2,       only: telsem2_atlas_data    ! Data type to hold TELSEM2 atlas info
  use mod_cnrm_mw_atlas,    only: cnrm_mw_atlas_data    ! Data type to hold CNRM atlas info
#endif

#if defined(_RTIFC_USE_MPI_DACE)
  use mo_mpi_dace,          only: p_bcast
#endif
#if defined(_RTIFC_USE_MPI_ICON)
  use mo_mpi,               only: p_bcast
#endif
#if defined(_RTIFC_DISTRIBCOEF) && defined(HAVE_MPI_MOD)
  use mpi     ! prefer MPI module over mpif.h
#endif

#if defined(_RTTOV_USE_OPENMP)
  use omp_lib,              only: omp_get_max_threads     !
#endif

  implicit none

  !-------------
  ! Public stuff
  !-------------

  private

  ! subroutines
  public :: rtifc_version          ! Version string
  public :: rtifc_set_opts_sub     ! Set RTTOV options
  public :: rtifc_get_opts_sub     ! Set RTTOV options
  public :: rtifc_init             ! Initialise RTTOV modules, read coeffs
  public :: rtifc_coef_index       ! Returns index of coeffs for given satid/instr
  public :: rtifc_cleanup          ! frees memory allocated by rtifc_init
  public :: rtifc_fill_input       ! fills the profile-dependent part for RTTOV
  public :: rtifc_direct           ! calls RTTOV direct routine
  public :: rtifc_k                ! calls RTTOV K routine
  public :: rtifc_l2c_god          ! god-corrected l2c
  public :: rtifc_print_profiles   ! print rttov profile structure
  public :: rtifc_coef_prop        ! get properties of coefs
#if defined(_RTTOV_ATLAS)
  public :: rtifc_init_atlas
  public :: rtifc_init_brdf_atlas
  public :: rtifc_emis_atlas
  public :: rtifc_emis_retrieve
  public :: rtifc_emis_sea
  public :: rtifc_brdf_atlas
  public :: rtifc_tskin_retrieve
#endif

  ! RTTOV options
  public :: rttov_options

  ! parameters/variables
  public :: rtifc_vers


  ! RTTOV constants
  public :: gas_unit_specconc
  public :: gas_unit_ppmv
  public :: gas_unit_ppmvdry
  public :: gas_id_ozone
  public :: gas_id_co2

  !-----------
  ! Interfaces
  !-----------
#include "rttov_direct.interface"
#include "rttov_k.interface"
#if defined(_RTTOV_USE_OPENMP)
#include "rttov_parallel_direct.interface"
#include "rttov_parallel_k.interface"
#endif
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_copy_prof.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_coeffname.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_init_transmission.interface"
#include "rttov_read_coefs.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_coefs.interface"
#include "rttov_convert_profile_units.interface"
#include "rttov_apply_reg_limits.interface"
#include "rttov_hdf_save.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_user_options_checkinput.interface"
#if defined(_RTTOV_ATLAS)
#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"
#endif
#include "rttov_nullify_coefs.interface"

#if defined(_RTIFC_DISTRIBCOEF)
#if !defined(HAVE_MPI_MOD)
! include "mpif.h"   ! already imported via "use mo_rtifc_base"
#endif
#endif

  interface rtifc_fill_input
    module procedure rtifc_fill_input_var
#if defined(_DACE_)
    module procedure rtifc_fill_input_rad
#endif
  end interface

  ! MPI routines
#if defined(_RTIFC_DISTRIBCOEF)
  interface p_bcast
    module procedure p_bcast_rttov_coefs
    module procedure p_bcast_rttov_coef
    module procedure p_bcast_rttov_coef_scatt_ir
    module procedure p_bcast_rttov_optpar_ir
    module procedure p_bcast_rttov_coef_pccomp
    module procedure p_bcast_rttov_coef_pccomp1
    module procedure p_bcast_rttov_coef_optpiclb
    module procedure p_bcast_rttov_phasefn_int
    module procedure p_bcast_rttov_fast_coef
    module procedure p_bcast_rttov_fast_coef_gas
    module procedure p_bcast_rttov_nlte_coef
    module procedure p_bcast_rttov_coef_pccomp2
#if defined(_RTTOV_ATLAS)
    module procedure p_bcast_rttov_atlas
    module procedure p_bcast_rttov_telsem
    module procedure p_bcast_rttov_cnrm
#endif
#if defined(_RTTOV_GOD)
    module procedure p_bcast_god_par
#endif
  end interface

  interface p_bcast_rttov_container
    module procedure p_bcast_rttov_cnt_coef
    module procedure p_bcast_rttov_cnt_coef_scatt_ir
    module procedure p_bcast_rttov_cnt_optpar_ir
    module procedure p_bcast_rttov_cnt_coef_pccomp
    module procedure p_bcast_rttov_cnt_coef_pccomp1
    module procedure p_bcast_rttov_cnt_coef_optpiclb
    module procedure p_bcast_rttov_cnt_phasefn_int
    module procedure p_bcast_rttov_cnt_fast_coef
    module procedure p_bcast_rttov_cnt_fast_coef_gas
    module procedure p_bcast_rttov_cnt_nlte_coef
    module procedure p_bcast_rttov_cnt_coef_pccomp2
#if defined(_RTTOV_ATLAS)
    module procedure p_bcast_rttov_cnt_emis_atlas
    module procedure p_bcast_rttov_cnt_telsem
    module procedure p_bcast_rttov_cnt_cnrm
#endif
  end interface
#endif

#if defined(_RTIFC_USE_MPI_DACE) && defined (__GFORTRAN__) && (__GNUC__ >= 10)
  !---------------------------------------------
  ! include interface for external bcast routine
  !---------------------------------------------
  interface p_bcast_derivedtype
     subroutine p_bcast_derivedtype (buffer, count, source, comm)
       type(*) ,intent(inout)     :: buffer     ! variable to bcast
       integer ,intent(in)        :: count      ! len(byte) of variable
       integer ,intent(in)        :: source     ! source processor index
       integer ,intent(in)        :: comm       ! communicator
     end subroutine p_bcast_derivedtype
  end interface p_bcast_derivedtype
#endif
  !----------------------------
  ! Module parameters/variables
  !----------------------------

  ! IFC version
  integer, parameter :: rtifc_vers = 12

  ! RTTOV (de)allocation switches
  integer(jpim), parameter :: rttov_alloc   = 1   ! allocation switch
  integer(jpim), parameter :: rttov_dealloc = 0   ! deallocation switch

  integer,       parameter :: m_atl = 2           ! max. number of atlases for each coef

  ! RTTOV options
  ! Type for organization of options/coeffs
  type t_rtopts
    character(len=100)        :: name    = '' ! Name/description (useful for debugging/err. messages)
    integer                   :: satid   = -1 ! WMO satellite ID
    integer                   :: instr   = -1 ! RTTOV instrument ID
    integer                   :: grid    = -1 ! RTTOV grid ID
    integer                   :: icoeff  = 0  ! index of corresponding rttov_coefs (in "coefs" array)
    type(rttov_options)       :: opts
    integer                   :: iatl(m_atl) = -1 ! index of corresponding atlases in "atlas" array
    integer                   :: natl        =  0 ! Number of associated atlases
   end type t_rtopts
  type(t_rtopts),            pointer            :: rt_opts(:) => null()
  integer                                       :: n_opts = 0
  integer,                   parameter          :: n_alloc_opts = 25
  ! Templates for t_rtopts
  ! (if new templated are added, a new set_opts_template call must be added)
  integer,                   parameter          :: n_tmpl     = 2
  type(t_rtopts),            target             :: rt_opts_tmpl(n_tmpl)
  integer                                       :: idef0      = -1

  ! general dimension variables
  integer(jpim) :: nprof        ! no.of profiles processed in one rttov call
  integer(jpim) :: nmax_chans   ! maximum number of channels of all instrs. wanted
  integer(jpim) :: nlevs        ! number of model levels of external grid
  integer(jpim) :: nlay         ! number of model levels of external grid

  ! RTTOV variables
  type(rttov_coefs),          pointer           :: coefs(:)    => NULL()   ! container for coeffs
  type(rttov_profile),        pointer           :: profiles(:) => NULL()   ! atm. data
  type(rttov_transmission),                save :: transmission            ! transmittances,layer optical depths
  type(rttov_radiance),                    save :: radiance                ! radiances, brightness temperatures
  type(rttov_radiance2),                   save :: radiance2               ! upwelling and downwelling radiances(RTTOV12)
  type(rttov_profile),        pointer           :: profiles_k(:) => NULL() ! atm. data
  type(rttov_transmission),                save :: transmission_k          ! transmittances,layer optical depths
  type(rttov_radiance),                    save :: radiance_k              ! radiances, brightness temperatures
#if defined(_RTTOV_ATLAS)
  type(rttov_emis_atlas_data), pointer          :: atlas(:)      => NULL()
  integer                                       :: n_atlas       =  0
#endif


  ! Check of regularization limits
  integer,            parameter :: nprof_aux = 2
  type(rttov_profile),pointer   :: profiles_aux(:) => NULL()     ! auxiliary profiles


#if !defined(_DACE_) || defined (__ICON__)
  ! we are not sure, that we use a DWD version of RTTOV -> use dummy variables
  integer :: ipr_deb, pe_rt
#endif

  interface construct
    module procedure construct_rttov_options
  end interface construct


contains

#if defined (_FTRACE) && !defined (DISABLE_FTRACE_REGION)
#define FTRACE_BEGIN(text) CALL FTRACE_REGION_BEGIN (text)
#define FTRACE_END(text)   CALL FTRACE_REGION_END   (text)
#else
#define FTRACE_BEGIN(text)
#define FTRACE_END(text)
#endif


  function rtifc_version() result(vers)
    character(len=17) :: vers

    vers = rttov_version()
    write(vers(13:),'("IFC",I2.2)') rtifc_vers
  end function rtifc_version


  subroutine rtifc_set_opts_sub(iopt,               &!
                                new,                &!
                                init,               &!
                                tmpl,               &!
                                iopt_tmpl,          &!
                                name,               &!
                                satid,              &!
                                instr,              &!
                                grid,               &!
                                rttov_opts,         &!
                                addinterp,          &!
                                interp_mode,        &!
                                addrefrac,          &!
                                addclouds,          &!
                                addaerosl,          &!
                                addsolar,           &!
                                addpc,              &!
                                apply_reg_lims,     &!
                                verbose_reg_lims,   &!
                                crop_k_reg_lims,    &!
                                switchrad,          &!
                                conv_overc,         &!
                                fix_hgpl,           &!
                                fastem_version,     &!
                                ir_sea_emis_model,  &!
                                use_q2m,            &!
                                do_lambertian,      &!
                                ozone_data,         &!
                                co2_data,           &!
                                n2o_data,           &!
                                co_data,            &!
                                ch4_data,           &!
                                so2_data,           &!
                                clw_data            &!
                               )
    integer,             intent(inout), optional         :: iopt       ! ID of options, index in rt_opts
    logical,             intent(in),    optional         :: new        ! Create new options, i.e. new entry
                                                                       ! in rt_opts
    logical,             intent(in),    optional         :: init       ! initialization
    character(len=*),    intent(in),    optional         :: tmpl       ! name of template
    integer,             intent(in),    optional         :: iopt_tmpl  ! index of "source" options
    character(len=*),    intent(in),    optional         :: name       ! meta data: name of options
    integer,             intent(in),    optional         :: satid      ! meta data: satellite ID
    integer,             intent(in),    optional         :: instr      ! meta data: instrument
    integer,             intent(in),    optional         :: grid       ! meta data: grid
    type(rttov_options), intent(in),    optional         :: rttov_opts ! full rttov options
    logical,             intent(in),    optional         :: addinterp
    integer,             intent(in),    optional         :: interp_mode
    logical,             intent(in),    optional         :: addrefrac
    logical,             intent(in),    optional         :: addclouds
    logical,             intent(in),    optional         :: addaerosl
    logical,             intent(in),    optional         :: addsolar
    logical,             intent(in),    optional         :: addpc
    logical,             intent(in),    optional         :: apply_reg_lims
    logical,             intent(in),    optional         :: verbose_reg_lims
    logical,             intent(in),    optional         :: crop_k_reg_lims
    logical,             intent(in),    optional         :: switchrad
    logical,             intent(in),    optional         :: conv_overc
    integer,             intent(in),    optional         :: fix_hgpl
    integer,             intent(in),    optional         :: fastem_version
    integer,             intent(in),    optional         :: ir_sea_emis_model
    logical,             intent(in),    optional         :: use_q2m
    logical,             intent(in),    optional         :: do_lambertian
    logical,             intent(in),    optional         :: ozone_data
    logical,             intent(in),    optional         :: co2_data
    logical,             intent(in),    optional         :: n2o_data
    logical,             intent(in),    optional         :: co_data
    logical,             intent(in),    optional         :: ch4_data
    logical,             intent(in),    optional         :: so2_data
    logical,             intent(in),    optional         :: clw_data
    !--------------------------------------------------------------------------
    ! Set options in rttov_options type for later use in rtifc_* routines.
    !   Explanantion of crop_k_reg_lims Option by DWD (RF):
    !   If the humidity input to rttov12 exceeds the regression limits, we get
    !   humi_k=0. above level_hum_dum, which is consistent. The humi_k values
    !   above level_hum_dum are not relevant for the assimilation. But they
    !   are written into the *RTOVP* files. In order to have realistic profiles in the
    !   *RTOVP* files, we prevent RTTOV12 from setting the results to zero, where
    !   the reg_limits were exceeded.
    !--------------------------------------------------------------------------
    character(len=18),   parameter :: proc = 'rtifc_set_opts_sub'
    type(rttov_options), pointer   :: ropts => null()
    logical                        :: new_  =  .false.
    type(t_rtopts),      pointer   :: rt_opts_(:) => null() ! auxiliary for reallocation of rt_opts
    type(t_rtopts),      pointer   :: rto         => null()
    integer                        :: mx_opts
    integer                        :: i, itmpl
    logical                        :: l_mod_tmpl  = .false.



    rto => null()
    l_mod_tmpl = .false.

    ! Initialize templates if necessary
    if (idef0 <= 0) then
      idef0 = 1
      call set_opts_template(rt_opts_tmpl(idef0), 'default0')
      call set_opts_template(rt_opts_tmpl(    2), 'default' )
    end if

    ! Find template
    itmpl = -1
    if (present(tmpl)) then
      do i = 1, n_tmpl
        if (trim(tmpl) == trim(rt_opts_tmpl(i)%name)) then
          itmpl = i
          exit
        end if
      end do
      if (itmpl <= 0) call finish(proc, 'Did not find rt_opts-template "'//trim(tmpl)//'"')
    end if

    if (present(new)) then
      new_ = new
    else
      new_ = .false.
    end if
    if (new_) then
      if (.not.present(iopt)) call finish(proc, 'arguments: new without iopt not supported.')

      if (associated(rt_opts)) then
        mx_opts = size(rt_opts)
      else
        mx_opts = 0
      end if
      if (n_opts + 1 > mx_opts) then
        ! Reallocate rt_opts
        allocate(rt_opts_(mx_opts + n_alloc_opts))
        if (associated(rt_opts)) then
          rt_opts_(1:n_opts) = rt_opts(1:n_opts)
          deallocate(rt_opts)
        end if
        rt_opts => rt_opts_
      end if
      ! Create and initialize new entry in rt_opts
      n_opts = n_opts + 1
      iopt = n_opts
      rto => rt_opts(iopt)
      if (present(tmpl)) then
        rto = rt_opts_tmpl(itmpl)
      elseif (present(iopt_tmpl)) then
        if (iopt_tmpl < 0 .or. iopt_tmpl >= n_opts) call finish(proc, 'invalid iopt_tmpl')
        rto = rt_opts(iopt_tmpl)
      else
        rto = rt_opts_tmpl(idef0)
      end if
    elseif (present(iopt)) then
      if (iopt < 1 .or. iopt > n_opts) call finish(proc, 'invalid iopt')
      rto => rt_opts(iopt)
    elseif (present(tmpl)) then
      ! Modify template
      rto => rt_opts_tmpl(itmpl)
      l_mod_tmpl = .true.
    else
      call finish(proc, 'Either iopt or tmpl must be present')
    end if

    if (.not.l_mod_tmpl) then
      if (present(name )) rto%name  = trim(name)
      if (present(satid)) rto%satid = satid
      if (present(instr)) rto%instr = instr
      if (present(grid )) rto%grid  = grid
    end if

    ropts => rto%opts

    if (present(init)) then
      if (init) then
        call set_opts_template(rto)
      end if
    end if

    if (present(rttov_opts)) ropts = rttov_opts

    if (present(addinterp        )) ropts%interpolation%addinterp  = addinterp
    if (present(interp_mode      )) ropts%interpolation%interp_mode= interp_mode
    if (present(addrefrac        )) ropts%rt_all%addrefrac         = addrefrac
    if (present(addclouds        )) ropts%rt_ir%addclouds          = addclouds
    if (present(addaerosl        )) ropts%rt_ir%addaerosl          = addaerosl
    if (present(addsolar         )) ropts%rt_ir%addsolar           = addsolar
    if (present(addpc            )) ropts%rt_ir%pc%addpc           = addpc
    if (present(apply_reg_lims   )) ropts%config%apply_reg_limits  = apply_reg_lims
    if (present(verbose_reg_lims )) ropts%config%verbose           = verbose_reg_lims
    if (present(switchrad        )) ropts%rt_all%switchrad         = switchrad
    if (present(fastem_version   )) ropts%rt_mw%fastem_version     = fastem_version
    if (present(ir_sea_emis_model)) ropts%rt_ir%ir_sea_emis_model  = ir_sea_emis_model
    if (present(use_q2m          )) ropts%rt_all%use_q2m           = use_q2m
    if (present(do_lambertian    )) ropts%rt_all%do_lambertian     = do_lambertian
#if defined(_DACE_) && !defined(__ICON__)
    if (present(crop_k_reg_lims  )) ropts%config%crop_k_reg_limits = crop_k_reg_lims
    if (present(conv_overc       )) ropts%config%conv_overc        = conv_overc
    if (present(fix_hgpl         )) ropts%config%fix_hgpl          = fix_hgpl
#endif
    if (present(ozone_data       )) ropts%rt_ir%ozone_data         = ozone_data
    if (present(co2_data         )) ropts%rt_ir%co2_data           = co2_data
    if (present(n2o_data         )) ropts%rt_ir%n2o_data           = n2o_data
    if (present(co_data          )) ropts%rt_ir%co_data            = co_data
    if (present(ch4_data         )) ropts%rt_ir%ch4_data           = ch4_data
    if (present(so2_data         )) ropts%rt_ir%so2_data           = so2_data
    if (present(clw_data         )) ropts%rt_mw%clw_data           = clw_data

    if (.not.ropts%config%apply_reg_limits) then
      ropts%config%do_checkinput = .false.
    end if

  end subroutine rtifc_set_opts_sub


  subroutine rtifc_get_opts_sub(iopt,               &!
                                tmpl,               &!
                                name,               &!
                                satid,              &!
                                instr,              &!
                                grid,               &!
                                rttov_opts,         &!
                                addinterp,          &!
                                interp_mode,        &!
                                addrefrac,          &!
                                addclouds,          &!
                                addaerosl,          &!
                                addsolar,           &!
                                addpc,              &!
                                apply_reg_lims,     &!
                                verbose_reg_lims,   &!
                                crop_k_reg_lims,    &!
                                switchrad,          &!
                                conv_overc,         &!
                                fix_hgpl,           &!
                                fastem_version,     &!
                                ir_sea_emis_model,  &!
                                use_q2m,            &!
                                do_lambertian,      &!
                                ozone_data,         &!
                                co2_data,           &!
                                n2o_data,           &!
                                co_data,            &!
                                ch4_data,           &!
                                so2_data,           &!
                                clw_data            &!
                                )
    integer,             intent(in),     optional :: iopt       ! ID of options, index in rt_opts
    character(len=*),    intent(in),     optional :: tmpl       ! name of template
    character(len=*),    intent(out),    optional :: name       ! meta data: name of options
    integer,             intent(out),    optional :: satid      ! meta data: satellite ID
    integer,             intent(out),    optional :: instr      ! meta data: instrument
    integer,             intent(out),    optional :: grid       ! meta data: grid
    type(rttov_options), intent(out),    optional :: rttov_opts ! full rttov options
    logical,             intent(out),    optional :: addinterp
    integer,             intent(out),    optional :: interp_mode
    logical,             intent(out),    optional :: addrefrac
    logical,             intent(out),    optional :: addclouds
    logical,             intent(out),    optional :: addaerosl
    logical,             intent(out),    optional :: addsolar
    logical,             intent(out),    optional :: addpc
    logical,             intent(out),    optional :: apply_reg_lims
    logical,             intent(out),    optional :: verbose_reg_lims
    logical,             intent(out),    optional :: crop_k_reg_lims
    logical,             intent(out),    optional :: switchrad
    logical,             intent(out),    optional :: conv_overc
    integer,             intent(out),    optional :: fix_hgpl
    integer,             intent(out),    optional :: fastem_version
    integer,             intent(out),    optional :: ir_sea_emis_model
    logical,             intent(out),    optional :: use_q2m
    logical,             intent(out),    optional :: do_lambertian
    logical,             intent(out),    optional :: ozone_data
    logical,             intent(out),    optional :: co2_data
    logical,             intent(out),    optional :: n2o_data
    logical,             intent(out),    optional :: co_data
    logical,             intent(out),    optional :: ch4_data
    logical,             intent(out),    optional :: so2_data
    logical,             intent(out),    optional :: clw_data
    !--------------------------------------------------------------------------
    ! Get options from the rt_opts or rt_opts_tmpl arrays
    !--------------------------------------------------------------------------
    character(len=18),   parameter :: proc = 'rtifc_get_opts_sub'
    type(rttov_options), pointer   :: ropts => null()
    type(t_rtopts),      pointer   :: rto         => null()
    integer                        :: i, io

    if (present(iopt) .and. present(tmpl)) then
      call finish(proc, 'Either iopt or tmpl option allowed, but not both.')
    elseif (present(tmpl)) then
      io = -1
      do i = 1, n_tmpl
        if (trim(tmpl) == trim(rt_opts_tmpl(i)%name)) then
          io = i
          exit
        end if
      end do
      if (io <= 0) call finish(proc, 'Did not find rt_opts-template "'//trim(tmpl)//'"')
      rto => rt_opts_tmpl(io)
    elseif (present(iopt)) then
      if (iopt > 0 .and. iopt <= n_opts) then
        rto => rt_opts(iopt)
      else
        call finish(proc, 'Options index iopt out of range.')
      end if
    else
      call finish(proc, 'Either iopt or tmpl option must be given.')
    end if

    if (present(name )) name  = trim(rto%name)
    if (present(satid)) satid = rto%satid
    if (present(instr)) instr = rto%instr
    if (present(grid )) grid  = rto%grid

    ropts => rto%opts

    if (present(rttov_opts)) rttov_opts = ropts

    if (present(addinterp        )) addinterp         = ropts%interpolation%addinterp
    if (present(interp_mode      )) interp_mode       = ropts%interpolation%interp_mode
    if (present(addrefrac        )) addrefrac         = ropts%rt_all%addrefrac
    if (present(addclouds        )) addclouds         = ropts%rt_ir%addclouds
    if (present(addaerosl        )) addaerosl         = ropts%rt_ir%addaerosl
    if (present(addsolar         )) addsolar          = ropts%rt_ir%addsolar
    if (present(addpc            )) addpc             = ropts%rt_ir%pc%addpc
    if (present(apply_reg_lims   )) apply_reg_lims    = ropts%config%apply_reg_limits
    if (present(verbose_reg_lims )) verbose_reg_lims  = ropts%config%verbose
    if (present(switchrad        )) switchrad         = ropts%rt_all%switchrad
    if (present(fastem_version   )) fastem_version    = ropts%rt_mw%fastem_version
    if (present(ir_sea_emis_model)) ir_sea_emis_model = ropts%rt_ir%ir_sea_emis_model
    if (present(use_q2m          )) use_q2m           = ropts%rt_all%use_q2m
    if (present(do_lambertian    )) do_lambertian     = ropts%rt_all%do_lambertian
#if defined(_DACE_) && !defined(__ICON__)
    if (present(crop_k_reg_lims  )) crop_k_reg_lims   = ropts%config%crop_k_reg_limits
    if (present(conv_overc       )) conv_overc        = ropts%config%conv_overc
    if (present(fix_hgpl         )) fix_hgpl          = ropts%config%fix_hgpl
#endif
    if (present(ozone_data       )) ozone_data        = ropts%rt_ir%ozone_data
    if (present(co2_data         )) co2_data          = ropts%rt_ir%co2_data
    if (present(n2o_data         )) n2o_data          = ropts%rt_ir%n2o_data
    if (present(co_data          )) co_data           = ropts%rt_ir%co_data
    if (present(ch4_data         )) ch4_data          = ropts%rt_ir%ch4_data
    if (present(so2_data         )) so2_data          = ropts%rt_ir%so2_data
    if (present(clw_data         )) clw_data          = ropts%rt_mw%clw_data

  end subroutine rtifc_get_opts_sub


  subroutine construct_rttov_options(ropts)
    type(rttov_options), intent(out) :: ropts
  end subroutine construct_rttov_options


  subroutine set_opts_template(rto, tmpl)
    type(t_rtopts),   intent(out), target   :: rto
    character(len=*), intent(in),  optional :: tmpl

    character(len=17), parameter :: proc = 'set_opts_template'
    type(rttov_options), pointer :: ropts

    ropts => rto%opts

    call construct(ropts)
    ! Overwrite some NWPSAF defaults with DWD defaults
    ropts%interpolation%addinterp  = .false.
    ropts%rt_all%addrefrac         = .true.
    ropts%rt_ir%addclouds          = .false.
    ropts%rt_ir%addaerosl          = .false.
    ropts%rt_ir%addsolar           = .false.
    ropts%rt_ir%pc%addpc           = .false.
    ropts%config%apply_reg_limits  = .false.
    ropts%config%verbose           = .false.
    ropts%rt_all%switchrad         = .true.
    ropts%rt_mw%fastem_version     = 5
    ropts%rt_ir%ir_sea_emis_model  = 2
    ropts%rt_all%use_q2m           = .false.
    ropts%rt_all%do_lambertian     = .false.
#if defined(_DACE_)
    ropts%config%crop_k_reg_limits = .false.
    ropts%config%conv_overc        = .false.
    ropts%config%fix_hgpl          = 2
#endif
    ropts%rt_ir%ozone_data         = .false.
    ropts%rt_ir%co2_data           = .false.
    ropts%rt_ir%n2o_data           = .false.
    ropts%rt_ir%co_data            = .false.
    ropts%rt_ir%ch4_data           = .false.
    ropts%rt_ir%so2_data           = .false.
    ropts%rt_mw%clw_data           = .false.

    if (present(tmpl)) then
      select case(trim(tmpl))
      case('default0')
        ! Old default for old (deprecated) namelists
      case('default')
        ! Do nothing
      case default
        ! TODO: do not crash here, but set a status, that is given back to the user?
        call finish(proc, 'Options template "'//trim(tmpl)//'" not defined.')
      end select
      rto%name = trim(tmpl)
    end if
  end subroutine set_opts_template


  subroutine match_opts_coefs(opts,chans,coefs,match)
    type(rttov_options), intent(in)  :: opts
    integer,             intent(in)  :: chans(:)
    type(rttov_coefs),   intent(in)  :: coefs
    logical,             intent(out) :: match
    ! tells us whether coefs match with opts, i.e. whether the coefs were
    ! read with different options.
    ! rttov_read_coeffs.f90 was used as a template.
    integer :: i

    match = .true.

    !RF: i think, that it is not possible to reconstruct the channel indices that were
    !    used to load the coefficients. So we have to load the coeffs again, if not all
    !    channels were loaded
    match = match .and. (coefs%coef%fmv_ori_nchn == coefs%coef%fmv_chn)
    do i = 1, size(chans)
      if (all(chans(i) /= coefs%coef%ff_ori_chn(1:coefs%coef%fmv_ori_nchn))) then
        match = .false.
        return
      end if
    end do

    if (opts%rt_ir%pc%addpc) then
      ! pccomp file read?
      match = match .and. (coefs%coef_pccomp%fmv_pc_comp_pc > 0)
      if (opts%rt_ir%addclouds) then
        match = match .and. (coefs%coef_pccomp%fmv_pc_cld > 0)
      end if
      if (opts%rt_ir%pc%addradrec) then
        match = match .and. (coefs%coef_pccomp%fmv_pc_nchn > 0)
      endif
    endif
  end subroutine match_opts_coefs


  subroutine merge_opts(optsm,opts)
    type(rttov_options), intent(inout) :: optsm
    type(rttov_options), intent(in)    :: opts
    ! Purpose: to obtain options, that are adequate for loading the coefficients or
    ! allocating the profiles in a way that they can be used with multiple different
    ! rttov_options

    optsm%rt_ir%ozone_data        = optsm%rt_ir%ozone_data .OR. opts%rt_ir%ozone_data
    optsm%rt_ir%co2_data          = optsm%rt_ir%co2_data   .OR. opts%rt_ir%co2_data
    optsm%rt_ir%n2o_data          = optsm%rt_ir%n2o_data   .OR. opts%rt_ir%n2o_data
    optsm%rt_ir%co_data           = optsm%rt_ir%co_data    .OR. opts%rt_ir%co_data
    optsm%rt_ir%ch4_data          = optsm%rt_ir%ch4_data   .OR. opts%rt_ir%ch4_data
    optsm%rt_ir%so2_data          = optsm%rt_ir%so2_data   .OR. opts%rt_ir%so2_data
    optsm%rt_mw %clw_data          = optsm%rt_mw %clw_data   .OR. opts%rt_mw %clw_data

    optsm%rt_ir%addaerosl          = optsm%rt_ir%addaerosl          .OR.  opts%rt_ir%addaerosl
    optsm%rt_ir%user_aer_opt_param = optsm%rt_ir%user_aer_opt_param .AND. opts%rt_ir%user_aer_opt_param

    optsm%rt_ir%addclouds          = optsm%rt_ir%addclouds          .OR.  opts%rt_ir%addclouds
    optsm%rt_ir%user_cld_opt_param = optsm%rt_ir%user_cld_opt_param .AND. opts%rt_ir%user_cld_opt_param

    optsm%rt_ir%pc%addpc     = optsm%rt_ir%pc%addpc     .OR. opts%rt_ir%pc%addpc
    optsm%rt_ir%pc%addradrec = optsm%rt_ir%pc%addradrec .OR. opts%rt_ir%pc%addradrec

  end subroutine merge_opts


  subroutine rtifc_init(instruments,channels,nchans_inst,ch_id,iopts,my_proc_id,n_proc,  &
                        io_proc_id,mpi_comm_type,status,path_coefs)
    integer,             intent(in)          :: instruments(:,:)  ! info on processed instruments:
                                                                  ! (1,1:nInstr): rttov platform IDs
                                                                  ! (2,1:nInstr): rttov satellite IDs
                                                                  ! (3,1:nInstr): rttov instrument IDs
    integer,             intent(in)          :: channels(:,:)     ! channels to extract:
                                                                  ! (1:nchansPerInst,1:nInstr);
                                                                  ! numbers in channels are the indices
                                                                  ! of the channels in the coeff. file
    integer,             intent(in)          :: nchans_inst(:)    ! number of channels for each instrument
                                                                  ! as defined in channels
    integer,             intent(out)         :: ch_id(:,:)        ! channel indices (for RTTOV calls)
    integer,             intent(in)          :: iopts(:)          ! rttov_options for each instrument
    integer,             intent(in)          :: my_proc_id        ! ID of local processors
    integer,             intent(in)          :: n_proc            ! no. of processors in mpi communication domain
    integer,             intent(in)          :: io_proc_id        ! ID of IO processor
    integer,             intent(in)          :: mpi_comm_type     ! mpi communicator type
    integer,             intent(out)         :: status            ! exit status
    character(*),        intent(in), optional:: path_coefs        ! path to coefficient files
    !--------------------------------------------------------------------------
    ! initialization routine which initializes the rttov modules
    ! and reads the instrument specific coefficients which are needed (wanted).
    !--------------------------------------------------------------------------
    character(len=10), parameter :: proc = 'rtifc_init'
    type(rttov_coefs), pointer   :: coefs_tmp(:) => null()
    type(t_rtopts),    pointer   :: rto          => null()
    type(rttov_options)          :: opts_
    integer                      :: chans(mx_chan), ch_id_(mx_chan)
    character(len=120)           :: msg = ''
    integer                      :: errstat(size(instruments,2)), stat
    integer                      :: ninstr, ncoefs, ncoefs_, nchans
    integer                      :: instr, stat_
    integer                      :: i, j, k, ic
    integer                      :: id_platform, id_sat, id_inst
    logical                      :: l_distrib, match

FTRACE_BEGIN('rtifc_init')

    status = NO_ERROR

    l_distrib = .false.
#if defined(_RTIFC_DISTRIBCOEF)
    l_distrib = (n_proc > 1) .and. read1pe
#endif

    ! initialize some module variables
    nprof      = 0
    nlevs      = 0
    nmax_chans = maxval(nchans_inst(:))

    if (verbosity >= production .and. my_proc_id == io_proc_id) &
         print*,'Initialize '//rtifc_version()

    pe_ifc = my_proc_id
    pe_rt  = pe_ifc

    ! Determine coefs to be loaded or link to existing coefs
    ninstr     = size(instruments, 2)
    if (size(iopts) < ninstr) call finish(proc,'iopts smaller than number of instruments.')
    ncoefs = 0
    if (associated(coefs)) ncoefs = size(coefs)
    ncoefs_ = ncoefs
    do i = 1, ninstr
      rto => rt_opts(iopts(i))
      if (rto%icoeff > 0) call finish(proc, 'rt_opts are associated with coeffs already')
      id_platform = instruments(1,i)
      id_sat      = instruments(2,i)
      id_inst     = instruments(3,i)
      do ic = 1, ncoefs
        if ( id_platform == coefs(ic)%coef%id_platform .and. &
             id_sat      == coefs(ic)%coef%id_sat      .and. &
             id_inst     == coefs(ic)%coef%id_inst     ) then
          ! Can we use these coeffs for the new options?
          call match_opts_coefs(rto%opts, channels(1:nchans_inst(i),i), coefs(ic), match)
          if (.not.match) cycle
          call rttov_user_options_checkinput(stat, rto%opts, coefs(ic))
          if (stat /= 0) call finish(proc, 'bad options or existing coeffs are incompatible')
          rto%icoeff = ic
          ch_id(1:nchans_inst(i),i) = channels(1:nchans_inst(i),i)
        end if
      end do
      if (rto%icoeff <= 0) then
        ! check whether this instrument occured already
        do j = i-1, 1, -1
          if ( id_platform == instruments(1,j) .and. &
               id_sat      == instruments(2,j) .and. &
               id_inst     == instruments(3,j) ) then
            rto%icoeff = rt_opts(iopts(j))%icoeff
            exit
          end if
        end do
        if (rto%icoeff <= 0) then
          ! we have to initialize new coefficients
          ncoefs_ = ncoefs_ + 1
          rto%icoeff = ncoefs_
        end if
      end if
    end do

    if (verbosity >= production .and. my_proc_id == io_proc_id) then
      write(msg,'("Read ",I2," RTTOV coefficient files")') ncoefs_ - ncoefs
      if (l_distrib) msg = trim(msg)//' on one pe and distribute them'
      print*,trim(msg)
    end if

    ! allocate/reallocate coefficients
    if (ncoefs_ > ncoefs) then
      allocate(coefs_tmp(ncoefs_), stat=stat_)
      if(stat_ /= 0) then
        status = ERR_ALLOC
        return
      end if
      if (ncoefs > 0) coefs_tmp(1:ncoefs) = coefs(1:ncoefs)
      if (associated(coefs)) deallocate(coefs)
      coefs => coefs_tmp
    end if


    errstat = 0
    do instr = 1, ninstr
      rto => rt_opts(iopts(instr))
      ic = rto%icoeff
      if (coefs(ic)%initialised) cycle
      FTRACE_BEGIN('read_coeffs')

      ! Merge options/channels with following options/channels for the same instrument
      opts_ = rto%opts
      nchans = nchans_inst(instr)
      chans(1:nchans_inst(instr)) = channels(1:nchans_inst(instr),instr)
      do i = instr+1, ninstr
        if (rt_opts(iopts(i))%icoeff == ic) then
          call merge_opts(opts_,rt_opts(iopts(i))%opts)
          do j = 1, nchans_inst(i)
            if (all(chans(1:nchans) /= channels(j,i))) then
              nchans = nchans + 1
              if (nchans > size(chans)) call finish(proc,'chans array too small. &
                   &Recompile with increased mx_chan in mo_rtifc ')
              chans(nchans) = channels(j,i)
            end if
          end do
        end if
      end do
      ! sort channel
      do i = nchans-1, 1, -1
        do j = 1, i
          if (chans(j) > chans(j+1)) then
            k          = chans(j)
            chans(j)   = chans(j+1)
            chans(j+1) = k
          end if
        end do
      end do
      ! Determine channel indices
      do i = 1, nchans
        if (chans(i) <= 0 .or. chans(i) > mx_chan) call finish(proc, 'Either invalid &
             &channel number of mx_chan too small.')
        ch_id_(chans(i)) = i
      end do
      ch_id(1:nchans_inst(instr),instr) = ch_id_(channels(1:nchans_inst(instr),instr))
      do i = instr+1, ninstr
        if (rt_opts(iopts(i))%icoeff == ic) then
          ch_id(1:nchans_inst(i),i) = ch_id_(channels(1:nchans_inst(i),i))
        end if
      end do

      ! Read or initialize coeffs
      if (io_proc_id == my_proc_id .or. .not.l_distrib) then
        call rttov_read_coefs(                  &
             errstat(instr),                    &! -> return code
             coefs(ic),                         &! -> coefficients
             opts_,                             &! <- options
             channels   = chans(1:nchans),      &! <- channels to load
             instrument = instruments(:,instr), &!
             path       = path_coefs)
      else
        call rttov_nullify_coefs(coefs(ic))
      endif
      FTRACE_END('read_coeffs')
#if defined(_RTIFC_DISTRIBCOEF)
      ! check errstat
      if (l_distrib) call p_bcast(errstat(instr:instr), io_proc_id, mpi_comm_type)
      if (errstat(instr) /= errorstatus_success) cycle
      ! distribute loaded coefficients to all PEs
      if (l_distrib) then
        FTRACE_BEGIN('distrib_coeffs')
        call p_bcast(coefs(ic),io_proc_id,mpi_comm_type)
        FTRACE_END('distrib_coeffs')
      end if
#else
      ! check errstat
      if (errstat(instr) /= errorstatus_success) cycle
#endif
      if (.not.coefs(ic)%initialised) call finish('proc', 'coefs not initialised after reading')

   enddo

    if(any(errstat(1:ninstr) /= errorstatus_success)) then
       print *, 'rtifc_init: stat =', errstat
       status = ERR_RTTOV_SETUP
       return
    endif

    if(associated(profiles)) then
      call dealloc_rttov_arrays(status, profs=profiles, rads=radiance, rads2=radiance2, transm=transmission)
      if (status /= NO_ERROR) return
    end if
    if(associated(profiles_k)) then
      call dealloc_rttov_arrays(status, profs=profiles_k, rads=radiance_k, transm=transmission_k)
      if (status /= NO_ERROR) return
    end if

    if (ncoefs_ > ncoefs) then
      ! Check loaded coeffs
      do i = 1, ninstr
        rto => rt_opts(iopts(i))
        if (rto%icoeff > ncoefs) then
          call rttov_user_options_checkinput(stat, rto%opts, coefs(rto%icoeff))
          if (stat /= 0) call finish(proc, 'bad options or existing coeffs are incompatible')
        end if
      end do
      nlevs = coefs(1)%coef%nlevels
    end if

#if defined(_RTTOV_GOD)
    FTRACE_BEGIN('read_god_par')
    if (god_par_file /= '') then
#if defined(_RTIFC_DISTRIBCOEF)
      if (io_proc_id == my_proc_id .or. .not.l_distrib) call read_god_par(status)
      if (l_distrib) call p_bcast(status, io_proc_id, mpi_comm_type)
      if (status /= NO_ERROR) return
      if (l_distrib) call bcast_god_par(my_proc_id,io_proc_id,mpi_comm_type)
#else
      call read_god_par(status)
      if (status /= NO_ERROR) return
#endif
    end if
    FTRACE_END('read_god_par')
#endif

    if(my_proc_id == io_proc_id .and. verbosity >= production) then
       print *, 'Successfully initialized RTTOV.'
    endif

    FTRACE_END('rtifc_init')

  end subroutine rtifc_init



#if defined(_RTTOV_GOD)
  subroutine read_god_par(stat)
    integer, intent(inout) :: stat
    !-------------------
    ! Reads god_par_file
    !-------------------
    type(rttov_coef),    pointer :: coef => null()
    integer                      :: iu, istat, i, icoef, i_entry, ichan, n_used
    logical                      :: ldum
    ! File content
    character(len=10000)         :: line
    character(len=100)           :: tag
    integer                      :: version
    integer                      :: nlevels   ! Number of rttov levels
    ! god parameters
    integer                      :: instr, chan, layer, mode
    real(kind=jprb)              :: opdep, par1
    type(rttov_god_par), pointer :: gp

#define ERR stat = ERR_GOD_FILE ; close(iu) ; return
    stat = 0
    if (god_par_file /= '') then
      inquire(file=trim(god_par_file), exist=ldum)
      if (.not.ldum) then
        stat = ERR_GOD_FILE ; return
      end if
      iu = get_unit_number()
      open(iu, file=trim(god_par_file))
      version = 1
      nlevels = 54
      i_entry = 0
      n_used  = 0
      read_god_loop: do
        read(iu, '(A)', iostat=istat) line
        !print*,'line',trim(line)
        if (istat /= 0) then
          if (istat == iostat_end) then
            if (verbosity >= production) then
              write(stdout,*)
              write(stdout,*) 'Read god-parameters from "'//trim(god_par_file)//'":'
              write(stdout,*) '  #Entries          : ',i_entry
              write(stdout,*) '  #Useful entries   : ',n_used
            end if
            close(iu)
            exit
          else
            ERR
          end if
        end if
        line = adjustl(line)

        if (line == ''                            ) cycle read_god_loop ! Empty line
        if (line(1:1) == '#' .or. line(1:1) == '!') cycle read_god_loop ! Comment line

        if (i_entry == 0) then
          read(line,*,iostat=istat) tag, i
          if (istat == 0) then
            if (tag == 'version' .or. tag == 'Version') then
              version = i
              cycle read_god_loop
            end if
            if (tag == 'nlevels') then
              nlevels = i
              cycle read_god_loop
            end if
          end if
        end if

        istat = -1
        select case(version)
        case(1)
          read(line,*,iostat=istat) instr, chan, layer, mode, opdep, par1
        case default
          write(0,*) 'Invalid version ',version,' of god_par_file "'//trim(god_par_file)//'"'
          ERR
        end select
        if (istat == 0) then
          i_entry = i_entry + 1
          ldum = .false.
          do icoef = 1, size(coefs)
            coef => coefs(icoef)%coef
            if (coef%id_inst == instr) then
              if (layer < 1 .or. layer > coef%nlayers) then
                write(0,*) 'Invalid layer in god_par_file entry ',i_entry,': '//trim(line)
                ERR
              end if
              if (coef%nlevels /= nlevels) then
                write(0,*) 'level number of god_par_file does not match with RTTOV coeffs:'
                write(0,*) 'nlevels (god_par): ', nlevels
                write(0,*) 'nlevels (coeffs.): ', coef%nlevels
                ERR
              end if
              if (.not.associated(coef%god)) then
                allocate(coef%god(coef%nlayers, coef%fmv_chn))
                coef%god%ntr = 0
              end if
              do ichan = 1, coef%fmv_chn
                if (coef%ff_ori_chn(ichan) == chan) then
                  gp => coef%god(layer,ichan)
                  gp%ntr = gp%ntr + 1
                  if (gp%ntr > size(gp%tr)) then
                    write(0,*) 'Too many entries for instr=',instr,' chan=',chan,' layer=',layer,&
                         &' in god_par_file.'
                    ERR
                  end if
                  gp%tr(gp%ntr)%mode  = mode
                  gp%tr(gp%ntr)%opdep = opdep
                  gp%tr(gp%ntr)%par1  = par1
                  ldum = .true.
                end if
              end do
            end if
          end do
          if (ldum) then
            n_used = n_used + 1
          else
            if (verbosity >= verbose) write(stdout,*) 'not used entry: '//trim(line)
          end if
        else
          write(0,*) 'Invalid entry "',trim(line)//'" in god_par_file.'
          ERR
        end if
      end do read_god_loop
      close(iu)

      do icoef = 1, size(coefs)
        coef => coefs(icoef)%coef
        if (associated(coef%god)) then
          do ichan = 1, coef%fmv_chn
            do layer = 1, coef%nlayers
              gp => coef%god(layer, ichan)
              if (gp%ntr > 0) then
                call rttov_god_init(p=gp)
                if (wr_god) call write_god(gp, coef%id_inst, coef%ff_ori_chn(ichan), layer)
              end if
            end do
          end do
        end if
      end do

    end if
#undef ERR
  end subroutine read_god_par

#if defined(_RTIFC_DISTRIBCOEF)

  subroutine bcast_god_par(my_proc_id,io_proc_id,mpi_comm_type)
    integer,          intent(in) :: my_proc_id             ! ID of processor
    integer,          intent(in) :: io_proc_id             ! ID of IO processor
    integer,          intent(in) :: mpi_comm_type          ! mpi communicator type
    type(rttov_coef), pointer    :: coef => null()
    integer(kind=jpim)           :: icoef, ichan, layer
    logical                      :: l_io, l_assoc
    l_io = (io_proc_id == my_proc_id)
    do icoef = 1, size(coefs)
      coef => coefs(icoef)%coef
      l_assoc = associated(coef%god)
      call p_bcast(l_assoc,io_proc_id,mpi_comm_type)
      if (l_assoc) then
        if (.not.l_io) allocate(coef%god(coef%nlayers, coef%fmv_chn))
        do ichan = 1, coef%fmv_chn
          do layer = 1, coef%nlayers
            call p_bcast(coef%god(layer,ichan),io_proc_id,mpi_comm_type)
          end do
        end do
      end if
    end do
  end subroutine bcast_god_par


  subroutine p_bcast_god_par(buffer,source,comm)
    type(rttov_god_par), intent(inout)   :: buffer
    integer,             intent(in)      :: source
    integer,    optional,intent(in)      :: comm
    !------------------------------------------------------------------
    ! Broadcast an rttov_coef container across all available processors
    !------------------------------------------------------------------
    character(len=15), parameter :: proc = 'p_bcast_god_par'
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish(proc, 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_god_par

#endif /* _RTIFC_DISTRIBCOEF */

  subroutine write_god(god, instr, chan, layer)
    type(rttov_god_par), intent(in) :: god
    integer,             intent(in) :: instr
    integer,             intent(in) :: chan
    integer,             intent(in) :: layer

    integer, parameter  :: n = 1000
    integer             :: i, iu
    real(kind=wp)       :: x, dx, od
    character(len=100)  :: fname

    if (god%ntr > 0) then
      od = god%tr(1)%opdep
      dx = 6*od / (n-1)
      iu = get_unit_number()
      write(fname, '("god_instr",I3.3,"_chan",I5.5,"_layer",I3.3,".dat")') instr, chan, layer
      if (out_path /= '') fname = trim(out_path)//'/'//trim(fname)
      open(iu, file=trim(fname))
      print*,trim(fname),god%ntr
      do i = 1, n
        x = -2*od + (i-1) * dx
        write(iu, *) x, rttov_god2o(p=god, god=x)
      end do
      write(iu, *) '&'
      do i = 1, n
        x = -2*od + (i-1) * dx
        write(iu, *) x, rttov_d_god2o(p=god, god=x)
      end do

      write(iu, *) '&'
      do i = 1, n
        x = -2*od + (i-1) * dx
        write(iu, *) x, rttov_d_god2o(p=god, o=rttov_god2o(x, god))
      end do

      close(iu)
    end if
  end subroutine write_god

  function get_unit_number() result(iu)
    integer :: iu
    logical :: opened
    do iu = 10, 100000
      inquire(iu, opened=opened)
      if (.not.opened) return
    end do
  end function get_unit_number

#endif /* _RTTOV_GOD */


  subroutine rtifc_l2c_god(iopt, chans, valid, l2c, l2c_corr, transm, opdep, p_l2c, ideb)
    integer,             intent(in)           :: iopt
    integer,             intent(in)           :: chans(:)
    logical,             intent(in)           :: valid(:,:)
    real,                intent(in)           :: l2c(:,:)
    real,                intent(out)          :: l2c_corr(:,:)
    real(kind=jprb),     intent(in), optional :: transm(:,:,:)
    real(kind=jprb),     intent(in), optional :: opdep(:,:,:)
    real(kind=jprb),     intent(in), optional :: p_l2c(:,:)
    integer,             intent(in), optional :: ideb(:)

    character(len=13),   parameter   :: proc   = 'rtifc_l2c_god'
    type(t_rtopts),      pointer     :: rto    => null()
    type(rttov_options), pointer     :: ropts  => null()
    real(kind=jprb),     allocatable :: p_rt(:)
    integer :: ic
    integer :: nprof, nchan, nlev, nlev_
    integer :: ipr, ich
    logical :: l_opdep, ldeb

    if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
    l_opdep = present(opdep)
    if (l_opdep) then
      if (.not.present(p_l2c)) call finish(proc, 'argument opdep requires p_l2c')
    elseif (.not.present(transm)) then
      call finish(proc, 'Either opdep or transm required')
    end if

    rto   => rt_opts(iopt)
    ropts => rto%opts
    ic    =  rto%icoeff

#if defined(_RTTOV_GOD)
    if (associated(coefs(ic)%coef%god)) then
      if (l_opdep) then
        nprof = size(opdep, 3)
        nchan = size(opdep, 2)
        nlev  = size(opdep, 1) + 1
        allocate(p_rt(nlev))
        call rtifc_coef_prop(iopt, preslev=p_rt, nlevs=nlev_)
        if (nlev /= nlev_) call finish(proc, 'inconsistent nlev')
      else
        nprof = size(transm, 3)
        nchan = size(transm, 2)
      end if

      do ipr = 1, nprof
        if (present(ideb)) then
          ldeb = any(ipr == ideb)
        else
          ldeb = .false.
        end if
        do ich = 1, nchan
          if (ldeb) print*,'rtifc_l2c_god',ipr,ich,chans(ich), &
               any(coefs(ic)%coef%god(:,chans(ich))%ntr > 0)
          if (any(coefs(ic)%coef%god(:,chans(ich))%ntr > 0) .and. valid(ich,ipr)) then
            if (l_opdep) then
              call l2c_god_opd(opdep(:,ich,ipr), coefs(ic)%coef%god(:,chans(ich)), &
                               log(p_l2c(:,ipr)), log(p_rt(:)),                    &
                               l2c(ich,ipr), l2c_corr(ich,ipr), debug=ldeb)
            else
              call l2c_god_tau(transm(:,ich,ipr), coefs(ic)%coef%god(:,chans(ich)), &
                               l2c(ich,ipr), l2c_corr(ich,ipr), debug=ldeb)
            end if
            if (ldeb) print*,'rtifc_l2c_god result',ipr,ich,chans(ich), &
                 l2c(ich,ipr), l2c_corr(ich,ipr),l2c(ich,ipr)==l2c_corr(ich,ipr)
          else
            l2c_corr(ich,ipr) = l2c(ich,ipr)
          end if
        end do
      end do
    else
      l2c_corr = l2c
    end if
#else
    l2c_corr = l2c
#endif

  end subroutine rtifc_l2c_god


  function rtifc_coef_index(satid, platf, instr) result(idx)
    integer             :: idx
    integer, intent(in) :: satid
    integer, intent(in) :: platf
    integer, intent(in) :: instr

    integer :: k, nmatch

    idx = 0
    nmatch = 0
    do k = 1, size(coefs)
      if ( satid == coefs(k)%coef%id_sat      .and. &
           platf == coefs(k)%coef%id_platform .and. &
           instr == coefs(k)%coef%id_inst ) then
        idx = k
        nmatch = nmatch + 1
      end if
    end do
    select case (nmatch)
    case(0)
      idx = -ERR_NO_COEFS
    case(1)
      ! Fine: found exactly one coef structure
    case default
      ! Multiple coefs match. In this case you should not call this function. You should
      ! use the index given by the rtifc_init subroutine.
      idx = -ERR_MULTIPLE_COEFS
    end select

  end function rtifc_coef_index


  subroutine rtifc_cleanup(lprof, lcoef, latlas)
    logical, intent(in), optional :: lprof
    logical, intent(in), optional :: lcoef
    logical, intent(in), optional :: latlas
    !------------------------------------------------------
    ! frees memory allocated by rtifc_init/rtifc_fill_input
    !------------------------------------------------------
    character(len=13), parameter :: proc = 'rtifc_cleanup'
    integer :: i, k, stat
    logical :: l

    ! deallocation of rttov permanent arrays
    if (present(lprof)) then ; l = lprof ; else ; l = .true. ; endif
    if (l) then
      if (associated(profiles)) then
        call dealloc_rttov_arrays(stat, profs=profiles,rads=radiance,rads2=radiance2,transm=transmission)
        if (stat /= NO_ERROR) call finish(proc, 'failed to deallocate RTTOV profiles/rad./transm.')
      endif
      if (associated(profiles_k)) then
        call dealloc_rttov_arrays(stat, profs=profiles_k,rads=radiance_k,transm=transmission_k  )
        if (stat /= NO_ERROR) call finish(proc, 'failed to deallocate RTTOV-K profiles/rad./transm.')
      endif
    end if

    ! deallocate coefficients
    if (present(lcoef)) then ; l = lcoef ; else ; l = .true. ; endif
    if (l) then
      if ((associated(coefs))) then
        do i=1,size(coefs)
          call rttov_dealloc_coefs(stat, coefs(i))
          if (stat /= NO_ERROR) call finish(proc, 'failed to deallocate RTTOV coefficients')
        enddo
        deallocate(coefs)
      endif
    end if

   ! Cleanup mw emis atlases
#if defined(_RTTOV_ATLAS)
    if (present(latlas)) then ; l = latlas ; else ; l = .true. ; endif
    if (l) then
      if (associated(atlas)) then
        do k = 1, size(atlas)
          if (.not. atlas(k)% init) cycle
          call rttov_deallocate_emis_atlas(atlas(k))
        end do
        deallocate(atlas)
      end if
    endif
#endif
  end subroutine rtifc_cleanup


#if defined(_DACE_)
  subroutine rtifc_fill_input_rad (status,rad,iopts,ivect,istart,iend, pe, &
                                   wr_profs, wr_profs_fmt)
    integer,             intent(out)          :: status
    type(t_radv),        intent(in)           :: rad             ! derived type to store all the above info
    integer,             intent(in), target   :: iopts(:)        ! options index
    integer,             intent(in), optional :: ivect           ! Vectorization
                                                                 ! ivect =1 : vectorize profiles
                                                                 ! ivect/=1 : vectorize levels
    integer,             intent(in), optional :: istart          ! first profile (in rad) to be used
    integer,             intent(in), optional :: iend            ! last  profile (in rad) to be used
    integer,             intent(in), optional :: pe
    integer,             intent(in), optional :: wr_profs(:)
    character(len=*),    intent(in), optional :: wr_profs_fmt
    !------------------------------------------------
    ! this routine fills the RTTOV profile structures
    !------------------------------------------------
    character(len=20),   parameter :: proc   = 'rtifc_fill_input_rad'
    type(t_rtopts),      pointer   :: rto    => null()
    type(rttov_options), pointer   :: ropts  => null()
    type(rttov_options), target    :: opts_
    integer(jpim)                  :: stat
    integer(jpim)                  :: iprof, jprof, mpr
    integer(jpim)                  :: ivect_loc
    integer(jpim)                  :: i, j, iopt
    integer                        :: is, ie
    ! writing of profiles
    character(len=300):: fname             =  ""

    status = NO_ERROR

    if (present(pe)) then
      pe_ifc = pe
      pe_rt  = pe_ifc
    end if

    do i = 1, size(iopts)
      iopt = iopts(i)
      if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
      rto => rt_opts(iopt)
      ropts => rto%opts
      if (i == 1) then
        opts_ = ropts
      else
        call merge_opts(opts_, ropts)
      end if
    end do
    ropts => opts_

#define DIM_ERROR(text) status = ERR_DIM ; write(0,*) '*** dim_err '//text ; return

FTRACE_BEGIN('rtifc_fill_input_rad')

    ! Check whether required variables are available
    if (.not.associated(rad%t2m   )) then ; DIM_ERROR('rad%t2m'   ) ; endif
    if (.not.associated(rad%q2m   )) then ; DIM_ERROR('rad%q2m'   ) ; endif
    if (.not.associated(rad%ps_fg )) then ; DIM_ERROR('rad%ps_fg' ) ; endif
    if (.not.associated(rad%u10_fg)) then ; DIM_ERROR('rad%u10_fg') ; endif
    if (.not.associated(rad%v10_fg)) then ; DIM_ERROR('rad%v10_fg') ; endif
    if (.not.associated(rad%stype )) then ; DIM_ERROR('rad%stype' ) ; endif
    if (.not.associated(rad%ts_fg )) then ; DIM_ERROR('rad%ts_fg' ) ; endif
    if (.not.associated(rad%shgt  )) then ; DIM_ERROR('rad%shgt'  ) ; endif
    if (.not.associated(rad%stzen )) then ; DIM_ERROR('rad%stzen' ) ; endif
    if (.not.associated(rad%dlat  )) then ; DIM_ERROR('rad%dlat'  ) ; endif
    if (.not.associated(rad%dlon  )) then ; DIM_ERROR('rad%dlon'  ) ; endif
    if (.not.associated(rad%p     )) then ; DIM_ERROR('rad%p'     ) ; endif
    if (.not.associated(rad%t_fg  )) then ; DIM_ERROR('rad%t_fg'  ) ; endif
    if (.not.associated(rad%q_fg  )) then ; DIM_ERROR('rad%q_fg'  ) ; endif

    ! Derive dimensions
    if (present(istart)) then
      is = istart
    else
      is = lbound(rad%t_fg, 2)
    end if
    if (present(iend)) then
      ie = iend
    else
      ie = ubound(rad%t_fg, 2)
    end if
    nprof = ie - is + 1
    nlevs = size(rad%t_fg, 1)
    call check_nlevs(nlevs, nlevs_top, status)
    if (status /= NO_ERROR) then
      if (.not.opts_%interpolation%addinterp) then
        return
      else
        ! Model levels as input
        nlevs_top = 0
      end if
    end if
    nlay = nlevs + nlevs_top - 1

    ! Check dimensions
    !    required vars
    if (size(rad%t2m    ) < nprof) then ; DIM_ERROR('t2m (nprof)'   ) ; end if
    if (size(rad%q2m    ) < nprof) then ; DIM_ERROR('q2m (nprof)'   ) ; end if
    if (size(rad%ps_fg  ) < nprof) then ; DIM_ERROR('ps_fg (nprof)' ) ; end if
    if (size(rad%u10_fg ) < nprof) then ; DIM_ERROR('u10_fg (nprof)') ; end if
    if (size(rad%v10_fg ) < nprof) then ; DIM_ERROR('v10_fg (nprof)') ; end if
    if (size(rad%ts_fg  ) < nprof) then ; DIM_ERROR('ts_fg (nprof)' ) ; end if
    if (size(rad%stype,2) < nprof) then ; DIM_ERROR('stype (nprof)' ) ; end if
    if (size(rad%shgt, 2) < nprof) then ; DIM_ERROR('shgt (nprof)'  ) ; end if
    if (size(rad%stzen  ) < nprof) then ; DIM_ERROR('stzen (nprof)' ) ; end if
    if (size(rad%dlat   ) < nprof) then ; DIM_ERROR('dlat (nprof)'  ) ; end if
    if (size(rad%dlon   ) < nprof) then ; DIM_ERROR('dlon (nprof)'  ) ; end if
    if (size(rad%p,    1) < nlevs) then ; DIM_ERROR('p (nlevs)'     ) ; end if
    if (size(rad%t_fg, 1) < nlevs) then ; DIM_ERROR('t_fg (nlevs)'  ) ; end if
    if (size(rad%q_fg, 1) < nlevs) then ; DIM_ERROR('q_fg (nlevs)'  ) ; end if
    if (size(rad%t_fg, 2) < nprof) then ; DIM_ERROR('t_fg (nprof)'  ) ; end if
    if (size(rad%q_fg, 2) < nprof) then ; DIM_ERROR('q_fg (nprof)'  ) ; end if
    !   optional vars
    if (associated(rad% cld_top) .and. associated(rad% cld_frc)) then
      if (size(rad% cld_top) < nprof) then ; DIM_ERROR('cld_top (nprof)') ; end if
      if (size(rad% cld_frc) < nprof) then ; DIM_ERROR('cld_frc (nprof)') ; end if
    endif
    if (associated(rad% stazi)) then
      if (size(rad%stazi   ) < nprof) then ; DIM_ERROR('stazi (nprof)'  ) ; end if
    end if
    if (associated(rad% sunazi)) then
      if (size(rad% sunazi ) < nprof) then ; DIM_ERROR('sunazi (nprof)' ) ; end if
    endif
    if (associated(rad% sunzen)) then
      if (size(rad% sunzen ) < nprof) then ; DIM_ERROR('sunzen (nprof)' ) ; end if
    endif
    if (associated(rad% cld_fg) .and. associated(rad% cld_frc)) then
      if (size(rad% cld_fg, 3) < nprof) then ; DIM_ERROR('cld_fg (nprof)' ) ; end if
      if (size(rad% cld_fg, 1) < nlay ) then ; DIM_ERROR('cld_fg (nlay)'  ) ; end if
      if (size(rad% cld_frc,2) < nprof) then ; DIM_ERROR('cld_frc (nprof)') ; end if
      if (size(rad% cld_frc,1) < nlay ) then ; DIM_ERROR('cld_frc (nlay)' ) ; end if
    endif
    if (rad% i_o3 > 0) then
      if (.not.associated(rad%trg)   ) then ; DIM_ERROR('rad%trg (O3)'   ) ; end if
      if (size(rad%trg,1) < rad%i_o3) then ; DIM_ERROR('rad%trg (i_o3)' ) ; end if
      if (size(rad%trg,2) < nlevs  ) then ; DIM_ERROR('rad%trg (nlevs)') ; end if
      if (size(rad%trg,3) < nprof  ) then ; DIM_ERROR('rad%trg (nprof)') ; end if
    endif

    status = 0

#undef DIM_ERROR

    ! Check for inconsistent input
    if (associated(rad% cld_fg) .and. associated(rad% cld_frc)) then
       status = ERR_CLOUD_INCONS
       return
    endif
    if (rad%i_o3 > 0 .neqv. ropts%rt_ir%ozone_data) then
      status = ERR_TRACEGAS_INCONS
      return
    end if

    ! allocate the direct structures if not yet done
    call realloc_rttov_arrays(status, nprof, nlevs+nlevs_top, nprof*nmax_chans,&
                              ropts, profs=profiles)
    if (status /= NO_ERROR) return

    ! fill RTTOV profile input
    profiles(1:nprof)% nlevels    = nlevs + nlevs_top
    profiles(1:nprof)% nlayers    = nlevs + nlevs_top - 1
    profiles(1:nprof)% gas_units  = default_gas_units
    profiles(1:nprof)% mmr_cldaer = .false. ! g/cm^3 instead of kg/kg

    ! Scalar variables
    do iprof = is, ie
      jprof = iprof - is + 1
      if (associated(rad% obsnum)) then
        write(profiles(jprof)% id,*) rad% obsnum(iprof)
      end if

      ! 2 meter air variables:
      profiles(jprof)% s2m% t          = min(tmax_ifc,max(tmin_ifc,rad% t2m(iprof)))
      profiles(jprof)% s2m% q          = min(qmax_ifc,max(qmin_ifc,rad% q2m(iprof)))
      profiles(jprof)% s2m% p          = rad% ps_fg(iprof)
      profiles(jprof)% s2m% u          = rad% u10_fg(iprof)
      profiles(jprof)% s2m% v          = rad% v10_fg(iprof)
      profiles(jprof)% s2m% o          = default_o3_surf
      profiles(jprof)% s2m% wfetc      = default_wfetch

      ! skin variables
      profiles(jprof)% skin% t         = rad% ts_fg(iprof)
      profiles(jprof)% skin% snow_fraction = rad% snw_frc(iprof)
      profiles(jprof)% skin% surftype  = int(rad% stype(1,iprof), jpim)
      profiles(jprof)% skin% salinity  = default_salinity
      profiles(jprof)% skin% watertype = default_watertype
      profiles(jprof)% skin% fastem    = default_fastem

      ! cloud stuff
      if (associated(rad% cld_top) .and. associated(rad% cld_frc)) then
        profiles(jprof)% ctp           = rad% cld_top(iprof)
        profiles(jprof)% cfraction     = rad% cld_frc(1,iprof)
      else
        profiles(jprof)% ctp           = default_ctp
        profiles(jprof)% cfraction     = default_cfraction
      end if
      profiles(jprof)% idg             = default_idg
      profiles(jprof)% ice_scheme      = default_ice_scheme

      ! geometric variables
      profiles(jprof)% elevation  = rad% shgt (1,iprof)
      profiles(jprof)% zenangle   = abs(rad% stzen(iprof))
      profiles(jprof)% latitude   = rad% dlat (iprof)
      profiles(jprof)% longitude  = rad% dlon (iprof)
      if (associated(rad% stazi)) then
        profiles(jprof)% azangle = rad% stazi(iprof)
      else
        profiles(jprof)% azangle = default_satazim
      endif
      if (associated(rad% sunazi)) then
        if (abs(rad% sunazi(iprof)) <= 360._wp) then
          profiles(jprof)% sunazangle = rad% sunazi(iprof)
        else
          profiles(jprof)% sunazangle = default_sunazangle
        end if
      else
        profiles(jprof)% sunazangle = default_sunazangle
      endif
      if (associated(rad% sunzen)) then
        if (abs(rad% sunzen(iprof)) <= 360._wp) then
          profiles(jprof)% sunzenangle = rad% sunzen(iprof)
        else
          profiles(jprof)% sunzenangle = default_sunzenangle
        end if
      else
        profiles(jprof)% sunzenangle = default_sunzenangle
      endif
    enddo

    mpr = ubound(rad% p, 2)

    ! vectorization
    ivect_loc = 0
    if (present(ivect)) ivect_loc = ivect

    if (ivect_loc == 1) then ! Vectorize profiles
      if (nlevs_top==1) then
        do iprof = is, ie
          jprof = iprof - is + 1
          profiles(jprof)% p(1) = 0.005_jprb
          profiles(jprof)% t(1) = min(tmax_ifc,max(tmin_ifc,rad% t_fg(1,iprof)))
          profiles(jprof)% q(1) = min(qmax_ifc,max(qmin_ifc,rad% q_fg(1,iprof)))
        enddo
      end if

      do i=1,nlevs
        do iprof = is, ie
          j = min(iprof, mpr)
          jprof = iprof - is + 1
          profiles(jprof)% p(i+nlevs_top) = rad% p(i,j)
          profiles(jprof)% t(i+nlevs_top) = min(tmax_ifc,max(tmin_ifc,rad% t_fg(i,iprof)))
          profiles(jprof)% q(i+nlevs_top) = min(qmax_ifc,max(qmin_ifc,rad% q_fg(i,iprof)))
        enddo
      enddo

      if (associated(rad% cld_fg) .and. associated(rad% cld_frc)) then
!NEC$ novector
        do i=1,nlay
!NEC$ novector
          do j=1,size(rad% cld_fg,1)
            do iprof = is, ie
              jprof = iprof - is + 1
              profiles(jprof)% cloud(j,i) = rad% cld_fg(j,i,iprof)
            enddo
          enddo
        enddo
!NEC$ novector
        do i=1,nlay
          do iprof = is, ie
            jprof = iprof - is + 1
            profiles(jprof)% cfrac(i) = rad% cld_frc(i,iprof)
          enddo
        enddo
      endif

      if (associated(rad% clw)) then
        do iprof = is, ie
          profiles(iprof-is+1)% clw(1) = rad% clw(1,iprof)
        enddo
!NEC$ novector
        do i=1,nlevs
          do iprof = is, ie
            profiles(iprof-is+1)% clw(i+nlevs_top) = rad% clw(i,iprof)
          enddo
        enddo
      end if

      if (rad%i_o3 > 0) then
        do iprof = is, ie
          jprof = iprof - is + 1
          profiles(jprof)% o3(1) = rad%trg(rad%i_o3,1,iprof)
        enddo
!NEC$ novector
        do i=1,nlevs
          do iprof = is, ie
            jprof = iprof - is + 1
            profiles(jprof)% o3(i+nlevs_top) = rad%trg(rad%i_o3,i,iprof)
          enddo
        enddo
      endif

      if (rad%i_co2 > 0) then
        do iprof = is, ie
         jprof = iprof - is + 1
          profiles(jprof)% co2(1) = rad%trg(rad%i_co2,1,iprof)
        enddo
!NEC$ novector
        do i=1,nlevs
          do iprof = is, ie
           jprof = iprof - is + 1
            profiles(jprof)% co2(i+nlevs_top) = rad%trg(rad%i_co2,i,iprof)
          enddo
        enddo
      endif

!       ... similar for n2o, co, ch4, ...
    else ! Vectorize levels

      do iprof = is, ie
        jprof = iprof - is + 1
        if (associated(rad%obsnum)) then
          write(profiles(jprof)% id,*) rad%obsnum(iprof)
        end if
        if (nlevs_top==1) then
          profiles(jprof)% p(1)            = 0.005_jprb
          profiles(jprof)% t(1)            = min(tmax_ifc,max(tmin_ifc,rad%t_fg(1,iprof)))
          profiles(jprof)% q(1)            = min(qmax_ifc,max(qmin_ifc,rad%q_fg(1,iprof)))
        end if
        j = min(iprof, mpr)
        profiles(jprof)% p(1+nlevs_top:) = rad% p(1:nlevs,j)
        profiles(jprof)% t(1+nlevs_top:) = min(tmax_ifc,max(tmin_ifc,rad% t_fg(1:nlevs,iprof)))
        profiles(jprof)% q(1+nlevs_top:) = min(qmax_ifc,max(qmin_ifc,rad% q_fg(1:nlevs,iprof)))

        if (associated(rad% cld_fg) .and. associated(rad% cld_frc)) then
          profiles(jprof)% cloud(1:size(rad%cld_fg,1),1:nlay) = rad% cld_fg(:,1:nlay,iprof)
          profiles(jprof)% cfrac(1:nlay) = rad% cld_frc(1:nlay,iprof)
        end if

        if (associated(profiles(jprof)% clw)) then
          if (nlevs_top==1) profiles(jprof)% clw  (1) = rad% clw  (1,iprof)
          profiles(jprof)% clw(1+nlevs_top:) = rad% clw  (1:nlevs,iprof)
        end if
        if (rad%i_o3 > 0) then
          profiles(jprof)% o3(1)               = rad%trg(rad%i_o3,1,iprof)
          profiles(jprof)% o3(1+nlevs_top:)    = rad%trg(rad%i_o3,1:nlevs,iprof)
        endif
        if (rad%i_co2 > 0) then
          profiles(jprof)% co2(1)               = rad%trg(rad%i_co2,1,iprof)
          profiles(jprof)% co2(1+nlevs_top:)    = rad%trg(rad%i_co2,1:nlevs,iprof)
        endif
!       ... similar for n2o, co, ch4, ...
      end do
    endif

    ! Missing t_surf/T_G in input to var3d/mec results in skin%t == 0. This will crash
    ! in RTTOV if the do_checkinput option (apply_reg_lims option in var3d/mec) is not set.
    ! Thus, we check here for invalid skin temp. in order to get a useful error message.
    if (any(profiles(1:nprof)% skin% t < tmin_ifc)) then
      status = ERR_INVALID_TSKIN
      return
    end if

    if (present(wr_profs)) then
      if (size(wr_profs) > 0) then
        do iprof = 1, nprof
          if (any(rad%obsnum(iprof) == wr_profs(:))) then
            write(fname,trim(wr_profs_fmt)) rad%obsnum(iprof)
            call rttov_write_profile(profiles(iprof:iprof), ropts, trim(fname), stat)
            if (stat /= 0) then
              status = ERR_WR_PROF
              return
            end if
          end if
        end do
      end if
    end if

FTRACE_END('rtifc_fill_input_rad')

  end subroutine rtifc_fill_input_rad
#endif


  subroutine rtifc_fill_input_var (status,iopts,press,temp,humi,t2m,q2m,psurf,hsurf,&
                                   u10m,v10m,stemp,stype,lat,lon,sat_zen,sun_zen,   &
                                   sat_azi,sun_azi,cloud,cfrac,ctp,cfraction,       &
                                   ice_scheme, idg, watertype,snw_frc,id,           &
                                   ivect,istart,iend,pe,wr_profs, wr_profs_fmt)
    integer,             intent(out)          :: status
    integer,             intent(in)           :: iopts       (:) ! options index
    real(wp),            intent(in)           :: press     (:,:) ! press. grid (nlevs,nprofs) [hPa]
    real(wp),            intent(in)           :: temp      (:,:) ! temperature profiles (nlevs,nprofs)[K]
    real(wp),            intent(in)           :: humi      (:,:) ! water vapour profile (nlevs,nprofs) [ppmv]
    real(wp),            intent(in)           :: t2m         (:) ! 2m temperature (nprofs) [K]
    real(wp),            intent(in)           :: q2m         (:) ! 2m humidity (nprofs) [ppmv]
    real(wp),            intent(in)           :: psurf       (:) ! surface pressure (nprofs) [hPa]
    real(wp),            intent(in)           :: hsurf       (:) ! surface height (nprofs) [km]
    real(wp),            intent(in)           :: u10m        (:) ! 10m U wind component (nprofs) [m/s]
    real(wp),            intent(in)           :: v10m        (:) ! 10m V wind component (nprofs) [m/s]
    real(wp),            intent(in)           :: stemp       (:) ! radiative skin temperature (nprofs) [K]
    integer,             intent(in)           :: stype       (:) ! surface type (nprofs):0,1,2=land,sea,sea-ice
    real(wp),            intent(in)           :: lat         (:) ! latitude of ground point (nprofs) [deg]
    real(wp),            intent(in)           :: lon         (:) ! longitude (nprofs) [deg]
    real(wp),            intent(in)           :: sat_zen     (:) ! local satellite zenith angle (nprofs) [deg]
    real(wp),            intent(in), optional :: sun_zen     (:) ! local solar zenith angle (nprofs) [deg]
    real(wp),            intent(in), optional :: ctp         (:) ! cloud top pressure (of a black body cloud)
                                                                 !   (nprofs) [hPa] default: 500
    real(wp),            intent(in), optional :: cfraction   (:) ! cloud fraction (of a black body cloud)
    real(wp),            intent(in), optional :: sat_azi     (:) ! local satellite azimuth (nprofs) [deg]
                                                                 !   range: 0-360; east=90
    real(wp),            intent(in), optional :: sun_azi     (:) ! local solar azimuth angle (nprofs) [deg]
                                                                 ! range: 0-360; east=90
    real(wp),            intent(in), optional :: cloud   (:,:,:) ! cloud water/ice (IR only)
                                                                 ! (nCloudPar,nlevs,nprofs) [g m^-3]
    real(wp),            intent(in), optional :: cfrac     (:,:) ! cloud fractional cover (IR only)
                                                                 ! (nCloudPar,nlevs,nprofs) range (0-1)
    integer,             intent(in), optional :: ice_scheme  (:) !
    integer,             intent(in), optional :: idg         (:) !
    integer,             intent(in), optional :: watertype   (:) !
    real(wp),            intent(in), optional :: snw_frc     (:) ! snow fraction
    integer,             intent(in), optional :: id          (:) ! Profile id
    integer,             intent(in), optional :: ivect           ! Vectorization
                                                                 ! ivect =1 : vectorize profiles
                                                                 ! ivect/=1 : vectorize levels
    integer,             intent(in), optional :: istart          ! first profile (in rad) to be used
    integer,             intent(in), optional :: iend            ! last  profile (in rad) to be used
    integer,             intent(in), optional :: pe
    integer,             intent(in), optional :: wr_profs(:)
    character(len=*),    intent(in), optional :: wr_profs_fmt
    !------------------------------------------------
    ! this routine fills the RTTOV profile structures
    !------------------------------------------------
    character(len=20),   parameter :: proc   = 'rtifc_fill_input_var'
    type(t_rtopts),      pointer   :: rto    => null()
    type(rttov_options), pointer   :: ropts  => null()
    type(rttov_options), target    :: opts_
    integer(jpim)     :: stat
    integer(jpim)     :: iprof, jprof, mpr
    integer(jpim)     :: ivect_loc
    integer(jpim)     :: i, j, iopt
    integer           :: is, ie
    ! writing of profiles
    character(len=300):: fname             =  ""

    status = NO_ERROR

    if (present(pe)) then
      pe_ifc = pe
      pe_rt  = pe_ifc
    end if

    do i = 1, size(iopts)
      iopt = iopts(i)
      if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
      rto => rt_opts(iopt)
      ropts => rto%opts
      if (i == 1) then
        opts_ = ropts
      else
        call merge_opts(opts_, ropts)
      end if
    end do
    ropts => opts_

#define DIM_ERROR(text) status = ERR_DIM ; write(0,*) 'dim_error: '//text ; return

FTRACE_BEGIN('rtifc_fill_input_var')

    ! Derive dimensions
    if (present(istart)) then
      is = istart
    else
      is = lbound(temp, 2)
    end if
    if (present(iend)) then
      ie = iend
    else
      ie = ubound(temp, 2)
    end if
    nprof     = ie - is + 1
    nlevs     = size(temp, 1)
    call check_nlevs(nlevs, nlevs_top, status)
    if (status /= NO_ERROR) then
      if (.not.ropts%interpolation%addinterp) then
        return
      else
        ! Model levels as input
        nlevs_top = 0
      end if
    end if
    nlay = nlevs + nlevs_top - 1

    ! Check dimensions
    !    required vars
    if (size(t2m    ) < nprof) then ; DIM_ERROR('t2m (nprof)'    ) ; end if
    if (size(q2m    ) < nprof) then ; DIM_ERROR('q2m (nprof)'    ) ; end if
    if (size(psurf  ) < nprof) then ; DIM_ERROR('psurf (nprof)'  ) ; end if
    if (size(u10m   ) < nprof) then ; DIM_ERROR('u10m (nprof)'   ) ; end if
    if (size(v10m   ) < nprof) then ; DIM_ERROR('v10m (nprof)'   ) ; end if
    if (size(stemp  ) < nprof) then ; DIM_ERROR('stemp (nprof)'  ) ; end if
    if (size(stype  ) < nprof) then ; DIM_ERROR('stype (nprof)'  ) ; end if
    if (size(hsurf  ) < nprof) then ; DIM_ERROR('hsurf(nprof)'   ) ; end if
    if (size(sat_zen) < nprof) then ; DIM_ERROR('sat_zen (nprof)') ; end if
    if (size(lat    ) < nprof) then ; DIM_ERROR('lat (nprof)'    ) ; end if
    if (size(lon    ) < nprof) then ; DIM_ERROR('lon (nprof)'    ) ; end if
    if (size(press,1) < nlevs) then ; DIM_ERROR('press (nlevs)'  ) ; end if
    if (size(temp, 1) < nlevs) then ; DIM_ERROR('temp (nlevs)'   ) ; end if
    if (size(humi, 1) < nlevs) then ; DIM_ERROR('humi (nlevs)'   ) ; end if
    if (size(temp, 2) < nprof) then ; DIM_ERROR('temp (nprof)'   ) ; end if
    if (size(humi, 2) < nprof) then ; DIM_ERROR('humi (nprof)'   ) ; end if
    !   optional vars
    if (present(ctp) .and. present(cfraction)) then
      if (size(ctp)       < nprof) then ; DIM_ERROR('ctp (nprof)'      ) ; end if
      if (size(cfraction) < nprof) then ; DIM_ERROR('cfraction (nprof)') ; end if
    endif
    if (present(sat_azi)) then
      if (size(sat_azi) < nprof) then ; DIM_ERROR('sat_azi (nprof)') ; end if
    end if
    if (present(sun_azi)) then
      if (size(sun_azi) < nprof) then ; DIM_ERROR('sun_azi (nprof)') ; end if
    endif
    if (present(sun_zen)) then
      if (size(sun_zen) < nprof) then ; DIM_ERROR('sun_zen (nprof)') ; end if
    endif
    if (present(cloud) .and. present(cfrac)) then
      if (size(cloud, 3) < nprof) then ; DIM_ERROR('cloud (nprof)') ; end if
      if (size(cloud, 2) < nlay ) then ; DIM_ERROR('cloud (nlay)' ) ; end if
      if (size(cloud, 1) < maxval(coefs(:)%coef_scatt_ir%fmv_wcl_comp+1)) &
                                   then ; DIM_ERROR('cloud (ncld)' ) ; end if
      if (size(cfrac, 2) < nprof) then ; DIM_ERROR('cfrac (nprof)') ; end if
      if (size(cfrac, 1) < nlay ) then ; DIM_ERROR('cfrac (nlay)' ) ; end if
    endif
    if (present(ice_scheme)) then
      if (size(ice_scheme) < nprof) then ; DIM_ERROR('ice_scheme (nprof)') ; end if
    endif
    if (present(idg)) then
      if (size(idg) < nprof) then ; DIM_ERROR('idg (nprof)') ; end if
    endif
    if (present(watertype)) then
      if (size(watertype) < nprof) then ; DIM_ERROR('watertype (nprof)') ; end if
    endif

    status = 0

#undef DIM_ERROR

    if (present(cloud) .neqv. present(cfrac)) then
       status = ERR_CLOUD_INCONS
       return
    endif

    ! allocate the direct structures if not yet done
    call realloc_rttov_arrays(status, nprof, nlevs+nlevs_top, nprof*nmax_chans,&
                              ropts, profs=profiles)
    if (status /= NO_ERROR) return

    ! fill RTTOV profile input
    profiles(1:nprof)% nlevels = nlevs + nlevs_top
    profiles(1:nprof)% nlayers = nlevs + nlevs_top - 1
    profiles(1:nprof)% gas_units  = default_gas_units
    profiles(1:nprof)% mmr_cldaer = .false. ! g/cm^3 instead of kg/kg

    ! Scalar variables
    do iprof = is, ie
      jprof = iprof - is + 1
      if (present(id)) then
        write(profiles(jprof)% id,*) id(iprof)
      end if
      ! 2 meter air variables:
      profiles(jprof)% s2m% t          = min(tmax_ifc,max(tmin_ifc,t2m(iprof)))
      profiles(jprof)% s2m% q          = min(qmax_ifc,max(qmin_ifc,q2m(iprof)))
      profiles(jprof)% s2m% p          = psurf(iprof)
      profiles(jprof)% s2m% u          = u10m(iprof)
      profiles(jprof)% s2m% v          = v10m(iprof)
      profiles(jprof)% s2m% o          = default_o3_surf
      profiles(jprof)% s2m% wfetc      = default_wfetch

      ! skin variables
      profiles(jprof)% skin% t         = stemp(iprof)
      profiles(jprof)% skin% surftype  = int(stype(iprof), jpim)
      profiles(jprof)% skin% salinity  = default_salinity
      profiles(jprof)% skin% fastem    = default_fastem
      if (present(snw_frc)) then
        profiles(jprof)% skin% snow_fraction = snw_frc(iprof)
      else
        profiles(jprof)% skin% snow_fraction = 0._wp
      end if
       if (present(watertype)) then
        profiles(jprof)%skin%watertype = watertype(iprof)
      else
        profiles(jprof)%skin%watertype = default_watertype
      end if

      ! cloud stuff
      if (present(ctp) .and. present(cfraction)) then
        profiles(jprof)% ctp           = ctp(iprof)
        profiles(jprof)% cfraction     = cfraction(iprof)
      else
        profiles(jprof)% ctp           = default_ctp
        profiles(jprof)% cfraction     = default_cfraction
      end if
      if (present(idg)) then
        profiles(jprof)% idg           = idg(iprof)
      else
        profiles(jprof)% idg           = default_idg
      end if
      if (present(ice_scheme)) then
        profiles(jprof)% ice_scheme    = ice_scheme(iprof)
      else
        profiles(jprof)% ice_scheme    = default_ice_scheme
      end if

      ! geometric variables
      profiles(jprof)% elevation  = hsurf(iprof)
      profiles(jprof)% zenangle   = abs(sat_zen(iprof))
      profiles(jprof)% latitude   = lat(iprof)
      profiles(jprof)% longitude  = lon(iprof)
      if (present(sat_azi)) then
        profiles(jprof)% azangle = sat_azi(iprof)
      else
        profiles(jprof)% azangle = default_satazim
      endif
      if (present(sun_azi)) then
        if (abs(sun_azi(iprof)) <= 360._wp) then
          profiles(jprof)% sunazangle = sun_azi(iprof)
        else
          profiles(jprof)% sunazangle = default_sunazangle
        end if
      else
        profiles(jprof)% sunazangle = default_sunazangle
      endif
      if (present(sun_zen)) then
        if (abs(sun_zen(iprof)) <= 360._wp) then
          profiles(jprof)% sunzenangle = sun_zen(iprof)
        else
          profiles(jprof)% sunzenangle = default_sunzenangle
        end if
      else
        profiles(jprof)% sunzenangle = default_sunzenangle
      endif
    enddo

    mpr = ubound(press, 2)

    ! vectorization
    ivect_loc = 0
    if (present(ivect)) ivect_loc = ivect

    if (ivect_loc == 1) then ! Vectorize profiles
      if (nlevs_top==1) then
        do iprof = is, ie
          jprof = iprof - is + 1
          profiles(jprof)% p(1) = 0.005_jprb
          profiles(jprof)% t(1) = min(tmax_ifc,max(tmin_ifc,temp(1,iprof)))
          profiles(jprof)% q(1) = min(qmax_ifc,max(qmin_ifc,humi(1,iprof)))
        enddo
      end if

      do i=1,nlevs
        do iprof = is, ie
          j = min(iprof, mpr)
          jprof = iprof - is + 1
          profiles(jprof)% p(i+nlevs_top) = press(i,j)
          profiles(jprof)% t(i+nlevs_top) = min(tmax_ifc,max(tmin_ifc,temp(i,iprof)))
          profiles(jprof)% q(i+nlevs_top) = min(qmax_ifc,max(qmin_ifc,humi(i,iprof)))
        enddo
      enddo

      if (present(cloud) .and. present(cfrac)) then
!NEC$ novector
        do i=1,nlay
!NEC$ novector
          do j=1,size(cloud,1)
            do iprof = is, ie
              jprof = iprof - is + 1
              profiles(jprof)% cloud(j,i) = cloud(j,i,iprof)
            enddo
          enddo
        enddo
!NEC$ novector
        do i=1,nlay
          do iprof = is, ie
            jprof = iprof - is + 1
            profiles(jprof)% cfrac(i) = cfrac(i,iprof)
          enddo
        enddo
      endif

    else ! Vectorize levels

      do iprof = is, ie
        jprof = iprof - is + 1
        if (nlevs_top==1) then
          profiles(jprof)% p(1)            = 0.005_jprb
          profiles(jprof)% t(1)            = min(tmax_ifc,max(tmin_ifc,temp(1,iprof)))
          profiles(jprof)% q(1)            = min(qmax_ifc,max(qmin_ifc,humi(1,iprof)))
        end if
        j = min(iprof, mpr)
        profiles(jprof)% p(1+nlevs_top:) = press(1:nlevs,j)
        profiles(jprof)% t(1+nlevs_top:) = min(tmax_ifc,max(tmin_ifc,temp(1:nlevs,iprof)))
        profiles(jprof)% q(1+nlevs_top:) = min(qmax_ifc,max(qmin_ifc,humi(1:nlevs,iprof)))
        if (present(cloud) .and. present(cfrac)) then
          profiles(jprof)% cloud(1:size(cloud,1),1:nlay) = cloud(:,1:nlay,iprof)
          profiles(jprof)% cfrac(1:nlay) = cfrac(1:nlay,iprof)
        end if
      end do
    endif

    ! Missing t_surf/T_G in input to var3d/mec results in skin%t == 0. This will crash
    ! in RTTOV if the do_checkinput option (apply_reg_lims option in var3d/mec) is not set.
    ! Thus, we check here for invalid skin temp. in order to get a useful error message.
    if (any(profiles(1:nprof)% skin% t < tmin_ifc)) then
      status = ERR_INVALID_TSKIN
      return
    end if

    if (present(wr_profs) .and. present(id)) then
      if (size(wr_profs) > 0) then
        do iprof = 1, nprof
          if (any(id(iprof) == wr_profs(:))) then
            write(fname,trim(wr_profs_fmt)) id(iprof)
            call rttov_write_profile(profiles(iprof:iprof), ropts, trim(fname), stat)
            if (stat /= 0) then
              status = ERR_WR_PROF
              return
            end if
          end if
        end do
      end if
    end if

FTRACE_END('rtifc_fill_input_var')

  end subroutine rtifc_fill_input_var


  subroutine rttov_write_profile(prof, ropts, name, stat)
    use hdf5
    use rttov_hdf_mod
    type(rttov_profile), intent(in)  :: prof(:)
    type(rttov_options), intent(in)  :: ropts
    character(len=*),    intent(in)  :: name
    integer,             intent(out) :: stat

    call open_hdf(.true., stat)
    if (stat /= 0) return
    call rttov_hdf_save(stat, trim(name), '/PROFILES', create=.true., profiles=prof(:))
    if (stat /= 0) return
    call rttov_hdf_save(stat, trim(name), '/OPTIONS', create=.false., options=ropts)
    if (stat /= 0) return
    call close_hdf(stat)
    if (stat /= 0) return

  end subroutine rttov_write_profile


  subroutine rtifc_direct (iopt,lprofs,chans,emissiv,t_b,status,                   &
                           sat_zen, sat_azi, t_b_clear,rad,radclear,radupclear,    &
                           raddnclear,refdnclear,radovercast,radtotal,             &
                           transm,transmtotal,opdep,height, istore,reg_lim, rflag, &
                           dealloc,iprint,rad_out_flg,pe, l_pio)
    integer,             intent(in)            :: iopt               ! options index
    integer,             intent(in)            :: lprofs         (:) ! list of profile indices (nchans*nprof)
    integer,             intent(in)            :: chans          (:) ! list of channel indices (nchans*nprof)
    real(wp),            intent(inout)         :: emissiv      (:,:) ! emissivities (nchans,nprof),
                                                                     ! if < 0.01, they will be calculated,
                                                                     ! else they will be used
    real(wp),            intent(out)           :: t_b          (:,:) ! brightness temp. (nchans,nprof) [K]
    integer,             intent(out)           :: status             ! exit status
    real(wp),            intent(in),  optional :: sat_zen        (:) ! satellite zenith angle
    real(wp),            intent(in),  optional :: sat_azi        (:) ! satellite azimuth angle
                                                                     ! (The above sat_* variables are useful here in order to use the same
                                                                     !  profile for different satellites, e.g. for synsat calc.)
    real(wp),            intent(out), optional :: t_b_clear    (:,:) ! calc. clear sky b.t. (nchans,nprof) [K]
    real(wp),            intent(out), optional :: rad          (:,:) ! calculated total radiance
                                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: radclear     (:,:) ! calculated clear sky radiances
                                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: radupclear   (:,:) ! calculated clear sky upwelling radiances
                                                                     ! includes surface emission term but not downwelling reflected
                                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: raddnclear   (:,:) ! calculated clear sky downwelling radiances at surface
                                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: refdnclear   (:,:) ! reflected clear sky downwelling radiances at TOA
                                                                     ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: radovercast(:,:,:) ! overcast radiances (nlevs,nchans,nprof)
                                                                     ! [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: radtotal     (:,:) ! total radiances (nchans,nprof)
                                                                     ! [mW/cm^-1/ster/m^2]
    real(wp),            intent(out), optional :: transm     (:,:,:) ! transmission (nlevs,nchans,nprof)
    real(wp),            intent(out), optional :: transmtotal  (:,:) ! total surface to TOA transmission (nchans,nprof)
    real(wp),            intent(out), optional :: opdep      (:,:,:) ! original optical depth (nlevs,nchans,nprof)
    real(wp),            intent(out), optional :: height       (:,:) ! height estimated on basis of weighting function
    integer,             intent(in),  optional :: istore       (:,:) ! put result i into array position (istore(i,1),istore(i,2))
    integer,             intent(out), optional :: reg_lim    (:,:,:) ! result of apply_reg_limits (nlevs,nvars,nprof)
    integer,             intent(out), optional :: rflag        (:,:) ! results of other test like e.g. the chk_god test
    integer,             intent(in),  optional :: iprint         (:) ! profile to printed in RTTOV (debug)
    logical,             intent(in),  optional :: dealloc
    integer,             intent(in),  optional :: rad_out_flg        ! bit field (see OUT_* parameters)
    integer,             intent(in),  optional :: pe
    logical,             intent(in),  optional :: l_pio

    !-----------------------------------------------------------------------
    ! this subroutine renews the instrument dependent part
    ! of a call of rttov_direct in the case that there are more
    ! than one instrument looking at the same ground point
    ! from the same satellite (or at least if their measurements
    ! are interpolated to the same ground point by e.g., aapp)
    !-----------------------------------------------------------------------
    character(len=12), parameter :: proc   = 'rtifc_direct'
    type(t_rtopts),     pointer  :: rto    => null()
    type(rttov_options),pointer  :: ropts  => null()
    integer                      :: ic
    integer                      :: ipr(5),npr
    integer                      :: rad_out
    integer(jpim)                :: iprof
    integer(jpim)                :: nlevs_coef
    integer(jpim)                :: nchans
    integer(jpim)                :: nchansprofs
    integer(jpim)                :: nprof_store
    integer(jpim)                :: nprof_calc
    integer(jpim)                :: iprof_s
    integer(jpim)                :: iprof_e
    integer(jpim)                :: sind
    integer(jpim)                :: eind
    integer(jpim)                :: i, k, minchan, maxchan, minprof, maxprof
    integer(jpim)                :: stat
    integer(jpim)                :: lprofs_aux    (size(lprofs))
    integer(jpim),   allocatable :: index_prof    (:)
    logical(jplm),   allocatable :: prof_used     (:)
    logical                      :: l_dealloc
    logical(jplm)                :: calcemis      (size(chans))
    type(rttov_emissivity)       :: emissivity    (size(chans))
    integer(jpim)                :: errstat
    type(rttov_chanprof)         :: chanprof(size(chans))
    integer                      :: lims_flag(nlevs,2)

    real(kind=jprb),     pointer :: height_aux(:) => NULL() ! height of weighting function

    logical                      :: lpio
    logical                      :: l_opdep, l_transm, l_radoverc
    logical                      :: l_chk_god

    type(rttov_options)          :: ropts_
    type(rttov_profile), pointer :: p => null()
    character(len=1000)          :: msg  = ''
    character(len=1000)          :: msg_ = ''

FTRACE_BEGIN('rtifc_direct')

    status = NO_ERROR

    ! Process arguments
    lpio = .false.
    if (present(l_pio)) lpio=l_pio

    if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
    rto   => rt_opts(iopt)
    ropts => rto%opts
    ic    =  rto%icoeff

    l_transm   = .false. ; if (present(transm     )) l_transm   = (size(transm)      > 0)
    l_opdep    = .false. ; if (present(opdep      )) l_opdep    = (size(opdep )      > 0)
    l_radoverc = .false. ; if (present(radovercast)) l_radoverc = (size(radovercast) > 0)

    npr = 0 ; ipr = 0
    if (present(iprint)) then
      npr = count(iprint > 0)
      ipr(1:npr) = pack(iprint, mask=(iprint>0))
    end if
    if (npr > 0) then
      ipr_deb = ipr(1)
      do i = 1, npr
        write(msg,*) 'rttov input profile ',ipr(i),trim(profiles(ipr(i))%id)
        call rttov_print_profile(profiles(ipr(i)), usd, trim(msg))
        call rttov_print_opts(ropts, usd, 'ropts (options for rttov_direct call):')
      end do
    else
      ipr_deb = 0
    end if

    if (present(pe)) then
      pe_ifc = pe
      pe_rt  = pe_ifc
    end if

    if (present(rad_out_flg)) then
      rad_out = rad_out_flg
    else
      rad_out = 0
      if (present(t_b_clear)) rad_out = ibset(rad_out,OUT_CSB)
      if (present(rad      )) rad_out = ibset(rad_out,OUT_ASR)
      if (present(radclear ) .or. present(radupclear) .or. &
          present(raddnclear) .or. present(refdnclear)) rad_out = ibset(rad_out,OUT_CSR)
    end if

    if (present(dealloc)) then
      l_dealloc = dealloc
    else
      l_dealloc = .true.
    end if

    ! Dimensions
    nchansprofs    = size(chans)
    if (present(istore)) then
       nchans      = maxval(istore(1:nchansprofs,1))
       nprof_store = maxval(istore(1:nchansprofs,2))
    else
       nchans      = nchansprofs/nprof
       nprof_store = nprof
    end if
    nlay       = profiles(1)% nlayers
    nlevs_coef = coefs(ic)%coef%nlevels

    ! Prepare profiles array
    ! RTTOV(10) crashes, if the number of profiles supplied to it is bigger than
    ! the number of forward-calculations. This can happen if e.g.
    ! chans = (/ 1, 2, 1, 2 /) and lprofs = (/ 1, 1, 7, 7 /)  (7 profiles and 4
    ! forward calculations. Therefore the profiles array has to be packed in the
    ! RTTOV call: pack(profiles(iprof_s:iprof_e), mask=prof_used)
    ! In the following an appropriate lprofs-array "lprofs_aux" is derived.
    ! 1. Determine the used profiles
    iprof_s     = minval(lprofs(:))
    iprof_e     = maxval(lprofs(:))
    allocate(prof_used(iprof_s:iprof_e))
    prof_used(:) = .false.
    do i = 1, nchansprofs
      prof_used(lprofs(i)) = .true.
    end do
    nprof_calc  = count(prof_used(:))
    if (nchansprofs < nprof_calc) then
      write(0,*) 'nchansprofs',nchansprofs
      write(0,*) 'nprof_calc',nprof_calc
      write(0,*) 'iprof_s, iprof_e', iprof_s, iprof_e
      write(0,*) 'lprofs', lprofs(:)
      call finish(proc, 'Less RTTOV-calculations than profiles! RTTOV will crash')
    end if
    ! 2. Determine an array to translate lprofs into lprofs_aux
    allocate(index_prof(iprof_s:iprof_e))
    k = 1
    do i = iprof_s, iprof_e
      if (prof_used(i)) then
        index_prof(i) = k
        k = k + 1
        if (i==ipr_deb) ipr_deb = index_prof(i)
      end if
    end do
    ! 3. Translate lprofs into lprofs_aux
    do i = 1, nchansprofs
      lprofs_aux(i) = index_prof(lprofs(i))
    end do
    deallocate(index_prof)

    ! Prepare storing
    if (.not.present(istore)) then
      sind = 1
      eind = nchans
      do iprof=1,nprof
        emissivity(sind:eind)%emis_in  = emissiv(1:nchans,iprof)
        emissivity(sind:eind)%emis_out = 0._jprb
        sind = eind + 1
        eind = eind + nchans
      enddo
    else
!NEC$ ivdep
      do i = 1, nchansprofs
        emissivity(i)%emis_in  = emissiv(istore(i,1),istore(i,2))
        emissivity(i)%emis_out = 0._jprb
      end do
    end if

    ! Check dimensions
#define DIM_ERROR(text) status = ERR_DIM ; write(0,*) 'dim_error: '//text ; return
    if (.not.present(istore)) then
      ! check dimensions of input and output arrays:
      if (size(lprofs ) /= nchansprofs) then
        DIM_ERROR('lprofs')
      endif
      if (size(emissiv) /= nchansprofs) then
        DIM_ERROR('emissiv')
      endif
      if ((size(t_b,1) /= nchans)  .or. &
           (size(t_b,2) /= nprof_store )) then
        DIM_ERROR('t_b')
      endif
      if (btest(rad_out,OUT_CSB)) then
        if (present(t_b_clear)) then
          if ((size(t_b_clear,1) /= nchans) .or. &
               (size(t_b_clear,2) /= nprof_store )) then
            DIM_ERROR('t_b_clear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_ASR)) then
        if (present(rad   )) then
          if ((size(rad,1) /= nchans ) .or. &
               (size(rad,2) /= nprof_store  )) then
            DIM_ERROR('rad')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radclear)) then
          if ((size(radclear,1) /= nchans) .or. &
               (size(radclear,2) /= nprof_store )) then
            DIM_ERROR('radclear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radupclear)) then
          if ((size(radupclear,1) /= nchans) .or. &
               (size(radupclear,2) /= nprof_store )) then
            DIM_ERROR('radupclear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(raddnclear)) then
          if ((size(raddnclear,1) /= nchans) .or. &
               (size(raddnclear,2) /= nprof_store )) then
            DIM_ERROR('raddnclear')
          endif
        endif
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(refdnclear)) then
          if ((size(refdnclear,1) /= nchans) .or. &
               (size(refdnclear,2) /= nprof_store )) then
            DIM_ERROR('refdnclear')
          endif
        endif
      end if
      if (l_radoverc) then
        if ( (size(radovercast,1) /= nlay       ) .or. &
             (size(radovercast,2) <  nchans     ) .or. &
             (size(radovercast,3) <  nprof_store)) then
          DIM_ERROR('radovercast')
        endif
      endif
      if (l_opdep) then
        if ( (size(opdep,1) /= nlevs_coef-1) .or. &
             (size(opdep,2) <  nchans     ) .or. &
             (size(opdep,3) <  nprof_store)) then
          DIM_ERROR('opdep')
        endif
      endif
      if (l_transm) then
        if ( (size(transm,1) /= nlevs      ) .or. &
             (size(transm,2) <  nchans     ) .or. &
             (size(transm,3) <  nprof_store)) then
          DIM_ERROR('transm')
        endif
      endif
      if (present(radtotal)) then
        if ( (size(radtotal,1) < nchans     ) .or. &
             (size(radtotal,2) < nprof_store)) then
          DIM_ERROR('radtotal')
        endif
      endif
       if (present(transmtotal)) then
        if ((size(transmtotal,1) /= nchans)  .or. &
             (size(transmtotal,2) /= nprof_store )) then
          DIM_ERROR('transmtotal')
        endif
      endif
      if (present(height)) then
        if ((size(height,1) /= nchans)  .or. &
             (size(height,2) /= nprof_store )) then
          DIM_ERROR('height')
        endif
      endif
    else
      minchan = minval(istore(1:nchansprofs,1))
      maxchan = maxval(istore(1:nchansprofs,1))
      minprof = minval(istore(1:nchansprofs,2))
      maxprof = maxval(istore(1:nchansprofs,2))
      if ((minchan < lbound(t_b, 1)) .or. (maxchan > ubound(t_b, 1))) then
        DIM_ERROR('t_b(1)')
      endif
      if ((minprof < lbound(t_b, 2)) .or. (maxprof > ubound(t_b, 2))) then
        DIM_ERROR('t_b(2)')
      endif
      if (btest(rad_out,OUT_CSB)) then
        if (present(t_b_clear)) then
          if ((minchan < lbound(t_b_clear, 1)) .or. (maxchan > ubound(t_b_clear, 1))) then
            DIM_ERROR('t_b_clear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_ASR)) then
        if (present(rad)) then
          if ((minchan < lbound(rad, 1)) .or. (maxchan > ubound(rad, 1))) then
            DIM_ERROR('rad')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radclear)) then
          if ((minchan < lbound(radclear, 1)) .or. (maxchan > ubound(radclear, 1))) then
            DIM_ERROR('radclear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(radupclear)) then
          if ((minchan < lbound(radupclear, 1)) .or. (maxchan > ubound(radupclear, 1))) then
            DIM_ERROR('radupclear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(raddnclear)) then
          if ((minchan < lbound(raddnclear, 1)) .or. (maxchan > ubound(raddnclear, 1))) then
            DIM_ERROR('raddnclear')
          endif
        end if
      end if
      if (btest(rad_out,OUT_CSR)) then
        if (present(refdnclear)) then
          if ((minchan < lbound(refdnclear, 1)) .or. (maxchan > ubound(refdnclear, 1))) then
            DIM_ERROR('refdnclear')
          endif
        end if
      end if
      if (l_radoverc) then
        if ((minchan < lbound(radovercast, 2)) .or. (maxchan > ubound(radovercast, 2))) then
          DIM_ERROR('radovercast')
        endif
      end if
      if (l_transm) then
        if ((minchan < lbound(transm, 2)) .or. (maxchan > ubound(transm, 2))) then
          DIM_ERROR('transm')
        endif
      end if
      if (l_opdep) then
        if ((minchan < lbound(opdep, 2)) .or. (maxchan > ubound(opdep, 2))) then
          DIM_ERROR('opdep')
        endif
      end if
      if (present(radtotal)) then
        if ((minchan < lbound(radtotal, 1)) .or. (maxchan > ubound(radtotal, 1))) then
          DIM_ERROR('radtotal')
        endif
      end if
      if (present(transmtotal)) then
        if ((minchan < lbound(transmtotal, 1)) .or. (maxchan > ubound(transmtotal, 1))) then
          DIM_ERROR('transmtotal')
        endif
      end if
      if (present(height)) then
        if ((minchan < lbound(height, 1)) .or. (maxchan > ubound(height, 1))) then
          DIM_ERROR('height')
        endif
      end if
    endif
    if (chk_reg_lims /= 0) then
      if (present(reg_lim)) then
        if ( size(reg_lim,1) < nlevs       .or. &
             size(reg_lim,2) < 2           .or. &
             size(reg_lim,3) < nprof_store) then
          DIM_ERROR('reg_lim')
        endif
      else
        DIM_ERROR('reg_lim not present')
      endif
    end if

    ! Put satellite direction into profiles
    if (present(sat_azi)) then
       if (size(sat_azi(:)) < nprof) then
         DIM_ERROR('sat_azi')
       endif
       do i = 1, nprof
          profiles(i)% azangle  = sat_azi(i)
       enddo
    endif
    if (present(sat_zen)) then
       if (size(sat_zen(:)) < nprof) then
         DIM_ERROR('sat_zen')
       endif
       do i = 1, nprof
          profiles(i)% zenangle = sat_zen(i)
       enddo
    endif

#undef DIM_ERROR


    ! Allocate/initialize arrays
#if defined(_DACE_)
    l_chk_god = (chk_god /= 0 .and. nprof > 0 .and. &
                 present(rflag) .and. associated(coefs(ic)%coef%god))
    ! for l_chk_god the original opdep (on rt-levels) is only required if addinterp=T,
    ! otherwise the transmission is used.
    transmission%l_opdep = l_opdep .or. (l_chk_god .and. ropts%interpolation%addinterp)
#endif
    call realloc_rttov_arrays(status, nprof, nlevs + nlevs_top, nchansprofs, &
                              ropts, rads=radiance, rads2=radiance2, transm=transmission,&
                              n_levels_coef=nlevs_coef)
    if (status /= NO_ERROR) return

    ! clear variable output part of rttov structures
    call rttov_clear_rad_var    (nchansprofs,radiance)
    call rttov_clear_rad2_var   (nchansprofs,radiance2)
    call rttov_init_transmission(transmission)

    ! Prepare calcemis.
    calcemis(1:nchansprofs) = emissivity(1:nchansprofs)%emis_in < 0.01_jprb

    ! fill chanprof array
    do i=1,nchansprofs
      chanprof(i)% chan = chans(i)
      chanprof(i)% prof = lprofs_aux(i)
    enddo

    ! Debug
    if (npr > 0) then
      do i = 1, nchansprofs
        if (any(chanprof(i)% prof == ipr(1:npr))) then
          print*,'rttov emissivity input: ',i,chanprof(i)%chan,emissivity(i),calcemis(i)
        end if
      end do
    end if

    ! call RTTOV
    errstat = errorstatus_success
#if defined(_RTTOV_USE_OPENMP)
    call rttov_parallel_direct (                            &
           errstat,                                         & ! --> error flag
           chanprof(1:nchansprofs),                         & ! <-- channels and profiles to calculate
           ropts,                                           & ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used), & ! <--  profile array
           coefs(ic),                                       & ! <--  coefs array
           transmission,                                    & ! --> array of transmittances
           radiance,                                        & ! <-> computed radiance array
           radiance2,                                       & ! <--> computed secondary radiance array
           calcemis=calcemis(1:nchansprofs),                & ! <-- flag for internal emissivity calc
           emissivity=emissivity(1:nchansprofs),            & ! <-- input emissivities per channel
           nthreads=omp_get_max_threads()                   )
#else
    call rttov_direct (                                     &
           errstat,                                         & ! --> error flag
           chanprof(1:nchansprofs),                         & ! <-- channels and profiles to calculate
           ropts,                                           & ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used), & ! <--  profile array
           coefs(ic),                                       & ! <--  coefs array
           transmission,                                    & ! --> array of transmittances
           radiance,                                        & ! <-> computed radiance array
           radiance2,                                       & ! <--> computed secondary radiance array
           calcemis=calcemis(1:nchansprofs),                & ! <-- flag for internal emissivity calc
           emissivity=emissivity(1:nchansprofs)            )  ! <-- input emissivities per channel
#endif
    if (errstat == ERRORSTATUS_FATAL) then
       print *, 'after rttov direct -> fatal error'
       status = ERR_RTTOV_CALL
       return
    endif

    ! Calculate height
    if (present(height)) then
#if defined(_DACE_)
      call realloc_rttov_arrays(status, nprof, nlevs+nlevs_top, nchansprofs, &
                                ropts, height=height_aux)
      if (status /= NO_ERROR) return
      do i = 1, nchansprofs
        call get_weighting_function(profiles(chanprof(i)%prof)%p(1:nlevs), &
                                    profiles(chanprof(i)%prof)%t(1:nlevs), &
                                    transmission%tau_levels(:,i),          &
                                    height=height_aux(i))
        if (stat /= NO_ERROR) then
          status = stat
          return
        endif
      end do
#else
      status = ERR_INPUT
      return
#endif
    end if

    ! fill result into output arrays
    if (.not.present(istore)) then
      sind = 1
      eind = nchans
      do iprof=1,nprof
        emissiv(1:nchans,iprof) = dble(emissivity(sind:eind)%emis_out)
        t_b(1:nchans,iprof) = dble(radiance% bt(sind:eind))
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear)) t_b_clear(1:nchans,iprof)=dble(radiance%bt_clear(sind:eind))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad      )) rad      (1:nchans,iprof)=dble(radiance%total   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear )) radclear (1:nchans,iprof)=dble(radiance%clear   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radupclear )) radupclear (1:nchans,iprof)=dble(radiance2%upclear   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(raddnclear )) raddnclear (1:nchans,iprof)=dble(radiance2%dnclear   (sind:eind))
          if (present(refdnclear )) refdnclear (1:nchans,iprof)=dble(radiance2%refldnclear   (sind:eind))
        end if
        if (present(radtotal))  radtotal (1:nchans,iprof)=dble(radiance%clear   (sind:eind)) !clear -> total ???
        if (l_radoverc) &
             radovercast(:,1:nchans,iprof) = dble(radiance% overcast(:,sind:eind))
        if (present(height))    height   (1:nchans,iprof)=height_aux(sind:eind)
        if (l_transm) &
             transm(:,1:nchans,iprof) = dble(transmission%tau_levels(1+nlevs_top:,sind:eind))
        if (present(transmtotal)) transmtotal(1:nchans,iprof) = dble(transmission%tau_total(sind:eind))
#if defined(_DACE_) && !defined(__ICON__)
        if (l_opdep) &
             opdep (:,1:nchans,iprof) = dble(transmission%opdep_ref (1+nlevs_top:,sind:eind))
#endif
        sind = eind + 1
        eind = eind + nchans
      enddo
    else
!NEC$ ivdep
      do i = 1, nchansprofs
        emissiv(istore(i,1),istore(i,2)) = dble(emissivity(i)%emis_out)
        t_b(istore(i,1),istore(i,2)) = dble(radiance% bt(i))
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear  )) t_b_clear  (istore(i,1),istore(i,2))   = dble(radiance% bt_clear(i))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad        )) rad        (istore(i,1),istore(i,2))   = dble(radiance% total   (i))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear   )) radclear   (istore(i,1),istore(i,2))   = dble(radiance% clear   (i))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radupclear   )) radupclear   (istore(i,1),istore(i,2))   = dble(radiance2% upclear  (i))
          if (present(raddnclear   )) raddnclear   (istore(i,1),istore(i,2))   = dble(radiance2% dnclear  (i))
          if (present(refdnclear   )) refdnclear   (istore(i,1),istore(i,2))   = dble(radiance2% refldnclear  (i))
        end if
        if (present(radtotal   )) radtotal    (istore(i,1),istore(i,2))   = dble(radiance% clear   (i))
        if (l_radoverc          ) radovercast (:,istore(i,1),istore(i,2)) = dble(radiance% overcast(:,i))
        if (present(height     )) height      (istore(i,1),istore(i,2))   = height_aux(i)
        if (l_transm            ) transm(:,istore(i,1),istore(i,2)) = dble(transmission%tau_levels(1+nlevs_top:,i))
        if (present(transmtotal)) transmtotal(istore(i,1),istore(i,2)) = dble(transmission%tau_total(i))
#if defined(_DACE_) && !defined(__ICON__)
        if (l_opdep) &
             opdep(:,istore(i,1),istore(i,2)) = dble(transmission%opdep_ref(1+nlevs_top:,i))
#endif
      end do
    end if

    ! Check profile limits
    if (chk_reg_lims /= 0 .and. nprof > 0) then
      call realloc_rttov_arrays(status, nprof_aux,profiles(1)%nlevels,0,ropts,profs=profiles_aux)
      if (status /= NO_ERROR) return
      ropts_ = ropts
      do iprof = iprof_s, iprof_e
        if (prof_used(iprof)) then
          call check_prof(iprof, coefs(ic), lims_flag, ropts_, lpio)
          if (chk_reg_lims > 1 .and. present(reg_lim)) then
            if (.not.present(istore)) then
              do i = 1, nprof
                if (lprofs((i-1)*nchans+1) ==iprof) then
                  reg_lim(:,:,i) = lims_flag(:,:)
                  exit
                end if
              end do
            else
              do i = 1, nchansprofs
                if (lprofs(i)==iprof) then
                  reg_lim(:,:,istore(i,2)) = lims_flag(:,:)
                  exit
                end if
              end do
            end if
          end if
        end if
      end do
      call dealloc_rttov_arrays(status, profs=profiles_aux)
      if (status /= NO_ERROR) return
    end if

    ! Check impact of god-smoothing
#if defined(_RTTOV_GOD)
    if (l_chk_god) then
      do i = 1, nchansprofs
        stat = 0
        if (any(lprofs(i) == ipr(1:npr))) then
          print*,'check_god_ifc prof',trim(profiles(lprofs(i))%id),chans(i),coefs(ic)%coef%ff_ori_chn(chans(i))
          print*,'check_god_ifc ntr',chans(i),any(coefs(ic)%coef%god(:,chans(i))%ntr > 0)
        end if
        if (any(coefs(ic)%coef%god(:,chans(i))%ntr > 0)) then
          if (transmission%l_opdep) then
            call check_god_infl(coefs(ic)%coef%god(:,chans(i)), stat, msg_, &
                                od_ref=transmission%opdep_ref(:,i), debug=any(lprofs(i)==ipr(1:npr)))
          else
            call check_god_infl(coefs(ic)%coef%god(:,chans(i)), stat, msg_, &
                                tau=transmission%tau_levels(:,i), debug=any(lprofs(i)==ipr(1:npr)))
          end if
          if (stat /= 0) then
            if (btest(chk_god, 0)) then
              p => profiles(lprofs(i))
              write(msg,'(2x,A,I5.5," god-smoothing too large in ",A," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p%id)//' chan ',coefs(ic)%coef%ff_ori_chn(chans(i)),&
                   trim(msg_),p%latitude,p%longitude
              WR
            end if
          end if
        end if
        if ((btest(chk_god,1) .or. btest(chk_god,2)) .and. all(shape(rflag) > 0)) then
          if (present(istore)) then
            if (any(lprofs(i) == ipr(1:npr))) then
              print*,'check_god_ifc stat',stat,istore(i,:)
            end if
            rflag(istore(i,1), istore(i,2)) = stat
          else
            ! TODO !!!
          end if
        end if
      end do
    elseif (present(rflag)) then
      rflag(:,:) = 0
    end if
#endif /* _RTTOV_GOD */

    ! Trailer
    if (l_dealloc) then
      call dealloc_rttov_arrays(status,rads=radiance, rads2=radiance2,transm=transmission)
      if (status /= NO_ERROR) return
    end if

FTRACE_END('rtifc_direct')

  end subroutine rtifc_direct


  subroutine rtifc_k (iopt,lprofs,chans,emissiv,emissiv_k,temp_k,humi_k,          &
                      t2m_k,q2m_k,stemp_k,t_b,status,sat_zen,sat_azi,t_b_clear,   &
                      rad,radclear,radtotal,radovercast,transm,opdep,             &
                      psurf_k,u10m_k,v10m_k,o3_surf_k,wfetc_k,ctp_k,cfraction_k,  &
                      clw_k,o3_k,co2_k,n2o_k,co_k,ch4_k,istore,reg_lim,rflag,     &
                      dealloc,iprint,rad_out_flg,pe,l_pio)
    integer,             intent(in)          :: iopt               ! options index
    integer,             intent(in)          :: lprofs         (:) ! list of profile indices(nchans*nprof)
    integer,             intent(in)          :: chans          (:) ! list of channel indices(nchans*nprof)
    real(wp),            intent(inout)       :: emissiv      (:,:) ! emissivities - calculated if < 0.01,
                                                                   !   else they will be used (nchans*nprof)
    real(wp),            intent(inout)       :: emissiv_k    (:,:) ! k-matrix of surf. emiss. (nchans*nprof)
    real(wp),            intent(out)         :: temp_k     (:,:,:) ! grad. of temp. (nlevs,nchans,nprof) [K/K]
    real(wp),            intent(out)         :: humi_k     (:,:,:) ! grad. of w.v. (nlevs,nchans,nprof) [K/ppmv]
    real(wp),            intent(out)         :: t2m_k        (:,:) ! grad. of 2m temp. (nchans,nprof) [K/K]
    real(wp),            intent(out)         :: q2m_k        (:,:) ! grad. of 2m hum. (nchans,nprof) [K/ppmv]
    real(wp),            intent(out)         :: stemp_k      (:,:) ! grad. of surf. skin t. (nchans,nprof) [K/K]
    real(wp),            intent(out)         :: t_b          (:,:) ! calculated b.t. (nchans,nprof) [K]
    integer,             intent(out)         :: status             ! exit status
    real(wp),            intent(in),optional :: sat_zen        (:) ! satellite zenith angle
    real(wp),            intent(in),optional :: sat_azi        (:) ! satellite azimuth angle
                                                                   ! (These variables are useful here in order to use the same
                                                                   !  profile for different satellites, e.g. for synsat calc.)
    real(wp),            intent(out),optional:: t_b_clear    (:,:) ! calc. clear sky b.t. (nchans,nprof) [K]
    real(wp),            intent(out),optional:: rad          (:,:) ! calculated total radiance
                                                                   ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out),optional:: radclear     (:,:) ! calculated clear sky radiances
                                                                   ! (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out),optional:: radtotal     (:,:) ! calculated total radiances
                                                                   !   (nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out),optional:: radovercast(:,:,:) ! calculated overcast radiances
                                                                   !   (nlevs,nchans,nprof) [mW/cm^-1/ster/m^2]
    real(wp),            intent(out),optional:: transm     (:,:,:) ! transmission (nlevs,nchans,nprof)
    real(wp),            intent(out),optional:: opdep      (:,:,:) ! transmission (nlevs,nchans,nprof)
    real(wp),            intent(out),optional:: psurf_k      (:,:) ! grad. of surf. press. (nchans,nprof) [K/hPa]
    real(wp),            intent(out),optional:: u10m_k       (:,:) ! grad. of 10m U wind (nchans,nprof) [K/(m/s)]
    real(wp),            intent(out),optional:: v10m_k       (:,:) ! grad. of 10m V wind (nchans,nprof) [K/(m/s)]
    real(wp),            intent(out),optional:: o3_surf_k    (:,:) ! k  of surf. o3 conc. (nchans,nprof) [K/ppmv]
    real(wp),            intent(out),optional:: wfetc_k      (:,:) ! grad. of wind fetch (nchans,nprof) [K/m]
    real(wp),            intent(out),optional:: ctp_k        (:,:) ! grad. w.resp. cloud top pressure   [K/hPa]
    real(wp),            intent(out),optional:: cfraction_k  (:,:) ! grad. w.resp. cloud fraction       [K]
    real(wp),            intent(out),optional:: clw_k      (:,:,:) ! grad. of cloud liquid water (microwave only)
                                                                   !   (nlevs,nchans,nprof) [K/(kg/kg)]
    real(wp),            intent(out),optional:: o3_k       (:,:,:) ! grad. of o3  (nlevs,nchans,nprof) [K/ppmv]
    real(wp),            intent(out),optional:: co2_k      (:,:,:) ! grad. of co2 (nlevs,nchans,nprof) [K/ppmv]
    real(wp),            intent(out),optional:: n2o_k      (:,:,:) ! grad. of n2o (nlevs,nchans,nprof) [K/ppmv]
    real(wp),            intent(out),optional:: co_k       (:,:,:) ! grad. of co  (nlevs,nchans,nprof) [K/ppmv]
    real(wp),            intent(out),optional:: ch4_k      (:,:,:) ! grad. of ch4 (nlevs,nchans,nprof) [K/ppmv]
    integer,             intent(in), optional:: istore       (:,:) ! put result i into array position (istore(i,1),istore(i,2))
    integer,             intent(out),optional:: reg_lim    (:,:,:) ! result of apply_reg_limits (nlevs,nvars,nprof)
    integer,             intent(out),optional:: rflag        (:,:) ! results of other test like e.g. the chk_god test
    logical,             intent(in), optional:: dealloc
    integer,             intent(in), optional:: rad_out_flg        ! bit field (see OUT_* parameters)
    integer,             intent(in), optional:: iprint         (:)
    integer,             intent(in), optional:: pe
    logical,             intent(in), optional:: l_pio

    !-----------------------------------------------------------------------
    ! this subroutine renews the instrument dependent part
    ! of a call of rttov_k in the case that there are more
    ! than one instrument looking at the same ground point
    ! from the same satellite (or at least if their measurements
    ! are interpolated to the same ground point by e.g., aapp)
    !----------------------------------------------------------------------
    character(len=7),  parameter :: proc   = 'rtifc_k'
    type(t_rtopts),     pointer  :: rto    => null()
    type(rttov_options),pointer  :: ropts  => null()
    integer                      :: ic
    integer                      :: rad_out
    integer                      :: pe_loc
    integer                      :: ipr(5),npr
    integer(jpim)                :: iprof
    integer(jpim)                :: ichan
    integer(jpim)                :: i, k
    integer(jpim)                :: stat
    integer(jpim)                :: nlevs_coef
    integer(jpim)                :: nchans
    integer(jpim)                :: nchansprofs
    integer(jpim)                :: nprof_store
    integer(jpim)                :: nprof_calc
    integer(jpim)                :: iprof_s
    integer(jpim)                :: iprof_e
    integer(jpim)                :: sind
    integer(jpim)                :: eind
    integer(jpim)                :: count1
    integer(jpim)                :: nusedchans
    integer(jpim)                :: lprofs_aux      (size(lprofs))
    integer(jpim),   allocatable :: index_prof      (:)
    logical(jplm),   allocatable :: prof_used       (:)
    logical                      :: l_dealloc
    logical(jplm)                :: calcemis        (size(chans ))
    type(rttov_emissivity)       :: emissivity      (size(chans ))
    type(rttov_emissivity)       :: emissivity_k    (size(chans ))
    integer(jpim)                :: errstat
    type(rttov_chanprof)         :: chanprof(size(chans))
    integer                      :: lims_flag(nlevs,2)
    logical                      :: lpio
    logical                      :: l_opdep, l_transm, l_radoverc
    logical                      :: l_chk_god
    type(rttov_options)          :: ropts_
    type(rttov_profile), pointer :: p => null()
    character(len=1000)          :: msg  = ''
    character(len=1000)          :: msg_ = ''


FTRACE_BEGIN('rtifc_k')

    status = NO_ERROR

    ! Process arguments
    lpio = .false.
    if (present(l_pio)) lpio=l_pio

    if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
    rto   => rt_opts(iopt)
    ropts => rto%opts
    ic    =  rto%icoeff

    l_transm   = .false. ; if (present(transm     )) l_transm   = (size(transm)      > 0)
    l_opdep    = .false. ; if (present(opdep      )) l_opdep    = (size(opdep )      > 0)
    l_radoverc = .false. ; if (present(radovercast)) l_radoverc = (size(radovercast) > 0)

    if (present(pe)) then
      pe_loc = pe
    else
      pe_loc = -1
    end if

    if (present(rad_out_flg)) then
      rad_out = rad_out_flg
    else
      rad_out = ibset(0,OUT_ASB)
      if (present(t_b_clear)) rad_out = ibset(rad_out,OUT_CSB)
      if (present(rad      )) rad_out = ibset(rad_out,OUT_ASR)
      if (present(radclear )) rad_out = ibset(rad_out,OUT_CSR)
     end if

    npr = 0 ; ipr = 0
    if (present(iprint)) then
      npr = count(iprint > 0)
      ipr(1:npr) = pack(iprint, mask=(iprint>0))
    end if
    if (npr > 0) then
      ipr_deb = ipr(1)
      do i = 1, npr
        write(msg,*) 'rttov input profile ',ipr(i),trim(profiles(ipr(i))%id)
        call rttov_print_profile(profiles(ipr(i)), usd, trim(msg))
        call rttov_print_opts(ropts, usd, 'ropts (options for rttov_direct call):')
      end do
    else
      ipr_deb = 0
    end if

    if (present(dealloc)) then
      l_dealloc = dealloc
    else
      l_dealloc = .true.
    end if

    ! Dimensions
    nchansprofs  = size(chans)
    if (present(istore)) then
      nchans      = maxval(istore(1:nchansprofs,1))
      nprof_store = maxval(istore(1:nchansprofs,2))
    else
      nchans      = nchansprofs / nprof
      nprof_store = nprof
    end if
    nlay       = profiles(1)% nlayers
    nlevs_coef = coefs(ic)%coef%nlevels

    ! Prepare lprofs array
    ! RTTOV(10) crashes, if the number of profiles supplied to it is bigger than
    ! the number of forward-calculations. This can happen if e.g.
    ! chans = (/ 1, 2, 1, 2 /) and lprofs = (/ 1, 1, 7, 7 /)  (7 profiles and 4
    ! forward calculations. Therefore the profiles array has to be packed in the
    ! RTTOV call: pack(profiles(iprof_s:iprof_e), mask=prof_used)
    ! In the following an appropriate lprofs-array "lprofs_aux" is derived.
    ! 1. Determine the used profiles
    iprof_s     = minval(lprofs(:))
    iprof_e     = maxval(lprofs(:))
    allocate(prof_used(iprof_s:iprof_e))
    prof_used(:) = .false.
    do i = 1, nchansprofs
      prof_used(lprofs(i)) = .true.
    end do
    nprof_calc  = count(prof_used(:))
    if (nchansprofs < nprof_calc) then
      write(0,*) 'nchansprofs',nchansprofs
      write(0,*) 'nprof_calc',nprof_calc
      write(0,*) 'iprof_s, iprof_e', iprof_s, iprof_e
      write(0,*) 'lprofs', lprofs(:)
      call finish('proc', 'Less RTTOV-calculations than profiles! RTTOV will crash')
    end if
    ! 2. Determine an array to translate lprofs into lprofs_aux
    allocate(index_prof(iprof_s:iprof_e))
    k = 1
    do i = iprof_s, iprof_e
      if (prof_used(i)) then
        index_prof(i) = k
        k = k + 1
        if (i==ipr_deb) ipr_deb = index_prof(i)
      end if
    end do
    ! 3. Translate lprofs into lprofs_aux
    do i = 1, nchansprofs
      lprofs_aux(i) = index_prof(lprofs(i))
    end do
    deallocate(index_prof)

    ! Prepare emissivity
    if (.not.present(istore)) then
      sind = 1
      eind = nchans
      do iprof=1,nprof
        emissivity  (sind:eind)%emis_in  = real (emissiv  (1:nchans,iprof),jprb)
        emissivity_k(sind:eind)%emis_in  = 0.0_jprb
        emissivity  (sind:eind)%emis_out = 0.0_jprb
        emissivity_k(sind:eind)%emis_out = 0.0_jprb
        sind = eind + 1
        eind = eind + nchans
      enddo
    else
!NEC$ ivdep
      do i = 1, nchansprofs
        emissivity  (i)%emis_in  = real (emissiv  (istore(i,1),istore(i,2)),jprb)
        emissivity_k(i)%emis_in  = 0.0_jprb
        emissivity  (i)%emis_out = 0.0_jprb
        emissivity_k(i)%emis_out = 0.0_jprb
      end do
    end if

    ! Check dimensions
#define DIM_ERROR(text) status = ERR_DIM ; write(0,*) 'dim_error: '//text ; return
    ! check dimensions of input and output arrays
    if ((size(temp_k,1) /= nlevs) .or. &
        (size(humi_k,1) /= nlevs)) then
      DIM_ERROR('temp_k/humi_k levels')
    endif
    if ((size(temp_k, 2) < nchans) .or. &
        (size(humi_k, 2) < nchans) .or. &
        (size(t2m_k,  1) < nchans) .or. &
        (size(q2m_k,  1) < nchans) .or. &
        (size(t_b,    1) < nchans) .or. &
        (size(stemp_k,1) < nchans)) then
      DIM_ERROR('*_k channels')
    endif
    if ((size(temp_k, 3) <  nprof_store) .or. &
        (size(humi_k, 3) <  nprof_store) .or. &
        (size(t2m_k,  2) <  nprof_store) .or. &
        (size(q2m_k,  2) <  nprof_store) .or. &
        (size(t_b,    2) <  nprof_store) .or. &
        (size(stemp_k,2) <  nprof_store)) then
      DIM_ERROR('*_k,t_b profiles')
    endif
    if ((size(lprofs   ) <  nchansprofs) .or. &
        (size(emissiv  ) <  nchansprofs) .or. &
        (size(emissiv_k) <  nchansprofs)) then
      DIM_ERROR('lprofs,emissiv* nchansprofs')
    endif
    if (present(radtotal)) then
       if ((size(radtotal,1) <  nchans) .or. &
           (size(radtotal,2) <  nprof_store )) then
         DIM_ERROR('radtotal')
       endif
    endif
    if (l_radoverc) then
      if ( (size(radovercast,1) /= nlay       ) .or. &
           (size(radovercast,2) <  nchans     ) .or. &
           (size(radovercast,3) <  nprof_store)) then
        DIM_ERROR('radovercast')
      endif
    endif
    if (l_opdep) then
      if ( (size(opdep,1) /= nlevs_coef-1) .or. &
           (size(opdep,2) <  nchans     ) .or. &
           (size(opdep,3) <  nprof_store)) then
        write(0,*) shape(opdep), nlevs_coef, nchans, nprof_store
        DIM_ERROR('opdep')
      endif
    endif
    if (l_transm) then
      if ( (size(transm,1) /= nlevs      ) .or. &
           (size(transm,2) <  nchans     ) .or. &
           (size(transm,3) <  nprof_store)) then
        DIM_ERROR('transm')
      endif
    endif
    if (present(psurf_k)) then
       if ((size(psurf_k,1) <  nchans) .or. &
           (size(psurf_k,2) <  nprof_store )) then
         DIM_ERROR('psurf_k')
       endif
    endif
    if (present(u10m_k)) then
       if ((size(u10m_k,1) <  nchans) .or. &
           (size(u10m_k,2) <  nprof_store )) then
         DIM_ERROR('u10m_k')
       endif
    endif
    if (present(v10m_k)) then
       if ((size(v10m_k,1) <  nchans) .or. &
           (size(v10m_k,2) <  nprof_store )) then
         DIM_ERROR('v10m_k')
       endif
    endif
    if (present(o3_surf_k)) then
       if ((size(o3_surf_k,1) <  nchans) .or. &
           (size(o3_surf_k,2) <  nprof_store )) then
         DIM_ERROR('o3_surf_k')
       endif
    endif
    if (present(wfetc_k)) then
       if ((size(wfetc_k,1) <  nchans) .or. &
           (size(wfetc_k,2) <  nprof_store )) then
         DIM_ERROR('wfetc_k')
       endif
    endif
    if (present(ctp_k)) then
       if ((size(ctp_k,1) <  nchans) .or. &
           (size(ctp_k,2) <  nprof_store )) then
         DIM_ERROR('ctp_k')
       endif
    endif
    if (present(cfraction_k)) then
       if ((size(cfraction_k,1) <  nchans) .or. &
           (size(cfraction_k,2) <  nprof_store )) then
         DIM_ERROR('cfraction_k')
       endif
    endif
    if (present(clw_k)) then
       if ((size(clw_k,1) /= nlevs ) .or. &
           (size(clw_k,2) <  nchans) .or. &
           (size(clw_k,3) <  nprof_store )) then
         DIM_ERROR('clw_k')
       endif
    endif
    if (present(o3_k)) then
       if ((size(o3_k,1) /= nlevs ) .or. &
           (size(o3_k,2) <  nchans) .or. &
           (size(o3_k,3) <  nprof_store )) then
         DIM_ERROR('o3_k')
       endif
    endif
    if (present(co2_k)) then
       if ((size(co2_k,1) /= nlevs ) .or. &
           (size(co2_k,2) <  nchans) .or. &
           (size(co2_k,3) <  nprof_store )) then
         DIM_ERROR('co2_k')
       endif
    endif
    if (present(n2o_k)) then
       if ((size(n2o_k,1) /= nlevs ) .or. &
           (size(n2o_k,2) <  nchans) .or. &
           (size(n2o_k,3) <  nprof_store )) then
         DIM_ERROR('n2o_k')
       endif
    endif
    if (present(co_k)) then
       if ((size(co_k,1) /= nlevs ) .or. &
           (size(co_k,2) <  nchans) .or. &
           (size(co_k,3) <  nprof_store )) then
         DIM_ERROR('co_k')
       endif
    endif
    if (present(ch4_k)) then
       if ((size(ch4_k,1) /= nlevs ) .or. &
           (size(ch4_k,2) <  nchans) .or. &
           (size(ch4_k,3) <  nprof_store )) then
         DIM_ERROR('ch4_k')
       endif
    endif
    if (chk_reg_lims /= 0) then
      if (present(reg_lim)) then
        if (size(reg_lim,1) < nlevs       .or. &
            size(reg_lim,2) < 2           .or. &
            size(reg_lim,3) < nprof_store) then
          DIM_ERROR('reg_lim')
        endif
      else
        DIM_ERROR('reg_lim not present')
      endif
    end if

    ! Put satellite direction into profiles
    if (present(sat_azi)) then
       if (size(sat_azi(:)) < nprof) then
         DIM_ERROR('sat_azi')
       endif
       do i = 1, nprof
          profiles(i)% azangle  = sat_azi(i)
       enddo
    endif
    if (present(sat_zen)) then
       if (size(sat_zen(:)) < nprof) then
         DIM_ERROR('sat_zen')
       endif
       do i = 1, nprof
          profiles(i)% zenangle = sat_zen(i)
       enddo
    endif

#undef DIM_ERROR

    ! Allocate/initialize arrays
#if defined(_RTTOV_GOD)
    l_chk_god = (chk_god /= 0 .and. nprof > 0 .and. &
                 present(rflag) .and. associated(coefs(ic)%coef%god))
    ! for l_chk_god the original opdep (on rt-levels) is only required if addinterp=T,
    ! otherwise the transmission is used.
    transmission%l_opdep = l_opdep .or. (l_chk_god .and. ropts%interpolation%addinterp)
    transmission_k%l_opdep = .false.
#endif
    ! Allocate/initilize arrays
    ! ATTENTION: k - matrix has requires as much profiles per
    !            instrument as channels exist
    call realloc_rttov_arrays(status, nchansprofs, nlevs+nlevs_top, nchansprofs, ropts, &
                              profs=profiles_k, rads=radiance_k, transm=transmission_k)
    if (status /= NO_ERROR) return
    call realloc_rttov_arrays(status, nchansprofs, nlevs+nlevs_top, nchansprofs, ropts, &
                              rads=radiance, transm=transmission,&
                              n_levels_coef=nlevs_coef)
    if (status /= NO_ERROR) return
    call rttov_clear_rad_var(nchansprofs,radiance)
    call rttov_clear_prof   (profiles_k(1:nchansprofs))
    call rttov_init_transmission(transmission)
    call rttov_init_transmission(transmission_k)
    call rttov_clear_rad_var(nchansprofs,radiance_k)

    profiles_k(1:nchansprofs)% gas_units = default_gas_units
    calcemis(1:nchansprofs) = emissivity(1:nchansprofs)%emis_in < 0.01_jprb
    if (ropts%rt_all%switchrad) then
      radiance_k% bt(1:nchansprofs)    = 1.0_jprb
      radiance_k% total(1:nchansprofs) = 0.0_jprb
    else
      radiance_k% bt(1:nchansprofs)    = 0.0_jprb
      radiance_k% total(1:nchansprofs) = 1.0_jprb
    endif

    ! fill chanprof array
    do i=1,nchansprofs
      chanprof(i)% chan = chans(i)
      chanprof(i)% prof = lprofs_aux(i)
    enddo

    ! Debug
    if (npr > 0) then
      do i = 1, nchansprofs
        if (any(chanprof(i)% prof == ipr(1:npr))) then
          print*,'rttov emissivity input: ',i,chanprof(i),emissivity(i),calcemis(i)
        end if
      end do
    end if

    ! call RTTOV
    errstat = errorstatus_success
#if defined(_RTTOV_USE_OPENMP)
    call rttov_parallel_k(                                  &
           errstat,                                         & !  --> error flag
           chanprof(1:nchansprofs),                         & ! <-- channels and profiles to calculate
           ropts,                                           & ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used), & ! <--  profile array
           profiles_k(1:nchansprofs),                       & ! <--> profile increments;
           coefs(ic),                                       & ! <--  coefs array
           transmission,                                    & ! <--> transmittances
           transmission_k,                                  & ! <--> K matrix of transmittances
           radiance,                                        & ! <--> direct model output radiances
           radiance_k,                                      & ! <--> K matrix of radiances
           calcemis(1:nchansprofs),                         & ! <--  flag for internal emissivity calc
           emissivity(1:nchansprofs),                       & ! <-- input emissivities per channel
           emissivity_k(1:nchansprofs),                     & ! <--> k matrix on input surface emissivity
           nthreads=omp_get_max_threads()                   )
#else
    call rttov_k(                                           &
           errstat,                                         & !  --> error flag
           chanprof(1:nchansprofs),                         & ! <-- channels and profiles to calculate
           ropts,                                           & ! <-- options
           pack(profiles(iprof_s:iprof_e), mask=prof_used), & ! <--  profile array
           profiles_k(1:nchansprofs),                       & ! <--> profile increments;
           coefs(ic),                                       & ! <--  coefs array
           transmission,                                    & ! <--> transmittances
           transmission_k,                                  & ! <--> K matrix of transmittances
           radiance,                                        & ! <--> direct model output radiances
           radiance_k,                                      & ! <--> K matrix of radiances
           calcemis(1:nchansprofs),                         & ! <--  flag for internal emissivity calc
           emissivity(1:nchansprofs),                       & ! <-- input emissivities per channel
           emissivity_k(1:nchansprofs)                      )  ! <--> k matrix on input surface emissivity
#endif
    if (errstat == errorstatus_fatal) then
      print *, 'after rttov_k -> fatal error'
      status = ERR_RTTOV_CALL
      return
    endif

    ! fill results into output arrays
    FTRACE_BEGIN('rtifc_k:store')
    if (.not.present(istore)) then
      sind   = 0
      eind   = 0
      count1 = 0
      do iprof=1,nprof
        nusedchans = count (lprofs == iprof)
        sind       = eind + 1
        eind       = eind + nusedchans
        do ichan=1,count(lprofs==iprof)
          count1 = count1 + 1
          temp_k(:,ichan,iprof)=dble(profiles_k(count1)%t  (1+nlevs_top:))
          humi_k(:,ichan,iprof)=dble(profiles_k(count1)%q  (1+nlevs_top:))
          if (present(clw_k)) clw_k(:,ichan,iprof)=dble(profiles_k(count1)%clw(1+nlevs_top:))
          if (present( o3_k))  o3_k(:,ichan,iprof)=dble(profiles_k(count1)%o3 (1+nlevs_top:))
          if (present(co2_k)) co2_k(:,ichan,iprof)=dble(profiles_k(count1)%co2(1+nlevs_top:))
          if (present(n2o_k)) n2o_k(:,ichan,iprof)=dble(profiles_k(count1)%n2o(1+nlevs_top:))
          if (present( co_k))  co_k(:,ichan,iprof)=dble(profiles_k(count1)%co (1+nlevs_top:))
          if (present(ch4_k)) ch4_k(:,ichan,iprof)=dble(profiles_k(count1)%ch4(1+nlevs_top:))
          if(nlevs_top==1) then
            !--------------------------
            ! rttov9 compatibility mode
            !--------------------------
            temp_k(1,ichan,iprof) = temp_k(1,ichan,iprof) + dble(profiles_k(count1)%t(1))
            humi_k(1,ichan,iprof) = humi_k(1,ichan,iprof) + dble(profiles_k(count1)%q(1))
            if (present(clw_k)) clw_k(1,ichan,iprof) = clw_k(1,ichan,iprof) + dble(profiles_k(count1)%clw(1))
            if (present( o3_k))  o3_k(1,ichan,iprof) =  o3_k(1,ichan,iprof) + dble(profiles_k(count1)%o3 (1))
            if (present(co2_k)) co2_k(1,ichan,iprof) = co2_k(1,ichan,iprof) + dble(profiles_k(count1)%co2(1))
            if (present(n2o_k)) n2o_k(1,ichan,iprof) = n2o_k(1,ichan,iprof) + dble(profiles_k(count1)%n2o(1))
            if (present( co_k))  co_k(1,ichan,iprof) =  co_k(1,ichan,iprof) + dble(profiles_k(count1)%co (1))
            if (present(ch4_k)) ch4_k(1,ichan,iprof) = ch4_k(1,ichan,iprof) + dble(profiles_k(count1)%ch4(1))
          endif
        enddo
        t2m_k    (1:nusedchans,iprof) = dble(profiles_k(sind:eind)% s2m% t)
        q2m_k    (1:nusedchans,iprof) = dble(profiles_k(sind:eind)% s2m% q)
        stemp_k  (1:nusedchans,iprof) = dble(profiles_k(sind:eind)% skin% t)
        t_b      (1:nusedchans,iprof) = dble(radiance%bt(sind:eind)        )
        emissiv  (1:nusedchans,iprof) = dble(emissivity  (sind:eind)%emis_out)
        emissiv_k(1:nusedchans,iprof) = dble(emissivity_k(sind:eind)%emis_out)
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear)) t_b_clear(1:nusedchans,iprof)=dble(radiance%bt_clear(sind:eind))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad      )) rad      (1:nusedchans,iprof)=dble(radiance%total   (sind:eind))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear )) radclear (1:nusedchans,iprof)=dble(radiance%clear   (sind:eind))
        end if
        if (present(   radtotal)) radtotal (  1:nchans,iprof)=dble(radiance%clear   (  sind:eind) )
        if (l_radoverc) then
          radovercast (:,1:nusedchans,iprof)=dble(radiance%overcast(:,sind:eind) )
        endif
        if (l_transm) &
             transm(:,1:nchans,iprof) = dble(transmission%tau_levels(1+nlevs_top:,sind:eind))
#if defined(_DACE_) && !defined(__ICON__)
        if (l_opdep) &
             opdep(:,1:nchans,iprof) = dble(transmission%opdep_ref(1+nlevs_top:,sind:eind))
#endif
        if (present(    psurf_k))     psurf_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%s2m%p    )
        if (present(     u10m_k))      u10m_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%s2m%u    )
        if (present(     v10m_k))      v10m_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%s2m%v    )
        if (present(  o3_surf_k))   o3_surf_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%s2m%o    )
        if (present(    wfetc_k))     wfetc_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%s2m%wfetc)
        if (present(      ctp_k))       ctp_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%ctp      )
        if (present(cfraction_k)) cfraction_k (  1:nusedchans,iprof)=dble(profiles_k(sind:eind)%cfraction)
      enddo
    else
      do i = 1, nchansprofs
        ! Forward calc
        t_b(istore(i,1),istore(i,2)) = dble(radiance% bt(i))
        if (btest(rad_out,OUT_CSB)) then
          if (present(t_b_clear  )) t_b_clear  (istore(i,1),istore(i,2)) = dble(radiance% bt_clear(i))
        end if
        if (btest(rad_out,OUT_ASR)) then
          if (present(rad        )) rad        (istore(i,1),istore(i,2)) = dble(radiance% total   (i))
        end if
        if (btest(rad_out,OUT_CSR)) then
          if (present(radclear   )) radclear   (istore(i,1),istore(i,2)) = dble(radiance% clear   (i))
        end if
        if (present(radtotal   )) radtotal   (istore(i,1),istore(i,2))   = dble(radiance% clear     (i))
        if (l_radoverc) radovercast(:,istore(i,1),istore(i,2)) = dble(radiance% overcast(:,i))
        if (l_transm) &
             transm(:,istore(i,1),istore(i,2)) = dble(transmission%tau_levels(1+nlevs_top:,i))
#if defined(_DACE_) && !defined(__ICON__)
        if (l_opdep) &
             opdep(:,istore(i,1),istore(i,2)) = dble(transmission%opdep_ref(1+nlevs_top:,i))
#endif
        ! K (profile variables)
        temp_k(:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%t (1+nlevs_top:))
        humi_k(:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%q (1+nlevs_top:))
        if (present(clw_k )) clw_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%clw(1+nlevs_top:))
        if (present(o3_k  )) o3_k  (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%o3 (1+nlevs_top:))
        if (present(co2_k )) co2_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%co2(1+nlevs_top:))
        if (present(n2o_k )) n2o_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%n2o(1+nlevs_top:))
        if (present(co_k  )) co_k  (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%co (1+nlevs_top:))
        if (present(ch4_k )) ch4_k (:,istore(i,1),istore(i,2)) = dble(profiles_k(i)%ch4(1+nlevs_top:))
        if(nlevs_top==1) then
          !--------------------------
          ! rttov9 compatibility mode
          !--------------------------
          temp_k(1,istore(i,1),istore(i,2)) = temp_k(1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%t(1))
          humi_k(1,istore(i,1),istore(i,2)) = humi_k(1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%q(1))
          if (present(clw_k )) clw_k (1,istore(i,1),istore(i,2)) = clw_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%clw(1))
          if (present(o3_k  )) o3_k  (1,istore(i,1),istore(i,2)) = o3_k  (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%o3 (1))
          if (present(co2_k )) co2_k (1,istore(i,1),istore(i,2)) = co2_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%co2(1))
          if (present(n2o_k )) n2o_k (1,istore(i,1),istore(i,2)) = n2o_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%n2o(1))
          if (present(co_k  )) co_k  (1,istore(i,1),istore(i,2)) = co_k  (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%co (1))
          if (present(ch4_k )) ch4_k (1,istore(i,1),istore(i,2)) = ch4_k (1,istore(i,1),istore(i,2)) + dble(profiles_k(i)%ch4(1))
        endif
        ! K (surface variables)
        t2m_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% t    )
        q2m_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% q    )
        stemp_k  (istore(i,1),istore(i,2)) = dble(profiles_k(i)% skin% t    )
        emissiv  (istore(i,1),istore(i,2)) = dble(emissivity  (i)%emis_out  )
        emissiv_k(istore(i,1),istore(i,2)) = dble(emissivity_k(i)%emis_out  )
        if (present(psurf_k    )) psurf_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% p    )
        if (present(u10m_k     )) u10m_k     (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% u    )
        if (present(u10m_k     )) v10m_k     (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% v    )
        if (present(o3_surf_k  )) o3_surf_k  (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% o    )
        if (present(wfetc_k    )) wfetc_k    (istore(i,1),istore(i,2)) = dble(profiles_k(i)% s2m% wfetc)
        if (present(ctp_k      )) ctp_k      (istore(i,1),istore(i,2)) = dble(profiles_k(i)% ctp        )
        if (present(cfraction_k)) cfraction_k(istore(i,1),istore(i,2)) = dble(profiles_k(i)% cfraction  )
      end do
    end if
    FTRACE_END('rtifc_k:store')

    ! Check profiles
    if (chk_reg_lims /= 0 .and. nprof > 0) then
      call realloc_rttov_arrays(status, nprof_aux,profiles(1)%nlevels,0,ropts,profs=profiles_aux)
      if (status /= NO_ERROR) return
      ropts_ = ropts
      do iprof = iprof_s, iprof_e
        if (prof_used(iprof)) then
          call check_prof(iprof, coefs(ic), lims_flag, ropts_, lpio)
          if (chk_reg_lims > 1 .and. present(reg_lim)) then
            if (.not.present(istore)) then
              do i = 1, nprof
                if (lprofs((i-1)*nchans+1) ==iprof) then
                  reg_lim(:,:,i) = lims_flag(:,:)
                  exit
                end if
              end do
            else
              do i = 1, nchansprofs
                if (lprofs(i)==iprof) then
                  reg_lim(:,:,istore(i,2)) = lims_flag(:,:)
                  exit
                end if
              end do
            end if
          end if
        end if
      end do
      call dealloc_rttov_arrays(status, profs=profiles_aux)
      if (status /= NO_ERROR) return
    end if

    ! Check impact of RTTOV-smoothing
#if defined(_RTTOV_GOD)
    if (l_chk_god) then
      do i = 1, nchansprofs
        stat = 0
        if (any(lprofs(i) == ipr(1:npr))) then
          print*,'check_god_ifc prof',trim(profiles(lprofs(i))%id),chans(i),coefs(ic)%coef%ff_ori_chn(chans(i))
          print*,'check_god_ifc ntr',chans(i),any(coefs(ic)%coef%god(:,chans(i))%ntr > 0)
        end if
        if (any(coefs(ic)%coef%god(:,chans(i))%ntr > 0)) then
          if (transmission%l_opdep) then
            call check_god_infl(coefs(ic)%coef%god(:,chans(i)), stat, msg_, &
                                od_ref=transmission%opdep_ref(:,i), debug=any(lprofs(i)==ipr(1:npr)))
          else
            call check_god_infl(coefs(ic)%coef%god(:,chans(i)), stat, msg_, &
                                tau=transmission%tau_levels(:,i), debug=any(lprofs(i)==ipr(1:npr)))
          end if
          if (stat /= 0) then
            if (btest(chk_god, 0)) then
              p => profiles(lprofs(i))
              write(msg,'(2x,A,I5.5," god-smoothing too large in ",A," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p%id)//' chan ',coefs(ic)%coef%ff_ori_chn(chans(i)),&
                   trim(msg_),p%latitude,p%longitude
              WR
            end if
          end if
        end if
        if ((btest(chk_god,1) .or. btest(chk_god,2)) .and. all(shape(rflag) > 0)) then
          if (present(istore)) then
            if (any(lprofs(i) == ipr(1:npr))) then
              print*,'check_god_ifc stat',stat,istore(i,:)
            end if
            rflag(istore(i,1), istore(i,2)) = stat
          else
            ! TODO !!!
          end if
        end if
      end do
    elseif (present(rflag)) then
      rflag(:,:) = 0
    end if
#endif /* _RTTOV_GOD */

    ! Trailer
    if (l_dealloc) then
       call dealloc_rttov_arrays(status,profs=profiles_k,rads=radiance_k,&
                                 transm=transmission_k)
       if (status /= NO_ERROR) return
       call dealloc_rttov_arrays(status,rads=radiance,transm=transmission)
       if (status /= NO_ERROR) return
    end if
FTRACE_END('rtifc_k')

  end subroutine rtifc_k


#if defined(_RTTOV_ATLAS)
 subroutine rtifc_emis_retrieve(iopt, lprofs, chans, obs, spec, emis, stat, pe, ldeb)
   integer,             intent(in)  :: iopt      ! options index
   integer,             intent(in)  :: lprofs(:) ! list of profile indices
   integer,             intent(in)  :: chans(:)  ! list of channel indices
   real(wp),            intent(in)  :: obs(:)    ! observed brightness temperature
   real(wp),            intent(in)  :: spec(:,:) ! specularity (ignored here, for compat. with RT13)
   real(wp),            intent(out) :: emis(:)   ! computed emissivities
   integer,             intent(out) :: stat      ! error status
   integer,             intent(in), optional:: pe
   logical,             intent(in), optional:: ldeb
    !--------------------------------------------------
   ! Dynamic retrieve of emissivity (following Karbou)
   !--------------------------------------------------
   character(len=19),  parameter   :: proc   = 'rtifc_emis_retrieve'
   type(t_rtopts),     pointer     :: rto    => null()
   type(rttov_options),pointer     :: ropts  => null()
   integer                         :: ic
   real(jprb), allocatable :: t_b(:,:)
   real(jprb), allocatable :: emis_in(:,:)
   real(jprb), allocatable :: radclear(:,:)
   real(jprb), allocatable :: radupclear(:,:)
   real(jprb), allocatable :: raddnclear(:,:)
   real(jprb), allocatable :: gamma(:,:)
   real(jprb), allocatable :: obsrad(:)
   real(jprb), allocatable :: radup(:)
   real(jprb), allocatable :: rademi(:)
   integer                 :: k

   if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
   rto   => rt_opts(iopt)
   ropts => rto%opts
   ic    =  rto%icoeff

   ! Check inputs consistency
   if (coefs(ic)% coef% id_sensor /= 2) then
     stat = ERR_INVALID_INSTR
     return
   end if
   if (size(lprofs) /= size(chans) .or. size(lprofs) /= size(chans) .or. &
        size(lprofs) /= size(obs) .or. size(lprofs) /= size(emis) ) then
     stat = ERR_DIM
     return
   end if
   if (any(lprofs(:) /= 1 )) then
     stat = ERR_DIM
     return
   end if

   ! Allocate auxiliary arrays
   allocate(t_b       (size(lprofs),1))
   allocate(emis_in   (size(lprofs),1))
   emis_in = 0._wp
   allocate(radclear  (size(lprofs),1))
   allocate(radupclear(size(lprofs),1))
   allocate(raddnclear(size(lprofs),1))
   allocate(gamma     (size(lprofs),1))
   allocate(obsrad    (size(lprofs)))
   allocate(radup     (size(lprofs)))
   allocate(rademi    (size(lprofs)))

   ! Call rtifc_direct for first guesses needed
   call rtifc_direct (                  &
        iopt,                           & ! <--  rttov instrument index
        lprofs,                         & ! <--  list of profile indices
        chans,                          & ! <--  list of channel indices
        emis_in,                        & ! <--> emissivities -
        t_b,                            & !  --> calculated brightness
        stat,                           & !  --> exit status
        radclear     = radclear,        & !  --> TOA radiance
        radupclear   = radupclear,      & !  --> TOA upweeling radiance(with emissivity term)
        raddnclear   = raddnclear,      & !  --> downelling radiance at surface
        transmtotal  = gamma )!,        & !  --> calculated transmittance
!        pe           = pe,             & ! <--  pe, debug
!        iprint       = iprint )          ! <--  debug
   if (stat /= NO_ERROR) return

   ! Compute emissivity
   do k = 1, size(chans)
     ! Compute radiance associated to observed brightness temperatures
     call planck(coefs(ic)%coef% planck1(chans(k)),coefs(ic)%coef% planck2(chans(k)), &
          obs(k),obsrad(k))
     ! Compute radiance associated to skin temperature
     call planck(coefs(ic)%coef% planck1(chans(k)),coefs(ic)%coef% planck2(chans(k)), &
          profiles(1)% skin% t, rademi(k))
     ! Compute upwelling radiance without emission term
     radup(k) = radupclear(k,1)-rademi(K)*emis_in(k,1)*gamma(k,1)
     ! Compute emissivity
     emis(k) = (obsrad(k)-radup(k)-raddnclear(k,1)*gamma(k,1))/(gamma(k,1)*(rademi(k)-raddnclear(k,1)))
   end do

   ! Clean up
   deallocate(t_b       )
   deallocate(emis_in   )
   deallocate(radclear  )
   deallocate(radupclear)
   deallocate(raddnclear)
   deallocate(gamma     )
   deallocate(obsrad    )
   deallocate(radup     )
   deallocate(rademi    )

 end subroutine rtifc_emis_retrieve

 subroutine rtifc_emis_sea(iopt, lprofs, chans, emis, stat, pe, ldeb)
   integer,             intent(in)  :: iopt      ! options index
   integer,             intent(in)  :: lprofs(:) ! list of profile indices
   integer,             intent(in)  :: chans(:)  ! list of channel indices
   real(wp),            intent(out) :: emis(:)   ! computed emissivities
   integer,             intent(out) :: stat      ! error status
   integer,             intent(in), optional:: pe
   logical,             intent(in), optional:: ldeb

   stat = ERR_NO_RTTOV_LIB

 end subroutine rtifc_emis_sea

 subroutine rtifc_init_atlas(iopts, atlas_id, angcorr, month, path, &
                             my_proc_id, n_proc, io_proc_id, mpi_comm_type)
   integer,            intent(in) :: iopts(:)        ! option indices
   integer,            intent(in) :: atlas_id(:)    ! atlas IDs
   logical,            intent(in) :: angcorr(:)     ! atlas IDs
   integer,            intent(in) :: month          ! Month number
   character(len=128), intent(in) :: path           ! path to atlases
   integer,            intent(in) :: my_proc_id     ! ID of local processors

   integer,            intent(in) :: n_proc         ! no. of processors in mpi communication domain
   integer,            intent(in) :: io_proc_id     ! ID of IO processor
   integer,            intent(in) :: mpi_comm_type  ! mpi communicator type
   !--------------------------------------------
   ! Loads MW emissivity atlas
   !--------------------------------------------
   character(len=16),       parameter   :: proc   = 'rtifc_init_atlas'
   integer,                 parameter   :: mcnrm  = 8
   type(t_rtopts),          pointer     :: rto
   character(len=300)                   :: msg
   integer,                 allocatable :: iatl_ir(:)
   integer,                 allocatable :: iatl_mw(:)
   type(rttov_emis_atlas_data), pointer :: atlas_(:)
   integer                              :: instr_cnrm(mcnrm), iatl_cnrm(mcnrm)
   integer                              :: ierr
   integer                              :: n_opt
   integer                              :: n_cnrm
   logical                              :: l_distrib
   integer                              :: ic, iopt
   integer                              :: i, j

   pe_ifc = my_proc_id
   pe_rt  = pe_ifc

   n_atlas = 0
   ierr    = 0

   l_distrib = .false.
#if defined(_RTIFC_DISTRIBCOEF)
   l_distrib = (n_proc > 1) .and. read1pe
#endif

   if (associated(atlas)) call finish(proc, 'Atlases initialized already.')

   ! Remove double entries
   n_opt = size(iopts)
   if (n_opt == 0) return

   ! do i = 1, n_opt
   !   write(0,*) 'init_atlas',i,iopts(i),atlas_id(i),rt_opts(iopts(i))%satid,rt_opts(iopts(i))%instr
   ! end do

   i = maxval(atlas_id)
   allocate(iatl_ir(i), iatl_mw(i))
   iatl_ir = 0 ; iatl_mw = 0

   if (io_proc_id == pe_ifc .or. .not.l_distrib) then
     allocate(atlas(n_opt))
     n_atlas = 0
     n_cnrm  = 0
     write(*,*)
     write(*,'(1x,A)') 'Initialize RTTOV emissivity atlases:'
     do i = 1, n_opt
       iopt = iopts(i)
       if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
       rto => rt_opts(iopt)
       ic = rto%icoeff
       ierr = 0
       msg  = ''
       select case(coefs(ic)%coef%id_sensor)
       case(sensor_id_ir,sensor_id_hi) !IR
         if (atlas_single_inst) then
           n_atlas = n_atlas + 1
           call rttov_setup_emis_atlas(ierr, rto%opts, month, atlas_type_ir, atlas(n_atlas), &
                atlas_id=atlas_id(i), path = path, coefs=coefs(ic), ir_atlas_ang_corr=angcorr(i))
           call add_iatl(n_atlas)
           msg = c_atlas(atlas_type_ir, atlas_id(i), coefs(ic)%coef%id_inst, rto%satid)
         else
           if (iatl_ir(atlas_id(i)) <= 0) then
             n_atlas = n_atlas + 1
             call rttov_setup_emis_atlas(ierr, rto%opts, month, atlas_type_ir, atlas(n_atlas), &
                  atlas_id=atlas_id(i), path = path, ir_atlas_ang_corr=angcorr(i))
             iatl_ir(atlas_id(i)) = n_atlas
             msg = c_atlas(atlas_type_ir, atlas_id(i))
           end if
           call add_iatl(iatl_ir(atlas_id(i)))
         end if
       case(sensor_id_mw,sensor_id_po) !MW
         if (atlas_id(i) == cnrm_mw_atlas_id) then
           if (any(instr_cnrm(1:n_cnrm) == coefs(ic)%coef%id_inst)) then
             do j = 1, n_cnrm
               if (instr_cnrm(j) == coefs(ic)%coef%id_inst) then
                 call add_iatl(iatl_cnrm(j))
                 exit
               endif
             end do
           else
             n_atlas = n_atlas + 1
             call rttov_setup_emis_atlas(ierr, rto%opts, month, atlas_type_mw, atlas(n_atlas), &
                  atlas_id=atlas_id(i), path = path, coefs=coefs(ic))
             n_cnrm = n_cnrm + 1
             if (n_cnrm > mcnrm) call finish(proc,'mcnrm too small')
             instr_cnrm(n_cnrm) = coefs(ic)%coef%id_inst
             iatl_cnrm (n_cnrm) = n_atlas
             call add_iatl(iatl_cnrm(n_cnrm))
             msg = c_atlas(atlas_type_mw, atlas_id(i), coefs(ic)%coef%id_inst)
           end if
         else
           if (iatl_mw(atlas_id(i)) <= 0) then
             n_atlas = n_atlas + 1
             call rttov_setup_emis_atlas(ierr, rto%opts, month, atlas_type_mw, atlas(n_atlas), &
                  atlas_id=atlas_id(i), path = path)
             iatl_mw(atlas_id(i)) = n_atlas
             msg = c_atlas(atlas_type_mw, atlas_id(i))
           end if
           call add_iatl(iatl_mw(atlas_id(i)))
         end if
       case default
         call finish(proc,'sensor_id currently not implemented')
       end select
       if (ierr == 0) then
         if (io_proc_id == pe_ifc .and. msg /= '') &
              write(stdout,'(3x,"atlas ",I2,": ",A)') n_atlas, trim(msg)
       else
         call finish(proc, trim(msg)//' emissivity atlas initialization FAILED!')
       end if
     end do
     ! Reduce atlas to correct size
     allocate(atlas_(n_atlas))
     atlas_(1:n_atlas) = atlas(1:n_atlas)
     deallocate(atlas)
     atlas => atlas_

   endif

#if defined(_RTIFC_DISTRIBCOEF)
   ! distribute atlases to all PEs
   if (l_distrib) then
     call p_bcast(n_atlas,io_proc_id,mpi_comm_type)
     if (io_proc_id == pe_ifc) then
       write(stdout,'(3x,A,I2,A)') 'Distribute ',n_atlas,' emissivity atlases.'
     else
       allocate(atlas(n_atlas))
     end if
     do i = 1, n_atlas
       if (io_proc_id /= pe_ifc) call rttov_deallocate_emis_atlas(atlas(i))
       call p_bcast(atlas(i),io_proc_id,mpi_comm_type)
     end do
     do i = 1, n_opt
       call p_bcast(rt_opts(iopts(i))%natl,io_proc_id,mpi_comm_type)
       call p_bcast(rt_opts(iopts(i))%iatl,io_proc_id,mpi_comm_type)
     end do
   end if
#endif

 contains

   subroutine add_iatl(iatl)
     integer, intent(in) :: iatl
     if (.not.any(rto%iatl(1:rto%natl) == iatl)) then
       rto%natl = rto%natl + 1
       if (rto%natl > m_atl) call finish('add_iatl@'//proc, 'rto%natl > m_atl')
       rto%iatl(rto%natl) = iatl
     end if
   end subroutine add_iatl

 end subroutine rtifc_init_atlas

 subroutine rtifc_emis_atlas(iopt, atlas_id, lprofs, chan, emis, stat, ldebug, max_dst)
   integer,  intent(in)           :: iopt     ! options index
   integer,  intent(in)           :: atlas_id ! 1: TELSEM, 2: CNRM, 3: CAMEL-Climatology
   integer,  intent(in)           :: lprofs(:)! list of profile indices
   integer,  intent(in)           :: chan(:)  ! list of channels
   real(wp), intent(out)          :: emis(:)  ! emissivities
   integer,  intent(out)          :: stat     ! error status
   logical,  intent(in), optional :: ldebug
   real(wp), intent(in), optional :: max_dst  ! ignored here, for compatibility with RT13
   !--------------------------
   ! Get emissivity from atlas
   !--------------------------
   character(len=16),  parameter :: proc   = 'rtifc_emis_atlas'
   type(t_rtopts),     pointer   :: rto    => null()
   type(rttov_options),pointer   :: ropts  => null()
   type(rttov_chanprof)          :: chanprof(size(chan))
   integer                       :: ic, iatl
   integer                       :: i,j
   logical                       :: ldeb

   ldeb = .false.
   if (present(ldebug)) ldeb = ldebug
   if (ldeb) then
     ipr_deb = lprofs(1)
   else
     ipr_deb = -1
   end if

   if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
   rto   => rt_opts(iopt)
   ropts => rto%opts
   ic    =  rto%icoeff
   iatl = 0
   do i = 1, rto%natl
     j = rto%iatl(i)
     if (j <= 0 .or. j > n_atlas) call finish(proc, 'invalid rto%iatl')
     if (atlas(j)% atlas_id == atlas_id) then
       iatl = j
       exit
     end if
   end do
   if (iatl <= 0) then
     write(0,*) 'sat,instr',rto%satid,rto%instr
     write(0,*) 'atlas',atlas_id
     call finish(proc, 'atlas not intialized')
   end if

   if (size(lprofs) /= size(chan) .or. &
        size(lprofs) /= size(emis) ) then
     stat = ERR_DIM
     return
   end if
   if (any(lprofs(:) /= 1 )) then
     ! Can only deal with one profile at a time
     stat = ERR_DIM
     return
   end if
   do i = 1,size(chan)
      chanprof(i)% chan = chan(i)
      chanprof(i)% prof = lprofs(i)
   enddo
   if (ldeb) call rttov_print_profile(profiles(lprofs(1)), stdout, trim(proc))
   call rttov_get_emis(stat, ropts, chanprof, profiles, coefs(ic), atlas(iatl), emis)
   if (stat /= 0) then
     write(0,*) 'stat=',stat
     call finish(proc, 'rttov_get_emis failed')
   end if

 end subroutine rtifc_emis_atlas


 function c_atlas(type, id, instr, satid) result(c)
   character(len=30)             :: c
   integer, intent(in)           :: type
   integer, intent(in)           :: id
   integer, intent(in), optional :: instr
   integer, intent(in), optional :: satid
   character(len=30) :: ci

   select case(type)
   case(atlas_type_ir)
     select case(id)
     case(uwiremis_atlas_id)
       c = 'UWIR'
     case(camel_atlas_id)
       c = 'CAMEL07'
     case default
       c = '??????'
     end select
   case(atlas_type_mw)
     select case(id)
     case(telsem2_atlas_id)
       c = 'TELSEM2'
     case(cnrm_mw_atlas_id)
       c = 'CNRM'
     case default
       c = '??????'
     end select
   case default
     c = '??????'
   end select
   if (present(instr)) then
     ci = ''
     if (present(satid)) then
       write(ci,'("(instr=",I3,",satid=",I3,")")') instr, satid
     else
       write(ci,'("(instr=",I3,")")') instr
     end if
   else
     ci = '(all instrs.)'
   end if
   c = trim(c)//trim(ci)

 end function c_atlas


 subroutine rtifc_init_brdf_atlas(month, path, my_proc_id, n_proc, io_proc_id, mpi_comm_type, stat)
   integer,            intent(in) :: month          ! Month number
   character(len=128), intent(in) :: path           ! path to atlases
   integer,            intent(in) :: my_proc_id     ! ID of local processors

   integer,            intent(in) :: n_proc         ! no. of processors in mpi communication domain
   integer,            intent(in) :: io_proc_id     ! ID of IO processor
   integer,            intent(in) :: mpi_comm_type  ! mpi communicator type
   integer,            intent(out):: stat           ! error status

   stat = ERR_NO_ATLAS
   return

 end subroutine rtifc_init_brdf_atlas

 subroutine rtifc_brdf_atlas(iopt, profs, chans, refl, stat,refl_flag)
   integer,            intent(in)  :: iopt        ! options index
   integer,            intent(in)  :: profs(:)    ! list of profile indices
   integer,            intent(in)  :: chans(:)    ! list of channels
   real(wp),           intent(out) :: refl(:)     ! emissivities
   integer,            intent(out) :: stat        ! error status
   integer, optional,  intent(out) :: refl_flag(:)! emissivities flags

   stat = ERR_NO_ATLAS
   return

 end subroutine rtifc_brdf_atlas

 subroutine rtifc_tskin_retrieve(iopt, lprofs, channum, chans, obs, &
      spec, emis, tskin, stat, tsfl, pe, ldeb, spt_hd_id)
   integer,  intent(in)            :: iopt      ! options index
   integer,  intent(in)            :: lprofs(:) ! list of profile indices
   integer,  intent(in)            :: channum(:)! list of channel numbers
   integer,  intent(in)            :: chans(:)  ! list of channel indices
   real(wp), intent(in)            :: obs(:)    ! observed brightness temperature
   real(wp), intent(in)            :: spec(:,:) ! specularity
   real(wp), intent(inout)         :: emis(:)   ! atlas emissivities
   real(wp), intent(out)           :: tskin(:)  ! computed skin temperature
   integer,  intent(out)           :: stat      ! error status
   logical,  intent(inout)         :: tsfl      ! set tskin flag
   integer,  intent(in),   optional:: pe
   logical,  intent(in),   optional:: ldeb
   integer,  intent(in),   optional:: spt_hd_id

   stat = ERR_NO_RTTOV_LIB

 end subroutine rtifc_tskin_retrieve

#endif /* _RTTOV_ATLAS */

!------------------------------------------------------------------------------

  subroutine rtifc_coef_prop(iopt, version, sign_opdep, preslev, nlevs, igas, gas_bkg, &
                             satid, grid, instr)
    integer,  intent(in)            :: iopt     ! options index
    integer,  intent(out), optional :: version
    integer,  intent(out), optional :: sign_opdep
    real(wp), intent(out), optional :: preslev(:)
    integer,  intent(out), optional :: nlevs
    integer,  intent(in),  optional :: igas
    real(wp), intent(out), optional :: gas_bkg(:)
    integer,  intent(out), optional :: satid
    integer,  intent(out), optional :: grid
    integer,  intent(out), optional :: instr
    !----------------------
    ! Get coef properties
    !----------------------
    character(len=15),  parameter :: proc   = 'rtifc_coef_prop'
    type(t_rtopts),     pointer   :: rto    => null()
    type(rttov_options),pointer   :: ropts  => null()
    integer                       :: ic
    integer                       :: nl, nlevs_user
    integer                       :: vers

    if (present(version   )) version    = -1
    if (present(sign_opdep)) sign_opdep =  0

    if (iopt<=0 .or. iopt>n_opts) call finish(proc, 'invalid option index')
    rto   => rt_opts(iopt)
    ropts => rto%opts
    ic    =  rto%icoeff

    vers = coefs(ic)%coef%fmv_model_ver

    if (present(version)) version = vers

    if (present(sign_opdep)) then
      select case(vers)
      case(7,8)
        sign_opdep = -1
      case(13)
        sign_opdep =  1
      case default
        call finish(proc, 'unknown version')
      end select
    end if

    if (present(satid)) satid = rto%satid
    if (present(instr)) instr = rto%instr
    if (present(grid )) grid  = rto%grid

    nl = coefs(ic)%coef%nlevels
    if (present(nlevs)) nlevs = nl
    if (present(preslev)) then
      nlevs_user = min(size(preslev), nl)
      if (nlevs_user < nl - 1) &
           call finish(proc, 'invalid size of preslev array')
      nlevs_top = nl - nlevs_user
      preslev(1:nlevs_user) = coefs(ic)%coef%ref_prfl_p(1+nlevs_top:)
    end if

    if (present(gas_bkg)) then
      if (.not.present(igas)) call finish(proc, 'gas_bkg requires option "igas"')
      if (igas < 1 .or. igas > coefs(ic)%coef%fmv_gas) &
           call finish(proc, 'invalid value for igas/gas not supported')
      nlevs_user = min(size(gas_bkg), nl)
      if (nlevs_user < nl - 1) &
           call finish(proc, 'invalid size of preslev array')
      nlevs_top = nl - nlevs_user
      gas_bkg(1:nlevs_user) = coefs(ic)%coef%bkg_prfl_mr(1+nlevs_top:,coefs(ic)%coef%fmv_gas_pos(igas))
    end if

  end subroutine rtifc_coef_prop



  subroutine check_prof(iprof, coefs, flg, ropts, lpio)
    integer,             intent(in)    :: iprof
    type(rttov_coefs),   intent(in)    :: coefs
    type(rttov_options), intent(inout) :: ropts
    integer,             intent(out)   :: flg(:,:)
    logical,             intent(in)    :: lpio

    character(len=256)           :: msg = ''
    type(rttov_profile), pointer :: p1, p2
    integer                      :: ilev
    logical                      :: apply_reg_lims_aux

    flg = 0
    if (chk_reg_lims /= 0) then
      call rttov_copy_prof(profiles_aux(1:1), profiles(iprof:iprof), larray=.true., lscalar=.true.)
      apply_reg_lims_aux            = ropts%config%apply_reg_limits
      ropts%config%apply_reg_limits = .true.
      call rttov_convert_profile_units(ropts, coefs, profiles(iprof:iprof), profiles_aux(1:1))
      call rttov_copy_prof(profiles_aux(2:2), profiles_aux(1:1), larray=.true., lscalar=.true.)
      call rttov_apply_reg_limits(ropts, profiles_aux(2:2), profiles_aux(1:1), coefs%coef, coefs%coef_pccomp)
      p1 => profiles_aux(2)
      p2 => profiles_aux(1)
      ropts%config%apply_reg_limits = apply_reg_lims_aux
      ! T
      do ilev = 1+nlevs_top, p1%nlevels
        flg(ilev-nlevs_top,1) = 0
        if (p1%p(ilev) < chk_plim_t) then
          if (p1%t(ilev) < p2%t(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",f7.3,"<",f7.3," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' t exceeds lower bound in level ',ilev,p1%t(ilev),p2%t(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,1) = 1
          elseif (p1%t(ilev) > p2%t(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",f7.3,">",f7.3," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' t exceeds upper bound in level ',ilev,p1%t(ilev),p2%t(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,1) = -1
          end if
        end if
      end do
      ! q
      do ilev = 1+nlevs_top, p1%nlevels
        flg(ilev-nlevs_top,2) = 0
        if (p1%p(ilev) < chk_plim_q) then
          if (p1%q(ilev) < p2%q(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",e13.6,"<",e13.6," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' q exceeds lower bound in level ',ilev,p1%q(ilev),p2%q(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,2) = 1
          elseif (p1%q(ilev) > p2%q(ilev)) then
            if (btest(chk_reg_lims,0)) then
              write(msg,'(2x,A,I3," : ",e13.6,">",e13.6," (",f7.3,"N",f8.3,"E)")') &
                   'profile '//trim(p1%id)//' q exceeds upper bound in level ',ilev,p1%q(ilev),p2%q(ilev),p1%latitude,p1%longitude
              WR
            end if
            flg(ilev-nlevs_top,2) = -1
          end if
        end if
      end do
      ! o3
!       assoc1 = associated(p1%o3)
!       assoc2 = associated(p2%o3)
!       if (.not.(assoc1 .eqv. assoc2)) flg = 1
!       if (assoc1 .and. assoc2) then
!         do ilev = 1, p1%nlevels
!           if (mask_lims_o3(ilev) .and. p1%o3(ilev) /= p2%o3(ilev)) then
!             print*,'profile '//trim(p1%id)//' o3 out of bounds in level ',ilev
!             flg = 1
!           end if
!         end do
!       end if
      ! TODO: other variables ...
    end if

  end subroutine check_prof


  subroutine rttov_clear_prof(profiles)
    type(rttov_profile),intent(inout) :: profiles(:) ! initialized profiles
    !------------------------------------------
    ! initialize rttov rttov_profile structures
    !------------------------------------------
    integer :: i

    ! Initialize scalar components
    do i=1,size(profiles)
       profiles(i)% skin% surftype  = -1       ! no meaning
       profiles(i)% skin% watertype = -1       ! no meaning
       profiles(i)% skin% t         =  0._jprb ! on temperature
       profiles(i)% skin% fastem    =  0._jprb ! Fastem
       profiles(i)% skin% salinity  =  0._jprb ! ?
       profiles(i)% skin% soil_moisture =  0._jprb ! ?
       profiles(i)% skin% snow_fraction =  0._jprb ! ?
       profiles(i)% skin% foam_fraction =  0._jprb ! ?
       profiles(i)% s2m% t          =  0._jprb ! temperature
       profiles(i)% s2m% q          =  0._jprb ! WV
       profiles(i)% s2m% o          =  0._jprb ! O3
       profiles(i)% s2m% p          =  0._jprb ! pressure
       profiles(i)% s2m% u          =  0._jprb ! wind components
       profiles(i)% s2m% v          =  0._jprb ! wind components
       profiles(i)% s2m% wfetc      =  0._jprb ! wind fetc
       profiles(i)% zenangle        =  0._jprb ! no meaning
       profiles(i)% azangle         =  0._jprb ! no meaning
       profiles(i)% sunzenangle     =  0._jprb ! no meaning
       profiles(i)% sunazangle      =  0._jprb ! no meaning
       profiles(i)% elevation       =  0._jprb ! no meaning
       profiles(i)% latitude        =  0._jprb ! no meaning
       profiles(i)% ctp             =  0._jprb ! cloud top pressure
       profiles(i)% cfraction       =  0._jprb ! cloud fraction
       profiles(i)% id              =  ''
       profiles(i)% date            =  -1
       profiles(i)% time            =  -1
       profiles(i)% gas_units       = default_gas_units
       profiles(i)% mmr_cldaer      = .false.
       profiles(i)% idg             = default_idg
       profiles(i)% ice_scheme      = default_ice_scheme
       profiles(i)% longitude       =  0._jprb ! no meaning
       profiles(i)% Be              =  0._jprb ! no meaning
       profiles(i)% cosbk           =  0._jprb ! no meaning
    enddo

    ! Initialize profile components
    do i=1,size(profiles)
!NEC$ shortloop
!CDIR ARRAYCOMB
       profiles(i)% p  (:)          =  0._jprb ! no AD on pressure levels
       profiles(i)% t  (:)          =  0._jprb ! temperature
       profiles(i)% q  (:)          =  0._jprb ! water vapour (ppmv)
       if (associated(profiles(i)%o3))  &
         profiles(i)% o3 (:)        =  0._jprb ! O3 ppmv
       if (associated(profiles(i)%co2)) &
         profiles(i)% co2(:)        =  0._jprb ! CO2 ppmv
       if (associated(profiles(i)%n2o)) &
         profiles(i)% n2o(:)        =  0._jprb ! N2O ppmv
       if (associated(profiles(i)%co))  &
         profiles(i)% co (:)         =  0._jprb ! CO ppmv
       if (associated(profiles(i)%ch4)) &
         profiles(i)% ch4(:)         =  0._jprb ! CH4 ppmv
       if (associated(profiles(i)%clw)) &
         profiles(i)% clw(:)         =  0._jprb ! cloud liquid water (kg/kg)
!CDIR END ARRAYCOMB
    end do

    do i=1,size(profiles)
      if (associated(profiles(i)% aerosols)) &
           profiles(i)% aerosols = 0._jprb ! (iaer,   nlevels)
      if (associated(profiles(i)% cloud)) &
          profiles(i)% cloud     =  0._jprb ! (ncldtyp,nlevels)
      if (associated(profiles(i)% cfrac)) &
          profiles(i)% cfrac     =  0._jprb ! (icld,   nlevels)
    end do

  end subroutine rttov_clear_prof


  subroutine rttov_clear_rad_var(nchannels,rad)
    integer(jpim),       intent(in)    :: nchannels
    type(rttov_radiance),intent(inout) :: rad
    !-----------------------------------------------------------------------
    ! subroutine clears the per event variable part
    ! of an rttov rttov_radiance structure.
    !-----------------------------------------------------------------------
    rad% clear     (1:nchannels)   = 0.0_jprb
    rad% cloudy    (1:nchannels)   = 0.0_jprb
    rad% total     (1:nchannels)   = 0.0_jprb
    rad% bt        (1:nchannels)   = 0.0_jprb
    rad% bt_clear  (1:nchannels)   = 0.0_jprb
    rad% overcast  (:,1:nchannels) = 0.0_jprb
    rad% refl_clear(1:nchannels)   = 0.0_jprb
    rad% refl      (1:nchannels)   = 0.0_jprb
  end subroutine rttov_clear_rad_var


  subroutine rttov_clear_rad2_var(nchannels,rad)
  integer(jpim),intent(in)           :: nchannels
  type(rttov_radiance2),intent(inout) :: rad
  !-----------------------------------------------------------------------
  ! subroutine clears the per event variable part
  ! of an rttov rttov_radiance structure.
  !-----------------------------------------------------------------------
    rad% upclear    (1:nchannels)   = 0.0_jprb
    rad% dnclear    (1:nchannels)   = 0.0_jprb
    rad% refldnclear(1:nchannels)   = 0.0_jprb
    rad% up         (:,1:nchannels) = 0.0_jprb
    rad% down       (:,1:nchannels) = 0.0_jprb
    rad% surf       (:,1:nchannels) = 0.0_jprb
  end subroutine rttov_clear_rad2_var


  subroutine alloc_rttov_arrays(status, n_profs,nlevs,n_channels,l_init,ropts,&
                                profs,rads,rads2,transm,n_levels_coef,height)
    integer,                  intent(out)            :: status
    integer(jpim),            intent(in)             :: n_profs        !Number of profiles to allocate
    integer(jpim),            intent(in)             :: nlevs
    integer(jpim),            intent(in)             :: n_channels     !Number of channels &
                                                                       !(e.g channels/prof * profiles)
    logical(jplm),            intent(in)             :: l_init
    type(rttov_options),      intent(in)             :: ropts
    type(rttov_profile),      pointer,      optional :: profs(:)
    type(rttov_radiance),     intent(inout),optional :: rads
    type(rttov_radiance2),    intent(inout),optional :: rads2
    type(rttov_transmission), intent(inout),optional :: transm
    integer(jpim),            intent(in),   optional :: n_levels_coef
    real(kind=jprb),          pointer,      optional :: height(:)
    !-------------------------------------------------------------
    ! subroutine allocates the RTTOV profile and radiance structures
    !-------------------------------------------------------------
    integer      :: stat

FTRACE_BEGIN('alloc_rttov_arrays')

    nlay    = nlevs - 1

    status  = NO_ERROR

    if (present(profs )) then
      call alloc_profs(profs )
      if (status /= NO_ERROR) return
    end if

    if (present(rads)) then
      ! allocate radiance structure
      call rttov_alloc_rad(stat, n_channels, rads, nlevs, rttov_alloc, &
                           radiance2=rads2, init=.true.)
      if(stat /= 0) then
        status = ERR_ALLOC
        return
      endif
    end if

    if (present(transm)) then
      ! allocate parts of transmission structure and error status
      call rttov_alloc_transmission(stat, transm, nlevs, n_channels, rttov_alloc, init=.true. &
#if defined(_DACE_) && !defined(__ICON__)
                                    , nlevels_coef=n_levels_coef &
#endif
                                    )
      if(stat /= 0) then
        status = ERR_ALLOC
        return
      end if
    end if

    if (present(height)) then
      allocate(height(n_channels), stat=stat)
      if(stat /= 0) then
        status = ERR_ALLOC
        return
      end if
    end if

FTRACE_END('alloc_rttov_arrays')

  contains

    subroutine alloc_profs(profs)
      type (rttov_profile),     pointer :: profs(:)
      integer(jpim) :: loc1 (1), loc2 (1)
      logical       :: l_req        ! Whether allocation of profiles is required

      ! allocate profile structure
      l_req = .not.associated(profs)
      if (.not.l_req) l_req = (size(profs) < n_profs)
      if (.not.l_req) l_req = (size(profs(1)%t) < nlevs)

      if (l_req) then
        if (associated(profs)) then
          call rttov_alloc_prof(stat,size(profs),profs,         &
               size(profs(1)%p),ropts,rttov_dealloc,coefs(1))
          if (stat/= 0) then ; status = ERR_ALLOC ; return ; end if
          deallocate(profs, stat = stat)
          if (stat/= 0) then ; status = ERR_ALLOC ; return ; end if
        END if
        allocate(profs(n_profs), stat = stat)
        if (stat/= 0) then ; status = ERR_ALLOC ; return ; end if

        loc1 = maxloc(coefs(:)%coef_scatt_ir%fmv_aer_comp)
        loc2 = maxloc(coefs(:)%coef_scatt_ir%fmv_wcl_comp)

        if (loc1(1) /= loc2(1)) then
          status = ERR_CLOUD_AERO_MISSM
          return
        endif

        call rttov_alloc_prof(stat, n_profs, profs, nlevs, ropts, rttov_alloc, &
                              coefs(loc1(1)), init      = l_init)
        if(stat /= 0) then
          status  = ERR_ALLOC
          return
        endif
      else
        call rttov_init_prof(profs(1:n_profs))
      end if

      profs(:)% gas_units  = default_gas_units
      profs(:)% mmr_cldaer = .false.
    end subroutine alloc_profs


  end subroutine alloc_rttov_arrays

  subroutine dealloc_rttov_arrays(status, profs,rads,rads2,transm,height)
    integer,                  intent(out)             :: status
    type(rttov_profile),      pointer,       optional :: profs(:)
    type(rttov_radiance),     intent(inout), optional :: rads
    type(rttov_radiance2),    intent(inout), optional :: rads2
    type(rttov_transmission), intent(inout), optional :: transm
    real(kind=jprb),          pointer,       optional :: height(:)
    !---------------------------------------------------------------
    ! subroutine deallocates the RTTOV profile and radiance structures
    !---------------------------------------------------------------
    integer                    :: status_rt
    logical                    :: lr2

FTRACE_BEGIN('dealloc_rttov_arrays')

    status = NO_ERROR

    if (idef0 <= 0) then
      status = ERR_NO_OPTS_TMPL
      return
    end if

    if (present(rads)) then
      if (associated(rads%clear)) then
        if (present(rads2)) then
          lr2 = associated(rads2%up)
        else
          lr2 = .false.
        end if
        if (lr2) then
          call rttov_alloc_rad(status_rt,size(rads%clear),             &
                               rads,size(rads% overcast, 1),rttov_dealloc, radiance2=rads2)
        else
          call rttov_alloc_rad(status_rt,size(rads%clear),             &
                               rads,size(rads% overcast, 1),rttov_dealloc)
        end if
        if(status_rt /= errorstatus_success) then
          status = ERR_ALLOC
          return
        endif
      end if
    endif
    if (present(profs)) then
      if (associated(profs)) then
        call rttov_alloc_prof(status_rt,size(profs),profs,         &
             size(profs(1)%p),rt_opts_tmpl(idef0)%opts,rttov_dealloc,coefs(1))
        if(status_rt /= errorstatus_success) then
          status = ERR_ALLOC
          return
        endif
        deallocate(profs)
        nullify(profs)
      endif
    end if

    if (present(transm)) then
      if (associated(transm%tau_levels)) then
        call rttov_alloc_transmission(status_rt,transm, &
                                      size(transm%tau_levels,1) - 1, &
                                      size(transm%tau_levels,2),     &
                                      rttov_dealloc)
        if(status_rt /= errorstatus_success) then
          status = ERR_ALLOC
          return
        endif
      end if

      if (present(height)) then
        deallocate(height, stat=status_rt)
        if (status_rt /= 0) then
          status = ERR_ALLOC
          return
        end if
      end if

    end if

FTRACE_END('dealloc_rttov_arrays')

  end subroutine dealloc_rttov_arrays


  subroutine realloc_rttov_arrays(status, n_profiles, n_levels, n_channels, &
       ropts, profs,rads,rads2,transm,n_levels_coef,height)
    integer,                  intent(out)             :: status
    integer(jpim),            intent(in)              :: n_profiles     ! requested number of profiles
    integer(jpim),            intent(in)              :: n_levels       ! requested levels
    integer(jpim),            intent(in)              :: n_channels     ! requested channels
    type(rttov_options),      intent(in)              :: ropts
    type(rttov_profile),      pointer,       optional :: profs(:)       ! profile array to realloc
    type(rttov_radiance),     intent(inout), optional :: rads           ! radiance array to realloc
    type(rttov_radiance2),    intent(inout), optional :: rads2          ! radiance array to realloc
    type(rttov_transmission), intent(inout), optional :: transm         ! transmission array to realloc
    integer(jpim),            intent(in),    optional :: n_levels_coef
    real(kind=jprb),          pointer,       optional :: height(:)
    !----------------------------------------------------------------------------------
    ! subroutine reallocates the profile, radiance and transmission structure if necessary
    !----------------------------------------------------------------------------------
    integer :: n_p
    logical :: l_req, l_req2

    status = NO_ERROR

    if (present(profs)) then
      l_req = .not.associated(profs)
      if (.not.l_req) then
        n_p = size(profs)
        l_req = (n_profiles > n_p)
        if (.not.l_req .and. n_p > 0) then
          l_req = n_levels /= profs(1)%nlevels
          ! RTTOV does not allocate gas profiles by default
          if (.not.l_req) l_req = ropts%rt_ir%ozone_data .neqv. associated(profs(1)%o3)
          if (.not.l_req) l_req = ropts%rt_ir%co2_data   .neqv. associated(profs(1)%co2)
          if (.not.l_req) l_req = ropts%rt_ir%n2o_data   .neqv. associated(profs(1)%n2o)
          if (.not.l_req) l_req = ropts%rt_ir%co_data    .neqv. associated(profs(1)%co)
          if (.not.l_req) l_req = ropts%rt_ir%ch4_data   .neqv. associated(profs(1)%ch4)
          if (.not.l_req) l_req = ropts%rt_mw%clw_data   .neqv. associated(profs(1)%clw)
          if (.not.l_req) l_req = ropts%rt_ir%addaerosl  .neqv. associated(profs(1)%aerosols)
          if (.not.l_req) l_req = ropts%rt_ir%addclouds  .neqv. associated(profs(1)%cloud)
        end if
      end if
      if (l_req) then
        call dealloc_rttov_arrays(status, profs =profs)
        if (status /= NO_ERROR) return
        call  alloc_rttov_arrays (status, n_profiles,n_levels, n_channels,&
                                  .true., ropts, profs =profs)
        if (status /= NO_ERROR) return
      end if
    end if

    if (present(transm)) then
      l_req = .not.associated(transm%tau_total)
      if (.not.l_req) l_req = (n_channels /= size(transm%tau_total) .and. rtifc_alloc_mode == 0) .or. &
                              (n_channels >  size(transm%tau_total) .and. rtifc_alloc_mode >= 1)
      if (.not.l_req) l_req = n_levels   /= size(transm%tau_levels,1)
      if (.not.l_req .and. transm%l_opdep) then
        if (associated(transm%opdep_ref)) then
          l_req = (size(transm%opdep_ref,1) < n_levels_coef-1)
        else
          l_req = .true.
        end if
      end if
      if (l_req) then
        call dealloc_rttov_arrays(status, transm=transm)
        if (status /= NO_ERROR) return
        call alloc_rttov_arrays (status, n_profiles, n_levels, n_channels,&
                                 .true., ropts, transm =transm, n_levels_coef=n_levels_coef)
        if (status /= NO_ERROR) return
      end if
    end if

    if (present(rads)) then
      l_req  = .not.associated(rads%clear)
      if (.not.l_req)  l_req  = (n_channels /= size(rads%clear) .and. rtifc_alloc_mode == 0) .or. &
                                (n_channels >  size(rads%clear) .and. rtifc_alloc_mode >= 1)
      if (.not.l_req)  l_req  = n_levels   /= size(rads%overcast,1) + 1
      if (present(rads2)) then
        l_req2 = .not.associated(rads2%upclear)
        if (.not.l_req2) l_req2 = (n_channels /= size(rads2%upclear) .and. rtifc_alloc_mode == 0) .or. &
                                  (n_channels >  size(rads2%upclear) .and. rtifc_alloc_mode >= 1)
        if (.not.l_req2) l_req2 = n_levels   /= size(rads2%up,1)+1
        if (l_req .and. l_req2) then
          call dealloc_rttov_arrays(status, rads=rads, rads2=rads2)
          if (status /= NO_ERROR) return
          call alloc_rttov_arrays (status, n_profiles, n_levels, n_channels,&
                                   .true., ropts, rads=rads, rads2=rads2)
          if (status /= NO_ERROR) return
        else if (l_req) then
          call dealloc_rttov_arrays(status, rads=rads)
          if (status /= NO_ERROR) return
          call alloc_rttov_arrays (status, n_profiles, n_levels, n_channels,&
                                  .true., ropts, rads=rads)
          if (status /= NO_ERROR) return
        end if
      else
        if (l_req) then
          call dealloc_rttov_arrays(status, rads=rads)
          if (status /= NO_ERROR) return
          call alloc_rttov_arrays (status, n_profiles, n_levels, n_channels,&
                                   .true., ropts, rads=rads)
          if (status /= NO_ERROR) return
        end if
      end if
    end if
    if (present(height)) then
      l_req = .not.associated(height)
      if (.not.l_req) l_req = (n_channels /= size(height) .and. rtifc_alloc_mode == 0) .or. &
                              (n_channels >  size(height) .and. rtifc_alloc_mode >= 1)
      if (l_req) then
        call dealloc_rttov_arrays(status, height=height)
        if (status /= NO_ERROR) return
        call alloc_rttov_arrays (status, n_profiles, n_levels, n_channels,&
                                 .true., ropts, height=height)
        if (status /= NO_ERROR) return
      end if
    end if

  end subroutine realloc_rttov_arrays


#if defined(_RTIFC_DISTRIBCOEF)
  !--------------------------------
  ! RTTOV IFC MPI Transfer routines
  !--------------------------------

  subroutine p_bcast_rttov_coefs (coefs, source, comm)
    type(rttov_coefs), intent(inout) :: coefs
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: comm
    !-------------------------------------------------------------------
    ! Broadcast an rttov_coefs structure across all available processors
    !-------------------------------------------------------------------
    call p_bcast(coefs% initialised,     source, comm)
    call p_bcast(coefs% coef,          source, comm)
    call p_bcast(coefs% coef_scatt_ir, source, comm)
    call p_bcast(coefs% optp,          source, comm)
    call p_bcast(coefs% coef_pccomp,   source, comm)
  end subroutine


  subroutine p_bcast_rttov_coef(coef, source, comm)
    type(rttov_coef),  intent(inout) :: coef
    integer,           intent(in)    :: source
    integer, optional, intent(in)    :: comm
    !------------------------------------------------------------------
    ! Broadcast an rttov_coef structure across all available processors
    !------------------------------------------------------------------
    integer :: i

    call p_bcast_rttov_container(coef, source, comm)

    if (pe_ifc /= source) then
      if (associated(coef%fmv_gas_id       )) allocate(coef%fmv_gas_id       (coef%fmv_gas))
      if (associated(coef%fmv_gas_pos      )) allocate(coef%fmv_gas_pos      (ngases_max))
      if (associated(coef%fmv_var          )) allocate(coef%fmv_var          (coef%fmv_gas))
      if (associated(coef%fmv_lvl          )) allocate(coef%fmv_lvl          (coef%fmv_gas))
      if (associated(coef%fmv_coe          )) allocate(coef%fmv_coe          (coef%fmv_gas))
      if (associated(coef%ff_ori_chn       )) allocate(coef%ff_ori_chn       (coef%fmv_chn))
      if (associated(coef%ff_val_chn       )) allocate(coef%ff_val_chn       (coef%fmv_chn))
      if (associated(coef%ff_cwn           )) allocate(coef%ff_cwn           (coef%fmv_chn))
      if (associated(coef%ff_bco           )) allocate(coef%ff_bco           (coef%fmv_chn))
      if (associated(coef%ff_bcs           )) allocate(coef%ff_bcs           (coef%fmv_chn))
      if (associated(coef%ff_gam           )) allocate(coef%ff_gam           (coef%fmv_chn))
      if (associated(coef%tt_val_chn       )) allocate(coef%tt_val_chn       (coef%fmv_chn))
      if (associated(coef%tt_a0            )) allocate(coef%tt_a0            (coef%fmv_chn))
      if (associated(coef%tt_a1            )) allocate(coef%tt_a1            (coef%fmv_chn))
      if (associated(coef%pw_val_chn       )) allocate(coef%pw_val_chn       (coef%fmv_chn))
      if (associated(coef%ss_val_chn       )) allocate(coef%ss_val_chn       (coef%fmv_chn))
      if (associated(coef%ss_solar_spectrum)) allocate(coef%ss_solar_spectrum(coef%fmv_chn))
      if (associated(coef%refl_visnir_ow   )) allocate(coef%refl_visnir_ow   (coef%fmv_chn))
      if (associated(coef%refl_visnir_fw   )) allocate(coef%refl_visnir_fw   (coef%fmv_chn))
      if (associated(coef%woc_waopc_ow     )) allocate(coef%woc_waopc_ow     (coef%fmv_chn))
      if (associated(coef%woc_waopc_fw     )) allocate(coef%woc_waopc_fw     (coef%fmv_chn))
      if (associated(coef%ws_npoint        )) allocate(coef%ws_npoint        (coef%ws_nomega))
      if (associated(coef%ws_k_omega       )) allocate(coef%ws_k_omega       (coef%ws_nomega))
      if (associated(coef%fastem_polar     )) allocate(coef%fastem_polar     (coef%fmv_chn))
      if (associated(coef%ssirem_a0        )) allocate(coef%ssirem_a0        (coef%fmv_chn))
      if (associated(coef%ssirem_a1        )) allocate(coef%ssirem_a1        (coef%fmv_chn))
      if (associated(coef%ssirem_a2        )) allocate(coef%ssirem_a2        (coef%fmv_chn))
      if (associated(coef%ssirem_xzn1      )) allocate(coef%ssirem_xzn1      (coef%fmv_chn))
      if (associated(coef%ssirem_xzn2      )) allocate(coef%ssirem_xzn2      (coef%fmv_chn))
      if (associated(coef%iremis_coef      )) allocate(coef%iremis_coef      (coef%iremis_ncoef,coef%fmv_chn))
    endif

    if (associated(coef%fmv_gas_id       )) call p_bcast(coef%fmv_gas_id,       source,comm)
    if (associated(coef%fmv_gas_pos      )) call p_bcast(coef%fmv_gas_pos,      source,comm)
    if (associated(coef%fmv_var          )) call p_bcast(coef%fmv_var,          source,comm)
    if (associated(coef%fmv_lvl          )) call p_bcast(coef%fmv_lvl,          source,comm)
    if (associated(coef%fmv_coe          )) call p_bcast(coef%fmv_coe,          source,comm)
    if (associated(coef%ff_ori_chn       )) call p_bcast(coef%ff_ori_chn,       source,comm)
    if (associated(coef%ff_val_chn       )) call p_bcast(coef%ff_val_chn,       source,comm)
    if (associated(coef%ff_cwn           )) call p_bcast(coef%ff_cwn,           source,comm)
    if (associated(coef%ff_bco           )) call p_bcast(coef%ff_bco,           source,comm)
    if (associated(coef%ff_bcs           )) call p_bcast(coef%ff_bcs,           source,comm)
    if (associated(coef%ff_gam           )) call p_bcast(coef%ff_gam,           source,comm)
    if (associated(coef%tt_val_chn       )) call p_bcast(coef%tt_val_chn,       source,comm)
    if (associated(coef%tt_a0            )) call p_bcast(coef%tt_a0,            source,comm)
    if (associated(coef%tt_a1            )) call p_bcast(coef%tt_a1,            source,comm)
    if (associated(coef%pw_val_chn       )) call p_bcast(coef%pw_val_chn,       source,comm)
    if (associated(coef%ss_val_chn       )) call p_bcast(coef%ss_val_chn,       source,comm)
    if (associated(coef%ss_solar_spectrum)) call p_bcast(coef%ss_solar_spectrum,source,comm)
    if (associated(coef%refl_visnir_ow   )) call p_bcast(coef%refl_visnir_ow,   source,comm)
    if (associated(coef%refl_visnir_fw   )) call p_bcast(coef%refl_visnir_fw,   source,comm)
    if (associated(coef%woc_waopc_ow     )) call p_bcast(coef%woc_waopc_ow,     source,comm)
    if (associated(coef%woc_waopc_fw     )) call p_bcast(coef%woc_waopc_fw,     source,comm)
    if (associated(coef%ws_npoint        )) call p_bcast(coef%ws_npoint,        source,comm)
    if (associated(coef%ws_k_omega       )) call p_bcast(coef%ws_k_omega,       source,comm)
    if (associated(coef%fastem_polar     )) call p_bcast(coef%fastem_polar,     source,comm)
    if (associated(coef%ssirem_a0        )) call p_bcast(coef%ssirem_a0,        source,comm)
    if (associated(coef%ssirem_a1        )) call p_bcast(coef%ssirem_a1,        source,comm)
    if (associated(coef%ssirem_a2        )) call p_bcast(coef%ssirem_a2,        source,comm)
    if (associated(coef%ssirem_xzn1      )) call p_bcast(coef%ssirem_xzn1,      source,comm)
    if (associated(coef%ssirem_xzn2      )) call p_bcast(coef%ssirem_xzn2,      source,comm)
    if (associated(coef%iremis_coef      )) call p_bcast(coef%iremis_coef,      source,comm)

    if (pe_ifc /= source) then
      if (associated(coef%bkg_prfl_mr   )) allocate(coef%bkg_prfl_mr   (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%env_prfl_tmax )) allocate(coef%env_prfl_tmax (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%env_prfl_tmin )) allocate(coef%env_prfl_tmin (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%env_prfl_gmax )) allocate(coef%env_prfl_gmax (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%env_prfl_gmin )) allocate(coef%env_prfl_gmin (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%ref_prfl_p    )) allocate(coef%ref_prfl_p    (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%ref_prfl_t    )) allocate(coef%ref_prfl_t    (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%ref_prfl_mr   )) allocate(coef%ref_prfl_mr   (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%lim_prfl_p    )) allocate(coef%lim_prfl_p    (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%lim_prfl_tmax )) allocate(coef%lim_prfl_tmax (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%lim_prfl_tmin )) allocate(coef%lim_prfl_tmin (coef%fmv_lvl(gas_id_mixed)              ))
      if (associated(coef%lim_prfl_gmax )) allocate(coef%lim_prfl_gmax (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%lim_prfl_gmin )) allocate(coef%lim_prfl_gmin (coef%fmv_lvl(gas_id_mixed), coef%fmv_gas))
      if (associated(coef%thermal       )) allocate(coef%thermal      (coef%fmv_chn))
      if (coef%solarcoef) then
         if (associated(coef%solar      )) allocate(coef%solar        (coef%fmv_chn))
      else
         coef%solar => coef%thermal
      endif
      if (associated(coef%nlte_coef     )) allocate(coef%nlte_coef                  )
      if (associated(coef%pmc_pnominal  )) allocate(coef%pmc_pnominal (coef%fmv_chn))
      if (associated(coef%pmc_coef      )) allocate(coef%pmc_coef     (coef%pmc_nlay, coef%fmv_chn, coef%pmc_nvar))
      if (associated(coef%pmc_ppmc      )) allocate(coef%pmc_ppmc     (coef%fmv_chn))
      if (associated(coef%planck1       )) allocate(coef%planck1      (coef%fmv_chn))
      if (associated(coef%planck2       )) allocate(coef%planck2      (coef%fmv_chn))
      if (associated(coef%frequency_ghz )) allocate(coef%frequency_ghz(coef%fmv_chn))
      if (associated(coef%dp            )) allocate(coef%dp           (coef%nlayers))
      if (associated(coef%dpp           )) allocate(coef%dpp          (0:coef%nlayers))
      if (associated(coef%tstar         )) allocate(coef%tstar        (coef%nlayers))
      if (associated(coef%to3star       )) allocate(coef%to3star      (coef%nlayers))
      if (associated(coef%wstar         )) allocate(coef%wstar        (coef%nlayers))
      if (associated(coef%ostar         )) allocate(coef%ostar        (coef%nlayers))
      if (associated(coef%co2star       )) allocate(coef%co2star      (coef%nlayers))
      if (associated(coef%n2ostar       )) allocate(coef%n2ostar      (coef%nlayers))
      if (associated(coef%costar        )) allocate(coef%costar       (coef%nlayers))
      if (associated(coef%ch4star       )) allocate(coef%ch4star      (coef%nlayers))
      if (associated(coef%so2star       )) allocate(coef%so2star      (coef%nlayers))
      if (associated(coef%bounds        )) allocate(coef%bounds       (2, coef%fmv_gas, coef%fmv_chn, 2))
    endif

    if (associated(coef%bkg_prfl_mr   )) call p_bcast(coef%bkg_prfl_mr,   source,comm)
    if (associated(coef%env_prfl_tmax )) call p_bcast(coef%env_prfl_tmax, source,comm)
    if (associated(coef%env_prfl_tmin )) call p_bcast(coef%env_prfl_tmin, source,comm)
    if (associated(coef%env_prfl_gmax )) call p_bcast(coef%env_prfl_gmax, source,comm)
    if (associated(coef%env_prfl_gmin )) call p_bcast(coef%env_prfl_gmin, source,comm)
    if (associated(coef%ref_prfl_p    )) call p_bcast(coef%ref_prfl_p,    source,comm)
    if (associated(coef%ref_prfl_t    )) call p_bcast(coef%ref_prfl_t,    source,comm)
    if (associated(coef%ref_prfl_mr   )) call p_bcast(coef%ref_prfl_mr,   source,comm)
    if (associated(coef%lim_prfl_p    )) call p_bcast(coef%lim_prfl_p,    source,comm)
    if (associated(coef%lim_prfl_tmax )) call p_bcast(coef%lim_prfl_tmax, source,comm)
    if (associated(coef%lim_prfl_tmin )) call p_bcast(coef%lim_prfl_tmin, source,comm)
    if (associated(coef%lim_prfl_gmax )) call p_bcast(coef%lim_prfl_gmax, source,comm)
    if (associated(coef%lim_prfl_gmin )) call p_bcast(coef%lim_prfl_gmin, source,comm)
    if (associated(coef%thermal       )) then
       do i = 1, size(coef%thermal)
          call p_bcast(coef%thermal(i), source,comm)
       end do
    end if
    if (coef%solarcoef) then
       if (associated(coef%solar      )) then
          do i = 1, size(coef%solar)
             call p_bcast(coef%solar(i), source,comm)
          end do
       end if
    end if

    if (associated(coef%nlte_coef     )) call p_bcast(coef%nlte_coef,     source,comm)
    if (associated(coef%pmc_pnominal  )) call p_bcast(coef%pmc_pnominal,  source,comm)
    if (associated(coef%pmc_coef      )) call p_bcast(coef%pmc_coef,      source,comm)
    if (associated(coef%pmc_ppmc      )) call p_bcast(coef%pmc_ppmc,      source,comm)
    if (associated(coef%planck1       )) call p_bcast(coef%planck1,       source,comm)
    if (associated(coef%planck2       )) call p_bcast(coef%planck2,       source,comm)
    if (associated(coef%frequency_ghz )) call p_bcast(coef%frequency_ghz, source,comm)
    if (associated(coef%dp            )) call p_bcast(coef%dp,            source,comm)
    if (associated(coef%dpp           )) call p_bcast(coef%dpp,           source,comm)
    if (associated(coef%tstar         )) call p_bcast(coef%tstar,         source,comm)
    if (associated(coef%to3star       )) call p_bcast(coef%to3star,       source,comm)
    if (associated(coef%wstar         )) call p_bcast(coef%wstar,         source,comm)
    if (associated(coef%ostar         )) call p_bcast(coef%ostar,         source,comm)
    if (associated(coef%co2star       )) call p_bcast(coef%co2star,       source,comm)
    if (associated(coef%n2ostar       )) call p_bcast(coef%n2ostar,       source,comm)
    if (associated(coef%costar        )) call p_bcast(coef%costar,        source,comm)
    if (associated(coef%ch4star       )) call p_bcast(coef%ch4star,       source,comm)
    if (associated(coef%so2star       )) call p_bcast(coef%so2star,       source,comm)
    if (associated(coef%bounds        )) call p_bcast(coef%bounds,        source,comm)

  end subroutine


  subroutine p_bcast_rttov_coef_scatt_ir(coef, source, comm)
    type(rttov_coef_scatt_ir),intent(inout) :: coef
    integer,                  intent(in)    :: source
    integer, optional,        intent(in)    :: comm
    !------------------------------------------------------------------
    ! Broadcast an rttov_coef structure across all available processors
    !------------------------------------------------------------------
    integer :: count
    integer :: dimensions(21)


    dimensions(:) = -1
    call p_bcast(dimensions,source,comm)

    call p_bcast_rttov_container (coef, source, comm)

    dimensions(:) = -1

    if (associated(coef%fmv_aer_rh_val)) dimensions(1:1)   = shape(coef%fmv_aer_rh_val)
    if (associated(coef%fmv_wcl_rh_val)) dimensions(2:2)   = shape(coef%fmv_wcl_rh_val)
    if (associated(coef%abs))            dimensions(4:5)   = shape(coef%abs)
    if (associated(coef%sca))            dimensions(6:7)   = shape(coef%sca)
    if (associated(coef%bpr))            dimensions(8:9)   = shape(coef%bpr)
    if (associated(coef%pha))            dimensions(10:12) = shape(coef%pha)
    if (associated(coef%fmv_icl_deff ))  dimensions(13:13) = shape(coef%fmv_icl_deff)
    if (associated(coef%aer_pha_index))  dimensions(14:14) = shape(coef%aer_pha_index)
    if (associated(coef%wcl_pha_index))  dimensions(15:15) = shape(coef%wcl_pha_index)
    if (associated(coef%icl_pha_index))  dimensions(16:16) = shape(coef%icl_pha_index)
    if (associated(coef%nmom         ))  dimensions(17:18) = shape(coef%nmom)
    if (associated(coef%legcoef      ))  dimensions(19:21) = shape(coef%legcoef)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef%fmv_aer_comp_name )) allocate(coef%fmv_aer_comp_name (coef%fmv_aer_comp))
      if (associated(coef%fmv_wcl_comp_name )) allocate(coef%fmv_wcl_comp_name (coef%fmv_wcl_comp))
      if (associated(coef%fmv_aer_rh        )) allocate(coef%fmv_aer_rh        (coef%fmv_aer_comp))
      if (associated(coef%fmv_wcl_rh        )) allocate(coef%fmv_wcl_rh        (coef%fmv_wcl_comp))
      if (associated(coef%fmv_aer_rh_val    )) allocate(coef%fmv_aer_rh_val    (dimensions(1)))
      if (associated(coef%fmv_wcl_rh_val    )) allocate(coef%fmv_wcl_rh_val    (dimensions(2)))
      if (associated(coef%fmv_icl_deff      )) allocate(coef%fmv_icl_deff      (dimensions(13)))
      if (associated(coef%aer_pha_chanlist  )) allocate(coef%aer_pha_chanlist  (coef%fmv_aer_pha_chn))
      if (associated(coef%wcl_pha_chanlist  )) allocate(coef%wcl_pha_chanlist  (coef%fmv_wcl_pha_chn))
      if (associated(coef%icl_pha_chanlist  )) allocate(coef%icl_pha_chanlist  (coef%fmv_icl_pha_chn))
      if (associated(coef%aer_pha_index     )) allocate(coef%aer_pha_index     (dimensions(14)))
      if (associated(coef%wcl_pha_index     )) allocate(coef%wcl_pha_index     (dimensions(15)))
      if (associated(coef%icl_pha_index     )) allocate(coef%icl_pha_index     (dimensions(16)))
      if (associated(coef%aer_phangle       )) allocate(coef%aer_phangle       (coef%aer_nphangle))
      if (associated(coef%wcl_phangle       )) allocate(coef%wcl_phangle       (coef%wcl_nphangle))
      if (associated(coef%icl_phangle       )) allocate(coef%icl_phangle       (coef%icl_nphangle))
      if (associated(coef%nmom              )) allocate(coef%nmom              (dimensions(17),dimensions(18)))
      if (associated(coef%legcoef           )) allocate(coef%legcoef           (dimensions(19),dimensions(20),dimensions(21)))
      if (associated(coef%aer_mmr2nd        )) allocate(coef%aer_mmr2nd        (coef%fmv_aer_comp))
      if (associated(coef%abs               )) allocate(coef%abs               (dimensions(4),dimensions(5)))
      if (associated(coef%sca               )) allocate(coef%sca               (dimensions(6),dimensions(7)))
      if (associated(coef%bpr               )) allocate(coef%bpr               (dimensions(8),dimensions(9)))
      if (associated(coef%pha               )) allocate(coef%pha               (dimensions(10),dimensions(11),dimensions(12)))
      if (associated(coef%confac            )) allocate(coef%confac            (coef%fmv_wcl_comp))
    endif

    if (associated(coef%fmv_aer_comp_name )) call p_bcast(coef%fmv_aer_comp_name, source,comm)
    if (associated(coef%fmv_wcl_comp_name )) call p_bcast(coef%fmv_wcl_comp_name, source,comm)
    if (associated(coef%fmv_aer_rh        )) call p_bcast(coef%fmv_aer_rh,        source,comm)
    if (associated(coef%fmv_wcl_rh        )) call p_bcast(coef%fmv_wcl_rh,        source,comm)
    if (associated(coef%fmv_aer_rh_val    )) call p_bcast(coef%fmv_aer_rh_val,    source,comm)
    if (associated(coef%fmv_wcl_rh_val    )) call p_bcast(coef%fmv_wcl_rh_val,    source,comm)
    if (associated(coef%fmv_icl_deff      )) call p_bcast(coef%fmv_icl_deff,      source,comm)
    if (associated(coef%aer_pha_chanlist  )) call p_bcast(coef%aer_pha_chanlist,  source,comm)
    if (associated(coef%wcl_pha_chanlist  )) call p_bcast(coef%wcl_pha_chanlist,  source,comm)
    if (associated(coef%icl_pha_chanlist  )) call p_bcast(coef%icl_pha_chanlist,  source,comm)
    if (associated(coef%aer_pha_index     )) call p_bcast(coef%aer_pha_index,     source,comm)
    if (associated(coef%wcl_pha_index     )) call p_bcast(coef%wcl_pha_index,     source,comm)
    if (associated(coef%icl_pha_index     )) call p_bcast(coef%icl_pha_index,     source,comm)
    if (associated(coef%aer_phangle       )) call p_bcast(coef%aer_phangle,       source,comm)
    if (associated(coef%wcl_phangle       )) call p_bcast(coef%wcl_phangle,       source,comm)
    if (associated(coef%icl_phangle       )) call p_bcast(coef%icl_phangle,       source,comm)
    if (associated(coef%nmom              )) call p_bcast(coef%nmom,              source,comm)
    if (associated(coef%legcoef           )) call p_bcast(coef%legcoef,           source,comm)
    if (associated(coef%aer_mmr2nd        )) call p_bcast(coef%aer_mmr2nd,        source,comm)
    if (associated(coef%abs               )) call p_bcast(coef%abs,               source,comm)
    if (associated(coef%sca               )) call p_bcast(coef%sca,               source,comm)
    if (associated(coef%bpr               )) call p_bcast(coef%bpr,               source,comm)
    if (associated(coef%pha               )) call p_bcast(coef%pha,               source,comm)
    if (associated(coef%confac            )) call p_bcast(coef%confac,            source,comm)

    call p_bcast(coef%aer_phfn_int,source,comm)
    call p_bcast(coef%wcl_phfn_int,source,comm)
    call p_bcast(coef%icl_phfn_int,source,comm)

  end subroutine


  subroutine p_bcast_rttov_optpar_ir(optp, source, comm)
    type(rttov_optpar_ir),intent(inout) :: optp
    integer,              intent(in)    :: source
    integer, optional,    intent(in)    :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir structure across all available processors
    !-----------------------------------------------------------------------
    integer :: run1
    integer :: dimensions(3)


    call p_bcast_rttov_container (optp, source, comm)

    dimensions(:) = -1
    if(associated(optp% optpaer)) dimensions(1:1) = shape(optp% optpaer)
    if(associated(optp% optpwcl)) dimensions(2:2) = shape(optp% optpwcl)
    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if(associated(optp% optpaer)) allocate(optp% optpaer(dimensions(1)))
      if(associated(optp% optpwcl)) allocate(optp% optpwcl(dimensions(2)))
      if(associated(optp% optpicl )) allocate(optp% optpicl )
      if(associated(optp% optpiclb)) allocate(optp% optpiclb)
    endif

    if(associated(optp% optpaer)) then
      do run1=1,dimensions(1)
        call p_bcast(optp% optpaer(run1),source,comm)
      enddo
    endif

    if(associated(optp% optpwcl)) then
      do run1=1,dimensions(2)
        call p_bcast(optp% optpwcl(run1),source,comm)
      enddo
    endif

    if(associated(optp% optpicl)) then
      call p_bcast(optp% optpicl, source,comm)
      call p_bcast(optp% optpiclb,source,comm)
    endif
  end subroutine

  subroutine p_bcast_rttov_coef_optpiclb(coef, source, comm)
  type(rttov_coef_optpiclb),intent(inout) :: coef
  integer,                  intent(in)    :: source
  integer,     optional,    intent(in)    :: comm

    integer :: dimensions(1)

    call p_bcast_rttov_container (coef, source, comm)

    dimensions = -1

    if (associated(coef%iwn2014)) dimensions(1:1) = shape(coef%iwn2014)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef%iwn2014   )) allocate(coef%iwn2014   (dimensions(1)))
      if (associated(coef%jwn2014   )) allocate(coef%jwn2014   (dimensions(1)))
      if (associated(coef%dx_dwn2014)) allocate(coef%dx_dwn2014(dimensions(1)))
      if (associated(coef%q         )) allocate(coef%q         (baran_ngauss))
      if (associated(coef%w         )) allocate(coef%w         (baran_ngauss))
   end if

    if (associated(coef%iwn2014   )) call p_bcast(coef%iwn2014,   source,comm)
    if (associated(coef%jwn2014   )) call p_bcast(coef%jwn2014,   source,comm)
    if (associated(coef%dx_dwn2014)) call p_bcast(coef%dx_dwn2014,source,comm)
    if (associated(coef%q         )) call p_bcast(coef%q,         source,comm)
    if (associated(coef%w         )) call p_bcast(coef%w,         source,comm)

   call p_bcast(coef%baran_phfn_int, source, comm)

  end subroutine p_bcast_rttov_coef_optpiclb


  subroutine p_bcast_rttov_phasefn_int(phfn, source, comm)
    type(rttov_phasefn_int),intent(inout) :: phfn
    integer,                intent(in)    :: source
    integer,     optional,  intent(in)    :: comm

    integer :: dimensions(2)

    call p_bcast_rttov_container (phfn, source, comm)

    dimensions = -1

    if (associated(phfn%cosphangle)) dimensions(1:1) = shape(phfn%cosphangle)
    if (associated(phfn%iphangle  )) dimensions(2:2) = shape(phfn%iphangle)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(phfn%cosphangle)) allocate(phfn%cosphangle(dimensions(1)))
      if (associated(phfn%iphangle  )) allocate(phfn%iphangle  (dimensions(2)))
    end if

    if (associated(phfn%cosphangle)) call p_bcast(phfn%cosphangle,source,comm)
    if (associated(phfn%iphangle  )) call p_bcast(phfn%iphangle,  source,comm)

  end subroutine p_bcast_rttov_phasefn_int


  subroutine p_bcast_rttov_fast_coef(coef, source, comm)
    type(rttov_fast_coef),intent(inout) :: coef
    integer,              intent(in)    :: source
    integer,   optional,  intent(in)    :: comm

    integer :: dimensions(19)
    integer :: i

    call p_bcast_rttov_container (coef, source, comm)

    dimensions = -1

    if (associated(coef%mixedgas   )) dimensions( 1: 2) = shape (coef%mixedgas)
    if (associated(coef%watervapour)) dimensions( 3: 4) = shape (coef%watervapour)
    if (associated(coef%ozone      )) dimensions( 5: 6) = shape (coef%ozone)
    if (associated(coef%wvcont     )) dimensions( 7: 8) = shape (coef%wvcont)
    if (associated(coef%co2        )) dimensions( 9:10) = shape (coef%co2)
    if (associated(coef%n2o        )) dimensions(11:12) = shape (coef%n2o)
    if (associated(coef%co         )) dimensions(13:14) = shape (coef%co)
    if (associated(coef%ch4        )) dimensions(15:16) = shape (coef%ch4)
    if (associated(coef%co2        )) dimensions(17:18) = shape (coef%co2)
    if (associated(coef%gasarray   )) dimensions(19:19) = shape (coef%gasarray)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef%mixedgas   )) allocate(coef%mixedgas   (dimensions(1),dimensions(2)))
      if (associated(coef%watervapour)) allocate(coef%watervapour(dimensions(3),dimensions(4)))
      if (associated(coef%ozone      )) allocate(coef%ozone      (dimensions(5),dimensions(6)))
      if (associated(coef%wvcont     )) allocate(coef%wvcont     (dimensions(7),dimensions(8)))
      if (associated(coef%co2        )) allocate(coef%co2        (dimensions(9),dimensions(10)))
      if (associated(coef%n2o        )) allocate(coef%n2o        (dimensions(11),dimensions(12)))
      if (associated(coef%co         )) allocate(coef%co         (dimensions(13),dimensions(14)))
      if (associated(coef%ch4        )) allocate(coef%ch4        (dimensions(15),dimensions(16)))
      if (associated(coef%co2        )) allocate(coef%co2        (dimensions(17),dimensions(18)))
      if (associated(coef%gasarray   )) allocate(coef%gasarray   (dimensions(19)))
   end if

    if (associated(coef%mixedgas   )) call p_bcast(coef%mixedgas,   source,comm)
    if (associated(coef%watervapour)) call p_bcast(coef%watervapour,source,comm)
    if (associated(coef%ozone      )) call p_bcast(coef%ozone,      source,comm)
    if (associated(coef%wvcont     )) call p_bcast(coef%wvcont,     source,comm)
    if (associated(coef%co2        )) call p_bcast(coef%co2,        source,comm)
    if (associated(coef%n2o        )) call p_bcast(coef%n2o,        source,comm)
    if (associated(coef%co         )) call p_bcast(coef%co,         source,comm)
    if (associated(coef%ch4        )) call p_bcast(coef%ch4,        source,comm)
    if (associated(coef%co2        )) call p_bcast(coef%co2,        source,comm)
    if (associated(coef%gasarray   )) then
       do i = 1, dimensions(19)
          call p_bcast(coef%gasarray(i),source,comm)
       end do
    end if

  end subroutine p_bcast_rttov_fast_coef


  subroutine p_bcast_rttov_fast_coef_gas(coef, source, comm)
    type(rttov_fast_coef_gas),intent(inout) :: coef
    integer,                  intent(in)    :: source
    integer,       optional,  intent(in)    :: comm

    integer :: dimensions(2)

    call p_bcast_rttov_container (coef, source, comm)

    dimensions = -1

    if (associated(coef%coef)) dimensions(1:2) = shape(coef%coef)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef%coef)) allocate(coef%coef(dimensions(1),dimensions(2)))
   end if

    if (associated(coef%coef)) call p_bcast(coef%coef,source,comm)

  end subroutine p_bcast_rttov_fast_coef_gas


  subroutine p_bcast_rttov_nlte_coef(coef, source, comm)
    type(rttov_nlte_coef),intent(inout) :: coef
    integer,              intent(in)    :: source
    integer,   optional,  intent(in)    :: comm

    call p_bcast_rttov_container (coef, source, comm)

    if (pe_ifc /= source) then
      if (associated(coef%coef         )) allocate(coef%coef         (coef%ncoef,coef%nsat,coef%nsol,coef%nchan))
      if (associated(coef%sol_zen_angle)) allocate(coef%sol_zen_angle(coef%nsol))
      if (associated(coef%sat_zen_angle)) allocate(coef%sat_zen_angle(coef%nsat))
      if (associated(coef%cos_sol      )) allocate(coef%cos_sol      (coef%nsol))
      if (associated(coef%sec_sat      )) allocate(coef%sec_sat      (coef%nsat))
   end if

   if (associated(coef%coef         )) call p_bcast(coef%coef,         source,comm)
   if (associated(coef%sol_zen_angle)) call p_bcast(coef%sol_zen_angle,source,comm)
   if (associated(coef%sat_zen_angle)) call p_bcast(coef%sat_zen_angle,source,comm)
   if (associated(coef%cos_sol      )) call p_bcast(coef%cos_sol,      source,comm)
   if (associated(coef%sec_sat      )) call p_bcast(coef%sec_sat,      source,comm)

  end subroutine p_bcast_rttov_nlte_coef


  subroutine p_bcast_rttov_coef_pccomp(coef_pccomp,source,comm)
    type(rttov_coef_pccomp),intent(inout) :: coef_pccomp
    integer,              intent(in)      :: source
    integer, optional,    intent(in)      :: comm
    !-------------------------------------------------------------------------
    ! Broadcast an rttov_coef_pccomp structure across all available processors
    !-------------------------------------------------------------------------
    integer :: run1, run2
    integer :: dimensions(4)

    call p_bcast_rttov_container (coef_pccomp, source, comm)

    dimensions = -1

    if (associated(coef_pccomp%co2_pc_ref)) dimensions(1:1) = shape(coef_pccomp%co2_pc_ref)
    if (associated(coef_pccomp%n2o_pc_ref)) dimensions(2:2) = shape(coef_pccomp%n2o_pc_ref)
    if (associated(coef_pccomp%co_pc_ref))  dimensions(3:3) = shape(coef_pccomp%co_pc_ref)
    if (associated(coef_pccomp%ch4_pc_ref)) dimensions(4:4) = shape(coef_pccomp%ch4_pc_ref)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef_pccomp%fmv_pc_sets      )) allocate(coef_pccomp%fmv_pc_sets      (coef_pccomp%fmv_pc_bands))
      if (associated(coef_pccomp%emiss_chn        )) allocate(coef_pccomp%emiss_chn        (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c1         )) allocate(coef_pccomp%emiss_c1         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c2         )) allocate(coef_pccomp%emiss_c2         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c3         )) allocate(coef_pccomp%emiss_c3         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c4         )) allocate(coef_pccomp%emiss_c4         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c5         )) allocate(coef_pccomp%emiss_c5         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c6         )) allocate(coef_pccomp%emiss_c6         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c7         )) allocate(coef_pccomp%emiss_c7         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c8         )) allocate(coef_pccomp%emiss_c8         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%emiss_c9         )) allocate(coef_pccomp%emiss_c9         (coef_pccomp%fmv_pc_nche))
      if (associated(coef_pccomp%ref_pc_prfl_p    )) allocate(coef_pccomp%ref_pc_prfl_p    (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%ref_pc_prfl_mr   )) allocate(coef_pccomp%ref_pc_prfl_mr   (coef_pccomp%fmv_pc_nlev,&
                                                                                            coef_pccomp%fmv_pc_gas))
      if (associated(coef_pccomp%lim_pc_prfl_tmin )) allocate(coef_pccomp%lim_pc_prfl_tmin (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_tmax )) allocate(coef_pccomp%lim_pc_prfl_tmax (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_qmin )) allocate(coef_pccomp%lim_pc_prfl_qmin (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_qmax )) allocate(coef_pccomp%lim_pc_prfl_qmax (coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_ozmin)) allocate(coef_pccomp%lim_pc_prfl_ozmin(coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%lim_pc_prfl_ozmax)) allocate(coef_pccomp%lim_pc_prfl_ozmax(coef_pccomp%fmv_pc_nlev))
      if (associated(coef_pccomp%co2_pc_ref       )) allocate(coef_pccomp%co2_pc_ref       (dimensions(1)))
      if (associated(coef_pccomp%n2o_pc_ref       )) allocate(coef_pccomp%n2o_pc_ref       (dimensions(2)))
      if (associated(coef_pccomp%co_pc_ref        )) allocate(coef_pccomp%co_pc_ref        (dimensions(3)))
      if (associated(coef_pccomp%ch4_pc_ref       )) allocate(coef_pccomp%ch4_pc_ref       (dimensions(4)))
      if (associated(coef_pccomp%noise            )) allocate(coef_pccomp%noise            (coef_pccomp%fmv_pc_nchn_noise))
      if (associated(coef_pccomp%noise_in         )) allocate(coef_pccomp%noise_in         (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_ori_chn_in    )) allocate(coef_pccomp%ff_ori_chn_in    (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_cwn_in        )) allocate(coef_pccomp%ff_cwn_in        (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_bco_in        )) allocate(coef_pccomp%ff_bco_in        (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%ff_bcs_in        )) allocate(coef_pccomp%ff_bcs_in        (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%planck1_in       )) allocate(coef_pccomp%planck1_in       (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%planck2_in       )) allocate(coef_pccomp%planck2_in       (coef_pccomp%fmv_pc_nchn))
      if (associated(coef_pccomp%pcreg            )) allocate(coef_pccomp%pcreg            (coef_pccomp%fmv_pc_bands,&
                                                                                            coef_pccomp%fmv_pc_msets))
      if (associated(coef_pccomp%eigen            )) allocate(coef_pccomp%eigen            (coef_pccomp%fmv_pc_bands))
    endif

    if (associated(coef_pccomp%fmv_pc_sets      )) call p_bcast(coef_pccomp%fmv_pc_sets,      source,comm)
    if (associated(coef_pccomp%emiss_chn        )) call p_bcast(coef_pccomp%emiss_chn,        source,comm)
    if (associated(coef_pccomp%emiss_c1         )) call p_bcast(coef_pccomp%emiss_c1,         source,comm)
    if (associated(coef_pccomp%emiss_c2         )) call p_bcast(coef_pccomp%emiss_c2,         source,comm)
    if (associated(coef_pccomp%emiss_c3         )) call p_bcast(coef_pccomp%emiss_c3,         source,comm)
    if (associated(coef_pccomp%emiss_c4         )) call p_bcast(coef_pccomp%emiss_c4,         source,comm)
    if (associated(coef_pccomp%emiss_c5         )) call p_bcast(coef_pccomp%emiss_c5,         source,comm)
    if (associated(coef_pccomp%emiss_c6         )) call p_bcast(coef_pccomp%emiss_c6,         source,comm)
    if (associated(coef_pccomp%emiss_c7         )) call p_bcast(coef_pccomp%emiss_c7,         source,comm)
    if (associated(coef_pccomp%emiss_c8         )) call p_bcast(coef_pccomp%emiss_c8,         source,comm)
    if (associated(coef_pccomp%emiss_c9         )) call p_bcast(coef_pccomp%emiss_c9,         source,comm)
    if (associated(coef_pccomp%ref_pc_prfl_p    )) call p_bcast(coef_pccomp%ref_pc_prfl_p,    source,comm)
    if (associated(coef_pccomp%ref_pc_prfl_mr   )) call p_bcast(coef_pccomp%ref_pc_prfl_mr,   source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_tmin )) call p_bcast(coef_pccomp%lim_pc_prfl_tmin, source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_tmax )) call p_bcast(coef_pccomp%lim_pc_prfl_tmax, source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_qmin )) call p_bcast(coef_pccomp%lim_pc_prfl_qmin, source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_qmax )) call p_bcast(coef_pccomp%lim_pc_prfl_qmax, source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_ozmin)) call p_bcast(coef_pccomp%lim_pc_prfl_ozmin,source,comm)
    if (associated(coef_pccomp%lim_pc_prfl_ozmax)) call p_bcast(coef_pccomp%lim_pc_prfl_ozmax,source,comm)
    if (associated(coef_pccomp%co2_pc_ref       )) call p_bcast(coef_pccomp%co2_pc_ref,        source,comm)
    if (associated(coef_pccomp%n2o_pc_ref       )) call p_bcast(coef_pccomp%n2o_pc_ref,        source,comm)
    if (associated(coef_pccomp%co_pc_ref        )) call p_bcast(coef_pccomp%co_pc_ref,         source,comm)
    if (associated(coef_pccomp%ch4_pc_ref       )) call p_bcast(coef_pccomp%ch4_pc_ref,        source,comm)
    if (associated(coef_pccomp%noise            )) call p_bcast(coef_pccomp%noise,             source,comm)
    if (associated(coef_pccomp%noise_in         )) call p_bcast(coef_pccomp%noise_in,          source,comm)
    if (associated(coef_pccomp%ff_ori_chn_in    )) call p_bcast(coef_pccomp%ff_ori_chn_in,     source,comm)
    if (associated(coef_pccomp%ff_cwn_in        )) call p_bcast(coef_pccomp%ff_cwn_in,         source,comm)
    if (associated(coef_pccomp%ff_bco_in        )) call p_bcast(coef_pccomp%ff_bco_in,         source,comm)
    if (associated(coef_pccomp%ff_bcs_in        )) call p_bcast(coef_pccomp%ff_bcs_in,         source,comm)
    if (associated(coef_pccomp%planck1_in       )) call p_bcast(coef_pccomp%planck1_in,        source,comm)
    if (associated(coef_pccomp%planck2_in       )) call p_bcast(coef_pccomp%planck2_in,        source,comm)

    if (associated(coef_pccomp%eigen)) then
      do run1 = 1,coef_pccomp%fmv_pc_bands
        call p_bcast(coef_pccomp%eigen(run1),source,comm)
      enddo
    endif
    if (associated(coef_pccomp%pcreg)) then
      do run1 = 1,coef_pccomp%fmv_pc_bands
         do run2 = 1,coef_pccomp%fmv_pc_msets
            call p_bcast(coef_pccomp%pcreg(run1,run2),source,comm)
         end do
      enddo
    endif

  end subroutine p_bcast_rttov_coef_pccomp


  subroutine p_bcast_rttov_coef_pccomp1(coef_pccomp1,source,comm)
    type(rttov_coef_pccomp1),intent(inout) :: coef_pccomp1
    integer,              intent(in)      :: source
    integer, optional,    intent(in)      :: comm
    !-------------------------------------------------------------------------
    ! Broadcast an rttov_coef_pccomp structure across all available processors
    !-------------------------------------------------------------------------
    integer :: run1
    integer :: dimensions(2)

    call p_bcast_rttov_container (coef_pccomp1, source, comm)

    if (associated(coef_pccomp1%coefficients)) &
      dimensions(1:2) = shape(coef_pccomp1%coefficients)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef_pccomp1%predictindex)) &
        allocate(coef_pccomp1%predictindex(coef_pccomp1%fmv_pc_npred))
      if (associated(coef_pccomp1%coefficients)) &
        allocate(coef_pccomp1%coefficients(dimensions(1),dimensions(2)))
    endif

    if (associated(coef_pccomp1%predictindex)) &
      call p_bcast(coef_pccomp1%predictindex, source, comm)
    if (associated(coef_pccomp1%coefficients)) &
      call p_bcast(coef_pccomp1%coefficients, source, comm)
  end subroutine p_bcast_rttov_coef_pccomp1


  subroutine p_bcast_rttov_coef_pccomp2(coef_pccomp2,source,comm)
    type(rttov_coef_pccomp2),intent(inout) :: coef_pccomp2
    integer,              intent(in)      :: source
    integer, optional,    intent(in)      :: comm
    !-------------------------------------------------------------------------
    ! Broadcast an rttov_coef_pccomp structure across all available processors
    !-------------------------------------------------------------------------
    integer :: run1
    integer :: dimensions(4)

    call p_bcast_rttov_container (coef_pccomp2, source, comm)

    if (associated(coef_pccomp2%eigenvectors  )) dimensions(1:2) = shape(coef_pccomp2%eigenvectors  )
    if (associated(coef_pccomp2%eigenvectors_t)) dimensions(3:4) = shape(coef_pccomp2%eigenvectors_t)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(coef_pccomp2%eigenvectors  )) allocate(coef_pccomp2%eigenvectors  (dimensions(1),dimensions(2)))
      if (associated(coef_pccomp2%eigenvectors_t)) allocate(coef_pccomp2%eigenvectors_t(dimensions(3),dimensions(4)))
    endif

    if (associated(coef_pccomp2%eigenvectors  )) call p_bcast(coef_pccomp2%eigenvectors,   source, comm)
    if (associated(coef_pccomp2%eigenvectors_t)) call p_bcast(coef_pccomp2%eigenvectors_t, source, comm)

  end subroutine p_bcast_rttov_coef_pccomp2


  subroutine p_bcast_rttov_cnt_coef(buffer,source,comm)
    type(rttov_coef), intent(inout)   :: buffer
    integer,          intent(in)      :: source
    integer, optional,intent(in)      :: comm
    !------------------------------------------------------------------
    ! Broadcast an rttov_coef container across all available processors
    !------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)

    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_coef@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_coef


  subroutine p_bcast_rttov_cnt_coef_scatt_ir(buffer,source,comm)
    type(rttov_coef_scatt_ir), intent(inout) :: buffer
    integer,          intent(in)             :: source
    integer, optional,intent(in)             :: comm
    !---------------------------------------------------------------------------
    ! Broadcast an rttov_coef_scatt_ir container across all available processors
    !---------------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_coef_scatt_ir@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_coef_scatt_ir


  subroutine p_bcast_rttov_cnt_optpar_ir(buffer,source,comm)
    type(rttov_optpar_ir), intent(inout) :: buffer
    integer,          intent(in)         :: source
    integer, optional,intent(in)         :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir container across all available processors
    !-----------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_optpar_ir@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_optpar_ir


  subroutine p_bcast_rttov_cnt_coef_optpiclb(buffer,source,comm)
    type(rttov_coef_optpiclb), intent(inout) :: buffer
    integer,                   intent(in)    :: source
    integer,         optional, intent(in)    :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir container across all available processors
    !-----------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_coef_optpiclb@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_coef_optpiclb


  subroutine p_bcast_rttov_cnt_phasefn_int(buffer,source,comm)
    type(rttov_phasefn_int),  intent(inout) :: buffer
    integer,                  intent(in)    :: source
    integer,         optional,intent(in)    :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir container across all available processors
    !-----------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_phasefn_int@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_phasefn_int


  subroutine p_bcast_rttov_cnt_fast_coef(buffer,source,comm)
    type(rttov_fast_coef),  intent(inout) :: buffer
    integer,                intent(in)    :: source
    integer,       optional,intent(in)    :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir container across all available processors
    !-----------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_fast_coef@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_fast_coef


  subroutine p_bcast_rttov_cnt_fast_coef_gas(buffer,source,comm)
    type(rttov_fast_coef_gas),  intent(inout) :: buffer
    integer,                    intent(in)    :: source
    integer,           optional,intent(in)    :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir container across all available processors
    !-----------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_fast_coef_gas@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_fast_coef_gas


  subroutine p_bcast_rttov_cnt_nlte_coef(buffer,source,comm)
    type(rttov_nlte_coef),  intent(inout) :: buffer
    integer,                intent(in)    :: source
    integer,       optional,intent(in)    :: comm
    !-----------------------------------------------------------------------
    ! Broadcast an rttov_optpar_ir container across all available processors
    !-----------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_nlte_coef@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_nlte_coef


  subroutine p_bcast_rttov_cnt_coef_pccomp(buffer,source,comm)
    type(rttov_coef_pccomp), intent(inout) :: buffer
    integer,          intent(in)           :: source
    integer, optional,intent(in)           :: comm
    !-------------------------------------------------------------------------
    ! Broadcast an rttov_coef_pccomp container across all available processors
    !-------------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_coef_pccomp@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_coef_pccomp


  subroutine p_bcast_rttov_cnt_coef_pccomp1(buffer,source,comm)
    type(rttov_coef_pccomp1), intent(inout) :: buffer
    integer,          intent(in)            :: source
    integer, optional,intent(in)            :: comm
    !--------------------------------------------------------------------------
    ! Broadcast an rttov_coef_pccomp1 container across all available processors
    !--------------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm

    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_coef_pccomp1@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_coef_pccomp1


  subroutine p_bcast_rttov_cnt_coef_pccomp2(buffer,source,comm)
    type(rttov_coef_pccomp2), intent(inout) :: buffer
    integer,          intent(in)            :: source
    integer, optional,intent(in)            :: comm
    !--------------------------------------------------------------------------
    ! Broadcast an rttov_coef_pccomp2 container across all available processors
    !--------------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_coef_pccomp2@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_coef_pccomp2


#if defined(_RTTOV_ATLAS)
  subroutine p_bcast_rttov_atlas (atlas, source, comm)
    type(rttov_emis_atlas_data), intent(inout) :: atlas
    integer,                     intent(in)    :: source
    integer,           optional, intent(in)    :: comm
    !-----------------------------------------------------------------------------
    ! Broadcast an rttov_emis_atlas_data structure across all available processors
    !-----------------------------------------------------------------------------
    character(len=19), parameter :: proc = 'p_bcast_rttov_atlas'

    call p_bcast_rttov_container(atlas, source, comm)
    if (atlas% is_mw) then
      if (atlas% atlas_id == 1) then
        call p_bcast(atlas% telsem2_atlas, source, comm)
      else
        call p_bcast(atlas% cnrm_mw_atlas, source, comm)
      endif
    else
      call finish(proc,'broadcasting of IR atlases not implemented so far.')
    end if
  end subroutine p_bcast_rttov_atlas


  subroutine p_bcast_rttov_telsem(telsem, source, comm)
    type(telsem2_atlas_data), intent(inout) :: telsem
    integer,                  intent(in)    :: source
    integer, optional,        intent(in)    :: comm
    !--------------------------------------------------------------------------
    ! Broadcast a telsem2_atlas_data structure across all available processors
    !--------------------------------------------------------------------------
    integer :: dimensions(13)

    call p_bcast_rttov_container (telsem, source, comm)

    dimensions(:) = -1

    if (associated(telsem% ncells))         dimensions(1:1)   = shape(telsem% ncells)
    if (associated(telsem% firstcell))      dimensions(2:2)   = shape(telsem% firstcell)
    if (associated(telsem% emis))           dimensions(3:4)   = shape(telsem% emis)
    if (associated(telsem% correl))         dimensions(5:7)   = shape(telsem% correl)
    if (associated(telsem% emis_err))       dimensions(8:9)   = shape(telsem% emis_err)
    if (associated(telsem% class1))         dimensions(10:10) = shape(telsem% class1)
    if (associated(telsem% class2))         dimensions(11:11) = shape(telsem% class2)
    if (associated(telsem% cellnum))        dimensions(12:12) = shape(telsem% cellnum)
    if (associated(telsem% correspondance)) dimensions(13:13) = shape(telsem% correspondance)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(telsem% ncells))         allocate(telsem% ncells    (dimensions(1)))
      if (associated(telsem% firstcell))      allocate(telsem% firstcell (dimensions(2)))
      if (associated(telsem% emis))           allocate(telsem% emis      (dimensions(3),dimensions(4)))
      if (associated(telsem% correl))         allocate(telsem% correl    (dimensions(5),dimensions(6),dimensions(7)))
      if (associated(telsem% emis_err))       allocate(telsem% emis_err  (dimensions(8),dimensions(9)))
      if (associated(telsem% class1))         allocate(telsem% class1    (dimensions(10)))
      if (associated(telsem% class2))         allocate(telsem% class2    (dimensions(11)))
      if (associated(telsem% cellnum))        allocate(telsem% cellnum   (dimensions(12)))
      if (associated(telsem% correspondance)) allocate(telsem% correspondance (dimensions(13)))
    endif

    if (associated(telsem% ncells))          call p_bcast(telsem% ncells,source,comm)
    if (associated(telsem% firstcell))       call p_bcast(telsem% firstcell,source,comm)
    if (associated(telsem% emis))            call p_bcast(telsem% emis,source,comm)
    if (associated(telsem% correl))          call p_bcast(telsem% correl,source,comm)
    if (associated(telsem% emis_err))        call p_bcast(telsem% emis_err,source,comm)
    if (associated(telsem% class1))          call p_bcast(telsem% class1,source,comm)
    if (associated(telsem% class2))          call p_bcast(telsem% class2,source,comm)
    if (associated(telsem% cellnum))         call p_bcast(telsem% cellnum,source,comm)
    if (associated(telsem% correspondance))  call p_bcast(telsem% correspondance,source,comm)

  end subroutine p_bcast_rttov_telsem


  subroutine p_bcast_rttov_cnrm(cnrm, source, comm)
    type(cnrm_mw_atlas_data), intent(inout) :: cnrm
    integer,                  intent(in)    :: source
    integer, optional,        intent(in)    :: comm
    !--------------------------------------------------------------------------
    ! Broadcast a cnrm2_atlas_data structure across all available processors
    !--------------------------------------------------------------------------
    integer :: dimensions(9)

    call p_bcast_rttov_container (cnrm, source, comm)

    dimensions(:) = -1

    if (associated(cnrm% emissivity))     dimensions(1:4)   = shape(cnrm% emissivity)
    if (associated(cnrm% frequencies))    dimensions(5:5)   = shape(cnrm% frequencies)
    if (associated(cnrm% polarisation))   dimensions(6:6)   = shape(cnrm% polarisation)
    if (associated(cnrm% freq_range_min)) dimensions(7:7)   = shape(cnrm% freq_range_min)
    if (associated(cnrm% freq_range_max)) dimensions(8:8)   = shape(cnrm% freq_range_max)
    if (associated(cnrm% zenith_angles))  dimensions(9:9)   = shape(cnrm% zenith_angles)

    call p_bcast(dimensions,source,comm)

    if (pe_ifc /= source) then
      if (associated(cnrm% emissivity))     allocate(cnrm% emissivity(dimensions(1),dimensions(2),dimensions(3),dimensions(4)))
      if (associated(cnrm% frequencies))    allocate(cnrm% frequencies(dimensions(5)))
      if (associated(cnrm% polarisation))   allocate(cnrm% polarisation(dimensions(6)))
      if (associated(cnrm% freq_range_min)) allocate(cnrm% freq_range_min(dimensions(7)))
      if (associated(cnrm% freq_range_max)) allocate(cnrm% freq_range_max(dimensions(8)))
      if (associated(cnrm% zenith_angles))  allocate(cnrm% zenith_angles(dimensions(9)))
    endif

    if (associated(cnrm% emissivity))      call p_bcast(cnrm% emissivity,source,comm)
    if (associated(cnrm% frequencies))     call p_bcast(cnrm% frequencies,source,comm)
    if (associated(cnrm% polarisation))    call p_bcast(cnrm% polarisation,source,comm)
    if (associated(cnrm% freq_range_min))  call p_bcast(cnrm% freq_range_min,source,comm)
    if (associated(cnrm% freq_range_max))  call p_bcast(cnrm% freq_range_max,source,comm)
    if (associated(cnrm% zenith_angles))   call p_bcast(cnrm% zenith_angles,source,comm)
  end subroutine p_bcast_rttov_cnrm


  subroutine p_bcast_rttov_cnt_emis_atlas(buffer,source,comm)
    type(rttov_emis_atlas_data), intent(inout) :: buffer
    integer,                     intent(in)    :: source
    integer,            optional,intent(in)    :: comm
    !-----------------------------------------------------------------------------
    ! Broadcast an rttov_emis_atlas_data container across all available processors
    !-----------------------------------------------------------------------------
    integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
    call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
    lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
    call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
                   source, lcom, errorcode)
    if (errorcode /= MPI_SUCCESS) &
         call finish('p_bcast_rttov_cnt_emis_atlas@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
  end subroutine p_bcast_rttov_cnt_emis_atlas


 subroutine p_bcast_rttov_cnt_telsem(buffer,source,comm)
   type(telsem2_atlas_data), intent(inout)   :: buffer
   integer,                  intent(in)      :: source
   integer,         optional,intent(in)      :: comm
   !------------------------------------------------------------------
   ! Broadcast an rttov_coef container across all available processors
   !------------------------------------------------------------------
   integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
   call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
   lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
   call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)
   if (errorcode /= MPI_SUCCESS) &
        call finish('p_bcast_rttov_cnt_telsem@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
 end subroutine p_bcast_rttov_cnt_telsem


 subroutine p_bcast_rttov_cnt_cnrm(buffer,source,comm)
   type(cnrm_mw_atlas_data), intent(inout)   :: buffer
   integer,                  intent(in)      :: source
   integer,         optional,intent(in)      :: comm
   !------------------------------------------------------------------
   ! Broadcast an rttov_coef container across all available processors
   !------------------------------------------------------------------
   integer :: lcom, errorcode

#if defined(_RTIFC_USE_MPI_DACE)
   call p_bcast_derivedtype(buffer,size(transfer(buffer,(/' '/))),source,comm)
#else
   lcom = MPI_COMM_WORLD ;if (present (comm)) lcom = comm
   call MPI_Bcast(buffer,size(transfer(buffer,(/' '/))), MPI_BYTE, &
        source, lcom, errorcode)
   if (errorcode /= MPI_SUCCESS) &
        call finish('p_bcast_rttov_cnt_cnrm@mo_rtifc_13', 'MPI ERROR in MPI_Bcast')
#endif
 end subroutine p_bcast_rttov_cnt_cnrm

#endif /* _RTTOV_ATLAS */

#endif /* _RTIFC_DISTRIBCOEF */


  subroutine rtifc_print_profiles(unit)
    integer, intent(in) :: unit
    integer             :: j
    do j = 1, size(profiles)
      call print_profile(j,unit=unit)
    end do
  end subroutine rtifc_print_profiles


  subroutine print_profile(i, unit, form)
    integer, intent(in)           :: i
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: form
    !--------------------------------------------
    ! subroutine print_profile for debugging only
    !--------------------------------------------

    real(jprb), parameter :: RDRD = mh2o/mair
    real(jprb) :: q, ppmv_dry
    integer :: iunit, j
    logical :: form_loc

    if (present(unit)) then
      iunit=unit
    else
      iunit=stdout
    end if
    form_loc=.true.
    if (present(form)) form_loc = form

    if (form_loc) then
      if (i > size(profiles)) then
        write(iunit,*) 'rttov input profile: invalid index'
        return
      end if
      write(iunit,*) 'rttov input profile ',i,trim(profiles(i)%id)
      write(iunit,*) profiles(i)%nlevels
      if (associated(profiles(i)%p) .and. associated(profiles(i)%t) .and. associated(profiles(i)%q)) then
        write(iunit,'(3x,9x,"p",9x,"t",13x,"q",4x,"ppmv(dry)")')
        do j = 1, size(profiles(i)%p)
          select case (profiles(i)%gas_units)
          case(gas_unit_specconc)
            q = profiles(i)%q(j)
            ppmv_dry = (1.d6 * q) / (RDRD * (1 - q))
          case(gas_unit_ppmvdry)
            ppmv_dry = profiles(i)%q(j)
            q = (RDRD * ppmv_dry) / (1.d6 + RDRD * ppmv_dry)
          case(gas_unit_ppmv)
            ppmv_dry = profiles(i)%q(j) / (1._jprb - 1.E-6_jprb * profiles(i)%q(j))
            q = (RDRD * ppmv_dry) / (1.d6 + RDRD * ppmv_dry)
          end select
          write(iunit,'(I3,2(1x,f9.4),2(1x,E13.6))') j, profiles(i)%p(j), profiles(i)%t(j), q, ppmv_dry
        end do
      end if
      write(iunit,*) 'cloud'
      !    if (associated(profiles(i)%cloud)) write(iunit,'(6(1x,e25.17))'),profiles(i)%cloud(:,:)
      if (associated(profiles(i)%cloud)) write(iunit,*) profiles(i)%cloud(:,:)
      write(iunit,*) 'cfrac'
      !    if (associated(profiles(i)%cfrac)) write(iunit,'(6(1x,e25.17))'),profiles(i)%cfrac(:,:)
      if (associated(profiles(i)%cfrac)) write(iunit,*) profiles(i)%cfrac(:)
      write(iunit,*) 'idg=',profiles(i)%idg
      write(iunit,*) 'ice_scheme=',profiles(i)%ice_scheme
      write(iunit,*) 'skin'
      write(iunit,*) 'surftype=',profiles(i)%skin%surftype
      write(iunit,*) 'watertype=',profiles(i)%skin%watertype
      write(iunit,*) 't=',profiles(i)%skin%t
      write(iunit,*) 'fastem=',profiles(i)%skin%fastem
      write(iunit,*) '2m'
      write(iunit,*) 't=',profiles(i)%s2m%t
      write(iunit,*) 'q=',profiles(i)%s2m%q
      write(iunit,*) 'o=',profiles(i)%s2m%o
      write(iunit,*) 'p=',profiles(i)%s2m%p
      write(iunit,*) 'u=',profiles(i)%s2m%u
      write(iunit,*) 'v=',profiles(i)%s2m%v
      write(iunit,*) 'wfetch =',profiles(i)%s2m%wfetc
      write(iunit,*) 'zenangle=',profiles(i)%zenangle
      write(iunit,*) 'azangle=',profiles(i)%azangle
      write(iunit,*) 'sunzenangle=',profiles(i)%sunzenangle
      write(iunit,*) 'sunazangle=',profiles(i)%sunazangle
      write(iunit,*) 'elevation=',profiles(i)%elevation
      write(iunit,*) 'latitude=',profiles(i)%latitude
      write(iunit,*) 'longitude=',profiles(i)%longitude
      write(iunit,*) 'Be=',profiles(i)%Be
      write(iunit,*) 'cosbk=',profiles(i)%cosbk
      write(iunit,*) 'ctp=',profiles(i)%ctp
      write(iunit,*) 'cfraction=',profiles(i)%cfraction
    else
      if (i > size(profiles)) then
        write(iunit) 'rttov input profile: invalid index'
        return
      end if
      write(iunit) profiles(i)%nlevels
      if (associated(profiles(i)%p) .and. associated(profiles(i)%t) .and. associated(profiles(i)%q)) then
        do j = 1, size(profiles(i)%p)
          write(iunit) j, profiles(i)%p(j), profiles(i)%t(j), profiles(i)%q(j)
        end do
      end if
      if (associated(profiles(i)%cloud)) write(iunit) profiles(i)%cloud(:,:)
      if (associated(profiles(i)%cfrac)) write(iunit) profiles(i)%cfrac(:)
      write(iunit) profiles(i)%idg
      write(iunit) profiles(i)%ice_scheme
      write(iunit) profiles(i)%skin%surftype
      write(iunit) profiles(i)%skin%watertype
      write(iunit) profiles(i)%skin%t
      write(iunit) profiles(i)%skin%fastem
      write(iunit) profiles(i)%s2m%t
      write(iunit) profiles(i)%s2m%q
      write(iunit) profiles(i)%s2m%o
      write(iunit) profiles(i)%s2m%p
      write(iunit) profiles(i)%s2m%u
      write(iunit) profiles(i)%s2m%v
      write(iunit) profiles(i)%s2m%wfetc
      write(iunit) profiles(i)%zenangle
      write(iunit) profiles(i)%azangle
      write(iunit) profiles(i)%sunzenangle
      write(iunit) profiles(i)%sunazangle
      write(iunit) profiles(i)%elevation
      write(iunit) profiles(i)%latitude
      write(iunit) profiles(i)%longitude
      write(iunit) profiles(i)%Be
      write(iunit) profiles(i)%cosbk
      write(iunit) profiles(i)%ctp
      write(iunit) profiles(i)%cfraction
    end if
  end subroutine print_profile

!==================================
#endif /* (_RTTOV_VERSION == 12) */
!==================================

end module mo_rtifc_12
