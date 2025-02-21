! Interface for the RTTOV library (version 10 and later)
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

MODULE mo_rtifc

!-------------------------------------------------------------------------------------
!
! Description:
!   This module and the associated mo_rtifc_*.f90 modules contain subroutines
!   to work with RTTOV. The purpose of this module is to make work with
!   different RTTOV versions more user friendly. The interfaces and options of
!   this module should not change from one RTTOV version to the next - in
!   contrast to the original RTTOV code.
!   The following is provided:
!    - routines to modify/check/print RTTOV configuration/options
!    - an initialization routine which initializes the RTTOV modules
!      and reads the required instrument specific coefficients.
!      If _RTIFC_DISTRIBCOEF is set during compilation, this routine has the
!      option to read the coefficients on one PE and distribute them to the
!      others.
!    - Routines to fill the RTTOV profiles structure.
!    - Routines to call rttov in forward and k mode
!    - Some auxiliary routines
!   This module (mo_rtifc.f90) is just a wrapper module, that uses routines,
!   variables and parameters from associated modules (depending on the Macros that
!   were set during compilation).
!   Basic routines and stuff, that should not change from one RTTOV version to
!   another (e.g. DWD specific stuff) is located in mo_rtifc_base.f90.
!   All stuff that depends somehow on the RTTOV version is located in version
!   specific modules mo_rtifc_${rttov_version}.f90 (rttov_version=nort,10,12,13...).
!
!-------------------------------------------------------------------------------------

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

  use mo_rtifc_base

  use mo_rtifc_nort,     only: rtifc_vers,             &! RTTOV version this interface is compiled for
                               rtifc_version,          &! version string
                               rtifc_init,             &! Initialize RTTOV, i.e. read coeffs
                               rtifc_coef_index,       &! Return index of coeffs for given satid/instr
                               rtifc_cleanup,          &! Cleanup RTTOV and mo_rtifc structures
                               rtifc_fill_input,       &! Fill rttov profile structure
                               rtifc_direct,           &! rttov_direct call
                               rtifc_k,                &! rttov_k call
                               rtifc_l2c_god,          &! god-corrected l2c
                               rtifc_print_profiles,   &! print rttov profile structure
                               rtifc_coef_prop,        &! get coef properties
                               rttov_options,          &! Options for RTTOV (type definition)
                               rt_gas_id_o3  => gas_id_ozone, &!
                               rt_gas_id_co2 => gas_id_co2


  implicit none

  !-------------
  ! Public stuff
  !-------------

  private

  ! RTTOV version
  public :: rtifc_vers           ! the rttov version this interface was compiled for

  ! subroutines
  public :: rtifc_check_config     ! Check RTTOV version and level number, set nlevs_top
  public :: rtifc_version          ! RTTOV and interface version string
  public :: rtifc_set_opts         ! set RTTOV options
  public :: rtifc_get_opts         ! get RTTOV options
  public :: rtifc_init             ! Initialise RTTOV modules, read coeffs
  public :: rtifc_coef_index       ! Returns index of coeffs for given satid/instr
  public :: rtifc_cleanup          ! frees memory allocated by rtifc_init
  public :: rtifc_fill_input       ! fills the profile-dependent part for RTTOV
  public :: rtifc_direct           ! calls RTTOV direct routine
  public :: rtifc_k                ! calls RTTOV K routine
  public :: rtifc_l2c_god          ! god-corrected l2c
  public :: rtifc_print_profiles   ! print profiles
  public :: rtifc_errmsg           ! gives error message corresponding to exit status
  public :: rtifc_coef_prop        ! get properties of coefs
  public :: rtifc_check_nlevs

  ! Optimization
  public :: read1pe
  public :: rtifc_alloc_mode

  ! error codes/messages
  public :: NO_ERROR             ! everything was ok.
  public :: WARN_RTTOV           ! warning

  ! options
  public :: rttov_options        ! type definition

  ! output flags
  public :: OUT_ASB
  public :: OUT_CSB
  public :: OUT_ASR
  public :: OUT_CSR
  public :: OUT_VIS

  ! RTTOV levels above user levels
  public :: nlevs_top

  ! default profile values
  public :: default_wfetch
  public :: default_fastem
  public :: default_watertype
  public :: default_salinity
  public :: default_o3_surf
  public :: default_satazim
  public :: default_sunzenangle
  public :: default_sunazangle
  public :: default_ctp
  public :: default_cfraction
  public :: default_idg
  public :: default_clw_scheme
  public :: default_gas_units

  ! hard limits on profile variables
  public :: qmin_ifc
  public :: qmax_ifc
  public :: tmin_ifc
  public :: tmax_ifc

  ! RTTOV "constants"
  public :: min_od
  public :: rt_gas_id_o3
  public :: rt_gas_id_co2

  ! Regularization limits
  public :: chk_reg_lims
  public :: chk_plim_t
  public :: chk_plim_q

  ! god (generalized optical depth) parameters
  public :: god_par_file
  public :: wr_god
  public :: out_path
  public :: chk_god
  public :: god_thresh

  ! atlas use
  public :: atlas_single_inst

  ! surftypes
  public :: rts_land, rts_sea, rts_ice, rts_name

  interface rtifc_check_nlevs
    module procedure check_nlevs
  end interface rtifc_check_nlevs

contains


  subroutine rtifc_set_opts(iopt,               &!
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
                            use_t2m_opdep,      &!
                            use_q2m,            &!
                            do_lambertian,      &!
                            cloud_overlap,      &!
                            do_checkinput,      &!
                            ozone_data,         &!
                            co2_data,           &!
                            n2o_data,           &!
                            co_data,            &!
                            ch4_data,           &!
                            so2_data,           &!
                            clw_data,           &!
                            dom_rayleigh,       &!
                            dom_nstreams,       &!
                            ir_scatt_model,     &!
                            vis_scatt_model,    &!
                            clip_gas_opdep      &!
                           )
    integer,             intent(inout), optional :: iopt  ! ID of options, index in rt_opts
    logical,             intent(in),    optional :: new   ! Create new options, i.e. new entry
                                                          ! in rt_opts
    logical,             intent(in),    optional :: init  ! initialization
    character(len=*),    intent(in),    optional :: tmpl       ! name of template
    integer,             intent(in),    optional :: iopt_tmpl  ! index of "source" options
    character(len=*),    intent(in),    optional :: name       ! meta data: name of options
    integer,             intent(in),    optional :: satid      ! meta data: satellite ID
    integer,             intent(in),    optional :: instr      ! meta data: instrument
    integer,             intent(in),    optional :: grid       ! meta data: grid
    type(rttov_options), intent(in),    optional :: rttov_opts ! full rttov options
    logical,             intent(in),    optional :: addinterp
    integer,             intent(in),    optional :: interp_mode
    logical,             intent(in),    optional :: addrefrac
    logical,             intent(in),    optional :: addclouds
    logical,             intent(in),    optional :: addaerosl
    logical,             intent(in),    optional :: addsolar
    logical,             intent(in),    optional :: addpc
    logical,             intent(in),    optional :: apply_reg_lims
    logical,             intent(in),    optional :: verbose_reg_lims
    logical,             intent(in),    optional :: crop_k_reg_lims
    logical,             intent(in),    optional :: switchrad
    logical,             intent(in),    optional :: conv_overc
    integer,             intent(in),    optional :: fix_hgpl
    integer,             intent(in),    optional :: fastem_version
    integer,             intent(in),    optional :: ir_sea_emis_model
    logical,             intent(in),    optional :: use_t2m_opdep
    logical,             intent(in),    optional :: use_q2m
    logical,             intent(in),    optional :: do_lambertian
    integer,             intent(in),    optional :: cloud_overlap
    logical,             intent(in),    optional :: do_checkinput
    logical,             intent(in),    optional :: ozone_data
    logical,             intent(in),    optional :: co2_data
    logical,             intent(in),    optional :: n2o_data
    logical,             intent(in),    optional :: co_data
    logical,             intent(in),    optional :: ch4_data
    logical,             intent(in),    optional :: so2_data
    logical,             intent(in),    optional :: clw_data
    logical,             intent(in),    optional :: dom_rayleigh
    integer,             intent(in),    optional :: dom_nstreams
    integer,             intent(in),    optional :: ir_scatt_model
    integer,             intent(in),    optional :: vis_scatt_model
    logical,             intent(in),    optional :: clip_gas_opdep

    character(len=14),   parameter :: proc = 'rtifc_set_opts'


  end subroutine rtifc_set_opts


  subroutine rtifc_get_opts(iopt,               &!
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
                            use_t2m_opdep,      &!
                            use_q2m,            &!
                            do_lambertian,      &!
                            cloud_overlap,      &!
                            do_checkinput,      &!
                            ozone_data,         &!
                            co2_data,           &!
                            n2o_data,           &!
                            co_data,            &!
                            ch4_data,           &!
                            so2_data,           &!
                            clw_data,           &!
                            dom_rayleigh,       &!
                            dom_nstreams,       &!
                            ir_scatt_model,     &!
                            vis_scatt_model,    &!
                            clip_gas_opdep      &!
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
    logical,             intent(out),    optional :: use_t2m_opdep
    logical,             intent(out),    optional :: use_q2m
    logical,             intent(out),    optional :: do_lambertian
    integer,             intent(out),    optional :: cloud_overlap
    logical,             intent(out),    optional :: do_checkinput
    logical,             intent(out),    optional :: ozone_data
    logical,             intent(out),    optional :: co2_data
    logical,             intent(out),    optional :: n2o_data
    logical,             intent(out),    optional :: co_data
    logical,             intent(out),    optional :: ch4_data
    logical,             intent(out),    optional :: so2_data
    logical,             intent(out),    optional :: clw_data
    logical,             intent(out),    optional :: dom_rayleigh
    integer,             intent(out),    optional :: dom_nstreams
    integer,             intent(out),    optional :: ir_scatt_model
    integer,             intent(out),    optional :: vis_scatt_model
    logical,             intent(out),    optional :: clip_gas_opdep

    character(len=14),   parameter :: proc = 'rtifc_get_opts'


  end subroutine rtifc_get_opts


  function rtifc_errmsg(code, v) result(msg)
    character(len=120)            :: msg
    integer, intent(in)           :: code
    logical, intent(in), optional :: v

    logical :: v_

    msg = trim(errmsg(code))

    if (present(v)) then
      v_ = v
    else
      v_ = .false.
    end if
    if (v_) then
      msg = rtifc_version()//' '//trim(msg)
    end if

  end function rtifc_errmsg

end module mo_rtifc
