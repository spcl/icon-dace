!
! Classes and functions for the turbulent mixing package (tmx)
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_vdf_atmo

  USE mo_kind,              ONLY: wp, sp, i1
  USE mo_exception,         ONLY: message, finish
  USE mtime,                ONLY: t_datetime => datetime
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain
  USE mo_variable,          ONLY: t_variable
  USE mo_variable_list,     ONLY: t_variable_list, t_variable_set
  USE mo_util_string,       ONLY: real2string
  USE mo_model_domain,      ONLY: t_patch
  USE mo_intp_data_strc,    ONLY: t_int_state, p_int_state
  USE mo_nonhydro_types,    ONLY: t_nh_metrics
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_aes_sfc_indices,   ONLY: nsfc_type, iwtr, iice, ilnd
  USE mo_physical_constants,ONLY: grav, rd, cpd, cpv, cvd, rd_o_cpd, &
    &                             vtmpc1, p0ref, rgrav, alvdcp, alv
  USE mo_aes_thermo,        ONLY: sat_pres_water, specific_humidity
  USE mo_turb_vdiff_params, ONLY: ckap
  USE mo_math_constants,    ONLY: pi_2, ln2
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_nh_vert_interp_les,ONLY: brunt_vaisala_freq, vert_intp_full2half_cell_3d
  USE mo_intp,              ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array,     &
  &                             sync_patch_array_mult
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_intp_rbf,          ONLY: rbf_vec_interpol_vertex, rbf_vec_interpol_edge
  USE mo_fortran_tools,     ONLY: init

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vdf_atmo, t_vdf_atmo_inputs, t_vdf_atmo_config, t_vdf_atmo_diagnostics, test !, &
    ! & compute_temp_from_static_energy

  !Parameters for surface layer parameterizations: From Zeng_etal 1997 J. Clim
  REAL(wp), PARAMETER :: bsm = 5.0_wp  !Businger Stable Momentum
  REAL(wp), PARAMETER :: bum = 16._wp  !Businger Untable Momentum
  REAL(wp), PARAMETER :: bsh = 5.0_wp  !Businger Stable Heat
  REAL(wp), PARAMETER :: buh = 16._wp  !Businger Untable Heat


  TYPE, EXTENDS(t_tmx_process) :: t_vdf_atmo
    ! TYPE(t_vdf_atmo_inputs) :: inputs
    ! PROCEDURE(i_temp_to_energy), POINTER :: temp_to_energy => NULL()
  CONTAINS
    PROCEDURE :: Compute
    PROCEDURE :: Compute_diagnostics
    ! PROCEDURE(i_temp_to_energy) :: temp_to_energy
    PROCEDURE :: temp_to_energy
    PROCEDURE :: energy_to_temp
    PROCEDURE :: compute_flux_x
    PROCEDURE :: Update_diagnostics
  END TYPE t_vdf_atmo
  
  ! ABSTRACT INTERFACE
  !   SUBROUTINE i_temp_to_energy(this, configs, inputs, temperature, energy)
  !     IMPORT t_vdf_atmo, t_variable_set
  !     CLASS(t_vdf_atmo), INTENT(in) :: this
  !     CLASS(t_variable_set), INTENT(in) :: configs
  !     CLASS(t_variable_set), INTENT(in) :: inputs
  !     REAL(wp), INTENT(in) :: temperature(:,:,:)
  !     REAL(wp), INTENT(out) :: energy(:,:,:)
  !   END SUBROUTINE
  ! END INTERFACE

  INTERFACE t_vdf_atmo
    MODULE PROCEDURE t_vdf_atmo_construct
  END INTERFACE

  TYPE, EXTENDS(t_variable_set) :: t_vdf_atmo_inputs
    REAL(wp), POINTER :: &
      ! exchange coefficient
      & ptvm1(:,:,:) => NULL() , &
      & rho(:,:,:) => NULL() , &
      & mair(:,:,:) => NULL() , &
      & cvair(:,:,:) => NULL() , &
      & zf(:,:,:) => NULL() , &
      & zh(:,:,:) => NULL() , &
      ! exchange coefficient + boundary condition
      & ptm1(:,:,:) => NULL() , &
      & pum1(:,:,:) => NULL() , &
      & pvm1(:,:,:) => NULL() , &
      & pwp1(:,:,:) => NULL() , &
      & pqm1(:,:,:) => NULL() , &
      & pxlm1(:,:,:) => NULL() , &
      & pxim1(:,:,:) => NULL() , &
      & pxrm1(:,:,:) => NULL() , &
      & pxsm1(:,:,:) => NULL() , &
      & pxgm1(:,:,:) => NULL() , &
      & papm1(:,:,:) => NULL() , &
      & paphm1(:,:,:) => NULL() , &
      & dz(:,:,:) => NULL() , &
      & inv_dzh(:,:,:) => NULL(), &
      & pfrc(:,:,:) => NULL()
#ifdef __MIXED_PRECISION
    REAL(sp), POINTER :: &
      & inv_dzf(:,:,:) => NULL()
#else
    REAL(wp), POINTER :: &
      & inv_dzf(:,:,:) => NULL()
#endif
  CONTAINS
    ! PROCEDURE :: Init => init_t_vdf_atmo_variable_set
    PROCEDURE :: Set_pointers => Set_pointers_inputs
  END TYPE t_vdf_atmo_inputs

  TYPE, EXTENDS(t_variable_set) :: t_vdf_atmo_config
    REAL(wp), POINTER :: &
      cpd => NULL(), &
      cvd => NULL(), &
      dissipation_factor => NULL(), &
      rturb_prandtl => NULL(), &
      turb_prandtl  => NULL(), &
      louis_constant_b => NULL(), &
      km_min => NULL(), &
      dtime => NULL(),  &
      k_s => NULL()
    INTEGER, POINTER :: &
      solver_type => NULL(), &
      energy_type => NULL()
    LOGICAL, POINTER :: &
      use_louis  => NULL()
    CONTAINS
    ! PROCEDURE :: Init => init_t_vdf_atmo_variable_set
    PROCEDURE :: Set_pointers => Set_pointers_config
  END TYPE t_vdf_atmo_config

  TYPE, EXTENDS(t_variable_set) :: t_vdf_atmo_diagnostics
    REAL(wp), POINTER :: &
      & ghf(:,:,:) => NULL(), &
      & ctgz(:,:,:) => NULL(), &
      & ctgzvi(:,:) => NULL(), &
      & div_c(:,:,:) => NULL(), &
      & theta_v(:,:,:) => NULL(), &
      & pprfac(:,:,:) => NULL(), &
      & rho_ic(:,:,:) => NULL(), &
      & bruvais(:,:,:) => NULL(), &
      & stability_function(:,:,:) => NULL(), &
      & vn_ie(:,:,:) => NULL(), &
      & vt_ie(:,:,:) => NULL(), &
      & w_ie(:,:,:)  => NULL(), &
      & km_ic(:,:,:) => NULL(), &
      & kh_ic(:,:,:) => NULL(), &
      & km_c(:,:,:) => NULL(), &
      & km_ie(:,:,:) => NULL(), &
      & km_iv(:,:,:) => NULL(), &
      & km(:,:,:)    => NULL(), &
      & kh(:,:,:)    => NULL(), &
      & vn(:,:,:)    => NULL(), &
      & shear(:,:,:) => NULL(), &
      & div_of_stress(:,:,:) => NULL(), &
      & mech_prod(:,:,:) => NULL(), &
      & u_vert(:,:,:) => NULL(), &
      & v_vert(:,:,:) => NULL(), &
      & w_vert(:,:,:) => NULL(), &
      !
      & dissip_kin_energy(:,:,:) => NULL(), &
      & dissip_kin_energy_vi(:,:) => NULL(), &
      & heating(:,:,:) => NULL(), &
      & scaling_factor_louis(:,:) => NULL(), &
      & internal_energy_vi(:,:) => NULL(), &
      & internal_energy_vi_tend(:,:) => NULL()
      ! boundary condition
      ! & pcfm_tile(:,:,:) => NULL(), &
      ! & pcfh_tile(:,:,:) => NULL(), &
      ! & pch_tile(:,:,:)  => NULL(), &
      ! & pbn_tile(:,:,:)  => NULL(), &
      ! & pbhn_tile(:,:,:) => NULL(), &
      ! & pbm_tile(:,:,:)  => NULL(), &
      ! & pbh_tile(:,:,:)  => NULL(), &
      ! & pcpt_tile(:,:,:) => NULL(), &
      ! & pqsat_tile(:,:,:) => NULL()
CONTAINS
    ! PROCEDURE :: Init => init_t_vdf_atmo_variable_set
    PROCEDURE :: Set_pointers => Set_pointers_diagnostics
  END TYPE t_vdf_atmo_diagnostics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf_atmo'

CONTAINS

  FUNCTION t_vdf_atmo_construct(name, dt, domain) RESULT(result)

    CHARACTER(len=*), INTENT(in) :: name
    REAL(wp),         INTENT(in) :: dt
    TYPE(t_domain),      POINTER :: domain
    TYPE(t_vdf_atmo),    POINTER :: result

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_vdf_atmo_construct'

    CALL message(routine, '')

    ALLOCATE(t_vdf_atmo::result)
    !$ACC ENTER DATA COPYIN(result)
    ! Call Init of abstract parent class
    CALL result%Init(dt=dt, name=name, domain=domain)
    __acc_attach(result%domain)

    ! Initialize variable sets
    ALLOCATE(t_vdf_atmo_config :: result%config)
    result%config%list = build_atmo_config_list(result%domain)
    !$ACC ENTER DATA COPYIN(result%config)
    ! CALL result%config%Init(build_atmo_config_list(result%domain))

    ALLOCATE(t_vdf_atmo_inputs :: result%inputs)
    result%inputs%list = build_atmo_input_list(result%domain)
    !$ACC ENTER DATA COPYIN(result%inputs)
    ! CALL result%inputs%Init(build_atmo_input_list(result%domain))

    ALLOCATE(t_vdf_atmo_diagnostics :: result%diagnostics)
    result%diagnostics%list = build_atmo_diagnostic_list(result%domain)
    !$ACC ENTER DATA COPYIN(result%diagnostics)
    ! CALL result%diagnostics%list%allocator()
    ! CALL result%diagnostics%Set_pointers()

  END FUNCTION t_vdf_atmo_construct

  SUBROUTINE Compute(this, datetime)

    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in),   POINTER :: datetime     !< date and time at beginning of time step

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute'

    CALL message(routine, '')

  END SUBROUTINE Compute

  SUBROUTINE Compute_diagnostics(this, datetime)

    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime

    TYPE(t_vdf_atmo_config),      POINTER :: conf
    TYPE(t_vdf_atmo_inputs),      POINTER :: ins
    TYPE(t_vdf_atmo_diagnostics), POINTER :: diags

    INTEGER :: jg, nlev, jsfc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute_diagnostics'


    jg = this%domain%patch%id
    jsfc = 1

    nlev = this%domain%nlev

    SELECT TYPE (set => this%config)
    TYPE IS (t_vdf_atmo_config)
      conf => set
    END SELECT
    __acc_attach(conf)
    SELECT TYPE (set => this%inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins => set
    END SELECT
    __acc_attach(ins)
    SELECT TYPE (set => this%diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags => set
    END SELECT
    __acc_attach(diags)

    CALL compute_geopotential_height_above_ground(this%domain, ins%zf(:,:,:), ins%zh(:,:,:) , diags%ghf(:,:,:))

    CALL compute_static_energy(this%domain, conf%cpd, ins%ptm1(:,:,:), diags%ghf(:,:,:), diags%ctgz(:,:,:))

    CALL compute_exchange_coefficient(this%domain, this%config, this%inputs, this%diagnostics)
    ! CALL compute_exchange_coefficient(this%domain, p_nh_state(jg)%metrics, p_int_state(jg))

  END SUBROUTINE Compute_diagnostics

  SUBROUTINE Update_diagnostics(this)

    CLASS(t_vdf_atmo), INTENT(inout), TARGET :: this

    TYPE(t_vdf_atmo_config),      POINTER :: conf
    TYPE(t_vdf_atmo_inputs),      POINTER :: ins
    TYPE(t_vdf_atmo_diagnostics), POINTER :: diags

    INTEGER :: jb, jk, jc
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & state_ta, state_qv, state_qc, state_qi, &
      & tend_ta, tend_qv, tend_qc, tend_qi, &
      & new_state_ta, new_state_qv, new_state_qc, new_state_qi

    REAL(wp), DIMENSION(this%domain%nproma,this%domain%nblks_c) :: &
      & internal_energy_vi_old

    CHARACTER(len=*), PARAMETER :: routine = modname//':Update_diagnostics'

    ! CALL message(routine, 'start')

    SELECT TYPE (set => this%config)
    TYPE IS (t_vdf_atmo_config)
      conf => set
    END SELECT
    __acc_attach(conf)
    SELECT TYPE (set => this%inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins => set
    END SELECT
    __acc_attach(ins)
    SELECT TYPE (set => this%diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags => set
    END SELECT
    __acc_attach(diags)

    state_ta     => this%states    %Get_ptr_r3d('temperature') ! old state
    tend_ta      => this%tendencies%Get_ptr_r3d('temperature') ! tendency
    new_state_ta => this%new_states%Get_ptr_r3d('temperature') ! new state
    state_qv     => this%states    %Get_ptr_r3d('water vapor') ! old state
    tend_qv      => this%tendencies%Get_ptr_r3d('water vapor') ! tendency
    new_state_qv => this%new_states%Get_ptr_r3d('water vapor') ! new state
    state_qc     => this%states    %Get_ptr_r3d('cloud water') ! old state
    tend_qc      => this%tendencies%Get_ptr_r3d('cloud water') ! tendency
    new_state_qc => this%new_states%Get_ptr_r3d('cloud water') ! new state
    state_qi     => this%states    %Get_ptr_r3d('cloud ice')   ! old state
    tend_qi      => this%tendencies%Get_ptr_r3d('cloud ice')   ! tendency
    new_state_qi => this%new_states%Get_ptr_r3d('cloud ice') ! new state

    ASSOCIATE( &
      domain => this%domain, &
      cpd => conf%cpd, &
      dtime => conf%dtime, &
      dz => ins%dz, &
      rho => ins%rho, &
      ctgz => diags%ctgz, &
      ctgzvi => diags%ctgzvi, &
      dissip_kin_energy => diags%dissip_kin_energy, &
      dissip_kin_energy_vi => diags%dissip_kin_energy_vi, &
      internal_energy_vi => diags%internal_energy_vi, &
      internal_energy_vi_tend => diags%internal_energy_vi_tend &
      & )

    !$ACC DATA CREATE(internal_energy_vi_old)

    CALL compute_static_energy( &
      & domain, &
      & cpd, new_state_ta(:,:,:), diags%ghf(:,:,:), &
      & ctgz(:,:,:) &
      & )

    ! Vertical integrals

    ! TODO: include hydrometeors from microphysics?
    CALL compute_internal_energy_vi( &
      & domain, &
      & rho(:,:,:), dz(:,:,:), ins%pxrm1(:,:,:), ins%pxsm1(:,:,:), ins%pxgm1(:,:,:), &
      & state_ta(:,:,:), state_qv(:,:,:), state_qc(:,:,:), state_qi(:,:,:), &
      & internal_energy_vi_old(:,:) &
      & )
    CALL compute_internal_energy_vi( &
      & domain, &
      & rho(:,:,:), dz(:,:,:), ins%pxrm1(:,:,:), ins%pxsm1(:,:,:), ins%pxgm1(:,:,:), &
      & new_state_ta(:,:,:), new_state_qv(:,:,:), new_state_qc(:,:,:), new_state_qi(:,:,:), &
      & internal_energy_vi(:,:) &
      & )

!$OMP PARALLEL
      CALL init(ctgzvi, lacc=.TRUE.)
      CALL init(dissip_kin_energy_vi, lacc=.TRUE.)
      CALL init(internal_energy_vi_tend, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb),domain%i_endidx_c(jb)

          ! static energy (reference air path)
          ctgzvi(jc,jb) = ctgzvi(jc,jb) + ctgz(jc,jk,jb) * rho(jc,jk,jb) * dz(jc,jk,jb)
          ! kinetic energy dissipation
          dissip_kin_energy_vi(jc,jb) = dissip_kin_energy_vi(jc,jb) + dissip_kin_energy(jc,jk,jb)

          ! internal energy tendency
          internal_energy_vi_tend(jc,jb) = (internal_energy_vi(jc,jb) - internal_energy_vi_old(jc,jb)) / dtime

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)
    !$ACC END DATA

    END ASSOCIATE

    ! CALL message(routine, 'end')

  END SUBROUTINE Update_diagnostics

  ! SUBROUTINE init_t_vdf_atmo_variable_set(this, varlist)

  !   CLASS(t_variable_set), INTENT(inout) :: this
  !   TYPE(t_variable_list), INTENT(in)    :: varlist

  !   INTEGER :: nproma, nblks_c, nlev, nlevp1
  !   INTEGER :: shape_2d(2), shape_3d(3)

  !   this%list = varlist

  ! END SUBROUTINE init_t_vdf_atmo_variable_set

  FUNCTION build_atmo_config_list(domain) RESULT(configlist)

    TYPE(t_variable_list) :: configlist
    TYPE(t_domain),        INTENT(in)    :: domain

    INTEGER :: nproma, nblks_c, nlev, nlevp1
    INTEGER :: shape_0d(1)

    configlist = t_variable_list('config')

    nproma  = domain%nproma
    nblks_c = domain%nblks_c
    nlev    = domain%nlev
    nlevp1  = nlev + 1

    shape_0d = [0]
    CALL configlist%append(t_variable('cpd', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('cvd', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('reverse prandtl number', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('minimum Km', shape_0d, "m2/s", type_id="real"))
    CALL configlist%append(t_variable('k_s', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('prandtl number', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('switch to activate Louis formula', shape_0d, "", type_id="logical"))
    CALL configlist%append(t_variable('Louis constant b', shape_0d, "", type_id="real"))
    CALL configlist%append(t_variable('time step', shape_0d, "s", type_id="real"))
    CALL configlist%append(t_variable('solver type', shape_0d, "", type_id="integer"))
    CALL configlist%append(t_variable('energy type', shape_0d, "", type_id="integer"))
    CALL configlist%append(t_variable('dissipation factor', shape_0d, "", type_id="real"))

  END FUNCTION build_atmo_config_list

  SUBROUTINE Set_pointers_config(this)

    CLASS(t_vdf_atmo_config), INTENT(inout) :: this

    SELECT TYPE (this)
    TYPE IS (t_vdf_atmo_config)
      this%cpd => this%list%Get_ptr_r0d('cpd')
      __acc_attach(this%cpd)
      this%cvd => this%list%Get_ptr_r0d('cvd')
      __acc_attach(this%cvd)
      this%rturb_prandtl => this%list%Get_ptr_r0d('reverse prandtl number')
      __acc_attach(this%rturb_prandtl)
      this%turb_prandtl  => this%list%Get_ptr_r0d('prandtl number')
      __acc_attach(this%turb_prandtl)
      this%use_louis  => this%list%Get_ptr_l0d('switch to activate Louis formula')
      __acc_attach(this%use_louis)
      this%louis_constant_b => this%list%Get_ptr_r0d('Louis constant b')
      __acc_attach(this%louis_constant_b)
      this%km_min => this%list%Get_ptr_r0d('minimum Km')
      __acc_attach(this%km_min)
      this%k_s => this%list%Get_ptr_r0d('k_s')
      __acc_attach(this%k_s)
      this%dtime => this%list%Get_ptr_r0d('time step')
      __acc_attach(this%dtime)
      this%solver_type => this%list%Get_ptr_i0d('solver type')
      __acc_attach(this%solver_type)
      this%energy_type => this%list%Get_ptr_i0d('energy type')
      __acc_attach(this%energy_type)
      this%dissipation_factor => this%list%Get_ptr_r0d('dissipation factor')
      __acc_attach(this%dissipation_factor)
    END SELECT

  END SUBROUTINE Set_pointers_config

  FUNCTION build_atmo_input_list(domain) RESULT(inlist)

    TYPE(t_variable_list) :: inlist
    TYPE(t_domain),        INTENT(in)    :: domain

    INTEGER :: nproma, nblks_c, nlev, nlevp1
    INTEGER :: shape_2d(2), shape_3d(3)

    inlist = t_variable_list('inputs')

    nproma  = domain%nproma
    nblks_c = domain%nblks_c
    nlev    = domain%nlev
    nlevp1  = nlev + 1

    shape_3d = [nproma,nlev,nblks_c]
    CALL inlist%append(t_variable('temperature',     shape_3d, "K", type_id="real"))
    CALL inlist%append(t_variable('virtual temperature',     shape_3d, "K", type_id="real"))
    CALL inlist%append(t_variable('air density',     shape_3d, "kg/m3", type_id="real"))
    CALL inlist%append(t_variable('zonal wind',      shape_3d, "m/s", type_id="real"))
    CALL inlist%append(t_variable('meridional wind', shape_3d, "m/s", type_id="real"))
    CALL inlist%append(t_variable('water vapor', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('cloud water', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('cloud ice', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('rain', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('snow', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('graupel', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('full level pressure', shape_3d, "Pa", type_id="real"))
    CALL inlist%append(t_variable('layer thickness', shape_3d, "m", type_id="real"))
#ifdef __MIXED_PRECISION
    CALL inlist%append(t_variable('inverse layer thickness full', shape_3d, "m", type_id="single"))
#else
    CALL inlist%append(t_variable('inverse layer thickness full', shape_3d, "m", type_id="real"))
#endif
    ! CALL inlist%append(t_variable('saturation specific humidity', shape_3d, "kg/kg", type_id="real"))
    CALL inlist%append(t_variable('moist air mass', shape_3d, "kg/m2", type_id="real"))
    ! CALL inlist%append(t_variable('specific heat of air at constant pressure', shape_3d, "J/kg/K", type_id="real"))
    CALL inlist%append(t_variable('specific heat of air at constant volume',   shape_3d, "J/kg/K", type_id="real"))
    ! CALL inlist%append(t_variable('conv. factor layer heating to temp. tendency', shape_3d, "(K/s)/(W/m2)", type_id="real"))
    CALL inlist%append(t_variable('geometric height full', shape_3d, "m", type_id="real"))

    shape_3d = [nproma,nlev+1,nblks_c]
    CALL inlist%append(t_variable('half level pressure', shape_3d, "Pa", type_id="real"))
    CALL inlist%append(t_variable('inverse layer thickness half', shape_3d, "m", type_id="real"))
    CALL inlist%append(t_variable('vertical wind', shape_3d, "m/s", type_id="real"))
    CALL inlist%append(t_variable('geometric height half', shape_3d, "m", type_id="real"))

  END FUNCTION build_atmo_input_list

  SUBROUTINE Set_pointers_inputs(this)

    CLASS(t_vdf_atmo_inputs), INTENT(inout) :: this

    SELECT TYPE (this)
    TYPE IS (t_vdf_atmo_inputs)
      this%ptm1    => this%list%Get_ptr_r3d('temperature')
      __acc_attach(this%ptm1)
      this%ptvm1   => this%list%Get_ptr_r3d('virtual temperature')
      __acc_attach(this%ptvm1)
      this%rho     => this%list%Get_ptr_r3d('air density')
      __acc_attach(this%rho)
      this%mair    => this%list%Get_ptr_r3d('moist air mass')
      this%cvair   => this%list%Get_ptr_r3d('specific heat of air at constant volume')
      __acc_attach(this%cvair)
      this%zf      => this%list%Get_ptr_r3d('geometric height full')
      __acc_attach(this%zf)
      this%zh      => this%list%Get_ptr_r3d('geometric height half')
      __acc_attach(this%zh)
      this%pum1    => this%list%Get_ptr_r3d('zonal wind')
      __acc_attach(this%pum1)
      this%pvm1    => this%list%Get_ptr_r3d('meridional wind')
      __acc_attach(this%pvm1)
      this%pwp1    => this%list%Get_ptr_r3d('vertical wind')
      __acc_attach(this%pwp1)
      this%pqm1    => this%list%Get_ptr_r3d('water vapor')
      __acc_attach(this%pqm1)
      this%pxlm1   => this%list%Get_ptr_r3d('cloud water')
      __acc_attach(this%pxlm1)
      this%pxim1   => this%list%Get_ptr_r3d('cloud ice')
      __acc_attach(this%pxim1)
      this%pxrm1   => this%list%Get_ptr_r3d('rain')
      __acc_attach(this%pxrm1)
      this%pxsm1   => this%list%Get_ptr_r3d('snow')
      __acc_attach(this%pxsm1)
      this%pxgm1   => this%list%Get_ptr_r3d('graupel')
      __acc_attach(this%pxgm1)
      this%papm1   => this%list%Get_ptr_r3d('full level pressure')
      __acc_attach(this%papm1)
      this%paphm1  => this%list%Get_ptr_r3d('half level pressure')
      __acc_attach(this%paphm1)
      this%dz      => this%list%Get_ptr_r3d('layer thickness')
      __acc_attach(this%dz)
#ifdef __MIXED_PRECISION
      this%inv_dzf => this%list%Get_ptr_s3d('inverse layer thickness full')
#else
      this%inv_dzf => this%list%Get_ptr_r3d('inverse layer thickness full')
#endif
      __acc_attach(this%inv_dzf)
      this%inv_dzh => this%list%Get_ptr_r3d('inverse layer thickness half')
      __acc_attach(this%inv_dzh)
  END SELECT

  END SUBROUTINE Set_pointers_inputs

  FUNCTION build_atmo_diagnostic_list(domain) RESULT(diaglist)

    TYPE(t_variable_list) :: diaglist
    TYPE(t_domain),        INTENT(in)    :: domain

    INTEGER :: nproma, nblks_c, nblks_e, nblks_v, nlev, nlevp1
    INTEGER :: shape_2d(2), shape_3d(3)

    diaglist = t_variable_list('diagnostics')

    nproma  = domain%nproma
    nblks_c = domain%nblks_c
    nblks_e = domain%nblks_e
    nblks_v = domain%nblks_v
    nlev    = domain%nlev
    nlevp1  = nlev + 1

    shape_3d = [nproma,nlev,nblks_v]
    CALL diaglist%append(t_variable('zonal wind vertice', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('meridional wind vertice', shape_3d, "m/s", type_id="real"))

    shape_3d = [nproma,nlevp1,nblks_v]
    CALL diaglist%append(t_variable('vertical wind vertice', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient momentum interface vertice', shape_3d, "m2/s", type_id="real"))

    shape_3d = [nproma,nlev,nblks_e]
    CALL diaglist%append(t_variable('normal wind', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('shear', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('stress div', shape_3d, "", type_id="real"))

    shape_3d = [nproma,nlevp1,nblks_e]
    CALL diaglist%append(t_variable('normal wind at edge', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('tangential wind at edge', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('vertical wind at edge', shape_3d, "m/s", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient momentum interface edge', shape_3d, "m2/s", type_id="real"))

    shape_3d = [nproma,nlev,nblks_c]
    CALL diaglist%append(t_variable('full level geopotential height above ground', shape_3d, "m", type_id="real"))
    CALL diaglist%append(t_variable('static energy', shape_3d, "m2 s-2", type_id="real"))
    CALL diaglist%append(t_variable('divergence cell', shape_3d, "", type_id="real"))
    ! CALL diaglist%append(t_variable('diff coefficient scalar', shape_3d, "m2/s", type_id="real"))
    ! CALL diaglist%append(t_variable('diff coefficient momentum', shape_3d, "m2/s", type_id="real"))
    CALL diaglist%append(t_variable('layer heating', shape_3d, "W/m2", type_id="real"))
    CALL diaglist%append(t_variable('dissipation of kinetic energy', shape_3d, "W/m2", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient momentum full', shape_3d, "m2/s", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient scalar full', shape_3d, "m2/s", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient momentum for hor. diff.', shape_3d, "m2/s", type_id="real"))
    CALL diaglist%append(t_variable('virtual potential temperature', shape_3d, "K", type_id="real"))
    CALL diaglist%append(t_variable('rho over dz', shape_3d, "", type_id="real"))

    shape_3d = [nproma,nlevp1,nblks_c]
    CALL diaglist%append(t_variable('air density interface', shape_3d, "kg/m3", type_id="real"))
    CALL diaglist%append(t_variable('brunt vaisal freq', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('stability function', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('mechanical production', shape_3d, "", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient momentum interface', shape_3d, "m2/s", type_id="real"))
    CALL diaglist%append(t_variable('exchange coefficient scalar interface', shape_3d, "m2/s", type_id="real"))

    shape_2d = [nproma,nblks_c]
    CALL diaglist%append(t_variable('latent heat flux surface', shape_2d, "W/m2", type_id="real"))
    CALL diaglist%append(t_variable('sensible heat flux surface', shape_2d, "W/m2", type_id="real"))
    CALL diaglist%append(t_variable('static energy, vert. int.', shape_2d, "m2 s-2", type_id="real"))
    CALL diaglist%append(t_variable('dissipation of kinetic energy, vert. int.', shape_2d, "W/m2", type_id="real"))
    CALL diaglist%append(t_variable('scaling factor for Louis constant b', shape_2d, "", type_id="real"))
    CALL diaglist%append(t_variable('moist internal energy after tmx, vert. int.', shape_2d, "J m-2", type_id="real"))
    CALL diaglist%append(t_variable('tendency of vert. int. moist internal energy', shape_2d, "J m-2 s-1", type_id="real"))

    shape_3d = [nproma,nblks_c,nsfc_type]
    ! CALL diaglist%append(t_variable('latent heat flux surface tile', shape_3d, "W/m2", type_id="real"))
    ! CALL diaglist%append(t_variable('sensible heat flux surface tile', shape_3d, "W/m2", type_id="real"))
    ! CALL diaglist%append(t_variable('density surface tile', shape_3d, "kg/m3", type_id="real"))
    ! CALL diaglist%append(t_variable('saturation specific humidity', shape_3d, "kg/kg", type_id="real"))
    ! CALL diaglist%append(t_variable('dry static energy', shape_3d, "m2/s2", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for momentum', shape_3d, "-", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for scalar', shape_3d, "-", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for diagnostics 1', shape_3d, "-", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for diagnostics 2', shape_3d, "-", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for diagnostics 3', shape_3d, "-", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for diagnostics 4', shape_3d, "-", type_id="real"))
    ! CALL diaglist%append(t_variable('drag coefficient for diagnostics 5', shape_3d, "-", type_id="real"))

  END FUNCTION build_atmo_diagnostic_list

  SUBROUTINE Set_pointers_diagnostics(this)

    CLASS(t_vdf_atmo_diagnostics), INTENT(inout) :: this

    SELECT TYPE (this)
    TYPE IS (t_vdf_atmo_diagnostics)
      this%ghf => this%list%Get_ptr_r3d('full level geopotential height above ground')
      __acc_attach(this%ghf)
      this%ctgz => this%list%Get_ptr_r3d('static energy')
      __acc_attach(this%ctgz)
      this%ctgzvi => this%list%get_ptr_r2d('static energy, vert. int.')
      __acc_attach(this%ctgzvi)
    ! this%kh => this%list%Get_ptr_r3d('diff coefficient scalar')
      ! this%km => this%list%Get_ptr_r3d('diff coefficient momentum')
      this%km   => this%list%Get_ptr_r3d('exchange coefficient momentum full')
      __acc_attach(this%km)
      this%kh   => this%list%Get_ptr_r3d('exchange coefficient scalar full')
      __acc_attach(this%kh)
      this%km_c => this%list%Get_ptr_r3d('exchange coefficient momentum for hor. diff.')
      __acc_attach(this%km_c)
      this%div_c => this%list%Get_ptr_r3d('divergence cell')
      __acc_attach(this%div_c)
      this%theta_v => this%list%Get_ptr_r3d('virtual potential temperature')
      __acc_attach(this%theta_v)
      this%pprfac => this%list%Get_ptr_r3d('rho over dz')
      __acc_attach(this%pprfac)
      this%rho_ic => this%list%Get_ptr_r3d('air density interface')
      __acc_attach(this%rho_ic)
      this%km_ic => this%list%Get_ptr_r3d('exchange coefficient momentum interface')
      __acc_attach(this%km_ic)
      this%kh_ic => this%list%Get_ptr_r3d('exchange coefficient scalar interface')
      __acc_attach(this%kh_ic)
      this%bruvais => this%list%Get_ptr_r3d('brunt vaisal freq')
      __acc_attach(this%bruvais)
      this%stability_function => this%list%Get_ptr_r3d('stability function')
      __acc_attach(this%stability_function)
      this%vn_ie => this%list%Get_ptr_r3d('normal wind at edge')
      __acc_attach(this%vn_ie)
      this%vt_ie => this%list%Get_ptr_r3d('tangential wind at edge')
      __acc_attach(this%vt_ie)
      this%w_ie  => this%list%Get_ptr_r3d('vertical wind at edge')
      __acc_attach(this%w_ie)
      this%km_ie => this%list%Get_ptr_r3d('exchange coefficient momentum interface edge')
      __acc_attach(this%km_ie)
      this%vn    => this%list%Get_ptr_r3d('normal wind')
      __acc_attach(this%vn)
      this%shear => this%list%Get_ptr_r3d('shear')
      __acc_attach(this%shear)
      this%div_of_stress => this%list%Get_ptr_r3d('stress div')
      __acc_attach(this%div_of_stress)
      this%mech_prod => this%list%Get_ptr_r3d('mechanical production')
      __acc_attach(this%mech_prod)
      this%u_vert => this%list%Get_ptr_r3d('zonal wind vertice')
      __acc_attach(this%u_vert)
      this%v_vert => this%list%Get_ptr_r3d('meridional wind vertice')
      __acc_attach(this%v_vert)
      this%w_vert => this%list%Get_ptr_r3d('vertical wind vertice')
      __acc_attach(this%w_vert)
      this%km_iv => this%list%Get_ptr_r3d('exchange coefficient momentum interface vertice')
      __acc_attach(this%km_iv)
      this%heating => this%list%Get_ptr_r3d('layer heating')
      __acc_attach(this%heating)
      this%dissip_kin_energy => this%list%Get_ptr_r3d('dissipation of kinetic energy')
      __acc_attach(this%dissip_kin_energy)
      this%dissip_kin_energy_vi => this%list%Get_ptr_r2d('dissipation of kinetic energy, vert. int.')
      __acc_attach(this%dissip_kin_energy_vi)
      this%scaling_factor_louis => this%list%Get_ptr_r2d('scaling factor for Louis constant b')
      __acc_attach(this%scaling_factor_louis)
      this%internal_energy_vi => this%list%Get_ptr_r2d('moist internal energy after tmx, vert. int.')
      __acc_attach(this%internal_energy_vi)
      this%internal_energy_vi_tend => this%list%Get_ptr_r2d('tendency of vert. int. moist internal energy')
      __acc_attach(this%internal_energy_vi_tend)
      ! boundary condition
      ! this%lhfl  => this%list%Get_ptr_r2d('latent heat flux surface')
      ! this%shfl  => this%list%Get_ptr_r2d('sensible heat flux surface')
      ! this%lhfl_tile => this%list%Get_ptr_r3d('latent heat flux surface tile')
      ! this%shfl_tile => this%list%Get_ptr_r3d('sensible heat flux surface tile')
      ! this%rho_tile => this%list%Get_ptr_r3d('density surface tile')
      ! this%pqsat_tile => this%list%Get_ptr_r3d('saturation specific humidity')
      ! this%pcpt_tile => this%list%Get_ptr_r3d('dry static energy')
      ! this%pcfm_tile => this%list%Get_ptr_r3d('drag coefficient for momentum')
      ! this%pcfh_tile => this%list%Get_ptr_r3d('drag coefficient for scalar')
      ! this%pch_tile => this%list%Get_ptr_r3d('drag coefficient for diagnostics 1')
      ! this%pbn_tile => this%list%Get_ptr_r3d('drag coefficient for diagnostics 2')
      ! this%pbhn_tile => this%list%Get_ptr_r3d('drag coefficient for diagnostics 3')
      ! this%pbm_tile => this%list%Get_ptr_r3d('drag coefficient for diagnostics 4')
      ! this%pbh_tile => this%list%Get_ptr_r3d('drag coefficient for diagnostics 5')
    END SELECT

  END SUBROUTINE Set_pointers_diagnostics

  SUBROUTINE temp_to_energy(this, temperature, energy, use_new_moisture_state)

    CLASS(t_vdf_atmo), INTENT(in) :: this
    LOGICAL,  INTENT(in), OPTIONAL :: use_new_moisture_state
    REAL(wp), INTENT(in) :: temperature(:,:,:)
    REAL(wp), INTENT(out) :: energy(:,:,:)

    INTEGER :: jb, jk, jc
    INTEGER,  POINTER :: energy_type
    LOGICAL :: use_updated_moisture
    REAL(wp), POINTER :: cpd
    REAL(wp), POINTER :: geo_height(:,:,:)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & qr, qs, qg, qv, qc, qi

    CHARACTER(len=*), PARAMETER :: routine = modname//':temp_to_energy'

    energy_type      => this%config%list%Get_ptr_i0d('energy type')

    ! No effect for this%energy_type=1
    use_updated_moisture = .TRUE.
    IF (PRESENT(use_new_moisture_state)) use_updated_moisture = use_new_moisture_state

    geo_height => this%diagnostics%list%get_ptr_r3d('full level geopotential height above ground')

    SELECT CASE(energy_type)
    CASE (1)
      cpd => this%config%list%Get_ptr_r0d('cpd')
      CALL compute_static_energy(this%domain, cpd, temperature, geo_height, energy)
    CASE (2)
      qr  => this%inputs%list%get_ptr_r3d('rain')
      qs  => this%inputs%list%get_ptr_r3d('snow')
      qg  => this%inputs%list%get_ptr_r3d('graupel')
      IF (use_updated_moisture) THEN
        qv  => this%new_states%get_ptr_r3d('water vapor')
        qc  => this%new_states%get_ptr_r3d('cloud water')
        qi  => this%new_states%get_ptr_r3d('cloud ice')
      ELSE
        qv  => this%states%get_ptr_r3d('water vapor')
        qc  => this%states%get_ptr_r3d('cloud water')
        qi  => this%states%get_ptr_r3d('cloud ice')
      END IF
      CALL compute_internal_energy( &
        & this%domain, &
        & geo_height,  &
        & qr, qs, qg,  &
        & temperature, &
        & qv, qc, qi,  &
        & energy       &
        & )
    END SELECT

  END SUBROUTINE temp_to_energy

  SUBROUTINE energy_to_temp(this, energy, temperature, use_new_moisture_state)

    CLASS(t_vdf_atmo), INTENT(in) :: this
    LOGICAL,  INTENT(in), OPTIONAL :: use_new_moisture_state
    REAL(wp), INTENT(in)  :: energy(:,:,:)
    REAL(wp), INTENT(out) :: temperature(:,:,:)

    INTEGER :: jb, jk, jc
    INTEGER,  POINTER :: energy_type
    LOGICAL :: use_updated_moisture
    REAL(wp), POINTER :: cpd
    REAL(wp), POINTER :: geo_height(:,:,:)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & qr, qs, qg, qv, qc, qi

    CHARACTER(len=*), PARAMETER :: routine = modname//':energy_to_temp'

    energy_type => this%config%list%Get_ptr_i0d('energy type')

    ! No effect for this%energy_type=1
    use_updated_moisture = .TRUE.
    IF (PRESENT(use_new_moisture_state)) use_updated_moisture = use_new_moisture_state

    geo_height => this%diagnostics%list%get_ptr_r3d('full level geopotential height above ground')

    SELECT CASE(energy_type)
    CASE (1)
      cpd => this%config%list%Get_ptr_r0d('cpd')
      CALL compute_temp_from_static_energy(this%domain, cpd, energy, geo_height, temperature)
    CASE (2)
      qr  => this%inputs%list%get_ptr_r3d('rain')
      qs  => this%inputs%list%get_ptr_r3d('snow')
      qg  => this%inputs%list%get_ptr_r3d('graupel')
      IF (use_updated_moisture) THEN
        qv  => this%new_states%get_ptr_r3d('water vapor')
        qc  => this%new_states%get_ptr_r3d('cloud water')
        qi  => this%new_states%get_ptr_r3d('cloud ice')
      ELSE
        qv  => this%states%get_ptr_r3d('water vapor')
        qc  => this%states%get_ptr_r3d('cloud water')
        qi  => this%states%get_ptr_r3d('cloud ice')
      END IF
      CALL compute_temperature_from_internal_energy( &
        & this%domain, &
        & geo_height,  &
        & qr, qs, qg,  &
        & energy,      &
        & qv, qc, qi,  &
        & temperature  &
        & )
      END SELECT

  END SUBROUTINE energy_to_temp

  SUBROUTINE compute_flux_x(this, shflx, ufts, ufvs, flux_x)

    CLASS(t_vdf_atmo), INTENT(in) :: this
    REAL(wp), INTENT(in) :: &
      & shflx(:,:), & !< sensible heat flux
      & ufts(:,:),  & !< energy flux at surface from thermal exchange
      & ufvs(:,:)     !< energy flux at surface from vapor exchange
    REAL(wp), INTENT(out) :: flux_x(:,:)

    INTEGER :: jb, jc
    INTEGER,  POINTER :: energy_type
    REAL(wp), POINTER :: cpd, cvd

    CHARACTER(len=*), PARAMETER :: routine = modname//':energy_flux_to_flux_x'

    energy_type      => this%config%list%Get_ptr_i0d('energy type')

    ASSOCIATE(domain => this%domain)

!$OMP PARALLEL
    CALL init(flux_x, lacc=.TRUE.)
!$OMP END PARALLEL

    SELECT CASE(energy_type)
    CASE (1)
      cpd => this%config%list%Get_ptr_r0d('cpd')
      cvd => this%config%list%Get_ptr_r0d('cvd')

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c, domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          flux_x(jc,jb) = shflx(jc,jb) * cpd / cvd
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    CASE (2)

!$OMP PARALLEL DO PRIVATE(jb, jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c, domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          flux_x(jc,jb) = ufts(jc,jb) + ufvs(jc,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    END SELECT

    END ASSOCIATE

    !$ACC WAIT(1)

  END SUBROUTINE compute_flux_x
  !
  !=================================================================
  !
  ! Subroutines for diagnostics
  !
  SUBROUTINE compute_static_energy(domain, spec_heat, temperature, geo_height, static_energy)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: spec_heat
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & temperature, geo_height
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & static_energy

    INTEGER :: jb, jk, jc

!$OMP PARALLEL
    CALL init(static_energy, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          static_energy(jc,jk,jb) = spec_heat * temperature(jc,jk,jb) + grav * geo_height(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_static_energy

  SUBROUTINE compute_temp_from_static_energy(domain, spec_heat, static_energy, geo_height, temperature)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), INTENT(in) :: spec_heat
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & static_energy, geo_height
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & temperature

    INTEGER :: jb, jk, jc

!$OMP PARALLEL
    CALL init(temperature, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          temperature(jc,jk,jb) = (static_energy(jc,jk,jb) - grav * geo_height(jc,jk,jb)) / spec_heat                                    
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_temp_from_static_energy

  SUBROUTINE compute_geopotential_height_above_ground(domain, zf, zh, ghf)

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & zf, zh
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & ghf

    INTEGER :: jb, jk, jc, nlevp1

    nlevp1 = domain%nlev + 1

!$OMP PARALLEL
    CALL init(ghf, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          ghf(jc,jk,jb) = zf(jc,jk,jb) - zh(jc,nlevp1,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_geopotential_height_above_ground

  ! Compute mass specific internal energy + geopotential from temperature and moisture state
  SUBROUTINE compute_internal_energy( &
    & domain,      &
    & geo_height,  &
    & qr, qs, qg,  &
    & temperature, &
    & qv, qc, qi,  &
    & energy       &
    & )

    USE mo_aes_thermo, ONLY: internal_energy

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & geo_height,  &
      & qr,          &
      & qs,          &
      & qg,          &
      & temperature, &
      & qv,          &
      & qc,          &
      & qi
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & energy

    INTEGER :: jb, jk, jc
    REAL(wp) :: q_liquid, q_solid

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_internal_energy'

!$OMP PARALLEL
    CALL init(energy, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, q_liquid, q_solid) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q_liquid, q_solid)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          q_liquid = qc(jc,jk,jb) + qr(jc,jk,jb)
          q_solid  = qi(jc,jk,jb) + qs(jc,jk,jb) + qg(jc,jk,jb)
          energy(jc,jk,jb) = &
            & internal_energy(           &
            &   temperature(jc,jk,jb), & ! temperature
            &   qv         (jc,jk,jb), & ! qv
            &   q_liquid,              & ! liquid
            &   q_solid,               & ! solid
            &   1._wp,                 & ! density
            &   1._wp                  & ! delta z
            &) + grav * geo_height(jc,jk,jb) * cvd/cpd
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_internal_energy

  ! Compute temperature from mass specific internal energy + geopotential and moisture state
  SUBROUTINE compute_temperature_from_internal_energy( &
    & domain,      &
    & geo_height,  &
    & qr, qs, qg,  &
    & energy,      &
    & qv, qc, qi,  &
    & temperature  &
    & )

    USE mo_aes_thermo, ONLY: T_from_internal_energy

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & geo_height, &
      & qr,         &
      & qs,         &
      & qg,         &
      & energy,     &
      & qv,         &
      & qc,         &
      & qi
    REAL(wp), DIMENSION(:,:,:), INTENT(out) :: &
      & temperature

    INTEGER :: jb, jk, jc
    REAL(wp) :: q_liquid, q_solid, u

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_temperature_from_internal_energy'

!$OMP PARALLEL
    CALL init(temperature, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, q_liquid, q_solid, u) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q_liquid, q_solid, u)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          q_liquid = qc(jc,jk,jb) + qr(jc,jk,jb)
          q_solid  = qi(jc,jk,jb) + qs(jc,jk,jb) + qg(jc,jk,jb)
          u        = energy(jc,jk,jb) - grav * geo_height(jc,jk,jb) * cvd/cpd
          temperature(jc,jk,jb) = &
            & T_from_internal_energy(  &
            &   u,                     & ! internal energy
            &   qv         (jc,jk,jb), & ! qv
            &   q_liquid,              & ! liquid
            &   q_solid,               & ! solid
            &   1._wp,                 & ! density
            &   1._wp                  & ! delta z
            &)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_temperature_from_internal_energy

  SUBROUTINE compute_internal_energy_vi( &
    & domain,      &
    & rho,         &
    & dz,          &
    & qr, qs, qg,  &
    & temperature, &
    & qv, qc, qi,  &
    & uvi          &
    & )

    USE mo_aes_thermo, ONLY: internal_energy

    TYPE(t_domain), INTENT(in), POINTER :: domain
    REAL(wp), DIMENSION(:,:,:), INTENT(in) :: &
      & rho, &
      & dz, &
      & qr, &
      & qs, &
      & qg, &
      & temperature, &
      & qv, &
      & qc, &
      & qi
    REAL(wp), DIMENSION(:,:), INTENT(out) :: &
      & uvi

    INTEGER :: jb, jk, jc
    REAL(wp) :: q_liquid, q_solid

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_internal_energy_vi'

!$OMP PARALLEL
    CALL init(uvi, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, q_liquid, q_solid) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c,domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1,domain%nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(q_liquid, q_solid)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          q_liquid = qc(jc,jk,jb) + qr(jc,jk,jb)
          q_solid  = qi(jc,jk,jb) + qs(jc,jk,jb) + qg(jc,jk,jb)
          uvi(jc,jb) = uvi(jc,jb) + &
            & internal_energy(           &
            &   temperature(jc,jk,jb), & ! temperature
            &   qv         (jc,jk,jb), & ! qv
            &   q_liquid,              & ! liquid
            &   q_solid,               & ! solid
            &   rho(jc,jk,jb),         & ! density
            &   dz (jc,jk,jb)          & ! delta z
            &)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)

  END SUBROUTINE compute_internal_energy_vi
  !
  !#########################################################################
  !## Smagorinsky_model
  !!------------------------------------------------------------------------
  !! Computes the sgs viscosity and diffusivity using Smagorinsky model
  !! \tau_ij = KD_ij where D_ij = du_i/dx_j + du_j/dx_i
  !!
  !! and  K = cs * \Delta**2 * D / sqrt(2), where D = sqrt(D_ijD_ij)
  !!
  !! and D**2 = D_11**2 + D_22**2 + D_33**2 + 2D_12**2 + 2D_13**2 + 2D_23**2
  !!
  !! where, D_11 = 2 * du_1/dx_1
  !!        D_22 = 2 * d_u2/dx_2
  !!        D_33 = 2 * d_u3/dx_3
  !!        D_12 = du_1/dx_2 + du_2/dx_1
  !!        D_13 = du_1/dx_3 + du_3/dx_1
  !!        D_23 = du_2/dx_3 + du_3/dx_2
  !! For triangles: 1=normal, 2=tangential, and 3 = z directions
  !!------------------------------------------------------------------------
  SUBROUTINE compute_exchange_coefficient(domain, config, inputs, diagnostics)

    TYPE(t_domain),        INTENT(in),   POINTER :: domain
    CLASS(t_variable_set), INTENT(in),    TARGET :: config
    CLASS(t_variable_set), INTENT(in),    TARGET :: inputs
    CLASS(t_variable_set), INTENT(inout), TARGET :: diagnostics

    TYPE(t_vdf_atmo_config),      POINTER :: conf
    TYPE(t_vdf_atmo_inputs),      POINTER :: ins
    TYPE(t_vdf_atmo_diagnostics), POINTER :: diags
    TYPE(t_int_state),            POINTER :: p_int         !< interpolation state
    TYPE(t_nh_metrics),           POINTER :: p_nh_metrics
    TYPE(t_patch),                POINTER :: patch

    INTEGER, DIMENSION(:,:,:), POINTER :: ividx, ivblk, iecidx, iecblk, ieidx, ieblk

    INTEGER :: jg, jl, jk
    INTEGER :: jb,jc,je
    INTEGER :: jb_max, jk_max, jc_max
    INTEGER :: jcn, jbn
    INTEGER :: nlev, nlevm1, nlevp1
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end

    REAL(wp), DIMENSION(:,:), POINTER :: area

    REAL(wp) :: zvn1, zvn2
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4
    REAL(wp) :: vt_vert1, vt_vert2, vt_vert3, vt_vert4
    REAL(wp) :: w_full_c1, w_full_c2, w_full_v1, w_full_v2
    REAL(wp) :: D_11, D_12, D_13, D_22, D_23, D_33
    REAL(wp) :: Ri

    REAL(wp), POINTER :: rturb_prandtl

    ! Global mean of cell area for R2B8 [m]
    REAL(wp), PARAMETER :: mean_area_R2B8 = 97294071.23714285_wp

    CHARACTER(len=*), PARAMETER :: routine = modname//':compute_exchange_coefficients'

    patch => domain%patch

    SELECT TYPE (config)
    TYPE IS (t_vdf_atmo_config)
      conf => config
    END SELECT
    __acc_attach(conf)
    SELECT TYPE (inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins => inputs
    END SELECT
    __acc_attach(ins)
    SELECT TYPE (diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags => diagnostics
    END SELECT
    __acc_attach(diags)

    ! Use ASSOCIATE for readability
    ! (not sure if this is ok with OpenACC!)
    ASSOCIATE ( &
      nproma  => domain%nproma,  &
      nblks_c => domain%nblks_c, &
      nlev    => domain%nlev,    &
      ! i_startblk_c => domain%i_startblk_c,    & ! Start block on cells
      ! i_endblk_c   => domain%i_endblk_c,      & ! End block on cells
      ! i_startidx_c => domain%i_startidx_c(:), & ! Start indices on cells (for each block)
      ! i_endidx_c   => domain%i_endidx_c(:),   & ! End indices on cells (for each block)
      ! i_startblk_e => domain%i_startblk_e,    & ! Start block on edges
      ! i_endblk_e   => domain%i_endblk_e,      & ! End block on edges
      ! i_startidx_e => domain%i_startidx_e(:), & ! Start indices on edges (for each block)
      ! i_endidx_e   => domain%i_endidx_e(:),   & ! End indices on edges (for each block)
      !
      dz       => ins%dz,        & !< layer thickness on full levels
      inv_dzf  => ins%inv_dzf,   & !< inverse of layer thickness on full levels
      inv_dzh  => ins%inv_dzh,   & !< inverse of layer thickness on half levels
      zf       => ins%zf,        & 
      zh       => ins%zh,        & 
      ptm1     => ins%ptm1,      &
      ptvm1    => ins%ptvm1,     &
      rho      => ins%rho,       &
      pum1     => ins%pum1,      &
      pvm1     => ins%pvm1,      &
      pwp1     => ins%pwp1,      &
      papm1    => ins%papm1,     &
      paphm1   => ins%paphm1,    &
      theta_v  => diags%theta_v, &
      scaling_factor_louis => diags%scaling_factor_louis, &
      pprfac   => diags%pprfac,  &
      bruvais  => diags%bruvais, &
      stability_function => diags%stability_function, &
      rho_ic   => diags%rho_ic,  &
      km_min   => conf%km_min,   &
      ! kh       => diags%kh,      &
      ! km       => diags%km,      &
      km_c     => diags%km_c,    &
      km_ic    => diags%km_ic,   &
      kh_ic    => diags%kh_ic,   &
      km_iv    => diags%km_iv,   &
      vn_ie    => diags%vn_ie,   &
      vt_ie    => diags%vt_ie,   &
      w_ie     => diags%w_ie,    &
      km_ie    => diags%km_ie,   &
      vn       => diags%vn,      &
      shear    => diags%shear,   &
      div_of_stress   => diags%div_of_stress,&
      mech_prod=> diags%mech_prod,&
      u_vert   => diags%u_vert,  &
      v_vert   => diags%v_vert,  &
      w_vert   => diags%w_vert,  &
      div_c    => diags%div_c,   &
      ! rturb_prandtl=> conf%rturb_prandtl,&
      use_louis => conf%use_louis, &
      louis_constant_b => conf%louis_constant_b, &
      turb_prandtl=> conf%turb_prandtl &
      )

    jg = patch%id             ! Only jg=1 supported currently!
    p_int  => p_int_state(jg)
    p_nh_metrics => p_nh_state(jg)%metrics

    nlevm1 = nlev-1
    nlevp1 = nlev+1

    ! Scaling factor for Louis constant b, designed to be 1 with R2B8
    area => patch%cells%area
!$OMP PARALLEL DO PRIVATE(jb,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        scaling_factor_louis(jc,jb) = mean_area_R2B8 / area(jc,jb)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    CALL sync_patch_array(SYNC_C, patch, pum1)
    CALL sync_patch_array(SYNC_C, patch, pvm1)

    rturb_prandtl=> conf%rturb_prandtl

    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, jcn, jbn, zvn1, zvn2) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk ! domain%i_startblk_e, domain%i_endblk_e
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2) PRIVATE(jcn, jbn, zvn1, zvn2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx ! domain%i_startidx_e(jb), domain%i_endidx_e(jb)
          jcn  = patch%edges%cell_idx(je,jb,1)
          jbn  = patch%edges%cell_blk(je,jb,1)

          zvn1 =   pum1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v1 &
            &    + pvm1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,1)%v2
          !
          jcn  =   patch%edges%cell_idx(je,jb,2)
          jbn  =   patch%edges%cell_blk(je,jb,2)
          zvn2 =   pum1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v1 &
            &    + pvm1(jcn,jk,jbn)*patch%edges%primal_normal_cell(je,jb,2)%v2
          !
          vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)*zvn1 &
            &          + p_int%c_lin_e(je,2,jb)*zvn2
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

   CALL sync_patch_array(SYNC_E, patch, vn)


!#########################################################################
!## initialize
!#########################################################################

!$OMP PARALLEL
    CALL init(km_iv, lacc=.TRUE.)
    CALL init(km_c, lacc=.TRUE.)
    CALL init(km_ie, lacc=.TRUE.)
    CALL init(kh_ic, lacc=.TRUE.)
    CALL init(km_ic, lacc=.TRUE.)
    CALL init(u_vert, lacc=.TRUE.)
    CALL init(v_vert, lacc=.TRUE.)
    CALL init(w_vert, lacc=.TRUE.)
!$OMP END PARALLEL

    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          theta_v(jc,jk,jb) = ptvm1(jc,jk,jb)*(p0ref/papm1(jc,jk,jb))**rd_o_cpd
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !Get rho at interfaces to be used later
    CALL vert_intp_full2half_cell_3d(patch, p_nh_metrics, rho, rho_ic, &
                                     2, min_rlcell_int-2, lacc=.TRUE.)

    CALL brunt_vaisala_freq(patch, p_nh_metrics, nproma, theta_v, bruvais, &
                            opt_rlstart=3, lacc=.TRUE.)

    !--------------------------------------------------------------------------
    !1) Interpolate velocities at desired locations- mostly around the quadrilateral
    !
    !<It assumes that prog values are all synced at this stage while diag values might not>
    !--------------------------------------------------------------------------


    CALL cells2verts_scalar(pwp1, patch, p_int%cells_aw_verts, w_vert,                   &
                            opt_rlend=min_rlvert_int, opt_acc_async=.TRUE.)
    CALL cells2edges_scalar(pwp1, patch, p_int%c_lin_e, w_ie, opt_rlend=min_rledge_int-2,&
                            lacc=.TRUE.)

    ! RBF reconstruction of velocity at vertices: include halos
    CALL rbf_vec_interpol_vertex(vn, patch, p_int, u_vert, v_vert,                       &
                                 opt_rlend=min_rlvert_int)

    !sync them
    CALL sync_patch_array_mult(SYNC_V, patch, 3, w_vert, u_vert, v_vert)

    !Get vn at interfaces and then get vt at interfaces
    !Boundary values are extrapolated like dynamics although
    !they are not required in current implementation


    rl_start   = 2
    rl_end     = min_rledge_int-3
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx
          vn_ie(je,jk,jb) = p_nh_metrics%wgtfac_e(je,jk,jb) * vn(je,jk,jb) +            &
                            ( 1._wp - p_nh_metrics%wgtfac_e(je,jk,jb) ) * vn(je,jk-1,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO je = i_startidx, i_endidx
        vn_ie(je,1,jb)      = p_nh_metrics%wgtfacq1_e(je,1,jb) * vn(je,1,jb) +          &
                              p_nh_metrics%wgtfacq1_e(je,2,jb) * vn(je,2,jb) +          &
                              p_nh_metrics%wgtfacq1_e(je,3,jb) * vn(je,3,jb)

        vn_ie(je,nlevp1,jb) = p_nh_metrics%wgtfacq_e(je,1,jb) * vn(je,nlev,jb)   +      &
                              p_nh_metrics%wgtfacq_e(je,2,jb) * vn(je,nlev-1,jb) +      &
                              p_nh_metrics%wgtfacq_e(je,3,jb) * vn(je,nlev-2,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO


    CALL rbf_vec_interpol_edge(vn_ie, patch, p_int, vt_ie, opt_rlstart=3, &
                               opt_rlend=min_rledge_int-2)

    !--------------------------------------------------------------------------
    !2) Compute horizontal strain rate tensor at full levels
    !--------------------------------------------------------------------------
!    ividx => patch%edges%vertex_idx
!    ivblk => patch%edges%vertex_blk

!    iecidx => patch%edges%cell_idx
!    iecblk => patch%edges%cell_blk

!    ieidx => patch%cells%edge_idx
!    ieblk => patch%cells%edge_blk

    rl_start   = 4
    rl_end     = min_rledge_int-2
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, vn_vert1, vn_vert2, vn_vert3, vn_vert4,  &
!$OMP                     vt_vert1, vt_vert2, vt_vert3, vt_vert4, w_full_c1, w_full_c2, w_full_v1, &
!$OMP                     w_full_v2, D_11, D_12, D_13, D_22, D_23, D_33) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk

      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(u_vert, v_vert)
      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(vn_vert1, vn_vert2, vn_vert3, vn_vert4, vt_vert1, vt_vert2, vt_vert3, vt_vert4) &
      !$ACC   PRIVATE(w_full_c1, w_full_c2, w_full_v1, w_full_v2, D_11, D_12, D_13, D_22, D_23, D_33)
#ifdef __LOOP_EXCHANGE
      DO je = i_startidx, i_endidx
        DO jk = 1, nlev
#else
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
#endif

          vn_vert1 = u_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1))     *   &
                     patch%edges%primal_normal_vert(je,jb,1)%v1 +   &
                     v_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1))     *   &
                     patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 = u_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2))     *   &
                     patch%edges%primal_normal_vert(je,jb,2)%v1 +   &
                     v_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2))     *   &
                     patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 = u_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3))     *   &
                     patch%edges%primal_normal_vert(je,jb,3)%v1 +   &
                     v_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3))     *   &
                     patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 = u_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4))     *   &
                     patch%edges%primal_normal_vert(je,jb,4)%v1 +   &
                     v_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4))     *   &
                     patch%edges%primal_normal_vert(je,jb,4)%v2

          vt_vert1 = u_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1))     *   &
                     patch%edges%dual_normal_vert(je,jb,1)%v1   +   &
                     v_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1))     *   &
                     patch%edges%dual_normal_vert(je,jb,1)%v2

          vt_vert2 = u_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2))     *   &
                     patch%edges%dual_normal_vert(je,jb,2)%v1   +   &
                     v_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2))     *   &
                     patch%edges%dual_normal_vert(je,jb,2)%v2

          vt_vert3 = u_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3))     *   &
                     patch%edges%dual_normal_vert(je,jb,3)%v1   +   &
                     v_vert(patch%edges%vertex_idx(je,jb,3),jk,patch%edges%vertex_blk(je,jb,3))     *   &
                     patch%edges%dual_normal_vert(je,jb,3)%v2

          vt_vert4 = u_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4))     *   &
                     patch%edges%dual_normal_vert(je,jb,4)%v1   +   &
                     v_vert(patch%edges%vertex_idx(je,jb,4),jk,patch%edges%vertex_blk(je,jb,4))     *   &
                     patch%edges%dual_normal_vert(je,jb,4)%v2

          ! W at full levels
          w_full_c1  = 0.5_wp *                                              &
                       ( pwp1(patch%edges%cell_idx(je,jb,1),jk,patch%edges%cell_blk(je,jb,1)) +   &
                         pwp1(patch%edges%cell_idx(je,jb,1),jk+1,patch%edges%cell_blk(je,jb,1)) )

          w_full_c2  = 0.5_wp *                                              &
                       ( pwp1(patch%edges%cell_idx(je,jb,2),jk,patch%edges%cell_blk(je,jb,2)) +   &
                         pwp1(patch%edges%cell_idx(je,jb,2),jk+1,patch%edges%cell_blk(je,jb,2)) )

          ! W at full levels vertices from w at vertices at interface levels
          w_full_v1  = 0.5_wp *                                              &
                       ( w_vert(patch%edges%vertex_idx(je,jb,1),jk,patch%edges%vertex_blk(je,jb,1)) +          &
                         w_vert(patch%edges%vertex_idx(je,jb,1),jk+1,patch%edges%vertex_blk(je,jb,1)) )

          w_full_v2  = 0.5_wp *                                              &
                       ( w_vert(patch%edges%vertex_idx(je,jb,2),jk,patch%edges%vertex_blk(je,jb,2)) +          &
                         w_vert(patch%edges%vertex_idx(je,jb,2),jk+1,patch%edges%vertex_blk(je,jb,2)) )


          ! Strain rates at edge center
          D_11       = 2._wp * ( vn_vert4 - vn_vert3 ) *                     &
                       patch%edges%inv_vert_vert_length(je,jb)

          D_12       = patch%edges%tangent_orientation(je,jb) *            &
                       ( vn_vert2 - vn_vert1 ) *                             &
                       patch%edges%inv_primal_edge_length(je,jb) +         &
                       ( vt_vert4-vt_vert3 ) *                               &
                       patch%edges%inv_vert_vert_length(je,jb)

          D_13       = ( vn_ie(je,jk,jb) - vn_ie(je,jk+1,jb) ) *             &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +           &
                       ( w_full_c2 - w_full_c1 ) *                           &
                       patch%edges%inv_dual_edge_length(je,jb)

          D_22       = 2._wp * ( vt_vert2-vt_vert1 ) *                       &
                       patch%edges%tangent_orientation(je,jb) *            &
                       patch%edges%inv_primal_edge_length(je,jb)

          D_23       = ( vt_ie(je,jk,jb) - vt_ie(je,jk+1,jb) ) *             &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)  +           &
                       patch%edges%tangent_orientation(je,jb) *            &
                       ( w_full_v2 - w_full_v1 ) *                           &
                       patch%edges%inv_primal_edge_length(je,jb)

          D_33       = 2._wp * ( w_ie(je,jk,jb) - w_ie(je,jk+1,jb) ) *       &
                       p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb)

          ! Mechanical prod is half of this value divided by km
          shear(je,jk,jb) = D_11**2 + D_22**2 + D_33**2 +                    &
                            2._wp * ( D_12**2 + D_13**2 + D_23**2 )

          ! calculate divergence to get the deviatoric part of stress tensor in
          ! diffusion: D_11 - 1/3 * (D_11 + D_22 + D_33)
          div_of_stress(je,jk,jb) = 0.5_wp * ( D_11 + D_22 + D_33 )
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END PARALLEL DO

    !Interpolate mech production term from mid level edge to interface level cell
    !except top and bottom boundaries
    rl_start   = grf_bdywidth_c+1
    rl_end     = min_rlcell_int-1
    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx, i_endidx
          div_c(jc,jk,jb) =                                                                       &
                  ( div_of_stress(patch%cells%edge_idx(jc,jb,1),jk,patch%cells%edge_blk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +  &
                    div_of_stress(patch%cells%edge_idx(jc,jb,2),jk,patch%cells%edge_blk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +  &
                    div_of_stress(patch%cells%edge_idx(jc,jb,3),jk,patch%cells%edge_blk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    ! Interpolate mech. production term from mid level edge to interface level cell
    ! except top and bottom boundaries
    rl_start   = 3
    rl_end     = min_rlcell_int-1
    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk,      &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 2, nlev
        DO jc = i_startidx, i_endidx
          mech_prod(jc,jk,jb) = p_nh_metrics%wgtfac_c(jc,jk,jb) * (                      &
                      shear(patch%cells%edge_idx(jc,jb,1),jk,patch%cells%edge_blk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb)   +      &
                      shear(patch%cells%edge_idx(jc,jb,2),jk,patch%cells%edge_blk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb)   +      &
                      shear(patch%cells%edge_idx(jc,jb,3),jk,patch%cells%edge_blk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) ) +      &
                      ( 1._wp - p_nh_metrics%wgtfac_c(jc,jk,jb) ) * (                             &
                      shear(patch%cells%edge_idx(jc,jb,1),jk-1,patch%cells%edge_blk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) +      &
                      shear(patch%cells%edge_idx(jc,jb,2),jk-1,patch%cells%edge_blk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) +      &
                      shear(patch%cells%edge_idx(jc,jb,3),jk-1,patch%cells%edge_blk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO


    !--------------------------------------------------------------------------
    ! 3) Classical Smagorinsky model with stability correction due to Lilly 1962
    !    at interface cell centers. At this point mech_prod is twice the actual
    !    mechanical production term.
    !--------------------------------------------------------------------------
    ! MP = Mechanical prod term calculated above
    ! visc = mixing_length_sq * SQRT(MP/2) * SQRT(1-Ri/Pr) where
    ! Ri = (g / theta) * d_theta_dz / (MP/2), where Brunt_vaisala_freq (byncy prod term/kh)
    !    = (g / theta) * d_theta_dz.
    ! After simplification: visc = mixing_length_sq/SQRT(2) * SQRT[MP/2 - (Brunt_vaisala_frq/Pr)]
    ! Note that the factor SQRT(2) with mixing_length_sq is considered into the Smag constant
    rl_start   = 3
    rl_end     = min_rlcell_int
    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

    IF (use_louis) THEN

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
          DO jc = i_startidx, i_endidx
#endif

            Ri  = 2._wp * bruvais(jc,jk,jb) / mech_prod(jc,jk,jb) 

            Stability_function(jc,jk,jb) =                                &
            &  MAX(1.0_wp - Ri*rturb_prandtl,                             &
            &      MIN(1._wp, 1._wp/(1._wp+louis_constant_b*scaling_factor_louis(jc,jb)*ABS(Ri))**4))

            kh_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * rturb_prandtl *          &
                              p_nh_metrics%mixing_length_sq(jc,jk,jb) *   &
                              SQRT( mech_prod(jc,jk,jb) * 0.5_wp ) *      &
                              SQRT( Stability_function(jc,jk,jb) )

            km_ic(jc,jk,jb) = kh_ic(jc,jk,jb) * turb_prandtl
          END DO
        END DO
        !$ACC END PARALLEL
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
          kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
          kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
          km_ic(jc,1,jb)      = km_ic(jc,2,jb)
          km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
        END DO
        !$ACC END PARALLEL
      END DO
!$OMP END PARALLEL DO

  ELSE

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
#ifdef __LOOP_EXCHANGE
      DO jc = i_startidx, i_endidx
        DO jk = 2 , nlev
#else
      DO jk = 2 , nlev
        DO jc = i_startidx, i_endidx
#endif
          kh_ic(jc,jk,jb) = rho_ic(jc,jk,jb) * rturb_prandtl *                           &
                            p_nh_metrics%mixing_length_sq(jc,jk,jb) *                    &
                            SQRT( MAX( 0._wp, mech_prod(jc,jk,jb) * 0.5_wp -             &
                            rturb_prandtl * bruvais(jc,jk,jb) ) )
          km_ic(jc,jk,jb) = kh_ic(jc,jk,jb) * turb_prandtl
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
        kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
        km_ic(jc,1,jb)      = km_ic(jc,2,jb)
        km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

END IF

    !$ACC WAIT

    CALL sync_patch_array(SYNC_C, patch, kh_ic)
    CALL sync_patch_array(SYNC_C, patch, km_ic)


    !--------------------------------------------------------------------------
    !4) Interpolate difusivity (viscosity) to different locations: calculate them for
    !   halos also because they will be used later in diffusion
    !--------------------------------------------------------------------------

    !4a) visc at cell center
    rl_start = grf_bdywidth_c
    rl_end   = min_rlcell_int-1
    i_startblk = patch%cells%start_block(rl_start)
    i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                          i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1 , nlev
        DO jc = i_startidx, i_endidx
          km_c(jc,jk,jb) = MAX( km_min,                                   &
                                ( kh_ic(jc,jk,jb) + kh_ic(jc,jk+1,jb) ) * &
                                0.5_wp * turb_prandtl )
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    !4b) visc at vertices
    CALL cells2verts_scalar(kh_ic, patch, p_int%cells_aw_verts, km_iv, &
                            opt_rlstart=5, opt_rlend=min_rlvert_int-1,   &
                            opt_acc_async=.TRUE.)

    jb_max=SIZE(km_iv, 3)
    jk_max=SIZE(km_iv, 2)
    jc_max=SIZE(km_iv, 1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jb = 1, jb_max
      DO jk = 1, jk_max
        DO jc = 1, jc_max
          km_iv(jc,jk,jb) = MAX( km_min,  km_iv(jc,jk,jb) * turb_prandtl )
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

!    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
!    km_iv = MAX( km_min, km_iv * turb_prandtl )
!    !$ACC END KERNELS

    !4c) Now calculate visc at half levels at edge
    CALL cells2edges_scalar(kh_ic, patch, p_int%c_lin_e, km_ie,                   &
                            opt_rlstart=grf_bdywidth_e, opt_rlend=min_rledge_int-1, &
                            lacc=.TRUE.)

    jb_max=SIZE(km_ie, 3)
    jk_max=SIZE(km_ie, 2)
    jc_max=SIZE(km_ie, 1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO jb = 1, jb_max
      DO jk = 1, jk_max
        DO jc = 1, jc_max
          km_ie(jc,jk,jb) = MAX( km_min,  km_ie(jc,jk,jb) * turb_prandtl )
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

!    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
!    km_ie = MAX( km_min, km_ie * turb_prandtl )
!    !$ACC END KERNELS

    END ASSOCIATE

  END SUBROUTINE compute_exchange_coefficient
  !
  !=================================================================
  !
  SUBROUTINE test(inputs)

    CLASS(t_variable_set), INTENT(in), TARGET :: inputs

    TYPE(t_vdf_atmo_inputs), POINTER :: ins

    SELECT TYPE (inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins => inputs
    END SELECT

    ! ins%kh = 6._wp

  END SUBROUTINE
  !
  !=================================================================
  !
END MODULE mo_vdf_atmo
