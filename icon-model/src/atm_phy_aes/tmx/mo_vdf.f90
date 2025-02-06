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

! Classes and functions for the turbulent mixing package (tmx)

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_vdf

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_fortran_tools,     ONLY: init, copy, insert_dimension
  USE mo_util_string,       ONLY: int2string, real2string
  USE mtime,                ONLY: t_datetime => datetime
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_domain
  USE mo_vdf_atmo,          ONLY: t_vdf_atmo, t_vdf_atmo_config, t_vdf_atmo_inputs, &
    &                             t_vdf_atmo_diagnostics, prepare_diffusion_matrix
  USE mo_vdf_sfc,           ONLY: t_vdf_sfc, t_vdf_sfc_config, t_vdf_sfc_inputs, t_vdf_sfc_diagnostics
  USE mo_tmx_numerics,      ONLY: t_time_scheme_explicit_euler, &
    &                             diffuse_vertical_explicit, diffuse_vertical_implicit
  USE mo_math_utilities,    ONLY: tdma_solver
  USE mo_nonhydro_types,    ONLY: t_nh_metrics
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_intp_data_strc,    ONLY: t_int_state, p_int_state
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
    &                             sync_patch_array_mult
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c
  USE mo_intp              ,ONLY: cells2edges_scalar
  USE mo_intp_rbf          ,ONLY: rbf_vec_interpol_edge, rbf_vec_interpol_cell

  ! Todo: refactor so that t_patch is not needed
  USE mo_model_domain      ,ONLY: t_patch

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_vdf, new_vdf
  PUBLIC :: heat_type, momentum_type

  ! Turbulence scheme
  ENUM, BIND(C)
    ENUMERATOR :: tte=1, smagorinsky, tke
  END ENUM

  ! Types for diffused fields
  ENUM, BIND(C)
    ENUMERATOR :: heat_type, momentum_type
  END ENUM

  TYPE, EXTENDS(t_tmx_process) :: t_vdf
    ! TYPE(t_vdf_atmo), ALLOCATABLE :: atmo
    ! TYPE(t_vdf_sfc),  ALLOCATABLE :: sfc
    TYPE(t_vdf_atmo), POINTER :: atmo
    TYPE(t_vdf_sfc),  POINTER :: sfc
    CHARACTER(len=:), ALLOCATABLE :: sfc_coupling_type ! explicit or implicit
    INTEGER                       :: turbulence_scheme ! Turbulence scheme
  CONTAINS
    PROCEDURE :: Init => Init_vdf
    PROCEDURE :: Compute !< Compute
    PROCEDURE :: Compute_diagnostics
    PROCEDURE :: Update_diagnostics
  END TYPE t_vdf
  
  ! INTERFACE t_vdf
  !   MODULE PROCEDURE t_vdf_construct
  ! END INTERFACE t_vdf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf'

CONTAINS
  !
  !============================================================================
  !
  SUBROUTINE Init_vdf(this)
    CLASS(t_vdf), INTENT(inout), TARGET :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init'

    !CALL message(routine, '')

  END SUBROUTINE Init_vdf
  !
  !============================================================================
  !
  FUNCTION new_vdf(patch, nproma, nlev, nsfc_tiles, sfc_types, dt, &
    &              sfc_coupling_type, turbulence_scheme) &
    &      RESULT(result)

    TYPE(t_patch), POINTER :: patch
    INTEGER, INTENT(in) :: &
      & nproma, &
      & nlev, &
      & nsfc_tiles, &
      & sfc_types(:)
    REAL(wp), INTENT(in) :: &
      & dt
    CHARACTER(len=*), INTENT(in), OPTIONAL :: &
      & sfc_coupling_type,                    &
      & turbulence_scheme
    TYPE(t_vdf), POINTER :: result

    TYPE(t_vdf_atmo), POINTER :: atmo

    TYPE(t_time_scheme_explicit_euler) :: time_scheme_explicit_euler

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_vdf_construct'

    CALL message(routine, 'Construct new vdf process: vdf')
    CALL message(routine, '  Time step = '//real2string(dt))

    ALLOCATE(t_vdf::result)

    IF (PRESENT(sfc_coupling_type)) THEN
      result%sfc_coupling_type = sfc_coupling_type
    ELSE
      result%sfc_coupling_type = 'explicit'
    END IF
    IF (result%sfc_coupling_type /= 'explicit') CALL finish(routine, 'Only explicit srf-atmo coupling supported')

    IF (PRESENT(turbulence_scheme)) THEN
      SELECT CASE (turbulence_scheme)
      CASE ('tke')
        result%turbulence_scheme = tke
      CASE ('tte')
        result%turbulence_scheme = tke
      CASE ('smag')
        result%turbulence_scheme = smagorinsky
      CASE DEFAULT
        CALL finish(routine, 'Unknown turbulence scheme '//turbulence_scheme)
      END SELECT
    ELSE
      result%turbulence_scheme = smagorinsky
    END IF
    CALL message(routine, 'Turbulence scheme = '//int2string(result%turbulence_scheme))
    IF (result%turbulence_scheme /= smagorinsky) CALL finish(routine, 'Only smagorinsky turbulence scheme supported')

    !$ACC ENTER DATA COPYIN(result)

    CALL result%Init_process(dt=dt, name='vdf', domain=t_domain(patch, nproma, nlev=nlev))

    ! Sub-process for atmosphere
    result%atmo => t_vdf_atmo('vdf atmo', dt=dt, domain=t_domain(patch, nproma, nlev=nlev))
    CALL result%atmo%Set_time_scheme(time_scheme_explicit_euler)
    __acc_attach(result%atmo)

    ! Sub-process for surface
    result%sfc  => t_vdf_sfc ('vdf sfc',  dt=dt, domain=t_domain(patch, nproma, nlev=1, ntiles=nsfc_tiles, sfc_types=sfc_types))
    __acc_attach(result%sfc)

    CALL result%Add_process(result%atmo)
    CALL result%Add_process(result%sfc)

  END FUNCTION new_vdf
  !
  !============================================================================
  !
  SUBROUTINE Compute(this, datetime)

    CLASS(t_vdf), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in),   POINTER :: datetime     !< date and time at beginning of time step

    INTEGER :: jg
    TYPE(t_vdf_atmo_config),      POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_config),       POINTER :: conf_sfc
    TYPE(t_vdf_sfc_inputs),       POINTER :: ins_sfc
    TYPE(t_vdf_sfc_diagnostics),  POINTER :: diags_sfc

    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state
    TYPE(t_patch),      POINTER :: patch

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute'

    SELECT TYPE (v => this%atmo%config)
    TYPE IS (t_vdf_atmo_config)
      conf_atmo => v
    END SELECT
    __acc_attach(conf_atmo)
    SELECT TYPE (v => this%atmo%inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins_atmo => v
    END SELECT
    __acc_attach(ins_atmo)
    SELECT TYPE (v => this%atmo%diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags_atmo => v
    END SELECT
    __acc_attach(diags_atmo)

    SELECT TYPE (v => this%sfc%config)
    TYPE IS (t_vdf_sfc_config)
      conf_sfc => v
    END SELECT
    __acc_attach(conf_sfc)
    SELECT TYPE (v => this%sfc%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins_sfc => v
    END SELECT
    __acc_attach(ins_sfc)
    SELECT TYPE (v => this%sfc%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags_sfc => v
    END SELECT
    __acc_attach(diags_sfc)

    patch => this%atmo%domain%patch

    jg = patch%id
    p_int => p_int_state(jg)
    p_nh_metrics => p_nh_state(jg)%metrics

    !----------------------------------------------------------------------------
    ! Compute diagnostics at start of time step (e.g. sfc and atmo exchange coefficients)
    ! CALL this%atmo%Compute_diagnostics()
    ! Possibly put needed variables at lowest atmo level into sfc inputs collection
    !----------------------------------------------------------------------------
    CALL this%Compute_diagnostics(datetime)

    !----------------------------------------------------------------------------
    ! Call surface model and compute fluxes (so far, only explicit land/atmo is used!)
    !----------------------------------------------------------------------------
    CALL this%sfc%Compute(datetime)

    !----------------------------------------------------------------------------
    ! Call diffusioin of hydrometeors
    !----------------------------------------------------------------------------
    CALL Compute_diffusion_hydrometeors(patch,p_int,this%atmo%domain,      &
                                        this%atmo,conf_atmo,ins_atmo,      &
                                        diags_atmo,diags_sfc)

    !----------------------------------------------------------------------------
    ! Call diffusioin of temperature
    !
    ! Note: uses the new state of moisture variables
    ! Note: new_state_ta and tend_ta will be updated later by 
    !       horizontal diffusion and by additional heating 
    !----------------------------------------------------------------------------
    CALL Compute_diffusion_temperature(patch,p_int,this%atmo%domain,      &
                                       this%atmo,conf_atmo,ins_atmo,      &
                                       diags_atmo,diags_sfc)

    !----------------------------------------------------------------------------
    ! Call diffusion of horizontal wind
    !----------------------------------------------------------------------------
    CALL Compute_diffusion_hor_wind(patch,p_int,p_nh_metrics,this%atmo%domain, &
                                   this%atmo,conf_atmo,ins_atmo,              &
                                   diags_atmo,diags_sfc)

    !----------------------------------------------------------------------------
    ! Call diffusion of vertical wind
    !----------------------------------------------------------------------------
    CALL Compute_diffusion_vert_wind(patch,p_int,p_nh_metrics,this%atmo%domain, &
                                    this%atmo,conf_atmo,ins_atmo, diags_atmo)

    !----------------------------------------------------------------------------
    ! Update energy/temperature tendencies
    !----------------------------------------------------------------------------
    CALL Update_energy_tendencies(patch,this%atmo%domain,this%atmo,conf_atmo,& 
                                  ins_atmo,diags_atmo,diags_sfc)

    !----------------------------------------------------------------------------
    ! Update diagnostics at end of time step
    !----------------------------------------------------------------------------
    CALL this%Update_diagnostics()


  END SUBROUTINE Compute
  !
  !============================================================================
  !
  SUBROUTINE Compute_diagnostics(this, datetime)

    CLASS(t_vdf), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in), POINTER :: datetime

    INTEGER :: iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute_diagnostics'

    ! CALL message(routine, '')

    DO iproc=1,SIZE(this%processes)

      CALL this%processes(iproc)%p%Compute_diagnostics(datetime)

    END DO

  END SUBROUTINE Compute_diagnostics
  !
  !============================================================================
  !
  SUBROUTINE Update_diagnostics(this)

    USE mo_tmx_surface_interface, ONLY: compute_2m_temperature, compute_2m_humidity_and_dewpoint, compute_10m_wind
    USE mo_vdf_sfc,               ONLY: average_tiles

    CLASS(t_vdf), INTENT(inout), TARGET :: this

    TYPE(t_vdf_atmo_config),      POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_config),       POINTER :: conf_sfc
    TYPE(t_vdf_sfc_inputs),       POINTER :: ins_sfc
    TYPE(t_vdf_sfc_diagnostics),  POINTER :: diags_sfc

    INTEGER :: iproc, jtile, nlev, nlevm1
    INTEGER :: jc, jb, jk
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      & new_ta, new_qv, new_qc, new_qi, new_ua, new_va, &
      & new_tsfc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Update_diagnostics'

    ! CALL message(routine, '')

    DO iproc=1,SIZE(this%processes)

      CALL this%processes(iproc)%p%Update_diagnostics()

    END DO

    SELECT TYPE (v => this%atmo%config)
    TYPE IS (t_vdf_atmo_config)
      conf_atmo => v
    END SELECT
    __acc_attach(conf_atmo)
    SELECT TYPE (v => this%atmo%inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins_atmo => v
    END SELECT
    __acc_attach(ins_atmo)
    SELECT TYPE (v => this%atmo%diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags_atmo => v
    END SELECT
    __acc_attach(diags_atmo)

    SELECT TYPE (v => this%sfc%config)
    TYPE IS (t_vdf_sfc_config)
      conf_sfc => v
    END SELECT
    __acc_attach(conf_sfc)
    SELECT TYPE (v => this%sfc%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins_sfc => v
    END SELECT
    __acc_attach(ins_sfc)
    SELECT TYPE (v => this%sfc%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags_sfc => v
    END SELECT
    __acc_attach(diags_sfc)

    ASSOCIATE( &
      domain => this%atmo%domain,    &
      domain_sfc => this%sfc%domain, &
      km     => diags_atmo%km,       &
      kh     => diags_atmo%kh,       &
      km_ic  => diags_atmo%km_ic,    &
      kh_ic  => diags_atmo%kh_ic,    &
      km_sfc => diags_sfc%km,        &
      kh_sfc => diags_sfc%kh         &
      & )

    new_tsfc => this%sfc %new_states%Get_ptr_r3d('surface temperature')
    new_ta   => this%atmo%new_states%Get_ptr_r3d('temperature')
    new_qv   => this%atmo%new_states%Get_ptr_r3d('water vapor')
    new_qc   => this%atmo%new_states%Get_ptr_r3d('cloud water')
    new_qi   => this%atmo%new_states%Get_ptr_r3d('cloud ice')
    new_ua   => this%atmo%new_states%Get_ptr_r3d('eastward wind')
    new_va   => this%atmo%new_states%Get_ptr_r3d('northward wind')

    nlev = domain%nlev
    nlevm1 = nlev - 1

    DO jtile = 1, domain_sfc%ntiles

      CALL compute_10m_wind( &
        & domain_sfc, domain_sfc%sfc_types(jtile), diags_sfc%nvalid(:,jtile), diags_sfc%indices(:,:,jtile), &
        & ins_sfc%zf(:,:), ins_sfc%zh(:,:), &
        & new_ua(:,nlev,:), new_va(:,nlev,:), ins_sfc%u_oce_current(:,:), ins_sfc%v_oce_current(:,:), &
        & diags_sfc%moist_rich_tile(:,:,jtile), diags_sfc%km_tile(:,:,jtile), diags_sfc%km_neutral_tile(:,:,jtile), &
        & diags_sfc%u10m_tile(:,:,jtile), diags_sfc%v10m_tile(:,:,jtile), diags_sfc%wind10m_tile(:,:,jtile) &
        & )
    
      CALL compute_2m_temperature( &
          & domain_sfc, domain_sfc%sfc_types(jtile), diags_sfc%nvalid(:,jtile), diags_sfc%indices(:,:,jtile), &
          & ins_sfc%zf(:,:), ins_sfc%zh(:,:), new_ta(:,nlev,:), new_tsfc(:,:,jtile), &
          & diags_sfc%moist_rich_tile(:,:,jtile), diags_sfc%kh_tile(:,:,jtile), diags_sfc%km_tile(:,:,jtile), &
          & diags_sfc%kh_neutral_tile(:,:,jtile), diags_sfc%km_neutral_tile(:,:,jtile), &
          & diags_sfc%t2m_tile(:,:,jtile) &
          & )

      CALL compute_2m_humidity_and_dewpoint( &
        & domain_sfc, diags_sfc%nvalid(:,jtile), diags_sfc%indices(:,:,jtile), &
        & ins_sfc%pa(:,:), ins_sfc%psfc(:,:), &
        & new_ta(:,nlev,:), diags_sfc%t2m_tile(:,:,jtile), &
        & new_qv(:,nlev,:), new_qc(:,nlev,:), new_qi(:,nlev,:), &
        & diags_sfc%hus2m_tile(:,:,jtile), &
        & diags_sfc%dew2m_tile(:,:,jtile) &
        & )

    END DO

    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%t2m_tile,   diags_sfc%t2m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%hus2m_tile, diags_sfc%hus2m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%dew2m_tile, diags_sfc%dew2m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%u10m_tile,  diags_sfc%u10m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%v10m_tile,  diags_sfc%v10m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, &
      & diags_sfc%wind10m_tile, diags_sfc%wind10m)

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = domain%i_startblk_c, domain%i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlevm1
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          km(jc,jk,jb) = km_ic(jc,jk+1,jb)
          kh(jc,jk,jb) = kh_ic(jc,jk+1,jb)
        END DO
      END DO
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
        km(jc,nlev,jb) = km_sfc(jc,jb)
        kh(jc,nlev,jb) = kh_sfc(jc,jb)
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    END ASSOCIATE

  END SUBROUTINE Update_diagnostics
  !
  !============================================================================
  !
  SUBROUTINE Compute_diffusion_hydrometeors(patch,p_int,domain,         &
                                            atmo,conf_atmo,ins_atmo,    &
                                            diags_atmo,diags_sfc)

    TYPE(t_patch),      INTENT(in), POINTER :: patch
    TYPE(t_int_state),  INTENT(in), POINTER :: p_int         !< interpolation state
    TYPE(t_domain),     INTENT(in), POINTER :: domain

    TYPE(t_vdf_atmo),             INTENT(in), POINTER :: atmo
    TYPE(t_vdf_atmo_config),      INTENT(in), POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      INTENT(in), POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), INTENT(in), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_diagnostics),  INTENT(in), POINTER :: diags_sfc

    INTEGER :: jb, jc, je, jk, itrac
    INTEGER :: nproma, nlev, nblks_c, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk, ividx, ivblk

    REAL(wp) :: rdtime
    REAL(wp) :: sfc_flx(domain%nproma,domain%nblks_c), top_flx(domain%nproma,domain%nblks_c)
    REAL(wp) :: nabla2_e(domain%nproma,atmo%domain%nlev,domain%nblks_e)

    REAL(wp) :: &
      hori_tend_c(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      a(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      b(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      c(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      rhs(domain%nproma,atmo%domain%nlev,domain%nblks_c)

    REAL(wp) :: inv_mair(domain%nproma,atmo%domain%nlev,domain%nblks_c)
    REAL(wp), POINTER :: state(:,:,:), tend(:,:,:), new_state(:,:,:)
  
    iecidx => patch%edges%cell_idx;     iecblk => patch%edges%cell_blk
    ieidx  => patch%cells%edge_idx;     ieblk  => patch%cells%edge_blk

    nproma  = domain%nproma
    nlev    = domain%nlev
    nblks_c = domain%nblks_c

    !$ACC DATA &
    !$ACC   CREATE(sfc_flx, top_flx) &
    !$ACC   CREATE(nabla2_e, hori_tend_c) &
    !$ACC   CREATE(a, b, c, rhs, inv_mair)

    ASSOCIATE ( &
      i_startblk_c => domain%i_startblk_c,    &
      i_endblk_c   => domain%i_endblk_c,      &
      i_startidx_c => domain%i_startidx_c(:), &
      i_endidx_c   => domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,      &
      solver_type  => conf_atmo%solver_type,&
      rturb_prandtl=> conf_atmo%rturb_prandtl,&
      kh_ic        => diags_atmo%kh_ic,        &
      km_ie        => diags_atmo%km_ie,        &
      evapotrans   => diags_sfc%evapotrans, &
      rho_ic    => diags_atmo%rho_ic,       &
      zf        => ins_atmo%zf,             &
      mair      => ins_atmo%mair,           &
      rho       => ins_atmo%rho,            &
      inv_dzh   => ins_atmo%inv_dzh         &
    )

!$OMP PARALLEL
    CALL init(top_flx, lacc=.TRUE.)
!$OMP END PARALLEL

    rdtime = 1._wp / dtime

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          inv_mair(jc,jk,jb) = 1._wp / mair(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      ! In comparison to the implicit version, no matrix
      ! operations are required in the explicit version.
      ! However, the same coefficients as for the implicit
      ! version are used in the explicit version to reduce
      ! code duplication.
      CALL prepare_diffusion_matrix(                &
        & ics=i_startidx_c(jb), ice=i_endidx_c(jb), & ! in
        & minlvl=1, maxlvl=nlev,                    & ! in
        & lhalflvl=.FALSE.,                         & ! in
        & inv_mair=inv_mair(:,:,jb),                & ! in
        & inv_dz=inv_dzh(:,:,jb),                   & ! in
        & zk=kh_ic(:,:,jb),                         & ! in
        & a=a(:,:,jb),                              & ! out
        & b=b(:,:,jb),                              & ! out
        & c=c(:,:,jb)                               & ! out
        & )

    END DO
!$OMP END PARALLEL DO

    DO itrac=1,3
      SELECT CASE(itrac)
      CASE (1)
        state => atmo%states%Get_ptr_r3d('water vapor')
        tend => atmo%tendencies%Get_ptr_r3d('water vapor')
        new_state => atmo%new_states%Get_ptr_r3d('water vapor')
!$OMP PARALLEL
        CALL copy(evapotrans, sfc_flx, lacc=.TRUE.)
!$OMP END PARALLEL
      CASE (2)
        state => atmo%states%Get_ptr_r3d('cloud water')
        tend => atmo%tendencies%Get_ptr_r3d('cloud water')
        new_state => atmo%new_states%Get_ptr_r3d('cloud water')
!$OMP PARALLEL
        CALL init(sfc_flx, lacc=.TRUE.)
!$OMP END PARALLEL
      CASE (3)
        state => atmo%states%Get_ptr_r3d('cloud ice')
        tend => atmo%tendencies%Get_ptr_r3d('cloud ice')
        new_state => atmo%new_states%Get_ptr_r3d('cloud ice')
!$OMP PARALLEL
        CALL init(sfc_flx, lacc=.TRUE.)
!$OMP END PARALLEL
      END SELECT

!$OMP PARALLEL
      CALL init(tend, lacc=.TRUE.)
      CALL init(rhs, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jc) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk_c,i_endblk_c

        ! Set the right hand side
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          rhs(jc,nlev,jb) = - sfc_flx(jc,jb) * inv_mair(jc,nlev,jb)
          rhs(jc,1   ,jb) = + top_flx(jc,jb) * inv_mair(jc,1   ,jb)
        END DO
        !$ACC END PARALLEL LOOP

      END DO
!$OMP END PARALLEL DO

      IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

!$OMP PARALLEL DO PRIVATE(jb) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk_c,i_endblk_c

          ! Compute the tendencies
          CALL diffuse_vertical_explicit( &
            & ics=i_startidx_c(jb),       & ! in
            & ice=i_endidx_c(jb),         & ! in
            & minlvl=1, maxlvl=nlev,      & ! in
            & a=a(:,:,jb),                & ! in
            & b=b(:,:,jb),                & ! in
            & c=c(:,:,jb),                & ! in
            & rhs=rhs(:,:,jb),            & ! in
            & var=state(:,:,jb),          & ! in
            & tend=tend(:,:,jb)           & ! inout
            & )

        END DO
!$OMP END PARALLEL DO

      ELSE !Implicit solver

!$OMP PARALLEL DO PRIVATE(jb) ICON_OMP_RUNTIME_SCHEDULE
        DO jb = i_startblk_c,i_endblk_c

          ! Compute the tendencies
          CALL diffuse_vertical_implicit( &
            & ics=i_startidx_c(jb),       & ! in
            & ice=i_endidx_c(jb),         & ! in
            & minlvl=1, maxlvl=nlev,      & ! in
            & a=a(:,:,jb),                & ! in
            & c=c(:,:,jb),                & ! in
            & bb=b(:,:,jb),               & ! in
            & rhs=rhs(:,:,jb),            & ! in
            & rdtime=rdtime,              & ! in
            & var=state(:,:,jb),          & ! in
            & tend=tend(:,:,jb)           & ! inout
            & )

        END DO
!$OMP END PARALLEL DO

      END IF

      !---------------------------------------------------------------
      ! Horizontal diffusion (conservative; following mo_nh_diffusion)
      !---------------------------------------------------------------

      !include halo points and boundary points because these values will be
      !used in next loop
      CALL sync_patch_array(SYNC_C, patch, state)

!$OMP PARALLEL PRIVATE(rl_start, rl_end, i_startblk, i_endblk)

      rl_start   = grf_bdywidth_e
      rl_end     = min_rledge_int-1
      i_startblk = patch%edges%start_block(rl_start)
      i_endblk   = patch%edges%end_block(rl_end)

!$OMP DO PRIVATE(jk,je,jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                          i_startidx, i_endidx, rl_start, rl_end)

        ! compute kh_ie * grad_horiz(state)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) COLLAPSE(2)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            nabla2_e(je,jk,jb) = 0.5_wp * rturb_prandtl * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) &
                                * patch%edges%inv_dual_edge_length(je,jb)                         &
                                * (  state(iecidx(je,jb,2),jk,iecblk(je,jb,2))                    &
                                   - state(iecidx(je,jb,1),jk,iecblk(je,jb,1))                    &
                                  )
          ENDDO
        ENDDO
      !$ACC END PARALLEL LOOP
      ENDDO
!$OMP END DO

    ! now compute the divergence of the quantity above
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk_c,i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx_c(jb), i_endidx_c(jb)
            ! horizontal tendency
            hori_tend_c(jc,jk,jb) = (  nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) &
                                     + nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) &
                                     + nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) &
                                    ) / rho(jc,jk,jb)
            tend(jc,jk,jb) = tend(jc,jk,jb) + hori_tend_c(jc,jk,jb)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk_c,i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1) COLLAPSE(2)
        DO jk = 1, nlev
          DO jc = i_startidx_c(jb), i_endidx_c(jb)
            new_state(jc,jk,jb) = state(jc,jk,jb) + tend(jc,jk,jb) * dtime
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

      NULLIFY(state)
      NULLIFY(tend)
      NULLIFY(new_state)

    END DO ! itrac

    END ASSOCIATE

    !$ACC END DATA

  END SUBROUTINE Compute_diffusion_hydrometeors
  !
  !============================================================================
  !
  SUBROUTINE Compute_diffusion_temperature(patch,p_int,domain,         &
                                           atmo,conf_atmo,ins_atmo,    &
                                           diags_atmo,diags_sfc)

    TYPE(t_patch),      INTENT(in), POINTER :: patch
    TYPE(t_int_state),  INTENT(in), POINTER :: p_int         !< interpolation state
    TYPE(t_domain),     INTENT(in), POINTER :: domain

    TYPE(t_vdf_atmo),             INTENT(in), POINTER :: atmo
    TYPE(t_vdf_atmo_config),      INTENT(in), POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      INTENT(in), POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), INTENT(in), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_diagnostics),  INTENT(in), POINTER :: diags_sfc

    INTEGER :: jb, jc, je, jk
    INTEGER :: nproma, nlev, nblks_c, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk

    REAL(wp) :: rdtime
    REAL(wp) :: sfc_flx(domain%nproma,domain%nblks_c), top_flx(domain%nproma,domain%nblks_c)
    REAL(wp) :: nabla2_e(domain%nproma,atmo%domain%nlev,domain%nblks_e)

    REAL(wp) :: &
      inv_mair   (domain%nproma,atmo%domain%nlev,domain%nblks_c), &
      energy     (domain%nproma,atmo%domain%nlev,domain%nblks_c), &
      new_energy (domain%nproma,atmo%domain%nlev,domain%nblks_c), &
      tend_energy(domain%nproma,atmo%domain%nlev,domain%nblks_c)

    REAL(wp) :: &
      hori_tend_c(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      a(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      b(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      c(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      rhs(domain%nproma,atmo%domain%nlev,domain%nblks_c)

    REAL(wp), POINTER :: state_ta(:,:,:), tend_ta(:,:,:), new_state_ta(:,:,:)

    iecidx => patch%edges%cell_idx;     iecblk => patch%edges%cell_blk
    ieidx  => patch%cells%edge_idx;     ieblk  => patch%cells%edge_blk

    nlev    = domain%nlev

    !$ACC DATA &
    !$ACC   CREATE(sfc_flx, top_flx) &
    !$ACC   CREATE(energy, new_energy, tend_energy) &
    !$ACC   CREATE(nabla2_e, hori_tend_c) &
    !$ACC   CREATE(a, b, c, rhs, inv_mair)

    ASSOCIATE ( &
      i_startblk_c => domain%i_startblk_c,    &
      i_endblk_c   => domain%i_endblk_c,      &
      i_startidx_c => domain%i_startidx_c(:), &
      i_endidx_c   => domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,      &
      solver_type  => conf_atmo%solver_type,&
      rturb_prandtl=> conf_atmo%rturb_prandtl,&
      kh_ic        => diags_atmo%kh_ic,        &
      km_ie        => diags_atmo%km_ie,        &
      heating      => diags_atmo%heating,   &
      shfl         => diags_sfc%shfl,          &
      ufts         => diags_sfc%ufts, &
      ufvs         => diags_sfc%ufvs, &
      q_snocpymlt  => diags_sfc%q_snocpymlt_lnd, &
      rho_ic       => diags_atmo%rho_ic,       &
      zf        => ins_atmo%zf,             &
      mair      => ins_atmo%mair,           &
      rho       => ins_atmo%rho,            &
      inv_dzh   => ins_atmo%inv_dzh         &
    )

    rdtime = 1._wp / dtime

    state_ta => atmo%states%Get_ptr_r3d('temperature')
    tend_ta  => atmo%tendencies%Get_ptr_r3d('temperature')
    new_state_ta => atmo%new_states%Get_ptr_r3d('temperature')

!$OMP PARALLEL
    CALL init(top_flx, lacc=.TRUE.)
    CALL init(tend_energy, lacc=.TRUE.)
    CALL init(rhs, lacc=.TRUE.)
!$OMP END PARALLEL

    CALL atmo%temp_to_energy(state_ta(:,:,:), energy(:,:,:), use_new_moisture_state=.FALSE.)

    ! sfc_flx(:,:) = shfl(:,:)
    CALL atmo%compute_flux_x(shfl(:,:), ufts(:,:), ufvs(:,:), sfc_flx(:,:))

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c

      !$ACC PARALLEL DEFAULT(PRESENT)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          inv_mair(jc,jk,jb) = 1._wp / mair(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      ! In comparison to the implicit version, no matrix
      ! operations are required in the explicit version.
      ! However, the same coefficients as for the implicit
      ! version are used in the explicit version to reduce
      ! code duplication.
      CALL prepare_diffusion_matrix(                &
        & ics=i_startidx_c(jb), ice=i_endidx_c(jb), & ! in
        & minlvl=1, maxlvl=nlev,                    & ! in
        & lhalflvl=.FALSE.,                         & ! in
        & inv_mair=inv_mair(:,:,jb),                & ! in
        & inv_dz=inv_dzh(:,:,jb),                   & ! in
        & zk=kh_ic(:,:,jb),                         & ! in
        & a=a(:,:,jb),                              & ! out
        & b=b(:,:,jb),                              & ! out
        & c=c(:,:,jb)                               & ! out
        & )

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
      ! Set the right hand side
      DO jc = i_startidx_c(jb), i_endidx_c(jb)
        rhs(jc,nlev,jb) = - sfc_flx(jc,jb) * inv_mair(jc,nlev,jb)
        rhs(jc,1   ,jb) = + top_flx(jc,jb) * inv_mair(jc,1   ,jb)
      END DO
      !$ACC END PARALLEL LOOP

    END DO
!$OMP END PARALLEL DO

    IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

!$OMP PARALLEL DO PRIVATE(jb) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk_c,i_endblk_c

        ! Compute the tendencies
        CALL diffuse_vertical_explicit( &
          & ics=i_startidx_c(jb),       & ! in
          & ice=i_endidx_c(jb),         & ! in
          & minlvl=1, maxlvl=nlev,      & ! in
          & a=a(:,:,jb),                & ! in
          & b=b(:,:,jb),                & ! in
          & c=c(:,:,jb),                & ! in
          & rhs=rhs(:,:,jb),            & ! in
          & var=energy(:,:,jb),         & ! in
          & tend=tend_energy(:,:,jb)    & ! inout
          & )

      END DO
!$OMP END PARALLEL DO

    ELSE !Implicit solver

!$OMP PARALLEL DO PRIVATE(jb) ICON_OMP_RUNTIME_SCHEDULE
      DO jb = i_startblk_c,i_endblk_c

        ! Compute the tendencies
        CALL diffuse_vertical_implicit( &
          & ics=i_startidx_c(jb),       & ! in
          & ice=i_endidx_c(jb),         & ! in
          & minlvl=1, maxlvl=nlev,      & ! in
          & a=a(:,:,jb),                & ! in
          & c=c(:,:,jb),                & ! in
          & bb=b(:,:,jb),               & ! in
          & rhs=rhs(:,:,jb),            & ! in
          & rdtime=rdtime,              & ! in
          & var=energy(:,:,jb),         & ! in
          & tend=tend_energy(:,:,jb)    & ! inout
          & )

      END DO
!$OMP END PARALLEL DO

    END IF

    !---------------------------------------------------------------
    ! Horizontal diffusion (conservative; following mo_nh_diffusion)
    !---------------------------------------------------------------

    !include halo points and boundary points because these values will be
    !used in next loop
    CALL sync_patch_array(SYNC_C, patch, energy)

    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL

!$OMP DO PRIVATE(jk, je, jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                        i_startidx, i_endidx, rl_start, rl_end)

      ! compute kh_ie * grad_horiz(energy)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          nabla2_e(je,jk,jb) = 0.5_wp * rturb_prandtl * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) &
                              * patch%edges%inv_dual_edge_length(je,jb)                         &
                              * (  energy(iecidx(je,jb,2),jk,iecblk(je,jb,2))                   &
                                 - energy(iecidx(je,jb,1),jk,iecblk(je,jb,1))                   &
                                )
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDDO
!$OMP END DO

    ! now compute the divergence of the quantity above
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c, i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          ! horizontal tendency
          hori_tend_c(jc,jk,jb) = (  nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) &
                                   + nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) &
                                   + nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) &
                                  ) / rho(jc,jk,jb)
          tend_energy(jc,jk,jb) = tend_energy(jc,jk,jb) + hori_tend_c(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1, nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          new_energy(jc,jk,jb) = energy(jc,jk,jb) + tend_energy(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    CALL atmo%energy_to_temp(new_energy(:,:,:), new_state_ta(:,:,:), use_new_moisture_state=.TRUE.)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1, nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          tend_ta(jc,jk,jb) = (new_state_ta(jc,jk,jb) - state_ta(jc,jk,jb)) * rdtime
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

    END ASSOCIATE

    !$ACC END DATA

  END SUBROUTINE Compute_diffusion_temperature
  !
  !============================================================================
  !
  SUBROUTINE Compute_diffusion_hor_wind(patch,p_int,p_nh_metrics,domain, &
                                           atmo,conf_atmo,ins_atmo,      &
                                           diags_atmo,diags_sfc)

    TYPE(t_patch),      INTENT(in), POINTER :: patch
    TYPE(t_int_state),  INTENT(in), POINTER :: p_int         !< interpolation state
    TYPE(t_nh_metrics), INTENT(in) ,POINTER :: p_nh_metrics
    TYPE(t_domain),     INTENT(in), POINTER :: domain

    TYPE(t_vdf_atmo),             INTENT(in), POINTER :: atmo
    TYPE(t_vdf_atmo_config),      INTENT(in), POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      INTENT(in), POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), INTENT(in), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_diagnostics),  INTENT(in), POINTER :: diags_sfc

    INTEGER :: jb, jc, je, jk, jcn, jbn, jvn
    INTEGER :: nproma, nlev, nblks_c, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ividx, ivblk

    REAL(wp) :: rdtime
    REAL(wp) :: flux_up_v, flux_dn_v, flux_up_c, flux_dn_c, flux_dn_e, stress_c1n, stress_c2n, dwdn
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt

    REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp
    REAL(wp),         PARAMETER :: z_2by3  = 2._wp/3._wp
    !$ACC DECLARE COPYIN(z_1by3, z_2by3)

    REAL(wp) :: &
      za(domain%nproma,atmo%domain%nlev,domain%nblks_e),  &
      zb(domain%nproma,atmo%domain%nlev,domain%nblks_e),  &
      zc(domain%nproma,atmo%domain%nlev,domain%nblks_e),  &
      zrhs(domain%nproma,atmo%domain%nlev,domain%nblks_e)


    ! horizontal diffusion for horizontal wind
    REAL(wp) :: &
      inv_rhoe( domain%nproma,atmo%domain%nlev,domain%nblks_e),    &
      inv_maire(domain%nproma,atmo%domain%nlev,domain%nblks_e),    &
      tot_tend( domain%nproma,atmo%domain%nlev,domain%nblks_e)

    REAL(wp), POINTER ::                                      &
      & state_u(:,:,:),  tend_u(:,:,:),  new_state_u(:,:,:),  &
      & state_v(:,:,:),  tend_v(:,:,:),  new_state_v(:,:,:)

    iecidx => patch%edges%cell_idx;     iecblk => patch%edges%cell_blk
    ividx  => patch%edges%vertex_idx;   ivblk  => patch%edges%vertex_blk

    nlev    = domain%nlev

    state_u     => atmo%states%Get_ptr_r3d('eastward wind')
    tend_u      => atmo%tendencies%Get_ptr_r3d('eastward wind')
    new_state_u => atmo%new_states%Get_ptr_r3d('eastward wind')
    state_v     => atmo%states%Get_ptr_r3d('northward wind')
    tend_v      => atmo%tendencies%Get_ptr_r3d('northward wind')
    new_state_v => atmo%new_states%Get_ptr_r3d('northward wind')

!$OMP PARALLEL
    CALL init(tend_u, lacc=.TRUE.)
    CALL init(tend_v, lacc=.TRUE.)
!$OMP END PARALLEL

    !$ACC DATA &
    !$ACC   CREATE(inv_rhoe, tot_tend) &
    !$ACC   CREATE(za, zb, zc, zrhs, inv_maire)

    ASSOCIATE ( &
      i_startblk_c => domain%i_startblk_c,    &
      i_endblk_c   => domain%i_endblk_c,      &
      i_startidx_c => domain%i_startidx_c(:), &
      i_endidx_c   => domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,             &
      solver_type  => conf_atmo%solver_type,       &
      rturb_prandtl=> conf_atmo%rturb_prandtl,&
      km_c         => diags_atmo%km_c,         &
      km_iv        => diags_atmo%km_iv,        &
      km_ic        => diags_atmo%km_ic,        &
      km_ie        => diags_atmo%km_ie,        &
      mflux_u      => diags_sfc%ustress,       &
      mflux_v      => diags_sfc%vstress,       &
      rho_ic       => diags_atmo%rho_ic,       &
      u_vert       => diags_atmo%u_vert,       &
      v_vert       => diags_atmo%v_vert,       &
      vn           => diags_atmo%vn,           &
      pwp1         => ins_atmo%pwp1,           &
      div_c        => diags_atmo%div_c,        &
      zf           => ins_atmo%zf,             &
      mair         => ins_atmo%mair,           &
      rho          => ins_atmo%rho,            &
      inv_dzh      => ins_atmo%inv_dzh,        &
      dissip_kin_energy => diags_atmo%dissip_kin_energy &
    )

    rdtime = 1._wp / dtime

    !---------------------------------------------------------------
    ! Horizontal diffusion for horizontal wind
    !---------------------------------------------------------------

!$OMP PARALLEL
    CALL init(tot_tend, lacc=.TRUE.)
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_C, patch, rho)

    !density at edge
    CALL cells2edges_scalar(rho, patch, p_int%c_lin_e, inv_rhoe,                  &
                            opt_rlstart=grf_bdywidth_e+1, opt_rlend=min_rledge_int, &
                            lacc=.TRUE.)

    rl_start   = grf_bdywidth_e+1
    rl_end     = min_rledge_int
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          inv_rhoe(je,jk,jb) = 1._wp / inv_rhoe(je,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    ! 1) First get the horizontal tendencies

!$OMP PARALLEL DO PRIVATE(jb,jk,je,i_startidx,i_endidx,vn_vert1,vn_vert2,vn_vert3,vn_vert4,&
!$OMP            dvt,jcn,jbn,flux_up_c,flux_dn_c,jvn,flux_up_v,flux_dn_v)  ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt, jcn, jbn, jvn) &
      !$ACC   PRIVATE(flux_up_c, flux_dn_c, flux_up_v, flux_dn_v)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          vn_vert1 =   u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * patch%edges%primal_normal_vert(je,jb,1)%v1 &
                     + v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1)) * patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 =   u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * patch%edges%primal_normal_vert(je,jb,2)%v1 &
                     + v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2)) * patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 =   u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * patch%edges%primal_normal_vert(je,jb,3)%v1 &
                     + v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 =   u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * patch%edges%primal_normal_vert(je,jb,4)%v1 &
                     + v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * patch%edges%primal_normal_vert(je,jb,4)%v2

          dvt      =   u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * patch%edges%dual_normal_vert(je,jb,4)%v1 &
                     + v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) * patch%edges%dual_normal_vert(je,jb,4)%v2 &
                     - u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * patch%edges%dual_normal_vert(je,jb,3)%v1 &   
                     - v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) * patch%edges%dual_normal_vert(je,jb,3)%v2 
                      

          ! tendency in normal direction:
          ! flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

          jcn       = iecidx(je,jb,2)
          jbn       = iecblk(je,jb,2)
          flux_up_c = km_c(jcn,jk,jbn)                                                                  &
                      * ( 4._wp * ( vn_vert4 - vn(je,jk,jb) ) * patch%edges%inv_vert_vert_length(je,jb) &
                          - z_2by3 * div_c(jcn,jk,jbn) )


          jcn       = iecidx(je,jb,1)
          jbn       = iecblk(je,jb,1)
          flux_dn_c = km_c(jcn,jk,jbn)                                                                  &
                      * ( 4._wp * ( vn(je,jk,jb) - vn_vert3 ) * patch%edges%inv_vert_vert_length(je,jb) &
                          - z_2by3 * div_c(jcn,jk,jbn) )

          ! tendency in tangential direction

          ! D_12 between edge center and the vertex: delta_v/(primal_edge_len/2) +
          ! ((vt4+vt2)/2-(vt3+vt2)/2)/(distance_opp_edges)
          ! flux = D_12 * visc

          ! Note that the tangential velocities at vertices are used in D_12 is an
          ! approximation for speed. Better way is to use vt reconstructed from vn at
          ! each edges. Also, visc at somewhere between edge mid point and the vertex
          ! should be used but this is a good approximation

          jvn       = ividx(je,jb,2)
          jbn       = ivblk(je,jb,2)
          flux_up_v = ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) )                              &
                      * ( patch%edges%tangent_orientation(je,jb) * ( vn_vert2 - vn(je,jk,jb) ) &    
                          * patch%edges%inv_primal_edge_length(je,jb)                          &
                          + 0.5_wp * dvt * patch%edges%inv_vert_vert_length(je,jb)             &
                        )

          jvn       = ividx(je,jb,1)
          jbn       = ivblk(je,jb,1)
          flux_dn_v = ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) )                              &
                      * ( patch%edges%tangent_orientation(je,jb) * ( vn(je,jk,jb) - vn_vert1 ) &    
                          * patch%edges%inv_primal_edge_length(je,jb)                          &
                          + 0.5_wp * dvt * patch%edges%inv_vert_vert_length(je,jb)             &
                        )

          tot_tend(je,jk,jb) = ( ( flux_up_c - flux_dn_c ) * patch%edges%inv_dual_edge_length(je,jb)          &
                                 + 2._wp * patch%edges%tangent_orientation(je,jb) * ( flux_up_v - flux_dn_v ) & 
                                   * patch%edges%inv_primal_edge_length(je,jb)                                &
                               ) * inv_rhoe(je,jk,jb)

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    ! 2) Vertical tendency

!$OMP PARALLEL
    CALL init(zrhs, lacc=.TRUE.)
!$OMP END PARALLEL

    ! Sync momentum fluxes (otherwise, MPI test fails)
    ! TODO: Can't we just use sync_patch_array directly?
    CALL sync_uvml_s(mflux_u, mflux_v, patch)

!$OMP PARALLEL DO PRIVATE(jb,jk,je,i_startidx,i_endidx, stress_c1n, stress_c2n, flux_dn_e) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 2, nlev-1
        DO je = i_startidx, i_endidx
          inv_maire(je,jk,jb) = p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)
          zrhs(je,jk,jb) = (  km_ie(je,jk,jb) * patch%edges%inv_dual_edge_length(je,jb)                                    & 
                              * ( pwp1(iecidx(je,jb,2),jk,iecblk(je,jb,2)) - pwp1(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )    &
                            - km_ie(je,jk+1,jb) * patch%edges%inv_dual_edge_length(je,jb)                                  &
                              * ( pwp1(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) - pwp1(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) )&
                           ) * inv_maire(je,jk,jb)
        END DO
      END DO

      !jk = 1
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO je = i_startidx, i_endidx
        inv_maire(je,1,jb) = inv_rhoe(je,1,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,1,jb)
        zrhs(je,1,jb) = ( - km_ie(je,2,jb) * patch%edges%inv_dual_edge_length(je,jb)                              &
                          * ( pwp1(iecidx(je,jb,2),2,iecblk(je,jb,2)) - pwp1(iecidx(je,jb,1),2,iecblk(je,jb,1)) ) &
                        ) * inv_maire(je,1,jb)
      END DO

      ! jk = nlev
      !$ACC LOOP GANG(STATIC: 1) VECTOR &
      !$ACC   PRIVATE(dwdn, stress_c1n, stress_c2n, flux_dn_e)
      DO je = i_startidx, i_endidx
        inv_maire(je,nlev,jb) = inv_rhoe(je,nlev,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb)

        ! term due to dwdn- goes to RHS
        dwdn = km_ie(je,nlev,jb) * patch%edges%inv_dual_edge_length(je,jb)                                   &
               * ( pwp1(iecidx(je,jb,2),nlev,iecblk(je,jb,2)) - pwp1(iecidx(je,jb,1),nlev,iecblk(je,jb,1)) ) &
               * inv_maire(je,nlev,jb)

        ! Get net shear stress in the direction of vn at surface
        ! shear stress in normal direction from cell 1
        stress_c1n =   mflux_u(iecidx(je,jb,1),iecblk(je,jb,1)) * patch%edges%primal_normal_cell(je,jb,1)%v1 &
                     + mflux_v(iecidx(je,jb,1),iecblk(je,jb,1)) * patch%edges%primal_normal_cell(je,jb,1)%v2

        ! shear stress in normal direction from cell 2
        stress_c2n =   mflux_u(iecidx(je,jb,2),iecblk(je,jb,2)) * patch%edges%primal_normal_cell(je,jb,2)%v1 &
                     + mflux_v(iecidx(je,jb,2),iecblk(je,jb,2)) * patch%edges%primal_normal_cell(je,jb,2)%v2

        ! Net stress at the edge
        flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + stress_c2n * p_int%c_lin_e(je,2,jb)

        zrhs(je,nlev,jb) = dwdn - flux_dn_e * inv_maire(je,nlev,jb)
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)

      ! In comparison to the implicit version, no matrix
      ! operations are required in the explicit version.
      ! However, the same coefficients as for the implicit
      ! version are used in the explicit version to reduce
      ! code duplication.
      CALL prepare_diffusion_matrix(                      &
        & ics=i_startidx, ice=i_endidx,                   & ! in
        & minlvl=1, maxlvl=nlev,                          & ! in
        & lhalflvl=.FALSE.,                               & ! in
        & inv_mair=inv_maire(:,:,jb),                     & ! in
        & inv_dz=p_nh_metrics%inv_ddqz_z_half_e(:,:,jb),  & ! in
        & zk=km_ie(:,:,jb),                               & ! in
        & a=za(:,:,jb),                                   & ! out
        & b=zb(:,:,jb),                                   & ! out
        & c=zc(:,:,jb)                                    & ! out
        & )

    END DO
!$OMP END PARALLEL DO

    IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)

      ! Compute the tendencies
      CALL diffuse_vertical_explicit( &
        & ics=i_startidx,             & ! in
        & ice=i_endidx,               & ! in
        & minlvl=1, maxlvl=nlev,      & ! in
        & a=za(:,:,jb),               & ! in
        & b=zb(:,:,jb),               & ! in
        & c=zc(:,:,jb),               & ! in
        & rhs=zrhs(:,:,jb),           & ! in
        & var=vn(:,:,jb),             & ! in
        & tend=tot_tend(:,:,jb)       & ! inout
        & )

    END DO
!$OMP END PARALLEL DO

    ELSE !Implicit solver

!$OMP PARALLEL DO PRIVATE(jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                            i_startidx, i_endidx, rl_start, rl_end)

        ! Compute the tendencies
        CALL diffuse_vertical_implicit( &
          & ics=i_startidx,             & ! in
          & ice=i_endidx,               & ! in
          & minlvl=1, maxlvl=nlev,      & ! in
          & a=za(:,:,jb),               & ! in
          & c=zc(:,:,jb),               & ! in
          & bb=zb(:,:,jb),              & ! in
          & rhs=zrhs(:,:,jb),           & ! in
          & rdtime=rdtime,              & ! in
          & var=vn(:,:,jb),             & ! in
          & tend=tot_tend(:,:,jb)       & ! inout
          & )

      END DO !jb
!$OMP END PARALLEL DO

    END IF

    CALL sync_patch_array(SYNC_E, patch, tot_tend)
    CALL rbf_vec_interpol_cell(tot_tend, patch, p_int, tend_u, tend_v, opt_rlend=min_rlcell_int)

!$OMP PARALLEL DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          new_state_u(jc,jk,jb) = state_u(jc,jk,jb) + tend_u(jc,jk,jb) * dtime
          new_state_v(jc,jk,jb) = state_v(jc,jk,jb) + tend_v(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    ! TODO: Are these necessary?
    CALL sync_patch_array_mult(SYNC_C, patch, 2, tend_u, tend_v)
    CALL sync_patch_array_mult(SYNC_C, patch, 2, new_state_u, new_state_v)

    END ASSOCIATE

    !$ACC END DATA
  END SUBROUTINE Compute_diffusion_hor_wind
  !
  !============================================================================
  !
  SUBROUTINE Compute_diffusion_vert_wind(patch,p_int,p_nh_metrics,domain, &
                                         atmo,conf_atmo,ins_atmo,         &
                                         diags_atmo)

    TYPE(t_patch),      INTENT(in), POINTER :: patch
    TYPE(t_int_state),  INTENT(in), POINTER :: p_int         !< interpolation state
    TYPE(t_nh_metrics), INTENT(in) ,POINTER :: p_nh_metrics
    TYPE(t_domain),     INTENT(in), POINTER :: domain

    TYPE(t_vdf_atmo),             INTENT(in), POINTER :: atmo
    TYPE(t_vdf_atmo_config),      INTENT(in), POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      INTENT(in), POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), INTENT(in), POINTER :: diags_atmo

    INTEGER :: jb, jc, je, jk, jcn, jbn, jvn
    INTEGER :: nlev, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk, ividx, ivblk

    REAL(wp) :: rdtime
    REAL(wp) :: flux_up_v, flux_dn_v, flux_up_c, flux_dn_c
    REAL(wp) :: dvn1, dvn2, dvt1, dvt2
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt

    REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp
    !$ACC DECLARE COPYIN(z_1by3)

    REAL(wp) :: &
      hori_tend_e(domain%nproma,atmo%domain%nlev,domain%nblks_e),   &
      a(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      b(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      c(domain%nproma,atmo%domain%nlev,domain%nblks_c),   &
      rhs(domain%nproma,atmo%domain%nlev,domain%nblks_c)
    
    REAL(wp) :: &
      vt_e   (domain%nproma,atmo%domain%nlev,domain%nblks_e),      &
      inv_mair_ic(domain%nproma,atmo%domain%nlev,domain%nblks_c),  &
      inv_rho_ic(domain%nproma,atmo%domain%nlev,domain%nblks_c)

    REAL(wp),                     POINTER :: &
      & state(:,:,:), tend(:,:,:), new_state(:,:,:)

    iecidx => patch%edges%cell_idx;     iecblk => patch%edges%cell_blk
    ividx  => patch%edges%vertex_idx;   ivblk  => patch%edges%vertex_blk
    ieidx  => patch%cells%edge_idx;     ieblk  => patch%cells%edge_blk

    nlev    = domain%nlev
    
    state => atmo%states%Get_ptr_r3d('vertical velocity')
    tend  => atmo%tendencies%Get_ptr_r3d('vertical velocity')
    new_state => atmo%new_states%Get_ptr_r3d('vertical velocity')

    !$ACC DATA &
    !$ACC   CREATE(hori_tend_e) &
    !$ACC   CREATE(a, b, c, rhs, vt_e, inv_rho_ic, inv_mair_ic)

    ASSOCIATE ( &
      i_startblk_c => domain%i_startblk_c,    &
      i_endblk_c   => domain%i_endblk_c,      &
      i_startidx_c => domain%i_startidx_c(:), &
      i_endidx_c   => domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,      &
      rturb_prandtl=> conf_atmo%rturb_prandtl,&
      km_ie     => diags_atmo%km_ie,        &
      km_c      => diags_atmo%km_c,         &
      km_ic     => diags_atmo%km_ic,        &
      km_iv     => diags_atmo%km_iv,        &
      w_ie      => diags_atmo%w_ie,         &
      rho_ic    => diags_atmo%rho_ic,       &
      u_vert    => diags_atmo%u_vert,       &
      v_vert    => diags_atmo%v_vert,       &
      w_vert    => diags_atmo%w_vert,       &
      vn        => diags_atmo%vn,           &
      div_c     => diags_atmo%div_c,        &
      pum1      => ins_atmo%pum1,           &
      pvm1      => ins_atmo%pvm1,           &
      pwp1      => ins_atmo%pwp1,           &
      zf        => ins_atmo%zf,             &
      rho       => ins_atmo%rho,            &
      inv_dzf   => ins_atmo%inv_dzf,        &
      inv_dzh   => ins_atmo%inv_dzh         &
      )

      rdtime = 1._wp / dtime

!$OMP PARALLEL
    CALL init(tend, lacc=.TRUE.)
    CALL init(new_state, lacc=.TRUE.)
    CALL init(rhs, lacc=.TRUE.)
!$OMP END PARALLEL

    !---------------------------------------------------------------
    ! Vertical diffusion for w-wind
    !---------------------------------------------------------------
    CALL rbf_vec_interpol_edge( vn, patch, p_int, vt_e, opt_rlend=min_rledge_int-1)

!$OMP PARALLEL DO PRIVATE(jb,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          inv_rho_ic(jc,jk,jb) = 1._wp / rho_ic(jc,jk,jb)
          inv_mair_ic(jc,jk,jb) = inv_rho_ic(jc,jk,jb) * inv_dzh(jc,jk,jb)
          rhs(jc,jk,jb) = 2._wp * inv_mair_ic(jc,jk,jb)                       &
                          * (   km_c(jc,jk,jb)   * z_1by3 * div_c(jc,jk,jb)   &
                              - km_c(jc,jk-1,jb) * z_1by3 * div_c(jc,jk-1,jb) & 
                            ) 
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC WAIT(1)

      ! Compute the coefficients of the matrix
      CALL prepare_diffusion_matrix(                &
        & ics=i_startidx_c(jb), ice=i_endidx_c(jb), & ! in
        & minlvl=2, maxlvl=nlev,                    & ! in
        & lhalflvl=.TRUE.,                          & ! in
        & inv_mair=inv_mair_ic(:,:,jb),             & ! in
        & inv_dz=inv_dzf(:,:,jb),                   & ! in
        & zk=km_c(:,:,jb),                          & ! in
        & zprefac=2._wp,                            & ! in
        & a=a(:,:,jb),                              & ! out
        & b=b(:,:,jb),                              & ! out
        & c=c(:,:,jb)                               & ! out
        & )

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR
      DO jc = i_startidx_c(jb), i_endidx_c(jb)     
        ! This results from the condition w=0 at the top and bottom boundary.
        b(jc,2,jb)    = b(jc,2,jb)    + 2._wp * km_c(jc,1,jb)    * inv_dzf(jc,1,jb)    * inv_mair_ic(jc,2,jb)
        b(jc,nlev,jb) = b(jc,nlev,jb) + 2._wp * km_c(jc,nlev,jb) * inv_dzf(jc,nlev,jb) * inv_mair_ic(jc,nlev,jb)
      END DO
      !$ACC END PARALLEL LOOP

    END DO

!$OMP END PARALLEL DO

    ! 2) Vertical tendency: evaluated at w point

!$OMP PARALLEL DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c

      ! Compute the tendencies
      CALL diffuse_vertical_implicit( &
        & ics=i_startidx_c(jb),       & ! in
        & ice=i_endidx_c(jb),         & ! in
        & minlvl=2, maxlvl=nlev,      & ! in
        & a=a(:,:,jb),                & ! in
        & c=c(:,:,jb),                & ! in
        & bb=b(:,:,jb),               & ! in
        & rhs=rhs(:,:,jb),            & ! in
        & rdtime=rdtime,              & ! in
        & var=pwp1(:,:,jb),           & ! in
        & tend=tend(:,:,jb)           & ! inout
        & )

    END DO !jb
!$OMP END PARALLEL DO

    !---------------------------------------------------------------
    ! Horizontal diffusion for w-wind
    !---------------------------------------------------------------

    ! 1) Get horizontal tendencies at half level edges
    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL
    CALL init(hori_tend_e, lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, jcn, jbn, dvn1, dvn2, flux_up_c, flux_dn_c, &
!$OMP                     jvn, dvt1, dvt2, flux_up_v, flux_dn_v) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(jcn, jbn, dvn1, dvn2, flux_up_c, flux_dn_c, jvn, dvt1, dvt2, flux_up_v, flux_dn_v)
      DO jk = 2, nlev
        DO je = i_startidx, i_endidx

          ! tendency in normal direction
          ! flux = visc_c * D_31_c where D_31(=D_13) is calculated at half level
          ! cell center

          jcn   = iecidx(je,jb,2)
          jbn   = iecblk(je,jb,2)

          dvn2  =   pum1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,2)%v1 &
                  + pvm1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,2)%v2 &
                  - pum1(jcn,jk,jbn)   * patch%edges%primal_normal_cell(je,jb,2)%v1 &
                  - pvm1(jcn,jk,jbn)   * patch%edges%primal_normal_cell(je,jb,2)%v2

          flux_up_c = km_ic(jcn,jk,jbn)                                                   &
                      * ( dvn2 * inv_dzh(jcn,jk,jbn)                                      & 
                          + ( w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) - w_ie(je,jk,jb) ) &
                            * 2.0_wp * patch%edges%inv_vert_vert_length(je,jb)            &
                        )

          jcn   = iecidx(je,jb,1)
          jbn   = iecblk(je,jb,1)

          dvn1  =   pum1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v1 &
                  + pvm1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v2 &
                  - pum1(jcn,jk,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v1   &
                  - pvm1(jcn,jk,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v2


          flux_dn_c = km_ic(jcn,jk,jbn)                                                   &
                      * ( dvn1 * inv_dzh(jcn,jk,jbn)                                      & 
                          + ( w_ie(je,jk,jb) - w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) ) &
                            * 2.0_wp * patch%edges%inv_vert_vert_length(je,jb)            &
                        )


         ! tendency in tangential direction
         ! flux = visc_v * D_32_v where D_32(= D_23) is calculated at half level
         ! between vertex and edge center

          jvn  = ividx(je,jb,2)
          jbn  = ivblk(je,jb,2)

          dvt2 =   0.5_wp * (  u_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v1 &
                             + v_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v2 &
                             + vt_e(je,jk-1,jb) )                                              &
                 - 0.5_wp * (  u_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v1   &
                             + v_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v2   &
                             + vt_e(je,jk,jb) ) 

          flux_up_v = km_iv(jvn,jk,jbn)                                       &                    
                      * (   dvt2 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) &                    
                          + patch%edges%tangent_orientation(je,jb)            &
                            * ( w_vert(jvn,jk,jbn) - w_ie(je,jk,jb) )         &
                              / patch%edges%edge_cell_length(je,jb,2)         &
                        )

          jvn  = ividx(je,jb,1)
          jbn  = ivblk(je,jb,1)

          dvt1 =  0.5_wp * (  u_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v1  &
                            + v_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v2  &
                            + vt_e(je,jk-1,jb) )                                               &
                - 0.5_wp * (  u_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v1    &
                            + v_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v2    &
                            + vt_e(je,jk,jb) )

          flux_dn_v = km_iv(jvn,jk,jbn)                                       &
                      * (   dvt1 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) &
                          + patch%edges%tangent_orientation(je,jb)            &
                            * ( w_ie(je,jk,jb) - w_vert(jvn,jk,jbn) )         &
                              / patch%edges%edge_cell_length(je,jb,1)         &
                        )

          hori_tend_e(je,jk,jb) =   ( flux_up_c - flux_dn_c ) * patch%edges%inv_dual_edge_length(je,jb) &
                                  + ( flux_up_v - flux_dn_v ) * patch%edges%tangent_orientation(je,jb)  &
                                    * 2._wp * patch%edges%inv_primal_edge_length(je,jb) 
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    ! Interpolate horizontal tendencies to w point: except top and bottom boundaries
    ! w==0 at these boundaries
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          tend(jc,jk,jb) = tend(jc,jk,jb)                                                                  &
                           + inv_rho_ic(jc,jk,jb)                                                          &            
                             * (  hori_tend_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%e_bln_c_s(jc,1,jb) & 
                                + hori_tend_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%e_bln_c_s(jc,2,jb) & 
                                + hori_tend_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%e_bln_c_s(jc,3,jb) & 
                               )
          ! new_state(jc,jk,jb) = state(jc,jk,jb) + tend(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END DO
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          new_state(jc,jk,jb) = state(jc,jk,jb) + tend(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END DO
!$OMP END PARALLEL

    CALL sync_patch_array(SYNC_C, patch, new_state)

    END ASSOCIATE

    !$ACC END DATA
    
  END SUBROUTINE Compute_diffusion_vert_wind
  !
  !============================================================================
  !
  SUBROUTINE sync_uvml_s(u, v, patch)
    REAL(wp), TARGET, INTENT(inout) :: u(:,:), v(:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch

    REAL(wp), POINTER :: pu(:,:,:), pv(:,:,:)
 
    CALL insert_dimension(pu, u, 2)
    CALL insert_dimension(pv, v, 2)
    CALL sync_patch_array_mult(SYNC_C, patch, 2, pu, pv)
  END SUBROUTINE sync_uvml_s
  !
  !============================================================================
  !
  SUBROUTINE Update_energy_tendencies(patch,domain, &
                                      atmo,conf_atmo,ins_atmo,      &
                                      diags_atmo,diags_sfc)

    TYPE(t_patch),      INTENT(in), POINTER :: patch
    TYPE(t_domain),     INTENT(in), POINTER :: domain

    TYPE(t_vdf_atmo),             INTENT(in), POINTER :: atmo
    TYPE(t_vdf_atmo_config),      INTENT(in), POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      INTENT(in), POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), INTENT(in), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_diagnostics),  INTENT(in), POINTER :: diags_sfc

    INTEGER :: jb, jk, jc, nlev
    REAL(wp) :: rdtime

    REAL(wp), POINTER ::                                      &
      & state_ta(:,:,:), tend_ta(:,:,:), new_state_ta(:,:,:), &
      & state_u(:,:,:),  new_state_u(:,:,:),  &
      & state_v(:,:,:),  new_state_v(:,:,:)

    state_ta      => atmo%states%Get_ptr_r3d('temperature')
    tend_ta       => atmo%tendencies%Get_ptr_r3d('temperature')
    new_state_ta  => atmo%new_states%Get_ptr_r3d('temperature')
    state_u       => atmo%states%Get_ptr_r3d('eastward wind')
    new_state_u   => atmo%new_states%Get_ptr_r3d('eastward wind')
    state_v       => atmo%states%Get_ptr_r3d('northward wind')
    new_state_v   => atmo%new_states%Get_ptr_r3d('northward wind')

    nlev    = domain%nlev


    ASSOCIATE ( &
      i_startblk_c => domain%i_startblk_c,    &
      i_endblk_c   => domain%i_endblk_c,      &
      i_startidx_c => domain%i_startidx_c(:), &
      i_endidx_c   => domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,      &
      heating      => diags_atmo%heating,   &
      dissipation_factor => conf_atmo%dissipation_factor, &
      dissip_kin_energy => diags_atmo%dissip_kin_energy, &
      q_snocpymlt       => diags_sfc%q_snocpymlt_lnd, &
      mair              => ins_atmo%mair,           &
      cvair             => ins_atmo%cvair           &
    ) 

    rdtime = 1._wp / dtime

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk=1,nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          heating(jc,jk,jb) = 0._wp
        END DO
      END DO
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx_c(jb), i_endidx_c(jb)
        heating(jc,nlev,jb) = - q_snocpymlt(jc,jb) ! non-zero only for land
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          dissip_kin_energy(jc,jk,jb) = 0.5_wp * mair(jc,jk,jb) * dissipation_factor * rdtime  &
                                        * (   state_u(jc,jk,jb)**2 - new_state_u(jc,jk,jb)**2  &
                                            + state_v(jc,jk,jb)**2 - new_state_v(jc,jk,jb)**2  &
                                          ) 
          heating(jc,jk,jb)      = heating(jc,jk,jb) + dissip_kin_energy(jc,jk,jb)
          tend_ta(jc,jk,jb)      = tend_ta(jc,jk,jb) + heating(jc,jk,jb) / cvair(jc,jk,jb)
          new_state_ta(jc,jk,jb) = state_ta(jc,jk,jb) + tend_ta(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    ! TODO: Are these necessary?
    CALL sync_patch_array_mult(SYNC_C, patch, 2, new_state_ta, tend_ta)

    END ASSOCIATE

  END SUBROUTINE Update_energy_tendencies
  !
  !============================================================================
  !
END MODULE mo_vdf 
