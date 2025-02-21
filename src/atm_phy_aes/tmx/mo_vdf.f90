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

MODULE mo_vdf

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_fortran_tools,     ONLY: init, copy, insert_dimension
  USE mo_util_string,       ONLY: int2string, real2string
  USE mtime,                ONLY: t_datetime => datetime
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_field_class,   ONLY: t_tmx_field, t_domain, &
    &                             isfc_oce, isfc_ice, isfc_lnd
  USE mo_vdf_atmo,          ONLY: t_vdf_atmo, t_vdf_atmo_config, t_vdf_atmo_inputs, t_vdf_atmo_diagnostics
  USE mo_vdf_sfc,           ONLY: t_vdf_sfc, t_vdf_sfc_config, t_vdf_sfc_inputs, t_vdf_sfc_diagnostics
  USE mo_tmx_numerics,      ONLY: t_time_scheme_explicit_euler, &
    &                             diffuse_scalar_vertical_explicit, diffuse_scalar_vertical_implicit
  USE mo_math_utilities,    ONLY: tdma_solver, tdma_solver_vec
  USE mo_nonhydro_types,    ONLY: t_nh_metrics
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_intp_data_strc,    ONLY: t_int_state, p_int_state
  USE mo_sync,              ONLY: SYNC_E, SYNC_C, SYNC_V, sync_patch_array, &
    &                             sync_patch_array_mult
  USE mo_impl_constants,    ONLY: min_rlcell, min_rledge_int, min_rlcell_int, min_rlvert_int
  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_loopindices,       ONLY: get_indices_e, get_indices_c
  USE mo_intp              ,ONLY: cells2verts_scalar, cells2edges_scalar
  USE mo_intp_rbf          ,ONLY: rbf_vec_interpol_edge, rbf_vec_interpol_cell

  ! Todo: refactor so that t_patch is not needed
  USE mo_model_domain      ,ONLY: t_patch


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
    PROCEDURE :: Compute !< Compute
    PROCEDURE :: Compute_diagnostics
    PROCEDURE :: Update_diagnostics
  END TYPE t_vdf
  
  ! INTERFACE t_vdf
  !   MODULE PROCEDURE t_vdf_construct
  ! END INTERFACE t_vdf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_vdf'

CONTAINS

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

    CALL result%Init(dt=dt, name='vdf', domain=t_domain(patch, nproma, nlev=nlev))

    ! Sub-process for atmosphere
    result%atmo => t_vdf_atmo('vdf atmo', dt=dt, domain=t_domain(patch, nproma, nlev=nlev))
    CALL result%atmo%Set_time_scheme(time_scheme_explicit_euler)
    

    ! Sub-process for surface
    result%sfc  => t_vdf_sfc ('vdf sfc',  dt=dt, domain=t_domain(patch, nproma, nlev=1, ntiles=nsfc_tiles, sfc_types=sfc_types))
    

    CALL result%Add_process(result%atmo)
    CALL result%Add_process(result%sfc)

  END FUNCTION new_vdf

  SUBROUTINE Compute(this, datetime)

    CLASS(t_vdf), INTENT(inout), TARGET :: this
    TYPE(t_datetime), OPTIONAL, INTENT(in),   POINTER :: datetime     !< date and time at beginning of time step

    INTEGER :: jg, jb, jc, je, jk
    TYPE(t_vdf_atmo_config),      POINTER :: conf_atmo
    TYPE(t_vdf_atmo_inputs),      POINTER :: ins_atmo
    TYPE(t_vdf_atmo_diagnostics), POINTER :: diags_atmo
    TYPE(t_vdf_sfc_config),       POINTER :: conf_sfc
    TYPE(t_vdf_sfc_inputs),       POINTER :: ins_sfc
    TYPE(t_vdf_sfc_diagnostics),  POINTER :: diags_sfc

    TYPE(t_nh_metrics) ,POINTER :: p_nh_metrics
    TYPE(t_int_state)  ,POINTER :: p_int         !< interpolation state
    TYPE(t_patch),      POINTER :: patch

    INTEGER :: nproma, nlev, nblks_c
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end

    INTEGER,  DIMENSION(:,:,:), POINTER :: iecidx, iecblk, ieidx, ieblk, ividx, ivblk

    REAL(wp) :: sfc_flx(this%domain%nproma,this%domain%nblks_c), top_flx(this%domain%nproma,this%domain%nblks_c)
    REAL(wp) :: nabla2_e(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e)

    REAL(wp) :: &
      hori_tend_c(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c),   &
      hori_tend_e(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),   &
      a(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c),   &
      b(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c),   &
      c(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c),   &
      rhs(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c), &
      za(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),   &
      zb(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),   &
      zc(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),   &
      zrhs(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e)


    REAL(wp),                     POINTER :: &
      & state(:,:,:), tend(:,:,:), new_state(:,:,:), &
      & state_ta(:,:,:), tend_ta(:,:,:), new_state_ta(:,:,:), &
      & state_u(:,:,:), tend_u(:,:,:), new_state_u(:,:,:), state_v(:,:,:), tend_v(:,:,:), new_state_v(:,:,:)

    REAL(wp) :: &
      energy     (this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c), &
      new_energy (this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c), &
      tend_energy(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c)

    ! horizontal diffusion for horizontal wind
    REAL(wp) :: &
      inv_rhoe(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),     &
      tot_tend(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e)
      ! ddt_u   (this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c),     &
      ! ddt_v   (this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c)

    ! vertical diffusion for w-wind
    REAL(wp) :: &
      var_new(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c),      &
      var_new_e(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),      &
      vt_e   (this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_e),      &
      inv_rho_ic(this%domain%nproma,this%atmo%domain%nlev,this%domain%nblks_c)

    REAL(wp) :: dvn1, dvn2, dvt1, dvt2

    REAL(wp) :: flux_up_v, flux_dn_v, flux_up_c, flux_dn_c, flux_up_e, flux_dn_e, stress_c1n, stress_c2n, dwdn
    REAL(wp) :: vn_vert1, vn_vert2, vn_vert3, vn_vert4, dvt
    INTEGER :: jcn, jbn, jvn
    INTEGER :: itrac

    REAL(wp),         PARAMETER :: z_1by3  = 1._wp/3._wp
    !$ACC DECLARE COPYIN(z_1by3)
    !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Compute'

    ! CALL message(routine, '')

    nproma  = this%domain%nproma
    nlev    = this%domain%nlev
    nblks_c = this%domain%nblks_c

    !$ACC DATA &
    !$ACC   CREATE(sfc_flx, top_flx) &
    !$ACC   CREATE(energy, new_energy, tend_energy) &
    !$noACC   CREATE(inv_rhoe, nabla2_e, tot_tend, ddt_u, ddt_v, hori_tend_c, hori_tend_e) &
    !$ACC   CREATE(inv_rhoe, nabla2_e, tot_tend, hori_tend_c, hori_tend_e) &
    !$ACC   CREATE(a, b, c, rhs, var_new, vt_e, inv_rho_ic)

    SELECT TYPE (v => this%atmo%config)
    TYPE IS (t_vdf_atmo_config)
      conf_atmo => v
    END SELECT
    
    SELECT TYPE (v => this%atmo%inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins_atmo => v
    END SELECT
    
    SELECT TYPE (v => this%atmo%diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags_atmo => v
    END SELECT
    

    SELECT TYPE (v => this%sfc%config)
    TYPE IS (t_vdf_sfc_config)
      conf_sfc => v
    END SELECT
    
    SELECT TYPE (v => this%sfc%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins_sfc => v
    END SELECT
    
    SELECT TYPE (v => this%sfc%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags_sfc => v
    END SELECT
    

    conf_sfc%cpd => conf_atmo%cpd
    conf_sfc%cvd => conf_atmo%cvd

    ! Compute diagnostics at start of time step (e.g. sfc and atmo exchange coefficients)
    ! CALL this%atmo%Compute_diagnostics()
    ! Possibly put needed variables at lowest atmo level into sfc inputs collection
    ! ...
    ! CALL this%sfc%Compute_diagnostics()
    CALL this%Compute_diagnostics(datetime)

    ! Call surface model and compute fluxes (so far, only explicit land/atmo is used!)
    CALL this%sfc%Compute(datetime)
  
    ! Pseudo code:
    ! allocate (tri-diagonal) matrix
    ! for state in atmo%states:
    !   compute atmo part of tri-diagonal matrix
    !   get boundary values from surface
    !   compute sfc part of tri-diagnoal matrix
    !   solve and compute tendency of state
    ! end

    ! IF ( conf_atmo%solver_type == 1 ) THEN
    !   CALL message(routine, 'solver is explicit')
    ! ELSE
    !   CALL message(routine, 'solver is implicit')
    ! END IF

    patch => this%atmo%domain%patch

    jg = patch%id
    p_int => p_int_state(jg)
    p_nh_metrics => p_nh_state(jg)%metrics

    iecidx => patch%edges%cell_idx
    iecblk => patch%edges%cell_blk

    ieidx => patch%cells%edge_idx
    ieblk => patch%cells%edge_blk

    ividx  => patch%edges%vertex_idx
    ivblk  => patch%edges%vertex_blk

    ! ---------------------------------------
    ! 1. Vertical and horizontal diffusion of hydrometeors (iqv, iqc, ici)
    ! ---------------------------------------
!$OMP PARALLEL
    CALL init(top_flx)
!$OMP END PARALLEL

    ASSOCIATE ( &
      i_startblk_c => this%domain%i_startblk_c,    &
      i_endblk_c   => this%domain%i_endblk_c,      &
      i_startidx_c => this%domain%i_startidx_c(:), &
      i_endidx_c   => this%domain%i_endidx_c(:),   &
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
      inv_dzf   => ins_atmo%inv_dzf         &
      )

    DO itrac=1,3
      SELECT CASE(itrac)
      CASE (1)
        state => this%atmo%states%Get_ptr_r3d('water vapor')
        tend => this%atmo%tendencies%Get_ptr_r3d('water vapor')
        new_state => this%atmo%new_states%Get_ptr_r3d('water vapor')
!$OMP PARALLEL
        CALL copy(evapotrans, sfc_flx)
!$OMP END PARALLEL
      CASE (2)
        state => this%atmo%states%Get_ptr_r3d('cloud water')
        tend => this%atmo%tendencies%Get_ptr_r3d('cloud water')
        new_state => this%atmo%new_states%Get_ptr_r3d('cloud water')
!$OMP PARALLEL
        CALL init(sfc_flx)
!$OMP END PARALLEL
      CASE (3)
        state => this%atmo%states%Get_ptr_r3d('cloud ice')
        tend => this%atmo%tendencies%Get_ptr_r3d('cloud ice')
        new_state => this%atmo%new_states%Get_ptr_r3d('cloud ice')
!$OMP PARALLEL
        CALL init(sfc_flx)
!$OMP END PARALLEL
      END SELECT

!$OMP PARALLEL
      CALL init(tend)
      CALL init(new_state)
!$OMP END PARALLEL

      IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

        CALL diffuse_scalar_vertical_explicit( &
          & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
          & mair(:,:,:), zf(:,:,:), &
          & rho_ic(:,:,:), kh_ic(:,:,:), &
          & state(:,:,:), sfc_flx(:,:), top_flx(:,:), &
          & tend(:,:,:) &
          & )

      ELSE !Implicit solver

        CALL diffuse_scalar_vertical_implicit( &
          & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
          & dtime, &
          & mair(:,:,:), zf(:,:,:), &
          & rho_ic(:,:,:), kh_ic(:,:,:), &
          & state(:,:,:), sfc_flx(:,:), top_flx(:,:), &
          & tend(:,:,:) &
          & )

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
            nabla2_e(je,jk,jb) = 0.5_wp * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) *                    &
                                  rturb_prandtl *                                                       &
                                  patch%edges%inv_dual_edge_length(je,jb) *                           &
                                  ( state(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -                           &
                                    state(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )
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
            hori_tend_c(jc,jk,jb) = (                                                  &
              nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) +  &
              nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) +  &
              nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) ) / rho(jc,jk,jb)
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

    ! ---------------------------------------
    ! 2. Vertical and horizontal diffusion of temperature resp. thermal energy
    !
    ! Note: uses the new state of moisture variables
    ! ---------------------------------------

    ASSOCIATE ( &
      i_startblk_c => this%domain%i_startblk_c,    &
      i_endblk_c   => this%domain%i_endblk_c,      &
      i_startidx_c => this%domain%i_startidx_c(:), &
      i_endidx_c   => this%domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,      &
      solver_type  => conf_atmo%solver_type,&
      rturb_prandtl=> conf_atmo%rturb_prandtl,&
      kh_ic        => diags_atmo%kh_ic,        &
      km_ie        => diags_atmo%km_ie,        &
      heating      => diags_atmo%heating,   &
      shfl         => diags_sfc%shfl,          &
      q_snocpymlt  => diags_sfc%q_snocpymlt_lnd, &
      rho_ic       => diags_atmo%rho_ic,       &
      zf        => ins_atmo%zf,             &
      mair      => ins_atmo%mair,           &
      rho       => ins_atmo%rho,            &
      inv_dzf   => ins_atmo%inv_dzf         &
      )

    state_ta => this%atmo%states%Get_ptr_r3d('temperature')
    tend_ta  => this%atmo%tendencies%Get_ptr_r3d('temperature')
    new_state_ta => this%atmo%new_states%Get_ptr_r3d('temperature')

!$OMP PARALLEL
    CALL init(top_flx)
    CALL init(tend_ta)
    CALL init(energy)
    CALL init(tend_energy)
    CALL init(new_energy)
!$OMP END PARALLEL

    CALL this%atmo%temp_to_energy(state_ta(:,:,:), energy(:,:,:), use_new_moisture_state=.FALSE.)

    ! sfc_flx(:,:) = shfl(:,:)
    CALL this%atmo%energy_flux_to_flux_x(shfl(:,:), sfc_flx(:,:))

    IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

      CALL diffuse_scalar_vertical_explicit( &
        & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
        & mair(:,:,:), zf(:,:,:), &
        & rho_ic(:,:,:), kh_ic(:,:,:), &
        & energy(:,:,:), sfc_flx(:,:), top_flx(:,:), &
        & tend_energy(:,:,:) &
        & )

    ELSE !Implicit solver

      CALL diffuse_scalar_vertical_implicit( &
        & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
        & dtime, &
        & mair(:,:,:), zf(:,:,:), &
        & rho_ic(:,:,:), kh_ic(:,:,:), &
        & energy(:,:,:), sfc_flx(:,:), top_flx(:,:), &
        & tend_energy(:,:,:) &
        & )

    END IF

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

    CALL this%atmo%energy_to_temp(new_energy(:,:,:), new_state_ta(:,:,:), use_new_moisture_state=.TRUE.)

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk = 1, nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          tend_ta(jc,jk,jb) = (new_state_ta(jc,jk,jb) - state_ta(jc,jk,jb)) / dtime
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

    ! Note: new_state_ta and tend_ta will be updated later by horizontal diffusion and by additional heating (see below)
    ! Additional heating
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk=1,nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          ! heating(jc,jk,jb) = tend_energy(jc,jk,jb) * mair(jc,jk,jb)
          heating(jc,jk,jb) = 0._wp
        END DO
      END DO
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx_c(jb), i_endidx_c(jb)
        ! heating(jc,nlev,jb) = heating(jc,nlev,jb) - q_snocpymlt(jc,jb)
        heating(jc,nlev,jb) = - q_snocpymlt(jc,jb) ! non-zero only for land
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

!$OMP END PARALLEL

    !---------------------------------------------------------------
    ! Horizontal diffusion (conservative; following mo_nh_diffusion)
    !---------------------------------------------------------------

    !include halo points and boundary points because these values will be
    !used in next loop
    CALL sync_patch_array(SYNC_C, patch, state_ta)

    rl_start   = grf_bdywidth_e
    rl_end     = min_rledge_int-1
    i_startblk = patch%edges%start_block(rl_start)
    i_endblk   = patch%edges%end_block(rl_end)

!$OMP PARALLEL

!$OMP DO PRIVATE(jk, je, jb, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                        i_startidx, i_endidx, rl_start, rl_end)

      ! compute kh_ie * grad_horiz(state_ta)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          nabla2_e(je,jk,jb) = 0.5_wp * ( km_ie(je,jk,jb) + km_ie(je,jk+1,jb) ) *                    &
                                rturb_prandtl *                                                       &
                                patch%edges%inv_dual_edge_length(je,jb) *                           &
                                ( state_ta(iecidx(je,jb,2),jk,iecblk(je,jb,2)) -                           &
                                  state_ta(iecidx(je,jb,1),jk,iecblk(je,jb,1)) )
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
          hori_tend_c(jc,jk,jb) = (                                                  &
            nabla2_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) * p_int%geofac_div(jc,1,jb) +  &
            nabla2_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) * p_int%geofac_div(jc,2,jb) +  &
            nabla2_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) * p_int%geofac_div(jc,3,jb) ) / rho(jc,jk,jb)
          tend_ta(jc,jk,jb) = tend_ta(jc,jk,jb) + hori_tend_c(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END DO

!$OMP END PARALLEL

    END ASSOCIATE

    ! ---------------------------------------
    ! 3. Diffusion of horizontal velocity
    ! ---------------------------------------

    state_u => this%atmo%states%Get_ptr_r3d('eastward wind')
    tend_u => this%atmo%tendencies%Get_ptr_r3d('eastward wind')
    new_state_u => this%atmo%new_states%Get_ptr_r3d('eastward wind')
    state_v => this%atmo%states%Get_ptr_r3d('northward wind')
    tend_v => this%atmo%tendencies%Get_ptr_r3d('northward wind')
    new_state_v => this%atmo%new_states%Get_ptr_r3d('northward wind')

!$OMP PARALLEL
    CALL init(tend_u)
    CALL init(new_state_u)
    CALL init(tend_v)
    CALL init(new_state_v)
!$OMP END PARALLEL

    ASSOCIATE ( &
      i_startblk_c => this%domain%i_startblk_c,    &
      i_endblk_c   => this%domain%i_endblk_c,      &
      i_startidx_c => this%domain%i_startidx_c(:), &
      i_endidx_c   => this%domain%i_endidx_c(:),   &
      dtime        => conf_atmo%dtime,      &
      solver_type  => conf_atmo%solver_type,&
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
      heating      => diags_atmo%heating,      &
      zf           => ins_atmo%zf,             &
      mair         => ins_atmo%mair,           &
      q2t_factor   => ins_atmo%q2t_factor,     &
      rho          => ins_atmo%rho,            &
      inv_dzf      => ins_atmo%inv_dzf,        &
      dissip_kin_energy => diags_atmo%dissip_kin_energy &
      )

    !---------------------------------------------------------------
    ! Vertical diffusion for u-wind
    !---------------------------------------------------------------

    ! IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

    !   CALL diffuse_scalar_vertical_explicit( &
    !     & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
    !     & mair(:,:,:), zf(:,:,:), &
    !     & rho_ic(:,:,:), km_ic(:,:,:), &
    !     & state_u(:,:,:), mflux_u(:,:), top_flx(:,:), &
    !     & tend_u(:,:,:) &
    !     & )

    ! ELSE !Implicit solver

    !   CALL diffuse_scalar_vertical_implicit( &
    !     & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
    !     & dtime, &
    !     & mair(:,:,:), zf(:,:,:), &
    !     & rho_ic(:,:,:), km_ic(:,:,:), &
    !     & state_u(:,:,:), mflux_u(:,:), top_flx(:,:), &
    !     & tend_u(:,:,:) &
    !     & )

    ! END IF

    ! !---------------------------------------------------------------
    ! ! Vertical diffusion for v-wind
    ! !---------------------------------------------------------------

    ! IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

    !   CALL diffuse_scalar_vertical_explicit( &
    !     & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
    !     & mair(:,:,:), zf(:,:,:), &
    !     & rho_ic(:,:,:), km_ic(:,:,:), &
    !     & state_v(:,:,:), mflux_v(:,:), top_flx(:,:), &
    !     & tend_v(:,:,:) &
    !     & )

    ! ELSE !Implicit solver

    !   CALL diffuse_scalar_vertical_implicit( &
    !     & i_startblk_c, i_endblk_c, i_startidx_c(:), i_endidx_c(:), &
    !     & dtime, &
    !     & mair(:,:,:), zf(:,:,:), &
    !     & rho_ic(:,:,:), km_ic(:,:,:), &
    !     & state_v(:,:,:), mflux_v(:,:), top_flx(:,:), &
    !     & tend_v(:,:,:) &
    !     & )

    ! END IF

    !---------------------------------------------------------------
    ! Horizontal diffusion for horizontal wind
    !---------------------------------------------------------------
!$OMP PARALLEL
    CALL init(tot_tend)
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
          vn_vert1 = u_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     patch%edges%primal_normal_vert(je,jb,1)%v1 +   &
                     v_vert(ividx(je,jb,1),jk,ivblk(je,jb,1))     *   &
                     patch%edges%primal_normal_vert(je,jb,1)%v2

          vn_vert2 = u_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     patch%edges%primal_normal_vert(je,jb,2)%v1 +   &
                     v_vert(ividx(je,jb,2),jk,ivblk(je,jb,2))     *   &
                     patch%edges%primal_normal_vert(je,jb,2)%v2

          vn_vert3 = u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     patch%edges%primal_normal_vert(je,jb,3)%v1 +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     patch%edges%primal_normal_vert(je,jb,3)%v2

          vn_vert4 = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     patch%edges%primal_normal_vert(je,jb,4)%v1 +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     patch%edges%primal_normal_vert(je,jb,4)%v2

          dvt      = u_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     patch%edges%dual_normal_vert(je,jb,4)%v1   +   &
                     v_vert(ividx(je,jb,4),jk,ivblk(je,jb,4))     *   &
                     patch%edges%dual_normal_vert(je,jb,4)%v2   -   &
                     (u_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))    *   &
                     patch%edges%dual_normal_vert(je,jb,3)%v1   +   &
                     v_vert(ividx(je,jb,3),jk,ivblk(je,jb,3))     *   &
                     patch%edges%dual_normal_vert(je,jb,3)%v2)

          ! tendency in normal direction:
          ! flux = visc*(D_11-2/3DIV) = visc*(2*delta_v/(vert_vert_len/2)-2/3*div_of_stress)

          jcn       = iecidx(je,jb,2)
          jbn       = iecblk(je,jb,2)
          flux_up_c = km_c(jcn,jk,jbn) * ( 4._wp * ( vn_vert4 - vn(je,jk,jb) ) *     &
                      patch%edges%inv_vert_vert_length(je,jb) - 2._wp * z_1by3 *   &
                      div_c(jcn,jk,jbn) )


          jcn       = iecidx(je,jb,1)
          jbn       = iecblk(je,jb,1)
          flux_dn_c = km_c(jcn,jk,jbn) * ( 4._wp * ( vn(je,jk,jb) - vn_vert3 ) *     &
                      patch%edges%inv_vert_vert_length(je,jb) -  2._wp * z_1by3 *  &
                      div_c(jcn,jk,jbn) )

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
          flux_up_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                      &
                      ( patch%edges%tangent_orientation(je,jb) *                 &
                        ( vn_vert2 - vn(je,jk,jb) )    *                                          &
                        patch%edges%inv_primal_edge_length(je,jb) * 2._wp +      &
                        dvt * patch%edges%inv_vert_vert_length(je,jb) )

          jvn       = ividx(je,jb,1)
          jbn       = ivblk(je,jb,1)
          flux_dn_v = 0.5_wp * ( km_iv(jvn,jk,jbn) + km_iv(jvn,jk+1,jbn) ) *                      &
                      ( patch%edges%tangent_orientation(je,jb) *                 &
                        ( vn(je,jk,jb) - vn_vert1 )    *                                          &
                        patch%edges%inv_primal_edge_length(je,jb) * 2._wp +      &
                        dvt * patch%edges%inv_vert_vert_length(je,jb) )

          tot_tend(je,jk,jb) = ( ( flux_up_c - flux_dn_c ) *                                      &
                                 patch%edges%inv_dual_edge_length(je,jb) +       &
                                 patch%edges%tangent_orientation(je,jb) *        &
                                 ( flux_up_v - flux_dn_v ) *                                      &
                                 patch%edges%inv_primal_edge_length(je,jb) * 2._wp ) * &
                               inv_rhoe(je,jk,jb)

        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    ! 2) Vertical tendency

    ! Sync momentum fluxes (otherwise, MPI test fails)
    ! TODO: Can't we just use sync_patch_array directly?
    CALL sync_uvml_s(mflux_u, mflux_v)

    IF ( SOLVER_TYPE == 1 ) THEN !Explicit solver

!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, flux_up_e,flux_dn_e, stress_c1n, stress_c2n) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(flux_up_e, flux_dn_e)
      DO jk = 2, nlev-1
        DO je = i_startidx, i_endidx
          flux_up_e = km_ie(je,jk,jb) *                                                         &
          ( ( vn(je,jk-1,jb) - vn(je,jk,jb) ) *                 &
            p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb)              +                 &
            ( pwp1(iecidx(je,jb,2),jk,iecblk(je,jb,2) )    -                 &
            pwp1(iecidx(je,jb,1),jk,iecblk(je,jb,1) ) )    *                 &
            patch%edges%inv_dual_edge_length(je,jb) )

          flux_dn_e = km_ie(je,jk+1,jb) *                                                       &
          ( ( vn(je,jk,jb) - vn(je,jk+1,jb) ) *                 &
            p_nh_metrics%inv_ddqz_z_half_e(je,jk+1,jb)            +                 &
            ( pwp1(iecidx(je,jb,2),jk+1,iecblk(je,jb,2))   -                 &
            pwp1(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) )   *                 &
            patch%edges%inv_dual_edge_length(je,jb) )

          tot_tend(je,jk,jb) = tot_tend(je,jk,jb) +  ( flux_up_e - flux_dn_e ) *                &
                   p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)
        END DO
      END DO

      ! 3) Boundary treatment in vertical

      !-----------------------------------------------------------------
      ! jk = 1
      !-----------------------------------------------------------------
      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(flux_dn_e)
      DO je = i_startidx, i_endidx
        flux_dn_e = km_ie(je,2,jb) *                                              &
                    ( ( vn(je,1,jb) - vn(je,2,jb) )  *        &
                       p_nh_metrics%inv_ddqz_z_half_e(je,2,jb)           +        &
                       ( pwp1(iecidx(je,jb,2),2,iecblk(je,jb,2) ) -        &
                       pwp1(iecidx(je,jb,1),2,iecblk(je,jb,1)) )  *        &
                       patch%edges%inv_dual_edge_length(je,jb) )

        tot_tend(je,1,jb) = tot_tend(je,1,jb) - flux_dn_e *                       &
                            p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)
      END DO
      !-----------------------------------------------------------------
      ! jk = nlev
      !-----------------------------------------------------------------
      !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
      !$ACC   PRIVATE(flux_up_e, flux_dn_e, stress_c1n, stress_c2n)
      DO je = i_startidx, i_endidx
        flux_up_e = km_ie(je,nlev,jb) *                                           &
                    ( ( vn(je,nlev-1,jb) - vn(je,nlev,jb) ) * &
                       p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb)               + &
                      ( pwp1(iecidx(je,jb,2),nlev,iecblk(je,jb,2))       - &
                      pwp1(iecidx(je,jb,1),nlev,iecblk(je,jb,1)) )       * &
                      patch%edges%inv_dual_edge_length(je,jb) )

        ! Get net shear stress in the direction of vn at surface

        ! shear stress in normal direction from cell 1
        stress_c1n = mflux_u(iecidx(je,jb,1),iecblk(je,jb,1))  *          &
          patch%edges%primal_normal_cell(je,jb,1)%v1      +          &
          mflux_v(iecidx(je,jb,1),iecblk(je,jb,1))  *          &
          patch%edges%primal_normal_cell(je,jb,1)%v2

        ! shear stress in normal direction from cell 2
        stress_c2n = mflux_u(iecidx(je,jb,2),iecblk(je,jb,2))  *          &
                      patch%edges%primal_normal_cell(je,jb,2)%v1     +          &
                      mflux_v(iecidx(je,jb,2),iecblk(je,jb,2)) *          &
                      patch%edges%primal_normal_cell(je,jb,2)%v2

        ! Net stress at the edge
        flux_dn_e  = stress_c1n * p_int%c_lin_e(je,1,jb) + stress_c2n * p_int%c_lin_e(je,2,jb)

        tot_tend(je,nlev,jb) = tot_tend(je,nlev,jb) + ( flux_up_e - flux_dn_e ) *               &
                              p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)
      END DO

      !$ACC END PARALLEL

    END DO
!$OMP END PARALLEL DO

    ELSE !Implicit solver

!$OMP PARALLEL
      CALL init(za)
      CALL init(zb)
      CALL init(zc)
      CALL init(zrhs)
      CALL init(var_new_e)
!$OMP END PARALLEL
  
!$OMP PARALLEL DO PRIVATE(jb, jk, je, i_startidx, i_endidx, dwdn, var_new,&
!$OMP                     flux_dn_e, stress_c1n, stress_c2n) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                            i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

        !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
        !$ACC   PRIVATE(dwdn)
        DO jk = 2, nlev-1
          DO je = i_startidx, i_endidx

            za(je,jk,jb)   = - km_ie(je,jk,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *           &
                            p_nh_metrics%inv_ddqz_z_half_e(je,jk,jb) * inv_rhoe(je,jk,jb)

            zc(je,jk,jb)   = - km_ie(je,jk+1,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) *         &
                            p_nh_metrics%inv_ddqz_z_half_e(je,jk+1,jb) * inv_rhoe(je,jk,jb)

            zb(je,jk,jb)   =  1._wp / dtime - za(je,jk,jb) - zc(je,jk,jb)

            ! term due to dwdn - goes to RHS
            dwdn       =  ( km_ie(je,jk,jb) * patch%edges%inv_dual_edge_length(je,jb)   *       &
                            ( pwp1(iecidx(je,jb,2),jk,iecblk(je,jb,2))   -                 &
                              pwp1(iecidx(je,jb,1),jk,iecblk(je,jb,1)) ) -                 &
                            km_ie(je,jk+1,jb) * patch%edges%inv_dual_edge_length(je,jb) *       &
                            ( pwp1(iecidx(je,jb,2),jk+1,iecblk(je,jb,2)) -                 &
                              pwp1(iecidx(je,jb,1),jk+1,iecblk(je,jb,1)) ) )       *       &
                          p_nh_metrics%inv_ddqz_z_full_e(je,jk,jb) * inv_rhoe(je,jk,jb)

            zrhs(je,jk,jb) =  vn(je,jk,jb) / dtime + dwdn
          END DO
        END DO

        ! Boundary treatment
        !--------------------------------------------------------
        ! jk = 1
        !--------------------------------------------------------
        !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
        !$ACC   PRIVATE(dwdn)
        DO je = i_startidx, i_endidx
          zc(je,1,jb)   = - km_ie(je,2,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) *                &
                        p_nh_metrics%inv_ddqz_z_half_e(je,2,jb) * inv_rhoe(je,1,jb)

          zb(je,1,jb)   =  1._wp / dtime - zc(je,1,jb)

          !term due to dwdn- goes to RHS
          dwdn      = - km_ie(je,2,jb) * patch%edges%inv_dual_edge_length(je,jb) *              &
                        ( pwp1(iecidx(je,jb,2),2,iecblk(je,jb,2))   -                      &
                          pwp1(iecidx(je,jb,1),2,iecblk(je,jb,1)) ) *                      &
                        p_nh_metrics%inv_ddqz_z_full_e(je,1,jb) * inv_rhoe(je,1,jb)

          zrhs(je,1,jb) =  vn(je,1,jb) / dtime + dwdn
        END DO

        !--------------------------------------------------------
        ! jk = nlev
        !--------------------------------------------------------

        !$ACC LOOP GANG(STATIC: 1) VECTOR TILE(32, 4) &
        !$ACC   PRIVATE(dwdn, stress_c1n, stress_c2n, flux_dn_e)
        DO je = i_startidx, i_endidx
          za(je,nlev,jb)  = - km_ie(je,nlev,jb) * p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) *        &
                          p_nh_metrics%inv_ddqz_z_half_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

          zb(je,nlev,jb)  = 1._wp / dtime - za(je,nlev,jb)

          ! term due to dwdn- goes to RHS
          dwdn        = km_ie(je,nlev,jb) * patch%edges%inv_dual_edge_length(je,jb) *           &
                        ( pwp1(iecidx(je,jb,2),nlev,iecblk(je,jb,2))   -                   &
                          pwp1(iecidx(je,jb,1),nlev,iecblk(je,jb,1)) ) *                   &
                        p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)

          ! Get net shear stress in the direction of vn at surface

          ! shear stress in normal direction from cell 1
          stress_c1n = mflux_u(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        patch%edges%primal_normal_cell(je,jb,1)%v1     + &
                        mflux_v(iecidx(je,jb,1),iecblk(je,jb,1)) * &
                        patch%edges%primal_normal_cell(je,jb,1)%v2

          ! shear stress in normal direction from cell 2
          stress_c2n = mflux_u(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        patch%edges%primal_normal_cell(je,jb,2)%v1     + &
                        mflux_v(iecidx(je,jb,2),iecblk(je,jb,2)) * &
                        patch%edges%primal_normal_cell(je,jb,2)%v2

          ! Net stress at the edge
          flux_dn_e    = stress_c1n * p_int%c_lin_e(je,1,jb) + stress_c2n * p_int%c_lin_e(je,2,jb)

          zrhs(je,nlev,jb) = vn(je,nlev,jb) / dtime + dwdn - flux_dn_e *                   &
                        p_nh_metrics%inv_ddqz_z_full_e(je,nlev,jb) * inv_rhoe(je,nlev,jb)
        END DO

        !$ACC END PARALLEL

        CALL tdma_solver_vec(za(:,:,jb),zb(:,:,jb),zc(:,:,jb),zrhs(:,:,jb),1,nlev,i_startidx,i_endidx,var_new_e(:,:,jb))

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
            tot_tend(je,jk,jb) = tot_tend(je,jk,jb) + (var_new_e(je,jk,jb) - vn(je,jk,jb)) / dtime
          END DO
        END DO
        !$ACC END PARALLEL LOOP
        
      END DO !jb
!$OMP END PARALLEL DO

    END IF

    ! 4) Update vn: it makes more sense to first apply diffusion on vn
    ! and then get ddt_u/v than to interpolate tot_tend directly to
    ! ddt_u/v. Proof: during the test phase it was found that the latter slowed
    ! down the computation by 15%, although the results look same.

!$OMP PARALLEL DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk,i_endblk
      CALL get_indices_e(patch, jb, i_startblk, i_endblk,       &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1, nlev
        DO je = i_startidx, i_endidx
          var_new_e(je,jk,jb) = vn(je,jk,jb) + dtime * tot_tend(je,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

    ! 5) Get turbulent tendency at cell center
    CALL sync_patch_array(SYNC_E, patch, var_new_e)
    CALL rbf_vec_interpol_cell(var_new_e, patch, p_int, new_state_u, new_state_v, opt_rlend=min_rlcell_int)

!$OMP PARALLEL DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          tend_u(jc,jk,jb) = ( new_state_u(jc,jk,jb) - state_u(jc,jk,jb) ) / dtime
          tend_v(jc,jk,jb) = ( new_state_v(jc,jk,jb) - state_v(jc,jk,jb) ) / dtime
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          dissip_kin_energy(jc,jk,jb) = 0.5_wp * ( state_u(jc,jk,jb)**2 - new_state_u(jc,jk,jb)**2 &
            &                                    + state_v(jc,jk,jb)**2 - new_state_v(jc,jk,jb)**2 &
            &                                    ) * mair(jc,jk,jb) / dtime
          heating(jc,jk,jb) = heating(jc,jk,jb) + dissip_kin_energy(jc,jk,jb)
          tend_ta(jc,jk,jb) = tend_ta(jc,jk,jb) + heating(jc,jk,jb) * q2t_factor(jc,jk,jb)
          new_state_ta(jc,jk,jb) = state_ta(jc,jk,jb) + tend_ta(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!$OMP END PARALLEL DO

    ! TODO: Are these necessary?
    CALL sync_patch_array_mult(SYNC_C, patch, 2, tend_u, tend_v)
    CALL sync_patch_array_mult(SYNC_C, patch, 2, new_state_u, new_state_v)
    CALL sync_patch_array_mult(SYNC_C, patch, 2, new_state_ta, tend_ta)

    END ASSOCIATE

    ! ---------------------------------------
    ! 4. Diffusion of vertical velocity
    ! ---------------------------------------

    state => this%atmo%states%Get_ptr_r3d('vertical velocity')
    tend  => this%atmo%tendencies%Get_ptr_r3d('vertical velocity')
    new_state => this%atmo%new_states%Get_ptr_r3d('vertical velocity')
    ASSOCIATE ( &
      i_startblk_c => this%domain%i_startblk_c,    &
      i_endblk_c   => this%domain%i_endblk_c,      &
      i_startidx_c => this%domain%i_startidx_c(:), &
      i_endidx_c   => this%domain%i_endidx_c(:),   &
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

!$OMP PARALLEL
    CALL init(tend)
    CALL init(new_state)
!$OMP END PARALLEL

    !---------------------------------------------------------------
    ! Vertical diffusion for w-wind
    !---------------------------------------------------------------

    p_nh_metrics => p_nh_state(jg)%metrics

    CALL rbf_vec_interpol_edge( vn, patch, p_int, vt_e, opt_rlend=min_rledge_int-1)

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          inv_rho_ic(jc,jk,jb) = 1._wp / rho_ic(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
END DO
!$OMP END PARALLEL DO

    ! 2) Vertical tendency: evaluated at w point

!$OMP PARALLEL
    CALL init(a)
    CALL init(b)
    CALL init(c)
    CALL init(rhs)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jc,jb,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c

      ! DO jc = i_startidx_c(jb), i_endidx_c(jb)
      !   a  (jc,1,jb) = 0._wp
      !   a  (jc,2,jb) = 0._wp
      !   c  (jc,1,jb) = 0._wp
      !   c  (jc,nlev,jb) = 0._wp
      ! END DO

     !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
     !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 3, nlev-1
      ! DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)

            a(jc,jk,jb)   = - 2._wp * km_c(jc,jk-1,jb) * inv_dzf(jc,jk-1,jb) * &
                            inv_dzh(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

            c(jc,jk,jb)   = - 2._wp * km_c(jc,jk,jb) * inv_dzf(jc,jk,jb)     * &
                            inv_dzh(jc,jk,jb) * inv_rho_ic(jc,jk,jb)

            b(jc,jk,jb)   =  1._wp/dtime - a(jc,jk,jb) - c(jc,jk,jb)

            rhs(jc,jk,jb) =  pwp1(jc,jk,jb) / dtime +                                       &
                          2._wp * ( km_c(jc,jk,jb)   * z_1by3 * div_c(jc,jk,jb)     -            &
                                    km_c(jc,jk-1,jb) * z_1by3 * div_c(jc,jk-1,jb) ) *            &
                          inv_dzh(jc,jk,jb) * inv_rho_ic(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL

      ! TODO: why this? The formulation is the same as in the loop above.

      ! Boundary treatment
      !--------------------------------------------------------
      ! jk = 2 (w == 0)
      !--------------------------------------------------------
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = i_startidx_c(jb), i_endidx_c(jb)
        c(jc,2,jb)   = - 2._wp * km_c(jc,2,jb) * inv_dzf(jc,2,jb) *         &
                      inv_dzh(jc,2,jb) * inv_rho_ic(jc,2,jb)

        b(jc,2,jb)   = 1._wp/dtime - c(jc,2,jb) + 2._wp * km_c(jc,1,jb) *                                &
                    inv_dzf(jc,1,jb) *                                   &
                    inv_dzh(jc,2,jb) * inv_rho_ic(jc,2,jb)

        rhs(jc,2,jb) = pwp1(jc,2,jb) / dtime +                                           &
                    2._wp * ( km_c(jc,2,jb) * z_1by3 * div_c(jc,2,jb) -                       &
                              km_c(jc,1,jb) * z_1by3 * div_c(jc,1,jb) ) *                     &
                    inv_dzh(jc,2,jb) * inv_rho_ic(jc,2,jb)
      END DO
      !$ACC END PARALLEL LOOP
      !--------------------------------------------------------
      ! jk = nlev (w == 0)
      !--------------------------------------------------------
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jc = i_startidx_c(jb), i_endidx_c(jb)
        a(jc,nlev,jb)   = - km_c(jc,nlev-1,jb) * inv_dzf(jc,nlev-1,jb) *    &
                          inv_dzh(jc,nlev,jb) * 2._wp *                   &
                          inv_rho_ic(jc,nlev,jb)

        b(jc,nlev,jb)   =   1._wp/dtime - a(jc,nlev,jb) + 2._wp * km_c(jc,nlev,jb) *                     &
                          inv_dzf(jc,nlev,jb) *                           &
                          inv_dzh(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)

        rhs(jc,nlev,jb) =   pwp1(jc,nlev,jb) / dtime +                                   &
                          2._wp * ( km_c(jc,nlev,jb) * z_1by3 * div_c(jc,nlev,jb) -            &
                                    km_c(jc,nlev-1,jb) * z_1by3 * div_c(jc,nlev-1,jb) ) *      &
                          inv_dzh(jc,nlev,jb) * inv_rho_ic(jc,nlev,jb)
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC WAIT
      CALL tdma_solver_vec(a(:,:,jb),b(:,:,jb),c(:,:,jb),rhs(:,:,jb),2,nlev,i_startidx_c(jb),i_endidx_c(jb),var_new(:,:,jb))

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          tend(jc,jk,jb) = ( var_new(jc,jk,jb) - pwp1(jc,jk,jb) ) / dtime 
        END DO
      END DO
      !$ACC END PARALLEL LOOP

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
    CALL init(hori_tend_e)
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

          dvn2  = pum1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,2)%v1 +  &
                  pvm1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,2)%v2 -  &
                  pum1(jcn,jk,jbn) * patch%edges%primal_normal_cell(je,jb,2)%v1   -  &
                  pvm1(jcn,jk,jbn) * patch%edges%primal_normal_cell(je,jb,2)%v2

          flux_up_c = km_ic(jcn,jk,jbn) * (                                            &
                      dvn2 * inv_dzh(jcn,jk,jbn) +                &
                      ( w_vert(ividx(je,jb,4),jk,ivblk(je,jb,4)) - w_ie(je,jk,jb) ) *  &
                      patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp )

          jcn   = iecidx(je,jb,1)
          jbn   = iecblk(je,jb,1)

          dvn1  = pum1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v1 +  &
                  pvm1(jcn,jk-1,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v2 -  &
                  pum1(jcn,jk,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v1   -  &
                  pvm1(jcn,jk,jbn) * patch%edges%primal_normal_cell(je,jb,1)%v2


          flux_dn_c = km_ic(jcn,jk,jbn) * (                                            &
                      dvn1 * inv_dzh(jcn,jk,jbn) +                &
                      ( w_ie(je,jk,jb) - w_vert(ividx(je,jb,3),jk,ivblk(je,jb,3)) ) *  &
                      patch%edges%inv_vert_vert_length(je,jb) * 2.0_wp )


         ! tendency in tangential direction
         ! flux = visc_v * D_32_v where D_32(= D_23) is calculated at half level
         ! between vertex and edge center

          jvn  = ividx(je,jb,2)
          jbn  = ivblk(je,jb,2)

          dvt2 = ( u_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v1 +            &
                   v_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v2 +            &
                   vt_e(je,jk-1,jb) ) * 0.5_wp  -                                                 &
                 ( u_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v1 +              &
                   v_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,2)%v2 +              &
                   vt_e(je,jk,jb) ) * 0.5_wp

          flux_up_v = km_iv(jvn,jk,jbn) * ( dvt2 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +   &
                      patch%edges%tangent_orientation(je,jb) * ( w_vert(jvn,jk,jbn) -           &
                      w_ie(je,jk,jb) ) / patch%edges%edge_cell_length(je,jb,2) )


          jvn  = ividx(je,jb,1)
          jbn  = ivblk(je,jb,1)

          dvt1 = ( u_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v1 +            &
                   v_vert(jvn,jk-1,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v2 +            &
                   vt_e(je,jk-1,jb) ) * 0.5_wp           - &
                 ( u_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v1 +              &
                   v_vert(jvn,jk,jbn) * patch%edges%dual_normal_vert(je,jb,1)%v2 +              &
                   vt_e(je,jk,jb) ) * 0.5_wp


          flux_dn_v = km_iv(jvn,jk,jbn) * ( dvt1 * p_nh_metrics%inv_ddqz_z_half_v(jvn,jk,jbn) +   &
                      patch%edges%tangent_orientation(je,jb) * ( w_ie(je,jk,jb) -               &
                      w_vert(jvn,jk,jbn) ) / patch%edges%edge_cell_length(je,jb,1) )

          hori_tend_e(je,jk,jb) = ( flux_up_c - flux_dn_c ) *                                        &
                                 patch%edges%inv_dual_edge_length(je,jb) +                      &
                                 patch%edges%tangent_orientation(je,jb)  *                      &
                               ( flux_up_v - flux_dn_v ) *                                        &
                                 patch%edges%inv_primal_edge_length(je,jb) * 2._wp

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
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
      DO jk = 2, nlev
        DO jc = i_startidx_c(jb), i_endidx_c(jb)
          tend(jc,jk,jb) = tend(jc,jk,jb) + inv_rho_ic(jc,jk,jb)             *                     &
                               ( hori_tend_e(ieidx(jc,jb,1),jk,ieblk(jc,jb,1)) *                     &
                                 p_int%e_bln_c_s(jc,1,jb) +                                       &
                                 hori_tend_e(ieidx(jc,jb,2),jk,ieblk(jc,jb,2)) *                     &
                                 p_int%e_bln_c_s(jc,2,jb) +                                       &
                                 hori_tend_e(ieidx(jc,jb,3),jk,ieblk(jc,jb,3)) *                     &
                                 p_int%e_bln_c_s(jc,3,jb) )
          ! new_state(jc,jk,jb) = state(jc,jk,jb) + tend(jc,jk,jb) * dtime
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!$OMP END DO
!$OMP DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk_c,i_endblk_c
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR ASYNC(1)
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

    ! Update diagnostics at end of time step
    CALL this%Update_diagnostics()

  CONTAINS

    SUBROUTINE sync_uvml_s(u, v)
      REAL(wp), TARGET, INTENT(inout) :: u(:,:), v(:,:)
      REAL(wp), POINTER :: pu(:,:,:), pv(:,:,:)
      CALL insert_dimension(pu, u, 2)
      CALL insert_dimension(pv, v, 2)
      CALL sync_patch_array_mult(SYNC_C, patch, 2, pu, pv)
    END SUBROUTINE sync_uvml_s

  END SUBROUTINE Compute

  ! SUBROUTINE Compute_explicit(this)

  !   CLASS(t_vdf), INTENT(inout) :: this

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':Compute_explicit'

  !   CALL this%sfc%Compute()
  !   CALL this%atmo%Compute()

  ! END SUBROUTINE Compute_explicit

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

  SUBROUTINE Update_diagnostics(this)

    USE mo_tmx_surface_interface, ONLY: compute_2m_temperature, compute_10m_wind
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
    
    SELECT TYPE (v => this%atmo%inputs)
    TYPE IS (t_vdf_atmo_inputs)
      ins_atmo => v
    END SELECT
    
    SELECT TYPE (v => this%atmo%diagnostics)
    TYPE IS (t_vdf_atmo_diagnostics)
      diags_atmo => v
    END SELECT
    

    SELECT TYPE (v => this%sfc%config)
    TYPE IS (t_vdf_sfc_config)
      conf_sfc => v
    END SELECT
    
    SELECT TYPE (v => this%sfc%inputs)
    TYPE IS (t_vdf_sfc_inputs)
      ins_sfc => v
    END SELECT
    
    SELECT TYPE (v => this%sfc%diagnostics)
    TYPE IS (t_vdf_sfc_diagnostics)
      diags_sfc => v
    END SELECT
    

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
      
    END DO

    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%t2m_tile,     diags_sfc%t2m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%u10m_tile,    diags_sfc%u10m)
    CALL average_tiles(domain_sfc, ins_sfc%fract_tile, diags_sfc%nvalid, diags_sfc%indices, diags_sfc%v10m_tile,    diags_sfc%v10m)
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

END MODULE mo_vdf 
