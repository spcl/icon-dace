! Computation of vertical tracer flux
!
! Vertical fluxes are calculated at triangle centers on half-levels.
! Possible options for vertical flux calculation include
! - first order Godunov method (UP1)
! - third order PPM method without CFL restriction
! - third order PSM method without CFL restriction
!
! Semi-monotone and monotone limiters are available for PPM
!
! These routines compute only the correct half level value of
! 'c*w'. The vertical divergence is computed in step_advection.
!
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
MODULE mo_advection_vflux

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int,   &
    &                               iup_v, ippm_v, ipsm_v,                      &
    &                               islopel_vsm, islopel_vm, ifluxl_vpd,        &
    &                               ino_flx, izero_grad, iparent_flx
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_math_utilities,      ONLY: tdma_solver_vec
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level, lvert_nest, timers_level, iqtke
  USE mo_advection_config,    ONLY: t_advection_config, advection_config,     &
    &                               lcompute, lcleanup, t_trList 
  USE mo_advection_vlimit,    ONLY: v_limit_parabola_mo, v_limit_parabola_sm,  &
   &                                vflx_limiter_pd,                           &
   &                                v_limit_slope_mo, v_limit_slope_sm,        &
   &                                v_limit_face_mc_mo, v_limit_face_mc_sm
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: global_max
  USE mo_mpi,                 ONLY: process_mpi_stdio_id, my_process_is_stdio, get_my_mpi_work_id, &
                                    get_glob_proc0, comm_lev
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device
  USE mo_timer,               ONLY: timer_adv_vflx, timer_start, timer_stop


  IMPLICIT NONE

  PRIVATE


  PUBLIC :: vert_upwind_flux
  PUBLIC :: upwind_vflux_up
  PUBLIC :: upwind_vflux_ppm
  PUBLIC :: upwind_vflux_ppm4gpu
  PUBLIC :: implicit_sedim_tracer

  CHARACTER(len=*), PARAMETER :: modname = 'mo_advection_vflux'

  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !
  !

  !! Calculation of vertical upwind flux at triangle centers on half levels
  !!
  !! Calculation of vertical upwind flux at triangle centers on half levels
  !! using either
  !! - the first order Godunov method (UP1)
  !! - the third order PPM method
  !! - the third order PSM method
  !!
  !
  ! !LITERATURE
  ! see below
  !
  SUBROUTINE vert_upwind_flux( p_patch, p_cc, p_mflx_contra_v,                &
    &                      p_dtime, p_cellhgt_mc_now, p_cellmass_now,         &
    &                      lprint_cfl, p_upflux, q_ubc, q_int,                &
    &                      opt_rlstart, opt_rlend  )

    CHARACTER(len=*), PARAMETER :: routine = modname//':vert_upwind_flux'

    TYPE(t_patch), INTENT(IN) ::  &  !< patch on which computation is 
      &  p_patch                             !< performed

    REAL(wp), INTENT(IN) ::  &      !< advected cell centered variable
      &  p_cc(:,:,:,:)              !< dim: (nproma,nlev,nblks_c,ntracer)

    REAL(wp), INTENT(INOUT) ::  &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)     !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) :: p_dtime !< time step

    REAL(wp), INTENT(IN) ::  &      !< cell height defined at full levels for
      &  p_cellhgt_mc_now(:,:,:)    !< time step n (either \Delta p or \Delta z)
                                    !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN)::  &       !< NH: density times cell thickness at cell center
      &  p_cellmass_now(:,:,:)      !< at time step n [kg/m**2]
                                    !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &      !< determines if vertical CFL number shall be printed
      &  lprint_cfl                 !< in routine upwind_vflux_ppm

    REAL(wp), INTENT(INOUT) :: &    !< variable in which the upwind flux is stored
      &  p_upflux(:,:,:,:)          !< dim: (nproma,nlevp1,nblks_c,ntracer)

    REAL(wp), INTENT(IN)    :: &    !< tracer mass fraction at (nest) upper boundary 
      &  q_ubc(:,:,:)               !< NH: [kg/kg]

    REAL(wp), INTENT(OUT)   :: &    !< tracer mass fraction at child nest interface level 
      &  q_int(:,:,:)               !< NH: [kg/kg]
                                        !< dim: (nproma,ntracer,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    TYPE(t_trList), POINTER ::  &      !< pointer to tracer sublist
      &  trAdvect

    INTEGER :: jt, nt                  !< tracer index and loop index
    INTEGER :: jc, jb                  !< cell and block loop index
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: i_rlstart_c, i_rlend_c
    INTEGER :: iadv_min_slev           !< scheme specific minimum slev

    TYPE(t_advection_config), POINTER :: advconf

    REAL(wp) :: z_mflx_contra_v !< auxiliary variable for computing vertical nest interface quantities

    !-----------------------------------------------------------------------

    IF (timers_level > 2) CALL timer_start(timer_adv_vflx)

    ! pointer to advection_config(p_patch%id) to save some paperwork
    advconf => advection_config(p_patch%id)

    trAdvect => advconf%trAdvect

    !
    ! Loop over different tracers
    !
    ! Note (DR): Since we define ivadv_tracer as a 1D array of dimension
    ! ntracer, we can select different flux calculation methods for
    ! different tracers. We may also decide not to advect special
    ! tracers. Furthermore, options regarding the desired limiter can be passed 
    ! via the argument list. The same is true for precomputed lcompute/lcleanup
    ! values.
    DO nt = 1, trAdvect%len

      jt = trAdvect%list(nt)

      IF (.NOT. PRESENT(opt_rlend) .OR. (jt == iqtke .AND. advconf%iadv_tke == 1)) THEN
        i_rlend_c = min_rlcell_int
      ELSE
        i_rlend_c = opt_rlend
      ENDIF

      ! Select desired flux calculation method
      SELECT  CASE( advconf%ivadv_tracer(jt) )

      CASE( iup_v )
        ! CALL first order upwind
        !$ACC WAIT
        CALL upwind_vflux_up(                                &
          &         p_patch        = p_patch,                & !in
          &         p_cc           = p_cc(:,:,:,jt),         & !in
          &         p_iubc_adv     = advconf%iubc_adv,       & !in
          &         p_mflx_contra_v= p_mflx_contra_v(:,:,:), & !in 
          &         p_upflux       = p_upflux(:,:,:,jt),     & !out
          &         opt_q_ubc      = q_ubc(:,jt,:),          & !in
          &         opt_slev       = advconf%iadv_slev(jt),  & !in
          &         opt_rlstart    = opt_rlstart,            & !in
          &         opt_rlend      = i_rlend_c               ) !in


      CASE( ippm_v, ipsm_v )

        iadv_min_slev = advconf%ppm_v%iadv_min_slev

        ! CALL third order PPM/PSM (unrestricted timestep-version) (i.e. CFL>1)
        CALL upwind_vflux_ppm(                                           &
          &         p_patch             = p_patch,                       & !in
          &         p_cc                = p_cc(:,:,:,jt),                & !in
          &         p_iubc_adv          = advconf%iubc_adv,              & !in
          &         p_mflx_contra_v     = p_mflx_contra_v(:,:,:),        & !in
          &         p_dtime             = p_dtime,                       & !in
          &         ld_compute          = lcompute%ppm_v(jt),            & !in
          &         ld_cleanup          = lcleanup%ppm_v(jt),            & !in
          &         p_itype_vlimit      = advconf%itype_vlimit(jt),      & !in
          &         p_ivlimit_selective = advconf%ivlimit_selective(jt), & !in
          &         p_cellhgt_mc_now    = p_cellhgt_mc_now(:,:,:),       & !in
          &         p_cellmass_now      = p_cellmass_now(:,:,:),         & !in
          &         lprint_cfl          = lprint_cfl,                    & !in
          &         ivadv_tracer        = advconf%ivadv_tracer(jt),      & !in
          &         p_upflux            = p_upflux(:,:,:,jt),            & !out
          &         opt_q_ubc           = q_ubc(:,jt,:),                 & !in
          &         opt_slev            = advconf%iadv_slev(jt),         & !in
          &         opt_ti_slev         = iadv_min_slev,                 & !in
          &         opt_rlstart         = opt_rlstart,                   & !in
          &         opt_rlend           = i_rlend_c                      ) !in

      END SELECT

    END DO  ! Tracer loop


    !
    ! get face value "q_int" at vertical boundary of child nest
    !
    ! determine if upper boundary values are needed
    IF (lvert_nest .AND. (p_patch%nshift_child > 0)) THEN 

      ! refinement control start/end level for cells
      i_rlstart_c = 1
      i_rlend_c   = min_rlcell_int

      i_startblk = p_patch%cells%start_block(i_rlstart_c)
      i_endblk   = p_patch%cells%end_block(i_rlend_c)

      !$ACC DATA PRESENT(p_mflx_contra_v, p_patch, q_int, p_upflux, trAdvect)

!$OMP PARALLEL DO PRIVATE(jb,jt,jc,nt,i_startidx,i_endidx,z_mflx_contra_v) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
          &                 i_startidx, i_endidx, i_rlstart_c, i_rlend_c )

        ! Be sure to avoid division by zero
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP SEQ
        DO nt = 1, trAdvect%len
          !$ACC LOOP GANG VECTOR PRIVATE(z_mflx_contra_v)
          DO jc = i_startidx, i_endidx
            jt = trAdvect%list(nt)
            z_mflx_contra_v = SIGN( MAX(ABS(p_mflx_contra_v(jc,p_patch%nshift_child,jb)),dbl_eps), &
            &                               p_mflx_contra_v(jc,p_patch%nshift_child,jb) )
            q_int(jc,jt,jb) = p_upflux(jc,p_patch%nshift_child,jb,jt) / z_mflx_contra_v
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!$OMP END PARALLEL DO

    !$ACC END DATA

    ENDIF

    IF (timers_level > 2) CALL timer_stop(timer_adv_vflx)

  END SUBROUTINE vert_upwind_flux




  !-------------------------------------------------------------------------
  !! The first order Godunov method
  !!
  !! Calculation of time averaged vertical tracer fluxes using the first
  !! order Godunov method.
  !!
  SUBROUTINE upwind_vflux_up( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v,    &
    &                         p_upflux, opt_q_ubc, opt_slev,                 &
    &                         opt_rlstart, opt_rlend )

!!$    CHARACTER(len=*), PARAMETER :: routine = modname//':upwind_vflux_up'

    TYPE(t_patch), INTENT(IN) ::  & !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::   &   !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(IN) ::   &   !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(INOUT) ::  & !< vertical tracer flux at half levels
      &  p_upflux(:,:,:)          !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< tracer mass fraction at (nest) upper boundary 
      &  opt_q_ubc(:,:)                 !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    REAL(wp) ::  &
      &  zq_ubc(nproma,p_patch%nblks_c)

    INTEGER  :: slev                   !< vertical start level
    INTEGER  :: nlev, nlevp1           !< number of full and half levels
    INTEGER  :: jc, jk, jb             !< index of cell, vertical level and block
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend
    !-------------------------------------------------------------------------


    !$ACC DATA CREATE(zq_ubc) PRESENT(p_cc, p_mflx_contra_v, p_upflux) &
    !$ACC   IF(i_am_accel_node)

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_q_ubc) ) THEN
      !$ACC KERNELS PRESENT(opt_q_ubc, zq_ubc) ASYNC(1) IF(i_am_accel_node)
      zq_ubc(:,:) = opt_q_ubc(:,:)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS PRESENT(zq_ubc) ASYNC(1) IF(i_am_accel_node)
      zq_ubc(:,:) = 0._wp
      !$ACC END KERNELS
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    !
    ! advection is done with 1st order upwind scheme,
    ! i.e. a piecewise constant approx. of the cell centered values
    ! is used.
    !
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev+1, nlev
        DO jc = i_startidx, i_endidx
          ! calculate vertical tracer flux   -- removed flaky laxfr macro
          p_upflux(jc,jk,jb) = p_mflx_contra_v(jc,jk,jb) *                    &
                               MERGE( p_cc(jc,jk,jb),p_cc(jc,jk-1,jb),        &
                                      p_mflx_contra_v(jc,jk,jb) .GE. 0.0_wp ) 
        END DO ! end loop over cells
      ENDDO ! end loop over vertical levels
      !$ACC END PARALLEL

      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(i_start      = i_startidx,                 & !in
        &              i_end        = i_endidx,                   & !in
        &              iubc_adv     = p_iubc_adv,                 & !in
        &              llbc_no_flux = .TRUE.,                     & !in
        &              mflx_top     = p_mflx_contra_v(:,slev,jb), & !in
        &              q_top        = zq_ubc(:,jb),               & !in
        &              upflx_top    = p_upflux(:,slev,jb),        & !out
        &              upflx_bottom = p_upflux(:,nlevp1,jb))        !out
      !$ACC WAIT

    ENDDO ! end loop over blocks

    !$ACC WAIT(1)
    !$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE upwind_vflux_up



  !-------------------------------------------------------------------------
  !! The third order PPM/PSM scheme for large time steps (CFL>1)
  !!
  !! Calculation of time averaged vertical tracer fluxes or tracer edge 
  !! values using the third order PPM/PSM scheme. This scheme can handle 
  !! large time steps (i.e. CFL>1)
  !!
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201 (PPM)
  ! - Carpenter et al. (1989), MWR, 118, 586-612  (PPM)
  ! - Zerroukat et al. (2006), Int. J. Numer. Meth. Fluids, 51, 1297-1318 (PSM)
  ! - Lin et al (1994), MWR, 122, 1575-1593 (filtered reconstruction)
  ! - Lin and Rood (1996), MWR, 124, 2046-2070 (CFL-independent version)
  !
  SUBROUTINE upwind_vflux_ppm( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v,     &
    &                      p_dtime,  ld_compute, ld_cleanup, p_itype_vlimit,   &
    &                      p_ivlimit_selective,                                &
    &                      p_cellhgt_mc_now, p_cellmass_now,                   &
    &                      lprint_cfl, ivadv_tracer,                           &
    &                      p_upflux, opt_lout_edge, opt_q_ubc, opt_slev,       &
    &                      opt_ti_slev, opt_rlstart, opt_rlend, opt_elev )

    CHARACTER(len=*), PARAMETER :: routine = modname//':upwind_vflux_ppm'

    TYPE(t_patch), INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::  &    !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(INOUT) ::  & !< contravariant vertical mass flux
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< time step
      &  p_dtime

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. compute geometrical terms
      &  ld_compute

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. clean up geometrical terms
      &  ld_cleanup

    INTEGER, INTENT(IN)  ::  &    !< parameter to select the limiter for
      &  p_itype_vlimit           !< vertical transport

    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< NH: density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)    !< at time step n [kg/m**2]
                                  !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &    !< determines if vertical CFL number shall be written out
      &  lprint_cfl

    INTEGER, INTENT(IN) ::   &    !< type of vertical transport (PPM or PSM)
      &  ivadv_tracer

    REAL(wp), INTENT(INOUT) :: &  !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    LOGICAL, INTENT(IN), OPTIONAL ::  & !< optional: output edge value (.TRUE.),
      &  opt_lout_edge                  !< or the flux across the edge 
                                        !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< tracer mass fraction at (nest) upper boundary 
      &  opt_q_ubc(:,:)                 !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end   level (for sedimentation)
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart                   !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    LOGICAL  :: l_out_edgeval     !< corresponding local variable; default 
                                  !< .FALSE. i.e. output flux across the edge

    REAL(wp) :: &                 !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< integer fluxes
      &  z_iflx(nproma,p_patch%nlevp1)

    REAL(wp) :: &                 !< difference between upper and lower face value times 0.5
      &  z_delta_q(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< 1/6 * a6,i (see Colella and Woodward (1984))
      &  z_a1(nproma,p_patch%nlev)

    REAL(wp) :: &                 !< p_cc of upwind cell
      &  q_up

    REAL(wp) :: &                 !< z_delta_q of upwind cell
      &  dq_up

    REAL(wp) :: &                 !< z_a1 of upwind cell, including sign flip  
      &  a_up                     !< for w < 0

    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: ikm1, ikp1               !< vertical level minus and plus one
    INTEGER  :: slev, slevp1             !< vertical start level and start level +1
    INTEGER  :: slev_ti, slevp1_ti       !< vertical start level (+1)  (tracer independent part)
    INTEGER  :: nlev, nlevp1             !< number of full and half levels

    ! JF: for treatment of sedimentation
    INTEGER  :: elev, elev_lim           !< vertical end level
    LOGICAL  :: llbc_no_flux             !< TRUE: apply 'no flux' lower boundary condition
    INTEGER  :: ik                       !< = MIN(jk,nlev)

    INTEGER  :: ji                       !< loop variable for index list
    INTEGER  :: ist                      !< status variable
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend

    INTEGER  :: nlist_max                !< maximum number of index lists
    INTEGER  :: nlist                    !< list loop variable

    REAL(wp) :: wsign                    !< wind direction: introduced, in order to merge flux formula  
                                         !< for w>0 and w<0.
                                         !< +1, if w >0
                                         !< -1, if w <0

    REAL(wp), ALLOCATABLE, SAVE :: &     !< fractional (mass weighted) Courant number 
      &  z_cflfrac(:,:,:)                !< always positive

    REAL(wp), ALLOCATABLE, SAVE ::    &  !< maximum vertical Courant number
      &  max_cfl_blk(:)                  !< per block

    INTEGER, ALLOCATABLE, SAVE  ::    &  !< Index lists, level lists and list dimensions 
      &  i_indlist(:,:,:),    &          !< for points with CFL>1,2,3
      &  i_levlist(:,:,:),    &
      &  i_listdim(:,:)

    INTEGER, ALLOCATABLE, SAVE  ::    &  !< upwind shifted vertical index
      &  jk_shifted(:,:,:)               !< jk +/- s, with positive shift s

    INTEGER  :: jks                      !< shifted vertical index
                                         !< can vary betwen jk and jk_shifted

    INTEGER  :: bot_bound                !< shifted index jk_shifted must fall within the 
                                         !< range [top_bound, bot_bound]. Note that the 
                                         !< permissible range depends on the sign of w.
                                         !< Note that the variable top_bound is currently 
                                         !< not needed. It can savely be replaced by 
                                         !< slevp1_ti (see below) 

    INTEGER  :: counter, counter_ji      !< check whether any of the points has 
                                         !< CFL>nlist

    INTEGER  :: jg                       !< patch ID

    REAL(wp) ::   &                      !< high order flux
      &  z_flx_frac_high

    REAL(wp) ::  &
      &  zq_ubc(nproma,p_patch%nblks_c)

    REAL(wp) ::   &                      !< maximum CFL within one layer, and domain-wide maximum
      &  max_cfl_lay(p_patch%nlevp1,p_patch%nblks_c), max_cfl_tot, max_cfl_lay_tot(p_patch%nlevp1)

    REAL(wp) ::   &                      !< auxiliary for fractional CFL number computation
      &  z_aux(nproma)

    REAL(wp) :: rdtime                   !< 1/dt







    !-----------------------------------------------------------------------

    ! inverse of time step for computational efficiency
    rdtime = 1._wp/p_dtime

    ! get patch ID
    jg = p_patch%id

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
      slevp1= opt_slev + 1
    ELSE
      slev  = 1
      slevp1= 2
    END IF

    ! check optional arguments
    IF ( PRESENT(opt_ti_slev) ) THEN
      slev_ti  = opt_ti_slev
      slevp1_ti= opt_ti_slev + 1
    ELSE
      slev_ti  = 1
      slevp1_ti= 2
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_q_ubc) ) THEN
      zq_ubc(:,:) = opt_q_ubc(:,:)
    ELSE
      zq_ubc(:,:) = 0._wp
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! maximum number of lists
    nlist_max = advection_config(jg)%ivcfl_max

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! check optional arguments
    llbc_no_flux = .TRUE.
    IF ( PRESENT(opt_elev) ) THEN
      IF ( opt_elev == nlevp1 ) THEN
        elev = nlevp1
        llbc_no_flux = .FALSE.  ! e.g. when used for sedimentation
      ELSE
        elev = nlev
      END IF
    ELSE
      elev = nlev
    END IF
    elev_lim = nlev

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    !
    ! advection is done with an upwind scheme where a piecwise parabolic
    ! approx. of the subgrid distribution is used.
    !
    IF ( ld_compute ) THEN
      ! allocate temporary arrays 
      ALLOCATE( i_indlist(nproma*nlevp1,nlist_max,p_patch%nblks_c),    &
        &       i_levlist(nproma*nlevp1,nlist_max,p_patch%nblks_c),    &
        &       i_listdim(nlist_max,p_patch%nblks_c),                  &
        &       jk_shifted(nproma,nlevp1,p_patch%nblks_c),             &
        &       z_cflfrac(nproma,nlevp1,p_patch%nblks_c),              &
        &       max_cfl_blk(p_patch%nblks_c), STAT=ist                 )
      IF (ist /= SUCCESS) THEN
        CALL finish(routine,                                           &
          &  'allocation for i_indlist, i_levlist, i_listdim, '    //  &
          &  'jk_shifted, z_cflfrac, max_cfl_blk  failed '    )
      ENDIF
    END IF

!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jc,ik,ikm1,i_startidx,i_endidx,nlist,                  &
!$OMP            counter,counter_ji,z_aux,ikp1,ji,wsign,jks,                  &
!$OMP            z_iflx,z_delta_q,z_a1,z_face,z_face_up,z_face_low,           &
!$OMP            z_flx_frac_high,q_up,dq_up,a_up,bot_bound) ICON_OMP_GUIDED_SCHEDULE
  DO jb = i_startblk, i_endblk

    CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
      &                 i_startidx, i_endidx, i_rlstart, i_rlend )


    ! The contravariant mass flux should never exactly vanish
    !
    IF (l_out_edgeval) THEN
      DO jk = slevp1, elev
        p_mflx_contra_v(i_startidx:i_endidx,jk,jb) =                            &
        &              p_mflx_contra_v(i_startidx:i_endidx,jk,jb)               &
        &              + SIGN(dbl_eps,p_mflx_contra_v(i_startidx:i_endidx,jk,jb))
      ENDDO
    ENDIF

    !
    ! 1. Compute density weighted (fractional) Courant number 
    !    for w<0 and w>0 and integer shift s
    !
    IF (ld_compute) THEN

      ! set start values for index list dimension
      i_listdim(1:nlist_max,jb) = 0

      ! initialize (fractional) Courant number at top
      z_cflfrac(i_startidx:i_endidx,slev_ti,jb) = 0._wp
      ! initialize (fractional) Courant number at bottom
      z_cflfrac(i_startidx:i_endidx,nlevp1,jb) = 0._wp

      !
      ! compute (fractional) Courant number
      !
      DO jk = slevp1_ti, elev

        ik   = MIN(jk,nlev)
        ikm1 = jk-1

        DO jc = i_startidx, i_endidx

          ! initialize shifted index with upwind index
          jk_shifted(jc,jk,jb) = MERGE(ik, ikm1, p_mflx_contra_v(jc,jk,jb) > 0._wp)

          z_aux(jc) = p_dtime * ABS(p_mflx_contra_v(jc,jk,jb))

          ! compute CFL number
          z_cflfrac(jc,jk,jb) = z_aux(jc) / MERGE(p_cellmass_now(jc,ik,jb),p_cellmass_now(jc,ikm1,jb), &
                                                  p_mflx_contra_v(jc,jk,jb) > 0._wp)

        ENDDO
        max_cfl_lay(jk,jb) = MAXVAL(z_cflfrac(i_startidx:i_endidx,jk,jb))


        ! If CFL>1 then split the CFL number into the fractional CFL number 
        ! and the index shift s.
        IF (max_cfl_lay(jk,jb) <= 1._wp) CYCLE


        ! initialize list number
        nlist = 0

        ! checks whether there exists any point with 'large CFL number'
        counter = 1

        ! loop until no point has CFL > 1, or nlist > nlist_max
        DO WHILE(counter > 0 .AND. nlist < nlist_max )

          ! get number of current list
          nlist     = nlist + 1
          ! re-initialize counter for CFL>nlist
          counter   = 0

          ! copy value from counter in vector form to counter in scalar form, 
          ! since otherwise the following DO-Loop will not vectorize.
          counter_ji = i_listdim(nlist,jb)

          DO jc = i_startidx, i_endidx

            ! jk_shifted must fall within the range [top_bound, bot_bound] in order 
            ! to pass the following if condition. Unfortunately, the range depends on 
            ! the sign of w. 
            ! Note that for the time being we can skip the computation of top_bound, 
            ! as the CFL computation only starts at slevp1_ti anyways.  
            bot_bound = MERGE(nlev-1 , nlev     , p_mflx_contra_v(jc,jk,jb) > 0._wp)
            !top_bound = MERGE(slev_ti, slevp1_ti, p_mflx_contra_v(jc,jk,jb) > 0._wp)

            IF ( z_aux(jc) > p_cellmass_now(jc,jk_shifted(jc,jk,jb),jb)  &
              &  .AND. jk_shifted(jc,jk,jb) <= bot_bound                 &
              &  .AND. jk_shifted(jc,jk,jb) >= slevp1_ti                 ) THEN


              z_aux(jc) = z_aux(jc) - p_cellmass_now(jc,jk_shifted(jc,jk,jb),jb)

              ! Index shift
              jk_shifted(jc,jk,jb) = jk_shifted(jc,jk,jb)  &
                &                    + MERGE(1, -1, p_mflx_contra_v(jc,jk,jb) > 0._wp )

              ! tests whether we need to loop once again
              counter = counter + 1

              ! Fill index lists with those points that need index shifts
              ! Note that we have to use a scalar counter instead of a vector, like
              ! i_listdim(nlist,jb). Otherwise this loop will not vectorize.  
              counter_ji = counter_ji + 1
              i_indlist(counter_ji,nlist,jb) = jc
              i_levlist(counter_ji,nlist,jb) = jk

              ! compute fractional Courant number
              z_cflfrac(jc,jk,jb) = z_aux(jc) / p_cellmass_now(jc,jk_shifted(jc,jk,jb),jb)
            ENDIF

          END DO ! end loop over cells

          ! store index of current index list
          ! after the last jk-loop this will be the list dimension
          i_listdim(nlist,jb) = counter_ji

        ENDDO  ! DO WHILE loop
 
      ENDDO ! end loop over vertical levels

      max_cfl_blk(jb) = MAXVAL(max_cfl_lay(slevp1_ti:elev,jb))

    END IF ! ld_compute



    !
    ! 2. Edge value reconstruction
    !
    SELECT CASE(ivadv_tracer)
    CASE (IPPM_V)

      !
      ! PPM Reconstruction following Colella and Woodward (1984)
      !
      CALL compute_face_values_ppm( i_startidx       = i_startidx,               & !in
        &                           i_endidx         = i_endidx,                 & !in
        &                           slev             = slev,                     & !in
        &                           elev             = nlev,                     & !in
        &                           p_itype_vlimit   = p_itype_vlimit,           & !in
        &                           p_cc             = p_cc(:,:,jb),             & !in
        &                           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,jb), & !in
        &                           p_face           = z_face(:,:)               ) !inout


    CASE (IPSM_V)


      !
      ! PSM Reconstruction following Zerroukat et al. (2006)
      !
      CALL compute_face_values_psm( i_startidx       = i_startidx,               & !in
        &                           i_endidx         = i_endidx,                 & !in
        &                           slev             = slev,                     & !in
        &                           elev             = nlev,                     & !in
        &                           p_itype_vlimit   = p_itype_vlimit,           & !in
        &                           p_cc             = p_cc(:,:,jb),             & !in
        &                           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,jb), & !in
        &                           p_face           = z_face(:,:)               ) !inout


    END SELECT  ! ivadv_tracer


      !
      ! 4. Limitation/filtering of first guess parabola (which is based on z_face)
      ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
      ! the limitation procedure.
      ! Therefore 2 additional fields z_face_up and z_face_low are introduced.
      !
      SELECT CASE (p_itype_vlimit)
      CASE(ISLOPEL_VSM)

        ! semi-monotonic (sm) filter
        CALL v_limit_parabola_sm( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in

      CASE(ISLOPEL_VM)

        ! monotonic (mo) filter
        CALL v_limit_parabola_mo( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in


      CASE default

        ! simply copy face values to 'face_up' and 'face_low' arrays
        !
        DO jk = slev, nlev
          ! index of bottom half level
          ikp1 = jk + 1
          z_face_up(i_startidx:i_endidx,jk)  = z_face(i_startidx:i_endidx,jk)
          z_face_low(i_startidx:i_endidx,jk) = z_face(i_startidx:i_endidx,ikp1)
        ENDDO

      END SELECT  ! p_itype_vlimit



      !
      ! 5. Computation of upwind fluxes. IF CFL > 1, the fluxes are the sum of
      !    integer-fluxes and a fractional flux. IF CFL <1 the fluxes are only
      !    comprised of the fractional flux. The fractional flux is calculated
      !    assuming a piecewise parabolic subgrid distribution.
      !

      ! 5a. Compute coefficients of reconstructed parabola as they are used at 
      !     various places below.
      !     Terminology follows Colella (1984)
      !     z_delta_q = 0.5*\Delta q
      !     z_a1 = 1/6*a_6
      !
      DO jk = slev, nlev
        DO jc = i_startidx, i_endidx
          z_delta_q(jc,jk) = 0.5_wp * (z_face_up(jc,jk) - z_face_low(jc,jk))
          z_a1(jc,jk)      = p_cc(jc,jk,jb) - 0.5_wp*(z_face_up(jc,jk) + z_face_low(jc,jk))
        ENDDO
      ENDDO


      !
      ! 5b. First compute fluxes for the CFL<1 case for all grid points
      ! On the grid points where CFL>1, they will be overwritten afterwards 
      ! This part has been adopted from the restricted time step PPM-scheme.
      !
      DO jk = slevp1, elev

        ik   = MIN(jk,nlev)
        ! index of top half level
        ikm1 = jk -1

        DO jc = i_startidx, i_endidx
          q_up  = MERGE(p_cc(jc,ik,jb)  , p_cc(jc,ikm1,jb)         , p_mflx_contra_v(jc,jk,jb) >= 0._wp)
          dq_up = MERGE(z_delta_q(jc,ik), -1._wp*z_delta_q(jc,ikm1), p_mflx_contra_v(jc,jk,jb) >= 0._wp)
          a_up  = MERGE(z_a1(jc,ik)     , z_a1(jc,ikm1)            , p_mflx_contra_v(jc,jk,jb) >= 0._wp)

          !
          ! full flux
          ! fluxes for CFL>1 are corrected lateron
          !
          p_upflux(jc,jk,jb) = p_mflx_contra_v(jc,jk,jb)                      &
            &                * ( q_up                                         &
            &                + (dq_up * (1._wp - z_cflfrac(jc,jk,jb)))        &
            &                - a_up * (1._wp - 3._wp*z_cflfrac(jc,jk,jb)      &
            &                + 2._wp*z_cflfrac(jc,jk,jb)*z_cflfrac(jc,jk,jb)) )


        END DO ! end loop over cells

      ENDDO ! end loop over vertical levels



      !
      ! 5c. Now execute the special computations needed for CFL>1:
      !     Computation of integer fluxes and a fractional flux
      !
      IF (max_cfl_blk(jb) > 1._wp) THEN

        z_iflx(i_startidx:i_endidx,slev:nlevp1) = 0._wp

        ! Loop over all lists (nlist will serve as index shift)
        ! i.e. List nlist=1 contains all points with an index shift of at least 1
        !      List nlist=2 contains all points with an index shift of at least 2 
        !      and so on
        DO nlist = 1, nlist_max

          IF (i_listdim(nlist,jb) == 0) CYCLE

          !
          ! loop over all cells in i_indlist
          !
          ! integer fluxes
          !
!$NEC ivdep
          DO ji=1,i_listdim(nlist,jb)

            ! get jc and jk index from precomputed list
            jc = i_indlist(ji,nlist,jb)
            jk = i_levlist(ji,nlist,jb)

            ! shifted index (depends on the sign of w)
            jks = MERGE(jk+nlist-1, jk-nlist, p_mflx_contra_v(jc,jk,jb) >= 0._wp)

            ! cycle if the model level is in a region where advection is 
            ! turned off for the present variable
            IF (jk < slevp1) CYCLE

            ! cycle if the source model level is in a region where advection is
            ! turned off for the present variable
            IF (jks < slevp1) CYCLE

            ! Integer flux (division by dtime is done at the end)
            z_iflx(jc,jk) = z_iflx(jc,jk) + p_cc(jc,jks,jb) &
                             * p_cellmass_now(jc,jks,jb)

          ENDDO  ! loop over cells in i_indlist_p

        ENDDO ! list loop


        ! Now use the first list (which contains all points with CFL>1)
        ! to compute the corrected full fluxes
        IF (i_listdim(1,jb) > 0) THEN
!$NEC ivdep
          DO ji=1,i_listdim(1,jb)

            ! get jc and jk index from precomputed list
            jc = i_indlist(ji,1,jb)
            jk = i_levlist(ji,1,jb)

            jks = jk_shifted(jc,jk,jb)

            ! cycle if the model level is in a region where advection is 
            ! turned off for the present variable
            IF (jk < slevp1) CYCLE

            ! this is needed in addition in order to avoid accessing non-existing (uninitalized)
            ! source levels for tracers that are not advected on all model levels
            IF (jks < slev) CYCLE

            wsign = MERGE(1._wp, -1._wp, p_mflx_contra_v(jc,jk,jb) >= 0._wp)

            ! fractional high order flux   
            z_flx_frac_high = p_cellmass_now(jc,jks,jb) * z_cflfrac(jc,jk,jb)       &
              &         * ( p_cc(jc,jks,jb)                                         &
              &         + wsign*(z_delta_q(jc,jks) * (1._wp - z_cflfrac(jc,jk,jb))) &
              &         - z_a1(jc,jks)*(1._wp - 3._wp*z_cflfrac(jc,jk,jb)           &
              &         + 2._wp*z_cflfrac(jc,jk,jb)*z_cflfrac(jc,jk,jb)) )

            ! full flux (integer- plus high order fractional flux)
            p_upflux(jc,jk,jb) = wsign*rdtime * (z_iflx(jc,jk) + z_flx_frac_high) 
          ENDDO

        ENDIF

      ENDIF  ! IF (max_cfl_blk(jb) > 1)


      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(i_start      = i_startidx,                 & !in
        &              i_end        = i_endidx,                   & !in
        &              iubc_adv     = p_iubc_adv,                 & !in
        &              llbc_no_flux = llbc_no_flux,               & !in
        &              mflx_top     = p_mflx_contra_v(:,slev,jb), & !in
        &              q_top        = zq_ubc(:,jb),               & !in
        &              upflx_top    = p_upflux(:,slev,jb),        & !out
        &              upflx_bottom = p_upflux(:,nlevp1,jb))        !out


      ! If desired, get edge value of advected quantity 
      IF ( l_out_edgeval ) THEN

        DO jk = slevp1, nlev
          DO jc = i_startidx, i_endidx
            p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb)/p_mflx_contra_v(jc,jk,jb)
          ENDDO
        ENDDO

      ENDIF

      !
      ! 6. If desired, apply positive-definite flux limiter to limit 
      !    computed fluxes (based on work by Zalesak (1979)).
      !
      IF (p_itype_vlimit == IFLUXL_VPD) THEN
        ! positive-definite (pd) flux limiter
        CALL vflx_limiter_pd( p_dtime         = p_dtime,                & !in
          &                   p_cc            = p_cc(:,:,jb),           & !in
          &                   p_rhodz_now     = p_cellmass_now(:,:,jb), & !in
          &                   p_mflx_tracer_v = p_upflux(:,:,jb),       & !inout
          &                   i_startidx      = i_startidx,             & !in
          &                   i_endidx        = i_endidx,               & !in
          &                   slev            = slev,                   & !in
          &                   elev            = elev_lim                ) !in
        !$ACC WAIT
      ENDIF

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL



    !
    ! If desired, print maximum vertical CFL number
    !
    IF ( ld_compute .AND. msg_level >= 10 .AND. lprint_cfl) THEN

      max_cfl_tot = MAXVAL(max_cfl_blk(i_startblk:i_endblk))

      ! Take maximum over all PEs
      IF (msg_level >= 13) THEN
        max_cfl_tot = global_max(max_cfl_tot)
      ELSE
        max_cfl_tot = global_max(max_cfl_tot, iroot=process_mpi_stdio_id)
      ENDIF
      IF (my_process_is_stdio() .OR. comm_lev>0 .AND. get_my_mpi_work_id() == get_glob_proc0() ) THEN
        ! otherwise it is possible that max_cfl_tot is undefined
        WRITE(message_text,'(a,e16.8)') 'maximum vertical CFL =',max_cfl_tot
        CALL message(routine, message_text)
      ENDIF

      ! Add layer-wise diagnostic if the maximum CFL value is close to the stability limit
      IF (msg_level >= 13) THEN
        IF (max_cfl_tot > (nlist_max-1)) THEN
          DO jk = slevp1_ti, nlev
            max_cfl_lay_tot(jk) = MAXVAL(max_cfl_lay(jk,i_startblk:i_endblk))
          ENDDO

          max_cfl_lay_tot(slevp1_ti:nlev) = global_max(max_cfl_lay_tot(slevp1_ti:nlev), iroot=process_mpi_stdio_id)
          DO jk = slevp1_ti,nlev
            WRITE(message_text,'(a,i4,a,e16.8)') 'maximum vertical CFL in layer', jk,' =', max_cfl_lay_tot(jk)
            CALL message(routine, message_text)
          ENDDO
        ENDIF
      ENDIF

    END IF

    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays
      DEALLOCATE( i_indlist, i_levlist, i_listdim,             &
        &         jk_shifted, z_cflfrac, max_cfl_blk, STAT=ist )

      IF (ist /= SUCCESS) THEN
        CALL finish(routine,                                            &
          &  'deallocation for i_indlist, i_levlist, i_listdim_p, ' //  &
          &  'jk_shifted, z_cflfrac, max_cfl_blk failed '      )
      ENDIF
    END IF

  END SUBROUTINE upwind_vflux_ppm





  !-------------------------------------------------------------------------
  !! The third order PPM/PSM scheme for large time steps (CFL>1)
  !! GPU-enabled version without index lists.
  !!
  !! Calculation of time averaged vertical tracer fluxes or tracer edge 
  !! values using the third order PPM/PSM scheme. This scheme can handle 
  !! large time steps (CFL>1).
  !
  ! !LITERATURE
  ! - Colella and Woodward (1984), JCP, 54, 174-201 (PPM)
  ! - Carpenter et al. (1989), MWR, 118, 586-612  (PPM)
  ! - Zerroukat et al. (2006), Int. J. Numer. Meth. Fluids, 51, 1297-1318 (PSM)
  ! - Lin et al (1994), MWR, 122, 1575-1593 (filtered reconstruction)
  ! - Lin and Rood (1996), MWR, 124, 2046-2070 (CFL-independent version)
  !
  SUBROUTINE upwind_vflux_ppm4gpu( p_patch, p_cc, p_iubc_adv, p_mflx_contra_v, &
    &                      p_dtime,  ld_compute, ld_cleanup, p_itype_vlimit,   &
    &                      p_ivlimit_selective,                                &
    &                      p_cellhgt_mc_now, p_cellmass_now, lprint_cfl,       &
    &                      ivadv_tracer,                                       &
    &                      p_upflux, opt_lout_edge, opt_q_ubc, opt_slev,       &
    &                      opt_ti_slev, opt_rlstart, opt_rlend, opt_elev )

    CHARACTER(len=*), PARAMETER :: routine = modname//':upwind_vflux_ppm4gpu'

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is performed
      &  p_patch

    REAL(wp), INTENT(IN) ::  &    !< advected cell centered variable
      &  p_cc(:,:,:)              !< dim: (nproma,nlev,nblks_c)

    INTEGER, INTENT(IN)  ::   &   !< selects upper boundary condition
      &  p_iubc_adv

    REAL(wp), INTENT(INOUT) ::  & !< contravariant vertical mass flux [kg/m**2/s]
      &  p_mflx_contra_v(:,:,:)   !< dim: (nproma,nlevp1,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< time step [s]
      &  p_dtime

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. compute geometric terms
      &  ld_compute

    LOGICAL, INTENT(IN)  ::  &    !< flag, if .TRUE. clean up geometric terms
      &  ld_cleanup

    INTEGER, INTENT(IN)  ::  &    !< parameter to select the limiter for
      &  p_itype_vlimit           !< vertical transport

    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &    !< layer thickness at cell center at time n [m]
      &  p_cellhgt_mc_now(:,:,:)  !< dim: (nproma,nlev,nblks_c)

    REAL(wp), INTENT(IN) ::  &    !< density weighted cell height at full levels
      &  p_cellmass_now(:,:,:)    !< at time step n [kg/m**2]
                                  !< dim: (nproma,nlev,nblks_c)

    LOGICAL, INTENT(IN) ::   &    !< determines if vertical CFL number shall be written out
      &  lprint_cfl

    INTEGER, INTENT(IN) ::   &    !< type of vertical transport (PPM or PSM)
      &  ivadv_tracer

    REAL(wp), INTENT(INOUT) :: &  !< output field, containing the tracer mass flux
      &  p_upflux(:,:,:)          !< or the reconstructed edge value
                                  !< dim: (nproma,nlevp1,nblks_c)

    LOGICAL, INTENT(IN), OPTIONAL ::  & !< optional: output edge value (.TRUE.),
      &  opt_lout_edge                  !< or the flux across the edge 
                                        !< (.FALSE./not specified)

    REAL(wp), INTENT(IN), OPTIONAL :: & !< tracer mass fraction at (nest) upper boundary 
      &  opt_q_ubc(:,:)                 !< dim: (nproma,nblks_c)

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical start level (tracer independent part)
      &  opt_ti_slev

    INTEGER, INTENT(IN), OPTIONAL ::  & !< optional vertical end level (for sedimentation)
      &  opt_elev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
     &  opt_rlstart                    !< only valid for calculation of 'cell value'

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
     &  opt_rlend                      !< (to avoid calculation of halo points)



    ! local vars

    REAL(wp), ALLOCATABLE, SAVE  ::   &  !< sum of integer and fractional Courant number 
      &  z_cfl(:,:,:)                    !< Sign equals sign of w.

    REAL(wp) :: &                        !< absolute value of fractional Courant number
      &  z_cflfrac                       !< i.e. always >=0

    REAL(wp) :: &                        !< face values of transported field
      &  z_face(nproma,p_patch%nlevp1)

    REAL(wp) :: &                        !< face value (upper face)
      &  z_face_up(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< face value (lower face)
      &  z_face_low(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< difference between upper and lower face value times 0.5
      &  z_delta_q(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< 1/6 * a6,i (see Colella and Woodward (1984))
      &  z_a1(nproma,p_patch%nlev)

    REAL(wp) :: &                        !< mass crossing cell face during \Delta t [kg/m**2]
      &  z_mass                          !< can be positive, or negative, depending on the sign of w

    REAL(wp) :: &                        !< integrated value of q (from 0 to z_cflfrac)
      &  z_q_int

    REAL(wp) :: &                        !< integer flux for w>0 or w<0  [kg/m**2/s]
      &  z_iflx

    INTEGER  :: jc, jk, jb               !< index of cell, vertical level and block
    INTEGER  :: ikp1                     !< vertical level plus one
    INTEGER  :: slev, slevp1             !< vertical start level and start level +1
    INTEGER  :: slev_ti, slevp1_ti       !< vertical start level (+1)  (tracer independent part)
    INTEGER  :: nlev, nlevp1             !< number of full and half levels

    INTEGER  :: jk_shift, jks            !< shifted vertical index

    INTEGER  :: js                       !< the shift itself (always positive), i.e. jks = jk \pm js

    ! JF: for treatment of sedimentation
    INTEGER  :: elev, elev_lim           !< vertical end level
    LOGICAL  :: llbc_no_flux             !< TRUE: apply 'no flux' lower boundary condition

    LOGICAL  :: l_out_edgeval            !< corresponding local variable; default 
                                         !< .FALSE. i.e. output flux across the edge

    INTEGER  :: ist                      !< status variable
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER  :: i_rlstart, i_rlend

    INTEGER  :: n                        !< loop index

    REAL(wp) :: wsign                    !< wind direction: introduced, in order to merge flux formula  
                                         !< for w>0 and w<0 into one.
                                         !< +1, if w >0
                                         !< -1, if w <0

    REAL(wp) ::   &                      !< absolute CFL at grid point
      &  abs_cfl(nproma)

    REAL(wp) ::   &                      !< maximum CFL for each layer (blockwise)
      &  max_cfl_lay(p_patch%nlevp1,p_patch%nblks_c)

    REAL(wp) ::   &                      !< maximum CFL for each layer  
      &  max_cfl_lay_tot(p_patch%nlevp1)

    REAL(wp) ::    &                     !< maximum vertical Courant number per block
      &  max_cfl_blk(p_patch%nblks_c)

    REAL(wp) ::   &                      !< domain-wide maximum CFL
      &  max_cfl_tot

    REAL(wp) ::  &
      &  zq_ubc(nproma,p_patch%nblks_c)

    REAL(wp) :: rdtime                   !< 1/dt


    !-----------------------------------------------------------------------

    ! inverse of time step for computational efficiency
    rdtime = 1._wp/p_dtime

    ! check optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev  = opt_slev
      slevp1= opt_slev + 1
    ELSE
      slev  = 1
      slevp1= 2
    END IF

    ! check optional arguments
    IF ( PRESENT(opt_ti_slev) ) THEN
      slev_ti  = opt_ti_slev
      slevp1_ti= opt_ti_slev + 1
    ELSE
      slev_ti  = 1
      slevp1_ti= 2
    END IF

    IF ( PRESENT(opt_lout_edge) ) THEN
      l_out_edgeval = opt_lout_edge
    ELSE
      l_out_edgeval = .FALSE.
    ENDIF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = grf_bdywidth_c-1
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rlcell_int
    ENDIF

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! check optional arguments
    llbc_no_flux = .TRUE.
    IF ( PRESENT(opt_elev) ) THEN
      IF ( opt_elev == nlevp1 ) THEN
        elev = nlevp1
        llbc_no_flux = .FALSE.  ! e.g. when used for sedimentation
      ELSE
        elev = nlev
      END IF
    ELSE
      elev = nlev
    END IF
    elev_lim = nlev

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    IF ( ld_compute ) THEN
      !
      ! allocate field for storing the density weighted Courant number 
      !
      ! DA: this is a performance issue
      ! TODO: figure out z_cfl lifetime
      ALLOCATE( z_cfl(nproma,nlevp1,p_patch%nblks_c), STAT=ist  )
      IF (ist /= SUCCESS) THEN
        CALL finish(routine, 'allocation for z_cfl failed')
      ENDIF
      !$ACC ENTER DATA CREATE(z_cfl) IF(i_am_accel_node)
    END IF

    !$ACC DATA CREATE(z_face, z_face_up, z_face_low, z_delta_q, z_a1, zq_ubc) &
    !$ACC   IF(i_am_accel_node)

    IF ( PRESENT(opt_q_ubc) ) THEN
      !$ACC KERNELS DEFAULT(PRESENT) PRESENT(opt_q_ubc) ASYNC(1) IF(i_am_accel_node)
      zq_ubc(:,:) = opt_q_ubc(:,:)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      zq_ubc(:,:) = 0._wp
      !$ACC END KERNELS
    ENDIF


!$OMP PARALLEL

!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ikp1,jks,z_mass, &
!$OMP            jk_shift,js,n,z_iflx,z_delta_q,z_a1,z_q_int,  &
!$OMP            wsign,z_cflfrac,z_face,z_face_up,z_face_low   ) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )


      ! The contravariant mass flux should never exactly vanish
      !
      IF (l_out_edgeval) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slevp1, elev
          DO jc = i_startidx, i_endidx
            p_mflx_contra_v(jc,jk,jb) =                            &
          &              p_mflx_contra_v(jc,jk,jb)               &
          &              + SIGN(dbl_eps,p_mflx_contra_v(jc,jk,jb))
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF

      !
      ! 1. Compute density weighted Courant number for w<0 and w>0. 
      !    It is the sum of the fractional Courant number and the integer shift s.
      !    Stored at cell faces
      !
      IF (ld_compute) THEN

        ! initialize Courant number
        !$ACC KERNELS DEFAULT(PRESENT) PRESENT(z_cfl) ASYNC(1) IF(i_am_accel_node)
        z_cfl(i_startidx:i_endidx,slev_ti:nlevp1,jb) = 0._wp
        !$ACC END KERNELS


        ! Split density-weighted Courant number into integer and fractional 
        ! part and store the sum in z_cfl (for w>0 and w<0)
        !
        !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(z_cfl) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR PRIVATE(z_mass, jks, z_cflfrac) COLLAPSE(2)
        DO jk = slevp1_ti, elev
          DO jc = i_startidx, i_endidx

            ! total mass crossing jk'th edge during \Delta t
            z_mass = p_dtime*p_mflx_contra_v(jc,jk,jb)

            !
            ! Case: w > 0
            !
            IF (z_mass > 0._wp) THEN

              jks = jk   ! initialize shifted index

              DO WHILE( (z_mass >= p_cellmass_now(jc,jks,jb)) .AND. (jks <= nlev-1) )
                z_mass = z_mass - p_cellmass_now(jc,jks,jb)
                jks = jks+1
                ! update Courant number
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + 1._wp
              ENDDO

              z_cflfrac = z_mass/p_cellmass_now(jc,jks,jb)

              ! now we add the fractional Courant number
              !
              IF (z_cflfrac < 1._wp) THEN
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + z_cflfrac
              ELSE  ! z_cflfrac >= 1._wp
                !
                ! z_cflfrac >= 1 indicates that the backward trajectory hit the
                ! lower boundary and the previous while loop was exited ahead of time.
                ! In this case we limit the fractional Courant number to 1-epsilon.
                !
                ! Subtracting epsilon is essential in order to avoid out-of-bounds
                ! memory access at the lower boundary during the computation of fractional
                ! fluxes. The underlying reason is that we do not store the shifted index
                ! jks for later use (see above), but use FLOOR(z_cfl) to recompute it.
                !
                ! Note that this limitation of the Courant number essentially limits
                ! the mass flux across the edge considered, and hence, results in
                ! a local (and potentially rare) violation of the tracer and air mass
                ! consistency.
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + (1._wp - dbl_eps)
              ENDIF

            ELSE
            !
            ! Case w < 0
            !
              jks = jk-1   ! initialize shifted index

              DO WHILE( (ABS(z_mass) >= p_cellmass_now(jc,jks,jb)) .AND. &
                &       (jks >= slevp1_ti) )
                z_mass = z_mass + p_cellmass_now(jc,jks,jb)
                jks = jks-1
                ! update Courant number
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) - 1._wp
              ENDDO


              z_cflfrac = z_mass/p_cellmass_now(jc,jks,jb)

              ! now we add the fractional Courant number
              !
              IF (ABS(z_cflfrac) < 1._wp) THEN
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + z_cflfrac
              ELSE  ! ABS(z_cflfrac) >= 1._wp
                !
                ! ABS(z_cflfrac) >= 1 indicates that the backward trajectory hit the
                ! upper boundary and the previous while loop was exited ahead of time.
                ! In this case we limit the fractional Courant number to -1+epsilon.
                !
                ! Adding epsilon is essential in order to avoid out-of-bounds
                ! memory access at the uppper boundary during the computation of fractional
                ! fluxes. The underlying reason is that we do not store the shifted index
                ! jks for later use (see above), but use FLOOR(z_cfl) to recompute it.
                !
                ! Note that this limitation of the Courant number essentially limits
                ! the mass flux across the edge considered, and hence, results in
                ! a local (and potentially rare) violation of the tracer and air mass
                ! consistency.
                z_cfl(jc,jk,jb) = z_cfl(jc,jk,jb) + (dbl_eps - 1._wp)
              ENDIF

            ENDIF

          ENDDO  ! jc
        ENDDO  ! jk
        !$ACC END PARALLEL

      END IF ! ld_compute



      !
      ! 2. Edge value reconstruction
      !
      SELECT CASE(ivadv_tracer)
      CASE (IPPM_V)

        !
        ! PPM Reconstruction following Colella and Woodward (1984)
        !
        CALL compute_face_values_ppm( i_startidx       = i_startidx,               & !in
          &                           i_endidx         = i_endidx,                 & !in
          &                           slev             = slev,                     & !in
          &                           elev             = nlev,                     & !in
          &                           p_itype_vlimit   = p_itype_vlimit,           & !in
          &                           p_cc             = p_cc(:,:,jb),             & !in
          &                           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,jb), & !in
          &                           p_face           = z_face(:,:)               ) !inout


      CASE (IPSM_V)


        !
        ! PSM Reconstruction following Zerroukat et al. (2006)
        !
        CALL compute_face_values_psm( i_startidx       = i_startidx,               & !in
          &                           i_endidx         = i_endidx,                 & !in
          &                           slev             = slev,                     & !in
          &                           elev             = nlev,                     & !in
          &                           p_itype_vlimit   = p_itype_vlimit,           & !in
          &                           p_cc             = p_cc(:,:,jb),             & !in
          &                           p_cellhgt_mc_now = p_cellhgt_mc_now(:,:,jb), & !in
          &                           p_face           = z_face(:,:)               ) !inout


      END SELECT  ! ivadv_tracer



      !
      ! 4. Limitation/filtering of first guess parabola (which is based on z_face)
      ! Note that z_face_up(k) does not need to equal z_face_low(k-1) after
      ! the limitation procedure.
      ! Therefore 2 additional fields z_face_up and z_face_low are introduced.
      !
      SELECT CASE (p_itype_vlimit)
      CASE(ISLOPEL_VSM)

        ! semi-monotonic (sm) filter
        CALL v_limit_parabola_sm( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in

      CASE(ISLOPEL_VM)

        ! monotonic (mo) filter
        CALL v_limit_parabola_mo( p_ivlimit_selective,              & !in
          &                   p_cc(:,:,jb), z_face(:,:),            & !in
          &                   z_face_up(:,:), z_face_low(:,:),      & !inout
          &                   i_startidx, i_endidx, slev, elev_lim  ) !in


      CASE default

        ! simply copy face values to 'face_up' and 'face_low' arrays
        !
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR PRIVATE(ikp1) COLLAPSE(2)
        DO jk = slev, nlev
          DO jc = i_startidx, i_endidx
            ! index of bottom half level
            ikp1 = jk + 1
            z_face_up(jc,jk)  = z_face(jc,jk)
            z_face_low(jc,jk) = z_face(jc,ikp1)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      END SELECT  ! p_itype_vlimit


      !
      ! 5. calculation of upwind fluxes. For CFL>1, the total flux is the sum of
      !    integer-fluxes and a fractional flux. IF CFL<=1 the fluxes are only
      !    comprised of the fractional flux. The fractional flux is calculated
      !    by assuming a piecewise parabolic approx. for the subgrid distribution.
      !

      ! 5a. Compute coefficients of reconstructed parabola as they are used below.
      !     Terminology follows Colella (1984)
      !     z_delta_q = 0.5*\Delta q
      !     z_a1 = 1/6*a_6
      !
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = slev, nlev
        DO jc = i_startidx, i_endidx
          z_delta_q(jc,jk) = 0.5_wp * (z_face_up(jc,jk) - z_face_low(jc,jk))
          z_a1(jc,jk)      = p_cc(jc,jk,jb) - 0.5_wp*(z_face_up(jc,jk) + z_face_low(jc,jk))
        ENDDO
      ENDDO
      !$ACC END PARALLEL


      !
      ! 5b. First compute the fractional fluxes for all cell faces.
      !     For cell faces with CFL>1, integer fluxes will be added lateron.
      !
      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(z_cfl) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR PRIVATE(js, z_cflfrac, jks, wsign, z_q_int) COLLAPSE(2)
      DO jk = slevp1, elev
        DO jc = i_startidx, i_endidx

          ! get integer shift (always non-negative)
          js = FLOOR(ABS(z_cfl(jc,jk,jb)))

          ! get fractional part of Courant number (always non-negative)
          z_cflfrac = ABS(z_cfl(jc,jk,jb)) - REAL(js,wp)

          ! compute shifted cell index
          IF (z_cfl(jc,jk,jb) > 0._wp) THEN
            jks = MIN(jk,nlev)+js
            wsign = 1._wp
          ELSE
            jks = (jk-1) - js
            wsign = -1._wp
          ENDIF

          ! this is needed in addition in order to avoid accessing non-existing (uninitalized)
          ! source levels for tracers that are not advected on all model levels
          IF (jks < slev) THEN
            p_upflux(jc,jk,jb) = 0._wp
            CYCLE
          ENDIF

          ! compute flux
          !
          ! flux formula differs between w>0 and w<0. 
          ! By using the coefficient 'wsign' we are able to merge 
          ! the two formula into one.
          z_q_int = p_cc(jc,jks,jb)                                 &
            &     + wsign*(z_delta_q(jc,jks) * (1._wp - z_cflfrac)) &
            &     - z_a1(jc,jks)*(1._wp - 3._wp*z_cflfrac + 2._wp*z_cflfrac*z_cflfrac)

          p_upflux(jc,jk,jb) = wsign * p_cellmass_now(jc,jks,jb)    &
            &                * z_cflfrac * z_q_int * rdtime
        ENDDO

      ENDDO ! end loop over vertical levels
      !$ACC END PARALLEL


      !
      ! 5c. Now compute the integer fluxes and add them to the fractional flux
      !
      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(z_cfl) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR PRIVATE(js, z_iflx, jk_shift) COLLAPSE(2)
      DO jk = slevp1, elev

        DO jc = i_startidx, i_endidx

          ! get integer shift (always non-negative)
          js = FLOOR(ABS(z_cfl(jc,jk,jb)))

          IF (js == 0) CYCLE   ! no work to do

          z_iflx = 0._wp

          ! case w > 0
          IF (z_cfl(jc,jk,jb) > 0._wp) THEN

            DO n = 1, js
              jk_shift = jk-1 + n
              ! Integer flux (division by p_dtime is done at the end)
              z_iflx = z_iflx + p_cc(jc,jk_shift,jb) * p_cellmass_now(jc,jk_shift,jb)
            ENDDO

          ! case w <= 0
          ELSE

            DO n = 1, js
              jk_shift = jk - n

              ! cycle if the source model level is in a region where advection is
              ! turned off for the present variable
              IF (jk_shift < slev) CYCLE

              ! Integer flux (division by p_dtime is done at the end)
              z_iflx = z_iflx - p_cc(jc,jk_shift,jb) * p_cellmass_now(jc,jk_shift,jb)
            ENDDO

          ENDIF

          ! compute full (integer- plus high order fractional) flux
          p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb) + z_iflx*rdtime
        ENDDO  ! jc
      ENDDO ! jk
      !$ACC END PARALLEL


      !
      ! set upper and lower boundary condition
      !
      CALL set_bc_vadv(i_start      = i_startidx,                 & !in
        &              i_end        = i_endidx,                   & !in
        &              iubc_adv     = p_iubc_adv,                 & !in
        &              llbc_no_flux = llbc_no_flux,               & !in
        &              mflx_top     = p_mflx_contra_v(:,slev,jb), & !in
        &              q_top        = zq_ubc(:,jb),               & !in
        &              upflx_top    = p_upflux(:,slev,jb),        & !out
        &              upflx_bottom = p_upflux(:,nlevp1,jb))        !out


      ! If desired, get edge value of advected quantity
      IF ( l_out_edgeval ) THEN

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = slevp1, nlev
          DO jc = i_startidx, i_endidx
            p_upflux(jc,jk,jb) = p_upflux(jc,jk,jb)/p_mflx_contra_v(jc,jk,jb)
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ENDIF


      !
      ! 6. If desired, apply positive-definite flux limiter to limit 
      !    computed fluxes (based on work by Zalesak (1979)).
      !
      IF (p_itype_vlimit == IFLUXL_VPD) THEN
        ! positive-definite (pd) flux limiter
        CALL vflx_limiter_pd( p_dtime         = p_dtime,                & !in
          &                   p_cc            = p_cc(:,:,jb),           & !in
          &                   p_rhodz_now     = p_cellmass_now(:,:,jb), & !in
          &                   p_mflx_tracer_v = p_upflux(:,:,jb),       & !inout
          &                   i_startidx      = i_startidx,             & !in
          &                   i_endidx        = i_endidx,               & !in
          &                   slev            = slev,                   & !in
          &                   elev            = elev_lim                ) !in
      ENDIF

    ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL




! These diagnostics are not viable on GPU
    !
    ! If desired, print maximum vertical CFL number
    !
    IF ( ld_compute .AND. msg_level >= 10 .AND. lprint_cfl ) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx,abs_cfl)
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,       &
          &                 i_startidx, i_endidx, i_rlstart, i_rlend )

        DO jk = slevp1_ti, nlev
          DO jc = i_startidx, i_endidx
            abs_cfl(jc) = ABS(z_cfl(jc,jk,jb))
          ENDDO  ! jc
          max_cfl_lay(jk,jb) = MAXVAL(abs_cfl(i_startidx:i_endidx))
        ENDDO  ! jk
        !
        max_cfl_blk(jb) = MAXVAL(max_cfl_lay(slevp1_ti:elev,jb))
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      max_cfl_tot = MAXVAL(max_cfl_blk(i_startblk:i_endblk))

      ! Take maximum over all PEs
      IF (msg_level >= 13) THEN
        max_cfl_tot = global_max(max_cfl_tot)
      ELSE
        max_cfl_tot = global_max(max_cfl_tot, iroot=process_mpi_stdio_id)
      ENDIF
      IF (my_process_is_stdio() .OR. comm_lev>0 .AND. get_my_mpi_work_id() == get_glob_proc0() ) THEN
        ! otherwise it is possible that max_cfl_tot is undefined
        WRITE(message_text,'(a,e16.8)') 'maximum vertical CFL =',max_cfl_tot
        CALL message(routine,message_text)
      ENDIF

      ! Add layer-wise diagnostic if the maximum CFL value is close to the stability limit
      IF (msg_level >= 13 .AND. max_cfl_tot > 4._wp) THEN
        DO jk = slevp1_ti, nlev
          max_cfl_lay_tot(jk) = MAXVAL(max_cfl_lay(jk,i_startblk:i_endblk))
        ENDDO

        max_cfl_lay_tot(slevp1_ti:nlev) = global_max(max_cfl_lay_tot(slevp1_ti:nlev), iroot=process_mpi_stdio_id)
        DO jk = slevp1_ti,nlev
          WRITE(message_text,'(a,i4,a,e16.8)') 'maximum vertical CFL in layer', jk,' =', max_cfl_lay_tot(jk)
          CALL message(routine,message_text)
        ENDDO
      ENDIF

    END IF

    !$ACC WAIT
    !$ACC END DATA

    IF ( ld_cleanup ) THEN
      ! deallocate temporary arrays
      !$ACC EXIT DATA DELETE(z_cfl) IF(i_am_accel_node)
      DEALLOCATE( z_cfl, STAT=ist )
      IF (ist /= SUCCESS) THEN
        CALL finish(routine, 'deallocation for z_cfl failed' )
      ENDIF
    END IF

  END SUBROUTINE upwind_vflux_ppm4gpu




  !---------------------------------------------------------------
  !! Description:
  !!   solve the vertical flux advection equation for sedimentation
  !!   of scalar variables (a purely downward directed transport)
  !!   with the 2-point implicit scheme described in
  !!   COSMO Sci. doc. II, section 5.2.4.
  !!
  !! Method:
  !!   index convention for the sedimentation velocity :
  !!   v_new(i,j,k) = v(i,j,k-1/2)
  !!   sign: v_new, v_old > 0 ! (i.e. directed downward)
  !!
  !!   negative values in phi_new are clipped; this destroys
  !!   mass conservation.
  !
  SUBROUTINE implicit_sedim_tracer( tracer,                    &
    &                        rho, rho_inv,                     &
    &                        v_new, v_old,                     &
    &                        dt,                               &
    &                        p_patch, p_metrics,               &
    &                        i_rlstart, i_rlend,               &
    &                        rhoS, lacc )

    USE mo_parallel_config,         ONLY: nproma
    USE mo_nonhydro_types ,         ONLY: t_nh_metrics

    IMPLICIT NONE

    REAL (wp),     INTENT(INOUT) :: tracer(:,:,:) !advected cell centered variable
                                                       !< dim: (nproma, nlev, nblks_c)
    REAL (wp),     INTENT(IN)  :: rho    (:,:,:)  ! mass density  (in kg/m^3)
    REAL (wp),     INTENT(IN)  :: rho_inv(:,:,:)  ! 1/rho  (in m^3/kg)
    REAL (wp),     INTENT(IN)  :: v_new  (:,:,:)  ! sedimentation velocity at time level n+1 (in m/s)
    REAL (wp),     INTENT(IN)  :: v_old  (:,:,:)  ! sedimentation velocity at time level n   (in m/s)
    REAL (wp),     INTENT(IN)  :: dt     ! time step

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch      !< Patch on which computation is performed
    TYPE(t_nh_metrics), TARGET,  INTENT(IN) :: p_metrics    !< Metrical fields

    INTEGER, INTENT(IN) :: i_rlstart, i_rlend

    REAL (wp), INTENT(IN), OPTIONAL :: rhoS(:,:,:)

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc

    INTEGER    :: jk, jb, jc
    REAL (wp)  :: h, dz, c, lambda_im

    REAL (wp) :: phi_old( nproma, p_patch%nlev )
    REAL (wp) :: phi_new( nproma, p_patch%nlev )

    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_startidx, i_endidx

    LOGICAL :: is_rhoS_present

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    is_rhoS_present = PRESENT( rhoS )   ! own variable may help in vectorisation for some compilers

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    !$ACC DATA CREATE(phi_new, phi_old) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,dz,c,lambda_im,h,phi_old,phi_new) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,        &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend )

      ! calculate densities (for the following flux advection scheme)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx
          phi_old(jc,jk) = tracer(jc,jk,jb) * rho(jc,jk,jb)
        ENDDO ! jc
      ENDDO ! jk
      !$ACC END PARALLEL

      IF ( is_rhoS_present ) THEN

        ! top level
        jk = 1
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(dz, c, lambda_im, h)
        DO jc = i_startidx, i_endidx

          dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
          c = dt / ( 2.0_wp * dz );
          lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

          h = phi_old(jc,jk) - c * ( v_old(jc,jk+1,jb) * phi_old(jc,jk) )

          phi_new(jc,jk) = MAX( lambda_im * ( h + rhoS(jc,jk,jb)*dt ), 0.0_wp)
        END DO ! jc
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP SEQ
        DO jk=2, p_patch%nlev
          !$ACC LOOP GANG VECTOR PRIVATE(dz, c, lambda_im, h)
          DO jc = i_startidx, i_endidx

            dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
            c = dt / ( 2.0_wp * dz );
            lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

            h = phi_old(jc,jk) + c *                      &
              &  ( v_new(jc,jk  ,jb) * phi_new(jc,jk-1)   &
              &  + v_old(jc,jk  ,jb) * phi_old(jc,jk-1)   &
              &  - v_old(jc,jk+1,jb) * phi_old(jc,jk  )  )

            phi_new(jc,jk) = MAX( lambda_im * ( h + rhoS(jc,jk,jb)*dt ), 0.0_wp)

          END DO ! jc
        END DO ! jk
        !$ACC END PARALLEL

      ELSE

        ! the same code as in the if-block before but without the source term rhoS:

        ! top level
        jk = 1
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR PRIVATE(dz, c, lambda_im, h)
        DO jc = i_startidx, i_endidx

          dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
          c = dt / ( 2.0_wp * dz );
          lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

          h = phi_old(jc,jk) - c * ( v_old(jc,jk+1,jb) * phi_old(jc,jk) )

          phi_new(jc,jk) = MAX( lambda_im * h, 0.0_wp)
        END DO ! jc
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP SEQ
        DO jk=2, p_patch%nlev
          !$ACC LOOP GANG VECTOR PRIVATE(dz, c, lambda_im, h)
          DO jc = i_startidx, i_endidx

            dz = p_metrics%z_ifc(jc,jk,jb) - p_metrics%z_ifc(jc,jk+1,jb)
            c = dt / ( 2.0_wp * dz );
            lambda_im = 1.0_wp / ( 1.0_wp + c * v_new(jc,jk,jb) )

            h = phi_old(jc,jk) + c *                      &
              &  ( v_new(jc,jk  ,jb) * phi_new(jc,jk-1)   &
              &  + v_old(jc,jk  ,jb) * phi_old(jc,jk-1)   &
              &  - v_old(jc,jk+1,jb) * phi_old(jc,jk  )  )

            phi_new(jc,jk) = MAX( lambda_im * h, 0.0_wp)
          END DO ! jc
        END DO ! jk
        !$ACC END PARALLEL

      END IF       ! IF ( is_rhoS_present )

      ! calculate back the specific mass:
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx
          tracer(jc,jk,jb) = phi_new(jc,jk) * rho_inv(jc,jk,jb)
        ENDDO ! jc
      ENDDO ! jk
      !$ACC END PARALLEL

    END DO   ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE implicit_sedim_tracer



  !-------------------------------------------------------------------------
  !! Set top and bottom  boundary condition for vertical transport
  !!
  SUBROUTINE set_bc_vadv(i_start, i_end, iubc_adv, llbc_no_flux, mflx_top, q_top, &
    &                        upflx_top, upflx_bottom )

!!$    CHARACTER(len=*), PARAMETER :: routine = modname//':set_ubc_adv'

    INTEGER, INTENT(IN)      :: & !< start and end index
      &  i_start, i_end
    INTEGER, INTENT(IN)      :: & !< selects upper boundary condition
      &  iubc_adv
    LOGICAL, INTENT(IN)      :: & !< TRUE: apply 'no flux' lower boundary condition
      &  llbc_no_flux
    REAL(wp), INTENT(IN)     :: & !< mass flux at upper boundary
      &  mflx_top(:)
    REAL(wp), INTENT(IN)     :: & !< tracer mass fraction at upper boundary
      &  q_top(:)
    REAL(wp), INTENT(OUT)    :: & !< upper boundary condition
      &  upflx_top(:)
    REAL(wp), INTENT(INOUT)  :: & !< lower boundary condition
      &  upflx_bottom(:)

    INTEGER:: jc
    !-------------------------------------------------------------------------

    !
    ! flux at top boundary
    !
    SELECT CASE (iubc_adv)
    CASE ( ino_flx )

      ! no flux condition

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = i_start,i_end
        upflx_top(jc) = 0._wp
      ENDDO
      !$ACC END PARALLEL
    CASE ( iparent_flx ) ! interpolated flux from parent grid

      ! multiply given face value q_ubc with time averaged mass flux
      ! at (nest) upper boundary

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = i_start,i_end
        upflx_top(jc) = q_top(jc)*mflx_top(jc)
      ENDDO
      !$ACC END PARALLEL
    END SELECT

    !
    ! flux at bottom boundary
    !
    IF ( llbc_no_flux ) THEN

      ! no flux condition

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = i_start,i_end
        upflx_bottom(jc) = 0._wp
      ENDDO
      !$ACC END PARALLEL
    END IF

  END SUBROUTINE set_bc_vadv


  !-------------------------------------------------------------------------
  !! PSM Face value reconstruction after Zerroukat et al (2006)
  !!
  !
  SUBROUTINE compute_face_values_psm( i_startidx, i_endidx, slev, elev, &
    &                                 p_itype_vlimit, p_cc, p_cellhgt_mc_now, p_face )

    INTEGER, INTENT(IN) ::    &    !< horizontal start and end indices
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN) ::    &    !< vertical start and end levels
      &  slev, elev

    INTEGER, INTENT(IN) ::    &    !< selects the vertical limiter
      &  p_itype_vlimit

    REAL(wp), INTENT(IN) ::   &    !< cell centered variable (cell average)
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::   &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:)     !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT):: &    !< face values of transported field
      &  p_face(:,:)               !< dim: (nproma,nlev)


    ! local variables
    INTEGER :: jc, jk
    INTEGER :: elevp1              !< end level + 1

    ! TDMA arrays
    REAL(wp) :: a(SIZE(p_face,1),SIZE(p_face,2))   !< sub-diagonal
    REAL(wp) :: b(SIZE(p_face,1),SIZE(p_face,2))   !< main diagonal
    REAL(wp) :: c(SIZE(p_face,1),SIZE(p_face,2))   !< super diagonal
    REAL(wp) :: rhs(SIZE(p_face,1),SIZE(p_face,2)) !< right hand side
    REAL(wp) :: dzfrac                             !< ratio of neighboring cell heights

    !-------------------------------------------------------------------------

    elevp1 = elev + 1

    !$ACC DATA CREATE(a, b, c, rhs) IF(i_am_accel_node)

    !
    ! 1. reconstruct face values at vertical half-levels using splines
    !
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx
      ! top BC
      a  (jc,slev) = 0._wp
      b  (jc,slev) = 2._wp
      c  (jc,slev) = 1._wp
      rhs(jc,slev) = 3._wp*p_cc(jc,slev)
      !
      ! bottom BC
      a  (jc,elevp1) = 1._wp
      b  (jc,elevp1) = 2._wp
      c  (jc,elevp1) = 0._wp
      rhs(jc,elevp1) = 3._wp*p_cc(jc,elev)
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR PRIVATE(dzfrac) COLLAPSE(2)
    DO jk=slev+1,elev
      DO jc = i_startidx, i_endidx
        dzfrac = p_cellhgt_mc_now(jc,jk-1)/p_cellhgt_mc_now(jc,jk)
        a  (jc,jk) = 1._wp
        c  (jc,jk) = dzfrac
        b  (jc,jk) = 2._wp*(1._wp + dzfrac)
        rhs(jc,jk) = 3._wp*(dzfrac*p_cc(jc,jk) + p_cc(jc,jk-1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    ! solve tri-diagonal system
    CALL tdma_solver_vec(a       = a,          &
      &                  b       = b,          &
      &                  c       = c,          &
      &                  d       = rhs,        &
      &                  slev    = slev,       &
      &                  elev    = elevp1,     &
      &                  startidx= i_startidx, &
      &                  endidx  = i_endidx,   &
      &                  varout  = p_face      )  ! out



    ! 2. OPTIONAL: Limit face values
    !
    SELECT CASE (p_itype_vlimit)
    CASE(ISLOPEL_VSM)

      ! make sure that face values are bounded by the neighbouring cell averages
      !
      CALL v_limit_face_mc_sm(p_cc             = p_cc(:,:),             & ! in
        &                     p_cellhgt_mc_now = p_cellhgt_mc_now(:,:), & ! in
        &                     p_face           = p_face(:,:),           & ! inout
        &                     i_startidx       = i_startidx,            & ! in
        &                     i_endidx         = i_endidx,              & ! in
        &                     slev             = slev,                  & ! in
        &                     elev             = elev                   ) ! in


      ! top and bottom face
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        p_face(jc,slev)   = MAX(p_cc(jc,slev),p_face(jc,slev))
        p_face(jc,elevp1) = MAX(p_cc(jc,elev),p_face(jc,elevp1))
      ENDDO
      !$ACC END PARALLEL


    CASE(ISLOPEL_VM)

      ! make sure that face values are bounded by the neighbouring cell averages
      !
      CALL v_limit_face_mc_mo(p_cc             = p_cc(:,:),             & ! in
        &                     p_cellhgt_mc_now = p_cellhgt_mc_now(:,:), & ! in
        &                     p_face           = p_face(:,:),           & ! inout
        &                     i_startidx       = i_startidx,            & ! in
        &                     i_endidx         = i_endidx,              & ! in
        &                     slev             = slev,                  & ! in
        &                     elev             = elev                   ) ! in
 

      ! top and bottom face
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        p_face(jc,slev)   = p_cc(jc,slev)
        p_face(jc,elevp1) = p_cc(jc,elev)
      ENDDO
      !$ACC END PARALLEL

    CASE default
      ! do nothing
    END SELECT

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE compute_face_values_psm


  !-------------------------------------------------------------------------
  !! PPM Face value reconstruction after Colella and Woodward (1984)
  !!
  !
  SUBROUTINE compute_face_values_ppm( i_startidx, i_endidx, slev, elev, &
    &                                 p_itype_vlimit, p_cc, p_cellhgt_mc_now, p_face )

    INTEGER, INTENT(IN) ::    &    !< horizontal start and end indices
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN) ::    &    !< vertical start and end levels
      &  slev, elev

    INTEGER, INTENT(IN) ::    &    !< selects the vertical limiter
      &  p_itype_vlimit

    REAL(wp), INTENT(IN) ::   &    !< cell centered variable (cell average)
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::   &    !< layer thickness at cell center at time n
      &  p_cellhgt_mc_now(:,:)     !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT):: &    !< face values of transported field
      &  p_face(:,:)               !< dim: (nproma,nlev)


    ! local variables
    INTEGER :: jc, jk
    INTEGER :: ikm1, ikp1, ikm2      !< vertical level minus and plus one, minus two
    INTEGER :: slevp1, elevp1        !< start/end level + 1

    REAL(wp) ::   &                  !< auxiliaries for optimization
      &  zfac, zfac_m1

    REAL(wp) ::   &                  !< auxiliary field for optimization
      &  zfac_n(nproma)

    REAL(wp) ::   &                  !< geometric factors
      &  zgeo1, zgeo2, zgeo3, zgeo4

    REAL(wp) :: &                     !< (monotonized) slope
      &  z_slope(SIZE(p_cc,1),SIZE(p_cc,2))

    !-------------------------------------------------------------------------

    slevp1 = slev + 1
    elevp1 = elev + 1

    !$ACC DATA CREATE(z_slope) IF(i_am_accel_node)

    !
    ! 1. Compute slope
    !
    ! Initialize z_slope for jk=slev
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    z_slope(i_startidx:i_endidx,slev) = 0._wp
    !$ACC END KERNELS

    ! Initialize zfac_n for jk=slev
    zfac_n(i_startidx:i_endidx) = (p_cc(i_startidx:i_endidx,slevp1) - p_cc(i_startidx:i_endidx,slev)) &
      &                         /( p_cellhgt_mc_now(i_startidx:i_endidx,slevp1)                       &
      &                          + p_cellhgt_mc_now(i_startidx:i_endidx,slev) )

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR PRIVATE(ikm1, ikp1, zfac, zfac_m1) COLLAPSE(2)
    DO jk = slevp1, elev
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, elev )

        zfac_m1 = zfac_n(jc)

        zfac = (p_cc(jc,ikp1) - p_cc(jc,jk)) &
          &  / (p_cellhgt_mc_now(jc,ikp1) + p_cellhgt_mc_now(jc,jk))

        z_slope(jc,jk) = ( p_cellhgt_mc_now(jc,jk)                                       &
          &  / (p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk)                      &
          &  + p_cellhgt_mc_now(jc,ikp1)) )                                              &
          &  * ( (2._wp * p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk)) * zfac    &
          &  + (p_cellhgt_mc_now(jc,jk) + 2._wp * p_cellhgt_mc_now(jc,ikp1)) * zfac_m1)

        zfac_n(jc) = zfac
      END DO  ! jc

    END DO  ! jk
    !$ACC END PARALLEL

    !
    ! 2. Optional: monotonize slope if necessary
    !     - only necessary, when using the monotonic or semi-monotonic 
    !       sub-grid scale filter (i.e. if the parabola is modified). 
    !
    IF (p_itype_vlimit == ISLOPEL_VSM) THEN
      CALL v_limit_slope_sm(p_cc(:,:), i_startidx, i_endidx, slevp1, elev, z_slope)
    ELSE IF (p_itype_vlimit == ISLOPEL_VM) THEN
      CALL v_limit_slope_mo(p_cc(:,:), i_startidx, i_endidx, slevp1, elev, z_slope)
    ENDIF



    !
    ! 3. face value reconstruction at vertical half-levels
    !

    ! Boundary values for two uppermost and lowermost half-levels
    !
    ! for faces k=slevp1 and k=nlevp1-1 face values are reconstructed by
    ! interpolating a quadratic (instead of quartic) polynomial through 3
    ! values of the indefinite integral A=\int_{z_{0}}^{z}q\,\mathrm{d}z
    !
    ! for faces k=slev and k=nlevp1 a zero gradient condition is assumed and the
    ! face values are set to the values of the corresponding cell centers
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR
    DO jc = i_startidx, i_endidx

      ! see transport documentation for derivation
      p_face(jc,slevp1) = p_cc(jc,slevp1)*(1._wp - (p_cellhgt_mc_now(jc,slevp1) &
        &       / p_cellhgt_mc_now(jc,slev))) + (p_cellhgt_mc_now(jc,slevp1)    &
        &       /(p_cellhgt_mc_now(jc,slev) + p_cellhgt_mc_now(jc,slevp1)))     &
        &       * ((p_cellhgt_mc_now(jc,slevp1) / p_cellhgt_mc_now(jc,slev))    &
        &       * p_cc(jc,slevp1) + p_cc(jc,slev))

      p_face(jc,elev) = p_cc(jc,elev)*( 1._wp                                &
        &       - (p_cellhgt_mc_now(jc,elev) / p_cellhgt_mc_now(jc,elev-1))) &
        &       + (p_cellhgt_mc_now(jc,elev)/(p_cellhgt_mc_now(jc,elev-1)    &
        &       + p_cellhgt_mc_now(jc,elev))) * ((p_cellhgt_mc_now(jc,elev)  &
        &       / p_cellhgt_mc_now(jc,elev-1)) * p_cc(jc,elev)               &
        &       + p_cc(jc,elev-1))

      p_face(jc,slev)   = p_cc(jc,slev)
      p_face(jc,elevp1) = p_cc(jc,elev)

    ENDDO
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(i_am_accel_node)
    !$ACC LOOP GANG VECTOR PRIVATE(ikm1, ikm2, ikp1, zgeo1, zgeo2, zgeo3, zgeo4) COLLAPSE(2)
    DO jk = slev+2, elev-1
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1 = jk - 1
        ikm2 = jk - 2
        ! index of bottom half level
        ikp1 = jk + 1

        zgeo1 = p_cellhgt_mc_now(jc,ikm1)                                      &
          &   / (p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk))
        zgeo2 = 1._wp / (p_cellhgt_mc_now(jc,ikm2) + p_cellhgt_mc_now(jc,ikm1) &
          &   + p_cellhgt_mc_now(jc,jk) + p_cellhgt_mc_now(jc,ikp1))
        zgeo3 = (p_cellhgt_mc_now(jc,ikm2) + p_cellhgt_mc_now(jc,ikm1))        &
          &   / (2._wp*p_cellhgt_mc_now(jc,ikm1) + p_cellhgt_mc_now(jc,jk))
        zgeo4 = (p_cellhgt_mc_now(jc,ikp1) + p_cellhgt_mc_now(jc,jk))          &
          &   / (2._wp*p_cellhgt_mc_now(jc,jk) + p_cellhgt_mc_now(jc,ikm1))


        p_face(jc,jk) = p_cc(jc,ikm1)                                  &
          &  + zgeo1 * (p_cc(jc,jk) - p_cc(jc,ikm1))                   &
          &  + zgeo2 * ( (2._wp * p_cellhgt_mc_now(jc,jk) * zgeo1)     &
          &  * ( zgeo3 - zgeo4 ) * (p_cc(jc,jk) - p_cc(jc,ikm1))       &
          &  - zgeo3 * p_cellhgt_mc_now(jc,ikm1)   * z_slope(jc,jk)    &
          &  + zgeo4 * p_cellhgt_mc_now(jc,jk) * z_slope(jc,ikm1) )

      END DO  !jc

    END DO  !jk
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE compute_face_values_ppm


END MODULE mo_advection_vflux
