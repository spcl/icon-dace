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

! This module contains utility routines needed for the initialization of the
! NH model

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_init_utils

  USE mo_kind,                  ONLY: wp, vp
  USE mo_model_domain,          ONLY: t_patch
  USE mo_parallel_config,       ONLY: nproma
  USE mo_grid_config,           ONLY: l_limited_area
  USE mo_physical_constants,    ONLY: grav, cpd, rd, cvd_o_rd, p0ref
  USE mo_vertical_coord_table,  ONLY: vct_b
  USE mo_impl_constants,        ONLY: min_rlcell
  USE mo_sync,                  ONLY: sync_patch_array, SYNC_C
  USE mo_intp_data_strc,        ONLY: t_int_state
  USE mo_intp,                  ONLY: edges2cells_scalar
  USE mo_math_gradients,        ONLY: grad_fd_norm
  USE mo_loopindices,           ONLY: get_indices_c, get_indices_e
  USE mo_initicon_types,        ONLY: t_init_state
  USE mo_util_phys,             ONLY: virtual_temp
  USE mo_ifs_coord,             ONLY: geopot
  USE mo_fortran_tools,         ONLY: init

  IMPLICIT NONE

  PRIVATE


  ! subroutines
  !
  PUBLIC :: interp_uv_2_vn
  PUBLIC :: init_w, adjust_w
  PUBLIC :: convert_thdvars
  PUBLIC :: convert_omega2w
  PUBLIC :: compute_input_pressure_and_height
  PUBLIC :: compute_exner_pert

CONTAINS

  !-------------
  !> Compute pressure and height of input data.
  !
  !  OUT: initicon%const%z_mc_in
  !       initicon%atm_in%pres
  !
  SUBROUTINE compute_input_pressure_and_height(p_patch, psfc, phi_sfc, initicon, opt_lmask)
    TYPE(t_patch),          INTENT(IN)       :: p_patch
    REAL(wp),               INTENT(INOUT)    :: psfc(:,:)
    REAL(wp),               INTENT(INOUT)    :: phi_sfc(:,:)
    CLASS(t_init_state),    INTENT(INOUT)    :: initicon
    LOGICAL, OPTIONAL,      INTENT(IN)       :: opt_lmask(:,:)
    ! LOCAL VARIABLES
    INTEGER :: jb, nlen, nlev_in
    INTEGER :: jc, jc1, jb1
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev  ) :: delp, rdelp, rdlnpr, rdalpha, geop_mc
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_c) :: temp_v_in
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev+1) :: pres_ic, lnp_ic, geop_ic

    nlev_in = initicon%atm_in%nlev

    ! Compute virtual temperature of input data
    CALL virtual_temp(p_patch, initicon%atm_in%temp, initicon%atm_in%qv, initicon%atm_in%qc, &
                      initicon%atm_in%qi, initicon%atm_in%qr, initicon%atm_in%qs,            &
                      temp_v=temp_v_in)


    ! 1. Compute pressure and height of input data, using the IFS routines

    ! If mask field is provided, fill data-void points (mask=.FALSE.) 
    ! with dummy value.
    IF (PRESENT(opt_lmask)) THEN
      !
      ! Detect first grid point for which the mask field is .true.
      outer: DO jb = 1, p_patch%nblks_c
        IF (jb /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
         nlen = p_patch%npromz_c
        ENDIF
        inner: DO jc = 1, nlen
          IF (opt_lmask(jc,jb)) THEN
            jc1 = jc
            jb1 = jb
            EXIT outer
          ENDIF
        ENDDO inner
      ENDDO outer

      ! Do filling for psfc, phi_sfc, temp_v_in
!$OMP PARALLEL DO PRIVATE(jb,jc,nlen)
      DO jb = 1,p_patch%nblks_c

        IF (jb /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        ENDIF

        DO jc = 1, nlen
          IF (.NOT. opt_lmask(jc,jb)) THEN
            psfc(jc,jb)                = psfc(jc1,jb1)
            phi_sfc(jc,jb)             = phi_sfc(jc1,jb1)
            temp_v_in(jc,1:nlev_in,jb) = temp_v_in(jc1,1:nlev_in,jb1)
          ENDIF
        ENDDO
      ENDDO  ! jb
!$OMP END PARALLEL DO

    ENDIF


!$OMP PARALLEL
!$OMP DO PRIVATE(jb, nlen, pres_ic, lnp_ic, geop_ic, delp, rdelp, rdlnpr, &
!$OMP            rdalpha, geop_mc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1,p_patch%nblks_c

      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF
      
      ! Check if psfc is really psfc or LOG(psfc)
      IF (MAXVAL(psfc(1:nlen,jb)) <= 100._wp) THEN
        psfc(1:nlen,jb) = EXP(psfc(1:nlen,jb))
      ENDIF
      
      CALL initicon%const%vct%half_level_pressure(psfc(:,jb), nproma, nlen, nlev_in, pres_ic)
      
      CALL initicon%const%vct%full_level_pressure(pres_ic,nproma, nlen, nlev_in, initicon%atm_in%pres(:,:,jb))
      
      CALL initicon%const%vct%auxhyb(pres_ic, nproma, nlen, nlev_in,     & ! in
        delp, rdelp, lnp_ic, rdlnpr, rdalpha) ! out
      
      CALL geopot(temp_v_in(:,:,jb), rdlnpr, rdalpha, phi_sfc(:,jb), & ! in
        nproma, 1, nlen, nlev_in, geop_mc, geop_ic ) ! inout
      
      ! Compute 3D height coordinate field
      initicon%const%z_mc_in(1:nlen,1:nlev_in,jb) = geop_mc(1:nlen,1:nlev_in)/grav
      
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE compute_input_pressure_and_height




  !-------------
  !>
  !! SUBROUTINE convert_thdvars
  !! Converts the hydrostatic set of thermodynamic variables into the nonhydrostatic one
  !!
  !! Required input fields: pressure, virtual temperature
  !! Output: density, Exner pressure, virtual potential temperature
  !!
  SUBROUTINE convert_thdvars(p_patch, pres, temp_v, &
                             rho, exner, theta_v    )


    TYPE(t_patch), INTENT(IN) :: p_patch

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN) :: pres  (:,:,:) ! pressure (Pa)
    REAL(wp), INTENT(IN) :: temp_v(:,:,:) ! virtual temperature (K)

    ! Output fields (prognostic model variables) - all defined at full model levels
    REAL(wp), INTENT(OUT) :: rho(:,:,:)        ! density (kg/m**3)
    REAL(wp), INTENT(OUT) :: exner(:,:,:)      ! Exner pressure
    REAL(wp), INTENT(OUT) :: theta_v(:,:,:)    ! virtual potential temperature (K)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, jc
    INTEGER :: nlen, nlev

    nlev = p_patch%nlev

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, p_patch%nblks_c
      IF (jb /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          exner(jc,jk,jb)   = (pres(jc,jk,jb)/p0ref)**(rd/cpd)
          theta_v(jc,jk,jb) = temp_v(jc,jk,jb)/exner(jc,jk,jb)
          rho(jc,jk,jb)     = exner(jc,jk,jb)**cvd_o_rd*p0ref/rd/theta_v(jc,jk,jb)
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE convert_thdvars

  !-------------
  !>
  !! SUBROUTINE convert_omega2w
  !! Converts the hydrostatic vertical velocity (omega, Pa/s)
  !! into physical vertical velocity (m/s)
  !! Note: this routine has to be called on the input grid,
  !! where omega, pressure and temperature are not vertically staggered
  !!
  !! Required input fields: omega, pressure, temperature
  !! Output: vertical wind speed
  !!
  SUBROUTINE convert_omega2w(omega, w, pres, temp, nblks, npromz, nlev, opt_lmask)


    ! Input fields
    REAL(wp), INTENT(IN) :: omega (:,:,:) ! omega (Pa/s)
    REAL(wp), INTENT(IN) :: pres  (:,:,:) ! pressure (Pa)
    REAL(wp), INTENT(IN) :: temp  (:,:,:) ! virtual temperature (K)

    ! Output
    REAL(wp), INTENT(OUT) :: w(:,:,:)  ! vertical velocity (m/s)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlev       ! Number of model levels
    
    LOGICAL , INTENT(IN), OPTIONAL :: opt_lmask(:,:) ! logical mask of points to process

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, jc
    INTEGER :: nlen

    IF(PRESENT(opt_lmask)) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks
        IF (jb /= nblks) THEN
          nlen = nproma
        ELSE
          nlen = npromz
        ENDIF

        DO jk = 1, nlev
          DO jc = 1, nlen
            IF (opt_lmask(jc,jb)) THEN
              w(jc,jk,jb) = -rd*omega(jc,jk,jb)*temp(jc,jk,jb)/(grav*pres(jc,jk,jb))
            ELSE ! fill with dummy value
              w(jc,jk,jb) = 0._wp
            ENDIF
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ELSE ! not present opt_lmask

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, nblks
        IF (jb /= nblks) THEN
          nlen = nproma
        ELSE
          nlen = npromz
        ENDIF

        DO jk = 1, nlev
          DO jc = 1, nlen
            w(jc,jk,jb) = -rd*omega(jc,jk,jb)*temp(jc,jk,jb)/(grav*pres(jc,jk,jb))
          ENDDO
        ENDDO

      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    ENDIF ! PRESENT(opt_lmask)

  END SUBROUTINE convert_omega2w


  !-------------
  !>
  !! SUBROUTINE interp_uv_2_vn
  !! Interpolates u and v on cell points to vn on edge points
  !!
  !! Required input fields: u and v on cell points
  !! Output: vn on edge points
  !!
  SUBROUTINE interp_uv_2_vn(p_patch, p_int, u, v, vn )


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields - all defined at full model levels
    REAL(wp), INTENT(IN) :: u(:,:,:) ! zonal wind component on cell points (m/s)
    REAL(wp), INTENT(IN) :: v(:,:,:) ! meridional wind component on cell points (m/s)

    ! Output field (prognostic model variable) - defined at full model levels
    ! Intent (INOUT) because lateral nest boundaries cannot be filled here
    REAL(wp), INTENT(INOUT) :: vn(:,:,:)  ! edge-normal wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je
    INTEGER :: nlev, nblks_e, i_startblk,i_endblk, i_startidx,i_endidx
    REAL(wp) :: z_u, z_v

    INTEGER, DIMENSION(:,:,:), POINTER :: iidx, iblk

    nlev = p_patch%nlev
    nblks_e = p_patch%nblks_e

    iidx => p_patch%edges%cell_idx
    iblk => p_patch%edges%cell_blk

!$OMP PARALLEL PRIVATE(i_startblk,i_endblk)

    IF (l_limited_area .OR. p_patch%id > 1) THEN ! Fill outermost nest boundary

      i_startblk = p_patch%edges%start_blk(1,1)
      i_endblk   = p_patch%edges%end_blk(1,1)

!$OMP DO PRIVATE(jb,i_startidx,i_endidx,je,jk,z_u,z_v) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 1, 1)

        DO je = i_startidx, i_endidx
          IF (iidx(je,jb,1) >= 1 .AND. iblk(je,jb,1) >= 1) THEN
            DO jk = 1, nlev
              z_u = u(iidx(je,jb,1),jk,iblk(je,jb,1))
              z_v = v(iidx(je,jb,1),jk,iblk(je,jb,1))
              vn(je,jk,jb) =  z_u*p_patch%edges%primal_normal(je,jb)%v1 + &
                              z_v*p_patch%edges%primal_normal(je,jb)%v2
            END DO
          ELSE IF (iidx(je,jb,2) >= 1 .AND. iblk(je,jb,2) >= 1) THEN
            DO jk = 1, nlev
              z_u = u(iidx(je,jb,2),jk,iblk(je,jb,2))
              z_v = v(iidx(je,jb,2),jk,iblk(je,jb,2))
              vn(je,jk,jb) =  z_u*p_patch%edges%primal_normal(je,jb)%v1 + &
                              z_v*p_patch%edges%primal_normal(je,jb)%v2
            END DO
          ENDIF
        END DO

      END DO
!$OMP END DO
    ENDIF

    i_startblk = p_patch%edges%start_blk(2,1)

!$OMP DO PRIVATE(jb,jk,je,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE

      DO jb = i_startblk, nblks_e

        CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                           i_startidx, i_endidx, 2)

#ifdef __LOOP_EXCHANGE
        DO je = i_startidx, i_endidx
          DO jk = 1, nlev
#else
        DO jk = 1, nlev
          DO je = i_startidx, i_endidx
#endif

            vn(je,jk,jb) = p_int%c_lin_e(je,1,jb)                                              &
              *(u(iidx(je,jb,1),jk,iblk(je,jb,1))*p_patch%edges%primal_normal_cell(je,jb,1)%v1 &
              + v(iidx(je,jb,1),jk,iblk(je,jb,1))*p_patch%edges%primal_normal_cell(je,jb,1)%v2)&
              +            p_int%c_lin_e(je,2,jb)                                              &
              *(u(iidx(je,jb,2),jk,iblk(je,jb,2))*p_patch%edges%primal_normal_cell(je,jb,2)%v1 &
              + v(iidx(je,jb,2),jk,iblk(je,jb,2))*p_patch%edges%primal_normal_cell(je,jb,2)%v2 )

          ENDDO
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE interp_uv_2_vn

  !-------------
  !>
  !! SUBROUTINE init_w
  !! Initializes the vertical wind field based on the lower boundary condition
  !! w = v grad h and an empirical vertical decay function
  !! The discretization used here is simpler than that used in the dynamical core
  !! but is sufficient to avoid excessive generation of sound waves during the start phase
  !!
  !! Required input fields: vn, z_ifc
  !! Output: w
  !!
  SUBROUTINE init_w(p_patch, p_int, vn, z_ifc, w)


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields
    REAL(wp), INTENT(IN) :: vn(:,:,:)    ! edge-normal wind component (m/s)
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) ! height of half levels (m)

    ! Output field - defined at half model levels
    ! Intent (INOUT) because lateral nest boundaries cannot be filled here
    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je, jc, ktop
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, nshift, i_startblk, i_startidx, i_endidx

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c
    nshift  = p_patch%nshift_total

    ! In order to initialize w(1) = 0 except for vertical nesting
    IF (nshift == 0) THEN
      ktop = 2
    ELSE
      ktop = 1
    ENDIF

    ! Compute slope at edges
    CALL grad_fd_norm (z_ifc, p_patch, z_slope_e, 1, nlevp1)

    ! slope cannot be computed at outer boundary edges
    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      ! Extrapolation of vn to half levels is neglected here
      DO je = i_startidx, i_endidx
        z_wsfc_e(je,1,jb) = vn(je,nlev,jb)*z_slope_e(je,nlevp1,jb)
      ENDDO
    ENDDO

    CALL edges2cells_scalar(z_wsfc_e,p_patch,p_int%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)

    i_startblk = p_patch%cells%start_blk(2,1)

!$OMP PARALLEL
    ! First, initialize w with zero in order to avoid undefined nest boundary points
    CALL init(w(:,:,:), lacc=.FALSE.)
!$OMP BARRIER

    ! specify a reasonable initial vertical wind speed
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO
      DO jk = nlev, ktop, -1
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = z_wsfc_c(jc,1,jb)*vct_b(jk+nshift)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_w


  !-------------
  !>
  !! SUBROUTINE adjust_w
  !! Computes the lower boundary condition for w in a similar way as init_w,
  !! but the result is then merged with an already available vertical wind
  !! field provided by an external data source.
  !!
  !! Required input fields: vn, w, z_ifc
  !! Output: w
  !!
  SUBROUTINE adjust_w(p_patch, p_int, vn, z_ifc, w)


    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch
    TYPE(t_int_state),     INTENT(IN)   :: p_int

    ! Input fields
    REAL(wp), INTENT(IN) :: vn(:,:,:)    ! edge-normal wind component (m/s)
    REAL(wp), INTENT(IN) :: z_ifc(:,:,:) ! height of half levels (m)

    ! INOUT field - defined at half model levels
    REAL(wp), INTENT(INOUT) :: w(:,:,:)  ! vertical wind component (m/s)

    ! LOCAL VARIABLES
    INTEGER :: jb, jk, je, jc, ktop
    INTEGER :: nlev, nlevp1, nblks_e, nblks_c, nshift, i_startblk, i_startidx, i_endidx
    REAL(wp):: wfac

    REAL(wp) :: z_wsfc_e(nproma,1,p_patch%nblks_e) ! w at surface (edge points)
    REAL(wp) :: z_wsfc_c(nproma,1,p_patch%nblks_c) ! w at surface (cell points)
    REAL(wp) :: z_slope_e(nproma,p_patch%nlevp1,p_patch%nblks_e) ! slope at edges

    nlev    = p_patch%nlev
    nlevp1  = p_patch%nlevp1
    nblks_e = p_patch%nblks_e
    nblks_c = p_patch%nblks_c
    nshift  = p_patch%nshift_total

    ! In order to initialize w(1) = 0 except for vertical nesting
    IF (nshift == 0) THEN
      ktop = 2
    ELSE
      ktop = 1
    ENDIF

    ! Compute slope at edges
    CALL grad_fd_norm (z_ifc, p_patch, z_slope_e, 1, nlevp1)

    ! slope cannot be computed at outer boundary edges
    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, &
                         i_startidx, i_endidx, 2)

      ! Extrapolation of vn to half levels is neglected here
      DO je = i_startidx, i_endidx
        z_wsfc_e(je,1,jb) = vn(je,nlev,jb)*z_slope_e(je,nlevp1,jb)
      ENDDO
    ENDDO

    CALL edges2cells_scalar(z_wsfc_e,p_patch,p_int%e_inn_c,z_wsfc_c,&
                            1,1,opt_rlstart=2)

    i_startblk = p_patch%cells%start_blk(2,1)

    ! specify lower boundary condition and merge with w field provided on input
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jc,jk,wfac) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                         i_startidx, i_endidx, 2)

      ! Lower boundary condition
      DO jc = i_startidx, i_endidx
        w(jc,nlevp1,jb) = z_wsfc_c(jc,1,jb)
      ENDDO

      ! Merging of lower boundary condition with interpolated data
      DO jk = nlev, ktop, -1
        wfac = vct_b(jk+nshift)**2
        DO jc = i_startidx, i_endidx
          w(jc,jk,jb) = (1._wp-wfac)*w(jc,jk,jb) + wfac*z_wsfc_c(jc,1,jb)
        ENDDO
      ENDDO

      IF (nshift == 0) THEN
        ! Upper boundary condition and smooth transition below
        ! if domain is not vertically nested
        DO jc = i_startidx, i_endidx
          w(jc,1,jb) = 0._wp
          w(jc,2,jb) = 0.33_wp*w(jc,2,jb)
          w(jc,3,jb) = 0.66_wp*w(jc,3,jb)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE adjust_w


  !-------------
  !>
  !! Computes perturbation exner pressure by subtracting the
  !! exner reference state from the actual exner pressure.
  !!
  SUBROUTINE compute_exner_pert(exner, exner_ref, exner_pr, use_acc)
    REAL(wp), INTENT(IN)    :: exner(:,:,:)         !< exner pressure
    REAL(vp), INTENT(IN)    :: exner_ref(:,:,:)     !< exner reference state
    REAL(wp), INTENT(INOUT) :: exner_pr(:,:,:)      !< perturbation exner pressure
    LOGICAL,  INTENT(IN)    :: use_acc              !< if True, use openACC

    INTEGER :: i,j,k,ie,je,ke

    ie = SIZE(exner_pr, 1)
    je = SIZE(exner_pr, 2)
    ke = SIZE(exner_pr, 3)
!$OMP PARALLEL
#if (defined(_CRAYFTN) || defined(__INTEL_COMPILER))
!$OMP DO PRIVATE(i,j,k)
#else
!$OMP DO COLLAPSE(3) PRIVATE(i,j,k)
#endif
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(use_acc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO k = 1, ke
      DO j = 1, je
        DO i = 1, ie
          exner_pr(i,j,k) = exner(i,j,k) - REAL(exner_ref(i,j,k), wp)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE compute_exner_pert

END MODULE mo_nh_init_utils
