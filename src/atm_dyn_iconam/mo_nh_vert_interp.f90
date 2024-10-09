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

! This module contains routines for the vertical interpolation of
! atmospheric data provided by external analyses to the ICON grid

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nh_vert_interp

  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_metrics
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp,                ONLY: edges2cells_scalar, &
    &                               cells2edges_scalar, cells2verts_scalar
  USE mo_parallel_config,     ONLY: nproma
  USE mo_physical_constants,  ONLY: grav, rd, rdv, o_m_rdv, dtdz_standardatm, p0sl_bg
  USE mo_grid_config,         ONLY: n_dom, grid_sphere_radius
  USE mo_run_config,          ONLY: num_lev, ntracer
  USE mo_io_config,           ONLY: itype_pres_msl
  USE mo_impl_constants,      ONLY: PRES_MSL_METHOD_GME, PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_DWD, &
    &                               PRES_MSL_METHOD_IFS_CORR, MODE_ICONVREMAP, MODE_IFSANA,        &
    &                               MODE_COMBINED, SUCCESS
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_initicon_config,     ONLY: zpbl1, zpbl2, l_coarse2fine_mode, init_mode, lread_vn, lvert_remap_fg
  USE mo_initicon_types,      ONLY: t_init_state, t_initicon_state
  USE mo_fortran_tools,       ONLY: init, set_acc_host_or_device, assert_acc_device_only, &
    &                               assert_acc_host_only, minval_1d, minval_2d
  USE mo_vertical_coord_table,ONLY: vct_a
  USE mo_nh_init_utils,       ONLY: interp_uv_2_vn, adjust_w, convert_thdvars
  USE mo_util_phys,           ONLY: virtual_temp, vap_pres
  USE mo_math_gradients,      ONLY: grad_fd_norm, grad_fd_tang
  USE mo_loopindices,         ONLY: get_indices_e, get_indices_c
  USE mo_grf_intp_data_strc,  ONLY: t_gridref_state
  USE mo_grf_bdyintp,         ONLY: interpol_scal_grf, interpol2_vec_grf
  USE mo_sync,                ONLY: sync_patch_array, SYNC_C, SYNC_E
  USE mo_satad,               ONLY: sat_pres_water
  USE mo_nwp_sfc_interp,      ONLY: process_sfcfields
  USE mo_upatmo_config,       ONLY: upatmo_config
  USE mo_nh_vert_extrap_utils,ONLY: t_expol_state
  USE mo_dynamics_config,     ONLY: ldeepatmo
  USE mo_deepatmo,            ONLY: deepatmo_htrafo

  IMPLICIT NONE
  PRIVATE


  ! Threshold for switching between analytical formulas for constant temperature and
  ! constant vertical gradient of temperature, respectively
  REAL(wp), PARAMETER :: dtvdz_thresh = 1.e-4_wp ! 0.1 K/km

  ! Artificial limits on temperature profile used for extrapolation below the ground
  REAL(wp), PARAMETER :: t_low  = 255.0_wp
  REAL(wp), PARAMETER :: t_high = 290.5_wp

  ! fill value for masked (data-void) grid points
  REAL(wp), PARAMETER :: fill_temp = 288.15_wp

  ! Height above ground level from where the temperature is used for downward extrapolation
  REAL(wp), PARAMETER :: zagl_extrap_1 = 10._wp
  REAL(wp), PARAMETER :: topo_extrap_1 = 1500._wp  ! 10 m AGL valid for grid points below 1500 m
  REAL(wp), PARAMETER :: zagl_extrap_2 = 150._wp
  REAL(wp), PARAMETER :: trans_depth   = 500._wp   ! depth of transition zone 500 m


  PUBLIC :: vert_interp
  PUBLIC :: vert_interp_atm
  PUBLIC :: vert_interp_sfc
  PUBLIC :: prepare_lin_intp
  PUBLIC :: prepare_extrap, prepare_extrap_ifspp
  PUBLIC :: prepare_cubic_intp
  PUBLIC :: lin_intp, uv_intp, qv_intp, temperature_intp, pressure_intp


CONTAINS

  !-------------
  !>
  !! FUNCTION start_idx_diff_threshold
  !! Computes the minimum index where the difference between a z-level and z-level
  !! reference is less or equal the threshold.
  !! This wrapper enables the use of OpenACC
  !!
  FUNCTION start_idx_diff_threshold(threshold, z2d_in, z_reference, nlen, nlevs, lacc)
    REAL(wp), INTENT(IN) :: threshold
    ! The CONTIGUOUS keyword is necessary to circumvent an inlining-openmp bug in nfort 3.2 and 3.5.1
    REAL(wp), CONTIGUOUS, INTENT(IN) :: z2d_in(:, :), z_reference(:)
    INTEGER, INTENT(IN) :: nlen, nlevs
    LOGICAL, INTENT(IN) :: lacc ! if true, use OpenACC
    INTEGER :: start_idx_diff_threshold, jc, jk

    !nlen = SIZE(z2d_in,1)
    !nlevs = SIZE(z2d_in,2)

    start_idx_diff_threshold = nlevs-1
    ! OpenACC requires a different approach, as EXIT would be illegal
    ! and MINVAL would be very inefficient within an ACC-loop.
#ifndef _OPENACC
    DO jk = 1, nlevs
      IF (MINVAL(z2d_in(1:nlen,jk)-z_reference(1:nlen)) <= threshold) THEN
        start_idx_diff_threshold = jk - 1
        EXIT
      ENDIF
    ENDDO
#else
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) ASYNC(1) DEFAULT(PRESENT) REDUCTION(MIN: start_idx_diff_threshold) IF(lacc)
    DO jk = 1, nlevs
      DO jc = 1, nlen
        IF ( z2d_in(jc,jk)-z_reference(jc) <= threshold ) THEN
          start_idx_diff_threshold = MIN(start_idx_diff_threshold, jk)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT ! required to sync result back to CPU

    start_idx_diff_threshold = start_idx_diff_threshold - 1
    ! OpenACC does its own initialization of reduction variables with a very high value.
    ! Thus we make sure that the final value is in a valid range
    IF ( start_idx_diff_threshold > nlevs-1 ) start_idx_diff_threshold = nlevs-1
#endif
  END FUNCTION start_idx_diff_threshold


  !-------------
  !>
  !! FUNCTION start_idx_threshold
  !! Computes the minimum index where zalml is less than the threshold.
  !! This wrapper enables the use of OpenACC
  !!
  FUNCTION start_idx_threshold(threshold, zalml, nlen, nlevs, lacc)
    REAL(wp), INTENT(IN) :: threshold
    ! The CONTIGUOUS keyword is necessary to circumvent an inlining-openmp bug in nfort 3.2 and 3.5.1
    REAL(wp), CONTIGUOUS, INTENT(IN) :: zalml(:, :)
    INTEGER, INTENT(IN) :: nlen, nlevs
    LOGICAL, INTENT(IN) :: lacc ! if true, use OpenACC
    INTEGER :: start_idx_threshold, jc, jk

    !nlen = SIZE(zalml,1)
    !nlevs = SIZE(zalml,2)

    start_idx_threshold = nlevs-1
    ! OpenACC requires a different approach, as EXIT would be illegal
    ! and MINVAL would be very inefficient within an ACC-loop.
#ifndef _OPENACC
    DO jk = 1, nlevs
      IF (MINVAL(zalml(1:nlen,jk)) < threshold) THEN
        start_idx_threshold = jk
        EXIT
      ENDIF
    ENDDO
#else
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) REDUCTION(MIN: start_idx_threshold) ASYNC(1) IF(lacc)
    DO jk = 1, nlevs
      DO jc = 1, nlen
        IF ( zalml(jc,jk) < threshold ) THEN
          start_idx_threshold = MIN(start_idx_threshold, jk)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT ! required to sync result back to CPU
#endif
  END FUNCTION start_idx_threshold


  !-------------
  !>
  !! SUBROUTINE vert_interp_atm
  !! Outer driver routine for vertical interpolation of analysis
  !! data (atmosphere only) interpolated horizontally by IFS2ICON to
  !! the ICON grid
  !!
  !! Only used during initialization and no OpenACC support.
  SUBROUTINE vert_interp_atm(p_patch, p_nh_state, p_int, p_grf, initicon)

    TYPE(t_patch),          INTENT(INOUT)    :: p_patch(:)
    TYPE(t_nh_state),       INTENT(IN)       :: p_nh_state(:)
    TYPE(t_int_state),      INTENT(IN)       :: p_int(:)
    TYPE(t_gridref_state),  INTENT(IN)       :: p_grf(:)
    TYPE(t_initicon_state), INTENT(INOUT)    :: initicon(:)

    ! LOCAL VARIABLES
    INTEGER :: jg, jn, jgc

!-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE ! skip model domains not active at initial time

      WRITE(message_text,'(a,i4)') 'Vertical interpolation of analysis data, domain ',jg
      CALL message('vert_interp_atm', TRIM(message_text))

      IF (p_patch(jg)%n_patch_cells==0) CYCLE ! skip empty patches

      CALL vert_interp(p_patch(jg), p_int(jg), p_nh_state(jg)%metrics, initicon(jg), &
        &              opt_inputonzgpot=(init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED))

      ! Apply boundary interpolation for u and v because the outer nest boundary
      ! points would remain undefined otherwise
      DO jn = 1, p_patch(jg)%n_childdom

        jgc = p_patch(jg)%child_id(jn)

        CALL interpol_scal_grf (p_pp=p_patch(jg), p_pc=p_patch(jgc), p_grf=p_grf(jg)%p_dom(jn), &
                                nfields=1, lacc=.FALSE., f3din1=initicon(jg)%atm%w, f3dout1=initicon(jgc)%atm%w )

        CALL interpol2_vec_grf (p_pp=p_patch(jg), p_pc=p_patch(jgc), p_grf=p_grf(jg)%p_dom(jn), &
                                nfields=1, lacc=.FALSE., f3din1=initicon(jg)%atm%vn, f3dout1=initicon(jgc)%atm%vn )

      ENDDO
    ENDDO

  END SUBROUTINE vert_interp_atm



  !-------------
  !>
  !! SUBROUTINE vert_interp_sfc
  !! Outer driver routine for vertical interpolation of analysis
  !! data (surface only) interpolated horizontally by ICONREMAP to
  !! the ICON grid
  !!
  SUBROUTINE vert_interp_sfc(p_patch, ext_data, initicon)

    TYPE(t_patch),         INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data), INTENT(IN)    :: ext_data(:)
    CLASS(t_init_state),   INTENT(INOUT) :: initicon(:)

    ! LOCAL VARIABLES
    INTEGER :: jg

    !-------------------------------------------------------------------------

    DO jg = 1, n_dom

      IF (p_patch(jg)%n_patch_cells==0) CYCLE ! skip empty patches
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE ! skip model domains not active at initial time

      !
      ! process surface fields
      !
      CALL process_sfcfields(p_patch(jg), ext_data(jg), initicon(jg))
    ENDDO

  END SUBROUTINE vert_interp_sfc


  !-------------
  !>
  !! SUBROUTINE vert_interp
  !! Domain-wise driver routine for vertical interpolation of analysis
  !! data (atmosphere only) interpolated horizontally by ICONREMAP to
  !! the ICON grid
  !!
  !! SOURCE data taken from initicon%atm_in: 
  !! z3d (psfc or log(psfc), phi_sfc), temp, vn (u,v), qv, qc, qi, qr, qs, pres, w (omega), tke
  !! Alternatives are given in brackets.
  !! - if z3d is not available, it is derived from (psfc or log(psfc), phi_sfc)
  !! - if omega rather than w is given, it is converted to w prior to the interpolation step.
  !!
  !!
  !! NOTE: This subroutine is usually MPI-collective, since it
  !!       contains several synchronization calls. It must therefore
  !!       be passed by all worker PEs. However, there is the common
  !!       situation where vertical interpolation shall performed on a
  !!       subset of PEs only, while no valid data is available on the
  !!       remaining PE. For this situation, the optional "opt_latbcmode"
  !!       allows skipping all computations requiring MPI synchronization
  !!
  !!
  !! Fields are vertically interpolated onto the ICON vertical grid 
  !! and stored in the initicon%atm state. They are finally converted into the following 
  !! set of prognostic variables:
  !! vn, w, qv, qc, qi, qr, qs, rho, exner, theta_v
  !!
  SUBROUTINE vert_interp(p_patch, p_int, p_metrics, initicon, opt_use_vn, opt_lmask_c, opt_lmask_e, opt_latbcmode, &
    &                    opt_inputonzgpot)

    TYPE(t_patch),          INTENT(INOUT)    :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int
    TYPE(t_nh_metrics),     INTENT(IN)       :: p_metrics

    ! (Note: This is a CLASS parameter, st. we can feed
    ! t_initicon_state into this subroutine)
    CLASS(t_init_state),    INTENT(INOUT), TARGET :: initicon

    LOGICAL, OPTIONAL,      INTENT(IN)       :: opt_use_vn, opt_latbcmode, opt_inputonzgpot
    LOGICAL, OPTIONAL,      INTENT(IN)       :: opt_lmask_c(:,:)
    LOGICAL, OPTIONAL,      INTENT(IN)       :: opt_lmask_e(:,:)

    ! LOCAL VARIABLES
    CHARACTER(LEN=*), PARAMETER :: routine = 'vert_interp'

    INTEGER :: jg
    INTEGER :: nlev, nlevp1, idx
    INTEGER :: nlev_in         ! number of vertical levels in source vgrid 
    LOGICAL :: lc2f, l_use_vn, latbcmode, linputonzgpot

    ! Auxiliary fields for input data
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_c) :: temp_v_in

    ! Auxiliary field for output data
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      z_tempv

    ! CELLS !
    !
    ! Auxiliary fields for coefficients and filled input height fields
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      wfac_lin, coef1, coef2, coef3

    REAL(wp), DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: &
      wfac_lin_w

    INTEGER , DIMENSION(nproma,p_patch%nlev,p_patch%nblks_c) :: &
      idx0_lin, idx0_cub

    INTEGER , DIMENSION(nproma,p_patch%nlevp1,p_patch%nblks_c) :: &
      idx0_lin_w

    ! extrapolation weights and indices for cells
    REAL(wp), DIMENSION(nproma,p_patch%nblks_c) :: &
      wfacpbl1, wfacpbl2, slope

    INTEGER , DIMENSION(nproma,p_patch%nblks_c) :: &
      bot_idx_lin, bot_idx_lin_w, bot_idx_cub, kpbl1, kpbl2

    ! array for filling input height field on masked points
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_c), TARGET :: &
      z_mc_fill

    ! EDGES !
    !
    ! Auxiliary fields for coefficients at edges
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: &
      wfac_lin_e, coef1_e, coef2_e, coef3_e

    INTEGER , DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: &
      idx0_lin_e, idx0_cub_e

    ! extrapolation weights
    REAL(wp), DIMENSION(nproma,p_patch%nblks_e) :: &
      wfacpbl1_e, wfacpbl2_e

    INTEGER , DIMENSION(nproma,p_patch%nblks_e) :: &
      bot_idx_lin_e, bot_idx_cub_e, kpbl1_e, kpbl2_e

    ! full level heights of initicon input fields at edge midpoints
    REAL(wp), DIMENSION(nproma,initicon%atm_in%nlev,p_patch%nblks_e) :: atm_in_z_me

    ! full level heights of ICON vertical grid at edge midpoints
    REAL(wp), DIMENSION(nproma,p_patch%nlev,p_patch%nblks_e) :: z_me

    REAL(wp), DIMENSION(:,:,:), POINTER :: z_mc_in

    ! upper-atmospher extrapolation type
    TYPE(t_expol_state) :: expol
    LOGICAL :: lexpol

    INTEGER  :: istat

!-------------------------------------------------------------------------

    ! Does not support OpenACC
    ! In order to support OpenACC, the lacc flags in the following code must be turned on.
    ! Additional porting might be necessary.

    jg = p_patch%id

    nlev_in = initicon%atm_in%nlev

    IF (nlev_in == 0) THEN
       CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
    END IF

    IF (PRESENT(opt_use_vn)) THEN
      l_use_vn = opt_use_vn
    ELSE
      l_use_vn = lread_vn ! use vn field if available as input from file
    ENDIF

    IF (PRESENT(opt_latbcmode)) THEN
      latbcmode = opt_latbcmode  ! true: mode for interpolating lateral boundary conditions
    ELSE
      latbcmode = .FALSE.
    ENDIF

    ! The switch opt_inputonzgpot is only relevant for ldeepatmo = .TRUE. (see below)
    IF (PRESENT(opt_inputonzgpot)) THEN
      linputonzgpot = opt_inputonzgpot  ! true: input data defined on geopotential height levels
    ELSE
      linputonzgpot = .FALSE.
    ENDIF

    ! (Deep atmosphere: 
    ! The 'z_mc_in'-field is the geopotential height in case of IFS-date. For the shallow atmosphere geopotential 
    ! and geometric height coincide. Unfortunately this does no longer hold for the deep atmosphere. 
    ! So a further step is required in which geopotential height is transformed into geometric height.)
    IF (ldeepatmo .AND. linputonzgpot) THEN
      CALL deepatmo_htrafo(z_inout=initicon%const%z_mc_in,                                             &
        & nblks_nproma_npromz=[p_patch%nblks_c,nproma,p_patch%npromz_c], start_end_levels=[1,nlev_in], &
        & radius=grid_sphere_radius, trafo_type='zgpot2z', ierror=istat, lacc=.FALSE.)
      IF (istat /= SUCCESS) CALL finish(routine, "deepatmo_htrafo failed.")
    ENDIF

    IF (PRESENT(opt_lmask_c)) THEN
      IF (.NOT. ANY(opt_lmask_c)) RETURN
      CALL fill_input_height(initicon%const%z_mc_in, z_mc_fill, p_patch%nblks_c, p_patch%npromz_c, nlev_in, opt_lmask_c)
      z_mc_in  => z_mc_fill
    ELSE
      z_mc_in  => initicon%const%z_mc_in
    ENDIF

    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ! Switch to determine if corrections for coarse-to-fine-mesh interpolation are to be used
    lc2f   = l_coarse2fine_mode(p_patch%id) .AND. .NOT. latbcmode

    ! Compute virtual temperature of input data
    CALL virtual_temp(p_patch, initicon%atm_in%temp, initicon%atm_in%qv, initicon%atm_in%qc, &
                      initicon%atm_in%qi, initicon%atm_in%qr, initicon%atm_in%qs,            &
                      temp_v=temp_v_in)

    ! Prepare interpolation coefficients for cells
    IF (.NOT. ASSOCIATED(initicon%const%z_mc)) THEN
       CALL finish(routine, "Internal error!")
    END IF
    IF (.NOT. ASSOCIATED(initicon%const%z_mc_in)) THEN
       CALL finish(routine, "Internal error: z_mc_in not associated!")
    END IF
    ! (Deep atmosphere: it is assumed here and in the following 
    ! that the weighting factors computed in 'prepare_...(_intp)' result from 
    ! a pure distance weighting (along coordinate lines etc.), 
    ! so no deep-atmosphere modification is applied, this would have to be reconsidered, 
    ! if the weighting would be regarded as a volume weighting)
    CALL prepare_lin_intp(z_mc_in, initicon%const%z_mc,                     &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                          wfac_lin, idx0_lin, bot_idx_lin, lacc=.FALSE.)

    CALL prepare_extrap(z_mc_in,                                    &
                        p_patch%nblks_c, p_patch%npromz_c, nlev_in, &
                        kpbl1, wfacpbl1, kpbl2, wfacpbl2, lacc=.FALSE.)



    CALL prepare_cubic_intp(z_mc_in, initicon%const%z_mc,                     &
                            p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                            coef1, coef2, coef3, idx0_cub, bot_idx_cub, lacc=.FALSE.)

    ! (Initialize upper-atmosphere extrapolation type.)
    lexpol = upatmo_config(jg)%exp%l_expol
    IF (lexpol) CALL expol%initialize(p_patch, latbcmode)


    ! Perform vertical interpolation

    ! Temperature
    IF (lc2f) CALL compute_slope(p_patch, p_int, initicon%const%topography_c, slope)

    CALL temperature_intp(initicon%atm_in%temp, initicon%atm%temp,                &
                          z_mc_in, initicon%const%z_mc,                           &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,       &
                          coef1, coef2, coef3, wfac_lin,                          &
                          idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,           &
                          wfacpbl1, kpbl1, wfacpbl2, kpbl2,                       &
                          l_restore_sfcinv=.TRUE., l_hires_corr=lc2f,             &
                          extrapol_dist=-1500._wp, l_pz_mode=.FALSE., slope=slope,&
                          opt_lmask=opt_lmask_c, lacc=.FALSE. )

    ! (Extrapolate temperature to upper atmosphere
    ! Note: the following subroutine is a post-processing routine, 
    ! it has to be positioned after 'temperature_intp'!
    ! This applies to the remaining extrapolations as well.)                     
    IF (lexpol) CALL expol%temp(p_patch, initicon%atm%temp)


    ! horizontal wind components
    IF (l_use_vn) THEN

      !
      ! prepare interpolation coefficients for edges
      !

      ! full level heights at edges for initicon input
      CALL cells2edges_scalar(z_mc_in, p_patch, p_int%c_lin_e, atm_in_z_me, opt_fill_latbc=.TRUE., lacc=.FALSE.)


      ! full level heights at edges for ICON vertical grid
      CALL cells2edges_scalar(p_metrics%z_mc, p_patch, p_int%c_lin_e, z_me, opt_fill_latbc=.TRUE., lacc=.FALSE.)


      ! compute extrapolation coefficients for edges
      CALL prepare_extrap(atm_in_z_me, p_patch%nblks_e, p_patch%npromz_e, &
        &                 nlev_in, kpbl1_e, wfacpbl1_e, kpbl2_e, wfacpbl2_e, lacc=.FALSE.)


      ! compute linear interpolation coefficients for edges
      CALL prepare_lin_intp(atm_in_z_me, z_me,                                &
                            p_patch%nblks_e, p_patch%npromz_e, nlev_in, nlev, &
                            wfac_lin_e, idx0_lin_e, bot_idx_lin_e, lacc=.FALSE.)

      ! compute cubic interpolation coefficients for edges
      CALL prepare_cubic_intp(atm_in_z_me, z_me,                                &
                            p_patch%nblks_e, p_patch%npromz_e, nlev_in, nlev,   &
                            coef1_e, coef2_e, coef3_e, idx0_cub_e, bot_idx_cub_e, lacc=.FALSE.)

      ! vertically interpolate vn to ICON grid
      CALL uv_intp(initicon%atm_in%vn, initicon%atm%vn,                  &
                   atm_in_z_me, z_me,                                    &
                   p_patch%nblks_e, p_patch%npromz_e, nlev_in, nlev,     &
                   coef1_e, coef2_e, coef3_e, wfac_lin_e,                &
                   idx0_cub_e, idx0_lin_e, bot_idx_cub_e, bot_idx_lin_e, &
                   wfacpbl1_e, kpbl1_e, wfacpbl2_e, kpbl2_e,             &
                   l_hires_intp=lc2f, lacc=.FALSE.                       )

    ELSE
      CALL uv_intp(initicon%atm_in%u, initicon%atm%u,                &
                   z_mc_in, initicon%const%z_mc,                     &
                   p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                   coef1, coef2, coef3, wfac_lin,                    &
                   idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                   wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                   l_hires_intp=lc2f, lacc=.FALSE.                   )
      CALL uv_intp(initicon%atm_in%v, initicon%atm%v,                &
                   z_mc_in, initicon%const%z_mc,                     &
                   p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev, &
                   coef1, coef2, coef3, wfac_lin,                    &
                   idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,     &
                   wfacpbl1, kpbl1, wfacpbl2, kpbl2,                 &
                   l_hires_intp=lc2f, lacc=.FALSE.                   )

      ! Convert u and v on cell points to vn at edge points
      CALL interp_uv_2_vn(p_patch, p_int, initicon%atm%u, initicon%atm%v, initicon%atm%vn)
    ENDIF

    ! This synchronization is executed in the calling routine in latbc mode
    ! (Postpone 'sync' to after extrapolation below in case 
    ! upper-atmosphere extrapolation has been switched on)
    IF (.NOT. latbcmode .AND. (.NOT. lexpol)) CALL sync_patch_array(SYNC_E,p_patch,initicon%atm%vn)


    ! Preliminary interpolation of QV: this is needed to compute virtual temperature
    ! below, which in turn is required to integrate the hydrostatic equation
    ! A lower limit of 2.5 ppm is imposed
    CALL lin_intp(initicon%atm_in%qv, initicon%atm%qv,                 &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,    &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,    &
                  wfacpbl2, kpbl2, l_loglin=.TRUE., l_extrapol=.TRUE., &
                  l_pd_limit=.TRUE., lower_limit=2.5e-7_wp, lacc=.FALSE.)

    ! Cloud and precipitation variables - linear interpolation only because cubic may
    ! cause negative values, and no-gradient condition for downward extrapolation
    ! Positive definite limiter is used in all cases
    CALL lin_intp(initicon%atm_in%qc, initicon%atm%qc,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE., lacc=.FALSE.)

    CALL lin_intp(initicon%atm_in%qi, initicon%atm%qi,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE., lacc=.FALSE.)

    CALL lin_intp(initicon%atm_in%qr, initicon%atm%qr,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE., lacc=.FALSE.)

    CALL lin_intp(initicon%atm_in%qs, initicon%atm%qs,                   &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,      &
                  wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.TRUE., lacc=.FALSE.)

    ! ... additional tracer variables
    DO idx=1, ntracer
      IF ( ASSOCIATED(initicon%atm_in%tracer(idx)%field) ) THEN
        IF ( .NOT. latbcmode ) THEN
          ! allocate target array for vertical interpolation
          ALLOCATE(initicon%atm%tracer(idx)%field(nproma,nlev,p_patch%nblks_c))
!$OMP PARALLEL
          CALL init(initicon%atm%tracer(idx)%field(:,:,:), lacc=.FALSE.)   !_jf: necessary?
!$OMP END PARALLEL
          ! set pointer to var_element of atm_in
          initicon%atm%tracer(idx)%var_element => initicon%atm_in%tracer(idx)%var_element
        ENDIF
        ! do vertical interpolation
        CALL lin_intp(initicon%atm_in%tracer(idx)%field, initicon%atm%tracer(idx)%field, &
                      p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,                  &
                      wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,                  &
                      wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE.,             &
                      l_pd_limit=.TRUE., lacc=.FALSE.)
      ENDIF
    ENDDO

    ! (Extrapolate water phases, except for water vapour 'qv', to upper atmosphere.
    ! Note: this is only done for subroutine-internal use, 
    ! see the final treatment of the 'qx' in 'mo_initicon_utils: copy_initicon2prog_atm')
    IF (lexpol) CALL expol%qx(p_patch, initicon%atm%qc, initicon%atm%qi, initicon%atm%qr, initicon%atm%qs)

    ! Compute virtual temperature with preliminary QV
    CALL virtual_temp(p_patch, initicon%atm%temp, initicon%atm%qv, initicon%atm%qc,    &
                      initicon%atm%qi, initicon%atm%qr, initicon%atm%qs, temp_v=z_tempv)

    ! Interpolate pressure on ICON grid
    CALL pressure_intp_initmode(initicon%atm_in%pres, temp_v_in, z_mc_in,              &
      &                        initicon%atm%pres, z_tempv, initicon%const%z_mc,        &
      &                        p_patch%nblks_c, p_patch%npromz_c, nlev, nlev_in,       &
      &                        wfac_lin, idx0_lin, bot_idx_lin, opt_lmask=opt_lmask_c)

    ! (Extrapolate pressure to upper atmosphere) 
    IF (lexpol) CALL expol%pres(p_patch, initicon%atm%pres, z_tempv, p_metrics)


    CALL qv_intp(initicon%atm_in%qv, initicon%atm%qv, z_mc_in,                      &
                 initicon%const%z_mc, initicon%atm_in%temp, initicon%atm_in%pres,   &
                 initicon%atm%temp, initicon%atm%pres,                              &
                 p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlev,                  &
                 coef1, coef2, coef3, wfac_lin,                                     &
                 idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,                      &
                 wfacpbl1, kpbl1, wfacpbl2, kpbl2, l_satlimit=.TRUE.,               &
                 lower_limit=2.5e-7_wp, l_restore_pbldev=.TRUE.,                    &
                 opt_hires_corr=lc2f, opt_qc=initicon%atm%qc, opt_lmask=opt_lmask_c,&
                 lacc=.FALSE. )

    ! Compute virtual temperature with final QV
    CALL virtual_temp(p_patch, initicon%atm%temp, initicon%atm%qv, initicon%atm%qc,    &
                      initicon%atm%qi, initicon%atm%qr, initicon%atm%qs, temp_v=z_tempv)

    ! Final interpolation of pressure on ICON grid
    CALL pressure_intp_initmode(initicon%atm_in%pres, temp_v_in, z_mc_in,              &
      &                        initicon%atm%pres, z_tempv, initicon%const%z_mc,        &
      &                        p_patch%nblks_c, p_patch%npromz_c, nlev, nlev_in,       &
      &                        wfac_lin, idx0_lin, bot_idx_lin, opt_lmask=opt_lmask_c) 

    ! (In case of an upper-atmosphere extrapolation, 'hydro_adjust' should not be called 
    ! in 'src/atm_dyn_iconam/mo_initicon_utils: copy_initicon2prog_atm', so it has to be done here. 
    ! In its current form, this can be covered by 'expol%pres', 
    ! but be careful, if you change something here or there.) 
    IF (lexpol) THEN
      CALL expol%pres(p_patch, initicon%atm%pres, z_tempv, p_metrics, opt_slev=1, opt_elev=nlev)
      ! (Below, the extrapolation of vn requires the horizontal pressure gradient, 
      ! but 'initicon%atm%pres' should have been processed on the entire grid, so no sync should be necessary) 
    ENDIF


    ! Convert thermodynamic variables into set of NH prognostic variables
    CALL convert_thdvars(p_patch, initicon%atm%pres, z_tempv,    &
                         initicon%atm%rho, initicon%atm%exner,   &
                         initicon%atm%theta_v)

    ! (The extrapolation of 'vn' makes use of the geostrophic balance, 
    ! so it should be computed after the thermodynamic state has been settled.)
    IF (lexpol) THEN
      CALL expol%vn(p_patch, initicon%atm%vn, initicon%atm%theta_v, initicon%atm%exner, p_metrics, p_int)
      IF (.NOT. latbcmode) CALL sync_patch_array(SYNC_E,p_patch,initicon%atm%vn)
    ENDIF    


    ! Compute coefficients for w interpolation
    CALL prepare_lin_intp(z_mc_in, initicon%const%z_ifc,       &
                          p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1, &
                          wfac_lin_w, idx0_lin_w, bot_idx_lin_w, lacc=.FALSE.)


    ! Perform linear interpolation of w
    ! Note: the coefficients for gradient computation (*pbl*) do not have to be recomputed
    ! because of l_extrapol=.FALSE.,
    CALL lin_intp(initicon%atm_in%w, initicon%atm%w,                     &
                  p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1,    &
                  wfac_lin_w, idx0_lin_w, bot_idx_lin_w, wfacpbl1, kpbl1,&
                  wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                  l_pd_limit=.FALSE., lacc=.FALSE.)

    ! (Extrapolation of the vertical wind to the upper atmosphere.)
    IF (lexpol) CALL expol%w(p_patch, initicon%atm%w)

    IF (.NOT. latbcmode) THEN

      ! Impose appropriate lower boundary condition on vertical wind field
      CALL adjust_w(p_patch, p_int, initicon%atm%vn, initicon%const%z_ifc, initicon%atm%w)

      CALL sync_patch_array(SYNC_C,p_patch,initicon%atm%w)
    ENDIF

    IF ((init_mode == MODE_ICONVREMAP .OR. lvert_remap_fg) .AND. ASSOCIATED(initicon%atm_in%tke)) THEN
      CALL lin_intp(initicon%atm_in%tke, initicon%atm%tke,                 &
                    p_patch%nblks_c, p_patch%npromz_c, nlev_in, nlevp1,    &
                    wfac_lin_w, idx0_lin_w, bot_idx_lin_w, wfacpbl1, kpbl1,&
                    wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.FALSE., &
                    l_pd_limit=.TRUE., lower_limit=1.e-5_wp, lacc=.FALSE.)
    ENDIF

    ! (Finalize upper-atmosphere extrapolation)
    IF (lexpol) CALL expol%finalize(p_patch)

  END SUBROUTINE vert_interp

  !-------------
  !>
  !! SUBROUTINE fill_input_height
  !!
  !! Fills input height field with dummy values if a mask field for data-void grid points is present
  !!
  SUBROUTINE fill_input_height(z3d_in, z3d_fill, nblks, npromz, nlev, lmask)

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlev       ! Number of input levels

    LOGICAL, INTENT(IN)  :: lmask(:,:) ! mask field, false on data-void points

    ! Output fields
    REAL(wp), INTENT(OUT) :: z3d_fill(:,:,:) ! height coordinate field filled with dummy values on data-void grid points

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jc1, jb1
    INTEGER :: nlen

!-------------------------------------------------------------------------

    ! Detect first grid point for which the mask field is .true.
    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
      ENDIF
      DO jc = 1, nlen
        IF (lmask(jc,jb)) THEN
          jc1 = jc
          jb1 = jb
          EXIT
        ENDIF
      ENDDO
    ENDDO

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen
          IF (lmask(jc,jb)) THEN
            z3d_fill(jc,jk,jb) = z3d_in(jc,jk,jb)
          ELSE
            z3d_fill(jc,jk,jb) = z3d_in(jc1,jk,jb1)
          ENDIF
        ENDDO
      ENDDO

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE fill_input_height


  !-------------
  !>
  !! SUBROUTINE prepare_lin_intp
  !! Computes coefficient fields for linear vertical interpolation
  !!
  !! Required input fields: 3D coordinate fields of input and output data
  !! Output: weighting coefficient, index of upper level,
  !! index of lowest level for which interpolation is possible
  !! (as opposed to extrapolation)
  !!
  !! It is assumed that the highest level of the input data is at least
  !! as high as the highest level of the output data
  !!
  SUBROUTINE prepare_lin_intp(z3d_in, z3d_out,                    &
                              nblks, npromz, nlevs_in, nlevs_out, &
                              wfac, idx0, bot_idx, lextrap, lacc )

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Switch for extrapolation beyond top of input data
    LOGICAL, INTENT(IN), OPTIONAL :: lextrap ! if true, apply linear extrapolation

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true use OpenACC

    ! Output fields
    REAL(wp), INTENT(OUT) :: wfac(:,:,:)       ! weighting factor of upper level
    INTEGER , INTENT(OUT) :: idx0(:,:,:)       ! index of upper level
    INTEGER , INTENT(OUT) :: bot_idx(:,:)      ! index of lowest level for which interpolation is possible

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk1, jk_start
    INTEGER :: nlen, ierror(nblks), nerror
    LOGICAL :: l_found(nproma), lfound_all, l_extrap, lzacc

!-------------------------------------------------------------------------

    IF (PRESENT(lextrap)) THEN
      l_extrap = lextrap
    ELSE
      l_extrap = .TRUE.
    ENDIF

    CALL set_acc_host_or_device(lzacc, lacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk1,jk_start,l_found,lfound_all) ICON_OMP_DEFAULT_SCHEDULE
    !$ACC DATA CREATE(l_found) IF(lzacc)

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = nlen+1, nproma
          bot_idx(jc,jb) = nlevs_out
          wfac(jc,:,jb)  = 0.5_wp
          idx0(jc,:,jb)  = nlevs_in-1
        ENDDO
        !$ACC END PARALLEL
      ENDIF
      ierror(jb) = 0

      jk_start = 1
      DO jk = 1, nlevs_out
        lfound_all = .FALSE.
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          l_found(jc) = .FALSE.
        ENDDO

        !$ACC LOOP SEQ
        DO jk1 = jk_start,nlevs_in-1
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            IF (.NOT. l_found(jc) ) THEN
              IF (z3d_out(jc,jk,jb) <= z3d_in(jc,jk1,jb) .AND. &
                  z3d_out(jc,jk,jb) >  z3d_in(jc,jk1+1,jb)) THEN
                idx0(jc,jk,jb) = jk1
                wfac(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb))/&
                                 (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
                bot_idx(jc,jb) = jk
                l_found(jc) = .TRUE.
              ELSE IF (z3d_out(jc,jk,jb) <= z3d_in(jc,nlevs_in,jb)) THEN
                l_found(jc) = .TRUE.
                idx0(jc,jk,jb) = nlevs_in
                IF (jk == 1)  bot_idx(jc,jb) = 0
              ELSE IF (z3d_out(jc,jk,jb) > z3d_in(jc,1,jb)) THEN ! linear extrapolation
                idx0(jc,jk,jb) = 1
                IF (l_extrap) THEN
                  wfac(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,2,jb))/&
                                   (z3d_in(jc,1,jb)-z3d_in(jc,2,jb))
                ELSE
                  wfac(jc,jk,jb) = 1._wp ! use constant values above top of input data
                ENDIF
                bot_idx(jc,jb) = jk
                l_found(jc) = .TRUE.
              ENDIF
            ENDIF
          ENDDO

#ifndef _OPENACC
          ! ACC: the following EXIT would be illegal within an OpenACC Kernel, thus we skip this CPU optimization
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
#endif

        ENDDO
        !$ACC END PARALLEL

#ifdef _OPENACC
        lfound_all = .TRUE.
        ! The following reduction must appear in its own small kernel as it did not work otherwise with Nvidia 21.2
        !$ACC PARALLEL LOOP VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc) REDUCTION(.AND.: lfound_all)
        DO jc = 1, nlen
          lfound_all = lfound_all .AND. l_found(jc)
        ENDDO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT
#endif

        IF (lfound_all) THEN
          jk_start = MIN(minval_1d(idx0(1:nlen,jk,jb), lacc=lzacc),nlevs_in-1)
        ELSE
          ierror(jb) = ierror(jb) + 1
          ! Write extra debug output before finishing
          WRITE(0,*) 'prepare_lin_intp',jb,jk,nlevs_out,jk_start,jk1,nlevs_in
          !$ACC UPDATE HOST(l_found, z3d_in, z3d_in, z3d_out) ASYNC(1) IF(lzacc)
          !$ACC WAIT(1)
          DO jc = 1, nlen
            IF(.NOT.l_found(jc)) THEN
              WRITE(0,*)'prepare_lin_intp',z3d_in(jc,jk_start:jk1,jb),&
                z3d_in(jc,nlevs_in,jb),z3d_out(jc,jk,jb)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      jk_start = minval_1d(bot_idx(1:nlen,jb), lacc=lzacc) + 1

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = jk_start, nlevs_out
        DO jc = 1, nlen
          IF (jk >= bot_idx(jc,jb)+1) THEN
            ! Store extrapolation distance on wfac if target point is below the
            ! surface of the input data (note: this is a negative quantity)
            idx0(jc,jk,jb) = nlevs_in
            wfac(jc,jk,jb) = z3d_out(jc,jk,jb)-z3d_in(jc,nlevs_in,jb)
          ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL


    ENDDO
    !$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    nerror = SUM(ierror)

    IF (nerror > 0) CALL finish("prepare_lin_intp:",&
      "Error in computing linear interpolation coefficients")

  END SUBROUTINE prepare_lin_intp


  !-------------
  !>
  !! SUBROUTINE prepare_extrap
  !! Computes coefficient fields for vertical extrapolation below the surface level
  !! of the input model
  !!
  !! Required input fields: 3D coordinate fields of input data
  !! Output: index and coefficient fields to compute field values at
  !! 500 m and 1000 m above ground
  !!
  SUBROUTINE prepare_extrap(z3d_in, nblks, npromz, nlevs_in, &
                            kpbl1, wfacpbl1, kpbl2, wfacpbl2, lacc )

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Output fields
    INTEGER , INTENT(OUT) :: kpbl1(:,:), & ! Indices of model levels lying immediately above
                             kpbl2(:,:)    ! (by default) 500 m / 1000 m AGL
    REAL(wp), INTENT(OUT) :: wfacpbl1(:,:), & ! Corresponding interpolation coefficients
                             wfacpbl2(:,:)

    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true use OpenACC
    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk_start
    INTEGER :: nlen
    REAL(wp):: kpbl1_min, kpbl2_min
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk_start) ICON_OMP_DEFAULT_SCHEDULE

    !$ACC DATA PRESENT(z3d_in, kpbl1, kpbl2, wfacpbl1, wfacpbl2) IF(lzacc)

    DO jb = 1, nblks

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = 1, nproma
        kpbl1(jc,jb) = -1
        kpbl2(jc,jb) = -1
      ENDDO
      !$ACC END PARALLEL LOOP

      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = nlen+1, nproma
          kpbl1(jc,jb) = nlevs_in
          kpbl2(jc,jb) = nlevs_in

          wfacpbl1(jc,jb) = 0.5_wp
          wfacpbl2(jc,jb) = 0.5_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      jk_start = start_idx_diff_threshold(zpbl2, z3d_in(:, : ,jb), z3d_in(:,nlevs_in,jb), nlen, nlevs_in, lzacc)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = jk_start, nlevs_in-1
        DO jc = 1, nlen

          IF (z3d_in(jc,jk,jb)  >= z3d_in(jc,nlevs_in,jb)+zpbl2 .AND. &
              z3d_in(jc,jk+1,jb) < z3d_in(jc,nlevs_in,jb)+zpbl2) THEN
            kpbl2(jc,jb) = jk
            wfacpbl2(jc,jb) = (z3d_in(jc,nlevs_in,jb)+zpbl2 - z3d_in(jc,jk+1,jb)) / &
                              (z3d_in(jc,jk,jb)            - z3d_in(jc,jk+1,jb))
          ENDIF

          IF (z3d_in(jc,jk,jb)  >= z3d_in(jc,nlevs_in,jb)+zpbl1 .AND. &
              z3d_in(jc,jk+1,jb) < z3d_in(jc,nlevs_in,jb)+zpbl1) THEN
            kpbl1(jc,jb) = jk
            wfacpbl1(jc,jb) = (z3d_in(jc,nlevs_in,jb)+zpbl1 - z3d_in(jc,jk+1,jb)) / &
                              (z3d_in(jc,jk,jb)            - z3d_in(jc,jk+1,jb))
          ENDIF

        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    kpbl1_min = minval_2d(kpbl1(:,1:nblks), lzacc)
    kpbl2_min = minval_2d(kpbl2(:,1:nblks), lzacc)

    !$ACC WAIT
    !$ACC END DATA

    IF ( kpbl1_min < 0 .OR. kpbl2_min < 0 ) THEN
      CALL finish("prepare_extrap:", &
        &         "No kpbl found, check vertical coordinate input data.")
    ENDIF

  END SUBROUTINE prepare_extrap


  !-------------
  !>
  !! SUBROUTINE prepare_extrap_ifspp
  !! Computes coefficient fields for vertical extrapolation below the surface level
  !! of the input model for IFS postprocessing method
  !!
  !! Required input fields: 3D coordinate fields of input data
  !! Output: index and coefficient fields to compute field values at
  !! 150 m above ground
  !!
  SUBROUTINE prepare_extrap_ifspp(z3d_h_in, z3d_in, nblks, npromz, nlevs_in, kextrap, zextrap, wfac_extrap, lacc)

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_h_in(:,:,:) ! half-level height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! full-level height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Output fields
    INTEGER , INTENT(OUT) :: kextrap(:,:) ! Indices of model levels lying immediately above zextrap
    REAL(wp), INTENT(OUT) :: zextrap(:,:) ! AGL height from which downward extrapolation starts (between 10 m and 150 m)
    REAL(wp), INTENT(OUT) :: wfac_extrap(:,:) ! Corresponding interpolation coefficients
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk_start
    INTEGER :: nlen
    REAL(wp) :: zagl_extrap


!-------------------------------------------------------------------------

    CALL assert_acc_device_only("prepare_extrap_ifspp", lacc)

    ! Use extrapolation from 150 m AGL if one of the IFS methods is selected, and 
    ! orography-height dependent blending between 10 m and 150 m AGL if the new DWD
    ! method is selected (which constitutes a mixture between the IFS method and the old GME method)
    SELECT CASE (itype_pres_msl)
    CASE (PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_IFS_CORR)
      zagl_extrap = zagl_extrap_2
    CASE (PRES_MSL_METHOD_DWD)
      zagl_extrap = zagl_extrap_1
    END SELECT

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk_start) ICON_OMP_DEFAULT_SCHEDULE

    !$ACC DATA PRESENT(z3d_h_in, z3d_in, kextrap, zextrap, wfac_extrap)

    DO jb = 1, nblks
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
      DO jc = 1, nproma
        kextrap(jc,jb) = -1
      ENDDO
      !$ACC END PARALLEL LOOP
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
        DO jc = nlen+1, nproma
          kextrap(jc,jb) = nlevs_in
          wfac_extrap(jc,jb) = 0.5_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      ! Compute start height above ground for downward extrapolation, depending on topography height
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = 1, nlen
        IF (z3d_h_in(jc,nlevs_in+1,jb) <= topo_extrap_1) THEN
          zextrap(jc,jb) = zagl_extrap
        ELSE IF (z3d_h_in(jc,nlevs_in+1,jb) >= topo_extrap_1 + trans_depth) THEN
          zextrap(jc,jb) = zagl_extrap_2
        ELSE
          zextrap(jc,jb) = zagl_extrap + (zagl_extrap_2 - zagl_extrap) *            &
                           (z3d_h_in(jc,nlevs_in+1,jb) - topo_extrap_1) / trans_depth
        ENDIF
      ENDDO
      !$ACC END PARALLEL

      jk_start = start_idx_diff_threshold(zagl_extrap_2, z3d_in(:,:,jb), z3d_h_in(:,nlevs_in+1,jb), nlen, nlevs_in, lacc=.TRUE.)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ ! These two loops must not be collapsed because of kextrap
      DO jk = jk_start, nlevs_in-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen

          IF (z3d_in(jc,jk,jb)  >= z3d_h_in(jc,nlevs_in+1,jb)+zextrap(jc,jb) .AND. &
              z3d_in(jc,jk+1,jb) < z3d_h_in(jc,nlevs_in+1,jb)+zextrap(jc,jb) .OR.  &
              kextrap(jc,jb) == -1 .AND. jk == nlevs_in-1) THEN 
            kextrap(jc,jb) = jk
            wfac_extrap(jc,jb) = (z3d_h_in(jc,nlevs_in+1,jb)+zextrap(jc,jb) - z3d_in(jc,jk+1,jb)) / &
                                 (z3d_in(jc,jk,jb)                          - z3d_in(jc,jk+1,jb))
          ENDIF

        ENDDO

      ENDDO
      !$ACC END PARALLEL

    ENDDO

    !$ACC END DATA
    !$ACC WAIT

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prepare_extrap_ifspp

  !>
  !! SUBROUTINE prepare_cubic_intp
  !! Computes coefficient fields for cubic vertical interpolation
  !!
  !! Required input fields: 3D coordinate fields of input and output data
  !! Output: weighting coefficients, index of upper level,
  !! index of lowest level for which cubic interpolation is possible
  !!
  !! It is assumed that the highest level of the input data is at least
  !! as high as the highest level of the output data
  !!
  SUBROUTINE prepare_cubic_intp(z3d_in, z3d_out,                    &
                                nblks, npromz, nlevs_in, nlevs_out, &
                                coef1, coef2, coef3, idx0, bot_idx, lacc)

    ! Input fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:) ! height coordinate field of input data (m)
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)! height coordinate field of input data (m)

    ! Input dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Output fields
    REAL(wp), INTENT(OUT) :: coef1(:,:,:)      ! coefficient for linear term
    REAL(wp), INTENT(OUT) :: coef2(:,:,:)      ! coefficient for quadratic term
    REAL(wp), INTENT(OUT) :: coef3(:,:,:)      ! coefficient for cubic term
    INTEGER , INTENT(OUT) :: idx0(:,:,:)       ! index of upper level
    INTEGER , INTENT(OUT) :: bot_idx(:,:)      ! index of lowest level for which interpolation is possible

    LOGICAL , INTENT(IN), OPTIONAL :: lacc      ! If true, use OpenACC

    ! LOCAL VARIABLES

    INTEGER :: jb, jk, jc, jk1, jk_start
    INTEGER :: nlen, ierror(nblks), nerror
    LOGICAL :: l_found(nproma),lfound_all
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (nlevs_in <= 3) CALL finish("prepare_cubic_intp:",&
      "Error, number of levels too small for cubic interpolation")

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,jk1,jk_start,l_found,lfound_all) ICON_OMP_DEFAULT_SCHEDULE

    !$ACC DATA CREATE(l_found) IF(lzacc)
    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = nlen+1, nproma
          bot_idx(jc,jb) = nlevs_out
          coef1(jc,:,jb) = 0.5_wp
          coef2(jc,:,jb) = 0.5_wp
          coef3(jc,:,jb) = 0.5_wp
          idx0(jc,:,jb)  = nlevs_in-1
        ENDDO
        !$ACC END PARALLEL
      ENDIF
      ierror(jb) = 0

      jk_start = 2
      DO jk = 1, nlevs_out
        lfound_all = .FALSE.
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          l_found(jc) = .FALSE.
        ENDDO

        !$ACC LOOP SEQ
        DO jk1 = jk_start,nlevs_in-2 ! cubic interpolation requires jk1-1, jk1, jk1+1 and jk1+2
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            IF (.NOT. l_found(jc) ) THEN
              IF (z3d_out(jc,jk,jb) <= z3d_in(jc,jk1,jb) .AND. &
                  z3d_out(jc,jk,jb) >  z3d_in(jc,jk1+1,jb)) THEN
                idx0(jc,jk,jb)  = jk1-1
                coef1(jc,jk,jb) = z3d_out(jc,jk,jb)-z3d_in(jc,jk1-1,jb)
                coef2(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1-1,jb))*&
                                  (z3d_out(jc,jk,jb)-z3d_in(jc,jk1  ,jb))
                coef3(jc,jk,jb) = (z3d_out(jc,jk,jb)-z3d_in(jc,jk1-1,jb))*&
                                  (z3d_out(jc,jk,jb)-z3d_in(jc,jk1  ,jb))*&
                                  (z3d_out(jc,jk,jb)-z3d_in(jc,jk1+1,jb))
                bot_idx(jc,jb) = jk
                l_found(jc) = .TRUE.
              ELSE IF (z3d_out(jc,jk,jb) <= z3d_in(jc,nlevs_in-1,jb)) THEN
                l_found(jc) = .TRUE.
                idx0(jc,jk,jb) = nlevs_in-1
              ELSE IF (z3d_out(jc,jk,jb) >  z3d_in(jc,2,jb)) THEN
                ! linear interpolation between two upper input levels or extrapolation beyond top level
                idx0(jc,jk,jb) = 1
                coef1(jc,jk,jb) = z3d_out(jc,jk,jb)-z3d_in(jc,1,jb)
                coef2(jc,jk,jb) = 0._wp
                coef3(jc,jk,jb) = 0._wp
                bot_idx(jc,jb) = jk
                l_found(jc) = .TRUE.
              ENDIF
            ENDIF
          ENDDO

#ifndef _OPENACC
          ! ACC: the following EXIT would be illegal within an OpenACC Kernel, thus we skip this CPU optimization
          IF (ALL(l_found(1:nlen))) THEN
            lfound_all = .TRUE.
            EXIT
          ENDIF
#endif

        ENDDO
        !$ACC END PARALLEL

#ifdef _OPENACC
        lfound_all = .TRUE.
        ! The following reduction must appear in its own small kernel as it did not work otherwise with Nvidia 21.2
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc) REDUCTION(.AND.: lfound_all)
        DO jc = 1, nlen
          lfound_all = lfound_all .AND. l_found(jc)
        ENDDO
        !$ACC END PARALLEL LOOP
#endif
        !$ACC WAIT

        IF (lfound_all) THEN
          jk_start = MIN(minval_1d(idx0(1:nlen,jk,jb), lacc=lzacc)+1,nlevs_in-2)
          jk_start = MAX(jk_start, 2) ! avoid out-of-bounds access
        ELSE
          ierror(jb) = ierror(jb) + 1
        ENDIF
      ENDDO

      jk_start = minval_1d(bot_idx(1:nlen,jb), lacc=lzacc) + 1

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = jk_start, nlevs_out
        DO jc = 1, nlen
          IF (jk >= bot_idx(jc,jb)+1) THEN
            coef1(jc,jk,jb) = 0._wp
            coef2(jc,jk,jb) = 0._wp
            coef3(jc,jk,jb) = 0._wp
            idx0 (jc,jk,jb) = nlevs_in
           ENDIF
        ENDDO
      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

    nerror = SUM(ierror)

    IF (nerror > 0) CALL finish("prepare_cubic_intp:",&
      "Error in computing cubic interpolation coefficients")

  END SUBROUTINE prepare_cubic_intp


  !-------------
  !>
  !! SUBROUTINE lin_intp
  !! Performs linear vertical interpolation of a 3D field
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Note: field values below the surface level of the input field are
  !! computed by linear extrapolation, using the gradient of the lowest (by default)
  !! 500 m above ground. Do not use for temperature!
  !!
  !! Setting l_loglin=.TRUE. activates logarithmic interpolation
  !!
  SUBROUTINE lin_intp(f3d_in, f3d_out,                      &
                      nblks, npromz, nlevs_in, nlevs_out,   &
                      wfac, idx0, bot_idx, wfacpbl1, kpbl1, &
                      wfacpbl2, kpbl2, l_loglin, l_pd_limit,&
                      l_extrapol, lower_limit, lacc         )


    ! Atmospheric fields
    REAL(wp), INTENT(IN)  :: f3d_in (:,:,:) ! input field
    REAL(wp), INTENT(OUT) :: f3d_out(:,:,:) ! output field (on ICON vertical grid)

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    ! Control switches
    LOGICAL,  INTENT(IN) :: l_loglin    ! switch for logarithmic interpolation
    LOGICAL,  INTENT(IN) :: l_pd_limit  ! switch for use of positive definite limiter
    LOGICAL,  INTENT(IN) :: l_extrapol  ! switch for use of downward extrapolation (no-gradient condition otherwise)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! if true use OpenACC

    REAL(wp), INTENT(IN), OPTIONAL :: lower_limit ! lower limit of variable

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jc
    INTEGER  :: nlen
    REAL(wp) :: zf_in_tr(nproma,nlevs_in), zf_in_lim(nproma,nlevs_in), z_limit, f3d_z1, f3d_z2, vgrad_f3d
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (PRESENT(lower_limit)) THEN
      z_limit = lower_limit
    ELSE
      z_limit = 0._wp
    ENDIF

    !$ACC DATA IF(lzacc) &
    !$ACC   PRESENT(f3d_in, f3d_out, wfac, idx0, bot_idx, kpbl1, wfacpbl1, kpbl2, wfacpbl2) &
    !$ACC   CREATE(zf_in_tr, zf_in_lim)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jc,zf_in_tr,zf_in_lim,f3d_z1,f3d_z2,vgrad_f3d) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jk = 1, nlevs_out
          DO jc = nlen+1, nproma
            f3d_out(jc,jk,jb)  = 0.0_wp
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk = 1, nlevs_in
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          IF (l_pd_limit) THEN
            zf_in_lim(jc, jk) = MAX(z_limit,f3d_in(jc, jk, jb))
          ELSE
            zf_in_lim(jc, jk) = f3d_in(jc, jk, jb)
          ENDIF
          IF (l_loglin) THEN
            zf_in_tr(jc, jk) = LOG(MAX(1.e-20_wp,zf_in_lim(jc, jk)))
          ELSE
            zf_in_tr(jc, jk) = zf_in_lim(jc, jk)
          ENDIF
        ENDDO
      ENDDO

      !$ACC LOOP SEQ
      DO jk = 1, nlevs_out
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(f3d_z1, f3d_z2, vgrad_f3d)
        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN

            ! linear interpolation
            f3d_out(jc,jk,jb) = wfac(jc,jk,jb)*zf_in_tr(jc,idx0(jc,jk,jb)) + &
              (1._wp-wfac(jc,jk,jb))*zf_in_tr(jc,idx0(jc,jk,jb)+1)

            IF (l_loglin) f3d_out(jc,jk,jb) = EXP(f3d_out(jc,jk,jb))

          ELSE IF (l_extrapol) THEN

            ! linear extrapolation, using the gradient between (by default) 1000 m and 500 m AGL
            ! Logarithmic computation is not used here because it would be numerically unstable for extrapolation

            ! Field value at height zpbl1
            f3d_z1 = wfacpbl1(jc,jb) *zf_in_lim(jc,kpbl1(jc,jb)  ) +  &
              (1._wp-wfacpbl1(jc,jb))*zf_in_lim(jc,kpbl1(jc,jb)+1)

            ! Field value at height zpbl2
            f3d_z2 = wfacpbl2(jc,jb) *zf_in_lim(jc,kpbl2(jc,jb)  ) +  &
              (1._wp-wfacpbl2(jc,jb))*zf_in_lim(jc,kpbl2(jc,jb)+1)

            ! vertical gradient
            vgrad_f3d = (f3d_z2-f3d_z1)/(zpbl2-zpbl1)

            ! wfac carries the (negative) extrapolation distance (in m) in this case
            f3d_out(jc,jk,jb) = zf_in_lim(jc,nlevs_in) + vgrad_f3d*wfac(jc,jk,jb)

          ELSE ! use no-gradient condition for extrapolation

            f3d_out(jc,jk,jb) = zf_in_lim(jc,nlevs_in)

          ENDIF
        ENDDO
      ENDDO

      IF (l_pd_limit) THEN
        !$ACC LOOP SEQ
        DO jk = 1, nlevs_out
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            f3d_out(jc, jk, jb) = MAX(z_limit,f3d_out(jc, jk, jb))
          ENDDO
        ENDDO
      ENDIF
      !$ACC END PARALLEL

    ENDDO
    !$ACC WAIT
    !$ACC END DATA
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE lin_intp


  !-------------
  !>
  !! SUBROUTINE pressure_intp
  !! Performs vertical interpolation of pressure
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! field of input data, virtual temperature and 3D height coordinate of output data
  !! Output: pressure field of output data
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  SUBROUTINE pressure_intp(pres_in, tempv_in, z3d_in, pres_out, z3d_out,                   &
                           nblks, npromz, nlevs_in, nlevs_out,                             &
                           wfac, idx0, bot_idx, wfacpbl1, kpbl1, wfacpbl2, kpbl2, zextrap, &
                           lacc                                                            )


    ! Input fields
    REAL(wp),         INTENT(IN)  :: pres_in  (:,:,:) ! pressure field of input data
    REAL(wp),         INTENT(IN)  :: tempv_in (:,:,:) ! virtual temperature of input data
    REAL(wp), TARGET, INTENT(IN)  :: z3d_in   (:,:,:) ! 3D height coordinate field of input data
    REAL(wp), TARGET, INTENT(IN)  :: z3d_out  (:,:,:) ! 3D height coordinate field of output data

    ! Output
    REAL(wp), INTENT(OUT) :: pres_out (:,:,:) ! pressure field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)     ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)  ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)     ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)  ! corresponding interpolation coefficient

    REAL(wp), OPTIONAL, INTENT(IN) :: zextrap(:,:)   ! AGL height from which downward extrapolation starts (in postprocesing mode)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! LOCAL VARIABLES

    CHARACTER(*), PARAMETER :: routine = 'mo_nh_vert_interp: pressure_intp'
    INTEGER  :: jb, jk, jc, jk1, nlen
    REAL(wp), DIMENSION(nproma,nlevs_out) :: dtvdz_down
    REAL(wp), DIMENSION(nproma)           :: tmsl, tsfc_mod, tempv1, tempv2, vtgrad_up, sfc_inv
    REAL(wp), ALLOCATABLE, TARGET         :: zgpot_out(:,:,:), zgpot_in(:,:,:)
    REAL(wp),              POINTER        :: z_in(:,:,:), z_out(:,:,:)
    REAL(wp) :: p_up, p_down, t_extr
    INTEGER  :: istat

!-------------------------------------------------------------------------

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (.NOT. PRESENT(zextrap) .AND. itype_pres_msl >= 3 ) THEN
      CALL finish(routine, "zextrap missing in argument list")
    ENDIF

    CALL assert_acc_device_only(routine, lacc)

    IF (.NOT. ldeepatmo) THEN
      z_in  => z3d_in
      z_out => z3d_out
    ELSE
      ALLOCATE(zgpot_in(nproma, nlevs_in,  nblks), zgpot_out(nproma, nlevs_out, nblks), STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, 'Allocation of zgpot failed')
      !$ACC ENTER DATA &
      !$ACC   CREATE(zgpot_in, zgpot_out)
      CALL deepatmo_htrafo(z_in=z3d_out, z_out=zgpot_out, nblks_nproma_npromz=[nblks,nproma,npromz], &
        & start_end_levels=[1,nlevs_out], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=lacc)
      IF (istat /= SUCCESS) CALL finish(routine, 'deepatmo_htrafo failed')
      CALL deepatmo_htrafo(z_in=z3d_in, z_out=zgpot_in, nblks_nproma_npromz=[nblks,nproma,npromz], &
        & start_end_levels=[1,nlevs_in], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=lacc)
      IF (istat /= SUCCESS) CALL finish(routine, 'deepatmo_htrafo failed')
      z_in  => zgpot_in
      z_out => zgpot_out
      ! Note: the heights above ground level, 'zpbl1', 'zpbl2', 'zextrap' and heights derived from them 
      ! have relatively low values (~ 1 km), so no deep-atmosphere modification is applied to them. 
      ! (Put another way: we regard 'zpbl1', 'zpbl2', and 'zextrap' to represent geopotential heights.)
    ENDIF

    !$ACC DATA PRESENT(pres_in, tempv_in, z3d_in, z3d_out, pres_out, wfac, idx0, bot_idx, kpbl1, wfacpbl1, kpbl2) &
    !$ACC   PRESENT(wfacpbl2, z_in, z_out, zextrap, num_lev, vct_a) &
    !$ACC   CREATE(dtvdz_down, tmsl, tsfc_mod, tempv1, tempv2, vtgrad_up, sfc_inv)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jk1,jc,dtvdz_down,p_up,p_down,tmsl,tsfc_mod,tempv1,tempv2,&
!$OMP            vtgrad_up,sfc_inv,t_extr) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
        DO jc = nlen+1, nproma
          pres_out(jc,:,jb)  = 0.0_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
      IF (itype_pres_msl == PRES_MSL_METHOD_GME) THEN

        ! Preparations for extrapolation below the ground: calculate temperature
        ! profile with artificial limits like in IFS. This temperature profile
        ! is NOT consistent with the temperature on pressure levels computed
        ! afterwards, implying that temperature and geopotential on pressure
        ! levels are not in hydrostatic balance. We are not happy with this
        ! solution, but it is used for the time being in order to be consistent
        ! with IFS and GME.

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          tsfc_mod(jc) = tempv_in(jc,nlevs_in,jb)
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z_in(jc,nlevs_in,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO

      ELSE IF ( itype_pres_msl >= 3 ) THEN

        ! Similar method to option 1, but we use the temperature at 150 m AGL
        ! for extrapolation (kpbl1 contains the required index in this case)

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen

          tsfc_mod(jc) = wfacpbl1(jc,jb) *tempv_in(jc,kpbl1(jc,jb),jb  ) + &
                  (1._wp-wfacpbl1(jc,jb))*tempv_in(jc,kpbl1(jc,jb)+1,jb) - &
                  dtdz_standardatm*(zextrap(jc,jb)-0.5_wp*vct_a(num_lev(1)))

          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z_in(jc,nlevs_in,jb)
          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF
        ENDDO

      ELSE

        ! Use a temperature profile that is (roughly) consistent with
        ! that used for temperature extrapolation

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          ! Virtual temperature at height zpbl1
          tempv1(jc) = wfacpbl1(jc,jb) *tempv_in(jc,kpbl1(jc,jb),jb  ) + &
                (1._wp-wfacpbl1(jc,jb))*tempv_in(jc,kpbl1(jc,jb)+1,jb)

          ! Virtual temperature at height zpbl2
          tempv2(jc) = wfacpbl2(jc,jb) *tempv_in(jc,kpbl2(jc,jb),jb  ) + &
                (1._wp-wfacpbl2(jc,jb))*tempv_in(jc,kpbl2(jc,jb)+1,jb)

          ! Vertical gradient between zpbl1 and zpbl2
          vtgrad_up(jc) = (tempv2(jc) - tempv1(jc))/(zpbl2 - zpbl1)

          ! Modified (extrapolated) surface temperature without inversion
          tsfc_mod(jc) = tempv1(jc) - zpbl1*vtgrad_up(jc)

          ! "surface inversion", defined by the difference between the extrapolated
          ! extrapolated temperature from above and the original input temperature
          sfc_inv(jc) = tsfc_mod(jc) - tempv_in(jc,nlevs_in,jb)

          ! Reduction of the surface inversion depending on the extrapolation
          ! distance. The surface inversion is fully restored for extrapolation distances
          ! up to zpbl1 and disregarded for distances larger than 3*zpbl1
          IF (z_in(jc,nlevs_in,jb) > 3._wp*zpbl1) THEN
            sfc_inv(jc) = 0._wp
          ELSE IF (z_in(jc,nlevs_in,jb) > zpbl1) THEN
            sfc_inv(jc) = sfc_inv(jc)*(1._wp - (z_in(jc,nlevs_in,jb)-zpbl1)/(2._wp*zpbl1))
          ENDIF

          tsfc_mod(jc) = tsfc_mod(jc) - sfc_inv(jc)

          ! Limitation of very cold temperatures according to GME method
          IF (tsfc_mod(jc) < t_low) tsfc_mod(jc) = 0.5_wp*(t_low+tsfc_mod(jc))

          ! Estimated temperature at mean sea level
          tmsl(jc) = tsfc_mod(jc) - dtdz_standardatm*z_in(jc,nlevs_in,jb)

          IF (tmsl(jc) > t_high) THEN
            IF (tsfc_mod(jc) > t_high) THEN
              tsfc_mod(jc) = 0.5_wp*(t_high+tsfc_mod(jc))
              tmsl(jc)     = tsfc_mod(jc)
            ELSE
              tmsl(jc)     = t_high
            ENDIF
          ENDIF

        ENDDO

      ENDIF

      ! Compute vertical gradients of virtual potential temperature
      !$ACC LOOP SEQ
      DO jk = 1, nlevs_out

        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen           ! The vertical gradient is needed for downward extrapolation only

          IF (z_in(jc,nlevs_in,jb) > 1._wp) THEN
            dtvdz_down(jc,jk) = (tsfc_mod(jc)-tmsl(jc))/z_in(jc,nlevs_in,jb)

          ELSE ! avoid pathological results at grid points below sea level
            dtvdz_down(jc,jk) = dtdz_standardatm

          ENDIF

        ENDDO

      ENDDO

      ! Now compute pressure on target grid
      !$ACC LOOP SEQ
      DO jk = 1, nlevs_out

        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk1, p_up, p_down, t_extr)
        DO jc = 1, nlen
          IF (jk <= bot_idx(jc,jb)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation:
            ! we use only the formula for constant (layer-averaged) temperature in order not to
            ! require an already interpolated (virtual) temperature on the output grid;
            ! the linear weighting between the upward and downward integrals nevertheless
            ! ensures negligibly small errors
            jk1 = idx0(jc,jk,jb)

            p_up = pres_in(jc,jk1,jb)*EXP(-grav*(z_out(jc,jk,jb)-z_in(jc,jk1,jb)) / &
              (rd*0.5_wp*(tempv_in(jc,jk1,jb)+tempv_in(jc,jk1+1,jb))) )

            p_down = pres_in(jc,jk1+1,jb)*EXP(-grav*(z_out(jc,jk,jb)-z_in(jc,jk1+1,jb)) / &
              (rd*0.5_wp*(tempv_in(jc,jk1,jb)+tempv_in(jc,jk1+1,jb))) )

            ! apply inverse-distance weighting between top-down and bottom-up integrated value
            ! except in the case of extrapolation above the top of the input data
            IF (jk1 == 1 .AND. wfac(jc,jk,jb) > 1._wp) THEN
              pres_out(jc,jk,jb) = p_up
            ELSE
              pres_out(jc,jk,jb) = wfac(jc,jk,jb)*p_up + (1._wp-wfac(jc,jk,jb))*p_down
            ENDIF

          ELSE ! downward extrapolation based on the vertical gradient computed above
            t_extr = tsfc_mod(jc) + (z_out(jc,jk,jb)-z_in(jc,nlevs_in,jb))*dtvdz_down(jc,jk)
            IF (ABS(dtvdz_down(jc,jk)) > dtvdz_thresh) THEN
              p_down = pres_in(jc,nlevs_in,jb)*EXP(-grav/(rd*dtvdz_down(jc,jk))* &
                LOG(t_extr/tsfc_mod(jc)) )
            ELSE
              p_down = pres_in(jc,nlevs_in,jb)*EXP(-grav*(z_out(jc,jk,jb)-z_in(jc,nlevs_in,jb)) / &
                (rd*0.5_wp*(t_extr+tsfc_mod(jc))) )
            ENDIF
            pres_out(jc,jk,jb) = p_down
          ENDIF
        ENDDO

      ENDDO
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

    NULLIFY(z_in, z_out)
    IF (ldeepatmo) THEN
      !$ACC EXIT DATA &
      !$ACC   DELETE(zgpot_in, zgpot_out)
      DEALLOCATE(zgpot_in, zgpot_out, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, 'Deallocation of zgpot failed')
    ENDIF

  END SUBROUTINE pressure_intp


  !-------------
  !>
  !! SUBROUTINE pressure_intp_initmode
  !! Performs vertical interpolation of pressure
  !!
  !! Required input fields: pressure, virtual temperature and 3D height coordinate
  !! field of input data, virtual temperature and 3D height coordinate of output data
  !! Output: pressure field of output data
  !!
  !! Method: piecewise analytical integration of the hydrostatic equation
  !!
  SUBROUTINE pressure_intp_initmode(pres_in, tempv_in, z3d_in, pres_out, tempv_out, z3d_out, &
    &                               nblks, npromz, nlevs_out, nlevs_in, wfac, idx0, bot_idx, &
    &                               opt_lmask                                                )


    ! Input fields
    REAL(wp),         INTENT(IN)  :: pres_in  (:,:,:) ! pressure field of input data
    REAL(wp),         INTENT(IN)  :: tempv_in (:,:,:) ! virtual temperature of input data
    REAL(wp), TARGET, INTENT(IN)  :: z3d_in   (:,:,:) ! 3D height coordinate field of input data
    REAL(wp),         INTENT(IN)  :: tempv_out(:,:,:) ! virtual temperature of output data
    REAL(wp), TARGET, INTENT(IN)  :: z3d_out  (:,:,:) ! 3D height coordinate field of output data

    ! Output
    REAL(wp), INTENT(OUT) :: pres_out (:,:,:) ! pressure field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: wfac(:,:,:)    ! weighting factor of upper level
    INTEGER , INTENT(IN) :: idx0(:,:,:)    ! index of upper level
    INTEGER , INTENT(IN) :: bot_idx(:,:)   ! index of lowest level for which interpolation is possible
    LOGICAL, OPTIONAL,  INTENT(IN) :: opt_lmask(:,:)

    ! LOCAL CONSTANTS

    REAL(wp), PARAMETER :: TOL = 1e-12_wp

    ! LOCAL VARIABLES

    REAL(wp), ALLOCATABLE, TARGET  :: zgpot_out(:,:,:), zgpot_in(:,:,:)
    REAL(wp),              POINTER :: z_in(:,:,:), z_out(:,:,:)

    INTEGER  :: jb, jk, jc, jk1, nlen
    REAL(wp), DIMENSION(nproma,nlevs_out) :: dtvdz_up, dtvdz_down
    REAL(wp) :: p_up, p_down, inv_scal_hgt
    LOGICAL  :: lmask(nproma)
    INTEGER  :: istat

!-------------------------------------------------------------------------
    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    ! inverse scale height for filling pressure on data-void grid points with artificial values
    inv_scal_hgt = grav/(rd*fill_temp)

    IF (.NOT. ldeepatmo) THEN
      z_in  => z3d_in
      z_out => z3d_out
    ELSE
      ALLOCATE(zgpot_in(nproma, nlevs_in,  nblks), zgpot_out(nproma, nlevs_out, nblks), STAT=istat)
      IF (istat /= SUCCESS) CALL finish('mo_nh_vert_interp: pressure_intp_initmode', 'Allocation of zgpot failed') 
      ! Compute geopotential heights in case of the deep atmosphere
      CALL deepatmo_htrafo(z_in=z3d_out, z_out=zgpot_out, nblks_nproma_npromz=[nblks,nproma,npromz], &
        & start_end_levels=[1,nlevs_out], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=.FALSE.)
      IF (istat /= SUCCESS) CALL finish('mo_nh_vert_interp: pressure_intp_initmode', 'deepatmo_htrafo failed')
      CALL deepatmo_htrafo(z_in=z3d_in, z_out=zgpot_in, nblks_nproma_npromz=[nblks,nproma,npromz], &
        & start_end_levels=[1,nlevs_in], radius=grid_sphere_radius, trafo_type='z2zgpot', ierror=istat, lacc=.FALSE.)
      IF (istat /= SUCCESS) CALL finish('mo_nh_vert_interp: pressure_intp_initmode', 'deepatmo_htrafo failed')
      z_in  => zgpot_in
      z_out => zgpot_out
    ENDIF

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,nlen,jk,jk1,jc,dtvdz_up,dtvdz_down,p_up,p_down,lmask) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        pres_out(nlen+1:nproma,:,jb)  = 0.0_wp
      ENDIF

      IF (PRESENT(opt_lmask)) THEN
        lmask(:) = opt_lmask(:,jb)
        IF (.NOT. ANY(lmask(:)) ) THEN
          pres_out(:,:,jb) = p0sl_bg*EXP(-z_out(:,:,jb)*inv_scal_hgt)
          CYCLE
        ENDIF
      ELSE
        lmask(:) = .TRUE.
      ENDIF

      ! Compute vertical gradients of virtual potential temperature
      DO jk = 1, nlevs_out

        DO jc = 1, nlen

          IF (jk <= bot_idx(jc,jb)) THEN
            jk1 = idx0(jc,jk,jb)

            IF (ABS(z_in  (jc,jk1,jb)-z_out  (jc,jk,jb)) > TOL) THEN
              dtvdz_up(jc,jk) = (tempv_in(jc,jk1,jb)-tempv_out(jc,jk,jb)) / &
                &               (z_in    (jc,jk1,jb)-z_out    (jc,jk,jb))
            ELSE
              dtvdz_up(jc,jk) = 0._wp
            ENDIF

            ! Paranoia
            IF (ABS(z_in  (jc,jk1+1,jb)-z_out  (jc,jk,jb)) > TOL) THEN
              dtvdz_down(jc,jk) = (tempv_in(jc,jk1+1,jb)-tempv_out(jc,jk,jb)) / &
                &                 (z_in    (jc,jk1+1,jb)-z_out    (jc,jk,jb))
            ELSE
              dtvdz_down(jc,jk) = 0._wp
            END IF

          ELSE ! downward extrapolation; only dtvdz_down is needed

            dtvdz_down(jc,jk) = (tempv_out(jc,jk-1,jb)-tempv_out(jc,jk,jb)) / &
              (z_out(jc,jk-1,jb)-z_out(jc,jk,jb))

          ENDIF
        ENDDO

      ENDDO

      ! Now compute pressure on target grid
      DO jk = 1, nlevs_out

        DO jc = 1, nlen

          IF (.NOT. lmask(jc)) CYCLE
          IF (jk <= bot_idx(jc,jb)) THEN

            ! interpolation based on piecewise analytical integration of the hydrostatic equation
            jk1 = idx0(jc,jk,jb)

            IF (ABS(dtvdz_up(jc,jk)) > dtvdz_thresh) THEN
              p_up = pres_in(jc,jk1,jb)*EXP(-grav/(rd*dtvdz_up(jc,jk)) * &
                LOG(tempv_out(jc,jk,jb)/tempv_in(jc,jk1,jb)) )
            ELSE
              p_up = pres_in(jc,jk1,jb)*EXP(-grav*(z_out(jc,jk,jb)-z_in(jc,jk1,jb)) / &
                (rd*0.5_wp*(tempv_out(jc,jk,jb)+tempv_in(jc,jk1,jb))) )
            ENDIF

            IF (ABS(dtvdz_down(jc,jk)) > dtvdz_thresh) THEN
              p_down = pres_in(jc,jk1+1,jb)*EXP(-grav/(rd*dtvdz_down(jc,jk))* &
                LOG(tempv_out(jc,jk,jb)/tempv_in(jc,jk1+1,jb)) )
            ELSE
              p_down = pres_in(jc,jk1+1,jb)*EXP(-grav*(z_out(jc,jk,jb)-z_in(jc,jk1+1,jb)) / &
                (rd*0.5_wp*(tempv_out(jc,jk,jb)+tempv_in(jc,jk1+1,jb))) )
            ENDIF

            ! Finally, apply inverse-distance weighting between top-down and bottom-up integrated value
            ! except in the case of extrapolation above the top of the input data
            IF (jk1 == 1 .AND. wfac(jc,jk,jb) > 1._wp) THEN
              pres_out(jc,jk,jb) = p_up
            ELSE
              pres_out(jc,jk,jb) = wfac(jc,jk,jb)*p_up + (1._wp-wfac(jc,jk,jb))*p_down
            ENDIF

          ELSE ! downward extrapolation

            IF (ABS(dtvdz_down(jc,jk)) > dtvdz_thresh) THEN
              p_down = pres_out(jc,jk-1,jb)*EXP(-grav/(rd*dtvdz_down(jc,jk))* &
                LOG(tempv_out(jc,jk,jb)/tempv_out(jc,jk-1,jb)) )
            ELSE
              p_down = pres_out(jc,jk-1,jb)*EXP(-grav*(z_out(jc,jk,jb)-z_out(jc,jk-1,jb)) / &
                (rd*0.5_wp*(tempv_out(jc,jk,jb)+tempv_out(jc,jk-1,jb))) )
            ENDIF

            pres_out(jc,jk,jb) = p_down
          ENDIF
        ENDDO

      ENDDO


      ! Fill data-void grid points
      IF (PRESENT(opt_lmask)) THEN
        DO jc = 1, nlen
          IF (.NOT. lmask(jc)) THEN
            pres_out(jc,:,jb) = p0sl_bg*EXP(-z_out(jc,:,jb)*inv_scal_hgt)
          ENDIF
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    NULLIFY(z_in, z_out)
    IF (ldeepatmo) THEN
      DEALLOCATE(zgpot_in, zgpot_out, STAT=istat)
      IF (istat /= SUCCESS) CALL finish('mo_nh_vert_interp: pressure_intp_initmode', 'Deallocation of zgpot failed') 
    ENDIF

  END SUBROUTINE pressure_intp_initmode


  !-------------
  !>
  !! SUBROUTINE compute_slope
  !! Computes the slope of the mass points, which is an optional input
  !! to the temperature interpolation routine
  !!
  !! Required input fields: topography fields, patch and interpolation state
  !! Output: slope
  !!
  SUBROUTINE compute_slope(p_patch, p_int, topo_c, slope_c)

    TYPE(t_patch),          INTENT(INOUT)    :: p_patch
    TYPE(t_int_state),      INTENT(IN)       :: p_int

    ! Topography data
    REAL(wp), INTENT(IN) :: topo_c(:,:)  ! topography height of mass points

    ! Output
    REAL(wp), INTENT(OUT) :: slope_c(:,:) ! surface slope of mass points

    ! LOCAL VARIABLES

    INTEGER  :: jb, jc, je
    INTEGER  :: i_startblk, i_startidx, i_endidx, nblks_c, nblks_e
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_c) :: z_topo_c, slope_abs_c
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_v) :: z_topo_v
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_e) :: slope_norm, slope_tang, slope_abs_e

!-------------------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! Initialization - no slope can be computed on lateral boundary points
    slope_c(:,:)    = 0._wp

    slope_norm (:,1,:) = 0._wp
    slope_tang (:,1,:) = 0._wp
    slope_abs_e(:,1,:) = 0._wp
    z_topo_v   (:,1,:) = 0._wp

    z_topo_c(:,1,:) = topo_c(:,:)

    ! Compute auxiliary topography at verices
    CALL cells2verts_scalar(z_topo_c, p_patch, p_int%cells_aw_verts, z_topo_v, 1, 1)

    ! Compute slopes
    CALL grad_fd_norm ( z_topo_c, p_patch, slope_norm, 1, 1)
    CALL grad_fd_tang ( z_topo_v, p_patch, slope_tang, 1, 1)

    i_startblk = p_patch%edges%start_blk(2,1)

    DO jb = i_startblk, nblks_e

      CALL get_indices_e(p_patch, jb, i_startblk, nblks_e, i_startidx, i_endidx, 2)

      DO je = i_startidx, i_endidx
        slope_abs_e(je,1,jb) = SQRT(slope_norm(je,1,jb)**2+slope_tang(je,1,jb)**2)
      ENDDO
    ENDDO

    ! Interpolate absolute slope to mass points
    CALL edges2cells_scalar(slope_abs_e, p_patch, p_int%e_bln_c_s, slope_abs_c, 1, 1)

    i_startblk = p_patch%cells%start_blk(2,1)

    DO jb = i_startblk, nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx
        slope_c(jc,jb) = slope_abs_c(jc,1,jb)
      ENDDO
    ENDDO

    CALL sync_patch_array(SYNC_C,p_patch,slope_c)

  END SUBROUTINE compute_slope


  !-------------
  !>
  !! SUBROUTINE temperature_intp
  !! Performs vertical interpolation and extrapolation of temperature with
  !! special treatment of surface inversions
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin/cubic_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Performs cubic interpolation where possible, turning to linear interpolation
  !! close to the surface
  !! The most important ingredient of the refined temperature extrapolation
  !! is to remove the surface inversion, if present, before the extrapolation
  !! and to add it again afterwards with a variable weighting coefficient,
  !! accounting for the slope of the target grid point and for the height
  !! difference between source and target grid
  !!
  SUBROUTINE temperature_intp(temp_in, temp_out, z3d_in, z3d_out,            &
                             nblks, npromz, nlevs_in, nlevs_out,             &
                             coef1, coef2, coef3, wfac_lin,                  &
                             idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,   &
                             wfacpbl1, kpbl1, wfacpbl2, kpbl2,               &
                             l_hires_corr, l_restore_sfcinv, extrapol_dist,  &
                             l_pz_mode, zextrap, slope, opt_lmask, lacc)


    ! Atmospheric fields
    REAL(wp), INTENT(IN)  :: temp_in (:,:,:) ! input temperature field
    REAL(wp), INTENT(OUT) :: temp_out(:,:,:) ! output temperature field

    ! Coordinate fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)  ! 3D height coordinate field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: coef1(:,:,:)    ! coefficient for linear term
    REAL(wp), INTENT(IN) :: coef2(:,:,:)    ! coefficient for quadratic term
    REAL(wp), INTENT(IN) :: coef3(:,:,:)    ! coefficient for cubic term
    REAL(wp), INTENT(IN) :: wfac_lin(:,:,:) ! weighting factor for linear interpolation

    INTEGER , INTENT(IN) :: idx0_cub(:,:,:) ! index of upper level, cubic interpolation
    INTEGER , INTENT(IN) :: idx0_lin(:,:,:) ! index of upper level, linear interpolation

    INTEGER , INTENT(IN) :: bot_idx_cub(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: bot_idx_lin(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    REAL(wp), OPTIONAL, INTENT(IN) :: zextrap(:,:)    ! AGL height from which downward extrapolation starts (in postprocesing mode)

    ! Logical switch if slope-based reduction of surface inversion is to be performed
    ! (recommended when target model has a much finer resolution than source model)
    LOGICAL, INTENT(IN) :: l_hires_corr

    ! Logical switch if surface inversion is to be restored after downward extrapolation
    ! (may be set to .FALSE. when interpolating on constant height/pressure levels)
    LOGICAL, INTENT(IN) :: l_restore_sfcinv

    ! Logical switch if interpolation mode to height/pressure levels is to be used
    LOGICAL, INTENT(IN) :: l_pz_mode

    REAL(wp), INTENT(IN) :: extrapol_dist ! Maximum extrapolation distance using the local vertical gradient
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    REAL(wp), INTENT(IN), OPTIONAL :: slope(:,:)  ! slope of mass points
    LOGICAL, OPTIONAL,  INTENT(IN) :: opt_lmask(:,:)

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, sfcinv, tmsl_max

    REAL(wp), DIMENSION(nproma) :: temp1, temp2, vtgrad_up, zdiff_inout, &
                                   redinv1, redinv2, tmsl, tmsl_mod
    LOGICAL , DIMENSION(nproma) :: l_found

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, sfc_inv, temp_mod, g1, g2, g3
    REAL(wp), DIMENSION(nproma,nlevs_in-1) :: zalml_in_d
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (l_hires_corr .AND. .NOT. PRESENT(slope)) CALL finish("temperature_intp:",&
      "slope correction requires slope data as input")

    IF (.NOT. PRESENT(zextrap) .AND. l_pz_mode .AND. itype_pres_msl >= 3) THEN
      CALL finish("temperature_intp:", "zextrap missing in argument list")
    ENDIF

    ! Artificial upper limit on sea-level temperature for so-called plateau correction
    tmsl_max = 298._wp

!$OMP PARALLEL private(jc, jk, zalml_in, zalml_out)

    !$ACC DATA IF(lzacc) &
    !$ACC   PRESENT(num_lev) &
    !$ACC   CREATE(ik1, temp1, temp2, vtgrad_up, zdiff_inout) &
    !$ACC   CREATE(redinv1, redinv2, tmsl, tmsl_mod) &
    !$ACC   CREATE(l_found, zalml_in, sfc_inv, temp_mod, g1, g2, g3) &
    !$ACC   CREATE(zalml_out, zalml_in_d) &
    !$ACC   PRESENT(temp_in, temp_out, z3d_in, z3d_out, coef1, coef2, coef3, wfac_lin, idx0_cub, idx0_lin) &
    !$ACC   PRESENT(bot_idx_cub, bot_idx_lin, kpbl1, wfacpbl1, kpbl2, wfacpbl2, zextrap, slope, opt_lmask)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO jk = 1, nlevs_in
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nproma
        zalml_in(jc, jk) = (jc - 1) * nlevs_in + jk
      END DO
    END DO
    !$ACC LOOP SEQ
    DO jk = 1, nlevs_out
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nproma
        zalml_out(jc, jk) = (jc - 1) * nlevs_out + jk
      END DO
    END DO
    !$ACC END PARALLEL

!-------------------------------------------------------------------------

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,sfcinv,    &
!$OMP            temp1,temp2,tmsl,tmsl_mod,vtgrad_up,zdiff_inout,redinv1,redinv2,l_found,&
!$OMP            zalml_in_d,temp_mod,sfc_inv,g1,g2,g3) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc=nlen+1, nproma
          temp_out(jc,:,jb)  = 0.0_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      IF (PRESENT(opt_lmask)) THEN
        CALL assert_acc_host_only("temperature_intp:opt_lmask", lacc)
        IF (.NOT. ANY(opt_lmask(:,jb))) THEN
          temp_out(:,:,jb) = fill_temp
          CYCLE
        ENDIF
      ENDIF

      IF (l_pz_mode .AND. itype_pres_msl < 3) THEN

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk1 = 1, nlevs_in
          DO jc = 1, nlen

            ! just copy temp_in to temp_mod
            temp_mod(jc,jk1) = temp_in(jc,jk1,jb)

          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ELSE IF (l_pz_mode) THEN ! new IFS method: extrapolate downward from 150 m AGL
                               ! with the standard atmosphere gradient

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          ! Temperature at height zpbl1 (150 m AGL in this case)
          temp1(jc) = wfacpbl1(jc,jb) *temp_in(jc,kpbl1(jc,jb),jb  ) + &
               (1._wp-wfacpbl1(jc,jb))*temp_in(jc,kpbl1(jc,jb)+1,jb)
        ENDDO

        !$ACC LOOP SEQ
        DO jk1 = 1, nlevs_in
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            IF (jk1 <= kpbl1(jc,jb)) THEN ! just use the input temperature
              temp_mod(jc,jk1) = temp_in(jc,jk1,jb)
            ELSE ! extrapolate downward from 150 m AGL with 6.5 K/km
              temp_mod(jc,jk1) = temp1(jc) - dtdz_standardatm*(zextrap(jc,jb) -       &
                0.5_wp*vct_a(num_lev(1)) - (z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)) )
            ENDIF
          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ELSE
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          ! Temperature at height zpbl1
          temp1(jc) = wfacpbl1(jc,jb) *temp_in(jc,kpbl1(jc,jb),jb  ) + &
               (1._wp-wfacpbl1(jc,jb))*temp_in(jc,kpbl1(jc,jb)+1,jb)

          ! Temperature at height zpbl2
          temp2(jc) = wfacpbl2(jc,jb) *temp_in(jc,kpbl2(jc,jb),jb  ) + &
               (1._wp-wfacpbl2(jc,jb))*temp_in(jc,kpbl2(jc,jb)+1,jb)

          ! Vertical gradient between zpbl1 and zpbl2
          vtgrad_up(jc) = (temp2(jc) - temp1(jc))/(zpbl2 - zpbl1)

          ! Set reasonable limits
          vtgrad_up(jc) = MAX(vtgrad_up(jc),-8.5e-3_wp)
          vtgrad_up(jc) = MIN(vtgrad_up(jc),-1.5e-3_wp)

          ! height distance between lowest input and output grid point
          ! (negative if extrapolation takes place)
          zdiff_inout(jc) = z3d_out(jc,nlevs_out,jb) - z3d_in(jc,nlevs_in,jb)
        ENDDO

        !$ACC LOOP SEQ
        DO jk1 = 1, nlevs_in
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen

            ! Height above lowest model level
            zalml_in(jc,jk1) = z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)

            ! Modified temperature with surface inversion removed
            IF (zalml_in(jc,jk1) < zpbl1) THEN
              temp_mod(jc,jk1) = temp1(jc)+vtgrad_up(jc)*(zalml_in(jc,jk1)-zpbl1)
            ELSE
              temp_mod(jc,jk1) = temp_in(jc,jk1,jb)
            ENDIF

            ! "surface inversion", defined by the difference between the extrapolated
            ! extrapolated temperature from above and the original input temperature
            !
            ! Note: the extrapolated temperature may also be slightly colder than the original
            ! one in the presence of (super-)adiabatic temperature gradients. In such
            ! cases, using vtgrad_up for downward extrapolation is still appropriate because
            ! the atmosphere in narrow mountain valleys usually does not become adiabatic
            ! because of slope heating
            !
            sfc_inv(jc,jk1) = temp_mod(jc,jk1) - temp_in(jc,jk1,jb)

          ENDDO
        ENDDO
        !$ACC END PARALLEL

        jk_start_in = start_idx_threshold(zpbl1, zalml_in(:,:), nlen, nlevs_in, lzacc)

      ENDIF

      ! Reduce the surface inversion strength depending on the slope of the target point and
      ! on the height difference between the source and target topography.
      !
      ! The following empirical modifications are designed for the case that the target
      ! model has a significantly finer spatial resolution than the source model,
      ! implying that nocturnal surface inversions should be removed on grid points lying
      ! higher than on the source grid (i.e. mountain ranges not resolved in the source model)
      ! and on slope points because cold air drainage in reality lets the cold air
      ! accumulate at the valley bottom.

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      IF (l_hires_corr) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen

          IF (slope(jc,jb) <= 5.e-2_wp) THEN ! 50 m/km
            redinv1(jc) = 1._wp
          ELSE IF (slope(jc,jb) <= 3.5e-1_wp) THEN ! 350 m/km
            redinv1(jc) = 1._wp - LOG(slope(jc,jb)/5.e-2_wp)/LOG(10._wp)
          ELSE
            redinv1(jc) = 1._wp - LOG(7._wp)/LOG(10._wp)
          ENDIF

          ! zdiff_inout > 0 means that the target grid point is higher than the source point
          IF (zdiff_inout(jc) <= 100._wp) THEN
            redinv2(jc) = 1._wp + MAX(0._wp,MIN(1.5_wp,-1.e-3_wp*zdiff_inout(jc)))
          ELSE IF (zdiff_inout(jc) <= 300._wp) THEN
            redinv2(jc) = 1._wp - (zdiff_inout(jc)-100._wp)/200._wp
          ELSE
            redinv2(jc) = 0._wp
          ENDIF
        ENDDO

        ! Reduce surface inversion if sfc_inv > 0 (i.e. there is really enhanced static stability)
        !$ACC LOOP SEQ
        DO jk1 = jk_start_in, nlevs_in
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            IF (sfc_inv(jc,jk1) > 0._wp) THEN
              sfc_inv(jc,jk1) = sfc_inv(jc,jk1)*MIN(1._wp, redinv1(jc)*redinv2(jc))
            ENDIF
          ENDDO
        ENDDO

      ENDIF

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nlen
        ! Compute vertical gradients of input data
        g1(jc,1) = (temp_mod(jc,1)-temp_mod(jc,2))/ &
             (z3d_in(jc,1,jb)-z3d_in(jc,2,jb))
        g1(jc,2) = (temp_mod(jc,2)-temp_mod(jc,3))/ &
             (z3d_in(jc,2,jb)-z3d_in(jc,3,jb))
        ! Compute vertical gradients of gradients
        g2(jc,2-1) = (g1(jc,1)-g1(jc,2))/(z3d_in(jc,1,jb)-z3d_in(jc,3,jb))
      ENDDO
      !$ACC LOOP SEQ
      DO jk1 = 3, nlevs_in-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          ! Compute vertical gradients of input data
          g1(jc,jk1) = (temp_mod(jc,jk1 )-temp_mod(jc,jk1+1 ))/ &
                       (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
          ! Compute vertical gradients of gradients
          g2(jc,jk1-1) = (g1(jc,jk1-1)-g1(jc,jk1))/(z3d_in(jc,jk1-1,jb)-z3d_in(jc,jk1+1,jb))
          ! Compute third-order vertical gradients
          g3(jc,jk1-2) = (g2(jc,jk1-2)-g2(jc,jk1-1))/(z3d_in(jc,jk1-2,jb)-z3d_in(jc,jk1+1,jb))
        ENDDO
      ENDDO

      ! Now perform vertical interpolation, based on the modified temperature field
      IF (.NOT. l_pz_mode) THEN ! mode for interpolating initial data

        !$ACC LOOP SEQ
        DO jk = 1, nlevs_out
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk1)
          DO jc = 1, nlen
            IF (jk <= bot_idx_cub(jc,jb)) THEN

              ! cubic interpolation
              jk1 = idx0_cub(jc,jk,jb)
              temp_out(jc,jk,jb) = temp_mod(jc,jk1)           + coef1(jc,jk,jb)*g1(jc,jk1) + &
                                   coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

            ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

              ! linear interpolation
              jk1 = idx0_lin(jc,jk,jb)
              temp_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*temp_mod(jc,jk1) + &
                (1._wp-wfac_lin(jc,jk,jb))*temp_mod(jc,jk1+1)

            ELSE

              ! linear extrapolation using the upper temperature gradient; wfac_lin
              ! carries the extrapolation distance
              IF (wfac_lin(jc,jk,jb) > extrapol_dist) THEN
                temp_out(jc,jk,jb) = temp_mod(jc,nlevs_in) + wfac_lin(jc,jk,jb)*vtgrad_up(jc)
              ELSE
                temp_out(jc,jk,jb) = temp_mod(jc,nlevs_in) + extrapol_dist*vtgrad_up(jc) +  &
                                     (wfac_lin(jc,jk,jb)-extrapol_dist)*dtdz_standardatm
              ENDIF

            ENDIF

            ! Height above lowest model level - needed for restoring the surface inversion
            zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

          ENDDO
        ENDDO

      ELSE ! mode for interpolation to pressure or height levels

        ! Some precomputations that do not depend on the model level
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen

          ! extrapolated mean sea level temperature
          tmsl(jc) = temp_mod(jc,nlevs_in) - dtdz_standardatm * z3d_in(jc,nlevs_in,jb)

          ! "Plateau correction"
          IF (z3d_in(jc,nlevs_in,jb) > 2000._wp .AND. tmsl(jc) > tmsl_max) THEN
            IF (z3d_in(jc,nlevs_in,jb) < 2500._wp) THEN
              tmsl_mod(jc) = 0.002_wp*((2500._wp-z3d_in(jc,nlevs_in,jb))*tmsl(jc)+  &
                (z3d_in(jc,nlevs_in,jb)-2000._wp)*tmsl_max)
            ELSE
              tmsl_mod(jc) = tmsl_max
            ENDIF
          ELSE
            tmsl_mod(jc) = tmsl(jc)
          ENDIF
        ENDDO

        !$ACC LOOP SEQ
        DO jk = 1, nlevs_out
          !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk1)
          DO jc = 1, nlen
            IF (jk <= bot_idx_cub(jc,jb)) THEN

              ! cubic interpolation
              jk1 = idx0_cub(jc,jk,jb)
              temp_out(jc,jk,jb) = temp_mod(jc,jk1)           + coef1(jc,jk,jb)*g1(jc,jk1) + &
                                   coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

            ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

              ! linear interpolation
              jk1 = idx0_lin(jc,jk,jb)
              temp_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*temp_mod(jc,jk1) + &
                (1._wp-wfac_lin(jc,jk,jb))*temp_mod(jc,jk1+1)

            ELSE ! apply empirical correction computed above for extrapolation below the ground

              temp_out(jc,jk,jb) = temp_mod(jc,nlevs_in) + wfac_lin(jc,jk,jb)* &
                (temp_mod(jc,nlevs_in)-tmsl_mod(jc))/z3d_in(jc,nlevs_in,jb)

            ENDIF

          ENDDO
        ENDDO

      ENDIF
      !$ACC END PARALLEL

      ! Finally, subtract surface inversion from preliminary temperature field
      IF (l_restore_sfcinv .AND. .NOT. l_pz_mode) THEN

        jk_start_out = start_idx_threshold(zpbl1, zalml_out(:,:), nlen, nlevs_out, lzacc)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlevs_in-1
          DO jc = 1, nlen
            zalml_in_d(jc,jk) = 1.0_wp / (zalml_in(jc,jk) - zalml_in(jc,jk+1))
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            l_found(jc) = .FALSE.
          ENDDO
          !$ACC LOOP SEQ
          DO jk1 = jk_start, nlevs_in-1
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wfac, sfcinv)
            DO jc = 1, nlen
              IF(.NOT. l_found(jc)) THEN
                IF (zalml_out(jc,jk) >= zpbl1) THEN
                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk_start
                ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                         zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                  wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1)) * zalml_in_d(jc,jk1)
                  sfcinv = wfac*sfc_inv(jc,jk1) + (1._wp-wfac)*sfc_inv(jc,jk1+1)

                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk1

                  temp_out(jc,jk,jb) =  temp_out(jc,jk,jb) - sfcinv

                ELSE IF (zalml_out(jc,jk) > zalml_in(jc,jk_start)) THEN
                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk_start
                ENDIF
              ENDIF
            ENDDO
#ifndef _OPENACC
            ! ACC: the following EXIT would be illegal within an OpenACC Kernel, thus we skip this CPU optimization
            IF (ALL(l_found(1:nlen))) EXIT
#endif
          ENDDO
          !$ACC END PARALLEL
          jk_start = minval_1d(ik1(1: nlen), lacc=lzacc)
        ENDDO
      ENDIF

      ! Fill data-void grid points
      IF (PRESENT(opt_lmask)) THEN
        CALL assert_acc_host_only("temperature_intp:opt_lmask", lacc)
        DO jc = 1, nlen
          IF (.NOT. opt_lmask(jc,jb)) THEN
            temp_out(jc,:,jb) = fill_temp
          ENDIF
        ENDDO
      ENDIF

    ENDDO
    !$ACC WAIT
    !$ACC END DATA
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE temperature_intp


  !-------------
  !>
  !! SUBROUTINE uv_intp
  !! Performs vertical interpolation and extrapolation of horizontal wind components
  !! with special treatment of boundary layer effects
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin/cubic_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Performs cubic interpolation where possible, turning to linear interpolation
  !! close to the surface
  !! Boundary-layer treatment follows a similar reasoning as for temperature,
  !! but there are different consistency checks and limitations.
  !! In particular, wind components are requested not to change sign when being
  !! extrapolated downward
  !!
  SUBROUTINE uv_intp(uv_in, uv_out, z3d_in, z3d_out,               &
                     nblks, npromz, nlevs_in, nlevs_out,           &
                     coef1, coef2, coef3, wfac_lin,                &
                     idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin, &
                     wfacpbl1, kpbl1, wfacpbl2, kpbl2,             &
                     l_hires_intp, l_restore_fricred,              &
                     extrap_limit, lacc)


    ! Atmospheric fields
    REAL(wp), INTENT(IN)  :: uv_in (:,:,:) ! input wind (component) field
    REAL(wp), INTENT(OUT) :: uv_out(:,:,:) ! output wind (component) field

    ! Coordinate fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)  ! 3D height coordinate field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: coef1(:,:,:)    ! coefficient for linear term
    REAL(wp), INTENT(IN) :: coef2(:,:,:)    ! coefficient for quadratic term
    REAL(wp), INTENT(IN) :: coef3(:,:,:)    ! coefficient for cubic term
    REAL(wp), INTENT(IN) :: wfac_lin(:,:,:) ! weighting factor for linear interpolation

    INTEGER , INTENT(IN) :: idx0_cub(:,:,:) ! index of upper level, cubic interpolation
    INTEGER , INTENT(IN) :: idx0_lin(:,:,:) ! index of upper level, linear interpolation

    INTEGER , INTENT(IN) :: bot_idx_cub(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: bot_idx_lin(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    LOGICAL,  INTENT(IN) :: l_hires_intp ! mode for interpolation to (much) finer grid
    LOGICAL,  INTENT(IN), OPTIONAL :: lacc ! If true, use OpenACC
    LOGICAL,  INTENT(IN), OPTIONAL :: l_restore_fricred ! subtract/restore frictional reduction of wind speed
    REAL(wp), INTENT(IN), OPTIONAL :: extrap_limit  ! multiplicative limit in case of downward extrapolation

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, fricred, mult_limit


    REAL(wp), DIMENSION(nproma) :: uv1, uv2, dudz_up
    LOGICAL , DIMENSION(nproma) :: l_found

    LOGICAL :: lrestore_fricred

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, fric_red, uv_mod, g1, g2, g3
    REAL(wp), DIMENSION(nproma,nlevs_in-1) :: zalml_in_d
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out, zdiff_inout, red_speed
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0)))  RETURN

    ! consistency check:
    IF (UBOUND(uv_out,2) < nlevs_out) CALL finish("Wrong size of output field!")

    IF (PRESENT(l_restore_fricred)) THEN
      lrestore_fricred = l_restore_fricred
    ELSE
      lrestore_fricred = .TRUE.
    ENDIF

    ! Multiplicative limit for downward extrapolation; not used in high-resolution mode
    ! (wind speed in narrow valleys not resolved in the source model may be reduced
    ! more strongly)
    IF (PRESENT(extrap_limit)) THEN
      mult_limit = extrap_limit
    ELSE IF (lrestore_fricred) THEN
      mult_limit = 0.5_wp
    ELSE
      mult_limit = 1.0_wp
    ENDIF

!$OMP PARALLEL private(jc, jk, zalml_in,zalml_out)

    !$ACC DATA IF(lzacc) &
    !$ACC   PRESENT(uv_in, uv_out, z3d_in, z3d_out, coef1, coef2, coef3) &
    !$ACC   PRESENT(wfac_lin, idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin, kpbl1, wfacpbl1, kpbl2, wfacpbl2) &
    !$ACC   CREATE(ik1, uv1, uv2, dudz_up, l_found, zalml_in, fric_red, uv_mod, g1, g2, g3) &
    !$ACC   CREATE(zalml_out, zalml_in_d, zdiff_inout, red_speed)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO jk = 1, nlevs_in
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nproma
        zalml_in(jc, jk) = (jc - 1) * nlevs_in + jk
      END DO
    END DO
    !$ACC LOOP SEQ
    DO jk = 1, nlevs_out
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nproma
        zalml_out(jc, jk) = (jc - 1) * nlevs_out + jk
      END DO
    END DO
    !$ACC END PARALLEL

!-------------------------------------------------------------------------

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,fricred,&
!$OMP            uv1,uv2,dudz_up,zdiff_inout,red_speed,l_found,    &
!$OMP            uv_mod,zalml_in_d, fric_red,g1,g2,g3) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = nlen+1, nproma
          uv_out(jc,:,jb)  = 0.0_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nlen
        ! Wind component at height zpbl1
        uv1(jc) = wfacpbl1(jc,jb) *uv_in(jc,kpbl1(jc,jb),jb  ) + &
           (1._wp-wfacpbl1(jc,jb))*uv_in(jc,kpbl1(jc,jb)+1,jb)

        ! Wind component at height zpbl2
        uv2(jc) = wfacpbl2(jc,jb) *uv_in(jc,kpbl2(jc,jb),jb  ) + &
           (1._wp-wfacpbl2(jc,jb))*uv_in(jc,kpbl2(jc,jb)+1,jb)

        ! Vertical gradient between zpbl1 and zpbl2
        dudz_up(jc) = (uv2(jc) - uv1(jc))/(zpbl2 - zpbl1)
      ENDDO

      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen

          ! Height above lowest model level
          zalml_in(jc,jk1) = z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)

          ! Modified wind speed with frictional reduction removed
          IF (lrestore_fricred .AND. zalml_in(jc,jk1) < zpbl1) THEN
            uv_mod(jc,jk1) = uv1(jc)+dudz_up(jc)*(zalml_in(jc,jk1)-zpbl1)
          ELSE
            uv_mod(jc,jk1) = uv_in(jc,jk1,jb)
          ENDIF

          ! frictional reduction of wind speed
          fric_red(jc,jk1) = uv_mod(jc,jk1) - uv_in(jc,jk1,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      jk_start_in = start_idx_threshold(zpbl1, zalml_in(:,:), nlen, nlevs_in, lzacc)

      ! Compute factor for artificial wind speed reduction when interpolating
      ! into deep valleys, whose presence is inferred from zdiff_inout
      ! (of course, this inference assumes that the target grid has a much finer
      ! resolution than the source grid; thus, the reduction is used for l_hires_intp only

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      IF (l_hires_intp) THEN

        !$ACC LOOP SEQ
        DO jk = 1, nlevs_out
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen

            ! height distance between lowest input and output grid points
            ! (negative if extrapolation takes place)
            zdiff_inout(jc,jk) = z3d_out(jc,jk,jb) - z3d_in(jc,nlevs_in,jb)

            IF (zdiff_inout(jc,jk) <= -500._wp) THEN
              red_speed(jc,jk) = 0.25_wp
            ELSE IF (zdiff_inout(jc,jk) <= -100._wp) THEN
              red_speed(jc,jk) = 1._wp + (zdiff_inout(jc,jk) + 100._wp)*1.875e-3_wp
            ELSE
              red_speed(jc,jk) = 1._wp
            ENDIF

          ENDDO
        ENDDO

      ELSE !  no artificial reduction

        !$ACC LOOP SEQ
        DO jk = 1, nlevs_out
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            red_speed(jc,jk) = 1._wp
          ENDDO
        ENDDO

      ENDIF


      ! Compute vertical gradients of input data
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          g1(jc,jk1) = (uv_mod(jc,jk1 )-uv_mod(jc,jk1+1 ))/ &
                       (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
        ENDDO
      ENDDO

      ! Compute vertical gradients of gradients
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in-2
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          g2(jc,jk1) = (g1(jc,jk1)-g1(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+2,jb))
        ENDDO
      ENDDO

      ! Compute third-order vertical gradients
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in-3
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          g3(jc,jk1) = (g2(jc,jk1)-g2(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+3,jb))
        ENDDO
      ENDDO

      ! Now perform vertical interpolation, based on the modified temperature field
      !$ACC LOOP SEQ
      DO jk = 1, nlevs_out
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk1)
        DO jc = 1, nlen
          IF (jk <= bot_idx_cub(jc,jb)) THEN

            ! cubic interpolation
            jk1 = idx0_cub(jc,jk,jb)
            uv_out(jc,jk,jb) = uv_mod(jc,jk1)             + coef1(jc,jk,jb)*g1(jc,jk1) + &
                               coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

          ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

            ! linear interpolation
            jk1 = idx0_lin(jc,jk,jb)
            uv_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*uv_mod(jc,jk1) + &
              (1._wp-wfac_lin(jc,jk,jb))*uv_mod(jc,jk1+1)

          ELSE IF (lrestore_fricred) THEN

            ! downward extrapolation
            uv_out(jc,jk,jb) = uv_mod(jc,nlevs_in) + wfac_lin(jc,jk,jb)*dudz_up(jc)

            ! ensure that the extrapolated wind does not change sign and
            ! stays within the specified limit of speed reduction
            IF (uv_out(jc,jk,jb)*uv_mod(jc,nlevs_in) >= 0._wp) THEN ! still correct sign

              uv_out(jc,jk,jb) = SIGN(MAX(ABS(uv_out(jc,jk,jb)),  &
                ABS(mult_limit*uv_mod(jc,nlevs_in))), uv_mod(jc,nlevs_in))

              ! For extrapolation to z-levels or p-levels, this limitation is needed in addition
              uv_out(jc,jk,jb) = SIGN(MIN(ABS(uv_out(jc,jk,jb)),  &
                ABS(1.25_wp*uv_mod(jc,nlevs_in))), uv_mod(jc,nlevs_in))

            ELSE  ! i.e. spurious sign change due to extrapolation

              uv_out(jc,jk,jb) = mult_limit*uv_mod(jc,nlevs_in)

            ENDIF

            ! Artificial reduction for large downward extrapolation distances
            ! if high-resolution mode is selected
            uv_out(jc,jk,jb) = uv_out(jc,jk,jb)*red_speed(jc,jk)

          ELSE ! in post-processing mode, just use constant profile below ground

            uv_out(jc,jk,jb) = uv_mod(jc,nlevs_in)

          ENDIF

          ! Height above lowest model level - needed for restoring the frictional layer
          zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! subtract frictional reduction from preliminary wind field
      IF (lrestore_fricred) THEN

        jk_start_out = start_idx_threshold(zpbl1, zalml_out(:,:), nlen, nlevs_out, lzacc)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlevs_in-1
          DO jc = 1, nlen
            zalml_in_d(jc,jk) = 1.0_wp / (zalml_in(jc,jk) - zalml_in(jc,jk+1))
          ENDDO
        ENDDO
        !$ACC END PARALLEL

        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            l_found(jc) = .FALSE.
          ENDDO
          !$ACC LOOP SEQ
          DO jk1 = jk_start, nlevs_in-1
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wfac, fricred)
            DO jc = 1, nlen
              IF (.NOT. l_found(jc) ) THEN
                IF (zalml_out(jc,jk) >= zpbl1) THEN
                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk_start
                ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                         zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                  wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1)) * zalml_in_d(jc,jk1)
                  fricred = wfac*fric_red(jc,jk1) + (1._wp-wfac)*fric_red(jc,jk1+1)

                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk1

                  uv_out(jc,jk,jb) = (uv_out(jc,jk,jb)/red_speed(jc,jk)-fricred)*red_speed(jc,jk)

                ELSE IF (zalml_out(jc,jk) > zalml_in(jc,jk_start)) THEN
                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk_start
                ENDIF
              ENDIF
            ENDDO
#ifdef _OPENACC
            ! ACC: the following EXIT would be illegal within an OpenACC Kernel, thus we skip this CPU optimization
            IF (ALL(l_found(1:nlen))) EXIT
#endif
          ENDDO
          !$ACC END PARALLEL

          jk_start = minval_1d(ik1(1:nlen), lacc=lzacc)
        ENDDO
      ENDIF

    ENDDO
!$OMP END DO NOWAIT
    !$ACC WAIT
    !$ACC END DATA
!$OMP END PARALLEL

  END SUBROUTINE uv_intp


  !-------------
  !>
  !! SUBROUTINE qv_intp
  !! Performs vertical interpolation of specific humidity,
  !! including checks for supersaturation and unrealistically small values
  !!
  !! Required input fields: 3D input field to be interpolated,
  !! coefficient fields from prepare_lin/cubic_intp and prepare_extrap
  !! Output: interpolated 3D field
  !!
  !! Performs cubic interpolation where possible, turning to linear interpolation
  !! close to the surface
  !!
  SUBROUTINE qv_intp(qv_in, qv_out, z3d_in, z3d_out,                  &
                     temp_in, pres_in, temp_out, pres_out,            &
                     nblks, npromz, nlevs_in, nlevs_out,              &
                     coef1, coef2, coef3, wfac_lin,                   &
                     idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin,    &
                     wfacpbl1, kpbl1, wfacpbl2, kpbl2,                &
                     lower_limit, l_satlimit, l_restore_pbldev,       &
                     opt_hires_corr, opt_qc, opt_lmask, lacc)


    ! Specific humidity fields
    REAL(wp), INTENT(IN)  :: qv_in (:,:,:) ! input field
    REAL(wp), INTENT(OUT) :: qv_out(:,:,:) ! output field

    ! Optional cloud water field: if provided, consistency checks between QV and QC will
    ! be performed at the end of this routine
    REAL(wp), INTENT(INOUT), OPTIONAL :: opt_qc(:,:,:) ! specific cloud water

    ! Additional atmospheric fields on input and output grids (needed for supersaturation check)
    ! Note: these are all input data
    REAL(wp), INTENT(IN) :: temp_in(:,:,:)  ! temperature on input grid
    REAL(wp), INTENT(IN) :: pres_in(:,:,:)  ! pressure on input grid
    REAL(wp), INTENT(IN) :: temp_out(:,:,:) ! temperature on output grid
    REAL(wp), INTENT(IN) :: pres_out(:,:,:) ! pressure on output grid

    ! Coordinate fields
    REAL(wp), INTENT(IN) :: z3d_in(:,:,:)   ! 3D height coordinate field of input data
    REAL(wp), INTENT(IN) :: z3d_out(:,:,:)  ! 3D height coordinate field of output data

    ! Dimension parameters
    INTEGER , INTENT(IN) :: nblks      ! Number of blocks
    INTEGER , INTENT(IN) :: npromz     ! Length of last block
    INTEGER , INTENT(IN) :: nlevs_in   ! Number of input levels
    INTEGER , INTENT(IN) :: nlevs_out  ! Number of output levels

    ! Coefficients
    REAL(wp), INTENT(IN) :: coef1(:,:,:)    ! coefficient for linear term
    REAL(wp), INTENT(IN) :: coef2(:,:,:)    ! coefficient for quadratic term
    REAL(wp), INTENT(IN) :: coef3(:,:,:)    ! coefficient for cubic term
    REAL(wp), INTENT(IN) :: wfac_lin(:,:,:) ! weighting factor for linear interpolation

    INTEGER , INTENT(IN) :: idx0_cub(:,:,:) ! index of upper level, cubic interpolation
    INTEGER , INTENT(IN) :: idx0_lin(:,:,:) ! index of upper level, linear interpolation

    INTEGER , INTENT(IN) :: bot_idx_cub(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: bot_idx_lin(:,:)! index of lowest level for which cubic interpolation is possible
    INTEGER , INTENT(IN) :: kpbl1(:,:)      ! index of model level immediately above (by default) 500 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl1(:,:)   ! corresponding interpolation coefficient
    INTEGER , INTENT(IN) :: kpbl2(:,:)      ! index of model level immediately above (by default) 1000 m AGL
    REAL(wp), INTENT(IN) :: wfacpbl2(:,:)   ! corresponding interpolation coefficient

    REAL(wp), INTENT(IN) :: lower_limit     ! lower limit of QV
    LOGICAL , INTENT(IN) :: l_satlimit       ! limit input field to water saturation
    LOGICAL , INTENT(IN) :: l_restore_pbldev ! restore PBL deviation of QV from extrapolated profile
    LOGICAL , INTENT(IN), OPTIONAL :: lacc ! If true, use OpenACC

    LOGICAL, OPTIONAL, INTENT(IN) :: opt_hires_corr   ! apply corrections / limits for coarse-to-fine grid interpolation
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_lmask(:,:)

    ! LOCAL VARIABLES

    INTEGER  :: jb, jk, jk1, jc, nlen, jk_start, jk_start_in, jk_start_out, ik1(nproma)
    REAL(wp) :: wfac, pbldev, rhum, qtot
    REAL(wp) :: e_vapor                     ! vapor pressure


    REAL(wp), DIMENSION(nproma) :: qv1, qv2, dqvdz_up
    LOGICAL , DIMENSION(nproma) :: l_found
    LOGICAL                     :: l_check_qv_qc, lhr_corr

    REAL(wp), DIMENSION(nproma,nlevs_in)  :: zalml_in, pbl_dev, qv_mod, g1, g2, g3, qsat_in, qv_in_lim
    REAL(wp), DIMENSION(nproma,nlevs_in-1) :: zalml_in_d
    REAL(wp), DIMENSION(nproma,nlevs_out) :: zalml_out, qsat_out
    LOGICAL :: lmask(nproma)
    LOGICAL :: lzacc ! non-optional version of lacc

!-------------------------------------------------------------------------
    CALL set_acc_host_or_device(lzacc, lacc)

    ! return, if nothing to do:
    IF ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0))) RETURN

    IF (PRESENT(opt_qc)) THEN
      l_check_qv_qc = .TRUE.
    ELSE
      l_check_qv_qc = .FALSE.
    ENDIF

    IF (PRESENT(opt_hires_corr)) THEN
      lhr_corr = opt_hires_corr
    ELSE
      lhr_corr = .FALSE.
    ENDIF

!$OMP PARALLEL private(jc, jk, zalml_in, zalml_out)

    !$ACC DATA IF(lzacc) &
    !$ACC   PRESENT(qv_in, qv_out, opt_qc, temp_in, pres_in, temp_out, pres_out, z3d_in, z3d_out, coef1, coef2) &
    !$ACC   PRESENT(coef3, wfac_lin, idx0_cub, idx0_lin, bot_idx_cub, bot_idx_lin, kpbl1, wfacpbl1, kpbl2) &
    !$ACC   PRESENT(wfacpbl2) &
    !$ACC   CREATE(qv1, qv2, dqvdz_up) &
    !$ACC   CREATE(l_found, lmask, ik1) &
    !$ACC   CREATE(zalml_in, pbl_dev, qv_mod, g1, g2, g3, qsat_in, qv_in_lim) &
    !$ACC   CREATE(zalml_out, qsat_out, zalml_in_d)

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO jk = 1, nlevs_in
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nproma
        zalml_in(jc, jk) = (jc - 1) * nlevs_in + jk
      END DO
    END DO
    !$ACC LOOP SEQ
    DO jk = 1, nlevs_out
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nproma
        zalml_out(jc, jk) = (jc - 1) * nlevs_out + jk
      END DO
    END DO
    !$ACC END PARALLEL

!-------------------------------------------------------------------------

!$OMP DO PRIVATE(jb,jk,jk1,jc,nlen,jk_start,jk_start_in,jk_start_out,ik1,wfac,pbldev,&
!$OMP            rhum,qtot,qv1,qv2,dqvdz_up,l_found,zalml_in_d,qv_in_lim,    &
!$OMP            qv_mod,pbl_dev,g1,g2,g3,e_vapor,qsat_in,qsat_out,lmask) ICON_OMP_DEFAULT_SCHEDULE

    DO jb = 1, nblks
      IF (jb /= nblks) THEN
        nlen = nproma
      ELSE
        nlen = npromz
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = nlen+1, nproma
          qv_out(jc,:,jb)  = 0.0_wp
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDIF

      IF (PRESENT(opt_lmask)) THEN
        CALL assert_acc_host_only("qv_intp:opt_lmask", lacc)
        lmask(:) =  opt_lmask(:,jb)
        IF (.NOT. ANY(lmask)) THEN
          qv_out(:,:,jb)  = 0.0_wp
          CYCLE
        ENDIF
      ELSE
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = 1, nproma
          lmask(jc) = .TRUE.
        ENDDO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(e_vapor)
        DO jc = 1, nlen
          IF (lmask(jc)) THEN

            ! Limit input data to the 'lower_limit' specified in the argument list
            qv_in_lim(jc,jk1) = MAX(qv_in(jc,jk1,jb),lower_limit)

            ! compute vapour pressure e_vapor=f(qv_in,pres_in)
            e_vapor = vap_pres(qv_in_lim(jc,jk1),pres_in(jc,jk1,jb))

            ! saturation specific humidity of input data
            qsat_in(jc,jk1) = rdv*sat_pres_water(temp_in(jc,jk1,jb)) /    &
                              (pres_in(jc,jk1,jb)-o_m_rdv*e_vapor)

            IF (l_satlimit) THEN
              ! limit input data to water saturation when processing interpolated IFS data:
              ! This is needed to remove supersaturations generated (primarily) by the interpolation
              ! from the spherical harmonics to the Gaussain grid; without this limitation, the
              ! QV-QC-adjustment at the end of this routine would generate nonsensically large
              ! cloud water peaks
              qv_in_lim(jc,jk1) = MIN(qv_in_lim(jc,jk1),qsat_in(jc,jk1))
            ENDIF
          ELSE
            qv_in_lim(jc,jk1) = MAX(qv_in(jc,jk1,jb),lower_limit)
            qsat_in(jc,jk1)   = qv_in_lim(jc,jk1)
          ENDIF
        ENDDO
      ENDDO



      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = 1, nlen
        ! QV at height zpbl1
        qv1(jc) = wfacpbl1(jc,jb) *qv_in_lim(jc,kpbl1(jc,jb)  ) + &
           (1._wp-wfacpbl1(jc,jb))*qv_in_lim(jc,kpbl1(jc,jb)+1)

        ! QV at height zpbl2
        qv2(jc) = wfacpbl2(jc,jb) *qv_in_lim(jc,kpbl2(jc,jb)  ) + &
           (1._wp-wfacpbl2(jc,jb))*qv_in_lim(jc,kpbl2(jc,jb)+1)

        ! Vertical gradient between zpbl1 and zpbl2
        dqvdz_up(jc) = MIN(5.e-6_wp, MAX(-5.e-6_wp, (qv2(jc) - qv1(jc))/(zpbl2 - zpbl1)))
      ENDDO

      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen

          ! Height above lowest model level
          zalml_in(jc,jk1) = z3d_in(jc,jk1,jb) - z3d_in(jc,nlevs_in,jb)

          ! Modified QV with boundary-layer deviation from the extrapolated profile removed
          IF (zalml_in(jc,jk1) < zpbl1) THEN
            qv_mod(jc,jk1) = qv1(jc)+dqvdz_up(jc)*(zalml_in(jc,jk1)-zpbl1)
            qv_mod(jc,jk1) = MIN(qv_mod(jc,jk1),1.25_wp*qv_in_lim(jc,jk1),qsat_in(jc,jk1))
            qv_mod(jc,jk1) = MAX(qv_mod(jc,jk1),0.75_wp*qv_in_lim(jc,jk1))
          ELSE
            qv_mod(jc,jk1) = qv_in_lim(jc,jk1)
          ENDIF

          ! boundary-layer deviation of QV, converted to RH
          pbl_dev(jc,jk1) = (qv_mod(jc,jk1) - qv_in_lim(jc,jk1)) / qsat_in(jc,jk1)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      jk_start_in = start_idx_threshold(zpbl1, zalml_in(:,:), nlen, nlevs_in, lzacc)

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      ! Compute vertical gradients of input data
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          g1(jc,jk1) = (qv_mod(jc,jk1 )-qv_mod(jc,jk1+1 ))/ &
                       (z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+1,jb))
        ENDDO
      ENDDO

      ! Compute vertical gradients of gradients
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in-2
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          g2(jc,jk1) = (g1(jc,jk1)-g1(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+2,jb))
        ENDDO
      ENDDO

      ! Compute third-order vertical gradients
      !$ACC LOOP SEQ
      DO jk1 = 1, nlevs_in-3
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nlen
          g3(jc,jk1) = (g2(jc,jk1)-g2(jc,jk1+1))/(z3d_in(jc,jk1,jb)-z3d_in(jc,jk1+3,jb))
        ENDDO
      ENDDO

      ! Now perform vertical interpolation, based on the modified temperature field
      !$ACC LOOP SEQ
      DO jk = 1, nlevs_out
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(jk1, e_vapor)
        DO jc = 1, nlen

          IF (jk <= bot_idx_cub(jc,jb)) THEN

            ! cubic interpolation
            jk1 = idx0_cub(jc,jk,jb)
            qv_out(jc,jk,jb) = qv_mod(jc,jk1)             + coef1(jc,jk,jb)*g1(jc,jk1) + &
                               coef2(jc,jk,jb)*g2(jc,jk1) + coef3(jc,jk,jb)*g3(jc,jk1)

          ELSE IF (jk <= bot_idx_lin(jc,jb)) THEN

            ! linear interpolation
            jk1 = idx0_lin(jc,jk,jb)
            qv_out(jc,jk,jb) = wfac_lin(jc,jk,jb)*qv_mod(jc,jk1) + &
              (1._wp-wfac_lin(jc,jk,jb))*qv_mod(jc,jk1+1)

          ELSE

            ! if extrapolation is needed, maintain relative humidity (with a slight
            ! approximation of the qsat expression to avoid iterations)
            qv_out(jc,jk,jb) = qv_mod(jc,nlevs_in)/qsat_in(jc,nlevs_in)* &
              rdv*sat_pres_water(temp_out(jc,jk,jb))/pres_out(jc,jk,jb)

            IF (lhr_corr) THEN ! prevent excessively high humidities in valleys not resolved in the source model
              qv_out(jc,jk,jb) = MIN(qv_out(jc,jk,jb),                                             &
                qv_mod(jc,nlevs_in) * (1._wp+4.e-4_wp*(z3d_in(jc,nlevs_in,jb)-z3d_out(jc,jk,jb))), &
                qv_mod(jc,nlevs_in) + 2.5e-6_wp*(z3d_in(jc,nlevs_in,jb)-z3d_out(jc,jk,jb))          )
            ENDIF

          ENDIF

          ! Height above lowest model level - needed for restoring the boundary-layer deviation
          zalml_out(jc,jk) = z3d_out(jc,jk,jb) - z3d_out(jc,nlevs_out,jb)

          ! compute vapour pressure e_vapor=f(qv_out,pres_out)
          e_vapor = vap_pres(qv_out(jc,jk,jb),pres_out(jc,jk,jb))

          ! saturation specific humidity of output data
          qsat_out(jc,jk) = rdv*sat_pres_water(temp_out(jc,jk,jb)) /    &
                            (pres_out(jc,jk,jb)-o_m_rdv*e_vapor)

        ENDDO
      ENDDO
      !$ACC END PARALLEL

      jk_start_out = start_idx_threshold(zpbl1, zalml_out(:,:), nlen, nlevs_out, lzacc)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlevs_in-1
        DO jc = 1, nlen
          zalml_in_d(jc,jk) = 1.0_wp / (zalml_in(jc,jk) - zalml_in(jc,jk+1))
        ENDDO
      ENDDO
      !$ACC END PARALLEL

      ! subtract boundary-layer deviation from preliminary QV
      IF (l_restore_pbldev) THEN
        jk_start = jk_start_in - 1
        DO jk = jk_start_out, nlevs_out
          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO jc = 1, nlen
            l_found(jc) = .FALSE.
          ENDDO
          !$ACC LOOP SEQ
          DO jk1 = jk_start, nlevs_in-1
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO jc = 1, nlen
              IF( .NOT. l_found(jc) ) THEN
                IF (zalml_out(jc,jk) >= zpbl1) THEN
                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk_start
                ELSE IF (zalml_out(jc,jk) <  zalml_in(jc,jk1) .AND. &
                         zalml_out(jc,jk) >= zalml_in(jc,jk1+1)) THEN

                  wfac = (zalml_out(jc,jk)-zalml_in(jc,jk1+1)) * zalml_in_d(jc,jk1)
                  pbldev = wfac*pbl_dev(jc,jk1) + (1._wp-wfac)*pbl_dev(jc,jk1+1)

                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk1

                  qv_out(jc,jk,jb) =  qv_out(jc,jk,jb) - pbldev*qsat_out(jc,jk)

                ELSE IF (zalml_out(jc,jk) > zalml_in(jc,jk_start)) THEN
                  l_found(jc) = .TRUE.
                  ik1(jc)     = jk_start
                ENDIF
              ENDIF
            ENDDO
#ifdef _OPENACC
            ! ACC: the following EXIT would be illegal within an OpenACC Kernel, thus we skip this CPU optimization
            IF (ALL(l_found(1:nlen))) EXIT
#endif
          ENDDO
          !$ACC END PARALLEL
          jk_start = minval_1d(ik1(1:nlen), lacc=lzacc)

        ENDDO
      ENDIF

      IF (l_check_qv_qc) THEN ! apply consistency checks between QV and QC
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(rhum, qtot)
        DO jk = 1, nlevs_out
          DO jc = 1, nlen

            ! Impose limiter to avoid negative or unreasonably small values,
            ! and ensure that RH > 1% below 200 hPa
            qv_out(jc,jk,jb) = MAX(lower_limit,qv_out(jc,jk,jb))

            IF (pres_out(jc,jk,jb) >= 20000._wp) THEN
              qv_out(jc,jk,jb) = MAX(qv_out(jc,jk,jb),0.01_wp*qsat_out(jc,jk))
            ENDIF

            rhum = qv_out(jc,jk,jb)/qsat_out(jc,jk)

            IF (rhum < 0.75_wp) THEN ! remove all the cloud water if present
              opt_qc(jc,jk,jb) = 0._wp
            ELSE ! conserve QV plus QC
              qtot = qv_out(jc,jk,jb) + opt_qc(jc,jk,jb)
              opt_qc(jc,jk,jb) = MAX(0._wp,qtot-qsat_out(jc,jk))
              qv_out(jc,jk,jb) = qtot - opt_qc(jc,jk,jb)
            ENDIF

          ENDDO
        ENDDO
        !$ACC END PARALLEL

      ELSE ! impose only some obvious constraints on QV
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO jk = 1, nlevs_out
          DO jc = 1, nlen

            ! Impose limiter to avoid negative or unreasonably small values
            qv_out(jc,jk,jb) = MAX(lower_limit,qv_out(jc,jk,jb))

            ! in addition, limit output data to water saturation,
            ! ensure that RH > 1% below 200 hPa
            qv_out(jc,jk,jb) = MIN(qv_out(jc,jk,jb),qsat_out(jc,jk))

            IF (pres_out(jc,jk,jb) >= 20000._wp) THEN
              qv_out(jc,jk,jb) = MAX(qv_out(jc,jk,jb),0.01_wp*qsat_out(jc,jk))
            ENDIF

          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDIF

    ENDDO
    !$ACC WAIT
    !$ACC END DATA

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE qv_intp

END MODULE mo_nh_vert_interp
