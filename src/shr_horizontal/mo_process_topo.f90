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

! Contains routine for topography smoothing

MODULE mo_process_topo

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text
  USE mo_model_domain,       ONLY: t_patch
  USE mo_parallel_config,    ONLY: nproma
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_sync,               ONLY: SYNC_C, sync_patch_array
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_extpar_config,      ONLY: fac_smooth_topo, n_iter_smooth_topo,           &
    &                              heightdiff_threshold, hgtdiff_max_smooth_topo, &
    &                              lrevert_sea_height, pp_sso
  USE mo_impl_constants,     ONLY: min_rlcell, min_rlcell_int
  USE mo_math_constants,     ONLY: rad2deg
  USE mo_math_laplace,       ONLY: nabla2_scalar, nabla4_scalar
  USE mo_intp_rbf,           ONLY: rbf_interpol_c2grad

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: compute_smooth_topo
  PUBLIC :: smooth_topo
  PUBLIC :: smooth_topo_real_data, postproc_sso
  PUBLIC :: smooth_frland, smooth_urbfrac

CONTAINS

  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !! Computes the smoothed topography needed e.g. for the SLEVE coordinate.
  !! May be bypassed once an option for reading the smooth topography from data
  !! is available
  !!
  SUBROUTINE compute_smooth_topo(p_patch, p_int, topo_c, niter, topo_smt_c)

    TYPE(t_patch),TARGET,INTENT(IN) :: p_patch
    TYPE(t_int_state), INTENT(IN) :: p_int

    ! Input fields: topography on cells
    REAL(wp), INTENT(IN) :: topo_c(:,:)

    ! number of iterations
    INTEGER,  INTENT(IN) :: niter

    ! Output fields: smooth topography on cells
    REAL(wp), INTENT(OUT) :: topo_smt_c(:,:)

    INTEGER  :: jb, jc, iter
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),nabla2_topo(nproma,1,p_patch%nblks_c)

    !-------------------------------------------------------------------------

    ! Initialize auxiliary fields for topography with data and nullify nabla2 field
    z_topo(:,1,:)      = topo_c(:,:)
    nabla2_topo(:,1,:) = 0._wp

    i_startblk = p_patch%cells%start_blk(2,1)
    nblks_c    = p_patch%nblks_c

    CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ! Apply nabla2-diffusion niter times to create smooth topography
    DO iter = 1, niter

      CALL nabla2_scalar(z_topo, p_patch, p_int, nabla2_topo, &
        &                 slev=1, elev=1, rl_start=2, rl_end=min_rlcell )

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 2)

        DO jc = i_startidx, i_endidx
          z_topo(jc,1,jb) = z_topo(jc,1,jb) + 0.125_wp*nabla2_topo(jc,1,jb) &
            &                               * p_patch%cells%area(jc,jb)
        ENDDO
      ENDDO

      CALL sync_patch_array(SYNC_C,p_patch,z_topo)

    ENDDO

    ! Store smooth topography on output fields
    topo_smt_c(:,:) = z_topo(:,1,:)

  END SUBROUTINE compute_smooth_topo



  !-----------------------------------------------------------------------
  !! Topography smoothing
  !!
  !! Topography smoothing by selectively applying nabla2 and nabla4 
  !! operators.
  !!
  SUBROUTINE smooth_topo (p_patch, p_int, topography_c)

    TYPE(t_patch)      , INTENT(IN)    :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(INOUT) :: topography_c(:,:) ! original topography on input
                                                            ! smoothed one on output

    ! local variables
    INTEGER  :: jg, jb, jc, iter, il, npts, niter, iter_sub
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo(nproma,1,p_patch%nblks_c),z_nabla4_topo(nproma,1,p_patch%nblks_c),    &
      &         z_topo_old(nproma,1,p_patch%nblks_c),z_nabla2_topo(nproma,1,p_patch%nblks_c),&
      &         z_hdiffmax(nproma,p_patch%nblks_c)
    REAL(wp) :: zmaxtop,zmintop,z_topo_new,zdcoeff,z_heightdiff_threshold,rms_hdiff,smooth_fac
    LOGICAL  :: lnabla2_mask(nproma,p_patch%nblks_c)


    jg = p_patch%id

    z_nabla4_topo(:,1,:) = 0._wp
    z_nabla2_topo(:,1,:) = 0._wp

    nblks_c    = p_patch%nblks_c

    zdcoeff = 0.05_wp ! diffusion coefficient for nabla2 diffusion
    niter   = 20      ! number of iterations for local nabla2 diffusion

    IF (n_iter_smooth_topo(jg) > 0) THEN
      WRITE(message_text,'(a,i3,a,i3)') 'number of topography smoothing steps in domain ', &
        jg, ': ', n_iter_smooth_topo(jg)
      CALL message('', TRIM(message_text))
    ENDIF

    z_topo(:,1,:)   = topography_c(:,:)

    DO iter = 1, n_iter_smooth_topo(jg)

      ! Step 1: local nabla2 diffusion at grid points where the height difference to the
      ! neighbors exceeds a certain threshold value

      DO iter_sub = 1, niter ! perform niter iterations
                             ! note: a variable number of iterations (with an exit condition) potentially
                             ! causes trouble with MPI reproducibility

        z_heightdiff_threshold = heightdiff_threshold(jg)
        IF (iter_sub >= niter-1 .AND. heightdiff_threshold(jg) < 2500._wp) &
          z_heightdiff_threshold = 0.75_wp*heightdiff_threshold(jg)

        z_hdiffmax(:,:)   = 0._wp
        lnabla2_mask(:,:) = .FALSE.

        CALL nabla2_scalar( z_topo, p_patch, p_int, z_nabla2_topo, &
          &                 slev=1, elev=1, rl_start=2, rl_end=min_rlcell )

        npts = 0

        i_startblk = p_patch%cells%start_blk(2,1)

        DO jb = i_startblk,nblks_c

          CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                             i_startidx, i_endidx, 2)

          DO il=1,p_patch%geometry_info%cell_type
            DO jc = i_startidx, i_endidx
              z_hdiffmax(jc,jb) = MAX(z_hdiffmax(jc,jb), ABS(z_topo(jc,1,jb) - &
                &  z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1,              &
                &         p_patch%cells%neighbor_blk(jc,jb,il))) )

              IF (z_hdiffmax(jc,jb) > z_heightdiff_threshold) lnabla2_mask(jc,jb) = .TRUE.

            ENDDO
          ENDDO

        ENDDO

        ! set diffusion mask also true if one of the neighboring grid points 
        ! fulfills the height difference criterion
        i_startblk = p_patch%cells%start_blk(3,1)

        DO jb = i_startblk,nblks_c

          CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                             i_startidx, i_endidx, 3)

          DO il=1,p_patch%geometry_info%cell_type
            DO jc = i_startidx, i_endidx
              IF (z_hdiffmax(p_patch%cells%neighbor_idx(jc,jb,il),                   &
                  p_patch%cells%neighbor_blk(jc,jb,il)) > z_heightdiff_threshold) THEN
                lnabla2_mask(jc,jb) = .TRUE.
              ENDIF
            ENDDO
          ENDDO

          DO jc = i_startidx, i_endidx
            IF (lnabla2_mask(jc,jb)) THEN
              npts = npts + 1
              topography_c(jc,jb) = z_topo(jc,1,jb) + zdcoeff *                      &
                                    p_patch%cells%area(jc,jb) * z_nabla2_topo(jc,1,jb)
            ENDIF
          ENDDO

        ENDDO

        z_topo(:,1,:) = topography_c(:,:)
        CALL sync_patch_array(SYNC_C, p_patch, z_topo)
        topography_c(:,:) = z_topo(:,1,:)

      ENDDO ! iteration of local nabla2 smoothing

      i_startblk = p_patch%cells%start_blk(3,1)

      ! Step 2: local nabla4 diffusion with monotonous limiter
      z_topo(:,1,:)   = topography_c(:,:)
      z_topo_old(:,1,:) = z_topo(:,1,:)

      CALL nabla4_scalar(z_topo, p_patch, p_int, z_nabla4_topo, &
        & slev=1, elev=UBOUND(z_topo,2), rl_start=3, rl_end=min_rlcell, &
        & p_nabla2=z_nabla2_topo )

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema

          ! compute maximum (zmaxtop) and minimum (zmintop) of neighbor cells topography
          !set zmaxtop and zmintop to first neighbor's value
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          rms_hdiff = (z_topo(jc,1,jb)-z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
                      p_patch%cells%neighbor_blk(jc,jb,1)))**2
          DO il=2,p_patch%geometry_info%cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            !zmaxtop and zmintop are now max resp min of all neighbors

            rms_hdiff = rms_hdiff + (z_topo(jc,1,jb)-z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                      p_patch%cells%neighbor_blk(jc,jb,il)))**2
          ENDDO

          rms_hdiff = SQRT(rms_hdiff/REAL(p_patch%geometry_info%cell_type,wp))

          IF (hgtdiff_max_smooth_topo(jg) < 1._wp) THEN
            smooth_fac = fac_smooth_topo
          ELSE
            smooth_fac = fac_smooth_topo * MIN(1._wp,rms_hdiff/hgtdiff_max_smooth_topo(jg))
          ENDIF

          z_topo_new = z_topo (jc,1,jb) -  &
              smooth_fac * z_nabla4_topo(jc,1,jb) * &
              p_patch%cells%area(jc,jb)*p_patch%cells%area(jc,jb)

          !If it was a local maximum in the old field, dont make it higher
          IF ( zmaxtop < z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) CYCLE
          ENDIF
          !If it was a local minimum in the old field, dont make it lower:
          IF ( zmintop > z_topo (jc,1,jb) ) THEN
             IF ( z_nabla4_topo(jc,1,jb) > 0.0_wp ) CYCLE
          ENDIF

          !If it became a local maximum in the new field with regard to old neighbors, avoid it:
          IF (( zmaxtop < z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp )) THEN
            topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
            CYCLE
          ENDIF
          !If it became a local minimum in the new field with regard to old neighbors, avoid it:
          IF (( zmintop > z_topo_new  ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
            topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
            CYCLE
          ENDIF

          topography_c(jc,jb)=z_topo_new

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)
        DO jc = i_startidx, i_endidx

          !Limiter to avoid amplification of local extrema
          zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,1),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,1))
          zmintop = zmaxtop
          DO il=2,p_patch%geometry_info%cell_type

            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) > &
              & zmaxtop ) THEN
              zmaxtop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
            IF ( z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
              &          p_patch%cells%neighbor_blk(jc,jb,il)) < &
              & zmintop ) THEN
              zmintop = z_topo(p_patch%cells%neighbor_idx(jc,jb,il),1, &
                &          p_patch%cells%neighbor_blk(jc,jb,il) )
            ENDIF
           !zmaxtop and zmintop are now mx resp min of all neighbors
          ENDDO

           !If it became a local maximum in the new field, avoid it
           IF ( ( zmaxtop < z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) < 0.0_wp ) ) THEN
             topography_c(jc,jb) = MAX(z_topo_old(jc,1,jb),zmaxtop)
           !If it became a local minimum in the new field, avoid it
           ELSEIF (( zmintop > z_topo(jc,1,jb) ) .AND. ( z_nabla4_topo(jc,1,jb) > 0.0_wp )) THEN
             topography_c(jc,jb) = MIN(z_topo_old(jc,1,jb),zmintop)
           ENDIF

        ENDDO

      ENDDO

      z_topo(:,1,:)   = topography_c(:,:)
      CALL sync_patch_array(SYNC_C, p_patch, z_topo)
      topography_c(:,:)=z_topo(:,1,:)

    ENDDO !iter

  END SUBROUTINE smooth_topo



  !-----------------------------------------------------------------------
  !!
  !! Topography smoothing for real-case runs. Apart from smoothing the 
  !! topography field, 
  !! - the height of sea-points is reset to that of the raw topography 
  !!   data set
  !! - the SSO standard deviation field is updated based on the smoothed 
  !!   topography
  !!
  SUBROUTINE smooth_topo_real_data (p_patch, p_int, fr_land, topography_c, fr_lake, sso_stdh)

    TYPE(t_patch)      , INTENT(IN)    :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(IN)    :: fr_land(:,:)
    REAL(wp)           , INTENT(INOUT) :: topography_c(:,:)
    REAL(wp), OPTIONAL , INTENT(IN)    :: fr_lake(:,:)
    REAL(wp), OPTIONAL , INTENT(INOUT) :: sso_stdh(:,:)

    ! local variables
    INTEGER  :: jb, jc
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: z_topo_c_sv(nproma,p_patch%nblks_c)
    REAL(wp) :: zhdiff
    LOGICAL  :: is_present_fr_lake

    nblks_c    = p_patch%nblks_c

    i_startblk = p_patch%cells%start_blk(3,1)

    ! save original raw topography
    z_topo_c_sv(:,:) = topography_c(:,:)

    ! topography smoothing
    !
    CALL smooth_topo(p_patch, p_int, topography_c)

    ! bring sea-points back to zero height (i.e. to original extpar values)
    !
    IF (lrevert_sea_height) THEN

      is_present_fr_lake = PRESENT(fr_lake)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        DO jc=i_startidx, i_endidx
          !
          ! bring grid cell back to zero height, if it is entirely covered by sea.
          IF (is_present_fr_lake) THEN
            ! NPW physics
            IF ((fr_land(jc,jb) + fr_lake(jc,jb)) == 0._wp) topography_c(jc,jb) = z_topo_c_sv(jc,jb)
          ELSE
            ! AES physics
            IF (fr_land(jc,jb) < EPSILON(1._wp)) topography_c(jc,jb) = z_topo_c_sv(jc,jb)
          ENDIF
        ENDDO  !jc

      ENDDO  !jb
    ENDIF

    ! re-compute SSO_STDH based on smoothed topography
    !
    IF (PRESENT(sso_stdh) .AND. pp_sso <= 1) THEN

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 3)

        DO jc = i_startidx, i_endidx

          zhdiff = topography_c(jc,jb) - z_topo_c_sv(jc,jb)
          sso_stdh(jc,jb) = SQRT(sso_stdh(jc,jb)**2 + zhdiff**2)

        ENDDO
      ENDDO

    ENDIF


  END SUBROUTINE smooth_topo_real_data

  !-----------------------------------------------------------------------
  !!
  !! Smoothed land fraction for adaptive tuning of sea ice bottom heat flux
  !!
  SUBROUTINE smooth_frland (p_patch, p_int, fr_land, fr_land_smt)

    TYPE(t_patch)      , INTENT(IN)    :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(IN)    :: fr_land(:,:)
    REAL(wp)           , INTENT(OUT)   :: fr_land_smt(:,:)

    ! local variables
    INTEGER  :: jb, jc, jb1, jc1, iter, j
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx
    REAL(wp) :: zaux(nproma,p_patch%nblks_c)

    nblks_c    = p_patch%nblks_c

    i_startblk = p_patch%cells%start_blk(3,1)

    fr_land_smt(:,:) = fr_land(:,:)
    zaux(:,:)        = fr_land(:,:)

    DO iter = 1, NINT(200.e3_wp/p_patch%geometry_info%mean_characteristic_length)

      DO jb = i_startblk,nblks_c

        CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, &
                           i_startidx, i_endidx, 3)

        zaux(i_startidx:i_endidx,jb) = 0._wp

        DO j=1, p_int%cell_environ%max_nmbr_nghbr_cells
          DO jc = i_startidx, i_endidx
            IF (j <= p_int%cell_environ%nmbr_nghbr_cells(jc,jb)) THEN
              jc1 = p_int%cell_environ%idx(jc,jb,j)
              jb1 = p_int%cell_environ%blk(jc,jb,j)
              zaux(jc,jb) = zaux(jc,jb) + fr_land_smt(jc1,jb1)/REAL(p_int%cell_environ%nmbr_nghbr_cells(jc,jb),wp)
            ENDIF
          ENDDO
        ENDDO

      ENDDO  !jb

      CALL sync_patch_array(SYNC_C, p_patch, zaux)
      fr_land_smt(:,:) = zaux(:,:)

    ENDDO

  END SUBROUTINE smooth_frland

  !-----------------------------------------------------------------------
  !!
  !! Smoothed urban fraction for a provisional tuning of the sealed fraction
  !!
  SUBROUTINE smooth_urbfrac (p_patch, p_int, fr_urb, fr_urb_smt)

    TYPE(t_patch)      , INTENT(IN)    :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(IN)    :: fr_urb(:,:)
    REAL(wp)           , INTENT(OUT)   :: fr_urb_smt(:,:)

    ! local variables
    INTEGER  :: jb, jc, jb1, jc1, j
    INTEGER  :: i_startblk, nblks_c, i_startidx, i_endidx

    i_startblk = p_patch%cells%start_blk(3,1)
    nblks_c    = p_patch%nblks_c

    fr_urb_smt(:,1:i_startblk) = 0._wp

    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 3)

      fr_urb_smt(i_startidx:i_endidx,jb) = 0._wp

      DO j=1, p_int%cell_environ%max_nmbr_nghbr_cells
        DO jc = i_startidx, i_endidx
          IF (j <= p_int%cell_environ%nmbr_nghbr_cells(jc,jb)) THEN
            jc1 = p_int%cell_environ%idx(jc,jb,j)
            jb1 = p_int%cell_environ%blk(jc,jb,j)
            fr_urb_smt(jc,jb) = fr_urb_smt(jc,jb) + fr_urb(jc1,jb1)/REAL(p_int%cell_environ%nmbr_nghbr_cells(jc,jb),wp)
          ENDIF
        ENDDO
      ENDDO

    ENDDO  !jb

    CALL sync_patch_array(SYNC_C, p_patch, fr_urb_smt)


  END SUBROUTINE smooth_urbfrac

  !-----------------------------------------------------------------------
  !! Reduce SSO stdh and slope over glaciers depending on the ratio between SSO slope and resolved slope
  !!
  !! This should preferably done in extpar, but the results with a globally modified calculation of 
  !! the SSO parameters (by removing the grid-scale slope before calculating stdh and SSO-slope) 
  !! were quite ambivalent
  !!
  SUBROUTINE postproc_sso (p_patch, p_int, fr_glac, topography_c, sso_stdh, sso_sigma)

    TYPE(t_patch)      , INTENT(IN)    :: p_patch
    TYPE(t_int_state)  , INTENT(IN)    :: p_int
    REAL(wp)           , INTENT(IN)    :: fr_glac(:,:)
    REAL(wp)           , INTENT(IN)    :: topography_c(:,:)
    REAL(wp)           , INTENT(INOUT) :: sso_stdh(:,:)
    REAL(wp)           , INTENT(INOUT) :: sso_sigma(:,:)

    ! local variables
    INTEGER  :: jb, jc
    INTEGER  :: i_startblk, i_endblk, nblks_c, i_startidx, i_endidx
    REAL(wp), DIMENSION(nproma,1,p_patch%nblks_c) :: z_topo_c, slope_x, slope_y, slope
    REAL(wp) :: slc_fac, slope_ratio, slred_fac, offset_fac, stdh_sv


    nblks_c    = p_patch%nblks_c

    i_startblk = p_patch%cells%start_block(2)
    i_endblk   = p_patch%cells%end_block(min_rlcell_int)

    z_topo_c(:,1,:) = topography_c(:,:)


    ! compute grid-scale slope
    !
    CALL rbf_interpol_c2grad(z_topo_c, p_patch, p_int, slope_x, slope_y)

    DO jb = i_startblk,i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, 2, min_rlcell_int)

      DO jc=i_startidx, i_endidx
        slope(jc,1,jb) = SQRT(slope_x(jc,1,jb)**2 + slope_y(jc,1,jb)**2)
      ENDDO

    ENDDO

    CALL sync_patch_array(SYNC_C, p_patch, slope)

    ! Reduce SSO stdh and slope depending on the glacier fration and the ratio between SSO slope and grid-scale slope
    ! With Merit/Rema raw data, it turned out to be beneficial to extend the slope correction the whole Arctic region
    DO jb = i_startblk,nblks_c

      CALL get_indices_c(p_patch, jb, i_startblk, nblks_c, i_startidx, i_endidx, 2)

      DO jc = i_startidx, i_endidx

        slope_ratio = MIN(1._wp, slope(jc,1,jb)/(1.e-4_wp + sso_sigma(jc,jb)))
        slc_fac = MAX(0._wp, 2._wp*(fr_glac(jc,jb) - 0.5_wp))
        IF (pp_sso == 1) THEN
          slred_fac = (1._wp - 0.6_wp*slc_fac*slope_ratio)
          sso_stdh(jc,jb)  = slred_fac*sso_stdh(jc,jb)
          sso_sigma(jc,jb) = slred_fac*sso_sigma(jc,jb)
        ELSE
          IF (p_patch%cells%center(jc,jb)%lat*rad2deg > 66.5_wp) slc_fac = 1._wp
          offset_fac = slc_fac
          slred_fac = MAX(0.05_wp,1._wp - slc_fac*slope_ratio)
          stdh_sv   = MAX(1.e-4_wp,sso_stdh(jc,jb))
          sso_stdh(jc,jb)  = MAX(0._wp, MIN(slred_fac*sso_stdh(jc,jb),sso_stdh(jc,jb)-offset_fac*10._wp) )
          sso_sigma(jc,jb) = (sso_stdh(jc,jb)/stdh_sv)**1.5_wp*sso_sigma(jc,jb)
        ENDIF

      ENDDO
    ENDDO


  END SUBROUTINE postproc_sso


END MODULE mo_process_topo

