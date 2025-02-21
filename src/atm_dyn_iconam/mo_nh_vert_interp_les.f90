!
! This module has useful routines for LES runs
!
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

MODULE mo_nh_vert_interp_les

  USE mo_kind,                ONLY: wp
  USE mo_nonhydro_types,      ONLY: t_nh_metrics
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_intp,                ONLY: cells2verts_scalar, edges2cells_scalar
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_sync,                ONLY: global_sum_array, SYNC_C, SYNC_V, &
                                    sync_patch_array_mult
  USE mo_impl_constants_grf,  ONLY: grf_bdywidth_c
  USE mo_impl_constants,      ONLY: SUCCESS, min_rlcell_int
  USE mo_parallel_config,     ONLY: nproma, p_test_run
  USE mo_physical_constants,  ONLY: grav
  USE mo_les_config,          ONLY: les_config
  USE mo_exception,           ONLY: finish
  USE mo_fortran_tools,       ONLY: init, set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vert_intp_full2half_cell_3d, vert_intp_linear_1d, global_hor_mean
  PUBLIC :: vertical_derivative, brunt_vaisala_freq, init_vertical_grid_for_les

  CONTAINS


  !!------------------------------------------------------------------------
  !! init_vertical_grid_for_les

  SUBROUTINE init_vertical_grid_for_les(jg, p_patch, p_int, p_metrics)
    INTEGER,                   INTENT(in)     :: jg
    TYPE(t_patch),             INTENT(inout)  :: p_patch
    TYPE(t_int_state),         INTENT(in)     :: p_int
    TYPE(t_nh_metrics),        INTENT(inout)  :: p_metrics

    ! local variables
    CHARACTER(*), PARAMETER :: routine = &
        "mo_les_utilities:init_vertical_grid_for_les"

    IF(.NOT.les_config(jg)%les_metric) &
      RETURN

    IF (p_test_run) THEN
!$OMP PARALLEL
      CALL init(p_metrics%ddxt_z_half_v)
      CALL init(p_metrics%ddxn_z_half_c)
      CALL init(p_metrics%ddxn_z_full_c)
      CALL init(p_metrics%ddxn_z_full_v)
      CALL init(p_metrics%ddxt_z_half_c)
      CALL init(p_metrics%ddxt_z_full_c)
      CALL init(p_metrics%ddxt_z_full_v)
      CALL init(p_metrics%inv_ddqz_z_full_v)
!$OMP END PARALLEL
    END IF

    ! half_c sync
    CALL edges2cells_scalar(p_metrics%ddxn_z_half_e, p_patch, p_int%e_bln_c_s, &
      &                     p_metrics%ddxn_z_half_c)

    CALL edges2cells_scalar(p_metrics%ddxt_z_half_e, p_patch, p_int%e_bln_c_s, &
      &                     p_metrics%ddxt_z_half_c)

    ! full_c sync
    CALL edges2cells_scalar(p_metrics%ddxn_z_full, p_patch, p_int%e_bln_c_s, &
      &                     p_metrics%ddxn_z_full_c)

    CALL edges2cells_scalar(p_metrics%ddxt_z_full, p_patch, p_int%e_bln_c_s, &
      &                     p_metrics%ddxt_z_full_c)
    CALL sync_patch_array_mult(SYNC_C, p_patch, 4, p_metrics%ddxn_z_half_c, &
      &                        p_metrics%ddxt_z_half_c, &
      &                        p_metrics%ddxn_z_full_c, &
      &                        p_metrics%ddxt_z_full_c)

    ! full_v sync
    CALL cells2verts_scalar(p_metrics%ddxn_z_full_c, p_patch, &
      &                     p_int%cells_aw_verts, p_metrics%ddxn_z_full_v)

    CALL cells2verts_scalar(p_metrics%ddxt_z_full_c, p_patch, &
      &                     p_int%cells_aw_verts, p_metrics%ddxt_z_full_v)

    CALL cells2verts_scalar(p_metrics%inv_ddqz_z_full, p_patch, &
      &                     p_int%cells_aw_verts, p_metrics%inv_ddqz_z_full_v)

    ! half_v sync
    CALL cells2verts_scalar(p_metrics%ddxt_z_half_c, p_patch, &
         p_int%cells_aw_verts, p_metrics%ddxt_z_half_v)
    CALL sync_patch_array_mult(SYNC_V, p_patch, 4, p_metrics%ddxn_z_full_v, &
      &                        p_metrics%ddxt_z_full_v, &
      &                        p_metrics%inv_ddqz_z_full_v, &
      &                        p_metrics%ddxt_z_half_v)


  END SUBROUTINE init_vertical_grid_for_les


  !!------------------------------------------------------------------------
  !! vert_intp_full2half_3d

  SUBROUTINE vert_intp_full2half_cell_3d(p_patch, p_metrics, varin, varout, rl_start, rl_end, lacc)

    TYPE(t_nh_metrics),INTENT(in) :: p_metrics
    TYPE(t_patch),     INTENT(in) :: p_patch
    REAL(wp), INTENT(in)                  :: varin(:,:,:)
    INTEGER,  INTENT(in)                  :: rl_start, rl_end 
    REAL(wp), INTENT(out)                 :: varout(:,:,:)                     

    INTEGER :: i_startblk, i_endblk
    INTEGER :: i_endidx, i_startidx, nlevp1, nlev
    INTEGER :: jk, jc, jb

    LOGICAL, OPTIONAL, INTENT(in)   :: lacc  !< GPU flag
    LOGICAL                         :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   PRESENT(p_metrics, varin, varout) IF(lzacc)

    nlev      = p_patch%nlev
    nlevp1    = p_patch%nlev+1

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, rl_end)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)




        DO jk = 2 , nlev  
         DO jc = i_startidx , i_endidx

          varout(jc,jk,jb) = p_metrics%wgtfac_c(jc,jk,jb)*varin(jc,jk,jb) + &
                        (1._wp-p_metrics%wgtfac_c(jc,jk,jb))*varin(jc,jk-1,jb)
         END DO
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx, i_endidx
           varout(jc,1,jb) =                                &
             p_metrics%wgtfacq1_c(jc,1,jb)*varin(jc,1,jb) + &
             p_metrics%wgtfacq1_c(jc,2,jb)*varin(jc,2,jb) + &
             p_metrics%wgtfacq1_c(jc,3,jb)*varin(jc,3,jb)

           varout(jc,nlevp1,jb) =                               &
             p_metrics%wgtfacq_c(jc,1,jb)*varin(jc,nlev,jb)   + &
             p_metrics%wgtfacq_c(jc,2,jb)*varin(jc,nlev-1,jb) + &
             p_metrics%wgtfacq_c(jc,3,jb)*varin(jc,nlev-2,jb)
        END DO     
        !$ACC END PARALLEL
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

  !$ACC WAIT
  !$ACC END DATA

  END SUBROUTINE vert_intp_full2half_cell_3d


  !!------------------------------------------------------------------------
  !! linear vertical interpolation: grid za to zb
  !! - It extrapolates if no data given for zb > za
  !! Taken from UCLA-LES

  SUBROUTINE vert_intp_linear_1d(za, xa, zb, xb) 
     REAL(wp), INTENT(IN)  :: za(:), zb(:), xa(:)
     REAL(wp), INTENT(OUT) :: xb(:)
  
     REAL(wp) :: wt
     INTEGER  :: l, k, na, nb

     na = SIZE(za)
     nb = SIZE(zb)

     l = na
     DO k = nb, 1, -1
       IF (zb(k) <= za(1)) THEN
          DO WHILE ( zb(k) > za(l-1) .AND. l > 1)
             l = l-1
          END DO
          wt=(zb(k)-za(l))/(za(l-1)-za(l))
          xb(k)=xa(l)+(xa(l-1)-xa(l))*wt    
       ELSE
          wt=(zb(k)-za(1))/(za(2)-za(1))
          xb(k)=xa(1)+(xa(2)-xa(1))*wt
       END IF
    END DO

  END SUBROUTINE vert_intp_linear_1d

  !!------------------------------------------------------------------------
  !! global_hor_mean: only called for interior points
  !! Calculates horizontally averaged vertically varying quantaties 

  SUBROUTINE global_hor_mean(p_patch, var, varout, inv_no_cells)

    TYPE(t_patch),     INTENT(in) :: p_patch
    REAL(wp), INTENT(in)                  :: var(:,:,:), inv_no_cells
    REAL(wp), INTENT(out)                 :: varout(:)                     

    REAL(wp) :: var_aux(SIZE(var,1),SIZE(var,2),SIZE(var,3))
    INTEGER  :: i_startblk, i_endblk, rl_start
    INTEGER  :: i_endidx, i_startidx
    INTEGER  :: jk, jc, jb, nz, nblk, nproma

    rl_start   = grf_bdywidth_c+1
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(min_rlcell_int)
    nproma     = SIZE(var,1)
    nz         = SIZE(var,2)
    nblk       = SIZE(var,3)

   !Now put values in interior nodes
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, nblk
      IF (jb >= i_startblk .AND. jb <= i_endblk) THEN
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                           i_startidx, i_endidx, rl_start, min_rlcell_int)
        DO jk = 1, nz
          DO jc = 1, nproma
            var_aux(jc,jk,jb) = MERGE(var(jc,jk,jb), 0._wp, &
                 jc >= i_startidx .AND. jc <= i_endidx)
          END DO
        END DO
      ELSE
        var_aux(:,:,jb) = 0.0_wp
      END IF
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

   DO jk = 1 , nz
    varout(jk) = global_sum_array(var_aux(:,jk,:)) * inv_no_cells
   END DO

  END SUBROUTINE global_hor_mean

  !!------------------------------------------------------------------------
  !! vertical_derivative

  FUNCTION vertical_derivative (var, inv_dz) RESULT(dvardz)

    REAL(wp), INTENT(in) :: var(:), inv_dz(:)
                     
    REAL(wp) :: dvardz(SIZE(inv_dz))                     
    INTEGER  :: jk

!$OMP PARALLEL
!$OMP DO PRIVATE(jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jk = 1 , SIZE(inv_dz)
      dvardz(jk) = ( var(jk) - var(jk+1) ) * inv_dz(jk)
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END FUNCTION vertical_derivative

  !!------------------------------------------------------------------------
  !! Brunt Vaisala Frequency: 
  !! Calculates BVF for unsaturated and saturated case based on Durran & Klemp 1982
  !! Eq. 4. and using moist lapse rate expression from Marshall and Plumb

  SUBROUTINE brunt_vaisala_freq(p_patch, p_metrics, thetav, bru_vais, lacc)

    TYPE(t_patch), INTENT(in) :: p_patch
    TYPE(t_nh_metrics), INTENT(in) :: p_metrics
    REAL(wp), INTENT(in):: thetav(:,:,:)
    REAL(wp), INTENT(INOUT)               :: bru_vais(:,:,:)

    REAL(wp) :: thetav_ic(nproma,p_patch%nlev+1,p_patch%nblks_c)
    INTEGER  :: i_startblk, i_endblk, rl_start, rl_end
    INTEGER  :: i_endidx, i_startidx, nlev
    INTEGER  :: jk, jc, jb

    LOGICAL, OPTIONAL, INTENT(in)   :: lacc  !< GPU flag
    LOGICAL                         :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA &
    !$ACC   PRESENT(thetav, bru_vais, p_metrics, p_metrics%inv_ddqz_z_half) &
    !$ACC   CREATE(thetav_ic) IF(lzacc)

    !To be calculated at all cells at interface levels, except top/bottom 
    !boundaries

    nlev      = p_patch%nlev

    rl_start   = 2
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    CALL vert_intp_full2half_cell_3d(p_patch, p_metrics, thetav, thetav_ic, rl_start, rl_end, lacc=lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, rl_start, rl_end)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)




      DO jk = 2 , nlev
        DO jc = i_startidx , i_endidx

          bru_vais(jc,jk,jb) = grav * ( thetav(jc,jk-1,jb) - thetav(jc,jk,jb) ) * &
                               p_metrics%inv_ddqz_z_half(jc,jk,jb)/thetav_ic(jc,jk,jb)    
        END DO
      END DO     
      !$ACC END PARALLEL
    END DO 
!$OMP END DO NOWAIT
!$OMP END PARALLEL     

  !$ACC WAIT
  !$ACC END DATA
   
  END SUBROUTINE brunt_vaisala_freq

!-------------------------------------------------------------------------------

END MODULE mo_nh_vert_interp_les



