! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.
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

MODULE mo_icon_interpolation_vector
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp, vp
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_model_domain,        ONLY: t_patch
  USE mo_run_config,          ONLY: timers_level
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_intp_data_strc,      ONLY: t_int_state
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_intp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: edges2cells_vector

CONTAINS



!------------------------------------------------------------------------
!!  Bilinear interpolation of normal and tangential velocity components
!!  at the edges to u and v at the cells
!!  Works only for triangles (bilinear interpolation weights are not implemented
!!  for hexagons)
!!
SUBROUTINE edges2cells_vector( p_vn_in, p_vt_in, ptr_patch, p_int, &
  &  p_u_out, p_v_out, opt_slev, opt_elev, opt_rlstart, opt_rlend  )
!


TYPE(t_patch), TARGET, INTENT(in) :: ptr_patch

! normal velocity component at edges
REAL(wp), INTENT(in) ::  p_vn_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! (reconstructed) tangential velocity component at edges
REAL(vp), INTENT(in) ::  p_vt_in(:,:,:)  ! dim: (nproma,nlev,nblks_e)

! Interpolation state
TYPE(t_int_state), INTENT(IN) :: p_int

! cell based output fields: u and v
REAL(wp), INTENT(inout) :: p_u_out(:,:,:), p_v_out(:,:,:) ! dim: (nproma,nlev,nblks_c)

INTEGER, INTENT(in), OPTIONAL ::  opt_slev ! optional vertical start level

INTEGER, INTENT(in), OPTIONAL ::  opt_elev ! optional vertical end level

! start and end values of refin_ctrl flag
INTEGER, INTENT(in), OPTIONAL ::  opt_rlstart, opt_rlend

INTEGER :: slev, elev     ! vertical start and end level
INTEGER :: jc, jk, jb
INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, i_nchdom
INTEGER :: rl_start, rl_end

INTEGER,  DIMENSION(:,:,:),   POINTER :: iidx, iblk

!-------------------------------------------------------------------------

! check optional arguments
IF ( PRESENT(opt_slev) ) THEN
  slev = opt_slev
ELSE
  slev = 1
END IF
IF ( PRESENT(opt_elev) ) THEN
  elev = opt_elev
ELSE
  elev = UBOUND(p_vn_in,2)
END IF

IF ( PRESENT(opt_rlstart) ) THEN
  rl_start = opt_rlstart
ELSE
  rl_start = 2
END IF
IF ( PRESENT(opt_rlend) ) THEN
  rl_end = opt_rlend
ELSE
  rl_end = min_rlcell_int
END IF

iidx => ptr_patch%cells%edge_idx
iblk => ptr_patch%cells%edge_blk

! values for the blocking
i_nchdom   = MAX(1,ptr_patch%n_childdom)
i_startblk = ptr_patch%cells%start_blk(rl_start,1)
i_endblk   = ptr_patch%cells%end_blk(rl_end,i_nchdom)


IF (timers_level > 10) CALL timer_start(timer_intp)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,jk,jc) ICON_OMP_DEFAULT_SCHEDULE
DO jb = i_startblk, i_endblk

  CALL get_indices_c(ptr_patch, jb, i_startblk, i_endblk, &
                     i_startidx, i_endidx, rl_start, rl_end)


!CDIR UNROLL=6
  DO jk = slev, elev
    DO jc = i_startidx, i_endidx

      p_u_out(jc,jk,jb) =  &
        p_int%e_bln_c_u(jc,1,jb)*p_vn_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_u(jc,2,jb)*p_vt_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_u(jc,3,jb)*p_vn_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_u(jc,4,jb)*p_vt_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_u(jc,5,jb)*p_vn_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
        p_int%e_bln_c_u(jc,6,jb)*p_vt_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

      p_v_out(jc,jk,jb) =  &
        p_int%e_bln_c_v(jc,1,jb)*p_vn_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_v(jc,2,jb)*p_vt_in(iidx(jc,jb,1),jk,iblk(jc,jb,1)) + &
        p_int%e_bln_c_v(jc,3,jb)*p_vn_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_v(jc,4,jb)*p_vt_in(iidx(jc,jb,2),jk,iblk(jc,jb,2)) + &
        p_int%e_bln_c_v(jc,5,jb)*p_vn_in(iidx(jc,jb,3),jk,iblk(jc,jb,3)) + &
        p_int%e_bln_c_v(jc,6,jb)*p_vt_in(iidx(jc,jb,3),jk,iblk(jc,jb,3))

     ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF (timers_level > 10) CALL timer_stop(timer_intp)


END SUBROUTINE edges2cells_vector


END MODULE mo_icon_interpolation_vector
