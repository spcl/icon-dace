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

! Some utilities which are specific to the transport algorithm.
!
! Module contains some functions and procedures which are specifically related
! to the transport schemes. These subroutines or functions are needed at
! various places within the transport scheme. Therefore outsourcing these
! routines protects from possible circular dependencies.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_quadrature

  USE mo_kind,                ONLY: wp, vp
  USE mo_advection_config,    ONLY: t_gauss_quad_2d, gaussq_2d_o1, gaussq_2d_o2
  USE mo_advection_utils,     ONLY: t_list2D
  USE mo_model_domain,        ONLY: t_patch
  USE mo_loopindices,         ONLY: get_indices_e
  USE mo_impl_constants,      ONLY: min_rledge_int
  USE mo_math_constants,      ONLY: dbl_eps, eps

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: prep_gauss_quadrature_l
  PUBLIC :: prep_gauss_quadrature_l_list
  PUBLIC :: prep_gauss_quadrature_q
  PUBLIC :: prep_gauss_quadrature_q_list
  PUBLIC :: prep_gauss_quadrature_c
  PUBLIC :: prep_gauss_quadrature_c_list
  PUBLIC :: prep_gauss_quadrature_q_miura3
  PUBLIC :: prep_gauss_quadrature_c_miura3


CONTAINS


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of linear tracer subgrid distribution
  !! over a quadrilateral element.
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 1.
  !! I.e. a single quadrature point in physical space and the product of weights
  !! and the determinant of the Jacobian for the quadrature point.
  !! This subroutine is specific to a linear polynomial. It needs to be called 
  !! only once per time step, independent of the number of advected fields.
  !!
  SUBROUTINE prep_gauss_quadrature_l( p_patch, p_coords_dreg_v,         &
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                             !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,3,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pt_x,         &    !< in physical space
      &  z_gauss_pt_y

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac               !< gaussian quadrature point.

    REAL(wp) :: z_x(4), z_y(4) !< storage for local coordinates

    INTEGER  :: jb, je, jk          !< loop index for blocks and edges, levels
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev, elev          !< vertical start and end level

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_x,z_y
#endif
  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! convenience pointer
    gq => gaussq_2d_o1


    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    !$ACC DATA PRESENT(p_coords_dreg_v, p_quad_vector_sum, p_dreg_area) &
    !$ACC   CREATE(z_x, z_y)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,i_startidx,i_endidx,z_gauss_pt_x,z_gauss_pt_y,wgt_t_detjac,z_x,z_y &
!$OMP ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)


      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq)
      !$ACC LOOP GANG VECTOR PRIVATE(z_x, z_y, z_gauss_pt_x, z_gauss_pt_y) COLLAPSE(2)
      DO jk = slev, elev

        DO je = i_startidx, i_endidx

          z_x(1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get coordinates of the quadrature points in physical space (mapping)
!WS: TODO:  make sure that DOT_PRODUCT is supported in this OpenACC context
          z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,1),z_x(1:4))
          z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,1),z_y(1:4))

          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value'
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac = ( jac(z_x(1:4),z_y(1:4),gq%zeta(1),gq%eta(1)) * gq%wgt(1)) + dbl_eps

          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac. No summation necessary, since a
          ! single integration point is used.
          !
          ! const
          p_quad_vector_sum(je,1,jk,jb) = wgt_t_detjac
          ! x
          p_quad_vector_sum(je,2,jk,jb) = wgt_t_detjac * z_gauss_pt_x
          ! y
          p_quad_vector_sum(je,3,jk,jb) = wgt_t_detjac * z_gauss_pt_y

          ! area of departure region
          p_dreg_area(je,jk,jb) = wgt_t_detjac

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_l


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of linear tracer subgrid distribution
  !! over a quadrilateral element.
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 1.
  !! I.e. a single quadrature point in physical space and the product of weights
  !! and the determinant of the Jacobian for the quadrature point.
  !! This subroutine is specific to a linear polynomial. It needs to be called 
  !! only once per time step, independent of the number of advected fields.
  !!
  !! Index-list based version. Otherwise identical to prep_gauss_quadrature_l
  !!
  SUBROUTINE prep_gauss_quadrature_l_list( p_patch, p_coords_dreg_v, falist, &
    &                                 p_quad_vector_sum, p_dreg_area,        &
    &                                 opt_rlstart, opt_rlend                 )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                             !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:)   !< in 2D cartesian coordinates
                                    !< dim: (npoints,4,2,nblks_e)

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard 
      &  falist                     !< Miura-type treatment of flux areas is 
                                    !< insufficient

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:)   !< dim: (npoints,3,nblks_e)

    REAL(vp), INTENT(INOUT) :: &    !< total area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pt_x,         &    !< in physical space
      &  z_gauss_pt_y

    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac               !< gaussian quadrature point.

    REAL(wp) :: z_x(4), z_y(4) !< storage for local coordinates

    INTEGER  :: jb, je, jk          !< loop index for blocks and edges, levels
    INTEGER  :: ie                  !< index list loop counter
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_x,z_y
#endif
  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    ! convenience pointer
    gq => gaussq_2d_o1

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    !$ACC DATA PRESENT(p_coords_dreg_v, falist, falist%len, falist%eidx, falist%elev, p_dreg_area) &
    !$ACC   PRESENT(p_quad_vector_sum) &
    !$ACC   CREATE(z_x, z_y)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,ie,z_gauss_pt_x,z_gauss_pt_y,wgt_t_detjac,z_x,z_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq)
      !$ACC LOOP GANG VECTOR PRIVATE(je, jk, z_gauss_pt_x, z_gauss_pt_y)
!$NEC ivdep
      DO ie = 1, falist%len(jb)

        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)

        z_x(1:4) = p_coords_dreg_v(ie,1:4,1,jb)
        z_y(1:4) = p_coords_dreg_v(ie,1:4,2,jb)

        ! get coordinates of the quadrature points in physical space (mapping)
        z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,1),z_x(1:4))
        z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,1),z_y(1:4))

        ! get Jacobian determinant for each quadrature point and multiply with
        ! corresponding weights
        ! Note: dbl_eps is added, in order to have a meaningful 'edge value'
        ! (better: area-average) even when the integration-area tends to zero.
        wgt_t_detjac = ( jac(z_x(1:4),z_y(1:4),gq%zeta(1),gq%eta(1)) * gq%wgt(1) ) + dbl_eps


        ! Get quadrature vector for each integration point and multiply by
        ! corresponding wgt_t_detjac. No summation necessary, since a
        ! single integration point is used.
        !
        ! const
        p_quad_vector_sum(ie,1,jb) = wgt_t_detjac
        ! x
        p_quad_vector_sum(ie,2,jb) = wgt_t_detjac * z_gauss_pt_x
        ! y
        p_quad_vector_sum(ie,3,jb) = wgt_t_detjac * z_gauss_pt_y

        ! Add contribution to total area of departure region
        p_dreg_area(je,jk,jb) = p_dreg_area(je,jk,jb) + wgt_t_detjac

      ENDDO ! ie: loop over index list
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_l_list

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of quadratic tracer subgrid distribution
  !! over a quadrilateral element.
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a reconstruction based on a quadratic
  !! polynomial. It needs to be called only once per time step, independent
  !! of the number of advected fields.
  !!
  SUBROUTINE prep_gauss_quadrature_q( p_patch, p_coords_dreg_v,         &
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,6,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

    ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pt_x,         &    !< in physical space
      &  z_gauss_pt_y
    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.

    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,6)

    REAL(wp) :: z_x(4), z_y(4)      !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients
    REAL(wp) :: z_area                       !< auxiliary for dreg area

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev, elev          !< vertical start and end level

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: wgt_t_detjac,z_quad_vector,z_x,z_y
!DIR$ ATTRIBUTES ALIGN :64 :: z_wgt,z_eta
#endif

  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    ! convenience pointer
    gq => gaussq_2d_o2

    z_wgt(1:4) = 0.0625_wp * gq%wgt(1:4)

    z_eta(1:4,1) = 1._wp - gq%eta(1:4)
    z_eta(1:4,2) = 1._wp + gq%eta(1:4)
    z_eta(1:4,3) = 1._wp - gq%zeta(1:4)
    z_eta(1:4,4) = 1._wp + gq%zeta(1:4)

    !$ACC DATA PRESENT(p_coords_dreg_v, p_quad_vector_sum, p_dreg_area) &
    !$ACC   COPYIN(z_wgt, z_eta) &
    !$ACC   CREATE(z_x, z_y, wgt_t_detjac, z_quad_vector)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pt_x,z_gauss_pt_y,&
!$OMP            wgt_t_detjac,z_quad_vector,z_x,z_y,z_area) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)


      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq)
      !$ACC LOOP SEQ
      DO jk = slev, elev

        !$ACC LOOP GANG VECTOR PRIVATE(jg, z_area)
        DO je = i_startidx, i_endidx

          z_x(1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value'
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1:4) = dbl_eps + z_wgt(1:4) * ( &
            &   (z_eta(1:4,1)*(z_x(2)-z_x(1)) + z_eta(1:4,2)*(z_x(3)-z_x(4))) &
            & * (z_eta(1:4,3)*(z_y(4)-z_y(1)) - z_eta(1:4,4)*(z_y(2)-z_y(3))) &
            & - (z_eta(1:4,1)*(z_y(2)-z_y(1)) + z_eta(1:4,2)*(z_y(3)-z_y(4))) &
            & * (z_eta(1:4,3)*(z_x(4)-z_x(1)) - z_eta(1:4,4)*(z_x(2)-z_x(3))) )

          DO jg=1,4
            ! get coordinates of the quadrature points in physical space (mapping)
            z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,jg),z_x(1:4))
            z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,jg),z_y(1:4))

            ! const
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            ! x
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pt_x
            ! y
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pt_y
            ! x**2
            z_quad_vector(jg,4) = z_quad_vector(jg,2) * z_gauss_pt_x
            ! y**2
            z_quad_vector(jg,5) = z_quad_vector(jg,3) * z_gauss_pt_y
            ! xy
            z_quad_vector(jg,6) = z_quad_vector(jg,2) * z_gauss_pt_y
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je,1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je,2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je,3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je,4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je,5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je,6,jk,jb) = SUM(z_quad_vector(:,6))

          ! area of departure region
          z_area = SUM(wgt_t_detjac(1:4))
          p_dreg_area(je,jk,jb) = SIGN(MAX(eps,ABS(z_area)),z_area)

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_q


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of a 2D quadratic polynomial
  !! over a parallelogram-shaped element.
  !!
  !! Provides tracer-independent vector of polynomial points x^{k}y^{l}
  !! for a 2D polynomial of degree 2 (quadratic). The polynomial points are
  !! evaluated at the quadrature points and summed up, using a
  !! Gauss-Legendre quadrature of order 2 (two quadrature points in
  !! each direction).
  !! This routine must be called only once per time step, as it is
  !! independent of the advected fields.
  !!
  !! Note: This routine is only applicable for parallelogram-shaped
  !!       integration areas! This is because explicit use was made
  !!       of the fact that the Jacobian of the mapping is constant
  !!       and equals one quarter of the parallelogram area.
  !!
  SUBROUTINE prep_gauss_quadrature_q_miura3( p_patch, p_coords_dreg_v,  &
    &                                 p_quad_vector_sum,                &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev                          )

    TYPE(t_patch), INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                     !< performed

    REAL(vp),      INTENT(IN) ::  &  !< vertices of the flux area
      &  p_coords_dreg_v(:,:,:,:,:)  !< in 2D cartesian coordinates
                                     !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp),      INTENT(OUT)::  & !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,6,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation on halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,6)
    REAL(wp) :: z_x(4), z_y(4)      !< storage for local coordinates
    REAL(wp) :: z_wgt(4)

    REAL(wp) :: z_gauss_pt_x, z_gauss_pt_y

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev, elev          !< vertical start and end level

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_quad_vector,z_x,z_y,z_wgt
#endif

  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)


    ! convenience pointer
    gq => gaussq_2d_o2

    z_wgt(1:4) = 0.25_wp * gq%wgt(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_quad_vector, &
!$OMP            z_x,z_y,z_gauss_pt_x,z_gauss_pt_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq)
      !$ACC LOOP GANG VECTOR PRIVATE(z_x, z_y, z_quad_vector) COLLAPSE(2)
      DO jk = slev, elev
!NEC$ nolstval
        DO je = i_startidx, i_endidx

          ! vertices of the flux area
          z_x(1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding quadrature weight gq%wgt.
          DO jg = 1,4
            ! get coordinates of the quadrature point jg in physical space (mapping)
            z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,jg),z_x(1:4))
            z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,jg),z_y(1:4))

            ! const
            z_quad_vector(jg,1) = z_wgt(jg)
            ! x
            z_quad_vector(jg,2) = z_wgt(jg) * z_gauss_pt_x
            ! y
            z_quad_vector(jg,3) = z_wgt(jg) * z_gauss_pt_y
            ! x**2
            z_quad_vector(jg,4) = z_quad_vector(jg,2) * z_gauss_pt_x
            ! y**2
            z_quad_vector(jg,5) = z_quad_vector(jg,3) * z_gauss_pt_y
            ! xy
            z_quad_vector(jg,6) = z_quad_vector(jg,2) * z_gauss_pt_y
          ENDDO

          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je, 1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je, 2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je, 3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je, 4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je, 5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je, 6,jk,jb) = SUM(z_quad_vector(:,6))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prep_gauss_quadrature_q_miura3


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of quadratic tracer subgrid distribution
  !! over a quadrilateral element.
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a reconstruction based on a quadratic
  !! polynomial. It needs to be called only once per time step, independent
  !! of the number of advected fields.
  !!
  !! Index-list based version. Otherwise identical to prep_gauss_quadrature_q
  !!
  SUBROUTINE prep_gauss_quadrature_q_list( p_patch, p_coords_dreg_v, falist, &
    &                                 p_quad_vector_sum, p_dreg_area,        &
    &                                 opt_rlstart, opt_rlend                 )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:)   !< in 2D cartesian coordinates
                                    !< dim: (npoints,4,2,nblks_e)

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard 
      &  falist                     !< Miura-type treatment of flux areas is
                                    !< insufficient

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:)   !< dim: (npoints,6,nblks_e)

    REAL(vp), INTENT(INOUT) :: &    !< total area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pt_x,         &    !< in physical space
      &  z_gauss_pt_y

    REAL(wp) ::                       &     !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)                    !< each gaussian quadrature point.

    REAL(wp) ::                          &  !< quadrature vector for single integration point
      &  z_quad_vector(4,6)

    REAL(wp) :: z_x(4), z_y(4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: ie                  !< index list loop counter
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: wgt_t_detjac,z_quad_vector,z_x,z_y
!DIR$ ATTRIBUTES ALIGN :64 :: z_wgt,z_eta
#endif
  !-----------------------------------------------------------------------

   ! Check for optional arguments
    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    ! convenience pointer
    gq => gaussq_2d_o2

    z_wgt(1:4) = 0.0625_wp * gq%wgt(1:4)

    z_eta(1:4,1) = 1._wp - gq%eta(1:4)
    z_eta(1:4,2) = 1._wp + gq%eta(1:4)
    z_eta(1:4,3) = 1._wp - gq%zeta(1:4)
    z_eta(1:4,4) = 1._wp + gq%zeta(1:4)

    !$ACC DATA PRESENT(p_coords_dreg_v, falist, falist%len, falist%eidx, falist%elev, p_dreg_area) &
    !$ACC   PRESENT(p_quad_vector_sum) COPYIN(z_wgt, z_eta) &
    !$ACC   CREATE(z_x, z_y, z_quad_vector, wgt_t_detjac)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,ie,jg,z_gauss_pt_x,z_gauss_pt_y,wgt_t_detjac, &
!$OMP            z_quad_vector,z_x,z_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq)
      !$ACC LOOP GANG VECTOR
!$NEC ivdep
      DO ie = 1, falist%len(jb)

        z_x(1:4) = p_coords_dreg_v(ie,1:4,1,jb)
        z_y(1:4) = p_coords_dreg_v(ie,1:4,2,jb)

        ! get Jacobian determinant for each quadrature point and multiply with
        ! corresponding weights
        ! Note: dbl_eps is added, in order to have a meaningful 'edge value'
        ! (better: area-average) even when the integration-area tends to zero.
        wgt_t_detjac(1:4) = dbl_eps + z_wgt(1:4) * ( &
          &   (z_eta(1:4,1)*(z_x(2)-z_x(1)) + z_eta(1:4,2)*(z_x(3)-z_x(4))) &
          & * (z_eta(1:4,3)*(z_y(4)-z_y(1)) - z_eta(1:4,4)*(z_y(2)-z_y(3))) &
          & - (z_eta(1:4,1)*(z_y(2)-z_y(1)) + z_eta(1:4,2)*(z_y(3)-z_y(4))) &
          & * (z_eta(1:4,3)*(z_x(4)-z_x(1)) - z_eta(1:4,4)*(z_x(2)-z_x(3))) )

        ! Get quadrature vector for each integration point and multiply by
        ! corresponding wgt_t_detjac
        DO jg=1,4
          ! get coordinates of the quadrature point in physical space (mapping)
          z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,jg),z_x(1:4))
          z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,jg),z_y(1:4))

          ! const
          z_quad_vector(jg,1) = wgt_t_detjac(jg)
          ! x
          z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pt_x
          ! y
          z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pt_y
          ! x**2
          z_quad_vector(jg,4) = z_quad_vector(jg,2) * z_gauss_pt_x
          ! y**2
          z_quad_vector(jg,5) = z_quad_vector(jg,3) * z_gauss_pt_y
          ! xy
          z_quad_vector(jg,6) = z_quad_vector(jg,2) * z_gauss_pt_y
        ENDDO


        ! Sum quadrature vectors over all integration points
        p_quad_vector_sum(ie,1,jb) = SUM(z_quad_vector(:,1))
        p_quad_vector_sum(ie,2,jb) = SUM(z_quad_vector(:,2))
        p_quad_vector_sum(ie,3,jb) = SUM(z_quad_vector(:,3))
        p_quad_vector_sum(ie,4,jb) = SUM(z_quad_vector(:,4))
        p_quad_vector_sum(ie,5,jb) = SUM(z_quad_vector(:,5))
        p_quad_vector_sum(ie,6,jb) = SUM(z_quad_vector(:,6))

        ! Add contribution to total area of departure region
        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)
        p_dreg_area(je,jk,jb) = p_dreg_area(je,jk,jb) + SUM(wgt_t_detjac(1:4))

      ENDDO ! ie: loop over index list
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC END DATA

  END SUBROUTINE prep_gauss_quadrature_q_list


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of cubic tracer subgrid distribution
  !! over a quadrilateral element.
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a cubic reconstruction.
  !! It needs to be called only once per time step, independent of the number
  !! of advected fields.
  !!
  SUBROUTINE prep_gauss_quadrature_c( p_patch, p_coords_dreg_v,         &
    &                                 p_quad_vector_sum, p_dreg_area,   &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev                          )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:,:) !< in 2D cartesian coordinates
                                    !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:) !< dim: (nproma,10,nlev,nblks_e)

    REAL(vp), INTENT(OUT) :: &      !< area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pt_x,         &    !< in physical space
      &  z_gauss_pt_y
    REAL(wp) ::                &    !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)            !< each gaussian quadrature point.
    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,10)
    REAL(wp) :: z_x(4), z_y(4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients
    REAL(wp) :: z_area                       !< auxiliary for dreg area

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev, elev          !< vertical start and end level

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: wgt_t_detjac,z_quad_vector,z_x,z_y
!DIR$ ATTRIBUTES ALIGN :64 :: z_wgt,z_eta
#endif

  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF


    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    ! convenience pointer
    gq => gaussq_2d_o2

    z_wgt(1:4) = 0.0625_wp * gq%wgt(1:4)

    z_eta(1:4,1) = 1._wp - gq%eta(1:4)
    z_eta(1:4,2) = 1._wp + gq%eta(1:4)
    z_eta(1:4,3) = 1._wp - gq%zeta(1:4)
    z_eta(1:4,4) = 1._wp + gq%zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_gauss_pt_x,z_gauss_pt_y,&
!$OMP            wgt_t_detjac,z_quad_vector,z_x,z_y,z_area) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq) COPYIN(z_eta, z_wgt)
      !$ACC LOOP GANG VECTOR PRIVATE(z_x, z_y, wgt_t_detjac, z_quad_vector) COLLAPSE(2)
      DO jk = slev, elev

        DO je = i_startidx, i_endidx
          z_x(1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)

          ! get Jacobian determinant for each quadrature point and multiply with
          ! corresponding weights
          ! Note: dbl_eps is added, in order to have a meaningful 'edge value'
          ! (better: area-average) even when the integration-area tends to zero.
          wgt_t_detjac(1:4) = dbl_eps + z_wgt(1:4) * ( &
            &   (z_eta(1:4,1)*(z_x(2)-z_x(1)) + z_eta(1:4,2)*(z_x(3)-z_x(4))) &
            & * (z_eta(1:4,3)*(z_y(4)-z_y(1)) - z_eta(1:4,4)*(z_y(2)-z_y(3))) &
            & - (z_eta(1:4,1)*(z_y(2)-z_y(1)) + z_eta(1:4,2)*(z_y(3)-z_y(4))) &
            & * (z_eta(1:4,3)*(z_x(4)-z_x(1)) - z_eta(1:4,4)*(z_x(2)-z_x(3))) )

          ! Get quadrature vector for each integration point and multiply by
          ! corresponding wgt_t_detjac
          DO jg=1,4
            ! get coordinates of the quadrature point in physical space (mapping)
            z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,jg),z_x(1:4))
            z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,jg),z_y(1:4))

            ! const
            z_quad_vector(jg,1) = wgt_t_detjac(jg)
            ! x
            z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pt_x
            ! y
            z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pt_y
            ! x**2
            z_quad_vector(jg,4) = z_quad_vector(jg,2) * z_gauss_pt_x
            ! y**2
            z_quad_vector(jg,5) = z_quad_vector(jg,3) * z_gauss_pt_y
            ! xy
            z_quad_vector(jg,6) = z_quad_vector(jg,2) * z_gauss_pt_y
            ! x**3
            z_quad_vector(jg,7) = z_quad_vector(jg,4) * z_gauss_pt_x
            ! y**3
            z_quad_vector(jg,8) = z_quad_vector(jg,5) * z_gauss_pt_y
            ! x**2y
            z_quad_vector(jg,9) = z_quad_vector(jg,4) * z_gauss_pt_y
            ! xy**2
            z_quad_vector(jg,10)= z_quad_vector(jg,5) * z_gauss_pt_x
          ENDDO


          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je, 1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je, 2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je, 3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je, 4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je, 5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je, 6,jk,jb) = SUM(z_quad_vector(:,6))
          p_quad_vector_sum(je, 7,jk,jb) = SUM(z_quad_vector(:,7))
          p_quad_vector_sum(je, 8,jk,jb) = SUM(z_quad_vector(:,8))
          p_quad_vector_sum(je, 9,jk,jb) = SUM(z_quad_vector(:,9))
          p_quad_vector_sum(je,10,jk,jb) = SUM(z_quad_vector(:,10))

          ! area of departure region
          z_area = SUM(wgt_t_detjac(1:4))
          p_dreg_area(je,jk,jb) = SIGN(MAX(eps,ABS(z_area)),z_area)

!!$IF (p_dreg_area(je,jk,jb) < 0._wp) THEN
!!$  WRITE(0,*) "ATTENTION: negative areas at je,jk,jb= ", je, jk, jb, p_dreg_area(je,jk,jb)
!!$  WRITE(0,*) "system orientation: ", p_patch%edges%tangent_orientation(je,jb)
!!$  ELSE IF ((p_dreg_area(je,jk,jb) >= 0._wp)) THEN
!!$  WRITE(0,*) "OK for system orientation= ", je, jk, jb, p_patch%edges%tangent_orientation(je,jb)
!!$ENDIF
        ENDDO ! loop over edges

      ENDDO  ! loop over levels
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prep_gauss_quadrature_c


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of a 2D cubic polynomial 
  !! over a parallelogram-shaped element.
  !!
  !! Provides tracer-independent vector of polynomial points x^{k}y^{l}
  !! for a 2D polynomial of degree 3 (cubic). The polynomial points are
  !! evaluated at the quadrature points and summed up, using a
  !! Gauss-Legendre quadrature of order 2 (two quadrature points in
  !! each direction).
  !! This routine must be called only once per time step, as it is
  !! independent of the advected fields.
  !!
  !! Note: This routine is only applicable for parallelogram-shaped
  !!       integration areas! This is because explicit use was made
  !!       of the fact that the Jacobian of the mapping is constant
  !!       and equals one quarter of the parallelogram area.
  !!
  SUBROUTINE prep_gauss_quadrature_c_miura3( p_patch, p_coords_dreg_v,  &
    &                                 p_quad_vector_sum,                &
    &                                 opt_rlstart, opt_rlend, opt_slev, &
    &                                 opt_elev                          )

    TYPE(t_patch), INTENT(IN) :: &   !< patch on which computation is
      &  p_patch                     !< performed

    REAL(vp),      INTENT(IN) :: &   !< vertices of the flux area
      &  p_coords_dreg_v(:,:,:,:,:)  !< in 2D cartesian coordinates
                                     !< dim: (nproma,4,2,nlev,nblks_e)

    REAL(vp),      INTENT(OUT):: &   !< quadrature vector
      &  p_quad_vector_sum(:,:,:,:)  !< dim: (nproma,10,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation on halo points)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical start level
      &  opt_slev

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional vertical end level
      &  opt_elev

   ! local variables
    REAL(wp) ::                &    !< quadrature vector for single integration point
      &  z_quad_vector(4,10)
    REAL(wp) :: z_x(4), z_y(4)      !< storage for local coordinates
    REAL(wp) :: z_wgt(4)

    REAL(wp) :: z_gauss_pt_x, z_gauss_pt_y

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: i_startidx, i_endidx, i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend
    INTEGER  :: slev, elev          !< vertical start and end level

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: z_quad_vector,z_x,z_y,z_wgt
#endif

  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_slev) ) THEN
      slev = opt_slev
    ELSE
      slev = 1
    END IF

    IF ( PRESENT(opt_elev) ) THEN
      elev = opt_elev
    ELSE
      elev = p_patch%nlev
    END IF

    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)


    ! convenience pointer
    gq => gaussq_2d_o2

    z_wgt(1:4) = 0.25_wp * gq%wgt(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,jg,i_startidx,i_endidx,z_quad_vector, &
!$OMP            z_x,z_y,z_gauss_pt_x,z_gauss_pt_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, i_rlstart, i_rlend)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq)
      !$ACC LOOP GANG VECTOR PRIVATE(z_x, z_y, z_quad_vector) COLLAPSE(2)
      DO jk = slev, elev
!NEC$ nolstval
        DO je = i_startidx, i_endidx

          ! vertices of the flux area
          z_x(1:4) = p_coords_dreg_v(je,1:4,1,jk,jb)
          z_y(1:4) = p_coords_dreg_v(je,1:4,2,jk,jb)


          ! Get quadrature vector for each integration point and multiply by
          ! corresponding quadrature weight gq%wgt.
          DO jg = 1,4
            ! get coordinates of the quadrature point jg in physical space (mapping)
            z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,jg),z_x(1:4))
            z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,jg),z_y(1:4))

            ! const
            z_quad_vector(jg,1) = z_wgt(jg)
            ! x
            z_quad_vector(jg,2) = z_wgt(jg) * z_gauss_pt_x
            ! y
            z_quad_vector(jg,3) = z_wgt(jg) * z_gauss_pt_y
            ! x**2
            z_quad_vector(jg,4) = z_quad_vector(jg,2) * z_gauss_pt_x
            ! y**2
            z_quad_vector(jg,5) = z_quad_vector(jg,3) * z_gauss_pt_y
            ! xy
            z_quad_vector(jg,6) = z_quad_vector(jg,2) * z_gauss_pt_y
            ! x**3
            z_quad_vector(jg,7) = z_quad_vector(jg,4) * z_gauss_pt_x
            ! y**3
            z_quad_vector(jg,8) = z_quad_vector(jg,5) * z_gauss_pt_y
            ! x**2y
            z_quad_vector(jg,9) = z_quad_vector(jg,4) * z_gauss_pt_y
            ! xy**2
            z_quad_vector(jg,10)= z_quad_vector(jg,5) * z_gauss_pt_x
          ENDDO

          ! Sum quadrature vectors over all integration points
          p_quad_vector_sum(je, 1,jk,jb) = SUM(z_quad_vector(:,1))
          p_quad_vector_sum(je, 2,jk,jb) = SUM(z_quad_vector(:,2))
          p_quad_vector_sum(je, 3,jk,jb) = SUM(z_quad_vector(:,3))
          p_quad_vector_sum(je, 4,jk,jb) = SUM(z_quad_vector(:,4))
          p_quad_vector_sum(je, 5,jk,jb) = SUM(z_quad_vector(:,5))
          p_quad_vector_sum(je, 6,jk,jb) = SUM(z_quad_vector(:,6))
          p_quad_vector_sum(je, 7,jk,jb) = SUM(z_quad_vector(:,7))
          p_quad_vector_sum(je, 8,jk,jb) = SUM(z_quad_vector(:,8))
          p_quad_vector_sum(je, 9,jk,jb) = SUM(z_quad_vector(:,9))
          p_quad_vector_sum(je,10,jk,jb) = SUM(z_quad_vector(:,10))

        ENDDO ! loop over edges

      ENDDO  ! loop over levels
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prep_gauss_quadrature_c_miura3


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Prepares integration of cubic tracer subgrid distribution
  !! over a quadrilateral element.
  !!
  !! Provides tracer-independent parts for a gauss quadrature of order 2.
  !! I.e. the 4 quadrature points in physical space and the product of weights
  !! and the determinant of the Jacobian for each quadrature point.
  !! This subroutine is specific to a cubic reconstruction.
  !! It needs to be called only once per time step, independent of the number
  !! of advected fields.
  !!
  !! Index-list based version. Otherwise identical to prep_gauss_quadrature_c
  !!
  SUBROUTINE prep_gauss_quadrature_c_list( p_patch, p_coords_dreg_v, falist, &
    &                                 p_quad_vector_sum, p_dreg_area,        &
    &                                 opt_rlstart, opt_rlend                 )

    TYPE(t_patch), TARGET, INTENT(IN) ::  &  !< patch on which computation is
      &  p_patch                           !< performed

    REAL(vp), INTENT(IN)  ::   &    !< vertices of departure regions
      &  p_coords_dreg_v(:,:,:,:)   !< in 2D cartesian coordinates
                                    !< dim: (npoints,4,2,nblks_e)

    TYPE(t_list2D), INTENT(IN) :: & !< index list with points for which the standard 
      &  falist                     !< Miura-type treatment of flux areas is 
                                    !< insufficient

    REAL(vp), INTENT(OUT) :: &      !< quadrature vector
      &  p_quad_vector_sum(:,:,:)   !< dim: (npoints,10,nblks_e)

    REAL(vp), INTENT(INOUT) :: &    !< total area of departure region  [m**2]
      &  p_dreg_area(:,:,:)         !< dim: (nproma,nlev,nblks_e)

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control start level
      &  opt_rlstart

    INTEGER, INTENT(IN), OPTIONAL :: & !< optional: refinement control end level
      &  opt_rlend                     !< (to avoid calculation of halo points)


   ! local variables
    REAL(wp) ::                &    !< coordinates of gaussian quadrature points
      &  z_gauss_pt_x,         &    !< in physical space
      &  z_gauss_pt_y

    REAL(wp) ::                       &     !< weights times determinant of Jacobian for
      &  wgt_t_detjac(4)     !< each gaussian quadrature point.

    REAL(wp) ::                           & !< quadrature vector for single integration point
      &  z_quad_vector(4,10)

    REAL(wp) :: z_x(4), z_y(4) !< storage for local coordinates
    REAL(wp) :: z_wgt(4), z_eta(4,4)         !< for precomputation of coefficients

    INTEGER  :: jb, je, jk, jg      !< loop index for blocks and edges, levels and
                                    !< integration points
    INTEGER  :: ie                  !< index list loop counter
    INTEGER  :: i_startblk, i_endblk
    INTEGER  :: i_rlstart, i_rlend

    TYPE(t_gauss_quad_2d), POINTER :: gq

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: wgt_t_detjac,z_quad_vector,z_x,z_y
!DIR$ ATTRIBUTES ALIGN :64 :: z_wgt,z_eta
#endif

  !-----------------------------------------------------------------------

    ! Check for optional arguments
    IF ( PRESENT(opt_rlstart) ) THEN
      i_rlstart = opt_rlstart
    ELSE
      i_rlstart = 4
    ENDIF

    IF ( PRESENT(opt_rlend) ) THEN
      i_rlend = opt_rlend
    ELSE
      i_rlend = min_rledge_int - 2
    ENDIF

    i_startblk = p_patch%edges%start_block(i_rlstart)
    i_endblk   = p_patch%edges%end_block(i_rlend)

    ! convenience pointer
    gq => gaussq_2d_o2

    z_wgt(1:4) = 0.0625_wp * gq%wgt(1:4)

    z_eta(1:4,1) = 1._wp - gq%eta(1:4)
    z_eta(1:4,2) = 1._wp + gq%eta(1:4)
    z_eta(1:4,3) = 1._wp - gq%zeta(1:4)
    z_eta(1:4,4) = 1._wp + gq%zeta(1:4)

!$OMP PARALLEL
!$OMP DO PRIVATE(je,jk,jb,ie,jg,z_gauss_pt_x,z_gauss_pt_y,&
!$OMP            wgt_t_detjac,z_quad_vector,z_x,z_y) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) PRESENT(gq) COPYIN(z_wgt, z_eta)
      !$ACC LOOP GANG VECTOR PRIVATE(wgt_t_detjac, z_quad_vector, z_x, z_y)
!$NEC ivdep
      DO ie = 1, falist%len(jb)

        z_x(1:4) = p_coords_dreg_v(ie,1:4,1,jb)
        z_y(1:4) = p_coords_dreg_v(ie,1:4,2,jb)

        ! get Jacobian determinant for each quadrature point and multiply with
        ! corresponding weights
        ! Note: dbl_eps is added, in order to have a meaningful 'edge value'
        ! (better: area-average) even when the integration-area tends to zero.
        wgt_t_detjac(1:4) = dbl_eps + z_wgt(1:4) * ( &
          &   (z_eta(1:4,1)*(z_x(2)-z_x(1)) + z_eta(1:4,2)*(z_x(3)-z_x(4))) &
          & * (z_eta(1:4,3)*(z_y(4)-z_y(1)) - z_eta(1:4,4)*(z_y(2)-z_y(3))) &
          & - (z_eta(1:4,1)*(z_y(2)-z_y(1)) + z_eta(1:4,2)*(z_y(3)-z_y(4))) &
          & * (z_eta(1:4,3)*(z_x(4)-z_x(1)) - z_eta(1:4,4)*(z_x(2)-z_x(3))) )

        ! Get quadrature vector for each integration point and multiply by
        ! corresponding wgt_t_detjac
        DO jg=1,4
          ! get coordinates of the quadrature point in physical space (mapping)
          z_gauss_pt_x = DOT_PRODUCT(gq%shape_func(1:4,jg),z_x(1:4))
          z_gauss_pt_y = DOT_PRODUCT(gq%shape_func(1:4,jg),z_y(1:4))

          ! const
          z_quad_vector(jg,1) = wgt_t_detjac(jg)
          ! x
          z_quad_vector(jg,2) = wgt_t_detjac(jg) * z_gauss_pt_x
          ! y
          z_quad_vector(jg,3) = wgt_t_detjac(jg) * z_gauss_pt_y
          ! x**2
          z_quad_vector(jg,4) = z_quad_vector(jg,2) * z_gauss_pt_x
          ! y**2
          z_quad_vector(jg,5) = z_quad_vector(jg,3) * z_gauss_pt_y
          ! xy
          z_quad_vector(jg,6) = z_quad_vector(jg,2) * z_gauss_pt_y
          ! x**3
          z_quad_vector(jg,7) = z_quad_vector(jg,4) * z_gauss_pt_x
          ! y**3
          z_quad_vector(jg,8) = z_quad_vector(jg,5) * z_gauss_pt_y
          ! x**2y
          z_quad_vector(jg,9) = z_quad_vector(jg,4) * z_gauss_pt_y
          ! xy**2
          z_quad_vector(jg,10)= z_quad_vector(jg,5) * z_gauss_pt_x
        ENDDO


        ! Sum quadrature vectors over all integration points
        p_quad_vector_sum(ie, 1,jb) = SUM(z_quad_vector(:,1))
        p_quad_vector_sum(ie, 2,jb) = SUM(z_quad_vector(:,2))
        p_quad_vector_sum(ie, 3,jb) = SUM(z_quad_vector(:,3))
        p_quad_vector_sum(ie, 4,jb) = SUM(z_quad_vector(:,4))
        p_quad_vector_sum(ie, 5,jb) = SUM(z_quad_vector(:,5))
        p_quad_vector_sum(ie, 6,jb) = SUM(z_quad_vector(:,6))
        p_quad_vector_sum(ie, 7,jb) = SUM(z_quad_vector(:,7))
        p_quad_vector_sum(ie, 8,jb) = SUM(z_quad_vector(:,8))
        p_quad_vector_sum(ie, 9,jb) = SUM(z_quad_vector(:,9))
        p_quad_vector_sum(ie,10,jb) = SUM(z_quad_vector(:,10))

        ! Add contribution to total area of departure region
        je = falist%eidx(ie,jb)
        jk = falist%elev(ie,jb)
        p_dreg_area(je,jk,jb) = p_dreg_area(je,jk,jb) + SUM(wgt_t_detjac(1:4))

      ENDDO ! ie: loop over index list
      !$ACC END PARALLEL

    ENDDO  ! loop over blocks
    !$ACC WAIT(1)

!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE prep_gauss_quadrature_c_list


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computes Jacobian determinant for gaussian quadrature
  !!
  !! Computes Jacobian determinant for gaussian quadrature
  !!
  FUNCTION jac(x, y, zeta, eta)  RESULT(det_jac)
    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN) :: x(1:4), y(1:4)  !< coordinates of vertices in x-y-system
    REAL(wp), INTENT(IN) :: zeta, eta       !< integration point in \zeta,\eta-system

    ! RETURN VALUE:
    REAL(wp) :: det_jac

    REAL(wp), DIMENSION(2,2) :: jacob
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: jacob
#endif

  !-----------------------------------------------------------------------

    jacob(1,1) = (1._wp - eta) * ( x(2) - x(1))   &
      &        + (1._wp + eta) * ( x(3) - x(4))
    jacob(1,2) = (1._wp - eta) * ( y(2) - y(1))   &
      &        + (1._wp + eta) * ( y(3) - y(4))
    jacob(2,1) = (1._wp - zeta)* ( x(4) - x(1))   &
      &        - (1._wp + zeta)* ( x(2) - x(3))
    jacob(2,2) = (1._wp - zeta)* ( y(4) - y(1))   &
      &        - (1._wp + zeta)* ( y(2) - y(3))

    det_jac = 0.0625_wp * (jacob(1,1)*jacob(2,2) - jacob(1,2)*jacob(2,1))

  END FUNCTION jac


END MODULE mo_advection_quadrature

