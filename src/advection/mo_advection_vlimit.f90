! Flux and reconstruction limiter for vertical tracer transport
!
! This module contains various filters and flux limiters for
! vertical tracer transport.
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
#include "omp_definitions.inc"
!----------------------------
MODULE mo_advection_vlimit

  USE mo_kind,                ONLY: wp
  USE mo_math_constants,      ONLY: dbl_eps
  USE mo_parallel_config,     ONLY: nproma, p_test_run

  IMPLICIT NONE

  PRIVATE


  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: vflx_limiter_pd
  PUBLIC :: v_limit_parabola_mo
  PUBLIC :: v_limit_parabola_sm
  PUBLIC :: v_limit_slope_mo
  PUBLIC :: v_limit_slope_sm
  PUBLIC :: v_limit_face_mc_mo
  PUBLIC :: v_limit_face_mc_sm


CONTAINS



  !-------------------------------------------------------------------------
  !>
  !! Positive definite flux limiter for vertical advection
  !!
  !! Positive definite Zalesak Flux-Limiter (Flux corrected transport).
  !! for the nonhydrostatic core. Only outward fluxes are re-scaled, in 
  !! order to maintain positive definiteness.
  !!
  !! @par Literature:
  !! - Zalesak, S.T. (1979): Fully Multidimensional Flux-corrected Transport
  !!   Algorithms for Fluids. JCP, 31, 335-362
  !! - Harris, L. M. and P. H. Lauritzen (2010): A flux-form version of the
  !!   Conservative Semi-Lagrangian Multi-tracer transport scheme (CSLAM) on
  !!   the cubed sphere grid.  J. Comput. Phys., 230, 1215-1237
  !! - Smolarkiewicz, P. K., 1989: Comment on "A positive definite advection 
  !!   scheme obtained by nonlinear renormalization of the advective fluxes.", 
  !!   Mon. Wea. Rev., 117, 2626-2632
  !!
  SUBROUTINE vflx_limiter_pd( p_dtime, p_cc, p_rhodz_now,             &
    &                         p_mflx_tracer_v, i_startidx, i_endidx,  &
    &                         slev, elev )

    REAL(wp), INTENT(IN) ::     &    !< time step [s]
      &  p_dtime

    REAL(wp), INTENT(IN) ::     &    !< advected cell centered variable at time (n)
      &  p_cc(:,:)                   !< dim: (nproma,nlev)
                                     !< [kg kg^-1]

    REAL(wp), INTENT(IN) ::     &    !< density times cell thickness at timestep n
      &  p_rhodz_now(:,:)            !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT) ::  &    !< calculated vertical tracer mass flux
      &  p_mflx_tracer_v(:,:)        !< dim: (nproma,nlevp1)
                                     !< [kg m^-2 s^-1]

    INTEGER, INTENT(IN) ::      &    !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN) ::      &    !< vertical start and end level
      &  slev, elev

    ! local
    !
    REAL(wp) ::                 &    !< fraction with which all outgoing fluxes 
      &  r_m(nproma,SIZE(p_cc,2))    !< of cell jk are multiplied 
                                     !< to guarantee positive definiteness

    REAL(wp) :: p_m                  !< sum of fluxes out of cell
                                     !< [kg m^-2]

    REAL(wp) :: z_signum             !< sign of mass flux
                                     !< >0: upward; <0: downward

    INTEGER  :: jk, jc               !< index of vert level, cell
    INTEGER  :: jkp1, jkm1

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN :64 :: r_m
#endif
  !-------------------------------------------------------------------------


    !$ACC DATA CREATE(r_m) PRESENT(p_cc, p_rhodz_now, p_mflx_tracer_v)

    IF (p_test_run) THEN
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1)
      r_m = 0._wp
      !$ACC END KERNELS
    ENDIF

    !
    ! 1. Compute total outward mass (loop over full levels)
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG COLLAPSE(2) PRIVATE(jkp1, p_m)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx
        jkp1 = jk+1

        ! Sum of all outgoing fluxes out of cell jk
        p_m = p_dtime                               &
          & * (MAX(0._wp,p_mflx_tracer_v(jc,jk))    &  ! upper half level
          &  - MIN(0._wp,p_mflx_tracer_v(jc,jkp1)) )   ! lower half level

        ! fraction with which all the fluxes out of cell jk must be multiplied 
        ! to guarantee no undershoot
        ! Nominator: maximum allowable mass loss \rho^n q^n
        r_m(jc,jk) = MIN(1._wp, (p_cc(jc,jk)*p_rhodz_now(jc,jk)) &
          &         /(p_m + dbl_eps) )

      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !
    ! 2. Limit outward fluxes (loop over half levels)
    !    Choose r_m depending on the sign of p_mflx_tracer_v
    !
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(jkm1, z_signum) COLLAPSE(2)
    DO jk = slev+1, elev
      DO jc = i_startidx, i_endidx
        jkm1 = jk-1

        ! NH:
        ! p_mflx_tracer_v(k-1/2) > 0: flux directed from cell k   -> k-1
        ! p_mflx_tracer_v(k-1/2) < 0: flux directed from cell k-1 -> k
        !
        z_signum = SIGN(1._wp,p_mflx_tracer_v(jc,jk))

        p_mflx_tracer_v(jc,jk) =  p_mflx_tracer_v(jc,jk)  * 0.5_wp      &
          &                       * ( (1._wp + z_signum) * r_m(jc,jk)   &
          &                       +   (1._wp - z_signum) * r_m(jc,jkm1) )
  
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE vflx_limiter_pd



  !-------------------------------------------------------------------------
  !>
  !! Filter for PPM/PSM (3rd order) vertical advection (monotonic version)
  !!
  !! Removes over- and undershoots in first guess parabola by resetting the
  !! upper or lower interface values.
  !! Avoids non-physical over/undershoots in advected fields.
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  SUBROUTINE v_limit_parabola_mo( p_ivlimit_selective, p_cc, p_face, &
    &                           p_face_up, p_face_low, i_startidx,   &
    &                           i_endidx, slev, elev )


    INTEGER, INTENT(IN) ::   &    !< avoids limiting of smooth extrema
      &  p_ivlimit_selective      !< if activated

    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::  &     !< reconstructed face values of the advected field
      &  p_face(:,:)               !< dim: (nproma,nlevp1)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (upper face, height based)
      &  p_face_up(:,:)            !< dim: (nproma,nlevp)

    REAL(wp), INTENT(INOUT) :: &   !< final face value (lower face, height based)
      &  p_face_low(:,:)           !< dim: (nproma,nlevp)

    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    ! local
    REAL(wp) :: z_delta(SIZE(p_cc,1),SIZE(p_cc,2)) !< lower minus upper face value
    REAL(wp) :: z_a6i(SIZE(p_cc,1),SIZE(p_cc,2))   !< curvature of parabola
    LOGICAL  :: l_limit(SIZE(p_cc,1),SIZE(p_cc,2)) !< is limiting of subgrid parabola 
                                                   !< in particular cell necessary [yes/no]
    LOGICAL  :: is_main_crit          !< is main criterion for limiter activation TRUE
    REAL(wp) :: q_face_up, q_face_low !< face values at upper and lower cell edge

    INTEGER  :: jc, jk             !< index of cell and vertical level
    INTEGER  :: jkm2, jkm1,    &   !< neighbour indices  
      &         jkp1, jkp2, jkp3
    INTEGER  :: elev_slim          !< end level for spurious extremum dectector
    !-----------------------------------------------------------------------

    !$ACC DATA PRESENT(p_cc, p_face, p_face_up, p_face_low) &
    !$ACC   CREATE(z_delta, z_a6i, l_limit)

    ! selective limitation yes or no
    IF (p_ivlimit_selective == 1) THEN
      ! 
      elev_slim = MIN(elev+1,UBOUND(p_face,2))
      !
      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(z_delta, z_a6i, l_limit) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jkp1, jkp2, jkp3, jkm1, jkm2, is_main_crit) COLLAPSE(2)
      DETECT_SEL:DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! index of bottom half level
          jkp1 = jk+1
          jkp2 = MIN(jk+2,elev_slim)
          jkp3 = MIN(jk+3,elev_slim)
          jkm1 = MAX(jk-1,slev)
          jkm2 = MAX(jk-2,slev)

          z_delta(jc,jk) = p_face(jc,jk) - p_face(jc,jkp1)        
          z_a6i(jc,jk)   = 6._wp * (p_cc(jc,jk)                      &
            &             - 0.5_wp * (p_face(jc,jk) + p_face(jc,jkp1)))

          ! main criterion upon which it is decided whether an undershoot 
          ! is spurious and whether the limiter is activated.
          is_main_crit = ABS(z_delta(jc,jk)) < ABS(z_a6i(jc,jk))

          ! final (more selective) criterion which decides upon limiter activation
          ! takes into account main criterion is_main_crit 
          l_limit(jc,jk) = isExtremumSpurious(is_main_crit, z_delta(jc,jk), z_a6i(jc,jk), p_cc(jc,jk), &
            &                                 p_face(jc,jkm2), p_face(jc,jkm1), p_face(jc,jk),         & 
            &                                 p_face(jc,jkp1), p_face(jc,jkp2), p_face(jc,jkp3) )
        ENDDO
      ENDDO DETECT_SEL
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(z_delta, z_a6i, l_limit) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jkp1) COLLAPSE(2)
      DETECT:DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! index of bottom half level
          jkp1 = jk+1

          z_delta(jc,jk) = p_face(jc,jk) - p_face(jc,jkp1)        
          z_a6i(jc,jk)   = 6._wp * (p_cc(jc,jk)                      &
            &             - 0.5_wp * (p_face(jc,jk) + p_face(jc,jkp1)))

          ! main criterion upon which it is decided whether an undershoot 
          ! is spurious and whether the limiter is activated. 
          l_limit(jc,jk) = ABS(z_delta(jc,jk)) < ABS(z_a6i(jc,jk))
        ENDDO
      ENDDO DETECT
      !$ACC END PARALLEL
      !
    ENDIF  ! p_ivlimit_selective


    !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(z_delta, z_a6i, l_limit) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(jkp1, q_face_up, q_face_low) COLLAPSE(2)
    LIMIT:DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! index of bottom half level
        jkp1 = jk+1

        !
        ! check if parabola must be modified to remove spurious extrema
        !
        IF ( l_limit(jc,jk) ) THEN

          ! if cell average presents a local extremum, replace parabola 
          ! by piecewise constant function
          IF ( ((p_cc(jc,jk) - p_face(jc,jkp1))*(p_cc(jc,jk)-p_face(jc,jk))) > 0._wp) THEN
            q_face_up  = p_cc(jc,jk)
            q_face_low = p_cc(jc,jk)

          ELSE 
            !
            ! monotonize parabola by modifying one of the edge values
            IF (z_delta(jc,jk) * z_a6i(jc,jk) > z_delta(jc,jk) * z_delta(jc,jk)) THEN
              q_face_up  = p_face(jc,jk)
              q_face_low = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jk)
              
            ELSE IF (z_delta(jc,jk) * z_a6i(jc,jk) < -1._wp * (z_delta(jc,jk) * z_delta(jc,jk))) THEN
              q_face_up  = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jkp1)
              q_face_low = p_face(jc,jkp1)
              
            ELSE
              ! necessary if z_delta and z_a6i become very tiny.
              q_face_up  = p_face(jc,jk)
              q_face_low = p_face(jc,jkp1)
            ENDIF
            !
          ENDIF
          !
        ELSE
          ! no monotonization required
          q_face_up  = p_face(jc,jk)
          q_face_low = p_face(jc,jkp1)
        ENDIF

        ! backcopy of eventually modified face values
        p_face_up(jc,jk)  = q_face_up
        p_face_low(jc,jk) = q_face_low

      END DO  ! jc
    END DO LIMIT
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE v_limit_parabola_mo 



  !-------------------------------------------------------------------------
  !>
  !! Filter for PPM/PSM (3rd order) vertical advection (semi-monotonic version)
  !!
  !! Removes undershoots in first guess parabola by resetting either the
  !! upper or lower interface value.
  !! Avoids non-physical undershoots in advected fields.
  !!
  !! Literature
  !! Lin and Rood (1996), MWR, 124, 2046-2070
  !!
  SUBROUTINE v_limit_parabola_sm( p_ivlimit_selective, p_cc, p_face, &
    &                           p_face_up, p_face_low, i_startidx,   &
    &                           i_endidx, slev, elev     )


    INTEGER, INTENT(IN) ::     &      !< avoids limiting of smooth extrema
      &  p_ivlimit_selective          !< if activated

    REAL(wp), INTENT(IN) ::    &      !< advected cell centered variable
      &  p_cc(:,:)                    !< dim: (nproma,nlev)

    REAL(wp), INTENT(IN) ::    &      !< reconstructed face values of the advected field
      &  p_face(:,:)                  !< dim: (nproma,nlevp1)

    REAL(wp), INTENT(INOUT) :: &      !< final face value (upper face, height based)
      &  p_face_up(:,:)               !< dim: (nproma,nlevp)

    REAL(wp), INTENT(INOUT) :: &      !< final face value (lower face, height based)
      &  p_face_low(:,:)              !< dim: (nproma,nlevp)

    INTEGER, INTENT(IN)     :: &      !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &      !< vertical start and end index of DO loop
      &  slev, elev

    ! local
    REAL(wp) :: z_delta               !< undivided cell gradient
    REAL(wp) :: z_a6i                 !< curvature of parabola
    LOGICAL  :: is_main_crit          !< main criterion for limiting
    LOGICAL  :: l_limit(SIZE(p_cc,1),SIZE(p_cc,2)) !< is limiting of subgrid parabola 
                                                   !< in particular cell necessary [yes/no]
    REAL(wp) :: q_face_up, q_face_low !< face values at upper and lower cell edge

    INTEGER  :: jc, jk                !< index of cell and vertical level
    INTEGER  :: jkm2, jkm1,    &      !< neighbour indices  
      &         jkp1, jkp2, jkp3
    INTEGER  :: elev_slim             !< end level for spurious extremum dectector

  !-------------------------------------------------------------------------


    !$ACC DATA PRESENT(p_cc, p_face, p_face_up, p_face_low) CREATE(l_limit)

    ! selective limitation yes or no
    IF (p_ivlimit_selective == 1) THEN
      ! 
      elev_slim = MIN(elev+1,UBOUND(p_face,2))
      !
      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(l_limit) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jkp1, jkp2, jkp3, jkm1, jkm2, z_delta, z_a6i, is_main_crit) COLLAPSE(2)
      DETECT_SEL:DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! index of bottom half level
          jkp1 = jk+1
          jkp2 = MIN(jk+2,elev_slim)
          jkp3 = MIN(jk+3,elev_slim)
          jkm1 = MAX(jk-1,slev)
          jkm2 = MAX(jk-2,slev)

          z_delta = p_face(jc,jk) - p_face(jc,jkp1)        
          z_a6i   = 6._wp * (p_cc(jc,jk)                      &
            &     - 0.5_wp * (p_face(jc,jk) + p_face(jc,jkp1)))

          ! main criterion upon which it is decided whether an undershoot 
          ! is spurious and whether the limiter is activated.
          is_main_crit = ABS(z_delta) < -1._wp*z_a6i

          ! final (more selective) criterion which decides upon limiter activation
          ! takes into account main criterion is_main_crit 
          l_limit(jc,jk) = isExtremumSpurious(is_main_crit, z_delta, z_a6i, p_cc(jc,jk),        &
            &                                 p_face(jc,jkm2), p_face(jc,jkm1), p_face(jc,jk),  & 
            &                                 p_face(jc,jkp1), p_face(jc,jkp2), p_face(jc,jkp3) )
        ENDDO
      ENDDO DETECT_SEL
      !$ACC END PARALLEL
    ELSE
      !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(l_limit) ASYNC(1)
      !$ACC LOOP GANG VECTOR PRIVATE(jkp1, z_delta, z_a6i) COLLAPSE(2)
      DETECT:DO jk = slev, elev
        DO jc = i_startidx, i_endidx

          ! index of bottom half level
          jkp1 = jk+1

          z_delta = p_face(jc,jk) - p_face(jc,jkp1)        
          z_a6i   = 6._wp * (p_cc(jc,jk)                      &
            &     - 0.5_wp * (p_face(jc,jk) + p_face(jc,jkp1)))

          ! main criterion upon which it is decided whether an undershoot 
          ! is spurious and whether the limiter is activated. 
          l_limit(jc,jk) = ABS(z_delta) < -1._wp*z_a6i
        ENDDO
      ENDDO DETECT
      !$ACC END PARALLEL
      !
    ENDIF  ! p_ivlimit_selective


    !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(l_limit) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(jkp1, q_face_up, q_face_low) COLLAPSE(2)
    LIMIT:DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! index of bottom half level
        jkp1 = jk+1

        !
        ! check if parabola must be modified to remove local undershoots
        !
        IF (l_limit(jc,jk)) THEN

          ! if cell average presents a local minimum, replace parabola 
          ! by piecewise constant function
          IF (p_cc(jc,jk) < MIN(p_face(jc,jk),p_face(jc,jkp1)) ) THEN
            q_face_up  = p_cc(jc,jk)
            q_face_low = p_cc(jc,jk)

          ELSE
            !
            ! monotonize parabola by modifying one of the edge values
            IF (p_face(jc,jk) > p_face(jc,jkp1)) THEN
              q_face_up  = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jkp1)
              q_face_low = p_face(jc,jkp1)

            ELSE
              q_face_up  = p_face(jc,jk)
              q_face_low = 3._wp*p_cc(jc,jk) - 2._wp*p_face(jc,jk)

           ENDIF
              !
          ENDIF
          !
        ELSE
          q_face_up  = p_face(jc,jk)
          q_face_low = p_face(jc,jkp1)
        ENDIF

        ! backcopy of eventually modified face values
        p_face_up(jc,jk)  = q_face_up
        p_face_low(jc,jk) = q_face_low

      END DO
    END DO LIMIT
    !$ACC END PARALLEL

    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE v_limit_parabola_sm



  !-------------------------------------------------------------------------
  !>
  !! Detect if the subgrid reconstruction has a spurious extremum
  !!
  !! Detect if the subgrid reconstruction has an extremum 
  !! in the cell interior, i.e. if 0<\zeta_ext<1. \zeta_ext denotes 
  !! the location of the extremum in dimnesionless coordinates.
  !! This is checked by the main criterion is_main_crit.
  !!
  !! The extremum is deemed spurious, only if at least one of 5 additional 
  !! constraints is fulfilled. This is checked by is_add_crit.
  !!
  !! Literature
  !! Zerroukat et al. (2005), QJRMS, 131, 2923-2936
  !!
  LOGICAL FUNCTION isExtremumSpurious(is_main_crit, z_delta, z_a6i, p_cc, q_face_jkm2,  &
    &                                 q_face_jkm1, q_face_jk, q_face_jkp1, q_face_jkp2, &
    &                                 q_face_jkp3)
    !$ACC ROUTINE SEQ

    LOGICAL,  INTENT(IN) :: is_main_crit     !< is main criterion for limiter activation TRUE
    REAL(wp), INTENT(IN) :: z_delta          !< undivided cell gradient
    REAL(wp), INTENT(IN) :: z_a6i            !< curvature of parabola
    REAL(wp), INTENT(IN) :: p_cc             !< cell average at given cell jk
    REAL(wp), INTENT(IN) ::   &
     & q_face_jkm2, q_face_jkm1, q_face_jk,& !< face values for jk-2, jk-1, jk, jk+1, jk+2, jk+3
     & q_face_jkp1, q_face_jkp2, q_face_jkp3 

    ! local
    LOGICAL  :: is_add_crit         ! additional, more selective criterion
    REAL(wp) :: z_a1                ! linear coefficient of parabolic interpolant
    REAL(wp) :: loc_extremum        ! location of extreme value in dimensionless coordinate
    REAL(wp) :: val_extremum        ! extreme value
    !-----------------------------------------------------------------------


    ! linear coefficient of parabolic interpolant
    z_a1  = 6._wp * p_cc - 2._wp*q_face_jk - 4._wp*q_face_jkp1

    ! determine location of extremum
    loc_extremum = 0.5_wp*(1._wp + z_delta/(SIGN(ABS(z_a6i)+dbl_eps,z_a6i)))
    ! determine extremum
    val_extremum = p_cc - z_delta*(0.5_wp - loc_extremum)  &
      &           - z_a6i*(1._wp/6._wp - loc_extremum + loc_extremum*loc_extremum)


    ! additional, more selective criterion
    ! Only if any of the following constraints is TRUE, the extremum is deemed spurious 
    is_add_crit = ((q_face_jkp1 - q_face_jkp2) * (q_face_jkm1 - q_face_jk)  ) >= 0._wp .OR. &
      &           ((q_face_jkp2 - q_face_jkp3) * (q_face_jkp1 - q_face_jkp2)) <= 0._wp .OR. &
      &           ((q_face_jkm1 - q_face_jk  ) * (q_face_jkm2 - q_face_jkm1)) <= 0._wp .OR. &
      &           ((q_face_jkp1 - q_face_jkp2) * z_a1 )                       <= 0._wp .OR. &
      &           val_extremum                                                <  0._wp

    ! result
    isExtremumSpurious = is_main_crit .AND. is_add_crit

  END FUNCTION isExtremumSpurious



  !-------------------------------------------------------------------------
  !>
  !! Monotonicity preserving slope limiter for PPM after Lin et al (1994)
  !!
  !! Monotonicity preserving slope limiter after Lin et al (1994).
  !! Guarantees e.g. that reconstructed PPM edge values lie in 
  !! the range of values defined by neighbouring cells. 
  !!
  !! Literature
  !! Lin et al. (1994), MWR, 122, 1575-1593
  !!
  SUBROUTINE v_limit_slope_mo( p_cc, i_startidx, i_endidx, slev, elev, slope )


    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)


    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp), INTENT(INOUT) :: &   !< advected cell centered variable
      &  slope(:,:)                !< dim: (nproma,nlev)

    REAL(wp) :: p_cc_min, p_cc_max    !< 3-point max/min values

    INTEGER  :: jc, jk                !< index of cell and vertical level
    INTEGER  :: ikp1, ikm1            !< vertical level plus-minus one

  !-------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(ikm1, ikp1, p_cc_min, p_cc_max) COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, elev )

        ! equivalent formulation of Colella and Woodward (1984) slope limiter 
        ! following Lin et al (1994).
        p_cc_min = MIN(p_cc(jc,ikm1),p_cc(jc,jk),p_cc(jc,ikp1))
        p_cc_max = MAX(p_cc(jc,ikm1),p_cc(jc,jk),p_cc(jc,ikp1))
        slope(jc,jk) = SIGN(                                                   &
          &            MIN( ABS(slope(jc,jk)), 2._wp*(p_cc(jc,jk)-p_cc_min),   &
          &                                    2._wp*(p_cc_max-p_cc(jc,jk)) ), &
          &            slope(jc,jk) )
      END DO  ! jc

    END DO  ! jk 
    !$ACC END PARALLEL

  END SUBROUTINE v_limit_slope_mo



  !-------------------------------------------------------------------------
  !>
  !! Semi-monotonic slope limiter for PPM after Lin et al (1994)
  !!
  !! Semi-monotonic slope limiter after Lin et al (1994).
  !! Guarantees e.g. that reconstructed PPM edge values do not 
  !! fall below the range of values defined by neighbouring cells. 
  !!
  !! Literature
  !! Lin et al. (1994), MWR, 122, 1575-1593
  !!
  SUBROUTINE v_limit_slope_sm( p_cc, i_startidx, i_endidx, slev, elev, slope )


    REAL(wp), INTENT(IN) ::  &     !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma,nlev)


    INTEGER, INTENT(IN)     :: &   !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     :: &   !< vertical start and end index of DO loop
      &  slev, elev

    REAL(wp), INTENT(INOUT) :: &   !< advected cell centered variable
      &  slope(:,:)                !< dim: (nproma,nlev)

    REAL(wp) :: p_cc_min           !< 3-point min values

    INTEGER  :: jc, jk                !< index of cell and vertical level
    INTEGER  :: ikp1, ikm1            !< vertical level plus-minus one

  !-------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(ikm1, ikp1, p_cc_min) COLLAPSE(2)
    DO jk = slev, elev
      DO jc = i_startidx, i_endidx

        ! index of top half level
        ikm1    = jk - 1
        ! index of bottom half level
        ikp1    = MIN( jk+1, elev )

        ! equivalent formulation of Colella and Woodward (1984) slope limiter 
        ! following Lin et al (1994).
        p_cc_min = MIN(p_cc(jc,ikm1),p_cc(jc,jk),p_cc(jc,ikp1))
        slope(jc,jk) = SIGN(                                                   &
          &            MIN( ABS(slope(jc,jk)), 2._wp*(p_cc(jc,jk)-p_cc_min) ), &
          &            slope(jc,jk) )
      END DO  ! jc

    END DO  ! jk 
    !$ACC END PARALLEL

  END SUBROUTINE v_limit_slope_sm



  !-------------------------------------------------------------------------
  !>
  !! Monotonic face value limiter for PSM
  !! Checks if a face value is bounded by the neighbouring cell averages. 
  !! If not, a linear reconstruction based on a monotoniced 
  !! centered-difference (mc) slope is computed for each of the two neighbouring 
  !! cells, in order to estimate a (bounded) edge value.
  !!
  !! Literature
  !! - White et al. (2008), JCP, 227, 7394-7422
  !! - Leveque (2002), Finite Volume Methods for Hyperbolic Problems, 
  !!                  Cambridge University Press
  !!
  SUBROUTINE v_limit_face_mc_mo( p_cc, p_cellhgt_mc_now, p_face, i_startidx, i_endidx, &
    &                            slev, elev )


    REAL(wp), INTENT(IN) ::     &  !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma, nlev)

    REAL(wp), INTENT(IN) ::     &  !< layer thickness at cell center
      &  p_cellhgt_mc_now(:,:)     !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT) ::  &  !< reconstructed face values
      &  p_face(:,:)               !< dim: (nproma, nlevp1)

    INTEGER, INTENT(IN)     ::  &  !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     ::  &  !< vertical start and end index of DO loop
      &  slev, elev


    ! local
    REAL(wp):: mc_slope_u, mc_slope_l !< monotonized central-difference slope 
                                      !< for adjacent upper (u) and lower (l) cell

    REAL(wp):: faceval_u, faceval_l   !< reconstructed face value
                                      !< based on the linear reconstruction for the 
                                      !< adjacent upper (u) and lower (l) cell 
                                
    INTEGER :: jc, jk                 !< loop indices
    INTEGER :: ikm2, ikm1, ikp1       !< vertical level minus two, minus/plus one

    LOGICAL :: l_limit                !< decides whether edge value needs to be limited or not

  !-------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(ikm2, ikm1, ikp1, l_limit, mc_slope_u, mc_slope_l, faceval_u, faceval_l) &
    !$ACC   COLLAPSE(2)
    DO jk= slev, elev
      DO jc = i_startidx, i_endidx

        ikm2 = MAX(jk-2, slev)
        ikm1 = MAX(jk-1, slev)
        ikp1 = MIN(jk+1, elev)


        ! monotonized central-difference slope for cell jk-1
        mc_slope_u  = mc_limiter(                                 &
          &             p_cc_u       = p_cc(jc,ikm2),             &
          &             p_cc_c       = p_cc(jc,ikm1),             &
          &             p_cc_l       = p_cc(jc,jk),               &
          &             cellhgt_mc_u = p_cellhgt_mc_now(jc,ikm2), &
          &             cellhgt_mc_c = p_cellhgt_mc_now(jc,ikm1), &
          &             cellhgt_mc_l = p_cellhgt_mc_now(jc,jk)    )

        ! monotonized central-difference slope for cell jk
        mc_slope_l  = mc_limiter(                                 &
          &             p_cc_u       = p_cc(jc,ikm1),             &
          &             p_cc_c       = p_cc(jc,jk),               &
          &             p_cc_l       = p_cc(jc,ikp1),             &
          &             cellhgt_mc_u = p_cellhgt_mc_now(jc,ikm1), &
          &             cellhgt_mc_c = p_cellhgt_mc_now(jc,jk),   &
          &             cellhgt_mc_l = p_cellhgt_mc_now(jc,ikp1)  )

        ! monotonized, reconstructed face value from adjacent upper cell
        faceval_u = p_cc(jc,ikm1) - 0.5_wp*p_cellhgt_mc_now(jc,ikm1) * mc_slope_u

        ! monotonized, reconstructed face value from adjacent lower cell
        faceval_l = p_cc(jc,jk) + 0.5_wp*p_cellhgt_mc_now(jc,jk) * mc_slope_l

        ! edge value must be limited, if it is not bounded by the neighbouring cell averages
        l_limit = ( (p_cc(jc,ikm1)-p_face(jc,jk))*(p_face(jc,jk)-p_cc(jc,jk)) ) < 0._wp

        !DR: We have replaced the 'IF (l_limit)' condition by a 'MERGE' statement, 
        !    which leads to a speedup of roughly a factor of 4 on the SX AURORA. 
        !    The previous 'IF (l_limit)' implementation encapsulated the entire computation, 
        !    i.e. mc_slope_u/l and faceval_u/l.
        !      
        !    The drawback of the current implementation is that it is performed 
        !    for all cells, no matter if the face value must be limited or not.
        ! 
        ! monotonized face value for jk-1/2
        ! take average of linear reconstructions from adjacent cells
        p_face(jc,jk) = MERGE( 0.5_wp * (faceval_u + faceval_l), p_face(jc,jk), l_limit)

      END DO  ! jc
    END DO  ! jk
    !$ACC END PARALLEL

  END SUBROUTINE v_limit_face_mc_mo



  !-------------------------------------------------------------------------
  !>
  !! Semi-monotonic face value limiter for PSM
  !! Checks if a face value is bounded by the neighbouring cell averages. 
  !! If not, a linear reconstruction based on a monotoniced 
  !! centered-difference (mc) slope is computed for each of the two neighbouring 
  !! cells, in order to estimate a (bounded) edge value.
  !!
  !! Note that this routine only differs from v_limit_face_mc_mo wrt the 
  !! detection criterion l_limit. In principle it would be possible to merge 
  !! both routines with the help of procedure pointers. However, this would 
  !! introduce some performance penalty, as the procedure pointer hinders 
  !! inlining of the detection citerion.
  !!
  !! Literature
  !! - White et al. (2008), JCP, 227, 7394-7422
  !! - Leveque (2002), Finite Volume Methods for Hyperbolic Problems, 
  !!                  Cambridge University Press
  !!
  SUBROUTINE v_limit_face_mc_sm( p_cc, p_cellhgt_mc_now, p_face, i_startidx, i_endidx, &
    &                            slev, elev )


    REAL(wp), INTENT(IN) ::     &  !< advected cell centered variable
      &  p_cc(:,:)                 !< dim: (nproma, nlev)

    REAL(wp), INTENT(IN) ::     &  !< layer thickness at cell center
      &  p_cellhgt_mc_now(:,:)     !< dim: (nproma,nlev)

    REAL(wp), INTENT(INOUT) ::  &  !< reconstructed face values
      &  p_face(:,:)               !< dim: (nproma, nlevp1)

    INTEGER, INTENT(IN)     ::  &  !< horizontal start and end index of DO loop
      &  i_startidx, i_endidx

    INTEGER, INTENT(IN)     ::  &  !< vertical start and end index of DO loop
      &  slev, elev


    ! local
    REAL(wp):: mc_slope_u, mc_slope_l !< monotonized central-difference slope 
                                      !< for adjacent upper (u) and lower (l) cell

    REAL(wp):: faceval_u, faceval_l   !< reconstructed face value
                                      !< based on the linear reconstruction for the 
                                      !< adjacent upper (u) and lower (l) cell 
                                
    INTEGER :: jc, jk                 !< loop indices
    INTEGER :: ikm2, ikm1, ikp1       !< vertical level minus two, minus/plus one

    LOGICAL :: l_limit                !< decides whether edge value needs to be limited or not

  !-------------------------------------------------------------------------

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(ikm2, ikm1, ikp1, l_limit, mc_slope_u, mc_slope_l, faceval_u, faceval_l) &
    !$ACC   COLLAPSE(2)
    DO jk= slev, elev
      DO jc = i_startidx, i_endidx

        ikm2 = MAX(jk-2, slev)
        ikm1 = MAX(jk-1, slev)
        ikp1 = MIN(jk+1, elev)


        ! monotonized central-difference slope for cell jk-1
        mc_slope_u  = mc_limiter(                                 &
          &             p_cc_u       = p_cc(jc,ikm2),             &
          &             p_cc_c       = p_cc(jc,ikm1),             &
          &             p_cc_l       = p_cc(jc,jk),               &
          &             cellhgt_mc_u = p_cellhgt_mc_now(jc,ikm2), &
          &             cellhgt_mc_c = p_cellhgt_mc_now(jc,ikm1), &
          &             cellhgt_mc_l = p_cellhgt_mc_now(jc,jk)    )

        ! monotonized central-difference slope for cell jk
        mc_slope_l  = mc_limiter(                                 &
          &             p_cc_u       = p_cc(jc,ikm1),             &
          &             p_cc_c       = p_cc(jc,jk),               &
          &             p_cc_l       = p_cc(jc,ikp1),             &
          &             cellhgt_mc_u = p_cellhgt_mc_now(jc,ikm1), &
          &             cellhgt_mc_c = p_cellhgt_mc_now(jc,jk),   &
          &             cellhgt_mc_l = p_cellhgt_mc_now(jc,ikp1)  )



        ! monotonized, reconstructed face value from adjacent upper cell
        faceval_u = p_cc(jc,ikm1) - 0.5_wp*p_cellhgt_mc_now(jc,ikm1) * mc_slope_u

        ! monotonized, reconstructed face value from adjacent lower cell
        faceval_l = p_cc(jc,jk) + 0.5_wp*p_cellhgt_mc_now(jc,jk) * mc_slope_l

        ! edge value must be limited, if it is not bounded by the neighbouring cell averages
        l_limit = p_face(jc,jk) < MIN(p_cc(jc,ikm1),p_cc(jc,jk))

        !DR: We have replaced the 'IF (l_limit)' condition by a 'MERGE' statement, 
        !    which leads to a speedup of roughly a factor of 4 on the SX AURORA. 
        !    The previous 'IF (l_limit)' implementation encapsulated the entire computation, 
        !    i.e. mc_slope_u/l and faceval_u/l.
        !      
        !    The drawback of the current implementation is that it is performed 
        !    for all cells, no matter if the face value must be limited or not.
        ! 
        ! monotonized face value for jk-1/2
        ! take average of linear reconstructions from adjacent cells
        p_face(jc,jk) = MERGE( 0.5_wp * (faceval_u + faceval_l), p_face(jc,jk), l_limit)

      END DO  ! jc
    END DO  ! jk
    !$ACC END PARALLEL

  END SUBROUTINE v_limit_face_mc_sm



  !-------------------------------------------------------------------------
  !>
  !! van Leer-type monotonized central-difference slope
  !!
  !! Literature
  !! - White et al. (2008), JCP, 227, 7394-7422
  !! - Leveque (2002), Finite Volume Methods for Hyperbolic Problems, 
  !!                 Cambridge University Press
  !!
  FUNCTION mc_limiter (p_cc_u, p_cc_c, p_cc_l, cellhgt_mc_u, cellhgt_mc_c, cellhgt_mc_l) &
    &                  RESULT(mc_slope)
    !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN)  :: p_cc_u, p_cc_c, p_cc_l  ! advected variable for 
                                                     ! upper (u) center (c) and lower (l) cell
    REAL(wp), INTENT(IN)  :: cellhgt_mc_u, cellhgt_mc_c, cellhgt_mc_l

    ! local
    REAL(wp) :: slope_u, slope_l   ! one-sided slopes
    REAL(wp) :: slope_c            ! centered-difference slope

    ! monotonized central-difference slope
    REAL(wp) :: mc_slope
  !-----------------------------------------------

    ! one-sided upper slope
    slope_u = 2._wp*(p_cc_u - p_cc_c)/cellhgt_mc_c
    ! one-sided lower slope
    slope_l = 2._wp*(p_cc_c - p_cc_l)/cellhgt_mc_c
    ! centered difference slope
    slope_c = 2._wp*(p_cc_u - p_cc_l)  &
      &     /(cellhgt_mc_u + 2._wp*cellhgt_mc_c + cellhgt_mc_l)

    IF ( slope_u*slope_l > 0._wp) THEN
      mc_slope = SIGN(MIN(ABS(slope_u),ABS(slope_l),ABS(slope_c)),slope_c)
    ELSE
      mc_slope = 0._wp
    ENDIF

  END FUNCTION mc_limiter

END MODULE mo_advection_vlimit

