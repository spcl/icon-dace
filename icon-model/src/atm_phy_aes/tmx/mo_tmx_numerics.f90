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

MODULE mo_tmx_numerics

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_fortran_tools,     ONLY: init
  USE mo_surrogate_class,   ONLY: t_surrogate
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_time_integration_class, ONLY: t_time_scheme
  USE mo_model_domain,      ONLY: t_patch
  ! USE mo_variable_list, ONLY: t_variable_list

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
    & t_time_scheme_explicit_euler, &
    & diffuse_vertical_explicit, diffuse_vertical_implicit,                     &
    & get_normal_velocity_vertex, get_tangential_velocity_vertex,               &
    & vertical_interpolation_scalar_cell, vertical_interpolation_scalar_vertex

  TYPE, EXTENDS(t_time_scheme) :: t_time_scheme_explicit_euler
  CONTAINS
    PROCEDURE, NOPASS :: Step_forward => step_forward_explicit_euler
  END TYPE t_time_scheme_explicit_euler

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_numerics'
  
CONTAINS
  !============================================================================
  SUBROUTINE step_forward_explicit_euler(process, dt)
    CLASS(t_surrogate), INTENT(inout) :: process
    REAL(wp), INTENT(in) :: dt

    CHARACTER(len=*), PARAMETER :: routine = modname//':step_forward_explicit_euler'

    SELECT TYPE (process)
    CLASS IS (t_tmx_process)
      CALL process%Compute()
      !TBD:
      ! For all fields/tendencies: process[field] = process[field] + process[tendency] * dt
    CLASS DEFAULT
      CALL finish(routine, 'Unkown class for process')
    END SELECT

  END SUBROUTINE step_forward_explicit_euler
  !============================================================================
  ! Explicit vertical diffusion of any physical
  ! quantity. The coefficients of the system of
  ! equations are set in prepare_diffusion_matrix
  ! in mo_vdf_atmo. They are the same as for the
  ! implicit treatment.
  !============================================================================
  SUBROUTINE diffuse_vertical_explicit( &
    & ics, ice,       & ! in
    & minlvl, maxlvl, & ! in
    & a, b, c, rhs,   & ! in
    & var,            & ! in
    & tend            & ! out
    & )

    ! Iteration boundaries for blocks, cells, and level
    INTEGER, INTENT(in) :: ics, ice, minlvl, maxlvl

    ! Solve system of equations of shape
    ! a*x_(k-1) + b*x_(k) + c*x_(k+1) = d,
    ! where x is the variable at time step t.
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & var, &
      & a,   &
      & b,   &
      & c,   &
      & rhs

    ! Tendency of the vertically diffused physical quantity
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: tend

    ! Loop iterators
    INTEGER  :: jc, jk

    CHARACTER(len=*), PARAMETER :: routine = modname//':diffuse_vertical_explicit'

    ! Calculate the tendency explicitly.
    ! Caution, the sign of the coefficients a, b, c must be negative, since they were
    ! defined in prepare_diffusion_matrix for an implicit scheme. In the explicit scheme
    ! they are transferred on the other side of the equation, so their sign changed.
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
    DO jk=minlvl+1,maxlvl-1
      DO jc = ics, ice
        tend(jc,jk) = tend(jc,jk) -              &
                       a(jc,jk) * var(jc,jk-1) - &
                       b(jc,jk) * var(jc,jk) -   &
                       c(jc,jk) * var(jc,jk+1) + &
                       rhs(jc,jk)
      END DO
    END DO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = ics, ice
      tend(jc,minlvl) = tend(jc,minlvl) -                  &
                         b(jc,minlvl) * var(jc,minlvl) -   &
                         c(jc,minlvl) * var(jc,minlvl+1) + &
                         rhs(jc,minlvl)
    END DO

    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = ics, ice
      tend(jc,maxlvl) = tend(jc,maxlvl) -                  &
                         a(jc,maxlvl) * var(jc,maxlvl-1) - &
                         b(jc,maxlvl) * var(jc,maxlvl) +   &
                         rhs(jc,maxlvl)
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT(1)

  END SUBROUTINE diffuse_vertical_explicit
  !============================================================================
  ! Implicit vertical diffusion of any physical
  ! quantity. The coefficients of the system of
  ! equations are set in prepare_diffusion_matrix
  ! in mo_vdf_atmo.
  !============================================================================
  SUBROUTINE diffuse_vertical_implicit( &
    & ics, ice,       & ! in
    & minlvl, maxlvl, & ! in
    & a, c,           & ! in
    & bb, rhs,        & ! in
    & rdtime,         & ! in
    & var,            & ! in
    & tend            & ! inout
    & )

    USE mo_math_utilities, ONLY: tdma_solver_vec

    ! Iteration boundaries for blocks, cells, and level
    INTEGER, INTENT(in) :: ics, ice, minlvl, maxlvl

    ! Solve system of equations of shape
    ! a*x_(k-1) + b*x_(k) + c*x_(k+1) = d,
    ! where x is the variable at time step t+dt.
    ! At this point, bb and rhs still lack a term
    ! depending on the time increment (see module mo_vdf_atmo;
    ! subroutine prepare_diffusion_matrix). Hence,
    ! this term is added to bb and rhs to get the
    ! correct coefficients b and d. 
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & var, &
      & a,   &
      & bb,  &
      & c,   &
      & rhs

    ! Reziprocal time increment
    REAL(wp), INTENT(in) :: rdtime

    ! Tendendy of the vertically diffused physical quantity.
    REAL(wp), INTENT(inout), DIMENSION(:,:) :: tend

    ! Loop iterators
    INTEGER  :: jc, jk

    ! Variable value at time step t+dt
    REAL(wp) :: new_var(SIZE(var,1),SIZE(var,2))

    ! Right hand side of the system of the equation
    REAL(wp) :: d(SIZE(var,1),SIZE(var,2))

    ! Coefficient b of the matrix including the time increment
    REAL(wp) :: b(SIZE(var,1),SIZE(var,2))

    CHARACTER(len=*), PARAMETER :: routine = modname//':diffuse_vertical_implicit'

    !$ACC DATA CREATE(new_var, b, d)

    ! The reciprocal time increment is added at this point to b as it could not be
    ! added earlier, for example in prepare_diffusion_matrix. The same applies for
    ! the product of the physical quantity var and the reciprocal time increment.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
    DO jk=minlvl,maxlvl
      DO jc = ics, ice
        b(jc,jk) = rdtime + bb(jc,jk) 
        d(jc,jk) = var(jc,jk) * rdtime + rhs(jc,jk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !ACC WAIT(1)

    ! Solve the system of equations.
    CALL tdma_solver_vec(a(:,:), b(:,:), c(:,:), d(:,:), &
                      minlvl, maxlvl, ics, ice, new_var(:,:))

    ! Calculate the tendency.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
    DO jk=minlvl,maxlvl
      DO jc = ics, ice
        tend(jc,jk) = tend(jc,jk) + (new_var(jc,jk) - var(jc,jk)) * rdtime
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    !$ACC END DATA

  END SUBROUTINE diffuse_vertical_implicit
  !============================================================================
  ! Determines horizontal velocity component at vertex. 
  ! jv: vertex index of target edge following the numbering in figure 1 in Zaengl et al. 
  ! 2015, Q. J. R. Meteorol. Soc..
  !============================================================================
  FUNCTION get_normal_velocity_vertex(                         &
    u_vert, v_vert, patch, je, jb, jk, jv                           &
    ) RESULT(vn_vert)

    REAL(wp), INTENT(in), POINTER :: u_vert(:,:,:), v_vert(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch

    INTEGER, INTENT(in) :: je, jb, jk, jv
    REAL(wp) :: vn_vert
   
    vn_vert =   u_vert(patch%edges%vertex_idx(je,jb,jv),jk,patch%edges%vertex_blk(je,jb,jv)) &    
                * patch%edges%primal_normal_vert(je,jb,jv)%v1                                &
              + v_vert(patch%edges%vertex_idx(je,jb,jv),jk,patch%edges%vertex_blk(je,jb,jv)) &    
                * patch%edges%primal_normal_vert(je,jb,jv)%v2

  END FUNCTION get_normal_velocity_vertex
  !============================================================================
  ! Determines tangential velocity component at vertex. 
  ! jv: vertex index of target edge following the numbering in figure 1 in Zaengl et al. 
  ! 2015, Q. J. R. Meteorol. Soc..
  !============================================================================
  FUNCTION get_tangential_velocity_vertex(                          &
    u_vert, v_vert, patch, je, jb, jk, jv                           &
    ) RESULT(vt_vert)

    REAL(wp), INTENT(in), POINTER :: u_vert(:,:,:), v_vert(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch

    INTEGER, INTENT(in) :: je, jb, jk, jv
    REAL(wp) :: vt_vert
   
    vt_vert =   u_vert(patch%edges%vertex_idx(je,jb,jv),jk,patch%edges%vertex_blk(je,jb,jv)) &    
                * patch%edges%dual_normal_vert(je,jb,jv)%v1                                  &
              + v_vert(patch%edges%vertex_idx(je,jb,jv),jk,patch%edges%vertex_blk(je,jb,jv)) &    
                * patch%edges%dual_normal_vert(je,jb,jv)%v2

  END FUNCTION get_tangential_velocity_vertex 
  !============================================================================
  ! Local vertical interpolation of pointer scalar variable at cell center. 
  ! Interpolates between jk and jk+1. 
  ! jc: cell index to interpolate (cell numbering see figure A1 in Zaengl et al. 
  ! 2015, Q. J. R. Meteorol. Soc.).
  !============================================================================
  FUNCTION vertical_interpolation_scalar_cell(                      &
    pc,patch,je,jb,jk,jc                                            &
    ) RESULT(pcint)

    REAL(wp), INTENT(in), POINTER :: pc(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch

    INTEGER, INTENT(in) :: je, jb, jk, jc
    REAL(wp) :: pcint
   
    pcint = 0.5_wp * ( & 
                pc(patch%edges%cell_idx(je,jb,jc),jk,  patch%edges%cell_blk(je,jb,jc)) &   
              + pc(patch%edges%cell_idx(je,jb,jc),jk+1,patch%edges%cell_blk(je,jb,jc)) &  
              )

  END FUNCTION vertical_interpolation_scalar_cell 
  !============================================================================
  ! Local vertical interpolation of pointer scalar variable at vertex. 
  ! Interpolates between jk and jk+1. 
  ! jv: vertex index to interpolate (vertex numbering see figure 1 in Zaengl et al. 
  ! 2015, Q. J. R. Meteorol. Soc.).
  !============================================================================
  FUNCTION vertical_interpolation_scalar_vertex(                      &
    pv,patch,je,jb,jk,jv                                            &
    ) RESULT(pvint)

    REAL(wp), INTENT(in), POINTER :: pv(:,:,:)
    TYPE(t_patch), INTENT(in), POINTER :: patch

    INTEGER, INTENT(in) :: je, jb, jk, jv
    REAL(wp) :: pvint
   
    pvint = 0.5_wp * ( & 
                pv(patch%edges%vertex_idx(je,jb,jv),jk,  patch%edges%vertex_blk(je,jb,jv)) &   
              + pv(patch%edges%vertex_idx(je,jb,jv),jk+1,patch%edges%vertex_blk(je,jb,jv)) &  
              )

  END FUNCTION vertical_interpolation_scalar_vertex
  !============================================================================
END MODULE mo_tmx_numerics
