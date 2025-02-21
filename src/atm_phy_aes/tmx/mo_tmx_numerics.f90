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

MODULE mo_tmx_numerics

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: finish
  USE mo_fortran_tools,     ONLY: init
  USE mo_surrogate_class,   ONLY: t_surrogate
  USE mo_tmx_process_class, ONLY: t_tmx_process
  USE mo_tmx_time_integration_class, ONLY: t_time_scheme
  ! USE mo_variable_list, ONLY: t_variable_list

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: &
    & t_time_scheme_explicit_euler, &
    & diffuse_scalar_vertical_explicit, diffuse_scalar_vertical_implicit

  TYPE, EXTENDS(t_time_scheme) :: t_time_scheme_explicit_euler
  CONTAINS
    PROCEDURE, NOPASS :: Step_forward => step_forward_explicit_euler
  END TYPE t_time_scheme_explicit_euler

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_numerics'
  
CONTAINS

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

  SUBROUTINE diffuse_scalar_vertical_explicit( &
    & ibs, ibe, ics, ice, &
    & dz, zf, &
    & rho_ic, &
    & k_ic, &
    & var, &
    & sfc_flx, &
    & top_flx, &
    & tend &
    & )

    INTEGER, INTENT(in) :: &
      & ibs, ibe, ics(:), ice(:)

    REAL(wp), INTENT(in), DIMENSION(:,:,:) :: &
      & var, &
      & rho_ic, &
      & k_ic, &
      & dz, &
      & zf

    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & sfc_flx, &
      & top_flx

    REAL(wp), INTENT(out), DIMENSION(:,:,:) :: &
      & tend

    INTEGER :: nlev
    INTEGER :: jb, jc, jk

    CHARACTER(len=*), PARAMETER :: routine = modname//':diffuse_scalar_vertical_explicit'

    nlev = SIZE(var,2)

!$OMP PARALLEL
    CALL init(tend)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = ibs, ibe

      DO jk=2,nlev-1
        DO jc = ics(jb), ice(jb)
          tend(jc,jk,jb) =  1._wp/dz(jc,jk,jb) * &
          & ( &
          &   rho_ic(jc,jk  ,jb) * k_ic(jc,jk  ,jb) * (var(jc,jk-1,jb)-var(jc,jk  ,jb)) / (zf(jc,jk-1,jb)-zf(jc,jk  ,jb))&
          & - rho_ic(jc,jk+1,jb) * k_ic(jc,jk+1,jb) * (var(jc,jk  ,jb)-var(jc,jk+1,jb)) / (zf(jc,jk  ,jb)-zf(jc,jk+1,jb))&
          & )
        END DO
      END DO

      jk=1
      DO jc = ics(jb), ice(jb)
        tend(jc,jk,jb) =  1._wp/dz(jc,jk,jb) * &
        & ( &
        &                      top_flx(jc,jb) &
        & - rho_ic(jc,jk+1,jb) * k_ic(jc,jk+1,jb) * (var(jc,jk  ,jb)-var(jc,jk+1,jb)) / (zf(jc,jk  ,jb)-zf(jc,jk+1,jb)) &
        & )
      END DO

      jk=nlev
      DO jc = ics(jb), ice(jb)
        tend(jc,jk,jb) =  1._wp/dz(jc,jk,jb) * &
        & ( &
        &   rho_ic(jc,jk  ,jb) * k_ic(jc,jk,  jb) * (var(jc,jk-1,jb)-var(jc,jk  ,jb)) / (zf(jc,jk-1,jb)-zf(jc,jk  ,jb)) &
        & - sfc_flx(jc,jb) &
        & )
      END DO

    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE diffuse_scalar_vertical_explicit

  SUBROUTINE diffuse_scalar_vertical_implicit( &
    & ibs, ibe, ics, ice, &
    & dtime, &
    & dz, zf, &
    & rho_ic, &
    & k_ic, &
    & var, &
    & sfc_flx, &
    & top_flx, &
    & tend &
    & )

    USE mo_math_utilities, ONLY: tdma_solver_vec

    INTEGER, INTENT(in) :: &
      & ibs, ibe, ics(:), ice(:)

    REAL(wp), INTENT(in) :: &
      & dtime

    REAL(wp), INTENT(in), DIMENSION(:,:,:) :: &
      & var, &
      & rho_ic, &
      & k_ic, &
      & dz, &
      & zf

    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & sfc_flx, &
      & top_flx

    REAL(wp), INTENT(out), DIMENSION(:,:,:) :: &
      & tend

    INTEGER  :: nlev
    INTEGER  :: jb, jc, jk
    REAL(wp) :: &
      & new_var(SIZE(var,1),SIZE(var,2),SIZE(var,3)), &
      & a      (SIZE(var,1),SIZE(var,2),SIZE(var,3)), &
      & b      (SIZE(var,1),SIZE(var,2),SIZE(var,3)), &
      & c      (SIZE(var,1),SIZE(var,2),SIZE(var,3)), &
      & rhs    (SIZE(var,1),SIZE(var,2),SIZE(var,3))

    CHARACTER(len=*), PARAMETER :: routine = modname//':diffuse_scalar_vertical_implicit'

    nlev = SIZE(var,2)

    !$ACC DATA CREATE(a, b, c, rhs, new_var)

!$OMP PARALLEL
    CALL init(new_var)
    CALL init(tend)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(jb,jc,jk) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = ibs, ibe

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO jk=2,nlev-1
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = ics(jb), ice(jb)
          a(jc,jk,jb) = - k_ic(jc,jk  ,jb) / dz(jc,jk,jb) * rho_ic(jc,jk  ,jb) / (zf(jc,jk-1,jb)-zf(jc,jk  ,jb))
          c(jc,jk,jb) = - k_ic(jc,jk+1,jb) / dz(jc,jk,jb) * rho_ic(jc,jk+1,jb) / (zf(jc,jk  ,jb)-zf(jc,jk+1,jb))
          b(jc,jk,jb) = 1._wp / dtime - a(jc,jk,jb) - c(jc,jk,jb)
          rhs(jc,jk,jb) = var(jc,jk,jb) / dtime
        END DO
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = ics(jb), ice(jb)
        a(jc,1,jb) = 0._wp 
        c(jc,1,jb) = - k_ic(jc,2,jb) / dz(jc,1,jb) * rho_ic(jc,2,jb) / (zf(jc,1,jb)-zf(jc,2,jb))
        b(jc,1,jb) = 1._wp / dtime - a(jc,1,jb) - c(jc,1,jb)
        rhs(jc,1,jb) = var(jc,1,jb) / dtime  ! TODO: only correct for top_flx=0 !
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = ics(jb), ice(jb)
        a(jc,nlev,jb) = - k_ic(jc,nlev  ,jb) / dz(jc,nlev,jb) * rho_ic(jc,nlev  ,jb) / (zf(jc,nlev-1,jb)-zf(jc,nlev,jb))
        c(jc,nlev,jb) = 0._wp 
        b(jc,nlev,jb) = 1._wp / dtime - a(jc,nlev,jb) - c(jc,nlev,jb)
        rhs(jc,nlev,jb) = var(jc,nlev,jb) / dtime - sfc_flx(jc,jb) / dz(jc,nlev,jb) 
      END DO

      !$ACC END PARALLEL

      ! DO jc = ics(jb), ice(jb)
      !   CALL tdma_solver(a(jc,:,jb), b(jc,:,jb), c(jc,:,jb), rhs(jc,:,jb), &
      !                   nlev, new_var(jc,:,jb))
      CALL tdma_solver_vec(a(:,:,jb), b(:,:,jb), c(:,:,jb), rhs(:,:,jb), &
                        1, nlev, ics(jb), ice(jb), new_var(:,:,jb))
      ! END DO

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG(STATIC: 1) VECTOR COLLAPSE(2) ASYNC(1)
      DO jk=1,nlev
        DO jc = ics(jb), ice(jb)
          tend(jc,jk,jb) = (new_var(jc,jk,jb)-var(jc,jk,jb)) / dtime
        END DO
      END DO
      !$ACC END PARALLEL LOOP

    ENDDO
!$OMP END PARALLEL DO

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE diffuse_scalar_vertical_implicit

END MODULE mo_tmx_numerics
