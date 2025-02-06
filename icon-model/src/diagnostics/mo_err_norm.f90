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

! Calculate error norms as in Williamson et al. (1992)
!
! General subroutines to be used for any variable (scalar, or
! vector with 2 or 3 components) and any test case

MODULE mo_err_norm
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!

USE mo_kind,               ONLY: wp

IMPLICIT NONE

PRIVATE

PUBLIC :: t_errors              ! errors TYPE definition
PUBLIC :: err_norm            ! calculate all errors
PUBLIC :: err_norm_scal
PUBLIC :: err_norm_vec2
PUBLIC :: err_norm_vec3
PRIVATE:: compute_err_norm

INTERFACE err_norm
  MODULE PROCEDURE err_norm_scal
  MODULE PROCEDURE err_norm_vec2
  MODULE PROCEDURE err_norm_vec3
END INTERFACE

! error type

TYPE t_errors
  REAL(wp):: abs_l1,       & ! absolute l1 error
    &           abs_l2,       & ! absolute l2 error
    &           abs_linf,     & ! absolute linf error
    &           rel_l1,       & ! relative l1 error
    &           rel_l2,       & ! relative l2 error
    &           rel_linf        ! relative linf error
END TYPE t_errors

CONTAINS

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!
!! Calculation of the normalized errors as in Willianson et al.
!!
!! (1992)
!! for a scalar field.
!! Calculate the "absolute difference field" and the "absolute reference
!! field" for a scalar field and then call compute_err_norm
!! Absolute difference field is the square root of the difference field to the
!! square. Absolute reference field is the square root of the reference field
!! to the square.
!!
SUBROUTINE err_norm_scal (k_npoints, p_xref, p_xmod, p_wgt, p_errors)
!
INTEGER,   INTENT(IN)    :: k_npoints                    ! number of grid pts
REAL(wp),  INTENT(IN)    :: p_xref( k_npoints), & ! reference field
  & p_xmod( k_npoints), & ! model output
  & p_wgt ( k_npoints)    ! integration weigths

TYPE(t_errors), INTENT(OUT) :: p_errors                     ! all errors

REAL(wp)                  :: z_diff( k_npoints), & ! mod.(diff. field)
  & z_reff( k_npoints)    ! mod.(ref. field)

!-----------------------------------------------------------------------

!$OMP PARALLEL WORKSHARE

z_diff = ABS ( p_xmod - p_xref )
z_reff = ABS ( p_xref )

!$OMP END PARALLEL WORKSHARE

CALL compute_err_norm (k_npoints, z_diff, z_reff, p_wgt, p_errors)

END SUBROUTINE  err_norm_scal

!-------------------------------------------------------------------------
!
!! Calculation of the normalized errors as in Williamson et al.
!!
!! (1992)
!! for a vector field with two components..
!! Calculate the "absolute difference field" and the "absolute reference field"
!! for a vector field and then call to compute_err_norm.
!! Absolute difference field is the square root of the difference vector field
!! to the square; absolute reference field is the square root of the reference
!! vector field to the square.
!!
SUBROUTINE err_norm_vec2 (k_npoints, p_uref, p_vref, p_umod, p_vmod, p_wgt,  &
  & p_errors)
!
! !INPUT PARAMETERS
INTEGER,   INTENT(IN)    :: k_npoints                    ! number of grid pts
REAL(wp),  INTENT(IN)    :: p_uref( k_npoints), & ! reference field...
  & p_vref( k_npoints), & ! ...components
  & p_umod( k_npoints), & ! model output...
  & p_vmod( k_npoints), & ! ...components
  & p_wgt ( k_npoints)    ! weigths

! !OUTPUT PARAMETERS
TYPE(t_errors), INTENT(OUT) :: p_errors                     ! all errors

REAL(wp)                  :: z_diff( k_npoints), & ! mod.(diff. field)
  & z_reff( k_npoints)    ! mod.(ref. field)

!-----------------------------------------------------------------------

!$OMP PARALLEL WORKSHARE

z_diff = SQRT ((p_umod - p_uref) * (p_umod - p_uref) +  &
  (p_vmod - p_vref) * (p_vmod - p_vref))

z_reff = SQRT ( p_uref * p_uref + p_vref * p_vref )

!$OMP END PARALLEL WORKSHARE

CALL compute_err_norm (k_npoints, z_diff, z_reff, p_wgt, p_errors)

END SUBROUTINE err_norm_vec2

!-------------------------------------------------------------------------
!
!! Calculation of the normalized errors as in Williamson et al.
!!
!! (1992)
!! for a vector field with three components..
!! Calculate the "absolute difference field" and the "absolute reference field"
!! for a vector field and then call to compute_err_norm.
!! Absolute difference field is the square root of the difference vector field
!! to the square; absolute reference field is the square root of the reference
!! vector field to the square.
!!
SUBROUTINE err_norm_vec3 (k_npoints, p_uref, p_vref, p_wref, p_umod, p_vmod, &
  p_wmod, p_wgt, p_errors)
!
! !INPUT PARAMETERS
INTEGER,   INTENT(IN)    :: k_npoints                    ! number of grid pts
REAL(wp),  INTENT(IN)    :: p_uref( k_npoints), & ! reference...
  & p_vref( k_npoints), & ! ...field...
  & p_wref( k_npoints), & ! ...components
  & p_umod( k_npoints), & ! model...
  & p_vmod( k_npoints), & ! ...output...
  & p_wmod( k_npoints), & ! ...components
  & p_wgt ( k_npoints)    ! weigths

! !OUTPUT PARAMETERS
TYPE(t_errors), INTENT(OUT) :: p_errors                     ! all errors

REAL(wp)                  :: z_diff( k_npoints), & ! mod.(diff. field)
  & z_reff( k_npoints)    ! mod.(ref. field)

!-----------------------------------------------------------------------

!$OMP PARALLEL WORKSHARE

z_diff = (p_umod - p_uref) * (p_umod - p_uref) +  &
  (p_vmod - p_vref) * (p_vmod - p_vref) +  &
  (p_wmod - p_wref) * (p_wmod - p_wref)
z_diff = SQRT (z_diff)

z_reff = SQRT ( p_uref * p_uref + p_vref * p_vref + p_wref * p_wref)

!$OMP END PARALLEL WORKSHARE

CALL compute_err_norm (k_npoints, z_diff, z_reff, p_wgt, p_errors)

END SUBROUTINE err_norm_vec3

!-------------------------------------------------------------------------
!
!! Calculate the normalized errors as in Williamson et al.
!!
!! (1992).
!! It needs the modulus of the difference field and the modulus of the reference
!! field, the number of grid points and the weights (areas) for each grid point.
!! In case the reference field is more or less zero, the relative error is
!! defined as absolute error.
!!
SUBROUTINE compute_err_norm (k_npoints, p_diff, p_reff, p_wgt, p_errors)
!
REAL(wp), PARAMETER       :: z_tol = 1.e-10_wp            ! smallest value...
! ...for modulus of reference field
! !INPUT PARAMETERS
INTEGER, INTENT(IN)       :: k_npoints                    ! number of grid pts
REAL(wp), INTENT(IN)      :: p_diff( k_npoints), & ! mod. difference
  & p_reff( k_npoints), & ! mod. reference
  & p_wgt ( k_npoints)    ! weights
! !OUTPUT PARAMETERS
TYPE(t_errors), INTENT(OUT) :: p_errors                     ! all errors

REAL(wp)                  :: z_diff2( k_npoints),& ! square(mod.diff.)
  & z_reff2( k_npoints)   ! square(mod.ref.)
REAL(wp)                  :: z_den                        ! denominator in...
! ...error calculation
REAL(wp)                  :: z_err                        ! error
REAL(wp)                  :: z_sum_wgt                    ! sum of all weights

!-----------------------------------------------------------------------

!$OMP PARALLEL

!$OMP WORKSHARE

z_diff2 = p_diff * p_diff
z_reff2 = p_reff * p_reff

z_err     = SUM( p_diff * p_wgt )
z_den     = SUM( p_reff * p_wgt )
z_sum_wgt = SUM( p_wgt )

!$OMP END WORKSHARE

!$OMP SINGLE

p_errors%abs_l1 = z_err / z_sum_wgt

IF ( z_den > z_tol) THEN
  p_errors%rel_l1 = z_err / z_den
ELSE
  p_errors%rel_l1 = p_errors%abs_l1
ENDIF

!$OMP END SINGLE

!$OMP WORKSHARE

z_err     = SQRT( SUM(z_diff2 * p_wgt) )
z_den     = SQRT( SUM(z_reff2 * p_wgt) )
z_sum_wgt = SQRT( z_sum_wgt )

!$OMP END WORKSHARE

!$OMP SINGLE

p_errors%abs_l2 = z_err / z_sum_wgt

IF ( z_den > z_tol) THEN
  p_errors%rel_l2 = z_err / z_den
ELSE
  p_errors%rel_l2 = p_errors%abs_l2
ENDIF

!$OMP END SINGLE

!$OMP WORKSHARE

z_err = MAXVAL(p_diff)
z_den = MAXVAL(p_reff)

!$OMP END WORKSHARE

!$OMP SINGLE

p_errors%abs_linf = z_err

IF ( z_den > z_tol) THEN
  p_errors%rel_linf = z_err / z_den
ELSE
  p_errors%rel_linf = z_err
ENDIF

!$OMP END SINGLE

!$OMP END PARALLEL

END SUBROUTINE  compute_err_norm

END MODULE mo_err_norm

