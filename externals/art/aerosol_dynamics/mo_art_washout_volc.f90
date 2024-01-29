!
! mo_art_washout_volc
! This module provides washout for one-moment aerosol
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

MODULE mo_art_washout_volc
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  PUBLIC   :: art_washout_volc 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_washout_volc(dtime, istart, iend, kstart, nlev, p_qr, prr_gsp, prr_con, prr_con3d, tracer, lacc)
!<
! SUBROUTINE art_washout_volc
! This subroutine provides washout for one-moment aerosol, for the derivation of the equations
! used see Rieger et al. (2015)
! Part of Module: mo_art_washout_volc
! Author: Kristina Lundgren, KIT
! Initial Release: 2012-06-15
! Modifications:
!! 2016-09-19: Daniel Rieger, KIT
!! - Passing of 2D fields by argument instead of types containing 3D fields
!! - Outsourced jb loop for OMP parallelization
!>
  REAL(wp), INTENT(IN)              :: &
    &  dtime,                          & !< time interval, fast physics
    &  p_qr(:,:),                      & !< rain water mass mixing ratio
    &  prr_gsp(:),                     & !< rain rate due to grid scale precipitation
    &  prr_con(:),                     & !< rain rate due to convective precipitation
    &  prr_con3d(:,:)                    !< 3D-rain rate due to convective precipitation
  INTEGER, INTENT(IN)               :: &
    &  istart, iend, kstart, nlev        !< Loop boundaries
  REAL(wp), INTENT(INOUT)           :: &
    &  tracer(:,:)                       !< Tracer field on which washout is calculated
  LOGICAL, OPTIONAL, INTENT(IN)     :: &
    &  lacc
!Local variables
  INTEGER  ::         &
    &  jk,jc            !< loop indices
  LOGICAL :: lzacc             ! OpenACC flag
  CALL set_acc_host_or_device(lzacc, lacc)

  !$ACC DATA PRESENT(prr_con3d, prr_con, prr_gsp, p_qr, tracer) IF(lzacc)

  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO jk= kstart, nlev
    DO jc= istart, iend
      IF ( p_qr(jc,jk)> 0.0_wp .OR. prr_con3d(jc,jk)> 0.0_wp) THEN
        IF ((prr_con(jc) + prr_gsp(jc)) > 0.0_wp) THEN
          ! To be consistent, multiplication and afterwards division with rho
          ! has to be done. But as this is not necessary it is not done to keep
          ! computational effort on a minimum
          tracer(jc,jk)= tracer(jc,jk)                       &
            &             - (3.*(prr_con(jc) + prr_gsp(jc))  &
            &             * tracer(jc,jk)                    &
            &             * dtime)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  !$ACC END PARALLEL

  !$ACC END DATA

END SUBROUTINE art_washout_volc
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_washout_volc
