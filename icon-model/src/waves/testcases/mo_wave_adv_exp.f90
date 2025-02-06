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

! Subroutine to initialize the wave test case

MODULE mo_wave_adv_exp

  USE mo_kind,                 ONLY: wp
  USE mo_model_domain,         ONLY: t_patch
  USE mo_wave_forcing_types,   ONLY: t_wave_forcing
  USE mo_wave_config,          ONLY: t_wave_config
  USE mo_math_constants,       ONLY: pi, rad2deg, deg2rad, dbl_eps
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, min_rlcell
  USE mo_loopindices,          ONLY: get_indices_c

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wind_adv_test

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_adv_exp'

  !--------------------------------------------------------------------

CONTAINS

  SUBROUTINE init_wind_adv_test(p_patch, wave_config, p_forcing)

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//'::init_wind_adv_test'

    TYPE(t_patch),        INTENT(IN)         :: p_patch
    TYPE(t_wave_config),  TARGET, INTENT(IN) :: wave_config
    TYPE(t_wave_forcing), INTENT(INOUT)      :: p_forcing

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: jc, jb
    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    REAL(wp):: sin_tmp, cos_tmp, zlat, zlon, d1, r

    REAL(wp):: lambda0, phi0  ! coordinates of center point
    REAL(wp), PARAMETER :: RR  = 1._wp/3._wp ! horizontal half width divided by 'a'

    wc => wave_config

    lambda0 = wc%peak_lon*deg2rad
    phi0    = wc%peak_lat*deg2rad

    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,zlon,zlat,sin_tmp,cos_tmp,r,d1)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,      &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      DO jc = i_startidx, i_endidx
        !test wind from NCAR_TESTCASE
        zlon = p_patch%cells%center(jc,jb)%lon
        zlat = p_patch%cells%center(jc,jb)%lat

        sin_tmp = SIN(zlat) * SIN(phi0)
        cos_tmp = COS(zlat) * COS(phi0)
        r  = ACOS (sin_tmp + cos_tmp*COS(zlon-lambda0))       ! great circle distance without 'a'
        d1 = MIN( 1._wp, (r/RR) )

        ! calculate U and V wind components and ensure nonzero values in order to
        ! avoid division by zero in the following ATAN2 function
        p_forcing%u10m(jc,jb) = MAX(0.5_wp*(1._wp + COS(pi*d1)) * wc%peak_u10, dbl_eps)
        p_forcing%v10m(jc,jb) = MAX(0.5_wp*(1._wp + COS(pi*d1)) * wc%peak_v10, dbl_eps)
        p_forcing%sp10m(jc,jb) = SQRT(p_forcing%u10m(jc,jb)**2 + p_forcing%v10m(jc,jb)**2)
        ! 45 degree towards NE
        p_forcing%dir10m(jc,jb) = ATAN2(p_forcing%v10m(jc,jb),p_forcing%u10m(jc,jb))*rad2deg
      END DO ! cell loop
    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE init_wind_adv_test

END MODULE mo_wave_adv_exp
