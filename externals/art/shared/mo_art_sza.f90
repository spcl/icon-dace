!
! mo_art_sza
! This module provides the sza calculation calculation
!
! original from atm_phy_aes/mo_aes_phy_main.f90
!
! elliptical seasonal orbit,
! with diurnal cycle
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

MODULE mo_art_sza

  ! ICON
  USE mo_kind,                 ONLY: wp
  USE mo_math_constants,       ONLY: pi, rad2deg

  USE mtime,                   ONLY: datetime, no_of_sec_in_a_day,   &
                                 &   getDayOfYearFromDateTime
  USE mo_fortran_tools,        ONLY: set_acc_host_or_device

  ! ART
  USE mo_art_data,               ONLY: p_art_data
  USE mo_art_atmo_data,          ONLY: t_art_atmo
  USE mo_art_wrapper_routines,   ONLY: art_get_indices_c


  IMPLICIT NONE


  PRIVATE

  PUBLIC   ::                             &
      &   art_calc_sza

  CONTAINS


  SUBROUTINE art_calc_sza(current_date, jg, lacc)


    TYPE(datetime), POINTER, INTENT(IN)  ::  current_date !< actual date
    INTEGER, INTENT(in)                  ::  jg           !< patch id
    LOGICAL, INTENT(in), OPTIONAL        ::  lacc ! If true, use openacc


    INTEGER               ::  year,hour,minute,second,dayinyear


    TYPE(t_art_atmo),POINTER    :: &
      &  art_atmo                     !< Pointer to ART diagnostic fields

    REAL(wp) :: zleapfrac
    REAL(wp) :: zyearfrac
    REAL(wp) :: zdeclination_sun
    REAL(wp) :: ztime_dateline

    INTEGER ::                     &
       &    jc, jb,                &                   !< loop indizes
       &    i_startidx, i_endidx

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    ! ----------------------------------
    ! --- Get the loop indices
    ! ----------------------------------

    art_atmo => p_art_data(jg)%atmo


    year   = current_date%date%year

    hour   = current_date%time%hour
    minute = current_date%time%minute
    second = current_date%time%second
    dayinyear = getDayOfYearFromDateTime(current_date)


!
    zleapfrac = 0.681_wp + 0.2422_wp * REAL(year - 1949,wp) - &
                REAL((year - 1949) / 4,wp)
    zyearfrac = 2._wp * pi * (REAL(dayinyear,wp) - 1.0_wp + zleapfrac) / 365.2422_wp
    zdeclination_sun = 0.006918_wp - 0.399912_wp * COS(zyearfrac) + &
                       0.070257_wp * SIN(zyearfrac) -               &
                       0.006758_wp * COS(2._wp * zyearfrac) +       &
                       0.000907_wp * SIN(2._wp * zyearfrac) -       &
                       0.002697_wp * COS(3._wp * zyearfrac) +       &
                       0.001480_wp * SIN(3._wp * zyearfrac)
    ztime_dateline = ((REAL(hour,wp) * 3600._wp + &
                      REAL(minute,wp) * 60._wp +  &
                      REAL(second,wp)) /          &
                      REAL(no_of_sec_in_a_day,wp)) - 0.5_wp
    ztime_dateline = ztime_dateline * 2._wp * pi + 0.000075_wp +              &
       0.001868_wp * COS(zyearfrac) - 0.032077_wp * SIN(zyearfrac) -          &
       0.014615_wp * COS(2._wp * zyearfrac) - 0.040849_wp * SIN(2._wp * zyearfrac)

    !$ACC DATA PRESENT(art_atmo) IF(lzacc)

    DO jb = art_atmo%i_startblk, art_atmo%i_endblk
      CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO jc = i_startidx, i_endidx
        art_atmo%sza(jc,jb)     = SIN(zdeclination_sun) * SIN(art_atmo%lat(jc,jb)) + &
                                  COS(zdeclination_sun) * COS(art_atmo%lat(jc,jb)) * &
                                  COS(ztime_dateline + art_atmo%lon(jc,jb))

!         art_atmo%sza_deg(jc,jb)  = ACOS(art_atmo%sza(jc,jb)) * 180._wp / rad2deg
        art_atmo%sza_deg(jc,jb)  = ACOS(art_atmo%sza(jc,jb)) * rad2deg
      ENDDO
      !$ACC END PARALLEL
    ENDDO

    !$ACC END DATA


END SUBROUTINE art_calc_sza

END MODULE mo_art_sza
