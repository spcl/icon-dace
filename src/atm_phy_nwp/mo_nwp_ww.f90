!  Module to calculate weather interpretation, WW.
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

MODULE mo_nwp_ww

#ifdef ONLYWW

  USE mo_wwonly, ONLY: wp, vct_a, &
                       rdv, O_m_rdv, rd_o_cpd, tmelt, alvdcp, b3, b1, b2w, b4w, &
                       max_dom, kstart_moist, sat_pres_water

#else

  USE mo_kind,               ONLY: wp

  USE mo_physical_constants, ONLY: rdv, O_m_rdv, tmelt,    &
                                   rd_o_cpd, alvdcp,          &
                                   b3    => tmelt
  USE mo_convect_tables,     ONLY: b1    => c1es  , & !! constants for computing the sat. vapour
                                   b2w   => c3les , & !! pressure over water (l) and ice (i)
                                   b4w   => c4les     !!               -- " --
  USE mo_impl_constants,     ONLY: max_dom

  USE mo_nonhydrostatic_config, ONLY: kstart_moist
  USE mo_vertical_coord_table,  ONLY: vct_a

  USE mo_exception,             ONLY: message, message_text
  USE mtime,                    ONLY: datetime, newDatetime, timeDelta

  USE mo_satad,                 ONLY: sat_pres_water
  USE mo_fortran_tools,         ONLY: set_acc_host_or_device

#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: configure_ww, ww_diagnostics
#ifdef ONLYWW
  PUBLIC :: rkgrenz1, rkgrenz2
  PUBLIC :: rf_fog, clc_fog, rain_l_m, rain_m_s, snow_l_m, snow_m_s, &
 &          rash_lm_s, rash_s_vs, snsh_l_ms, driz_l_m, driz_m_s,     &
 &          drif_l_ms, raif_l_ms, rgdiff_th1, rgdiff_th2
#else
  PUBLIC :: ww_datetime

  TYPE(datetime) :: ww_datetime(max_dom)
#endif

!  The following indices are used in the diagnostics of WW
  INTEGER :: ih_500hPa(max_dom),  &    ! index of layer next to 500 hPa for US standard atmosphere
             ih_700hPa(max_dom),  &    ! index of layer next to 700 hPa for US standard atmosphere
             ih_850hPa(max_dom),  &    ! index of layer next to 850 hPa for US standard atmosphere
             ih_950hPa(max_dom),  &    ! index of layer next to 950 hPa for US standard atmosphere
             ihb500hPa(max_dom),  &    ! index of next layer below 500 hPa
             ihb950hPa(max_dom)        ! index of next layer below 950 hPa
  !$ACC DECLARE CREATE(ih_500hPa, ih_700hPa, ih_850hPa, ih_950hPa, ihb500hpa, ihb950hpa)

! The following indices are used to calculate a fog index
  INTEGER :: ifog_temp(max_dom),  &    ! calculates minimum spread in this layer
             ifog_wind(max_dom)        ! calculate mean wind in this layer
  !$ACC DECLARE CREATE(ifog_temp, ifog_wind)

  REAL(wp) :: rkgrenz1,   &            ! limit for thunderstorm: convective precipitation in 1 hour
              rkgrenz2                 ! limit for strong thunderstorm: convective precip. in 1 hour
  REAL(wp) :: rf_fog                   ! threshold value for relative humidity indicating fog
  REAL(wp) :: clc_fog                  ! threshold value for cloudcover near ground indicating fog
  REAL(wp) :: rain_l_m, rain_m_s       ! tresholds for rain
  REAL(wp) :: snow_l_m, snow_m_s       ! tresholds for snow
  REAL(wp) :: rash_lm_s,rash_s_vs      ! tresholds for rain shower
  REAL(wp) :: snsh_l_ms                ! tresholds for snow shower
  REAL(wp) :: driz_l_m,driz_m_s        ! tresholds for drizzle
  REAL(wp) :: drif_l_ms                ! tresholds for freezing drizzle
  REAL(wp) :: raif_l_ms                ! tresholds for freezing rain
  REAL(wp) :: rgdiff_th1, rgdiff_th2   ! thresholds for weather with precipitation
                                       ! rgdiff_th1 (below) no precipitation
                                       ! rgdiff_th2 (below) only drizzle,freezing drizzle, snow griesel

CONTAINS

  SUBROUTINE configure_ww( ini_datetime, jg, nlev, nshift_total,  ymodel)
!
!   Set some constants for the calculation of WW
!
    CHARACTER(LEN=24), PARAMETER :: routine = 'atm_phy_nwp:configure_ww'

#ifdef ONLYWW
    INTEGER                                   :: ini_datetime  ! not used
#else
    TYPE(datetime),   INTENT(IN)              :: ini_datetime  ! init datetime (mtime)
#endif

    INTEGER,  INTENT(IN) :: jg           !< patch 
    INTEGER,  INTENT(IN) :: nlev         !< number of full vertical levels 
    INTEGER,  INTENT(IN) :: nshift_total 
    CHARACTER(LEN=4), INTENT(IN) :: ymodel   ! Use scalings factors for this NWP model

!   The following values are calculated as height difference between the height of the pressure level
!   and the height of 1000 hPa for a US standard atmosphere; h_1000hPa=110.88 m
    REAL(wp), PARAMETER :: h_500hPa = 5463.56_wp    ! height of 500 hPa level
    REAL(wp), PARAMETER :: h_700hPa = 2901.30_wp    ! height of 700 hPa level
    REAL(wp), PARAMETER :: h_850hPa = 1346.42_wp    ! height of 850 hPa level
    REAL(wp), PARAMETER :: h_950hPa =  429.46_wp    ! height of 950 hPa level
!   The followin heights are used to calculate a fox index
    REAL(wp), PARAMETER :: hfog_temp = 877.62_wp   ! h(900 hPa)-h(1000 hPa)
    REAL(wp), PARAMETER :: hfog_wind = 298.74_wp   ! h(965 hPa)-h(1000 hPa)

    INTEGER             :: jk, jk1

#ifndef ONLYWW
    ww_datetime(jg) = ini_datetime
#endif

!   Find next layer below height of 500 hPa, and next layer to height of 500 hPa
    ihb500hPa(jg) = kstart_moist(jg)
    DO jk = MAX( kstart_moist(jg),2) , nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_500hPa) THEN
        ihb500hPa(jg) = jk
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_500hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_500hPa) ) THEN
          ih_500hPa(jg) = jk
        ELSE
          ih_500hPa(jg) = jk-1
        END IF
        EXIT
      END IF
    ENDDO
    !$ACC UPDATE DEVICE(ih_500hPa(jg:jg), ihb500hPa(jg:jg)) ASYNC(1)

    ih_700hPa(jg) = ihb500hPa(jg)
    DO jk = ihb500hPa(jg), nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_700hPa) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_700hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_700hPa) ) THEN
          ih_700hPa(jg) = jk
        ELSE
          ih_700hPa(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE(ih_700hPa(jg:jg)) ASYNC(1)

    ih_850hPa(jg) = ih_700hPa(jg)+1
    DO jk = ih_700hPa(jg)+1, nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_850hPa) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_850hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_850hPa) ) THEN
          ih_850hPa(jg) = jk
        ELSE
          ih_850hPa(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE(ih_850hPa(jg:jg)) ASYNC(1)

    ifog_temp(jg) = ih_850hPa(jg)+1
    DO jk = ih_850hPa(jg)+1, nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < hfog_temp) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - hfog_temp) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - hfog_temp) ) THEN
          ifog_temp(jg) = jk
        ELSE
          ifog_temp(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE(ifog_temp(jg:jg)) ASYNC(1)

    ihb950hPa(jg) = ifog_temp(jg)
    DO jk = ifog_temp(jg), nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < h_950hPa) THEN
        ihb950hPa(jg) = jk
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - h_950hPa) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - h_950hPa) ) THEN
          ih_950hPa(jg) = jk
        ELSE
          ih_950hPa(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE(ih_950hPa(jg:jg), ihb950hPa(jg:jg)) ASYNC(1)

    ifog_wind(jg) = ih_950hPa(jg)
    DO jk = ih_950hPa(jg), nlev
      jk1 = jk + nshift_total
      IF  ( 0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) < hfog_wind) THEN
        IF ( ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1+1)) - hfog_wind) <  &
   &         ABS(0.5_wp*(vct_a(jk1)+vct_a(jk1-1)) - hfog_wind) ) THEN
          ifog_wind(jg) = jk
        ELSE
          ifog_wind(jg) = jk-1
        END IF
        EXIT
      ENDIF
    ENDDO
    !$ACC UPDATE DEVICE(ifog_wind(jg:jg)) ASYNC(1)

    SELECT CASE (ymodel)
      CASE( 'ICON')  ! these values are the same as for COSMO-EU (LME)
        rkgrenz1 = 0.75_wp
        rkgrenz2 = 2.0_wp
      CASE( 'ICEU')  ! these values are the same as for COSMO-EU (LME)
        rkgrenz1 = 0.75_wp
        rkgrenz2 = 2.0_wp
      CASE( 'GME ')
        rkgrenz1 = 0.2_wp
        rkgrenz2 = 1.0_wp
      CASE( 'GME6')         ! old GME with mesh width <= 60 km (NI<=128)
        rkgrenz1 = 0.25_wp
        rkgrenz2 = 1.0_wp
      CASE( 'LME ')
        rkgrenz1 = 0.75_wp
        rkgrenz2 = 2.0_wp
      CASE( 'LMK ')
        rkgrenz1 = 0.60_wp
        rkgrenz2 = 2.0_wp
      CASE( 'TEST')
        PRINT *, 'Test ww_diagnostics with rkgrenz1=', rkgrenz1, '  rkgrenz2=',rkgrenz2
      CASE DEFAULT
        PRINT *, 'Unknown model ', ymodel,' in configure_ww. Use ICON values!'
        rkgrenz1 = 0.75_wp
        rkgrenz2 = 2.0_wp
    END SELECT

    rf_fog    =  95.   ! %
    clc_fog   =  99.   ! %
    rain_l_m  =  2.5   ! mm/h
    rain_m_s  = 10.0   ! mm/h
    snow_l_m  =  1.00  ! mm/h
    snow_m_s  =  5.0   ! mm/h
    rash_lm_s =  2.5   ! mm/h
    rash_s_vs = 20.0   ! mm/h
    snsh_l_ms =  1.00  ! mm/h
    driz_l_m  =  0.1   ! mm/h
    driz_m_s  =  0.5   ! mm/h
    drif_l_ms =  0.1   ! mm/h
    raif_l_ms =  2.5   ! mm/h
    rgdiff_th1=  0.015 ! mm/h
    rgdiff_th2=  0.050 ! mm/h

#ifdef ONLYWW
    WRITE( *,'(a,i4,a,6i4)') 'Domain', jg,                                    &
  &    '; ih_500hPa, ihb500hPa, ih_700hPa, ih_850hPa, ih_950hPa, ihb950hPa ', &
  &       ih_500hPa(jg), ihb500hPa(jg), ih_700hPa(jg), ih_850hPa(jg),         &
  &       ih_950hPa(jg), ihb950hPa(jg)
    WRITE( *,'(2(a,i4),i4)') 'Domain', jg,                                    &
  &    '; ifog_temp, ifog_wind ', ifog_temp(jg), ifog_wind(jg)
#else
    WRITE( message_text,'(a,i4,a,6i4)') 'Domain', jg,                         &
  &    '; ih_500hPa, ihb500hPa, ih_700hPa, ih_850hPa, ih_950hPa, ihb950hPa ', &
  &       ih_500hPa(jg), ihb500hPa(jg), ih_700hPa(jg), ih_850hPa(jg),         &
  &       ih_950hPa(jg), ihb950hPa(jg)
    CALL message( TRIM(routine), message_text)
    WRITE( message_text,'(2(a,i4),i4)') 'Domain', jg,                         &
  &    '; ifog_temp, ifog_wind ', ifog_temp(jg), ifog_wind(jg)
    CALL message( TRIM(routine), message_text)
#endif

  END SUBROUTINE configure_ww


  SUBROUTINE ww_diagnostics( ie, ke, ke1, i_startidx, i_endidx, jg,           &
                             t   , qv   , qc , u   , v   , clc  , pf   , ph,  &
                             t_2m, td_2m, t_g, clct, clcm, u_10m, v_10m,      &
                             rain_gsp0, rain_gsp, rain_con0, rain_con,        &
                             snow_gsp0, snow_gsp, snow_con0, snow_con,        &
                             bas_con, top_con, time_diff, iww, lacc )

!
! Description:
! The subroutime evaluates for each grid point of the model a "ww-number";
! that means, a number corresponding to those numbers that
! would be used for weather observation (e.g. thunderstorm=95, 
! heavy rain=65 etc.)
! This subroutine is based on auswert_lm.f of the ww-library of DWD.


  INTEGER, INTENT(IN)  :: ie, ke, ke1, &        ! horizontal and vertical dimensions
                          i_startidx, i_endidx  ! start and end indices of loops in horizontal
  INTEGER, INTENT(IN)  :: jg                    !  patch 

  REAL(wp), INTENT(IN) :: t  (ie,ke),  &  ! temperature (K)
                          qv (ie,ke),  &  ! specific humidity (on model levels) (kg/kg)
                          qc (ie,ke),  &  ! specific cloud liquid water (on model levels) (kg/kg)
                          u  (ie,ke),  &  ! zonal comp. of wind vel. (m/s)
                          v  (ie,ke),  &  ! meridional comp. of wind vel. (m/s)
                          clc(ie,ke),  &  ! cloud cover (on model levels) (per cent)
                          pf (ie,ke),  &  ! pressure on full levels (Pa)
                          ph (ie,ke1)     ! pressure on half levels (Pa)

  REAL(wp), INTENT(IN) :: t_2m (ie),  &  ! temperature at 2m above the ground (K)
                          td_2m(ie),  &  ! dewpoint temperature at 2m above the ground (K)
                          t_g  (ie),  &  ! surface temperature (weighted from t_s and t_snow) (K)
                          clct (ie),  &  ! total cloud cover
                          clcm (ie),  &  ! medium cloud cover
                          u_10m(ie),  &  ! zonal wind component at 10m above the ground (m/s)
                          v_10m(ie),  &  ! meridional wind component at 10m above the ground (m/s)
                          rain_gsp0(ie), & ! grid scale rain, valid at previous call of ww_diagnostics
                          rain_gsp (ie), & ! grid scale rain, valid now
                          rain_con0(ie), & ! convective rain, valid at previous call of ww_diagnostics
                          rain_con (ie), & ! convective rain, valid now
                          snow_gsp0(ie), & ! grid scale snow, valid at previous call of ww_diagnostics
                          snow_gsp (ie), & ! grid scale snow, valid now
                          snow_con0(ie), & ! convective snow, valid at previous call of ww_diagnostics
                          snow_con (ie)    ! convective snow, valid now

  INTEGER,  INTENT(IN) :: bas_con(ie),   & ! base index of main convective cloud
                          top_con(ie)      ! top index of main convective cloud

  LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

#if ONLYWW
  REAL(wp)                        :: time_diff ! time since previous call in hours
#else
  TYPE(timeDelta),        POINTER :: time_diff ! time since previous call
  CHARACTER(LEN=26), PARAMETER    :: routine = 'atm_phy_nwp:ww_diagnostics'
#endif

  INTEGER, INTENT(OUT) :: iww(ie)          ! significant weather code

    REAL(wp) :: rgdiff         !  Difference total precip. in mm/hour
    REAL(wp) :: rkdiff         !  Difference convective precip. in mm/hour
    REAL(wp) :: rldiff         !  Difference grid scale precip. in mm/hour
    REAL(wp) :: rrdiff         !  Difference rain in mm/hour
    REAL(wp) :: rsdiff         !  Difference snow in mm/hour

    INTEGER :: i, k
    INTEGER :: iwolk, iwolkc

    REAL(wp) :: rko            ! KO-Index
    REAL(wp) :: vbetr          ! wind speed in 10 m a. g.
    REAL(wp) :: tblmax         ! max. value of t in lowest 100 hPa
    REAL(wp) :: tdbl, tdblmax  ! value/max. value of td in lowest 100 hPa
    REAL(wp) :: rfbl, rfblmax  ! value/max. value of rf in lowest 100 hPa
    REAL(wp) :: vbl            ! mean value of wind speed in lowest 35 hPa
    REAL(wp) :: neb_i          ! Nebulae-Index
    REAL(wp) :: zfrac, qvmin   ! variables used in determination of dew point temperature
    REAL(wp) :: dp             ! pressure difference
    REAL(wp) :: e              ! vapour pressure
    REAL(wp) :: e_s            ! saturation vapour pressure
    INTEGER  :: igfb, irrb, isprb
    LOGICAL  :: test_fog
    LOGICAL  :: lconvb, lsnb

  REAL(wp) :: t_t  (ke,ie),  &  ! temperature (K)
              ph_t (ke1,ie),  &  ! specific humidity (on model levels) (kg/kg)
              qc_t (ke,ie)  ! specific cloud liquid water (on model levels) (kg/kg)

    REAL(wp), PARAMETER :: rkogrenz = 1._wp          ! limit for rko
    REAL(wp), PARAMETER :: ms_kn = 1._wp/0.51444_wp  ! factor to convert numbers in m/s to knots

    LOGICAL :: lzacc ! non-optional version of lacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(t_t, ph_t, qc_t) PRESENT(t, qc, ph) IF(lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO i = 1, ie
      DO k = 1, ke
        t_t(k,i) = t(i,k)
        qc_t(k,i) = qc(i,k)
      end do
    end do
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO i = 1, ie
      DO k = 1, ke1
        ph_t(k,i) = ph(i,k)
      end do
    end do
    !$ACC END PARALLEL


    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(test_fog, rgdiff, lsnb, irrb, igfb, iwolk) &
    !$ACC   PRIVATE(iwolkc, isprb, rgdiff, rkdiff, rldiff, rrdiff, rsdiff, dp, lconvb) &
    !$ACC   PRIVATE(vbetr, tblmax, tdblmax, rfblmax, vbl, qvmin, zfrac, tdbl, rfbl) &
    !$ACC   PRIVATE(neb_i, e, e_s)
    DO i = i_startidx, i_endidx
      iww(i) = -9
      test_fog = .FALSE.
      rgdiff = -1._wp
      lsnb = .true.
      irrb =  9
      igfb =  9
      iwolk = 9
      isprb = 9
      rgdiff =  rain_gsp (i) + snow_gsp (i) + rain_con (i) + snow_con (i)     &
   &          - rain_gsp0(i) - snow_gsp0(i) - rain_con0(i) - snow_con0(i)

      !.... Total precipitation < 0.015 mm/h => weather without precipitation
WW_PRECIP: IF (rgdiff < rgdiff_th1) THEN
        test_fog = .TRUE.

      ELSE WW_PRECIP

        rkdiff = rain_con(i) + snow_con(i) - rain_con0(i) - snow_con0(i)
        rldiff = rain_gsp(i) + snow_gsp(i) - rain_gsp0(i) - snow_gsp0(i)
        rrdiff = rain_gsp(i) + rain_con(i) - rain_gsp0(i) - rain_con0(i)
        rsdiff = snow_gsp(i) + snow_con(i) - snow_gsp0(i) - snow_con0(i)

!.... Thunderstorm

        IF (top_con(i) > 0 .AND. top_con(i) <= ke1) THEN
          dp = ph(i,bas_con(i)) - pf(i,top_con(i))
        ELSE
          dp   = 0.0_wp
        ENDIF

!     Determination of KO-Index (rko)

        CALL ko_index( pf(i,ih_950hPa(jg)), pf(i,ih_850hPa(jg)), pf(i,ih_700hPa(jg)), pf(i,ih_500hPa(jg)),  &
     &                   t(i,ih_950hPa(jg)),  t(i,ih_850hPa(jg)),  t(i,ih_700hPa(jg)),  t(i,ih_500hPa(jg)),  &
     &                  qv(i,ih_950hPa(jg)), qv(i,ih_850hPa(jg)), qv(i,ih_700hPa(jg)), qv(i,ih_500hPa(jg)),  &
     &                  rko)

        IF ( (rkdiff > rkgrenz1 .AND. (top_con(i) > 0 .AND. top_con(i) <= ke1))     &
     &        .AND.( dp > 400.e2_wp .AND. t(i, MAX(1,top_con(i))) < tmelt-25._wp)   &
     &        .AND.(rko < rkogrenz) ) THEN

          IF (rko < -6._wp .AND. t(i, top_con(i)) < tmelt-45._wp .AND.  rkdiff > rkgrenz2) THEN
            iww(i)   = 96
          ELSE
            iww(i)   = 95
          ENDIF

        ELSE

          iwolk = 0
          iwolkc = 0
          !$ACC LOOP SEQ
          DO k = ihb950hPa(jg)+1, ke
            IF ( qc(i,k) > 1.e-8_wp) THEN
              iwolkc = 1
              EXIT
            ENDIF
          ENDDO
          IF (iwolkc == 1.AND. clcm(i) < 0.95_wp) iwolk = 1

          CALL gefr( jg, ke, ke1, ph_t(:,i), qc_t(:,i), t_t(:,i), t_2m(i), t_g(i), &
                    igfb, irrb, isprb )

!.... Convective precipitation

          lconvb  = ( rkdiff > rldiff)
          lsnb    =  ( ( (rsdiff > rrdiff .AND. t_2m(i) <= tmelt+5._wp) .OR.  &
       &                 t_2m(i) < tmelt-10._wp) )

          IF (lconvb) THEN
            IF (lsnb .OR. irrb == 0) THEN
!.... Snow shower (threshold value 0.05 due to convective drizzle problem in ICON/ICON-EU)
              IF ( rgdiff >= MAX(0.05_wp,rgdiff_th2)) THEN
                IF (rgdiff >  snsh_l_ms) THEN
                  iww(i) = 86 
                ELSE
                  iww(i) = 85 
                ENDIF
              ELSE
                test_fog = .true.
              ENDIF
            ELSE IF (( .NOT. lsnb) .AND. irrb == 1) THEN
              IF(rgdiff >= rgdiff_th2) THEN
                IF (igfb == 1) THEN
!....             Freezing rain (without drizzle)
                  IF (rgdiff > raif_l_ms) THEN
                    iww(i) = 67
                  ELSE
                    iww(i) = 66
                  ENDIF
                ELSE

!....             Rain shower (threshold value 0.2 due to convective drizzle problem in ICON/ICON-EU)
                  IF(rgdiff > 0.2_wp) THEN
                    IF (rgdiff <= rash_lm_s) THEN
                      iww(i) = 80
                    ELSE IF (rgdiff <= rash_s_vs) THEN
                      iww(i) = 81
                    ELSE
                      iww(i) = 82
                    END IF
                  ELSE
                    test_fog = .true.
                  ENDIF
                ENDIF
              ELSE
                test_fog = .true.
              ENDIF
            ENDIF


          ELSE ! IF(.not.lconvb)

!.... Large-scale precipitaion
            IF (iwolk == 0) THEN
!.... Snow/Rain/Freezing rain (no drizzle/snow griesel)         
              IF ((lsnb) .OR. (irrb == 0)) THEN
!.... Snow
                IF (rgdiff >= rgdiff_th2) THEN
                  IF (rgdiff <= snow_l_m) THEN
                    iww(i) = 71
                  ELSE IF (rgdiff <= snow_m_s) THEN
                    iww(i) = 73
                  ELSE
                    iww(i) = 75
                  END IF
                ELSE
                  test_fog = .true.
                ENDIF
              ELSE IF ((.not.lsnb) .AND. irrb == 1) THEN
                IF(rgdiff >= rgdiff_th2) THEN
                  IF (igfb > 0) THEN

!.... Freezing rain
                    IF (rgdiff <= raif_l_ms) THEN
                      iww(i) = 66
                    ELSE
                      iww(i) = 67
                    ENDIF
                  ELSE
!.... Rain
                    IF (rgdiff <= rain_l_m) THEN
                      iww(i) = 61
                    ELSE IF (rgdiff <= rain_m_s) THEN
                      iww(i) = 63
                    ELSE
                      iww(i) = 65
                    ENDIF
                  ENDIF
                ELSE
                  test_fog = .true.
                ENDIF
              ENDIF

            ELSE ! IF (iwolk == 1) THEN
!.... Drizzle/Freezing drizzle/Snow griesel 
              IF (( .NOT. lsnb .AND. igfb < 1)) THEN
!.... Drizzle (must be evaluated before Freezing drizzle/Snow griesel)
                IF (rgdiff <= driz_l_m) THEN
                  iww(i) = 51
                ELSE IF (rgdiff <= driz_m_s) THEN
                  iww(i) = 53
                ELSE
                  iww(i) = 55
                ENDIF
!.... Freezing drizzle/Snow griesel
              ELSE IF (.not.lsnb .AND. igfb > 0 .AND. isprb == 1) THEN
!.... Freezing drizzle
                iww(i) = 57
                IF (rgdiff <= drif_l_ms) iww(i) = 56
              ELSE
!.... Snow griesel
!     cases     IF (.not.lsnb .AND. igfb > 0 .AND. isprb == 0) THEN
!               or
!               IF ( lsnb ) THEN
                iww(i) = 77
              ENDIF !ok
            ENDIF ! iwolk=0/1
          ENDIF ! lconvb
        ENDIF ! thunderstorm
      ENDIF WW_PRECIP ! weather with preciptitation

!.... Weather without precipitation

!.... Fog
      IF (test_fog ) THEN

        vbetr   = SQRT (u_10m(i)**2 + v_10m(i)**2)
        tblmax  = t_2m(i)
        tdblmax = td_2m(i)
        rfblmax = 0._wp
        vbl = 0._wp
        DO k = ifog_wind(jg), ke
          vbl = vbl + SQRT(u(i,k)**2 + v(i,k)**2)*(ph(i,k)-ph(i,k+1))
        ENDDO
        vbl = vbl / ( ph(i,ifog_wind(jg))-ph(i,ke1) )

        DO k = ifog_temp(jg), ke
          tblmax  = MAX( t(i,k), tblmax)
          qvmin   = MAX( qv(i,k), 1.e-12_wp)
          zfrac   = LOG( pf(i,k)*qvmin/(b1*(rdv+O_m_rdv*qvmin)))
          tdbl    = MIN( (b2w*b3-b4w*zfrac)/(b2w-zfrac), t(i,k) )
          tdblmax = MAX( tdbl, tdblmax)
!.... calculate rf for k=ke-1 and k=ke
          IF (k >= ke-1) THEN
            e_s = sat_pres_water(t(i,k))
            e   = sat_pres_water(tdbl)
! compute vapor pressure from formula for specific humidity:
!DR           e   = pf(i,k)*qv(i,k) / (rdv + o_m_rdv*qv(i,k))
            rfbl = 100._wp * e/e_s
            rfblmax = MAX(rfbl,rfblmax)
          ENDIF
        ENDDO
        neb_i = ( t_2m(i)-tblmax) + ( t_2m(i)-tdblmax) + ms_kn*vbl

        IF ( clc(i,ke) > clc_fog .OR.  clc(i,ke-1) > clc_fog       &
     &  .OR. neb_i     < 5._wp   .AND. rfblmax     > rf_fog ) THEN
!           Fog
          IF ( (t_2m(i) < tmelt-1.0_wp .OR. t_g(i) < tmelt-1.0_wp)   &
        &     .AND. vbetr > 1.5_wp) THEN
             iww(i) = 48
          ELSE
             iww(i) = 45
          END IF


        ELSE

!     Clouds
           iww(i) = clct2ww( clct(i) )

        ENDIF

      ENDIF

!.... Control Printout
#ifndef _OPENACC
      IF ( iww(i) == -9 ) THEN
#ifdef ONLYWW
           WRITE(*,'(a,i7,a,i3,a,f6.3,2(a,f5.1),2(a,L1),4(a,i2),a,L1)') &
     &     'ww(',i,'):',iww(i), ' rg:', rgdiff,               &
     &     ' clc(ke):', clc(i,ke),' clc(ke-1):', clc(i,ke-1), &
     &     ' lconvb: ',lconvb,' lsnb: ',lsnb,' irrb: ',irrb,  &
     &     ' igfb: ',igfb,' iwolk: ',iwolk,' isprb: ',isprb,  &
     &     ' test_fog: ',test_fog
#else
           WRITE(message_text,'(a,i7,a,i3,a,f6.3,2(a,f5.1),2(a,L1),4(a,i2),a,L1)') &
     &     'ww(',i,'):',iww(i), ' rg:', rgdiff,               &
     &     ' clc(ke):', clc(i,ke),' clc(ke-1):', clc(i,ke-1), &
     &     ' lconvb: ',lconvb,' lsnb: ',lsnb,' irrb: ',irrb,  &
     &     ' igfb: ',igfb,' iwolk: ',iwolk,' isprb: ',isprb,  &
     &     ' test_fog: ',test_fog
    CALL message( TRIM(routine), message_text)
#endif
      ENDIF
#endif
   
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE ww_diagnostics


  SUBROUTINE ko_index ( p950,  p850,  p700,  p500,      &
     &                  t950,  t850,  t700,  t500,      &
     &                 qv950, qv850, qv700, qv500,      &
     &                 rko)
!
! Description:
! Calculation of thunderstorm indices KO as for COSMO-EU as in subroutine index.f in ww-library.
!
! Method:
! Calculation of the equivalent potential temperature: then,  
! using the appropriate formulas, calculation of KO index.
!

!=======================================================================
!
      REAL(wp), INTENT(IN) :: p950, p850, p700, p500       ! single values of pressure field p in hPa
      REAL(wp), INTENT(IN) :: t950, t850, t700, t500       ! single values of temperature field t in degree C
      REAL(wp), INTENT(IN) :: qv950, qv850, qv700, qv500   ! single values of specific humidity field qv

      REAL(wp), INTENT(OUT)::  rko                         ! KO-Index
!-----------------------------------------------------------------------
!
      REAL(wp), PARAMETER :: qvmin = 1.e-12_wp            ! minimum value for specific humidity
      REAL(wp)            :: ta950, ta850, ta700, ta500   ! equivalent potential temperature

      !$ACC ROUTINE SEQ
!
  
!     950 hPa
      ta950 = ( t950 + alvdcp*MAX( qvmin, qv950) )*(1.e5/p950)**rd_o_cpd

!     850 hPa
      ta850 = ( t850 + alvdcp*MAX( qvmin, qv850) )*(1.e5/p850)**rd_o_cpd

!     700 hPa
      ta700 = ( t700 + alvdcp*MAX( qvmin, qv700) )*(1.e5/p700)**rd_o_cpd

!     500 hPa
      ta500 = ( t500 + alvdcp*MAX( qvmin, qv500) )*(1.e5/p500)**rd_o_cpd

      rko = 0.5_wp *( ta700 + ta500 - ta950 - ta850 )

    END SUBROUTINE ko_index


    SUBROUTINE gefr ( jg, ke, ke1, p_gefr, qc_gefr, t_gefr, t2m, tg, &
                      igfb, irrb, isprb )
!
! Description:
! Determination of igfb, irrb and isprb (indicators for freezing)
!
! Method:
! Main part:
!    Set igfb, irrb and isprb to 0. If the main freezing condition is fulfilled,
!    set igfb to 1. From level ihb500hPa(jg) to ke (bottom level of the model) calculate
!    the products t*p, if t>0; build the sum of those products. Dependent on
!    this sum resp. the value of t_2m  and t_g set the value for irrb resp.
!    igfb.
!
! Additional conditions for drizzle:
!   Determine the index of the bottom layer of low clouds (kwug); 
!   dependent on the value of t at this level continue (or return to the
!   calling program):
!   If kwug equals ke set isprb to 1 and return to the calling program);
!   otherwise build the sum of the products t*p from kwug+1 to ke; dependent on
!   the value of this sum resp. the value of the temperature at level ke set
!   isprb and return to the calling program.
!
!=======================================================================
!
      INTEGER, INTENT(IN)  :: jg                ! patch index
      INTEGER, INTENT(IN)  :: ke, ke1           ! number of levels
      REAL(wp), INTENT(IN) :: p_gefr(ke1)       ! pressure at half levels
      REAL(wp), INTENT(IN) :: qc_gefr(ke)       ! profile of cloud water qc
      REAL(wp), INTENT(IN) :: t_gefr(ke)        ! profile of temperature t
      REAL(wp), INTENT(IN) :: t2m, tg           ! 2m temp. t_2m, surface temp. t_g

      INTEGER, INTENT(OUT) :: igfb,         &   ! indicator for freezing
                              irrb,         &   ! indicator for rain
                              isprb             ! indicator for drizzle
!
! Local scalars:

      INTEGER kk                         ! loop variable
      INTEGER :: kwug                    ! index of bottom level of low clouds

      REAL(wp) :: sumtpr                 ! sum of products t*p for t>0 from level ihb500hPa(jg) to ke
      REAL(wp) :: sumtpsp                ! sum of products t*p from level kwug+1 to ke
      REAL(wp) :: dp                     ! help variables to calculate
                                         ! sumtpr resp. sumtpsp
      REAL(wp), PARAMETER :: xtpsp = -20.e2_wp          ! limit for sumtpsp
      REAL(wp), PARAMETER :: xtsp  = tmelt - 2.0_wp     ! limit for t at level kwug
      REAL(wp), PARAMETER :: xtpr  = 100.e2_wp          ! limit for sumtpr

      !$ACC ROUTINE SEQ
!
!- End of header
!

!   **********************************************************************
!   * Limits for t und t*p (variable name ending with "r" for rain,      *
!   * ending with "sp" for drizzle)                                      *
!   **********************************************************************

      igfb = 0
      irrb = 0
      isprb = 0

!.... Freezing condition fulfilled

      IF ( t2m < tmelt-0.1_wp .OR. tg < tmelt-0.1_wp) igfb = 1

!.... Additional condition for rain

      sumtpr = 0._wp

      DO kk = ihb500hPa(jg), ke
         IF ( t_gefr(kk) > tmelt) THEN
            dp = p_gefr(kk+1)  - p_gefr(kk)
            sumtpr = sumtpr + (t_gefr(kk)-tmelt) * dp
         ENDIF
      ENDDO

      IF (sumtpr > xtpr .OR. (t2m > tmelt+0.5_wp .AND. tg >= tmelt-0.1_wp)) THEN
         irrb = 1
      ELSE
         IF (igfb == 1) igfb = 2
      ENDIF

      kwug = ihb950hPa(jg)

      DO kk = ihb950hPa(jg)+1, ke
         IF ( qc_gefr(kk) > 1.e-8) kwug = kk
      ENDDO   

!.... First additional condition

      IF ( t_gefr(kwug) <= xtsp) RETURN

!.... Return, if kwug = ke

      IF (kwug == ke) THEN
         isprb = 1
         RETURN
      ENDIF

!.... Second additional

      sumtpsp = 0._wp
      DO  kk = kwug+1, ke
         dp = p_gefr(kk+1)  - p_gefr(kk)
         sumtpsp = sumtpsp + (t_gefr(kk)-tmelt) * dp
      ENDDO

      IF (sumtpsp > xtpsp .OR. t_gefr(ke) > tmelt) THEN
         isprb = 1
      ENDIF

    END SUBROUTINE gefr

    INTEGER FUNCTION clct2ww( clct)
!
!   WW key (significant weather) for different cloud covers
!
      REAL(wp), INTENT(IN) :: clct
      !$ACC ROUTINE SEQ
      IF (      clct > 0.8125_wp) THEN
        clct2ww = 3
      ELSE IF ( clct > 0.4375_wp) THEN
        clct2ww = 2
      ELSE IF ( clct > 0.0625_wp) THEN
        clct2ww = 1
      ELSE
        clct2ww = 0
      END IF
    END  FUNCTION clct2ww
    

END MODULE mo_nwp_ww
