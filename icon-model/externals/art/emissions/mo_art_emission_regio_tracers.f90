!
! mo_art_emission_regio_tracers
! This module provides the emission subroutines for regional tracers
! It is based on Vogel et al., ACP, 2015
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

MODULE mo_art_emission_regio_tracers
  ! ICON
  USE mo_kind,                   ONLY: wp
  USE mo_math_constants,         ONLY: rad2deg

  ! ART
  USE mo_art_data,               ONLY: p_art_data
  USE mo_art_atmo_data,          ONLY: t_art_atmo
  USE mo_art_chem_data,          ONLY: t_art_chem_indices
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC   ::                        &
      &   art_emiss_TR_art,          &
      &   art_emiss_TR_bgs,          &
      &   art_emiss_TR_bgn,          &
      &   art_emiss_TR_ech,          &
      &   art_emiss_TR_eur,          &
      &   art_emiss_TR_naf,          &
      &   art_emiss_TR_saf,          &
      &   art_emiss_TR_sea,          &
      &   art_emiss_TR_mdg,          &
      &   art_emiss_TR_sam,          &
      &   art_emiss_TR_nam,          &
      &   art_emiss_TR_aus,          &
      &   art_emiss_TR_med,          &
      &   art_emiss_TR_sib,          &
      &   art_emiss_TR_sin,          &
      &   art_emiss_TR_nin,          &
      &   art_emiss_TR_tpo,          &
      &   art_emiss_TR_tio,          &
      &   art_emiss_TR_tao   

  REAL(wp), PARAMETER, DIMENSION(4) :: & !(lat_south, lat_north, lon_west, lon_east)
    &  area_ech = (/   20._wp,   40._wp,   90._wp,  125._wp /) / rad2deg,  &
    &  area_eur = (/   45._wp,   75._wp,  -20._wp,   55._wp /) / rad2deg,  &
    &  area_naf = (/    0._wp,   35._wp,  -20._wp,   55._wp /) / rad2deg,  &
    &  area_saf = (/  -36._wp,    0._wp,    7._wp,   42._wp /) / rad2deg,  &
    &  area_sea = (/  -12._wp,   20._wp,   90._wp,  155._wp /) / rad2deg,  &
    &  area_mdg = (/  -27._wp,  -12._wp,   42._wp,   52._wp /) / rad2deg,  &
    &  area_sam = (/  -55._wp,   15._wp,  -80._wp,  -35._wp /) / rad2deg,  &
    &  area_nam = (/   15._wp,   75._wp, -160._wp,  -50._wp /) / rad2deg,  &
    &  area_aus = (/  -40._wp,  -12._wp,  110._wp,  155._wp /) / rad2deg,  &
    &  area_med = (/   35._wp,   45._wp,  -20._wp,   55._wp /) / rad2deg,  &
    &  area_sib = (/   40._wp,   75._wp,   55._wp,  180._wp /) / rad2deg,  &
    &  area_sin = (/    0._wp,   20._wp,   55._wp,   90._wp /) / rad2deg,  &
    &  area_tpo1= (/  -20._wp,   20._wp,  155._wp,  180._wp /) / rad2deg,  &
    &  area_tpo2= (/  -20._wp,   20._wp, -180._wp,  -60._wp /) / rad2deg,  &
    &  area_tao = (/  -20._wp,   20._wp,  -60._wp,   20._wp /) / rad2deg,  &
    &  area_tio = (/  -20._wp,   20._wp,   20._wp,  110._wp /) / rad2deg,  &
    &  area_nin = (/   20._wp,   40._wp,   55._wp,   90._wp /) / rad2deg

  INTEGER, PARAMETER :: &
    &  num_levs_emiss = 3


  CONTAINS

FUNCTION is_in_area(lat, lon, area)
!<
! FUNCTION is_in_area
! This routine returns if lat lon are within the four corners given in "area"
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<
  REAL(wp), INTENT(IN) :: &
    &  lat, lon, area(4)  !< everything in radian
  LOGICAL :: &
    &  is_in_area

  is_in_area = (lat >=  area(1) .AND. lat <=  area(2) .AND. lon >=  area(3) .AND. lon <= area(4))
END FUNCTION is_in_area

SUBROUTINE art_emiss_regio(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer, area)
!<
! SUBROUTINE art_emiss_TR_art
! This subroutine sets the artificial tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer(:,:,:)
  REAL(wp), INTENT(IN) :: area(4)

  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

  art_atmo    => p_art_data(jg)%atmo

  DO i = jcs, jce
    IF (is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area)) THEN

      DO jk = nlev-num_levs_emiss,nlev
        p_tracer(i,jk,jb) = 1.0_wp
      END DO
    ENDIF
  END DO
END SUBROUTINE art_emiss_regio



SUBROUTINE art_emiss_TR_art(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_art)
!<
! SUBROUTINE art_emiss_TR_art
! This subroutine sets the artificial tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_art(:,:,:)

  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

  art_atmo    => p_art_data(jg)%atmo

  DO i = jcs, jce
    IF (ABS(art_atmo%lat(i,jb))*rad2deg <= 15.) THEN
      DO jk=1,nlev
        IF (art_atmo%pres(i,jk,jb) <= 200.*100.       &
          & .AND. art_atmo%pres(i,jk,jb) >= 150.*100. ) THEN

          p_tracer_art(i,jk,jb) = 1._wp
        ENDIF
      ENDDO
    ENDIF 
  ENDDO
END SUBROUTINE art_emiss_TR_art

SUBROUTINE art_emiss_TR_nin(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_nin)
!<
! SUBROUTINE art_emiss_TR_nin
! This subroutine sets the Northern Indian tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_nin(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_nin,area_nin)
END SUBROUTINE art_emiss_TR_nin

SUBROUTINE art_emiss_TR_sin(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_sin)
!<
! SUBROUTINE art_emiss_TR_sin
! This subroutine sets the Southern Indian tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_sin(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_sin,area_sin)
END SUBROUTINE art_emiss_TR_sin

SUBROUTINE art_emiss_TR_ech(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_ech)
!<
! SUBROUTINE art_emiss_TR_ech
! This subroutine sets the Eastern China tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_ech(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_ech,area_ech)
END SUBROUTINE art_emiss_TR_ech

SUBROUTINE art_emiss_TR_sea(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_sea)
!<
! SUBROUTINE art_emiss_TR_sea
! This subroutine sets the Southeast Asia tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_sea(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_sea,area_sea)
END SUBROUTINE art_emiss_TR_sea

SUBROUTINE art_emiss_TR_sib(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_sib)
!<
! SUBROUTINE art_emiss_TR_sib
! This subroutine sets the Siberia tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_sib(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_sib,area_sib)
END SUBROUTINE art_emiss_TR_sib

SUBROUTINE art_emiss_TR_eur(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_eur)
!<
! SUBROUTINE art_emiss_TR_eur
! This subroutine sets the Europe tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_eur(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_eur,area_eur)
END SUBROUTINE art_emiss_TR_eur

SUBROUTINE art_emiss_TR_naf(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_naf)
!<
! SUBROUTINE art_emiss_TR_naf
! This subroutine sets the Northern Africa tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_naf(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_naf,area_naf)
END SUBROUTINE art_emiss_TR_naf

SUBROUTINE art_emiss_TR_med(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_med)
!<
! SUBROUTINE art_emiss_TR_med
! This subroutine sets the Mediterranean tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_med(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_med,area_med)
END SUBROUTINE art_emiss_TR_med

SUBROUTINE art_emiss_TR_saf(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_saf)
!<
! SUBROUTINE art_emiss_TR_saf
! This subroutine sets the Southern Africa tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_saf(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_saf,area_saf)
END SUBROUTINE art_emiss_TR_saf

SUBROUTINE art_emiss_TR_mdg(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_mdg)
!<
! SUBROUTINE art_emiss_TR_mdg
! This subroutine sets the Madagascar tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_mdg(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_mdg,area_mdg)
END SUBROUTINE art_emiss_TR_mdg

SUBROUTINE art_emiss_TR_aus(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_aus)
!<
! SUBROUTINE art_emiss_TR_aus
! This subroutine sets the Australia tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_aus(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_aus,area_aus)
END SUBROUTINE art_emiss_TR_aus

SUBROUTINE art_emiss_TR_nam(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_nam)
!<
! SUBROUTINE art_emiss_TR_nam
! This subroutine sets the North America tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_nam(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_nam,area_nam)
END SUBROUTINE art_emiss_TR_nam


SUBROUTINE art_emiss_TR_sam(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_sam)
!<
! SUBROUTINE art_emiss_TR_sam
! This subroutine sets the South America tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_sam(:,:,:)

  CALL art_emiss_regio(jg,jb,jcs,jce,nproma,nlev,p_tracer_sam,area_sam)
END SUBROUTINE art_emiss_TR_sam


SUBROUTINE art_emiss_TR_tpo(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_tpo)
!<
! SUBROUTINE art_emiss_TR_tpo
! This subroutine sets the Tropical Pacific tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_tpo(:,:,:)


  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

  art_atmo    => p_art_data(jg)%atmo

  DO i = jcs, jce
    ! Tropical Pacific Ocean (tpo) emission tracer
    IF ((is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tpo1) .OR.   &
      &  is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tpo2)) .AND. &
      ! check sam region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_sam)) .AND. &
      ! check nam region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_nam)) ) THEN
  
      DO jk = nlev-num_levs_emiss,nlev
        p_tracer_tpo(i,jk,jb) = 1.0_wp
      END DO
    ENDIF
  END DO

END SUBROUTINE art_emiss_TR_tpo

SUBROUTINE art_emiss_TR_tao(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_tao)
!<
! SUBROUTINE art_emiss_TR_tao
! This subroutine sets the Tropical Atlantic tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_tao(:,:,:)


  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

  art_atmo    => p_art_data(jg)%atmo

  DO i = jcs, jce
    ! Tropical Atlantic Ocean (tao) emission tracer
    IF ((is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tao) .AND. &
      ! check sam region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_sam)) .AND. &
      ! check nam region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_nam))) .AND. &
      ! check saf region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_saf)) .AND. &
      ! check naf region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_naf))) THEN

      DO jk = nlev-num_levs_emiss,nlev
        p_tracer_tao(i,jk,jb) = 1.0_wp
      END DO
    ENDIF
  END DO

END SUBROUTINE art_emiss_TR_tao

SUBROUTINE art_emiss_TR_tio(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_tio)
!<
! SUBROUTINE art_emiss_TR_tio
! This subroutine sets the Tropical Indian Ocean tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_tio(:,:,:)


  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields

  art_atmo    => p_art_data(jg)%atmo

  DO i = jcs, jce
    ! Tropical Atlantic Ocean (tao) emission tracer
    IF (is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_tio)  .AND. &
      ! check saf region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_saf))  .AND. &
      ! check naf region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_naf))  .AND. &
      ! check mdg region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_mdg))  .AND. &
      ! check sin region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_sin))  .AND. &
      ! check sea region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_sea))  .AND. &
      ! check aus region
      &  (.NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_aus)) ) THEN

      DO jk = nlev-num_levs_emiss,nlev
        p_tracer_tio(i,jk,jb) = 1.0_wp
      END DO
    ENDIF
  END DO

END SUBROUTINE art_emiss_TR_tio


SUBROUTINE art_emiss_TR_bgs(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_bgs)
!<
! SUBROUTINE art_emiss_TR_bgs
! This subroutine sets the southern hemispheric background tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_bgs(:,:,:)

  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields
  TYPE(t_art_chem_indices),POINTER    :: &
     &  art_indices        !< Pointer to ART tracer indices

  art_atmo    => p_art_data(jg)%atmo
  art_indices => p_art_data(jg)%chem%indices

  DO i = jcs, jce
    DO jk = nlev-num_levs_emiss,nlev
      IF ((art_atmo%lat(i,jb) <=  0._wp) .AND. &
      ! check sam tracer
      &  ((art_indices%iTR_sam == 0) .OR.       &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_sam)) .AND. &
      ! check saf tracer
      &  ((art_indices%iTR_saf == 0) .OR.    &
      &   .NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_saf))  .AND. &
      ! check mdg tracer
      &  ((art_indices%iTR_mdg == 0) .OR.   &
      &   .NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_mdg))  .AND. &
      ! check sea tracer
      &  ((art_indices%iTR_sea == 0) .OR.   &
      &   .NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_sea))  .AND. &
      ! check aus tracer
      &  ((art_indices%iTR_aus == 0) .OR.   &
      &   .NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_aus)) .AND. &
      ! check tpo tracer
      &  ((art_indices%iTR_tpo == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tpo1) .AND. &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tpo2)) .AND. &
      ! check tao tracer
      &  ((art_indices%iTR_tao == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tao)) .AND. &
      ! check tio tracer
      &  ((art_indices%iTR_tio == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tio))  ) THEN
      
        p_tracer_bgs(i,jk,jb) = 1.0_wp
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE art_emiss_TR_bgs

SUBROUTINE art_emiss_TR_bgn(jg,jb,jcs,jce ,&
            &             nproma,nlev,p_tracer_bgn)
!<
! SUBROUTINE art_emiss_TR_bgn
! This subroutine sets the northern hemispheric background tracer
! Part of Module: mo_art_emission_regio_tracers
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!<

  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  REAL(wp), INTENT(INOUT) :: p_tracer_bgn(:,:,:)

  INTEGER :: &
    &  i, jk
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields
  TYPE(t_art_chem_indices),POINTER    :: &
     &  art_indices        !< Pointer to ART tracer indices

  art_atmo    => p_art_data(jg)%atmo
  art_indices => p_art_data(jg)%chem%indices

  DO i = jcs, jce
    DO jk = nlev-num_levs_emiss,nlev
      IF ((art_atmo%lat(i,jb) >=  0._wp) .AND.  &
      ! check sam tracer
      &  ((art_indices%iTR_sam == 0) .OR.       &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_sam)) .AND. &
      ! check sea tracer
      &  ((art_indices%iTR_sea == 0) .OR.   &
      &   .NOT. is_in_area(art_atmo%lat(i,jb),art_atmo%lon(i,jb), area_sea))  .AND. &
      ! check nam tracer
      &  ((art_indices%iTR_nam == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_nam)) .AND. &
      ! check eur tracer
      &  ((art_indices%iTR_eur == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_eur)) .AND. &
      ! check med tracer
      &  ((art_indices%iTR_med == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_med)) .AND. &
      ! check naf tracer
      &  ((art_indices%iTR_naf == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_naf)) .AND. &
      ! check sib tracer
      &  ((art_indices%iTR_sib == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_sib)) .AND. &
      ! check nin tracer
      &  ((art_indices%iTR_nin == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_nin)) .AND. &
      ! check sin tracer
      &  ((art_indices%iTR_sin == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_sin)) .AND. &
      ! check ech tracer
      &  ((art_indices%iTR_ech == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_ech)) .AND. &
      ! check tpo tracer
      &  ((art_indices%iTR_tpo == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tpo1) .AND. &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tpo2)) .AND. &
      ! check tao tracer
      &  ((art_indices%iTR_tao == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tao)) .AND. &
      ! check tio tracer
      &  ((art_indices%iTR_tio == 0) .OR.             &
      &  .NOT. is_in_area(art_atmo%lat(i,jb), art_atmo%lon(i,jb), area_tio)) ) THEN
        p_tracer_bgn(i,jk,jb) = 1.0_wp
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE art_emiss_TR_bgn


END MODULE mo_art_emission_regio_tracers
