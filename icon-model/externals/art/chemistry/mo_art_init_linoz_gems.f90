!
! mo_art_init_linoz_gems
! Contains subroutines to calculate ozone climatology depending on day of the year.
! Taken from DWD's GME.
! Also routines overtaken from AES providing vertical interpolation between
! obs-data and model levels (another time interpolation will be added later)
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

MODULE mo_art_init_linoz_gems
    !<
    ! SUBROUTINE calc_o3_gems_linoz
    ! Calculates GEMS ozone climatology.
    ! Taken and adapted from ECMWF's IFS (37r2).
    ! Part of Module: mo_art_init_linoz_gems
    ! Author: Thorsten Reinhardt, AGeoBw, Offenbach
    ! Initial Release: 2011-10-18
    ! Modifications: Christian Stassen (KIT)
    !                ext_data -> o3_profile	20.03.2015
    !                return o3 VMR instead of MMR 07.05.2015
    !>

  USE mo_exception,              ONLY: finish
  USE mo_parallel_config,        ONLY: nproma
  USE mo_kind,                   ONLY: wp
  USE mo_math_constants,         ONLY: deg2rad,rad2deg
  USE mo_model_domain,           ONLY: t_patch
  USE mo_physical_constants,     ONLY: amd,amo3
  USE mo_o3_gems_data,           ONLY: RGHG7

  USE mtime,                     ONLY: datetime

  USE mo_art_data,               ONLY: p_art_data
  USE mo_art_atmo_data,          ONLY: t_art_atmo
  USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: calc_o3_gems_linoz


CONTAINS

!=======================================================================

!>
!! Calculates GEMS ozone climatology.
!! Taken and adapted from ECMWF's IFS (37r2).
!!
!! @par Revision History
!! Initial Release by Thorsten Reinhardt, AGeoBw, Offenbach (2011-10-18)
!!
!! Changelog :: - Christian Stassen (CS) ext_data -> o3_profile	20.03.2015
!!              -                   (CS) return o3 VMR instead of MMR 07.05.2015
SUBROUTINE calc_o3_gems_linoz(pt_patch,current_date,o3_profile)


    

  CHARACTER(len=*), PARAMETER :: routine =  'calc_o3_gems'

  INTEGER , PARAMETER :: ilat=64,nlev_gems=91

  ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
  REAL(wp), PARAMETER     :: zrefp(nlev_gems) = (/ &
    & 0.10000000E+01_wp, 0.29900000E+01_wp, 0.56834998E+01_wp, 0.10147500E+02_wp, &
    & 0.17160500E+02_wp, 0.27682501E+02_wp, 0.42848999E+02_wp, 0.63956501E+02_wp, &
    & 0.92440987E+02_wp, 0.12985049E+03_wp, 0.17781149E+03_wp, 0.23799651E+03_wp, &
    & 0.31209000E+03_wp, 0.40175449E+03_wp, 0.50860199E+03_wp, 0.63416602E+03_wp, &
    & 0.77987903E+03_wp, 0.94705548E+03_wp, 0.11368750E+04_wp, 0.13503740E+04_wp, &
    & 0.15884360E+04_wp, 0.18517920E+04_wp, 0.21410149E+04_wp, 0.24565259E+04_wp, &
    & 0.27986001E+04_wp, 0.31673640E+04_wp, 0.35628110E+04_wp, 0.39848059E+04_wp, &
    & 0.44330962E+04_wp, 0.49073169E+04_wp, 0.54070068E+04_wp, 0.59314971E+04_wp, &
    & 0.64797832E+04_wp, 0.70505981E+04_wp, 0.76428970E+04_wp, 0.82572246E+04_wp, &
    & 0.88957646E+04_wp, 0.95614326E+04_wp, 0.10257570E+05_wp, 0.10988080E+05_wp, &
    & 0.11757620E+05_wp, 0.12571580E+05_wp, 0.13435160E+05_wp, 0.14352070E+05_wp, &
    & 0.15325200E+05_wp, 0.16357520E+05_wp, 0.17452211E+05_wp, 0.18612539E+05_wp, &
    & 0.19841900E+05_wp, 0.21144000E+05_wp, 0.22522650E+05_wp, 0.23981760E+05_wp, &
    & 0.25525529E+05_wp, 0.27158301E+05_wp, 0.28884641E+05_wp, 0.30709301E+05_wp, &
    & 0.32637240E+05_wp, 0.34673641E+05_wp, 0.36823859E+05_wp, 0.39093570E+05_wp, &
    & 0.41487949E+05_wp, 0.44010520E+05_wp, 0.46651148E+05_wp, 0.49386160E+05_wp, &
    & 0.52190051E+05_wp, 0.55035488E+05_wp, 0.57894629E+05_wp, 0.60745699E+05_wp, &
    & 0.63574230E+05_wp, 0.66367703E+05_wp, 0.69114523E+05_wp, 0.71801031E+05_wp, &
    & 0.74412922E+05_wp, 0.76938641E+05_wp, 0.79368109E+05_wp, 0.81689688E+05_wp, &
    & 0.83891563E+05_wp, 0.85965023E+05_wp, 0.87903594E+05_wp, 0.89700203E+05_wp, &
    & 0.91347422E+05_wp, 0.92841422E+05_wp, 0.94181656E+05_wp, 0.95367977E+05_wp, &
    & 0.96399797E+05_wp, 0.97280203E+05_wp, 0.98034609E+05_wp, 0.98677672E+05_wp, &
    & 0.99198242E+05_wp, 0.99594977E+05_wp, 0.99881500E+05_wp /)

  ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
  REAL(wp), PARAMETER     :: zmday(12) = (/ &
    &  31._wp,    59.25_wp,  90.25_wp, 120.25_wp, 151.25_wp, 181.25_wp, &
    & 212.25_wp, 243.25_wp, 273.25_wp, 304.25_wp, 334.25_wp, 365.25_wp /)

  ! Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).
  REAL(wp), PARAMETER     :: zytime(12) = (/ &
    &  22320._wp,  64980._wp, 107640._wp, 151560._wp, 195480._wp, 239400._wp, &
    & 283320._wp, 327960._wp, 371880._wp, 415800._wp, 459720._wp, 503640._wp /)    

  TYPE(t_patch),         INTENT(in) :: pt_patch    ! Patch
  TYPE(datetime),POINTER,INTENT(in) :: current_date

  REAL(wp), INTENT(inout)   :: o3_profile(:,:,:)  !! Profile of ozone mmr

  ! local fields
  INTEGER  :: idx0(nproma,0:pt_patch%nlev,pt_patch%nblks_c)
  REAL(wp) :: zlat(0:ilat+1)
  REAL(wp) :: zozn(0:ilat+1,1:nlev_gems)
  REAL(wp) :: zpresh(0:nlev_gems)
  REAL(wp) :: rclpr(0:nlev_gems)
  REAL(wp) :: zo3(nproma,1:nlev_gems,pt_patch%nblks_c)
  REAL(wp) :: zviozo(nproma,0:pt_patch%nlev)
  REAL(wp) :: zozovi(nproma,0:nlev_gems,pt_patch%nblks_c)
  LOGICAL  :: l_found(nproma)

  ! local scalars
  INTEGER  :: jk,jkk,jl,jc,jb !loop indices
  INTEGER  :: idy,im,imn,im1,im2,jk_start,i_startidx,i_endidx
  REAL(wp) :: ztimi,zxtime,zjl,zlatint,zint
  LOGICAL  :: lfound_all

  INTEGER :: jg
  TYPE(t_art_atmo),POINTER    :: &
     &  art_atmo           !< Pointer to ART atmospheric fields


  jg = pt_patch%id

  art_atmo => p_art_data(jg)%atmo
  

  !Time index. Taken from su_ghgclim.F90 of ECMWF's IFS (37r2).

  IDY = current_date%date%day - 1 !NDD(KINDAT)-1
  IMN = current_date%date%month ! NMM(KINDAT)
  IF (IMN == 1) THEN
    ZXTIME=REAL(IDY,wp)*1440._wp + current_date%time%minute !KMINUT
  ELSEIF (IMN == 2) THEN
    IF(IDY == 28) IDY=IDY-1
    ! A DAY IN FEB. IS 28.25*24*60/28=1452.8571min LONG.
    ZXTIME=44640._wp+REAL(IDY,wp)*1452.8571_wp + current_date%time%minute !KMINUT
  ELSE
    ZXTIME=(ZMDAY(IMN-1)+REAL(IDY,KIND(ZXTIME)))*1440._wp + current_date%time%minute !KMINUT
  ENDIF
  ! 525960=MINUTES IN A SIDERAL YEAR (365.25d)
  ZXTIME=MOD(ZXTIME,525960._wp)

  IM1=0
  IM2=0
  IF (ZXTIME <= ZYTIME(1)) THEN
    IM1=12
    IM2=1
    ZTIMI=(ZYTIME(1)-ZXTIME)/44640._wp
  ELSEIF(ZXTIME > ZYTIME(12)) THEN
    IM1=12
    IM2=1
    ZTIMI=(548280._wp-ZXTIME)/44640._wp
    ! 548280.=(365.25d + 15.5d)*24*60
  ELSE
    DO IM=1,11
      IF (ZXTIME > ZYTIME(IM) .AND. ZXTIME <= ZYTIME(IM+1)) THEN
        IM1=IM
        IM2=IM+1
        ZTIMI=(ZXTIME-ZYTIME(IM2))/(ZYTIME(IM1)-ZYTIME(IM2))
      ENDIF
    ENDDO
    IF ( IM1 == 0 .OR. IM2 == 0 ) THEN
      CALL finish(TRIM(routine),'Problem with time interpolation in suecozc!')
    ENDIF
  ENDIF

!==================================================================================

  !---------------------------------------------------------------------------
  ! 2023-11-20
  ! Quick hack for initializiation of linoz ozone
  ! Reason: pres_ifc and other diagnostic variables in new icon version not set at this stage
  ! Solution: assume a standard altitude pressure relation
  ! Possible problem: total ozone is not strictly conserved by the interpolation

  ! ini o3_profile 
  o3_profile = 1.e-8_wp

  ! time interpolation of mixing ratio
  DO jk=1,nlev_gems
    DO jl=1,ilat
      zozn(JL,JK) = (RGHG7(JL,JK,IM2)&
        & +ZTIMI*(RGHG7(JL,JK,IM1)-RGHG7(JL,JK,IM2)))
    ENDDO
  ENDDO

  !--- prepare latitude interpolation
  DO jk=1,nlev_gems
    zozn(0,JK)      = zozn(1,jk)
    zozn(ilat+1,jk) = zozn(ilat,jk)
  ENDDO

  zlatint=180._wp/REAL(ilat,wp)
  DO jl=0,ilat+1
    zlat(jl)=(-90._wp+0.5_wp*zlatint+(jl-1)*zlatint)*deg2rad
  ENDDO
  
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk

    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    !DO jc=i_startidx,i_endidx
    !  !--- check 
    !  IF ( minval( art_atmo%pres(jc,:,jb) ) < 0._wp  .OR. &
    !    & maxval( art_atmo%pres(jc,:,jb) ) <= 0._wp ) THEN
!
    !    print *, 'min/max pres:', jb,  jc, minval( art_atmo%pres(jc,:,jb) ) &
    !        & , maxval( art_atmo%pres(jc,:,jb) )
!
    !  ENDIF
    !ENDDO

    ! Latitude interpolation
    DO jkk=1,nlev_gems
      DO jc=i_startidx,i_endidx

        zjl=1._wp+(rad2deg*pt_patch%cells%center(jc,jb)%lat+90._wp-.5_wp*zlatint)/zlatint

        !linear interpolation
        zo3(jc,jkk,jb) = zozn(INT(zjl),jkk) &
          & + (zozn(INT(zjl)+1,jkk)-zozn(INT(zjl),jkk)) &
          &  /  (zlat(INT(zjl)) - zlat(INT(zjl)+1)) &
          &  * (zlat(INT(zjl)) - pt_patch%cells%center(jc,jb)%lat)

      ENDDO !jc
    ENDDO !jkk

    ! Vertical interpolation; p ordered by increasing p, not guarantedd for non-hydrostatic
    ! mass mixing ratio
    DO jc = i_startidx,i_endidx
      jkk = 1
      DO jk = 1,pt_patch%nlev
        ! for pressure > 100 mb use default value
        IF ( art_atmo%pres(jc,jk,jb) > 10000._wp ) cycle
        DO WHILE ( art_atmo%pres(jc,jk,jb) > zrefp(jkk+1) .AND. jkk < nlev_gems ) 
           jkk = jkk + 1
        ENDDO
        o3_profile(jc,jk,jb) = amo3/amd*max( zo3(jc,jkk,jb) + ( zo3(jc,jkk+1,jb) - zo3(jc,jkk,jb) ) &
                & * ( art_atmo%pres(jc,jk,jb) - zrefp(jkk) )/(zrefp(jkk+1) - zrefp(jkk) ), 0._wp )
      ENDDO
      !if ( minval(o3_profile(jc,:,jb)) < 0._wp ) then
      !  print *, 'Negative Ozone:', jc, jb
      !  do jk=1,pt_patch%nlev
      !    print *, jk, -7.*log( art_atmo%pres(jc,jk,jb)/101325._wp ), o3_profile(jc,jk,jb)
      !  enddo
      !endif
           
    ENDDO

  !CALL finish(TRIM(routine),'Stop in Test')

  ENDDO !jb


  RETURN

  !--- End of quick hack
  !------------------------------------------------------------------------------------------
!=============================================================================================

  ! Pressure levels of climatological ozone fields. From su_ghgclim.F90 of ECMWF's IFS (37r2).

  ZPRESH(0)=0.0_wp
  RCLPR(0) =0.0_wp
  DO JK=1,NLEV_GEMS-1
    ZPRESH(JK)=(ZREFP(JK)+ZREFP(JK+1))*0.5_wp
    RCLPR(JK) =ZPRESH(JK)
  ENDDO
  ZPRESH(NLEV_GEMS)=110000._wp
  RCLPR(nlev_gems) =ZPRESH(nlev_gems)


  ! volume mixing ratio to ozone pressure thickness
  DO jk=1,nlev_gems
    DO jl=1,ilat
      zozn(JL,JK) = amo3/amd * (RGHG7(JL,JK,IM2)&
        & +ZTIMI*(RGHG7(JL,JK,IM1)-RGHG7(JL,JK,IM2)))
      zozn(JL,JK) = zozn(JL,JK) * (ZPRESH(JK)-ZPRESH(JK-1))
    ENDDO
  ENDDO

  DO jk=1,nlev_gems
    zozn(0,JK)      = zozn(1,jk)
    zozn(ilat+1,jk) = zozn(ilat,jk)
  ENDDO

  !Preparations for latitude interpolations

  zlatint=180._wp/REAL(ilat,wp)

  DO jl=0,ilat+1
    zlat(jl)=(-90._wp+0.5_wp*zlatint+(jl-1)*zlatint)*deg2rad
  ENDDO
  

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,jkk,i_startidx,i_endidx,zjl,jk_start,l_found,lfound_all,&
!$OMP zint,zviozo) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk

    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    ! Latitude interpolation
    
    DO jkk=1,nlev_gems

      DO jc=i_startidx,i_endidx

        zjl=1._wp+(rad2deg*pt_patch%cells%center(jc,jb)%lat+90._wp-.5_wp*zlatint)/zlatint

        !just select nearest value
        !zo3(jc,jkk,jb) = zozn(NINT(zjl),jkk)

        !linear interpolation
        zo3(jc,jkk,jb) = zozn(INT(zjl),jkk) &
          & + (zozn(INT(zjl)+1,jkk)-zozn(INT(zjl),jkk)) &
          &  /  (zlat(INT(zjl)) - zlat(INT(zjl)+1)) &
          &  * (zlat(INT(zjl)) - pt_patch%cells%center(jc,jb)%lat)

      ENDDO !jc

    ENDDO !jkk

    ! ACCUMULATE FROM TOP TO BOTTOM THE LATITUDE INTERPOLATED FIELDS
    ! From radghg.F90 of ECMWF's IFS.

    zozovi(i_startidx:i_endidx,0,jb) = 0._wp
    DO jkk=1,nlev_gems
      zozovi(i_startidx:i_endidx,jkk,jb) = zozovi(i_startidx:i_endidx,jkk-1,jb) &
        &                                 + zo3(i_startidx:i_endidx,jkk,jb)
    ENDDO

    ! REDISTRIBUTE THE VERTIC. INTEGR. CLIM. O3 ON THE MODEL GRID
    ! Adapted from radghg.F90 of ECMWF's IFS.

    jk = 0
    zviozo(i_startidx:i_endidx,0) = 0._wp

    jk_start = 0

    DO jk = 0,pt_patch%nlev
      l_found(:) = .FALSE.
      lfound_all = .FALSE.
      DO jkk = jk_start,nlev_gems-1
        DO jc = i_startidx,i_endidx
          IF( art_atmo%pres_ifc(jc,jk+1,jb) >= RCLPR(jkk)  &
            & .AND. art_atmo%pres_ifc(jc,jk+1,jb) < RCLPR(jkk+1)) THEN
            ZINT=(art_atmo%pres_ifc(jc,jk+1,jb)-RCLPR(jkk))/(RCLPR(jkk+1)-RCLPR(jkk)) 
            ZVIOZO(jc,JK) = ZOZOVI(jc,jkk,jb) + ZINT * (ZOZOVI(jc,jkk+1,jb)-ZOZOVI(jc,jkk,jb))
            l_found(jc) = .TRUE.
            idx0(jc,jk,jb) = jkk
          ELSEIF ( art_atmo%pres_ifc(jc,jk+1,jb) > RCLPR(nlev_gems) ) THEN
            l_found(jc) = .TRUE.
            idx0(jc,jk,jb) = nlev_gems
          ENDIF
        ENDDO !jc
        IF (ALL(l_found(i_startidx:i_endidx))) THEN
          lfound_all = .TRUE.
          EXIT
        ENDIF
      ENDDO !jkk
      IF (lfound_all) THEN
        jk_start = MIN(MINVAL(idx0(i_startidx:i_endidx,jk,jb)),nlev_gems-1)
      ENDIF
    ENDDO !jk

    ! COMPUTE THE MASS MIXING RATIO

    DO jk = 1,pt_patch%nlev
      DO jc = i_startidx,i_endidx
        o3_profile(jc,jk,jb)=(ZVIOZO(jc,jk)-ZVIOZO(jc,jk-1)) / art_atmo%dpres_mc(jc,jk,jb)

! MK, 2014-01: tuning to increase ozone by maximum 50% at 15km decreasing linearly
!              towards 0% at 30km and 10km; goal: warmer lower stratosphere 0.5K
!
! >>> No tuning for init.
!          IF ( ltuning_ozone ) THEN
!            o3_profile(jc,jk,jb) = o3_profile(jc,jk,jb) &
!                                    & * atm_phy_nwp_config(pt_patch%id)%fac_ozone(jk)
!          ENDIF

        !IF ( jk > 30 .and. jk < 55 ) THEN   ! triangle in [30,55] with maximum at L50
        !  o3_profile(jc,jk,jb) = o3_profile(jc,jk,jb) &
        !    * ( 1.0_wp + 0.5_wp * min( (jk-30)/20._wp, (55-jk)/5._wp ) )
        !ENDIF
      ENDDO
    ENDDO

  ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL  


END SUBROUTINE calc_o3_gems_linoz
END MODULE mo_art_init_linoz_gems

