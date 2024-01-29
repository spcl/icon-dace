!
! mo_art_simnoy
! This module provides a simplified photochemistry for
! N2O and NOy
! first introduced by Olsen, 2001
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

MODULE mo_art_simnoy

  ! ICON
  USE mo_kind,                 ONLY: wp
  USE mo_math_constants,       ONLY: pi
  USE mtime,                   ONLY: datetime
  
  ! ART
  USE mo_art_read_simnoy,      ONLY:  tparm_noy, lparm,    &
                                 &    kparm,               & 
                                 &    lparm_min, lparm_dist
  USE mo_art_data,             ONLY: p_art_data
  USE mo_art_atmo_data,        ONLY: t_art_atmo
  USE mo_art_chem_types_param, ONLY: t_chem_meta_simnoy

  IMPLICIT NONE

  PRIVATE   

  PUBLIC  ::  art_noy_polarchem, art_noy_polarchem_sedi
  PUBLIC  ::  art_n2onoy_chemistry_wmo, art_n2onoy_chemistry_pres, art_n2onoy_chemistry_extp


CONTAINS


! -----------------------------------------
! --- simplified N2O / NOy chemistry
! --- Olsen et al., 2001
! --- by IMK-ASF DATE
! -----------------------------------------

SUBROUTINE art_simnoy_get_tab_values(current_date,lat,z_mc,nlev, jcs, jce, &
                      &              n2onoy_tab1, n2onoy_tab2,             &
                      &              n2onoy_tab3, n2onoy_tab4)
  ! inout variables
  TYPE(datetime), INTENT(in)    ::  &
    &  current_date           !< actual date
  REAL(wp),       INTENT(in)    ::  &
    &  lat(:)           !< [lat] = rad
  REAL(wp),       INTENT(in)    ::  &
    &  z_mc(:,:)          !< geometric height of full levels
  INTEGER,        INTENT(in)    ::  &
    &  nlev, jcs, jce
  REAL(wp), INTENT(inout)         ::  &
    &  n2onoy_tab1(:,:),              &
    &  n2onoy_tab2(:,:),              &
    &  n2onoy_tab3(:,:),              &
    &  n2onoy_tab4(:,:)

  ! local variables
  INTEGER              :: jk,ik,il,ilp, jc
  REAL(wp)             :: lp
  REAL(wp)             :: alf,alf1


  lp = lparm_min - lparm_dist

  DO jc = jcs,jce
    DO jk=1,nlev

      ! nearest neighbour for lat
      ik = INT( (lat(jc) + pi/2.) / pi * kparm ) + 1  

      ! linear interpolation for lev
      il = INT( (z_mc(jc,jk)/1000. - lp) / lparm_dist)
      ilp = il + 1

      IF (il >= lparm) THEN
        il = lparm
        ilp = lparm
      ENDIF

      IF (il < 1) THEN
        il = 1
        ilp = 1
      ENDIF

      alf = ( z_mc(jc,jk)/1000._wp - FLOAT(il*INT(lparm_dist) + INT(lp) )) &
            & / lparm_dist

      IF ((alf > 1.0) .OR. (alf <= 0.0)) THEN
        alf = 0.0
        ilp = il
      ENDIF

      alf1 = 1.0 - alf

      ! creating interpolated coefficient lists
      n2onoy_tab1(jc,jk) = alf1 * tparm_noy(il,ik,current_date%date%month,1) &
                        & + alf * tparm_noy(ilp,ik,current_date%date%month,1)
      n2onoy_tab2(jc,jk) = alf1 * tparm_noy(il,ik,current_date%date%month,2) &
                        & + alf * tparm_noy(ilp,ik,current_date%date%month,2)
      n2onoy_tab3(jc,jk) = alf1 * tparm_noy(il,ik,current_date%date%month,3) &
                        & + alf * tparm_noy(ilp,ik,current_date%date%month,3)
      n2onoy_tab4(jc,jk) = alf1 * tparm_noy(il,ik,current_date%date%month,4) &
                        & + alf * tparm_noy(ilp,ik,current_date%date%month,4)
    ENDDO
  END DO

END SUBROUTINE art_simnoy_get_tab_values

SUBROUTINE art_n2onoy_chemistry_extp(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
! SUBROUTINE n2onoy_chemistry
! based on Olsen 2001
! Author: Christopher Diekmann, KIT
!>
  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime

  ! local variables
  REAL(wp)             :: dn2odt, noyvmr, &
                        & dnoydt, noyloss
  INTEGER              :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo

  nlev = art_atmo%nlev

  ! -------------------------------------------------
  ! interpolation of coefficients onto ICON grid
  ! -------------------------------------------------

  CALL art_simnoy_get_tab_values(current_date, art_atmo%lat(:,jb),   &
          &                      art_atmo%z_mc(:,:,jb), nlev,        &
          &                      jcs, jce, tracer%n2onoy_tab1, tracer%n2onoy_tab2, &
          &                      tracer%n2onoy_tab3, tracer%n2onoy_tab4)


  ! ------------------------------------------------------------
  ! start chemistry
  ! calculate new mixing ratios for NOy and N2O
  ! ------------------------------------------------------------

  DO jc = jcs,jce
    DO jk=1,nlev
      ! N2O loss:
      dn2odt = (1.0 - exp(-tracer%n2onoy_tab1(jc,jk)*p_dtime)) * tracer%n2o_tracer(jc,jk,jb)
      
      tracer%tend_n2o(jc,jk,jb) = - dn2odt

      ! NOy production:
      tracer%tend(jc,jk,jb) = 2.0 * dn2odt * tracer%n2onoy_tab3(jc,jk)

      ! NOy loss:
      noyvmr = tracer%tracer(jc,jk,jb) / p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

      dnoydt = (1.0 - exp(-tracer%n2onoy_tab2(jc,jk)*p_dtime)) * tracer%tracer(jc,jk,jb)
      noyloss = 2.0 * noyvmr * tracer%n2onoy_tab4(jc,jk)       &
              & / (noyvmr * tracer%n2onoy_tab4(jc,jk) + 1.0)

      tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb)  - noyloss*dnoydt
    ENDDO
  END DO
END SUBROUTINE art_n2onoy_chemistry_extp

SUBROUTINE art_n2onoy_chemistry_wmo(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
! SUBROUTINE n2onoy_chemistry
! based on Olsen 2001
! Author: Christopher Diekmann, KIT
!>
  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime

  ! local variables
  REAL(wp)             :: dn2odt, noyvmr, &
                        & dnoydt, noyloss
  INTEGER              :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo

  nlev = art_atmo%nlev

  ! -------------------------------------------------
  ! interpolation of coefficients onto ICON grid
  ! -------------------------------------------------

  CALL art_simnoy_get_tab_values(current_date, art_atmo%lat(:,jb),   &
          &                      art_atmo%z_mc(:,:,jb), nlev,        &
          &                      jcs, jce, tracer%n2onoy_tab1, tracer%n2onoy_tab2, &
          &                      tracer%n2onoy_tab3, tracer%n2onoy_tab4)


  ! ------------------------------------------------------------
  ! start chemistry
  ! calculate new mixing ratios for NOy and N2O
  ! ------------------------------------------------------------
  DO jc = jcs,jce
    DO jk = 1, art_atmo%ktrpwmo(jc,jb)
      IF (jk <= art_atmo%ktrpwmo(jc,jb)) THEN
        ! N2O loss:
        dn2odt = (1.0 - exp(-tracer%n2onoy_tab1(jc,jk)*p_dtime)) * tracer%n2o_tracer(jc,jk,jb)
        
        tracer%tend_n2o(jc,jk,jb) = - dn2odt

        ! NOy production:
        tracer%tend(jc,jk,jb) = 2.0 * dn2odt * tracer%n2onoy_tab3(jc,jk)

        ! NOy loss:
        noyvmr = tracer%tracer(jc,jk,jb) / p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

        dnoydt = (1.0 - exp(-tracer%n2onoy_tab2(jc,jk)*p_dtime)) * tracer%tracer(jc,jk,jb)

        noyloss = 2.0 * noyvmr * tracer%n2onoy_tab4(jc,jk)     &
                & / (noyvmr * tracer%n2onoy_tab4(jc,jk) + 1.0)

        tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) - noyloss*dnoydt
      ELSE 
        tracer%tend_n2o(jc,jk,jb) = tracer%n2o_tracer(jc,jk,jb)           &
                &                    * exp(-1._wp*p_dtime*tracer%des_N2O) &
                &                    - tracer%n2o_tracer(jc,jk,jb)

        tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb)           &
                &                * exp(-1._wp*p_dtime*tracer%des) &
                &                - tracer%tracer(jc,jk,jb)
      ENDIF
    ENDDO
  END DO
END SUBROUTINE art_n2onoy_chemistry_wmo

SUBROUTINE art_n2onoy_chemistry_pres(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!<
! SUBROUTINE n2onoy_chemistry
! based on Olsen 2001
! Author: Christopher Diekmann, KIT
!>
  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime

  ! local variables
  REAL(wp)             :: dn2odt, noyvmr, &
                        & dnoydt, noyloss
  INTEGER              :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo

  nlev = art_atmo%nlev

  ! -------------------------------------------------
  ! interpolation of coefficients onto ICON grid
  ! -------------------------------------------------

  CALL art_simnoy_get_tab_values(current_date, art_atmo%lat(:,jb),   &
          &                      art_atmo%z_mc(:,:,jb), nlev,        &
          &                      jcs, jce, tracer%n2onoy_tab1, tracer%n2onoy_tab2, &
          &                      tracer%n2onoy_tab3, tracer%n2onoy_tab4)


  ! ------------------------------------------------------------
  ! start chemistry
  ! calculate new mixing ratios for NOy and N2O
  ! ------------------------------------------------------------

  DO jc = jcs,jce
    DO jk = 1,nlev
      IF (art_atmo%pres(jc,jk,jb) <= 900.*100.) THEN 
        ! N2O loss:
        dn2odt = (1.0 - exp(-tracer%n2onoy_tab1(jc,jk)*p_dtime)) * tracer%n2o_tracer(jc,jk,jb)
        
        tracer%tend_n2o(jc,jk,jb) =  - dn2odt

        ! NOy production:
        tracer%tend(jc,jk,jb) =  2.0 * dn2odt * tracer%n2onoy_tab3(jc,jk)

        ! NOy loss:
        noyvmr = tracer%tracer(jc,jk,jb) / p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)

        dnoydt = (1.0 - exp(-tracer%n2onoy_tab2(jc,jk)*p_dtime)) * tracer%tracer(jc,jk,jb)
        noyloss = 2.0 * noyvmr * tracer%n2onoy_tab4(jc,jk)      &
                & / (noyvmr * tracer%n2onoy_tab4(jc,jk) + 1.0)

        tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) - noyloss*dnoydt

      ELSE
        tracer%tend_n2o(jc,jk,jb) = tracer%n2o_tracer(jc,jk,jb)           &
                &                    * exp(-1._wp*p_dtime*tracer%des_N2O) &
                &                    - tracer%n2o_tracer(jc,jk,jb)

        tracer%tend(jc,jk,jb) = tracer%tracer(jc,jk,jb)               &
                &                * exp(-1._wp*p_dtime * tracer%des)   &
                &                - tracer%tracer(jc,jk,jb)
      
      ENDIF
    ENDDO
  END DO

END SUBROUTINE art_n2onoy_chemistry_pres

SUBROUTINE art_noy_polarchem_sedi(jg,jb,jcs,jce ,&
            &             tracer, current_date, p_dtime )
!< 
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Author: Christopher Diekmann, KIT
!---------------------------------------------------------
!>

  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date
  REAL(wp), INTENT(in) :: &
    &  p_dtime


  ! Local variables
  REAL(wp) :: diff_n_noy, n_noy, n_noy_new
  REAL(wp) :: rate

  INTEGER         :: jk, nlev, jc

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo
  nlev = art_atmo%nlev

  !--------------------------------------------------------------------

  DO jc = jcs, jce
    diff_n_noy = 0._wp

    DO jk=1,nlev
      n_noy =  tracer%tracer(jc,jk,jb) * 1.e6_wp

      n_noy = n_noy + diff_n_noy

      IF (jk /= nlev) THEN 
          rate = tracer%des_noysed * tracer%cold_tracer(jc,jk,jb)  &
              &   /  p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb)
          n_noy_new = n_noy * exp(-1._wp*p_dtime * rate)

          diff_n_noy = n_noy - n_noy_new
      END IF
          
      tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) +  n_noy_new * 1.e-6_wp          &
           &                - tracer%tracer(jc,jk,jb)
    END DO
  END DO

END SUBROUTINE art_noy_polarchem_sedi

SUBROUTINE art_noy_polarchem(jg,jb,jcs,jce ,&
            &             nproma, nlev, tracer)
!< 
!---------------------------------------------------------
!--- Calculate non-linear heterogenous ozone loss for polar areas
!--- Author: Christopher Diekmann, KIT
!---------------------------------------------------------
!>

  ! inout variables
  INTEGER, INTENT(in) :: &
    &  jg,jb,jcs,jce,    &
    &  nproma, nlev
  TYPE(t_chem_meta_simnoy), INTENT(inout) :: &
    &  tracer

  ! Local variables
  REAL(wp) :: diff_n_noy, n_noy, n_noy_new
  REAL(wp) :: diff_vmr_cold

  INTEGER  :: jk, jc

  !--------------------------------------------------------------------
  DO jc = jcs, jce
    diff_n_noy  = 0._wp
  
    DO jk=1,nlev
      ! vertical redistribution of noy
      n_noy = tracer%tracer(jc,jk,jb) * 1.e6_wp

      n_noy = n_noy + diff_n_noy
  
      IF (jk /= nlev) THEN 
        !diff_vmr_cold = tracer%p_cold_sed(jc,jk,jb) * 1.e-6_wp       &
        !              &  /  p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb) 
        ! The upper part is not restartable, due to p_cold_sed not being
        ! properly recreated at the restart (needs some code restructuring)
        ! AS A PRELIMARY SOLUTION P_COLD_SED IS SET TO A CONSTANT 0.0
        ! THIS BASICALLY TURNS THIS PROCESS OFF. THIS NEEDS PROPER
        ! TREATMENT BEFORE IT CAN BE ACTIVATED AGAIN.
        ! IF YOU CAN READ THIS MESSAGE, THIS WAS NOT YET DONE.
        diff_vmr_cold = 0.0_wp * 1.e-6_wp       &
                      &  /  p_art_data(jg)%chem%vmr2Nconc(jc,jk,jb) 
        n_noy_new     = n_noy * (1._wp - diff_vmr_cold)
  
        diff_n_noy    = n_noy - n_noy_new
      END IF
      
      tracer%tend(jc,jk,jb) = tracer%tend(jc,jk,jb) + n_noy_new * 1.e-6_wp   &
           &                  - tracer%tracer(jc,jk,jb)
    ENDDO
  END DO

END SUBROUTINE art_noy_polarchem

END MODULE mo_art_simnoy
