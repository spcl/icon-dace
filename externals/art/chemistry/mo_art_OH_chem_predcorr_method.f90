!
! mo_art_OH_chem_predcorr_method
! This module provides routines for the execution of the predictor-corrector
! method for a sublist of tracers
! Based on Seinfeld and Pandis - Atmospheric Chemistry and Physics:
!                       (2012)   From Air Pollution to Climate Change
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

MODULE mo_art_OH_chem_predcorr_method
  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_physical_constants,            ONLY: amd, amo3
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  ! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN, IART_QV
  USE mo_art_chem_types_param,          ONLY: t_chem_meta_OH, OH_get_destruct
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_chem_data,                 ONLY: t_art_chem
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_OH_chem_types,             ONLY: t_art_OH_chem_meta
  USE mo_art_kinetic_constants,         ONLY: art_determine_bimolecular_kinetic_constant, &
                                          &   art_CH4_CO_kinetic_constants
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  
  PRIVATE

  PUBLIC ::                                   &
    &  art_calculate_OH_number_concentration, &
    &  art_predictor_corrector_method,        &
    &  art_OH_chem_set_tdep_meta,             &
    &  luse_gems_ozone

  LOGICAL, PARAMETER :: &
    &  luse_gems_ozone = .TRUE.
  CHARACTER(len=*), PARAMETER ::   &
    &  routine = 'mo_art_OH_chem_predcorr_method'

CONTAINS


!!
!!-----------------------------------------------------------------------------
!!
SUBROUTINE art_OH_chem_set_tdep_meta(jg,OH_chem_meta,temp,p_ozone,CH4_Nconc)
!<
! SUBROUTINE art_predcorr_set_tdep_meta
! This subroutine sets the time-dependent metadata of the predcorr_list
! For the user: if new tracer should be added to this mechanism, please include herein 
!               its reaction rate with OH
! Based on: 
! Part of Module: mo_art_OH_chem_predcorr_method
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-16
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  &
    &  jg     !< patch id on which computation is performed
  TYPE(t_art_OH_chem_meta), INTENT(inout) ::   &
    &  OH_chem_meta             !< general meta data of OH chemistry
  REAL(wp) , INTENT(in) ::    &
    &  temp(:,:,:)              !< current temperature
  REAL(wp), POINTER ::      &
    &  p_ozone(:,:,:)           !< pointer to ozone (either linoz or GEMS mmr)
  REAL(wp), INTENT(in), TARGET :: &
    &  CH4_Nconc(:,:,:)         !< number concentration of CH4 (# / cm3)

  !local variables
  INTEGER ::  &
    &  jc,jk,jb, i_startidx, i_endidx   !< loop indices
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  TYPE(t_art_chem), POINTER :: &
    &  art_chem


  art_atmo => p_art_data(jg)%atmo
  art_chem => p_art_data(jg)%chem

  IF (.NOT. ASSOCIATED(p_ozone) ) THEN
    CALL finish('mo_art_OH_chem_predcorr_method:art_predcorr_set_tdep_meta',  &
       &        'Ozone not found which is necessarily needed'                 &
       &      //' for predictor-corrector method.')
  END IF


  IF (luse_gems_ozone) THEN
    ! ozone in GEMS climatology is given as mass mixing ratio 
    ! and has to be converted to number concentration
    DO jb = art_atmo%i_startblk, art_atmo%i_endblk
      CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
      DO jk = 1,art_atmo%nlev
        DO jc = i_startidx,i_endidx
          OH_chem_meta%ozone_Nconc(jc,jk,jb) = p_ozone(jc,jk,jb) * amd / amo3           &
                &                              * art_chem%vmr2Nconc(jc,jk,jb)
        END DO
      END DO
    END DO
  ELSE
    DO jb = art_atmo%i_startblk, art_atmo%i_endblk
      CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
      DO jk = 1,art_atmo%nlev
        DO jc = i_startidx,i_endidx
          OH_chem_meta%ozone_Nconc(jc,jk,jb) = p_ozone(jc,jk,jb)
        END DO
      END DO
    END DO
  END IF

  
  CALL art_CH4_CO_kinetic_constants(jg,OH_chem_meta%k_CH4_OH,    &
                    &               OH_chem_meta%k_CO_OH,        &
                    &               art_chem%vmr2Nconc, temp)

  CALL art_calc_level_CH4_gt_1ppm(jg,OH_chem_meta%level_CH4_gt_1ppm,  &
                       &          CH4_Nconc, art_chem%vmr2Nconc)
  


END SUBROUTINE art_OH_chem_set_tdep_meta

!!
!!-----------------------------------------------------------------------------
!!

SUBROUTINE art_calc_level_CH4_gt_1ppm(jg,level_CH4_gt_1ppm,CH4_Nconc,vmr2Nconc)
!<
! SUBROUTINE art_calc_level_CH4_gt_1ppm
! This subroutine calculates the level index where CH4 gets lower than 1ppmv the
! first time
! Part of Module: mo_art_OH_chem_predcorr_method
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
!>
  INTEGER, INTENT(in) :: &
    &  jg                      !< patch%id
  INTEGER, INTENT(out) :: &
    &  level_CH4_gt_1ppm(:,:)  !< level where CH4 gets greater than 1ppmv
  REAL(wp), INTENT(in) :: &
    &  CH4_Nconc(:,:,:),  &    !< number concentration of CH4 (cm-3)
    & vmr2Nconc(:,:,:)         !< number concentration (cm-3) to volume mixing ratio (mol/mol)

  ! local variables
  INTEGER :: &
    &  i_startidx,i_endidx, jc,jk,jb
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
    
    DO jc = i_startidx,i_endidx
      jk = art_atmo%nlev
      DO WHILE ((jk > 0) .AND. (CH4_Nconc(jc,jk,jb) / vmr2Nconc(jc,jk,jb) > 1.e-6_wp))
        jk = jk - 1
      END DO

      IF (jk < art_atmo%nlev) THEN
        level_CH4_gt_1ppm(jc,jb) = jk + 1
      ELSE
        level_CH4_gt_1ppm(jc,jb) = art_atmo%nlev
      END IF
    END DO
  END DO

END SUBROUTINE art_calc_level_CH4_gt_1ppm
!!
!!-----------------------------------------------------------------------------
!!

SUBROUTINE art_calculate_OH_number_concentration(jg,OH_chem_meta,temp,J_O3, &
                             &                   CH4_Nconc, CO_Nconc)
!<
! SUBROUTINE art_calculate_OH_number_concentration
! This subroutine calculates diagnostically the number concentration of OH via photolysis of 
! ozone, reaction of O(1D) with water vapour. Depletion is performed with reaction of OH with
! CH4 and CO, respectively.
! Based on: e.g. Jacob (1999), Introduction to atmospheric chemistry
! Part of Module: mo_art_OH_chem_predcorr_method
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
!>
  IMPLICIT NONE
  INTEGER,       INTENT(in) :: &
    & jg                    !< patch id on which computation is performed
  REAL(wp), INTENT(in) :: &
    & temp(:,:,:)           !< current temperature 
  REAL(wp), INTENT(in) :: &
    & J_O3(:,:,:)           !< photolysis rate of O3 to O(1D)  in [1 / s]
  TYPE(t_art_OH_chem_meta), INTENT(inout) :: &
    &  OH_chem_meta        !< meta data of OH chemistry
  REAL(wp), INTENT(in) :: &
    &  CH4_Nconc(:,:,:),    & !< number concentrations of CH4 and CO (# / cm3)
    &  CO_Nconc(:,:,:)

  ! local variables
  REAL(wp) ::  &
    &  H2O_Nconc, O3_Nconc, & !< number concentrations of the substances
    &  N2_Nconc, O2_Nconc,  &
    &  k_O1D_H2O, k_O1D_O2, & !< kinetic constants (cm3/s)
    &  k_O1D_N2 
  INTEGER :: &
    &  jc,jk,jb, i_startidx, i_endidx
  TYPE(t_art_chem), POINTER :: &
    &  art_chem
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo => p_art_data(jg)%atmo
  art_chem => p_art_data(jg)%chem

  
  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
    DO jk=1,art_atmo%nlev
      DO jc=i_startidx,i_endidx

        N2_Nconc   =  0.78084_wp * art_chem%vmr2Nconc(jc,jk,jb)
        O2_Nconc   =  0.20946_wp * art_chem%vmr2Nconc(jc,jk,jb)
        H2O_Nconc  =  art_chem%water_tracers(jc,jk,jb,IART_QV)
        O3_Nconc   =  OH_chem_meta%ozone_Nconc(jc,jk,jb)

        CALL art_determine_bimolecular_kinetic_constant  &
               &  (k_O1D_H2O,temp(jc,jk,jb),1.63e-10_wp, -60._wp)
        CALL art_determine_bimolecular_kinetic_constant  &
               &  (k_O1D_O2, temp(jc,jk,jb),3.30e-11_wp, -55._wp)
        CALL art_determine_bimolecular_kinetic_constant  &
               &  (k_O1D_N2, temp(jc,jk,jb),2.15e-11_wp,-110._wp)


        OH_chem_meta%OH_Nconc(jc,jk,jb) = 2._wp * J_O3(jc,jk,jb)  &
                                                      &  * k_O1D_H2O * O3_Nconc * H2O_Nconc   &
                     &  / (k_O1D_O2 * O2_Nconc + k_O1D_N2 * N2_Nconc + k_O1D_H2O * H2O_Nconc) &
                     &  / (OH_chem_meta%k_CH4_OH(jc,jk,jb) * CH4_Nconc(jc,jk,jb)              &
                     &   + OH_chem_meta%k_CO_OH(jc,jk,jb) * CO_Nconc(jc,jk,jb))
      END DO
    END DO
  END DO

END SUBROUTINE art_calculate_OH_number_concentration  

!
! -------------------------------------------------------------------------------------------
!

SUBROUTINE art_predictor_corrector_method(jg,OH_chem_meta,p_prog_list,dtime)
!<
! SUBROUTINE art_predictor_corrector_method
! This subroutine performes the predictor-corrector method for the tracers with
! <c_solve>OH</c_solve> in the XML
! Based on: Seinfeld and Pandis (2012), equations 25.118 and 25.119
! Part of Module: mo_art_OH_chem_predcorr_method
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
!>
  IMPLiCIT NONE
  INTEGER, INTENT(in)       ::  &
    &  jg                    !< patch id on which computation is performed
  TYPE(t_art_OH_chem_meta), INTENT(inout) ::   &
    &  OH_chem_meta          !< general meta data of OH chemistry
  TYPE(t_var_list_ptr), INTENT(inout) :: &
    &  p_prog_list           !< list of prognistic variables
  REAL(wp), INTENT(in) ::     &
    &  dtime                 !< model time step
  
  ! local variables
  REAL(wp) ::                    &
    &  lifetime_star,            &
    &  lt_tr                      !< temporal saves of tracer lifetime and tracer production
  
  INTEGER ::                         &
    &  jb, jk, jc, iv,               &  !< loop indices
    &  i_startidx, i_endidx             !< loop indices
  INTEGER, POINTER ::    &
    &  jsp                 !< index of the tracer
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                     !< Pointer to ART atmo fields
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem                     !< Pointer to ART chem fields
  CHARACTER(:), ALLOCATABLE   :: &
    &  tracer_name                  !< name of the tracer

  TYPE(t_chem_meta_OH),POINTER ::  meta_ptr


  art_chem => p_art_data(jg)%chem
  art_atmo => p_art_data(jg)%atmo


  ! Calculate lifetimes from previous determined destruction rate 
  ! and take care of values near zero 
  DO iv = 1, p_prog_list%p%nvars
    jsp => p_prog_list%p%vl(iv)%p%info%ncontained
    SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      TYPE IS (t_chem_meta_OH)
        meta_ptr => meta
        DO jb = art_atmo%i_startblk,art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
      
          DO jk=1,art_atmo%nlev        
            DO jc=i_startidx,i_endidx
              ! if lifetime is too low, predictor/corrector method gets instable
              ! --> use the implicit method in this case
              lt_tr = 1._wp / meta%des_3d(jc,jk,jb)
              IF ( lt_tr < dtime) THEN
                  
                meta%tracer(jc,jk,jb)                                          &
                  & = meta%prod(jc,jk,jb) * lt_tr                              &
                  & + ( meta%tracer(jc,jk,jb) - meta%prod(jc,jk,jb) * lt_tr )  &
                  & * exp( - dtime / lt_tr)
              ELSE
                meta%tracer_star(jc,jk,jb)                               &
                  & = (meta%tracer(jc,jk,jb) * (2._wp * lt_tr - dtime)   &
                  &     +  2._wp * dtime * lt_tr * meta%prod(jc,jk,jb))  &
                  & / (2._wp * lt_tr + dtime)
              END IF
            ENDDO
          ENDDO
        ENDDO

        CALL meta%get_tracer_name(tracer_name)

        IF (TRIM(tracer_name) == 'TRCO') THEN
          OH_chem_meta%CO_star => meta_ptr%tracer_star
        ELSE IF (TRIM(tracer_name) == 'TRCH4') THEN
          OH_chem_meta%CH4_star => meta_ptr%tracer_star
        END IF
    END SELECT

    NULLIFY(jsp)
  

  ENDDO ! tracers
  
  ! calculate OH number concentration with star values
  CALL art_calculate_OH_number_concentration(jg,OH_chem_meta,art_atmo%temp,  &
                &                             art_chem%photo%rate(:,:,:,4),  &
                &                             OH_chem_meta%CH4_star(:,:,:),  &
                &                             OH_chem_meta%CO_star(:,:,:))
  
  DO iv = 1, p_prog_list%p%nvars
    SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      TYPE IS (t_chem_meta_OH)
        CALL OH_get_destruct(meta,jg,OH_chem_meta, .TRUE.)
    END SELECT

  END DO

  
  DO iv = 1, p_prog_list%p%nvars
    SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      TYPE IS (t_chem_meta_OH)
        CALL meta%get_prod_star()

        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
  
          DO jk = 1,art_atmo%nlev
            DO jc = i_startidx,i_endidx
              lt_tr = 1._wp / meta%des_3d(jc,jk,jb)
              lifetime_star = 1._wp / meta%des_star(jc,jk,jb)


              ! Corrector step
              IF (lt_tr >= dtime) THEN
                meta%tracer(jc,jk,jb) = ( meta%tracer(jc,jk,jb)                                 &
                                       &  * ( lifetime_star + lt_tr - dtime)                    &
                                       &  + 0.5_wp * dtime                                      &
                                       &    * (meta%prod_star(jc,jk,jb) + meta%prod(jc,jk,jb))  &
                                       &  * (lifetime_star + lt_tr)        )                    &
                                       &  / (lifetime_star + lt_tr + dtime )
              END IF
            END DO
          END DO
        END DO
    END SELECT


  END DO

END SUBROUTINE art_predictor_corrector_method
END MODULE mo_art_OH_chem_predcorr_method
