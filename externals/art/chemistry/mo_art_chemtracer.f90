!
! mo_art_chemtracer
! This module provides the subroutines for chemical tracer processes.
! Loss processes are based on simple parametrizations regarding to life time
! and / or linearization
! This Module is only envoked if lart_chemtracer == .TRUE.
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

MODULE mo_art_chemtracer

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic, t_var_metadata

  USE mtime,                            ONLY: datetime

! ART
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_chem_data,                 ONLY: t_art_chem,         &
                                          &   t_art_chem_indices, &
                                          &   t_art_chem_param
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_impl_constants,            ONLY: IART_LINOZ_ANA,        &
                                          &   IART_LINOZ_LT,         &
                                          &   IART_POLARCHEM,        &
                                          &   IART_CHEM_NO,          &
                                          &   IART_SIMNOY_PRES,      &
                                          &   IART_SIMNOY_WMO,       &
                                          &   IART_SIMNOY_EXTP,      &
                                          &   IART_SIMNOY_SEDI
  USE mo_art_OH_chem_types,             ONLY: t_art_OH_chem_meta
  USE mo_art_OH_chem_predcorr_method,   ONLY: art_calculate_OH_number_concentration, &
                                          &   art_OH_chem_set_tdep_meta,             &
                                          &   art_predictor_corrector_method,        &
                                          &   luse_gems_ozone
  USE mo_art_chem_types,                ONLY: t_chem_meta_param
  USE mo_art_chem_types_param,          ONLY: t_chem_meta_lt, t_chem_meta_OH,        &
                                          &   t_chem_meta_linoz, t_chem_meta_simnoy, &
                                          &   t_chem_meta_cold, OH_get_destruct,     &
                                          &   linoz_integrate
  USE mo_art_chem_utils,                ONLY: art_get_CO2_des, art_get_CO_des,     &
                                          &   art_get_CH4_des, art_get_CH3CN_des,  &
                                          &   art_get_H2O_des, art_calc_col
  USE mo_art_linoz,                     ONLY: art_calc_linoz_polarchem_lt,   &
                                          &   art_calc_linoz_polarchem,      &
                                          &   art_calc_linoz, art_calc_linoz_ana
  USE mo_art_simnoy,                    ONLY: art_n2onoy_chemistry_wmo,      &
                                          &   art_n2onoy_chemistry_extp,     &
                                          &   art_n2onoy_chemistry_pres,     &
                                          &   art_noy_polarchem,             &
                                          &   art_noy_polarchem_sedi
  USE mo_art_cell_loop,                 ONLY: art_loop_cell_tracer,         &
                                          &   art_loop_cell_tracer_td,      &
                                          &   art_loop_cell_array
  USE mo_art_diagnostics,               ONLY: art_calc_o2_column, art_chem_calc_column
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c


  IMPLICIT NONE

  PRIVATE

  PUBLIC   ::                        &
      &   art_loss_chemtracer

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_chemtracer'

CONTAINS

!!
!!-----------------------------------------------------------------------------
!!

SUBROUTINE art_loss_chemtracer(jg,current_date,p_dtime,p_prog_list,p_tracer_now)
!<
! SUBROUTINE art_loss_chemtracer                 
! This subroutine calculates loss terms for the tracers using either implicit Euler method 
! or predictor-corrector method. Cold tracer, N2O-NOy
! chemistry and linoz tracer are treated separately.
! Part of Module: mo_art_chemtracer
! Author: Roland Ruhnke, KIT
! Initial Release: 2013-05-15                
! Modifications: 
! YYYY-MM-DD: Jennifer Schroeter, KIT
! - 
! 2015-03-04: Christian Stassen, KIT
! - included linoz (?)
! until 2016-04-13: Michael Weimer, KIT
! - included OH chemistry mechanism with predictor-corrector method
! around 2016-03-10: Roland Ruhnke, KIT
! - included artificial and regional tracers (??)
! 2016-04-15: Michael Weimer, KIT
! - added (simple) CO2 deposition into the ocean
!   based on Jacob (1999): Introduction to atmospheric chemistry
! 2016-11-22: Michael Weimer, KIT
! - included reading PSC climatologies
! 2017-03-15: Christopher Diekmann, KIT
! - added Noy/N2O mechanism
!>
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &
    &  jg                  !< patch id
  TYPE(datetime),POINTER,INTENT(IN) :: &
    &  current_date        !< current date and time
  REAL(wp), INTENT(IN)              ::  &
    &  p_dtime             !< time step
  TYPE(t_var_list_ptr), INTENT(INOUT)   ::  &
    &  p_prog_list         !< current prognostic state list

  REAL(wp), INTENT(INOUT), TARGET ::  &
    &  p_tracer_now(:,:,:,:)  !< tracer mixing ratios (specific concentrations)
                              !< at current time level n (before transport)
                              !< [kg/kg]
                              !< dim: (nproma,nlev,nblks_c,ntracer)
  
  ! local variables
  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element
  TYPE(t_art_OH_chem_meta), POINTER  ::  &
    &  OH_chem_meta     !< meta data for OH chemistry
  REAL(wp), POINTER ::   &
    &  ozone(:,:,:)         !< pointer to ozone for OH chemistry
  INTEGER, POINTER ::   &
    &  jsp               !< returns index of element
  INTEGER ::            &
    &  jb,              &!< loop indice
    &  i_startidx, i_endidx
  TYPE(t_art_chem_indices),POINTER    :: &
    &  art_indices          !< Pointer to ART chem indices
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem          !< Pointer to ART chem fields
  TYPE(t_art_chem_param),POINTER    :: &
    &  art_param          !< Pointer to ART paretrised chem fields
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo          !< Pointer to ART atmo fields
  INTEGER                     :: &
    &  iv                 !< loop index

  ! ----------------------------------
  ! --- Init pointer
  ! ----------------------------------

  art_chem => p_art_data(jg)%chem
  art_param => p_art_data(jg)%chem%param
  art_indices => p_art_data(jg)%chem%indices
  art_atmo => p_art_data(jg)%atmo
    

  ! ------------------------------------------------------------
  ! --- Start routine with some calculations of the OH chemistry
  ! ------------------------------------------------------------

  IF (art_param%OH_chem_meta%is_init) THEN
    OH_chem_meta => art_param%OH_chem_meta

    IF (luse_gems_ozone) THEN
      ozone => art_atmo%o3_clim
    ELSE
        IF (art_indices%iTRO3 /= 0) THEN
          ozone => p_tracer_now(:,:,:,art_indices%iTRO3)
        ELSE
          CALL finish('mo_art_chemtracer',                                 &
                &     'luse_gems_ozone is FALSE but O3 is not found in tracer list')
        ENDIF
    END IF

    CALL art_OH_chem_set_tdep_meta(jg,OH_chem_meta,art_atmo%temp, ozone, &
                                &  p_tracer_now(:,:,:,art_indices%iTRCH4))

    CALL art_calculate_OH_number_concentration(jg,OH_chem_meta,                                   &
           &  art_atmo%temp,art_chem%photo%rate(:,:,:,4), p_tracer_now(:,:,:,art_indices%iTRCH4), &
           &  p_tracer_now(:,:,:,art_indices%iTRCO))

  END IF

  ! ----------------------------------
  ! --- Start
  ! ----------------------------------


  DO iv = 1, p_prog_list%p%nvars

    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    IF (info_dyn%tracer%lis_tracer) THEN

      jsp=>info%ncontained

      SELECT TYPE(tracer => info_dyn%tracer)
        TYPE IS (t_chem_meta_param)
          CALL tracer%set_tracer(p_tracer_now(:,:,:,jsp))

        TYPE IS (t_chem_meta_lt)

          CALL tracer%set_tracer(p_tracer_now(:,:,:,jsp))

          IF (jsp == art_indices%iTRCH4) THEN
            CALL art_loop_cell_tracer(jg, tracer, art_get_CH4_des)
          ELSE IF (jsp == art_indices%iTRCO) THEN
            CALL art_loop_cell_tracer(jg, tracer, art_get_CO_des)
          ELSE IF (jsp == art_indices%iTRCO2) THEN
            CALL art_loop_cell_array(jg, art_param%o2_column, art_calc_o2_column)
            CALL art_loop_cell_tracer(jg, tracer, art_get_CO2_des)
          ELSE IF (jsp == art_indices%iTRCH3CN) THEN
            CALL art_loop_cell_tracer(jg, tracer, art_get_CH3CN_des)
          ELSE IF (jsp == art_indices%iTRH2O) THEN
            CALL art_loop_cell_tracer(jg, tracer, art_get_H2O_des)
          END IF

        TYPE IS (t_chem_meta_cold)

          CALL tracer%set_tracer(p_tracer_now(:,:,:,jsp))
          CALL tracer%get_destruct(jg)

        TYPE IS (t_chem_meta_linoz)

          IF (tracer%polarchem == IART_POLARCHEM) THEN
            CALL tracer%set_tracer_linoz(p_tracer_now(:,:,:,jsp),  &
                           &             p_tracer_now(:,:,:,art_indices%iTR_cold))
          ELSE
            CALL tracer%set_tracer(p_tracer_now(:,:,:,jsp))
          END IF

          CALL art_loop_cell_tracer(jg, tracer, art_calc_col)

          IF (tracer%O3_paramet == IART_LINOZ_ANA) THEN
            CALL art_loop_cell_tracer_td(jg, tracer,current_date,p_dtime, art_calc_linoz_ana)
          ELSE
            CALL art_loop_cell_tracer_td(jg, tracer,current_date,p_dtime, art_calc_linoz)
          END IF

          SELECT CASE(tracer%polarchem)
            CASE(IART_POLARCHEM)
              CALL art_loop_cell_tracer_td(jg, tracer, current_date,p_dtime,  &
                           &               art_calc_linoz_polarchem)
            CASE(IART_LINOZ_LT)
              CALL art_loop_cell_tracer_td(jg, tracer, current_date,p_dtime,   &
                          &                art_calc_linoz_polarchem_lt)
          END SELECT

        TYPE IS (t_chem_meta_simnoy)
          IF (jsp == art_indices%iTRNOy) THEN
            IF (tracer%polarchem /= IART_CHEM_NO) THEN
              CALL tracer%set_tracer_simnoy(p_tracer_now(:,:,:,jsp),                &
                               &            p_tracer_now(:,:,:,art_indices%iTRN2O),   &
                               &            p_tracer_now(:,:,:,art_indices%iTR_cold), &
                               &            p_prog_list, jg)
            ELSE
              CALL tracer%set_tracer_simnoy(p_tracer_now(:,:,:,jsp),  &
                               &            p_tracer_now(:,:,:,art_indices%iTRN2O))
            END IF

            SELECT CASE (tracer%cnoyn2o_tropo)
              CASE (IART_SIMNOY_EXTP)
                CALL art_loop_cell_tracer_td(jg, tracer, current_date,p_dtime,  &
                              &              art_n2onoy_chemistry_extp)
              CASE (IART_SIMNOY_WMO)
                CALL art_loop_cell_tracer_td(jg, tracer, current_date,p_dtime,  &
                              &              art_n2onoy_chemistry_wmo)
              CASE (IART_SIMNOY_PRES)
                CALL art_loop_cell_tracer_td(jg, tracer, current_date,p_dtime,  &
                              &              art_n2onoy_chemistry_pres)
            END SELECT


            SELECT CASE (tracer%polarchem)
              CASE (IART_POLARCHEM)
                CALL art_loop_cell_tracer(jg, tracer,  &
                              &           art_noy_polarchem)
              CASE(IART_SIMNOY_SEDI)
                CALL art_loop_cell_tracer_td(jg, tracer,current_date,p_dtime,  &
                              &              art_noy_polarchem_sedi)
            END SELECT
          END IF


        TYPE IS (t_chem_meta_OH)
          CALL tracer%set_tracer(p_tracer_now(:,:,:,jsp))
          CALL OH_get_destruct(tracer,jg,OH_chem_meta,.FALSE.)

      END SELECT

  ENDIF

  
  ENDDO


  ! Missing
  ! iTRH2O_feed


  DO iv = 1, p_prog_list%p%nvars
    SELECT TYPE(tracer =>  p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      CLASS IS (t_chem_meta_param)
        CALL tracer%get_prod()
    END SELECT
    
  END DO



  IF (p_art_data(jg)%chem%param%OH_chem_meta%is_init) THEN
    CALL art_predictor_corrector_method  &
      & (jg,OH_chem_meta,p_prog_list,  &
      &  p_dtime)
  END IF



  DO iv = 1, p_prog_list%p%nvars

    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    IF (info_dyn%tracer%lis_tracer) THEN

      jsp=>info%ncontained

      SELECT TYPE(tracer => info_dyn%tracer)
        TYPE IS (t_chem_meta_param)
          CALL tracer%add_prod(p_dtime)

        TYPE IS (t_chem_meta_lt)
          CALL tracer%integrate(p_dtime)

        TYPE IS (t_chem_meta_cold)
          CALL tracer%integr_and_prescribe(jg, p_dtime, p_art_data(jg)%chem%vmr2Nconc)

        TYPE IS (t_chem_meta_linoz)
        
          CALL linoz_integrate(tracer)

        TYPE IS (t_chem_meta_simnoy)
          IF (jsp == art_indices%iTRNOy) THEN
            CALL tracer%integrate()
          END IF

      END SELECT

    ENDIF 

    ! ----------------------------------
    ! --- Total column
    ! ----------------------------------

    IF (art_config(jg)%lart_diag_out) THEN
      IF (jsp == art_indices%iTRSO2) THEN
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
          CALL art_chem_calc_column(art_chem%so2_column(:,:,jb), art_atmo%pres(:,:,jb),               &
            &                       art_atmo%temp(:,:,jb), art_atmo%z_ifc(:,:,jb),                    &
            &                       p_tracer_now(:,:,jb, jsp), i_startidx, i_endidx, art_atmo%nlev)
          ! Transform Column values from mol/m2 into DU
          art_chem%so2_column(:,:,jb) = art_chem%so2_column(:,:,jb) / 0.44615e-3_wp
        ENDDO
      ENDIF
    ENDIF
  
  ENDDO !loop elements


!  ! ----------------------------------
!  ! --- TR_vortex (Vortexsplit tracer)
!  ! ----------------------------------
!
!  IF (art_chem%iTR_vortex /= 0) THEN
!    ! ----------------------------------
!    ! --- Should artificial vortex
!    ! --- tracer (TR_vortex) be set ?
!    ! ----------------------------------
!    ALLOCATE(ptr_vortex_log(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
!
!    CALL art_define_polar_vortex_boundary(current_date,ptr_vortex_log, select_hemisphere,  &
!                  &                       p_patch(jg), 65, 10 )
!    ! --- Set the tracer to zero on the desired half of the earth before setting new vortextracer
!    SELECT CASE ( select_hemisphere )
!      CASE ( 'SH' )
!        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
!          CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
!          DO jk = 10,65
!            DO jc = i_startidx,i_endidx
!              IF ( art_atmo%lat(jc,jb) <= 0._wp ) THEN
!                 p_tracer_now(jc,jk,jb,art_chem%iTR_vortex) = 0._wp
!
!                 IF ( ptr_vortex_log(jc,jk,jb) ) THEN
!                   p_tracer_now(jc,jk,jb,art_chem%iTR_vortex) = 1._wp
!                 END IF
!               END IF
!             END DO
!          END DO
!        END DO
!      CASE ( 'NH' )
!        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
!          CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
!          DO jk = 10,65
!            DO jc = i_startidx,i_endidx
!              IF ( art_atmo%lat(jc,jb) >= 0._wp ) THEN
!                p_tracer_now(jc,jk,jb,art_chem%iTR_vortex) = 0._wp
!
!                IF ( ptr_vortex_log(jc,jk,jb) ) THEN
!                  !RR p_tracer_now(jc,jk,jb,art_chem%iTR_vortex)   &
!                         &   = abs(p_patch%cells%center(jc,jb)%lat)*180._wp/pi
!
!                  p_tracer_now(jc,jk,jb,art_chem%iTR_vortex) = 1._wp
!                END IF
!              END IF
!            ENDDO
!          ENDDO
!        ENDDO
!    END SELECT !select_hemisphere
!  ENDIF


END SUBROUTINE art_loss_chemtracer

!!
!!-----------------------------------------------------------------------------
!!




END MODULE mo_art_chemtracer
