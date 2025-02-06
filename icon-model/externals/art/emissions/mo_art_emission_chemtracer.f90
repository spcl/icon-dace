!
! mo_art_emission_chemtracer
! This module provides the emission subroutines for chemical tracer processes.
! This Module is only envoked if lart_chemtracer = .TRUE.
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

MODULE mo_art_emission_chemtracer
    ! ICON
    USE mo_math_constants,              ONLY: pi
    USE mo_kind,                        ONLY: wp
    USE mo_physical_constants,          ONLY: amw, amd
    USE mo_var_list,                    ONLY: t_var_list_ptr
    USE mo_var,                         ONLY: t_var
    USE mo_var_metadata_types,          ONLY: t_var_metadata_dynamic, t_var_metadata

    USE mtime,                          ONLY: datetime,                     &
                                          &   timedelta,newTimedelta,       &
                                          &   getTimeDeltaFromDateTime,     &
                                          &   getTotalMillisecondsTimedelta,&
                                          &   deallocateTimedelta
    USE mo_time_config,                 ONLY: time_config
    ! ART
    USE mo_art_cell_loop,               ONLY: art_loop_cell_tracer_td,   &
                                          &   art_loop_cell_tracer,      &
                                          &   art_loop_cell_array
    USE mo_art_chem_types,              ONLY: t_chem_meta_passive
    USE mo_art_chem_types_param,        ONLY: t_chem_meta_lt
    USE mo_art_external_types,          ONLY: t_art_online_dms
    USE mo_art_emission_online_dms,     ONLY: art_add_dms_emission_to_tracer
    USE mo_art_data,                    ONLY: p_art_data
    USE mo_art_chem_data,               ONLY: t_art_chem_indices,  &
                                          &   t_art_chem_param
    USE mo_art_atmo_data,               ONLY: t_art_atmo
    USE mo_art_emission_regio_tracers,  ONLY: art_emiss_TR_art,          &
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


  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC   ::                        &
      &   art_emiss_chemtracer

  CHARACTER(LEN=*), PARAMETER :: routine = 'mo_art_emission_chemtracer'
         
CONTAINS

!!-----------------------------------------------------------------------------
!! art_emiss_chemtracer subroutine provides the boundary conditions 
!! of the chemical and transport tracer
!!-----------------------------------------------------------------------------

SUBROUTINE art_emiss_chemtracer(current_date,dtime,p_tracer_now,jg,p_prog_list)
!<
! SUBROUTINE art_emiss_chemtracer
! This subroutine sets the lower boundary condition for chemical tracers
! Passive tracers based on: Vogel (2015)
! Part of Module: mo_art_emission_chemtracer
! Author: Roland Ruhnke, KIT (et al. see above)
! Initial Release: 2013-05-15
! Modifications:
! 2015-08-06: Jennifer Schroeter, KIT
! -
! 2015-08-06: Michael Weimer, KIT
! - Included time-dependent prescribed emissions
! 2016-03-15: Roland Ruhnke, KIT
! - included passive tracers (?)
! 2016-07-27: Jennifer Schroeter, KIT
! - Envoked lpassive switch
! 2016-10-28: Michael Weimer, KIT
! - removed prescribed emission (moved to mo_art_emission_interface)
!>
  IMPLICIT NONE

  INTEGER, INTENT(IN)               :: jg              !< patch on which computation is performed

  REAL(wp), INTENT(INOUT), TARGET ::  &  !< tracer mixing ratios (specific concentrations)
    &  p_tracer_now(:,:,:,:)             !< at current time level n (before transport)
                                         !< [kg/kg]
                                         !< dim: (nproma,nlev,nblks_c,ntracer)
  REAL(wp), INTENT(IN) :: dtime  ! time step    

  TYPE(datetime), POINTER, INTENT(IN)  :: current_date
  TYPE(t_var_list_ptr), INTENT(INOUT)   ::  &
    &  p_prog_list         !< current prognostic state list

  !local variables
  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element
  INTEGER, POINTER ::   &
    &  jsp               !< returns index of element
  INTEGER          ::   &
    &  iv

  TYPE(t_art_online_dms), POINTER :: &
    &  p_online_dms           !< online dms
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                !< Pointer to ART atmo fields
  TYPE(t_art_chem_indices),POINTER    :: &
    &  art_indices             !< Pointer to ART chem indices
  TYPE(t_art_chem_param),POINTER    :: &
    &  art_param               !< Pointer to ART parametrised chem fields

  !!------------------------------------------------------------------------------
  !! --- Allocation
  !!------------------------------------------------------------------------------


  art_indices    => p_art_data(jg)%chem%indices
  art_param      => p_art_data(jg)%chem%param
  art_atmo       => p_art_data(jg)%atmo
  p_online_dms   => p_art_data(jg)%ext%online_dms
  


  DO iv = 1, p_prog_list%p%nvars

    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    IF (info_dyn%tracer%lis_tracer) THEN

      jsp=>info%ncontained

      SELECT TYPE(tracer => info_dyn%tracer)
        CLASS IS (t_chem_meta_lt)

          IF (art_indices%iTRCHBr3 == jsp )   THEN
            CALL tracer%set_tracer(p_tracer_now(:,:,:,art_indices%iTRCHBr3))
            CALL art_loop_cell_tracer(jg, tracer, art_emiss_CHBr3)
          ELSE IF (art_indices%iTRCH2Br2 == jsp )  THEN
            CALL tracer%set_tracer(p_tracer_now(:,:,:,art_indices%iTRCH2Br2))
            CALL art_loop_cell_tracer(jg, tracer, art_emiss_CH2Br2)
          ELSE IF (art_indices%iTRH2O == jsp ) THEN 
            CALL tracer%set_tracer(p_tracer_now(:,:,:,art_indices%iTRH2O))
            CALL art_loop_cell_tracer(jg, tracer, art_emiss_H2O)
          ELSE IF (art_indices%iTRDMS == jsp .AND. p_online_dms%lcalc_onl) THEN
            !!----------------------
            !! --- Online DMS Calculation
            !!----------------------
            CALL tracer%set_tracer(p_tracer_now(:,:,:,art_indices%iTRDMS))
            CALL art_loop_cell_tracer_td(jg, tracer, current_date, dtime,   &
                         &               art_add_dms_emission_to_tracer)
          ENDIF

        CLASS IS (t_chem_meta_passive)

           !-------------------------------------------------------------------------------
           ! --- Prescribing artificial age of air tracer near the surface (TRAGE)
           !-------------------------------------------------------------------------------
           IF (art_indices%iTRAGE == jsp ) THEN
             CALL tracer%set_tracer(p_tracer_now(:,:,:,art_indices%iTRAGE))
             CALL art_loop_cell_tracer_td(jg, tracer, current_date, dtime, art_emiss_AGE)
           END IF

           ! only if regional tracers are being calculated, this subroutine is called
           IF (art_param%lregio_tracers) THEN
             IF (jsp == art_indices%iTR_art) THEN
               CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_art),   &
                             &          art_emiss_TR_art)
             ENDIF

             IF (current_date%time%hour == 0) THEN
               IF (jsp == art_indices%iTR_ech) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_ech),   &
                               &          art_emiss_TR_ech)
               ELSE IF (jsp == art_indices%iTR_eur) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_eur),   &
                               &          art_emiss_TR_eur)
               ELSE IF (jsp == art_indices%iTR_naf) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_naf),   &
                               &          art_emiss_TR_naf)
               ELSE IF (jsp == art_indices%iTR_saf) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_saf),   &
                               &          art_emiss_TR_saf)
               ELSE IF (jsp == art_indices%iTR_sea) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_sea),   &
                               &          art_emiss_TR_sea)
               ELSE IF (jsp == art_indices%iTR_mdg) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_mdg),   &
                               &          art_emiss_TR_mdg)
               ELSE IF (jsp == art_indices%iTR_sam) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_sam),   &
                               &          art_emiss_TR_sam)
               ELSE IF (jsp == art_indices%iTR_nam) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_nam),   &
                               &          art_emiss_TR_nam)
               ELSE IF (jsp == art_indices%iTR_aus) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_aus),   &
                               &          art_emiss_TR_aus)
               ELSE IF (jsp == art_indices%iTR_med) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_med),   &
                               &          art_emiss_TR_med)
               ELSE IF (jsp == art_indices%iTR_sib) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_sib),   &
                               &          art_emiss_TR_sib)
               ELSE IF (jsp == art_indices%iTR_sin) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_sin),   &
                               &          art_emiss_TR_sin)
               ELSE IF (jsp == art_indices%iTR_nin) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_nin),   &
                               &          art_emiss_TR_nin)
               ELSE IF (jsp == art_indices%iTR_tpo) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_tpo),   &
                               &          art_emiss_TR_tpo)
               ELSE IF (jsp == art_indices%iTR_tio) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_tio),   &
                               &          art_emiss_TR_tio)
               ELSE IF (jsp == art_indices%iTR_tao) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_tao),   &
                               &          art_emiss_TR_tao)
               ELSE IF (jsp == art_indices%iTR_bgs) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_bgs),   &
                               &          art_emiss_TR_bgs)
               ELSE IF (jsp == art_indices%iTR_bgn) THEN
                 CALL art_loop_cell_array(jg, p_tracer_now(:,:,:,art_indices%iTR_bgn),   &
                               &          art_emiss_TR_bgn)
               ENDIF
             END IF ! hour == 0

           END IF  ! art_param%lregio_tracers
         
         
      END SELECT
    END IF

  END DO

  NULLIFY(art_atmo)

END SUBROUTINE art_emiss_chemtracer
!
! ----------------------------------------------------------
!
SUBROUTINE art_emiss_AGE(jg,jb,jcs,jce ,&
            &             tracer,current_date, p_dtime)
!<
! SUBROUTINE art_emiss_AGE
! This subroutine sets the lower boundary condition for TRAGE
! Part of Module: mo_art_emission_chemtracer
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!>
   INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
   TYPE(t_chem_meta_passive), INTENT(INOUT) :: tracer
   TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
   REAL(wp), INTENT(IN) :: p_dtime
         

  ! local variables
  INTEGER :: &
    &  i ,jk
  TYPE(timedelta),POINTER :: &
    &  time_diff               !< mtime object: Elapsed time since experiment start
  REAL(wp)                :: &
    &  p_sim_time              !< elapsed simulation time on this grid level
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo
 

  time_diff  => newTimedelta("PT0S")
  time_diff  =  getTimeDeltaFromDateTime(current_date, time_config%tc_exp_startdate)
  p_sim_time =  getTotalMillisecondsTimedelta(time_diff, current_date)*1.e-3_wp
  CALL deallocateTimedelta(time_diff)

  DO jk = 1, art_atmo%nlev
    DO i = jcs, jce
      IF (art_atmo%pres(i,jk,jb) >= 950._wp*100._wp) THEN
        tracer%tracer(i,jk,jb) = 7._wp+(p_sim_time/86400._wp)/365.2425_wp
        ! sim_time : [sec] seconds elapsed since the model start
        ! linear increase of TRAGE in troposphere +1 per year; Offset = +7
      END IF
    END DO
  END DO
END SUBROUTINE art_emiss_AGE
!
! -----------------------------------------------------------
!
SUBROUTINE art_emiss_CHBr3(jg,jb,jcs,jce ,&
            &             nproma,nlev,tracer)
!<
! SUBROUTINE art_emiss_AGE
! This subroutine sets the lower boundary condition for CHBr3
! Part of Module: mo_art_emission_chemtracer
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!>
  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  TYPE(t_chem_meta_lt), INTENT(INOUT) :: tracer

  ! local variables
  INTEGER :: &
    &  i ,jk
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  DO i = jcs, jce
    ! only water points in lat <= abs (15 degrees) 
    IF (.NOT. art_atmo%llsm(i,jb) .AND. ABS(art_atmo%lat(i,jb))*180._wp/pi <= 15._wp) THEN

      DO jk = 1,art_atmo%nlev
        IF (art_atmo%pres(i,jk,jb) >= 900._wp*100._wp) THEN
          tracer%tracer(i,jk,jb)  = 1.60e-12_wp * tracer%mol_weight / amd * 1000._wp
        END IF
      END DO

    END IF
  END DO
END SUBROUTINE art_emiss_CHBr3
!
! -----------------------------------------------------------
!
SUBROUTINE art_emiss_CH2Br2(jg,jb,jcs,jce ,&
            &             nproma,nlev,tracer)
!<
! SUBROUTINE art_emiss_AGE
! This subroutine sets the lower boundary condition for CH2Br2
! Part of Module: mo_art_emission_chemtracer
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!>
  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  TYPE(t_chem_meta_lt), INTENT(INOUT) :: tracer

  ! local variables
  INTEGER :: &
    &  i ,jk
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  DO i = jcs, jce
    ! only water points in lat <= abs (15 degrees) 
    IF (.NOT. art_atmo%llsm(i,jb) .AND. ABS(art_atmo%lat(i,jb))*180._wp/pi <= 15._wp) THEN

      DO jk = 1,art_atmo%nlev
        IF (art_atmo%pres(i,jk,jb) >= 900._wp*100._wp) THEN
          tracer%tracer(i,jk,jb)  = 1.10e-12_wp * tracer%mol_weight / amd * 1000._wp
        END IF
      END DO

    END IF
  END DO
END SUBROUTINE art_emiss_CH2Br2
!
! -----------------------------------------------------------
!
SUBROUTINE art_emiss_H2O(jg,jb,jcs,jce ,&
            &             nproma,nlev,tracer)
!<
! SUBROUTINE art_emiss_AGE
! This subroutine sets the lower boundary condition for CH2Br2
! Part of Module: mo_art_emission_chemtracer
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-29
! Modifications:
!>
  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  INTEGER        ,INTENT(IN) :: nproma,nlev
  TYPE(t_chem_meta_lt), INTENT(INOUT) :: tracer

  ! local variables
  INTEGER :: &
    &  i ,jk
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  DO jk = 1,art_atmo%nlev
    DO i = jcs, jce
      IF (art_atmo%pres(i,jk,jb) >= 300._wp*100._wp) THEN
        tracer%tracer(i,jk,jb)   = 1.00e-04_wp * amw / amd
      END IF
    END DO
  END DO
END SUBROUTINE art_emiss_H2O
!
! -----------------------------------------------------------
!
END MODULE mo_art_emission_chemtracer
