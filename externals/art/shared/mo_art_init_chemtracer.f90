!
! mo_art_init_chemtracer
! This module provides initialization routines for chemistry
! tracer with fixed values in ICON-ART
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

MODULE mo_art_init_chemtracer
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_math_constants,                ONLY: rad2deg
  USE mo_exception,                     ONLY: finish
  USE mo_master_config,                 ONLY: isRestart
  USE mo_run_config,                    ONLY: ntracer
  USE mo_art_config,                    ONLY: art_config
  USE mo_io_units,                      ONLY: filename_max
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_tracer_metadata_types,         ONLY: t_chem_meta
  USE mo_impl_constants,                ONLY: SUCCESS
! ART
  USE mo_art_vinterp ,                  ONLY: art_prepare_vinterp
  USE mo_art_vinterp_chem_init,         ONLY: art_vinterp_chem_init,art_read_in_chem_init 
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_data,                 ONLY: t_art_chem_indices
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_chem_init_types,           ONLY: t_chem_init_state
  USE mo_physical_constants,            ONLY: ppmv2gg
  USE mo_util_string,                   ONLY: int2string
  USE mo_exception,                     ONLY: message
  USE mo_art_chem_init_utils,           ONLY: deallocate_chem_init_chem, &
                                          &   deallocate_chem_init_atm
  
                                         
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_init_chemtracer, art_init_chemtracer_ext, art_get_chemtracer_idx
  PUBLIC :: art_set_init_tracer_trop_strat

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_init_chemtracer'

CONTAINS

SUBROUTINE art_init_chemtracer(jg, tracer)
!<
! SUBROUTINE art_init_chemtracer
! This subroutine sets initial values for tracers
! Part of Module: mo_art_init_chemtracer
! Author: Jennifer Schroeter, KIT
! Initial Release: YYYY-MM-DD
! Modifications:
! 2016-11-08: Michael Weimer, KIT
! - removed emission initialisation
!>
  INTEGER, INTENT(in)               :: &
    &  jg                                 !< patch on which computation is performed
  REAL(wp), INTENT(inout)           :: &  !< tracer mixing ratios (specific concentrations)
    &  tracer(:,:,:,:)                    !< at current time level n (before transport)
  ! local variables
  TYPE(t_art_atmo),POINTER    :: &
      &  art_atmo                     !< Pointer to ART atmo fields
  TYPE(t_art_chem_indices),POINTER    :: &
      &  art_indices                     !< Pointer to ART indicesetrised chem fields
  INTEGER ::                  &
    &  nlev,                  &
    &  jk,                    &           !< loop index
    &  jb, jc, i_startidx, i_endidx
  REAL(wp) :: zlat        !< latitude (deg)


  art_atmo => p_art_data(jg)%atmo
  art_indices => p_art_data(jg)%chem%indices
 
  IF (.NOT. isRestart()) THEN

    ! ----------------------------------
    ! CO initialisation
    ! ----------------------------------
    ! Values are valid for 90 level version of ICON 
    
    IF (art_indices%iTRCO   /=0) THEN 
      CALL art_set_init_tracer_trop_strat(jg,                               &
                      &                   tracer(:,:,:,art_indices%iTRCO),  &
                      &                   1.09e-07_wp,                      &
                      &                   1.8e-08_wp,                       &
                      &                   5)
    ENDIF 
    ! --------------------------------------------
    ! Initialisation of propane: tropospheric value estimated from
    ! Pozzer (2010), ACP, 10, 4403-4422
    ! --------------------------------------------
    IF (art_indices%iTRC3H8  /= 0) THEN
      CALL art_set_init_tracer_trop_strat(jg,                                &
                      &                   tracer(:,:,:,art_indices%iTRC3H8), &
                      &                   1.0e-10_wp,                        &
                      &                   0.0_wp,                            &
                      &                   5)
    ENDIF

    !--------------------------------
    ! OCS Initialisation: estimation according to grafic in Turco (1979)
    !--------------------------------
    IF (art_indices%iTROCS /=0)THEN
      CALL art_set_init_tracer_trop_strat(jg,                               &
                      &                   tracer(:,:,:,art_indices%iTROCS), &
                      &                   4.0e-10_wp,                       &
                      &                   2.14e-13_wp,                      &
                      &                   5)
    ENDIF

    !-------------------------------
    ! Acetone initialisation: According to Sprung and Zahn, JGR (2010)
    !------------------------------
    IF (art_indices%iTRCH3COCH3 /= 0) THEN
      CALL art_set_init_tracer_trop_strat(jg,                                    &
                      &                   tracer(:,:,:,art_indices%iTRCH3COCH3), &
                      &                   4.e-10_wp,                             &
                      &                   7.e-11_wp,                             &
                      &                   5)
    ENDIF
      
    ! -------------------------------
    ! Initialisation of other species, just first approx.
    ! -------------------------------

    IF (art_indices%iTRCHBR3        /= 0)     tracer(:,:,:,art_indices%iTRCHBR3)      = 1.00e-20_wp
    IF (art_indices%iTRCH2BR2       /= 0)     tracer(:,:,:,art_indices%iTRCH2BR2)     = 1.00e-20_wp
    IF (art_indices%iTRC2H6         /= 0)     tracer(:,:,:,art_indices%iTRC2H6)       = 1.00e-20_wp
    IF (art_indices%iTRC5H8         /= 0)     tracer(:,:,:,art_indices%iTRC5H8)       = 1.00e-20_wp
    IF (art_indices%iTRCH3CN        /= 0)     tracer(:,:,:,art_indices%iTRCH3CN)      = 1.00e-12_wp
    IF (art_indices%iTRCO2          /= 0)     tracer(:,:,:,art_indices%iTRCO2)        = 3.60e-04_wp 
    IF (art_indices%iTRH2O          /= 0)     tracer(:,:,:,art_indices%iTRH2O)        = 3.50e-06_wp
    IF (art_indices%iTRN2O          /= 0)     tracer(:,:,:,art_indices%iTRN2O)        = 3.00e-07_wp
    IF (art_indices%iTRNOy          /= 0)     tracer(:,:,:,art_indices%iTRNOy)        = 2.50e-10_wp
    IF (art_indices%iTRNO2          /= 0)     tracer(:,:,:,art_indices%iTRNO2)        = 1.00e-08_wp
    IF (art_indices%iTRN2O5         /= 0)     tracer(:,:,:,art_indices%iTRN2O5)       = 5.03e-10_wp
    IF (art_indices%iTRHNO3         /= 0)     tracer(:,:,:,art_indices%iTRHNO3)       = 1.00e-09_wp
    IF (art_indices%iTRNH3          /= 0)     tracer(:,:,:,art_indices%iTRNH3)        = 1.00e-11_wp
    IF (art_indices%iTRSO2          /= 0)     tracer(:,:,:,art_indices%iTRSO2)        = 1.00e-25_wp !1.00e-11_wp
    IF (art_indices%iTRDMS          /= 0)     tracer(:,:,:,art_indices%iTRDMS)        = 3.00e-10_wp
    IF (art_indices%iTRH2SO4        /= 0)     tracer(:,:,:,art_indices%iTRH2SO4)      = 1.00e-25_wp !1.00e-11_wp
    IF (art_indices%iTRCFCl3        /= 0)     tracer(:,:,:,art_indices%iTRCFCl3)      = 2.80e-10_wp
    IF (art_indices%iTRCF2Cl2       /= 0)     tracer(:,:,:,art_indices%iTRCF2Cl2)     = 5.03e-10_wp
    IF (art_indices%iTRHCl          /= 0)     tracer(:,:,:,art_indices%iTRHCl)        = 150.e-12_wp
    IF (art_indices%iTRHOCl         /= 0)     tracer(:,:,:,art_indices%iTRHOCl)       = 5.03e-10_wp
    IF (art_indices%iTRClONO2       /= 0)     tracer(:,:,:,art_indices%iTRClONO2)     = 150.e-12_wp
    IF (art_indices%iTRHBr          /= 0)     tracer(:,:,:,art_indices%iTRHBr)        = 5.03e-10_wp
    IF (art_indices%iTRHOBr         /= 0)     tracer(:,:,:,art_indices%iTRHOBr)       = 5.03e-10_wp
    IF (art_indices%iTRBrONO2       /= 0)     tracer(:,:,:,art_indices%iTRBrONO2)     = 5.03e-10_wp
    IF (art_indices%iTRAGE          /= 0)     tracer(:,:,:,art_indices%iTRAGE)        = 0.00_wp

    ! ----------------------------------
    ! --- TRAGE (Age of Air tracer)
    ! ----------------------------------

    IF (art_indices%iTRAGE /= 0) THEN
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg,jb,i_startidx,i_endidx)

        DO jk = 1,art_atmo%nlev
          DO jc = i_startidx,i_endidx
            tracer(jc,jk,jb,art_indices%iTRAGE) = 7._wp - (MIN( 7._wp,MAX( 0._wp, ( (-7000._wp *  &
                            &        LOG(art_atmo%pres(jc,jk,jb)/101325._wp))/1.e3_wp - 15._wp )  &
                            &        / 25._wp * 7._wp )))
          ENDDO
        ENDDO
      END DO
    END IF

    nlev = art_atmo%nlev

    IF (art_indices%iTR_stn    /= 0)   tracer(:,:,:,art_indices%iTR_stn)    = 0.00_wp
    IF (art_indices%iTR_stt    /= 0)   tracer(:,:,:,art_indices%iTR_stt)    = 0.00_wp
    IF (art_indices%iTR_sts    /= 0)   tracer(:,:,:,art_indices%iTR_sts)    = 0.00_wp
    IF (art_indices%iTR_trn    /= 0)   tracer(:,:,:,art_indices%iTR_trn)    = 0.00_wp
    IF (art_indices%iTR_trt    /= 0)   tracer(:,:,:,art_indices%iTR_trt)    = 0.00_wp
    IF (art_indices%iTR_trs    /= 0)   tracer(:,:,:,art_indices%iTR_trs)    = 0.00_wp
    IF (art_indices%iTR_tiln   /= 0)   tracer(:,:,:,art_indices%iTR_tiln)   = 0.00_wp
    IF (art_indices%iTR_tils   /= 0)   tracer(:,:,:,art_indices%iTR_tils)   = 0.00_wp
    IF (art_indices%iTR_nin    /= 0)   tracer(:,:,:,art_indices%iTR_nin)    = 0.00_wp
    IF (art_indices%iTR_sin    /= 0)   tracer(:,:,:,art_indices%iTR_sin)    = 0.00_wp
    IF (art_indices%iTR_ech    /= 0)   tracer(:,:,:,art_indices%iTR_ech)    = 0.00_wp
    IF (art_indices%iTR_sea    /= 0)   tracer(:,:,:,art_indices%iTR_sea)    = 0.00_wp
    IF (art_indices%iTR_sib    /= 0)   tracer(:,:,:,art_indices%iTR_sib)    = 0.00_wp
    IF (art_indices%iTR_eur    /= 0)   tracer(:,:,:,art_indices%iTR_eur)    = 0.00_wp
    IF (art_indices%iTR_med    /= 0)   tracer(:,:,:,art_indices%iTR_med)    = 0.00_wp
    IF (art_indices%iTR_naf    /= 0)   tracer(:,:,:,art_indices%iTR_naf)    = 0.00_wp
    IF (art_indices%iTR_saf    /= 0)   tracer(:,:,:,art_indices%iTR_saf)    = 0.00_wp
    IF (art_indices%iTR_mdg    /= 0)   tracer(:,:,:,art_indices%iTR_mdg)    = 0.00_wp
    IF (art_indices%iTR_aus    /= 0)   tracer(:,:,:,art_indices%iTR_aus)    = 0.00_wp
    IF (art_indices%iTR_nam    /= 0)   tracer(:,:,:,art_indices%iTR_nam)    = 0.00_wp
    IF (art_indices%iTR_sam    /= 0)   tracer(:,:,:,art_indices%iTR_sam)    = 0.00_wp
    IF (art_indices%iTR_tpo    /= 0)   tracer(:,:,:,art_indices%iTR_tpo)    = 0.00_wp
    IF (art_indices%iTR_tao    /= 0)   tracer(:,:,:,art_indices%iTR_tao)    = 0.00_wp
    IF (art_indices%iTR_tio    /= 0)   tracer(:,:,:,art_indices%iTR_tio)    = 0.00_wp
    IF (art_indices%iTR_bgn    /= 0)   tracer(:,:,:,art_indices%iTR_bgn)    = 0.00_wp
    IF (art_indices%iTR_bgs    /= 0)   tracer(:,:,:,art_indices%iTR_bgs)    = 0.00_wp
    IF (art_indices%iTR_art    /= 0)   tracer(:,:,:,art_indices%iTR_art)    = 0.00_wp
    IF (art_indices%iTR_nest   /= 0)   tracer(:,:,:,art_indices%iTR_nest)   = 0.00_wp

    DO jb = art_atmo%i_startblk, art_atmo%i_endblk
      CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

      IF (art_indices%iTR_nest /= 0) THEN
        IF (jg == 2) THEN
          DO jc = i_startidx, i_endidx
            DO jk=1,nlev
              tracer(jc,jk,jb,art_indices%iTR_nest) = 1._wp
            END DO
          END DO
        ELSE
          DO jc = i_startidx, i_endidx
            DO jk=1,nlev
              tracer(jc,jk,jb,art_indices%iTR_nest) = 0._wp
            END DO
          END DO
        END IF
      ENDIF

      ! ----------------------------------
      ! --- Set artificial tracer 
      ! --- in stratosphere, TTL and troposphere 
      ! ----------------------------------
      DO jc=i_startidx,i_endidx
        zlat = art_atmo%lat(jc,jb)*rad2deg
        IF (zlat <= -23._wp) THEN
          DO jk=1,art_atmo%ktrpwmo(jc,jb)
            IF (art_atmo%theta(jc,jk,jb) >= 380._wp) THEN
                IF ( art_indices%iTR_sts /= 0)    &
                    &  tracer(jc,jk,jb,art_indices%iTR_sts)  = 1._wp
                IF ( art_indices%iTR_tils /= 0)   &
                    &  tracer(jc,jk,jb,art_indices%iTR_tils) = 0._wp
                IF ( art_indices%iTR_trs /= 0)    &
                    &  tracer(jc,jk,jb,art_indices%iTR_trs)  = 0._wp
            ELSE
              IF ( art_indices%iTR_sts  /= 0)     &
                  &  tracer(jc,jk,jb,art_indices%iTR_sts)  = 0._wp
              IF ( art_indices%iTR_tils /= 0)     &
                  &  tracer(jc,jk,jb,art_indices%iTR_tils) = 1._wp
              IF ( art_indices%iTR_trs  /= 0)     &
                  &  tracer(jc,jk,jb,art_indices%iTR_trs)  = 0._wp
            ENDIF
          ENDDO
          IF ( art_indices%iTR_sts  /= 0)   &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_sts)  = 0._wp
          IF ( art_indices%iTR_tils /= 0)   &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_tils) = 0._wp
          IF ( art_indices%iTR_trs  /= 0)   &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_trs)  = 1._wp
        ELSE IF (zlat >= -23._wp .and. &
          &      zlat <=  23._wp ) THEN
          !Strat-Tracer
          IF ( art_indices%iTR_stt  /= 0)   &
              &  tracer(jc,:art_atmo%ktrpwmo(jc,jb),jb,art_indices%iTR_stt) = 1._wp
          IF ( art_indices%iTR_trt  /= 0)   &
              &  tracer(jc,:art_atmo%ktrpwmo(jc,jb),jb,art_indices%iTR_trt) = 0._wp
          !Trop-Tracer
          IF ( art_indices%iTR_stt  /= 0)   &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_stt) = 0._wp
          IF ( art_indices%iTR_trt  /= 0)   &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_trt) = 1._wp
        ELSE IF (zlat >= 23._wp) THEN
          DO jk=1,art_atmo%ktrpwmo(jc,jb)
            IF (art_atmo%theta(jc,jk,jb) >= 380._wp) THEN
              IF ( art_indices%iTR_stn  /= 0)   &
                  &  tracer(jc,jk,jb,art_indices%iTR_stn)  = 1._wp
              IF ( art_indices%iTR_tiln /= 0)   &
                  &  tracer(jc,jk,jb,art_indices%iTR_tiln) = 0._wp
              IF ( art_indices%iTR_trn  /= 0)   &
                  &  tracer(jc,jk,jb,art_indices%iTR_trn)  = 0._wp
            ELSE
              IF ( art_indices%iTR_stn  /= 0)   &
                  &  tracer(jc,jk,jb,art_indices%iTR_stn)  = 0._wp
              IF ( art_indices%iTR_tiln /= 0)   &
                  &  tracer(jc,jk,jb,art_indices%iTR_tiln) = 1._wp
              IF ( art_indices%iTR_trn  /= 0)   &
                  &  tracer(jc,jk,jb,art_indices%iTR_trn)  = 0._wp
            ENDIF
          ENDDO
          IF ( art_indices%iTR_stn  /= 0)    &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_stn)  = 0._wp
          IF ( art_indices%iTR_tiln /= 0)    &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_tiln) = 0._wp
          IF ( art_indices%iTR_trn  /= 0)    &
              &  tracer(jc,art_atmo%ktrpwmo(jc,jb):,jb,art_indices%iTR_trn)  = 1._wp

        ENDIF
      ENDDO
      ! ----------------------------------
      ! --- Set artIFicial tracer 
      ! --- in northern hemisphere (TR_nh) and southern hemisphere (TR_sh)
      ! ----------------------------------

      DO jc = i_startidx, i_endidx
        zlat = art_atmo%lat(jc,jb)*rad2deg
        DO jk=1,nlev
          IF (zlat <= 0._wp) THEN
            IF ( art_indices%iTR_sh  /= 0)   &
                &  tracer(jc,jk,jb,art_indices%iTR_sh) = 1._wp
            IF ( art_indices%iTR_nh  /= 0)   &
                &  tracer(jc,jk,jb,art_indices%iTR_nh) = 0._wp
          ELSE
            IF ( art_indices%iTR_sh  /= 0)   &
                &  tracer(jc,jk,jb,art_indices%iTR_sh) = 0._wp
            IF ( art_indices%iTR_nh  /= 0)   &
                &  tracer(jc,jk,jb,art_indices%iTR_nh) = 1._wp

          ENDIF
        ENDDO
      ENDDO

    END DO ! jb
  END IF
END SUBROUTINE art_init_chemtracer

!
!----------------------------------------------------------------------------------
!

SUBROUTINE art_set_init_tracer_trop_strat(jg,tracer,tracer_trop,tracer_strat,levels_trans)
!<
! SUBROUTINE art_set_init_tracer_trop_strat
! Sets tropospheric and stratospheric values and interpolates above the
! tropopause linearly in levels_trans levels.
! Part of Module: mo_art_init_chemtracer
! Author: Michael Weimer, KIT
! Initial Release: 2020-03-05
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg              !< patch id
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:)   !< tracer volume mixing ratio (mol / mol)
  REAL(wp), INTENT(in) :: &
    &  tracer_trop,       &  !< tropospheric value
    &  tracer_strat          !< stratospheric value
  INTEGER, INTENT(in)  :: &
    &  levels_trans          !< number of levels of the transistion
  !local variables
  INTEGER ::         &
    &  jb, jc, jk,   &  !< loop indices
    &  level_strat,  &  !< level from which stratospheric value is used
    &  level_trop,   &  !< level from which tropospheric value is used
    &  i_startidx,   &  !< loop indices
    &  i_endidx
  REAL(wp)           &
    &  diff_tp          !< difference of level strat and trop
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo         !< Pointer to ART atmo fields

  art_atmo => p_art_data(jg)%atmo

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg,jb,i_startidx,i_endidx)

    DO jc = i_startidx, i_endidx
      level_strat = MAX(art_atmo%ktrpwmo(jc,jb) - levels_trans, 1)
      level_trop  = art_atmo%ktrpwmop1(jc,jb)

      ! Stratosphere
      DO jk = 1, level_strat
        tracer(jc,jk,jb) = tracer_strat
      END DO

      ! Transition
      diff_tp = level_strat - level_trop

      DO jk = level_strat + 1, level_trop - 1
        tracer(jc,jk,jb) = tracer_trop + (jk - level_trop) / diff_tp * (tracer_strat - tracer_trop)
      END DO

      ! Tropopshere
      DO jk = level_trop, art_atmo%nlev
        tracer(jc,jk,jb) = tracer_trop
      END DO
    END DO
  END DO
END SUBROUTINE art_set_init_tracer_trop_strat

!
!-----------------------------------------------------------------------------------
!


SUBROUTINE art_init_chemtracer_ext(p_patch, p_prog_list,tracer)
!<
! SUBROUTINE art_init_chemtracer_ext
! <Description>
! Part of Module: mo_art_init_chemtracer
! Author: Jennifer Schroeter, KIT
! Initial Release: YYYY-MM-DD
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  TYPE(t_patch), TARGET, INTENT(in) ::  &
    &  p_patch  !< patch on which computation is performed

  TYPE(t_var_list_ptr), INTENT(in)       :: &
    &  p_prog_list      !< current prognostic state list
  REAL(wp), INTENT(inout)           :: &  !< tracer mixing ratios (specific concentrations)
    &  tracer(:,:,:,:)                    !< at current time level n (before transport)

  TYPE(t_var_metadata), POINTER      :: & 
    &  info                              !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata_dynamic), POINTER :: & 
    &  info_dyn                          !< returns reference to tracer metadata of current element
  INTEGER, POINTER ::  jsp                           !< returns index of element
  INTEGER          :: init_mode, iv

  CHARACTER(LEN=filename_max) :: coord_path
  CHARACTER(LEN=filename_max) :: chem_init_path
  INTEGER :: jg, ierror

  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                         !< Pointer to ART atmo fields
  TYPE(t_art_chem_indices),POINTER    :: &
    &  art_indices                         !< Pointer to ART indicesetrised chem fields
  TYPE(t_chem_init_state),POINTER    :: &
    &  chem_init                        !< Pointer to ART diagnostic fields

  REAL(wp):: conversion_factor        !< factor to convert e.g. kg/kg -> mol/mol

  INTEGER, TARGET :: &
    &  ntracer_target

  jg = p_patch%id

  chem_init    =>  p_art_data(jg)%chem%init
  art_atmo     =>  p_art_data(jg)%atmo
  art_indices  =>  p_art_data(jg)%chem%indices
  
  coord_path     = TRIM(art_config(jg)%cart_cheminit_coord)
  chem_init_path = TRIM(art_config(jg)%cart_cheminit_file)

  chem_init%chem_init_in%model_name = TRIM(art_config(jg)%cart_cheminit_type)
    
  conversion_factor = 1.0_wp
 
  DO iv = 1, p_prog_list%p%nvars

    info_dyn => p_prog_list%p%vl(iv)%p%info_dyn
    info => p_prog_list%p%vl(iv)%p%info
    IF (info_dyn%tracer%lis_tracer) THEN

      jsp => info%ncontained
      SELECT TYPE(meta => info_dyn%tracer)
        CLASS IS (t_chem_meta)
          CALL meta%opt_meta%get('init_mode',init_mode,ierror)

          IF (ierror/= SUCCESS) CALL finish(TRIM(routine)//':art_init_chemtracer_xml', &
               &         'Metadata init_mode not available for tracer '//TRIM(meta%name)//'.')

          IF (init_mode == 1 ) THEN

            IF ( .NOT. chem_init%chem_init_in%linitialized ) THEN
                ntracer_target = ntracer
                chem_init%chem_init_chem%n_spec_readin => ntracer_target

                CALL art_prepare_vinterp(chem_init,p_patch, coord_path, chem_init_path)
                CALL art_read_in_chem_init(chem_init,p_patch,p_prog_list, chem_init_path)
                CALL art_vinterp_chem_init(chem_init,ntracer,jg, p_prog_list)
            ENDIF
            
            IF (TRIM(ADJUSTL(art_config(jg)%cart_cheminit_type)) == 'ERA') THEN

                IF (jsp == art_indices%iTRO3_pas .or. jsp == art_indices%iTRO3) THEN
                    conversion_factor = 1.0_wp/( ppmv2gg * 1.e6_wp)
                ENDIF

            ENDIF
            
            tracer(:,:,:,jsp) = chem_init%chem_init_chem%spec_interp(:,:,:,jsp)*conversion_factor

            CALL message('init chemtracer', 'init tracer ' //  int2string(jsp))

          ENDIF
      END SELECT
    ENDIF

  ENDDO
 
  IF (chem_init%chem_init_in%linitialized) THEN
    CALL deallocate_chem_init_atm(chem_init)
    CALL deallocate_chem_init_chem(chem_init)
  END IF

  NULLIFY(art_atmo)

END SUBROUTINE art_init_chemtracer_ext
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_get_chemtracer_idx(jg)
!<
! SUBROUTINE art_get_chemtracer_idx
! This routine sets the indices for chemtracers.
! Part of Module: mo_art_init_chemtracer
! Author: Jennifer Schroeter, KIT
! Initial Release: around 2018-10
! Modifications:
!>
  INTEGER, INTENT(in) ::   &
    &  jg  !< patch id
  ! local variables
  INTEGER :: ierror
  TYPE(t_art_chem_indices),POINTER    :: &
      &  art_indices    !< Pointer to ART indicesetrised chem fields

  art_indices => p_art_data(jg)%chem%indices

  CALL p_art_data(jg)%dict_tracer%get('TRCHBR3',art_indices%iTRCHBR3,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCHBR3    = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCH2BR2',art_indices%iTRCH2BR2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCH2BR2   = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCH4',art_indices%iTRCH4,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCH4      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRC2H6',art_indices%iTRC2H6,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRC2H6     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRC5H8',art_indices%iTRC5H8,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRC5H8     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCH3COCH3',art_indices%iTRCH3COCH3,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCH3COCH3 = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCO',art_indices%iTRCO,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCO       = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCO2',art_indices%iTRCO2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCO2      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRC3H8',art_indices%iTRC3H8,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRC3H8     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRH2O',art_indices%iTRH2O,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRH2O      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRH2O_feed',art_indices%iTRH2O_feed,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRH2O_feed      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRO3',art_indices%iTRO3,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRO3       = 0
  CALL p_art_data(jg)%dict_tracer%get('TRN2O',art_indices%iTRN2O,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRN2O      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRNOy',art_indices%iTRNOy,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRNOy      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRNH3',art_indices%iTRNH3,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRNH3      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRSO2',art_indices%iTRSO2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRSO2      = 0
  CALL p_art_data(jg)%dict_tracer%get('TRH2SO4',art_indices%iTRH2SO4,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRH2SO4    = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCH3CN',art_indices%iTRCH3CN,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCH3CN     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRHNO3',art_indices%iTRHNO3,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRHNO3     = 0
  CALL p_art_data(jg)%dict_tracer%get('TROCS',art_indices%iTROCS,ierror)
    IF (ierror /= SUCCESS) art_indices%iTROCS     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRDMS',art_indices%iTRDMS,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRDMS     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRNO2',art_indices%iTRNO2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRNO2     = 0
  CALL p_art_data(jg)%dict_tracer%get('TRO3_pas',art_indices%iTRO3_pas,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRO3_pas  = 0
  CALL p_art_data(jg)%dict_tracer%get('TRAGE',art_indices%iTRAGE,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRAGE      = 0
  CALL p_art_data(jg)%dict_tracer%get('TR_vortex',art_indices%iTR_vortex,ierror)
    IF (ierror /= SUCCESS) art_indices%iTR_vortex  = 0
  CALL p_art_data(jg)%dict_tracer%get('TRN2O5',art_indices%iTRN2O5,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRN2O5 = 0
  CALL p_art_data(jg)%dict_tracer%get('TRHCl',art_indices%iTRHCl,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRHCl = 0
  CALL p_art_data(jg)%dict_tracer%get('TRHBr',art_indices%iTRHBr,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRHBr = 0
  CALL p_art_data(jg)%dict_tracer%get('TRHOCl',art_indices%iTRHOCl,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRHOCl = 0
  CALL p_art_data(jg)%dict_tracer%get('TRHOBr',art_indices%iTRHOBr,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRHOBr = 0
  CALL p_art_data(jg)%dict_tracer%get('TRClONO2',art_indices%iTRClONO2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRClONO2 = 0
  CALL p_art_data(jg)%dict_tracer%get('TRBrONO2',art_indices%iTRBrONO2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRBrONO2 = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCFCl3',art_indices%iTRCFCl3,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCFCl3 = 0
  CALL p_art_data(jg)%dict_tracer%get('TRCF2Cl2',art_indices%iTRCF2Cl2,ierror)
    IF (ierror /= SUCCESS) art_indices%iTRCF2Cl2 = 0
  CALL p_art_data(jg)%dict_tracer%get('TR_cold',art_indices%iTR_cold,ierror)
    IF (ierror /= SUCCESS) art_indices%iTR_cold = 0

  CALL get_regio_tracer_idx(jg,'TR_stn',art_indices%iTR_stn)
  CALL get_regio_tracer_idx(jg,'TR_stt',art_indices%iTR_stt)
  CALL get_regio_tracer_idx(jg,'TR_sts',art_indices%iTR_sts)
  CALL get_regio_tracer_idx(jg,'TR_trn',art_indices%iTR_trn)
  CALL get_regio_tracer_idx(jg,'TR_trt',art_indices%iTR_trt)
  CALL get_regio_tracer_idx(jg,'TR_trs',art_indices%iTR_trs)
  CALL get_regio_tracer_idx(jg,'TR_tiln',art_indices%iTR_tiln)
  CALL get_regio_tracer_idx(jg,'TR_tils',art_indices%iTR_tils)
  CALL get_regio_tracer_idx(jg,'TR_nh',art_indices%iTR_nh)
  CALL get_regio_tracer_idx(jg,'TR_sh',art_indices%iTR_sh)
  CALL get_regio_tracer_idx(jg,'TR_nin',art_indices%iTR_nin)
  CALL get_regio_tracer_idx(jg,'TR_sin',art_indices%iTR_sin)
  CALL get_regio_tracer_idx(jg,'TR_ech',art_indices%iTR_ech)
  CALL get_regio_tracer_idx(jg,'TR_sea',art_indices%iTR_sea)
  CALL get_regio_tracer_idx(jg,'TR_sib',art_indices%iTR_sib)
  CALL get_regio_tracer_idx(jg,'TR_eur',art_indices%iTR_eur)
  CALL get_regio_tracer_idx(jg,'TR_med',art_indices%iTR_med)
  CALL get_regio_tracer_idx(jg,'TR_naf',art_indices%iTR_naf)
  CALL get_regio_tracer_idx(jg,'TR_saf',art_indices%iTR_saf)
  CALL get_regio_tracer_idx(jg,'TR_mdg',art_indices%iTR_mdg)
  CALL get_regio_tracer_idx(jg,'TR_aus',art_indices%iTR_aus)
  CALL get_regio_tracer_idx(jg,'TR_nam',art_indices%iTR_nam)
  CALL get_regio_tracer_idx(jg,'TR_sam',art_indices%iTR_sam)
  CALL get_regio_tracer_idx(jg,'TR_tpo',art_indices%iTR_tpo)
  CALL get_regio_tracer_idx(jg,'TR_tao',art_indices%iTR_tao)
  CALL get_regio_tracer_idx(jg,'TR_tio',art_indices%iTR_tio)
  CALL get_regio_tracer_idx(jg,'TR_bgs',art_indices%iTR_bgs)
  CALL get_regio_tracer_idx(jg,'TR_bgn',art_indices%iTR_bgn)
  CALL get_regio_tracer_idx(jg,'TR_art',art_indices%iTR_art)
  CALL get_regio_tracer_idx(jg,'TR_nest',art_indices%iTR_nest)

END SUBROUTINE art_get_chemtracer_idx

SUBROUTINE get_regio_tracer_idx(jg, tracer_name, iTR)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: &
    &  jg
  CHARACTER(LEN = *), INTENT(in) :: &
    &  tracer_name
  INTEGER, INTENT(inout) :: &
    &  iTR
  ! local variables
  INTEGER :: &
    &  ierror

  CALL p_art_data(jg)%dict_tracer%get(TRIM(tracer_name),iTR,ierror)

  IF (ierror == SUCCESS)  THEN
    p_art_data(jg)%chem%param%lregio_tracers = .TRUE.
  ELSE
    iTR    = 0
  END IF


END SUBROUTINE get_regio_tracer_idx

END MODULE mo_art_init_chemtracer
