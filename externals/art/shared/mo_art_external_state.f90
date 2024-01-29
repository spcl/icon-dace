!
! mo_art_external_state
! This module provides initialization procedures for data from
! external sources (e.g. soil types, volcano characteristics,
! radioactive source characteristics, etc.)
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

MODULE mo_art_external_state
! ICON
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_exception,                     ONLY: message,message_text
  USE mo_art_config,                    ONLY: t_art_config
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_art_emiss_types,               ONLY: t_art_aero_emiss
  USE mo_var_list,                      ONLY: t_var_list_ptr
! ART
  USE mo_art_external_types,            ONLY: t_art_external
  USE mo_art_read_extdata,              ONLY: art_read_ext_soil, art_read_pollendata, &
    &                                         art_read_ext_dms, art_read_biomBurndata
  USE mo_art_external_init_soilhwsd,    ONLY: art_extinit_hwsdsoil_prepare, &
    &                                         art_extinit_hwsdsoil_finalize
  USE mo_art_emission_online_dms,       ONLY: art_extinit_dms
  USE mo_art_external_init_volc,        ONLY: art_extinit_volcanoes
  USE mo_art_external_init_pollen,      ONLY: art_extinit_pollen
  USE mo_art_external_init_biomBurn,    ONLY: art_extinit_biomBurn
  USE mo_art_read_pol_oper,             ONLY: art_read_pol_oper
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_modes_linked_list,         ONLY: t_mode
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_external_state'

  PUBLIC :: art_external_init

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_external_init(p_patch, artconf, ext_data, dims2d, dict_tracer, &
  &                          art_extdata, tracer2aeroemiss, p_prog_list)
!<
! SUBROUTINE art_external_init
! This subroutine provides an initialization structures for external
! data.
! Based on: -
! Part of Module: mo_art_external_state
! Author: Daniel Rieger, KIT
! Initial Release: 2016-12-14
! Modifications:
! 2016-12-21: Jonas Straub, KIT
! - included initilization for pollen
! YYYY-MM-DD: Carmen Ullwer, KIT
! - included initilization for DMS
!>
  TYPE(t_patch), TARGET, INTENT(in)      :: &
    &  p_patch                                !< patch on which computation is performed
  TYPE(t_art_config),INTENT(IN)          :: &
    &  artconf                                !< ART configuration state
  TYPE(t_external_data), INTENT(in)      :: &
    &  ext_data                               !< Atmosphere external data (ICON)
  INTEGER, INTENT(IN)                    :: &
    &  dims2d(2)                              !< 2D-Dimensions used to allocate fields
  TYPE(t_key_value_store), INTENT(in)    :: &
    &  dict_tracer                            !< Tracer index dictionary
  TYPE(t_art_external),INTENT(INOUT)     :: &
    &  art_extdata                            !< External data container for ART
  TYPE(t_art_aero_emiss),INTENT(IN)      :: &
    &  tracer2aeroemiss                       !< Aerosol tracer to emission dictionary
  TYPE(t_var_list_ptr), INTENT(in)           :: &
    &  p_prog_list                            !< current prognostic state list
  INTEGER ::  jg              !< patch id

  !local variable
  TYPE(t_mode),POINTER :: current
  NULLIFY(current)
  jg        = p_patch%id
  
  IF (artconf%lart_aerosol) THEN
    IF (tracer2aeroemiss%lisinit) THEN
      current=>tracer2aeroemiss%e2t_list%p%first_mode
      DO WHILE(ASSOCIATED(current))
        SELECT TYPE(this=>current%fields)
          TYPE IS(t_art_emiss2tracer)
            ! Soil initializations
              IF (this%name=='dust' .AND. this%lcalcemiss) THEN

                WRITE (message_text,*) 'ART: Initialization of HWSD soil data.'
                CALL message (TRIM(routine)//':art_external_init', message_text)

                ! Preparations (ALLOCATE, Definitions, ...)
                CALL art_extinit_hwsdsoil_prepare(dims2d, art_extdata%soil_prop)

                ! Read external data
                CALL art_read_ext_soil(jg,art_extdata%soil_prop,TRIM(artconf%cart_input_folder))

                ! Other calculations
                CALL art_extinit_hwsdsoil_finalize(p_patch,ext_data,art_extdata%soil_prop)
              ENDIF !dust%lcalcemiss

              ! Volcano initializations (get source locations)
  !!          IF (artconf%iart_volcano == 1 .OR. artconf%iart_volcano == 2) THEN
              IF (this%name=='volc' .AND. this%lcalcemiss) THEN
                WRITE (message_text,*) 'ART: Initialization of volcano data.'
                CALL message (TRIM(routine)//':art_external_init', message_text)

                CALL art_extinit_volcanoes(p_patch, artconf%cart_volcano_file, dict_tracer,  &
                  &                        this, art_extdata%volc_data, artconf%iart_volcano)
  !!          ENDIF !iart_volcano
              ENDIF !volc%lcalcemiss
        END SELECT
        current=>current%next_mode
      END DO
    END IF
  
    ! Pollen initializations (get source locations)
    IF (artconf%iart_pollen == 1) THEN
      WRITE (message_text,*) 'ART: Initialization of pollen tracer data.'
      CALL message (TRIM(routine)//':art_external_init', message_text)
    
      ! ALLOCATE, inital calculations
      CALL art_extinit_pollen(p_patch,p_patch%nblks_c, art_extdata%pollen_prop, dict_tracer, &
        &                     TRIM(artconf%cart_input_folder), artconf%lart_diag_out)
    
      ! Read external data
      ! art_read_pol_oper was called in art_diagnostics_interface in previous
      ! code versions - moved here to avoid initialization related misbehaviour
      CALL art_read_pol_oper(jg)
      CALL art_read_pollendata(jg, art_extdata%pollen_prop, TRIM(artconf%cart_input_folder))
    ENDIF !iart_pollen

    ! Biomass burning initializations ( and get source locations)
    IF (artconf%iart_fire == 1) THEN
      WRITE (message_text,*) 'ART: Initialization of biomass burning tracer data.'
      CALL message (TRIM(routine)//':art_external_init', message_text)

      ! ALLOCATE, inital calculations
      CALL art_extinit_biomBurn(p_patch%nblks_c, art_extdata%biomBurn_prop)

      ! Read external data
      CALL art_read_biomBurndata(jg ,art_extdata%biomBurn_prop, TRIM(artconf%cart_input_folder))
    ENDIF !iart_pollen
    
  ENDIF !lart_aerosol

  ! PFT initializations
  ! ...missing here as pointers are set earlier 
  ! and therefore this initialization has to be performed earlier
  
  IF (artconf%lart_chem) THEN
    WRITE (message_text,*) 'ART: Initialization of DMS tracer data.'
    CALL message (TRIM(routine)//':art_external_init', message_text)

    ! ALLOCATE, inital calculations
    IF (artconf%lart_chemtracer) THEN
      CALL art_extinit_dms(art_extdata%online_dms, jg, 'TRDMS', p_prog_list)
    END IF

    IF (artconf%lart_mecca) THEN
#ifdef __ART_GPL
      CALL art_extinit_dms(art_extdata%online_dms, jg, 'DMS', p_prog_list)
#endif
    END IF

    ! Read external data
    IF (p_art_data(jg)%ext%online_dms%lcalc_onl) THEN
      CALL art_read_ext_dms(jg, art_extdata%online_dms,TRIM(artconf%cart_input_folder))
    ELSE
      WRITE (message_text,*) 'ART: No online DMS was selected.'
      CALL message (TRIM(routine)//':art_external_init', message_text)
    ENDIF
  ENDIF  !lart_chem


END SUBROUTINE art_external_init
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_external_state
