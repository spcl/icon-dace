!
! mo_art_init
! This module provides initialization routines for ICON-ART
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

MODULE mo_art_init
! ICON
  USE mo_kind,                          ONLY: wp, dp
  USE mo_model_domain,                  ONLY: p_patch
  USE mo_grid_config,                   ONLY: start_time
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mo_run_config,                    ONLY: number_of_grid_used, ntracer
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_art_config,                    ONLY: t_art_config,art_config,IART_PATH_LEN
!++un 12.3.20
!  irad_aero of aes has to be taken into account when running aes physics
!--un  USE mo_radiation_config,              ONLY: irad_aero
  USE mo_radiation_config,              ONLY: irad_aero_nwp => irad_aero
  USE mo_aes_rad_config,                ONLY: aes_rad_config
  USE mo_impl_constants,                ONLY: inwp, iaes
  USE mo_run_config,                    ONLY: iforcing
!+un 12.3.20

  USE mo_ext_data_state,                ONLY: ext_data ! temporary!
  USE mtime,                            ONLY: timedelta, datetime
  ! ART
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom
  USE mo_art_data,                      ONLY: p_art_data,t_art_air_properties,                 &
    &                                         t_art_turb_fields, t_art_tend
  USE mo_art_chem_data,                 ONLY: t_art_chem
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_external_state,            ONLY: art_external_init
  USE mo_art_fplume_init,               ONLY: art_fplume_init
  USE mo_art_init_radiation,            ONLY: art_init_radiation
#ifdef __ART_GPL
  USE mo_art_init_full_chemistry,       ONLY: art_init_mapping_kpp_indices
#endif
!  USE mo_art_define_polar_vortex,       ONLY: art_init_polar_vortex_boundary
  USE mo_art_create_filenames,          ONLY: art_set_io_suffix
  USE mo_art_photo_init,                ONLY: art_photo_init
#ifdef __ART_GPL
  USE mo_art_psc_init,                  ONLY: art_psc_init_arrays
#endif

  USE mo_art_pntSrc_state,              ONLY: art_emiss_init_pntSrc
  USE mo_art_emiss_state,               ONLY: art_aerosol_assign_emission2tracer,  &
                                          &   art_init_emissions_from_xml
  USE mo_art_prescribed_state,          ONLY: art_create_prescr_list
#ifdef __ART_GPL
  USE messy_mecca_kpp_global,           ONLY: c
#endif

  ! --------------------------------------------
  ! ---------  online emission module ----------
  ! --------------------------------------------
  USE mo_art_oem_init,                  ONLY: art_oem_init

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_init'

  PUBLIC :: art_init

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init(jg, tc_dt_model, tc_exp_refdate, p_prog_list, tracer)
!<
! SUBROUTINE art_init
! This subroutine provides an interface to initialization routines in ICON-ART.
! Part of Module: mo_art_init
! Author: Daniel Rieger, KIT
! Initial Release: 2013-05-16
! Modifications:
! 2016-08-10: Daniel Rieger, KIT
! - Update to new tracer structure
! 2016-11-22: Michael Weimer, KIT
! - included PSC data structure allocation
! 2018-07-  : Lukas Muser, KIT
! - initialization of aerosol tracers with fix values
!>
  INTEGER,INTENT(in)                        :: &
    &  jg                          !< patch id
  TYPE(timedelta), POINTER, INTENT(in)      :: &
    &  tc_dt_model                 !< Model timestep
  TYPE(datetime), POINTER, INTENT(in)       :: &
    &  tc_exp_refdate              !< Experiment reference date
  TYPE(t_var_list_ptr),INTENT(in) :: &
    &  p_prog_list                 !< current list: prognostic
  REAL(wp), POINTER, INTENT(inout),OPTIONAL :: &
    &  tracer(:,:,:,:)             !< Tracer mixing ratios [# kg-1], [amug kg-1] or [mol mol-1]
! Local variables
  INTEGER                  ::    &
    &  idims(3)                    !< Allocation dimensions
  CHARACTER(LEN=2)         ::    &
    &  cjg                         !< patch id as character
  TYPE(t_mode),POINTER     ::    &
    &  current_mode                !< pointer to loop through mode structure
  TYPE(t_art_config), POINTER :: &
    &  artconf                     !< ART configuration (namelist setting)
  TYPE(t_art_air_properties),POINTER :: &
    &  airprop                     !< pointer to air properties
  TYPE(t_art_turb_fields),   POINTER :: &
    &  turbfields                  !< pointer to turbulence fields
  TYPE(t_art_chem),POINTER    :: &
    &  art_chem                    !< Pointer to ART chem fields
  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo                    !< Pointer to ART atmo fields
  TYPE(t_art_tend),POINTER    :: &
    &  art_tend                     !< Pointer to ART tendency fields

!++un 12.3.20 
!  irad_aero of aes has to be taken into account when running aes physics
   INTEGER :: irad_aero
  !OEM
  INTEGER :: ierror
  CHARACTER (LEN= 255) :: yerrmsg

  ! ----------------------------------
  ! --- 0.0 Some Definitions
  ! ----------------------------------

  artconf    => art_config(jg)
  airprop    => p_art_data(jg)%air_prop
  turbfields => p_art_data(jg)%turb_fields
  art_chem   => p_art_data(jg)%chem
  art_atmo   => p_art_data(jg)%atmo
  art_tend   => p_art_data(jg)%tend

  IF (jg <= 9) THEN
    write(cjg,"(I1)") jg
    cjg = '0'//cjg
  ELSE
    write(cjg,"(I2)") jg
  ENDIF

  idims(1) = art_atmo%nproma
  idims(2) = art_atmo%nlev
  idims(3) = art_atmo%nblks

  ! ----------------------------------
  ! --- 1.0 Initialize data structures
  ! ----------------------------------

  CALL message('','ART: Constructing ART Data structures for domain '//TRIM(cjg))
  
  ALLOCATE (airprop%art_free_path(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ALLOCATE (airprop%art_dyn_visc(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  
  ALLOCATE (turbfields%sv(art_atmo%nproma,art_atmo%nblks,artconf%nturb_tracer))
  turbfields%sv   = 0.0_wp
  ALLOCATE (turbfields%vdep(art_atmo%nproma,art_atmo%nblks,artconf%nturb_tracer))
  turbfields%vdep = 0.0_wp

  !$ACC ENTER DATA CREATE(airprop%art_dyn_visc, airprop%art_free_path, turbfields%sv, turbfields%vdep)

  IF (artconf%iart_aci_cold == 6 .OR. artconf%iart_aci_cold == 7) THEN
    ALLOCATE(p_art_data(jg)%diag%ndust_tot(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ENDIF


  ! initialize FPLUME structure for output
  IF (art_config(jg)%iart_fplume /= 0) THEN
    CALL art_fplume_init(p_patch(jg),p_art_data(jg)%fplume_init_all,art_config(jg)%iart_volc_numb,&
                       & art_config(jg)%cart_fplume_inp,tc_exp_refdate,tc_dt_model)
  ENDIF


  ! ----------------------------------
  ! --- 4.0 Initializations
  ! ----------------------------------
  ! This actually changes the original namelist parameter cart_io_suffix. This
  ! cannot be done within the check of the namelist parameters since
  ! number_of_grid_used is not initialised at this point
  CALL art_set_io_suffix(art_config(jg)%cart_io_suffix, number_of_grid_used(jg))

  IF (art_config(jg)%lart_chem) THEN
    ALLOCATE(p_art_data(jg)%chem%vmr2Nconc(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
    ALLOCATE(art_chem%water_tracers(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks,3))
    ALLOCATE(art_chem%CO2_mmr_depos(art_atmo%nlev))

    p_art_data(jg)%chem%vmr2Nconc = 1._wp
    art_chem%water_tracers = 0._wp

    IF ((p_art_data(jg)%chem%param%OH_chem_meta%is_init) &
      &  .OR. (art_config(jg)%lart_mecca)) THEN
      CALL art_photo_init(p_art_data(jg)%chem%photo, &
                  &       art_atmo%nproma,           &
                  &       art_atmo%nlev,             &
                  &       art_atmo%nblks)
    END IF

    IF (art_config(jg)%lart_chemtracer) THEN
      ALLOCATE(art_chem%param%o2_column(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))

      ! LK: if one does not understand, how things work ...
      !     art_checm%so2_column is allocated two times one time per add_ref
      !     used for I/O and so on than once more here.
      !     as so2_column is a pointer it gets reassociated here and the memory in
      !     add_ref lost its connection and susequently is not connected anymore.
      !
      !     THIS IS A SEVERE BUG !!!!!!!

      ALLOCATE(art_tend%chem(art_atmo%nproma, art_atmo%nlev, art_atmo%nblks, ntracer))
      art_tend%chem(:,:,:,:) = 0.0_wp

!      CALL art_init_polar_vortex_boundary(jg)
      IF (.NOT. PRESENT(tracer))  THEN
        CALL finish('mo_art_init:art_init',   &
          &         'No tracer fields passed for chemtracer initialization.')
      END IF
    END IF

    IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
      CALL art_init_mapping_kpp_indices(jg, p_prog_list, TRIM(artconf%cart_mecca_xml))

      IF (.NOT. PRESENT(tracer))  THEN
        CALL finish('mo_art_init:art_init',   &
          &         'No tracer fields passed for MECCA initialization.')
      END IF
      IF (start_time(jg) <= 0._wp) THEN
        c(:) = 0._dp
      END IF
#else
      CALL finish('mo_art_init', &
          &     'You have activated MECCA which is using software under GPL licence. Deactivate MECCA or activate GPL software for ART.')
#endif
    END IF


    IF (art_config(jg)%lart_psc) THEN
#ifdef __ART_GPL
      CALL art_psc_init_arrays(p_art_data(jg)%chem%PSC_meta, jg, tracer, p_prog_list)
#else
      CALL finish('mo_art_init', &
          &     'You have activated PSCs which depends on software using GPL licence. Deactivate PSCs or activate GPL software for ART.')
#endif
    END IF

  ENDIF !lart_chem

  IF ((artconf%lart_chem) .OR. (artconf%lart_aerosol)) THEN
    CALL art_create_prescr_list(p_art_data(jg)%prescr_list,p_prog_list, jg)
  END IF

  WRITE (message_text,*) 'ART: Initializing external properties'
  CALL message (TRIM(routine)//':art_init', message_text)
  CALL art_external_init(p_patch(jg), art_config(jg),     &
    &                    ext_data(jg), (/art_atmo%nproma, art_atmo%nblks/), &
    &                    p_art_data(jg)%dict_tracer,      &
    &                    p_art_data(jg)%ext,              &
    &                    p_art_data(jg)%tracer2aeroemiss, &
    &                    p_prog_list)

  IF(art_config(jg)%lart_pntSrc) THEN
    CALL art_emiss_init_pntSrc(p_patch(jg),                             &
      &                        art_atmo%z_ifc,                          &
      &                        tc_dt_model,                             &
      &                        tc_exp_refdate,                          &
      &                        TRIM(art_config(jg)%cart_pntSrc_xml),    &
      &                        art_config(jg)%iart_radioact,            &
      &                        TRIM(art_config(jg)%cart_radioact_file), &
      &                        p_prog_list,                             &
      &                        p_art_data(jg))
  ENDIF

  IF (art_config(jg)%lart_aerosol) THEN
!++un 12.3.20 for aerosol radiation coupling with AES physics
! irad_aero of aes has to be taken into account when running aes physics
! otherwise, the default irad_aer == 2 would used within ART and
! the value given in the namlist be ignored
    IF (iforcing == inwp) THEN
      irad_aero = irad_aero_nwp
    ELSEIF (iforcing == iaes) THEN
      irad_aero = aes_rad_config(1)%irad_aero
    ELSE
      CALL finish('','iforcing not supported with ART')
    ENDIF
!+un  12.3.20
    IF ( irad_aero == 9 .AND. artconf%iart_ari>0) THEN

        current_mode => p_mode_state(jg)%p_mode_list%p%first_mode
        DO WHILE(ASSOCIATED(current_mode))
          SELECT TYPE (fields=>current_mode%fields)
            CLASS IS (t_fields_2mom)
              CALL art_init_radiation(jg,                              &
              &                     TRIM(fields%name),               &
              &                     fields%opt_props%ext_coeff,      &
              &                     fields%opt_props%ssa_coeff,      & 
              &                     fields%opt_props%asy_coeff,      &
              &                     fields%opt_props%ext_param,      &
              &                     fields%opt_props%ssa_param,      &
              &                     fields%opt_props%asy_param,      &
              &                     fields%opt_props%ext_default,    &
              &                     fields%opt_props%ssa_default,    &
              &                     fields%opt_props%asy_default,    &
              &                     fields%opt_props%dia_min_factor, &
              &                     fields%opt_props%dia_max_factor)
            CLASS DEFAULT
          END SELECT
          current_mode=>current_mode%next_mode
      ENDDO
    ENDIF

  ENDIF !lart_aerosol

  ! --------------------------------------------
  ! ---  4.1 prescribed emissions
  ! --------------------------------------------

  IF ((artconf%lart_chem) .OR. (artconf%lart_aerosol))  THEN

    CALL art_init_emissions_from_xml(jg,p_prog_list,             &
                          &          artconf%cart_emiss_xml_file)
  END IF

  ! --------------------------------------------
  ! ---------  online emission module ----------
  ! --------------------------------------------

  IF (artconf%lart_chem) THEN

    CALL art_oem_init(jg,p_prog_list,ierror,yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF

  END IF


  ! ----------------------------------
  ! --- 5.0 Clean Up
  ! ----------------------------------

  NULLIFY(current_mode)
END SUBROUTINE art_init
!!
!!-------------------------------------------------------------------------
!!

END MODULE mo_art_init
