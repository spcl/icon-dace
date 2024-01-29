!
! Provides interface to ART-routines dealing with emissions
!
! This module provides an interface to the ART emission routines.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
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

MODULE mo_art_emission_interface
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_lnd_nwp_config,                ONLY: dzsoil
  USE mo_exception,                     ONLY: finish, message
  USE mo_parallel_config,               ONLY: nproma
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_nonhydro_state,                ONLY: p_nh_state_lists
  USE mo_ext_data_types,                ONLY: t_external_data
  USE mo_nwp_lnd_types,                 ONLY: t_lnd_diag
  USE mo_run_config,                    ONLY: lart,ntracer
  USE mo_time_config,                   ONLY: time_config
  USE mtime,                            ONLY: datetime, getDayOfYearFromDateTime
  USE mo_util_mtime,                    ONLY: mtime_utils, FMT_HHH
  USE mo_mpi,                           ONLY: p_max,p_comm_work
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_emissInt
  USE mo_fortran_tools,                 ONLY: assert_acc_device_only

  ! ART
! Infrastructure Routines
  USE mo_art_modes_linked_list,         ONLY: p_mode_state,t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom,t_fields_radio, &
                                          &   t_fields_pollen, t_fields_volc
  USE mo_art_aerosol_utilities,         ONLY: art_air_properties, art_calc_number_from_mass
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_config,                    ONLY: art_config
! Emission Routines
  USE mo_art_emission_volc_1mom,        ONLY: art_organize_emission_volc
  USE mo_art_emission_volc_2mom,        ONLY: art_prepare_emission_volc,    &
                                          &   art_calculate_emission_volc
  USE mo_art_emission_seas,             ONLY: art_seas_emiss_martensson,    &
                                          &   art_seas_emiss_monahan,       &
                                          &   art_seas_emiss_smith,         &
                                          &   art_seas_emiss_mode1,         &
                                          &   art_seas_emiss_mode2,         &
                                          &   art_seas_emiss_mode3
  USE mo_art_emission_dust,             ONLY: art_emission_dust,art_prepare_emission_dust
  USE mo_art_emission_chemtracer,       ONLY: art_emiss_chemtracer
#ifdef __ART_GPL
  USE mo_art_emission_full_chemistry,   ONLY: art_emiss_full_chemistry
#endif
  USE mo_art_emission_pollen,           ONLY: art_emiss_pollen, art_pollen_get_nstns, &
                                          &   art_prepare_tsum, art_prepare_sdes,     &
                                          &   art_prepare_saisl
  USE mo_art_emission_pntSrc,           ONLY: art_emission_pntSrc
  USE mo_art_emission_biomBurn,         ONLY: art_emission_biomBurn_prepare,art_emission_biomBurn
  USE mo_art_read_emissions,            ONLY: art_add_emission_to_tracers
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_art_prescribed_state,          ONLY: art_prescribe_tracers
#ifdef _OPENMP
  USE omp_lib
#endif
  USE mo_sync,                          ONLY: sync_patch_array_mult, SYNC_C
  USE mo_art_fplume_emission,           ONLY: art_fplume_emission
  USE mo_run_config,                    ONLY: iqv 
  USE mo_art_chem_deposition,           ONLY: art_CO2_deposition
  USE mo_art_diagnostics,               ONLY: art_save_aerosol_emission
  USE mo_art_read_extdata,              ONLY: art_read_sdes_ambrosia

  ! OEM
  USE mo_art_oem_emission,              ONLY: art_oem_compute_emissions
  USE mo_art_oem_vprm,                  ONLY: art_oem_compute_biosphere_fluxes
  USE mo_art_oem_types,                 ONLY: p_art_oem_data,  &
                                          &   t_art_oem_config

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_emission_interface'

  PUBLIC  :: art_emission_interface

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
  SUBROUTINE art_emission_interface(p_prog_list,ext_data,p_patch,dtime, &
       &                            p_diag_lnd, current_date,tracer,lacc)
  !! Interface for ART: Emissions
  !! @par Revision History
  !! Initial revision by Daniel Reinert, DWD (2012-01-27)
  !! Modification by Kristina Lundgren, KIT (2012-01-30)
  !! Rewritten by Daniel Rieger, KIT (2013-09-30)
  TYPE(t_var_list_ptr), INTENT(inout) :: &
    &  p_prog_list             !< list of prognostic variables
  TYPE(t_external_data), INTENT(in) ::  &
    &  ext_data                !< Atmosphere external data
  TYPE(t_patch), TARGET, INTENT(in) ::  &
    &  p_patch                 !< Patch on which computation is performed
  REAL(wp), INTENT(in)    :: &
    &  dtime                   !< Time step (advection)
  TYPE(t_lnd_diag), INTENT(in)      :: &
    &  p_diag_lnd              !< List of diagnostic fields (land)
  TYPE(datetime), INTENT(in), POINTER :: &
    &  current_date            !< Date and time information
  REAL(wp), INTENT(inout) :: &
    &  tracer(:,:,:,:)         !< Tracer mixing ratios [kg kg-1]
  LOGICAL, OPTIONAL, INTENT(IN) :: lacc
  ! Local variables
  TYPE(t_art_emiss2tracer),POINTER:: &
    &  this                    !< Current emiss2tracer dictionary item
  INTEGER                 :: &
    &  jg, jb, ijsp, jc,     & !< Patch id, counter for block loop, jsp loop, vertical loop
    &  istart, iend, nlev,   & !< Start and end of nproma loop, number of levels
    &  nblks,                & !< Number of blocks
    &  n_stns,               & !< variables needed for atab-readout
    &  doy_dec1, ierr,       & !< days since 1st December (for initial time of run)
    &  current_doy_dec1,     & !< days since 1st December (for current date)
    &  doy_start_season,     & !< days since 1st December for start of pollen season
    &  doy_end_season,       & !< days since 1st December for end   of pollen season
    &  ipoll,                &
    &  ierror,               &
    &  inumb_volc              !< loop counter for number of fplume volcanoes
  LOGICAL                 :: &
    &  lvolc_block
  REAL(wp),ALLOCATABLE    :: &
    &  emiss_rateM(:,:,:),   & !< Mass emission rates [UNIT m-3 s-1], UNIT might be mug, kg
    &  emiss_rate0(:,:,:),   & !< Number emission rates [m-3 s-1]
    &  saisl_stns(:)        
  CHARACTER(LEN=3)  :: hhh    !< hours since model start, e.g., "002"
  TYPE(t_mode), POINTER   :: &
    &  this_mode               !< pointer to current aerosol mode
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo                !< pointer to ART atmo fields

  ! OEM
  TYPE(t_art_oem_config), POINTER :: &
    &  oem_config    !< OEM data structure -> config
  CHARACTER (LEN= 255) :: yerrmsg

  ! Get the number of blocks
  nblks  = p_patch%nblks_c

  ! --- Get the loop indizes
  jg         = p_patch%id

  ! --- Initialize local variables
  n_stns           = 0
  doy_dec1         = 0
  current_doy_dec1 = 0
  doy_start_season = 0
  doy_end_season   = 0
  ipoll            = 0
  this_mode        => NULL()
  art_atmo         => NULL()

  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_emissInt)

    art_atmo => p_art_data(jg)%atmo
    nlev = art_atmo%nlev

    IF (art_config(jg)%lart_pntSrc) THEN
      ! Point sources
      CALL art_emission_pntSrc(p_art_data(jg)%pntSrc, current_date, dtime, art_atmo%rho,  &
        &                      art_atmo%cell_area, art_atmo%dz, tracer)
    ENDIF

    IF (art_config(jg)%lart_aerosol .OR. art_config(jg)%lart_chem) THEN

      CALL art_prescribe_tracers(tracer, p_art_data(jg)%prescr_list,     &
               &                 p_patch, current_date)

      IF (p_art_data(jg)%emiss%is_init) THEN
        CALL art_add_emission_to_tracers(tracer,p_art_data(jg)%emiss,p_patch,   &
                                    &  dtime,                                   &
                                    &  current_date)
      ENDIF

      ! --------------------------------------------
      ! ---------  online emission module ----------
      ! --------------------------------------------
      oem_config => p_art_oem_data%configure
      IF (oem_config%emis_tracer>0) THEN
        CALL art_oem_compute_emissions(tracer,p_patch,dtime,current_date,ierror,yerrmsg)
        IF (ierror /= 0) THEN
          CALL finish ('art_emission_interface', yerrmsg)
        ENDIF
      ENDIF
      IF (oem_config%vprm_tracer>0) THEN
        CALL art_oem_compute_biosphere_fluxes(tracer,p_patch,dtime,current_date,ierror,yerrmsg)
        IF (ierror /= 0) THEN
          CALL finish ('art_vprm_interface', yerrmsg)
        ENDIF
      ENDIF

    ENDIF

    IF (art_config(jg)%lart_aerosol) THEN
!$omp parallel do default (shared) private(jb, istart, iend)
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, istart, iend)

        CALL art_air_properties(art_atmo%pres(:,:,jb),art_atmo%temp(:,:,jb), &
!          &                     istart,iend,1,art_atmo%nlev,jb,p_art_data(jg))
          &                     istart,iend,1,nlev,p_art_data(jg)%air_prop%art_free_path(:,:,jb), &
          &                     p_art_data(jg)%air_prop%art_dyn_visc(:,:,jb),lacc=lacc)
        
        ! ----------------------------------
        ! --- Preparations for emission routines
        ! ----------------------------------
        IF (p_art_data(jg)%tracer2aeroemiss%lisinit) THEN
          this_mode=>p_art_data(jg)%tracer2aeroemiss%e2t_list%p%first_mode
          DO WHILE(ASSOCIATED(this_mode))
            SELECT TYPE(this=>this_mode%fields)
              TYPE IS(t_art_emiss2tracer)
                IF(this%lcalcemiss) THEN
                  SELECT CASE(this%name)
                    CASE('dust')
                      CALL art_prepare_emission_dust(jb, istart, iend, art_atmo%u(:,nlev,jb),&
                        &                            art_atmo%v(:,nlev,jb),                  &
                        &                            art_atmo%rho(:,nlev,jb), art_atmo%tcm(:,jb),   &
                        &                            dzsoil(1), p_diag_lnd%w_so(:,1,jb),            &
                        &                            p_diag_lnd%w_so_ice(:,1,jb),                   &
                        &                            p_art_data(jg)%ext%soil_prop,                  &
                        &                            p_diag_lnd%h_snow(:,jb),                       &
                        &                            this%dg3(1,1), this%dg3(2,1), this%dg3(3,1),   &
                        &                            art_config(jg)%lart_diag_out,                  &
                        &                            p_art_data(jg)%diag%ustar_threshold(:,jb),     &
                        &                            p_art_data(jg)%diag%ustar(:,jb))
                    CASE('volc')
                      CALL art_prepare_emission_volc(current_date,jb,nlev,                          &
                        &                            art_atmo%z_ifc(:,:,jb),              &
                        &                            p_art_data(jg)%dict_tracer,                    &
                        &                            p_art_data(jg)%ext%volc_data,                  &
                        &                            this%itr3(1,1), this%itr3(2,1), this%itr3(3,1))
                    CASE DEFAULT
                      !nothing to do
                  END SELECT
                END IF !this%lcalcemiss
            END SELECT
            this_mode=>this_mode%next_mode
          END DO
          NULLIFY(this_mode)
        END IF
        ! TODO: Incorporate this into emiss2tracer-scheme above 
        SELECT CASE(art_config(jg)%iart_fire)
          CASE(0)
            ! Nothing to do, no biomass burning emissions
          CASE(1)
            CALL art_emission_biomBurn_prepare(                                               &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_crop_irrig),   &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_crop_rain),    &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_crop_mos),     &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_veg_mos),      &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_b_eg),  &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_b_d),   &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_woodland),     &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_n_eg),  &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_n_d),   &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_bn),    &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_mos),    &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),     &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),        &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),        &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),       &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_rf),    &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_forest_pf),    &
              &          ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass_rf),     &
              &          p_art_data(jg)%ext%biomBurn_prop%dc_hflux_min_res(:,:,:),            &
              &          p_art_data(jg)%ext%biomBurn_prop%dc_hflux_max_res(:,:,:),            &
              &          p_art_data(jg)%ext%biomBurn_prop%dc_burnt_area_res(:,:,:),           &
              &          p_art_data(jg)%ext%biomBurn_prop%dc_emis_res(:,:,:),                 &
              &          jb, istart, iend )
          CASE default
            CALL finish('mo_art_emission_interface:art_emission_interface', &
                 &      'ART: Unknown biomass burning emissions configuration')
        END SELECT 
      ENDDO !jb   
!$omp end parallel do

        ! ----------------------------------
        ! --- Call the emission routines
        ! ----------------------------------
!! NOT !$omp parallel do default (shared) private (jb, jc, istart, iend, emiss_rateM, emiss_rate0, dz)
!$omp parallel do default (shared) private (jb, jc, istart, iend, emiss_rateM, emiss_rate0)
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)

          IF (p_art_data(jg)%tracer2aeroemiss%lisinit) THEN
            this_mode=>p_art_data(jg)%tracer2aeroemiss%e2t_list%p%first_mode
            p_art_data(jg)%diag%emiss(:,jb,:) = 0.0_wp
            DO WHILE(ASSOCIATED(this_mode))
              SELECT TYPE(this=>this_mode%fields)
                TYPE IS(t_art_emiss2tracer)
                  IF(this%lcalcemiss) THEN
                    ALLOCATE(emiss_rateM(istart:iend,nlev,this%nmodes))
                    ALLOCATE(emiss_rate0(istart:iend,nlev,this%nmodes))
                    SELECT CASE(this%name)
                      CASE('seas_martensson') !CASE('seas')
                        ! Sea salt emission Martensson et al. (sodium and chloride)
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_seas_emiss_martensson(art_atmo%u_10m(:,jb), art_atmo%v_10m(:,jb),  &
                          &                            art_atmo%dz(:,nlev,jb), p_diag_lnd%t_s(:,jb),&
                          &                            ext_data%atm%fr_land(:,jb),                  &
                          &                            p_diag_lnd%fr_seaice(:,jb),                  &
                          &                            ext_data%atm%fr_lake(:,jb),                  &
                          &                            istart, iend, emiss_rateM(:,nlev,1))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,     &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend,   &
                          &                            nlev, nlev, jb)
                      CASE('seas_monahan')
                        ! Sea salt emission Monahan et al. (sodium and chloride)
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_seas_emiss_monahan(art_atmo%u_10m(:,jb), art_atmo%v_10m(:,jb),     &
                          &                         art_atmo%dz(:,nlev,jb),                         &
                          &                         ext_data%atm%fr_land(:,jb),                     &
                          &                         p_diag_lnd%fr_seaice(:,jb),                     &
                          &                         ext_data%atm%fr_lake(:,jb),                     &
                          &                         istart, iend, emiss_rateM(:,nlev,1))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,     &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend,   &
                          &                            nlev, nlev, jb)
                      CASE('seas_smith')
                        ! Sea salt emission Smith et al. (unspecified)
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_seas_emiss_smith(art_atmo%u_10m(:,jb), art_atmo%v_10m(:,jb),       &
                          &                       art_atmo%dz(:,nlev,jb),ext_data%atm%fr_land(:,jb),&
                          &                       p_diag_lnd%fr_seaice(:,jb),                       &
                          &                       ext_data%atm%fr_lake(:,jb),                       &
                          &                       istart, iend, emiss_rateM(:,nlev,1))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,     &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend,   &
                          &                            nlev, nlev, jb)
                      CASE('seas_mode1')
                        ! Sea salt emission  (sodium and chloride)
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_seas_emiss_mode1(art_atmo%u_10m(:,jb),                           &
                          &                       art_atmo%v_10m(:,jb),art_atmo%dz(:,art_atmo%nlev,jb),     &
                          &                       p_diag_lnd%t_s(:,jb),ext_data%atm%fr_land(:,jb),          &
                          &                       p_diag_lnd%fr_seaice(:,jb), ext_data%atm%fr_lake(:,jb),   &
                          &                       istart,iend,emiss_rateM(:,nlev,1))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,     &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend,   &
                          &                            nlev, nlev, jb)
                      CASE('seas_mode2')
                        ! Sea salt emission  (sodium and chloride)
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_seas_emiss_mode2(art_atmo%u_10m(:,jb),                           &
                          &                       art_atmo%v_10m(:,jb),art_atmo%dz(:,art_atmo%nlev,jb),     &
                          &                       p_diag_lnd%t_s(:,jb),ext_data%atm%fr_land(:,jb),          &
                          &                       p_diag_lnd%fr_seaice(:,jb), ext_data%atm%fr_lake(:,jb),   &
                          &                       istart,iend,emiss_rateM(:,nlev,1))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,     &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend,   &
                          &                            nlev, nlev, jb)
                      CASE('seas_mode3')
                        ! Sea salt emission  (sodium and chloride)
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_seas_emiss_mode3(art_atmo%u_10m(:,jb),                           &
                          &                       art_atmo%v_10m(:,jb),art_atmo%dz(:,art_atmo%nlev,jb),     &
                          &                       p_diag_lnd%t_s(:,jb),ext_data%atm%fr_land(:,jb),          &
                          &                       p_diag_lnd%fr_seaice(:,jb), ext_data%atm%fr_lake(:,jb),   &
                          &                       istart,iend,emiss_rateM(:,nlev,1))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,     &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend,   &
                          &                            nlev, nlev, jb)
                      CASE('dust')
                        ! Mineral dust emission Rieger et al. 2017
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_emission_dust(art_atmo%dz(:,nlev,jb),                             &
                          &  ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub_eg),      &
                          &  ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_shrub),         &
                          &  ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_grass),         &
                          &  ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_bare_soil),     &
                          &  ext_data%atm%lu_class_fraction(:,jb,ext_data%atm%i_lc_sparse),        &
                          &  jb, istart, iend, p_art_data(jg)%ext%soil_prop, emiss_rateM(:,nlev,:))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend,    &
                          &                             nlev, nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:), &
                          &                   art_atmo%rho(:,:,jb), dtime, istart, iend, nlev,    &
                          &                   nlev, jb)

                      CASE('volc')
                        ! Volcanic ash emission Rieger et al. 2015 extended by 2mom-PSD
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
            
                        CALL art_calculate_emission_volc(jb, art_atmo%dz(:,:,jb),                   &
                          &                              p_patch%cells%area(:,jb),                  &
                          &                              p_art_data(jg)%ext%volc_data, nlev,        &
                          &                              this%nmodes, this%itr3(:,1),               &
                          &                              emiss_rateM(:,:,:))
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend, 1,  &
                          &                             nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend, 1,&
                          &                            nlev,jb)
                      CASE('volc_fplume')
                        ! Is there a volcano at all and is it in the current block?
                        lvolc_block = .FALSE.
                        DO inumb_volc = 1,art_config(jg)%iart_volc_numb
                          IF(p_art_data(jg)%fplume_init_all%p(inumb_volc)%ithis_nlocal_pts>=1) THEN
                            lvolc_block=.TRUE.
                          ENDIF
                        ENDDO

                        IF( .NOT. lvolc_block) THEN
                          this_mode => this_mode%next_mode
                          IF(ALLOCATED(emiss_rateM)) DEALLOCATE(emiss_rateM)
                          IF(ALLOCATED(emiss_rate0)) DEALLOCATE(emiss_rate0)
                          CYCLE
                        END IF
                        DO inumb_volc = 1,art_config(jg)%iart_volc_numb
                          IF (p_art_data(jg)%fplume_init_all%p(inumb_volc)%tri_iblk_loc == jb) THEN
                            emiss_rateM(:,:,:) = 0.0_wp
                            emiss_rate0(:,:,:) = 0.0_wp
                            CALL art_fplume_emission(p_art_data(jg),                                   &
                              &                      inumb_volc,  &
                              &                      art_atmo%z_mc(:,:,jb), art_atmo%z_ifc(:,:,jb),    &
                              &                      art_atmo%rho(:,:,jb), art_atmo%pres(:,:,jb),      &
                              &                      art_atmo%temp(:,:,jb), art_atmo%u(:,:,jb),        &
                              &                      art_atmo%v(:,:,jb), iqv, jg, art_atmo%nlev,       &
                              &                      current_date, p_art_data(jg)%ext%volc_fplume,     &
                              &                      dtime, art_atmo%dz(:,:,jb),p_patch%cells%area(:,jb),&
                              &                      tracer, emiss_rateM(:,:,:))
                            CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend, 1,  &
                              &                             nlev)
                            CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),  &
                              &                            art_atmo%rho(:,:,jb), dtime, istart, iend, 1,&
                              &                            nlev, jb)
                          ENDIF
                        ENDDO

                      CASE('biomass_burning')
                        emiss_rateM(:,:,:) = 0.0_wp
                        emiss_rate0(:,:,:) = 0.0_wp
                        CALL art_emission_biomBurn(                                               &
                                     !dimensions:(nproma,nlev,nblks)
                          &          art_atmo%temp(:,:,jb),                                       &
                          &          art_atmo%pres(:,:,jb),                                       &
                          &          art_atmo%u(:,:,jb),                                          &
                          &          art_atmo%v(:,:,jb),                                          &
                          &          tracer(:,:,jb,iqv),                                          &
                          &          art_atmo%z_mc(:,:,jb),                                       &
                          &          art_atmo%z_ifc(:,:,jb),                                      &
                          &          art_atmo%lon(:,jb),                                          &
                          &          art_atmo%cell_area(:,jb),                                    &
                          &          art_atmo%dz(:,:,jb),                                         &
                          &          current_date,                                                &
                          &          p_art_data(jg)%ext%biomBurn_prop%dc_hflux_min_res(:,:,:),    &
                          &          p_art_data(jg)%ext%biomBurn_prop%dc_hflux_max_res(:,:,:),    &
                          &          p_art_data(jg)%ext%biomBurn_prop%dc_burnt_area_res(:,:,:),   &
                          &          p_art_data(jg)%ext%biomBurn_prop%dc_emis_res(:,:,:),         &
                          &          p_art_data(jg)%ext%biomBurn_prop%flux_bc(:,jb),              &
                          &          emiss_rateM(:,:,1),                                          &
                          &          art_atmo%nlev, art_atmo%nlevp1, jb, istart, iend )
                        CALL this%calc_number_from_mass(emiss_rateM, emiss_rate0, istart, iend, 1,  &
                          &                             nlev)
                        CALL this%distribute_emissions(emiss_rateM, emiss_rate0, tracer(:,:,:,:),   &
                          &                            art_atmo%rho(:,:,jb), dtime, istart, iend, 1,&
                          &                            nlev, jb)
                      CASE DEFAULT
                        !nothing to do
                    END SELECT
                    CALL art_save_aerosol_emission(this, p_art_data(jg),                           &
                      & emiss_rateM(:,:,:), emiss_rate0(:,:,:), art_atmo%dz(:,:,:),                &
                      & dtime,  jb, istart, iend, 1, art_atmo%nlev)
                    IF(ALLOCATED(emiss_rateM)) DEALLOCATE(emiss_rateM)
                    IF(ALLOCATED(emiss_rate0)) DEALLOCATE(emiss_rate0)
                  END IF !this%lcalcemiss
              END SELECT
              this_mode=>this_mode%next_mode
            END DO
            NULLIFY(this_mode)
          END IF
        ENDDO !jb
!$omp end parallel do

        this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

        DO WHILE(ASSOCIATED(this_mode))
          ! Check how many moments the mode has
          SELECT TYPE (fields=>this_mode%fields)

            TYPE is (t_fields_2mom)
              ! handled above

            TYPE is (t_fields_pollen)

              CALL assert_acc_device_only('art_emission_interface', lacc)

              IF (art_config(jg)%iart_pollen == 0) THEN
                CALL message('art_emission_interface:art_emission_interface',  &
                  &          'pollen provided in XML, but iart_pollen is set to 0')
                this_mode => this_mode%next_mode
                CYCLE
              END IF

              ! days since 1st December (= first day) for current run (initial time)
              IF (time_config%tc_exp_startdate%date%month == 12) THEN
                doy_dec1 = time_config%tc_exp_startdate%date%day
              ELSE
                doy_dec1 = getDayOfYearFromDateTime(time_config%tc_exp_startdate, ierr) + 31
              END IF

              ! define day of year of start and end of pollen season
              ! JF: numbers are now defined in the modes.xml file
! JF:               SELECT CASE(TRIM(fields%name))
! JF:               CASE ('pollalnu', 'pollbetu')
! JF:                 doy_start_season = 1
! JF:                 doy_end_season   = 146
! JF:               CASE ('pollpoac')
! JF:                 doy_start_season = 60
! JF:                 doy_end_season   = 305
! JF:               CASE ('pollambr')
! JF:                 doy_start_season = 213
! JF:                 doy_end_season   = 280
! JF:               END SELECT

              ! define start and end of pollen season, counting days from 1st December
              doy_start_season = fields%doy_start_season + 31
              doy_end_season   = fields%doy_end_season   + 31

              ! date check: day of year of current run (initial time) is in pollen season
              IF ( doy_dec1 >= doy_start_season .AND. doy_dec1 <= doy_end_season ) THEN

                !----------------------------------------------------------------------------------
                !--   Run calculations once a day: at 12 UTC --------------------------------------
                !----------------------------------------------------------------------------------
                ! time check
                IF ( current_date%time%hour == 12 .AND.  &
                  &  (current_date%time%minute * 60 + current_date%time%second) < INT(dtime) ) THEN
                 ! Get n_stns
                 CALL art_pollen_get_nstns( p_art_data(jg)%ext%pollen_prop, fields%name, n_stns )

                 DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                    CALL art_get_indices_c(jg, jb, istart, iend)

                    !$ACC UPDATE HOST(art_atmo%t_2m)

                    ! maybe has to be moved to preparation-block
                    CALL art_prepare_tsum( current_date,                   &
                      &                    fields%name,                    &
                      &                    p_art_data(jg)%ext%pollen_prop, &
                      &                    art_atmo%t_2m(:,jb),            &
                      &                    jb, istart, iend                )
                  ENDDO !jb

                  IF(.NOT.ALLOCATED(saisl_stns)) ALLOCATE(saisl_stns(n_stns))
                  saisl_stns = 0._wp

                  IF ( TRIM(fields%name) /= 'pollpoac' .AND.  &
                    &  TRIM(fields%name) /= 'pollambr') THEN

                    CALL art_prepare_saisl( p_art_data(jg)%ext%pollen_prop, &
                      &                     current_date,                   &
                      &                     fields%name,                    &
                      &                     saisl_stns )

                    !synchronization across processors
                    !using p_max-routine since only one processor will provide a useful value
                    !that is greater than 0 (per station)
                    saisl_stns = p_max(saisl_stns,comm=p_comm_work)

                  END IF

                  IF (TRIM(fields%name) /= 'pollambr') THEN

                    DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                      CALL art_get_indices_c(jg, jb, istart, iend)

                      CALL art_prepare_sdes( p_art_data(jg)%ext%pollen_prop, &
                        &                    p_patch,                        &
                        &                    jb, istart, iend,               &
                        &                    fields%name,                    &
                        &                    saisl_stns )
                    END DO ! jb

                  END IF

                  DEALLOCATE(saisl_stns)

                END IF ! time check 12 UTC

                ! read daily static ambrosia sdes files each day at 0 UTC
                ! during model start, the AMBRsdes is read from the POV file.
                ! For the following days, the code below is used
                IF ( TRIM(fields%name) == 'pollambr' ) then
                  ! read AMBRsdes at 00 UTC but not after the start of the run
                  ! therefore we check, that the hours since model start are > 0
                  hhh = TRIM(mtime_utils%ddhhmmss(time_config%tc_exp_startdate, &
                         &                        current_date, FMT_HHH))
                  IF (hhh /= "000" .AND. current_date%time%hour == 0 .AND.  &
                  &  (current_date%time%minute * 60 + current_date%time%second) < INT(dtime) ) THEN

                    ! check if current_date is within the pollen saison. 
                    ! Otherwise no sdes file exists.
                    IF (current_date%date%month == 12) THEN
                      current_doy_dec1 = current_date%date%day
                    ELSE
                      current_doy_dec1 = getDayOfYearFromDateTime(current_date, ierr) + 31
                    END IF

                    IF ( current_doy_dec1 >= doy_start_season .AND. &
                      &  current_doy_dec1 <= doy_end_season ) THEN
                      CALL art_read_sdes_ambrosia(jg,  p_art_data(jg)%ext%pollen_prop, &
                        &      TRIM(art_config(jg)%cart_input_folder), current_date)
                    END IF

                    ! set sum of radiation to zero at midnight
                    CALL p_art_data(jg)%ext%pollen_prop%dict_pollen%get(fields%name, ipoll, ierr)
                    IF(ierr /= SUCCESS) CALL finish ('mo_art_emission_interface:'//               &
                      &   'art_emission_interface', 'ipoll not found in pollen table dictionary.')

                    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
                    !$ACC LOOP GANG VECTOR COLLAPSE(2)
                    DO jb = 1, nblks
                      DO jc = 1, nproma
                        p_art_data(jg)%ext%pollen_prop%pollen_type(ipoll)%sobs_sum(jc,jb) = 0._wp
                        p_art_data(jg)%ext%pollen_prop%pollen_type(ipoll)%rh_sum(jc,jb)   = 0._wp
                      END DO
                    END DO
                    !$ACC END PARALLEL

                  END IF ! time check 00 UTC but not model start
                END IF ! Ambrosia

                DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                  CALL art_get_indices_c(jg, jb, istart, iend)

                  CALL art_emiss_pollen(dtime, current_date,                &
                    &                   art_atmo%rho(:,art_atmo%nlev,jb),   &
                    &                   fields%name,                        &
                    &                   p_art_data(jg)%ext%pollen_prop,     &
                    &                   tracer(:,art_atmo%nlev,jb,:),       &
                    &                   p_art_data(jg)%dict_tracer,         &
                    &                   art_atmo%temp(:,art_atmo%nlev,jb),  &
                    &                   art_atmo%tke(:,art_atmo%nlev,jb),   &
                    &                   art_atmo%rain_gsp_rate(:,jb),       &
                    &                   art_atmo%rain_con_rate(:,jb),       &
                    &                   art_atmo%rh_2m(:,jb),               &
                    &                   art_atmo%swflxsfc(:,jb),            &
                    &                   art_atmo%dz(:,art_atmo%nlev,jb),    &
                    &                   art_atmo%llsm(:,jb),                &
                    &                   jb, istart, iend, lacc=lacc)

                END DO !jb

              ENDIF ! date check

            CLASS is (t_fields_radio)
              ! handled via pntSrc structure

            CLASS is (t_fields_volc)
              ! handled below
            CLASS default
              CALL finish('mo_art_emission_interface:art_emission_interface',                     &
                &         'ART: Unknown mode field type')
          END SELECT
          this_mode => this_mode%next_mode
        ENDDO ! while(associated)

      ! volcano emissions
      IF (art_config(jg)%iart_volcano == 1) THEN
        CALL art_organize_emission_volc(p_patch, current_date, dtime,art_atmo%rho,  &
          &                             p_art_data(jg)%dict_tracer,                 &
          &                             p_art_data(jg)%ext%volc_data,tracer)
      ENDIF


    ENDIF !lart_aerosol

    ! ----------------------------------
    ! --- emissions of chemical tracer
    ! ----------------------------------

    IF (art_config(jg)%lart_chem) THEN


      IF (p_art_data(jg)%chem%indices%iTRCO2 /= 0) THEN
        CALL art_CO2_deposition(jg,tracer(:,:,:,p_art_data(jg)%chem%indices%iTRCO2),  &
                   &            dtime, p_art_data(jg)%atmo)
      END IF
     

      IF (art_config(jg)%lart_chemtracer) THEN
        CALL art_emiss_chemtracer(current_date,                   &
          &                       dtime,                          &
          &                       tracer,                         &
          &                       jg,                             &
          &                       p_prog_list)
      END IF


      IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
        CALL art_emiss_full_chemistry(current_date,                   &
          &                           dtime,                          &
          &                           tracer,                         &
          &                           jg,                             &
          &                           p_prog_list)
#endif
      END IF
    ENDIF !lart_chem
    CALL sync_patch_array_mult(SYNC_C, p_patch, ntracer,  f4din=tracer(:,:,:,:))

    IF (timers_level > 3) CALL timer_stop(timer_art_emissInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)

  ENDIF !lart

END SUBROUTINE art_emission_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_interface
