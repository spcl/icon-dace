! This module checks the read-in namelist parameters and, in case of
! inconsistencies, it tries to correct these.
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

MODULE mo_nml_crosscheck


  USE, INTRINSIC :: iso_c_binding, ONLY: c_int64_t
  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_impl_constants,           ONLY: inwp, tracer_only,                                &
    &                                    iaes, RAYLEIGH_CLASSIC, INOFORCING,               &
    &                                    icosmo, MODE_IAU, MODE_IAU_OLD,                   &
    &                                    max_echotop, max_wshear, max_srh,                 &
    &                                    LSS_JSBACH, ivdiff, IHELDSUAREZ, ILDF_DRY
  USE mo_time_config,              ONLY: time_config, dt_restart
  USE mo_extpar_config,            ONLY: itopo
  USE mo_io_config,                ONLY: dt_checkpoint, lnetcdf_flt64_output, echotop_meta,&
    &                                    wshear_uv_heights, n_wshear, srh_heights, n_srh
  USE mo_parallel_config,          ONLY: check_parallel_configuration,                     &
    &                                    ignore_nproma_use_nblocks_c,                      &
    &                                    ignore_nproma_use_nblocks_e,                      &
    &                                    num_io_procs,                                     &
    &                                    num_prefetch_proc, use_dp_mpi2io, num_io_procs_radar
  USE mo_limarea_config,           ONLY: latbc_config, LATBC_TYPE_CONST, LATBC_TYPE_EXT
  USE mo_master_config,            ONLY: isRestart
  USE mo_run_config,               ONLY: nsteps, dtime, iforcing, output_mode,             &
    &                                    ltransport, ltestcase, ltimer,                    &
    &                                    activate_sync_timers, timers_level, lart,         &
    &                                    msg_level, luse_radarfwo
  USE mo_dynamics_config,          ONLY: ldeepatmo, lmoist_thdyn
  USE mo_advection_config,         ONLY: advection_config
  USE mo_nonhydrostatic_config,    ONLY: itime_scheme_nh => itime_scheme,                  &
    &                                    rayleigh_type, ivctype, iadv_rhotheta
  USE mo_atm_phy_nwp_config,       ONLY: atm_phy_nwp_config, icpl_aero_conv, iprog_aero,   &
    &                                    icpl_aero_ice
  USE mo_lnd_nwp_config,           ONLY: ntiles_lnd, lsnowtile, sstice_mode, llake
  USE mo_aes_phy_config,           ONLY: aes_phy_config
  USE mo_aes_vdf_config,           ONLY: aes_vdf_config
  USE mo_radiation_config,         ONLY: irad_aero, iRadAeroNone, iRadAeroConst,           &
    &                                    iRadAeroTegen, iRadAeroART, iRadAeroConstKinne,   &
    &                                    iRadAeroCAMSclim, iRadAeroCAMStd,                 &
    &                                    iRadAeroKinne, iRadAeroVolc, iRadAeroKinneVolc,   &
    &                                    iRadAeroKinneVolcSP, iRadAeroKinneSP,             &
    &                                    irad_o3, irad_h2o, irad_co2, irad_ch4,            &
    &                                    irad_n2o, irad_o2, irad_cfc11, irad_cfc12,        &
    &                                    icld_overlap, ecrad_llw_cloud_scat, isolrad,      &
    &                                    ecrad_iliquid_scat, ecrad_iice_scat,              &
    &                                    ecrad_isnow_scat, ecrad_irain_scat,               &
    &                                    ecrad_igraupel_scat,                              & 
    &                                    ecrad_isolver, ecrad_igas_model,                  &
    &                                    ecrad_use_general_cloud_optics
  USE mo_turbdiff_config,          ONLY: turbdiff_config
  USE mo_initicon_config,          ONLY: init_mode, dt_iau, ltile_coldstart, timeshift,    &
    &                                    itype_vert_expol, iterate_iau
  USE mo_nh_testcases_nml,         ONLY: nh_test_name, layer_thickness
  USE mo_meteogram_config,         ONLY: meteogram_output_config, check_meteogram_configuration
  USE mo_grid_config,              ONLY: lplane, n_dom, l_limited_area, start_time,        &
    &                                    nroot, is_plane_torus, n_dom_start, l_scm_mode,   &
    &                                    vct_filename
  USE mo_art_config,               ONLY: art_config
  USE mo_time_management,          ONLY: compute_timestep_settings,                        &
    &                                    compute_restart_settings,                         &
    &                                    compute_date_settings
  USE mo_event_manager,            ONLY: initEventManager
  USE mtime,                       ONLY: getTotalMilliSecondsTimeDelta, datetime,          &
    &                                    newDatetime, deallocateDatetime
  USE mo_sleve_config,             ONLY: itype_laydistr, flat_height, top_height
  USE mo_nudging_config,           ONLY: nudging_config, indg_type
  USE mo_nwp_tuning_config,        ONLY: itune_gust_diag
  USE mo_nudging_nml,              ONLY: check_nudging
  USE mo_upatmo_config,            ONLY: check_upatmo
  USE mo_name_list_output_config,  ONLY: is_variable_in_output_dom
  USE mo_coupling_config,          ONLY: is_coupled_to_ocean, is_coupled_to_waves, is_coupled_to_hydrodisc

  USE mo_assimilation_config,      ONLY: assimilation_config
  USE mo_scm_nml,                  ONLY: i_scm_netcdf, scm_sfc_temp, scm_sfc_qv, scm_sfc_mom
#ifndef __NO_ICON_LES__
  USE mo_ls_forcing_nml,           ONLY: is_ls_forcing
#endif
#ifdef __ICON_ART
  USE mo_grid_config,              ONLY: lredgrid_phys
#endif

#ifdef HAVE_RADARFWO
  USE radar_data,            ONLY: ndoms_max_radar => ndoms_max
#endif

  USE mo_sppt_config,              ONLY: sppt_config, crosscheck_sppt


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atm_crosscheck

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_nml_crosscheck"


CONTAINS

  SUBROUTINE atm_crosscheck

    INTEGER  :: jg, i
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::atm_crosscheck'
    REAL(wp) :: secs_restart, secs_checkpoint, secs_iau_end
    TYPE(datetime), POINTER :: reference_dt
    LOGICAL  :: l_global_nudging
    INTEGER(c_int64_t) :: msecs_restart

    !--------------------------------------------------------------------
    ! Compute date/time/time step settings
    ! and initialize the event manager
    !--------------------------------------------------------------------
    !
    ! Note that the ordering of the following three calls must not be
    ! changed, since they rely on previous results:
    !
    CALL compute_timestep_settings()
    CALL compute_restart_settings()
    CALL compute_date_settings("atm", dt_restart, nsteps)
    !
    ! Create an event manager, ie. a collection of different events
    !
    CALL initEventManager(time_config%tc_exp_refdate)

    !--------------------------------------------------------------------
    ! Parallelization
    !--------------------------------------------------------------------
    CALL check_parallel_configuration()

    !
    ! nblocks_c or nblocks_e does not work with nesting
    ! It crashes and would be a waste of memory, if the nest were significant smaller than the parent.
    !
    IF ((ignore_nproma_use_nblocks_c .OR. ignore_nproma_use_nblocks_e) .AND. (n_dom > 1)) CALL finish(routine, &
      'Currently nblocks_c or nblocks_e (>0) is not supported for nested domains.')

    !--------------------------------------------------------------------
    ! Limited Area Mode and LatBC read-in:
    !--------------------------------------------------------------------


    IF (lplane) CALL finish(routine,&
     'Currently a plane version is not available')

    ! Reset num_prefetch_proc to zero if the model does not run in limited-area mode
    ! or in global nudging mode or if there are no lateral boundary data to be read
    l_global_nudging = ANY(nudging_config(1:n_dom)%nudge_type == indg_type%globn)
    IF (.NOT. (l_limited_area .OR. l_global_nudging) .OR. latbc_config%itype_latbc == 0) THEN
      IF (num_prefetch_proc /=0) CALL message(routine,' WARNING! num_prefetch_proc reset to 0 !')
      num_prefetch_proc = 0
    ENDIF

    ! If LatBC data is unavailable: Idle-wait-and-retry
    !
    IF (latbc_config%nretries > 0) THEN
      !
      ! ... only supported for prefetching LatBC mode
      IF (.NOT. (l_limited_area .OR. l_global_nudging) .OR. (num_prefetch_proc == 0)) THEN
        CALL finish(routine, "LatBC: Idle-wait-and-retry only supported for prefetching LatBC mode!")
      END IF
      !
      ! ... only supported for MPI parallel configuration
#ifdef NOMPI
      CALL finish(routine, "LatBC: Idle-wait-and-retry requires MPI!")
#endif
    END IF

    ! Limited area mode must not be enabled for torus grid:
    IF (is_plane_torus .AND. l_limited_area) THEN
      CALL finish(routine, 'Plane torus grid requires l_limited_area = .FALSE.!')
    END IF

    ! Root bisection "0" does not make sense for limited area mode; it
    ! is more likely that the user tried to use a torus grid here:
    IF (l_limited_area .AND. (nroot == 0)) THEN
      CALL finish(routine, "Root bisection 0 does not make sense for limited area mode; did you try to use a torus grid?")
    END IF


    !--------------------------------------------------------------------
    ! Grid and dynamics
    !--------------------------------------------------------------------

    IF (lplane) CALL finish( routine,&
      'Currently a plane version is not available')

    !--------------------------------------------------------------------
    ! If ltestcase is set to .FALSE. in run_nml set testcase name to empty
    ! (in case it is still set in the run script)
    IF (.NOT. ltestcase) THEN
      nh_test_name = ''
    END IF
    !--------------------------------------------------------------------

    ! Vertical grid
    IF (ivctype==1) THEN
      IF (TRIM(vct_filename) == "") THEN
        IF (itype_laydistr /= 3 .AND. layer_thickness < 0) THEN
          CALL finish(routine, &
            & "ivctype=1 requires layer_thickness>0 or itype_laydistr=3 or vct_filename/=''")
        ENDIF
      ENDIF
    ENDIF



    !--------------------------------------------------------------------
    ! Testcases (nonhydrostatic)
    !--------------------------------------------------------------------
    IF (.NOT. ltestcase .AND. rayleigh_type == RAYLEIGH_CLASSIC) THEN
      CALL finish(routine, &
        & 'rayleigh_type = RAYLEIGH_CLASSIC not applicable to real case runs.')
    ENDIF

    IF ( ( nh_test_name=='APE_nwp'.OR. nh_test_name=='dcmip_tc_52' ) .AND.  &
      &  ( ANY(atm_phy_nwp_config(:)%inwp_surface > 0) ) ) THEN
      CALL finish(routine, &
        & 'surface scheme must be switched off, when running the APE test')
    ENDIF

    IF ( lmoist_thdyn ) THEN
      SELECT CASE (iforcing)
      CASE(IHELDSUAREZ,INOFORCING,ILDF_DRY)     
        CALL message(routine, &
           'lmoist_thdyn is reset to false .FALSE. because a dry model configuration is used')
        lmoist_thdyn = .FALSE.
      END SELECT
    ENDIF

    !--------------------------------------------------------------------
    ! SCM single column model
    !--------------------------------------------------------------------
    IF (l_scm_mode) THEN
      ! data read from 0: ASCII, 1: normal netcdf file, 2: DEPHY unified format
      IF ( (i_scm_netcdf /= 1) .AND. (i_scm_netcdf /= 2) ) &
        CALL finish(routine, 'i_scm_netcdf not valid for SCM, only 1 or 2 allowed')
      ltestcase      = .TRUE.
      is_plane_torus = .TRUE.
      CALL message( routine, 'l_scm_mode: ltestcase and is_plane_torus has been set to TRUE')

      !if time dependent surface BD conditions then use ls_forcing
      !routines for interpolation in time
      IF ( (scm_sfc_temp .GE. 1) .OR. (scm_sfc_qv .GE. 1) .OR. (scm_sfc_mom .GE. 1) ) THEN
#ifndef __NO_ICON_LES__
        is_ls_forcing = .TRUE.
        CALL message( routine, 'scm_sfc_... requires is_ls_forcing=TRUE. is_ls_forcing has been set to TRUE')
#else
        CALL finish( routine, 'scm_sfc_... requires is_ls_forcing=TRUE, but --disable-les has been set')
#endif
      END IF
    ELSE
      i_scm_netcdf   = 0
    END IF


    !--------------------------------------------------------------------
    ! Nonhydrostatic atm
    !--------------------------------------------------------------------

    IF (ldeepatmo) THEN
      IF (.NOT. ANY([inoforcing, inwp, iaes] == iforcing)) THEN
        CALL finish(routine, 'Deep-atmosphere configuration: incompatible iforcing')
      ELSEIF (ltestcase .AND. TRIM(nh_test_name) /= 'dcmip_bw_11') THEN
        CALL finish(routine, 'Deep-atmosphere configuration: the only supported testcase is "dcmip_bw_11"')
      ELSEIF (lplane .OR. is_plane_torus) THEN
        CALL finish(routine, 'Deep-atmosphere configuration is incompatible with plane or torus modes')
      ELSEIF (iadv_rhotheta /= 2) THEN
        CALL finish(routine, 'Deep-atmosphere configuration requires iadv_rhotheta = 2')
      ENDIF
    ENDIF ! IF (ldeepatmo)

    !--------------------------------------------------------------------
    ! Atmospheric physics, general
    !--------------------------------------------------------------------

#ifdef __NO_AES__
    IF ( iforcing==iaes ) &
      CALL finish( routine, 'AES physics desired, but compilation with --disable-aes' )
#endif

#ifdef __NO_NWP__
    IF ( iforcing==inwp ) &
      CALL finish( routine, 'NWP physics desired, but compilation with --disable-nwp' )
#endif

    !--------------------------------------------------------------------
    ! NWP physics
    !--------------------------------------------------------------------
    IF (iforcing==inwp) THEN

      DO jg =1,n_dom

        IF (atm_phy_nwp_config(jg)%inwp_gscp /= 8) THEN
          IF( atm_phy_nwp_config(jg)%inwp_satad == 0       .AND. &
          & ((atm_phy_nwp_config(jg)%inwp_convection >0 ) .OR. &
          &  (atm_phy_nwp_config(jg)%inwp_gscp > 0 )   ) ) &
          &  CALL finish( routine,'satad has to be switched on')
        ENDIF

        IF( (atm_phy_nwp_config(jg)%inwp_gscp==0) .AND. &
          & (atm_phy_nwp_config(jg)%inwp_convection==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_radiation==0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_sso==0)  .AND. &
          & (atm_phy_nwp_config(jg)%inwp_surface == 0) .AND.&
          & (atm_phy_nwp_config(jg)%inwp_turb> 0) )   &
        CALL message(routine,' WARNING! NWP forcing set but '//&
                    'only turbulence selected!')

        IF( atm_phy_nwp_config(jg)%inwp_surface == 1 .AND. &
        &   atm_phy_nwp_config(jg)%inwp_gscp == 0 ) &
        ! Perhaps it would be easy to implement this combination if needed.
        &  CALL finish( routine,'Surface model TERRA requires a gscp scheme at the moment.')


        IF (( atm_phy_nwp_config(jg)%inwp_turb == icosmo ) .AND. &
          & (turbdiff_config(jg)%lconst_z0) ) THEN
          CALL message(routine,' WARNING! NWP forcing set but '//  &
                      'idealized (horizontally homogeneous) roughness '//&
                      'length z0 selected!')
        ENDIF

        IF (.NOT. ltestcase .AND. atm_phy_nwp_config(jg)%inwp_surface == 0) THEN
          CALL finish( routine,'Real-data applications require using a surface scheme!')
        ENDIF

        IF ( (atm_phy_nwp_config(jg)%icpl_rad_reff == 2)  .AND.  & 
              (atm_phy_nwp_config(jg)%inwp_radiation/= 4) ) THEN
          CALL finish( routine, 'Wrong value for: icpl_rad_reff. Coupling effective radius for all '//  &  
                                 'hydrometeors only works with ECRAD!')
        ENDIF


        ! check radiation scheme in relation to chosen ozone and irad_aero=iRadAeroTegen to itopo

        IF ( (atm_phy_nwp_config(jg)%inwp_radiation > 0) )  THEN

          SELECT CASE (irad_o3)
          CASE (0) ! ok
            CALL message(routine,'radiation is used without ozone')
          CASE (2,4,5,6,7,9,11,79,97) ! ok
            CALL message(routine,'radiation is used with ozone')
          CASE (10) ! ok
            CALL message(routine,'radiation is used with ozone calculated from ART')
            IF ( .NOT. lart ) THEN
              CALL finish(routine,'irad_o3 currently is 10 but lart is false.')
            ENDIF
          CASE default
            CALL finish(routine,'irad_o3 currently has to be 0, 2, 4, 5, 6, 7, 9, 10, 11, 79 or 97.')
          END SELECT

          ! Tegen aerosol and itopo (Tegen aerosol data have to be read from external data file)
          IF ( ( irad_aero == iRadAeroTegen ) .AND. ( itopo /=1 ) ) THEN
            CALL finish(routine,'irad_aero=6 (Tegen) requires itopo=1')
          ENDIF

          IF ( .NOT. ANY ( irad_aero ==  (/iRadAeroTegen, iRadAeroCAMSclim , iRadAeroCAMStd ,  &
                         &  iRadAeroART, iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc,      &
                         &  iRadAeroKinneVolc, iRadAeroKinneVolcSP, iRadAeroKinneSP/) ) .AND.  &
            &  ( atm_phy_nwp_config(jg)%icpl_aero_gscp > 0 .OR. icpl_aero_conv > 0 ) ) THEN
            CALL finish(routine,'aerosol-precipitation coupling requires irad_aero=6,7,8,9,12,13,14,15,18 or 19')
          ENDIF

          ! reset lscale_cdnc to .false. if the SP scheme or icpl_aero_gscp = 3 (or both) are not set
          IF ( atm_phy_nwp_config(jg)%lscale_cdnc .AND. atm_phy_nwp_config(jg)%icpl_aero_gscp /= 3 ) THEN
            IF ( .NOT. ANY ( irad_aero ==  (/iRadAeroKinneVolcSP, iRadAeroKinneSP/) ) ) THEN
              atm_phy_nwp_config(jg)%lscale_cdnc = .false.
              CALL message(routine,'cdnc scaling is only effective in combination with the simple plumes &
                                   &(irad_aero=18,19) and icpl_aero_gscp = 3; reset lscale_cdnc to .false.')
            ENDIF
          ENDIF

          ! check if CAMS/Tegen aerosols are available for DeMott ice nucleation scheme
          IF (icpl_aero_ice == 1 .AND. .NOT. ANY(irad_aero == (/iRadAeroTegen, iRadAeroCAMSclim, iRadAeroCAMStd/) ) ) &
            & CALL finish(routine,'icpl_aero_ice = 1 requires irad_aero= 6,7 or 8')

#ifdef _OPENACC
          IF ( icpl_aero_ice == 1 ) THEN
            CALL finish(routine,'DeMott ice nucleation icpl_aero_ice > 0 is currently not supported on GPU.')
          END IF
#endif

          ! Kinne, CMIP6 volcanic aerosol only work with ecRad
          IF ( ANY( irad_aero == (/iRadAeroConstKinne,iRadAeroKinne,iRadAeroVolc,            &
            &                      iRadAeroKinneVolc,iRadAeroKinneVolcSP,iRadAeroKinneSP/) ) &
            &  .AND. atm_phy_nwp_config(jg)%inwp_radiation /= 4 ) THEN
            CALL finish(routine,'irad_aero = 12, 13, 14, 15, 18 or 19 requires inwp_radiation=4')
          ENDIF

          IF ( irad_aero == 5 ) THEN
            CALL finish(routine,'irad_aero=5 (Tanre climatology) has been removed')
          ENDIF

          ! Transient solar radiation only works with ecRad
          IF ( ANY( isolrad == (/2/) ) .AND. atm_phy_nwp_config(jg)%inwp_radiation /= 4 ) THEN
            CALL finish(routine,'isolrad = 2 requires inwp_radiation = 4')
          ENDIF

          ! ecRad specific checks
          IF ( (atm_phy_nwp_config(jg)%inwp_radiation == 4) )  THEN
            IF (.NOT. ANY( irad_h2o     == (/0,1/)         ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_h2o has to be 0 or 1')
            IF (.NOT. ANY( irad_co2     == (/0,2,4/)       ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_co2 has to be 0, 2 or 4')
            IF (.NOT. ANY( irad_ch4     == (/0,2,3,4/)     ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_ch4 has to be 0, 2, 3 or 4')
            IF (.NOT. ANY( irad_n2o     == (/0,2,3,4/)     ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_n2o has to be 0, 2, 3 or 4')
            IF (.NOT. ANY( irad_o3      == (/0,5,7,9,10,11,79,97/) ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_o3 has to be 0, 5, 7, 9, 10, 11, 79 or 97')
            IF (.NOT. ANY( irad_o2      == (/0,2/)         ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_o2 has to be 0 or 2')
            IF (.NOT. ANY( irad_cfc11   == (/0,2,4/)       ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_cfc11 has to be 0, 2 or 4')
            IF (.NOT. ANY( irad_cfc12   == (/0,2,4/)       ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, irad_cfc12 has to be 0, 2 or 4')
            IF (.NOT. ANY( irad_aero    == (/iRadAeroNone, iRadAeroConst, iRadAeroTegen, iRadAeroART, &
              &             iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc, iRadAeroCAMSclim,        &
              &             iRadAeroCAMStd, iRadAeroKinneVolc, iRadAeroKinneVolcSP, iRadAeroKinneSP/) ) ) THEN
              WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' is invalid for inwp_radiation=4'
              CALL finish(routine,message_text)
            ENDIF
            IF (.NOT. ANY( icld_overlap == (/1,2,5/)       ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, icld_overlap has to be 1, 2 or 5')
            IF ( ecrad_igas_model   == 1  .AND. .NOT. ecrad_use_general_cloud_optics) THEN
              ecrad_use_general_cloud_optics = .TRUE.
              CALL message(routine,'Warning: Reset ecrad_use_general_cloud_optics = T. &
                                    &It is the only valid option for ECCKD')
            ENDIF
            IF ( ecrad_use_general_cloud_optics )THEN
              IF (.NOT. ANY( ecrad_iliquid_scat == (/0/)   ) ) &
                &  CALL finish(routine,'For inwp_radiation = 4 ecrad_use_general_cloud_optics = T, &
                                       &ecrad_iliquid_scat has to be 0')
              IF (.NOT. ANY( ecrad_iice_scat    == (/0,10,11/) ) ) &
                &  CALL finish(routine,'For inwp_radiation = 4 ecrad_use_general_cloud_optics = T, &
                                       &ecrad_iice_scat has to be 0, 10 or 11')
            ELSE
              IF (.NOT. ANY( ecrad_iliquid_scat == (/0,1/)   ) ) &
                &  CALL finish(routine,'For inwp_radiation = 4 ecrad_use_general_cloud_optics = F, &
                                       &ecrad_iliquid_scat has to be 0 or 1')
              IF (.NOT. ANY( ecrad_iice_scat    == (/0,1,2/) ) ) &
                &  CALL finish(routine,'For inwp_radiation = 4 ecrad_use_general_cloud_optics = F, & 
                                       &ecrad_iice_scat has to be 0, 1 or 2')
            ENDIF
            IF (.NOT. ANY( ecrad_isnow_scat    == (/-1,0,10/) ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, ecrad_isnow_scat has to be -1, 0 or 10')
            IF (.NOT. ANY( ecrad_irain_scat    == (/-1,0/) ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, ecrad_isnow_scat has to be -1 or 0')
            IF (.NOT. ANY( ecrad_igraupel_scat    == (/-1,0,10/) ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, ecrad_igraupel_scat has to be -1, 0 or 10')
            IF (.NOT. ANY( ecrad_isolver  == (/0,1,2,3/)       ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, ecrad_isolver has to be 0, 1, 2 or 3')
            IF (.NOT. ANY( ecrad_igas_model   == (/0,1/)   ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, ecrad_igas_model has to be 0 or 1')
            IF (ecrad_igas_model == 1 .AND. .NOT. ANY(irad_aero == (/iRadAeroNone, iRadAeroConst, iRadAeroTegen/) ) )&
              &  CALL finish(routine,'For ecrad_igas_model=1, only Tegen aerosol implemented (irad_aero=0,2,6)')
            IF (ecrad_igas_model == 1 .AND. isolrad /= 1) THEN
              isolrad = 1
              CALL message(routine,'Warning: For ecrad_igas_model = 1, only Coddington scaling is available. Setting isolrad=1')
            ENDIF
            IF (.NOT. ANY( isolrad      == (/0,1,2/)       ) ) &
              &  CALL finish(routine,'For inwp_radiation = 4, isolrad has to be 0, 1 or 2')
            IF ( .NOT. ecrad_use_general_cloud_optics .AND. &
              & (ecrad_isnow_scat > -1 .OR. ecrad_igraupel_scat > -1 .OR. ecrad_irain_scat > -1 )) &
              &  CALL finish(routine,'Ecrad with qr/qs/qg can only be used with ecrad_use_general_cloud_optics = T')
            IF ( ecrad_use_general_cloud_optics .AND. & 
              & (ecrad_isnow_scat > -1 .OR. ecrad_igraupel_scat > -1 .OR. ecrad_irain_scat > -1 ) .AND. &
              &  atm_phy_nwp_config(jg)%icpl_rad_reff /= 2 ) &
              &  CALL finish(routine,'Ecrad with qr/qs/qg can only be used with icpl_rad_reff= 2')
            IF ( ecrad_igraupel_scat > -1 .AND. ANY(atm_phy_nwp_config(jg)%inwp_gscp == (/1,3/) ) ) THEN
              ecrad_igraupel_scat = -1
              CALL message(routine,'Warning: Reset ecrad_igraupel_scat = -1 because of no graupel in microphysics')
            ENDIF
          ELSE
            IF ( ecrad_llw_cloud_scat ) &
              &  CALL message(routine,'Warning: ecrad_llw_cloud_scat is set to .true., but ecRad is not used')
            IF ( ecrad_iliquid_scat /= 0 ) &
              &  CALL message(routine,'Warning: ecrad_iliquid_scat is explicitly set, but ecRad is not used')
            IF ( ecrad_iice_scat /= 0 ) &
              &  CALL message(routine,'Warning: ecrad_iice_scat is explicitly set, but ecRad is not used')
            IF ( ecrad_isolver /= 0 ) &
              &  CALL message(routine,'Warning: ecrad_isolver is explicitly set, but ecRad is not used')
            IF ( ecrad_igas_model /= 0 ) &
              &  CALL message(routine,'Warning: ecrad_igas_model is explicitly set, but ecRad is not used')
            IF ( ecrad_use_general_cloud_optics ) &
              &  CALL message(routine,'Warning: ecrad_use_general_cloud_optics is explicitly set, but ecRad is not used')
          ENDIF

        ELSE

          SELECT CASE (irad_o3)
          CASE(0) ! ok
          CASE default
            irad_o3 = 0
            CALL message(routine,'running without radiation => irad_o3 reset to 0')
          END SELECT

        ENDIF !inwp_radiation

        !! Checks for simple prognostic aerosol scheme
        IF (iprog_aero > 0) THEN
#ifndef __NO_ICON_LES__
          IF (atm_phy_nwp_config(jg)%is_les_phy) &
            & CALL finish(routine,'iprog_aero > 0 can not be combined with LES physics')
#endif
          IF (irad_aero /= iRadAeroTegen) &
            & CALL finish(routine,'iprog_aero > 0 currently only available for irad_aero=6 (Tegen)')
        ENDIF

        !! check microphysics scheme
        IF (   ANY(atm_phy_nwp_config(1:n_dom)%inwp_gscp == 2) .AND. &
             & ANY(atm_phy_nwp_config(1:n_dom)%inwp_gscp == 1) ) THEN
          CALL finish(routine,'combining inwp_gscp=1 and inwp_gscp=2 in nested runs is not allowed')
        END IF

        IF (  atm_phy_nwp_config(jg)%mu_rain < 0.0   .OR. &
          &   atm_phy_nwp_config(jg)%mu_rain > 5.0)  THEN
          CALL finish(routine,'mu_rain requires: 0 < mu_rain < 5')
        END IF

        IF (  atm_phy_nwp_config(jg)%mu_snow < 0.0   .OR. &
          &   atm_phy_nwp_config(jg)%mu_snow > 5.0)  THEN
          CALL finish(routine,'mu_snow requires: 0 < mu_snow < 5')
        END IF ! microphysics

        IF (  atm_phy_nwp_config(jg)%inwp_turb /= icosmo  .AND. &
          &  atm_phy_nwp_config(jg) % cfg_2mom % lturb_enhc ) THEN
          CALL finish(routine,' Turbulence enhancement of collisions '//  &
                      'in two-moment scheme (lturb_enhc) only applicable for inwp_turb = 1')
        ENDIF

        IF (  iforcing==iaes  .AND. &
          &  atm_phy_nwp_config(jg) % cfg_2mom % lturb_enhc ) THEN
          CALL finish(routine,' Turbulence enhancement of collisions '//  &
                      'in two-moment scheme (lturb_enhc) not applicable for aes physics.')
        ENDIF

        SELECT CASE (atm_phy_nwp_config(jg)%inwp_surface)
        CASE (0)
          IF (ntiles_lnd > 1) THEN
            ntiles_lnd = 1
            CALL message(routine,'Warning: ntiles reset to 1 because the surface scheme is turned off')
          ENDIF

        CASE (LSS_JSBACH)
          IF (ntiles_lnd > 1) THEN
            ntiles_lnd = 1
            CALL message(routine,'Warning: ntiles reset to 1 because JSBACH handles tiles internally')
          ENDIF
          IF (llake) THEN
            llake = .FALSE.
            CALL message(routine,'Warning: llake=.FALSE. because JSBACH handles lakes internally')
          END IF
        END SELECT

      ENDDO

#ifdef _OPENACC
    IF ( irad_aero == iRadAeroCAMSclim) THEN
        CALL finish(routine,'CAMS 3D climatology irad_aero=7 is currently not supported on GPU.')
    END IF
    IF ( irad_aero == iRadAeroCAMStd) THEN
        CALL finish(routine,'CAMS forecast irad_aero=8 is currently not supported on GPU.')
    END IF
    IF ( atm_phy_nwp_config(jg)%icpl_aero_gscp == 3 ) THEN
        CALL finish(routine,'Using cloud-droplet number climatology icpl_aero_gscp = 3 is currently not supported on GPU.')
    END IF
#endif

    END IF

    !--------------------------------------------------------------------
    ! Tracers and diabatic forcing
    !--------------------------------------------------------------------

    ! General
    IF ((itime_scheme_nh==tracer_only) .AND. (.NOT.ltransport)) THEN
      WRITE(message_text,'(A,i2,A)') &
        'nonhydrostatic_nml:itime_scheme set to ', tracer_only, &
        '(TRACER_ONLY), but ltransport to .FALSE.'
      CALL finish( routine,message_text)
    END IF



#ifdef _OPENACC
    IF (ltransport) THEN
      DO jg =1,n_dom
        IF ( .not. advection_config(jg)%llsq_svd ) THEN
          ! The .FALSE. version is not supported by OpenACC. However,
          ! both versions do the same. The .TRUE. version is faster
          ! during runtime on most platforms but involves a more
          ! expensive calculation of coefficients.
          CALL message(routine, 'WARNING: llsq_svd has been set to .true. for this OpenACC GPU run.')
          advection_config(jg)%llsq_svd = .TRUE.
        END IF
      ENDDO
    ENDIF
#endif


    !--------------------------------------------------------------------
    ! checking the meanings of the io settings
    !--------------------------------------------------------------------


    IF (lnetcdf_flt64_output) THEN
       CALL message(routine,'NetCDF output of floating point variables will be in 64-bit accuracy')
       IF (.NOT. use_dp_mpi2io) THEN
          use_dp_mpi2io = .TRUE.
          CALL message(routine,'--> use_dp_mpi2io is changed to .TRUE. to allow 64-bit accuracy in the NetCDF output.')
       END IF
    ELSE
       CALL message(routine,'NetCDF output of floating point variables will be in 32-bit accuracy')
    END IF

    IF (activate_sync_timers .AND. .NOT. ltimer) THEN
      activate_sync_timers = .FALSE.
      CALL message(routine, "namelist parameter 'activate_sync_timers' has &
        &been set to .FALSE., because global 'ltimer' flag is disabled.")
    END IF
    IF (timers_level > 9 .AND. .NOT. activate_sync_timers) THEN
      activate_sync_timers = .TRUE.
      CALL message(routine, "namelist parameter 'activate_sync_timers' has &
        &been set to .TRUE., because global 'timers_level' is > 9.")
    END IF

    DO jg =1,n_dom
      echotop_meta(jg)%nechotop = 0
      DO i=1, max_echotop
        IF (echotop_meta(jg)%dbzthresh(i) >= -900.0_wp) THEN
          echotop_meta(jg)%nechotop = echotop_meta(jg)%nechotop + 1
        END IF
      END DO
      IF ( is_variable_in_output_dom(var_name="echotop" , jg=jg) .AND. &
           echotop_meta(jg)%nechotop == 0 ) THEN
        WRITE (message_text, '(a,i2,a,i2.2,a)') 'output of "echotop" in ml_varlist on domain ', jg, &
             ' not possible due to invalid echotop_meta(', jg, ')%dbzthresh specification'
        CALL finish(routine, message_text)
      END IF
      IF ( is_variable_in_output_dom(var_name="echotopinm" , jg=jg) .AND. &
           echotop_meta(jg)%nechotop == 0 ) THEN
        WRITE (message_text, '(a,i2,a,i2.2,a)') 'output of "echotopinm" in ml_varlist on domain ', jg, &
             ' not possible due to invalid echotop_meta(', jg, ')%dbzthresh specification'
        CALL finish(routine, message_text)
      END IF
      IF (echotop_meta(jg)%time_interval < 0.0_wp) THEN
        WRITE (message_text, '(a,i2.2,a,f0.1,a)') 'invalid echotop_meta(', jg, &
             ')%time_interval = ', echotop_meta(jg)%time_interval, ' [seconds] given in namelist /io_nml/. Must be >= 0.0!'
        CALL finish(routine, message_text)
      END IF
    END DO

    n_wshear = 0
    DO i=1, max_wshear
      IF (wshear_uv_heights(i) > 0.0_wp) THEN
        n_wshear = n_wshear + 1
        ! shift valid values towards the start of the vector:
        wshear_uv_heights(n_wshear) = wshear_uv_heights(i)
      END IF
    END DO
    IF (n_wshear < max_wshear) wshear_uv_heights(n_wshear+1:) = -999.99_wp
    DO jg=1, n_dom
      IF ( ( is_variable_in_output_dom(var_name="wshear_u" , jg=jg) .OR. &
             is_variable_in_output_dom(var_name="wshear_v" , jg=jg) ) .AND. &
           n_wshear == 0 ) THEN
        message_text(:) = ' '
        WRITE (message_text, '(a)') 'output of "wshear_u" and/or "wshear_v" in ml_varlist'// &
             ' not possible because nml-parameter "wshear_uv_heights" (io_nml) contains no heights > 0.0!'
        CALL finish(routine, message_text)
      END IF
    END DO

    n_srh = 0
    DO i=1, max_srh
      IF (srh_heights(i) > 0.0_wp) THEN
        n_srh = n_srh + 1
        ! shift valid values towards the start of the vector:
        srh_heights(n_srh) = srh_heights(i)
      END IF
    END DO
    IF (n_srh < max_srh) srh_heights(n_srh+1:) = -999.99_wp
    DO jg=1, n_dom
      IF ( ( is_variable_in_output_dom(var_name="srh" , jg=jg) ) .AND. &
           n_srh == 0 ) THEN
        message_text(:) = ' '
        WRITE (message_text, '(a)') 'output of "srh" in ml_varlist'// &
             ' not possible because nml-parameter "srh_heights" (io_nml) contains no heights > 0.0!'
        CALL finish(routine, message_text)
      END IF
    END DO


    !--------------------------------------------------------------------
    ! Realcase runs
    !--------------------------------------------------------------------

    IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN  ! start from dwd analysis with incremental update

      ! check if the appropriate physics package has been selected
      IF (iforcing /= inwp) THEN
        CALL finish(routine,"(iterative) IAU is available for NWP physics only")
      ENDIF

      ! check analysis update window
      !
      IF ( (dt_iau > 0._wp) .AND. (dt_iau < dtime)) THEN
        ! If dt_iau is chosen to be larger than 0, it must be >= dtime at least.
        dt_iau = dtime
        WRITE (message_text,'(a,a,f6.2)') "Wrong value for dt_iau. ", &
          &   "If >0 then at least equal to advective/phys tstep ",dtime
        CALL finish('initicon_nml:', message_text)
      ENDIF

      IF (.NOT. isRestart()) THEN
        reference_dt => newDatetime("1980-06-01T00:00:00.000")

        msecs_restart   = getTotalMilliSecondsTimeDelta(time_config%tc_dt_restart, reference_dt)
        secs_restart    = 0.001_wp * REAL(msecs_restart,wp)
        secs_checkpoint = dt_checkpoint
        secs_iau_end    = dt_iau+timeshift%dt_shift
        ! ignore restart and/or checkpoint intervals if zero:
        IF (secs_restart    <= 0._wp)  secs_restart    = secs_iau_end
        IF (secs_checkpoint <= 0._wp)  secs_checkpoint = secs_restart
        IF (MIN(secs_checkpoint, secs_restart) < secs_iau_end) THEN
          CALL finish('atm_crosscheck:', "Restarting is not allowed within the IAU phase")
        ENDIF

        CALL deallocateDatetime(reference_dt)

        IF (l_limited_area) THEN
          ! For a negative IAU shift, no extra boundary file can be read. So it has to be taken
          ! from the first guess file.
          IF (timeshift%dt_shift < 0._wp .AND. .NOT. latbc_config%init_latbc_from_fg) THEN
            CALL finish('atm_crosscheck:', "For dt_shift<0, latbc has &
              &to be taken from first guess (init_latbc_from_fg)")
          ENDIF
        ENDIF
      ENDIF

      DO jg = 2, n_dom
        IF (start_time(jg) > timeshift%dt_shift .AND. start_time(jg) < dt_iau+timeshift%dt_shift) THEN
          CALL finish('atm_crosscheck:', "Starting a nest is not allowed within the IAU phase")
        ENDIF
      ENDDO

      ! IAU modes MODE_IAU_OLD cannot be combined with snowtiles
      ! when performing snowtile warmstart.
      IF ((ntiles_lnd > 1) .AND. (.NOT. ltile_coldstart) .AND. (lsnowtile)) THEN
        IF ( init_mode == MODE_IAU_OLD ) THEN
          WRITE (message_text,'(a,i2)') "lsnowtile=.TRUE. not allowed for IAU-Mode ", init_mode
          CALL finish(routine, message_text)
        ENDIF
      ENDIF

    ENDIF

    IF (itune_gust_diag == 3 .AND. ntiles_lnd == 1) THEN
      WRITE (message_text,'(a)') "itune_gust_diag = 3 requires ntiles > 1"
      CALL finish(routine, message_text)
    ENDIF

    ! check meteogram configuration
    IF (ANY(meteogram_output_config(:)%lenabled) .AND. .NOT. output_mode%l_nml) THEN
      CALL finish(routine, "Meteograms work only for run_nml::output='nml'!")
    END IF
    CALL check_meteogram_configuration(num_io_procs)

    IF (ANY(iforcing == [iaes, inwp])) CALL land_crosscheck()

    IF (iforcing==inwp) CALL coupled_crosscheck()

    CALL art_crosscheck()

    CALL check_nudging( n_dom, iforcing, ivctype, top_height,                            &
      &                 l_limited_area, num_prefetch_proc, latbc_config%lsparse_latbc,   &
      &                 latbc_config%itype_latbc, latbc_config%nudge_hydro_pres,         &
      &                 latbc_config%latbc_varnames_map_file, LATBC_TYPE_CONST,          &
      &                 LATBC_TYPE_EXT, is_plane_torus, lart, ltransport  )

    CALL check_upatmo( n_dom_start, n_dom, iforcing, atm_phy_nwp_config(:)%lupatmo_phy,   &
      &                l_limited_area, ivctype, flat_height, itype_vert_expol, init_mode, &
      &                atm_phy_nwp_config(:)%inwp_turb, atm_phy_nwp_config(:)%inwp_radiation)

    !--------------------------------------------------------------------
    ! assimiliation
    !--------------------------------------------------------------------
    IF ( MODE_IAU == init_mode ) THEN
      IF( ANY(assimilation_config(:)%dace_coupling) .AND. .NOT. iterate_iau ) THEN 
        ! The MEC in dace_coupling needs the fully initialized state to compute the
        ! model equivalents for vv=0.
        CALL finish(routine,                                        &
        &  'assimilation: dace_coupling requires iterate_iau')
        ! One might loosen this check slightly by allowing non-iterative IAU if
        ! dace_time_ctrl(0) > 0. However this case should be tested before
        ! modifying this crosscheck.
      ENDIF
    ENDIF

    ! ********************************************************************************
    !
    ! Cross checks for EMVORADO-related namelist parameters:
    ! *  EMVORADO may be switched on only if nonhydro and inwp enabled.
    ! *  EMVORADO may be switched on only if preprocessor flag enabled (see below).
    !
    ! ********************************************************************************

    CALL emvorado_crosscheck()



    ! ********************************************************************************
    !
    !  Cross checks for SPPT (Stochastic Perturbation of Physics Tendencies)
    !
    ! ********************************************************************************

    IF( ANY(sppt_config(1:n_dom)%lsppt) ) THEN
      CALL crosscheck_sppt()
    ENDIF

  END  SUBROUTINE atm_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE land_crosscheck

    INTEGER  :: jg
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::land_crosscheck'

#ifdef __NO_JSBACH__
    IF (ANY(aes_phy_config(:)%ljsb)) THEN
      CALL finish(routine, "This version was compiled without jsbach. Compile with __JSBACH__, or set ljsb=.FALSE.")
    ENDIF

    IF (ANY(atm_phy_nwp_config(1:n_dom)%inwp_surface == LSS_JSBACH)) &
        & CALL finish(routine, "This version was compiled without jsbach. Compile with __JSBACH__, or set inwp_surface to a different value.")
#else
    DO jg=1,n_dom
      IF (.NOT.aes_phy_config(jg)%ljsb) THEN
         IF (aes_phy_config(jg)%llake) THEN
            CALL message(routine, 'Setting llake = .FALSE. since ljsb = .FALSE.')
            aes_phy_config(jg)%llake = .FALSE.
         END IF
      ELSE
        IF (aes_vdf_config(jg)%use_tmx) THEN
          CALL message(routine, 'Setting llake = .FALSE. since using tmx (lakes are treated inside JSBACH)')
          aes_phy_config(jg)%llake = .FALSE.
        END IF
      END IF
    END DO
#endif

  END SUBROUTINE land_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE coupled_crosscheck

    CHARACTER(len=*), PARAMETER :: routine =  modname//'::coupled_crosscheck'

    IF ( ntiles_lnd == 1 .AND. ( is_coupled_to_ocean() .OR. is_coupled_to_hydrodisc() ) .AND. .NOT. &
        & (iforcing == inwp .AND. ALL(atm_phy_nwp_config(1:n_dom)%inwp_turb == ivdiff)) ) THEN
      CALL finish(routine, "Coupled atm/hydrodisc/ocean runs not supported with ntiles=1 when not using VDIFF")
    ENDIF

    IF ( sstice_mode /= 1 .AND. is_coupled_to_ocean() ) THEN
      CALL finish(routine, "Coupled atm/ocean runs only supported with sstice_mode=1 named SSTICE_ANA")
    ENDIF

    IF ( is_coupled_to_waves() .AND. (.NOT. iforcing == inwp) ) THEN
      CALL finish(routine, "Coupled atm/wave runs only supported with NWP physics (iforcing=3)")
    ENDIF

    IF ( is_coupled_to_waves() .AND. (ntiles_lnd == 1) ) THEN
      CALL finish(routine, "Coupled atm/wave runs require ntiles_lnd>1.")
    ENDIF

#ifdef _OPENACC
    IF ( is_coupled_to_hydrodisc() ) THEN
      CALL finish(routine, "Coupled atm/hydrodisc/ocean runs are not available on GPU")
    END IF

    IF ( is_coupled_to_waves() ) THEN
      CALL finish(routine, "Coupled atm/wave runs are not available on GPU")
    END IF

    IF ( is_coupled_to_ocean() ) THEN
      CALL finish(routine, "Coupled atm/ocean runs are not available on GPU")
    END IF
#endif

  END SUBROUTINE coupled_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE art_crosscheck

    CHARACTER(len=*), PARAMETER :: routine =  modname//'::art_crosscheck'
#ifdef __ICON_ART
    INTEGER  :: &
      &  jg
#endif

#ifndef __ICON_ART
    IF (lart) THEN
        CALL finish( routine,'run_nml: lart is set .TRUE. but ICON was compiled without -D__ICON_ART')
    ENDIF
#endif

    IF (.NOT. lart .AND. irad_aero == iRadAeroART ) THEN
      CALL finish(routine,'irad_aero=9 (ART) needs lart = .TRUE.')
    END IF

    IF ( ( irad_aero == iRadAeroART ) .AND. ( iprog_aero /= 0 ) ) THEN
      CALL finish(routine,'irad_aero=9 (ART) requires iprog_aero=0')
    ENDIF

#ifdef __ICON_ART
    IF ( ( irad_aero == iRadAeroART ) .AND. ( itopo /=1 ) ) THEN
      CALL finish(routine,'irad_aero=9 (ART) requires itopo=1')
    ENDIF

    DO jg= 1,n_dom
      IF(lredgrid_phys(jg) .AND. irad_aero == iRadAeroART .AND. atm_phy_nwp_config(jg)%inwp_radiation /= 4) THEN
        CALL finish(routine,'irad_aero=9 (ART) and a reduced radiation grid only works with ecRad (inwp_radiation=4)')
      ENDIF
      IF(art_config(jg)%iart_ari == 0 .AND. irad_aero == iRadAeroART) THEN
        CALL finish(routine,'irad_aero=9 (ART) needs iart_ari > 0')
      ENDIF
      IF(art_config(jg)%iart_ari > 0  .AND. irad_aero /= iRadAeroART) THEN
        CALL finish(routine,'iart_ari > 0 requires irad_aero=9 (ART)')
      ENDIF
    ENDDO

#ifdef _OPENACC
    DO jg = 1, n_dom
      IF (    art_config(jg)%iart_aci_cold  >  0  .OR.  &
          &   art_config(jg)%iart_aci_warm  >  0  .OR.  &
          &   art_config(jg)%iart_ari       >  0  .OR.  &      
          &   art_config(jg)%iart_dust      >  0  .OR.  &
          &   art_config(jg)%iart_init_aero >  0  .OR.  &
          &   art_config(jg)%iart_radioact  >  0  .OR.  &
          &   art_config(jg)%iart_seasalt   >  0  .OR.  &
          &   art_config(jg)%iart_volcano   >  0  ) THEN
        CALL finish(routine,  &
          &  'mo_nml_crosscheck: art_crosscheck: some activated art switches are currently not supported on GPU.')
      ELSEIF (    art_config(jg)%lart_chem            .OR.  &
              &   art_config(jg)%lart_chemtracer            ) THEN 
        CALL message(routine, 'WARNING: The switches lartchem and lart_chemtracer are not supported on GPU. Use them at your own risk. However, using this with OEM-specific cases is safe at the moment.') 
      END IF
    ENDDO
#endif

#endif
  END SUBROUTINE art_crosscheck
  !---------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------
  SUBROUTINE emvorado_crosscheck
    CHARACTER(len=*), PARAMETER :: routine =  'mo_nml_crosscheck:emvorado_crosscheck'

    INTEGER  :: jg, ndoms_radaractive
    CHARACTER(len=255) :: errstring

#ifndef HAVE_RADARFWO
    IF ( ANY(luse_radarfwo) ) THEN
        CALL finish( routine,'run_nml: luse_radarfwo is set .TRUE. in some domains but ICON was compiled without -DHAVE_RADARFWO')
    ENDIF
#endif

    ndoms_radaractive = 0
    DO jg = 1, n_dom
      IF (luse_radarfwo(jg)) ndoms_radaractive = ndoms_radaractive + 1
    END DO
#ifdef HAVE_RADARFWO
    IF ( ndoms_radaractive > ndoms_max_radar ) THEN
      errstring(:) = ' '
      WRITE (errstring, '(a,i3,a,i2,a)') 'luse_radarfwo is enabled (.true.) for ', ndoms_radaractive, &
           ' ICON domains, but EMVORADO supports max. ', ndoms_max_radar, &
           ' domains. You may increase  parameter ''ndoms_max'' in radar_data.f90.'
      CALL finish(routine, 'run_nml: '//TRIM(errstring))
    END IF
#endif

    IF ( num_io_procs_radar > 0 .AND. .NOT.ANY(luse_radarfwo) ) THEN
      CALL message(routine, 'Setting num_io_procs_radar = 0 because luse_radarfwo(:) = .FALSE.')
      num_io_procs_radar = 0
    END IF

    IF ( iforcing /= INWP .AND. ANY(luse_radarfwo) ) THEN
      errstring(:) = ' '
      WRITE (errstring, '(a,i2)') 'luse_radarfwo = .true. is only possible for NWP physics iforcing = ', INWP
      CALL finish(routine, 'run_nml: '//TRIM(errstring))
    END IF

  END SUBROUTINE emvorado_crosscheck

END MODULE mo_nml_crosscheck
