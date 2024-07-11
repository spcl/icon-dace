! Namelist for NWP physics
!
! these Subroutines are called by control model and construct the
! physics composition
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

MODULE mo_nwp_phy_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: max_dom, ivdiff, LSS_JSBACH
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,      ONLY: use_restart_namelists

  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,    &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config,                        &
    &                               config_lrtm_filename   => lrtm_filename,   &
    &                               config_cldopt_filename => cldopt_filename, &
    &                               config_icpl_aero_conv  => icpl_aero_conv,  &
    &                               config_iprog_aero      => iprog_aero,      &
    &                               config_icpl_o3_tp      => icpl_o3_tp,      &
    &                               config_icpl_aero_ice   => icpl_aero_ice,   &
    &                               config_lcuda_graph_turb_tran => lcuda_graph_turb_tran

  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_cuparameters,        ONLY: icapdcycl

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_nwp_phy_namelist

   !
   ! user defined calling intervals
   !
  REAL(wp) :: dt_conv(max_dom)   !> field element for convection
  REAL(wp) :: dt_ccov(max_dom)   !> field element for cloud cover
  REAL(wp) :: dt_rad(max_dom)    !! "-"                     radiation
  REAL(wp) :: dt_sso(max_dom)    !! "-"  for subscale orographic gravity waves
  REAL(wp) :: dt_gwd(max_dom)    !! "-"  for subscale gravity waves

  ! switches defining physics packages
  INTEGER  :: inwp_convection(max_dom)    !! convection
  LOGICAL  :: lshallowconv_only(max_dom)  !! use shallow convection only
  LOGICAL  :: lstoch_expl(max_dom)        !! use explicit stochastic shallow convection
  !!! NOTE! Explicit stochastic scheme is experimental, cloud ensemble variables not saved for restart!!!
  !!! --> Will fail restart test. Use SDE version if restart capability is required.
  LOGICAL  :: lstoch_sde(max_dom)         !! use stochastic differential eqns for convection
  LOGICAL  :: lstoch_deep(max_dom)        !! use stochastic deep convection parameterization
  LOGICAL  :: lvvcouple(max_dom)          !! use vertical velocity at 650hPa as criterion to couple
                                          !! shallow convection with resolved deep convection 
  LOGICAL  :: lvv_shallow_deep(max_dom)   !! use vertical velocity at 650hPa to distinguish between shallow and 
                                          !! deep convection within convection routines (instead of cloud depth)
  LOGICAL  :: lstoch_spinup(max_dom)      !! spin up cloud ensemble to equilibrium, in shallow stochastic convection
  LOGICAL  :: lrestune_off(max_dom)       !! switch off all resolution-dependent tuning in convection setup
  LOGICAL  :: lmflimiter_off(max_dom)     !! switch off MF limiters in convection setup
  INTEGER  :: nclds(max_dom)              !! max number of clouds in stochastic cloud ensemble
  LOGICAL  :: lgrayzone_deepconv(max_dom) !! use grayzone tuning for deep convection
  LOGICAL  :: ldetrain_conv_prec(max_dom) !! detrain convective rain and snow
  INTEGER  :: inwp_cldcover(max_dom)      !! cloud cover
  LOGICAL  :: lsgs_cond(max_dom)          !! subgrid-scale condensation related to cloud cover
  INTEGER  :: inwp_radiation(max_dom)     !! radiation
  INTEGER  :: inwp_sso(max_dom)           !! sso
  INTEGER  :: inwp_gwd(max_dom)           !! non-orographic gravity wave drag
  INTEGER  :: inwp_gscp(max_dom)          !! microphysics
  INTEGER  :: inwp_satad(max_dom)         !! saturation adjustment
  INTEGER  :: inwp_turb(max_dom)          !! turbulence
  INTEGER  :: inwp_surface(max_dom)       !! surface including soil, ocean, ice,lake

  INTEGER  :: itype_z0           !! type of roughness length data
  INTEGER  :: icpl_aero_gscp     !! type of aerosol-microphysics coupling
  LOGICAL  :: lscale_cdnc        !! switch to activate the scaling of external CDNCs
  INTEGER  :: icpl_aero_ice      !! type of aerosol-ice nucleation coupling
  INTEGER  :: icpl_aero_conv     !! type of coupling between aerosols and convection scheme
  INTEGER  :: iprog_aero         !! type of prognostic aerosol
  INTEGER  :: icpl_o3_tp         !! type of ozone-tropopause coupling
  REAL(wp) :: qi0, qc0           !! variables for hydci_pp
  REAL(wp) :: ustart_raylfric    !! velocity at which extra Rayleigh friction starts
  REAL(wp) :: efdt_min_raylfric  !! e-folding time corresponding to maximum relaxation coefficient
  LOGICAL  :: latm_above_top(max_dom) !! use extra layer above model top for radiation (reduced grid only)
  LOGICAL  :: lupatmo_phy(max_dom)    !! switch on/off upper-atmosphere physics in domains
  ! parameter for cloud microphysics
  real(wp) :: mu_rain            !! shape parameter in gamma distribution for rain
  real(wp) :: rain_n0_factor     !! tuning factor for intercept parameter of raindrop size distribution
  real(wp) :: mu_snow            !! ...for snow
  LOGICAL  :: lsbm_warm_full     !! false: Piggy Backing with 2M, true: full warm-phase SBM

  INTEGER  :: icalc_reff(max_dom)    !! type of effective radius calculation
  INTEGER  :: icpl_rad_reff(max_dom) !! coupling radiation and effective radius
  INTEGER  :: ithermo_water(max_dom) !! thermodynamic of water

  LOGICAL  :: lcuda_graph_turb_tran  !! Activate CUDA GRAPH in turbulent transfer

  !> NetCDF file containing longwave absorption coefficients and other data
  !> for RRTMG_LW k-distribution model ('rrtmg_lw.nc')
  CHARACTER(LEN=filename_max) :: lrtm_filename

  !> NetCDF file with RRTM Cloud Optical Properties for ECHAM6
  CHARACTER(LEN=filename_max) :: cldopt_filename

  NAMELIST /nwp_phy_nml/ inwp_convection, inwp_cldcover, lsgs_cond,  &
    &                    inwp_radiation, inwp_sso, inwp_gwd,         &
    &                    inwp_gscp, inwp_satad,                      &
    &                    inwp_turb, inwp_surface,                    &
    &                    dt_conv, dt_rad, dt_sso, dt_gwd, dt_ccov,   &
    &                    qi0, qc0, icpl_aero_gscp, icpl_aero_ice,    &
    &                    ustart_raylfric, efdt_min_raylfric,         &
    &                    latm_above_top, itype_z0, mu_rain,          &
    &                    mu_snow, icapdcycl, icpl_aero_conv,         &
    &                    lrtm_filename, cldopt_filename, icpl_o3_tp, &
    &                    iprog_aero, lshallowconv_only,lstoch_expl,  &
    &                    lvvcouple, lvv_shallow_deep,lstoch_spinup,  &
    &                    lrestune_off,nclds,                         &
    &                    lmflimiter_off, lstoch_sde,lstoch_deep,     &
    &                    ldetrain_conv_prec, rain_n0_factor,         &
    &                    icalc_reff, lupatmo_phy, icpl_rad_reff,     &
    &                    lgrayzone_deepconv, ithermo_water,          &
    &                    lsbm_warm_full, lcuda_graph_turb_tran,      &
    &                    lscale_cdnc

CONTAINS

  !-------------------------------------------------------------------------
  !
  !! Read Namelist for NWP physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - performs sanity checks
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  SUBROUTINE read_nwp_phy_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, jg
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = 'mo_nwp_phy_nml:read_nwp_phy_namelist'
    INTEGER  :: param_def
    REAL(wp) :: dt_conv_def, dt_rad_def, dt_sso_def, dt_gwd_def
    INTEGER  :: icalc_reff_def, icpl_rad_reff_def, ithermo_water_def

    !-----------------------
    ! 1a. default settings for domain-specific parmeters; will be rest to dummy values afterwards
    !     to determine what is actually specified in the namelist
    !-----------------------
    param_def   = 1  ! default for all parameterizations is currently '1'
    dt_conv_def = 600._wp
    dt_rad_def  = 1800._wp
    dt_sso_def  = 1200._wp
    dt_gwd_def  = 1200._wp

    icalc_reff_def = 0    ! Default is no calculation of effectives radius
    icpl_rad_reff_def = 0 ! Default is no coupling of effective radius and radiation
    ithermo_water_def = 0 ! Default is latent heat as a function of temperature in saturation adjustment 
                          ! but constant in microphysics. 

    inwp_gscp(:)       = param_def
    inwp_satad(:)      = param_def
    inwp_convection(:) = param_def
    inwp_radiation(:)  = param_def
    inwp_sso(:)        = param_def
    inwp_gwd(:)        = param_def
    inwp_cldcover(:)   = param_def
    inwp_turb(:)       = param_def
    inwp_surface(:)    = param_def

    dt_conv (:) = dt_conv_def
    dt_ccov (:) = dt_conv_def 
    dt_rad  (:) = dt_rad_def
    dt_sso  (:) = dt_sso_def
    dt_gwd  (:) = dt_gwd_def

    !------------------------------------------------------------------
    ! 1b. Defaults for some other variables
    !------------------------------------------------------------------

    nclds(:)              = 5000
    lshallowconv_only(:)  = .FALSE.
    lstoch_expl(:)        = .FALSE. ! default: no stochastic convection
    lstoch_sde(:)         = .FALSE. ! default: no stochastic convection with differential equations
    lstoch_deep(:)        = .FALSE. ! default: no stochastic deep convection is used
    lvvcouple(:)          = .FALSE. ! default: no vertical velocity used to couple shallow and resolve deep convection
    lvv_shallow_deep(:)   = .FALSE. ! default: shallow/deep distinguished by cloud depth (rdepths), not by vertical velocity
    lstoch_spinup(:)      = .FALSE. ! default: no spinup of stochastic shallow cloud ensemble
    lrestune_off(:)       = .FALSE. ! default: all tunings as for default master branch
    lmflimiter_off(:)     = .FALSE. ! default: mass flux limiters on
    lgrayzone_deepconv(:) = .FALSE.
    ldetrain_conv_prec(:) = .FALSE.
    lsgs_cond(:)          = .TRUE.  ! activate subgrid-scale condensation in cloud cover scheme


    lrtm_filename   = 'rrtmg_lw.nc'  
    cldopt_filename = 'ECHAM6_CldOptProps.nc'

    itype_z0 = 2  !  2 = land-cover related roughness lenght only (i.e. no orographic contrib)
    qi0      = 0.0_wp 
    qc0      = 0.0_wp 

    ! shape parameter for gamma distribution for rain and snow
    mu_rain = 0.0_wp
    mu_snow = 0.0_wp
    rain_n0_factor = 1.0_wp

    lsbm_warm_full = .TRUE. ! false: Piggy Backing with 2M, true: full warm-phase SBM

    ustart_raylfric    = 160._wp
    efdt_min_raylfric  = 10800._wp

    latm_above_top(:)  = .FALSE.  ! no extra layer above model top for radiation computation

    lupatmo_phy(:)     = .TRUE.   ! switch on upper-atmosphere physics
    lupatmo_phy(1)     = .FALSE.  ! switch off upper-atmosphere physics on dom 1 (please, do not touch this)

    ! CAPE correction to improve diurnal cycle of convection (moved from mo_cuparameters)
    icapdcycl = 0  ! 0= no CAPE diurnal cycle correction (IFS default prior to cy40r1, i.e. 2013-11-19)
                   ! 1=    CAPE - surface buoyancy flux (intermediate testing option)
                   ! 2=    CAPE - subcloud CAPE (IFS default starting with cy40r1)
                   ! 3=    Apply CAPE modification of (2) over land only, with additional restriction to the tropics

    ! coupling between aersols and cloud microphysics
    icpl_aero_gscp = 0  ! 0 = none
                        ! 1 = simple coupling with aerosol climatology disregarding the dependency of aerosol activation on vertical wind speed
                        ! 2 = more accurate coupling with aerosol climatology as a function of vertical wind speed
                        ! 3 = like 1 but using the cdnc from external parameter 

    ! scaling of external CDNCs (only for icpl_aero_gscp = 3), mainly for climate projections
    lscale_cdnc = .FALSE.  ! FALSE - no scaling
                           ! TRUE  - scaling according to the temporal evolution of the simple plumes

    ! coupling between aersols and ice nucleation
    icpl_aero_ice = 0   ! 0 = Cooper (1986)
                        ! 1 = DeMott (2015) with CAMS/Tegen dust concentration

    ! coupling between aersols and convection scheme
    icpl_aero_conv = 0  ! 0 = none
                        ! 1 = specify thresholds (QC and cloud thickness) for precip initiation depending on aerosol climatology instead of land-sea mask

    ! type of prognostic aerosol
    iprog_aero = 0  ! 0 = pure climatology
                    ! 1 = very simple prognostic scheme based on advection of and relaxation towards climatology for mineral dust
                    ! 2 = very simple prognostic scheme based on advection of and relaxation towards climatology for all species
                    ! 3 = very simple prognostic scheme based on advection of and relaxation towards climatology for all species including wildfires

    ! coupling between ozone and the tropopause
    icpl_o3_tp = 1      ! 0 = none
                        ! 1 = take climatological values from 100/350 hPa above/below the tropopause in the extratropics

    
    ! Calculation of effective radius
    icalc_reff(:)   =  icalc_reff_def ! 0      = no calculation (current default)
                       ! 1,2,4,5,6,7 = corresponding to the microphysics scheme 1,....,7 (same terminology as inwp_gscp)
                       ! 11     = corresponding to the microphysics scheme  inwp_gscp(jg)
                       ! 12     = RRTM parameterization

    icpl_rad_reff(:)=  icpl_rad_reff_def ! 0    = no coupling (using old RRTM Parameterization)
                       ! 1      = coupling RRTM/ECRAD with the effective radius parameterization (through two phases ice and frozen)
                       ! 2      = coupling ECRAD with the effective radius parameterization (each phase alone)

    ithermo_water(:)=  ithermo_water_def ! 0   = Latent heats (LH) constant in microphysics
                                         ! 1   = LH as function of temperature in microphysics

    lcuda_graph_turb_tran = .FALSE.   ! cuda graph deactivated by default

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, nwp_phy_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('nwp_phy_nml')
      READ(funit,NML=nwp_phy_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('nwp_phy_nml', status=istat)

    SELECT CASE (istat)
    CASE (POSITIONED)

      ! Set array parameters to dummy values to determine which ones are actively set in the namelist
      inwp_gscp(:)       = -1
      inwp_satad(:)      = -1
      inwp_convection(:) = -1
      inwp_radiation(:)  = -1
      inwp_sso(:)        = -1
      inwp_gwd(:)        = -1
      inwp_cldcover(:)   = -1
      inwp_turb(:)       = -1
      inwp_surface(:)    = -1

      dt_conv (:) = -999._wp
      dt_ccov (:) = -999._wp
      dt_rad  (:) = -999._wp
      dt_sso  (:) = -999._wp
      dt_gwd  (:) = -999._wp

      icalc_reff(:)      = -1
      icpl_rad_reff(:)   = -1
      ithermo_water(:)   = -1

      READ (nnml, nwp_phy_nml)   ! overwrite default settings

      ! Restore default values for global domain where nothing at all has been specified

      ! Physics packages
      IF (inwp_gscp(1)       < 0) inwp_gscp(1)       = param_def  !> 1 = hydci (COSMO-EU microphysics)
      IF (inwp_satad(1)      < 0) inwp_satad(1)      = param_def  !> 1 = saturation adjustment on
      IF (inwp_convection(1) < 0) inwp_convection(1) = param_def  !> 1 = Tiedtke/Bechthold convection
      IF (inwp_radiation(1)  < 0) inwp_radiation(1)  = param_def  !> 1 = RRTM radiation
      IF (inwp_sso(1)        < 0) inwp_sso(1)        = param_def  !> 1 = Lott and Miller scheme (COSMO)
      IF (inwp_gwd(1)        < 0) inwp_gwd(1)        = param_def  !> 1 = Orr-Ern-Bechthold scheme (IFS)
      IF (inwp_cldcover(1)   < 0) inwp_cldcover(1)   = param_def  !> 1 = diagnostic cloud cover (by Martin Koehler)
      IF (inwp_turb(1)       < 0) inwp_turb(1)       = param_def  !> 1 = turbdiff (COSMO diffusion and transfer)
      IF (inwp_surface(1)    < 0) inwp_surface(1)    = param_def  !> 1 = TERRA

      ! Time steps
      IF (dt_conv (1) < 0._wp) dt_conv (1) = dt_conv_def
      IF (dt_ccov (1) < 0._wp) dt_ccov (1) = dt_conv (1) ! this sets the previous default of dt_ccov=dt_conv
      IF (dt_sso  (1) < 0._wp) dt_sso  (1) = dt_sso_def
      IF (dt_gwd  (1) < 0._wp) dt_gwd  (1) = dt_gwd_def
      IF (dt_rad  (1) < 0._wp) dt_rad  (1) = dt_rad_def

      ! Extra calculation
      IF (icalc_reff(1)      < 0) icalc_reff(1)      = icalc_reff_def    ! Default no calculation of effective radius
      IF (icpl_rad_reff(1)   < 0) icpl_rad_reff(1)   = icpl_rad_reff_def ! Default no coupling of radiation with effective radius
      IF (ithermo_water(1)   < 0) ithermo_water(1)   = ithermo_water_def ! Default is constant LH in micro. 
      
      ! Copy values of parent domain (in case of linear nesting) to nested domains where nothing has been specified

      DO jg = 2, max_dom

        ! Physics packages
        IF (inwp_gscp(jg)       < 0) inwp_gscp(jg)       = inwp_gscp(jg-1)
        IF (inwp_satad(jg)      < 0) inwp_satad(jg)      = inwp_satad(jg-1)
        IF (inwp_convection(jg) < 0) inwp_convection(jg) = inwp_convection(jg-1)
        IF (inwp_radiation(jg)  < 0) inwp_radiation(jg)  = inwp_radiation(jg-1)
        IF (inwp_sso(jg)        < 0) inwp_sso(jg)        = inwp_sso(jg-1)
        IF (inwp_gwd(jg)        < 0) inwp_gwd(jg)        = inwp_gwd(jg-1)
        IF (inwp_cldcover(jg)   < 0) inwp_cldcover(jg)   = inwp_cldcover(jg-1)
        IF (inwp_turb(jg)       < 0) inwp_turb(jg)       = inwp_turb(jg-1)
        IF (inwp_surface(jg)    < 0) inwp_surface(jg)    = inwp_surface(jg-1)

        ! Time steps
        IF (dt_conv (jg) < 0._wp) dt_conv (jg) = dt_conv (jg-1) 
        IF (dt_ccov (jg) < 0._wp) dt_ccov (jg) = dt_conv (jg) ! this sets the previous default of dt_ccov=dt_conv
        IF (dt_sso  (jg) < 0._wp) dt_sso  (jg) = dt_sso  (jg-1)
        IF (dt_gwd  (jg) < 0._wp) dt_gwd  (jg) = dt_gwd  (jg-1)
        IF (dt_rad  (jg) < 0._wp) dt_rad  (jg) = dt_rad  (jg-1)
        
        ! Extra calculations
        IF (icalc_reff(jg)      < 0) icalc_reff(jg)       = icalc_reff(jg-1)
        IF (icpl_rad_reff(jg)   < 0) icpl_rad_reff(jg)    = icpl_rad_reff(jg-1)
        IF (ithermo_water(jg)   < 0) ithermo_water(jg)    = ithermo_water(jg-1)


        ! Upper-atmosphere physics
        IF (lupatmo_phy(jg)) lupatmo_phy(jg) = lupatmo_phy(jg-1)

      ENDDO


      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, nwp_phy_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    
    ! check for valid parameters in namelists

    DO jg = 1, max_dom

      IF ( ALL((/0,1,2,3,4,5,6,7,8,9/) /= inwp_gscp(jg)) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting for inwp_gscp. Must be 0,1,2,3,4,5,6,7,8 or 9.')
      END IF
      
      IF ( ALL((/0,1,2,4,5,6,7,8,100,101/) /= icalc_reff(jg)) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting for icalc_reff. Must be 0,1,2,4,5, 6,7,8, 100 or 101.')
      END IF

      IF ( ALL((/0,1,2/) /= icpl_rad_reff(jg)) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting for icpl_rad_reff. Must be 0,1,2')
      END IF
      
      IF ( (icpl_rad_reff(jg) > 0) .AND. (icalc_reff(jg) == 0) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting for icpl_rad_reff. It must be 0 if no reff is defined (icalc_reff =0)')
      END IF

      IF ( ALL((/0,1/) /= ithermo_water(jg)) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting for ithermo_water. Must be 0 or 1')
      END IF

      IF ( (ithermo_water(jg) == 1) .AND. ALL((/1,2,4,5,7/) /= inwp_gscp(jg)) ) THEN
        CALL finish( TRIM(routine), 'Incorrect setting: ithermo_water = 1 is not developed for chosen inwp_gscp. ')
      END IF

#ifndef __ICON_ART
      IF (inwp_gscp(jg) == 6) THEN
        CALL finish( TRIM(routine),'inwp_gscp == 6, but ICON was compiled with --disable-art')
      ENDIF
#endif
      
      IF (inwp_surface(jg) == 0 .AND. itype_z0 > 1) THEN
        CALL message(TRIM(routine), 'Warning: itype_z0 is reset to 1 because surface scheme is turned off')
        itype_z0 = 1
      ENDIF

      IF (lsgs_cond(jg) .AND. inwp_cldcover(jg) /= 1) THEN
        CALL message(TRIM(routine),'Subgrid-scale condensation only available for cloud cover 1.')
        lsgs_cond(jg) = .FALSE.
      ENDIF

      IF (icpl_aero_gscp > 0 .AND. inwp_gscp(jg) > 3) THEN
        CALL finish( TRIM(routine), 'Aerosol-microphysics coupling currently available only for inwp_gscp=1,2,3')
      ENDIF

#ifdef _OPENACC
      IF ( ALL((/0,1,5/) /= inwp_cldcover(jg)) ) THEN
        CALL finish(routine,'GPU version only available for cloud cover 0, 1 and 5')
      ENDIF

      IF (ALL((/0,1/) /= inwp_sso(jg))) THEN
        CALL finish(routine,'GPU version only available for inwp_sso == 0 or 1.')
      ENDIF

      IF (inwp_gscp(jg) == 8) THEN
        CALL finish(routine,'GPU version not available for Stochastic Bin Microphysics (inwp_gscp=8).')
      ENDIF
#endif

      IF (inwp_surface(jg) == LSS_JSBACH .AND. inwp_turb(jg) /= ivdiff) THEN
        CALL finish( TRIM(routine), 'JSBACH land-surface scheme has to be used with VDIFF turbulence, set inwp_turb=6')
      ENDIF

      IF (inwp_surface(jg) /= LSS_JSBACH .AND. inwp_turb(jg) == ivdiff) THEN
        CALL finish( TRIM(routine), 'VDIFF turbulence scheme has to be used with JSBACH land-surface scheme, set inwp_surface=2')
      ENDIF

      ! For backward compatibility, do not throw an error message, if inwp_turb=10,11 or 12 
      ! is chosen. reset inwp_turb to 1, instead
      IF ( ANY((/10,11,12/) == inwp_turb(jg)) ) THEN
        inwp_turb(jg) = 1
        WRITE(message_text,'(a,i2)') 'Reset inwp_turb to 1 for domain ', jg
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF

      IF (lstoch_expl(jg) .and. lstoch_sde(jg) ) THEN
         lstoch_sde(jg)=.FALSE.
         WRITE(message_text,'(a,i2)') 'lstoch_expl and lstoch_sde cannot both be true. SDE has been switched off for domain ', jg
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF

    ENDDO

    ! deactivate cuda graph if no cpp key => make sure ACC WAIT is activated where needed
#ifndef ICON_USE_CUDA_GRAPH
    lcuda_graph_turb_tran = .FALSE.
#endif


    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg=1,max_dom
      atm_phy_nwp_config(jg)%inwp_convection = inwp_convection(jg)
      atm_phy_nwp_config(jg)%inwp_cldcover   = inwp_cldcover(jg)
      atm_phy_nwp_config(jg)%inwp_radiation  = inwp_radiation(jg)
      atm_phy_nwp_config(jg)%inwp_sso        = inwp_sso(jg)
      atm_phy_nwp_config(jg)%inwp_gwd        = inwp_gwd(jg) 
      atm_phy_nwp_config(jg)%inwp_gscp       = inwp_gscp(jg)
      atm_phy_nwp_config(jg)%inwp_satad      = inwp_satad(jg)
      atm_phy_nwp_config(jg)%inwp_turb       = inwp_turb(jg)
      atm_phy_nwp_config(jg)%inwp_surface    = inwp_surface(jg)
      atm_phy_nwp_config(jg)%itype_z0        = itype_z0

      atm_phy_nwp_config(jg)%nclds              = nclds(jg)
      atm_phy_nwp_config(jg)%lshallowconv_only  = lshallowconv_only(jg)
      atm_phy_nwp_config(jg)%lstoch_expl        = lstoch_expl(jg)
      atm_phy_nwp_config(jg)%lstoch_sde         = lstoch_sde(jg)
      atm_phy_nwp_config(jg)%lstoch_deep        = lstoch_deep(jg)
      atm_phy_nwp_config(jg)%lvvcouple          = lvvcouple(jg)
      atm_phy_nwp_config(jg)%lvv_shallow_deep   = lvv_shallow_deep(jg)
      atm_phy_nwp_config(jg)%lstoch_spinup      = lstoch_spinup(jg)
      atm_phy_nwp_config(jg)%lrestune_off       = lrestune_off(jg)
      atm_phy_nwp_config(jg)%lmflimiter_off     = lmflimiter_off(jg)
      atm_phy_nwp_config(jg)%lgrayzone_deepconv = lgrayzone_deepconv(jg)
      atm_phy_nwp_config(jg)%ldetrain_conv_prec = ldetrain_conv_prec(jg)
      atm_phy_nwp_config(jg)%lsgs_cond          = lsgs_cond(jg)

      atm_phy_nwp_config(jg)%dt_conv         = dt_conv (jg)
      atm_phy_nwp_config(jg)%dt_ccov         = dt_ccov (jg)
      atm_phy_nwp_config(jg)%dt_rad          = dt_rad  (jg)
      atm_phy_nwp_config(jg)%dt_sso          = dt_sso  (jg)
      atm_phy_nwp_config(jg)%dt_gwd          = dt_gwd  (jg)
      atm_phy_nwp_config(jg)%qi0             = qi0 
      atm_phy_nwp_config(jg)%qc0             = qc0 
      atm_phy_nwp_config(jg)%ustart_raylfric = ustart_raylfric 
      atm_phy_nwp_config(jg)%efdt_min_raylfric = efdt_min_raylfric
      atm_phy_nwp_config(jg)%latm_above_top  = latm_above_top(jg)
      atm_phy_nwp_config(jg)%lupatmo_phy     = lupatmo_phy(jg)
      atm_phy_nwp_config(jg)%mu_rain         = mu_rain
      atm_phy_nwp_config(jg)%rain_n0_factor  = rain_n0_factor
      atm_phy_nwp_config(jg)%mu_snow         = mu_snow
      atm_phy_nwp_config(jg)%icpl_aero_gscp  = icpl_aero_gscp
      atm_phy_nwp_config(jg)%lscale_cdnc     = lscale_cdnc
      atm_phy_nwp_config(jg)%icalc_reff      = icalc_reff (jg)
      atm_phy_nwp_config(jg)%icpl_rad_reff   = icpl_rad_reff (jg)
      atm_phy_nwp_config(jg)%ithermo_water   = ithermo_water(jg)
      atm_phy_nwp_config(jg)%lsbm_warm_full  = lsbm_warm_full
    ENDDO

    config_lrtm_filename         = TRIM(lrtm_filename)
    config_cldopt_filename       = TRIM(cldopt_filename)
    config_icpl_aero_conv        = icpl_aero_conv
    config_iprog_aero            = iprog_aero
    config_icpl_o3_tp            = icpl_o3_tp
    config_lcuda_graph_turb_tran = lcuda_graph_turb_tran
    config_icpl_aero_ice         = icpl_aero_ice

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=nwp_phy_nml)                    
      CALL store_and_close_namelist(funit, 'nwp_phy_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=nwp_phy_nml)

  END SUBROUTINE read_nwp_phy_namelist

END MODULE mo_nwp_phy_nml

