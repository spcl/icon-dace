! Namelist for surface physics
!
! these Subroutines are called by control model and construct the
! surface scheme composition
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

MODULE mo_lnd_nwp_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish
  USE mo_impl_constants,      ONLY: SSTICE_ANA, max_nsoil
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,  &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_lnd_nwp_config,      ONLY: config_nlev_snow          => nlev_snow         , &
    &                               config_ntiles             => ntiles_lnd        , &
    &                               config_frlnd_thrhld       => frlnd_thrhld      , &
    &                               config_frlndtile_thrhld   => frlndtile_thrhld  , &
    &                               config_frlake_thrhld      => frlake_thrhld     , &
    &                               config_frsea_thrhld       => frsea_thrhld      , &
    &                               config_hice_min           => hice_min          , &
    &                               config_hice_max           => hice_max          , &
    &                               config_lseaice            => lseaice           , &
    &                               config_lprog_albsi        => lprog_albsi       , &
    &                               config_lbottom_hflux      => lbottom_hflux     , &
    &                               config_llake              => llake             , &
    &                               config_lmelt              => lmelt             , &
    &                               config_lmelt_var          => lmelt_var         , &
    &                               config_lmulti_snow        => lmulti_snow       , &
    &                               config_l2lay_rho_snow     => l2lay_rho_snow    , &
    &                               config_max_toplaydepth    => max_toplaydepth   , &
    &                               config_idiag_snowfrac     => idiag_snowfrac    , &
    &                               config_itype_snowevap     => itype_snowevap    , &
    &                               config_cwimax_ml          => cwimax_ml         , &
    &                               config_c_soil             => c_soil            , &
    &                               config_c_soil_urb         => c_soil_urb        , &
    &                               config_cr_bsmin           => cr_bsmin          , &
    &                               config_itype_trvg         => itype_trvg        , &
    &                               config_itype_evsl         => itype_evsl        , &
    &                               config_itype_lndtbl       => itype_lndtbl      , &
    &                               config_itype_root         => itype_root        , &
    &                               config_itype_canopy       => itype_canopy      , &
    &                               config_cskinc             => cskinc            , &
    &                               config_tau_skin           => tau_skin          , &
    &                               config_lterra_urb         => lterra_urb        , &
    &                               config_lurbalb            => lurbalb           , &
    &                               config_lurbahf            => lurbahf           , &
    &                               config_itype_kbmo         => itype_kbmo        , &
    &                               config_itype_eisa         => itype_eisa        , &
    &                               config_lstomata           => lstomata          , &
    &                               config_l2tls              => l2tls             , &
    &                               config_itype_heatcond     => itype_heatcond    , &
    &                               config_itype_interception => itype_interception, &
    &                               config_itype_hydbound     => itype_hydbound    , &
    &                               config_lana_rho_snow      => lana_rho_snow     , &
    &                               config_lsnowtile          => lsnowtile         , &
    &                               config_sstice_mode        => sstice_mode       , &
    &                               config_sst_td_filename    => sst_td_filename   , &
    &                               config_ci_td_filename     => ci_td_filename    , &
    &                               config_zml_soil           => zml_soil          , &
    &                               config_nlev_soil          => nlev_soil         , &
    &                               config_czbot_w_so         => czbot_w_so

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_nwp_lnd_namelist

CONTAINS


  !-------------------------------------------------------------------------
  !
  !! Read Namelist for NWP land physics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP land physics
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  SUBROUTINE read_nwp_lnd_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: js            ! loop indices
    INTEGER :: iunit

    ! Variable to set the sst and seaice fraction mode
    INTEGER ::  sstice_mode
    !> Action Variables for physical schemes
    ! --------------------------------------
    INTEGER ::  nlev_snow         !< number of snow layers
    INTEGER ::  nlev_soil         !< number of soil layers
    REAL(wp):: zml_soil(max_nsoil)!< Soil full levels
    REAL(wp)::  czbot_w_so        !< thickness of the hydraulical active soil layer [m]
    INTEGER ::  ntiles            !< number of static tiles
    REAL(wp)::  frlnd_thrhld      !< fraction threshold for creating a land grid point
    REAL(wp)::  frlndtile_thrhld  !< fraction threshold for retaining the respective
    !! tile for a grid point
    REAL(wp)::  frlake_thrhld     !< fraction threshold for creating a lake grid point
    REAL(wp)::  frsea_thrhld      !< fraction threshold for creating a sea grid point
    REAL(wp)::  hice_min          !< minimum sea-ice thickness [m]
    REAL(wp)::  hice_max          !< maximum sea-ice thickness [m]
    LOGICAL ::  lbottom_hflux     !< use simple parameterization for heat flux through sea ice bottom
    REAL(wp)::  max_toplaydepth   !< maximum depth of uppermost snow layer for multi-layer snow scheme
    INTEGER ::  itype_trvg        !< type of vegetation transpiration parameterization
    INTEGER ::  itype_evsl        !< type of parameterization of bare soil evaporation (see Schulz and Vogel 2020)
    INTEGER ::  itype_lndtbl      !< choice of table for associating surface parameters to land-cover classes
    INTEGER ::  itype_root        !< type of root density distribution
    INTEGER ::  itype_heatcond    !< type of soil heat conductivity (see Schulz et al. 2016)
    INTEGER ::  itype_interception!< type of plant interception
    REAL(wp)::  cwimax_ml         !< scaling parameter for maximum interception storage
    REAL(wp)::  c_soil            !< surface area density of the (evaporative) soil surface
    REAL(wp)::  c_soil_urb        !< surface area density of the (evaporative) soil surface, urban areas
    REAL(wp)::  cr_bsmin          !< minimum bare soil evaporation resistance (see Schulz and Vogel 2020)
    INTEGER ::  itype_canopy      !< type of canopy parameterisation with respect to the surface energy balance
                                  !< (see Schulz and Vogel 2020)
    REAL(wp)::  cskinc            !< skin conductivity (W/m**2/K)
    REAL(wp)::  tau_skin          !< relaxation time scale for the computation of the skin temperature
    LOGICAL ::  lterra_urb        !< activate urban model TERRA_URB (see Schulz et al. 2023)
    LOGICAL ::  lurbalb           !< use urban albedo and emissivity
    LOGICAL ::  lurbahf           !< use urban anthropogenic heat flux
    INTEGER ::  itype_kbmo        !< type of bluff-body thermal roughness length parameterisation
    INTEGER ::  itype_eisa        !< type of evaporation from impervious surface area
    INTEGER ::  itype_hydbound    !< type of hydraulic lower boundary condition
    INTEGER ::  idiag_snowfrac    !< method for diagnosis of snow-cover fraction
    INTEGER ::  itype_snowevap    !< treatment of snow evaporation in the presence of vegetation

    CHARACTER(LEN=filename_max) :: sst_td_filename, ci_td_filename


    LOGICAL ::           &
         lseaice,        & !> forecast with sea ice model
         lprog_albsi,    & !> sea-ice albedo is computed prognostically 
         llake,          & !> forecast with lake model FLake
         lmelt     ,     & !> soil model with melting process
         lmelt_var ,     & !> freezing temperature dependent on water content
         lmulti_snow,    & !> run the multi-layer snow model
         l2lay_rho_snow, & !> use two-layer snow density for single-layer snow scheme
         lstomata   ,    & !> map of minimum stomata resistance
         l2tls      ,    & !> forecast with 2-TL integration scheme
         lana_rho_snow,  & !> if .TRUE., take rho_snow-values from analysis file
         lsnowtile         !> if .TRUE., snow is considered as a separate tile
    !--------------------------------------------------------------------
    ! nwp forcing (right hand side)
    !--------------------------------------------------------------------

    NAMELIST/lnd_nml/ nlev_snow, zml_soil, ntiles                             , &
         &               frlnd_thrhld, lseaice, lprog_albsi, llake, lmelt     , &
         &               frlndtile_thrhld, frlake_thrhld                      , &
         &               frsea_thrhld, lmelt_var, lmulti_snow                 , &
         &               hice_min, hice_max, lbottom_hflux                    , &
         &               itype_trvg, idiag_snowfrac, max_toplaydepth          , &
         &               itype_evsl                                           , &
         &               itype_lndtbl                                         , &
         &               itype_root                                           , &
         &               itype_heatcond                                       , &
         &               itype_interception                                   , &
         &               itype_hydbound                                       , &
         &               itype_canopy, cskinc, tau_skin                       , &
         &               lterra_urb, lurbalb, lurbahf, itype_kbmo, itype_eisa , &
         &               lstomata                                             , &
         &               l2tls                                                , &
         &               lana_rho_snow, l2lay_rho_snow                        , &
         &               lsnowtile, itype_snowevap                            , &
         &               sstice_mode                                          , &
         &               sst_td_filename                                      , &
         &               ci_td_filename, cwimax_ml, c_soil, c_soil_urb        , &
         &               czbot_w_so, cr_bsmin

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_lnd_nwp_nml:read_nwp_lnd_namelist'

    !-----------------------
    ! 1. default settings
    !-----------------------

    sstice_mode  = SSTICE_ANA  ! forecast mode, sst and sea ice fraction is read from
                               ! the analysis, sst ist kept constant, sea ice fraction
                               ! is modified by the sea ice model
                               ! default names for the time dependent SST and CI ext param files
                               ! if sstice=SSTICE_CLIM, <year> is substituted by "CLIM"
    sst_td_filename = "<path>SST_<year>_<month>_<gridfile>"
    ci_td_filename = "<path>CI_<year>_<month>_<gridfile>"


    nlev_snow      = 2       ! 2 = default value for number of snow layers
    zml_soil(:)    = -1._wp
    czbot_w_so     = 2.5_wp  ! default thickness of 2.5m for hydraulical active soil layer 
    ntiles         = 1       ! 1 = default value for number of static surface types
    frlnd_thrhld   = 0.05_wp ! fraction threshold for creating a land grid point

    frlake_thrhld  = 0.05_wp ! fraction threshold for creating a lake grid point

    frsea_thrhld   = 0.05_wp ! fraction threshold for creating a sea grid point
    frlndtile_thrhld = 0.05_wp ! fraction threshold for retaining the respective 
                             ! tile for a grid point
    hice_min       = 0.05_wp ! minimum sea-ice thickness [m]
    hice_max       = 3.0_wp  ! maximum sea-ice thickness [m]
    lbottom_hflux  = .FALSE. ! true: use simple parameterization for heat flux through sea ice bottom
    lmelt          = .TRUE.  ! soil model with melting process
    lmelt_var      = .TRUE.  ! freezing temperature dependent on water content
    lmulti_snow    = .FALSE. ! .TRUE. = run the multi-layer snow model, .FALSE. = use single-layer scheme
    l2lay_rho_snow = .FALSE. ! use two-layer snow density for single-layer snow model
    max_toplaydepth = 0.25_wp ! maximum depth of uppermost snow layer for multi-layer snow scheme (25 cm)
                              ! (also used for simplified two-layer snow density scheme)
    lsnowtile      = .FALSE. ! if .TRUE., snow is considered as a separate tile
    idiag_snowfrac = 1       ! 1: old method based on SWE, 2: more advanced method used in operational system
    itype_snowevap = 2       ! 1: old method, 2: empirical correction, 3: more advanced empirical correction 
    !
    itype_trvg     = 2       ! type of vegetation transpiration parameterization
                             ! Note that this is currently the only available option!
    itype_evsl     = 2       ! type of bare soil evaporation parameterization
                             !  2: based on BATS (Dickinson 1984)
                             !  4: resistance-based formulation by Schulz and Vogel (2020)
                             !  5: same as 4, but uses cr_bsmin instead of c_soil for tuning, and c_soil is set to 2
    itype_lndtbl   = 3       ! choice of look-up table for associating surface parameters to land-cover classes
    itype_root     = 2       ! type of root density distribution
                             !  1: uniform
                             !  2: exponential
    itype_heatcond = 2       ! type of soil thermal conductivity (see Schulz et al. 2016)
                             !  1: vertically constant, representing a mean soil water content
                             !  2: dependent on soil water content in each layer (Johansen 1975)
    itype_interception = 1   ! type of plant interception
    cwimax_ml      = 1.e-6_wp ! scaling parameter for maximum interception storage. Almost turned off by default;
                              ! the recommended value to activate interception storage is 5.e-4
    c_soil         = 1._wp   ! surface area density of the (evaporative) soil surface
    c_soil_urb     = 1._wp   ! surface area density of the (evaporative) soil surface, urban areas
    cr_bsmin       = 110._wp ! minimum bare soil evaporation resistance (s/m) (see Schulz and Vogel 2020)
    itype_hydbound = 1       ! type of hydraulic lower boundary condition
    !
    itype_canopy   = 1       ! type of canopy parameterisation with respect to the surface energy balance
                             !  1: surface energy balance equation solved at the ground surface,
                             !     canopy energetically not represented
                             !  2: skin temperature formulation by Schulz and Vogel (2020),
                             !     based on Viterbo and Beljaars (1995)
    cskinc         = -1._wp  ! use map of skin conductivity (W/m**2/K)
    tau_skin      = 3600._wp ! relaxation time scale for the computation of the skin temperature
    !
    lterra_urb     = .FALSE. ! if .TRUE., activate urban model TERRA_URB by Wouters et al. (2016, 2017)
                             ! (see Schulz et al. 2023)
    lurbalb        = .TRUE.  ! if .TRUE., use urban albedo and emissivity (Wouters et al. 2016)
    lurbahf        = .TRUE.  ! if .TRUE., use urban anthropogenic heat flux (Wouters et al. 2016)
    itype_kbmo     = 2       ! type of bluff-body thermal roughness length parameterisation
                             !  1: standard SAI-based turbtran (Raschendorfer 2001)
                             !  2: Brutsaert-Kanda parameterisation for bluff-body elements (kB-1)
                             !     (Kanda et al. 2007)
                             !  3: Zilitinkevich (1970)
    itype_eisa     = 3       ! type of evaporation from impervious surface area
                             !  1: evaporation like bare soil (see Schulz and Vogel 2020)
                             !  2: no evaporation
                             !  3: PDF-based puddle evaporation (Wouters et al. 2015)
    !
    lstomata       = .TRUE.  ! map of minimum stomata resistance
    l2tls          = .TRUE.  ! forecast with 2-TL integration scheme
    lana_rho_snow  = .TRUE.  ! if .TRUE., take rho_snow-values from analysis file 

    lseaice        = .TRUE.  ! .TRUE.: sea-ice model is used
    lprog_albsi    = .FALSE. ! .TRUE.: sea-ice albedo is computed prognostically 
                             ! (only takes effect if "lseaice=.TRUE.")
    llake          = .TRUE.  ! .TRUE.: lake model is used

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('lnd_nml')
      READ(funit,NML=lnd_nml)
      CALL close_tmpfile(funit)
    END IF

    !-------------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !-------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('lnd_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, lnd_nml)   ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, lnd_nml)                                       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, lnd_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check (if necessary)
    !----------------------------------------------------

    !Multi-layer snow model
    !
    IF ( (lmulti_snow) .AND. (nlev_snow <= 1) ) THEN
      CALL finish( TRIM(routine),                                   &
        &  'nlev_snow must be >1 when running the multi-layer snow model')
    ENDIF

    IF ( lmulti_snow .AND. l2lay_rho_snow ) THEN
      CALL finish( TRIM(routine), 'multi-layer snow model cannot be combined with l2lay_rho_snow option')
    ENDIF

    IF ( czbot_w_so <= 0._wp) THEN
      CALL finish(TRIM(routine), 'thickness of hydrological active soil layer czbot_w_wo must be > 0.')
    ENDIF

    ! For TERRA_URB: Range of allowed values for itype_kbmo
    IF ((itype_kbmo < 1) .OR. (itype_kbmo > 3)) THEN
      CALL finish(TRIM(routine), 'Value of itype_kbmo must be 1, 2 or 3.')
    ENDIF

    ! For TERRA_URB: Range of allowed values for itype_eisa
    IF ((itype_eisa < 1) .OR. (itype_eisa > 3)) THEN
      CALL finish(TRIM(routine), 'Value of itype_eisa must be 1, 2 or 3.')
    ENDIF

    ! For simplicity, in order to avoid further case discriminations
    IF (l2lay_rho_snow) nlev_snow = 2

    ! Reset prognostic sea-ice albedo switch if the sea-ice scheme is not used
    IF ( .NOT.lseaice ) THEN
      lprog_albsi = .FALSE.  
    ENDIF


    ! Number of actual soil layers

    nlev_soil = count( zml_soil(:) > 0.0_wp )

    ! Check if namelist parameter zml_soil defined in ascending order

    DO js = 1, nlev_soil -1
      IF ( zml_soil(js+1) <= zml_soil(js) ) THEN
        CALL finish(routine, "ERROR namelist parameter zml_soil was not defined in ascending order.")
      ENDIF
    ENDDO

    ! Actual soil layer thickness
    IF (nlev_soil == 0) THEN ! set default values:
      nlev_soil = 8
      ALLOCATE(config_zml_soil(nlev_soil))
      config_zml_soil = (/ 0.005_wp,0.02_wp,0.06_wp,0.18_wp,0.54_wp,1.62_wp,4.86_wp,14.58_wp /)
    ELSE ! use actual values that are defined in namelist setting:
      ALLOCATE(config_zml_soil(nlev_soil))
      config_zml_soil = zml_soil(1:nlev_soil)
    ENDIF
    !$ACC ENTER DATA COPYIN(config_zml_soil)

    !Check if target GPU configuration is supported




    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    config_nlev_snow          = nlev_snow
    config_ntiles             = ntiles
    config_frlnd_thrhld       = frlnd_thrhld
    config_frlndtile_thrhld   = frlndtile_thrhld
    config_frlake_thrhld      = frlake_thrhld
    config_frsea_thrhld       = frsea_thrhld
    config_hice_min           = hice_min
    config_hice_max           = hice_max
    config_lbottom_hflux      = lbottom_hflux
    config_lseaice            = lseaice
    config_lprog_albsi        = lprog_albsi 
    config_llake              = llake
    config_lmelt              = lmelt
    config_lmelt_var          = lmelt_var
    config_lmulti_snow        = lmulti_snow
    config_max_toplaydepth    = max_toplaydepth
    config_idiag_snowfrac     = idiag_snowfrac
    config_itype_snowevap     = itype_snowevap
    config_itype_trvg         = itype_trvg
    config_itype_evsl         = itype_evsl
    config_itype_lndtbl       = itype_lndtbl
    config_itype_root         = itype_root
    config_itype_canopy       = itype_canopy
    config_cskinc             = cskinc
    config_tau_skin           = tau_skin
    config_lterra_urb         = lterra_urb
    config_lurbalb            = lurbalb
    config_lurbahf            = lurbahf
    config_itype_kbmo         = itype_kbmo
    config_itype_eisa         = itype_eisa
    config_lstomata           = lstomata
    config_l2tls              = l2tls
    config_itype_heatcond     = itype_heatcond
    config_itype_interception = itype_interception
    config_cwimax_ml          = cwimax_ml
    config_c_soil             = c_soil
    config_c_soil_urb         = c_soil_urb
    config_cr_bsmin           = cr_bsmin
    config_itype_hydbound     = itype_hydbound
    config_lana_rho_snow      = lana_rho_snow
    config_l2lay_rho_snow     = l2lay_rho_snow
    config_lsnowtile          = lsnowtile
    config_sstice_mode        = sstice_mode
    config_sst_td_filename    = sst_td_filename
    config_ci_td_filename     = ci_td_filename
    config_nlev_soil          = nlev_soil
    config_czbot_w_so         = czbot_w_so
    !$ACC UPDATE DEVICE(config_itype_interception) ASYNC(1)

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=lnd_nml)                    
      CALL store_and_close_namelist(funit, 'lnd_nml') 
    ENDIF


    ! 7. write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=lnd_nml)

  END SUBROUTINE read_nwp_lnd_namelist


END MODULE mo_lnd_nwp_nml

