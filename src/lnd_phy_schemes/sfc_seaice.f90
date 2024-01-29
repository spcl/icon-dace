! The main program unit of the sea-ice parameterization scheme for NWP.
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
!
!
! The main program unit of the sea-ice parameterization scheme for NWP.
!
! It contains a number of procedures.
! SUBROUTINE seaice_init_nwp
! initializes the sea-ice scheme and performs some consistency checks.
! SUBROUTINE seaice_timestep_nwp
! advances prognostic variables of the sea-ice scheme one time step.
! SUBROUTINE seaice_coldinit_nwp
! performs a cold start of the sea-ice scheme.
! SUBROUTINE seaice_coldinit_albsi_nwp
! performs a cold start initialization of prognostic sea-ice albedo.
!
! The present sea-ice parameterization scheme is a bulk thermodynamic (no rheology) scheme
! intended for use in NWP and similar applications.
! The scheme is based on a self-similar parametric representation (assumed shape) of the
! evolving temperature profile within the ice and on the integral heat budget of the ice slab.
! The scheme carries ordinary differential equations (in time)
! for the ice surface temperature and the ice thickness.
! An explicit Euler scheme is used for time advance.
! In the current configuration of the scheme, snow over sea ice is not treated explicitly.
! The effect of snow above the ice is accounted for implicitly (parametrically).
! To this end, a rate equation for the sea-ice surface albedo
! with respect to (diffuse) solar radiation is used.
! The rate equation contains
! the relaxation terms that drive sea-ice albedo towards its equilibrium value,
! and "albedo source term" due to precipitation that accounts for the increase
! of albedo after snowfalls.
! The equilibrium albedo is a function of the sea-ice surface temperature.
! Optionally, the sea-ice albedo may be treated diagnostically using
! a temperature-dependent equilibrium albedo.
! For the "sea water" type ICON grid boxes, the snow thickness is set to zero and
! the snow surface temperature is set equal to the ice surface temperature
! (both temperatures are set equal to the fresh-water freezing point if the ice is absent).
!
! Prognostic equations for the ice thickness and the ice surface temperature are solved
! for the ICON grid boxes with the ice fraction
! (area fraction of a given model grid box of the type "sea water" that is covered by ice)
! that exceeds a threshold value of 0.0015. Otherwise, the grid box is treated as ice-free.
! The ice fraction is determined on the basis of observational data
! by the data assimilation scheme and is kept constant over the entire model forecast period.
! However, if the ice melts away during the forecast, the ice fraction is reset to zero.
! This is done within SUBROUTINE update_idx_lists_sea.
! If the ICON grid box is set ice-free during the initialization,
! no ice is created over the forecast period.
! If observational data indicate open water conditions for a given ICON grid box,
! residual ice from the previous model run is removed,
! i.e., the ice thickness is set to zero and
! the ice surface temperature is set to the fresh-water freezing point.
! The newly formed ice has the surface temperature equal to the salt-water freezing point
! and the thickness from 0.1 m to 0.5 m depending on the ice fraction.
! The new ice is formed instantaneously if the data assimilation scheme
! indicates the presence of ice in a given ICON grid box
! but there was no ice in that grid box at the end of the previous model run.
! Prognostic ice thickness is limited by a maximum value of 3 m and a minimum value of 0.05 m.
! Constant values of the ice density, ice molecular heat conductivity, specific heat of ice,
! the latent heat of fusion, and the salt-water freezing point are used.
!
! In uncoupled runs (no interaction of the ice slab with the ocean boundary layer beneath the ice),
! a simple ad hoc formulation of the heat flux from water to ice is used.
! Optionally, adaptive tuning of the parameter(s) of the temperature profile within the ice
! and of the heat flux from water to ice can be used that is based on the assimilation increments
! of the atmopsheric surface layer quantities.
!
! A detailed description of the sea-ice scheme is given in
! Mironov, D., B. Ritter, J.-P. Schulz, M. Buchhold, M. Lange, and E. Machulskaya, 2012:
! Parameterization of sea and lake ice in numerical weather prediction models
! of the German Weather Service.
! Tellus A, 64, 17330. doi:10.3402/tellusa.v64i0.17330
!
! The present sea-ice scheme (with minor modifications)
! is also implemented into the NWP models GME and COSMO
! (see Mironov et al. 2012, for details).
!
! Adaptive parameter tuning in NWP models is described in
! Zaengl, G., 2023:
! Adaptive tuning of uncertain parameters in a numerical weather prediction model
! based upon data assimilation.
! Q. J. R. Meteorol. Soc., 1-20. doi:https://doi.org/10.1002/qj.4535
!
!
!
! Lines embraced with "!_tmp>" and "!_tmp<" contain temporary parts of the code.
! Lines embraced/marked with "!_dev>" and "!_dev<" may be replaced
! as improved formulations are developed and tested.
! Lines embraced/marked with "!_cdm>" and "!_cdm<" are DM's comments that may be helpful to a user.
! Lines embraced/marked with "!_dbg>" and "!_dbg<" are used for debugging purposes only.
! Lines starting with "!_nu" are not used.

MODULE sfc_seaice

  USE mo_kind, ONLY:      &
                   &  wp  !< KIND-type parameter for real variables

  USE mo_exception, ONLY:           &          
                        &  finish , &  !< external procedure, finishes model run and reports the reason  
                        &  message     !< external procedure, sends a message (error, warning, etc.)

  USE mo_impl_constants, ONLY :     &
                              &  ALB_SI_MISSVAL             !< missing value for prognostic seaice albedo

!_cdm>
! Note that ki is equal to 2.1656 in ICON, but is 2.29 in COSMO and GME. 
!_cdm<
  USE mo_physical_constants, ONLY:                      &
                                 & tf_fresh => tmelt  , &  !< fresh-water freezing point [K]
                                 &             tf_salt, &  !< salt-water freezing point [K]
                                 &             alf    , &  !< latent heat of fusion [J/kg]
                                 &             rhoi   , &  !< density of ice [kg/m^3]
                                 &             ci     , &  !< specific heat of ice [J/(kg K)]
                                 &             ki          !< molecular heat conductivity of ice [J/(m s K)]  

  USE mo_lnd_nwp_config,     ONLY:                      &
                                 & lprog_albsi        , &  !< sea-ice albedo is computed prognostically 
                                 & lbottom_hflux      , &  !< use parameterization for heat flux through sea ice bottom
                                 & frsi_min           , &  !< minimum sea-ice fraction  [-]
                                 & hice_min           , &  !< minimum sea-ice thickness [m]
                                 & hice_max                !< maximum sea-ice thickness [m]

  USE sfc_terra_data,        ONLY:                      &
                                 & csalb              , &  !< solar albedo for different soil types
                                 & ist_seaice              !< ID of soiltype "sea ice"

  USE mo_coupling_config,    ONLY: is_coupled_to_ocean     !< TRUE for coupled ocean-atmosphere runs 

  IMPLICIT NONE

  PRIVATE

  REAL (wp), PARAMETER ::                             &
                       &  hice_ini_min = 0.1_wp     , &  !< minimum thickness of the newly formed sea ice [m]
                       &  hice_ini_max = 0.5_wp     , &  !< maximum thickness of the newly formed sea ice [m]
                       &  csi_lin      = 0.5_wp     , &  !< shape factor for linear temperature profile [-]
                       &  phiipr0_lin  = 1.0_wp     , &  !< derivative (at zeta=0) of the linear 
                                                         !< temperature profile shape function [-]
                       &  csidp_nlin   = 2.0_wp     , &  !< disposable parameter for non-linear 
                                                         !< temperature profile shape factor [-]
                       &  csidp_nlin_d =              &  !< derived disposable parameter for non-linear 
                                                         !< temperature profile shape factor [-]
                       &  (1._wp+csidp_nlin)/12._wp , &
                       &  cmaxearg     = 1.0E+02_wp , &  !< maximum value of the EXP function argument (security constant) [-]
                       &  csmall       = 1.0E-05_wp      !< small number (security constant) [-]
 
  REAL (wp), PARAMETER ::                             &
    &  taualbsi_min   = 7.0_wp*86400._wp            , &  !< minimum relaxation time scale for sea-ice albedo [s]
    &  taualbsi_max   = 21.0_wp*86400._wp           , &  !< maximum relaxation time scale for sea-ice albedo [s]
    &  t_taualbsi_min = 268.15_wp                   , &  !< lower bound of the temperature range in the interpolation
                                                         !< formula for the sea-ice albedo relaxation time scale [K]
    &  rdelt_taualbsi =                               &  !< reciprocal of the temperature range in the interpolation
    &    1._wp/(tf_fresh-t_taualbsi_min)            , &  !< formula for the sea-ice albedo relaxation time scale [K^{-1}]
    &  albsi_snow_max = 0.80_wp                     , &  !< maximum albedo of snow over sea ice [-]
    &  albsi_snow_min = 0.50_wp                     , &  !< minimum albedo of snow over sea ice [-]
    &  c1_albsi_snow  =                               &  !< constant in the expression for albedo 
    &    1._wp-albsi_snow_min/albsi_snow_max        , &  !< of snow over sea ice [-] 
    &  c2_albsi_snow  = 136.6_wp/tf_fresh           , &  !< constant in the expression for albedo
                                                         !< of snow over sea ice [K^{-1}] 
    &  c_tausi_snow   = 1._wp/5._wp                 , &  !< constant used to define the relaxation time scale towards
                                                         !< the equilibrium albedo of snow over sea ice [(kg/m^2)^{-1}]
                                                         !< (corresponds to 5 mm snow water equaivalent precipitated over 
                                                         !< an e-folding time scale)
    &  t_albsi_snow_max  = 272.95_wp                     !< upper bound of the temperature range over which relaxation 
                                                         !< towards snow-overice albedo is applied [K]

  !>
  !! Optical characteristics of sea ice. 
  !!
  !! A storage for an n-band approximation of the exponential decay law 
  !! for the flux of solar radiation is allocated.
  !! A maximum value of the wave-length bands is currently set equal to two.
  !!
  
  INTEGER, PARAMETER ::                      &
                     &  nband_optic_max = 2  !< maximum number of wave-length bands in the decay law
                                             !< for the solar radiation flux [-]

  TYPE opticpar_seaice
    INTEGER ::                                       &
            &  nband_optic                              !< number of wave-length bands [-]
    REAL (wp) ::                                     &
              &  frac_optic(nband_optic_max)       , &  !< fractions of total solar radiation flux for different bands [-]
              &  extincoef_optic(nband_optic_max)       !< extinction coefficients for different bands [1/m]
  END TYPE opticpar_seaice

  ! One-band approximation for opaque sea ice.
  ! The use of large extinction coefficient prevents the penetration of solar radiation
  ! into the ice interior, i.e. the volumetric character of the solar radiation heating is ignored.
  TYPE (opticpar_seaice), PARAMETER ::                             &
                                    &  opticpar_seaice_opaque =    & 
                                    &  opticpar_seaice(1,          &
                                    &  (/1._wp, 0._wp/),           &
                                    &  (/1.0E+07_wp, 1.E+10_wp/)) 

  ! Minimum values of the sea-ice fraction and of the sea-ice thickness are used outside 
  ! "sfc_seaice" to compose an index list of grid boxes where sea ice is present. 
  PUBLIC ::                       &
         &  hice_ini_min             , & ! parameter
         &  hice_ini_max             , & ! parameter
         &  seaice_init_nwp          , & ! procedure     
         &  seaice_coldinit_nwp      , & ! procedure
         &  seaice_coldinit_albsi_nwp, & ! procedure
         &  seaice_timestep_nwp      , & ! procedure
         &  alb_seaice_equil

#ifdef ICON_USE_CUDA_GRAPH
  LOGICAL, PARAMETER :: using_cuda_graph = .TRUE.
#else
  LOGICAL, PARAMETER :: using_cuda_graph = .FALSE.
#endif

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

CONTAINS

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Prognostic variables of the sea-ice parameterization scheme are initialized 
  !! and some consistency checks are performed. 
  !!
  !! The procedure arguments are arrays (vectors) of the sea-ice fraction, and 
  !! of the sea-ice scheme prognostic variables.
  !! The vector length is equal to the number of grid boxes (within a given block)
  !! where sea ice is present, i.e. where the sea-ice fraction exceeds its minimum value.
  !! First, the sea-ice fraction is checked. 
  !! If a value less than a minimum threshold value is found 
  !! (indicating that a grid box is declared as partially ice-covered but no sea ice 
  !! should be present), an error message is sent and the model abort is called.
  !! Next, "new" ice is formed in the grid boxes where "old" value of the sea-ice thickness
  !! is less than a minimum threshold value (i.e. there was no ice at the end 
  !! of the the previous model run). The newly formed ice has the surface temperature equal 
  !! to the salt-water freezing point and the thickness of 0.1 m to 0.5 m  
  !! depending on the ice fraction.
  !! For the grid boxes where new ice is created, prognostic sea-ice albedo is set equal 
  !! to its equilibrium value (function of sea-ice surface temperature).
  !! Then, the ice thickness is limited from above and below 
  !! by the maximum and minimum threshold values, 
  !! and the ice surface temperature is limited from above by the fresh-water freezing point. 
  !! These security measures (taken for each grid box where sea ice is present)
  !! are required to avoid non-allowable values of prognostic variables 
  !! that may occure due to a loss of accuracy during the model IO 
  !! (e.g. due to GRIB encoding and decoding).  
  !! Finally, the snow thickness is set to zero, and
  !! the snow surface temperature is set equal to the ice surface temperature 
  !! (recall that snow over sea ice is not treated explicitly).
  !! 

  SUBROUTINE seaice_init_nwp (                                  &
                          &  linit_hice,                        &
                          &  nsigb,                             &
                          &  frsi,                              &
                          &  tice_p, hice_p, tsnow_p, hsnow_p,  &
                          &  albsi_p,                           &
                          &  tice_n, hice_n, tsnow_n, hsnow_n,  &
                          &  albsi_n                            &
                          &  )

    IMPLICIT NONE

    ! Procedure arguments 

    LOGICAL, INTENT(IN) ::        &
                        &  linit_hice     !< TRUE if initialization of hice is required 
                                          !< for new seaice points.
                                          !< FALSE if updated information on hice is provided 
                                          !< (e.g. in coupled atmosphere-ocean runs) 

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< Array (vector) dimension
                                  !< (equal to the number of grid boxes within a block 
                                  !< where the sea ice is present) 

    REAL(wp), DIMENSION(:), INTENT(IN)    ::         &
                                          &  frsi       !< sea-ice fraction [-]
 
    REAL(wp), DIMENSION(:), INTENT(INOUT) ::           &
                                          &  tice_p  , &  !< temperature of ice upper surface at previous time level [K]
                                          &  hice_p  , &  !< ice thickness at previous time level [m] 
                                          &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K]
                                          &  hsnow_p , &  !< snow thickness at previous time level [m] 
                                          &  albsi_p , &  !< sea-ice albedo at previous time level [-] 
                                          &  tice_n  , &  !< temperature of ice upper surface at new time level [K] 
                                          &  hice_n  , &  !< ice thickness at new time level [m] 
                                          &  tsnow_n , &  !< temperature of snow upper surface at new time level [K] 
                                          &  hsnow_n , &  !< snow thickness at new time level [m] 
                                          &  albsi_n      !< sea-ice albedo at new time level [-] 

    ! Local variables 

    INTEGER ::      &
            &  isi  !< DO loop index

    CHARACTER(len=256) ::           &
                       &  nameerr , &  !< name of procedure where an error occurs
                       &  texterr      !< error/warning message text

    LOGICAL ::             &
            &  lcallabort  !< logical switch, set .TRUE. if errors are encountered 
                           !< (used to call modell abort outside a DO loop)


    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Logical switch, default value is .FALSE. 
    lcallabort = .FALSE.

    ! Loop over grid boxes where sea ice is present
    GridBoxesWithSeaIce: DO isi=1, nsigb

      ! Check sea-ice fraction
      IF( frsi(isi) < frsi_min ) THEN 
        ! Sea-ice fraction less than a minimum threshold value is found
        ! Set logical switch 
        lcallabort = .TRUE.  
        ! Exit DO loop to call model abort
        EXIT GridBoxesWithSeaIce 
      END IF 

      IF ( linit_hice ) THEN
        ! Here we assume that new seaice points are characterized by (fr_seaice>0, hice_p<hice_min)
        ! Create new ice
        IF( hice_p(isi) < (hice_min-csmall) ) THEN
          hice_p(isi) = hice_ini_min + frsi(isi) * (hice_ini_max-hice_ini_min)
          tice_p(isi) = tf_salt
          ! Set sea-ice albedo to its equilibrium value
          ! (only required if sea-ice albedo is treated prognostically)
          IF ( lprog_albsi ) THEN
            albsi_p(isi) = alb_seaice_equil( tice_p(isi) )
          ENDIF
        ENDIF

      ELSE ! linit_hice=.FALSE.
        !
        ! Here we assume that an educated guess for ice thickness is provided, 
        ! with hice>=hice_min (e.g. in coupled mode) in addition to the ice fraction. 
        ! Hence, we skip the initialization of ice thickness but keep 
        ! the initialization of ice temperature and albedo for new ice points.
        !
        ! In this case the identification of new ice points by the above condition 
        ! (fr_seaice>0, hice_p<hice_min) fails and so does any identification via tice. 
        ! Hence, new ice points must be identified before entering 
        ! this routine by comparing the old and new seaice fraction fields and 
        ! by creating a "NEW seaice points only" index list.
        !
        ! For this IF condition to work correctly, seaice_init_nwp must be called for 
        ! NEW seaice points only rather than for ALL seaice points 
        ! (see example in mo_nwp_sfc_utils:process_sst_and_seaice).
        ! Otherwise ice temperature and albedo will erroneously be reset 
        ! to default values for all seaice points.
        !
        tice_p(isi) = tf_salt
        ! Set sea-ice albedo to its equilibrium value
        ! (only required if sea-ice albedo is treated prognostically)
        IF ( lprog_albsi ) THEN
          albsi_p(isi) = alb_seaice_equil( tice_p(isi) )
        ENDIF
      ENDIF


      ! In general we assume that new seaice points are characterized by 
      ! ( fr_seaice>0, h_ice_p=0 ). However, it may happen that h_ice_p 
      ! has already been adjusted consistently by the data assimilation process. 
      ! In that case, we fail to identify new sea-ice points by the above condition, 
      ! and we miss the initialization of the prognostic seaice albedo. Thus, 
      ! the following statement is added. 
      IF ( lprog_albsi .AND. albsi_p(isi) <= 0._wp) THEN
        albsi_p(isi) = alb_seaice_equil( tice_p(isi) )
      ENDIF

      ! Take security measures 
      hice_p(isi) = MAX(MIN(hice_p(isi), hice_max), hice_min) 
      tice_p(isi) = MIN(tice_p(isi), tf_fresh)

      ! Set temperature of snow upper surface and snow thickness
      tsnow_p(isi) = tice_p(isi)
      hsnow_p(isi) = 0._wp

      ! Set variables at new time level
      tice_n(isi)  = tice_p(isi)     
      hice_n(isi)  = hice_p(isi)     
      tsnow_n(isi) = tsnow_p(isi)    
      hsnow_n(isi) = hsnow_p(isi)    
      IF ( lprog_albsi ) THEN
        albsi_n(isi) = albsi_p(isi)
      ENDIF

    END DO GridBoxesWithSeaIce  

    ! Call model abort if errors are encountered
    IF( lcallabort ) THEN 
      ! Send an error message 
      WRITE(nameerr,*) "MODULE sfc_seaice, SUBROUTINE seaice_init_nwp" 
      WRITE(texterr,*) "Sea-ice fraction ", frsi(isi),                        & 
                    &  " is less than a minimum threshold value ", frsi_min,  & 
                    &  " Call model abort."
      CALL message(TRIM(nameerr), TRIM(texterr))

      ! Call model abort
      WRITE(nameerr,*) "sfc_seaice:seaice_init_nwp" 
      WRITE(texterr,*) "error in sea-ice fraction"
      CALL finish(TRIM(nameerr), TRIM(texterr))
    END IF 

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_init_nwp


!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Prognostic variables of the sea-ice scheme are advanced one time step.
  !!
  !! Ordinary differential equations (in time) for the ice surface temperature and 
  !! the ice thickness are solved using an explicit Euler scheme for time advance.
  !! The sea-ice surface albedo with respect to solar radiation is determined  
  !! by solving a rate equation (if the sea-ice albedo is treated diagnostically, 
  !! no albedo calculations are performed in the present routine).
  !! The shape factor for the temperature profile within the ice and the derivative of the 
  !! temperature profile shape function at the underside of the ice are functions of the ice 
  !! thickness (see Mironov et al. 2012, for details).  
  !! Optionally, constant values of the shape factor and of the shape-function derivative
  !! corresponding to the linear temperature profile within the ice 
  !! (cf. the sea-ice schemes of GME and COSMO) can be used 
  !! (as there is no logical switch to activate this option, changes in the code should be made).
  !! In the regime of ice growth or melting from below, the solution may become spurious
  !! when the ice thickness is small and/or the model time step is large.
  !! In such a case, a quasi-steady heat transfer through the ice is assumed 
  !! and the differential equation for the ice surface temperature 
  !! is reduced to an algebraic relation. 
  !! In the current configuration of the sea-ice scheme, snow over ice is not treated explicitly. 
  !! The effect of snow is accounted for implicitely through changes
  !! of the sea-ice albedo with respect to solar radiation.
  !! For the "sea water" type grid boxes, the snow thickness is set to zero and
  !! the snow surface temperature is set equal to the ice surface temperature.
  !! Prognostic ice thickness is limited by a maximum value of 3 m and a minimum value of 0.05 m. 
  !! In "uncoupled" runs, i.e., where the interaction between the sea ice and the sea water 
  !! beneath is not considered, an ad hoc formulation for the heat flux from water to ice is used
  !! that serves to quench the ice growth rate as the ice thickness approaches its maximum value.
  !! No ice is created during the forecast period. 
  !! If the ice melts away during the forecast (i.e., the ice becomes thinner than 0.05 m), 
  !! the ice tickness is set to zero and 
  !! the ice surface temperatures is set to the fresh-water freezing point.
  !! The procedure arguments are arrays (vectors) 
  !! of the sea-ice scheme prognostic variables,
  !! of the components of the heat balance at the ice upper surface  
  !! (i.e., the fluxes of sensible and latent heat, the net flux of long-wave radiation,
  !! and the net flux of solar radiation with due regard for the ice surface albedo),
  !  of the precipitation rates of snow and rain,
  !! and of the sea-ice surface albedo with respect to solar radiation.
  !! Fluxes are positive when directed downward.
  !! The vector length is equal to the number of model grid boxes 
  !! (within a given block) where sea ice is present  
  !! (i.e., where the sea-ice fraction exceeds its minimum threshold value).
  !! The time tendecies of the ice thickness, the ice surface temperature, 
  !! the snow thickness and the snow surface temperature are also computed. 
  !! These are optional arguments of the procedure. 
  !!

  SUBROUTINE seaice_timestep_nwp (                                      &
                              &  dtime,                                 &
                              &  nsigb,                                 &
                              &  qsen, qlat, qlwrnet, qsolnet,          &
                              &  snow_rate, rain_rate, fac_bottom_hflx, &
                              &  tice_p, hice_p, tsnow_p, hsnow_p,      &
                              &  albsi_p,                               &
                              &  tice_n, hice_n, tsnow_n, hsnow_n,      &
                              &  condhf, albsi_n,                       &
                              &  opt_dticedt, opt_dhicedt, opt_dtsnowdt,&
                              &  opt_dhsnowdt                           )

    IMPLICIT NONE

    ! Procedure arguments 

    REAL(wp), INTENT(IN) ::        & 
                         &  dtime  !< model time step [s]

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< number of grid boxes within a block 
                                  !< where the sea ice is present (<=nproma)

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                       &  qsen    , &  !< sensible heat flux at the surface [W/m^2]
                                       &  qlat    , &  !< latent heat flux at the surface [W/m^2]
                                       &  qlwrnet , &  !< net long-wave radiation flux at the surface [W/m^2] 
                                       &  qsolnet      !< net solar radiation flux at the surface [W/m^2] 

    REAL(wp), DIMENSION(:), INTENT(IN) ::             &
                                       &  snow_rate , &  !< snow rate (convecive + grid-scale) [kg/(m^2 s)]
                                       &  rain_rate      !< rain rate (convecive + grid-scale) [kg/(m^2 s)]

    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(IN) ::   &
                                       &  fac_bottom_hflx  !< tuning factor for heat flux from water to ice [-]

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                       &  tice_p  , &  !< temperature of ice upper surface at previous time level [K] 
                                       &  hice_p  , &  !< ice thickness at previous time level [m] 
                                       &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K] 
                                       &  hsnow_p , &  !< snow thickness at previous time level [m] 
                                       &  albsi_p      !< sea-ice albedo at previous time level [-] 

    REAL(wp), DIMENSION(:), INTENT(OUT) ::          &
                                       &  tice_n  , &  !< temperature of ice upper surface at new time level [K] 
                                       &  hice_n  , &  !< ice thickness at new time level [m] 
                                       &  tsnow_n , &  !< temperature of snow upper surface at new time level [K] 
                                       &  hsnow_n , &  !< snow thickness at new time level [m] 
                                       &  condhf  , &  !< conductive heat flux within the sea ice
                                                       !< just above the ice lower boundary [W/m^2]
                                       &  albsi_n      !< sea-ice albedo at new time level [-] 

    REAL(wp), DIMENSION(:), INTENT(OUT), OPTIONAL ::    &
                                       &  opt_dticedt , &  !< time tendency of ice surface temperature [K/s] 
                                       &  opt_dhicedt , &  !< time tendency of ice thickness [m/s] 
                                       &  opt_dtsnowdt, &  !< time tendency of snow surface temperature [K/s] 
                                       &  opt_dhsnowdt     !< time tendency of snow thickness [m/s] 

    ! Derived parameters 
    ! (combinations of physical constants encountered several times in the code)

    REAL (wp), PARAMETER ::                                   &  
                         &  r_rhoici     = 1._wp/(rhoi*ci)  , &  
                         &  ki_o_rhoici  = ki*r_rhoici      , & 
                         &  r_rhoialf    = 1._wp/(rhoi*alf) , &
                         &  ki_o_rhoialf = ki*r_rhoialf     , & 
                         &  ci_o_alf     = ci/alf        

    ! Local variables 
    REAL(wp), DIMENSION(nsigb) ::            &
                                &  dticedt , &  !< time tendency of ice surface temperature [K/s] 
                                &  dhicedt , &  !< time tendency of ice thickness [m/s] 
                                &  dtsnowdt, &  !< time tendency of snow surface temperature [K/s] 
                                &  dhsnowdt     !< time tendency of snow thickness [m/s] 

    INTEGER ::      &
            &  isi  !< DO loop index

    REAL (wp) ::                &
              &  qsoliw       , &  !< solar radiation flux at the ice-water interface (positive downward) [W/m^2]
              &  qatm         , &  !< "total atmospheric heat flux" for the ice slab  (positive downward) [W/m^2]
                                   !< (sum of the sensible heat flux, latent heat flux and the net flux 
                                   !< of long-waver radiation at the ice upper surface, 
                                   !< and the difference of solar radiation fluxes 
                                   !< at the upper and lower surfaces of the ice slab)
              &  qwat         , &  !< heat flux from water to ice (positive downward) [W/m^2]
              &  csice        , &  !< shape factor for the temperature profile within the ice [-]
              &  phiipr0      , &  !< derivative (at zeta=0) of the temperature profile shape function [-]
              &  hice_thrshld , &  !< threshold value of the ice tickness to switch between a quasi-equilibrium
                                   !< model of heat transfer through the ice and a complete model [m]
              &  rti          , &  !< dimensionless parameter [-]
              &  r_dtime      , &  !< reciprocal of the time step [1/s]
              &  strg_1       , &  !< temporary storage variable 
              &  strg_2            !< temporary storage variable 

    REAL (wp) ::                &
              &  albsi_e      , &  !< equilibrium sea-ice albedo [-]
              &  albsi_snow_e , &  !< equilibrium albedo of snow over sea ice [-]
              &  taualbsi     , &  !< relaxation time scale for sea-ice albedo [s]
              &  rtaualbsisn  , &  !< reciprocal of relaxation time scale for snow-over-ice albedo [s^{-1}]
              &  albsi_e_wghtd     !< weighted equilibrium albedo (storage variable) [-] 

    LOGICAL ::   lis_coupled_run   !< TRUE for coupled ocean-atmosphere runs (copy for vectorisation)
    LOGICAL ::   lpres_fac_hflux   !< TRUE if fac_bottom_hflx is provided as input (from NWP interface only)

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Reciprocal of the time step 
    r_dtime = 1._wp/dtime

    ! If fac_bottom_hflx is provided, adaptive tuning of the parameter(s) 
    ! of the temperature profile within the ice and of the heat flux from water to ice is used
    lpres_fac_hflux = PRESENT(fac_bottom_hflx)

    ! for vectorisation
    lis_coupled_run = is_coupled_to_ocean()
    !$ACC DATA CREATE(dticedt, dhicedt, dtsnowdt, dhsnowdt) &
    !$ACC   PRESENT(qsen, qlat, qlwrnet, qsolnet, snow_rate, rain_rate, tice_p, hice_p, tsnow_p, hsnow_p) &
    !$ACC   PRESENT(albsi_p, tice_n, hice_n, tsnow_n, hsnow_n, condhf, albsi_n)
    !$ACC DATA PRESENT(fac_bottom_hflx) IF(PRESENT(fac_bottom_hflx))

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR &
    !$ACC   PRIVATE(qsoliw, qatm, qwat, csice, phiipr0, hice_thrshld, rti, strg_1, strg_2) &
    !$ACC   PRIVATE(albsi_e, albsi_snow_e, taualbsi, rtaualbsisn, albsi_e_wghtd)
    ! Loop over grid boxes where sea ice is present
    GridBoxesWithSeaIce: DO isi=1, nsigb

      ! Compute solar radiation flux at the ice-water interface (positive downward) 
      qsoliw = qsolnet(isi)*(                                                                 &
             & opticpar_seaice_opaque%frac_optic(1)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(1)*hice_p(isi), cmaxearg)) +  &
             & opticpar_seaice_opaque%frac_optic(2)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(2)*hice_p(isi), cmaxearg)) )

      ! Compute total atmospheric heat flux for the ice slab  (positive downward) 
      qatm = qsen(isi) + qlat(isi) + qlwrnet(isi) + qsolnet(isi) - qsoliw 

      ! Provision is made to account for the heat flux from water to ice (upward flux is negative).
      ! This is the ocean heat flux just below the ice, 
      ! not including the latent flux from melting and freezing.
      ! In the case of coupled icon atmosphere-ocean runs, 
      ! this flux has to be passed from the ocean if available.
      ! For uncoupled runs, an ad hoc formulation for the heat flux from water to ice is used 
      ! that basically serves to quench the ice growth rate 
      ! as the ice thickness approaches its maximum value.

      IF (lbottom_hflux) THEN
        ! Derivative (at zeta=0) of the temperature profile shape function.
        IF (lpres_fac_hflux) THEN 
          ! Adaptive tuning of phiipr0 and of bottom heat flux
          phiipr0 = phiipr0_lin - fac_bottom_hflx(isi)*hice_p(isi)/hice_max
        ELSE
          ! A constant value is used (cf. linear temperature profile).
          phiipr0 = phiipr0_lin
        END IF
        ! Heat flux from water to ice is limited from above (no flux from ice to water is allowed).
        qwat = (1._wp-hice_p(isi)/hice_max-phiipr0)*ki*(tf_salt-tice_p(isi))/MAX(hice_p(isi),hice_min)
        qwat = MIN(qwat, 0._wp)
      ELSE
        ! Derivative (at zeta=0) of the temperature profile shape function 
        phiipr0 = 1._wp - hice_p(isi)/hice_max
        ! Heat flux from water to ice is set to zero.
        qwat = 0._wp 
      END IF

      ! Compute temperature profile shape factor and temporary storage variables
      ! (to recover linear temperature profile, set csice=csi_lin)
      csice = csi_lin - csidp_nlin_d*hice_p(isi)/hice_max
      rti = ci_o_alf*(tice_p(isi)-tf_salt)
      strg_1 = rti*(1.5_wp-2.0_wp*csice)
      strg_2 = 1._wp + strg_1 

      FreezingMeltingRegime: IF( tice_p(isi)>=(tf_fresh-csmall) .AND. qatm>0._wp ) THEN 

        ! Melting from above

        ! Set the ice surface temperature equal to the fresh-water freezing point
        tice_n(isi) = tf_fresh

        IF ( lis_coupled_run ) THEN
          ! Coupling to icon-o: 
          ! Heat flux within the ice just above the ice-water intrface, phiipr0 is computed above 
          ! Note that for coupled runs, lbottom_hflux=.FALSE. and the heat flux 
          ! from water to ice qwat is set to zero (qwat should be provided by the ocean model).
          condhf(isi) = ki*phiipr0*(tice_n(isi)-tf_salt)/hice_p(isi)
          ! No change of ice thickness
          hice_n(isi) = hice_p(isi)
        ELSE
          ! Compute the rate of ice melting (note the sign of heat fluxes)
          dhicedt(isi) = -(qatm-qwat)*r_rhoialf/strg_2 
          ! Update the ice thickness
          hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)
        END IF

      ELSE FreezingMeltingRegime 

        ! Freezing or melting from below

        ! Compute threshold value of the ice thickness  
        ! Note that an expression in parentheses should be multiplied with 
        ! MAX(1._wp, ABS(2._wp*csice*rti)) (cf. the code of the lake parameterization scheme FLake).
        ! However, |2*csice*(ci/alf)*(tice-tf_salt)| < 1 
        ! at all conceivable values of the ice surface temperature.
        hice_thrshld = SQRT(phiipr0*ki_o_rhoici*dtime/csice)

        IF( hice_p(isi)<hice_thrshld) THEN 

          ! Use a quasi-equilibrium model of heat transfer through the ice 

          IF ( lis_coupled_run ) THEN
            ! Coupling to icon-o: 
            ! Heat flux (phiipr0 is computed above) 
            condhf(isi) = ki*phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)
            ! No change of ice thickness
            hice_n(isi) = hice_p(isi)
          ELSE

            ! Compute the time-rate-of-change of the ice thickness (note the sign of heat fluxes)
            dhicedt(isi) = -(qatm-qwat)*r_rhoialf/strg_2
            ! Update the ice thickness
            hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)

          END IF

          ! Compute the (updated) ice surface temperature (note the sign of heat fluxes)
          tice_n(isi) = tf_salt + (qatm+qwat*strg_1)*hice_n(isi)/(phiipr0*ki*strg_2)

        ELSE

          ! Use a complete model of heat transfer through the ice 

          ! Use dhsnowdt as a temporary storage
          dhsnowdt(isi) = phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)

          ! Compute the time-rate-of-change of the ice surface temperature 
          ! (note the sign of heat fluxes)
          dticedt(isi) = ( (qatm+strg_1*qwat)*r_rhoici - ki_o_rhoici*dhsnowdt(isi)*strg_2 )  & 
                     & /(csice*hice_p(isi))  
          ! Update the ice surface temperature
          tice_n(isi) = tice_p(isi) + dtime*dticedt(isi)

          IF ( lis_coupled_run ) THEN
            ! Coupling to icon-o: 
            ! Heat flux (phiipr0 is computed above) 
            condhf(isi) = ki*phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)
            ! No change of ice thickness
            hice_n(isi) = hice_p(isi)
          ELSE
            ! Compute the time-rate-of-change of the ice thickness (note the sign of heat fluxes)
            dhicedt(isi) = -ki_o_rhoialf*dhsnowdt(isi) + qwat*r_rhoialf
            ! Update the ice thickness
            hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)
          END IF

        END IF 

      END IF FreezingMeltingRegime 

      ! Remove too thin ice or impose security constraints
      IF( hice_n(isi)<hice_min ) THEN
        ! Remove too thin ice
        hice_n(isi) = 0._wp
        ! Set the ice surface temperature equal to the fresh water freezing point
        tice_n(isi) = tf_fresh
      ELSE 
        ! Limit the ice thickness from above (security)
        hice_n(isi) = MIN(hice_n(isi), hice_max)
        ! Limit the ice surface temperature from above (security)
        tice_n(isi) = MIN(tice_n(isi), tf_fresh)
      END IF 

      ! Currently, snow over sea ice is not treated explicitly
      ! Set the snow thickness to zero 
      hsnow_n(isi) = 0._wp
      ! Set the snow surface temperature equal to the ice surface temperature
      tsnow_n(isi) = tice_n(isi)

      ! Compute tendencies (for eventual use outside the sea-ice scheme program units)
      dticedt(isi)  = (tice_n(isi)-tice_p(isi))*r_dtime
      dhicedt(isi)  = (hice_n(isi)-hice_p(isi))*r_dtime
      dtsnowdt(isi) = (tsnow_n(isi)-tsnow_p(isi))*r_dtime
      dhsnowdt(isi) = (hsnow_n(isi)-hsnow_p(isi))*r_dtime

    END DO GridBoxesWithSeaIce

    ! Compute sea-ice albedo through a rate equation  
    PrognosticSeaIceAlbedo: IF ( lprog_albsi ) THEN

      ! Loop over grid boxes where sea ice is present
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(albsi_e, albsi_snow_e, taualbsi, rtaualbsisn, albsi_e_wghtd)
      DO isi=1, nsigb

        ! Equilibrium sea-ice albedo (function of sea-ice surface temperature)
        albsi_e = alb_seaice_equil( tice_n(isi) )

        ! Equilibrium albedo of snow over sea ice (function of sea-ice surface temperature)
        albsi_snow_e = albsi_snow_max * ( 1.0_wp - c1_albsi_snow                   &
          &                           * EXP(-c2_albsi_snow*(tf_fresh-tice_n(isi))) )

        ! Relaxation time scale for sea-ice albedo
        ! Interpolate linearly between maximum and minimum relaxation time scales
        ! over a given temperaure range
        taualbsi = taualbsi_max + (taualbsi_min-taualbsi_max)  &
                 * ((tice_n(isi)-t_taualbsi_min)*rdelt_taualbsi)
        ! Limit relaxation time scale from below and from above
        taualbsi = MIN(taualbsi_max,MAX(taualbsi,taualbsi_min))
        ! Use temperature-dependent relaxation time scale
        ! if sea-ice albedo tends to decrease,
        ! and a maximum time scale otherwise
        taualbsi = MERGE( taualbsi, taualbsi_max, (albsi_p(isi)>albsi_e) )

        ! Reciprocal of the relaxation time scale for snow-over-ice albedo
        ! Relaxation towards snow-over-ice albedo is only applied if
        ! albedo tends to increase, and
        ! the sea-ice surface temperature is not too close to the freezing point
        rtaualbsisn = snow_rate(isi)*c_tausi_snow*MERGE( 1._wp, 0._wp,  &
          &           ((albsi_p(isi)<albsi_snow_e).AND.(tice_n(isi)<t_albsi_snow_max)) )

        ! Weighted equilibrium albedo 
        albsi_e_wghtd = (albsi_e+taualbsi*rtaualbsisn*albsi_snow_e)  &
          &           / (1._wp+taualbsi*rtaualbsisn)

        ! Relax sea-ice albedo towards equilibrium value 
        albsi_n(isi) = albsi_e_wghtd+(albsi_p(isi)-albsi_e_wghtd)  &
          &          * EXP(-dtime*(1._wp/taualbsi+rtaualbsisn))

      END DO

    ENDIF PrognosticSeaIceAlbedo
    !$ACC END PARALLEL

    ! Store time tendencies (optional)
    IF (PRESENT(opt_dticedt)) THEN
#ifdef _OPENACC
        CALL finish ('seaice_timestep_nwp', 'OpenACC version currently does not support the optional argument opt_dticedt')
#endif
      opt_dticedt(1:nsigb)  = dticedt(1:nsigb)
      IF (nsigb < SIZE(opt_dticedt)) THEN
        opt_dticedt(nsigb+1:) = 0._wp 
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhicedt)) THEN
#ifdef _OPENACC
        CALL finish ('seaice_timestep_nwp', 'OpenACC version currently does not support the optional argument opt_dhicedt')
#endif
      opt_dhicedt(1:nsigb)  = dhicedt(1:nsigb)
      IF (nsigb < SIZE(opt_dhicedt)) THEN
        opt_dhicedt(nsigb+1:) = 0._wp
      ENDIF
    ENDIF
    IF (PRESENT(opt_dtsnowdt)) THEN
#ifdef _OPENACC
        CALL finish ('seaice_timestep_nwp', 'OpenACC version currently does not support the optional argument opt_dtsnowdt')
#endif
      opt_dtsnowdt(1:nsigb) = dtsnowdt(1:nsigb)
      IF (nsigb < SIZE(opt_dtsnowdt)) THEN
        opt_dtsnowdt(nsigb+1:)= 0._wp
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhsnowdt)) THEN
#ifdef _OPENACC
        CALL finish ('seaice_timestep_nwp', 'OpenACC version currently does not support the optional argument opt_dhsnowdt')
#endif
      opt_dhsnowdt(1:nsigb) = dhsnowdt(1:nsigb)
      IF (nsigb < SIZE(opt_dhsnowdt)) THEN
        opt_dhsnowdt(nsigb+1:)= 0._wp
      ENDIF
    ENDIF
 
    IF (.NOT. using_cuda_graph) THEN
      !$ACC WAIT(1)
    END IF
    !$ACC END DATA
    !$ACC END DATA
    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_timestep_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Coldstart for sea-ice parameterization scheme. 
  !!
  !! Coldstart for sea-ice parameterization scheme. Sea-ice surface temperature and sea-ice 
  !! thickness are initialized with meaningful values.
  !! Note that an estimate of the sea-ice temperature is required for the cold start and is 
  !! assumed to be available. The only option at the time being is to use the IFS skin 
  !! temperature for the cold start initialization of t_ice. Since an estimate of the 
  !! ice thickness h_ice is generally not available, h_ice is initialized with a 
  !! meaningful constant value (1m). 
  !! Note that only "sea" grid boxes are initialized;
  !! the "lake" grid boxes are left intact. 
  !! 

  SUBROUTINE seaice_coldinit_nwp (                              & 
                          &  nswgb,                             &
                          &  frice_thrhld,                      &
                          &  frsi,                              &
                          &  temp_in,                           &
                          &  tice_p, hice_p, tsnow_p, hsnow_p,  &
                          &  albsi_p,                           &
                          &  tice_n, hice_n, tsnow_n, hsnow_n,  &
                          &  albsi_n                            &
                          &  )

    IMPLICIT NONE

    ! Procedure arguments 

    INTEGER, INTENT(IN)  ::                &
                        &  nswgb      !< number of "sea" grid boxes within a block (<=nproma)

    REAL(wp), INTENT(IN) :: frice_thrhld     !< fraction threshold for creating a sea grid point

    REAL(wp), DIMENSION(:), INTENT(IN)    ::           &
                                          &  frsi    , &  !< sea-ice fraction [-]
                                          &  temp_in      !< meaningfull guess of ice surface temperature [K] 
                                                          !  e.g. tskin from IFS 


    REAL(wp), DIMENSION(:), INTENT(INOUT) ::           &
                                          &  tice_p  , &  !< temperature of ice upper surface at previous time level [K] 
                                          &  hice_p  , &  !< ice thickness at previous time level [m] 
                                          &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K] 
                                          &  hsnow_p , &  !< snow thickness at previous time level [m] 
                                          &  albsi_p , &  !< sea-ice albedo at previous time level [-] 
                                          &  tice_n  , &  !< temperature of ice upper surface at new time level [K] 
                                          &  hice_n  , &  !< ice thickness at new time level [m] 
                                          &  tsnow_n , &  !< temperature of snow upper surface at new time level [K] 
                                          &  hsnow_n , &  !< snow thickness at new time level [m] 
                                          &  albsi_n      !< sea-ice albedo at new time level [-] 

    ! Local variables 

    INTEGER ::      &
            &  isi  !< DO loop index


    REAL(wp), PARAMETER :: h_ice_coldstart = 1.0_wp   ! sea-ice thickness for cold start [m]

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Loop over all grid boxes
    DO isi=1, nswgb

      ! Note that we make use of >= instead of > in order to be consistent 
      ! with the seaice index list generation routine
      IF ( frsi(isi) >= frice_thrhld ) THEN  ! ice point

        hice_p(isi)  = h_ice_coldstart             ! constant ice thickness of 1m 
        tice_p(isi)  = temp_in(isi)                ! some proper estimate (here: tskin from IFS)
        tice_p(isi)  = MIN(tice_p(isi), tf_fresh)  ! security
        tsnow_p(isi) = tice_p(isi)                ! snow temperature is equal to ice temperature
        hsnow_p(isi) = 0._wp                      ! snow over ice is not treated explicitly 
        IF ( lprog_albsi ) THEN                   ! set sea-ice albedo to its equilibrium value
          albsi_p(isi) = alb_seaice_equil( tice_p(isi) )
        ENDIF

        ! Set variables at new time level
        tice_n(isi)  = tice_p(isi)     
        hice_n(isi)  = hice_p(isi)     
        tsnow_n(isi) = tsnow_p(isi)    
        hsnow_n(isi) = hsnow_p(isi)
        IF ( lprog_albsi ) THEN                   ! set sea-ice albedo to its equilibrium value
          albsi_n(isi) = albsi_p(isi)
        ENDIF
      ENDIF

    END DO ! isi  

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_coldinit_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Cold start initialization of prognostic sea-ice albedo.
  !!
  !! This routine is used when 
  !! cold start initialization of prognostic sea-ice albedo is necessary,
  !! whereas the sea-ice surface temperature and the sea-ice thickness 
  !! should not be (re-)initilaized.
  !! This occurs, for example, if the sea-ice scheme has already been used,
  !! but the sea-ice albedo from previous runs is not available.
  !! The sea-ice albedo is initialized with its equilibrium value 
  !! that is a function of sea-ice surface temperature.
  !! Note that only "sea" grid boxes are initialized;
  !! the "lake" grid boxes are left intact. 
  !!

  SUBROUTINE seaice_coldinit_albsi_nwp (              &
                                    &  nswgb,         &
                                    &  frice_thrhld,  &
                                    &  frsi,          &
                                    &  tice_p,        &
                                    &  albsi_p,       &
                                    &  albsi_n        &
                                    &  )

    IMPLICIT NONE

    ! Procedure arguments 

    INTEGER, INTENT(IN) ::          &
                        &  nswgb      !< number of "see" grid boxes within a block (<=nproma)

    REAL(wp), INTENT(IN) :: frice_thrhld     !< fraction threshold for creating a sea grid point

    REAL(wp), DIMENSION(:), INTENT(IN)    ::           &
                                          &  frsi    , &  !< sea-ice fraction [-]
                                          &  tice_p       !< temperature of ice upper surface at previous time level [K] 


    REAL(wp), DIMENSION(:), INTENT(INOUT) ::           &
                                          &  albsi_p , &  !< sea-ice albedo at previous time level [-] 
                                          &  albsi_n      !< sea-ice albedo at new time level [-] 
    ! Local variables 

    INTEGER ::      &
            &  isi  !< DO loop index

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Loop over sea-water grid boxes
    DO isi=1, nswgb

      ! Note that we make use of >= instead of > in order to be consistent 
      ! with the seaice index list generation routine
      IF ( frsi(isi) >= frice_thrhld ) THEN  ! ice point

        ! set sea-ice albedo to its equilibrium value
        albsi_p(isi) = alb_seaice_equil( tice_p(isi) )

        ! set albedo at new time level
        albsi_n(isi) = albsi_p(isi)

      ENDIF

    END DO ! isi  

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_coldinit_albsi_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Equilibrium sea-ice albedo is computed 
  !! as function of the sea-ice surface temperature.
  !!

  REAL (wp) FUNCTION alb_seaice_equil ( t_ice )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    ! Procedure arguments 

    REAL(wp), INTENT(IN) :: t_ice       !< temperature of ice upper surface [K]


    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    alb_seaice_equil = csalb(ist_seaice) * ( 1.0_wp - 0.3143_wp * EXP(-0.35_wp*(tf_fresh-t_ice)) )

    ! A derived constant 0.35 is equal to 95.6/tf_fresh, where tf_fresh=273.15 K, 
    ! and has a dimensions of K^{-1}.
    ! A derived constant 0.3143 is equal to (albsi_max-albsi_min)/albsi_max, 
    ! where albsi_max=csalb(ist_seaice)=0.7 and albsi_min=0.48.

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END FUNCTION alb_seaice_equil

END MODULE sfc_seaice
