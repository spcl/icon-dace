! This module configures:
! - Upper-atmosphere physics
! - Upper-atmosphere extrapolation
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

MODULE mo_upatmo_config

  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: message, message_text, finish
  USE mo_impl_constants,           ONLY: max_dom,                    &
    &                                    MODE_IFSANA, MODE_COMBINED, &
    &                                    MODE_DWDANA, inoforcing,    &
    &                                    SUCCESS, inwp, iaes,        &
    &                                    inh_atmosphere, ivexpol
  USE mo_model_domain,             ONLY: t_patch
  USE mo_upatmo_impl_const,        ONLY: iUpatmoStat, imsg_thr, itmr_thr, &
    &                                    iUpatmoPrcStat, iUpatmoGasStat,  &
    &                                    iUpatmoGrpId,  iUpatmoExtdatStat
  USE mo_upatmo_phy_config,        ONLY: t_upatmo_aes_phy,         &
    &                                    t_upatmo_nwp_phy,         &
    &                                    t_upatmo_phy_config,      &
    &                                    configure_upatmo_physics, &
    &                                    print_config_upatmo_physics
  USE mo_upatmo_utils,             ONLY: init_logical_1d
  USE mo_name_list_output_types,   ONLY: t_output_name_list
  USE mo_name_list_output_config,  ONLY: is_variable_in_output
  USE mtime,                       ONLY: datetime
  USE mo_timer,                    ONLY: timers_level, timer_start, timer_stop,   &
    &                                    timer_upatmo_constr, timer_upatmo_destr, &
    &                                    timer_upatmo


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: upatmo_exp_config
  PUBLIC :: upatmo_phy_config
  PUBLIC :: t_upatmo_config
  PUBLIC :: upatmo_config
  PUBLIC :: configure_upatmo
  PUBLIC :: destruct_upatmo
  PUBLIC :: check_upatmo

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_config'

  ! Note: upper-atmosphere physics and upper-atmosphere extrapolation 
  ! are regarded as two elements of the upper-atmosphere extension of ICON. 
  ! This is why we treat them in one module (as well as code economy).

  !------------------------------------------------------------
  !                    Configuration types
  !------------------------------------------------------------

  ! Namelist parameters that control the vertical extrapolation 
  ! of the initial atmospheric state into the upper atmosphere.
  !
  TYPE t_upatmo_exp_config
    REAL(wp) :: expol_start_height    ! [m] Height above which extrapolation (blending) of initial data starts
    REAL(wp) :: expol_blending_scale  ! [m] Blending scale height
    REAL(wp) :: expol_vn_decay_scale  ! [m] Scale height for (exponential) decay of extrapolated 
                                      ! horizontal wind component (for stability reasons)
    REAL(wp) :: expol_temp_infty      ! [K] Climatological temperature of exosphere (for z -> infinity)
    LOGICAL  :: lexpol_sanitycheck    ! .TRUE. -> Apply sanity check to extrapolated fields
    ! Status
    LOGICAL :: lset = .FALSE.         ! .TRUE. after assignment of namelist entries
  END TYPE t_upatmo_exp_config

  !------------------------------------------------------------

  ! For reasons of easier generalizability 
  ! the variables of the above types become domain-dependent. 
  ! Since they are used, to store the namelist input, 
  ! i.e. at a point in the program sequence before the actual number of domains is known, 
  ! they have to be allocated for the maximum permissible domain range '0:max_dom'.
  ! A significant number of parameters are derived from the namelist settings 
  ! (especially in the context of the upper-atmosphere physics),  
  ! so they are not gathered in the above, but in the following type(s), 
  ! which can be allocated for the actual number of domains 'n_dom_start:n_dom', 
  ! in order to save some memory in standard simulations without upper-atmosphere.
  !
  TYPE t_exp
    LOGICAL :: l_expol       ! .TRUE. -> upper-atmosphere-specific vertical extrapolation 
                             ! of initial data towards climatological values is switched on    
    INTEGER :: nexpollev     ! Index of grid layer for and above which the extrapolation 
                             ! of initial data will take place
    LOGICAL :: l_initicon_config = .FALSE. ! .TRUE. -> first step of configuration of extrapolation 
                             ! after namelist read-in has taken place
    ! Status variables
    LOGICAL :: l_status(iUpatmoStat%nitem)
  END TYPE t_exp
  !
  ! Collector
  !
  TYPE t_upatmo_config
    TYPE(t_upatmo_aes_phy)   :: aes_phy
    TYPE(t_upatmo_nwp_phy)   :: nwp_phy
    TYPE(t_exp) :: exp
    ! Miscellaneous
    REAL(wp) :: dt_fastphy             ! On several occasions we need the domain-specific dtime 
                                       ! and 'atm_phy_nwp_config(jg)%dt_fastphy' may not be available
    REAL(wp) :: dt_dyn_nom             ! Domain-specific nominal dynamics time step
    ! Status variables
    ! (Although desirable, the creation 
    ! of a status object is currently too complicated 
    ! and not worth the effort.)
    LOGICAL :: l_status(iUpatmoStat%nitem)
  END TYPE t_upatmo_config

  !------------------------------------------------------------

  TYPE(t_upatmo_exp_config), TARGET              :: upatmo_exp_config(0:max_dom)
  TYPE(t_upatmo_phy_config), TARGET              :: upatmo_phy_config(0:max_dom)
  TYPE(t_upatmo_config),     TARGET, ALLOCATABLE :: upatmo_config(:) 

  !------------------------------------------------------------

CONTAINS !..................................................................................

  !------------------------------------------------------------
  !                       Subroutines
  !------------------------------------------------------------

  !>
  !! This subroutine configures:
  !! - Upper-atmosphere physics
  !! - Upper-atmosphere extrapolation
  !!
  !!
  SUBROUTINE configure_upatmo( n_dom_start,            & !in
    &                          n_dom,                  & !in
    &                          p_patch,                & !in
    &                          lrestart,               & !in
    &                          lupatmo_phy,            & !in
    &                          init_mode,              & !in
    &                          iforcing,               & !in
    &                          tc_exp_startdate,       & !in
    &                          tc_exp_stopdate,        & !in
    &                          start_time,             & !in
    &                          end_time,               & !in
    &                          dtime,                  & !in
    &                          dt_rad_nwp,             & !in
    &                          ndyn_substeps,          & !in
    &                          flat_height,            & !in
    &                          l_orbvsop87,            & !in
    &                          cecc,                   & !in
    &                          cobld,                  & !in
    &                          clonp,                  & !in
    &                          lyr_perp,               & !in
    &                          yr_perp,                & !in
    &                          model_base_dir,         & !in
    &                          msg_level,              & !in
    &                          vct_a                   ) !(opt)in

    ! In/out variables
    INTEGER,            INTENT(IN) :: n_dom_start            ! Start index of domains
    INTEGER,            INTENT(IN) :: n_dom                  ! End index of domains
    TYPE(t_patch),      INTENT(IN) :: p_patch(n_dom_start:)  ! Domain properties
    LOGICAL,            INTENT(IN) :: lrestart               ! Switch for restart mode
    LOGICAL,            INTENT(IN) :: lupatmo_phy(:)         ! (max_dom) Switch for upper-atmosphere physics in nwp-mode
    INTEGER,            INTENT(IN) :: init_mode              ! Initialization mode
    INTEGER,            INTENT(IN) :: iforcing               ! Switch for physics package (NWP, AES etc.)
    TYPE(datetime),     INTENT(IN) :: tc_exp_startdate       ! Experiment start date
    TYPE(datetime),     INTENT(IN) :: tc_exp_stopdate        ! Experiment end date
    REAL(wp),           INTENT(IN) :: start_time(:)          ! (max_dom) Time at which execution of domain starts
    REAL(wp),           INTENT(IN) :: end_time(:)            ! (max_dom) Time at which execution of domain ends
    REAL(wp),           INTENT(IN) :: dtime                  ! Fast-physics/advective time step for primary domain
    REAL(wp),           INTENT(IN) :: dt_rad_nwp(:)          ! Tendency update period for radiation under NWP forcing
    INTEGER,            INTENT(IN) :: ndyn_substeps          ! Number of dynamics' substeps per fast-physics time step
    REAL(wp),           INTENT(IN) :: flat_height            ! Below 'flat_height' grid layer interfaces follow topography
    LOGICAL,            INTENT(IN) :: l_orbvsop87(:)         ! (max_dom) .TRUE. for VSOP87 orbit, .FALSE. for Kepler orbit
    REAL(wp),           INTENT(IN) :: cecc(:)                ! (max_dom) Eccentricity  of  orbit
    REAL(wp),           INTENT(IN) :: cobld(:)               ! (max_dom) Obliquity of Earth axis
    REAL(wp),           INTENT(IN) :: clonp(:)               ! (max_dom) Long. of the perihelion
    LOGICAL,            INTENT(IN) :: lyr_perp(:)            ! (max_dom) Switch for perpetuation of Earth orbit
    INTEGER,            INTENT(IN) :: yr_perp(:)             ! (max_dom) Year, for which Earth orbit is perpetuated
    CHARACTER(LEN=*),   INTENT(IN) :: model_base_dir         ! Path for input files
    INTEGER,            INTENT(IN) :: msg_level              ! Message level
    REAL(wp), OPTIONAL, INTENT(IN) :: vct_a(:)               ! (nlev+1) Nominal heights of grid layer interfaces

    ! Local variables
    REAL(wp) :: dt_grp_prevdom(iUpatmoGrpId%nitem)
    INTEGER  :: jg, jg_ordered, jg_ref, jg_aux
    INTEGER  :: nlev, nlevp1, nshift_total, n_dom_shift
    LOGICAL  :: l_upatmo_phy, ltimer
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':configure_upatmo'

    !---------------------------------------------------------

    ltimer   = timers_level > itmr_thr%med
    IF (ltimer) THEN
      CALL timer_start(timer_upatmo)
      CALL timer_start(timer_upatmo_constr)
    ENDIF

    !-----------------------------------------------------
    !             Domain-independent checks
    !-----------------------------------------------------

    ! (We assume that the setup of 'p_patch' took already place.)

    IF (.NOT. ANY((/0, 1/) == n_dom_start)) THEN
      ! For domain-dependent fields the max. range of allocation 
      ! with which we can reckon is '0:max_dom'
      !                              -
      CALL finish(routine, 'Something has changed regarding n_dom_start.' )
    ELSEIF (.NOT. PRESENT(vct_a)) THEN
      ! For the setup of the upper-atmosphere extrapolation, we need 'vct_a', 
      ! and we take its absence as an indicator that it has not been allocated yet
      CALL finish(routine, 'vct_a still uninitialized.')
    ENDIF

    ! Some of the settings below (especially for the upper-atmosphere extrapolation) 
    ! are more conveniently done, if the domain loop follows the order: 
    ! jg = 1, 2, 3, ..., n_dom, n_dom_start, 
    ! provided that 'n_dom_start = 0'. 
    n_dom_shift = 1 - n_dom_start

    ! Domain loop
    DO jg_ordered = n_dom_start, n_dom

      ! First, 'jg = 1' to 'n_dom' ... 
      jg = jg_ordered + n_dom_shift
      ! ... and then 'jg = 0'
      IF (jg > n_dom) jg = n_dom + 1 - jg

      ! Number of grid layers
      nlev   = p_patch(jg)%nlev
      nlevp1 = p_patch(jg)%nlevp1
      
      ! Shift of grid layer index, 
      ! to account for vertical nesting
      IF (jg > 0) THEN
        nshift_total = p_patch(jg)%nshift_total
      ELSE
        ! 'nshift_total' is not initialized for 'jg = n_dom_start', 
        ! if 'n_dom_start = 0', so we take the value from the primary domain 
        nshift_total = p_patch(1)%nshift_total
      ENDIF

      !-----------------------------------------------------
      !               Domain-dependent checks
      !-----------------------------------------------------
      
      IF ((jg < n_dom_start) .OR. (jg > n_dom)) THEN
        ! Some rudimentary checks of the reordering 
        ! of the domain sequence may be in order
        CALL finish(routine, "Domain index jg has unexpected value.")
      ELSEIF ((jg_ordered == n_dom_start) .AND. (jg /= 1)) THEN
        CALL finish(routine, "First domain to be configured has to be jg=1.")
      ELSEIF (.NOT. upatmo_config(jg)%l_status(iUpatmoStat%checked)) THEN
        ! (Just to make sure. Actually it should have been set to .true. 
        ! right after the allocation of 'upatmo_config')
        CALL finish(routine, "Check calling sequence: check_upatmo -> configure_upatmo.")
      ELSEIF (.NOT. upatmo_config(jg)%exp%l_initicon_config) THEN
        ! For the final configuration of the upper-atmosphere extrapolation it is required 
        ! that the preliminary configuration in 'src/configure_model/mo_initicon_config' has taken place
        CALL finish(routine, "Check calling sequence: configure_initicon -> configure_upatmo.")
      ENDIF

      !-----------------------------------------------------
      !                   Configuration
      !-----------------------------------------------------

      ! On several occasions we need the domain-specific value of dtime, 
      ! but 'iforcing /= inwp', so that 'atm_phy_nwp_config(jg)%dt_fastphy' is not available. 
      ! So we compute it here following its computation in 
      ! 'src/configure_model/mo_atm_phy_nwp_config: configure_atm_phy_nwp'. 
      ! A potential rescaling of 'dtime' with 'src/configure_model/mo_grid_config: grid_rescale_factor' 
      ! took already place in 'src/shared/mo_time_management: compute_timestep_settings', 
      ! which is called in 'src/configure_model/mo_nml_crosscheck: atm_crosscheck'.
      upatmo_config(jg)%dt_fastphy = dtime / 2._wp**(p_patch(jg)%level-p_patch(1)%level)
      
      ! We take the opportunity, to compute the nominal dynamics time step
      ! ("nominal", because 'ndyn_substeps' may change its value during runtime)
      upatmo_config(jg)%dt_dyn_nom = upatmo_config(jg)%dt_fastphy / REAL(ndyn_substeps, wp)

      !---------------
      !    Physics
      !---------------

      IF (jg > 0) THEN 
        l_upatmo_phy = lupatmo_phy(jg)
        jg_aux       = jg
      ELSE
        ! Upper-atmosphere physics are not required 
        ! on the radiation grid 'jg = 0'.
        ! In addition, there is no entry available, since 'lupatmo_phy(1:max)'. 
        l_upatmo_phy = .FALSE.
        ! Likewise, input for 'l_orbvsop87', 'cecc', 'cobld', 'clonp', 'lyr_perp', 'yr_perp', 
        ! 'start_time' and 'end_time' is only available for domains >= 1, so for 'jg = 0' 
        ! we simply take the values from the primary domain (that should do no harm)
        jg_aux = 1
      ENDIF

      ! Initialize update period for physics tendencies on previous domain for domain 1 (and 0). 
      ! It is such that it should have no effect on the update period modification in 
      ! 'src/upper_atmosphere/mo_upatmo_phy_config: configure_nwp_event' for domain 1.
      ! (For the subsequent domains 'configure_upatmo_physics' will return 
      ! the corresponding updates of 'dt_grp_prevdom'.)
      IF (jg < 2) THEN
        dt_grp_prevdom(:) = (2._wp + AINT(upatmo_phy_config(jg)%nwp_grp(:)%dt / upatmo_config(jg)%dt_fastphy)) * &
          &                 upatmo_config(jg)%dt_fastphy
      ENDIF

      CALL configure_upatmo_physics( jg                      = jg,                           & !in
        &                            lupatmo_phy             = l_upatmo_phy,                 & !in
        &                            iforcing                = iforcing,                     & !in
        &                            l_orbvsop87             = l_orbvsop87(jg_aux),          & !in
        &                            cecc                    = cecc(jg_aux),                 & !in
        &                            cobld                   = cobld(jg_aux),                & !in
        &                            clonp                   = clonp(jg_aux),                & !in
        &                            lyr_perp                = lyr_perp(jg_aux),             & !in
        &                            yr_perp                 = yr_perp(jg_aux),              & !in
        &                            nlev                    = nlev,                         & !in
        &                            nshift_total            = nshift_total,                 & !in
        &                            tc_exp_startdate        = tc_exp_startdate,             & !in
        &                            tc_exp_stopdate         = tc_exp_stopdate,              & !in
        &                            start_time              = start_time(jg_aux),           & !in
        &                            end_time                = end_time(jg_aux),             & !in
        &                            dtime                   = dtime,                        & !in
        &                            dt_fastphy              = upatmo_config(jg)%dt_fastphy, & !in
        &                            dt_rad_nwp              = dt_rad_nwp(jg_aux),           & !in
        &                            dt_grp_prevdom          = dt_grp_prevdom,               & !inout
        &                            model_base_dir          = model_base_dir,               & !in
        &                            msg_level               = msg_level,                    & !in
        &                            timers_level            = timers_level,                 & !in
        &                            upatmo_phy_config       = upatmo_phy_config(jg),        & !inout                     
        &                            upatmo_aes_phy_config   = upatmo_config(jg)%aes_phy,    & !inout
        &                            upatmo_nwp_phy_config   = upatmo_config(jg)%nwp_phy,    & !inout
        &                            vct_a                   = vct_a                         ) !(opt)in

      !---------------
      ! Extrapolation
      !---------------
      
      ! The primary domain is our reference
      ! (-> One reason for the reordering of the domain loop)
      jg_ref = 1

      ! Because of the consistency check within 'configure_upatmo_extrapolation', 
      ! the initialization of 'l_status' has to be done here
      CALL init_logical_1d( variable=upatmo_config(jg)%exp%l_status, value=.FALSE., &
        &                   opt_ilist=(/iUpatmoStat%checked/), opt_mask="list"      )

      IF (jg == jg_ref) THEN
        CALL configure_upatmo_extrapolation( jg                = jg,                    & !in
          &                                  jg_ref            = jg_ref,                & !in
          &                                  nlev              = nlev,                  & !in
          &                                  nlevp1            = nlevp1,                & !in
          &                                  nshift_total      = nshift_total,          & !in
          &                                  flat_height       = flat_height,           & !in
          &                                  vct_a             = vct_a,                 & !(opt)in
          &                                  upatmo_exp_config = upatmo_exp_config(jg), & !in
          &                                  upatmo_config     = upatmo_config(jg)      ) !inout
      ELSE
        CALL configure_upatmo_extrapolation( jg                    = jg,                    & !in
          &                                  jg_ref                = jg_ref,                & !in
          &                                  nlev                  = nlev,                  & !in
          &                                  nlevp1                = nlevp1,                & !in
          &                                  nshift_total          = nshift_total,          & !in
          &                                  flat_height           = flat_height,           & !in
          &                                  vct_a                 = vct_a,                 & !(opt)in
          &                                  upatmo_exp_config     = upatmo_exp_config(jg), & !in
          &                                  upatmo_config         = upatmo_config(jg),     & !inout
          &                                  opt_upatmo_config_ref = upatmo_config(jg_ref)  ) !optin
      ENDIF

      !---------------
      ! Miscellaneous
      !---------------

      ! General purpose message output desired?
      upatmo_config(jg)%l_status(iUpatmoStat%message) = msg_level >= imsg_thr%high

      ! Timer monitoring desired?
      upatmo_config(jg)%l_status(iUpatmoStat%timer) = timers_level > itmr_thr%med
      
      ! 'upatmo_config' is allocated in any case, but not necessarily required
      upatmo_config(jg)%l_status(iUpatmoStat%required) = upatmo_config(jg)%aes_phy%l_status(iUpatmoStat%required) .OR. &
        &                                                upatmo_config(jg)%nwp_phy%l_status(iUpatmoStat%required) .OR. &
        &                                                upatmo_config(jg)%exp%l_status(iUpatmoStat%required)
      
      ! Indicate that configuration has taken place
      upatmo_config(jg)%l_status(iUpatmoStat%configured) = upatmo_config(jg)%aes_phy%l_status(iUpatmoStat%configured) .OR. &
        &                                                  upatmo_config(jg)%nwp_phy%l_status(iUpatmoStat%configured) .OR. &
        &                                                  upatmo_config(jg)%exp%l_status(iUpatmoStat%configured)

      !---------------
      !   Messages
      !---------------
      
      IF (upatmo_config(jg)%l_status(iUpatmoStat%required)) THEN
        CALL print_config( jg,                    & !in
          &                n_dom,                 & !in
          &                iforcing,              & !in
          &                nshift_total,          & !in
          &                lrestart,              & !in
          &                l_upatmo_phy,          & !in
          &                msg_level,             & !in
          &                upatmo_phy_config(jg), & !in
          &                upatmo_config(jg),     & !in
          &                vct_a                  ) !(opt)in
      ENDIF
      
    ENDDO  !jg_ordered

    IF (ltimer) THEN
      CALL timer_stop(timer_upatmo_constr)
      CALL timer_stop(timer_upatmo)
    ENDIF
    
  END SUBROUTINE configure_upatmo

  !====================================================================================

  !>
  !! Configure upper-atmosphere extrapolation.
  !!
  SUBROUTINE configure_upatmo_extrapolation( jg,                   & !in
    &                                        jg_ref,               & !in
    &                                        nlev,                 & !in
    &                                        nlevp1,               & !in
    &                                        nshift_total,         & !in
    &                                        flat_height,          & !in
    &                                        vct_a,                & !in
    &                                        upatmo_exp_config,    & !in
    &                                        upatmo_config,        & !inout
    &                                        opt_upatmo_config_ref ) !optin

    ! In/out variables
    INTEGER,                             INTENT(IN)    :: jg, jg_ref            ! Domain indices
    INTEGER,                             INTENT(IN)    :: nlev, nlevp1          ! Number of grid layers and interfaces
    INTEGER,                             INTENT(IN)    :: nshift_total          ! Index shift for vertical nesting
    REAL(wp),                            INTENT(IN)    :: flat_height           ! Below 'flat_height' grid layer interfaces 
                                                                                ! follow topography
    REAL(wp),                  OPTIONAL, INTENT(IN)    :: vct_a(:)              ! Nominal heights of grid layer interfaces
    TYPE(t_upatmo_exp_config),           INTENT(IN)    :: upatmo_exp_config     ! Namelist parameters
    TYPE(t_upatmo_config),               INTENT(INOUT) :: upatmo_config         ! Upper-atmosphere configuration
    TYPE(t_upatmo_config),     OPTIONAL, INTENT(IN)    :: opt_upatmo_config_ref ! Configuration of reference domain 

    ! Local variables
    INTEGER :: jk, jks
    LOGICAL :: l_found
    LOGICAL :: l_present_ref, l_configured_ref, l_expol_ref
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':configure_upatmo_extrapolation'

    !---------------------------------------------------------

    ! The initialization of the extrapolation switch 'l_expol' 
    ! takes place in 'src/configure_model/mo_initicon_config', 
    ! depending on 'initicon_nml: itype_vert_expol'.
    ! Here, the final setting takes place.

    !-----------------------------------------------------
    !                 Consistency check
    !-----------------------------------------------------

    ! Actually, vct_a is only optional, because it is already optional one level higher
    IF (.NOT. PRESENT(vct_a)) CALL finish(routine, 'vct_a has to be present.')

    IF (PRESENT(opt_upatmo_config_ref)) THEN
      l_present_ref    = .TRUE.
      l_configured_ref = opt_upatmo_config_ref%exp%l_status(iUpatmoStat%configured)
      l_expol_ref      = opt_upatmo_config_ref%exp%l_expol
    ELSE
      l_present_ref    = .FALSE.
      l_configured_ref = .FALSE.
      l_expol_ref      = .FALSE.
    ENDIF

    IF (l_present_ref .EQV. (jg == jg_ref)) THEN
      ! The reference configuration has to be present, 
      ! if the current domain is not the reference domain, but only then
      CALL finish(routine, "Presence/absence of optional reference has to coincide with jg/=jg_ref/jg=jg_ref.")
    ELSEIF (l_present_ref .AND. (.NOT. l_configured_ref)) THEN
      ! The reference domain should be the first one, 
      ! which enters this subroutine
      CALL finish(routine, "Optional reference is present, but not yet configured.")
    ENDIF

    !-----------------------------------------------------
    !         Initialization with default values
    !-----------------------------------------------------

    ! Please, do not set 'upatmo_config%exp%l_expol' here, 
    ! otherwise, you overwrite the assignment 
    ! in 'src/configure_model/mo_initicon_config: configure_initicon'
    upatmo_config%exp%nexpollev = -1         

    !-----------------------------------------------------
    !                   Configuration
    !-----------------------------------------------------

    IF ((jg /= jg_ref) .AND. (.NOT. l_expol_ref)) THEN
      ! If extrapolation has been switched off for domain 1 (reference),    
      ! it will be switched off for all other domains
      upatmo_config%exp%l_expol = .FALSE. 
    ELSEIF (upatmo_exp_config%expol_start_height < flat_height .OR.     &
      &     upatmo_exp_config%expol_start_height > vct_a(1+nshift_total)) THEN
      ! If extrapolation start heigth 'expol_start_height' is below 'flat_height' or
      ! above the domain top, no extrapolation of the upper-atmosphere type should be necessary
      upatmo_config%exp%l_expol = .FALSE. 
      ! In all other cases the initial setting of 'l_expol' remains
    ENDIF
    
    ! Determine vertical index of grid layer for and above which 
    ! the extrapolation of the upper-atmosphere type of initial data will take place 
    IF (upatmo_config%exp%l_expol .AND. (jg==jg_ref)) THEN
      l_found = .FALSE.
      DO jk=1, nlevp1
        jks = jk + nshift_total  ! (Not really necessary, but for completeness)
        IF (vct_a(jks) < upatmo_exp_config%expol_start_height) THEN
          ! The grid cell layer within which 'expol_start_height' is located 
          ! shall be the lowermost layer to which the upper-atmosphere extrapolation is applied, 
          ! and the first interface above 'expol_start_height' is the lowermost to which 
          ! the extrapolation is applied, i.e.'nexpollev' applies to full and half(!) levels
          upatmo_config%exp%nexpollev = jk - 1
          ! Indicate the find
          l_found = .TRUE.
          EXIT
        ENDIF
      ENDDO  !jk
      IF (.NOT. l_found) CALL finish(routine, 'Could not find nexpollev.')
      ! Just to make sure
      upatmo_config%exp%nexpollev = MIN( MAX( 1, upatmo_config%exp%nexpollev ), nlev )
    ELSEIF (upatmo_config%exp%l_expol .AND. (jg/=jg_ref)) THEN
      ! In the current configuration of the upper-atmosphere namelist, 
      ! no domain-wise specification of the extrapolation parameters is possible. 
      ! So the settings for jg /= jg_ref should be derivable 
      ! from the setting for the primary domain
      upatmo_config%exp%nexpollev = opt_upatmo_config_ref%exp%nexpollev - nshift_total
    ENDIF  !IF (upatmo_config%exp%l_expol .AND. (jg==jg_ref))
    
    ! 'upatmo_config' is allocated in any case, but not necessarily required
    upatmo_config%exp%l_status(iUpatmoStat%required) = upatmo_config%exp%l_expol
    
    ! Indicate that configuration has taken place
    upatmo_config%exp%l_status(iUpatmoStat%configured) = .TRUE.
    
  END SUBROUTINE configure_upatmo_extrapolation

  !====================================================================================

  !>
  !! Print message on upper-atmosphere configuration.
  !!
  SUBROUTINE print_config( jg,                & !in
    &                      n_dom,             & !in
    &                      iforcing,          & !in
    &                      nshift_total,      & !in
    &                      lrestart,          & !in
    &                      lupatmo_phy,       & !in
    &                      msg_level,         & !in
    &                      upatmo_phy_config, & !in
    &                      upatmo_config,     & !in
    &                      vct_a              ) !(opt)in

    ! In/out variables
    INTEGER,                         INTENT(IN)    :: jg                ! Domain index
    INTEGER,                         INTENT(IN)    :: n_dom             ! Number of domains
    INTEGER,                         INTENT(IN)    :: iforcing          ! Switch for physics package 
                                                                        ! (NWP, AES etc.)
    INTEGER,                         INTENT(IN)    :: nshift_total      ! Shift of vertical grid index 
    LOGICAL,                         INTENT(IN)    :: lrestart          ! Switch for restart mode
    LOGICAL,                         INTENT(IN)    :: lupatmo_phy       ! Switch for upper-atmosphere physics (NWP)
                                                                        ! for vertical nesting
    INTEGER,                         INTENT(IN)    :: msg_level         ! Message level
    TYPE(t_upatmo_phy_config),       INTENT(IN)    :: upatmo_phy_config ! Upper-atmosphere physics configuration 
                                                                        ! with namelist settings
    TYPE(t_upatmo_config),           INTENT(INOUT) :: upatmo_config     ! Upper-atmosphere configuration
    REAL(wp),              OPTIONAL, INTENT(IN)    :: vct_a(:)          ! Nominal heights of
                                                                        ! grid layer interfaces

    ! Local variables
    LOGICAL :: l_onlyPrimDom
    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':print_config'

    !---------------------------------------------------------

    ! Some print is only required for the primary domain
    l_onlyPrimDom = (jg==1)

    !-----------------------------------------------------
    !                "Any-case" output
    !-----------------------------------------------------

    ! Miscellaneous:
    IF (l_onlyPrimDom .AND. upatmo_config%l_status(iUpatmoStat%required)) THEN
      WRITE (message_text, '(a,i0,a)') &
        "(Info: most upper-atmosphere-related message output requires &
        &msg_level >= ", imsg_thr%high, ")"
      CALL message(routine, message_text)
      WRITE (message_text, '(a,i0,a)') &
        "(Info: most upper-atmosphere-related timers require timers_level >= ",&
        itmr_thr%med, ")"
      CALL message(routine, message_text)
    ENDIF

    IF (msg_level >= imsg_thr%low) THEN

      !-----------------------------------------------------
      !              "High-priority" output
      !-----------------------------------------------------

      IF (msg_level >= imsg_thr%med) THEN

        !-----------------------------------------------------
        !             "Medium-priority" output
        !-----------------------------------------------------
        
        IF (msg_level >= imsg_thr%high) THEN 

          !-----------------------------------------------------
          !               "Low-priority" output
          !-----------------------------------------------------

          ! Upper-atmosphere extrapolation:
          IF (upatmo_config%exp%l_status(iUpatmoStat%required) .AND. PRESENT(vct_a)) THEN
#define msg_prefix 'upatmo-expol(', jg, '): '
            WRITE (message_text, '(a,i0,2a,i0)') &
              msg_prefix, 'nexpollev: ', upatmo_config%exp%nexpollev
            CALL message(routine, message_text)
            WRITE (message_text, '(a,i0,2a,g32.5)') &
              msg_prefix, 'interface height above which extrapolation &
              &potentially takes place: ', &
              vct_a(upatmo_config%exp%nexpollev + nshift_total)
            CALL message(routine, message_text)
#undef msg_prefix
          ENDIF

        ENDIF  !imsg_thr%high
      ENDIF  !imsg_thr%med
    ENDIF  !imsg_thr%low

    ! Physics:
    CALL print_config_upatmo_physics( jg                      = jg,                      & !in
      &                               n_dom                   = n_dom,                   & !in
      &                               iforcing                = iforcing,                & !in
      &                               lrestart                = lrestart,                & !in
      &                               lupatmo_phy             = lupatmo_phy,             & !in
      &                               msg_level               = msg_level,               & !in
      &                               upatmo_phy_config       = upatmo_phy_config,       & !in
      &                               upatmo_aes_phy_config   = upatmo_config%aes_phy,   & !in
      &                               upatmo_nwp_phy_config   = upatmo_config%nwp_phy    ) !inout

  END SUBROUTINE print_config

  !====================================================================================

  !>
  !! Destruct upper-atmosphere configuration.
  !!
  SUBROUTINE destruct_upatmo( n_dom_start, & !in
    &                         n_dom        ) !in

    ! In/out variables
    INTEGER, INTENT(IN) :: n_dom_start ! Start index of domains
    INTEGER, INTENT(IN) :: n_dom       ! End index of domains

    ! Local variables
    INTEGER :: jg
    INTEGER :: istat
    LOGICAL :: ltimer
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':destruct_upatmo'

    !---------------------------------------------------------

    ltimer   = timers_level > itmr_thr%med
    IF (ltimer) THEN
      CALL timer_start(timer_upatmo)
      CALL timer_start(timer_upatmo_destr)
    ENDIF

    ! Deallocate upatmo_config
    ! (Allocated in 'src/namelists/mo_upatmo_nml: check_upatmo')
    IF (ALLOCATED(upatmo_config)) THEN
      ! Deallocate the event management variables, if required
      DO jg = n_dom_start, n_dom
        IF (upatmo_config(jg)%nwp_phy%l_phy_stat(iUpatmoPrcStat%enabled)) THEN
          CALL upatmo_config(jg)%nwp_phy%event_mgmt_grp%finalize()
        ENDIF
        IF (upatmo_config(jg)%nwp_phy%l_extdat_stat(iUpatmoExtdatStat%required)) THEN
          CALL upatmo_config(jg)%nwp_phy%event_mgmt_extdat%finalize()
        ENDIF
      ENDDO  !jg
      DEALLOCATE(upatmo_config, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, "Deallocation of upatmo_config failed.")
    ENDIF

    IF (ltimer) THEN
      CALL timer_stop(timer_upatmo_destr)
      CALL timer_stop(timer_upatmo)
    ENDIF

  END SUBROUTINE destruct_upatmo


  !====================================================================================

  !>
  !! Check for conflicts with other namelist settings.
  !!
  SUBROUTINE check_upatmo( n_dom_start,      & !in
    &                      n_dom,            & !in
    &                      iforcing,         & !in
    &                      lupatmo_phy,      & !in
    &                      l_limited_area,   & !in
    &                      ivctype,          & !in
    &                      flat_height,      & !in
    &                      itype_vert_expol, & !in
    &                      init_mode,        & !in
    &                      inwp_turb,        & !in
    &                      inwp_radiation)     !in

    ! In/out variables
    INTEGER,  INTENT(IN) :: n_dom_start       ! Start index of domains
    INTEGER,  INTENT(IN) :: n_dom             ! End index of domains
    INTEGER,  INTENT(IN) :: iforcing          ! Switch for physics package 
                                              ! (nwp, aes etc.)
    LOGICAL,  INTENT(IN) :: lupatmo_phy(:)    ! Switch for upper-atmosphere physics
                                              ! in nwp-mode
    LOGICAL,  INTENT(IN) :: l_limited_area    ! Switch for limited-area mode
    INTEGER,  INTENT(IN) :: ivctype           ! Type of vertical grid (SLEVE etc.)
    REAL(wp), INTENT(IN) :: flat_height       ! Below 'flat_height' grid layer 
                                              ! interfaces follow topography
    INTEGER,  INTENT(IN) :: itype_vert_expol  ! Type of vertical extrapolation 
                                              ! of initial atmosphere state
    INTEGER,  INTENT(IN) :: init_mode         ! Initialization mode 
    INTEGER,  INTENT(IN) :: inwp_turb(:)      ! Switch for turbulence scheme (NWP)
    INTEGER,  INTENT(IN) :: inwp_radiation(:) ! Switch for radiation scheme (NWP)

    ! Local variables
    INTEGER :: jg, istat
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':check_upatmo'

    !------------------------------------------------

    !---------------
    !    Physics
    !---------------

    IF (lupatmo_phy(1)) THEN 
      
      ! If upper-atmosphere physics are switched on for the NWP-mode ...
      
      IF (iforcing /= inwp) THEN
        ! ... the NWP-mode should be switched on
        WRITE (message_text, '(a,i0)') "nwp_phy_nml: lupatmo_phy only &
          &available, if run_nml: iforcing = inwp = ", inwp
        CALL finish(routine, message_text)
      ELSEIF (.NOT. ANY((/MODE_DWDANA, MODE_IFSANA, MODE_COMBINED/) == init_mode)) THEN
        ! ... only initialization with IFS or DWD analyses is allowed
        ! (it has not yet been figured out, how to include the upper-atmosphere physics into the IAU-infrastructure)
        WRITE (message_text, '(3(a,i0))') "NWP + upper-atmosphere physics &
          &exclusively allowed for MODE_DWDANA = ", MODE_DWDANA, &
          ", MODE_IFSANA = ", MODE_IFSANA, " and MODE_COMBINED = ",MODE_COMBINED
        CALL finish(routine, message_text)
      ENDIF

      ! Domain-wise checks
      DO jg = 1, n_dom
        ! A turbulence model has to be switched on (otherwise upper-atmosphere wind tendencies will not be accumulated)
        IF (lupatmo_phy(jg) .AND. inwp_turb(jg) == 0) THEN
          WRITE (message_text, '(a,i0)') &
            "Upper-atmosphere physics require inwp_turb > 0 in dom ", jg
          CALL finish(routine, message_text)
        ENDIF
        ! Some diagnostic variables, such as the cosine of the solar zenith angle, are computed by 
        ! the 'standard' radiation schemes, and they are required by the upper-atmosphere radiation schemes as well, 
        ! so we have to make sure that some standard scheme is switched on
        IF (lupatmo_phy(jg) .AND. inwp_radiation(jg) == 0) THEN
          WRITE (message_text, '(a,i0)') &
            "Upper-atmosphere physics require inwp_radiation > 0 in dom ", jg
          CALL finish(routine, message_text)
        ENDIF
      ENDDO  !jg

      ! ... OpenACC-parallelization is not available in combination with upper-atmosphere physics in NWP mode
#ifdef _OPENACC
      CALL finish(routine, "Upper-atmosphere physics (NWP) are not available in combination with Open-ACC.")
#endif

    ELSEIF (ANY(lupatmo_phy(2:max_dom))) THEN

      ! Just to be on the safe side:
      ! if 'lupatmo_phy' is switched off on dom 1, this should hold for all other doms > 1. 
      ! Skipping domains - for instance switch on domain 1, switch off domain 2, switch on domain 3 - 
      ! is not possible for the time being. This is because 'lupatmo_phy' has only two states, 
      ! so that domain skipping cannot be implemented into 'src/namelists/mo_nwp_phy_nml' 
      ! in a way comparable to the integer switches, such as 'inwp_satad'.
      CALL finish(routine, "Something is wrong with setting of lupatmo_phy in mo_nwp_phy_nml.")

    ENDIF    

    !---------------
    ! Extrapolation
    !---------------

    IF (itype_vert_expol == ivexpol%upatmo) THEN

      IF (ivctype == 2 .AND. &
        & ANY(upatmo_exp_config(:)%expol_start_height < flat_height)) THEN
        ! Upper-atmosphere extrapolation: start height above which extrapolation takes place 
        ! should not lie below 'flat_height'
        CALL finish(routine, &
          & "Upper-atmosphere extrapolation: start height has to be above flat_height.")
      ELSEIF (l_limited_area) THEN
        ! This type of extrapolation is not intended for the limited-area mode
        WRITE (message_text, '(a,i0)') &
          "The limited-area mode requires: itype_vert_expol = ", ivexpol%lin
        CALL finish(routine, message_text)
      ENDIF

    ENDIF  !IF (itype_vert_expol == ivexpol%upatmo)

    !---------------
    ! Miscellaneous
    !---------------

    IF (ALLOCATED(upatmo_config)) THEN
      CALL finish(routine, "Error in calling sequence: upatmo_config is already allocated.")
    ELSE
      ! Some of the variables in 'upatmo_config' might be necessary for a coarser radiation grid as well, 
      ! so the index range starts with 'n_dom_start' (which is zero, if a coarser radiation grid is used).
      ALLOCATE(upatmo_config(n_dom_start:n_dom), STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, "Allocation of upatmo_config failed.")
    ENDIF

    !---------------
    ! Status update
    !---------------
    
    ! Indicate that crosscheck took place
    upatmo_config(n_dom_start:n_dom)%aes_phy%l_status(iUpatmoStat%checked)   = .TRUE.
    upatmo_config(n_dom_start:n_dom)%nwp_phy%l_status(iUpatmoStat%checked)   = .TRUE.
    upatmo_config(n_dom_start:n_dom)%exp%l_status(iUpatmoStat%checked)       = .TRUE.
    upatmo_config(n_dom_start:n_dom)%l_status(iUpatmoStat%checked)           = .TRUE.

  END SUBROUTINE check_upatmo

END MODULE mo_upatmo_config
