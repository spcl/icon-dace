! This module supports src/configure_model/mo_upatmo_config
!
! In order to unburden mo_upatmo_config,
! we moved most of the type definitions and subroutines
! related to the upper-atmosphere physics here.
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

MODULE mo_upatmo_phy_config

  USE mo_kind,                     ONLY: wp, i8
  USE mo_exception,                ONLY: finish, message, message_text
  USE mo_impl_constants,           ONLY: MAX_CHAR_LENGTH, SUCCESS, &
    &                                    inwp, iaes
  USE mo_physical_constants,       ONLY: amd, amco2, amo2, amo3, amo, amno, amn2, avo, &
    &                                    cpd, cvd
  USE mo_math_constants,           ONLY: dbl_eps
  USE mo_upatmo_impl_const,        ONLY: iUpatmoStat, imsg_thr, iUpatmoPrcId, iUpatmoGrpId,  &
    &                                    iUpatmoGasId, iUpatmoExtdatId, iUpatmoTendId,       &
    &                                    iUpatmoPrcStat, iUpatmoGasStat, iUpatmoPrcMode,     &
    &                                    iUpatmoGasMode, iUpatmoExtdatStat, iorbit,          &
    &                                    startHeightDef, itmr_thr, iThermdynCoupling
  USE mo_upatmo_utils,             ONLY: init_logical_1d, is_variable_in_output_cond
  USE mo_util_string,              ONLY: real2string, &
    &                                    t_keyword_list, associate_keyword,       &
    &                                    with_keywords
  USE mtime,                       ONLY: MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, &
    &                                    datetime, timedelta, deallocateTimedelta,    &
    &                                    deallocateDatetime, getPTStringFromMS,       &
    &                                    newTimedelta, newDatetime,                   &
    &                                    getTotalMilliSecondsTimeDelta,               &
    &                                    OPERATOR(+), OPERATOR(-), OPERATOR(>)
  USE mo_phy_events,               ONLY: t_phyProcSlow, t_phyProcGroup
  USE mo_io_units,                 ONLY: filename_max
  USE mo_mpi,                      ONLY: my_process_is_stdio
  USE mo_name_list_output_types,   ONLY: t_output_name_list
  USE mo_grid_config,              ONLY: DEFAULT_ENDTIME
  USE mo_cdi,                      ONLY: FILETYPE_GRB, FILETYPE_GRB2
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: t_upatmo_aes_phy
  PUBLIC :: t_upatmo_nwp_phy
  PUBLIC :: t_upatmo_phy_config
  PUBLIC :: t_nwp_gas
  PUBLIC :: configure_upatmo_physics
  PUBLIC :: print_config_upatmo_physics

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_phy_config'

  REAL(wp), PARAMETER :: eps = 10._wp * ABS(dbl_eps)

  !------------------------------------------------------------
  !                  Configuration types
  !------------------------------------------------------------

  ! Types for the namelist entries:

  ! Control type for the groups in which the single upper-atmosphere physics processes are clustered, 
  ! following the example set by 'src/configre_model/mo_aes_phy_config' (only required for NWP-physics)  
  TYPE t_nwp_prc
    INTEGER                             :: imode         ! Selector for mode/type of group
    REAL(wp)                            :: dt            ! Time step for process group
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: t_start       ! Start time of process group
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: t_end         ! End time of process group
    REAL(wp)                            :: start_height  ! Start height of process group
  END TYPE t_nwp_prc

  !-----------------------------------------------------------------------

  ! Type for the radiatively active gase in the upper atmosphere
  TYPE t_nwp_gas
    INTEGER  :: imode   ! Gas mode
    REAL(wp) :: vmr     ! Volume mixing ratio ((m3/m3), should be equal to mole fraction (mol/mol))
    REAL(wp) :: mmr     ! Mass mixing ratio (kg/kg)
    REAL(wp) :: fscale  ! Scaling factor for mixing ratio ('frad_<gas>' in 'src/configre_model/mo_aes_rad_config')
  END TYPE t_nwp_gas

  !-----------------------------------------------------------------------

  ! Type for external data
  TYPE t_nwp_extdat
    REAL(wp)                    :: dt        ! Update period for time interpolation
    CHARACTER(LEN=filename_max) :: filename  ! Name of file containing external data
  END TYPE t_nwp_extdat

  !-----------------------------------------------------------------------

  ! The main typer for the namelist entries
  TYPE t_upatmo_phy_config
    INTEGER  :: orbit_type   ! Orbit model 
    INTEGER  :: solvar_type  ! Solar activity
    INTEGER  :: solvar_data  ! Solar activity data
    INTEGER  :: solcyc_type  ! Solar cycle 
    REAL(wp) :: cecc         ! Eccentricity of orbit
    REAL(wp) :: cobld        ! Obliquity of Earth axis
    REAL(wp) :: clonp        ! Longitude of perihelion
    LOGICAL  :: lyr_perp     ! Switch for perpetuation of Earth orbit for year 'yr_perp'
    INTEGER  :: yr_perp      ! Year, for which Earth orbit is perpetuated
    LOGICAL  :: lsanitycheck ! Switch for applying sanity checks
    ! AES-specific
    REAL(wp) :: aes_start_height(iUpatmoPrcId%nitem)         ! Start heights, above which 
                                                             ! processes compute tendencies
    ! NWP-specific
    TYPE(t_nwp_prc)    :: nwp_grp(iUpatmoGrpId%nitem)        ! Control for physics groups 
    TYPE(t_nwp_gas)    :: nwp_gas(iUpatmoGasId%nitem)        ! Radiatively active gases
    TYPE(t_nwp_extdat) :: nwp_extdat(iUpatmoExtdatId%nitem)  ! External data
    !
    INTEGER :: nwp_thermdyn_cpl         ! Type of thermodynamic coupling of physics & dynamics
    LOGICAL :: nwp_ldiss_from_heatdiff  ! Switch for considering heat source from heat diffusion
    ! Status
    LOGICAL :: lset = .FALSE.  ! .TRUE. after assignment of namelist entries
  END TYPE t_upatmo_phy_config

  !====================================================================================

  !--------------------------------

  TYPE t_upatmo_aes_phy
    REAL(wp) :: start_height(iUpatmoPrcId%nitem)  ! Start heights, above which 
                                                  ! processes compute tendencies
    REAL(wp) :: end_height(iUpatmoPrcId%nitem)    ! End heights, below which 
                                                  ! processes compute tendencies
    INTEGER  :: istartlev(iUpatmoPrcId%nitem)     ! Grid layer index corresponding to end_height(!)
    INTEGER  :: iendlev(iUpatmoPrcId%nitem)       ! Grid layer index corresponding to start_height(!)
    ! Status variables
    LOGICAL  :: l_enabled = .FALSE.               ! .TRUE.: upper-atmosphere physics are switched on
    LOGICAL  :: l_status(iUpatmoStat%nitem)
  END TYPE t_upatmo_aes_phy
 
  !==================================================================================== 

  ! Please note:
  ! * The NWP-specific configuration type contains a confusingly large number of variables,   
  !   but most of them were introduced to meet the requirements of NWP on memory 
  !   and runtime efficiency and neutrality at least to some extent.
  ! * A significant number of the variables is not domain-dependent, 
  !   so when the collector type 't_upatmo_nwp_phy' becomes part of 'upatmo_config(n_dom_start:n_dom)'
  !   in 'src/configure_model/mo_upatmo_config', we store redundant information. 
  !   However, for reasons of convenience, simplicity and coherence, we prefer to keep it this way.
  ! * One might get the impression that the flow control of the upper atmosphere,  
  !   currently covered by the following variables, would be in better hands at object-oriented programming.
  !   However, we think that the additional level of abstraction and complexity this involves 
  !   is only justified, if that flow control tool would be universal, 
  !   i.e. usable not only for the upper atmosphere, but for other elements of ICON as well (in principle).
  !   Creating such tool is beyond our means.
  ! * Currently, the flow control of the upper atmosphere assumes a "linear" program sequence. 
  !   It is not integrated into the IAU iterations. If this should once be desired, 
  !   the flow control needs a thorough revision! E.g., a proper reset of the status switches 'l_stat(:), ' 
  !   whose settings are called up at several places in the upper-atmosphere code, has to take place. 

  ! Process-related variables
  TYPE t_phy_prc
    INTEGER                        :: id           ! Process identifier
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: name         ! Process name
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: longname     ! Long process name
    INTEGER                        :: imode        ! Process mode
    INTEGER                        :: igrp         ! Identifier for group in which process is clustered
    REAL(wp)                       :: dt           ! Update period ("time step") for tendencies
    TYPE(datetime)                 :: t_start      ! Start time for parameterization
    TYPE(datetime)                 :: t_end        ! End time for parameterization
    TYPE(t_phyProcSlow)            :: event        ! Event management variable for single processes
                                                   ! (all upper-atmosphere physics are assumed to be 'slow' processes)
    REAL(wp)                       :: start_height ! Height above which tendencies are computed
    REAL(wp)                       :: end_height   ! Height below which tendencies are computed
    INTEGER                        :: istartlev    ! Grid layer index corresponding to end_height(!)
    INTEGER                        :: iendlev      ! Grid layer index corresponding to start_height(!)
    LOGICAL                        :: l_update(iUpatmoTendId%nitem)  ! Switches to indicate which vars will be updated
    LOGICAL                        :: l_gas(iUpatmoGasId%nitem)      ! Required gases
    LOGICAL                        :: l_stat(iUpatmoPrcStat%nitem)   ! Process status
    !
    LOGICAL, PRIVATE               :: lBeforeOpPhase ! Aux. switch for 'isBeforeOpPhase'
    LOGICAL, PRIVATE               :: lAfterOpPhase  ! Aux. Switch for 'isAfterOpPhase'
  CONTAINS
    !
    PROCEDURE :: isBeforeOpPhase => prc_isBeforeOpPhase
    PROCEDURE :: isAfterOpPhase  => prc_isAfterOpPhase
    PROCEDURE :: isInOpPhase     => prc_isInOpPhase
  END TYPE t_phy_prc

  !--------------------------------

  ! Gas-related variables
  TYPE t_phy_gas
    INTEGER                        :: id          ! Gas identifier
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: name        ! Gas name (e.g., 'co2')
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: longname    ! Long gas name (e.g., 'carbon dioxide')
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit        ! Gas unit (e.g., 'kg kg-1')
    INTEGER                        :: imode       ! Gas mode (corresponds to 'irad_<gas>' in 'mo_radiation_config')
    REAL(wp)                       :: vmr2mmr     ! Auxiliary factor to convert from volume mixing ratio to mass mixing ratio
    REAL(wp)                       :: mmr2vmr     ! Auxiliary factor to convert from mass mixing ratio to volume mixing ratio
    REAL(wp)                       :: mass2mol    ! Auxiliary factor to convert from mass to molecules
    LOGICAL                        :: l_prc(iUpatmoPrcId%nitem)     ! Processes making use of the gas
    LOGICAL                        :: l_grp(iUpatmoGrpId%nitem)     ! Process group --,,--
    LOGICAL                        :: l_stat(iUpatmoGasStat%nitem)  ! Gas status
  END TYPE t_phy_gas

  !--------------------------------

  ! Variables related to the external data
  TYPE t_phy_extdat
    INTEGER                        :: id       ! Identifier of external data type
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: name     ! Name of external data type (e.g., 'gases')
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: longname ! Long name (e.g., 'radiatively active gases')
    REAL(wp)                       :: dt       ! Update period for time interpolation of external data
    TYPE(t_phyProcSlow)            :: event    ! Event management variable
    LOGICAL                        :: l_stat(iUpatmoExtdatStat%nitem) ! Status
  END TYPE t_phy_extdat

  !--------------------------------

  TYPE t_upatmo_nwp_phy
    TYPE(t_phy_prc)    :: prc(iUpatmoPrcId%nitem)        ! Control units for single processes
    TYPE(t_phy_prc)    :: grp(iUpatmoGrpId%nitem)        ! --,,-- for process groups
    TYPE(t_phy_gas)    :: gas(iUpatmoGasId%nitem)        ! --,,-- for gases
    TYPE(t_phy_extdat) :: extdat(iUpatmoExtdatId%nitem)  ! --,,-- for external data
    !
    TYPE(t_phyProcGroup) :: event_mgmt_grp          ! Event management variable for all process groups
    TYPE(t_phyProcGroup) :: event_mgmt_extdat       ! --,,-- for all types of external data
    !
    LOGICAL :: l_phy_stat(iUpatmoPrcStat%nitem)        ! Total (summarizing) status of processes
    LOGICAL :: l_gas_stat(iUpatmoGasStat%nitem)        ! Total (summarizing) status of gases
    LOGICAL :: l_extdat_stat(iUpatmoExtdatStat%nitem)  ! Total (summarizing) status of external data
    LOGICAL :: l_any_update(iUpatmoTendId%nitem)       ! Switches to indicate which prog vars will be updated at all
    !
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: vname_prefix  ! Prefix for upper-atmosphere variable names 
                                                    ! (i.e., "vname_prefixVariable name", 
                                                    ! has to be entered to a varlist of 'output_nml', 
                                                    ! compare: 'src/atm_dyn_iconam/mo_nonhydro_state')
    !
    REAL(wp) :: thermdyn_cpl_fac                    ! Factor for thermodynamic coupling of physics 
    !
    LOGICAL :: l_status(iUpatmoStat%nitem)  ! Status variables
    !
    TYPE(datetime), PRIVATE :: t_start        ! Greatest lower bound of start times of groups
    TYPE(datetime), PRIVATE :: t_end          ! Least upper bound of end times of groups
    LOGICAL,        PRIVATE :: lBeforeOpPhase ! Aux. switch for 'isBeforeOpPhase'
    LOGICAL,        PRIVATE :: lAfterOpPhase  ! Aux. Switch for 'isAfterOpPhase'
  CONTAINS
    !
    PROCEDURE :: isBeforeOpPhase => grp_isBeforeOpPhase
    PROCEDURE :: isAfterOpPhase  => grp_isAfterOpPhase
    PROCEDURE :: isInOpPhase     => grp_isInOpPhase
  END TYPE t_upatmo_nwp_phy

CONTAINS !..................................................................................

  !>
  !! Configure (upper-atmosphere) physics.
  !!
  SUBROUTINE configure_upatmo_physics( jg,                      & !in
    &                                  lupatmo_phy,             & !in
    &                                  iforcing,                & !in
    &                                  l_orbvsop87,             & !in
    &                                  cecc,                    & !in
    &                                  cobld,                   & !in
    &                                  clonp,                   & !in
    &                                  lyr_perp,                & !in
    &                                  yr_perp,                 & !in
    &                                  nlev,                    & !in
    &                                  nshift_total,            & !in
    &                                  tc_exp_startdate,        & !in
    &                                  tc_exp_stopdate,         & !in
    &                                  start_time,              & !in
    &                                  end_time,                & !in
    &                                  dtime,                   & !in
    &                                  dt_fastphy,              & !in
    &                                  dt_rad_nwp,              & !in
    &                                  dt_grp_prevdom,          & !inout
    &                                  model_base_dir,          & !in
    &                                  msg_level,               & !in
    &                                  timers_level,            & !in
    &                                  upatmo_phy_config,       & !inout
    &                                  upatmo_aes_phy_config,   & !inout
    &                                  upatmo_nwp_phy_config,   & !inout
    &                                  vct_a                    ) !(opt)in

    ! In/out variables
    INTEGER,                   INTENT(IN)    :: jg                      ! Domain index
    LOGICAL,                   INTENT(IN)    :: lupatmo_phy             ! Switch for upper-atmosphere physics in nwp-mode
    INTEGER,                   INTENT(IN)    :: iforcing                ! Switch for physics package (NWP, AES etc.) 
    LOGICAL,                   INTENT(IN)    :: l_orbvsop87             ! .TRUE. for VSOP87 orbit, 
                                                                        ! .FALSE. for Kepler orbit
    REAL(wp),                  INTENT(IN)    :: cecc                    ! Eccentricity  of  orbit
    REAL(wp),                  INTENT(IN)    :: cobld                   ! Obliquity of Earth axis
    REAL(wp),                  INTENT(IN)    :: clonp                   ! Long. of the perihelion
    LOGICAL,                   INTENT(IN)    :: lyr_perp                ! Switch for perpetuation of Earth orbit
    INTEGER,                   INTENT(IN)    :: yr_perp                 ! Year, for which Earth orbit is perpetuated
    INTEGER,                   INTENT(IN)    :: nlev                    ! Number of vertical grid layers
    INTEGER,                   INTENT(IN)    :: nshift_total            ! Shift of vertical grid index for vertical nesting
    TYPE(datetime),            INTENT(IN)    :: tc_exp_startdate        ! Experiment start date
    TYPE(datetime),            INTENT(IN)    :: tc_exp_stopdate         ! Experiment end date
    REAL(wp),                  INTENT(IN)    :: start_time              ! Time at which execution of domain starts
    REAL(wp),                  INTENT(IN)    :: end_time                ! Time at which execution of domain ends
    REAL(wp),                  INTENT(IN)    :: dtime                   ! Fast-physics/advective time step for primary domain
    REAL(wp),                  INTENT(IN)    :: dt_fastphy              ! Fast physics' time step for this domain
    REAL(wp),                  INTENT(IN)    :: dt_rad_nwp              ! Tendency update period for radiation under NWP forcing
    REAL(wp),                  INTENT(INOUT) :: dt_grp_prevdom(:)       ! (ngrp) Tendency update period from previous domain
    CHARACTER(LEN=*),          INTENT(IN)    :: model_base_dir          ! Path for input files
    INTEGER,                   INTENT(IN)    :: msg_level               ! Message level
    INTEGER,                   INTENT(IN)    :: timers_level            ! Control parameter for timer
    TYPE(t_upatmo_phy_config), INTENT(INOUT) :: upatmo_phy_config       ! Upper-atmosphere configuration 
                                                                        ! with namelist settings
    TYPE(t_upatmo_aes_phy),    INTENT(INOUT) :: upatmo_aes_phy_config   ! Upper-atmosphere configuration for AES
    TYPE(t_upatmo_nwp_phy),    INTENT(INOUT) :: upatmo_nwp_phy_config   ! Upper-atmosphere configuration for NWP
    REAL(wp),        OPTIONAL, INTENT(IN)    :: vct_a(:)                ! (nlev+1) Nominal heights of grid layer interfaces

    ! Local variables
    TYPE (t_keyword_list), POINTER  :: keywords
    TYPE(timedelta), POINTER :: domTime
    TYPE(datetime) :: domStartDate, domEndDate
    REAL(wp) :: start_time_sggstn, dtime_sggstn
    INTEGER  :: jgrp, jgas, jprc, jtnd, jext
    INTEGER  :: imode, igrp, iext, cjgadj
    INTEGER  :: istat
    LOGICAL  :: l_belowtop, l_on, l_offline, l_exist, l_enabled, l_first
    CHARACTER(LEN=MAX_CHAR_LENGTH)       :: vname_prefix, filename
    CHARACTER(LEN=11) :: cjg
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: domTimeString
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':configure_upatmo_physics'

    !---------------------------------------------------------

    !-----------------------------------------------------
    !                 Consistency check
    !-----------------------------------------------------

    ! Actually, vct_a is only optional, because it is already optional one level higher
    IF (.NOT. PRESENT(vct_a)) CALL finish(routine, 'vct_a has to be present.')

    !-----------------------------------------------------
    !         Initialization with default values
    !-----------------------------------------------------

    CALL init_upatmo_phy_config( upatmo_aes_phy_config   = upatmo_aes_phy_config,   & !inout
      &                          upatmo_nwp_phy_config   = upatmo_nwp_phy_config    ) !inout

    !-----------------------------------------------------
    !                   Configuration
    !-----------------------------------------------------

    WRITE (cjg, '(i0)') jg
    cjgadj = VERIFY(cjg, " ")

    keywords => NULL()

    ! The configuration depends on what kind of forcing is used

    SELECT CASE(iforcing)
      
    CASE(iaes)
      
      !
      !*******************************************************************************
      !                                AES forcing
      !
      !*******************************************************************************

      ! Extreme-ultraviolet heating requires psrad orbit parameters.
      ! Only in case of AES-forcing, we overwrite the default settings 
      ! in 'src/namelists/mo_upatmo_nml' with potential namelist input to 'aes_rad_nml'.
      ! (In all other cases the input variables: 
      ! 'l_orbvsop87', 'cecc', 'cobld', 'clonp', 'lyr_perp' and 'yr_perp'
      ! might contain compiler-dependent default values and should not be used.) 
      upatmo_phy_config%orbit_type = MERGE(iorbit%vsop87, iorbit%kepler, l_orbvsop87)
      upatmo_phy_config%cecc       = cecc
      upatmo_phy_config%cobld      = cobld
      upatmo_phy_config%clonp      = clonp
      upatmo_phy_config%lyr_perp   = lyr_perp
      upatmo_phy_config%yr_perp    = yr_perp

      ! Process-wise configuration
      DO jprc = 1, iUpatmoPrcId%nitem

        ! Currently, the end height, above which the processes 
        ! would compute no more tendencies, is always the domain top
        upatmo_aes_phy_config%end_height( jprc ) = vct_a( 1 + nshift_total )
        upatmo_aes_phy_config%istartlev( jprc )  = 1

        ! Start height, above which processes start to compute tendencies.
        ! (Do not use any namelist input for the start height in case of NLTE.)
        CALL configure_start_height( start_height_nml = upatmo_phy_config%aes_start_height( jprc ),   & !in
          &                          nlev             = nlev,                                         & !in
          &                          nshift_total     = nshift_total,                                 & !in
          &                          start_height     = upatmo_aes_phy_config%start_height( jprc ),   & !inout
          &                          iendlev          = upatmo_aes_phy_config%iendlev( jprc ),        & !inout
          &                          vct_a            = vct_a,                                        & !(opt)in
          &                          opt_ldiscardnml  = jprc == iUpatmoPrcId%nlte                     ) !optin

      ENDDO  !jprc

      ! Status changes:
      ! * Required?
      upatmo_aes_phy_config%l_status(iUpatmoStat%required) = .FALSE.
      ! * Message output desired?
      upatmo_aes_phy_config%l_status(iUpatmoStat%message)  = msg_level >= imsg_thr%high
      ! * Timer monitoring desired?
      upatmo_aes_phy_config%l_status(iUpatmoStat%timer)    = timers_level > itmr_thr%med

    CASE(inwp)
      
      !
      !*******************************************************************************
      !                                 NWP forcing
      !
      !*******************************************************************************
      
      ! The following settings are only required, 
      ! if upper-atmosphere physics are switched on in 'nwp_phy_nml'
      IF (lupatmo_phy) THEN

        ! Loop over groups in which physics processes are clustered
        DO jgrp = 1, iUpatmoGrpId%nitem

          ! A process group is enabled, if its mode is either 'on' 
          ! or  'offline', and if its start height is below the model top
          imode      = upatmo_phy_config%nwp_grp( jgrp )%imode
          l_on       = imode == iUpatmoPrcMode%on
          l_offline  = imode == iUpatmoPrcMode%offline
          ! Please note that 'l_belowtop' is relatively useless, 
          ! if there is no explicit namelist entry, because its default value is < 0
          l_belowtop = upatmo_phy_config%nwp_grp( jgrp )%start_height < &
            & 0.5_wp * ( vct_a( 1 + nshift_total ) + vct_a( 2 + nshift_total ) )
          
          upatmo_nwp_phy_config%grp( jgrp )%imode                            = imode
          upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled ) = l_belowtop .AND.      &
            &                                                                  ( l_on .OR. l_offline )
          upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%offline ) = l_offline

        ENDDO  !jgrp

        ! Transfer information to the single processes
        DO jprc = 1, iUpatmoPrcId%nitem

          ! Process belongs to which group?
          igrp = upatmo_nwp_phy_config%prc( jprc )%igrp

          ! Mode
          upatmo_nwp_phy_config%prc( jprc )%imode = upatmo_nwp_phy_config%grp( igrp )%imode

          ! For the currently available group modes a group, which is enabled 
          ! means that all sub-processes are enabled, too
          upatmo_nwp_phy_config%prc( jprc )%l_stat( iUpatmoPrcStat%enabled ) = &
            & upatmo_nwp_phy_config%grp( igrp )%l_stat( iUpatmoPrcStat%enabled )

          ! The same holds for the offline-mode
          upatmo_nwp_phy_config%prc( jprc )%l_stat( iUpatmoPrcStat%offline ) = &
            & upatmo_nwp_phy_config%grp( igrp )%l_stat( iUpatmoPrcStat%offline )

        ENDDO  !jprc

      ELSE

        ! Some of the namelist parameters required special defaults 
        ! for the handling of secondary domains etc. 
        ! We reset them to reasonable values for 'lupatmo_phy = .FALSE.'
        DO jgrp = 1, iUpatmoGrpId%nitem
          upatmo_phy_config%nwp_grp( jgrp )%imode = iUpatmoPrcMode%off
        ENDDO
        DO jgas = 1, iUpatmoGasId%nitem
          upatmo_phy_config%nwp_gas( jgas )%imode = iUpatmoGasMode%zero
        ENDDO

      ENDIF  !lupatmo_phy

      ! Any parameterization switched on at all?
      upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled ) = & 
        & ANY( upatmo_nwp_phy_config%grp( : )%l_stat( iUpatmoPrcStat%enabled ) )

      ! Any parameterization in offline-mode?
      upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%offline ) = & 
        & ANY( upatmo_nwp_phy_config%grp( : )%l_stat( iUpatmoPrcStat%offline ) )

      ! Determine which variables will experience tendencies 
      ! from physics groups
      DO jprc = 1, iUpatmoPrcId%nitem

        ! Group, process is classified with
        igrp = upatmo_nwp_phy_config%prc( jprc )%igrp

        ! The information in 'l_update' on the process level
        ! is on whether a variable might experience tendencies 
        ! from the process in general.
        ! On the group level, 'l_update' is required for the handling
        ! of the accumulated tendencies within the NWP interface. 
        ! This means that we have to incorporate the information 
        ! whether the group is enabled and whether it is not in the offline-mode.
        ! In the latter case, the accumulated tendencies are not required.
        IF ( upatmo_nwp_phy_config%grp( igrp )%l_stat( iUpatmoPrcStat%enabled ) .AND. &       
          & .NOT. upatmo_nwp_phy_config%grp( igrp )%l_stat( iUpatmoPrcStat%offline )  ) THEN

          ! Loop over tendencies
          DO jtnd = 1, iUpatmoTendId%nitem

            upatmo_nwp_phy_config%grp( igrp )%l_update( jtnd ) =        &
              & upatmo_nwp_phy_config%grp( igrp )%l_update( jtnd ) .OR. &
              & upatmo_nwp_phy_config%prc( jprc )%l_update( jtnd )

          ENDDO  !jtnd

        ENDIF  !Group enabled and not in offline-mode?

      ENDDO  !jprc    

      ! Now, we transfer this information one level higher, 
      ! to see for which variables we have to care about an update at all.
      ! (Here, we spare us queries on 'enabled' and 'offline', 
      ! since defaulting all 'l_update' with .FALSE. in 'init_upatmo_phy_config' below 
      ! should lead to the correct result in any case.)
      DO jgrp = 1, iUpatmoGrpId%nitem
        
        DO jtnd = 1, iUpatmoTendId%nitem
          
          upatmo_nwp_phy_config%l_any_update( jtnd ) =         & 
            & upatmo_nwp_phy_config%l_any_update( jtnd ) .OR.  &
            & upatmo_nwp_phy_config%grp( jgrp )%l_update( jtnd )
          
        ENDDO  !jtnd
        
      ENDDO  !jgrp      

      ! Currently, the end height, above which the processes 
      ! would compute no more tendencies, is always the domain top
      DO jprc = 1, iUpatmoPrcId%nitem
        IF (upatmo_nwp_phy_config%prc( jprc )%l_stat( iUpatmoPrcStat%enabled )) THEN          
          upatmo_nwp_phy_config%prc( jprc )%end_height = vct_a( 1 + nshift_total )
          upatmo_nwp_phy_config%prc( jprc )%istartlev  = 1
        ENDIF
      ENDDO
      DO jgrp = 1, iUpatmoGrpId%nitem        
        IF (upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled )) THEN
          upatmo_nwp_phy_config%grp( jgrp )%end_height = vct_a( 1 + nshift_total )
          upatmo_nwp_phy_config%grp( jgrp )%istartlev  = 1
          ! Initialize the start height of a group
          upatmo_nwp_phy_config%grp( jgrp )%start_height = vct_a( 1 + nshift_total )
        ENDIF
      ENDDO

      ! Start heights, above which the processes compute tendencies
      DO jprc = 1, iUpatmoPrcId%nitem

        IF (upatmo_nwp_phy_config%prc( jprc )%l_stat( iUpatmoPrcStat%enabled )) THEN          

          ! Process belongs to which group?
          igrp = upatmo_nwp_phy_config%prc( jprc )%igrp

          ! Please note that a switch-off of single processes of a group, 
          ! which has been switched on, is currently not envisaged. 
          ! This means that single processes may be switched off implicitly, 
          ! if their default start height is above the domain top 
          ! and there is no or no valid entry for the start height 
          ! of the entire group!

          ! (Do not use any namelist input for the start height in case of NLTE.)
          CALL configure_start_height( start_height_nml = upatmo_phy_config%nwp_grp( igrp )%start_height, & !in
            &                          nlev             = nlev,                                           & !in
            &                          nshift_total     = nshift_total,                                   & !in
            &                          start_height     = upatmo_nwp_phy_config%prc( jprc )%start_height, & !inout
            &                          iendlev          = upatmo_nwp_phy_config%prc( jprc )%iendlev,      & !inout
            &                          vct_a            = vct_a,                                          & !(opt)in
            &                          opt_ldiscardnml  = jprc == iUpatmoPrcId%nlte                       ) !optin

          ! The start height for a group is the lowest start height among its processes
          ! (see its initialization above)
          IF (upatmo_nwp_phy_config%prc( jprc )%start_height < upatmo_nwp_phy_config%grp( igrp )%start_height) THEN
            upatmo_nwp_phy_config%grp( igrp )%start_height = upatmo_nwp_phy_config%prc( jprc )%start_height
            upatmo_nwp_phy_config%grp( igrp )%iendlev      = upatmo_nwp_phy_config%prc( jprc )%iendlev
          ENDIF

        ENDIF  !Enabled?

      ENDDO  !jprc

      ! Processing of the  radiatively active gases
      iext = iUpatmoExtdatId%gases

      ! Determine from single processes, which gases are required by a group.
      ! (This is comparable to the treatment of 'l_update' above.)
      DO jprc = 1, iUpatmoPrcId%nitem

        ! Group, process is classified with
        igrp = upatmo_nwp_phy_config%prc( jprc )%igrp

        IF (upatmo_nwp_phy_config%grp( igrp )%l_stat( iUpatmoPrcStat%enabled )) THEN

          ! Loop over gases
          DO jgas = 1, iUpatmoGasId%nitem

            ! Group requires which gas?
            upatmo_nwp_phy_config%grp( igrp )%l_gas( jgas ) =        &
              & upatmo_nwp_phy_config%grp( igrp )%l_gas( jgas ) .OR. &
              & upatmo_nwp_phy_config%prc( jprc )%l_gas( jgas ) 

            ! Gas is required by which group?
            upatmo_nwp_phy_config%gas( jgas )%l_grp( igrp ) =        &
              & upatmo_nwp_phy_config%gas( jgas )%l_grp( igrp ) .OR. &
              & upatmo_nwp_phy_config%gas( jgas )%l_prc( jprc )

          ENDDO  !jgas

        ENDIF  !Group enabled?

      ENDDO  !jprc

      ! Gas-wise evaluation
!NEC$ novector
      DO jgas = 1, iUpatmoGasId%nitem
        
        ! Is gas required by some of the process groups at all?
        DO jgrp = 1, iUpatmoGrpId%nitem

          upatmo_nwp_phy_config%gas( jgas )%l_stat( iUpatmoGasStat%enabled ) =        &
            & upatmo_nwp_phy_config%gas( jgas )%l_stat( iUpatmoGasStat%enabled ) .OR. &
            & upatmo_nwp_phy_config%gas( jgas )%l_grp( jgrp ) 

        ENDDO  !jgrp

        IF (upatmo_nwp_phy_config%gas( jgas )%l_stat( iUpatmoGasStat%enabled )) THEN

          imode = upatmo_phy_config%nwp_gas( jgas )%imode
          upatmo_nwp_phy_config%gas( jgas )%imode = imode

          ! External data for radiatively active gases required?
          upatmo_nwp_phy_config%extdat( iext )%l_stat( iUpatmoExtdatStat%required ) =        &
            & upatmo_nwp_phy_config%extdat( iext )%l_stat( iUpatmoExtdatStat%required ) .OR. &
            & (imode == iUpatmoGasMode%extdat)

        ENDIF
        
      ENDDO  !jgas

      ! Any gas required at all?
      upatmo_nwp_phy_config%l_gas_stat( iUpatmoGasStat%enabled ) = &
        & ANY(upatmo_nwp_phy_config%gas( : )%l_stat( iUpatmoGasStat%enabled ))

      ! External data for chemical heating tendencies required?
      iext = iUpatmoExtdatId%chemheat

      upatmo_nwp_phy_config%extdat( iext )%l_stat( iUpatmoExtdatStat%required ) = &
        & upatmo_nwp_phy_config%prc( iUpatmoPrcId%chemheat )%l_stat( iUpatmoPrcStat%enabled )

      ! Any external data required at all?
      upatmo_nwp_phy_config%l_extdat_stat( iUpatmoExtdatStat%required ) = &
        & ANY(upatmo_nwp_phy_config%extdat( : )%l_stat( iUpatmoExtdatStat%required ))

      ! Finalize name of external data files
      DO jext = 1, iUpatmoExtdatId%nitem

        CALL associate_keyword("<path>", TRIM(model_base_dir), keywords)
        filename = with_keywords(keywords, upatmo_phy_config%nwp_extdat( jext )%filename)
        upatmo_phy_config%nwp_extdat( jext )%filename = filename
        keywords => NULL()
   
      ENDDO  !jext

      ! Check, if files with external data exist
      IF (my_process_is_stdio()) THEN

        ! Loop over external data types
        DO jext = 1, iUpatmoExtdatId%nitem

          IF (upatmo_nwp_phy_config%extdat( jext )%l_stat( iUpatmoExtdatStat%required )) THEN

            INQUIRE(file = upatmo_phy_config%nwp_extdat( jext )%filename, exist=l_exist)
            
            IF (.NOT. l_exist) THEN
              message_text = 'The external data file: '                          &
                & //TRIM(upatmo_phy_config%nwp_extdat( jext )%filename)//' for ' &
                & //TRIM(upatmo_nwp_phy_config%extdat( jext )%longname)          &
                & //' is required, but it cannot be found.'
              CALL finish(routine, message_text)
            ENDIF

          ENDIF  !Required?

        ENDDO  !jext

      ENDIF  !IF (my_process_is_stdio())

      ! Check the varlists of 'output_nml':
      ! Some upper-atmosphere fields are only allocated, 
      ! if the corresponding process is switched on, 
      ! and only in this case is an output possible. 
      ! The name of any upatmo field that can be output, 
      ! has the prefix 'vname_prefix' (i.e. "upatmo_<name>", see initialization below), 
      ! so we search for this prefix.
      ! (The multiple invocation of 'is_variable_in_output_cond' is likely expensive 
      ! and definitely inefficient. However, we assume this to be bearable, 
      ! since it is only done once, during model setup.  
      ! In addition, the construction of a sophisticated inquiry tool 
      ! on output namelist settings would be an overkill for our purposes.)
      vname_prefix = upatmo_nwp_phy_config%vname_prefix
      IF ( .NOT. upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled ) .AND. &
        &  (LEN_TRIM(vname_prefix) > 0)                                     .AND. &
        &  is_variable_in_output_cond( var_name=vname_prefix,                     &
        &                              opt_dom=(/jg/) )                           ) THEN

        ! Provided the prefix for upper-atmosphere physics variable names has been chosen 
        ! distinguishable enough, it seems that the output of such a variable is requested, 
        ! which is not allowed under these circumstances
        message_text = 'Current namelist settings do not provide output of upatmo variables '    &
          & //'on domain '//cjg(cjgadj:)//', check output_nml-varlists for variables with prefix: ' &
          & //vname_prefix
        CALL finish(routine, message_text)

      ELSEIF ( .NOT. upatmo_nwp_phy_config%l_gas_stat( iUpatmoGasStat%enabled ) .AND. &
        &  is_variable_in_output_cond( var_name="group:upatmo_rad_gases",         &
        &                              opt_dom=(/jg/) )                           ) THEN

        ! Output of all radiatively active gases (varlist-group-prefix: 'group:') is only possible, 
        ! if radiation is switched on
        message_text = 'Current namelist settings do not provide output of upatmo gases ' &
          & //'on domain '//cjg(cjgadj:)//', check output_nml-varlists for: group:upatmo_rad_gases'
        CALL finish(routine, message_text)

      ELSEIF ( .NOT. upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled ) .AND. &
        &  is_variable_in_output_cond( var_name="group:upatmo_tendencies",        &
        &                              opt_dom=(/jg/) )                           ) THEN

        ! Output of all physics tendencies is only possible, if they are switched on
        message_text = 'Current namelist settings do not provide output of upatmo tendencies ' &
          & //'on domain '//cjg(cjgadj:)//', check output_nml-varlists for: group:upatmo_tendencies'
        CALL finish(routine, message_text)

      ELSEIF ( is_variable_in_output_cond( var_name=vname_prefix,                          &
        &                                  opt_dom=(/jg/),                                 &
        &                                  opt_filetype =(/FILETYPE_GRB, FILETYPE_GRB2/) ) ) THEN

        ! For the time being, we disable an output in the GRIB format 
        ! for the following reason:
        ! * Currently, we cannot guarantee that upper-atmosphere output variables 
        !   might overlap with output variables from standard ICON to the effect that 
        !   they would be indistinguishable with regard to their GRIB metadata 
        !   (at least without the specification of additional "non-standard" GRIB keys)
        message_text = 'The output of upatmo variables ' &
          & //'(desired for domain '//cjg(cjgadj:)//') in the GRIB format is not possible'
        CALL finish(routine, message_text)

      ELSEIF ( is_variable_in_output_cond( var_name="group:upatmo_rad_gases",              &
        &                                  opt_dom=(/jg/),                                 &
        &                                  opt_filetype =(/FILETYPE_GRB, FILETYPE_GRB2/) ) ) THEN

        message_text = 'The output of upatmo gases ' &
          & //'(desired for domain '//cjg(cjgadj:)//') in the GRIB format is not possible'
        CALL finish(routine, message_text)
  
      ELSEIF ( is_variable_in_output_cond( var_name="group:upatmo_tendencies",             &
        &                                  opt_dom=(/jg/),                                 &
        &                                  opt_filetype =(/FILETYPE_GRB, FILETYPE_GRB2/) ) ) THEN

        message_text = 'The output of upatmo tendencies ' &
          & //'(desired for domain '//cjg(cjgadj:)//') in the GRIB format is not possible'
        CALL finish(routine, message_text)
        
      ENDIF

      ! Set up flow control of events.

      ! Some notes:
      !
      ! * All upper-atmosphere physics processes are assumed 
      !   to be "slow" processes and so is 
      !   the handling of external data.
      !
      ! * Most of the settings below, are copied from 
      !   'src/configure_model/mo_atm_phy_nwp_config' 
      !   and 'src/configure_model/mo_aes_phy_config',
      !   so please see there for any descriptions. 
      !
      ! * Due to the extreme complexity of the NWP-events 
      !   the following modifications are applied:         
      !   - No IAU-time-shift, since currently the upper-atmosphere physics 
      !     are not intended for the extremely involved IAU-mode
      !     (for that reason alone, that no data for assimilation 
      !     are available for the upper atmosphere (although one could  
      !     let the latter evolve 'freely'))
      !   - The fast physics' time step 'dt_fastphy' is not added 
      !     to the event start time in the hope that 
      !     the physics tendencies would be computed as soon as 
      !     the model time exceeds the event start time
      !   - It might be desirable to borrow the AES-functionality of 
      !     starting physics processes not right away at model start, 
      !     but possibly later on (e.g., if a simulation is initialized 
      !     with IFS-data, the model atmosphere needs some spin-up time 
      !     to "calm down", and it might be advantageous to switch on 
      !     a process leading to strong tendencies only after this 
      !     turbulent spin-up phase to reduce the risk of model crashes).
      !   - Because of the previous point we make no use of an initial call, 
      !     how it is done for the "standard" NWP-physics
      !   - For all above mentioned points we have to be very careful 
      !     that the upper-atmosphere physics time stepping does not get 
      !     out of phase with 'dt_fastphy', otherwise newly computed 
      !     wind tendencies would not be accumulated!
      !
      ! * The event management infrastructure can be found in
      !   'src/atm_phy_nwp/mo_phy_events'
      !
      ! * The "heartbeat" period of ICON 
      !   is the fast physics' time step 'dt_fastphy'. 
      !   For this reason all update periods 'dt' have to be 
      !   multiples of 'dt_fastphy'. The multiple closest to 'dt'
      !   is computed according to:
      !
      !   dt_new = Max(1, Rounding(dt / dt_fastphy)) * dt_fasphy.
      ! 
      !   'dt_new' may be greater or lower than 'dt', 
      !   which could have undesirable effects in multi-domain simulations, 
      !   as demonstrated by the following example:
      !
      !   A simulation comprises the global domain plus one nest, 
      !   with the following namelist settings: 
      !   - dtime = 180 s
      !   - dt    = 600 s for all domains
      !   This means for domain 1:
      !   - dt_fastphy = dtime = 180 s
      !   - dt                 = 600 s
      !   - dt_new             = 540 s
      !   and for domain 2:
      !   - dt_fastphy = dtime / 2 =  90 s
      !   - dt                     = 600 s
      !   - dt_new                 = 630 s
      !
      !   That 'dt_new(nest) > dt_new(global domain)' may happen 
      !   is probably undesirable. In order to avoid this, 
      !   we set the result for 'dt_new' on the previous domain 
      !   (stored in 'dt_grp_prevdom') as upper limit for 'dt_new' 
      !   on the current domain (so 'dt_new = 630 s' would be reset 
      !   to 'dt_new = 540 s' in the above example).
      !   This applies only to the update periods of the physics tendencies, 
      !   since they can be set domain-wise in the namelist. 
      !   For the update periods of the time interpolation of the external data 
      !   only one value can be set in the namelist, which applies to all domains.
      !
      ! * To consider: the event management subroutines in 'mo_aes_phy_config', 
      !   'mo_atm_phy_nwp_config' and this one might be merged, 
      !   in order to reduce the code overhead and to facilitate maintenance.

      IF (upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled )) THEN

        domTime => NULL()

        ! Set up start date of domain
        IF (jg > 1) THEN
          IF (MOD(start_time, dtime) > eps) THEN
            start_time_sggstn = MAX(1._wp, ANINT(start_time / dtime)) * dtime
            dtime_sggstn      = start_time / MAX(1._wp, ANINT(start_time / dtime))
            message_text      = "start_time(dom"//cjg(cjgadj:)//") = "                &
              & //TRIM(ADJUSTL(real2string(start_time, opt_fmt="(F20.3)")))        &
              & //" s is no multiple of dtime = "                                  &
              & //TRIM(ADJUSTL(real2string(dtime, opt_fmt="(F20.3)")))             & 
              & //" s! Please, either adjust start_time (choose, e.g., "           &
              & //TRIM(ADJUSTL(real2string(start_time_sggstn, opt_fmt="(F20.3)"))) &
              & //" s) or adjust dtime (choose, e.g., "                            &
              & //TRIM(ADJUSTL(real2string(dtime_sggstn, opt_fmt="(F20.3)")))//" s)."
            CALL finish(routine, message_text)
          ENDIF
          CALL getPTStringFromMS(INT(start_time * 1000._wp, i8), domTimeString)
          domTime => newTimedelta(domTimeString, istat)
          IF (istat /= SUCCESS) CALL finish(routine, "Invalid domain start time string.")
          domStartDate = tc_exp_startdate + domTime
          CALL deallocateTimedelta(domTime)
        ELSE
          domStartDate = tc_exp_startdate
        ENDIF

        ! Set up end date of domain
        IF (jg > 1) THEN
          IF (end_time /= DEFAULT_ENDTIME) THEN
            CALL getPTStringFromMS(INT(end_time * 1000._wp, i8), domTimeString)
            domTime => newTimedelta(domTimeString, istat)
            IF (istat /= SUCCESS) CALL finish(routine, "Invalid domain end time string.")
            domEndDate = tc_exp_startdate + domTime
            CALL deallocateTimedelta(domTime)
            IF (domEndDate > tc_exp_stopdate) domEndDate = tc_exp_stopdate
          ELSE
            domEndDate = tc_exp_stopdate
          ENDIF
        ELSE
          domEndDate = tc_exp_stopdate
        ENDIF

        ! The tendency update period for upper-atmosphere radiation 
        ! has to divide the update period for NWP radiation evenly. 
        ! (Please see 'src/upper_atmosphere/mo_nwp_upatmo_interface' for an explanation.)
        IF ( upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%enabled ) .AND. &
          &  (MOD(dt_rad_nwp, upatmo_phy_config%nwp_grp( iUpatmoGrpId%rad )%dt) > eps)            ) THEN
          message_text = "WARNING, update period dt for group "                         &
            & //TRIM(upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%name)//" on domain " &
            & //cjg(cjgadj:)//" is adjusted to evenly divide dt_rad for NWP forcing..."
          CALL message(routine, message_text)
          message_text = "... its old value was " &
            & //TRIM(ADJUSTL(real2string(upatmo_phy_config%nwp_grp( iUpatmoGrpId%rad )%dt, opt_fmt="(F20.3)")))//" s"
          CALL message(routine, message_text)
          upatmo_phy_config%nwp_grp( iUpatmoGrpId%rad )%dt = &
            & dt_rad_nwp / MAX(1._wp, ANINT(dt_rad_nwp / upatmo_phy_config%nwp_grp( iUpatmoGrpId%rad )%dt))
          message_text = "... its new value is " &
            & //TRIM(ADJUSTL(real2string(upatmo_phy_config%nwp_grp( iUpatmoGrpId%rad )%dt, opt_fmt="(F20.3)")))//" s"
          CALL message(routine, message_text)
        ENDIF

        ! Construct event group
        CALL upatmo_nwp_phy_config%event_mgmt_grp%construct( grpName = 'phyNwpUpatmoEventGroupGrp', & 
          &                                                  pid     = jg,                          & 
          &                                                  grpSize = iUpatmoGrpId%nitem           )

        l_first = .TRUE.
        upatmo_nwp_phy_config%lBeforeOpPhase = .TRUE.
        upatmo_nwp_phy_config%lAfterOpPhase  = .FALSE.

        DO jgrp = 1, iUpatmoGrpId%nitem

          l_enabled = upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled )

          ! Set up process group events 
          ! and integrate them into event group
          CALL configure_nwp_event( eventStartDateIn       = upatmo_phy_config%nwp_grp( jgrp )%t_start, & !in
            &                       eventEndDateIn         = upatmo_phy_config%nwp_grp( jgrp )%t_end,   & !in
            &                       eventIntervalIn        = upatmo_phy_config%nwp_grp( jgrp )%dt,      & !in
            &                       domainStartDate        = domStartDate,                              & !in
            &                       domainEndDate          = domEndDate,                                & !in
            &                       basicInterval          = dt_fastphy,                                & !in
            &                       eventName              = upatmo_nwp_phy_config%grp( jgrp )%name,    & !in
            &                       domainName             = cjg(cjgadj:),                              & !in
            &                       eventId                = upatmo_nwp_phy_config%grp( jgrp )%id,      & !in
            &                       eventEnabled           = l_enabled,                                 & !in
            &                       eventObject            = upatmo_nwp_phy_config%event_mgmt_grp,      & !inout
            &                       eventSubobject         = upatmo_nwp_phy_config%grp( jgrp )%event,   & !inout
            &                       optEventReqInit        = .FALSE.,                                   & !optin
            &                       optEventInclStart      = .TRUE.,                                    & !optin
            &                       optLowerBound4Interval = dt_fastphy,                                & !optin
            &                       optUpperBound4Interval = dt_grp_prevdom( jgrp ),                    & !optin
            &                       optEventStartDateOut   = upatmo_nwp_phy_config%grp( jgrp )%t_start, & !optout
            &                       optEventEndDateOut     = upatmo_nwp_phy_config%grp( jgrp )%t_end,   & !optout
            &                       optEventIntervalOut    = upatmo_nwp_phy_config%grp( jgrp )%dt       ) !optout
          
          ! Update 'dt_grp_prevdom' for the next call of this subroutine
          ! on the subsequent domain (if there is one)
          dt_grp_prevdom( jgrp ) = upatmo_nwp_phy_config%grp( jgrp )%dt

          upatmo_nwp_phy_config%grp( jgrp )%lBeforeOpPhase = .TRUE.
          upatmo_nwp_phy_config%grp( jgrp )%lAfterOpPhase  = .FALSE.

          IF (l_enabled .AND. l_first) THEN
            upatmo_nwp_phy_config%t_start = upatmo_nwp_phy_config%grp( jgrp )%t_start
            upatmo_nwp_phy_config%t_end   = upatmo_nwp_phy_config%grp( jgrp )%t_end
            l_first = .FALSE.
          ELSEIF (l_enabled) THEN
            IF (upatmo_nwp_phy_config%t_start > upatmo_nwp_phy_config%grp( jgrp )%t_start) &
              & upatmo_nwp_phy_config%t_start = upatmo_nwp_phy_config%grp( jgrp )%t_start
            IF (upatmo_nwp_phy_config%grp( jgrp )%t_end > upatmo_nwp_phy_config%t_end) &
              & upatmo_nwp_phy_config%t_end = upatmo_nwp_phy_config%grp( jgrp )%t_end
          ENDIF

        ENDDO  !jgrp

        ! Repeat (almost) the same for the update periods 
        ! for the time interpolation of the external data.
        
        IF (upatmo_nwp_phy_config%l_extdat_stat( iUpatmoExtdatStat%required )) THEN

          CALL upatmo_nwp_phy_config%event_mgmt_extdat%construct( grpName = 'phyNwpUpatmoEventGroupExtdat', & 
            &                                                     pid     = jg,                             & 
            &                                                     grpSize = iUpatmoExtdatId%nitem           )

          DO jext = 1, iUpatmoExtdatId%nitem

            l_enabled = upatmo_nwp_phy_config%extdat( jext )%l_stat( iUpatmoExtdatStat%required )

            CALL configure_nwp_event( eventStartDateIn       = ' ',                                       & !in
              &                       eventEndDateIn         = ' ',                                       & !in
              &                       eventIntervalIn        = upatmo_phy_config%nwp_extdat( jext )%dt,   & !in
              &                       domainStartDate        = domStartDate,                              & !in
              &                       domainEndDate          = domEndDate,                                & !in
              &                       basicInterval          = dt_fastphy,                                & !in
              &                       eventName              = upatmo_nwp_phy_config%extdat( jext )%name, & !in
              &                       domainName             = cjg(cjgadj:),                              & !in
              &                       eventId                = upatmo_nwp_phy_config%extdat( jext )%id,   & !in
              &                       eventEnabled           = l_enabled,                                 & !in
              &                       eventObject            = upatmo_nwp_phy_config%event_mgmt_extdat,   & !inout
              &                       eventSubobject         = upatmo_nwp_phy_config%extdat( jext )%event,& !inout
              &                       optEventReqInit        = .FALSE.,                                   & !optin
              &                       optEventInclStart      = .TRUE.,                                    & !optin
              &                       optLowerBound4Interval = dt_fastphy,                                & !optin
              &                       optEventIntervalOut    = upatmo_nwp_phy_config%extdat( jext )%dt    ) !optout

          ENDDO  !jext

        ENDIF  !External data required?

        ! Thermodynamic coupling of physics:
        ! Notation:
        ! * T, p, rho -> Temperature, pressure, density
        ! * "D"       -> Parcel-fixed, Lagrangian system ("convenience" system of theoretical considerations)
        ! * "d"       -> Cell-fixed, Eulerian system (the actual system at hand in ICON)
        ! * Q         -> Heat sources and heat flux divergences from parameterizations multiplied by temperature
        ! * e         -> Mass-specific energy
        ! * s         -> Mass-specific heat (entropy)
        ! * mu        -> Chemical potential
        ! * theta     -> Potential temperature
        ! * Exner     -> Exner pressure
        ! * p00       -> Const. reference pressure

        ! Please note: a single coupling factor works only 
        ! as long as cp, cv and R are const.
        
        SELECT CASE(upatmo_phy_config%nwp_thermdyn_cpl)
        CASE (iThermdynCoupling%isobaric)
          ! 1st Isobaric coupling:
          !
          ! Gibbs' eq. expanded according to:
          ! rho * cp * DT/Dt - Dp/Dt = Q 
          !
          ! The state change due to the parameterized process 
          ! is assumed to be isobaric:
          ! 
          ! Dp/Dt = 0  =>  DT/Dt = Q / rho / cp
          !
          ! Practically all parameterizations compute the temperature tendency accordingly, 
          ! since they were developed in the framework of hydrostatic models.
          ! Consequently:
          upatmo_nwp_phy_config%thermdyn_cpl_fac = 1._wp
        CASE (iThermdynCoupling%isochoric)
          ! 2nd Isochoric (isodensic) coupling:
          !
          ! Gibbs' eq. expanded according to:
          ! rho * cv * DT/Dt - p / rho * Drho/Dt = Q 
          !
          ! In ICON the temperature tendencies are computed 
          ! according to the isochoric approximation 
          ! (or isodensic from the cell-fixed, Eulerian point of view):
          !
          ! Drho/Dt = 0  =>  DT/Dt = Q / rho / cv
          !
          ! Since cv < cp, |DT/Dt|_isobaric < |DT/Dt|_isochoric. 
          ! This excess is required to go into the volume work, 
          ! when the tendencies DT/dt are dribbled into the dynamics (for the slow physics).
          ! As a consequence the complete state change (physics + dynamics) 
          ! is neither isochoric nor isobaric. 
          ! Since the parameterizations compute DT/Dt in the isobaric way, 
          ! the coupling factor has to compensate that:
          upatmo_nwp_phy_config%thermdyn_cpl_fac = cpd / cvd
        CASE (iThermdynCoupling%entropic)
          ! 3rd Entropic coupling:
          !
          ! For the cell-fixed systems of ICON 
          ! Gibbs' eq. (without momentum and gravity) reads:
          ! d(rho*e)/dt = T * d(rho*s)/dt + mu * drho/dt,
          !
          ! where drho/dt is assumed to be exclusively treated by the dynamics:
          !
          ! drho/dt = (drho/dt)_dyn + (drho/dt)_phy, with (drho/dt)_phy = 0.   (1)
          !
          ! d(rho*s)/dt is partly covered by the dynamics and partly covered by the physics:
          !
          ! d(rho*s)/dt = (d(rho*s)/dt)_dyn + (d(rho*s)/dt)_phy,   (2)
          !
          ! with (d(rho*s)/dt)_phy = ... + Q / T,
          !
          ! where Q represents sources of heat from irreversible subgrid-scale processes
          ! plus a divergence of irreversible heat fluxes.
          ! Whether the state change takes place constrained (e.g., isobarically or isochorically), 
          ! makes no difference at this place.
          ! In ICON eq. (2) is used implicitly via the potential temperature 
          ! or Exner pressure equations: 
          ! 
          ! (2) * theta / cp  =>  d(rho*theta)/dt = ... + Q / Exner / cp  | * R * Exner / rho / theta / cv  (3)
          !                   =>  dExner/dt = ... + R * Q / rho / theta / cp / cv                           (4)
          !
          ! Eq. (4), however, is not used directly within the physics. 
          ! Given that Exner = (p / p00)**(R/cp) and the equation of state p = rho * R * T, 
          ! the following expansion is used instead:
          !
          ! dExner/dt = (R/cp) * (p / p00)**(R/cp-1) * dp/dt / p00 
          !           = R * Exner / cp * (dT/dt / T + drho/dt / rho).   (5)
          !
          ! Eq. (5) is used to transform the temperature tendencies into Exner pressure tendencies. 
          ! In case of the isochoric coupling, the physics part of Eq. (5) reads:
          !
          ! (dExner/dt)_phy = R * Exner / cp * (dT/dt)_phy / T = R * Q / rho / theta / cp / cv, 
          !
          ! where eq. (1) is applied and (dT/dt)_phy = (DT/dt)_phy = Q / rho / cv is assumed.
          ! Under these circumstances eq. (4) becomes equal to eq. (3), 
          ! so that the entropic coupling factor is the same as with the isochoric coupling:
          upatmo_nwp_phy_config%thermdyn_cpl_fac = cpd / cvd
        CASE DEFAULT
          CALL finish(routine, "Invalid thermdyn_cpl")
        END SELECT
        
      ENDIF  !Physics enabled?
      
      ! Status changes:
      ! * Required?
      upatmo_nwp_phy_config%l_status(iUpatmoStat%required) = &
        & upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled )
      ! * Message output desired?
      upatmo_nwp_phy_config%l_status(iUpatmoStat%message)  = msg_level >= imsg_thr%high
      ! * Timer monitoring desired?
      upatmo_nwp_phy_config%l_status(iUpatmoStat%timer)    = timers_level > itmr_thr%med
      
    END SELECT  !SELECT CASE(iforcing)

    !
    !*******************************************************************************
    !                                 Miscellaneous
    !
    !*******************************************************************************
    
    ! Status changes:
    ! * Configured?
    upatmo_aes_phy_config%l_status(iUpatmoStat%configured)   = .TRUE.
    upatmo_nwp_phy_config%l_status(iUpatmoStat%configured)   = .TRUE.
    
  END SUBROUTINE configure_upatmo_physics

  !==================================================================================== 

  !>
  !! Initialize configuration types with default values.
  !!
  SUBROUTINE init_upatmo_phy_config( upatmo_aes_phy_config,   & !inout
    &                                upatmo_nwp_phy_config    ) !inout

    ! In/out variables
    TYPE(t_upatmo_aes_phy),   INTENT(INOUT) :: upatmo_aes_phy_config   ! Upper-atmosphere configuration for AES
    TYPE(t_upatmo_nwp_phy),   INTENT(INOUT) :: upatmo_nwp_phy_config   ! Upper-atmosphere configuration for NWP

    ! Local variables
    INTEGER  :: jprc, jgrp, jgas, jext

    !---------------------------------------------------------

    ! Please note:
    ! * Many of the string variables in the following will become part 
    !   of variable names, descriptions and other metadata that enter NetCDF (and GRIB) output. 
    !   Please, be careful with special characters (even stuff like '.', ':').  
    !   Not everything might be digestible by the two formats (an internal check 
    !   for format compliance would be too complicated currently). 
    ! * Be careful, if you change the defaults for the switches 'l_...', 
    !   the above configuration often builds on them implicitly.
    ! * Be careful, if you change the order. Some settings depend on 
    !   preceding settings.

    !
    !*******************************************************************************
    !                                AES forcing
    !
    !*******************************************************************************

    ! Initialize the status switches, 
    ! but skip 'iUpatmoStat%checked', so as not to overwrite 
    ! its assignment in 'src/namelists/mo_upatmo_nml: check_upatmo'
    CALL init_logical_1d( variable=upatmo_aes_phy_config%l_status,                          &
      &                   value=.FALSE., opt_ilist=(/iUpatmoStat%checked/), opt_mask="list" )

    ! Start height, above which processes compute tendencies
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%srbc )     = startHeightDef%srbc
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%nlte )     = startHeightDef%nlte
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%euv )      = startHeightDef%euv
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%vdfmol )   = startHeightDef%vdfmol
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%fric )     = startHeightDef%fric
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%iondrag )  = startHeightDef%iondrag
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%joule )    = &  ! Joule heating gets same start height as ion drag
      & upatmo_aes_phy_config%start_height( iUpatmoPrcId%iondrag )
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%no )       = startHeightDef%no
    upatmo_aes_phy_config%start_height( iUpatmoPrcId%chemheat ) = startHeightDef%chemheat
    ! End height, below which processes compute tendencies,  
    ! and grid layer indices
    DO jprc = 1, iUpatmoPrcId%nitem
      upatmo_aes_phy_config%end_height( jprc ) = -999._wp
      upatmo_aes_phy_config%istartlev( jprc )  = 0
      upatmo_aes_phy_config%iendlev( jprc )    = 0
    ENDDO  !jprc    
    
    !
    !*******************************************************************************
    !                                 NWP forcing
    !
    !*******************************************************************************

    ! Initialize the status switches, 
    ! but skip 'iUpatmoStat%checked', so as not to overwrite 
    ! its assignment in 'src/namelists/mo_upatmo_nml: check_upatmo'
    CALL init_logical_1d( variable=upatmo_nwp_phy_config%l_status,                          &
      &                   value=.FALSE., opt_ilist=(/iUpatmoStat%checked/), opt_mask="list" )

    ! Some overall switches
    upatmo_nwp_phy_config%l_phy_stat(:)          = .FALSE.
    upatmo_nwp_phy_config%l_gas_stat(:)          = .FALSE.
    upatmo_nwp_phy_config%l_extdat_stat(:)       = .FALSE.
    upatmo_nwp_phy_config%l_any_update(:)        = .FALSE.
    !
    upatmo_nwp_phy_config%vname_prefix = 'upatmo_'  ! Prefix for variable names in output_nml-varlists.
                                                    ! (Please, choose something 'unique' enough to be very unlikely 
                                                    ! used by other namelists, and please update 
                                                    ! Namelist_overview: upatmo_nml, if you would change it. Thank you!)

    ! Control units for single processes
    DO jprc = 1, iUpatmoPrcId%nitem
      upatmo_nwp_phy_config%prc( jprc )%id           = jprc
      upatmo_nwp_phy_config%prc( jprc )%imode        = iUpatmoPrcMode%off
      upatmo_nwp_phy_config%prc( jprc )%dt           = -999._wp
      upatmo_nwp_phy_config%prc( jprc )%start_height = -999._wp
      upatmo_nwp_phy_config%prc( jprc )%end_height   = -999._wp
      upatmo_nwp_phy_config%prc( jprc )%istartlev    = 0
      upatmo_nwp_phy_config%prc( jprc )%iendlev      = 0
      upatmo_nwp_phy_config%prc( jprc )%l_update(:)  = .FALSE.
      upatmo_nwp_phy_config%prc( jprc )%l_gas(:)     = .FALSE.
      upatmo_nwp_phy_config%prc( jprc )%l_stat(:)    = .FALSE.
    ENDDO  !jprc
    ! Individual settings
    !
    ! Heating due to Schumann-Runge bands and continuum of O2
    jprc = iUpatmoPrcId%srbc
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'srbc'
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'heating due to Schumann-Runge bands and continuum of O2'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%rad
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%o2 )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%srbc
    !
    ! Heating due to processes in non-local-thermodynamic-equilibrium
    jprc = iUpatmoPrcId%nlte
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'nlte'
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'heating due to processes in non-local-thermodynamic-equilibrium'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%rad
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%o3 )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%o2 )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%o )          = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%co2 )        = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%n2 )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%nlte
    !
    ! Extreme-ultraviolet heating
    jprc = iUpatmoPrcId%euv
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'euv' 
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'extreme-ultraviolet heating'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%rad
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%o2 )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%o )          = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%n2 )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%euv
    !
    ! Molecular diffusion
    jprc = iUpatmoPrcId%vdfmol
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'vdfmol' 
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'molecular diffusion'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%imf
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%wind_h ) = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%qx )     = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%vdfmol
    !
    ! Frictional heating
    jprc = iUpatmoPrcId%fric
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'fric' 
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'frictional heating'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%imf
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE. 
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%fric
    !
    ! Ion drag
    jprc = iUpatmoPrcId%iondrag
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'iondrag' 
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'ion drag'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%imf
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%wind_h ) = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%iondrag
    !
    ! Joule heating
    jprc = iUpatmoPrcId%joule
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'joule' 
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'Joule heating'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%imf
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = &  ! The same as ion drag
      & upatmo_nwp_phy_config%prc( iUpatmoPrcId%iondrag )%start_height
    !
    ! Near-infrared heating by NO
    jprc = iUpatmoPrcId%no
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'no' 
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'near-infrared heating by NO'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%rad
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%l_gas( iUpatmoGasId%no )         = .TRUE.
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%no
    !
    ! Chemical heating
    jprc = iUpatmoPrcId%chemheat
    upatmo_nwp_phy_config%prc( jprc )%name                             = 'chemheat'
    upatmo_nwp_phy_config%prc( jprc )%longname                         = &
      & 'chemical heating'
    upatmo_nwp_phy_config%prc( jprc )%igrp                             = iUpatmoGrpId%rad
    upatmo_nwp_phy_config%prc( jprc )%l_update( iUpatmoTendId%temp )   = .TRUE. 
    upatmo_nwp_phy_config%prc( jprc )%start_height                     = startHeightDef%chemheat

    ! Control units for groups in which processes are clustered
    DO jgrp = 1, iUpatmoGrpId%nitem
      upatmo_nwp_phy_config%grp( jgrp )%id           = jgrp
      upatmo_nwp_phy_config%grp( jgrp )%imode        = iUpatmoPrcMode%off
      upatmo_nwp_phy_config%grp( jgrp )%igrp         = -1                  ! does not apply
      upatmo_nwp_phy_config%grp( jgrp )%dt           = -999._wp
      upatmo_nwp_phy_config%grp( jgrp )%start_height = -999._wp
      upatmo_nwp_phy_config%grp( jgrp )%end_height   = -999._wp
      upatmo_nwp_phy_config%grp( jgrp )%istartlev    = 0
      upatmo_nwp_phy_config%grp( jgrp )%iendlev      = 0
      upatmo_nwp_phy_config%grp( jgrp )%l_update(:)  = .FALSE.
      upatmo_nwp_phy_config%grp( jgrp )%l_gas(:)     = .FALSE.
      upatmo_nwp_phy_config%grp( jgrp )%l_stat(:)    = .FALSE.
    ENDDO  !jgrp
    ! Individual settings
    !
    ! Ion drag (I), molecular diffusion (M) and frictional heating (F)
    upatmo_nwp_phy_config%grp( iUpatmoGrpId%imf )%name     = 'imf'
    upatmo_nwp_phy_config%grp( iUpatmoGrpId%imf )%longname = & 
      & 'ion drag, molecular diffusion and frictional heating'
    !
    ! Radiation and chemical heating 
    upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%name     = 'rad'
    upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%longname = &
      & 'radiation and chemical heating '

    ! Control units for gases
    DO jgas = 1, iUpatmoGasId%nitem
      upatmo_nwp_phy_config%gas( jgas )%id        = jgas
      upatmo_nwp_phy_config%gas( jgas )%unit      = 'kg kg-1'
      upatmo_nwp_phy_config%gas( jgas )%imode     = iUpatmoGasMode%zero
      upatmo_nwp_phy_config%gas( jgas )%l_prc(:)  = .FALSE.
      upatmo_nwp_phy_config%gas( jgas )%l_grp(:)  = .FALSE.
      upatmo_nwp_phy_config%gas( jgas )%l_stat(:) = .FALSE.
    ENDDO  !jgas
    ! Individual settings
    !
    ! Ozone
    jgas = iUpatmoGasId%o3
    upatmo_nwp_phy_config%gas( jgas )%name                       = 'o3'
    upatmo_nwp_phy_config%gas( jgas )%longname                   = 'ozone'
    upatmo_nwp_phy_config%gas( jgas )%vmr2mmr                    = amo3 / amd  ! Conversion factor
    upatmo_nwp_phy_config%gas( jgas )%mmr2vmr                    = amd / amo3  ! Conversion factor
    upatmo_nwp_phy_config%gas( jgas )%mass2mol                   = avo / amo3  ! Conversion factor
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%nlte ) = .TRUE.
    !
    ! Dioxygen
    jgas = iUpatmoGasId%o2
    upatmo_nwp_phy_config%gas( jgas )%name                       = 'o2'
    upatmo_nwp_phy_config%gas( jgas )%longname                   = 'dioxygen'
    upatmo_nwp_phy_config%gas( jgas )%vmr2mmr                    = amo2 / amd   
    upatmo_nwp_phy_config%gas( jgas )%mmr2vmr                    = amd / amo2
    upatmo_nwp_phy_config%gas( jgas )%mass2mol                   = avo / amo2
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%srbc ) = .TRUE.
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%nlte ) = .TRUE.
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%euv )  = .TRUE.
    !
    ! Atomic oxygen 
    jgas = iUpatmoGasId%o
    upatmo_nwp_phy_config%gas( jgas )%name                       = 'o'
    upatmo_nwp_phy_config%gas( jgas )%longname                   = 'atomic oxygen'
    upatmo_nwp_phy_config%gas( jgas )%vmr2mmr                    = amo / amd
    upatmo_nwp_phy_config%gas( jgas )%mmr2vmr                    = amd / amo
    upatmo_nwp_phy_config%gas( jgas )%mass2mol                   = avo / amo
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%nlte ) = .TRUE.
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%euv )  = .TRUE.
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%no )   = .TRUE.
    !
    ! Carbon dioxide 
    jgas = iUpatmoGasId%co2
    upatmo_nwp_phy_config%gas( jgas )%name                       = 'co2'
    upatmo_nwp_phy_config%gas( jgas )%longname                   = 'carbon dioxide'
    upatmo_nwp_phy_config%gas( jgas )%vmr2mmr                    = amco2 / amd
    upatmo_nwp_phy_config%gas( jgas )%mmr2vmr                    = amd / amco2
    upatmo_nwp_phy_config%gas( jgas )%mass2mol                   = avo / amco2
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%nlte ) = .TRUE.
    ! 
    ! Nitric oxide
    jgas = iUpatmoGasId%no
    upatmo_nwp_phy_config%gas( jgas )%name                     = 'no'
    upatmo_nwp_phy_config%gas( jgas )%longname                 = 'nitric oxide'
    upatmo_nwp_phy_config%gas( jgas )%vmr2mmr                  = amno / amd
    upatmo_nwp_phy_config%gas( jgas )%mmr2vmr                  = amd / amno
    upatmo_nwp_phy_config%gas( jgas )%mass2mol                 = avo / amno
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%no ) = .TRUE.
    !
    ! Dinitrogen => diagnostic mode!
    jgas = iUpatmoGasId%n2
    upatmo_nwp_phy_config%gas( jgas )%name                       = 'n2'
    upatmo_nwp_phy_config%gas( jgas )%longname                   = 'dinitrogen'
    upatmo_nwp_phy_config%gas( jgas )%vmr2mmr                    = amn2 / amd
    upatmo_nwp_phy_config%gas( jgas )%mmr2vmr                    = amd / amn2
    upatmo_nwp_phy_config%gas( jgas )%mass2mol                   = avo / amn2 
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%nlte ) = .TRUE.
    upatmo_nwp_phy_config%gas( jgas )%l_prc( iUpatmoPrcId%euv )  = .TRUE.

    ! Control units for external data
    DO jext = 1, iUpatmoExtdatId%nitem
      upatmo_nwp_phy_config%extdat( jext )%id        = jext
      upatmo_nwp_phy_config%extdat( jext )%dt        = -999._wp
      upatmo_nwp_phy_config%extdat( jext )%l_stat(:) = .FALSE.
    ENDDO  !jext
    ! Individual settings
    !
    ! Gas concentrations
    jext = iUpatmoExtdatId%gases
    upatmo_nwp_phy_config%extdat( jext )%name     = 'gases'
    upatmo_nwp_phy_config%extdat( jext )%longname = 'radiatively active gases'
    !
    ! Chemical heating tendencies
    jext = iUpatmoExtdatId%chemheat
    upatmo_nwp_phy_config%extdat( jext )%name     = 'chemheat'
    upatmo_nwp_phy_config%extdat( jext )%longname = 'chemical heating tendencies'

    ! Thermodynamic coupling factor
    upatmo_nwp_phy_config%thermdyn_cpl_fac = 1._wp
    
  END SUBROUTINE init_upatmo_phy_config

  !==================================================================================== 

  !>
  !! Configure start height.
  !!
  SUBROUTINE configure_start_height( start_height_nml,  & !in
    &                                nlev,              & !in
    &                                nshift_total,      & !in
    &                                start_height,      & !inout
    &                                iendlev,           & !inout
    &                                vct_a,             & !(opt)in
    &                                opt_ldiscardnml    ) !optin
    
    ! In/out variables
    REAL(wp),           INTENT(IN)    :: start_height_nml  ! Start height from namelist
    INTEGER,            INTENT(IN)    :: nlev              ! Number of grid layers
    INTEGER,            INTENT(IN)    :: nshift_total      ! Shift of vertical grid index for vertical nesting
    REAL(wp),           INTENT(INOUT) :: start_height      ! Configured start height
    INTEGER,            INTENT(INOUT) :: iendlev           ! Grid layer corresponding to start_height
    REAL(wp), OPTIONAL, INTENT(IN)    :: vct_a(:)          ! Nominal heights of grid layer interfaces
    LOGICAL,  OPTIONAL, INTENT(IN)    :: opt_ldiscardnml   ! .TRUE.: do not use start_height_nml, 
                                                           ! but start_height in any case
    
    ! Local variables
    INTEGER  :: jk, jks
    LOGICAL  :: ldiscardnml
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':configure_start_height'

    !---------------------------------------------------------

    IF (.NOT. PRESENT(vct_a)) THEN 
      CALL finish(routine, 'vct_a still uninitialized.')
    ELSEIF (start_height < 0._wp) THEN
      ! A valid in-value for start_height is required
      CALL finish(routine, 'start_height requires non-negative in-value.')
    ENDIF
    
    IF (PRESENT(opt_ldiscardnml)) THEN
      ldiscardnml = opt_ldiscardnml
    ELSE
      ldiscardnml = .FALSE.
    ENDIF

    ! A valid namelist entry overwrites the default start height (if opt_ldiscardnml = .FALSE.)
    IF (.NOT. (start_height_nml < 0._wp .OR. ldiscardnml)) start_height = start_height_nml

    ! Derive the grid layer index corresponding to start height
    IF (.NOT. ((start_height - 0.5_wp * (vct_a(1+nshift_total)+vct_a(2+nshift_total))) < 0._wp)) THEN
      ! Start height lies at or above the center of the uppermost grid layer:
      ! process implicitly switched off!
      iendlev = 0
    ELSEIF (.NOT. ((0.5_wp * (vct_a(nlev+nshift_total)+vct_a(nlev+1+nshift_total)) - start_height) < 0._wp)) THEN
      ! Start height lies at or below the center of the lowermost grid layer:
      ! process acts on entire grid column
      start_height = 0._wp
      iendlev      = nlev
    ELSE
      ! Start height lies somewhere in between
      LEVEL_LOOP: DO jk = 1, nlev  ! (Loops downward)
        jks = jk + nshift_total
        IF ( 0.5_wp * (vct_a(jks)+vct_a(jks+1)) < start_height ) THEN
          iendlev = jk - 1
          EXIT LEVEL_LOOP
        ENDIF
      ENDDO LEVEL_LOOP
      iendlev      = MIN(MAX(1, iendlev), nlev)
      start_height = 0.5_wp * (vct_a(iendlev+nshift_total)+vct_a(iendlev+nshift_total+1))
    ENDIF

  END SUBROUTINE configure_start_height

  !==================================================================================== 

  !>
  !! Configure flow control of events.
  !!
  SUBROUTINE configure_nwp_event( eventStartDateIn,       & !in
    &                             eventEndDateIn,         & !in
    &                             eventIntervalIn,        & !in
    &                             domainStartDate,        & !in
    &                             domainEndDate,          & !in
    &                             basicInterval,          & !in
    &                             eventName,              & !in
    &                             domainName,             & !in
    &                             eventId,                & !in
    &                             eventEnabled,           & !in
    &                             eventObject,            & !inout
    &                             eventSubobject,         & !inout
    &                             optEventReqInit,        & !optin
    &                             optEventInclStart,      & !optin
    &                             optLowerBound4Interval, & !optin
    &                             optUpperBound4Interval, & !optin
    &                             optEventStartDateOut,   & !optout
    &                             optEventEndDateOut,     & !optout
    &                             optEventIntervalOut     ) !optout

    ! In/out variables
    CHARACTER(LEN=*),         INTENT(IN)    :: eventStartDateIn, eventEndDateIn
    REAL(wp),                 INTENT(IN)    :: eventIntervalIn
    TYPE(datetime),           INTENT(IN)    :: domainStartDate, domainEndDate
    REAL(wp),                 INTENT(IN)    :: basicInterval
    CHARACTER(LEN=*),         INTENT(IN)    :: eventName, domainName
    INTEGER,                  INTENT(IN)    :: eventId
    LOGICAL,                  INTENT(IN)    :: eventEnabled
    TYPE(t_phyProcGroup),     INTENT(INOUT) :: eventObject
    TYPE(t_phyProcSlow),      INTENT(INOUT) :: eventSubobject
    LOGICAL,        OPTIONAL, INTENT(IN)    :: optEventReqInit, optEventInclStart
    REAL(wp),       OPTIONAL, INTENT(IN)    :: optLowerBound4Interval, optUpperBound4Interval
    TYPE(datetime), OPTIONAL, INTENT(OUT)   :: optEventStartDateOut, optEventEndDateOut
    REAL(wp),       OPTIONAL, INTENT(OUT)   :: optEventIntervalOut

    ! Local variables
    REAL(wp) :: eventInterval, startDateDiff
    INTEGER  :: istat, tlen
    LOGICAL  :: presentLowerBound, presentUpperBound
    TYPE(datetime) :: eventStartDate, eventEndDate
    TYPE(datetime),  POINTER :: eventDate
    TYPE(timedelta), POINTER :: eventDelta, plusSlack, eventDateCorr
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: eventDeltaString
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':configure_nwp_event'
    CHARACTER(len=*), PARAMETER :: its_adjusted = "It is adjusted to the"

    !---------------------------------------------------------

    !-------------
    ! Preparation
    !-------------
    
    eventDate     => NULL()
    eventDelta    => NULL()
    plusSlack     => NULL()
    eventDateCorr => NULL()



    presentLowerBound = PRESENT(optLowerBound4Interval)
    presentUpperBound = PRESENT(optUpperBound4Interval)
    IF ( presentLowerBound .AND. presentUpperBound) THEN
      IF  (optLowerBound4Interval > optUpperBound4Interval) THEN
        message_text = "Lower bound for interval is greater than upper bound " &
          &       // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      END IF
    ENDIF
    IF (presentLowerBound) THEN
      IF (optLowerBound4Interval < 0._wp) THEN
        message_text = "optLowerBound4Interval has to be positive " &
          &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      ELSEIF (MOD(optLowerBound4Interval, basicInterval) > eps) THEN
        message_text = "optLowerBound4Interval has to be a multiple of basicInterval " &
          &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      ENDIF
    ENDIF
    IF (presentUpperBound) THEN
      IF (optUpperBound4Interval < 0._wp) THEN
        message_text = "optUpperBound4Interval has to be positive " &
          &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      ELSEIF (MOD(optUpperBound4Interval, basicInterval) > eps) THEN
        message_text = "optUpperBound4Interval has to be a multiple of basicInterval " &
          &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      ENDIF
    ENDIF

    !--------------------
    ! Configure interval
    !--------------------

    IF (eventIntervalIn < 0._wp) THEN
      ! A negative time interval is not allowed
      message_text = "Invalid eventIntervalIn " &
        &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
      CALL finish(routine, message_text)
    ELSEIF (eventEnabled .AND. MOD(eventIntervalIn, basicInterval) > eps) THEN
      ! The event interval is not a multiple of the basic interval, 
      ! so we adjust it
      message_text = "WARNING, event interval = "                          &
        & //TRIM(ADJUSTL(real2string(eventIntervalIn, opt_fmt="(F20.3)"))) &
        & //" s is no multiple of basic interval = "                       &
        & //TRIM(ADJUSTL(real2string(basicInterval, opt_fmt="(F20.3)")))//" s "&
        & //"(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
      CALL message(routine, message_text)
      ! Compute a new event interval, which is the multiple of the basic interval 
      ! closest to the original event interval 
      ! (it can be lower or greater than the original value!)
      eventInterval = MAX(1._wp, ANINT(eventIntervalIn / basicInterval)) * basicInterval
      message_text = its_adjusted
      tlen = LEN(its_adjusted)
      IF (eventInterval > eventIntervalIn) THEN
        message_text(tlen+1:) = " GREATER"
        tlen = tlen + 8
      ELSE
        message_text(tlen+1:) = " LOWER"
        tlen = tlen + 6
      ENDIF
      message_text(tlen+1:) = " value of " &
        & //TRIM(ADJUSTL(real2string(eventInterval, opt_fmt="(F20.3)")))//" s!"
      CALL message(routine, message_text)
    ELSE
      eventInterval = eventIntervalIn
    ENDIF
    ! Check for lower and upper bounds for the interval
    IF (presentLowerBound) THEN
      IF (eventInterval < optLowerBound4Interval) THEN 
        message_text = "WARNING, event interval = "                                 &
          & //TRIM(ADJUSTL(real2string(eventInterval, opt_fmt="(F20.3)")))          &
          & //" s is lower than optLowerBound4Interval = "                          &
          & //TRIM(ADJUSTL(real2string(optLowerBound4Interval, opt_fmt="(F20.3)"))) &
          & //" s "//"(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL message(routine, message_text)
        message_text = "It is adjusted to event interval = optLowerBound4Interval"
        CALL message(routine, message_text)
        eventInterval = optLowerBound4Interval
      ENDIF
    ENDIF
    IF (presentUpperBound) THEN
      IF (eventInterval > optUpperBound4Interval) THEN 
        message_text = "WARNING, event interval = "                                 &
          & //TRIM(ADJUSTL(real2string(eventInterval, opt_fmt="(F20.3)")))          &
          & //" s is greater than optUpperBound4Interval = "                        &
          & //TRIM(ADJUSTL(real2string(optUpperBound4Interval, opt_fmt="(F20.3)"))) &
          & //" s "//"(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL message(routine, message_text)
        message_text = "It is adjusted to event interval = optUpperBound4Interval"
        CALL message(routine, message_text)
        eventInterval = optUpperBound4Interval
      ENDIF
    ENDIF
    ! Transform interval into date-time variable
    CALL getPTStringFromMS(INT(eventInterval * 1000._wp, i8), eventDeltaString)
    eventDelta => newTimedelta(eventDeltaString, istat)
    IF (istat /= SUCCESS) THEN
      message_text = "Could not transform event interval into date-time variable " &
        &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
      CALL finish(routine, message_text)
    END IF

    !----------------------
    ! Configure start date
    !----------------------

    IF (LEN_TRIM(eventStartDateIn) == 0) THEN
      ! If the string for the start date of the triggering of events is empty, 
      ! we take the start date of the domain
      eventStartDate = domainStartDate
    ELSE
      ! We try to transform the string into a date-time variable
      eventDate => newDatetime(eventStartDateIn, istat)
      IF (istat /= SUCCESS) THEN
        message_text = "Invalid eventStartDateIn " &
         &          // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      END IF
      eventStartDate = eventDate
      ! Difference between domainStartDate and eventStartDate in seconds
      startDateDiff = 1.0e-3_wp * REAL(getTotalMilliSecondsTimeDelta(eventStartDate - domainStartDate, &
        & domainStartDate), wp)
      IF (startDateDiff > 0._wp .AND. MOD(startDateDiff, basicInterval) > eps) THEN
        startDateDiff = MAX(1._wp, ANINT(startDateDiff / basicInterval)) * basicInterval - startDateDiff
        CALL getPTStringFromMS(INT(startDateDiff * 1000._wp, i8), eventDeltaString)
        eventDateCorr => newTimedelta(eventDeltaString, istat)
        IF (istat /= SUCCESS) THEN
          message_text = "Could not transform event date correction into date-time variable " &
            &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
          CALL finish(routine, message_text)
        END IF
        eventStartDate = eventStartDate + eventDateCorr
        CALL deallocateTimedelta(eventDateCorr)
        message_text = "WARNING, eventStartDate is adjusted to be in phase with cycle time of ICON " &
          & // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL message(routine, message_text)
      ENDIF
      ! If the domain start date happens to be after the event start date, 
      ! the event start date is replaced by the domain start date
      IF (domainStartDate > eventStartDate) THEN 
        eventStartDate = domainStartDate
        message_text = "WARNING, domainStartDate > eventStartDate => eventStartDate = domainStartDate " &
          & // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL message(routine, message_text)
      ENDIF
      CALL deallocateDatetime(eventDate)
    ENDIF

    !--------------------
    ! Configure end date
    !--------------------

    IF(.NOT. eventEnabled) THEN
      ! If the event is disabled, 
      ! we take as end date the start date
      eventEndDate = eventStartDate
    ELSEIF (LEN_TRIM(eventEndDateIn) == 0) THEN
      ! If the string for the end date of the triggering of events is empty, 
      ! we take the end date of the domain
      eventEndDate = domainEndDate
    ELSE
      ! We try to transform the string into a date-time variable
      eventDate => newDatetime(eventEndDateIn, istat)
      IF (istat /= SUCCESS) THEN
        message_text = "Invalid eventEndDateIn " // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL finish(routine, message_text)
      END IF
      eventEndDate = eventDate
      ! If the domain end date happens to be before the event end date, 
      ! the latter is replaced by the former
      IF (eventEndDate > domainEndDate) THEN 
        eventEndDate = domainEndDate
        message_text = "WARNING, eventEndDate > domainEndDate => eventEndDate = domainEndDate " &
          &         // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
        CALL message(routine, message_text)
      ENDIF
      CALL deallocateDatetime(eventDate)
    ENDIF

    ! If after all the start date happens to be after the end date, 
    ! we set the start date equal to the end date
    IF (eventStartDate > eventEndDate) THEN
      eventStartDate = eventEndDate
      message_text = "WARNING, eventStartDate > eventEndDate => eventStartDate = eventEndDate " &
        & // "(event: "//TRIM(eventName)//", domain: "//TRIM(domainName)//")"
      CALL message(routine, message_text)
    ENDIF

    !-----------------
    ! Configure event
    !-----------------

    plusSlack => newTimedelta("PT0S")

    ! Set up event
    ! (see 'src/atm_phy_nwp/mo_phy_events: phyProcBase_initialize' about meaning of input)
    CALL eventSubobject%initialize( name         = eventName,        & !in
      &                             id           = eventId,          & !in
      &                             is_enabled   = eventEnabled,     & !in
      &                             startDate    = eventStartDate,   & !in
      &                             endDate      = eventEndDate,     & !in
      &                             dt           = eventDelta,       & !in
      &                             plusSlack    = plusSlack,        & !in
      &                             optReqInit   = optEventReqInit,  & !optin
      &                             optInclStart = optEventInclStart ) !optin

    ! Add set up event to event group
    CALL eventObject%addToGroup(phyProc = eventSubobject)

    !-------------------
    ! Store event times
    !-------------------

    IF (PRESENT(optEventIntervalOut))  optEventIntervalOut  = eventInterval
    IF (PRESENT(optEventStartDateOut)) optEventStartDateOut = eventStartDate
    IF (PRESENT(optEventEndDateOut))   optEventEndDateOut   = eventEndDate

    !----------
    ! Clean-up
    !----------

    CALL deallocateTimedelta(eventDelta)
    CALL deallocateTimedelta(plusSlack)



  END SUBROUTINE configure_nwp_event

  !==================================================================================== 

  !>
  !! Print message on upper-atmosphere physics configuration.
  !!
  SUBROUTINE print_config_upatmo_physics( jg,                      & !in
    &                                     n_dom,                   & !in
    &                                     iforcing,                & !in
    &                                     lrestart,                & !in
    &                                     lupatmo_phy,             & !in
    &                                     msg_level,               & !in
    &                                     upatmo_phy_config,       & !in
    &                                     upatmo_aes_phy_config,   & !in
    &                                     upatmo_nwp_phy_config    ) !inout

    ! In/out variables
    INTEGER,                   INTENT(IN)    :: jg                      ! Domain index
    INTEGER,                   INTENT(IN)    :: n_dom                   ! Number of domains
    INTEGER,                   INTENT(IN)    :: iforcing                ! Switch for physics package (nwp, aes etc.)
    LOGICAL,                   INTENT(IN)    :: lrestart                ! Switch for restart mode
    LOGICAL,                   INTENT(IN)    :: lupatmo_phy             ! Switch for upper-atmosphere physics (NWP)
    INTEGER,                   INTENT(IN)    :: msg_level               ! Message level
    TYPE(t_upatmo_phy_config), INTENT(IN)    :: upatmo_phy_config       ! Upper-atmosphere configuration 
                                                                        ! with namelist settings
    TYPE(t_upatmo_aes_phy),    INTENT(IN)    :: upatmo_aes_phy_config   ! Upper-atmosphere configuration for AES
    TYPE(t_upatmo_nwp_phy),    INTENT(INOUT) :: upatmo_nwp_phy_config   ! Upper-atmosphere configuration for NWP

    ! Local variables
    INTEGER  :: jgrp, jprc, jgas, cjgadj, tlen, ladj_s, ladj_e, ladj_h
    CHARACTER(LEN=13) :: cremark
    CHARACTER(LEN=80) :: msg_prefix
    CHARACTER(LEN=11) :: cjg
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':print_config_upatmo_physics'
    CHARACTER(len=*), PARAMETER :: &
      hght_msg_sfx = ', => switched off effectively!)'

    !---------------------------------------------------------

    WRITE (cjg, '(i0)') jg
    cjgadj = VERIFY(cjg, " ")

    !-----------------------------------------------------
    !                "Any-case" output
    !-----------------------------------------------------

    IF (msg_level >= imsg_thr%low) THEN
      
      !-----------------------------------------------------
      !              "High-priority" output
      !-----------------------------------------------------
      
      SELECT CASE(iforcing)
        
      CASE(iaes)
        
        !
        !*******************************************************************************
        !                                AES forcing
        !
        !*******************************************************************************
        
      CASE(inwp)
        
        !
        !*******************************************************************************
        !                                 NWP forcing
        !
        !*******************************************************************************
        
        IF (jg == 1 .AND. upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
          
          CALL message(routine, 'Notes on configuration of upper-atmosphere physics:')
          message_text = '* upatmo_nml: nwp_grp_<procgroup>%dt ' &
            & //'(default or NAMELIST input) is not rescaled, if grid_nml: grid_rescale_factor /= 1!'
          CALL message(' ', message_text)
          ! The upper-atmosphere physics are not and cannot be made restart-safe,  
          ! if more than 2 domains are in use 
          ! for reasons discussed in 'src/upper_atmosphere/mo_upatmo_flowevent_utils'
          IF (lrestart .AND. n_dom > 2) THEN
            message_text = '* upper-atmosphere physics are not restart-safe, if more than 2 domains are in use!'
            CALL message(' ', message_text)
          ENDIF
          !
          CALL message(routine, 'Notes on output of (upper-atmosphere) variables:')
          CALL message(' ', '* only variables associated with switched on process groups &
            &can be selected in output_nml!')
          message_text = '* every variable name requires the following prefix: ' &
            & // upatmo_nwp_phy_config%vname_prefix
          CALL message(' ', message_text)
          CALL message(' ', '* heating rates are isobaric, not isochoric &
            &(no matter how they are processed internally)!')
          IF (upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%enabled )) THEN
            CALL message(' ', '* efficiency and scaling factors sclrlw and effrsw are multiplied &
              &to the temperature tendencies from the standard radiation,')
            CALL message(' ', '  so the forcing is NOT ddt_temp_radsw and ddt_temp_radlw, &
              &but effrsw * ddt_temp_radsw and sclrlw * ddt_temp_lw')
          ENDIF  !Radiation enabled?
          IF (n_dom > 1) THEN
            CALL message(' ', '* halo cells typically contain 0s, which might be visible &
              &on the lateral boundary of lon-lat-regridded output of nests')
          ENDIF
          
        ENDIF  !IF (jg == 1 .AND. upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled ))
        
      END SELECT  !SELECT CASE(iforcing)
      
      IF (msg_level >= imsg_thr%med) THEN
        
        !-----------------------------------------------------
        !             "Medium-priority" output
        !-----------------------------------------------------
        
        IF (msg_level >= imsg_thr%high) THEN 

          !-----------------------------------------------------
          !               "Low-priority" output
          !-----------------------------------------------------

          SELECT CASE(iforcing)
            
          CASE(iaes)
            
            !
            !*******************************************************************************
            !                                AES forcing
            !
            !*******************************************************************************

            IF (jg == 1 .AND. upatmo_aes_phy_config%l_enabled) THEN
              ! (For the time being, we take the configuration for NWP, 
              ! where a corresponding configuration for AES is not yet implemented)
              CALL message(routine, 'Info on processes:')
              DO jprc = 1, iUpatmoPrcId%nitem
                message_text = '* '//TRIM(upatmo_nwp_phy_config%prc( jprc )%longname) &
                  & // ' (start height = ' &
                  & //TRIM(ADJUSTL(real2string(upatmo_aes_phy_config%start_height( jprc ), opt_fmt='(F20.1)')))//' m'
                tlen = LEN_TRIM(message_text)
                IF (upatmo_aes_phy_config%iendlev( jprc ) > 0) THEN
                  message_text(tlen+1:) = ')'
                ELSE
                  message_text(tlen+1:) = ', => switched off effectively!)'
                ENDIF
                CALL message(' ', message_text)
              ENDDO  !jprc
            ENDIF  !IF (jg == 1 .AND. upatmo_aes_phy_config%l_enabled)
            
          CASE(inwp)

            !
            !*******************************************************************************
            !                                 NWP forcing
            !
            !*******************************************************************************

            IF (upatmo_nwp_phy_config%l_status( iUpatmoStat%required )) THEN
              msg_prefix = 'upatmo_config('//cjg(cjgadj:)//')%nwp_phy%'
              tlen = LEN_TRIM(msg_prefix)
              !
              WRITE (message_text, '(3a,l7)') 'lupatmo_phy(dom', cjg(cjgadj:), '): ', &
                & lupatmo_phy
              CALL message(routine, message_text)
              !
              WRITE (message_text, '(3a,l7)') 'l_phy_stat(dom', cjg(cjgadj:), '): ', &
                & upatmo_nwp_phy_config%l_phy_stat(iUpatmoPrcStat%enabled)
              CALL message(routine, message_text)
            ENDIF

            IF (upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
              !
              message_text = 'Switched on process groups on dom '//cjg(cjgadj:)//':'
              CALL message(routine, message_text)
              DO jgrp = 1, iUpatmoGrpId%nitem
                IF (upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled )) THEN
                  IF (upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%offline )) THEN
                    cremark = ' (offline'
                  ELSE
                    cremark = ' (interactive'
                  ENDIF
                  WRITE (message_text, '(8a)') '* ', TRIM(upatmo_nwp_phy_config%grp( jgrp )%longname),      &
                       & TRIM(cremark),                                                                       &
                       & ', dt(namelist) = ',                                                                 &
                       & TRIM(ADJUSTL(real2string(upatmo_phy_config%nwp_grp( jgrp )%dt, opt_fmt='(F20.3)'))), &
                       & ' s, dt(used) = ',                                                                   &
                       & TRIM(ADJUSTL(real2string(upatmo_nwp_phy_config%grp( jgrp )%dt, opt_fmt='(F20.3)'))), &
                       & ' s)'
                  CALL message(' ', message_text)
                  CALL message(' ', '   Info on single processes:')
                  DO jprc = 1, iUpatmoPrcId%nitem
                    IF (upatmo_nwp_phy_config%prc( jprc )%igrp == jgrp) THEN
                      WRITE (msg_prefix, '(2(f20.1))') &
                        upatmo_nwp_phy_config%prc( jprc )%start_height, &
                        upatmo_nwp_phy_config%prc( jprc )%end_height
                      ladj_s = 1 + VERIFY(msg_prefix(1:20), ' ')
                      ladj_e = 21 + VERIFY(msg_prefix(21:40), ' ')
                      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc(jprc)%longname)
                      ladj_h = MERGE(1, LEN(hght_msg_sfx), &
                        &            upatmo_nwp_phy_config%prc(jprc)%iendlev>0)
                      WRITE (message_text, '(7a,i0,a,i0,a)') &
                        & '   * ', upatmo_nwp_phy_config%prc(jprc)%longname(1:tlen),                                      &
                        & ' (start height = ',  msg_prefix(ladj_s:20),                                                    &
                        & ' m, end height = ',  msg_prefix(ladj_e:40),                                                    &
                        & ' m, start level = ',                                                                           &
                        & upatmo_nwp_phy_config%prc( jprc )%istartlev,                                                    &
                        & ', end level = ',                                                                               &
                        & upatmo_nwp_phy_config%prc( jprc )%iendlev,                                                      &
                        & hght_msg_sfx(ladj_h:)
                      CALL message(' ', message_text)
                    ENDIF  !IF (upatmo_nwp_phy_config%prc( jprc )%igrp == jgrp)
                  ENDDO  !jprc
                ENDIF  !IF (upatmo_nwp_phy_config%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled ))
              ENDDO  !jgrp
              ! Info on event management (time control) of process groups.
              ! (The call of 'printSetup' is the reason why 
              ! 'upatmo_nwp_phy_config' got the attribute 'INTENT(INOUT)'.)
              CALL upatmo_nwp_phy_config%event_mgmt_grp%printSetup()
            ENDIF  !IF (upatmo_nwp_phy_config%l_phy_stat( iUpatmoPrcStat%enabled ))
            !
            IF (jg == 1 .AND. upatmo_nwp_phy_config%l_gas_stat( iUpatmoGasStat%enabled )) THEN
              !
              CALL message(routine, 'Required radiatively active gases:')
              DO jgas = 1, iUpatmoGasId%nitem
                IF (upatmo_nwp_phy_config%gas( jgas )%l_stat( iUpatmoGasStat%enabled )) THEN
                  tlen = LEN_TRIM(upatmo_nwp_phy_config%gas( jgas )%longname)
                  message_text = '* '//upatmo_nwp_phy_config%gas( jgas )%longname(1:tlen)
                  IF (upatmo_nwp_phy_config%gas( jgas )%imode == iUpatmoGasMode%extdat) THEN
                    message_text(tlen+3:) = ' (external data from file)'
                  END IF
                  CALL message(' ', message_text)
                ENDIF  !IF (upatmo_nwp_phy_config%gas( jgas )%l_stat( iUpatmoGasStat%enabled ))
              ENDDO  !jgas
              !
            ENDIF  !IF (jg == 1 .AND. upatmo_nwp_phy_config%l_gas_stat( iUpatmoGasStat%enabled ))
            !
            ! Info on event management of external data handling
            IF (upatmo_nwp_phy_config%l_extdat_stat( iUpatmoExtdatStat%required )) THEN
              CALL upatmo_nwp_phy_config%event_mgmt_extdat%printSetup()
            ENDIF
            !
            ! Thermodynamic coupling factor
            message_text = 'Thermodynamic coupling factor = ' &
              & //TRIM(ADJUSTL(real2string(upatmo_nwp_phy_config%thermdyn_cpl_fac, opt_fmt='(F6.2)')))
            CALL message(routine, message_text)
            WRITE (message_text, '(a,l7)') 'Consider heat source from heat diffusion: ', &
              & upatmo_phy_config%nwp_ldiss_from_heatdiff
            CALL message(routine, message_text)

          END SELECT  !SELECT CASE(iforcing)

        ENDIF  !imsg_thr%high
      ENDIF  !imsg_thr%med
    ENDIF  !imsg_thr%low

  END SUBROUTINE print_config_upatmo_physics

  !====================================================================================

  !>
  !! Checks, whether current time is before the operational phase.
  !!
  !!
  LOGICAL FUNCTION prc_isBeforeOpPhase( prc,          & !class
    &                                   mtime_current ) !in

    ! In/out variables
    CLASS(t_phy_prc), INTENT(INOUT) :: prc
    TYPE(datetime), INTENT(IN)      :: mtime_current

    !---------------------------------------------------------

    IF (prc%l_stat( iUpatmoPrcStat%enabled )) THEN
      IF (prc%lBeforeOpPhase) THEN
        ! 'lBeforeOpPhase' is initialized with .true.
        ! Once we are after the start date ('lBeforeOpPhase -> .false.'), 
        ! this holds for the rest of the time 
        ! and the following operation becomes superfluous
        prc%lBeforeOpPhase = prc%t_start > mtime_current
      ENDIF
      prc_isBeforeOpPhase = prc%lBeforeOpPhase
    ELSE
      prc_isBeforeOpPhase = .FALSE.
    ENDIF

  END FUNCTION prc_isBeforeOpPhase

  !====================================================================================

  !>
  !! Checks, whether current time is after the operational phase.
  !!
  LOGICAL FUNCTION prc_isAfterOpPhase( prc,          & !class
    &                                  mtime_current ) !in

    ! In/out variables
    CLASS(t_phy_prc), INTENT(INOUT) :: prc
    TYPE(datetime), INTENT(IN)      :: mtime_current

    !---------------------------------------------------------

    IF (prc%l_stat( iUpatmoPrcStat%enabled )) THEN
      IF (.NOT. prc%lAfterOpPhase) THEN
        ! 'lAfterOpPhase' is initialized with .false.
        ! Once we are after the end date ('lAfterOpPhase -> .true.'), 
        ! this holds for the rest of the time 
        ! and the following operation becomes superfluous
        prc%lAfterOpPhase = mtime_current > prc%t_end
      ENDIF
      prc_isAfterOpPhase = prc%lAfterOpPhase
    ELSE
      prc_isAfterOpPhase = .FALSE.
    ENDIF    

  END FUNCTION prc_isAfterOpPhase

  !====================================================================================

  !>
  !! Checks, whether current time is in the operational phase.
  !!
  LOGICAL FUNCTION prc_isInOpPhase( prc,          & !class
    &                               mtime_current ) !in

    ! In/out variables
    CLASS(t_phy_prc),        INTENT(INOUT) :: prc
    TYPE(datetime), POINTER, INTENT(IN)    :: mtime_current

    !---------------------------------------------------------

    IF (prc%l_stat( iUpatmoPrcStat%enabled )) THEN
      ! If the current time is not before 't_start' and not after 't_end', 
      ! it should be within the (closed) start-date-end-date-interval
      prc_isInOpPhase = .NOT. ( prc%isBeforeOpPhase(mtime_current) .OR. &
        &                       prc%isAfterOpPhase(mtime_current)       )
    ELSE
      prc_isInOpPhase = .FALSE.
    ENDIF    

  END FUNCTION prc_isInOpPhase

  !====================================================================================

  !>
  !! Checks, whether current time is before the greatest lower bound of start dates.
  !! (Compare 'prc_isBeforeOpPhase' above.)
  !!
  LOGICAL FUNCTION grp_isBeforeOpPhase( grp,          & !class
    &                                   mtime_current ) !in

    ! In/out variables
    CLASS(t_upatmo_nwp_phy), INTENT(INOUT) :: grp
    TYPE(datetime), INTENT(IN)             :: mtime_current

    !---------------------------------------------------------

    IF (grp%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
      IF (grp%lBeforeOpPhase) THEN
        grp%lBeforeOpPhase = grp%t_start > mtime_current
      ENDIF
      grp_isBeforeOpPhase = grp%lBeforeOpPhase
    ELSE
      grp_isBeforeOpPhase = .FALSE.
    ENDIF

  END FUNCTION grp_isBeforeOpPhase

  !====================================================================================

  !>
  !! Checks, whether current time is after the least upper bound of the end dates.
  !! (Compare 'prc_isAfterOpPhase' above.)
  !!
  LOGICAL FUNCTION grp_isAfterOpPhase( grp,          & !class
    &                                  mtime_current ) !in

    ! In/out variables
    CLASS(t_upatmo_nwp_phy), INTENT(INOUT) :: grp
    TYPE(datetime), POINTER, INTENT(IN)    :: mtime_current

    !---------------------------------------------------------

    IF (grp%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
      IF (.NOT. grp%lAfterOpPhase) THEN
        grp%lAfterOpPhase = mtime_current > grp%t_end
      ENDIF
      grp_isAfterOpPhase = grp%lAfterOpPhase
    ELSE
      grp_isAfterOpPhase = .FALSE.
    ENDIF

  END FUNCTION grp_isAfterOpPhase

  !====================================================================================

  !>
  !! Checks, whether current time is in the operational phase.
  !! (Compare 'prc_isInOpPhase' above.)
  !!
  LOGICAL FUNCTION grp_isInOpPhase( grp,          & !class
    &                               mtime_current ) !in

    ! In/out variables
    CLASS(t_upatmo_nwp_phy), INTENT(INOUT) :: grp
    TYPE(datetime), POINTER, INTENT(IN)    :: mtime_current

    !---------------------------------------------------------

    IF (grp%l_phy_stat( iUpatmoPrcStat%enabled )) THEN
      grp_isInOpPhase = .NOT. ( grp%isBeforeOpPhase(mtime_current) .OR. &
        &                       grp%isAfterOpPhase(mtime_current)       )
    ELSE
      grp_isInOpPhase = .FALSE.
    ENDIF

  END FUNCTION grp_isInOpPhase

END MODULE mo_upatmo_phy_config
