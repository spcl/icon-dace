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

! Interface between 'src/atm_phy_nwp/mo_nh_interface_nwp: nwp_nh_interface'
! and the schemes of upper-atmosphere physics in 'src/upper_atmosphere'

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_upatmo_interface

  USE mo_kind,                   ONLY: wp, vp
  USE mo_exception,              ONLY: finish, message, message_text
  USE mo_impl_constants,         ONLY: SUCCESS, MAX_CHAR_LENGTH, &
    &                                  min_rlcell_int, min_rledge_int
  USE mo_physical_constants,     ONLY: rd_o_cpd, vtmpc1, cpd, cvd
  USE mo_math_constants,         ONLY: dbl_eps
  USE mo_upatmo_impl_const,      ONLY: iUpatmoStat, iUpatmoPrcStat, iUpatmoGrpId,   &
    &                                  iUpatmoTendId, iUpatmoGasStat, iUpatmoGasId, & 
    &                                  iUpatmoTracerId, iUpatmoPrcId
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c, grf_bdywidth_e
  USE mo_model_domain,           ONLY: t_patch
  USE mo_intp_data_strc,         ONLY: t_int_state
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_upatmo_types,           ONLY: t_upatmo, t_upatmo_tend
  USE mo_upatmo_state,           ONLY: prm_upatmo_vec => prm_upatmo   ! WS 2022-02-07, no longer passed as arg
  USE mo_upatmo_config,          ONLY: upatmo_config, upatmo_phy_config
  USE mo_upatmo_phy_config,      ONLY: t_upatmo_nwp_phy
  USE mo_parallel_config,        ONLY: nproma
  USE mo_run_config,             ONLY: iqv
  USE mo_loopindices,            ONLY: get_indices_c, get_indices_e
  USE mo_sync,                   ONLY: sync_patch_array, sync_patch_array_mult, &
    &                                  SYNC_C, SYNC_C1
  USE mtime,                     ONLY: MAX_DATETIME_STR_LEN, &
    &                                  datetime, datetimeToString
  USE mo_phy_events,             ONLY: mtime_ctrl_physics
  USE mo_timer,                  ONLY: timer_start, timer_stop,                    & 
    &                                  timer_upatmo_phy, timer_upatmo_phy_tend,    &
    &                                  timer_upatmo_phy_imf, timer_upatmo_phy_rad, &
    &                                  timer_upatmo_phy_acc
  USE mo_upatmo_phy_iondrag,     ONLY: iondrag
  USE mo_upatmo_phy_vdfmol,      ONLY: vdf_mol
  USE mo_upatmo_phy_fric,        ONLY: fric_heat
  USE mo_upatmo_phy_srbc,        ONLY: srbc_heating
  USE mo_upatmo_phy_euv,         ONLY: euv_heating
  USE mo_upatmo_phy_nlte,        ONLY: nlte_heating
  USE mo_upatmo_phy_no,          ONLY: no_heating
  USE mo_upatmo_phy_diag,        ONLY: update_diagnostic_variables
  USE mo_nh_vert_extrap_utils,   ONLY: sanity_check
  USE mo_util_string,            ONLY: int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nwp_upatmo_interface
  PUBLIC :: nwp_upatmo_update

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_upatmo_interface'

  ! For convenience
  INTEGER, PARAMETER :: ngrp       = iUpatmoGrpId%nitem    ! Number of upper-atmosphere physics groups,
                                                           ! in which single processes are clusterd
  INTEGER, PARAMETER :: nprc       = iUpatmoPrcId%nitem    ! Number of single processes
  INTEGER, PARAMETER :: ntnd       = iUpatmoTendId%nitem   ! Number of variables,
                                                           ! for which tendencies are computed 
                                                           ! excluding the Exner pressure
  INTEGER, PARAMETER :: ntnd_2     = iUpatmoTendId%nitem_2 ! --,,-- including the Exner pressure
  INTEGER, PARAMETER :: ngas       = iUpatmoGasId%nitem    ! Number of radiatively active gases 
  INTEGER, PARAMETER :: ntrc       = iUpatmoTracerId%nitem ! Number of upatmo-affected tracers
  INTEGER, PARAMETER :: igrpIMF    = iUpatmoGrpId%imf      ! Identifier for physics group IMF
  INTEGER, PARAMETER :: igrpRAD    = iUpatmoGrpId%rad      ! Identifier for physics group RAD
  INTEGER, PARAMETER :: igasO3     = iUpatmoGasId%o3       ! Identifier for ozone
  INTEGER, PARAMETER :: igasO2     = iUpatmoGasId%o2       ! Identifier for dioxygen
  INTEGER, PARAMETER :: igasO      = iUpatmoGasId%o        ! Identifier for atomar oxygen
  INTEGER, PARAMETER :: igasCO2    = iUpatmoGasId%co2      ! Identifier for carbon dioxide
  INTEGER, PARAMETER :: igasNO     = iUpatmoGasId%no       ! Identifier for nitric oxide
  INTEGER, PARAMETER :: igasN2     = iUpatmoGasId%n2       ! Identifier for dinitrogen
  INTEGER, PARAMETER :: itendTemp  = iUpatmoTendId%temp    ! Identifier for temperature tendency
  INTEGER, PARAMETER :: itendExner = iUpatmoTendId%exner   ! Identifier for Exner pressure tendency
  INTEGER, PARAMETER :: itendWind  = iUpatmoTendId%wind_h  ! Identifier for u- and v-wind tendencies
  INTEGER, PARAMETER :: itendQx    = iUpatmoTendId%qx      ! Identifier for tracer tendencies
  INTEGER, PARAMETER :: itracerQv  = iUpatmoTracerId%qv    ! Identifier for water vapor

CONTAINS

  !>
  !! Upper-atmosphere-to-NWP interface.
  !!
  SUBROUTINE nwp_upatmo_interface( dt_loc,         &  !in
    &                              mtime_datetime, &  !in
    &                              p_patch,        &  !in 
    &                              p_int_state,    &  !in
    &                              p_metrics,      &  !in
    &                              p_prog,         &  !in
    &                              p_prog_rcf,     &  !in
    &                              p_diag,         &  !in
    &                              prm_nwp_diag,   &  !in
    &                              prm_nwp_tend,   &  !in
    &                              kstart_moist    )  !in

    ! In/out variables
    REAL(wp),                      INTENT(IN)    :: dt_loc          ! Advective time step
    TYPE(datetime),       POINTER, INTENT(IN)    :: mtime_datetime  ! Date/time information
    TYPE(t_patch),        TARGET,  INTENT(IN)    :: p_patch         ! Grid/patch info
    TYPE(t_int_state),             INTENT(IN)    :: p_int_state     ! Horizontal interpolation weights
    TYPE(t_nh_metrics),            INTENT(IN)    :: p_metrics       ! Vertical grid variables
    TYPE(t_nh_prog),               INTENT(IN)    :: p_prog          ! Prognostic variables
    TYPE(t_nh_prog),      TARGET,  INTENT(IN)    :: p_prog_rcf      ! Prognostic variables 
                                                                    ! with reduced calling frequency
    TYPE(t_nh_diag),      TARGET,  INTENT(IN)    :: p_diag          ! Diagnostic variables
    TYPE(t_nwp_phy_diag),          INTENT(IN)    :: prm_nwp_diag    ! Diag vars of NWP-physics
    TYPE(t_nwp_phy_tend),          INTENT(IN)    :: prm_nwp_tend    ! NWP tendencies from slow physics
    INTEGER,                       INTENT(IN)    :: kstart_moist    ! Index of grid layer above which 
                                                                    ! no condensed water phases exist in the model

    ! Local variables
    TYPE(t_upatmo), POINTER :: prm_upatmo           ! WS: Convenience pointer to keep previous code structure

    REAL(wp), POINTER :: ddt_temp_tot(:,:,:), &     ! (nproma,nlev,nblks_c) Pointer to accumulative temperature tendency
      &                  ddt_vn_tot(:,:,:),   &     ! (nproma,nlev,nblks_e) Pointer to accumulative vn-wind tendency
      &                  ddt_qx_tot(:,:,:),   &     ! (nproma,nlev,nblks_c) Pointer to accumulative tracer tendency
      &                  ddt_exner_tot(:,:,:)       ! (nproma,nlev,nblks_c) Pointer to accumulative Exner pressure tendency
    REAL(wp), ALLOCATABLE :: ddt_u_tot(:,:,:), &    ! (nproma,nlev,nblks_c) Accumulative u-wind tendency
      &                      ddt_v_tot(:,:,:)       ! (nproma,nlev,nblks_c) Accumulative v-wind tendency
    ! (The 'gas_...' are allocatable, because they are 
    ! not necessarily required in every call of this subroutine)
    REAL(wp), ALLOCATABLE :: gas_vmr(:,:,:,:),   &  ! (nproma,nlev,nblks_c,ngas) Gas volume mixing ratio
      &                      gas_col(:,:,:,:),   &  ! (nproma,nlev,nblks_c,ngas) Gas column number density
      &                      gas_cumcol(:,:,:,:)    ! (nproma,nlev,nblks_c,ngas) Gas accumulated column density
    REAL(wp) :: cecc, cobld, clonp
    REAL(wp) :: thermdyn_cpl_fac

    INTEGER, POINTER :: iidx(:,:,:), iblk(:,:,:)
    INTEGER, ALLOCATABLE, TARGET :: sunlit_idx( :, : )
    INTEGER, ALLOCATABLE         :: nsunlit( : )
    INTEGER  :: iendlev_prc( nprc ), istartlev_prc( nprc )
    INTEGER  :: iendlev_grp( ngrp ), istartlev_grp( ngrp )
    INTEGER  :: ierror( nprc, p_patch%nblks_c )
    INTEGER  :: iendlev_exner, istartlev_exner
    INTEGER  :: iendlev_vn, istartlev_vn
    INTEGER  :: jg, jb, jk, jc, je, jgrp, jprc, jtnd, jtrc
    INTEGER  :: inewTemp, inewWind, inewQx, inewExner
    INTEGER  :: nlev, nblks_c
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: orbit_type, solvar_type, solvar_data, solcyc_type, yr_perp
    INTEGER  :: istat, error

    LOGICAL  :: lcall_phy( ngrp )
    LOGICAL  :: lgrp_enabled( ngrp ) 
    LOGICAL  :: lgrp_offline( ngrp )
    LOGICAL  :: lgrp_accumulate( ngrp)
    LOGICAL  :: lupdate( ntnd_2 )
    LOGICAL  :: lddt_tot( ntnd_2 )
    LOGICAL  :: lmessage, ltimer, lupdate_gas, lyr_perp
    LOGICAL  :: ldiss_from_heatdiff

    TYPE(t_upatmo_nwp_phy), POINTER :: upatmo_nwp

    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datetime_str
    CHARACTER(LEN=2) :: dom_str

    REAL(wp), PARAMETER :: inv_vtmpc1 = 1._wp / vtmpc1
    REAL(wp), PARAMETER :: cvd_o_cpd  = cvd / cpd
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':nwp_upatmo_interface'
    
    !--------------------------------------------------------------

    ! Domain index
    jg = p_patch%id
    prm_upatmo => prm_upatmo_vec(jg)   ! WS: set convenience pointer

    ltimer = upatmo_config(jg)%l_status( iUpatmoStat%timer )

    IF (ltimer) THEN 
      CALL timer_start(timer_upatmo_phy)
      CALL timer_start(timer_upatmo_phy_tend)
    ENDIF

    !---------------------------------------------------------------------
    !                              Checks
    !---------------------------------------------------------------------

    IF (.NOT. upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%initialized )) THEN
      CALL finish(TRIM(routine), 'Initialization took not yet place')
    ELSEIF (ntrc > 1) THEN
      CALL finish(TRIM(routine), 'Not prepared for more than one tracer (water vapor)')
    ENDIF

    !---------------------------------------------------------------------
    !                             Preparation
    !---------------------------------------------------------------------

    ! Number of vertical levels
    nlev = p_patch%nlev

    ! Number of blocks
    nblks_c = p_patch%nblks_c

    ! Loop boundaries for prognostic domain.
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! Message output desired?
    lmessage = upatmo_config(jg)%l_status( iUpatmoStat%message )

    ! For messages and error handling
    CALL datetimeToString(mtime_datetime, datetime_str)
    dom_str = TRIM(int2string(jg))

    ddt_temp_tot  => NULL()
    ddt_vn_tot    => NULL()
    ddt_qx_tot    => NULL()
    ddt_exner_tot => NULL()
    iidx          => NULL()
    iblk          => NULL()

    ! For brevity
    upatmo_nwp => upatmo_config(jg)%nwp_phy

    DO jgrp = 1, ngrp
      ! Which physics group is enabled for the simulation?
      ! (Required for flow point "Inactivation" below, too.)
      lgrp_enabled( jgrp ) = upatmo_nwp%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled )
      ! Offline-mode? 
      ! (Compute, but not apply tendencies.)
      lgrp_offline( jgrp ) = upatmo_nwp%grp( jgrp )%l_stat( iUpatmoPrcStat%offline )
      ! Accumulate tendencies from group?
      lgrp_accumulate( jgrp ) = lgrp_enabled( jgrp ) .AND. (.NOT. lgrp_offline( jgrp )) .AND. &
        &                       upatmo_nwp%grp( jgrp )%isInOpPhase( mtime_datetime )
      ! Get start and end indices of grid layer range, 
      ! for which tendencies are computed
      iendlev_grp( jgrp )   = upatmo_nwp%grp( jgrp )%iendlev
      istartlev_grp( jgrp ) = upatmo_nwp%grp( jgrp )%istartlev
    ENDDO  !jgrp

    ! Which accumulative tendencies are available in general
    ! (e.g., not in offline-mode)?
    lddt_tot( 1:ntnd ) = upatmo_nwp%l_any_update( 1:ntnd )
    ! 'l_any_update' is available for the basic tendencies only 
    ! (i.e., for wind, temperature and water vapor tendencies). 
    ! The Exner pressure tendency is derived form the temperature 
    ! and water vapor tendencies.
    lddt_tot( itendExner ) = lddt_tot( itendTemp ) .OR. lddt_tot( itendQx )

    ! Check if model time is within the overall start-date-end-date interval 
    ! of the upper-atmosphere physics.
    ! (Rough naming convention for time management switches: 
    !  * 'enabled' -> The process will be executed sometime during the simulation. 
    !                 (A "static" state that does not change during the simulation.)
    !  * 'active'  -> The current time is within the start-date-end-date interval 
    !                 of the process. Within this interval, calls of the process
    !                 are triggered (a "dynamic" state that changes during the simulation:
    !                 'before' -> 'in' -> 'after' the 'active' phase)
    !  * 'event'   -> A call of the process is triggered within its 'active' phase
    !                 (a "dynamic" state).)
    IF (upatmo_nwp%isInOpPhase( mtime_datetime )) THEN

      ! Check if a call of a physics group is due.
      CALL mtime_ctrl_physics( phyProcs      = upatmo_nwp%event_mgmt_grp, & !in
        &                      mtime_current = mtime_datetime,            & !in
        &                      isInit        = .FALSE.,                   & !in
        &                      lcall_phy     = lcall_phy                  ) !inout

      ! The following should only be necessary, if any group is to be called at all
      IF (ANY(lcall_phy( : ))) THEN

        IF (lmessage) CALL message(TRIM(routine), &
          & 'Start computation of tendencies from upper-atmosphere physics on domain '//dom_str)

        ! Determine which variables will have to be updated due to the current call (if any). 
        ! (This information is required for the handling of the accumulative tendencies 
        ! 'prm_upatmo%tend%ddt%<variable>%tot'.)
        lupdate( 1:ntnd_2 ) = .FALSE.  ! Include Exner pressure (-> 'ntnd_2')
        DO jgrp = 1, ngrp
          IF (lcall_phy( jgrp )) THEN
            DO jtnd = 1, ntnd  ! Exclude Exner pressure (-> 'ntnd')
              ! 'l_update' contains the information on the offline mode
              lupdate( jtnd ) = lupdate( jtnd ) .OR. upatmo_nwp%grp( jgrp )%l_update( jtnd )
            ENDDO  !jtnd
          ENDIF  !IF (lcall_phy( jgrp ))
        ENDDO  !jgrp
        ! Update of Exner pressure necessary?
        lupdate( itendExner ) = lupdate( itendTemp ) .OR. lupdate( itendQx ) 

        !---------------------------------------------------------------------
        !                      Update diagnostic variables
        !---------------------------------------------------------------------

        ! Please note that the diagnostic variables are updated only for the prognostic domain. 
        ! A value of zero will likely persist in their halo cells 
        ! throughout the simulation (set initially when they were allocated by 'add_var'). 
        ! (Remember that, if you do a lat-lon output of diagnostic variables on nests.)

        ! Update of radiatively active gases due?
        ! (They are only required for radiation, currently.)
        lupdate_gas = upatmo_nwp%l_gas_stat( iUpatmoGasStat%enabled ) .AND. & ! Gases enabled?
          &           lcall_phy( igrpRAD )                                    ! Call of radiation due?

        IF (lupdate_gas) THEN

          ! Allocate the gas volume mixing ratio and gas-column-number-density-related variables
          ALLOCATE( gas_vmr( nproma, nlev, nblks_c, ngas ),    &
            &       gas_col( nproma, nlev, nblks_c, ngas ),    &
            &       gas_cumcol( nproma, nlev, nblks_c, ngas ), &
            &       STAT=istat                                 )
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of gas_vmr/col/cumcol failed')

          ! Compute the diagnostic variables      
          CALL update_diagnostic_variables( mtime_datetime    = mtime_datetime,        & !in
            &                               lupdate_gas       = lupdate_gas,           & !in
            &                               p_patch           = p_patch,               & !in
            &                               p_prog_rcf        = p_prog_rcf,            & !in
            &                               p_diag            = p_diag,                & !in
            &                               prm_upatmo_diag   = prm_upatmo%diag,       & !inout
            &                               prm_upatmo_tend   = prm_upatmo%tend,       & !inout
            &                               prm_upatmo_extdat = prm_upatmo%extdat,     & !inout
            &                               upatmo_config     = upatmo_config(jg),     & !in
            &                               upatmo_phy_config = upatmo_phy_config(jg), & !in
            &                               nproma            = nproma,                & !in
            &                               opt_kstart_moist  = kstart_moist,          & !optin
            &                               opt_gas_vmr       = gas_vmr,               & !optinout
            &                               opt_gas_col       = gas_col,               & !optinout
            &                               opt_gas_cumcol    = gas_cumcol             ) !optinout 

        ELSE
          
          CALL update_diagnostic_variables( mtime_datetime    = mtime_datetime,        & !in
            &                               lupdate_gas       = lupdate_gas,           & !in
            &                               p_patch           = p_patch,               & !in
            &                               p_prog_rcf        = p_prog_rcf,            & !in
            &                               p_diag            = p_diag,                & !in
            &                               prm_upatmo_diag   = prm_upatmo%diag,       & !inout
            &                               prm_upatmo_tend   = prm_upatmo%tend,       & !inout
            &                               prm_upatmo_extdat = prm_upatmo%extdat,     & !inout
            &                               upatmo_config     = upatmo_config(jg),     & !in
            &                               upatmo_phy_config = upatmo_phy_config(jg), & !in
            &                               nproma            = nproma,                & !in
            &                               opt_kstart_moist  = kstart_moist           ) !optin
          
        ENDIF  !Update of radiatively active gases due?

        !---------------------------------------------------------------------
        !                     Computation of tendencies
        !--------------------------------------------------------------------- 

        ! Loop over single processes
        DO jprc = 1, nprc
          ! Initialize error indicator
          ierror( jprc, : ) = SUCCESS
          ! Get start and end indices of grid layer range, 
          ! for which tendencies are computed.
          ! Please note that the start height of some processes 
          ! might lie above the top of the current domain. 
          ! In this case 'istartlev=1 > iendlev=0', 
          ! and the single process subroutines would return 
          ! right after initializing the tendencies with zero.
          istartlev_prc( jprc ) = upatmo_nwp%prc( jprc )%istartlev
          iendlev_prc( jprc )   = upatmo_nwp%prc( jprc )%iendlev
        ENDDO  !jprc

        ! Loop over groups
        DO jgrp = 1, ngrp

          ! Is group due?
          IF (lcall_phy( jgrp )) THEN

            IF (lmessage) CALL message(TRIM(routine), & 
              & 'Start computation of '//TRIM(upatmo_nwp%grp( jgrp )%longname))

            SELECT CASE( jgrp )

            CASE( igrpIMF )

              !---------------------------------------------------------------------
              !                  Tendencies from physics group:
              !       Ion drag, Molecular diffusion and Frictional heating
              !---------------------------------------------------------------------  

              IF (ltimer) CALL timer_start(timer_upatmo_phy_imf)

              ! Solar activity for ion drag
              solvar_type = upatmo_phy_config(jg)%solvar_type

              ! The original formulation of 'fric_heat' considers "only" the heat source 
              ! from frictional heating. If the following switch is set to .true.
              ! we compute the heat source from the heat diffusion, too.
              ldiss_from_heatdiff = upatmo_phy_config(jg)%nwp_ldiss_from_heatdiff
              
              ! (Index boundaries ('rl_start' etc.) from above should still apply)
              
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, i_startidx, i_endidx, jprc, error) ICON_OMP_GUIDED_SCHEDULE
              DO jb = i_startblk, i_endblk
                
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
                
                !----------
                ! Ion drag
                !----------
                
                jprc = iUpatmoPrcId%iondrag
                
                CALL iondrag( jcs           = i_startidx,                             & !in
                  &           jce           = i_endidx,                               & !in
                  &           kbdim         = nproma,                                 & !in
                  &           klev          = nlev,                                   & !in
                  &           solvar_type   = solvar_type,                            & !in
                  &           psteplen      = dt_loc,                                 & !in
                  &           lat           = p_patch%cells%center(:,jb)%lat,         & !in
                  &           pum1          = p_diag%u(:,:,jb),                       & !in
                  &           pvm1          = p_diag%v(:,:,jb),                       & !in
                  &           pqm1          = p_prog_rcf%tracer(:,:,jb,iqv),          & !in
                  &           grav          = prm_upatmo%diag%grav(:,:,jb),           & !in
                  &           pgeom1        = p_metrics%geopot_agl(:,:,jb),           & !in
                  &           pcp           = prm_upatmo%diag%cpair(:,:,jb),          & !in
                  &           pvom          = prm_upatmo%tend%ddt_u_iondrag(:,:,jb),  & !out 
                  &           pvol          = prm_upatmo%tend%ddt_v_iondrag(:,:,jb),  & !out
                  &           ptte          = prm_upatmo%tend%ddt_temp_joule(:,:,jb), & !out
                  &           opt_istartlev = istartlev_prc( jprc ),                  & !optin
                  &           opt_iendlev   = iendlev_prc( jprc ),                    & !optin
                  &           opt_error     = error                                   ) !optout
                
                IF (error /= SUCCESS) ierror( jprc, jb ) = error
                
                !---------------------
                ! Molecular diffusion
                !---------------------
                
                ! Please note that the upper-atmosphere physics consider only one tracer, 
                ! namely water vapor ('ntrc = 1'). 'vdf_mol' assumes that the 'iqx' 
                ! and 'iUpatmoTracerId%qx' are identical. 
                ! Once more than one tracer has to be treated by upatmo, 
                ! 'iUpatmoTracerId%qx' has to be harmonized with the 'iqx' 
                ! or some index map has to be implemented.
                
                jprc = iUpatmoPrcId%vdfmol
                
                CALL vdf_mol( jcs           = i_startidx,                                                & !in
                  &           jce           = i_endidx,                                                  & !in
                  &           kbdim         = nproma,                                                    & !in
                  &           klev          = nlev,                                                      & !in
                  &           ktracer       = ntrc,                                                      & !in
                  &           psteplen      = dt_loc,                                                    & !in
                  &           ptvm1         = p_diag%tempv(:,:,jb),                                      & !in
                  &           ptm1          = p_diag%temp(:,:,jb),                                       & !in
                  &           pqm1          = p_prog_rcf%tracer(:,:,jb,iqv:iqv),                         & !in
                  &           pum1          = p_diag%u(:,:,jb),                                          & !in
                  &           pvm1          = p_diag%v(:,:,jb),                                          & !in
                  &           papm1         = p_diag%pres(:,:,jb),                                       & !in
                  &           paphm1        = p_diag%pres_ifc(:,:,jb),                                   & !in
                  &           grav          = prm_upatmo%diag%grav(:,:,jb),                              & !in
                  &           amu           = prm_upatmo%diag%amd(:,:,jb),                               & !in
                  &           ptte          = prm_upatmo%tend%ddt_temp_vdfmol(:,:,jb),                   & !out
                  &           pvom          = prm_upatmo%tend%ddt_u_vdfmol(:,:,jb),                      & !out
                  &           pvol          = prm_upatmo%tend%ddt_v_vdfmol(:,:,jb),                      & !out
                  &           pqte          = prm_upatmo%tend%ddt_qx_vdfmol(:,:,jb,itracerQv:itracerQv), & !out
                  &           opt_istartlev = istartlev_prc( jprc ),                                     & !optin
                  &           opt_iendlev   = iendlev_prc( jprc ),                                       & !optin
                  &           opt_error     = error                                                      ) !optout
                
                IF (error /= SUCCESS) ierror( jprc, jb ) = error
                
                !--------------------
                ! Frictional heating
                !--------------------
                
                jprc = iUpatmoPrcId%fric

                CALL fric_heat( jcs                     = i_startidx,                            & !in
                  &             jce                     = i_endidx,                              & !in
                  &             kbdim                   = nproma,                                & !in
                  &             klev                    = nlev,                                  & !in
                  &             ptm1                    = p_diag%temp(:,:,jb),                   & !in
                  &             ptvm1                   = p_diag%tempv(:,:,jb),                  & !in
                  &             pum1                    = p_diag%u(:,:,jb),                      & !in
                  &             pvm1                    = p_diag%v(:,:,jb),                      & !in
                  &             papm1                   = p_diag%pres(:,:,jb),                   & !in
                  &             paphm1                  = p_diag%pres_ifc(:,:,jb),               & !in
                  &             grav                    = prm_upatmo%diag%grav(:,:,jb),          & !in 
                  &             amu                     = prm_upatmo%diag%amd(:,:,jb),           & !in
                  &             cp                      = prm_upatmo%diag%cpair(:,:,jb),         & !in
                  &             ptte_fc                 = prm_upatmo%tend%ddt_temp_fric(:,:,jb), & !out
                  &             opt_istartlev           = istartlev_prc( jprc ),                 & !optin
                  &             opt_iendlev             = iendlev_prc( jprc ),                   & !optin
                  &             opt_ldiss_from_heatdiff = ldiss_from_heatdiff                    ) !optin

              ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

              IF (ltimer) CALL timer_stop(timer_upatmo_phy_imf)
              
            CASE ( igrpRAD )
              
              !---------------------------------------------------------------------
              !                  Tendencies from physics group:
              !                  Radiation and Chemical heating 
              !---------------------------------------------------------------------  

              IF (ltimer) CALL timer_start(timer_upatmo_phy_rad)

              ! Orbit model 
              orbit_type = upatmo_phy_config(jg)%orbit_type
              ! Solar activity
              solvar_type = upatmo_phy_config(jg)%solvar_type 
              ! Solar activity data
              solvar_data = upatmo_phy_config(jg)%solvar_data
              ! Solar cycle 
              solcyc_type = upatmo_phy_config(jg)%solcyc_type
              ! Eccentricity of orbit
              cecc = upatmo_phy_config(jg)%cecc
              ! Obliquity of Earth axis
              cobld = upatmo_phy_config(jg)%cobld
              ! Longitude of perihelion
              clonp = upatmo_phy_config(jg)%clonp
              ! Switch for perpetuation of Earth orbit for year 'yr_perp'
              lyr_perp = upatmo_phy_config(jg)%lyr_perp
              ! Year, for which Earth orbit is perpetuated
              yr_perp = upatmo_phy_config(jg)%yr_perp

              ALLOCATE(sunlit_idx( nproma, nblks_c ), nsunlit( nblks_c ), STAT=istat)
              IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of sunlit_idx and nsunlit failed')

              ! Initialization
              sunlit_idx( :, : ) = 0
              nsunlit( : )       = 0

              ! (Index boundaries ('rl_start' etc.) from above should still apply)
              
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
              DO jb = i_startblk, i_endblk
                
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

                !-------------------------------------------------
                !          Determine sunlit grid columns
                !-------------------------------------------------

                ! Some notes:
                !
                ! * Currently, the only advantage of this step is 
                !   that the following schemes do not have to do this themselves.
                !
                ! * The full potential of this could be exploited 
                !   only if the upper-atmosphere radiation would be computed 
                !   on the radiation grid (if there is a separate one), 
                !   because on such grid sunlit and dark grid columns 
                !   are evenly distributed among the MPI processes.
                !
                ! * The current load balance is bad, because there might be MPI processes
                !   with dark grid columns only, which have to do almost nothing 
                !   in the following (apart from thermal radiation), 
                !   whereas processes in full sunshine would have to do 
                !   the expensive computations for each of their columns.
                !
                ! * Nevertheless, we refrain from computing the radiative processes 
                !   on a radiation grid for the following reasons:
                !   - The required infrastructure is extremely complex 
                !     and beyond our means currently.
                !   - For a proper simulation of, for instance, gravity wave dynamics 
                !     on each domain, it might be necessary to cover processes 
                !     such as "radiative damping" (see, e.g., S. B. Fels (1984) 
                !     "The radiative damping of short vertical scale waves in the mesosphere",
                !     JAS, 41, 1755-1764). 
                !     With the above measure, only a poor coverage of such processes 
                !     would be possible, if at all, since we do not divide the group RAD 
                !     into short-wave and long-wave raditive processes.
                !
                ! * Even if it turns out that the thread has no sunlit grid columns ('nsunlit = 0'),
                !   the following schemes have to be called, because the initialization 
                !   of the temperature tendencies with zero takes place within them.
                !
                ! * In our first implementation, 'sunlit_idx' was allocated for '1:nproma' 
                !   and not for '1:nproma,1:nblks_c', 'nsunlit' was no field at all, 
                !   both were added to the OMP-PRIVATE list
                !   and were computed directly before the call of 'srbc_heating'. 
                !   For yet unidentified reasons that did not work, unfortunately.
                !   It seemed that at least within 'srbc_heating' 'sunlit_idx' 
                !   contained values for 'jc' that belonged to halo cells, 
                !   something which should not have happened. 
                !   For that reason, we decided for the current solution 
                !   and precompute 'sunlit_idx' and 'nsunlit' in a loop 
                !   before the actual loop of the radiation processes.
                !   
                ! * We use 'p_nwp_diag%cosmu0' as the cosine of the solar zenith angle. 
                !   However, in 'src/atm_phy_nwp/mo_nh_interface_nwp: nwp_nh_interface' it seems that 
                !   'zcosmu0' is used instead. Comments on why this is the case are unfortunately sparse, 
                !   so we hope for the best when using 'p_nwp_diag%cosmu0'.

                DO jc = i_startidx, i_endidx
                  IF (prm_nwp_diag%cosmu0(jc,jb) > 0._wp) THEN
                    nsunlit(jb)                = nsunlit(jb) + 1
                    sunlit_idx(nsunlit(jb),jb) = jc
                  ENDIF
                ENDDO  !jc
              ENDDO  !jb
!$OMP END DO

!$OMP DO PRIVATE(jb, jc, i_startidx, i_endidx, jprc, error) ICON_OMP_GUIDED_SCHEDULE
              DO jb = i_startblk, i_endblk
                
                CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

                ! Please note that most of the following radiation processes pose a challenge for restart tests. 
                ! If "gentle" conditions are chosen (a relatively small time step and short simulation period),
                ! most of them can pass the test, except for NLTE 
                ! (at least for the typical compiler settings on the CRAY).
                ! It contains a number of computations that are prone to accumulate round-off errors. 

                !-------------------------------------------------
                ! Heating from Schumann-Runge bands and continuum
                !-------------------------------------------------

                jprc = iUpatmoPrcId%srbc

                CALL srbc_heating( jcs            = i_startidx,                            & !in
                  &                jce            = i_endidx,                              & !in
                  &                kbdim          = nproma,                                & !in
                  &                klev           = nlev,                                  & !in
                  &                ppf            = p_diag%pres(:,:,jb),                   & !in
                  &                prmu0          = prm_nwp_diag%cosmu0(:,jb),             & !in
                  &                am             = prm_upatmo%diag%amd(:,:,jb),           & !in
                  &                cp             = prm_upatmo%diag%cpair(:,:,jb),         & !in
                  &                zo2            = gas_vmr(:,:,jb,igasO2),                & !in
                  &                tto2           = gas_col(:,:,jb,igasO2),                & !in
                  &                heato2         = prm_upatmo%tend%ddt_temp_srbc(:,:,jb), & !out
                  &                solvar_type    = solvar_type,                           & !in
                  &                solvar_data    = solvar_data,                           & !in
                  &                opt_sunlit_idx = sunlit_idx(:,jb),                      & !optin
                  &                opt_nsunlit    = nsunlit(jb),                           & !optin
                  &                opt_istartlev  = istartlev_prc( jprc ),                 & !optin
                  &                opt_iendlev    = iendlev_prc( jprc ),                   & !optin
                  &                opt_error      = error                                  ) !optout

                IF (error /= SUCCESS) ierror( jprc, jb ) = error

                !----------------------------------------------
                ! Heating due to extreme ultraviolet radiation
                !----------------------------------------------            

                jprc = iUpatmoPrcId%euv

                CALL euv_heating( jcs            = i_startidx,                           & !in
                  &               jce            = i_endidx,                             & !in
                  &               kbdim          = nproma,                               & !in
                  &               klev           = nlev,                                 & !in
                  &               prmu0          = prm_nwp_diag%cosmu0(:,jb),            & !in
                  &               zo2            = gas_vmr(:,:,jb,igasO2),               & !in
                  &               zn2            = gas_vmr(:,:,jb,igasN2),               & !in
                  &               zo             = gas_vmr(:,:,jb,igasO),                & !in
                  &               tto2           = gas_cumcol(:,:,jb,igasO2),            & !in
                  &               ttn2           = gas_cumcol(:,:,jb,igasN2),            & !in
                  &               ttox           = gas_cumcol(:,:,jb,igasO),             & !in
                  &               amu            = prm_upatmo%diag%amd(:,:,jb),          & !in
                  &               cp             = prm_upatmo%diag%cpair(:,:,jb),        & !in
                  &               ptte           = prm_upatmo%tend%ddt_temp_euv(:,:,jb), & !out
                  &               this_datetime  = mtime_datetime,                       & !in
                  &               orbit_type     = orbit_type,                           & !in
                  &               solvar_type    = solvar_type,                          & !in
                  &               solcyc_type    = solcyc_type,                          & !in
                  &               cecc           = cecc,                                 & !in
                  &               cobld          = cobld,                                & !in
                  &               clonp          = clonp,                                & !in
                  &               lyr_perp       = lyr_perp,                             & !in
                  &               yr_perp        = yr_perp,                              & !in
                  &               opt_sunlit_idx = sunlit_idx(:,jb),                     & !optin
                  &               opt_nsunlit    = nsunlit(jb),                          & !optin
                  &               opt_istartlev  = istartlev_prc( jprc ),                & !optin
                  &               opt_iendlev    = iendlev_prc( jprc ),                  & !optin
                  &               opt_error      = error                                 ) !optout
                
                IF (error /= SUCCESS) ierror( jprc, jb ) = error

                !--------------------------------------------
                ! Non-LTE infrared cooling due to CO2 and O3
                !--------------------------------------------

                ! Please note that in addition to the temperature tendency 'ddt_temp_nlte', 
                ! a scale factor 'sclrlw' is computed, the temperature tendency 
                ! from the "standard" long-wave radiation 'prm_nwp_tend%ddt_temp_radlw'
                ! has to be multiplied with. In case of the offline-mode this factor should be 1. 

                jprc = iUpatmoPrcId%nlte

                CALL nlte_heating( jcs            = i_startidx,                            & !in
                  &                jce            = i_endidx,                              & !in
                  &                kbdim          = nproma,                                & !in
                  &                klev           = nlev,                                  & !in
                  &                ppf            = p_diag%pres(:,:,jb),                   & !in
                  &                ptf            = p_diag%temp(:,:,jb),                   & !in
                  &                pco2           = gas_vmr(:,:,jb,igasCO2),               & !in
                  &                pco2col        = gas_cumcol(:,:,jb,igasCO2),            & !in
                  &                po3            = gas_vmr(:,:,jb,igasO3),                & !in 
                  &                po2            = gas_vmr(:,:,jb,igasO2),                & !in
                  &                po             = gas_vmr(:,:,jb,igasO),                 & !in
                  &                pn2            = gas_vmr(:,:,jb,igasN2),                & !in
                  &                pam            = prm_upatmo%diag%amd(:,:,jb),           & !in
                  &                pcp            = prm_upatmo%diag%cpair(:,:,jb),         & !in 
                  &                prmu0          = prm_nwp_diag%cosmu0(:,jb),             & !in
                  &                phnlte         = prm_upatmo%tend%ddt_temp_nlte(:,:,jb), & !out
                  &                sclrlw         = prm_upatmo%diag%sclrlw(:,:,jb),        & !out
                  &                opt_sunlit_idx = sunlit_idx(:,jb),                      & !optin
                  &                opt_nsunlit    = nsunlit(jb),                           & !optin
                  &                opt_loffline   = lgrp_offline( jgrp )                   ) !optin

                !------------
                ! NO-heating
                !------------

                jprc = iUpatmoPrcId%no

                CALL no_heating( jcs           = i_startidx,                          & !in
                  &              jce           = i_endidx,                            & !in
                  &              kbdim         = nproma,                              & !in
                  &              klev          = nlev,                                & !in
                  &              zo            = gas_vmr(:,:,jb,igasO),               & !in
                  &              zno           = gas_vmr(:,:,jb,igasNO),              & !in
                  &              cp            = prm_upatmo%diag%cpair(:,:,jb),       & !in
                  &              tm1           = p_diag%temp(:,:,jb),                 & !in
                  &              apm1          = p_diag%pres(:,:,jb),                 & !in
                  &              amu           = prm_upatmo%diag%amd(:,:,jb),         & !in
                  &              ptte          = prm_upatmo%tend%ddt_temp_no(:,:,jb), & !out
                  &              opt_istartlev = istartlev_prc( jprc ),               & !optin
                  &              opt_iendlev   = iendlev_prc( jprc )                  ) !optin

                !------------------
                ! Chemical heating
                !------------------

                ! The temperature tendency 'ddt_temp_chemheat' from chemical heating 
                ! is determined in 'src/upper_atmosphere/mo_upatmo_phy_diag: update_diagnostic_variables' 
                ! called above. The efficiency factor 'effrsw',  the temperature tendency 
                ! from the "standard" short-wave radiation 'prm_nwp_tend%ddt_temp_radsw' 
                ! has to be multiplied with, has been computed once in 
                ! 'src/upper_atmosphere/mo_upatmo_phy_setup: init_upatmo_phy_nwp'.

              ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

              DEALLOCATE(sunlit_idx, nsunlit, STAT=istat)
              IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of sunlit_idx and nsunlit failed')

              IF (ltimer) CALL timer_stop(timer_upatmo_phy_rad)

            CASE DEFAULT

              CALL finish (TRIM(routine), 'Please implement the new physics group into ' &
                & //'src/upper_atmosphere/mo_nwp_upatmo_interface: nwp_upatmo_interface. Thank you!') 

            END SELECT  !Physics groups

            IF (lmessage) CALL message(TRIM(routine), & 
              & 'Finish computation of '//TRIM(upatmo_nwp%grp( jgrp )%longname))

          ENDIF  !Is group due?

        ENDDO  !jgrp

        !---------------------------------------------------------------------
        !                         Check for errors
        !---------------------------------------------------------------------         

        DO jb = 1, nblks_c
          DO jprc = 1, nprc
            IF (ierror( jprc, jb ) /= SUCCESS) THEN
              CALL finish(TRIM(routine), 'Process '//TRIM(upatmo_nwp%prc( jprc )%longname)            &
                & //' on domain '//TRIM(dom_str)//' on block '//TRIM(int2string(jb))//' at time '     &
                & //TRIM(datetime_str)//' reportet error code '//TRIM(int2string(ierror( jprc, jb ))) )
            ENDIF
          ENDDO  !jprc
        ENDDO  !jb

        !---------------------------------------------------------------------
        !                       Accumulate tendencies
        !--------------------------------------------------------------------- 

        ! Some notes:
        !
        ! * The tendencies of groups that operate in offline mode are not accumulated.
        !
        ! * The accumulative tendencies 'prm_upatmo%tend%ddt%<variable>%tot' 
        !   are allocated only if 'upatmo_nwp%l_any_update( iUpatmoTendId%<variable> ) = .TRUE.'. 
        !   (See 'src/upper_atmosphere/mo_upatmo_state: new_upatmo_tend_list'.)
        !
        ! * The conversion of the accumulative temperature and water vapor tendencies 
        !   into an Exner pressure tendency has to take place at the end, 
        !   because we need the most recent state of these tendencies.
        !
        ! * The temperature tendencies from the parameterizations are isobaric, not isochoric, 
        !   and we leave them unchanged! Please remember that, if they go into data output.
        !   Changing the single tendencies would be too computationally expensive 
        !   (or, alternatively, a too invasive intervention within the schemes). 
        !   A transform from isobaric to isochoric (if desired) takes place 
        !   via the factor 'thermdyn_cpl_fac' in the computation of the Exner pressure tendency 
        !   from the accumulative temperature tendency.

        !------------------------------------
        ! New wind tendencies to care about?
        !------------------------------------

        IF (lupdate( itendWind )) THEN

          ALLOCATE( ddt_u_tot( nproma, nlev, nblks_c ), &
            &       ddt_v_tot( nproma, nlev, nblks_c ), &
            &       STAT=istat                          )
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Allocation of ddt_u/v_tot failed')
          
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
            
            ! Initialize accumulative wind tendencies.
            ! (Actually, this is not really necessary 
            ! for the vertical grid layer range, for which tendencies are computed,
            ! because IMF is the only physics group that produces wind tendencies. 
            ! Nevertheless, we prefer the more general procedure, 
            ! for reasons of safety and an easier and less error-prone extensibility.)
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                ! Zonal wind component
                ddt_u_tot(jc,jk,jb) = 0._wp
                ! Meridional wind component
                ddt_v_tot(jc,jk,jb) = 0._wp
              ENDDO  !jc
            ENDDO  !jk  
            
            ! From IMF
            IF (lgrp_accumulate( igrpIMF )) THEN
              
              DO jk = istartlev_grp( igrpIMF ), iendlev_grp( igrpIMF )
                DO jc = i_startidx, i_endidx
                  ! Zonal wind component
                  ddt_u_tot(jc,jk,jb) = ddt_u_tot(jc,jk,jb)                     &
                    &                 + prm_upatmo%tend%ddt_u_vdfmol(jc,jk,jb)  &
                    &                 + prm_upatmo%tend%ddt_u_iondrag(jc,jk,jb) 
                  ! Meridional wind component
                  ddt_v_tot(jc,jk,jb) = ddt_v_tot(jc,jk,jb)                     &
                    &                 + prm_upatmo%tend%ddt_v_vdfmol(jc,jk,jb)  &
                    &                 + prm_upatmo%tend%ddt_v_iondrag(jc,jk,jb)
                ENDDO  !jc
              ENDDO  !jk
              
            ENDIF  !IMF enabled?
            
            ! No wind tendencies from RAD
            
          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
          
          ! For the accumulation of the new wind tendencies 
          ! into the corresponding NWP tendencies in 'nwp_upatmo_update' below, 
          ! the edge-normal component has to be interpolated from du/dt and dv/dt. 
          ! In preparation for this interpolation, we update some halo cells here.
          CALL sync_patch_array_mult(SYNC_C1, p_patch, 2, ddt_u_tot, ddt_v_tot)
          
          ! Swap new and old states and get index of new state
          CALL prm_upatmo%tend%ddt%state( itendWind )%swap()
          inewWind = prm_upatmo%tend%ddt%state( itendWind )%inew()
          ddt_vn_tot => prm_upatmo%tend%ddt%vn( inewWind )%tot
          
          ! Get start and end indices for the grid layer range, 
          ! for which wind tendencies are computed
          istartlev_vn = prm_upatmo%tend%ddt%info( itendWind )%istartlev
          iendlev_vn   = prm_upatmo%tend%ddt%info( itendWind )%iendlev
          
          ! Indices of cells in the neighborhood of an edge
          iidx => p_patch%edges%cell_idx
          iblk => p_patch%edges%cell_blk
          
          ! Loop boundaries for prognostic edges
          rl_start   = grf_bdywidth_e + 1
          rl_end     = min_rledge_int
          i_startblk = p_patch%edges%start_block(rl_start)
          i_endblk   = p_patch%edges%end_block(rl_end)
          
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,     &
              &                 i_startidx, i_endidx, rl_start, rl_end )
            
            ! Interpolate edge-normal component of wind tendency 
            ! from ddt_u_tot and ddt_v_tot.
            ! (The vector interpolation is quite expensive, 
            ! so we limit it to the vertical grid layer range 
            ! for which wind tendencies are computed.)
            DO jk = istartlev_vn, iendlev_vn
              DO je = i_startidx, i_endidx
                ddt_vn_tot(je,jk,jb) = p_int_state%c_lin_e(je,1,jb) *                 &
                  &                  ( ddt_u_tot(iidx(je,jb,1),jk,iblk(je,jb,1)) *    &
                  &                    p_patch%edges%primal_normal_cell(je,jb,1)%v1   &
                  &                  + ddt_v_tot(iidx(je,jb,1),jk,iblk(je,jb,1)) *    &
                  &                    p_patch%edges%primal_normal_cell(je,jb,1)%v2 ) &
                  &                  + p_int_state%c_lin_e(je,2,jb) *                 &
                  &                  ( ddt_u_tot(iidx(je,jb,2),jk,iblk(je,jb,2)) *    &
                  &                    p_patch%edges%primal_normal_cell(je,jb,2)%v1   &
                  &                  + ddt_v_tot(iidx(je,jb,2),jk,iblk(je,jb,2)) *    &
                  &                    p_patch%edges%primal_normal_cell(je,jb,2)%v2   )
              ENDDO  !je
            ENDDO  !jk
            
            DO jk = 1, istartlev_vn - 1
              DO je = i_startidx, i_endidx
                ddt_vn_tot(je,jk,jb) = 0._wp
              ENDDO  !je
            ENDDO  !jk
            
            DO jk = iendlev_vn + 1, nlev
              DO je = i_startidx, i_endidx
                ddt_vn_tot(je,jk,jb) = 0._wp
              ENDDO  !je
            ENDDO  !jk        
            
          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL          
          
          DEALLOCATE(ddt_u_tot, ddt_v_tot, STAT=istat)
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of ddt_u/v_tot failed')

          ddt_vn_tot => NULL()
          iidx       => NULL()
          iblk       => NULL()

          ! Reset loop boundaries to prognostic cells
          rl_start   = grf_bdywidth_c + 1
          rl_end     = min_rlcell_int
          i_startblk = p_patch%cells%start_block(rl_start)
          i_endblk   = p_patch%cells%end_block(rl_end)

        ENDIF  !IF (lupdate( itendWind ))

        !-------------------------------------------
        ! New temperature tendencies to care about?
        !-------------------------------------------

        IF (lupdate( itendTemp )) THEN

          CALL prm_upatmo%tend%ddt%state( itendTemp )%swap()
          inewTemp = prm_upatmo%tend%ddt%state( itendTemp )%inew()
          ddt_temp_tot => prm_upatmo%tend%ddt%temp( inewTemp )%tot

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
            
            ! Initialize accumulative temperature tendency
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                ddt_temp_tot(jc,jk,jb) = 0._wp
              ENDDO  !jc
            ENDDO  !jk
            
            ! From IMF
            IF (lgrp_accumulate( igrpIMF )) THEN
              
              ! Add tendencies from single processes
              DO jk = istartlev_grp( igrpIMF ), iendlev_grp( igrpIMF )
                DO jc = i_startidx, i_endidx
                  ddt_temp_tot(jc,jk,jb) = ddt_temp_tot(jc,jk,jb)                    &
                    &                    + prm_upatmo%tend%ddt_temp_vdfmol(jc,jk,jb) &
                    &                    + prm_upatmo%tend%ddt_temp_fric(jc,jk,jb)   &
                    &                    + prm_upatmo%tend%ddt_temp_joule(jc,jk,jb)
                ENDDO  !jc
              ENDDO  !jk
              
            ENDIF  !IMF enabled?
            
            ! From RAD
            IF (lgrp_accumulate( igrpRAD )) THEN
              
              ! Please note that the upper-atmosphere radiation processing 
              ! includes to modify the temperature tendencies 
              ! from the "standard" radiation in the following way:
              ! * prm_nwp_tend%ddt_temp_radsw  -> effrsw * prm_nwp_tend%ddt_temp_radsw
              ! * prm_nwp_tend%ddt_temp_radlw  -> sclrlw * prm_nwp_tend%ddt_temp_radlw
              ! For this purpose, we add here:
              ! * + (effrsw - 1) * (cv/cp) * prm_nwp_tend%ddt_temp_radsw
              ! * + (sclrlw - 1) * (cv/cp) * prm_nwp_tend%ddt_temp_radlw
              ! This makes life a bit easier for 'nwp_upatmo_update' below, 
              ! but it is an approximation(!), since the update of 'ddt_temp_radsw/radlw' 
              ! in 'src/atm_phy_nwp/mo_nh_interface_nwp: nwp_nh_interface'
              ! takes place every 'dt_fastphy' (named 'dt_loc' here), 
              ! so potentially much more often than the following summation 
              ! of 'ddt_temp_tot' takes place.
              ! However, we accept this inaccuracy for the following reasons:
              ! * The total temperature tendency is transformed into 
              !   an Exner pressure tendency. That is done more straightforwardly 
              !   in this subroutine than in 'nwp_upatmo_update' below.
              ! * We assume that the variations of 'ddt_temp_radsw/radlw' 
              !   are relatively moderate within the update period of RAD. 
              ! A change in value of 'ddt_temp_radsw/radlw' might be strongest 
              ! after calls of the radiative transfer scheme. 
              ! In order to resolve at least this variation in 'ddt_temp_radsw/radlw', 
              ! we enforce that the tendency update period of RAD divides
              ! the calling period of the radiative transfer scheme of NWP evenly
              ! (see 'src/upper_atmosphere/mo_upatmo_phy_config: configure_upatmo_physics').
              ! Since 'ddt_temp_radsw/radlw' are isochoric tendencies 
              ! and not isobaric ones like all other, we convert them with the factor '(cv/cp)'. 
              ! The conversion from isobaric to isochoric tendencies (if desired)
              ! will take place during the transformation of 'ddt_temp_tot' 
              ! into an Exner pressure tendency below.

              DO jk = istartlev_grp( igrpRAD ), iendlev_grp( igrpRAD )
                DO jc = i_startidx, i_endidx
                  ddt_temp_tot(jc,jk,jb) = ddt_temp_tot(jc,jk,jb)                            &
                    &                    + prm_upatmo%tend%ddt_temp_srbc(jc,jk,jb)           &
                    &                    + prm_upatmo%tend%ddt_temp_nlte(jc,jk,jb)           &
                    &                    + prm_upatmo%tend%ddt_temp_euv(jc,jk,jb)            &
                    &                    + prm_upatmo%tend%ddt_temp_no(jc,jk,jb)             &
                    &                    + prm_upatmo%tend%ddt_temp_chemheat(jc,jk,jb)       &
                    &                    + ( prm_upatmo%diag%effrsw(jc,jk,jb) - 1._wp ) *    &
                    &                      cvd_o_cpd * prm_nwp_tend%ddt_temp_radsw(jc,jk,jb) &
                    &                    + ( prm_upatmo%diag%sclrlw(jc,jk,jb) - 1._wp ) *    &
                    &                      cvd_o_cpd * prm_nwp_tend%ddt_temp_radlw(jc,jk,jb)
                ENDDO  !jc
              ENDDO  !jk
              
            ENDIF  !RAD enabled?
            
          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

          ddt_temp_tot => NULL()

        ENDIF  !IF (lupdate( itendTemp ))

        !--------------------------------------
        ! New tracer tendencies to care about?
        !--------------------------------------

        IF (lupdate( itendQx )) THEN

          CALL prm_upatmo%tend%ddt%state( itendQx )%swap()
          inewQx = prm_upatmo%tend%ddt%state( itendQx )%inew()
          
          DO jtrc = 1, ntrc

            ddt_qx_tot => prm_upatmo%tend%ddt%qx( inewQx )%tot(:,:,:,jtrc)
            
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
            DO jb = i_startblk, i_endblk
              
              CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

              ! Initialize accumulative tracer tendency
              DO jk = 1, nlev
                DO jc = i_startidx, i_endidx
                  ddt_qx_tot(jc,jk,jb) = 0._wp
                ENDDO  !jc
              ENDDO  !jk
              
              ! From IMF
              IF (lgrp_accumulate( igrpIMF )) THEN
                
                DO jk = istartlev_grp( igrpIMF ), iendlev_grp( igrpIMF )
                  DO jc = i_startidx, i_endidx
                    ddt_qx_tot(jc,jk,jb) = ddt_qx_tot(jc,jk,jb) &
                      &                  + prm_upatmo%tend%ddt_qx_vdfmol(jc,jk,jb,jtrc)
                  ENDDO  !jc
                ENDDO  !jk
                
              ENDIF  !IMF enabled?

              ! No tracer tendencies from RAD

            ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

            ddt_qx_tot => NULL()
            
          ENDDO  !jtrc

        ENDIF  !IF (lupdate( itendQx ))

        !--------------------------------
        ! New Exner pressure tendencies?
        !--------------------------------
        
        IF (lupdate( itendExner )) THEN

          CALL prm_upatmo%tend%ddt%state( itendExner )%swap()
          inewExner = prm_upatmo%tend%ddt%state( itendExner )%inew()
          ddt_exner_tot => prm_upatmo%tend%ddt%exner( inewExner )%tot
          IF (lddt_tot( itendTemp )) THEN
            inewTemp = prm_upatmo%tend%ddt%state( itendTemp )%inew()
            ddt_temp_tot => prm_upatmo%tend%ddt%temp( inewTemp )%tot
          ENDIF
          IF (lddt_tot( itendQx )) THEN
            inewQx = prm_upatmo%tend%ddt%state( itendQx )%inew()
            ddt_qx_tot => prm_upatmo%tend%ddt%qx( inewQx )%tot(:,:,:,itracerQv)
          ENDIF
          istartlev_exner = prm_upatmo%tend%ddt%info( itendExner )%istartlev
          iendlev_exner   = prm_upatmo%tend%ddt%info( itendExner )%iendlev

          ! For the conversion into an Exner pressure tendency, 
          ! we use the formula:
          !
          ! dexner/dt = R / (cp * theta_v) * [dtemp/dt * (1 + alpha) + temp * dalpha/dt] <=>
          !
          ! dexner/dt = R * (1 + alpha) * temp / (cp * theta_v) * [dtemp/dt / temp + dalpha/dt / (1 + alpha)] <=>
          !
          ! dexner/dt =  R * exner / cp * [dtemp/dt / temp + dalpha/dt / (1 + alpha)],
          !
          ! where we have used that temp * (1 + alpha) / theta_v = temp_v / theta_v = exner.
          ! In addition, we assume that upper-atmosphere tendencies become significant 
          ! only far above the tropopause. So we neglect the contribution of the condensed water phases to alpha:
          !
          ! alpha ~ (Rv / R - 1) * qv.
          !         |__________|
          !               |
          !            = beta
          !
          ! So we can write:
          !
          ! dalpha/dt / (1 + alpha) = beta * dqv/dt / (1 + beta * qv) = dqv/dt / (1 / beta + qv).
          !
          ! An additional reduction of the computational costs could be achieved by applying the approximation:
          !
          ! 1 / (1 + alpha) ~ 1 - alpha
          ! => dqv/dt / (1 / beta + qv) ~ beta^2 * (1 / beta - qv) * dqv/dt.
          !
          ! However, if we can expect alpha << 1 to always hold in the upper atmosphere is not yet clear.

          ! Factor for thermodynamic coupling (either isobaric, isochoric or entropic).
          ! Please note that this factor is only applied to the total temperature tendency. 
          ! The temperature tendencies from the single parameterizations remain the isobaric ones! 
          ! They are not modified, in order to avoid additional computational costs.
          ! In addition, this way of coupling will work only as long as 'prm_upatmo%diag%cpair' 
          ! is a const. field with values equal to "cpd". 
          ! Once this should be changed, the coupling factor has to become a field, too, 
          ! in order to account for the variation of 'cpair', or 'cpair' has to be replaced 
          ! by 'cvair' in the interfaces, where possible!
          thermdyn_cpl_fac = upatmo_nwp%thermdyn_cpl_fac

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

            IF (lddt_tot( itendTemp ) .AND. lddt_tot( itendQx )) THEN

              DO jk = istartlev_exner, iendlev_exner
                DO jc = i_startidx, i_endidx
                  ddt_exner_tot(jc,jk,jb) = rd_o_cpd * p_prog%exner(jc,jk,jb) *              &
                    &                     ( thermdyn_cpl_fac * ddt_temp_tot(jc,jk,jb) /      &
                    &                       p_diag%temp(jc,jk,jb)                            &
                    &                     + ddt_qx_tot(jc,jk,jb) /                           &
                    &                       ( inv_vtmpc1 + p_prog_rcf%tracer(jc,jk,jb,iqv) ) )
                ENDDO  !jc
              ENDDO  !jk

            ELSEIF (lddt_tot( itendTemp )) THEN

              ! Water vapor tendencies are in offline-mode, 
              ! so we do not have to consider them here.

              DO jk = istartlev_exner, iendlev_exner
                DO jc = i_startidx, i_endidx
                  ddt_exner_tot(jc,jk,jb) = rd_o_cpd * p_prog%exner(jc,jk,jb) * thermdyn_cpl_fac * &
                    &                       ddt_temp_tot(jc,jk,jb) / p_diag%temp(jc,jk,jb)
                ENDDO  !jc
              ENDDO  !jk
              
            ELSEIF (lddt_tot( itendQx )) THEN

              ! Temperature tendencies are in offline-mode, 
              ! so we do not have to consider them here.
              
              DO jk = istartlev_exner, iendlev_exner
                DO jc = i_startidx, i_endidx
                  ddt_exner_tot(jc,jk,jb) = rd_o_cpd * p_prog%exner(jc,jk,jb) * &
                    &                       ddt_qx_tot(jc,jk,jb) /              &
                    &                       ( inv_vtmpc1 + p_prog_rcf%tracer(jc,jk,jb,iqv) ) 
                ENDDO  !jc
              ENDDO  !jk
              
            ENDIF
            
          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL          
          
          ddt_exner_tot => NULL()
          IF (lddt_tot( itendTemp )) ddt_temp_tot => NULL()
          IF (lddt_tot( itendQx ))   ddt_qx_tot   => NULL()

        ENDIF  !IF (lupdate( itendExner ))

        !----------
        ! Clean-up
        !----------

        IF (ALLOCATED(gas_vmr)) THEN 
          DEALLOCATE(gas_vmr, STAT=istat)
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of gas_vmr failed')
        ENDIF
        IF (ALLOCATED(gas_col)) THEN 
          DEALLOCATE(gas_col, STAT=istat)
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of gas_col failed')
        ENDIF
        IF (ALLOCATED(gas_cumcol)) THEN 
          DEALLOCATE(gas_cumcol, STAT=istat)
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), 'Deallocation of gas_cumcol failed')
        ENDIF

        IF (lmessage) CALL message(TRIM(routine), &
          & 'Finish computation of tendencies from upper-atmosphere physics on domain '//dom_str)

      ENDIF  !Any physics group to be called at all?

    ENDIF  !In overall start-date-end-date interval?

    ! The following information will be required for 'nwp_upatmo_update'
    DO jtnd = 1, ntnd  ! Exclude Exner pressure!
      prm_upatmo%tend%ddt%info( jtnd )%linActivePhase = .FALSE.
      DO jgrp = 1, ngrp
        ! (Both 'l_update' and 'lgrp_accumulate' contain 
        ! the information on the offline mode)
        IF (upatmo_nwp%grp( jgrp )%l_update( jtnd )) THEN
          prm_upatmo%tend%ddt%info( jtnd )%linActivePhase =        &
            & prm_upatmo%tend%ddt%info( jtnd )%linActivePhase .OR. &
            & lgrp_accumulate( jgrp )
        ENDIF
      ENDDO  !jgrp
    ENDDO  !jtnd
    ! Exner pressure
    prm_upatmo%tend%ddt%info( itendExner )%linActivePhase =       &
      & prm_upatmo%tend%ddt%info( itendTemp )%linActivePhase .OR. &
      & prm_upatmo%tend%ddt%info( itendQx )%linActivePhase
    
    !---------------------------------------------------------------------
    !                           Inactivation
    !--------------------------------------------------------------------- 

    ! Some notes: 
    !
    ! * This point is required, because the active phase of a physics group 
    !   may only be a subset of the time interval of the current experiment, 
    !   depending on the settings for 'upatmo_nml: nwp_grp_<physics group>%t_start/end'. 
    !   Although such functionality is uncommon for the "standard" NWP physics, 
    !   we adopted it from ECHAM-upatmo, in order to keep the namelist control parameters 
    !   for the upper-atmosphere physics relatively uniform. 
    !   In addition, 'upatmo_nml: nwp_grp_<physics group>%t_start' may be of use 
    !   under NWP forcing as a potential means of avoiding crashes during 
    !   the "turbulent" initial adjustment time of the model atmosphere, 
    !   if tendencies from upper-atmosphere physics contribute significantly to the crash.
    !
    ! * So if the current model time 'mtime_datetime' exceeds the namelist-selectable end time 
    !   of a process group, we set the corresponding tendencies to zero (once)  
    !   for reasons of a proper output.

    ! Initialize tendency-wise evaluation of information on inactivation
    prm_upatmo%tend%ddt%info( : )%lafterActivePhase = lddt_tot( : )

    ! Loop over groups
    DO jgrp = 1, ngrp

      IF (lgrp_enabled( jgrp )) THEN

        ! Check if we have left the start-date-end-date intervall 
        IF ( upatmo_nwp%grp( jgrp )%isAfterOpPhase( mtime_datetime ) .AND. &
          &  .NOT. upatmo_nwp%grp( jgrp )%l_stat( iUpatmoPrcStat%afterActivePhase ) ) THEN

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, jtrc) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

            SELECT CASE( jgrp )

            CASE( igrpIMF )

              !---------------------------------------------------------------------
              !                     Inactivate physics group:
              !       Ion drag, Molecular diffusion and Frictional heating
              !---------------------------------------------------------------------

              DO jk = 1, nlev
                DO jc = i_startidx, i_endidx
                  ! Temperature
                  prm_upatmo%tend%ddt_temp_vdfmol(jc,jk,jb) = 0._wp
                  prm_upatmo%tend%ddt_temp_fric(jc,jk,jb)   = 0._wp
                  prm_upatmo%tend%ddt_temp_joule(jc,jk,jb)  = 0._wp
                  ! Zonal wind component
                  prm_upatmo%tend%ddt_u_vdfmol(jc,jk,jb)  = 0._wp
                  prm_upatmo%tend%ddt_u_iondrag(jc,jk,jb) = 0._wp
                  ! Meridional wind component
                  prm_upatmo%tend%ddt_v_vdfmol(jc,jk,jb)  = 0._wp
                  prm_upatmo%tend%ddt_v_iondrag(jc,jk,jb) = 0._wp
                ENDDO  !jc
              ENDDO  !jk
              ! Tracer
              ! (The current loop order is inefficient, 
              ! but a special treatment of the tracers 
              ! involves too much code overhead. 
              ! In addition, this takes place only once, 
              ! so the inefficiency should be bearable.)
              DO jtrc = 1, ntrc
                DO jk = 1, nlev
                  DO jc = i_startidx, i_endidx
                    prm_upatmo%tend%ddt_qx_vdfmol(jc,jk,jb,jtrc) = 0._wp
                  ENDDO  !jc
                ENDDO  !jk
              ENDDO  !jtrc
              
            CASE ( igrpRAD ) 

              !---------------------------------------------------------------------
              !                     Inactivate physics group:
              !                  Radiation and Chemical heating 
              !---------------------------------------------------------------------                
                    
              DO jk = 1, nlev
                DO jc = i_startidx, i_endidx
                  prm_upatmo%tend%ddt_temp_srbc(jc,jk,jb)     = 0._wp
                  prm_upatmo%tend%ddt_temp_nlte(jc,jk,jb)     = 0._wp
                  prm_upatmo%tend%ddt_temp_euv(jc,jk,jb)      = 0._wp
                  prm_upatmo%tend%ddt_temp_no(jc,jk,jb)       = 0._wp
                  prm_upatmo%tend%ddt_temp_chemheat(jc,jk,jb) = 0._wp
                  ! The scale factor 'sclrlw' from non-LTE infrared cooling 
                  ! and the efficiency factor 'effrsw' from chemical heating are reset to 1 
                  prm_upatmo%diag%sclrlw(jc,jk,jb) = 1._wp
                  prm_upatmo%diag%effrsw(jc,jk,jb) = 1._wp
                ENDDO  !jc
              ENDDO  !jk
              
            END SELECT !SELECT CASE( jgrp )

          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

          ! Indicate that we have left the active phase
          upatmo_nwp%grp( jgrp )%l_stat( iUpatmoPrcStat%afterActivePhase ) = .TRUE.

        ENDIF  !Left active phase?

      ENDIF  !IF (lgrp_enabled( jgrp ))

      ! Evaluate information on inactivation tendency-wise
      ! (required for 'nwp_upatmo_update')
      DO jtnd = 1, ntnd  ! Exclude Exner pressure!
        IF (upatmo_nwp%grp( jgrp )%l_update( jtnd )) THEN
          prm_upatmo%tend%ddt%info( jtnd )%lafterActivePhase =         &
            & prm_upatmo%tend%ddt%info( jtnd )%lafterActivePhase .AND. &
            & upatmo_nwp%grp( jgrp )%l_stat( iUpatmoPrcStat%afterActivePhase )
        ENDIF
      ENDDO  !jtnd

    ENDDO  !jgrp

    ! Evaluate information on inactivation for Exner pressure
    prm_upatmo%tend%ddt%info( itendExner )%lafterActivePhase =        &
      & prm_upatmo%tend%ddt%info( itendTemp )%lafterActivePhase .AND. &
      & prm_upatmo%tend%ddt%info( itendQx )%lafterActivePhase 

    DO jtnd = 1, ntnd_2  ! Include Exner pressure
      IF ( prm_upatmo%tend%ddt%info( jtnd )%lafterActivePhase .AND. &
        &  prm_upatmo%tend%ddt%state( jtnd )%lunlockable()          ) THEN
        ! If a tendency left its active phase, we do a final swap, 
        ! in order to simplify a proper treatment of this event 
        ! in 'nwp_upatmo_update' below.
        ! (Effects of the argument 'optFinal = .TRUE.': 
        !  * 'lfinal()' will return .true. until the next call of 'clear()'.
        !  * The next call of 'clear()' in 'nwp_upatmo_update' 
        !    does not unlock 'swap()' for a next call, 
        !    but locks 'swap()' permanently 
        !    and a call to 'lunlockable()' will return .false.
        ! Effect of the argument 'optCountAsUpdate = .FALSE.':
        !  * In general, 'lupdated()' returns .true. after a call of 'swap()'
        !    and until the next call of 'clear()'. 
        !    The abovementioned argument, however, makes 'lupdated()' return .false.
        CALL prm_upatmo%tend%ddt%state( jtnd )%swap(optFinal = .TRUE., optCountAsUpdate = .FALSE.)
      ENDIF
    ENDDO  !jtnd    
      
    !---------------------------------------------------------------------
    !                 Inactivate accumulative tendencies
    !                 after all processes went to sleep
    !---------------------------------------------------------------------

    IF (upatmo_nwp%l_phy_stat( iUpatmoPrcStat%enabled )) THEN 

      IF ( upatmo_nwp%isAfterOpPhase( mtime_datetime ) .AND. &
        &  .NOT. upatmo_nwp%l_phy_stat( iUpatmoPrcStat%afterActivePhase ) ) THEN

        ! The accumulative tendencies cannot be output 
        ! and there is no more call of 'nwp_upatmo_update' 
        ! after the overall active phase, 
        ! so currently there is no reason to set the accumulative tendencies to zero. 

        upatmo_nwp%l_phy_stat( iUpatmoPrcStat%afterActivePhase ) = .TRUE.

        IF (lmessage) CALL message(TRIM(routine), &
          & 'Upper-atmosphere physics shut down on domain '//dom_str)

      ENDIF  !Left active phase?

    ENDIF  !Any upatmo physics enabled?

    ! Clean-up
    upatmo_nwp => NULL()

    IF (ltimer) THEN 
      CALL timer_stop(timer_upatmo_phy_tend)
      CALL timer_stop(timer_upatmo_phy)
    ENDIF

  END SUBROUTINE nwp_upatmo_interface

  !==================================================================================== 

  !>
  !! Interface to accumulate the upper-atmosphere physics tendencies.
  !!
  SUBROUTINE nwp_upatmo_update( lslowphys,       &  !in
    &                           lradheat,        &  !in
    &                           lturb,           &  !in
    &                           dt_loc,          &  !in
    &                           p_patch,         &  !inout
    &                           p_prog_rcf,      &  !inout
    &                           p_diag           )  !inout

    ! In/out variables
    LOGICAL,                     INTENT(IN)    :: lslowphys        ! Any slow physics called in NWP interface?
    LOGICAL,                     INTENT(IN)    :: lradheat         ! Radiative heating called?
    LOGICAL,                     INTENT(IN)    :: lturb            ! Turbulence scheme called?
    REAL(wp),                    INTENT(IN)    :: dt_loc           ! Advective time step (fast-physics time step)
    TYPE(t_patch),       TARGET, INTENT(INOUT) :: p_patch          ! Grid/patch info
    TYPE(t_nh_prog),             INTENT(INOUT) :: p_prog_rcf       ! Prog vars (with red. calling frequency for tracers)
    TYPE(t_nh_diag),             INTENT(INOUT) :: p_diag           ! Diagnostic variables

    ! Local variables
    TYPE(t_upatmo_tend), POINTER               :: prm_upatmo_tend  ! Upper-atmosphere physics tendencies
    REAL(wp) :: mv( nproma, p_patch%nlev), dmv( nproma, p_patch%nlev)
    REAL(wp) :: sum_dmv( nproma ), mv_tot( nproma )

    INTEGER  :: istartlev( ntnd_2 ), iendlev( ntnd_2 )
    INTEGER  :: jg, jb, jk, jc, je, jtnd
    INTEGER  :: nswap, inew, iold
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx 

    LOGICAL  :: lupdated_upatmo( ntnd_2 ), linActivePhase( ntnd_2 )
    LOGICAL  :: ladd( ntnd_2 ), lsubtract( ntnd_2 ), lafterActivePhase( ntnd_2 )
    LOGICAL  :: lqvconstrained( p_patch%nblks_c )
    LOGICAL  :: lmessage, ltimer, lupdated_nwp, lpassed

    REAL(wp), POINTER :: ddt_exner_tot_old(:,:,:), ddt_exner_tot_new(:,:,:), &
      &                  ddt_vn_tot_old(:,:,:), ddt_vn_tot_new(:,:,:),       &
      &                  ddt_qv_tot(:,:,:)

    CHARACTER(LEN=2)               :: dom_str
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: sanity_check_msg

    REAL(wp), PARAMETER :: eps = ABS(dbl_eps) * 1000._wp
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':nwp_upatmo_update'
    
    !--------------------------------------------------------------

    ! Some notes:
    !
    ! * Within the overall start-date-end-date interval of the upper-atmosphere physics 
    !   this subroutine has to be called every fast-physics time step for the following reasons: 
    !   - We may otherwise miss switches of 'lslowphys' and 'lradheat', 
    !     which would lead to misinterpretations of the content of the NWP tendencies.
    !   - Specific humidity has to be updated every fast-physics time step.
    !
    ! * The reason for the tendencies to be accumulated here and not directly 
    !   in 'mo_nh_interface: nwp_nh_interface', from where this subroutine is called, 
    !   is mainly that the tendency fields for the upper-atmosphere physics, 
    !   'prm_upatmo_tend%ddt%<variable>%tot' would be allocated only if the respective 
    !   physics group ('imf' or 'rad') would be switched on. 
    !   Adding the upatmo tendencies in 'nwp_nh_interface' under these circumstances, 
    !   would have meant many expensive additional if-queries and clearness-reducing 
    !   code doublings. With this subroutine, we can keep most of the costs on 
    !   the upatmo side, although this comes at a price:
    !   - If new upper-atmosphere tendencies have to be integrated into the NWP tendencies, 
    !     but the NWP tendencies have not been updated in 'nwp_nh_interface' 
    !     since the last call of this subroutine, we have to subtract 
    !     the old upper-atmosphere tendencies before we can add the new ones.
    !   - Even if there are no new upper-atmosphere tendencies to be integrated 
    !     into the NWP tendencies at this call, but the NWP tendencies have been updated 
    !     in 'nwp_nh_interface' the contribution from the upper-atmosphere physics is lost.
    !     So we have to reintegrate the most recent state of the upper-atmosphere tendencies.

    !---------------------------------------------------------------------
    !                           Preparation
    !---------------------------------------------------------------------

    ! Domain index
    jg = p_patch%id
    prm_upatmo_tend => prm_upatmo_vec(jg)%tend

    ltimer = upatmo_config(jg)%l_status( iUpatmoStat%timer )

    IF (ltimer) THEN 
      CALL timer_start(timer_upatmo_phy)
      CALL timer_start(timer_upatmo_phy_acc)
    ENDIF

    ! Message output desired?
    lmessage = upatmo_config(jg)%l_status( iUpatmoStat%message )

    ! For messages
    dom_str = TRIM(int2string(jg))

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Start integration of upatmo tendencies into NWP tendencies on domain '//dom_str)

    ddt_exner_tot_old => NULL()
    ddt_exner_tot_new => NULL()
    ddt_vn_tot_old    => NULL()
    ddt_vn_tot_new    => NULL()
    ddt_qv_tot        => NULL()

    DO jtnd = 1, ntnd_2

      ! Is current time in time period, within which tendency experiences updates?
      linActivePhase( jtnd ) = prm_upatmo_tend%ddt%info( jtnd )%linActivePhase

      ! If the event of inactivation occured for a tendency, 
      ! we may have to subtract previously integrated upper-atmosphere tendencies, 
      ! but no more (re)integration is required. 
      ! The case that 'lafterActivePhase = .true.' should occur only once:
      ! indeed, once 'ddt%info%lafterActivePhase = .true.' it stays like this, 
      ! but 'ddt%state%lfinal()' will return .true. only once,
      ! because if it does so, the call of 'ddt%state%clear()' at the end of this subroutine 
      ! will make 'ddt%state%lfinal()' return .false. again and will lock  
      ! 'ddt%state%swap()' permanently (i.e., any further call of 'swap' 
      ! will lead to a call of 'finish').
      ! (Please note that this information is not equivalent to '.NOT. linActivePhase')
      lafterActivePhase( jtnd ) = prm_upatmo_tend%ddt%info( jtnd )%lafterActivePhase .AND. &
        &                         prm_upatmo_tend%ddt%state( jtnd )%lfinal()

      ! Number of swaps up to now
      nswap = prm_upatmo_tend%ddt%state( jtnd )%nswap()

      ! The total tendencies of which variables have been updated by upper-atmosphere physics?
      ! ('state%lupdated' has been set by the call of 'state%swap' in 'nwp_upatmo_interface' above.)
      lupdated_upatmo( jtnd ) = prm_upatmo_tend%ddt%state( jtnd )%lupdated()

      ! The total tendencies of which variables have been updated by NWP physics?
      lupdated_nwp = lslowphys
      IF ((jtnd == itendTemp) .OR. (jtnd == itendExner)) lupdated_nwp = lupdated_nwp .OR. lradheat
      IF (jtnd == itendWind) lupdated_nwp = lupdated_nwp .AND. lturb

      ! Which tendencies have to be integrated or reintegrated into the total NWP tendencies?
      ladd( jtnd ) = lupdated_upatmo( jtnd ) .OR. &
        &            (linActivePhase( jtnd ) .AND. (nswap > 0) .AND. lupdated_nwp)

      ! Do we have to subtract upper-atmosphere tendencies 
      ! previously integrated into the total NWP tendencies 
      ! and still contained therein?
      lsubtract( jtnd ) = (lupdated_upatmo( jtnd ) .OR. lafterActivePhase( jtnd )) .AND. &
        &                 (nswap > 1) .AND. (.NOT. lupdated_nwp)

      ! Start and end indices of the grid layer range 
      ! for which tendencies are computed
      istartlev( jtnd ) = prm_upatmo_tend%ddt%info( jtnd )%istartlev
      iendlev( jtnd )   = prm_upatmo_tend%ddt%info( jtnd )%iendlev

    ENDDO  !jtnd_2

    !---------------------------------------------------------------------
    !                         Wind tendencies
    !---------------------------------------------------------------------

    IF (ladd( itendWind ) .OR. lsubtract( itendWind )) THEN

      IF (lsubtract( itendWind )) THEN
        iold = prm_upatmo_tend%ddt%state( itendWind )%iold()
        ddt_vn_tot_old => prm_upatmo_tend%ddt%vn( iold )%tot
      ENDIF
      IF (ladd( itendWind )) THEN
        inew = prm_upatmo_tend%ddt%state( itendWind )%inew()
        ddt_vn_tot_new => prm_upatmo_tend%ddt%vn( inew )%tot     
      ENDIF

      ! Loop boundaries for prognostic domain.
      rl_start   = grf_bdywidth_e + 1
      rl_end     = min_rledge_int
      i_startblk = p_patch%edges%start_block(rl_start)
      i_endblk   = p_patch%edges%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, je, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_e( p_patch, jb, i_startblk, i_endblk,     &
          &                 i_startidx, i_endidx, rl_start, rl_end )

        !---------------------------------------------------------------------
        !                   Subtract old wind tendencies
        !---------------------------------------------------------------------

        IF (lsubtract( itendWind )) THEN 

          DO jk = istartlev( itendWind ), iendlev( itendWind )
            DO je = i_startidx, i_endidx
              p_diag%ddt_vn_phy(je,jk,jb) = p_diag%ddt_vn_phy(je,jk,jb) &
                &                         - REAL(ddt_vn_tot_old(je,jk,jb), vp)
            ENDDO  !je
          ENDDO  !jk

        ENDIF

        !---------------------------------------------------------------------
        !                      Add new wind tendencies
        !---------------------------------------------------------------------

        IF (ladd( itendWind )) THEN

          DO jk = istartlev( itendWind ), iendlev( itendWind )
            DO je = i_startidx, i_endidx
              p_diag%ddt_vn_phy(je,jk,jb) = p_diag%ddt_vn_phy(je,jk,jb) &
                &                         + REAL(ddt_vn_tot_new(je,jk,jb), vp)
            ENDDO  !je
          ENDDO  !jk

        ENDIF

      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      IF (lsubtract( itendWind )) ddt_vn_tot_old => NULL()
      IF (ladd( itendWind ))      ddt_vn_tot_new => NULL()

    ENDIF  !Process wind tendencies?

    ! The rest of the tendencies live in cells
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !---------------------------------------------------------------------
    !                     Exner pressure tendencies
    !---------------------------------------------------------------------

    IF (ladd( itendExner ) .OR. lsubtract( itendExner )) THEN

      IF (lsubtract( itendExner )) THEN
        iold = prm_upatmo_tend%ddt%state( itendExner )%iold()
        ddt_exner_tot_old => prm_upatmo_tend%ddt%exner( iold )%tot
      ENDIF
      IF (ladd( itendExner )) THEN
        inew = prm_upatmo_tend%ddt%state( itendExner )%inew()
        ddt_exner_tot_new => prm_upatmo_tend%ddt%exner( inew )%tot
      ENDIF
      
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk
        
        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
          &                 i_startidx, i_endidx, rl_start, rl_end )

        !---------------------------------------------------------------------
        !               Subtract old Exner pressure tendencies
        !---------------------------------------------------------------------
        
        IF (lsubtract( itendExner )) THEN 

          DO jk = istartlev( itendExner ), iendlev( itendExner )
            DO jc = i_startidx, i_endidx
              p_diag%ddt_exner_phy(jc,jk,jb) = p_diag%ddt_exner_phy(jc,jk,jb) &
                &                            - REAL(ddt_exner_tot_old(jc,jk,jb), vp)
            ENDDO  !jc
          ENDDO  !jk

        ENDIF

        !---------------------------------------------------------------------
        !                  Add new Exner pressure tendencies
        !---------------------------------------------------------------------

        IF (ladd( itendExner )) THEN

          DO jk = istartlev( itendExner ), iendlev( itendExner )
            DO jc = i_startidx, i_endidx
              p_diag%ddt_exner_phy(jc,jk,jb) = p_diag%ddt_exner_phy(jc,jk,jb) &
                &                            + REAL(ddt_exner_tot_new(jc,jk,jb), vp)
            ENDDO  !jc
          ENDDO  !jk

        ENDIF
        
      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      IF (lsubtract( itendExner )) ddt_exner_tot_old => NULL()
      IF (ladd( itendExner ))      ddt_exner_tot_new => NULL()
      
    ENDIF  !Process Exner pressure tendencies?

    !---------------------------------------------------------------------
    !                       Water vapor tendencies
    !---------------------------------------------------------------------

    ! Water vapor tendencies are not accumulated. 
    ! Rather the tracer mass fractions are updated in 
    ! 'src/atm_phy_nwp/mo_util_phys: tracer_add_phytend' 
    ! with the individual tendencies from the slow physics  
    ! (i.e., 'prm_nwp_tend%ddt_tracer_pconv' from convection)
    ! every fast-physics time step. 
    ! This means that we can proceed in the same way here. 
    ! We do not have to care about a "reintegration".
    ! Every call of this subroutine, we update specific humidity 
    ! by 'ddt_qv_tot * dt_fasphys'.

    ! Initialize switch that indicates, 
    ! if one of the following constraints on qv took effect
    lqvconstrained( : ) = .FALSE.

    IF (linActivePhase( itendQx )) THEN

      ! Currently, water vapor tendencies may result from molecular diffusion
      ! (see 'vdf_mol' above). Neglecting advective (reversible) fluxes, 
      ! the diffusion equation for the water vapor mass mv in a grid cell may read:
      !    dmv/dt = -divergence_of_diffusive_fluxes <=> d(qv*m)/dt = ... | * 1/V
      ! => d(qv*rho)/dt = -divergence_of_diffusive_fluxes/V | - qv*drho/dt = 0
      ! => rho*dqv/dt = ... | * 1/rho 
      ! => dqv/dt = -divergence_of_diffusive_fluxes/(V*rho), 
      ! where m denotes the total air mass in a grid cell, V is the grid cell volume, 
      ! and rho = m/V. We applied the common modeling assumption that there are 
      ! no diffusive (irreversible) fluxes of the total air mass, which means here:
      ! dm/dt = 0 | * 1/V => drho/dt = 0. 
      ! Unlike dmv/dt = ..., dqv/dt = ... is not in budget form, 
      ! so that its discretized formulation does not per se satisfy mass conservation 
      ! unless it is appropriately adjusted to the discretized formulation of 
      ! the continuity equation for the total air mass. 
      ! In addition, 'ddt_qv_tot * dt_fasphys' might lead to negative values for 'qv'.
      ! To overcome these problems, we proceed as in 'tracer_add_phytend':
      ! * We write: 
      !   dqv = ddt_qv * dt_fasphys = d(mv/m) = d[mv/(md + mv)] = dmv/m - mv*(dmv + dmd)/m^2 <=>
      !       = dmv/m * (1 - qv) - qv*dmd/m,
      !   where md denotes the mass of dry air.
      ! * From our above assumption that there is no diffusive flux of total air mass, 
      !   it follows that:
      !   dm = dmd + dmv = 0  <=>  dmd = -dmv  <=>  dqv = dmv/m.
      ! * The scheme of molecular diffusion assumes a vanishing vertical diffusive flux 
      !   of water vapor at the model top and at the grid layer interface, 
      !   below which no more tendencies are computed. 
      !   So the water vapor mass in this column should be conserved:
      !   Sum_jk=1_to_iendlev(dmv_jk) =: Dmv =! 0.
      ! * In addition, negative values for the water vapor mass in a grid cell are not allowed:
      !   mv + dmv >=! 0.
      ! * In order to satisfy the latter constraint, we set a lower limit to dmv:
      !   dmv' = Max(dmv, -mv)
      ! * Whenever this limitation takes effect, the mass -(mv + dmv) is created. 
      !   We add this mass to Dmv:
      !   Dmv' = Dmv - Sum_jk=1_to_iendlev(Min(mv + dmv, 0))
      ! * Finally, in order to satisfy the mass conservation, 
      !   we distribute the negative excess mass -Dmv' among the grid cells:
      !   dmv'' = dmv' - mv*Dmv'/mv_tot,
      !   where mv_tot = Sum_jk=1_to_iendlev(mv_jk).
      !   The factor mv/mv_tot effects that the more water vapor mass a grid cell contains, 
      !   the larger the fraction of the excess mass that is subtracted. 
      !   By this measure, we try to reduce the likelihood for mv + dmv'' < 0 to occur (again). 
      !   Now, at least Sum_jk=1_to_iendlev(dmv''_jk) = 0 should be satisfied.
      ! * The new specific humidity qv_new follows from:
      !   qv_new = (mv + dmv'')/m .
      ! * We stress that this approach is a brute-force method 
      !   without any physical background whatsoever! It is completely arbitrary. 
      !   Any other approach that achieves the global mass conservation
      !   and positive semidefiniteness of the local mass 
      !   is equally justified/unjustified.
      ! * The computational costs of the current approach are relatively high. 
      !   So if a less costly approach comes to your mind, 
      !   please replace this approach by yours.
      !
      ! IMPORTANT: if one of the above constraints would take effect 
      ! and would change ddt_qv_tot effectively, we do not update ...
      ! - ... ddt_qv_vdfmol, a field which might have been selected for data output, 
      !   since we would be unable to disentangle the contribution from single processes 
      !   (without undue effort) once further processes, in addition to vdfmol,
      !   would contribute to ddt_qv_tot, 
      ! - ... the Exner pressure tendency for efficiency reasons.

      inew = prm_upatmo_tend%ddt%state( itendQx )%inew()
      ddt_qv_tot => prm_upatmo_tend%ddt%qx( inew )%tot(:,:,:,itracerQv)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, mv, dmv, sum_dmv, mv_tot) ICON_OMP_GUIDED_SCHEDULE
      DO jb = i_startblk, i_endblk

        CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
          &                 i_startidx, i_endidx, rl_start, rl_end )

        ! Initialize sum of changes in water vapor mass
        ! and total water vapor mass
        DO jc = i_startidx, i_endidx
          sum_dmv(jc) = 0._wp
          mv_tot(jc)  = 0._wp
        ENDDO  !jc

        ! Change in water vapor mass
        DO jk = istartlev( itendQx ), iendlev( itendQx )
          DO jc = i_startidx, i_endidx
            mv(jc,jk)   = p_prog_rcf%tracer(jc,jk,jb,iqv) * p_diag%airmass_new(jc,jk,jb)
            dmv(jc,jk)  = dt_loc * ddt_qv_tot(jc,jk,jb) * p_diag%airmass_new(jc,jk,jb)
            sum_dmv(jc) = sum_dmv(jc) + dmv(jc,jk)
            mv_tot(jc)  = mv_tot(jc)  + mv(jc,jk)
            IF (mv(jc,jk) < -dmv(jc,jk)) THEN 
              sum_dmv(jc)        = sum_dmv(jc) - ( mv(jc,jk) + dmv(jc,jk) )
              dmv(jc,jk)         = -mv(jc,jk)
              lqvconstrained(jb) = .TRUE.
            ENDIF
          ENDDO  !jc
        ENDDO  !jk

        ! Avoid division by zero
        DO jc = i_startidx, i_endidx
          IF (mv_tot(jc) < eps) mv_tot(jc) = eps
        ENDDO  !jc

        ! Distribute the excess mass among the cells 
        ! and compute the new value of specific humidity
        DO jk = istartlev( itendQx ), iendlev( itendQx )
          DO jc = i_startidx, i_endidx
            dmv(jc,jk) = dmv(jc,jk) + mv(jc,jk) * sum_dmv(jc) / mv_tot(jc)
            p_prog_rcf%tracer(jc,jk,jb,iqv) = ( mv(jc,jk) + dmv(jc,jk) ) / p_diag%airmass_new(jc,jk,jb)
          ENDDO  !jc
        ENDDO  !jk       

      ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

      ! 'lsanitycheck' can be set in 'src/namelists/mo_upatmo_nml'
      IF (upatmo_phy_config(jg)%lsanitycheck) THEN
        
        ! Is updated specific humidity still positive-definite?
        CALL sanity_check ( p_patch     = p_patch,                      &  !in
          &                 state       = p_prog_rcf%tracer(:,:,:,iqv), &  !in
          &                 bound       = 0._wp,                        &  !in
          &                 keys        = 'cell,lower',                 &  !in
          &                 lpassed     = lpassed,                      &  !out
          &                 opt_slev    = istartlev( itendQx ),         &  !optin
          &                 opt_elev    = iendlev( itendQx ),           &  !optin
          &                 opt_message = sanity_check_msg              )  !optout

        IF (.NOT. lpassed) THEN
          message_text = 'Updated specific humidity is not positive-definite: ' &
            & //TRIM(sanity_check_msg)//' Please reconsider its computation.'
          CALL finish(TRIM(routine), TRIM(message_text))
        ENDIF
        
      ENDIF  !IF (upatmo_phy_config(jg)%lsanitycheck)

      ! Update halo cells.
      ! (This call is the reason why 'p_patch' got the attribute 'INTENT(INOUT)' 
      ! instead of 'INTENT(IN)', which would otherwise have been sufficient.)
      CALL sync_patch_array(SYNC_C, p_patch, p_prog_rcf%tracer(:,:,:,iqv))

    ENDIF  !Process water vapor tendencies?

    !---------------------------------------------------------------------
    !                          Postprocessing
    !---------------------------------------------------------------------

    DO jtnd = 1, ntnd_2  ! Include Exner pressure ('ntnd_2' instead of 'ntnd')

      ! Unlock 'swap' for next call in 'nwp_upatmo_interface' 
      ! - a step, without which the next call of 'swap' would lead to a call of 'finish' - 
      ! or lock 'swap' permanently, if the active phase of a tendency is over
      IF (lupdated_upatmo( jtnd ) .OR. lafterActivePhase( jtnd )) &
        & CALL prm_upatmo_tend%ddt%state( jtnd )%clear()

    ENDDO

    IF (lmessage) THEN
      message_text = 'Finish integration of upatmo tendencies into NWP tendencies on domain' &
        & //TRIM(dom_str)
      IF (ANY(lqvconstrained(:))) THEN
        message_text = TRIM(message_text)//' (constraint: qv + dqv_upatmo >=! 0 took effect!)'
      ENDIF
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

    IF (ltimer) THEN 
      CALL timer_stop(timer_upatmo_phy_acc)
      CALL timer_stop(timer_upatmo_phy)
    ENDIF

  END SUBROUTINE nwp_upatmo_update

END MODULE mo_nwp_upatmo_interface
