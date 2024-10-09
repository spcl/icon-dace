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

! Update diagnostic variables of upper-atmosphere physics.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_upatmo_phy_diag

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message, message_text
  USE mo_impl_constants,         ONLY: MAX_CHAR_LENGTH, &
    &                                  min_rlcell_int
  USE mo_physical_constants,     ONLY: cpd, cpv
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
  USE mo_upatmo_impl_const,      ONLY: iUpatmoStat, iUpatmoGasId,         &
    &                                  iUpatmoGasMode, iUpatmoExtdatStat, &
    &                                  iUpatmoExtdatId
  USE mo_run_config,             ONLY: iqv, iqc, iqi 
  USE mo_model_domain,           ONLY: t_patch
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag
  USE mo_upatmo_types,           ONLY: t_upatmo_diag, t_upatmo_tend, &
    &                                  t_upatmo_extdat
  USE mo_upatmo_config,          ONLY: t_upatmo_config
  USE mo_upatmo_phy_config,      ONLY: t_upatmo_phy_config, t_nwp_gas, &
    &                                  t_upatmo_nwp_phy
  USE mo_loopindices,            ONLY: get_indices_c
  USE mtime,                     ONLY: datetime
  USE mo_phy_events,             ONLY: mtime_ctrl_physics
  USE mo_o3_util,                ONLY: o3_pl2ml
  USE mo_upatmo_extdat,          ONLY: update_upatmo_extdat_nwp
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_upatmo_phy_diag
  USE mo_util_string,            ONLY: int2string
  USE mo_nh_vert_extrap_utils,   ONLY: sanity_check

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: update_diagnostic_variables

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_phy_diag'

  ! For convenience
  INTEGER, PARAMETER :: next        = iUpatmoExtdatId%nitem ! Number of external data types
  INTEGER, PARAMETER :: ngas        = iUpatmoGasId%nitem    ! Number of radiatively active gases 
  INTEGER, PARAMETER :: igasDiag    = iUpatmoGasId%diag     ! Identifier for diagnostic gas
  INTEGER, PARAMETER :: imodeConst  = iUpatmoGasMode%const  ! Identifier for constant gas mode
  INTEGER, PARAMETER :: imodeExtdat = iUpatmoGasMode%extdat ! Identifier for external-data gas mode
  INTEGER, PARAMETER :: iextGases   = iUpatmoExtdatId%gases ! Identifier for external gas data
  
CONTAINS

  !>
  !! Subroutine to update the upper-atmosphere diagnositc variables 
  !! under NWP forcing.
  !!
  SUBROUTINE update_diagnostic_variables( mtime_datetime,    &  !in
    &                                     lupdate_gas,       &  !in
    &                                     p_patch,           &  !in
    &                                     p_prog_rcf,        &  !in
    &                                     p_diag,            &  !in
    &                                     prm_upatmo_diag,   &  !inout
    &                                     prm_upatmo_tend,   &  !inout
    &                                     prm_upatmo_extdat, &  !inout
    &                                     upatmo_config,     &  !in
    &                                     upatmo_phy_config, &  !in
    &                                     nproma,            &  !in
    &                                     opt_linit,         &  !optin
    &                                     opt_lrestart,      &  !optin
    &                                     opt_kstart_moist,  &  !optin
    &                                     opt_gas_vmr,       &  !optinout
    &                                     opt_gas_col,       &  !optinout
    &                                     opt_gas_cumcol     )  !optinout 

    ! In/out variables
    TYPE(datetime),                      INTENT(IN)    :: mtime_datetime          ! Date/time information
    LOGICAL,                             INTENT(IN)    :: lupdate_gas             ! Is update of gases due?
    TYPE(t_patch),             TARGET,   INTENT(IN)    :: p_patch                 ! Grid/patch info
    TYPE(t_nh_prog),           TARGET,   INTENT(IN)    :: p_prog_rcf              ! Prog vars with reduced calling frequency
    TYPE(t_nh_diag),                     INTENT(IN)    :: p_diag                  ! Diagnostic variables
    TYPE(t_upatmo_diag),                 INTENT(INOUT) :: prm_upatmo_diag         ! Upper-atmosphere diagnostic fields
    TYPE(t_upatmo_tend),       TARGET,   INTENT(INOUT) :: prm_upatmo_tend         ! Upper-atmosphere physics tendencies
    TYPE(t_upatmo_extdat),     TARGET,   INTENT(INOUT) :: prm_upatmo_extdat       ! Upper-atmosphere external data
    TYPE(t_upatmo_config),     TARGET,   INTENT(IN)    :: upatmo_config           ! General upper-atmosphere configuration
    TYPE(t_upatmo_phy_config), TARGET,   INTENT(IN)    :: upatmo_phy_config       ! Upper-atmosphere physics configuration 
                                                                                  ! with namelist settings
    INTEGER,                             INTENT(IN)    :: nproma                  ! Blocking length
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: opt_linit               ! Switch for initialization mode
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: opt_lrestart            ! Switch for resumed simulation
    INTEGER,                   OPTIONAL, INTENT(IN)    :: opt_kstart_moist        ! Index of grid layer above which 
                                                                                  ! no condensed water phases exist in the model
    REAL(wp),                  OPTIONAL, INTENT(INOUT) :: opt_gas_vmr(:,:,:,:)    ! Gas volume mixing ratio
                                                                                  ! (nproma,nlev,nblks_c,ngas)
    REAL(wp),                  OPTIONAL, INTENT(INOUT) :: opt_gas_col(:,:,:,:)    ! Gas column number density
                                                                                  ! (nproma,nlev,nblks_c,ngas) 
    REAL(wp),                  OPTIONAL, INTENT(INOUT) :: opt_gas_cumcol(:,:,:,:) ! Gas accumulated column density
                                                                                  ! (nproma,nlev,nblks_c,ngas) 
  
    ! Local variables    
    REAL(wp), POINTER :: qv(:,:,:), qc(:,:,:), qi(:,:,:)
    REAL(wp) :: mmr, mmr2vmr, mass2mol
        
    INTEGER  :: jg, jb, jk, jc, jgas, jext
    INTEGER  :: nlev, kstart_moist, nlev_gas
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx 
    INTEGER  :: igas_ext

    LOGICAL  :: lupdate_extdat( next )
    LOGICAL  :: linit, lrestart
    LOGICAL  :: lgas_vmr, lgas_col, lgas_cumcol
    LOGICAL  :: lmessage, ltimer, lpassed

    CHARACTER(LEN=2)               :: dom_str
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: sanity_check_msg

    TYPE(t_nwp_gas),        POINTER :: upatmo_nwp_gas(:)
    TYPE(t_upatmo_nwp_phy), POINTER :: upatmo_nwp

    TYPE(datetime), TARGET  :: lastactive
    TYPE(datetime), POINTER :: lastactive_ptr

    REAL(wp), PARAMETER :: cpv_m_cpd   = cpv - cpd
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':update_diagnostic_variables'
    
    !--------------------------------------------------------------

    !---------------------------------------------------------------------
    !                     Checks and preparation
    !---------------------------------------------------------------------

    ltimer = upatmo_config%l_status( iUpatmoStat%timer )
    IF (ltimer) CALL timer_start(timer_upatmo_phy_diag)

    IF (PRESENT(opt_linit)) THEN
      linit = opt_linit
    ELSE
      linit = .FALSE.
    ENDIF

    IF (PRESENT(opt_lrestart)) THEN
      ! A restart is only considered, 
      ! if this is an initialization call
      lrestart = opt_lrestart .AND. linit
    ELSE
      lrestart = .FALSE.
    ENDIF

    lgas_vmr    = PRESENT(opt_gas_vmr)
    lgas_col    = PRESENT(opt_gas_col)
    lgas_cumcol = PRESENT(opt_gas_cumcol)

    IF (lupdate_gas .AND. (.NOT. linit) .AND. (.NOT. (lgas_vmr .AND. lgas_col .AND. lgas_cumcol))) THEN
      CALL finish(TRIM(routine), 'Optional arguments are missing')
    ELSEIF (linit .AND. (lgas_vmr .OR. lgas_col .OR. lgas_cumcol)) THEN
      CALL finish(TRIM(routine), 'No input of optional arguments allowed for initialization call')
    ENDIF

    ! Domain index
    jg = p_patch%id
    
    ! Message output desired?
    lmessage = upatmo_config%l_status( iUpatmoStat%message )
    dom_str  = TRIM(int2string(jg))
    IF (lmessage .AND. linit) THEN 
      CALL message(TRIM(routine), &
        & 'Start initialization of diagnostic variables for upper-atmosphere physics on domain '//dom_str)
    ELSEIF (lmessage) THEN
      CALL message(TRIM(routine), &
        & 'Start update of diagnostic variables for upper-atmosphere physics on domain '//dom_str)
    ENDIF

    ! For convenience
    upatmo_nwp => upatmo_config%nwp_phy

    IF (upatmo_nwp%l_extdat_stat( iUpatmoExtdatStat%required )) THEN

      !-------------------------------------
      ! Time interpolation of external data
      !-------------------------------------
      
      ! * For chemical heating this means that this 
      !   is the actual place where the heating tendencies 
      !   'prm_upatmo_tend%ddt_temp_chemheat' are determined. 
      
      ! * The external gas data are interpolated in time 
      !   and onto the horizontal grid of ICON, 
      !   but they are still on the original pressure levels 
      !   of the external data. 
      !   Their vertical interpolation is done below.

      ! * In case of a restart things get ugly. 
      !   The fields contained by 'prm_upatmo_extdat' 
      !   are not (and cannot be) handeled by 'add_var'. 
      !   As a consequence they cannot be saved in the restart file. 
      !   So if this is a resumed simulation, 
      !   we have to interpolate the external data to that datetime 
      !   when the last interpolation took place 
      !   in the simulation resumed.

      IF (.NOT. lrestart) THEN

        ! Any update of external data due?
        CALL mtime_ctrl_physics( phyProcs      = upatmo_nwp%event_mgmt_extdat, & !in
          &                      mtime_current = mtime_datetime,               & !in
          &                      isInit        = .FALSE.,                      & !in
          &                      lcall_phy     = lupdate_extdat                ) !inout
        
        ! If this is the initialization call, 
        ! we force a first interpolation, 
        ! in order to fill the gas fields 'prm_upatmo_diag%gas' 
        ! with reasonable values, for they might have been selected
        ! for data output at time = 0.
        IF (linit) lupdate_extdat(:) = upatmo_nwp%extdat(:)%l_stat( iUpatmoExtdatStat%required )
        
        IF (ANY(lupdate_extdat(:))) THEN
          
          CALL update_upatmo_extdat_nwp( mtime_datetime    = mtime_datetime,    & !in
            &                            p_patch           = p_patch,           & !in
            &                            prm_upatmo_extdat = prm_upatmo_extdat, & !inout
            &                            prm_upatmo_tend   = prm_upatmo_tend,   & !inout
            &                            upatmo_config     = upatmo_config,     & !in
            &                            lupdate           = lupdate_extdat,    & !in
            &                            lmessage          = lmessage           ) !in
        
        ENDIF  !IF (ANY(lupdate_extdat(:)))

      ELSE  ! Resumed simulation

        ! The last interpolation of the different types of external data 
        ! may have taken place at different times. 
        ! That's why we have to loop over the types of external data.
        lupdate_extdat(:) = .FALSE.
        DO jext = 1, next

          lupdate_extdat(jext) = upatmo_nwp%extdat(jext)%l_stat( iUpatmoExtdatStat%required )

          IF (lupdate_extdat(jext)) THEN
          
            ! We do not have to worry about the case 
            ! that the restart point in time coincides 
            ! with a trigger time for one of the external data, 
            ! because this subroutine will be called again 
            ! before the computation of the physics tendencies.
            lastactive     = upatmo_nwp%extdat(jext)%event%getLastActive()
            lastactive_ptr => lastactive

            CALL update_upatmo_extdat_nwp( mtime_datetime    = lastactive_ptr,    & !in
              &                            p_patch           = p_patch,           & !in
              &                            prm_upatmo_extdat = prm_upatmo_extdat, & !inout
              &                            prm_upatmo_tend   = prm_upatmo_tend,   & !inout
              &                            upatmo_config     = upatmo_config,     & !in
              &                            lupdate           = lupdate_extdat,    & !in
              &                            lmessage          = lmessage           ) !in

            lastactive_ptr => NULL()

          ENDIF  !IF (lupdate_extdat(jext))

        ENDDO  !jext

      ENDIF  !IF (.NOT. lrestart)

    ENDIF  !External data required?

    ! In case of a resumed simulation the computations below
    ! are unnecessary, because the fields contained by 'prm_upatmo_diag' 
    ! would have been stored in the restart file.
    IF (lrestart) THEN
      upatmo_nwp => NULL()
      IF (lmessage .AND. linit) &
        & CALL message(TRIM(routine), &
        & 'Finish initialization of diagnostic variables for upper-atmosphere physics on domain '//dom_str)
      IF (ltimer) CALL timer_stop(timer_upatmo_phy_diag)
      RETURN
    ENDIF

    ! Number of vertical levels
    nlev = p_patch%nlev

    IF (PRESENT(opt_kstart_moist)) THEN
      kstart_moist = MAX(1, MIN(opt_kstart_moist, nlev))
    ELSE
      kstart_moist = 0
    ENDIF

    ! For convenience
    qv => p_prog_rcf%tracer(:,:,:,iqv)
    qc => p_prog_rcf%tracer(:,:,:,iqc)
    qi => p_prog_rcf%tracer(:,:,:,iqi)

    ! Loop boundaries for prognostic domain. 
    ! (The upper-atmosphere physics make no use 
    ! of the halo cells.)
    rl_start   = grf_bdywidth_c + 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
    DO jb = i_startblk, i_endblk
      
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
        &                 i_startidx, i_endidx, rl_start, rl_end )

      !---------------------------------------------------------------------
      !                    Some quantities of the air
      !--------------------------------------------------------------------- 
      
      ! For grid layers without condensed water phases ...
      DO jk = 1, kstart_moist
        DO jc = i_startidx, i_endidx      
          ! Mass of dry air [kg/m2]:
          ! mass of air minus water vapor mass
          prm_upatmo_diag%mdry(jc,jk,jb) = ( 1._wp - qv(jc,jk,jb) ) * p_diag%airmass_new(jc,jk,jb)
          ! Heat capacity of moist air at constant pressure [J/K/kg]
          prm_upatmo_diag%cpair(jc,jk,jb) = cpd + cpv_m_cpd * qv(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk

      ! ... and for grid layers with condensed water phases
      DO jk = kstart_moist + 1, nlev
        DO jc = i_startidx, i_endidx      
          prm_upatmo_diag%mdry(jc,jk,jb) = ( 1._wp - qv(jc,jk,jb)     &
            &                                      - qc(jc,jk,jb)     &
            &                                      - qi(jc,jk,jb) ) * & 
            &                              p_diag%airmass_new(jc,jk,jb)
          prm_upatmo_diag%cpair(jc,jk,jb) = cpd + cpv_m_cpd * qv(jc,jk,jb)
        ENDDO  !jc
      ENDDO  !jk

      !---------------------------------------------------------------------
      !                     Initialize diagnostic gas
      !--------------------------------------------------------------------- 

      IF (lupdate_gas) THEN

        ! The diagnostic gas has to be initialized with 1, 
        ! because of its iterative computation below 
        ! (diagnostic gas mass = air mass - all other gas masses)
        DO jk = 1, nlev
          DO jc = i_startidx, i_endidx              
            prm_upatmo_diag%gas(jc,jk,jb,igasDiag) = 1._wp
          ENDDO  !jc
        ENDDO  !jk            
        
      ENDIF  !IF (lupdate_gas)
     
    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

    !---------------------------------------------------------------------
    !                  Update radiatively active gases
    !--------------------------------------------------------------------- 

    IF (lupdate_gas) THEN

      ! For brevity
      upatmo_nwp_gas => upatmo_phy_config%nwp_gas

      ! Loop over gases
      DO jgas = 1, ngas

        ! Skip diagnostic gas
        IF (jgas /= igasDiag) THEN

          !---------------------------------------------------------------------
          !                   Treat the different gas modes
          !--------------------------------------------------------------------- 

          IF (upatmo_nwp_gas(jgas)%imode == imodeExtdat) THEN

            !----------------------------------------
            ! Preparation for external-data gas mode
            !----------------------------------------

            ! Get extdat-internal gas index
            igas_ext = prm_upatmo_extdat%mapgasid2indx( jgas )
            IF (igas_ext < 1) CALL finish(TRIM(routine), 'Invalid extdat gas index')

            ! Number of grid layers in external data file
            nlev_gas = prm_upatmo_extdat%gas( igas_ext )%nlev

          ELSEIF (linit .AND. upatmo_nwp_gas(jgas)%imode == imodeConst) THEN

            !-------------------
            ! Constant gas mode
            !-------------------
            
            ! This has to be done only once during initialization.

            ! Mass mixing ratio of gas from settings in upatmo_nml [kg/kg]
            mmr = upatmo_nwp_gas( jgas )%mmr
            
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
            DO jb = i_startblk, i_endblk
            
              CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
                &                 i_startidx, i_endidx, rl_start, rl_end )
              
              DO jk = 1, nlev
                DO jc = i_startidx, i_endidx  
                  ! Amount of gas in a grid cell in terms of kg/kg
                  prm_upatmo_diag%gas(jc,jk,jb,jgas) = mmr
                ENDDO  !jc
              ENDDO  !jk               
            ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL            
            
          ENDIF  !Gas mode

!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
              &                 i_startidx, i_endidx, rl_start, rl_end )
            
            IF (upatmo_nwp_gas(jgas)%imode == imodeExtdat) THEN
              
              !------------------------
              ! External-data gas mode
              !------------------------
              
              ! The external gas data are provided on pressure levels. 
              ! Here, we interpolate them to the grid layers of ICON. 
              CALL o3_pl2ml( jcs            = i_startidx,                                         & !in
                &            jce            = i_endidx,                                           & !in
                &            kbdim          = nproma,                                             & !in
                &            nlev_pres      = nlev_gas,                                           & !in
                &            klev           = nlev,                                               & !in
                &            pfoz           = prm_upatmo_extdat%gas( igas_ext )%lev,              & !in
                &            phoz           = prm_upatmo_extdat%gas( igas_ext )%lev_half,         & !in
                &            ppf            = p_diag%pres(:,:,jb),                                & !in
                &            pph            = p_diag%pres_ifc(:,:,jb),                            & !in
                &            o3_time_int    = prm_upatmo_extdat%gas_interm( igas_ext )%p(:,:,jb), & !in
                &            o3_clim        = prm_upatmo_diag%gas(:,:,jb,jgas),                   & !out
                &            opt_o3_initval = 0._wp                                               ) !optin              
              
            ENDIF  !IF (upatmo_nwp_gas(jgas)%imode == imodeExtdat)
            
            !---------------------------
            ! Update the diagnostic gas 
            !---------------------------

            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx
                ! We assume that the amount of all gases sums up to the amount of dray air:
                ! 
                ! mass(dry air) = sum[ mass(jgas) ] = sum\[ mass(jgas) ] + mass(igasDiag) <=>
                !
                ! mass(igasDiag) = mass(dry air) - sum\[ mass(jgas) ],
                !
                ! where sum\ denotes the sum over all gases except for the diagnostic gas igasDiag.
                prm_upatmo_diag%gas(jc,jk,jb,igasDiag) = prm_upatmo_diag%gas(jc,jk,jb,igasDiag) &
                  &                                    - prm_upatmo_diag%gas(jc,jk,jb,jgas)                
              ENDDO  !jc
            ENDDO  !jk

          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

        ENDIF  !IF (jgas /= igasDiag)

      ENDDO  !jgas
      
      !--------------------------------
      ! Sanity check of diagnostic gas
      !--------------------------------
      
      ! 'lsanitycheck' can be set in 'src/namelists/mo_upatmo_nml'
      IF (upatmo_phy_config%lsanitycheck) THEN
        
        ! Is mass of diagnostic gas positive-definite?
        CALL sanity_check ( p_patch     = p_patch,                             &  !in
          &                 state       = prm_upatmo_diag%gas(:,:,:,igasDiag), &  !in
          &                 bound       = 0._wp,                               &  !in
          &                 keys        = 'cell,lower',                        &  !in
          &                 lpassed     = lpassed,                             &  !out
          &                 opt_slev    = 1,                                   &  !optin
          &                 opt_elev    = nlev,                                &  !optin
          &                 opt_message = sanity_check_msg                     )  !optout

        IF (.NOT. lpassed) THEN
          message_text = 'Mass of diagnostic gas N2 is not positive-definite: '//TRIM(sanity_check_msg) &
            & //' Please reconsider your settings for vmr and fscale in upatmo_nml '                    &
            & //'and/or check the external gas data.'
          CALL finish(TRIM(routine), TRIM(message_text))
        ENDIF
        
      ENDIF  !IF (upatmo_phy_config%lsanitycheck)
      
      !---------------------------------------------------------------------
      !                     Derived gas quantities
      !--------------------------------------------------------------------- 
      
      ! If this subroutine is called for initialization, 
      ! 'lupdate_gas = .TRUE.' may hold, but 'opt_gas_vmr', 
      ! 'opt_gas_col' and 'opt_gas_cumcol' are not present
      IF (.NOT. linit) THEN
        
        ! The following computations 
        ! have to be done for all gases
        DO jgas = 1, ngas
          
          ! Conversion factors
          mmr2vmr  = upatmo_nwp%gas( jgas )%mmr2vmr
          mass2mol = upatmo_nwp%gas( jgas )%mass2mol / 10._wp
          
!$OMP PARALLEL
!$OMP DO PRIVATE(jb, jk, jc, i_startidx, i_endidx) ICON_OMP_GUIDED_SCHEDULE
          DO jb = i_startblk, i_endblk
            
            CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,     &
              &                 i_startidx, i_endidx, rl_start, rl_end )
            
            ! Depending on the gas mode, some of the following computations 
            ! could be done computationally more efficiently, 
            ! but to simplify matters we do without a case differentiation for the time being.
            
            !--------------------------------------------------
            ! Volume mixing ratios and column number densities
            !--------------------------------------------------
            
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx 
                ! Volume mixing ratio ([kg/kg] -> [m3/m3] (= [mol/mol]))
                opt_gas_vmr(jc,jk,jb,jgas) = mmr2vmr * prm_upatmo_diag%gas(jc,jk,jb,jgas)
                ! Column number density ([kg/m2] -> [molecules/cm2])
                opt_gas_col(jc,jk,jb,jgas) = mass2mol * prm_upatmo_diag%gas(jc,jk,jb,jgas) * &
                  &                          prm_upatmo_diag%mdry(jc,jk,jb)              
              ENDDO  !jc
            ENDDO  !jk  
            
            !-------------------------------------
            ! Accumulated column number densities
            !-------------------------------------
            
            ! Initialization
            DO jc = i_startidx, i_endidx 
              opt_gas_cumcol(jc,1,jb,jgas) = opt_gas_col(jc,1,jb,jgas)
            ENDDO
            
            ! Integrate from model top downwards
            DO jk = 2, nlev
              DO jc = i_startidx, i_endidx  
                opt_gas_cumcol(jc,jk,jb,jgas) = opt_gas_cumcol(jc,jk-1,jb,jgas) &
                  &                           + opt_gas_col(jc,jk,jb,jgas)
              ENDDO  !jc
            ENDDO  !jk  
            
          ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL

        ENDDO  !jgas
        
      ENDIF  !IF (.NOT. linit)
      
      upatmo_nwp_gas => NULL()
      
    ENDIF  !IF (lupdate_gas)
    
    NULLIFY(upatmo_nwp, qv, qc, qi)

    IF (lmessage .AND. linit) THEN
      CALL message(TRIM(routine), &
        & 'Finish initialization of diagnostic variables for upper-atmosphere physics on domain '//dom_str)
    ELSEIF (lmessage) THEN 
      CALL message(TRIM(routine), &
        & 'Finish update of diagnostic variables for upper-atmosphere physics on domain '//dom_str)
    ENDIF

    IF (ltimer) CALL timer_stop(timer_upatmo_phy_diag)

  END SUBROUTINE update_diagnostic_variables

END MODULE mo_upatmo_phy_diag
