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

! Initialization and clean-up of fields for upper-atmosphere physics.

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_upatmo_phy_setup

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_impl_constants,       ONLY: min_rlcell
  USE mo_physical_constants,   ONLY: amd, grav
  USE mo_upatmo_impl_const,    ONLY: iUpatmoStat, iUpatmoPrcStat, iUpatmoGasStat, &
    &                                iUpatmoGrpId,  iUpatmoExtdatStat
  USE mo_master_config,        ONLY: isRestart
  USE mo_upatmo_config,        ONLY: upatmo_config, upatmo_phy_config
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_metrics, t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag, t_nwp_phy_tend
  USE mo_upatmo_types,         ONLY: t_upatmo
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_upatmo_state,         ONLY: prm_upatmo_vec => prm_upatmo
  USE mo_upatmo_phy_diag,      ONLY: update_diagnostic_variables
  USE mo_upatmo_phy_chemheat,  ONLY: zeroz, onez, effrswmin
  USE mo_upatmo_phy_nlte,      ONLY: nlte_std_co2, nlte_set_co2_pre
  USE mtime,                   ONLY: datetime
  USE mo_timer,                ONLY: timer_start, timer_stop,   & 
    &                                timer_upatmo_phy, timer_upatmo_phy_init
  USE mo_util_string,          ONLY: int2string

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: init_upatmo_phy_nwp
  PUBLIC :: finalize_upatmo_phy_nwp

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_phy_setup'
  
CONTAINS

  !>
  !! Initialize upper-atmosphere physics state variables for NWP forcing.
  !!  
  SUBROUTINE init_upatmo_phy_nwp( mtime_datetime, &  !in
    &                             p_patch,        &  !in
    &                             p_metrics,      &  !in
    &                             p_prog,         &  !in
    &                             p_diag,         &  !in
    &                             prm_nwp_diag,   &  !in
    &                             prm_nwp_tend,   &  !in
    &                             nproma          )  !in

    ! In/out-variables
    TYPE(datetime),       POINTER,INTENT(IN)    :: mtime_datetime
    TYPE(t_patch),        TARGET, INTENT(IN)    :: p_patch 
    TYPE(t_nh_metrics),           INTENT(IN)    :: p_metrics
    TYPE(t_nh_prog),              INTENT(IN)    :: p_prog
    TYPE(t_nh_diag),              INTENT(IN)    :: p_diag
    TYPE(t_nwp_phy_diag),         INTENT(IN)    :: prm_nwp_diag
    TYPE(t_nwp_phy_tend),         INTENT(IN)    :: prm_nwp_tend
    INTEGER,                      INTENT(IN)    :: nproma

    ! Local variables
    TYPE(t_upatmo),       POINTER               :: prm_upatmo  ! WS: convenience pointer
    REAL(wp) :: scale4effrsw, fac4effrsw
    INTEGER  :: jg, jb, jk, jc, jgrp
    INTEGER  :: nlev
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk 
    INTEGER  :: i_startidx, i_endidx 
    LOGICAL  :: lmessage, ltimer, lupdate_gas, lchemheat, lrestart

    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':init_upatmo_phy_nwp'

    !---------------------------------------------------------

    ! Domain index
    jg  = p_patch%id
    prm_upatmo => prm_upatmo_vec(jg)

    ltimer = upatmo_config(jg)%l_status( iUpatmoStat%timer )

    IF (ltimer) THEN 
      CALL timer_start(timer_upatmo_phy)
      CALL timer_start(timer_upatmo_phy_init)
    ENDIF

    IF (.NOT. upatmo_config(jg)%l_status( iUpatmoStat%configured )) THEN 
      CALL finish(routine, 'Upper atmosphere not yet configured.')
    ELSEIF (.NOT. prm_upatmo%diag%linitialized) THEN
      CALL finish(routine, 'prm_upatmo%diag not yet initialized.')
    ELSEIF (.NOT. prm_upatmo%tend%linitialized) THEN
      CALL finish(routine, 'prm_upatmo%tend not yet initialized.')
    ELSEIF (.NOT. prm_upatmo%extdat%linitialized) THEN
      CALL finish(routine, 'prm_upatmo%extdat not yet initialized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%initialized )) THEN
      CALL finish(routine, 'Gases already initialized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%initialized )) THEN
      CALL finish(routine, 'Parameterizations already initialized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%initialized )) THEN
      CALL finish(routine, 'External data already initialized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%finalized )) THEN
      CALL finish(routine, 'Gases already finalized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%finalized )) THEN
      CALL finish(routine, 'Parameterizations already finalized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%finalized )) THEN
      CALL finish(routine, 'External data already finalized.')
    ! Check allocation status of some (but not all!) fields 
    ! that we need for the upper-atmosphere physics
    ELSEIF (.NOT. ASSOCIATED(p_prog%tracer)) THEN
      CALL finish(routine, 'p_prog%tracer is not allocated.')
    ELSEIF (.NOT. ASSOCIATED(p_diag%ddt_vn_phy)) THEN
      CALL finish(routine, 'p_diag%ddt_vn_phy is not allocated.')
    ELSEIF (.NOT. ASSOCIATED(p_diag%ddt_exner_phy)) THEN
      CALL finish(routine, 'p_diag%ddt_exner_phy is not allocated.')
    ELSEIF (.NOT. ASSOCIATED(prm_nwp_diag%cosmu0)) THEN
      CALL finish(routine, 'prm_nwp_diag%cosmu0 is not allocated.')
    ELSEIF (.NOT. ASSOCIATED(prm_nwp_tend%ddt_temp_radsw)) THEN
      CALL finish(routine, 'prm_nwp_tend%ddt_temp_radsw is not allocated.')
    ELSEIF (.NOT. ASSOCIATED(prm_nwp_tend%ddt_temp_radlw)) THEN
      CALL finish(routine, 'prm_nwp_tend%ddt_temp_radlw is not allocated.')
    ENDIF

    lmessage = upatmo_config(jg)%l_status( iUpatmoStat%message )

    IF (lmessage) THEN
      WRITE (message_text, '(a,i0)') 'Initialization of upper-atmosphere &
        &physics for NWP forcing started on domain ', jg
      CALL message(routine, message_text)
    END IF

    lrestart = isRestart()

    !---------------------------------------------------------------------
    !                  Initialize diagnostic variables
    !---------------------------------------------------------------------

    ! Any upper-atmosphere physics switched on?
    ! All fields in 'prm_upatmo' (except for 'prm_upatmo%extdat') 
    ! are stored in the restart file. 
    ! So the following should not be done in case of a restart.
    IF (upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled )) THEN

      !--------------------------
      ! Time-dependent variables
      !--------------------------

      lupdate_gas = upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%enabled )

      CALL update_diagnostic_variables( mtime_datetime    = mtime_datetime,        & !in
        &                               lupdate_gas       = lupdate_gas,           & !in
        &                               p_patch           = p_patch,               & !in
        &                               p_prog_rcf        = p_prog,                & !in
        &                               p_diag            = p_diag,                & !in
        &                               prm_upatmo_diag   = prm_upatmo%diag,       & !inout
        &                               prm_upatmo_tend   = prm_upatmo%tend,       & !inout
        &                               prm_upatmo_extdat = prm_upatmo%extdat,     & !inout
        &                               upatmo_config     = upatmo_config(jg),     & !in
        &                               upatmo_phy_config = upatmo_phy_config(jg), & !in
        &                               nproma            = nproma,                & !in
        &                               opt_linit         = .TRUE.,                & !optin
        &                               opt_lrestart      = lrestart               ) !optin

      !----------------------------
      ! Time-independent variables
      !----------------------------

      ! The fields contained by 'prm_upatmo%diag' would be stored 
      ! in the restart file as the case may be. 
      ! So we can skip most of the following computations 
      ! in case of a resumed simulation.
      IF (.NOT. lrestart) THEN
        
        ! Number of vertical grid levels
        nlev = p_patch%nlev
        
        ! Some switches 
        lchemheat    = upatmo_config(jg)%nwp_phy%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%enabled ) .AND. &
          &            .NOT. upatmo_config(jg)%nwp_phy%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%offline )
        
        ! Auxiliary factors for efficiency factor from chemical heating
        IF (lchemheat) THEN
          scale4effrsw = 1._wp / (onez - zeroz)
          fac4effrsw   = 1._wp - effrswmin
        ELSE
          scale4effrsw = 0._wp
          fac4effrsw   = 0._wp
        ENDIF
        
        ! Do initialization for all cells
        rl_start   = 1
        rl_end     = min_rlcell
        i_startblk = p_patch%cells%start_block(rl_start)
        i_endblk   = p_patch%cells%end_block(rl_end)
        
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx) ICON_OMP_GUIDED_SCHEDULE
        DO jb = i_startblk, i_endblk
          
          CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
          
          DO jk = 1, nlev
            DO jc = i_startidx, i_endidx  
              
              ! Gravitational acceleration
              ! (includes deep-atmosphere modification)
              prm_upatmo%diag%grav(jc,jk,jb) = grav * p_metrics%deepatmo_gradh_mc(jk)**2
              
              ! Molar mass of dry air 
              ! (constant for the time being)
              prm_upatmo%diag%amd(jc,jk,jb) = amd
              
            ENDDO  !jc
          ENDDO  !jk
          
          IF (upatmo_config(jg)%nwp_phy%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%enabled )) THEN
            
            ! The scaling and efficiency factors are actually only required, 
            ! if the radiation group is switched on.    
            
            DO jk = 1, nlev
              DO jc = i_startidx, i_endidx  
                
                ! Scaling factor for heating rate from "standard" long-wave radiation
                ! (just and initialization, it will be diagnosed in the nwp-upatmo-interface)
                prm_upatmo%diag%sclrlw(jc,jk,jb) = 1._wp
                
                ! Time-independent efficiency factor for heating rate 
                ! from "standard" short-wave radiation
                prm_upatmo%diag%effrsw(jc,jk,jb) = 1._wp - fac4effrsw * &
                  & MIN(1._wp, MAX(0._wp, scale4effrsw * (p_metrics%z_mc(jc,jk,jb) - zeroz)))
                
              ENDDO  !jc
            ENDDO  !jk
            
          ENDIF  !Radiation enabled?
          
        ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL
        
      ENDIF  !IF (.NOT. lrestart)

    ENDIF  !Any upper-atmosphere physics switched on?

    ! Indicate that initialization took place:
    ! * Gases
    upatmo_config(jg)%nwp_phy%gas( : )%l_stat( iUpatmoGasStat%initialized ) = &
      & upatmo_config(jg)%nwp_phy%gas( : )%l_stat( iUpatmoGasStat%enabled )
    upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%initialized ) = .TRUE.
    ! * External data
    upatmo_config(jg)%nwp_phy%extdat( : )%l_stat( iUpatmoExtdatStat%initialized ) = .TRUE.
    upatmo_config(jg)%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%initialized )      = .TRUE.

    !---------------------------------------------------------------------------------
    !      Initialize sub-packets of upper-atmosphere physics parameterizations
    !---------------------------------------------------------------------------------

    DO jgrp = 1, iUpatmoGrpId%nitem

      ! Physics group enabled?
      IF (upatmo_config(jg)%nwp_phy%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled )) THEN

        IF (jgrp == iUpatmoGrpId%imf) THEN
          
          !---------------------------------------------------------------------------------
          !   Initialize ion drag (I), molecular diffusion (M) and frictional heating (F)
          !---------------------------------------------------------------------------------
          
        ELSEIF (jgrp == iUpatmoGrpId%rad) THEN
          
          !---------------------------------------------------------------------------------
          !                   Initialize radiation and chemical heating 
          !---------------------------------------------------------------------------------

          ! Non-LTE infrared cooling due to CO2 and O3.
          ! This setup has to be done only once for all domains. 
          ! There are switches in 'src/upper_atmosphere/mo_upatmo_phy_nlte' to guarantee 
          ! a onetime setup. Although less efficient, we prefer to keep the check 
          ! for the switch position within the subroutines, in order to retain
          ! the safer private status of the switches.
          CALL nlte_std_co2()
          CALL nlte_set_co2_pre()
          
        ELSE
          
          CALL finish(routine, 'Please implement initialization of new physics group ' &
            & //'into upper_atmosphere/mo_upatmo_phy_setup: init_upatmo_phy_nwp. Thank you!'  ) 
          
        ENDIF  !IF (jgrp == iUpatmoGrpId%...)
        
      ENDIF  !Physics group enabled?
      
      ! Indicate initialization
      upatmo_config(jg)%nwp_phy%grp( jgrp )%l_stat( iUpatmoPrcStat%initialized ) = .TRUE.

    ENDDO  !jgrp

    ! Indicate initialization for single processes
    upatmo_config(jg)%nwp_phy%prc( : )%l_stat( iUpatmoPrcStat%initialized ) = .TRUE.

    ! Indicate that initialization took place
    upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%initialized ) = .TRUE.

    IF (lmessage) THEN
      WRITE (message_text, '(a,i0)') 'Initialization of upper-atmosphere &
        &physics for NWP forcing finished on domain ', jg
      CALL message(routine, message_text)
    END IF

    IF (ltimer) THEN
      CALL timer_stop(timer_upatmo_phy_init)
      CALL timer_stop(timer_upatmo_phy)
    ENDIF

  END SUBROUTINE init_upatmo_phy_nwp

  !====================================================================================

  !>
  !! Finalize upper-atmosphere physics for NWP forcing.
  !!
  SUBROUTINE finalize_upatmo_phy_nwp( p_patch )

    ! In/out-variables
    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch 

    ! Local variables
    INTEGER  :: jg, jgrp
    LOGICAL  :: lmessage

    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':finalize_upatmo_phy_nwp'

    !---------------------------------------------------------

    ! Domain index
    jg  = p_patch%id

    IF (.NOT. upatmo_config(jg)%l_status( iUpatmoStat%configured )) THEN
      CALL finish(routine, 'Error: upper atmosphere is not configured.')
    ELSEIF (.NOT. upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%initialized )) THEN
      CALL finish(routine, 'Error: gases not initialized.')
    ELSEIF (.NOT. upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%initialized )) THEN
      CALL finish(routine, 'Error: parameterizations not initialized.')
    ELSEIF (.NOT. upatmo_config(jg)%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%initialized )) THEN
      CALL finish(routine, 'Error: external data not initialized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%finalized )) THEN
      CALL finish(routine, 'Error: gases already finalized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%finalized )) THEN
      CALL finish(routine, 'Error: parameterizations already finalized.')
    ELSEIF (upatmo_config(jg)%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%finalized )) THEN
      CALL finish(routine, 'Error: external data already finalized.')
    ENDIF

    lmessage = upatmo_config(jg)%l_status( iUpatmoStat%message )

    IF (lmessage) THEN
      WRITE (message_text, '(a,i0)') 'Finalization of upper-atmosphere &
           &physics for NWP forcing started on domain ', jg
      CALL message(routine, message_text)
    END IF

    !---------------------------------------------------------------------------------
    !      Finalize sub-packets of upper-atmosphere physics parameterizations
    !---------------------------------------------------------------------------------

    DO jgrp = 1, iUpatmoGrpId%nitem

      ! Physics group enabled?
      IF (upatmo_config(jg)%nwp_phy%grp( jgrp )%l_stat( iUpatmoPrcStat%enabled )) THEN
        
        IF (jgrp == iUpatmoGrpId%imf) THEN
          
          !---------------------------------------------------------------------------------
          !    Finalize ion drag (I), molecular diffusion (M) and frictional heating (F)
          !---------------------------------------------------------------------------------
          
        ELSEIF (jgrp == iUpatmoGrpId%rad) THEN
          
          !---------------------------------------------------------------------------------
          !                     Finalize radiation and chemical heating 
          !---------------------------------------------------------------------------------

          ! The destruction of the auxiliary fields in 'src/upper_atmosphere/mo_upatmo_phy_chemheat'
          ! is triggered in 'src/upper_atmosphere/mo_upatmo_extdat_state: destruct_upatmo_extdat_nwp'.
          
        ELSE
          
          CALL finish(routine, 'Please implement finalization of new physics group '      &
            & //'into upper_atmosphere/mo_upatmo_phy_setup: finalize_upatmo_phy_nwp. Thank you!' ) 
          
        ENDIF  !IF (jgrp == iUpatmoGrpId%...)

      ENDIF  !Physics group enabled?
        
      ! Indicate finalization
      upatmo_config(jg)%nwp_phy%grp( jgrp )%l_stat( iUpatmoPrcStat%finalized ) = .TRUE.

    ENDDO  !jgrp

    ! Indicate finalization for single processes
    upatmo_config(jg)%nwp_phy%prc( : )%l_stat( iUpatmoPrcStat%finalized ) = .TRUE.

    !---------------------------------------------------------------------
    !                 Finalize radiatively active gases
    !---------------------------------------------------------------------

    ! Indicate finalization:
    ! * Gases
    upatmo_config(jg)%nwp_phy%gas( : )%l_stat( iUpatmoGasStat%finalized ) = .TRUE.
    upatmo_config(jg)%nwp_phy%l_gas_stat( iUpatmoGasStat%finalized )      = .TRUE.
    ! * External data
    upatmo_config(jg)%nwp_phy%extdat( : )%l_stat( iUpatmoExtdatStat%finalized ) = .TRUE.
    upatmo_config(jg)%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%finalized )      = .TRUE.

    ! Indicate that finalization took place
    upatmo_config(jg)%nwp_phy%l_phy_stat( iUpatmoPrcStat%finalized ) = .TRUE.

    IF (lmessage) THEN
      WRITE (message_text, '(a,i0)') 'Finalization of upper-atmosphere &
        &physics for NWP forcing finished on domain ', jg
      CALL message(routine, message_text)
    END IF

  END SUBROUTINE finalize_upatmo_phy_nwp

END MODULE mo_upatmo_phy_setup
