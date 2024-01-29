!
! Construction and destruction of the external data
! for the upper atmosphere.
!
!----------------------------
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_upatmo_extdat_state

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_impl_constants,         ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_upatmo_impl_const,      ONLY: iUpatmoStat, iUpatmoGasId, iUpatmoGasMode, &
    &                                  iUpatmoExtdatId, iUpatmoExtdatLatId,       &
    &                                  iUpatmoExtdatLevId, iUpatmoExtdatTimeId,   &
    &                                  iUpatmoExtdatStat
  USE mo_model_domain,           ONLY: t_patch
  USE mo_upatmo_config,          ONLY: t_upatmo_config
  USE mo_upatmo_phy_config,      ONLY: t_upatmo_phy_config
  USE mo_upatmo_types,           ONLY: t_upatmo_extdat
  USE mo_util_string,            ONLY: toupper, int2string
  USE mo_upatmo_extdat_utils,    ONLY: read_extdat_gas, read_extdat_chemheat, &
    &                                  construct_interpolation_lat,           &
    &                                  construct_interpolation_lev
  USE mo_upatmo_phy_chemheat,    ONLY: chem_heat_clean

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_upatmo_extdat_nwp
  PUBLIC :: destruct_upatmo_extdat_nwp

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_extdat_state'

CONTAINS

  !>
  !! Initialize external data for the upper atmosphere 
  !! under NWP forcing.
  !!
  SUBROUTINE construct_upatmo_extdat_nwp( jg,                &  !in
    &                                     nproma,            &  !in
    &                                     p_patch,           &  !in
    &                                     prm_upatmo_extdat, &  !inout
    &                                     upatmo_config,     &  !in
    &                                     upatmo_phy_config, &  !in
    &                                     vct_a              )  !(opt)in

    ! In/out variables
    INTEGER,                            INTENT(IN)    :: jg                 ! Domain index
    INTEGER,                            INTENT(IN)    :: nproma             ! For horizontal blocks
    TYPE(t_patch),             TARGET,  INTENT(IN)    :: p_patch            ! Grid/patch info
    TYPE(t_upatmo_extdat),              INTENT(INOUT) :: prm_upatmo_extdat  ! External data to be constructed
    TYPE(t_upatmo_config),              INTENT(IN)    :: upatmo_config      ! General upper-atmosphere configuration
    TYPE(t_upatmo_phy_config),          INTENT(IN)    :: upatmo_phy_config  ! Upper-atmosphere physics configuration 
                                                                            ! with namelist settings
    REAL(wp),                 OPTIONAL, INTENT(IN)    :: vct_a(:)           ! Nominal heights of grid layer interfaces

    ! Local variables
    REAL(wp) :: vmr2mmr

    INTEGER  :: mapgasid2indx( iUpatmoGasId%nitem )
    INTEGER  :: mapgasindx2id( iUpatmoGasId%nitem )
    INTEGER  :: jgas, ngas, nblks_c, nlev, nplev
    INTEGER  :: istat

    LOGICAL  :: lmessage

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: dom_str, filename, gasname, error_str

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':construct_upatmo_extdat_nwp'
    
    !-------------------------------------------------------------- 

    ! Some info:
    !
    ! * To store the external data domain-wise  
    !   means to store a lot of redundant information, 
    !   since there is only one globally valid external data set 
    !   for the time being. However, for the handling of external data  
    !   within ICON, it makes life a bit easier. 
    !   In addition, the domain-wise treatment means a lot of unnecessary 
    !   procedural overhead in the following, but since this setup 
    !   is only done once, we assume this to be bearable.

    !---------------------------------------------------------------------
    !                               Checks
    !---------------------------------------------------------------------

    ! Message output desired?
    lmessage  = upatmo_config%l_status( iUpatmoStat%message )
    dom_str   = TRIM(int2string(jg))
    error_str = 'Allocation of prm_upatmo_extdat('//TRIM(dom_str)//')'

    IF (prm_upatmo_extdat%linitialized) THEN
      CALL finish(TRIM(routine), 'External data are already initialized for domain '//TRIM(dom_str))
    ELSEIF (.NOT. PRESENT(vct_a)) THEN 
      CALL finish(TRIM(routine), 'vct_a has to be present.')
    ENDIF
    
    IF (lmessage) CALL message(TRIM(routine), &
      & 'Start construction of external data on domain '//TRIM(dom_str))
    
    !---------------------------------------------------------------------
    !                            Preparation
    !---------------------------------------------------------------------
    
    ! Number of levels
    nlev = p_patch%nlev
    
    ! Number of blocks
    nblks_c = p_patch%nblks_c

    !---------------------------------------------------------------------
    !                       Initialize external data:
    !                       Radiatively active gases
    !---------------------------------------------------------------------
    
    IF (upatmo_config%nwp_phy%extdat( iUpatmoExtdatId%gases )%l_stat( iUpatmoExtdatStat%required ) ) THEN
      
      ! Scan gases for mode 'extdat'
      ngas               = 0
      mapgasid2indx( : ) = 0
      mapgasindx2id( : ) = 0
      ! Loop over all radiatively active gases
      DO jgas = 1, iUpatmoGasId%nitem
        IF (upatmo_phy_config%nwp_gas( jgas )%imode == iUpatmoGasMode%extdat) THEN
          ngas                  = ngas + 1
          mapgasindx2id( ngas ) = jgas
          mapgasid2indx( jgas ) = ngas
        ENDIF
      ENDDO  !jgas
      
      prm_upatmo_extdat%ngas = ngas
      
      ! We can skip the following, if we do not have to read gases 
      ! from the external data file
      IF (ngas > 0) THEN
        
        ALLOCATE( prm_upatmo_extdat%mapgasindx2id( ngas ),               &
          &       prm_upatmo_extdat%mapgasid2indx( iUpatmoGasId%nitem ), &
          &       prm_upatmo_extdat%gas( ngas ),                         &
          &       prm_upatmo_extdat%gas_interm( ngas ),                  &
          &       STAT=istat                                             )
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%gas/%mapgasindx2id/%mapgasid2indx/gas_interm failed.')
        
        prm_upatmo_extdat%mapgasindx2id( 1:ngas )               = mapgasindx2id( 1:ngas )
        prm_upatmo_extdat%mapgasid2indx( 1:iUpatmoGasId%nitem ) = mapgasid2indx( 1:iUpatmoGasId%nitem )
        
        ! Name of file with external gas data
        filename = TRIM(upatmo_phy_config%nwp_extdat( iUpatmoExtdatId%gases )%filename)
        
        ! Loop over extdat gases
        DO jgas = 1, ngas
          
          prm_upatmo_extdat%gas( jgas )%data_id = prm_upatmo_extdat%mapgasindx2id( jgas )
          ! Latitudes are measured in degree north
          prm_upatmo_extdat%gas( jgas )%lat_id = iUpatmoExtdatLatId%deg
          ! Levels are pressure levels
          prm_upatmo_extdat%gas( jgas )%lev_id = iUpatmoExtdatLevId%p
          ! Times are months
          prm_upatmo_extdat%gas( jgas )%time_id = iUpatmoExtdatTimeId%month
          
          ! Gas name
          gasname = toupper(upatmo_config%nwp_phy%gas( jgas )%name)

          ! The external gas data are provided in units 
          ! of volume mixing ratio (mole fraction) [mol mol-1], 
          ! but we need them in units of mass mixing ration [kg kg-1]. 
          ! In addition, a constant scaling of the gas concentration 
          ! might be desired. We can incorporate the corresponding 
          ! scaling factor from the namelist here.
          vmr2mmr = upatmo_config%nwp_phy%gas( jgas )%vmr2mmr * &
            &       upatmo_phy_config%nwp_gas( jgas )%fscale
          
          ! For reasons of simplicity and convenience, 
          ! the read-in of the external data is done gas-wise
          ! (i.e. open and close the file, check dimensions etc.). 
          ! This is not very efficient, since all gases are stored 
          ! in one file for the time being. 
          ! But since this is done only once during model setup, 
          ! it should be bearable.
          CALL read_extdat_gas( gas          = prm_upatmo_extdat%gas( jgas ), &  !inout
            &                   gasname      = gasname,                       &  !in
            &                   vmr2mmr      = vmr2mmr,                       &  !in
            &                   filename     = filename,                      &  !in
            &                   opt_lmessage = lmessage                       )  !optin
          
          ALLOCATE( prm_upatmo_extdat%gas( jgas )%intrpl%lat%idx( 2, nproma, nblks_c ), &
            &       prm_upatmo_extdat%gas( jgas )%intrpl%lat%wgt( 2, nproma, nblks_c ), &
            &       STAT=istat                                                          )
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
            & TRIM(error_str)//'%gas('//TRIM(gasname)//')%intrpl%lat%idx/wgt failed.')
          
          ! Get auxiliary quantities for meridional interpolation during runtime
          CALL construct_interpolation_lat( p_patch    = p_patch,                                      &  !in
            &                               lat_stzstl = prm_upatmo_extdat%gas( jgas )%lat,            &  !in
            &                               istart     = prm_upatmo_extdat%gas( jgas )%istartlat,      &  !in
            &                               iend       = prm_upatmo_extdat%gas( jgas )%iendlat,        &  !in
            &                               istep      = prm_upatmo_extdat%gas( jgas )%isteplat,       &  !in
            &                               intrpl_idx = prm_upatmo_extdat%gas( jgas )%intrpl%lat%idx, &  !out
            &                               intrpl_wgt = prm_upatmo_extdat%gas( jgas )%intrpl%lat%wgt  )  !out
          
          ! In order to reduce the interpolation work during runtime, 
          ! we introduce a field, which contains the gas concentrations 
          ! already on the horizontal grid of ICON, 
          ! but still on the original pressure levels of the external data
          nplev = prm_upatmo_extdat%gas( jgas )%nlev
          ALLOCATE( prm_upatmo_extdat%gas_interm( jgas )%p( nproma, nplev, nblks_c ), &
            &       STAT=istat                                                         )
          IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
            & TRIM(error_str)//'%gas_interm('//TRIM(gasname)//')%p failed.')
          
        ENDDO  !jgas
        
      ENDIF  !IF (ngas > 0)
      
    ENDIF  !External gas data required?
    
    !---------------------------------------------------------------------
    !                      Initialize external data:
    !                     Chemical heating tendencies
    !---------------------------------------------------------------------
    
    IF (upatmo_config%nwp_phy%extdat( iUpatmoExtdatId%chemheat )%l_stat( iUpatmoExtdatStat%required ) ) THEN

      ! Name of file with chemical heating tendencies
      filename = TRIM(upatmo_phy_config%nwp_extdat( iUpatmoExtdatId%chemheat )%filename)
      
      ! Latitudes are measured in degree north
      prm_upatmo_extdat%chemheat%lat_id = iUpatmoExtdatLatId%deg
      ! Levels are geometric heights
      prm_upatmo_extdat%chemheat%lev_id = iUpatmoExtdatLevId%z
      ! Times are months
      prm_upatmo_extdat%chemheat%time_id = iUpatmoExtdatTimeId%month
      
      ! Read external data
      CALL read_extdat_chemheat( chemheat     = prm_upatmo_extdat%chemheat, &  !inout
        &                        filename     = filename,                   &  !in
        &                        opt_lmessage = lmessage                    )  !optin
      
      ALLOCATE( prm_upatmo_extdat%chemheat%intrpl%lat%idx( 2, nproma, nblks_c ), &
        &       prm_upatmo_extdat%chemheat%intrpl%lat%wgt( 2, nproma, nblks_c ), &
        &       prm_upatmo_extdat%chemheat%intrpl%lev%idx( 2, nlev ),            &
        &       prm_upatmo_extdat%chemheat%intrpl%lev%wgt( 2, nlev ),            &
        &       STAT=istat                                                       )
      IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
        & TRIM(error_str)//'%chemheat%intrpl%lat/lev%idx/wgt failed.')
      
      ! Get auxiliary quantities for meridional interpolation during runtime
      CALL construct_interpolation_lat( p_patch    = p_patch,                                   &  !in
        &                               lat_stzstl = prm_upatmo_extdat%chemheat%lat,            &  !in
        &                               istart     = prm_upatmo_extdat%chemheat%istartlat,      &  !in
        &                               iend       = prm_upatmo_extdat%chemheat%iendlat,        &  !in
        &                               istep      = prm_upatmo_extdat%chemheat%isteplat,       &  !in
        &                               intrpl_idx = prm_upatmo_extdat%chemheat%intrpl%lat%idx, &  !out
        &                               intrpl_wgt = prm_upatmo_extdat%chemheat%intrpl%lat%wgt  )  !out
      
      ! Get auxiliary quantities for vertical interpolation during runtime
      ! Please note that the chemical heating tendencies 
      ! are interpolated with respect to the nominal grid layer heights (stored in 'vct_a'), 
      ! not with respect to the actual height of each cell (stored in 'p_nh_state%metrics%z_mc') 
      ! for the following reasons:
      ! * Chemical heating tendencies are assumed 
      !   to become significant only far above 'flat_height'.
      ! * The external chemical heating tendencies
      !   provide only zonal averages, 
      !   so any terrain-imprint on the data 
      !   is significantly degraded.
      ! * A storage of the auxiliary quantities, 
      !   computed in the following, on the full 3d-grid 
      !   means a considerable consumption of memory.
      CALL construct_interpolation_lev( p_patch    = p_patch,                                   &  !in
        &                               lev_stzstl = prm_upatmo_extdat%chemheat%lev,            &  !in
        &                               istart     = prm_upatmo_extdat%chemheat%istartlev,      &  !in
        &                               iend       = prm_upatmo_extdat%chemheat%iendlev,        &  !in
        &                               istep      = prm_upatmo_extdat%chemheat%isteplev,       &  !in
        &                               intrpl_idx = prm_upatmo_extdat%chemheat%intrpl%lev%idx, &  !out
        &                               intrpl_wgt = prm_upatmo_extdat%chemheat%intrpl%lev%wgt, &  !out
        &                               vct_a      = vct_a                                      )  !(opt)in
      
    ENDIF  !External chemical heating tendencies required?

    IF (lmessage) CALL message(TRIM(routine), &
      & 'Finish construction of external data on domain '//TRIM(dom_str))
    
  END SUBROUTINE construct_upatmo_extdat_nwp

  !==================================================================================== 

  !>
  !! Destruct external data for the upper atmosphere 
  !! under NWP forcing.
  !!
  SUBROUTINE destruct_upatmo_extdat_nwp( jg,                &  !in
    &                                    prm_upatmo_extdat, &  !inout
    &                                    upatmo_config      )  !in

    ! In/out variables
    INTEGER,               INTENT(IN)    :: jg                 ! Domain index
    TYPE(t_upatmo_extdat), INTENT(INOUT) :: prm_upatmo_extdat  ! External data to be constructed
    TYPE(t_upatmo_config), INTENT(IN)    :: upatmo_config      ! General upper-atmosphere configuration

    ! Local variables
    INTEGER  :: istat, jgas

    LOGICAL  :: lmessage

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: dom_str, gasname, error_str

    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':destruct_upatmo_extdat_nwp'
    
    !-------------------------------------------------------------- 

    !---------------------------
    ! Destruct chemical heating 
    !---------------------------

    ! This has to be done only once for all domains. 
    ! There is a switch in 'src/upper_atmosphere/mo_upatmo_phy_chemheat' 
    ! to accomplish this. Yet still calling 'chem_heat_clean'
    ! for each domain is not very efficient, but the control switch 
    ! is private to the chemical heating module and we prefer to keep it this way.
    CALL chem_heat_clean()
    
    ! Message output desired?
    lmessage  = upatmo_config%l_status( iUpatmoStat%message )
    dom_str   = TRIM(int2string(jg))
    error_str = 'Deallocation of prm_upatmo_extdat('//TRIM(dom_str)//')'
    
    IF (lmessage) CALL message(TRIM(routine), &
      & 'Start destruction of external data on domain '//TRIM(dom_str))
    
    !---------------------------------------------------------------------
    !                      prm_upatmo_extdat%gas
    !---------------------------------------------------------------------
    
    IF (upatmo_config%nwp_phy%extdat( iUpatmoExtdatId%gases )%l_stat( iUpatmoExtdatStat%required ) ) THEN
      
      IF (ALLOCATED(prm_upatmo_extdat%gas)) THEN
        DO jgas = 1, prm_upatmo_extdat%ngas
          gasname = toupper(upatmo_config%nwp_phy%gas( jgas )%name)
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%data)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%data, STAT=istat)
            IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%data failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%lat)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%lat, STAT=istat)
            IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%lat failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%lev)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%lev, STAT=istat)
            IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%lev failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%lev_half)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%lev_half, STAT=istat)
            IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%lev_half failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%time)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%time, STAT=istat)
            IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%time failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%intrpl%lat%idx)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%intrpl%lat%idx, STAT=istat)
            IF(istat /= SUCCESS) & 
              & CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%intrpl%lat%idx failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%intrpl%lat%wgt)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%intrpl%lat%wgt, STAT=istat)
            IF(istat /= SUCCESS) &
              & CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%intrpl%lat%wgt failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%intrpl%lev%idx)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%intrpl%lev%idx, STAT=istat)
            IF(istat /= SUCCESS) & 
              & CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%intrpl%lev%idx failed.')
          ENDIF
          IF (ALLOCATED(prm_upatmo_extdat%gas( jgas )%intrpl%lev%wgt)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas( jgas )%intrpl%lev%wgt, STAT=istat)
            IF(istat /= SUCCESS) &
              & CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas('//TRIM(gasname)//')%intrpl%lev%wgt failed.')
          ENDIF
        ENDDO  !jgas
        DEALLOCATE(prm_upatmo_extdat%gas, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%gas failed.')
      ENDIF
      
      !---------------------------------------------------------------------
      !            rm_upatmo_extdat%mapgasid2indx/mapgasindx2id
      !---------------------------------------------------------------------
      
      IF (ALLOCATED(prm_upatmo_extdat%mapgasid2indx)) THEN
        DEALLOCATE(prm_upatmo_extdat%mapgasid2indx, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%mapgasid2indx failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%mapgasindx2id)) THEN
        DEALLOCATE(prm_upatmo_extdat%mapgasindx2id, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%mapgasindx2id failed.')
      ENDIF
      
      !---------------------------------------------------------------------
      !                    prm_upatmo_extdat%gas_interm
      !---------------------------------------------------------------------
      
      IF (ALLOCATED(prm_upatmo_extdat%gas_interm)) THEN
        DO jgas = 1, prm_upatmo_extdat%ngas
          gasname = toupper(upatmo_config%nwp_phy%gas( jgas )%name)
          IF (ASSOCIATED(prm_upatmo_extdat%gas_interm( jgas )%p)) THEN
            DEALLOCATE(prm_upatmo_extdat%gas_interm( jgas )%p, STAT=istat)
            IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
              & TRIM(error_str)//'%gas_interm('//TRIM(gasname)//')%p failed.')
          ENDIF
        ENDDO  !jgas
        DEALLOCATE(prm_upatmo_extdat%gas_interm, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%gas_interm failed.')
      ENDIF
      
    ENDIF  !External gas data required?
      
    !---------------------------------------------------------------------
    !                     prm_upatmo_extdat%chemheat
    !---------------------------------------------------------------------
    
    IF (upatmo_config%nwp_phy%extdat( iUpatmoExtdatId%chemheat )%l_stat( iUpatmoExtdatStat%required ) ) THEN
      
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%data)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%data, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%data failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%lat)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%lat, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%lat failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%lev)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%lev, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%lev failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%lev_half)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%lev_half, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%lev_half failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%time)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%time, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%time failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%intrpl%lat%idx)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%intrpl%lat%idx, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%intrpl%lat%idx failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%intrpl%lat%wgt)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%intrpl%lat%wgt, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%intrpl%lat%wgt failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%intrpl%lev%idx)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%intrpl%lev%idx, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%intrpl%lev%idx failed.')
      ENDIF
      IF (ALLOCATED(prm_upatmo_extdat%chemheat%intrpl%lev%wgt)) THEN
        DEALLOCATE(prm_upatmo_extdat%chemheat%intrpl%lev%wgt, STAT=istat)
        IF(istat /= SUCCESS) CALL finish(TRIM(routine), &
          & TRIM(error_str)//'%chemheat%intrpl%lev%wgt failed.')
      ENDIF
      
    ENDIF  !External chemical heating tendencies required?
    
    IF (lmessage) CALL message(TRIM(routine), &
      & 'Finish destruction of external data on domain '//TRIM(dom_str))
    
  END SUBROUTINE destruct_upatmo_extdat_nwp

END MODULE mo_upatmo_extdat_state

