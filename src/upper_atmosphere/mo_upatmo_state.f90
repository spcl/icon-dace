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
! Conceptual copy of:
! * src/atm_dyn_iconam/mo_nonhydro_state
! * src/atm_phy_nwp/mo_nwp_phy_state
! * src/atm_phy_echam/mo_echam_phy_memory
! for the upper-atmosphere variables.

MODULE mo_upatmo_state

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH, VINTP_METHOD_LIN
  USE mo_upatmo_impl_const,    ONLY: iUpatmoStat, iUpatmoPrcStat, iUpatmoGasStat, &
    &                                iUpatmoTendId, iUpatmoGrpId, iUpatmoGasId,   &
    &                                iUpatmoTracerId, iUpatmoPrcId,               &
    &                                iUpatmoExtdatStat, itmr_thr
  USE mo_master_control,       ONLY: get_my_process_name
  USE mo_model_domain,         ONLY: t_patch
  USE mo_upatmo_types,         ONLY: t_upatmo_diag, t_upatmo_tend, t_upatmo
  USE mo_upatmo_config,        ONLY: t_upatmo_config
  USE mo_upatmo_phy_config,    ONLY: t_upatmo_nwp_phy, t_upatmo_phy_config
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var
  USE mo_io_config,            ONLY: lnetcdf_flt64_output
  USE mo_var_list,             ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_cdi,                  ONLY: DATATYPE_PACK16, DATATYPE_FLT32, DATATYPE_FLT64, & 
    &                                GRID_UNSTRUCTURED
  USE mo_cdi_constants,        ONLY: GRID_CELL, GRID_UNSTRUCTURED_CELL, &
    &                                GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_zaxis_type,           ONLY: ZA_REFERENCE
  USE mo_var_metadata,         ONLY: create_vert_interp_metadata, vintp_types
  USE mo_var_groups,           ONLY: groups
  USE mo_upatmo_extdat_state,  ONLY: construct_upatmo_extdat_nwp, &
    &                                destruct_upatmo_extdat_nwp
  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop,   &
    &                                timer_upatmo_constr, timer_upatmo_destr, &
    &                                timer_upatmo


  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: prm_upatmo
  PUBLIC :: construct_upatmo_state
  PUBLIC :: destruct_upatmo_state

  TYPE(t_upatmo),   ALLOCATABLE, TARGET :: prm_upatmo(:)           ! Shape: (n_dom)

  TYPE(t_var_list_ptr), ALLOCATABLE         :: prm_upatmo_diag_list(:) ! Shape: (n_dom)

  TYPE(t_var_list_ptr), ALLOCATABLE         :: prm_upatmo_tend_list(:) ! Shape: (n_dom)

  ! Please note that we cannot use 'mo_impl_constants: TIMELEVEL_SUFFIX' 
  ! for the different time levels of the total tendencies, 
  ! because of its use in 'mo_var_list'
  CHARACTER(LEN=*), PARAMETER :: STATELEVEL_SUFFIX = '.SL'

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_state'
  
CONTAINS

  !>
  !! Construct upper-atmosphere state variables.
  !!
  SUBROUTINE construct_upatmo_state( n_dom,             & !in
    &                                nproma,            & !in
    &                                p_patch,           & !in
    &                                upatmo_config,     & !in
    &                                upatmo_phy_config, & !in
    &                                vct_a              ) !(opt)in

    ! In/out variables
    INTEGER,                       INTENT(IN) :: n_dom
    INTEGER,                       INTENT(IN) :: nproma
    TYPE(t_patch),         TARGET, INTENT(IN) :: p_patch(:)
    TYPE(t_upatmo_config),         INTENT(IN) :: upatmo_config(:)
    TYPE(t_upatmo_phy_config),     INTENT(IN) :: upatmo_phy_config(:)
    REAL(wp),            OPTIONAL, INTENT(IN) :: vct_a(:)

    ! Local variables 
    INTEGER  :: jg, istat, nblks_c, nblks_e, nlev
    LOGICAL  :: lmessage, ltimer
    CHARACTER(LEN=80) :: listname
    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':construct_upatmo_state'

    !---------------------------------------------------------

    ltimer   = timers_level > itmr_thr%med
    IF (ltimer) THEN
      CALL timer_start(timer_upatmo)
      CALL timer_start(timer_upatmo_constr)
    ENDIF

    IF (.NOT. ANY(upatmo_config( : )%l_status( iUpatmoStat%configured ))) THEN
      CALL finish (routine, 'Information required is not yet available')
    ELSEIF (ALLOCATED(prm_upatmo)) THEN
      CALL finish (routine, 'prm_upatmo is already allocated')
    ELSEIF (.NOT. PRESENT(vct_a)) THEN 
      CALL finish(routine, 'vct_a has to be present')
    ENDIF

    lmessage = ANY(upatmo_config( : )%l_status( iUpatmoStat%message ))
    IF (lmessage) CALL message (routine, 'Construction of upatmo state started')

    !------------------------------------
    !     Upper-atmosphere physics
    !------------------------------------

    ! Unfortunately some compilers require the domain allocation 
    ! of 'prm_upatmo' to always take place, 
    ! so we had to move the query for the enabling of the physics 
    ! to after their allocation
      
    ! Cummulative upper-atmosphere type:
    ! * Diagnostic fields
    ! * Tendencies from upper-atmosphere physics parameterizations
    ! * External data
    ALLOCATE(prm_upatmo( n_dom ), STAT=istat)
    IF(istat/=SUCCESS) CALL finish (routine, 'Allocation of prm_upatmo failed')

    ! Any upper-atmosphere physics switched on?
    IF (ANY(upatmo_config( : )%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled ))) THEN

      ! List for diagnostic fields
      ALLOCATE(prm_upatmo_diag_list( n_dom ), STAT=istat)
      IF(istat/=SUCCESS) CALL finish (routine, 'Allocation of prm_upatmo_diag_list failed')
      
      ! List for tendencies
      ALLOCATE(prm_upatmo_tend_list( n_dom ), STAT=istat)
      IF(istat/=SUCCESS) CALL finish (routine, 'Allocation of prm_upatmo_tend_list failed')

      DO jg = 1, n_dom

        ! Physics switched on on domain?
        IF (upatmo_config( jg )%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled )) THEN

          ! Determine size of arrays
          nblks_c = p_patch( jg )%nblks_c
          nblks_e = p_patch( jg )%nblks_e
          
          ! Number of vertical levels
          nlev   = p_patch( jg )%nlev
          
          ! Prefix for variable names

          WRITE(listname,'(a,i2.2)') 'prm_upatmo_diag_of_domain_', jg

          ! Allocate diagnostic upper-atmosphere fields
          CALL new_upatmo_diag_list( jg, nlev, nblks_c, nproma, listname,               &
            &                        upatmo_config( jg )%nwp_phy%vname_prefix,          &
            &                        upatmo_config( jg )%nwp_phy,                       &
            &                        prm_upatmo_diag_list( jg ), prm_upatmo( jg )%diag  )
          
          WRITE(listname,'(a,i2.2)') 'prm_upatmo_tend_of_domain_', jg
          
          ! Allocate tendencies from upper-atmosphere physics parameterizations
          CALL new_upatmo_tend_list( jg, nlev, nblks_c, nblks_e, nproma, listname,      &
            &                        upatmo_config( jg )%nwp_phy%vname_prefix,          &
            &                        upatmo_config( jg )%nwp_phy,                       &
            &                        prm_upatmo_tend_list( jg ), prm_upatmo( jg )%tend, &
            &                        lmessage                                           )
          
        ENDIF  !Physics switched on on domain?

        ! Indicate that setup of 'prm_upatmo_diag' and 'prm_upatmo_tend' took place
        prm_upatmo( jg )%diag%linitialized = .TRUE. 
        prm_upatmo( jg )%tend%linitialized = .TRUE.

        ! The external upper-atmosphere data
        ! are not allocated with the add_var procedure
        IF (upatmo_config( jg )%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%required )) THEN
          CALL construct_upatmo_extdat_nwp( jg                = jg,                      & !in
            &                               nproma            = nproma,                  & !in
            &                               p_patch           = p_patch( jg ),           & !in
            &                               prm_upatmo_extdat = prm_upatmo( jg )%extdat, & !inout
            &                               upatmo_config     = upatmo_config( jg ),     & !in
            &                               upatmo_phy_config = upatmo_phy_config( jg ), & !in
            &                               vct_a             = vct_a                    ) !(opt)in
        ENDIF

        ! Indicate initialization
        ! (should be set to .TRUE. even if external date are not required)
        prm_upatmo( jg )%extdat%linitialized = .TRUE.

      ENDDO  !jg

    ENDIF !Any upper-atmosphere physics switched on?

    IF (lmessage) CALL message (routine, 'Upatmo state construction completed')

    IF (ltimer) THEN
      CALL timer_stop(timer_upatmo_constr)
      CALL timer_stop(timer_upatmo)
    ENDIF

  END SUBROUTINE construct_upatmo_state

  !==================================================================================== 

  !>
  !! Destruct upper atmosphere physics state variables.
  !!
  SUBROUTINE destruct_upatmo_state( n_dom,        & !in
    &                               upatmo_config ) !in

    ! In/out variables
    INTEGER,               INTENT(IN) :: n_dom
    TYPE(t_upatmo_config), INTENT(IN) :: upatmo_config(:)

    ! Local variables 
    INTEGER  :: jg, jst, istat
    INTEGER  :: nstate
    LOGICAL  :: lmessage, ltimer
    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':destruct_upatmo_state'

    !---------------------------------------------------------

    ltimer   = timers_level > itmr_thr%med
    IF (ltimer) THEN
      CALL timer_start(timer_upatmo)
      CALL timer_start(timer_upatmo_destr)
    ENDIF

    lmessage = ANY(upatmo_config( : )%l_status( iUpatmoStat%message ))

    IF (lmessage) CALL message (routine, 'Destruction of upatmo state started')

    !------------------------------------
    !     Upper-atmosphere physics
    !------------------------------------

    ! Any upper-atmosphere physics switched on?
    IF (ANY(upatmo_config( : )%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled ))) THEN

      DO jg = 1, n_dom

        ! Physics switched on on domain?
        IF (upatmo_config( jg )%nwp_phy%l_phy_stat( iUpatmoPrcStat%enabled )) THEN

          CALL vlr_del(prm_upatmo_diag_list(jg))
          
          CALL vlr_del(prm_upatmo_tend_list(jg))     

          ! Deallocate the fields, which have not been allocated via add_var/add_ref.
          ! Currently these are:
          ! * prm_upatmo%diag%gas_ptr
          ! * prm_upatmo%tend%ddt_qx_vdfmol_ptr
          ! * prm_upatmo%tend%ddt%temp
          ! * prm_upatmo%tend%ddt%exner
          ! * prm_upatmo%tend%ddt%vn
          ! * prm_upatmo%tend%ddt%qx%tot_ptr
          ! * prm_upatmo%tend%ddt%qx
          ! * prm_upatmo%tend%ddt%info
          ! * prm_upatmo%tend%ddt%state
  
          ! prm_upatmo%diag%gas_ptr
          IF (ALLOCATED(prm_upatmo( jg )%diag%gas_ptr)) THEN
            DEALLOCATE(prm_upatmo( jg )%diag%gas_ptr, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%diag%gas_ptr failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          ! prm_upatmo%tend%ddt_qx_vdfmol_ptr
          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt_qx_vdfmol_ptr)) THEN
            DEALLOCATE(prm_upatmo( jg )%tend%ddt_qx_vdfmol_ptr, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt_qx_vdfmol_ptr failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          ! prm_upatmo%tend%ddt%temp
          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%temp)) THEN
            DEALLOCATE(prm_upatmo( jg )%tend%ddt%temp, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt%temp failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          ! prm_upatmo%tend%ddt%exner
          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%exner)) THEN
            DEALLOCATE(prm_upatmo( jg )%tend%ddt%exner, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt%exner failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          ! prm_upatmo%tend%ddt%vn
          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%vn)) THEN
            DEALLOCATE(prm_upatmo( jg )%tend%ddt%vn, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt%vn failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%qx)) THEN
            ! prm_upatmo%tend%ddt%qx%tot_ptr
            nstate = prm_upatmo( jg )%tend%ddt%info( iUpatmoTendId%qx )%nstate
            DO jst = 1, nstate
              ! prm_upatmo%tend%ddt%qx%tot_ptr
              IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%qx( jst )%tot_ptr)) THEN
                DEALLOCATE(prm_upatmo( jg )%tend%ddt%qx( jst )%tot_ptr, STAT=istat)
                IF(istat/=SUCCESS) THEN
                  WRITE (message_text, '(a,i0,a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                       ')%tend%ddt%qx(', jst, ')%tot_ptr failed'
                  CALL finish (routine, message_text)
                ENDIF
              ENDIF
            ENDDO  !jst
            ! prm_upatmo%tend%ddt%qx
            DEALLOCATE(prm_upatmo( jg )%tend%ddt%qx, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt%qx failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          ! prm_upatmo%tend%ddt%info
          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%info)) THEN
            DEALLOCATE(prm_upatmo( jg )%tend%ddt%info, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt%info failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

          ! prm_upatmo%tend%ddt%state
          IF (ALLOCATED(prm_upatmo( jg )%tend%ddt%state)) THEN
            DEALLOCATE(prm_upatmo( jg )%tend%ddt%state, STAT=istat)
            IF(istat/=SUCCESS) THEN
              WRITE (message_text, '(a,i0,a)') 'Deallocation prm_upatmo(', jg, &
                   ')%tend%ddt%state failed'
              CALL finish(routine, message_text)
            ENDIF
          ENDIF

        ENDIF  !Physics switched on on domain?

        ! Destruct the external data
        IF (upatmo_config( jg )%nwp_phy%l_extdat_stat( iUpatmoExtdatStat%required )) THEN
          CALL destruct_upatmo_extdat_nwp( jg                = jg,                      &  !in
            &                              prm_upatmo_extdat = prm_upatmo( jg )%extdat, &  !inout
            &                              upatmo_config     = upatmo_config( jg )      )  !in
        ENDIF

      ENDDO  !jg

      DEALLOCATE(prm_upatmo_diag_list, STAT=istat)
      IF(istat/=SUCCESS) CALL finish (routine, 'Deallocation prm_upatmo_diag_list failed')
      
      DEALLOCATE(prm_upatmo_tend_list, STAT=istat)
      IF(istat/=SUCCESS) CALL finish (routine, 'Deallocation prm_upatmo_tend_list failed')

    ENDIF !Any upper-atmosphere physics switched on?

    DEALLOCATE(prm_upatmo, STAT=istat)
    IF(istat/=SUCCESS) CALL finish (routine, 'Deallocation prm_upatmo failed')

    IF (lmessage) CALL message (routine, 'Upatmo state destruction completed')

    IF (ltimer) THEN
      CALL timer_stop(timer_upatmo_destr)
      CALL timer_stop(timer_upatmo)
    ENDIF

  END SUBROUTINE destruct_upatmo_state

  !==================================================================================== 

  !>
  !! Allocation of diagnostic upper-atmosphere physics state variables.
  !!
  SUBROUTINE new_upatmo_diag_list( jg, nlev, nblks_c, nproma, listname, vname_prefix, &
    &                              upatmo_nwp_phy_config, diag_list, diag             )

    ! In/out variables
    INTEGER,                INTENT(IN)    :: jg, nlev, nblks_c, nproma
    CHARACTER(LEN=*),       INTENT(IN)    :: listname
    CHARACTER(LEN=*),       INTENT(IN)    :: vname_prefix
    TYPE(t_upatmo_nwp_phy), INTENT(IN)    :: upatmo_nwp_phy_config
    TYPE(t_var_list_ptr),       INTENT(INOUT) :: diag_list
    TYPE(t_upatmo_diag),    INTENT(INOUT) :: diag

    ! Local variables 
    INTEGER  :: shape3d(3), shape4d(4)
    INTEGER  :: ibits
    INTEGER  :: datatype_flt
    INTEGER  :: istat
    INTEGER  :: jgas
    INTEGER  :: vn_pfx_len

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    CHARACTER(len=LEN(vname_prefix)+10) :: var_name
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: var_name_ref

    INTEGER, PARAMETER :: ngas = iUpatmoGasId%nitem
    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':new_upatmo_diag_list'

    !---------------------------------------------------------

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 

    shape3d = (/nproma, nlev, nblks_c/)

    ! Ensure that all pointers have a defined association status
    ! (should be covered by the '=> NULL()'-statement at their 
    ! declaration, but to be on the safe side ...)
    NULLIFY( diag%gas,     &
      &      diag%mdry,    &
      &      diag%amd,     &
      &      diag%cpair,   &
      &      diag%grav,    &
      &      diag%sclrlw,  &
      &      diag%effrsw   )

    ! Register the field list and apply default settings
    ! The fields in 'prm_upatmo_diag' are required for restart, in general
    CALL vlr_add(diag_list, listname, patch_id=jg, lrestart=.TRUE., &
      &          model_type=get_my_process_name())

    ! Please note that apart from a few exceptions (e.g., ozone) 
    ! GRIB2 triplets (discipline, category, number) do net (yet) exist 
    ! for the fields that we allocate here. 
    ! For that reason, we hand over a triplet with the missing valuses (255,255,255). 
    ! In addition, 'src/upper_atmosphere/mo_upatmo_phy_config: configure_upatmo_physics' 
    ! contains a check that should make the program stop, if the user desires 
    ! the output of upper-atmosphere variables in the GRIB format.

    !------------------------------------
    !          Allocate fields
    !------------------------------------

    !------------------------------------
    !          General fields
    !   (practically always required)
    !------------------------------------

    ! &      diag%mdry(nproma,nlev,nblks_c) 
    !--------------------------------------
    ! Construct variable name
    vn_pfx_len = LEN_TRIM(vname_prefix)
    var_name   = vname_prefix(1:vn_pfx_len)//'mdry'
    cf_desc    = t_cf_var(var_name, 'kg m-2', 'mass of dry air', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, var_name, diag%mdry,                             &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
      &           ldims=shape3d,                                              &
      &           vert_interp=create_vert_interp_metadata(                    &
      &                       vert_intp_type=vintp_types("P","Z","I"),        &
      &                       vert_intp_method=VINTP_METHOD_LIN),             &
      &           loutput=.TRUE., lrestart=.TRUE.                             )         
    
    ! &      diag%amd(nproma,nlev,nblks_c) 
    !-------------------------------------
    var_name   = vname_prefix(1:vn_pfx_len)//'amd'
    cf_desc    = t_cf_var(var_name, 'g mol-1', 'molar mass of dry air', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, var_name, diag%amd,                              &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
      &           ldims=shape3d,                                              &
      &           vert_interp=create_vert_interp_metadata(                    &
      &                       vert_intp_type=vintp_types("P","Z","I"),        &
      &                       vert_intp_method=VINTP_METHOD_LIN),             &
      &           loutput=.TRUE., lrestart=.TRUE.                             )         
    
    ! &      diag%cpair(nproma,nlev,nblks_c) 
    !---------------------------------------
    var_name   = vname_prefix(1:vn_pfx_len)//'cpair'
    cf_desc    = t_cf_var(var_name, 'J K-1 kg-1',                               &
      &                   'heat capacity of (moist) air at constant pressure',  &
      &                   datatype_flt                                      )
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, var_name, diag%cpair,                            &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
      &           ldims=shape3d,                                              &
      &           vert_interp=create_vert_interp_metadata(                    &
      &                       vert_intp_type=vintp_types("P","Z","I"),        &
      &                       vert_intp_method=VINTP_METHOD_LIN),             &
      &           loutput=.TRUE., lrestart=.TRUE.                             )
    
    ! &      diag%grav(nproma,nlev,nblks_c) 
    !--------------------------------------
    var_name   = vname_prefix(1:vn_pfx_len)//'grav'
    cf_desc    = t_cf_var(var_name, 'm s-2',                                  &
      &                   'gravitational acceleration of Earth',              &
      &                   datatype_flt                                        )
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, var_name, diag%grav,                             &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
      &           ldims=shape3d,                                              &
      &           vert_interp=create_vert_interp_metadata(                    &
      &                       vert_intp_type=vintp_types("P","Z","I"),        &
      &                       vert_intp_method=VINTP_METHOD_LIN),             &
      &           loutput=.TRUE., lrestart=.TRUE.                             )     

    !------------------------------------
    !  Radiation and chemical heating
    !------------------------------------

    IF (upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%enabled )) THEN

      ! The scaling and efficiency factors are actually only required, 
      ! if the radiation group is switched on.  

      ! &      diag%sclrlw(nproma,nlev,nblks_c) 
      !----------------------------------------
      var_name   = vname_prefix(1:vn_pfx_len)//'sclrlw'
      cf_desc    = t_cf_var(var_name, '1',                                          &
        &                   'scaling factor for long-wave radiation heating rate',  &
        &                   datatype_flt                                            )
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, var_name, diag%sclrlw,                           &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d,                                              &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           loutput=.TRUE., lrestart=.TRUE.                             )         

      ! &      diag%effrsw(nproma,nlev,nblks_c) 
      !----------------------------------------
      var_name   = vname_prefix(1:vn_pfx_len)//'effrsw'
      cf_desc    = t_cf_var(var_name, '1',                                              &
        &                   'efficiency factor for short-wave radiation heating rate',  &
        &                   datatype_flt                                                )
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, var_name, diag%effrsw,                           &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d,                                              &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           loutput=.TRUE., lrestart=.TRUE.                             )         

    ENDIF !Radiation switched on?

    !------------------------------------
    !     Radiatively active gases
    !------------------------------------

    IF (upatmo_nwp_phy_config%l_gas_stat( iUpatmoGasStat%enabled )) THEN

      shape4d = (/nproma, nlev, nblks_c, ngas/)

      ALLOCATE( diag%gas_ptr( ngas ), STAT=istat )
      IF (istat/=SUCCESS) THEN
        WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_diag(', jg, &
             ')%gas_ptr failed'
        CALL finish(routine, message_text)
      ENDIF

      ! &      diag%gas(nproma,nlev,nblks_c,ngas)
      !------------------------------------------
      var_name   = vname_prefix(1:vn_pfx_len)//'rad_gases'
      cf_desc    = t_cf_var(var_name, ' ',                                      &
        &                   'radiatively active gases', datatype_flt            )
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, var_name,                                        &
        &           diag%gas,                                                   &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape4d,                                              &
        &           lcontainer=.TRUE., loutput=.FALSE., lrestart=.FALSE.        )   

      var_name_ref(1:vn_pfx_len) = vname_prefix(1:vn_pfx_len)
      DO jgas = 1, ngas
        
        ! &      diag%gas(nproma,nlev,nblks_c,ngas)
        !------------------------------------------
        var_name_ref(vn_pfx_len+1:) = upatmo_nwp_phy_config%gas(jgas)%name
        cf_desc    = t_cf_var(var_name_ref, &
             upatmo_nwp_phy_config%gas(jgas)%unit, &
             upatmo_nwp_phy_config%gas(jgas)%longname, datatype_flt)
        grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( diag_list, var_name, var_name_ref,                          & 
          &           diag%gas_ptr( jgas )%p_3d,                                  &
          &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
          &           ldims=shape3d,                                              &
          &           vert_interp=create_vert_interp_metadata(                    &
          &                       vert_intp_type=vintp_types("P","Z","I"),        &
          &                       vert_intp_method=VINTP_METHOD_LIN),             &
          &           loutput=.TRUE., lrestart=.TRUE.,                            &
          &           in_group=groups("upatmo_rad_gases"),                        &
          &           opt_var_ref_pos=4, ref_idx=jgas                             )   

      ENDDO  !jgas
      
    ENDIF !Gases required?

  END SUBROUTINE new_upatmo_diag_list

  !==================================================================================== 

  !>
  !! Allocation of tendencies from physics parameterizations.
  !!
  SUBROUTINE new_upatmo_tend_list( jg, nlev, nblks_c, nblks_e, nproma, listname, vname_prefix, &
    &                              upatmo_nwp_phy_config, tend_list, tend, lmessage            )

    ! In/out variables
    INTEGER,                INTENT(IN)    :: jg, nlev, nblks_c, nblks_e, nproma
    CHARACTER(LEN=*),       INTENT(IN)    :: listname
    CHARACTER(LEN=*),       INTENT(IN)    :: vname_prefix
    TYPE(t_upatmo_nwp_phy), INTENT(IN)    :: upatmo_nwp_phy_config
    TYPE(t_var_list_ptr),       INTENT(INOUT) :: tend_list
    TYPE(t_upatmo_tend),    INTENT(INOUT) :: tend
    LOGICAL,                INTENT(IN)    :: lmessage
    TARGET :: upatmo_nwp_phy_config

    ! Local variables 
    INTEGER :: shape3d_c(3), shape3d_e(3), shape4d_c(4)
    INTEGER :: ibits
    INTEGER :: datatype_flt
    INTEGER :: jtrc, jtnd, jst, jgrp, nstate
    INTEGER :: istartlev, iendlev
    INTEGER :: istat
    INTEGER :: vn_pfx_len, tlen

    LOGICAL :: loutput, ltend( iUpatmoTendId%nitem_2 )

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    CHARACTER(len=2) :: cjst
    CHARACTER(len=9) :: ctrc
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: &
      & var_name_prefix, var_name, var_dscrptn, var_name_ref

    CHARACTER(len=*), PARAMETER :: unit_m_over_s2 = 'm s-2', &
         unit_k_over_s = 'K s-1', &
         unit_kg_over_kg_s = 'kg kg-1 s-1'
    CHARACTER(LEN = :), POINTER :: prc_name_vdfmol, &
      &                            prc_name_fric, &
      &                            prc_name_iondrag, &
      &                            prc_name_joule, &
      &                            prc_name_srbc, &
      &                            prc_name_nlte, &
      &                            prc_name_euv, &
      &                            prc_name_no, &
      &                            prc_name_chemheat

    INTEGER, PARAMETER :: ntrc   = iUpatmoTracerId%nitem
    INTEGER, PARAMETER :: ngrp   = iUpatmoGrpId%nitem
    INTEGER, PARAMETER :: ntnd   = iUpatmoTendId%nitem    ! (Excludes Exner pressure)
    INTEGER, PARAMETER :: ntnd_2 = iUpatmoTendId%nitem_2  ! (Includes Exner pressure)
    CHARACTER(LEN=*), PARAMETER ::  &
      &  routine = modname//':new_upatmo_tend_list'

    !---------------------------------------------------------

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 

    shape3d_c = (/nproma, nlev, nblks_c/)
    shape3d_e = (/nproma, nlev, nblks_e/)
    shape4d_c = (/nproma, nlev, nblks_c, ntrc/)

    ! Ensure that all pointers have a defined association status
    ! (should be covered by the '=> NULL()'-statement at their 
    ! declaration, but to be on the safe side ...)
    NULLIFY( tend%ddt_temp_srbc,     &
      &      tend%ddt_temp_nlte,     &
      &      tend%ddt_temp_euv,      &
      &      tend%ddt_temp_vdfmol,   &
      &      tend%ddt_temp_fric,     &
      &      tend%ddt_temp_no,       &
      &      tend%ddt_temp_chemheat, &
      &      tend%ddt_temp_joule,    &
      &      tend%ddt_u_vdfmol,      &
      &      tend%ddt_u_iondrag,     &
      &      tend%ddt_v_vdfmol,      &
      &      tend%ddt_v_iondrag,     &
      &      tend%ddt_qx_vdfmol      )

    ! Register the field list and apply default settings
    ! The tendencies will be written in the restart file as the case may be by default, 
    ! for the special way the total tendencies are computed in 
    ! 'mo_nwp_upatmo_interface: nwp_upatmo_interface' alone

    CALL vlr_add(tend_list, listname, patch_id=jg, lrestart=.TRUE., &
      &          model_type=get_my_process_name())

    ! Please note that for most of the tendencies, which we allocate here, 
    ! either a GRIB2 triplet (discipline, category, number) does not (yet) exist 
    ! or is already occupied by tendencies of the "standard" configuration of ICON.
    ! For that reason, we hand over a triplet with the missing valuses (255,255,255). 
    ! In addition, 'src/upper_atmosphere/mo_upatmo_phy_config: configure_upatmo_physics' 
    ! contains a check that should make the program stop, if the user desires 
    ! the output of upper-atmosphere variables in the GRIB format.

    !------------------------------------
    !          Allocate fields
    !------------------------------------

    ! Note: we could save a significant amount of memory, 
    ! if we would limit the allocation of the tendencies 
    ! to only those vertical grid layers,
    ! for which the processes compute tendencies. 
    ! However, this is not straightforward for at least 
    ! the following reason:
    ! * For output purposes, the vertical grid type 
    !   has to be specified in 'add_var'. 
    !   The standard entry for the 'nlev' model levels is: 
    !   - CALL add_var(..., vgrid = ZA_REFERENCE, ...)
    !   So for each process we would have to register 
    !   a new vertical grid type 'ZA_REFERENCE_<process>', say, 
    !   in 'src/shared/mo_zaxis_type', 
    !   and add a definition in 
    !   'src/io/shared/mo_name_list_output_zaxes: setup_ml_axes_atmo'.
    !   Now, 'setup_ml_axes_atmo' seems to be called in 
    !   'src/io/shared/mo_name_list_output_init: create_vertical_axes', 
    !   which in turn seems to be called in 
    !   'src/io/shared/mo_name_list_output: name_list_io_main_proc', 
    !   which in turn seems to be called in
    !   'src/drivers/mo_atmo_model: construct_atmo_model', 
    !   which in turn is calle in 
    !   'src/drivers/mo_atmo_model: atmo_model' 
    !   and now we finally come to the point, BEFORE the call of 'atmo_nonhydrostatic'!
    !   That is, the information, which we would need for the definition 
    !   of the new vertical grid types is not yet available. 
    !   Any modification to this infrastructure is too delicate, 
    !   so we refrain from it. 
    !   To cut a long story short, we have no choice but to allocate the tendencies 
    !   for the entire range of model levels.

    ! The tendencies from the single processes can be selected for output
    loutput = .TRUE.

    ! actual useful length of vname_prefix
    vn_pfx_len = LEN_TRIM(vname_prefix)
    !------------------------------------
    !  Radiation and chemical heating
    !------------------------------------

    IF (upatmo_nwp_phy_config%grp( iUpatmoGrpId%rad )%l_stat( iUpatmoPrcStat%enabled )) THEN

      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%srbc )%name)
      prc_name_srbc     => upatmo_nwp_phy_config%prc( iUpatmoPrcId%srbc )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%nlte )%name)
      prc_name_nlte     => upatmo_nwp_phy_config%prc( iUpatmoPrcId%nlte )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%euv )%name)
      prc_name_euv      => upatmo_nwp_phy_config%prc( iUpatmoPrcId%euv )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%no )%name)
      prc_name_no       => upatmo_nwp_phy_config%prc( iUpatmoPrcId%no )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%chemheat )%name)
      prc_name_chemheat => upatmo_nwp_phy_config%prc( iUpatmoPrcId%chemheat )%name(1:tlen)

      !------------------------------------
      !      Temperature tendencies
      !------------------------------------

      ! Prefix of variable name
      var_name_prefix = vname_prefix(1:vn_pfx_len)//'ddt_temp_'

      ! &      tend%ddt_temp_srbc(nproma,nlev,nblks_c) 
      !-----------------------------------------------
      ! Construc variable name
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_srbc
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                             &
        &                   'temperature tendency due to absorbtion by O2 in SRB and SRC',  &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_srbc,                    &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )   

      ! &      tend%ddt_temp_nlte(nproma,nlev,nblks_c) 
      !--------------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_nlte
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to Non-LTE heating',      &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_nlte,                    &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )   

      ! &      tend%ddt_temp_euv(nproma,nlev,nblks_c) 
      !----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_euv
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to EUV heating',          &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_euv,                     &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )  

      ! &      tend%ddt_temp_no(nproma,nlev,nblks_c) 
      !---------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_no
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to NO heating at NIR',    &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_no,                      &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )  

      ! &  tend%ddt_temp_chemheat(nproma,nlev,nblks_c) 
      !-----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_chemheat
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to chemical heating',     &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_chemheat,                &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             ) 

    ENDIF !RAD-group switched on?

    !------------------------------------
    !           Ion drag (I), 
    !       molecular diffusion (M) 
    !     and frictional heating (F)
    !------------------------------------

    IF (upatmo_nwp_phy_config%grp( iUpatmoGrpId%imf )%l_stat( iUpatmoPrcStat%enabled )) THEN

      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%vdfmol )%name)
      prc_name_vdfmol  => upatmo_nwp_phy_config%prc( iUpatmoPrcId%vdfmol )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%fric )%name)
      prc_name_fric    => upatmo_nwp_phy_config%prc( iUpatmoPrcId%fric )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%iondrag )%name)
      prc_name_iondrag => upatmo_nwp_phy_config%prc( iUpatmoPrcId%iondrag )%name(1:tlen)
      tlen = LEN_TRIM(upatmo_nwp_phy_config%prc( iUpatmoPrcId%joule )%name)
      prc_name_joule   => upatmo_nwp_phy_config%prc( iUpatmoPrcId%joule )%name(1:tlen)

      !------------------------------------
      !      Temperature tendencies
      !------------------------------------

      var_name_prefix = vname_prefix(1:vn_pfx_len)//'ddt_temp_'

      ! &      tend%ddt_temp_vdfmol(nproma,nlev,nblks_c) 
      !-------------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_vdfmol
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to molecular diffusion',  &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_vdfmol,                  &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )  

      ! &      tend%ddt_temp_fric(nproma,nlev,nblks_c) 
      !-----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_fric
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to frictional heating',   &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_fric,                    &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &  
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )

      ! &      tend%ddt_temp_joule(nproma,nlev,nblks_c) 
      !------------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+9)//prc_name_joule
      cf_desc    = t_cf_var(var_name, unit_k_over_s,                                 &
        &                   'temperature tendency due to Joule heating',        &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_temp_joule,                   &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )  

      !------------------------------------
      !         U-wind tendencies
      !------------------------------------

      var_name_prefix = vname_prefix(1:vn_pfx_len)//'ddt_u_'

      ! &      tend%ddt_u_vdfmol(nproma,nlev,nblks_c) 
      !----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+6)//prc_name_vdfmol
      cf_desc    = t_cf_var(var_name, unit_m_over_s2,                                 &
        &                   'u-wind tendency due to molecular diffusion',       &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_u_vdfmol,                     &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )   

      ! &      tend%ddt_u_iondrag(nproma,nlev,nblks_c) 
      !-----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+6)//prc_name_iondrag
      cf_desc    = t_cf_var(var_name, unit_m_over_s2,                                 &
        &                   'u-wind tendency due to ion drag',                  &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_u_iondrag,                    &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )  

      !------------------------------------
      !         V-wind tendencies
      !------------------------------------

      var_name_prefix = vname_prefix(1:vn_pfx_len)//'ddt_v_'

      ! &      tend%ddt_v_vdfmol(nproma,nlev,nblks_c) 
      !----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+6)//prc_name_vdfmol
      cf_desc    = t_cf_var(var_name, unit_m_over_s2,                                 &
        &                   'v-wind tendency due to molecular diffusion',       &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_v_vdfmol,                     &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )  

      ! &      tend%ddt_v_iondrag(nproma,nlev,nblks_c) 
      !-----------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+6)//prc_name_iondrag
      cf_desc    = t_cf_var(var_name, unit_m_over_s2,                                 &
        &                   'v-wind tendency due to ion drag',                  &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_v_iondrag,                    &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape3d_c,                                            &
        &           vert_interp=create_vert_interp_metadata(                    &
        &                       vert_intp_type=vintp_types("P","Z","I"),        &
        &                       vert_intp_method=VINTP_METHOD_LIN),             &
        &           in_group=groups("upatmo_tendencies"),                       &
        &           loutput=loutput                                             )

      !------------------------------------
      !         Tracer tendencies
      !------------------------------------

      ! Currently only specific humidity
      ALLOCATE(tend%ddt_qx_vdfmol_ptr( ntrc ), STAT=istat)
      IF (istat/=SUCCESS) THEN
        WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', jg, &
             ')%ddt_qx_vdfmol_ptr failed'
        CALL finish(routine, message_text)
      ENDIF

      var_name_prefix = vname_prefix(1:vn_pfx_len)//'ddt_q'

      ! &     tend%ddt_qx_vdfmol(nproma,nlev,nblks_c,ntrc) 
      !---------------------------------------------------
      var_name   = var_name_prefix(1:vn_pfx_len+5)//'x_'//prc_name_vdfmol
      cf_desc    = t_cf_var(var_name, unit_kg_over_kg_s,                                 &
        &                   'tendencies of mass mixing ratio of tracers '//     &
        &                   'due to molecular diffusion',                       &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( tend_list, var_name, tend%ddt_qx_vdfmol,                    &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
        &           ldims=shape4d_c,                                            &
        &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.        )  

      jtrc = iUpatmoTracerId%qv

      ! &      tend%ddt_qv_vdfmol(nproma,nlev,nblks_c) 
      !-----------------------------------------------
      var_name_ref = var_name_prefix(1:vn_pfx_len+5)//'v_'//prc_name_vdfmol
      cf_desc    = t_cf_var(var_name_ref, unit_kg_over_kg_s,                                   &
        &                   'tendency of specific humidity '//                        &
        &                   'due to molecular diffusion',                             &
        &                   datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_ref( tend_list, var_name, var_name_ref,                                & 
        &           tend%ddt_qx_vdfmol_ptr( jtrc )%p_3d,                              &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
        &           ldims=shape3d_c,                                                  &
        &           vert_interp=create_vert_interp_metadata(                          &
        &                       vert_intp_type=vintp_types("P","Z","I"),              &
        &                       vert_intp_method=VINTP_METHOD_LIN),                   &
        &           in_group=groups("upatmo_tendencies"),                             &
        &           loutput=loutput, lrestart=.TRUE., opt_var_ref_pos=4, ref_idx=jtrc )

    ENDIF !IMF-group switched on?

    !------------------------------------
    !      Accumulative tendencies
    !------------------------------------

    !-------------
    ! Preparation
    !-------------

    ! For convenience, 'ddt%info' and 'ddt%state' are allocated, 
    ! even if the accumulative tendencies are not required 
    ! (due to an offline-mode).

    ALLOCATE(tend%ddt%info( ntnd_2 ), STAT=istat)
    IF (istat/=SUCCESS) THEN
      WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', jg, &
           ')%ddt%info failed'
      CALL finish(routine, message_text)
    ENDIF
    
    ALLOCATE(tend%ddt%state( ntnd_2 ), STAT=istat)
    IF (istat/=SUCCESS) THEN
      WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', jg, &
           ')%ddt%state failed'
      CALL finish(routine, message_text)
    ENDIF

    ! Name and number of states
    ! For temperature only the current state is required
    jtnd = iUpatmoTendId%temp
    tend%ddt%info( jtnd )%name     = 'ddt_temp_tot'
    tend%ddt%info( jtnd )%longname = 'accumulative temperature tendency'
    tend%ddt%info( jtnd )%unit     = 'K s-1'
    tend%ddt%info( jtnd )%nstate   = 1
    ! For Exner pressure the previous and current states are required
    jtnd = iUpatmoTendId%exner
    tend%ddt%info( jtnd )%name     = 'ddt_exner_tot'
    tend%ddt%info( jtnd )%longname = 'accumulative Exner pressure tendency'
    tend%ddt%info( jtnd )%unit     = '1'
    tend%ddt%info( jtnd )%nstate   = 2
    ! For horizontal wind the previous and current states are required
    jtnd = iUpatmoTendId%wind_h
    tend%ddt%info( jtnd )%name     = 'ddt_vn_tot'
    tend%ddt%info( jtnd )%longname = 'accumulative wind tendency'
    tend%ddt%info( jtnd )%unit     = 'm s-2'
    tend%ddt%info( jtnd )%nstate   = 2
    ! For tracers only the current state is required
    jtnd = iUpatmoTendId%qx
    tend%ddt%info( jtnd )%name     = 'ddt_qx_tot'
    tend%ddt%info( jtnd )%longname = 'accumulative tracer tendency'
    tend%ddt%info( jtnd )%unit     = 'kg kg-1 s-1'
    tend%ddt%info( jtnd )%nstate   = 1
    
    ! Initialize state transition type
    DO jtnd = 1, ntnd_2
      CALL tend%ddt%state( jtnd )%init( nstate                  = tend%ddt%info( jtnd )%nstate, & !in
        &                               optDefaultFinishOnError = .TRUE.,                       & !optin
        &                               optDefaultLocking       = .TRUE.                        ) !optin
    ENDDO  !jtnd

    ! The following should only be necessary, 
    ! if any of the upper-atmosphere physics groups
    ! is in interactive mode. 

    IF (ANY(upatmo_nwp_phy_config%l_any_update( : ))) THEN

      ! Output is disabled for the accumulative tendencies 
      loutput = .FALSE.

      ltend( : ) = .FALSE.

      DO jtnd = 1, ntnd  ! Exclude Exner pressure!
        ! Determine vertical start and end level
        DO jgrp = 1, ngrp
          IF (upatmo_nwp_phy_config%grp( jgrp )%l_update( jtnd )) THEN
            istartlev = upatmo_nwp_phy_config%grp( jgrp )%istartlev
            iendlev   = upatmo_nwp_phy_config%grp( jgrp )%iendlev
            IF (jgrp == 1) THEN
              tend%ddt%info( jtnd )%istartlev = istartlev
              tend%ddt%info( jtnd )%iendlev   = iendlev
            ELSE
              IF (istartlev < tend%ddt%info( jtnd )%istartlev) THEN 
                tend%ddt%info( jtnd )%istartlev = istartlev
              ENDIF
              IF (iendlev > tend%ddt%info( jtnd )%iendlev) THEN
                tend%ddt%info( jtnd )%iendlev = iendlev
              ENDIF
            ENDIF  !IF (jgrp == 1)
          ENDIF  !IF (upatmo_nwp_phy_config%grp( jgrp )%l_update( ntnd ))
        ENDDO  !jgrp
        ! Will tendency experience updates?
        ltend( jtnd ) = upatmo_nwp_phy_config%l_any_update( jtnd )
      ENDDO  !jtnd
      ! For Exner pressure:
      ! (Assumes proper initialization of 'istartlev' and 'iendlev' in 'src/upper_atmosphere/mo_upatmo_types'!)
      tend%ddt%info( iUpatmoTendId%exner )%istartlev = MIN( tend%ddt%info( iUpatmoTendId%temp )%istartlev, &
        &                                                   tend%ddt%info( iUpatmoTendId%qx )%istartlev    )
      tend%ddt%info( iUpatmoTendId%exner )%iendlev   = MAX( tend%ddt%info( iUpatmoTendId%temp )%iendlev, &
        &                                                   tend%ddt%info( iUpatmoTendId%qx )%iendlev    )
      ! 'l_any_update' and 'l_update' contain only information 
      ! on the basic tendencies: temp, u, v, qx, 
      ! but not on the derived tendency of the Exner pressure
      ltend( iUpatmoTendId%exner ) = ltend( iUpatmoTendId%temp ) .OR. ltend( iUpatmoTendId%qx )

      IF (lmessage) THEN
        DO jtnd = 1, ntnd_2
          IF (ltend( jtnd )) THEN
            WRITE (message_text, '(3a,i0)') 'Start level of ', &
                 TRIM(tend%ddt%info( jtnd )%longname), ': ', tend%ddt%info( jtnd )%istartlev
            CALL message (routine, message_text)
            WRITE (message_text, '(3a,i0)') 'End level of ', &
                 TRIM(tend%ddt%info( jtnd )%longname), ': ', tend%ddt%info( jtnd )%iendlev
            CALL message (routine, message_text)
          ENDIF  !IF (ltend( jtnd )
        ENDDO  !jtnd
      ENDIF  !IF (lmessage)

      !-------------
      ! Temperature
      !-------------

      jtnd = iUpatmoTendId%temp

      IF (ltend( jtnd )) THEN
        
        nstate = tend%ddt%info( jtnd )%nstate

        ALLOCATE(tend%ddt%temp( nstate ), STAT=istat)
        IF (istat/=SUCCESS) THEN
          WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', &
               jg, ')%ddt%temp failed'
          CALL finish(routine, message_text)
        ENDIF
        
        ! Loop over states
        DO jst = 1, nstate

          WRITE (cjst, '(i2.2)') jst

          NULLIFY(tend%ddt%temp( jst )%tot)

          ! &      tend%ddt%temp(jst)%tot(nproma,nlev,nblks_c)
          !---------------------------------------------------
          ! Construct variable name for output
          var_name = vname_prefix(1:vn_pfx_len)//TRIM(tend%ddt%info( jtnd )%name)//STATELEVEL_SUFFIX//cjst
          ! Variable description
          var_dscrptn = TRIM(tend%ddt%info( jtnd )%longname)//' of time level with index '//cjst
          ! Variable unit
          cf_desc    = t_cf_var(var_name, tend%ddt%info( jtnd )%unit, &
               var_dscrptn, datatype_flt)
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, var_name, tend%ddt%temp( jst )%tot,              &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
            &           ldims=shape3d_c, loutput=loutput                            )  
          
        ENDDO  !jst

      ENDIF !Any update of temperature?        

      !----------------
      ! Exner pressure
      !----------------

      jtnd = iUpatmoTendId%exner

      IF (ltend( jtnd )) THEN

        nstate = tend%ddt%info( jtnd )%nstate

        ALLOCATE(tend%ddt%exner( nstate ), STAT=istat)
        IF (istat/=SUCCESS) THEN
          WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', &
               jg, ')%ddt%exner failed'
          CALL finish(routine, message_text)
        ENDIF

        DO jst = 1, nstate

          WRITE (cjst, '(i2.2)') jst

          NULLIFY(tend%ddt%exner( jst )%tot)

          ! &      tend%ddt%exner(jst)%tot(nproma,nlev,nblks_c)
          !----------------------------------------------------
          var_name   = vname_prefix(1:vn_pfx_len)//TRIM(tend%ddt%info( jtnd )%name)//STATELEVEL_SUFFIX//cjst
          var_dscrptn = TRIM(tend%ddt%info( jtnd )%longname)//' of time level with index '//cjst
          cf_desc    = t_cf_var(var_name, tend%ddt%info( jtnd )%unit, var_dscrptn, datatype_flt)
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, var_name, tend%ddt%exner( jst )%tot,             &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
            &           ldims=shape3d_c, loutput=loutput                            )  
          
        ENDDO  !jst

      ENDIF !Any update of Exner pressure? 

      !-----------------
      ! Horizontal wind
      !-----------------

      jtnd = iUpatmoTendId%wind_h

      IF (ltend( jtnd )) THEN

        nstate = tend%ddt%info( jtnd )%nstate

        ALLOCATE(tend%ddt%vn( nstate ), STAT=istat)
        IF (istat/=SUCCESS) THEN
          WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', &
               jg, ')%ddt%vn failed'
          CALL finish(routine, message_text)
        ENDIF

        DO jst = 1, nstate

          WRITE (cjst, '(i2.2)') jst

          NULLIFY(tend%ddt%vn( jst )%tot)    
          
          ! &      tend%ddt%vn(jst)%tot(nproma,nlev,nblks_e)
          !-------------------------------------------------
          var_name   = vname_prefix(1:vn_pfx_len)//TRIM(tend%ddt%info( jtnd )%name)//STATELEVEL_SUFFIX//cjst
          var_dscrptn = TRIM(tend%ddt%info( jtnd )%longname)//' of time level with index '//cjst
          cf_desc    = t_cf_var(var_name, tend%ddt%info( jtnd )%unit, var_dscrptn, datatype_flt)
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
          CALL add_var( tend_list, var_name, tend%ddt%vn( jst )%tot,                &
            &           GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
            &           ldims=shape3d_e, loutput=loutput                            )  
          
        ENDDO  !jst
        
      ENDIF !Any update of horizontal wind component?

      !---------
      ! Tracers
      !---------

      jtnd = iUpatmoTendId%qx

      IF (ltend( jtnd )) THEN

        nstate = tend%ddt%info( jtnd )%nstate

        ALLOCATE(tend%ddt%qx( nstate ), STAT=istat)
        IF (istat/=SUCCESS) THEN
          WRITE (message_text, '(a,i0,a)') 'Allocation of prm_upatmo_tend(', &
               jg, ')%ddt%qx failed'
          CALL finish(routine, message_text)
        ENDIF

        DO jst = 1, nstate

          WRITE (cjst, '(i2.2)') jst

          ALLOCATE(tend%ddt%qx( jst )%tot_ptr( ntrc ), STAT=istat)
          IF (istat/=SUCCESS) THEN
            WRITE (message_text, '(a,i0,a,i0,a)') 'Allocation of prm_upatmo_tend(', &
                 jg, ')%dqx%dt(', jst, ')%tot_ptr failed'
            CALL finish (routine, message_text)
          ENDIF

          NULLIFY(tend%ddt%qx( jst )%tot)

          ! &     tend%ddt%qx(jst)%tot(nproma,nlev,nblks_c,ntrc)
          !-----------------------------------------------------
          var_name   = vname_prefix(1:vn_pfx_len)//TRIM(tend%ddt%info( jtnd )%name)//STATELEVEL_SUFFIX//cjst
          var_dscrptn = TRIM(tend%ddt%info( jtnd )%longname)//' of time level with index '//cjst
          cf_desc    = t_cf_var(var_name, tend%ddt%info( jtnd )%unit, &
               var_dscrptn, datatype_flt)
          grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( tend_list, var_name, tend%ddt%qx( jst )%tot,                  &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
            &           ldims=shape4d_c,                                              &
            &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.          )

          DO jtrc = 1, ntrc

            IF (jtrc == iUpatmoTracerId%qv) THEN
              ctrc = 'v'
            ELSE
              WRITE (ctrc, '(i0)') jtrc
            ENDIF

            ! &     tend%ddt%qx(jst)%tot(nproma,nlev,nblks_c)
            !------------------------------------------------
            var_name_ref = vname_prefix(1:vn_pfx_len)//'ddt_q'//TRIM(ctrc)//STATELEVEL_SUFFIX//cjst
            var_dscrptn  = 'accumulative tendency of tracer q'//TRIM(ctrc) &
              & //' of time level with index '//cjst
            cf_desc    = t_cf_var(var_name_ref, tend%ddt%info( jtnd )%unit, &
                 var_dscrptn, datatype_flt)
            grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
            CALL add_ref( tend_list, var_name, var_name_ref,                                &
              &           tend%ddt%qx( jst )%tot_ptr( jtrc )%p_3d,                          &
              &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
              &           ldims=shape3d_c,                                                  &
              &           loutput=loutput, lrestart=.TRUE., opt_var_ref_pos=4, ref_idx=jtrc )   
            
          ENDDO  !jtrc
          
        ENDDO  !jst
        
      ENDIF !Any update of tracers?

    ENDIF  !Accumulative tendencies required?

  END SUBROUTINE new_upatmo_tend_list
    
END MODULE mo_upatmo_state
