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

MODULE mo_comin_adapter

#ifndef __NO_ICON_COMIN__
  USE mo_kind,                    ONLY : wp
  USE mo_impl_constants,          ONLY : SUCCESS, vname_len, TLEV_NNOW_RCF, MIURA, ippm_v, &
    &                                    max_ntracer,  ifluxl_sm, islopel_vsm
  USE mo_cdi_constants,           ONLY : GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT,   &
    &                                    GRID_UNSTRUCTURED_EDGE
  USE mo_zaxis_type,              ONLY : ZA_REFERENCE, ZA_REFERENCE_HALF, ZA_REFERENCE_HALF_HHL, &
    & ZA_SURFACE, zaxisTypeList, t_zaxisType
  USE mo_master_control,          ONLY : get_my_process_name
  USE mo_grid_config,             ONLY : n_dom
  USE mo_run_config,              ONLY : msg_level
  USE mo_exception,               ONLY : message, message_text, finish
  USE mo_model_domain,            ONLY : t_patch, p_patch
  USE mo_parallel_config,         ONLY : nproma
  USE mo_var_list_register,       ONLY : t_vl_register_iter, vlr_add
  USE mo_var_metadata,            ONLY : get_var_timelevel, get_var_name, &
    &                                    get_timelevel_string
  USE mo_var_metadata_types,      ONLY : t_var_metadata_dynamic, t_var_metadata
  USE mo_var,                     ONLY : t_var, level_type_ml
  USE mo_var_list,                ONLY : add_var, add_ref, t_var_list_ptr, find_list_element
  USE mo_var_groups,              ONLY : groups
  USE mo_cf_convention,           ONLY : t_cf_var
  USE mo_grib2,                   ONLY : t_grib2_var
  USE mo_nonhydro_types,          ONLY : t_nh_state, t_nh_state_lists
  USE mo_nwp_phy_state,           ONLY : prm_nwp_tend_list, prm_nwp_tend
  USE mo_advection_config,        ONLY : t_advection_config, advection_config
  USE mo_advection_utils,         ONLY : add_tracer_ref
  USE mo_tracer_metadata,         ONLY : create_tracer_metadata
  USE mo_name_list_output_types,  ONLY : t_patch_info
  USE mo_util_vgrid_types,        ONLY : t_vgrid_buffer
  USE mo_decomposition_tools,     ONLY : get_local_index
  USE mo_time_config,             ONLY : t_time_config
  USE mo_comin_config,            ONLY : comin_config, t_comin_tracer_info
  USE mo_coupling_config,         ONLY : is_coupled_run
  USE mo_coupling_utils,          ONLY : cpl_get_instance_id
  USE comin_host_interface,       ONLY : t_comin_var_ptr,                         &
    &                                    t_comin_var_descriptor,                  &
    &                                    t_var_request_list_item,                 &
    &                                    comin_request_get_list_head,             &
    &                                    comin_var_list_append,                   &
    &                                    COMIN_ZAXIS_2D, COMIN_ZAXIS_3D,          &
    &                                    COMIN_ZAXIS_3D_HALF, COMIN_ZAXIS_UNDEF,  &
    &                                    t_comin_descrdata_global,                &
    &                                    t_comin_descrdata_domain,                &
    &                                    comin_descrdata_set_global,              &
    &                                    comin_descrdata_set_domain,              &
    &                                    t_comin_descrdata_simulation_interval,   &
    &                                    comin_descrdata_set_simulation_interval, &
    &                                    comin_descrdata_set_fct_glb2loc_cell,    &
    &                                    comin_current_set_datetime,              &
    &                                    comin_metadata_set,comin_metadata_get_or,&
    &                                    t_comin_var_metadata_iterator,           &
    &                                    comin_descrdata_set_timesteplength,      &
    &                                    comin_callback_context_call,             &
    &                                    comin_setup_set_verbosity_level,         &
    &                                    comin_plugin_primaryconstructor,         &
    &                                    COMIN_HGRID_UNSTRUCTURED_CELL,           &
    &                                    COMIN_HGRID_UNSTRUCTURED_EDGE,           &
    &                                    COMIN_HGRID_UNSTRUCTURED_VERTEX,         &
    &                                    comin_var_set_sync_device_mem,           &
    &                                    COMIN_SUCCESS,                           &
    &                                    COMIN_METADATA_TYPEID_INTEGER,           &
    &                                    COMIN_METADATA_TYPEID_REAL,              &
    &                                    COMIN_METADATA_TYPEID_CHARACTER,         &
    &                                    COMIN_METADATA_TYPEID_LOGICAL

  USE mo_timer,                   ONLY : timers_level, timer_comin_init, &
                                         timer_comin_primary_constructors, &
                                         timer_comin_callbacks, &
                                         timer_start, timer_stop
  USE mo_fortran_tools,           ONLY: set_acc_host_or_device
  USE mtime,                      ONLY : datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_util_vcs,                ONLY: get_remote_url, &
       &                                get_local_branch,    &
       &                                get_revision,   &
       &                                get_icon_version
  USE iso_c_binding,              ONLY: C_LOC, C_PTR, C_NULL_PTR, C_DOUBLE
#ifdef _OPENACC
  USE openacc,                    ONLY: acc_get_device_type, acc_get_device_num,    &
    &                                   acc_get_property_string, acc_property_name, &
    &                                   acc_property_vendor, acc_property_driver,   &
    &                                   acc_device_kind, acc_device_none
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: icon_append_comin_tracer_variables
  PUBLIC :: icon_append_comin_tracer_phys_tend
  PUBLIC :: icon_append_comin_variables
  PUBLIC :: icon_expose_variables
  PUBLIC :: icon_expose_descrdata
  PUBLIC :: icon_update_current_datetime
  PUBLIC :: icon_expose_timesteplength_domain
  PUBLIC :: icon_update_expose_variables
  PUBLIC :: icon_call_callback

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_comin_adapter'

  ! variable lists (domain-wise) for additional ComIn variables
  TYPE(t_var_list_ptr), TARGET, ALLOCATABLE :: p_comin_varlist(:)

  TYPE :: t_exposed_timedep_vars_list
        TYPE(t_var), POINTER :: icon_var => NULL()
        TYPE(t_comin_var_ptr), POINTER :: comin_var_ptr => NULL()
        TYPE(t_exposed_timedep_vars_list), POINTER :: next => NULL()
  END TYPE t_exposed_timedep_vars_list
  TYPE(t_exposed_timedep_vars_list), POINTER :: exposed_timedep_vars_head => NULL()

#ifdef _OPENACC
 INTERFACE
    FUNCTION acc_deviceptr(h_ptr) BIND(C) RESULT(device_ptr)
      USE iso_c_binding
      IMPLICIT NONE
      TYPE(C_PTR), VALUE, INTENT(IN) :: h_ptr
      TYPE(C_PTR) :: device_ptr
    END FUNCTION acc_deviceptr
 END INTERFACE
#endif

CONTAINS

  !> loop over the total list of additional requested tracer variables
  !  and perform `add_tracer_ref` operations needed.
  !
  ! TODOs:
  ! - flag "lrestart"
  ! - vertical axis
  ! - OpenACC handling
  ! - GRIB2
  !
  SUBROUTINE icon_append_comin_tracer_variables(p_patch, p_nh_state, p_nh_state_lists)
    TYPE(t_patch),             INTENT(IN)       :: p_patch(:)
    TYPE(t_nh_state), POINTER, INTENT(IN)       :: p_nh_state(:)
    TYPE(t_nh_state_lists), POINTER, INTENT(IN) :: p_nh_state_lists(:)
    ! Local variables
    CHARACTER(*), PARAMETER    :: routine = modname//"::icon_append_comin_tracer_variables"
    INTEGER                    :: jg, shape3d_c(3), jt, ntl, tracer_idx
    TYPE(t_var_request_list_item), POINTER :: ptr
    TYPE(t_cf_var)             :: cf_desc
    TYPE(t_grib2_var)          :: grib2_desc
    CHARACTER(LEN=vname_len)   :: tracer_container_name
    CHARACTER(len=4)           :: suffix
    CHARACTER(LEN=vname_len+LEN(suffix)) :: tracer_name
    TYPE(t_var_list_ptr)       :: p_prog_list !< current prognostic state list
    TYPE(t_var_list_ptr)       :: p_tracer_list
    TYPE(t_advection_config), POINTER :: advconf
    LOGICAL                    :: tracer, tracer_turb, tracer_conv
    INTEGER                    :: ivadv, ihadv, ihlimit, ivlimit
    TYPE(t_comin_tracer_info), POINTER :: this_info => NULL()

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    DOM_LOOP : DO jg=1,n_dom
      shape3d_c = [ nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c ]

      ptr => comin_request_get_list_head()
      VAR_LOOP : DO WHILE (ASSOCIATED(ptr))
        ! since we are dealing with a generic linked list, we have to
        ! perform a "dynamic type-cast" into `t_comin_request_item`:
        ASSOCIATE (comin_request_item => ptr%item_value)
          ! skip if variable was not requested for this domain
          IF (comin_request_item%descriptor%id /= jg) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP
          END IF

          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer", tracer, .FALSE.)
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_turb", tracer_turb, .FALSE.)
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_conv", tracer_conv, .FALSE.)
          IF (.NOT. tracer) THEN
            IF (tracer_turb .OR. tracer_conv) THEN
              WRITE (message_text,'(A,A,A,L1,A,L1,A,L1)') 'Requested variable ', comin_request_item%descriptor%name, &
                &                                  'has tracer = ',tracer,', but tracer_turb = ',tracer_turb, &
                &                                  'and tracer_conv = ',tracer_conv
              CALL finish(routine, message_text)
            ENDIF
            ptr => ptr%next()
            CYCLE VAR_LOOP
          END IF

          IF ( .NOT. ASSOCIATED(this_info) ) THEN
            ALLOCATE(comin_config%comin_icon_domain_config(jg)%tracer_info_head)
            this_info => comin_config%comin_icon_domain_config(jg)%tracer_info_head
          ELSE
            ALLOCATE(this_info%next)
            this_info => this_info%next
          END IF

          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_hadv",   ihadv,       MIURA      )
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_vadv",   ivadv,       ippm_v     )
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_hlimit", ihlimit,     ifluxl_sm  )
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_vlimit", ivlimit,     islopel_vsm)

          ! add reference for every timelevel
          ntl = SIZE(p_nh_state_lists(jg)%tracer_list)
          DO jt = 1, ntl
            suffix = get_timelevel_string(timelevel=jt)
            tracer_container_name = 'tracer'//suffix

            ! add tracer tracer_name to container tracer_container_name
            p_prog_list   = p_nh_state_lists(jg)%prog_list(jt)
            p_tracer_list = p_nh_state_lists(jg)%tracer_list(jt)
            advconf      => advection_config(jg)
            tracer_name   = comin_request_item%descriptor%name//TRIM(suffix)
            CALL add_tracer_ref(p_prog_list, tracer_container_name,                &
              &          tracer_name, tracer_idx,                                  &
              &          p_nh_state(jg)%prog(jt)%tracer_ptr, cf_desc, grib2_desc,  &
              &          advconf,                                                  &
              &          ldims=shape3d_c,                                          &
              &          loutput=.TRUE., tlev_source=TLEV_NNOW_RCF,                &
              &          in_group=groups("comin_vars"),                            &
              &          tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,     &
              &                        name        = tracer_name,                  &
              &                        lfeedback   = .FALSE.,                      &
              &                        ihadv_tracer=ihadv,                         &
              &                        ivadv_tracer=ivadv,                         &
              &                        lturb_tracer=tracer_turb,                   &
              &                        lconv_tracer=tracer_conv)  )
            ! update metadata as set via ComIn
            ! note that ihadv_tracer and ivadv_tracer are also set in add_tracer_ref
            advconf%ihadv_tracer(tracer_idx) = ihadv
            advconf%ivadv_tracer(tracer_idx) = ivadv
            advconf%itype_hlimit(tracer_idx) = ihlimit
            advconf%itype_vlimit(tracer_idx) = ivlimit
          END DO

          this_info%name       = TRIM(comin_request_item%descriptor%name)
          this_info%idx_tracer = tracer_idx ! Assumes that time levels have the same number of tracers

          IF (msg_level >= 10) THEN
            WRITE (message_text,'(A,A,A,I2)') "Add tracer variable '", &
                 &  TRIM(comin_request_item%descriptor%name),             &
                 &  "', requested by third party plugins for domain ", jg
            CALL message(routine, message_text)
          ENDIF

        END ASSOCIATE
        ptr => ptr%next()
      END DO VAR_LOOP

      NULLIFY(this_info)
    END DO DOM_LOOP
    IF (timers_level > 2) CALL timer_stop(timer_comin_init)

  END SUBROUTINE icon_append_comin_tracer_variables

  !> Append physical tendencies to the respective list of ICON
  !  (Necessary for getting convective and turbulent tendencies)
  !
  ! TODOs:
  ! - vertical axis
  ! - OpenACC handling
  ! - GRIB2, CF descriptors
  !
  SUBROUTINE icon_append_comin_tracer_phys_tend(p_patch)
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    ! Local variables
    CHARACTER(*), PARAMETER   :: routine = modname//"::icon_append_comin_tracer_phys_tend"
    CHARACTER(LEN=vname_len)  :: name_tend
    INTEGER ::            &
      &  tracer_tend_idx, & !< Index of current tendency in container
      &  nturb_tracer,    & !< Number of turbulent tracers from ComIn (for crosscheck)
      &  shape3d_c(3),    & !< 3D (hor.+vert.) variable shape
      &  jg                 !< Current domain
    LOGICAL ::            &
      &  tracer,          & !< Is the variable a tracer?
      &  tracer_turb,     & !< Does the tracer participate in turbulent transport?
      &  tracer_conv        !< Does the tracer participate in convective transport?
    TYPE(t_var), POINTER :: &
      &  target_var         !< Pointer to tendency container in phys. tend. list
    TYPE(t_var_metadata), POINTER :: &
      &  target_info        !< Pointer to metadata of target_var
    TYPE(t_var_request_list_item), POINTER :: &
      &  ptr                !< Current comin request item
    TYPE(t_cf_var) :: &
      &  cf_desc            !< NETCDF variable description
    TYPE(t_grib2_var) :: &
      &  grib2_desc         !< GRIB2 variable description
    TYPE(t_comin_tracer_info), POINTER :: &
      &  this_info          !< ICON-internal ComIn variable information

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    DOM_LOOP : DO jg=1,n_dom
      shape3d_c = [ nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c ]
      ptr => comin_request_get_list_head()
      VAR_LOOP : DO WHILE (ASSOCIATED(ptr))
        ASSOCIATE (comin_request_item => ptr%item_value)
          ! skip if variable was not requested for this domain
          IF (comin_request_item%descriptor%id /= jg) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP
          END IF

          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer", tracer, .FALSE.)
          IF (.NOT. tracer) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP
          END IF

          ! Adding tendencies of tracers due to turbulence
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_turb", tracer_turb, .FALSE.)
          IF (tracer_turb) THEN
            name_tend   = "ddt_"//TRIM(comin_request_item%descriptor%name)//"_turb"

            target_var      => find_list_element(prm_nwp_tend_list(jg), 'ddt_tracer_turb')
            target_info     => target_var%info
            tracer_tend_idx =  target_info%ncontained+1

            CALL add_ref( prm_nwp_tend_list(jg), 'ddt_tracer_turb', name_tend,      &
              &           prm_nwp_tend(jg)%tracer_turb_ptr(tracer_tend_idx)%p_3d,   &
              &           target_info%hgrid, target_info%vgrid,                     &
              &           cf_desc, grib2_desc, ref_idx=tracer_tend_idx,             &
              &           ldims=shape3d_c, lrestart=.FALSE.,                        &
              &           in_group=groups("comin_tendencies") )

            this_info => comin_config%comin_icon_domain_config(jg)%tracer_info_head
            TRACER_INFO_LOOP_TURB: DO WHILE (ASSOCIATED(this_info))
              IF ( TRIM(this_info%name) == TRIM(comin_request_item%descriptor%name) ) THEN
                this_info%idx_turb = tracer_tend_idx
                EXIT TRACER_INFO_LOOP_TURB
              END IF
              this_info => this_info%next
            END DO TRACER_INFO_LOOP_TURB

            IF (msg_level >= 10) THEN
              WRITE (message_text,'(A,A,A,A,A,I2)') "Tendency due to turbulence '",TRIM(name_tend), &
                   &                                   "' added for tracer '", TRIM(comin_request_item%descriptor%name), &
                   &                                   "' in domain ",jg
              CALL message(routine, message_text)
            ENDIF

          END IF !comin_request_item%metadata%tracer_turb

          ! Adding tendencies of tracers due to convection
          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer_conv", tracer_conv, .FALSE.)
          IF (tracer_conv) THEN
            name_tend   = "ddt_"//TRIM(comin_request_item%descriptor%name)//"_conv"

            target_var      => find_list_element(prm_nwp_tend_list(jg), 'ddt_tracer_pconv')
            target_info     => target_var%info
            tracer_tend_idx =  target_info%ncontained+1

            CALL add_ref( prm_nwp_tend_list(jg), 'ddt_tracer_pconv', name_tend,     &
              &           prm_nwp_tend(jg)%tracer_conv_ptr(tracer_tend_idx)%p_3d,   &
              &           target_info%hgrid, target_info%vgrid,                     &
              &           cf_desc, grib2_desc, ref_idx=tracer_tend_idx,             &
              &           ldims=shape3d_c, lrestart=.FALSE.,                        &
              &           in_group=groups("comin_tendencies") )

            this_info => comin_config%comin_icon_domain_config(jg)%tracer_info_head
            TRACER_INFO_LOOP_CONV: DO WHILE (ASSOCIATED(this_info))
              IF ( TRIM(this_info%name) == TRIM(comin_request_item%descriptor%name) ) THEN
                this_info%idx_conv = tracer_tend_idx
                EXIT TRACER_INFO_LOOP_CONV
              END IF
              this_info => this_info%next
            END DO TRACER_INFO_LOOP_CONV

            IF (msg_level >= 10) THEN
              WRITE (message_text,'(A,A,A,A,A,I2)') "Tendency due to convection '",TRIM(name_tend), &
                   &                                   "' added for tracer '", TRIM(comin_request_item%descriptor%name), &
                   &                                   "' in domain ",jg
              CALL message(routine, message_text)
            ENDIF

          END IF !comin_request_item%metadata%tracer_conv

        END ASSOCIATE
        ptr => ptr%next()
      END DO VAR_LOOP

      ! Consistency check for the number of turbulent tracer
      nturb_tracer = 0
      this_info => comin_config%comin_icon_domain_config(jg)%tracer_info_head
      DO WHILE (ASSOCIATED(this_info))
        IF (this_info%idx_turb > 0) nturb_tracer = nturb_tracer + 1
        this_info => this_info%next
      END DO
      IF (nturb_tracer /= comin_config%comin_icon_domain_config(jg)%nturb_tracer &
           & .AND. msg_level >= 10) THEN
        WRITE (message_text,'(A,A,I3,A,I3)') 'Inconsistent number of turb tracers, ', &
          &                                   'nturb_tracer: ',nturb_tracer,          &
          &                                   ', comin_config: ',                     &
          &                                   comin_config%comin_icon_domain_config(jg)%nturb_tracer
        CALL finish(routine, message_text)
      END IF
    END DO DOM_LOOP
    IF (timers_level > 2) CALL timer_stop(timer_comin_init)

  END SUBROUTINE icon_append_comin_tracer_phys_tend

  !> loop over the total list of additional requested variables and
  !  perform `add_var` operations needed.
  !
  ! TODOs:
  ! - flag "lrestart"
  ! - vertical axis
  ! - OpenACC handling
  ! - GRIB2
  !
  SUBROUTINE icon_append_comin_variables(p_patch)
    TYPE(t_patch), INTENT(IN)  :: p_patch(:)
    !
    CHARACTER(*), PARAMETER    :: routine = modname//"::icon_append_comin_variables"
    INTEGER                    :: jg, ist, shape2d_c(2), shape3d_c(3), shape3d_c_half(3), zaxis_id
    CHARACTER(LEN=2)           :: dom_str
    CLASS(t_var_request_list_item), POINTER :: ptr
    REAL(wp),          POINTER :: tmp2d(:,:), tmp3d(:,:,:)
    TYPE(t_cf_var)             :: cf_desc
    TYPE(t_grib2_var)          :: grib2_desc
    LOGICAL                    :: restart, tracer
    CHARACTER(LEN=:), ALLOCATABLE :: charbuff

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    ! ComIn variables are added to a separate variable list.
    ALLOCATE(p_comin_varlist(n_dom), STAT=ist)
    IF (ist /= SUCCESS)  CALL finish (routine, "allocation of ComIn variable list failed")

    CALL comin_var_set_sync_device_mem(icon_sync_variable)

    DOM_LOOP : DO jg=1,n_dom

      WRITE(dom_str, "(i2.2)") jg
      CALL vlr_add(p_comin_varlist(jg), 'comin__'//dom_str,        &
        & patch_id=jg, vlevel_type=level_type_ml, lrestart=.TRUE., &
        & model_type=get_my_process_name())

      shape2d_c      = [ nproma,                     p_patch(jg)%nblks_c ]
      shape3d_c      = [ nproma, p_patch(jg)%nlev,   p_patch(jg)%nblks_c ]
      shape3d_c_half = [ nproma, p_patch(jg)%nlev+1, p_patch(jg)%nblks_c ]

      ptr => comin_request_get_list_head()
      VAR_LOOP : DO WHILE (ASSOCIATED(ptr))
        ! Note: only single-time-level variables are added
        !
        ! since we are dealing with a generic linked list, we have to
        ! perform a "dynamic type-cast" into `t_comin_request_item`:
        ASSOCIATE (comin_request_item => ptr%item_value)
          ! skip if variable was not requested for this domain
          IF (comin_request_item%descriptor%id /= jg) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP
          END IF

          CALL comin_metadata_get_or(comin_request_item%metadata, "tracer", tracer, .FALSE.)
          IF (tracer) THEN
            ptr => ptr%next()
            CYCLE VAR_LOOP
          END IF

          CALL comin_metadata_get_or(comin_request_item%metadata, "restart", restart, .FALSE.)
          CALL comin_metadata_get_or(comin_request_item%metadata, "units", charbuff, "")
          cf_desc%units         = charbuff
          CALL comin_metadata_get_or(comin_request_item%metadata, "standard_name", charbuff, "")
          cf_desc%standard_name = charbuff
          CALL comin_metadata_get_or(comin_request_item%metadata, "long_name", charbuff, "")
          cf_desc%long_name     = charbuff
          CALL comin_metadata_get_or(comin_request_item%metadata, "short_name", charbuff, "")
          cf_desc%short_name    = charbuff

          CALL comin_metadata_get_or(comin_request_item%metadata, "zaxis_id", zaxis_id, COMIN_ZAXIS_3D)

          IF (zaxis_id == COMIN_ZAXIS_3D) THEN
            CALL add_var(p_comin_varlist(jg),                                       &
              &          comin_request_item%descriptor%name, tmp3d,           &
              &          GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
              &          ldims=shape3d_c, lrestart=restart, in_group=groups("comin_vars"),            &
              &          lopenacc = .TRUE. )
          ELSE IF(zaxis_id == COMIN_ZAXIS_3D_HALF) THEN
            CALL add_var(p_comin_varlist(jg),                                       &
              &          comin_request_item%descriptor%name, tmp3d,           &
              &          GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
              &          ldims=shape3d_c_half, lrestart=restart, in_group=groups("comin_vars"),            &
              &          lopenacc = .TRUE. )
          ELSE
            CALL add_var(p_comin_varlist(jg),                                       &
              &          comin_request_item%descriptor%name, tmp2d,           &
              &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
              &          ldims=shape2d_c, lrestart=restart, in_group=groups("comin_vars"),            &
              &          lopenacc = .TRUE. )
          END IF

          IF (msg_level >= 10) THEN
            WRITE (message_text,*) "add variable '", &
                 &  TRIM(comin_request_item%descriptor%name), &
                 &  "', requested by third party plugins."
            CALL message(routine, message_text)
          ENDIF
        END ASSOCIATE ! comin_request_item
        ptr => ptr%next()

      END DO VAR_LOOP
    END DO DOM_LOOP
    IF (timers_level > 2) CALL timer_stop(timer_comin_init)

  END SUBROUTINE icon_append_comin_variables

  !> Expose ICON's variables to the ComIn infrastructure.
  !
  SUBROUTINE icon_expose_variables()
    CHARACTER(*), PARAMETER          :: routine = modname//"::icon_expose_variables"
    CHARACTER(:), ALLOCATABLE        :: temp_name
    TYPE(t_vl_register_iter)         :: vl_iter
    LOGICAL                          :: is_2d_field,  multi_timelevel_logical
    INTEGER                          :: zaxis_id, comin_zaxis_id, pos_jcjkjb(3), &
      &                                 ierr, jg, iv, ic, itrac, comin_hgrid_id, &
      &                                 iref_pos
    TYPE(t_var),            POINTER  :: elem
    TYPE(t_comin_var_ptr),  POINTER  :: comin_var_ptr
    REAL(wp), POINTER                :: ptr(:,:,:,:,:)
    TYPE(t_comin_var_descriptor)     :: descriptor
    TYPE(t_advection_config),POINTER :: advconf
    TYPE(t_zaxisType)                :: zaxisType
    TYPE(t_exposed_timedep_vars_list), POINTER :: exposed_timedep_vars_temp => NULL()
    TYPE(t_exposed_timedep_vars_list), POINTER :: temp_ptr => NULL()

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    IF (msg_level >= 10) THEN
      WRITE (message_text,*) "expose ICON's variables to the ComIn infrastructure"
      CALL message(routine, message_text)
    ENDIF

    ! loop of variable lists
    DO WHILE(vl_iter%next())
      jg = vl_iter%cur%p%patch_id
      IF(vl_iter%cur%p%vlevel_type /= level_type_ml)  CYCLE
      ! loop over variables
      DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p

        ! register only double-precision floating-point variables
        IF (.NOT. ASSOCIATED(elem%r_ptr))  CYCLE

        ! for time-level dependent variables: register only once
        IF (get_var_timelevel(elem%info%name) > 1)  CYCLE

        zaxis_id = elem%info%vgrid

        IF (zaxis_id == ZA_REFERENCE) THEN
          comin_zaxis_id = COMIN_ZAXIS_3D
        ELSE IF (zaxis_id == ZA_REFERENCE_HALF .OR. zaxis_id == ZA_REFERENCE_HALF_HHL ) THEN
          comin_zaxis_id = COMIN_ZAXIS_3D_HALF
        ELSE IF (zaxis_id == ZA_SURFACE) THEN
          comin_zaxis_id = COMIN_ZAXIS_2D
        ELSE
          comin_zaxis_id = COMIN_ZAXIS_UNDEF
        END IF

        IF (elem%info%hgrid == GRID_UNSTRUCTURED_CELL) THEN
          comin_hgrid_id = COMIN_HGRID_UNSTRUCTURED_CELL
        ELSE IF (elem%info%hgrid == GRID_UNSTRUCTURED_EDGE) THEN
          comin_hgrid_id = COMIN_HGRID_UNSTRUCTURED_EDGE
        ELSE IF (elem%info%hgrid == GRID_UNSTRUCTURED_VERT) THEN
          comin_hgrid_id = COMIN_HGRID_UNSTRUCTURED_VERTEX
        ELSE
          CYCLE
        END IF

        ! Provide data pointer, together with index positions (for
        ! correct interpretation).
        zaxisType = zaxisTypeList%getEntry(zaxis_id, ierr)
        IF (ierr == 0) THEN
          is_2d_field = zaxisType%is_2d
        ELSE
          is_2d_field = .FALSE.
        END IF

        ! In ICON, arrays have the implicit ordering: `jc`, `jk`, `jb`.
        !
        ! - in the case of 2D arrays, `jk` is omitted
        ! - the dimensions start with position 1 (i.e. usually:
        !   dimension position 1 = `jc`, position 2 = `jk`, position 3 = `jb`)
        ! - tracer variables correspond to 4D slices of the 5D array,
        !   where the position of the slicing dimension, and the slice
        !   index are stored in the variable's meta-data (`info%var_ref_pos`)

        IF (is_2d_field) THEN
          pos_jcjkjb =  [1, -1, 2]
        ELSE
          pos_jcjkjb =  [1, 2, 3]
        END IF

        iref_pos = elem%info%var_ref_pos
        IF (elem%info%ncontained > 0) THEN
          ! set reference position for containers, default: ndim + 1
          IF (elem%info%var_ref_pos < 0) iref_pos = MERGE(3, 4, is_2d_field)
          IF (iref_pos <= pos_jcjkjb(1)) THEN
            pos_jcjkjb(1:3) = pos_jcjkjb(1:3) + 1
          ELSEIF (iref_pos <= pos_jcjkjb(2)) THEN
            pos_jcjkjb(2:3) = pos_jcjkjb(2:3) + 1
          ELSEIF (iref_pos <= pos_jcjkjb(3)) THEN
            pos_jcjkjb(3:3) = pos_jcjkjb(3:3) + 1
          END IF
        END IF

        ! ComIn expects only the slice of the tracer
        ptr => elem%r_ptr
        IF (.NOT. elem%info%lcontainer) THEN
          ic = elem%info%ncontained
          SELECT CASE(iref_pos)
          CASE (1)
            ptr => elem%r_ptr(ic:ic,:,:,:,:)
          CASE (2)
            ptr => elem%r_ptr(:,ic:ic,:,:,:)
          CASE (3)
            ptr => elem%r_ptr(:,:,ic:ic,:,:)
          CASE (4)
            ptr => elem%r_ptr(:,:,:,ic:ic,:)
          CASE (5)
            ptr => elem%r_ptr(:,:,:,:,ic:ic)
          END SELECT
        END IF
        ALLOCATE(comin_var_ptr)
        comin_var_ptr%ptr => ptr
#ifdef _OPENACC
        comin_var_ptr%device_ptr = acc_deviceptr(C_LOC(ptr))
#endif
        comin_var_ptr%pos_jc = pos_jcjkjb(1)
        comin_var_ptr%pos_jk = pos_jcjkjb(2)
        comin_var_ptr%pos_jb = pos_jcjkjb(3)
        comin_var_ptr%pos_jn = iref_pos
        comin_var_ptr%lcontainer = elem%info%lcontainer
        comin_var_ptr%ncontained=elem%info%ncontained


        comin_var_ptr%descriptor%name=get_var_name(elem%info)
        comin_var_ptr%descriptor%id=jg
        CALL comin_var_list_append(       &
          &   p          = comin_var_ptr)

        ! Copy all metadata set by plugins from request_list to var_list
        descriptor = comin_var_ptr%descriptor
        CALL icon_expose_metadata_from_request_list(descriptor)

        ! Set metadata based on ICON settings. If these are different from metadata requested by a plugin
        !   this overwrites plugin settings without warning
        advconf => advection_config(jg)
        CALL comin_metadata_set(descriptor, "restart", LOGICAL(elem%info%lrestart))
        CALL comin_metadata_set(descriptor, "tracer", elem%info_dyn%tracer%lis_tracer)
        IF (elem%info_dyn%tracer%lis_tracer) THEN
          CALL comin_metadata_set(descriptor, "tracer_turb", elem%info_dyn%tracer%lturb_tracer)
          CALL comin_metadata_set(descriptor, "tracer_conv", elem%info_dyn%tracer%lconv_tracer)
          CALL comin_metadata_set(descriptor, "tracer_hadv", elem%info_dyn%tracer%ihadv_tracer)
          CALL comin_metadata_set(descriptor, "tracer_vadv", elem%info_dyn%tracer%ivadv_tracer)
          DO itrac = 1, max_ntracer
             IF (TRIM(ADJUSTL(get_var_name(elem%info))) == TRIM(ADJUSTL(advconf%tracer_names(itrac)))) THEN
               CALL comin_metadata_set(descriptor, "tracer_hlimit", advconf%itype_hlimit(itrac))
               CALL comin_metadata_set(descriptor, "tracer_vlimit", advconf%itype_vlimit(itrac))
             END IF
          END DO
        END IF

        CALL comin_metadata_set(descriptor, "zaxis_id",      comin_zaxis_id)
        CALL comin_metadata_set(descriptor, "hgrid_id",      comin_hgrid_id)
        CALL comin_metadata_set(descriptor, "units",         elem%info%cf%units)
        CALL comin_metadata_set(descriptor, "standard_name", elem%info%cf%standard_name)
        CALL comin_metadata_set(descriptor, "long_name",     elem%info%cf%long_name)
        CALL comin_metadata_set(descriptor, "short_name",    elem%info%cf%short_name)
        IF (get_var_timelevel(elem%info%name) == -1) THEN
           multi_timelevel_logical = .FALSE.
        ELSE
           multi_timelevel_logical = .TRUE.
        END IF
        CALL comin_metadata_set(descriptor, "multi_timelevel", multi_timelevel_logical)

        !here, a linked list is filled, which consists of icon variables that are exposed to comin and are time-dependent
        ALLOCATE(exposed_timedep_vars_temp)
        !for time-dependent variables (only time level 1 is considered), add to the beginning of the linked list
        IF (get_var_timelevel(elem%info%name) /= -1) THEN
           exposed_timedep_vars_temp%icon_var => elem
           exposed_timedep_vars_temp%comin_var_ptr => comin_var_ptr
           exposed_timedep_vars_temp%next => exposed_timedep_vars_head
           exposed_timedep_vars_head => exposed_timedep_vars_temp
        END IF
     END DO
    END DO
   !fill the linked list for the time level 2
    DO WHILE(vl_iter%next())
       jg = vl_iter%cur%p%patch_id
       IF(vl_iter%cur%p%vlevel_type /= level_type_ml)  CYCLE
       ! loop over variables
       DO iv = 1, vl_iter%cur%p%nvars
          elem => vl_iter%cur%p%vl(iv)%p
          ! register only double-precision floating-point variables
          IF (.NOT. ASSOCIATED(elem%r_ptr))  CYCLE
          !consider only variables for time level 2
          IF (get_var_timelevel(elem%info%name) <= 1)  CYCLE
          ALLOCATE(exposed_timedep_vars_temp)
          exposed_timedep_vars_temp%icon_var => elem
          temp_ptr => exposed_timedep_vars_head
          temp_name = TRIM(ADJUSTL(get_var_name(elem%info)))
          search_loop: DO WHILE(ASSOCIATED(temp_ptr))
             IF (temp_ptr%comin_var_ptr%descriptor%id==jg .and. temp_name == TRIM(ADJUSTL(temp_ptr%comin_var_ptr%descriptor%name))) THEN
                exposed_timedep_vars_temp%comin_var_ptr => temp_ptr%comin_var_ptr
                EXIT search_loop
             END IF
             temp_ptr => temp_ptr%next
          END DO search_loop
          exposed_timedep_vars_temp%next => exposed_timedep_vars_head
          exposed_timedep_vars_head => exposed_timedep_vars_temp
       END DO
    END DO
    IF (timers_level > 2) CALL timer_stop(timer_comin_init)
  END SUBROUTINE icon_expose_variables

  !> Expose all metadata coming from the plugins via request list
  !    to plugins via var list.
  !    Some of those metadata may be overwritten in the following
  !    by icon_expose_variables
  SUBROUTINE icon_expose_metadata_from_request_list(descriptor)
    TYPE(t_comin_var_descriptor)           :: descriptor
    ! Local variables
    TYPE(t_var_request_list_item), POINTER :: ptr
    TYPE(t_comin_var_metadata_iterator)    :: metadata_it
    CHARACTER(LEN=:), ALLOCATABLE          :: current_metadata_key, charbuff
    LOGICAL        :: logbuff
    INTEGER        :: intbuff
    REAL(C_DOUBLE) :: realbuff


    ptr => comin_request_get_list_head()
    VAR_LOOP : DO WHILE (ASSOCIATED(ptr))
      ASSOCIATE (comin_request_item => ptr%item_value)
        IF (TRIM(descriptor%name) == TRIM(comin_request_item%descriptor%name) &
          & .AND. descriptor%id == comin_request_item%descriptor%id) THEN

          ! Iterate through request list metadata, forward metadata to var list
          CALL ptr%item_value%metadata%get_iterator(metadata_it)
          DO WHILE(.NOT. metadata_it%is_end())
            current_metadata_key = metadata_it%key()
            SELECT CASE(ptr%item_value%metadata%query(current_metadata_key))
            CASE (COMIN_METADATA_TYPEID_INTEGER)
              CALL ptr%item_value%metadata%get(current_metadata_key, intbuff)
              CALL comin_metadata_set(descriptor, current_metadata_key, intbuff)
            CASE (COMIN_METADATA_TYPEID_REAL)
              CALL ptr%item_value%metadata%get(current_metadata_key, realbuff)
              CALL comin_metadata_set(descriptor, current_metadata_key, realbuff)
            CASE (COMIN_METADATA_TYPEID_CHARACTER)
              CALL ptr%item_value%metadata%get(current_metadata_key, charbuff)
              CALL comin_metadata_set(descriptor, current_metadata_key, charbuff)
            CASE (COMIN_METADATA_TYPEID_LOGICAL)
              CALL ptr%item_value%metadata%get(current_metadata_key, logbuff)
              CALL comin_metadata_set(descriptor, current_metadata_key, logbuff)
            END SELECT
            CALL metadata_it%next()
          ENDDO
          CALL metadata_it%delete()
        ENDIF ! descriptor names and ids
      END ASSOCIATE
      ptr => ptr%next()

    END DO VAR_LOOP

  END SUBROUTINE icon_expose_metadata_from_request_list

  !> fill ComIn descriptive data structures from ICON
  !
  SUBROUTINE icon_expose_descrdata(msg_level, n_dom, max_dom, nproma, min_rlcell_int,    &
    &                              min_rlcell, max_rlcell, min_rlvert_int, min_rlvert,   &
    &                              max_rlvert, min_rledge_int, min_rledge, max_rledge,   &
    &                              grf_bdywidth_c, grf_bdywidth_e, lrestart,             &
    &                              vct_a, p_patch, number_of_grid_used, vgrid_buffer,    &
    &                              start_time, end_time, time_config, dtime)
    INTEGER,              INTENT(IN) :: msg_level, n_dom, max_dom, nproma
    INTEGER,              INTENT(IN) :: grf_bdywidth_c, grf_bdywidth_e
    INTEGER,              INTENT(IN) :: min_rlcell_int, min_rlcell, max_rlcell
    INTEGER,              INTENT(IN) :: min_rlvert_int, min_rlvert, max_rlvert
    INTEGER,              INTENT(IN) :: min_rledge_int, min_rledge, max_rledge
    LOGICAL,              INTENT(IN) :: lrestart
    REAL(wp),             INTENT(IN) :: vct_a(:), start_time(:), end_time(:)
    INTEGER,              INTENT(IN) :: number_of_grid_used(:)
    TYPE(t_patch),        INTENT(IN), TARGET :: p_patch(:)
    TYPE(t_vgrid_buffer), INTENT(IN), TARGET :: vgrid_buffer(:)
    TYPE(t_time_config),  INTENT(IN) :: time_config
    REAL(wp),             INTENT(IN) :: dtime
    !
    CHARACTER(*), PARAMETER :: routine = modname//"::icon_expose_descrdata"
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: sim_start, sim_end, sim_current, &
      &                                    run_start, run_stop

    CALL comin_setup_set_verbosity_level(msg_level)
    ! expose descriptive data structures
    CALL icon_expose_descrdata_global(n_dom, max_dom, nproma, min_rlcell_int, min_rlcell, &
      &                               max_rlcell, min_rlvert_int, min_rlvert, max_rlvert, &
      &                               min_rledge_int, min_rledge, max_rledge,             &
      &                               grf_bdywidth_c, grf_bdywidth_e, lrestart, vct_a)
    ! p_patch is used from index 1 to exclude potential coarse radiation grid
    CALL icon_expose_descrdata_domain(p_patch(1:), number_of_grid_used, vgrid_buffer, &
         &                     start_time, end_time)
    CALL datetimeToString(time_config%tc_exp_startdate, sim_start)
    CALL datetimeToString(time_config%tc_exp_stopdate, sim_end)
    CALL datetimeToString(time_config%tc_exp_startdate, sim_current)
    CALL datetimeToString(time_config%tc_startdate, run_start)
    CALL datetimeToString(time_config%tc_stopdate, run_stop)
    CALL icon_expose_descrdata_state(sim_start, sim_end, sim_current, run_start, run_stop)
    CALL icon_expose_timesteplength_domain(1, dtime)
    ! - call primary constructors
    IF (timers_level > 2) CALL timer_start(timer_comin_primary_constructors)
    CALL comin_plugin_primaryconstructor(comin_config%plugin_list(1:comin_config%nplugins))
    IF (timers_level > 2) CALL timer_stop(timer_comin_primary_constructors)
  END SUBROUTINE icon_expose_descrdata


  !> fill ComIn descriptive data structures from ICON: global information
  !
  SUBROUTINE icon_expose_descrdata_global(n_dom, max_dom, nproma, min_rlcell_int, min_rlcell, &
    &                                     max_rlcell, min_rlvert_int, min_rlvert, max_rlvert, &
    &                                     min_rledge_int, min_rledge, max_rledge,             &
    &                                     grf_bdywidth_c, grf_bdywidth_e, lrestart, vct_a)
    INTEGER,  INTENT(IN) :: n_dom
    INTEGER,  INTENT(IN) :: max_dom
    INTEGER,  INTENT(IN) :: nproma
    INTEGER,  INTENT(IN) :: min_rlcell_int
    INTEGER,  INTENT(IN) :: min_rlcell
    INTEGER,  INTENT(IN) :: max_rlcell
    INTEGER,  INTENT(IN) :: min_rlvert_int
    INTEGER,  INTENT(IN) :: min_rlvert
    INTEGER,  INTENT(IN) :: max_rlvert
    INTEGER,  INTENT(IN) :: min_rledge_int
    INTEGER,  INTENT(IN) :: min_rledge
    INTEGER,  INTENT(IN) :: max_rledge

    INTEGER,  INTENT(IN) :: grf_bdywidth_c
    INTEGER,  INTENT(IN) :: grf_bdywidth_e
    LOGICAL,  INTENT(IN) :: lrestart
    REAL(wp), INTENT(IN) :: vct_a(:)
    ! local variables
    TYPE(t_comin_descrdata_global) :: comin_descrdata_global_data
    INTEGER :: host_url_length, host_branch_length, host_revision_length, host_tag_length
#ifdef _OPENACC
    INTEGER(kind=acc_device_kind)  :: devicetype
    INTEGER                        :: devicenum, nullidx
    CHARACTER(64)                  :: buffer_str
#endif

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    comin_descrdata_global_data%n_dom            = n_dom
    comin_descrdata_global_data%max_dom          = max_dom
    comin_descrdata_global_data%nproma           = nproma
    comin_descrdata_global_data%min_rlcell_int   = min_rlcell_int
    comin_descrdata_global_data%min_rlcell       = min_rlcell
    comin_descrdata_global_data%max_rlcell       = max_rlcell
    comin_descrdata_global_data%min_rlvert_int   = min_rlvert_int
    comin_descrdata_global_data%min_rlvert       = min_rlvert
    comin_descrdata_global_data%max_rlvert       = max_rlvert
    comin_descrdata_global_data%min_rledge_int   = min_rledge_int
    comin_descrdata_global_data%min_rledge       = min_rledge
    comin_descrdata_global_data%max_rledge       = max_rledge
    comin_descrdata_global_data%grf_bdywidth_c   = grf_bdywidth_c
    comin_descrdata_global_data%grf_bdywidth_e   = grf_bdywidth_e
    comin_descrdata_global_data%lrestartrun      = lrestart
    ALLOCATE(comin_descrdata_global_data%vct_a(SIZE(vct_a)))
    comin_descrdata_global_data%vct_a(:)         = vct_a(:)

    comin_descrdata_global_data%yac_instance_id  = -1
#ifdef YAC_coupling
    IF (is_coupled_run()) THEN
      comin_descrdata_global_data%yac_instance_id  = cpl_get_instance_id()
    END IF
#endif


    comin_descrdata_global_data%host_git_remote_url = get_remote_url('icon')
    comin_descrdata_global_data%host_git_branch = get_local_branch('icon')
    comin_descrdata_global_data%host_git_tag = get_icon_version()
    comin_descrdata_global_data%host_revision = get_revision('icon')

    ! extract information from openacc
    ! NVHPC's openacc impl. fills the entire buffer with \x00
#if _OPENACC
    devicetype = acc_get_device_type()
    comin_descrdata_global_data%has_device = (devicetype /= acc_device_none)
    IF (comin_descrdata_global_data%has_device) THEN
       devicenum = acc_get_device_num(devicetype)

       buffer_str(:) = " "
       CALL acc_get_property_string( devicenum, devicetype, &
            acc_property_name, buffer_str )
       nullidx = INDEX(buffer_str, CHAR(0))
       IF (nullidx == 0) nullidx = LEN_TRIM(buffer_str)
       comin_descrdata_global_data%device_name = buffer_str(1:nullidx-1)

       buffer_str(:) = " "
       CALL acc_get_property_string( devicenum, devicetype, &
            acc_property_vendor, buffer_str )
       nullidx = INDEX(buffer_str, CHAR(0))
       IF (nullidx == 0) nullidx = LEN_TRIM(buffer_str)
       comin_descrdata_global_data%device_vendor = buffer_str(1:nullidx-1)

       buffer_str(:) = " "
       CALL acc_get_property_string( devicenum, devicetype, &
            acc_property_driver, buffer_str )
       nullidx = INDEX(buffer_str, CHAR(0))
       IF (nullidx == 0) nullidx = LEN_TRIM(buffer_str)
       comin_descrdata_global_data%device_driver = buffer_str(1:nullidx-1)
    ENDIF
#else
    comin_descrdata_global_data%has_device = .FALSE.
#endif

    ! register global info in ComIn
    CALL comin_descrdata_set_global(comin_descrdata_global_data)
    IF (timers_level > 2) CALL timer_stop(timer_comin_init)

  END SUBROUTINE icon_expose_descrdata_global

  !> fill ComIn descriptive data structures from ICON: domain information
  !
  SUBROUTINE icon_expose_descrdata_domain(patch, number_of_grid_used, vgrid_buffer, start_time, end_time)
    CHARACTER(*), PARAMETER                  :: routine = modname//"::icon_expose_descrdata_domain"
    TYPE(t_patch), TARGET,        INTENT(IN) :: patch(:)
    INTEGER,                      INTENT(IN) :: number_of_grid_used(:)
    TYPE(t_vgrid_buffer), TARGET, INTENT(IN) :: vgrid_buffer(:)
    REAL(wp),                     INTENT(IN) :: start_time(:), end_time(:)
    !local var
    INTEGER :: jg
    TYPE(t_comin_descrdata_domain) :: comin_descrdata_domain(SIZE(patch))

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    DO jg=1,SIZE(patch)
      ! cell, vertex and edge independent variables
      comin_descrdata_domain(jg)%grid_filename => patch(jg)%grid_filename
      comin_descrdata_domain(jg)%grid_uuid => patch(jg)%grid_uuid%DATA
      ALLOCATE(comin_descrdata_domain(jg)%number_of_grid_used(SIZE(number_of_grid_used)))
      ! Note: number_of_grid_used is not part of a pointer structure
      comin_descrdata_domain(jg)%number_of_grid_used = number_of_grid_used
      comin_descrdata_domain(jg)%id => patch(jg)%id
      comin_descrdata_domain(jg)%n_childdom => patch(jg)%n_childdom
      comin_descrdata_domain(jg)%dom_start = start_time(jg)
      comin_descrdata_domain(jg)%dom_end = end_time(jg)
      comin_descrdata_domain(jg)%nlev => patch(jg)%nlev
      comin_descrdata_domain(jg)%nshift => patch(jg)%nshift
      comin_descrdata_domain(jg)%nshift_total => patch(jg)%nshift_total

      ! cell components
      comin_descrdata_domain(jg)%cells%ncells => patch(jg)%n_patch_cells
      comin_descrdata_domain(jg)%cells%ncells_global => patch(jg)%n_patch_cells_g
      comin_descrdata_domain(jg)%cells%nblks => patch(jg)%nblks_c
      comin_descrdata_domain(jg)%cells%refin_ctrl => patch(jg)%cells%refin_ctrl
      comin_descrdata_domain(jg)%cells%max_connectivity => patch(jg)%cells%max_connectivity
      comin_descrdata_domain(jg)%cells%num_edges => patch(jg)%cells%num_edges
      comin_descrdata_domain(jg)%cells%start_index => patch(jg)%cells%start_index
      comin_descrdata_domain(jg)%cells%end_index => patch(jg)%cells%end_index
      comin_descrdata_domain(jg)%cells%start_block => patch(jg)%cells%start_block
      comin_descrdata_domain(jg)%cells%end_block => patch(jg)%cells%end_block
      comin_descrdata_domain(jg)%cells%child_id => patch(jg)%cells%child_id
      comin_descrdata_domain(jg)%cells%child_idx => patch(jg)%cells%child_idx
      comin_descrdata_domain(jg)%cells%child_blk => patch(jg)%cells%child_blk
      comin_descrdata_domain(jg)%cells%parent_glb_idx => patch(jg)%cells%parent_glb_idx
      comin_descrdata_domain(jg)%cells%parent_glb_blk => patch(jg)%cells%parent_glb_blk
      comin_descrdata_domain(jg)%cells%vertex_idx => patch(jg)%cells%vertex_idx
      comin_descrdata_domain(jg)%cells%vertex_blk => patch(jg)%cells%vertex_blk
      comin_descrdata_domain(jg)%cells%neighbor_idx => patch(jg)%cells%neighbor_idx
      comin_descrdata_domain(jg)%cells%neighbor_blk => patch(jg)%cells%neighbor_blk
      comin_descrdata_domain(jg)%cells%edge_idx => patch(jg)%cells%edge_idx
      comin_descrdata_domain(jg)%cells%edge_blk => patch(jg)%cells%edge_blk
      ! Note: no pointer access to lon and lat since their structure is changed
      comin_descrdata_domain(jg)%cells%clon = patch(jg)%cells%center%lon
      comin_descrdata_domain(jg)%cells%clat = patch(jg)%cells%center%lat
      comin_descrdata_domain(jg)%cells%area => patch(jg)%cells%area
      ! Note: at the primary constructor p_nh_state(jg)%metrics%z_ifc is not filled yet
      ! vgrid_buffer (from which p_nh_state is filled at mo_vertical_grid::set_nh_metrics)
      ! needs to be copied here since it is deallocated after set_nh_metrics
      ALLOCATE(comin_descrdata_domain(jg)%cells%hhl(nproma,patch(jg)%nlevp1,patch(jg)%nblks_c))
      comin_descrdata_domain(jg)%cells%hhl = vgrid_buffer(jg)%z_ifc

      comin_descrdata_domain(jg)%cells%glb_index => patch(jg)%cells%decomp_info%glb_index
      comin_descrdata_domain(jg)%cells%decomp_domain => patch(jg)%cells%decomp_info%decomp_domain

      ! vertex components
      comin_descrdata_domain(jg)%verts%nverts => patch(jg)%n_patch_verts
      comin_descrdata_domain(jg)%verts%nverts_global => patch(jg)%n_patch_verts_g
      comin_descrdata_domain(jg)%verts%nblks => patch(jg)%nblks_v
      comin_descrdata_domain(jg)%verts%refin_ctrl => patch(jg)%verts%refin_ctrl
      comin_descrdata_domain(jg)%verts%start_index => patch(jg)%verts%start_index
      comin_descrdata_domain(jg)%verts%end_index => patch(jg)%verts%end_index
      comin_descrdata_domain(jg)%verts%start_block => patch(jg)%verts%start_block
      comin_descrdata_domain(jg)%verts%end_block => patch(jg)%verts%end_block
      comin_descrdata_domain(jg)%verts%neighbor_blk => patch(jg)%verts%neighbor_blk
      comin_descrdata_domain(jg)%verts%neighbor_idx => patch(jg)%verts%neighbor_idx
      comin_descrdata_domain(jg)%verts%cell_idx => patch(jg)%verts%cell_idx
      comin_descrdata_domain(jg)%verts%cell_blk => patch(jg)%verts%cell_blk
      comin_descrdata_domain(jg)%verts%edge_idx => patch(jg)%verts%edge_idx
      comin_descrdata_domain(jg)%verts%edge_blk => patch(jg)%verts%edge_blk
      comin_descrdata_domain(jg)%verts%num_edges => patch(jg)%verts%num_edges
      ! Note: no pointer access to lon and lat since their structure is changed
      comin_descrdata_domain(jg)%verts%vlon = patch(jg)%verts%vertex%lon
      comin_descrdata_domain(jg)%verts%vlat = patch(jg)%verts%vertex%lat

      ! edge components
      comin_descrdata_domain(jg)%edges%nedges => patch(jg)%n_patch_edges
      comin_descrdata_domain(jg)%edges%nedges_global => patch(jg)%n_patch_edges_g
      comin_descrdata_domain(jg)%edges%nblks => patch(jg)%nblks_e
      comin_descrdata_domain(jg)%edges%refin_ctrl => patch(jg)%edges%refin_ctrl
      comin_descrdata_domain(jg)%edges%start_index => patch(jg)%edges%start_index
      comin_descrdata_domain(jg)%edges%end_index => patch(jg)%edges%end_index
      comin_descrdata_domain(jg)%edges%start_block => patch(jg)%edges%start_block
      comin_descrdata_domain(jg)%edges%end_block => patch(jg)%edges%end_block
      comin_descrdata_domain(jg)%edges%child_id => patch(jg)%edges%child_id
      comin_descrdata_domain(jg)%edges%child_idx => patch(jg)%edges%child_idx
      comin_descrdata_domain(jg)%edges%child_blk => patch(jg)%edges%child_blk
      comin_descrdata_domain(jg)%edges%parent_glb_idx => patch(jg)%edges%parent_glb_idx
      comin_descrdata_domain(jg)%edges%parent_glb_blk => patch(jg)%edges%parent_glb_blk
      comin_descrdata_domain(jg)%edges%cell_idx => patch(jg)%edges%cell_idx
      comin_descrdata_domain(jg)%edges%cell_blk => patch(jg)%edges%cell_blk
      comin_descrdata_domain(jg)%edges%vertex_idx => patch(jg)%edges%vertex_idx
      comin_descrdata_domain(jg)%edges%vertex_blk => patch(jg)%edges%vertex_blk
      ! Note: no pointer access to lon and lat since their structure is changed
      comin_descrdata_domain(jg)%edges%elon = patch(jg)%edges%center%lon
      comin_descrdata_domain(jg)%edges%elat = patch(jg)%edges%center%lat
    END DO

    ! register domain info in ComIn
    CALL comin_descrdata_set_domain(comin_descrdata_domain)

    CALL comin_descrdata_set_fct_glb2loc_cell(glb2loc_index)

    IF (timers_level > 2) CALL timer_stop(timer_comin_init)

  END SUBROUTINE icon_expose_descrdata_domain

  !> fill ComIn descriptive data structure from ICON: simulation status
  !
  SUBROUTINE icon_expose_descrdata_state(sim_time_start, sim_time_end, sim_time_current, &
    &                                    run_time_start, run_time_stop)
    CHARACTER(LEN=*), INTENT(IN) :: sim_time_start
    CHARACTER(LEN=*), INTENT(IN) :: sim_time_end
    CHARACTER(LEN=*), INTENT(IN) :: sim_time_current
    CHARACTER(LEN=*), INTENT(IN) :: run_time_start
    CHARACTER(LEN=*), INTENT(IN) :: run_time_stop
    !local variables
    TYPE(t_comin_descrdata_simulation_interval) :: comin_descrdata_state

    IF (timers_level > 2) CALL timer_start(timer_comin_init)
    comin_descrdata_state%exp_start = sim_time_start
    comin_descrdata_state%exp_stop  = sim_time_end
    comin_descrdata_state%run_start = run_time_start
    comin_descrdata_state%run_stop  = run_time_stop

    ! register time info in ComIn
    CALL comin_descrdata_set_simulation_interval(comin_descrdata_state)
    CALL comin_current_set_datetime(sim_time_current)
    IF (timers_level > 2) CALL timer_stop(timer_comin_init)

  END SUBROUTINE icon_expose_descrdata_state

  ! Expose ICON's physics/advection time step of each domain
  !
  RECURSIVE SUBROUTINE icon_expose_timesteplength_domain(jg, dt_current)
    INTEGER,  INTENT(IN) :: jg
    REAL(wp), INTENT(IN) :: dt_current
    !local variables
    INTEGER :: jn

    CALL comin_descrdata_set_timesteplength(jg, dt_current)

    DO jn=1,p_patch(jg)%n_childdom
       CALL icon_expose_timesteplength_domain(p_patch(jg)%child_id(jn), dt_current/2.0_wp)
    END DO

  END SUBROUTINE icon_expose_timesteplength_domain

  INTEGER FUNCTION glb2loc_index(jg, glb) RESULT(loc)
    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(IN) :: glb

    loc = get_local_index(p_patch(jg)%cells%decomp_info%glb2loc_index,glb)
  END FUNCTION glb2loc_index

  !> Update ComIn descriptive data structure from ICON: state, sim_current
  !
  SUBROUTINE icon_update_current_datetime(sim_time_current)
    CHARACTER(LEN=*), INTENT(IN) :: sim_time_current

    ! register time info in ComIn
    CALL comin_current_set_datetime(sim_time_current)

  END SUBROUTINE icon_update_current_datetime

  !> Update pointers to current timelevel of exposed ICON variables.
  !
  SUBROUTINE icon_update_expose_variables(tlev_source, tlev)
    INTEGER, INTENT(IN)     :: tlev_source, tlev
    CHARACTER(*), PARAMETER :: routine = modname//"::icon_update_expose_variables"
    REAL(wp), POINTER       :: ptr(:,:,:,:,:)
    INTEGER                 :: ic, iref_pos
    TYPE(t_exposed_timedep_vars_list), POINTER  :: current => NULL()

    IF (msg_level >= 10) THEN
      WRITE (message_text,*) "Update pointers to current timelevel of exposed ICON variables."
      CALL message(routine, message_text)
    ENDIF
    current => exposed_timedep_vars_head
    DO WHILE (ASSOCIATED(current))
       IF (.NOT. ASSOCIATED(current%icon_var%r_ptr)) THEN
            current => current%next
            CYCLE
       END IF
       IF (current%icon_var%info%tlev_source /= tlev_source) THEN   ! select between TLEV_NNOW and TLEV_NNOW_RCF
          current => current%next
          CYCLE
       END IF
       IF (get_var_timelevel(current%icon_var%info%name) /= tlev) THEN !select the current timelevel to update comin pointer
          current => current%next
          CYCLE
       END IF
       ptr => current%icon_var%r_ptr
       iref_pos = current%comin_var_ptr%pos_jn
       
       IF (.NOT. current%icon_var%info%lcontainer) THEN
          ic = current%icon_var%info%ncontained
          SELECT CASE(iref_pos)
          CASE (1)
             ptr => current%icon_var%r_ptr(ic:ic,:,:,:,:)
          CASE (2)
             ptr => current%icon_var%r_ptr(:,ic:ic,:,:,:)
          CASE (3)
             ptr => current%icon_var%r_ptr(:,:,ic:ic,:,:)
          CASE (4)
             ptr => current%icon_var%r_ptr(:,:,:,ic:ic,:)
          CASE (5)
             ptr => current%icon_var%r_ptr(:,:,:,:,ic:ic)
          END SELECT
       END IF   
       current%comin_var_ptr%ptr => ptr
#ifdef _OPENACC
       current%comin_var_ptr%device_ptr = acc_deviceptr(C_LOC(ptr))
#endif

       current => current%next

    END DO
  END SUBROUTINE icon_update_expose_variables

  SUBROUTINE icon_call_callback(ep, jg, lacc)
    INTEGER, INTENT(IN) :: ep
    INTEGER, INTENT(IN) :: jg
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

    IF (timers_level > 2) CALL timer_start(timer_comin_callbacks)
    CALL comin_callback_context_call(ep, jg, lzacc)
    IF (timers_level > 2) CALL timer_stop(timer_comin_callbacks)
  END SUBROUTINE icon_call_callback

  SUBROUTINE icon_sync_variable(var_ptr, direction)
    TYPE(t_comin_var_ptr), POINTER, INTENT(IN) :: var_ptr
    LOGICAL, INTENT(IN) :: direction

    IF (direction) THEN ! to device
       !$ACC UPDATE DEVICE(var_ptr%ptr)
    ELSE
       !$ACC UPDATE HOST(var_ptr%ptr)
    END IF
  END SUBROUTINE icon_sync_variable
#endif
END MODULE mo_comin_adapter
