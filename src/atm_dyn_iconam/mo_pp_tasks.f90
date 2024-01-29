!
! Tasks for internal post-processing.
!
! The subroutines in this module can be inserted into a dynamic "job queue".
! See module "mo_pp_scheduler" for detailed info.
!
! -----------------------------------------------------------------------------------
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

MODULE mo_pp_tasks

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message, finish, warning
  USE mo_impl_constants,          ONLY: SUCCESS,                      &
    & VINTP_METHOD_VN, VINTP_METHOD_LIN, VINTP_METHOD_QV,             &
    & VINTP_METHOD_LIN_NLEVP1,                                        &
    & TASK_INIT_VER_Z, TASK_INIT_VER_P, TASK_INIT_VER_I,              &
    & TASK_FINALIZE_IPZ,                                              &
    & TASK_INTP_HOR_LONLAT, TASK_INTP_VER_PLEV,                       &
    & TASK_COMPUTE_RH, TASK_COMPUTE_PV, TASK_COMPUTE_SMI,             &
    & TASK_COMPUTE_SDI2, TASK_COMPUTE_LPI, TASK_COMPUTE_CEILING,      &
    & TASK_COMPUTE_HBAS_SC, TASK_COMPUTE_HTOP_SC,                     &
    & TASK_COMPUTE_INVERSION,                                         &
    & TASK_COMPUTE_TWATER, TASK_COMPUTE_Q_SEDIM,                      &
    & TASK_COMPUTE_DBZCMAX, TASK_COMPUTE_DBZ850,                      &
    & TASK_COMPUTE_DBZLMX_LOW, TASK_COMPUTE_SRH, TASK_COMPUTE_VIS,    &
    & TASK_COMPUTE_WSHEAR_U, TASK_COMPUTE_WSHEAR_V,                   &
    & TASK_COMPUTE_LAPSERATE,                                         &
    & TASK_INTP_VER_ZLEV,                                             &
    & TASK_INTP_VER_ILEV,                                             &
    & PRES_MSL_METHOD_SAI, PRES_MSL_METHOD_GME, max_dom,              &
    & ALL_TIMELEVELS, PRES_MSL_METHOD_IFS, PRES_MSL_METHOD_DWD,       &
    & PRES_MSL_METHOD_IFS_CORR, RH_METHOD_WMO, RH_METHOD_IFS,         &
    & RH_METHOD_IFS_CLIP, TASK_COMPUTE_OMEGA, HINTP_TYPE_LONLAT_BCTR, &
    & TLEV_NNOW, TLEV_NNOW_RCF, HINTP_TYPE_LONLAT_RBF
  USE mo_model_domain,            ONLY: t_patch, p_patch_local_parent
  USE mo_var,                     ONLY: t_var
  USE mo_var_metadata_types,      ONLY: t_var_metadata, t_vert_interp_meta
  USE mo_intp,                    ONLY: cell_avg, cells2edges_scalar
  USE mo_intp_data_strc,          ONLY: t_int_state, p_int_state,     &
    &                                   p_int_state_local_parent
  USE mo_intp_lonlat_types,       ONLY: t_lon_lat_intp, lonlat_grids
  USE mo_intp_rbf,                ONLY: rbf_vec_interpol_cell
  USE mo_nh_vert_interp,          ONLY: lin_intp, uv_intp, qv_intp,         &
    &                                   prepare_extrap, prepare_extrap_ifspp
  USE mo_nh_vert_interp_ipz,      ONLY: prepare_vert_interp_z,              &
    &                                   prepare_vert_interp_p,              &
    &                                   prepare_vert_interp_i
  USE mo_nh_diagnose_pmsl,        ONLY: diagnose_pmsl, diagnose_pmsl_gme,   &
    &                                   diagnose_pmsl_ifs
  USE mo_nonhydro_types,          ONLY: t_nh_state, t_nh_prog, t_nh_diag,   &
    &                                   t_nh_metrics
  USE mo_opt_diagnostics,         ONLY: t_nh_diag_pz, t_nh_opt_diag, t_vcoeff, &
    &                                   vcoeff_deallocate, t_vcoeff_lin, t_vcoeff_cub
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_nh_pzlev_config,         ONLY: t_nh_pzlev_config
  USE mo_parallel_config,         ONLY: nproma
  USE mo_dynamics_config,         ONLY: nnow, nnow_rcf
  USE mo_zaxis_type,              ONLY: zaxisTypeList
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_CELL,                  &
    &                                   GRID_UNSTRUCTURED_EDGE
  USE mo_sync,                    ONLY: sync_patch_array,                        &
    &                                   SYNC_C, SYNC_E,                          &
    &                                   cumulative_sync_patch_array,             &
    &                                   complete_cumulative_sync
  USE mo_util_phys,               ONLY: compute_field_rel_hum_wmo,               &
    &                                   compute_field_rel_hum_ifs
  USE mo_opt_nwp_diagnostics,     ONLY: compute_field_omega,                     &
    &                                   compute_field_pv,                        &
    &                                   compute_field_sdi,                       &
    &                                   compute_field_lpi,                       &
    &                                   compute_field_ceiling,                   &
    &                                   compute_field_visibility,                &
    &                                   compute_field_hbas_sc, compute_field_htop_sc, &
    &                                   compute_field_twater, compute_field_q_sedim,  &
    &                                   compute_field_inversion_height,          &
    &                                   compute_field_dbz850,                    &
    &                                   compute_field_dbzlmx,                    &
    &                                   compute_field_dbzcmax,                   &
    &                                   compute_field_smi,                       &
    &                                   compute_field_lapserate,                 &
    &                                   compute_field_srh,                       &
    &                                   compute_field_wshear
  USE mo_io_config,               ONLY: itype_pres_msl, itype_rh,                &
    &                                   n_wshear, wshear_uv_heights, n_srh, srh_heights
  USE mo_grid_config,             ONLY: l_limited_area, n_dom_start
  USE mo_interpol_config,         ONLY: support_baryctr_intp
  USE mo_advection_config,        ONLY: advection_config
  USE mo_fortran_tools,           ONLY: init, copy, assert_acc_device_only, assert_acc_host_only, &
   &                                    assert_lacc_equals_i_am_accel_node
#ifdef _OPENACC
  USE mo_mpi,                     ONLY: i_am_accel_node
#endif

  ! Workaround for SMI computation. Not nice, however by making 
  ! direct use of the states below, we avoid enhancing the type t_data_input.
  USE mo_nwp_lnd_state,           ONLY: p_lnd_state
  USE mo_ext_data_state,          ONLY: ext_data

  IMPLICIT NONE

  ! interface definition
  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_pp_tasks'

  ! max. name string length
  INTEGER, PARAMETER, PUBLIC :: MAX_NAME_LENGTH   =   256

  ! priority levels for tasks (smaller is earlier):
  INTEGER, PARAMETER, PUBLIC  :: HIGH_PRIORITY     =    0  
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY0 =    8  
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY1 =    9   
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY2 =   10  
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY3 =   11  
  INTEGER, PARAMETER, PUBLIC  :: DEFAULT_PRIORITY4 =   12  
  INTEGER, PARAMETER, PUBLIC  :: LOW_PRIORITY      =  100  

  ! level of output verbosity
  INTEGER, PUBLIC :: dbg_level = 0

  ! functions and subroutines
  PUBLIC :: pp_task_lonlat
  PUBLIC :: pp_task_sync
  PUBLIC :: pp_task_ipzlev_setup
  PUBLIC :: pp_task_ipzlev
  PUBLIC :: pp_task_intp_msl
  PUBLIC :: pp_task_compute_field
  PUBLIC :: pp_task_edge2cell
  ! variables
  PUBLIC :: job_queue
  ! data types
  PUBLIC :: t_data_input
  PUBLIC :: t_data_output
  PUBLIC :: t_simulation_status
  PUBLIC :: t_activity_status
  PUBLIC :: t_job_queue

  !--- JOB QUEUE DEFINITION ----------------------------------------------------------

  !> data necessary for job input.
  !
  !  This type is likely to contain more data than really needed for
  !  your specific post-processing job (we are lacking polymorphic
  !  data structures).
  !
  !  @note Please avoid using any non-static data from other places
  !  for your post-processing tasks!  The only exceptions are the use of
  !  the global, volatile "nnow" value and "nproma"!
  !
  !  @note Elements of this type are COPIED. Therefore avoid large
  !  data structures or use POINTERs.
  TYPE t_data_input
    ! pointer for model variable (array)
    TYPE (t_var),            POINTER :: var
    INTEGER                          :: jg ! domain ID
    TYPE(t_patch),           POINTER :: p_patch         
    TYPE(t_int_state),       POINTER :: p_int_state     
    TYPE(t_nh_state),        POINTER :: p_nh_state      
    TYPE(t_nwp_phy_diag),    POINTER :: prm_diag        
    TYPE(t_nh_opt_diag),     POINTER :: p_nh_opt_diag   
    TYPE(t_nh_pzlev_config), POINTER :: nh_pzlev_config 
  END TYPE t_data_input


  !> data necessary for job output.
  !  See also @p t_data_input.
  TYPE t_data_output
    ! pointer for model variable (array)
    TYPE (t_var), POINTER :: var    => NULL()
    ! (optional) pointer for second component of model variable.
    ! necessary for lon-lat interpolation of edge-based fields.
    TYPE (t_var), POINTER :: var_2  => NULL()
  END TYPE t_data_output


  !> Definition of simulation status, a list of LOGICAL flags like
  !  "first_step", "last_step", "output_time"
  !
  !  Based on these values, we determine if a post-processing
  !  task is "active" and will be processed.
  !
  !  Flags are stored in a contiguous LOGICAL array to make
  !  Boolean comparisons more convenient.
  !
  !  @note There might be better places in the code for such a
  !  variable!
  TYPE t_simulation_status
    LOGICAL :: status_flags(4)           !< l_output_step, l_first_step, l_last_step, l_accumulation_step
    LOGICAL :: ldom_active(max_dom)      !< active domains
    INTEGER :: i_timelevel_dyn(max_dom)  !< active time level (for dynamics output variables related to nnow)
    INTEGER :: i_timelevel_phy(max_dom)  !< active time level (for physics output variables related to nnow_rcf)
  END TYPE t_simulation_status


  !> Definition of task activity, i.e. settings when a task should be
  !> triggered.
  TYPE t_activity_status
    LOGICAL :: status_flags(4)           !< l_output_step, l_first_step, l_last_step
    LOGICAL :: check_dom_active          !< check if this task's domain is active
    INTEGER :: i_timelevel               !< time level for this task
  END TYPE t_activity_status


  !> A variable of type @p t_job_queue defines a single post-processing task.
  !
  !  Jobs with smaller priority values are processed first.
  TYPE t_job_queue
    
    INTEGER                         :: job_priority   !< Task priority.
    CHARACTER(len=MAX_NAME_LENGTH)  :: job_name       !< job name string (for status output)
    INTEGER                         :: job_type       !< task type (quasi function pointer)
    TYPE(t_activity_status)         :: activity       !< "under which conditions does this task run?"

    TYPE(t_data_input)              :: data_input     !< input of post-processing task
    TYPE(t_data_output)             :: data_output    !< result of post-processing task

    TYPE(t_job_queue), POINTER      :: next => NULL() !< pointer to next element in list

  END TYPE t_job_queue


  !--- MODULE DATA -------------------------------------------------------------------
  TYPE(t_job_queue), POINTER   :: job_queue  =>  NULL() !< head of (ordered) job queue

#ifndef _OPENACC
  LOGICAL, PARAMETER :: i_am_accel_node=.FALSE.
#endif

CONTAINS

  !--- POST-PROCESSING TASKS ---------------------------------------------------------

  !---------------------------------------------------------------
  !> Performs interpolation of a variable onto a regular grid.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_lonlat(ptr_task)
    TYPE(t_job_queue), TARGET :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_lonlat"
    INTEGER                            ::        &
      &  lonlat_id, jg,                          &
      &  in_var_idx, out_var_idx, out_var_idx_2, &
      &  ierrstat, dim1, dim2, dim3, hintp_type
    TYPE (t_var), POINTER :: in_var, out_var, out_var_2
    TYPE (t_var_metadata),     POINTER :: p_info
    TYPE (t_lon_lat_intp),     POINTER :: ptr_int_lonlat
    REAL(wp), ALLOCATABLE, TARGET      :: tmp_var(:,:,:)
    INTEGER,  ALLOCATABLE, TARGET      :: tmp_int_var(:,:,:)
    REAL(wp), POINTER                  :: tmp_ptr(:,:,:)
    INTEGER,  POINTER                  :: tmp_int_ptr(:,:,:)
    TYPE(t_patch),             POINTER :: p_patch
    INTEGER                            :: var_ref_pos

    p_patch        => ptr_task%data_input%p_patch      ! patch
    p_info         => ptr_task%data_input%var%info
    in_var         => ptr_task%data_input%var
    out_var        => ptr_task%data_output%var
    lonlat_id      =  ptr_task%data_output%var%info%hor_interp%lonlat_id
    jg             =  ptr_task%data_input%jg
    ptr_int_lonlat => lonlat_grids%list(lonlat_id)%intp(jg)
    hintp_type     = p_info%hor_interp%hor_intp_type

#ifdef _OPENACC
      IF(i_am_accel_node .AND. .NOT. p_info%lopenacc) THEN
        CALL message(routine, "WARNING: " // TRIM(p_info%name) // " is temporarily copied onto accelerator. If " // &
          "this warning appears at every output step, consider making the variable present permanently. However " // &
          "as this variable is not present on the device it is unlikely that it changes over time as it is not " // &
          "part of any physics. Therefore one might want to limit the output frequency of " // TRIM(p_info%name))
        ! MJ: Some fields that are constant over time are not necessary to be available on the devices. Their memory
        ! can be saved. In principle it would be enough to output them only at the first time step. Thus this warning
        ! should only appear once.
        IF(ASSOCIATED(in_var%r_ptr)) THEN
          !$ACC ENTER DATA COPYIN(in_var%r_ptr)
        ELSEIF(ASSOCIATED(in_var%s_ptr)) THEN
          !$ACC ENTER DATA COPYIN(in_var%s_ptr)
        ELSEIF(ASSOCIATED(in_var%i_ptr)) THEN
          !$ACC ENTER DATA COPYIN(in_var%i_ptr)
        ENDIF
      ENDIF
#endif

    ! --------------------------------------------------------------------------
    !
    ! IMPORTANT: Currently, barycentric interpolation supported only
    !            - if the namelist parameter "interpol_nml/support_baryctr_intp"
    !              has been set to .TRUE.
    !
    ! If these two prerequisites are not met, then we choose a different
    ! interpolation algorithm as a fallback option. This algorithm is specified
    ! by "hor_interp%fallback_type".
    ! --------------------------------------------------------------------------
    IF (hintp_type == HINTP_TYPE_LONLAT_BCTR) THEN
      IF (.NOT. support_baryctr_intp) THEN
        hintp_type = p_info%hor_interp%fallback_type
      END IF
    END IF
    in_var_idx        = 1
    IF (in_var%info%lcontained)  in_var_idx  = in_var%info%ncontained
    out_var_idx       = 1
    IF (out_var%info%lcontained) out_var_idx = out_var%info%ncontained
    ! For edge-based interpolation: retrieve data on Y-component:
    IF (ASSOCIATED(ptr_task%data_output%var_2)) THEN
      out_var_2      => ptr_task%data_output%var_2
      out_var_idx_2  =  1
      IF (out_var_2%info%lcontained) out_var_idx_2 = out_var_2%info%ncontained
    END IF
    IF (zaxisTypeList%is_2d(p_info%vgrid) .AND. (p_info%ndims /= 2)) THEN
      CALL finish(routine, "Inconsistent dimension info!")
    END IF

    SELECT CASE (p_info%hgrid)
    CASE (GRID_UNSTRUCTURED_CELL)
      IF (ASSOCIATED(in_var%r_ptr) .OR. ASSOCIATED(in_var%s_ptr)) THEN

        ! --------------------------------------
        ! REAL and SINGLE PRECISION FLOAT fields
        ! --------------------------------------

        IF (zaxisTypeList%is_2d(p_info%vgrid)) THEN
          ! A 2D variable (nproma, nblks) is copied a to 1-level 3D variable 
          ! (nproma, nlevs=1, nblks). This requires a temporary variable:

          var_ref_pos = 3
          IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos

          IF (var_ref_pos /= 2 .OR. ASSOCIATED(in_var%s_ptr)) THEN
            dim1 = p_info%used_dimensions(1)
            dim2 = p_info%used_dimensions(2)
            ALLOCATE(tmp_var(dim1, 1, dim2), STAT=ierrstat)
            IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
            !$ACC ENTER DATA CREATE(tmp_var) IF(i_am_accel_node)
          ENDIF

          IF (ASSOCIATED(in_var%r_ptr)) THEN

            SELECT CASE(var_ref_pos)
            CASE (1)
              !$OMP PARALLEL
              CALL copy(in_var%r_ptr(in_var_idx,:,:,1,1), tmp_var(:,1,:))
              !$OMP END PARALLEL
              tmp_ptr => tmp_var(:,:,:)
            CASE (2)
              ! no need to copy in this particular case (the second dim has already length 1)
              tmp_ptr => in_var%r_ptr(:,in_var_idx:in_var_idx,:,1,1)
            CASE (3)
              !$OMP PARALLEL
              CALL copy(in_var%r_ptr(:,:,in_var_idx,1,1), tmp_var(:,1,:))
              !$OMP END PARALLEL
              tmp_ptr => tmp_var(:,:,:)
            CASE default
              CALL finish(routine, "internal error!")
            END SELECT

          ELSE IF (ASSOCIATED(in_var%s_ptr)) THEN

            ! A SP variable has to be copied to a temporary DP array.

            SELECT CASE(var_ref_pos)
            CASE (1)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(in_var_idx,:,:,1,1), tmp_var(:,1,:))
              !$OMP END PARALLEL
            CASE (2)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(:,in_var_idx,:,1,1), tmp_var(:,1,:))
              !$OMP END PARALLEL
            CASE (3)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(:,:,in_var_idx,1,1), tmp_var(:,1,:))
              !$OMP END PARALLEL
            CASE default
              CALL finish(routine, "internal error!")
            END SELECT
            tmp_ptr => tmp_var(:,:,:)
          ELSE
            CALL finish(routine, "internal error!")
          ENDIF

        ELSE

          var_ref_pos = 4
          IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos

          IF (ASSOCIATED(in_var%r_ptr)) THEN
            SELECT CASE(var_ref_pos)
            CASE (1)
              tmp_ptr => in_var%r_ptr(in_var_idx,:,:,:,1)
            CASE (2)
              tmp_ptr => in_var%r_ptr(:,in_var_idx,:,:,1)
            CASE (3)
              tmp_ptr => in_var%r_ptr(:,:,in_var_idx,:,1)
            CASE (4)
              tmp_ptr => in_var%r_ptr(:,:,:,in_var_idx,1)
            CASE default
              CALL finish(routine, "internal error!")
            END SELECT
          ELSE  IF (ASSOCIATED(in_var%s_ptr)) THEN
            ! A SP variable has to be copied to a temporary DP array.
            SELECT CASE(var_ref_pos)
            CASE (1)
              dim1 = SIZE(in_var%s_ptr,2)
              dim2 = SIZE(in_var%s_ptr,3)
              dim3 = SIZE(in_var%s_ptr,4)
              ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
              IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
              !$ACC ENTER DATA CREATE(tmp_var) IF(i_am_accel_node)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(in_var_idx,:,:,:,1), tmp_var)
              !$OMP END PARALLEL
            CASE (2)
              dim1 = SIZE(in_var%s_ptr,1)
              dim2 = SIZE(in_var%s_ptr,3)
              dim3 = SIZE(in_var%s_ptr,4)
              ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
              IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
              !$ACC ENTER DATA CREATE(tmp_var) IF(i_am_accel_node)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(:,in_var_idx,:,:,1), tmp_var)
              !$OMP END PARALLEL
            CASE (3)
              dim1 = SIZE(in_var%s_ptr,1)
              dim2 = SIZE(in_var%s_ptr,2)
              dim3 = SIZE(in_var%s_ptr,4)
              ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
              IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
              !$ACC ENTER DATA CREATE(tmp_var) IF(i_am_accel_node)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(:,:,in_var_idx,:,1), tmp_var)
              !$OMP END PARALLEL
            CASE (4)
              dim1 = SIZE(in_var%s_ptr,1)
              dim2 = SIZE(in_var%s_ptr,2)
              dim3 = SIZE(in_var%s_ptr,3)
              ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
              IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
              !$ACC ENTER DATA CREATE(tmp_var) IF(i_am_accel_node)
              !$OMP PARALLEL
              CALL copy(in_var%s_ptr(:,:,:,in_var_idx,1), tmp_var)
              !$OMP END PARALLEL
            CASE default
              CALL finish(routine, "internal error!")
            END SELECT
            tmp_ptr => tmp_var(:,:,:)
          ELSE
            CALL finish(routine, "internal error!")
          ENDIF
        END IF ! 2D

        ! for cell-based variables: interpolate gradients (finite
        ! differences) and reconstruct
        CALL ptr_int_lonlat%interpolate(          &
          &   TRIM(p_info%name), tmp_ptr, nproma, &
          &   out_var%r_ptr(:,:,:,out_var_idx,1), &
          &   hintp_type)

      ELSE IF (ASSOCIATED(in_var%i_ptr)) THEN

        ! --------------
        ! INTEGER fields
        ! --------------

        IF (zaxisTypeList%is_2d(p_info%vgrid)) THEN
          ! A 2D variable (nproma, nblks) is copied a to 1-level 3D variable 
          ! (nproma, nlevs=1, nblks). This requires a temporary variable:
          var_ref_pos = 3
          IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos

          IF (var_ref_pos /= 2) THEN
            dim1 = p_info%used_dimensions(1)
            dim2 = p_info%used_dimensions(2)
            ALLOCATE(tmp_int_var(dim1, 1, dim2), STAT=ierrstat)
            IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_int_var failed')
            !$ACC ENTER DATA CREATE(tmp_int_var) IF(i_am_accel_node)
          ENDIF

          SELECT CASE(var_ref_pos)
          CASE (1)
            CALL copy(in_var%i_ptr(in_var_idx,:,:,1,1), tmp_int_var(:,1,:))
            tmp_int_ptr => tmp_int_var
          CASE (2)
            ! no need to copy in this particular case (the second dim has already length 1)
            tmp_int_ptr => in_var%i_ptr(:,in_var_idx:in_var_idx,:,1,1)
          CASE (3)
            CALL copy(in_var%i_ptr(:,:,in_var_idx,1,1), tmp_int_var(:,1,:))
            tmp_int_ptr => tmp_int_var
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT

        ELSE

          var_ref_pos = 4
          IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos
          SELECT CASE(var_ref_pos)
          CASE (1)
            tmp_int_ptr => in_var%i_ptr(in_var_idx,:,:,:,1)
          CASE (2)
            tmp_int_ptr => in_var%i_ptr(:,in_var_idx,:,:,1)
          CASE (3)
            tmp_int_ptr => in_var%i_ptr(:,:,in_var_idx,:,1)
          CASE (4)
            tmp_int_ptr => in_var%i_ptr(:,:,:,in_var_idx,1)
          CASE default
            CALL finish(routine, "internal error!")
          END SELECT

        END IF ! 2D

        ! for cell-based variables: interpolate gradients (finite
        ! differences) and reconstruct
        CALL ptr_int_lonlat%interpolate(               &
          &   TRIM(p_info%name), tmp_int_ptr, nproma,  &
          &   out_var%i_ptr(:,:,:,out_var_idx,1),      &
          &   hintp_type)

        IF (ALLOCATED(tmp_int_var)) THEN
          ! clean up:
          !$ACC WAIT IF(i_am_accel_node)
          !$ACC EXIT DATA DELETE(tmp_int_var) IF(i_am_accel_node)
          DEALLOCATE(tmp_int_var, STAT=ierrstat)
          IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation failed')
        END IF

      END IF

      ! --------------------------------------------------------------
      !
    CASE (GRID_UNSTRUCTURED_EDGE)
      ! throw error message, if this variable is not a REAL field:
      IF (.NOT. ASSOCIATED(in_var%r_ptr)) THEN
        CALL finish(routine, TRIM(p_info%name)//": Interpolation not implemented.")
      END IF

      IF (zaxisTypeList%is_2d(p_info%vgrid)) THEN
        ! For 2D variables (nproma, nblks) we use a 1-level 3D pointer (nproma, nlevs=1, nblks).

        var_ref_pos = 3
        IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos

        IF (var_ref_pos /= 2) THEN
          dim1 = p_info%used_dimensions(1)
          dim2 = p_info%used_dimensions(2)
          ALLOCATE(tmp_var(dim1, 1, dim2), STAT=ierrstat)
          IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
          !$ACC ENTER DATA CREATE(tmp_var) IF(i_am_accel_node)
        ENDIF

        SELECT CASE(var_ref_pos)
        CASE (1)
          !$OMP PARALLEL
          CALL copy(in_var%r_ptr(in_var_idx,:,:,1,1), tmp_var(:,1,:))
          !$OMP END PARALLEL
          tmp_ptr => tmp_var(:,:,:)
        CASE (2)
          ! no need to copy in this particular case (the second dim has already length 1)
          tmp_ptr => in_var%r_ptr(:,in_var_idx:in_var_idx,:,1,1)
        CASE (3)
          !$OMP PARALLEL
          CALL copy(in_var%r_ptr(:,:,in_var_idx,1,1), tmp_var(:,1,:))
          !$OMP END PARALLEL
          tmp_ptr => tmp_var(:,:,:)
        CASE default
          CALL finish(routine, "internal error!")
        END SELECT

      ELSE

        var_ref_pos = 4
        IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos
        SELECT CASE(var_ref_pos)
        CASE (1)
          tmp_ptr => in_var%r_ptr(in_var_idx,:,:,:,1)
        CASE (2)
          tmp_ptr => in_var%r_ptr(:,in_var_idx,:,:,1)
        CASE (3)
          tmp_ptr => in_var%r_ptr(:,:,in_var_idx,:,1)
        CASE (4)
          tmp_ptr => in_var%r_ptr(:,:,:,in_var_idx,1)
        CASE default
          CALL finish(routine, "internal error!")
        END SELECT
      END IF ! 2D

      ! for edge-based variables: simple interpolation
      CALL ptr_int_lonlat%interpolate( tmp_ptr, nproma,                            &
        &                              out_var%r_ptr(:,:,:,out_var_idx,1),         &
        &                              out_var_2%r_ptr(:,:,:,out_var_idx_2,1),     &
        &                              HINTP_TYPE_LONLAT_RBF )

    CASE DEFAULT
      CALL finish(routine, 'Unknown grid type.')
    END SELECT

    !$ACC WAIT IF(i_am_accel_node)

    ! clean up
    IF (ALLOCATED(tmp_var)) THEN
      !$ACC EXIT DATA DELETE(tmp_var) IF(i_am_accel_node)
      DEALLOCATE(tmp_var, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation of tmp_var failed')
    ENDIF
    IF(i_am_accel_node .AND. .NOT. p_info%lopenacc) THEN
      IF(ASSOCIATED(in_var%r_ptr)) THEN
        !$ACC EXIT DATA DELETE(in_var%r_ptr)
      ELSEIF(ASSOCIATED(in_var%s_ptr)) THEN
        !$ACC EXIT DATA DELETE(in_var%s_ptr)
      ELSEIF(ASSOCIATED(in_var%i_ptr)) THEN
        !$ACC EXIT DATA DELETE(in_var%i_ptr)
      ENDIF
    ENDIF
  END SUBROUTINE pp_task_lonlat


  !---------------------------------------------------------------
  !> Performs synchronization of halo regions.
  !  This is necessary before starting the lon-lat interpolation.
  !  All variables that are part of an interpolation task are
  !  synchronized.
  !
  ! To avoid unnecessary overhead in the synchronization,
  ! several improvements are possible:
  !
  !   - Copy 2D fields into a 3D field which is
  !     synchronized. Afterwards, 2D fields are extracted again from
  !     the temporary 3D field.
  !
  !   - Introduce some kind of meta information of each variable list.
  !     For example, prognostic fields must not be synchronized,
  !     therefore the corresponding variable lists can be marked as
  !     "skip_sync".
  !
  SUBROUTINE pp_task_sync(sim_status)
    TYPE(t_simulation_status),  INTENT(IN) :: sim_status
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_sync"
    TYPE(t_job_queue),         POINTER :: ptr_task
    INTEGER :: in_var_idx, jg, sync_mode, var_ref_pos, timelevel
    TYPE (t_var), POINTER :: in_var
    TYPE (t_var_metadata), POINTER :: p_info
    TYPE(t_patch), POINTER :: p_patch

    ptr_task => job_queue
    ! loop over job queue
    LOOP_JOB : DO
      IF (.NOT. ASSOCIATED(ptr_task)) EXIT
      IF (ptr_task%job_type == TASK_INTP_HOR_LONLAT) THEN
        p_patch     => ptr_task%data_input%p_patch
        p_info      => ptr_task%data_input%var%info
        jg          =  p_patch%id
        SELECT CASE (p_info%tlev_source)
        CASE(TLEV_NNOW);     timelevel = sim_status%i_timelevel_dyn(jg)
        CASE(TLEV_NNOW_RCF); timelevel = sim_status%i_timelevel_phy(jg)
        CASE DEFAULT
          CALL finish(routine, 'Unsupported tlev_source')
        END SELECT
        IF ((ptr_task%activity%i_timelevel == timelevel) .OR.  &
          & (ptr_task%activity%i_timelevel == ALL_TIMELEVELS))  THEN
          in_var      => ptr_task%data_input%var
          in_var_idx  =  1
          IF (in_var%info%lcontained) in_var_idx = in_var%info%ncontained
          IF (zaxisTypeList%is_2d(p_info%vgrid) .AND. (p_info%ndims /= 2)) &
            &  CALL finish(routine, "Inconsistent dimension info!")
          IF (dbg_level >= 10) & 
               CALL message(routine, "synchronize variable "//TRIM(p_info%name))
          SELECT CASE (p_info%hgrid)
          CASE (GRID_UNSTRUCTURED_CELL)
            sync_mode = SYNC_C
          CASE (GRID_UNSTRUCTURED_EDGE)
            sync_mode = SYNC_E
          CASE DEFAULT
            CALL finish(routine, 'Unknown grid type.')
          END SELECT
          IF (zaxisTypeList%is_2d(p_info%vgrid)) THEN
            var_ref_pos = 3
            IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos
            IF (ASSOCIATED(in_var%r_ptr)) THEN
              SELECT CASE(var_ref_pos)
              CASE (1)
                CALL sync_patch_array(sync_mode, p_patch, in_var%r_ptr(in_var_idx,:,:,1,1) )
              CASE (2)
                CALL sync_patch_array(sync_mode, p_patch, in_var%r_ptr(:,in_var_idx,:,1,1) )
              CASE (3)
                CALL sync_patch_array(sync_mode, p_patch, in_var%r_ptr(:,:,in_var_idx,1,1) )
              CASE default
                CALL finish(routine, "internal error!")
              END SELECT
            END IF
            IF (ASSOCIATED(in_var%i_ptr)) THEN
              SELECT CASE(var_ref_pos)
              CASE (1)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(in_var_idx,:,:,1,1) )
              CASE (2)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(:,in_var_idx,:,1,1) )
              CASE (3)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(:,:,in_var_idx,1,1) )
              CASE default
                CALL finish(routine, "internal error!")
              END SELECT
            END IF
          ELSE
            var_ref_pos = 4
            IF (in_var%info%lcontained)  var_ref_pos = in_var%info%var_ref_pos
            IF (ASSOCIATED(in_var%r_ptr)) THEN
              SELECT CASE(var_ref_pos)
              CASE (1)
                CALL cumulative_sync_patch_array(sync_mode, p_patch, in_var%r_ptr(in_var_idx,:,:,:,1))
              CASE (2)
                CALL cumulative_sync_patch_array(sync_mode, p_patch, in_var%r_ptr(:,in_var_idx,:,:,1))
              CASE (3)
                CALL cumulative_sync_patch_array(sync_mode, p_patch, in_var%r_ptr(:,:,in_var_idx,:,1))
              CASE (4)
                CALL cumulative_sync_patch_array(sync_mode, p_patch, in_var%r_ptr(:,:,:,in_var_idx,1))
              CASE default
                CALL finish(routine, "internal error!")
              END SELECT
            END IF
            IF (ASSOCIATED(in_var%i_ptr)) THEN
              SELECT CASE(var_ref_pos)
              CASE (1)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(in_var_idx,:,:,:,1))
              CASE (2)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(:,in_var_idx,:,:,1))
              CASE (3)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(:,:,in_var_idx,:,1))
              CASE (4)
                CALL sync_patch_array(sync_mode, p_patch, in_var%i_ptr(:,:,:,in_var_idx,1))
              CASE default
                CALL finish(routine, "internal error!")
              END SELECT
            END IF
          END IF
        END IF
      END IF
      ptr_task => ptr_task%next
    END DO LOOP_JOB
    ! complete pending syncs:
    CALL complete_cumulative_sync()
  END SUBROUTINE pp_task_sync

  !---------------------------------------------------------------
  !> Performs setup of vertical interpolation.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  SUBROUTINE pp_task_ipzlev_setup(ptr_task, lacc)
    TYPE(t_job_queue), POINTER :: ptr_task
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_ipzlev_setup"
    INTEGER                            :: jg, nzlev, nplev, nilev
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    

    ! prognostic state: note that we only use p_prog(nnow(jg))
    TYPE(t_nh_prog),           POINTER :: p_prog
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config
    TYPE(t_int_state),         POINTER :: intp_hrz

    CALL assert_acc_device_only(routine, lacc)

    ! patch, state, and metrics
    jg             =  ptr_task%data_input%jg
    p_patch        => ptr_task%data_input%p_patch
    p_metrics      => ptr_task%data_input%p_nh_state%metrics
    p_prog         => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_diag         => ptr_task%data_input%p_nh_state%diag
    p_diag_pz      => ptr_task%data_input%p_nh_opt_diag%diag_pz
    prm_diag       => ptr_task%data_input%prm_diag
    intp_hrz       => ptr_task%data_input%p_int_state

    ! ipz-level interpolation data
    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config

    nzlev          =  nh_pzlev_config%zlevels%nvalues
    nplev          =  nh_pzlev_config%plevels%nvalues
    nilev          =  nh_pzlev_config%ilevels%nvalues

    ! build data structure "vcoeff" containing coefficient tables                      
    SELECT CASE ( ptr_task%job_type )
    CASE ( TASK_INIT_VER_Z )
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_Z")
      CALL prepare_vert_interp_z(p_patch, p_diag, p_metrics, intp_hrz, nzlev,          &
        &                        p_diag_pz%z_temp, p_diag_pz%z_pres,                   &
        &                        nh_pzlev_config%z3d, p_diag_pz%vcoeff_z, lacc=.TRUE.)
      !
    CASE ( TASK_INIT_VER_P )
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_P")
      CALL prepare_vert_interp_p(p_patch, p_diag, p_metrics, intp_hrz, nplev,          &
        &                        p_diag_pz%p_gh, p_diag_pz%p_temp,                     &
        &                        nh_pzlev_config%p3d, p_diag_pz%vcoeff_p, lacc=.TRUE.)
      !
    CASE ( TASK_INIT_VER_I )
      IF (dbg_level >= 10)  CALL message(routine, "TASK_INIT_VER_I")
      CALL prepare_vert_interp_i(p_patch, p_prog, p_diag, p_metrics, intp_hrz, nilev,  &
        &                        p_diag_pz%i_gh, p_diag_pz%i_temp,                     &
        &                        nh_pzlev_config%i3d, p_diag_pz%vcoeff_i, lacc=.TRUE.)
      !
    CASE ( TASK_FINALIZE_IPZ )
      ! deallocate coefficient tables:
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_z, lacc=.TRUE.)
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_p, lacc=.TRUE.)
      CALL vcoeff_deallocate(p_diag_pz%vcoeff_i, lacc=.TRUE.)
      !
    CASE DEFAULT
      CALL finish(routine, "Internal error.")
    END SELECT ! vert_intp_method

  END SUBROUTINE pp_task_ipzlev_setup


  !---------------------------------------------------------------
  !> Performs vertical interpolation of a variable onto i/p/z-levels.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  !
  SUBROUTINE pp_task_ipzlev(ptr_task, lacc)
    TYPE(t_job_queue), POINTER, INTENT(INOUT) :: ptr_task
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_ipzlev"
    INTEGER                            :: &
      &  vert_intp_method,                        &
      &  in_var_idx, out_var_idx, nlev, nlevp1,   &
      &  n_ipzlev, npromz, nblks, ierrstat,       &
      &  in_var_ref_pos, out_var_ref_pos,         &
      &  dim1, dim2, dim3
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nh_diag_pz),        POINTER :: p_diag_pz
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag
    TYPE(t_vert_interp_meta),  POINTER :: pzlev_flags

    TYPE (t_var), POINTER :: in_var, out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    TYPE(t_vcoeff),            POINTER :: vcoeff
    TYPE(t_nh_pzlev_config),   POINTER :: nh_pzlev_config
    REAL(wp),                  POINTER :: p_z3d(:,:,:), p_pres(:,:,:), p_temp(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET      :: z_me(:,:,:), p_z3d_edge(:,:,:)
    REAL(wp), ALLOCATABLE, TARGET      :: tmp_var(:,:,:)
    REAL(wp),                  POINTER :: in_z3d(:,:,:), in_z_mc(:,:,:)
    TYPE(t_int_state),         POINTER :: intp_hrz

    LOGICAL                            :: &
      &  l_hires_intp, l_restore_fricred, l_loglin, &
      &  l_extrapol, l_satlimit, l_restore_pbldev,  &
      &  l_pd_limit
    REAL(wp)                           :: lower_limit
    TYPE (t_vcoeff_lin), POINTER       :: vcoeff_lin, vcoeff_lin_nlevp1
    TYPE (t_vcoeff_cub), POINTER       :: vcoeff_cub

    REAL(wp), POINTER :: in_ptr(:,:,:), out_ptr(:,:,:)

    CALL assert_acc_device_only(routine, lacc)
    CALL assert_lacc_equals_i_am_accel_node(routine, lacc)

    ! input/output field for this task
    p_info            => ptr_task%data_input%var%info
    in_var            => ptr_task%data_input%var
    out_var           => ptr_task%data_output%var

    in_var_idx        = 1
    in_var_ref_pos    = 4
    IF (ptr_task%data_input%var%info%lcontained) THEN
      in_var_idx      = ptr_task%data_input%var%info%ncontained
      in_var_ref_pos  = ptr_task%data_input%var%info%var_ref_pos
    END IF
    out_var_idx       = 1
    out_var_ref_pos   = 4
    IF (ptr_task%data_output%var%info%lcontained) THEN
      out_var_idx     = ptr_task%data_output%var%info%ncontained
      out_var_ref_pos = ptr_task%data_output%var%info%var_ref_pos
    END IF

    !--- load some items from input/output data structures
    vert_intp_method = p_info%vert_interp%vert_intp_method

    ! patch, state, and metrics
    p_patch           => ptr_task%data_input%p_patch
    p_metrics         => ptr_task%data_input%p_nh_state%metrics
    p_diag            => ptr_task%data_input%p_nh_state%diag
    p_diag_pz         => ptr_task%data_input%p_nh_opt_diag%diag_pz
    prm_diag          => ptr_task%data_input%prm_diag
    intp_hrz          => ptr_task%data_input%p_int_state

    nh_pzlev_config   => ptr_task%data_input%nh_pzlev_config
    nlev              = p_patch%nlev
    nlevp1            = p_patch%nlevp1

    ! pz-level interpolation data
    SELECT CASE ( ptr_task%job_type )
    CASE ( TASK_INTP_VER_ZLEV )
      ! vertical levels for z-level interpolation
      n_ipzlev  =  nh_pzlev_config%zlevels%nvalues
      vcoeff  =>  p_diag_pz%vcoeff_z
      p_z3d   =>  nh_pzlev_config%z3d
      p_pres  =>  p_diag_pz%z_pres
      p_temp  =>  p_diag_pz%z_temp
    CASE ( TASK_INTP_VER_PLEV )
      ! vertical levels for p-level interpolation
      n_ipzlev  =  nh_pzlev_config%plevels%nvalues
      vcoeff  =>  p_diag_pz%vcoeff_p
      p_z3d   =>  p_diag_pz%p_gh
      p_pres  =>  nh_pzlev_config%p3d
      p_temp  =>  p_diag_pz%p_temp
    CASE ( TASK_INTP_VER_ILEV )
      ! vertical levels for isentropic-level interpolation
      n_ipzlev  =   nh_pzlev_config%ilevels%nvalues
      vcoeff  =>  p_diag_pz%vcoeff_i
      p_z3d   =>  p_diag_pz%i_gh
      p_pres  =>  nh_pzlev_config%p3d ! ** this still needs to be fixed! we either need i_pres here
      p_temp  =>  p_diag_pz%i_temp    !    or have to turn off the saturation adjustment for theta levels
                                      !    or have to use log-linear interpolation in this case **
      IF (vert_intp_method == VINTP_METHOD_QV) & ! Let's stop with an error message for the time being
        CALL finish(routine, "QV interpolation to isentropic levels not available.")
    CASE DEFAULT
      CALL finish(routine, "Unknown post-processing job.")
    END SELECT
                     
    ! interpolation flags + parameters
    pzlev_flags => in_var%info%vert_interp
    l_hires_intp      = pzlev_flags%l_hires_intp      
    l_restore_fricred = pzlev_flags%l_restore_fricred 
    l_loglin          = pzlev_flags%l_loglin          
    l_extrapol        = pzlev_flags%l_extrapol        
    l_satlimit        = pzlev_flags%l_satlimit        
    l_restore_pbldev  = pzlev_flags%l_restore_pbldev  
    l_pd_limit        = pzlev_flags%l_pd_limit
    lower_limit       = pzlev_flags%lower_limit       

    !-- perform some consistency checks
    IF (p_info%ndims /= 3) &
      & CALL finish(routine, "Wrong number of variables dimensions!")
    IF (.NOT. vcoeff%l_initialized) &
      CALL finish(routine, "Interpolation coefficients not yet initialized!")

    SELECT CASE ( p_info%hgrid )
    CASE (GRID_UNSTRUCTURED_CELL) 
      nblks  = p_patch%nblks_c
      npromz = p_patch%npromz_c
      ! 
      vcoeff_lin        => vcoeff%lin_cell
      vcoeff_lin_nlevp1 => vcoeff%lin_cell_nlevp1
      vcoeff_cub        => vcoeff%cub_cell
      in_z3d            => p_z3d
      in_z_mc           => p_metrics%z_mc

    CASE (GRID_UNSTRUCTURED_EDGE) 
      nblks  = p_patch%nblks_e
      npromz = p_patch%npromz_e
      !
      vcoeff_lin        => vcoeff%lin_edge
      vcoeff_lin_nlevp1 => NULL()
      vcoeff_cub        => vcoeff%cub_edge

      ! Compute geometric height at edge points (temporary variable)
      ALLOCATE(p_z3d_edge(nproma,n_ipzlev,nblks), z_me(nproma,p_patch%nlev,nblks), STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'ALLOCATE failed')
      !$ACC ENTER DATA CREATE(p_z3d_edge, z_me)

      CALL cells2edges_scalar(p_metrics%z_mc, p_patch, intp_hrz%c_lin_e,    &
        &                     z_me, opt_fill_latbc=.TRUE., lacc=.TRUE.)
      CALL cells2edges_scalar(p_z3d, p_patch, intp_hrz%c_lin_e, p_z3d_edge, &
        &                     opt_fill_latbc=.TRUE., lacc=.TRUE.)

      in_z3d            => p_z3d_edge
      in_z_mc           => z_me
    END SELECT

    IF (ASSOCIATED(in_var%r_ptr)) THEN
      SELECT CASE(in_var_ref_pos)
      CASE (1)
        in_ptr => in_var%r_ptr(in_var_idx,:,:,:,1)
      CASE (2)
        in_ptr => in_var%r_ptr(:,in_var_idx,:,:,1)
      CASE (3)
        in_ptr => in_var%r_ptr(:,:,in_var_idx,:,1)
      CASE (4)
        in_ptr => in_var%r_ptr(:,:,:,in_var_idx,1)
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
    ELSE IF (ASSOCIATED(in_var%s_ptr)) THEN
      SELECT CASE(in_var_ref_pos)
      CASE (1)
        dim1 = SIZE(in_var%s_ptr,2)
        dim2 = SIZE(in_var%s_ptr,3)
        dim3 = SIZE(in_var%s_ptr,4)
        ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
        !$ACC ENTER DATA CREATE(tmp_var)
!$OMP PARALLEL
        CALL copy(in_var%s_ptr(in_var_idx,:,:,:,1), tmp_var(:,:,:))
!$OMP END PARALLEL
      CASE (2)
        dim1 = SIZE(in_var%s_ptr,1)
        dim2 = SIZE(in_var%s_ptr,3)
        dim3 = SIZE(in_var%s_ptr,4)
        ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
        !$ACC ENTER DATA CREATE(tmp_var)
!$OMP PARALLEL
        CALL copy(in_var%s_ptr(:,in_var_idx,:,:,1), tmp_var(:,:,:))
!$OMP END PARALLEL
      CASE (3)
        dim1 = SIZE(in_var%s_ptr,1)
        dim2 = SIZE(in_var%s_ptr,2)
        dim3 = SIZE(in_var%s_ptr,4)
        ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
        !$ACC ENTER DATA CREATE(tmp_var)
!$OMP PARALLEL
        CALL copy(in_var%s_ptr(:,:,in_var_idx,:,1), tmp_var(:,:,:))
!$OMP END PARALLEL
      CASE (4)
        dim1 = SIZE(in_var%s_ptr,1)
        dim2 = SIZE(in_var%s_ptr,2)
        dim3 = SIZE(in_var%s_ptr,3)
        ALLOCATE(tmp_var(dim1, dim2, dim3), STAT=ierrstat)
        IF (ierrstat /= SUCCESS)  CALL finish (routine, 'allocation of tmp_var failed')
        !$ACC ENTER DATA CREATE(tmp_var)
!$OMP PARALLEL
        CALL copy(in_var%s_ptr(:,:,:,in_var_idx,1), tmp_var(:,:,:))
!$OMP END PARALLEL
      CASE default
        CALL finish(routine, "internal error!")
      END SELECT
      in_ptr => tmp_var(:,:,:) 
    ELSE
      CALL finish (routine, 'internal error!')
    ENDIF

    SELECT CASE(out_var_ref_pos)
    CASE (1)
      out_ptr => out_var%r_ptr(out_var_idx,:,:,:,1)
    CASE (2)
      out_ptr => out_var%r_ptr(:,out_var_idx,:,:,1)
    CASE (3)
      out_ptr => out_var%r_ptr(:,:,out_var_idx,:,1)
    CASE (4)
      out_ptr => out_var%r_ptr(:,:,:,out_var_idx,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT

    SELECT CASE ( p_info%hgrid )
    CASE (GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE) 
      ! consistency check:
      IF ((UBOUND(in_ptr,1) > nproma)              .OR.  &
        & (UBOUND(in_ptr,2) > UBOUND(in_ptr, 2))   .OR.  &
        & (UBOUND(in_ptr,3) > nblks)) THEN
        CALL finish(routine, "Inconsistent array dimensions")
      END IF
    CASE DEFAULT
      CALL finish(routine, "Internal error!")
    END SELECT

    !--- actually perform vertical interpolation task
    IF (.NOT. ((nblks == 0) .OR. ((nblks == 1) .AND. (npromz == 0)))) THEN

      !$ACC DATA PRESENT(in_ptr, out_ptr)
    
      SELECT CASE ( vert_intp_method )
      CASE ( VINTP_METHOD_VN )
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_VN")
        IF (.NOT. ASSOCIATED(vcoeff_lin) .OR. .NOT. ASSOCIATED(vcoeff_cub)) &
          CALL finish(routine, "Internal error!")
        CALL uv_intp(in_ptr,                                                        & !in
          &          out_ptr,                                                       & !out
          &          in_z_mc, in_z3d,                                               & !in
          &          nblks, npromz, nlev, n_ipzlev,                                 & !in
          &          vcoeff_cub%coef1, vcoeff_cub%coef2,                            & !in
          &          vcoeff_cub%coef3, vcoeff_lin%wfac_lin,                         & !in
          &          vcoeff_cub%idx0_cub, vcoeff_lin%idx0_lin,                      & !in
          &          vcoeff_cub%bot_idx_cub, vcoeff_lin%bot_idx_lin,                & !in
          &          vcoeff_lin%wfacpbl1, vcoeff_lin%kpbl1,                         & !in
          &          vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl2,                         & !in
          &          l_hires_intp=l_hires_intp,                                     & !in
          &          l_restore_fricred=l_restore_fricred, lacc=.TRUE. )               !in
        !
      CASE ( VINTP_METHOD_LIN )        
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN")
        IF (.NOT. ASSOCIATED(vcoeff_lin)) CALL finish(routine, "Internal error!")
        CALL lin_intp(in_ptr,                                                       & !in
          &           out_ptr,                                                      & !out
          &           nblks, npromz, nlev, n_ipzlev,                                & !in
          &           vcoeff_lin%wfac_lin, vcoeff_lin%idx0_lin,                     & !in
          &           vcoeff_lin%bot_idx_lin,                                       & !in
          &           vcoeff_lin%wfacpbl1, vcoeff_lin%kpbl1,                        & !in
          &           vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl2,                        & !in
          &           l_loglin=l_loglin,                                            & !in
          &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,                 & !in
          &           lower_limit=lower_limit, lacc=.TRUE. )                          !in
        !
      CASE ( VINTP_METHOD_LIN_NLEVP1 )        
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_LIN_NLEVP1")
        IF (.NOT. ASSOCIATED(vcoeff_lin_nlevp1)) CALL finish(routine, "Internal error!")
        CALL lin_intp(in_ptr,                                                       & !in
          &           out_ptr,                                                      & !out
          &           nblks, npromz, nlevp1, n_ipzlev,                              & !in
          &           vcoeff_lin_nlevp1%wfac_lin,                                   & !in
          &           vcoeff_lin_nlevp1%idx0_lin,                                   & !in
          &           vcoeff_lin_nlevp1%bot_idx_lin,                                & !in
          &           vcoeff_lin_nlevp1%wfacpbl1, vcoeff_lin_nlevp1%kpbl1,          & !in
          &           vcoeff_lin_nlevp1%wfacpbl2, vcoeff_lin_nlevp1%kpbl2,          & !in
          &           l_loglin=l_loglin,                                            & !in
          &           l_extrapol=l_extrapol, l_pd_limit=l_pd_limit,                 & !in
          &           lower_limit=lower_limit, lacc=.TRUE. )                          !in
        !
      CASE (VINTP_METHOD_QV )
        IF (dbg_level > 15)  CALL message(routine, "VINTP_METHOD_QV")
        IF (.NOT. ASSOCIATED(vcoeff_lin) .OR. .NOT. ASSOCIATED(vcoeff_cub)) &
          CALL finish(routine, "Internal error!")
        CALL qv_intp(in_ptr,                                                        & !in
          &          out_ptr,                                                       & !out
          &          in_z_mc, in_z3d, p_diag%temp,                                  & !in
          &          p_diag%pres, p_temp, p_pres,                                   & !in
          &          nblks, npromz, nlev, n_ipzlev,                                 & !in
          &          vcoeff_cub%coef1, vcoeff_cub%coef2,                            & !in
          &          vcoeff_cub%coef3,                                              & !in
          &          vcoeff_lin%wfac_lin, vcoeff_cub%idx0_cub,                      & !in
          &          vcoeff_lin%idx0_lin,                                           & !in
          &          vcoeff_cub%bot_idx_cub, vcoeff_lin%bot_idx_lin,                & !in
          &          vcoeff_lin%wfacpbl1, vcoeff_lin%kpbl1,                         & !in
          &          vcoeff_lin%wfacpbl2, vcoeff_lin%kpbl2,                         & !in
          &          l_satlimit=l_satlimit, lower_limit=lower_limit,                & !in
          &          l_restore_pbldev=l_restore_pbldev, lacc=.TRUE. )                 !in
      END SELECT ! vert_intp_method

      !$ACC END DATA

   END IF

    ! clean up
    IF (p_info%hgrid == GRID_UNSTRUCTURED_EDGE) THEN
      !$ACC WAIT
      !$ACC EXIT DATA DELETE(p_z3d_edge, z_me)
      DEALLOCATE(p_z3d_edge, z_me, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'DEALLOCATE failed')
    END IF
    IF (ALLOCATED(tmp_var)) THEN
      !$ACC WAIT
      !$ACC EXIT DATA DELETE(tmp_var)
      DEALLOCATE(tmp_var, STAT=ierrstat)
      IF (ierrstat /= SUCCESS)  CALL finish (routine, 'deallocation of tmp_var failed')
    ENDIF

  END SUBROUTINE pp_task_ipzlev


  !---------------------------------------------------------------
  !> Performs interpolation of a 2D field onto mean sea level, z=0.
  !
  !  This routine is completely independent from the data structures
  !  used for pz-level interpolation.
  !
  SUBROUTINE pp_task_intp_msl(ptr_task, lacc)
    TYPE(t_job_queue), POINTER, INTENT(INOUT) :: ptr_task
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc
    ! local variables    
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_intp_msl"
    INTEGER,  PARAMETER :: nzlev         =        1     ! just a single z-level... 
    REAL(wp), PARAMETER :: ZERO_HEIGHT   =    0._wp, &
      &                    EXTRAPOL_DIST = -500._wp

    INTEGER                            :: nblks_c, npromz_c, nblks_e, jg, jc, jb, &
      &                                   out_var_idx, nlev, i_endblk
    TYPE (t_var), POINTER :: in_var, out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_nh_metrics),        POINTER :: p_metrics    
    TYPE(t_nh_diag),           POINTER :: p_diag

    REAL(wp) :: pmsl_aux(nproma,1,ptr_task%data_input%p_patch%nblks_c), &
                pmsl_avg(nproma,1,ptr_task%data_input%p_patch%nblks_c)

    INTEGER,  DIMENSION(nproma, ptr_task%data_input%p_patch%nblks_c) :: &
      &  kpbl1, kpbl2
    REAL(wp), DIMENSION(nproma, ptr_task%data_input%p_patch%nblks_c) :: &
      &  zextrap, wfacpbl1, wfacpbl2

    CALL assert_acc_device_only(routine, lacc)
    CALL assert_lacc_equals_i_am_accel_node(routine, lacc)

    ! patch, state, and metrics
    jg             =  ptr_task%data_input%jg
    p_patch        => ptr_task%data_input%p_patch
    p_metrics      => ptr_task%data_input%p_nh_state%metrics
    p_diag         => ptr_task%data_input%p_nh_state%diag

    ! input/output field for this task
    p_info            => ptr_task%data_input%var%info
    in_var            => ptr_task%data_input%var
    out_var           => ptr_task%data_output%var

    IF (TRIM(p_info%name) /= "pres")  CALL message(routine, "Invalid input field!")

    nlev     = p_patch%nlev
    nblks_c  = p_patch%nblks_c
    npromz_c = p_patch%npromz_c

    out_var_idx = 1
    IF (out_var%info%lcontained) out_var_idx = out_var%info%ncontained

    !$ACC DATA CREATE(pmsl_aux, pmsl_avg)

    SELECT CASE (itype_pres_msl)
    CASE (PRES_MSL_METHOD_SAI) ! stepwise analytical integration 

      IF (dbg_level >= 10)  CALL message(routine, "PRES_MSL_METHOD_SAI: stepwise analytical integration")

      !$ACC DATA CREATE(kpbl1, wfacpbl1, kpbl2, wfacpbl2)
      CALL init(kpbl1, opt_acc_async=.TRUE.)
      CALL init(wfacpbl1, opt_acc_async=.TRUE.)
      CALL init(kpbl2, opt_acc_async=.TRUE.)
      CALL init(wfacpbl2, opt_acc_async=.TRUE.)

      ! compute extrapolation coefficients:
      CALL prepare_extrap(p_metrics%z_mc,                                     & !in
        &                 nblks_c, npromz_c, nlev,                            & !in
        &                 kpbl1, wfacpbl1, kpbl2, wfacpbl2,                   & !out
        &                 lacc=.TRUE.                                         & !in
        & )

      ! Interpolate pressure on z-level "0": 
      CALL diagnose_pmsl(p_diag%pres, p_diag%tempv, p_metrics%z_mc,           &
        &                pmsl_aux(:,1,:),                                     &
        &                nblks_c, npromz_c, p_patch%nlev,                     &
        &                wfacpbl1, kpbl1, wfacpbl2, kpbl2,                    &
        &                ZERO_HEIGHT, EXTRAPOL_DIST,                          &
        &                lacc=.TRUE.)
      !$ACC WAIT
      !$ACC END DATA

    CASE (PRES_MSL_METHOD_GME) ! GME-type extrapolation
      CALL assert_acc_host_only(routine//':PRES_MSL_METHOD_GME', lacc) ! ACC not tested

      IF (dbg_level >= 10)  CALL message(routine, "PRES_MSL_METHOD_GME")
      ! Interpolate pressure on z-level "0":

      CALL diagnose_pmsl_gme(p_diag%pres, p_diag%pres_sfc, p_diag%temp, &  ! in
        &                    p_metrics%z_ifc,                           &  ! in
        &                    pmsl_aux(:,1,:),                           &  ! out
        &                    nblks_c, npromz_c, p_patch%nlev,           &  ! in
        &                    lacc=.FALSE.)                                 ! in

    CASE (PRES_MSL_METHOD_IFS,PRES_MSL_METHOD_IFS_CORR,PRES_MSL_METHOD_DWD) ! IFS or new DWD extrapolation method

      IF (dbg_level >= 10)  THEN
        IF (itype_pres_msl == PRES_MSL_METHOD_DWD) THEN
          CALL message(routine, "PRES_MSL_METHOD_DWD")
        ELSE
          CALL message(routine, "PRES_MSL_METHOD_IFS")
        ENDIF
      ENDIF
      !$ACC DATA CREATE(kpbl1, wfacpbl1, zextrap)
      CALL init(kpbl1, opt_acc_async=.TRUE.)
      CALL init(wfacpbl1, opt_acc_async=.TRUE.)
      CALL init(zextrap, opt_acc_async=.TRUE.)
      ! compute extrapolation coefficients:
      CALL prepare_extrap_ifspp(p_metrics%z_ifc, p_metrics%z_mc,              & !in
        &                 nblks_c, npromz_c, nlev,                            & !in
        &                 kpbl1, zextrap, wfacpbl1,                           & !out
        &                 lacc=.TRUE.)

      ! Interpolate pressure on z-level "0":
      CALL diagnose_pmsl_ifs(p_diag%pres_sfc, p_diag%temp, p_metrics%z_ifc,   & ! in
        &                    pmsl_aux(:,1,:),                                 & ! out
        &                    nblks_c, npromz_c, p_patch%nlev,                 & ! in
        &                    wfacpbl1, kpbl1, zextrap, itype_pres_msl,        & ! in
        &                    lacc=.TRUE.)                                       ! in
      !$ACC WAIT
      !$ACC END DATA

    CASE DEFAULT
      CALL finish(routine, 'Internal error!')
    END SELECT

    IF (l_limited_area .OR. jg > 1) THEN ! copy outermost nest boundary row in order to avoid missing values
      i_endblk = p_patch%cells%end_blk(1,1)
      !$OMP PARALLEL
      CALL copy(pmsl_aux(:,1,1:i_endblk), pmsl_avg(:,1,1:i_endblk), opt_acc_async=.TRUE.)
      !$OMP END PARALLEL
    ENDIF

    CALL cell_avg(pmsl_aux, p_patch, p_int_state(jg)%c_bln_avg, pmsl_avg, lacc=.TRUE.)
    i_endblk = ptr_task%data_input%p_patch%nblks_c
    !$OMP PARALLEL
    CALL copy(pmsl_avg(:,1,1:i_endblk), out_var%r_ptr(:,1:i_endblk,out_var_idx,1,1), opt_acc_async=.TRUE.)
    !$OMP END PARALLEL
    !$ACC WAIT
    !$ACC END DATA

  END SUBROUTINE pp_task_intp_msl


  !---------------------------------------------------------------
  !> Performs computation of optional diagnostic fields.
  !
  !  Selects subroutines for field computation based on the variable
  !  name.
  !
  !  @note This could be easily replaced by a procedure pointer,
  !        alas, this is an F2003 feature.
  !
  !  @todo Change order of processing: First, interpolate input fields
  !        onto z-levels, then compute rel_hum.
  !
  SUBROUTINE pp_task_compute_field(ptr_task)

    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    INTEGER                            :: jg, out_var_idx
    TYPE (t_var), POINTER :: out_var
    TYPE(t_var_metadata),      POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch
    !TYPE(t_int_state),         POINTER :: p_int_state
    TYPE(t_nh_prog),           POINTER :: p_prog, p_prog_rcf
    TYPE(t_nh_diag),           POINTER :: p_diag
    TYPE(t_nwp_phy_diag),      POINTER :: prm_diag
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_compute_field"
    LOGICAL :: lclip                   ! limit rh to MAX(rh,100._wp)
    
    ! output field for this task
    out_var   => ptr_task%data_output%var
    p_info    => out_var%info    
    out_var_idx = 1
    IF (out_var%info%lcontained)  out_var_idx = out_var%info%ncontained

    ! input data required for computation:
    jg          =  ptr_task%data_input%jg
    p_patch     => ptr_task%data_input%p_patch
    !p_int_state => ptr_task%data_input%p_int_state
    p_prog      => ptr_task%data_input%p_nh_state%prog(nnow(jg))
    p_prog_rcf  => ptr_task%data_input%p_nh_state%prog(nnow_rcf(jg))
    p_diag      => ptr_task%data_input%p_nh_state%diag
    prm_diag    => ptr_task%data_input%prm_diag

    SELECT CASE(ptr_task%job_type)
    CASE (TASK_COMPUTE_RH)

      SELECT CASE (itype_rh)
      CASE (RH_METHOD_WMO)
        CALL compute_field_rel_hum_wmo(p_patch, p_prog, p_diag, &
          &                        out_var%r_ptr(:,:,:,out_var_idx,1), lacc=i_am_accel_node)
      CASE (RH_METHOD_IFS, RH_METHOD_IFS_CLIP)
        IF (itype_rh == RH_METHOD_IFS_CLIP) THEN
          lclip = .TRUE.
        ELSE
          lclip = .FALSE.
        ENDIF
#ifdef _OPENACC
        CALL finish(routine, 'not yet ported postproc RH_METHOD_IFS, RH_METHOD_IFS_CLIP for variable '//TRIM(p_info%name) )
#endif
        CALL compute_field_rel_hum_ifs(p_patch, p_prog, p_diag,        &
          &                        out_var%r_ptr(:,:,:,out_var_idx,1), &
          &                        opt_lclip=lclip)

      CASE DEFAULT
        CALL finish(routine, 'Internal error!')
      END SELECT

    CASE (TASK_COMPUTE_OMEGA)
      CALL compute_field_omega(p_patch, p_prog, &
        &                      out_var%r_ptr(:,:,:,out_var_idx,1), lacc=i_am_accel_node)
    
    CASE (TASK_COMPUTE_PV)
      CALL compute_field_pv(p_patch, p_int_state(jg),                  &
        &   ptr_task%data_input%p_nh_state%metrics, p_prog, p_diag,    &  
        &   out_var%r_ptr(:,:,:,out_var_idx,1), lacc=i_am_accel_node)

    CASE (TASK_COMPUTE_SDI2)
      IF ( jg >= n_dom_start+1 ) THEN
        ! p_patch_local_parent(jg) seems to exist
        CALL compute_field_sdi( p_patch, jg, p_patch_local_parent(jg), p_int_state_local_parent(jg),     &
          &   ptr_task%data_input%p_nh_state%metrics, p_prog, p_diag,    &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1
      ELSE
        CALL message( routine, "WARNING: SDI2 cannot be computed since no reduced grid is available" )
      END IF

    CASE (TASK_COMPUTE_LPI)
      IF ( jg >= n_dom_start+1 ) THEN
        ! p_patch_local_parent(jg) seems to exist
        CALL compute_field_lpi( p_patch, jg, p_patch_local_parent(jg), p_int_state_local_parent(jg),     &
          &   ptr_task%data_input%p_nh_state%metrics, p_prog, p_prog_rcf, p_diag,    &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1
      ELSE
        CALL message( routine, "WARNING: LPI cannot be computed since no reduced grid is available" )
      END IF

    CASE (TASK_COMPUTE_CEILING)
      CALL compute_field_ceiling( p_patch, jg,                                       &
          &   ptr_task%data_input%p_nh_state%metrics, prm_diag,                      &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_VIS)
      CALL compute_field_visibility( p_patch, p_prog, p_diag, prm_diag, jg,          &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_INVERSION)
      CALL compute_field_inversion_height( p_patch, jg, ptr_task%data_input%p_nh_state%metrics, p_prog, p_diag,prm_diag,   &
          &   out_var%r_ptr(:,:,out_var_idx,1,1))   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_HBAS_SC)
      CALL compute_field_hbas_sc( p_patch,                                           &
          &   ptr_task%data_input%p_nh_state%metrics, prm_diag,                      &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_HTOP_SC)
      CALL compute_field_htop_sc( p_patch,                                           &
          &   ptr_task%data_input%p_nh_state%metrics, prm_diag,                      &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_TWATER)
      CALL compute_field_twater( p_patch, ptr_task%data_input%p_nh_state%metrics%ddqz_z_full, &
          &                      p_prog%rho, p_prog_rcf%tracer,                               &
          &                      advection_config(jg)%trHydroMass%list,                       &
          &                      out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)

    CASE (TASK_COMPUTE_Q_SEDIM)
      CALL compute_field_q_sedim( p_patch, jg, p_prog,                               &
          &   out_var%r_ptr(:,:,:,out_var_idx,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_DBZ850)
      CALL compute_field_dbz850( p_patch, prm_diag%k850(:,:), prm_diag%dbz3d_lin(:,:,:), &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_DBZLMX_LOW)
      ! NOTE: The layer bounds 1000 m and 2000 m were found more appropriate than the fixed bounds
      !       500 m and 2500 m in the eccodes definition of DBZLMX_LOW. Because of possible
      !       further adaptions in the near future and to avoid several consecutive adaptions of eccodes
      !       until consolidation of the layer bounds, we for now set the bounds to 1000 and 2000 here without
      !       changing the fixed bounds in the eccodes definitions.
      CALL compute_field_dbzlmx( p_patch, jg, 1000.0_wp, 2000.0_wp, &
          &   ptr_task%data_input%p_nh_state%metrics, prm_diag%dbz3d_lin(:,:,:), &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_DBZCMAX)
      CALL compute_field_dbzcmax( p_patch, jg, prm_diag%dbz3d_lin(:,:,:),            &
          &   out_var%r_ptr(:,:,out_var_idx,1,1), lacc=i_am_accel_node)   ! unused dimensions are filled up with 1

    CASE (TASK_COMPUTE_SMI)
      CALL compute_field_smi(p_patch, p_lnd_state(jg)%diag_lnd, &
           &                 ext_data(jg), out_var%r_ptr(:,:,:,out_var_idx,1), lacc=i_am_accel_node)

    CASE (TASK_COMPUTE_WSHEAR_U)
#ifdef _OPENACC
      CALL finish(routine, 'not yet ported postproc TASK_COMPUTE_WSHEAR_U for variable '//TRIM(p_info%name) )
#endif
      CALL compute_field_wshear( p_patch, ptr_task%data_input%p_nh_state%metrics, &
           &                     p_diag%u, wshear_uv_heights(1:n_wshear), out_var%r_ptr(:,:,:,out_var_idx,1) )
      
    CASE (TASK_COMPUTE_WSHEAR_V)
#ifdef _OPENACC
      CALL finish(routine, 'not yet ported postproc TASK_COMPUTE_WSHEAR_V for variable '//TRIM(p_info%name) )
#endif
      CALL compute_field_wshear( p_patch, ptr_task%data_input%p_nh_state%metrics, &
           &                     p_diag%v, wshear_uv_heights(1:n_wshear), out_var%r_ptr(:,:,:,out_var_idx,1) )

    CASE (TASK_COMPUTE_LAPSERATE)
#ifdef _OPENACC
      CALL finish(routine, 'not yet ported postproc TASK_COMPUTE_LAPSERATE for variable '//TRIM(p_info%name) )
#endif
      CALL compute_field_lapserate( p_patch, ptr_task%data_input%p_nh_state%metrics, &
           &                        p_diag, 500e2_wp, 850e2_wp, out_var%r_ptr(:,:,out_var_idx,1,1) )

    CASE (TASK_COMPUTE_SRH)
#ifdef _OPENACC
      CALL finish(routine, 'not yet ported postproc TASK_COMPUTE_SRH for variable '//TRIM(p_info%name) )
#endif
      CALL compute_field_srh( ptr_patch     = p_patch,   &
           &                  p_metrics     = ptr_task%data_input%p_nh_state%metrics, &
           &                  p_diag        = p_diag,    &
           &                  z_up_srh      = srh_heights(1:n_srh), &
           &                  z_up_meanwind = 6000.0_wp, &
           &                  z_low_shear   = 250.0_wp,  &
           &                  z_up_shear    = 5750.0_wp, &
           &                  dz_shear      = 500.0_wp,  &
           &                  srh           = out_var%r_ptr(:,:,:,out_var_idx,1) )
     
    CASE DEFAULT
      CALL finish(routine, 'Internal error!')
    END SELECT

  END SUBROUTINE pp_task_compute_field


  !---------------------------------------------------------------
  !> Performs interpolation of a edge-based variable onto cell centers.
  !
  !  This is only a wrapper for the corresponding routines from the
  !  interpolation module.
  !  This routine should be GPU-capable as rbf_vec_interpol_cell
  !  has been ported with OpenACC
  SUBROUTINE pp_task_edge2cell(ptr_task)
    TYPE(t_job_queue), POINTER :: ptr_task
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::pp_task_edge2cell"
    INTEGER :: &
      &  in_var_idx, out_var_idx_1, out_var_idx_2, &
      &  in_var_ref_pos, out_var_ref_pos_1,        &
      &  out_var_ref_pos_2
    TYPE (t_var), POINTER :: in_var, out_var_1, out_var_2
    TYPE (t_var_metadata),     POINTER :: p_info
    TYPE(t_patch),             POINTER :: p_patch
    TYPE(t_int_state),         POINTER :: intp_hrz
    REAL(wp),                  POINTER :: in_ptr(:,:,:), out_ptr_1(:,:,:), &
      &                                   out_ptr_2(:,:,:)

    p_patch        => ptr_task%data_input%p_patch      ! patch
    intp_hrz       => ptr_task%data_input%p_int_state
    in_var         => ptr_task%data_input%var
    p_info         => ptr_task%data_input%var%info

    ! Consistency check: We make the following assumptions:
    ! - This is a 3D variable
    IF (zaxisTypeList%is_2d(p_info%vgrid))  CALL finish(routine, "Internal error!")
    ! - We have two output components:
    IF (.NOT. ASSOCIATED(ptr_task%data_output%var) .OR.  &
      & .NOT. ASSOCIATED(ptr_task%data_output%var_2)) THEN
      CALL finish(routine, "Internal error!")
    END IF

    in_var_idx        = 1
    in_var_ref_pos    = 4
    IF (in_var%info%lcontained) THEN
      in_var_idx     = in_var%info%ncontained
      in_var_ref_pos = in_var%info%var_ref_pos
    END IF

    out_var_1         => ptr_task%data_output%var
    out_var_idx_1     = 1
    out_var_ref_pos_1 = 4
    IF (out_var_1%info%lcontained) THEN
      out_var_idx_1     = out_var_1%info%ncontained
      out_var_ref_pos_1 = out_var_1%info%var_ref_pos
    END IF
    out_var_2         => ptr_task%data_output%var_2
    out_var_idx_2     =  1
    out_var_ref_pos_2 = 4
    IF (out_var_2%info%lcontained) THEN
      out_var_idx_2     = out_var_2%info%ncontained
      out_var_ref_pos_2 = out_var_2%info%var_ref_pos
    END IF

    ! throw error message, if this variable is not a REAL field:
    IF (.NOT. ASSOCIATED(in_var%r_ptr)) THEN
      CALL finish(routine, TRIM(p_info%name)//": Implemented for REAL fields only.")
    END IF

    SELECT CASE(in_var_ref_pos)
    CASE (1)
      in_ptr => in_var%r_ptr(in_var_idx,:,:,:,1)
    CASE (2)
      in_ptr => in_var%r_ptr(:,in_var_idx,:,:,1)
    CASE (3)
      in_ptr => in_var%r_ptr(:,:,in_var_idx,:,1)
    CASE (4)
      in_ptr => in_var%r_ptr(:,:,:,in_var_idx,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    SELECT CASE(out_var_ref_pos_1)
    CASE (1)
      out_ptr_1 => out_var_1%r_ptr(out_var_idx_1,:,:,:,1)
    CASE (2)
      out_ptr_1 => out_var_1%r_ptr(:,out_var_idx_1,:,:,1)
    CASE (3)
      out_ptr_1 => out_var_1%r_ptr(:,:,out_var_idx_1,:,1)
    CASE (4)
      out_ptr_1 => out_var_1%r_ptr(:,:,:,out_var_idx_1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    SELECT CASE(out_var_ref_pos_2)
    CASE (1)
      out_ptr_2 => out_var_2%r_ptr(out_var_idx_2,:,:,:,1)
    CASE (2)
      out_ptr_2 => out_var_2%r_ptr(:,out_var_idx_2,:,:,1)
    CASE (3)
      out_ptr_2 => out_var_2%r_ptr(:,:,out_var_idx_2,:,1)
    CASE (4)
      out_ptr_2 => out_var_2%r_ptr(:,:,:,out_var_idx_2,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT

    CALL rbf_vec_interpol_cell(in_ptr,                                  &   !< normal wind comp.
      &                        p_patch, intp_hrz,                       &   !< patch, interpolation state
      &                        out_ptr_1, out_ptr_2 )                       !< reconstr. u,v wind

  END SUBROUTINE pp_task_edge2cell

END MODULE mo_pp_tasks
