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

MODULE mo_output_coupling

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_var                 ,ONLY: t_var_ptr
  USE mo_var_groups          ,ONLY: MAX_GROUPS, var_groups_dyn
  USE mo_run_config          ,ONLY: nlev, msg_level
  USE mo_run_config          ,ONLY: ltimer
  USE mo_timer               ,ONLY: timer_start, timer_stop, &
       &                            timer_coupling_put, timer_coupling_init
  USE mo_util_string         ,ONLY: int2string
  USE mo_exception           ,ONLY: message, finish
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_zaxis_type          ,ONLY: zaxisTypeList

#ifdef _OPENACC
  USE mo_mpi                 ,ONLY: i_am_accel_node
  USE openacc
#endif

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_output_coupling' ! Output of module for debug
  CHARACTER(len=2), PARAMETER :: newline = ACHAR(13) // ACHAR(10)

  TYPE, PRIVATE :: t_exposed_var
     INTEGER :: yac_field_id
     TYPE(t_var_ptr) :: var(1:3)
     INTEGER :: tlev_source, var_size
     TYPE(t_exposed_var), POINTER :: next => NULL()
  END TYPE t_exposed_var

  PUBLIC :: construct_output_coupling
  PUBLIC :: winnow_field_list
  PUBLIC :: output_coupling
  PUBLIC :: destruct_output_coupling

  TYPE(t_exposed_var), POINTER :: exposed_vars_head => NULL()
  INTEGER :: max_collection_size = 0, max_hor_size = 0


CONTAINS

  !>
  !! SUBROUTINE construct_output_coupling -- the initialisation for
  !! the coupling of atmosphere and output components This routine
  !! iterates over all variables in all variablelists and defines proper
  !! variables as a fields in the coupler.

  SUBROUTINE construct_output_coupling ( &
    p_patch, comp_id, cell_point_ids, vertex_point_ids, timestepstring)

    USE mo_var_list_register,   ONLY: t_vl_register_iter
    USE mo_var_metadata,        ONLY: get_var_timelevel, get_var_name
    USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT
    USE mo_yac_finterface      ,ONLY: yac_fdef_field, YAC_TIME_UNIT_ISO_FORMAT, &
         yac_fdef_field_metadata, yac_fget_component_name, yac_fget_grid_name

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_ids(1), vertex_point_ids(1)
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    TYPE(t_vl_register_iter), ALLOCATABLE :: vl_iter
    TYPE(t_exposed_var), POINTER :: exposed_var
    CHARACTER(len=:), ALLOCATABLE :: var_name, metadata, comp_name, grid_name
    INTEGER :: iv, tl, collection_size, key_notl, count = 0, var_size, grpi, nblks, pos(3)
    INTEGER :: point_ids(1), var_ref_pos

    TYPE t_tmp_timelevel_var
       INTEGER :: key_notl
       TYPE(t_exposed_var), POINTER :: exposed_var
       TYPE(t_tmp_timelevel_var), POINTER :: next => NULL()
    END type t_tmp_timelevel_var

    TYPE(t_tmp_timelevel_var), POINTER :: tmp_timelevel_var_head => NULL(), tmp_timelevel_var => NULL()

    max_hor_size = MAX(p_patch(1)%n_patch_cells, p_patch(1)%n_patch_verts)

    ALLOCATE(vl_iter)
    VARLIST_LOOP: DO WHILE(vl_iter%next())
       IF (vl_iter%cur%p%patch_id .NE. 1) CYCLE
       DO iv = 1, vl_iter%cur%p%nvars

          ASSOCIATE( elem => vl_iter%cur%p%vl(iv)%p )

            IF (.NOT. elem%info%loutput) THEN
               IF (msg_level >= 15) &
                  CALL message(str_module, "Omitted due to missing output flag: " // TRIM(elem%info%name))
               CYCLE
            END IF

            ! donot expose container variables
            IF (elem%info%lcontainer) THEN
               IF (msg_level >= 15) &
                  CALL message(str_module, "Omitted due to lcontainer: " // TRIM(elem%info%name))
               CYCLE
            END IF

            ! expose only vars on the icon grid (cells or vertices)
            IF (vl_iter%cur%p%hgrid(iv) .EQ. GRID_UNSTRUCTURED_CELL) THEN
               point_ids = cell_point_ids
               var_size = p_patch(1)%n_patch_cells
               nblks = p_patch(1)%nblks_c
            ELSEIF (vl_iter%cur%p%hgrid(iv) .EQ. GRID_UNSTRUCTURED_VERT) THEN
               point_ids = vertex_point_ids
               var_size = p_patch(1)%n_patch_verts
               nblks = p_patch(1)%nblks_v
            ELSE
               IF (msg_level >= 15) &
                  CALL message(str_module, "Omitted due to invalid grid: " // TRIM(elem%info%name))
               CYCLE
            END IF

            var_ref_pos = MERGE(elem%info%var_ref_pos, 4, elem%info%lcontained)
            IF (zaxisTypeList%is_2d(elem%info%vgrid)) THEN
               pos = (/ 1, -99, 2 /)
            ELSE
               pos = (/ 1, 2, 3 /)
            END IF
            pos = MERGE(pos, pos+1, pos <= var_ref_pos)

            ! ensure correct dimension
            IF (elem%info%used_dimensions(pos(1)) .NE. nproma .OR. &
               elem%info%used_dimensions(pos(3)) .LT. nblks)  THEN !the blk dimension could be larger than nblks
               IF (msg_level >= 15) &
                  CALL message(str_module, "Omitted due to invalid dimensions: " // TRIM(elem%info%name) // &
                       "nproma: " // int2string(elem%info%used_dimensions(pos(1))) // ", " // &
                       "nblks: " // int2string(elem%info%used_dimensions(pos(3))))
               CYCLE
            END IF

            IF (zaxisTypeList%is_2d(elem%info%vgrid)) THEN
               collection_size = 1
            ELSE
               collection_size = elem%info%used_dimensions(pos(2))
            END IF

            tl = get_var_timelevel(elem%info%name)
            key_notl = vl_iter%cur%p%key_notl(iv)
            var_name = get_var_name(elem%info)
            exposed_var => NULL()
            IF(tl /= -1) THEN
               ! check if we already have a timelevel val registered
               tmp_timelevel_var => tmp_timelevel_var_head
               TL_VAR_LOOP: DO WHILE(ASSOCIATED(tmp_timelevel_var))
                  IF (tmp_timelevel_var%key_notl == key_notl) THEN
                     exposed_var => tmp_timelevel_var%exposed_var
                     EXIT TL_VAR_LOOP
                  ENDIF
                  tmp_timelevel_var => tmp_timelevel_var%next
               END DO TL_VAR_LOOP
            END IF
            IF (.NOT. ASSOCIATED(exposed_var)) THEN
               ALLOCATE(exposed_var)
               exposed_var%var_size = var_size
               IF(tl /= -1) THEN
                  exposed_var%var(tl) = vl_iter%cur%p%vl(iv)
                  exposed_var%tlev_source = elem%info%tlev_source
               ELSE
                  exposed_var%var(1) = vl_iter%cur%p%vl(iv)
                  exposed_var%tlev_source = -1
               ENDIF
               CALL yac_fdef_field(             &
                    & TRIM(var_name),           &
                    & comp_id,                  &
                    & point_ids,                &
                    & 1,                        &
                    & collection_size,          &
                    & timestepstring,           &
                    & YAC_TIME_UNIT_ISO_FORMAT, &
                    & exposed_var%yac_field_id )
               count = count + 1
               max_collection_size = MAX(max_collection_size, collection_size)
               ! add element to list

               exposed_var%next => exposed_vars_head
               exposed_vars_head => exposed_var
               ! if it has a nontrival timelevel we add it into the temporary list
               IF (tl /= -1) THEN
                  ALLOCATE(tmp_timelevel_var)
                  tmp_timelevel_var%key_notl = key_notl
                  tmp_timelevel_var%exposed_var => exposed_var
                  tmp_timelevel_var%next => tmp_timelevel_var_head
                  tmp_timelevel_var_head => tmp_timelevel_var
               END IF

               ! add some metadata
               metadata = "standard_name: !!str '" // TRIM(elem%info%cf%standard_name) // "'" // newline &
                    //    "units: !!str '" // TRIM(elem%info%cf%units) // "'" // newline &
                    //    "short_name: !!str '" // TRIM(elem%info%cf%short_name) // "'" // newline &
                    //    "long_name: !!str '" // TRIM(elem%info%cf%long_name) // "'" // newline &
                    //    "groups:" // newline
               DO grpi = 1,SIZE(var_groups_dyn%gname)
                  IF (elem%info%in_group(grpi)) metadata = metadata // "  - " // TRIM(var_groups_dyn%gname(grpi)) // newline
               END DO
               CALL yac_fdef_field_metadata( &
                    yac_fget_component_name(exposed_var%yac_field_id),&
                    yac_fget_grid_name(exposed_var%yac_field_id), &
                    var_name, metadata)
            ELSE
               exposed_var%var(tl) = vl_iter%cur%p%vl(iv)
            END IF

          END ASSOCIATE

       ENDDO
    ENDDO VARLIST_LOOP
    comp_name = yac_fget_component_name(exposed_vars_head%yac_field_id)
    grid_name = yac_fget_grid_name(exposed_vars_head%yac_field_id)
    IF (msg_level >= 15) &
       CALL message(str_module, int2string(count) // &
                  " variables exposed on " // comp_name // " / " // grid_name)
    tmp_timelevel_var => tmp_timelevel_var_head
    DO WHILE(ASSOCIATED(tmp_timelevel_var))
       tmp_timelevel_var_head => tmp_timelevel_var
       tmp_timelevel_var => tmp_timelevel_var%next
       DEALLOCATE(tmp_timelevel_var_head)
    END DO
    DEALLOCATE(vl_iter)
  END SUBROUTINE construct_output_coupling


  !>
  !! SUBROUTINE winnow_field_list -- sort out all non-coupled fields
  !! from the field_list

  SUBROUTINE winnow_field_list()
    USE mo_yac_finterface, ONLY: yac_fget_role_from_field_id, &
         YAC_EXCHANGE_TYPE_NONE, YAC_EXCHANGE_TYPE_SOURCE
    TYPE(t_exposed_var), POINTER :: exposed_var, tmp
    INTEGER :: role, count = 1

    ! find new head and delete
    role = yac_fget_role_from_field_id(exposed_vars_head%yac_field_id)
    DO WHILE(role /= YAC_EXCHANGE_TYPE_SOURCE)
       IF (role /= YAC_EXCHANGE_TYPE_NONE) CALL finish(str_module, "invalid role for exposed variable")
       tmp => exposed_vars_head
       exposed_vars_head => exposed_vars_head%next
       DEALLOCATE(tmp)
       IF (ASSOCIATED(exposed_vars_head)) THEN
          role = yac_fget_role_from_field_id(exposed_vars_head%yac_field_id)
       ELSE
          count = 0
          IF (msg_level >= 15) &
             CALL message(str_module, "No exposed fields are coupled")
          RETURN
       END IF
    END DO

    ! remove other fields without role
    exposed_var => exposed_vars_head
    DO WHILE(ASSOCIATED(exposed_var%next))
       role = yac_fget_role_from_field_id(exposed_var%next%yac_field_id)
       IF (role == YAC_EXCHANGE_TYPE_SOURCE) THEN
          count = count + 1
          exposed_var => exposed_var%next
       ELSE
          IF (role /= YAC_EXCHANGE_TYPE_NONE) &
               CALL finish(str_module, "invalid role for exposed variable")
          tmp => exposed_var%next
          exposed_var%next => exposed_var%next%next
          DEALLOCATE(tmp)
       END IF
    END DO

    IF (msg_level >= 15) &
       CALL message(str_module, int2string(count) // " exposed vars left after clean up")
  END SUBROUTINE winnow_field_list

  !>
  !! SUBROUTINE output_coupling -- Exchange fields between
  !! atmosphere and output components.
  SUBROUTINE output_coupling (valid_mask)

    USE, INTRINSIC :: ieee_arithmetic
    USE mo_impl_constants      ,ONLY: TLEV_NNOW, TLEV_NNEW, TLEV_NNOW_RCF, TLEV_NNEW_RCF
    USE mo_dynamics_config,     ONLY: nnow, nnow_rcf, nnew, nnew_rcf
    USE mo_yac_finterface,      ONLY: yac_fget_field_collection_size, yac_fput, yac_fget_action, &
      &                               yac_fupdate, YAC_ACTION_NONE, yac_dble_ptr

    INTEGER               :: info, ierror, collection_size, nn, now, ncontained, var_size, var_ref_pos
    REAL(wp), ALLOCATABLE, TARGET, SAVE :: buffer(:,:)
    REAL(wp), CONTIGUOUS, POINTER :: tmp_buffer(:,:)
    TYPE(yac_dble_ptr), ALLOCATABLE :: buffer_ptr(:, :)
    REAL(wp), OPTIONAL :: valid_mask(:,:,:)
    TYPE(t_exposed_var), POINTER :: cur_field
    TYPE(t_var_ptr) :: var_now

    cur_field => exposed_vars_head
    IF (.NOT. ALLOCATED(buffer)) ALLOCATE(buffer(max_hor_size, max_collection_size))
    IF (.NOT. ALLOCATED(buffer_ptr)) ALLOCATE(buffer_ptr(1, max_collection_size))

    DO WHILE (ASSOCIATED(cur_field))
       IF (cur_field%tlev_source == -1) THEN
          now = 1
       ELSE
          SELECT CASE (cur_field%tlev_source)
          CASE(TLEV_NNOW);     now = nnow(1)
          CASE(TLEV_NNOW_RCF); now = nnow_rcf(1)
          CASE(TLEV_NNEW);     now = nnew(1)
          CASE(TLEV_NNEW_RCF); now = nnew_rcf(1)
          CASE DEFAULT
             CALL finish(str_module,'Unsupported tlev_source')
          END SELECT
       ENDIF
       var_now = cur_field%var(now)

       collection_size = yac_fget_field_collection_size(cur_field%yac_field_id)

       CALL yac_fget_action(cur_field%yac_field_id, info)
       IF ( info == YAC_ACTION_NONE ) THEN
          CALL yac_fupdate(cur_field%yac_field_id)
          cur_field => cur_field%next
          IF (msg_level >= 15) &
             CALL message(str_module, " skipping field " // TRIM(var_now%p%info%name))
          CYCLE
       ENDIF
       IF (msg_level >= 15) &
          CALL message(str_module, " sending field " // TRIM(var_now%p%info%name))

!$ACC UPDATE HOST(var_now%p%r_ptr) IF(i_am_accel_node .AND. acc_is_present(var_now%p%r_ptr))
       var_ref_pos = MERGE(var_now%p%info%var_ref_pos, 4, var_now%p%info%lcontained)
       ncontained = MERGE(var_now%p%info%ncontained, 1, var_now%p%info%lcontained)
       var_size = cur_field%var_size

       IF (zaxisTypeList%is_2d(var_now%p%info%vgrid)) THEN
          SELECT CASE (var_ref_pos)
          CASE (1)
             buffer(1:var_size,1:1) = var_now%p%r_ptr(ncontained, :, :, 1, 1)
             buffer_ptr(1, 1)%p => buffer(:,1)
          CASE (2)
             buffer(1:var_size,1:1) = var_now%p%r_ptr(:, ncontained, :, 1, 1)
             buffer_ptr(1, 1)%p => buffer(:,1)
          CASE (3)
             tmp_buffer => var_now%p%r_ptr(:, :, ncontained, 1, 1)
             buffer_ptr(1, 1)%p(1:var_size) => tmp_buffer
          CASE (4)
             tmp_buffer => var_now%p%r_ptr(:, :, 1, ncontained, 1)
             buffer_ptr(1, 1)%p(1:var_size) => tmp_buffer
          CASE (5)
             tmp_buffer => var_now%p%r_ptr(:, :, 1, 1, ncontained)
             buffer_ptr(1, 1)%p(1:var_size) => tmp_buffer
          CASE DEFAULT
             CALL finish(str_module, "Unsupported var_ref_pos " // int2string(var_ref_pos) // &
                  " for variable " // TRIM(var_now%p%info%name))
          END SELECT
       ELSE
          DO nn = 1 , collection_size
             SELECT CASE (var_ref_pos)
             CASE (1)
                buffer(:,nn) = RESHAPE(var_now%p%r_ptr(ncontained, :, nn, :, 1), &
                     (/ var_size /))
             CASE (2)
                buffer(:,nn) = RESHAPE(var_now%p%r_ptr(:, ncontained, nn, :, 1), &
                     (/ var_size /))
             CASE (3)
                buffer(:,nn) = RESHAPE(var_now%p%r_ptr(:, nn, ncontained, :, 1), &
                     (/ var_size /))
             CASE (4)
                buffer(:,nn) = RESHAPE(var_now%p%r_ptr(:, nn, :, ncontained, 1), &
                     (/ var_size /))
             CASE (5)
                buffer(:,nn) = RESHAPE(var_now%p%r_ptr(:, nn, :, 1, ncontained), &
                     (/ var_size /))
             CASE DEFAULT
                CALL finish(str_module, "Unsupported var_ref_pos " // int2string(var_ref_pos) // &
                     " for variable " // TRIM(var_now%p%info%name))
             END SELECT
             buffer_ptr(1, nn)%p => buffer(:,nn)
          ENDDO
       END IF

       ! The ocean model does not mask land cells hence we set them to NaN manually before coupling to YAC.
       IF ( PRESENT(valid_mask) ) THEN
          DO nn = 1 , collection_size
             ! dont overwrite the variable itself
             IF (.NOT. ASSOCIATED(buffer_ptr(1, nn)%p, buffer(:,nn))) THEN
                buffer(:,nn) = buffer_ptr(1, nn)%p
                buffer_ptr(1, nn)%p => buffer(:,nn)
             ENDIF
             ! Duplicate first level of ocean-mask for half-depth fields.
             WHERE ( RESHAPE(valid_mask(:, MAX(1, nn - MAX(0, collection_size - nlev)), :), &
                  (/ var_size /)) .LT. 0.5_wp )
                buffer_ptr(1, nn)%p = ieee_value(buffer_ptr(1, nn)%p, ieee_quiet_nan)
             ENDWHERE
          ENDDO
       ENDIF

       IF (ltimer) CALL timer_start(timer_coupling_put)
       CALL yac_fput(cur_field%yac_field_id, 1, &
            collection_size, buffer_ptr(:, 1:collection_size), info, ierror)
       IF (ltimer) CALL timer_stop(timer_coupling_put)
       cur_field => cur_field%next
    ENDDO
  END SUBROUTINE output_coupling

  !>
  !! SUBROUTINE destruct_output_coupling -- destructs the fields list
  SUBROUTINE destruct_output_coupling()
    TYPE(t_exposed_var), POINTER :: exposed_var, tmp
    exposed_var => exposed_vars_head
    DO WHILE(ASSOCIATED(exposed_var))
       tmp => exposed_var
       exposed_var => exposed_var%next
       DEALLOCATE(tmp)
    END DO
  END SUBROUTINE destruct_output_coupling


END MODULE mo_output_coupling
