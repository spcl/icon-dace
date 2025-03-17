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

MODULE mo_output_coupling

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_var                 ,ONLY: t_var_ptr
  USE mo_var_groups          ,ONLY: MAX_GROUPS, var_groups_dyn
  USE mo_run_config          ,ONLY: nlev, msg_level
  USE mo_run_config          ,ONLY: ltimer
  USE mo_timer               ,ONLY: timer_start, timer_stop, &
       &                            timer_coupling_output_put, timer_coupling_output_1stput, &
       &                            timer_coupling_output, timer_coupling_output_buf_prep
  USE mo_util_string         ,ONLY: int2string
  USE mo_exception           ,ONLY: message, finish
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_zaxis_type          ,ONLY: zaxisTypeList






  CHARACTER(len=*), PARAMETER :: str_module = 'mo_output_coupling' ! Output of module for debug
  CHARACTER(len=2), PARAMETER :: newline = ACHAR(13) // ACHAR(10)

  TYPE, PRIVATE :: t_exposed_var
     INTEGER :: yac_field_id
     TYPE(t_var_ptr) :: var(1:3)
     INTEGER :: tlev_source, var_size
     TYPE(t_exposed_var), POINTER :: next => NULL()
  END TYPE t_exposed_var

  PUBLIC :: construct_output_coupling
  PUBLIC :: construct_output_coupling_finalize
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
    p_patch, comp_id, cell_point_id, vertex_point_id, timestepstring)

    USE mo_var_list_register,   ONLY: t_vl_register_iter
    USE mo_var_metadata,        ONLY: get_var_timelevel, get_var_name
    USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT
    USE mo_var,                 ONLY: level_type_ml
    USE mo_coupling_utils,      ONLY: cpl_get_instance_id





    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id, vertex_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    TYPE(t_vl_register_iter), ALLOCATABLE :: vl_iter
    TYPE(t_exposed_var), POINTER :: exposed_var
    CHARACTER(len=:), ALLOCATABLE :: var_name, metadata, comp_name, grid_name
    INTEGER :: iv, tl, collection_size, key_notl, count = 0, var_size, grpi, nblks, pos(3)
    INTEGER :: point_id, var_ref_pos, instance_id

    TYPE t_tmp_timelevel_var
       INTEGER :: key_notl
       TYPE(t_exposed_var), POINTER :: exposed_var
       TYPE(t_tmp_timelevel_var), POINTER :: next => NULL()
    END type t_tmp_timelevel_var

    TYPE(t_tmp_timelevel_var), POINTER :: tmp_timelevel_var_head => NULL(), tmp_timelevel_var => NULL()


    CALL finish(str_module // 'construct_output_coupling', &
                'built without coupling support.')
  END SUBROUTINE construct_output_coupling


  !>
  !! SUBROUTINE construct_output_coupling_finalize -- sort out all non-coupled fields
  !! from the field_list (has to be called after the enddef operation)

  SUBROUTINE construct_output_coupling_finalize()

   CALL finish(str_module // 'construct_output_coupling_finalize', &
               "built without coupling support.")
  END SUBROUTINE construct_output_coupling_finalize

  !>
  !! SUBROUTINE output_coupling -- Exchange fields between
  !! atmosphere and output components.
  SUBROUTINE output_coupling (valid_mask)

    USE, INTRINSIC :: ieee_arithmetic
    USE mo_impl_constants      ,ONLY: TLEV_NNOW, TLEV_NNEW, TLEV_NNOW_RCF, TLEV_NNEW_RCF
    USE mo_dynamics_config,     ONLY: nnow, nnow_rcf, nnew, nnew_rcf

    REAL(wp), OPTIONAL :: valid_mask(:,:,:)

   CALL finish(str_module // 'output_coupling', &
               'built without coupling support')
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


  SUBROUTINE construct_output_nml_coupling(comp_id, cell_point_id, vertex_point_id)

    USE mo_cdi_constants,          ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT
    USE mo_exception,              ONLY: message_text
    USE mo_name_list_output_init,  ONLY: output_file, nlevs_of_var
    USE mo_name_list_output_types, ONLY: FILETYPE_YAC, t_output_name_list
    USE mo_var_metadata,           ONLY: get_var_name
    USE mo_var_metadata_types,     ONLY: t_var_metadata

    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: cell_point_id, vertex_point_id

    CALL finish(str_module // 'construct_output_coupling', &
                'built without coupling support.')
  END SUBROUTINE construct_output_nml_coupling

END MODULE mo_output_coupling
