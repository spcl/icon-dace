!> @file comin_plugin_interface.F90
!! @brief The module comin_plugin_interface exports all procedures, variables
!!        and constants that are exposed to third party plugins.
!!
!! @defgroup plugin_interface Plugin Interface
!! @{
!!
!! **Entities that are exposed to both, the host interface and the
!!   plugin interface are listed in the group** @ref common.
!!
!! @}
!!
!! - This module does not provide any actual implementation.
!! - The convention is that from the third party plugin side, no other
!!   module than comin_plugin_interface must be used.
!!
!  @authors 01/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_plugin_interface

  USE comin_setup,             ONLY: comin_current_get_plugin_info
  USE comin_callback,          ONLY: comin_callback_register,         &
    &                                comin_callback_get_ep_name
  USE comin_variable_types,    ONLY: t_comin_var_descriptor, t_comin_var_ptr,                         &
    &                                t_comin_var_descr_list_item
  USE comin_variable,          ONLY: comin_var_get, comin_var_request_add,                            &
    &                                comin_var_get_descr_list_head, comin_var_to_3d
  USE comin_metadata,          ONLY: comin_metadata_set => comin_metadata_set_request,                &
    &                                comin_metadata_get, comin_metadata_get_typeid,                   &
    &                                comin_metadata_get_iterator
  USE comin_metadata_types,    ONLY: t_comin_var_metadata_iterator
  USE comin_descrdata_types,   ONLY: t_comin_descrdata_domain, t_comin_descrdata_global,              &
    &                                t_comin_descrdata_simulation_interval
  USE comin_descrdata,         ONLY: comin_descrdata_get_domain, comin_descrdata_get_global,          &
    &                                comin_descrdata_get_simulation_interval, &
    &                                comin_descrdata_get_index, comin_descrdata_get_block,            &
    &                                comin_descrdata_get_cell_indices,                                &
    &                                comin_descrdata_get_cell_npromz, comin_descrdata_get_edge_npromz, &
    &                                comin_descrdata_get_vert_npromz,                                  &
    &                                comin_descrdata_index_lookup_glb2loc_cell,                        &
    &                                comin_current_get_datetime, comin_descrdata_get_timesteplength
  USE comin_setup_utils,       ONLY: t_comin_setup_version_info,     &
    &                                comin_setup_get_version
  USE comin_setup_constants,   ONLY: wp,                             &
    &                                EP_SECONDARY_CONSTRUCTOR,       &
    &                                EP_ATM_INIT_FINALIZE,           &
    &                                EP_ATM_TIMELOOP_BEFORE,         &
    &                                EP_ATM_TIMELOOP_START,          &
    &                                EP_ATM_TIMELOOP_END,            &
    &                                EP_ATM_TIMELOOP_AFTER,          &
    &                                EP_ATM_INTEGRATE_BEFORE,        &
    &                                EP_ATM_INTEGRATE_START,         &
    &                                EP_ATM_INTEGRATE_END,           &
    &                                EP_ATM_INTEGRATE_AFTER,         &
    &                                EP_ATM_WRITE_OUTPUT_BEFORE,     &
    &                                EP_ATM_WRITE_OUTPUT_AFTER,      &
    &                                EP_ATM_CHECKPOINT_BEFORE,       &
    &                                EP_ATM_CHECKPOINT_AFTER,        &
    &                                EP_ATM_ADVECTION_BEFORE,        &
    &                                EP_ATM_ADVECTION_AFTER,         &
    &                                EP_ATM_PHYSICS_BEFORE,          &
    &                                EP_ATM_PHYSICS_AFTER,           &
    &                                EP_ATM_NUDGING_BEFORE,          &
    &                                EP_ATM_NUDGING_AFTER,           &
    &                                EP_ATM_SURFACE_BEFORE,          &
    &                                EP_ATM_SURFACE_AFTER,           &
    &                                EP_ATM_TURBULENCE_BEFORE,       &
    &                                EP_ATM_TURBULENCE_AFTER,        &
    &                                EP_ATM_MICROPHYSICS_BEFORE,     &
    &                                EP_ATM_MICROPHYSICS_AFTER,      &
    &                                EP_ATM_CONVECTION_BEFORE,       &
    &                                EP_ATM_CONVECTION_AFTER,        &
    &                                EP_ATM_RADIATION_BEFORE,        &
    &                                EP_ATM_RADIATION_AFTER,         &
    &                                EP_ATM_RADHEAT_BEFORE,          &
    &                                EP_ATM_RADHEAT_AFTER,           &
    &                                EP_ATM_GWDRAG_BEFORE,           &
    &                                EP_ATM_GWDRAG_AFTER,            &
    &                                EP_FINISH,                      &
    &                                EP_DESTRUCTOR,                  &
    &                                COMIN_FLAG_READ,                &
    &                                COMIN_FLAG_WRITE,               &
    &                                COMIN_FLAG_DEVICE,              &
    &                                COMIN_ZAXIS_NONE,               &
    &                                COMIN_ZAXIS_2D, COMIN_ZAXIS_3D, &
    &                                COMIN_ZAXIS_3D_HALF,            &
    &                                COMIN_ZAXIS_UNDEF,              &
    &                                COMIN_DOMAIN_OUTSIDE_LOOP,      &
    &                                COMIN_HGRID_UNSTRUCTURED_CELL,  &
    &                                COMIN_HGRID_UNSTRUCTURED_EDGE,  &
    &                                COMIN_HGRID_UNSTRUCTURED_VERTEX
  USE comin_state,             ONLY: comin_current_get_ep, comin_current_get_domain_id, &
    &                                comin_setup_get_verbosity_level
  USE comin_plugin_types,      ONLY: t_comin_plugin_info
  USE comin_parallel,          ONLY: comin_parallel_get_plugin_mpi_comm,                     &
    &                                comin_parallel_get_host_mpi_comm,                       &
    &                                comin_parallel_get_host_mpi_rank
  USE comin_errhandler,        ONLY: comin_plugin_finish, comin_error_get_message,           &
    &                                comin_error_check, comin_error_set_errors_return
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS

  IMPLICIT NONE

  PRIVATE

#include "comin_global.inc"

  ! From comin_callback:
  PUBLIC :: comin_callback_register, comin_current_get_ep

  ! From comin_setup:
  PUBLIC :: comin_current_get_plugin_info,   &
    &       COMIN_ZAXIS_NONE, COMIN_ZAXIS_2D, COMIN_ZAXIS_3D, COMIN_ZAXIS_3D_HALF, COMIN_ZAXIS_UNDEF

  ! From comin_variable:
  PUBLIC :: comin_var_get, t_comin_var_descriptor, t_comin_var_ptr
  PUBLIC :: comin_var_to_3d
  PUBLIC :: comin_var_request_add, comin_var_get_descr_list_head
  PUBLIC :: t_comin_var_descr_list_item
  ! From comin_metadata:
  PUBLIC :: comin_metadata_set
  PUBLIC :: comin_metadata_get
  PUBLIC :: comin_metadata_get_typeid
  PUBLIC :: comin_metadata_get_iterator
  ! From comin_metadata_types:
  PUBLIC :: t_comin_var_metadata_iterator

  ! From comin_descrdata:
  PUBLIC :: comin_descrdata_get_domain, t_comin_descrdata_domain, &
            comin_descrdata_get_global, t_comin_descrdata_global
  PUBLIC :: comin_descrdata_get_simulation_interval, t_comin_descrdata_simulation_interval
  PUBLIC :: comin_descrdata_get_index, comin_descrdata_get_block, comin_descrdata_get_cell_indices
  PUBLIC :: comin_descrdata_get_cell_npromz, comin_descrdata_get_edge_npromz, comin_descrdata_get_vert_npromz
  PUBLIC :: comin_descrdata_index_lookup_glb2loc_cell
  PUBLIC :: comin_current_get_datetime
  PUBLIC :: comin_descrdata_get_timesteplength

  ! From comin_setup_utils:
  PUBLIC :: t_comin_setup_version_info, comin_setup_get_version

  ! From comin_setup_constants:
  PUBLIC :: wp
  PUBLIC :: EP_SECONDARY_CONSTRUCTOR,     &
    &       EP_ATM_INIT_FINALIZE,         &
    &       EP_ATM_TIMELOOP_BEFORE,       &
    &       EP_ATM_TIMELOOP_START,        &
    &       EP_ATM_TIMELOOP_END,          &
    &       EP_ATM_TIMELOOP_AFTER,        &
    &       EP_ATM_INTEGRATE_BEFORE,      &
    &       EP_ATM_INTEGRATE_START,       &
    &       EP_ATM_INTEGRATE_END,         &
    &       EP_ATM_INTEGRATE_AFTER,       &
    &       EP_ATM_WRITE_OUTPUT_BEFORE,   &
    &       EP_ATM_WRITE_OUTPUT_AFTER,    &
    &       EP_ATM_CHECKPOINT_BEFORE,     &
    &       EP_ATM_CHECKPOINT_AFTER,      &
    &       EP_ATM_ADVECTION_BEFORE,      &
    &       EP_ATM_ADVECTION_AFTER,       &
    &       EP_ATM_PHYSICS_BEFORE,        &
    &       EP_ATM_PHYSICS_AFTER,         &
    &       EP_ATM_NUDGING_BEFORE,        &
    &       EP_ATM_NUDGING_AFTER,         &
    &       EP_ATM_SURFACE_BEFORE,        &
    &       EP_ATM_SURFACE_AFTER,         &
    &       EP_ATM_TURBULENCE_BEFORE,     &
    &       EP_ATM_TURBULENCE_AFTER,      &
    &       EP_ATM_MICROPHYSICS_BEFORE,   &
    &       EP_ATM_MICROPHYSICS_AFTER,    &
    &       EP_ATM_CONVECTION_BEFORE,     &
    &       EP_ATM_CONVECTION_AFTER,      &
    &       EP_ATM_RADIATION_BEFORE,      &
    &       EP_ATM_RADIATION_AFTER,       &
    &       EP_ATM_RADHEAT_BEFORE,        &
    &       EP_ATM_RADHEAT_AFTER,         &
    &       EP_ATM_GWDRAG_BEFORE,         &
    &       EP_ATM_GWDRAG_AFTER,          &
    &       EP_FINISH,                    &
    &       EP_DESTRUCTOR,                &
    &       COMIN_HGRID_UNSTRUCTURED_CELL,  &
    &       COMIN_HGRID_UNSTRUCTURED_EDGE,  &
    &       COMIN_HGRID_UNSTRUCTURED_VERTEX
  PUBLIC :: COMIN_FLAG_READ, COMIN_FLAG_WRITE, COMIN_FLAG_DEVICE
  PUBLIC :: COMIN_DOMAIN_OUTSIDE_LOOP
  PUBLIC :: comin_callback_get_ep_name

  ! From comin_plugin_types
  PUBLIC :: t_comin_plugin_info

  ! From comin_parallel:
  PUBLIC :: comin_parallel_get_plugin_mpi_comm, comin_parallel_get_host_mpi_comm
  PUBLIC :: comin_parallel_get_host_mpi_rank

  ! From comin_errhandler:
  PUBLIC :: comin_plugin_finish, comin_error_get_message, comin_error_check

  ! From comin_errhandler_constants
  PUBLIC :: COMIN_SUCCESS

  ! From comin_state:
  PUBLIC :: comin_current_get_domain_id
  PUBLIC :: comin_setup_get_verbosity_level

  ! From comin_global.inc
  PUBLIC :: COMIN_METADATA_TYPEID_UNDEFINED
  PUBLIC :: COMIN_METADATA_TYPEID_INTEGER
  PUBLIC :: COMIN_METADATA_TYPEID_REAL
  PUBLIC :: COMIN_METADATA_TYPEID_CHARACTER
  PUBLIC :: COMIN_METADATA_TYPEID_LOGICAL

  ! Create fortran integer parameters from values defined in comin_global.inc
  INTEGER, PARAMETER :: COMIN_METADATA_TYPEID_UNDEFINED = COMIN_TYPEID_UNDEFINED
  INTEGER, PARAMETER :: COMIN_METADATA_TYPEID_INTEGER   = COMIN_TYPEID_INTEGER
  INTEGER, PARAMETER :: COMIN_METADATA_TYPEID_REAL      = COMIN_TYPEID_REAL
  INTEGER, PARAMETER :: COMIN_METADATA_TYPEID_CHARACTER = COMIN_TYPEID_CHARACTER
  INTEGER, PARAMETER :: COMIN_METADATA_TYPEID_LOGICAL   = COMIN_TYPEID_LOGICAL

END MODULE comin_plugin_interface
