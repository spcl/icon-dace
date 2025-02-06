!> comin_host_interface.F90
!! @brief ComIn entities exposed to the host model (e.g. ICON).
!!
!!
!! @defgroup host_interface Host Interface
!! @{
!!
!! **Entities that are exposed to both, the host interface and the
!!   plugin interface are listed in the group** @ref common.
!!
!! @}
!
!  - No actual implementations within this module
!
!  @authors 01/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_host_interface

  USE comin_setup,             ONLY: comin_setup_check,              &
    &                                comin_plugin_primaryconstructor,&
    &                                comin_setup_finalize,           &
    &                                comin_descrdata_set_fct_glb2loc_cell, &
    &                                comin_setup_init, &
    &                                comin_setup_errhandler
  USE comin_callback,          ONLY: comin_callback_context_call,    &
    &                                comin_callback_get_ep_name

  USE comin_setup_constants,   ONLY: EP_SECONDARY_CONSTRUCTOR,       &
    &                                EP_ATM_YAC_DEFCOMP_BEFORE,      &
    &                                EP_ATM_YAC_DEFCOMP_AFTER,       &
    &                                EP_ATM_YAC_SYNCDEF_BEFORE,      &
    &                                EP_ATM_YAC_SYNCDEF_AFTER,       &
    &                                EP_ATM_YAC_ENDDEF_BEFORE,       &
    &                                EP_ATM_YAC_ENDDEF_AFTER,        &
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
    &                                COMIN_ZAXIS_NONE,               &
    &                                COMIN_ZAXIS_2D, COMIN_ZAXIS_3D, &
    &                                COMIN_ZAXIS_3D_HALF,            &
    &                                COMIN_ZAXIS_UNDEF,              &
    &                                COMIN_DOMAIN_OUTSIDE_LOOP,      &
    &                                COMIN_HGRID_UNSTRUCTURED_CELL,  &
    &                                COMIN_HGRID_UNSTRUCTURED_EDGE,  &
    &                                COMIN_HGRID_UNSTRUCTURED_VERTEX
  USE comin_setup_utils,       ONLY: t_comin_setup_version_info, &
    &                                comin_setup_get_version
  USE comin_variable_types,    ONLY: t_comin_var_ptr, t_comin_var_descriptor,                 &
    &                                t_comin_var_descr_list_item, t_var_request_list_item
  USE comin_variable,          ONLY: comin_var_list_finalize, comin_var_get_descr_list_head,  &
    &                                comin_var_list_append,  comin_request_get_list_head,     &
    &                                comin_var_update, comin_var_set_sync_device_mem,         &
    &                                comin_var_request_list_finalize
  USE comin_metadata,          ONLY: comin_metadata_set => comin_metadata_set_host, &
    &                                comin_metadata_get_typeid, comin_metadata_get_iterator, &
    &                                comin_metadata_get_or
  USE comin_metadata_types,    ONLY: t_comin_var_metadata_iterator
  USE comin_parallel,          ONLY: comin_parallel_get_host_mpi_comm,   &
    &                                comin_parallel_mpi_handshake
  USE comin_plugin_types,      ONLY: t_comin_plugin_description
  USE comin_descrdata_types,   ONLY: t_comin_descrdata_global, t_comin_descrdata_domain,  &
    &                                t_comin_descrdata_simulation_interval
  USE comin_descrdata,         ONLY: comin_descrdata_finalize,                                     &
    &                                comin_descrdata_set_global, comin_descrdata_set_domain,       &
    &                                comin_descrdata_set_simulation_interval, &
    &                                comin_current_set_datetime, comin_descrdata_set_timesteplength
  USE mo_mpi_handshake,        ONLY: mpi_handshake, mpi_handshake_dummy
  USE comin_state,             ONLY: comin_setup_set_verbosity_level, &
    &                                comin_setup_get_verbosity_level
  USE comin_errhandler_constants, ONLY: COMIN_SUCCESS

  IMPLICIT NONE

  PRIVATE

#include "comin_global.inc"

  ! From comin_state:
  PUBLIC :: comin_setup_set_verbosity_level
  PUBLIC :: comin_setup_get_verbosity_level
  PUBLIC :: comin_descrdata_set_fct_glb2loc_cell
  PUBLIC :: comin_setup_errhandler

  ! From comin_setup:
  PUBLIC :: comin_setup_init
  PUBLIC :: comin_setup_check, comin_plugin_primaryconstructor, comin_setup_finalize

  ! From comin_plugin_types:
  PUBLIC :: t_comin_plugin_description

  ! From comin_callback:
  PUBLIC :: comin_callback_context_call

  ! From comin_setup_utils:
  PUBLIC :: t_comin_setup_version_info,   &
    &       comin_setup_get_version
  PUBLIC :: EP_SECONDARY_CONSTRUCTOR,     &
    &       EP_ATM_YAC_DEFCOMP_BEFORE,    &
    &       EP_ATM_YAC_DEFCOMP_AFTER,     &
    &       EP_ATM_YAC_SYNCDEF_BEFORE,    &
    &       EP_ATM_YAC_SYNCDEF_AFTER,     &
    &       EP_ATM_YAC_ENDDEF_BEFORE,     &
    &       EP_ATM_YAC_ENDDEF_AFTER,      &
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
    &       EP_DESTRUCTOR
  PUBLIC :: COMIN_ZAXIS_NONE, COMIN_ZAXIS_2D, COMIN_ZAXIS_3D, COMIN_ZAXIS_3D_HALF, COMIN_ZAXIS_UNDEF
  PUBLIC :: COMIN_DOMAIN_OUTSIDE_LOOP
  PUBLIC :: comin_callback_get_ep_name
  PUBLIC :: COMIN_HGRID_UNSTRUCTURED_CELL, COMIN_HGRID_UNSTRUCTURED_EDGE,  &
    &                                COMIN_HGRID_UNSTRUCTURED_VERTEX

  ! From comin_variable:
  PUBLIC :: comin_var_list_finalize, comin_var_get_descr_list_head
  PUBLIC :: comin_var_request_list_finalize
  PUBLIC :: t_comin_var_ptr, t_comin_var_descriptor,                       &
    &       comin_var_list_append
  PUBLIC :: comin_request_get_list_head, comin_var_update
  PUBLIC :: t_comin_var_descr_list_item, t_var_request_list_item
  PUBLIC :: comin_var_set_sync_device_mem
  ! From comin_metadata:
  PUBLIC :: comin_metadata_set
  PUBLIC :: comin_metadata_get_typeid
  PUBLIC :: comin_metadata_get_iterator
  PUBLIC :: comin_metadata_get_or
  ! From comin_metadata_types:
  PUBLIC :: t_comin_var_metadata_iterator
  ! From comin_parallel:
  PUBLIC :: comin_parallel_get_host_mpi_comm
  PUBLIC :: comin_parallel_mpi_handshake

  ! From comin_descrdata:
  PUBLIC :: comin_descrdata_finalize
  PUBLIC :: t_comin_descrdata_global, comin_descrdata_set_global, &
            t_comin_descrdata_domain, comin_descrdata_set_domain
  PUBLIC :: t_comin_descrdata_simulation_interval, comin_descrdata_set_simulation_interval, &
    &       comin_current_set_datetime, comin_descrdata_set_timesteplength

  ! From mo_mpi_handshake
  PUBLIC :: mpi_handshake, mpi_handshake_dummy

  ! From comin_errhandler_constants
  PUBLIC :: COMIN_SUCCESS

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

END MODULE comin_host_interface
