!> @file comin_run_recorder_plugin.F90
!! @brief A plugin that records all information to a netcdf file to replay the run.
! Variables need to be handled separately by `comin_var_recorder_plugin`.
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

MODULE comin_run_recorder_plugin

  USE mpi,                    ONLY: MPI_COMM_SIZE
  USE comin_plugin_interface, ONLY: t_comin_plugin_info, t_comin_setup_version_info,                 &
       &                            t_comin_descrdata_global, t_comin_descrdata_simulation_interval, &
       &                            comin_setup_get_version, comin_current_get_plugin_info,          &
       &                            comin_parallel_get_host_mpi_rank, comin_callback_register,       &
       &                            comin_plugin_finish, comin_descrdata_get_global,                 &
       &                            comin_descrdata_get_domain, comin_descrdata_get_timesteplength,  &
       &                            comin_descrdata_get_simulation_interval, comin_current_get_ep,   &
       &                            comin_current_get_domain_id, comin_current_get_datetime,         &
       &                            EP_DESTRUCTOR, comin_parallel_get_host_mpi_comm,                 &
       &                            comin_error_check
  USE netcdf,                 ONLY: nf90_create, nf90_put_att, nf90_def_dim, nf90_def_var,           &
       &                            nf90_def_grp, nf90_put_var, nf90_close, NF90_CHAR, NF90_GLOBAL,  &
       &                            NF90_INT, NF90_NETCDF4, NF90_UNLIMITED
  USE netcdf_utils,           ONLY: nf90_utils_def_var, nf90_utils_get_shape, nf90
  USE comin_descrdata_save,   ONLY: comin_descrdata_save_global, comin_descrdata_save_domain
  USE utils,                  ONLY: int2string

  IMPLICIT NONE

  INTEGER :: ierr, host_rank, host_comm, host_comm_size
  INTEGER :: ep_counter = 1
  INTEGER :: ncid, group_ncid
  INTEGER :: callback_call_dimid, datetime_len_dimid
  INTEGER :: current_ep_varid, current_domain_id_varid, current_datetime_varid

  TYPE(t_comin_plugin_info)               :: plugin_info
  TYPE(t_comin_setup_version_info)        :: comin_version
  TYPE(t_comin_descrdata_global), POINTER :: global

CONTAINS

  SUBROUTINE comin_main() BIND(C)

    INTEGER :: ep, jg
    TYPE(t_comin_descrdata_simulation_interval), POINTER :: simulation_interval

    comin_version = comin_setup_get_version()
    CALL comin_current_get_plugin_info(plugin_info)
    host_rank = comin_parallel_get_host_mpi_rank()
    host_comm = comin_parallel_get_host_mpi_comm()
    CALL MPI_COMM_SIZE(host_comm, host_comm_size, ierr)

    ! register callbacks
    DO ep=1,EP_DESTRUCTOR-1
      CALL comin_callback_register(ep, comin_run_recorder_plugin_ep)
    END DO
    CALL comin_callback_register(EP_DESTRUCTOR, comin_run_recorder_plugin_destructor)

    CALL nf90(nf90_create(TRIM(plugin_info%options) // int2string(host_rank) // ".nc", &
                          NF90_NETCDF4, ncid))

    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "comin_version", &
         & [comin_version%version_no_major, &
         &  comin_version%version_no_minor, &
         &  comin_version%version_no_patch]))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "host_comm_rank", host_rank))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "host_comm_size", host_comm_size))

    ! create variables for ep data
    CALL nf90(nf90_def_dim(ncid, "callback", NF90_UNLIMITED, callback_call_dimid))
    CALL nf90(nf90_def_dim(ncid, "datetime_len", 32, datetime_len_dimid))

    CALL nf90(nf90_def_var(ncid, "current_ep", NF90_INT, [ callback_call_dimid ], &
                           current_ep_varid))
    CALL nf90(nf90_def_var(ncid, "current_domain_id", NF90_INT, [ callback_call_dimid ], &
                           current_domain_id_varid))
    CALL nf90(nf90_def_var(ncid, "current_datetime", NF90_CHAR, &
                           [ datetime_len_dimid, callback_call_dimid ], current_datetime_varid))

    global => comin_descrdata_get_global()
    CALL comin_descrdata_save_global(ncid, global)

    DO jg=1, global%n_dom
      CALL nf90(nf90_def_grp(ncid, "domain_"//int2string(jg), group_ncid))
      CALL comin_descrdata_save_domain(group_ncid, comin_descrdata_get_domain(jg))
      CALL nf90(nf90_put_att(group_ncid, NF90_GLOBAL, "timestep_length", &
           & comin_descrdata_get_timesteplength(jg)))
    END DO

    simulation_interval => comin_descrdata_get_simulation_interval()
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "exp_start", simulation_interval%exp_start))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "exp_stop", simulation_interval%exp_stop))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "run_start", simulation_interval%run_start))
    CALL nf90(nf90_put_att(ncid, NF90_GLOBAL, "run_stop", simulation_interval%run_stop))

  END SUBROUTINE comin_main

  RECURSIVE SUBROUTINE comin_run_recorder_plugin_ep() BIND(C)
    INTEGER :: current_ep, current_domain_id
    CHARACTER(LEN=:), ALLOCATABLE :: current_datetime
    current_ep = comin_current_get_ep()
    current_domain_id = comin_current_get_domain_id()
    CALL comin_current_get_datetime(current_datetime)

    CALL nf90(nf90_put_var(ncid, current_ep_varid, current_ep, [ ep_counter ]))
    CALL nf90(nf90_put_var(ncid, current_domain_id_varid, current_domain_id, [ ep_counter ]))
    CALL nf90(nf90_put_var(ncid, current_datetime_varid, TRIM(current_datetime), [ 1, ep_counter ]))
    ep_counter = ep_counter +1
  END SUBROUTINE comin_run_recorder_plugin_ep

  SUBROUTINE comin_run_recorder_plugin_destructor() BIND(C)
    CALL comin_run_recorder_plugin_ep
    CALL nf90(nf90_close(ncid))
  END SUBROUTINE comin_run_recorder_plugin_destructor

END MODULE comin_run_recorder_plugin
