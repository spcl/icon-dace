! Contains the main stepping routine the 3-dim hydrostatic ocean model.
!
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

MODULE mo_ocean_output
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp
  USE mo_impl_constants,         ONLY: max_char_length
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d, t_subset_range, t_patch_vert
  USE mo_grid_config,            ONLY: n_dom
  USE mo_sync,                   ONLY: sync_patch_array, sync_e, sync_c !, sync_v
  USE mo_ocean_nml,              ONLY: iswm_oce, n_zlev, no_tracer, &
    & diagnostics_level, &
    & eos_type, i_sea_ice, gibraltar
  USE mo_dynamics_config,        ONLY: nold, nnew
  USE mo_io_config,              ONLY: timeSteps_per_outputStep
  USE mo_run_config,             ONLY: nsteps, dtime, ltimer, output_mode
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_ocean_types, ONLY: t_hydro_ocean_state, t_hydro_ocean_diag, &
    & t_hydro_ocean_prog, t_operator_coeff
  USE mo_ocean_state, ONLY: ocean_restart_list
  USE mo_sea_ice_types,          ONLY: t_atmos_fluxes, t_sea_ice
  USE mo_ocean_surface_types,    ONLY: t_ocean_surface
  !USE mo_ocean_physics,            ONLY: t_ho_params
  USE mo_name_list_output,       ONLY: write_name_list_output, istime4name_list_output
  USE mo_ocean_diagnostics, ONLY: destruct_oce_diagnostics, calc_moc, calc_psi
  USE mo_mpi,                    ONLY: my_process_is_stdio
  USE mo_statistics
  USE mo_sea_ice_nml,            ONLY: i_ice_dyn
  USE mo_util_dbg_prnt,          ONLY: dbg_print
  USE mtime,                     ONLY: datetime, MAX_DATETIME_STR_LEN, datetimeToPosixString
  USE mo_mpi,                    ONLY: i_am_accel_node
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: output_ocean

  INTEGER :: nsteps_since_last_output = 0

  !-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !
!<Optimize:inUse>
  SUBROUTINE output_ocean( &
    & patch_3d,        &
    & ocean_state,     &
    & this_datetime,   &
    & surface_fluxes,  &
    & sea_ice,         &
    & jstep, jstep0,   &
    & force_output)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)          :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(inout) :: ocean_state(n_dom)
    TYPE(datetime), POINTER                          :: this_datetime
    TYPE(t_ocean_surface)                            :: surface_fluxes
    TYPE (t_sea_ice),         INTENT(inout)          :: sea_ice
    INTEGER,   INTENT(in)                            :: jstep, jstep0
    LOGICAL, OPTIONAL                                :: force_output
   
    ! local variables
    LOGICAL :: use_force_output
    INTEGER :: jg, jtrc, out_step
    INTEGER :: ocean_statistics
    !LOGICAL                         :: l_outputtime
    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_patch_vert), POINTER :: patch_1d
    INTEGER, POINTER :: dolic(:,:)
    REAL(wp), POINTER :: prism_thickness(:,:,:)
    
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: datestring
    CHARACTER(len=32) :: fmtstr
    
    !CHARACTER(LEN=filename_max)  :: outputfile, gridfile
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_output:output_ocean'
    !------------------------------------------------------------------

    patch_2D      => patch_3d%p_patch_2d(1)
    jg = 1
    nsteps_since_last_output = nsteps_since_last_output + 1
    !------------------------------------------------------------------
    use_force_output = .false.
    IF (PRESENT(force_output)) &
      use_force_output = force_output

    out_step = jstep
    IF (use_force_output) THEN
      out_step = jstep - nsteps_since_last_output + timeSteps_per_outputStep
    ENDIF

   !write(0,*) "out_step=", jstep, nsteps_since_last_output, timeSteps_per_outputStep, out_step
   IF (.not. istime4name_list_output(out_step) )  RETURN
   !write(0,*) "write ....."

    !------------------------------------------------------------------
    fmtstr = '%Y-%m-%d %H:%M:%S'
    call datetimeToPosixString(this_datetime, datestring, fmtstr)
    WRITE(message_text,'(a,a)') 'Write output at:', TRIM(datestring)
    CALL message (TRIM(routine),message_text)




    IF (output_mode%l_nml) CALL write_name_list_output(out_step)


  END SUBROUTINE output_ocean
  !-------------------------------------------------------------------------

END MODULE mo_ocean_output
