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

MODULE mo_ser_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist

  IMPLICIT NONE
  PUBLIC

  INTEGER, TARGET ::  ser_initialization(3)
  INTEGER, TARGET ::  ser_output_diag_dyn(3)
  INTEGER, TARGET ::  ser_output_diag(3)
  INTEGER, TARGET ::  ser_output_opt(3)
  INTEGER, TARGET ::  ser_latbc_data(3)
  INTEGER, TARGET ::  ser_nesting_save_progvars(3)
  INTEGER, TARGET ::  ser_dynamics(3)
  INTEGER, TARGET ::  ser_diffusion(3)
  INTEGER, TARGET ::  ser_nesting_compute_tendencies(3)
  INTEGER, TARGET ::  ser_nesting_boundary_interpolation(3)
  INTEGER, TARGET ::  ser_nesting_relax_feedback(3)
  INTEGER, TARGET ::  ser_step_advection(3)
  INTEGER, TARGET ::  ser_physics(3)
  INTEGER, TARGET ::  ser_physics_init(3)
  INTEGER, TARGET ::  ser_lhn(3)
  INTEGER, TARGET ::  ser_nudging(3)
  INTEGER, TARGET ::  ser_surface(3)
  INTEGER, TARGET ::  ser_microphysics(3)
  INTEGER, TARGET ::  ser_turbdiff(3)
  INTEGER, TARGET ::  ser_turbtrans(3)
  INTEGER, TARGET ::  ser_convection(3)
  INTEGER, TARGET ::  ser_cover(3)
  INTEGER, TARGET ::  ser_set_reff(3)
  INTEGER, TARGET ::  ser_radiation(3)
  INTEGER, TARGET ::  ser_radheat(3)
  INTEGER, TARGET ::  ser_gwdrag(3)
  INTEGER, TARGET ::  ser_time_loop_end(3)
  INTEGER, TARGET ::  ser_reset_to_initial_state(3)
  INTEGER, TARGET ::  ser_all_debug(3)            !! serialize statements using ser_all anywhere for debug purposes

  LOGICAL  ::  ser_rbf
  REAL(wp) ::  ser_nfail
  INTEGER  ::  ser_nreport
  LOGICAL  ::  ser_debug                          !! serialize the debug calls from mo_ser_debug

  NAMELIST /ser_nml/ ser_initialization, ser_output_diag, ser_output_opt, ser_output_diag_dyn, ser_latbc_data, &
  &                  ser_nesting_save_progvars, ser_dynamics, &
  &                  ser_debug, ser_rbf, ser_diffusion, ser_nesting_compute_tendencies, &
  &                  ser_nesting_boundary_interpolation, ser_nesting_relax_feedback, ser_physics_init, &
  &                  ser_step_advection, ser_physics, ser_lhn, ser_nudging, ser_all_debug, ser_surface, &
  &                  ser_microphysics, ser_convection, ser_cover, ser_set_reff, ser_radiation, ser_radheat, &
  &                  ser_turbtrans, ser_turbdiff, ser_gwdrag, ser_time_loop_end, ser_reset_to_initial_state, &
  &                  ser_nfail, ser_nreport

  CONTAINS

  SUBROUTINE read_ser_namelist( filename )

   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER :: istat, funit 
   INTEGER :: iunit
   INTEGER :: param_def(3)

   ! turn serialization off by default
   !             # of times to serialize, relative threshold, absolute threshold
   param_def = (/0,                       12,                 12/)
   ser_initialization = param_def
   ser_output_diag_dyn = param_def
   ser_output_diag = param_def
   ser_output_opt = param_def
   ser_latbc_data = param_def
   ser_nesting_save_progvars = param_def
   ser_dynamics = param_def
   ser_diffusion = param_def
   ser_nesting_compute_tendencies = param_def
   ser_nesting_boundary_interpolation = param_def
   ser_nesting_relax_feedback = param_def
   ser_step_advection = param_def
   ser_physics = param_def
   ser_physics_init = param_def
   ser_lhn = param_def
   ser_nudging = param_def
   ser_surface = param_def
   ser_microphysics = param_def
   ser_turbtrans = param_def
   ser_turbdiff = param_def
   ser_convection = param_def
   ser_cover = param_def
   ser_set_reff = param_def
   ser_radiation = param_def
   ser_radheat = param_def
   ser_gwdrag = param_def
   ser_time_loop_end = param_def
   ser_reset_to_initial_state = param_def
   ser_all_debug = param_def
   ser_nfail = 1._wp
   ser_nreport = 10
   ser_debug = .FALSE.
   ser_rbf = .FALSE.

   IF (my_process_is_stdio()) THEN
     iunit = temp_defaults()
     WRITE(iunit, ser_nml)   ! write defaults to temporary text file
   END IF

   CALL open_nml(TRIM(filename))
   CALL position_nml ('ser_nml', status=istat)

   SELECT CASE (istat)
   CASE (POSITIONED)
      READ (nnml, ser_nml)   ! overwrite default settings

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, ser_nml)   ! write settings to temporary text file
      END IF
   END SELECT

   IF(my_process_is_stdio())  THEN
     funit = open_tmpfile()
     WRITE(funit,NML=ser_nml)                    
     CALL store_and_close_namelist(funit, 'ser_nml') 
   ENDIF

   IF(my_process_is_stdio()) WRITE(nnml_output,nml=ser_nml)

  END SUBROUTINE read_ser_namelist

END MODULE mo_ser_nml
