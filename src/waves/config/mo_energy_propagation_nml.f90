! Namelist for wave energy propagation
!
! The following routine is called by read_wave_namelist and controls
! specifics of the wave energy propagation.
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
MODULE mo_energy_propagation_nml

  USE mo_kind,                      ONLY: wp
  USE mo_exception,                 ONLY: finish, message, print_value
  USE mo_io_units,                  ONLY: nnml, nnml_output
  USE mo_impl_constants,            ONLY: max_dom, VNAME_LEN, MAX_NTRACER
  USE mo_master_control,            ONLY: use_restart_namelists
  USE mo_namelist,                  ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_restart_nml_and_att,       ONLY: open_tmpfile, store_and_close_namelist,     &
    &                                     open_and_restore_namelist, close_tmpfile
  USE mo_mpi,                       ONLY: my_process_is_stdio
  USE mo_nml_annotate,              ONLY: temp_defaults, temp_settings
  USE mo_energy_propagation_config, ONLY: energy_propagation_config

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_energy_propagation_nml

CONTAINS

  !>
  !! Read namelist for wave energy propagation
  !!
  !! This subroutine
  !! - reads the Namelist for wave energy propagation
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_energy_propagation_nml (filename)

    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: istat, funit
    INTEGER :: jg          !< patch loop index
    INTEGER :: iunit
    INTEGER :: it

    CHARACTER(len=VNAME_LEN) ::  &    !< tracer-specific name suffixes
      &  tracer_names(MAX_NTRACER)    !< these are only required for
    CHARACTER(len=VNAME_LEN) :: tname !< tracer name

    INTEGER :: itype_limit    !< selects the numerical flux limiter
                              !< 0: no limiter
                              !< 3: monotonous flux limiter
                              !< 4: positive definite flux limiter

    INTEGER :: igrad_c_miura  !< selects the reconstruction method for
                              !< the gradient defined at the cell center.
                              !< Used by the transport scheme.

    REAL(wp):: beta_fct       !< global boost factor for range of permissible values in
                              !< (semi-) monotonous flux limiter. A value larger than
                              !< 1 allows for (small) over and undershoots, while a value
                              !< of 1 gives strict monotonicity (at the price of increased
                              !< diffusivity).

    LOGICAL :: lgrid_refr     !< if .TRUE., calculate grid refraction

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_energy_propagation_nml: read_energy_propagation_nml'

    NAMELIST /energy_propagation_nml/ itype_limit, igrad_c_miura, beta_fct, lgrid_refr

    !-----------------------
    ! 1. default settings
    !-----------------------

    itype_limit   = 0          ! no flux limiter
    igrad_c_miura = 1          ! linear least squares reconstruction
    beta_fct      = 1.005_wp   ! factor of allowed over-/undershooting in flux limiter
    lgrid_refr    = .TRUE.     ! if .TRUE., calculate grid refraction

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, energy_propagation_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('energy_propagation_nml')
      READ(funit,NML=energy_propagation_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('energy_propagation_nml', STATUS=istat)
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, energy_propagation_nml)       ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, energy_propagation_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------

    IF ( ALL((/0,3,4/) /= itype_limit) ) THEN
      CALL finish(routine,                                     &
        &  'incorrect settings for itype_limit. Must be 0, 3, or 4')
    ENDIF

    ! setting default names for tracers
    !
    CALL message(routine,'Setting default names for tracers: q<no>')
    DO it=1, MAX_NTRACER
      WRITE(tname,'(i3)') it
      tracer_names(it) = 'q'//TRIM(ADJUSTL(tname))
      CALL print_value('tracer name '//TRIM(ADJUSTL(tracer_names(it)))//' in position',it)
    END DO
    CALL message(' ',' ')


    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------
    DO jg = 1,max_dom
      energy_propagation_config(jg)%itype_limit     = itype_limit
      energy_propagation_config(jg)%igrad_c_miura   = igrad_c_miura
      energy_propagation_config(jg)%beta_fct        = beta_fct
      energy_propagation_config(jg)%lgrid_refr      = lgrid_refr
      energy_propagation_config(jg)%tracer_names(:) = tracer_names(:)
    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=energy_propagation_nml)
      CALL store_and_close_namelist(funit, 'energy_propagation_nml')
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=energy_propagation_nml)

  END SUBROUTINE read_energy_propagation_nml

END MODULE mo_energy_propagation_nml
