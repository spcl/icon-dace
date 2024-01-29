! Contains the variables to set up the grid configuration
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

MODULE mo_grid_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,               ONLY: wp
! USE mo_exception,          ONLY: message, message_text, finish
  USE mo_io_units,           ONLY: nnml, nnml_output,filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_math_constants,     ONLY: rad2deg
!  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_physical_constants, ONLY: earth_angular_velocity
!  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
!    &                               open_and_restore_namelist, close_tmpfile

  USE mo_grid_config,        ONLY:                                          &
    & config_lfeedback                    => lfeedback,                     &
    & config_ifeedback_type               => ifeedback_type,                &
    & config_start_time                   => start_time,                    &
    & config_end_time                     => end_time,                      &
    & config_lplane                       => lplane,                        &
    & config_is_plane_torus               => is_plane_torus,                &
    & config_corio_lat                    => corio_lat,                     &
    & config_l_limited_area               => l_limited_area,                &
    & config_l_scm_mode                   => l_scm_mode,                    &
    & config_patch_weight                 => patch_weight,                  &
    & config_lredgrid_phys                => lredgrid_phys,                 &
    & config_nexlevs_rrg_vnest            => nexlevs_rrg_vnest,             &
    & config_dynamics_grid_filename       => dynamics_grid_filename,        &
    & config_dynamics_parent_grid_id      => dynamics_parent_grid_id,       &
    & config_radiation_grid_filename      => radiation_grid_filename,       &
    & config_grid_rescale_factor          => grid_rescale_factor,           &
    & config_grid_angular_velocity        => namelist_grid_angular_velocity,&
    & config_lrescale_timestep            => lrescale_timestep,             &
    & config_lrescale_ang_vel             => lrescale_ang_vel,              &
    & config_use_duplicated_connectivity  => use_duplicated_connectivity,   &
    & config_use_dummy_cell_closure       => use_dummy_cell_closure,        &
    & DEFAULT_ENDTIME,                                         &
    & config_create_vgrid                 => create_vgrid,                  &
    & config_vertical_grid_filename       => vertical_grid_filename,        &
    & config_vct_filename                 => vct_filename

  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_grid_namelist
  

  ! ------------------------------------------------------------------------
  ! 1.0 Namelist variables and auxiliary variables
  ! ------------------------------------------------------------------------


  CONTAINS

 !!  Initialization of grid namelist variables
 !!
  SUBROUTINE read_grid_namelist( filename )
    
    CHARACTER(LEN=*), INTENT(IN) :: filename                                           
    INTEGER  :: i_status, i
!    INTEGER  :: funit

    LOGICAL    :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
    INTEGER    :: ifeedback_type           ! type of feedback (incremental or relaxation)
    REAL(wp)   :: start_time(max_dom)      ! time at which execution of a (nested) model domain starts
    REAL(wp)   :: end_time(max_dom)        ! time at which execution of a (nested) model domain terminates
    LOGICAL    :: lredgrid_phys(max_dom)   ! If set to .true. radiation is calculated on a reduced grid
    INTEGER    :: nexlevs_rrg_vnest        ! Maximum number of extra model layers used for calculating radiation if 
                                           ! a reduced radiation grid is combined with vertical nesting
    LOGICAL    :: l_limited_area           ! limited area setup where forcing comes in from sides            
    LOGICAL    :: l_scm_mode               ! SCM mode is designed for tests where all columns are identical
    LOGICAL    :: lsep_grfinfo             ! If .true., read fields related to grid refinement from separate 
                                           ! grid files
    LOGICAL    :: use_duplicated_connectivity ! if true, the zero connectivity is replaced by the last non-zero value
    LOGICAL    :: use_dummy_cell_closure   ! if true then create a dummy cell and connect it to cells and edges with no neigbor

    LOGICAL    :: lplane                   ! f-plane option
    LOGICAL    :: is_plane_torus           ! f-plane with doubly periodic boundary==> like a plane torus
    REAL(wp)   :: corio_lat                ! Latitude, where the f-plane is located if lplane=.true.
  
    REAL(wp)   :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                           ! for any of the first level child patches,
                                           ! processor splitting will be performed

    CHARACTER(LEN=filename_max) :: dynamics_grid_filename(max_dom)
    INTEGER                     :: dynamics_parent_grid_id(max_dom)
    CHARACTER(LEN=filename_max) :: radiation_grid_filename
    INTEGER                     :: dynamics_radiation_grid_link(max_dom)
        
    REAL(wp) :: grid_rescale_factor, grid_angular_velocity
    LOGICAL  :: lrescale_timestep, lrescale_ang_vel
    INTEGER                    :: iunit

    LOGICAL    :: create_vgrid   ! switch if files containing vct_a, vct_b, z_ifc shall be created

    !> files containing vct_a, vct_b, z_ifc
    CHARACTER(LEN=filename_max) :: vertical_grid_filename(max_dom)
    CHARACTER(LEN=filename_max) :: vct_filename ! for reading in vct_a, vct_b
                                                ! we implicitly assume that all other domains 
                                                ! share their levels wih dom01.



    NAMELIST /grid_nml/ lfeedback, ifeedback_type,                 &
      &  lplane, is_plane_torus, corio_lat, l_limited_area,        &
      &  l_scm_mode, grid_rescale_factor, lsep_grfinfo,            &
      &  lrescale_timestep, lrescale_ang_vel,                      &
      &  patch_weight, lredgrid_phys, start_time, end_time,        &
      &  dynamics_grid_filename,  dynamics_parent_grid_id,         &
      &  radiation_grid_filename, nexlevs_rrg_vnest,               &
      &  grid_angular_velocity, use_duplicated_connectivity,       &
      &  use_dummy_cell_closure, create_vgrid,                     &
      &  vertical_grid_filename, vct_filename


!    INTEGER  :: funit
!    INTEGER  :: jg, jlev
!    CHARACTER(filename_max) :: patch_file, gridtype
!    INTEGER  ::  patch_level(max_dom)
!    LOGICAL :: l_exist

    !------------------------------------------------------------
    !  set up the default values for grid_nml
    !------------------------------------------------------------
    DO i = 1, max_dom
      dynamics_grid_filename(i)   = ""

      ! "dynamics_parent_grid_id": This namelist parameter is
      ! necessary only for old ICON grids without a proper UUID
      ! attribute.
      dynamics_parent_grid_id(i)  = -1 
    ENDDO

    ! For the coarsest grid, we set the default "parent ID" to 0. This
    ! is necessary, since we cannot deduce this from another grid's
    ! UUID later on:
    dynamics_parent_grid_id(1) = 0

    radiation_grid_filename  = ""
    vertical_grid_filename   = " "
    vct_filename             = ""
      
    lfeedback   = .TRUE.
    ifeedback_type = 2
    start_time(:) = 0._wp
    end_time(:)   = DEFAULT_ENDTIME
    lplane         = .FALSE.
    is_plane_torus = .FALSE.
    l_limited_area = .FALSE.
    l_scm_mode     = .FALSE.
    lsep_grfinfo   = .FALSE.
    corio_lat   = 0.0_wp
    patch_weight= 0.0_wp
    lredgrid_phys = .FALSE.
    nexlevs_rrg_vnest = 8
    use_duplicated_connectivity = config_use_duplicated_connectivity
    use_dummy_cell_closure      = config_use_dummy_cell_closure

    create_vgrid = .FALSE.

    !----------------------------------------------------------------
    grid_rescale_factor   = 1.0_wp
    grid_angular_velocity = earth_angular_velocity
    lrescale_timestep     = .FALSE.
    lrescale_ang_vel      = .FALSE.

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
 !   IF (use_restart_namelists()) THEN
 !     funit = open_and_restore_namelist('grid_nml')
 !     READ(funit,NML=grid_nml)
 !     CALL close_tmpfile(funit)
 !   END IF

    !------------------------------------------------------------
    ! Read the namelist (done so far by all MPI processes)
    !------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('grid_nml', status=i_status)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, grid_nml)  ! write defaults to temporary text file
    END IF
    IF (i_status == POSITIONED) THEN
      READ (nnml, grid_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, grid_nml)  ! write settings to temporary text file
      END IF
    ENDIF
    CALL close_nml

    ! convert degrees in radiant for the Coriolis latitude
    corio_lat =  corio_lat/rad2deg

    ! Reset start and end times for global domain in order to avoid 
    ! interferences with output flow control
    start_time(1) = 0._wp
    end_time(1) = DEFAULT_ENDTIME

    !-----------------------------------------------------
    !  Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
    ! funit = open_tmpfile()
    ! WRITE(funit,NML=grid_nml)
    ! CALL store_and_close_namelist(funit, 'grid_nml')

      ! write the contents of the namelist to an ASCII file
      WRITE(nnml_output,nml=grid_nml)
    ENDIF


    !-----------------------------------------------------
    config_lfeedback         = lfeedback
    config_ifeedback_type    = ifeedback_type
    config_start_time        = start_time
    config_end_time          = end_time
    config_lplane            = lplane
    config_is_plane_torus    = is_plane_torus
    config_corio_lat         = corio_lat
    config_l_limited_area    = l_limited_area
    config_l_scm_mode        = l_scm_mode
    config_patch_weight      = patch_weight
    config_lredgrid_phys     = lredgrid_phys
    config_nexlevs_rrg_vnest = nexlevs_rrg_vnest
    config_use_duplicated_connectivity = use_duplicated_connectivity
    config_use_dummy_cell_closure      = use_dummy_cell_closure
    config_dynamics_grid_filename  = dynamics_grid_filename
    config_dynamics_parent_grid_id = dynamics_parent_grid_id
    config_radiation_grid_filename = radiation_grid_filename
    config_grid_rescale_factor     = grid_rescale_factor
    config_grid_angular_velocity   = grid_angular_velocity
    config_lrescale_timestep       = lrescale_timestep
    config_lrescale_ang_vel        = lrescale_ang_vel
    config_create_vgrid            = create_vgrid
    config_vertical_grid_filename  = vertical_grid_filename
    config_vct_filename            = vct_filename

    ! Throw a warning for deprecated parameters: ---------

    IF (lsep_grfinfo .EQV. .TRUE.) THEN
      ! default value has been modified; inform the user that this
      ! switch has no effect:
      IF (my_process_is_stdio()) &
        WRITE (0,*) "WARNING: Namelist switch 'lsep_grfinfo' is deprecated and will soon be removed!"
    END IF
          
  END SUBROUTINE read_grid_namelist
  !-----------------------------------------------------------------------
  
  
END MODULE mo_grid_nml
