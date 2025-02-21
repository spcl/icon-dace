! Contains the variables to set up the grid configuration
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

MODULE mo_grid_config
!-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message_text, finish, message
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_io_units,           ONLY: filename_max 
  USE mo_physical_constants, ONLY: earth_radius, earth_angular_velocity
  USE mo_parallel_config,    ONLY: division_method, division_file_name
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf, ONLY: &
      & nf_nowrite, nf_global, nf_noerr, nf_strerror, nf_open, nf_close, &
      & nf_get_att_int => nfx_get_att, &
      & nf_get_att_double => nfx_get_att

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_grid_configuration, get_grid_rescale_factor
  PUBLIC :: nroot, start_lev, n_dom, lfeedback, l_scm_mode,                  &
    &       lplane, is_plane_torus, corio_lat, l_limited_area, patch_weight, &
    &       lredgrid_phys, ifeedback_type, start_time, end_time
  PUBLIC :: grid_rescale_factor, grid_length_rescale_factor, &
     & grid_sphere_radius, grid_angular_velocity
  PUBLIC :: lrescale_timestep, lrescale_ang_vel
  PUBLIC :: namelist_grid_angular_velocity
  PUBLIC :: dynamics_grid_filename,  dynamics_parent_grid_id,     &
    &       radiation_grid_filename, nexlevs_rrg_vnest
  PUBLIC :: vertical_grid_filename, vct_filename, create_vgrid
  PUBLIC :: set_patches_grid_filename

! !   PUBLIC :: radiation_grid_distribution
  
  PUBLIC :: n_dom_start, max_childdom     
  PUBLIC :: n_phys_dom
  PUBLIC :: no_of_dynamics_grids
  PUBLIC :: use_duplicated_connectivity, use_dummy_cell_closure
  PUBLIC :: DEFAULT_ENDTIME
  ! ------------------------------------------------------------------------

  !> module name string
  CHARACTER(*), PARAMETER :: modname = 'mo_grid_config'

  REAL(wp), PARAMETER :: DEFAULT_ENDTIME = 1.e30_wp

  ! ------------------------------------------------------------------------
  !Configuration variables
  ! ------------------------------------------------------------------------
  INTEGER  :: nroot                    ! root division of initial edges
  INTEGER  :: start_lev                ! coarsest bisection level
  INTEGER  :: n_dom                    ! number of model domains, 1=global domain only 
  INTEGER  :: n_dom_start=1 
  INTEGER  :: max_childdom             ! type of feedback (incremental or relaxation)
  INTEGER  :: ifeedback_type
  REAL(wp) :: start_time(max_dom)      ! Time at which execution of a (nested) model domain starts
  REAL(wp) :: end_time(max_dom)        ! Time at which execution of a (nested) model domain terminates
  INTEGER  :: n_phys_dom=1             ! Number of physical domains, computed when reading the patches

  LOGICAL  :: lfeedback(max_dom)       ! specifies if feedback to parent grid is performed
  LOGICAL  :: lredgrid_phys(max_dom)   ! If set to .true. radiation is calculated on a reduced grid
  INTEGER  :: nexlevs_rrg_vnest        ! Maximum number of extra model layers used for calculating radiation if 
                                       ! a reduced radiation grid is combined with vertical nesting
  LOGICAL  :: l_limited_area           ! limited area setup where forcing comes in from sides 
  LOGICAL  :: l_scm_mode               ! SCM mode is designed for tests where all columns are identical
  LOGICAL  :: use_duplicated_connectivity  = .true.  ! if true, the zero connectivity is replaced by the last non-zero value
  LOGICAL  :: use_dummy_cell_closure = .false.  ! if true then create a dummy cell and connect it to cells and edges with no neigbor
   
!   INTEGER  :: radiation_grid_distribution   ! 0=do nothing
                                       ! 1=redistribute for radiaiton reading from file

  LOGICAL  :: lplane                   ! f-plane option
  LOGICAL  :: is_plane_torus           ! f-plane with doubly periodic boundary==> like a plane torus
  REAL(wp) :: corio_lat                ! Latitude, where the f-plane is located if 
                                       ! lplane or is_plane_torus=.true.

  REAL(wp) :: patch_weight(max_dom)    ! If patch_weight is set to a value > 0
                                       ! for any of the first level child patches,
                                       ! processor splitting will be performed

  REAL(wp) :: grid_rescale_factor = 1.0_wp
  REAL(wp) :: grid_length_rescale_factor = 1.0_wp
!   REAL(wp) :: grid_area_rescale_factor = 1.0_wp
  REAL(wp) :: grid_sphere_radius  = earth_radius
  LOGICAL  :: lrescale_timestep   = .FALSE.
  LOGICAL  :: lrescale_ang_vel    = .FALSE.
  REAL(wp) :: grid_angular_velocity  = earth_angular_velocity
  REAL(wp) :: namelist_grid_angular_velocity  = earth_angular_velocity

  CHARACTER(LEN=filename_max) :: dynamics_grid_filename(max_dom)
  INTEGER                     :: dynamics_parent_grid_id(max_dom)
  CHARACTER(LEN=filename_max) :: radiation_grid_filename

  LOGICAL    :: create_vgrid   ! switch if files containing vct_a, vct_b, z_ifc shall be created
  !> files containing vct_a, vct_b, z_ifc
  CHARACTER(LEN=filename_max) :: vertical_grid_filename(max_dom)
  CHARACTER(LEN=filename_max) :: vct_filename !< for reading in vct_a, vct_b

  INTEGER :: no_of_dynamics_grids  = 0

  ! -----------------------------------------------------------------------
  ! 2.0 Declaration of dependent variables
  ! -----------------------------------------------------------------------

CONTAINS

  !!  Initialization of variables that contain grid configuration
  SUBROUTINE init_grid_configuration
                                               
    !local variables
    INTEGER  :: jg, ncid
    LOGICAL  :: file_exists, lradiation_grid
    CHARACTER(*), PARAMETER :: routine = modname//"::init_grid_configuration"

    IF (no_of_dynamics_grids /= 0) &
      CALL finish( routine, 'should not be called twice')
    
    !-----------------------------------------------------------------------
    ! find out how many grids we have
    ! and check if they exist
    no_of_dynamics_grids  = 0

!     write(0,*) routine, TRIM(dynamics_grid_filename(1))
    
    jg=1
    DO WHILE (dynamics_grid_filename(jg) /= "")
      IF (my_process_is_stdio()) THEN
        INQUIRE (FILE=dynamics_grid_filename(jg), EXIST=file_exists)
        IF (.NOT. file_exists)   THEN
          WRITE (message_text,'(a,a)')  TRIM(dynamics_grid_filename(jg)), &
            " file does not exist"
          CALL finish( routine, TRIM(message_text))
        ENDIF
      ENDIF
      jg=jg+1
      IF (jg > max_dom) EXIT
    END DO
    no_of_dynamics_grids  = jg-1

    lradiation_grid = (LEN_TRIM(radiation_grid_filename) > 0)    
    IF (lradiation_grid .AND. (my_process_is_stdio())) THEN
      INQUIRE (FILE=radiation_grid_filename, EXIST=file_exists)
      IF (.NOT. file_exists)   THEN
        WRITE (message_text,'(a,a)')  TRIM(radiation_grid_filename), &
          " file does not exist"
        CALL finish( routine, TRIM(message_text))
      ENDIF
      DO jg = 1, no_of_dynamics_grids
        IF (TRIM(radiation_grid_filename) == TRIM(dynamics_grid_filename(jg))) THEN
          CALL finish( routine, "radiation_grid_filename must not be equal to dynamics_grid_filename!")
        END IF
      END DO
    END IF
    n_dom = no_of_dynamics_grids

    IF (no_of_dynamics_grids < 1) &
      CALL finish( routine, 'no dynamics grid is defined')
    
    IF (no_of_dynamics_grids > 1) CALL message(routine,'Warning: You are using one nproma for multiple grids!')

    ! get here the nroot, eventually it should be moved into the patch info
    CALL nf(nf_open(dynamics_grid_filename(1), nf_nowrite, ncid), routine)
    CALL get_gridfile_root_level(ncid, nroot, start_lev)
!     CALL get_gridfile_sphere_radius(ncid, grid_sphere_radius)
    grid_sphere_radius = earth_radius ! the grid-based radious is not used
    CALL nf(nf_close(ncid), routine)

    ! domain geometric properties
    IF ( grid_rescale_factor <= 0.0_wp ) grid_rescale_factor = 1.0_wp
    grid_sphere_radius = grid_sphere_radius * grid_rescale_factor
    grid_length_rescale_factor = grid_rescale_factor
    IF ( lrescale_ang_vel) THEN
       grid_angular_velocity = namelist_grid_angular_velocity / grid_rescale_factor
    ELSE
       grid_angular_velocity = namelist_grid_angular_velocity
    END IF

    IF (lradiation_grid) THEN
      n_dom_start = 0
    ELSE
      n_dom_start = 1
      lredgrid_phys(1) = .FALSE.    ! lredgrid_phys requires presence of patch0 => reset to false
    
      ! the division method starts from 0, shift if there's no 0 grid (ie no reduced radiation)
      DO jg = no_of_dynamics_grids-1, 0, -1
        division_method(jg+1)              = division_method(jg)
        division_file_name(jg+1)           = division_file_name(jg)
      ENDDO

    ENDIF
    
    !------------------------------------------------------------
    ! Reset lfeedback to false for all model domains if lfeedback(1) = false
    IF (.NOT. lfeedback(1)) lfeedback(2:max_dom) = .FALSE.
    
  END SUBROUTINE init_grid_configuration
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE get_gridfile_root_level( ncid, grid_root, grid_level )
    INTEGER,    INTENT(in)     :: ncid
    INTEGER,    INTENT(inout)  :: grid_root, grid_level
    CHARACTER(*), PARAMETER :: routine = modname//":get_gridfile_root_level"

    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_root', grid_root), routine)
    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_level', grid_level), routine)

  END SUBROUTINE get_gridfile_root_level
  !-------------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  SUBROUTINE get_gridfile_sphere_radius( ncid, sphere_radius )
    INTEGER,    INTENT(in)     :: ncid
    REAL(wp),   INTENT(out)    ::  sphere_radius

    INTEGER :: netcd_status

    netcd_status = nf_get_att_double(ncid, nf_global,'sphere_radius', &
        & sphere_radius )
    IF (netcd_status /= nf_noerr) sphere_radius = earth_radius

  END SUBROUTINE get_gridfile_sphere_radius
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  FUNCTION get_grid_rescale_factor( ) result(rescale_factor)
    REAL(wp) ::  rescale_factor

    IF (no_of_dynamics_grids < 1) &
      CALL finish( "get_grid_rescale_factor", 'no dynamics grid is defined')
    rescale_factor = grid_rescale_factor

  END FUNCTION get_grid_rescale_factor
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  INTEGER FUNCTION get_grid_root( patch_file )
    CHARACTER(*),    INTENT(in)  ::  patch_file   ! name of grid file
    CHARACTER(*), PARAMETER :: routine = modname//":get_grid_root"
    INTEGER :: ncid, grid_root

    CALL nf(nf_open(TRIM(patch_file), nf_nowrite, ncid), routine)
    CALL nf(nf_get_att_int(ncid, nf_global, 'grid_root', grid_root), routine)
    CALL nf(nf_close(ncid), routine)
    get_grid_root = grid_root
  END FUNCTION get_grid_root
  !-------------------------------------------------------------------------

  SUBROUTINE set_patches_grid_filename( grid_filename, grid_filename_grfinfo )

    CHARACTER(LEN=filename_max), INTENT(OUT) :: grid_filename(n_dom_start:)
    CHARACTER(LEN=filename_max), INTENT(OUT) :: grid_filename_grfinfo(n_dom_start:)
    ! local variables
    INTEGER :: jg, iind
    TYPE (t_keyword_list), POINTER :: keywords => NULL()
    CHARACTER(LEN=filename_max) :: grid_name

    !-----------------------------------------------------------------------
    DO jg = n_dom_start, n_dom

      CALL associate_keyword("<path>", TRIM(getModelBaseDir()), keywords)
      grid_name = ""
      IF (jg==0) THEN
        grid_name = TRIM(with_keywords(keywords, radiation_grid_filename))
      ELSE
        grid_name = TRIM(with_keywords(keywords, dynamics_grid_filename(jg)))
      ENDIF
      iind = INDEX(TRIM(grid_name),'.nc')
      grid_filename(jg)         = grid_name
      grid_filename_grfinfo(jg) = grid_name(1:iind-1)//"-grfinfo.nc"
    ENDDO

  END SUBROUTINE set_patches_grid_filename
  !-------------------------------------------------------------------------

END MODULE mo_grid_config
