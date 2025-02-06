!> Contains JSBACH grid structure and methods
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_grid_class
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_exception,          ONLY: message_text, message, finish
  USE mo_jsb_domain_iface,   ONLY: t_patch, get_host_patch_id,                    &
                                   get_grid_filename_domain => get_grid_filename, &
                                   get_ntotal_domain        => get_ntotal,        &
                                   get_ntotal_g_domain      => get_ntotal_g,      &
                                   get_dims_g_domain        => get_dims_g,        &
                                   get_nlat_g_domain        => get_nlat_g,        &
                                   get_nproma_domain        => get_nproma,        &
                                   get_nblks_domain         => get_nblks
  USE mo_jsb_grid_iface,     ONLY: get_lon_domain           => get_lon,           &
                                   get_lat_domain           => get_lat,           &
                                   get_area_domain          => get_area
  USE mo_jsb_impl_constants, ONLY: SHORT_NAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_grid, t_jsb_grid_ptr, t_jsb_vgrid, t_jsb_vgrid_ptr
  PUBLIC :: new_grid, delete_grid, new_vgrid, delete_vgrid

  TYPE t_jsb_grid
    INTEGER                       :: id = 0
    CHARACTER(LEN=SHORT_NAME_LEN) :: name = ''
    CHARACTER(LEN=filename_max)   :: filename = ''
    INTEGER                       :: ntotal = 0          !< Number of cells on this PE
    INTEGER                       :: nland = 0
    REAL(wp), POINTER             :: lon(:,:) => NULL()  !< Longitudes of grid in degrees
    REAL(wp), POINTER             :: lat(:,:) => NULL()  !< Latitudes of grid in degrees
    REAL(wp), POINTER             :: coslat(:,:) => NULL()
    REAL(wp), POINTER             :: sinlat(:,:) => NULL()
    REAL(wp), POINTER             :: coslon(:,:) => NULL()
    REAL(wp), POINTER             :: sinlon(:,:) => NULL()
    REAL(wp), POINTER             :: area(:,:) => NULL()
    LOGICAL,  ALLOCATABLE         :: lsm(:,:)
    REAL(wp), ALLOCATABLE         :: lsf(:,:)
    INTEGER                       :: nblks = 1
    INTEGER                       :: nproma = 0
    INTEGER                       :: npromz = 0
    INTEGER                       :: ntotal_g              !< Total number of cells in grid
    INTEGER                       :: dims_g(2)             !< Dimensions of global grid
    INTEGER                       :: nlat_g = 0            !< Number of latitudes in global grid (effective for ICON)
    LOGICAL                       :: configured = .FALSE.
    CHARACTER(LEN=:), ALLOCATABLE :: type
    TYPE(t_patch), POINTER        :: patch => NULL()
    INTEGER                       :: host_patch_id
  CONTAINS
    PROCEDURE :: Get_nblks
    PROCEDURE :: Get_nproma
    PROCEDURE :: Get_npromz
    PROCEDURE :: Get_blk_start
    PROCEDURE :: Get_blk_end
    PROCEDURE :: Get_col_start
    PROCEDURE :: Get_col_end
    PROCEDURE :: Print             => Print_grid
  END TYPE t_jsb_grid

  TYPE t_jsb_grid_ptr
    TYPE(t_jsb_grid), POINTER :: p => NULL()
  END TYPE t_jsb_grid_ptr

  TYPE t_jsb_vgrid
    INTEGER ::           id           = 0
    CHARACTER(len=20) :: name         = ''
    CHARACTER(len=50) :: longname     = ''
    CHARACTER(len=10) :: units        = ''
    INTEGER :: n_levels               = 0
    REAL(wp), ALLOCATABLE :: levels(:)
    REAL(wp), ALLOCATABLE :: lbounds(:)
    REAL(wp), ALLOCATABLE :: ubounds(:)
    REAL(wp), ALLOCATABLE :: dz(:)
    INTEGER :: cdi_axis_type          = 0
    INTEGER :: ZaxisID                = 0   ! ICON: zaxis type (ZA_*), ECHAM: cdi zaxis id
    INTEGER :: echamZaxisIdx          = 0   ! ICON: N/A, ECHAM: levelindx for dimension
  END TYPE t_jsb_vgrid

  TYPE t_jsb_vgrid_ptr
    TYPE(t_jsb_vgrid), POINTER :: p => NULL()
  END TYPE t_jsb_vgrid_ptr

  INTERFACE new_grid
    MODULE PROCEDURE new_grid_from_file
    MODULE PROCEDURE new_grid_from_host
    ! MODULE PROCEDURE new_grid_from_lonlat ! Not implemented, yet!
  END INTERFACE new_grid

  INTEGER, SAVE :: ngrid_count = 0

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_grid_class'

CONTAINS

  !>
  !! Grid constructor

  FUNCTION new_grid_from_host(patch, type) RESULT(grid)

    USE mo_jsb_math_constants, ONLY: pi

    TYPE(t_patch),    TARGET, INTENT(in) :: patch
    CHARACTER(LEN=*),         INTENT(in) :: type
    TYPE(t_jsb_grid),         POINTER    :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_grid_from_host'

    ! CALL message(routine, 'starting construction of new JSBACH grid instance from host model')

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
!!$    new%name = name

    SELECT CASE (TRIM(type))
    CASE ('icon')
      new%patch => patch
    CASE ('echam')
      ALLOCATE(new%patch)
      new%patch = patch
    END SELECT
    new%type = TRIM(type)
    new%host_patch_id = get_host_patch_id(patch)

    new%filename = get_grid_filename_domain(patch)
    new%ntotal   = get_ntotal_domain       (patch)
    new%ntotal_g = get_ntotal_g_domain     (patch)
    new%dims_g(1:2) = get_dims_g_domain    (patch)
    new%nlat_g   = get_nlat_g_domain       (patch)
    new%nblks    = get_nblks_domain        (patch)
    new%nproma   = get_nproma_domain       (patch)
    new%npromz   = new%ntotal - new%nproma * (new%nblks - 1)
    new%lon     => get_lon_domain          (patch)
    new%lat     => get_lat_domain          (patch)
    ALLOCATE(new%coslat(new%nproma,new%nblks), new%sinlat(new%nproma,new%nblks), &
      &      new%coslon(new%nproma,new%nblks), new%sinlon(new%nproma,new%nblks))
    new%coslat(:,:) = COS(pi * new%lat(:,:) / 180._wp)
    new%sinlat(:,:) = SIN(pi * new%lat(:,:) / 180._wp)
    new%coslon(:,:) = COS(pi * new%lon(:,:) / 180._wp)
    new%sinlon(:,:) = SIN(pi * new%lon(:,:) / 180._wp)
    new%area    => get_area_domain         (patch)
    ALLOCATE(new%lsm(new%nproma,new%nblks), new%lsf(new%nproma,new%nblks))
    new%lsm(:,:) = .TRUE.
    new%lsf(:,:) = 1.0_wp

    new%configured = .TRUE.

    grid => new

    !$ACC ENTER DATA COPYIN(grid, grid%lon, grid%lat, grid%area, grid%lsf, grid%lsm, grid%sinlat, grid%coslat)

    ! CALL message(routine, 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_host

  FUNCTION new_grid_from_file(filename, name) RESULT(grid)

    ! USE mo_jsb_io_netcdf,      ONLY: ncfile_has_dim
    ! USE mo_jsb_parallel, ONLY: my_process_is_stdio

    CHARACTER(LEN=*), INTENT(in)    :: filename
    CHARACTER(LEN=*), INTENT(in)    :: name
    TYPE(t_jsb_grid), POINTER       :: grid

    TYPE(t_jsb_grid), POINTER       :: new

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_grid_from_file'

    CALL message(routine, 'starting construction of new JSBACH grid instance from file:')
    CALL message('     ', TRIM(filename))

    ALLOCATE(new)

    ngrid_count = ngrid_count + 1
    new%id = ngrid_count
    new%name = name
    new%filename = filename
    new%host_patch_id = 0

    ! IF (my_process_is_stdio()) THEN

    !   IF (ncfile_has_dim(TRIM(filename), 'lon_cell_centre')) THEN
    !     ! ICON grid
    !   ELSE
    !     ! ECHAM grid
    !   END IF

    ! END IF

    new%configured = .TRUE.

    grid => new

    CALL message(routine, 'construction of JSBACH grid completed.')

  END FUNCTION new_grid_from_file

  ! Not implemented, yet!
  ! FUNCTION new_grid_from_lonlat(lon, lat) RESULT(grid)

  !   REAL(wp)                        :: lon
  !   REAL(wp)                        :: lat
  !   TYPE(t_jsb_grid), POINTER       :: grid

  !   TYPE(t_jsb_grid), POINTER       :: new

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':new_grid_from_lonlat'

  !   ALLOCATE(new)

  !   grid => new

  ! END FUNCTION new_grid_from_lonlat

  SUBROUTINE delete_grid(self)

#ifdef HAVE_F2003
    TYPE(t_jsb_grid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_grid), POINTER :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':delete_grid'

    CALL message(routine, 'starting destruction of JSBACH grid instance.')

    IF (ALLOCATED(self%lsm)) DEALLOCATE(self%lsm)
    IF (ALLOCATED(self%lsf)) DEALLOCATE(self%lsf)
    IF (ASSOCIATED(self%lon)) DEALLOCATE(self%lon)
    IF (ASSOCIATED(self%lat)) DEALLOCATE(self%lat)
    IF (ASSOCIATED(self%coslat)) DEALLOCATE(self%coslat)
    IF (ASSOCIATED(self%sinlat)) DEALLOCATE(self%sinlat)
    IF (ASSOCIATED(self%coslon)) DEALLOCATE(self%coslon)
    IF (ASSOCIATED(self%sinlon)) DEALLOCATE(self%sinlon)

    DEALLOCATE(self)

    CALL message(routine, 'destruction of JSBACH grid completed.')

  END SUBROUTINE delete_grid

  !>
  !! Print method for grid instance.
  !!
  !! Subroutine to print description of grid instance.
  !!
  SUBROUTINE Print_grid(this)

    CLASS(t_jsb_grid), INTENT(in) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':Print_grid'

    CALL message(routine, 'print JSBACH grid description')

    WRITE(message_text,'(A,I2,A,A)') 'Grid id: ', this%id, ', grid name: ', TRIM(this%name)
    CALL message('  ', message_text)

    IF (.NOT. this%configured) THEN
      CALL message(routine, 'Grid not configured yet!')
      RETURN
    END IF

    IF (this%filename /= '') THEN
      WRITE(message_text,'(A)') 'Grid filename: '//TRIM(this%filename)
      CALL message('  ', message_text)
    END IF

    IF (this%host_patch_id > 0) THEN
      WRITE(message_text,'(A,I2)') 'Patch ID in host model: ', this%host_patch_id
      CALL message('  ', message_text)
    END IF

  END SUBROUTINE Print_grid

  INTEGER FUNCTION get_nblks(self)

#ifdef HAVE_F2003
    CLASS(t_jsb_grid), INTENT(in) :: self
#else
    CLASS(t_jsb_grid)             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nblks'

    get_nblks = self%nblks

  END FUNCTION get_nblks

  INTEGER FUNCTION get_nproma(self)

#ifdef HAVE_F2003
    CLASS(t_jsb_grid), INTENT(in) :: self
#else
    CLASS(t_jsb_grid)             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_nproma'

    get_nproma = self%nproma

  END FUNCTION get_nproma

  INTEGER FUNCTION get_npromz(self)

#ifdef HAVE_F2003
    CLASS(t_jsb_grid), INTENT(in) :: self
#else
    CLASS(t_jsb_grid)             :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_npromz'

    get_npromz = self%npromz

  END FUNCTION get_npromz

  INTEGER FUNCTION get_blk_start(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_start_domain => get_blk_start

    CLASS(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_start'

    get_blk_start = get_blk_start_domain(grid%patch)

  END FUNCTION get_blk_start

  INTEGER FUNCTION get_blk_end(grid)

    USE mo_jsb_domain_iface, ONLY: get_blk_end_domain => get_blk_end

    CLASS(t_jsb_grid), INTENT(in) :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_blk_end'

    get_blk_end = get_blk_end_domain(grid%patch)

  END FUNCTION get_blk_end

  INTEGER FUNCTION get_col_start(grid, iblk)

    USE mo_jsb_domain_iface, ONLY: get_col_start_domain => get_col_start

    CLASS(t_jsb_grid), INTENT(in) :: grid
    INTEGER,          INTENT(in) :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_start'

    get_col_start = get_col_start_domain(iblk, grid%patch)

  END FUNCTION get_col_start

  INTEGER FUNCTION get_col_end(grid, iblk)

    USE mo_jsb_domain_iface, ONLY: get_col_end_domain => get_col_end

    CLASS(t_jsb_grid), INTENT(in) :: grid
    INTEGER,          INTENT(in) :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_col_end'

    get_col_end = get_col_end_domain(iblk, grid%patch)

  END FUNCTION get_col_end

  FUNCTION new_vgrid(name, cdi_axis_type, length, longname, units, levels, lbounds, ubounds, dz) RESULT(new)


    CHARACTER(len=*),           INTENT(in) :: name
    INTEGER,                    INTENT(in) :: cdi_axis_type
    INTEGER,                    INTENT(in) :: length
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: longname
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: units
    REAL(wp),         OPTIONAL, INTENT(in) :: levels(length)
    REAL(wp),         OPTIONAL, INTENT(in) :: lbounds(length)
    REAL(wp),         OPTIONAL, INTENT(in) :: ubounds(length)
    REAL(wp),         OPTIONAL, INTENT(in) :: dz(length)

    TYPE(t_jsb_vgrid), POINTER :: new

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':new_vgrid'

    IF (TRIM(name) == '') CALL finish(routine, 'Must specify axis name')
    IF (length < 1) CALL finish(routine, TRIM(name)//' - length < 1')

    ALLOCATE(new)

    new%cdi_axis_type = cdi_axis_type
    new%name = TRIM(name)
    new%n_levels = length

    IF (PRESENT(longname)) THEN
      new%longname = TRIM(longname)
    ELSE
      new%longname = ''
    END IF
    IF (PRESENT(units)) THEN
      new%units  = TRIM(units)
    ELSE
      new%units = ''
    END IF

    ALLOCATE(new%levels(length), new%lbounds(length), new%ubounds(length), new%dz(length))
    new%levels (:) = 0._wp
    new%lbounds(:) = 0._wp
    new%ubounds(:) = 0._wp
    new%dz     (:) = 0._wp
    !$ACC ENTER DATA PCREATE(new%levels, new%lbounds, new%ubounds, new%dz)

    IF (PRESENT(levels))  new%levels (:) = levels(:)
    IF (PRESENT(lbounds)) new%lbounds(:) = lbounds(:)
    IF (PRESENT(ubounds)) new%ubounds(:) = ubounds(:)
    IF (PRESENT(dz)) THEN
      IF (ANY(dz(:) <= 0._wp)) CALL finish(routine, 'Dz must be larger than zero')
      new%dz     (:) = dz(:)
    END IF

    IF (PRESENT(lbounds) .AND. PRESENT(ubounds)) THEN
      IF (ANY(ubounds(:) < lbounds(:))) CALL finish(routine, 'Ubounds must be larger than Lbounds')
      IF (PRESENT(dz)) THEN
        IF (ANY(ABS((ubounds(:)-lbounds(:))-dz(:)) > EPSILON(1._wp))) CALL finish(routine, 'Dz != Ubounds-Lbounds')
      ELSE
        new%dz(:) = ubounds(:) - lbounds(:)
      END IF
      IF (.NOT. PRESENT(levels)) THEN
        new%levels(:) = 0.5_wp * (lbounds(:) + ubounds(:))
      END IF
    ELSE IF (PRESENT(lbounds)) THEN
      IF (PRESENT(dz)) THEN
        new%ubounds(:) = lbounds(:) + dz(:)
        IF (.NOT. PRESENT(levels)) THEN
          new%levels(:) = 0.5_wp * (lbounds(:) + new%ubounds(:))
        END IF
      ELSE IF (PRESENT(levels)) THEN
        new%ubounds(:) = levels(:) + (levels(:) - lbounds(:))
        new%dz(:) = new%ubounds(:) - lbounds(:)
      ELSE
        CALL finish(routine, 'Lbounds present ... specify either Levels or Dz!')
      END IF
    ELSE IF (PRESENT(ubounds)) THEN
      IF (PRESENT(dz)) THEN
        new%lbounds(:) = ubounds(:) - dz(:)
        IF (.NOT. PRESENT(levels)) THEN
          new%levels(:) = 0.5_wp * (new%lbounds(:) + ubounds(:))
        END IF
      ELSE IF (PRESENT(levels)) THEN
        new%lbounds(:) = levels(:) - (ubounds(:) - levels(:))
        new%dz(:) = ubounds(:) - new%lbounds(:)
      ELSE
        CALL finish(routine, 'Ubounds present ... specify either Levels or Dz!')
      END IF
    ELSE
      IF (PRESENT(levels) .AND. PRESENT(dz)) THEN
        CALL finish(routine, 'Levels and dz given ... please specify lbounds and/or ubounds')
      ELSE IF (PRESENT(dz)) THEN
        new%lbounds(1) = 0._wp
        DO i=2,length
          new%lbounds(i) = new%lbounds(i-1) + dz(i-1)
        END DO
        new%ubounds(:) = new%lbounds(:) + dz(:)
        new%levels(:) = 0.5_wp * (new%lbounds(:) + new%ubounds(:))
      ELSE IF (PRESENT(levels)) THEN
        IF (length > 1) THEN ! Multi-level axis
          new%lbounds(1) = 0._wp
          DO i=2,length
            new%lbounds(i) = levels(i-1) + (levels(i-1) - new%lbounds(i-1))
          END DO
          new%ubounds(:) = levels(:) + (levels(:) - new%lbounds(:))
          new%dz(:) = new%ubounds(:) - new%lbounds(:)
        ELSE ! Single level axis (level specified)
          new%lbounds(1) = levels(1)
          new%ubounds(1) = levels(1)
          new%dz(1)      = 0._wp
        END IF
      ELSE
        IF (length > 1) THEN ! Multi-level axis
          DO i=1,length
            new%levels(i) = REAL(i,wp)
            new%lbounds(i) = REAL(i,wp)
            new%ubounds(i) = REAL(i+1,wp)
            new%dz(i) = 1._wp
          END DO
        ELSE ! Single level axis (level not specified)
          new%levels(1) = 0._wp
          new%lbounds(1) = 0._wp
          new%ubounds(1) = 0._wp
          new%dz(1)      = 0._wp
        END IF
      END IF
    END IF
    !$ACC UPDATE DEVICE(new%levels) ASYNC(1) IF(present(levels))
    !$ACC UPDATE DEVICE(new%lbounds) ASYNC(1) IF(present(lbounds))
    !$ACC UPDATE DEVICE(new%ubounds) ASYNC(1) IF(present(ubounds))
    !$ACC UPDATE DEVICE(new%dz) ASYNC(1) IF(present(dz))
  END FUNCTION new_vgrid

  SUBROUTINE delete_vgrid(self)

    ! USE mo_jsb_io_iface, ONLY: Destroy_zaxis

#ifdef HAVE_F2003
    TYPE(t_jsb_vgrid), POINTER, INTENT(inout) :: self
#else
    TYPE(t_jsb_vgrid), POINTER :: self
#endif

    CHARACTER(len=*), PARAMETER :: routine = modname//':delete_vgrid'

    IF (self%id == 0) RETURN

    CALL message(routine, 'starting destruction of JSBACH vgrid instance.')

    !$ACC EXIT DATA DELETE(self%levels) IF(allocated(self%levels))
    !$ACC EXIT DATA DELETE(self%lbounds) IF(allocated(self%lbounds))
    !$ACC EXIT DATA DELETE(self%ubounds) IF(allocated(self%ubounds))
    !$ACC EXIT DATA DELETE(self%dz) IF(allocated(self%dz))
    IF (ALLOCATED(self%levels))  DEALLOCATE(self%levels)
    IF (ALLOCATED(self%lbounds)) DEALLOCATE(self%lbounds)
    IF (ALLOCATED(self%ubounds)) DEALLOCATE(self%ubounds)
    IF (ALLOCATED(self%dz))      DEALLOCATE(self%dz)

    DEALLOCATE(self)

    CALL message(routine, 'destruction of JSBACH vgrid completed.')

  END SUBROUTINE delete_vgrid

#endif
END MODULE mo_jsb_grid_class
