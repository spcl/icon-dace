!> Interface between HD and the ocean, through a coupler
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
MODULE mo_interface_hd_ocean
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_coupling_config,     ONLY: is_coupled_run
#ifdef YAC_coupling
  USE yac,                    ONLY: yac_fdef_field, yac_fput, &
       &                            yac_fdef_mask,            &
       &                            YAC_LOCATION_CELL,        &
       &                            YAC_TIME_UNIT_ISO_FORMAT, &
       &                            YAC_ACTION_PUT_FOR_RESTART
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN
#endif

  USE mo_jsb_task_class,     ONLY: t_jsb_task_options

  dsl4jsb_Use_processes HYDRO_, HD_
  dsl4jsb_Use_config(HD_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(HD_)

  IMPLICIT NONE
  PRIVATE

#ifdef YAC_coupling
  PUBLIC :: jsb_fdef_hd_fields
#endif
  PUBLIC :: interface_hd_ocean

  CHARACTER(len=*), PARAMETER :: modname = 'mo_interface_hd_ocean'

  INTEGER :: yac_river_runoff_field_id

CONTAINS

#ifdef YAC_coupling
  SUBROUTINE jsb_fdef_hd_fields(comp_id, cell_point_ids, grid_id, n_patch_cells)

    USE mo_jsb_model_class,    ONLY: t_jsb_model
    USE mo_jsb_class,          ONLY: Get_model
    USE mo_jsb_grid_class,     ONLY: t_jsb_grid
    USE mo_jsb_grid,           ONLY: Get_grid

    INTEGER, INTENT(in) :: comp_id           ! yac component id
    INTEGER, INTENT(in) :: cell_point_ids(:) ! yac point is
    INTEGER, INTENT(in) :: grid_id           ! yac grid id
    INTEGER, INTENT(in) :: n_patch_cells     ! number of cells (compute + halo)
                                            ! n_patch_cells is <= nproma*nblks
    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid
    TYPE(t_jsb_model), POINTER    :: model

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    INTEGER :: cell_mask_ids(1)
    INTEGER :: model_id
    INTEGER :: jb, jc
    INTEGER :: nproma, nblks
    LOGICAL,  ALLOCATABLE :: is_valid(:)
    CHARACTER(len=MAX_TIMEDELTA_STR_LEN) :: modelTimeStep

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HD_)

    dsl4jsb_Real2D_onDomain :: hd_mask

    ! It is assumed that only one domain is active (jg=1) i.e. no nesting.
    ! Currently, this requirement is fulfilled since ocean coupling is only set up for one domain.
    model_id = 1

    model => Get_model(model_id)

    dsl4jsb_Get_config(HD_)
    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    CALL model%Get_top_tile(tile)
    dsl4jsb_Get_memory(HD_)

    grid  => Get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    CALL timedeltaToString(time_config%tc_dt_model, modelTimeStep)

    ! Define the HD mask for YAC.
    ! It shall contain ocean coast points only for
    ! source point mapping (source_to_target_map)

    ALLOCATE(is_valid(nproma*nblks))

    is_valid(:) = .FALSE.

    dsl4jsb_Get_var2D_onDomain(HD_, hd_mask)

    DO jb = 1, nblks
      DO jc = 1, nproma
        IF ( hd_mask(jc, jb) .EQ. 0.0 ) THEN
          is_valid((jb-1)*nproma+jc) = .TRUE.
        END IF
      END DO
    END DO

    CALL yac_fdef_mask (       &
      & grid_id,               &
      & n_patch_cells,         &
      & YAC_LOCATION_CELL,     &
      & is_valid,              &
      & cell_mask_ids(1) )

    CALL yac_fdef_field_mask ( &
      & 'river_runoff',        &
      & comp_id,               &
      & cell_point_ids,        &
      & cell_mask_ids(1),      &
      & 1,                     &
      & 1,                     &
      & modelTimeStep,         &
      & YAC_TIME_UNIT_ISO_FORMAT, &
      & yac_river_runoff_field_id )

    DEALLOCATE(is_valid)

  END SUBROUTINE jsb_fdef_hd_fields
#endif

  SUBROUTINE interface_hd_ocean(tile, options)

    USE mo_fortran_tools, ONLY: init

    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile
    TYPE(t_jsb_task_options),   INTENT(in) :: options

    ! Local variables
    !
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HD_)
    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':interface_hd_ocean'

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onDomain ::discharge_ocean

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells
    REAL(wp), ALLOCATABLE :: buffer(:,:)
    INTEGER               :: collection_size, nbr_pointsets, nbr_subdomains
    INTEGER               :: info, ierror !< return values form cpl_put/get calls
    INTEGER               :: iblk, n, nn, nlen

    n =options%nc ! only to avoid compiler warnings

    IF ( .NOT. is_coupled_run() ) RETURN

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HD_)
    IF (.NOT. dsl4jsb_Config(HD_)%active) RETURN

    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'HD model works on root tile only!')

    grid => Get_grid(model%grid_id)

    ! Get reference to variables in memory
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onDomain(HYDRO_, discharge_ocean) ! IN

    ! Send river discharge
    ! --------------------

    write_coupler_restart = .FALSE.

    ALLOCATE(buffer(grid%nproma * grid%nblks, 1))
    !$ACC DATA CREATE(buffer)
!$OMP PARALLEL
    CALL init(buffer(:,:), lacc=.TRUE.)
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(iblk, nn, nlen, n)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO iblk = 1, grid%nblks
      nn = (iblk - 1) * grid%nproma
      IF (iblk /= grid%nblks) THEN
       nlen = grid%nproma
      ELSE
       nlen = grid%npromz
      END IF
      !$ACC LOOP GANG VECTOR
      DO n = 1, nlen
        buffer(nn+n,1) = discharge_ocean(n,iblk)
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC UPDATE HOST(buffer) ASYNC(1)
    !$ACC WAIT(1)
!$OMP END PARALLEL DO

    nbr_hor_cells = grid%ntotal
    collection_size = 1
    nbr_pointsets   = 1
    nbr_subdomains  = 1
#ifdef YAC_coupling
    CALL yac_fput(yac_river_runoff_field_id, nbr_hor_cells, collection_size, buffer, info, ierror)
    IF ( info == YAC_ACTION_PUT_FOR_RESTART ) write_coupler_restart = .TRUE.
#endif

    IF ( write_coupler_restart ) THEN
       CALL message(TRIM(routine), 'YAC says it is put for restart')
    ENDIF

    !$ACC END DATA
    DEALLOCATE(buffer)

  END SUBROUTINE interface_hd_ocean

#endif
END MODULE mo_interface_hd_ocean
