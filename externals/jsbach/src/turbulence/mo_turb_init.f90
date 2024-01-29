!> Initialization of the the turbulence memory.
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
MODULE mo_turb_init
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_grid_class,    ONLY: t_jsb_grid
  USE mo_jsb_grid,          ONLY: Get_grid
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on, jsbach_runs_standalone
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,            ONLY: missval

  dsl4jsb_Use_processes TURB_
  dsl4jsb_Use_config(TURB_)
  dsl4jsb_Use_memory(TURB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: turb_init

  TYPE t_turb_init_vars
    REAL(wp), POINTER ::  &
      & rough_m(:,:)
  END TYPE t_turb_init_vars

  TYPE(t_turb_init_vars) :: turb_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_init'

CONTAINS

  !
  !> Intialize turbulence process (after memory has been set up)
  !
  SUBROUTINE turb_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init'

    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL turb_read_init_vars(tile)
    END IF

    CALL turb_init_bc(tile)
    CALL turb_init_ic(tile)

    IF (tile%Is_last_process_tile(TURB_)) THEN
      CALL turb_finalize_init_vars()
    END IF

  END SUBROUTINE turb_init

  SUBROUTINE turb_finalize_init_vars

    DEALLOCATE( &
      & turb_init_vars%rough_m &
      & )

  END SUBROUTINE turb_finalize_init_vars

  SUBROUTINE turb_read_init_vars(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    dsl4jsb_Def_config(TURB_)

    REAL(wp), POINTER :: ptr_2D(:,:)       !< temporary pointer
    TYPE(t_input_file) :: input_file

    INTEGER :: nproma, nblks

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_read_init_vars'

    model => get_model(tile%owner_model_id)
    grid  => get_grid(model%grid_id)
    nproma = grid%Get_nproma()
    nblks  = grid%Get_nblks()

    ALLOCATE(turb_init_vars%rough_m(nproma, nblks))
    turb_init_vars%rough_m(:,:) = missval

    ! Get turb config
    dsl4jsb_Get_config(TURB_)

    IF (.NOT. dsl4jsb_Config(TURB_)%active) RETURN

    IF (debug_on()) CALL message(TRIM(routine), &
         'Reading turb init vars for tile from '//TRIM(dsl4jsb_Config(TURB_)%bc_filename))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(TURB_)%bc_filename), model%grid_id)

    ! Surface roughness length
    ptr_2D => input_file%Read_2d(         &
      & variable_name='roughness_length',         &
      & fill_array = turb_init_vars%rough_m)
    ptr_2D = MERGE(ptr_2D, 1._wp, ptr_2D > 0._wp)
    IF (dsl4jsb_Config(TURB_)%max_ini_rough_m > 0._wp) THEN
      ptr_2D = MIN(ptr_2D, dsl4jsb_Config(TURB_)%max_ini_rough_m)
                                  ! TODO: rough_m should not be larger than lowest model layer thickness, or even half of that !?
    END IF

    CALL input_file%Close()

  END SUBROUTINE turb_read_init_vars

  SUBROUTINE turb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init_bc'

    IF (tile%level > 0) CONTINUE ! only to avoid compiler warnings
    !model => get_model(tile%owner_model_id)

  END SUBROUTINE turb_init_bc

  SUBROUTINE turb_init_ic(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(TURB_)
    dsl4jsb_Def_memory(TURB_)

    CHARACTER(len=*), PARAMETER :: routine = modname//':turb_init_ic'

    model => get_model(tile%owner_model_id)

    ! Get turb config
    dsl4jsb_Get_config(TURB_)

    IF (.NOT. dsl4jsb_Config(TURB_)%active) RETURN

    ! Get turb memory of the tile
    dsl4jsb_Get_memory(TURB_)

    ! Initialize parameters
    !
    ! rough_m from ini file is not used in update_roughness, it's computed there (except for jsbach_lite). However,
    ! the atmosphere also reads this variable from this file and uses it to compute the turbulent exchange parameters for
    ! the first time step which are then passed to JSBACH.
    ! In the standalone model, these exchange parameters are computed in mo_forcing using rough_m and then
    ! passed to the first call of JSBACH. Therefore, we need rough_m from the ini file for the standalone model.
    ! We also need rough_m for runs without PFTs, e.g. jsbach_lite, on vegetation tiles that are not PFTs (e.g. veg_tile)
    IF (jsbach_runs_standalone() .OR. (tile%contains_vegetation .AND. tile%lcts(1)%lib_id == 0)) THEN

      IF (debug_on()) CALL message(TRIM(routine), &
          'Initializing turb memory for tile '//TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(TURB_)%bc_filename))

      ! Surface roughness length
      dsl4jsb_var2D_onDomain(TURB_,rough_m) = turb_init_vars%rough_m
      dsl4jsb_var2D_onDomain(TURB_,rough_m_star) = 1._wp / LOG(dsl4jsb_Config(TURB_)%blending_height / turb_init_vars%rough_m)**2

      ! Needed for tmx
      IF (.NOT. jsbach_runs_standalone()) & ! (TODO: doesn't work for standalone ... why?)
        & dsl4jsb_var2D_onDomain(TURB_,rough_h) = turb_init_vars%rough_m

      !$ACC UPDATE &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(TURB_,rough_m), dsl4jsb_var2D_onDomain(TURB_,rough_m_star)) &
      !$ACC   DEVICE(dsl4jsb_var2D_onDomain(TURB_,rough_h))

    END IF

  END SUBROUTINE turb_init_ic

#endif
END MODULE mo_turb_init
