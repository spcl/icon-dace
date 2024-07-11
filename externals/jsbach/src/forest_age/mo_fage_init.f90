!> fage (forest age) memory initialisation
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
!>#### Initialisation of fage (forest age) variables
!>
MODULE mo_fage_init
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_vgrid, t_jsb_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input

  dsl4jsb_Use_processes FAGE_, PHENO_
  dsl4jsb_Use_config(PHENO_)
  dsl4jsb_Use_config(FAGE_)
  dsl4jsb_Use_memory(FAGE_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fage_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_init'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize fage process
  !
  SUBROUTINE fage_init(tile)
    !-----------------------------------------------------------------------
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    !-----------------------------------------------------------------------
    CHARACTER(len=*), PARAMETER :: routine = modname//':fage_init'
    !-----------------------------------------------------------------------

    CALL fage_init_ic(tile)

  END SUBROUTINE fage_init

  ! ====================================================================================================== !
  !
  !> Initialize initial conditions for fage process
  !
  SUBROUTINE fage_init_ic(tile)

    USE mo_fage_util,      ONLY: get_mean_age
    USE mo_jsb_parallel,   ONLY: Is_enabled_openmp
    !-----------------------------------------------------------------------
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    !-----------------------------------------------------------------------
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    REAL(wp), POINTER :: tile_fract(:,:)

    TYPE(t_input_file) :: input_file
    TYPE(t_jsb_vgrid), POINTER :: forest_age

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    dsl4jsb_Real3D_onDomain :: fract_per_age

    CHARACTER(len=*), PARAMETER :: routine = modname//':fage_init_ic'
    !-----------------------------------------------------------------------

    model => get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    dsl4jsb_Get_config(PHENO_)
    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var3D_onDomain(FAGE_, fract_per_age)

    ! Assert: does currently not run with openmp
    IF (Is_enabled_openmp()) THEN
      CALL finish(TRIM(routine), "Currently forest age class implementation does not run with openmp!")
    END IF

    ! Assert: l_forestRegrowth
    IF (.NOT. dsl4jsb_Config(PHENO_)%l_forestRegrowth) THEN
      CALL finish(TRIM(routine), "Current forest age class implementation requires l_forestRegrowth to be true!")
    END IF

    ! Distribution of age fractions - needs to comply with the fractions in the age classes (mo_jsb_tile - Init_tile)
    IF (dsl4jsb_Config(FAGE_)%init_ac_fracts_scheme == 'allFirst') THEN
      ! all fracts to youngest age
      dsl4jsb_var3D_onDomain(FAGE_, fract_per_age) = 0._wp

      ALLOCATE(tile_fract(grid%nproma, grid%nblks))
      CALL tile%Get_fraction(fract=tile_fract(:,:))
      dsl4jsb_var_ptr(FAGE_, fract_per_age) (:,1,:) = tile_fract(:,:)
      DEALLOCATE(tile_fract)

    ELSE
      CALL finish(TRIM(routine), 'Currently only "allFirst" implemented for the distribution of age fractions and age classes.')
    END IF

    ! input_file = jsb_netcdf_open_input(dsl4jsb_Config(FAGE_)%ic_filename, model%grid_id)

    IF (tile%Has_children()) THEN
      dsl4jsb_var2D_onDomain(FAGE_, mean_age) = get_mean_age(dsl4jsb_var3D_onDomain(FAGE_, fract_per_age))
    ELSE
      !Only one age -> for now initialise age with zero later e.g. from file
      dsl4jsb_var2D_onDomain(FAGE_, mean_age) = 0._wp
    END IF

  END SUBROUTINE fage_init_ic

#endif
END MODULE mo_fage_init
