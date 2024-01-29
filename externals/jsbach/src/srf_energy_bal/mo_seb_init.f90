!> Initialization of the the seb memory.
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
MODULE mo_seb_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_class,         ONLY: get_model
  USE mo_jsb_model_class,   ONLY: t_jsb_model

  !USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  !USE mo_jsb_class,         ONLY: get_model
  !USE mo_jsb_io_netcdf,     ONLY: t_input_file !, jsb_netcdf_open_input

  dsl4jsb_Use_processes SEB_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_memory(SEB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: seb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_init'

CONTAINS

  !
  !> Intialize soil process (after memory has been set up)
  !
  SUBROUTINE seb_init(tile)
    USE mo_jsb_time,        ONLY: timestep_in_days
    USE mo_pheno_constants, ONLY: PhenoParam

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init'

    dsl4jsb_Real2D_onDomain :: &
      & F_pseudo_soil_temp,    &
      & N_pseudo_soil_temp

    dsl4jsb_Def_memory(SEB_)

    model => get_model(tile%owner_model_id)

    CALL seb_init_bc(tile)
    !CALL seb_init_ic(tile)

    ! Initial values for the calculation of pseudo_soil_temp in SUBROUTINE calc_pseudo_soil_temp
    dsl4jsb_Get_memory(SEB_)


    IF (tile%contains_land .OR. model%config%use_tmx) THEN
      dsl4jsb_Get_var2D_onDomain(SEB_, N_pseudo_soil_temp) ! OUT
      dsl4jsb_Get_var2D_onDomain(SEB_, F_pseudo_soil_temp) ! OUT

      F_pseudo_soil_temp = EXP(- timestep_in_days(model_id=tile%owner_model_id)/PhenoParam%EG_SG%tau_pseudo_soil)
      N_pseudo_soil_temp = 1._wp / (1._wp - F_pseudo_soil_temp)

      !$ACC UPDATE DEVICE(F_pseudo_soil_temp, N_pseudo_soil_temp)
    END IF

  END SUBROUTINE seb_init

  SUBROUTINE seb_init_bc(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    !TYPE(t_jsb_model),       POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_bc'

    IF (tile%level > 0) CONTINUE ! only to avoid compiler warnings
    !model => get_model(tile%owner_model_id)

  END SUBROUTINE seb_init_bc

  SUBROUTINE seb_init_ic(tile)

    !USE mo_jsb_time, ONLY: get_time_interpolation_weights

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    ! TYPE(t_jsb_model),       POINTER :: model

    ! dsl4jsb_Def_config(SEB_)
    !dsl4jsb_Def_memory(SEB_)

    !REAL(wp), POINTER :: &
    !  & ptr_2D(:,  :)       !< temporary pointer

    !TYPE(t_input_file) :: input_file

    !INTEGER  :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':seb_init_ic'

    !model => get_model(tile%owner_model_id)

    ! Get seb config
    ! dsl4jsb_Get_config(SEB_)

    ! IF (.NOT. dsl4jsb_Config(SEB_)%active) RETURN

    ! Get seb memory of the tile
    !dsl4jsb_Get_memory(SEB_)

    ! Initialize parameters
    !
    ! IF (debug_on()) CALL message(TRIM(routine), 'Setting initial state of seb memory for tile '// &
      ! &                          TRIM(tile%name)//' from '//TRIM(dsl4jsb_Config(SEB_)%ic_filename))

    ! Initial surface temperature
    ! Is initialized in sse_init_ic together with soil temperatures

  END SUBROUTINE seb_init_ic

#endif
END MODULE mo_seb_init
