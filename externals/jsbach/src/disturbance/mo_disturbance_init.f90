!> Initialization of the the disturbance memory.
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
MODULE mo_disturb_init
#ifndef __NO_JSBACH__

  !USE mo_kind,              ONLY: wp
  !USE mo_exception,         ONLY: message

  !USE mo_jsb_model_class,   ONLY: t_jsb_model
  !USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  !USE mo_jsb_io_netcdf,     ONLY: t_input_file !, jsb_netcdf_open_input
  !USE mo_jsb_time_iface,    ONLY: is_time_experiment_start, is_time_restart

  !dsl4jsb_Use_processes DISTURB_
  !dsl4jsb_Use_config(DISTURB_)
  !dsl4jsb_Use_memory(DISTURB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: disturb_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_init'

CONTAINS

  !
  !> Intialize disturb process (after memory has been set up)
  !
  SUBROUTINE disturb_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init'

    IF (tile%level > 0) CONTINUE ! avoid compiler warnings about dummy arguments not being used

    !dsl4jsb_Def_memory(DISTURB_)
    !dsl4jsb_Get_memory(DISTURB_)

    !CALL disturb_init_ic(tile)

    !CALL disturb_init_bc(tile)

  END SUBROUTINE disturb_init
!
!   SUBROUTINE disturb_init_bc(tile)
!     !-----------------------------------------------------------------------
!     !  ARGUMENTS
!     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!
!     !-----------------------------------------------------------------------
!     !  LOCAL VARIABLES
!     TYPE(t_jsb_model),       POINTER :: model
!     dsl4jsb_Def_memory(DISTURB_)
!
!     CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init_bc'
!
!     !-----------------------------------------------------------------------
!     ! CONTENT
!
!     model => get_model(tile%owner_model_id)
!
!     ! Get disturb memory of the tile
!     dsl4jsb_Get_memory(DISTURB_)
!
!     IF (.NOT. ASSOCIATED(disturb__mem)) RETURN
!
!     CALL message(TRIM(routine), 'Setting  disturb boundary conditions for tile '//TRIM(tile%name))
!
!     ! R: nothing to do in the moment ...
!
!   END SUBROUTINE disturb_init_bc
!
!   SUBROUTINE disturb_init_ic(tile)
!
!     !USE mo_disturb_constants,  ONLY: AlbedoCanopySnow_age, AlbedoCanopySnow_temp, AlbedoVisInitial, AlbedoNirInitial
!
!     CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
!
!     TYPE(t_jsb_model),         POINTER :: model
!     dsl4jsb_Def_config(DISTURB_)
!     dsl4jsb_Def_memory(DISTURB_)
!
!     !INTEGER :: idisturb
!     REAL(wp), POINTER :: &
!       & return_pointer  (:,  :) !, & !< temporary pointer
!       ! & return_pointer3d(:,:,:)
!     !LOGICAL :: &
!     !  read_cpools
!
!     TYPE(t_input_file) :: input_file
!     !INTEGER :: i, n1, n2
!
!     CHARACTER(len=*), PARAMETER :: routine = modname//':disturb_init_ic'
!
!     model => get_model(tile%owner_model_id)
!
!     ! Get disturb memory and config  of the tile
!     dsl4jsb_Get_memory(DISTURB_)
!     dsl4jsb_Get_config(DISTURB_)
!
!     IF (.NOT. dsl4jsb_Config(DISTURB_)%active) RETURN
!
!     CALL message(TRIM(routine), 'Reading disturb pools in memory for tile '//TRIM(tile%name))
!
!     ! Initialize physical state
!     !
!
!     input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(DISTURB_)%ic_filename), model%grid_id)
!
!     ! TBD: adapt names of variables in input files to names in memory
!
!     ! R: These variables should be 3D as they are tile dependent. However for THIS tile they would be 2D but in the input file
!     !    they are 3D!
!     ! R: jedenfalls wurden 3D variablen bisher nirgendwo über einlesen vom File initialisiert!
!     return_pointer => input_file%Read_3d(         &
!       & variable_name='c_green',           &
!       & fill_array = dsl4jsb_var_ptr(DISTURB_,c_green))
!     dsl4jsb_var3D_onDomain(DISTURB_,c_green) = &
!       MERGE(dsl4jsb_var3D_onDomain(DISTURB_,c_green), 0._wp, dsl4jsb_var3D_onDomain(DISTURB_,c_green) >= 0._wp)
!
!
!     CALL input_file%Close()
!
!   END SUBROUTINE disturb_init_ic
!
#endif
END MODULE mo_disturb_init
