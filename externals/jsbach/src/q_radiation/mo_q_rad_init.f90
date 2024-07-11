!> QUINCY radiation variables init
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### initialization of radiation memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_rad_init
#ifndef __NO_QUINCY__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  USE mo_jsb_grid_class,    ONLY: t_jsb_grid
  USE mo_jsb_grid,          ONLY: Get_grid
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_io_netcdf,     ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io,            ONLY: missval

  dsl4jsb_Use_processes Q_RAD_
  dsl4jsb_Use_memory(Q_RAD_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: q_rad_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_rad_init'

CONTAINS

  !
  !> Intialize radiation process (after memory has been set up)
  !
  SUBROUTINE q_rad_init(tile)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':q_rad_init'

    IF (tile%contains_vegetation) THEN
        CALL q_rad_init_ic(tile)
    END IF
  END SUBROUTINE q_rad_init

  !-----------------------------------------------------------------------------------------------------
  !> Intialize radiation process (after memory has been set up)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE q_rad_init_ic(tile)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_q_rad_parameters,      ONLY: rfr_ratio_toc

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct

    CHARACTER(len=*), PARAMETER :: routine = modname//':q_rad_init_ic'

    dsl4jsb_Def_memory(Q_RAD_)

    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc
    dsl4jsb_Real2D_onDomain      :: rfr_ratio_boc_tvegdyn_mavg
    dsl4jsb_Real3D_onDomain      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onDomain      :: ppfd_shaded_tcnl_mavg_cl

    IF (.NOT. tile%Is_process_active(Q_RAD_)) RETURN

    IF (debug_on()) CALL message(routine, 'Setting initial conditions of radiation memory (quincy) for tile '// &
      &                                   TRIM(tile%name))

    dsl4jsb_Get_memory(Q_RAD_)

    dsl4jsb_Get_var2D_onDomain(Q_RAD_, rfr_ratio_boc)
    dsl4jsb_Get_var2D_onDomain(Q_RAD_, rfr_ratio_boc_tvegdyn_mavg)
    dsl4jsb_Get_var3D_onDomain(Q_RAD_, ppfd_sunlit_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(Q_RAD_, ppfd_sunlit_tcnl_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(Q_RAD_, ppfd_shaded_tfrac_mavg_cl)
    dsl4jsb_Get_var3D_onDomain(Q_RAD_, ppfd_shaded_tcnl_mavg_cl)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    rfr_ratio_boc(:,:)                  = rfr_ratio_toc
    rfr_ratio_boc_tvegdyn_mavg(:,:)     = rfr_ratio_boc(:,:)

    ppfd_sunlit_tfrac_mavg_cl(:,:,:)    = 150.0_wp
    ppfd_sunlit_tcnl_mavg_cl(:,:,:)     = 150.0_wp

    ppfd_shaded_tfrac_mavg_cl(:,:,:)    = 10.0_wp
    ppfd_shaded_tcnl_mavg_cl(:,:,:)     = 10.0_wp

  END SUBROUTINE q_rad_init_ic

#endif
END MODULE mo_q_rad_init
