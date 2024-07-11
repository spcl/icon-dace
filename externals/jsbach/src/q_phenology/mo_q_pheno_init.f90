!> QUINCY phenology variables init
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
!>#### initialization of phenology memory variables using, e.g., ic & bc input files
!>
MODULE mo_q_pheno_init
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_grid,            ONLY: Get_grid
  USE mo_jsb_class,           ONLY: get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract

  dsl4jsb_Use_processes Q_PHENO_
  dsl4jsb_Use_config(Q_PHENO_)
  dsl4jsb_Use_memory(Q_PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: q_pheno_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_init'

CONTAINS

  ! ======================================================================================================= !
  !> run quincy phenology init
  !>
  !>
  SUBROUTINE q_pheno_init(tile)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':q_pheno_init'

    model => Get_model(tile%owner_model_id)

    ! init boundary conditions
    CALL q_pheno_init_bc(tile)
  END SUBROUTINE q_pheno_init

  ! ======================================================================================================= !
  !> Intialize pheno variables for boundary conditions (bc)
  !>
  !>
  SUBROUTINE q_pheno_init_bc(tile)

    ! see also the USE statements in the module header
    USE mo_q_pheno_constants,     ONLY: ievergreen
    USE mo_jsb_impl_constants,    ONLY: false, true

    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),       POINTER  :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':q_pheno_init_bc'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_PHENO_)
    dsl4jsb_Def_memory(Q_PHENO_)
    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onDomain    :: growing_season
    dsl4jsb_Real2D_onDomain    :: lai_max
    ! ---------------------------
    ! 0.4 Debug Option
    IF (debug_on()) CALL message(routine, 'Setting boundary conditions of pheno memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model => get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(Q_PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(Q_PHENO_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onDomain(Q_PHENO_, growing_season)
    dsl4jsb_Get_var2D_onDomain(Q_PHENO_, lai_max)

    !> 1.0 init
    !>
    !>   differs between QS and IQ
    !>
    ! init growing_season, lai_max, lai
    ! work with lctlib only if the present tile is a pft
    IF (tile%lcts(1)%lib_id /= 0) THEN
#ifdef __QUINCY_STANDALONE__
      !> QS - QUINCY standalone
      !>
      ! set lai_max from the site specific LAI from all_site_list.dat (saved in pheno_config:lai_max) and check for sensible values
      IF (dsl4jsb_Config(Q_PHENO_)%lai_max > 0.5_wp .AND. dsl4jsb_Config(Q_PHENO_)%lai_max < 8.0_wp) THEN
        lai_max(:,:) = dsl4jsb_Config(Q_PHENO_)%lai_max
      ELSE
        lai_max(:,:) = 6.0_wp
      END IF
#else
      !> IQ - QUINCY in ICON-Land
      !>
      ! lai_max
      ! for jsbach4 a lctlib parameter 'dsl4jsb_Lctlib_param(MaxLAI)' is available
      lai_max(:,:) = 6.0_wp
#endif
      ! growing season flag, set depending on whether this is a evergreen plant (TRUE) or not (FALSE)
      ! NOTE: false and true are of type REAL(wp)
      IF(dsl4jsb_Lctlib_param(phenology_type) == ievergreen) THEN
        growing_season(:,:) = true
      ELSE
        growing_season(:,:) = false
      ENDIF
    ELSE
      ! default for non-PFT tiles
      lai_max(:,:)        = 1.0_wp
      growing_season(:,:) = false
    ENDIF
  END SUBROUTINE q_pheno_init_bc

#endif
END MODULE mo_q_pheno_init
