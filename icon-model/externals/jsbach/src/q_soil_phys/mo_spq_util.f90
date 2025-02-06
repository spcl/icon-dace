!> QUINCY helper routines for soil-physics-quincy
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
!>#### various helper routines for the soil-physics-quincy process
!>
MODULE mo_spq_util
#ifndef __NO_QUINCY__

  USE mo_kind,          ONLY: wp
  USE mo_exception,     ONLY: finish, message, message_text
  USE mo_jsb_control,   ONLY: debug_on

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_qmax_texture
  PUBLIC :: reset_spq_fluxes           ! called in the land model interface prior to all Tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION calc_qmax_texture
  !!
  !! Output: sorption capacity for organic material
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_qmax_texture(qmax_fine_particle, silt_sl, clay_sl) RESULT(qmax_texture)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: qmax_fine_particle    !< maximum sorption capacity of fine soil particle [mol / kg fine particle (silt or clay)]
    REAL(wp), INTENT(in) :: silt_sl               !< silt content of the mineral soil [kg silt / kg mineral soil]
    REAL(wp), INTENT(in) :: clay_sl               !< clay content of the mineral soil [kg clay / kg mineral soil]
    REAL(wp)             :: qmax_texture          !< maximum sorption capacity of the soil [mol / kg]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_texture'

    ! maximum sorption capacity per kg mineral soil (mol/kg)
    qmax_texture = qmax_fine_particle * (silt_sl + clay_sl)

  END FUNCTION calc_qmax_texture


  !-----------------------------------------------------------------------------------------------------
  !> init/reset soil fluxes (with/to zero)
  !!
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_spq_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SPQ_, HYDRO_, VEG_, HD_, Q_RAD_, Q_ASSIMI_, Q_PHENO_, SB_
    USE mo_jsb_math_constants,    ONLY: zero

    ! Use of process memories
    dsl4jsb_Use_memory(SPQ_)        !  USE mo_spq_memory_class,      ONLY: t_spq_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    !TYPE(t_lctlib_element), POINTER       :: lctlib         !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_spq_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SPQ_)                          ! TYPE(t_spq_memory),      POINTER :: spq__mem
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_calculated(SPQ_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process memories
    dsl4jsb_Get_memory(SPQ_)           ! spq__mem => spq_memory
    model  => Get_model(tile%owner_model_id)
    !lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    dsl4jsb_var2D_onChunk(SPQ_, transpiration) = zero        !< transpiration in kg m-2 s-1
    dsl4jsb_var3D_onChunk(SPQ_, w_soil_freeze_flux) = zero
    dsl4jsb_var3D_onChunk(SPQ_, w_soil_melt_flux) = zero

  END SUBROUTINE reset_spq_fluxes

#endif
END MODULE mo_spq_util
