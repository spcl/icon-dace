!> Contains conserved quantity types, definitions and methods for to-be-conserved quantities
!>
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
!>#### Contains conserved quantity types (CQTs) and the type for cqt collection (used on tiles)
!>
!> NOTE: currently only 2D variables are dealt with -- if also 3D variables need to be conserved,
!>       these will have to specially be dealt with in the concept and thus need to be
!>       identified with an onw cqt e.g. distinguish WATER_2D_CQ_TYPE and WATER_3D_CQ_TYPE?
!>       Furthermore, integration of pool-structure might need some additional work.
!>
MODULE mo_jsb_cqt_class
#ifndef __NO_JSBACH__

  USE mo_util,                ONLY: int2string
  USE mo_exception,           ONLY: finish
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: WATER_CQ_TYPE
  PUBLIC :: FLUX_C_CQ_TYPE, LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE
  PUBLIC :: Get_number_of_types, Get_cqt_name, t_jsb_consQuan, t_jsb_consQuan_p, max_cqt_name_length

  ENUM, BIND(C)
    ENUMERATOR ::                 &
      & WATER_CQ_TYPE = 0,        &
      & LIVE_CARBON_CQ_TYPE,      &
      & AG_DEAD_C_CQ_TYPE,        &
      & BG_DEAD_C_CQ_TYPE,        &
      & PRODUCT_CARBON_CQ_TYPE,   &
      & FLUX_C_CQ_TYPE,           &
      & LAST_CQ_TYPE ! needs to always be the last -- it is only used to determine the max number of types
  END ENUM

  !> Type used on tiles to collect conserved quantities from the process memories of the tile
  !
  TYPE :: t_jsb_consQuan
    INTEGER :: type_id = -1     !< one of the CQ_TYPEs found in mo_jsb_cqt_class
    INTEGER :: no_of_vars = 0   !< number of vars - required for allocation of cq_vars_2D
    INTEGER :: last_index_used = 0 !< last index used - required for filling of cq_vars_2D
    INTEGER, ALLOCATABLE :: associated_process(:) !< process to which each CQ in cq_vars_2D belongs
    TYPE(t_jsb_var_p), ALLOCATABLE :: cq_vars_2D(:)
        !< collection of all variables of a certain CQT potentially from different process memories
  END TYPE t_jsb_consQuan
  TYPE :: t_jsb_consQuan_p
    TYPE(t_jsb_consQuan), POINTER :: p => NULL()
  END TYPE t_jsb_consQuan_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_cqt_class'

  INTEGER, PARAMETER :: max_cqt_name_length = 25

CONTAINS

  ! ====================================================================================================== !
  !
  !> Returns the number of CQTs defined in the enumerator
  !
  FUNCTION Get_number_of_types() RESULT(last_type_id)
    INTEGER :: last_type_id

    last_type_id = LAST_CQ_TYPE

  END FUNCTION Get_number_of_types

  ! ====================================================================================================== !
  !
  !> Returns the name for this id - restricted length (max_cqt_name_length)
  !
  FUNCTION Get_cqt_name(id) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)  :: id
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_cqt_name'
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id)
      CASE (WATER_CQ_TYPE)
        return_value = 'water'
      CASE (LIVE_CARBON_CQ_TYPE)
        return_value = 'living_carbon'
      CASE (AG_DEAD_C_CQ_TYPE)
        return_value = 'ag_dead_carbon'
      CASE (BG_DEAD_C_CQ_TYPE)
        return_value = 'bg_dead_carbon'
      CASE (PRODUCT_CARBON_CQ_TYPE)
        return_value = 'product_carbon'
      CASE (FLUX_C_CQ_TYPE)
        return_value = 'carbon_flux'
      CASE DEFAULT
        CALL finish(TRIM(routine), 'No name specified for cq type of id '//int2string(id))
    END SELECT

  END FUNCTION Get_cqt_name

#endif
END MODULE mo_jsb_cqt_class
