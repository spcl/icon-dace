!> QUINCY model modes
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
!>#### defines QUINCY model modes and related helper routines
!>
MODULE mo_quincy_model_config
#ifndef __NO_QUINCY__

  USE mo_util,             ONLY: tolower
  USE mo_exception,        ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Get_quincy_model_config_name, Get_quincy_model_config_id
  PUBLIC :: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION

  !> Available quincy model ids
  !>
  ENUM, BIND(C)
    ENUMERATOR :: QLAND=1, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = ' mo_quincy_model_config'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Returns the QUINCY model name for this id
  !
  FUNCTION Get_quincy_model_config_name(id) RESULT(return_value)

    INTEGER, INTENT(in) :: id
    CHARACTER(len=:), ALLOCATABLE :: return_value

    SELECT CASE(id)
    CASE (QLAND)
      return_value = 'land'
    CASE (QPLANT)
      return_value = 'plant'
    CASE (QSOIL)
      return_value = 'soil'
    CASE (QCANOPY)
      return_value = 'canopy'
    CASE (Q_TEST_CANOPY)
      return_value = 'test_canopy'
    CASE (Q_TEST_RADIATION)
      return_value = 'test_radiation'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION Get_quincy_model_config_name

  ! ====================================================================================================== !
  !
  !> Returns the id for this QUINCY model name
  !
  FUNCTION Get_quincy_model_config_id(name) RESULT(return_value)

    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_quincy_model_config_id'

    SELECT CASE(tolower(TRIM(name)))
    CASE ('land')
      return_value = QLAND
    CASE ('plant')
      return_value = QPLANT
    CASE ('soil')
      return_value = QSOIL
    CASE ('canopy')
      return_value = QCANOPY
    CASE ('test_canopy')
      return_value = Q_TEST_CANOPY
    CASE ('test_radiation')
      return_value = Q_TEST_RADIATION
    CASE DEFAULT
      CALL finish(routine,'Invalid QUINCY model name. Try either land, plant, soil, canopy, test_canopy or test_radiation.')
    END SELECT

  END FUNCTION Get_quincy_model_config_id

#endif
END MODULE mo_quincy_model_config
