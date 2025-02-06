!> Contains process class
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
MODULE mo_jsb_process_class
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish
  USE mo_util_string,         ONLY: tolower
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_process_task_p
  USE mo_jsb_config_class,    ONLY: t_jsb_config
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process, t_jsb_process_p, Get_process_name, Get_process_id
  PUBLIC :: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, HD_, RAD_, ASSIMI_, PHENO_, CARBON_ , DISTURB_, FUEL_, &
    &       PPLCC_, FAGE_, ALCC_, FLCC_, WLCC_, NLCC_, VEG_, SB_, SPQ_, &
    &       Q_PHENO_, Q_ASSIMI_, Q_RAD_
  PUBLIC :: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_, Get_action_name

  TYPE t_jsb_process
    INTEGER                       :: id = 0
    CHARACTER(len=SHORT_NAME_LEN) :: name = ''
    CLASS(t_jsb_config), POINTER  :: config => NULL()
    TYPE(t_jsb_process_task_p), ALLOCATABLE :: task_list(:)
    INTEGER                       :: owner_model_id
    LOGICAL                       :: l_changes_fractions = .FALSE.  !< Does this process change tile fractions?
    LOGICAL                       :: has_memory = .TRUE.  !< Does this process has an own memory module
  CONTAINS
    PROCEDURE :: Configure => Configure_process
    PROCEDURE :: Register_task
    PROCEDURE :: Get_task
  END TYPE t_jsb_process

  ! The constructor for t_jsb_process is Create_process in mo_jsb_process_factory.

  TYPE t_jsb_process_p
    CLASS(t_jsb_process), POINTER :: p
  END TYPE t_jsb_process_p

  ENUM, BIND(C)
    ENUMERATOR :: A2L_=1, L2A_, SEB_, TURB_, SSE_, HYDRO_, HD_, RAD_, ASSIMI_, PHENO_,  CARBON_ , DISTURB_, &
      &           FUEL_, PPLCC_, FAGE_, ALCC_, FLCC_, WLCC_, NLCC_, VEG_, SB_, SPQ_, &
      &           Q_PHENO_, Q_ASSIMI_, Q_RAD_
  END ENUM

  ! Process actions
  ! It is important that AGGREGATE comes before the other actions
  ENUM, BIND(C)
    ENUMERATOR :: SKIP_=0, INHERIT_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_class'

CONTAINS

  FUNCTION Get_process_name(iproc) RESULT(return_value)

    INTEGER, INTENT(in) :: iproc
    CHARACTER(len=:), ALLOCATABLE :: return_value

    SELECT CASE(iproc)
    CASE (A2L_)
      return_value = 'a2l'
    CASE (L2A_)
      return_value = 'l2a'
    CASE (SEB_)
      return_value = 'seb'
    CASE (TURB_)
      return_value = 'turb'
    CASE (SSE_)
      return_value = 'sse'
    CASE (HYDRO_)
      return_value = 'hydro'
    CASE (RAD_)
      return_value = 'rad'
    CASE (Q_RAD_)
      return_value = 'q_rad'
    CASE (ASSIMI_)
      return_value = 'assimi'
    CASE (Q_ASSIMI_)
      return_value = 'q_assimi'
    CASE (PHENO_)
      return_value = 'pheno'
    CASE (Q_PHENO_)
      return_value = 'q_pheno'
#ifndef __NO_JSBACH_HD__
    CASE (HD_)
      return_value = 'hd'
#endif
    CASE (CARBON_)
      return_value = 'carbon'
    CASE (DISTURB_)
      return_value = 'disturb'
    CASE (FUEL_)
      return_value = 'fuel'
    CASE (PPLCC_)
      return_value = 'pplcc'
    CASE (FAGE_)
      return_value = 'fage'
    CASE (ALCC_)
      return_value = 'alcc'
    CASE (FLCC_)
      return_value = 'flcc'
    CASE (WLCC_)
      return_value = 'wlcc'
    CASE (NLCC_)
      return_value = 'nlcc'
    CASE (VEG_)
      return_value = 'veg'
    CASE (SB_)
      return_value = 'sb'
    CASE (SPQ_)
      return_value = 'spq'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION Get_process_name

  FUNCTION Get_process_id(name) RESULT(return_value)

    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value

    SELECT CASE(tolower(TRIM(name)))
    CASE ('a2l')
      return_value = A2L_
    CASE ('l2a')
      return_value = L2A_
    CASE ('seb')
      return_value = SEB_
    CASE ('turb')
      return_value = TURB_
    CASE ('sse')
      return_value = SSE_
    CASE ('hydro')
      return_value = HYDRO_
    CASE ('rad')
      return_value = RAD_
    CASE ('q_rad')
      return_value = Q_RAD_
    CASE ('assimi')
      return_value = ASSIMI_
    CASE ('q_assimi')
      return_value = Q_ASSIMI_
    CASE ('pheno')
      return_value = PHENO_
    CASE ('q_pheno')
      return_value = Q_PHENO_
#ifndef __NO_JSBACH_HD__
    CASE ('hd')
      return_value = HD_
#endif
    CASE ('carbon')
      return_value = CARBON_
    CASE ('disturb')
      return_value = DISTURB_
    CASE ('fuel')
      return_value = FUEL_
    CASE ('pplcc')
      return_value = PPLCC_
    CASE ('fage')
      return_value = FAGE_
    CASE ('alcc')
      return_value = ALCC_
    CASE ('flcc')
      return_value = FLCC_
    CASE ('wlcc')
      return_value = WLCC_
    CASE ('nlcc')
      return_value = NLCC_
    CASE ('veg')
      return_value = VEG_
    CASE ('sb')
      return_value = SB_
    CASE ('spq')
      return_value = SPQ_
    CASE DEFAULT
      return_value = 0
    END SELECT

  END FUNCTION Get_process_id

  SUBROUTINE Register_task(this, task)

    CLASS(t_jsb_process),     INTENT(inout) :: this
    CLASS(t_jsb_process_task), TARGET, INTENT(in)    :: task

    TYPE(t_jsb_process_task_p), ALLOCATABLE :: task_temp(:)
    INTEGER :: n

    CHARACTER(len=*), PARAMETER :: routine = modname//':Register_task'

    IF (this%id /= task%process_id) &
      CALL finish(TRIM(routine), 'Task is for another process')

    IF (.NOT. ALLOCATED(this%task_list)) THEN
      ALLOCATE(this%task_list(1))
      this%task_list(1) = t_jsb_process_task_p(task)
    ELSE
      n = SIZE(this%task_list)
      ALLOCATE(task_temp(n+1))
      task_temp(1:n) = this%task_list
      task_temp(n+1) = t_jsb_process_task_p(task)
      CALL MOVE_ALLOC(task_temp, this%task_list)
    END IF

    IF (debug_on()) CALL message(TRIM(routine), 'Registered task "'//TRIM(task%name)// &
      &                                         '" for process "'//TRIM(this%name)//'".')

  END SUBROUTINE Register_task

  FUNCTION Get_task(this, name) RESULT(return_ptr)

    CLASS(t_jsb_process), INTENT(inout) :: this
    CHARACTER(LEN=*),              INTENT(in)    :: name
    CLASS(t_jsb_process_task), POINTER           :: return_ptr

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_task'

    return_ptr => NULL()
    DO i=1,SIZE(this%task_list)
      !IF (.NOT. ASSOCIATED())
      IF (TRIM(this%task_list(i)%p%name) == TRIM(name) ) THEN
        return_ptr => this%task_list(i)%p
        EXIT
      END IF
    END DO

    IF (.NOT.(ASSOCIATED(return_ptr))) &
      & CALL finish(TRIM(routine), 'Task "'//TRIM(name)//'" not found for process "'// &
      &                            TRIM(this%name)//'".')

  END FUNCTION Get_task

  SUBROUTINE Configure_process(this)

    CLASS(t_jsb_process), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':Configure_process'

    IF (.NOT. ASSOCIATED(this%config)) RETURN

    CALL this%config%Init()

  END SUBROUTINE Configure_process

  FUNCTION Get_action_name(iact) RESULT(return_value)

    INTEGER, INTENT(in) :: iact
    CHARACTER(len=:), ALLOCATABLE :: return_value

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_action_name'

    SELECT CASE(iact)
    CASE (SKIP_)
      return_value = 'SKIP'
    CASE (AGGREGATE_)
      return_value = 'AGGREGATE'
    CASE (ON_LEAFS_)
      return_value = 'ON_LEAFS'
    CASE(ON_TILE_)
      return_value = 'ON_TILE'
    CASE (ON_SUBTREE_)
      return_value = 'ON_SUBTREE'
    CASE (INHERIT_)
      return_value = 'INHERIT'
    END SELECT

  END FUNCTION Get_action_name

#endif
END MODULE mo_jsb_process_class
