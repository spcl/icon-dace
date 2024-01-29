!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_dynamics_config

  USE mo_exception,             ONLY: message, print_value
  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: MAX_DOM
  USE mo_restart_nml_and_att,   ONLY: getAttributesForRestarting
  USE mo_key_value_store,       ONLY: t_key_value_store

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iequations, divavg_cntrwgt
  PUBLIC :: lcoriolis
  PUBLIC :: lmoist_thdyn
  PUBLIC :: ldeepatmo
  PUBLIC :: nold, nnow, nnew, nsav1, nsav2, nnow_rcf, nnew_rcf
  PUBLIC :: configure_dynamics

  !--------------------------------------------------------------------------
  ! Basic settings for the dynamical core 
  !--------------------------------------------------------------------------
  !TYPE :: t_dynamics_config

    ! namelist variables

    INTEGER  :: iequations      !< Choice of governing equation set
    REAL(wp) :: divavg_cntrwgt  !< Weight of central cell for divergence averaging
    LOGICAL  :: lcoriolis       !< if .TRUE., Coriolis force is switched on   
    LOGICAL  :: lmoist_thdyn    !< if .TRUE., include moisture terms in first law
    LOGICAL  :: ldeepatmo       !< if .TRUE., dynamical core assumes a deep atmosphere
                                !< instead of a shallow atmosphere

    ! derived variables

    INTEGER :: nold(MAX_DOM)      !< variables denoting time levels
    INTEGER :: nnow(MAX_DOM)      !< variables denoting time levels
    INTEGER :: nnew(MAX_DOM)      !< variables denoting time levels

    INTEGER :: nsav1(MAX_DOM)     !< Extra 'time levels' of prognostic variables
    INTEGER :: nsav2(MAX_DOM)     !< needed to compute boundary tendencies and
                                          !< feedback increments

    INTEGER :: nnow_rcf(MAX_DOM)  !< Extra time levels for reduced
    INTEGER :: nnew_rcf(MAX_DOM)  !< calling frequency (rcf)

  !END TYPE t_dynamics_config
  !>
  !!
  !TYPE(t_dynamics_config) :: dynamics_config(MAX_DOM)

CONTAINS
  !>
  !!
  SUBROUTINE configure_dynamics( ndom, ldynamics, ltransport )

    INTEGER,INTENT(IN) :: ndom
    LOGICAL,INTENT(IN) :: ldynamics, ltransport

    INTEGER :: jdom
    CHARACTER(2) :: sdom
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(LEN=*),PARAMETER :: routine='mo_dynamics_config:setup_dynamics_config'

    !------------------------
    ! Set time level indices

    CALL message(routine,'Set time level indices')

    CALL getAttributesForRestarting(restartAttributes)
    IF (restartAttributes%is_init) THEN
      ! Read time level indices from restart file.
      ! NOTE: this part will be modified later for a proper handling
      ! of multiple domains!!!
      DO jdom = 1,ndom
        WRITE (sdom, "(i2.2)") jdom
        CALL restartAttributes%get('nold_DOM'//sdom, nold(jdom))
        CALL restartAttributes%get('nnow_DOM'//sdom, nnow(jdom))
        IF (ldynamics) THEN
           ! This run is  with dynamics --> nnew must be different from nnow
           CALL restartAttributes%get('nnew_DOM'//sdom, nnew(jdom))
           IF (nnew(jdom) == nnow(jdom)) THEN
              ! previous run was with ldynamics = .FALSE. and therefore nnew was the same as nnow
              ! but this run is  with ldynamics = .TRUE.  and therefore nnew must be different from nnow
              IF (nnow(jdom) == 1) nnew(jdom) = 2
              IF (nnow(jdom) == 2) nnew(jdom) = 1
           END IF
        ELSE
           ! This run is  without dynamics --> nnew takes the same value as nnow
           nnew(jdom) = nnow(jdom)
        END IF
        CALL restartAttributes%get('nnow_rcf_DOM'//sdom, nnow_rcf(jdom))
        IF (ltransport) THEN
           ! This run is  with transport --> nnew_rcf must be different from nnow_rcf
           CALL restartAttributes%get('nnew_rcf_DOM'//sdom, nnew_rcf(jdom))
           IF (nnew_rcf(jdom) == nnow_rcf(jdom)) THEN
              ! previous run was with ltransport = .FALSE. and therefore nnew_rcf was the same as nnow_rcf
              ! but this run is  with ltransport = .TRUE.  and therefore nnew_rcf must be different from nnow_rcf
              IF (nnow_rcf(jdom) == 1) nnew_rcf(jdom) = 2
              IF (nnow_rcf(jdom) == 2) nnew_rcf(jdom) = 1
           END IF
        ELSE
           ! This run is  without transport --> nnew_rcf takes the same value as nnow_rcf
           nnew_rcf(jdom) = nnow_rcf(jdom)
        END IF
      END DO

    ELSE ! not isRestart

      ! time level indices for prognostic fields of the dynamical core
      IF (ldynamics) THEN
         ! if dynamics is active, then all necessary time levels are needed and
         ! the time level indices are swapped before the end of each time step
         nnow(:) = 1
         nnew(:) = 2
         nold(:) = 3
      ELSE
         ! if dynamics is not used, then use only a single time level
         nnow(:) = 1
         nnew(:) = 1
         nold(:) = 1
      END IF

      ! time level indices for prognostic fields of the transport scheme
      IF (ltransport) THEN
         ! if transport is active, then time levels "now" and "new" are needed and
         ! the time level indices are swapped before the end of each time step
         nnow_rcf(:) = 1
         nnew_rcf(:) = 2
      ELSE
         ! if transport is not used, then use only a single time level
         nnow_rcf(:) = 1
         nnew_rcf(:) = 1
      END IF

    END IF

    IF (ldynamics) THEN
       CALL message('ldynamics  =  .TRUE.','2 time levels used for progn. dynamics  variables -> nnew /= nnow')
    ELSE
       CALL message('ldynamics  = .FALSE.','1 time level  used for progn. dynamics  variables -> nnew = nnow')
    END IF
    IF (ltransport) THEN
       CALL message('ltransport =  .TRUE.','2 time levels used for progn. transport variables -> nnew_rcf /= nnow_rcf')
    ELSE
       CALL message('ltransport = .FALSE.','1 time level  used for progn. transport variables -> nnew_rcf = nnow_rcf')
    END IF
    DO jdom = 1,ndom
      WRITE (sdom, "(i2.2)") jdom
      CALL print_value('nnow    ('//sdom//')', nnow(jdom))
      CALL print_value('nnew    ('//sdom//')', nnew(jdom))
      CALL print_value('nold    ('//sdom//')', nold(jdom))
      CALL print_value('nnow_rcf('//sdom//')', nnow_rcf(jdom))
      CALL print_value('nnew_rcf('//sdom//')', nnew_rcf(jdom))
    END DO

    nsav1(:) = 0
    nsav2(:) = 4

  END SUBROUTINE configure_dynamics

END MODULE mo_dynamics_config
