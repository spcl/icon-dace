!
! Classes and functions for the turbulent mixing package (tmx)
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

MODULE mo_tmx_process_class

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: message, finish
  USE mo_fortran_tools, ONLY: init
  USE mtime,        ONLY: t_datetime => datetime
  USE mo_variable, ONLY: t_variable, bind_variable, allocate_variable
  USE mo_variable_list, ONLY: t_variable_list, t_variable_set
  USE mo_surrogate_class, ONLY: t_surrogate
  USE mo_tmx_field_class, ONLY: t_tmx_field, t_tmx_field_list, t_domain !bind_tmx_field
  USE mo_tmx_time_integration_class, ONLY: t_time_scheme








  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tmx_process, t_tmx_process_p

  TYPE, ABSTRACT, EXTENDS(t_surrogate) :: t_tmx_process
    ! PRIVATE
    CHARACTER(LEN=:),      ALLOCATABLE :: name         !< Process name
    TYPE(t_domain),        POINTER     :: domain       !< Spatial domain
    REAL(wp)                           :: dt           !< Time step
    TYPE(t_tmx_field_list)             :: states       !< State variables
    TYPE(t_variable_list)              :: tendencies   !< Tendency variables
    TYPE(t_variable_list)              :: new_states   !< New state variables
    CLASS(t_variable_set), ALLOCATABLE :: diagnostics  !< Diagnostic variables
    CLASS(t_variable_set), ALLOCATABLE :: config       !< Config variables
    CLASS(t_variable_set), ALLOCATABLE :: inputs       !< Input variables
    CLASS(t_time_scheme),  ALLOCATABLE :: time_scheme  !< Time integration scheme
    TYPE(t_tmx_process_p), POINTER     :: processes(:) => NULL() !< Subprocesses
  CONTAINS
    PROCEDURE                          :: Init => Init_process
    ! PROCEDURE(init_iface),     DEFERRED :: Init
    PROCEDURE                          :: Lock_variable_sets
    PROCEDURE                          :: Add_process
    PROCEDURE                          :: Add_state_r2d
    PROCEDURE                          :: Add_state_r3d
    PROCEDURE                          :: Add_state_shape_real
    GENERIC                            :: Add_state => Add_state_r2d, Add_state_r3d, Add_state_shape_real
    PROCEDURE                          :: Set_time_scheme
    PROCEDURE                          :: Step_forward
    PROCEDURE(compute_iface),      DEFERRED :: Compute !< Compute
    PROCEDURE(compute_iface),      DEFERRED :: Compute_diagnostics !< Compute diagnostics
    PROCEDURE(compute_diag_iface), DEFERRED :: Update_diagnostics !< Update diagnostics
    PROCEDURE                          :: Get_tendency_r2d
    PROCEDURE                          :: Get_tendency_r3d
    PROCEDURE                          :: Get_diagnostic_r2d
    PROCEDURE                          :: Get_diagnostic_r3d
  END TYPE t_tmx_process

  TYPE t_tmx_process_p
    CLASS(t_tmx_process), POINTER :: p
  END TYPE t_tmx_process_p

  ABSTRACT INTERFACE
    ! SUBROUTINE init_iface(this, config)
    !   IMPORT :: t_tmx_process, t_variable_list
    !   CLASS(t_tmx_process),  INTENT(inout) :: this
    !   TYPE(t_variable_list), INTENT(in)    :: config
    ! END SUBROUTINE
    SUBROUTINE compute_iface(this, datetime)
      IMPORT :: t_tmx_process, t_datetime
      CLASS(t_tmx_process),     INTENT(inout), TARGET :: this
      TYPE(t_datetime), OPTIONAL, INTENT(in),   POINTER :: datetime     !< date and time at beginning of time step
    END SUBROUTINE compute_iface
    SUBROUTINE compute_diag_iface(this)
      IMPORT :: t_tmx_process
      CLASS(t_tmx_process),     INTENT(inout), TARGET :: this
    END SUBROUTINE compute_diag_iface
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_process_class'
  
CONTAINS

  SUBROUTINE Init_process(this, dt, name, domain)

    CLASS(t_tmx_process), INTENT(inout)        :: this
    REAL(wp),             INTENT(in)           :: dt
    CHARACTER(len=*),     INTENT(in), OPTIONAL :: name 
    TYPE(t_domain),       POINTER,    OPTIONAL :: domain

    this%dt = dt

    IF (PRESENT(name)) THEN
      this%name = name
    ELSE
      this%name = 'Unnamed'
    END IF

    IF (PRESENT(domain)) THEN
      this%domain => domain
      
    END IF

    this%states      = t_tmx_field_list('states')
    !$ACC ENTER DATA COPYIN(this%states)
    this%new_states  = t_variable_list ('new states')
    !$ACC ENTER DATA COPYIN(this%new_states)
    ! this%tendencies  = t_tmx_field_list('tendencies')
    this%tendencies  = t_variable_list ('tendencies')
    !$ACC ENTER DATA COPYIN(this%tendencies)
    ! this%config      = t_variable_list ('config')
    ! this%inputs      = t_variable_list ('inputs')
    ! this%diagnostics = t_variable_list ('diagnostics')

  END SUBROUTINE Init_process

  RECURSIVE SUBROUTINE Lock_variable_sets(this)

    CLASS(t_tmx_process), INTENT(inout):: this

    INTEGER :: iproc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Lock_variable_sets'

    IF (ALLOCATED(this%config))      CALL this%config%Set_pointers()
    IF (ALLOCATED(this%inputs))      CALL this%inputs%Set_pointers()
    IF (ALLOCATED(this%diagnostics)) CALL this%diagnostics%list%allocator()
    IF (ALLOCATED(this%diagnostics)) CALL this%diagnostics%Set_pointers()

    IF (ASSOCIATED(this%processes)) THEN
      DO iproc=1,SIZE(this%processes)
        CALL this%processes(iproc)%p%Lock_variable_sets()
      END DO
    END IF

  END SUBROUTINE Lock_variable_sets

  SUBROUTINE Add_process(this, process)

    CLASS(t_tmx_process), INTENT(inout) :: this
    CLASS(t_tmx_process), TARGET, INTENT(in) :: process

    TYPE(t_tmx_process_p), POINTER :: processes_(:)
    INTEGER :: n_processes, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_process'

    IF (.NOT. ASSOCIATED(this%processes)) THEN
      ALLOCATE(this%processes(1))
      this%processes(1)%p => process
    ELSE
      n_processes = SIZE(this%processes)
      ALLOCATE(processes_(n_processes+1))
      DO i=1,n_processes
        processes_(i)%p => this%processes(i)%p
      END DO
      processes_(n_processes+1)%p => process
      DEALLOCATE(this%processes)
      this%processes => processes_
    END IF

    CALL message(routine, 'Added sub-process '//process%name//' to '//this%name)

  END SUBROUTINE Add_process

  SUBROUTINE Add_state_shape_real(this, name, type, dims)

    USE mo_variable,        ONLY: t_variable, bind_variable

    CLASS(t_tmx_process), INTENT(inout), TARGET :: this
    CHARACTER(len=*),     INTENT(in)            :: name
    INTEGER,              INTENT(in)            :: type
    INTEGER,              INTENT(in)            :: dims(:)

    TYPE(t_variable), POINTER :: tv
    INTEGER :: ndims

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_state_shape_real'

    ndims = SIZE(dims)

    IF (ndims < 2 .OR. ndims > 3) CALL finish(routine, 'Only 2d or 3d real states supported at this time')

    IF (ndims == 2) THEN
      IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nblks_c) THEN
        CALL finish(routine, 'Dimension mismatch')
      END IF
    ELSE IF (ndims == 3) THEN
      IF (this%domain%ntiles > 1) THEN
        IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nblks_c .OR. dims(3) /= this%domain%ntiles) THEN
          CALL finish(routine, 'Dimension mismatch')
        END IF
      ELSE IF (this%domain%ntiles == 1) THEN
        IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nblks_c .OR. dims(3) /= this%domain%ntiles) THEN
          CALL finish(routine, 'Dimension mismatch')
        END IF
      ELSE IF (this%domain%nlev > 1) THEN
        IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nlev .OR. dims(3) /= this%domain%nblks_c) THEN
          CALL finish(routine, 'Dimension mismatch')
        END IF
      ELSE
        CALL finish(routine, 'Dimension mismatch')
      END IF
    ELSE IF (ndims == 4) THEN
      IF (dims(1) /= this%domain%nproma .OR. dims(2) /= this%domain%nlev .OR. &
        & dims(3) /= this%domain%nblks_c .OR. dims(4) /= this%domain%ntiles) THEN
        CALL finish(routine, 'Dimension mismatch')
      END IF
    END IF

    CALL this%states%append(t_tmx_field(name, dims, type))
    tv => this%states%Search(name)
    IF (tv%bound) CALL finish(routine, 'State variable for '//name//' already bound to outside variable')
    CALL allocate_variable(tv)

    ! Add state variable to inputs list of process
    ! tv => this%states%Search(name)
    ! CALL this%inputs%list%append(tv)

    ! Add variable to new_state varlist with the same attributes as state
    CALL this%new_states%append(t_variable(name, dims, "", type_id="real"))
    tv => this%new_states%Search(name)
    IF (tv%bound) CALL finish(routine, 'new_state variable for '//name//' should not be bound to outside variable')
    CALL allocate_variable(tv)
!$OMP PARALLEL
    IF (ndims == 2) THEN
      CALL init(tv%r2d)
    ELSE IF (ndims == 3) THEN
      CALL init(tv%r3d)
    ELSE IF (ndims == 4) THEN
      CALL init(tv%r4d)
    ELSE
      CALL finish(routine, 'ERROR')
    END IF
!$OMP END PARALLEL
    
    ! Add variable to tendencies varlist with the same attributes as state
    CALL this%tendencies%append(t_variable(name, dims, "", type_id="real"))
    tv => this%tendencies%Search(name)
    IF (tv%bound) CALL finish(routine, 'Tendency variable for '//name//' should not be bound to outside variable')
    CALL allocate_variable(tv)
!$OMP PARALLEL
    IF (ndims == 2) THEN
      CALL init(tv%r2d)
    ELSE IF (ndims == 3) THEN
      CALL init(tv%r3d)
    ELSE IF (ndims == 4) THEN
      CALL init(tv%r4d)
    ELSE
      CALL finish(routine, 'ERROR')
    END IF
!$OMP END PARALLEL

    tv => this%states%search(name)
    CALL message(routine, 'New state: '//tv%name//' for process '//this%name)

  END SUBROUTINE Add_state_shape_real

  SUBROUTINE Add_state_r2d(this, name, type, field)

    USE mo_variable,        ONLY: t_variable, bind_variable

    CLASS(t_tmx_process), INTENT(inout), TARGET :: this
    CHARACTER(len=*),     INTENT(in)            :: name
    INTEGER,              INTENT(in)            :: type
    REAL(wp),             POINTER               :: field(:,:)

    TYPE(t_variable), POINTER :: tv

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_state_r2d'

    IF (SIZE(field,1) /= this%domain%nproma .OR. SIZE(field,2) /= this%domain%nblks_c) THEN
      CALL finish(routine, 'Dimension mismatch')
    END IF

    CALL this%states%append(t_tmx_field(name, SHAPE(field), type))
    ! CALL bind_tmx_field(this%states%search_field(name), field)
    CALL bind_variable(this%states%search(name), field)
    tv => this%states%Search(name)

    ! Add state variable to inputs list of process
    ! tv => this%states%Search(name)
    ! CALL this%inputs%list%append(tv)

    ! Add variable to new_states varlist with the same attributes as state
    CALL this%new_states%append(t_variable(name, SHAPE(field), "", type_id="real"))
    tv => this%new_states%Search(name)
    IF (tv%bound) CALL finish(routine, 'new_state variable for '//name//' should not be bound to outside variable')
    CALL allocate_variable(tv)
!$OMP PARALLEL
    CALL init(tv%r2d)
!$OMP END PARALLEL

    ! Add variable to tendencies varlist with the same attributes as state
    CALL this%tendencies%append(t_variable(name, SHAPE(field), "", type_id="real"))
    tv => this%tendencies%Search(name)
    IF (tv%bound) CALL finish(routine, 'Tendency variable for '//name//' should not be bound to outside variable')
    CALL allocate_variable(tv)
!$OMP PARALLEL
    CALL init(tv%r2d)
!$OMP END PARALLEL

    tv => this%states%search(name)
    CALL message(routine, 'New state: '//tv%name//' for process '//this%name)

  END SUBROUTINE Add_state_r2d

  SUBROUTINE Add_state_r3d(this, name, type, field)

    USE mo_variable,        ONLY: t_variable, bind_variable

    CLASS(t_tmx_process), INTENT(inout), TARGET :: this
    CHARACTER(len=*),     INTENT(in)            :: name
    INTEGER,              INTENT(in)            :: type
    REAL(wp),             POINTER               :: field(:,:,:)

    TYPE(t_variable), POINTER :: tv

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_state_r3d'

    IF (this%domain%ntiles > 1) THEN
      IF (SIZE(field,1) /= this%domain%nproma .OR. SIZE(field,2) /= this%domain%nblks_c) THEN
        CALL finish(routine, 'Dimension mismatch for '//name//' in '//this%name)
      END IF
    ELSE IF (this%domain%nlev > 1) THEN
      IF (SIZE(field,1) /= this%domain%nproma .OR. SIZE(field,3) /= this%domain%nblks_c) THEN
        CALL finish(routine, 'Dimension mismatch for '//name//' in '//this%name)
      END IF
    END IF

    CALL this%states%append(t_tmx_field(name, SHAPE(field), type))
    ! CALL bind_tmx_field(this%states%search_field(name), field)
    CALL bind_variable(this%states%search(name), field)
    ! tv => this%states%Search(name)

    ! Add state variable to inputs list of process
    ! tv => this%states%Search(name)
    ! CALL this%inputs%list%append(tv)

    ! Add variable to new_states varlist with the same attributes as state
    CALL this%new_states%append(t_variable(name, SHAPE(field), "", type_id="real"))
    tv => this%new_states%Search(name)
    IF (tv%bound) CALL finish(routine, 'new_states variable for '//name//' should not be bound to outside variable')
    CALL allocate_variable(tv)
!$OMP PARALLEL
    CALL init(tv%r3d)
!$OMP END PARALLEL

    ! Add variable to tendencies varlist with the same attributes as state
    CALL this%tendencies%append(t_variable(name, SHAPE(field), "", type_id="real"))
    tv => this%tendencies%Search(name)
    IF (tv%bound) CALL finish(routine, 'Tendency variable for '//name//' should not be bound to outside variable')
    CALL allocate_variable(tv)
!$OMP PARALLEL
    CALL init(tv%r3d)
!$OMP END PARALLEL

    tv => this%states%search(name)
    CALL message(routine, 'New state: '//tv%name//' for process '//this%name)

  END SUBROUTINE Add_state_r3d

  FUNCTION Get_tendency_r2d(this, name) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    CHARACTER(len=*),     INTENT(in) :: name
    REAL(wp), POINTER                :: result(:,:)

    TYPE(t_variable), POINTER :: tv

    tv => this%tendencies%Search(name)
    result => tv%r2d

  END FUNCTION Get_tendency_r2d

  FUNCTION Get_tendency_r3d(this, name) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    CHARACTER(len=*),     INTENT(in) :: name
    REAL(wp), POINTER                :: result(:,:,:)

    TYPE(t_variable), POINTER :: tv

    tv => this%tendencies%Search(name)
    result => tv%r3d

  END FUNCTION Get_tendency_r3d

  FUNCTION Get_diagnostic_r2d(this, name) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    CHARACTER(len=*),     INTENT(in) :: name
    REAL(wp), POINTER                :: result(:,:)

    TYPE(t_variable), POINTER :: tv

    tv => this%diagnostics%Search(name)
    result => tv%r2d

  END FUNCTION Get_diagnostic_r2d

  FUNCTION Get_diagnostic_r3d(this, name) RESULT(result)

    CLASS(t_tmx_process), INTENT(in) :: this
    CHARACTER(len=*),     INTENT(in) :: name
    REAL(wp), POINTER                :: result(:,:,:)

    TYPE(t_variable), POINTER :: tv

    tv => this%diagnostics%Search(name)
    result => tv%r3d

  END FUNCTION Get_diagnostic_r3d

  SUBROUTINE Set_time_scheme(this, time_scheme)

    CLASS(t_tmx_process), INTENT(inout) :: this
    CLASS(t_time_scheme), INTENT(in)    :: time_scheme

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_time_scheme'

    IF (ALLOCATED(this%time_scheme)) DEALLOCATE(this%time_scheme)

    ALLOCATE(this%time_scheme, source=time_scheme)

  END SUBROUTINE Set_time_scheme

  SUBROUTINE Step_forward(this)

    CLASS(t_tmx_process) :: this
    ! REAL(wp), INTENT(in) :: dt

    CHARACTER(len=*), PARAMETER :: routine = modname//':Step_forward'

    IF (ALLOCATED(this%time_scheme)) THEN
      CALL this%time_scheme%Step_forward(this, this%dt)
    ELSE
      CALL finish(routine, 'No time scheme defined')
    END IF

    END SUBROUTINE Step_forward

  END MODULE mo_tmx_process_class
