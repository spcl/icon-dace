! Collection of some tools for the upper-atmosphere configuration.
!
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

MODULE mo_upatmo_utils

  USE mo_kind,                      ONLY: wp
  USE mo_exception,                 ONLY: finish
  USE mo_impl_constants,            ONLY: SUCCESS
  USE mo_name_list_output_types,    ONLY: t_output_name_list
  USE mo_name_list_output_config,   ONLY: first_output_name_list, is_variable_in_output_nml
  USE mo_util_string,               ONLY: int2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_logical_1d
  PUBLIC :: is_variable_in_output_cond
  PUBLIC :: isInInterval
  PUBLIC :: t_varstate_set
  PUBLIC :: t_varstate

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_utils'

  !====================================================================================

  INTEGER, PARAMETER :: IERR_INIT   = SUCCESS + 100
  INTEGER, PARAMETER :: IERR_NSTATE = SUCCESS + 200
  INTEGER, PARAMETER :: IERR_LOCKED = SUCCESS + 300
  INTEGER, PARAMETER :: IERR_ISET   = SUCCESS + 400

  INTEGER, PARAMETER :: I4USE   = 1
  INTEGER, PARAMETER :: I4RESET = 2

  ! If you would add further variables to t_varstate_set, 
  ! please add them to src/upper_atmosphere/mo_upatmo_flowevent_utils as well. Thank you!
  TYPE t_varstate_set
    INTEGER :: i_old             = 1
    INTEGER :: i_new             = 2
    INTEGER :: n_swap            = 0
    INTEGER :: n_state           = 2
    INTEGER :: n_statep1         = 3
    LOGICAL :: l_swapped         = .FALSE.
    LOGICAL :: l_updated         = .FALSE.
    LOGICAL :: l_locked          = .FALSE.
    LOGICAL :: l_locking         = .TRUE.
    LOGICAL :: l_unlockable      = .TRUE.
    LOGICAL :: l_final           = .FALSE.
    LOGICAL :: l_finish_on_error = .TRUE.
    LOGICAL :: l_initialized     = .FALSE.
  END TYPE t_varstate_set

  TYPE t_varstate
    TYPE(t_varstate_set), PRIVATE :: set(2)
  CONTAINS
    PROCEDURE, PUBLIC  :: init        => t_varstate_init
    PROCEDURE, PUBLIC  :: swap        => t_varstate_swap
    PROCEDURE, PRIVATE :: lock        => t_varstate_lock
    PROCEDURE, PRIVATE :: unlock      => t_varstate_unlock
    PROCEDURE, PUBLIC  :: reset       => t_varstate_reset
    PROCEDURE, PUBLIC  :: clear       => t_varstate_clear
    PROCEDURE, PUBLIC  :: iold        => t_varstate_iold
    PROCEDURE, PUBLIC  :: inew        => t_varstate_inew
    PROCEDURE, PUBLIC  :: nstate      => t_varstate_nstate
    PROCEDURE, PUBLIC  :: nswap       => t_varstate_nswap
    PROCEDURE, PRIVATE :: lswapped    => t_varstate_lswapped
    PROCEDURE, PUBLIC  :: lupdated    => t_varstate_lupdated
    PROCEDURE, PUBLIC  :: llocked     => t_varstate_llocked
    PROCEDURE, PUBLIC  :: lfinal      => t_varstate_lfinal
    PROCEDURE, PUBLIC  :: lunlockable => t_varstate_lunlockable
    PROCEDURE, PUBLIC  :: lswappable  => t_varstate_lswappable
    PROCEDURE, PRIVATE :: iget        => t_varstate_iget
    PROCEDURE, PRIVATE :: lget        => t_varstate_lget
    PROCEDURE, PUBLIC  :: getSet      => t_varstate_getSet
  END TYPE t_varstate

  !====================================================================================

  INTERFACE isInInterval
    MODULE PROCEDURE isInInterval_integer
    MODULE PROCEDURE isInInterval_real
  END INTERFACE isInInterval

CONTAINS !..................................................................................

  !>
  !! Initialize logical 1d-array.
  !! (Introduced, because 'externals/fortran-support/src/mo_fortran_tools: init_contiguous_l' 
  !! does not suit our purposes.)
  !!
  SUBROUTINE init_logical_1d( variable,   & !inout
    &                         value,      & !in
    &                         opt_ilist,  & !optin
    &                         opt_mask    ) !optin
    ! In/out variables
    LOGICAL,                    INTENT(INOUT) :: variable(:)  ! Logical array to be assigned with 'value'
    LOGICAL,                    INTENT(IN)    :: value         
    INTEGER,          OPTIONAL, INTENT(IN)    :: opt_ilist(:) ! Optional list with indices of 'variable' 
                                                              ! that shall or shall not be assigned with 'value'. 
                                                              ! The indices are assumed to be 
                                                              ! in '[1, size(variable)]', 
                                                              ! if present, or in '[1, SIZE(variable)]' otherwise.
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: opt_mask     ! "list" -> those indices of 'variable' stored 
                                                              ! in 'opt_ilist' are not assigned with 'value'
                                                              ! "complement" -> those indices of 'variable' 
                                                              ! not stored in 'opt_ilist' are not assigned with 'value'. 
                                                              ! The case that 'opt_ilist' is present, 
                                                              ! while 'opt_mask' is absent, is interpreted 
                                                              ! as 'opt_mask = "list"'.

    ! Local variables
    LOGICAL, ALLOCATABLE :: mask(:)
    LOGICAL :: lmask
    INTEGER :: varsize, jloop, istat
    CHARACTER(LEN=*), PARAMETER :: mask_list       = "list"
    CHARACTER(LEN=*), PARAMETER :: mask_complement = "complement"
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':init_logical_1d'

    !---------------------------------------------------------

    varsize = SIZE(variable)

    IF (PRESENT(opt_ilist)) THEN
      lmask = .TRUE.
      IF ( MINVAL(opt_ilist) < 1 .OR. &
        &  MAXVAL(opt_ilist) > varsize        ) THEN
        CALL finish(routine, "Index in opt_ilist outside index range of variable.")
      ENDIF
      ALLOCATE(mask(varsize), STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, "Allocation of mask failed.")

      IF (PRESENT(opt_mask)) THEN
        SELECT CASE(TRIM(opt_mask))
        CASE(mask_list)
          mask(:) = .FALSE.
        CASE(mask_complement)
          mask(:) = .TRUE.
        CASE default
          CALL finish(routine, "Invalid opt_mask.")
        END SELECT
      ELSE
        mask(:) = .FALSE.
      ENDIF

      DO jloop = 1, SIZE(opt_ilist)
        mask(opt_ilist(jloop)) = .NOT. mask(opt_ilist(jloop))
      ENDDO
    ELSE
      lmask = .FALSE.
    ENDIF
    
    IF (lmask) THEN
      DO jloop = 1, varsize
        IF (.NOT. mask(jloop)) variable(jloop) = value
      ENDDO
    ELSE
      DO jloop = 1, varsize
        variable(jloop) = value
      ENDDO
    ENDIF

    IF (lmask) THEN
      DEALLOCATE(mask, STAT=istat)
      IF (istat /= SUCCESS) CALL finish(routine, "Deallocation of mask failed.")
    ENDIF

  END SUBROUTINE init_logical_1d

  !====================================================================================

  !>
  !! Copy of 'src/configure_model/mo_name_list_output_config: is_variable_in_output', 
  !! where (optional) conditions have to be met.
  !!
  FUNCTION is_variable_in_output_cond( var_name,               &
    &                                  opt_dom,                &
    &                                  opt_filetype            ) RESULT(retval)

    ! In/out variables
    LOGICAL                                        :: retval
    CHARACTER(LEN=*),                   INTENT(IN) :: var_name               ! Variable name
    ! Optional conditions:
    INTEGER,                  OPTIONAL, INTENT(IN) :: opt_dom(:)             ! Domain
    INTEGER,                  OPTIONAL, INTENT(IN) :: opt_filetype(:)        ! File type (GRIB1/2, NetCDF2/4)

    ! Local variables
    TYPE (t_output_name_list), POINTER :: p_onl
    LOGICAL :: l_dom, l_filetype
    LOGICAL :: l_anycond, l_met

    !---------------------------------------------------------

    retval = .FALSE.

    ! Check, if variable name is non-empty
    IF (LEN_TRIM(var_name) == 0) RETURN

    ! Check for optional conditions
    !
    l_dom      = PRESENT(opt_dom)
    l_filetype = PRESENT(opt_filetype)
    l_anycond  = l_dom .OR. l_filetype

    p_onl  => first_output_name_list

    ! If there is no optional condition, 
    ! this function becomes 'is_variable_in_output' effectively
    IF (.NOT. l_anycond) THEN
      DO WHILE (ASSOCIATED(p_onl) .AND. .NOT. retval)
        retval = is_variable_in_output_nml(p_onl, var_name=var_name)
        p_onl  => p_onl%next
      END DO
    ELSE
      DO WHILE (ASSOCIATED(p_onl) .AND. .NOT. retval)
        ! Check, if current list element satisfies ALL optional conditions, ...
        l_met = .TRUE.
        IF (l_dom)      l_met = l_met .AND. (ANY(opt_dom(:) == p_onl%dom))
        IF (l_filetype) l_met = l_met .AND. (ANY(opt_filetype(:) == p_onl%filetype))
        ! ... and apply check for variable only in this case
        IF(l_met) retval = is_variable_in_output_nml(p_onl, var_name=var_name)
        p_onl => p_onl%next
      END DO
    ENDIF  !IF (.NOT. l_anycond) 

  END FUNCTION is_variable_in_output_cond

  !====================================================================================

  FUNCTION isInInterval_integer( number,    &
    &                            opt_clbnd, &
    &                            opt_olbnd, &
    &                            opt_cubnd, &
    &                            opt_oubnd  )

    ! In/out variables
    LOGICAL                       :: isInInterval_integer
    INTEGER,           INTENT(IN) :: number
    INTEGER, OPTIONAL, INTENT(IN) :: opt_clbnd
    INTEGER, OPTIONAL, INTENT(IN) :: opt_olbnd
    INTEGER, OPTIONAL, INTENT(IN) :: opt_cubnd
    INTEGER, OPTIONAL, INTENT(IN) :: opt_oubnd

    ! Local variables
    LOGICAL  :: l_present_clbnd, l_present_cubnd
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':isInInterval_integer'
    !---------------------------------------------------------

    ! Covered (mutually exclusive) cases:
    ! (Open bounds might be regarded as questionable 
    ! for the integer version of this function.  
    ! However, for reasons of simplicity and coherence 
    ! with the real-number version of this function, 
    ! we keep it this way)
    ! * number in [opt_clbnd, opt_cubnd]
    ! * number in (opt_olbnd, opt_cubnd]
    ! * number in [opt_clbnd, opt_oubnd)
    ! * number in (opt_olbnd, opt_oubnd)

    ! Closed lower bound present?
    l_present_clbnd = PRESENT(opt_clbnd)
    ! Closed upper bound present?
    l_present_cubnd = PRESENT(opt_cubnd)

    IF (l_present_clbnd .EQV. PRESENT(opt_olbnd)) THEN
      ! Either both a closed lower boundary and an open lower boundary 
      ! are present or none is present
      CALL finish(routine, "Invalid argument for lower bound.")
    ELSEIF (l_present_cubnd .EQV. PRESENT(opt_oubnd)) THEN
      ! Either both a closed upper boundary and an open upper boundary 
      ! are present or none is present
      CALL finish(routine, "Invalid argument for upper bound.")
    ENDIF
    
    ! Is number equal to or greater than lower bound?
    IF (l_present_clbnd) THEN
      isInInterval_integer = number >= opt_clbnd
    ELSE
      isInInterval_integer = number > opt_olbnd
    ENDIF

    ! Is number equal to less than upper bound?
    IF (l_present_cubnd) THEN
      isInInterval_integer = isInInterval_integer .AND. number <= opt_cubnd
    ELSE
      isInInterval_integer = isInInterval_integer .AND. number < opt_oubnd
    ENDIF    

  END FUNCTION isInInterval_integer

  !====================================================================================

  FUNCTION isInInterval_real( number,    &
    &                         opt_clbnd, &
    &                         opt_olbnd, &
    &                         opt_cubnd, &
    &                         opt_oubnd  )

    ! In/out variables
    LOGICAL                        :: isInInterval_real
    REAL(wp),           INTENT(IN) :: number
    REAL(wp), OPTIONAL, INTENT(IN) :: opt_clbnd
    REAL(wp), OPTIONAL, INTENT(IN) :: opt_olbnd
    REAL(wp), OPTIONAL, INTENT(IN) :: opt_cubnd
    REAL(wp), OPTIONAL, INTENT(IN) :: opt_oubnd

    ! Local variables
    LOGICAL  :: l_present_clbnd, l_present_cubnd
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':isInInterval_real'
    !---------------------------------------------------------

    ! Covered (mutually exclusive) cases:
    ! * number in [opt_clbnd, opt_cubnd]
    ! * number in (opt_olbnd, opt_cubnd]
    ! * number in [opt_clbnd, opt_oubnd)
    ! * number in (opt_olbnd, opt_oubnd)

    ! Closed lower bound present?
    l_present_clbnd = PRESENT(opt_clbnd)
    ! Closed upper bound present?
    l_present_cubnd = PRESENT(opt_cubnd)

    IF (l_present_clbnd .EQV. PRESENT(opt_olbnd)) THEN
      ! Either both a closed lower boundary and an open lower boundary 
      ! are present or none is present
      CALL finish(routine, "Invalid argument for lower bound.")
    ELSEIF (l_present_cubnd .EQV. PRESENT(opt_oubnd)) THEN
      ! Either both a closed upper boundary and an open upper boundary 
      ! are present or none is present
      CALL finish(routine, "Invalid argument for upper bound.")
    ENDIF
    
    ! Is number equal to or greater than lower bound?
    IF (l_present_clbnd) THEN
      isInInterval_real = .NOT. (number < opt_clbnd)
    ELSE
      isInInterval_real = number > opt_olbnd
    ENDIF

    ! Is number equal to or less than upper bound?
    IF (l_present_cubnd) THEN
      isInInterval_real = isInInterval_real .AND. .NOT. (number > opt_cubnd)
    ELSE
      isInInterval_real = isInInterval_real .AND. number < opt_oubnd
    ENDIF    

  END FUNCTION isInInterval_real

  !====================================================================================

  SUBROUTINE t_varstate_init( varstate,                 & !class
    &                         nstate,                   & !in
    &                         optDefaultFinishOnError,  & !optin
    &                         optDefaultLocking,        & !optin
    &                         optError                  ) !optout
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    INTEGER,                   INTENT(IN)    :: nstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optDefaultFinishOnError
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optDefaultLocking
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError

    TYPE(t_varstate_set), POINTER :: set4use, set4reset
    INTEGER :: error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_init'

    !----------------------------------------------

    error = SUCCESS
    set4use   => varstate%set(I4USE)
    set4reset => varstate%set(I4RESET)
    IF (PRESENT(optDefaultFinishOnError)) THEN
      set4use%l_finish_on_error = optDefaultFinishOnError
    ELSE
      set4use%l_finish_on_error = .TRUE.
    ENDIF
    IF (PRESENT(optDefaultLocking)) THEN
      set4use%l_locking = optDefaultLocking
    ELSE
      set4use%l_locking = .TRUE.
    ENDIF
    IF (.NOT. set4use%l_initialized) THEN
      IF (nstate == 1 .OR. nstate == 2) THEN
        set4use%n_state       = nstate
        set4use%n_statep1     = nstate + 1
        set4use%i_old         = 1
        set4use%i_new         = set4use%n_statep1 - set4use%i_old
        set4use%n_swap        = 0
        set4use%l_swapped     = .FALSE.
        set4use%l_updated     = .FALSE.
        set4use%l_locked      = .FALSE.
        set4use%l_unlockable  = .TRUE.
        set4use%l_final       = .FALSE.
        set4use%l_initialized = .TRUE.
      ELSE
        error = IERR_NSTATE
      ENDIF
    ELSE
      error = IERR_INIT
    ENDIF
    set4reset = set4use
    set4reset => NULL()
    IF (PRESENT(optError)) optError = error
    IF (set4use%l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF
    set4use   => NULL()

  END SUBROUTINE t_varstate_init

  !***************************************************************

  ! Conceptual copy of 'externals/fortran-support/src/mo_fortran_tools: swap_init':

  SUBROUTINE t_varstate_swap( varstate,         & !class
    &                         optFinal,         & !optin
    &                         optCountAsUpdate, & !optin
    &                         optFinishOnError, & !optin
    &                         optError          ) !optout
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinal
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optCountAsUpdate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error
    LOGICAL :: l_finish_on_error, l_final, l_count_as_update
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_swap'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (PRESENT(optFinal)) THEN
      l_final = optFinal
    ELSE
      l_final = .FALSE.
    ENDIF
    IF (PRESENT(optCountAsUpdate)) THEN
      l_count_as_update = optCountAsUpdate
    ELSE
      l_count_as_update = .TRUE.
    ENDIF
    IF (set4use%l_initialized) THEN
      IF (.NOT. (set4use%l_locking .AND. set4use%l_locked)) THEN
        set4use%i_old     = set4use%n_statep1 - set4use%i_old
        set4use%i_new     = set4use%n_statep1 - set4use%i_new
        set4use%n_swap    = set4use%n_swap + 1
        set4use%l_swapped = .TRUE.
        set4use%l_updated = l_count_as_update
        IF (set4use%l_locking) THEN
          CALL varstate%lock(optFinishOnError, optError)
          set4use%l_final = l_final
        ENDIF
      ELSEIF (set4use%l_locking .AND. set4use%l_locked) THEN
        error = IERR_LOCKED
      ENDIF
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END SUBROUTINE t_varstate_swap

  !***************************************************************

  SUBROUTINE t_varstate_clear( varstate,         & !class
    &                          optFinishOnError, & !optin
    &                          optError          ) !optout
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_clear'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (set4use%l_initialized) THEN
      set4use%l_swapped = .FALSE.
      set4use%l_updated = .FALSE.
      IF (set4use%l_locking .AND. set4use%l_unlockable .AND. .NOT. set4use%l_final) THEN
        CALL varstate%unlock(optFinishOnError, optError)
      ELSEIF (set4use%l_locking .AND. set4use%l_final) THEN
        CALL varstate%lock(optFinishOnError, optError)
        set4use%l_unlockable = .FALSE.
        set4use%l_final      = .FALSE.
      ENDIF
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END SUBROUTINE t_varstate_clear

  !***************************************************************

  SUBROUTINE t_varstate_lock( varstate,         & !class
    &                         optFinishOnError, & !optin
    &                         optError          ) !optout
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_lock'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (set4use%l_initialized) THEN
      set4use%l_locked = .TRUE.
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END SUBROUTINE t_varstate_lock

  !***************************************************************

  SUBROUTINE t_varstate_unlock( varstate,         & !class
    &                           optFinishOnError, & !optin
    &                           optError          ) !optout
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_unlock'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (set4use%l_initialized) THEN
      IF (.NOT. (set4use%l_locking .AND. .NOT. set4use%l_unlockable)) THEN
        set4use%l_locked = .FALSE.
      ELSE
        error = IERR_LOCKED
      ENDIF
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END SUBROUTINE t_varstate_unlock

  !***************************************************************

  SUBROUTINE t_varstate_reset( varstate,          & !class
    &                          optSet4Reset,      & !optin
    &                          optFinishOnError,  & !optin
    &                          optError           ) !optout
    
    CLASS(t_varstate),              TARGET, INTENT(INOUT) :: varstate
    TYPE(t_varstate_set), OPTIONAL, TARGET, INTENT(IN)    :: optSet4Reset
    LOGICAL,              OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER,              OPTIONAL,         INTENT(OUT)   :: optError

    TYPE(t_varstate_set), POINTER :: set4use, set4reset
    INTEGER :: error
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_reset'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    IF (PRESENT(optSet4Reset)) THEN
      set4reset => optSet4Reset
    ELSE
      set4reset => varstate%set(I4RESET)
    ENDIF
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (set4use%l_initialized) THEN
      set4use = set4reset
    ELSE
      error = IERR_INIT
    ENDIF
    set4use   => NULL()
    set4reset => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END SUBROUTINE t_varstate_reset

  !***************************************************************

  FUNCTION t_varstate_iold( varstate,           & !class
    &                       optFinishOnError,   & !optin
    &                       optError          ) & !optout
    &                       RESULT(iold)
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    INTEGER                                  :: iold

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    iold = varstate%iget(set4use%i_old, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_iold

  !***************************************************************

  FUNCTION t_varstate_inew( varstate,           & !class
    &                       optFinishOnError,   & !optin
    &                       optError          ) & !optout
    &                       RESULT(inew)
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    INTEGER                                  :: inew

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    inew = varstate%iget(set4use%i_new, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_inew

  !***************************************************************

  FUNCTION t_varstate_nswap( varstate,           & !class
    &                        optFinishOnError,   & !optin
    &                        optError          ) & !optout
    &                        RESULT(nswap)
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    INTEGER                                  :: nswap

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    nswap = varstate%iget(set4use%n_swap, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_nswap

  !***************************************************************

  FUNCTION t_varstate_nstate( varstate,           & !class
    &                         optFinishOnError,   & !optin
    &                         optError          ) & !optout
    &                         RESULT(nstate)
    
    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    INTEGER                                  :: nstate

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    nstate = varstate%iget(set4use%n_state, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_nstate

  !***************************************************************

  FUNCTION t_varstate_iget( varstate,           & !class
    &                       iwant,              & !in
    &                       optFinishOnError,   & !optin
    &                       optError          ) & !optout
    &                       RESULT(iget)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    INTEGER,                   INTENT(IN)    :: iwant
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    INTEGER                          :: iget

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_iget'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    iget = -999
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (set4use%l_initialized) THEN
      iget = iwant
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END FUNCTION t_varstate_iget

  !***************************************************************

  FUNCTION t_varstate_lswapped( varstate,           & !class
    &                           optFinishOnError,   & !optin
    &                           optError          ) & !optout
    &                           RESULT(lswapped)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: lswapped

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    lswapped = varstate%lget(set4use%l_swapped, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_lswapped

  !***************************************************************

  FUNCTION t_varstate_lupdated( varstate,           & !class
    &                           optFinishOnError,   & !optin
    &                           optError          ) & !optout
    &                           RESULT(lupdated)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: lupdated

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    lupdated = varstate%lget(set4use%l_updated, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_lupdated

  !***************************************************************

  FUNCTION t_varstate_llocked( varstate,           & !class
    &                          optFinishOnError,   & !optin
    &                          optError          ) & !optout
    &                          RESULT(llocked)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: llocked

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    llocked = varstate%lget(set4use%l_locked, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_llocked

  !***************************************************************

  FUNCTION t_varstate_lfinal( varstate,           & !class
    &                         optFinishOnError,   & !optin
    &                         optError          ) & !optout
    &                         RESULT(lfinal)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: lfinal

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    lfinal = varstate%lget(set4use%l_final, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_lfinal

  !***************************************************************

  FUNCTION t_varstate_lunlockable( varstate,           & !class
    &                              optFinishOnError,   & !optin
    &                              optError          ) & !optout
    &                              RESULT(lunlockable)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: lunlockable

    TYPE(t_varstate_set), POINTER :: set4use

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    lunlockable = varstate%lget(set4use%l_unlockable, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_lunlockable

  !***************************************************************

  FUNCTION t_varstate_lswappable( varstate,           & !class
    &                             optFinishOnError,   & !optin
    &                             optError          ) & !optout
    &                             RESULT(lswappable)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: lswappable

    TYPE(t_varstate_set), POINTER :: set4use
    LOGICAL :: l_swappable

    !----------------------------------------------

    set4use => varstate%set(I4USE)
    IF (set4use%l_locking) THEN
      l_swappable = .NOT. (set4use%l_locked .OR. set4use%l_final)
    ELSE
      l_swappable = .TRUE.
    ENDIF
    lswappable = varstate%lget(l_swappable, optFinishOnError, optError)
    set4use => NULL()

  END FUNCTION t_varstate_lswappable

  !***************************************************************

  FUNCTION t_varstate_lget( varstate,           & !class
    &                       lwant,              & !in
    &                       optFinishOnError,   & !optin
    &                       optError          ) & !optout
    &                       RESULT(lget)

    CLASS(t_varstate), TARGET, INTENT(INOUT) :: varstate
    LOGICAL,                   INTENT(IN)    :: lwant
    LOGICAL, OPTIONAL,         INTENT(IN)    :: optFinishOnError
    INTEGER, OPTIONAL,         INTENT(OUT)   :: optError
    LOGICAL                                  :: lget

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_lget'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    lget = .FALSE.
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (set4use%l_initialized) THEN
      lget = lwant
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END FUNCTION t_varstate_lget

  !***************************************************************

  FUNCTION t_varstate_getSet( varstate,           & !class
    &                         optWhichSet,        & !optin
    &                         optFinishOnError,   & !optin
    &                         optError          ) &!optout
    &                         RESULT(set)
    
    CLASS(t_varstate),  TARGET, INTENT(INOUT) :: varstate
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: optWhichSet
    LOGICAL,          OPTIONAL, INTENT(IN)    :: optFinishOnError
    INTEGER,          OPTIONAL, INTENT(OUT)   :: optError
    TYPE(t_varstate_set)                      :: set

    TYPE(t_varstate_set), POINTER :: set4use
    INTEGER :: error, iset
    LOGICAL :: l_finish_on_error
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':t_varstate_getSet'

    !----------------------------------------------

    error = SUCCESS
    set4use => varstate%set(I4USE)
    IF (PRESENT(optFinishOnError)) THEN
      l_finish_on_error = optFinishOnError
    ELSE
      l_finish_on_error = set4use%l_finish_on_error
    ENDIF
    IF (PRESENT(optWhichSet)) THEN
      SELECT CASE (optWhichSet)
      CASE("set4use")
        iset = I4USE
      CASE("set4reset")
        iset = I4RESET
      CASE DEFAULT
        error = IERR_ISET
      END SELECT
    ELSE
      iset = I4USE
    ENDIF
    IF (set4use%l_initialized) THEN
      set = varstate%set(iset)
    ELSE
      error = IERR_INIT
    ENDIF
    set4use => NULL()
    IF (PRESENT(optError)) optError = error
    IF (l_finish_on_error .AND. (error /= SUCCESS)) THEN
      CALL finish (routine, 'Error code: '//TRIM(int2string(error)))
    ENDIF

  END FUNCTION t_varstate_getSet

END MODULE mo_upatmo_utils
