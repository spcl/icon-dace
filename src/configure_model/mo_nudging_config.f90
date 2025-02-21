! Processing of nudging parameters for the nudging types:
! - Upper boundary nudging
! - Global nudging
! For the lateral boundary nudging, please see:
! - src/configure_model/mo_limarea_config
! - src/configure_model/mo_interpol_config
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

MODULE mo_nudging_config

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message, message_text
  USE mo_impl_constants,        ONLY: MAX_CHAR_LENGTH, max_dom
  USE mo_vertical_coord_table,  ONLY: vct_a
  USE mo_util_string,           ONLY: int2string, real2string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: indg_type, indg_profile, ithermdyn_type
  PUBLIC :: NDG_VAR_STR_LEN, NDG_VARLIST_STR_LEN
  PUBLIC :: indg_var
  PUBLIC :: cndg_var
  PUBLIC :: t_nudging_config
  PUBLIC :: nudging_config
  PUBLIC :: configure_nudging
  PUBLIC :: l_global_nudging

  ! Module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nudging_config'

  !---------------------------------------------------------
  !                       Parameters
  !---------------------------------------------------------

  INTEGER, PARAMETER :: NDG_VAR_STR_LEN = 10  ! For variable identifiers
  ! - NDG_VARLIST_STR_LEN -> see a bit further below

  !---------------------------------------------------------
  !                     Parameter types
  !---------------------------------------------------------

  ! Note: the following parameter types are introduced for convenience: 
  ! * combine thematically related identifiers, to unburden the USE-areas in modules, where they are required
  ! * use the number of list box entries 'ixyz%nitem' for allocation purposes and loop boundaries
  ! * most changes with regard to the identifiers can be confined to this module
  ! * avoid optical overload of 'src/shared/mo_impl_constants'

  ! Their construction generally follows the rule:
  !
  ! TYPE t_ixyz
  !   INTEGER :: a      ! Identifier 1
  !   INTEGER :: b      ! Identifier 2
  !   INTEGER :: c      ! Identifier 3
  !   INTEGER :: d      ! Identifier 4
  !   ! ...
  !   !
  !   INTEGER :: nitem  ! Number of identifiers in t_ixyz
  ! END TYPE t_ixyz
  ! TYPE(t_ixyz), PARAMETER :: ixyz = t_ixyz( 1, & ! a
  !                                           2, & ! b = a + 1
  !                                           3, & ! c = b + 1
  !                                           4, & ! d = c + 1
  !
  !                                           4  ) ! nitem = d, since d is the last list element
  !
  ! Why this construction rule?
  ! -> We may want to loop over whatever the items of ixyz identify: 
  !
  !    DO jitem = 1, ixyz%nitem
  !      ! jitem successivly takes the values of ixyz%a, ixyz%b, ixyz%c, ... 
  !      ! and can be compared with namelist entries etc.
  !    ENDDO
  ! 
  ! -> The items of ixyz shall be parameters, so a setup during runtime 
  !    by means of some object is ineligible. 
  !    We have to count the items of ixyz manually and store the value in nitem. 
  !
  ! Please, be careful, if you modify the following types, 
  ! since a reasonable check for the adherence to the above construction rule 
  ! during runtime is not possible.

  ! IDENTIFIERS for nudging type:.....................................................................................
  !
  TYPE t_indg_type
    INTEGER :: off      ! Switched off
    INTEGER :: ubn      ! Upper boundary nudging 
    INTEGER :: globn    ! Global nudging
    !   
    INTEGER :: nitem    ! Number of preceding list elements
  END TYPE t_indg_type
  TYPE(t_indg_type), PARAMETER :: indg_type = t_indg_type( 0, & ! off
    &                                                      1, & ! ubn
    &                                                      2, & ! globn
    !
    &                                                      2  ) ! (not 3) nitem 

  ! IDENTIFIERS for nudging profile:..................................................................................
  !
  ! For a visualization of the global nudging profiles  
  ! please see the end of this file.
  !
  TYPE t_indg_profile
    INTEGER :: sqrddist  ! Squared scaled vertical distance
    INTEGER :: const     ! Constant profile for nudging strength (for global nudging)
    INTEGER :: tanh      ! Hyperbolic tangent profile for nudging strength (for global nudging)
    INTEGER :: trapid    ! Trapezoidal profile (for global nudging)
    ! 
    INTEGER :: nitem     ! Number of list elements
  END TYPE t_indg_profile
  TYPE(t_indg_profile), PARAMETER :: indg_profile = t_indg_profile( 1, & ! sqrddist
    &                                                               2, & ! const
    &                                                               3, & ! tanh
    &                                                               4, & ! trapid
    ! 
    &                                                               4  ) ! nitem 

  ! IDENTIFIERS for type of thermodynamic variables:..................................................................
  !  
  TYPE t_ithermdyn_type
    INTEGER :: hydrostatic      ! Hydrostatically balanced
    INTEGER :: nonhydrostatic   ! Non-hydrostatic
    ! 
    INTEGER :: nitem            ! Number of list elements
  END TYPE t_ithermdyn_type
  TYPE(t_ithermdyn_type), PARAMETER ::  ithermdyn_type = t_ithermdyn_type( 1, & ! hydrostatic
    &                                                                      2, & ! nonhydrostatic
    !
    &                                                                      2  ) ! nitem 

  ! IDENTIFIERS for variables/variable groups subject to nudging:.....................................................
  !  
  TYPE t_indg_var
    INTEGER :: vn         ! Horizontal wind
    INTEGER :: thermdyn   ! Thermodynamic variables
    INTEGER :: qv         ! Water vapour
    ! 
    INTEGER :: nitem      ! Number of preceding list elements
    ! Abstracts
    INTEGER :: all        ! All variables
    !
    INTEGER :: nitem_2    ! Number of preceding list elements, except for nitem
  END TYPE t_indg_var
  TYPE(t_indg_var), PARAMETER :: indg_var = t_indg_var( 1, & ! vn
    &                                                   2, & ! thermdyn
    &                                                   3, & ! qv
    !
    &                                                   3, & ! nitem 
    ! Abstracts
    &                                                   4, & ! all
    !
    &                                                   4  ) ! nitem_2

  ! Max string length for variable lists
  INTEGER, PARAMETER :: NDG_VARLIST_STR_LEN = 2 * indg_var%nitem * NDG_VAR_STR_LEN  ! Factor 2 is margin for commas etc.

  ! CHARACTER IDENTIFIERS for variables/variable groups subject to nudging:...........................................
  !  
  CHARACTER(LEN=NDG_VAR_STR_LEN), PARAMETER :: cndg_var(indg_var%nitem_2) = &
    & (/ "vn        ", & ! vn
    &    "thermdyn  ", & ! thermdyn
    &    "qv        ", & ! qv
    &    "all       " /) ! all

  !--------------------------------------------------------
  !         Derived types for nudging parameters
  !--------------------------------------------------------

  ! Thresholds
  TYPE t_thr
    INTEGER :: low
    INTEGER :: med
    integer :: high
  END TYPE t_thr

  ! MAIN CONFIGURATION TYPE:..........................................................................................
  !
  TYPE t_nudging_config

    ! Namelist parameters
    INTEGER  :: nudge_type                           ! Nudging type identifier:
                                                     ! - 0 -> switched off
                                                     ! - 1 -> upper boundary nudging
                                                     ! - 2 -> global nudging 
    REAL(wp) :: max_nudge_coeff_vn                   ! Max. nudging coefficient for horizontal wind
    REAL(wp) :: max_nudge_coeff_thermdyn             ! Max. nudging coefficient for thermodynamic variables
    REAL(wp) :: nudge_start_height                   ! Nudging start height

    ! Non-namelist parameters 
    ! in case of upper boundary nudging, 
    ! namelist parameters 
    ! in case of global nudging
    INTEGER  :: thermdyn_type                        ! Type of thermodynamic variables (from the perspective
                                                     ! of the nudging-target variables in ICON):
                                                     ! - 1 -> hydrostatic: hydrostatic pressure & temperature are nudged
                                                     !        (diagnostic variables)
                                                     ! - 2 -> nonhydrostatic: density & virtual potential temperature are nudged
                                                     !        (prognostic variables)
    INTEGER  :: nudge_profile                        ! Profile of nudging strength:
                                                     ! - 1 -> inverse squared scaled vertical distance from nudging start height
                                                     ! - 2 -> constant profile (for global nudging only)
                                                     ! - 3 -> hyperbolic tangent profile (for global nudging only)
                                                     ! - 4 -> trapezoidal profile (for global nudging only)
    REAL(wp) :: nudge_end_height                     ! Nudging end height
    REAL(wp) :: nudge_scale_height                   ! Scale height for profiles
    CHARACTER(LEN=NDG_VARLIST_STR_LEN) :: nudge_var  ! Nudging variable(s) (for global nudging only): 
                                                     ! - 'vn'       -> horizontal wind 
                                                     ! - 'thermdyn' -> thermodynamic variables
                                                     ! - 'qv'       -> water vapour
                                                     ! - 'all'      -> all variables (i.e., 'vn' & 'thermdyn' & 'qv')
                                                     ! - string with a comma-separated list of variables, e.g., 'vn,thermdyn'
    REAL(wp) :: max_nudge_coeff_qv                   ! Max. nudging coefficient for water vapour (for global nudging only)

    ! Derived parameters
    LOGICAL  :: lnudging = .FALSE.                   ! Overall switch for nudging
    LOGICAL  :: ltype(indg_type%nitem)               ! And nudging-type-specific
    INTEGER  :: ilev_start                           ! Index of grid layer corresponding to 'nudge_end_height'
    INTEGER  :: ilev_end                             ! Index of grid layer corresponding to 'nudge_start_height'
    LOGICAL  :: lvar(indg_var%nitem)                 ! .TRUE. -> variable is subject to nudging
    
    ! Auxiliary variables
    TYPE(t_thr)      :: msg_thr                      ! Message thresholds for print
    TYPE(t_thr)      :: tmr_thr                      ! Timer thresholds
    LOGICAL          :: ltimer = .FALSE.             ! .TRUE.  -> global nudging timer is switched on
    LOGICAL          :: lmessage = .FALSE.           ! .TRUE.  -> message output during integration is switched on
    LOGICAL          :: llatbc_type_const = .FALSE.  ! .TRUE.  -> constant-in-time lateral boundary condition data
    LOGICAL          :: lqv_nudgable = .TRUE.        ! .TRUE.  -> model setup allows nudging of water vapour

    ! Parameters for type status
    LOGICAL  :: lchecked    = .FALSE.                ! .TRUE. -> crosschecks of namelist took place
    LOGICAL  :: lconfigured = .FALSE.                ! .TRUE. -> configuration of this type took place
 
  END TYPE t_nudging_config

  ! Important: 'nudging_config' is uninitialized until the call of 'configure_nudging'. 
  ! So whenever you use 'nudging_config' outside this module at a relatively early point of the program sequence, 
  ! please check its status by: 
  ! 
  !  IF (nudging_config%lconfigured) THEN 
  !    ! Do what you would like to do
  !  ELSE 
  !    CALL finish(routine, "nudging_config still unconfigured.")
  !  ENDIF
  ! 
  ! to be on the safe side.
  TYPE(t_nudging_config), TARGET :: nudging_config(max_dom)

  ! Convenience substitute for 'nudging_config%ltype(indg_type%globn)'
  ! for use in 'src/atm_dyn_iconam/mo_nh_stepping' 
  ! (and 'src/drivers/mo_atmo_nonhydrostatic: destruct_atmo_nonhydrostatic').
  LOGICAL :: l_global_nudging = .FALSE.

CONTAINS !............................................................................................................

  !>
  !! Subroutine to configure:
  !! - Upper boundary nudging
  !! - Global nudging
  !!
  SUBROUTINE configure_nudging( nlev, nshift_total, n_dom, msg_level, timers_level )

    ! In/out variables
    INTEGER, INTENT(IN) :: nlev            ! Number of grid layers of primary domain
    INTEGER, INTENT(IN) :: nshift_total(:) ! total shift of model top with respect to base domain
    INTEGER, INTENT(IN) :: n_dom           ! total number of domains
    INTEGER, INTENT(IN) :: msg_level       ! Control parameter for message output
    INTEGER, INTENT(IN) :: timers_level    ! Control parameter for timer

    ! Local variables
    REAL(wp) :: height, start_height, end_height
    INTEGER  :: jg
    INTEGER  :: jk, nlevp1, jvar
    LOGICAL  :: lfound, lremove_qv
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: c4print
    CHARACTER(LEN=MAX_CHAR_LENGTH), PARAMETER :: &
      & routine = modname//":configure_nudging"

    !----------------------------------------

    ! Crosscheck of namelist entries should have taken place beforehand 
    ! in 'src/namelists/mo_nudging_nml: check_nudging'.
    IF (.NOT. ALL(nudging_config(1:n_dom)%lchecked)) THEN
      CALL finish(routine, "Crosscheck of nudging_nml still pending. " &
        & //"Please, check the program sequence.")
    ENDIF

    !---------------------------------------------------
    ! Set auxiliary variables in 'nudging_config'
    !---------------------------------------------------

    DO jg = 1, n_dom
      ! Set thresholds for message output
      nudging_config(jg)%msg_thr%low  = 4   ! For important output
      nudging_config(jg)%msg_thr%med  = 10  ! For less important output
      nudging_config(jg)%msg_thr%high = 12  ! For least important output
    
      ! Set thresholds for timers
      nudging_config(jg)%tmr_thr%low  = 1
      nudging_config(jg)%tmr_thr%med  = 5
      nudging_config(jg)%tmr_thr%high = 10
    ENDDO

    !---------------------------------------------------
    ! Initialize derived parameters in 'nudging_config'
    !---------------------------------------------------

    DO jg = 1, n_dom
      nudging_config(jg)%ltype(:)          = .FALSE.
      nudging_config(jg)%ilev_start        = -1
      nudging_config(jg)%ilev_end          = -1
      nudging_config(jg)%lvar(:)           = .FALSE.
    ENDDO

    lremove_qv = .FALSE.

    !---------------------------------------------------
    ! Determine derived parameters
    !---------------------------------------------------

    DO jg = 1,n_dom

      IF (nudging_config(jg)%lnudging) THEN

        ! Set switch for selected nudging type
        nudging_config(jg)%ltype(nudging_config(jg)%nudge_type) = .TRUE.


        ! Vertical start and end levels are determined for DOM 1 only, and lateron 
        ! adapted for eventually vertically nested domains by taking into account 
        ! the vertical shift parameter nshift_total(jg).
        IF (jg==1) THEN

          ! Evaluate nudging start and end heights:
          ! Check status of 'vct_a'
          IF (.NOT. ALLOCATED(vct_a)) CALL finish(TRIM(routine), 'vct_a still uninitialized.')
          nlevp1    = nlev + 1
          lfound    = .FALSE.

          ! 1) Start height:
          start_height = nudging_config(jg)%nudge_start_height
          DO jk = nlevp1, 1, -1
            ! (No shift of vertical index necessary, because this is done for the primary domain, only)
            ! Nominal height of half level
            ! (Note: we (arbitrarily) defined that those vertical grid layers, 
            ! in which 'start/end_height' lie, should be subject to nudging. 
            ! This is, why we do not use 'height = 0.5_wp * (vct_a(jk) + vct_a(jk+1))' here.)
            height = vct_a(jk)
            IF (height > start_height) THEN
              ! We have found the lowermost grid layer, where nudging is applied
              ! (Because vertical looping starts from model top, 'ilev_end' corresponds 
              ! to the start height, and 'ilev_start' corresponds to the end height)
              nudging_config(jg)%ilev_end = jk
              ! Indicate the find
              lfound = .TRUE.
              EXIT 
            ENDIF
          ENDDO  !jk
          IF (.NOT. lfound) CALL finish(TRIM(routine), 'Could not find ilev_end corresponding to nudge_start_height.')
          ! Just to make sure
          nudging_config(jg)%ilev_end = MIN( MAX( 1, nudging_config(jg)%ilev_end ), nlev )
          !
          ! 2) End height
          end_height = nudging_config(jg)%nudge_end_height
          lfound     = .FALSE.
          DO jk = 1, nudging_config(jg)%ilev_end + 1
            height = vct_a(jk)
            IF (height < end_height) THEN
              ! We have found the uppermost grid layer, where nudging is applied
              nudging_config(jg)%ilev_start = jk - 1
              ! Indicate the find
              lfound = .TRUE.
              EXIT
            ENDIF
          ENDDO
          IF (.NOT. lfound) CALL finish(TRIM(routine), 'Could not find ilev_start corresponding to nudge_end_height.')
          ! Just to make sure
          nudging_config(jg)%ilev_start = MIN( MAX( 1, nudging_config(jg)%ilev_start ), nlev )
          ! We replace the namelist entry for end height by the effective end height
          nudging_config(jg)%nudge_end_height   = 0.5_wp * ( vct_a(nudging_config(jg)%ilev_start)     &
            &                                          + vct_a(nudging_config(jg)%ilev_start + 1) )
          ! Finally, we can compute the nudging scale height 
          ! (if there is no explicit namelist input for it, or 'nudge_profile = 1' is selected)
          IF ( nudging_config(jg)%nudge_scale_height < 0._wp .OR.        & 
            &  nudging_config(jg)%nudge_profile == indg_profile%sqrddist ) THEN
            nudging_config(jg)%nudge_scale_height = nudging_config(jg)%nudge_end_height - nudging_config(jg)%nudge_start_height
          ENDIF
          !
        ELSE  ! for jg > 1
          !
          ! copy from DOM 1
          nudging_config(jg)%nudge_end_height   = nudging_config(1)%nudge_end_height
          nudging_config(jg)%nudge_scale_height = nudging_config(1)%nudge_scale_height
          !
          IF (nudging_config(1)%ilev_end <= nshift_total(jg)) THEN
            ! switch off nudging implicitly, as the domain top height is below the nudging start height.
            nudging_config(jg)%ilev_start = 0
            nudging_config(jg)%ilev_end   = -1
          ELSE
            nudging_config(jg)%ilev_start = MAX(1,nudging_config(1)%ilev_start-nshift_total(jg))
            nudging_config(jg)%ilev_end   = MAX(1,nudging_config(1)%ilev_end  -nshift_total(jg))
          ENDIF
        ENDIF


        ! For global nudging only

        IF (nudging_config(jg)%nudge_type == indg_type%globn) THEN
          ! Evaluate the namelist input for the variables that are to be nudged
          ! (please note that the evaluation is very rudimentary without any sophisticated error handling)
          IF (LEN_TRIM(nudging_config(jg)%nudge_var) == 0) THEN 
            ! The string list with the variables that shall be nudged is empty
            CALL finish(routine, "nudge_var is empty.")
          ELSEIF (INDEX(TRIM(nudging_config(jg)%nudge_var), TRIM(cndg_var(indg_var%all))) > 0) THEN
            ! All variables, for which nudging has been implemented, are subject to nudging
            nudging_config(jg)%lvar(:) = .TRUE.
          ELSE
            ! In this case, the string should contain either the character identifier of a single variable, 
            ! or a comma-separated list with the character identifiers of several variables
            lfound = .FALSE.
            DO jvar = 1, indg_var%nitem
              IF (INDEX(TRIM(nudging_config(jg)%nudge_var), TRIM(cndg_var(jvar))) > 0) THEN
                ! Current variable 'jvar' shall be nudged
                nudging_config(jg)%lvar(jvar) = .TRUE.
                lfound                        = .TRUE.
              ENDIF
            ENDDO  !jvar
            IF (.NOT. lfound) CALL finish(routine, "No valid variable identifier found in nudge_var.")
          ENDIF  !Sift variables
          ! The model setup might not allow to nudge water vapour, although the user requests it
          IF (nudging_config(jg)%lvar(indg_var%qv) .AND. .NOT. nudging_config(jg)%lqv_nudgable) THEN
            nudging_config(jg)%lvar(indg_var%qv) = .FALSE.
            lremove_qv                           = .TRUE.
          ENDIF
          ! Include global nudging in timer monitoring?
          nudging_config(jg)%ltimer = timers_level > nudging_config(jg)%tmr_thr%med
          ! Message output during the integration stepping?
          nudging_config(jg)%lmessage = msg_level >= nudging_config(jg)%msg_thr%high
        ENDIF  !IF (nudging_config(jg)%nudge_type == indg_type%globn)

        !---------------------------------------------------
        ! Print some info
        !---------------------------------------------------

        IF (msg_level >= nudging_config(jg)%msg_thr%med) THEN 

          ! Nudging type
          SELECT CASE(nudging_config(jg)%nudge_type)
          CASE(indg_type%ubn)
            c4print = "upper boundary nudging for limited-area mode"
          CASE(indg_type%globn)
            c4print = "global nudging"
          END SELECT
          WRITE(message_text,'(a,i2)') 'Nudging is switched on for DOM', jg
          CALL message(TRIM(routine), message_text)
          WRITE(message_text,'(a)') 'Selected nudging type: '//TRIM(c4print)
          CALL message(TRIM(routine), message_text)

          ! Profile of nudging strength
          SELECT CASE(nudging_config(jg)%nudge_profile)
          CASE(indg_profile%sqrddist)
            c4print = "inverse squared scaled vertical distance from nudging start height"
          CASE(indg_profile%const)
            c4print = "constant"
          CASE(indg_profile%tanh)
            c4print = "hyperbolic tangent decreasing in magnitude from nudging end height downwards"
          CASE(indg_profile%trapid) 
            c4print = "trapezoidal"
          END SELECT
          WRITE(message_text,'(a)') 'Selected profile of nudging strength between start and end height: ' &
            & //TRIM(c4print)
          CALL message(TRIM(routine), message_text)

          ! Nudging start and end heights and the scale height
          WRITE(message_text,'(a)') 'Nudging includes the full levels between:'
          CALL message(TRIM(routine), message_text)
          jk     = nudging_config(jg)%ilev_start
          height = 0.5_wp * ( vct_a(jk) + vct_a(jk+1) )  
          WRITE(message_text,'(a)') ' - Level jk = '//TRIM(int2string(jk))//', height z(jk) = ' &
            & //TRIM(real2string(height))//' m'
          CALL message(' ', message_text)
          jk     = nudging_config(jg)%ilev_end
          height = 0.5_wp * ( vct_a(jk) + vct_a(jk+1) )  
          WRITE(message_text,'(a)') ' - Level jk = '//TRIM(int2string(jk))//', height z(jk) = ' &
            & //TRIM(real2string(height))//' m'
          CALL message(' ', message_text)
          IF (ANY((/indg_profile%sqrddist, indg_profile%tanh, indg_profile%trapid/) &
            & == nudging_config(jg)%nudge_profile)) THEN 
            WRITE(message_text,'(a)') 'Scale height for profile of nudging strength: ' &
              & //TRIM(real2string(nudging_config(jg)%nudge_scale_height))//' m'
            CALL message(TRIM(routine), message_text)
          ENDIF

          ! For upper boundary nudging only

          IF (nudging_config(jg)%nudge_type == indg_type%ubn) THEN

            ! Max. nudging coefficients
            WRITE(message_text,'(a)') 'Max. nudging coefficient for horizontal wind: ' &
              & //TRIM(real2string(nudging_config(jg)%max_nudge_coeff_vn))
            CALL message(TRIM(routine), message_text)
            WRITE(message_text,'(a)') 'Max. nudging coefficient for thermodynamic variables: ' &
              & //TRIM(real2string(nudging_config(jg)%max_nudge_coeff_thermdyn))
            CALL message(TRIM(routine), message_text)

          ! For global nudging only

          ELSEIF (nudging_config(jg)%nudge_type == indg_type%globn) THEN

            ! Nudging variables and max. nudging coefficients
            WRITE(message_text,'(a)') 'Active nudging variables:'
            CALL message(TRIM(routine), message_text)
            DO jvar = 1, indg_var%nitem
              IF (nudging_config(jg)%lvar(jvar)) THEN
                IF (jvar == indg_var%vn) THEN
                  c4print = "horizontal wind (max. nudging coefficient: " &
                    & //TRIM(real2string(nudging_config(jg)%max_nudge_coeff_vn))//")"
                ELSEIF (jvar == indg_var%thermdyn .AND. &
                  &     nudging_config(jg)%thermdyn_type == ithermdyn_type%hydrostatic) THEN
                  c4print = "hydrostatically balanced pressure & temperature (max. nudging coefficient: " &
                    & //TRIM(real2string(nudging_config(jg)%max_nudge_coeff_thermdyn))//")"
                ELSEIF (jvar == indg_var%thermdyn .AND. &
                  &     nudging_config(jg)%thermdyn_type == ithermdyn_type%nonhydrostatic) THEN
                  c4print = "density & virtual potential temperature (max. nudging coefficient: " &
                    & //TRIM(real2string(nudging_config(jg)%max_nudge_coeff_thermdyn))//")"
                ELSEIF (jvar == indg_var%qv) THEN
                  c4print = "water vapour (max. nudging coefficient: " &
                    & //TRIM(real2string(nudging_config(jg)%max_nudge_coeff_qv))//")"
                ENDIF
                WRITE(message_text,'(a)') '- '//TRIM(c4print)
                CALL message(' ', message_text)
              ENDIF  !IF (nudging_config(jg)%lvar(jvar))
            ENDDO  !jvar
            IF (lremove_qv) THEN
              ! Although requested, water vapour may not be available for nudging
              WRITE(message_text,'(a)') 'Note: qv has been removed from the list of nudging variables, ...' 
              CALL message(TRIM(routine), message_text)
              WRITE(message_text,'(a)') '... since the current model setup does not allow to nudge it.'
              CALL message(TRIM(routine), message_text)
            ENDIF  !IF (lremove_qv)

          ENDIF  !IF (nudging_config(jg)%nudge_type == indg_type%...)

        ENDIF  !IF (msg_level >= nudging_config(jg)%msg_thr%med)

      ENDIF  !IF (nudging_config(jg)%lnudging)

      ! Indicate that the setup of 'nudging_config(jg)' took place
      nudging_config(jg)%lconfigured = .TRUE.

    ENDDO  ! jg

    ! Convenience substitute for 'nudging_config(1)%ltype(indg_type%globn)'
    ! for use in 'src/atm_dyn_iconam/mo_nh_stepping'
    ! (and 'src/drivers/mo_atmo_nonhydrostatic: destruct_atmo_nonhydrostatic')
    l_global_nudging = nudging_config(1)%ltype(indg_type%globn)


  END SUBROUTINE configure_nudging

END MODULE mo_nudging_config

!............................................................................................................

  !--------------------------------------------
  !   Visualization of the nudging profiles
  !--------------------------------------------

  ! For upper-boundary nudging 
  ! and global nudging:

  !--------------------------------------------
  ! 1st nudge_profile = indg_profile%sqrddist
  !--------------------------------------------

  !
  !               nominal height z
  ! 
  !                     /|  !                      |

  !                      |
  !         top_height ---    *                         . 
  !                      |    *                         . 
  !                      |    *                         . 
  !                      |    *                         .
  !   nudge_end_height --- .............................*..  --- (upper-boundary nudging: nudge_end_height = top_height)
  !                      |                        *     .     |
  !                      |                   *          .     |
  !                      |               *              .     |
  !                      |            *                 .     |
  !                      |          *                   .     |
  !                      |        *                     .     | < nudge_scale_height = nudge_end_height - nudge_start_height
  !                      |       *                      .     | 
  !                      |      *                       .     |  
  !                      |     *                        .     |
  !                      |     *                        .     |
  !  nudge_start_height --- ..*............................  ---
  !                      |    *                         .
  !                      |    *                         .
  !                      |    *                         .
  !                      |    *                         .
  !                   0 ---   *                         .
  !                      |____|_________________________|____\  nudging coefficient /
  !                           |                         |    /    max_nudge_coeff
  !                           0                         1

  ! For global nudging:

  !--------------------------------------------
  ! 2nd nudge_profile = indg_profile%const
  !--------------------------------------------

  !
  !               nominal height z
  ! 
  !                     /|  !                      |

  !                      |
  !         top_height ---    *                         .
  !                      |    *                         . 
  !                      |    *                         . 
  !                      |    *                         .
  !   nudge_end_height --- .............................*..  
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !  nudge_start_height --- ............................*..  
  !                      |    *                         .
  !                      |    *                         .
  !                      |    *                         .
  !                      |    *                         .
  !                   0 ---   *                         .
  !                      |____|_________________________|____\  nudging coefficient /
  !                           |                         |    /    max_nudge_coeff
  !                           0                         1

  !--------------------------------------------
  ! 3rd nudge_profile = indg_profile%tanh
  !--------------------------------------------

  !
  !               nominal height z
  ! 
  !                     /|  !                      |

  !                      |
  !         top_height ---    *                         .
  !                      |    *                         . 
  !                      |    *                         . 
  !                      |    *                         .
  !   nudge_end_height --- .............................*..  ---
  !                      |                         *    .     |
  !                      |                     *        .     |
  !                      |                  *           .     |
  !                      |                *             .     |
  !                      |              *               .     |
  !                      |             *                .     |
  !                      |            *                 .     |      
  !                      |           *                  .     | < nudge_scale_height = namelist input
  !                      |          *                   .     |
  !                      |          *                   .     |
  !  nudge_start_height --- ........*.......................  |
  !                      |    *                         .     |
  !                      |    *                         .     |
  !                      |    *                         .     |
  !                      |    *                         .    ...
  !                   0 ---   *                         .
  !                      |____|_________________________|____\  nudging coefficient /
  !                           |                         |    /    max_nudge_coeff
  !                           0                         1

  !--------------------------------------------
  ! 4th nudge_profile = indg_profile%trapid
  !--------------------------------------------

  !
  !               nominal height z
  ! 
  !                     /|  !                      |

  !                      |
  !         top_height ---    *                         .
  !                      |    *                         . 
  !                      |    *                         . 
  !                      |    *                         .
  !   nudge_end_height --- ...*............................  ---  
  !                      |     *                        .     |
  !                      |       *                      .     |
  !                      |          *                   .     |
  !                      |               *              .     | < nudge_scale_height = namelist input
  !                      |                     *        .     |
  !                      |                         *    .     |
  !                      |                            * .     | 
  !                      |                             *.    ~~~          
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                              *     
  !                      |                             *.    ~~~      
  !                      |                            * .     |     
  !                      |                         *    .     |      
  !                      |                     *        .     |      
  !                      |               *              .     | < nudge_scale_height = namelist input     
  !                      |          *                   .     |      
  !                      |       *                      .     |      
  !                      |     *                        .     |      
  !  nudge_start_height --- ..*............................  ---  
  !                      |    *                         .
  !                      |    *                         .
  !                      |    *                         .
  !                      |    *                         .
  !                   0 ---   *                         .
  !                      |____|_________________________|____\  nudging coefficient /
  !                           |                         |    /    max_nudge_coeff
  !                           0                         1

!............................................................................................................
