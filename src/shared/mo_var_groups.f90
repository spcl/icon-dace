! Group information for ICON variables.
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

MODULE mo_var_groups

  USE mo_impl_constants,        ONLY: vname_len
  USE mo_util_string,           ONLY: toupper
  USE mo_key_value_store,       ONLY: t_key_value_store
  USE mo_exception,             ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MAX_GROUPS, var_groups_dyn, groups

  !> module name string
  CHARACTER(*), PARAMETER :: modname = 'mo_var_groups'

  ! maximum number of variable groups supported by a single info state
  INTEGER, PARAMETER :: MAX_GROUPS = 200

  ! ---------------------------------------------------------------
  ! STATICALLY DEFINED VARIABLE GROUPS
  ! ---------------------------------------------------------------

  ! List of *static* variable groups.
  ! 
  ! A variable can have any combination of this which means that it is
  ! part of each of these different variable sets.
  ! A variable is added to an existing group by setting the meta-data
  ! information "in_group" as follows
  !
  !   CALL add_var( p_prog_list, ..., in_group=groups("nh_prog_vars") )
  !
  ! It is also possible to add a variable to more than one group:
  !
  !   CALL add_var( diag_list, ...,   &
  !                 in_group=groups("multisnow_vars", "snow_vars"))
  !
  ! New groups can be added by extending the VAR_GROUPS list.
  !
  ! Note that the statically defined group list "var_groups" is
  ! non-public. Its contents are copied to a dynamically growing list
  ! "var_groups_dyn".
  INTEGER, PARAMETER :: N_VAR_GROUPS_STATIC = 73

  CHARACTER(LEN=vname_len), PARAMETER :: VAR_GROUPS_STATIC(N_VAR_GROUPS_STATIC) = &
     [ "ALL                   ",  &
    &  "ATMO_ML_VARS          ",  &
    &  "ATMO_PL_VARS          ",  &
    &  "ATMO_ZL_VARS          ",  &
    &  "NH_PROG_VARS          ",  &
    &  "ATMO_DERIVED_VARS     ",  &
    &  "RAD_VARS              ",  &
    &  "PRECIP_VARS           ",  &
    &  "CLOUD_DIAG            ",  &
    &  "COMIN_TENDENCIES      ",  &
    &  "COMIN_VARS            ",  &
    &  "PBL_VARS              ",  &
    &  "PHYS_TENDENCIES       ",  &
    &  "PROG_TIMEMEAN         ",  &
    &  "ECHAM_TIMEMEAN        ",  &
    &  "TRACER_TIMEMEAN       ",  &
    &  "ATMO_TIMEMEAN         ",  &
    &  "LAND_VARS             ",  &
    &  "LAND_TILE_VARS        ",  &
    &  "MULTISNOW_VARS        ",  &
    &  "ADDITIONAL_PRECIP_VARS",  &
    &  "SNOW_VARS             ",  &
    &  "DWD_FG_ATM_VARS       ",  &  ! DWD First Guess (atmosphere) 
    &  "DWD_FG_SFC_VARS       ",  &  ! DWD First Guess (surface/soil)
    &  "DWD_FG_SFC_VARS_T     ",  &  ! DWD First Guess (surface/soil) tiles
    &  "MODE_DWD_FG_IN        ",  &  ! Input first guess fields for MODE_DWD
    &  "MODE_DWD_ANA_IN       ",  &  ! Input analysis fields for MODE_DWD
    &  "MODE_IAU_FG_IN        ",  &  ! First guess input for IAU
    &  "MODE_IAU_ANA_IN       ",  &  ! Analysis input for IAU
    &  "MODE_IAU_ANAATM_IN    ",  &  ! Atmospheric analysis input for (old/new) IAU
    &  "MODE_IAU_OLD_FG_IN    ",  &  ! First guess input for old IAU mode
    &  "MODE_IAU_OLD_ANA_IN   ",  &  ! Analysis input for old IAU mode
    &  "MODE_COMBINED_IN      ",  &  ! Input fields for MODE_COMBINED
    &  "MODE_COSMO_IN         ",  &  ! Input fields for MODE_COSMO
    &  "OCE_PROG              ",  &
    &  "OCE_DIAG              ",  &
    &  "OCE_EDDY              ",  &
    &  "OCE_DEFAULT           ",  &
    &  "OCEAN_MOC             ",  &  ! meant o hold all kinds of overturning fields (atl, pac, global,...)
    &  "OCEAN_FLOWS           ",  &  ! meant o hold all through flows
    &  "HAMOCC_BASE           ",  &
    &  "HAMOCC_TEND           ",  &
    &  "HAMOCC_MONI           ",  &
    &  "HAMOCC_SED            ",  &
    &  "oce_essentials        ",  &
    &  "oce_force_essentials  ",  &
    &  "OCE_AUX               ",  &
    &  "OCEAN_MONITOR         ",  &
    &  "OCE_GEOMETRY          ",  &
    &  "OCE_PHYSICS           ",  &
    &  "OCE_COEFFS            ",  &
    &  "ICE_DEFAULT           ",  &
    &  "ICE_BUDGETS           ",  &
    &  "ICE_DIAG              ",  &
    &  "LATBC_PREFETCH_VARS   ",  &
    &  "mode_iniana           ",  &  ! Variable set needed to initialize ICON (MODE_ICONVREMAP)
    &  "icon_lbc_vars         ",  &  ! Variable set needed for ICON-LAM lateral boundary conditions
    &  "ART_AEROSOL           ",  &  ! ICON-ART fields for aerosol particles
    &  "ART_CHEMISTRY         ",  &  ! ICON-ART fields for chemical tracers
    &  "ART_FPLUME            ",  &  ! ICON-ART fields for FPlume output
    &  "ART_DIAGNOSTICS       ",  &  ! ICON-ART fields for diagnostic fields
    &  "ART_ROUTINE_DIAG      ",  &  ! ICON-ART fields for routine diagnostic fields
    &  "RTTOV                 ",  &
    &  "UPATMO_TENDENCIES     ",  &  ! Upper-atmosphere physics tendencies
    &  "UPATMO_RAD_GASES      ",  &  ! Upper-atmosphere radiatively active gases
    &  "IAU_RESTORE_VARS      ",  &  ! Dynamics + NWP Physics fields for the IAU save/restore mechanism
    &  "IAU_INIT_VARS         ",  &  ! Dynamics + NWP Physics fields for the IAU reinitialization mechanism
    &  "WAVE_SPECTRUM         ",  &  ! ICON-waves: full wave spectrum output
    &  "WIND_10M              ",  &  ! ICON-waves: wind 10m
    &  "WAVE_PHY              ",  &  ! ICON-waves: wave physics
    &  "WAVE_PHY_EXT          ",  &  ! ICON-waves: extended wave physics
    &  "WAVE_SHORT            ",  &  ! ICON-waves: short/"standart"
    &  "WAVE_FORCING          "   ]  ! ICON-waves: wave forcing

  ! ---------------------------------------------------------------
  ! DYNAMICALLY DEFINED VARIABLE GROUPS
  ! ---------------------------------------------------------------
  TYPE t_var_groups
    ! variable group names. The list ordering is important, since the
    ! entry #i (= ID #i) corresponds to the i'th LOGICAL entry in the
    ! "in_group" of the t_var_metadata type below.
    CHARACTER(LEN=vname_len) :: gname(MAX_GROUPS) = ""
    CHARACTER(LEN=vname_len) :: gname_upper(MAX_GROUPS) = ""
    INTEGER :: gname_len(MAX_GROUPS) = 0
    INTEGER, PRIVATE :: n_grps = 0
    TYPE(t_key_value_store), PRIVATE :: map
  CONTAINS
    PROCEDURE :: group_id
    PROCEDURE :: get_n_grps
  END TYPE t_var_groups

  ! List of variable groups - including the statically defined group
  ! names from VAR_GROUPS as well as dynamic variable groups, used,
  ! e.g., for tiles.
  TYPE(t_var_groups) :: var_groups_dyn

  INTERFACE groups
    MODULE PROCEDURE groups_arg
    MODULE PROCEDURE groups_vec
  END INTERFACE

CONTAINS

  INTEGER FUNCTION get_n_grps(this)
    CLASS(t_var_groups), INTENT(IN) :: this

    get_n_grps = this%n_grps
  END FUNCTION get_n_grps

  !----------------------------------------------------------------------------------------
  !> Implements a (somewhat randomly chosen) one-to-one mapping
  !  between a string and an integer ID number between 1 and
  !  size(var_groups%name)+1.
  INTEGER FUNCTION group_id(this, in_str)
    CLASS(t_var_groups), INTENT(INOUT) :: this
    CHARACTER(*), INTENT(IN) :: in_str
    INTEGER :: ierr
    CHARACTER(*), PARAMETER :: routine = modname//":group_id"

    IF (this%n_grps .EQ. 0) CALL init()
    ! search the variable groups (which includes the statically
    ! defined groups and the dynamically defined ones):
    CALL this%map%get(in_str, group_id, opt_err=ierr)
    ! If the group does not exist, create it.
    IF (ierr .NE. 0) CALL append() 
  CONTAINS

    SUBROUTINE init()
      INTEGER :: i

      this%n_grps                            = N_VAR_GROUPS_STATIC
      CALL this%map%init(.FALSE.)
      DO i = 1, N_VAR_GROUPS_STATIC
        this%gname_len(i) = LEN_TRIM(VAR_GROUPS_STATIC(i))
        this%gname(i) = VAR_GROUPS_STATIC(i)(1:this%gname_len(i))
        this%gname_upper(i) = toupper(VAR_GROUPS_STATIC(i))
        CALL this%map%put(this%gname(i), i)
      END DO
    END SUBROUTINE init

    SUBROUTINE append()

      IF (this%n_grps .EQ. MAX_GROUPS) CALL finish(routine, "too many groups")
      this%n_grps = this%n_grps + 1
      this%gname(this%n_grps) = toupper(in_str)
      this%gname_upper(this%n_grps) = this%gname(this%n_grps)
      this%gname_len(this%n_grps) = LEN_TRIM(in_str)
      CALL this%map%put(in_str, this%n_grps)
      group_id = this%n_grps
    END  SUBROUTINE append
  END FUNCTION group_id

  !----------------------------------------------------------------------------------------
  !> Utility function with *a lot* of optional string parameters g1,
  !  g2, g3, g4, ...; mapping those onto a
  !  LOGICAL(DIMENSION=MAX_GROUPS) according to the "group_id"
  !  function.
  FUNCTION groups_arg(g01, g02, g03, g04, g05, g06, g07, g08, g09, g10, g11, g12, g13, groups_in)
    LOGICAL :: groups_arg(MAX_GROUPS)
    CHARACTER(*), INTENT(IN), OPTIONAL :: &
      &   g01, g02, g03, g04, g05, g06, g07, g08, g09, g10, g11, g12, g13
    LOGICAL, INTENT(IN), OPTIONAL :: groups_in(MAX_GROUPS)

    groups_arg(1) = .TRUE. ! this is "ALL" - obviously true
    groups_arg(2:) = .FALSE.
    IF (PRESENT(groups_in)) groups_arg(2:) = groups_in(2:)
    IF (PRESENT(g01)) groups_arg(var_groups_dyn%group_id(g01)) = .TRUE.
    IF (PRESENT(g02)) groups_arg(var_groups_dyn%group_id(g02)) = .TRUE.
    IF (PRESENT(g03)) groups_arg(var_groups_dyn%group_id(g03)) = .TRUE.
    IF (PRESENT(g04)) groups_arg(var_groups_dyn%group_id(g04)) = .TRUE.
    IF (PRESENT(g05)) groups_arg(var_groups_dyn%group_id(g05)) = .TRUE.
    IF (PRESENT(g06)) groups_arg(var_groups_dyn%group_id(g06)) = .TRUE.
    IF (PRESENT(g07)) groups_arg(var_groups_dyn%group_id(g07)) = .TRUE.
    IF (PRESENT(g08)) groups_arg(var_groups_dyn%group_id(g08)) = .TRUE.
    IF (PRESENT(g09)) groups_arg(var_groups_dyn%group_id(g09)) = .TRUE.
    IF (PRESENT(g10)) groups_arg(var_groups_dyn%group_id(g10)) = .TRUE.
    IF (PRESENT(g11)) groups_arg(var_groups_dyn%group_id(g11)) = .TRUE.
    IF (PRESENT(g12)) groups_arg(var_groups_dyn%group_id(g12)) = .TRUE.
    IF (PRESENT(g13)) groups_arg(var_groups_dyn%group_id(g13)) = .TRUE.
  END FUNCTION groups_arg

  !----------------------------------------------------------------------------------------
  !> The same, but provide list of groups as one character vector of group names.
  !  Attention: the strings passed in group_list must be of length vname_len !
  FUNCTION groups_vec(group_list)
    LOGICAL :: groups_vec(MAX_GROUPS)
    CHARACTER(LEN=vname_len), INTENT(IN) :: group_list(:)
    INTEGER :: i

    groups_vec(1) = .TRUE. ! this is "ALL" - obviously true
    groups_vec(2:) = .FALSE.
    DO i = 1, SIZE(group_list)
      groups_vec(var_groups_dyn%group_id(group_list(i))) = .TRUE.
    END DO
  END FUNCTION groups_vec

END MODULE mo_var_groups
