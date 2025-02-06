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

MODULE mo_var

  USE mo_kind,               ONLY: dp, sp
  USE mo_var_metadata_types, ONLY: t_var_metadata, t_var_metadata_dynamic, &
                                   VINTP_TYPE_LIST
  USE mo_exception,          ONLY: message, finish, message_text
  USE mo_impl_constants,     ONLY: REAL_T, SINGLE_T, BOOL_T, INT_T, &
    &                              STR_HINTP_TYPE
  USE mo_var_groups,         ONLY: var_groups_dyn

  IMPLICIT NONE
  PRIVATE

  ! export index constants: model/pressure/height levels
  PUBLIC :: level_type_ml, level_type_pl, level_type_hl, level_type_il, lev_type_str
  PUBLIC :: t_var, t_var_ptr

  ! constants defining level type:
  INTEGER, PARAMETER             :: level_type_ml = 1
  INTEGER, PARAMETER             :: level_type_pl = 2
  INTEGER, PARAMETER             :: level_type_hl = 3
  INTEGER, PARAMETER             :: level_type_il = 4
  ! string defining level type:
  CHARACTER(LEN=2), PARAMETER    :: lev_type_str(4) = (/ 'ML', 'PL', 'HL', 'IL' /)

  TYPE t_var
    REAL(dp), CONTIGUOUS, POINTER :: r_ptr(:,:,:,:,:) => NULL()
    REAL(sp), CONTIGUOUS, POINTER :: s_ptr(:,:,:,:,:) => NULL()
    INTEGER,  CONTIGUOUS, POINTER :: i_ptr(:,:,:,:,:) => NULL()
    LOGICAL,  CONTIGUOUS, POINTER :: l_ptr(:,:,:,:,:) => NULL()
    INTEGER                      :: var_base_size      ! generic size in bytes of variable used
    TYPE(t_var_metadata)         :: info               ! meta data for this entry
    TYPE(t_var_metadata_dynamic) :: info_dyn           ! dynamic meta data for this entry (see type description)
    TYPE(t_var), POINTER :: ref_to => NULL()
  CONTAINS
    PROCEDURE :: print_short => print_var_short
    PROCEDURE :: print_rigorous => print_var_rigorous
  END TYPE t_var

  TYPE t_var_ptr
    TYPE(t_var), POINTER :: p => NULL()
  END TYPE t_var_ptr

  CHARACTER(*), PARAMETER :: modname = "mo_var"

CONTAINS

  SUBROUTINE print_var_short(var)
    CLASS(t_var), INTENT(IN), TARGET :: var
    CHARACTER(LEN=1), PARAMETER :: lm3(3) = ['i', 'm', 'a'], &
                                   lm4(3) = ['c', 'v', 'e']
    CHARACTER(LEN=4) :: localMode
    TYPE(t_var_metadata), POINTER :: info

    info => var%info
    localMode = '----'
    IF (info%lrestart) localMode(1:1) = 'r'
    IF (info%lcontained) localMode(2:2) = 't'
    SELECT CASE (info%isteptype)
    CASE (1,2,3)
      localMode(3:3) = lm3(info%isteptype)
    END SELECT
    SELECT CASE (info%hgrid)
    CASE (1,2,3)
      localMode(4:4) = lm4(info%hgrid)
    CASE (45)
      localMode(4:4) = 'L'
    END SELECT
    WRITE(message_text, '(a4,3i4,a24,a40)') localMode, &
      & info%grib2%discipline, info%grib2%category,   &
      & info%grib2%number, TRIM(info%name), TRIM(info%cf%standard_name)
    CALL message('', message_text)
  END SUBROUTINE print_var_short

  SUBROUTINE print_var_rigorous(var)
    CLASS(t_var), INTENT(IN), TARGET :: var
    TYPE(t_var_metadata), POINTER :: info
    TYPE(t_var_metadata_dynamic), POINTER :: info_dyn
    CHARACTER(*), PARAMETER :: routine = modname//":print_var_rigorous", &
      & category(22) = [ &
      & 'Table entry name                            : ', &
      & 'Pointer status                              : ', &
      & 'Local field dimensions                      : ', &
      & 'Assigned GRIB discipline/category/parameter : ', &
      & 'CF convention standard name/unit            : ', &
      & 'CF convention long name                     : ', &
      & 'Field is in a container                     : ', &
      & ' Index in container                         : ', &
      & ' horizontal grid type used (C=1,V=2,E=3)    : ', &
      & ' vertical grid type used (see cdilib.c)     : ', &
      & ' stat. processing type (I=1,AVG=2,ACC=3...) : ', &
      & 'Missing value                               : ', &
      & 'Added to restart                            : ', &
      & 'Tracer field                                : ', &
      & '   Child-to-parent feedback                 : ', &
      & '   Horizontal transport method              : ', &
      & '   Vertical transport method                : ', &
      & '   Turbulent transport                      : ', &
      & 'Variable class/species                      : ', &
      & 'Variable group(s)                           : ', &
      & 'Horizontal interpolation                    : ', &
      & 'Vertical interpolation                      : ' ]
    CHARACTER(LEN=2) :: istr
    CHARACTER(LEN=20) :: mvstr
    CHARACTER(LEN=9) :: ipstr
    CHARACTER(LEN=32) :: dimstr
    INTEGER :: i, mlen, gnlen

    info => var%info
    info_dyn => var%info_dyn
    CALL message('', category(1)//TRIM(info%name))
    IF (ANY([ASSOCIATED(var%r_ptr), ASSOCIATED(var%s_ptr), &
             ASSOCIATED(var%i_ptr), ASSOCIATED(var%l_ptr)])) THEN
      CALL message ('',category(2)//'in use.')
      dimstr = '('
      DO i = 1, info%ndims
        WRITE(dimstr,"(a,i0,a1)") TRIM(dimstr), info%used_dimensions(i), &
          & MERGE(')', ',', info%ndims .EQ. i)
      ENDDO
      CALL message('', category(3)//TRIM(dimstr))
    ELSE
      CALL message('', category(2)//'not in use.')
    ENDIF
    WRITE (message_text,'(a,3i4)') category(4), &
      & info%grib2%discipline, info%grib2%category, info%grib2%number
    CALL message('', message_text)
    WRITE (message_text,'(a,a,a,a)') category(5), &
         TRIM(info%cf%standard_name), '     ', TRIM(info%cf%units)
    CALL message('', message_text)
    CALL message('', category(6)//TRIM(info%cf%long_name))
    CALL message('', category(7)//MERGE('yes.', 'no. ', info%lcontained))
    istr = '--'
    IF (info%lcontained) WRITE(istr, "(i2)") info%ncontained
    CALL message('', category(8)//istr)
    WRITE(istr, "(i2)") info%hgrid
    CALL message('', category(9)//istr)
    WRITE(istr, "(i2)") info%vgrid
    CALL message('', category(10)//istr)
    WRITE(istr, "(i2)") info%isteptype
    CALL message('', category(11)//istr)
    IF (info%lmiss) THEN
      SELECT CASE(info%data_type)
      CASE (REAL_T)
        WRITE(mvstr, "(e20.12)") info%missval%rval
      CASE (SINGLE_T)
        WRITE(mvstr, "(e20.12)") info%missval%sval
      CASE (INT_T)
        WRITE(mvstr, "(i8)") info%missval%ival
      CASE (BOOL_T)
        WRITE(mvstr, "(l8)") info%missval%lval
      END SELECT
    ELSE
      WRITE(mvstr, "(a)") 'off.'
    ENDIF
    CALL message('', category(12)//mvstr)
    CALL message('', category(13)//MERGE("yes.", " no.", info%lrestart))
    CALL message('', category(14)//MERGE("yes.", " no.", info_dyn%tracer%lis_tracer))
    IF (info_dyn%tracer%lis_tracer) THEN
      CALL message('', category(15)//MERGE("yes.", " no.", info_dyn%tracer%lfeedback))
      WRITE(ipstr, "(3i3)") info_dyn%tracer%ihadv_tracer
      CALL message('', category(16)//ipstr)
      WRITE(ipstr, "(3i3)") info_dyn%tracer%ivadv_tracer
      CALL message('', category(17)//ipstr)
      CALL message('', category(18)//MERGE("yes.", " no.", info_dyn%tracer%lturb_tracer))
    ENDIF !lis_tracer
    WRITE(istr, "(i2)") info%var_class
    CALL message('', category(19)//istr)
    IF (ANY(info%in_group(:))) THEN
      CALL message('', category(20))
      message_text = ' + '
      mlen = 3
      DO i = 1, var_groups_dyn%get_n_grps()
        IF (.NOT.info%in_group(i)) CYCLE
        gnlen = var_groups_dyn%gname_len(i)
        IF (mlen + 2 + gnlen .GT. LEN(message_text)) CALL finish(routine, "message_text overflow")
        WRITE(message_text, "(3a)") message_text(1:mlen), var_groups_dyn%gname(i)(1:gnlen), ", "
        mlen = mlen + 2 + gnlen
      END DO
      CALL message('', message_text(1:mlen-2))
    END IF
    CALL message('', category(21)//STR_HINTP_TYPE(info%hor_interp%hor_intp_type))
    DO i = 1, SIZE(VINTP_TYPE_LIST)
      IF (info%vert_interp%vert_intp_type(i)) &
        CALL message('', category(22)//VINTP_TYPE_LIST(i))
    END DO
    CALL message('', '')
  END SUBROUTINE print_var_rigorous

END MODULE mo_var
