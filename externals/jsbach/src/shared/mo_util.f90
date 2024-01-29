!> Contains utility functions
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
MODULE mo_util
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_depth_to_layers_2d, ifs2soil, one_of, toupper, int2string, real2string, logical2string, &
    &       report_memory_usage, soil_init_from_texture

  INTERFACE one_of
    MODULE PROCEDURE one_of_str
    MODULE PROCEDURE one_of_int
  END INTERFACE one_of

  CHARACTER(len=*), PARAMETER :: modname = 'mo_util'

CONTAINS

  !!
  ! !>      Compute soil layer depths from (fixed) layer thicknesses and total depth
  !
  ! @param      depth    Total soil depth
  ! @param      dz       Fixed layer thicknesses
  !
  ! @return     Layer thicknesses until soil depth
  !!
  FUNCTION soil_depth_to_layers_2d(depth, dz) RESULT(depth_l)

    REAL(wp), INTENT(in) :: depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth,1), SIZE(dz), SIZE(depth,2))

    REAL(wp) :: bottom, bottom_m1
    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':soil_depth_to_layer_2d'

    nsoil = SIZE(dz)

    depth_l(:,:,:) = 0._wp
    bottom = 0._wp
    bottom_m1 = 0._wp
    DO i=1,nsoil
      bottom = bottom + dz(i)
      WHERE (depth(:,:) >= bottom)
        depth_l(:,i,:) = dz(i)
      ELSEWHERE (depth(:,:) > bottom_m1)
        depth_l(:,i,:) = depth(:,:) - bottom_m1
      END WHERE
      bottom_m1 = bottom_m1 + dz(i)
    END DO

  END FUNCTION soil_depth_to_layers_2d

  SUBROUTINE ifs2soil(ifs_3d, ifs_depth, soil_3d, soil_depth, ifs_sfc)

    REAL(wp), INTENT(in) :: &
      & ifs_3d (:,:,:),     &
      & ifs_depth(:),       &
      & soil_depth(:)
    REAL(wp), OPTIONAL, INTENT(in) :: &
      & ifs_sfc(:,:)

    REAL(wp), INTENT(out) :: &
      & soil_3d(:,:,:)

    REAL(wp) :: ifs(SIZE(ifs_3d,1),0:SIZE(ifs_depth)+1,SIZE(ifs_3d,3)), zifs(0:SIZE(ifs_depth)+1), zsoil(0:SIZE(soil_depth))
    REAL(wp) :: zwt, zwt_sum
    INTEGER  :: nsoil_in, nsoil_out, isoil_in, isoil_out, ik1

    CHARACTER(len=*), PARAMETER :: routine = modname//':ifs2soil'

    nsoil_in  = SIZE(ifs_depth)
    nsoil_out = SIZE(soil_depth)

    zsoil(1:nsoil_out) = soil_depth(:)
    zsoil(0)           = 0._wp
    zifs(1:nsoil_in) = ifs_depth(:)
    zifs(0)          = 0._wp
    zifs(nsoil_in+1) = soil_depth(nsoil_out)
    ifs(:,1:nsoil_in,:) = ifs_3d(:,:,:)
    IF (PRESENT(ifs_sfc)) THEN
      ifs(:,0         ,:) = ifs_sfc(:,:)
    ELSE
      ifs(:,0         ,:) = ifs_3d(:,1,:)
    END IF
    ifs(:,nsoil_in+1,:) = ifs_3d(:,nsoil_in,:)

    soil_3d(:,:,:) = 0._wp

    IF (PRESENT(ifs_sfc)) THEN
      soil_3d(:,1,:) = ifs_sfc(:,:)
      ik1 = 2
    ELSE
      ik1 = 1
    END IF
    DO isoil_out = ik1, nsoil_out
      zwt_sum        = 0._wp
      DO isoil_in = 1,nsoil_in+1
        IF (zifs(isoil_in) <= zsoil(isoil_out-1)) THEN
          ! Do nothing
          zwt = 0._wp
        ELSE IF (zifs(isoil_in-1) > zsoil(isoil_out)) THEN
          ! Do nothing
          zwt = 0._wp
        ELSE IF (zifs(isoil_in-1) > zsoil(isoil_out-1) .AND. zifs(isoil_in) <= zsoil(isoil_out)) THEN
          zwt = 1._wp
          zwt_sum = zwt_sum + zwt
          soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) + zwt * ifs(:,isoil_in,:)
        ELSE IF (zifs(isoil_in) > zsoil(isoil_out-1) .AND. zifs(isoil_in) <= zsoil(isoil_out)) THEN
          zwt = (zifs(isoil_in) - zsoil(isoil_out-1)) / (zsoil(isoil_out) - zsoil(isoil_out-1))
          zwt_sum = zwt_sum + zwt
          soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) + zwt * ifs(:,isoil_in,:)
        ELSE IF (zifs(isoil_in) > zsoil(isoil_out) .AND. zifs(isoil_in-1) <= zsoil(isoil_out)) THEN
          zwt = (zsoil(isoil_out) - zifs(isoil_in-1)) / (zsoil(isoil_out) - zsoil(isoil_out-1))
          zwt_sum = zwt_sum + zwt
          soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) + zwt * ifs(:,isoil_in,:)
        END IF
        ! print*, 'AAA isoil_out, isoil_in, zwt, zwt_sum', isoil_out, isoil_in, zwt, zwt_sum
      END DO
      IF (zwt_sum < EPSILON(1._wp)) CALL finish(routine, 'zwt_sum = 0 (should not happen)!')
      soil_3d(:,isoil_out,:) = soil_3d(:,isoil_out,:) / zwt_sum
    END DO

  END SUBROUTINE ifs2soil

  !> Calculate soil parameters fom soil texture.
  !!
  !! Input to this routine should conform:
  !! fraction_sand + fraction_silt + fraction_clay = 1
  !! fraction_sand_deep + fraction_silt_deep + fraction_clay_deep = 1
  !! Properties of upper soil and deep soil are identical. Their composition is averaged.
  !! The organic fraction of the soil is the non-mineral fraction of the soil.
  PURE ELEMENTAL SUBROUTINE soil_init_from_texture(                 &
    & parameter_sand, parameter_silt, parameter_clay, parameter_oc, &
    & fraction_sand, fraction_silt, fraction_clay,                  &
    & fraction_sand_deep, fraction_silt_deep, fraction_clay_deep,   &
    & fraction_oc, fraction_oc_deep,                                &
    & parameter_from_texture)

    REAL(wp), INTENT(in)  ::  parameter_sand, parameter_silt, parameter_clay, parameter_oc
    REAL(wp), INTENT(in)  ::                                      &
      fraction_sand, fraction_silt, fraction_clay,                &
      fraction_sand_deep, fraction_silt_deep, fraction_clay_deep, &
      fraction_oc, fraction_oc_deep
    REAL(wp), INTENT(out)  ::  parameter_from_texture

    parameter_from_texture = MAX(1._wp - ((fraction_oc + fraction_oc_deep) &
                             * 0.5_wp),0._wp) *                            &
                             MAX(((fraction_sand + fraction_sand_deep)     &
                             * 0.5_wp * parameter_sand) +                  &
                             ((fraction_silt + fraction_silt_deep)         &
                             * 0.5_wp * parameter_silt) +                  &
                             ((fraction_clay + fraction_clay_deep)         &
                             * 0.5_wp * parameter_clay),                   &
                             MIN(parameter_sand,parameter_silt,parameter_clay)) &
                             + ((fraction_oc + fraction_oc_deep)           &
                             * 0.5_wp) * parameter_oc

  END SUBROUTINE soil_init_from_texture

  !
  !> Conversion: Lowercase -> Uppercase
  !
  !  Copied from ICON module shared/mo_util_string
  !
  FUNCTION toupper (lowercase)
    CHARACTER(len=*), INTENT(in) :: lowercase
    CHARACTER(len=LEN_TRIM(lowercase)) :: toupper
    !
    INTEGER, PARAMETER :: idel = ICHAR('A')-ICHAR('a')
    INTEGER :: i
    !
    DO i = 1, LEN_TRIM(lowercase)
      IF (ICHAR(lowercase(i:i)) >= ICHAR('a') .AND. ICHAR(lowercase(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lowercase(i:i)) + idel )
      ELSE
        toupper(i:i) = lowercase(i:i)
      ENDIF
    ENDDO
    !
  END FUNCTION toupper

  !> Function for convenience
  !
  !  If "in_str" is matching one of the arguments "arg(i)" return the
  !  index "i". Returns "-1" if none of the strings matches.
  !
  !  Copied from ICON module shared/mo_util_string
  !
  FUNCTION one_of_str(in_str, arg)
    INTEGER :: one_of_str
    CHARACTER(len=*), INTENT(IN)           :: in_str    ! input string
    CHARACTER(len=*), INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i, ipos, len_in_str

    len_in_str = LEN_TRIM(in_str)
    one_of_str = -1
    DO i=1,SIZE(arg)
      ipos = SCAN(arg(i), '*')
      IF (ipos == 0) THEN
        ipos = LEN_TRIM(arg(i))
      ELSE
        ipos = ipos - 1
      END IF
      IF (ipos > 0) THEN
        IF ( ADJUSTL(toupper(in_str(1:MIN(ipos,len_in_str)))) == ADJUSTL(toupper(arg(i)(1:ipos))) ) THEN
          one_of_str=i
          EXIT
        END IF
      ELSE
        one_of_str=i
        EXIT
      END IF
    END DO
  END FUNCTION one_of_str

  FUNCTION one_of_int(in_int, arg)
    INTEGER :: one_of_int
    INTEGER, INTENT(IN)           :: in_int    ! input integer
    INTEGER, INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i

    one_of_int = -1
    DO i=1,SIZE(arg)
      IF (in_int == arg(i)) THEN
        one_of_int=i
        EXIT
      END IF
    END DO
  END FUNCTION one_of_int

  !
  ! returns integer n as a string (often needed in printing messages)
  !
  ! Copied from ICON module shared/mo_util_string
  !
  FUNCTION int2string(n, opt_fmt)
    CHARACTER(len=:), ALLOCATABLE :: int2string ! result
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=:), ALLOCATABLE :: fmt
    CHARACTER(len=30) :: string

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(I10)'
    END IF
    WRITE(string,fmt) n
    int2string = TRIM(ADJUSTL(string))
    !
  END FUNCTION int2string

  FUNCTION real2string(x, opt_fmt)
    CHARACTER(len=:), ALLOCATABLE :: real2string ! result
    REAL(wp), INTENT(in) :: x
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=:), ALLOCATABLE :: fmt
    CHARACTER(len=30) :: string

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(E10.4)'
    END IF
    WRITE(string,fmt) x
    real2string = TRIM(ADJUSTL(string))
    !
  END FUNCTION real2string

  FUNCTION logical2string(value)
    CHARACTER(len=:), ALLOCATABLE :: logical2string ! result
    LOGICAL, INTENT(in) :: value
    !
    IF (value) THEN
      logical2string = '.TRUE.'
    ELSE
      logical2string = '.FALSE.'
    END IF
    !
  END FUNCTION logical2string

  SUBROUTINE report_memory_usage()

    INTEGER           :: size=0, resident=0, share=0, text=0, lib=0, data=0, dt=0, oom_score=0
    INTEGER           :: ierr=0

    OPEN(unit=1, file='/proc/self/statm', form='formatted', status='OLD', action='READ', iostat=ierr)
    READ(unit=1, fmt=*, iostat=ierr) size, resident, share, text, lib, data, dt
    CLOSE(unit=1)
    OPEN(unit=1, file='/proc/self/oom_score', form='formatted', status='OLD', action='READ', iostat=ierr)
    READ(unit=1, fmt=*, iostat=ierr) oom_score
    CLOSE(unit=1)

    IF (ierr < 0) THEN
      CALL message('', 'Problem reading /proc/self/statm')
    ELSE
      CALL message('INFO -- Memory usage in pages and OOM score', ' size='//TRIM(int2string(size))// &
        &                             ' resident='//TRIM(int2string(resident))// &
        &                             ' data='//TRIM(int2string(data))// &
        &                             ' oom_score='//TRIM(int2string(oom_score)))
    END IF

  END SUBROUTINE report_memory_usage

#endif
END MODULE mo_util
