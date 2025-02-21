!> @file ppm_extents.f90
!! @brief describe rectilinear partitions and partitionings
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
! Keywords: partition partitioning descriptor
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!> contains definition of extent and interval types and associated
!! functions

MODULE ppm_extents
  USE ppm_std_type_kinds, ONLY: i4, i8, dp, sp
  USE ppm_math_extensions, ONLY: de_g_sp, de_g_sp_width, de_g_dp, de_g_dp_width
  IMPLICIT NONE
  PRIVATE
  !> Describes range as start and (directed) size
  TYPE, PUBLIC :: extent
    SEQUENCE
    !> range is anchored by this value
    INTEGER(i4) :: first
    !> range has this size
    INTEGER(i4) :: size
  END TYPE extent

  !> canonical representation of empty range
  TYPE(extent), PUBLIC, PARAMETER :: empty_extent = extent(0, 0)

  !> string representation of extent size/position takes
  !! at most 11 decimal places (10 + sign)
  INTEGER(i4), PARAMETER :: ext_i2s_len=11
  !> string representation of extent takes at most ext2s_len characters
  INTEGER(i4), PUBLIC, PARAMETER :: ext2s_len=2*ext_i2s_len+4

  !> Describes range as start and (directed) size
  TYPE, PUBLIC :: extent_i8
    SEQUENCE
    !> range is anchored by this value
    INTEGER(i8) :: first
    !> range has this size
    INTEGER(i8) :: size
  END TYPE extent_i8

  !> canonical representation of empty range
  TYPE(extent_i8), PUBLIC, PARAMETER :: empty_extent_i8 = extent_i8(0_i8, 0_i8)

  !> string representation of 8-byte extent size/position takes
  !! at most 20 decimal places (10 + sign)
  INTEGER(i4), PARAMETER :: ext_i2s_len_i8=20
  !> string representation of i8 extent takes at most ext2s_len_i8 characters
  INTEGER(i4), PUBLIC, PARAMETER :: ext2s_len_i8=2*ext_i2s_len_i8+4

  !> interval including both limits, i.e. [first,last] where first and last are
  !! integral numbers, this exhibits a problem for 0-size intervals, because
  !! we also want to allow negative stride intervals like .e.g. 2..-3
  !! of size -6, meaning there is no iinterval representation equal to
  !! the empty set
  TYPE, PUBLIC :: iinterval
    SEQUENCE
    !> first delimiter of range
    INTEGER(i4) :: first
    !> last delimiter of range
    INTEGER(i4) :: last
  END TYPE iinterval

  !> canonical representation of empty interval
  TYPE(iinterval), PUBLIC, PARAMETER :: empty_iinterval = iinterval(1, 0)

  !> interval including both limits, i.e. [first,last] where first and last are
  !! integral numbers, this exhibits a problem for 0-size intervals, because
  !! we also want to allow negative stride intervals like .e.g. 2..-3
  !! of size -6, meaning there is no iinterval representation equal to
  !! the empty set
  TYPE, PUBLIC :: iinterval_i8
    SEQUENCE
    !> first delimiter of range
    INTEGER(i8) :: first
    !> last delimiter of range
    INTEGER(i8) :: last
  END TYPE iinterval_i8

  !> canonical representation of empty interval
  TYPE(iinterval_i8), PUBLIC, PARAMETER :: empty_iinterval_i8 &
       = iinterval_i8(1_i8, 0_i8)

  TYPE, PUBLIC :: iinterval_sp
    SEQUENCE
    !> first delimiter of range
    REAL(sp) :: first
    !> last delimiter of range
    REAL(sp) :: last
  END TYPE iinterval_sp

  !> string representation of iinterval_sp takes at most iinterval_sp2s_len
  !! characters
  INTEGER, PUBLIC, PARAMETER :: iinterval_sp2s_len = 2*de_g_sp_width+4


  CHARACTER(len=*), PARAMETER :: iinterval_sp_fmt = &
           '("[",'//de_g_sp//',",",'//de_g_sp//',"]")'
  CHARACTER(len=*), PARAMETER :: iinterval_dp_fmt = &
           '("[",'//de_g_dp//',",",'//de_g_dp//',"]")'
  TYPE, PUBLIC :: iinterval_dp
    SEQUENCE
    !> first delimiter of range
    REAL(dp) :: first
    !> last delimiter of range
    REAL(dp) :: last
  END TYPE iinterval_dp

  !> string representation of iinterval_dp takes at most iinterval_dp2s_len
  !! characters
  INTEGER, PUBLIC, PARAMETER :: iinterval_dp2s_len = 2*de_g_dp_width+4

  PUBLIC :: extent_size, extent_shape, extent_start, extent_end
  PUBLIC :: extent_size_i8
  PUBLIC :: rebased_extent
  PUBLIC :: extent_from_iinterval, extent_set_iinterval, char, ASSIGNMENT(=)
  PUBLIC :: iinterval_from_extent
  PUBLIC :: extent_intersect, extents_do_intersect, is_contained_in
  PUBLIC :: sprint
  PUBLIC :: init_io_formats
  PUBLIC :: OPERATOR(==), OPERATOR(/=)
  INTERFACE extent_size
    MODULE PROCEDURE extent_size_1d
    MODULE PROCEDURE extent_size_nd
    MODULE PROCEDURE extent_i8_size_1d
    MODULE PROCEDURE extent_i8_size_nd
    MODULE PROCEDURE iinterval_size_1d
    MODULE PROCEDURE iinterval_size_nd
    MODULE PROCEDURE iinterval_i8_size_1d
    MODULE PROCEDURE iinterval_i8_size_nd
  END INTERFACE
  INTERFACE extent_size_i8
    MODULE PROCEDURE extent_size_1d_i8
    MODULE PROCEDURE extent_size_nd_i8
    MODULE PROCEDURE iinterval_size_1d_i8
    MODULE PROCEDURE iinterval_size_nd_i8
  END INTERFACE
  INTERFACE extent_shape
    MODULE PROCEDURE extent_size_1d
    MODULE PROCEDURE extent_i8_size_1d
    MODULE PROCEDURE iinterval_size_1d
    MODULE PROCEDURE iinterval_i8_size_1d
  END INTERFACE extent_shape
  INTERFACE extent_start
    MODULE PROCEDURE extent_start_1d
    MODULE PROCEDURE extent_i8_start_1d
    MODULE PROCEDURE extent_start_nd
    MODULE PROCEDURE extent_i8_start_nd
    MODULE PROCEDURE iinterval_start_1d
    MODULE PROCEDURE iinterval_i8_start_1d
    MODULE PROCEDURE iinterval_start_nd
    MODULE PROCEDURE iinterval_i8_start_nd
  END INTERFACE
  INTERFACE extent_end
    MODULE PROCEDURE extent_end_1d
    MODULE PROCEDURE extent_i8_end_1d
    MODULE PROCEDURE extent_end_nd
    MODULE PROCEDURE extent_i8_end_nd
    MODULE PROCEDURE iinterval_end_1d
    MODULE PROCEDURE iinterval_i8_end_1d
    MODULE PROCEDURE iinterval_end_nd
    MODULE PROCEDURE iinterval_i8_end_nd
  END INTERFACE
  INTERFACE rebased_extent
    MODULE PROCEDURE rebased_extent_1d
    MODULE PROCEDURE rebased_extent_idx
  END INTERFACE
  INTERFACE char
    MODULE PROCEDURE char_auto_e
    MODULE PROCEDURE char_auto_e_i8
    MODULE PROCEDURE char_auto_i
    MODULE PROCEDURE char_auto_i_i8
    MODULE PROCEDURE char_auto_i_sp
    MODULE PROCEDURE char_auto_i_dp
  END INTERFACE
  INTERFACE sprint
    ELEMENTAL SUBROUTINE ppm_sprint_extent(buf, ext)
      IMPORT :: ext2s_len, extent
      CHARACTER(len=ext2s_len), INTENT(out) :: buf
      TYPE(extent), INTENT(in) :: ext
    END SUBROUTINE ppm_sprint_extent
    ELEMENTAL SUBROUTINE ppm_sprint_extent_i8(buf, ext)
      IMPORT :: ext2s_len_i8, extent_i8
      CHARACTER(len=ext2s_len_i8), INTENT(out) :: buf
      TYPE(extent_i8), INTENT(in) :: ext
    END SUBROUTINE ppm_sprint_extent_i8
    ELEMENTAL SUBROUTINE ppm_sprint_iinterval(buf, rng)
      IMPORT :: ext2s_len, iinterval
      CHARACTER(len=ext2s_len), INTENT(out) :: buf
      TYPE(iinterval), INTENT(in) :: rng
    END SUBROUTINE ppm_sprint_iinterval
    ELEMENTAL SUBROUTINE ppm_sprint_iinterval_i8(buf, rng)
      IMPORT :: ext2s_len_i8, iinterval_i8
      CHARACTER(len=ext2s_len_i8), INTENT(out) :: buf
      TYPE(iinterval_i8), INTENT(in) :: rng
    END SUBROUTINE ppm_sprint_iinterval_i8
    ELEMENTAL SUBROUTINE ppm_sprint_iinterval_sp(buf, rng)
      IMPORT :: iinterval_sp2s_len, iinterval_sp
      CHARACTER(len=iinterval_sp2s_len), INTENT(out) :: buf
      TYPE(iinterval_sp), INTENT(in) :: rng
    END SUBROUTINE ppm_sprint_iinterval_sp
    ELEMENTAL SUBROUTINE ppm_sprint_iinterval_dp(buf, rng)
      IMPORT :: iinterval_dp2s_len, iinterval_dp
      CHARACTER(len=iinterval_dp2s_len), INTENT(out) :: buf
      TYPE(iinterval_dp), INTENT(in) :: rng
    END SUBROUTINE ppm_sprint_iinterval_dp
    MODULE PROCEDURE sprint_e
    MODULE PROCEDURE sprint_e_i8
    MODULE PROCEDURE sprint_i
  END INTERFACE sprint
  INTERFACE extent_set_iinterval
    MODULE PROCEDURE extent_set_iinterval_2i
    MODULE PROCEDURE extent_set_iinterval_i4
    MODULE PROCEDURE extent_set_iinterval_i8
  END INTERFACE
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE extent_set_iinterval_i4
    MODULE PROCEDURE extent_set_iinterval_i8
    MODULE PROCEDURE iinterval_set_extent_i4
    MODULE PROCEDURE iinterval_set_extent_i8
  END INTERFACE
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE extent_equality
    MODULE PROCEDURE extent_equality_i8
    MODULE PROCEDURE iinterval_equality
    MODULE PROCEDURE iinterval_equality_i8
  END INTERFACE
  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE extent_inequality
    MODULE PROCEDURE extent_inequality_i8
    MODULE PROCEDURE iinterval_inequality
    MODULE PROCEDURE iinterval_inequality_i8
  END INTERFACE
  INTERFACE extent_intersect
    MODULE PROCEDURE extent_intersect_1d
    MODULE PROCEDURE extent_i8_intersect_1d
    MODULE PROCEDURE iinterval_intersect_1d
    MODULE PROCEDURE iinterval_i8_intersect_1d
  END INTERFACE extent_intersect
  INTERFACE extents_do_intersect
    MODULE PROCEDURE extents_do_intersect_1d
    MODULE PROCEDURE extents_do_intersect_nd
    MODULE PROCEDURE extents_i8_do_intersect_1d
    MODULE PROCEDURE extents_i8_do_intersect_nd
    MODULE PROCEDURE iintervals_do_intersect_1d
    MODULE PROCEDURE iintervals_do_intersect_nd
    MODULE PROCEDURE iintervals_i8_do_intersect_1d
    MODULE PROCEDURE iintervals_i8_do_intersect_nd
  END INTERFACE extents_do_intersect
  INTERFACE is_contained_in
    MODULE PROCEDURE is_contained_in_e
    MODULE PROCEDURE is_contained_in_e_nd
  END INTERFACE is_contained_in
CONTAINS
  !> size of one-dimensional range
  !! @param ext extent to query
  !! @return size of extent
  ELEMENTAL FUNCTION extent_size_1d(ext) RESULT(esize)
    TYPE(extent), INTENT(in) :: ext
    INTEGER(i4) :: esize
    esize = ext%size
  END FUNCTION extent_size_1d

  !> size of multi-dimensional extent,
  !! i.e. product of sizes in individual dimensions
  !! @param ext extents to query
  !! @return product of sizes of extents
  PURE FUNCTION extent_size_nd(ext) RESULT(esize)
    TYPE(extent), INTENT(in) :: ext(:)
    INTEGER(i4) :: esize, i, n
    n = SIZE(ext)
    IF (n /= 0) THEN
      esize = ext(1)%size
      DO i = 2, n
        esize = esize * ext(i)%size
      END DO
    ELSE
      esize = 0
    END IF
  END FUNCTION extent_size_nd

  !> size of one-dimensional range
  !! @param ext extent to query
  !! @return size of extent
  ELEMENTAL FUNCTION extent_i8_size_1d(ext) RESULT(esize)
    TYPE(extent_i8), INTENT(in) :: ext
    INTEGER(i8) :: esize
    esize = ext%size
  END FUNCTION extent_i8_size_1d

  !> size of multi-dimensional extent,
  !! i.e. product of sizes in individual dimensions
  !! @param ext extents to query
  !! @return product of sizes of extents
  PURE FUNCTION extent_i8_size_nd(ext) RESULT(esize)
    TYPE(extent_i8), INTENT(in) :: ext(:)
    INTEGER(i8) :: esize
    INTEGER :: i, n
    n = SIZE(ext)
    IF (n /= 0) THEN
      esize = ext(1)%size
      DO i = 2, n
        esize = esize * ext(i)%size
      END DO
    ELSE
      esize = 0_i8
    END IF
  END FUNCTION extent_i8_size_nd

  !> size of one-dimensional range
  !! @param rng iinterval to query
  !! @return size of interval
  ELEMENTAL FUNCTION iinterval_size_1d(rng) RESULT(isize)
    TYPE(iinterval), INTENT(in) :: rng
    INTEGER(i4) :: isize
    isize = rng%last - rng%first + SIGN(1_i4, rng%last - rng%first)
  END FUNCTION iinterval_size_1d

  !> size of one-dimensional range
  !! @param rng iinterval to query
  !! @return size of interval
  ELEMENTAL FUNCTION iinterval_i8_size_1d(rng) RESULT(isize)
    TYPE(iinterval_i8), INTENT(in) :: rng
    INTEGER(i8) :: isize
    isize = rng%last - rng%first + SIGN(1_i8, rng%last - rng%first)
  END FUNCTION iinterval_i8_size_1d

  !> size of multi-dimensional range,
  !! i.e. product of sizes in individual dimensions
  !! @param rng multi-dimensional range to query
  !! @return product of sizes of ranges
  PURE FUNCTION iinterval_size_nd(rng) RESULT(isize)
    TYPE(iinterval), INTENT(in) :: rng(:)
    INTEGER(i4) :: isize
    INTEGER :: i, n
    n = SIZE(rng)
    IF (n /= 0) THEN
      isize = iinterval_size_1d(rng(1))
      DO i = 2, n
        isize = isize * iinterval_size_1d(rng(i))
      END DO
    ELSE
      isize = 0
    END IF
  END FUNCTION iinterval_size_nd

  !> size of multi-dimensional range,
  !! i.e. product of sizes in individual dimensions
  !! @param rng multi-dimensional range to query
  !! @return product of sizes of ranges
  PURE FUNCTION iinterval_i8_size_nd(rng) RESULT(isize)
    TYPE(iinterval_i8), INTENT(in) :: rng(:)
    INTEGER(i8) :: isize
    INTEGER :: i, n
    n = SIZE(rng)
    IF (n /= 0) THEN
      isize = iinterval_i8_size_1d(rng(1))
      DO i = 2, n
        isize = isize * iinterval_i8_size_1d(rng(i))
      END DO
    ELSE
      isize = 0_i8
    END IF
  END FUNCTION iinterval_i8_size_nd

  !> size of one-dimensional range
  !! @param ext extent to query
  !! @return size of extent
  ELEMENTAL FUNCTION extent_size_1d_i8(ext) RESULT(extent_size)
    TYPE(extent), INTENT(in) :: ext
    INTEGER(i8) :: extent_size
    extent_size = INT(ext%size, i8)
  END FUNCTION extent_size_1d_i8

  !> size of multi-dimensional extent,
  !! i.e. product of sizes in individual dimensions
  !! @param ext extents to query
  !! @return product of sizes of extents
  PURE FUNCTION extent_size_nd_i8(ext) RESULT(esize)
    TYPE(extent), INTENT(in) :: ext(:)
    INTEGER(i8) :: esize
    INTEGER :: i, n
    n = SIZE(ext)
    IF (n /= 0) THEN
      esize = INT(ext(1)%size, i8)
      DO i = 2, n
        esize = esize * INT(ext(i)%size, i8)
      END DO
    ELSE
      esize = 0_i8
    END IF
  END FUNCTION extent_size_nd_i8

  !> size of one-dimensional range
  !! @param rng iinterval to query
  !! @return size of interval
  ELEMENTAL FUNCTION iinterval_size_1d_i8(rng) RESULT(isize)
    TYPE(iinterval), INTENT(in) :: rng
    INTEGER(i8) :: isize
    isize = INT(rng%last - rng%first + SIGN(1, rng%last - rng%first), i8)
  END FUNCTION iinterval_size_1d_i8

  !> size of multi-dimensional range,
  !! i.e. product of sizes in individual dimensions
  !! @param rng multi-dimensional range to query
  !! @return product of sizes of ranges
  PURE FUNCTION iinterval_size_nd_i8(rng) RESULT(isize)
    TYPE(iinterval), INTENT(in) :: rng(:)
    INTEGER(i8) :: isize
    INTEGER(i4) :: i, n
    n = SIZE(rng)
    IF (n /= 0) THEN
      isize = INT(iinterval_size_1d(rng(1)), i8)
      DO i = 2, n
        isize = isize * iinterval_size_1d_i8(rng(i))
      END DO
    ELSE
      isize = 0_i8
    END IF
  END FUNCTION iinterval_size_nd_i8

  !> start index of one-dimensional range
  !! @param ext extent to query
  !! @return first index of interval
  ELEMENTAL FUNCTION extent_start_1d(ext) RESULT(start_idx)
    TYPE(extent), INTENT(in) :: ext
    INTEGER(i4) :: start_idx
    start_idx = ext%first
  END FUNCTION extent_start_1d

  !> start index of one-dimensional range
  !! @param ext extent to query
  !! @return first index of interval
  ELEMENTAL FUNCTION extent_i8_start_1d(ext) RESULT(start_idx)
    TYPE(extent_i8), INTENT(in) :: ext
    INTEGER(i8) :: start_idx
    start_idx = ext%first
  END FUNCTION extent_i8_start_1d

  !> start indices of multi-dimensional ranges
  !! @param ext extents to query
  !! @return first indices of interval
  PURE FUNCTION extent_start_nd(ext) RESULT(starts)
    TYPE(extent), INTENT(in) :: ext(:)
    INTEGER(i4) :: starts(SIZE(ext))
    starts = ext%first
  END FUNCTION extent_start_nd

  !> start indices of multi-dimensional ranges
  !! @param ext extents to query
  !! @return first indices of interval
  PURE FUNCTION extent_i8_start_nd(ext) RESULT(starts)
    TYPE(extent_i8), INTENT(in) :: ext(:)
    INTEGER(i8) :: starts(SIZE(ext))
    starts = ext%first
  END FUNCTION extent_i8_start_nd

  !> start index of one-dimensional range
  !! @param rng iinterval to query
  !! @return first index of interval
  ELEMENTAL FUNCTION iinterval_start_1d(rng) RESULT(start_idx)
    TYPE(iinterval), INTENT(in) :: rng
    INTEGER(i4) :: start_idx
    start_idx = rng%first
  END FUNCTION iinterval_start_1d

  !> start index of one-dimensional range
  !! @param rng iinterval to query
  !! @return first index of interval
  ELEMENTAL FUNCTION iinterval_i8_start_1d(rng) RESULT(start_idx)
    TYPE(iinterval_i8), INTENT(in) :: rng
    INTEGER(i8) :: start_idx
    start_idx = rng%first
  END FUNCTION iinterval_i8_start_1d

  !> start indices of multi-dimensional ranges
  !! @param ranges iintervals to query
  !! @return first indices of intervals
  PURE FUNCTION iinterval_start_nd(ranges) RESULT(starts)
    TYPE(iinterval), INTENT(in) :: ranges(:)
    INTEGER(i4) :: starts(SIZE(ranges))
    starts(:) = ranges(:)%first
  END FUNCTION iinterval_start_nd

  !> start indices of multi-dimensional ranges
  !! @param ranges iintervals to query
  !! @return first indices of intervals
  PURE FUNCTION iinterval_i8_start_nd(ranges) RESULT(starts)
    TYPE(iinterval_i8), INTENT(in) :: ranges(:)
    INTEGER(i8) :: starts(SIZE(ranges))
    starts = ranges%first
  END FUNCTION iinterval_i8_start_nd

  ELEMENTAL FUNCTION extent_end_1d(ext) RESULT(end_idx)
    TYPE(extent), INTENT(in) :: ext
    INTEGER(i4) :: end_idx
    end_idx = ext%first + ext%size - SIGN(1, ext%size)
  END FUNCTION extent_end_1d

  ELEMENTAL FUNCTION extent_i8_end_1d(ext) RESULT(end_idx)
    TYPE(extent_i8), INTENT(in) :: ext
    INTEGER(i8) :: end_idx
    end_idx = ext%first + ext%size - SIGN(1_i8, ext%size)
  END FUNCTION extent_i8_end_1d

  PURE FUNCTION extent_end_nd(ext) RESULT(ext_end)
    TYPE(extent), INTENT(in) :: ext(:)
    INTEGER(i4) :: ext_end(SIZE(ext))
    ext_end = extent_end_1d(ext)
  END FUNCTION extent_end_nd

  PURE FUNCTION extent_i8_end_nd(ext) RESULT(ext_end)
    TYPE(extent_i8), INTENT(in) :: ext(:)
    INTEGER(i8) :: ext_end(SIZE(ext))
    ext_end = extent_i8_end_1d(ext)
  END FUNCTION extent_i8_end_nd

  ELEMENTAL FUNCTION iinterval_end_1d(rng) RESULT(r_end)
    TYPE(iinterval), INTENT(in) :: rng
    INTEGER(i4) :: r_end
    r_end = rng%last
  END FUNCTION iinterval_end_1d

  ELEMENTAL FUNCTION iinterval_i8_end_1d(rng) RESULT(r_end)
    TYPE(iinterval_i8), INTENT(in) :: rng
    INTEGER(i8) :: r_end
    r_end = rng%last
  END FUNCTION iinterval_i8_end_1d

  PURE FUNCTION iinterval_end_nd(ranges) RESULT(r_ends)
    TYPE(iinterval), INTENT(in) :: ranges(:)
    INTEGER(i4) :: r_ends(SIZE(ranges))
    r_ends = iinterval_end_1d(ranges)
  END FUNCTION iinterval_end_nd

  PURE FUNCTION iinterval_i8_end_nd(ranges) RESULT(r_ends)
    TYPE(iinterval_i8), INTENT(in) :: ranges(:)
    INTEGER(i8) :: r_ends(SIZE(ranges))
    r_ends = iinterval_i8_end_1d(ranges)
  END FUNCTION iinterval_i8_end_nd

  ELEMENTAL FUNCTION extent_iinterval(ext)
    TYPE(extent), INTENT(in) :: ext
    TYPE(iinterval) :: extent_iinterval
    extent_iinterval%first = extent_start(ext)
    extent_iinterval%last = extent_end(ext)
  END FUNCTION extent_iinterval

  ELEMENTAL FUNCTION rebased_extent_1d(ext, base) RESULT(norm)
    TYPE(extent) :: norm
    TYPE(extent), INTENT(in) :: ext, base
    norm = extent(ext%first - base%first, ext%size)
  END FUNCTION rebased_extent_1d

  ELEMENTAL FUNCTION rebased_extent_idx(ext, base) RESULT(norm)
    TYPE(extent) :: norm
    TYPE(extent), INTENT(in) :: ext
    INTEGER(i4), INTENT(in) :: base
    norm = extent(ext%first - base, ext%size)
  END FUNCTION rebased_extent_idx

  SUBROUTINE extent_set_iinterval_2i(ext, first, last)
    TYPE(extent), INTENT(out) :: ext
    INTEGER(i4), INTENT(in) :: first, last
    ext%first = first
    ext%size = last - first + SIGN(1, last - first)
  END SUBROUTINE extent_set_iinterval_2i

  ELEMENTAL SUBROUTINE extent_set_iinterval_i4(ext, rng)
    TYPE(extent), INTENT(out) :: ext
    TYPE(iinterval), INTENT(in) :: rng
    ext%first = rng%first
    ext%size = rng%last - rng%first + SIGN(1, rng%last - rng%first)
  END SUBROUTINE extent_set_iinterval_i4

  ELEMENTAL SUBROUTINE extent_set_iinterval_i8(ext, rng)
    TYPE(extent_i8), INTENT(out) :: ext
    TYPE(iinterval_i8), INTENT(in) :: rng
    ext%first = rng%first
    ext%size = rng%last - rng%first + SIGN(1_i8, rng%last - rng%first)
  END SUBROUTINE extent_set_iinterval_i8

  !> assign extent to range,
  !! caution: this is only well if defined either last < first denotes
  !! the empty range or zero-size extents are not part of the input
  ELEMENTAL SUBROUTINE iinterval_set_extent_i4(rng, ext)
    TYPE(iinterval), INTENT(out) :: rng
    TYPE(extent), INTENT(in) :: ext
    rng%first = extent_start(ext); rng%last = extent_end(ext)
  END SUBROUTINE iinterval_set_extent_i4

  !> assign extent to range,
  !! caution: this is only well if defined either last < first denotes
  !! the empty range or zero-size extents are not part of the input
  ELEMENTAL SUBROUTINE iinterval_set_extent_i8(rng, ext)
    TYPE(iinterval_i8), INTENT(out) :: rng
    TYPE(extent_i8), INTENT(in) :: ext
    rng%first = extent_start(ext); rng%last = extent_end(ext)
  END SUBROUTINE iinterval_set_extent_i8

  !> construct range from extent,
  !! caution: this is only well if defined either last < first denotes
  !! the empty range or zero-size extents are not part of the input
  ELEMENTAL FUNCTION iinterval_from_extent(ext) RESULT(rng)
    TYPE(iinterval) :: rng
    TYPE(extent), INTENT(in) :: ext
    rng = ext
  END FUNCTION iinterval_from_extent

  !> convert inclusive interval to extent representing same range
  ELEMENTAL FUNCTION extent_from_iinterval(first, last) RESULT(ext)
    TYPE(extent) :: ext
    INTEGER(i4), INTENT(in) :: first, last
    ext%first = first
    ext%size = last - first + SIGN(1, last - first)
  END FUNCTION extent_from_iinterval

  ELEMENTAL FUNCTION char_auto_e(ext) RESULT(str)
    CHARACTER(len=ext2s_len) :: str
    TYPE(extent), INTENT(in) :: ext
    IF (ext%size /= 0) THEN



      WRITE (str, '(a,i0,a,i0,a)') '[', extent_start(ext), ',', &
           extent_end(ext), ']'

    ELSE
      str = '{}'
    END IF
  END FUNCTION char_auto_e

  ELEMENTAL FUNCTION char_auto_e_i8(ext) RESULT(str)
    CHARACTER(len=ext2s_len_i8) :: str
    TYPE(extent_i8), INTENT(in) :: ext
    IF (ext%size /= 0_i8) THEN



      WRITE (str, '(a,i0,a,i0,a)') '[', extent_start(ext), ',', &
           extent_end(ext), ']'

    ELSE
      str = '{}'
    END IF
  END FUNCTION char_auto_e_i8

  SUBROUTINE sprint_e(s, ranges)
    CHARACTER(*), INTENT(out) :: s
    TYPE(extent), INTENT(in) :: ranges(:)

    CHARACTER(len=ext2s_len) :: estr
    INTEGER :: i, n, slen

    slen = LEN(s)
    n = 1
    IF (SIZE(ranges) > 0) THEN
      estr = CHAR(ranges(1))
      n = LEN_TRIM(estr)
      s(1:MIN(n, slen)) = TRIM(estr)
      n = n + 1
    END IF
    DO i = 2, SIZE(ranges)
      estr = CHAR(ranges(i))
      IF (n <= slen) s(n:slen) = ', ' // estr
      n = n + LEN_TRIM(estr) + 2
      IF (n > slen) EXIT
    END DO
    IF (n <= slen) s(n:slen) = ' '
  END SUBROUTINE sprint_e

  SUBROUTINE sprint_e_i8(s, ranges)
    CHARACTER(*), INTENT(out) :: s
    TYPE(extent_i8), INTENT(in) :: ranges(:)

    CHARACTER(len=ext2s_len_i8) :: estr
    INTEGER :: i, n, slen

    slen = LEN(s)
    n = 1
    IF (SIZE(ranges) > 0) THEN
      estr = CHAR(ranges(1))
      n = LEN_TRIM(estr)
      s(1:MIN(n, slen)) = TRIM(estr)
      n = n + 1
    END IF
    DO i = 2, SIZE(ranges)
      estr = CHAR(ranges(i))
      IF (n <= slen) s(n:slen) = ', ' // estr
      n = n + LEN_TRIM(estr) + 2
      IF (n > slen) EXIT
    END DO
    IF (n <= slen) s(n:slen) = ' '
  END SUBROUTINE sprint_e_i8

  SUBROUTINE sprint_i(s, ranges)
    CHARACTER(*), INTENT(out) :: s
    TYPE(iinterval), INTENT(in) :: ranges(:)

    CHARACTER(len=ext2s_len) :: estr
    INTEGER :: i, n, slen

    slen = LEN(s)
    n = 1
    IF (SIZE(ranges) > 0) THEN
      estr = CHAR(ranges(1))
      n = LEN_TRIM(estr)
      s(1:MIN(n, slen)) = TRIM(estr)
      n = n + 1
    END IF
    DO i = 2, SIZE(ranges)
      estr = CHAR(ranges(i))
      IF (n <= slen) s(n:slen) = ', ' // estr
      n = n + LEN_TRIM(estr) + 2
      IF (n > slen) EXIT
    END DO
    IF (n <= slen) s(n:slen) = ' '
  END SUBROUTINE sprint_i

  ELEMENTAL FUNCTION char_auto_i(rng) RESULT(str)
    CHARACTER(len=ext2s_len) :: str
    TYPE(iinterval), INTENT(in) :: rng



    WRITE (str, '(a,i0,a,i0,a)') '[', rng%first, ',', &
         rng%last, ']'

  END FUNCTION char_auto_i

  ELEMENTAL FUNCTION char_auto_i_i8(rng) RESULT(str)
    CHARACTER(len=ext2s_len_i8) :: str
    TYPE(iinterval_i8), INTENT(in) :: rng



    WRITE (str, '(a,i0,a,i0,a)') '[', rng%first, ',', &
         rng%last, ']'

  END FUNCTION char_auto_i_i8

  SUBROUTINE init_io_formats






  END SUBROUTINE init_io_formats

  ELEMENTAL FUNCTION char_auto_i_sp(rng) RESULT(str)
    CHARACTER(len=iinterval_sp2s_len) :: str
    TYPE(iinterval_sp), INTENT(in) :: rng



    WRITE (str, fmt=iinterval_sp_fmt) rng%first, rng%last

  END FUNCTION char_auto_i_sp

  ELEMENTAL FUNCTION char_auto_i_dp(rng) RESULT(str)
    CHARACTER(len=iinterval_dp2s_len) :: str
    TYPE(iinterval_dp), INTENT(in) :: rng



    WRITE (str, fmt=iinterval_dp_fmt) rng%first, rng%last

  END FUNCTION char_auto_i_dp

  ELEMENTAL FUNCTION extent_equality(a, b) RESULT(l)
    TYPE(extent), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first == b%first &
         .AND. (a%size == b%size .OR. ABS(a%size) == 1 .AND. ABS(b%size) == 1)
  END FUNCTION extent_equality

  ELEMENTAL FUNCTION extent_inequality(a, b) RESULT(l)
    TYPE(extent), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first /= b%first .OR. a%size /= b%size
  END FUNCTION extent_inequality

  ELEMENTAL FUNCTION extent_equality_i8(a, b) RESULT(l)
    TYPE(extent_i8), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first == b%first &
         .AND. (a%size == b%size &
         &      .OR. ABS(a%size) == 1_i8 .AND. ABS(b%size) == 1_i8)
  END FUNCTION extent_equality_i8

  ELEMENTAL FUNCTION extent_inequality_i8(a, b) RESULT(l)
    TYPE(extent_i8), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first /= b%first .OR. a%size /= b%size
  END FUNCTION extent_inequality_i8

  ELEMENTAL FUNCTION iinterval_equality(a, b) RESULT(l)
    TYPE(iinterval), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first == b%first .AND. a%last == b%last
  END FUNCTION iinterval_equality

  ELEMENTAL FUNCTION iinterval_inequality(a, b) RESULT(l)
    TYPE(iinterval), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first /= b%first .OR. a%last /= b%last
  END FUNCTION iinterval_inequality

  ELEMENTAL FUNCTION iinterval_equality_i8(a, b) RESULT(l)
    TYPE(iinterval_i8), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first == b%first .AND. a%last == b%last
  END FUNCTION iinterval_equality_i8

  ELEMENTAL FUNCTION iinterval_inequality_i8(a, b) RESULT(l)
    TYPE(iinterval_i8), INTENT(in) :: a, b
    LOGICAL :: l
    l = a%first /= b%first .OR. a%last /= b%last
  END FUNCTION iinterval_inequality_i8

  ELEMENTAL FUNCTION extent_intersect_1d(a, b) RESULT(intersection)
    TYPE(extent) :: intersection
    TYPE(extent), INTENT(in) :: a, b
    TYPE(iinterval) :: a_i, b_i, temp_i
    IF (a%size == 0 .OR. b%size == 0) THEN
      intersection = empty_extent
    ELSE
      a_i = MERGE(iinterval(a%first, a%first + a%size - 1), &
           &      iinterval(a%first + a%size + 1, a%first), a%size > 0)
      b_i = MERGE(iinterval(b%first, b%first + b%size - 1), &
           &      iinterval(b%first + b%size + 1, b%first), b%size > 0)
      IF (a_i%first > b_i%first) THEN
        temp_i = a_i
        a_i = b_i
        b_i = temp_i
      END IF
      IF (b_i%first > a_i%last) THEN
        intersection = empty_extent
      ELSE
        intersection = extent(b_i%first, MIN(a_i%last - b_i%first + 1, &
             b_i%last - b_i%first + 1))
      END IF
    END IF
  END FUNCTION extent_intersect_1d

  ELEMENTAL FUNCTION extent_i8_intersect_1d(a, b) RESULT(intersection)
    TYPE(extent_i8) :: intersection
    TYPE(extent_i8), INTENT(in) :: a, b
    TYPE(iinterval_i8) :: a_i, b_i, temp_i
    IF (a%size == 0_i8 .OR. b%size == 0_i8) THEN
      intersection = empty_extent_i8
    ELSE
      a_i = MERGE(iinterval_i8(a%first, a%first + a%size - 1_i8), &
           &      iinterval_i8(a%first + a%size + 1_i8, a%first), a%size > 0_i8)
      b_i = MERGE(iinterval_i8(b%first, b%first + b%size - 1_i8), &
           &      iinterval_i8(b%first + b%size + 1_i8, b%first), b%size > 0_i8)
      IF (a_i%first > b_i%first) THEN
        temp_i = a_i
        a_i = b_i
        b_i = temp_i
      END IF
      IF (b_i%first > a_i%last) THEN
        intersection = empty_extent_i8
      ELSE
        intersection = extent_i8(b_i%first, MIN(a_i%last - b_i%first + 1_i8, &
             b_i%last - b_i%first + 1_i8))
      END IF
    END IF
  END FUNCTION extent_i8_intersect_1d

  ELEMENTAL FUNCTION iinterval_intersect_1d(a, b) RESULT(intersection)
    TYPE(iinterval) :: intersection
    TYPE(iinterval), INTENT(in) :: a, b
    IF ((a%first > a%last .OR. b%first > b%last) &
         .OR. (a%first < b%first .AND. b%first > a%last) &
         .OR. (a%first > b%first .AND. a%first > b%last)) THEN
      intersection = empty_iinterval
    ELSE IF (a%first > b%first) THEN
      intersection = iinterval(a%first, a%first + MIN(b%last - a%first, &
           a%last - a%first))
    ELSE
      intersection = iinterval(b%first, b%first + MIN(a%last - b%first, &
           b%last - b%first))
    END IF
  END FUNCTION iinterval_intersect_1d

  ELEMENTAL FUNCTION iinterval_i8_intersect_1d(a, b) RESULT(intersection)
    TYPE(iinterval_i8) :: intersection
    TYPE(iinterval_i8), INTENT(in) :: a, b
    IF ((a%first > a%last .OR. b%first > b%last) &
         .OR. (a%first < b%first .AND. b%first > a%last) &
         .OR. (a%first > b%first .AND. a%first > b%last)) THEN
      intersection = empty_iinterval_i8
    ELSE IF (a%first > b%first) THEN
      intersection = iinterval_i8(a%first, a%first + MIN(b%last - a%first, &
           a%last - a%first))
    ELSE
      intersection = iinterval_i8(b%first, b%first + MIN(a%last - b%first, &
           b%last - b%first))
    END IF
  END FUNCTION iinterval_i8_intersect_1d

  PURE FUNCTION extents_do_intersect_1d(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(extent), INTENT(in) :: a, b
    p = extent_size(extent_intersect(a, b)) /= 0
  END FUNCTION extents_do_intersect_1d

  PURE FUNCTION extents_do_intersect_nd(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(extent), INTENT(in) :: a(:), b(:)
    INTEGER :: i, n
    n = SIZE(a)
    p = n >= 1
    DO i = 1, n
      p = p .AND. extent_size(extent_intersect(a(i), b(i))) /= 0
    END DO
  END FUNCTION extents_do_intersect_nd

  PURE FUNCTION extents_i8_do_intersect_1d(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(extent_i8), INTENT(in) :: a, b
    p = extent_size(extent_intersect(a, b)) /= 0_i8
  END FUNCTION extents_i8_do_intersect_1d

  PURE FUNCTION extents_i8_do_intersect_nd(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(extent_i8), INTENT(in) :: a(:), b(:)
    INTEGER :: i, n
    n = SIZE(a)
    p = n >= 1
    DO i = 1, n
      p = p .AND. extent_size(extent_intersect(a(i), b(i))) /= 0
    END DO
  END FUNCTION extents_i8_do_intersect_nd

  PURE FUNCTION iintervals_do_intersect_1d(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(iinterval), INTENT(in) :: a, b
    p = extent_size(extent_intersect(a, b)) /= 0
  END FUNCTION iintervals_do_intersect_1d

  PURE FUNCTION iintervals_do_intersect_nd(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(iinterval), INTENT(in) :: a(:), b(:)
    INTEGER :: i, n
    n = SIZE(a)
    p = n >= 1
    DO i = 1, n
      p = p .AND. extent_size(extent_intersect(a(i), b(i))) /= 0
    END DO
  END FUNCTION iintervals_do_intersect_nd

  PURE FUNCTION iintervals_i8_do_intersect_1d(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(iinterval_i8), INTENT(in) :: a, b
    p = extent_size(extent_intersect(a, b)) /= 0_i8
  END FUNCTION iintervals_i8_do_intersect_1d

  PURE FUNCTION iintervals_i8_do_intersect_nd(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(iinterval_i8), INTENT(in) :: a(:), b(:)
    INTEGER :: i, n
    n = SIZE(a)
    p = n >= 1
    DO i = 1, n
      p = p .AND. extent_size(extent_intersect(a(i), b(i))) /= 0_i8
    END DO
  END FUNCTION iintervals_i8_do_intersect_nd

  ELEMENTAL FUNCTION is_contained_in_e(i, rng) RESULT(p)
    INTEGER(i4), INTENT(in) :: i
    TYPE(extent), INTENT(in) :: rng
    LOGICAL :: p
    p = i >= rng%first .AND. i < rng%first + rng%size
  END FUNCTION is_contained_in_e

  FUNCTION is_contained_in_e_nd(i, rng) RESULT(p)
    INTEGER(i4), INTENT(in) :: i(:)
    TYPE(extent), INTENT(in) :: rng(:)
    LOGICAL :: p
    p = ALL(i >= rng%first .AND. i < rng%first + rng%size)
  END FUNCTION is_contained_in_e_nd

END MODULE ppm_extents
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
