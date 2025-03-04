!
! Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!
! Version: 1.0
! Keywords: NetCDF
! Author: Thomas Jahns <jahns@dkrz.de>
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
!! @file ppm_ncdf_dump.f90
!> write a single array to a NetCDF file, for debugging mostly

MODULE ppm_ncdf_dump
  USE ppm_std_type_kinds, ONLY: dp, i4, sp
  USE ppm_base, ONLY: assertion, abort_ppm
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  PRIVATE
  PUBLIC :: dump_ncdf_single_array

  INTERFACE dump_ncdf_single_array
    MODULE PROCEDURE dump_ncdf_single_array_sp_asize
    MODULE PROCEDURE dump_ncdf_single_array_sp_1d
    MODULE PROCEDURE dump_ncdf_single_array_sp_2d
    MODULE PROCEDURE dump_ncdf_single_array_sp_3d
    MODULE PROCEDURE dump_ncdf_single_array_dp_asize
    MODULE PROCEDURE dump_ncdf_single_array_dp_1d
    MODULE PROCEDURE dump_ncdf_single_array_dp_2d
    MODULE PROCEDURE dump_ncdf_single_array_dp_3d
    MODULE PROCEDURE dump_ncdf_single_array_i4_asize
    MODULE PROCEDURE dump_ncdf_single_array_i4_1d
    MODULE PROCEDURE dump_ncdf_single_array_i4_2d
    MODULE PROCEDURE dump_ncdf_single_array_i4_3d
  END INTERFACE

  CHARACTER(*), PARAMETER :: longitude = 'longitude', latitude = 'latitude', &
       degrees_east = 'degrees_east', degrees_north = 'degrees_north'

  CHARACTER(len=*), PARAMETER :: filename = 'ppm_ncdf_dump.f90'
CONTAINS
  !> write array to file in NetCDF format
  !! @param dump_fname file name to use
  !! @param a array of data items to write
  !! @param a_shape a is shaped as given in a_shape (i.e. a_shape=shape(a))
  !! @param a_name metadata name for a
  !! @param a_valid_range valid entries in a are assumed to be within
  !! [a_valid_range(1),a_valid_range(2)]
  !! @param a_fill_value entries of a with this value are assumed undefined
  !! @param title description of data dump
  !! @param source description source of the data
  !! @param institution name of originating institution
  !! @param a_long_name long display name (is not subject to same
  !! conventions as a_name and may contain e.g. spaces)
  !! @param alon longitude array of same shape as a, if given, a is
  !! assumed to contain data on planetary grid, where each entry
  !! contains the longitude of the corresponding entry of a
  !! @param alat same as alon for latitudes
  SUBROUTINE dump_ncdf_single_array_sp_asize(dump_fname, a, ndims, a_shape, &
       a_name, a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    INTEGER, INTENT(in) :: ndims
    REAL(sp), INTENT(in) :: a(*)
    INTEGER, INTENT(in) :: a_shape(ndims)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(sp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER :: ncid, varid_a, varid_lat, varid_lon
    INTEGER :: a_dimids(ndims)

    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_real, ndims, &
         a_dimids, varid_a), 133)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_sp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 140)

    CALL handle_ncdf_err(nf_put_var_real(ncid, varid_a, &
         a(1:PRODUCT(a_shape))), 143)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 148)
  END SUBROUTINE dump_ncdf_single_array_sp_asize

  SUBROUTINE dump_ncdf_single_array_sp_1d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(sp), INTENT(in) :: a(:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(sp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(1), a_shape(1)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_real, 1, &
         a_dimids, varid_a), 172)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_sp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 179)

    CALL handle_ncdf_err(nf_put_var_real(ncid, varid_a, a), 181)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 186)
  END SUBROUTINE dump_ncdf_single_array_sp_1d

  SUBROUTINE dump_ncdf_single_array_sp_2d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(sp), INTENT(in) :: a(:,:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(sp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 2
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(2), a_shape(2)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_real, 2, &
         a_dimids, varid_a), 210)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_sp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 217)

    CALL handle_ncdf_err(nf_put_var_real(ncid, varid_a, a), 219)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 224)
  END SUBROUTINE dump_ncdf_single_array_sp_2d

  SUBROUTINE dump_ncdf_single_array_sp_3d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(sp), INTENT(in) :: a(:,:,:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(sp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(3), a_shape(3)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_real, 3, &
         a_dimids, varid_a), 248)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_sp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 255)

    CALL handle_ncdf_err(nf_put_var_real(ncid, varid_a, a), 257)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 262)
  END SUBROUTINE dump_ncdf_single_array_sp_3d

  SUBROUTINE dump_ncdf_single_array_dp_asize(dump_fname, a, ndims, a_shape, &
       a_name, a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(dp), INTENT(in) :: a(*)
    INTEGER, INTENT(in) :: ndims
    INTEGER, INTENT(in) :: a_shape(ndims)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(dp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER :: ncid, varid_a, varid_lat, varid_lon
    INTEGER :: a_dimids(ndims)

    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_double, ndims, &
         a_dimids, varid_a), 287)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_dp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 294)

    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_a, &
         a(1:PRODUCT(a_shape))), 297)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 302)
  END SUBROUTINE dump_ncdf_single_array_dp_asize

  SUBROUTINE dump_ncdf_single_array_dp_1d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(dp), INTENT(in) :: a(:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(dp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(1), a_shape(1)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_double, 1, &
         a_dimids, varid_a), 326)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_dp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 333)

    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_a, a), 335)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 340)
  END SUBROUTINE dump_ncdf_single_array_dp_1d

  SUBROUTINE dump_ncdf_single_array_dp_2d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(dp), INTENT(in) :: a(:,:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(dp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 2
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(2), a_shape(2)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_double, 2, &
         a_dimids, varid_a), 364)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_dp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 371)

    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_a, a), 373)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 378)
  END SUBROUTINE dump_ncdf_single_array_dp_2d

  SUBROUTINE dump_ncdf_single_array_dp_3d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(dp), INTENT(in) :: a(:,:,:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(dp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(3), a_shape(3)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_double, 3, &
         a_dimids, varid_a), 402)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_dp(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 409)

    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_a, a), 411)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 416)
  END SUBROUTINE dump_ncdf_single_array_dp_3d

  SUBROUTINE dump_ncdf_single_array_i4_asize(dump_fname, a, ndims, a_shape, &
       a_name, a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    INTEGER(i4), INTENT(in) :: a(*)
    INTEGER, INTENT(in) :: ndims
    INTEGER, INTENT(in) :: a_shape(ndims)
    CHARACTER(len=*), INTENT(in) :: a_name
    INTEGER(i4), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER :: ncid, varid_a, varid_lat, varid_lon
    INTEGER :: a_dimids(ndims)

    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_int, ndims, &
         a_dimids, varid_a), 441)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_i4(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 448)

    CALL handle_ncdf_err(nf_put_var_int(ncid, varid_a, a(1:PRODUCT(a_shape))), &
         451)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 456)
  END SUBROUTINE dump_ncdf_single_array_i4_asize

  SUBROUTINE dump_ncdf_single_array_i4_1d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    INTEGER(i4), INTENT(in) :: a(:)
    CHARACTER(len=*), INTENT(in) :: a_name
    INTEGER(i4), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 1
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(1), a_shape(1)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_int, 1, &
         a_dimids, varid_a), 480)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_i4(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 487)

    CALL handle_ncdf_err(nf_put_var_int(ncid, varid_a, a), &
         490)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 495)
  END SUBROUTINE dump_ncdf_single_array_i4_1d

  SUBROUTINE dump_ncdf_single_array_i4_2d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    INTEGER(i4), INTENT(in) :: a(:,:)
    CHARACTER(len=*), INTENT(in) :: a_name
    INTEGER(i4), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 2
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(2), a_shape(2)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_int, 2, &
         a_dimids, varid_a), 519)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_i4(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 526)

    CALL handle_ncdf_err(nf_put_var_int(ncid, varid_a, a), &
         529)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 534)
  END SUBROUTINE dump_ncdf_single_array_i4_2d

  SUBROUTINE dump_ncdf_single_array_i4_3d(dump_fname, a, a_name, &
       a_valid_range, a_fill_value, &
       title, source, institution, a_long_name, alon, alat)
    CHARACTER(len=*), intent(in) :: dump_fname
    INTEGER(i4), TARGET, INTENT(in) :: a(:,:,:)
    CHARACTER(len=*), INTENT(in) :: a_name
    INTEGER(i4), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alon(:,:), alat(:,:)

    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: ncid, varid_a, varid_lat, varid_lon, a_dimids(3), a_shape(3)

    a_shape = SHAPE(a)
    CALL cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
         title, source, institution)

    CALL handle_ncdf_err(nf_def_var(ncid, a_name, nf_int, 3, &
         a_dimids, varid_a), 558)

    CALL put_cf_array_attributes_common(ncid, varid_a, &
         ndims, a_dimids, a_shape, &
         varid_lat, varid_lon, alon, alat, a_long_name)
    CALL put_cf_array_attributes_i4(ncid, varid_a, a_valid_range, a_fill_value)

    CALL handle_ncdf_err(nf_enddef(ncid), 565)

    CALL handle_ncdf_err(nf_put_var_int(ncid, varid_a, a), &
         568)

    IF (PRESENT(alon) .AND. PRESENT(alat)) &
         CALL put_cf_geometry(ncid, varid_lat, varid_lon, alon, alat)

    CALL handle_ncdf_err(nf_close(ncid), 573)
  END SUBROUTINE dump_ncdf_single_array_i4_3d

  SUBROUTINE cf_common_setup(dump_fname, ncid, a_shape, a_dimids, &
       title, source, institution)
    CHARACTER(len=*), INTENT(in) :: dump_fname
    INTEGER, INTENT(inout) :: ncid
    INTEGER, INTENT(in) :: a_shape(:)
    INTEGER, INTENT(inout) :: a_dimids(:)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: title
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: institution

    CHARACTER(len=3), PARAMETER :: dim_name(7) = (/ 'x  ', 'y  ', 'z  ', &
         'xx ', 'yy ', 'zz ', 'xxx' /)
    INTEGER :: ndims, i

    ndims = SIZE(a_shape)
    CALL handle_ncdf_err(nf_create(dump_fname, nf_clobber, ncid), &
         592)
    DO i = 1, ndims
      CALL handle_ncdf_err(nf_def_dim(ncid, TRIM(dim_name(i)), a_shape(i), &
           a_dimids(i)), 595)
    END DO
    ! write global cf information
    IF (PRESENT(title)) &
         CALL handle_ncdf_err(nf_put_att_text(ncid, nf_global, &
         'title', LEN(title), title), 600)
    IF (PRESENT(institution)) &
         CALL handle_ncdf_err(nf_put_att_text(ncid, nf_global, &
         'institution', LEN(institution), institution), 603)
    IF (PRESENT(source)) &
         CALL handle_ncdf_err(nf_put_att_text(ncid, nf_global, &
         'source', LEN(source), source), 606)
  END SUBROUTINE cf_common_setup

  SUBROUTINE put_cf_array_attributes_common(ncid, varid_a, &
       ndims, a_dimids, a_shape, &
       varid_lat, varid_lon, alat, alon, a_long_name)
    INTEGER, INTENT(in) :: ndims
    INTEGER, INTENT(in) :: ncid, varid_a, a_shape(ndims)
    INTEGER, INTENT(inout) :: a_dimids(ndims)
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: a_long_name
    REAL(dp), OPTIONAL, INTENT(in) :: alat(:,:), alon(:,:)
    INTEGER, INTENT(inout) :: varid_lat, varid_lon

    IF (PRESENT(a_long_name)) &
         CALL handle_ncdf_err(nf_put_att_text(ncid, varid_a, &
         'long_name', LEN(a_long_name), a_long_name), 621)
    IF (PRESENT(alat) .AND. PRESENT(alon)) THEN
      CALL assertion(ndims >= 2, filename, 623)
      CALL assertion(SIZE(alon, 1) == SIZE(alat, 1) &
           .AND. SIZE(alon, 2) == SIZE(alat, 2) &
           .AND. SIZE(alon, 1) == a_shape(1) &
           .AND. SIZE(alon, 2) == a_shape(2), filename, 627)
      CALL cf_latlon_geometry_setup(ncid, varid_a, a_shape(1:2), &
           a_dimids(1:2), varid_lat, varid_lon)
    END IF
  END SUBROUTINE put_cf_array_attributes_common

  SUBROUTINE put_cf_array_attributes_sp(ncid, varid_a, &
       a_valid_range, a_fill_value)
    INTEGER, INTENT(in) :: ncid, varid_a
    REAL(sp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value

    IF (PRESENT(a_valid_range)) &
         CALL handle_ncdf_err(nf_put_att_real(ncid, varid_a, &
         'valid_range', NF_REAL, 2, a_valid_range(1)), 640)
    IF (PRESENT(a_fill_value)) &
         CALL handle_ncdf_err(nf_put_att_real(ncid, varid_a, &
         '_FillValue', NF_REAL, 1, a_fill_value), 643)
  END SUBROUTINE put_cf_array_attributes_sp

  SUBROUTINE put_cf_array_attributes_dp(ncid, varid_a, &
       a_valid_range, a_fill_value)
    INTEGER, INTENT(in) :: ncid, varid_a
    REAL(dp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value

    IF (PRESENT(a_valid_range)) &
         CALL handle_ncdf_err(nf_put_att_double(ncid, varid_a, &
         'valid_range', nf_double, 2, a_valid_range(1)), 653)
    IF (PRESENT(a_fill_value)) &
         CALL handle_ncdf_err(nf_put_att_double(ncid, varid_a, &
         '_FillValue', nf_double, 1, a_fill_value), 656)
  END SUBROUTINE put_cf_array_attributes_dp

  SUBROUTINE put_cf_array_attributes_i4(ncid, varid_a, &
       a_valid_range, a_fill_value)
    INTEGER, INTENT(in) :: ncid, varid_a
    INTEGER(i4), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value

    IF (PRESENT(a_valid_range)) &
         CALL handle_ncdf_err(nf_put_att_int(ncid, varid_a, &
         'valid_range', NF_INT, 2, a_valid_range(1)), 666)
    IF (PRESENT(a_fill_value)) &
         CALL handle_ncdf_err(nf_put_att_int(ncid, varid_a, &
         '_FillValue', NF_INT, 1, a_fill_value), 669)
  END SUBROUTINE put_cf_array_attributes_i4

  SUBROUTINE cf_latlon_geometry_setup(ncid, varid_a, a_shape, a_dimids, &
       varid_lat, varid_lon)
    INTEGER, INTENT(in) :: ncid, varid_a, a_shape(2)
    INTEGER, INTENT(inout) :: varid_lat, varid_lon
    INTEGER, INTENT(in) :: a_dimids(2)

    CHARACTER(*), PARAMETER :: lat_lon_coord = 'lat lon'
    INTEGER :: ndims

    ndims = SIZE(a_shape)
    CALL handle_ncdf_err(nf_put_att_text(ncid, varid_a, &
         'coordinates', LEN(lat_lon_coord), lat_lon_coord), 683)
    CALL handle_ncdf_err(nf_def_var(ncid, 'lat', NF_DOUBLE, ndims, &
         a_dimids, varid_lat), 685)
    CALL handle_ncdf_err(nf_def_var(ncid, 'lon', NF_DOUBLE, ndims, &
         a_dimids, varid_lon), 687)
  END SUBROUTINE

  SUBROUTINE put_cf_geometry(ncid, varid_lat, varid_lon, &
         alon, alat)
    INTEGER, INTENT(inout) :: ncid, varid_lat, varid_lon
    REAL(dp), INTENT(in) :: alon(:,:), alat(:,:)

    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_lat, alat), 695)
    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_lon, alon), 696)
  END SUBROUTINE put_cf_geometry

  SUBROUTINE handle_ncdf_err(errcode, line)
    INTEGER, INTENT(in) :: errcode, line
    IF (errcode .NE. nf_noerr) THEN
      CALL abort_ppm(nf_strerror(errcode), filename, line)
    END IF
  END SUBROUTINE handle_ncdf_err

END MODULE ppm_ncdf_dump
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
