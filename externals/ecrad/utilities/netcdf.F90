MODULE netcdf
  IMPLICIT NONE
  ! REF: https://github.com/Unidata/netcdf-fortran/blob/main/fortran/netcdf_constants.F90
  INTEGER, PARAMETER :: nf90_noerr = 0
  INTEGER, PARAMETER :: nf90_nowrite = 0
  INTEGER, PARAMETER :: nf90_clobber = 0
  INTEGER, PARAMETER :: nf90_max_var_dims = 1024
  INTEGER, PARAMETER :: nf90_enotvar = -49
  INTEGER, PARAMETER :: nf90_global = 0
  INTEGER, PARAMETER :: nf90_double = 6
  INTEGER, PARAMETER :: nf90_float = 5
  INTEGER, PARAMETER :: nf90_byte = 1
  INTEGER, PARAMETER :: nf90_int = 4
  INTEGER, PARAMETER :: nf90_short = 3
  INTERFACE nf90_def_var
    MODULE PROCEDURE stub_0
    MODULE PROCEDURE stub_1
    MODULE PROCEDURE stub_2
    MODULE PROCEDURE stub_3
    MODULE PROCEDURE stub_4
    MODULE PROCEDURE stub_5
    MODULE PROCEDURE stub_6
    MODULE PROCEDURE stub_7
    MODULE PROCEDURE stub_8
    MODULE PROCEDURE stub_9
    MODULE PROCEDURE stub_10
    MODULE PROCEDURE stub_11
    MODULE PROCEDURE stub_12
    MODULE PROCEDURE stub_13
    MODULE PROCEDURE stub_14
    MODULE PROCEDURE stub_15
    MODULE PROCEDURE stub_16
    MODULE PROCEDURE stub_17
    MODULE PROCEDURE stub_18
    MODULE PROCEDURE stub_19
    MODULE PROCEDURE stub_20
    MODULE PROCEDURE stub_21
    MODULE PROCEDURE stub_22
    MODULE PROCEDURE stub_23
    MODULE PROCEDURE stub_24
  END INTERFACE nf90_def_var
  INTERFACE nf90_get_att
    MODULE PROCEDURE stub_0
    MODULE PROCEDURE stub_1
    MODULE PROCEDURE stub_2
    MODULE PROCEDURE stub_3
    MODULE PROCEDURE stub_4
    MODULE PROCEDURE stub_5
    MODULE PROCEDURE stub_6
    MODULE PROCEDURE stub_7
    MODULE PROCEDURE stub_8
    MODULE PROCEDURE stub_9
    MODULE PROCEDURE stub_10
    MODULE PROCEDURE stub_11
    MODULE PROCEDURE stub_12
    MODULE PROCEDURE stub_13
    MODULE PROCEDURE stub_14
    MODULE PROCEDURE stub_15
    MODULE PROCEDURE stub_16
    MODULE PROCEDURE stub_17
    MODULE PROCEDURE stub_18
    MODULE PROCEDURE stub_19
    MODULE PROCEDURE stub_20
    MODULE PROCEDURE stub_21
    MODULE PROCEDURE stub_22
    MODULE PROCEDURE stub_23
    MODULE PROCEDURE stub_24
  END INTERFACE nf90_get_att
  INTERFACE nf90_get_var
    MODULE PROCEDURE stub_0
    MODULE PROCEDURE stub_1
    MODULE PROCEDURE stub_2
    MODULE PROCEDURE stub_3
    MODULE PROCEDURE stub_4
    MODULE PROCEDURE stub_5
    MODULE PROCEDURE stub_6
    MODULE PROCEDURE stub_7
    MODULE PROCEDURE stub_8
    MODULE PROCEDURE stub_9
    MODULE PROCEDURE stub_10
    MODULE PROCEDURE stub_11
    MODULE PROCEDURE stub_12
    MODULE PROCEDURE stub_13
    MODULE PROCEDURE stub_14
    MODULE PROCEDURE stub_15
    MODULE PROCEDURE stub_16
    MODULE PROCEDURE stub_17
    MODULE PROCEDURE stub_18
    MODULE PROCEDURE stub_19
    MODULE PROCEDURE stub_20
    MODULE PROCEDURE stub_21
    MODULE PROCEDURE stub_22
    MODULE PROCEDURE stub_23
    MODULE PROCEDURE stub_24
  END INTERFACE nf90_get_var
  INTERFACE nf90_put_var
    MODULE PROCEDURE stub_0
    MODULE PROCEDURE stub_1
    MODULE PROCEDURE stub_2
    MODULE PROCEDURE stub_3
    MODULE PROCEDURE stub_4
    MODULE PROCEDURE stub_5
    MODULE PROCEDURE stub_6
    MODULE PROCEDURE stub_7
    MODULE PROCEDURE stub_8
    MODULE PROCEDURE stub_9
    MODULE PROCEDURE stub_10
    MODULE PROCEDURE stub_11
    MODULE PROCEDURE stub_12
    MODULE PROCEDURE stub_13
    MODULE PROCEDURE stub_14
    MODULE PROCEDURE stub_15
    MODULE PROCEDURE stub_16
    MODULE PROCEDURE stub_17
    MODULE PROCEDURE stub_18
    MODULE PROCEDURE stub_19
    MODULE PROCEDURE stub_20
    MODULE PROCEDURE stub_21
    MODULE PROCEDURE stub_22
    MODULE PROCEDURE stub_23
    MODULE PROCEDURE stub_24
  END INTERFACE nf90_put_var
  INTERFACE nf90_put_att
    MODULE PROCEDURE stub_0
    MODULE PROCEDURE stub_1
    MODULE PROCEDURE stub_2
    MODULE PROCEDURE stub_3
    MODULE PROCEDURE stub_4
    MODULE PROCEDURE stub_5
    MODULE PROCEDURE stub_6
    MODULE PROCEDURE stub_7
    MODULE PROCEDURE stub_8
    MODULE PROCEDURE stub_9
    MODULE PROCEDURE stub_10
    MODULE PROCEDURE stub_11
    MODULE PROCEDURE stub_12
    MODULE PROCEDURE stub_13
    MODULE PROCEDURE stub_14
    MODULE PROCEDURE stub_15
    MODULE PROCEDURE stub_16
    MODULE PROCEDURE stub_17
    MODULE PROCEDURE stub_18
    MODULE PROCEDURE stub_19
    MODULE PROCEDURE stub_20
    MODULE PROCEDURE stub_21
    MODULE PROCEDURE stub_22
    MODULE PROCEDURE stub_23
    MODULE PROCEDURE stub_24
  END INTERFACE nf90_put_att
  CONTAINS
  FUNCTION stub_0(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    REAL, INTENT(IN) :: a_2
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_0
    stub_0 = nf90_noerr
  END FUNCTION stub_0
  FUNCTION stub_1(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2(:)
    CHARACTER(LEN = *), INTENT(IN) :: a_3(:)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_1
    stub_1 = nf90_noerr
  END FUNCTION stub_1
  FUNCTION stub_2(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    INTEGER(KIND = SELECTED_INT_KIND(1)), INTENT(IN) :: a_3
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_2
    stub_2 = nf90_noerr
  END FUNCTION stub_2
  FUNCTION stub_3(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    REAL, INTENT(IN) :: a_3
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_3
    stub_3 = nf90_noerr
  END FUNCTION stub_3
  FUNCTION stub_4(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    DOUBLE PRECISION, INTENT(IN) :: a_2(:, :, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_4
    stub_4 = nf90_noerr
  END FUNCTION stub_4
  FUNCTION stub_5(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    DOUBLE PRECISION, INTENT(IN) :: a_2(:)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_5
    stub_5 = nf90_noerr
  END FUNCTION stub_5
  FUNCTION stub_6(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    INTEGER, INTENT(IN) :: a_2(:, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_6
    stub_6 = nf90_noerr
  END FUNCTION stub_6
  FUNCTION stub_7(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    INTEGER, INTENT(IN) :: a_2(:, :, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_7
    stub_7 = nf90_noerr
  END FUNCTION stub_7
  FUNCTION stub_8(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    REAL, INTENT(IN) :: a_2(:, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_8
    stub_8 = nf90_noerr
  END FUNCTION stub_8
  FUNCTION stub_9(a_0, a_1, a_2, a_3, a_4, start, count)
    INTEGER, INTENT(IN) :: a_0
    CHARACTER(LEN = *), INTENT(IN) :: a_1
    INTEGER, INTENT(IN) :: a_2
    INTEGER, INTENT(IN) :: a_3(:)
    INTEGER, INTENT(IN) :: a_4
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_9
    stub_9 = nf90_noerr
  END FUNCTION stub_9
  FUNCTION stub_10(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    DOUBLE PRECISION, INTENT(IN) :: a_2(:, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_10
    stub_10 = nf90_noerr
  END FUNCTION stub_10
  FUNCTION stub_11(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2(:)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_11
    stub_11 = nf90_noerr
  END FUNCTION stub_11
  FUNCTION stub_12(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_12
    stub_12 = nf90_noerr
  END FUNCTION stub_12
  FUNCTION stub_13(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    INTEGER(KIND = SELECTED_INT_KIND(4)), INTENT(IN) :: a_3
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_13
    stub_13 = nf90_noerr
  END FUNCTION stub_13
  FUNCTION stub_14(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    INTEGER, INTENT(IN) :: a_2(:)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_14
    stub_14 = nf90_noerr
  END FUNCTION stub_14
  FUNCTION stub_15(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    REAL, INTENT(IN) :: a_2(:, :, :, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_15
    stub_15 = nf90_noerr
  END FUNCTION stub_15
  FUNCTION stub_16(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    DOUBLE PRECISION, INTENT(IN) :: a_2
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_16
    stub_16 = nf90_noerr
  END FUNCTION stub_16
  FUNCTION stub_17(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2(:, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_17
    stub_17 = nf90_noerr
  END FUNCTION stub_17
  FUNCTION stub_18(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    CHARACTER(LEN = *), INTENT(IN) :: a_3
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_18
    stub_18 = nf90_noerr
  END FUNCTION stub_18
  FUNCTION stub_19(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    DOUBLE PRECISION, INTENT(IN) :: a_3
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_19
    stub_19 = nf90_noerr
  END FUNCTION stub_19
  FUNCTION stub_20(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    INTEGER, INTENT(IN) :: a_2
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_20
    stub_20 = nf90_noerr
  END FUNCTION stub_20
  FUNCTION stub_21(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    DOUBLE PRECISION, INTENT(IN) :: a_2(:, :, :, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_21
    stub_21 = nf90_noerr
  END FUNCTION stub_21
  FUNCTION stub_22(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    REAL, INTENT(IN) :: a_2(:, :, :)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_22
    stub_22 = nf90_noerr
  END FUNCTION stub_22
  FUNCTION stub_23(a_0, a_1, a_2, a_3, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    CHARACTER(LEN = *), INTENT(IN) :: a_2
    INTEGER, INTENT(IN) :: a_3
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_23
    stub_23 = nf90_noerr
  END FUNCTION stub_23
  FUNCTION stub_24(a_0, a_1, a_2, start, count)
    INTEGER, INTENT(IN) :: a_0
    INTEGER, INTENT(IN) :: a_1
    REAL, INTENT(IN) :: a_2(:)
    INTEGER, OPTIONAL :: start(:), count(:)
    INTEGER :: stub_24
    stub_24 = nf90_noerr
  END FUNCTION stub_24
  FUNCTION nf90_open(path, mode, ncid, chunksize)
    CHARACTER(LEN = *), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: mode
    INTEGER, INTENT(OUT) :: ncid
    INTEGER, OPTIONAL, INTENT(INOUT) :: chunksize
    INTEGER :: nf90_open
    nf90_open = nf90_noerr
  END FUNCTION nf90_open
  FUNCTION nf90_create(path, cmode, ncid, initialsize, chunksize)
    CHARACTER(LEN = *), INTENT(IN) :: path
    INTEGER, INTENT(IN) :: cmode
    INTEGER, INTENT(OUT) :: ncid
    INTEGER, OPTIONAL, INTENT(IN) :: initialsize
    INTEGER, OPTIONAL, INTENT(INOUT) :: chunksize
    INTEGER :: nf90_create
    nf90_create = nf90_noerr
  END FUNCTION nf90_create
  FUNCTION nf90_strerror(ncerr)
    CHARACTER(LEN = 80) :: nf90_strerror
    INTEGER, INTENT(IN) :: ncerr
    nf90_strerror = "!"
  END FUNCTION nf90_strerror
  FUNCTION nf90_close(ncid)
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: nf90_close
    nf90_close = nf90_noerr
  END FUNCTION nf90_close
  FUNCTION nf90_enddef(ncid, h_minfree, v_align, v_minfree, r_align)
    INTEGER, INTENT(IN) :: ncid
    INTEGER, OPTIONAL, INTENT(IN) :: h_minfree, v_align, v_minfree, r_align
    INTEGER :: nf90_enddef
    nf90_enddef = nf90_noerr
  END FUNCTION nf90_enddef
  FUNCTION nf90_inq_varid(ncid, name, varid)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN = *), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: varid
    INTEGER :: nf90_inq_varid
    nf90_inq_varid = nf90_noerr
  END FUNCTION nf90_inq_varid
  FUNCTION nf90_inquire_variable(ncid, varid, name, xtype, ndims, dimids, nAtts)
    INTEGER, INTENT(IN) :: ncid, varid
    CHARACTER(LEN = *), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT) :: xtype, ndims
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: dimids
    INTEGER, OPTIONAL, INTENT(OUT) :: nAtts
    INTEGER :: nf90_inquire_variable
    nf90_inquire_variable = nf90_noerr
  END FUNCTION nf90_inquire_variable
  FUNCTION nf90_inquire_dimension(ncid, dimid, name, len)
    INTEGER, INTENT(IN) :: ncid, dimid
    CHARACTER(LEN = *), OPTIONAL, INTENT(OUT) :: name
    INTEGER, OPTIONAL, INTENT(OUT) :: len
    INTEGER :: nf90_inquire_dimension
    nf90_inquire_dimension = nf90_noerr
  END FUNCTION nf90_inquire_dimension
  FUNCTION nf90_inquire_attribute(ncid, varid, name, xtype, len, attnum)
    INTEGER, INTENT(IN) :: ncid, varid
    CHARACTER(LEN = *), INTENT(IN) :: name
    INTEGER, INTENT(OUT), OPTIONAL :: xtype, len, attnum
    INTEGER :: nf90_inquire_attribute
    nf90_inquire_attribute = nf90_noerr
  END FUNCTION nf90_inquire_attribute
  FUNCTION nf90_inq_dimid(ncid, name, dimid)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN = *), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: dimid
    INTEGER :: nf90_inq_dimid
    nf90_inq_dimid = nf90_noerr
  END FUNCTION nf90_inq_dimid
  FUNCTION nf90_inq_attname(ncid, varid, attnum, name)
    INTEGER, INTENT(IN) :: ncid, varid, attnum
    CHARACTER(LEN = *), INTENT(OUT) :: name
    INTEGER :: nf90_inq_attname
    nf90_inq_attname = nf90_noerr
  END FUNCTION nf90_inq_attname
  FUNCTION nf90_def_dim(ncid, name, len, dimid)
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN = *), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: len
    INTEGER, INTENT(OUT) :: dimid
    INTEGER :: nf90_def_dim
    nf90_def_dim = nf90_noerr
  END FUNCTION nf90_def_dim
  FUNCTION nf90_copy_att(ncid_in, varid_in, name, ncid_out, varid_out)
    INTEGER, INTENT(IN) :: ncid_in, varid_in
    CHARACTER(LEN = *), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: ncid_out, varid_out
    INTEGER :: nf90_copy_att
    nf90_copy_att = nf90_noerr
  END FUNCTION nf90_copy_att
END MODULE netcdf
