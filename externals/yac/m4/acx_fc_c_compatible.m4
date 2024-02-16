# ACX_FC_C_COMPATIBLE([ACTION-IF-SUCCESS],
#                     [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the Fortran compiler can link objects compiled with the C
# compiler. Tries to compile a simple C code with the C compiler and to link the
$ resulting object into a Fortran program with the Fortran compiler.
#
# The implementation implies that the Fortran compiler supports the BIND(C)
# attribute.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_c_compatible variable.
#
AC_DEFUN([ACX_FC_C_COMPATIBLE],
  [AC_CACHE_CHECK(
     [whether Fortran compiler can link objects compiled with C compiler],
     [acx_cv_fc_c_compatible],
     [_ACX_FC_C_COMPATIBLE
      acx_cv_fc_c_compatible=$acx_fc_c_compatiable])
   AS_VAR_IF([acx_cv_fc_c_compatible], [yes], [$1],
     [m4_default([$2],
        [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
C compiler])])])])

# ACX_FC_C_COMPATIBLE_MPI([MPIRUN = true],
#                         [ACTION-IF-SUCCESS],
#                         [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether Fortran and C MPI libraries are compatible. Tries to compile
# a C function that make calls to MPI functions and to link the resulting
# object into a Fortran program. To make the result of the test more reliable,
# the Fortran program is run using the MPI launch command MPIRUN (defaults to
# "true", which unconditionally exits with a zero exit status) appended with
# two additional arguments: "-n 1" and the name of the linked executable. The
# result is either "yes" or "no".
#
# The implementation implies that the Fortran compiler supports the BIND(C)
# attribute and can link object files compiled with the C compiler.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_c_compatible_mpi variable.
#
AC_DEFUN([ACX_FC_C_COMPATIBLE_MPI],
  [AC_CACHE_CHECK([whether Fortran and C MPI libraries are compatible],
     [acx_cv_fc_c_compatible_mpi],
     [acx_cv_fc_c_compatible_mpi=no
      _ACX_FC_C_COMPATIBLE([[
#include <mpi.h>
void conftest_foo()
{ int world_size, world_rank, name_len;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Init(0, 0);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Get_processor_name(processor_name, &name_len);
  MPI_Finalize(); }]],
        [AC_TRY_COMMAND([m4_default([$1], [true]) -n 1 ./conftest$ac_exeext dnl
>&AS_MESSAGE_LOG_FD])])
      acx_cv_fc_c_compatible_mpi=$acx_fc_c_compatiable])
   AS_VAR_IF([acx_cv_fc_c_compatible_mpi], [yes], [$2],
     [m4_default([$3],
        [AC_MSG_FAILURE([Fortran and C MPI libraries are not compatible])])])])

# ACX_FC_C_COMPATIBLE_OPENMP([ACTION-IF-SUCCESS],
#                            [ACTION-IF-FAILURE = FAILURE])
# -----------------------------------------------------------------------------
# Checks whether the Fortran compiler can link C code that uses OpenMP. Tries
# to compile a C function that make calls to OpenMP functions and to link the
# resulting object into a Fortran program.
#
# The implementation implies that the Fortran compiler supports the BIND(C)
# attribute and can link object files compiled with the C compiler.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# The result is cached in the acx_cv_fc_c_compatible_openmp variable.
#
AC_DEFUN([ACX_FC_C_COMPATIBLE_OPENMP],
  [AC_CACHE_CHECK(
     [whether Fortran compiler can link C code that uses OpenMP],
     [acx_cv_fc_c_compatible_openmp],
     [_ACX_FC_C_COMPATIBLE([[
#include <omp.h>
void conftest_foo () {omp_get_num_threads();}]])
      acx_cv_fc_c_compatible_openmp=$acx_fc_c_compatiable])
   AS_VAR_IF([acx_cv_fc_c_compatible_openmp], [yes], [$1],
     [m4_default([$2],
        [AC_MSG_FAILURE(
           [Fortran compiler cannot link C code that uses OpenMP])])])])

# _ACX_FC_C_COMPATIBLE([FOO-C-CODE = "void conftest_foo(){}"]
#                      [EXTRA-ACTIONS])
# -----------------------------------------------------------------------------
# Checks whether the Fortran compiler can link a program that makes a call to
# a parameterless C function "conftest_foo", which returns void, using the
# ISO_C_BINDING module. First, tries to compile FOO-C-CODE (defaults to a dummy
# function "conftest_foo") with the C compiler and to link the resulting object
# into a Fortran program with the Fortran compiler. The result is either "yes"
# or "no". If you need to run extra commands upon successful linking (e.g. you
# need to run the result of the linking, i.e. "./conftest$ac_exeext"), you can
# put them as EXTRA-ACTIONS argument. In that case the result of the macro will
# be "yes" only if the exit code of the last command listed in EXTRA-ACTIONS is
# zero.
#
# The result is stored in the acx_fc_c_compatiable variable.
#
m4_define([_ACX_FC_C_COMPATIBLE],
  [AC_LANG_ASSERT([Fortran])dnl
   AC_REQUIRE([AC_PROG_CC])dnl
   AC_LANG_PUSH([C])
   acx_fc_c_compatiable=no
   AC_COMPILE_IFELSE([AC_LANG_SOURCE(
     [m4_default([$1], [[void conftest_foo(){}]])])],
     [AC_LANG_POP([C])
      mv ./conftest.$ac_objext ./conftest_c.$ac_objext
      acx_save_LIBS=$LIBS; LIBS="./conftest_c.$ac_objext $LIBS"
      AC_LINK_IFELSE(
        [AC_LANG_SOURCE(
[[      program conftest
      implicit none
      interface
      subroutine conftest_foo() bind(c)
      end subroutine
      end interface
      call conftest_foo()
      end program]])],
        [m4_ifval([$2],
           [$2
            AS_IF([test $? -eq 0], [acx_fc_c_compatiable=yes])],
           [acx_fc_c_compatiable=yes])])
      LIBS=$acx_save_LIBS
      rm -f conftest_c.$ac_objext
      AC_LANG_PUSH([C])])
   AC_LANG_POP([C])])
