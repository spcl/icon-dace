# ACX_PROG_MPIRUN([VARIABLE],
#                 [ACTION-IF-SUCCESS],
#                 [ACTION-IF-FAILURE = FAILURE],
#                 [CANDIDATES = mpirun mpiexec srun],
#                 [MPI-JOB-COUNT = 2])
# -----------------------------------------------------------------------------
# Searches for the MPI launcher command. If the argument VARIABLE is given and
# the shell variable VARIABLE is set, checks whether the value it stores is a
# valid MPI launcher command. If VARIABLE is not set, iterates over the values
# of the blank-separated list CANDIDATES and stops when the first valid command
# is found. The value of VARIABLE is never set or changed. The number of MPI
# jobs in each test is equal to MPI-JOB-COUNT (defaults to 2). See
# _ACX_PROG_MPIRUN for more details.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# A positive result of this test is cached in the acx_cv_prog_mpirun variable.
#
AC_DEFUN([ACX_PROG_MPIRUN],
  [_ACX_PROG_MPIRUN([_ACX_PROG_MPIRUN_CHECK_PROGRAM], $@)])

# ACX_PROG_MPIRUN_FC_HEADER([VARIABLE],
#                           [ACTION-IF-SUCCESS],
#                           [ACTION-IF-FAILURE = FAILURE],
#                           [CANDIDATES = mpirun mpiexec srun],
#                           [MPI-JOB-COUNT = 2])
# -----------------------------------------------------------------------------
# The same as ACX_PROG_MPIRUN but the tests (see _ACX_PROG_MPIRUN) are based on
# a Fortran program that includes the header file "mpif.h" to interface the MPI
# functionality, whereas the default Fortran program in ACX_PROG_MPIRUN uses
# the Fortran module "MPI".
#
AC_DEFUN([ACX_PROG_MPIRUN_FC_HEADER],
  [AC_LANG_ASSERT([Fortran])dnl
   _ACX_PROG_MPIRUN([AC_LANG_PROGRAM([], [[      implicit none
      include "mpif.h"
      integer :: s, r, e1, e2
      s = -1
      r = -1
      call MPI_Init(e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif
      call MPI_Comm_size(MPI_COMM_WORLD, s, e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif
      call MPI_Comm_rank(MPI_COMM_WORLD, r, e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif
      if (r == 0) write(*, "(a,i0)") 'conftest: ', s
      call MPI_Finalize(e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif]])], $@)])

# _ACX_PROG_MPIRUN(CHECK-PROGRAM,
#                  [VARIABLE],
#                  [ACTION-IF-SUCCESS],
#                  [ACTION-IF-FAILURE = FAILURE],
#                  [CANDIDATES = mpirun mpiexec srun],
#                  [MPI-JOB-COUNT = 2])
# -----------------------------------------------------------------------------
# Searches for the MPI launcher command as described in ACX_PROG_MPIRUN. Each
# tested command is run appended with two additional arguments:
# "-n MPI-JOB-COUNT" (MPI-JOB-COUNT defaults to 2) and the name of the
# executable of the program CHECK-PROGRAM. The program is expected to print a
# line 'conftest: N', where N is the total number of tasks in the
# MPI_COMM_WORLD communicator. The test is successful if it exits with a zero
# exit status and prints a line 'conftest: MPI-JOB-COUNT'.
#
# If successful, runs ACTION-IF-SUCCESS, otherwise runs ACTION-IF-FAILURE
# (defaults to failing with an error message).
#
# A positive result of this test is cached in the acx_cv_prog_mpirun variable.
#
m4_define([_ACX_PROG_MPIRUN],
  [AC_MSG_CHECKING([for MPI launch program])
   acx_tmp=
   AC_CACHE_VAL([acx_cv_prog_mpirun],
     [AC_LINK_IFELSE([$1],
        [m4_ifval([$2],
           [AS_VAR_SET_IF([$2],
             [set dummy "$$2"],
             [set dummy m4_default([$5], [mpirun mpiexec srun])])],
           [set dummy m4_default([$5], [mpirun mpiexec srun])])
         shift
         for acx_candidate in "$[@]"; do
            _AS_ECHO_LOG(
              [acx_exec_result=`$acx_candidate -n m4_default([$6], [2]) dnl
./conftest$ac_exeext 2>/dev/null`])
            acx_exec_result=dnl
`$acx_candidate -n m4_default([$6], [2]) ./conftest$ac_exeext 2>/dev/null`
            acx_status=$?
            _AS_ECHO_LOG([\$? = $acx_status])
            AS_IF([test $acx_status -eq 0],
              [_AS_ECHO_LOG([\$acx_exec_result = $acx_exec_result])
               _AS_ECHO_LOG(
                 [acx_exec_result=`AS_ECHO(["\$acx_exec_result"]) | dnl
sed -n '/^conftest: m4_default([$6], [2])$/p'`])
               acx_exec_result=`AS_ECHO(["$acx_exec_result"]) | dnl
sed -n '/^conftest: m4_default([$6], [2])$/p'`
               _AS_ECHO_LOG([\$acx_exec_result = $acx_exec_result])
               AS_VAR_IF([acx_exec_result],
                 ['conftest: m4_default([$6], [2])'],
                 [acx_cv_prog_mpirun=$acx_candidate
                  break])])
         done],
        [acx_tmp='failed to link MPI test program'])])
   AS_VAR_SET_IF([acx_cv_prog_mpirun],
     [AC_MSG_RESULT([$acx_cv_prog_mpirun])
      $3],
     [AC_MSG_RESULT([unknown])
      AS_IF([test -n "$acx_tmp"], [AC_MSG_WARN([$acx_tmp])])
      m4_default([$4],
        [AC_MSG_FAILURE([unable to find a valid MPI launch program])])])])

# _ACX_PROG_MPIRUN_CHECK_PROGRAM()
# -----------------------------------------------------------------------------
# Expands into a program in the current language that prints a line
# 'conftest: N', where N is the size of the group associated with the MPI
# communicator MPI_COMM_WORLD. By default, expands to m4_fatal with the message
# saying that _AC_LANG is not supported.
#
AC_DEFUN([_ACX_PROG_MPIRUN_CHECK_PROGRAM],
  [m4_ifdef([$0(]_AC_LANG[)],
     [m4_indir([$0(]_AC_LANG[)], $@)],
     [m4_fatal([the MPI check program is not defined for ]dnl
_AC_LANG[ language])])])

# _ACX_PROG_MPIRUN_CHECK_PROGRAM(C)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_PROG_MPIRUN_CHECK_PROGRAM for C language.
#
m4_define([_ACX_PROG_MPIRUN_CHECK_PROGRAM(C)],
  [AC_LANG_SOURCE([[#include<stdio.h>
#include<mpi.h>
int main(int argc, char **argv) {
  int s = -1, r = -1, e;
  if ((e = MPI_Init(&argc, &argv)) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, e);
  if ((e = MPI_Comm_size(MPI_COMM_WORLD, &s)) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, e);
  if ((e = MPI_Comm_rank(MPI_COMM_WORLD, &r)) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, e);
  if (r == 0)
    printf("conftest: %i\n", s);
  if ((e = MPI_Finalize()) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, e);
  return e; }]])])

# _ACX_PROG_MPIRUN_CHECK_PROGRAM(Fortran)()
# -----------------------------------------------------------------------------
# Implementation of _ACX_PROG_MPIRUN_CHECK_PROGRAM for Fortran language.
#
m4_define([_ACX_PROG_MPIRUN_CHECK_PROGRAM(Fortran)],
  [AC_LANG_PROGRAM([], [[      use mpi
      implicit none
      integer :: s, r, e1, e2
      s = -1
      r = -1
      call MPI_Init(e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif
      call MPI_Comm_size(MPI_COMM_WORLD, s, e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif
      call MPI_Comm_rank(MPI_COMM_WORLD, r, e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif
      if (r == 0) write(*, "(a,i0)") 'conftest: ', s
      call MPI_Finalize(e1)
      if (e1 /= MPI_SUCCESS) then
        call MPI_Abort(MPI_COMM_WORLD, e1, e2)
      endif]])])
