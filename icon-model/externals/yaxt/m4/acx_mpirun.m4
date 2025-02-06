dnl acx_mpirun.m4 --- check whether launching MPI programs works
dnl
dnl Copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl
dnl ACX_MPIRUN([NUM_TASKS=4],[ACTION-IF-WORKING],
dnl            [ACTION-IF-FAILED=AC_MSG_FAILURE])
dnl
dnl First determines mpi launcher program (and sets MPI_LAUNCH to its path),
dnl then builds and runs simple program with 4 tasks
dnl (or NUM_TASKS if specified).
dnl MPI_LAUNCH is either set to a valid MPI launcher program path (unless
dnl cross-compiling, because it can't be tested then) or true (when it
dnl doesn't work).
dnl
dnl TODO: instead of setting C language, perform test for active AC_LANG
AC_DEFUN([ACX_MPIRUN],
  [AC_PATH_PROGS([MPI_LAUNCH],[mpirun mpiexec],[true])
   AC_ARG_VAR([MPI_LAUNCH],[absolute path to launcher for MPI programs, must be working unless configuring in cross-compilation mode])
   AS_IF([test x"$MPI_LAUNCH" = xtrue],
     [MPI_LAUNCH_failMsg="Failed to find MPI launcher"],
     [AS_IF([test x"$cross_compiling" = xno],
        [AC_MSG_CHECKING([if $MPI_LAUNCH works])
         AS_IF([test -x `echo "$MPI_LAUNCH" | sed -e 's/@<:@ 	@:>@.*//'`],
           [AC_LANG_PUSH([C])
            AC_LINK_IFELSE([AC_LANG_SOURCE([
@%:@include <stdio.h>
@%:@include <stdlib.h>

@%:@include <mpi.h>

@%:@define xmpi(ret)           \\
  do {                      \\
    if (ret != MPI_SUCCESS) \\
      exit(EXIT_FAILURE);   \\
  } while (0)

int
main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  char *numarg = argv@<:@1@:>@;
  int cmdnum = atoi(numarg);
  int procnum = 1;
  xmpi(MPI_Allreduce(MPI_IN_PLACE, &procnum, 1, MPI_INT, MPI_SUM,
                     MPI_COMM_WORLD));
  xmpi(MPI_Finalize());
  return (procnum == cmdnum)?EXIT_SUCCESS:EXIT_FAILURE;
}
])],
              [acx_mpirun_test="$MPI_LAUNCH -n m4_ifval([$1],[$1],[4]) ./conftest$EXEEXT m4_ifval([$1],[$1],[4])"
               AS_IF([expr "$ac_link" : '.*/libtool --mode=link' >/dev/null],
                 [acx_mpirun_test=`echo "$ac_link" | sed -e 's@\(.*/libtool --mode=\)link.*@\1@'`"execute $acx_mpirun_test"])
               _AC_RUN_LOG([LIBC_FATAL_STDERR_=1 $acx_mpirun_test >&2],[echo "running $acx_mpirun_test"])
               AS_IF([test $ac_status -eq 0],,
                 [MPI_LAUNCH=true ; MPI_LAUNCH_failMsg="Failed to run MPI programs"])],
              [MPI_LAUNCH=true ; MPI_LAUNCH_failMsg="Cannot compile or link simple MPI program"])
            AC_LANG_POP([C])],
           [MPI_LAUNCH_failMsg="Cannot execute $MPI_LAUNCH" ; MPI_LAUNCH=true])
         AS_IF([test x"$MPI_LAUNCH" = xtrue],
           [AC_MSG_RESULT([no])],
           [AC_MSG_RESULT([yes])])])])
   AS_IF([test x"$MPI_LAUNCH" = xtrue],
     [m4_ifval([$3],[$3],[AC_MSG_FAILURE([$MPI_LAUNCH_failMsg])])],
     [m4_ifval([$2],[$2],[:])])
])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
