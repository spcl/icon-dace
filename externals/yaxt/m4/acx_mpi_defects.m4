dnl acx_mpi_defects.m4 --- check whether MPI has one or more of
dnl                        several known defects
dnl
dnl Copyright  (C)  2014  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
dnl ACX_MPI_DEFECTS([TEST-SOURCE-DIR=$srcdir/config/checksrc],
dnl                 [ACTION-IF-CHECK-SUCCEEDS],
dnl                 [ACTION-IF-CHECK-FAILED=AC_MSG_FAILURE],
dnl                 [LIST-OF-TEST-GLOBS=[*]])
dnl
dnl Requires MPI_LAUNCH program. Also CC/FC must be setup to
dnl build MPI programs.
dnl Builds and runs simple programs from TEST-SOURCE-DIR, each of which
dnl should represent a test for a known defect that affects the library
dnl code.
dnl Each test is built according to it's file suffix as either Fortran
dnl or C MPI program.
dnl Within ACTION-IF-CHECK-SUCCEEDS and ACTION-IF-CHECK-FAILED,
dnl the following variables are set to test-specific values:
dnl acx_subtestname = base file name of the test
dnl acx_mpi_check_src = path to check source file
dnl acx_suffix = file suffix of source file
dnl
dnl Each test source may contain zero or more of the following stanzas
dnl acx_mpirun_num_tasks = N
dnl   specify number of tasks N (positive integer) to run this test with
dnl TODO: extend for F77 and C++
AC_DEFUN([ACX_MPI_DEFECTS],
  [AC_ARG_ENABLE([cross-mpi-defect-checks],
     AS_HELP_STRING([--enable-cross-mpi-defect-checks],
       [Run MPI defect tests even in cross-compilation mode. This is]dnl
[ typically useful in a setup where the mpi-starter (see MPI_LAUNCH)]dnl
[ executes programs on another system than the one running configure]dnl
[ which matches the host specification.]dnl
[ @<:@default=no@:>@]),
     [AS_IF([test x"$enable_cross_mpi_defect_checks" != xno],
        [enable_cross_mpi_defect_checks=yes])],
     [enable_cross_mpi_defect_checks=no])
   AS_IF([test x"$MPI_LAUNCH" = xtrue],
     [AC_MSG_NOTICE([Skipping tests for known MPI defects: MPI launcher unavailable])],
     [test $cross_compiling = yes -a $enable_cross_mpi_defect_checks = no],
     [AC_MSG_NOTICE([Skipping tests for known MPI defects in cross-]dnl
[compilation mode, consider enabling them via ]dnl
[--enable-cross-mpi-defect-checks, if possible])],
     [AC_LANG_PUSH([C])
      AC_MSG_CHECKING([MPI for known defects])
      AC_MSG_RESULT([])
m4_ifval([$4],
  [m4_pushdef([acx_test_pat],
     m4_dquote(m4_foreach([acx_temp],[$4],
       ["]m4_dquote(m4_default([$1],[$srcdir/config/checksrc]))[/"]acx_temp)))],
  [m4_pushdef([acx_test_pat],
     ["]m4_dquote(m4_default([$1],[$srcdir/config/checksrc]))[/"*])])dnl
      AS_FOR([acx_mpi_check_src_],[acx_mpi_check_src],acx_test_pat,[
	AS_IF([test ! -r "]acx_mpi_check_src_[" \
            && test x"]acx_mpi_check_src_[" = x]m4_bpatsubst(m4_dquote(m4_defn([acx_test_pat])),[[*?]],[\\\&]),
	  [break])
	acx_suffix=`echo "acx_mpi_check_src_" | sed 's/^.*\.\(@<:@^.@:>@*\)$/\1/'`
	acx_subtestname=`echo "acx_mpi_check_src_" | sed 's/^.*\/\(@<:@^\/@:>@*\)\.@<:@^.@:>@*/\1/'`
	AC_MSG_CHECKING([$acx_subtestname])
	acx_mpirun_num_tasks=`sed -n '/acx_mpirun_num_tasks *= *\(@<:@0-9@:>@*\)/{
s/.*acx_mpirun_num_tasks *= *\(@<:@0-9@:>@*\).*/\1/
p
q
}
' "acx_mpi_check_src_"`
	AS_IF([test `expr "$acx_mpirun_num_tasks" : "@<:@0-9@:>@@<:@0-9@:>@*$"` -gt 0 \
	       && test "$acx_mpirun_num_tasks" -gt 0],,
	  [acx_mpirun_num_tasks=1])
        acx_mpirun_expected_exitcode=`sed -n '/acx_mpirun_expected_exitcode *= *\(@<:@0-9@:>@*\)/{
s/.*acx_mpirun_expected_exitcode *= *\(@<:@0-9@:>@*\).*/\1/
p
q
}
' "acx_mpi_check_src_"`
	AS_IF([test `expr "$acx_mpirun_expected_exitcode" : "@<:@0-9@:>@@<:@0-9@:>@*$"` -gt 0 \
	       && test "$acx_mpirun_expected_exitcode" -gt 0],,
	  [acx_mpirun_expected_exitcode=0])
dnl check if test should only run if MPI version matches
        acx_mpi_expected_version=`sed -n '/acx_mpi_expected_version *\(==\|>=\|<=\|>\|<\|!=\) *\(@<:@0-9@:>@@<:@0-9.@:>@*\)/{
s/.*acx_mpi_expected_version *\(\(==\|>=\|<=\|>\|<\|!=\) @<:@0-9@:>@@<:@0-9.@:>@*\).*/\1/
p
q
}
' "acx_mpi_check_src_"`
	AS_IF([expr "$acx_mpi_expected_version" : "\(==\|>=\|<=\|>\|<\|!=\) @<:@0-9@:>@@<:@0-9.@:>@*$" >/dev/null],
          [acx_mpi_expected_version=`echo "$acx_mpi_expected_version" \
            | sed -e 's/^==/eq/;s/^>=/ge/;s/^<=/le/;s/^>/gt/;s/^</lt/;]dnl
[s/!=/ne/;s/  */ v/'`],
	  [ASX_VAR_UNSET([acx_mpi_expected_version])])
        AS_IF([test x${acx_mpi_expected_version+set} = xset],
          [AC_COMPUTE_INT([MPI_VERSION], [MPI_VERSION],
             [@%:@include <mpi.h>],
             [AC_MSG_ERROR([cannot determine MPI version!])])
           AC_COMPUTE_INT([MPI_SUBVERSION], [MPI_SUBVERSION],
             [@%:@include <mpi.h>],
             [AC_MSG_ERROR([cannot determine MPI sub-version!])])])
        AS_IF([test x${acx_mpi_expected_version+set} = x \
          || $PERL -e "use strict; use warnings; exit(!(v${MPI_VERSION}.${MPI_SUBVERSION} ${acx_mpi_expected_version}))"],
	  [AS_CASE([$acx_suffix],
	     [c],
	     [cat confdefs.h "acx_mpi_check_src_" >conftest."$acx_suffix"],
	     [f90|F90],[cat "acx_mpi_check_src_" >conftest."$acx_suffix"
	      _AC_LANG_SET(,[Fortran])],
	     [AC_MSG_FAILURE([Unexpected language in MPI check: ${acx_subtestname}.${acx_suffix}])])
           AC_LINK_IFELSE(,
	     [acx_mpirun_num_tasks="$MPI_LAUNCH -n $acx_mpirun_num_tasks ./conftest$EXEEXT"
             AS_IF([expr "$ac_link" : '.*/libtool --mode=link' >/dev/null],
               [acx_mpirun_test=`echo "$ac_link" | sed -e 's@\(.*/libtool --mode=\)link.*@\1@'`"execute $acx_mpirun_test"])
	     _AC_RUN_LOG_LIMIT([LIBC_FATAL_STDERR_=1 $acx_mpirun_num_tasks >&2],[echo "running $acx_mpirun_num_tasks"])
	     AS_IF([test $ac_status -eq $acx_mpirun_expected_exitcode],
	       [acx_mpi_defects_result=okay; acx_mpi_defects_fail=no]m4_ifval([$2],[
	        $2]),
	     [acx_mpi_defects_result=error ; acx_mpi_defects_fail=yes])],
	     [acx_mpi_defects_result=error ; acx_mpi_defects_fail=yes])
           AS_CASE([$acx_suffix],
	     [f90|F90],[_AC_LANG_SET(,[C])])],
        [acx_mpi_defects_result=skipped ; acx_mpi_defects_fail=no])
        AC_MSG_RESULT([$acx_mpi_defects_result])
        AS_IF([test x"$acx_mpi_defects_fail" = xyes],
          [m4_default([$3],[AC_MSG_FAILURE([chosen MPI has known error $acx_subtestname])])])])
m4_popdef([acx_test_pat])dnl
      ASX_VAR_UNSET([acx_mpirun_num_tasks])
      ASX_VAR_UNSET([acx_mpi_check_src])
      ASX_VAR_UNSET([acx_suffix])
      ASX_VAR_UNSET([acx_subtestname])
      AC_LANG_POP([C])
     ])])
dnl dump text documentation of defect test to stderr
dnl ACX_MPI_DEFECTS_DOCUMENT([TEST-DOC-DIR=config/checkdoc])
AC_DEFUN([ACX_MPI_DEFECTS_DOCUMENT],
  [AS_IF([test -r "$srcdir/m4_default([$1],[config/checkdoc])/${acx_subtestname}.txt"],
             [cat "$srcdir/m4_default([$1],[config/checkdoc])/${acx_subtestname}.txt" >&2])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
