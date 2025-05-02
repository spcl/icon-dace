dnl acx_libxml2_workaround.m4 --- test whether libxml2 can be
dnl run-time-patched to work around a known problem executing division
dnl by zero during xmlXPathInit
dnl
dnl Copyright  (C)  2023  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf libxml2 xmlXPathInit SIGFPE
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
dnl ACX_LIBXML2_WORKAROUND([TEST-SOURCE-DIR=config/checksrc],
dnl                 [LIST-OF-WORK-AROUND-TESTS],
dnl                 [ACTION-IF-WORKAROUND-SUCCEEDS],
dnl                 [ACTION-IF-WORKAROUND-FAILED=AC_MSG_FAILURE],
dnl                 [WORKAROUND-DIR=src/mpi-workarounds])
dnl
dnl Assumes CC and CFLAGS are setup for MPI.
dnl Then builds files needed for work-around in addition to test programs
dnl and checks if linking with those fixes the issue.
dnl See ACX_MPI_DEFECTS for further description of tests.
AC_DEFUN([ACX_LIBXML2_WORKAROUND],
  [AC_LANG_PUSH([C])
   dnl 1. detect libxml2 version
   dnl currently skipped, only one implementation known to be needed
   dnl so far
   acx_saved_LIBS=$LIBS
   acx_saved_CPPFLAGS=$CPPFLAGS
   acx_is_affected_libxml2=:
   dnl 2. copy workaround source into place for configure logic
   acx_workaround_dir="$srcdir/m4_default([$5],[src/mpi-workarounds])/"
   for acx_workaround in xt_xmlXPathInit xt_xmlInitParser \
     xt_xmlInitParser.c_def
   do
     AS_IF([test "$acx_workaround" = xt_xmlInitParser.c_def],
       [acx_workaround=xt_xmlInitParser
        CPPFLAGS="$CPPFLAGS -DXT_LIBXML_INCLUDE_SUBDIR"])
     cp "$acx_workaround_dir/${acx_workaround}.c" conftest.c
     dnl 3. build fixed xt_xmlXPathInit.c
     AC_COMPILE_IFELSE(,
       [ACX_MV_OBJ([conftest],[xt_workaround])],
       [rm -f core conftest.err conftest.$ac_objext conftest.$ac_ext
        continue])
     dnl 4. verify listed tests now succeed
     LIBS="xt_workaround.$ac_objext -lxml2 $LIBS"
     acx_libxml2_check_dir=`echo "acx_mpi_check_src_" | sed 's@/@<:@^/@:>@*$[]@@'`
     ACX_MPI_DEFECTS([$acx_libxml2_check_dir],, [acx_is_affected_libxml2=false], [$2])
     rm -f xt_workaround.*
   done
   LIBS=$acx_saved_LIBS
   dnl 5. execute action corresponding to outcome of tests
   AS_IF([$acx_is_affected_libxml2],
     [$3
      AS_IF([test "x$CPPFLAGS" != "x$acx_saved_CPPFLAGS"],
        [AC_DEFINE([XT_LIBXML_INCLUDE_SUBDIR],,
          [are libxml header files buried in sub-directory?])])],
     [m4_default([$4],
        [AC_MSG_FAILURE(
           [Cannot apply libxml2 SIGFPE bug work-around.])])])
   CPPFLAGS=$acx_saved_CPPFLAGS
   AC_LANG_POP([C])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
