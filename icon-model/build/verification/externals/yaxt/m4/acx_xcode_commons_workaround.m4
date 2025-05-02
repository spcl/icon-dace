dnl acx_xcode_commons_workaround.m4 --- test whether -Wl,-commons,use_dylibs is needed
dnl                                     for correct linking
dnl
dnl Copyright  (C)  2024  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf commons xcode
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
dnl ACX_XCODE_COMMONS_WORKAROUND([TEST-SOURCE-DIR=config/checksrc],
dnl                 [LIST-OF-WORK-AROUND-TESTS],
dnl                 [ACTION-IF-WORKAROUND-SUCCEEDS],
dnl                 [ACTION-IF-WORKAROUND-FAILED=AC_MSG_FAILURE])
dnl
dnl Assumes CC and CFLAGS are setup for MPI.
dnl Then flag combinations possibly constituting a work-around
dnl and checks if linking with those fixes the issue.
dnl See ACX_MPI_DEFECTS for further description of tests.
AC_DEFUN([ACX_XCODE_COMMONS_WORKAROUND],
  [AC_LANG_PUSH([C])
   dnl 1. detect ld version
   dnl currently skipped, only one implementation known to be needed
   dnl so far
   acx_saved_LDFLAGS=$LDFLAGS
   acx_is_affected_xcode=false
   dnl 2. test amended set of LDFLAGS
   for acx_workaround in -Wl,-commons,use_dylibs \
     -Wl,-commons,use_dylibs\ -Wl,-ld_classic
   do
     LDFLAGS="${LDFLAGS+$LDFLAGS }${acx_workaround}"
     dnl 3. verify listed tests now succeed
     ACX_MPI_DEFECTS([$1],[acx_is_affected_xcode=:],[:],[$2])
     LDFLAGS=$acx_saved_LDFLAGS
     AS_IF([$acx_is_affected_xcode],
       [break])
   done
   dnl 5. execute action corresponding to outcome of tests
   AS_IF([$acx_is_affected_xcode],
     [$3],
     [m4_default([$4],
        [AC_MSG_FAILURE(
           [Cannot apply Xcode common linking workaround.])])])
   AC_LANG_POP([C])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
