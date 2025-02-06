dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords: configure configure.ac autotools
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
dnl Check for doxygen properties
dnl
dnl ACX_CHECK_DOXYGEN
dnl
dnl check if doxygen executable is available
dnl
AC_DEFUN([ACX_CHECK_DOXYGEN],
  [AC_CHECK_PROGS([DOXYGEN], [doxygen], [:])
   AS_IF([test x"$DOXYGEN" = x:],
     [acx_doxygen_version=0.0.0],
     [acx_doxygen_version=`"$DOXYGEN" --version`])])dnl
dnl
dnl ACX_CHECK_DOXYGEN_LAYOUT
dnl
AC_DEFUN([ACX_CHECK_DOXYGEN_LAYOUT],
  [AC_REQUIRE([ACX_CHECK_DOXYGEN])
   AS_IF([test x"$DOXYGEN" = x:],[DOXYFILE_HEADER=Doxyfile-1.5.6],
     [AC_MSG_CHECKING([if doxygen supports layouting])
      AS_IF(["$DOXYGEN" -l /dev/null > /dev/null 2>&1],
        [AC_MSG_RESULT([yes]); DOXYFILE_HEADER=Doxyfile-1.5.8],
        [AC_MSG_RESULT([no]) ; DOXYFILE_HEADER=Doxyfile-1.5.6])])])dnl
dnl
dnl ACX_CHECK_DOXYGEN_HTML_EXTRA_STYLESHEET
dnl
AC_DEFUN([ACX_CHECK_DOXYGEN_HTML_EXTRA_STYLESHEET],
  [AC_REQUIRE([ACX_CHECK_DOXYGEN])
   AC_MSG_CHECKING([if doxygen supports an extra HTML CSS stylesheet])
   AS_VERSION_COMPARE(["$acx_doxygen_version"],[1.8.2],
     [ACX_DOXYFILE_HTML_EXTRA_STYLESHEET=no],
     [ACX_DOXYFILE_HTML_EXTRA_STYLESHEET=yes],
     [ACX_DOXYFILE_HTML_EXTRA_STYLESHEET=yes])
   AC_MSG_RESULT([$ACX_DOXYFILE_HTML_EXTRA_STYLESHEET])])dnl
dnl Local Variables:
dnl mode: autoconf
dnl End:
