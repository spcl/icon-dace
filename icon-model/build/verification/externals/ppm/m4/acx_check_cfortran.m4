dnl acx_check_cfortran.m4 --- test if a program compiled from
dnl                           a main Fortran program and
dnl                           C functions gives expected results
dnl
dnl Copyright  (C)  2013  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
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
dnl helper function to emit Fortran test code
AC_DEFUN([_ACX_CHECK_CFORTRAN_FC],
  [AC_LANG_PUSH([Fortran])
         AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      MODULE conftest_data
      IMPLICIT NONE
      PRIVATE
      REAL :: ri
      PUBLIC :: ri
      END MODULE conftest_data

      FUNCTION cftstf(v) RESULT(r)
      USE conftest_data, ONLY: ri
      REAL, INTENT(in) :: v
      REAL :: r
      r = v * 100.0
      ri = 1.0 / v
      END FUNCTION cftstf])],
           [ACX_MV_OBJ([conftest],[conftest_f])
            save_LIBS=$LIBS
            LIBS="conftest_c.$ac_objext conftest_f.$ac_objext $LIBS"
            AS_CASE([$FC_FPP_FLAG],[-x\ f??-cpp-input],
              [LIBS="-x none $LIBS"])
            AC_RUN_IFELSE([AC_LANG_PROGRAM(,
[      USE conftest_data, ONLY: ri
      IMPLICIT NONE
      INTERFACE
       FUNCTION cftstc(i, v, p, q) RESULT(f)
         INTEGER, INTENT(in) :: i
         REAL, INTENT(in) :: v
         INTEGER, INTENT(out) :: p
         REAL, INTENT(out) :: q
         REAL :: f
       END FUNCTION cftstc
       FUNCTION chtst() result(s)
         CHARACTER(99) :: s
       END FUNCTION chtst
      END INTERFACE
      REAL, PARAMETER :: eps = 10e-6
      REAL :: foo, boo, too
      INTEGER :: bar, baz, i
      CHARACTER(99) :: aaaaaa
      bar = 5
      foo = 0.3
      too = cftstc(bar, foo, baz, boo)
      IF (ABS(baz - NINT(bar * foo)) /= 0) THEN
        WRITE (0, '(2(a,i0))') "error checking, when baz, baz=", baz, &
             ", NINT(bar * foo) =", NINT(bar * foo)
        FLUSH(0)
        CALL err_exit
      END IF
      IF (ABS((ri - 1.0 / (bar * foo)) / ABS(ri)) > eps)  THEN
        WRITE (0, '(2(a,g24.15))') "error checking ri, ri=", ri, ", 1.0 / &
             &(bar * foo) = ", 1.0 / (bar * foo)
        FLUSH(0)
        CALL err_exit
      END IF
      IF (ABS((boo - (bar * foo * 100.0))/ABS(boo)) > eps)  THEN
        WRITE (0, '(2(a,g24.15))') "error checking boo, boo=", boo, &
             ", bar * foo * 100.0 = ", bar * foo * 100.0
        FLUSH(0)
        CALL err_exit
      END IF
      IF (too /= boo) THEN
        WRITE (0, '(2(a,g24.15))') "error checking too vs. boo, too=", too, &
             ", boo = ", boo
        FLUSH(0)
        CALL err_exit
      END IF
      aaaaaa = chtst()
      DO i = 1, 99
        IF (aaaaaa(i:i) /= 'A') THEN
          WRITE (0, '(a,i0,a)') "error checking aaaaaa(", i, ")=", &
              aaaaaa(i:i)
          FLUSH(0)
          CALL err_exit
        END IF
      END DO])],
              [acx_cv_cfortran_works=yes],
              [acx_cv_cfortran_works="error"],
              [AC_MSG_NOTICE([Skipping run test for cfortran.h in cross-compilation mode,])
	       AC_MSG_NOTICE([cfortran.h link test succeeded for $FC.])
               acx_cv_cfortran_works=yes])
	       rm -f "conftest_f.$ac_objext" "conftest_f.$OBJEXT"
dnl Some Fortran compilers create module files not in the current
dnl working directory but in the directory with the object file,
dnl therefore we try to delete everything:
	       rm -f conftest_data* CONFTEST_DATA* 2>/dev/null
	    LIBS=$save_LIBS],
           [acx_cv_cfortran_works="error compiling Fortran subroutine"])
         AC_LANG_POP([Fortran])])

dnl helper function to emit Fortran 77 test code
AC_DEFUN([_ACX_CHECK_CFORTRAN_F77],
  [AC_LANG_PUSH([Fortran 77])
         AC_COMPILE_IFELSE([AC_LANG_SOURCE(
[      REAL FUNCTION CFTSTF(v)
      REAL RI
      COMMON /CFTSTD/ RI
      REAL V
      REAL R
      CFTSTF = V * 100.0
      RI = 1.0 / V
      END FUNCTION CFTSTF])],
           [ACX_MV_OBJ([conftest],[conftest_f])
            save_LIBS=$LIBS
            LIBS="conftest_c.$ac_objext conftest_f.$ac_objext $LIBS"
            AS_CASE([$FC_FPP_FLAG],[-x\ f??-cpp-input],
              [LIBS="-x none $LIBS"])
            AC_RUN_IFELSE([AC_LANG_PROGRAM(,
[      REAL RI
      COMMON /CFTSTD/ RI
      REAL EPS
      PARAMETER(EPS=10E-6)
      REAL FOO, BOO, TOO
      INTEGER BAR, BAZ, I
      CHARACTER(99) AAAAAA
      EXTERNAL CFTSTC, CFTSTF, CHTST, ERR_EXIT
      REAL CFTSTC, CFTSTF
      CHARACTER(99) CHTST
      BAR = 5
      FOO = 0.3
      TOO = CFTSTC(BAR, FOO, BAZ, BOO)
      IF (ABS(BAZ - NINT(BAR * FOO)) /= 0) THEN
        WRITE (0, '(2(A,I0))') "ERROR CHECKING, WHEN BAZ, BAZ=", BAZ,
     &       ", NINT(BAR * FOO) =", NINT(BAR * FOO)
        CALL ERR_EXIT
      END IF
      IF (ABS((RI - 1.0 / (BAR * FOO)) / ABS(RI)) > EPS)  THEN
        WRITE (0, '(2(A,F24.15))') "ERROR CHECKING RI, RI=", RI, ",
     &       1.0 / (BAR * FOO) = ", 1.0 / (BAR * FOO)
        CALL err_exit
      END IF
      IF (ABS((BOO - (BAR * FOO * 100.0))/ABS(BOO)) > EPS)  THEN
        WRITE (0, '(2(A,F24.15))') "ERROR CHECKING BOO, BOO=", BOO,
     &       ", BAR * FOO * 100.0 = ", BAR * FOO * 100.0
        CALL ERR_EXIT
      END IF
      IF (TOO /= BOO) THEN
        WRITE (0, '(2(A,F24.15))') "ERROR CHECKING TOO VS. BOO, TOO=",
     &       TOO, ", BOO = ", BOO
        CALL ERR_EXIT
      END IF
      AAAAAA = CHTST()
      DO i = 1, 99
        IF (AAAAAA(I:I) /= 'A') THEN
          WRITE (0, '(A,I0,2A)') "ERROR CHECKING AAAAAA(", I, ")=",
     &        AAAAAA(I:I)
          CALL ERR_EXIT
        END IF
      END DO])],
              [acx_cv_cfortran_works=yes],
              [acx_cv_cfortran_works="error"],
              [AC_MSG_NOTICE([Skipping run test for cfortran.h in cross-compilation mode,])
	       AC_MSG_NOTICE([cfortran.h link test succeeded for $F77.])
               acx_cv_cfortran_works=yes])
	       rm -f "conftest_f.$ac_objext" "conftest_f.$OBJEXT"
	    LIBS=$save_LIBS],
           [acx_cv_cfortran_works="error compiling Fortran subroutine"])
         AC_LANG_POP([Fortran 77])])

dnl ACX_CHECK_CFORTRAN([OPTIONAL-CFORTRAN-INC-DIR],[ACTION-IF-SUCCESS],
dnl                    [ACTION-IF-FAILED])
dnl Test if C compiler can produce objects that link fine to Fortran programs
dnl when using cfortran.h.
dnl
AC_DEFUN([ACX_CHECK_CFORTRAN],
  [AC_CACHE_CHECK([if C externals constructed with cfortran.h work],
     [acx_cv_cfortran_works],
     [acx_cv_cfortran_works=no
      save_CPPFLAGS=$CPPFLAGS
      CPPFLAGS="-I]m4_ifval([$1],[$1],[$srcdir/include])[ $CPPFLAGS"
dnl build C function
      AC_LANG_PUSH([C])
      AC_COMPILE_IFELSE([AC_LANG_SOURCE([@%:@include "cfortran.h"
@%:@include <math.h>
@%:@include <stdio.h>
@%:@include <stdlib.h>

PROTOCCALLSFFUN1(FLOAT,CFTSTF,cftstf,FLOAT)
#define conftest_F(v) \
  CCALLSFFUN1(CFTSTF,cftstf,FLOAT,(v));

static float
cftstC(int i, float v, int *p, float *q)
{
  float f;
  *p = (int)roundf(v * i);
  *q = f = conftest_F(v * i);
  return f;
}

FCALLSCFUN4(FLOAT,cftstC,CFTSTC,cftstc,INT,FLOAT,PINT,PFLOAT)

/* test string returns */
static const char *
conftest_str_C(void)
{
  static const char msg@<:@100@:>@ = "AAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
  return msg;
}

FCALLSCFUN0(STRING,conftest_str_C,CHTST,chtst)

/* This function is required simply because some Fortran compilers
 * won't stop with exit code n when encountering STOP n */
static void
errExit(void)
{
  exit(1);
}

FCALLSCSUB0(errExit,ERR_EXIT,err_exit)

])],
        [ACX_MV_OBJ([conftest],[conftest_c])
         AC_PROVIDE_IFELSE([AC_PROG_FC],
           [AS_IF([test -n "$FC" -a X"$FC" != Xno],
              [_ACX_CHECK_CFORTRAN_FC],
              [acx_cv_cfortran_works=${acx_cv_cfortran_works-yes}])])
         AC_PROVIDE_IFELSE([AC_PROG_F77],
           [AS_IF([test -n "$F77" -a X"$F77" != Xno AC_PROVIDE_IFELSE(
              [AC_PROG_FC],[-a \( x"$acx_cv_cfortran_works" = xyes dnl
-o -z "$FC" -o X"$FC" = Xno \)])],[_ACX_CHECK_CFORTRAN_F77])])
         rm -f "conftest_c.$ac_objext" "conftest_c.$OBJEXT"
        ],
        [acx_cv_cfortran_works="compiling with cfortran.h failed"])
      AC_LANG_POP([C])
      CPPFLAGS=$save_CPPFLAGS
     ])
   m4_ifval([$3],
     [AS_IF([test x"$acx_cv_cfortran_works" = xyes],[$2],[$3])],
     [AS_CASE([x"$acx_cv_cfortran_works"],
       [x"error"],
       [AC_MSG_FAILURE([Linking/Running with C EXTERNAL built with cfortran.h does not work!])],
       [x"compiling with cfortran.h failed"],
       [AC_MSG_FAILURE([Compilation with cfortran.h is not working!])],
       [x"error compiling Fortran subroutine"],
       [AC_MSG_FAILURE([compilation of simple Fortran source failed!])],
       [xyes],[$2],
       [AC_MSG_FAILURE([Unexpected error when linking C and Fortran via cfortran.h!])])])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
