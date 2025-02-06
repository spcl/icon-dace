dnl configure.ac --- configure template for ScalES-PPM
dnl
dnl Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
dnl Process this file with autoconf to produce a configure script.
dnl
dnl we need 2.64 because the LANG macros use autom4te properties only
dnl available since sometime in 2.63b
AC_PREREQ([2.64])
AC_INIT([ppm],[1.0.8.1],[jahns@dkrz.de])
AC_CONFIG_SRCDIR([src/core/core.c])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([1.11 foreign subdir-objects])
AM_MAINTAINER_MODE
dnl
dnl Checks for package options.
dnl ######################################################################


dnl ######################################################################

AC_CONFIG_HEADER([include/config.h])

AC_ARG_ENABLE([mh-file],
  [AS_HELP_STRING([--enable-mh-file=<name>],
   [select specific build and host configuration @<:@default: off@:>@])],
   ,[enable_mh_file=off])
AC_ARG_WITH([setup],
  [AS_HELP_STRING([--with-setup=<name>],
     [use user specific setup (mh_setup)])],,)
AS_IF([test x"$enable_mh_file" = xauto],
  [host_frag=""
   AS_CASE([$host],
     [alpha-*-*], [host_frag=config/mh-alpha],
     [5000-fujitsu-*], [host_frag=config/mh-fujitsu-5000],
     [*-fujitsu-*], [host_frag=config/mh-fujitsu],
     [sr8000-hitachi-*], [host_frag=config/mh-hitachi-sr8000],
     [*90-cray-*], [host_frag=config/mh-cray-c90],
     [sv1-cray-*], [host_frag=config/mh-cray-c90],
     [nv1-cray-*], [host_frag=config/mh-cray-x1],
     [t3e-cray-*], [host_frag=config/mh-cray-t3e],
     [xd1-cray-*], [host_frag=config/mh-cray-xd1],
     [xt3-cray-*], [host_frag=config/mh-cray-xt3],
     [mips*-sgi-irix6.*], [host_frag=config/mh-irix6],
     [*-*-solaris2*], [host_frag=config/mh-solaris],
     [sx*-*-*], [host_frag=config/mh-nec],
     [es*-*-*], [host_frag=config/mh-es],
     [*-ibm-aix*], [host_frag=config/mh-aix],
     [ia64-*-linux-*], [host_frag=config/mh-linux-ia64 ],
     [x86_64-*-linux-*], [host_frag=config/mh-linux-x64],
     [i*86-*-linux-*], [host_frag=config/mh-linux],
     [*-apple-darwin*], [host_frag=config/mh-darwin],
     [*-*-hpux*], [host_frag=config/mh-hpux],
     [AC_MSG_WARN([This configuration is not supported.
Please create a valid config/mh-* file.])])],
  [host_frag="config/$enable_mh_file"])
dnl
AS_IF([test x"$enable_mh_file" != xno -a x"$enable_mh_file" != xoff],
  [AS_IF([test x"$host_frag" != x -a ! -f $srcdir/$host_frag],
     [AC_MSG_FAILURE([machine dependent configuration file $host_frag \
does not exist!], 1)])

   AS_IF([test "$ARCH" = CRAY_XT3],
     [cross_compiling=yes])

   mh_setup=""
   AS_IF([test x"$with_setup" != x],
     [AC_MSG_CHECKING([for user specific setup])
      AC_MSG_RESULT($with_setup)
      mh_setup=$with_setup
   ])
   AC_GET_MH([$srcdir/$host_frag])])
dnl ######################################################################
dnl Checks for programs.
dnl ######################################################################
CHECK_GNU_MAKE
dnl we prefer GNU make because the System V make on many commercial
dnl platforms is slow
AC_CHECK_PROGS([MAKE], [$_cv_gnu_make_command gmake gnumake make])
dnl maintainers need GNU make for the regeneration rules, so enforce that
dnl
AS_IF([test "x${_cv_gnu_make_command}x$USE_MAINTAINER_MODE" = xxyes],
  [AC_MSG_FAILURE([maintainers need to use GNU make or compatible])])
AC_PROG_MAKE_SET
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MKDIR_P
m4_pushdef([_AC_PROG_CC_C89],
  [rm -rf conftest.dSYM # needed when configured on MacOS with CFLAGS='-g'
   _AC_PROG_CC_C99(,[AC_MSG_FAILURE([CC=$CC is not an ISO C99 compiler])])])dnl
AC_PROG_CC
m4_popdef([_AC_PROG_CC_C89])dnl
AC_PROG_CPP
AC_PROG_GREP
if test "X$FC" != Xno ; then
AC_PROG_FC
m4_define([test_with_fortran],[test -n "$FC" && test "X$FC" != Xno])dnl
dnl check that FC is working just as well as CC
AC_LANG_PUSH([Fortran])
_AC_COMPILER_EXEEXT
rm -rf conftest.dSYM # needed when configured on MacOS without FCFLAGS or with FCFLAGS='-g'
AC_LANG_POP([Fortran])
AC_FC_MODULE_OUTPUT_FLAG
AC_FC_SRCEXT([f90])
AC_FC_FREEFORM
AC_LANG_PUSH([Fortran])
ACX_PROG_FC_CHECK_FPP
ac_fcflags_srcext="${ac_fcflags_srcext+$ac_fcflags_srcext }$FC_FPP_FLAG"
FCFLAGS_f90="${FCFLAGS_f90+$FCFLAGS_f90 }$FC_FPP_FLAG"
AC_PROG_FPP
AC_LANG_POP([Fortran])
fi
AM_CONDITIONAL([with_fortran], [test_with_fortran])

AM_PROG_AS
AS_IF([test x"$USE_MAINTAINER_MODE" = xyes],
  [AC_PROG_GNU_M4])
# FIXME: is this still needed?
AS_CASE([$host],
  [sx*-*-*], NEC_ES_SX_HIGH_RESOLUTION_TIMER=yes,
  [es*-*-*], NEC_ES_SX_HIGH_RESOLUTION_TIMER=yes)
AM_CONDITIONAL([NEC_ES_SX_HIGH_RESOLUTION_TIMER],
  [test "x$NEC_ES_SX_HIGH_RESOLUTION_TIMER" = xyes])
AX_WITH_PERL([/bin/false], $PATH:$prefix/bin:/usr/bin:/usr/local/bin)
AX_PROG_PERL_VERSION([5.0005],[],
        AC_MSG_ERROR([Perl (version>=5.0005) is required]))
AX_WITH_PYTHON([/bin/false])
dnl ######################################################################
dnl use libtool
dnl before switching on libtool, identify compilers that prevent us from
dnl certain build configurations
AC_MSG_NOTICE([setting up libtool])
LT_PREREQ([2.4])
ACX_USE_LIBTOOL_CONFIGURATION
dnl put information on shared/static archives into ppm.settings
AC_SUBST([enable_shared])
AM_SUBST_NOTMAKE([enable_shared])
AC_SUBST([enable_static])
AM_SUBST_NOTMAKE([enable_static])

dnl ######################################################################
dnl Properties of the compiler
dnl ######################################################################
m4_define([PPM_ADD_FC_PP_DEFINITION],
  [PPM_FC_FEATURE_DEFS="$PPM_FC_FEATURE_DEFS
@%:@define $1[]m4_ifval([$2],[ $2])"])
AM_PROG_CC_C_O
if test_with_fortran ; then
AC_PROG_FC_C_O
ACX_SL_FC_MOD_SUFFIX(,[AC_MSG_ERROR(dnl
[Fortran module file extension must be known to install F90 modules.])])
AC_LANG_PUSH([Fortran])
ACX_SL_FC_CHECK_MOD_PATH_FLAG
AC_LANG_POP([Fortran])
ACX_FC_VOLATILE([PPM_ADD_FC_PP_DEFINITION([HAVE_VOLATILE])])
m4_pushdef([_AC_PROG_FC_V],m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_FC_V])),
  [-v -verbose --verbose -V],[-v -Wl,-v -verbose --verbose -V]))dnl
AC_FC_DUMMY_MAIN(,[:])
AC_FC_WRAPPERS
m4_popdef([_AC_PROG_FC_V])dnl

ACX_FIND_CFORTRAN_DEF
ACX_XLF_QEXTNAME_ADD_APPENDUS
ACX_CHECK_CFORTRAN
fi
AM_CONDITIONAL([FC_MOD_UPPERCASE],
  [test x$ac_fc_mod_uppercase = xyes])

dnl ######################################################################
dnl Checks for header files.
dnl ######################################################################
dnl
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h \
 unistd.h time.h sys/time.h sys/unistd.h \
 inttypes.h stdint.h sys/types.h setjmp.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_CHECK_TYPES([u_char])
dnl
dnl Checks for library functions.
AC_CHECK_FUNCS([gettimeofday sysconf uname])
dnl
dnl Check for byte ordering
AC_C_BIGENDIAN
dnl
dnl Check whether char is unsigned by default
AC_C_CHAR_UNSIGNED
AS_IF([test x$ac_cv_c_char_unsigned = xyes],
  [C_CHAR_IS_UNSIGNED=1],
  [C_CHAR_IS_UNSIGNED=0])
AC_SUBST([C_CHAR_IS_UNSIGNED])
dnl
dnl check for C datatype sizes
AC_CHECK_SIZEOF([char])
AC_CHECK_SIZEOF([short])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([float])
AC_CHECK_SIZEOF([double])
AC_CHECK_SIZEOF([long double])
AC_CHECK_SIZEOF([int *])
AC_CHECK_SIZEOF([size_t],[AC_INCLUDES_DEFAULT
@%:@include <stdlib.h>
])
dnl
dnl preprocessor implicit type conversions depends on the size of these
AC_SUBST([C_SIZEOF_INT],[$ac_cv_sizeof_int])
AC_SUBST([C_SIZEOF_LONG],[$ac_cv_sizeof_long])
AC_SUBST([C_SIZEOF_LONG_LONG],[$ac_cv_sizeof_long_long])
TJ_FIND_INTEGRAL_TYPE([int64_t],[DEFINED_INT64_T],[@%:@include <inttypes.h>])
AS_IF([test x"$DEFINED_INT64_T" = x"long long"],
  [AC_DEFINE([INT64_T_IS_LONG_LONG],[1],
     [defined if int64_t is long long in disguise])],
  [test x"$DEFINED_INT64_T" = x"long"],
  [AC_DEFINE([INT64_T_IS_LONG],[1],
     [defined if int64_t is long in disguise])])
dnl
if test_with_fortran ; then
ACX_FC_OPENMP_TASK([PPM_ADD_FC_PP_DEFINITION([HAVE_FC_OPENMP_TASK])])
dnl Check for Fortran types:
dnl first find wether we have iso_c_binding module
ACX_FORTRAN_CHECK_MODULE([ISO_C_BINDING], [integer(c_int) :: a])
dnl Defines FORT_INT1=1, FORT_INT2=1, FORT_INT4=1, FORT_INT8=1, FORT_INT16=1
dnl         FORT_REAL4=1, FORT_REAL8=1, FORT_REAL16=1
dnl         FORT_COMPLEX8=1, FORT_COMPLEX16=1, FORT_COMPLEX32=1
dnl         FORT_DOUBLE_COMPLEX=1 if the types are available.
ACX_FORTRAN_USUAL_TYPE_KINDS
dnl ACX_FORTRAN_C_INT
ACX_FORTRAN_C_INT64_T
ACX_FORTRAN_C_FLOAT
ACX_FORTRAN_TYPE_KIND([integer],,[void *],,,[PPM_ADDRESS_KIND=$acx_fortran_kind_subst])
fi
AC_SUBST([PPM_ADDRESS_KIND])

ACX_C_FLOAT_DIGITS

if test_with_fortran ; then
ACX_FC_REAL_EDIT_DESCRIPTORS([PPM])
AC_SUBST([PPM_de_g_sp])
AM_SUBST_NOTMAKE([PPM_de_g_sp])
AC_SUBST([PPM_de_g_sp_width])
AM_SUBST_NOTMAKE([PPM_de_g_sp_width])
AC_SUBST([PPM_de_g_dp])
AM_SUBST_NOTMAKE([PPM_de_g_dp])
AC_SUBST([PPM_de_g_dp_width])
AM_SUBST_NOTMAKE([PPM_de_g_dp_width])
fi

dnl Find size of system limits and Fortran intrinsic datatypes
ACX_C_LIMIT([PATH_MAX],,[AC_MSG_FAILURE([Cannot determine PATH_MAX.])])
ACX_C_CHAR_BITS
CHAR_BIT=$acx_cv_c_char_bits
AC_SUBST([CHAR_BIT])
if test_with_fortran ; then
ACX_FORTRAN_CHECK_SIZEOF_INTEGRAL_TYPE([INTEGER],,subst)
fi
dnl Set maximal alignment, this is needed in cross-host exchange of
dnl memory buffers via MPI_{Get|Put} to pre-compute the memory layout
dnl without "asking" the communication partner
m4_define([default_maximum_alignment],[64])
AC_ARG_WITH([maximum-alignment],
  [AS_HELP_STRING([--with-maximum-alignment=ALIGN],
     [Align all relevant addresses to multiples of ALIGN bytes @<:@default=]default_maximum_alignment[ unless __BIGGEST_ALIGNMENT__ is defined@:>@.])],,
  [with_maximum_alignment=default_maximum_alignment])
AS_IF([echo "$with_maximum_alignment" | grep '^[0-9][0-9]*$' >/dev/null \
   && test "$with_maximum_alignment" -gt 0],,
  [AC_CHECK_DECL([__BIGGEST_ALIGNMENT__],
     [AC_COMPUTE_INT([with_maximum_alignment],[__BIGGEST_ALIGNMENT__])],
     [with_maximum_alignment=default_maximum_alignment])])
AC_SUBST([PPM_MAXIMUM_ALIGNMENT],[$with_maximum_alignment])
AC_DEFINE_UNQUOTED([PPM_MAXIMUM_ALIGNMENT],[$with_maximum_alignment],
  [Align data of unknown alignment to multiples of this value.])
AC_MSG_NOTICE([Aligning unknown data to multiples of $with_maximum_alignment bytes.])

dnl make setjmp/longjmp available to Fortran
AC_CHECK_SIZEOF([jmp_buf],,[AC_INCLUDES_DEFAULT
@%:@ifdef HAVE_SETJMP_H
@%:@include <setjmp.h>
@%:@endif
])
JMP_BUF_ISIZE=`expr \( "$ac_cv_sizeof_jmp_buf" + "$ac_cv_sizeof_int" - 1 \) / "$ac_cv_sizeof_int"`
AC_SUBST([JMP_BUF_ISIZE])

dnl check if the compiler supports large bounds values
if test_with_fortran ; then
ACX_FC_LARGE_ARRAY_BOUNDS(,
  [AS_CASE(["$large_bounds_arg"],
     [yes\ assumed\ for\ cross\ compilation],
     [AC_MSG_NOTICE([If large array bounds work cannot be tested in cross compilation mode! You better have a good compiler])],
     [FCFLAGS="${FCFLAGS+$FCFLAGS }$large_bounds_arg"])
   AS_UNSET([large_bounds_arg])],
  [AC_MSG_NOTICE([large array bounds are unsupported!])],[large_bounds_arg],["yes assumed for cross compilation"])
dnl check if the compiler supports pointer bounds remapping
ACX_FC_POINTER_BOUNDS_REMAP([HAVE_POINTER_REMAP=yes],[HAVE_POINTER_REMAP=no])
AS_IF([test $HAVE_POINTER_REMAP = yes],
  [ACX_FC_POINTER_HUGE_BOUNDS_REMAP(,
     [AC_MSG_NOTICE([Disabling native low-bounds mapping, $FC has a known defective implementation.])
      HAVE_POINTER_REMAP=no])])
AC_SUBST([HAVE_POINTER_REMAP])
AM_SUBST_NOTMAKE([HAVE_POINTER_REMAP])
fi
AM_CONDITIONAL([HAVE_POINTER_REMAP],[test x"$HAVE_POINTER_REMAP" = xyes])

dnl ######################################################################
dnl Optional feature arguments
dnl ######################################################################
dnl
dnl Let the user choose an alternative math library
AS_IF([test x${PLATFORM_OPT_LIBM+set} = xset],
  [ACX_ALTERNATIVE_LIBM([[$PLATFORM_OPT_LIBM]])],
  [ACX_ALTERNATIVE_LIBM])

dnl checking for library routines and intrinsics that might require
dnl emulation
if test_with_fortran ; then
ACX_FORTRAN_CHECK_LIBFUNC([get_environment_variable],
[integer :: home_length
      character(len=10), parameter :: home_varname = 'HOME'
      call get_environment_variable(home_varname, length=home_length)],
[], [AC_LIBOBJ([get_environment_variable])])

ACX_F90_PACKAGE([ieee_arithmetic],[ieee_arithmetic],,
  [$IEEE_ARITHMETIC_FCFLAGS],
  [AC_MSG_WARN([Optional module ieee_arithmetic not found or not compilable.])
   have_ieee_arithmetic_fc_bindings=no],
  [ieee_is_normal],,[[$IEEE_FC_LIB]],,
  [AC_MSG_WARN([Cannot link programs using ieee_arithmetic ieee_is_normal.])
   have_ieee_arithmetic_fc_bindings=no],,
  [logical :: l],,[l = ieee_is_normal(0.0)])
AS_IF([test x"$have_ieee_arithmetic_fc_bindings" != xno],
  [PPM_ADD_FC_PP_DEFINITION([HAVE_IEEE_SIGNALING_NAN])],
  [USE_IEEE_ARITHMETIC=
   AC_CHECK_HEADER([math.h],,
     [AC_MSG_FAILURE([math.h header not found but required for ieee_arithmetic emulation.])])
   AC_CHECK_DECL([isnormal],[],
     [AC_MSG_FAILURE([isnormal macro not found but required for ieee_arithmetic emulation.])],[[@%:@include <math.h>]])
   ACX_C_SIGNALING_NAN([PPM_ADD_FC_PP_DEFINITION([HAVE_IEEE_SIGNALING_NAN])])
   ACX_C_IEEE_COPYSIGN
   AC_MSG_WARN([Emulating ieee_arithmetic module with C99 functionality.])
  ])
fi
AM_CONDITIONAL([GETENV_EMULATION],
  [test x$acx_fortran_routine_get_environment_variable = xno])
AM_CONDITIONAL([IEEE_ARITHMETIC_EMULATION],
  [test x$have_ieee_arithmetic_fc_bindings = xno])

dnl the implicit FPU precision of 80bits on x86 and DAZ/FTZ behaviour
dnl on SSE2-enabled systems are a potential problem, check if we can
dnl adjust for this
ACX_CHECK_KAHAN_SUMMATION_SETTINGS(
  [PPM_ADD_FC_PP_DEFINITION([NEED_PRECISION_CONTROL])],
  [PPM_ADD_FC_PP_DEFINITION([NEED_UNDERFLOW_CONTROL])])

if test_with_fortran ; then
AC_LANG_PUSH([Fortran])
m4_pushdef([fail_action],
  [PPM_ADD_FC_PP_DEFINITION([UNRELIABLE_DENORMAL])])
ACX_FC_CHECK_PRINT_DENORMAL_IFELSE(
  [ACX_LANG_CHECK_COMPARE_DENORMAL_IFELSE(
     [AC_LANG_PUSH([C])
      ACX_LANG_CHECK_COMPARE_DENORMAL_IFELSE(
        [ACX_C_STRTOD_DENORMAL(,
           [fail_action])],
        [fail_action])
      AC_LANG_POP([C])],
     [fail_action])],
  [fail_action])
m4_popdef([fail_action])
AC_LANG_POP([Fortran])
fi

dnl
dnl ######################################################################
dnl Checks for MPI.
dnl ######################################################################
have_MPI_c_bindings=no
have_MPI_fc_bindings=no
AC_ARG_ENABLE([MPI],
  [AS_HELP_STRING([--enable-MPI],
   [compile for and use MPI environment @<:@default: auto@:>@])],
   [AS_IF([test x"$enable_MPI" != xno],
      [force_enable_MPI=yes])])
AS_IF([test x"$enable_MPI" != xno],
  [enable_MPI=yes
   if test_with_fortran ; then
   ACX_F90_PACKAGE([MPI],[mpi],,,
     [AC_MSG_WARN([Optional mpi module not found or not compilable.])
      enable_MPI=no],[mpi_waitall],
     [mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
     [AC_MSG_WARN([Cannot link Fortran 90 MPI programs.])
      enable_MPI=no],,
[INTEGER :: cnt, req(2), stati(MPI_STATUS_SIZE, 2), ierror],
[(cnt, req, stati, ierror)])
   AS_IF([test $enable_MPI = no],
     [enable_MPI=yes
      ACX_FORTRAN_PACKAGE([MPI], [mpif.h],,,
        [AC_MSG_WARN([Required include mpif.h not found or not compilable.])
         enable_MPI=no],[mpi_waitall],
        [mpi mpi_f90 mpi_f77 mpich],[[-lmpi_f77 -lmpi],[-lmpi]],,
        [AC_MSG_WARN([Cannot link Fortran MPI programs.])
         enable_MPI=no])
      AS_IF([test $enable_MPI = yes],
        [FC_IMPORT_MPI="INCLUDE 'mpif.h'"])],
     [FC_IMPORT_MPI="USE mpi"
      MPI_FC_INCLUDE=$MPI_FC_MOD])
   fi
   ACX_C_PACKAGE([MPI],[mpi.h],,,
     [AC_MSG_WARN([Required header mpi.h not found or not compilable.])
      enable_MPI=no],[MPI_Waitall],[mpi mpich],,,
     [AC_MSG_WARN([Cannot link C MPI programs.])
      enable_MPI=no])
  ])
enable_mpi_type_dup_rma=$enable_MPI
save_CFLAGS=$CFLAGS
save_FCFLAGS=$FCFLAGS
save_LIBS=$LIBS
CFLAGS="$CFLAGS $MPI_C_INCLUDE"
FCFLAGS="$FCFLAGS $MPI_FC_INCLUDE"
LIBS="$LIBS $MPI_C_LIB"
AC_ARG_WITH([regard-for-quality],
  [AS_HELP_STRING([--without-regard-for-quality],
     [ignore results of checks for known defects @<:@default: abort on error@:>@])],
  [AS_IF([test x"$withval" != xno],
     [with_regard_for_quality=yes])],[with_regard_for_quality=yes])
AS_IF([test $enable_MPI = yes],
  [ACX_MPIRUN(,[saved_MPI_LAUNCH=$MPI_LAUNCH],
     [AC_MSG_WARN([$MPI_LAUNCH_failMsg])
      AS_VAR_IF([with_regard_for_quality], [yes],
        [AC_MSG_FAILURE([unable to find a working MPI launch program, which is
required for checks for known MPI defects (see --without-regard-for-quality)])],
        [saved_MPI_LAUNCH=$MPI_LAUNCH ; MPI_LAUNCH=true])])])
AC_ARG_VAR([NOMPI_LAUNCH],
  [optional command prefix for running test programs not using MPI
This will most often be empty, but can be used to e.g. forward commands to
another cluster node, adjust CPU frequency governor, or bind the task.])
dnl do not remove conftest.err we need to inspect the output
m4_pushdef([_AC_RUN_LOG_LIMIT],
  m4_bpatsubst(m4_dquote(m4_defn([_AC_RUN_LOG_LIMIT])),
    [rm -f conftest.er1 conftest.err],[rm -f conftest.er1]))dnl
m4_pushdef([_AC_LINK_IFELSE],
  m4_bpatsubst(m4_dquote(m4_defn([_AC_LINK_IFELSE])),
    [rm -f core conftest.err],[rm -f core]))dnl
AS_IF([test $enable_MPI = yes],
  [ACX_MPI_DEFECTS(,,
     [ACX_MPI_DEFECTS_DOCUMENT
      AS_CASE([$acx_subtestname],[openmpi-typedup-rma],
        [AC_MSG_WARN([known OpenMPI pre-1.7.x problem detected that might affect distributed array!])],
        [incorrect_mpi_abort_exitcode],
        [AC_MSG_WARN([Your installation drops the exit code of MPI_Abort!
Consider installing a better MPI/batch scheduler.])
         AS_IF([grep -v -i 'mpi_abort(@<:@^,()@:>@*, 2)' conftest.err >/dev/null],
           [AC_MSG_NOTICE([Implementing fragile work-around for tests!])
            AC_DEFINE([XT_NEED_MPI_ABORT_WORK_AROUND],[1],
              [To pass the value for MPI_Abort to the surrounding execution environment, a hack is needed.])])],
        [AS_VAR_IF([with_regard_for_quality],[yes],
           [AC_MSG_FAILURE([test for known defect $acx_subtestname failed,
re-configure with --without-regard-for-quality in case you must use the
given MPI implementation, but some tests might fail])],
           [AC_MSG_NOTICE([ignoring failed test for known defect $acx_subtestname,
because configure option --without-regard-for-quality is active])])])])])
m4_popdef([_AC_RUN_LOG_LIMIT])m4_popdef([_AC_LINK_IFELSE])dnl
rm -f conftest.err
AS_IF([test $enable_MPI = no],
  [MPI_C_INCLUDE= MPI_C_LIB= MPI_FC_LIB= MPI_FC_INCLUDE= MPI_FINT_FC_KIND=-1],
  [AC_DEFINE([USE_MPI],[1],[use MPI parallel process environment])
   if test_with_fortran ; then
   AC_LANG_PUSH([Fortran])
   $MKDIR_P include/f77
   AS_IF([test x"$FC_IMPORT_MPI" = x"USE mpi"],
     [acx_temp="      $FC_IMPORT_MPI
      IMPLICIT NONE"],
     [acx_temp="      IMPLICIT NONE
           $FC_IMPORT_MPI"])
   AS_FOR([MPI_INT2TEST],[mpi_int2test],
     [MPI_INTEGER4 MPI_INTEGER8 MPI_REAL4 MPI_REAL8],
     [AC_MSG_CHECKING([if Fortran MPI bindings define $mpi_int2test])
      AC_COMPILE_IFELSE(AC_LANG_PROGRAM([$acx_temp],
[      INTEGER :: i
      i = MPI_INT2TEST]),[AC_MSG_RESULT([yes])
echo '@%:@define HAVE_'"$mpi_int2test" >>include/f77/mpi_fc_conf.inc],
     [AC_MSG_RESULT([no])])])
   AC_LANG_POP([Fortran])
   TJ_FIND_INTEGRAL_TYPE([MPI_Fint],,[@%:@include <mpi.h>])
   acx_temp=AS_TR_SH([c_$tj_cv_c_type_MPI_Fint])
   ACX_FORTRAN_TYPE_KIND([integer],[$acx_temp],[MPI_Fint],
    [      use iso_c_binding],[@%:@include <mpi.h>],
    [MPI_FINT_FC_KIND=$acx_fortran_kind_subst],
    [MPI_FINT_FC_KIND=$acx_fortran_kind_subst])
   fi])
MPI_LAUNCH=$saved_MPI_LAUNCH
CFLAGS=$save_CFLAGS
FCFLAGS=$save_FCFLAGS
LIBS=$save_LIBS
AS_IF([test $enable_MPI != yes -a x"$force_enable_MPI" = xyes],
  [AC_MSG_ERROR([user-requested MPI support unavailable.])])
AM_CONDITIONAL([USE_MPI],
  [test $enable_MPI = yes])
AM_CONDITIONAL([USE_MPI_MOD],
  [test x"$FC_IMPORT_MPI" = x"USE mpi"])
AC_SUBST([USE_MPI],[$enable_MPI])
AM_SUBST_NOTMAKE([USE_MPI])
AM_CONDITIONAL([ENABLE_DISTRIBUTED_ARRAY_C],
  [test x"$enable_mpi_type_dup_rma" = xyes])
AC_SUBST([ENABLE_DISTRIBUTED_ARRAY_C],[$enable_mpi_type_dup_rma])
AC_SUBST([HAVE_MPI_C_BINDINGS],["$have_MPI_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_MPI_C_BINDINGS])
AC_SUBST([HAVE_MPI_FC_BINDINGS],["$have_MPI_fc_bindings"])
AM_SUBST_NOTMAKE([HAVE_MPI_fc_BINDINGS])
AC_SUBST([MPI_FINT_FC_KIND])
AM_SUBST_NOTMAKE([MPI_FINT_FC_KIND])

dnl ######################################################################
dnl Checks for NETCDF.
dnl ######################################################################
dnl   NetCDF might use HDF5 on some platforms, which in turn might
dnl   require szip and almost definitely requires zlib
if test_with_fortran ; then
AC_ARG_ENABLE([netcdf],
  [AS_HELP_STRING([--disable-netcdf],
     [disable checks for netCDF library @<:@default: auto@:>@])],,
  [enable_netcdf=auto])

AC_ARG_ENABLE([hdf5],
  [AS_HELP_STRING([--disable-hdf5],
     [disable checks for HDF5 library @<:@default: auto@:>@])],,
  [enable_hdf5=auto])

AS_IF([test x"$enable_hdf5" != xno && test x"$enable_netcdf" != xno],
  [dnl
dnl ----------------------------------------------------------------------
dnl   test for ZLIB compression library
dnl
   ACX_C_PACKAGE([zlib],[zlib.h],,,
     [AC_MSG_NOTICE([zlib library header not found, zlib compression not supported])],[deflate],[z],,,
     [AC_MSG_NOTICE([cannot link zlib library, zlib compression not supported])])
   AS_VAR_IF([have_zlib_c_bindings], [no],
     [ASX_VAR_UNSET([ZLIB_C_INCLUDE])
      ASX_VAR_UNSET([ZLIB_C_LIB])])
dnl ----------------------------------------------------------------------
dnl   test for SZIP compression library
dnl
   ACX_C_PACKAGE([szip],[szlib.h],[@%:@ifdef STDC_HEADERS
@%:@include <stdlib.h>
@%:@endif],,
     [AC_MSG_NOTICE([szip library header not found, szip compression not supported])],[SZ_BufftoBuffCompress],[sz],,,
     [AC_MSG_NOTICE([cannot link szip library, szip compression not generally supported])])
   AS_VAR_IF([have_szip_c_bindings], [no],
     [ASX_VAR_UNSET([SZIP_C_INCLUDE])
      ASX_VAR_UNSET([SZIP_C_LIB])])
dnl ----------------------------------------------------------------------
dnl   test for HDF5 file format library, netCDF 4.x usually requires it
dnl
   ACX_C_PACKAGE([hdf5],[hdf5.h],,ACX_M4_GENERATE_SUBSETS([[$SZIP_C_INCLUDE],[$ZLIB_C_INCLUDE],[$MPI_C_INCLUDE]],[ ]),
     [AC_MSG_NOTICE([HDF5 library header not found, might not support netCDF 4.x])],[H5open],[hdf5],ACX_M4_GENERATE_SUBSETS([[$SZIP_C_LIB],[$ZLIB_C_LIB],[$LIBM]],[ ]),,
     [AC_MSG_NOTICE([cannot link HDF5 library])])
   AS_VAR_IF([have_hdf5_c_bindings], [yes],
     [AC_DEFINE([HAVE_LIBHDF5],[1],[HDF 5 library is present if defined to 1])],
     [ASX_VAR_UNSET([HDF5_C_INCLUDE])
      ASX_VAR_UNSET([HDF5_C_LIB])
      AS_VAR_IF([enable_hdf5], [yes],
        [AC_MSG_FAILURE([Requested package HDF5 unavailable.])])])],
  [have_zlib_c_bindings=no
   ASX_VAR_UNSET([ZLIB_C_INCLUDE])
   ASX_VAR_UNSET([ZLIB_C_LIB])
   have_szip_c_bindings=no
   ASX_VAR_UNSET([SZIP_C_INCLUDE])
   ASX_VAR_UNSET([SZIP_C_LIB])
   have_hdf5_c_bindings=no
   ASX_VAR_UNSET([HDF5_C_INCLUDE])
   ASX_VAR_UNSET([HDF5_C_LIB])])

AC_SUBST([HAVE_ZLIB_C_BINDINGS],["$have_zlib_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_ZLIB_C_BINDINGS])
AC_SUBST([HAVE_SZIP_C_BINDINGS],["$have_szip_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_SZIP_C_BINDINGS])
AC_SUBST([HAVE_HDF5_C_BINDINGS],["$have_hdf5_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_HDF5_C_BINDINGS])

AS_IF([test x"$enable_netcdf" != xno],
  [dnl
dnl ----------------------------------------------------------------------
dnl   test for curl library, possibly needed by netCDF DAP
dnl
   ACX_C_PACKAGE([curl],[curl/curl.h],,,
     [AC_MSG_NOTICE([curl library header not found, netCDF remote protocols might be unsupported])],
     [curl_easy_init],[curl],,,
     [AC_MSG_NOTICE([cannot link curl library, netCDF remote protocols might be unsupported])])
   AS_VAR_IF([have_curl_c_bindings], [no],
     [ASX_VAR_UNSET([CURL_C_INCLUDE])
      ASX_VAR_UNSET([CURL_C_LIB])])
dnl ----------------------------------------------------------------------
dnl   test for netCDF file format library
dnl
   ACX_F90_PACKAGE([NetCDF],[netcdf],,,
     [AC_MSG_WARN([Optional netcdf module not found or not compilable.])],
     [nf90_close],
     [netcdf netcdff],ACX_M4_GENERATE_SUBSETS([[-lnetcdf],[$CURL_C_LIB],[-lhdf5_hl],[$HDF5_C_LIB]],[ ]),,
     [AC_MSG_WARN([Cannot link Fortran 90 netcdf programs.])],,
     [INTEGER :: file_ncid, error_id],,
     [      error_id = nf90_close(file_ncid)])
   AS_VAR_IF([have_NetCDF_fc_bindings], [no],
     [ACX_FORTRAN_PACKAGE([NetCDF],[netcdf.inc],,,
       [AC_MSG_NOTICE([Include netcdf.inc not found or not compilable.])],
       [nf_open],[netcdf netcdff],ACX_M4_GENERATE_SUBSETS([[-lnetcdf],[$CURL_C_LIB],[-lhdf5_hl],[$HDF5_C_LIB]],[ ]),,
       [AC_MSG_NOTICE([Cannot link with NetCDF library.])],[])
      AS_VAR_IF([have_NetCDF_fc_bindings], [yes],
        [FC_IMPORT_NETCDF="INCLUDE 'netcdf.inc'"])],
     [FC_IMPORT_NETCDF='USE netcdf' NETCDF_FC_INCLUDE=$NETCDF_FC_MOD])
   AS_VAR_IF([have_NetCDF_fc_bindings], [yes],
     [AC_DEFINE([HAVE_FORTRAN_NETCDF],[1],[netCDF Fortran bindings present if defined to 1])],
     [ASX_VAR_UNSET([NETCDF_FC_INCLUDE])
      ASX_VAR_UNSET([NETCDF_FC_LIB])
      AS_VAR_IF([enable_netcdf], [yes],
        [AC_MSG_FAILURE([Requested package NetCDF unavailable.])])])],
  [have_curl_c_bindings=no
   ASX_VAR_UNSET([CURL_C_INCLUDE])
   ASX_VAR_UNSET([CURL_C_LIB])
   have_NetCDF_fc_bindings=no
   ASX_VAR_UNSET([NETCDF_FC_INCLUDE])
   ASX_VAR_UNSET([NETCDF_FC_LIB])])
else
  have_curl_c_bindings=no
  ASX_VAR_UNSET([CURL_C_INCLUDE])
  ASX_VAR_UNSET([CURL_C_LIB])
  have_NetCDF_fc_bindings=no
  ASX_VAR_UNSET([NETCDF_FC_INCLUDE])
  ASX_VAR_UNSET([NETCDF_FC_LIB])
fi

AC_SUBST([HAVE_CURL_C_BINDINGS],["$have_curl_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_CURL_C_BINDINGS])
AM_CONDITIONAL([USE_FC_NETCDF],[test x"$have_NetCDF_fc_bindings" = xyes])
AM_CONDITIONAL([USE_FC_NETCDF_MOD],[test "x$FC_IMPORT_NETCDF" = 'xUSE netcdf'])
AC_SUBST([HAVE_NETCDF_FC_BINDINGS],["$have_NetCDF_fc_bindings"])
AM_SUBST_NOTMAKE([HAVE_NETCDF_FC_BINDINGS])
dnl ######################################################################
dnl Checks for ParMetis
dnl ######################################################################
AC_ARG_ENABLE([parmetis],
  [AS_HELP_STRING([--enable-parmetis],
   [Provide wrapper for ParMETIS graph partitioners @<:@default: yes@:>@])])
AS_IF([test $enable_MPI = yes -a x"$enable_parmetis" != xno],
  [enable_parmetis=yes
   ACX_C_PACKAGE([parmetis],[parmetis.h],[@%:@include <mpi.h>],
     [[],[$MPI_C_INCLUDE]],
     [AC_MSG_WARN([Header for package ParMETIS not found.])
      enable_parmetis=no],
[ParMETIS_V3_PartKway],[parmetis],ACX_M4_GENERATE_SUBSETS([[-lmetis],[$MPI_C_LIB],[$LIBM]],[ ]),,
     [AC_MSG_WARN([Cannot link with package ParMETIS library.])
      enable_parmetis=no])],
     [have_parmetis_c_bindings=no])
save_CFLAGS=$CFLAGS
save_LIBS=$LIBS
CFLAGS="$MPI_C_INCLUDE $PARMETIS_C_INCLUDE $CFLAGS"
LIBS="$PARMETIS_C_LIB $MPI_C_LIB $LIBS"
PARMETIS_MAJOR_VERSION=-1
AS_IF([test x"$enable_parmetis" = xyes],
  [AC_COMPUTE_INT([PARMETIS_MAJOR_VERSION],[PARMETIS_MAJOR_VERSION],
     [@%:@include <parmetis.h>],
     [AC_MSG_FAILURE([Cannot detect ParMetis version!])])])
AM_CONDITIONAL([USE_PARMETIS], [test x"$enable_parmetis" = xyes])
AM_CONDITIONAL([HAVE_PARMETIS_V3], [test $PARMETIS_MAJOR_VERSION -eq 3])
AC_SUBST([HAVE_PARMETIS_C_BINDINGS],["$have_parmetis_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_PARMETIS_C_BINDINGS])
dnl
AC_ARG_ENABLE([metis],
  [AS_HELP_STRING([--enable-metis],
     [Provide wrapper for serial graph partitioners from METIS library @<:@default: auto@:>@])])
CFLAGS=$save_CFLAGS
LIBS=$save_LIBS
AS_IF([test x"$enable_metis" != xno],
  [enable_metis=yes
   AS_IF([test x"$enable_parmetis" = xyes -a $PARMETIS_MAJOR_VERSION -eq 3],
# search metis/metis.h if parmetis version 3 is also used
     [METIS_HEADER='metis/metis.h'],
     [METIS_HEADER='metis.h'])
   ACX_C_PACKAGE([metis],[$METIS_HEADER],,ACX_M4_GENERATE_SUBSETS([[$PARMETIS_C_INCLUDE],[$MPI_C_INCLUDE]],[ ]),
     [AC_MSG_WARN([Header for package METIS not found.])
      enable_metis=no],
     [METIS_PartGraphKway],[metis],[[$PARMETIS_C_LIB],[],[$LIBM]],,
     [AC_MSG_WARN([Cannot link with METIS library.])
      enable_metis=no])
  ],[have_metis_c_bindings=no])
CFLAGS="$MPI_C_INCLUDE $PARMETIS_C_INCLUDE $METIS_C_INCLUDE $save_CFLAGS"
LIBS="$PARMETIS_C_LIB $METIS_C_LIB $MPI_C_LIB $save_LIBS"
METIS_VER_MAJOR=-1
AS_IF([test x"$enable_metis" = xyes],
  [AC_CHECK_DECL([METIS_VER_MAJOR],
     [AC_COMPUTE_INT([METIS_VER_MAJOR],[METIS_VER_MAJOR],
       [@%:@include <$METIS_HEADER>],
       [AC_MSG_FAILURE([Cannot detect Metis version!])])],
     [AC_CHECK_FUNC([METIS_mCPartGraphKway],[METIS_VER_MAJOR=4],
        [AC_MSG_FAILURE([Cannot detect Metis version!])])
      AC_DEFINE_UNQUOTED([METIS_VER_MAJOR],[$METIS_VER_MAJOR],
        [Older versions of Metis don't define this themselves.])],
     [@%:@include <$METIS_HEADER>])],
  [METIS_C_INCLUDE= ; METIS_C_LIB=])
AM_CONDITIONAL([USE_METIS], [test x"$enable_metis" = xyes])
AM_CONDITIONAL([HAVE_METIS_V4], [test $METIS_VER_MAJOR -eq 4])
AC_SUBST([HAVE_METIS_C_BINDINGS],["$have_metis_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_METIS_C_BINDINGS])
dnl
dnl determine the exact type used by parmetis/metis to represent node
dnl indices and weights
PARMETIS_C_REAL_T=float
AS_IF([test x"$enable_parmetis" = xyes],
  [acx_idxtype_include="@%:@include <mpi.h>
@%:@include <parmetis.h>
@%:@if METIS_VER_MAJOR >= 5
typedef idx_t idxtype;
@%:@endif"
   TJ_FIND_INTEGRAL_TYPE([idxtype],[PARMETIS_C_IDXTYPE],[$acx_idxtype_include])
   AC_CHECK_TYPE([real_t],
     [TJ_FIND_TYPE([real_t], [PARMETIS_C_REAL_T],
        [$acx_idxtype_include], [float double])],,[$acx_idxtype_include])])

AS_IF([test x"$enable_metis" = xyes],
  [acx_metis_idxtype_include="@%:@include <$METIS_HEADER>
@%:@if METIS_VER_MAJOR >= 5
typedef idx_t idxtype;
@%:@endif"
   ASX_VAR_UNSET([ac_cv_type_real_t tj_cv_c_type_idxtype])
   TJ_FIND_INTEGRAL_TYPE([idxtype],[METIS_C_IDXTYPE],
     [$acx_metis_idxtype_include])
dnl next determine corresponding Fortran type kind
   AS_IF([test x"$enable_parmetis" != xyes],
     [acx_idxtype_include=$acx_metis_idxtype_include])
   AC_CHECK_TYPE([real_t],
     [TJ_FIND_TYPE([real_t], [METIS_C_REAL_T], [@%:@include <metis.h>],
        [float double])],[METIS_C_REAL_T=float],[$acx_metis_idxtype_include])],
  [METIS_C_REAL_T=$PARMETIS_C_REAL_T
   METIS_C_IDXTYPE=$PARMETIS_C_IDXTYPE])

# we require compatible types for METIS and ParMETIS
AS_IF([test x"$enable_parmetis" = xyes -a x"$enable_metis" = xyes],
  [AS_IF([test "$METIS_C_IDXTYPE" != "$PARMETIS_C_IDXTYPE" \
          || test "$PARMETIS_C_REAL_T" != "$METIS_C_REAL_T"],
     [_AS_ECHO_LOG([METIS_C_IDXTYPE=$METIS_C_IDXTYPE PARMETIS_C_IDXTYPE=$PARMETIS_C_IDXTYPE PARMETIS_C_REAL_T=$PARMETIS_C_REAL_T METIS_C_REAL_T=$METIS_C_REAL_T])
      AC_MSG_FAILURE([Must use compatible versions of METIS and ParMETIS])])])

AS_IF([test_with_fortran && { test x"$enable_parmetis" = xyes || test x"$enable_metis" = xyes; }],
  [acx_temp=AS_TR_SH([c_$METIS_C_IDXTYPE])
   ACX_FORTRAN_TYPE_KIND([integer],[$acx_temp],[idxtype],
    [      use iso_c_binding],,[METIS_FC_IDXTYPE_KIND=$acx_fortran_kind_subst],
    [METIS_FC_IDXTYPE_KIND=$acx_fortran_kind_subst])
   AS_IF([test x${METIS_FC_IDXTYPE_KIND+set} != xset],
     [AC_MSG_FAILURE([Cannot determine type kind of ParMETIS index type.])])
   acx_temp=AS_TR_SH([c_$METIS_C_REAL_T])
   ACX_FORTRAN_TYPE_KIND([real],[$acx_temp],[$METIS_C_REAL_T],
     [      use iso_c_binding],,[METIS_FC_REAL_T_KIND=$acx_fortran_kind_subst],
     [METIS_FC_REAL_T_KIND=$acx_fortran_kind_subst])
   AS_IF([test x${METIS_FC_REAL_T_KIND+set} != xset],
     [AC_MSG_FAILURE([Cannot determine type kind of ParMETIS real type.])])],
  [METIS_FC_IDXTYPE_KIND=-1
   METIS_FC_REAL_T_KIND=-1])
AC_SUBST([METIS_FC_IDXTYPE_KIND])
AM_SUBST_NOTMAKE([METIS_FC_IDXTYPE_KIND])
AC_SUBST([METIS_FC_REAL_T_KIND])
AM_SUBST_NOTMAKE([METIS_FC_REAL_T_KIND])
METIS_IDX_KIND_MATCH_I4=0
METIS_IDX_KIND_MATCH_I8=0
AS_IF([test_with_fortran && { test x"$enable_parmetis" = xyes || test x"$enable_metis" = xyes; }],
  [AC_LANG_PUSH([Fortran])
   AC_MSG_CHECKING([whether $METIS_FC_IDXTYPE_KIND matches i4 or i8])
   AC_COMPILE_IFELSE(
     [ACX_LANG_CHECK_INTEGER_EXPRESSION_SOURCE(
       [integer, parameter :: i4=SELECTED_INT_KIND(9)], [i4],
       [$METIS_FC_IDXTYPE_KIND], [use iso_c_binding])],
     [AC_COMPILE_IFELSE(
       [ACX_LANG_CHECK_INTEGER_EXPRESSION_SOURCE(
          [integer, parameter :: i8=SELECTED_INT_KIND(14)], [i8],
          [$METIS_FC_IDXTYPE_KIND], [use iso_c_binding])],
       [AC_MSG_FAILURE([Metis index type does not match i4 or i8])],
       [METIS_IDX_KIND_MATCH_I8=1
        AC_MSG_RESULT([i8])])],
     [METIS_IDX_KIND_MATCH_I4=1
        AC_MSG_RESULT([i4])])
   AC_LANG_POP([Fortran])])
AC_SUBST([METIS_IDX_KIND_MATCH_I4])
AM_SUBST_NOTMAKE([METIS_IDX_KIND_MATCH_I4])
AC_SUBST([METIS_IDX_KIND_MATCH_I8])
AM_SUBST_NOTMAKE([METIS_IDX_KIND_MATCH_I8])
CFLAGS=$save_CFLAGS
LIBS=$save_LIBS

dnl
dnl adjust library paths for Fortran compiler
have_parmetis_fc_bindings=no
AS_IF([test_with_fortran && test x"$enable_parmetis" = xyes],
  [ACX_FC_C_LINK([PARMETIS_C_LIB],[PARMETIS_FC_LIB],[AC_LANG_PROGRAM(,
[      CALL ParMETIS_V3_PartKway])])
   AC_SUBST([PARMETIS_FC_LIB])
   have_parmetis_fc_bindings=yes
   ])
AC_SUBST([HAVE_PARMETIS_FC_BINDINGS],["$have_parmetis_fc_bindings"])
AM_SUBST_NOTMAKE([HAVE_PARMETIS_FC_BINDINGS])
dnl
have_metis_fc_bindings=no
AS_IF([test_with_fortran && test x"$enable_metis" = xyes],
  [ACX_FC_C_LINK([METIS_C_LIB],[METIS_FC_LIB],[AC_LANG_PROGRAM(,
[      CALL METIS_PartGraphKway])])
   AC_SUBST([METIS_FC_LIB])
   have_metis_fc_bindings=yes])
AC_SUBST([HAVE_METIS_FC_BINDINGS],["$have_metis_fc_bindings"])
AM_SUBST_NOTMAKE([HAVE_METIS_FC_BINDINGS])
dnl
dnl ######################################################################
dnl ######################################################################

dnl Check for cryptographic library
AC_ARG_ENABLE([crypto],
  [AS_HELP_STRING([--disable-crypto],
   [wrap SSL crypto algorithms @<:@default: auto@:>@])],,[enable_crypto=auto])
AS_IF([test x"$enable_crypto" != xno],
  [ACX_C_PACKAGE([crypto],[openssl/evp.h],,,
     [AC_MSG_NOTICE([crypto evp library header not found, SSL digests not supported])],
     [EVP_Digest],[crypto],,,
     [AC_MSG_NOTICE([cannot link crypto library, SSL digests not supported])])
   AS_IF([test $have_crypto_c_bindings = yes],
     [if test_with_fortran ; then
      ACX_FC_C_LINK([CRYPTO_C_LIB],[CRYPTO_FC_LIB])
      AC_SUBST([CRYPTO_FC_LIB])
      fi
      AC_DEFINE([USE_CRYPTO],,[SSL crypto functions are available])],
     [ASX_VAR_UNSET([CRYPTO_C_INCLUDE])
      ASX_VAR_UNSET([CRYPTO_C_LIB])])
  ],[have_crypto_c_bindings=no])
dnl
AM_CONDITIONAL([USE_CRYPTO],[test $have_crypto_c_bindings = yes])
AC_SUBST([HAVE_CRYPTO_C_BINDINGS],["$have_crypto_c_bindings"])
AM_SUBST_NOTMAKE([HAVE_CRYPTO_C_BINDINGS])
AC_SUBST([HAVE_HASH_MD5],[yes])
AM_SUBST_NOTMAKE([HAVE_HASH_MD5])
AC_SUBST([HAVE_HASH_SHA1],
  [`test x"$have_crypto_c_bindings" = xyes && echo yes || echo no`])
AM_SUBST_NOTMAKE([HAVE_HASH_SHA1])
dnl
dnl ######################################################################
dnl ######################################################################
dnl
AC_ARG_ENABLE([low-mem],
  [AS_HELP_STRING([--enable-low-mem],
   [save memory by reducing precomputation tables @<:@default: off@:>@])])
AS_CASE(["$enable_low_mem"],
  [yes|on], [:],
  [off|no|''], [:],
  [*,*|bitvect], [flags="$enable_low_mem"
     while test -n "$flags" ; do
       flag=`echo "$flags" | sed 's/,.*//'`
       flags=`echo "$flags" | sed 's/@<:@^,@:>@*\(,\|$\)//'`
       AS_CASE(["$flag"],
         [bitvect], [:],
         [*], [AC_MSG_ERROR([Unrecognized option passed to --enable-low-mem])])
     done],
  [*], [AC_MSG_ERROR([Unrecognized option passed to --enable-low-mem])])
dnl
dnl
dnl add flag to disable automatic creation of test programs
AC_ARG_WITH([on-demand-check-programs],
  [AS_HELP_STRING([--with-on-demand-check-programs],
     [only build test programs when calling make check (default: build during make all)])],
  [AS_CASE([$withval],
     [yes|no],,
     [n],[with_on_demand_check_programs=no],
     [y],[with_on_demand_check_programs=yes],
     [AC_MSG_FAILURE([unexpected value for --with-on-demand-check-programs: $withval])])],
  [with_on_demand_check_programs=no])
AM_CONDITIONAL([with_on_demand_check_programs],
  [test x$with_on_demand_check_programs = xyes])
dnl
dnl
dnl add flag to disable automatic creation of example programs
AC_ARG_WITH([example-programs],
  [AS_HELP_STRING([--without-example-programs],
     [only build example programs when calling make examples in examples directory (default: build during make all)])],
  [AS_CASE([$withval],
     [yes|no],,
     [n],[with_example_programs=no],
     [y],[with_example_programs=yes],
     [AC_MSG_FAILURE([unexpected value for --with-example-programs: '$withval'])])],
  [with_example_programs=yes])
AM_CONDITIONAL([with_example_programs],
  [test x$with_example_programs = xyes])
dnl
dnl
dnl
AR="${AR:-ar}"
ARFLAGS="${ARFLAGS:-cru}"
AC_SUBST([AR])
AC_SUBST([ARFLAGS])

AS_IF([test x${FCLDFLAGS+set} != xset],
  [FCLDFLAGS=$FCFLAGS])
m4_foreach([flags_var],[[CFLAGS],[FCFLAGS],[LDFLAGS],[FCLDFLAGS]],
  [AC_ARG_VAR([BUILD_]flags_var,
    [append to ]flags_var[ during build but not in configure phase])
   AC_CONFIG_COMMANDS_PRE(flags_var[="$]flags_var[${BUILD_]flags_var[+ $BUILD_]flags_var[}"])])


# translate cpp DEFS into corresponding directives for the Fortran
# compiler to pass to the preprocessor
AC_CONFIG_COMMANDS_PRE([FCDEFS=`for ac_def in $DEFS ; do echo "$ac_def" ; done |
  sed 's/^-D/'"$FPP_DEFOPT"'/'`])
AC_SUBST([FCDEFS])

AC_SUBST([FCFLAGS])
AC_SUBST([FCLDFLAGS])

AC_CONFIG_FILES([Makefile include/Makefile include/f77/ppm.inc \
   include/f77/ftype_size.inc src/core/ppm_real_sp_dp_edit_descriptor.inc \
   test/atlocal])

AS_IF([test x = x"$PERL"],
  [AC_MSG_ERROR([PERL is undefined!])])

AM_MISSING_PROG([AUTOM4TE], [autom4te])

ACX_CHECK_DOXYGEN_LAYOUT
ACX_CHECK_DOXYGEN_HTML_EXTRA_STYLESHEET
AM_CONDITIONAL([HAVE_DOXYGEN_HTML_EXTRA_STYLESHEET],
  [test x"$ACX_DOXYFILE_HTML_EXTRA_STYLESHEET" = xyes])
AC_CONFIG_FILES([doc/unitdoc/scales-ppm.doxygen \
  src/pkgconfig/scales-ppm-core.pc ppm.settings \
  doc/Makefile doc/unitdoc/Makefile])
AS_IF([test_with_fortran],
  [AC_CONFIG_FILES([src/pkgconfig/scales-ppm.pc])])
AC_CONFIG_FILES([example/Graphvis.py],
  [chmod a+x "$ac_file"])
m4_define([f90makefiles],
[src/Makefile \
 include/f90/Makefile \
 include/f77/Makefile \
 test/Makefile \
 example/Makefile \
])dnl
AC_CONFIG_FILES(f90makefiles,
  [sed '
/FCCOMPILE =/s/(FCFLAGS)/(FCFLAGS) \$(FCDEFS)/
/mode=link/{
s/(FCFLAGS)/(FCLDFLAGS)/g
n
s/(FCFLAGS)/(FCLDFLAGS)/g
}
/@<:@	 @:>@\$(FC) \$(AM_FCFLAGS)/s/(FCFLAGS)/(FCFLAGS) \$(FCDEFS)/' \
   $ac_file  \
    >conftest.mk.fclinkrewrite;
   cmp $ac_file conftest.mk.fclinkrewrite >/dev/null \
    || mv conftest.mk.fclinkrewrite $ac_file; \
    rm -f conftest.mk.fclinkrewrite;
   $MKDIR_P `dirname $ac_file`/.deps
  ])
AC_SUBST([PPM_FC_FEATURE_DEFS])
AM_SUBST_NOTMAKE([PPM_FC_FEATURE_DEFS])
AC_CONFIG_FILES([include/f77/fc_feature_defs.inc])

AC_CONFIG_TESTDIR([test])
# edit out the am--include-marker that's required for automake 1.16
# and later on older automake versions
AS_VERSION_COMPARE([$am__api_version],[1.16],[extrasub="${extrasub+$extrasub$as_nl}s/\# am--include-marker\$//"
])

dnl Create dummy Fortran dependency files for Automake version 1.16 or later,
dnl which don't do so automatically anymore, older versions
dnl created the files by the time this command is executed. Having
dnl created the dummy Fortran dependency files, we enforce their
dnl regeneration at the build time by updating the timestamps of the
dnl corresponding Makefiles, which the Fortran dependency file have as
dnl prerequisites.
# force make to regenerate FC.deps file (because it depends on Makefile)
AC_CONFIG_COMMANDS([force_dep_regeneration],
  [if test x"$with_fortran" = xyes ; then
$MKDIR_P "test/${DEPDIR}" "src/${DEPDIR}" "example/${DEPDIR}"
touch "test/${DEPDIR}/FC.deps" "src/${DEPDIR}/FC.deps" "example/${DEPDIR}/FC.deps"
sleep 1
touch f90makefiles
fi],
[with_fortran=`test -z "${with_fortran_TRUE}" && echo yes || echo no`
DEPDIR="${DEPDIR}"])
AC_OUTPUT

echo ""
echo "Configuration done. Run $MAKE to build ScalES PPM."
echo ""
dnl remove a.out file and module files that some tests leave lying around
/bin/rm -f a.out *.${FCMODEXT}
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
