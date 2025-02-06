# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# AC_LANG([HIP])
# -----------------------------------------------------------------------------
# Define language HIP with _AC_LANG_ABBREV set to "hip", _AC_LANG_PREFIX set to
# "HIP" (i.e. _AC_LANG_PREFIX[]FLAGS equals to 'HIPFLAGS'), _AC_CC set to
# "HIPCXX", and various language-specific Autoconf macros (e.g. AC_LANG_SOURCE)
# copied from C++ language.
#
AC_LANG_DEFINE([HIP], [hip], [HIP], [HIPCXX], [C++],
[ac_ext=cc
ac_compile='$HIPCXX -c $HIPFLAGS conftest.$ac_ext >&AS_MESSAGE_LOG_FD'
ac_link='$HIPCXX -o conftest$ac_exeext $HIPFLAGS $LDFLAGS conftest.$ac_ext $LIBS >&AS_MESSAGE_LOG_FD'
ac_compiler_gnu=$acx_cv_hip_compiler_gnu
])
dnl Avoid mixing macro definitions that are relevant for C or C++ compiler and
dnl might be irrelevant for HIP:
m4_copy_force([AC_LANG_CONFTEST()], [AC_LANG_CONFTEST(HIP)])
AC_DEFUN([AC_LANG_COMPILER(HIP)], [AC_REQUIRE([ACX_PROG_HIPCXX])])

# ACX_LANG_HIP()
# ACX_LANG_PUSH_HIP()
# ACX_LANG_POP_HIP()
# -----------------------------------------------------------------------------
# A set of macros wrapping AC_LANG([HIP]), AC_LANG_PUSH([HIP]), and
# AC_LANG_POP([HIP]), respectively, that can be used in the "configure.ac" to
# switch to/from HIP language without the need to call "aclocal/autoreconf"
# with additional flag "-I /path/to/dir/with/this/file". Useful for projects
# that do not use Automake for the makefile generation and, therefore, cannot
# specify "ACLOCAL_AMFLAGS" in the top-level "Makefile.am".
#
AC_DEFUN([ACX_LANG_HIP], [AC_LANG([HIP])])
AC_DEFUN([ACX_LANG_PUSH_HIP], [AC_LANG_PUSH([HIP])])
AC_DEFUN([ACX_LANG_POP_HIP], [AC_LANG_POP([HIP])])

# ACX_PROG_HIPCXX([LIST-OF-COMPILERS = hipcc])
# -----------------------------------------------------------------------------
# Searches for the HIP C++ compiler command among the values of the
# blank-separated list LIST-OF-COMPILERS (defaults to a single value "hipcc").
# Declares precious variables HIPCXX and HIPFLAGS to be set to the compiler
# command and the compiler flags, respectively. If the environment variable
# HIPCXX is set, values of the list LIST-OF-COMPILERS are ignored.
#
# Checks whether the compiler supports an ISO C++ standard and/or whether it
# can compile a basic HIP program. The result of the check is stored in the
# shell variable acx_prog_hipcxx_works, which can be either "cxx14", (ISO C++
# 2014, ISO C++ 2011 and the basic HIP features are supported), "cxx11" (ISO
# C++ 2011 and the basic HIP features are supported), "basic" (only the basic
# HIP features are supported) and "no" (the compiler cannot compile even the
# basic HIP program). If the compiler supports the C++ standard, the flag that
# is required to enable it is cached in the
# acx_cv_prog_hipcxx_${acx_prog_hipcxx_works}_flag variable. The users are
# expected to either append the flag to HIPCXX or prepended it to HIPFLAGS
# based on the value of the acx_prog_hipcxx_works shell variable.
#
AC_DEFUN_ONCE([ACX_PROG_HIPCXX],
  [AC_LANG_PUSH([HIP])dnl
   AC_ARG_VAR([HIPCXX], [HIP C++ compiler command])dnl
   AC_ARG_VAR([HIPFLAGS], [HIP C++ compiler flags])dnl
   _AC_ARG_VAR_LDFLAGS()dnl
   _AC_ARG_VAR_LIBS()dnl
   AS_IF([test -z "$HIPCXX"],
     [AC_CHECK_PROGS([HIPCXX], [m4_default([$1], [hipcc])])])
   _AS_ECHO_LOG([checking for HIP C++ compiler version])
   set dummy $ac_compile
   ac_compiler=$[2]
   _AC_DO_LIMIT([$ac_compiler --version >&AS_MESSAGE_LOG_FD])
   m4_expand_once([_AC_COMPILER_EXEEXT])[]dnl
   m4_expand_once([_AC_COMPILER_OBJEXT])[]dnl
   _AC_LANG_COMPILER_GNU
   acx_test_HIPFLAGS=${HIPFLAGS+set}
   acx_save_HIPFLAGS=$HIPFLAGS
   AC_CACHE_CHECK([whether $HIPCXX accepts -g], [acx_cv_prog_hipc_g],
     [acx_save_hip_werror_flag=$ac_hip_werror_flag
      ac_hip_werror_flag=yes
      acx_cv_prog_hipc_g=no
      HIPFLAGS='-g'
      _AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
        [acx_cv_prog_hipc_g=yes],
        [HIPFLAGS=''
         _AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
           [],
           [ac_hip_werror_flag=$acx_save_hip_werror_flag
            HIPFLAGS='-g'
            _AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
              [acx_cv_prog_hipc_g=yes])])])
      ac_hip_werror_flag=$acx_save_hip_werror_flag])
   AS_IF(
     [test "$acx_test_HIPFLAGS" = set],
     [HIPFLAGS=$acx_save_HIPFLAGS],
     [test "$acx_cv_prog_hipc_g" = yes],
     [AS_VAR_IF([ac_cv_hip_compiler_gnu], [yes],
        [HIPFLAGS='-g -O2'], [HIPFLAGS='-g'])],
     [AS_VAR_IF([ac_cv_hip_compiler_gnu], [yes],
        [HIPFLAGS='-O2'], [HIPFLAGS=''])])
   acx_prog_hipcxx_works=no
   m4_map_sep([_ACX_HIP_CXX_TEST], [
], [[14],[11]])
   _ACX_HIP_BASIC_TEST
   AC_LANG_POP([HIP])])

# _ACX_HIP_CXX_TEST(STANDARD)
# -----------------------------------------------------------------------------
# Checks whether the HIP C++ compiler supports the ISO C++ STANDARD (only
# checks for ISO C++ 2014 and ISO C++ 2011 are currently supported) and can
# compile a basic HIP program. The value of STANDARD is expected to be the last
# two digits of the standard's year (e.g. 11). The check is skipped if the
# acx_prog_hipcxx_works shell variable is set to a value other than "no" (in a
# normal scenario that means that a more recent ISO C++ standard is supported).
#
# If successful, sets the shell variable acx_prog_hipcxx_works to
# "cxx[]STANDARD" (e.g. "cxx11") and the cache variable
# acx_cv_prog_hipcxx_cxx[]STANDARD[]_flag to the flag that is required to
# enable the standard.
#
AC_DEFUN([_ACX_HIP_CXX_TEST],
  [AC_LANG_ASSERT([HIP])dnl
   AC_REQUIRE([_ACX_HIP_CXX$1_TEST_PROGRAM])dnl
   AS_VAR_IF([acx_prog_hipcxx_works], [no],
     [AC_MSG_CHECKING([for $HIPCXX option to enable C++$1 features])
      m4_pushdef([acx_cache_var], [acx_cv_prog_hipcxx_cxx$1_flag])dnl
      AC_CACHE_VAL([acx_cache_var],
        [acx_cache_var=unsupported
         AC_LANG_CONFTEST([$acx_hip_conftest_cxx$1_program])
         acx_save_HIPFLAGS=$HIPFLAGS
         for acx_flag in '' m4_normalize(m4_defn([_ACX_HIP_CXX$1_OPTIONS]))
         do
           HIPFLAGS="$acx_flag $acx_save_HIPFLAGS"
           _AC_COMPILE_IFELSE([], [acx_cache_var=$acx_flag])
           test "x$acx_cache_var" != xunsupported && break
         done
         HIPFLAGS=$acx_save_HIPFLAGS
         rm -f conftest.$ac_ext])
      AS_IF([test -n "$acx_cache_var"],
        [AC_MSG_RESULT([$acx_cache_var])],
        [AC_MSG_RESULT([none needed])])
      AS_IF([test "x$acx_cache_var" != xunsupported],
        [acx_prog_hipcxx_works=cxx$1])
      m4_popdef([acx_cache_var])])])

# _ACX_HIP_BASIC_TEST()
# -----------------------------------------------------------------------------
# Checks whether the HIP C++ compiler can compile a basic HIP program. The
# check is skipped if the acx_prog_hipcxx_works shell variable is set to a
# value other than "no" (in a normal scenario that means that the compiler
# supports an ISO C++ standard, as well as the basic HIP features).
#
# If successful, sets the shell variable acx_prog_hipcxx_works to "basic".
#
# The result is cached in the acx_cv_prog_hipxx_basic variable.
#
AC_DEFUN([_ACX_HIP_BASIC_TEST],
  [AC_LANG_ASSERT([HIP])dnl
   AC_REQUIRE([_ACX_HIP_BASIC_TEST_PROGRAM])dnl
   AS_VAR_IF([acx_prog_hipcxx_works], [no],
     [AC_CACHE_CHECK([whether $HIPCXX can compile basic HIP code],
        [acx_cv_prog_hipxx_basic],
        [_AC_COMPILE_IFELSE([$acx_hip_conftest_basic_program],
           [acx_cv_prog_hipxx_basic=yes],
           [acx_cv_prog_hipxx_basic=no])])
      AS_VAR_IF([acx_cv_prog_hipxx_basic], [yes],
        [acx_prog_hipcxx_works=basic])])])

# _ACX_HIP_BASIC_TEST_GLOBALS()
# _ACX_HIP_BASIC_TEST_MAIN()
# _ACX_HIP_BASIC_TEST_PROGRAM()
# _ACX_HIP_CXX11_TEST_GLOBALS()
# _ACX_HIP_CXX11_TEST_MAIN()
# _ACX_HIP_CXX11_TEST_PROGRAM()
# _ACX_HIP_CXX14_TEST_GLOBALS()
# _ACX_HIP_CXX14_TEST_MAIN()
# _ACX_HIP_CXX14_TEST_PROGRAM()
# -----------------------------------------------------------------------------
# A set of macros that expand to the shell code in the INIT_PREPARE section of
# the configure script that assigns respective shell variables to the source
# code of the test programs that are used when checking for the HIP C++
# compiler features.
#
AC_DEFUN([_ACX_HIP_BASIC_TEST_GLOBALS],
  [m4_divert_once([INIT_PREPARE],
[[# Basic HIP test code (global declarations):
acx_hip_conftest_basic_globals='#include <hip/hip_runtime.h>
__global__ void conftest_foo() {}'
]])])

AC_DEFUN([_ACX_HIP_BASIC_TEST_MAIN],
  [m4_divert_once([INIT_PREPARE],
[[# Basic HIP test code (body of main):
acx_hip_conftest_basic_main='conftest_foo<<<1, 1>>>();'
]])])

AC_DEFUN([_ACX_HIP_BASIC_TEST_PROGRAM],
  [AC_REQUIRE([_ACX_HIP_BASIC_TEST_GLOBALS])dnl
   AC_REQUIRE([_ACX_HIP_BASIC_TEST_MAIN])dnl
   m4_divert_once([INIT_PREPARE],
[[# Basic HIP test code (complete):
acx_hip_conftest_basic_program="$acx_hip_conftest_basic_globals
int main (int argc, char **argv)
{
  $acx_hip_conftest_basic_main
  return 0;
}
"
]])])

AC_DEFUN([_ACX_HIP_CXX11_TEST_GLOBALS],
  [m4_divert_once([INIT_PREPARE],
[[# C++11 HIP test code (global declarations):
acx_hip_conftest_cxx11_globals='
#if !defined __cplusplus || __cplusplus < 201103L
# error "Compiler does not advertise C++11 conformance"
#endif

namespace cxx11test {
  template<typename T>
  class ConftestClassCXX11 {
  public:
    void operator() (T elem) {}
  };
}'
]])])

AC_DEFUN([_ACX_HIP_CXX11_TEST_MAIN],
  [m4_divert_once([INIT_PREPARE],
[[# C++11 HIP test code (body of main):
acx_hip_conftest_cxx11_main='
int* pcxx11 = nullptr;
cxx11test::ConftestClassCXX11<int> some_class =
  cxx11test::ConftestClassCXX11<int>();
int icxx11 = 1;
pcxx11 = &icxx11;
some_class(*pcxx11);'
]])])

AC_DEFUN([_ACX_HIP_CXX11_TEST_PROGRAM],
  [AC_REQUIRE([_ACX_HIP_BASIC_TEST_GLOBALS])dnl
   AC_REQUIRE([_ACX_HIP_BASIC_TEST_MAIN])dnl
   AC_REQUIRE([_ACX_HIP_CXX11_TEST_GLOBALS])dnl
   AC_REQUIRE([_ACX_HIP_CXX11_TEST_MAIN])dnl
   m4_divert_once([INIT_PREPARE],
[[# C++11 HIP test code (complete):
acx_hip_conftest_cxx11_program="$acx_hip_conftest_basic_globals
$acx_hip_conftest_cxx11_globals
int main (int argc, char **argv)
{
  $acx_hip_conftest_basic_main
  $acx_hip_conftest_cxx11_main
  return 0;
}
"
]])])

AC_DEFUN([_ACX_HIP_CXX14_TEST_GLOBALS],
  [m4_divert_once([INIT_PREPARE],
[[# C++14 HIP test code (global declarations):
acx_hip_conftest_cxx14_globals='
#if !defined __cplusplus || __cplusplus < 201402L
# error "Compiler does not advertise C++14 conformance"
#endif
'
]])])

AC_DEFUN([_ACX_HIP_CXX14_TEST_MAIN],
  [m4_divert_once([INIT_PREPARE],
[[# C++14 HIP test code (body of main):
acx_hip_conftest_cxx14_main=''
]])])

AC_DEFUN([_ACX_HIP_CXX14_TEST_PROGRAM],
  [AC_REQUIRE([_ACX_HIP_BASIC_TEST_GLOBALS])dnl
   AC_REQUIRE([_ACX_HIP_BASIC_TEST_MAIN])dnl
   AC_REQUIRE([_ACX_HIP_CXX11_TEST_GLOBALS])dnl
   AC_REQUIRE([_ACX_HIP_CXX11_TEST_MAIN])dnl
   AC_REQUIRE([_ACX_HIP_CXX14_TEST_GLOBALS])dnl
   AC_REQUIRE([_ACX_HIP_CXX14_TEST_MAIN])dnl
   m4_divert_once([INIT_PREPARE],
[[# C++14 HIP test code (complete):
acx_hip_conftest_cxx14_program="$acx_hip_conftest_basic_globals
$acx_hip_conftest_cxx11_globals
$acx_hip_conftest_cxx14_globals
int main (int argc, char **argv)
{
  $acx_hip_conftest_basic_main
  $acx_hip_conftest_cxx11_main
  $acx_hip_conftest_cxx14_main
  return 0;
}
"
]])])

# _ACX_HIP_CXX11_OPTIONS()
# -----------------------------------------------------------------------------
# Expands into a space-separated list of known flags needed to enable the
# support for the ISO C++ 2011 standard when running the HIP C++ compiler.
#
m4_define([_ACX_HIP_CXX11_OPTIONS], [-std=gnu++11 -std=c++11])

# _ACX_HIP_CXX14_OPTIONS()
# -----------------------------------------------------------------------------
# Expands into a space-separated list of known flags needed to enable the
# support for the ISO C++ 2014 standard when running the HIP C++ compiler.
#
m4_define([_ACX_HIP_CXX14_OPTIONS], [-std=gnu++14 -std=c++14])
