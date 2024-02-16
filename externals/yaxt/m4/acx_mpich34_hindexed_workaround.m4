dnl acx_mpich34_hindexed_workaround.m4 --- test whether MPICH can be
dnl run-time-patched to work around known problem with zero-stride
dnl hindexed datatype creation
dnl
dnl Copyright  (C)  2021  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec OpenMPI
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
dnl ACX_MPICH34_INDEXED_WORKAROUND([TEST-SOURCE-DIR=config/checksrc],
dnl                 [LIST-OF-WORK-AROUND-TESTS],
dnl                 [ACTION-IF-WORKAROUND-SUCCEEDS],
dnl                 [ACTION-IF-WORKAROUND-FAILED=AC_MSG_FAILURE],
dnl                 [PATCHED-FILES-DIR=config/workarounds/mpich34_yaksa_patch])
dnl
dnl Tries to detect MPICH installation from CC/mpirun.
dnl Assumes CC and CFLAGS are setup for MPICH.
dnl Then builds files needed for work-around in addition to test programs
dnl and checks if linking with those fixes the issue.
dnl See ACX_MPI_DEFECTS for further description of tests.
AC_DEFUN([ACX_MPICH34_DT_WORKAROUND],
  [AC_LANG_PUSH([C])
   dnl 1. detect Open MPI version
   is_affected_mpich=:
   acx_saved_CPPFLAGS=$CPPFLAGS
   acx_saved_LIBS=$LIBS
   while : ; do
   AC_COMPUTE_INT([acx_mpich_numversion],[MPICH_NUMVERSION],[AC_INCLUDES_DEFAULT
@%:@include <mpi.h>],
     [is_affected_mpich=false ; break])
   $is_affected_mpich || break
   AS_IF([expr "$acx_mpich_version >= 30400300 & $acx_mpich_version <= 40000000" >/dev/null],,[is_affected_mpich=false ; break])
   dnl 2.a. download header files
   acx_yaksa_add_workaround=m4_default([$5],[config/workarounds/mpich34_yaksa_patch])
   AS_IF([test -r "$acx_yaksa_add_workaround/xt_yaksa_indexed.c"],,
     [AC_CHECK_PROGS([DL_CMD], [wget curl],[/bin/false])
      AS_CASE([`$DL_CMD --version`],
        [*curl\ *], [
acx_fn_downloader()
{
  acx_temp=1
  for arg ; do
    if test "$acx_temp" -eq 1 ; then
      set -- ; acx_temp=0
    fi
    set -- "$][@" -O "$arg"
  done
  $DL_CMD -s -C - "$][@"
}
],
        [*GNU\ Wget*], [
acx_fn_downloader()
{
  $DL_CMD -q -c "$][@"
}
],
        [is_affected_mpich=false ; break])
      $MKDIR_P "$acx_yaksa_add_workaround"
      cd "$acx_yaksa_add_workaround"
      acx_temp='https://raw.githubusercontent.com/pmodels/yaksa/78739ed61c835a68991a34cff581b59b35be310e/'
      _AC_RUN_LOG([acx_fn_downloader \
     "$acx_temp"'src/frontend/types/yaksa_indexed.c' \
     "$acx_temp"'src/frontend/include/yaksi.h' \
     "$acx_temp"'src/frontend/include/yaksa.h.in' \
     "$acx_temp"'src/util/yaksu.h' \
     "$acx_temp"'src/util/yaksu_atomics.h' \
     "$acx_temp"'src/util/yaksu_base.h' \
     "$acx_temp"'src/util/yaksu_buffer_pool.h' \
     "$acx_temp"'src/util/yaksu_handle_pool.h' \
     "$acx_temp"'src/external/yuthash.h' \
     "$acx_temp"'src/backend/src/yaksur_pre.h' \
     "$acx_temp"'src/backend/src/yaksur_post.h' \
     "$acx_temp"'src/backend/seq/include/yaksuri_seq_pre.h' \
     "$acx_temp"'src/backend/seq/include/yaksuri_seq_post.h' \
     "$acx_temp"'src/backend/cuda/include/yaksuri_cuda_pre.h' \
     "$acx_temp"'src/backend/cuda/include/yaksuri_cuda_post.h' \
     "$acx_temp"'src/backend/ze/include/yaksuri_ze_pre.h' \
     "$acx_temp"'src/backend/ze/include/yaksuri_ze_post.h' \
         >&AS_MESSAGE_LOG_FD], [echo "downloading sources to patch"])
   dnl 2.b. and patch version into header template
   sed -e 's/@YAKSA_NUMVERSION@/unreleased/;s/@YAKSA_VERSION@/unreleased/' \
     yaksa.h.in >yaksa.h
      cd ../../..])
   CPPFLAGS="-I$acx_yaksa_add_workaround $CPPFLAGS"
   dnl 2.c. detect whether to use C11 atomics
   AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <pthread.h>

extern pthread_mutex_t yaksui_atomic_mutex;],
     [pthread_mutex_lock(&yaksui_atomic_mutex);])],
     [echo "" >"$acx_yaksa_add_workaround/yaksa_config.h"],
     [echo "@%:@define HAVE_C11_ATOMICS 1" >"$acx_yaksa_add_workaround/yaksa_config.h"])
   dnl 3. build fixed yaksa_indexed.c
   patch -p5 -d "$acx_yaksa_add_workaround" -o "$ac_pwd/conftest.c" \
     < "$srcdir/config/checkpatch/mpich_3.4_yaksa_hindexed.patch" \
     || { is_affected_mpich=false ; break ; }
   AC_COMPILE_IFELSE(,
     [acx_yaksa_add_workaround="$acx_yaksa_add_workaround/xt_yaksa_indexed.c"
      mv conftest.c "$acx_yaksa_add_workaround"
      ACX_MV_OBJ([conftest],[mpich_workaround])],
     [is_affected_mpich=false
      rm -f core conftest.err conftest.$ac_objext conftest.$ac_ext
      break])
   dnl 4. verify listed tests now succeed
   CPPFLAGS=$acx_saved_CPPFLAGS
   LIBS="mpich_workaround.$ac_objext $LIBS"
   acx_mpich_check_dir=`echo "acx_mpi_check_src_" | sed -e 's@/@<:@^/@:>@*$[]@@'`
   ACX_MPI_DEFECTS([$acx_mpich_check_dir],, [is_affected_mpich=false], [$2])
   rm -f "mpich_workaround."*
   break
   done
   CPPFLAGS=$acx_saved_CPPFLAGS
   LIBS=$acx_saved_LIBS
   dnl 5. execute action corresponding to outcome of tests
   AS_IF([$is_affected_mpich],[$3],
     [ASX_VAR_UNSET([acx_yaksa_add_workaround])
      m4_default([$4],
        [AC_MSG_FAILURE([Cannot apply MPICH 3.4-3.4.2 datatype bug work-around.])])])

   AC_LANG_POP([C])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
