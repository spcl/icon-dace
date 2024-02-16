dnl acx_ompi_dt_workaround.m4 --- test whether OpenMPI can be
dnl run-time-patched to work around known problem with zero-extent datatype
dnl
dnl Copyright  (C)  2019  Thomas Jahns <jahns@dkrz.de>
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
dnl ACX_MPI_DT_WORKAROUND([TEST-SOURCE-DIR=config/checksrc],
dnl                 [LIST-OF-WORK-AROUND-TESTS],
dnl                 [ACTION-IF-WORKAROUND-SUCCEEDS],
dnl                 [ACTION-IF-WORKAROUND-FAILED=AC_MSG_FAILURE],
dnl                 [PATCHED-FILES-DIR=config/workarounds/openmpi])
dnl
dnl Tries to detect Open MPI installation from CC/mpirun.
dnl Assumes CC and CFLAGS are setup for Open MPI.
dnl Then builds files needed for work-around in addition to test programs
dnl and checks if linking with those fixes the issue.
dnl See ACX_MPI_DEFECTS for further description of tests.
AC_DEFUN([ACX_OMPI_DT_WORKAROUND],
  [AC_LANG_PUSH([C])
   dnl 1. detect Open MPI version
   is_affected_ompi=:
   acx_saved_CPPFLAGS=$CPPFLAGS
   acx_saved_LIBS=$LIBS
   acx_saved_CFLAGS=$CFLAGS
   while : ; do
   AS_FOR([vvar_macro],[acx_temp],
     [OMPI_MAJOR_VERSION OMPI_MINOR_VERSION OMPI_RELEASE_VERSION],
     [AC_COMPUTE_INT([vvar_macro],[vvar_macro],[AC_INCLUDES_DEFAULT
@%:@include <mpi.h>],[is_affected_ompi=false ; break])])
   $is_affected_ompi || break
   acx_ompi_version="$OMPI_MAJOR_VERSION.$OMPI_MINOR_VERSION.$OMPI_RELEASE_VERSION"
   AS_CASE([$acx_ompi_version],
     [2.0.3|2.0.4|2.1.*|3.*|4.0.0|4.0.1], ,
     [1.*|2.0.0|2.0.1|2.0.2],       [is_affected_ompi=false ; break],
     [is_affected_ompi=false ; break])
   dnl todo: download sources and patch automatically as for mpich3.4
   dnl 2. deduce header location and flags
   ASX_VAR_UNSET([acx_ompi_path])
   AC_ARG_VAR([OMPI_MCA_C_INCLUDE],[flags needed to compile Open MPI MCA modules])
   AS_IF([test x"${OMPI_MCA_C_INCLUDE+set}" != xset],
     [dnl 2.a. find header directory candidates if not set by user
      AC_CACHE_CHECK([Open MPI MCA headers],[acx_cv_ompi_mca_c_include],
        [# Extract the first word of "$CC", so it can be a program name with args.
         set dummy $CC; ac_word=$[2]
         AS_CASE([$ac_word],
           [/*mpicc],
           [AS_IF([expr "`$ac_word --showme:version`" : '.* Open MPI '"$acx_ompi_version" >/dev/null],
              [acx_ompi_path_candidates=`echo "$ac_word" | sed -e 's@/@<:@^/@:>@*mpicc$][@@'`])],
           [*mpicc],
           [acx_ompi_path_candidates=`which "$ac_word" 2>/dev/null | sed -e 's@/@<:@^/@:>@*mpicc$][@@'`],
           [acx_ompi_path_candidates=`echo " $CPPFLAGS $CFLAGS " \
 | sed -e 's/ -I \(@<:@^ 	@:>@*\)/ -I\1/g
s/ \(@<:@^- 	@:>@\|-@<:@^I 	@:>@\)@<:@^ 	@:>@*//g
s/ -I/ /g
s/^ *//;s/ *$//'`])
         AS_FOR([incdir],[acx_temp],[$acx_ompi_path_candidates],
           [AS_IF([test -r "incdir/mpi.h" -a -d "incdir/openmpi/opal"],
              [acx_cv_ompi_mca_c_include="-I]incdir[/openmpi" ; break],
              [test -r "incdir/../include/mpi.h" -a -d "incdir/../include/openmpi/opal"],
              [acx_temp=`echo "incdir" | sed -e 's@/@<:@^/@:>@*$][@@'`/include
               AS_IF([test -r "$acx_temp/mpi.h" -a -d "$acx_temp/openmpi/opal"],
                 [acx_cv_ompi_mca_c_include="-I$acx_temp/openmpi"; break])])])])
      AS_IF([test x"${acx_cv_ompi_mca_c_include+set}" = xset],
        [OMPI_MCA_C_INCLUDE=$acx_cv_ompi_mca_c_include],
        [is_affected_ompi=false ; break])])
   dnl 2.b. establish Open MPI build CFLAGS if possible
   AC_ARG_VAR([OMPI_BUILD_CFLAGS],[C compiler flags matching used Open MPI])
   AS_IF([test x"${OMPI_BUILD_CFLAGS+set}" != xset],
     [AC_CACHE_CHECK([Open MPI C build flags],[acx_cv_ompi_build_cflags],
        [# Extract the first word of "$CC", so it can be a program name with args.
         set dummy $CC; ac_word=$[2]
         AS_CASE([$ac_word],
           [/*mpicc],
           [AS_IF([expr "`$ac_word --showme:version`" : '.* Open MPI '"$acx_ompi_version" >/dev/null],
              [acx_ompi_path_candidates=`echo "$ac_word" | sed -e 's@/@<:@^/@:>@*mpicc$][@@'`])],
           [*mpicc],
           [acx_ompi_path_candidates=`which "$ac_word" 2>/dev/null | sed -e 's@/@<:@^/@:>@*mpicc$][@@'`],
           [acx_ompi_path_candidates=`echo " $CPPFLAGS $CFLAGS " \
 | sed -e 's/ -I \(@<:@^ 	@:>@*\)/ -I\1/g
s/ \(@<:@^- 	@:>@\|-@<:@^I 	@:>@\)@<:@^ 	@:>@*//g
s/ -I/ /g
s/^ *//;s/ *$//'`])
         # from the path candidates derived above, try to find ompi_info
         AS_FOR([bindir],[acx_temp],[$acx_ompi_path_candidates],
           [AS_IF([test -x "bindir/ompi_info"],
              [acx_cv_ompi_build_cflags="bindir/ompi_info" ; break],
              [test -x "bindir/../bin/ompi_info"],
              [acx_temp=`echo "bindir" | sed -e 's@/@<:@^/@:>@*$][@@'`
               AS_IF([test -x "$acx_temp/bin/ompi_info"],
                 [acx_cv_ompi_build_cflags="$acx_temp/bin/ompi_info"; break])])])
         AS_IF([test x"${acx_cv_ompi_build_cflags+set}" = xset && expr "`$acx_cv_ompi_build_cflags --version | head -n 1`" : 'Open MPI v'"$acx_ompi_version" >/dev/null],
          [acx_cv_ompi_build_cflags=`$acx_cv_ompi_build_cflags -c | sed -n -e '/^ *Build CFLAGS:/s/^ *Build CFLAGS: *//p'`])])
      AS_IF([test x"${acx_cv_ompi_build_cflags+set}" = xset],
        [OMPI_BUILD_CFLAGS=$acx_cv_ompi_build_cflags],
        [is_affected_ompi=false ; break])],
     [AC_MSG_CHECKING([Open MPI C build flags])
      AC_MSG_RESULT([$OMPI_BUILD_CFLAGS])])
   dnl 2.c. make sure header matches what compiler uses
   CPPFLAGS="$OMPI_MCA_C_INCLUDE $CPPFLAGS"
   CFLAGS=$OMPI_BUILD_CFLAGS
   AS_FOR([vvar_macro],[acx_temp],
     [OMPI_MAJOR_VERSION OMPI_MINOR_VERSION OMPI_RELEASE_VERSION],
     [AC_COMPUTE_INT([OPAL_]vvar_macro,[vvar_macro],[AC_INCLUDES_DEFAULT
@%:@include <opal_config.h>],[is_affected_ompi=false ; break])
      AS_IF([eval test x\"\$OPAL_]vvar_macro[\" != x\"\$vvar_macro\"],
        [is_affected_ompi=false ; break])])
   $is_affected_ompi || break
   dnl 3. build fixed opal/datatype/opal_datatype_add.c,
   dnl corresponding to Open MPI version
   AS_CASE([$OMPI_MAJOR_VERSION],
     [2],
     [acx_opal_datatype_add_workaround=opal_datatype_add.c],
     [3|4],
     [acx_opal_datatype_add_workaround=opal_datatype_optimize.c])
   acx_temp="m4_default([$5],[config/workarounds/openmpi])/v$acx_ompi_version"
   AS_IF([test -r "$ac_pwd/$acx_temp/$acx_opal_datatype_add_workaround"],
     [acx_opal_datatype_add_workaround="$ac_pwd/$acx_temp/$acx_opal_datatype_add_workaround"
      cp "$acx_opal_datatype_add_workaround" conftest.c \
        || { is_affected_ompi=false ; break ; }],
     [test -r "$srcdir/$acx_temp/$acx_opal_datatype_add_workaround"],
     [acx_opal_datatype_add_workaround="$srcdir/$acx_temp/$acx_opal_datatype_add_workaround"
      cp "$acx_opal_datatype_add_workaround" conftest.c \
        || { is_affected_ompi=false ; break ; }],
     [dnl patched source not availabe, download and patch
      AC_CHECK_PROGS([DL_CMD], [wget curl],[/bin/false])
      AS_CASE([`$DL_CMD --version`],
        [*curl\ *], [
acx_fn_downloader()
{
  acx_tmp=1
  for arg ; do
    if test "$acx_tmp" -eq 1 ; then
      set -- ; acx_tmp=0
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
        [is_affected_ompi=false ; break])
      _AC_RUN_LOG([$MKDIR_P "$ac_pwd/$acx_temp/stage" && \
        cd "$ac_pwd/$acx_temp/stage" && \
        acx_fn_downloader \
          "https://raw.githubusercontent.com/open-mpi/ompi/v$acx_ompi_version/opal/datatype/$acx_opal_datatype_add_workaround" && \
        cd "$ac_pwd" && \
        patch -p3 -d "$acx_temp/stage" -o "$ac_pwd/conftest.c" \
          <"$srcdir/config/checkpatch/openmpi_datatype.$acx_opal_datatype_add_workaround.patch"],
        [echo "downloading and patching Open MPI source file $acx_opal_datatype_add_workaround"]) \
        || { cd "$ac_pwd" ; is_affected_ompi=false ; break ; }
      acx_opal_datatype_add_workaround="$ac_pwd/$acx_temp/stage/$acx_opal_datatype_add_workaround"])
   as_dir=src as_fn_mkdir_p
   AC_COMPILE_IFELSE(,
     [ACX_MV_OBJ([conftest],[ompi_workaround])
      AS_CASE([$acx_opal_datatype_add_workaround],
        ["$ac_pwd/$acx_temp/stage/"*],
        [acx_temp=`echo "$acx_opal_datatype_add_workaround" | sed -e 's/\/stage\(\/@<:@^/@:>@*\)$/\\1/'`
         mv conftest.c "$acx_temp"
         acx_opal_datatype_add_workaround=$acx_temp])],
     [is_affected_ompi=false
      rm -f core conftest.err conftest.$ac_objext conftest.$ac_ext
      AS_CASE([$acx_opal_datatype_add_workaround],
        ["$acx_temp/stage"*],
        [rm -f "$acx_opal_datatype_add_workaround"])
      break])
   dnl 4. verify listed tests now succeed
   CPPFLAGS=$acx_saved_CPPFLAGS
   CFLAGS=$acx_saved_CFLAGS
   dnl prepend -lopen-pal to LIBS if opal symbol cannot already be resolved
   AC_SEARCH_LIBS([opal_output],[open-pal],,[is_affected_ompi=false ; break])
   LIBS="ompi_workaround.$ac_objext $LIBS"
   acx_ompi_check_dir=`echo "acx_mpi_check_src_" | sed -e 's@/@<:@^/@:>@*$[]@@'`
   ACX_MPI_DEFECTS([$acx_ompi_check_dir],, [is_affected_ompi=false], [$2])
   rm -f "ompi_workaround."*
   break
   done
   CFLAGS=$acx_saved_CFLAGS
   CPPFLAGS=$acx_saved_CPPFLAGS
   LIBS=$acx_saved_LIBS
   dnl 5. execute action corresponding to outcome of tests
   AS_IF([$is_affected_ompi],[$3],
     [ASX_VAR_UNSET([acx_opal_datatype_add_workaround])
      m4_default([$4],
        [AC_MSG_FAILURE([Cannot apply Open MPI datatype bug work-around.])])])

   AC_LANG_POP([C])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
