# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------
AC_PREREQ([2.69])

ACX_POLYMORPHIC_OPTION([enable], [fcgroup-*], [icon_fcgroups])

AC_INIT([icon], [2024.01], [], [icon-2024.01])

AC_CONFIG_MACRO_DIR([m4])
m4_pattern_forbid([^A[SC]?X_])dnl
AC_CONFIG_AUX_DIR([utils])
AC_CONFIG_SRCDIR([src/drivers/icon.f90])

icon_usr_note="
********************************************************************************
The users are recommended to run site-specific wrappers instead of this generic
configure script.

You can find the wrappers in '${srcdir}/config/<organization>/'.

If you are running a wrapper and still see this message than either you are
using a wrong one or the wrapper is broken, or there is really a problem with
the configure script.
********************************************************************************
"
m4_define([AC_MSG_FAILURE],
[{ AS_MESSAGE([error: in `$ac_pwd':], 2)
AC_MSG_ERROR([$1
See `config.log' for more details
${icon_usr_note}], [$2]); }])
dnl End of the temporary message injection. Delete up until this line once the
dnl message is not needed anymore.

dnl We need to make sure that the source directory is clean when building
dnl out-of-source for two reasons:
dnl   1. Building systems of most of the bundled libraries are Automake-based
dnl      and, therefore, they rely on the standard non-selective search path
dnl      (i.e., VPATH) feature of make for out-of-source builds. This makes
dnl      having a clean source directory a hard requirement for them. Otherwise,
dnl      instead of generating object files in the build directory, make would
dnl      simply take them from the dirty source directory. For that reason, the
dnl      respective configure scripts will run a check similar to the one below
dnl      and fail if their source directories are dirty. It is better that we
dnl      run the check much earlier and fail the configuration with a more
dnl      informative error message.
dnl   2. The makefiles of ICON rely on the selective search path (i.e. vpath)
dnl      feature of GNU Make. Therefore, the problem with the object files that
dnl      is described above is irrelevant. However, icon.mk supports Fortran
dnl      source file preprocessing. The input and output files of the most of
dnl      the preprocessing steps are valid Fortran files. Thus, both the
dnl      original files in the source directory and the preprocessed files in
dnl      the build directory have the same suffix '.f90', which makes it
dnl      impossible to declare selective search paths that would work for the
dnl      original files but wouldn't work for the preprocessed files. This might
dnl      lead to the preprocessed files being taken from the dirty source
dnl      directory instead of being generated in the build directory. A
dnl      workaround for this, which had been implemented before, is to refer to
dnl      the preprocessed files by their absolute paths. However, it makes the
dnl      makefiles overcomplicated and error-prone. Instead, we simply forbid
dnl      the source directory to be dirty in the case of out-of-source build.
AS_IF([test "`cd $srcdir && pwd`" != "`pwd`" && dnl
test -f "$srcdir/config.status"],
  [AC_MSG_ERROR([source directory already configured; dnl
run "make distclean" there first])])

############################# OPTION DECLARATIONS ##############################

AC_ARG_ENABLE([rpaths],
  [AS_HELP_STRING([--enable-rpaths],
     [add directories specified with -L flags in LDFLAGS and LIBS to the
runtime library search paths (RPATH) @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rpaths=yes])],
  [enable_rpaths=yes])

AC_ARG_ENABLE([silent-rules],
  [AS_HELP_STRING([--enable-silent-rules],
     [less verbose build output (undo: "make V=1") @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_silent_rules=yes])],
  [enable_silent_rules=yes])

AC_ARG_ENABLE([delayed-config],
  [AS_HELP_STRING([--enable-delayed-config],
     [configure the bundled libraries at the build time @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_delayed_config=yes])],
  [enable_delayed_config=yes])

m4_divert_once([HELP_ENABLE], [[
Model Features:]])

AC_ARG_ENABLE([atmo],
  [AS_HELP_STRING([--enable-atmo],
     [enable the atmosphere component of ICON @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_atmo=yes])],
  [enable_atmo=yes])

AC_ARG_ENABLE([les],
  [AS_HELP_STRING([--enable-les],
     [enable the Large-Eddy Simulation component of ICON @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_les=yes])],
  [enable_les=yes])

AC_ARG_ENABLE([upatmo],
  [AS_HELP_STRING([--enable-upatmo],
     [enable the upper atmosphere component of ICON @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_upatmo=yes])],
  [enable_upatmo=yes])

AC_ARG_ENABLE([ocean],
  [AS_HELP_STRING([--enable-ocean],
     [enable the ocean component of ICON @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_ocean=yes])],
  [enable_ocean=yes])

AC_ARG_ENABLE([jsbach],
  [AS_HELP_STRING([--enable-jsbach],
     [enable the land component JSBACH @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_jsbach=yes])],
  [enable_jsbach=yes])

AC_ARG_ENABLE([jsbach-hd],
  [AS_HELP_STRING([--enable-jsbach-hd],
     [enable the hydrological discharge scheme of the JSBACH land component
(JSBACH HD) @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_jsbach_hd=yes])],
  [enable_jsbach_hd=auto])

AS_VAR_IF([enable_jsbach_hd], [auto],
  [AS_VAR_COPY([enable_jsbach_hd], [enable_jsbach])])

AC_ARG_ENABLE([waves],
  [AS_HELP_STRING([--enable-waves],
     [enable the ocean surface wave component of ICON @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_waves=yes])],
    [enable_waves=no])

AC_ARG_ENABLE([coupling],
  [AS_HELP_STRING([--enable-coupling],
     [enable the coupling @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_coupling=yes])],
  [enable_coupling=yes])

AC_ARG_ENABLE([aes],
  [AS_HELP_STRING([--enable-aes],
     [enable the AES physics package @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_aes=yes])],
  [enable_aes=yes])

AC_ARG_ENABLE([nwp],
  [AS_HELP_STRING([--enable-nwp],
     [enable the NWP physics package @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_nwp=yes])],
  [enable_nwp=yes])

AC_ARG_ENABLE([ecrad],
  [AS_HELP_STRING([--enable-ecrad],
     [enable usage of the ECMWF radiation scheme (ECRAD) @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_ecrad=yes])],
  [enable_ecrad=no])

AC_ARG_ENABLE([rte-rrtmgp],
  [AS_HELP_STRING([--enable-rte-rrtmgp],
     [enable usage of the RTE+RRTMGP toolbox for radiation calculations
@<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rte_rrtmgp=yes])],
  [enable_rte_rrtmgp=yes])

AC_ARG_ENABLE([rttov],
  [AS_HELP_STRING([--enable-rttov],
     [enable usage of the radiative transfer model for TOVS
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_rttov=yes])],
  [enable_rttov=no])

AC_ARG_ENABLE([dace],
  [AS_HELP_STRING([--enable-dace],
     [enable the DACE modules for data assimilation @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_dace=yes])],
  [enable_dace=no])

AC_ARG_ENABLE([emvorado],
  [AS_HELP_STRING([--enable-emvorado],
     [enable the radar forward operator EMVORADO @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_emvorado=yes])],
  [enable_emvorado=no])

AC_ARG_ENABLE([art],
  [AS_HELP_STRING([--enable-art],
     [enable the aerosols and reactive trace component ART
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_art=yes])],
  [enable_art=no])

AC_ARG_ENABLE([art-gpl],
  [AS_HELP_STRING([--enable-art-gpl],
     [enable GPL-licensed code parts of the ART component @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_art_gpl=yes])],
  [enable_art_gpl=no])

AC_ARG_ENABLE([comin],
  [AS_HELP_STRING([--enable-comin],
     [enable the ICON community interfaces @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_comin=yes])],
  [enable_comin=no])

AC_ARG_ENABLE([hd],
[  --enable-hd             enable the Hydrological Discharge (HD) model. The
                          value of the argument must be one of the following
                          @<:@default=no@:>@:
                            05deg|yes  enable for the 0.5-degree global domain
                            5min       enable for the 5min global domain
                            no         disable the HD model],
[AS_CASE([$enableval],
   [yes], [enable_hd=05deg],
   [05deg|5min|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-hd='$enableval'; valid values are '05deg', 'yes' (same dnl
as '05deg'), '5min', 'no'])])],
[enable_hd=no])

AC_ARG_ENABLE([acm-license],
  [AS_HELP_STRING([--enable-acm-license],
     [enable code parts that require accepting the ACM Software License ]dnl
[Agreement @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_acm_license=yes])],
  [enable_acm_license=no])

m4_divert_once([HELP_ENABLE], [[
Infrastructural Features:]])

AC_ARG_ENABLE([mpi],
  [AS_HELP_STRING([--enable-mpi],
     [enable MPI (parallelization) support @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mpi=yes])],
  [enable_mpi=yes])

AC_ARG_ENABLE([mpi-checks],
  [AS_HELP_STRING([--enable-mpi-checks],
     [enable configure-time checks of MPI library for known defects ]dnl
[@<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mpi_checks=yes])],
  [enable_mpi_checks=yes])

AC_ARG_ENABLE([active-target-sync],
  [AS_HELP_STRING([--enable-active-target-sync],
     [enable MPI active target mode (otherwise, passive target mode is used)
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_active_target_sync=yes])],
  [enable_active_target_sync=no])

dnl switch to disable MPI_Rget for systems where it doesn't work
AC_ARG_ENABLE([mpi-rget],
  [AS_HELP_STRING([--enable-mpi-rget],
     [enable usage of the MPI_Rget routine @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_mpi_rget=yes])],
  [enable_mpi_rget=auto])

AC_ARG_ENABLE([mpi-gpu],
  [AS_HELP_STRING([--enable-mpi-gpu],
     [enable usage of the GPU-aware MPI features @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mpi_gpu=yes])],
  [enable_mpi_gpu=no])

AC_ARG_ENABLE([async-io-rma],
  [AS_HELP_STRING([--enable-async-io-rma],
     [enable MPI RMA for asynchronous I/O @<:@default=yes@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_async_io_rma=yes])],
  [enable_async_io_rma=yes])

AC_ARG_ENABLE([openmp],
  [AS_HELP_STRING([--enable-openmp],
     [enable OpenMP support @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_openmp=yes])],
  [enable_openmp=no])

AC_ARG_ENABLE([gpu],
[  --enable-gpu            enable GPU support. The value of the argument must
                          be one of the following @<:@default=no@:>@:
                            openacc+cuda  enable support with OpenACC and CUDA
                            openacc+hip   enable support with OpenACC and HIP
                            openacc       alias for 'openacc+cuda'
                            yes           alias for 'openacc+cuda'
                            no            disable GPU support],
[AS_CASE([$enableval],
   [openacc|yes], [enable_gpu='openacc+cuda'],
   [openacc+cuda|openacc+hip|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-gpu='$enableval'; valid values are 'openacc+cuda', 'openacc+hip', dnl
'openacc' (same as 'openacc+cuda'), 'yes' (same as 'openacc+cuda'), 'no'])])],
[enable_gpu=no])

AS_CASE([$enable_gpu],
  [openacc+*], [icon_gpu_primary=openacc],
  [icon_gpu_primary=no])
AS_CASE([$enable_gpu],
  [*+cuda], [icon_gpu_secondary=cuda],
  [*+hip], [icon_gpu_secondary=hip],
  [icon_gpu_secondary=no])

AC_ARG_ENABLE([grib2],
  [AS_HELP_STRING([--enable-grib2],
     [enable GRIB2 I/O @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_grib2=yes])],
  [enable_grib2=no])

AC_ARG_ENABLE([parallel-netcdf],
  [AS_HELP_STRING([--enable-parallel-netcdf],
     [enable usage of the parallel features of NetCDF
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_parallel_netcdf=yes])],
  [enable_parallel_netcdf=no])

AC_ARG_ENABLE([cdi-pio],
  [AS_HELP_STRING([--enable-cdi-pio],
     [enable usage of the parallel features of CDI @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_cdi_pio=yes])],
  [enable_cdi_pio=no])

AC_ARG_ENABLE([sct],
  [AS_HELP_STRING([--enable-sct],
     [enable the SCT timer @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_sct=yes])],
  [enable_sct=no])

AC_ARG_ENABLE([yaxt],
  [AS_HELP_STRING([--enable-yaxt],
     [enable the YAXT data exchange @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_yaxt=yes])],
  [enable_yaxt=no])

AC_ARG_ENABLE([explicit-fpp],
  [AS_HELP_STRING([--enable-explicit-fpp],
     [enable explicit Fortran preprocessing @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_explicit_fpp=yes])],
  [enable_explicit_fpp=auto])

AC_ARG_ENABLE([serialization],
[  --enable-serialization  enable the Serialbox2 serialization. The value of
                          the argument must be one of the following
                          @<:@default=no@:>@:
                            read     enable READ mode
                            perturb  enable READ & PERTURB mode
                            create   enable CREATE mode
                            yes      alias for 'read'
                            no       disable serialization],
[AS_CASE([$enableval],
   [yes], [enable_serialization=read],
   [read|perturb|create|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-serialization='$enableval'; valid values are 'read', 'perturb', dnl
'create', 'yes' (same as 'read'), 'no'])])],
[enable_serialization=no])

AS_VAR_IF([enable_explicit_fpp], [auto],
  [AS_IF([test x"$enable_serialization" != xno],
     [AC_MSG_NOTICE([explicit Fortran preprocessing is enabled because the dnl
Serialbox2 serialization is requested: disable the explicit preprocessing dnl
if required (--disable-explicit-fpp)])
      enable_explicit_fpp=yes],
     [enable_explicit_fpp=no])])

AC_ARG_ENABLE([testbed],
  [AS_HELP_STRING([--enable-testbed],
     [enable ICON Testbed infrastructure @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_testbed=yes])],
  [enable_testbed=no])

AC_ARG_ENABLE([memory-tracing],
[  --enable-memory-tracing enable native dynamic memory tracing facility. The
                          value of the argument must be one of the following
                          @<:@default=no@:>@:
                            mtrace  enable tracing with mtrace (glibc)
                            yes     alias for 'mtrace'
                            no      disable memory tracing],
[AS_CASE([$enableval],
   [yes], [enable_memory_tracing=mtrace],
   [mtrace|no], [],
   [AC_MSG_ERROR([unexpected value for the argument dnl
--enable-memory-tracing='$enableval'; valid values are 'mtrace', 'yes' (same dnl
as 'mtrace'), 'no'])])],
[enable_memory_tracing=no])

m4_divert_once([HELP_ENABLE], [[
Optimization Features:]])

AC_ARG_ENABLE([loop-exchange],
  [AS_HELP_STRING([--enable-loop-exchange],
     [enable loop exchange @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_loop_exchange=yes])],
  [enable_loop_exchange=auto])

AS_VAR_IF([enable_loop_exchange], [auto],
  [AS_VAR_IF([icon_gpu_primary], [openacc],
     [enable_loop_exchange=no],
     [AC_MSG_NOTICE([loop exchange is enabled because no GPU support is dnl
requested: disable the loop exchange if required (--disable-loop-exchange)])
      enable_loop_exchange=yes])])

AC_ARG_ENABLE([dim-swap],
  [AS_HELP_STRING([--enable-dim-swap],
     [enable dimension swap @<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_dim_swap=yes])],
  [enable_dim_swap=auto])

AS_VAR_IF([enable_dim_swap], [auto],
  [AS_VAR_IF([icon_gpu_primary], [openacc],
     [AC_MSG_NOTICE([dimension swap is enabled because the GPU support is dnl
requested: disable the dimension swap if required (--disable-dim-swap)])
      enable_dim_swap=yes],
     [enable_dim_swap=no])])

AC_ARG_ENABLE([realloc-buf],
  [AS_HELP_STRING([--enable-realloc-buf],
     [enable reallocatable buffer in the communication @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno ], [enable_realloc_buf=yes])],
  [enable_realloc_buf=no])

AC_ARG_ENABLE([vectorized-lrtm],
  [AS_HELP_STRING([--enable-vectorized-lrtm],
     [enable the parallelization-invariant version of LRTM
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_vectorized_lrtm=yes])],
  [enable_vectorized_lrtm=no])

AC_ARG_ENABLE([mixed-precision],
  [AS_HELP_STRING([--enable-mixed-precision],
     [enable mixed precision dycore @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_mixed_precision=yes])],
  [enable_mixed_precision=no])

AC_ARG_ENABLE([intel-consistency],
  [AS_HELP_STRING([--enable-intel-consistency],
     [enable Intel compiler directives enforcing consistency
@<:@default=auto@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_intel_consistency=yes])],
  [enable_intel_consistency=auto])

AC_ARG_ENABLE([pgi-inlib],
  [AS_HELP_STRING([--enable-pgi-inlib],
     [enable PGI/NVIDIA cross-file function inlining via an inline library
@<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno && test x"$enableval" != xauto],
     [enable_pgi_inlib=yes])],
  [enable_pgi_inlib=no])

AC_ARG_ENABLE([nccl],
  [AS_HELP_STRING([--enable-nccl],
     [enable NCCL for communication @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_nccl=yes])],
  [enable_nccl=no])

AC_ARG_ENABLE([cuda-graphs],
  [AS_HELP_STRING([--enable-cuda-graphs],
     [enable CUDA graphs @<:@default=no@:>@])],
  [AS_IF([test x"$enableval" != xno], [enable_cuda_graphs=yes])],
  [enable_cuda_graphs=no])

m4_divert_once([HELP_ENABLE],
[  --enable-fcgroup-<NAME> enable Fortran compile group <NAME>, which can be
                          compiled with an extra set of compiler flags. <NAME>
                          can be any non-empty string consisting of capital
                          letters, digits and underscores. The Fortran
                          compiler flags to be appended to FCFLAGS when
                          compiling files of the group are read from the
                          environment variable ICON_<NAME>_FCFLAGS. The paths
                          to the directories and files included in the group
                          must be provided relative to the root source
                          directory. The value of the argument must be one of
                          the following @<:@default=no@:>@:
                            path1:path2  a colon-separated list of paths to
                                         the directories and files included in
                                         the group
                            yes          read paths to the directories and
                                         files of the group from the
                                         environment variable ICON_<NAME>_PATH
                            no           disable the group])
AC_SUBST([ICON_FCGROUP_ASSIGNMENTS], [''])
icon_fcgroups_unique=
icon_fcgroup_all_entries=
for icon_fcgroup in $icon_fcgroups; do
dnl Skip if the group has already been processed:
  AS_CASE([" $icon_fcgroups_unique "],
    [*"$icon_fcgroup"*], [continue],
    [AS_VAR_APPEND([icon_fcgroups_unique], [" $icon_fcgroup"])])
dnl Check that the name of the group consists of capital letters, digits and
dnl underscores:
  expr "X$icon_fcgroup" : "[.*[^_$as_cr_LETTERS$as_cr_digits]]" >/dev/null &&
    AC_MSG_ERROR([invalid compile group name '$icon_fcgroup': only capital dnl
letters and digits are allowed])
dnl Resolve trivial option values:
  enableval=AS_VAR_GET([enable_fcgroup_$icon_fcgroup])
  AS_CASE([$enableval],
    [no], [continue],
    [yes], [AS_VAR_COPY([enableval], [ICON_${icon_fcgroup}_PATH])])
  enableval=`AS_ECHO(["$enableval"]) | sed 's/::*/:/g;s/^://;s/:*$//'`
dnl
  icon_save_IFS=$IFS
  IFS=':'
  for icon_fcgroup_entry in $enableval; do
    IFS=$icon_save_IFS
    AS_CASE([$icon_fcgroup_entry],
      [[[\\/]]* | ?:[[\\/]]*],
      [AC_MSG_WARN([entry '$icon_fcgroup_entry' of the Fortran compile group dnl
'$icon_fcgroup' is ignored: not a relative path])
       continue])
    icon_fcgroup_path="$srcdir/$icon_fcgroup_entry"
    AS_IF(
      [test -d "$icon_fcgroup_path"],
      [icon_fcgroup_entry_dirname=$icon_fcgroup_path
       icon_fcgroup_entry_basename=],
      [test -f "$icon_fcgroup_path"],
      [icon_fcgroup_entry_dirname=`AS_DIRNAME(["$icon_fcgroup_path"])`
       icon_fcgroup_entry_basename=`AS_BASENAME(["$icon_fcgroup_path"])`],
      [AC_MSG_WARN([entry '$icon_fcgroup_entry' of the Fortran compile group dnl
'$icon_fcgroup' is ignored: not recognized as an existing directory or a file])
       continue])
    icon_fcgroup_entry_dirname=`cd "$icon_fcgroup_entry_dirname" >/dev/null dnl
2>&1 && pwd` || { AC_MSG_WARN([entry '$icon_fcgroup_entry' of the Fortran dnl
compile group '$icon_fcgroup' is ignored: unable to detect absolute path to dnl
'$icon_fcgroup_entry_dirname']); continue; }
    AS_CASE([$icon_fcgroup_entry_dirname],
      ["$ac_abs_confdir"[[\\/]]*], [],
      [AC_MSG_WARN([entry '$icon_fcgroup_entry' of the Fortran compile group dnl
'$icon_fcgroup' is ignored: not a subdirectory of the root source directory])
       continue])
    icon_fcgroup_entry_dirname=`AS_ECHO(["$icon_fcgroup_entry_dirname"]) | dnl
sed "s:^$ac_abs_confdir[[\\/]]::"`
    icon_fcgroup_entry_basename=`AS_ECHO(["$icon_fcgroup_entry_basename"]) | dnl
sed 's:[[^\.]]*$:%:'`
    icon_fcgroup_pattern=dnl
"${icon_fcgroup_entry_dirname}/${icon_fcgroup_entry_basename}"
    AS_CASE(["$icon_fcgroup_all_entries:"],
      [*":$icon_fcgroup_pattern:"*],
      [AC_MSG_WARN([entry '$icon_fcgroup_entry' of the Fortran compile group dnl
'$icon_fcgroup' is ignored: an entry with the same relative path (modulo dnl
file extension) has already been registered])
       continue],
      [AS_VAR_APPEND([icon_fcgroup_all_entries], [":$icon_fcgroup_pattern"])])
    AS_IF([test -n "$ICON_FCGROUP_ASSIGNMENTS"],
      [AS_VAR_APPEND([ICON_FCGROUP_ASSIGNMENTS], ["
"])])
    AS_VAR_APPEND([ICON_FCGROUP_ASSIGNMENTS],
      ["${icon_fcgroup_pattern}: dnl
ICON_FCFLAGS:= \$(ICON_${icon_fcgroup}_FCFLAGS)"])
  done
  IFS=$icon_save_IFS
done
dnl Newer versions of GNU Make (starting at least 4.1) give preference to the
dnl MOST SPECIALIZED matching pattern when performing pattern-specific variable
dnl assignment. Older versions (at least up to 3.81), however, prefer the LAST
dnl matching pattern. This means that, as long as we support the older version
dnl of GNU Make, we have to declare the most general assignments first. We can
dnl achieve that by sorting the list of assignments:
AS_IF([test -n "$ICON_FCGROUP_ASSIGNMENTS"],
  [ICON_FCGROUP_ASSIGNMENTS=`AS_ECHO(["$ICON_FCGROUP_ASSIGNMENTS"]) | sort`])

dnl Optional Packages:

AC_ARG_WITH([external-rte-rrtmgp],
     [AS_HELP_STRING([--with-external-rte-rrtmgp],
        [use external RTE+RRTMGP library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_rte_rrtmgp=yes])],
  [with_external_rte_rrtmgp=no])

AC_ARG_WITH([external-ecrad],
  [AS_HELP_STRING([--with-external-ecrad],
     [use external ECRAD library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_ecrad=yes])],
  [with_external_ecrad=no])

AC_ARG_WITH([external-sct],
  [AS_HELP_STRING([--with-external-sct],
     [use external SCT library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_sct=yes])],
  [with_external_sct=no])

AC_ARG_WITH([external-yaxt],
  [AS_HELP_STRING([--with-external-yaxt],
     [use external YAXT library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_yaxt=yes])],
  [with_external_yaxt=no])

AC_ARG_WITH([external-ppm],
  [AS_HELP_STRING([--with-external-ppm],
     [use external PPM library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_ppm=yes])],
  [with_external_ppm=no])

AC_ARG_WITH([external-cdi],
  [AS_HELP_STRING([--with-external-cdi],
     [use external CDI library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_cdi=yes])],
  [with_external_cdi=no])

AC_ARG_WITH([external-mtime],
  [AS_HELP_STRING([--with-external-mtime],
     [use external MTIME library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_mtime=yes])],
  [with_external_mtime=no])

AC_ARG_WITH([external-yac],
  [AS_HELP_STRING([--with-external-yac],
     [use external YAC library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_yac=yes])],
  [with_external_yac=no])

AC_ARG_WITH([external-tixi],
  [AS_HELP_STRING([--with-external-tixi],
     [use external TIXI library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno], [with_external_tixi=yes])],
  [with_external_tixi=no])

AC_ARG_WITH([external-fortran-support],
  [AS_HELP_STRING([--with-external-fortran-support],
     [use external FORTRAN-SUPPORT library @<:@default=no@:>@])],
  [AS_IF([test x"$withval" != xno],
     [with_external_fortran_support=yes])],
  [with_external_fortran_support=no])

############################# CODEBASE EXTENSIONS ##############################

AS_VAR_IF([enable_jsbach], [yes],
  [icon_subdir_srcdir="$srcdir/externals/jsbach"
   icon_subdir_unique_file="src/base/mo_jsb_base.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the land component JSBACH dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable the component dnl
(--disable-jsbach)])])])

AS_VAR_IF([enable_dace], [yes],
  [icon_subdir_srcdir="$srcdir/externals/dace_icon"
   icon_subdir_unique_file="src_for_icon/mo_fdbk_tables.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the DACE modules for data dnl
assimilation ($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable dnl
the component (--disable-dace)])])])

AS_VAR_IF([enable_emvorado], [yes],
  [icon_subdir_srcdir="$srcdir/externals/emvorado"
   icon_subdir_unique_file="src_emvorado/radar_organize.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the radar forward operator dnl
EMVORADO ($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable the dnl
operator (--disable-emvorado)])])
   AS_VAR_IF([enable_dace], [no],
     [AC_MSG_WARN([the radar forward operator EMVORADO is enabled but the dnl
DACE modules for data assimilation are disabled: EMVORADO will not be able dnl
to produce the FEEDBACK (fdbk) files for data assimilation])])])

AS_VAR_IF([enable_art], [yes],
  [icon_subdir_srcdir="$srcdir/externals/art"
   icon_subdir_unique_file="chem_init/mo_art_chem_init_utils.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the aerosols and reactive dnl
trace component ART ($icon_subdir_unique_file) in '$icon_subdir_srcdir': dnl
disable the component (--disable-art)])])])

############################## ENVIRONMENT CHECKS ##############################

ACX_BUILD_ENVIRONMENT

AC_CANONICAL_HOST
ACX_HOST_FQDN([], [:])
AC_SUBST([host_fqdn], [$acx_cv_host_fqdn])

###################### FORTRAN COMPILER COMMAND AND FLAGS ######################
AC_LANG([Fortran])
dnl Override the default flags set by AC_PROG_FC:
FCFLAGS=${FCFLAGS-}
dnl Save the flags to pass them later to the bundled libraries:
subdir_FCFLAGS=$FCFLAGS

dnl Make sure the configure checks are run with ICON_*FLAGS:
test -n "${ICON_LDFLAGS}" && AS_VAR_APPEND([LDFLAGS], [" ${ICON_LDFLAGS}"])
dnl We cannot literally append ICON_FCFLAGS to FCFLAGS now because it will be
dnl difficult to replace ICON_FCFLAGS with ICON_OCEAN_FCFLAGS afterwards.
dnl Instead, we make sure that all the checks are run with
dnl '$FCFLAGS $ICON_FCFLAGS':
m4_pushdef([AC_LANG(Fortran)],
  m4_bpatsubst(m4_dquote(m4_defn([AC_LANG(Fortran)])),
    [FCFLAGS],[FCFLAGS $ICON_FCFLAGS]))dnl

dnl Delay the cross-compilation check:
ACX_COMPILER_CROSS_CHECK_DELAY
AC_PROG_FC
dnl Make sure that the bundled libraries pick up the same FC if not set by the
dnl user on the command line:
export FC

dnl We put the following declarations here (i.e. after AC_PROG_FC) for better
dnl readability of the help message:
AC_ARG_VAR([ICON_FCFLAGS],
  [Fortran compiler flags to be appended to FCFLAGS when configuring ICON])dnl
AC_ARG_VAR([ICON_LDFLAGS],
  [Fortran compiler flags to be appended to LDFLAGS when configuring ICON])dnl

m4_divert_push([HELP_VAR])dnl
AS_HELP_STRING([ICON_<NAME>_FCFLAGS],
  [Fortran compiler flags to be appended to FCFLAGS instead of ICON_FCFLAGS when
compiling files of the compile group <NAME> (see --enable-fcgroup-<NAME>)
@<:@default=$ICON_FCFLAGS@:>@], [14])
AS_HELP_STRING([ICON_<NAME>_PATH],
  [colon-separated list of paths (relative to the root source directory) to
files and directories included in the compile group <NAME> (see
--enable-fcgroup-<NAME>)], [14])
m4_divert_pop([HELP_VAR])dnl

AC_SUBST([ICON_FCGROUP_FLAGS], [''])dnl
for icon_fcgroup in $icon_fcgroups_unique; do
  eval "icon_fcgroup_flags=\${ICON_${icon_fcgroup}_FCFLAGS-\$ICON_FCFLAGS}"
  AS_IF([test -n "$ICON_FCGROUP_FLAGS"],
    [AS_VAR_APPEND([ICON_FCGROUP_FLAGS], ["
"])])
  AS_VAR_APPEND([ICON_FCGROUP_FLAGS],
    ["ICON_${icon_fcgroup}_FCFLAGS= ${icon_fcgroup_flags}"])
done

AC_ARG_VAR([ICON_BUNDLED_FCFLAGS],
  [Fortran compiler flags to be appended to FCFLAGS when configuring the bundled
libraries @<:@default=$ICON_FCFLAGS@:>@])dnl
ICON_BUNDLED_FCFLAGS=${ICON_BUNDLED_FCFLAGS-$ICON_FCFLAGS}

ACX_FC_PP_SRCEXT([f90])
ACX_COMPILER_FC_VENDOR_SIMPLE
AC_SUBST([FC_VENDOR], [$acx_cv_fc_compiler_vendor])

AS_VAR_IF([enable_intel_consistency], [auto],
  [AS_VAR_IF([FC_VENDOR], [intel],
     [AC_MSG_NOTICE([Intel compiler directives enforcing consistency are dnl
enabled because $FC is recognized as the Intel Fortran compiler: disable dnl
the directives if required (--disable-intel-consistency)])
      enable_intel_consistency=yes],
     [enable_intel_consistency=no])])

AS_VAR_IF([enable_pgi_inlib], [auto],
  [AS_VAR_IF([FC_VENDOR], [portland],
     [AC_MSG_NOTICE([PGI/NVIDIA cross-file function inlining via an inline dnl
library is enabled because $FC is recognized as the PGI/NVIDIA compiler: dnl
disable the inlining if required (--disable-pgi-inlib)])
      enable_pgi_inlib=yes],
     [enable_pgi_inlib=no])])

ACX_COMPILER_FC_VERSION_SIMPLE
AC_SUBST([FC_VERSION], [$acx_cv_fc_compiler_version])

dnl Generate RPATH flags:
AS_VAR_IF([enable_rpaths], [yes],
  [ACX_SHLIB_FC_RPATH_FLAG
   ASX_EXTRACT_ARGS([icon_L_args], ["$LDFLAGS $LIBS"], ['-L@<:@ @:>@*'])
   for icon_L_arg in $icon_L_args; do
     ASX_VAR_APPEND_UNIQ([LDFLAGS], ["$acx_cv_fc_rpath_flag$icon_L_arg"], [' '])
   done
   ACX_SHLIB_RPATH_FLAGS_CHECK])

dnl Run the delayed cross-compilation check:
ACX_COMPILER_CROSS_CHECK_NOW

ACX_FC_LINE_LENGTH([unlimited])

ACX_FC_MODULE_IN_FLAG([AC_SUBST([FC_MOD_IN], [$acx_cv_fc_module_in_flag])])
ACX_FC_MODULE_OUT_FLAG([AC_SUBST([FC_MOD_OUT], [$acx_cv_fc_module_out_flag])])

ACX_FC_MODULE_NAMING(
  [AC_SUBST([FC_MOD_FILE_UPPER], [$acx_cv_fc_module_naming_upper])
   AC_SUBST([FC_MOD_FILE_EXT], [$acx_cv_fc_module_naming_ext])])

ACX_FC_INCLUDE_FLAG(
  [AC_SUBST([FC_INC_FLAG], [$acx_cv_fc_ftn_include_flag])])
ACX_FC_INCLUDE_ORDER(
  [FC_INC_ORDER=$acx_cv_fc_ftn_include_order],
  [FC_INC_ORDER=src,flg])
AC_SUBST([FC_INC_ORDER])

ACX_FC_INCLUDE_FLAG_PP(
  [AC_SUBST([FC_INC_FLAG_PP_f90], [$acx_cv_fc_pp_include_flag])])
ACX_FC_INCLUDE_ORDER_PP(
  [FC_INC_ORDER_PP_f90=$acx_cv_fc_pp_include_order],
  [FC_INC_ORDER_PP_f90=inc,flg])
AC_SUBST([FC_INC_ORDER_PP_f90])

ACX_LANG_MACRO_FLAG([AC_SUBST([FC_PP_DEF], [$acx_cv_fc_macro_flag])])

ACX_FC_ATTRIBUTE_CONTIGUOUS(
  [AS_VAR_APPEND([FCFLAGS],
  [" ${FC_PP_DEF}HAVE_FC_ATTRIBUTE_CONTIGUOUS"])], [:])

AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__ICON__"])

AS_VAR_IF([enable_atmo], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_ATMO__"])])

AS_VAR_IF([enable_les], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_LES__"])])

AS_VAR_IF([enable_upatmo], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_UPATMO__"])])

AS_VAR_IF([enable_ocean], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_OCEAN__"])])

AS_VAR_IF([enable_jsbach], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_JSBACH__"])])

AS_VAR_IF([enable_jsbach_hd], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_JSBACH_HD__"])],
  [AS_VAR_IF([enable_jsbach], [no],
     [AC_MSG_ERROR([unsupported combination of the configure options: you dnl
are trying to configure with the JSBACH HD scheme but the land component dnl
JSBACH is disabled: either disable the hydrological discharge scheme of the dnl
JSBACH land component (--disable-jsbach-hd) or enable the land component dnl
JSBACH (--enable-jsbach)])])])

AS_VAR_IF([enable_waves], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_WAVES__"])])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}YAC_coupling"])])

AS_VAR_IF([enable_aes], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_AES__"])])

AS_VAR_IF([enable_nwp], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_NWP__"])])

AS_VAR_IF([enable_ecrad], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__ECRAD"])])

AS_VAR_IF([enable_rte_rrtmgp], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_RTE_RRTMGP__"])])

AS_VAR_IF([enable_rttov], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__USE_RTTOV"])])

AS_VAR_IF([enable_dace], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__DACE__"])])

AS_VAR_IF([enable_emvorado], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}HAVE_RADARFWO ${FC_PP_DEF}TWOMOM_SB_NEW"])
   AS_VAR_IF([enable_dace], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NUDGING"])])])

AS_VAR_IF([enable_art], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}__ICON_ART ${FC_PP_DEF}__ART_TMP_IFDEF"])])

AS_VAR_IF([enable_art_gpl], [yes],
  [AS_VAR_IF([enable_art], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__ART_GPL"])],
     [AC_MSG_ERROR([unsupported combination of the configure options: you dnl
are trying to configure with the GPL-licensed code parts of the ART dnl
component but the ART component is disabled: either disable the the dnl
GPL-licensed code of the component (--disable-art-gpl) or enable the ART dnl
component (--enable-art)])])])

AS_VAR_IF([enable_comin], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_COMIN__"])])

AS_VAR_IF([enable_acm_license], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_ACM_LICENSE"])])

AS_VAR_IF([enable_mpi], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NOMPI"])])

AS_VAR_IF([enable_active_target_sync], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}HAVE_SLOW_PASSIVE_TARGET_ONESIDED"])])

AS_VAR_IF([enable_mpi_gpu], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__USE_G2G"])])

AS_VAR_IF([icon_gpu_secondary], [hip],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__HIP__"])])

AS_VAR_IF([enable_parallel_netcdf], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_PARALLEL_NETCDF"])])

AS_VAR_IF([enable_async_io_rma], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NO_ASYNC_IO_RMA"])])

AS_VAR_IF([enable_cdi_pio], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_CDI_PIO"])])

AS_VAR_IF([enable_sct], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__SCT__"])])

AS_VAR_IF([enable_yaxt], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_YAXT"])])

AS_IF([test x"$enable_serialization" != xno],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}SERIALIZE"])
   AS_CASE([$enable_serialization],
     [read], [AS_VAR_APPEND([FCFLAGS],
                [" ${FC_PP_DEF}SERIALIZE_READ_REFERENCE"])],
     [perturb], [AS_VAR_APPEND([FCFLAGS],
                   [" ${FC_PP_DEF}SERIALIZE_PERTURB_REFERENCE"])],
     [create], [AS_VAR_APPEND([FCFLAGS],
                  [" ${FC_PP_DEF}SERIALIZE_CREATE_REFERENCE"])])])

AS_VAR_IF([enable_testbed], [no],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__NO_ICON_TESTBED__"])])

AS_VAR_IF([enable_loop_exchange], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__LOOP_EXCHANGE"])])

AS_VAR_IF([enable_dim_swap], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__SWAPDIM"])])

AS_VAR_IF([enable_realloc_buf], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}__REALLOC_BUF"])])

AS_VAR_IF([enable_vectorized_lrtm], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}LRTM_FULL_VECTORIZATION"])])

AS_VAR_IF([enable_mixed_precision], [yes],
  [AS_VAR_APPEND([FCFLAGS],
     [" ${FC_PP_DEF}__MIXED_PRECISION"])])

AS_VAR_IF([enable_intel_consistency], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}IFORT_CONSISTENCY_ENFORCE"])])

AS_VAR_IF([enable_nccl], [yes],
  [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}USE_NCCL"])])

AS_VAR_IF([enable_cuda_graphs], [yes],
  [icon_warn=yes
   AS_CASE(["$FC_VENDOR,$FC_VERSION"],
     [portland,nv:*],
     [AX_COMPARE_VERSION([$FC_VERSION" | cut -d: -f2"], [ge], [23.3],
        [icon_warn=no])])
   AS_VAR_IF([icon_warn], [yes],
     [AC_MSG_WARN([Fortran compiler is not recognized as NVIDIA HPC Fortran dnl
23.3 or newer: it is strongly recommended to disable CUDA graphs dnl
(--disable-cuda-graphs)])])
   AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}ICON_USE_CUDA_GRAPH"])])

######################### C COMPILER COMMAND AND FLAGS #########################
AC_LANG([C])
dnl Override the default flags set by AC_PROG_CC:
CFLAGS=${CFLAGS-}
dnl Save the flags to pass them later to the bundled libraries:
subdir_CFLAGS=$CFLAGS
subdir_CPPFLAGS=$CPPFLAGS

dnl Make sure the configure checks are run with ICON_CFLAGS:
test -n "${ICON_CFLAGS}" && AS_VAR_APPEND([CFLAGS], [" ${ICON_CFLAGS}"])

dnl Fix a bug in Autoconf 2.70 (see https://savannah.gnu.org/support/?110396):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [0],
  [m4_pushdef([_AC_PROG_CC_C99],
        m4_bpatsubst(m4_dquote(m4_defn([_AC_PROG_CC_C99])),
          [\[ac_c_conftest_c89_program\]], [[ac_c_conftest_c99_program]]))dnl
   AC_PROG_CC
   m4_popdef([_AC_PROG_CC_C99])],
  [AC_PROG_CC])
dnl Expand the legacy macro for older version of Autoconf (2.69 and older):
m4_if(m4_version_compare(m4_defn([AC_AUTOCONF_VERSION]), [2.70]), [-1],
  [AC_PROG_CC_C99])
AS_VAR_IF([ac_cv_prog_cc_c99], [no],
  [AC_MSG_FAILURE([C compiler does not support ISO C99])])
dnl Make sure that the bundled libraries pick up the same CC if not set by the
dnl user on the command line:
export CC

dnl We put the following declaration here (i.e. after AC_PROG_CC) for better
dnl readability of the help message:
AC_ARG_VAR([ICON_CFLAGS],
  [C compiler flags to be appended to CFLAGS when configuring ICON])dnl

AC_ARG_VAR([ICON_BUNDLED_CFLAGS],
  [C compiler flags to be appended to CFLAGS when configuring the bundled
libraries @<:@default=$ICON_CFLAGS@:>@])dnl
ICON_BUNDLED_CFLAGS=${ICON_BUNDLED_CFLAGS-$ICON_CFLAGS}

ACX_COMPILER_CC_VENDOR_SIMPLE
ACX_COMPILER_CC_VERSION_SIMPLE

dnl AC_PROG_CC, AC_PROG_CC_C99, and AC_INCLUDES_DEFAULT do check for some of the
dnl following header files but the actual list of files they check for might
dnl vary from version to version. Therefore, we check for all the headers that
dnl are used in the code unconditionally:
AC_CHECK_HEADERS([stdlib.h stdio.h math.h time.h ctype.h string.h assert.h ]dnl
  [stdbool.h dirent.h errno.h fcntl.h sys/stat.h sys/types.h unistd.h],
  [], [AC_MSG_FAILURE([unable to find header file $ac_header])],
  [[// override default includes]])

######################## C++ COMPILER COMMAND AND FLAGS ########################

dnl Override the default flags set by AC_PROG_CXX:
CXXFLAGS=${CXXFLAGS-}
dnl Save the flags to pass them later to the bundled libraries:
subdir_CXXFLAGS=$CXXFLAGS

dnl We do not expand AC_PROG_CXX because ICON itself does not need a C++
dnl compiler yet. However, some of the bundled libraries might need it.
dnl Therefore we simply declare the following variables to document them in the
dnl help message:
AC_ARG_VAR([CXX], [C++ compiler command])dnl
AC_ARG_VAR([CXXFLAGS], [C++ compiler flags])dnl
dnl Make sure that the bundled libraries pick up the same CXX if not set by the
dnl user on the command line:
export CXX

AC_ARG_VAR([ICON_BUNDLED_CXXFLAGS],
  [C++ compiler flags to be appended to CXXFLAGS when configuring the bundled
libraries])dnl

###################### SECONDARY GPU COMPILERS AND FLAGS #######################

dnl We do not use AS_CASE to prevent checking for CUDACXX and HIPCXX when GPU
dnl support is disabled:
case $icon_gpu_secondary in #(
  cuda) :
dnl Honor NVCC and NVCFLAGS for backwards compatibility:
    AS_IF([${CUDACXX+false} : && ${NVCC+:} false],
      [CUDACXX=$NVCC
       AC_MSG_WARN(
         [variable NVCC is deprecated: use variable CUDACXX instead])])
    AS_IF([${CUDAFLAGS+false} : && ${NVCFLAGS+:} false],
      [CUDAFLAGS=$NVCFLAGS
       AC_MSG_WARN(
         [variable NVCFLAGS is deprecated: use variable CUDAFLAGS instead])])
    ACX_PROG_CUDACXX
    AS_CASE([$acx_prog_cudacxx_works],
      [no],
      [AC_MSG_FAILURE([unable to find a working CUDA Compiler: either switch dnl
the secondary GPU compiler (e.g. --enable-gpu=${icon_gpu_primary}+hip) or dnl
disable GPU support (--disable-gpu), or set CUDACXX accordingly])],
      [basic],
      [AC_MSG_FAILURE([unable to detect CUDA compiler flag needed to enable
support for the C++11 standard: either switch the secondary GPU compiler dnl
(e.g. --enable-gpu=${icon_gpu_primary}+hip) or disable GPU support dnl
(--disable-gpu), or set CUDAFLAGS accordingly])],
      [CUDAFLAGS=dnl
"AS_VAR_GET([acx_cv_prog_cudacxx_${acx_prog_cudacxx_works}_flag]) $CUDAFLAGS"])
    AC_LANG_PUSH([Fortran])
    ACX_LANG_CUDA_COMPATIBLE([],
      [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
CUDA compiler: a possible reason is that it cannot find CUDA runtime dnl
libraries (e.g. '-lcudart'): set LDFLAGS and LIBS accordingly])])
    ACX_LANG_CUDA_COMPATIBLE_STDCXX([],
      [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
CUDA compiler that require C++ standard library: a possible reason is that dnl
it cannot find C++ standard library (i.e. '-lstdc++'): set LDFLAGS and LIBS dnl
accordingly])])
    AC_LANG_POP([Fortran])
    HIPCXX=; HIPFLAGS= ;; #(
  hip) :
    ACX_PROG_HIPCXX
    AS_CASE([$acx_prog_hipcxx_works],
      [no],
      [AC_MSG_FAILURE([unable to find a working HIP Compiler: either switch dnl
the secondary GPU compiler (e.g. --enable-gpu=${icon_gpu_primary}+cuda) or dnl
disable GPU support (--disable-gpu), or set HIPCXX accordingly])],
      [basic],
      [AC_MSG_FAILURE([unable to detect HIP compiler flag needed to enable dnl
support for the C++11 standard: either switch the secondary GPU compiler dnl
(e.g. --enable-gpu=${icon_gpu_primary}+cuda) or disable GPU support dnl
(--disable-gpu), or set HIPFLAGS accordingly])],
      [HIPFLAGS=dnl
"AS_VAR_GET([acx_cv_prog_hipcxx_${acx_prog_hipcxx_works}_flag]) $HIPFLAGS"])
    AC_LANG_PUSH([Fortran])
    ACX_LANG_HIP_COMPATIBLE([],
      [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
HIP compiler: a possible reason is that it cannot find HIP runtime dnl
libraries (e.g. '-lamdhip64'): set LDFLAGS and LIBS accordingly])])
    ACX_LANG_HIP_COMPATIBLE_STDCXX([],
      [AC_MSG_FAILURE([Fortran compiler cannot link objects compiled with dnl
HIP compiler that require C++ standard library: a possible reason is that dnl
it cannot find C++ standard library (i.e. '-lstdc++'): set LDFLAGS and LIBS dnl
accordingly])])
    AC_LANG_POP([Fortran])
    CUDACXX=; CUDAFLAGS= ;; #(
  *) :
    CUDACXX=; CUDAFLAGS=
    HIPCXX=; HIPFLAGS= ;;
esac

################################ LIBRARY CHECKS ################################

dnl One of the topological orderings of the ICON dependency graph
dnl (run ./utils/icon_sorted_deps.sh):
dnl
dnl '-licon_tixi', '-lyac', '-lmtime', '-lSerialboxFortranShared',
dnl '-lcdi_f2003 -lcdipio -lcdi', '-lscalesppmcore', '-leccodes_f90',
dnl '-lyaxt -lyaxt_c', '-lsct', '-lrttov_other -lrttov_emis_atlas (etc.)',
dnl '-llapack -lblas (or other)', '-lradiation -lifsrrtm -lutilities -lifsaux',
dnl '-lrrtmgp -lrte', '-lnetcdff', '-lhdf5hl_fortran -lhdf5_fortran', '-lz',
dnl '-lgomp (or other)', '-lmpifort -lmpi (or other)'
dnl
dnl We check these dependencies in the reverse order.

dnl Currently, we assume that the Fortran compiler implicitly links to libm and
dnl libc, which are required by ICON itself and some of the bundled libraries
dnl (including the mandatory ones).

dnl Most of the following checks are for Fortran compiler:
AC_LANG([Fortran])

dnl Additional macro definitions to aid the dependency generator,
dnl e.g. compiler intrinsic macros:
AC_SUBST([DEPGEN_FCFLAGS], [''])

dnl A list of library files (i.e. *.a files) of the enabled bundled libraries:
AC_SUBST([BUNDLED_LIBFILES], [''])

dnl A list of include search flags targeting directories of the bundled
dnl libraries:
bundled_includes=

dnl Get absolute paths to the source and build directories:
ASX_SRCDIRS([.])
icon_abs_top_srcdir=$ac_abs_top_srcdir
icon_abs_top_builddir=$ac_abs_top_builddir

AS_VAR_IF([enable_memory_tracing], [mtrace],
  [AC_CACHE_CHECK([whether C mtrace function is available via the Fortran dnl
ISO_C_BINDING interface], [icon_cv_have_mtrace],
     [icon_cv_have_mtrace=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      implicit none
      interface
      subroutine mtrace() bind(c, name="mtrace")
      end subroutine mtrace
      end interface
      call mtrace]])],
        [icon_cv_have_mtrace=yes])])
   AS_VAR_IF([icon_cv_have_mtrace], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}ICON_MEMORY_TRACING=1"])],
     [AC_MSG_FAILURE([unable to link a test program using mtrace function dnl
via the Fortran ISO_C_BINDING interface: either disable the native dynamic dnl
memory tracing facility (--disable-memory-tracing) or set LDFLAGS and LIBS dnl
accordingly])])])

AS_IF([test "x$enable_mpi" = xyes || test "x$enable_yaxt" = xyes],
  [AC_CACHE_CHECK([whether the Fortran 90 interface of MPI library is dnl
available], [icon_cv_have_mpi],
     [icon_cv_have_mpi=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
       use mpi, only: mpi_init
       implicit none
       integer :: error
       call mpi_init(error)]])],
        [icon_cv_have_mpi=yes])])
   AS_VAR_IF([icon_cv_have_mpi], [yes],
     [AS_IF([test x"$enable_dace" = xyes || test x"$enable_rttov" = xyes],
        [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HAVE_MPI_MOD"])])],
     [AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interface of MPI library: either disable MPI support (--disable-mpi) and the dnl
YAXT data exchange (--disable-yaxt) or set FC to the MPI Fortran compiler dnl
wrapper, or set FCFLAGS, LDFLAGS and LIBS accordingly])])
dnl ######################################################################
dnl test for MPI standard implementation version
   AC_MSG_CHECKING([for the version of the supported MPI standard])
   AS_IF([AS_VAR_TEST_SET([icon_cv_mpi_version]) &&
AS_VAR_TEST_SET([icon_cv_mpi_subversion])],
     [AS_ECHO_N(["(cached) "]) >&AS_MESSAGE_FD],
     [AS_VAR_SET_IF([icon_cv_mpi_version], [],
        [icon_cv_mpi_version=unknown
         ACX_FORTRAN_CHECK_INTEGER_EXPRESSION([], [MPI_VERSION], [use mpi],
           [0], [100], [1], [icon_cv_mpi_version], [:])])
      AS_VAR_SET_IF([icon_cv_mpi_subversion], [],
        [icon_cv_mpi_subversion=unknown
         AS_IF([test x"$icon_cv_mpi_version" != xunknown],
           [ACX_FORTRAN_CHECK_INTEGER_EXPRESSION([], [MPI_SUBVERSION],
              [use mpi], [0], [100], [1], [icon_cv_mpi_subversion], [:])])])])
   AS_IF([test x"$icon_cv_mpi_version" != xunknown && dnl
test x"$icon_cv_mpi_subversion" != xunknown],
     [AC_MSG_RESULT([$icon_cv_mpi_version.$icon_cv_mpi_subversion])
      AS_VAR_APPEND([FCFLAGS],
        [" ${FC_PP_DEF}ICON_MPI_VERSION=$icon_cv_mpi_version dnl
${FC_PP_DEF}ICON_MPI_SUBVERSION=$icon_cv_mpi_subversion"])
      icon_mpi_version=$icon_cv_mpi_version
      icon_mpi_subversion=$icon_cv_mpi_subversion],
     [AC_MSG_RESULT([unknown])
      icon_mpi_version=unknown; icon_mpi_subversion=unknown])
dnl ######################################################################
   AS_VAR_IF([enable_mpi_rget], [auto],
     [AS_IF([test x"$icon_mpi_version" != xunknown && dnl
test "$icon_mpi_version" -ge 3 >/dev/null 2>&1],
        [AC_MSG_NOTICE([usage of the MPI_Rget routine is enabled because the dnl
provided MPI library supports standard 3.0 or higher: disable the usage if dnl
required (--disable-mpi-rget)])
         enable_mpi_rget=yes],
        [enable_mpi_rget=no])])
   AS_VAR_IF([enable_mpi_rget], [no],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NO_MPI_RGET"])])])
dnl ######################################################################
   AC_CACHE_CHECK(
     [for interfaces for MPI routines with choice arguments],
     [icon_cv_have_mpi_choice],
     [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
       use mpi, only: mpi_send, mpi_integer, mpi_comm_world
       implicit none
       integer :: ierror
       call mpi_send(0, 1, mpi_integer, 1, 0, mpi_comm_world, ierror)]])],
       [icon_cv_have_mpi_choice=yes],
       [icon_cv_have_mpi_choice=no])])
    AS_VAR_IF([icon_cv_have_mpi_choice], [no],
      [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NO_MPI_CHOICE_ARG"])])
dnl ######################################################################
   AC_CACHE_CHECK(
     [whether the MPI_Alloc_mem interface accepts a TYPE(C_PTR) argument],
     [icon_cv_have_mpi_cptr],
     [AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
       use, intrinsic :: iso_c_binding, only: c_ptr
       use mpi, only: mpi_alloc_mem, mpi_address_kind
       implicit none
       type(c_ptr) :: ptr
       integer :: ierror
       call mpi_alloc_mem(4_mpi_address_kind, 0, ptr, ierror)]])],
       [icon_cv_have_mpi_cptr=yes],
       [icon_cv_have_mpi_cptr=no])])
    AS_VAR_IF([icon_cv_have_mpi_cptr], [no],
      [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NO_MPI_CPTR_ARG"])])
dnl ######################################################################

dnl The configure scripts of the most of the bundled libraries expect
dnl MPI_LAUNCH to be set in the form of an absolute path, which we try to
dnl accommodate when required and if possible:
AS_IF([test -n "$MPI_LAUNCH"],
  [AS_CASE([$MPI_LAUNCH],
     [[[\\/]]* | ?:[[\\/]]*], [],
     [ACX_PROG_SEARCH_ABSPATH([$MPI_LAUNCH],
        [MPI_LAUNCH="$acx_prog_search_abspath"],
        [AC_MSG_WARN([unable to find the absolute path to $MPI_LAUNCH: the dnl
bundled libraries might pick other MPI launcher commands for their tests])])])])

AS_VAR_IF([enable_openmp], [yes],
  [ACX_LANG_OPENMP_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openmp_flag"])
      ACX_LANG_MACRO_CHECK_VALUE([_OPENMP],
         [201511 201307 201107 200805 200505 200011 199911 199710])
      AS_VAR_IF([acx_macro_value], [unknown],
        [AC_MSG_WARN([unable to detect supported OpenMP version])
         AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENMP"])],
        [AS_IF([test "$acx_macro_value" -lt 201307 >&AS_MESSAGE_LOG_FD 2>&1],
           [AC_MSG_WARN([Fortran compiler does not support OpenMP 4.0 dnl
(201307): it is strongly recommended to disable OpenMP support dnl
(--disable-openmp)])])
         AS_VAR_APPEND([DEPGEN_FCFLAGS],
           [" ${FC_PP_DEF}_OPENMP=$acx_macro_value"])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENMP])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenMP support by default: set dnl
FCFLAGS accordingly to disable it])])])

AS_VAR_IF([icon_gpu_primary], [openacc],
  [ACX_LANG_OPENACC_FLAG(
     [AS_VAR_APPEND([FCFLAGS], [" $acx_cv_fc_openacc_flag"])
      ACX_LANG_MACRO_CHECK_VALUE([_OPENACC],
         [201711 201510 201308 201306 201111])
      AS_VAR_IF([acx_macro_value], [unknown],
        [AC_MSG_WARN([unable to detect supported OpenACC version])
         AS_VAR_APPEND([DEPGEN_FCFLAGS], [" ${FC_PP_DEF}_OPENACC"])],
        [AS_IF([test "$acx_macro_value" -lt 201711 >&AS_MESSAGE_LOG_FD 2>&1],
           [AC_MSG_WARN([Fortran compiler does not support OpenACC 2.6 dnl
(201711): it is strongly recommended to disable GPU support (--disable-gpu)])])
      AS_VAR_APPEND([DEPGEN_FCFLAGS],
        [" ${FC_PP_DEF}_OPENACC=$acx_macro_value"])])])],
  [ACX_LANG_MACRO_CHECK_DEFINED([_OPENACC])
   AS_VAR_IF([acx_macro_defined], [yes],
     [AC_MSG_ERROR([Fortran compiler enables OpenACC support by default: dnl
either enable GPU support (--enable-gpu) or set FCFLAGS accordingly to dnl
disable OpenACC])])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_CACHE_CHECK([whether C ZLIB library is available via the Fortran dnl
ISO_C_BINDING interface], [icon_cv_have_zlib],
     [icon_cv_have_zlib=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      implicit none
      interface
      subroutine gzOpen() bind(c, name="gzopen")
      end subroutine gzOpen
      end interface
      call gzOpen]])],
        [icon_cv_have_zlib=yes])])
   AS_VAR_IF([icon_cv_have_zlib], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}WITH_ZLIB"])],
     [AC_MSG_FAILURE([unable to link a test program using ZLIB library via dnl
the Fortran ISO_C_BINDING interface: either disable the radar forward dnl
operator EMVORADO (--disable-emvorado) or set LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_CACHE_CHECK([whether the Fortran 90 low- and high-level interfaces of dnl
HDF5 library are available], [icon_cv_have_hdf5],
     [icon_cv_have_hdf5=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use hdf5, only: h5open_f, &
                    & h5fopen_f, &
                    & hid_t, &
                    & H5F_ACC_RDONLY_F
      use h5lt, only: h5ltget_attribute_string_f
      implicit none
      integer :: ierr
      integer(kind=hid_t) :: file_id
      character(len=500) :: catt_value
      call h5open_f(ierr)
      call h5fopen_f('conftest.h5', H5F_ACC_RDONLY_F, &
                   & file_id, ierr)
      call h5ltget_attribute_string_f(file_id, &
                                    & 'group', 'attr', &
                                    & catt_value, ierr)]])],
     [icon_cv_have_hdf5=yes])])
   AS_VAR_IF([icon_cv_have_hdf5], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}HDF5_RADAR_INPUT"])],
     [AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
low- and high-level interfaces of HDF5 library: set FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])])

AC_CACHE_CHECK([whether the Fortran 77 interface of NetCDF library is dnl
available], [icon_cv_have_netcdf77],
  [icon_cv_have_netcdf77=no
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      implicit none
      include 'netcdf.inc'
      integer :: status, ncid
      status = nf_open('conftest.nc', nf_nowrite, ncid)]])],
  [icon_cv_have_netcdf77=yes])])
AS_VAR_IF([icon_cv_have_netcdf77], [no],
  [AC_MSG_FAILURE([unable to link a test program using the Fortran 77 dnl
interface of NetCDF library: set FCFLAGS, LDFLAGS and LIBS accordingly])])
dnl
AS_VAR_IF([enable_parallel_netcdf], [yes],
  [AC_CACHE_CHECK([whether the provided version of NetCDF is built with dnl
parallel support], [icon_cv_have_netcdf_parallel],
     [icon_cv_have_netcdf_parallel=unknown
dnl The only known way to check whether a NetCDF-Fortran is built with parallel
dnl support is to build and run a test program. However, we need to run an MPI
dnl program. There is no support for running an MPI program without an MPI
dnl launcher on some platforms. Moreover, MPI launchers often return invalid
dnl exit codes. Therefore, we build a test program that prints
dnl 'conftest_parallel_yes' if NetCDF-Fortran library is built with parallel
dnl support, otherwise, the output is 'conftest_parallel_no'. We launch the
dnl program only if MPI_LAUNCH is set to a non-empty value and fail the
dnl configuration only if the output of the test program contains
dnl 'conftest_parallel_no' (some MPI implementations print warning messages to
dnl the standard output stream, therefore we cannot compare the output
dnl literally). If the output contains, 'conftest_parallel_yes', the test is
dnl considered successful. In all other cases, we issue a warning. Note that we
dnl do not run a proper check whether MPI_LAUNCH is a valid interactive MPI
dnl launcher command (as some of the bundled libraries do) because the main
dnl purpose of the test is to bring the user's attention to the fact that
dnl something is wrong, and if MPI_LAUNCH command is invalid, we will not get
dnl the expected output from the test program and issue a warning.
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mpi, only: mpi_init, mpi_finalize, &
                   & MPI_COMM_WORLD, MPI_INFO_NULL
      implicit none
      include 'netcdf.inc'
      integer :: result, tmp
      call mpi_init(tmp)
      result = nf_create_par('conftest.nc', &
                           & NF_MPIIO, &
                           & MPI_COMM_WORLD, &
                           & MPI_INFO_NULL, &
                           & tmp)
      if (result == NF_ENOPAR) then
        write(*,"(a)") "conftest_parallel_no"
      else
        write(*,"(a)") "conftest_parallel_yes"
      endif
      call mpi_finalize(tmp)]])],
        [icon_exec_result=dnl
`${MPI_LAUNCH:-false} -n 1 ./conftest$ac_exeext 2>/dev/null`
         rm -f conftest.nc
         AS_IF(
           [echo "$icon_exec_result" | dnl
grep 'conftest_parallel_no' >/dev/null 2>&1],
           [icon_cv_have_netcdf_parallel=no],
           [echo "$icon_exec_result" | dnl
grep 'conftest_parallel_yes' >/dev/null 2>&1],
           [icon_cv_have_netcdf_parallel=yes])])])
   AS_IF(
     [test x"$icon_cv_have_netcdf_parallel" = xno],
     [AC_MSG_FAILURE([the provided version of NetCDF does not support the dnl
parallel features: either disable usage of the parallel features of NetCDF dnl
(--disable-parallel-netcdf) or set FCFLAGS, LDFLAGS and LIBS accordingly])],
     [test x"$icon_cv_have_netcdf_parallel" = xunknown],
     [AC_MSG_WARN([unable to check whether the provided version of NetCDF dnl
supports the parallel features])])])
dnl
AC_CACHE_CHECK([whether the Fortran 90 interface of NetCDF library is dnl
available], [icon_cv_have_netcdf90],
  [icon_cv_have_netcdf90=no
   AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use netcdf, only: nf90_open, nf90_nowrite
      implicit none
      integer :: status, ncid
      status = nf90_open('conftest.nc', nf90_nowrite, ncid)]])],
  [icon_cv_have_netcdf90=yes])])
AS_VAR_IF([icon_cv_have_netcdf90], [yes],
  [AS_VAR_IF([enable_emvorado], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}NETCDF"])])],
  [AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interface of NetCDF library: set FCFLAGS, LDFLAGS and LIBS accordingly])])

AS_VAR_IF([enable_rte_rrtmgp], [yes],
  [AS_VAR_IF([with_external_rte_rrtmgp], [yes],
     [AC_CACHE_CHECK([whether Fortran RTE+RRTMGP library is available],
        [icon_cv_have_rte_rrtmgp],
        [icon_cv_have_rte_rrtmgp=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_optical_props, only: ty_optical_props_1scl, &
                                & ty_optical_props_2str
      use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
      use mo_source_functions, only: ty_source_func_lw
      use mo_rte_kind, only: wp
      use mo_fluxes, only: ty_fluxes_broadband
      use mo_rte_lw, only: rte_lw
      use mo_rte_sw, only: rte_sw
      use mo_gas_concentrations, only: ty_gas_concs
      implicit none
      type(ty_optical_props_1scl) :: optical_props_lw
      type(ty_optical_props_2str) :: optical_props_sw
      logical :: top_at_1
      type(ty_source_func_lw) :: source
      real(wp), dimension(:,:), allocatable :: sfc_emis_spec, &
                                             & toa_flux, &
                                             & sfc_alb_spec
      type(ty_fluxes_broadband) :: fluxes
      real(wp), dimension(:), allocatable :: mu0
      type(ty_gas_concs), dimension(:), allocatable :: gas_conc_array
      real(wp), dimension(:,:,:), allocatable :: gas_conc_temp_3d
      character(len=128) :: error_msg

      error_msg = rte_lw(optical_props_lw, top_at_1, source, &
                       & sfc_emis_spec, fluxes, n_gauss_angles = 1)
      error_msg = rte_sw(optical_props_sw, top_at_1, mu0, &
                       & toa_flux, sfc_alb_spec, sfc_alb_spec, fluxes)
      error_msg = gas_conc_array(1)%set_vmr('h2o',gas_conc_temp_3d(:,:,1))]])],
           [icon_cv_have_rte_rrtmgp=yes])])
      AS_VAR_IF([icon_cv_have_rte_rrtmgp], [no],
        [AC_MSG_FAILURE([unable to link a test program using RTE+RRTMGP dnl
library: either disable usage of the RTE+RRTMGP toolbox for radiation dnl
calculations (--disable-rte-rrtmgp) or use the bundled version of RTE+RRTMGP dnl
library (--without-external-rte-rrtmgp), or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/rte-rrtmgp"
      icon_subdir_unique_file="rte/mo_rte_sw.F90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the RTE+RRTMGP library dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable usage of dnl
the RTE+RRTMGP toolbox for radiation calculations (--disable-rte-rrtmgp) or dnl
provide an external version of RTE+RRTMGP library dnl
(--with-external-rte-rrtmgp) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      AC_ARG_VAR([ICON_RTE_RRTMGP_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of RTE+RRTMGP library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_RTE_RRTMGP_FCFLAGS=${ICON_RTE_RRTMGP_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/rte-rrtmgp], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/rte-rrtmgp],
        [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
        [NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([gpu])])
      ACX_SUBDIR_APPEND_ARGS([externals/rte-rrtmgp],
        [FCFLAGS="$subdir_FCFLAGS $ICON_RTE_RRTMGP_FCFLAGS"],
        [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
        [--enable-gpu="$icon_gpu_primary"],
        [--enable-silent-rules="$enable_silent_rules"])
      AS_VAR_APPEND([bundled_includes],
        ["${FC_MOD_IN}externals/rte-rrtmgp/mod "])
      BUNDLED_LIBFILES="externals/rte-rrtmgp/librte.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/rte-rrtmgp/librrtmgp.a $BUNDLED_LIBFILES"])])

AS_VAR_IF([enable_ecrad], [yes],
  [AS_VAR_IF([with_external_ecrad], [yes],
     [AC_CACHE_CHECK([whether Fortran ECRAD library is available],
        [icon_cv_have_ecrad],
        [icon_cv_have_ecrad=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use radiation_config, only: config_type
      use radiation_interface, only: setup_radiation
      implicit none
      type(config_type) :: conf
      call setup_radiation(conf)]])],
           [icon_cv_have_ecrad=yes])])
      AS_VAR_IF([icon_cv_have_ecrad], [no],
        [AC_MSG_FAILURE([unable to link a test program using ECRAD library: dnl
either disable usage of the ECMWF radiation scheme (--disable-ecrad) or use dnl
the bundled version of ECRAD library (--without-external-ecrad), or set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/ecrad"
      icon_subdir_unique_file="radiation/radiation_flux.F90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the the ECRAD library dnl
($icon_subdir_unique_file) in '$icon_subdir_srcdir': either disable usage of dnl
the ECMWF radiation scheme (--disable-ecrad) or provide an external version dnl
of ECRAD library (--with-external-ecrad) by setting FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])
      AC_ARG_VAR([ICON_ECRAD_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of ECRAD library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_ECRAD_FCFLAGS=${ICON_ECRAD_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/ecrad], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/ecrad],
        [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
        [NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([single-precision])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([pgi-inlib])])
      ACX_SUBDIR_APPEND_ARGS([externals/ecrad],
        [FCFLAGS="$subdir_FCFLAGS $ICON_ECRAD_FCFLAGS"],
        [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
        [--enable-silent-rules="$enable_silent_rules"],
        [--enable-pgi-inlib="$enable_pgi_inlib"])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/ecrad/mod "])
      BUNDLED_LIBFILES="externals/ecrad/libifsaux.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/ecrad/libutilities.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/ecrad/libifsrrtm.a $BUNDLED_LIBFILES"
      BUNDLED_LIBFILES="externals/ecrad/libradiation.a $BUNDLED_LIBFILES"])])

AC_CACHE_CHECK([whether Fortran BLAS library is available], [icon_cv_have_blas],
  [icon_cv_have_blas=no
   AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([dgemm])], [icon_cv_have_blas=yes])])
AS_VAR_IF([icon_cv_have_blas], [no],
  [AC_MSG_FAILURE([unable to link a test program using BLAS library: set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])

AC_CACHE_CHECK([whether Fortran LAPACK library is available],
  [icon_cv_have_lapack],
  [icon_cv_have_lapack=no
   AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY([dgeev])], [icon_cv_have_lapack=yes])])
AS_VAR_IF([icon_cv_have_lapack], [no],
  [AC_MSG_FAILURE([unable to link a test program using LAPACK library: set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])

AS_VAR_IF([enable_rttov], [yes],
  [AC_CACHE_CHECK([whether Fortran RTTOV library is available],
     [icon_cv_have_rttov],
     [icon_cv_have_rttov=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use rttov_types, only: rttov_coefs
      implicit none
#include "rttov_nullify_coef.interface"
      type(rttov_coefs),target:: coefs(1)
      call rttov_nullify_coef(coefs(1)%coef)]])],
        [icon_cv_have_rttov=yes])])
   AS_VAR_IF([icon_cv_have_rttov], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}_RTTOV_VERSION=13"])],
     [AC_MSG_FAILURE([unable to link a test program using RTTOV library: dnl
either disable usage of the radiative transfer model for TOVS dnl
(--disable-rttov) or set FCFLAGS, LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([enable_sct], [yes],
  [AS_VAR_IF([with_external_sct], [yes],
     [AC_CACHE_CHECK([whether the Fortran interface of SCT library is dnl
available], [icon_cv_have_sct],
        [icon_cv_have_sct=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use sct, only: sct_init
      implicit none
      call sct_init(timer_max=512)]])],
           [icon_cv_have_sct=yes])])
      AS_VAR_IF([icon_cv_have_sct], [no],
        [AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of SCT library: either disable the SCT timer (--disable-sct) or dnl
use the bundled version of SCT library (--without-external-sct), or set dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/sct"
      icon_subdir_unique_file="src/sct_reporter.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
SCT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the SCT timer (--disable-sct) or provide an external dnl
version of SCT library (--with-external-sct) by setting FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])
      AC_ARG_VAR([ICON_SCT_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of SCT library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_SCT_FCFLAGS=${ICON_SCT_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_SCT_CFLAGS],
        [C compiler flags to be appended to CFLAGS when configuring the bundled
version of SCT library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_SCT_CFLAGS=${ICON_SCT_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/sct], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/sct],
        [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
        [MPICC=*|MPIFC=*|MPI_LAUNCH=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([static])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([shared])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([mpi])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([mpirun-check])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([openmp])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([rpath])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([papi])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([hdf5])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([libpapi-prefix])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([libhdf5-prefix])])
      ACX_SUBDIR_APPEND_ARGS([externals/sct],
        [CFLAGS="$subdir_CFLAGS $ICON_SCT_CFLAGS"],
        [FCFLAGS="$subdir_FCFLAGS $ICON_SCT_FCFLAGS"],
        [MPICC=], [MPIFC=], [--enable-static], [--disable-shared],
        [--disable-maintainer-mode], [--enable-mpi="$enable_mpi"],
        [--enable-openmp="$enable_openmp"], [--disable-rpath],
        [--disable-papi], [--enable-hdf5],
        [--enable-silent-rules="$enable_silent_rules"],
        [--without-libpapi-prefix], [--without-libhdf5-prefix],
        [MPI_LAUNCH="$MPI_LAUNCH"])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/sct/src "])
      BUNDLED_LIBFILES="externals/sct/src/.libs/libsct.a $BUNDLED_LIBFILES"])])

AS_VAR_IF([with_external_yaxt], [yes],
  [AS_IF([test x"$enable_yaxt" = xyes || test x"$enable_cdi_pio" = xyes],
     [AC_CACHE_CHECK([whether the Fortran interface of YAXT library is dnl
available], [icon_cv_have_yaxt],
        [icon_cv_have_yaxt=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use yaxt, only: xt_initialized, xt_initialize
      implicit none
      integer :: a
      if (.not. xt_initialized()) call xt_initialize(a)]])],
           [icon_cv_have_yaxt=yes])])
      AS_VAR_IF([icon_cv_have_yaxt], [no],
        [AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of YAXT library: either disable the YAXT data exchange dnl
(--disable-yaxt) and usage of the parallel features of CDI dnl
(--disable-cdi-pio) or use the bundled version of YAXT library dnl
(--without-external-yaxt), or set FCFLAGS, LDFLAGS and LIBS accordingly])])])],
dnl
  [icon_bundled_yaxt_fc=no
   icon_bundled_yaxt_c=no
   AS_IF(
     [test x"$enable_yaxt" = xyes || test x"$enable_cdi_pio" = xyes],
     [icon_bundled_yaxt_fc=direct],
     [test x"$enable_coupling" = xyes && test x"$with_external_yac" = xno],
     [icon_bundled_yaxt_c=transitive])
dnl Ensure consistency:
   AS_CASE([$icon_bundled_yaxt_fc],
dnl Promote icon_bundled_yaxt_c (if it's not higher yet) to transitive if we
dnl need to link to the Fortran libraries (this is correct because we link
dnl statically):
     [direct|transitive],
     [AS_IF([test x"$icon_bundled_yaxt_c" != xdirect],
        [icon_bundled_yaxt_c=transitive])],
dnl Promote icon_bundled_yaxt_c (if it's not higher yet) to dummy if we need the
dnl Fortran interface:
     [dummy],
     [AS_VAR_IF([icon_bundled_yaxt_c], [no],
        [icon_bundled_yaxt_c=dummy])])
   AS_IF([test x"$icon_bundled_yaxt_c" != xno],
     [icon_subdir_srcdir="$srcdir/externals/yaxt"
      icon_subdir_unique_file="src/xt_idxlist.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
YAXT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the YAXT data exchange (--disable-yaxt), usage of the parallel dnl
features of CDI (--disable-cdi-pio) and usage of the bundled YAC library dnl
(--with-external-yac or --disable-coupling) or provide an external version dnl
of YAXT library (--with-external-yaxt) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      AC_ARG_VAR([ICON_YAXT_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of YAXT library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_YAXT_FCFLAGS=${ICON_YAXT_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_YAXT_CFLAGS],
        [C compiler flags to be appended to CFLAGS when configuring the bundled
version of YAXT library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_YAXT_CFLAGS=${ICON_YAXT_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/yaxt], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/yaxt],
        [BUILD_*=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|MPI_LAUNCH=*],
        [MPIROOT=*|MPI_C_INCLUDE=*|MPI_C_LIB=*|MPI_FC_MOD=*|MPI_FC_LIB=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([shared])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([static])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-root])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-include])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-lib])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-fc-mod])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-fc-lib])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([regard-for-quality])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([on-demand-check-programs])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([example-programs])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([perf-programs])])
      ACX_SUBDIR_APPEND_ARGS([externals/yaxt],
        [CFLAGS="$subdir_CFLAGS $ICON_YAXT_CFLAGS"],
        [FCFLAGS="$subdir_FCFLAGS $ICON_YAXT_FCFLAGS"],
        [MPI_LAUNCH="$MPI_LAUNCH"], [MPIROOT=], [MPI_C_INCLUDE=], [MPI_C_LIB=],
        [MPI_FC_MOD=], [MPI_FC_LIB=], [BUILD_CFLAGS=], [BUILD_FCFLAGS=],
        [BUILD_LDFLAGS=], [BUILD_LIBS=],
        [--enable-silent-rules="$enable_silent_rules"],
        [--disable-maintainer-mode], [--disable-shared], [--enable-static],
        [--with-regard-for-quality="$enable_mpi_checks"],
        [--with-on-demand-check-programs], [--without-example-programs],
        [--without-perf-programs])
dnl Currently, we cannot properly override the environment variables BUILD_CC
dnl and BUILD_FC on the command line when configuring YAXT because the variables
dnl are used even if they are set to empty strings. We could set the variables
dnl on the command line to "$CC" and "$FC", respectively, but that would discard
dnl any compiler flags that the configure script of YAXT might append to CC and
dnl FC instead of CFLAGS and FCFLAGS. As a temporary solution we unset BUILD_CC
dnl and BUILD_FC here (and hope that they will not be unexpectedly set anywhere
dnl below) since it should not affect any other bundled libraries.
      AS_UNSET([BUILD_CC]); AS_UNSET([BUILD_FC])
      AS_CASE([$icon_bundled_yaxt_c],
        [direct],
        [dnl TODO: extend CPPFLAGS when this becomes a real case
         AS_VAR_SET([BUNDLED_LIBFILES],
           ["externals/yaxt/src/.libs/libyaxt_c.a $BUNDLED_LIBFILES"])],
        [transitive],
        [AS_VAR_SET([BUNDLED_LIBFILES],
           ["externals/yaxt/src/.libs/libyaxt_c.a $BUNDLED_LIBFILES"])])
      AS_CASE([$icon_bundled_yaxt_fc],
        [direct],
        [AS_VAR_APPEND([bundled_includes],
           ["${FC_MOD_IN}externals/yaxt/inst_headers/f90 "])
         AS_VAR_SET([BUNDLED_LIBFILES],
           ["externals/yaxt/src/.libs/libyaxt.a $BUNDLED_LIBFILES"])],
        [transitive],
        [AS_VAR_SET([BUNDLED_LIBFILES],
           ["externals/yaxt/src/.libs/libyaxt.a $BUNDLED_LIBFILES"])],
        [no],
        [ACX_SUBDIR_APPEND_ARGS([externals/yaxt], [FC=no])])])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_CACHE_CHECK([whether the Fortran 90 interfaces of ECCODES library is dnl
available], [icon_cv_have_eccodes],
     [icon_cv_have_eccodes=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use eccodes, only: codes_open_file, GRIB_SUCCESS
      implicit none
      integer :: igrib, ierr
      call codes_open_file(igrib, 'conftest.grb', 'w', ierr)
      if (ierr  /= GRIB_SUCCESS) igrib = 0]])],
        [icon_cv_have_eccodes=yes])])
   AS_VAR_IF([icon_cv_have_eccodes], [yes],
     [AS_VAR_APPEND([FCFLAGS], [" ${FC_PP_DEF}GRIBAPI"])],
     [AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interfaces of ECCODES library: set FCFLAGS, LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([with_external_ppm], [yes],
dnl ICON currently does not use PPM directly, therefore we do not check
dnl 'whether the C interface of PPM library is available' when an external
dnl version of PPM library is used. It is the job of the configure script of
dnl CDI.
  [],
dnl
  [icon_bundled_ppm_needed=no
   AS_IF(
     [test x"$enable_cdi_pio" = xyes && test x"$with_external_cdi" = xno],
     [icon_bundled_ppm_needed=transitive_c])
   AS_IF([test x"$icon_bundled_ppm_needed" != xno],
     [icon_subdir_srcdir="$srcdir/externals/ppm"
   icon_subdir_unique_file="src/core/core.c"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
PPM library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable usage of the parallel features of CDI (--disable-cdi-pio) or provide dnl
an external version of PPM library (--with-external-ppm), or an external dnl
version of CDI library (--with-external-cdi) by setting FCFLAGS, LDFLAGS and dnl
LIBS accordingly])])
   AC_ARG_VAR([ICON_PPM_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of PPM library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_PPM_FCFLAGS=${ICON_PPM_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_PPM_CFLAGS],
     [C compiler flags to be appended to CFLAGS when configuring the bundled
version of PPM library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_PPM_CFLAGS=${ICON_PPM_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_SUBDIR_INIT_CONFIG([externals/ppm], [no-run])
   ACX_SUBDIR_REMOVE_ARGS([externals/ppm],
     [BUILD_*=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|MPI_LAUNCH=*],
     [PLATFORM_OPT_LIBM=*|IEEE_ARITHMETICROOT=*|IEEE_ARITHMETIC_FC_MOD=*],
     [IEEE_ARITHMETIC_FC_LIB=*|IEEE_ARHMETIC_FCFLAGS=*|MPIROOT=*],
     [MPI_FC_MOD=*|MPI_FC_LIB=*|MPI_FC_INCLUDE=*|MPI_C_INCLUDE=*|MPI_C_LIB=*],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([mh-file])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([shared])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([static])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([MPI])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([netcdf])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([hdf5])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([parmetis])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([metis])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([crypto])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([setup])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([ieee-arithmetic-root])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([ieee-arithmetic-fc-mod])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([ieee-arithmetic-fc-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-root])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-fc-mod])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-fc-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-include])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([regard-for-quality])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([on-demand-check-programs])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([example-programs])])
   ACX_SUBDIR_APPEND_ARGS([externals/ppm],
     [CFLAGS="$subdir_CFLAGS $ICON_PPM_CFLAGS"],
     [FCFLAGS="$subdir_FCFLAGS $ICON_PPM_FCFLAGS"],
     [MPI_LAUNCH="$MPI_LAUNCH"], [PLATFORM_OPT_LIBM=], [IEEE_ARITHMETICROOT=],
     [IEEE_ARITHMETIC_FC_MOD=], [IEEE_ARITHMETIC_FC_LIB=],
     [IEEE_ARHMETIC_FCFLAGS=], [MPIROOT=], [MPI_C_INCLUDE=], [MPI_C_LIB=],
     [MPI_FC_MOD=], [MPI_FC_LIB=], [BUILD_CFLAGS=], [BUILD_FCFLAGS=],
     [BUILD_LDFLAGS=], [BUILD_FCLDFLAGS=],
     [--enable-silent-rules="$enable_silent_rules"],
     [--disable-maintainer-mode], [--disable-shared], [--enable-static],
     [--enable-MPI], [--disable-netcdf], [--disable-hdf5],
     [--disable-parmetis], [--disable-metis], [--disable-crypto],
     [--with-regard-for-quality="$enable_mpi_checks"],
     [--with-on-demand-check-programs], [--without-example-programs],
dnl There is currently no need in the Fortran interface of the library:
     [FC=no])
   BUNDLED_LIBFILES="externals/ppm/src/.libs/libscalesppmcore.a dnl
$BUNDLED_LIBFILES"])])

AS_IF([test "x$enable_mpi" = xno && test "x$enable_cdi_pio" = xyes],
  [dnl
dnl Fail the configuration because we get the following at the linking stage:
dnl lib/libicon.a(mo_name_list_output.o): In function `data_write_cdipio':
dnl src/io/shared/mo_name_list_output.f90:2445: undefined reference to `var2buf_'
dnl src/io/shared/mo_name_list_output.f90:2451: undefined reference to `var2buf_'
  AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with the parallel features of CDI but MPI dnl
(parallelization) support is disabled: either disable usage of the parallel dnl
features of CDI (--disable-cdi-pio) or enable MPI (parallelization) support dnl
(--enable-mpi)])])

AS_VAR_IF([with_external_cdi], [yes],
dnl Check whether an external version of CDI is provided:
  [AC_CACHE_CHECK([whether the Fortran 90 interface of CDI 1.8.x library is dnl
available], [icon_cv_have_cdi],
     [icon_cv_have_cdi=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: TSTEP_CONSTANT, cdiInqNatts, &
                        vlistInqVarTypeOfGeneratingProcess
      implicit none
      integer :: a, b, c, d
! Appeared in 1.8.x
      a = cdiInqNatts(b, c, d)
! Disappeared in 1.9.x
      a = vlistInqVarTypeOfGeneratingProcess(b, c)]])],
        [icon_cv_have_cdi=yes])])
   AS_VAR_IF([icon_cv_have_cdi], [no],
     [AC_MSG_FAILURE([unable to link a test program using the Fortran 90 dnl
interface of CDI 1.8.x library: either use the bundled version of CDI dnl
library (--without-external-cdi) or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
dnl Check whether the provided external version is compatible with ICON:
   AC_CACHE_CHECK([whether the provided version of CDI 1.8.x is patched],
     [icon_cv_have_cdi_patched],
     [icon_cv_have_cdi_patched=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: t_CdiIterator, cdiIterator_inqFiletype
      implicit none
      integer :: a
      type(t_CdiIterator) :: b
      a = cdiIterator_inqFiletype(b)]])],
        [icon_cv_have_cdi_patched=yes])])
   AS_VAR_IF([icon_cv_have_cdi_patched], [no],
     [AC_MSG_FAILURE([the provided version of CDI 1.8.x does not have the dnl
required patches: use the bundled version of CDI library dnl
(--without-external-cdi)])])
dnl Check whether the provided external version supports the required file
dnl formats:
   icon_cdi_formats="NetCDF-2:CDI_FILETYPE_NC2 NetCDF-4:CDI_FILETYPE_NC4 dnl
GRIB1:CDI_FILETYPE_GRB"
   AS_VAR_IF([enable_grib2], [yes],
     [AS_VAR_APPEND([icon_cdi_formats], [" GRIB2:CDI_FILETYPE_GRB2"])])
   for icon_cdi_format in $icon_cdi_formats; do
     icon_cdi_format_name=`echo "$icon_cdi_format" | cut -d: -f1`
     m4_pushdef([icon_cache_var],
       [icon_cv_have_cdi_[]AS_TR_SH([$icon_cdi_format_name])])dnl
     AC_CACHE_CHECK([whether the provided version of CDI is built with dnl
$icon_cdi_format_name support], [icon_cache_var],
       [AS_VAR_SET([icon_cache_var], [unknown])
        icon_cdi_format_ID=`echo "$icon_cdi_format" | cut -d: -f2`
        AS_VAR_IF([cross_compiling], [no],
          [AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_cdi, only: cdiHaveFiletype, $icon_cdi_format_ID
      write(*, "(i0)") cdiHaveFiletype($icon_cdi_format_ID)]])],
           [icon_exec_result=`./conftest$ac_exeext 2>/dev/null`
            AS_IF([test $? -eq 0],
              [AS_IF([test $icon_exec_result -eq 1],
                 [AS_VAR_SET([icon_cache_var], [yes])],
                 [AS_VAR_SET([icon_cache_var], [no])])])])])])
     AS_IF(
       [test x"AS_VAR_GET(icon_cache_var)" = xno],
       [AC_MSG_FAILURE([the provided version of CDI does not support dnl
$icon_cdi_format_name: use the bundled version of CDI library dnl
(--without-external-cdi)])],
       [test x"AS_VAR_GET(icon_cache_var)" = xunknown],
       [AC_MSG_WARN([unable to check whether the provided version of CDI dnl
supports $icon_cdi_format_name])])
     m4_popdef([icon_cache_var])dnl
   done
dnl Check whether the provided external version supports the parallel features:
   AS_VAR_IF([enable_cdi_pio], [yes],
     [AS_VAR_IF([with_external_yaxt], [no],
        [AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with the parallel features of an external CDI dnl
library and the bundled YAXT library: either disable usage of the parallel dnl
features of CDI (--disable-cdi-pio) or use the bundled version of CDI dnl
library (--without-external-cdi), or provide the version of YAXT library dnl
(--with-external-yaxt) that the provided CDI library has been built with by dnl
setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      AC_CACHE_CHECK([whether the Fortran 77 parallel I/O interface of CDI dnl
library is available], [icon_cv_have_cdi_parallel],
        [icon_cv_have_cdi_parallel=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use yaxt, only: xt_idxlist
      implicit none
      include 'cdipio.inc'
      integer :: cdiCellGridID, gridtype, n_glb, &
               & max_cell_connectivity, grid_deco_part(2)
      type(xt_idxlist), allocatable :: reorder_idxlst_xt(:)
      type(xt_idxlist) :: null_idxlist
      cdiCellGridID = &
      &  cdiPioDistGridCreate(gridtype, n_glb, -1, -1, &
                            & max_cell_connectivity, &
                            & grid_deco_part, &
                            & reorder_idxlst_xt(1), &
                            & null_idxlist, null_idxlist)]])],
           [icon_cv_have_cdi_parallel=yes])])
      AS_VAR_IF([icon_cv_have_cdi_parallel], [no],
        [AC_MSG_FAILURE([unable to link a test program using the Fortran 77 dnl
parallel I/O interface of CDI library: either disable usage of the parallel dnl
features of CDI (--disable-cdi-pio) or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])])],
dnl Configure the bundled version of CDI:
  [icon_subdir_srcdir="$srcdir/externals/cdi"
   icon_subdir_unique_file="src/cdi.h"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of CDI dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of CDI library (--with-external-cdi) by setting FCFLAGS, dnl
LDFLAGS and LIBS accordingly])])
   AC_ARG_VAR([ICON_CDI_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of CDI library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_CDI_FCFLAGS=${ICON_CDI_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_CDI_CFLAGS],
     [C compiler flags to be appended to CFLAGS when configuring the bundled
version of CDI library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_CDI_CFLAGS=${ICON_CDI_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_SUBDIR_INIT_CONFIG([externals/cdi], [no-run])
   ACX_SUBDIR_REMOVE_ARGS([externals/cdi],
     [BUILD_*=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|F77=*|CXX=*],
     [UTIL_LINUX_UUID_C_INCLUDE=*|UTIL_LINUX_UUID_C_LIB=*],
     [OSSP_UUID_C_INCLUDE=*|OSSP_UUID_C_LIB=*|DCE_UUID_C_INCLUDE=*],
     [DCE_UUID_C_LIB=*|MPIROOT=*|MPI_C_INCLUDE=*|MPI_C_LIB=*|MPI_LAUNCH=*],
     [MPI_FC_MOD=*|MPI_FC_LIB=*|PKG_CONFIG=*|YAXT_C_INCLUDE=*|YAXT_C_LIB=*],
     [YAXT_FC_MOD=*|YAXT_FC_LIB=*|PPM_CORE_C_INCLUDE=*|PPM_CORE_C_LIB=*],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([mpi])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([iso-c-interface])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([cf-interface])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([ruby-interface])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([python-interface])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([openmp])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([shared])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([static])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([grib])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([across])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([cgribex])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([cdi-app])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([ppm-dist-array])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([util-linux-uuid-include])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([util-linux-uuid-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([ossp-uuid-include])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([ossp-uuid-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([dce-uuid-include])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([dce-uuid-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([threads])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([fdb5])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([szlib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([netcdf])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([eccodes])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([grib_api])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-root])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-include])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-fc-mod])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([mpi-fc-lib])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([on-demand-check-programs])],
     [ACX_SUBDIR_CONFIG_PATTERN_WITH([example-programs])])
   ACX_SUBDIR_APPEND_ARGS([externals/cdi],
     [CFLAGS="$subdir_CFLAGS $ICON_CDI_CFLAGS"],
     [FCFLAGS="$subdir_FCFLAGS $ICON_CDI_FCFLAGS"],
dnl No need to check for compatibility with F77.
     [F77=no],
dnl Avoid possible inconsistency in the generation of the 'postdeps' variable
dnl for Fortran (!!!) compiler inside the *.la files. We care about those files
dnl because they are used when linking the test suite of CDI (i.e. when running
dnl 'make check'). If you want to know more, take a look at the macro
dnl _LT_LANG_FC_CONFIG in libtool.m4 and note that it does not initialize the
dnl variable 'output_verbose_link_cmd' as it should but uses the value that is
dnl set for the C++ compiler. Therefore, by disabling CXX we get the same
dnl behaviour (still questionable but satisfactory for us) as if CDI didn't use
dnl the C++ compiler at all (e.g. like in YAXT).
     [CXX=no],
     [UTIL_LINUX_UUID_C_INCLUDE=], [UTIL_LINUX_UUID_C_LIB=],
     [OSSP_UUID_C_INCLUDE=], [OSSP_UUID_C_LIB=], [DCE_UUID_C_INCLUDE=],
     [DCE_UUID_C_LIB=], [MPIROOT=], [MPI_C_INCLUDE=], [MPI_C_LIB=],
     [MPI_LAUNCH="$MPI_LAUNCH"], [MPI_FC_MOD=], [MPI_FC_LIB=], [PKG_CONFIG=],
     [YAXT_C_INCLUDE=], [YAXT_C_LIB=], [YAXT_FC_MOD=], [YAXT_FC_LIB=],
     [PPM_CORE_C_INCLUDE=], [PPM_CORE_C_LIB=], [BUILD_CFLAGS=],
     [BUILD_FCFLAGS=], [BUILD_LDFLAGS=], [BUILD_LIBS=], [BUILD_MPI_C_LIB=],
     [BUILD_MPI_FC_LIB=], [BUILD_CC=], [BUILD_CXX=], [BUILD_FC=], [BUILD_F77=],
     [--enable-silent-rules="$enable_silent_rules"],
     [--disable-maintainer-mode], [--enable-mpi="$enable_cdi_pio"],
     [--enable-iso-c-interface], [--enable-cf-interface="$enable_cdi_pio"],
     [--disable-ruby-interface], [--disable-python-interface],
     [--disable-openmp], [--disable-shared], [--enable-static],
     [--enable-grib], [--disable-across], [--enable-cgribex],
     [--disable-cdi-app], [--enable-ppm-dist-array="$enable_cdi_pio"],
     [--without-util-linux-uuid], [--without-ossp-uuid], [--without-dce-uuid],
     [--without-threads], [--without-fdb5], [--without-szlib], [--with-netcdf],
dnl Currently, we don't have full control on choosing the GRIB2 backend (e.g.
dnl grib_api, eccodes) due to the limitations of the configure script of
dnl libcdi. Therefore, we set --with-eccodes=$enable_grib2 --without-grib_api
dnl unconditionally. This way, the configure script of libcdi enables support
dnl for grib2 regardless of whether grib_api or eccodes library is provided
dnl with the LDFLAGS and LIBS variables.
     [--with-eccodes="$enable_grib2"], [--without-grib_api],
     [--with-on-demand-check-programs], [--without-example-programs],
dnl The configure script of CDI runs nc-config to check for the features of the
dnl provided NetCDF library. The script might not be available or might belong
dnl to another NetCDF installation. Therefore, we enforce enabling of NetCDF2
dnl and NetCDF4 features by setting the corresponding cache variables on the
dnl command line.
     [acx_cv_have_netcdf2=yes], [acx_cv_have_netcdf4=yes],
dnl For the same reasons as above, we tell CDI not to use the Parallel-NetCDF
dnl features by default. If the version of the provided NetCDF library is not
dnl very old (i.e. 4.3.3 or later) and the library is configured with the
dnl Parallel-NetCDF features enabled, CDI will override the default value.
     [acx_cv_have_pnetcdf=no],
dnl Disable unwanted extra dependencies.
     [acx_cv_have_libnc_dap=no], [acx_cv_have_nc4hdf5=no])
   AS_VAR_IF([enable_cdi_pio], [yes],
dnl To enable the test suite of CDI, we need to provide the configure script of
dnl the library with information on how to link to the bundled YAXT and PPM
dnl (when they are used). The general approach that we currently follow is to
dnl avoid building the shared versions of the bundled libraries and link to the
dnl static ones without the -L and -l flags. Instead, to avoid possible
dnl ambiguity, the linking is done by passing paths to static libraries as
dnl positional arguments for the linker. That, however, does not work well with
dnl libtool scripts as they inject all files referred by the positional
dnl arguments into the generated static libraries (archives). For example,
dnl libcdipio.a will contain libyaxt_c.a if we tell the configure script of CDI
dnl that linking to YAXT should be done with
dnl "${icon_abs_top_builddir}/externals/yaxt/src/.libs/libyaxt_c.a". Although
dnl that is usually fine with the GNU linker, the MacOS one might ignore the
dnl whole library file:
dnl   ld: warning: ignoring file libcdipio.a, building for macOS-x86_64 but attempting to link with file built for macOS-x86_64
dnl To circumvent all those problems, we tell the configure of CDI to link the
dnl test programs to the bundled YAXT and PPM using libtool .la files.
     [AS_VAR_IF([with_external_yaxt], [no],
        [ACX_SUBDIR_APPEND_ARGS([externals/cdi],
           [acx_cv_c_include_yaxt_h=dnl
"-I${icon_abs_top_builddir}/externals/yaxt/inst_headers"],
           [acx_cv_option_search_xt_initialized_c=dnl
"${icon_abs_top_builddir}/externals/yaxt/src/libyaxt_c.la"],
           [acx_cv_fortran_mod_yaxt=dnl
"${FC_MOD_IN}${icon_abs_top_builddir}/externals/yaxt/inst_headers/f90"],
           [acx_cv_option_search_xt_initialized_fc=dnl
"${icon_abs_top_builddir}/externals/yaxt/src/libyaxt.la dnl
${icon_abs_top_builddir}/externals/yaxt/src/libyaxt_c.la"])])
      AS_VAR_IF([with_external_ppm], [no],
        [ACX_SUBDIR_APPEND_ARGS([externals/cdi],
           [acx_cv_c_include_ppm_ppm_h=dnl
"-I${icon_abs_top_builddir}/externals/ppm/include"],
           [ac_cv_header_ppm_dist_array_h=yes],
           [acx_cv_option_search_PPM_initialize_c=dnl
"${icon_abs_top_builddir}/externals/ppm/src/libscalesppmcore.la"])])
      AS_VAR_APPEND([bundled_includes],
        ["${FC_INC_FLAG}${srcdir}/externals/cdi/src "])
      BUNDLED_LIBFILES="externals/cdi/src/.libs/libcdipio.a $BUNDLED_LIBFILES"])
   AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/cdi/src "])
   BUNDLED_LIBFILES="externals/cdi/src/.libs/libcdi_f2003.a dnl
externals/cdi/src/.libs/libcdi.a $BUNDLED_LIBFILES"])

AS_IF([test x"$enable_serialization" != xno],
  [AC_CACHE_CHECK(
     [whether the Fortran interface of SERIALBOX2 library is available],
     [icon_cv_have_serialbox2],
     [icon_cv_have_serialbox2=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use utils_ppser, only: ppser_set_mode
      implicit none
      call ppser_set_mode(1)]])],
        [icon_cv_have_serialbox2=yes])])
   AS_VAR_IF([icon_cv_have_serialbox2], [no],
     [AC_MSG_FAILURE([unable to link a test program using the Fortrain dnl
interface of SERIALBOX2 library: either disable the serialization dnl
(--disable-serialization) or set FCFLAGS, LDFLAGS and LIBS accordingly])])])

AS_VAR_IF([with_external_mtime], [yes],
  [AC_CACHE_CHECK([whether the Fortran interface of MTIME library is dnl
available], [icon_cv_have_mtime],
     [icon_cv_have_mtime=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mtime, only: max_date_str_len, date, &
                     & newdate, dateToString
      implicit none
      character(len=max_date_str_len) :: test_date_string
      type(date), pointer :: test_date
      test_date => newdate(0,0,0)
      call dateToString(test_date, test_date_string)]])],
        [icon_cv_have_mtime=yes])])
   AS_VAR_IF([icon_cv_have_mtime], [no],
     [AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of MTIME library: either use the bundled version of MTIME dnl
(--without-external-mtime) or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
  [icon_subdir_srcdir="$srcdir/externals/mtime"
   icon_subdir_unique_file="include/mtime_calendar.h"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of MTIME dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': provide an dnl
external version of MTIME library (--with-external-mtime) by setting CFLAGS, dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])
   AC_ARG_VAR([ICON_MTIME_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of MTIME library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_MTIME_FCFLAGS=${ICON_MTIME_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_MTIME_CFLAGS],
     [C compiler flags to be appended to CFLAGS when configuring the bundled
version of MTIME library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_MTIME_CFLAGS=${ICON_MTIME_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_SUBDIR_INIT_CONFIG([externals/mtime], [no-run])
   ACX_SUBDIR_REMOVE_ARGS([externals/mtime],
     [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([static])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([shared])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([examples])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([check])],
     [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([openmp])])
   ACX_SUBDIR_APPEND_ARGS([externals/mtime],
     [CFLAGS="$subdir_CFLAGS $ICON_MTIME_CFLAGS"],
     [FCFLAGS="$subdir_FCFLAGS $ICON_MTIME_FCFLAGS"],
     [--enable-static], [--disable-shared], [--disable-maintainer-mode],
     [--enable-silent-rules="$enable_silent_rules"], [--disable-examples],
dnl Disable the broken unit test:
     [--disable-check],
     [--enable-openmp="$enable_openmp"])
   AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/mtime/src "])
   BUNDLED_LIBFILES="externals/mtime/src/.libs/libmtime.a $BUNDLED_LIBFILES"])

AS_VAR_IF([enable_coupling], [yes],
  [AS_VAR_IF([with_external_yac], [yes],
     [AS_VAR_IF([with_external_mtime], [no],
        [AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with an external YAC library and the bundled dnl
MTIME library: either disable the coupling (--disable-coupling) or use the dnl
bundled version of YAC library (--without-external-yac), or provide the dnl
version of MTIME library (--with-external-mtime) that the provided YAC dnl
library has been built with by setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      AS_VAR_IF([with_external_yaxt], [no],
        [AC_MSG_ERROR([unsupported combination of the configure options: dnl
you are trying to configure with an external YAC library and the bundled dnl
YAXT library: either disable the coupling (--disable-coupling) or use the dnl
bundled version of YAC library (--without-external-yac), or provide the dnl
version of YAXT library (--with-external-yaxt) that the provided YAC dnl
library has been built with by setting FCFLAGS, LDFLAGS and LIBS accordingly])])
      AC_CACHE_CHECK([whether the Fortran interface of YAC 2.x library is dnl
available], [icon_cv_have_yac],
        [icon_cv_have_yac=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_yac_finterface, only: yac_fdef_grid
      implicit none
      integer, parameter :: v = 3, c = 1, vpc = 3
      integer :: c2v(vpc, c), id
      real :: x(v), y(v)
      call yac_fdef_grid('', v, c, vpc, x, y, c2v, id)]])],
           [icon_cv_have_yac=yes])])
      AS_VAR_IF([icon_cv_have_yac], [no],
        [AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of YAC 2.x library: either disable the coupling dnl
(--disable-coupling) or use the bundled version of YAC library dnl
(--without-external-yac), or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/yac"
      icon_subdir_unique_file="src/grid.c"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of YAC dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the coupling (--disable-coupling) or provide an external version of dnl
YAC library (--with-external-yac) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
      AC_ARG_VAR([ICON_YAC_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of YAC library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_YAC_FCFLAGS=${ICON_YAC_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_YAC_CFLAGS],
        [C compiler flags to be appended to CFLAGS when configuring the bundled
version of YAC library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_YAC_CFLAGS=${ICON_YAC_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/yac], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/yac],
        [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*|MPI_LAUNCH=*],
        [MPI_CFLAGS=*|MPI_CLIBS=*|MPI_FCFLAGS=*|MPI_FCLIBS=*],
        [YAXT_CFLAGS=*|YAXT_CLIBS=*|YAXT_FCFLAGS=*|YAXT_FCLIBS=*],
        [NETCDF_CFLAGS=*|NETCDF_CLIBS=*|FYAML_CFLAGS=*|FYAML_CLIBS=*],
        [MTIME_CFLAGS=*|MTIME_CLIBS=*|FORTRAN_LAPACK_CLIBS=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([lib-only])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([netcdf])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([external-lapack])],
        [ACX_SUBDIR_CONFIG_PATTERN_WITH([external-mtime])])
      ACX_SUBDIR_APPEND_ARGS([externals/yac],
        [CFLAGS="$subdir_CFLAGS $ICON_YAC_CFLAGS"],
        [FCFLAGS="$subdir_FCFLAGS $ICON_YAC_FCFLAGS"],
        [MPI_CFLAGS=], [MPI_CLIBS=], [MPI_FCFLAGS=], [MPI_FCLIBS=],
        [YAXT_CLIBS=], [YAXT_FCLIBS=], [NETCDF_CFLAGS=], [NETCDF_CLIBS=],
        [FYAML_CFLAGS=], [FYAML_CLIBS=], [MTIME_CLIBS=],
        [FORTRAN_LAPACK_CLIBS=], [--enable-silent-rules="$enable_silent_rules"],
        [--disable-maintainer-mode], [--enable-lib-only], [--enable-netcdf],
        [--with-external-lapack=fortran], [--with-external-mtime])
      AS_VAR_SET_IF([MPI_LAUNCH],
        [ACX_SUBDIR_APPEND_ARGS([externals/yac], [MPI_LAUNCH="$MPI_LAUNCH"])])
      AS_VAR_IF([with_external_yaxt], [yes],
        [ACX_SUBDIR_APPEND_ARGS([externals/yac],
           [YAXT_CFLAGS=], [YAXT_FCFLAGS=])],
        [ACX_SUBDIR_RUN_EXTRA([externals/yac],
           [icon_try=dnl
"${MAKE-make} -C '${icon_abs_top_builddir}/externals/yaxt/inst_headers' yaxt.h"
            AC_MSG_NOTICE([running $icon_try])
            AC_TRY_EVAL([icon_try]) >&AS_MESSAGE_LOG_FD])
         ACX_SUBDIR_APPEND_ARGS([externals/yac],
           [YAXT_CFLAGS=dnl
"-I${icon_abs_top_builddir}/externals/yaxt/inst_headers"],
           [acx_cv_c_lib_func_xt_initialized=dnl
"${icon_abs_top_builddir}/externals/yaxt/src/.libs/libyaxt_c.a"])
         AS_IF([test x"${icon_bundled_yaxt_fc-no}" != xno],
           [ACX_SUBDIR_APPEND_ARGS([externals/yac],
              [YAXT_FCFLAGS=dnl
"${FC_MOD_IN}${icon_abs_top_builddir}/externals/yaxt/inst_headers/f90"],
              [acx_cv_fc_module_YAXT=yes],
              [acx_cv_fc_lib_func_XT_INITIALIZED=dnl
"${icon_abs_top_builddir}/externals/yaxt/src/.libs/libyaxt.a dnl
${icon_abs_top_builddir}/externals/yaxt/src/.libs/libyaxt_c.a"])],
           [ACX_SUBDIR_APPEND_ARGS([externals/yac], [YAXT_FCFLAGS=])])])
      AS_VAR_IF([with_external_mtime], [yes],
        [ACX_SUBDIR_APPEND_ARGS([externals/yac], [MTIME_CFLAGS=])],
        [ACX_SUBDIR_APPEND_ARGS([externals/yac],
           [MTIME_CFLAGS="-I${icon_abs_top_srcdir}/externals/mtime/include"],
           [acx_cv_c_lib_func_datetimeToString=dnl
"${icon_abs_top_builddir}/externals/mtime/src/.libs/libmtime.a"])])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/yac/src "])
      BUNDLED_LIBFILES="externals/yac/src/libyac.a $BUNDLED_LIBFILES"])])

AS_VAR_IF([enable_art], [yes],
  [AS_VAR_IF([with_external_tixi], [yes],
     [AC_CACHE_CHECK([whether the Fortran interface of TIXI library is dnl
available], [icon_cv_have_tixi],
        [icon_cv_have_tixi=no
         AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use tixi, only: tixiOpenDocument
      implicit none
      character(len=200) :: a
      integer :: b
      b = tixiOpenDocument(a, b)]])],
           [icon_cv_have_tixi=yes])])
      AS_VAR_IF([icon_cv_have_tixi], [no],
        [AC_MSG_FAILURE([unable to link a test program using the Fortran dnl
interface of TIXI library: either disable the aerosols and reactive trace dnl
component (--disable-art) or use the bundled version of TIXI library dnl
(--without-external-tixi), or set FCFLAGS, LDFLAGS and LIBS accordingly])])],
dnl
     [icon_subdir_srcdir="$srcdir/externals/tixi"
      icon_subdir_unique_file="src/tixi.f90"
      AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
        [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
TIXI library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': either dnl
disable the aerosols and reactive trace component (--disable-art) or provide dnl
an external version of TIXI library (--with-external-tixi) by setting dnl
FCFLAGS, LDFLAGS and LIBS accordingly])])
      AC_ARG_VAR([ICON_TIXI_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of TIXI library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_TIXI_FCFLAGS=${ICON_TIXI_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      AC_ARG_VAR([ICON_TIXI_CFLAGS],
        [C compiler flags to be appended to CFLAGS when configuring the bundled
version of TIXI library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
      ICON_TIXI_CFLAGS=${ICON_TIXI_CFLAGS-$ICON_BUNDLED_CFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/tixi], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/tixi],
        [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
        [XML2_CFLAGS=*|XML2_CLIBS=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([maintainer-mode])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])])
      ACX_SUBDIR_APPEND_ARGS([externals/tixi],
        [CFLAGS="$subdir_CFLAGS $ICON_TIXI_CFLAGS"],
        [FCFLAGS="$subdir_FCFLAGS $ICON_TIXI_FCFLAGS"],
        [XML2_CFLAGS=], [XML2_CLIBS=], [--disable-maintainer-mode],
        [--enable-silent-rules="$enable_silent_rules"])
      AS_VAR_APPEND([bundled_includes], ["${FC_MOD_IN}externals/tixi/src "])
      BUNDLED_LIBFILES="externals/tixi/src/libicon_tixi.a dnl
$BUNDLED_LIBFILES"])])

AS_CASE([$icon_gpu_secondary],
  [cuda],
  [ACX_LANG_PUSH_CUDA
   AC_CHECK_HEADER([cub/device/device_select.cuh], [],
     [AC_MSG_FAILURE([unable to find header files of CUB library: either dnl
switch the secondary GPU compiler (e.g. dnl
--enable-gpu=${icon_gpu_primary}+hip) or disable GPU support
(--disable-gpu), or set CUDAFLAGS accordingly])],
     [[// override default includes]])
   ACX_LANG_POP_CUDA],
  [hip],
  [ACX_LANG_PUSH_HIP
   AC_CHECK_HEADER([hipcub/device/device_select.hpp], [],
     [AC_MSG_FAILURE([unable to find header files of HIPCUB library: either dnl
switch the secondary GPU compiler (e.g. dnl
--enable-gpu=${icon_gpu_primary}+cuda) or disable GPU support dnl
(--disable-gpu), or set HIPFLAGS accordingly])],
     [[// override default includes]])
   ACX_LANG_POP_HIP])

AS_VAR_IF([enable_comin], [yes],
  [icon_subdir_srcdir="$srcdir/externals/comin"
   icon_subdir_unique_file="src/comin_setup.F90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of COMIN dnl
library ($icon_subdir_unique_file) in '$icon_subdir_srcdir'])])
   AC_ARG_VAR([ICON_COMIN_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of COMIN library @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
   ICON_COMIN_FCFLAGS=${ICON_COMIN_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_COMIN_CFLAGS],
     [C compiler flags to be appended to CFLAGS when configuring the bundled
version of COMIN library @<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   AC_ARG_VAR([ICON_COMIN_CXXFLAGS],
     [C++ compiler flags to be appended to CXXFLAGS when configuring the bundled
version of COMIN library @<:@default=$ICON_BUNDLED_CXXFLAGS@:>@])dnl
   ICON_COMIN_CFLAGS=${ICON_COMIN_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_SUBDIR_INIT_CMAKE([externals/comin], [no-run no-adjust-compilers], [],
     [AS_ESCAPE(["${CMAKE}"])])
   ACX_SUBDIR_REMOVE_ARGS([externals/comin],
dnl The compiler flags are replaced below:
     [-DCMAKE_C_FLAGS=*|-DCMAKE_Fortran_FLAGS=*|-DCMAKE_CXX_FLAGS=*],
dnl The interprocedural optimization is not enabled in the project:
     [-DCMAKE_*_COMPILER_AR=*|-DCMAKE_*_COMPILER_RANLIB=*])
   ACX_SUBDIR_APPEND_ARGS([externals/comin],
dnl Extend the compiler flags:
     [-DCMAKE_C_FLAGS="$subdir_CFLAGS $ICON_COMIN_CFLAGS $subdir_CPPFLAGS"],
     [-DCMAKE_Fortran_FLAGS="$subdir_FCFLAGS $ICON_COMIN_FCFLAGS"],
     [-DCMAKE_CXX_FLAGS=dnl
"$subdir_CXXFLAGS $ICON_COMIN_CXXFLAGS $subdir_CPPFLAGS"],
dnl Build the static version of the library:
     [-DBUILD_SHARED_LIBS:BOOL=OFF],
dnl Turn off the examples (only build the static lib)
     [-DCOMIN_ENABLE_EXAMPLES:BOOL=OFF],
dnl Set non-existing build type to avoid extra compiler flags:
     [-DCMAKE_BUILD_TYPE:STRING=NOEXTRAFLAGS])
   AS_VAR_APPEND([bundled_includes],
     ["${FC_MOD_IN}externals/comin/build/mod "])
   BUNDLED_LIBFILES="externals/comin/build/libcomin.a $BUNDLED_LIBFILES"])

AS_VAR_IF([with_external_fortran_support], [yes],
  [AC_CACHE_CHECK([whether Fortran support library is available],
     [icon_cv_have_fortran_support],
     [icon_cv_have_fortran_support=no
      AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[
      use mo_util_backtrace, only: ftn_util_backtrace
      implicit none
      call ftn_util_backtrace()]])],
        [icon_cv_have_fortran_support=yes])])
   AS_VAR_IF([icon_cv_have_fortran_support], [no],
     [AC_MSG_FAILURE([unable to link a test program using Fortran support dnl
library: either use the bundled version of FORTRAN-SUPPORT dnl
(--without-external-fortran-support) or set FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])],
dnl
  [icon_subdir_srcdir="$srcdir/externals/fortran-support"
   icon_subdir_unique_file="src/mo_util_backtrace.F90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
FORTRAN-SUPPORT library ($icon_subdir_unique_file) in '$icon_subdir_srcdir': dnl
provide an external version of FORTRAN-SUPPORT library dnl
(--with-external-fortran-support) by setting FCFLAGS, LDFLAGS and LIBS dnl
accordingly])])
   AC_ARG_VAR([ICON_FORTRAN_SUPPORT_FCFLAGS],
     [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of FORTRAN-SUPPORT library
@<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
ICON_FORTRAN_SUPPORT_FCFLAGS=dnl
${ICON_FORTRAN_SUPPORT_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
   AC_ARG_VAR([ICON_FORTRAN_SUPPORT_CFLAGS],
     [Fortran compiler flags to be appended to CFLAGS when configuring the
bundled version of FORTRAN-SUPPORT library
@<:@default=$ICON_BUNDLED_CFLAGS@:>@])dnl
   ICON_FORTRAN_SUPPORT_CFLAGS=dnl
${ICON_FORTRAN_SUPPORT_CFLAGS-$ICON_BUNDLED_CFLAGS}
   ACX_SUBDIR_INIT_CMAKE([externals/fortran-support],
     [no-run no-adjust-compilers], [],
     [AS_ESCAPE(["${CMAKE}"])])
   ACX_SUBDIR_REMOVE_ARGS([externals/fortran-support],
dnl C and Fortran compiler flags are replaced below:
     [-DCMAKE_C_FLAGS=*|-DCMAKE_Fortran_FLAGS=*],
dnl C++ compiler is used for tests only, which we currently disable:
     [-DCMAKE_CXX_*=*],
dnl We build static libraries only:
     [-DCMAKE_SHARED_LINKER_FLAGS=*|-DCMAKE_MODULE_LINKER_FLAGS=*],
dnl The interprocedural optimization is not enabled in the project:
     [-DCMAKE_*_COMPILER_AR=*|-DCMAKE_*_COMPILER_RANLIB=*])
   ACX_SUBDIR_APPEND_ARGS([externals/fortran-support],
dnl Extend C and Fortran compiler flags:
     [-DCMAKE_C_FLAGS=dnl
"$subdir_CFLAGS $ICON_FORTRAN_SUPPORT_CFLAGS $subdir_CPPFLAGS"],
     [-DCMAKE_Fortran_FLAGS="$subdir_FCFLAGS $ICON_FORTRAN_SUPPORT_FCFLAGS"],
dnl Build the static version of the library and disable the tests:
     [-DBUILD_SHARED_LIBS:BOOL=OFF], [-DBUILD_TESTING:BOOL=OFF],
dnl Set non-existing build type to avoid extra compiler flags:
     [-DCMAKE_BUILD_TYPE:STRING=NOEXTRAFLAGS])
   AS_VAR_APPEND([bundled_includes],
     ["${FC_MOD_IN}externals/fortran-support/build/src/mod "])
   BUNDLED_LIBFILES=dnl
"externals/fortran-support/build/src/libfortran-support.a $BUNDLED_LIBFILES"])

dnl Prepend the include search flags to FCFLAGS:
AS_IF([test -n "$bundled_includes"], [FCFLAGS="${bundled_includes}${FCFLAGS}"])

############################# BUNDLED EXECUTABLES ##############################

AS_IF([test x"$enable_hd" != xno],
  [icon_subdir_srcdir="$srcdir/externals/hd"
   icon_subdir_unique_file="code/src/hd_driver.f90"
   AS_IF([test ! -r "$icon_subdir_srcdir/$icon_subdir_unique_file"],
     [AC_MSG_ERROR([unable to find sources of the bundled version of dnl
HD model ($icon_subdir_unique_file) in '$icon_subdir_srcdir': disable dnl
building of the bundled Hydrological Discharge (HD) model dnl
(--without-bundled-hd)])])
      AC_ARG_VAR([ICON_HD_FCFLAGS],
        [Fortran compiler flags to be appended to FCFLAGS when configuring the
bundled version of HD model @<:@default=$ICON_BUNDLED_FCFLAGS@:>@])dnl
      ICON_HD_FCFLAGS=${ICON_HD_FCFLAGS-$ICON_BUNDLED_FCFLAGS}
      ACX_SUBDIR_INIT_CONFIG([externals/hd], [no-run])
      ACX_SUBDIR_REMOVE_ARGS([externals/hd],
        [BUILD_ENV=*|ICON_*FLAGS=*|CFLAGS=*|FCFLAGS=*],
        [MPI_FCFLAGS=*|MPI_FCLIBS=*|NETCDF_FCFLAGS=*|NETCDF_FCLIBS=*],
        [YAC_FCFLAGS=*|YAC_FCLIBS=*|PKG_CONFIG=*],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([rpaths])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([silent-rules])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([mpi])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([resolution])],
        [ACX_SUBDIR_CONFIG_PATTERN_ENABLE([coupling])])
      ACX_SUBDIR_APPEND_ARGS([externals/hd],
        [FCFLAGS="$subdir_FCFLAGS $ICON_HD_FCFLAGS"],
        [MPI_FCFLAGS=], [MPI_FCLIBS=], [NETCDF_FCFLAGS=], [NETCDF_FCLIBS=],
        [YAC_FCLIBS=], [PKG_CONFIG=], [--enable-rpaths="$enable_rpaths"],
        [--enable-silent-rules="$enable_silent_rules"],
        [--enable-mpi="$enable_mpi"], [--enable-resolution="$enable_hd"],
        [--enable-coupling="$enable_coupling"])
      AS_VAR_IF([enable_coupling], [yes],
        [AS_VAR_IF([with_external_yac], [yes],
           [ACX_SUBDIR_APPEND_ARGS([externals/hd], [YAC_FCFLAGS=])],
           [icon_hd_yac_lib=dnl
"${icon_abs_top_builddir}/externals/yac/src/libyac.a"
            AS_VAR_IF([with_external_mtime], [no],
              [AS_VAR_APPEND([icon_hd_yac_lib],
                 [dnl
" ${icon_abs_top_builddir}/externals/mtime/src/.libs/libmtime.a"])])
            AS_VAR_IF([with_external_yaxt], [no],
              [AS_VAR_APPEND([icon_hd_yac_lib],
                 [dnl
" ${icon_abs_top_builddir}/externals/yaxt/src/.libs/libyaxt_c.a"])])
            ACX_SUBDIR_APPEND_ARGS([externals/hd],
              [YAC_FCFLAGS=dnl
"${FC_MOD_IN}${icon_abs_top_builddir}/externals/yac/src"],
              [acx_cv_fc_module_MO_YAC_FINTERFACE=yes],
              [acx_cv_fc_lib_func_YAC_FGET_VERSION="${icon_hd_yac_lib}"])])])])

######################## BUNDLED PACKAGES FINALIZATION #########################

dnl Generate lists of the build directories of all active Autoconf-based and
dnl CMake-based bundled libraries:
AC_SUBST([icon_bundled_config_subdirs], [''])
AC_SUBST([icon_bundled_cmake_subdirs], [''])
for icon_bundled_srcdir in $extra_src_subdirs; do
  ACX_SUBDIR_GET_BUILD_DIR([icon_bundled_builddir], [$icon_bundled_srcdir])
  ACX_SUBDIR_GET_BUILD_TYPE([icon_bundled_type], [$icon_bundled_srcdir])
  AS_VAR_IF([icon_bundled_${icon_bundled_type}_subdirs], [],
    [AS_VAR_COPY([icon_bundled_${icon_bundled_type}_subdirs],
       [icon_bundled_builddir])],
    [AS_VAR_APPEND([icon_bundled_${icon_bundled_type}_subdirs],
       [" $icon_bundled_builddir"])])
done

dnl Either generate rules for the delayed configuration of the bundled
dnl libraries or register all the libraries for the configuration by this
dnl top-level configure script:
AC_SUBST([ICON_DELAYED_CONFIG_RULES], [''])
AS_VAR_IF([enable_delayed_config], [yes],
  [dnl
dnl Make sure that the bundled libraries pick up the same FC and CC if not set
dnl by the user on the command line when configuring ICON:
   ICON_DELAYED_CONFIG_RULES="export FC CC CXX$as_nl$as_nl"
   icon_rule_sep=''
   for icon_bundled_srcdir in $extra_src_subdirs; do
     ACX_SUBDIR_GET_BUILD_DIR([icon_bundled_builddir], [$icon_bundled_srcdir])
     ACX_SUBDIR_GET_RUN_CMD([icon_bundled_config_cmd], [$icon_bundled_srcdir])
     ACX_SUBDIR_GET_BUILD_TYPE([icon_bundled_type], [$icon_bundled_srcdir])
     AS_VAR_IF([icon_bundled_type], [cmake],
       [icon_bundled_silent_var='silent_CMAKE'
        icon_bundled_silent_suffix='>/dev/null'],
       [icon_bundled_silent_var='silent_CONFIG'
        icon_bundled_silent_suffix='--silent'])
     AS_VAR_APPEND([ICON_DELAYED_CONFIG_RULES], [dnl
"$icon_rule_sep$icon_bundled_builddir: dnl
\$(if \$(filter-out mostlyclean clean distclean,\$(or \$(MAKECMDGOALS),all)),dnl
$icon_bundled_builddir/Makefile)
$icon_bundled_builddir/Makefile: icon.mk
	\$($icon_bundled_silent_var)\$(MKDIR_P) \$(@D) && cd \$(@D) && dnl
$icon_bundled_config_cmd \$(if \$($icon_bundled_silent_var),dnl
$icon_bundled_silent_suffix)"])
     icon_rule_sep="$as_nl$as_nl"
   done],
  [for icon_bundled_srcdir in $extra_src_subdirs; do
     ACX_SUBDIR_RUN_RESET([$icon_bundled_srcdir], [yes])
   done])

##################### PREPROCESSORS AND AUXILIARY COMMANDS #####################

AC_PROG_INSTALL
AC_PROG_MKDIR_P
dnl The message generated by AC_PROG_LN_S is misleading in our case, therefore
dnl we do what the standard macro does but without the message:
AC_SUBST([LN_S], [$as_ln_s])

AC_ARG_VAR([AR], [archiver command])dnl
AC_ARG_VAR([ARFLAGS], [archiver flags])dnl
ARFLAGS=${ARFLAGS-cr}
ACX_PROG_SEARCH([AR], [ar],
  [AC_COMPILE_IFELSE([AC_LANG_PROGRAM],
     [icon_prog_try=dnl
'$acx_candidate ${ARFLAGS} libconftest.a conftest.$ac_objext >/dev/null'
      AC_TRY_EVAL([icon_prog_try])
      rm -f libconftest.a])
   test $ac_status -eq 0],
  [AR=$acx_cv_prog_AR])
AC_ARG_VAR([RANLIB], [archive indexer command])dnl
AC_PROG_RANLIB

AC_ARG_VAR([PYTHON], [Python interpreter command])dnl
ACX_PROG_SEARCH([PYTHON], [python python3],
  [icon_prog_try=dnl
'$acx_candidate $srcdir/utils/mkhelper/depgen.py -h >/dev/null'
   AC_TRY_EVAL([icon_prog_try])],
  [PYTHON=$acx_cv_prog_PYTHON])

AC_ARG_VAR([PERL], [Perl interpreter command])dnl
ACX_PROG_SEARCH([PERL], [perl],
  [icon_prog_try='$acx_candidate $srcdir/utils/pvcs.pl --help >/dev/null'
   AC_TRY_EVAL([icon_prog_try])
   test $ac_status -eq 0 || test $ac_status -eq 1],
  [PERL=$acx_cv_prog_PERL])

AC_ARG_VAR([CMAKE], [CMake command])dnl
AS_IF([test -n "$icon_bundled_cmake_subdirs"],
  [ACX_PROG_SEARCH([CMAKE], [cmake],
     [icon_prog_try='$acx_candidate --version >&AS_MESSAGE_LOG_FD'
      AC_TRY_EVAL([icon_prog_try])],
     [CMAKE=$acx_cv_prog_CMAKE])],
  [CMAKE=])

AC_ARG_VAR([FPP], [Fortran preprocessor command (optional)])dnl
AS_VAR_IF([enable_explicit_fpp], [yes],
  [ACX_PROG_FPP([FPP],
     [FPP=$acx_cv_prog_fpp],
     [AC_MSG_FAILURE([unable to find a valid Fortran preprocessor command: dnl
either disable explicit Fortran preprocessing (--disable-explicit-fpp) or dnl
set FPP accordingly])],
     [$FC_INC_FLAG_PP_f90], [$FC_PP_DEF], [utils/fpp-wrappers])],
  [FPP=])

AC_ARG_VAR([SB2PP], [Serialbox2 preprocessor command (optional)])dnl
AS_IF([test x"$enable_serialization" != xno],
  [AC_LANG_ASSERT([Fortran])
   AC_MSG_CHECKING([for Serialbox2 preprocessor command])
   icon_sb2pp_result=unknown
   SB2PP=${SB2PP-'${PYTHON} -m pp_ser'}
dnl pp_ser.py can not process a piece of code like this:
dnl   program main
dnl   end
dnl Therefore, we need to generate the source code manually:
   AC_LANG_CONFTEST([AC_LANG_SOURCE(
[[      program main
      implicit none
      !\$ser mode write
      call ppser_set_mode(1)
      end program]])])
   icon_prog_try=dnl
"$SB2PP -o conftest.$ac_ext conftest.$ac_ext >&AS_MESSAGE_LOG_FD"
   AC_TRY_EVAL([icon_prog_try])
   AS_IF([test $? -eq 0], [AC_LINK_IFELSE([], [icon_sb2pp_result=$SB2PP])])
   rm -f conftest.$ac_ext
   AC_MSG_RESULT([$icon_sb2pp_result])
   AS_VAR_IF([icon_sb2pp_result], [unknown],
     [AC_MSG_FAILURE([unable to find a valid Serialbox2 preprocessor dnl
command: either disable the Serialbox2 serialization dnl
(--disable-serialization) or make sure that path to the 'pp_ser' module is dnl
in the PYTHONPATH environment variable, or set SB2PP accordingly])])],
   [SB2PP=])

dnl We have checked this command earlier but put the declaration here so
dnl it would be listed in the help message of the configure script at the very
dnl end, together with other optional commands.
AC_ARG_VAR([MPI_LAUNCH], [interactive MPI launcher command (optional)])dnl

############################ MAKEFILE CONDITIONALS #############################

AS_VAR_IF([enable_silent_rules], [yes],
  [AC_SUBST([DEFAULT_VERBOSITY], [0])],
  [AC_SUBST([DEFAULT_VERBOSITY], [1])])

AS_VAR_IF([enable_delayed_config], [yes],
  [AC_SUBST([DELAYED_CONFIG_ENABLED], [''])],
  [AC_SUBST([DELAYED_CONFIG_ENABLED], ['#'])])

AS_VAR_IF([enable_jsbach], [yes],
  [AC_SUBST([JSBACH_ENABLED], [''])],
  [AC_SUBST([JSBACH_ENABLED], ['#'])])

AS_VAR_IF([enable_dace], [yes],
  [AC_SUBST([DACE_ENABLED], [''])],
  [AC_SUBST([DACE_ENABLED], ['#'])])

AS_VAR_IF([enable_emvorado], [yes],
  [AC_SUBST([EMVORADO_ENABLED], [''])],
  [AC_SUBST([EMVORADO_ENABLED], ['#'])])

AS_VAR_IF([enable_art], [yes],
  [AC_SUBST([ART_ENABLED], [''])],
  [AC_SUBST([ART_ENABLED], ['#'])])

AS_VAR_IF([enable_hd], [no],
  [AC_SUBST([HD_ENABLED], ['#'])],
  [AC_SUBST([HD_ENABLED], [''])])

AS_VAR_IF([enable_mpi], [yes],
  [AC_SUBST([MPI_ENABLED], [''])],
  [AC_SUBST([MPI_ENABLED], ['#'])])

AS_VAR_IF([enable_openmp], [yes],
  [AC_SUBST([OPENMP_ENABLED], [''])],
  [AC_SUBST([OPENMP_ENABLED], ['#'])])

AS_IF([test x"$icon_gpu_primary" != xno],
  [AC_SUBST([GPU_ENABLED], [''])],
  [AC_SUBST([GPU_ENABLED], ['#'])])

AS_VAR_IF([icon_gpu_secondary], [cuda],
  [AC_SUBST([CUDA_ENABLED], [''])],
  [AC_SUBST([CUDA_ENABLED], ['#'])])

AS_VAR_IF([icon_gpu_secondary], [hip],
  [AC_SUBST([HIP_ENABLED], [''])],
  [AC_SUBST([HIP_ENABLED], ['#'])])

AS_VAR_IF([enable_explicit_fpp], [yes],
  [AC_SUBST([FPP_ENABLED], [''])],
  [AC_SUBST([FPP_ENABLED], ['#'])])

AS_IF([test x"$enable_serialization" != xno],
  [AC_SUBST([SB2_ENABLED], [''])],
  [AC_SUBST([SB2_ENABLED], ['#'])])

AS_VAR_IF([enable_pgi_inlib], [yes],
  [AC_SUBST([PGI_INLIB_ENABLED], [''])
   AC_CONFIG_FILES([inlib.mk])],
  [AC_SUBST([PGI_INLIB_ENABLED], ['#'])])

#################################### OUTPUT ####################################

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
  Makefile
  depgen.c.config
  depgen.f90.config
  deplist.config
  icon.mk
])

AC_CONFIG_FILES([
  collect.extra-libs
  run/collect.set-up.info
], [chmod +x "$ac_file"])

dnl The make_runscripts generator should call run/collect.set-up.info and
dnl generate run/set-up.info. This is currently not the case, therefore we
dnl delegate this task to the configure script:
AC_CONFIG_COMMANDS([run/set-up.info],
  [run/collect.set-up.info "$ac_file" 2>/dev/null || dnl
   AC_MSG_ERROR([failed to create $ac_file])])

dnl The mkexp generator should not use run/set-up.info. Ideally, it should not
dnl use run/collect.set-up.info either but have its own template file. The file
dnl to be generated by the configure script using that template could still be a
dnl shell or even a Python script but it should be run by mkexp itself, not by
dnl the configure script. This way, both generators would be independent of each
dnl other. Until then, we generate run/SETUP.config as a symlink to
dnl run/set-up.info. Note that we have to use AC_CONFIG_COMMANDS and not
dnl AC_CONFIG_LINKS because the configure script generates all symlinks before
dnl running any of the commands, which in our case results into a failed attempt
dnl to create a symlink (i.e. run/SETUP.config) to a non-existent file (i.e.
dnl run/set-up.info). Therefore, we copy the standard logic of the symlink
dnl generation here:
AC_CONFIG_COMMANDS([run/SETUP.config],
  [rm -f "$ac_file"
   ac_rel_source='set-up.info'
   ac_source="run/${ac_rel_source}"
   ln -s "$ac_rel_source" "$ac_file" 2>/dev/null ||
     ln "$ac_source" "$ac_file" 2>/dev/null ||
     cp -p "$ac_source" "$ac_file" ||
     AC_MSG_ERROR([cannot link or copy $ac_source to $ac_file])])

dnl An additional step to generate dummy dependencies:
AC_CONFIG_COMMANDS([depfiles],
  [${MAKE-make} -f icon.mk dummy-depend >/dev/null 2>&AS_MESSAGE_LOG_FD
   sleep 1
   touch icon.mk])

AC_OUTPUT

AS_IF([test "x$enable_delayed_config" = xyes && test -n "$extra_src_subdirs"],
  [AS_ECHO(["
The following bundled libraries will be configured at the build time:dnl
`AS_ECHO(["$extra_src_subdirs"]) | tr ' ' '
'`
"]) >&AS_MESSAGE_FD])
