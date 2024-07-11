#/bin/bash

# Checks that the input files contain the license header (prefixed with the
# respective language-specific comment character):
license=' ICON

 ---------------------------------------------------------------
 Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
 Contact information: icon-model.org

 See AUTHORS.TXT for a list of authors
 See LICENSES/ for license information
 SPDX-License-Identifier: BSD-3-Clause
 ---------------------------------------------------------------'
# The command line arguments can be paths either to regular files or to
# directories. In the latter case, the directories are searched recursively and
# all files with names that match known patterns are considered as input files
# for checking. If no arguments are provided, the script runs the standard
# Buildbot test, i.e. checks files in ICON source directories.

# Currently, we support only Fortran:
comment='!'

# The known patterns are (newline-separated list of name patterns):
known_patterns='
*.F90
*.f90
*.inc
*.incf
'

# ICON source directories (newline-separated list of paths relative to the root
# source directory):
icon_directories='
src
support
test
'

# ICON ignored patterns (newline-separated list of path patterns relative to the
# root source directory):
icon_ignored_patterns='
src/sea_ice/dynamics_fem/mo_ice_fem_advection.f90
src/sea_ice/dynamics_fem/mo_ice_fem_evp.f90
src/sea_ice/dynamics_fem/mo_ice_fem_init.f90
src/sea_ice/dynamics_fem/mo_ice_fem_mesh.f90
src/diagnostics/atmosphere/mo_diag_hailcast.f90
src/atm_phy_schemes/cloud_random_numbers.f90
src/atm_phy_schemes/data_gwd.f90
src/atm_phy_schemes/mo_adjust.f90
src/atm_phy_schemes/mo_cuascn.f90
src/atm_phy_schemes/mo_cudescn.f90
src/atm_phy_schemes/mo_cuflxtends.f90
src/atm_phy_schemes/mo_cufunctions.f90
src/atm_phy_schemes/mo_cuinit.f90
src/atm_phy_schemes/mo_cumaster.f90
src/atm_phy_schemes/mo_cuparameters.f90
src/atm_phy_schemes/mo_gwd_wms.f90
src/atm_phy_schemes/mo_o3_gems_data.f90
src/atm_phy_schemes/mo_o3_macc_data.f90
src/atm_phy_schemes/mo_o3_util.f90
src/atm_phy_schemes/mo_lrtm_coeffs.f90
src/atm_phy_schemes/mo_lrtm_driver.f90
src/atm_phy_schemes/mo_lrtm_kgs.f90
src/atm_phy_schemes/mo_lrtm_netcdf.f90
src/atm_phy_schemes/mo_lrtm_par.f90
src/atm_phy_schemes/mo_lrtm_rtrnmr.f90
src/atm_phy_schemes/mo_lrtm_setup.f90
src/atm_phy_schemes/mo_lrtm_taumol.f90
src/atm_phy_schemes/mo_srtm_config.f90
src/atm_phy_schemes/mo_srtm.f90
src/atm_phy_schemes/mo_srtm_kgb16.f90
src/atm_phy_schemes/mo_srtm_kgb17.f90
src/atm_phy_schemes/mo_srtm_kgb18.f90
src/atm_phy_schemes/mo_srtm_kgb19.f90
src/atm_phy_schemes/mo_srtm_kgb20.f90
src/atm_phy_schemes/mo_srtm_kgb21.f90
src/atm_phy_schemes/mo_srtm_kgb22.f90
src/atm_phy_schemes/mo_srtm_kgb23.f90
src/atm_phy_schemes/mo_srtm_kgb24.f90
src/atm_phy_schemes/mo_srtm_kgb25.f90
src/atm_phy_schemes/mo_srtm_kgb26.f90
src/atm_phy_schemes/mo_srtm_kgb27.f90
src/atm_phy_schemes/mo_srtm_kgb28.f90
src/atm_phy_schemes/mo_srtm_kgb29.f90
src/atm_phy_schemes/mo_srtm_kgs.f90
src/atm_phy_schemes/mo_srtm_taumol.f90
src/atm_phy_schemes/mo_sso_cosmo.f90
src/atm_phy_schemes/mo_sso_ifs.f90
support/mo_util_uuid.f90
support/mo_util_uuid_types.f90
src/atm_phy_schemes/mo_voskin.f90
src/atm_phy_schemes/mo_vdftofdc.f90
src/atm_phy_schemes/random_rewrite.f90
src/atm_phy_rte_rrtmgp/mo_cloud_optics.f90
src/atm_phy_rte_rrtmgp/mo_load_cloud_coefficients.f90
src/atm_phy_rte_rrtmgp/mo_simple_netcdf.f90
src/atm_phy_rte_rrtmgp/mo_load_coefficients.f90
'

# Number of parallel jobs:
job_num=8

set -eu
set -o pipefail

quote_value ()
{
  eval "qv_x=\$${2-${1}}"
  case ${qv_x} in
    *\'*) qv_x=`echo "${qv_x}" | sed "s/'/'\\\\\\\\''/g"` ;;
  esac
  eval "${1}=\"'\${qv_x}'\""
  unset qv_x
}

flatten_and_quote_values ()
{
  eval "faqv_input=\$${2-${1}}; ${1}=''"
  while read faqv_x; do
    test -n "${faqv_x}" || continue
    quote_value faqv_x
    eval "${1}=\"\$${1}${faqv_x} \""
  done <<_EOF
${faqv_input}
_EOF
  unset faqv_input faqv_x
}

check_exist()
{
  for ce_input in "${@}"; do
    if test ! -e "${ce_input}"; then
      echo "ERROR: '${ce_input}' does not exist" >&2
      exit 2
    fi
  done
  unset ce_input
}

list_files()
{
  lf_issue_warn=${1}; shift
  lf_dirs=''
  for lf_input in "${@}"; do
    if test -d "${lf_input}"; then
      quote_value lf_input
      lf_dirs="${lf_dirs} ${lf_input}"
    elif test -f "${lf_input}"; then
      echo "${lf_input}"
    elif test x"${issue_warn}" = xyes; then
      echo "WARNING: input argument '${lf_input}' is neither a directory nor a file" >&2
    fi
  done
  if test -n "${lf_dirs}"; then
    lf_find_args="${lf_dirs} -type f -a \\("
    while read lf_x; do
      test -n "${lf_x}" || continue
      quote_value lf_x
      lf_find_args="${lf_find_args} -name ${lf_x} -o"
    done <<_EOF
${known_patterns}
_EOF
    lf_find_args="${lf_find_args} -false \\) -a \\! \\("
    eval "set dummy ${ignored_patterns}; shift"
    for lf_x in "${@}"; do
      quote_value lf_x
      lf_find_args="${lf_find_args} -path ${lf_x} -o"
    done
    lf_find_args="${lf_find_args} -false \\)"
    eval "find${lf_find_args}"
    unset lf_find_args
  fi
  unset lf_issue_warn lf_dirs lf_input
}

ignored_patterns=''
if test ${#} -eq 0; then
  icon_prefix="$(unset CDPATH; cd "$(dirname "${0}")/../.."; pwd)/"
  quote_value icon_prefix
  flatten_and_quote_values icon_ignored_patterns
  eval "set dummy ${icon_ignored_patterns}; shift"
  for pattern in "${@}"; do
    quote_value pattern
    ignored_patterns="${ignored_patterns} ${icon_prefix}${pattern}"
  done
  flatten_and_quote_values icon_directories
  eval "set dummy ${icon_directories}; shift"
  icon_directories=''
  for dir in "${@}"; do
    quote_value dir
    icon_directories="${icon_directories} ${icon_prefix}${dir}"
  done
  eval "set dummy ${icon_directories}; shift"
fi

check_exist "${@}"

# Prepend the comment character to each file of the ${license}:
license=$(echo -n "${license}" | sed "s/^/${comment}/")

# Prepare the regular expression based on the ${license}:
#   - escape the metacharacters;
#   - replace all new lines with \01.
regexp=$(echo -n "${license}" | sed 's/[][()\.^$?*+]/\\&/g' | tr '\n' '\01')

# Prepend the regular expression so that the source files are allowed to have
# only commented, preprocessor (start with #) and space-only lines before the
# license text:
regexp='^([[:space:]]*(('"${comment}"'|#)[^\01]*)?\01)*?'"${regexp}"
quote_value regexp

exitcode=0

list_files yes "${@}" | xargs -P ${job_num} -I{} -- ${SHELL-${BASH}} -c "{ cat '{}' | tr '\\n' '\\01' | grep -qPz ${regexp}; } || { echo '{}'; exit 1; }" || {
  {
    echo "ERROR: some of the input files (see above) do not contain the following license"
    echo "       header):"
    echo "${license}"
  } >&2
  exitcode=1
}

# Check also for forbidden Doxygen entries (the following regular expression is
# used in the case-insensitive mode):
regexp='^[[:space:]]*'"${comment}"'.*@(par[[:space:]]+revision[[:space:]]+history|author)'
quote_value regexp

list_files no "${@}" | xargs -P ${job_num} -I{} -- ${SHELL-${BASH}} -c "grep --color='auto' -HniP ${regexp} '{}' >&2; test \${?} -eq 1" || {
  {
    echo "ERROR: some of the input files contain forbidden Doxygen commands"
    echo "       '@par Revision history' and '@author' (see above)"
    echo "       remove the forbidden Doxygen commands"
  } >&2
  exitcode=1
}

exit ${exitcode}
