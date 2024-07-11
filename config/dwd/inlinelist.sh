# ICON
#
# ------------------------------------------
# Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ------------------------------------------

# Join all arguments using the separator given by the first argument, prefixing each with the second.
# Usage: join_arr <sep> <prefix> <item0> <item1>...
function join_arr {
  sep="$1"
  prefix="$2"
  shift 2
  while (( $# > 1 )); do
    echo -n "$prefix$1$sep"
    shift
  done
  echo "$prefix$1"
}


INLINE_LIST_ICON=(
  src/advection/mo_advection_utils.f90
  src/atm_phy_nwp/mo_util_phys.f90
  src/atm_phy_schemes/cloud_random_numbers.f90
  src/atm_phy_schemes/mo_2mom_mcrph_driver.f90
  src/atm_phy_schemes/mo_2mom_mcrph_processes.f90
  src/atm_phy_schemes/mo_2mom_mcrph_util.f90
  src/atm_phy_schemes/mo_aerosol_sources.f90
  src/atm_phy_schemes/mo_cpl_aerosol_microphys.f90
  src/atm_phy_schemes/mo_albedo.f90
  src/atm_phy_schemes/mo_cufunctions.f90
  src/atm_phy_schemes/mo_satad.f90
  src/atm_phy_schemes/mo_turb_vdiff.f90
  src/atm_phy_schemes/random_rewrite.f90
  src/atm_phy_schemes/turb_utilities.f90
  src/configure_model/mo_parallel_config.f90
  src/lnd_phy_nwp/mo_nwp_sfc_interp.f90
  src/parallel_infrastructure/mo_extents.f90
  src/shared/mo_math_utilities.f90
  src/shared/mo_statistics.f90
)


INLINE_LIST_ART=(
  externals/art/aerosol_dynamics/mo_art_aerosol_utilities.f90
  externals/art/shared/mo_art_modes.f90
  externals/art/tools/mo_art_clipping.f90
)


INLINE_LIST_DACE=(
  externals/dace_icon/src_for_icon/mo_physics.f90
)


INLINE_LIST_ECRAD=(
  externals/ecrad/radiation/radiation_two_stream.F90
  externals/ecrad/radiation/radiation_liquid_optics_socrates.F90
  externals/ecrad/radiation/radiation_ice_optics_fu.F90
)


INLINE_LIST_EMVORADO=(
  externals/emvorado/src_emvorado/radar_gamma_functions_vec.f90
  externals/emvorado/src_emvorado/radar_mie_meltdegree.f90
  externals/emvorado/src_emvorado/radar_mie_specint.f90
  externals/emvorado/src_emvorado/radar_mie_utils.f90
  externals/emvorado/src_emvorado/radar_mielib_vec.f90
  externals/emvorado/src_emvorado/radar_model2rays.f90
  externals/emvorado/src_emvorado/radar_utilities.f90
  externals/emvorado/src_emvorado/radar_dmin_wetgrowth.f90
  externals/emvorado/src_iface_icon/radar_interface.f90
)
