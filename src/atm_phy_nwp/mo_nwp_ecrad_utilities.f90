!
! The radiation scheme ecRad expects information from the host model (i.e., ICON) to be
!   copied to different ecRad data structures: t_ecrad_single_level_type, t_ecrad_gas_type,
!   t_ecrad_thermodynamics_type and t_ecrad_cloud_type. Similarly, the output of ecRad, i.e.
!   the radiative fluxes, are stored in a data structure t_ecrad_flux_type.
! This module offers subroutines that transfer the data from ICON variables into the
!   correct ecRad data structure. The subroutines are written in a way that they can be used
!   on the reduced radiation grid as well as on the full radiation grid. This ensures
!   consistency between the two modes.
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_nwp_ecrad_utilities

  USE mo_kind,                   ONLY: wp
  USE mo_math_constants,         ONLY: rad2deg, pi
  USE mo_exception,              ONLY: finish
  USE mo_fortran_tools,          ONLY: set_acc_host_or_device, assert_acc_device_only
  USE mo_math_types,             ONLY: t_geographical_coordinates
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_physical_constants,     ONLY: rd, grav
  USE mo_radiation_config,       ONLY: vmr_co2, vmr_n2o, vmr_o2, vmr_ch4,        &
                                   &   vmr_cfc11, vmr_cfc12,                     &
                                   &   irad_h2o, irad_o3, irad_co2,              &
                                   &   irad_n2o, irad_ch4,                       &
                                   &   irad_o2, irad_cfc11, irad_cfc12,          &
                                   &   vpp_ch4, vpp_n2o, decorr_pole, decorr_equator
  USE mo_nwp_tuning_config,      ONLY: tune_difrad_3dcont
  USE mtime,                     ONLY: datetime
  USE mo_bc_greenhouse_gases,    ONLY: ghg_co2mmr, ghg_ch4mmr, ghg_n2ommr, ghg_cfcmmr
  USE mo_parallel_config,        ONLY: nproma, nproma_sub

  USE mo_exception,              ONLY: message
  USE mo_grid_config,            ONLY: l_scm_mode
  USE mo_scm_nml,                ONLY: lon_scm, lat_scm 

  IMPLICIT NONE

  PRIVATE
END MODULE mo_nwp_ecrad_utilities
