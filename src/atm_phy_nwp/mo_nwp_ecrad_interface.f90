!
! This module is the interface between ICON:nwp_radiation to the radiation scheme ecRad
!
! - There are two interfaces within this module: nwp_ecrad_radiation and
!   nwp_ecrad_radiation_reduced. The former provides the interface to ecrad on the full
!   ICON dynamics grid. The latter one provides the interface on a grid with a reduced
!   spatial resolution.
! - The decision which of the two interface routines is used, is done via the namelist
!   switch lredgrid_phys. Based on the value of lredgrid_phys, the correct interface
!   is called by mo_nwp_rad_interface:nwp_radiation.
! - The interfaces have to fill the different ecRad input types (ecrad_aerosol,
!   ecrad_single_level, ecrad_thermodynamics, ecrad_gas, ecrad_cloud) with data from
!   ICON variables. Then, the ecRad radiation code is called. At the end, the fluxes
!   calculated by ecRad and stored in the ecrad_flux structure are copied to ICON variables.
! - The difference between nwp_ecrad_radiation and nwp_ecrad_radiation_reduced is mostly
!   an upscaling at the beginning and a downscaling at the end of the interface.
! - The transfer of data from ICON to ecRad and vice versa is performed within
!   routines from mo_nwp_ecrad_utilities and mo_nwp_ecrad_prep_aerosol, independent of
!   the choice to use a reduced radiation grid or not.
!
!
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

!----------------------------
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

!----------------------------

MODULE mo_nwp_ecrad_interface

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message
  USE mo_math_constants,         ONLY: pi
  USE mo_model_domain,           ONLY: t_patch, p_patch_local_parent
  USE mo_impl_constants,         ONLY: min_rlcell_int, n_camsaermr
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c, grf_ovlparea_start_c, grf_fbk_start_c
  USE mo_fortran_tools,          ONLY: init, assert_acc_device_only,t_ptr_2d
  USE mo_parallel_config,        ONLY: nproma, nproma_sub, nblocks_sub
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_grid_config,            ONLY: l_limited_area, nexlevs_rrg_vnest
  USE mo_ext_data_types,         ONLY: t_external_data
  USE mo_nwp_lnd_types,          ONLY: t_lnd_prog
  USE mo_nonhydro_types,         ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,          ONLY: t_nwp_phy_diag
  USE mo_physical_constants,     ONLY: rhoh2o
  USE mo_run_config,             ONLY: msg_level, iqv, iqi, iqc, iqr, iqs, iqg
  USE mo_atm_phy_nwp_config,     ONLY: atm_phy_nwp_config
  USE mo_radiation_config,       ONLY: irad_aero, ssi_radt,                                   &
                                   &   iRadAeroNone, iRadAeroConst, iRadAeroTegen,            &
                                   &   iRadAeroART, iRadAeroConstKinne, iRadAeroKinne,        &
                                   &   iRadAeroVolc, iRadAeroKinneVolc,  iRadAeroKinneVolcSP, &
                                   &   iRadAeroKinneSP, iRadAeroCAMSclim
  USE mo_phys_nest_utilities,    ONLY: t_upscale_fields, upscale_rad_input, downscale_rad_output
  USE mtime,                     ONLY: datetime


  IMPLICIT NONE

  PRIVATE
END MODULE mo_nwp_ecrad_interface
