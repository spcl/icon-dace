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

! In this module the configuration state for the ecRad radiation code is being set up.
!
! - Setup information is stored inside the object ecrad_conf of the derived type t_ecrad_conf
!   containing all the configuration information needed to run the radiation scheme.
! - The intention is that this part of the configuration is fixed for a given model run.
! - ICON namelist settings are translated to ecRad conform settings, if unsupported values
!   are provided via namelist, the user gets an error. (These values should already be
!   checked by the nml_crosscheck)
! - Please note that only a subset of the available configuration options of ecRad is
!   filled by this routine. For a full list of ecRad settings, please have a look at
!   externals/ecrad/radiation/radiation_config.F90

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

MODULE mo_nwp_ecrad_init

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: finish, message
  USE mo_radiation_config,     ONLY: icld_overlap, irad_aero, ecrad_data_path,             &
                                 &   ecrad_isolver, ecrad_igas_model, isolrad,             &
                                 &   ecrad_llw_cloud_scat, ecrad_iliquid_scat,             &
                                 &   ecrad_iice_scat, ecrad_isnow_scat,                    &
                                 &   ecrad_irain_scat, ecrad_igraupel_scat,                &
                                 &   ecrad_use_general_cloud_optics,                       &
                                 &   iRadAeroConst, iRadAeroTegen, iRadAeroART,            &
                                 &   iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc,      &
                                 &   iRadAeroKinneVolc,  iRadAeroKinneVolcSP,              &
                                 &   iRadAeroKinneSP, iRadAeroNone,                        &
                                 &   iRadAeroCAMSclim, iRadAeroCAMStd


  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_init'

END MODULE mo_nwp_ecrad_init
