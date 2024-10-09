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

MODULE mo_oem_config

  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: filename_max
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Basic configuration setup for the online emission module
  !--------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! oemctrl_nml:
    !------------------------------------------------------------------------
    CHARACTER(LEN=filename_max) :: vertical_profile_nc,   & !< name of the oae vertical profile
                              &    hour_of_day_nc,        & !< name of the oae hour of day file
                              &    day_of_week_nc,        & !< name of the oae day of week file
                              &    month_of_year_nc,      & !< name of the oae month of year file
                              &    hour_of_year_nc,       & !< name of the oae hour of year file
                              &    gridded_emissions_nc,  & !< name of the oae gridded emission file
                              &    ens_reg_nc,            & !< name of file with ensemble-regions
                              &    ens_lambda_nc,         & !< name of file with ensemble-lambdas
                              &    vegetation_indices_nc    !< name of file with MODIS reflectances
    REAL(wp), DIMENSION(8) ::      vprm_par,              & !< VPRM parameter values for PAR_0
                              &    vprm_lambda,           & !< VPRM parameter values for lambda
                              &    vprm_alpha,            & !< VPRM parameter values for alpha
                              &    vprm_beta,             & !< VPRM parameter values for beta
                              &    vprm_tmin,             & !< VPRM parameter values for T_min
                              &    vprm_tmax,             & !< VPRM parameter values for T_max
                              &    vprm_topt,             & !< VPRM parameter values for T_opt
                              &    vprm_tlow                !< VPRM parameter values for T_low
    LOGICAL ::                     lcut_area                !< Switch to turn on/off to select an
                                                            !< area where no fluxes are applied
    REAL(wp) ::                    lon_cut_start,         & !< longitude start coordinate 
                              &    lon_cut_end,           & !< longitude end coordinate
                              &    lat_cut_start,         & !< latitude start coordinate
                              &    lat_cut_end              !< latitude end coordinate
 



END MODULE mo_oem_config

