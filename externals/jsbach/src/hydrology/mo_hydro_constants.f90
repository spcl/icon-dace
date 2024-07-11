!> Contains constants for the radiation processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_hydro_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER ::                    &
    ! Identifiers used for the soil hydrology model
    & BrooksCorey_     = 1,                &
    & ClappHornberger_ = 2,                &
    & VanGenuchten_    = 3,                &
    ! Identifiers used for the interpolation scheme
    & Upstream_        = 1,                &
    & Arithmetic_      = 2,                &
    ! Identifiers used for the pond dynamics scheme
    & Quad_            = 1,                &
    & Tanh_            = 2

  REAL(wp), PARAMETER ::                   &
    ! Parameters used for snow cover fraction
    & wsn2fract_eps    = 1.E-12_wp,        &
    & wsn2fract_sigfac = 0.15_wp,          &
    & wsn2fract_const  = 0.95_wp,          & !< Inverse of snow height when snow is considered to completely cover the ground

    & InterceptionEfficiency = 0.25_wp,    & !< Efficiency of interception of precipitation (rain and snow)

    & Epar                   = 2.2E5_wp,   & ! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1.  R: for update_par
    & SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region.  R: for update_par
!    & LaiMax = 8._wp,                      & ! Maximum LAI (used for nitrogen scaling)
    & FcMax = 1.0_wp,                      & ! Maximum fractional vegetation cover. R: for calc_fapar_lai_cl
    & FcMin = 1.E-3_wp,                    & ! Minimum fractional vegetation cover. R: for calc_fapar_lai_cl
    & ZenithMinPar = 1.E-3_wp,             &  ! Minimum cos of zenith angle for which PAR is calculated. R: for calc_fapar_lai_cl
!    & ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
!    & SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region
!
!    & minStomaConductance = 0.0_wp,        & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]

    ! Parameters describing orographic steepness in infiltration computation
    & oro_var_min = 100._wp,               &
    & oro_var_max = 1000._wp,              &

    ! Parameters used for pond computations
    & oro_crit = 100._wp,                   & !< Reference topo standard deviation for pond scaling factor [m]

    ! Parameters for the computation of canopy conductance/resistance using Eq. 3.3.2.12 in ECHAM3 manual
    & conductance_k = 0.9_wp,              & !<
    & conductance_a = 5000._wp,            & !< Jm-3
    & conductance_b = 10._wp,              & !< Wm-2
    & conductance_c = 100._wp,             & !< ms-1

    ! Parameters for organic soil component (top layer)
    & vol_porosity_org_top    = 0.95_wp,       & !< volumetric porosity of organic part of soil layer [m/m]
    & vol_field_cap_org_top   = 0.95_wp,       & !< volumetric field capacity of organic layer [m/m]
    & vol_p_wilt_org_top      = 0.255_wp,      & !< volumetricwilting point of organic layer [m/m]
    & hyd_cond_sat_org_top    = 0.0001_wp,     & !< Saturated hydraulic conductivity of organic part of soil layer [m/s]
    & bclapp_org_top          = 4._wp,         & !< Exponent b in Clapp and Hornberger of organic part of soil layer  []
    & matrix_pot_org_top      = -0.1_wp,       & !< Matrix potential of organic part of soil layer [m]
                                                 ! (values (roughy) from Chadburn 2022)
    & pore_size_index_org_top = 0.7_wp,        & !< Pore size distribution index of organic part of soil layer []
    ! Parameters for organic soil component (deep layers)
    & vol_porosity_org_below  = 0.82_wp,       & !< volumetric porosity of organic part of soil layer [m/m]
    & vol_field_cap_org_below = 0.82_wp,       & !< volumetric field capacity of organic layer [m/m]
    & vol_p_wilt_org_below    = 0.255_wp,      & !< volumetricwilting point of organic layer [m/m]
    & hyd_cond_sat_org_below  = 0.00000001_wp, & !< Saturated hydraulic conductivity of organic part of soil layer [m/s]
    & bclapp_org_below        = 8._wp,         & !< Exponent b in Clapp and Hornberger of organic part of soil layer  []
    & matrix_pot_org_below    = -1.0_wp,       & !< Matrix potential of organic part of soil layer [m]
                                                 ! (values (roughy) from Chadburn 2022)
    & pore_size_index_org_below = 0.7_wp,      & !< Pore size distribution index of organic part of soil layer []
    & thresh_org              = 0.5_wp,        & !< Threshold of organic material above which connected flow pathways form
    & beta_perc               = 0.139_wp         !< Parameter from percolation theory

#endif
END MODULE mo_hydro_constants
