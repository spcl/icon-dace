!> QUINCY assimilation parameters
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and init assimilation parameters
!>
MODULE mo_q_assimi_parameters
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
  REAL(wp), SAVE :: &
    ! Values from Niinemets 1997 (SC 01/16)
    & jmax2n = def_parameters, &           !< electron-transport limited carboxylation rate per micro-mol CO2 / mmol N in ET
    & vcmax2n = def_parameters, &          !< Rubisco limited carboxylation rate per micro-mol CO2 / mmol N in RUB
    ! PEP C limited carboxylation rate per micro-mol CO2 / mmol N in PEP C, derived from k = 0.7 mol m-2 s-1 (Collatz et al. 1992)
    ! assuming 157 mmol N / m-2 leaf N and 4.5 % of leaf N in PS
    & pepc2n = def_parameters, &           !< PEP C limited carboxylation rate per micro-mol CO2 / mmol N in PEP C
    ! Value from Friend et al. 2009
    ! values from Evans 1989, Oecologia, 78: 9-19 (Table 2)
    & chl2n = def_parameters, &            !< micro-mol Chloroplast per mmol N (light harvesting complex)
    ! slope of the relationship of structural leaf N against leaf N in fraction/mmol (Friend 1997)
    & k1_fn_struc = def_parameters         !< slope of the fraction of leaf N not associated with PS against leaf N (--/mmolN)

  ! Parameters used in the photosynthesis routine (temperature sensitivies of PS processes)
  REAL(wp), SAVE :: &
    & E0kc = def_parameters, &             !< Scaling constant for Kc (unitless)
    & E1kc = def_parameters, &             !< Activation energy of Kc (J mol-1)
    & E0ko = def_parameters, &             !< Scaling constant for Ko (unitless)
    & E1ko = def_parameters, &             !< Activation energy of Ko (J mol-1)
    & E0pcp = def_parameters, &            !< Scaling constant for the photosynthetic compensation point (unitless)
    & E1pcp = def_parameters, &            !< Activation energy of the photosynthetic compensation point (J mol-1)
    & E0v = def_parameters, &              !< Scaling constant for Rubisco (unitless)
    & E1v = def_parameters, &              !< Activation energy of Rubisco (J mol-1)
    ! Temperature sensitivity of Electron transport
    & t_jmax_offset = def_parameters, &    !< Offset of the Tjmax~Tair relationship
    & t_jmax_slope = def_parameters, &     !< slope of the Tjmax~Tair relationship (candidate for PFT specific parameterisation?)
    & t_jmax_opt_min = def_parameters, &   !< Minimum of optimum Tjmax
    & t_jmax_opt_max = def_parameters, &   !< Maximum of optimum Tjmax
    ! Temperature sensitivity of PeP C4 photosynthesis
    & Tref_pepc = def_parameters, &        !< Temperature sensitivity of PeP C4 photosynthesis
    & Tbase_pepc = def_parameters          !< Temperature sensitivity of PeP C4 photosynthesis

  REAL(wp), SAVE :: &
    ! intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I))) (micro-mol CO2 / mol quanta)
    ! accounting for the absorptance of leaves (0.85), the maximum quantum yield of PS II (0.7)
    ! the fraction of total light that reaches PS II (0.5) and the NADPH limited RubGeneration (jointly leading to 0.066)
    & alpha_i = def_parameters, &            !< intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I)) (micro-mol CO2 / mol quanta)
    ! Partial Pressure of O2 in kPa
    & pO2 = def_parameters, &                !< Partial Pressure of O2 in kPa
    ! initial guess of the leaf Ci to Ca ratio
    & CiCa_default_C3 = def_parameters, &    !< default Ci:Ca for C3 plants
    & CiCa_default_C4 = def_parameters, &    !< default Ci:Ca for C4 plants
    ! maximum number of iteration in photosynthesis calculation (should probably be in mo_X_constants)
    & ps_it_max = def_parameters, &          !< maximum iteration of PS calculation
    & ci_max = def_parameters                !< saturating Ci in Pa C4 plants

  ! parameters describing the isotopic discrimination of PS (fractionation due to photosynthesis), derived from Drake 2014, Radiocarbon, 56, 29-38
  REAL(wp), SAVE :: &
    & discr_ps_a_C13 = def_parameters, &     !< discrimination of C13 due to stomatal diffusion
    & discr_ps_b_C13 = def_parameters, &     !< discrimination of C13 due to Rubisco
    & discr_ps_c_C13 = def_parameters, &     !< discrimination of C13 due to PEP C
    & discr_ps_a_C14 = def_parameters, &     !< discrimination of C14 due to stomatal diffusion
    & discr_ps_b_C14 = def_parameters, &     !< discrimination of C14 due to Rubisco
    & discr_ps_c_C14 = def_parameters, &     !< discrimination of C14 due to PEP C
    & discr_ps_phi   = def_parameters        !< leakage rate of bundle sheath cells; a typical value cf Hatch et al. 1995, Plant Phys. 108, 173-181

  ! parameters for canopy profile and canopy light extinction (parameters describing the canopy profile of N)
  REAL(wp), SAVE :: &
    & ka = def_parameters, &                 !< Extinction coefficient for PAR on chlorophyll
    & kn = def_parameters                    !< extinction coefficient to describe decline of N within the canopy

  REAL(wp), SAVE :: &
    & soa_b = def_parameters, &              !< parameter in state of acclimation calculation used for evergreen forests
    & soa_t_s = def_parameters               !< parameter in state of acclimation calculation used  for evergreen forests

  ! parameters for calculation of chlorophyll fluorescence
  REAL(wp), SAVE:: &
    & heat_dissipation_enrg_dep_param_1               = def_parameters, &   !< parameter for kn (energy dependent heat dissipation) calculation, source: drought dataset (Flexas 2002), for cotton kno=2.48 (Weiss and Berry 1987)
    & heat_dissipation_enrg_dep_param_2               = def_parameters, &   !< parameter for kn calculation, source: drought dataset (Flexas 2002), for cotton alpha_kn=2.83 (Weiss and Berry 1987)
    & heat_dissipation_enrg_dep_param_3               = def_parameters, &   !< parameter for kn calculation, source: drought dataset (Flexas 2002), for cotton beta_kn=0.114 (Weiss and Berry 1987)
    & k_fluorescence_rate                             = def_parameters, &   !< rate coefficient of fluorescence
    & thermal_dissipation_enrg_constitutive_max_rate  = def_parameters, &   !< maximum value for rate constant of energy constitutive thermal dissipation
    & thermal_dissipation_enrg_constitutive_param_1   = def_parameters, &   !< parameter for calculating the rate coefficient for thermal dissipation
    & thermal_dissipation_enrg_constitutive_param_2   = def_parameters, &   !< parameter for calculating the rate coefficient for thermal dissipation
    & k_photochemistry_rate                           = def_parameters      !< rate coefficient for photochemistry


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_assimi_parameters'

CONTAINS

  ! ======================================================================================================= !
  !> initialize parameters for the process: assimilation
  !>
  SUBROUTINE init_q_assimi_parameters
    USE mo_jsb_physical_constants,  ONLY: molar_mass_N    !< molar mass of nitrogen
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_q_assimi_parameters'

    ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
    jmax2n        = 17.6_wp / 4.0_wp                      !< .. reference
    vcmax2n       = 1.8_wp                                !< .. reference
    pepc2n        = 98777.97_wp                           !< .. reference
    chl2n         = 1000._wp / 39.8_wp                    !< .. reference
    k1_fn_struc   = 71.4_wp * molar_mass_N / 1000000._wp  !< .. reference

    ! Parameters used in the photosynthesis routine (temperature sensitivies of PS processes)
    E0kc            = 38.05_wp                          !< Bernacchi 2001
    E1kc            = 79430.0_wp                        !< Bernacchi 2001
    E0ko            = 20.3_wp                           !< Bernacchi 2001
    E1ko            = 36380.0_wp                        !< Bernacchi 2001
    E0pcp           = 19.02_wp                          !< Bernacchi 2001
    E1pcp           = 37830.0_wp                        !< Bernacchi 2001
    E0v             = 26.35_wp                          !< Bernacchi 2001
    E1v             = 65330.0_wp                        !< Bernacchi 2001
    t_jmax_offset   = 17.0_wp                           !< Friend 2009
    t_jmax_slope    = 0.35_wp                           !< Friend 2009
    t_jmax_opt_min  = t_jmax_offset                     !<
    t_jmax_opt_max  = 38.0_wp                           !< Friend 2009
    Tref_pepc       = 25._wp                            !< Friend 2009
    Tbase_pepc      = 10._wp                            !< Friend 2009
    alpha_i         = 0.85_wp * 0.066_wp                !< .. reference
    pO2             = 20.9_wp                           !< .. reference
    CiCa_default_C3 = 0.7_wp                            !< .. reference
    CiCa_default_C4 = 0.4_wp                            !< .. reference
    ps_it_max       = 10._wp                            !< .. reference
    ci_max          = 7800.0_wp                         !< .. reference

    ! parameters describing the isotopic discrimination of PS (fractionation due to photosynthesis)
    discr_ps_a_C13 = 4.4_wp                             !< derived from Drake 2014, Radiocarbon, 56, 29-38
    discr_ps_b_C13 = 27.0_wp                            !< derived from Drake 2014, Radiocarbon, 56, 29-38
    discr_ps_c_C13 = 5.7_wp                             !< derived from Drake 2014, Radiocarbon, 56, 29-38
    discr_ps_a_C14 = 1.97_wp * discr_ps_a_C13           !< derived from Drake 2014, Radiocarbon, 56, 29-38
    discr_ps_b_C14 = 1.89_wp * discr_ps_b_C13           !< derived from Drake 2014, Radiocarbon, 56, 29-38
    discr_ps_c_C14 = 1.89_wp * discr_ps_c_C13           !< derived from Drake 2014, Radiocarbon, 56, 29-38
    discr_ps_phi   = 0.16_wp                            !< derived from Drake 2014, Radiocarbon, 56, 29-38

    ! parameters for canopy profile and canopy light extinction (parameters describing the canopy profile of N)
    ka            = 0.005_wp                            !< .. reference
    kn            = 0.11_wp                             !< .. reference

    ! ...
    soa_b         = -0.5_wp                             !< .. reference
    soa_t_s       = 5.0_wp                              !< .. reference

    ! parameters for chlorophyll fluorescence
    heat_dissipation_enrg_dep_param_1               = 5.01_wp     !< Flexas 2002
    heat_dissipation_enrg_dep_param_2               = 1.93_wp     !< Flexas 2002
    heat_dissipation_enrg_dep_param_3               = 10.0_wp     !< Flexas 2002
    k_fluorescence_rate                             = 0.05_wp     !< van der Tol 2014
    thermal_dissipation_enrg_constitutive_max_rate  = 0.8738_wp   !< van der Tol 2014
    thermal_dissipation_enrg_constitutive_param_1   = 0.0301_wp   !< van der Tol 2014
    thermal_dissipation_enrg_constitutive_param_2   = 0.0773_wp   !< van der Tol 2014
    k_photochemistry_rate                           = 4.0_wp      !< van der Tol 2014

  END SUBROUTINE init_q_assimi_parameters

#endif
END MODULE mo_q_assimi_parameters
