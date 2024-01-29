!> Contains constants for the assimilation processes
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
MODULE mo_assimi_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                                    &
         ! C3 PLANTS: FARQHUAR, G.D., S. VON CAEMMERER AND J.A. BERRY, 1980.
         !            A BIOCHEMICAL MODEL OF PHOTOYNTHESIS IN LEAVES OF C3 SPECIES.
         !            PLANTA 149, 78-90.
                         & ALPHA = 0.28_wp,                 & ! EFFICIENCY OF OF PHOTON CAPTURE
                         & OX    = 0.21_wp,                 & ! OXYGEN CONCENTRATION [MOL(O2) / MOL(AIR)]
                         & KC0   = 460.E-6_wp,              & ! MICHAELIS-MENTEN CONSTANT FOR CO2 AT 25C [MOL(CO2) / MOL(AIR)]
                         & KO0   = 330.E-3_wp,              & ! MICHAELIS-MENTEN CONSTANT FOR O2 AT 25C [MOL(O2) / MOL(AIR)]
                         & EC    = 59356._wp,               & ! ACTIVATION ENERGY FOR KC [J / MOL]
                         & EO    = 35948._wp,               & ! ACTIVATION ENERGY FOR KO [J / MOL]
                         & EV    = 58520._wp,               & ! ACTIVATION ENERGY FOR VCMAX [J / MOL]
                         & ER    = 45000._wp,               & ! ACTIVATION ENERGY FOR DARK RESPIRATION [J / MOL]
                         & EK    = 50967._wp,               & !  = Q10=2 (Collatz et al. 1992)
                         & FRDC3 = 0.011_wp,                & ! RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C3
         ! C4 PLANTS: COLLATZ, G.J., M. RIBAS-CARBO AND J.A. BERRY, 1992.
         !            COUPLED PHOTOSYNTHESIS-STOMATAL CONDUCTANCE MODEL FOR LEAVES
         !            OF C4 PLANTS. AUST. J. PLANT PHYSIOL. 19, 519-538.
                         & FRDC4 = 0.031_wp,                & !RATIO OF DARK RESPIRATION TO "PVM" AT 25C for C4
                         & ALC4  = 0.04_wp,                 & !EFFECTIVE QUANTUM EFFICIENCY
                         & THETA = 0.83_wp,                 & !CURVATURE PARAMETER
                         & minOfMaxCarboxrate = 1.0e-12_wp, & ! Minimum of maximum carboxylation rate [10^(-6) mol/(m^2 s)].
         ! both
                         & minStomaConductance = 0.0_wp,    & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]
         ! Factors that relates leaf internal CO2-concentration to CO2-concentration of ambient air:
                         & FCI1C3        = 0.87_wp,         & ! For C3 plants
                         & FCI1C4        = 0.67_wp,         & ! For C4 plants
         ! Factors needed to calculate NPP_pot_rate:
                         & f_aut_leaf    = 0.40_wp,         & ! leaf fraction of plant-total (autotrophic) respiration
                         & cCost         = 1.25_wp    ! relative costs (measured in carbon) to produce 1 unit of carbon

!! quincy
  ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
  REAL(wp), SAVE :: &      ! In Fortran 2003 and Fortran 2008, a maximum of 255 continuation lines are allowed
      ! Values from Niinemets 1997 (SC 01/16)
      jmax2n        ,     &  !< electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
      vcmax2n       ,     &  !< Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
      ! PEP C limited carboxylation rate per µmol CO2 / mmol N in PEP C, derived from k = 0.7 mol m-2 s-1 (Collatz et al. 1992)
      ! assuming 157 mmol N / m-2 leaf N and 4.5 % of leaf N in PS 
      pepc2n        ,     & !< PEP C limited carboxylation rate per µmol CO2 / mmol N in PEP C
      ! Value from Friend et al. 2009 
      ! values from Evans 1989, Oecologia, 78: 9-19 (Table 2)
      chl2n         ,     & !< µmol Chloroplast per mmol N (light harvesting complex)
      ! slope of the relationship of structural leaf N against leaf N in fraction/mmol (Friend 1997)
      k1_fn_struc           !< slope of the fraction of leaf N not associated with PS against leaf N (--/mmolN) 

  ! Parameters used in the photosynthesis routine (temperature sensitivies of PS processes)
  REAL(wp), SAVE :: &
      E0kc    ,     & !< Scaling constant for Kc (unitless) 
      E1kc    ,     & !< Activation energy of Kc (J mol-1) 
      E0ko    ,     & !< Scaling constant for Ko (unitless)
      E1ko    ,     & !< Activation energy of Ko (J mol-1)
      E0pcp   ,     & !< Scaling constant for the photosynthetic compensation point (unitless)
      E1pcp   ,     & !< Activation energy of the photosynthetic compensation point (J mol-1)
      E0v     ,     & !< Scaling constant for Rubisco (unitless)
      E1v     ,     & !< Activation energy of Rubisco (J mol-1)
      ! Temperature sensitivity of Electron transport
      t_jmax_offset, & !< Offset of the Tjmax~Tair relationship
      t_jmax_slope, &  !< slope of the Tjmax~Tair relationship (candidate for PFT specific parameterisation?)
      t_jmax_opt_min, & !< Minimum of optimum Tjmax 
      t_jmax_opt_max, & !< Maximum of optimum Tjmax
      ! Temperature sensitivity of PeP C4 photosynthesis
      Tref_pepc   ,     & !< Temperature sensitivity of PeP C4 photosynthesis 
      Tbase_pepc         !< Temperature sensitivity of PeP C4 photosynthesis

  REAL(wp), SAVE :: &
      ! intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I))) (µmol CO2 / mol quanta)
      ! accounting for the absorptance of leaves (0.85), the maximum quantum yield of PS II (0.7)
      ! the fraction of total light that reaches PS II (0.5) and the NADPH limited RubGeneration (jointly leading to 0.066)
      alpha_i     ,     & !< intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I)) (µmol CO2 / mol quanta)

      ! Partial Pressure of O2 in kPa
      pO2             ,     & !< Partial Pressure of O2 in kPa
      ! initial guess of the leaf Ci to Ca ratio
      CiCa_default_C3 ,     & !< default Ci:Ca for C3 plants
      CiCa_default_C4 ,     & !< default Ci:Ca for C4 plants
      ! maximum number of iteration in photosynthesis calculation (should probably be in mo_X_constants)
      ps_it_max       ,     & !< maximum iteration of PS calculation
      ci_max          ,     & !< saturating Ci in Pa C4 plants
      ! ratio of the molecular diffusion constants of water and CO2 in air (Nobel, Plant Phys. 2005, p 380)
      Dwv2co2_air     ,     & !< ratio of diffusion coefficient for H2O and CO2 in air 
      ! ratio of the diffusion constant of water and CO2 in air, accounting for turbulent transfer (Nobel, Plant Phys. 2005, p 380) 
      Dwv2co2_turb            !< ratio of diffusion coefficient for H2O and CO2 in turbulent air 

  ! parameters describing the isotopic discrimination of PS (fractionation due to photosynthesis), derived from Drake 2014, Radiocarbon, 56, 29-38
  REAL(wp), SAVE :: &     
      discr_ps_a_C13  ,     &  !< discrimination of C13 due to stomatal diffusion
      discr_ps_b_C13  ,     &  !< discrimination of C13 due to Rubisco
      discr_ps_c_C13  ,     &  !< discrimination of C13 due to PEP C 
      discr_ps_a_C14  ,     &  !< discrimination of C14 due to stomatal diffusion
      discr_ps_b_C14  ,     &  !< discrimination of C14 due to Rubisco
      discr_ps_c_C14  ,     &  !< discrimination of C14 due to PEP C
      discr_ps_phi             !< leakage rate of bundle sheath cells; a typical value cf Hatch et al. 1995, Plant Phys. 108, 173-181

  ! parameters for canopy profile and canopy light extinction (parameters describing the canopy profile of N)
  REAL(wp), SAVE :: &     
      ka             ,   &  !< Extinction coefficient for PAR on chlorophyll
      kn                    !< extinction coefficient to describe decline of N within the canopy

  REAL(wp), SAVE :: &
      soa_b         ,   &  !< parameter in state of acclimation calculation used for evergreen forests
      soa_t_s              !< parameter in state of acclimation calculation used  for evergreen forests

  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
      INTEGER, PARAMETER ::    &
      ic3phot       = 1  ,     &
      ic4phot       = 2  

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_assimi_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the assimilation process
  !! called by jsbach_setup_models()
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_assimi_constants

    USE mo_jsb_physical_constants,  ONLY: molar_mass_N, &   !< parameter from mo_jsb_phys_constants
                                          Dwv, &            !< parameter from mo_jsb_phys_constants - diffusion coefficient for water vapour in air
                                          Dco2              !< parameter from mo_jsb_phys_constants - diffusion coefficient for CO2 in air

    IMPLICIT NONE
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_assimi_constants'

    ! Parameters describing the PS ~ An relationship (describing photosynthetic activity )
      jmax2n        = 17.6_wp/4.0_wp
      vcmax2n       = 1.8_wp       
      pepc2n        = 98777.97_wp 
      chl2n         = 1000._wp/39.8_wp 
      k1_fn_struc   = 71.4_wp * molar_mass_N / 1000000._wp 

    ! Parameters used in the photosynthesis routine (temperature sensitivies of PS processes)
      E0kc            = 38.05_wp    ! Bernacchi 2001 
      E1kc            = 79430.0_wp  ! Bernacchi 2001
      E0ko            = 20.3_wp     ! Bernacchi 2001
      E1ko            = 36380.0_wp  ! Bernacchi 2001
      E0pcp           = 19.02_wp    ! Bernacchi 2001
      E1pcp           = 37830.0_wp  ! Bernacchi 2001
      E0v             = 26.35_wp    ! Bernacchi 2001
      E1v             = 65330.0_wp  ! Bernacchi 2001
      t_jmax_offset   = 17.0_wp     ! Friend 2009
      t_jmax_slope    = 0.35_wp     ! Friend 2009
      t_jmax_opt_min  = t_jmax_offset 
      t_jmax_opt_max  = 38.0_wp     ! Friend 2009
      Tref_pepc       = 25._wp      ! Friend 2009
      Tbase_pepc      = 10._wp      ! Friend 2009
      alpha_i         = 0.85_wp * 0.066_wp  
      pO2             = 20.9_wp 
      CiCa_default_C3 = 0.7_wp
      CiCa_default_C4 = 0.4_wp
      ps_it_max       = 10._wp 
      ci_max          = 7800.0_wp  
      Dwv2co2_air     = Dwv / Dco2  
      Dwv2co2_turb    = ( Dwv / Dco2  ) ** (2._wp/3._wp) 

    ! parameters describing the isotopic discrimination of PS (fractionation due to photosynthesis), derived from Drake 2014, Radiocarbon, 56, 29-38
      discr_ps_a_C13 = 4.4_wp   
      discr_ps_b_C13 = 27.0_wp 
      discr_ps_c_C13 = 5.7_wp 
      discr_ps_a_C14 = 1.97_wp * discr_ps_a_C13
      discr_ps_b_C14 = 1.89_wp * discr_ps_b_C13
      discr_ps_c_C14 = 1.89_wp * discr_ps_c_C13
      discr_ps_phi   = 0.16_wp 

    ! parameters for canopy profile and canopy light extinction (parameters describing the canopy profile of N)
      ka            = 0.005_wp 
      kn            = 0.11_wp

      soa_b         = -0.5_wp
      soa_t_s       = 5.0_wp

  END SUBROUTINE init_assimi_constants

#endif
END MODULE mo_assimi_constants
