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

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_assimi_constants'

#endif
END MODULE mo_assimi_constants
