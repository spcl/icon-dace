!> QUINCY assimilation constants
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
!>#### declare and define assimilation constants
!>
MODULE mo_q_assimi_constants
#ifndef __NO_QUINCY__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_physical_constants,  ONLY: Dwv, &            !< diffusion coefficient for water vapour in air
                                        Dco2              !< diffusion coefficient for CO2 in air

  IMPLICIT NONE
  PUBLIC

  REAL(wp), SAVE :: &
    ! ratio of the molecular diffusion constants of water and CO2 in air (Nobel, Plant Phys. 2005, p 380)
    & Dwv2co2_air  = Dwv / Dco2, &                      !< ratio of diffusion coefficient for H2O and CO2 in air
    ! ratio of the diffusion constant of water and CO2 in air, accounting for turbulent transfer (Nobel, Plant Phys. 2005, p 380)
    ! ACCWA (nvhpc 21 on daint; fractional exponent doesn't work in module data initialization)
#if (__NVCOMPILER_MAJOR__ <= 21)
    & Dwv2co2_turb = EXP((2.0_wp / 3.0_wp) * LOG(Dwv / Dco2)) !< ratio of diffusion coefficient for H2O and CO2 in turbulent air
#else
    & Dwv2co2_turb = (Dwv / Dco2) ** (2.0_wp / 3.0_wp)        !< ratio of diffusion coefficient for H2O and CO2 in turbulent air
#endif

  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
  ! TODO re-write as ENUM
  INTEGER, PARAMETER ::    &
    & ic3phot   = 1, &
    & ic4phot   = 2

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_assimi_constants'

#endif
END MODULE mo_q_assimi_constants
