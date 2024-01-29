! @brief variables for marine aggregate settling
!
! Definition of variables and allocation of memory
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

MODULE mo_memory_agg

  USE mo_kind, ONLY : wp
  USE mo_control_bgc, ONLY: bgc_nproma, bgc_zlevs

  IMPLICIT NONE

  PUBLIC

  REAL(wp) :: AJ1, AJ2, AJ3, BJ1, BJ2, BJ3 ! constants for CD 

! for fix b (numbers distribution slope of aggregates)
  ! primary particle diameter for POM & PIM species involved in parametrized aggregation (m) 
  REAL(wp) :: dp_dust ! primary particle diameter dust
  REAL(wp) :: dp_det  ! primary particle diameter detritus
  REAL(wp) :: dp_calc ! primary particle diameter calc
  REAL(wp) :: dp_opal ! primary particle diameter opal
  REAL(wp) :: stickiness_tep  ! stickiness of TEP (related to opal frustules)
  REAL(wp) :: stickiness_det  ! normal detritus stickiness
  REAL(wp) :: stickiness_opal ! stickiness of opal (without TEP - just normal coating)
  REAL(wp) :: stickiness_calc ! stickiness of calc particles (coated with organics)
  REAL(wp) :: stickiness_dust ! stickiness of dust particles (coated with organics)
  REAL(wp) :: agg_df_max      ! maximum fractal dimension of aggregates (~2.5)
  REAL(wp) :: agg_df_min      ! minimum fractal dimension of aggregates (~1.2 - 1.6)
  REAL(wp) :: rho_tep         ! density of TEP particles
  REAL(wp) :: agg_Re_crit ! critical particle Reynolds number for nr-distribution limiting
  ! organic detritus density (alternative to orgdens to avoid negative ws)
  REAL(wp) :: agg_org_dens
  REAL(wp) :: det_mol2mass ! mol detritus P/m^3 to kg POM /m^3 (according to stoichiometry)

  INTEGER, PARAMETER :: &
       kavdp               =  1, &
       kavrhop             =  2, &
       kdfagg              =  3, &
       ksticka             =  4, &
       kLmaxagg            =  5, &
       kavrhof             =  6, &
!       kwsagg              =  7, &
!       kbagg               =  8, &
!       kstickf             =  9, &
!       kdynvis             = 10, &
!       kavdc               = 11, &
!       kavpor              = 12, &
       naggdiag            = 6

CONTAINS


END MODULE mo_memory_agg
