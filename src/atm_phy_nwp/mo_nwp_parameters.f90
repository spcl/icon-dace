!
!  Declares parameters computed during the initialization of the physics
!  parameterizations that have to be domain-dependent
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

MODULE mo_nwp_parameters

!-------------------------------------------------------------------------

  USE mo_kind,               ONLY: wp


  IMPLICIT NONE

  PRIVATE


  TYPE t_phy_params
    !
    ! Parameters which are only computed if convection is switched on
    !
    ! Level parameters for convection scheme
    INTEGER  :: kcon1, kcon2
    ! resolution-dependent parameters for convection scheme
    REAL(wp) :: tau, mfcfl, tau0
    ! relative humidity below which sub-cloud evaporation of rain starts over land and water, respectively
    REAL(wp) :: rhebc_land, rhebc_ocean, rhebc_land_trop, rhebc_ocean_trop
    ! 'excess values' of temperature and QV used for convection triggering (test parcel ascent)
    REAL(wp) :: texc, qexc
    ! fractional area covered by convective precipitation
    REAL(wp) :: rcucov, rcucov_trop
    ! tuning coefficient for organized entrainment of deep convection
    REAL(wp) :: entrorg
    ! detrainment rate for deep convection updrafts and downdrafts
    REAL(wp) :: detrpen, entrdd
    ! entrainment constants for test parcel ascent
    REAL(wp) :: entstpc1, entstpc2
    ! coefficient for conversion of cloud water into precipitation
    REAL(wp) :: rprcon
    ! maximum allowed depth of shallow convection (hPa)
    REAL(wp) :: rdepths
    ! critical stability threshold for stratocumulus (K)
    REAL(wp) :: eiscrit    
    ! switches for activation of shallow, midlevel and deep convection
    LOGICAL :: lmfscv, lmfmid, lmfpen
    ! switch for detrainment of rain and snow to gridscale scheme
    LOGICAL :: lmfdsnow
    ! switch for grayzone tuning for deep convection
    LOGICAL :: lgrayzone_deepconv
    ! switches on explicit stochastic shallow convection
    LOGICAL :: lstoch_expl
    ! switches on SDE stochastic shallow convection    
    LOGICAL :: lstoch_sde
    ! switches on SDE stochastic deep convection    
    LOGICAL :: lstoch_deep
    ! use 650hPa vertical velocity to switch off conv param at points with rising motion     
    LOGICAL :: lvvcouple
    ! use 650hPa vertical velocity to distinguish shallow vs deep convection
    LOGICAL :: lvv_shallow_deep
    !
    ! Parameters which are only computed if Gravity wave drag scheme is switched on
    !
    ! launch level for GWD scheme
    INTEGER  :: klaunch
    !
    ! Parameters which are only computed if Sub-grid Scale Orography (SSO) scheme is switched on
    !
    INTEGER  :: ngwdlim, ngwdtop, nktopg
    REAL(wp) :: gkwake, gkdrag, gfrcrit, grcrit, minsso, minsso_gwd, blockred, gkdrag_enh, grcrit_enh
    !
    ! Parameters which are always computed
    !
    ! characteristic horizontal length scale (grid-scale) for 
    ! turbulence scheme and convection scheme
    REAL(wp) :: mean_charlen
    ! level index corresponding to the HAG of the 60hPa level (identical to kcon2, if computed)
    INTEGER :: k060
  END TYPE t_phy_params


  PUBLIC :: t_phy_params


END MODULE mo_nwp_parameters
