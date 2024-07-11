! @brief configuration setup for turbulent diffusion (turbdiff)
!
! configuration setup for turbulent diffusion
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

MODULE mo_turbdiff_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_NTRACER, MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC



  !!--------------------------------------------------------------------------
  !! Basic configuration setup for turbulent diffusion (turbdiff)
  !!--------------------------------------------------------------------------
  TYPE :: t_turbdiff_config

    ! namelist variables

    REAL(wp), DIMENSION(:), POINTER :: &
      &  impl_weight    ! implicit weights for tridiagonal solver

    INTEGER :: &   ! mode of surface-atmosphere transfer
      &  imode_tran 
    INTEGER :: &   ! mode of cloud representation in transfer parametr.
      &  icldm_tran
    INTEGER :: &   ! mode of turbulent diffusion parametrization
      &  imode_turb
    INTEGER :: &   ! mode of cloud representation in turbulence parametr.
      &  icldm_turb
    INTEGER :: &   ! type of shear production for TKE
      &  itype_sher
    INTEGER :: &   ! mode of vertical smoothing of TKE source terms
      &  imode_frcsmot
    INTEGER :: &   ! mode of the separated horizontal shear mode 
      &  imode_shshear
    INTEGER :: &   ! mode of the SSO-turbulence coupling
      &  imode_tkesso

    LOGICAL :: &   ! calculation SSO-wake turbulence production for TKE
      &  ltkesso
    LOGICAL :: &   ! consider convective buoyancy production for TKE
      &  ltkecon
    LOGICAL :: &   ! calculation separ. horiz. shear production for TKE
      &  ltkeshs
    LOGICAL :: &   ! explicit corrections of the implicit calculated turbul. diff.
      &  lexpcor
    LOGICAL :: &   ! consideration of thermal TKE-sources in the enthalpy budget
      &  ltmpcor
    LOGICAL :: &   ! using the profile values of the lowest main level instead of
      &  lprfcor   ! the mean value of the lowest layer for surface flux calulations

    LOGICAL :: &   ! nonlocal calculation of vertical gradients used for turbul. diff.
      &  lnonloc 
    LOGICAL :: &   ! free-slip lower boundary condition (use for idealized runs only!)
      &  lfreeslip 
    LOGICAL :: &   ! consideration of fluctuations of the heat capacity of air
      &  lcpfluc
    LOGICAL :: &   ! lower flux condition for vertical diffusion calculation
      &  lsflcnd

    LOGICAL :: &   ! TRUE: horizontally homogeneous roughness length 
      &  lconst_z0 ! (for idealized testcases)

    REAL(wp):: &   ! turbulent master length scale 
      &  tur_len   ! (devided by roughness length)
    REAL(wp):: &   ! effectiv length scale of surface patterns
      &  pat_len   !
    REAL(wp):: &   ! scaling factor for stability correction of 'tur_len'
      &  a_stab    !
    REAL(wp):: &   ! mode of estimating the Charnock-Parameter
      & imode_charpar
    REAL(wp):: &   ! lower limit of velocity-dependent Charnock-parameter
      &  alpha0    !
    REAL(wp):: &   ! upper limit of velocity-dependent Charnock-parameter
      &  alpha0_max !
    REAL(wp):: &   ! additive ensemble perturbation of Charnock-parameter
      &  alpha0_pert !
    REAL(wp):: &   ! scaling factor for molecular roughness length of ocean waves
      &  alpha1    !
    REAL(wp):: &   ! minimal diffusion coefficient for scalars (heat)
      &  tkhmin    !
    REAL(wp):: &   ! minimal diffusion coefficient for momentum
      &  tkmmin    !
    REAL(wp):: &   ! enhanced minimal diffusion coefficient for scalars (heat) in the stratosphere
      &  tkhmin_strat    !
    REAL(wp):: &   ! enhanced minimal diffusion coefficient for momentum in the stratosphere
      &  tkmmin_strat    !
    REAL(wp):: &   ! constant diffusion coefficient for TKE
      &  c_diff    !
    REAL(wp):: &   ! scaling factor of laminar layer for scalars (heat)
      &  rlam_heat !
    REAL(wp):: &   ! scaling factor of laminar layer for momentum
      &  rlam_mom  !
    REAL(wp):: &   ! scaling correction factor for laminar layers of water surfaces
      &  rat_sea   ! 
    REAL(wp):: &   ! ratio of laminar scaling factors for vapour and heat
      &  rat_lam   !
    REAL(wp):: &   ! scaling correction factor for laminar layers of glacier surfaces
      &  rat_glac  ! 
    REAL(wp):: &   ! time smoothing factor for TKE
      &  tkesmot   ! 
    REAL(wp):: &   ! vertical smoothing factor of TKE forcing terms
      &  frcsmot   ! 
    REAL(wp):: &   ! implicit weight near the surface (maximal value)
      &  impl_s    ! 
    REAL(wp):: &   ! implicit weight near top of the atmosphere (maximal value)
      &  impl_t    ! 
    REAL(wp):: &   ! Length scale factor for the separated horizontal shear mode
      &  a_hshr    ! 
    REAL(wp):: &   ! Parameter for pdf width in turbulent cloud cover scheme
      &  q_crit    ! 

    REAL(wp):: &   ! horizontally homogeneous roughness length 
      &  const_z0  ! (for idealized testcases)

    LOGICAL :: &   ! turbulent diffusion of cloud ice QI
      &  ldiff_qi  ! .TRUE.: ON

    LOGICAL :: &   ! turbulent diffusion of snow QS
      &  ldiff_qs  ! .FALSE.: OFF

    !
    ! Switches controlling other physical parameterizations:
    !
    INTEGER :: &  ! type of water cloud diagnosis
      &  itype_wcld

    !
    ! derived variables
    !

                   
  END TYPE t_turbdiff_config

  !>
  !!
  TYPE(t_turbdiff_config), TARGET :: turbdiff_config(0:max_dom)


!CONTAINS


END MODULE mo_turbdiff_config
