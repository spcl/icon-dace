!
! mo_art_psc_types
! This module provides datastructures for PSCs
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

MODULE mo_art_psc_types

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_physical_constants,            ONLY: amw
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: t_art_psc

  PUBLIC :: Pa2torr
  PUBLIC :: mol_weight_HNO3
  PUBLIC :: mol_weight_NAT
  PUBLIC :: rho_NAT
  PUBLIC :: pi4_3_rhoNAT
  PUBLIC :: mol_weight_H2SO4



  REAL(wp), PARAMETER :: &
    &  Pa2torr = 1._wp / 133.322_wp, &     !< Pascal to Torr
    &  mol_weight_H2SO4 = 98.078_wp, &     !< g / mol; from NIST
    &  mol_weight_HNO3 = 63.0128_wp        !< g / mol; from NIST
  REAL(wp), PARAMETER   ::  &
    &  mol_weight_NAT  = 3.0_wp*amw + mol_weight_HNO3  !< g / mol
  REAL(wp), PARAMETER   ::   &
    &  rho_NAT     = 1.626e6_wp        !<  crystal mass density of NAT in g/m**3
  REAL(wp), PARAMETER   ::   &
    &  pi4_3_rhoNAT = 4._wp / 3._wp * pi * rho_NAT !< prefactor for calculation of particle mass 
                                                   !  (to be multiplied by radius**3)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_bin_ptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: &
      &  ptr
  END TYPE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_psc
    INTEGER, ALLOCATABLE ::      &
      &  tracer_indices_bins(:)       !< array of the tracer indices of NAT bins (NSB)
    INTEGER  ::                         &
      &  iTRN2O5, iTRClONO2, iTRBrONO2, &  !< indices of the reactive tracers
      &  iTRHBr, iTRHCl, iTRHOCl, iTRHOBr, &
      &  iTRHNO3, iTRH2SO4
    REAL(wp), POINTER      :: &
      &  liqsur(:,:,:),       &  !< total area of liquid aerosols (nproma, nlev, nblks) 
                                 !  (cm2/cm3 air)
      &  cgaml(:,:,:,:),      &  !< uptake coefficients (nproma, nlev, ihs_MAX)  (-)
      &  k_het(:,:,:,:)          !< heterogeneous reaction rates on PSCs (nproma,nlev,nblks,ihs_MAX)
                                 !  (cm3 s-1)
    INTEGER ::     &
      &  kstart,   &
      &  kend,     &
      &  NSB           !< number of size bins
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: &  !< (NSB)
      &  rbin_min, rbin_max, rbin_av,      &  !< minimum, maximum and average radius of the bins [m]
      &  no_density_limit                     !< limit of the number density in the bin (NSB) [cm-3]

    REAL(wp), POINTER, DIMENSION(:,:) :: & !< (nproma, nlev)
      &  H2O_Nconc_g,      &  !< number concentration of gaseous water (#/cm3)
      &  H2O_Nconc_l,      &  !< number concentration of liquid  water (#/cm3)
      &  H2O_Nconc_s          !< number concentration of solid   water (#/cm3)

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: & !< (nproma, nlev)
      &  HNO3_Nconc_g         !< number concentration of gaseous HNO3 (#/cm3)

    REAL(wp), POINTER, DIMENSION(:,:,:) :: & !< (nproma,nlev,nblks)
      &  ice_vmr_Marti,  &  !< volume mixing ratio of solid water by Marti param. (mol/mol)
      &  HNO3_Nconc_l,   &  !< number concentration of liquid  HNO3 (#/cm3)
      &  HNO3_Nconc_s       !< number concentration of solid   HNO3 (#/cm3)
                            !  (equals number concentration of NAT since 1 molecule NAT 
                            !   includes 1 HNO3 molecule)

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &  !< (nproma,nlev,nblks)
      &  radius_STS,     &  !< radius of STS particles (m)
      &  dens_ice,       &  !< number density of ice particles (m-3)
      &  radius_ice         !< radius of ice particles (m)

    REAL(wp), POINTER, DIMENSION(:,:,:,:) :: &  !< (nproma,nlev,nblks,NSB)
      &  v_sed_NAT,         &  !< sedimentation velocity of NAT particles (m s-1)
      &  v_sed_NAT_out,     &  !< sedimentation velocity of NAT particles for output (m s-1)
      &  dens_NAT,          &  !< number density of NAT molecules (1 / m3)
      &  radius_NAT            !< mean NAT radius (m)

    REAL(wp), POINTER, DIMENSION(:,:,:) :: &  !< (nproma,nblks,NSB)
      &  NAT_sedi_rel_diff     !< relative mass difference of NAT before and after sedimentation (-)

    REAL(wp)             ::             &
      &  total_max_density = 2.3e-4_wp, & !< "global" maximum number density
                                          !  integrated over all size bins
                                          !  (cm-3)
      &  NatFormThreshold = -3._wp,     & !< temperature threshold below which NAT should
                                          !  form initially (K)
      &  r_min = 1.e-7_wp                 !< min. radius which is used for thermodynamic NAT param. 
                                          !  to get number density in this case (m)

    REAL(wp)             ::    &
      &  khet_max = 1.e-13_wp, &  !< boundaries for heterogeneous reaction rates (cm3 s-1 molec-1)
      &  khet_min = 0.0_wp        !  to avoid unrealistically high values

    REAL(wp), ALLOCATABLE ::    & !< Variables for chemistry
      &  v_th(:,:),             & !< mean molecular speed without molar weight (m/s*sqrt(g/mol))
      &  mean_free_path(:,:),   & !< mean free path of the molecules in air (m)
      &  Nconc_higher(:,:)        !< number concentration of the reactant with higher concentration,
                                  !  i.e. H2O, HCl or HBr (molec / cm3)

  
    INTEGER, ALLOCATABLE ::    &
       &  iexception(:,:)         !< flag for computation of STS particles (see Eckstein, 2013)
    REAL (KIND=wp), ALLOCATABLE, DIMENSION(:,:) ::  &
       &    sqrtt,       &     !< square root of temperature
       &    pw,          &     !< pw     = water partial pressure [atm]
       &    logpw,       &     !< logarithm of pw
       &    ns,          &     !< ns     = total moles sulphate/m3 air 
                               !  (assumed to be totally in liquid phase)
       &    pn0,         &     !< pn0    = initial hno3 partial pressure [atm]
       &    phocl0,      &     !< phocl0 = initial hocl partial pressure [atm]
       &    phobr0,      &     !< phobr0 = initial hobr partial pressure [atm]
       &    phcl0,       &     !< phcl0  = initial hcl partial pressure [atm]
       &    phbr0,       &     !< phbr0  = initial hbr partial pressure [atm]
       &    parthno3,    &     !< parthno3 = fraction of hno3 remaining in gas phase after 
                               !             partitioning into liquid (1 = no removal to aerosols)
       &    parthcl,     &     !< parthcl = fraction of hcl remaining in gas phase 
                               !   after partitioning into liquid (1 = no removal to aerosols)
       &    parthbr,     &     !< parthbr = fraction of hbr remaining in gas phase 
                               !   after partitioning into liquid (1 = no removal to aerosols)
       &    hhocl,       &     !< hhocl  = effective hocl henry's law constant 
                               !           in hno3-h2so4-h2o solution [mol / kg water / atm]
       &    hhobr,       &     !< hhobr  = effective hobr henry's law constant 
                               !           in hno3-h2so4-h2o solution [mol / kg water / atm]
       &    hhbr,        &     !< hhbr   = effective hbr henry's law constant 
                               !           in hno3-h2so4-h2o solution [mol / kg water / atm]
       &    mn,          &     !< mn     = concentration of hno3 [mole / kg water]      
       &    ms,          &     !< ms     = concentration of h2so4 [mole / kg water]
       &    wn,          &     !< wn     = weight fraction of hno3 in the aerosol
       &    ws,          &     !< ws     = weight fraction of h2so4 in the aerosol
       &    whocl,       &     !< whocl  = weight fraction of hocl in the aerosol
       &    whobr,       &     !< whobr  = weight fraction of hobr in the aerosol
       &    wcl,         &     !< wcl    = weight fraction of hcl in the aerosol
       &    tice,        &     !< tice   = ice frost point temperature [k]
       &    density,     &     !< density of ternary solution in g/cm3
       &    diff,        &     !< diff: diffusivity of hocl in cm2 s-1
       &    wnen,        &     !< sum of molalities of substances in ternary solution multiplied 
                               !  by molar weights
       &    internal_temp      !< air temperature with lower limit 185 K or tice - 3 K

  END TYPE t_art_psc
  
END MODULE mo_art_psc_types
