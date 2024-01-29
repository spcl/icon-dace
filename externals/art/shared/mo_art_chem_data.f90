!
! mo_art_chem_data
! This module provides data storage structures and constants for chemistry.
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

MODULE mo_art_chem_data
! ICON
  USE mo_kind,                          ONLY: wp
! ART
  USE mo_art_OH_chem_types,             ONLY: t_art_OH_chem_meta
  USE mo_art_chem_init_types,           ONLY: t_chem_init_state
#ifdef __ART_GPL
  USE mo_art_psc_types,                 ONLY: t_art_psc
  USE mo_art_mecicon_data,              ONLY: t_art_mecicon
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: nphot
  PUBLIC :: t_art_chem
  PUBLIC :: t_art_chem_indices
  PUBLIC :: t_art_linoz
  PUBLIC :: t_art_photolysis
  PUBLIC :: t_art_chem_param

  INTEGER, PARAMETER :: &
    &  nphot = 72           !< number of photolysis rates


  TYPE t_art_photolysis
    REAL(wp), POINTER :: &
      &  rate(:,:,:,:)             !< photolysis rates in 1/s
    REAL(wp), ALLOCATABLE ::             &
      &  heating_rates(:,:,:,:),         &  !< heating rates in K/day
      &  input_fastjx_clf(:,:,:),        &  !< ozone vmr in ppmv
      &  input_fastjx_rh(:,:,:),         &  !< relative humidity in [0-1]
      &  input_fastjx_lwp(:,:,:),        &  !< liquid water path in g/m2
      &  input_fastjx_iwp(:,:,:),        &  !< liquid water path in g/m2
      &  input_fastjx_iwc(:,:,:),        &  !< ice water content in g/m3
      &  input_fastjx_lwc(:,:,:),        &  !< liquid water path in g/m3
      &  input_fastjx_reice(:,:,:),      &  !< effective radius [microns] ice clouds
      &  input_fastjx_reliq(:,:,:),      &  !< effective radius [microns] liq clouds
      &  fjx_zkap(:,:),                  &  !< breath parameter for scaling effective radius
      &  fjx_zland(:,:),                 &  !< land sea mask (0,1)
      &  fjx_zglac(:,:),                 &  !< fraction of land covered by glaciers
      &  output_fastjx_photo_new(:,:,:,:)   !< photolysis rates in 1/s

    REAL(wp),  ALLOCATABLE  ::  &
      &  PPP(:),                &
      &  ZZZ(:),                &
      &  TTT(:),                &
      &  DDD(:),                &
      &  RRR(:),                &
      &  OOO(:),                &
      &  LWP(:),                &
      &  IWP(:),                &
      &  REFFL(:),              &
      &  REFFI(:),              &
      &  CLF(:),                &
      &  AERSP(:,:),            &
      &  VALJXX(:,:),           &
      &  SKPERD(:,:),           &
      &  fjx_cdnc(:,:,:)          !< cloud droplet number concentration [1/m**3]

    INTEGER, ALLOCATABLE  ::  &
      &  CLDIW(:),            &
      &  NDXAER(:,:)

  END TYPE t_art_photolysis

  TYPE t_art_linoz
    REAL(wp), ALLOCATABLE :: &
      &  tparm(:,:,:,:)
    REAL(wp), ALLOCATABLE :: & 
      &  linoz_tab1(:,:), &
      &  linoz_tab2(:,:), &
      &  linoz_tab3(:,:), &
      &  linoz_tab4(:,:), &
      &  linoz_tab5(:,:), &
      &  linoz_tab6(:,:), &
      &  linoz_tab7(:,:)
    LOGICAL :: &
      &  is_init = .FALSE.

  END TYPE t_art_linoz

  TYPE t_art_chem_param
    REAL(wp), ALLOCATABLE ::  &
      &  o2_column(:,:,:)         !< Units: #/cm2 J
    INTEGER ::                 &
      &  number_param_tracers     !< number of tracers of class t_chem_meta_param and its extensions
    LOGICAL ::  &
      &  lregio_tracers = .FALSE. !< are there regional tracers?
    TYPE(t_art_linoz) ::  &
      &  linoz
    TYPE(t_art_OH_chem_meta)  :: &
      &  OH_chem_meta
  END TYPE t_art_chem_param

  TYPE t_art_chem_indices
    ! chemical tracers
    INTEGER ::  &
       &    iTRCHBR3 = 0,     & !< Indices in tracer container
       &    iTRCH2BR2 = 0,    &  
       &    iTRCH4 = 0,       &
       &    iTRC2H6 = 0,      &
       &    iTRC5H8 = 0,      &
       &    iTRC3H8 = 0,      &
       &    iTRCH3COCH3 = 0,  &
       &    iTRCH3CN = 0,     &
       &    iTRCO = 0,        &
       &    iTRCO2 = 0,       &
       &    iTRH2O = 0,       &
       &    iTRH2O_feed = 0,  &
       &    iTRO3 = 0,        &
       &    iTRO3_pas = 0,    &
       &    iTRN2O = 0,       &
       &    iTRNOy = 0,       &
       &    iTRNO2 = 0,       &
       &    iTRN2O5 = 0,      &
       &    iTRHNO3 = 0,      &
       &    iTRNH3 = 0,       &
       &    iTRSO2 = 0,       &
       &    iTROCS = 0,       &
       &    iTRDMS = 0,       &
       &    iTRH2SO4 = 0,     &
       &    iTRHCl = 0,       &
       &    iTRHOCl = 0,      &
       &    iTRClONO2 = 0,    &
       &    iTRHBr = 0,       &
       &    iTRHOBr = 0,      &
       &    iTRBrONO2 = 0,    &
       &    iTRCFCl3 = 0,     &
       &    iTRCF2Cl2 = 0,    &
       &    iTR_cold = 0,     &
       &    iTRAGE = 0,       &
       &    iTR_vortex = 0
    ! regional tracers
    INTEGER ::                 &
      &  iTR_stn  = 0, iTR_stt  = 0,  & !< indices in tracer container
      &  iTR_sts  = 0, iTR_trn  = 0,  & !< indices in tracer container
      &  iTR_trt  = 0, iTR_trs  = 0,  & !< indices in tracer container
      &  iTR_tiln = 0, iTR_tils = 0,  & !< indices in tracer container
      &  iTR_nh   = 0, iTR_sh   = 0,  & !< indices in tracer container
      &  iTR_nin  = 0, iTR_sin  = 0,  & !< indices in tracer container
      &  iTR_ech  = 0, iTR_sea  = 0,  & !< indices in tracer container
      &  iTR_sib  = 0, iTR_eur  = 0,  & !< indices in tracer container
      &  iTR_med  = 0, iTR_naf  = 0,  & !< indices in tracer container
      &  iTR_saf  = 0, iTR_mdg  = 0,  & !< indices in tracer container
      &  iTR_aus  = 0, iTR_nam  = 0,  & !< indices in tracer container
      &  iTR_sam  = 0, iTR_tpo  = 0,  & !< indices in tracer container
      &  iTR_tao  = 0, iTR_tio  = 0,  & !< indices in tracer container
      &  iTR_bgn  = 0, iTR_bgs  = 0,  & !< indices in tracer container
      &  iTR_nest = 0,                & !< indices in tracer container
      &  iTR_art  = 0
  END TYPE t_art_chem_indices

  TYPE t_art_chem
    REAL(wp), ALLOCATABLE :: &
      &  water_tracers(:,:,:,:)    !< water tracers (converted from kg/kg to #/cm3 in reaction
                                   !  interface
    REAL(wp), ALLOCATABLE :: &
      &  CO2_mmr_depos(:)
    REAL(wp), ALLOCATABLE :: &
      &  vmr2Nconc(:,:,:)          !< conversion from mol/mol to cm-3
    TYPE(t_art_photolysis) :: &
      &  photo
    TYPE(t_art_chem_indices) ::  &
      &  indices
    TYPE(t_art_chem_param) :: &
      &  param
    TYPE(t_chem_init_state) :: &
      &  init
#ifdef __ART_GPL
    TYPE(t_art_psc) ::  &
      &  PSC_meta
    TYPE(t_art_mecicon) :: &
      &  mecicon
#endif
    REAL(wp), POINTER :: &
      &  so2_column(:,:,:)           !< SO2 column in DU
  END TYPE t_art_chem


END MODULE mo_art_chem_data
