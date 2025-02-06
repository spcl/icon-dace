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

! Definition of variables and allocation of memory

MODULE mo_memory_bgc

  USE mo_kind, ONLY   : wp
  USE mo_param1_bgc, ONLY: natm
  USE mo_control_bgc, ONLY: rmasko, bgc_nproma, bgc_zlevs
  USE mo_hamocc_nml, ONLY: l_cpl_co2 
 

  IMPLICIT NONE

  PUBLIC

  INTEGER :: n90depth,n1000depth,n2000depth

  REAL(wp), DIMENSION(6):: dmsp

  REAL(wp) :: calcinp, orginp, silinp
  REAL(wp) :: roc13
  REAL(wp) :: ozkoa, totalarea
  REAL(wp) :: globalmean_co2, globalmean_o2, globalmean_n2
  REAL(wp) :: atmacon, atmacmol
  REAL(wp) :: ems_per_step
  REAL(wp) :: totarea
  REAL(wp) :: phytomi, grami, grazra, rrrcl,docmin
  !$ACC DECLARE CREATE(rrrcl)
  REAL(wp) :: remido, dyphy, zinges, epsher, spemor, gammap, gammaz, ecan
  REAL(wp) :: ro2ut, rcar, rnit, rnoi, rnit23, rnit13, ropal, rn2, p2gtc
  REAL(wp) :: bkphy, bkzoo, bkopal, bifr13, bifr14, plafr13, plafr14
  REAL(wp) :: dremdoc, dremn2o
  REAL(wp) :: thresh_aerob ! O2 threshold for aerob POC remineralization
  REAL(wp) :: thresh_sred ! O2 threshold for sulfate reduction
  REAL(wp) :: thresh_o2, prodn2o 
  REAL(wp) :: mc_fac, mc_eu_speed, mc_depth
  REAL(wp) :: n2_fixation
  REAL(wp) :: sulfate_reduction
  REAL(wp) :: wcya  ! daily boyancy speed of cyanos
  REAL(wp) :: gutc
  REAL(wp) :: psedi, csedi, ssedi
  REAL(wp) :: perc_diron, riron, fesoly, relaxfe, sinkspeed_dust,bolaymin
  REAL(wp) :: pi_alpha, fPAR, bkh2sox, rh2sox
  REAL(wp) :: nitdem, n2prod, ro2nitri
  REAL(wp) :: dremn3o, ro2bal
  REAL(wp) :: dustd1, dustd2, dustd3, dustsink

  REAL(wp) :: pi_alpha_cya,cya_growth_max,Topt_cya,T1_cya,T2_cya ! (namelist parameter)
  REAL(wp) :: bkcya_N, doccya_fac           ! (namelist parameter)
  REAL(wp) :: buoyancyspeed_cya      
  REAL(wp) :: ralk, cyamin, ro2ut_cya

  ! Extended Nitrogen cycle variables
  REAL(wp) :: bkno3_cya, bknh4_cya, ro2ammo, bkno3, bknh4
  REAL(wp) :: rmm, kg_denom, bkpo4, o2thresh, o2den_lim
  REAL(wp) :: no2denit, anamoxra, bkno2, bkrad, nitriox, nitrira
  REAL(wp) :: bkfe, rno3nh4, rno3no2, rnh4no2, rno2no3, alk_nrn2
  REAL(wp) :: rno2n2


  REAL(wp), PARAMETER :: ppm2con=0.35e-3_wp     !> ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
                                                !> --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2
  REAL(wp), PARAMETER :: contppm=1._wp/ppm2con  !> inverse of conversion factor molC/m**2 to ppm (in AT)

  REAL(wp), PARAMETER :: molw_co2=44.011_wp     !> Total Molecular Weight of CO2 (44.011 g/mol)
 
  REAL(wp), PARAMETER :: molw_dry_air=28.970_wp !> mean pseudo-molecular weight of dry air



CONTAINS

END MODULE mo_memory_bgc
