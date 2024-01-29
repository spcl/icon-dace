!
! Setup for 2-moment cloud microphysics scheme
!
! default configuration setup for synthetic radar data on the model grid
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

MODULE mo_2mom_mcrph_config_default
  
  USE mo_kind, ONLY: wp
  USE mo_2mom_mcrph_config, ONLY: t_cfg_2mom

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Namelist parameters
  !--------------------------------------------------------------------------

  !.. Type instance to hold the defaults for the config params:
  TYPE(t_cfg_2mom), PARAMETER :: cfg_2mom_default = t_cfg_2mom ( &

       !-----------------------
       ! .. General parameters:
       !-----------------------
       &            1, &          ! i2mom_solver: 0) explicit (1) semi-implicit solver

       !----------------------------------
       ! .. Parameters for cloud droplets:
       !----------------------------------
       &           -1, &          ! ccn_type: 6,7,8,9; if not set by namelist, the ccn_type_gscp4 or ccn_type_gscp5 will win

       ! .. Parameters for rain:
       !------------------------
       &            .FALSE.,  &   ! luse_mu_Dm_rain
       &          -999.99_wp, &   ! nu for rain, N(x) = N0 * D^nu * exp(-lambda*x^mu)
       &            6.0_wp, &     ! rain_cmu0
       &            30.0_wp, &    ! rain_cmu1
       &            1.1e-3_wp, &  ! rain_cmu3 = D_br
       &            1.0_wp, &     ! rain_cmu4

       !-----------------------------
       ! .. Parameters for cloud ice:
       !-----------------------------

       !------------------------
       ! .. Parameters for snow:
       !------------------------
       &            -999.99_wp, & ! cap_snow capacitance for snow deposition/sublimation - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! vsedi_max_s max fallspeed limit for snow - in this case the background value in mo_2mom_mcrph_main.f90 will win

       !---------------------------
       ! .. Parameters for graupel:
       !---------------------------
       &            -999.99_wp, & ! nu for graupel, N(x) = N0 * D^nu * exp(-lambda*x^mu) - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! mu for graupel, N(x) = N0 * D^nu * exp(-lambda*x^mu) - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! ageo for graupel, D = ageo*x^bgeo - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! bgeo for graupel, D = ageo*x^bgeo - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! avel for graupel, v = avel*x^bvel - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! bvel for graupel, v = avel*x^bvel - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            270.16_wp, &  ! Tmax_gr_rime
       &            1.0_wp,   &   ! melt_g_tune_fac: factor multiplying melting rate of graupel 

       !------------------------
       ! .. Parameters for hail:
       !------------------------
       &            -999.99_wp, & ! nu for hail, N(x) = N0 * D^nu * exp(-lambda*x^mu) - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! mu for hail, N(x) = N0 * D^nu * exp(-lambda*x^mu) - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! ageo for hail, D = ageo*x^bgeo - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! bgeo for hail, D = ageo*x^bgeo - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! avel for hail, v = avel*x^bvel - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            -999.99_wp, & ! bvel for hail, v = avel*x^bvel - in this case the background value in mo_2mom_mcrph_main.f90 will win
       &            1.0_wp, &     ! melt_h_tune_fac: factor multiplying melting rate of hail

       !------------------------------------------
       ! .. Parameters for conversions/collisions:
       !------------------------------------------
       &            .FALSE.,   &   ! lturb_enhc: Turbulent enhancement of collisons
       &            300.0_wp, &   ! lturb_len:  Turbulent lenght scale (dummy, later overtaken by the TKE scheme, NOT a namelist parameter)
       &            10,        &  ! iice_stick: sticking efficiency of cloud ice
       &            5,        &   ! isnow_stick: sticking efficiency of snow
       &            5,        &   ! iparti_stick: sticking efficiency of frozen inter-categorical collisions
       &            0.01_wp,  &   ! alpha_spacefilling
       &            75.0e-6_wp, & ! D_conv_ii: D-threshold for conversion to snow ice_selfcollection: newly created snowflakes have at least this mean mass diameter
       &            0.50e-3_wp, & ! D_rainfrz_ig
       &            1.25e-3_wp, & ! D_rainfrz_gh
       &            0.10_wp,  &   ! Collision efficiency for graupel autoconversion (dry graupel) 
       &            0.40_wp,  &   ! Collision efficiency for graupel autoconversion (wet graupel)
       &            270.16_wp, &  ! Temperature threshold for switching to wet graupel autoconversion
       &            1 &           ! iicephase: (0) warm-phase 2M (1) mixed-phase 2M
       &            )

END MODULE mo_2mom_mcrph_config_default
