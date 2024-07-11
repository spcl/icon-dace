!> QUINCY routines vertical transport in soil for JSM
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
!>#### multiple routines for the vertical transport along the soil column used with JSM (jena soil model)
!>
MODULE mo_q_sb_jsm_transport
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: one_day, dtime, eps8, eps4
  USE mo_lnd_bgcm_idx

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calc_particle_fluxrate, &
    &       calc_bioturbation_rate, &
    &       calc_bioturbation_transport_wrapper_jsm, &
    &       calc_liquid_phase_transport_wrapper , &
    &       calc_particle_transport_wrapper, &
    &       calc_bioturbation_transport, &
    &       calc_liquid_phase_transport

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_jsm_processes'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_jsm
  !
  !-----------------------------------------------------------------------------------------------------
  !> Interface to the liquid transport routine, owing to the lack of a suitable operator overload
  !!
  !! Input: leachable pools percolation rate, and soil depth
  !!
  !! Output: vertical transport rate and lateral loss (mol/m3/timestep)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_liquid_phase_transport_wrapper( &
    & nc, nsoil_sb, num_sl_above_bedrock, soil_depth_sl, &                            ! in
    & rtm_gasdiffusion_act, &
    & rmm_gasdiffusion_act, percolation_sl, frac_w_lat_loss_sl, &
    & nh4_solute, no3_solute, po4_solute, nh4_n15_solute, no3_n15_solute, &
    & noy, noy_n15, n2o, n2o_n15, n2, n2_n15, &
    & pool_mt_dom, &                                                                  ! in
    & transport_nh4_solute, transport_no3_solute, &                                   ! inout
    & transport_po4_solute, transport_nh4_n15_solute, transport_no3_n15_solute, &
    & leaching_dom_carbon, leaching_dom_nitrogen, leaching_dom_phosphorus, &
    & leaching_dom_carbon13, leaching_dom_carbon14, leaching_dom_nitrogen15, &
    & leaching_nh4_solute, leaching_no3_solute, leaching_po4_solute, &
    & leaching_nh4_n15_solute, leaching_no3_n15_solute, &
    & transport_mt_dom, &                                                             ! inout
    & lateral_loss_dom_carbon_sl, lateral_loss_dom_nitrogen_sl, &                     ! out
    & lateral_loss_dom_phosphorus_sl, &
    & lateral_loss_dom_carbon13_sl, lateral_loss_dom_carbon14_sl, &
    & lateral_loss_dom_nitrogen15_sl, &
    & lateral_loss_nh4_solute_sl, lateral_loss_no3_solute_sl, &
    & lateral_loss_po4_solute_sl, &
    & lateral_loss_nh4_n15_solute_sl, lateral_loss_no3_n15_solute_sl, &
    & transport_noy, transport_noy_n15, transport_n2o, &
    & transport_n2o_n15, transport_n2, transport_n2_n15, &
    & emission_noy, emission_noy_n15, emission_n2o, &
    & emission_n2o_n15, emission_n2, emission_n2_n15)                                 ! out

    INTEGER,                           INTENT(in)    :: nc                       !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                 !< number of soil layers
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: num_sl_above_bedrock     !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl            !< depth of each soil layer [m]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: rtm_gasdiffusion_act     !< temperature modifier for gas diffusion
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: rmm_gasdiffusion_act     !< moisture modifier for gas diffusion
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: percolation_sl           !< fraction of pool transported up or down per second
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: frac_w_lat_loss_sl       !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep)
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_solute               !< NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: no3_solute               !< NO3 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_solute               !< PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_n15_solute           !< N15H4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: no3_n15_solute           !< N15O3 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: noy                      !< NOy concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: noy_n15                  !< 15NOy concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: n2o                      !< N2O concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: n2o_n15                  !< 15N2O concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: n2                       !< N2 concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: n2_n15                   !< 15N2 concentration [mol / m3]
    REAL(wp),                          INTENT(in)    :: pool_mt_dom(:,:,:)       !< bgcm: DOM pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_solute     !< rate of NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_no3_solute     !< rate of NO3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_solute     !< rate of PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_n15_solute !< rate of N15H4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_no3_n15_solute !< rate of N15O3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_dom_carbon      !< rate of DOM-C leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_dom_nitrogen    !< rate of DOM-N leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_dom_phosphorus  !< rate of DOM-P leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_dom_carbon13    !< rate of DOM-C13 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_dom_carbon14    !< rate of DOM-C14 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_dom_nitrogen15  !< rate of DOM-N15 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_nh4_solute      !< rate of NH4 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_no3_solute      !< rate of NO3 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_po4_solute      !< rate of PO4 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_nh4_n15_solute  !< rate of 14NH4 leaching [mol/m2/timestep]
    REAL(wp), DIMENSION(nc),           INTENT(inout) :: leaching_no3_n15_solute  !< rate of 15NO3 leaching [mol/m2/timestep]
    REAL(wp),                          INTENT(inout) :: transport_mt_dom(:,:,:)  !< bgcm: rate of DOM transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_dom_carbon_sl     !< rate of DOM-C lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_dom_nitrogen_sl   !< rate of DOM-N lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_dom_phosphorus_sl !< rate of DOM-P lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_dom_carbon13_sl   !< rate of DOM-C13 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_dom_carbon14_sl   !< rate of DOM-C14 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_dom_nitrogen15_sl !< rate of DOM-N15 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_nh4_solute_sl     !< rate of NH4 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_no3_solute_sl     !< rate of NO3 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_po4_solute_sl     !< rate of PO4 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_nh4_n15_solute_sl !< rate of 14NH4 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: lateral_loss_no3_n15_solute_sl !< rate of 15NO3 lateral loss [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_noy                  !< vertical transport of NOy [micro-mol/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_noy_n15              !< vertical transport of 15NOy [micro-mol/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_n2o                  !< vertical transport of N2O [micro-mol/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_n2o_n15              !< vertical transport of 15N2O [micro-mol/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_n2                   !< vertical transport of N2 [micro-mol/m3/s]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_n2_n15               !< vertical transport of 15N2 [micro-mol/m3/s]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: emission_noy                   !< soil efflux of NOy [micro-mol/m2/s]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: emission_noy_n15               !< soil efflux of 15NOy [micro-mol/m2/s]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: emission_n2o                   !< soil efflux of N2O [micro-mol/m2/s]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: emission_n2o_n15               !< soil efflux of 15N2O [micro-mol/m2/s]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: emission_n2                    !< soil efflux of N2 [micro-mol/m2/s]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: emission_n2_n15                !< soil efflux of 15N2 [micro-mol/m2/s]
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_phase_transport_wrapper'
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    lateral_loss_dom_carbon_sl(:,:)     = 0.0_wp
    lateral_loss_dom_nitrogen_sl(:,:)   = 0.0_wp
    lateral_loss_dom_phosphorus_sl(:,:) = 0.0_wp
    lateral_loss_dom_carbon13_sl(:,:)   = 0.0_wp
    lateral_loss_dom_carbon14_sl(:,:)   = 0.0_wp
    lateral_loss_dom_nitrogen15_sl(:,:) = 0.0_wp
    lateral_loss_nh4_solute_sl(:,:)     = 0.0_wp
    lateral_loss_no3_solute_sl(:,:)     = 0.0_wp
    lateral_loss_po4_solute_sl(:,:)     = 0.0_wp
    lateral_loss_nh4_n15_solute_sl(:,:) = 0.0_wp
    lateral_loss_no3_n15_solute_sl(:,:) = 0.0_wp
    transport_noy(:,:)                  = 0.0_wp
    transport_noy_n15(:,:)              = 0.0_wp
    transport_n2o(:,:)                  = 0.0_wp
    transport_n2o_n15(:,:)              = 0.0_wp
    transport_n2(:,:)                   = 0.0_wp
    transport_n2_n15(:,:)               = 0.0_wp
    emission_noy(:)                     = 0.0_wp
    emission_noy_n15(:)                 = 0.0_wp
    emission_n2o(:)                     = 0.0_wp
    emission_n2o_n15(:)                 = 0.0_wp
    emission_n2(:)                      = 0.0_wp
    emission_n2_n15(:)                  = 0.0_wp


    !>1.0 transport of DOM
    !>
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              pool_mt_dom(ixC, :, :), transport_mt_dom(ixC, :, :), &
      &                              leaching_dom_carbon(:), lateral_loss_dom_carbon_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              pool_mt_dom(ixN, :, :), transport_mt_dom(ixN, :, :), &
      &                              leaching_dom_nitrogen(:), lateral_loss_dom_nitrogen_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              pool_mt_dom(ixP, :, :), transport_mt_dom(ixP, :, :), &
      &                              leaching_dom_phosphorus(:), lateral_loss_dom_phosphorus_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              pool_mt_dom(ixC13, :, :), transport_mt_dom(ixC13, :, :), &
      &                              leaching_dom_carbon13(:), lateral_loss_dom_carbon13_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              pool_mt_dom(ixC14, :, :), transport_mt_dom(ixC14, :, :), &
      &                              leaching_dom_carbon14(:), lateral_loss_dom_carbon14_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              pool_mt_dom(ixN15, :, :), transport_mt_dom(ixN15, :, :), &
      &                              leaching_dom_nitrogen15(:), lateral_loss_dom_nitrogen15_sl(:,:))

    !>2.0 transport of other solutes (NH4, NO3, PO4, and their isotopes)
    !>
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              nh4_solute(:,:), transport_nh4_solute(:,:), &
      &                              leaching_nh4_solute(:), lateral_loss_nh4_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              no3_solute(:,:), transport_no3_solute(:,:), &
      &                              leaching_no3_solute(:), lateral_loss_no3_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              po4_solute(:,:), transport_po4_solute(:,:), &
      &                              leaching_po4_solute(:), lateral_loss_po4_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              nh4_n15_solute(:,:), transport_nh4_n15_solute(:,:), &
      &                              leaching_nh4_n15_solute(:), lateral_loss_nh4_n15_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
      &                              percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
      &                              no3_n15_solute(:,:), transport_no3_n15_solute(:,:), &
      &                              leaching_no3_n15_solute(:), lateral_loss_no3_n15_solute_sl(:,:))

    !>3.0 gaseous loss of NOy, N2O, and N2, for now following Xu-Ri and Prentice 2008
    !>
    !>  thus ignoring any actual gas-transport and instead assume
    !>  that from each layer a certain fraction excapes directly to the atmosphere. [umol m-3 s-1]
    !>
    ! transport
    transport_noy(:,:)     = -rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * noy(:,:)     / one_day * 1.e6_wp
    transport_noy_n15(:,:) = -rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * noy_n15(:,:) / one_day * 1.e6_wp
    transport_n2o(:,:)     = -rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2o(:,:)     / one_day * 1.e6_wp
    transport_n2o_n15(:,:) = -rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2o_n15(:,:) / one_day * 1.e6_wp
    transport_n2(:,:)      = -rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2(:,:)      / one_day * 1.e6_wp
    transport_n2_n15(:,:)  = -rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2_n15(:,:)  / one_day * 1.e6_wp
    ! emission rates as positive flux to the atmosphere, integrated over the soil column [umol m-2 s-1]
    emission_noy(:)     = SUM(-transport_noy(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2o(:)     = SUM(-transport_n2o(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2(:)      = SUM(-transport_n2(:,:)      * soil_depth_sl(:,:), DIM=2)
    emission_noy_n15(:) = SUM(-transport_noy_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    emission_n2o_n15(:) = SUM(-transport_n2o_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    emission_n2_n15(:)  = SUM(-transport_n2_n15(:,:)  * soil_depth_sl(:,:), DIM=2)
  END SUBROUTINE calc_liquid_phase_transport_wrapper

  ! ======================================================================================================= !
  !>calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !>
  !>  Input:  pools subject to bioturbation, percolation rate, and soil depth
  !>         This explicitly excludes mycorrhiza, as these die when eaten
  !>
  !>  Output: vertical transport rate (mol/m2/timestep)
  !>
  SUBROUTINE calc_bioturbation_transport_wrapper_jsm( &
    & nc, &
    & nsoil_sb, &
    & num_sl_above_bedrock, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & k_bioturb, &
    & nh4_assoc, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & nh4_n15_assoc, &
    & sb_pool_mt, &
    & transport_nh4_assoc, &
    & transport_po4_assoc_fast, &
    & transport_po4_assoc_slow, &
    & transport_po4_occluded, &
    & transport_nh4_n15_assoc, &
    & sb_transport_mt)

    INTEGER,                            INTENT(in)    :: nc                         !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb                   !< number of soil layers
    REAL(wp), DIMENSION(nc),            INTENT(in)    :: num_sl_above_bedrock       !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl              !< depth of each soil layer [m]
    INTEGER,                            INTENT(in)    :: elements_index_map(:)      !< map bgcm element ID -> IDX
    LOGICAL,                            INTENT(in)    :: is_element_used(:)         !< is element in 'elements_index_map' used
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: k_bioturb                  !< diffusion factor for bioturbation
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_assoc                  !< minerally associated NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_assoc_fast             !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_assoc_slow             !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_occluded               !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_n15_assoc              !< N15H4 pool [mol/m3]
    REAL(wp),                           INTENT(in)    :: sb_pool_mt(:,:,:,:)        !< bgcm sb_pool
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_nh4_assoc        !< rate of NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_assoc_fast   !< rate of fast associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_assoc_slow   !< rate of slow associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_occluded     !< rate of occluded PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_nh4_n15_assoc    !< rate of N15H4 transport [mol/m3/timestep]
    REAL(wp),                           INTENT(inout) :: sb_transport_mt(:,:,:,:)   !< bgcm flux: sb_transport flux
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem             !< loop over bgcm elements
    INTEGER                     :: ix_elem           !< index of element in bgcm, used for looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport_wrapper_jsm'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 transport of organic pools
    !>
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        !>  1.1 polymeric litter
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_polymeric_litter, ix_elem, :, :), &
          &                              sb_transport_mt(ix_polymeric_litter, ix_elem, :, :))
        !>  1.2 minerally associated DOM
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_dom_assoc, ix_elem, :, :), &
          &                              sb_transport_mt(ix_dom_assoc, ix_elem, :, :))
        !>  1.3 fungi biomass
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_fungi, ix_elem, :, :), &
          &                              sb_transport_mt(ix_fungi, ix_elem, :, :))
        !>  1.4 microbial biomass
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_microbial, ix_elem, :, :), &
          &                              sb_transport_mt(ix_microbial, ix_elem, :, :))
        !>  1.5 microbial residue
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_residue, ix_elem, :, :), &
          &                              sb_transport_mt(ix_residue, ix_elem, :, :))
        !>  1.6 minerally associated microbial residue
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), &
          &                              sb_pool_mt(ix_residue_assoc, ix_elem, :, :), &
          &                              sb_transport_mt(ix_residue_assoc, ix_elem, :, :))
      END IF
    END DO

    !>2.0 inorganic pools
    !>
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), k_bioturb(:,:), soil_depth_sl(:,:), &
    !                                  nh4_assoc(:,:), &
    !                                  transport_nh4_assoc(:,:))
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), k_bioturb(:,:), soil_depth_sl(:,:), &
    !                                  po4_assoc_fast(:,:), &
    !                                  transport_po4_assoc_fast(:,:))
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), k_bioturb(:,:), soil_depth_sl(:,:), &
    !                                  po4_assoc_slow(:,:), &
    !                                  transport_po4_assoc_slow(:,:))
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), k_bioturb(:,:), soil_depth_sl(:,:), &
    !                                  po4_occluded(:,:), &
    !                                  transport_po4_occluded(:,:))
    ! CALL calc_bioturbation_transport(nc, nsoil_sb, num_sl_above_bedrock(:), k_bioturb(:,:), soil_depth_sl(:,:), &
    !                                  nh4_n15_assoc(:,:), &
    !                                  transport_nh4_n15_assoc(:,:))
  END SUBROUTINE calc_bioturbation_transport_wrapper_jsm

  ! ======================================================================================================= !
  !>calculates particile flux rate given the change in organic pools through formation and loss
  !>  corresponds to the sum of the omega term in Ahrens et al. 2015, eq. S1-S3
  !>
  !>  Input: change of soil pools through formation and loss
  !>
  !>  Output: omega_transport_term
  !>
  FUNCTION calc_particle_fluxrate( &
    & nc, &
    & nsoil_sb, &
    & ixC, &
    & soil_depth_sl, &
    & bulk_dens_corr_sl, &
    & sb_formation_mt, &
    & sb_loss_mt) &
    & RESULT (particle_fluxrate)

    USE mo_sb_constants,           ONLY: carbon_per_dryweight_SOM, rho_bulk_org
    USE mo_jsb_physical_constants, ONLY: molar_mass_C
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                          INTENT(in) :: nc                        !< dimensions
    INTEGER,                                          INTENT(in) :: nsoil_sb                  !< number of soil layers
    INTEGER,                                          INTENT(in) :: ixC                       !< index of C elem in elem map
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in) :: soil_depth_sl             !< depth of soil layer [m]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in) :: bulk_dens_corr_sl         !< bulk density of mineral+organic soil [kg/m3]
    REAL(wp),                                         INTENT(in) :: sb_formation_mt(:,:,:,:)  !< bgcm flux: formation of pools [mol/m3/timestep]
    REAL(wp),                                         INTENT(in) :: sb_loss_mt(:,:,:,:)       !< bgcm flux: loss from pools [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb)                            :: particle_fluxrate         !< particile flux rate [m/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: isoil
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_fluxrate'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate net change of organic matter cummulative down the profile [kg/m2/dtime]
    !>
    !>  NOTE omission of woody litter in the first layer, which is assumed not to be part of the profile
    !>
    particle_fluxrate(:,1) = soil_depth_sl(:,1)                           &
      & * ( sb_formation_mt(ix_soluable_litter, ixC, :, 1)  - sb_loss_mt(ix_soluable_litter, ixC, :, 1)   &
      &   + sb_formation_mt(ix_polymeric_litter, ixC, :, 1) - sb_loss_mt(ix_polymeric_litter, ixC, :, 1)  &
      &   + sb_formation_mt(ix_dom, ixC, :, 1)              - sb_loss_mt(ix_dom, ixC, :, 1)               &
      &   + sb_formation_mt(ix_dom_assoc, ixC, :, 1)        - sb_loss_mt(ix_dom_assoc, ixC, :, 1)         &
      &   + sb_formation_mt(ix_fungi, ixC, :, 1)            - sb_loss_mt(ix_fungi, ixC, :, 1)             &
      &   + sb_formation_mt(ix_mycorrhiza, ixC, :, 1)       - sb_loss_mt(ix_mycorrhiza, ixC, :, 1)        &
      &   + sb_formation_mt(ix_microbial, ixC, :, 1)        - sb_loss_mt(ix_microbial, ixC, :, 1)         &
      &   + sb_formation_mt(ix_residue, ixC, :, 1)          - sb_loss_mt(ix_residue, ixC, :, 1)           &
      &   + sb_formation_mt(ix_residue_assoc, ixC, :, 1)    - sb_loss_mt(ix_residue_assoc, ixC, :, 1)  )  &
      & / (carbon_per_dryweight_SOM / molar_mass_C * 1000._wp)

    ! For lower soil layers all compartments (including woody litter) are considered
    DO isoil = 2,nsoil_sb
      particle_fluxrate(:,isoil) = particle_fluxrate(:, isoil - 1) + soil_depth_sl(:, isoil)                         &
      & * ( sb_formation_mt(ix_soluable_litter, ixC, :, isoil)     - sb_loss_mt(ix_soluable_litter, ixC, :, isoil)   &
      &   + sb_formation_mt(ix_polymeric_litter, ixC, :, isoil)    - sb_loss_mt(ix_polymeric_litter, ixC, :, isoil)  &
      &   + sb_formation_mt(ix_woody_litter, ixC, :, isoil)        - sb_loss_mt(ix_woody_litter, ixC, :, isoil)      &
      &   + sb_formation_mt(ix_dom, ixC, :, isoil)                 - sb_loss_mt(ix_dom, ixC, :, isoil)               &
      &   + sb_formation_mt(ix_dom_assoc, ixC, :, isoil)           - sb_loss_mt(ix_dom_assoc, ixC, :, isoil)         &
      &   + sb_formation_mt(ix_fungi, ixC, :, isoil)               - sb_loss_mt(ix_fungi, ixC, :, isoil)             &
      &   + sb_formation_mt(ix_mycorrhiza, ixC, :, isoil)          - sb_loss_mt(ix_mycorrhiza, ixC, :, isoil)        &
      &   + sb_formation_mt(ix_microbial, ixC, :, isoil)           - sb_loss_mt(ix_microbial, ixC, :, isoil)         &
      &   + sb_formation_mt(ix_residue, ixC, :, isoil)             - sb_loss_mt(ix_residue, ixC, :, isoil)           &
      &   + sb_formation_mt(ix_residue_assoc, ixC, :, isoil)       - sb_loss_mt(ix_residue_assoc, ixC, :, isoil) )   &
      & / (carbon_per_dryweight_SOM / molar_mass_C * 1000._wp)
      ! should be the same -- but is not bitidentical...:
      ! particle_fluxrate(:,isoil) = particle_fluxrate(:,isoil-1) + soil_depth_sl(:,isoil) &
      !   & * SUM((sb_formation_mt(: ,ixC,:,isoil) - sb_loss_mt(: ,ixC,:,isoil)),DIM=1)      &
      !   & / (carbon_per_dryweight_SOM / molar_mass_C * 1000._wp)
    ENDDO

    !>2.0 convert change in pools in density change and thus transport term (m per dtime)
    !>
    WHERE(bulk_dens_corr_sl(:,:) > eps8)
      particle_fluxrate(:,:) = particle_fluxrate(:,:) / rho_bulk_org
    ENDWHERE
  END FUNCTION calc_particle_fluxrate

  ! ======================================================================================================= !
  !>calculates matter transport due to particle flux, by calling 'calc_particle_transport()'
  !>
  !>  Input:  pools subject to particle flux, and soil depth
  !>
  !>  Output: vertical transport rate [mol/m2/timestep]
  !>
  SUBROUTINE calc_particle_transport_wrapper( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & particle_fluxrate, &
    & nh4_solute, &
    & nh4_assoc, &
    & no3_solute, &
    & po4_solute, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & po4_primary, &
    & nh4_n15_solute, &
    & nh4_n15_assoc, &
    & no3_n15_solute, &
    & sb_pool_mt, &
    & transport_nh4_solute, &
    & transport_nh4_assoc, &
    & transport_no3_solute, &
    & transport_po4_solute, &
    & transport_po4_assoc_fast, &
    & transport_po4_assoc_slow, &
    & transport_po4_occluded, &
    & transport_po4_primary, &
    & transport_nh4_n15_solute, &
    & transport_nh4_n15_assoc, &
    & transport_no3_n15_solute, &
    & sb_transport_mt)

    INTEGER,                            INTENT(in)    :: nc                         !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb                   !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl              !< thickness of each soil layer [m]
    INTEGER,                            INTENT(in)    :: elements_index_map(:)      !< map bgcm element ID -> IDX
    LOGICAL,                            INTENT(in)    :: is_element_used(:)         !< is element in 'elements_index_map' used
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: particle_fluxrate          !< particle fluxrate (m/timestep)
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_solute                 !< NH4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_assoc                  !< minerally associated NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: no3_solute                 !< NO3 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_solute                 !< PO4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_assoc_fast             !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_assoc_slow             !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_occluded               !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: po4_primary                !< occluded PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_n15_solute             !< N15H4 pool in soil solution [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_n15_assoc              !< minerally associated N15H4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: no3_n15_solute             !< N15O3 pool in soil solution [mol/m3]
    REAL(wp),                           INTENT(in)    :: sb_pool_mt(:,:,:,:)        !< bgcm sb_pool: organic pools [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_nh4_solute       !< rate of solute NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_nh4_assoc        !< rate of associated NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_no3_solute       !< rate of solute NO3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_solute       !< rate of solute PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_assoc_fast   !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_assoc_slow   !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_occluded     !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_po4_primary      !< rate of associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_nh4_n15_solute   !< rate of solute N15H4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_no3_n15_solute   !< rate of solute N15O3 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport_nh4_n15_assoc    !< rate of associated N15H4 transport [mol/m3/timestep]
    REAL(wp),                           INTENT(inout) :: sb_transport_mt(:,:,:,:)   !< bgcm flux: rate of sb pool transport [mol/m3/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                           :: isoil                        !< loop over soil layers
    INTEGER                           :: ielem                        !< loop over bgcm elements
    INTEGER                           :: ix_elem                      !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc, nsoil_sb) :: particle_fluxrate_woody_sl   !< particle fluxrate for woody litter
    CHARACTER(len=*), PARAMETER       :: routine = TRIM(modname)//':calc_particle_transport_wrapper'
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 calc particle fluxrate for woody litter
    !>
    ! first layer of woody litter does not belong to profile
    DO isoil = 1,nsoil_sb
      particle_fluxrate_woody_sl(:,isoil) = particle_fluxrate(:,isoil) - particle_fluxrate(:,1)
    ENDDO

    !>1.0 transport of organic pools
    !>
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        !>  1.1 soluable litter
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:),  &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_soluable_litter, ix_elem, :, :), &
          & sb_transport_mt(ix_soluable_litter, ix_elem, :, :) )
        !>  1.2 polymeric litter
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_polymeric_litter, ix_elem, :, :), &
          & sb_transport_mt(ix_polymeric_litter, ix_elem, :, :) )
        !>  1.3 woody litter
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate_woody_sl(:,:), &
          & sb_pool_mt(ix_woody_litter, ix_elem, :, :), &
          & sb_transport_mt(ix_woody_litter, ix_elem, :, :) )
        !>  1.4 DOM
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_dom, ix_elem, :, :), &
          & sb_transport_mt(ix_dom, ix_elem, :, :) )
        !>  1.5 minerally associated DOM
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_dom_assoc, ix_elem, :, :), &
          & sb_transport_mt(ix_dom_assoc, ix_elem, :, :) )
        !>  1.6 fungi biomass
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_fungi, ix_elem, :, :), &
          & sb_transport_mt(ix_fungi, ix_elem, :, :) )
        !>  1.7 mycorrhiza biomass
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_mycorrhiza, ix_elem, :, :), &
          & sb_transport_mt(ix_mycorrhiza, ix_elem, :, :) )
        !>  1.8 microbial biomass
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_microbial, ix_elem, :, :), &
          & sb_transport_mt(ix_microbial, ix_elem, :, :) )
        !>  1.9 microbial residue
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_residue, ix_elem, :, :), &
          & sb_transport_mt(ix_residue, ix_elem, :, :) )
        !>  1.10 minerally associated microbial residue
        !>
        CALL calc_particle_transport( &
          & nc, nsoil_sb, soil_depth_sl(:,:), &
          & particle_fluxrate(:,:), &
          & sb_pool_mt(ix_residue_assoc, ix_elem, :, :), &
          & sb_transport_mt(ix_residue_assoc, ix_elem, :, :) )
      END IF
    END DO

    !>2.0 inorganic pools
    !>
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & nh4_solute(:,:), transport_nh4_solute(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & nh4_assoc(:,:), transport_nh4_assoc(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & no3_solute(:,:), transport_no3_solute(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & po4_solute(:,:), transport_po4_solute(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & po4_assoc_fast(:,:), transport_po4_assoc_fast(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & po4_assoc_slow(:,:), transport_po4_assoc_slow(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & po4_occluded(:,:), transport_po4_occluded(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & po4_primary(:,:), transport_po4_primary(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & nh4_n15_solute(:,:), transport_nh4_n15_solute(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & nh4_n15_assoc(:,:), transport_nh4_n15_assoc(:,:))
    CALL calc_particle_transport( &
      & nc, nsoil_sb, soil_depth_sl(:,:), &
      & particle_fluxrate(:,:), &
      & no3_n15_solute(:,:), transport_no3_n15_solute(:,:))
  END SUBROUTINE calc_particle_transport_wrapper

  ! ======================================================================================================= !
  !>Advective Transport in the vertical dimension
  !>  corresponds to the particule flux transport term in Ahrens et al. 2015, eq 4-7
  !>
  !> Input: pool to be transported for a specific element, percolation rate, and soil depth
  !>
  !> Output: vertical transport rate (mol/m2/s)
  !>
  SUBROUTINE calc_particle_transport( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & particle_fluxrate, &
    & pool, &
    & transport)

    INTEGER,                                          INTENT(in)    :: nc                     !< dimensions
    INTEGER,                                          INTENT(in)    :: nsoil_sb               !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: particle_fluxrate      !< fraction of pool transported up or down per timestep
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)    :: pool                   !< current state of pool to be transported
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(inout) :: transport              !< transport between layers per timestep
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)     :: hlp1
    REAL(wp), DIMENSION(nc)     :: hlp2
    INTEGER                     :: ic, isoil, isoil_hlp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_particle_transport'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 advective transport in the vertical axis of the soil (i.e. down)
    !>
    ! loop over: soil layers & points in chunk & soil layers limited by the 1st soil layer loop
    ! we can do '1,n-1' here, because by definition, there cannot be any upward transport from the uppermost layer
    DO isoil = 1,nsoil_sb-1
      DO ic = 1,nc
        IF (particle_fluxrate(ic, isoil) < 0.0_wp) THEN
          hlp1(ic) = particle_fluxrate(ic, isoil) * pool(ic, isoil + 1)
        ELSE
          ! when particle transport pushes down material, check if the fluxrate exceeds the soil layer depth
          isoil_hlp = isoil
          hlp1(ic)  = 0._wp
          hlp2(ic)  = particle_fluxrate(ic, isoil)
          ! DO loop: only when fluxrate is more than the soil layer depth,
          !          push the whole layer downwards and store the rest in hlp2
          DO WHILE (hlp2(ic) > soil_depth_sl(ic, isoil_hlp) .AND. isoil_hlp > 0)
            hlp1(ic)  = hlp1(ic) + soil_depth_sl(ic, isoil_hlp) * pool(ic, isoil_hlp)
            hlp2(ic)  = hlp2(ic) - soil_depth_sl(ic, isoil_hlp)
            isoil_hlp = isoil_hlp - 1
          ENDDO
          ! IF 1. fluxrate does not exceeds depth of layer isoil, material are transported from the layer isoil
          !    2. not all isoil layers are pushed down, extra material are transported from layer i
          ! ELSE (i<eps4) all isoil layers are pushed down (fluxrate is more than the total depth of isoil layers)
          IF (isoil_hlp > 0) THEN
            hlp1(ic) = hlp1(ic) + hlp2(ic) * pool(ic, isoil_hlp)
          ENDIF
        ENDIF
        transport(ic, isoil)     = transport(ic, isoil)     - hlp1(ic) / soil_depth_sl(ic, isoil)
        transport(ic, isoil + 1) = transport(ic, isoil + 1) + hlp1(ic) / soil_depth_sl(ic, isoil + 1)
      ENDDO
    ENDDO

    !>2.0 lower boundary condition, ... ??  @TODO
    !>
    !  based on the assumption that there is no upward flux, and that all other
    !  material is lost
    hlp1(:)               = 0.0_wp
    transport(:,nsoil_sb) = transport(:, nsoil_sb) - hlp1(:) / soil_depth_sl(:, nsoil_sb)
  END SUBROUTINE calc_particle_transport

  ! ======================================================================================================= !
  !>Advective Transport in the vertical dimension
  !>  corresponds to the water flux transport term in Ahrens et al. 2015, eq 5
  !>
  !>  Input: pool to be transported, percolation rate, lateral loss fraction, and soil depth
  !>
  !>  Output: vertical transport rate and lateral loss rate (mol/m3/timestep)
  !>
  SUBROUTINE calc_liquid_phase_transport( &
    & nc, &
    & nsoil_sb, &
    & num_sl_above_bedrock, &
    & percolation_sl, &
    & frac_w_lat_loss_sl, &
    & soil_depth_sl, &
    & pool, &
    & transport, &
    & leaching, &
    & lateral_loss_material)

    INTEGER,                            INTENT(in)    :: nc                     !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb               !< number of soil layers
    REAL(wp), DIMENSION(nc),            INTENT(in)    :: num_sl_above_bedrock   !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: percolation_sl         !< fraction of pool transported up or down (vertical) per second
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: frac_w_lat_loss_sl     !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep)
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: pool                   !< current state of pool to be transported [mol m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport              !< transport between layers per timestep [mol m-3]
    REAL(wp), DIMENSION(nc),            INTENT(out)   :: leaching               !< transport through the lower surface of the lowest
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(out)   :: lateral_loss_material  !< lateral loss flux of material
                                                                                !! layer per timestep [mol m-3 timestep-1]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)         :: hlp1
    INTEGER                         :: ic
    INTEGER                         :: is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_phase_transport'
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    leaching(:)                 = 0.0_wp
    lateral_loss_material(:,:)  = 0.0_wp

    !>1.0 advective transport in the vertical axis of the soil (i.e. down)
    !>
    ! loop over soil layers
    ! we can do '1,n-1' here, because by definition, there cannot be any upward transport from the uppermost layer
    DO ic = 1,nc
      DO is = 1,(INT(num_sl_above_bedrock(ic)) - 1)
        IF (percolation_sl(ic,is) > 0.0_wp) THEN
          hlp1(ic) = percolation_sl(ic,is) * dtime / 1000._wp * pool(ic,is) * soil_depth_sl(ic,is)
        !! @TODO BLARPP for now no upward transport - need to properly think about this in update_spq_physics
        ELSE
          hlp1(ic) = 0.0_wp ! percolation_sl(:,is) * dtime/1000._wp * pool(:,is+1) * soil_depth_sl(:,is+1)
        END IF
        transport(ic,is)   = transport(ic,is)   - hlp1(ic) / soil_depth_sl(ic,is)
        transport(ic,is+1) = transport(ic,is+1) + hlp1(ic) / soil_depth_sl(ic,is+1)
      ENDDO
    ENDDO

    !>2.0 lower boundary condition,
    !>
    !  based on the assumption that there is no upward flux from below,
    !  and that all other material is lost
    DO ic = 1,nc
      IF (percolation_sl(ic, INT(num_sl_above_bedrock(ic))) > 0.0_wp) THEN
        hlp1(ic) = percolation_sl(ic, INT(num_sl_above_bedrock(ic))) * dtime / 1000._wp &
          &        * pool(ic, INT(num_sl_above_bedrock(ic))) * soil_depth_sl(ic, INT(num_sl_above_bedrock(ic)))
      ELSE
        hlp1(ic) = 0.0_wp
      END IF
      transport(ic, INT(num_sl_above_bedrock(ic))) = transport(ic, INT(num_sl_above_bedrock(ic))) - hlp1(ic) &
        &                                            / soil_depth_sl(ic, INT(num_sl_above_bedrock(ic)))
      leaching(ic)                                 = hlp1(ic)
    END DO

    !>3.0 lateral loss
    !>
    !>  calculation of the the flux of lateral loss of material per soil layer
    !>  lateral_loss material flux =
    !>    "fraction water lateral_loss" * material_pool (mol m-3) * soil_depth (m) * dtime / 1000._wp
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        IF (frac_w_lat_loss_sl(ic, is) > 0.0_wp) THEN
          hlp1(ic) = frac_w_lat_loss_sl(ic, is) * pool(ic, is) * soil_depth_sl(ic, is) * dtime / 1000._wp
        ELSE
          hlp1(ic) = 0.0_wp
        END IF
        lateral_loss_material(ic, is) = hlp1(ic) / soil_depth_sl(ic, is)
      ENDDO
    ENDDO
  END SUBROUTINE calc_liquid_phase_transport

  ! ======================================================================================================= !
  !>calculates rate of bioturbation
  !>
  !>  Input: bulk density
  !>
  !>  Output: vertical transport rate (1/s)
  !>
  FUNCTION calc_bioturbation_rate( &
    & nc, &
    & nsoil_sb, &
    & num_sl_above_bedrock, &
    & soil_depth_sl, &
    & bulk_dens_corr_sl, &
    & root_fraction_sl) &
    & RESULT(k_bioturb)

    USE mo_sb_constants,          ONLY: k_diff_org
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc                     !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb               !< number of soil layers
    REAL(wp), DIMENSION(nc),            INTENT(in)    :: num_sl_above_bedrock   !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl          !< thickness of soil layers [m]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: bulk_dens_corr_sl      !< OM corrected bulk density [kg m-3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: root_fraction_sl       !< root fraction in soil layer
    REAL(wp), DIMENSION(nc, nsoil_sb)                 :: k_bioturb              !< diffusion factor for bioturbation
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb)               :: hlp1
    INTEGER                                         :: ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_rate'

    !> init out var
    !>
    k_bioturb(:,:) = 0.0_wp

    !>1.0 calculate root density relative to top-soil root density
    !>
    hlp1(:,:) = 0.0_wp
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        hlp1(ic, is) = root_fraction_sl(ic, is) / soil_depth_sl(ic, is)
      ENDDO
    ENDDO

    DO ic = 1,nc
      DO is = INT(num_sl_above_bedrock(ic)), 1, -1
        IF (hlp1(ic, 1) > eps8) THEN
          hlp1(ic, is) = hlp1(ic, is) / hlp1(ic, 1)
        ELSE
          hlp1(ic, is) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO

    !>2.0 adjust bioturbation coefficient by bulk density and root density
    !>
    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        IF (bulk_dens_corr_sl(ic, is) > eps8) THEN
          k_bioturb(ic, is) = hlp1(ic, is) * k_diff_org / bulk_dens_corr_sl(ic, is)
        ELSE
          k_bioturb(ic, is) = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
  END FUNCTION calc_bioturbation_rate

  ! ======================================================================================================= !
  !>calculates matter transport by bioturbation
  !>
  !>  corresponds to the diffusion via bioturbation term in Ahrens et al. 2015, eq. 4, 6 and 7
  !>
  SUBROUTINE calc_bioturbation_transport( &
    & nc, &
    & nsoil_sb, &
    & num_sl_above_bedrock, &
    & k_bioturb, &
    & soil_depth_sl, &
    & pool, &
    & transport)

    INTEGER,                            INTENT(in)    :: nc                   !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb             !< number of soil layers
    REAL(wp), DIMENSION(nc),            INTENT(in)    :: num_sl_above_bedrock !< number of soil layers above bedrock, i.e., with layer thickness > eps8
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: k_bioturb            !< diffusion factor for bioturbation [m2/s]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl        !< depth of each soil layer [m]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: pool                 !< bgcm pool: current state of pool to be transported
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: transport            !< bgcm flux: transport between layers per timestep [mol/m2/timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)     :: hlp1
    REAL(wp), DIMENSION(nc)     :: hlp2
    INTEGER                     :: ic, is
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport'
    ! ----------------------------------------------------------------------------------------------------- !

    !> 1.0 diffusion transport in the vertical axis of the soil (i.e. down)
    !>
    !! we can do 1,n-1 here, because by definition, there cannot be any upward transport from the uppermost layer
    !!
    !! we don't need soil_layer (nsoil_sb) here, because the assumption is made that the gradient at the bottom of the
    !! soil column is zero
    DO ic = 1,nc
      DO is = 1,(INT(num_sl_above_bedrock(ic)) - 1)
        ! distance between the centres of the soil layers
        hlp1(ic) = (soil_depth_sl(ic,is) / 2.0_wp + soil_depth_sl(ic,is + 1) / 2.0_wp)

        IF (hlp1(ic) > eps8) THEN
          hlp2(ic) = k_bioturb(ic,is) * dtime * (pool(ic,is) - pool(ic,is + 1)) / hlp1(ic)
        ELSE
          hlp2(ic) = 0.0_wp
        ENDIF
        transport(ic,is)     = transport(ic,is)     - hlp2(ic) / soil_depth_sl(ic,is)
        transport(ic,is + 1) = transport(ic,is + 1) + hlp2(ic) / soil_depth_sl(ic,is + 1)
      ENDDO
    ENDDO
  END SUBROUTINE calc_bioturbation_transport

#endif
END MODULE mo_q_sb_jsm_transport
