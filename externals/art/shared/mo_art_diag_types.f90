!
! mo_art_diag_types
! This module provides datastructures for diagnostic variables in ART
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

MODULE mo_art_diag_types
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_fortran_tools,                 ONLY: t_ptr_2d3d
! ART

  IMPLICIT NONE

  INTEGER, ALLOCATABLE  :: &
    & art_diag_tracer_index(:,:)

  PRIVATE

  PUBLIC :: t_art_aeronet, t_art_max_layer
  PUBLIC :: t_art_diag, art_diag_tracer_index
  PUBLIC :: art_deallocate_diag

  TYPE t_art_max_layer
    REAL(wp)             ::  &
      &  pres_bot,           & !< bottom pressure [Pa]
      &  pres_top              !< top pressure    [Pa]
    REAL(wp), POINTER    ::  &
      &  maximum(:,:)          !< maximum value in pressure layer pres_bot to pres_top
  END TYPE t_art_max_layer

  TYPE t_art_cloud_height
    REAL(wp)             ::  &
      &  threshold             !< threshold in [mg/m3] that defines volc ash cloud
    REAL(wp), POINTER    ::  &
      &   base(:,:),         & !< lowest height in [m] where threshold exceeded
      &   top(:,:)             !< greatest height           -"-
  ENDTYPE t_art_cloud_height

  TYPE t_art_radio_diag
    REAL(wp), POINTER    ::  &
      &  wetdepo(:,:),       & !< Accumulated wet deposition of radionuclides [Bq/m2]
      &  drydepo(:,:),       & !< Accumulated dry deposition of radionuclides [Bq/m2]
      &  avg(:,:,:),         & !< Averaged air concentration of radionuclides [Bq/m3]
      &  maxtint(:,:,:)        !< Maximum air concentration of radionuclides in time interval [Bq/m3]
    TYPE(t_art_max_layer), POINTER :: &
      &  maxtint_layer(:) => NULL() !< Maximal maximum air concentration of radionuclides
                                    !  in time interval between given pressure levels [Bq/m3]
  END TYPE t_art_radio_diag

  TYPE t_art_aeronet
    REAL(wp), POINTER    ::  &
      &  tau(:,:,:),         & !< aerosol optical depth at Aeronet wavelengths
      &  tau_vi(:,:)           !< total column of aerosol optical depth at Aeronet wavelengths
  END TYPE t_art_aeronet

  TYPE t_art_ceilo
    REAL(wp), POINTER    ::  &
      &  bsc(:,:,:),         & !< aerosol backscatter at Lidar wavelengths
      &  ceil_bsc(:,:,:),    & !< aerosol attenuated backscatter at Lidar wavelengths seen from ground
      &  sat_bsc(:,:,:)        !< aerosol attenuated backscatter at Lidar wavelengths seen from space
  END TYPE t_art_ceilo

  ! ----------------------------------
  ! --- Diagnostic fields for ART
  ! ----------------------------------

  TYPE t_art_diag
! Cloud Microphysics and Aerosol-Cloud-Interactions
    ! Warm Phase
    REAL(wp),POINTER          ::  &
      &  ncalls_warm(:,:,:)      => NULL(), & !< Number of calls of activation routine (accumulated)
      &  aci_nnuctot_warm(:,:,:) => NULL(), & !< Number of activated particles (accumulated, warm-phase)
      &  aci_nnucfhh_warm(:,:,:) => NULL(), & !< Number of activated particles according to FHH theory
                                              !   (accumulated, warm-phase)
      &  smax_water(:,:,:)       => NULL()    !< Maximum supersaturation over liquid water
! Cold Phase
    REAL(wp),POINTER          ::  &
      &  ncalls_cold(:,:,:)      => NULL(), & !< Number of calls of ice nucleation routine (accumulated)
      &  aci_nnuctot_cold(:,:,:) => NULL(), & !< Number of total nucleated (hom. freez. + het. nuc.)
                                              !   particles (accumulated, cold-phase)
      &  aci_nnuchet_cold(:,:,:) => NULL(), & !< Number of het. nucleated particles (accumulated, cold-phase)
      &  smax_ice(:,:,:)         => NULL()    !< Maximum supersaturation over ice
    ! Independent of  Phase
! Aerosol properties
    ! Input for KL06 nucleation scheme
    REAL(wp),POINTER          ::  &
      &  ndust_tot(:,:,:)        => NULL()           !< Total number of mineral dust particles [m-3]
! Debug Variables
    REAL(wp),POINTER          ::  &
      &  ART_3D_DBG01(:,:,:)     => NULL() , &
      &  art_3d_dbg_init (:,:,:) => NULL()
    REAL(wp),POINTER          ::  &
      &  ART_2D_DBG01(:,:)       => NULL()
! Multiple Fields Variables
    TYPE(t_art_aeronet),ALLOCATABLE :: &
      &  dust_aeronet(:),              & !< Optical thickness dust (9 spectral bands)
      &  so4_sol_aeronet(:),           & !< Optical thickness so4
      &  ash_insol_aeronet(:),         & !< Optical thickness ash insol
      &  ash_mixed_aeronet(:),         & !< Optical thickness ash mixed
      &  ash_giant_aeronet(:),         & !< Optical thickness ash giant
      &  seas_aeronet(:),              & !< Optical thickness seasalt (9 spectral bands)
      &  volc_aeronet(:),              & !< Optical thickness volcanic ash (9 spectral bands)
      &  soot_aeronet(:)                 !< Optical thickness soot (9 spectral bands)
    TYPE(t_art_ceilo),ALLOCATABLE ::   &
      &  dust_ceilo(:),                & !< Backscatter coefficient dust (3 spectral bands)
      &  seas_ceilo(:),                & !< Backscatter coefficient seasalt (3 spectral bands)
      &  volc_ceilo(:),                & !< Backscatter coefficient volcanic ash (3 spectral bands)
      &  so4_sol_ceilo(:),             & !< Backscatter coefficient volcanic ash(3 spectral bands)
      &  ash_insol_ceilo(:),           &
      &  ash_mixed_ceilo(:),           &
      &  ash_giant_ceilo(:),           &
      &  soot_ceilo(:),                & !< Backscatter coefficient soot (3 spectral bands)
!     &  dust_ext(:),                  & !< Extinciton coefficient dust (3 spectral bands)
!     &  seas_ext(:),                  & !< Extinction coefficient seasalt (3 spectral bands)
!     &  volc_ext(:),                  & !< Extinction coefficient volcanic ash (3 spectral bands)
      &  dust_att(:),                  & !< Attenuated Backscatter dust (3 spectral bands)
      &  seas_att(:),                  & !< Attenuated Backscatter seasalt (3 spectral bands)
      &  volc_att(:),                  & !< Attenuated Backscatter volcanic ash (3 spectral bands)
      &  soot_att(:),                  & !< Attenuated Backscatter soot (3 spectral bands)
      &  so4_sol_att(:),               & !< Attenuated Backscatter so4
      &  ash_insol_att(:),             &
      &  ash_mixed_att(:),             &
      &  ash_giant_att(:),             &
      &  dust_sat(:),                  & !< Attenuated Backscatter from satellite dust (3 spectral bands)
      &  ash_insol_sat(:),             &
      &  ash_mixed_sat(:),             &
      &  ash_giant_sat(:),             &
      &  so4_sat(:),                   & !< Attenuated Backscatter from satellite so4
      &  seas_sat(:),                  & !< Attenuated Backscatter from satellite seasalt (3 spectral bands)
      &  volc_sat(:),                  & !< Attenuated Backscatter from satellite volcanic ash (3 spectral bands)
      &  soot_sat(:)                     !< Attenuated Backscatter from satllite soot (3 spectral bands)
    REAL(wp),POINTER      ::           &
      &  art_photolysis(:,:,:,:) => NULL()
    TYPE(t_ptr_2d3d),ALLOCATABLE ::    &
      &  art_cgaml_ptr(:),             &
      &  art_reac_rates_ptr(:),        &
      &  art_NAT_rel_diff_ptr(:),      &
      &  art_NAT_sedi_vel_ptr(:),      &
      &  art_photolysis_ptr(:)
! Volcanic ash diagnostics
    REAL(wp),POINTER          ::         &
      &  ash_total_mc(:,:,:)  => NULL(), & !< Total volcanic ash mass concentration [kg/m3]
      &  ash_total_mc_vi(:,:) => NULL(), & !< Total column of volcanic ash mass concentration [kg/m2]
      &  ash_hml_max(:,:)     => NULL()    !< Height of maximal total volcanic ash mass concentration [m]
    TYPE(t_art_max_layer),POINTER ::   &
      &  ash_max_total_mc(:)  => NULL()     !< Maximal total ash mass between given pressure levels [kg/m3]
    TYPE(t_art_cloud_height),POINTER :: &
      &  ash_cloud(:)         => NULL()     !< base and top in [m] of volcanic ash cloud for different
                                           !  threshold conentrations in [ug/m3] as cloud definitions
! Sea salt diagnostics
    REAL(wp),POINTER          ::         &
      &  seas_total_mc(:,:,:) => NULL()    !< Total sea salt mass concentration [kg/m3]

! Mineral dust diagnostics
    REAL(wp),POINTER          ::          &
      &  dust_total_mc(:,:,:)  => NULL(), & !< Total mineral dust mass concentration [kg/m3]
      &  dust_total_mc_vi(:,:) => NULL()    !< Total column of mineral dust mass concentration [kg/m2]
    TYPE(t_art_max_layer),POINTER ::      &
      &  dust_max_total_mc(:)  => NULL()    !< Maximal total mineral dust mass between given pressure levels [kg/m3]
! Ustar threshold (for dust emission)
    REAL(wp),POINTER          ::         &
      &  ustar_threshold(:,:) => NULL(), & !< threshold friction velocity [m/s] (see Vogel et.al 2006, eq.3.4)
      &  ustar(:,:)           => NULL()    !< friction velocity [m/s]

! Tracer emission and deposition diagnostics
    REAL(wp),POINTER             ::    &  !< new diagnostics container for accumulated tracers (nproma, nblks_c, ntracer_dry_depo)
      &  acc_drydepo(:,:,:)       => NULL(), &  !< dry deposition 
      &  acc_sedim(:,:,:)         => NULL(), &  !< sedimentation
      &  acc_wetdepo_gscp(:,:,:)  => NULL(), &  !< wet deposition by grid scale precipitation
      &  acc_wetdepo_con(:,:,:)   => NULL(), &  !< wet deposition by convective precipitation
      &  acc_wetdepo_rrsfc(:,:,:) => NULL(), &  !< wet deposition when precipitation reaches the surface
      &  emiss(:,:,:)             => NULL(), &  !< emission of tracer
      &  acc_emiss(:,:,:)         => NULL()     !< accumulated emission of tracer
   TYPE(t_ptr_2d3d),ALLOCATABLE ::     &
      &  acc_drydepo_ptr(:),           &  !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]
      &  acc_sedim_ptr(:),             &  !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]
      &  acc_wetdepo_gscp_ptr(:),      &  !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]
      &  acc_wetdepo_con_ptr(:),       &  !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]
      &  acc_wetdepo_rrsfc_ptr(:),     &  !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]
      &  emiss_ptr(:),                 &  !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]
      &  acc_emiss_ptr(:)                 !< pointer array (one pointer for each tracer) [tracer-unit m-2], eg [mug m-2]

! Radioactive tracers diagnostics
    TYPE(t_art_radio_diag),POINTER ::  &
      &  radioact(:) => NULL()                     !< Diagnostic fields for radionuclides
  END TYPE t_art_diag

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_deallocate_diag(this)
!<
! SUBROUTINE art_deallocate_diag
! This subroutine does the clean up for the ART diag structure
! Based on: -
! Part of Module: mo_art_diag_types
! Author: Sven Werchner, KIT
! Initial Release: 2020-10-05
! Modifications:
!>

  TYPE(t_art_diag), INTENT(INOUT) ::     this


  IF(ALLOCATED(this%dust_aeronet))          DEALLOCATE(this%dust_aeronet)
  IF(ALLOCATED(this%so4_sol_aeronet))       DEALLOCATE(this%so4_sol_aeronet)
  IF(ALLOCATED(this%ash_insol_aeronet))     DEALLOCATE(this%ash_insol_aeronet)
  IF(ALLOCATED(this%ash_mixed_aeronet))     DEALLOCATE(this%ash_mixed_aeronet)
  IF(ALLOCATED(this%ash_giant_aeronet))     DEALLOCATE(this%ash_giant_aeronet)
  IF(ALLOCATED(this%seas_aeronet))          DEALLOCATE(this%seas_aeronet)
  IF(ALLOCATED(this%volc_aeronet))          DEALLOCATE(this%volc_aeronet)

  IF(ALLOCATED(this%dust_ceilo))            DEALLOCATE(this%dust_ceilo)
  IF(ALLOCATED(this%seas_ceilo))            DEALLOCATE(this%seas_ceilo)
  IF(ALLOCATED(this%volc_ceilo))            DEALLOCATE(this%volc_ceilo)
  IF(ALLOCATED(this%so4_sol_ceilo))         DEALLOCATE(this%so4_sol_ceilo)
  IF(ALLOCATED(this%ash_insol_ceilo))       DEALLOCATE(this%ash_insol_ceilo)
  IF(ALLOCATED(this%ash_mixed_ceilo))       DEALLOCATE(this%ash_mixed_ceilo)
  IF(ALLOCATED(this%ash_giant_ceilo))       DEALLOCATE(this%ash_giant_ceilo)
  IF(ALLOCATED(this%dust_att))              DEALLOCATE(this%dust_att)
  IF(ALLOCATED(this%seas_att))              DEALLOCATE(this%seas_att)
  IF(ALLOCATED(this%volc_att))              DEALLOCATE(this%volc_att)
  IF(ALLOCATED(this%so4_sol_att))           DEALLOCATE(this%so4_sol_att)
  IF(ALLOCATED(this%ash_insol_att))         DEALLOCATE(this%ash_insol_att)
  IF(ALLOCATED(this%ash_mixed_att))         DEALLOCATE(this%ash_mixed_att)
  IF(ALLOCATED(this%ash_giant_att))         DEALLOCATE(this%ash_giant_att)
  IF(ALLOCATED(this%dust_sat))              DEALLOCATE(this%dust_sat)
  IF(ALLOCATED(this%ash_insol_sat))         DEALLOCATE(this%ash_insol_sat)
  IF(ALLOCATED(this%ash_mixed_sat))         DEALLOCATE(this%ash_mixed_sat)
  IF(ALLOCATED(this%ash_giant_sat))         DEALLOCATE(this%ash_giant_sat)
  IF(ALLOCATED(this%so4_sat))               DEALLOCATE(this%so4_sat)
  IF(ALLOCATED(this%seas_sat))              DEALLOCATE(this%seas_sat)
  IF(ALLOCATED(this%volc_sat))              DEALLOCATE(this%volc_sat)

  IF(ALLOCATED(this%art_cgaml_ptr))         DEALLOCATE(this%art_cgaml_ptr)
  IF(ALLOCATED(this%art_reac_rates_ptr))    DEALLOCATE(this%art_reac_rates_ptr)
  IF(ALLOCATED(this%art_NAT_rel_diff_ptr))  DEALLOCATE(this%art_NAT_rel_diff_ptr)
  IF(ALLOCATED(this%art_NAT_sedi_vel_ptr))  DEALLOCATE(this%art_NAT_sedi_vel_ptr)
  IF(ALLOCATED(this%art_photolysis_ptr))    DEALLOCATE(this%art_photolysis_ptr)

  IF(ALLOCATED(this%acc_drydepo_ptr))       DEALLOCATE(this%acc_drydepo_ptr)
  IF(ALLOCATED(this%acc_sedim_ptr))         DEALLOCATE(this%acc_sedim_ptr)
  IF(ALLOCATED(this%acc_wetdepo_gscp_ptr))  DEALLOCATE(this%acc_wetdepo_gscp_ptr)
  IF(ALLOCATED(this%acc_wetdepo_con_ptr))   DEALLOCATE(this%acc_wetdepo_con_ptr)
  IF(ALLOCATED(this%acc_wetdepo_rrsfc_ptr)) DEALLOCATE(this%acc_wetdepo_rrsfc_ptr)
  IF(ALLOCATED(this%emiss_ptr))             DEALLOCATE(this%emiss_ptr)
  IF(ALLOCATED(this%acc_emiss_ptr))         DEALLOCATE(this%acc_emiss_ptr)

END SUBROUTINE art_deallocate_diag

END MODULE mo_art_diag_types
