!
! Data types for the two-moment bulk microphysics by Seifert and Beheng (2006)
!                  with prognostic cloud droplet number parameterization
!
! This module provides the data types for the variables used to
! configure the parameterization and to store the input and output
! for the parameterization.
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

MODULE mo_cloud_two_types

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cloud_two_config, t_cloud_two_input, t_cloud_two_output

  !-----------------------------------------------------------------------------
  
  TYPE t_cloud_two_config
     !
     ! configuration parameters
     ! ------------------------
     ! no parameters available.....
     !
     ! thresholds
     !
     ! grid scale microphysics
     !
  END TYPE t_cloud_two_config

  !-----------------------------------------------------------------------------
  
  TYPE t_cloud_two_input
     !
     ! Input arguments: parameters
     ! ---------------------------
     !
     INTEGER , POINTER :: jcs       (:,  :)=>NULL() !< column start index
     INTEGER , POINTER :: jce       (:,  :)=>NULL() !< column end   index
     INTEGER , POINTER :: msg_level (:,  :)=>NULL() !< message level
     REAL(wp), POINTER :: pdtime    (:,  :)=>NULL() !< [s] physics time step
     !
     ! Input arguments: fields
     ! -----------------------
     !
     ! atmospheric state
     REAL(wp), POINTER :: dz        (:,:,:)=>NULL() !< [m]      cell thickness
     REAL(wp), POINTER :: zh        (:,:,:)=>NULL() !< [m]      height of half levels
     REAL(wp), POINTER :: rho       (:,:,:)=>NULL() !< [kg/m3]  air density
     REAL(wp), POINTER :: pf        (:,:,:)=>NULL() !< [Pa]     air pressure
!inout
     REAL(wp), POINTER :: qv        (:,:,:)=>NULL() !< [kg/kg]  specific humidity
     REAL(wp), POINTER :: qc        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud water in air
     REAL(wp), POINTER :: qnc       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud droplet number
     REAL(wp), POINTER :: qr        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of rain        in air
     REAL(wp), POINTER :: qnr       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of rain droplet number
     REAL(wp), POINTER :: qi        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud ice   in air
     REAL(wp), POINTER :: qni       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud ice number
     REAL(wp), POINTER :: qs        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of snow        in air
     REAL(wp), POINTER :: qns       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of snow number
     REAL(wp), POINTER :: qg        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of graupel     in air
     REAL(wp), POINTER :: qng       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of graupel number
     REAL(wp), POINTER :: qh        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of hail        in air
     REAL(wp), POINTER :: qnh       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of hail number
     REAL(wp), POINTER :: ninact    (:,:,:)=>NULL() !< [kg/kg]  activated ice nuclei
     REAL(wp), POINTER :: ta        (:,:,:)=>NULL() !< [K]      air temperature
     REAL(wp), POINTER :: w         (:,:,:)=>NULL() !<          w
     ! fluxes at the surface
     REAL(wp), POINTER :: pr_rain     (:,  :)=>NULL() !< [kg/m2/s] sfc rain    flux
     REAL(wp), POINTER :: pr_ice      (:,  :)=>NULL() !< [kg/m2/s] sfc ice     flux
     REAL(wp), POINTER :: pr_snow     (:,  :)=>NULL() !< [kg/m2/s] sfc snow    flux
     REAL(wp), POINTER :: pr_grpl     (:,  :)=>NULL() !< [kg/m2/s] sfc graupel flux
     REAL(wp), POINTER :: pr_hail     (:,  :)=>NULL() !< [kg/m2/s] sfc hail    flux
     !
     !
  END TYPE t_cloud_two_input

  !-----------------------------------------------------------------------------
  
  TYPE t_cloud_two_output
     !
     ! Output arguments: fields
     ! ------------------------
     !
     ! tendencies in the atmosphere
     REAL(wp), POINTER :: tend_ta_two  (:,:,:)=>NULL() !< [K/s] tendency of temperature (cp)
     REAL(wp), POINTER :: tend_qv_two  (:,:,:)=>NULL() !< [1/s] tendency of specific humidity
     REAL(wp), POINTER :: tend_qc_two  (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of cloud water
     REAL(wp), POINTER :: tend_qnc_two (:,:,:)=>NULL() !< [1/s] tendency of number of cloud water
     REAL(wp), POINTER :: tend_qi_two  (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of cloud ice
     REAL(wp), POINTER :: tend_qni_two (:,:,:)=>NULL() !< [1/s] tendency of number of cloud ice
     REAL(wp), POINTER :: tend_qr_two  (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of rain
     REAL(wp), POINTER :: tend_qnr_two (:,:,:)=>NULL() !< [1/s] tendency of number of rain
     REAL(wp), POINTER :: tend_qs_two  (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of snow
     REAL(wp), POINTER :: tend_qns_two (:,:,:)=>NULL() !< [1/s] tendency of number of snow
     REAL(wp), POINTER :: tend_qg_two  (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of graupel
     REAL(wp), POINTER :: tend_qng_two (:,:,:)=>NULL() !< [1/s] tendency of number of graupel
     REAL(wp), POINTER :: tend_qh_two  (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of hail
     REAL(wp), POINTER :: tend_qnh_two (:,:,:)=>NULL() !< [1/s] tendency of number of hail
     REAL(wp), POINTER :: tend_ninact_two  (:,:,:)=>NULL() !< [1/s] tendency of activated ice nuclei

!inout
     REAL(wp), POINTER :: qv        (:,:,:)=>NULL() !< [kg/kg]  specific humidity
     REAL(wp), POINTER :: qc        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud water in air
     REAL(wp), POINTER :: qnc       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud droplet number
     REAL(wp), POINTER :: qr        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of rain        in air
     REAL(wp), POINTER :: qnr       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of rain droplet number
     REAL(wp), POINTER :: qi        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud ice   in air
     REAL(wp), POINTER :: qni       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud ice number
     REAL(wp), POINTER :: qs        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of snow        in air
     REAL(wp), POINTER :: qns       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of snow number
     REAL(wp), POINTER :: qg        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of graupel     in air
     REAL(wp), POINTER :: qng       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of graupel number
     REAL(wp), POINTER :: qh        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of hail        in air
     REAL(wp), POINTER :: qnh       (:,:,:)=>NULL() !< [kg/kg]  mass fraction of hail number
     REAL(wp), POINTER :: ninact    (:,:,:)=>NULL() !< [kg/kg]  activated ice nuclei
     REAL(wp), POINTER :: ta        (:,:,:)=>NULL() !< [K]      air temperature
     REAL(wp), POINTER :: w         (:,:,:)=>NULL() !<          w
     !
!inout
     ! fluxes at the surface
     REAL(wp), POINTER :: pr_rain     (:,  :)=>NULL() !< [kg/m2/s] sfc rain    flux
     REAL(wp), POINTER :: pr_ice      (:,  :)=>NULL() !< [kg/m2/s] sfc ice     flux
     REAL(wp), POINTER :: pr_snow     (:,  :)=>NULL() !< [kg/m2/s] sfc snow    flux
     REAL(wp), POINTER :: pr_grpl     (:,  :)=>NULL() !< [kg/m2/s] sfc graupel flux
     REAL(wp), POINTER :: pr_hail     (:,  :)=>NULL() !< [kg/m2/s] sfc hail    flux
     !
  END TYPE t_cloud_two_output

  !-----------------------------------------------------------------------------
  
END MODULE mo_cloud_two_types
