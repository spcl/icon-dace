!
! Data types for the cloud microphsics "graupel" parameterization
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

MODULE mo_cloud_mig_types

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_cloud_mig_config, t_cloud_mig_input, t_cloud_mig_output

  !-----------------------------------------------------------------------------
  
  TYPE t_cloud_mig_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! no parameters available.....
     !
     ! thresholds
     !
     ! grid scale microphysics
     !
  END TYPE t_cloud_mig_config

  !-----------------------------------------------------------------------------
  
  TYPE t_cloud_mig_input
     !
     ! Input arguments: parameters
     ! ---------------------------
     !
     INTEGER , POINTER :: jcs       (:,  :)=>NULL() !< column start index
     INTEGER , POINTER :: jce       (:,  :)=>NULL() !< column end   index
     REAL(wp), POINTER :: pdtime    (:,  :)=>NULL() !< [s] physics time step
     !
     ! Input arguments: fields
     ! -----------------------
     !
     ! atmospheric state
     REAL(wp), POINTER :: dz        (:,:,:)=>NULL() !< [m]      cell thickness
     REAL(wp), POINTER :: rho       (:,:,:)=>NULL() !< [kg/m3]  air density
     REAL(wp), POINTER :: pf        (:,:,:)=>NULL() !< [Pa]     air pressure
     REAL(wp), POINTER :: cpair     (:,:,:)=>NULL() !< [J/K/kg] isobaric specific heat of air
     REAL(wp), POINTER :: cvair     (:,:,:)=>NULL() !< [J/K/kg] isometric specific heat of air
     REAL(wp), POINTER :: ta        (:,:,:)=>NULL() !< [K]      air temperature
     REAL(wp), POINTER :: qv        (:,:,:)=>NULL() !< [kg/kg]  specific humidity
     REAL(wp), POINTER :: qc        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud water in air
     REAL(wp), POINTER :: qi        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of cloud ice   in air
     REAL(wp), POINTER :: qr        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of rain        in air
     REAL(wp), POINTER :: qs        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of snow        in air
     REAL(wp), POINTER :: qg        (:,:,:)=>NULL() !< [kg/kg]  mass fraction of graupel     in air
     !
  END TYPE t_cloud_mig_input

  !-----------------------------------------------------------------------------
  
  TYPE t_cloud_mig_output
     !
     ! Output arguments: fields
     ! ------------------------
     !
     ! tendencies in the atmosphere
     REAL(wp), POINTER :: tend_ta_mig (:,:,:)=>NULL() !< [K/s] tendency of air temperature (cp)
     REAL(wp), POINTER :: tend_qv_mig (:,:,:)=>NULL() !< [1/s] tendency of specific humidity
     REAL(wp), POINTER :: tend_qc_mig (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of cloud water in air
     REAL(wp), POINTER :: tend_qi_mig (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of cloud ice   in air
     REAL(wp), POINTER :: tend_qr_mig (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of rain        in air
     REAL(wp), POINTER :: tend_qs_mig (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of snow        in air
     REAL(wp), POINTER :: tend_qg_mig (:,:,:)=>NULL() !< [1/s] tendency of mass fraction of graupel     in air
     !
     ! fluxes at the surface
     REAL(wp), POINTER :: pr_rain     (:,  :)=>NULL() !< [kg/m2/s] sfc rain    flux
     REAL(wp), POINTER :: pr_ice      (:,  :)=>NULL() !< [kg/m2/s] sfc ice     flux
     REAL(wp), POINTER :: pr_snow     (:,  :)=>NULL() !< [kg/m2/s] sfc snow    flux
     REAL(wp), POINTER :: pr_grpl     (:,  :)=>NULL() !< [kg/m2/s] sfc graupel flux
     !
  END TYPE t_cloud_mig_output

  !-----------------------------------------------------------------------------
  
END MODULE mo_cloud_mig_types
