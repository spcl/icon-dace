!
! mo_art_data
! This module provides data storage structures and constants.
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

MODULE mo_art_data
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_diag_types,                ONLY: t_art_diag
  USE mo_art_external_types,            ONLY: t_art_external
  USE mo_art_emiss_types,               ONLY: t_art_aero_emiss,  &
                                          &   t_art_emiss_storage
  USE mo_art_pntSrc_types,              ONLY: t_art_all_pntSrc
  USE mo_art_prescribed_types,          ONLY: t_art_prescr_list
  USE mo_art_fplume_types,              ONLY: t_art_all_volc_fplume
  USE mo_art_chem_data,                 ONLY: t_art_chem
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_emission_pollen_atab,      ONLY: t_art_all_stns
  
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: akt, rakt
  PUBLIC :: p_art_data, t_art_data
  PUBLIC :: t_art_air_properties
  PUBLIC :: t_art_turb_fields
  PUBLIC :: t_art_tend
  PUBLIC :: t_art_aero

  REAL(wp), PARAMETER ::               &
    &  akt      =  0.4_wp,             &  !< von Karman-constant
    &  rakt     =  (1._wp / akt)          !< 1 / von Karman-constant


  ! ----------------------------------
  ! --- Type structure to store parameters characterizing the state of the air
  ! ----------------------------------
  TYPE t_art_air_properties
    REAL(wp),POINTER      ::  &
      &  art_free_path(:,:,:)   !< Mean free path in the air
    REAL(wp),POINTER      ::  &
      &  art_dyn_visc(:,:,:)    !< Dynamic viscosity
  END TYPE t_art_air_properties
  
  
  ! ----------------------------------
  ! --- Type structure to store fields in a suitable way for the turbulence routine
  ! ---  as the turbulence routine requires consecutive fields which have the same 
  ! ---  structure as the atmospheric values
  ! ----------------------------------
  TYPE t_art_turb_fields
    REAL(wp),POINTER      ::  &
      &  sv(:,:,:)              !< Surface values of species with deposition 
                                !   (nproma,nblocks,ntracer_turb)
    REAL(wp),POINTER      ::  &
      &  vdep(:,:,:)             !< deposition velocity (nproma,nblocks,ntracer_turb)
  END TYPE t_art_turb_fields

  ! ----------------------------------
  ! Structure for ART tendencies
  ! ----------------------------------


  TYPE t_art_tend
    
    REAL(wp), POINTER, DIMENSION (:,:,:,:)     :: chem

  END TYPE t_art_tend


  ! ----------------------------------
  ! --- General ART data structure
  ! ----------------------------------



  TYPE t_art_aero
    INTEGER ::                 &
      &  itr_start,            & !< first index of aerosol tracer in tracer container
      &  itr_end                 !< last  index of aerosol tracer in tracer container
  END TYPE


  TYPE t_art_data
    TYPE(t_art_air_properties)          :: air_prop
    TYPE(t_art_turb_fields)             :: turb_fields
    TYPE(t_art_diag)                    :: diag
    TYPE(t_art_external)                :: ext
    TYPE(t_art_chem)                    :: chem
    TYPE(t_art_emiss_storage)           :: emiss
    TYPE(t_art_aero_emiss)              :: tracer2aeroemiss
    TYPE(t_art_prescr_list)             :: prescr_list
    TYPE(t_art_all_pntSrc)              :: pntSrc
    TYPE(t_art_all_volc_fplume)         :: fplume_init_all
    TYPE(t_key_value_store)             :: dict_tracer
    TYPE(t_art_tend)                    :: tend
    TYPE(t_art_atmo)                    :: atmo
    TYPE(t_art_aero)                    :: aero
  END TYPE t_art_data

  TYPE(t_art_data),TARGET, ALLOCATABLE  :: &
    &  p_art_data(:)            !< dim[n_dom]

END MODULE mo_art_data
