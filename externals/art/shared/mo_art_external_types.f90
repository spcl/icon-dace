!
! mo_art_external_types
! This module provides data storage structures for data from
! external sources (e.g. soil types, volcano characteristics,
! radioactive source characteristics, etc.)
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

MODULE mo_art_external_types
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
! ART
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_art_emission_pollen_atab,      ONLY: t_art_pol_atab

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_art_external
  PUBLIC :: t_art_soil_table, t_art_soil_properties, t_art_land_properties
  PUBLIC :: t_art_volc_table, t_art_volcdata, t_art_volc_fplume
  PUBLIC :: t_art_pollen_table, t_art_pollen_properties
  PUBLIC :: t_art_online_dms  
  PUBLIC :: t_art_biomBurn_properties

  ! ----------------------------------
  ! --- Type structure to store soil parameters needed for mineral dust emission schemes
  ! ----------------------------------
  TYPE t_art_soil_table
    CHARACTER(LEN=4)      ::  &
      &  sname                   !< short name of this soil type
    INTEGER               ::  &
      &  nsoil_modes             !< number of soil modes, 4 for shao
    REAL(wp),POINTER      ::  &
      &  diam_med(:,:),       &  !< Median diameter, dim= (nsoilmodes,2) 1 = minimal dispersed,
                                 !      2= fully dispersed [m]
      &  std_dev(:,:),        &  !< standard deviation, dim= (nsoilmodes,2) 1 = minimal dispersed,
                                 !      2= fully dispersed
      &  fr_mode(:,:)            !< massfraction of this mode of total modal distribution,
                                 !    dim= (nsoilmodes,2) 1 = minimal dispersed, 2= fully dispersed
    REAL(wp),POINTER      ::  &
      &  fr_soil(:,:)            !< Fraction of soil type in grid box (dims: nproma, nblocks)
    REAL(wp)              ::  &
      &  stot_min,            &  !< total cross sectional area of soil particles for min. disp.
      &  stot_ful,            &  !< total cross sectional area of soil particles for full. disp.
      &  z0s                     !< total roughness length for smooth conditions [m]
  END TYPE t_art_soil_table
  
  ! ----------------------------------
  ! --- Type structure to store soil parameters needed for mineral dust emission schemes
  ! ----------------------------------
  TYPE t_art_soil_properties
    INTEGER               ::  &
      &  nsoil_types             !< Number of different soil types (=14 for HWSD data)
    REAL(wp),POINTER      ::  &
      &  f_z0(:,:)               !< impact of varying soil roughness length (dims: nproma, nblocks)
    REAL(wp),POINTER      ::  &
      &  wstrich(:,:)            !< Minimal value of gravimetric water content
    REAL(wp),POINTER      ::  & 
      &  emiss_rate_a(:,:),   &  !< emission rate mode a
      &  emiss_rate_b(:,:),   &  !< emission rate mode b
      &  emiss_rate_c(:,:)       !< emission rate mode c
    LOGICAL,POINTER       ::  &
      & dust_mask(:,:)           !< dust emission mask; indicates whether dust emissions take place
    TYPE(t_art_soil_table),POINTER :: &
      &  soil_type(:)
  END TYPE t_art_soil_properties

  ! ----------------------------------
  ! --- Type structure to store land parameters
  ! ----------------------------------
  TYPE t_art_land_properties
    REAL(wp),POINTER      ::  &
      pft(:,:,:) => NULL()       !< Plant functional type for online biogenic emission
  END TYPE t_art_land_properties

  ! ----------------------------------
  ! --- Volcanic ash
  ! ----------------------------------
  TYPE t_art_volc_table
    CHARACTER(LEN=8)  ::             id         !< = RN//SN//VN
    CHARACTER(LEN=2)  ::             rn
    CHARACTER(LEN=2)  ::             sn
    CHARACTER(LEN=4)  ::             vn
    CHARACTER(LEN=40) ::             name       !< name of volcano
    CHARACTER(LEN=40) ::             location   !< name of country
    CHARACTER(LEN=20) ::             status     !< status of volcano: Historical (= documented
                                                !   proof of eruption), Holocene, ...
    REAL(wp)          ::             latitude   !< latitude of volcano
    CHARACTER(LEN=1)  ::             ns         !< north or south hemisphere
    CHARACTER(LEN=1)  ::             vf         !< ?
    REAL(wp)          ::             longitude  !< longitude of volcano
    CHARACTER(LEN=1)  ::             ew         !< east or west
    REAL(wp)          ::             elev       !< elevation of volcano
    CHARACTER(LEN=20) ::             type       !< type of volcano: Stratovolcano,Shield_volcano,..
    CHARACTER(LEN=2)  ::             timeframe  !< ?
    CHARACTER(LEN=2)  ::             eruption_type !< eruption type defining source parameters
    REAL(wp) ::                      height_above_vent !< [m]
    REAL(wp) ::                      erup_duration     !< [h]
    REAL(wp) ::                      eruption_rate     !< [kg/s]
    REAL(wp) ::                      eruption_volume   !< [km^3]
    REAL(wp) ::                      ash_fraction  !< [] mass fraction of ash smaller than 63 mu m
    REAL(wp), ALLOCATABLE ::         cfactor(:) !< distribution factors for individual size classes
    REAL(wp) ::                      mastin_fac !< pre-factor in Mastin et al. 2009 eruption rate /
                                                !     height above vent relation
    REAL(wp) ::                      mastin_exp !< exponent   in Mastin et al. 2009 eruption rate /
                                                !     height above vent relation
    INTEGER  ::                      tri_iidx_loc  !< info about triangle:index
    INTEGER  ::                      tri_iblk_loc  !< info about triangle:blk
  END TYPE t_art_volc_table
  
  TYPE t_art_volcdata
    TYPE(t_art_volc_table),POINTER :: &
      &  info(:)               !< Table with information about current volcano
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: &
      &  volcanofile_path      !< Path of volcanofile
    INTEGER, POINTER      :: &
      &  owner(:)              !< global_id of local points
    REAL(wp), POINTER     :: &
      &  height_factor(:,:), & !< factor for each layer
      &  lat_idx(:,:),       & !< factor for each layer
      &  lon_idx(:,:)          !< factor for each layer
    INTEGER               :: &
      &  ithis_nlocal_pts      !< Local (PE) volcanoes
    INTEGER               :: &
      &  maxvolcs=1000         !< global number of volcanoes
    LOGICAL               :: &
      &  lprint_once = .TRUE.  !< flag to print some messages only once
  END TYPE t_art_volcdata

  ! ----------------------------------
  ! --- Volcanic ash from FPLUME
  ! ----------------------------------
  TYPE t_art_volc_fplume
    REAL(wp), POINTER      :: & 
      & MER_transport(:,:),           &  !< fraction of total MER with particle smaller than 32mu
      & plume_H(:,:),                   &
      & plume_MER(:,:)
  END TYPE t_art_volc_fplume

  ! ----------------------------------
  ! --- Pollen type structure for plant specific parameters
  ! ----------------------------------
  TYPE t_art_pollen_table
    LOGICAL               :: &
      &  linit = .FALSE.             !< set to .TRUE. if pollen type was successfully initialized
    CHARACTER(LEN=10)       ::  &
      &  sname                       !< short name of this pollen type
    CHARACTER(LEN=4)        ::  &
      &  shortname                   !< short name of this pollen type for grib output
    INTEGER                 ::  &
      &  jul_days_excl               !< pollen-specific start day for temperature sum minus 1         
    REAL(wp)                ::  &
      &  psi_random,            &    !< Loss of pollen due to random processes
      &  xi_r_precip,           &    !< Maximum capacity of the precipitation reservoir
      &  frac_xi_evap,          &    !< Fraction of xi_r_precip evaporating per timestep
      &  t_base                      !< pollen-specific base temperature for temperature sum
    REAL(wp),POINTER            ::  &
      &  fr_cov(:,:),           &    !< fraction of plant coverage in grid box/
                                     !  (dims: nproma, nblocks)      
      &  res_new_sum(:,:),      &    !< Daily number of pollen released into the reservoir
      &  res_old(:,:),          &    !< Number of pollen in the reservoir (previous timestep)
      &  f_q_alt(:,:),          &    !< Height correction for emission. Decreasing emission with
                                     !  height
      &  r_precip(:,:),         &    !< Precipitation reservoir
      &  ctsum(:,:),            &    !< cumulated temperature sum
      &  saisn(:,:),            &    !< number of days since the start of the pollen
                                     !  season if present day is during the season; zero outside 
                                     !  the season
      &  saisa(:,:),            &    !< as saisn, but contains length of season if after the 
                                     !  season and not zeros (contains all days of the season)
      &  saisl(:,:),            &    !< length of pollen seasons
      &  tthrs(:,:),            &    !< cumulated temperature sum threshold for the start of the
                                     !  pollen season
      &  tthre(:,:),            &    !< cumulated temperature sum threshold for the end of the
                                     !  pollen season
      &  tthrs_red(:,:),        &    !< reduction of the threshold to account for the fact that
                                     !  there are already pollen emitted before Pollen>30
                                     !  (dims:nproma, nblocks)									 
      &  f_q_seas(:,:),         &    !< state of pollen season
      &  fe_plant(:,:),         &    !< emission flux of pollen
      &  sobs_sum(:,:),         &    !< sum of radiation since midnight for ambrosia emissions
      &  rh_sum(:,:),           &    !< sum of relative humidity since midnight for ambrosia emissions
      &  tune(:,:),             &    !< tuning factor for pollen emission
      &  no_max_day(:,:),       &    !< Maximum number of pollen that can be released per day
                                     !  (per m2)
      &  no_max_timestep(:,:)        !< Maximum number of pollen that can be released per timestep
                                     !  (per m2)

    TYPE(t_art_pol_atab)        &
      &  pol_atab                    !< information read out from atab-file
  END TYPE t_art_pollen_table

  ! ----------------------------------
  ! --- Type structure to store ALL pollen parameters while emission
  ! ----------------------------------
  TYPE t_art_pollen_properties
    INTEGER                 ::  &
      &  npollen_types = 4,     &      !< Number of different pollen types (default = 1, birch)
                                       !  further Develpment: set this number dynamically
      &  npollen_used  = 0             !< actual number of different pollen types used
    TYPE(t_key_value_store) ::  &
      &  dict_pollen
  
    TYPE(t_art_pollen_table),POINTER :: &
      &  pollen_type(:)            !< store specific information for each pollen type
  END TYPE t_art_pollen_properties

  ! ----------------------------------
  ! --- Online DMS
  ! ----------------------------------
  TYPE t_art_online_dms
    INTEGER                 ::  ndms_months = 12       !< Number of months
    REAL(wp), ALLOCATABLE   ::  dms_month(:,:,:)       !< Data for each month
    REAL(wp), ALLOCATABLE   ::  mmr_onl_dms(:)         !< mass mixing ratio added to the tracer
    LOGICAL                 ::  lcalc_onl = .FALSE.    !< switch for online calculation 
  END TYPE t_art_online_dms

  
  ! ----------------------------------
  ! --- Biomass burning type structure 
  ! ----------------------------------  
  TYPE t_art_biomBurn_properties
    ! FRP diurnal cycle basic profiles
    REAL(wp),POINTER        :: &
      &  dc_hflux_min_res(:,:,:),  &  !< diurnal cycle for heat flux min
      &  dc_hflux_max_res(:,:,:),  &  !< diurnal cycle for heat flux max 
      &  dc_burnt_area_res(:,:,:), &  !< diurnal cycle for burnt area
      &  dc_emis_res(:,:,:),       &  !< (basic) diurnal cycle for emissions
                                      !  (multiplied later with GFAS FRP emissions)
      &  flux_bc(:,:)                 !< flux of black carbon (kg /s m2), GFAS product                                
  
  END TYPE t_art_biomBurn_properties
  
  
  
  ! Type containing all external data
  TYPE t_art_external
    TYPE(t_art_soil_properties)         :: soil_prop
    TYPE(t_art_land_properties)         :: land
    TYPE(t_art_volcdata)                :: volc_data
    TYPE(t_art_volc_fplume)             :: volc_fplume
    TYPE(t_art_pollen_properties)       :: pollen_prop
    TYPE(t_art_online_dms)              :: online_dms
    TYPE(t_art_biomBurn_properties)     :: biomBurn_prop
  END TYPE t_art_external

END MODULE mo_art_external_types
