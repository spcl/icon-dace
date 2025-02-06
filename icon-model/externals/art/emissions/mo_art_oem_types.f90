!
! mo_art_oem_types
! This module provides datastructures for the
! online emisison module OEM
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

MODULE mo_art_oem_types

! ICON
  USE mo_kind,                          ONLY: wp
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: p_art_oem_data, &
            t_art_oem, &
            t_art_oem_data, &
            t_art_oem_config, &
            t_art_oem_ensemble

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem_data

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
      & country_ids                 ! EMEP country code for each domain (nproma,nblks)

    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
      & gridded_emissions,     & ! "2D" emissions fields (nproma,nblks,ncat)
      & tp_dayofweek,          & ! day-of-week scaling factor (ndays,ncat,ncountries)
      & tp_monthofyear,        & ! seasonal scaling factor
      & tp_hourofday,          & ! diurnal scaling factor
      & tp_hourofyear            ! hourly scaling factor

    REAL(KIND=wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      & lambda_mat               ! lambdas for ensembles

    REAL(KIND=wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      & vert_scaling_fact        ! vertical scale factors on model grid

    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
      & lswi,                  & ! LSWI field (jc,jb,nday)
      & evi,                   & ! EVI field (jc,jb,nday)
      & vprm_lu_class_fraction   ! VPRM land-use class fractions (jc,jb,num_vprm_lu_classes)

    REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
      & lswi_min,              & ! Minimum of LSWI field (jc,jb)
      & evi_min,               & ! Minimum of EVI field (jc,jb)
      & lswi_max,              & ! Minimum of LSWI field (jc,jb)
      & evi_max                  ! Minumum of EVI field (jc,jb)

    INTEGER :: &
      & i_vprm_lc_evergreen,   & !< VPRM land-use class 'evergreen'
      & i_vprm_lc_deciduous,   & !< VPRM land-use class 'deciduous'
      & i_vprm_lc_mixed,       & !< VPRM land-use class 'mixed forest'
      & i_vprm_lc_shrub,       & !< VPRM land-use class 'shrubland'
      & i_vprm_lc_savanna,     & !< VPRM land-use class 'savanna'
      & i_vprm_lc_crop,        & !< VPRM land-use class 'cropland'
      & i_vprm_lc_grass,       & !< VPRM land-use class 'grassland'
      & i_vprm_lc_urban          !< VPRM land-use class 'urban area'

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
      & reg_map                  ! region masks for grid cells for ensembles (nreg,nrpoma,nblks)

    INTEGER, DIMENSION(:), ALLOCATABLE :: &
      & tp_countryid             ! EMEP country code

  END TYPE t_art_oem_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem_config

    INTEGER :: tp_ncountry, & ! number of countries in dataset
      &        emis_tracer, & ! number of emisison tracer
      &        ens_tracer,  & ! number of ensemble tracer
      &        vprm_tracer    ! number of VPRM tracer

    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE   :: &
      & gridded_emissions_idx,  & ! name of the annual mean emissions fields
      & vp_category,            & ! category names of vertical profiles
      & tp_category,            & ! category names of temporal profiles
      & emis_name,              & ! name of the emission tracer
      & vprm_name,              & ! name of the VPRM tracer
      & vprm_flux_type            ! name of the VPRM flux type ('resp' or 'gpp')

    CHARACTER(LEN=20), DIMENSION(:,:), ALLOCATABLE :: &
      & ycatl_l, yvpl_l, ytpl_l   ! categories, vertical- and temporal profiles (ntracer,ncat)

    INTEGER, DIMENSION(:), ALLOCATABLE             :: &
      & emis_idx,               & ! indices of OEM-tracers with emissions within the ICON container
      & vprm_idx,               & ! indices of OEM-tracers with VPRM within the ICON container
      & itype_tscale_l            ! type of temporal scaling


  END TYPE t_art_oem_config

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem_ensemble

    CHARACTER(LEN=20), DIMENSION(200)  :: &
      & ens_name                  ! name of the reference tracer, where the ensemble belongs to

    INTEGER, DIMENSION(2,200)          :: &
      & ens_table                 ! indices of OEM-tracers with emissions within the ICON container


  END TYPE t_art_oem_ensemble

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem

    TYPE(t_art_oem_data)               :: data_fields
    TYPE(t_art_oem_config)             :: configure
    TYPE(t_art_oem_ensemble)           :: ensemble

  END TYPE t_art_oem

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  TYPE(t_art_oem),TARGET  :: &
    &  p_art_oem_data
 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


END MODULE mo_art_oem_types

