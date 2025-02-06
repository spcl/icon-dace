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

! This module provides derived types for mo_aerosol_source which is used
! to calculate source terms for different aerosol species.
!
! Implemented up to now:
! - Mineral dust look-up tables for the clay content of different soil types
!   and the bare soil fraction of different land use classes

MODULE mo_aerosol_sources_types

  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_impl_constants,                ONLY: max_dom
  USE mo_io_units,                      ONLY: filename_max

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_dust_source_const, p_dust_source_const
  PUBLIC :: t_fire_source_info, p_fire_source_info

  ! Type to hold constant information needed for the mineral dust source
  TYPE t_dust_source_const
    LOGICAL                 :: &
      &  is_initialized=.FALSE.  !< Has this object been initialized?
    REAL(wp), ALLOCATABLE   :: &
      &  f_bare(:),            & !< Bare soil fraction, deducted from several land use classes (-)
      &  f_clay(:)               !< Clay fraction (-)
    ! ---------------------------------------------------------------------------------------------
    CONTAINS
      PROCEDURE :: init     => init_aerosol_dust_aod_source
      PROCEDURE :: finalize => finalize_aerosol_dust_aod_source
  END TYPE t_dust_source_const

  TYPE t_fire_species_info
    REAL(wp), POINTER            :: var(:,:)
    CHARACTER(LEN=10)            :: varname
    CHARACTER(LEN=filename_max)  :: current_filename
  END TYPE t_fire_species_info

  TYPE t_fire_source_info
    LOGICAL                   :: &
      &  is_initialized=.FALSE.    !< Has this object been initialized?
    INTEGER                   :: &
      &  nspecies                  !< Species: SO2, BC, OC
    TYPE(t_fire_species_info), ALLOCATABLE :: &
      &  species(:)
    ! ---------------------------------------------------------------------------------------------
    CONTAINS
      PROCEDURE :: init     => init_aerosol_fire_source_info
      PROCEDURE :: finalize => finalize_aerosol_fire_source_info
  END TYPE t_fire_source_info

  TYPE(t_dust_source_const) :: p_dust_source_const(max_dom)
  TYPE(t_fire_source_info)  :: p_fire_source_info(max_dom)

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_aerosol_sources_types'

CONTAINS

  !>
  !! SUBROUTINE init_aerosol_dust_aod_source
  !!
  !! * Allocates and initializes several arrays needed for mineral dust emission
  !!   which are constant in time.
  !! * Provide emission scaling factors for different soil types and land use classes
  !!
  SUBROUTINE init_aerosol_dust_aod_source (this_source,                             &
    &                                      i_lc_shrub_eg,  i_lc_shrub,  i_lc_grass, & 
    &                                      i_lc_bare_soil, i_lc_sparse,             &
    &                                      i_st_sand, i_st_sandyloam, i_st_loam,    &
    &                                      i_st_clayloam, i_st_clay,                &
    &                                      nlu_classes, soiltype_sidx, soiltype_eidx)

    CLASS(t_dust_source_const), INTENT(inout) :: &
      &  this_source
    INTEGER, INTENT(in) :: &
      &  i_lc_shrub_eg,    & !< Land use class index for shrub cover grassland/forest evergreen
      &  i_lc_shrub,       & !< Land use class index for closed to open shrubland (deciduous)
      &  i_lc_grass,       & !< Land use class index for grassland/herbaceous
      &  i_lc_bare_soil,   & !< Land use class index for bare soil
      &  i_lc_sparse,      & !< Land use class index for sparse vegetation
      &  i_st_sand,        & !< Soil type index for sand
      &  i_st_sandyloam,   & !< Soil type index for sandy loam
      &  i_st_loam,        & !< Soil type index for loam
      &  i_st_clayloam,    & !< Soil type index for clay loam
      &  i_st_clay,        & !< Soil type index for clay
      &  nlu_classes,      & !< Number of land use classes
      &  soiltype_sidx,    & !< Soil type array start index
      &  soiltype_eidx       !< Soil type array end index

    ! Local variables
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':init_aerosol_dust_aod_source'

    IF (this_source%is_initialized) THEN
      CALL finish(TRIM(routine),'this_source is already initialized.')
    ELSE

      IF(.NOT. ALLOCATED(this_source%f_bare)) THEN
        ALLOCATE(this_source%f_bare(nlu_classes) )
      ELSE
        CALL finish(TRIM(routine),'this_source%f_clay already allocated.')
      ENDIF

      IF(.NOT. ALLOCATED(this_source%f_clay)) THEN
        ALLOCATE(this_source%f_clay(soiltype_sidx:soiltype_eidx) )
      ELSE
        CALL finish(TRIM(routine),'this_source%f_clay already allocated.')
      ENDIF

      ! First set f_bare for all land use classes to 0. Then set it to specific values
      ! for some land use classes. These are tuning parameters, the emission at specific land use
      ! tiles can be linearly scaled with these parameters
      this_source%f_bare(:)              = 0._wp
      this_source%f_bare(i_lc_shrub_eg)  = 0._wp
      this_source%f_bare(i_lc_shrub)     = 0.25_wp
      this_source%f_bare(i_lc_grass)     = 0._wp
      this_source%f_bare(i_lc_bare_soil) = 1._wp
      this_source%f_bare(i_lc_sparse)    = 0.5_wp

      ! First set f_clay for all soil types to 0. Then set it to specific values for specfic 
      ! soil types. The values are chosen from the USDA soil texture triangle. The emissions
      ! for specific soil types can be linearly scaled with this parameter
      this_source%f_clay(:)              = 0._wp
      this_source%f_clay(i_st_sand)      = 0.075_wp  ! Terra: 5%
      this_source%f_clay(i_st_sandyloam) = 0.10_wp   ! Terra: 10%
      this_source%f_clay(i_st_loam)      = 0.125_wp  ! Terra: 20%
      this_source%f_clay(i_st_clayloam)  = 0.275_wp  ! Terra: 35%
      this_source%f_clay(i_st_clay)      = 0.40_wp   ! Terra: 70%

      ! Finally, set initialized flag to true
      this_source%is_initialized = .TRUE.

    ENDIF !is_initialized
  END SUBROUTINE init_aerosol_dust_aod_source

  !>
  !! SUBROUTINE finalize_aerosol_dust_aod_source
  !!
  !! * Cleanup for the dust source constant fields
  !!
  SUBROUTINE finalize_aerosol_dust_aod_source(this_source)

    CLASS(t_dust_source_const), INTENT(inout) :: &
      &  this_source

    IF(ALLOCATED(this_source%f_bare)) DEALLOCATE(this_source%f_bare)
    IF(ALLOCATED(this_source%f_clay)) DEALLOCATE(this_source%f_clay)
    this_source%is_initialized = .FALSE.

  END SUBROUTINE finalize_aerosol_dust_aod_source

  !>
  !! SUBROUTINE init_aerosol_fire_source_info
  !!
  !! * Initializes the t_fire_source_info derived type. 
  !! * Sets pointer to ICON memory where the wildfire data is stored
  !!
  SUBROUTINE init_aerosol_fire_source_info(this_info, bcfire, ocfire, so2fire)

    CLASS(t_fire_source_info), INTENT(inout)  :: &
      &  this_info
    REAL(wp), POINTER, INTENT(in) :: &
      & bcfire(:,:), ocfire(:,:), so2fire(:,:)
    ! Local variables
    INTEGER :: &
      &  js

    this_info%nspecies = 3
    ALLOCATE(this_info%species(this_info%nspecies))

    this_info%species(1)%var     => bcfire
    this_info%species(1)%varname =  'bcfire'
    this_info%species(2)%var     => ocfire
    this_info%species(2)%varname =  'ocfire'
    this_info%species(3)%var     => so2fire
    this_info%species(3)%varname =  'so2fire'

    DO js = 1, this_info%nspecies
      this_info%species(js)%current_filename =  ''
    ENDDO

    this_info%is_initialized = .TRUE.

  END SUBROUTINE init_aerosol_fire_source_info

  !>
  !! SUBROUTINE finalize_aerosol_fire_source_info
  !!
  !! * Cleanup for the fire source info type
  !!
  SUBROUTINE finalize_aerosol_fire_source_info(this_info)

    CLASS(t_fire_source_info), INTENT(inout)  :: &
      &  this_info
    ! Local variables
    INTEGER :: &
      &  js

    DO js = 1, this_info%nspecies
      this_info%species(js)%var              => NULL()
      this_info%species(js)%varname          =  ''
      this_info%species(js)%current_filename =  ''
    ENDDO

    DEALLOCATE(this_info%species)

    this_info%nspecies       = 0
    this_info%is_initialized = .FALSE.

  END SUBROUTINE finalize_aerosol_fire_source_info

END MODULE mo_aerosol_sources_types
