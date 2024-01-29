!
! mo_art_external_init_biomBurn
! This module provides initialization procedures for data from
! GFAS and for plume rise needed by the biomass burning emission routines
!
!
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

MODULE mo_art_external_init_biomBurn
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_parallel_config,               ONLY: nproma
! ART
  USE mo_art_external_types,            ONLY: t_art_biomBurn_properties
  USE mo_art_impl_constants,            ONLY: UNDEF_REAL_ART


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_extinit_biomBurn

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_biomBurn(nblks, p_biomBurn_prop)
!>
!! SUBROUTINE art_extinit_biomBurn
!! - 
!! - 
!! Based on Walter et al. - The importance of plume rise on the concentrations and
!!                          atmospheric impacts of biomass burning aerosol
!! Part of Module: mo_art_external_init_biomBurn
!! Author: Jonas Straub, KIT
!! Initial Release: 2018-02-12
!!
!! Modifications: 
!! YYYY-MM-DD: <name>, <institution>
!! - ...  
!<
 
! Within the first call of art_frp the lower and upper limits for the heat flux of 
! vegetation fires are attributed to the landuse classes in the domain. Freitas et al.
! (2006) (Impact of including the plume rise of vegetation fires in numerical simu-
! lations of associated atmospheric pollutants, Table 1) defined three biome types
! with corresponding lower and upper limits of heat flux. The variable frp_veg finally
! includes the fraction of these biome types (frp(:,:,1) - Tropical forest, frp(:,:,2) 
! - Woody savanna - cerrado, frp(:,:,3) - Grassland - pasture - cropland) to be able
! to use the suitable heat flux range for every gridpoint.
! Attention: The remote sensing frp obs defines the fire location! This landuse class 
! attribution only gives the link to the corresponding heat flux. 
! JS 2018: The flux of black carbon is proportionally to FRP, thus bc-flux is used for both
!          fire location and mass emission flux data.
! ======================================================================================== 
 
  TYPE(t_art_biomBurn_properties),INTENT(inout) :: &
    &  p_biomBurn_prop          !< temporal storage structure for emission calculation
  INTEGER,INTENT(in)       :: &
    &  nblks                    !< number of blocks
    
  ! ----------------------------------
  ! --- Allocate local storage
  ! ----------------------------------
  ! diurnal cycle for frp (hflux, burnt area and emissions)
  ALLOCATE( p_biomBurn_prop%dc_hflux_min_res (nproma,nblks,24) )
  ALLOCATE( p_biomBurn_prop%dc_hflux_max_res (nproma,nblks,24) )
  ALLOCATE( p_biomBurn_prop%dc_burnt_area_res(nproma,nblks,24) )
  ALLOCATE( p_biomBurn_prop%dc_emis_res      (nproma,nblks,24) )
  ALLOCATE( p_biomBurn_prop%flux_bc          (nproma,nblks)    )
 
  ! ----------------------------------
  ! initialize
  ! ----------------------------------
  p_biomBurn_prop%dc_hflux_min_res  = 0._wp
  p_biomBurn_prop%dc_hflux_max_res  = 0._wp
  p_biomBurn_prop%dc_burnt_area_res = 0._wp
  p_biomBurn_prop%dc_emis_res       = 0._wp
  p_biomBurn_prop%flux_bc           = UNDEF_REAL_ART
   
  
END SUBROUTINE art_extinit_biomBurn 
 
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_external_init_biomBurn
