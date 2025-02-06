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

! Cleanup/Destruction wrapper for NWP physics suite
!
! This routine destructs all variable lists related to the NWP physics suite,
! and calls deallocation routines for specific parameterizations,
! if necessary.

MODULE mo_nwp_phy_cleanup

  USE mo_nwp_phy_state,        ONLY: destruct_nwp_phy_state
  USE mo_nwp_lnd_state,        ONLY: destruct_nwp_lnd_state
  USE mo_nwp_reff_interface,   ONLY: reff_calc_dom
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_lnd_nwp_config,       ONLY: tile_list
  USE mo_grid_config,          ONLY: n_dom
  USE mo_aerosol_sources_types,ONLY: p_dust_source_const, p_fire_source_info
  USE mo_aerosol_util,         ONLY: tegen_scal_factors

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: cleanup_nwp_phy


CONTAINS

  !>
  !! Wrapper routine for NWP physics cleanup
  !!
  !! Performs destruction of NWP-specific variable lists and calls 
  !! deallocation routines of individual parameterizations (if necessary).
  !!
  SUBROUTINE cleanup_nwp_phy()

    ! local
    INTEGER :: jg
  !------------------------------------------------------

    ! destruct NWP physics and land variable lists
    !
    CALL destruct_nwp_phy_state()
    CALL destruct_nwp_lnd_state()

    ! destruct surface tile list
    CALL tile_list%destruct()

    DO jg = 1, n_dom
      IF ( atm_phy_nwp_config(jg)%icalc_reff > 0 ) CALL reff_calc_dom(jg)%destruct()
      !
      CALL atm_phy_nwp_config(jg)%finalize()
      
      IF ( iprog_aero > 0 ) CALL p_dust_source_const(jg)%finalize()
      IF ( iprog_aero > 2 ) CALL p_fire_source_info(jg)%finalize()
    ENDDO
    
    CALL tegen_scal_factors%finalize()

  END SUBROUTINE cleanup_nwp_phy

END MODULE mo_nwp_phy_cleanup

