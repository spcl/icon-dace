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

MODULE mo_ocean_hamocc_communication






  USE iso_c_binding,               ONLY: c_ptr
  USE mo_exception,                ONLY: message, finish

!  USE mo_kind,                     ONLY: wp
  USE mo_master_control,           ONLY: ocean_process, hamocc_process
  USE mo_parallel_config,          ONLY: nproma
  USE mo_model_domain,             ONLY: t_patch
  USE mo_grid_config,              ONLY: n_dom, n_dom_start
  USE mo_run_config,               ONLY: num_lev



  IMPLICIT NONE

  PRIVATE


  PUBLIC :: setup_ocean_2_hamocc_communication, &
            setup_hamocc_2_ocean_communication, &
            exchange_data_ocean_2_hamocc, &
            exchange_data_hamocc_2_ocean, &
            free_ocean_hamocc_communication

CONTAINS



  SUBROUTINE setup_ocean_2_hamocc_communication(p_patch, no_of_levels)
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(in) :: no_of_levels


  END SUBROUTINE setup_ocean_2_hamocc_communication

  SUBROUTINE setup_hamocc_2_ocean_communication(p_patch, no_of_levels)
    TYPE(t_patch), INTENT(IN)       :: p_patch
    INTEGER, INTENT(in) :: no_of_levels


  END SUBROUTINE setup_hamocc_2_ocean_communication

  SUBROUTINE free_ocean_hamocc_communication
 

  END SUBROUTINE free_ocean_hamocc_communication

  SUBROUTINE exchange_data_ocean_2_hamocc( &
    &  top_dilution_coeff, &
    &  h_old, &
    &  h_new, &
    &  h_old_withIce, &
    &  ice_concentration_sum, &
    &  temperature, &
    &  salinity, &
    &  ver_diffusion_coeff, &
    &  short_wave_flux, &
    &  wind10m, &
    &  co2_mixing_ratio, &
    &  mass_flux_e, &
    &  vn, &
    &  w, &
    &  press_hyd, &
    &  stretch_c, &
    &  stretch_c_new, &
    &  draftave)


    !INTEGER,INTENT(INOUT)  :: &
    TYPE(c_ptr), INTENT(IN) :: &
      &  top_dilution_coeff, &
      &  h_old, &
      &  h_new, &
      &  h_old_withIce, &
      &  ice_concentration_sum, &
      &  temperature, &
      &  salinity, &
      &  ver_diffusion_coeff, &
      &  short_wave_flux, &
      &  wind10m, &
      &  co2_mixing_ratio, &
      &  mass_flux_e, &
      &  vn, &
      &  w , & 
      &  press_hyd, &
      &  stretch_c, &
      &  stretch_c_new, &
      &  draftave
 

    TYPE(c_ptr) :: src_data_cptr(18), dst_data_cptr(18)

    src_data_cptr( 1) =  top_dilution_coeff
    src_data_cptr( 2) =  h_old
    src_data_cptr( 3) =  h_new
    src_data_cptr( 4) =  h_old_withIce
    src_data_cptr( 5) =  ice_concentration_sum
    src_data_cptr( 6) =  temperature
    src_data_cptr( 7) =  salinity
    src_data_cptr( 8) =  ver_diffusion_coeff
    src_data_cptr( 9) =  short_wave_flux
    src_data_cptr(10) =  wind10m
    src_data_cptr(11) =  co2_mixing_ratio
    src_data_cptr(12) =  mass_flux_e
    src_data_cptr(13) =  vn
    src_data_cptr(14) =  w
    src_data_cptr(15) =  press_hyd
    src_data_cptr(16) =  stretch_c
    src_data_cptr(17) =  stretch_c_new
    src_data_cptr(18) =  draftave

    dst_data_cptr = src_data_cptr

    CALL finish("exchange_data_ocean_2_hamocc", " Requires the YAXT library")


  END SUBROUTINE exchange_data_ocean_2_hamocc


  SUBROUTINE exchange_data_hamocc_2_ocean(co2_flux, swr_fraction)

    TYPE(c_ptr), INTENT(IN)  :: co2_flux
    TYPE(c_ptr), INTENT(IN)  :: swr_fraction

    TYPE(c_ptr) :: src_data_cptr(2), dst_data_cptr(2)
   
    src_data_cptr( 1) = co2_flux
    src_data_cptr( 2) = swr_fraction
    dst_data_cptr = src_data_cptr

    CALL finish("exchange_data_ocean_2_hamocc", " Requires the YAXT library")

  END SUBROUTINE exchange_data_hamocc_2_ocean

END MODULE mo_ocean_hamocc_communication
