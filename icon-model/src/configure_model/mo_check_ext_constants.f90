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

! This module checks whether the constants once copied from `mo_impl_constants`
! to `mo_gridman_constants` located in the external math-library sill match
! their original values. If there is any change in the original value, the
! subroutine here will send an error message and finish the calculation

MODULE mo_check_ext_constants


  USE mo_gridman_constants, ONLY: modelname_ext => modelname
  USE mo_gridman_constants, ONLY: modelversion_ext => modelversion
  USE mo_gridman_constants, ONLY: MAX_CHAR_LENGTH_ext => MAX_CHAR_LENGTH
  USE mo_gridman_constants, ONLY: SUCCESS_ext => SUCCESS
  USE mo_gridman_constants, ONLY: CELLS_ext => CELLS
  USE mo_gridman_constants, ONLY: EDGES_ext => EDGES
  USE mo_gridman_constants, ONLY: VERTS_ext => VERTS
  USE mo_gridman_constants, ONLY: ON_CELLS_ext => ON_CELLS
  USE mo_gridman_constants, ONLY: ON_EDGES_ext => ON_EDGES
  USE mo_gridman_constants, ONLY: ON_VERTICES_ext => ON_VERTICES
  USE mo_gridman_constants, ONLY: HALO_LEVELS_CEILING_ext => HALO_LEVELS_CEILING
  USE mo_gridman_constants, ONLY: grf_bdywidth_c_ext => grf_bdywidth_c
  USE mo_gridman_constants, ONLY: grf_bdywidth_e_ext => grf_bdywidth_e
  USE mo_gridman_constants, ONLY: max_hw_ext => max_hw
  USE mo_gridman_constants, ONLY: min_rlcell_int_ext => min_rlcell_int
  USE mo_gridman_constants, ONLY: min_rlcell_ext => min_rlcell
  USE mo_gridman_constants, ONLY: max_rlcell_ext => max_rlcell
  USE mo_gridman_constants, ONLY: min_rlvert_int_ext => min_rlvert_int
  USE mo_gridman_constants, ONLY: min_rlvert_ext => min_rlvert
  USE mo_gridman_constants, ONLY: max_rlvert_ext => max_rlvert
  USE mo_gridman_constants, ONLY: min_rledge_int_ext => min_rledge_int
  USE mo_gridman_constants, ONLY: min_rledge_ext => min_rledge
  USE mo_gridman_constants, ONLY: max_rledge_ext => max_rledge
  USE mo_gridman_constants, ONLY: start_prog_cells_ext => start_prog_cells
  USE mo_gridman_constants, ONLY: start_prog_edges_ext => start_prog_edges
  USE mo_gridman_constants, ONLY: end_prog_cells_ext => end_prog_cells
  USE mo_gridman_constants, ONLY: end_halo_lev1_cells_ext => end_halo_lev1_cells
  USE mo_gridman_constants, ONLY: end_all_cells_ext => end_all_cells
  USE mo_gridman_constants, ONLY: end_prog_edges_ext => end_prog_edges
  USE mo_gridman_constants, ONLY: end_edges_of_prog_cells_ext => end_edges_of_prog_cells
  USE mo_gridman_constants, ONLY: end_halo_lev1_edges_ext => end_halo_lev1_edges
  USE mo_gridman_constants, ONLY: end_all_edges_ext => end_all_edges
  USE mo_gridman_constants, ONLY: end_prog_verts_ext => end_prog_verts
  USE mo_gridman_constants, ONLY: end_verts_of_prog_cells_ext => end_verts_of_prog_cells
  USE mo_gridman_constants, ONLY: end_all_verts_ext => end_all_verts
  USE mo_gridman_constants, ONLY: max_dom_ext => max_dom
  USE mo_gridman_constants, ONLY: max_dom_dig10_ext => max_dom_dig10
  USE mo_gridman_constants, ONLY: max_phys_dom_ext => max_phys_dom
  USE mo_gridman_constants, ONLY: max_ntracer_ext => max_ntracer
  USE mo_gridman_constants, ONLY: max_echotop_ext => max_echotop
  USE mo_gridman_constants, ONLY: max_wshear_ext => max_wshear
  USE mo_gridman_constants, ONLY: max_srh_ext => max_srh
  USE mo_gridman_constants, ONLY: TORUS_MAX_LAT_ext => TORUS_MAX_LAT

  USE mo_impl_constants,    ONLY: modelname, modelversion, MAX_CHAR_LENGTH, SUCCESS, CELLS, EDGES, &
                                  VERTS, ON_CELLS, ON_EDGES, ON_VERTICES, HALO_LEVELS_CEILING, max_hw, &
                                  min_rlcell_int, min_rlcell, max_rlcell, min_rlvert_int, min_rlvert, max_rlvert,&
                                  min_rledge_int, min_rledge, max_rledge, start_prog_cells, start_prog_edges, &
                                  end_prog_cells, end_halo_lev1_cells, end_all_cells, end_prog_edges, &
                                  end_edges_of_prog_cells, end_halo_lev1_edges, end_all_edges, end_prog_verts, &
                                  end_verts_of_prog_cells, end_all_verts, max_dom, max_dom_dig10, max_phys_dom, &
                                  max_ntracer, max_echotop, max_wshear, max_srh, TORUS_MAX_LAT

  USE mo_impl_constants_grf,ONLY: grf_bdywidth_c, grf_bdywidth_e

  USE mo_exception,         ONLY: finish
!-------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: check_ext_constants

  CHARACTER(LEN = *), PARAMETER :: modname = "mo_check_ext_constants"

CONTAINS

!-------------------------------------------------------------------------

  SUBROUTINE check_ext_constants()

    CHARACTER(len=*), PARAMETER :: routine =  modname//'::check_ext_constants'

    IF (modelname_ext /= modelname) CALL finish(routine, 'MODELNAME has deviated from its original value')

    IF (modelversion_ext /= modelversion) CALL finish(routine, 'MODELVERSION has deviated from its original value')

    IF (MAX_CHAR_LENGTH_ext /= MAX_CHAR_LENGTH) CALL finish(routine, 'MAX_CHAR_LENGTH has deviated from its &
       &                                                              original value')

    IF (SUCCESS_ext /= SUCCESS) CALL finish(routine, 'SUCCESS has deviated from its original value')

    IF (CELLS_ext /= CELLS) CALL finish(routine, 'CELLS has deviated from its original value')

    IF (EDGES_ext /= EDGES) CALL finish(routine, 'EDGES has deviated from its original value')

    IF (VERTS_ext /= VERTS) CALL finish(routine, 'VERTS has deviated from its original value')

    IF (ON_CELLS_ext /= ON_CELLS) CALL finish(routine, 'ON_CELLS has deviated from its original value')

    IF (ON_EDGES_ext /= ON_EDGES) CALL finish(routine, 'ON_EDGES has deviated from its original value')

    IF (ON_VERTICES_ext /= ON_VERTICES) CALL finish(routine, 'ON_VERTICES has deviated from its original value')

    IF (HALO_LEVELS_CEILING_ext /= HALO_LEVELS_CEILING) CALL finish(routine, 'HALO_LEVELS_CEILING has deviated &
       &                                                                      from its original value')

    IF (grf_bdywidth_c_ext /= grf_bdywidth_c) CALL finish(routine, 'GRF_BDYWIDTH_C has deviated from its original &
       &                                                            value')

    IF (grf_bdywidth_e_ext /= grf_bdywidth_e) CALL finish(routine, 'GRF_BDYWIDTH_E has deviated from its original &
       &                                                            value')

    IF (max_hw_ext /= max_hw) CALL finish(routine, 'MAX_HW has deviated from its original value')

    IF (min_rlcell_int_ext /= min_rlcell_int) CALL finish(routine, 'MIN_RLCELL_INT has deviated from its  &
       &                                                            original value')

    IF (min_rlcell_ext /= min_rlcell) CALL finish(routine, 'MIN_RLCELL has deviated from its original value')

    IF (max_rlcell_ext /= max_rlcell) CALL finish(routine, 'MAX_RLCELL has deviated from its original value')

    IF (min_rlvert_int_ext /= min_rlvert_int) CALL finish(routine, 'MIN_RLVERT_INT has deviated from its  &
       &                                                            original value')

    IF (min_rlvert_ext /= min_rlvert) CALL finish(routine, 'MIN_RLVERT has deviated from its original value')

    IF (max_rlvert_ext /= max_rlvert) CALL finish(routine, 'MAX_RLVERT has deviated from its original value')

    IF (min_rledge_int_ext /= min_rledge_int) CALL finish(routine, 'MIN_RLEDGE_INT has deviated from its  &
       &                                                            original value')

    IF (min_rledge_ext /= min_rledge) CALL finish(routine, 'MIN_RLEDGE has deviated from its original value')

    IF (max_rledge_ext /= max_rledge) CALL finish(routine, 'MAX_RLEDGE has deviated from its original value')

    IF (start_prog_cells_ext /= start_prog_cells) CALL finish(routine, 'START_PROG_CELLS has deviated from &
       &                                                                its original value')

    IF (end_prog_cells_ext /= end_prog_cells) CALL finish(routine, 'END_PROG_CELLS has deviated from &
       &                                                            its original value')

    IF (end_halo_lev1_cells_ext /= end_halo_lev1_cells) CALL finish(routine, 'END_HALO_LEV1_CELLS has deviated &
       &                                                                      from its original value')

    IF (end_all_cells_ext /= end_all_cells) CALL finish(routine, 'END_ALL_CELLS has deviated from its & 
       &                                                          original value')

    IF (start_prog_edges_ext /= start_prog_edges) CALL finish(routine, 'START_PROG_EDGES has deviated from &
       &                                                                its original value')

    IF (end_prog_edges_ext /= end_prog_edges) CALL finish(routine, 'END_PROG_EDGES has deviated from &
       &                                                            its original value')

    IF (end_halo_lev1_edges_ext /= end_halo_lev1_edges) CALL finish(routine, 'END_HALO_LEV1_EDGES has deviated &
       &                                                                      from its original value')

    IF (end_all_edges_ext /= end_all_edges) CALL finish(routine, 'END_ALL_EDGES has deviated from its & 
       &                                                          original value')

    IF (end_prog_verts_ext /= end_prog_verts) CALL finish(routine, 'END_PROG_VERTS has deviated from &
       &                                                            its original value')

    IF (end_verts_of_prog_cells_ext /= end_verts_of_prog_cells) CALL finish(routine, 'END_VERTS_OF_PROG_CELLS &
       &                                                                    has deviated from its original value')

    IF (end_all_verts_ext /= end_all_verts) CALL finish(routine, 'END_ALL_VERTS has deviated from its &
       &                                                          original value')

    IF (max_dom_ext /= max_dom) CALL finish(routine, 'MAX_DOM has deviated from its original value')

    IF (max_dom_dig10_ext /= max_dom_dig10) CALL finish(routine, 'MAX_DOM_dig10 has deviated from & 
       &                                                          its original value')

    IF (max_phys_dom_ext /= max_phys_dom) CALL finish(routine, 'MAX_PHYS_DOM has deviated from its original value')

    IF (max_ntracer_ext /= max_ntracer) CALL finish(routine, 'MAX_NTRACER has deviated from its original value')

    IF (max_echotop_ext /= max_echotop) CALL finish(routine, 'MAX_ECHOTOP has deviated from its original value')

    IF (max_wshear_ext /= max_wshear) CALL finish(routine, 'MAX_WSHEAR has deviated from its original value')

    IF (max_srh_ext /= max_srh) CALL finish(routine, 'MAX_SRH has deviated from its original value')

    IF (TORUS_MAX_LAT_ext /= TORUS_MAX_LAT) CALL finish(routine, 'TORUS_MAX_LAT has deviated from its &
       &                                                         original value')

  END SUBROUTINE check_ext_constants

END MODULE mo_check_ext_constants

