! Contains the variables to set up the ocean model.
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
#include "iconfor_dsl_definitions.inc"
!=============================================================================================
MODULE mo_ocean_tracer_transport_types

  USE mo_kind,                ONLY: wp, sp
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_model_domain,        ONLY: t_patch_3d
  USE mo_impl_constants,      ONLY: max_char_length

  PUBLIC :: t_ocean_tracer, t_tracer_collection, t_ocean_transport_state
  
  !----------------------------------------------
  TYPE t_ocean_tracer
    onCells :: concentration

    onEdges :: hor_diffusion_coeff
    onCells_HalfLevels :: ver_diffusion_coeff

    onCells_2D :: top_bc, bottom_bc
    
    onCells_3D :: vertical_trasnport_tendencies

    LOGICAL :: is_advected

    TYPE(t_tracer_metadata), POINTER :: metadata

  END TYPE t_ocean_tracer
  !----------------------------------------------
    
  !-------------------------------
  TYPE t_tracer_metadata
    CHARACTER(LEN=max_char_length) :: tracer_longnames
    CHARACTER(LEN=max_char_length) :: tracer_stdnames
    CHARACTER(LEN=max_char_length) :: tracer_shortnames
    CHARACTER(LEN=max_char_length) :: tracer_units
    INTEGER                        :: tracer_codes
  END TYPE t_tracer_metadata
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  TYPE t_tracer_collection
    TYPE(t_patch_3d ),POINTER :: patch_3d
    CHARACTER(LEN=max_char_length) :: typeOfTracers
    INTEGER :: no_of_tracers
    TYPE(t_ocean_tracer), POINTER :: tracer(:)
  END TYPE t_tracer_collection
  !-------------------------------------------------------------------------------

  !----------------------------------------------
  TYPE t_ocean_transport_state
    TYPE(t_patch_3d ),POINTER :: patch_3d
    
    onCells_2D :: h_new, h_old

    onEdges    :: mass_flux_e
    onEdges    :: vn

    onCells_HalfLevels :: w

  END TYPE t_ocean_transport_state
  !----------------------------------------------
    
    
END MODULE mo_ocean_tracer_transport_types

