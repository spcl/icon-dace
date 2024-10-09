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

MODULE mo_sleve_config

  USE mo_kind,                ONLY: wp
 !USE mo_impl_constants,      ONLY: max_dom

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: itype_laydistr, min_lay_thckn, max_lay_thckn, htop_thcknlimit, nshift_above_thcklay, stretch_fac
  PUBLIC :: decay_scale_1, decay_scale_2, decay_exp, flat_height, top_height
  PUBLIC :: lread_smt
  !>
  !!--------------------------------------------------------------------------
  !! Type definition 
  !!--------------------------------------------------------------------------
  !TYPE :: t_sleve_config

    ! a) Parameters specifying the distrubution of the coordinate surfaces

    INTEGER :: itype_laydistr ! Type of analytical function used for computing the coordinate surface distribution
    REAL(wp):: min_lay_thckn  ! Layer thickness of lowermost level
    REAL(wp):: max_lay_thckn  ! Maximum layer thickness below htop_thcknlimit
    REAL(wp):: htop_thcknlimit! Height below which the layer thickness must not exceed max_lay_thckn
    REAL(wp):: stretch_fac    ! Factor for stretching/squeezing the model layer distribution
    REAL(wp):: top_height     ! Height of model top

    INTEGER :: nshift_above_thcklay ! Shift above constant-thickness layer for further calculation of layer distribution

    ! b) Parameters for SLEVE definition

    REAL(wp):: decay_scale_1  ! Decay scale for large-scale topography component
    REAL(wp):: decay_scale_2  ! Decay scale for small-scale topography component
    REAL(wp):: decay_exp      ! Exponent for decay function
    REAL(wp):: flat_height    ! Height above which the coordinate surfaces are exactly flat
                              ! additional feature not available in the standard
                              ! SLEVE definition

    ! c) Parameter for reading in smoothed topography
    LOGICAL :: lread_smt

  !END TYPE t_sleve_config
  !>
  !!
  !TYPE(t_sleve_config) :: sleve_config(max_dom)

END MODULE mo_sleve_config
