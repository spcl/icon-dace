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

! Configuration setup for wave energy transport

MODULE mo_energy_propagation_config

  USE mo_kind,                      ONLY: wp
  USE mo_impl_constants,            ONLY: max_dom, VNAME_LEN, MAX_NTRACER

  IMPLICIT NONE

  PRIVATE

  ! types
  PUBLIC :: t_energy_propagation_config

  ! variables
  PUBLIC :: energy_propagation_config


  CHARACTER(LEN = *), PARAMETER :: modname = "mo_energy_propagation_config"


  !!--------------------------------------------------------------------------
  !! Basic configuration setup for energy propagation
  !!--------------------------------------------------------------------------
  TYPE :: t_energy_propagation_config

    INTEGER :: &           !< parameter used to select the limiter
      &  itype_limit       !< for horizontal transport
                           !< 0: no limiter
                           !< 3: monotonous flux limiter
                           !< 4: positive definite flux limiter

    INTEGER :: &           !< parameter used to select the gradient
      &  igrad_c_miura     !< reconstruction method at cell center
                           !< for second order miura scheme


    REAL(wp):: &           !< global boost factor for range of permissible values in
      &  beta_fct          !< (semi-) monotonous flux limiter. A value larger than
                           !< 1 allows for (small) over and undershoots, while a value
                           !< of 1 gives strict monotonicity (at the price of increased
                           !< diffusivity).

    LOGICAL :: &           !< if .TRUE., calculate grid refraction
      &  lgrid_refr

    CHARACTER(len=VNAME_LEN) ::  &  !< tracer-specific name suffixes
      &  tracer_names(MAX_NTRACER)  !< set by namelist, e.g. 'hus' for specific humidity
                                    !< default: 'q<tracer index>'.

  END TYPE t_energy_propagation_config

  !>
  !!
  TYPE(t_energy_propagation_config), TARGET :: energy_propagation_config(1:max_dom)


CONTAINS

END MODULE mo_energy_propagation_config
