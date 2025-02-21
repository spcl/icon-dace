! Provide an implementation of the ocean tracer nuding functionality.
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

!----------------------------
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

! DSL definitions 


























!---------------------
! block definitions




!---------------------
! mappings



















!---------------------
! connectivity









!---------------------
! generic types




!---------------------
! shortcuts







































!---------------------
! Upper-lower case
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!



!----------------------------
MODULE mo_ocean_nudging_types
!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  !USE mo_math_utilities,            ONLY: t_cartesian_coordinates
  !USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  !USE mo_math_constants,            ONLY: pi
  !USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,              &
  !  & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
  !  & type_3dimrelax_temp, para_3dimrelax_temp,                           &
  !  & type_3dimrelax_salt, para_3dimrelax_salt
  !USE mo_util_dbg_prnt,             ONLY: dbg_print
  !USE mo_parallel_config,           ONLY: nproma
  !USE mo_dynamics_config,           ONLY: nold, nnew
  !USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
!!  USE mo_ocean_types,               ONLY: t_hydro_ocean_state, t_ocean_tracer
!  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
!  USE mo_exception,                 ONLY: finish !, message_text, message
!  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
!  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range


  IMPLICIT NONE

  TYPE t_ocean_nudge

    ! Variables for 3-dim tracer relaxation:
    REAL(wp),  POINTER, DIMENSION(:,:,:) :: &
      & data_3dimRelax_Temp, & ! 3-dim temperature relaxation data (T*)
      & forc_3dimRelax_Temp, & ! 3-dim temperature relaxation forcing (1/tau*(T-T*))
      & data_3dimRelax_Salt, & ! 3-dim salinity relaxation data (T*)
      & forc_3dimRelax_Salt, &    ! 3-dim salinity relaxation forcing (1/tau*(T-T*))
      & relax_3dim_coefficient ! 3-dim relaxation coefficient when the relaxation varies

  END TYPE t_ocean_nudge

  PRIVATE

  PUBLIC :: t_ocean_nudge

END MODULE mo_ocean_nudging_types


