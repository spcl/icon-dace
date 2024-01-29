! Contains the implementation of the semi-implicit Adams-Bashforth timestepping
! for the ICON ocean model.
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

MODULE mo_ocean_ab_timestepping
  USE mo_ocean_nml,                      ONLY: discretization_scheme
  USE mo_dynamics_config,                ONLY: nold, nnew
  USE mo_ocean_surface_types,            ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_model_domain,                   ONLY: t_patch_3D !, t_patch
  USE mo_ext_data_types,                 ONLY: t_external_data
  USE mo_ocean_ab_timestepping_mimetic,  ONLY: solve_free_sfc_ab_mimetic, &
    &  calc_normal_velocity_ab_mimetic, &
    &  calc_vert_velocity_mim_bottomup
  USE mo_ocean_physics_types,            ONLY: t_ho_params
  USE mo_ocean_types,                    ONLY: t_hydro_ocean_state, t_operator_coeff, t_solverCoeff_singlePrecision
  USE mo_exception,                      ONLY: finish!, message_text
  USE mo_fortran_tools,                  ONLY: set_acc_host_or_device

IMPLICIT NONE

PRIVATE

  INTEGER, PARAMETER :: MIMETIC_TYPE = 1
  INTEGER, PARAMETER :: RBF_TYPE     = 2

  PUBLIC :: solve_free_surface_eq_ab
  PUBLIC :: calc_normal_velocity_ab
  PUBLIC :: calc_vert_velocity
  PUBLIC :: update_time_indices

CONTAINS
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! !  Solves the free surface equation.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE solve_free_surface_eq_ab(patch_3D, ocean_state, external_data, p_as, p_oce_sfc, &
    & physics_parameters, timestep, op_coeffs, solverCoeff_sp, return_status, lacc)
    TYPE(t_patch_3D ), INTENT(IN), POINTER :: patch_3D
    TYPE(t_hydro_ocean_state) :: ocean_state
    TYPE(t_external_data) :: external_data
    TYPE(t_ocean_surface), INTENT(INOUT) :: p_oce_sfc
    TYPE(t_atmos_for_ocean), INTENT(INOUT) :: p_as
    TYPE (t_ho_params) :: physics_parameters
    INTEGER, INTENT(IN) :: timestep
    TYPE(t_operator_coeff), INTENT(IN), TARGET :: op_coeffs
    TYPE(t_solverCoeff_singlePrecision), INTENT(in), TARGET :: solverCoeff_sp
    INTEGER :: return_status
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL solve_free_sfc_ab_mimetic( patch_3D, ocean_state, external_data, p_as, p_oce_sfc, &
        & physics_parameters, timestep, op_coeffs, solverCoeff_sp, return_status, lacc=lzacc)

    ELSE
      CALL finish ('solve_free_surface_eq_ab: ',' Discretization type not supported !!')
    ENDIF

  END SUBROUTINE solve_free_surface_eq_ab
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computation of new velocity in Adams-Bashforth timestepping.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_normal_velocity_ab(patch_3D, ocean_state, operators_coefficients, &
    & solverCoeff_sp, external_data, physics_parameters, lacc)
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state), TARGET    :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(in) :: solverCoeff_sp
    TYPE(t_external_data), TARGET        :: external_data
    TYPE (t_ho_params)                   :: physics_parameters

    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    IF(discretization_scheme==MIMETIC_TYPE)THEN
      CALL calc_normal_velocity_ab_mimetic(patch_3D, ocean_state, operators_coefficients, lacc=lzacc)
    ELSE
      CALL finish ('calc_normal_velocity_ab: ',' Discreization type not supported !!')
    ENDIF

  END SUBROUTINE calc_normal_velocity_ab
  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Computation of new vertical velocity using continuity equation

  !! Calculate diagnostic vertical velocity from horizontal velocity using the
  !! incommpressibility condition in the continuity equation.
  !! For the case of the semi-implicit-AB scheme the land-sea-mask may be applied
  !! at least after collecting the whole explicit term.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE calc_vert_velocity(patch_3D, ocean_state, operators_coefficients, lacc)
    TYPE(t_patch_3D ),TARGET, INTENT(IN) :: patch_3D
    TYPE(t_hydro_ocean_state)            :: ocean_state
    TYPE(t_operator_coeff), INTENT(IN) :: operators_coefficients
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: lzacc
    !-----------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF(discretization_scheme==MIMETIC_TYPE)THEN

      CALL calc_vert_velocity_mim_bottomup( patch_3D,       &
                                  & ocean_state,            &
                                  & operators_coefficients, &
                                  & lacc=lzacc )

    ELSE
      CALL finish ('calc_vert_velocity: ',' Discretization type not supported !!')
    END IF

  END SUBROUTINE calc_vert_velocity
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE update_time_indices(jg)
    INTEGER, INTENT(IN) :: jg
    INTEGER             :: n_temp
    ! Step 7: Swap time indices before output
    !         half time levels of semi-implicit Adams-Bashforth timestepping are
    !         stored in auxiliary arrays g_n and g_nimd of p_diag%aux
    n_temp    = nold(jg)
    nold(jg)  = nnew(jg)
    nnew(jg)  = n_temp
  END SUBROUTINE update_time_indices
END MODULE mo_ocean_ab_timestepping
