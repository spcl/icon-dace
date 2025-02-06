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

! @brief Initialisation of atmosphere coupling

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_atmo_coupling_frame

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config
  USE mo_run_config          ,ONLY: iforcing, ltimer
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_coupling_init
  USE mo_impl_constants      ,ONLY: MAX_CHAR_LENGTH, inwp, LSS_JSBACH
  USE mo_ext_data_types      ,ONLY: t_external_data

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__) && defined(YAC_coupling)
  USE mo_interface_hd_ocean  ,ONLY: jsb_fdef_hd_fields
#endif

  USE mo_coupling_config     ,ONLY: is_coupled_run, is_coupled_to_ocean, &
    &                               is_coupled_to_hydrodisc, &
    &                               is_coupled_to_waves,     &
    &                               is_coupled_to_output,    &
    &                               is_coupled_to_aero,      &
    &                               is_coupled_to_o3
  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_time_config         ,ONLY: time_config

  USE mo_atmo_wave_coupling  ,ONLY: construct_atmo_wave_coupling
  USE mo_atmo_ocean_coupling ,ONLY: construct_atmo_ocean_coupling
  USE mo_atmo_o3_provider_coupling,ONLY: &
    construct_atmo_o3_provider_coupling_post_sync
  USE mo_atmo_aero_provider_coupling,ONLY: &
    construct_atmo_aero_provider_coupling_post_sync
  USE mo_nwp_hydrodisc_coupling,ONLY: construct_nwp_hydrodisc_coupling

  USE mo_exception           ,ONLY: finish, message

  USE mo_coupling_utils      ,ONLY: cpl_def_main, cpl_sync_def, cpl_enddef

  USE mtime                  ,ONLY: timedeltaToString, MAX_TIMEDELTA_STR_LEN

#ifndef __NO_ICON_COMIN__
  USE comin_host_interface, ONLY: EP_ATM_YAC_DEFCOMP_BEFORE,       &
       &                          EP_ATM_YAC_DEFCOMP_AFTER,        &
       &                          EP_ATM_YAC_SYNCDEF_BEFORE,       &
       &                          EP_ATM_YAC_SYNCDEF_AFTER,        &
       &                          EP_ATM_YAC_ENDDEF_BEFORE,        &
       &                          EP_ATM_YAC_ENDDEF_AFTER,         &
       &                          COMIN_DOMAIN_OUTSIDE_LOOP
  USE mo_comin_adapter,     ONLY: icon_call_callback
#endif

  USE mo_output_coupling     ,ONLY: construct_output_coupling, &
       &                            construct_output_coupling_finalize

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_atmo_coupling_frame' ! Output of module for debug

  PUBLIC :: construct_atmo_coupling, destruct_atmo_coupling
  PUBLIC :: nbr_inner_cells

  INTEGER, SAVE :: nbr_inner_cells

CONTAINS

  !>
  !! SUBROUTINE construct_atmo_coupling -- the initialisation for the coupling
  !! of the atmosphere
  SUBROUTINE construct_atmo_coupling (p_patch, ext_data)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    TYPE(t_external_data), INTENT(IN) :: ext_data(:)

    TYPE(t_patch), POINTER :: patch_horz

    !---------------------------------------------------------------------
    ! 11. Do the setup for the coupled run
    !
    ! For the time being this could all go into a subroutine which is
    ! common to atmo and ocean. Does this make sense if the setup deviates
    ! too much in future.
    !---------------------------------------------------------------------

    INTEGER :: comp_id, output_comp_id
    INTEGER :: grid_id
    INTEGER :: cell_point_id, vertex_point_id

    INTEGER :: jg

    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: timestepstring

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':construct_atmo_coupling'

    IF ( .NOT. is_coupled_run() ) RETURN

    IF (ltimer) CALL timer_start (timer_coupling_init)

    CALL message(str_module, 'Constructing the atmosphere coupling frame.')

    jg = 1
    patch_horz => p_patch(jg)

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_DEFCOMP_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, lacc=.FALSE.)
#endif

    ! Do basic initialisation of the component
    IF( is_coupled_to_output() ) THEN
      CALL cpl_def_main(routine,           & !in
                        patch_horz,        & !in
                        "icon_atmos_grid", & !in
                        comp_id,           & !out
                        output_comp_id,    & !out
                        grid_id,           & !out
                        cell_point_id,     & !out
                        vertex_point_id,   & !out
                        nbr_inner_cells)     !out
    ELSE
      CALL cpl_def_main(routine,           & !in
                        patch_horz,        & !in
                        "icon_atmos_grid", & !in
                        comp_id,           & !out
                        grid_id,           & !out
                        cell_point_id,     & !out
                        nbr_inner_cells)     !out
    ENDIF

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_DEFCOMP_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP, lacc=.FALSE.)
#endif

    ! get model timestep
    CALL timedeltaToString(time_config%tc_dt_model, timestepstring)

    IF( is_coupled_to_output() ) THEN

      CALL message(str_module, 'Constructing the coupling frame atmosphere-output.')

      CALL construct_output_coupling ( &
        p_patch, output_comp_id, cell_point_id, vertex_point_id, timestepstring)

    END IF

    IF ( is_coupled_to_ocean() ) THEN

      CALL message(str_module, 'Constructing the coupling frame atmosphere-ocean.')

      CALL construct_atmo_ocean_coupling( &
        p_patch, ext_data, comp_id, grid_id, cell_point_id, timestepstring)

#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__) && defined(YAC_coupling)

      ! Define coupling of runoff if HD model is present and interface is coded
      !  - discrimination between Proto2 (no HD) and Proto3 (with HD) is needed
      ! preliminary: coupling to jsbach/hd is active
      IF ( (iforcing /= INWP .AND. aes_phy_config(jg)%ljsb) .OR. ( atm_phy_nwp_config(jg)%inwp_surface == LSS_JSBACH .AND. .NOT. is_coupled_to_hydrodisc() ) ) THEN

        ! Construct coupling frame for atmosphere/JSBACH-hydrological discharge
        CALL message(str_module, 'Constructing the coupling frame atmosphere/JSBACH-hydrological discharge.')

        CALL jsb_fdef_hd_fields(comp_id, (/cell_point_id/), grid_id, patch_horz%n_patch_cells)

      ENDIF

#endif
    ENDIF   ! Construct coupling frame for atmosphere-ocean

    IF ( is_coupled_to_hydrodisc() ) THEN

      ! Construct coupling frame for atmosphere-hydrological discharge
      CALL message(str_module, 'Constructing the coupling frame atmosphere-hydrological discharge.')


      CALL construct_nwp_hydrodisc_coupling( &
        p_patch, ext_data, comp_id, grid_id, cell_point_id, timestepstring)

    ENDIF ! Construct coupling frame for atmosphere-hydrological discharge

    IF ( is_coupled_to_waves() ) THEN

      CALL message(str_module, 'Constructing the coupling frame atmosphere-wave.')

      CALL construct_atmo_wave_coupling( &
        comp_id, cell_point_id, timestepstring)

    END IF

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_SYNCDEF_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, lacc=.FALSE.)
#endif

    ! Synchronize all definitions until this point with other components
    CALL cpl_sync_def(str_module)

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_SYNCDEF_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP, lacc=.FALSE.)
#endif

    ! add Ozone data field if needed
    IF  ( is_coupled_to_o3() ) THEN

      ! Construct coupling frame for atmosphere-o3 provider
      CALL message(str_module, 'Constructing the coupling frame atmosphere-o3 provider.')

      CALL construct_atmo_o3_provider_coupling_post_sync( &
        comp_id, cell_point_id, TRIM(aes_phy_config(jg)%dt_rad))

    END IF

    IF ( is_coupled_to_aero() ) THEN

      ! Construct coupling frame for atmosphere-aero provider
      CALL message(str_module, 'Constructing the coupling frame atmosphere-aero provider.')

      CALL construct_atmo_aero_provider_coupling_post_sync( &
        comp_id, cell_point_id, TRIM(aes_phy_config(jg)%dt_rad))

    END IF

    ! End definition of coupling fields and search

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_ENDDEF_BEFORE, COMIN_DOMAIN_OUTSIDE_LOOP, lacc=.FALSE.)
#endif

    CALL cpl_enddef(str_module)

#ifndef __NO_ICON_COMIN__
    CALL icon_call_callback(EP_ATM_YAC_ENDDEF_AFTER, COMIN_DOMAIN_OUTSIDE_LOOP, lacc=.FALSE.)
#endif

    ! finalizes the output coupling
    IF( is_coupled_to_output() ) CALL construct_output_coupling_finalize()

    IF (ltimer) CALL timer_stop(timer_coupling_init)

  END SUBROUTINE construct_atmo_coupling

  !>
  !! SUBROUTINE destruct_atmo_coupling -- the finalization for the coupling
  !! of the atmosphere
  SUBROUTINE destruct_atmo_coupling()

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':destruct_atmo_coupling'

    IF ( is_coupled_run() ) THEN

      CALL message(str_module, 'Destructing the atmosphere coupling frame.')

    END IF

  END SUBROUTINE destruct_atmo_coupling

END MODULE mo_atmo_coupling_frame

!vim:fdm=marker
