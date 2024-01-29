!
! Subroutine aes_phy_main calls all the parameterization schemes
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

#if defined __xlC__ && !defined NOXLFPROCESS
@PROCESS HOT
@PROCESS SPILLSIZE(5000)
#endif
!OCL NOALIAS

MODULE mo_aes_phy_main

  USE mo_kind                ,ONLY: wp
  USE mo_exception           ,ONLY: message

  USE mtime                  ,ONLY: t_datetime => datetime, isCurrentEventActive, &
       &                            OPERATOR(<=), OPERATOR(>)

  USE mo_model_domain        ,ONLY: t_patch

  USE mo_omp_block_loop      ,ONLY: omp_block_loop_cell

  USE mo_aes_phy_config      ,ONLY: aes_phy_tc, dt_zero
  USE mo_aes_vdf_config      ,ONLY: aes_vdf_config
  USE mo_aes_phy_diag        ,ONLY: surface_fractions, &
    &                               droplet_number,    &
    &                               cpair_cvair_qconv, &
    &                               initialize,        &
    &                               finalize

  USE mo_aes_diagnostics     ,ONLY: aes_global_diagnostics
#if defined( _OPENACC )
  USE mo_exception           ,ONLY: warning
  USE mo_var_list_gpu        ,ONLY: gpu_update_var_list
#endif

  USE mo_diagnose_cov        ,ONLY: diagnose_cov

  USE mo_interface_aes_wmo   ,ONLY: interface_aes_wmo
  USE mo_interface_aes_rad   ,ONLY: interface_aes_rad
  USE mo_interface_aes_rht   ,ONLY: interface_aes_rht
  USE mo_interface_aes_vdf   ,ONLY: interface_aes_vdf
  USE mo_interface_aes_tmx   ,ONLY: interface_aes_tmx
  USE mo_interface_aes_car   ,ONLY: interface_aes_car
  USE mo_interface_aes_art   ,ONLY: interface_aes_art
  !
  USE mo_interface_cloud_mig ,ONLY: interface_cloud_mig
  !
  USE mo_interface_cloud_two ,ONLY: interface_cloud_two

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: aes_phy_main

CONTAINS

  !>
  !!
  SUBROUTINE aes_phy_main(patch, datetime, pdtime)


    ! Arguments
    !
    TYPE(t_patch)  ,TARGET ,INTENT(INOUT) :: patch
    TYPE(t_datetime)       ,POINTER       :: datetime
    REAL(wp)               ,INTENT(IN)    :: pdtime

    ! Local variables
    !
    INTEGER  :: jg                                         !< grid level/domain index

    jg = patch%id

    ! store grid specific time parameters for physics
    !
    aes_phy_tc(jg)%dt_phy_sec =  pdtime
    aes_phy_tc(jg)%datetime   => datetime

    !-------------------------------------------------------------------
    ! Prepare for physics
    !-------------------------------------------------------------------
    !
    CALL omp_block_loop_cell(patch, initialize)        ! initialize q_phy and q_phy_vi
    CALL omp_block_loop_cell(patch, cpair_cvair_qconv) ! cp, cv and W/m2 -> K/s
    CALL omp_block_loop_cell(patch, surface_fractions) ! surface fractions

    !-------------------------------------------------------------------
    ! single moment cloud microphysics "Graupel" (mig)
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_mig > dt_zero ) THEN
       !
       aes_phy_tc(jg)%is_in_sd_ed_interval_mig = (aes_phy_tc(jg)%sd_mig <= datetime) .AND. (aes_phy_tc(jg)%ed_mig > datetime)
       aes_phy_tc(jg)%is_active_mig            = isCurrentEventActive(aes_phy_tc(jg)%ev_mig, datetime)
       !
       CALL message_forcing_action('graupel microphysics (mig)',            &
            &                      aes_phy_tc(jg)%is_in_sd_ed_interval_mig, &
            &                      aes_phy_tc(jg)%is_active_mig)
       !
       CALL omp_block_loop_cell(patch, interface_cloud_mig)
       !
    END IF

    !--------------------------------------------------------------------
    ! two-moment cloud microphysics (two) by Seifert and Beheng (2006)
    !--------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_two > dt_zero ) THEN
#if defined( _OPENACC )
       CALL warning('GPU:aes_art_main','GPU host synchronization should be removed when port is done!')
       CALL gpu_update_var_list('prm_field_D', .false., jg, lacc=.TRUE.)
       CALL gpu_update_var_list('prm_tend_D' , .false., jg, lacc=.TRUE.)
#endif
       !
       aes_phy_tc(jg)%is_in_sd_ed_interval_two = (aes_phy_tc(jg)%sd_two <= datetime) .AND. (aes_phy_tc(jg)%ed_two > datetime)
       aes_phy_tc(jg)%is_active_two            = isCurrentEventActive(aes_phy_tc(jg)%ev_two, datetime)
       !
       CALL message_forcing_action('two-moment bulk microphysics (two)',    &
            &                      aes_phy_tc(jg)%is_in_sd_ed_interval_two, &
            &                      aes_phy_tc(jg)%is_active_two)
       !
       CALL omp_block_loop_cell(patch, interface_cloud_two)
       !
#if defined( _OPENACC )
       CALL warning('GPU:aes_art_main','GPU device synchronization should be removed when port is done!')
       CALL gpu_update_var_list('prm_field_D', .true., jg, lacc=.TRUE.)
       CALL gpu_update_var_list('prm_tend_D' , .true., jg, lacc=.TRUE.)
#endif
    END IF

    !-------------------------------------------------------------------
    ! Radiation (LW+SW)
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_rad > dt_zero ) THEN
       !
       aes_phy_tc(jg)%is_in_sd_ed_interval_rad = (aes_phy_tc(jg)%sd_rad <= datetime) .AND. (aes_phy_tc(jg)%ed_rad > datetime)
       aes_phy_tc(jg)%is_active_rad            = isCurrentEventActive(aes_phy_tc(jg)%ev_rad, datetime)
       !
       CALL message_forcing_action('LW and SW radiation (rad:fluxes )',     &
            &                      aes_phy_tc(jg)%is_in_sd_ed_interval_rad, &
            &                      aes_phy_tc(jg)%is_active_rad)
       !
       ! cloud droplet number concentration
       CALL omp_block_loop_cell(patch, droplet_number)
       !
       ! radiative fluxes
       CALL omp_block_loop_cell(patch, interface_aes_rad)
       !
       ! radiative heating is always active
       aes_phy_tc(jg)%is_active_rad = .TRUE.
       !
       CALL message_forcing_action('LW and SW radiation (rht:heating)',     &
            &                      aes_phy_tc(jg)%is_in_sd_ed_interval_rad, &
            &                      aes_phy_tc(jg)%is_active_rad)
       !
       ! radiative heating
       CALL omp_block_loop_cell(patch, interface_aes_rht)
       !
    END IF

    !-------------------------------------------------------------------
    ! Vertical diffusion, boundary layer and surface
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_vdf > dt_zero ) THEN
       !
       aes_phy_tc(jg)%is_in_sd_ed_interval_vdf = (aes_phy_tc(jg)%sd_vdf <= datetime) .AND. (aes_phy_tc(jg)%ed_vdf > datetime)
       aes_phy_tc(jg)%is_active_vdf            = isCurrentEventActive(aes_phy_tc(jg)%ev_vdf, datetime)
       !
       IF (aes_vdf_config(jg)%use_tmx) THEN
         CALL message_forcing_action('vertical diffusion (tmx)',              &
              &                      aes_phy_tc(jg)%is_in_sd_ed_interval_vdf, &
              &                      aes_phy_tc(jg)%is_active_vdf)

         CALL interface_aes_tmx(patch,                                        &
            &                   aes_phy_tc(jg)%is_in_sd_ed_interval_vdf,      &
            &                   aes_phy_tc(jg)%is_active_vdf,                 &
            &                   datetime, pdtime)
       ELSE
         CALL message_forcing_action('vertical diffusion (vdf)',              &
              &                      aes_phy_tc(jg)%is_in_sd_ed_interval_vdf, &
              &                      aes_phy_tc(jg)%is_active_vdf)
         !
         CALL omp_block_loop_cell(patch, diagnose_cov) ! cloud cover before turbulent diffusion (only for TKE scheme)
         !
         CALL interface_aes_vdf(patch)
       END IF
       !
    END IF

    !-------------------------------------------------------------------
    ! Linearized ozone chemistry of Cariolle
    !-------------------------------------------------------------------
    !
    IF ( aes_phy_tc(jg)%dt_car > dt_zero ) THEN
       !
       aes_phy_tc(jg)%is_in_sd_ed_interval_car = (aes_phy_tc(jg)%sd_car <= datetime) .AND. (aes_phy_tc(jg)%ed_car > datetime)
       aes_phy_tc(jg)%is_active_car            = isCurrentEventActive(aes_phy_tc(jg)%ev_car, datetime)
       !
       CALL message_forcing_action('lin. Cariolle ozone chem. (car)',       &
            &                      aes_phy_tc(jg)%is_in_sd_ed_interval_car, &
            &                      aes_phy_tc(jg)%is_active_car)
       !
       CALL omp_block_loop_cell(patch, interface_aes_car)
       !
    END IF

    !-------------------------------------------------------------------
    ! Atmospheric chemistry of ART
    !-------------------------------------------------------------------
    !
    IF (aes_phy_tc(jg)%dt_art > dt_zero) THEN
#if defined( _OPENACC )
       CALL warning('GPU:aes_art_main','GPU host synchronization should be removed when port is done!')
       CALL gpu_update_var_list('prm_field_D', .false., jg, lacc=.TRUE.)
       CALL gpu_update_var_list('prm_tend_D' , .false., jg, lacc=.TRUE.)
#endif
       !
       aes_phy_tc(jg)%is_in_sd_ed_interval_art = (aes_phy_tc(jg)%sd_art <= datetime) .AND. (aes_phy_tc(jg)%ed_art > datetime)
       aes_phy_tc(jg)%is_active_art            = isCurrentEventActive(aes_phy_tc(jg)%ev_art, datetime)
       !
       CALL message_forcing_action('ART (art)',                             &
            &                      aes_phy_tc(jg)%is_in_sd_ed_interval_art, &
            &                      aes_phy_tc(jg)%is_active_art)
       !
       ! OMP loops are hidden inside the ART routines. Hence the full patch needs
       ! to be passed to the ART routines and is it not possible to call the
       ! ART reaction interface inside the standard omp block loop.
       ! This should be reprogrammed.
       !
       CALL interface_aes_art(patch)
       !
#if defined( _OPENACC )
       CALL warning('GPU:aes_art_main','GPU device synchronization should be removed when port is done!')
       CALL gpu_update_var_list('prm_field_D', .true., jg, lacc=.TRUE.)
       CALL gpu_update_var_list('prm_tend_D' , .true., jg, lacc=.TRUE.)
#endif
    END IF

    !-------------------------------------------------------------------
    ! Finish physics
    !-------------------------------------------------------------------
    !
    CALL omp_block_loop_cell(patch,finalize)           ! dT/dt|phy,const.pressure -> dT/dt|phy,const.volume

    !-------------------------------------------------------------------
    ! Output diagnostics
    !-------------------------------------------------------------------
    !
    CALL omp_block_loop_cell(patch, diagnose_cov)      ! cloud cover
    CALL omp_block_loop_cell(patch, interface_aes_wmo) ! WMO tropopause height
    CALL aes_global_diagnostics(patch)                 ! global mean diagnostics

  END SUBROUTINE aes_phy_main
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  SUBROUTINE message_forcing_action(process, is_in_sd_ed_interval, is_active)
    CHARACTER(LEN=*) ,INTENT(in) :: process
    LOGICAL          ,INTENT(in) :: is_in_sd_ed_interval
    LOGICAL          ,INTENT(in) :: is_active

    IF (is_in_sd_ed_interval) THEN
       IF (is_active) THEN
          CALL message('aes_phy_main','compute forcing by '//process)
       ELSE
          CALL message('aes_phy_main','recycle forcing by '//process)
       END IF
    ELSE
       CALL    message('aes_phy_main','no      forcing by '//process)
    END IF

  END SUBROUTINE message_forcing_action
  !---------------------------------------------------------------------


END MODULE mo_aes_phy_main
