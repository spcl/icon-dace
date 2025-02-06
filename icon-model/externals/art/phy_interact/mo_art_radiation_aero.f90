!
! mo_art_radiation_aero
! Replaces aerosol optical properties of Tegen climatology by values
! obtained from current aerosol concentrations from the ART modules.
! Based on: Philipp Gasch (2016): Numerical simulations of an exceptional dust event
!                                 in the Eastern Mediterranean including the mineral
!                                 dust radiative feedback.
!                                 Master Thesis at the Fakultaet fuer Physik of
!                                 the Karlsruhe Institute of Technology (KIT)
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

MODULE mo_art_radiation_aero
  USE mo_kind,                          ONLY: wp
  USE mo_dynamics_config,               ONLY: nnew_rcf, nnew
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_nonhydro_state,                ONLY: p_nh_state
  USE mo_parallel_config,               ONLY: nproma
  USE mo_exception,                     ONLY: message
  ! ART ROUTINES
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_2mom
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_data,                      ONLY: p_art_data


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_radiation_aero

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_radiation_aero(zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,    &
  &                           zaea_rrtm,zaes_rrtm,zaeg_rrtm,    &
  &                           jg,jb,ks,ke,jcs,jce,nlong,nshort, &
  &                           aer_tau_lw_vr,                    &
  &                           aer_tau_sw_vr,                    &
  &                           aer_piz_sw_vr,                    &
  &                           aer_cg_sw_vr)
!<
! SUBROUTINE art_radiation_aero
! This SR replaces aerosol optical properties of Tegen climatology by values
! obtained from current aerosol concentrations from the ART modules.
! Based on: Philipp Gasch (2016): Numerical simulations of an exceptional dust event
!                                 in the Eastern Mediterranean including the mineral
!                                 dust radiative feedback.
!                                 Master Thesis at the Fakultaet fuer Physik of 
!                                 the Karlsruhe Institute of Technology (KIT)
! Part of Module: mo_art_radiation_aero
! Author: Carolin Walter, KIT, Philipp Gasch, KIT, Daniel Rieger, KIT
! Initial Release: 2015-11-12
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  REAL(wp),INTENT(in) ::     & !< Note that these values are only valid for Tegen climatology!
    &  zaeq1(:,:),           & !< continental aerosol
    &  zaeq2(:,:),           & !< maritime aerosol
    &  zaeq3(:,:),           & !< mineral dust aerosol
    &  zaeq4(:,:),           & !< urban aerosol
    &  zaeq5(:,:),           & !< stratospheric background aerosol
    &  zaea_rrtm(:,:),       & !< absorption coefficient
    &  zaes_rrtm(:,:),       & !< scatter coefficient
    &  zaeg_rrtm(:,:)          !< asymmetry coefficient
  INTEGER,INTENT(in) ::      &
    &  jg,jb,                & !< Domain ID, block index
    &  ks, ke,               & !< loop indices jk loop
    &  jcs, jce,             & !< loop indices jc loop
    &  nlong,nshort            !< number of bands long/shortwave. Sorting of arrays: 
                               !      longwave first, then shortwave bands

  REAL(wp), INTENT(out) ::  & !< DIMS: (nproma,nlev,nlong+nshort)
    &  aer_tau_lw_vr(:,:,:), & !< longwave aerosol optical depth (absorption only) 
                               !      [layer-1], vertically reversed
    &  aer_tau_sw_vr(:,:,:), & !< shortwave aerosol optical depth (absorption+scattering)
                               !      [layer-1], vertically reversed
    &  aer_piz_sw_vr(:,:,:), & !< shortwave aerosol single scattering albedo
                               !      [-], vertically reversed
    &  aer_cg_sw_vr(:,:,:)     !< shortwave aerosol asymmetry factor 
                               !      [-], vertically reversed
! Local variables
  TYPE(t_mode), POINTER ::   &
    &  this_mode               !< Pointer to current mode
  REAL(wp)              ::   & !< NOTE: No switch for volcanic ash as it is not part of climatology
    &  dust_clim,            & !< Use climatological values for dust optical properties
                               !      (1 = yes, 0 = no)
    &  seas_clim               !< Use climatological values for seasalt optical properties
                               !      (1 = yes, 0 = no)
  REAL(wp), POINTER     ::   &
    &  tracer(:,:,:,:),      & !< Tracer mixing ratio container [mug kg-1]
    &  rho(:,:,:),           & !< Air density [kg m-3]
    &  dz(:,:,:)               !< Model layer height [m]
  REAL(wp),ALLOCATABLE  ::   & 
    &  tau_vr(:,:,:),        & !< mode specific extinction optical depth from ART aerosols
    &  tau_s_vr(:,:,:),      & !< mode specific scattering optical depth from ART aerosols
    &  tauasy_vr(:,:,:),     & !< mode specific SUM(scattering optical depth * asymmetry parameter)
                               !      from ART aerosols
    &  tau_tot_vr(:,:,:),    & !< Combined extinction optical depth from ART aerosols First: [m-1],
                               !      Later: [layer-1]
    &  tau_s_tot_vr(:,:,:),  & !< Combined scattering optical depth from ART aerosols First: [m-1],
                               !      Later: [layer-1]
    &  tau_a_tot_vr(:,:,:),  & !< Combined absorption optical depth from ART aerosols First: [m-1],
                               !      Later: [layer-1]
    &  tauasy_tot_vr(:,:,:), & !< SUM(scattering optical depth * asymmetry parameter) from ART
                               !      aerosols [m-1]
    &  asy_tot_vr(:,:,:)       !< Mean asymmetry parameter from ART aerosols [-]
  REAL(wp) ::                &
    &  z_sum_aea, z_sum_aes    !< absorption/scattering optical depth from climatological 
                               !      + ART aerosols
  INTEGER ::                 &
    &  jspec, jk, jk_vr, jc, & !< loop indices
    &  var_med_dia             !< control variable for varying median parametrization 
                               !  (1=varying median dia,0=constant median dia)

  
  LOGICAL ::                 &
    &  lwarning_flag
  
  ALLOCATE(tau_vr(nproma,ke,(nshort+nlong)))
  ALLOCATE(tau_s_vr(nproma,ke,(nshort+nlong)))
  ALLOCATE(tauasy_vr(nproma,ke,(nshort+nlong)))
  ALLOCATE(tau_tot_vr(nproma,ke,(nshort+nlong)))
  ALLOCATE(tau_s_tot_vr(nproma,ke,(nshort+nlong)))
  ALLOCATE(tau_a_tot_vr(nproma,ke,(nshort+nlong)))
  ALLOCATE(tauasy_tot_vr(nproma,ke,(nshort+nlong))) !< This should be optimized 
                                                    !   -> Only allocate for shortwave bands
  ALLOCATE(asy_tot_vr(nproma,ke,(nshort+nlong)))    !< This should be optimized
                                                    !   -> Only allocate for shortwave bands
 
  lwarning_flag = .FALSE.

  ! Init
  tau_tot_vr    = 0.0_wp
  tau_s_tot_vr  = 0.0_wp
  tau_a_tot_vr  = 0.0_wp
  tauasy_tot_vr = 0.0_wp
  asy_tot_vr    = 0.0_wp

  ! Shortcuts (tracer container, rho and dz)
  tracer    => p_nh_state(jg)%prog(nnew_rcf(jg))%tracer
  rho       => p_nh_state(jg)%prog(nnew(jg))%rho
  dz        => p_nh_state(jg)%metrics%ddqz_z_full
  
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

  DO WHILE(ASSOCIATED(this_mode))
    ! Select type of mode
    SELECT TYPE (fields=>this_mode%fields)
      CLASS IS (t_fields_2mom)
        ! Before radiation, the modal parameters have to be calculated
        CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),  &
          &                     jcs, jce, ks, ke, jb, tracer(:,:,jb,:))
        ! TEMPORARY: hard coded names should be avoided
        SELECT CASE (TRIM(fields%name))
          CASE('dusta','dustb','dustc','asha','ashb','ashc')
            var_med_dia = 1
          CASE DEFAULT
            var_med_dia = 0
        END SELECT
        CALL get_tau_vr(tracer,rho,jb,nlong,nshort,ks,ke,jcs,jce,jg,fields,tau_vr,var_med_dia)
        CALL get_tau_s_vr(jb,nlong,nshort,ks,ke,jcs,jce,fields,tau_vr,tau_s_vr,var_med_dia)
        CALL get_tauasy_vr(jb,nlong,nshort,ks,ke,jcs,jce,fields,tau_s_vr,tauasy_vr,var_med_dia)
        !!CALL get_tau_s_vr(nlong,nshort,ks,ke,jcs,jce,fields,tau_vr,tau_s_vr)    !Original
        !!CALL get_tauasy_vr(nlong,nshort,ks,ke,jcs,jce,fields,tau_s_vr,tauasy_vr) !Original
        ! Calculate total optical depth for all modes
        CALL calc_sum_vr(tau_vr,nlong,nshort,ks,ke,jcs,jce,tau_tot_vr)
        ! Calculate total scattering optical depth for all modes
        CALL calc_sum_vr(tau_s_vr,nlong,nshort,ks,ke,jcs,jce,tau_s_tot_vr)
        ! Calculate total scattering optical depth for all modes
        CALL calc_sum_vr(tauasy_vr,nlong,nshort,ks,ke,jcs,jce,tauasy_tot_vr)
      CLASS DEFAULT
        ! No ARI for monodisperse particles
    END SELECT
    this_mode => this_mode%next_mode
  ENDDO !associated(this_mode)
  
  DO jspec=1,nlong+nshort
    DO jk = ks, ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce
        tau_a_tot_vr(jc,jk_vr,jspec) = tau_tot_vr(jc,jk_vr,jspec) - tau_s_tot_vr(jc,jk_vr,jspec)
        IF (tau_s_tot_vr(jc,jk_vr,jspec)> 1.e-20_wp) THEN
          asy_tot_vr(jc,jk_vr,jspec)   = tauasy_tot_vr(jc,jk_vr,jspec)/tau_s_tot_vr(jc,jk_vr,jspec)
        ELSE
          asy_tot_vr(jc,jk_vr,jspec)   = 0.0_wp
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  
  ! Height factors
  DO jspec=1,nlong+nshort
    DO jk = ks, ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce
        tau_tot_vr(jc,jk_vr,jspec)   = tau_tot_vr(jc,jk_vr,jspec)   * dz(jc,jk,jb) !< tau_EL
        tau_s_tot_vr(jc,jk_vr,jspec) = tau_s_tot_vr(jc,jk_vr,jspec) * dz(jc,jk,jb) !< tau_SL
        tau_a_tot_vr(jc,jk_vr,jspec) = tau_a_tot_vr(jc,jk_vr,jspec) * dz(jc,jk,jb) !< tau_AL
      ENDDO
    ENDDO
  ENDDO

  ! ----------------------------------
  ! --- Longwave
  ! ----------------------------------
  dust_clim   = 1.0_wp
  seas_clim   = 1.0_wp
  
  IF(art_config(jg)%iart_dust > 0)    dust_clim = 0.0_wp
  IF(art_config(jg)%iart_seasalt > 0) seas_clim = 0.0_wp
  ! TODO create similar if statement to check for control run?
  
  DO jspec=1,nlong
    DO jk=ks,ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce
        ! LW opt thickness of aerosols
        aer_tau_lw_vr(jc,jk,jspec) =  zaeq1(jc,jk_vr) * zaea_rrtm(jspec,1) &
          &                         + zaeq2(jc,jk_vr) * zaea_rrtm(jspec,2) * seas_clim &
          &                         + zaeq3(jc,jk_vr) * zaea_rrtm(jspec,3) * dust_clim &
          &                         + zaeq4(jc,jk_vr) * zaea_rrtm(jspec,4) &
          &                         + zaeq5(jc,jk_vr) * zaea_rrtm(jspec,5) &
          &                         + tau_a_tot_vr(jc,jk,jspec)
      ENDDO
    ENDDO
  ENDDO
  
  ! ----------------------------------
  ! --- Shortwave
  ! ----------------------------------
  
  DO jspec=1+nlong,nlong+nshort
    DO jk=ks,ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce

        z_sum_aea = zaeq1(jc,jk_vr) * zaea_rrtm(jspec,1) &
          &       + zaeq2(jc,jk_vr) * zaea_rrtm(jspec,2) * seas_clim &
          &       + zaeq3(jc,jk_vr) * zaea_rrtm(jspec,3) * dust_clim &
          &       + zaeq4(jc,jk_vr) * zaea_rrtm(jspec,4) &
          &       + zaeq5(jc,jk_vr) * zaea_rrtm(jspec,5) &
          &       + tau_a_tot_vr(jc,jk,jspec)
        
        z_sum_aes = zaeq1(jc,jk_vr) * zaes_rrtm(jspec,1) &
          &       + zaeq2(jc,jk_vr) * zaes_rrtm(jspec,2) * seas_clim &
          &       + zaeq3(jc,jk_vr) * zaes_rrtm(jspec,3) * dust_clim &
          &       + zaeq4(jc,jk_vr) * zaes_rrtm(jspec,4) &
          &       + zaeq5(jc,jk_vr) * zaes_rrtm(jspec,5) &
          &       + tau_s_tot_vr(jc,jk,jspec)
        
        IF (z_sum_aes < 1.e-20_wp) lwarning_flag = .TRUE.
        
        ! sw aerosol optical thickness
        aer_tau_sw_vr(jc,jk,jspec-nlong) = z_sum_aea + z_sum_aes
       
        ! sw aerosol single scattering albedo
        aer_piz_sw_vr(jc,jk,jspec-nlong) = z_sum_aes / ( z_sum_aea + z_sum_aes ) 
        
        ! sw aerosol asymmetry factor
        aer_cg_sw_vr(jc,jk,jspec-nlong) =                                             &
          & (   zaeq1(jc,jk_vr) * zaes_rrtm(jspec,1) * zaeg_rrtm(jspec,1)             &
          &   + zaeq2(jc,jk_vr) * zaes_rrtm(jspec,2) * zaeg_rrtm(jspec,2) * seas_clim &
          &   + zaeq3(jc,jk_vr) * zaes_rrtm(jspec,3) * zaeg_rrtm(jspec,3) * dust_clim &
          &   + zaeq4(jc,jk_vr) * zaes_rrtm(jspec,4) * zaeg_rrtm(jspec,4)             &
          &   + zaeq5(jc,jk_vr) * zaes_rrtm(jspec,5) * zaeg_rrtm(jspec,5)             &
          &   + tau_s_tot_vr(jc,jk,jspec)  * asy_tot_vr(jc,jk,jspec) )                &
          &   / z_sum_aes
          
      ENDDO
    ENDDO
  ENDDO
  
  IF (lwarning_flag) THEN
    CALL message('mo_art_radiation_aero:art_radiation_aero','WARNING: z_sum_aes very small')
  ENDIF
  
  DEALLOCATE(tau_vr)
  DEALLOCATE(tau_s_vr)
  DEALLOCATE(tauasy_vr)
  DEALLOCATE(tau_tot_vr)
  DEALLOCATE(tau_s_tot_vr)
  DEALLOCATE(tau_a_tot_vr)
  DEALLOCATE(tauasy_tot_vr)
  DEALLOCATE(asy_tot_vr)

END SUBROUTINE art_radiation_aero
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_tau_vr(tracer,rho,jb,nlong,nshort,ks,ke,jcs,jce,jg,fields,tau_vr,var_med_dia)
!<
! SUBROUTINE get_tau_vr
! Calculate vertically reversed (i.e. k-loop) optical thickness
! Based on: Philipp Gasch (2016): Numerical simulations of an exceptional dust event
!                                 in the Eastern Mediterranean including the mineral
!                                 dust radiative feedback.
!                                 Master Thesis at the Fakultaet fuer Physik of 
!                                 the Karlsruhe Institute of Technology (KIT)
! Part of Module: mo_art_radiation_aero
! Author: Carolin Walter, KIT, Philipp Gasch, KIT, Daniel Rieger, KIT
! Initial Release: 2015-11-12
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  REAL(wp),INTENT(in)    ::  &
    &  tracer(:,:,:,:),      & !< Tracer mixing ratios [mug kg-1]
    &  rho(:,:,:)              !< Air Density [kg m-3]
  INTEGER,INTENT(in)     ::  &
    &  jb,                   & !< Block index
    &  nlong, nshort,        & !< Number of longwave/shortwave bands
    &  ks, ke, jcs, jce ,    & !< vertical loop start/end index, jc start/end index
    &  jg,                   &        
    &  var_med_dia             !< control variable for varying median parametrization 
                               !  (1=varying median dia,0=constant median dia)
  TYPE(t_fields_2mom),INTENT(in) :: &
    &  fields                  !< All fields associated to 2-mom mode
  REAL(wp),INTENT(inout) ::  &
    &  tau_vr(:,:,:)           !< 

  ! Local variables
  REAL(wp),ALLOCATABLE   ::  &
    & loc_ext_coeff(:)
  INTEGER ::                 &
    &  jspec, jk, jk_vr,     & !< Loop indices
    &  jc, jsp,              & !< Loop indices
    &  tr_idx,               &
    &  ierror
  CHARACTER(LEN=IART_VARNAMELEN) ::  &
    &  tracer_str

  ALLOCATE(loc_ext_coeff(nproma))

  tau_vr(:,:,:) = 0._wp
  ierror = 1

  SELECT CASE(TRIM(fields%name))
    CASE('sol_ait')
      tracer_str = 'so4_sol_ait'
    CASE('sol_acc')
      tracer_str = 'so4_sol_acc'
    CASE('insol_acc')
      tracer_str = 'ash_insol_acc'
    CASE('insol_coa')
      tracer_str = 'ash_insol_coa'
    CASE('mixed_acc')
      tracer_str = 'ash_mixed_acc'
    CASE('mixed_coa')
      tracer_str = 'ash_mixed_coa'
    CASE('giant')
      tracer_str = 'ash_giant'
    CASE DEFAULT
      tracer_str = ''
  END SELECT

  IF(tracer_str/='') THEN
    CALL p_art_data(jg)%dict_tracer%get(tracer_str,tr_idx,ierror)
  END IF

  DO jspec=1,nlong+nshort
    DO jk = ks, ke
      jk_vr = ke+1-jk
      IF (var_med_dia==1) THEN
!NEC$ ivdep
        DO jc = jcs,jce
          IF (fields%diameter(jc,jk,jb) > fields%info%diameter_ini_nmb * &
            &                             fields%opt_props%dia_max_factor(1)) THEN
            loc_ext_coeff(jc) = fields%opt_props%ext_default(jspec,1)           
          ELSEIF (fields%diameter(jc,jk,jb) < fields%info%diameter_ini_nmb * &
              &                               fields%opt_props%dia_min_factor(1)) THEN          
            loc_ext_coeff(jc) = fields%opt_props%ext_default(jspec,2)               
          ELSE
            loc_ext_coeff(jc) = fields%opt_props%ext_param(jspec,1) * (fields%diameter(jc,jk,jb)  &
              &               * 1.e9_wp)**3                                                       &
              &               + fields%opt_props%ext_param(jspec,2) * (fields%diameter(jc,jk,jb)  &
              &               * 1.e9_wp)**2                                                       &
              &               + fields%opt_props%ext_param(jspec,3) * (fields%diameter(jc,jk,jb)  &
              &               * 1.e9_wp)**1                                                       &
              &               + fields%opt_props%ext_param(jspec,4)                    
          ENDIF
        ENDDO
      ELSE !> var_med_dia==0
!NEC$ ivdep
        DO jc = jcs,jce
          loc_ext_coeff(jc) = fields%opt_props%ext_coeff(jspec) 
        ENDDO
      ENDIF

      IF(ierror==SUCCESS) THEN  ! IF one of above mentioned modes/tracers are used
!NEC$ ivdep
        DO jc = jcs,jce
          tau_vr(jc,jk_vr,jspec) = tau_vr(jc,jk_vr,jspec)                   &
            &                    + (tracer(jc,jk,jb,tr_idx)                 &
            &                    * loc_ext_coeff(jc)                        &
            &                    * rho(jc,jk,jb) * 1.e-6_wp)           
        END DO 
      ELSE
        DO jsp = 1, fields%ntr-1
!NEC$ ivdep
          DO jc = jcs,jce
            tau_vr(jc,jk_vr,jspec) = tau_vr(jc,jk_vr,jspec)                  &
              &                    + (tracer(jc,jk,jb,fields%itr3(jsp))      &
              &                    * loc_ext_coeff(jc)                       &
              &                    * rho(jc,jk,jb) * 1.e-6_wp)
          ENDDO
        ENDDO
      END IF
      
    ENDDO
  ENDDO
  
  DEALLOCATE(loc_ext_coeff)

END SUBROUTINE get_tau_vr
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_tau_s_vr(jb,nlong,nshort,ks,ke,jcs,jce,fields,tau_vr,tau_s_vr, var_med_dia)
!<
! SUBROUTINE get_tau_s_vr
! Calculate vertically reversed (i.e. k-loop) scattering optical thickness
! Based on: Philipp Gasch (2016): Numerical simulations of an exceptional dust event
!                                 in the Eastern Mediterranean including the mineral
!                                 dust radiative feedback.
!                                 Master Thesis at the Fakultaet fuer Physik of 
!                                 the Karlsruhe Institute of Technology (KIT)
! Part of Module: mo_art_radiation_aero
! Author: Carolin Walter, KIT, Philipp Gasch, KIT, Daniel Rieger, KIT
! Initial Release: 2015-11-12
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  INTEGER,INTENT(in)     ::  &
    &  jb,                   & !< Block index
    &  nlong, nshort,        & !< Number of longwave/shortwave bands
    &  ks, ke, jcs, jce,     & !< vertical loop start/end index, jc start/end index
    &  var_med_dia             !< control variable for varying median parametrization 
                               !  (1=varying median dia,0=constant median dia)
  TYPE(t_fields_2mom),INTENT(in) :: &
    &  fields                  !< All fields associated to 2-mom mode
  REAL(wp),INTENT(in)    ::  &
    &  tau_vr(:,:,:)
  REAL(wp),INTENT(inout) ::  &
    &  tau_s_vr(:,:,:)         !<
  ! Local variables
  REAL(wp)               ::  &
    & loc_ssa_coeff
  INTEGER ::                 &
    &  jspec, jc, jk, jk_vr    !< Loop indices
  
  tau_s_vr(:,:,:) = 0._wp

  DO jspec=1,nlong+nshort
    DO jk = ks, ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce
        IF (var_med_dia==1) THEN
          IF (fields%diameter(jc,jk,jb) > fields%info%diameter_ini_nmb                   &
            &                           * fields%opt_props%dia_max_factor(1)) THEN
            loc_ssa_coeff = fields%opt_props%ssa_default(jspec,1)
          ELSEIF (fields%diameter(jc,jk,jb) < fields%info%diameter_ini_nmb               &
            &                               * fields%opt_props%dia_min_factor(1)) THEN
            loc_ssa_coeff = fields%opt_props%ssa_default(jspec,2)
          ELSE
            loc_ssa_coeff = fields%opt_props%ssa_param(jspec,1) * (fields%diameter(jc,jk,jb)       &
              &           * 1.e9_wp)**3                                                            &
              &           + fields%opt_props%ssa_param(jspec,2) * (fields%diameter(jc,jk,jb)       &
              &           * 1.e9_wp)**2                                                            &
              &           + fields%opt_props%ssa_param(jspec,3) * (fields%diameter(jc,jk,jb)       &
              &           * 1.e9_wp)**1                                                            &
              &           +  fields%opt_props%ssa_param(jspec,4)
          ENDIF
        ELSE ! var_med_dia == 0
          loc_ssa_coeff=fields%opt_props%ssa_coeff(jspec)
        ENDIF
        tau_s_vr(jc,jk_vr,jspec) = tau_vr(jc,jk_vr,jspec) * loc_ssa_coeff
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE get_tau_s_vr
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_tauasy_vr(jb,nlong,nshort,ks,ke,jcs,jce,fields,tau_s_vr,tauasy_vr,var_med_dia)
!<
! SUBROUTINE get_tauasy_vr
! Calculate vertically reversed (i.e. k-loop) scattering opt. thickness * assy. param.
! Based on: Philipp Gasch (2016): Numerical simulations of an exceptional dust event
!                                 in the Eastern Mediterranean including the mineral
!                                 dust radiative feedback.
!                                 Master Thesis at the Fakultaet fuer Physik of 
!                                 the Karlsruhe Institute of Technology (KIT)
! Part of Module: mo_art_radiation_aero
! Author: Carolin Walter, KIT, Philipp Gasch, KIT, Daniel Rieger, KIT
! Initial Release: 2015-11-12
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  INTEGER,INTENT(in)     ::  &
    &  jb,                   & !< Block index
    &  nlong, nshort,        & !< Number of longwave/shortwave bands
    &  ks, ke, jcs, jce ,    & !< vertical loop start/end index, jc start/end index
    &  var_med_dia             !< control variable for varying median parametrization 
                               !  (1=varying median dia,0=constant median dia)
  TYPE(t_fields_2mom),INTENT(in) :: &
    &  fields                  !< All fields associated to 2-mom mode
  REAL(wp),INTENT(in)    ::  &
    &  tau_s_vr(:,:,:)         !<
  REAL(wp),INTENT(inout) ::  &
    &  tauasy_vr(:,:,:)        !<
  ! Local variables
  REAL(wp)               :: &
    & loc_asy_coeff
  INTEGER ::                 &
    &  jspec, jk, jk_vr, jc    !< Loop indices

  tauasy_vr(:,:,:) = 0._wp

  DO jspec=nlong+1,nlong+nshort !< Only shortwave bands
    DO jk = ks, ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce
        IF (var_med_dia==1) THEN
          IF (fields%diameter(jc,jk,jb) > fields%info%diameter_ini_nmb                   &
            &                           * fields%opt_props%dia_max_factor(1)) THEN
            loc_asy_coeff = fields%opt_props%asy_default(jspec,1)
          ELSEIF (fields%diameter(jc,jk,jb) < fields%info%diameter_ini_nmb               &
            &                               * fields%opt_props%dia_min_factor(1)) THEN
            loc_asy_coeff = fields%opt_props%asy_default(jspec,2)
          ELSE
            loc_asy_coeff = fields%opt_props%asy_param(jspec,1) * (fields%diameter(jc,jk,jb)      &
              &           * 1.e9_wp)**3                                                           &
              &           + fields%opt_props%asy_param(jspec,2) * (fields%diameter(jc,jk,jb)      &
              &           * 1.e9_wp)**2                                                           &
              &           + fields%opt_props%asy_param(jspec,3) * (fields%diameter(jc,jk,jb)      &
              &           * 1.e9_wp)**1                                                           &
              &           + fields%opt_props%asy_param(jspec,4)
          ENDIF
        ELSE ! var_med_dia == 0
          loc_asy_coeff=fields%opt_props%asy_coeff(jspec)
        ENDIF
        tauasy_vr(jc,jk_vr,jspec) = tau_s_vr(jc,jk_vr,jspec) * loc_asy_coeff
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE get_tauasy_vr
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE calc_sum_vr(factor_mode_vr,nlong,nshort,ks,ke,jcs,jce,factor_tot_vr)
!<
! SUBROUTINE calc_sum_vr
! Calculate vertical sum of a given field
! Based on: Philipp Gasch (2016): Numerical simulations of an exceptional dust event
!                                 in the Eastern Mediterranean including the mineral
!                                 dust radiative feedback.
!                                 Master Thesis at the Fakultaet fuer Physik of 
!                                 the Karlsruhe Institute of Technology (KIT)
! Part of Module: mo_art_radiation_aero
! Author: Carolin Walter, KIT, Philipp Gasch, KIT, Daniel Rieger, KIT
! Initial Release: 2015-11-12
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  REAL(wp),INTENT(in)    ::  &
    &  factor_mode_vr(:,:,:)   !< Factor to add to sum
  INTEGER,INTENT(in)     ::  &
    &  nlong, nshort,        & !< Number of longwave/shortwave bands
    &  ks, ke, jcs, jce        !< vertical loop start/end index, jc start/end index
  REAL(wp),INTENT(inout) ::  &
    &  factor_tot_vr(:,:,:)    !< Sum
  ! Local variables
  INTEGER ::                 &
    &  jspec, jk, jk_vr, jc    !< Loop indices
  
  DO jspec=1,nlong+nshort
    DO jk = ks, ke
      jk_vr = ke+1-jk
!NEC$ ivdep
      DO jc = jcs,jce
        factor_tot_vr(jc,jk_vr,jspec) = factor_tot_vr(jc,jk_vr,jspec) &
          &                           + factor_mode_vr(jc,jk_vr,jspec)
      ENDDO
    ENDDO
  ENDDO
  
END SUBROUTINE calc_sum_vr
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_radiation_aero
