!
! mo_art_diagnostics
! This module provides routines to calculate diagnostic output variables.
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

MODULE mo_art_diagnostics
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_impl_constants,                ONLY: itfastphy, SUCCESS, MAX_CHAR_LENGTH
  USE mo_statistics,                    ONLY: time_avg
  USE mo_tropopause,                    ONLY: WMO_tropopause
  USE mo_aes_wmo_config,                ONLY: aes_wmo_config
  USE mo_fortran_tools,                 ONLY: set_acc_host_or_device

! ART
  USE mo_art_config,                    ONLY: npreslay
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio
  USE mo_art_data,                      ONLY: t_art_data, p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_clipping,                  ONLY: art_clip_lt
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_diag_types,                ONLY: art_diag_tracer_index, &
                                          &   t_art_aeronet, t_art_max_layer
  USE mo_art_impl_constants,            ONLY: IART_ACC_EMISS, IART_EMISS, IART_ACC_WETDEPO_GSCP, &
                                          &   IART_ACC_WETDEPO_CON, IART_ACC_WETDEPO_RRSFC
  USE mo_physical_constants,             ONLY: argas, avo, p0sl_bg
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_dust_diagnostics
  PUBLIC :: art_seas_diagnostics
  PUBLIC :: art_volc_diagnostics
  PUBLIC :: art_soot_diagnostics
  PUBLIC :: art_radio_diagnostics
  PUBLIC :: art_radio_diagnostics_dt_phy
  PUBLIC :: art_chem_calc_column
  PUBLIC :: art_calc_tropopause
  PUBLIC :: art_calc_o2_column
    
  PUBLIC :: art_save_aerosol_wet_deposition
  PUBLIC :: art_save_aerosol_emission

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_dust_diagnostics( rho, pres, p_trac, dz,                &
  &                              istart, iend, nlev, jb, p_art_data )

  !<
  ! SUBROUTINE art_dust_diagnostics
  ! integrates the 3d AOD for art-tracers
  ! vertically and saves the 2d AOD of art-tracers
  ! in addition it calculates the
  ! total mineral dust mass concentration [kg/m3]
  !
  ! Based on: -
  ! Part of Module: mo_art_diagnostics
  ! Author: Vanessa Bachmann, DWD
  ! Initial Release: 2022-06-09
  ! Modifications:
  ! 2022-06-09: -
  ! -
  !>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density  [kg/m3]
    &  pres(:,:),           & !< Air pressure [Pa]
    &  p_trac(:,:,:),       & !< Tracer mixing ratios [mug/kg]
    &  dz(:,:)                !< Layer thickness
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data), INTENT(inout), TARGET :: &
    &  p_art_data             !< Data container for ART
! Local variables
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_dust_diagnostics"

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of aerosol optical depth
  ! ----------------------------------------------------------------------
  CALL art_calc_tau_vi( istart, iend, nlev, jb, &
    &                   p_art_data%diag%dust_aeronet )

  ! -------------------------------------------------------------
  ! --- Calculate total mineral dust mass concentration [kg/m3]
  ! -------------------------------------------------------------
  CALL art_calc_total_mc( rho, p_trac, istart, iend, nlev, jb, p_art_data, &
    &                     'dust', p_art_data%diag%dust_total_mc )

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of mineral dust mass concentration [kg/m2]
  ! ----------------------------------------------------------------------
  CALL art_calc_total_mc_vi( dz, istart, iend, nlev, jb,        &
    &                        p_art_data%diag%dust_total_mc,     &
    &                        p_art_data%diag%dust_total_mc_vi )
  
  ! -----------------------------------------------------------------------------------------------
  ! -- Diagnose maximum total mineral dust mass concentration between given pressure levels [kg/m3]
  ! -----------------------------------------------------------------------------------------------
  CALL art_diag_max_preslay( pres, istart, iend, nlev, jb,              &
    &                        p_art_data, p_art_data%diag%dust_total_mc, &
    &                        p_art_data%diag%dust_max_total_mc )

END SUBROUTINE art_dust_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_seas_diagnostics( rho, p_trac,     &
  &                              istart, iend, nlev, jb, p_art_data )

  !<
  ! SUBROUTINE art_seas_diagnostics
  ! integrates the 3d AOD for art-tracers
  ! vertically and saves the 2d AOD of art-tracers
  ! in addition it calculates the
  ! total sea salt mass concentration [kg/m3]
  !
  ! Based on: -
  ! Part of Module: mo_art_diagnostics
  ! Author: Vanessa Bachmann, DWD
  ! Initial Release: 2022-09-01
  ! Modifications:
  ! 2022-09-01: -
  ! -
  !>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density  [kg/m3]
    &  p_trac(:,:,:)          !< Tracer mixing ratios [mug/kg]
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data), INTENT(inout), TARGET :: &
    &  p_art_data             !< Data container for ART
! Local variables
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_seas_diagnostics"

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of aerosol optical depth
  ! ----------------------------------------------------------------------
  CALL art_calc_tau_vi( istart, iend, nlev, jb, &
    &                   p_art_data%diag%seas_aeronet )

  ! -------------------------------------------------------------
  ! --- Calculate total sea salt mass concentration [kg/m3]
  ! -------------------------------------------------------------
  CALL art_calc_total_mc( rho, p_trac, istart, iend, nlev, jb, p_art_data, &
    &                     'seas', p_art_data%diag%seas_total_mc )

END SUBROUTINE art_seas_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_soot_diagnostics( rho, p_trac,     &
  &                              istart, iend, nlev, jb, p_art_data )

  !<
  ! SUBROUTINE art_soot_diagnostics
  ! integrates the 3d AOD for art-tracers
  ! vertically and saves the 2d AOD of art-tracers
  ! in addition it calculates the
  ! total biomass burning aerosol mass concentration [kg/m3]
  !
  ! Based on: -
  ! Part of Module: mo_art_diagnostics
  ! Author: Nikolas Porz, DWD
  ! Initial Release: 2023-12-13
  ! Modifications:
  ! 2023-12-13: -
  ! -
  !>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density  [kg/m3]
    &  p_trac(:,:,:)          !< Tracer mixing ratios [mug/kg]
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data), INTENT(inout), TARGET :: &
    &  p_art_data             !< Data container for ART
! Local variables
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_soot_diagnostics"

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of aerosol optical depth
  ! ----------------------------------------------------------------------
  CALL art_calc_tau_vi( istart, iend, nlev, jb, &
    &                   p_art_data%diag%soot_aeronet )

  ! -------------------------------------------------------------
  ! --- Calculate total biomass burning aerosol mass concentration [kg/m3]
  ! -------------------------------------------------------------
  CALL art_calc_total_mc( rho, p_trac, istart, iend, nlev, jb, p_art_data, &
    &                     'soot', p_art_data%diag%soot_total_mc )
  ! note the difference between 'biomass_burning' and 'soot' for emission scheme
END SUBROUTINE art_soot_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_volc_diagnostics(rho, pres, p_trac, dz, hml,           &
  &                             istart, iend, nlev, jb, p_art_data,   &
  &                             ash_scheme)
!<
! SUBROUTINE art_volc_diagnostics
! Diagnostics for volcanic ash
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2015-09-16
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density  [kg/m3]
    &  pres(:,:),           & !< Air pressure [Pa]
    &  p_trac(:,:,:),       & !< Tracer mixing ratios [mug/kg]
    &  dz(:,:),             & !< Layer thickness
    &  hml(:,:)               !< Height of main levels
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data), INTENT(inout), TARGET :: &
    &  p_art_data             !< Data container for ART
  INTEGER, INTENT(in)    :: &
    &  ash_scheme             !< 1 for bins, 2 for modes
! Local variables
  REAL(wp)               :: &
    &  ash_total_mix,       & !< total ash mixing ratio [ug/kg]
    &  pres_bot, pres_top
  INTEGER     ::            &
    &  jc, jk, jt             !< Loop indices (nproma, nlev, n_ash_cld_thr)
  LOGICAL ::                &
    &  l_ash_mc_mask(istart:iend,1:nlev)
  INTEGER ::                &
    &  i_maxloc(1:nlev),    &
    &  iash1, iash2,        & !< Tracer container indices
    &  iash3, iash4,        & !< Tracer container indices
    &  iash5, iash6,        & !< Tracer container indices
    &  ierror                 !< Error return value
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_volc_diagnostics"

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of aerosol optical depth
  ! ----------------------------------------------------------------------
  CALL art_calc_tau_vi( istart, iend, nlev, jb, &
    &                   p_art_data%diag%volc_aeronet )

  ! -------------------------------------------------------------
  ! --- Calculate total volcanic ash mass concentration [kg/m3]
  ! -------------------------------------------------------------  
  SELECT CASE(ash_scheme)
    CASE(1) ! for bins
      CALL p_art_data%dict_tracer%get('ash1',iash1,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash1 not found in dictionary.')
      CALL p_art_data%dict_tracer%get('ash2',iash2,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash2 not found in dictionary.')
      CALL p_art_data%dict_tracer%get('ash3',iash3,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash3 not found in dictionary.')
      CALL p_art_data%dict_tracer%get('ash4',iash4,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash4 not found in dictionary.')
      CALL p_art_data%dict_tracer%get('ash5',iash5,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash5 not found in dictionary.')
      CALL p_art_data%dict_tracer%get('ash6',iash6,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash6 not found in dictionary.')

      DO jk = 1, nlev
!NEC$ ivdep
        DO jc = istart, iend
          ash_total_mix =  &
            &    p_trac(jc,jk,iash1) + p_trac(jc,jk,iash2) + p_trac(jc,jk,iash3)  &
            &  + p_trac(jc,jk,iash4) + p_trac(jc,jk,iash5) + p_trac(jc,jk,iash6)
          p_art_data%diag%ash_total_mc(jc,jk,jb) =  &
            &    ash_total_mix * rho(jc,jk) * 1.e-9_wp
        END DO !jc
      END DO !jk
      ! clipping of negative mass conc. values
      CALL art_clip_lt(p_art_data%diag%ash_total_mc(:,:,jb),0.0_wp)

    CASE(2) ! for modes
      CALL art_calc_total_mc(rho, p_trac, istart, iend, nlev, jb, p_art_data, &
        &                    'volc', p_art_data%diag%ash_total_mc)

    CASE DEFAULT
      CALL finish (thisroutine, 'iart_volcano needs to be 1 or 2 for volcanic ash diagnostics.')
  END SELECT

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of volcanic ash mass concentration [kg/m2]
  ! ----------------------------------------------------------------------
  CALL art_calc_total_mc_vi( dz, istart, iend, nlev, jb,       &
    &                        p_art_data%diag%ash_total_mc,     &
    &                        p_art_data%diag%ash_total_mc_vi )
  
  ! -----------------------------------------------------------------------------------------------
  ! -- Diagnose maximum total volcanic ash mass concentration between given pressure levels [kg/m3]
  ! -----------------------------------------------------------------------------------------------
  CALL art_diag_max_preslay( pres, istart, iend, nlev, jb,             &
    &                        p_art_data, p_art_data%diag%ash_total_mc, &
    &                        p_art_data%diag%ash_max_total_mc )

  ! --------------------------------------------------------------------------
  ! --- Diagnose height of maximal total volcanic ash mass concentration [m]
  ! --------------------------------------------------------------------------
  DO jk = 1, nlev
!NEC$ ivdep
    DO jc = istart, iend
      IF ( p_art_data%diag%ash_total_mc(jc,jk,jb) >= 200.0_wp * 1.e-9_wp ) THEN
        l_ash_mc_mask(jc,jk) = .TRUE.
      ELSE
        l_ash_mc_mask(jc,jk) = .FALSE.
      ENDIF
    END DO !jc
  END DO !jk
!NEC$ ivdep
  DO jc = istart, iend
    IF ( ANY( l_ash_mc_mask(jc,:) ) ) THEN
      i_maxloc(:) = MAXLOC( p_art_data%diag%ash_total_mc(jc,1:nlev,jb), l_ash_mc_mask(jc,:) )
      p_art_data%diag%ash_hml_max(jc,jb) = hml(jc,i_maxloc(1))
    ELSE
      p_art_data%diag%ash_hml_max(jc,jb) = 0.0_wp
    ENDIF
  END DO !jc

  ! --------------------------------------------------------------------------
  ! --- Diagnose height of volcanic ash cloud base and cloud top  [m] 
  ! --------------------------------------------------------------------------
  DO jt=1,2 ! loop over different ash conc. thresholds for cloud def.
    DO jk = 1, nlev
!NEC$ ivdep
      DO jc = istart, iend
        IF ( p_art_data%diag%ash_total_mc(jc,jk,jb) >=  &
          &  p_art_data%diag%ash_cloud(jt)%threshold * 1.e-9_wp ) THEN
          l_ash_mc_mask(jc,jk) = .TRUE.
        ELSE
          l_ash_mc_mask(jc,jk) = .FALSE.
        ENDIF
      END DO !jc
    END DO !jk
!NEC$ ivdep
    DO jc = istart, iend
      IF ( ANY( l_ash_mc_mask(jc,:) ) ) THEN
        p_art_data%diag%ash_cloud(jt)%base(jc,jb) = MINVAL(hml(jc,:),l_ash_mc_mask(jc,:))
        p_art_data%diag%ash_cloud(jt)%top(jc,jb)  = MAXVAL(hml(jc,:),l_ash_mc_mask(jc,:))
      ELSE
        p_art_data%diag%ash_cloud(jt)%base(jc,jb) = 0.0_wp
        p_art_data%diag%ash_cloud(jt)%top(jc,jb)  = 0.0_wp
      ENDIF
    END DO !jc
  END DO !jt

END SUBROUTINE art_volc_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_radio_diagnostics( jg, pres, p_trac, istart, iend, nlev, jb, p_art_data )
!<
! SUBROUTINE art_radio_diagnostics
! Diagnostics for radionuclides called for output
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2023-08-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  INTEGER,INTENT(in)     :: &
    &  jg                     !< patch id
  REAL(wp), INTENT(in)   :: &
    &  pres(:,:),           & !< Air pressure [Pa]
    &  p_trac(:,:,:)          !< Tracer mixing ratios [mug/kg]
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data)       :: &
    &  p_art_data             !< Data container for ART

  !Local variables
  ! ---------------------------------------
  TYPE(t_mode), POINTER  :: &
    &  this_mode              !< Pointer to current mode
  INTEGER                :: &
    &  jc, jk, iradioact      !< Loop indices (art_atmo%nproma, nlev, tracer)

  iradioact = 0

  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

  DO WHILE(ASSOCIATED(this_mode))
    ! Select type of mode
    SELECT TYPE (fields=>this_mode%fields)
      TYPE IS (t_fields_radio)
        iradioact = iradioact + 1
        ! ----------------------------------------------------------------------------------------------------
        ! -- Diagnose maximal maximum air concentration in time interval between given pressure levels [Bq/m3]
        ! ----------------------------------------------------------------------------------------------------
        CALL art_diag_max_preslay( pres, istart, iend, nlev, jb,                            &
          &                        p_art_data, p_art_data%diag%radioact(iradioact)%maxtint, &
          &                        p_art_data%diag%radioact(iradioact)%maxtint_layer )
    END SELECT
    this_mode => this_mode%next_mode
  END DO

END SUBROUTINE art_radio_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_radio_diagnostics_dt_phy( jg, rho, p_trac, istart, iend, nlev, jb, p_art_data, &
  &                                      dt_phy_jg, p_sim_time )
!<
! SUBROUTINE art_radio_diagnostics_dt_phy
! Diagnostics for radionuclides called every physics time step
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2017-06-28
! Modifications:
! 2023-08-30: Jochen Foerstner, DWD
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  INTEGER,INTENT(in)     :: &
    &  jg                           !< patch id
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density  [kg/m3]
    &  p_trac(:,:,:)          !< Tracer mixing ratios [mug/kg]
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data)       :: &
    &  p_art_data             !< Data container for ART
  REAL(wp), INTENT(IN)   :: &
    &  dt_phy_jg(:)           !< time interval for all physics packages on domain jg
  REAL(wp), INTENT(IN)   :: &
    &  p_sim_time

  !Local variables
  ! ---------------------------------------
  TYPE(t_mode), POINTER  :: &
    &  this_mode              !< Pointer to current mode
  REAL(wp)               :: &
    &  t_wgt                  !< weight for running time average
  INTEGER                :: &
    &  jc, jk, iradioact      !< Loop indices (art_atmo%nproma, nlev, tracer)

  ! ------------------------------------------------------------------
  ! --- Calculate averaged air concentration of radionuclides [Bq/m3]
  ! ------------------------------------------------------------------

  ! time average weight
  t_wgt = dt_phy_jg(itfastphy) / MAX(1.e-6_wp, p_sim_time)

  iradioact = 0

  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

  IF ( p_sim_time <= 1.e-6_wp) THEN

    DO WHILE(ASSOCIATED(this_mode))
      ! Select type of mode
      SELECT TYPE (fields=>this_mode%fields)
        TYPE IS (t_fields_radio)
          iradioact = iradioact + 1
          DO jk = 1, nlev
!NEC$ ivdep
            DO jc = istart, iend
              p_art_data%diag%radioact(iradioact)%avg(jc,jk,jb) = &
                &  p_trac(jc,jk,fields%itr)*rho(jc,jk)
              p_art_data%diag%radioact(iradioact)%maxtint(jc,jk,jb) = &
                &  p_trac(jc,jk,fields%itr)*rho(jc,jk)
            END DO
          END DO
      END SELECT
      this_mode => this_mode%next_mode
    END DO

  ELSE

    DO WHILE(ASSOCIATED(this_mode))
      ! Select type of mode
      SELECT TYPE (fields=>this_mode%fields)
        TYPE IS (t_fields_radio)
          iradioact = iradioact + 1
          DO jk = 1, nlev
!NEC$ ivdep
            DO jc = istart, iend
              p_art_data%diag%radioact(iradioact)%avg(jc,jk,jb) = &
                &  time_avg( p_art_data%diag%radioact(iradioact)%avg(jc,jk,jb), &
                &            p_trac(jc,jk,fields%itr)*rho(jc,jk), t_wgt )
              p_art_data%diag%radioact(iradioact)%maxtint(jc,jk,jb) = &
                &  MAX( p_art_data%diag%radioact(iradioact)%maxtint(jc,jk,jb),  &
                &       p_trac(jc,jk,fields%itr)*rho(jc,jk) )
            END DO
          END DO
      END SELECT
      this_mode => this_mode%next_mode
    END DO

  ENDIF
  
END SUBROUTINE art_radio_diagnostics_dt_phy
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_chem_calc_column(column, pres, temp, z_ifc, p_tracer, istart, iend, nlev)
!<
! SUBROUTINE art_chem_calc_column
! Diagnostics for chemical species
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jennifer Schroeter, KIT
! Initial Release: 2017-05-11
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>

  REAL(wp),              INTENT(INOUT) ::  column(:,:)          !< calculated column [mol/m2]
  REAL(wp),              INTENT(IN)    ::  pres(:,:)            !< pressure
  REAL(wp),              INTENT(IN)    ::  temp(:,:)            !< temperatur
  REAL(wp),              INTENT(IN)    ::  z_ifc(:,:)           !< geometric height
  REAL(wp),              INTENT(IN)    ::  p_tracer(:,:)        !< tracer vmr

  INTEGER,               INTENT(IN)    :: istart, iend, nlev    !< loop parameters 
  INTEGER                              :: jc, jk

!NEC$ ivdep
  DO jc = istart, iend
    column(jc, 1) = (p_tracer(jc,1) * pres(jc, 1) / (argas * temp(jc, 1)))                     &
      &           * (z_ifc(jc, 1) - z_ifc(jc, 2))
  END DO

  DO jk=2,nlev
!NEC$ ivdep
    DO jc = istart, iend

      column(jc,jk) = (p_tracer(jc,jk) * pres(jc, jk) / (argas * temp(jc, jk)))                &
        &           * (z_ifc(jc, jk) - z_ifc(jc,jk+1))                                               &
        &           + column(jc, jk-1)
    ENDDO
  ENDDO


END SUBROUTINE
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_save_aerosol_wet_deposition(p_art_data, wash_rate,    &
  &           iart_aero_washout, rr_sfc, dz, dtime, idx_tracer, jb,  &
  &           istart, iend, kstart, nlev, third_mom, tracer)
  !<
  ! SUBROUTINE save_art_aerosol_wet_deposition
  ! integrates the 3d washout rates for aerosols vertically
  ! and saves the 2d wet deposition of aerosols
  !
  ! acc_wetdepo: is the total washout of aerosols in the column
  ! acc_wetdepo_rrsfc: is the total washout in the column but 
  !                    only if precipitation reaches the surface
  !
  ! Based on: -
  ! Part of Module: mo_art_diagnostics
  ! Author: Andrea Steiner, DWD
  ! Initial Release: 2018-01-02
  ! Modifications:
  ! 2018-01-02: -
  ! -
  !>
  TYPE(t_art_data),    INTENT(inout) :: &
    & p_art_data                          !< ART data container
  REAL(wp),            INTENT(in)    :: &
    & wash_rate(:,:),                   & !< Washout rate with respect to mass
                                          !  concentration [m3 m-3 s-1]
    & rr_sfc(:),                        & !< grid scale and convective precipitation rate 
                                          !  at the surface
    & dtime,                            & !< Time step (fast physics)
    & dz(:,:,:),                        & !< Vertical layer thickness (m)
    & third_mom(:,:),                   & !< third moment of considered mode (jc,jk) [m3 kg-1]
    & tracer(:,:)                         !< individual tracer concentration (jc,jk) [ug kg-1]
  INTEGER,             INTENT(in)    :: &
    & iart_aero_washout,                & !< Treatment of washout through grid scale and convective
                                          !  precipitation: 0:gscp+con;1:gscp,con;2:gscp,rcucov*con
    & idx_tracer,                       & !< Index of tracer in tracer container
    & jb,                               & !< Counter for block loop
    & istart, iend,                     & !< Start and end of art_atmo%nproma loop
    & kstart,                           & !< Index of level to start vertical loop
    & nlev                                !< Number of levels (equals index of lowest full level)

  !Local variables
  INTEGER :: idx_diag   !< Index of tracer in diagnostics container
  INTEGER :: jc,jk      !< Loop indices




  IF (iart_aero_washout < 100) THEN   ! "first call"

     !------ ACC_WETDEPO_GSCP:

     ! get index in diagnostics container (given the index in the tracer container)
     idx_diag = art_diag_tracer_index(IART_ACC_WETDEPO_GSCP, idx_tracer)

     ! save wet deposition
     IF ( idx_diag > 0 ) THEN
        ! compute accumulated wet deposition of art-tracer
        ! wet deposition shall be positive, thus substract negative washout rate for accumulation
       DO jk = kstart, nlev !start from kstart, as in art_aerosol_washout kstart = 15
!NEC$ ivdep
         DO jc = istart, iend
           p_art_data%diag%acc_wetdepo_gscp(jc,jb,idx_diag) =      &
             & p_art_data%diag%acc_wetdepo_gscp(jc,jb,idx_diag)    &
             & - wash_rate(jc,jk) * dtime / third_mom(jc,jk)       &
             & * tracer(jc,jk) * 1.e-9_wp * dz(jc,jk,jb)
         !Note: *rho not necessray anymore (see e.g. art_integrate_explicit)
         ENDDO !jc
       ENDDO !jk
     END IF

  ELSE 

     !------ ACC_WETDEPO_CON:

     ! get index in diagnostics container (given the index in the tracer container)
     idx_diag = art_diag_tracer_index(IART_ACC_WETDEPO_CON, idx_tracer)

     ! save wet deposition
     IF ( idx_diag > 0 ) THEN
        ! compute accumulated wet deposition of art-tracer
        ! wet deposition shall be positive, thus substract negative washout rate for accumulation
       DO jk = kstart, nlev !start from kstart, as in art_aerosol_washout kstart = 15
!NEC$ ivdep
         DO jc = istart, iend
           p_art_data%diag%acc_wetdepo_con(jc,jb,idx_diag) =       &
             & p_art_data%diag%acc_wetdepo_con(jc,jb,idx_diag)     &
             & - wash_rate(jc,jk) * dtime / third_mom(jc,jk)       &
             & * tracer(jc,jk) * 1.e-9_wp * dz(jc,jk,jb)
         !Note: *rho not necessray anymore (see e.g. art_integrate_explicit)
         ENDDO !jc
       ENDDO !jk
     END IF

  ENDIF

  !------ ACC_WETDEPO_RRSFC:

  ! get index in diagnostics container (given the index in the tracer container)
  idx_diag = art_diag_tracer_index(IART_ACC_WETDEPO_RRSFC, idx_tracer)

  ! save wet deposition only if precipitation reached the surface
  IF ( idx_diag > 0 ) THEN
    DO jk = kstart,nlev !start from kstart, as in art_aerosol_washout kstart = 15
!NEC$ ivdep
      DO jc = istart, iend
        IF (rr_sfc(jc) > 0.0_wp) THEN
          p_art_data%diag%acc_wetdepo_rrsfc(jc,jb,idx_diag) =     &
            & p_art_data%diag%acc_wetdepo_rrsfc(jc,jb,idx_diag)   &
            & - wash_rate(jc,jk) * dtime / third_mom(jc,jk)       &
            & * tracer(jc,jk) * 1.e-9_wp * dz(jc,jk,jb)
        !Note: *rho not necessray anymore (see e.g. art_integrate_explicit)
        END IF !(rr_sfc(jc) > 0.0_wp)
      ENDDO !jc
    ENDDO !jk
  END IF

END SUBROUTINE art_save_aerosol_wet_deposition
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_save_aerosol_emission(this,p_art_data, emiss_rate_mass, emiss_rate_numb,   &
  &                                  dz, dtime, jb, istart, iend, kstart, nlev, opt_fac)

  !<
  ! SUBROUTINE save_art_aerosol_emission
  ! integrates the 3d emission rates for art-tracers
  ! vertically and saves the 2d emisison of art-tracers
  !
  ! Based on: -
  ! Part of Module: mo_art_emission_interface
  ! Author: Andrea Steiner, DWD
  ! Initial Release: 2018-01-02
  ! Modifications:
  ! 2022-09-09: rewritten by Anika Rohde, KIT
  ! -
  !>

  CLASS(t_art_emiss2tracer), INTENT(INOUT) :: &
    &  this
  TYPE(t_art_data),    INTENT(inout) :: &
    & p_art_data                          !< ART data container
  INTEGER,             INTENT(in)    :: &
    & jb,                               & !< Counter for block loop
    & istart, iend,                     & !< Start and end of art_atmo%nproma loop
    & kstart,                           & !< Index of level to start vertical loop
    & nlev                                !< Number of levels (equals index of lowest full level)
  REAL(wp),            INTENT(in)    :: &
    & emiss_rate_mass(istart:iend,1:nlev,1:this%nmodes), & !< Mass emission rate   [kg m-3 s-1]
    & emiss_rate_numb(istart:iend,1:nlev,1:this%nmodes), & !< Number emission rate [m-3 s-1]
    & dtime,                            & !< Time step (fast physics)
    & dz(:,:,:)                           !< Vertical layer thickness (m)
  REAL(wp),            OPTIONAL      :: &
    & opt_fac                             !< Optional scalar factor, to transform mass emission
                                          !  into number (c.f. art_integrate_explicit)

  !Local variables
  INTEGER :: idx_diagM, idx_diag_accM   !< Index of tracer in diagnostics container
  INTEGER :: idx_diag0, idx_diag_acc0   !< Index of tracer in diagnostics container
  INTEGER :: jc, jk, imod, itr          !< Loop indices
  REAL(wp):: factor, conv_fac           !< Local instance of opt_fac if present, conversion factor to SI units


  IF(PRESENT(opt_fac)) THEN
    factor = opt_fac
  ELSE
    factor = 1.0_wp
  ENDIF


  DO imod = 1, this%nmodes
    DO itr = 1, this%ntr
      idx_diagM     = art_diag_tracer_index(IART_EMISS,     this%itr3(imod,itr))
      idx_diag_accM = art_diag_tracer_index(IART_ACC_EMISS, this%itr3(imod,itr))

      ! compute emission of art-tracer
      IF ( idx_diagM > 0 ) THEN
        DO jc = istart, iend
          DO jk = kstart, nlev
            conv_fac  = dtime*this%weight(imod,itr) * 1.e-9_wp
            p_art_data%diag%emiss(jc,jb,idx_diagM) =                  &
              &   p_art_data%diag%emiss(jc,jb,idx_diagM)              &
              & + emiss_rate_mass(jc,jk,imod) * conv_fac * factor * dz(jc,jk,jb)
          ENDDO !jk
        ENDDO ! jc
      ENDIF ! idx_diagM

      ! compute accumulated emission of art-tracer
      IF ( idx_diag_accM  > 0 ) THEN
        DO jc = istart, iend
          DO jk = kstart, nlev
            conv_fac  = dtime*this%weight(imod,itr) * 1.e-9_wp
            p_art_data%diag%acc_emiss(jc,jb,idx_diag_accM) =          &
              &   p_art_data%diag%acc_emiss(jc,jb,idx_diag_accM)      &
              & + emiss_rate_mass(jc,jk,imod) * conv_fac * factor * dz(jc,jk,jb)
          ENDDO !jk
        ENDDO ! jc
      ENDIF !dx_diag_accM
    ENDDO  !itr

    idx_diag0 = art_diag_tracer_index(IART_EMISS, this%itr0(imod))
    IF ( idx_diag0 > 0 ) THEN
      DO jc = istart, iend
        DO jk = kstart, nlev
          p_art_data%diag%emiss(jc,jb,idx_diag0) =                    &
            &   p_art_data%diag%emiss(jc,jb,idx_diag0)                &
            & + emiss_rate_numb(jc,jk,imod) * factor * dtime * dz(jc,jk,jb)
        ENDDO !jk
      ENDDO !jc
    ENDIF ! idx_diag0

    idx_diag_acc0 = art_diag_tracer_index(IART_ACC_EMISS, this%itr0(imod))
    IF ( idx_diag_acc0 > 0 ) THEN
      DO jc = istart, iend
        DO jk = kstart, nlev
          p_art_data%diag%acc_emiss(jc,jb,idx_diag_acc0) =            &
            &   p_art_data%diag%acc_emiss(jc,jb,idx_diag_acc0)        &
            & + emiss_rate_numb(jc,jk,imod) * factor * dtime * dz(jc,jk,jb)
        ENDDO !jk
      ENDDO !jc
    ENDIF !idx_diag_acc0
  ENDDO !imod

END SUBROUTINE art_save_aerosol_emission
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_tropopause(jg, nccbot_in, ncctop_in, lacc)

    IMPLICIT NONE 

    INTEGER, INTENT(in) :: jg
    INTEGER, INTENT(in), optional :: nccbot_in
    INTEGER, INTENT(in), optional :: ncctop_in   !< min. and max. level to search for tropopause
    LOGICAL, INTENT(in), optional :: lacc

    INTEGER :: jb, jc, jk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: nccbot, ncctop
    LOGICAL :: lzacc
 
    !Variables for tropopause and TTL determination
    LOGICAL    :: lo1, lo2

    TYPE(t_art_atmo),POINTER    :: &
       &  art_atmo                     !< Pointer to ART atmonostic fields

    CALL set_acc_host_or_device(lzacc, lacc)

    art_atmo   => p_art_data(jg)%atmo

    IF (PRESENT(nccbot_in)) THEN
      nccbot = nccbot_in
    ELSE
      nccbot = aes_wmo_config(jg)%jkewmo
    END IF

    IF (PRESENT(ncctop_in)) THEN
      ncctop = ncctop_in
    ELSE
      ncctop = aes_wmo_config(jg)%jkswmo+2
    END IF


    DO jb = art_atmo%i_startblk, art_atmo%i_endblk

      CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

      IF (PRESENT(nccbot_in)) THEN
        CALL WMO_tropopause(jg, 1, i_endidx, art_atmo%nproma, art_atmo%nlev,      &
            &              art_atmo%temp(:,:,jb), art_atmo%pres(:,:,jb),          &
            &              art_atmo%ptropo(:,jb), nccbot, ncctop, lacc=lzacc)
      ELSE
        CALL WMO_tropopause(jg, 1, i_endidx, art_atmo%nproma, art_atmo%nlev,      &
            &              art_atmo%temp(:,:,jb), art_atmo%pres(:,:,jb),          &
            &              art_atmo%ptropo(:,jb), lacc=lzacc)
      END IF
    
!NEC$ ivdep
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx,i_endidx
        art_atmo%ktrpwmop1(jc,jb) = ncctop
      END DO
    
      !$ACC LOOP SEQ
      DO jk = ncctop-1, nccbot+1
!NEC$ ivdep
        !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(lo1, lo2)
        DO jc = i_startidx, i_endidx
    
          lo1 =           art_atmo%pres(jc,jk,jb)       < art_atmo%ptropo(jc,jb)
          lo2 = lo1 .AND. art_atmo%pres(jc,jk+1,jb)     > art_atmo%ptropo(jc,jb)
    
          art_atmo%ktrpwmop1(jc,jb)  = MERGE(jk+1, art_atmo%ktrpwmop1(jc,jb), lo2)
          art_atmo%ktrpwmo(jc,jb)    = art_atmo%ktrpwmop1(jc,jb) - 1
    
          art_atmo%ktrpwmop1_real(jc,jb)  =  art_atmo%ktrpwmop1(jc,jb)
    
        END DO
      END DO
      !$ACC END PARALLEL

    END DO

END SUBROUTINE art_calc_tropopause

SUBROUTINE art_calc_o2_column(jg,jb,jcs,jce,nproma,nlev,o2_column)

  IMPLICIT NONE 

  INTEGER,  INTENT(in)    :: jg,jb,jcs,jce
  INTEGER,  INTENT(in)    :: nproma,nlev
  REAL(wp), INTENT(inout) :: o2_column(:,:,:)

  INTEGER :: jc, jk

  REAL(wp)              ::  press_alt_u              !< pressure altitude in m
  REAL(wp)              ::  press_alt_l              !< pressure altitude in m

  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo           !< Pointer to ART atmospheric fields 

  art_atmo   => p_art_data(jg)%atmo
  
!NEC$ ivdep
  DO jc = jcs, jce
    press_alt_u = 100000._wp
    press_alt_l =  -7000._wp * LOG(art_atmo%pres(jc,art_atmo%nlev,jb)/p0sl_bg)
    o2_column(jc,art_atmo%nlev,jb) = art_atmo%pres(jc,art_atmo%nlev,jb)    &
      &             / argas/art_atmo%temp(jc,art_atmo%nlev,jb)      &
      &             * avo/1000000._wp                              &
      &             * (press_alt_u-press_alt_l)*100._wp*0.2041_wp
  END DO

  DO jk = art_atmo%nlev-1,1,-1
!NEC$ ivdep
    DO jc = jcs, jce
      press_alt_u = -7000._wp * LOG(art_atmo%pres(jc,jk+1,jb)/p0sl_bg)
      press_alt_l = -7000._wp * LOG(art_atmo%pres(jc,jk  ,jb)/p0sl_bg)
      o2_column(jc,jk,jb) = art_atmo%pres(jc,jk,jb)                    &
        &              / argas/art_atmo%temp(jc,jk,jb)          &
        &              * avo/1000000._wp                       &
        &              * (press_alt_u-press_alt_l)*100._wp*0.2041_wp   &
        &              + o2_column(jc,jk+1,jb)
    ENDDO
  ENDDO

END SUBROUTINE art_calc_o2_column
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_tau_vi( istart, iend, nlev, jb, aeronet )
!<
! SUBROUTINE art_calc_tau_vi
! integrates the 3d AOD for art-tracers vertically
! and saves the 2d AOD of art-tracers
!
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2023-08-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! -
!>
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_aeronet), INTENT(inout) :: &
    &  aeronet(:)
! Local variables
  INTEGER     ::            &
    &  jc, jk, i_wavel        !< Loop indices (art_atmo%nproma, nlev, wavelength)
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_calc_tau_vi"

  ! ----------------------------------------------------------------------
  ! --- Calculate total column of aerosol optical depth
  ! ----------------------------------------------------------------------
  DO i_wavel = 1,9
     aeronet(i_wavel)%tau_vi(:,jb) = 0.0_wp
     DO jk = 1, nlev
!NEC$ ivdep
        DO jc = istart, iend
           aeronet(i_wavel)%tau_vi(jc,jb) =       &
             &    aeronet(i_wavel)%tau_vi(jc,jb)  &
             &  + aeronet(i_wavel)%tau(jc,jk,jb)
        END DO !jc
     END DO !jk
  END DO !i_wavel

END SUBROUTINE art_calc_tau_vi
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_total_mc( rho, p_trac, istart, iend, nlev, jb, &
  &                           p_art_data, emiss_name, total_mc)
!<
! SUBROUTINE art_calc_total_mc
! it calculates the total (sum over all modes) mass concentration [kg/m3]
!
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2023-08-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! -
!>
  REAL(wp), INTENT(in)   :: &
    &  rho(:,:),            & !< Air density  [kg/m3]
    &  p_trac(:,:,:)          !< Tracer mixing ratios [mug/kg]
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data), INTENT(inout), TARGET :: &
    &  p_art_data             !< Data container for ART
  CHARACTER(LEN=*), INTENT(in) :: &
    &  emiss_name
  REAL(wp), INTENT(out)   :: &
    &  total_mc(:,:,:)
! Local variables
  INTEGER     ::            &
    &  jc, jk,              & !< Loop indices (art_atmo%nproma, nlev)
    &  imod, itr              !< Loop indices (modes, tracers)
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_calc_total_mc"
  TYPE(t_mode), POINTER  :: &
    & current

  ! -------------------------------------------------------------
  ! --- Calculate total mass concentration [kg/m3]
  ! -------------------------------------------------------------
  total_mc(:,:,jb) = 0.0_wp
  IF (p_art_data%tracer2aeroemiss%lisinit) THEN
    current=>p_art_data%tracer2aeroemiss%e2t_list%p%first_mode
    DO WHILE(ASSOCIATED(current))
      SELECT TYPE(this=>current%fields)
      TYPE IS(t_art_emiss2tracer)
        IF (this%name(1:4) == emiss_name(1:4)) THEN
          IF (this%lcalcemiss) THEN
            DO imod = 1, this%nmodes
              DO itr = 1, this%ntr
                DO jk = 1, nlev
!NEC$ ivdep
                  DO jc = istart, iend
                    total_mc(jc,jk,jb) = total_mc(jc,jk,jb)             &
                      &  + p_trac(jc,jk,this%itr3(imod,itr))*rho(jc,jk) &
                      &    * 1.e-9_wp
                  ENDDO !jc
                ENDDO !jk
              ENDDO ! itr
            ENDDO ! imod
          ENDIF ! this%lcalcemiss
        ENDIF
      END SELECT
      current=>current%next_mode
    END DO
  END IF

  ! clipping of negative mass conc. values
  CALL art_clip_lt(total_mc(:,:,jb),0.0_wp)
END SUBROUTINE art_calc_total_mc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_total_mc_vi( dz, istart, iend, nlev, jb, &
  &                              total_mc, total_mc_vi )
!<
! SUBROUTINE art_calc_total_mc_vi
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2023-08-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)   :: &
    &  dz(:,:)                !< Layer thickness
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  REAL(wp), INTENT(in)   :: &
    &  total_mc(:,:,:)        !< total mass concentration (sum over all modes)
  REAL(wp), INTENT(out)  :: &
    &  total_mc_vi(:,:)       !< total column mass concentration
! Local variables
  INTEGER     ::            &
    &  jc, jk                 !< Loop indices (nproma, nlev)
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_calc_total_mc_vi"

  ! ----------------------------------------------------------------------
  ! --- Calculate total column mass concentration [kg/m2]
  ! ----------------------------------------------------------------------
  total_mc_vi(:,jb) = 0.0_wp
  DO jk = 1, nlev
!NEC$ ivdep
    DO jc = istart, iend
      total_mc_vi(jc,jb) = total_mc_vi(jc,jb)  &
        &                + dz(jc,jk) * total_mc(jc,jk,jb)
    END DO !jc
  END DO !jk
  
END SUBROUTINE art_calc_total_mc_vi
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_diag_max_preslay( pres, istart, iend, nlev, jb,        &
  &                              p_art_data, concentr, max_preslay )
!<
! SUBROUTINE art_diag_max_preslay
! Based on: -
! Part of Module: mo_art_diagnostics
! Author: Jochen Foerstner, DWD
! Initial Release: 2023-08-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  REAL(wp), INTENT(in)   :: &
    &  pres(:,:)              !< Air pressure [Pa]
  INTEGER, INTENT(in)    :: &
    &  istart, iend,        & !< Start and end index of art_atmo%nproma loop
    &  nlev, jb               !< Number of verical levels, Block index
  TYPE(t_art_data), INTENT(inout), TARGET :: &
    &  p_art_data             !< Data container for ART
  REAL(wp), INTENT(in)   :: &
    &  concentr(:,:,:)        !< mass or air concentration
  TYPE(t_art_max_layer), INTENT(inout) :: &
    &  max_preslay(:)
! Local variables
  REAL(wp)               :: &
    &  pres_bot, pres_top
  INTEGER     ::            &
    &  jc, jk, jp             !< Loop indices (nproma, nlev, npreslay)
  REAL(wp)               :: &
    &  masked_field(istart:iend,1:nlev)
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diagnostics:art_diag_max_preslay"

  ! -----------------------------------------------------------------------------------------------
  ! -- Diagnose maximum mass or air concentration between given pressure levels [kg/m3] / [Bq/m3]
  ! -----------------------------------------------------------------------------------------------
  DO jp = 1, npreslay
    pres_bot = max_preslay(jp)%pres_bot
    pres_top = max_preslay(jp)%pres_top
    DO jk = 1, nlev
!NEC$ ivdep
      DO jc = istart, iend
        masked_field(jc,jk) = concentr(jc,jk,jb)
        IF (pres_bot > 0.0_wp .AND. pres(jc,jk) >  pres_bot) masked_field(jc,jk) = 0.0_wp
        IF (pres_top > 0.0_wp .AND. pres(jc,jk) <= pres_top) masked_field(jc,jk) = 0.0_wp
      END DO !jc
    END DO !jk
!NEC$ ivdep
    DO jc = istart, iend
      max_preslay(jp)%maximum(jc,jb) = MAXVAL( masked_field(jc,:) )
    END DO !jc
  END DO !jp

END SUBROUTINE art_diag_max_preslay
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_diagnostics
