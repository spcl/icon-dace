! Provide an implementation of the zero-layer sea ice thermodynamics (Semtner, 1976).
! Used by both the atmopshere and ocean models.
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
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ice_zerolayer

  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_model_domain,        ONLY: t_patch
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, ki, ks, ci, alf, zemiss_def, stbo, tmelt
  USE mo_sea_ice_nml,         ONLY: hci_layer, use_no_flux_gradients
  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ice_growth_zerolayer
  PUBLIC :: set_ice_temp_zerolayer
  PUBLIC :: set_ice_temp_zerolayer_analytical

  CHARACTER(len=12)           :: str_module    = 'IceZeroLayer'  ! Output of module for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------------
  !>
  !! ! ice_growth_zerolayer:: change ice and snow thickness (Semtner 1976, Appendix)
  !!
  !! This function changes:
  !! ice % hs       new snow thickness for each ice category                [m]
  !! ice % hi       new ice  thickness for each ice category                [m]
  !! ice % Qbot     heat flux available for freezing/melting at ice bottom  [W/m^2]
  !! ice % heatOceI heat flux to the surface ocean layer (mixed layer)      [W/m^2]
  !! --- not currently ---
  !! ice % evapwi   amount of evaporated water from the mixed layer
  !!                in previously ice covered areas if all ice is gone      [kg/m^3]
  !!
  !! The counterpart to the ice_growth subroutine in MPIOM.
  !!
  SUBROUTINE ice_growth_zerolayer(p_patch, ice, lacc)
    TYPE(t_patch),             INTENT(IN), TARGET    :: p_patch
    TYPE(t_sea_ice),           INTENT(INOUT)         :: ice
    LOGICAL, INTENT(IN), OPTIONAL                    :: lacc

    ! Local variables
    REAL(wp), DIMENSION(nproma, ice%kice, p_patch%alloc_cell_blocks) ::        &
      & Q_surplus,  & ! Energy surplus during ice growth
      & hiold,      & ! Ice thickness before thermodynamic effects
      & hsold         ! Snow thickness before thermodynamic effects

    ! Loop indices
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: k, jb, jc, i_startidx_c, i_endidx_c
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(Q_surplus, hiold, hsold) IF(lzacc)

    !-------------------------------------------------------------------------------------------
    all_cells            => p_patch%cells%all


    !---------DEBUG DIAGNOSTICS-----------------------------------------------------------------
    CALL dbg_print('GrowZero bef.: Qtop'     , ice%Qtop     , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero bef.: Qbot'     , ice%Qbot     , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero bef.: Qbot_slow', ice%Qbot_slow, str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero bef.: hi'       , ice%hi       , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero bef.: hs'       , ice%hs       , str_module, 4, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------------------

!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, k, jc) SCHEDULE(dynamic)
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) IF(lzacc)
      DO k=1,ice%kice
        DO jc = i_startidx_c,i_endidx_c
        ! initialization
          Q_surplus   (jc, k, jb)    = 0.0_wp
          ! Save ice and snow thickness before thermodynamic effects
          hiold (jc, k, jb) = ice%hi(jc, k, jb)
          hsold (jc, k, jb) = ice%hs(jc, k, jb)
    !      ice%heatOceI(:,:,:) = 0.0_wp ! initialized in ice_zero
          IF (ice%hi(jc,k,jb) > 0._wp) THEN
            !     ------------------------------------------------
            ! (1) --------------- Update hi and hs ---------------
            !     ------------------------------------------------

            ! Add snowfall ice%totalsnowfall (in meters)
            ! Note, that liquid liquid precipitation (rprecw) is put into the ocean in ice_open_ocean
            ice%hs(jc,k,jb) = ice%hs(jc,k,jb) + ice%totalsnowfall(jc,jb)*rho_ref/rhos

            ! Qtop>0, if there is atmospheric heat flux available for melting
            IF ( ice%Qtop(jc,k,jb) > 0.0_wp ) THEN
              IF  ( ice%hs(jc,k,jb) > 0.0_wp )  THEN ! first melt any snow that is present

                ice%hs (jc,k,jb) =  ice%hs(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhos)
                ! put remaining heat, into melting ice below
                IF (ice%hs(jc,k,jb) < 0.0_wp) THEN
                  ice%hi(jc,k,jb) = ice%hi(jc,k,jb) + ice%hs(jc,k,jb) * (rhos/rhoi) ! snow thickness loss in ice equivalents
                  ice%hs(jc,k,jb) = 0.0_wp
                ENDIF

              ELSE ! where there's no snow

                ice%hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qtop(jc,k,jb) * dtime / (alf*rhoi)

              ENDIF
            ENDIF

            ! Heat flux at ice-ocean interface ("Qbot_slow") is sum of
            !  - Conductive heat flux calculated in ice_fast ("Qbot"):
            !    Qbot<0 is upward (Qbot=-F_S, negative F_S is upward) -> heat release to the atmosphere, ice growth
            !  - ocean-ice heat flux ("zHeatOceI", positive upward is melting ice) to Qbot
            ice%Qbot_slow(jc,k,jb) = ice%Qbot(jc,k,jb) + ice%zHeatOceI(jc,k,jb)

            !  Resulting flux Qbot_slow>0 melts ice from below, while Qbot_slow<0 grows ice.
            ice%hi(jc,k,jb) = ice%hi(jc,k,jb) - ice%Qbot_slow(jc,k,jb) * dtime / (alf*rhoi)

            !     ------------------------------------------------
            ! (2) --------------- Update heatOceI ----------------
            !     ------------------------------------------------

            !  hi < 0 - all ice melted, calculate residual heat flux into ocean mixed layer
            IF (ice%hi (jc,k,jb) < 0.0_wp) THEN

              ! heat flux is available for melting snow on top of melted ice and/or heating the ocean
              ! Q_surplus checks whether is enough to melt all the snow which was on top of melted ice
              Q_surplus(jc,k,jb) = ( -ice%hi(jc,k,jb)*rhoi - ice%hs(jc,k,jb)*rhos )*alf/dtime
              ice%hi   (jc,k,jb) = 0.0_wp

              ! melt snow
              IF ( Q_surplus(jc,k,jb) < 0 ) THEN ! not enough
                 ice%hs(jc,k,jb) = ice%hs(jc,k,jb) - (-ice%hi(jc,k,jb)*rhoi)/rhos ! guaranteed to be positive
                 Q_surplus(jc,k,jb) = 0
                 ! note: remaining snow will be converted into ice in the flooding routine.
              ELSE
                 ice%hs   (jc,k,jb) = 0.0_wp            ! all snow was melted
                 ice%Tsurf(jc,k,jb) = ice%Tfw(jc,jb)
              ENDIF
            ENDIF

            ! heatOceI - positive into ocean (positive=downward), i.e. same sign convention as HeatFlux_Total
            ! put the surplus after complete melting into heating of the ocean
            ice%heatOceI(jc,k,jb) = - ice%zHeatOceI(jc,k,jb) + Q_surplus(jc,k,jb)

            ! rescale with concentration, since later we calculate: HeatFlux_Total = heatOceI + heatOceW
            ice%heatOceI(jc,k,jb) = ice%heatOceI(jc,k,jb)*ice%conc(jc,k,jb)

            ! rescale with concentration in order to get output over whole cell-area
            ice%zHeatOceI(jc,k,jb) = ice%zHeatOceI(jc,k,jb)*ice%conc(jc,k,jb)
            ice%Qbot_slow(jc,k,jb) = ice%Qbot_slow(jc,k,jb)*ice%conc(jc,k,jb)  !! Attention !!
            Q_surplus    (jc,k,jb) =     Q_surplus(jc,k,jb)*ice%conc(jc,k,jb)

            ! Calculate mean change in ice and snow thickness due to thermodynamic effects
            ice%delhi(jc,k,jb) = ( ice%hi(jc,k,jb) - hiold(jc,k,jb) ) * ice%conc(jc,k,jb)
            ice%delhs(jc,k,jb) = ( ice%hs(jc,k,jb) - hsold(jc,k,jb) ) * ice%conc(jc,k,jb)

          ELSE  ! hi<=0
            ice%Qbot_slow(jc,k,jb) = 0.0_wp
            ice%Qbot     (jc,k,jb) = 0.0_wp
            ice%zHeatOceI(jc,k,jb) = 0.0_wp
            ice%heatOceI (jc,k,jb) = 0.0_wp
          ENDIF
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
!ICON_OMP_END_PARALLEL_DO

    !$ACC END DATA

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('GrowZero aft.: hi'         , ice%hi         , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: hs'         , ice%hs         , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: Qbot_slow'  , ice%Qbot_slow  , str_module, 3, in_subset=p_patch%cells%owned)

    CALL dbg_print('GrowZero aft.: zHeatOceI ' , ice%zHeatOceI  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: heatOceI '  , ice%heatOceI   , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: Qbot'       , ice%Qbot       , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: Q_surplus'  , Q_surplus      , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: delhi'      , ice%delhi      , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: delhs '     , ice%delhs      , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: conc  '     , ice%conc       , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('GrowZero aft.: concSum'    , ice%concSum    , str_module, 4, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

!---------------------------------------------------------------------

  END SUBROUTINE ice_growth_zerolayer

  !-------------------------------------------------------------------------------
  !>
  !! ! set_ice_temp_zerolayer:: calculate new ice and snow temperatures (Semtner 1976, Appendix)
  !!
  !! This function changes:
  !! ice % Tsurf    the new surface temperature for each ice category       [C]
  !! ice % Qbot     heat flux available for freezing/melting at ice bottom  [W/m^2]
  !! ice % Qtop     heat flux available for melting at ice surface          [W/m^2]
  !!
  !! Routines set_ice_temp_zero_nogradients is now a part of set_ice_temp_zerolayer.
  !! Switching between nogradients and full version is controlled by nfg_flag.
  !
  SUBROUTINE set_ice_temp_zerolayer(i_startidx_c, i_endidx_c, nbdim, kice, pdtime, &
            &   Tsurf,          &  ! surface temperature
            &   hi,             &  ! ice thickness
            &   hs,             &  ! snow thickness
            &   Qtop,           &  ! heat flux available for melting at ice surface
            &   Qbot,           &  ! heat flux available for freezing/melting at ice bottom
            &   SWnet,          &  ! net shortwave flux (includes albedo effect)
            &   nonsolar,       &  ! net nonsolar fluxes = lat + sens + LWnet
            &   dnonsolardT,    &  ! gradient of nonsolar fluxes = dlatdT + dsensdT + dLWdT
            &   Tfw,            &  ! sea surface freezing temperature
            &   lacc)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(IN)    :: pdtime
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: SWnet      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: nonsolar   (nbdim,kice)
    REAL(wp),INTENT(IN)    :: dnonsolardT(nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! Local variables
    REAL(wp) ::             &
      & F_A         ,       &  ! net atmospheric heat flux                  (positive=upward)
      & F_S         ,       &  ! conductive heat flux through snow and ice  (positive=upward)
      & k_effective ,       &  ! effective heat conductivity of ice/snow
      & deltaTdenominator,  &  ! prefactor of deltaT in sfc. flux balance formula
      & deltaT      ,       &  ! temperature increment from the prev timestep
      & c_icelayer  ,       &  ! contribution from heat capacity of a thin ice layer
      & nfg_flag               ! switches on/off the contribution of dnonsolardT

    ! & I_0 -- fraction of the net SW radiation penetrating the ice (not in use)

    INTEGER :: k, jk, jc ! loop indices
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! initialization of the output
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO k = 1,kice
      !$ACC LOOP GANG VECTOR
      DO jk = 1,nbdim
        Qbot(jk,k) = 0._wp
        Qtop(jk,k) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL

    ! declaration of constants
    ! ToDo: move into the sea-ice initialization step, doesn't belong here
    c_icelayer = rhoi*hci_layer*ci/pdtime

    IF (use_no_flux_gradients) THEN
        nfg_flag = 0._wp
    ELSE
        nfg_flag = 1._wp
    ENDIF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=1,kice
      !$ACC LOOP GANG VECTOR PRIVATE(k_effective, F_A, F_S, deltaT) &
      !$ACC   PRIVATE(deltaTdenominator)
      DO jc = i_startidx_c,i_endidx_c
        IF (hi(jc,k) > 0._wp) THEN

          k_effective = ki*ks/(ks*hi(jc,k) + ki*hs(jc,k))   ! conductivity of the ice-snow system
          F_S = - k_effective * (Tsurf(jc,k) - Tfw(jc))     ! conductive flux
          F_A = -SWnet(jc,k) - nonsolar(jc,k)               ! net atmospheric flux
          ! c_icelayer is added to stabilize the atmosphere
          deltaTdenominator = k_effective  - nfg_flag*dnonsolardT(jc,k) + c_icelayer

          deltaT = (F_S - F_A) / deltaTdenominator          ! temperature increment

          ! when melting, Tsurf = 0 deg C
          IF (Tsurf(jc,k) + deltaT > 0.0_wp)  deltaT = -Tsurf(jc,k)
          ! now assign new Tsurf
          Tsurf(jc,k) = Tsurf(jc,k) + deltaT

          ! Heat flux available for surface melting Qtop = -(F_A - F_S) evaluated at the new Tsurf.
          ! 1) if new Tsurf < 0, then surface fluxes are balanced, and Qtop == 0.
          ! 2) if new Tsurf = 0, then Qtop > 0 (goes into surface melting in ice_growth_zerolayer).
          Qtop(jc,k) = - F_A + F_S - deltaTdenominator * deltaT

          ! Heat flux available for bottom melting/freezing Qbot = -F_S evaluated at the new Tsurf.
          ! Note that Qbot is missing the flux from ocean to ice (added in ice_growth_zerolayer).
          Qbot(jc,k) = - F_S + k_effective * deltaT

        ELSE ! hi(jc,k) <= 0

            Tsurf(jc,k) = Tfw(jc) ! is it necessary?

        END IF
      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE set_ice_temp_zerolayer
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !>
  !! ! set_ice_temp_zerolayer_analytical: a copy of set_ice_temp_zerolayer above
  !!    with analytical atmospheric fluxes. Retained for historical reasons.
  !!
  !! Initial release by Vladimir Lapin, MPI (2016-11)
  !
  SUBROUTINE set_ice_temp_zerolayer_analytical(i_startidx_c, i_endidx_c, nbdim, kice, &
            &   Tsurf, hi, hs, Qtop, Qbot, Tfw, doy, lacc)

    INTEGER, INTENT(IN)    :: i_startidx_c, i_endidx_c, nbdim, kice
    REAL(wp),INTENT(INOUT) :: Tsurf      (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hi         (nbdim,kice)
    REAL(wp),INTENT(IN)    :: hs         (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qtop       (nbdim,kice)
    REAL(wp),INTENT(OUT)   :: Qbot       (nbdim,kice)
    REAL(wp),INTENT(IN)    :: Tfw        (nbdim)
    INTEGER, INTENT(IN)    :: doy
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    ! Local variables
    REAL(wp) ::             &
      & F_A         ,       &  ! net atmospheric heat flux                  (positive=upward)
      & F_S         ,       &  ! conductive heat flux through snow and ice  (positive=upward)
      & k_effective ,       &  ! effective heat conductivity of ice/snow
      & deltaTdenominator,  &  ! prefactor of deltaT in sfc. flux balance formula
      & deltaT                 ! temperature increment from the prev timestep
    ! Loop indices
    INTEGER :: k, jk, jc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! initialization of output variables
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO k = 1,kice
      !$ACC LOOP GANG VECTOR
      DO jk = 1,nbdim
        Qbot(jk,k) = 0._wp
        Qtop(jk,k) = 0._wp
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=1,kice
      !$ACC LOOP GANG VECTOR PRIVATE(k_effective, F_A, F_S, deltaT) &
      !$ACC   PRIVATE(deltaTdenominator)
      DO jc = i_startidx_c,i_endidx_c
        IF (hi(jc,k) > 0._wp) THEN

          k_effective = ki*ks/(ks*hi(jc,k) + ki*hs(jc,k))
          F_S = k_effective * (Tfw(jc) - Tsurf(jc,k))

          ! #achim: hard-coding simpler form of atmospheric fluxes (from Dirk's thesis, p.193)
          F_A =   zemiss_def * StBo * (Tsurf(jc,k) +tmelt)**4   &
              &            - ( 118.0_wp * EXP(-0.5_wp*((doy-206)/53.1_wp)**2) + 179.0_wp ) &
              &            - 314.0_wp * EXP(-0.5_wp*((doy-164)/47.9_wp)**2)  & !SW, NO 1-I_0 factor
              &                * ( 0.431_wp / (1.0_wp+((doy-207)/44.5_wp)**2) - 0.086_wp)!1-albedo

          deltaTdenominator = k_effective + 4.0_wp*zemiss_def*StBo*(Tsurf(jc,k)+tmelt)**3
          deltaT = (F_S - F_A) / deltaTdenominator

          IF (Tsurf(jc,k) + deltaT > 0.0_wp) THEN
            deltaT = -Tsurf(jc,k)
          END IF
          Tsurf(jc,k) = Tsurf(jc,k) + deltaT

          Qtop(jc,k) = -F_A + F_S - deltaTdenominator * deltaT
          Qbot(jc,k) = - F_S + k_effective * deltaT

        ELSE ! hi(jc,k) <= 0

            Tsurf(jc,k) = Tfw(jc)

        END IF
      END DO
    END DO
    !$ACC END PARALLEL

  END SUBROUTINE set_ice_temp_zerolayer_analytical
  !-------------------------------------------------------------------------------

END MODULE mo_ice_zerolayer
