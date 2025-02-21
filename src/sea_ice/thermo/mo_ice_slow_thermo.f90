! Provide an implementation of ice thermodynamic growth model: ice_slow_thermo.
! Calculates ice growth rate over thick ice and open ocean. Updates heat
! and fresh-water fluxes which are passed to the ocean for coupling.
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

!----------------------------
MODULE mo_ice_slow_thermo
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2007
  !
  !-------------------------------------------------------------------------
  !
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime, ltimer
  USE mo_dynamics_config,     ONLY: nold
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_impl_constants,      ONLY: sea_boundary
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, alf, clw
  USE mo_sea_ice_nml,         ONLY: i_ice_therm, hmin, hnull, leadclose_1, leadclose_2n, sice
  USE mo_ocean_nml,           ONLY: limit_seaice, limit_seaice_type, seaice_limit, &
    &                               vert_cor_type,                                 &
    &                               seaice_limit_abs
  USE mo_ocean_state,         ONLY: v_base
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_ocean_surface_types, ONLY: t_ocean_surface
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_timer,               ONLY: timer_start, timer_stop, timer_ice_slow

  USE mo_ice_parameterizations, ONLY: oce_ice_heatflx, ice_draft_and_flooding
  USE mo_ice_diagnostics,     ONLY: energy_in_surface
  USE mo_sea_ice,             ONLY: ice_conc_change
  USE mo_ice_zerolayer,       ONLY: ice_growth_zerolayer
  USE mo_ice_winton,          ONLY: ice_growth_winton
  USE mo_exception,           ONLY: finish
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC  :: ice_slow_thermo

  PRIVATE :: ice_growth_init
  PRIVATE :: ice_open_ocean
  PRIVATE :: update_ice_ocean_fluxes
  PRIVATE :: ice_update_conc
  PRIVATE :: ice_thickness_limiter

  CHARACTER(len=12)           :: str_module    = 'IceSlowTherm'  ! Output of module for 1 line debug

CONTAINS

  !--------------------------------------------------------------------------------
  !>
  !! ! ice_slow_thermo: Calculates ice growth rates and updates ice structure accordingly.
  !!             Updates heat and fresh-water fluxes which are passed to the ocean.
  !!
  !! This function changes:
  !! ice            thermodynamics-related fields in the sea ice data structure
  !! p_oce_sfc      heat and fresh-water fluxes passed to the ocean
  !!
  SUBROUTINE ice_slow_thermo(p_patch_3D, p_os, atmos_fluxes, ice, p_oce_sfc, lacc)
    TYPE(t_patch_3D), TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state),INTENT(IN)    :: p_os         ! sst, sss, ssh only
    TYPE(t_atmos_fluxes),     INTENT(IN)    :: atmos_fluxes
    TYPE(t_sea_ice),          INTENT(INOUT) :: ice
    TYPE(t_ocean_surface),    INTENT(INOUT) :: p_oce_sfc
    LOGICAL, INTENT(IN), OPTIONAL           :: lacc

    ! Local variables
    TYPE(t_patch),  POINTER :: p_patch
    LOGICAL :: lzacc

    REAL(wp), POINTER :: sss    (:,:) ! sea surface salinity (input only)       [psu]
    REAL(wp), POINTER :: sst    (:,:) ! sea surface temperature (input only)    [C]
    REAL(wp), POINTER :: ssh    (:,:) ! sea surface height (input only)         [m]

    REAL(wp), DIMENSION (nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks) :: energyCheck

    ! energy_in_surface below only returns something useful if energyCheck_dbg_lev
    ! is greater than 3. We therefore set the debug level for the following
    ! dbg_print calls to 4. To get output set energyCheck_dbg_lev to 4 and set the
    ! debug level in the namelist accordingly.

    INTEGER, PARAMETER :: energyCheck_dbg_lev   = 3
    INTEGER, PARAMETER :: energyCheck_dbg_print = 4
    INTEGER            :: jc, jb, i_startidx_c, i_endidx_c

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch => p_patch_3D%p_patch_2D(1)

    sst     => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    sss     => p_os%p_prog(nold(1))%tracer(:,1,:,2)
    ssh     => p_os%p_prog(nold(1))%h(:,:)
    !-----------------------------------------------------------------------

    IF (ltimer) CALL timer_start(timer_ice_slow)

    ! initialize growth-related variables with zeros (before ice_growth_*)
    CALL ice_growth_init (ice, p_patch, lacc=lzacc)

    ! heat flux from ocean into ice: ice%zHeatOceI. Will be applied in ice_growth_* routine
    CALL oce_ice_heatflx (p_patch, p_os, ice, lacc=lzacc)

    ! totalsnowfall is applied in ice_growth_*
    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        ice%totalsnowfall(jc, jb) =  atmos_fluxes%rpreci(jc, jb) * dtime
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    ! thick ice growth/melt (K-classes): calculates ice%hs, ice%hi, ice%heatOceI
    !-------------------------------------------------------------------------------
    CALL dbg_print('IceSlow: befZero: totalSnF', ice%totalsnowfall, str_module, 4, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------
    SELECT CASE (i_ice_therm)
    CASE (1,3)
        CALL ice_growth_zerolayer (p_patch, ice, lacc=lzacc)
    CASE (2)
        CALL ice_growth_winton    (p_patch, ice, lacc=lzacc)
    END SELECT

    ! for historical reasons, ice%totalsnowfall represents cell-average, when applied in mo_ocean_surface*
    ! ToDo: should not be done this way

    DO jb = 1,p_patch%nblks_c
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        ice%totalsnowfall(jc, jb) =  ice%totalsnowfall(jc, jb) * ice%concSum(jc, jb)
      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !-------------------------------------------------------------------------------
    CALL dbg_print('IceSlow: aftZero: totalSnF', ice%totalsnowfall, str_module, 3, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------

    ! Flooding (snow to ice conversion)
    ! Calculates ice%snow_to_ice and updates ice%draft, ice%draftave
    ! Note that such conversion is assumed to occur with no heat or salt exchange with the ocean.
    ! In particular, the routine does not change ice%heatOceI
    CALL ice_draft_and_flooding (p_patch, ice, lacc=lzacc)

    ! fluxes and ice growth over open water
    CALL ice_open_ocean(p_patch, ice, atmos_fluxes, sst, lacc=lzacc)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    energyCheck = energy_in_surface(p_patch, ice, ssh(:,:), sst(:,:), computation_type=0, &
         &                          dbg_lev=energyCheck_dbg_lev)

    CALL dbg_print('IceSlow: energy aft. Growth', energyCheck, str_module, &
         &          energyCheck_dbg_print, in_subset=p_patch%cells%all)

    !---------------------------------------------------------------------

    ! updates fluxes that ocean will receive
    CALL update_ice_ocean_fluxes(p_patch, ice, atmos_fluxes, sss, p_oce_sfc, lacc=lzacc)

    ! Ice Concentration Change
    IF ( i_ice_therm >= 1 ) THEN
!        CALL ice_update_conc(p_patch,ice)
        CALL ice_conc_change(p_patch, ice, p_os, lacc=lzacc)
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    energyCheck = energy_in_surface(p_patch, ice, ssh(:,:), sst(:,:), computation_type=0, &
         &                          dbg_lev=energyCheck_dbg_lev)

    CALL dbg_print('IceSlow: energy aft. ConcChange', energyCheck, str_module, &
         &         energyCheck_dbg_print, in_subset= p_patch%cells%all)

    !---------------------------------------------------------------------

    ! limits ice thinkness, adjust p_oce_sfc fluxes and calculates the final freeboard
    ! ToDo: limit ice thickness in the ice_growth_*. Then correction to p_oce_sfc will not be needed
    IF (limit_seaice) THEN
        IF (limit_seaice_type .eq. 1) CALL ice_thickness_limiter( p_patch, ice, p_oce_sfc, p_os, lacc=lzacc )
        IF (limit_seaice_type .eq. 2) CALL ice_thickness_limiter_hh( p_patch_3d, ice, p_oce_sfc, p_os, lacc=lzacc )
    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    energyCheck = energy_in_surface(p_patch, ice, ssh(:,:), sst(:,:), computation_type=1, &
                                    dbg_lev=energyCheck_dbg_lev)

    CALL dbg_print('IceSlow: energy aft. ThickLimiter', energyCheck, str_module, &
         &         energyCheck_dbg_print, in_subset= p_patch%cells%all)

    !---------------------------------------------------------------------

    IF (ltimer) CALL timer_stop(timer_ice_slow)

  END SUBROUTINE ice_slow_thermo

  !-------------------------------------------------------------------------------
  !>
  !! ! ice_open_ocean: Calculates ice growth in open ocean and corresponding heat flux
  !!
  !! This function changes:
  !! ice % newice       new ice growth in open water                            [m]
  !! ice % heatOceW     heat flux to the ocean surface layer (via open surface) [W/m^2]
  !!
  SUBROUTINE ice_open_ocean(p_patch, ice, atmos_fluxes, sst, lacc)
    TYPE(t_patch),TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),           INTENT(INOUT) :: ice
    TYPE(t_atmos_fluxes),      INTENT(IN)    :: atmos_fluxes
    REAL(wp),                  INTENT(IN)    :: sst(:,:)
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    REAL(wp), DIMENSION (nproma, p_patch%alloc_cell_blocks) ::   &
      & AvailMLHeat    ! heat content of the ocean's surface layer availabe for cooling    [J/m^2]

    ! Loop indices
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: jb, jc, i_startidx_c, i_endidx_c
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-------------------------------------------------------------------------------------------
    all_cells   => p_patch%cells%all

    !$ACC DATA CREATE(AvailMLHeat) IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    AvailMLHeat(:,:) = 0.0_wp
    !$ACC END KERNELS
    !-------------------------------------------------------------------------------
    CALL dbg_print('IceOpenOcean bef: SST',     sst,            str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceOpenOcean bef: newice',  ice%newice,     str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceOpenOcean bef: heatOceW',ice%heatOceW,   str_module, 4, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------

    !TODOram: openmp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c
        IF (all_cells%vertical_levels(jc,jb) < 1) CYCLE ! Ocean points only

        ! aggregate atm fluxes over open ocean (positive downward, into the ocean)
        ice%heatOceW(jc,jb) = ( atmos_fluxes%SWnetw(jc,jb) + atmos_fluxes%LWnetw(jc,jb) &
          &                   + atmos_fluxes%sensw (jc,jb) + atmos_fluxes%latw  (jc,jb))&
                              * (1.0_wp-ice%concSum(jc,jb))

        ! Calculate availabe heat content of the surface layer
        ! Note, that the total energy of warm water has already been (partly) depleeted by heatOceI
        AvailMLHeat(jc,jb) = (sst(jc,jb)-ice%Tfw(jc,jb))*ice%zUnderIce(jc,jb)*clw*rho_ref &
          &                   + dtime*sum(ice%heatOceI(jc,:,jb))

        ! Calculate possible super-cooling of the surface layer by ice%heatOceW
!dn        ice%newice(jc,jb) = ( -(sst(jc,jb)-ice%Tfw(jc,jb)) * ice%zUnderIce(jc,jb)*clw*rho_ref &
!dn             &                     - dtime*ice%heatOceW(jc,jb) )/( alf*rhoi )
        ice%newice(jc,jb) = ( -AvailMLHeat(jc,jb) - dtime*ice%heatOceW(jc,jb)  )/( alf*rhoi )

        ! New ice forms only if the ML heat content is depleeted, newice > 0
        ! ice%heatOceW is reduced by *(1-conc), therefore newice is value averages over whole grid area
        ice%newice(jc,jb) = MAX(ice%newice(jc,jb), 0.0_wp)

        ! If no new ice was formed, then ice%heatOceW(jc,jb) will be passed to the ocean as is
        ! Otherwise, AvailMLHeat was FULLY depleted and heatOceW = flux required to cool the ocean to the freezing point
        IF ( ice%newice(jc,jb) > 0.0_wp ) THEN
            ice%heatOceW(jc,jb) = - AvailMLHeat(jc,jb) / dtime
          ! #slo# 2015-01: bugfix for supercooled ocean: whole grid area is cooled below freezing point,
          !                no multiplication by 1-conc is required to calculate the part of ice%heatOceW due to supercooling
          !                but check: doubled supercooling, since it is already taken into account in oce_ice_heatflx?
!dn            ice%heatOceW(jc,jb) = (ice%Tfw(jc,jb)-sst(jc,jb)) * ice%zUnderIce(jc,jb)*clw*rho_ref / dtime
        ENDIF

      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !$ACC END DATA

    !-------------------------------------------------------------------------------
    CALL dbg_print('IceOpenOcean aft: newice',  ice%newice,     str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceOpenOcean aft: heatOceW',ice%heatOceW,   str_module, 3, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------

  END SUBROUTINE ice_open_ocean

  !-------------------------------------------------------------------------------
  !>
  !! ! update_ice_ocean_fluxes: Updates ocean fluxes
  !!   ToDo: reformulate in a way, that doesn't depend on ice%conc
  !!         then it can be applied after all ice-slow routines (and not before ice_conc_change, as now)
  !!
  !! This function changes:
  !! p_oce_sfc % HeatFlux_Total      sum of forcing surface heat flux                          [W/m2]
  !! p_oce_sfc % HeatFlux_ShortWave  shortwave heat flux for penetration into deeper layers    [W/m2]
  !! p_oce_sfc % HeatFlux_LongWave   surface long wave heat flux                               [W/m2]
  !! p_oce_sfc % HeatFlux_Sensible   surface sensible heat flux                                [W/m2]
  !! p_oce_sfc % HeatFlux_Latent     surface latent heat flux                                  [W/m2]
  !! p_oce_sfc % FrshFlux_VolumeIce  forcing volume flux for height equation under sea ice     [m/s]
  !! p_oce_sfc % FrshFlux_TotalIce   forcing surface freshwater flux due to sea ice change     [m/s]
  !!
  SUBROUTINE update_ice_ocean_fluxes(p_patch, ice, atmos_fluxes, sss, p_oce_sfc, lacc)
    TYPE(t_patch),TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),           INTENT(IN)    :: ice
    TYPE(t_atmos_fluxes),      INTENT(IN)    :: atmos_fluxes
    REAL(wp),                  INTENT(IN)    :: sss(:,:)
    TYPE(t_ocean_surface),     INTENT(INOUT) :: p_oce_sfc
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    !Local Variables
    REAL(wp), DIMENSION (nproma, p_patch%alloc_cell_blocks) ::   &
      & fi1,       &  ! fr water flux due to ice growth/melt plus snow_to_ice conversion    [m/s]
      & fi2,       &  ! fr water flux due to snow melt only                                 [m/s]
      & fi3,       &  ! fr water flux due to new ice formation                              [m/s]
      & snowiceave,&  ! average snow to ice conversion within a grid cell                   [m]
      & icegrowave,&  ! grid cell averaged thickness of grown/melted ice (water equiv.)     [m]
      & snowmelted    ! average change in snow thickness within a grid cell                 [m]

    ! Loop indices
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: jb, jc, i_startidx_c, i_endidx_c
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !-------------------------------------------------------------------------------------------
    all_cells   => p_patch%cells%all

    !$ACC DATA CREATE(fi1, fi2, fi3, snowiceave, icegrowave, snowmelted) IF(lzacc)

    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    fi1         (:,:) = 0.0_wp
    fi2         (:,:) = 0.0_wp
    fi3         (:,:) = 0.0_wp
    snowiceave  (:,:) = 0.0_wp
    icegrowave  (:,:) = 0.0_wp
    snowmelted  (:,:) = 0.0_wp
    !$ACC END KERNELS
    !-------------------------------------------------------------------------------

  !TODOram: openmp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c
        IF (all_cells%vertical_levels(jc,jb) < 1) CYCLE ! Ocean points only

        ! Collect the atm fluxes into p_oce_sfc
        p_oce_sfc%HeatFlux_ShortWave(jc,jb) = atmos_fluxes%SWnetw(jc,jb)*(1.0_wp-ice%concSum(jc,jb))
        p_oce_sfc%HeatFlux_LongWave (jc,jb) = atmos_fluxes%LWnetw(jc,jb)*(1.0_wp-ice%concSum(jc,jb))
        p_oce_sfc%HeatFlux_Sensible (jc,jb) = atmos_fluxes%sensw (jc,jb)*(1.0_wp-ice%concSum(jc,jb))
        p_oce_sfc%HeatFlux_Latent   (jc,jb) = atmos_fluxes%latw  (jc,jb)*(1.0_wp-ice%concSum(jc,jb))

        ! heatOceI is calculated in ice_growth_*; heatOceW is calculated in ice_open_ocean
        p_oce_sfc%HeatFlux_Total(jc,jb) = SUM(ice%heatOceI(jc,:,jb)) + ice%heatOceW(jc,jb)

        ! Volmue flux (rain that goes through ice)
        p_oce_sfc%FrshFlux_VolumeIce(jc,jb) = atmos_fluxes%rprecw(jc,jb)*ice%concSum(jc,jb)

        ! FrshFlux_TotalIce contains fluxes that do not change the total fresh-water volume in the cell
        IF (v_base%lsm_c(jc,1,jb) <= sea_boundary ) THEN

          snowiceave(jc,jb) = SUM( ice%snow_to_ice(jc,:,jb)*ice%conc(jc,:,jb) )
          ! thickness of ice growth/melt + new ice in open water + snow_to_ice conversion (in cell-average water equiv.)
          icegrowave(jc,jb) = ( (SUM(ice%delhi(jc,:,jb)) + ice%newice(jc,jb))*rhoi + snowiceave(jc,jb)*rhos ) / rho_ref
          ! total salt flux from ocean to sea ice:
          p_oce_sfc%FrshFlux_IceSalt(jc,jb) = - sice * icegrowave(jc,jb) / dtime

          ! The treatment of snow-to-ice conversion is done as follows: The amount of snow that is transformed to ice is melted
          ! This meltwater flux is included in "snowmelted" and thus in fi2.
          ! The energy for the snow melting is then used to form an equivalent amount of sea ice from sea water.
          ! Hence, energy is conserved. The related salt flux of the sea-ice formation from sea water is included in fi1.
          ! Growing ice from ocean water in fi1 and melted snow in fi2 which are not identical due to sea ice salinity.

          ! ice growth/melt plus snow_to_ice conversion proportional to salt difference of water and ice
          IF (vert_cor_type == 0) THEN
          fi1(jc,jb) = - (1._wp-sice/sss(jc,jb)) * ( SUM( ice%delhi(jc,:,jb) )*rhoi + snowiceave(jc,jb)*rhos )/(rho_ref*dtime)
          ELSEIF (vert_cor_type == 1) THEN
            fi1(jc,jb) = - ( SUM( ice%delhi(jc,:,jb) )*rhoi + snowiceave(jc,jb)*rhos )/(rho_ref*dtime)
          ENDIF

          ! Total snow melt (negative is net melt) is given by real snow melt and implied snow melt from snow-to ice conversion.
          ! Real snow melt is given by change in snow thickness, minus the new snow from snow fall.

          snowmelted(jc,jb) = ( SUM( ice%delhs(jc,:,jb) )- snowiceave(jc,jb))*rhos/rho_ref-ice%totalsnowfall(jc,jb)
          fi2(jc,jb) = - MIN(snowmelted(jc,jb), 0.0_wp)/dtime

          ! new ice formation
          IF (vert_cor_type == 0) THEN
          fi3(jc,jb) =-(1._wp-sice/sss(jc,jb))*ice%newice(jc,jb)*rhoi/(rho_ref*dtime)
          ELSEIF (vert_cor_type == 1) THEN
           fi3(jc,jb) = - ice%newice(jc,jb)*rhoi/(rho_ref*dtime)
         ENDIF

          ! total freshwater flux from sea ice thermodynamics:
          p_oce_sfc%FrshFlux_TotalIce(jc,jb) = fi1(jc,jb) + fi2(jc,jb) + fi3(jc,jb)

!          p_oce_sfc%FrshFlux_TotalIce(jc,jb) = -( ( SUM( ice%delhi(jc,:,jb))*rhoi &
!            & + SUM( ice%delhs(jc,:,jb) )*rhos + ice%newice(jc,jb)*rhoi)/rho_ref &
!            & - ice%totalsnowfall(jc,jb) )/dtime


        ENDIF

      END DO
      !$ACC END PARALLEL LOOP
    END DO

    !$ACC END DATA

    !-------------------------------------------------------------------------------
    CALL dbg_print('Update_IO_FL: SSS       ', sss           ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: heatOceI  ', ice%heatOceI  ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: heatOceW  ', ice%heatOceW  ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: fi1       ', fi1          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: fi2       ', fi2          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: fi3       ', fi3          ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: snowmelted',snowmelted    ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: snowiceave',snowiceave    ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: icegrowave',icegrowave    ,str_module, 4, in_subset=p_patch%cells%owned)

    CALL dbg_print('Update_IO_FL: TotalHeat', p_oce_sfc%HeatFlux_Total,     str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: VolumeIce', p_oce_sfc%FrshFlux_VolumeIce, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: TotalIce ', p_oce_sfc%FrshFlux_TotalIce , str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('Update_IO_FL: IceSalt  ', p_oce_sfc%FrshFlux_IceSalt  , str_module, 3, in_subset=p_patch%cells%owned)
    !-------------------------------------------------------------------------------

  END SUBROUTINE update_ice_ocean_fluxes

  !-------------------------------------------------------------------------------
  !>
  !! !! ice_update_conc: Calculates changes in concentration, update ice variables
  !! Result strongly depends on "lead-close" parameters (see Eq. 16 of Hibler 1979)
  !!
  !! This function changes:
  !! ice % conc     new snow thickness for each ice category                [ ]
  !! ice % concSum  new ice  thickness for each ice category                [ ]
  !!
  SUBROUTINE ice_update_conc(p_patch,ice)

    TYPE(t_patch),             INTENT(IN)    :: p_patch
    TYPE (t_sea_ice),          INTENT(INOUT) :: ice

    INTEGER  :: k

    ! TODO ram - remove all instances of p_patch%cells%area(:,:) and test
    ! change in concentration does not change ice and snow volumes:

    ! new snow and ice volumes, should be conserved
    DO k=1,ice%kice
      ! newice must not be multiplied by 1-conc, which is already included
      ice%vol (:,k,:) = ( ice%hi(:,k,:)*ice%conc(:,k,:) + ice%newice(:,:) ) * p_patch%cells%area(:,:)
      ice%vols(:,k,:) = ice%hs(:,k,:)*ice%conc(:,k,:) * p_patch%cells%area(:,:)
    ENDDO

    CALL dbg_print('IceConcCh: vol  at beg' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at beg' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

!ICON_OMP_PARALLEL
!ICON_OMP_WORKSHARE
    WHERE ( v_base%lsm_c(:,1,:) <= sea_boundary )

      ice%conc(:,1,:) = ice%conc(:,1,:) + &
      ! Hibler's way to change the concentration: conc=conc+newice/hnull
      ! New formulation: leadclose parameters 2 and 3 of MPIOM are combined into leadclose_2n:
      ! leadclose_2n = mpiom_leadclose(3)/mpiom_leadclose(2). According to mpiom default leadclose_2n=2/3.
      ! due to new ice formation (newice >= 0)
        &               ice%newice(:,:)/(hnull+leadclose_2n*(ice%hi(:,1,:)-hnull))      +    &
      ! due to lateral ice melting (delhi < 0)
        &               leadclose_1 * MIN(ice%delhi(:,1,:), 0._wp) / MAX(ice%hi(:,1,:), hmin)

      ! 0 <= ice%conc <= 1
      ice%conc(:,1,:) = MIN( 1._wp, MAX( 0._wp, ice%conc(:,1,:) ) )
    ENDWHERE
!ICON_OMP_END_WORKSHARE

!ICON_OMP_WORKSHARE
      ! New ice and snow thickness
    WHERE ( ice%conc(:,1,:) > 0._wp )
      ice%hi  (:,1,:) = ice%vol (:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
      ice%hs  (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
    ENDWHERE
      !TODO: Re-calculate temperatures to conserve energy when we change the ice thickness
!ICON_OMP_END_WORKSHARE

    ! Ice cannot grow thinner than hmin
!ICON_OMP_WORKSHARE
    WHERE ( ice%hi(:,1,:) < hmin .AND. ice%hi(:,1,:) > 0._wp )
      ice%hi  (:,1,:) = hmin
      ice%conc(:,1,:) = ice%vol(:,1,:) /( ice%hi  (:,1,:)*p_patch%cells%area(:,:) )
      ice%hs  (:,1,:) = ice%vols(:,1,:)/( ice%conc(:,1,:)*p_patch%cells%area(:,:) )
    ENDWHERE
!ICON_OMP_END_WORKSHARE

!ICON_OMP_WORKSHARE
!   ! should not be necessary
!    WHERE (ice%hi(:,1,:) <= 0._wp)
!      ice%Tsurf(:,1,:) = ice%Tfw(:,:)
!      ice%hi   (:,1,:) = 0.0_wp
!      ice%hs   (:,1,:) = 0.0_wp
!      ice%conc (:,1,:) = 0.0_wp
!      ice%vol  (:,1,:) = 0.0_wp
!      ice%vols  (:,1,:) = 0.0_wp
!    ENDWHERE
!ICON_OMP_END_WORKSHARE
!ICON_OMP_END_PARALLEL

    ice%concSum(:,:)  = SUM(ice%conc(:,:,:),2)

    CALL dbg_print('IceConcCh: IceConc end' ,ice%conc, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hi   at end' ,ice%hi  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: hs   at end' ,ice%hs  , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vol  at end' ,ice%vol , str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceConcCh: vols at end' ,ice%vols, str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_update_conc

   !-------------------------------------------------------------------------
  !>
  !! !  ice_thickness_limiter: Limit ice thinkness. Adjust fluxes and calculates the final freeboard
  !! !
  SUBROUTINE ice_thickness_limiter( p_patch, p_ice, p_oce_sfc, p_os, lacc )
    TYPE(t_patch),TARGET,      INTENT(IN)    :: p_patch
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_ice
    TYPE(t_ocean_surface),     INTENT(INOUT) :: p_oce_sfc
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    TYPE(t_subset_range), POINTER                           :: all_cells
    INTEGER                                                 :: k, jb, jc, i_startidx_c, i_endidx_c
    REAL(wp)                                                :: z_smax, prism_thick_flat
    REAL(wp), DIMENSION (nproma, p_patch%alloc_cell_blocks) :: sss ! Sea surface salinity
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    all_cells    => p_patch%cells%all

    !$ACC DATA CREATE(sss) IF(lzacc)

    ! surface layer thickness, same as patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,1,:)
    prism_thick_flat = v_base%del_zlev_m(1)
    ! sea surface salinity
    !$ACC KERNELS DEFAULT(PRESENT) IF(lzacc)
    sss(:,:)  =  p_os%p_prog(nold(1))%tracer(:,1,:,2)
    !$ACC END KERNELS

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c

        DO k = 1, p_ice%kice

        ! limit sea ice thickness to seaice_limit of surface layer depth, without elevation
            z_smax = seaice_limit * prism_thick_flat

            IF (vert_cor_type .eq. 1) THEN
              z_smax = z_smax*p_os%p_prog(nold(1))%stretch_c(jc, jb)
            END IF
            IF (seaice_limit_abs > 0._wp .AND. z_smax > seaice_limit_abs) THEN
              z_smax = seaice_limit_abs
            ENDIF

            IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary  .AND.  p_ice%hi(jc,k,jb) > z_smax ) THEN
                ! Tracer flux due to removal
                IF (vert_cor_type .EQ. 1 ) THEN  ! FrshFlx_Total_ice is used differently here
                  p_oce_sfc%FrshFlux_TotalIce (jc,jb) = p_oce_sfc%FrshFlux_TotalIce (jc,jb)                      &
                        & + (p_ice%hi(jc,k,jb)-z_smax)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime)  ! Ice
                ELSE
                  p_oce_sfc%FrshFlux_TotalIce (jc,jb) = p_oce_sfc%FrshFlux_TotalIce (jc,jb)                      &
                      & + (1._wp-sice/sss(jc,jb))*(p_ice%hi(jc,k,jb)-z_smax)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime)  ! Ice

                ENDIF

                p_oce_sfc%FrshFlux_IceSalt(jc,jb) = p_oce_sfc%FrshFlux_IceSalt(jc,jb) &
                  & + sice * (p_ice%hi(jc,k,jb)-z_smax)*p_ice%conc(jc,k,jb)*rhoi/(rho_ref*dtime)

                ! Heat flux due to removal
                !  #slo# 2015-02: - this heat did not come from the ocean, but from atmosphere, heating ocean is wrong
                !                 - check if conc must enter here as well, check energy for coupling
                p_oce_sfc%HeatFlux_Total(jc,jb) = p_oce_sfc%HeatFlux_Total(jc,jb)   &
                      & + (p_ice%hi(jc,k,jb)-z_smax)*p_ice%conc(jc,k,jb)*alf*rhoi/dtime           ! Ice
                p_ice%hi  (jc,k,jb) = z_smax
                p_ice%vol (jc,k,jb) = p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)*p_patch%cells%area(jc,jb)
            ENDIF

            p_ice%draft(jc,k,jb)   = (rhos * p_ice%hs(jc,k,jb) + rhoi * p_ice%hi(jc,k,jb)) / rho_ref

        ENDDO

        p_ice%draftave (jc,jb) = sum(p_ice%draft(jc,:,jb) * p_ice%conc(jc,:,jb))

      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO

    !$ACC END DATA

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('iceClUp: hi aft. limiter'     ,p_ice%hi       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('iceClUp: hs aft. limiter'     ,p_ice%hs       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: HeatTotal aft. limit',p_oce_sfc%HeatFlux_Total   ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: TotalIce aft. limit' ,p_oce_sfc%FrshFlux_TotalIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: draft aft. limit'    ,p_ice%draft    ,str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_thickness_limiter

  !-------------------------------------------------------------------------
  !>
  !! !  ice_thickness_limiter: Limit ice thinkness. Adjust fluxes and calculates the final freeboard
  !! !
  SUBROUTINE ice_thickness_limiter_hh( p_patch_3d, p_ice, p_oce_sfc, p_os, lacc )
    TYPE(t_patch_3d),TARGET,      INTENT(IN) :: p_patch_3d
    TYPE(t_sea_ice),           INTENT(INOUT) :: p_ice
    TYPE(t_ocean_surface),     INTENT(INOUT) :: p_oce_sfc
    TYPE(t_hydro_ocean_state), INTENT(IN)    :: p_os
    LOGICAL, INTENT(IN), OPTIONAL            :: lacc

    ! Local variables
    TYPE(t_patch), POINTER                   :: p_patch
    TYPE(t_subset_range), POINTER            :: all_cells
    INTEGER                                  :: k, jb, jc, i_startidx_c, i_endidx_c
    REAL(wp)                                 :: z_smax ,rhoicwa,zunderice_old,zunderice_new,ddice &
                                               ,svnew,svold,draft_old
    LOGICAL :: lzacc
    CHARACTER(len=*), PARAMETER :: routine = 'ice_thickness_limiter_hh'

    CALL set_acc_host_or_device(lzacc, lacc)


    p_patch => p_patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    rhoicwa = rhoi / rho_ref

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c

        DO k = 1, p_ice%kice

          ! limit sea ice thickness to seaice_limit of surface layer depth
          z_smax = seaice_limit * p_patch_3d%p_patch_1d(1)%prism_thick_c(jc,1,jb)
          IF ( v_base%lsm_c(jc,1,jb) <= sea_boundary  &
               .AND.  (rhoicwa*p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb) ) > z_smax ) THEN

            draft_old   = (rhos * p_ice%hs(jc,k,jb) + rhoi * p_ice%hi(jc,k,jb)) / rho_ref
            zunderice_old = p_patch_3d%p_patch_1d(1)%prism_thick_c(jc,1,jb) &
                 - draft_old*p_ice%conc(jc,k,jb)
            svold =sice*rhoicwa*(p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)) &
                       +p_os%p_prog(nold(1))%tracer(jc,1,jb,2)*zunderice_old

            ddice = MIN(  (draft_old*p_ice%conc(jc,k,jb)-z_smax)/rhoicwa &
                 , p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb) )

            p_ice%hi(jc,k,jb) = p_ice%hi(jc,k,jb)-ddice/p_ice%conc(jc,k,jb)
            p_ice%vol(jc,k,jb) = p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)
            p_ice%draft(jc,k,jb) = (rhos * p_ice%hs(jc,k,jb) + rhoi * p_ice%hi(jc,k,jb))/rho_ref
            p_ice%zunderice(jc,jb) = p_patch_3d%p_patch_1d(1)%prism_thick_c(jc,1,jb) &
                 - p_ice%draft(jc,k,jb)*p_ice%conc(jc,k,jb)

            !update surface salinity
            p_os%p_prog(nold(1))%tracer(jc,1,jb,2) =(rhoicwa*ddice*sice &
                             +p_os%p_prog(nold(1))%tracer(jc,1,jb,2)*zunderice_old)/p_ice%zunderice(jc,jb)

            svnew =sice*rhoicwa*(p_ice%hi(jc,k,jb)*p_ice%conc(jc,k,jb)) &
                       +p_os%p_prog(nold(1))%tracer(jc,1,jb,2)*zunderice_new

            print*,'seaice-limiter',svnew-svold,svnew,svold,ddice,&
                            zunderice_new,zunderice_old

!            p_oce_sfc%FrshFlux_TotalIce (jc,jb) = p_oce_sfc%FrshFlux_TotalIce (jc,jb)                      &
!                 & + ddice*rhoi/(rho_ref*dtime)  !
!            p_oce_sfc%HeatFlux_Total(jc,jb) = p_oce_sfc%HeatFlux_Total(jc,jb)   &
!                 & + ddice*alf*rhoi/dtime           ! Ice

          ENDIF

       ENDDO

       p_ice%draftave (jc,jb) = sum(p_ice%draft(jc,:,jb) * p_ice%conc(jc,:,jb))

      ENDDO
      !$ACC END PARALLEL LOOP
    ENDDO

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('iceClUp: hi aft. limiter'     ,p_ice%hi       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('iceClUp: hs aft. limiter'     ,p_ice%hs       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: HeatTotal aft. limit',p_oce_sfc%HeatFlux_Total   ,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: TotalIce aft. limit' ,p_oce_sfc%FrshFlux_TotalIce,str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('IceClUp: draft aft. limit'    ,p_ice%draft    ,str_module, 4, in_subset=p_patch%cells%owned)

  END SUBROUTINE ice_thickness_limiter_hh

  !-------------------------------------------------------------------------------
  !>
  !! ! ice_growth_init: save ice thickness before ice_growth_* and initialize variables
  !!
  SUBROUTINE ice_growth_init (ice, p_patch, lacc)

    TYPE (t_sea_ice), INTENT(INOUT) :: ice
    LOGICAL, INTENT(IN), OPTIONAL   :: lacc
    TYPE(t_patch),  POINTER         :: p_patch
    INTEGER                         :: jb, k, jc, i_startidx_c, i_endidx_c

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! initialize ice-growth-related variables with zeros
    DO jb = 1, p_patch%nblks_c
      CALL get_index_range(p_patch%cells%all, jb, i_startidx_c, i_endidx_c)
      DO k=1,ice%kice
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          ice % zHeatOceI   (jc, k, jb) = 0._wp
          ice % heatOceI    (jc, k, jb) = 0._wp
          ice % snow_to_ice (jc, k, jb) = 0._wp
          ice % delhi       (jc, k, jb) = 0._wp         ! change in mean ice thickness due to thermodynamic effects
          ice % delhs       (jc, k, jb) = 0._wp         ! change in mean snow thickness due to melting
          ice % hiold (jc, k, jb) = ice%hi(jc, k, jb)
!    ice % hsold (:,:,:) = ice%hs(:,:,:)
          IF (k == 1) THEN
            ice % heatOceW      (jc, jb) = 0._wp
            ice % newice        (jc, jb) = 0._wp
            ice % totalsnowfall (jc, jb) = 0._wp
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END DO
    END DO

  END SUBROUTINE ice_growth_init

END MODULE mo_ice_slow_thermo
