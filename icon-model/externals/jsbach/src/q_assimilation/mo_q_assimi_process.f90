!> QUINCY routines for photosynthesis related processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### routines for photosynthesis related processes, e.g., calculation of photosynthesis
!>
MODULE mo_q_assimi_process
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_math_constants,    ONLY: eps4, eps8
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, message_text, finish

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_canopy_fluxes, calc_photosynthesis, calc_soil_moisture_stress
  PUBLIC :: discrimination_ps

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_assimi_process'

CONTAINS

  ! ======================================================================================================= !
  !>update canopy fluxes, includes calculation of photosynthesis
  !>
  SUBROUTINE update_canopy_fluxes(tile, options)

    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,       ONLY: t_lctlib_element
    USE mo_jsb_process_class,      ONLY: A2L_, Q_ASSIMI_, VEG_, Q_RAD_, SPQ_
    USE mo_jsb_grid_class,         ONLY: t_jsb_vgrid
    USE mo_jsb_grid,               ONLY: Get_vgrid
    USE mo_atmland_constants,      ONLY: min_wind, max_wind
    USE mo_jsb_physical_constants, ONLY: r_gas, Tzero
    USE mo_isotope_util,           ONLY: calc_fractionation
    USE mo_phy_schemes             ! e.g., heat_transfer_coef
    USE mo_q_assimi_parameters,    ONLY: jmax2n, vcmax2n, E0v, E1v, t_jmax_offset, t_jmax_slope, t_jmax_opt_min, t_jmax_opt_max
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_ASSIMI_)
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_RAD_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model                    !< the model
    TYPE(t_lctlib_element), POINTER       :: lctlib                   !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_canopy_q_assimi    !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb            !< Vertical grid for soil layers
    INTEGER                               :: nsoil_sb                 !< number of soil layers
    INTEGER                               :: ncanopy                  !< number of canopy layers
    INTEGER                               :: icanopy                  !< for looping across ncanopy layers
    REAL(wp)                              :: lctlib_g1                !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    REAL(wp), DIMENSION(options%nc)       :: hlp1, hlp2, hlp3         !< dummy variables
    INTEGER                               :: iblk, ics, ice, nc       !< dimensions
    INTEGER                               :: ic, isoil                !< loop counter
    REAL(wp)                              :: steplen                  !< step length
    REAL(wp)                              :: alpha                    !< implicitness factor
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_canopy_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: q_air
    dsl4jsb_Real2D_onChunk      :: press_srf
    dsl4jsb_Real2D_onChunk      :: wind_10m
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_C13
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_C14
    dsl4jsb_Real2D_onChunk      :: swvis_srf_down
    dsl4jsb_Real2D_onChunk      :: swnir_srf_down
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: net_assimilation
    dsl4jsb_Real2D_onChunk      :: net_assimilation_boc
    dsl4jsb_Real2D_onChunk      :: maint_respiration_leaf
    dsl4jsb_Real2D_onChunk      :: aerodyn_cond
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: co2_conc_leaf
    dsl4jsb_Real2D_onChunk      :: beta_air
    dsl4jsb_Real2D_onChunk      :: beta_soa
    dsl4jsb_Real2D_onChunk      :: soa_tsoa_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt
    ! Q_ASSIMI_ 3D
    dsl4jsb_Real3D_onChunk      :: net_assimilation_cl
    dsl4jsb_Real3D_onChunk      :: gross_assimilation_cl
    dsl4jsb_Real3D_onChunk      :: maint_respiration_leaf_cl
    dsl4jsb_Real3D_onChunk      :: canopy_cond_cl
    dsl4jsb_Real3D_onChunk      :: ftranspiration_sl
    dsl4jsb_Real3D_onChunk      :: co2_conc_leaf_cl
    dsl4jsb_Real3D_onChunk      :: jmax_cl
    dsl4jsb_Real3D_onChunk      :: vcmax_cl
    dsl4jsb_Real3D_onChunk      :: chlfl_yield_cl
    ! Q_RAD_ 3D
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: spq_drag_srf
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot
    dsl4jsb_Real3D_onChunk      :: t_soil_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_pwp_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk      :: w_soil_pot_sl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tacclim_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen
    alpha   = options%alpha
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model                  => Get_model(tile%owner_model_id)
    lctlib                 => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy_q_assimi  => Get_vgrid('q_canopy_layer')
    ncanopy                =  vgrid_canopy_q_assimi%n_levels
    vgrid_soil_sb          => Get_vgrid('soil_layer_sb')
    nsoil_sb               =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)                          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)                      ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, wind_10m)                       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)               ! in     ! ICON-Land: calculated from jsb4 co2_air_mol
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio_C13)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio_C14)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swvis_srf_down)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swnir_srf_down)                 ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation)        ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C13)    ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C14)    ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, net_assimilation)          ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, net_assimilation_boc)      ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, maint_respiration_leaf)    ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, aerodyn_cond)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, canopy_cond)               ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, ftranspiration_sl)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, co2_conc_leaf)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air)                  ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soa)                  ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, soa_tsoa_mavg)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps)              ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, t_jmax_opt)                ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, net_assimilation_cl)       ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, gross_assimilation_cl)     ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, maint_respiration_leaf_cl) ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, canopy_cond_cl)            ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, co2_conc_leaf_cl)          ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, jmax_cl)                   ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, vcmax_cl)                  ! out
    dsl4jsb_Get_var3D_onChunk(Q_ASSIMI_, chlfl_yield_cl)            ! out
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_cl)               ! in
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_cl)               ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_pot)                ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, spq_drag_srf)                   ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, t_soil_sl)                      ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_pwp_sl)                  ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_fc_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_, w_soil_pot_sl)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tacclim_mavg)   ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                         ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)                ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                      ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                       ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                      ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                     ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)               ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)               ! in
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9.1 set output to zero
    !>
    canopy_cond(:)            = 0.0_wp
    co2_conc_leaf(:)          = 0.0_wp

    !>0.9.2 set g1 according to canopy_conductance_scheme
    !>
    SELECT CASE(TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme))
      CASE ("medlyn")
        lctlib_g1 = lctlib%g1_medlyn
      CASE ("ballberry")
        lctlib_g1 = lctlib%g1_bberry
    END SELECT

    !>1.0 calculate the moisture scaling factor for photosynthesis
    !>
    beta_air(:)     = calc_atm_moisture_stress(TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme), &
                                               t_air(:), q_air(:), press_srf(:))
    beta_soil_ps(:) = calc_soil_moisture_stress(w_soil_root_pot(:), lctlib%phi_leaf_min)
    beta_soil_gs(:) = calc_soil_moisture_stress(w_soil_root_pot(:), lctlib%phi_leaf_min)
    ! catch beta_soil_ps getting zero
    !  with "beta_soil_ps=0" calc_photosynthesis() does get a runtime error, because n1 & n2 become zero,
    !  and then msat & nco become zero, resulting in an error with the log() operation "nlim = -log(nco/..."
    WHERE(beta_soil_ps(:) < eps8) beta_soil_ps(:) = eps8

    !> 1.1 determine fraction of transpiration per soil layer
    !>
    ftranspiration_sl(:,:) = 0.0_wp
    DO ic = 1, nc
      DO isoil = 1, nsoil_sb
        IF ((w_soil_fc_sl(ic,isoil) - w_soil_pwp_sl(ic,isoil)) > eps8 .AND. t_soil_sl(ic,isoil) > Tzero) THEN
          hlp1(ic) = calc_soil_moisture_stress(w_soil_pot_sl(ic, isoil), lctlib%phi_leaf_min)
          ftranspiration_sl(ic,isoil) = root_fraction_sl(ic,isoil) * hlp1(ic)
        END IF
      END DO
      IF (SUM(ftranspiration_sl(ic,:)) > eps8) THEN
        ftranspiration_sl(ic,:) = ftranspiration_sl(ic,:) / SUM(ftranspiration_sl(ic,:))
      ELSE
        ftranspiration_sl(ic,:) = 0.0_wp
      END IF
    END DO

    !> 1.2 determine aerodynamic conductance
    !>
    DO ic = 1, nc
      aerodyn_cond(ic) = heat_transfer_coef(spq_drag_srf(ic), steplen, alpha) * MIN( MAX(wind_10m(ic), min_wind), max_wind)
    END DO

    ! calc scaling factor to account for the 'state of acclimation' of evergreen conifers in spring
    beta_soa(:) = calc_beta_soa(soa_tsoa_mavg(:),lctlib%phenology_type)

    !>2.0 calculate photosynthesis per layer
    !>
    ! calc photosynthesis
    DO icanopy = 1,ncanopy
      CALL calc_photosynthesis(                                       &
        & lctlib%gmin                                               , & ! lctlib
        & lctlib%g0                                                 , &
        & lctlib_g1                                                 , &
        & lctlib%t_jmax_omega                                       , &
        & lctlib%ps_pathway                                         , &
        & TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme) , & ! config
        & t_air(:)                                                  , & ! in
        & press_srf(:)                                              , &
        & co2_mixing_ratio(:)                                       , &
        & aerodyn_cond(:)                                           , &
        & t_air_tacclim_mavg(:)                                     , &
        & ppfd_sunlit_cl(:,icanopy)                                 , &
        & ppfd_shaded_cl(:,icanopy)                                 , &
        & fleaf_sunlit_cl(:,icanopy)                                , &
        & beta_air(:)                                               , &
        & beta_soa(:)                                               , &
        & beta_soil_ps(:)                                           , &
        & beta_sinklim_ps_tacclim_mavg(:)                           , &
        & beta_soil_gs(:)                                           , &
        & fn_chl_cl(:,icanopy)                                      , &
        & fn_et_cl(:,icanopy)                                       , &
        & fn_rub_cl(:,icanopy)                                      , &
        & fn_pepc_cl(:,icanopy)                                     , &
        & lai_cl(:,icanopy)                                         , & ! in - used for: IF statement whether this routine may run or not
        & leaf_nitrogen_cl(:,icanopy)                               , & ! in
        & t_jmax_opt_mavg(:)                                        , & ! in
        & gross_assimilation_cl(:,icanopy)                          , & ! out
        & net_assimilation_cl(:,icanopy)                            , & ! out
        & maint_respiration_leaf_cl(:,icanopy)                      , & ! out
        & canopy_cond_cl(:,icanopy)                                 , & ! out
        & co2_conc_leaf_cl(:,icanopy)                               , & ! out
        & hlp1(:)                                                   , & ! out  (output not needed here)
        & hlp2(:)                                                   , & ! out  (output not needed here)
        & hlp3(:)                                                   , & ! out  (output not needed here)
        & chlfl_yield_cl(:,icanopy)                                    )! out

      ! gross_assimilation incl. C13 & C14
      gross_assimilation(:) = gross_assimilation(:) + gross_assimilation_cl(:,icanopy) * lai_cl(:,icanopy)
      WHERE(gross_assimilation_cl(:,icanopy) > eps8)
        ! C13 flux associated with photosynthesis
        hlp1(:)                   = discrimination_ps(co2_conc_leaf_cl(:,icanopy), &
                                                      co2_mixing_ratio(:)        , &
                                                      'C13'                      , &
                                                      lctlib%ps_pathway)
        gross_assimilation_C13(:) = gross_assimilation_C13(:) + &
                                    gross_assimilation_cl(:,icanopy) * &
                                    lai_cl(:,icanopy) * &
                                    calc_fractionation(co2_mixing_ratio(:), co2_mixing_ratio_C13(:), hlp1(:))
        ! C14 flux associated with photosynthesis
        hlp1(:)                   = discrimination_ps(co2_conc_leaf_cl(:,icanopy), &
                                                      co2_mixing_ratio(:)        , &
                                                      'C14'                      , &
                                                      lctlib%ps_pathway)
        gross_assimilation_C14(:) = gross_assimilation_C14(:) + &
                                    gross_assimilation_cl(:,icanopy) * &
                                    lai_cl(:,icanopy) * &
                                    calc_fractionation(co2_mixing_ratio(:), co2_mixing_ratio_C14(:), &
                                                       hlp1(:), co2_mixing_ratio_C13(:))
      END WHERE

      net_assimilation(:)       = net_assimilation(:) + net_assimilation_cl(:,icanopy) * lai_cl(:,icanopy)
      maint_respiration_leaf(:) = maint_respiration_leaf(:) + maint_respiration_leaf_cl(:,icanopy) * lai_cl(:,icanopy)
      canopy_cond(:)            = canopy_cond(:) + canopy_cond_cl(:,icanopy) * lai_cl(:,icanopy)
      co2_conc_leaf(:)          = co2_conc_leaf(:) + co2_conc_leaf_cl(:,icanopy) * lai_cl(:,icanopy)

      ! get the net photosynthesis of the lowest canopy layer at this timestep
      WHERE(lai_cl(:,icanopy) > 0.0_wp)
        net_assimilation_boc(:) = net_assimilation_cl(:,icanopy)
      END WHERE

      ! diagnose jmax and vcmax
      !> N specific rate of electron transport rate
      jmax_cl(:,icanopy) = jmax2n * fn_et_cl(:,icanopy) * &
                           EXP(-(((t_air(:) - Tzero - t_jmax_opt_mavg(:)) / lctlib%t_jmax_omega) ** 2.0_wp)) * &
                           beta_soil_ps(:) * beta_sinklim_ps_tacclim_mavg(:) * beta_soa(:) * &
                           leaf_nitrogen_cl(:,icanopy)

      !> N specific rate of Rubisco turnover
      vcmax_cl(:,icanopy) = vcmax2n * fn_rub_cl(:,icanopy) * EXP(E0v - E1v / (r_gas * t_air(:))) * &
                            beta_soil_ps(:) * beta_sinklim_ps_tacclim_mavg(:) * beta_soa(:) * &
                            leaf_nitrogen_cl(:,icanopy)

    ENDDO ! DO icanopy = 1,ncanopy

    !> 3.0 update mean canopy CO2 concentration in the leaf
    !>   was accumulated from the CO2 concentration per canopy layer above ('co2_conc_leaf_cl(:,icanopy) * lai_cl(:,icanopy)')
    !>
    WHERE(lai(:) > eps8)
      co2_conc_leaf(:) = co2_conc_leaf(:) / lai(:)
    END WHERE

    !> 4.0 uptake optimum temperature for electron transport
    !>
    IF(dsl4jsb_Config(Q_ASSIMI_)%flag_t_jmax_acclimation) THEN
      t_jmax_opt(:) = t_jmax_offset + t_jmax_slope * ( t_air(:) - Tzero )
      t_jmax_opt(:) = MIN(MAX(t_jmax_opt(:), t_jmax_opt_min), t_jmax_opt_max)
    ELSE
      t_jmax_opt(:) = lctlib%t_jmax_opt
    ENDIF

  END SUBROUTINE update_canopy_fluxes

  ! ======================================================================================================= !
  !>Calculation of net photosynthesis and canopy conductance
  !>
  !>  according to Friend & Kiang 2005, J Climate
  !>  accounts for temperature acclimation of JMAX according to Friend, 2010 J Exp Bot
  ELEMENTAL SUBROUTINE calc_photosynthesis( gmin                , &
                                            g0                  , &
                                            g1                  , &
                                            t_jmax_omega        , &
                                            ps_pathway          , &
                                            canopy_cond_scheme  , &
                                            t_air               , &
                                            press_srf           , &
                                            co2_mixing_ratio    , &
                                            aerodyn_cond        , &
                                            temp_acclim         , &
                                            ppfd_sunlit_cl      , &
                                            ppfd_shaded_cl      , &
                                            fleaf_sunlit_cl     , &
                                            beta_air            , &
                                            beta_soa            , &
                                            beta_soil_ps        , &
                                            beta_sinklim_ps     , &
                                            beta_soil_gs        , &
                                            fn_chl_cl           , &
                                            fn_et_cl            , &
                                            fn_rub_cl           , &
                                            fn_pepc_cl          , &
                                            lai_cl              , &   ! used for: IF statement whether this routine may run or not
                                            leaf_nitrogen_cl    , &
                                            t_jmax_opt          , &
                                            ag_cl               , &
                                            an_cl               , &
                                            maint_resp_cl       , &
                                            gs_cl               , &
                                            ci_cl               , &
                                            m_rub               , &
                                            m_et                , &
                                            m_pepc              , &
                                            chlfl_yield_cl         )

  USE mo_veg_constants,             ONLY: fmaint_rate_base
  USE mo_jsb_physical_constants,    ONLY: r_gas, Tzero
  USE mo_q_assimi_parameters,       ONLY: E0kc, E1kc, E0ko, E1ko, E0pcp, E1pcp, E0v, E1v, CiCa_default_C3, CiCa_default_C4, &
                                          ps_it_max, ci_max, alpha_i, Tref_pepc, Tbase_pepc, pO2, jmax2n, vcmax2n, ka, pepc2n, &
                                          chl2n
  USE mo_q_assimi_constants,        ONLY: Dwv2co2_air, Dwv2co2_turb, ic3phot, ic4phot
  USE mo_q_veg_respiration,         ONLY: temperature_response_respiration

  IMPLICIT NONE
  ! ---------------------------
  ! 0.1 InOut
  REAL(wp),                 INTENT(IN) :: gmin            , &     !< land-cover-type library parameter
                                          g0              , &     !< land-cover-type library parameter
                                          g1              , &     !< land-cover-type library parameter
                                          t_jmax_omega            !< land-cover-type library parameter
  INTEGER,                  INTENT(IN) :: ps_pathway              !< land-cover-type library parameter
  CHARACTER(len=*),         INTENT(IN) :: canopy_cond_scheme      !< canopy_conductance_scheme: medlyn / ballberry
  REAL(wp),                 INTENT(IN) :: t_air                   !< 2m air temperature [K]
  REAL(wp),                 INTENT(IN) :: press_srf               !< air pressure (Pa)
  REAL(wp),                 INTENT(IN) :: aerodyn_cond            !< ga -- aerodynamic conductance (m s-1)
  REAL(wp),                 INTENT(IN) :: temp_acclim             !< veg\%t_air_tacclim_mavg  average temperature for respiration acclimation
  REAL(wp),                 INTENT(IN) :: ppfd_sunlit_cl          !< PAR of sunlit foliage par per canopy layer (micro-mol m-2 s-1)
  REAL(wp),                 INTENT(IN) :: ppfd_shaded_cl          !< PAR of shaded foliage par per canopy layer (micro-mol m-2 s-1)
  REAL(wp),                 INTENT(IN) :: fleaf_sunlit_cl         !< fraction of canopy layer that is sunlit
  REAL(wp),                 INTENT(IN) :: co2_mixing_ratio        !< CO2 concentration of the free air (ppm)
  REAL(wp),                 INTENT(IN) :: beta_air                !< scaling factor to account for air humidity effects on conductance
  REAL(wp),                 INTENT(IN) :: beta_soa                !< scaling factor to account for spring in coniferous forests
  REAL(wp),                 INTENT(IN) :: beta_soil_ps            !< scaling factor to account for soil moisture constraints on photosynthesis
  REAL(wp),                 INTENT(IN) :: beta_sinklim_ps         !< scaling factor to account for sink limitation constraints on photosynthesis
  REAL(wp),                 INTENT(IN) :: beta_soil_gs            !< scaling factor to account for soil moisture constraints on conductance
  REAL(wp),                 INTENT(IN) :: fn_et_cl                !< fraction of foliar N allocated to electron transport
  REAL(wp),                 INTENT(IN) :: fn_rub_cl               !< fraction of foliar N allocated to Rubisco
  REAL(wp),                 INTENT(IN) :: fn_pepc_cl              !< fraction of foliar N allocated to PEPC complex (C4 only)
  REAL(wp),                 INTENT(IN) :: fn_chl_cl               !< fraction of foliar N allocated to Chlorophyll
  REAL(wp),                 INTENT(IN) :: lai_cl                  !< lai of the particular canopy layer
  REAL(wp),                 INTENT(IN) :: leaf_nitrogen_cl        !< total foliar N content
  REAL(wp),                 INTENT(IN) :: t_jmax_opt              !< temperature optimum of electron transport
  REAL(wp),                 INTENT(OUT) :: ag_cl                  !< gross photsynthesis of layer (micro-mol CO2 / m2 /s)
  REAL(wp),                 INTENT(OUT) :: an_cl                  !< net photsynthesis of layer (micro-mol CO2 / m2 /s)
  REAL(wp),                 INTENT(OUT) :: maint_resp_cl          !< net canopy layer respiration (micro-mol CO2 / m2 /s)
  REAL(wp),                 INTENT(OUT) :: gs_cl                  !< canopy conductance of layer (m / s)
  REAL(wp),                 INTENT(OUT) :: ci_cl                  !< leaf internal CO2 concentration of layer (ppm)
  REAL(wp),                 INTENT(OUT) :: m_rub                  !< CO2 limitation of Rubisco factor for calculating optimal N fraction
  REAL(wp),                 INTENT(OUT) :: m_et                   !< CO2 limitation of ET factor for calculating optimal N fraction
  REAL(wp),                 INTENT(OUT) :: m_pepc                 !< CO2 limitation for PEPC factor for calculating optimal N fraction
  REAL(wp),                 INTENT(OUT) :: chlfl_yield_cl         !< leaf level chlorophyll fluorescence yield (unitless)

  ! ---------------------------
  ! 0.2 Local
  REAL(wp) :: kc, ko, km                  !< Michaelis-Menten Kinetic parameters for Rubisco (Pa)
  REAL(wp) :: pcp                         !< Photosynthetic compensation point (Pa)
  REAL(wp) :: f_resp                      !< N-specific layer respiration rate (micro-mol CO2 / mmol N /s)
  REAL(wp) :: nsat,nlim                   !< foliar N that is light saturated / limited (mmol N)
  REAL(wp) :: n1,n2,n3                    !< N specific rates of electron-transport/Rubisco/PEPc limited PS
  REAL(wp) :: m1,m2,m3                    !< ci sensitivity of electron-transport/Rubisco/PEPc limited PS
  REAL(wp) :: msat                        !< N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
  REAL(wp) :: nco                         !< co-limitation factor
  REAL(wp) :: agsh,agsl                   !< gross photosynthesis on sunlit/shaded leaves (micro-mol / m2 / s)
  REAL(wp) :: fabs,fabsb                  !< fractional light absorption on chloroplast
  REAL(wp) :: ppm2Pa                      !< factor to convert gas concentrations from ppm to Pa
  REAL(wp) :: rtck                        !< dry static energy
  REAL(wp) :: an_cl_old,dan_cl            !< dummy variables needed for iteration loop of CI
  REAL(wp) :: agsf1,agsf2                 !< dummy variables for PS calculation
  INTEGER  :: iterate                     !< Iteration loop counter
  LOGICAL  :: last_call                   !< boolean to stop iteration at low conductance
  REAL(wp) :: ag_max_shaded_illumination  !< maximum photosynthesis, all light limited, for only shaded illumination (for ChlF)
  REAL(wp) :: ag_max_sunlit               !< maximum photosynthesis, all light limited, for sunlit leaves (for ChlF)
  REAL(wp) :: ag_max_shaded               !< maximum photosynthesis, all light limited, for shaded leaves (for ChlF)
  REAL(wp) :: frac_tot_sunlit             !< fraction of how much photosynthesis taking place from maximum for sunlit (for ChlF)
  REAL(wp) :: frac_tot_shaded             !< fraction of how much photosynthesis taking place from maximum for shaded (for ChlF)
  REAL(wp) :: chlfl_yield_sunlit          !< chlorophyll fluorescence yield sunlit leafs
  REAL(wp) :: chlfl_yield_shaded          !< chlorophyll fluorescence yield shaded leafs
  REAL(wp) :: chlfl_yield                 !< chlorophyll fluorescence yield sunlit/shaded leafs [unitless]
  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_photosynthesis'

  !> init out variables
  !>
  ag_cl           = 0.0_wp
  an_cl           = 0.0_wp
  maint_resp_cl   = 0.0_wp
  gs_cl           = 0.0_wp
  ci_cl           = 0.0_wp
  m_rub           = 0.0_wp
  m_et            = 0.0_wp
  m_pepc          = 0.0_wp
  chlfl_yield_cl  = 0.0_wp

  !> this routine is processed only for canopy layers with leafs
  !>
  IF(lai_cl > eps8) THEN

    !> static energy (R_gas: J mol-1 K-1 -> rtck: J mol-1)
    rtck = r_gas * t_air
    !> conversion of ppm to Pa
    ppm2Pa = press_srf * 1.e-06_wp
    !>  Net canopy layer respiration (umol[CO2]/m2/s)
    maint_resp_cl = fmaint_rate_base / 1000._wp * beta_soil_ps * temperature_response_respiration(t_air, temp_acclim) * &
                      leaf_nitrogen_cl

    !> calculate photosynthesis only if light is available at this canopy layer
    !>
    !> else: initialise output for no-light condition
    !>
    IF (ppfd_sunlit_cl + ppfd_shaded_cl <= eps8) THEN
      ag_cl    = 0.0_wp
      an_cl    = ag_cl - maint_resp_cl
      gs_cl    = MAX(gmin,g0 + g1 * an_cl * beta_air * beta_soil_gs / co2_mixing_ratio) * rtck / press_srf
      ci_cl    = co2_mixing_ratio * ppm2Pa
      m_et   = 0.0_wp
      m_rub  = 0.0_wp
      m_pepc = 0.0_wp
    ELSE
      !> Michaelis-Menten constant of Rubisco for CO2 (Pa)
      kc  = exp(E0kc - E1kc/rtck) * ppm2Pa
      !> Michaelis-Menten constant of Rubisco for O2 (kPa)
      ko  = exp(E0ko - E1ko/rtck) * ppm2Pa
      !> Combined Rubisco Km for CO2 (Pa); kc (1 + pO2/ko) is part of equation 5 in Kull and Kruijt 1998
      km  = kc * (1.0_wp + pO2/ko)
      !> Photorespiratory compensation point (Pa[CO2])
      pcp = exp(E0pcp - E1pcp/rtck) * ppm2Pa

      !> N specific rate of electron transport rate
      n1  = jmax2n * fn_et_cl * exp(-(((t_air - Tzero - t_jmax_opt) / t_jmax_omega) ** 2.0_wp)) * &
              beta_soil_ps * beta_sinklim_ps * beta_soa

      !> N specific rate of Rubisco turnover
      n2  = vcmax2n * fn_rub_cl * exp(E0v - E1v/rtck) * beta_soil_ps * beta_sinklim_ps * beta_soa

      !> Initial guess of leaf internal CO2 concentration
      SELECT CASE (ps_pathway)
        CASE (ic3phot)
          ci_cl = co2_mixing_ratio * CiCa_default_C3 * ppm2Pa
        CASE (ic4phot)
          ci_cl = co2_mixing_ratio * CiCa_default_C4 * ppm2Pa
      ENDSELECT

      !> iteration to find correct combination of ci_cl, Anl and gs_cl
      an_cl_old = 1000.0_wp
      dan_cl    = 1.0_wp
      iterate        = 1
      last_call = .FALSE.
      !> @par DO WHILE (iterate < ps_it_max .AND. dan_cl > eps4 )
      DO WHILE (iterate < ps_it_max .AND. dan_cl > eps4 )

        SELECT CASE (ps_pathway)
        CASE (ic3phot)
          !> ci sensitivity of Electron-transport limited photosynthesis
          m1 = ci_cl/(ci_cl + 2.0_wp * pcp)

          !> ci sensitivity of Rubisco-limited photosynthesis
          m2 = ci_cl/(ci_cl + km)

          !> N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
          msat = MIN(m1 * n1, m2 * n2)

          !> nlim factor (part of equation 12 (Kull and Kruijt 1998))
          nco = msat/(alpha_i * ka * m1)

          !> A couple of factors to calculate photosynthesis
          !!  parts of equations 14-16 of Kull and Kruijt 1998
          agsf1 = MAX(1.0_wp - pcp / ci_cl, 0.0_wp)
          agsf2 = agsf1 * alpha_i * m1
          agsf1 = agsf1 * msat

        CASE (ic4phot)
          !> N specific rate of PePc limited photosynthessis
          n3 = pepc2n * fn_pepc_cl * 2.0_wp ** ((t_air - Tzero - Tref_pepc )/ Tbase_pepc)

          !> ci sensitivity of Electron-transport limited photosynthesis
          m1 = ci_max/(ci_max + 2.0_wp * pcp)

          !> ci sensitivity of Rubisco-limited photosynthesis
          m2 = ci_max/(ci_max + km)

          !> ci sensitivity of PePc-limited photosynthesis
          m3 = ci_cl / press_srf

          !> N-specific light saturated carboxylation (umol[CO2]/mmol[N]/s)
          msat = MIN(MIN(m1 * n1, m2 * n2), m3 * n3)

          !> nlim factor (part of equation 12 (Kull and Kruijt 1998))
          nco = msat/(alpha_i * ka * m1)

          !> A couple of factors to calculate photosynthesis
          !!  parts of equations 14-16 of Kull and Kruijt 1998
          agsf1 = MAX(1.0_wp - pcp / ci_max, 0.0_wp)
          agsf2 = agsf1 * alpha_i * m1
          agsf1 = agsf1 * msat

        END SELECT

        ! output factors for optimal fraction calculation
        m_et  = m1 * n1 / fn_et_cl
        m_rub = m2 * n2 / fn_rub_cl
        IF (fn_pepc_cl>eps8) THEN
          m_pepc = m3 * n3 / fn_pepc_cl
        ELSE
          m_pepc = 0._wp
        ENDIF

        !> fabsb: fraction of possible PAR absorption in present light conditions
        fabsb = exp(-ka * chl2n * fn_chl_cl * leaf_nitrogen_cl)

        !> @par IF: Case of no direct radiation
        IF (.NOT. fleaf_sunlit_cl > eps8) THEN
          !>  Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in sunlit foliage (mmol[N]/m2))
          nlim = -log(nco/(ppfd_shaded_cl * chl2n * fn_chl_cl + eps8))/(ka * chl2n * fn_chl_cl)
          nsat = MAX(MIN(nlim, leaf_nitrogen_cl), 0.0_wp)

          !> Relative PAR absorption in light-limited region (fraction)
          fabs = exp(-ka * chl2n * fn_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis (umol/m2/s)
          !!  sum of photosynthesis in colimited, light saturated PS and PS in the
          !!  light-limited region (assumed to scale linearily with light)
          ag_cl = agsf1 * nsat + agsf2 * ppfd_shaded_cl * fabs

          ! calc maximum photosynthesis, all light limited, for only shaded illumination
          ag_max_shaded_illumination = agsf2 * ppfd_shaded_cl * (1.0_wp - fabsb)
          ! calc fraction of how much photosynthesis taking place from photosynthesis maximum for shaded
          IF (ag_cl < eps8) THEN
            frac_tot_shaded = 0.0_wp
          ELSE
            frac_tot_shaded = ag_cl / ag_max_shaded_illumination
          END IF

        !> @par ELSE: Case of direct and diffuse radiation
        ELSE
          !>  Theoretical cumulative N at co-limitation by Rubisco and
           !!  light absorption in sunlit foliage (mmol[N]/m2))
          nlim = -log(nco/(ppfd_sunlit_cl * chl2n * fn_chl_cl + eps8))/(ka * chl2n * fn_chl_cl)
          nsat = MAX(MIN(nlim, leaf_nitrogen_cl), 0.0_wp)

          !>  Relative PAR absorption in light-limited region of sunlit foliage (fraction)
          fabs = exp(-ka * chl2n * fn_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis in sunlit foliage (umol/m2/s)
          agsl = agsf1 * nsat + agsf2 * ppfd_sunlit_cl * fabs

          ! Calculation for the chlorophyll fluorescence
          ag_max_sunlit =  agsf2 * ppfd_sunlit_cl * (1.0_wp - fabsb )
          ! sunlit part
          IF (agsl < eps8) THEN
            frac_tot_sunlit = 0.0_wp
          ELSE
            frac_tot_sunlit = agsl / ag_max_sunlit
          END IF

          !> Theoretical cumulative N at co-limitation by Rubisco and
          !!  light absorption in shaded foliage (mmol[N]/m2))
          nlim = -log(nco/(ppfd_shaded_cl * chl2n * fn_chl_cl + eps8))/(ka * chl2n * fn_chl_cl)
          nsat = MAX(MIN(nlim, leaf_nitrogen_cl), 0.0_wp)

          !> Relative PAR absorption in light-limited region of shaded foliage (fraction)
          fabs = exp(-ka * chl2n * fn_chl_cl * nsat) - fabsb

          !> Rate of gross photosynthesis in shaded foliage (umol/m2/s)
          agsh = agsf1 * nsat + agsf2 * ppfd_shaded_cl* fabs

          ! calc maximum photosynthesis, all light limited, for shaded leaves
          ag_max_shaded = agsf2 * ppfd_shaded_cl* (1.0_wp - fabsb)
          ! calc fraction of how much photosynthesis taking place from photosynthesis maximum for shaded
          IF (ag_max_shaded < eps8) THEN
            frac_tot_shaded = 0.0_wp
          ELSE
            frac_tot_shaded = agsh / ag_max_shaded
          END IF

          !> Mean gross photosynthesis over sunlit and shaded fractions (umol/m2/s)
          ag_cl = fleaf_sunlit_cl * agsl + (1.0_wp - fleaf_sunlit_cl) * agsh
        END IF

        !> Net layer photosynthesis (umol[CO2]/m2/s)
        an_cl = ag_cl - maint_resp_cl

        !> calculate canopy conductance (m/s)
        SELECT CASE(canopy_cond_scheme)
          CASE ("medlyn")
            !following Medlyn et al. 2011 [default]
            gs_cl = MAX(g0, 1.6_wp * ( 1._wp + g1 * beta_air * beta_soil_gs ) * an_cl / ( co2_mixing_ratio )) * rtck / press_srf
          CASE ("ballberry")
            ! following Ball & Berry
            gs_cl = MAX(gmin, g0 + g1 * an_cl * beta_air * beta_soil_gs / ( co2_mixing_ratio )) * rtck / press_srf
        END SELECT

        !> New internal CO2 from balance of supply and demand
        !!
        !! ci = cs - An * r_tot (Friend 2005 eq. 7)
        !! Friend 2005 eq. 8: r_tot = ( Dh2oco2_air / gs_cl + Dh2oco2_turb / (speed * q_cdrag))
        !! Corrected to Bonan et al. 2008, Ecol. Clim: accounting for the molecular diffusion differences
        !! at the stomatal interface (1.6) and in the boundary layer, where turbulence plays a role (1.37)
        ci_cl = co2_mixing_ratio * ppm2Pa - an_cl * rtck * 1e-06_wp * (Dwv2co2_air/gs_cl + Dwv2co2_turb/aerodyn_cond)
        ci_cl = MAX (ci_cl, pcp)

        !> change in An since last iteration and save new an_cl
        dan_cl = abs(an_cl-an_cl_old)
        an_cl_old = an_cl

        !> break-off point for numerical stability at low conductance
        IF(last_call) THEN
          iterate    = ps_it_max + 1   ! ps_it_max is of TYPE REAL; iterate is of TYPE INTEGER
          ci_cl = pcp
        ENDIF
!        IF( gs_cl <= g0 * rtck / press_srf .OR. ci_cl <= pcp ) THEN
!          last_call = .TRUE.
!          ci_cl     = pcp
!        ENDIF
        iterate =  iterate + 1
      END DO  ! closing: [DO WHILE (iterate < ps_it_max .AND. dan_cl > eps4 ) ]
      !> END DO WHILE
    ENDIF     ! light / no light condition
    ! for output purposes, convert ci to ppm
    ci_cl = ci_cl / ppm2Pa

    ! ChlF (chlorophyll fluorescence yield) calculation
    ! if light levels above eps8
    IF (ppfd_sunlit_cl + ppfd_shaded_cl > eps8) THEN
      ! sunlit area fluorescence
      IF (fleaf_sunlit_cl > eps8) THEN
        CALL calc_chlorophyll_fluorescence_yield(frac_tot_sunlit, t_air, chlfl_yield)
        chlfl_yield_sunlit = chlfl_yield
      ELSE
        chlfl_yield_sunlit = 0.0_wp
      ENDIF
      ! shaded area fluorescence
      CALL calc_chlorophyll_fluorescence_yield(frac_tot_shaded, t_air, chlfl_yield)
      chlfl_yield_shaded = chlfl_yield
      ! Add in the sunlit and shaded parts of the fluorescence yield
      chlfl_yield_cl = fleaf_sunlit_cl * chlfl_yield_sunlit + (1.0_wp - fleaf_sunlit_cl) * chlfl_yield_shaded
    ELSE
      chlfl_yield_cl = 0.0_wp
    ENDIF

  ! lai_cl < eps8
  ELSE
    !> set all output var to zero
    ag_cl          = 0.0_wp
    an_cl          = 0.0_wp
    maint_resp_cl  = 0.0_wp
    gs_cl          = 0.0_wp
    ci_cl          = 0.0_wp
    m_rub          = 0.0_wp
    m_et           = 0.0_wp
    m_pepc         = 0.0_wp

  END IF  ! end of IF(lai_cl > eps8)

  END SUBROUTINE calc_photosynthesis

  ! ======================================================================================================= !
  !>calculate atmospheric drought stress factor
  !>
  !>  Literature:
  !>    Medlyn et al. 2011: beta = MIN(5._wp, 1._wp / SQRT(vpd + eps8))
  !>    Ball & Berry: beta = q_air / q_sat
  !>
  PURE ELEMENTAL FUNCTION calc_atm_moisture_stress(canopy_cond_scheme, t_air, q_air, press_srf) RESULT(beta)

    USE mo_atmland_constants,  ONLY: eps_vpd
    USE mo_atmland_util,       ONLY: calc_spec_humidity_sat
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), INTENT(IN)  :: canopy_cond_scheme     !< canopy_conductance_scheme: medlyn / ballberry
    REAL(wp),INTENT(in)           :: t_air                  !< air temperature in K
    REAL(wp),INTENT(in)           :: q_air                  !< specific air humidity (g g-1)
    REAL(wp),INTENT(in)           :: press_srf              !< air pressure (Pa)
    REAL(wp)                      :: beta                   !< atmospheric moisture scalar for use in conductance calculation
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: q_sat                         !< saturating specific humidity
    REAL(wp) :: vpd                           !< vapour pressure deficit
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_atm_moisture_stress'

    !>1.0 saturating specific humidity [kg kg-1]
    !>
    q_sat = calc_spec_humidity_sat(t_air, press_srf)

    !>  1.1 vapour pressure deficit [kPa]
    !>
    vpd = MAX(q_sat - q_air, 0.0_wp) * press_srf / eps_vpd / 1000._wp

    !>  1.2 scalar for atmospheric humidity as relative humidity (fractions)
    !>
    SELECT CASE(canopy_cond_scheme)
      CASE ("medlyn")
        !following Medlyn et al. 2011 [default]
        beta = MIN(5._wp, 1._wp / SQRT(vpd + eps8))
      CASE ("ballberry")
        ! following Ball & Berry
        beta = q_air / q_sat
    END SELECT
  END FUNCTION calc_atm_moisture_stress

  ! ======================================================================================================= !
  !> calculate chlorophyll fluorescence - diagnostic output
  !>
  !> does not affect photsynthesis calculation
  !>
  !> re-implementation of the model as published in:
  !>    van der Tol et al, JGR, 2014. SCOPE v1.70
  !>    https://doi.org/10.1002/2014JG002713
  !>
  !> source of parameter values: van der Tol et al, JGR, 2014;
  !>   or papers cited within, e.g., Flexas 2002 & Weis and Berry 1987
  !>
  ELEMENTAL SUBROUTINE calc_chlorophyll_fluorescence_yield( &
    & frac_tot, &
    & t_air, &
    & chlfl_yield)

    USE mo_jsb_physical_constants, ONLY: Tzero
    USE mo_q_assimi_parameters,    ONLY: heat_dissipation_enrg_dep_param_1, &
      &                                  heat_dissipation_enrg_dep_param_2, &
      &                                  heat_dissipation_enrg_dep_param_3, &
      &                                  k_fluorescence_rate, &
      &                                  thermal_dissipation_enrg_constitutive_max_rate, &
      &                                  thermal_dissipation_enrg_constitutive_param_1, &
      &                                  thermal_dissipation_enrg_constitutive_param_2, &
      &                                  k_photochemistry_rate

    REAL(wp),     INTENT(IN)  :: frac_tot        !< fraction of actual photosynthesis to max
    REAL(wp),     INTENT(IN)  :: t_air           !< air temperature [K]
    REAL(wp),     INTENT(OUT) :: chlfl_yield     !< chlorophyll fluorescence yield sunlit/shaded leafs [unitless]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp) :: thermal_dissipation_enrg_constitutive_rate  !< rate constant for energy constitutive thermal dissipation
    REAL(wp) :: heat_dissipation_enrg_dep_rate  !< rate constant for energy dependent heat dissipation
    REAL(wp) :: photochemical_yield_rel         !< x in eq. 14 (van der Tol et al., 2014)
    REAL(wp) :: frac_dark_photochemistry        !< dark photochemistry fraction (Genty et al., 1989)
    REAL(wp) :: photochemical_yield             !< photochemical yield
    REAL(wp) :: light_saturation_rel            !< relative light saturation
    REAL(wp) :: fluorescence_yield_fo_dark      !< dark adapted fluorescence yield Fo
    REAL(wp) :: fluorescence_yield_fo_light     !< light adapted fluorescence yield Fo
    REAL(wp) :: fluorescence_yield_fm_dark      !< dark adapted fluorescence yield Fm
    REAL(wp) :: fluorescence_yield_fm_light     !< light adapted fluorescence yield Fm
    REAL(wp) :: quenching_phot                  !< photochemical quenching
    REAL(wp) :: quenching_nonphot               !< non-photochemical quenching
    REAL(wp) :: chlfl_yield_scaled              !< scaled chlorophyll fluorescence yield
    REAL(wp) :: photochemistry_rate             !< actual rate for photochemistry
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':chlorophyll fluorescence yield'
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.9.1 init local variables
    !>
    !> 0.9.2 parameter values are not PFT-specific (van der Tol et al. 2014)
    !>
    thermal_dissipation_enrg_constitutive_rate  = &
      &   MAX(thermal_dissipation_enrg_constitutive_max_rate, thermal_dissipation_enrg_constitutive_param_1 &
      &   * (t_air - Tzero) + thermal_dissipation_enrg_constitutive_param_2)
    frac_dark_photochemistry = &
      &   k_photochemistry_rate / (k_fluorescence_rate + thermal_dissipation_enrg_constitutive_rate + k_photochemistry_rate)
    photochemical_yield     = &
      &   MIN(frac_dark_photochemistry, frac_dark_photochemistry * frac_tot)
    photochemical_yield_rel = &
      &   MAX(0.0_wp, 1.0_wp - photochemical_yield / frac_dark_photochemistry)

    !> 1.0 Calculation of the relative light saturation
    !>
    IF (photochemical_yield_rel > -eps8 .AND. photochemical_yield_rel < eps8) THEN
      light_saturation_rel = 0.0_wp
    ELSE
      light_saturation_rel = EXP(LOG(photochemical_yield_rel) * heat_dissipation_enrg_dep_param_2)
    ENDIF
    heat_dissipation_enrg_dep_rate  = heat_dissipation_enrg_dep_param_1 * (1.0_wp + heat_dissipation_enrg_dep_param_3) &
      &                               * light_saturation_rel / (heat_dissipation_enrg_dep_param_3 + light_saturation_rel)
    fluorescence_yield_fo_dark      = k_fluorescence_rate / (k_fluorescence_rate + k_photochemistry_rate &
      &                               + thermal_dissipation_enrg_constitutive_rate)
    fluorescence_yield_fo_light     = k_fluorescence_rate / (k_fluorescence_rate + k_photochemistry_rate &
      &                               + thermal_dissipation_enrg_constitutive_rate + heat_dissipation_enrg_dep_rate)
    fluorescence_yield_fm_light     = k_fluorescence_rate / (k_fluorescence_rate + thermal_dissipation_enrg_constitutive_rate &
      &                               + heat_dissipation_enrg_dep_rate)
    fluorescence_yield_fm_dark      = k_fluorescence_rate / (k_fluorescence_rate + thermal_dissipation_enrg_constitutive_rate)

    !> 2.0 calc chlorophyll fluorescence yield
    !>
    chlfl_yield  = fluorescence_yield_fm_light * (1.0_wp - photochemical_yield)

    ! @NOTE the below calculated variables are not yet used
    !       but useful for future development (of this routine)
    ! !> 3.0 additional calculations
    ! !>
    ! chlfl_yield_scaled = chlfl_yield / fluorescence_yield_fo_dark
    ! quenching_phot     = 1.0_wp - (chlfl_yield - fluorescence_yield_fo_light) &
    !   &                  / (fluorescence_yield_fm_light - fluorescence_yield_fo_light)
    ! quenching_nonphot  = 1.0_wp - (fluorescence_yield_fm_light - fluorescence_yield_fo_light) &
    !   &                  / (fluorescence_yield_fm_dark - fluorescence_yield_fo_dark )
    ! ! calc actual rate for photochemistry
    ! photochemistry_rate = photochemical_yield / chlfl_yield * k_fluorescence_rate
  END SUBROUTINE calc_chlorophyll_fluorescence_yield

  ! ======================================================================================================= !
  !>calculate soil moisture stress on photosynthesis or conductance given current soil moisture
  !>
  PURE ELEMENTAL FUNCTION calc_soil_moisture_stress(phi, phi_leaf_min) RESULT(beta)

    REAL(wp),INTENT(in) :: phi               !< soil plant available water potential
    REAL(wp),INTENT(in) :: phi_leaf_min      !< minimum leaf water potential from lctlib
    REAL(wp)            :: beta              !< veg\%beta_soil_ps / veg\%beta_soil_gs
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_soil_moisture_stress'!

    !>1.0 Moisture stress as a function of soil water potential and minimum leaf water potential
    !>
    beta = MAX(0.0_wp, 1._wp - phi / phi_leaf_min)

  END FUNCTION calc_soil_moisture_stress

  ! ======================================================================================================= !
  !>Calculate spring recovery period for evergreen conifers
  !>
  !>  Based on Maekelae et al. (2004) https://doi.org/10.1093/treephys/24.4.369
  !>  returns beta_soa with a value [0.01,1.0]
  !>
  PURE ELEMENTAL FUNCTION calc_beta_soa(soa_tsoa_mavg,phenology_type) RESULT(beta_soa)

    USE mo_q_pheno_constants,          ONLY: ievergreen
    USE mo_q_assimi_parameters,        ONLY: soa_b, soa_t_s
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(in)  :: soa_tsoa_mavg    !< state of acclimation
    INTEGER,  INTENT(in)  :: phenology_type   !< phenology type of the plant
    REAL(wp)              :: beta_soa         !< beta value of soa
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_beta_soa'

    !>1.0 calc beta_soa
    !>
    IF (phenology_type == ievergreen) THEN
      beta_soa = 1.0_wp / (1.0_wp + EXP(soa_b * ( soa_tsoa_mavg - soa_t_s ) ))
      beta_soa = MIN ( MAX(beta_soa,0.01_wp), 1.0_wp)
    ELSE
      beta_soa = 1.0_wp
    END IF

  END FUNCTION calc_beta_soa

  ! ======================================================================================================= !
  !>Calculate discrimination of C13 and C14 due to photosynthesis
  !>
  !>  Follows the orignial ideas of Farquhar, as described in
  !>  Drake (2014), Radiocarbon, 56, 29-38
  !>
  !>  Input: leaf-level Ci and Ca, photosynthetic pathway (C3/C4) and name of the isotope
  !>
  !>  Output: Discrimination (per mill)
  !>
  PURE ELEMENTAL FUNCTION discrimination_ps(ci, ca, isotope_name, ps_pathway) RESULT(discrimination)

    USE mo_q_assimi_parameters,    ONLY: discr_ps_a_C13, discr_ps_b_C13, discr_ps_c_C13, &
                                         discr_ps_a_C14, discr_ps_b_C14, discr_ps_c_C14, &
                                         discr_ps_phi
    USE mo_q_assimi_constants,     ONLY: ic3phot, ic4phot
    !------------------------------------------------------------------------------------------------------ !
    REAL(wp),     INTENT(in) :: ci                     !< leaf internal CO2 concentration (ppm)
    REAL(wp),     INTENT(in) :: ca                     !< atmospheric CO2 concentration (ppm)
    CHARACTER(3), INTENT(in) :: isotope_name           !< isotope C13 or C14
    INTEGER,      INTENT(in) :: ps_pathway             !< photosynthetic pathway (C3 or C4) from lctlib
    REAL(wp)                 :: discrimination         !< isotopic discrimination (per mill)
    !------------------------------------------------------------------------------------------------------ !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':discrimination_ps'

    ! SELECT PS pathway C3 C4
    SELECT CASE(ps_pathway)
    CASE (ic3phot)
      ! SELECT isotope C13 C14
      SELECT CASE(isotope_name)
        CASE ('C13')
          discrimination = discr_ps_a_C13 + (discr_ps_b_C13 - discr_ps_a_C13) * ci / ca
        CASE ('C14')
          discrimination = discr_ps_a_C14 + (discr_ps_b_C14 - discr_ps_a_C14) * ci / ca
      END SELECT
    CASE (ic4phot)
      ! SELECT isotope C13 C14
      SELECT CASE(isotope_name)
        CASE ('C13')
          discrimination = discr_ps_a_C13 + (discr_ps_c_C13 + discr_ps_phi * discr_ps_b_C13 - discr_ps_a_C13) * ci / ca
        CASE ('C14')
          discrimination = discr_ps_a_C14 + (discr_ps_c_C14 + discr_ps_phi * discr_ps_b_C14 - discr_ps_a_C14) * ci / ca
      END SELECT
    END SELECT

  END FUNCTION discrimination_ps

#endif
END MODULE mo_q_assimi_process
