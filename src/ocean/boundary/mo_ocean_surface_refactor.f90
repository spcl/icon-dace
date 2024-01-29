! Provide an implementation of the ocean surface module.
!
! Provide an implementation of the parameters used for surface forcing
! of the hydrostatic ocean model.
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

MODULE mo_ocean_surface_refactor
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2007
!
!-------------------------------------------------------------------------
!
  USE mo_kind,                ONLY: wp
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: dtime
  USE mo_dynamics_config,     ONLY: nold
  USE mo_exception,           ONLY: finish, message
  USE mo_util_dbg_prnt,       ONLY: dbg_print

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range

  USE mo_ocean_nml,           ONLY: iforc_oce, no_tracer, type_surfRelax_Temp, type_surfRelax_Salt, &
    &  No_Forcing, Analytical_Forcing, OMIP_FluxFromFile, Coupled_FluxFromAtmo,                     &
    &  i_sea_ice, zero_freshwater_flux, atmos_flux_analytical_type, atmos_precip_const, &  ! atmos_evap_constant
    &  limit_elevation, lhamocc, lswr_jerlov, lhamocc, lfb_bgc_oce, lcheck_salt_content, &
    &  lfix_salt_content, ice_flux_type, heatflux_forcing_on_sst, &
    &  lfwflux_enters_with_sst, vert_cor_type

  USE mo_sea_ice_nml,        ONLY: sice

  USE mo_ocean_nml,           ONLY: atmos_flux_analytical_type, relax_analytical_type, &
    &  n_zlev, para_surfRelax_Salt, para_surfRelax_Temp, atmos_precip_const, &  ! atmos_evap_constant
    &  atmos_SWnet_const, atmos_LWnet_const, atmos_lat_const, atmos_sens_const, &
    &  atmos_SWnetw_const, atmos_LWnetw_const, atmos_latw_const, atmos_sensw_const, &
    &  relax_width, forcing_HeatFlux_amplitude, forcing_HeatFlux_base,              &
    &  basin_center_lat, basin_center_lon, basin_width_deg, basin_height_deg

  USE mo_math_constants,      ONLY: pi, deg2rad, rad2deg
  USE mo_physical_constants,  ONLY: rhoi, rhos, rho_ref, alv, tmelt, tf, clw, stbo, zemiss_def
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, MIN_DOLIC

  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_sea_ice_types,       ONLY: t_sea_ice, t_atmos_fluxes
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_operator_ocean_coeff_3d,ONLY: t_operator_coeff

  USE mtime,                  ONLY: datetime, getNoOfSecondsElapsedInDayDateTime
  USE mo_ice_interface,       ONLY: ice_fast_interface, ice_thermodynamics, ice_dynamics
  USE mo_ocean_bulk_forcing,  ONLY: update_surface_relaxation, apply_surface_relaxation, &
                                &   update_flux_fromFile, calc_omip_budgets_ice, calc_omip_budgets_oce, &
                                &   update_ocean_surface_stress, balance_elevation, &
                                &   balance_elevation_zstar

  USE mo_ocean_diagnostics, ONLY : diag_heat_salt_tendency
  USE mo_name_list_output_init, ONLY: isRegistered

  USE mo_swr_absorption
  USE mo_ocean_check_total_content,       ONLY: check_total_salt_content, check_total_si_volume

  USE mo_mpi, only: get_my_mpi_work_id
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE

  ! public interface
  PUBLIC  :: update_ocean_surface_refactor
  ! private routine

  PUBLIC  :: apply_surface_fluxes_slo
  PUBLIC  :: update_atmos_fluxes
  PRIVATE :: update_atmos_fluxes_analytical

  CHARACTER(len=12)           :: str_module    = 'OceanSurfaceRefactor'  ! Output of module for 1 line debug
  INTEGER                     :: idt_src       = 1               ! Level of detail for 1 line debug

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !
  SUBROUTINE apply_surface_fluxes_slo(p_patch_3D, p_os, p_ice, p_oce_sfc, lacc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    LOGICAL, INTENT(IN), OPTIONAL               :: lacc
    !
    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:apply_surface_fluxes_slo'
    INTEGER               :: jc, jb, trac_no
    INTEGER               :: i_startidx_c, i_endidx_c
    REAL(wp)              :: sss_inter, zUnderIceOld, zUnderIceIni, zUnderIceArt

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer
    REAL(wp)  :: zunderice_ini
    LOGICAL   :: lzacc

    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    CALL dbg_print('UpdSfc: h-old',p_os%p_prog(nold(1))%h,         str_module, 1, in_subset=p_patch%cells%owned)

    !  ******  (Thermodynamic Eq. 1)  ******
    ! Apply net surface heat flux to ocean surface (new p_oce_flx%SST)
    IF (no_tracer > 0) THEN
      IF ( heatflux_forcing_on_sst ) THEN
        ! sst-change in surface module after sea-ice thermodynamics using HeatFlux_Total and old freeboard
        ! freeboard before sea ice model (used for thermal boundary condition (Eq.1))
        ! by construction, is stored in p_oce_sfc%cellThicknessUnderIce
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(heatflux_surface_layer) ASYNC(1) IF(lzacc)
          DO jc = i_startidx_c, i_endidx_c
            IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

              p_os%p_diag%heatflux_rainevaprunoff(jc,jb) =                                   &
                  ( p_oce_sfc%FrshFlux_TotalOcean(jc,jb)+ p_oce_sfc%FrshFlux_Runoff(jc,jb) )    &
                    * ( p_oce_sfc%sst(jc,jb) - tf ) * clw * rho_ref

              ! fw flux enters with "SST" so the  flux is implicitly added

              ! substract the fraction of heatflux used for subsurface heating
              heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)
              p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
                &                    heatflux_surface_layer*dtime/(clw*rho_ref*p_oce_sfc%cellThicknessUnderIce(jc,jb))

            ENDIF
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO
        !$ACC WAIT(1)
      ENDIF
    END IF

    ! apply volume flux to surface elevation
    !  - add to h_old before explicit term
    !  - change in salt concentration applied here
    !    i.e. for salinity relaxation only, no volume flux is applied
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(zUnderIce_ini) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c
        IF (p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb) > 0) THEN

          !!  Provide total ocean forcing:
          !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
          !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
          !    - total freshwater volume forcing
          p_oce_sfc%FrshFlux_VolumeTotal(jc,jb) = p_oce_sfc%FrshFlux_Runoff    (jc,jb) &
            &                                   + p_oce_sfc%FrshFlux_VolumeIce (jc,jb) &
            &                                   + p_oce_sfc%FrshFlux_TotalOcean(jc,jb)
          ! provide total salinity forcing flux for diagnostics only
          p_oce_sfc%FrshFlux_TotalSalt(jc,jb)   = p_oce_sfc%FrshFlux_Runoff    (jc,jb) &
            &                                   + p_oce_sfc%FrshFlux_TotalIce  (jc,jb) &
            &                                   + p_oce_sfc%FrshFlux_TotalOcean(jc,jb)

          !******  (Thermodynamic Eq. 2)  ******
          !! Calculate the new freeboard caused by changes in ice thermodynamics
          !!  zUnderIce = z_surf + h_old - (z_draft - z_snowfall)
          !  #slo# 2015-01: totalsnowfall is needed for correct salt update (in surface module)
          !                 since draft was increased by snowfall but water below ice is not affected by snowfall
          !                 snow to ice conversion does not effect draft
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + &
          &                      p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave(jc,jb) + p_ice%totalsnowfall(jc,jb)

          !******  (Thermodynamic Eq. 3)  ******
          !! First, calculate internal salinity change caused by melting of snow and melt or growth of ice:
          !!   SSS_new * zUnderIce = SSS_old * zUnderIceArt
          !!   artificial freeboard zUnderIceArt is used for internal Salinity change only:
          !!   - melt/growth of ice and snow to ice conversion imply a reduced water flux compared to saltfree water
          !!   - reduced water flux is calculated in FrshFlux_TotalIce by the term  (1-Sice/SSS)
          !!   - respective zUnderIceArt for calculating salt change is derived from these fluxes
          !!     which are calculated in sea ice thermodynamics (upper_ocean_TS)
          !    - for i_sea_ice=0 it is FrshFlux_TotalIce=0 and no change here
          zUnderIceArt = p_ice%zUnderIce(jc,jb) - p_oce_sfc%FrshFlux_TotalIce(jc,jb)*dtime
          sss_inter    = p_oce_sfc%sss(jc,jb) * zUnderIceArt / p_ice%zUnderIce(jc,jb)

              !******  (Thermodynamic Eq. 4)  ******
          !! Next, calculate salinity change caused by rain and runoff without snowfall by adding their freshwater to zUnderIce
          zUnderIceOld           = p_ice%zUnderIce(jc,jb)
          p_ice%zUnderIce(jc,jb) = zUnderIceOld + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb) * dtime
          p_oce_sfc%SSS(jc,jb)   = sss_inter * zUnderIceOld / p_ice%zUnderIce(jc,jb)

         zUnderIce_ini=  p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
              &                    + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave_old(jc,jb)
          !******  (Thermodynamic Eq. 5)  ******
          !! Finally, let sea-level change from P-E+RO plus snow fall on ice, net total volume forcing to ocean surface
          p_os%p_prog(nold(1))%h(jc,jb) = p_os%p_prog(nold(1))%h(jc,jb)               &
            &                           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
            &                           + p_ice%totalsnowfall(jc,jb)

          !! update zunderice
          p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) + p_os%p_prog(nold(1))%h(jc,jb) &
            &                    - p_ice%draftave(jc,jb)

          p_oce_sfc%top_dilution_coeff(jc,jb) = zUnderIce_ini / p_ice%zUnderIce(jc,jb)
          
          !! set correct cell thickness under ice
          p_oce_sfc%cellThicknessUnderIce(jc,jb) = p_ice%zUnderIce(jc,jb)
        ENDIF  !  dolic>0
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    CALL dbg_print('UpdSfc: oce_sfc%HFTot ', p_oce_sfc%HeatFlux_Total,       str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%VolTot', p_oce_sfc%FrshFlux_VolumeTotal, str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: oce_sfc%TotIce', p_oce_sfc%FrshFlux_TotalIce,    str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: ice%totalsnowf', p_ice%totalsnowfall,            str_module, 4, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfc: zUnderIce   ',   p_ice%zUnderIce,                str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS,                  str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEnd: h-old+fwfVol',p_os%p_prog(nold(1))%h,         str_module, 2, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE apply_surface_fluxes_slo

  SUBROUTINE close_salt_budget(step, p_patch_3D, p_os, p_ice, p_oce_sfc)

    TYPE(t_patch_3D ),TARGET,   INTENT(IN)      :: p_patch_3D
    TYPE(t_hydro_ocean_state),  INTENT(INOUT)   :: p_os
    TYPE(t_sea_ice),            INTENT(IN)   :: p_ice
    TYPE(t_ocean_surface),      INTENT(INOUT)   :: p_oce_sfc

    REAL(wp)     :: zUnderIce_old
    !
    ! local variables
    TYPE(t_patch), POINTER          :: p_patch
    TYPE(t_subset_range), POINTER   :: all_cells
    INTEGER                         :: step, jc, jb, i_startidx_c, i_endidx_c



    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------


    IF (no_tracer > 0) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(i_startidx_c, i_endidx_c, jc,zunderice_old) SCHEDULE(dynamic)
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
            IF (step .EQ. 1 ) THEN

              zUnderIce_old = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
                   + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave_old(jc,jb)

              ! new salt from ice after the ice advection ;  liquid water part  before ice advection
              p_oce_sfc%surface_salt_content(jc,jb) =sice*rhoi/rho_ref*SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb)) &
                   +p_os%p_prog(nold(1))%tracer(jc,1,jb,2)*zUnderIce_old

            ENDIF

            IF (step .EQ. 2 ) THEN

              p_ice%draft(jc,:,jb) = (rhos * p_ice%hs(jc,:,jb) + rhoi * p_ice%hi(jc,:,jb))/rho_ref
              p_ice%draftave(jc,jb) = SUM(p_ice%draft(jc,:,jb) * p_ice%conc(jc,:,jb))
              p_ice%zUnderIce(jc,jb) = p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,1,jb) &
                   + p_os%p_prog(nold(1))%h(jc,jb) - p_ice%draftave(jc,jb)

              p_oce_sfc%sss(jc,jb)   = ( p_oce_sfc%surface_salt_content(jc,jb) - sice * rhoi/rho_ref * &
                 & SUM(p_ice%hi(jc,:,jb)*p_ice%conc(jc,:,jb))) &
                 &  /p_ice%zUnderIce(jc,jb)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO
    END IF


  END SUBROUTINE close_salt_budget


  !-------------------------------------------------------------------------
  !>
  !! Update ocean surface by applying flux forcing for hydrostatic ocean
  !!
  !! This function changes:
  !! p_ice      thermodynamical and dynamical fields of sea ice
  !! p_oce_sfc  surface fluxes and stress, passed to the ocean
  !! p_os       SSH, SST, SSS and HAMMOC tracers (dilution)
  !!
  !  Adapted for zstar
  !
  SUBROUTINE update_ocean_surface_refactor(p_patch_3D, p_os, p_as, p_ice, &
      & atmos_fluxes, p_oce_sfc, this_datetime, p_op_coeff, eta_c, stretch_c, &
      & lacc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    TYPE(datetime), POINTER                     :: this_datetime
    TYPE(t_operator_coeff),   INTENT(IN)        :: p_op_coeff
    REAL(wp), INTENT(INOUT),OPTIONAL :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht
    REAL(wp), INTENT(IN   ),OPTIONAL :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL,  INTENT(IN   ),OPTIONAL :: lacc
    !
    ! local variables
    TYPE(t_patch), POINTER                      :: p_patch
    TYPE(t_subset_range), POINTER               :: all_cells
    INTEGER                                     :: trac_no
    REAL(wp)                                    :: dsec
    LOGICAL                                     :: lzacc
    INTEGER                                     :: jc, jb, i_startidx_c, i_endidx_c

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:update_ocean_surface'

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    ! subset range pointer
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c
        p_oce_sfc%top_dilution_coeff(jc,jb) = 1.0_wp ! Initilize dilution factor to one
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)

    IF (no_tracer>=1) p_oce_sfc%sst => p_os%p_prog(nold(1))%tracer(:,1,:,1)
    IF (no_tracer>=2) p_oce_sfc%sss => p_os%p_prog(nold(1))%tracer(:,1,:,2)

    IF (iforc_oce == No_Forcing) RETURN  !  forcing for ocean not defined

    ! save values of ice, snow and temperature for the tendendy and hfbasin diagnostic
    IF ( isRegistered('delta_ice') .OR. isRegistered('delta_snow') .OR. &
           isRegistered('delta_thetao') .OR. &
           isRegistered('delta_so') .OR. &
           isRegistered('global_sltbasin') .OR. isRegistered('atlant_sltbasin') .OR. &
           isRegistered('pacind_sltbasin') .OR. &
           isRegistered('global_hfbasin') .OR. isRegistered('atlant_hfbasin') .OR. &
           isRegistered('pacind_hfbasin') ) THEN

      CALL diag_heat_salt_tendency(p_patch_3d, 1, p_ice,     &
         p_os%p_prog(nold(1))%tracer(:,:,:,1),               &
         p_os%p_prog(nold(1))%tracer(:,:,:,2),               &
         p_os%p_diag%delta_ice,                              &
         p_os%p_diag%delta_snow, p_os%p_diag%delta_thetao,   &
         p_os%p_diag%delta_so, lacc=lzacc)

    END IF

    !---------------------------------------------------------------------
    ! (1) Apply relaxation to surface temperature and salinity
    !---------------------------------------------------------------------
    IF (type_surfRelax_Temp >= 1) THEN
      trac_no = 1   !  tracer no 1: temperature

      IF (vert_cor_type .EQ. 0 ) THEN
        CALL update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no, lacc=lzacc)
      ELSE
        CALL update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no, stretch_c, lacc=lzacc)
      ENDIF

      !  apply restoring to surface temperature directly
      CALL apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, trac_no, lacc=lzacc)

    END IF

    IF (type_surfRelax_Salt >= 1 .AND. no_tracer >1) THEN
      trac_no = 2   !  tracer no 2: salinity

      IF (vert_cor_type .EQ. 0) THEN
        CALL update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no, lacc=lzacc)
      ELSE
        CALL update_surface_relaxation(p_patch_3D, p_os, p_ice, p_oce_sfc, trac_no, stretch_c, lacc=lzacc)
      ENDIF

      !  apply restoring to surface salinity directly
      CALL apply_surface_relaxation(p_patch_3D, p_os, p_oce_sfc, trac_no, lacc=lzacc)

    ENDIF

    !---------------------------------------------------------------------
    ! (2) Receive/calculate surface fluxes and wind stress
    !---------------------------------------------------------------------
    CALL update_atmos_fluxes(p_patch_3D, p_as, atmos_fluxes, p_oce_sfc, p_os, p_ice, this_datetime, lacc=lzacc)

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c,i_endidx_c
        ! copy atmospheric variables into new forcing variables for diagnostics
        p_oce_sfc%Wind_Speed_10m(jc,jb) = p_as%fu10(jc,jb)

    !---------------------------------------------------------------------
    ! (3) Sea ice thermodynamics & dynamics (at ocean time-step)
    !---------------------------------------------------------------------
        p_oce_sfc%cellThicknessUnderIce(jc,jb) = p_ice%zUnderIce(jc,jb) ! neccessary, because is not yet in restart

        p_ice%draftave_old(jc,jb) = p_ice%draftave(jc,jb)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)

    IF ( i_sea_ice > 0 ) THEN ! sea ice is on

        !  (3a) Fast sea ice thermodynamics (Analytical or OMIP cases only. Otherwise, done in the atmosphere)
        IF (iforc_oce == Analytical_Forcing .OR. iforc_oce == OMIP_FluxFromFile)  THEN
            CALL ice_fast_interface(p_patch, p_ice, atmos_fluxes, this_datetime, lacc=lzacc)
        ENDIF

        CALL ice_dynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff, lacc=lzacc)

        CALL ice_thermodynamics(p_patch_3D, p_ice, p_oce_sfc, atmos_fluxes, p_os, p_as, p_op_coeff, lacc=lzacc)

    ELSE !  sea ice is off

      ! for the setup without sea ice the SST is set to freezing temperature Tf
      ! should not be done here! Move to apply_surface_fluxes
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          IF(p_oce_sfc%SST(jc,jb) .LT. Tf) p_oce_sfc%SST(jc,jb) = Tf
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

    ENDIF

    !---------------------------------------------------------------------
    ! (4) Ocean surface stress boundary condition (atm-ocean + ice-ocean)
    !---------------------------------------------------------------------
    ! atm-ocean stress is either from file, bulk formula, or coupling
    ! if ice dynamics is off, ocean-ice stress is set to zero (no friction with ice)

    CALL update_ocean_surface_stress(p_patch_3D, p_ice, p_os, atmos_fluxes, p_oce_sfc, lacc=lzacc)

    !---------------------------------------------------------------------
    ! (5) Apply thermal and haline fluxes to the ocean surface layer
    !---------------------------------------------------------------------
    ! calculate the sw flux used for subsurface heating

    ! include hamoccs chlorophylls effect sw absorption
    ! FIXME zstar: Haven't checked for zstar requirements in hamocc
    IF ( lhamocc .AND. lfb_bgc_oce ) CALL dynamic_swr_absorption(p_patch_3d, p_os, lacc=lzacc)


    IF ( lswr_jerlov ) THEN

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          p_os%p_diag%heatabs(jc,jb)=(p_os%p_diag%swsum(jc,jb)  &
                  *p_oce_sfc%HeatFlux_ShortWave(jc,jb)*(1.0_wp-p_ice%concsum(jc,jb)))
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

    ELSE

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          p_os%p_diag%heatabs(jc,jb)=0.0_wp
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

    END IF

    IF (vert_cor_type .EQ. 0 ) THEN
       CALL apply_surface_fluxes_slo(p_patch_3D, p_os, p_ice, p_oce_sfc, lacc=lzacc)
    ELSE
       CALL apply_surface_fluxes_zstar_v13(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c, lacc=lzacc)
    ENDIF

    ! apply subsurface heating
    IF ( lswr_jerlov ) THEN

      IF (vert_cor_type .EQ. 0) THEN
        CALL subsurface_swr_absorption(p_patch_3d, p_os, lacc=lzacc)
      ELSE
        CALL subsurface_swr_absorption_zstar(p_patch_3d, p_os, stretch_c, lacc=lzacc)
      ENDIF

    ENDIF

    !---------------------------------------------------------------------
    ! (6) Apply volume flux correction
    !---------------------------------------------------------------------
    !  - sea level is balanced to zero over ocean surface
    !  - correction applied daily
    !  - calculate time
    dsec  = REAL(getNoOfSecondsElapsedInDayDateTime(this_datetime), wp)
    IF (limit_elevation) THEN

      IF (vert_cor_type .EQ. 0) THEN
        CALL balance_elevation(p_patch_3D, p_os%p_prog(nold(1))%h, p_oce_sfc, p_ice, lacc=lzacc)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
        CALL dbg_print('UpdSfc: h-old+BalElev',p_os%p_prog(nold(1))%h  ,str_module, 3, in_subset=p_patch%cells%owned)
      ELSE
        CALL balance_elevation_zstar(p_patch_3D, eta_c, p_oce_sfc, stretch_c, lacc=lzacc)
      !---------DEBUG DIAGNOSTICS-------------------------------------------
        CALL dbg_print('UpdSfc: h-old+BalElev', eta_c, routine, 2, in_subset=p_patch%cells%owned)
      !---------------------------------------------------------------------
      ENDIF

    END IF

  END SUBROUTINE update_ocean_surface_refactor

 !-------------------------------------------------------------------------
  !>
  !! Apply Thermodynamic Equations for Thermal and Haline Boundary Conditions
  !!
  !! Adapted for zstar
  !
  SUBROUTINE apply_surface_fluxes_zstar_v13(p_patch_3D, p_os, p_ice, p_oce_sfc, eta_c, stretch_c, lacc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)        :: p_patch_3D
    TYPE(t_hydro_ocean_state)                   :: p_os
    TYPE(t_sea_ice)                             :: p_ice
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    !
    REAL(wp), INTENT(INOUT) :: eta_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks) !! sfc ht
    REAL(wp), INTENT(IN   ) :: stretch_c(nproma, p_patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    ! local variables
    INTEGER               :: jc, jb, jt
    INTEGER               :: i_startidx_c, i_endidx_c
#ifdef __LVECTOR__
    REAL(wp)              :: dz_old
    REAL(wp)              :: dz_new(nproma,n_zlev), z_change(nproma,n_zlev)
    REAL(wp)              :: temp_stretch(nproma)
    INTEGER               :: lev, max_lev
#else
    REAL(wp)              :: dz_old(n_zlev), dz_new(n_zlev), z_change(n_zlev)
    REAL(wp)              :: temp_stretch
#endif
    REAL(wp)              :: tdc_i
    LOGICAL               :: lzacc

    REAL(wp) :: min_h

    INTEGER  :: bt_lev, jk, adj_lev
    REAL(wp) :: d_c
    REAL(wp) :: dz_ratio, adj_lev_ht

    REAL(wp) :: old_sss, new_s1, new_s2, new_sss, thresh_sss
    INTEGER  :: flag

    REAL(wp)  :: heatflux_surface_layer ! heatflux into the surface layer

    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    REAL(wp), POINTER :: prism_thick_flat_sfc_c(:,:,:)
    REAL(wp), POINTER :: tracer(:,:,:,:)
    INTEGER, POINTER :: dolic_c(:,:)

    CHARACTER(LEN=max_char_length), PARAMETER :: str_module = 'apply_surface_fluxes_zstar_v13'

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    tracer          => p_os%p_prog(nold(1))%tracer
    !-----------------------------------------------------------------------
    prism_thick_flat_sfc_c => p_patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(:,:,:)
    dolic_c => p_patch_3d%p_patch_1d(1)%dolic_c(:,:)
    tracer => p_os%p_prog(nold(1))%tracer(:,:,:,:)

    CALL dbg_print('UpdSfcSTART: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcSTART: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c
        !!  Provide total ocean forcing:
        !    - total heat fluxes are aggregated for ice/ocean in ice thermodynamics
        !    - total internal salt flux p_oce_sfc%FrshFlux_TotalIce is calculated in sea ice model
        !    - total freshwater volume forcing
        p_oce_sfc%FrshFlux_VolumeTotal(jc,jb) = p_oce_sfc%FrshFlux_Runoff    (jc,jb) &
          &                                   + p_oce_sfc%FrshFlux_VolumeIce (jc,jb) &
          &                                   + p_oce_sfc%FrshFlux_TotalOcean(jc,jb)
        ! provide total salinity forcing flux for diagnostics only
        p_oce_sfc%FrshFlux_TotalSalt(jc,jb)   = p_oce_sfc%FrshFlux_Runoff    (jc,jb) &
          &                                   + p_oce_sfc%FrshFlux_TotalIce  (jc,jb) &
          &                                   + p_oce_sfc%FrshFlux_TotalOcean(jc,jb)
      END DO
    END DO
    !$ACC WAIT(1)

    IF (no_tracer > 0) THEN
      IF ( heatflux_forcing_on_sst ) THEN
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(heatflux_surface_layer) ASYNC(1) IF(lzacc)
          DO jc = i_startidx_c, i_endidx_c
            IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN

              p_os%p_diag%heatflux_rainevaprunoff(jc,jb) =                                   &
                  ( p_oce_sfc%FrshFlux_TotalOcean(jc,jb)+ p_oce_sfc%FrshFlux_Runoff(jc,jb) )    &
                    * ( p_oce_sfc%sst(jc,jb) - tf ) * clw * rho_ref

              ! substract the fraction of heatflux used for subsurface heating

              heatflux_surface_layer=p_oce_sfc%HeatFlux_Total(jc,jb)-p_os%p_diag%heatabs(jc,jb)


              p_oce_sfc%sst(jc,jb) = p_oce_sfc%sst(jc,jb) + &
                &                    heatflux_surface_layer*dtime/(clw*rho_ref*p_oce_sfc%cellThicknessUnderIce(jc,jb))

            ENDIF
          ENDDO
          !$ACC END PARALLEL LOOP
        ENDDO
        !$ACC WAIT(1)
      ENDIF

      CALL dbg_print('UpdSfc: eta-old', eta_c,    str_module, 1, in_subset=p_patch%cells%owned)
      
      ! apply volume flux to surface elevation
#ifdef __LVECTOR__
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

        max_lev = MAXVAL(dolic_c(i_startidx_c:i_endidx_c, jb))

        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(dz_new, z_change, temp_stretch) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c, i_endidx_c
          IF (dolic_c(jc,jb) > 0) THEN

            bt_lev = dolic_c(jc, jb)
            d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)

            eta_c(jc,jb) = eta_c(jc,jb)               &
              &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
              &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime

            !! Only change the stretching parameter if it is above a certain threshold
            !! This avoids divide by 0
            min_h = prism_thick_flat_sfc_c(jc,1,jb)

            !! Update only if height is atleast dz
            IF (d_c > min_h) THEN
              temp_stretch(jc) = ( eta_c(jc, jb) + d_c)/( d_c )
            ELSE
              temp_stretch(jc) = stretch_c(jc, jb)
            END IF

            !! set dilution coefficient for HAMOCC
            p_oce_sfc%top_dilution_coeff(jc,jb) = stretch_c(jc,jb)/temp_stretch(jc)

            !! update zunderice
            p_ice%zUnderIce(jc,jb) = prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch(jc)

            dz_old = prism_thick_flat_sfc_c(jc,bt_lev,jb) * stretch_c(jc,jb)
            dz_new(jc,bt_lev) = prism_thick_flat_sfc_c(jc,bt_lev,jb) * temp_stretch(jc)
            z_change(jc,bt_lev) = dz_new(jc,bt_lev) - dz_old
          END IF
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        !$ACC LOOP SEQ
        DO jk = max_lev-1,1,-1
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx_c, i_endidx_c
            IF (jk <= dolic_c(jc,jb) - 1) THEN
              dz_old = prism_thick_flat_sfc_c(jc,jk,jb) * stretch_c(jc,jb)
              dz_new(jc,jk) = prism_thick_flat_sfc_c(jc,jk,jb) * temp_stretch(jc)
              z_change(jc,jk) = z_change(jc,jk+1) + dz_new(jc,jk) - dz_old
            END IF
          END DO
        END DO

        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx_c, i_endidx_c
          IF (p_oce_sfc%top_dilution_coeff(jc,jb) > 1.0_wp .AND. 1 <= dolic_c(jc,jb)) THEN
            dz_old = prism_thick_flat_sfc_c(jc,1,jb) * stretch_c(jc,jb)

            !NEC$ unroll_complete
            DO jt = 1, 2
              IF (jt /= 1 .OR. .NOT. lfwflux_enters_with_sst) THEN
                tracer(jc, 1, jb, jt) = tracer(jc, 1, jb, jt) * dz_old / ( dz_old + z_change(jc,1))
              ENDIF
            END DO
          END IF
        END DO

        !$ACC LOOP SEQ
        DO jk = 1, max_lev-1
          lev = max_lev + 1 - jk ! reverse iteration (max_lev -> 2) for growing layers

          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx_c, i_endidx_c
            IF (p_oce_sfc%top_dilution_coeff(jc,jb) > 1.0_wp .AND. jk <= dolic_c(jc,jb)-1) THEN

              !NEC$ unroll_complete
              DO jt = 1, 2
                tracer(jc, jk, jb, jt) = (tracer(jc, jk, jb, jt) &
                  &   * (dz_new(jc,jk) + z_change(jc,jk+1)) &
                  &   - z_change(jc,jk+1) * tracer(jc, jk+1, jb, jt))  &
                  &   / dz_new(jc,jk)
              END DO

            ELSE IF (p_oce_sfc%top_dilution_coeff(jc,jb) < 1.0_wp .AND. lev <= dolic_c(jc,jb)) THEN

              !NEC$ unroll_complete
              DO jt = 1, 2
                tracer(jc, lev, jb, jt) = (tracer(jc, lev, jb, jt) &
                  &   * (dz_new(jc,lev) - z_change(jc,lev)) &
                  &   + z_change(jc,lev) * tracer(jc, lev-1, jb, jt))  &
                  &   / dz_new(jc,lev)
              END DO

            END IF
          END DO ! jc
        END DO ! jk, lev

        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx_c, i_endidx_c
          IF (p_oce_sfc%top_dilution_coeff(jc,jb) < 1.0_wp .AND. 1 <= dolic_c(jc,jb)) THEN
            !NEC$ unroll_complete
            DO jt = 1, 2
              IF (jt /= 1 .OR. .NOT. lfwflux_enters_with_sst) THEN
                tracer(jc, 1, jb, jt) = tracer(jc, 1, jb, jt) * &
                    &  (dz_new(jc,1) - z_change(jc,1)) / dz_new(jc,1)
              ENDIF
            END DO
          END IF
        END DO

        !! Now put salt from ice into the surface layer (sss aliases the first level of the tracer field!)
        !$ACC LOOP GANG VECTOR
        DO jc = i_startidx_c, i_endidx_c
          IF (1 <= dolic_c(jc,jb)) THEN
            p_oce_sfc%sss(jc,jb) = p_oce_sfc%sss(jc,jb) + p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime / p_ice%zUnderIce(jc,jb)
          END IF
        END DO
        !$ACC END PARALLEL
        !$ACC WAIT(1)
      END DO ! jb
#else
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL DEFAULT(PRESENT) PRIVATE(dz_old, dz_new, z_change) ASYNC(1) IF(lzacc)
        !$ACC LOOP SEQ
        DO jc = i_startidx_c, i_endidx_c
          IF (dolic_c(jc,jb) > 0) THEN

            bt_lev = dolic_c(jc, jb)
            d_c    = p_patch_3d%p_patch_1d(1)%depth_CellInterface(jc, bt_lev + 1, jb)

            eta_c(jc,jb) = eta_c(jc,jb)               &
              &           + p_oce_sfc%FrshFlux_VolumeTotal(jc,jb)*dtime &
              &           + p_oce_sfc%FrshFlux_TotalIce(jc, jb)*dtime

            !! Only change the stretching parameter if it is above a certain threshold
            !! This avoids divide by 0
            temp_stretch = stretch_c(jc, jb)
            min_h        = prism_thick_flat_sfc_c(jc,1,jb)

            !! Update only if height is atleast dz
            if ( d_c  .GT.  min_h ) &
              & temp_stretch = ( eta_c(jc, jb) + d_c)/( d_c )

            !! set dilution coefficient for HAMOCC
            p_oce_sfc%top_dilution_coeff(jc,jb) = stretch_c(jc,jb)/temp_stretch

            !! update zunderice
            p_ice%zUnderIce(jc,jb) = prism_thick_flat_sfc_c(jc,1,jb) * temp_stretch

            !! Now interpolate T and S from old levels to new levels

            !! Get old and new level thicknesses and accumulated change of levels
            !! from bottom to jk
            dz_old(:) = prism_thick_flat_sfc_c(jc, :, jb) * stretch_c(jc, jb)
            !! dz_new = dz_old / tdc
            dz_new(:) = prism_thick_flat_sfc_c(jc, :, jb) * temp_stretch

            z_change(bt_lev) = dz_new(bt_lev) - dz_old(bt_lev)
            !$ACC LOOP SEQ
            DO jk = bt_lev-1,1,-1
              z_change(jk) = z_change(jk+1) + dz_new(jk) - dz_old(jk)
            ENDDO

            !! Loop through T and S
            !$ACC LOOP SEQ
            DO jt = 1,2

              IF (p_oce_sfc%top_dilution_coeff(jc,jb) > 1.0_wp) THEN
                ! add surface fwflux with temperature of sst, dilute sss
                IF ( jt .EQ. 1) THEN
                  IF (.NOT. lfwflux_enters_with_sst ) THEN
                    tracer(jc, 1, jb, jt) = tracer(jc, 1, jb, jt) * &
                        &  dz_old(1) / (dz_old(1) + z_change(1))
                  ENDIF
                ELSE
                  tracer(jc, 1, jb, jt) = tracer(jc, 1, jb, jt) * &
                      &  dz_old(1) / (dz_old(1) + z_change(1))
                ENDIF

                ! If tdc > 1 then bottom layer conc is unaffected (levels have
                ! become thinner)
                !$ACC LOOP SEQ
                DO jk = 1,bt_lev-1
                  tracer(jc, jk, jb, jt) = ( tracer(jc, jk, jb, jt) &
                      &   * (dz_new(jk) + z_change(jk+1)) &
                      &   - z_change(jk+1) * tracer(jc, jk+1, jb, jt))  &
                      &   / dz_new(jk)
                ENDDO
              ELSEIF (p_oce_sfc%top_dilution_coeff(jc,jb) < 1.0_wp) THEN
                !$ACC LOOP SEQ
                DO jk = bt_lev,2,-1
                  tracer(jc, jk, jb, jt) = ( tracer(jc, jk, jb, jt) &
                      &   * (dz_new(jk) - z_change(jk)) &
                      &   + z_change(jk) * tracer(jc, jk-1, jb, jt))  &
                      &   / dz_new(jk)
                ENDDO

                IF ( jt .EQ. 1) THEN
                  IF (.NOT. lfwflux_enters_with_sst ) THEN
                    tracer(jc, 1, jb, jt) = tracer(jc, 1, jb, jt) * &
                        &  (dz_new(1) - z_change(1)) / dz_new(1)
                  ENDIF
                ELSE
                    tracer(jc, 1, jb, jt) = tracer(jc, 1, jb, jt) * &
                        &  (dz_new(1) - z_change(1)) / dz_new(1)
                ENDIF

              ENDIF

          ENDDO ! T and S

            !! Now put salt from ice into the surface layer
            p_oce_sfc%sss(jc,jb) = p_oce_sfc%sss(jc,jb) + p_oce_sfc%FrshFlux_IceSalt(jc,jb) * dtime / p_ice%zUnderIce(jc,jb)

          ENDIF  !  dolic>0
        END DO
        !$ACC END PARALLEL
      END DO
      !$ACC WAIT(1)
#endif

    ENDIF

    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c
        !! set correct cell thickness under ice
        p_oce_sfc%cellThicknessUnderIce   (jc,jb) = p_ice%zUnderIce(jc,jb)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)

    CALL dbg_print('UpdSfcEND: oce_sfc%SST ',p_oce_sfc%SST, str_module, 2, in_subset=p_patch%cells%owned)
    CALL dbg_print('UpdSfcEND: oce_sfc%SSS ',p_oce_sfc%SSS, str_module, 2, in_subset=p_patch%cells%owned)

  END SUBROUTINE apply_surface_fluxes_zstar_v13



  !-------------------------------------------------------------------------
  !
  !>
  !!  Update surface fluxes for ocean forcing. Analytical, OMIP, coupled.
  !!
  !!  OMIP: atmos_fluxes over ice and open water are calculated with the help of
  !!        bulk formulas in calc_omip_budgets_ice and calc_omip_budgets_oce
  !!  Coupled: passes E/P fluxes and heat fluxes over open ocean. Heat fluxes over
  !!           the ice surface are not updated, because they are only used in ice_fast,
  !!           and ice_fast is called within the atmosphere.
  !!
  !
!<Optimize_Used>
  SUBROUTINE update_atmos_fluxes(p_patch_3D, p_as, atmos_fluxes, p_oce_sfc, p_os, p_ice, this_datetime, lacc)

    TYPE(t_patch_3D ),TARGET,  INTENT(IN)       :: p_patch_3D
    TYPE(t_atmos_for_ocean)                     :: p_as
    TYPE(t_atmos_fluxes)                        :: atmos_fluxes
    TYPE(t_ocean_surface)                       :: p_oce_sfc
    TYPE (t_hydro_ocean_state),INTENT(IN)       :: p_os
    TYPE (t_sea_ice),          INTENT(IN)       :: p_ice
    TYPE(datetime), POINTER                     :: this_datetime
    LOGICAL, INTENT(IN), OPTIONAL               :: lacc

    ! local variables
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface_refactor:update_atmos_fluxes'
    TYPE(t_patch), POINTER:: p_patch
    TYPE(t_subset_range), POINTER :: all_cells
    LOGICAL                       :: lzacc
    INTEGER                       :: jc, jb, i_startidx_c, i_endidx_c
    REAL(wp)                      :: ftdew_in(SIZE(p_as%ftdew,1), SIZE(p_as%ftdew,2))

    CALL set_acc_host_or_device(lzacc, lacc)

    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    all_cells       => p_patch%cells%all
    !-----------------------------------------------------------------------

    SELECT CASE (iforc_oce)

    CASE (Analytical_Forcing)        !  11      !  Driving the ocean with analytical fluxes

      CALL update_atmos_fluxes_analytical(p_patch_3D, p_os, p_ice, atmos_fluxes, p_oce_sfc, lacc=lzacc)

    CASE (OMIP_FluxFromFile)         !  12      !  Driving the ocean with OMIP fluxes

      !   a) read OMIP data into p_as
      CALL update_flux_fromFile(p_patch_3D, p_as, this_datetime, lacc=lzacc)

      !   b) calculate heat fluxes from p_as
      CALL calc_omip_budgets_oce(p_patch_3d, p_as, p_os, p_ice, atmos_fluxes, lacc=lzacc)

      IF (i_sea_ice >= 1) THEN ! sea ice is on

          !$ACC DATA CREATE(ftdew_in) IF(lzacc)

          !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          ftdew_in(:,:) = p_as%ftdew(:,:)-tmelt
          !$ACC END KERNELS
          !$ACC WAIT(1)

          CALL calc_omip_budgets_ice(                            &
            &                        p_patch_3d,                                     &  !  input parameter
            &                        p_as%tafo(:,:), ftdew_in(:,:),                  &  !  input parameter
            &                        p_as%fu10(:,:), p_as%fclou(:,:), p_as%pao(:,:), &  !  input parameter
            &                        p_as%fswr(:,:), p_ice%kice, p_ice%Tsurf(:,:,:), &  !  input parameter
            &                        p_ice%hi(:,:,:),                                &  !  input parameter
            &                        atmos_fluxes%albvisdir(:,:,:),                  &  !  input parameter
            &                        atmos_fluxes%albvisdif(:,:,:),                  &  !  input parameter
            &                        atmos_fluxes%albnirdir(:,:,:),                  &  !  input parameter
            &                        atmos_fluxes%albnirdif(:,:,:),                  &  !  input parameter
            &                        atmos_fluxes%LWnet    (:,:,:),                  &  !  output parameter
            &                        atmos_fluxes%SWnet    (:,:,:),                  &  !  output parameter
            &                        atmos_fluxes%sens     (:,:,:),                  &  !  output parameter
            &                        atmos_fluxes%lat      (:,:,:),                  &  !  output parameter
            &                        atmos_fluxes%dLWdT    (:,:,:),                  &  !  output parameter
            &                        atmos_fluxes%dsensdT  (:,:,:),                  &  !  output parameter
            &                        atmos_fluxes%dlatdT   (:,:,:),                  &  !  output parameter
            &                        lacc=lzacc)

          !$ACC END DATA

      ELSE   !  no sea ice

#ifdef _OPENACC
        IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

        ! apply net surface heat flux in W/m2 for OMIP case, since these fluxes are calculated in calc_omip_budgets_oce
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO jc = i_startidx_c,i_endidx_c
            IF (p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary) THEN
              p_oce_sfc%HeatFlux_ShortWave(jc,jb) = atmos_fluxes%SWnetw(jc,jb) ! net SW radiation flux over water
              p_oce_sfc%HeatFlux_LongWave (jc,jb) = atmos_fluxes%LWnetw(jc,jb) ! net LW radiation flux over water
              p_oce_sfc%HeatFlux_Sensible (jc,jb) = atmos_fluxes%sensw (jc,jb) ! Sensible heat flux over water
              p_oce_sfc%HeatFlux_Latent   (jc,jb) = atmos_fluxes%latw  (jc,jb) ! Latent heat flux over water
              ! sum of ocean heat fluxes for ocean boundary condition without ice, generally aggregated in ice thermodynamics
              p_oce_sfc%HeatFlux_Total(jc,jb) = atmos_fluxes%SWnetw(jc,jb) + atmos_fluxes%LWnetw(jc,jb) &
                &                              + atmos_fluxes%sensw(jc,jb)  + atmos_fluxes%latw(jc,jb)
            ELSE
              p_oce_sfc%HeatFlux_ShortWave(jc,jb) = 0.0_wp
              p_oce_sfc%HeatFlux_LongWave (jc,jb) = 0.0_wp
              p_oce_sfc%HeatFlux_Sensible (jc,jb) = 0.0_wp
              p_oce_sfc%HeatFlux_Latent   (jc,jb) = 0.0_wp
              p_oce_sfc%HeatFlux_Total    (jc,jb) = 0.0_wp
            ENDIF
          END DO
          !$ACC END PARALLEL LOOP
        END DO
        !$ACC WAIT(1)

      ENDIF

      ! c) wind stress is assigned in calc_omip_budgets_oce
      !    over ice: stress_x, stress_y; and over open water: stress_xw, stress_yw
      ! d) freshwater fluxes from p_as
      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c

          ! provide evaporation from latent heat flux for OMIP case
          ! under sea ice evaporation is neglected, atmos_fluxes%latw is flux in the absence of sea ice
          p_oce_sfc%FrshFlux_Evaporation(jc,jb) = atmos_fluxes%latw(jc,jb) / (alv*rho_ref)

          !  copy variables into atmos_fluxes
          p_oce_sfc%FrshFlux_Runoff(jc,jb)      = p_as%FrshFlux_Runoff(jc,jb)
          p_oce_sfc%FrshFlux_Precipitation(jc,jb) = p_as%FrshFlux_Precipitation(jc,jb)

          ! Precipitation on ice is snow when tsurf is below the freezing point
          !  - no snowfall from OMIP data
          !  - rprecw, rpreci are water equivalent over whole grid-area
          IF ( ALL( p_ice%Tsurf(jc,:,jb) < 0._wp) )  THEN!  Tsurf is -1.8 over open water, incorrect specification
            atmos_fluxes%rpreci(jc,jb) = p_as%FrshFlux_Precipitation(jc,jb)
            atmos_fluxes%rprecw(jc,jb) = 0._wp
          ELSE
            ! not considered in ice_growth_zero
            atmos_fluxes%rpreci(jc,jb) = 0._wp
            atmos_fluxes%rprecw(jc,jb) = p_as%FrshFlux_Precipitation(jc,jb)
          END IF

          ! evaporation and runoff not used in sea ice but in VolumeTotal, evaporation used for TotalOcean only
          p_oce_sfc%FrshFlux_TotalOcean(jc,jb) = p_patch_3d%wet_c(jc,1,jb)*( 1.0_wp-p_ice%concSum(jc,jb) ) * &
            &  (p_as%FrshFlux_Precipitation(jc,jb) + p_oce_sfc%FrshFlux_Evaporation(jc,jb))
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

    CASE (Coupled_FluxFromAtmo)

#ifdef _OPENACC
      IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

      !  Driving the ocean in a coupled mode:
      !  nothing to be done, atmospheric fluxes are provided at the end of time stepping
      !  atmospheric fluxes drive the ocean; fluxes are calculated by atmospheric model
      !  use atmospheric fluxes directly, i.e. no bulk formula as for OMIP is applied

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c

          ! HAMOCC uses p_as to get SW radiation and wind, so we need to copy
          ! the SW radiation onto it in the coupled case
          IF (lhamocc) p_as%fswr(jc,jb) = p_oce_sfc%HeatFlux_ShortWave(jc,jb)

          ! heatflux_total(jc,jb) is provided by coupling interface

          ! these 4 fluxes over open ocean are used in sea ice thermodynamics
          atmos_fluxes%SWnetw (jc,jb)   = p_oce_sfc%HeatFlux_ShortWave(jc,jb)
          atmos_fluxes%LWnetw (jc,jb)   = p_oce_sfc%HeatFlux_LongWave (jc,jb)
          atmos_fluxes%sensw  (jc,jb)   = p_oce_sfc%HeatFlux_Sensible (jc,jb)
          atmos_fluxes%latw   (jc,jb)   = p_oce_sfc%HeatFlux_Latent   (jc,jb)

          IF ( p_ice%concSum(jc,jb) > 0._wp) THEN !  corresponding to (1-concSum)*Precip in TotalOcean
          !  WHERE ( ALL( p_ice%hi   (jc,jb,:) > 0._wp, 2 ) )  !  corresponding to hi>0 in ice_growth_zero
          !  WHERE ( ALL( p_ice%Tsurf(jc,jb,:) < 0._wp, 2 ) )  !  Tsurf is -1.8 over open water, incorrect specification
            ! SnowFall and liquid rain over ice-covered part of ocean are taken from the atmosphere model
            atmos_fluxes%rpreci(jc,jb) = p_oce_sfc%FrshFlux_SnowFall(jc,jb)
            atmos_fluxes%rprecw(jc,jb) = p_oce_sfc%FrshFlux_Precipitation(jc,jb) - p_oce_sfc%FrshFlux_SnowFall(jc,jb)
          ELSE
            ! not considered in ice_growth_zero
            atmos_fluxes%rpreci(jc,jb) = 0._wp
            atmos_fluxes%rprecw(jc,jb) = p_oce_sfc%FrshFlux_Precipitation(jc,jb)
          END IF

          ! copy flux for use in TotalOcean, since analytical/omip use p_as:
          !p_as%FrshFlux_Precipitation      = p_oce_sfc%FrshFlux_Precipitation

          ! total water flux over ice-free ocean water: P*(1-C)+E
          !  - whole evaporation over grid-box enters open ocean, this includes evaporation over sea ice covered part
          !  - snowfall is included as (melted) water equivalent
          !  - runoff is added to VolumeTotal below
          p_oce_sfc%FrshFlux_TotalOcean(jc,jb) = p_patch_3d%wet_c(jc,1,jb)* &
            &  (( 1.0_wp-p_ice%concSum(jc,jb) ) * p_oce_sfc%FrshFlux_Precipitation(jc,jb) + p_oce_sfc%FrshFlux_Evaporation(jc,jb))
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

    CASE DEFAULT

      CALL message(TRIM(routine), 'STOP: Ocean Forcing option not implemented' )
      CALL finish(TRIM(routine), 'CHOSEN FORCING OPTION DOES NOT EXIST - TERMINATE')

    END SELECT


    IF (zero_freshwater_flux) THEN
#ifdef _OPENACC
      IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = i_startidx_c,i_endidx_c
          ! since latw<>0. we must set evap and TotalOcean again to zero:
          p_oce_sfc%FrshFlux_Evaporation  (jc,jb) = 0.0_wp
          p_oce_sfc%FrshFlux_TotalOcean   (jc,jb) = 0.0_wp
          p_oce_sfc%FrshFlux_Precipitation(jc,jb) = 0.0_wp
          p_oce_sfc%FrshFlux_SnowFall     (jc,jb) = 0.0_wp
          p_oce_sfc%FrshFlux_Evaporation  (jc,jb) = 0.0_wp
          p_oce_sfc%FrshFlux_Runoff       (jc,jb) = 0.0_wp
          p_oce_sfc%FrshFlux_TotalOcean   (jc,jb) = 0.0_wp
          atmos_fluxes%rpreci             (jc,jb) = 0.0_wp
          atmos_fluxes%rprecw             (jc,jb) = 0.0_wp
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

    ENDIF

    !---------DEBUG DIAGNOSTICS-------------------------------------------
    !idt_src=5  ! output print level (1-5, fix)
    !  these fluxes are always zero - fluxes over ice-covered area are Qbot, Qtop only
    !CALL dbg_print('aftAtmFluxUpd:atmflx%LWnetIce', atmos_fluxes%LWnet   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%SensIce',  atmos_fluxes%sens    ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%LatentIce',atmos_fluxes%lat     ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%dsensdT'  ,atmos_fluxes%dsensdT ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%dlatdT'   ,atmos_fluxes%dlatdT  ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:atmflx%dLWdT'    ,atmos_fluxes%dLWdt   ,str_module,idt_src, in_subset=p_patch%cells%owned)
    !CALL dbg_print('aftAtmFluxUpd:stress_x'        ,atmos_fluxes%stress_x,str_module,idt_src, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------
    CALL dbg_print('aftAtmFluxUpd: Precipitation', p_oce_sfc%FrshFlux_Precipitation,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: Evaporation'  , p_oce_sfc%FrshFlux_Evaporation  ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: SnowFall'     , p_oce_sfc%FrshFlux_SnowFall     ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: Runoff'       , p_oce_sfc%FrshFlux_Runoff       ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: TotalOcean'   , p_oce_sfc%FrshFlux_TotalOcean   ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: rprecw'       , atmos_fluxes%rprecw                ,str_module, 3, in_subset=p_patch%cells%owned)
    CALL dbg_print('aftAtmFluxUpd: rpreci'       , atmos_fluxes%rpreci                ,str_module, 3, in_subset=p_patch%cells%owned)
    !---------------------------------------------------------------------

  END SUBROUTINE update_atmos_fluxes

!**********************************************************************
!------------------------------ Analytical ----------------------------
!**********************************************************************

  !-------------------------------------------------------------------------
  !
  !>
  !! Update surface flux forcing for hydrostatic ocean
  !!
  !!
  !
!<Optimize_Used>
  SUBROUTINE update_atmos_fluxes_analytical(p_patch_3D, p_os, p_ice, atmos_fluxes, p_oce_sfc, lacc)

    TYPE(t_patch_3D ),TARGET, INTENT(IN)    :: p_patch_3D
    TYPE(t_hydro_ocean_state), INTENT(IN)   :: p_os
    TYPE(t_sea_ice), INTENT(IN)             :: p_ice
    TYPE(t_atmos_fluxes)                    :: atmos_fluxes
    TYPE(t_ocean_surface)                   :: p_oce_sfc
    LOGICAL, INTENT(IN), OPTIONAL           :: lacc
    !
    ! local variables
    INTEGER :: jc, jb
    INTEGER :: i_startblk_c, i_endblk_c, start_cell_index, end_cell_index
    !INTEGER :: i_startblk_c, i_endblk_c, i_startidx_c, i_endidx_c
    !INTEGER :: rl_start_c, rl_end_c

    REAL(wp) :: z_lat, z_lon, z_lat_deg
    !REAL(wp) :: y_length               !basin extension in y direction in degrees
    REAL(wp) :: z_T_init(nproma,p_patch_3D%p_patch_2D(1)%alloc_cell_blocks)
    REAL(wp) :: z_perlat, z_perlon, z_permax, z_perwid, z_relax, z_dst
    INTEGER  :: z_dolic
    REAL(wp) :: z_temp_max, z_temp_min, z_temp_incr
    REAL(wp) :: center, length, zonal_waveno,amplitude
    REAL(wp) :: no_flux_length, south_bound, max_flux_y
    LOGICAL  :: lzacc

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_bulk_forcing:update_atmos_fluxes_analytical'
    !-------------------------------------------------------------------------
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_patch), POINTER :: p_patch
    !-----------------------------------------------------------------------
    p_patch         => p_patch_3D%p_patch_2D(1)
    !-------------------------------------------------------------------------
    all_cells => p_patch%cells%all

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish(routine, 'OpenACC version currently not tested/validated')
#endif

    ! atmosphere fluxes for analytical testcased similar to mo_ocean_initial_conditions:
    SELECT CASE (atmos_flux_analytical_type)

    CASE(0)  !  all fluxes are unchanged, zero as default
      CONTINUE

    CASE(101,102,103)  !  constant fluxes for test of sea-ice processes

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_cell_index, end_cell_index
          ! set LW/SW/sensible/latent heat fluxes over ice to constant
          atmos_fluxes%SWnet(jc,1,jb) = atmos_SWnet_const
          atmos_fluxes%LWnet(jc,1,jb) = atmos_LWnet_const
          atmos_fluxes%sens (jc,1,jb) = atmos_sens_const
          atmos_fluxes%lat  (jc,1,jb) = atmos_lat_const
          ! set LW/SW/sensible/latent heat fluxes over water to constant
          atmos_fluxes%SWnetw(jc,jb)  = atmos_SWnetw_const
          atmos_fluxes%LWnetw(jc,jb)  = atmos_LWnetw_const
          atmos_fluxes%sensw(jc,jb)   = atmos_sensw_const
          atmos_fluxes%latw(jc,jb)    = atmos_lat_const
          ! set water fluxes over water to constant
          p_oce_sfc%FrshFlux_Precipitation(jc,jb) = atmos_precip_const
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

      CASE(200)

        IF(no_tracer>=1.AND.type_surfRelax_Temp==0)THEN

          center = basin_center_lat * deg2rad
          length = basin_height_deg * deg2rad
          zonal_waveno = 2.0_wp
          amplitude    =10.0_wp

          DO jb = all_cells%start_block, all_cells%end_block
            CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
            !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            DO jc = start_cell_index, end_cell_index

              IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

                z_lat    = p_patch%cells%center(jc,jb)%lat
                z_lon    = p_patch%cells%center(jc,jb)%lon

                p_oce_sfc%data_surfRelax_Temp(jc,jb) =0.0_wp

                p_oce_sfc%TopBC_Temp_vdiff(jc,jb) &
                &= amplitude * COS(zonal_waveno*pi*(z_lat-center)/length)

                !IF(z_lat<= center-0.5_wp*length)p_oce_sfc%TopBC_Temp_vdiff(jc,jb)=0.0_wp
                !IF(z_lat>= center+0.0_wp*length)p_oce_sfc%TopBC_Temp_vdiff(jc,jb)=10.0_wp
                !IF(z_lat>= center+0.75_wp*length)p_oce_sfc%TopBC_Temp_vdiff(jc,jb)=&
                !& -10.0_wp!amplitude * (COS(zonal_waveno*pi*(z_lat-center)/length))
              ENDIF
            END DO
            !$ACC END PARALLEL LOOP
          END DO
          !$ACC WAIT(1)
        ENDIF

      CASE(201) ! Abernathey 2011

        IF(no_tracer>=1.AND.type_surfRelax_Temp==0)THEN

          no_flux_length = 3.0_wp * relax_width * deg2rad
          south_bound = (basin_center_lat - 0.5_wp * basin_height_deg) * deg2rad
          length      = basin_height_deg * deg2rad
          max_flux_y  = length - no_flux_length
          zonal_waveno =  2.5_wp
          amplitude    = forcing_HeatFlux_amplitude

          DO jb = all_cells%start_block, all_cells%end_block
            CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
            !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            DO jc = start_cell_index, end_cell_index
              p_oce_sfc%data_surfRelax_Temp(jc,jb)     = 0.0_wp
              p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = 0.0_wp

              z_lat = p_patch%cells%center(jc,jb)%lat - south_bound

              IF(p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary .AND. z_lat < max_flux_y) THEN

                p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = &
                  & - amplitude * COS(zonal_waveno * pi * z_lat/max_flux_y) &
                  & + forcing_HeatFlux_base

              ENDIF
            END DO
            !$ACC END PARALLEL LOOP
          END DO
          !$ACC WAIT(1)
        ENDIF

      CASE default

        CALL finish(routine, "unknown atmos_flux_analytical_type")

    END SELECT

    SELECT CASE (relax_analytical_type)

    CASE(27,30,32)

     IF(no_tracer>=1.AND.type_surfRelax_Temp/=0)THEN

       !y_length = basin_height_deg * deg2rad
       DO jb = all_cells%start_block, all_cells%end_block
         CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
         !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
         DO jc = start_cell_index, end_cell_index

           IF(p_patch_3D%lsm_c(jc,1,jb)<=sea_boundary)THEN

             z_T_init(jc,jb) = 20.0_wp- p_patch_3D%p_patch_1D(1)%zlev_m(1)*15.0_wp/4000.0_wp

             z_lat    = p_patch%cells%center(jc,jb)%lat
             z_lon    = p_patch%cells%center(jc,jb)%lon

             ! Add temperature perturbation at new values
             z_perlat = basin_center_lat + 0.1_wp*basin_height_deg
             z_perlon = basin_center_lon + 0.1_wp*basin_width_deg
             z_permax = 0.1_wp
             z_perwid = 10.0_wp

             z_relax  = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

             z_dolic  = p_patch_3D%p_patch_1D(1)%dolic_c(jc,jb)
             IF (z_dolic > MIN_DOLIC) THEN

               z_dst = sqrt((z_lat-z_perlat*deg2rad)**2+(z_lon-z_perlon*deg2rad)**2)

               IF(z_dst <= 5.0_wp*deg2rad)THEN
                 z_T_init = z_T_init &
                 &        + z_permax*exp(-(z_dst/(z_perwid*deg2rad))**2) &
                 &        * sin(pi*p_patch_3D%p_patch_1D(1)%zlev_m(1)/4000.0_wp)
               ENDIF
               ! up to here z_init is identically initialized than temperature

               !add local cold perturbation
               IF(z_dst <= 10.5_wp*deg2rad)THEN
                 z_T_init(jc,jb) = z_T_init(jc,jb) - exp(-(z_dst/(z_perwid*deg2rad))**2)
               ENDIF

               p_oce_sfc%data_surfRelax_Temp(jc,jb)     = z_T_init(jc,jb)

               p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = z_relax * &
                 &  ( p_oce_sfc%data_surfRelax_Temp(jc,jb)-p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )

             END IF
           ELSE
             atmos_fluxes%topBoundCond_windStress_cc(jc,jb)%x(:) = 0.0_wp
!            atmos_fluxes%topBoundCond_windStress_u(jc,jb)       = 0.0_wp
!            atmos_fluxes%topBoundCond_windStress_v(jc,jb)       = 0.0_wp
           ENDIF
         END DO
         !$ACC END PARALLEL LOOP
       END DO
       !$ACC WAIT(1)

    ENDIF

    CASE (33)
      IF(type_surfRelax_Temp>=1)THEN
        z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO jc = start_cell_index, end_cell_index

            p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = z_relax*( p_oce_sfc%data_surfRelax_Temp(jc,jb) &
              &                                               -p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )
          END DO
          !$ACC END PARALLEL LOOP
        END DO
        !$ACC WAIT(1)

      END IF

    CASE(51)

      IF(type_surfRelax_Temp>=1)THEN

        z_relax = para_surfRelax_Temp/(30.0_wp*24.0_wp*3600.0_wp)

        z_temp_max  = 30.5_wp
        z_temp_min  = 0.5_wp
        z_temp_incr = (z_temp_max-z_temp_min)/(n_zlev-1.0_wp)

        !Add horizontal variation
        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          DO jc = start_cell_index, end_cell_index
            z_lat = p_patch%cells%center(jc,jb)%lat
            z_lat_deg = z_lat*rad2deg

            IF ( p_patch_3D%lsm_c(jc,1,jb) <= sea_boundary ) THEN

              z_temp_max     =0.01_wp*(z_lat_deg-basin_center_lat)*(z_lat_deg-basin_center_lat)
              z_T_init(jc,jb)=30.5_wp

              z_T_init(jc,jb)&
              &=z_T_init(jc,jb)*exp(-z_temp_max/basin_height_deg)
            ELSE
              z_T_init(jc,jb)=0.0_wp
            ENDIF

            p_oce_sfc%data_surfRelax_Temp(jc,jb)=z_T_init(jc,jb)

            p_oce_sfc%TopBC_Temp_vdiff(jc,jb) = z_relax*( p_oce_sfc%data_surfRelax_Temp(jc,jb) &
            &                                               -p_os%p_prog(nold(1))%tracer(jc,1,jb,1) )
          END DO
          !$ACC END PARALLEL LOOP
        END DO
        !$ACC WAIT(1)

      END IF

    END SELECT

    !-----------------------------
    !   varios ad-hoc fixes
    !-----------------------------
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_cell_index, end_cell_index
        !  needed for old heat flux BC
        p_oce_sfc%HeatFlux_Total(jc,jb)=p_oce_sfc%TopBC_Temp_vdiff(jc,jb)

        ! provide dLWdt for ice_fast as for OMIP
        atmos_fluxes%dLWdT (jc,jb,:)  = -4._wp*zemiss_def*stbo*(p_ice%tsurf(jc,jb,:)+tmelt)**3

        ! provide evaporation from latent heat flux
        p_oce_sfc%FrshFlux_Evaporation(jc,jb) = atmos_fluxes%latw(jc,jb) / (alv*rho_ref)
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)

  END SUBROUTINE update_atmos_fluxes_analytical

END MODULE mo_ocean_surface_refactor
