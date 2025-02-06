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

! Provide an implementation of the ocean forcing.
!
! Provide an implementation of the parameters used for surface forcing
! of the hydrostatic ocean model.

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_forcing
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_io_units,            ONLY: filename_max
  USE mo_grid_config,         ONLY: nroot
  USE mo_parallel_config,     ONLY: nproma
  USE mo_coupling_config,     ONLY: is_coupled_to_atmo
  USE mo_hamocc_nml,          ONLY: l_cpl_co2
  USE mo_ocean_nml,           ONLY: no_tracer,                                &
    & basin_height_deg, basin_width_deg, basin_center_lat, basin_center_lon,  &
    & forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno,       &
    & init_oce_relax, type_3dimRelax_Salt, type_3dimRelax_Temp,               &
    & type_surfRelax_Salt, type_surfRelax_Temp,                               &
    & forcing_windStress_u_amplitude, forcing_windStress_v_amplitude,         &
    & forcing_windstress_u_type, forcing_windstress_v_type,                   &
    & forcing_windspeed_type, forcing_windspeed_amplitude,                    &
    & relax_temperature_min, relax_temperature_max, initial_temperature_type, &
    & initial_temperature_bottom,initial_temperature_north,initial_temperature_south,&
    & relax_width, &
    & forcing_windstress_zonalWavePhase,                                      &
    & forcing_temperature_poleLat
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_math_constants,      ONLY: pi, deg2rad, pi_2
  USE mo_physical_constants,  ONLY: sfc_press_pascal
  USE mo_impl_constants,      ONLY: max_char_length, sea_boundary, success
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_ocean_surface_types, ONLY: t_ocean_surface, t_atmos_for_ocean
  USE mo_ocean_types,           ONLY: t_hydro_ocean_state
  USE mo_ocean_nudging,       ONLY: ocean_nudge
  USE mo_dynamics_config,     ONLY: nold

  USE mo_ocean_state,         ONLY: ocean_restart_list, ocean_default_list
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_var_list,            ONLY: add_var, add_ref
  USE mo_var_groups,          ONLY: groups
  USE mo_fortran_tools,       ONLY: assign_if_present
  USE mo_cf_convention
  USE mo_grib2
  USE mo_cdi,                ONLY: DATATYPE_FLT32, DATATYPE_FLT64, DATATYPE_PACK16, GRID_UNSTRUCTURED
  USE mo_zaxis_type,         ONLY: ZA_SURFACE
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_read_interface,     ONLY: openInputFile, closeFile, t_stream_id, &
    &                              on_cells, read_2D_1lev_1time
  USE mo_ocean_initial_conditions,     ONLY: tracer_ConstantSurface, &
    & varyTracerVerticallyExponentially
  USE mo_io_config,          ONLY: lnetcdf_flt64_output
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=12)           :: str_module    = 'oceForcing  '  ! Output of module for 1 line debug
  INTEGER :: idt_src       = 1               ! Level of detail for 1 line debug

  ! Public interface
  PUBLIC :: construct_ocean_surface, destruct_ocean_forcing
  PUBLIC :: construct_atmos_for_ocean, destruct_atmos_for_ocean
  PUBLIC :: init_ocean_forcing

CONTAINS


  !-------------------------------------------------------------------------
  !
  !> Constructor of ocean surface model, allocates all components and assigns zero.
  !!
  !! Combine with construct_ocean_forcing
  !
  SUBROUTINE construct_ocean_surface(p_patch_3D, p_oce_sfc)

    TYPE(t_patch_3D),             TARGET,INTENT(IN)    :: p_patch_3D
    TYPE (t_ocean_surface),              INTENT(INOUT) :: p_oce_sfc

    !Local variables
    INTEGER :: alloc_cell_blocks, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_surface:construct_ocean_surface'

    TYPE(t_patch),POINTER    :: p_patch
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    p_patch           => p_patch_3D%p_patch_2D(1)
    alloc_cell_blocks =  p_patch%alloc_cell_blocks
  ! nblks_e           =  p_patch%nblks_e

    !$ACC ENTER DATA COPYIN(p_oce_sfc)

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_default_list, 'top_dilution_coeff', p_oce_sfc%top_dilution_coeff, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('top_dilution_coeff', '', 'top_dilution_coeff', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), lrestart_cont=.false., in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%top_dilution_coeff)

    IF (is_coupled_to_atmo() .AND. l_cpl_co2) THEN
      CALL add_var(ocean_restart_list, 'co2_mixing_ratio', p_oce_sfc%CO2_Mixing_Ratio, &
        &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
        &        t_cf_var('Co2_Mixing_Ratio', 'ppmv', 'CO2 Mixing Ratio', datatype_flt),&
        &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
        &        ldims=(/nproma,alloc_cell_blocks/), lrestart_cont=.True., in_group=groups("oce_force_essentials"))
    END IF

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'sea_level_pressure', p_oce_sfc%sea_level_pressure, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('air_pressure_at_mean_sea_level', 'Pa', 'Sea Level Pressure', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'Wind_Speed_10m', p_oce_sfc%Wind_Speed_10m, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('wind_speed', 'm s-1', 'Wind Speed at 10m height', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! heat fluxes
    CALL add_var(ocean_restart_list, 'HeatFlux_Total', p_oce_sfc%HeatFlux_Total, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('surface_downward_heat_flux', 'W m-2', 'Total Heat Flux', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'HeatFlux_Shortwave', p_oce_sfc%HeatFlux_Shortwave, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('surface_downwelling_shortwave_flux', 'W m-2', 'Shortwave Heat Flux', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'HeatFlux_LongWave', p_oce_sfc%HeatFlux_LongWave , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('surface_downwelling_longwave_flux', 'W m-2', 'Longwave Heat Flux', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'HeatFlux_Sensible', p_oce_sfc%HeatFlux_Sensible , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('surface_downward_sensible_heat_flux', 'W m-2', 'Sensible Heat Flux', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'HeatFlux_Latent', p_oce_sfc%HeatFlux_Latent , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('surface_downward_latent_heat_flux', 'W m-2', 'Latent Heat Flux', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! freshwater fluxes
    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'FrshFlux_Precipitation', p_oce_sfc%FrshFlux_Precipitation , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('lwe_precipitation_rate', 'm s-1', 'FrshFlux_Precipitation', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'FrshFlux_Evaporation', p_oce_sfc%FrshFlux_Evaporation , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('lwe_water_evaporation_rate', 'm s-1', 'FrshFlux_Evaporation', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'FrshFlux_SnowFall', p_oce_sfc%FrshFlux_SnowFall , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('FrshFlux_SnowFall', 'm s-1', 'FrshFlux_SnowFall', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    ! Coupling fluxes must go into restart file:
    CALL add_var(ocean_restart_list, 'FrshFlux_Runoff', p_oce_sfc%FrshFlux_Runoff , &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('FrshFlux_Runoff', 'm s-1', 'FrshFlux_Runoff', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), &
      &        lrestart=is_coupled_to_atmo(), lrestart_cont=.TRUE., lopenacc=.NOT.is_coupled_to_atmo())

    IF (.NOT. is_coupled_to_atmo()) THEN
      __acc_attach(p_oce_sfc%sea_level_pressure)
      __acc_attach(p_oce_sfc%Wind_Speed_10m)
      __acc_attach(p_oce_sfc%HeatFlux_Total)
      __acc_attach(p_oce_sfc%HeatFlux_Shortwave)
      __acc_attach(p_oce_sfc%HeatFlux_LongWave)
      __acc_attach(p_oce_sfc%HeatFlux_Sensible)
      __acc_attach(p_oce_sfc%HeatFlux_Latent)
      __acc_attach(p_oce_sfc%FrshFlux_Precipitation)
      __acc_attach(p_oce_sfc%FrshFlux_Evaporation)
      __acc_attach(p_oce_sfc%FrshFlux_SnowFall)
      __acc_attach(p_oce_sfc%FrshFlux_Runoff)
    END IF

    CALL add_var(ocean_default_list, 'TopBC_WindStress_u', p_oce_sfc%TopBC_WindStress_u, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('TopBC_WindStress_u', 'Pa', 'Zonal Wind Stress', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%TopBC_WindStress_u)

    CALL add_var(ocean_default_list, 'TopBC_WindStress_v', p_oce_sfc%TopBC_WindStress_v, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('TopBC_WindStress_v', 'Pa', 'Meridional Wind Stress', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%TopBC_WindStress_v)

    CALL add_var(ocean_default_list, 'SST', p_oce_sfc%SST, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('SST', 'C', 'Sea Surface Temperature', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"))

    CALL add_var(ocean_default_list, 'SSS', p_oce_sfc%SSS, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('SSS', 'psu', 'Sea Surface Salinity', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"))

    CALL add_var(ocean_default_list,'SurfaceCellThicknessUnderIce', p_oce_sfc%cellThicknessUnderIce, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('SurfaceCellThicknessUnderIce', 'm', 'Cell Thickness at Surface under Ice', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%cellThicknessUnderIce)

    CALL add_var(ocean_default_list, 'surface_sal_content', p_oce_sfc%surface_salt_content, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('surface_salt_content', 'kg', 'surface_salt_content', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"))

    ! auxillary ice salt/freshwater fluxes,
    CALL add_var(ocean_default_list, 'FrshFlux_TotalIce', p_oce_sfc%FrshFlux_TotalIce, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('FrshFlux_TotalIce', 'm s-1', 'Freshwater Flux due to Sea Ice Change', datatype_flt),&
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_TotalIce)

    CALL add_var(ocean_default_list, 'FrshFlux_VolumeTotal', p_oce_sfc%FrshFlux_VolumeTotal, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('FrshFlux_VolumeTotal', 'm s-1', 'Freshwater Flux due to Volume Change', datatype_flt), &
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_VolumeTotal)

    CALL add_var(ocean_default_list, 'FrshFlux_IceSalt', p_oce_sfc%FrshFlux_IceSalt, &
      &        GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &        t_cf_var('FrshFlux_IceSalt', 'psu m s-1', 'Salt volume flux due to sea ice change', datatype_flt), &
      &        grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &        ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_IceSalt)

    ! auxillary freshwater fluxes, previously belonged to t_sfc_flx
    CALL add_var(ocean_default_list, 'FrshFlux_TotalSalt', p_oce_sfc%FrshFlux_TotalSalt , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_TotalSalt', 'm s-1', 'FrshFlux_TotalSalt', datatype_flt),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_TotalSalt)

    CALL add_var(ocean_default_list, 'FrshFlux_TotalOcean', p_oce_sfc%FrshFlux_TotalOcean , &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('surface_downward_water_flux', 'm s-1', 'FrshFlux_TotalOcean', datatype_flt),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_TotalOcean)

    CALL add_var(ocean_default_list, 'FrshFlux_VolumeIce', p_oce_sfc%FrshFlux_VolumeIce, &
    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
    &          t_cf_var('FrshFlux_VolumeIce', 'm s-1', 'FrshFlux_VolumeIce', datatype_flt),&
    &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
    &          ldims=(/nproma,alloc_cell_blocks/), in_group=groups("oce_force_essentials"), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_VolumeIce)

    ! relaxation fields {{{
!    IF(no_tracer>=1) THEN
!      ! there are four tracer related fields: tracer focing, tracer relaxation
!      ! and both accumulated
    CALL add_var(ocean_default_list, 'data_surfRelax_Temp', p_oce_sfc%data_surfRelax_Temp, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('data_surfRelax_Temp', 'C', 'Data to Relax Temperature to', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%data_surfRelax_Temp)

    CALL add_var(ocean_default_list, 'data_surfRelax_Salt', p_oce_sfc%data_surfRelax_Salt, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('data_surfRelax_Salt', 'psu', 'Data to Relax Salinity to', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%data_surfRelax_Salt)

    CALL add_var(ocean_default_list, 'TopBC_Temp_vdiff', p_oce_sfc%TopBC_Temp_vdiff, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('TopBC_Temp_vdiff', 'C m s-1', 'Forcing of temperature in vertical diffusion equation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.true., lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%TopBC_Temp_vdiff)

    CALL add_var(ocean_default_list, 'TopBC_Salt_vdiff', p_oce_sfc%TopBC_Salt_vdiff, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('TopBC_Salt_vdiff', 'psu m s-1', 'Forcing of salinity in vertical diffusion equation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), &
      &          lcontainer=.TRUE., lrestart=.FALSE., loutput=.true., lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%TopBC_Salt_vdiff)

    CALL add_var(ocean_default_list, 'TempFlux_Relax', p_oce_sfc%TempFlux_Relax , &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('TempFlux_Relax', 'K s-1', 'Temperature tracer flux due to relaxation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%TempFlux_Relax)

    CALL add_var(ocean_default_list, 'SaltFlux_Relax', p_oce_sfc%SaltFlux_Relax, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('SaltFlux_Relax', 'psu s-1', 'Salinity tracer flux due to relaxation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%SaltFlux_Relax)

    CALL add_var(ocean_default_list, 'HeatFlux_Relax', p_oce_sfc%HeatFlux_Relax, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('HeatFlux_Relax', 'W m-2', 'Surface heat flux due to relaxation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%HeatFlux_Relax)

    CALL add_var(ocean_default_list, 'FrshFlux_Relax', p_oce_sfc%FrshFlux_Relax, &
      &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, &
      &          t_cf_var('FrshFlux_Relax', 'm s-1', 'Surface freshwater flux due to relaxation', datatype_flt),&
      &          grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL),&
      &          ldims=(/nproma,alloc_cell_blocks/), lopenacc=.TRUE.)
    __acc_attach(p_oce_sfc%FrshFlux_Relax)
!    ENDIF
    ! }}}

    ! cartesians
    ALLOCATE(p_oce_sfc%TopBC_WindStress_cc(nproma,alloc_cell_blocks), STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for forcing wind stress cc  failed')
    END IF

    ! init of cartesian coordinates:
    p_oce_sfc%TopBC_WindStress_cc(:,:)%x(1) = 0.0_wp
    p_oce_sfc%TopBC_WindStress_cc(:,:)%x(2) = 0.0_wp
    p_oce_sfc%TopBC_WindStress_cc(:,:)%x(3) = 0.0_wp

    !$ACC ENTER DATA COPYIN(p_oce_sfc%TopBC_WindStress_cc)

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE construct_ocean_surface
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! Destructor surface flux forcing for hydrostatic ocean
  !!
  !
!<Optimize:inUse>
  SUBROUTINE destruct_ocean_forcing(p_oce_sfc)
    TYPE(t_ocean_surface), INTENT(INOUT) :: p_oce_sfc
    !
    ! Local variables

    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_ocean_forcing'
    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    ! forcing fields are handled by the ocean_default_list

    !$ACC EXIT DATA DELETE(p_oce_sfc%TopBC_WindStress_cc)

    DEALLOCATE(p_oce_sfc%TopBC_WindStress_cc, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for forcing wind stress cc failed')
    END IF

    CALL message(TRIM(routine), 'end' )

  END SUBROUTINE destruct_ocean_forcing
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !
  !>
  !! Constructor of atmospheric reprsentation  in ocean.
  !!
  !
  SUBROUTINE construct_atmos_for_ocean(p_patch, p_as)
    !
    TYPE(t_patch),                INTENT(IN):: p_patch
    TYPE(t_atmos_for_ocean ), INTENT(INOUT) :: p_as

    ! Local variables
    INTEGER :: alloc_cell_blocks, ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:construct_atmos_for_ocean'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    alloc_cell_blocks = p_patch%alloc_cell_blocks

    ALLOCATE(p_as%tafo(nproma,alloc_cell_blocks),                       &
             p_as%ftdew(nproma,alloc_cell_blocks),                      &
             p_as%fclou(nproma,alloc_cell_blocks),                      &
             p_as%fu10(nproma,alloc_cell_blocks),                       &
             p_as%co2(nproma,alloc_cell_blocks),                        &
             p_as%co2flx(nproma,alloc_cell_blocks),                     &
             p_as%fswr(nproma,alloc_cell_blocks),                       &
             p_as%pao(nproma,alloc_cell_blocks),                        &
             p_as%u(nproma,alloc_cell_blocks),                          &
             p_as%v(nproma,alloc_cell_blocks),                          &
!             p_as%precip(nproma,alloc_cell_blocks),                     &
!             p_as%evap(nproma,alloc_cell_blocks),                       &
!             p_as%runoff(nproma,alloc_cell_blocks),                     &
             p_as%topBoundCond_windStress_u(nproma,alloc_cell_blocks),  &
             p_as%topBoundCond_windStress_v(nproma,alloc_cell_blocks),  &
             p_as%FrshFlux_Precipitation(nproma,alloc_cell_blocks),     &
             p_as%FrshFlux_Runoff(nproma,alloc_cell_blocks),            &
             p_as%data_surfRelax_Temp(nproma,alloc_cell_blocks),        &
             p_as%data_surfRelax_Salt(nproma,alloc_cell_blocks), STAT=ist)

    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'allocation for p_as fields failed')
    END IF

    p_as%tafo  (:,:)                    = 0.0_wp
    p_as%ftdew (:,:)                    = 0.0_wp
    p_as%fclou (:,:)                    = 0.0_wp
    p_as%fu10  (:,:)                    = 0.0_wp
    p_as%co2  (:,:)                     = 0.0_wp
    p_as%co2flx  (:,:)                  = 0.0_wp
    p_as%fswr  (:,:)                    = 0.0_wp
    p_as%pao   (:,:)                    = sfc_press_pascal  ! initialize with reference pressure
    p_as%u     (:,:)                    = 0.0_wp
    p_as%v     (:,:)                    = 0.0_wp
!    p_as%precip(:,:)                    = 0.0_wp
!    p_as%evap  (:,:)                    = 0.0_wp
!    p_as%runoff(:,:)                    = 0.0_wp
    p_as%topBoundCond_windStress_u(:,:) = 0.0_wp
    p_as%topBoundCond_windStress_v(:,:) = 0.0_wp
    p_as%FrshFlux_Precipitation(:,:)    = 0.0_wp
    p_as%FrshFlux_Runoff(:,:)           = 0.0_wp
    p_as%data_surfRelax_Temp(:,:)       = 0.0_wp
    p_as%data_surfRelax_Salt(:,:)       = 0.0_wp

    !$ACC ENTER DATA &
    !$ACC   COPYIN(p_as, p_as%topBoundCond_windStress_u, p_as%topBoundCond_windStress_v) &
    !$ACC   COPYIN(p_as%tafo, p_as%ftdew, p_as%fu10, p_as%fclou, p_as%pao) &
    !$ACC   COPYIN(p_as%fswr, p_as%u, p_as%v, p_as%FrshFlux_Precipitation) &
    !$ACC   COPYIN(p_as%FrshFlux_Runoff, p_as%data_surfRelax_Temp, p_as%data_surfRelax_Salt) &
    !$ACC   COPYIN(p_as%co2, p_as%co2flx)

    CALL message(TRIM(routine), 'end')

  END SUBROUTINE construct_atmos_for_ocean

  !-------------------------------------------------------------------------
  !
  !>
  !!  Destructor of atmospheric reprsentation  in ocean.
  !!
  !
  SUBROUTINE destruct_atmos_for_ocean(p_as)
    !
    TYPE(t_atmos_for_ocean ), INTENT(INOUT) :: p_as

    ! Local variables
    INTEGER :: ist
    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_sea_ice:destruct_atmos_for_ocean'

    !-------------------------------------------------------------------------
    CALL message(TRIM(routine), 'start' )

    !$ACC EXIT DATA &
    !$ACC   DELETE(p_as%tafo, p_as%ftdew, p_as%fclou, p_as%fu10, p_as%co2) &
    !$ACC   DELETE(p_as%co2flx, p_as%fswr, p_as%pao, p_as%u, p_as%v)

    DEALLOCATE(p_as%tafo, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for tafo failed')
    END IF
    DEALLOCATE(p_as%ftdew, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for ftdew failed')
    END IF
    DEALLOCATE(p_as%fclou, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fclou failed')
    END IF

    DEALLOCATE(p_as%fu10, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fu10 failed')
    END IF

    DEALLOCATE(p_as%co2, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for co2 failed')
    END IF

    DEALLOCATE(p_as%co2flx, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for co2flx failed')
    END IF

    DEALLOCATE(p_as%fswr, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for fswr failed')
    END IF

    DEALLOCATE(p_as%pao, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for pao failed')
    END IF

    DEALLOCATE(p_as%u, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for u failed')
    END IF
    DEALLOCATE(p_as%v, STAT=ist)
    IF (ist/=SUCCESS) THEN
      CALL finish(TRIM(routine),'deallocation for v failed')
    END IF

!    DEALLOCATE(p_as%precip, STAT=ist)
!    IF (ist/=SUCCESS) THEN
!      CALL finish(TRIM(routine),'deallocation for precip failed')
!    END IF
!    DEALLOCATE(p_as%evap, STAT=ist)
!    IF (ist/=SUCCESS) THEN
!      CALL finish(TRIM(routine),'deallocation for evap failed')
!    END IF
!    DEALLOCATE(p_as%runoff, STAT=ist)
!    IF (ist/=SUCCESS) THEN
!      CALL finish(TRIM(routine),'deallocation for runoff failed')
!    END IF

    CALL message(TRIM(routine), 'end')

  END SUBROUTINE destruct_atmos_for_ocean

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE init_ocean_forcing(patch_2d, patch_3d, ocean_state, p_oce_sfc, fu10)
    !
    TYPE(t_patch),TARGET, INTENT(in)        :: patch_2d
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_ocean_surface), INTENT(INOUT)    :: p_oce_sfc
    REAL(wp), INTENT(INOUT)                 :: fu10(:,:)      !  windspeed: p_as%fu10

    TYPE(t_subset_range), POINTER :: all_cells, owned_cells


    CALL init_ocean_WindForcing(patch_3d, ocean_state, fu10)

    IF (init_oce_relax > 0) THEN
      CALL init_ho_relaxation(patch_3d, ocean_state, p_oce_sfc)
    END IF
  END SUBROUTINE init_ocean_forcing
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Initialization of temperature and salinity relaxation for the hydrostatic ocean model.
  !! Temperature and salinity relaxation data are read from external data
  !
  !
  !-------------------------------------------------------------------------
  !
!<Optimize:inUse>
  SUBROUTINE init_ho_relaxation(patch_3d, ocean_state, p_oce_sfc)

    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    TYPE(t_ocean_surface), INTENT(INOUT)    :: p_oce_sfc

    ! Local Variables

    CHARACTER(LEN=max_char_length), PARAMETER :: routine = 'mo_ocean_forcing:init_ho_relaxation'
    CHARACTER(filename_max) :: relax_init_file   !< file name for reading in

    LOGICAL :: l_exist
    INTEGER :: i_lev, no_cells, no_levels, jb, jc, jk
    INTEGER :: start_cell_index, end_cell_index
    INTEGER :: ncid, dimid
    TYPE(t_stream_id) :: stream_id

    REAL(wp):: z_c(nproma,1,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp):: z_surfRelax(nproma,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    REAL(wp) :: temperature_difference, poleLat, waveNo

    REAL(wp):: distan
    REAL(wp):: perturbation_lat, perturbation_lon,  max_perturbation, perturbation_width
    REAL(wp) :: basin_NorthBoundary, basin_SouthBoundary, lat_diff


    TYPE(t_patch),POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    !-------------------------------------------------------------------------
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL
    basin_NorthBoundary    = (basin_center_lat + 0.5_wp*basin_height_deg) * deg2rad
    basin_SouthBoundary    = (basin_center_lat - 0.5_wp*basin_height_deg) * deg2rad

    ! Read relaxation data from file
    IF (init_oce_relax == 1) THEN

       ! sphere_radius = grid_sphere_radius
      ! u0 =(2.0_wp*pi*sphere_radius)/(12.0_wp*24.0_wp*3600.0_wp)

      CALL message (TRIM(routine), 'start')

      i_lev        = patch_2d%level

      ! Relaxation variables are read from relax_init_file
!       WRITE (relax_init_file,'(a,i0,a,i2.2,a)') 'iconR',nroot,'B',i_lev, '-relax.nc'
      relax_init_file='ocean-relax.nc'

      IF (my_process_is_stdio()) THEN

        INQUIRE (FILE=relax_init_file, EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          WRITE(message_text,'(3a)') 'netcdf file named ', TRIM(relax_init_file),' not found!'
          CALL message(TRIM(routine),TRIM(message_text))
          CALL finish(TRIM(routine),'netcdf file for reading T/S relax. input not found - ABORT')
        ENDIF

        WRITE(message_text,'(3a)') 'netcdf file named ', TRIM(relax_init_file), &
          & ' opened for reading'
        CALL message(TRIM(routine),TRIM(message_text))

        !
        ! open file
        !
        CALL nf(nf90_open(TRIM(relax_init_file), nf90_nowrite, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf90_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf90_inquire_dimension(ncid, dimid, len = no_cells), routine)

        !
        ! check the number of cells
        !
        WRITE(message_text,'(a,i6)') 'No of cells =', no_cells
        CALL message(TRIM(routine),TRIM(message_text))
        IF (patch_2d%n_patch_cells_g /= no_cells) THEN
          CALL finish(TRIM(routine),&
            & 'Number of patch cells and cells in T/S relaxation input file do not match - ABORT')
        ENDIF
        !
        ! get number of levels
        !
        CALL nf(nf90_inq_dimid(ncid, 'level', dimid), routine)
        CALL nf(nf90_inquire_dimension(ncid, dimid, len = no_levels), routine)

        !
        ! check the number of cells
        !
        WRITE(message_text,'(a,i6)') 'No of vertical levels =', no_levels
        CALL message(TRIM(routine),TRIM(message_text))
        IF (no_levels /= 1) THEN
          CALL finish(TRIM(routine),'Number of vertical levels is not equal 1 - ABORT')
        ENDIF

        CALL nf(nf90_close(ncid), routine)

      ENDIF  !  stdio

      !-------------------------------------------------------
      !
      ! Read ocean relaxation data at cells
      !
      !-------------------------------------------------------

      CALL openinputfile(stream_id, relax_init_file, patch_2d)

      ! triangle center and edges

      ! read temperature
      !  - read one data set, annual mean only
      !  - "T": annual mean temperature
      CALL read_2D_1lev_1time(stream_id, on_cells, 'T', z_surfRelax)

      IF (no_tracer>=1) THEN
        p_oce_sfc%data_surfRelax_Temp(:,:) = z_surfRelax(:,:)
      ELSE
        CALL message( TRIM(routine),'WARNING: no tracer used, but init relaxation attempted')
      END IF

      ! read salinity
      !  - "S": annual mean salinity
      IF (no_tracer > 1) THEN
        CALL read_2D_1lev_1time(stream_id, on_cells, 'S', z_surfRelax)
        p_oce_sfc%data_surfRelax_Salt(:,:) = z_surfRelax(:,:)
      END IF

      ! close file
      CALL closeFile(stream_id)

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          IF ( patch_3d%lsm_c(jc,1,jb) > sea_boundary ) THEN
            p_oce_sfc%data_surfRelax_Temp(jc,jb) = 0.0_wp
            IF (no_tracer>1) p_oce_sfc%data_surfRelax_Salt(jc,jb) = 0.0_wp
          ENDIF
        END DO
      END DO

      CALL message( TRIM(routine),'Ocean T/S relaxation reading finished' )

    ENDIF  !  read relaxation data from file

    !-------------------------------------------------------
    !
    ! use initialized temperature/salinity, assigned to tracer, for 2-dim/3-dim relaxation
    !  - relaxation switch equals 3
    !
    !-------------------------------------------------------

    SELECT CASE(type_surfRelax_Temp)
    CASE(3)
      p_oce_sfc%data_surfRelax_Temp(:,:) = ocean_state%p_prog(nold(1))%tracer(:,1,:,1)

    CASE(4)
      ! smooth ape relaxation, as in temperature_smoothAPE in mo_cean_initial_conditions
      temperature_difference = relax_temperature_max - relax_temperature_min
      poleLat = ABS(forcing_temperature_poleLat * deg2rad)
      waveNo = pi_2 / poleLat

      DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
        DO jc = start_cell_index, end_cell_index
          DO jk=1, MIN(1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb))

            p_oce_sfc%data_surfRelax_Temp(jc,jb) = relax_temperature_min    + &
              & (COS(waveNo * MIN(ABS(patch_2d%cells%center(jc,jb)%lat), poleLat))**2) * temperature_difference

          END DO
        END DO
      END DO
      CASE(5)
      !Only valid for this specific testcase: temperature is resored to initial tempertaure without perturbation
      IF(initial_temperature_type==215.OR.initial_temperature_type==214)THEN

        p_oce_sfc%data_surfRelax_Temp(:,:) = 0.0_wp!initial_temperature_south

       temperature_difference = initial_temperature_north - initial_temperature_south
       lat_diff               = basin_NorthBoundary - basin_SouthBoundary  !  basin_height_deg*deg2rad

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            jk=1
            p_oce_sfc%data_surfRelax_Temp(jc,jb) = initial_temperature_north &
            &- temperature_difference*((basin_NorthBoundary-patch_2d%cells%center(jc,jb)%lat)/lat_diff)

            p_oce_sfc%data_surfRelax_Temp(jc,jb) = MERGE(ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1), &
            &initial_temperature_north, patch_2d%cells%center(jc,jb)%lat>basin_SouthBoundary)

            p_oce_sfc%data_surfRelax_Temp(jc,jb) = MERGE(&
            &initial_temperature_north-1.0_wp, &
                              &ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1), &
                              &patch_2d%cells%center(jc,jb)%lat>(basin_center_lat + 1.0_wp*basin_height_deg) * deg2rad)!&
                        !&.AND.lat(jc,jb)<(basin_center_lat + 1.25_wp*basin_height_deg) * deg2rad)
          END DO
        END DO

        !CALL SST_LinearMeridional(patch_3d, ocean_state%p_prog(nold(1))%tracer(:,:,:,1))
        !p_oce_sfc%data_surfRelax_Temp(:,:)=ocean_state%p_prog(nold(1))%tracer(:,1,:,1)
      !ENDIF
      ELSEIF(initial_temperature_type==203)THEN
        perturbation_lat = basin_center_lat + 0.1_wp * basin_height_deg
        perturbation_lon = basin_center_lon + 0.1_wp * basin_width_deg
        max_perturbation  = 0.1_wp!20.1_wp
        perturbation_width  = 10.0_wp!1.5_wp

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
          DO jc = start_cell_index, end_cell_index
            distan = SQRT((patch_2d%cells%center(jc,jb)%lat - perturbation_lat * deg2rad)**2 + &
            & (patch_2d%cells%center(jc,jb)%lon - perturbation_lon * deg2rad)**2)

            p_oce_sfc%data_surfRelax_Temp(jc,jb)=ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1)
            !Local cold perturbation
            IF(distan<=5.0_wp*deg2rad)THEN
              p_oce_sfc%data_surfRelax_Temp(jc,jb) =          &
              & ocean_state%p_prog(nold(1))%tracer(jc,1,jb,1)        &
              & - 2.0_wp*max_perturbation*EXP(-(distan/(perturbation_width*deg2rad))**2) &
              !                &   * sin(pi*v_base%zlev_m(jk)/4000.0_wp)!&
              & * SIN(pi*patch_3d%p_patch_1d(1)%zlev_m(1) / patch_3d%p_patch_1d(1)%zlev_i(2))
             ENDIF !Local cold perturbation
           END DO
         END DO
      ENDIF

    END SELECT

    IF (type_surfRelax_Salt == 3) THEN
      IF (no_tracer > 1) THEN
        p_oce_sfc%data_surfRelax_Salt(:,:) = ocean_state%p_prog(nold(1))%tracer(:,1,:,2)
      END IF
    END IF

    IF (type_3dimRelax_Temp >= 3) THEN
      ocean_nudge%data_3dimRelax_Temp(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,1)
      SELECT CASE (type_3dimRelax_Temp)
      CASE (4)
        ! 3D-relax the north and south boundary
        CALL init_3Drelax_coefficient_NS_boundaries( &
          & patch_3D = patch_3d, &
          & relax_coefficient=ocean_nudge%relax_3dim_coefficient, &
          & SouthBoundary=basin_SouthBoundary, &
          & NorthBoundary=basin_NorthBoundary, &
          & relaxWidth=relax_width * deg2rad)

      CASE (5)
        ! 3D-relax the north boundary (Abernathey)
        CALL init_3Drelax_coefficient_NS_boundaries( &
          & patch_3D = patch_3d, &
          & relax_coefficient=ocean_nudge%relax_3dim_coefficient, &
          & NorthBoundary=basin_NorthBoundary, &
          & relaxWidth=relax_width * deg2rad)

      CASE (6)
        ! as above 3D-relax the north boundary
        ! but with explicit relaxation temperature (Abernathey)
        CALL init_3Drelax_coefficient_NS_boundaries( &
          & patch_3D = patch_3d, &
          & relax_coefficient=ocean_nudge%relax_3dim_coefficient, &
          & NorthBoundary=basin_NorthBoundary, &
          & relaxWidth=relax_width * deg2rad)

        CALL tracer_ConstantSurface(patch_3d=patch_3d, ocean_tracer=ocean_nudge%data_3dimRelax_Temp, &
          & top_value=initial_temperature_north)

        CALL varyTracerVerticallyExponentially(patch_3d=patch_3d, &
          & ocean_tracer=ocean_nudge%data_3dimRelax_Temp,   &
          & bottom_value=initial_temperature_bottom,              &
          & scale_depth=1000.0_wp)

      END SELECT

    END IF
    IF (type_3dimRelax_Salt == 3) THEN
      IF (no_tracer > 1) THEN
        ocean_nudge%data_3dimRelax_Salt(:,:,:) = ocean_state%p_prog(nold(1))%tracer(:,:,:,2)
      ELSE
        CALL finish(TRIM(routine),' type_3dimRelax_Salt=3 and no_tracer<2 - ABORT')
      END IF
    END IF

    !---------Debug Diagnostics-------------------------------------------
    IF (type_surfRelax_Temp > 0) THEN
      idt_src=0  ! output print level - 0: print in any case
      ! z_c(:,1,:) = p_oce_sfc%data_surfRelax_Temp(:,:)
      CALL dbg_print('init relaxation - T'       ,p_oce_sfc%data_surfRelax_Temp(:,:), &
        & str_module,idt_src, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
    END IF
    IF (type_surfRelax_Salt > 0) THEN
      IF (no_tracer > 1) THEN
        z_c(:,1,:) = p_oce_sfc%data_surfRelax_Salt(:,:)
        CALL dbg_print('init relaxation - S'       ,z_c   ,str_module,idt_src, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
      ELSE
        CALL finish(TRIM(routine),' type_surfRelax_Salt>0 and no_tracer<2 - ABORT')
      END IF
    END IF

    CALL message( TRIM(routine),'end' )

  END SUBROUTINE init_ho_relaxation
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE  init_3Drelax_coefficient_NS_boundaries(patch_3D, relax_coefficient, SouthBoundary, NorthBoundary, relaxWidth)

    TYPE(t_patch_3d ),TARGET :: patch_3D
    REAL(wp) :: relax_coefficient(:,:,:)
    REAL(wp), OPTIONAL :: SouthBoundary, NorthBoundary
    REAL(wp) :: relaxWidth

    INTEGER :: jb, jc, start_cell_index, end_cell_index
    REAL(wp) :: lat_diff, south_boundary, north_boundary
    TYPE(t_patch),POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    !-------------------------------------------------------------------------
    patch_2d  => patch_3d%p_patch_2d(1)
    all_cells => patch_2d%cells%ALL

    IF (PRESENT(SouthBoundary)) THEN
      south_boundary = SouthBoundary
    ELSE
      south_boundary = -200.0_wp - relaxWidth
    ENDIF
    IF (PRESENT(NorthBoundary)) THEN
      north_boundary = NorthBoundary
    ELSE
      north_boundary = 200.0_wp + relaxWidth
    ENDIF

    relax_coefficient(:,:,:) = 0.0_wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_cell_index, end_cell_index)
      DO jc = start_cell_index, end_cell_index
!         write(0,*) start_cell_index, end_cell_index, jc, jb
        lat_diff = patch_2d%cells%center(jc,jb)%lat - south_boundary
        IF (lat_diff > relaxWidth) & ! check the north boundary
          lat_diff = north_boundary - patch_2d%cells%center(jc,jb)%lat

        IF (lat_diff >= 0.0_wp .AND. lat_diff < relaxWidth) THEN
          relax_coefficient(jc,:,jb) = (relaxWidth - lat_diff) / relaxWidth
        ENDIF

      END DO
    END DO

  END SUBROUTINE  init_3Drelax_coefficient_NS_boundaries
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !<Optimize:inUse>
  SUBROUTINE init_ocean_WindForcing(patch_3d, ocean_state, fu10)

    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET       :: ocean_state
    REAL(wp), INTENT(INOUT)                 :: fu10(:,:)      !  windspeed: p_as%fu10

    TYPE(t_patch), POINTER        :: patch_2d
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    REAL(wp):: windStress_u(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    REAL(wp):: windStress_v(nproma,patch_3D%p_patch_2D(1)%nblks_e)
    INTEGER :: start_edge_index, end_edge_index !,i_startblk_c, i_endblk_c,
    INTEGER :: je,blockNo, level

    patch_2d    => patch_3d%p_patch_2d(1)
    all_edges   => patch_2d%edges%all
    owned_edges => patch_2d%edges%owned

    CALL set_windstress(patch_2d, windStress_u, forcing_windstress_u_type, &
      & forcing_windStress_u_amplitude, forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno)

    CALL set_windstress(patch_2d, windStress_v, forcing_windstress_v_type, &
      & forcing_windStress_v_amplitude, forcing_windstress_zonal_waveno, forcing_windstress_merid_waveno)

    ! calculate normal edge stress from u,v
!ICON_OMP_PARALLEL_DO PRIVATE(start_edge_index,end_edge_index, je, level) ICON_OMP_DEFAULT_SCHEDULE
    DO blockNo = all_edges%start_block,all_edges%end_block
      CALL get_index_range(all_edges, blockNo, start_edge_index, end_edge_index)
      ocean_state%p_aux%bc_top_WindStress(:,blockNo) = 0.0_wp
      DO je =  start_edge_index, end_edge_index
        DO level=1,MAX(patch_3D%p_patch_1D(1)%dolic_e(je,blockNo),1)

          ocean_state%p_aux%bc_top_WindStress(je,blockNo) = &
            & windStress_u(je,blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v1 + &
            & windStress_v(je,blockNo) * patch_2d%edges%primal_normal(je,blockNo)%v2

          ENDDO
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO


!     field_2d(:,:) = MERGE(amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length),0.0_wp,mask(:,:) <= threshold)
!     CALL set_windspeed(all_cells, patch_3d%lsm_c(:,1,:), sea_boundary, fu10,                                      &
!       & atmos_fluxes%topBoundCond_windStress_u, atmos_fluxes%topBoundCond_windStress_v, forcing_windspeed_type,   &
!       & forcing_windspeed_amplitude)

    idt_src=0  ! output print level - 0: print in any case
    CALL dbg_print('init zonal wind stress'    ,windStress_u,str_module,idt_src,in_subset=owned_edges)
    CALL dbg_print('init merid. wind stress'   ,windStress_v,str_module,idt_src,in_subset=owned_edges)
!     CALL dbg_print('init wind speed'           ,fu10                                  ,str_module,idt_src,in_subset=owned_cells)

  END SUBROUTINE init_ocean_WindForcing
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE set_windstress(patch_2D, windstress, &
      &                     control, amplitude, zonal_waveno, meridional_waveno)
    TYPE(t_patch), POINTER           :: patch_2d
    REAL(wp), INTENT(INOUT)          :: windstress(:,:)
    INTEGER,  INTENT(IN)             :: control
    REAL(wp), INTENT(IN)             :: amplitude,zonal_waveno,meridional_waveno

!     REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp), TARGET                 :: lon(nproma,patch_2D%nblks_e), lat(nproma,patch_2D%nblks_e)
    REAL(wp)                         :: center, length

    lat(:,:) = patch_2D%edges%center(:,:)%lat
    lon(:,:) = patch_2D%edges%center(:,:)%lon
    center = basin_center_lat * deg2rad
    length = basin_height_deg * deg2rad

    SELECT CASE (control)
    CASE (0) ! NO FORCING, SET TO ZERO ========================================
      windstress = 0.0_wp
    CASE (1:100)      ! FILE INPUT, DONE ELSEWHERE ============================
      windstress = 0.0_wp
      CALL message('windstress forcing','file input')
    CASE (101:200)    ! ANALYTIC SETUP ========================================
      SELECT CASE (control)
      CASE(101)       ! constant amplitude
        windstress = amplitude
      CASE(102)       ! basin setup, zonally changed
        CALL basin_zonal(lon, lat, windstress,amplitude,zonal_waveno, length)
      CASE(103)       ! basin setup, meridionally changed
        CALL basin_meridional(lon, lat, windstress,amplitude,meridional_waveno, length)
      CASE(104)       ! zonally periodic, nonzero at pols, meridionally constant
        CALL zonal_periodic_nonzero_around_center_zero_at_pols(lon, lat,  windstress, amplitude)
      CASE(105)
        CALL meridional_periodic_around_center_zero_at_pols(lon, lat, windstress, amplitude)
      CASE(106)       ! zonally periodic around a given center, zero at pols, meridionally constant
        CALL zonal_periodic_zero_at_pols(lon, lat, windstress,amplitude,zonal_waveno)
      CASE(107)       ! latteral cells, zonal period only
        CALL cells_zonal_periodic(lon, lat, windstress,amplitude,zonal_waveno)
      CASE(108)       ! latteral cells, zonally and meridionally periodic
        CALL cells_zonal_and_meridional_periodic(lon, lat, windstress,amplitude,zonal_waveno,meridional_waveno)
      CASE(109)
        CALL cells_zonal_and_meridional_periodic_constant_amplitude_sin(lon, lat, windstress, amplitude)
      CASE(110)
        CALL cells_zonal_and_meridional_periodic_constant_amplitude_cosin(lon, lat, windstress, amplitude)
      CASE(111)
        CALL Wolfe_Cessi_TestCase(lon, lat, windstress, amplitude)
      CASE(112)       ! basin setup, zonally changed for Abernathey test case
        CALL basin_zonal_zeroOutside(lon, lat, windstress, amplitude, zonal_waveno, center, length)
      CASE(113)
        CALL zentral_jet(lon, lat, windstress, amplitude)
      CASE(114)
        CALL wind_shear_u(patch_2D, lon, lat, windstress, amplitude)
      CASE(115)
        CALL wind_shear_v(patch_2D, lon, lat, windstress, amplitude)
      END SELECT
    END SELECT

  END SUBROUTINE set_windstress

  SUBROUTINE basin_zonal(lon, lat, field_2d, amplitude, zonal_waveno, length)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN)             :: zonal_waveno, length

    field_2d(:,:) = amplitude * COS(zonal_waveno*pi*(lat(:,:)-length)/length)

  END SUBROUTINE basin_zonal

  SUBROUTINE basin_zonal_zeroOutside(lon, lat, field_2d, amplitude, zonal_waveno, center, length)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN)             :: zonal_waveno, center, length

    field_2d(:,:) = amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length)
    field_2d(:,:) = MERGE(field_2d(:,:),0.0_wp, lat(:,:) > center-0.5_wp*length)
    field_2d(:,:) = MERGE(field_2d(:,:),0.0_wp, lat(:,:) < center+0.5_wp*length)

  END SUBROUTINE basin_zonal_zeroOutside

  SUBROUTINE basin_meridional(lon, lat, field_2d, amplitude, meridional_waveno, length)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN)             :: meridional_waveno, length


    field_2d(:,:) = amplitude * COS(meridional_waveno*pi*(lon(:,:)-length)/length)

  END SUBROUTINE basin_meridional

  SUBROUTINE zonal_periodic_nonzero_around_center_zero_at_pols(lon, lat, field_2d, amplitude,&
      & center_opt,length_opt,zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: center_opt, length_opt, zonal_waveno_opt

    REAL(wp) :: center, length, zonal_waveno

    length       = 180.0_wp * deg2rad
    center       = -60.0_wp * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CaLL assign_if_present(zonal_waveno, zonal_waveno_opt)

    field_2d(:,:) = amplitude * COS(zonal_waveno*pi*(lat(:,:)-center)/length)
  END SUBROUTINE zonal_periodic_nonzero_around_center_zero_at_pols

  SUBROUTINE meridional_periodic_around_center_zero_at_pols(lon, lat, field_2d, amplitude, &
      & center_opt,length_opt,meridional_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: center_opt, length_opt, meridional_waveno_opt

    REAL(wp) :: center, length, meridional_waveno

    length                        =  90.0_wp * deg2rad
    center                        = -20.0_wp * deg2rad
    meridional_waveno             = forcing_windstress_merid_waveno

    CALL assign_if_present(center,center_opt)
    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(meridional_waveno, meridional_waveno_opt)

    field_2d(:,:) = amplitude * COS(meridional_waveno*pi*(lon(:,:)-center)/length)
  END SUBROUTINE meridional_periodic_around_center_zero_at_pols

  SUBROUTINE zonal_periodic_zero_at_pols(lon, lat, field_2d, amplitude, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt

    REAL(wp) :: zonal_waveno

    zonal_waveno = forcing_windstress_zonal_waveno

    CaLL assign_if_present(zonal_waveno, zonal_waveno_opt)

    field_2d(:,:) = amplitude * COS(lat(:,:)) * &
      & COS(forcing_windstress_zonalWavePhase + zonal_waveno * ABS(lat(:,:)))

  END SUBROUTINE zonal_periodic_zero_at_pols

  SUBROUTINE Wolfe_Cessi_TestCase(lon, lat, field_2d, amplitude)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    field_2d(:,:) = amplitude *                                   &
      & (0.8_wp*EXP(-1*lat(:,:)*lat(:,:)/0.01) - COS(6.0_wp*lat(:,:)))

  END SUBROUTINE Wolfe_Cessi_TestCase

  SUBROUTINE zentral_jet(lon, lat, field_2d, amplitude)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    field_2d(:,:) = amplitude *                                   &
      & (0.8_wp*EXP(-1*lat(:,:)*lat(:,:)/0.5) )

  END SUBROUTINE zentral_jet

  SUBROUTINE wind_shear_u(patch_2d, lon, lat, field_2d, amplitude)
    TYPE(t_patch), POINTER           :: patch_2d
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    INTEGER :: edge_index,edge_block, start_edges_index, end_edges_index
    TYPE(t_subset_range), POINTER :: all_edges
    REAL(wp) :: point_lon, point_lat, uu, vv

    all_edges => patch_2d%edges%ALL
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index

        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat

        IF(point_lat>=basin_center_lat)THEN
          uu=tanh((point_lat)*30.0_wp)
        ELSEIF (point_lat<basin_center_lat) THEN
          uu=tanh((-point_lat)*30.0_wp)
        ENDIF

        vv=0.1_wp*sin(2.0_wp*pi*point_lon)

        field_2d(edge_index,edge_block) =amplitude*uu * patch_2d%edges%primal_normal(edge_index,edge_block)%v1

        ENDDO
      ENDDO

  END SUBROUTINE wind_shear_u

  SUBROUTINE wind_shear_v(patch_2d, lon, lat, field_2d, amplitude)
    TYPE(t_patch), POINTER           :: patch_2d
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    INTEGER :: edge_index,edge_block, start_edges_index, end_edges_index
    TYPE(t_subset_range), POINTER :: all_edges
    REAL(wp) :: point_lon, point_lat, uu, vv

    all_edges => patch_2d%edges%ALL
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_edges_index, end_edges_index)
      DO edge_index = start_edges_index, end_edges_index

        point_lon = patch_2d%edges%center(edge_index,edge_block)%lon
        point_lat = patch_2d%edges%center(edge_index,edge_block)%lat
        vv=0.1_wp*sin(2.0_wp*pi*point_lon)
        field_2d(edge_index,edge_block) =amplitude*vv * patch_2d%edges%primal_normal(edge_index,edge_block)%v2

      ENDDO
    ENDDO



  END SUBROUTINE wind_shear_v


  SUBROUTINE cells_zonal_periodic(lon, lat, field_2d, amplitude, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt

    REAL(wp) :: zonal_waveno

    zonal_waveno = forcing_windstress_zonal_waveno
    CALL assign_if_present(zonal_waveno, zonal_waveno_opt)

    field_2d(:,:) = amplitude &
                & * COS(lat) &
                & * COS(zonal_waveno * lat) &
                & * COS(lon)
  END SUBROUTINE cells_zonal_periodic

  SUBROUTINE cells_zonal_and_meridional_periodic(lon, lat, field_2d, amplitude, &
      & zonal_waveno_opt,meridional_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp),INTENT(INOUT)           :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN),OPTIONAL    :: zonal_waveno_opt, meridional_waveno_opt

    REAL(wp) :: zonal_waveno, meridional_waveno

    zonal_waveno      = forcing_windstress_zonal_waveno
    meridional_waveno = forcing_windstress_merid_waveno

    CALL assign_if_present(zonal_waveno, zonal_waveno_opt)
    CALL assign_if_present(meridional_waveno, meridional_waveno_opt)

    field_2d(:,:) = amplitude &
      & * COS(lat(:,:)) &
      & * COS(zonal_waveno * lat(:,:)) &
      & * SIN(meridional_waveno * lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic

  SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_sin(lon, lat, field_2d, &
      & amplitude,length_opt, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp)                         :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno

    length       = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    ! -length/length ???
!     strength = amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)

    field_2d(:,:) = amplitude*amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)*sin(lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_sin

  SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_cosin(lon, lat, field_2d, &
      & amplitude,length_opt, zonal_waveno_opt)
    REAL(wp)                         :: lon(:,:), lat(:,:)
    REAL(wp)                         :: field_2d(:,:)
    REAL(wp), INTENT(IN)             :: amplitude
    REAL(wp), INTENT(IN) , OPTIONAL  :: length_opt,zonal_waveno_opt

    REAL(wp) :: length, zonal_waveno

    length       = basin_height_deg * deg2rad
    zonal_waveno = forcing_windstress_zonal_waveno

    CALL assign_if_present(length,length_opt)
    CALL assign_if_present(zonal_waveno,zonal_waveno_opt)

    ! -length/length ???
!     strength = amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)

    field_2d(:,:) = amplitude*amplitude*cos(zonal_waveno*pi*lat(:,:)-length/length)*cos(lon(:,:))

  END SUBROUTINE cells_zonal_and_meridional_periodic_constant_amplitude_cosin
  !-------------------------------------------------------------------------

!<Optimize:inUse>
  SUBROUTINE set_windspeed( windspeed, windstress_u, windstress_v, control, amplitude)
    REAL(wp), INTENT(INOUT)          :: windspeed(:,:)
    REAL(wp), INTENT(INOUT)          :: windstress_u(:,:)
    REAL(wp), INTENT(INOUT)          :: windstress_v(:,:)
    INTEGER,  INTENT(IN)             :: control
    REAL(wp), INTENT(IN)             :: amplitude

    SELECT CASE (control)
    CASE (0) ! NO FORCING, SET TO ZERO ========================================
      windspeed = 0.0_wp
    CASE (1:100)      ! FILE INPUT, DONE ELSEWHERE ============================
      CALL message('windspeed forcing','file input')
    CASE (101:200)    ! ANALYTIC SETUP ========================================
      SELECT CASE (control)
      CASE(101)       ! constant amplitude
        windspeed(:,:) = amplitude
      CASE(102)       ! windspeed backward calculated from windstress amplitude
        CALL calc_windspeed_fromwindstress( windspeed, windstress_u, windstress_v, amplitude)
   !  CASE(103)       ! basin setup, zonally changed
   !    CALL basin_zonal(subset,mask,threshold,windstress,amplitude,length)
      END SELECT
    END SELECT

  END SUBROUTINE set_windspeed
  !-------------------------------------------------------------------------
  !>
  !! Initialization of wind speed if wind stress is available
  !!  Following common relation of wind stress to wind speed (e.g. Smith 1980)
  !!
  !
  !
  !-------------------------------------------------------------------------
  !
  SUBROUTINE calc_windspeed_fromwindstress( windspeed, windstress_u, windstress_v, amplitude)
  !!
  !!
    REAL(wp)                         :: windspeed(:,:)
    REAL(wp)                         :: windstress_u(:,:)
    REAL(wp)                         :: windstress_v(:,:)
    REAL(wp), INTENT(IN)             :: amplitude

    ! |tau|  = Rho_air * C_D * U_10**2
    ! |U_10| = ( |tau| / (Rho_air*C_D) )**(1/2) or
    ! |U_10| = (Rho_air*C_D)**(-1/2) * |tau|**(1/2)
    !  with Rho_air=1.22 kg/m3, drag coefficient C_D ~ 1.3e-3, U_10 the wind speed in 10m height
    !  (Rho_air * C_D)**(-1/2) ~ (1.5 e-3 )**(-1/2) ~ 26
    !  or: |U_10| ~ 26 m/s for tau=1 Pa, or 8.2 m/s for tau=0.1 Pa
    !  |tau| = (tau_x**2 + tau_y**2)**(1/2)

    windspeed(:,:) = amplitude*26.0_wp*                                                                 &
      &              SQRT(SQRT(windstress_u(:,:)*windstress_u(:,:) + windstress_v(:,:)*windstress_v(:,:)))

  END SUBROUTINE calc_windspeed_fromwindstress

END MODULE mo_ocean_forcing
