!
! Interface between NWP physics and the ocean, through a coupler
!
! Notes on openMP parallelisation:
!   Most jb-loops cannot be optimized easily because converting a 2-D field into a 1-D field
!   requires remembering the 1-D index "ncount", which is incremented over the loop.  This
!   technique is used in the current code and openMP is not used.
!   A solution suggested by Rene Redler would be to calculate that index (nn) would look like this:
!      !ICON_OMP_PARALLEL_DO PRIVATE(jb, ic, jc, nn) ICON_OMP_RUNTIME_SCHEDULE
!      DO jb = i_startblk, i_endblk
!        nn = (jb-1)*nproma                               ! translation to 1-d buffer fields
!        DO ic = 1, ext_data%atm%list_sea%ncount(jb)      ! number of ocean points (open water & sea ice)
!          jc = ext_data%atm%list_sea%idx(ic,jb)
!          prm_field(jg)%ocv(n,i_blk) = buffer(nn+jc,1)
!  It might also be necessary to synch the data before each loop passing data to the ocean.
!      CALL sync_patch_array(sync_c, p_patch, prm_diag%swflxsfc_t (:,:,isub_water) )
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

MODULE mo_nwp_ocean_coupling

  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2vmr
  USE mo_ccycle_config       ,ONLY: ccycle_config,                          &
       & CCYCLE_MODE_NONE, CCYCLE_MODE_PRESCRIBED, CCYCLE_MODE_INTERACTIVE, &
       & CCYCLE_CO2CONC_CONST, CCYCLE_CO2CONC_FROMFILE
  USE mo_dbg_nml             ,ONLY: idbg_mxmn, idbg_val
  USE mo_exception           ,ONLY: finish
  USE mo_ext_data_types      ,ONLY: t_external_data
  USE mo_fortran_tools       ,ONLY: assert_acc_host_only, init
  USE mo_idx_list            ,ONLY: t_idx_list_blocked
  USE mo_impl_constants      ,ONLY: start_prog_cells, end_prog_cells
  USE mo_kind                ,ONLY: wp
  USE mo_lnd_nwp_config      ,ONLY: isub_water, isub_seaice, hice_max
  USE mo_loopindices         ,ONLY: get_indices_c
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nonhydro_types      ,ONLY: t_nh_diag
  USE mo_nwp_phy_types       ,ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types       ,ONLY: t_wtr_prog, t_lnd_diag
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_physical_constants  ,ONLY: vmr_to_mmr_co2, tmelt
  USE mo_util_dbg_prnt       ,ONLY: dbg_print

  USE mo_atmo_ocean_coupling ,ONLY: &
    & field_id_co2_flx, field_id_co2_vmr, field_id_freshflx, &
    & field_id_heatflx, field_id_oce_u, field_id_oce_v, field_id_pres_msl, &
    & field_id_seaice_atm, field_id_seaice_oce, field_id_sp10m, field_id_sst, &
    & field_id_umfl, field_id_vmfl
  USE mo_coupling_utils      ,ONLY: cpl_put_field, cpl_get_field

  IMPLICIT NONE

  PRIVATE

  !>
  !! Package of fields passed to the ocean via YAC.
  TYPE t_nwp_ocean_fields_tx

    !> Sea-water fraction [m2(ocean)/m2(gridcell)].
    REAL(wp), CONTIGUOUS, POINTER :: frac_w(:,:) => NULL()
    !> Sea-ice fraction [m2(seaice)/m2(gridcell)].
    REAL(wp), CONTIGUOUS, POINTER :: frac_i(:,:) => NULL()

    !> Resolved zonal surface stress over sea water [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: umfl_s_w(:,:) => NULL()
    !> Resolved zonal surface stress over sea ice [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: umfl_s_i(:,:) => NULL()
    !> Resolved meridional surface stress over sea water [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: vmfl_s_w(:,:) => NULL()
    !> Resolved meridional surface stress over sea ice [N/m2].
    REAL(wp), CONTIGUOUS, POINTER :: vmfl_s_i(:,:) => NULL()

    !> Moisture flux at sea water surface [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: qhfl_s_w(:,:) => NULL()
    !> Moisture flux at sea ice surface [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: qhfl_s_i(:,:) => NULL()

    !> Sensible heat flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: shfl_s_w(:,:) => NULL()
    !> Sensible heat flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: shfl_s_i(:,:) => NULL()

    !> Latent heat flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lhfl_s_w(:,:) => NULL()
    !> Latent heat flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lhfl_s_i(:,:) => NULL()

    !> Conductive heat flux at water-ice interface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: chfl_i(:,:) => NULL()
    !> Melt potential at atmosphere-ice interface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: meltpot_i(:,:) => NULL()

    !> 10m wind speed [m/s].
    REAL(wp), CONTIGUOUS, POINTER :: sp_10m(:,:) => NULL()

    !> Surface pressure [Pa].
    REAL(wp), CONTIGUOUS, POINTER :: pres_sfc(:,:) => NULL()

    !> Net short-wave flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: swflxsfc_w(:,:) => NULL()
    !> Net short-wave flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: swflxsfc_i(:,:) => NULL()

    !> Net long-wave flux at sea water surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lwflxsfc_w(:,:) => NULL()
    !> Net long-wave flux at sea ice surface [W/m2].
    REAL(wp), CONTIGUOUS, POINTER :: lwflxsfc_i(:,:) => NULL()

    !> Total rain rate [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: rain_rate(:,:) => NULL()
    !> Total snow rate [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: snow_rate(:,:) => NULL()

    !> Surface CO2 concentration [kg(CO2)/kg(air)].
    !! Not contiguous because it is the lowest level of a 3D field.
    REAL(wp), POINTER :: q_co2(:,:) => NULL()

  END TYPE t_nwp_ocean_fields_tx


  !>
  !! Package of fields received from the ocean via YAC.
  TYPE t_nwp_ocean_fields_rx

    !> Sea-surface temperature [K].
    REAL(wp), CONTIGUOUS, POINTER :: t_seasfc(:,:) => NULL()

    !> Sea-ice fraction [m2(seaice)/m2(ocean)].
    REAL(wp), CONTIGUOUS, POINTER :: fr_seaice(:,:) => NULL()

    !> Sea-ice thickness [m].
    REAL(wp), CONTIGUOUS, POINTER :: h_ice(:,:) => NULL()

    !> Zonal ocean surface velocity (optional) [m/s].
    REAL(wp), CONTIGUOUS, POINTER :: ocean_u(:,:) => NULL()

    !> Meridional ocean surface velocity (optional) [m/s].
    REAL(wp), CONTIGUOUS, POINTER :: ocean_v(:,:) => NULL()

    !> CO2 surface flux [kg/m2/s].
    REAL(wp), CONTIGUOUS, POINTER :: flx_co2(:,:) => NULL()

  END TYPE t_nwp_ocean_fields_rx


  PUBLIC :: nwp_couple_ocean
  PUBLIC :: couple_ocean
  PUBLIC :: t_nwp_ocean_fields_rx
  PUBLIC :: t_nwp_ocean_fields_tx

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_nwp_ocean_coupling' ! Output of module for debug

CONTAINS

  !>
  !! SUBROUTINE nwp_couple_ocean -- the interface between
  !! NWP physics and the ocean, through a coupler
  !!
  !! This subroutine is called from nwp_nh_interface.
  SUBROUTINE nwp_couple_ocean ( &
      & p_patch, pt_diag, lnd_diag, wtr_prog_new, prm_diag, ext_data, lacc &
    )

    TYPE(t_patch),                INTENT(INOUT)  :: p_patch
    TYPE(t_nh_diag),      TARGET, INTENT(INOUT)  :: pt_diag
    TYPE(t_wtr_prog),     TARGET, INTENT(INOUT)  :: wtr_prog_new
    TYPE(t_lnd_diag),     TARGET, INTENT(INOUT)  :: lnd_diag
    TYPE(t_nwp_phy_diag), TARGET, INTENT(INOUT)  :: prm_diag
    TYPE(t_external_data),        INTENT(IN)     :: ext_data
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    TYPE(t_nwp_ocean_fields_tx) :: tx
    TYPE(t_nwp_ocean_fields_rx) :: rx

    REAL(wp), TARGET :: rain_rate(nproma, p_patch%nblks_c)
    REAL(wp), TARGET :: snow_rate(nproma, p_patch%nblks_c)

    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    INTEGER :: jb, jc

    LOGICAL :: have_ice, have_hail, have_graupel

    CALL assert_acc_host_only('nwp_couple_ocean', lacc)

    ! include boundary interpolation zone of nested domains and halo points
    i_startblk = p_patch%cells%start_block(start_prog_cells)
    i_endblk   = p_patch%cells%end_block(end_prog_cells)

    have_ice = ASSOCIATED(prm_diag%ice_gsp_rate)
    have_hail = ASSOCIATED(prm_diag%hail_gsp_rate)
    have_graupel = ASSOCIATED(prm_diag%graupel_gsp_rate)

    ! Prepare rain and snow rates.
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, start_prog_cells, end_prog_cells)

      DO jc = i_startidx, i_endidx
        rain_rate(jc,jb) = prm_diag%rain_con_rate_corr(jc,jb) + prm_diag%rain_gsp_rate(jc,jb)
        snow_rate(jc,jb) = prm_diag%snow_con_rate_corr(jc,jb) + prm_diag%snow_gsp_rate(jc,jb)

        IF (have_ice) snow_rate(jc,jb) = snow_rate(jc,jb) + prm_diag%ice_gsp_rate(jc,jb)
        IF (have_hail) snow_rate(jc,jb) = snow_rate(jc,jb) + prm_diag%hail_gsp_rate(jc,jb)
        IF (have_graupel) snow_rate(jc,jb) = snow_rate(jc,jb) + prm_diag%graupel_gsp_rate(jc,jb)
      END DO
    END DO

    ! Send fields:
    ! 1. prm_diag%umfl_s_t(:,:,:)                       zonal resolved surface stress and  [N/m2]
    !    prm_diag%vmfl_s_t(:,:,:)                       meridional resolved surface stress [N/m2]
    !    prm_diag%sp_10m  (:,:)                         10m wind speed [m/s]
    !    pt_diag%pres_sfc (:,:)                         surface pressure [Pa]
    !
    ! 2. prm_diag%rain_con_rate_corr(:,:)               convective surface rain rate    [kg/m2/s]
    !    prm_diag%rain_gsp_rate(:,:)                    grid-scale surface rain rate    [kg/m2/s]
    !
    !    prm_diag%snow_con_rate_corr(:,:)               convective surface snow_rate    [kg/m2/s]
    !    prm_diag%snow_gsp_rate    (:,:)                grid_scale surface snow rate    [kg/m2/s]
    !    prm_diag%ice_gsp_rate     (:,:)                grid_scale surface ice rate     [kg/m2/s]
    !    prm_diag%graupel_gsp_rate (:,:)                grid_scale surface graupel rate [kg/m2/s]
    !    prm_diag%hail_gsp_rate    (:,:)                grid_scale surface hail rate    [kg/m2/s]
    !
    !    prm_diag%qhfl_s_t(:,:,isub_water/isub_seaice)  moisture flux (surface) aka evaporation rate at surface [kg/m2/s]
    !
    ! 3. prm_diag%swflxsfc_t (:,:,:)                    tile-based shortwave net flux at surface [W/m2]
    !    prm_diag%lwflxsfc_t (:,:,:)                    tile-based longwave net flux at surface  [W/m2]
    !    prm_diag%shfl_s_t   (:,:,:)                    tile-based sensible heat flux at surface [W/m2]
    !    prm_diag%lhfl_s_t   (:,:,:)                    tile-based latent   heat flux at surface [W/m2]
    !    lnd_diag%condhf_ice (:,:)                      conductive heat flux at ice-ocean interface [W/m2]

    tx%frac_w => ext_data%atm%frac_t(:,:,isub_water)
    tx%frac_i => ext_data%atm%frac_t(:,:,isub_seaice)

    tx%umfl_s_w => prm_diag%umfl_s_t(:,:,isub_water)
    tx%umfl_s_i => prm_diag%umfl_s_t(:,:,isub_seaice)
    tx%vmfl_s_w => prm_diag%vmfl_s_t(:,:,isub_water)
    tx%vmfl_s_i => prm_diag%vmfl_s_t(:,:,isub_seaice)

    tx%qhfl_s_w => prm_diag%qhfl_s_t(:,:,isub_water)
    tx%qhfl_s_i => prm_diag%qhfl_s_t(:,:,isub_seaice)

    tx%shfl_s_w => prm_diag%shfl_s_t(:,:,isub_water)
    tx%shfl_s_i => prm_diag%shfl_s_t(:,:,isub_seaice)

    tx%lhfl_s_w => prm_diag%lhfl_s_t(:,:,isub_water)
    tx%lhfl_s_i => prm_diag%lhfl_s_t(:,:,isub_seaice)

    tx%chfl_i => lnd_diag%condhf_ice(:,:)
    tx%meltpot_i => lnd_diag%meltpot_ice(:,:)

    tx%sp_10m => prm_diag%sp_10m(:,:)
    tx%pres_sfc => pt_diag%pres_sfc(:,:)

    tx%swflxsfc_w => prm_diag%swflxsfc_t(:,:,isub_water)
    tx%swflxsfc_i => prm_diag%swflxsfc_t(:,:,isub_seaice)

    tx%lwflxsfc_w => prm_diag%lwflxsfc_t(:,:,isub_water)
    tx%lwflxsfc_i => prm_diag%lwflxsfc_t(:,:,isub_seaice)

    tx%rain_rate => rain_rate(:,:)
    tx%snow_rate => snow_rate(:,:)

    tx%q_co2 => NULL()

    ! Receive fields
    ! 1. lnd_diag%t_seasfc (:,:)   SST [K]
    ! 2. lnd_diag%fr_seaice(:,:)   Sea-ice fraction [m2(ice)/m2(ocean)]
    ! 3. wtr_prog_new%h_ice(:,:)   Sea-ice height [m]

    rx%t_seasfc => lnd_diag%t_seasfc(:,:)
    rx%fr_seaice => lnd_diag%fr_seaice(:,:)
    rx%h_ice => wtr_prog_new%h_ice(:,:)
    rx%ocean_u => NULL()
    rx%ocean_v => NULL()
    rx%flx_co2 => NULL()

    CALL couple_ocean(p_patch, ext_data%atm%list_sea, tx, rx, lacc)

  END SUBROUTINE nwp_couple_ocean


  !>
  !! SUBROUTINE couple_ocean -- the actual interface between
  !! NWP physics and the ocean, through a coupler
  !!
  SUBROUTINE couple_ocean( p_patch, list_sea, tx, rx, lacc )

    ! Arguments

    TYPE(t_patch), INTENT(IN) :: p_patch
    TYPE(t_idx_list_blocked), INTENT(IN) :: list_sea
    TYPE(t_nwp_ocean_fields_tx), INTENT(IN) :: tx
    TYPE(t_nwp_ocean_fields_rx), INTENT(INOUT) :: rx
    LOGICAL, INTENT(IN), OPTIONAL :: lacc ! If true, use openacc

    ! Local variables

    INTEGER               :: jg                    ! grid index
    INTEGER               :: jb                    ! block loop count
    INTEGER               :: jc                    ! nproma loop count
    INTEGER               :: ic                    ! nproma loop count
    INTEGER               :: i_startblk, i_endblk  ! blocks
    INTEGER               :: i_startidx, i_endidx  ! slices

    REAL(wp), TARGET      :: buf(nproma, p_patch%nblks_c)

    REAL(wp), PARAMETER   :: csmall = 1.0E-5_wp    ! small number (security constant)

    LOGICAL :: received_data

    REAL(wp) :: co2conc

    CHARACTER(LEN=*), PARAMETER :: routine = str_module // ':couple_ocean'

    CALL assert_acc_host_only('couple_ocean', lacc)

    ! As YAC does not touch masked data an explicit initialisation with zero
    ! is required as some compilers are asked to initialise with NaN
    ! and as we loop over the full array.

    !$OMP PARALLEL
    CALL init(buf(:,:), lacc=.FALSE.)
    !$OMP END PARALLEL

    jg = p_patch%id

    i_startblk = p_patch%cells%start_block(start_prog_cells)
    i_endblk   = p_patch%cells%end_block(end_prog_cells)

    !-------------------------------------------------------------------------
    ! Send fields to ocean:
    !  "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !  "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !  "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !  "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !  "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !  "10m_wind_speed"                           - atmospheric wind speed
    !  "qtrc(nlev,co2)"                           - co2 mixing ratio
    !  "pres_msl"                                 - sea level pressure
    !
    ! Receive fields from ocean:
    !  "sea_surface_temperature"                  - SST
    !  "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !  "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !  "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !  "co2_flux"                                 - ocean co2 flux
    !-------------------------------------------------------------------------

    IF (.NOT. (SIZE(tx%umfl_s_w, 1) == nproma .AND. &
               SIZE(tx%umfl_s_w, 2) >= p_patch%nblks_c)) THEN
      CALL finish ('couple_ocean', 'first field extent must be nproma and &
          &field size has to be at least nproma*nblocks_c')
    END IF


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to ocean
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    !------------------------------------------------
    !  Send zonal wind stress bundle
    !    "surface_downward_eastward_stress" bundle
    !    - zonal wind stress component over ice and water
    !------------------------------------------------

    CALL cpl_put_field( &
      routine, field_id_umfl, 'u-stress', p_patch%n_patch_cells, &
      field_1=tx%umfl_s_w, field_2=tx%umfl_s_i)

    !------------------------------------------------
    !  Send meridional wind stress bundle
    !    "surface_downward_northward_stress" bundle
    !    - meridional wind stress component over ice and water
    !------------------------------------------------

    CALL cpl_put_field( &
      routine, field_id_vmfl, 'v-stress', p_patch%n_patch_cells, &
      field_1=tx%vmfl_s_w, field_2=tx%vmfl_s_i)

    !------------------------------------------------
    !  Send surface fresh water flux bundle
    !    "surface_fresh_water_flux" bundle
    !    - liquid rain, snowfall, evaporation
    !
    !    Note: the evap_tile should be properly updated and added;
    !          as long as evaporation over sea-ice is not used in ocean thermodynamics, the evaporation over the
    !          whole ocean part of grid-cell is passed to the ocean
    !          for pre04 a preliminary solution for evaporation in ocean model is to exclude the land fraction
    !          evap.oce = (evap.wtr*frac.wtr + evap.ice*frac.ice)/(1-frac.lnd)
    !------------------------------------------------

    DO jb = i_startblk, i_endblk
      CALL get_indices_c ( &
          & p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
          & start_prog_cells, end_prog_cells &
        )

      DO jc = i_startidx, i_endidx
        ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
        IF (tx%frac_w(jc,jb) + tx%frac_i(jc,jb) <= 0.0_wp) THEN
          ! ocean part is zero
          buf(jc,jb) = 0.0_wp
        ELSE
          buf(jc,jb) = tx%qhfl_s_w(jc,jb) * tx%frac_w(jc,jb) / (tx%frac_w(jc,jb) + tx%frac_i(jc,jb)) + &
            &          tx%qhfl_s_i(jc,jb) * tx%frac_i(jc,jb) / (tx%frac_w(jc,jb) + tx%frac_i(jc,jb))
        ENDIF
      ENDDO
    ENDDO

    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN
      CALL dbg_print('NWPOce: evapo-cpl', buf, str_module, 3, in_subset=p_patch%cells%owned)
    END IF

    CALL cpl_put_field( &
      routine, field_id_freshflx, 'fresh water flux', p_patch%n_patch_cells, &
      field_1=tx%rain_rate, field_2=tx%snow_rate, field_3=buf)

    !------------------------------------------------
    !  Send total heat flux bundle
    !    "total heat flux" bundle
    !    - short wave, long wave, sensible, latent heat flux
    !------------------------------------------------

    CALL cpl_put_field( &
      routine, field_id_heatflx, 'heat flux', p_patch%n_patch_cells, &
      field_1=tx%swflxsfc_w, field_2=tx%lwflxsfc_w, &
      field_3=tx%shfl_s_w, field_4=tx%lhfl_s_w)

    !------------------------------------------------
    !  Send sea ice flux bundle
    !    "atmosphere_sea_ice_bundle"
    !    - sea ice surface and bottom melt potentials Qtop, Qbot (conductive heat flux)
    !------------------------------------------------

    CALL cpl_put_field( &
      routine, field_id_seaice_atm, 'atmos sea ice', p_patch%n_patch_cells, &
      field_1=tx%meltpot_i, field_2=tx%chfl_i)

    !------------------------------------------------
    !  Send 10m wind speed
    !    "10m_wind_speed"
    !    - atmospheric wind speed
    !------------------------------------------------

    CALL cpl_put_field( &
      routine, field_id_sp10m, 'wind speed', p_patch%n_patch_cells, &
      tx%sp_10m)

    !------------------------------------------------
    !  Send sea level pressure
    !    "pres_msl"
    !    - atmospheric sea level pressure
    !    - pres_sfc is used insted of pres_msl because
    !      * it is available at each fast physics timestep
    !      * it calculated the hydrostatic surface pressure with less noise
    !------------------------------------------------

    CALL cpl_put_field( &
      routine, field_id_pres_msl, 'sea level pressure', p_patch%n_patch_cells, &
      tx%pres_sfc)

    !------------------------------------------------
    !  Send co2 mixing ratio
    !    "co2_mixing_ratio"
    !    - CO2 mixing ratio in ppmv
    !------------------------------------------------

#ifndef __NO_ICON_OCEAN__
    IF (ccycle_config(jg)%iccycle /= CCYCLE_MODE_NONE) THEN

      IF (ccycle_config(jg)%iccycle == CCYCLE_MODE_INTERACTIVE .AND. .NOT. ASSOCIATED(tx%q_co2)) THEN
        CALL finish('nwp_couple_ocean', 'Interactive carbon cycle is active but no CO2 concentration passed.')
      END IF

      SELECT CASE (ccycle_config(jg)%iccycle)

      CASE (CCYCLE_MODE_INTERACTIVE)
        DO jb = i_startblk, i_endblk
          CALL get_indices_c ( &
              & p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
              & start_prog_cells, end_prog_cells &
            )
          DO jc = i_startidx, i_endidx
            buf(jc,jb) = 1.0e6_wp * tx%q_co2(jc,jb) / vmr_to_mmr_co2
          END DO
        END DO

      CASE (CCYCLE_MODE_PRESCRIBED)
        SELECT CASE (ccycle_config(jg)%ico2conc)
        CASE (CCYCLE_CO2CONC_CONST)
          co2conc = 1e6_wp * ccycle_config(jg)%vmr_co2
        CASE (CCYCLE_CO2CONC_FROMFILE)
          co2conc = 1e6_wp * ghg_co2vmr
        CASE DEFAULT
          co2conc = 0._wp
          CALL finish('nwp_couple_ocean', 'unknown value of ccycle_config(jg)%ico2conc')
        END SELECT

        DO jb = i_startblk, i_endblk
          CALL get_indices_c ( &
              & p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, &
              & start_prog_cells, end_prog_cells &
            )
          DO jc = i_startidx, i_endidx
            buf(jc,jb) = co2conc
          END DO
        END DO

      END SELECT

      CALL cpl_put_field( &
        routine, field_id_co2_vmr, 'co2 vmr', p_patch%n_patch_cells, buf)

    ENDIF
#endif /* ifndef __NO_ICON_OCEAN__ */


    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from ocean to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !  Receive fields, only assign values if something was received
    !   - ocean fields have undefined values on land, which are not sent to the atmosphere,
    !     therefore buffer is set to zero to avoid unintended usage of ocean values over land

    !$OMP PARALLEL
    CALL init(buf(:,:), lacc=.FALSE.)
    !$OMP END PARALLEL

    !------------------------------------------------
    !  Receive SST
    !    "sea_surface_temperature"
    !    - SST
    !------------------------------------------------

    CALL cpl_get_field( &
      routine, field_id_sst, 'sst', p_patch%n_patch_cells, &
      rx%t_seasfc, first_get=.TRUE.)

    !------------------------------------------------
    !  Receive zonal velocity
    !    "eastward_sea_water_velocity"
    !    - zonal velocity, u component of ocean surface current
    !    RR: not used in NWP so far, not activated for exchange in coupling.xml
    !------------------------------------------------

    IF (ASSOCIATED(rx%ocean_u)) &
      CALL cpl_get_field( &
        routine, field_id_oce_u, 'u velocity', p_patch%n_patch_cells, &
        rx%ocean_u)

    !------------------------------------------------
    !  Receive meridional velocity
    !    "northward_sea_water_velocity"
    !    - meridional velocity, v component of ocean surface current
    !    RR: not used in NWP so far, not activated for exchange in
    !        YAC configuration file
    !------------------------------------------------

    IF (ASSOCIATED(rx%ocean_v)) &
      CALL cpl_get_field( &
        routine, field_id_oce_v, 'v velocity', p_patch%n_patch_cells, &
        rx%ocean_v)

    !------------------------------------------------
    !  Receive sea ice bundle
    !    "ocean_sea_ice_bundle"
    !    - ice thickness, snow thickness, ice concentration
    !------------------------------------------------

    CALL cpl_get_field( &
      routine, field_id_seaice_oce, 'sea ice', p_patch%n_patch_cells, &
      field_1=rx%h_ice, field_2=buf, field_3=rx%fr_seaice, &
      received_data=received_data)

    IF (received_data) THEN

      ! --- Here we loop only over ocean points, because h_ice is used by both oceans and lakes.

      DO jb = i_startblk, i_endblk
!$NEC ivdep
        DO ic = 1, list_sea%ncount(jb)                   ! number of ocean points (open water & sea ice)
          jc = list_sea%idx(ic,jb)                       ! index list of ocean points
          ! --- Limiting h_ice from ocean to hice_max used in atmospheric sea-ice.  Note that the ocean uses
          !     seaice_limit*dzlev_m(1), a fractional thickness of the ocean top layer.  The user is responsible
          ! --- for consistency of these three namelist parameters.  There is no ocean/atmo consistency check in ICON.
          rx%h_ice(jc,jb) = MIN (rx%h_ice(jc,jb), hice_max-csmall)
        END DO

      ENDDO

    END IF


    !------------------------------------------------
    !  Receive co2 flux
    !    "co2_flux"
    !    - ocean co2 flux
    !------------------------------------------------

    IF (ccycle_config(jg)%iccycle /= CCYCLE_MODE_NONE .AND. &
        ASSOCIATED(rx%flx_co2)) &
      CALL cpl_get_field( &
        routine, field_id_co2_flx, 'CO2 flux', p_patch%n_patch_cells, &
        rx%flx_co2)

    !------------------------------------------------
    ! Debug outputs
    !------------------------------------------------

    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN
      ! u/v-stress on ice and water
      CALL dbg_print('NWPOce: u_stress wtr', tx%umfl_s_w(:,:),   str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: u_stress ice', tx%umfl_s_i(:,:),   str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: v_stress wtr', tx%vmfl_s_w(:,:),   str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: v_stress ice', tx%vmfl_s_i(:,:),   str_module, 4, in_subset=p_patch%cells%owned)

      CALL dbg_print('NWPOce: sp_10m      ', tx%sp_10m(:,:),     str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: pres_msl    ', tx%pres_sfc(:,:),   str_module, 3, in_subset=p_patch%cells%owned)

      CALL dbg_print('NWPOce: rain_rate   ', tx%rain_rate(:,:),  str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: snow_rate   ', tx%snow_rate(:,:),  str_module, 3, in_subset=p_patch%cells%owned)

      IF (ASSOCIATED(tx%q_co2)) THEN
        CALL dbg_print('NWPOce: q_co2       ', tx%q_co2(:,:),    str_module, 3, in_subset=p_patch%cells%owned)
      END IF

      buf(:,:) = tx%qhfl_s_w(:,:) * tx%frac_w(:,:) + tx%qhfl_s_i(:,:) * tx%frac_i(:,:)
      CALL dbg_print('NWPOce: evaporation', buf(:,:),            str_module, 2, in_subset=p_patch%cells%owned)

      buf(:,:) = tx%swflxsfc_w(:,:) + tx%lwflxsfc_w(:,:) + tx%shfl_s_w(:,:) + tx%lhfl_s_w(:,:)
      CALL dbg_print('NWPOce: totalhfx.wtr', buf(:,:),           str_module, 2, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: swflxsfc.wtr', tx%swflxsfc_w(:,:), str_module, 3, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: lwflxsfc.wtr', tx%lwflxsfc_w(:,:), str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: shflx.wtr   ', tx%shfl_s_w(:,:),   str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: lhflx.wtr   ', tx%lhfl_s_w(:,:),   str_module, 4, in_subset=p_patch%cells%owned)

      ! Qtop and Qbot
      buf(:,:) = tx%shfl_s_i(:,:) + tx%swflxsfc_i + tx%lhfl_s_i(:,:) + tx%lwflxsfc_i(:,:)
      CALL dbg_print('NWPOce: ice-Qtop    ', buf,                str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: ice-Qbot    ', tx%chfl_i(:,:),     str_module, 3, in_subset=p_patch%cells%owned)

      ! SST, sea ice, ocean velocity received
      CALL dbg_print('NWPOce: t_seasfc    ', rx%t_seasfc(:,:),   str_module, 2, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: h_ice       ', rx%h_ice(:,:),      str_module, 4, in_subset=p_patch%cells%owned)

      IF (ASSOCIATED(rx%ocean_u) .AND. ASSOCIATED(rx%ocean_v)) THEN
        CALL dbg_print('NWPOce: ocu         ', rx%ocean_u(:,:),  str_module, 3, in_subset=p_patch%cells%owned)
        CALL dbg_print('NWPOce: ocv         ', rx%ocean_v(:,:),  str_module, 4, in_subset=p_patch%cells%owned)
      END IF

      IF (ASSOCIATED(rx%flx_co2)) THEN
        CALL dbg_print('NWPOce: flx_co2     ', rx%flx_co2(:,:),  str_module, 3, in_subset=p_patch%cells%owned)
      END IF

      ! Fraction of tiles:
      CALL dbg_print('NWPOce: fr_seaice   ', rx%fr_seaice(:,:),  str_module, 2, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: frac.si     ', tx%frac_i(:,:),     str_module, 4, in_subset=p_patch%cells%owned)
      CALL dbg_print('NWPOce: frac.wtr    ', tx%frac_w(:,:),     str_module, 4, in_subset=p_patch%cells%owned)

    END IF

  END SUBROUTINE couple_ocean


END MODULE mo_nwp_ocean_coupling
