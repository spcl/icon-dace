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

! @brief Interface between AES physics and the ocean, through a coupler

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_aes_ocean_coupling

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nonhydro_types      ,ONLY: t_nh_diag
  USE mo_aes_phy_memory      ,ONLY: prm_field
  USE mo_ccycle_config       ,ONLY: ccycle_config

  USE mo_parallel_config     ,ONLY: nproma

  USE mo_run_config          ,ONLY: ico2
  USE mo_aes_sfc_indices     ,ONLY: iwtr, iice, ilnd, nsfc_type
  USE mo_aes_phy_config      ,ONLY: aes_phy_config
  USE mo_aes_vdf_config      ,ONLY: aes_vdf_config

  USE mo_sync                ,ONLY: sync_c, sync_patch_array

  USE mo_bc_greenhouse_gases ,ONLY: ghg_co2vmr

  USE mo_coupling_config     ,ONLY: is_coupled_to_ocean
  USE mo_atmo_coupling_frame ,ONLY: nbr_inner_cells
  USE mo_atmo_ocean_coupling ,ONLY: mask_checksum, &
    & field_id_co2_flx, field_id_co2_vmr, field_id_freshflx, &
    & field_id_heatflx, field_id_oce_u, field_id_oce_v, field_id_pres_msl, &
    & field_id_seaice_atm, field_id_seaice_oce, field_id_sp10m, field_id_sst, &
    & field_id_umfl, field_id_vmfl
  USE mo_exception           ,ONLY: finish

  USE mo_coupling_utils      ,ONLY: cpl_put_field, cpl_get_field

  USE mo_util_dbg_prnt       ,ONLY: dbg_print
  USE mo_dbg_nml             ,ONLY: idbg_mxmn, idbg_val
  USE mo_physical_constants  ,ONLY: amd, amco2
  USE mo_physical_constants  ,ONLY: cvd, cpd
  USE mo_fortran_tools       ,ONLY: init

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: interface_aes_ocean

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_aes_ocean_coupling'  ! Output of module for 1 line debug

CONTAINS

  !>
  !! SUBROUTINE interface_aes_ocean -- the interface between
  !! AES physics and the ocean, through a coupler
  !!
  !! This subroutine is called in the time loop of the ICONAM model.
  !! It takes the following as input:
  !! <ol>
  !! <li> prognostic and diagnostic variables of the dynamical core;
  !! <li> tendency of the prognostic varibles induced by adiabatic dynamics;
  !! <li> time step;
  !! <li> information about the dynamics grid;
  !! <li> interplation coefficients.
  !! </ol>
  !!
  !! The output includes tendencies of the prognostic variables caused by
  !! the parameterisations.
  !!
  !! Note that each call of this subroutine deals with a single grid level
  !! rather than the entire grid tree.

  SUBROUTINE interface_aes_ocean( p_patch , pt_diag)

    ! Arguments

    TYPE(t_patch), TARGET, INTENT(INOUT)    :: p_patch
    TYPE(t_nh_diag), TARGET, INTENT(INOUT)  :: pt_diag

    ! Local variables

    INTEGER               :: nbr_hor_cells  ! = inner and halo points
    INTEGER               :: jg             ! grid index
    INTEGER               :: nlev           ! number of levels
    INTEGER               :: nblks_c        ! number of blocks
    INTEGER               :: n              ! nproma loop count
    INTEGER               :: nn             ! block offset
    INTEGER               :: i_blk          ! block loop count
    INTEGER               :: nlen           ! nproma/npromz
    INTEGER               :: no_arr         !  no of arrays in bundle for put/get calls

    REAL(wp), PARAMETER   :: dummy = 0.0_wp

    REAL(wp)              :: shflx_adjustment_factor

    REAL(wp)              :: scr(nproma,p_patch%alloc_cell_blocks)
    REAL(wp)              :: frac_oce(nproma,p_patch%alloc_cell_blocks)

    REAL(wp), ALLOCATABLE :: put_buffer(:,:,:)
    REAL(wp), ALLOCATABLE :: get_buffer(:,:)
    LOGICAL :: received_data

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':interface_aes_ocean'

    IF ( .NOT. is_coupled_to_ocean() ) RETURN

    ! adjust size if larger bundles are used (no_arr > 2 below)

    ALLOCATE(put_buffer(nproma,p_patch%nblks_c,2))

    ! adjust size if larger bundles are used (no_arr > 4 below)

    nbr_hor_cells = p_patch%n_patch_cells
    ALLOCATE(get_buffer(nbr_hor_cells,4))

    jg   = p_patch%id
    nlev = p_patch%nlev
    nblks_c = p_patch%nblks_c

    !-------------------------------------------------------------------------
    ! If running in atm-oce coupled mode, exchange information
    !-------------------------------------------------------------------------
    !
    ! Possible fields that contain information to be sent to the ocean include
    !
    ! 1. prm_field(jg)% u_stress_tile(:,:,iwtr/iice)  and
    !    prm_field(jg)% v_stress_tile(:,:,iwtr/iice)  which are the wind stress components over water and ice respectively
    !
    ! 2. prm_field(jg)% evap_tile(:,:,iwtr/iice)  evaporation rate over ice-covered and open ocean/lakes, no land;
    !
    ! 3. prm_field(jg)%rsfl and prm_field(jg)%ssfl
    !    which gives the liquid and solid precipitation rates, respectively
    !
    ! 4. prm_field(jg)% ta(:,nlev,:)  temperature at the lowest model level, or
    !    prm_field(jg)% tas(:,:)      2-m temperature, not available yet, or
    !    prm_field(jg)% shflx_tile(:,:,iwtr) sensible heat flux
    !    ... tbc
    !
    ! 5  prm_field(jg)% lhflx_tile(:,:,iwtr) latent heat flux
    ! 6. shortwave radiation flux at the surface
    !
    ! Possible fields to receive from the ocean include
    !
    ! 1. prm_field(jg)% ts_tile(:,:,iwtr)   SST
    ! 2. prm_field(jg)% ocu(:,:) and ocv(:,:) ocean surface current
    ! 3. ... tbc

    !  Send fields to ocean:
    !   "surface_downward_eastward_stress" bundle  - zonal wind stress component over ice and water
    !   "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water
    !   "surface_fresh_water_flux" bundle          - liquid rain, snowfall, evaporation
    !   "total heat flux" bundle                   - short wave, long wave, sensible, latent heat flux
    !   "atmosphere_sea_ice_bundle"                - sea ice surface and bottom melt potentials
    !   "10m_wind_speed"                           - atmospheric wind speed
    !   "qtrc_phy(nlev,co2)"                       - co2 mixing ratio
    !   "pres_msl"                                 - sea level pressure
    !
    !  Receive fields from ocean:
    !   "sea_surface_temperature"                  - SST
    !   "eastward_sea_water_velocity"              - zonal velocity, u component of ocean surface current
    !   "northward_sea_water_velocity"             - meridional velocity, v component of ocean surface current
    !   "ocean_sea_ice_bundle"                     - ice thickness, snow thickness, ice concentration
    !   "co2_flux"                                 - ocean co2 flux
    !
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Send fields from atmosphere to ocean
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****

    ! Calculate fractional ocean mask
    ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
    !  - lake part is included in land part, must be subtracted as well
    !  - if no lake part is present, subtract land part only
    !  - if no jsbach is present (aquaplanet), frac_oce is 1.

    !$ACC DATA CREATE(frac_oce, get_buffer, put_buffer)

    ! As YAC does not touch masked data an explicit initialisation
    ! is required as some compilers are asked to initialise with NaN
    ! and as we loop over the full array.

!ICON_OMP_PARALLEL
    CALL init(put_buffer(:,:,:))
!ICON_OMP_END_PARALLEL

    IF ( mask_checksum > 0 .AND. aes_phy_config(jg)%ljsb ) THEN
      IF ( aes_phy_config(jg)%llake ) THEN
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
        DO i_blk = 1, p_patch%nblks_c
          IF (i_blk /= p_patch%nblks_c) THEN
            nlen = nproma
          ELSE
            nlen = p_patch%npromz_c
          END IF
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
          DO n = 1, nlen
            frac_oce(n,i_blk) = 1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd) - prm_field(jg)%alake(n,i_blk)
          ENDDO
        ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ELSE
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
        DO i_blk = 1, p_patch%nblks_c
          IF (i_blk /= p_patch%nblks_c) THEN
            nlen = nproma
          ELSE
            nlen = p_patch%npromz_c
          END IF
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
          DO n = 1, nlen
            frac_oce(n,i_blk) = 1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd)
          ENDDO
        ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
      ENDIF
    ELSE
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen
          frac_oce(n,i_blk) = 1.0
        ENDDO
      ENDDO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL
    ENDIF

    ! ------------------------------
    !  Send zonal wind stress bundle
    !   field_id(1) represents "surface_downward_eastward_stress" bundle - zonal wind stress component over ice and water

    !$ACC UPDATE HOST(prm_field(jg)%u_stress_tile(:,:,iwtr)) ASYNC(1)
    !$ACC UPDATE HOST(prm_field(jg)%u_stress_tile(:,:,iice)) ASYNC(1)
    !$ACC WAIT(1)
    CALL cpl_put_field( &
      routine, field_id_umfl, 'u-stress', nbr_hor_cells, &
      field_1=prm_field(jg)%u_stress_tile(:,1:nblks_c,iwtr), &
      field_2=prm_field(jg)%u_stress_tile(:,1:nblks_c,iice))

    ! ------------------------------
    !  Send meridional wind stress bundle
    !   "surface_downward_northward_stress" bundle - meridional wind stress component over ice and water

    !$ACC UPDATE HOST(prm_field(jg)%v_stress_tile(:,:,iwtr)) ASYNC(1)
    !$ACC UPDATE HOST(prm_field(jg)%v_stress_tile(:,:,iice)) ASYNC(1)
    !$ACC WAIT(1)
    CALL cpl_put_field( &
      routine, field_id_vmfl, 'v-stress', nbr_hor_cells, &
      field_1=prm_field(jg)%v_stress_tile(:,1:nblks_c,iwtr), &
      field_2=prm_field(jg)%v_stress_tile(:,1:nblks_c,iice))

    ! ------------------------------
    !  Send surface fresh water flux bundle
    !   "surface_fresh_water_flux" bundle - liquid rain, snowfall, evaporation
    !
    !   Note: the evap_tile should be properly updated and added;
    !         as long as evaporation over sea-ice is not used in ocean thermodynamics, the evaporation over the
    !         whole ocean part of grid-cell is passed to the ocean

    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) scr(:,:) = 0.0_wp

    ! total rates of rain and snow over whole cell
    !$ACC UPDATE HOST(prm_field(jg)%rsfl, prm_field(jg)%ssfl) ASYNC(1)
    !$ACC WAIT(1)

    ! Aquaplanet coupling: surface types ocean and ice only
    IF (nsfc_type == 2) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen

          ! evaporation over ice-free and ice-covered water fraction - of whole ocean part
          put_buffer(n,i_blk,1) = &
            prm_field(jg)%evap_tile(n,i_blk,iwtr) * &
            prm_field(jg)%frac_tile(n,i_blk,iwtr) + &
            prm_field(jg)%evap_tile(n,i_blk,iice) * &
            prm_field(jg)%frac_tile(n,i_blk,iice)
        ENDDO
      ENDDO
      !$ACC UPDATE HOST(put_buffer(:,:,1)) ASYNC(1)
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    ! We zero out explicitly the put_buffer after the end of the last (possibly shorter) block
    ! in case it contains NaNs, Infs or other problematic stuff after transfer from the GPU
    ! This may not be necessary when running on CPUs, but the cost is negligble so we do it
    ! anyway. 
    put_buffer(p_patch%npromz_c+1:nproma, p_patch%nblks_c,1) = 0.0_wp
    ! Full coupling including jsbach: surface types ocean, ice, land
    ELSE IF (nsfc_type == 3) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      !$ACC DATA COPYOUT(scr) IF(idbg_mxmn >= 1 .OR. idbg_val >=1)
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) NO_CREATE(scr)
        DO n = 1, nlen

          ! evaporation over ice-free and ice-covered water fraction, of whole ocean part, without land part
          !  - lake part is included in land part, must be subtracted as well
          !    frac_oce(n,i_blk)= 1.0_wp-prm_field(jg)%frac_tile(n,i_blk,ilnd)-prm_field(jg)%alake(n,i_blk)

          IF (frac_oce(n,i_blk) <= 0.0_wp) THEN
            ! land part is zero
            put_buffer(n,i_blk,1) = 0.0_wp
          ELSE
            put_buffer(n,i_blk,1) = &
              (prm_field(jg)%evap_tile(n,i_blk,iwtr) * &
               prm_field(jg)%frac_tile(n,i_blk,iwtr) + &
               prm_field(jg)%evap_tile(n,i_blk,iice) * &
               prm_field(jg)%frac_tile(n,i_blk,iice))/frac_oce(n,i_blk)
          ENDIF
          IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) &
            scr(n,i_blk) = put_buffer(n,i_blk,1)
        ENDDO
      ENDDO 
      !$ACC UPDATE HOST(put_buffer(:,:,1)) ASYNC(1)
      !$ACC WAIT(1)
      !$ACC END DATA
!ICON_OMP_END_PARALLEL_DO
      put_buffer(p_patch%npromz_c+1:nproma, p_patch%nblks_c,1) = 0.0_wp

      IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
        &  CALL dbg_print('AESOce: evapo-cpl',scr,str_module,3,in_subset=p_patch%cells%owned)
    ELSE
      CALL finish( &
        routine, 'coupling only for nsfc_type equals 2 or 3. ' // &
        'Check your code/configuration!')
    ENDIF  !  nsfc_type

    CALL cpl_put_field( &
      routine, field_id_freshflx, 'fresh water flux', nbr_hor_cells, &
      field_1=prm_field(jg)%rsfl, &
      field_2=prm_field(jg)%ssfl, &
      field_3=put_buffer(:,:,1))

    ! ------------------------------
    !  Send total heat flux bundle
    !   "total heat flux" bundle - short wave, long wave, sensible, latent heat flux

    !$ACC UPDATE HOST(prm_field(jg)%swflxsfc_tile(:,:,iwtr)) ASYNC(1)
    !$ACC UPDATE HOST(prm_field(jg)%lwflxsfc_tile(:,:,iwtr)) ASYNC(1)
    !$ACC UPDATE HOST(prm_field(jg)%lhflx_tile(:,:,iwtr)) ASYNC(1)
    !$ACC WAIT(1)
    IF (aes_phy_config(jg)%use_shflx_adjustment .AND. &
        .NOT. aes_vdf_config(jg)%use_tmx) THEN

      shflx_adjustment_factor = cvd/cpd

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen
          put_buffer(n,i_blk,1) = &
            shflx_adjustment_factor*prm_field(jg)%shflx_tile(n,i_blk,iwtr)
        ENDDO
      ENDDO 
      !$ACC UPDATE HOST(put_buffer(:,:,1)) ASYNC(1)
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
      put_buffer(p_patch%npromz_c+1:nproma, p_patch%nblks_c,1) = 0.0_wp

      CALL cpl_put_field( &
        routine, field_id_heatflx, 'heat flux', nbr_hor_cells, &
        field_1=prm_field(jg)%swflxsfc_tile(:,:,iwtr), &
        field_2=prm_field(jg)%lwflxsfc_tile(:,:,iwtr), &
        field_3=put_buffer(:,:,1), &
        field_4=prm_field(jg)%lhflx_tile(:,:,iwtr))

    ELSE ! .NOT. use_shflx_adjustment .OR. use_tmx

      !$ACC UPDATE HOST(prm_field(jg)%shflx_tile(:,:,iwtr)) ASYNC(1)
      !$ACC WAIT(1)

      CALL cpl_put_field( &
        routine, field_id_heatflx, 'heat flux', nbr_hor_cells, &
        field_1=prm_field(jg)%swflxsfc_tile(:,:,iwtr), &
        field_2=prm_field(jg)%lwflxsfc_tile(:,:,iwtr), &
        field_3=prm_field(jg)%shflx_tile(:,:,iwtr), &
        field_4=prm_field(jg)%lhflx_tile(:,:,iwtr))

    ENDIF ! use_shflx_adjustment

    ! ------------------------------
    !  Send sea ice flux bundle
    !   "atmosphere_sea_ice_bundle" - sea ice surface and bottom melt potentials Qtop, Qbot

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
    DO i_blk = 1, p_patch%nblks_c
      IF (i_blk /= p_patch%nblks_c) THEN
        nlen = nproma
      ELSE
        nlen = p_patch%npromz_c
      END IF
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
      DO n = 1, nlen
        put_buffer(n,i_blk,1) = prm_field(jg)%Qtop(n,1,i_blk)
        put_buffer(n,i_blk,2) = prm_field(jg)%Qbot(n,1,i_blk)
      ENDDO
    ENDDO

    !$ACC UPDATE HOST(put_buffer(:,:,1:2)) ASYNC(1)
    !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
    put_buffer(p_patch%npromz_c+1:nproma, p_patch%nblks_c,1:2) = 0.0_wp

    CALL cpl_put_field( &
      routine, field_id_seaice_atm, 'atmos sea ice', nbr_hor_cells, &
      field_1=put_buffer(:,:,1), &
      field_2=put_buffer(:,:,2))

    ! ------------------------------
    !  Send 10m wind speed
    !   "10m_wind_speed" - atmospheric wind speed

    !$ACC UPDATE HOST(prm_field(jg)%sfcWind) ASYNC(1)
    !$ACC WAIT(1)

    CALL cpl_put_field( &
      routine, field_id_sp10m, 'wind speed', nbr_hor_cells, &
      prm_field(jg)%sfcWind(:,1:nblks_c))

    ! ------------------------------
    !  Send sea level pressure
    !   "pres_msl" - atmospheric sea level pressure

    !$ACC UPDATE HOST(pt_diag%pres_msl) ASYNC(1)
    !$ACC WAIT(1)

    CALL cpl_put_field( &
      routine, field_id_pres_msl, 'sea level pressure', nbr_hor_cells, &
      pt_diag%pres_msl(:,1:nblks_c))

#ifndef __NO_ICON_OCEAN__
    IF (ccycle_config(jg)%iccycle /= 0) THEN

       ! ------------------------------
       !  Send co2 mixing ratio
       !  "co2_mixing_ratio" - CO2 mixing ratio in ppmv

!!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
       DO i_blk = 1, p_patch%nblks_c
          IF (i_blk /= p_patch%nblks_c) THEN
             nlen = nproma
          ELSE
             nlen = p_patch%npromz_c
          END IF
          SELECT CASE (ccycle_config(jg)%iccycle)
          CASE (1) ! c-cycle with interactive atm. co2 concentration, qtrc_phy in kg/kg
             !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
             DO n = 1, nlen
                put_buffer(n,i_blk,1) = amd/amco2 * 1.0e6_wp * prm_field(jg)%qtrc_phy(n,nlev,i_blk,ico2)
             END DO
          CASE (2) ! c-cycle with prescribed  atm. co2 concentration
             SELECT CASE (ccycle_config(jg)%ico2conc)
             CASE (2) ! constant  co2 concentration, vmr_co2 in m3/m3
                !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
                DO n = 1, nlen
                   put_buffer(n,i_blk,1) = 1.0e6_wp * ccycle_config(jg)%vmr_co2
                END DO
             CASE (4) ! transient co2 concentration, ghg_co2vmr in m3/m3
                !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
                DO n = 1, nlen
                   put_buffer(n,i_blk,1) = 1.0e6_wp * ghg_co2vmr
                END DO
             END SELECT
          END SELECT
       ENDDO
       !$ACC UPDATE HOST(put_buffer(:,:,1)) ASYNC(1)
       !$ACC WAIT(1)
!!ICON_OMP_END_PARALLEL_DO
       put_buffer(p_patch%npromz_c+1:nproma, p_patch%nblks_c,1) = 0.0_wp

       CALL cpl_put_field( &
        routine, field_id_co2_vmr, 'co2 mr', nbr_hor_cells, &
        put_buffer(:,:,1))

    ENDIF
#endif

    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !  Receive fields from ocean to atmosphere
    !  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****  *****
    !
    !  Receive fields, only assign values if something was received
    !   - ocean fields have undefined values on land, which are not sent to the atmosphere,
    !     therefore get_buffer is set to zero to avoid unintended usage of ocean values over land

    ! ------------------------------
    !  Receive SST
    !   "sea_surface_temperature" - SST

!ICON_OMP_PARALLEL
    CALL init(get_buffer(:,:))
!ICON_OMP_END_PARALLEL

    no_arr = 1
    CALL cpl_get_field( &
      routine, field_id_sst, 'SST', &
      get_buffer(1:nbr_hor_cells,1:no_arr), &
      first_get=.TRUE., received_data=received_data)

    IF (received_data) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      !$ACC UPDATE DEVICE(get_buffer(:,1)) ASYNC(1)
      !$ACC DATA COPYOUT(scr) IF(idbg_mxmn >= 1 .OR. idbg_val >=1)
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) NO_CREATE(scr)
        DO n = 1, nlen
          !  - lake part is included in land part, must be subtracted as well, see frac_oce

          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ts_tile(n,i_blk,iwtr) = dummy
          ELSE
            IF ( frac_oce(n,i_blk) > EPSILON(1.0_wp) ) &
              prm_field(jg)%ts_tile(n,i_blk,iwtr) = get_buffer(nn+n,1)
            IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN
              IF ( frac_oce(n,i_blk) > 0.0_wp ) THEN
                scr(n,i_blk) = get_buffer(nn+n,1)
              ELSE
                scr(n,i_blk) = 285.0_wp  !  value over land - for dbg_print
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      !$ACC WAIT(1)
      !$ACC END DATA
!ICON_OMP_END_PARALLEL_DO
      IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 )  &
        &  CALL dbg_print('AESOce: SSToce-cpl',scr,str_module,4,in_subset=p_patch%cells%owned)

      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ts_tile(:,:,iwtr))
    END IF
    !
    ! ------------------------------
    !  Receive zonal velocity
    !   "eastward_sea_water_velocity" - zonal velocity, u component of ocean surface current
    !
    no_arr = 1
    CALL cpl_get_field( &
      routine, field_id_oce_u, 'u velocity', &
      get_buffer(1:nbr_hor_cells,1:no_arr), received_data=received_data)

    IF (received_data) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      !$ACC UPDATE DEVICE(get_buffer(:,1)) ASYNC(1)
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ocu(n,i_blk) = dummy
          ELSE
            prm_field(jg)%ocu(n,i_blk) = get_buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocu(:,:))
    END IF

    ! ------------------------------
    !  Receive meridional velocity
    !   "northward_sea_water_velocity" - meridional velocity, v component of ocean surface current
    !
    no_arr = 1
    CALL cpl_get_field( &
      routine, field_id_oce_v, 'v velocity', &
      get_buffer(1:nbr_hor_cells,1:no_arr), received_data=received_data)

    IF (received_data) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      !$ACC UPDATE DEVICE(get_buffer(:,1)) ASYNC(1)
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%ocv(n,i_blk) = dummy
          ELSE
            prm_field(jg)%ocv(n,i_blk) = get_buffer(nn+n,1)
          ENDIF
        ENDDO
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%ocv(:,:))
    END IF

    ! ------------------------------
    !  Receive sea ice bundle
    !   "ocean_sea_ice_bundle" - ice thickness, snow thickness, ice concentration
    !
    no_arr = 3
    CALL cpl_get_field( &
      routine, field_id_seaice_oce, 'sea ice', &
      get_buffer(1:nbr_hor_cells,1:no_arr), received_data=received_data)

    IF (received_data) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
      !$ACC UPDATE DEVICE(get_buffer(:,1:3)) ASYNC(1)
      DO i_blk = 1, p_patch%nblks_c
        nn = (i_blk-1)*nproma
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen
          IF ( nn+n > nbr_inner_cells ) THEN
            prm_field(jg)%hi  (n,1,i_blk) = dummy
            prm_field(jg)%hs  (n,1,i_blk) = dummy
            prm_field(jg)%conc(n,1,i_blk) = dummy
          ELSE
            prm_field(jg)%hi  (n,1,i_blk) = get_buffer(nn+n,1)
            prm_field(jg)%hs  (n,1,i_blk) = get_buffer(nn+n,2)
            prm_field(jg)%conc(n,1,i_blk) = get_buffer(nn+n,3)
          ENDIF
        ENDDO
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hi  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%hs  (:,1,:))
      CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%conc(:,1,:))

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nlen) ICON_OMP_RUNTIME_SCHEDULE
      DO i_blk = 1, p_patch%nblks_c
        IF (i_blk /= p_patch%nblks_c) THEN
          nlen = nproma
        ELSE
          nlen = p_patch%npromz_c
        END IF
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
        DO n = 1, nlen
          prm_field(jg)%seaice(n,i_blk) = prm_field(jg)%conc(n,1,i_blk)
          prm_field(jg)%siced(n,i_blk)  = prm_field(jg)%hi(n,1,i_blk)
        ENDDO
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO

    END IF

    IF (ccycle_config(jg)%iccycle /= 0) THEN
      !
      ! ------------------------------
      !  Receive co2 flux
      !   "co2_flux" - ocean co2 flux
      !
      get_buffer(:,1) = 0.0_wp ! needs to be checked if this is necessary
      no_arr = 1
      CALL cpl_get_field( &
        routine, field_id_co2_flx, 'CO2 flux', &
        get_buffer(1:nbr_hor_cells,1:no_arr), received_data=received_data)

       IF (received_data) THEN

!ICON_OMP_PARALLEL_DO PRIVATE(i_blk, n, nn, nlen) ICON_OMP_RUNTIME_SCHEDULE
          !$ACC UPDATE DEVICE(get_buffer(:,1)) ASYNC(1)
          DO i_blk = 1, p_patch%nblks_c
             nn = (i_blk-1)*nproma
             IF (i_blk /= p_patch%nblks_c) THEN
                nlen = nproma
             ELSE
                nlen = p_patch%npromz_c
             END IF
             !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1)
             DO n = 1, nlen
                IF ( nn+n > nbr_inner_cells ) THEN
                   prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = dummy
                ELSE
                   prm_field(jg)%co2_flux_tile(n,i_blk,iwtr) = get_buffer(nn+n,1)
                ENDIF
             ENDDO
          ENDDO
          !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
          !
          CALL sync_patch_array(sync_c, p_patch, prm_field(jg)%co2_flux_tile(:,:,iwtr))
        ENDIF

    END IF

!---------DEBUG DIAGNOSTICS-------------------------------------------

    ! calculations for debug print output for namelist debug-values >0 only
    IF ( idbg_mxmn >= 1 .OR. idbg_val >=1 ) THEN

      ! u/v-stress on ice and water sent
      scr(:,:) = prm_field(jg)%u_stress_tile(:,:,iwtr)
      CALL dbg_print('AESOce: u_stress.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%u_stress_tile(:,:,iice)
      CALL dbg_print('AESOce: u_stress.ice',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%v_stress_tile(:,:,iwtr)
      CALL dbg_print('AESOce: v_stress.wtr',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%v_stress_tile(:,:,iice)
      CALL dbg_print('AESOce: v_stress.ice',scr,str_module,4,in_subset=p_patch%cells%owned)

      ! rain, snow, evaporation
      scr(:,:) = prm_field(jg)%rsfl(:,:)
      CALL dbg_print('AESOce: total rain  ',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%ssfl(:,:)
      CALL dbg_print('AESOce: total sn/grp',scr,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: evaporation ',prm_field(jg)%evap   ,str_module,4,in_subset=p_patch%cells%owned)

      ! total: short wave, long wave, sensible, latent heat flux sent
      scr(:,:) = prm_field(jg)%swflxsfc_tile(:,:,iwtr) + &
        &        prm_field(jg)%lwflxsfc_tile(:,:,iwtr) + &
        &        shflx_adjustment_factor*prm_field(jg)%shflx_tile(:,:,iwtr)    + &
        &        prm_field(jg)%lhflx_tile(:,:,iwtr)
      CALL dbg_print('AESOce: totalhfx.wtr',scr,str_module,2,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%swflxsfc_tile(:,:,iwtr)
      CALL dbg_print('AESOce: swflxsfc.wtr',scr,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%lwflxsfc_tile(:,:,iwtr)
      CALL dbg_print('AESOce: lwflxsfc.wtr',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = shflx_adjustment_factor*prm_field(jg)%shflx_tile(:,:,iwtr)
      CALL dbg_print('AESOce: shflx.wtr   ',scr,str_module,4,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%lhflx_tile(:,:,iwtr)
      CALL dbg_print('AESOce: lhflx.wtr   ',scr,str_module,4,in_subset=p_patch%cells%owned)

      ! Qtop and Qbot, windspeed sent
      CALL dbg_print('AESOce: ice-Qtop    ',prm_field(jg)%Qtop        ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: ice-Qbot    ',prm_field(jg)%Qbot        ,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: sfcWind     ',prm_field(jg)%sfcWind     ,str_module,3,in_subset=p_patch%cells%owned)

      ! SST, sea ice, ocean velocity received
      scr(:,:) = prm_field(jg)%ts_tile(:,:,iwtr)
      CALL dbg_print('AESOce: ts_tile.iwtr',scr                       ,str_module,2,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: hi(1)       ',prm_field(jg)%hi(:,1,:)   ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: hs(1)       ',prm_field(jg)%hs(:,1,:)   ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: conc(1)     ',prm_field(jg)%conc(:,1,:) ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: siced       ',prm_field(jg)%siced       ,str_module,3,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: seaice      ',prm_field(jg)%seaice      ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: ocu         ',prm_field(jg)%ocu         ,str_module,4,in_subset=p_patch%cells%owned)
      CALL dbg_print('AESOce: ocv         ',prm_field(jg)%ocv         ,str_module,4,in_subset=p_patch%cells%owned)

      ! Fraction of tiles:
      !$ACC UPDATE HOST(frac_oce) ASYNC(1)
      !$ACC WAIT(1)
      CALL dbg_print('AESOce: frac_oce     ',frac_oce                 ,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%frac_tile(:,:,iwtr)
      CALL dbg_print('AESOce: frac_tile.wtr',scr                      ,str_module,3,in_subset=p_patch%cells%owned)
      scr(:,:) = prm_field(jg)%frac_tile(:,:,iice)
      CALL dbg_print('AESOce: frac_tile.ice',scr                      ,str_module,3,in_subset=p_patch%cells%owned)
      IF ( aes_phy_config(jg)%ljsb ) THEN
      scr(:,:) = prm_field(jg)%frac_tile(:,:,ilnd)
      CALL dbg_print('AESOce: frac_tile.lnd',scr                      ,str_module,4,in_subset=p_patch%cells%owned)
        IF ( aes_phy_config(jg)%llake ) &
          & CALL dbg_print('AESOce: frac_alake   ',prm_field(jg)%alake,str_module,4,in_subset=p_patch%cells%owned)
      ENDIF
    ENDIF
    !$ACC WAIT
    !$ACC END DATA ! frac_oce

    !---------------------------------------------------------------------

    DEALLOCATE(put_buffer, get_buffer)

  END SUBROUTINE interface_aes_ocean

END MODULE mo_aes_ocean_coupling
