! Data type defintion for wave forcing data state
!
! Defines the data type for storing wave-specific forcing
! fields.
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

MODULE mo_wave_forcing_types

  USE mo_kind,               ONLY: wp


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: t_wave_forcing

  ! wave-specific forcing data type
  !
  TYPE :: t_wave_forcing

    REAL(wp), POINTER, CONTIGUOUS ::  &
      & u10m (:,:),       & ! zonal wind in 10m (nproma,nblks_c)      ( m/s )
      & v10m (:,:),       & ! meridional wind in 10m (nproma,nblks_c) ( m/s )
      & sp10m(:,:),       & ! wind speed in 10m  (nproma,nblks_c)           ( m/s )
      & dir10m(:,:),      & ! wind direction in 10m (nproma,nblks_c)        ( deg )
      & sea_ice_c(:,:),   & ! sea ice concentration at centers (fraction of 1)
      & sea_ice_e(:,:),   & ! sea ice concentration at edges (fraction of 1)
      & sea_level_c(:,:), & ! sea level height at centers (nproma,nblks_c) ( m )
      & sea_level_e(:,:), & ! sea level height at edges (nproma,nblks_e)   ( m )
      & usoce_c(:,:),     & ! zonal ocean surface current at centers    (nproma,nblks_c) ( m/s )
      & vsoce_c(:,:),     & ! meridional ocean surface current at centers (nproma,nblks_c) ( m/s )
      & sp_soce_c(:,:),   & ! ocean surface current velocity at centers  (nproma,nblks_c) ( m/s )
      & dir_soce_c(:,:),  & ! ocean surface current direction at centers (nproma,nblks_c) ( m/s )
      & usoce_e(:,:),     & ! zonal ocean surface current at edges (nproma,nblks_e)    ( m/s )
      & vsoce_e(:,:)      & ! meridional ocean surface current at edges (nproma,nblks_e) ( m/s )
      & => NULL()

    INTEGER, POINTER, CONTIGUOUS :: &
      & ice_free_mask_c(:,:) & ! ice-free mask (nproma,nblks_c) 1 - no ice, 0 - ice
      & => NULL()

  END TYPE t_wave_forcing

END MODULE mo_wave_forcing_types
