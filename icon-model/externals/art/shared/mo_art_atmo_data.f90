!
! mo_art_data
! This module comprises all variables that are used as input from the host
! model.
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

MODULE mo_art_atmo_data
! ICON
  USE mo_kind,                          ONLY: wp
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: t_art_atmo


  ! ----------------------------------
  ! --- General ART data structure
  ! ----------------------------------


  TYPE t_art_atmo
    REAL(wp), POINTER ::       &
      &  temp(:,:,:) => NULL(),          &  !< air temperature (K)
      &  tempv(:,:,:) => NULL(),         &  !< virtual air temperature (K)
      &  temp_ifc(:,:,:) => NULL(),      &  !< temperature on interfaces (K)
      &  pres(:,:,:) => NULL(),          &  !< air pressure (Pa)
      &  exner(:,:,:) => NULL(),         &  !< Exner pressure (-)
      &  vor(:,:,:) => NULL(),           &  !< Vorticity (s-1)
      &  rho(:,:,:) => NULL(),           &  !< air density (kg/m3)
      &  theta_v(:,:,:) => NULL(),       &  !< virtual potential temperature (K)
      &  z_mc(:,:,:) => NULL(),          &  !< model level altitude (m)
      &  z_ifc(:,:,:) => NULL(),         &  !< model interface altitude (m)
      &  dz(:,:,:) => NULL(),            &  !< model layer heights (m)
      &  pres_ifc(:,:,:) => NULL(),      &  !< air pressure at interfaces (Pa)
      &  dpres_mc(:,:,:) => NULL(),      &  !< pressure difference between 
                                            !  two full model levels (Pa)
      &  o3_field_icon(:,:,:) => NULL(), &  !< O3 mixing ratio in ICON (kg/kg)
      &  o3_ext(:,:,:) => NULL()


    REAL(wp), POINTER ::  &
      &   fr_land(:,:) => NULL()  !< fraction that is covered by land

    LOGICAL, ALLOCATABLE ::  &
      &  llsm(:,:)                !< Land sea mask


    REAL(wp), ALLOCATABLE :: &
      &  o3_clim(:,:,:)           !< climatologic ozone calc. in ART (kg/kg)
     


    REAL(wp), POINTER ::               &
      &  cell_area(:,:) => NULL(),     & !< area of cells (m2)
      &  fr_glac(:,:) => NULL()          !< fraction that is covered by ice

    REAL(wp), ALLOCATABLE ::  &
      &  theta(:,:,:),        &  !< potential temperature (K)
      &  lat(:,:),            & !< latitude (rad)
      &  lon(:,:)               !< longitude (rad)

    ! CLoud and surface properties
    REAL(wp), POINTER ::             &
      &  clc(:,:,:) => NULL(),       &  !< cloud cover
      &  tot_cld(:,:,:,:) => NULL(), &  !< mass mixing ratio in cloud (kg/kg?)
      &  acdnc(:,:,:) => NULL(),     &
      &  albedo(:,:) => NULL()          !< surface albedo (-)

    ! for emissions and sedimentation
    REAL(wp), ALLOCATABLE ::   &
      &  swflx_par_sfc(:,:),   &  !< Photosynthetic active radiation at the surface (W m-2)
      &  swflxsfc(:,:),        &  !< shortwave net flux at surface (W m-2)
      &  t_2m(:,:),            &  !< air temperature in an altitude of 2 m (K)
      &  rh_2m(:,:),           &  !< relative humidity in an altitude of 2 m (%)
      &  gz0(:,:)                 !< gravity acceleration times surface altitude (m2/s2)

    REAL(wp), POINTER ::               &
      &  lai(:,:) => NULL(),           &  !< leaf area index (-)
      &  u(:,:,:) => NULL(),           &  !< eastward horizontal wind speed (m/s)
      &  v(:,:,:) => NULL(),           &  !< northward horizontal wind speed (m/s)
      &  tke(:,:,:) => NULL(),         &  !< turbulent kinetic energy (m2/s2)
      &  pres_sfc(:,:) => NULL(),      &  !< surface pressure (Pa)
      &  u_10m(:,:) => NULL(),         &  !< eastward horizontal wind speed in 10 m height (m/s)
      &  v_10m(:,:) => NULL(),         &  !< northward horizontal wind speed in 10 m height (m/s)
      &  rain_gsp_rate(:,:) => NULL(), &  !< rain rate due to grid scale rain (??)
      &  rain_con_rate(:,:) => NULL(), &  !< rain rate due to convective precipitation (??)
      &  tch(:,:) => NULL(),           &  !< ??
      &  tcm(:,:) => NULL()               !< ??

    REAL(wp),POINTER ::                & !< tropop. diagnostics
      &  ptropo(:,:),                  & !< level index below tropopause
      &  ktrpwmop1_real(:,:)
    INTEGER, POINTER ::                &
      &  ktrpwmop1(:,:),               & !< level index below tropopause
      &  ktrpwmo(:,:)                    !< level index above tropopause

    INTEGER ::                 &
      &  nroot,                & !< number of bisection (the "n" of R<n>B<k>)
      &  nbisect,              & !< number of bisection (the "k" of R<n>B<k>)
      &  nproma, nlev, nblks,  & !< dimensions (nproma, number of vertical levels, 
                                 !  and number of blocks)
      &  npromz,               & !< number of grid points in last block
      &  nlevp1,               & !< number of vertical levels plus 1
      &  i_startblk, i_endblk    !< loop indices for blocks

    REAL(wp), POINTER ::       &
      &  sza(:,:) => NULL(),   &   !< cosine of solar zenith angle
      &  sza_deg(:,:) => NULL()    !< solar zenith angle in degrees

    INTEGER ::                 &
      &  idust_insol_acc,      &   !< tracer index for insoluble dust in accumulation mode
      &  idust_insol_coa,      &   !< tracer index for insoluble dust in coarse       mode
      &  idust_giant               !< tracer index for           dust in giant        mode

  END TYPE t_art_atmo


END MODULE mo_art_atmo_data
