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

! Conceptual copy of:
! * src/atm_dyn_iconam/mo_nonhydro_types
! * src/atm_phy_nwp/mo_nwp_phy_types
! * src/atm_phy_echam/mo_echam_phy_memory
! for the upper-atmosphere variables:
! * External data
! * Concentration of radiatively active gases
! * Gas properties
! * Tendencies from the upper-atmosphere physics parameterizations

MODULE mo_upatmo_types

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
  USE mo_fortran_tools,        ONLY: t_ptr_2d3d, t_ptr_3d
  USE mo_upatmo_utils,         ONLY: t_varstate

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_upatmo_diag
  PUBLIC :: t_upatmo_tend
  PUBLIC :: t_upatmo_extdat
  PUBLIC :: t_extdat_latlevtime
  PUBLIC :: t_upatmo

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_upatmo_types'

  ! Some notes:
  ! * The following fields are only required for NWP forcing
  ! * They will be allocated only if the upper-atmosphere
  !   is switched on

  !------------------------------------------------------
  !                Diagnostic variables
  !------------------------------------------------------  

  TYPE t_upatmo_diag

    REAL(wp), POINTER, CONTIGUOUS :: &

     gas(:,:,:,:),     &  ! Gas mass mixing ratio
                          ! (nproma,nlev,nblks_c,ngas) [kg/kg]
     mdry(:,:,:),      &  ! Dry air mass
                          ! (nproma,nlev,nblks_c) [kg/m2]  
     amd(:,:,:),       &  ! Molar mass of dry air 
                          ! (nproma,nlev,nblks_c) [g/mol]  
     cpair(:,:,:),     &  ! Heat capacity of (moist) air at constant pressure 
                          ! (nproma,nlev,nblks_c) [J/K/kg] 
     grav(:,:,:),      &  ! Gravitational acceleration 
                          ! (nproma,nlev,nblks_c) [m/s2] 
     sclrlw(:,:,:),    &  ! Scaling factor for heating rate 
                          ! from "standard" long-wave radiation
                          ! (nproma,nlev,nblks_c) [1] 
     effrsw(:,:,:)     &  ! Efficiency factor for heating rate 
                          ! from "standard" short-wave radiation
                          ! (nproma,nlev,nblks_c) [1] 
     => NULL()

    TYPE(t_ptr_2d3d), ALLOCATABLE :: gas_ptr(:)   ! Pointer for gas container

    ! Status
    LOGICAL :: linitialized = .FALSE. 

  END TYPE t_upatmo_diag  

  !------------------------------------------------------
  !      Tendencies from physics parameterizations
  !------------------------------------------------------  

  ! The large number of fields in the following might be regarded 
  ! as a significant waste of memory. 
  ! However, for at least the following reasons, 
  ! we would prefer to leave it this way for the time being:
  ! * The implementation of the upper-atmosphere physics into ICON 
  !   is still in its evaluation phase, so that the possibility 
  !   to output the tendencies of each single process is probably desirable
  ! * As mentioned above, allocation happens only if upper-atmosphere physics
  !   are switched on, so it should not do much harm to standard simulations
  ! Once the 1st reason should become less important, 
  ! the current implementation should allow to store tendencies 
  ! only group-wise instead of process-wise, i.e.: 
  ! * IMF
  ! * RAD
  ! instead of:
  ! * SRBC
  ! * NLTE
  ! * EUV
  ! * NO
  ! * CHEMHEAT
  ! * VDFMOL
  ! * FRIC
  ! * JOULE
  ! * IONDRAG

  !---------------------------

  TYPE t_ddt_tot_3d
    REAL(wp), POINTER, CONTIGUOUS :: tot(:,:,:) => NULL()
  END TYPE t_ddt_tot_3d

  !---------------------------

  TYPE t_ddt_tot_4d
    REAL(wp), POINTER, CONTIGUOUS :: tot(:,:,:,:) => NULL()

    TYPE(t_ptr_2d3d), ALLOCATABLE :: tot_ptr(:)
  END TYPE t_ddt_tot_4d

  !---------------------------

  TYPE t_ddt_tot_info
    INTEGER  :: nstate            = -999
    INTEGER  :: istartlev         = 1
    INTEGER  :: iendlev           = 0
    LOGICAL  :: linActivePhase    = .FALSE.
    LOGICAL  :: lafterActivePhase = .FALSE.
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: name     = " "
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: longname = " "
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit     = " "
  END TYPE t_ddt_tot_info

  !---------------------------

  TYPE t_upatmo_tend_tot
    TYPE(t_ddt_tot_3d),   ALLOCATABLE :: temp(:)  ! Accumulative temperature tendency (nstate)
    TYPE(t_ddt_tot_3d),   ALLOCATABLE :: exner(:) ! Accumulative Exner pressure tendency (nstate)
    TYPE(t_ddt_tot_3d),   ALLOCATABLE :: vn(:)    ! Accumulative edge-normal component of wind tendency (nstate)
    TYPE(t_ddt_tot_4d),   ALLOCATABLE :: qx(:)    ! Accumulative tracer tendencies (nstate)
    !
    TYPE(t_ddt_tot_info), ALLOCATABLE :: info(:)  ! Some metadata (ntendency)
    TYPE(t_varstate),     ALLOCATABLE :: state(:) ! State (ntendency)
  END TYPE t_upatmo_tend_tot

  !---------------------------

  TYPE t_upatmo_tend

    REAL(wp), POINTER , CONTIGUOUS ::                           &  ! Tendencies of ...

     ddt_temp_srbc(:,:,:),        &  ! temperature due to SRBC heating by O2
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_nlte(:,:,:),        &  ! temperature due to Non-LTE heating
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_euv(:,:,:),         &  ! temperature due to EUV heating
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_vdfmol(:,:,:),      &  ! temperature due to molecular diffusion
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_fric(:,:,:),        &  ! temperature due to frictional heating
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_no(:,:,:),          &  ! temperature due to NO NIR heating
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_chemheat(:,:,:),    &  ! temperature due to chemical heating
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_temp_joule(:,:,:),       &  ! temperature due to Joule heating
                                     ! (nproma,nlev,nblks_c) [K/s]
     ddt_u_vdfmol(:,:,:),         &  ! zonal wind component due to molecular diffusion
                                     ! (nproma,nlev,nblks_c) [m/s2]
     ddt_u_iondrag(:,:,:),        &  ! zonal wind component due to ion drag
                                     ! (nproma,nlev,nblks_c) [m/s2]
     ddt_v_vdfmol(:,:,:),         &  ! meridional wind component due to molecular diffusion
                                     ! (nproma,nlev,nblks_c) [m/s2]
     ddt_v_iondrag(:,:,:),        &  ! meridional wind component due to ion drag
                                     ! (nproma,nlev,nblks_c) [m/s2]
     ddt_qx_vdfmol(:,:,:,:)       &  ! tracer due to molecular diffusion 
                                     ! (currently, only specific humidity, [:,:,:,iqv]) 
                                     ! (nproma,nlev,nblks_c,1)
     => NULL()

    TYPE(t_ptr_2d3d), ALLOCATABLE :: ddt_qx_vdfmol_ptr(:) ! Pointer for tracer tendencies

    TYPE(t_upatmo_tend_tot)       :: ddt                  ! Accumulative tendencies

    ! Status
    LOGICAL :: linitialized = .FALSE.

  END TYPE t_upatmo_tend

  !------------------------------------------------------
  !                   External data
  !------------------------------------------------------ 

  TYPE t_extdat_intrpl_1d
    INTEGER,  ALLOCATABLE :: idx(:,:) ! Indices for interpolation (2,nlev)
    REAL(wp), ALLOCATABLE :: wgt(:,:) ! Interpolation weights (2,nlev)
  END TYPE t_extdat_intrpl_1d

  !---------------------------

  TYPE t_extdat_intrpl_2d
    INTEGER,  ALLOCATABLE :: idx(:,:,:) ! Indices for interpolation (2,nproma,nblks)
    REAL(wp), ALLOCATABLE :: wgt(:,:,:) ! Interpolation weights (2,nproma,nblks)
  END TYPE t_extdat_intrpl_2d

  !---------------------------

  TYPE t_extdat_intrpl
    TYPE(t_extdat_intrpl_2d) :: lat
    TYPE(t_extdat_intrpl_1d) :: lev
  END TYPE t_extdat_intrpl

  !---------------------------

  TYPE t_extdat_latlevtime
    REAL(wp), ALLOCATABLE :: data(:,:,:) ! Data (nlat,nlev,ntime)
    REAL(wp), ALLOCATABLE :: lat(:)      ! Latitudes (nlat)
    REAL(wp), ALLOCATABLE :: lev(:)      ! Full levels (nlev)
    REAL(wp), ALLOCATABLE :: lev_half(:) ! Half levels (nlev+1)
    REAL(wp), ALLOCATABLE :: time(:)     ! Times (ntime)
    !
    INTEGER               :: nlat        ! Number of latitudes
    INTEGER               :: nlev        ! Number of levels
    INTEGER               :: ntime       ! Number of times
    !
    INTEGER               :: istartlat   ! Start index of latitudes
    INTEGER               :: iendlat     ! End index of latitudes
    INTEGER               :: isteplat    ! Loop step for latitudes
    INTEGER               :: istartlev   ! Start index of full levels
    INTEGER               :: iendlev     ! End index of full levels
    INTEGER               :: isteplev    ! Loop step for full levels
    INTEGER               :: istarttime  ! Start index of times
    INTEGER               :: iendtime    ! End index of times
    INTEGER               :: isteptime   ! Loop step for times
    !
    INTEGER               :: data_id     ! Identifier for data
    INTEGER               :: lat_id      ! Identifier for type of latitudes
    INTEGER               :: lev_id      ! Identifier for type of levels
    INTEGER               :: time_id     ! Identifier for type of times
    !
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit_data = " " ! Unit of external data
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit_lat  = " " ! Unit of latitudes
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit_lev  = " " ! Unit of levels
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit_mass = " " ! Unit of mass
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: unit_time = " " ! Unit of times
    !
    TYPE(t_extdat_intrpl)          :: intrpl          ! For interpolation
  END TYPE t_extdat_latlevtime

  !---------------------------

  TYPE t_upatmo_extdat

    TYPE(t_extdat_latlevtime), ALLOCATABLE :: gas(:)           ! External data for radiatively active gases (ngas)
    TYPE(t_extdat_latlevtime)                 chemheat         ! External data for chemical heating tendencies

    INTEGER                                :: ngas             ! Number of gases
    INTEGER,                   ALLOCATABLE :: mapgasid2indx(:) ! Map global gas id to local gas index
    INTEGER,                   ALLOCATABLE :: mapgasindx2id(:) ! Map local gas index to global gas id

    TYPE(t_ptr_3d),            ALLOCATABLE :: gas_interm(:)    ! Gas on horizontal grid of ICON, 
                                                               ! but still on pressure levels of external data (ngas)
    ! Status
    LOGICAL :: linitialized = .FALSE.

  END TYPE t_upatmo_extdat

  !------------------------------------------------------
  !                  Collective type
  !------------------------------------------------------ 

  TYPE t_upatmo

    ! Diagnostic variables
    TYPE(t_upatmo_diag)   :: diag

    ! Tendencies
    TYPE(t_upatmo_tend)   :: tend

    ! External data
    type(t_upatmo_extdat) :: extdat

  END TYPE t_upatmo

END MODULE mo_upatmo_types
