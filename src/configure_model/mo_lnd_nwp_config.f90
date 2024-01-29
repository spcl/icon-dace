! @brief configuration setup for NWP land scheme TERRA
!
! configuration setup for NWP land scheme TERRA
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

MODULE mo_lnd_nwp_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: GLOBCOVER2009, GLC2000, vname_len
  USE mo_io_units,           ONLY: filename_max
  USE mo_nwp_sfc_tiles,      ONLY: t_tile_list, setup_tile_list
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_coupling_config,    ONLY: is_coupled_to_ocean


  IMPLICIT NONE

  PRIVATE

  ! FUNCTIONS/SUBROUTINES
  PUBLIC :: configure_lnd_nwp
  PUBLIC :: convert_luc_ICON2GRIB


  ! VARIABLES
  PUBLIC :: dzsoil, zml_soil, nlev_soil, nlev_snow, ibot_w_so, ntiles_total, ntiles_lnd, ntiles_water
  PUBLIC :: frlnd_thrhld, frlndtile_thrhld, frlake_thrhld, frsea_thrhld, frsi_min, hice_min, hice_max
  PUBLIC :: lseaice, lprog_albsi, lbottom_hflux, llake, lmelt, lmelt_var, lmulti_snow, lsnowtile, max_toplaydepth
  PUBLIC :: itype_trvg, itype_evsl, itype_lndtbl, l2lay_rho_snow
  PUBLIC :: itype_root, itype_heatcond, itype_interception, &
            itype_hydbound, idiag_snowfrac, itype_snowevap, cwimax_ml, c_soil, c_soil_urb, cr_bsmin
  PUBLIC :: itype_canopy, cskinc, tau_skin
  PUBLIC :: lterra_urb, lurbalb, lurbahf, itype_kbmo, itype_eisa
  PUBLIC :: lstomata, l2tls, lana_rho_snow
  PUBLIC :: isub_water, isub_lake, isub_seaice
  PUBLIC :: sstice_mode, sst_td_filename, ci_td_filename
  PUBLIC :: tile_list
  PUBLIC :: groups_smi
  PUBLIC :: czbot_w_so


  !--------------------------------------------------------------------------
  ! Basic configuration setup for NWP land
  !--------------------------------------------------------------------------

!  TYPE t_nwp_lnd_config

  ! namelist variables
  INTEGER ::  nlev_snow          !< number of snow layers
  INTEGER ::  ntiles_lnd         !< number of static land surface types
  REAL(wp)::  frlnd_thrhld       !< fraction threshold for creating a land grid point
  REAL(wp)::  frlndtile_thrhld   !< fraction threshold for retaining the respective 
                                 !< tile for a grid point
  REAL(wp)::  frlake_thrhld      !< fraction threshold for creating a lake grid point
  REAL(wp)::  frsea_thrhld       !< fraction threshold for creating a sea grid point
  REAL(wp)::  frsi_min           !< minimum sea-ice fraction  [-]
  REAL(wp)::  hice_min           !< minimum sea-ice thickness [m]
  REAL(wp)::  hice_max           !< maximum sea-ice thickness [m]
  LOGICAL ::  lbottom_hflux      !< use simple parameterization for heat flux through sea ice bottom
  INTEGER ::  itype_trvg         !< type of vegetation transpiration parameterization
  INTEGER ::  itype_evsl         !< type of parameterization of bare soil evaporation (see Schulz and Vogel 2020)
  INTEGER ::  itype_lndtbl       !< choice of table for associating surface parameters to land-cover classes
  INTEGER ::  itype_root         !< type of root density distribution
  INTEGER ::  itype_heatcond     !< type of soil heat conductivity (see Schulz et al. 2016)
  INTEGER ::  itype_interception !< type of plant interception
  !$ACC DECLARE CREATE(itype_interception)
  REAL(wp)::  cwimax_ml          !< scaling parameter for maximum interception storage
  REAL(wp)::  c_soil             !< surface area density of the (evaporative) soil surface
  REAL(wp)::  c_soil_urb         !< surface area density of the (evaporative) soil surface, urban areas
  REAL(wp)::  cr_bsmin           !< minimum bare soil evaporation resistance (see Schulz and Vogel 2020)
  INTEGER ::  itype_canopy       !< type of canopy parameterisation with respect to the surface energy balance
                                 !< (see Schulz and Vogel 2020)
  REAL(wp)::  cskinc             !< skin conductivity (W/m**2/K)
  REAL(wp)::  tau_skin           !< relaxation time scale for the computation of the skin temperature
  LOGICAL ::  lterra_urb         !< activate urban model TERRA_URB (see Schulz et al. 2023)
  LOGICAL ::  lurbalb            !< use urban albedo and emissivity
  LOGICAL ::  lurbahf            !< use urban anthropogenic heat flux
  INTEGER ::  itype_kbmo         !< type of bluff-body thermal roughness length parameterisation
  INTEGER ::  itype_eisa         !< type of evaporation from impervious surface area
  INTEGER ::  itype_hydbound     !< type of hydraulic lower boundary condition
  INTEGER ::  idiag_snowfrac     !< method for diagnosis of snow-cover fraction
  INTEGER ::  itype_snowevap     !< treatment of snow evaporation in the presence of vegetation      

  LOGICAL ::  lseaice     !> forecast with sea-ice model
  LOGICAL ::  lprog_albsi !> sea-ice albedo is computed prognostically from a rate equation
  LOGICAL ::  llake       !! forecast with lake model FLake
  LOGICAL ::  lmelt       !! soil model with melting process
  LOGICAL ::  lmelt_var   !! freezing temperature dependent on water content
  LOGICAL ::  lmulti_snow !! run the multi-layer snow model
  LOGICAL ::  l2lay_rho_snow !! use two-layer snow density for single-layer snow model
  REAL(wp)::  max_toplaydepth !< maximum depth of uppermost snow layer for multi-layer snow scheme
  LOGICAL ::  lstomata    !! map of minimum stomata resistance
  LOGICAL ::  l2tls       !! forecast with 2-TL integration scheme
  LOGICAL ::  lana_rho_snow !! if .TRUE., take rho_snow-values from analysis file 
  LOGICAL ::  lsnowtile   !! if .TRUE., snow is considered as a separate tile
 
  INTEGER ::  sstice_mode      !< set if SST and sea ice cover are read from the analysis
                               !< and kept constant or read from external data files 
                               !< and updated regularly in run time
  REAL(wp)::  czbot_w_so       !< thickness of the hydraulical active soil layer [m]

  CHARACTER(LEN=filename_max) :: sst_td_filename, ci_td_filename

  ! derived variables
  INTEGER ::  ibot_w_so    !< number of hydrological active soil layers 
  INTEGER ::  isub_water   !< (open) water points tile number
  INTEGER ::  isub_lake    !< lake points tile number
  INTEGER ::  isub_seaice  !< seaice tile number
  INTEGER ::  ntiles_total !< total number of TILES
  INTEGER ::  ntiles_water !< number of extra tiles for ocean, seaice and lakes
  INTEGER ::  nlev_soil    !< number of soil layers
  REAL(wp), ALLOCATABLE :: zml_soil(:)   !< soil layer full level heights read in and used by ICON
  REAL(wp), ALLOCATABLE :: dzsoil(:)     !< soil layer thickness
  REAL(wp), ALLOCATABLE :: depth_hl(:)   !< depths of half levels

!  END TYPE t_nwp_lnd_config

  !>
  !!
!  TYPE(t_nwp_lnd_config) :: nwp_lnd_config(max_dom)



   TYPE(t_tile_list), TARGET :: tile_list  ! list of tiles


   CHARACTER(LEN = *), PARAMETER :: modname = "mo_lnd_nwp_config"

   CHARACTER(LEN=vname_len),DIMENSION(1) :: groups_smi = (/"mode_iniana"/)

CONTAINS

  !! setup components of the NWP land scheme
  !!
  !! Setup of additional nwp-land control variables depending on the 
  !! land-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  SUBROUTINE configure_lnd_nwp(ljsbach)
  !
    LOGICAL, INTENT(IN) :: ljsbach !< JSBACH is in use. Don't adjust land/lake/sea thresholds.
    ! local variables
    INTEGER  :: js              ! soil loop index

    CHARACTER(len=*), PARAMETER::  &
      &  routine = modname//"::configure_lnd_nwp"

    !-----------------------------------------------------------------------

    ! calculate the soil layer thickness based on depth of the soil full layers

    ALLOCATE(dzsoil(nlev_soil))
    ALLOCATE(depth_hl(nlev_soil))     ! half level depth

    depth_hl(1) = 2.0_wp*zml_soil(1)  ! depth of first half level
    dzsoil(1)=depth_hl(1)             ! layer thickness betw. half levels of uppermost layer
    DO js=2,nlev_soil
      ! depth of half levels
      depth_hl(js) = depth_hl(js-1) + 2.0_wp*(zml_soil(js) - depth_hl(js-1))
      ! layer thickness betw. half levels
      dzsoil(js)   = depth_hl(js) - depth_hl(js-1)
    ENDDO
    !$ACC ENTER DATA COPYIN(dzsoil)

    ! Determine ibot_w_so
    DO js=1,nlev_soil
      IF (depth_hl(js) <= czbot_w_so) ibot_w_so=js
    ENDDO
    ! make sure that 2 <= ibot_w_so <= nlev_soil-1
    ibot_w_so = MAX(2, MIN(ibot_w_so,nlev_soil-1))
    !
    WRITE(message_text,'(a,I2)') 'Number of hydrological active soil layers ibot_w_so is set to: ',ibot_w_so
    CALL message(TRIM(routine),message_text)

    ! seaice fraction limit
    IF ( is_coupled_to_ocean() ) THEN
       frsi_min = 1.0E-10_wp  ! ICON coupled with ocean, epsilon because ICON-O determines seaice fraction
    ELSE
       frsi_min = 0.015_wp    ! ICON uncoupled, limit at 1.5% seaice fraction (Dmitrii Mironov)
    END IF
 
    !
    ! settings dealing with surface tiles
    !

    IF (ljsbach) THEN
      ! switch off snowtiles
      lsnowtile = .FALSE.
    ELSE  ! TERRA
      ! when running without land tiles, switch off snowtiles
      ! and reset thresholds
      IF (ntiles_lnd==1) THEN
        lsnowtile     = .FALSE.
        frlnd_thrhld  = 0.5_wp
        frlake_thrhld = 0.5_wp
        frsea_thrhld  = 0.5_wp
      ENDIF
    ENDIF

    ! number of tiles; ntiles_lnd is set by the namelist variable ntiles
    IF(lsnowtile) THEN
      ntiles_total = 2*ntiles_lnd
    ELSE
      ntiles_total = ntiles_lnd
    END IF

    IF (ntiles_total == 1) THEN 
      ! no tile approach, thus no extra tile index for water points
      ntiles_water = 0
    ELSE
      ! extra tiles for water points
      ntiles_water = 3 ! one for open sea points
                       ! another one for lake points
                       ! another one for seaice
    ENDIF

    ! (open) water points tile number
    isub_water  = MAX(1,ntiles_total + ntiles_water - 2)

    ! lake points tile number
    isub_lake   = MAX(1,ntiles_total + ntiles_water - 1)

    ! sea-ice tile number
    isub_seaice = ntiles_total + ntiles_water


    ! this results in the following model internal tile structure
    !
    ! ASCII art for ICON-internal tile structure, depicting the configuration for:
    ! - ntiles_lnd = 3
    ! - lsnowtile  = .TRUE.
    !
    !------------------------------------------------------------------!
    !   1  |   2  |   3  ||   4  |   5  |   6  ||   7  |   8  |    9   !
    ! land | land | land || snow | snow | snow || open | lake | seaice !
    !      |      |      ||      |      |      || sea  |      |        !
    !------------------------------------------------------------------!
    !
    !<---  ntiles_lnd --->                      <---- ntiles_water ---->
    !
    !<-------------- ntiles_total ------------>
    !
    !<---------------- ntiles_total + ntiles_water ------------------->


    ! Setup tile meta-information required for tile I/O.
    !
    ! I.e. the mapping between the internal tile structure and the output 
    ! structure (according to GRIB2 Product Definition Template (PDT) 4.55) 
    ! is set up 
    CALL setup_tile_list (tile_list, ntiles_lnd, lsnowtile, isub_water, isub_lake, isub_seaice)


  END SUBROUTINE configure_lnd_nwp


  !! Given the internal land use class index, provide the official GRIB2 index.
  !!
  !! Given the ICON-internal land use class index, provide the official GRIB2 
  !! index according to GRIB2 table 4.2.43 (or vice versa).
  !!
  FUNCTION convert_luc_ICON2GRIB(lc_datbase,iluc_in,opt_linverse)  RESULT (iluc_out)

    INTEGER          , INTENT(IN) :: lc_datbase       !< landuse class database
    INTEGER          , INTENT(IN) :: iluc_in          !< landuse class for grid point i
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_linverse     !< given the GRIB2 landsuse class, 
                                                      !< provide the internal one

    ! Local
    LOGICAL :: linverse
    INTEGER :: iluc_out          !< landuse class index output
                                 !< Default: according to GRIB2 table 4.243

    INTEGER :: iluc_GLOBCOVER2009(23)  ! num_lcc
    INTEGER :: iluc_GLC2000(23)        ! num_lcc
    !
    ! maps GLOBCOVER2009 landuse classes onto tile class table 4.243
    DATA iluc_GLOBCOVER2009 / 24, &  ! 1 :irrigated croplands
                            & 25, &  ! 2 :rainfed croplands 
                            & 26, &  ! 3 :mosaic cropland (50-70%) - vegetation (20-50%)
                            & 27, &  ! 4 :mosaic vegetation (50-70%) - cropland (20-50%)
                            & 28, &  ! 5 :closed broadleaved evergreen forest           
                            & 2 , &  ! 6 :closed broadleaved deciduous forest           
                            & 3 , &  ! 7 :open broadleaved deciduous forest             
                            & 29, &  ! 8 :closed needleleaved evergreen forest          
                            & 30, &  ! 9 :open needleleaved deciduous forest            
                            & 31, &  ! 10:mixed broadleaved and needleleaved forest     
                            & 32, &  ! 11:mosaic shrubland (50-70%) - grassland (20-50%)
                            & 33, &  ! 12:mosaic grassland (50-70%) - shrubland (20-50%)
                            & 34, &  ! 13:closed to open shrubland                      
                            & 25, &  ! 14:closed to open herbaceous vegetation          
                            & 35, &  ! 15:sparse vegetation                             
                            & 36, &  ! 16:closed to open forest regularly flooded        
                            & 37, &  ! 17:closed forest or shrubland permanently flooded
                            & 38, &  ! 18:closed to open grassland regularly flooded    
                            & 22, &  ! 19:artificial surfaces                           
                            & 19, &  ! 20:bare areas                                     
                            & 20, &  ! 21:water bodies                                   
                            & 21, &  ! 22:permanent snow and ice                         
                            & 39  /  ! 23:undefined
    !
    ! maps GLC2000 landuse classes onto tile class table 4.243
    DATA iluc_GLC2000       / 1 , &  ! 1 :evergreen broadleaf forest
                            & 2 , &  ! 2 :deciduous broadleaf closed forest 
                            & 3 , &  ! 3 :deciduous broadleaf open forest
                            & 4 , &  ! 4 :evergreen needleleaf forest
                            & 5 , &  ! 5 :deciduous needleleaf forest           
                            & 6 , &  ! 6 :mixed leaf trees           
                            & 7 , &  ! 7 :fresh water flooded trees             
                            & 8 , &  ! 8 :saline water flooded trees          
                            & 9 , &  ! 9 :mosaic tree / natural vegetation            
                            & 10, &  ! 10:burnt tree cover     
                            & 11, &  ! 11:evergreen shrubs closed-open
                            & 12, &  ! 12:decidous shrubs closed-open
                            & 13, &  ! 13:herbaceous vegetation closed-open                      
                            & 14, &  ! 14:sparse herbaceous or grass          
                            & 15, &  ! 15:flooded shrubs or herbaceous                             
                            & 16, &  ! 16:cultivated & managed areas        
                            & 17, &  ! 17:mosaic crop / tree / natural vegetation
                            & 18, &  ! 18:mosaic crop / shrub / grass    
                            & 19, &  ! 19:bare areas                    
                            & 20, &  ! 20:water                                     
                            & 21, &  ! 21:snow & ice                                   
                            & 22, &  ! 22:artificial surface                         
                            & 39  /  ! 23:undefined
    !-----------------------------------------------------------------------

    IF (PRESENT(opt_linverse)) THEN
      linverse = opt_linverse
    ELSE
      linverse = .FALSE.
    ENDIF


    ! Special treatment for undefined points
    IF (iluc_in == -1) THEN
      iluc_out = 39                !undefined
      RETURN
    ENDIF
    IF (iluc_in == 0) THEN
      iluc_out = -9999       ! warning
      RETURN
    ENDIF


    SELECT CASE(lc_datbase)
    CASE(GLOBCOVER2009)
      IF (.NOT. linverse) THEN
        ! internal to GRIB2
        iluc_out = iluc_GLOBCOVER2009(iluc_in)
      ELSE
        ! GRIB2 to internal
!!$        iluc_out = find_loc(iluc_in,iluc_GLOBCOVER2009)
      ENDIF

    CASE(GLC2000)

      IF (.NOT. linverse) THEN
        ! internal to GRIB2
        iluc_out = iluc_GLC2000(iluc_in)
      ELSE
        ! GRIB2 to internal
!!$        iluc_out = find_loc(iluc_in,iluc_GLC2000)
      ENDIF

    END SELECT 


  END FUNCTION convert_luc_ICON2GRIB

END MODULE mo_lnd_nwp_config
