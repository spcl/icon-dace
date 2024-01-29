! @brief Setup for synthetic satellite images
!
! configuration setup for synthetic satellite images
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

MODULE mo_synsat_config

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: max_dom
  USE mo_exception,          ONLY: message_text, finish

  IMPLICIT NONE

  PUBLIC
  PUBLIC :: get_synsat_name
  PUBLIC :: get_synsat_grib_triple

  ! named constants for better readability
  INTEGER, PARAMETER :: CHAN_IR3_9   = 1
  INTEGER, PARAMETER :: CHAN_WV6_2   = 2
  INTEGER, PARAMETER :: CHAN_WV7_3   = 3
  INTEGER, PARAMETER :: CHAN_IR8_7   = 4
  INTEGER, PARAMETER :: CHAN_IR9_7   = 5
  INTEGER, PARAMETER :: CHAN_IR10_8  = 6
  INTEGER, PARAMETER :: CHAN_IR12_1  = 7
  INTEGER, PARAMETER :: CHAN_IR13_4  = 8


  ! Code imported from the COSMO model
  ! ---------------------------------------------------------------------------- !
  ! ---------------------------------------------------------------------------- !
  ! Global (i.e. public) Declarations:

  ! 1. Global parameters from the RTTOV-module MOD_CPARAM
  !------------------------------------------------------

  !------------------------------------------------------------------------------
  ! These are global parameters which users may want to edit to optimise
  ! for their application
  !------------------------------------------------------------------------------

  INTEGER, PARAMETER ::    &
     jppf   =  1,    & ! maximal number of profiles per RTTOV call
     jpch   = 10,    & ! maximal number of channels
     jpchus = 10,    & ! maximal number of channels computed/call
     jpnsat =  9,    & ! maximal number of sensors to be used
     jplev  = 43,    & ! number of pressure levels
     jpnav  =  4,    & ! number of profile variables
     jpnsav =  5,    & ! number of surface air variables
     jpnssv =  6       ! number of skin variables

  REAL(wp), PARAMETER ::     &
     rcnv = 6.03504E5_wp,    & ! kg/kg--> ppmv ozone
     rcnw = 1.60771704E6_wp    ! kg/kg--> ppmv water vapour


  ! 2. Variables for the organisation of the satellite computations
  !----------------------------------------------------------------

  TYPE sat_check_type
    CHARACTER(LEN= 8)  :: satellite          ! platform, e.g. METEOSAT, MSG
    CHARACTER(LEN=12)  :: sensor             ! Name of sensor used
    INTEGER            :: nsat_id            ! Satellite identification within rttov
    INTEGER            :: nsat_id_min        ! for range of satellite ids
    INTEGER            :: nsat_id_max        ! for range of satellite ids
    INTEGER            :: nrttov_id          ! sensor identification within rttov
    INTEGER            :: nchan_min          ! for range of channels
    INTEGER            :: nchan_max          ! for range of channels  
    INTEGER            :: num_chan_max       ! Max. Number of channels for that sensor
  END TYPE sat_check_type

  TYPE sat_org_type
    CHARACTER(LEN= 8)  :: satellite          ! platform, e.g. METEOSAT, MSG
    CHARACTER(LEN=12)  :: sensor             ! Name of sensor used
    CHARACTER(LEN=10)  :: chan_name(jpch)    ! Names of channels used (IRx.y, WVx.y)
    INTEGER            :: nsatell_table_id   ! entry in rttov satellite table
    INTEGER            :: nsat_id            ! Satellite identification
    INTEGER            :: wmo_satid          ! WMO satellite identification
    INTEGER            :: nsensor_table_id   ! entry in rttov sensor table
    INTEGER            :: num_chan           ! Number of channels used
    INTEGER            :: nchan_list(jpch)   ! List of channels used
    INTEGER            :: ngrib_chan(4*jpch) ! list of channels for grib output
    INTEGER            :: ngrib_aees(4*jpch) ! list of additional element numbers for grib output
    REAL(wp)           :: longitude          ! position of geost. satellite
    REAL(wp)           :: emissivity(jpch)   ! emissivities for all channels
    LOGICAL            :: lclear_rad         ! for clear sky radiance 
    LOGICAL            :: lcloud_rad         ! for cloudy sky radiance 
    LOGICAL            :: lclear_tem         ! for clear sky temperature
    LOGICAL            :: lcloud_tem         ! for cloudy sky temperature
  END TYPE sat_org_type

  TYPE (sat_org_type)  :: sat_compute(jpnsat)

  ! 3. Variables for the organisation of the Namelist Input
  !--------------------------------------------------------

  TYPE sat_input_type
    CHARACTER(LEN= 8)  :: satellite         ! platform, e.g. METEOSAT, MSG
    INTEGER            :: nsat_id            ! Satellite identification
    CHARACTER(LEN=12)  :: sensor            ! Name of sensor used
    INTEGER            :: num_chan           ! Number of channels used
    LOGICAL            :: lclear_rad         ! for clear sky radiance 
    LOGICAL            :: lcloud_rad         ! for cloudy sky radiance 
    LOGICAL            :: lclear_tem         ! for clear sky temperature
    LOGICAL            :: lcloud_tem         ! for cloudy sky temperature
  END TYPE sat_input_type

  LOGICAL :: lcon_clw        ! if .TRUE.: convective liquid water used in rttov

  ! 4. Additional control variables
  !--------------------------------

  INTEGER   ::                    &
   itype_rttov,     & ! Version of RTTOV model
   nlev_rttov,      & ! Number of RTTOV levels
   num_sensors,     & ! No. of sensors used
   num_images,      & ! Total number of images
   nsat_next          ! for the organization of the computations

  INTEGER, ALLOCATABLE   :: nsat_steps(:)
         ! list of time steps for which satellite computations must be done


  ! the following are the RTTOV satellite and sensor tables 
  ! (Table 3 from the RTTOV documentation)

  CHARACTER(LEN= 8) :: rttov_satell_table(  13)
  CHARACTER(LEN=12) :: rttov_sensor_table(0:26)

  LOGICAL ::          &
    lsynsat(max_dom), & ! main switch to produce the synthetic satellite images
    lobsrad,          & ! to process satellite observations and compute radiances
    linterp             ! do interpolation of p, t, q to half-levels

  ! 5. Variables that are initialized by RTTVI
  !-------------------------------------------

  ! These variables are the same for all instruments, but are initialized
  ! by the call to RTTVI, which is done only once for the whole program

  INTEGER  ::                &
    maxknchpf,               & ! maximum number of output radiances
    total_numchans(jpnsat),  & ! Number of valid channels (not all of them are necessarily computed in ICON)
    kiu1(jpnsat)               ! for input-file unit number

  REAL(wp), ALLOCATABLE ::               &
    o3_ref     (:)          ! default ozone values on the prescribed levels

  ! These variables can be different for every instrument and are initialized
  ! by the call to RTTVI, which is done only once for the whole program

  REAL(wp), ALLOCATABLE ::               &
    ppres_d    (:,:),     & ! default pressure on the prescribed levels
    utmx       (:,:),     & ! maximum temperature
    utmn       (:,:),     & ! minimum temperature
    uqmx       (:,:),     & ! maximum humidity
    uqmn       (:,:),     & ! minimum humidity
    uomx       (:,:),     & ! maximum ozone
    uomn       (:,:),     & ! minimum ozone
    ivch       (:,:)

  ! 6. Some constant variables
  !---------------------------

  REAL  (KIND=wp) ::           &
    const_aheat, r_sat

  ! for output of MSG-variables
  INTEGER , PARAMETER  :: nmsgchan = 8


  ! 7. Data structures and Variables for RTTOV9 and higher
  !-------------------------------------------------------

  ! Required for initialization of RTTOV9/10
  ! (function rttov_init of mo_rttov_ifc)

  INTEGER, ALLOCATABLE   :: instruments(:,:)
    ! for every sensor (2nd dimension):
    !    nsatell_table_id
    !    nsat_id
    !    nsensor_table_id

  INTEGER, ALLOCATABLE   :: total_channels(:,:)
    ! for every sensor (2nd dimension):
    !    wmo_satid
    !    nsatell_table_id

  INTEGER, ALLOCATABLE   :: total_n_chans(:)
    ! for every sensor the number of channels

  LOGICAL, ALLOCATABLE   :: addclouds(:)
    ! for every sensor if it shall use ir cloud scattering

  INTEGER                :: mchans
    ! maximum number of channels for one sensor (MAX (total_numchans))

  ! maximum satellite zenith angles
  REAL(wp), PARAMETER :: zenmax9  = 86.5_wp        ! deg
  REAL(wp), PARAMETER :: zenmax10 = 75.0_wp        ! deg


  REAL(wp), PARAMETER ::            &
    p_top = 9.9_wp,                 & ! Highest pressure level (in Pa) that is used for input to RTTOV
    t_top = 231.6_wp,               & ! Standard temperature  [K]
    w_top = 0.349555E-05_wp,        & ! Standard mixing ratio [Kg/Kg]
    q_top = w_top / (1._wp + w_top)   ! all at highest pressure level

  ! Bits for extrapolation of RTTOV input profiles above model top:
  INTEGER, PARAMETER :: extrp_logp    = 1 !extrapolate log(p) linearly instead of p
  INTEGER, PARAMETER :: extrp_const   = 2 !extrapolate t, q with constant values
  INTEGER, PARAMETER :: extrp_lin     = 4 !extrapolate t, q linearly
  INTEGER, PARAMETER :: extrp_clim    = 8 !extrapolate t, q to climatological value
  INTEGER            :: extrp_type    = extrp_const

  INTEGER :: iwc2effdiam   = 4
     ! 1: Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
     ! 2: Scheme by Wyser et al. (see McFarquhar et al. (2003))
     ! 3: Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
     ! 4: Scheme by McFarquhar et al. (2003)
     !RF  Don't use 1-3, these values might cause floating point exceptions!!
  INTEGER :: iceshape = 1
     ! 1: hexagonal
     ! 2: ice aggregates

  ! 8. Namelist variables for reading SEVIRI NWCSAF cloud products
  !---------------------------------------------------------------------------
  LOGICAL             :: lread_ct              ! namelist switch: read cloud type?
  CHARACTER (LEN=100) :: clouddir             ! directory of cloud data

  !==============================================================================


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_synsat_config'

  CONTAINS


  !! Configuration for synthetic satellite images
  !!
  SUBROUTINE configure_synsat
    CHARACTER(LEN=*), PARAMETER    :: routine = modname//"::configure_synsat"

    INTEGER :: i

    total_numchans(1) = 8
    total_numchans(2:) = 0

    num_sensors = 1

    sat_compute(1)%satellite         = 'MSG'
    sat_compute(1)%nsat_id            = 2
    sat_compute(1)%wmo_satid          = 56
    sat_compute(1)%longitude          = 0.
    sat_compute(1)%sensor            = 'SEVIRI'
    sat_compute(1)%num_chan           = 8
    sat_compute(1)%nchan_list(1:total_numchans(1)) = (/ (i, i=1,8) /)
    sat_compute(1)%chan_name(1:total_numchans(1))  = &
                (/ 'IR3.9 ', &
                   'WV6.2 ', &
                   'WV7.3 ', &
                   'IR8.7 ', &
                   'IR9.7 ', &
                   'IR10.8', &
                   'IR12.1', &
                   'IR13.4' /)
    sat_compute(1)%emissivity(1:jpch) = 0.
    sat_compute(1)%lclear_rad         = .true.
    sat_compute(1)%lcloud_rad         = .true.
    sat_compute(1)%lclear_tem         = .true.
    sat_compute(1)%lcloud_tem         = .true.

    num_images = 4*SUM(total_numchans(:))
    mchans = MAXVAL(total_numchans)

    ALLOCATE(instruments(3,num_sensors),addclouds(num_sensors))
    ALLOCATE(total_channels(mchans, num_sensors))
    ALLOCATE(total_n_chans(num_sensors))

    instruments(1:3,1) = (/ 12, 2,21/)
    total_channels(1:8,1) = (/ (i, i=1,8) /)
    total_n_chans(:) = 8
    addclouds(:) = .TRUE.

    ! --- consistency check
    
    ! Since "n_chans", "numchans" are used interchangeably (more or
    ! less) under the implicit assumption that we have only one
    ! satellite, we make place an assertion here:
    IF (ANY(total_numchans(2:) > 0) .OR. (num_sensors > 1)) THEN
      CALL finish(routine, "RTTOV interface not implemented for more than one satellite!")
    END IF

  END SUBROUTINE configure_synsat


  !> Auxiliary routine: Generate GRIB2 shortname for a given synthetic
  !> satellite image product.
  !
  SUBROUTINE get_synsat_name(lradiance, lcloudy, ichan, &
    &                        shortname, longname)
    LOGICAL, INTENT(IN) :: lradiance        ! shortname for radiance or brightness temp
    LOGICAL, INTENT(IN) :: lcloudy          ! shortname for cloudy or clear sky
    INTEGER, INTENT(IN) :: ichan            ! channel index
    CHARACTER(LEN=*), INTENT(OUT) :: shortname ! result: shortname
    CHARACTER(len=*), INTENT(OUT) :: longname  ! result: long name

    shortname = "SYNMSG_"
    IF (lradiance) THEN
      shortname = TRIM(shortname)//"RAD_"
    ELSE
      shortname = TRIM(shortname)//"BT_"
    END IF
    IF (lcloudy) THEN
      shortname = TRIM(shortname)//"CL_"
    ELSE
      shortname = TRIM(shortname)//"CS_"
    END IF
    shortname = TRIM(shortname)//sat_compute(1)%chan_name(ichan)

    longname = "synth. sat."
    IF (lradiance) THEN
      longname = TRIM(longname)//" radiance"
    ELSE
      longname = TRIM(longname)//" brightness temperature"
    END IF
    IF (lcloudy) THEN
      longname = TRIM(longname)//" cloudy"
    ELSE
      longname = TRIM(longname)//" clear sky"
    END IF
  END SUBROUTINE get_synsat_name


  !> Auxiliary routine: Returns GRIB2 category, discipline and number
  !> for a given synthetic satellite image product.
  !
  SUBROUTINE get_synsat_grib_triple(lradiance, lcloudy, ichan,       &
    &                               idiscipline, icategory, inumber, &
    &                               scaledValueOfCentralWaveNumber,  &
    &                               scaleFactorOfCentralWaveNumber)
    LOGICAL, INTENT(IN)  :: lradiance        ! shortname for radiance or brightness temp
    LOGICAL, INTENT(IN)  :: lcloudy          ! shortname for cloudy or clear sky
    INTEGER, INTENT(IN)  :: ichan            ! channel index
    INTEGER, INTENT(OUT) :: idiscipline, icategory, inumber ! result: GRIB2 triple
    INTEGER, INTENT(OUT) :: scaledValueOfCentralWaveNumber  ! result: GRIB2 meta-data
    INTEGER, INTENT(OUT) :: scaleFactorOfCentralWaveNumber  ! result: GRIB2 meta-data
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::get_synsat_grib_triple'

    ! common for all products
    idiscipline =  3
    icategory   =  1
    IF (lradiance) THEN
      IF (lcloudy) THEN
        inumber = 16
      ELSE
        inumber = 17
      END IF
    ELSE
      IF (lcloudy) THEN
        inumber = 14
      ELSE
        inumber = 15
      END IF
    END IF

    IF (lradiance) THEN
      ! radiances
      IF (lcloudy) THEN
        SELECT CASE(ichan)
        CASE (CHAN_IR3_9)
          scaledValueOfCentralWaveNumber = 256410
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV6_2)
          scaledValueOfCentralWaveNumber = 161290
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV7_3)
          scaledValueOfCentralWaveNumber = 136986
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR8_7)
          scaledValueOfCentralWaveNumber = 114942
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR9_7) 
          scaledValueOfCentralWaveNumber = 103092
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR10_8)
          scaledValueOfCentralWaveNumber =  92592
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR12_1)
          scaledValueOfCentralWaveNumber =  82644
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR13_4)
          scaledValueOfCentralWaveNumber =  74626
          scaleFactorOfCentralWaveNumber = 0
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT
      ELSE
        SELECT CASE(ichan)
        CASE (CHAN_IR3_9)
          scaledValueOfCentralWaveNumber = 256410
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV6_2)
          scaledValueOfCentralWaveNumber = 161290
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV7_3)
          scaledValueOfCentralWaveNumber = 136986
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR8_7)
          scaledValueOfCentralWaveNumber = 114942
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR9_7) 
          scaledValueOfCentralWaveNumber = 103092
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR10_8)
          scaledValueOfCentralWaveNumber =  92592
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR12_1)
          scaledValueOfCentralWaveNumber =  82644
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR13_4)
          scaledValueOfCentralWaveNumber =  74626
          scaleFactorOfCentralWaveNumber = 0
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT
      END IF
    ELSE
      ! brightness temperatures
      IF (lcloudy) THEN
        SELECT CASE(ichan)
        CASE (CHAN_IR3_9)
          scaledValueOfCentralWaveNumber = 256410
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV6_2)
          scaledValueOfCentralWaveNumber = 161290
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV7_3)
          scaledValueOfCentralWaveNumber = 136986
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR8_7)
          scaledValueOfCentralWaveNumber = 114942
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR9_7) 
          scaledValueOfCentralWaveNumber = 103092
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR10_8)
          scaledValueOfCentralWaveNumber =  92592
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR12_1)
          scaledValueOfCentralWaveNumber =  82644
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR13_4)
          scaledValueOfCentralWaveNumber =  74626
          scaleFactorOfCentralWaveNumber = 0
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT
      ELSE
        SELECT CASE(ichan)
        CASE (CHAN_IR3_9)
          scaledValueOfCentralWaveNumber = 256410
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV6_2)
          scaledValueOfCentralWaveNumber = 161290
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_WV7_3)
          scaledValueOfCentralWaveNumber = 136986
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR8_7)
          scaledValueOfCentralWaveNumber = 114942
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR9_7) 
          scaledValueOfCentralWaveNumber = 103092
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR10_8)
          scaledValueOfCentralWaveNumber =  92592
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR12_1)
          scaledValueOfCentralWaveNumber =  82644
          scaleFactorOfCentralWaveNumber = 0
        CASE (CHAN_IR13_4)
          scaledValueOfCentralWaveNumber =  74626
          scaleFactorOfCentralWaveNumber = 0
        CASE DEFAULT
          CALL finish(routine, "Internal error!")
        END SELECT
      END IF
    END IF
  END SUBROUTINE get_synsat_grib_triple

END MODULE mo_synsat_config
