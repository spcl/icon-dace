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

! Namelist setup for waves.
! The content of namelists is mostly adopted from the WAM 4.5.4.

MODULE mo_wave_nml

  USE mo_kind,                ONLY: wp
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_exception,           ONLY: finish
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_wave_config,         ONLY: wave_config


  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: read_wave_namelist

CONTAINS
  !>
  !!
  SUBROUTINE read_wave_namelist(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, iunit, funit
    INTEGER :: jg

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_wave_nml:read_wave_namelist'

    INTEGER  :: ndirs    ! NUMBER OF DIRECTIONS.
    INTEGER  :: nfreqs   ! NUMBER OF FREQUENCIES.
    INTEGER  :: IREF     ! FREQUENCY BIN NUMBER OF REFERENCE FREQUENCY

    REAL(wp) :: fr1      ! FIRST FREQUENCY [HZ].
    REAL(wp) :: CO       ! FREQUENCY RATIO

    REAL(wp) :: ALPHA      ! PHILLIPS' PARAMETER  (NOT USED IF IOPTI = 1)
    REAL(wp) :: FM         ! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY
    REAL(wp) :: GAMMA_wave ! OVERSHOOT FACTOR
    REAL(wp) :: SIGMA_A    ! LEFT PEAK WIDTH
    REAL(wp) :: SIGMA_B    ! RIGHT PEAK WIDTH
    REAL(wp) :: FETCH      ! FETCH IN METRES (IF ZERO THEN 0.5 OF THE LATITUDE INCREMENT IS USED.).

    REAL(wp) :: roair   ! AIR DENSITY
    REAL(wp) :: RNUAIR  ! KINEMATIC AIR VISCOSITY
    REAL(wp) :: RNUAIRM ! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER
    REAL(wp) :: ROWATER ! WATER DENSITY
    REAL(wp) :: XEPS
    REAL(wp) :: XINVEPS

    REAL(wp) :: ALPHA_CH ! minimum charnock constant (ecmwf cy45r1).
                         ! 0.0060, if le 30 frequencies changed !@waves todo
                         ! to 0.0075 in subroutine initmdl !@waves todo

    REAL(wp) :: depth     ! ocean depth (m) if not 0, then constant depth
    REAL(wp) :: depth_min ! allowed minimum of model depth (m)
    REAL(wp) :: depth_max ! allowed maximum of model depth (m)

    INTEGER  :: niter_smooth ! number of smoothing iterations for wave bathymetry

    REAL(wp) :: XKAPPA  ! VON KARMAN CONSTANT.
    REAL(wp) :: XNLEV   ! WINDSPEED REF. LEVEL.
    REAL(wp) :: BETAMAX ! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
    REAL(wp) :: ZALP    ! SHIFTS GROWTH CURVE (ECMWF CY45R1).
    !  REAL(wp)  :: ALPHA      ! MINIMUM CHARNOCK CONSTANT (ECMWF CY45R1).
    ! if LE 30 frequencies changed
    ! to 0.0075 in subroutine INITMDL

    INTEGER :: jtot_tauhf          ! dimension of wave_config%wtauhf. it must be odd !!!

    CHARACTER(LEN=filename_max) :: forc_file_prefix ! prefix of forcing file name
                                           ! the real file name will be constructed as:
                                           ! forc_file_prefix+'_wind' for U and V 10 meter wind (m/s)
                                           ! forc_file_prefix+'_ice'  for sea ice concentration (fraction of 1)
                                           ! forc_file_prefix+'_slh'  for sea level height (m)
                                           ! forc_file_prefix+'_osc'  for U and V ocean surface currents (m/s)

    LOGICAL :: linput_sf1      ! if .TRUE., calculate wind input source function term, first call
    LOGICAL :: linput_sf2      ! if .TRUE., calculate wind input source function term, second call
    LOGICAL :: ldissip_sf      ! if .TRUE., calculate dissipation source function term
    LOGICAL :: lwave_brk_sf    ! if .TRUE., calculate wave breaking dissipation source function term
    LOGICAL :: lnon_linear_sf  ! if .TRUE., calculate non linear source function term
    LOGICAL :: lbottom_fric_sf ! if .TRUE., calculate bottom_friction source function term
    LOGICAL :: lwave_stress1   ! if .TRUE., calculate wave stress, first call
    LOGICAL :: lwave_stress2   ! if .TRUE., calculate wave stress, second call

    ! for test case
    REAL(wp) :: peak_u10, peak_v10 ! peak values (m/s) of 10 m U and V wind speed for test case
    REAL(wp) :: peak_lat, peak_lon ! geographical location (deg) of wind peak value

    ! source function time integration
    REAL(wp) :: impl_fac       ! implicitness factor for total source function time integration
                               ! impl_fac=0.5 : second order Crank-Nicholson/trapezoidal scheme
                               ! impl_fac=1   : first order Euler backward scheme
                               ! valid range: 0.5 <= impl_fac <= 1

    NAMELIST /wave_nml/ &
         forc_file_prefix,          &
         ndirs, nfreqs, fr1, CO, IREF,                      &
         ALPHA, FM, GAMMA_wave, SIGMA_A, SIGMA_B, FETCH,    &
         roair, RNUAIR, RNUAIRM, ROWATER, XEPS, XINVEPS, &
         XKAPPA, XNLEV, BETAMAX, ZALP, jtot_tauhf, ALPHA_CH, &
         depth, depth_min, depth_max, niter_smooth, &
         linput_sf1, linput_sf2, ldissip_sf, lwave_brk_sf, lnon_linear_sf, lbottom_fric_sf, &
         lwave_stress1, lwave_stress2, peak_u10, peak_v10, peak_lat, peak_lon, &
         impl_fac

    !-----------------------------------------------------------
    ! 1. default settings
    !-----------------------------------------------------------
    ndirs      = 24             !! NUMBER OF DIRECTIONS.
    nfreqs     = 25             !! NUMBER OF FREQUENCIES.
    fr1        = 0.04177248_wp  !! FIRST FREQUENCY [HZ].
    CO         = 1.1_wp         !! FREQUENCY RATIO
    IREF       = 1              !! FREQUENCY BIN NUMBER OF REFERENCE FREQUENCY

    ALPHA      = 0.018_wp       !! PHILLIPS PARAMETER.
    FM         = 0.2_wp         !! PEAK FREQUENCY (HZ) AND/OR MAXIMUM FREQUENCY.
    GAMMA_wave = 3.0_wp         !! OVERSHOOT FACTOR.
    SIGMA_A    = 0.07_wp        !! LEFT PEAK WIDTH.
    SIGMA_B    = 0.09_wp        !! RIGHT PEAK WIDTH.

    FETCH      = 300000._wp     !! FETCH IN METRES.

    roair      = 1.225_wp       !! AIR DENSITY
    RNUAIR     = 1.5E-5_wp      !! KINEMATIC AIR VISCOSITY
    RNUAIRM    = 0.11_wp*RNUAIR !! KINEMATIC AIR VISCOSITY FOR MOMENTUM TRANSFER

    ROWATER    = 1000._wp       !! WATER DENSITY
    XEPS       = roair/ROWATER
    XINVEPS    = 1._wp/XEPS

    BETAMAX    = 1.20_wp        !! PARAMETER FOR WIND INPUT (ECMWF CY45R1).
    ZALP       = 0.0080_wp      !! SHIFTS GROWTH CURVE (ECMWF CY45R1).
    jtot_tauhf = 19             !! dimension of wtauhf. it must be odd
    ALPHA_CH   = 0.0075_wp      !! minimum charnock constant (ecmwf cy45r1).

    depth        = 0._wp        !! ocean depth (m) if not 0, then constant depth
    depth_min    = 0.2_wp       !! allowed minimum of model depth (m)
    depth_max    = 999.0_wp     !! allowed maximum of model depth (m)
    niter_smooth = 1            !! number of smoothing iterations for wave bathymetry
                                !! if 0 then no smoothing

    XKAPPA     = 0.40_wp        !! VON KARMAN CONSTANT.
    XNLEV      = 10.0_wp        !! WINDSPEED REF. LEVEL.

    forc_file_prefix = ''

    linput_sf1 =       .TRUE. !< if .TRUE., calculate wind input source function term, first call
    linput_sf2 =       .TRUE. !< if .TRUE., calculate wind input source function term, second call
    ldissip_sf =       .TRUE. !< if .TRUE., calculate dissipation source function term
    lwave_brk_sf =     .TRUE. !< if .TRUE., calculate wave breaking dissipation source function term
    lnon_linear_sf =   .TRUE. !< if .TRUE., calculate non linear source function term
    lbottom_fric_sf =  .TRUE. !< if .TRUE., calculate bottom_friction source function term
    lwave_stress1  =   .TRUE. !< if .TRUE., calculate wave stress, first call
    lwave_stress2  =   .TRUE. !< if .TRUE., calculate wave stress, second call

    peak_u10 = 9.0_wp         !! peak value (m/s) of 10 m U wind component for test case
    peak_v10 = 9.0_wp         !! peak value (m/s) of 10 m V wind component for test case
    peak_lat = -60.0_wp       !! latitude (deg) of wind peak value
    peak_lon = -140.0_wp      !! longitude (deg) of wind peak value

    impl_fac = 1.0_wp         !! first order Euler backward time integration scheme
                              !! for total source function

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('wave_nml')
      READ(funit,NML=wave_nml)
      CALL close_tmpfile(funit)
    END IF


    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, wave_nml)   ! write defaults to temporary text file
    END IF

    CALL open_nml(TRIM(filename))
    CALL position_nml ('wave_nml', STATUS=istat)

    SELECT CASE (istat)
    CASE (positioned)
       READ (nnml, wave_nml)                                        ! overwrite default settings
       IF (my_process_is_stdio()) THEN
          iunit = temp_settings()
          WRITE(iunit, wave_nml)    ! write settings to temporary text file
       END IF
    END SELECT
    CALL close_nml


    !----------------------------------------------------
    ! 4. Sanity checks
    !----------------------------------------------------

    IF (MOD(jtot_tauhf,2).eq.0) THEN
      CALL finish(TRIM(routine),'Error: jtot_tauhf must be odd')
    END IF

    IF ( (impl_fac<0.5_wp) .OR. (impl_fac>1.0_wp)) THEN
      CALL finish(TRIM(routine),'impl_fac outside permissible range 0.5<=impl_fac<=1.0')
    ENDIF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg=1,max_dom

      wave_config(jg)%ndirs             = ndirs
      wave_config(jg)%nfreqs            = nfreqs
      wave_config(jg)%fr1               = fr1
      wave_config(jg)%CO                = CO
      wave_config(jg)%IREF              = IREF
      wave_config(jg)%ALPHA             = ALPHA
      wave_config(jg)%FM                = FM
      wave_config(jg)%GAMMA_wave        = GAMMA_wave
      wave_config(jg)%SIGMA_A           = SIGMA_A
      wave_config(jg)%SIGMA_B           = SIGMA_B
      wave_config(jg)%FETCH             = FETCH
      wave_config(jg)%roair             = roair
      wave_config(jg)%RNUAIR            = RNUAIR
      wave_config(jg)%RNUAIRM           = RNUAIRM
      wave_config(jg)%ROWATER           = ROWATER
      wave_config(jg)%XEPS              = XEPS
      wave_config(jg)%XINVEPS           = XINVEPS
      wave_config(jg)%XKAPPA            = XKAPPA
      wave_config(jg)%XNLEV             = XNLEV
      wave_config(jg)%BETAMAX           = BETAMAX
      wave_config(jg)%ZALP              = ZALP
      wave_config(jg)%jtot_tauhf        = jtot_tauhf
      wave_config(jg)%ALPHA_CH          = ALPHA_CH
      wave_config(jg)%depth             = depth
      wave_config(jg)%depth_min         = depth_min
      wave_config(jg)%depth_max         = depth_max
      wave_config(jg)%niter_smooth      = niter_smooth
      wave_config(jg)%forc_file_prefix  = forc_file_prefix
      wave_config(jg)%linput_sf1        = linput_sf1
      wave_config(jg)%linput_sf2        = linput_sf2
      wave_config(jg)%ldissip_sf        = ldissip_sf
      wave_config(jg)%lwave_brk_sf      = lwave_brk_sf
      wave_config(jg)%lnon_linear_sf    = lnon_linear_sf
      wave_config(jg)%lbottom_fric_sf   = lbottom_fric_sf
      wave_config(jg)%lwave_stress1     = lwave_stress1
      wave_config(jg)%lwave_stress2     = lwave_stress2
      wave_config(jg)%peak_u10          = peak_u10
      wave_config(jg)%peak_v10          = peak_v10
      wave_config(jg)%peak_lat          = peak_lat
      wave_config(jg)%peak_lon          = peak_lon
      wave_config(jg)%impl_fac          = impl_fac
    ENDDO


    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=wave_nml)
      CALL store_and_close_namelist(funit, 'wave_nml')
    ENDIF

    !------------------------------------------------------------------
    ! 7. Write the namelist to an ASCII file
    !------------------------------------------------------------------
    IF ( my_process_is_stdio() ) WRITE(nnml_output,nml=wave_nml)

  END SUBROUTINE read_wave_namelist

END MODULE mo_wave_nml
