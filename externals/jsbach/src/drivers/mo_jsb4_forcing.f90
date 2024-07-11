!> Module to read and convert external forcing data
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> expects annual files with year in filename
!>
!> current implementation only for daily, timestep or constant data
!>
MODULE mo_jsb4_forcing
#ifndef __NO_JSBACH__

  USE mo_jsb_time,        ONLY: is_time_experiment_start, is_time_restart, get_year_day, &
    &                           get_date_components, t_datetime, is_newyear, is_newday,  &
    &                           get_year_length, get_year_at_experiment_start

  ! Apparently in echam and icon
  USE mo_math_constants,         ONLY: pi
  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util_string,            ONLY: tolower, int2string, real2string
  USE mo_physical_constants,     ONLY: p0sl_bg, vtmpc1, rdv
  USE mo_jsb_physical_constants, ONLY: grav, tmelt, rd, amco2, amd, rvd1, molar_mass_N, molar_mass_P
  USE mo_orbit_solar,            ONLY: inquire_declination

  USE mo_jsb_model_class,        ONLY: t_jsb_model, MODEL_QUINCY, MODEL_JSBACH
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
  USE mo_jsb_parallel,           ONLY: my_process_is_stdio, my_process_is_mpi_parallel, p_io, mpi_comm, p_bcast
  USE mo_jsb_io_netcdf,          ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io_netcdf_iface,    ONLY:           &
    & nf, nf90_inq_varid, nf90_get_var,          &
    & nf90_inq_dimid, nf90_inquire_dimension,    &
    & nf90_get_att, NF90_MAX_NAME, NF90_NOERR

  USE mo_jsb_class,              ONLY: get_model
  USE mo_isotope_util,           ONLY: calc_mixing_ratio_N15N14

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_standalone_driver, setup_forcing, init_forcing, finalize_external_forcing, &
    &       forcing_options, forcing_options_quincy

  REAL(wp), SAVE :: delta_time, time_step_len
  INTEGER,  SAVE :: time_steps_per_day

  LOGICAL, SAVE :: lstart, lresume

  ! We have to define time_days here as there is no equivalent TYPE in ICON-A
  TYPE, PUBLIC :: time_days ! relative calendar date and time format
    ! time_days [structure]
    !   day    [integer]  (day in calendar, -2147483648 ..... 2147483647
    !                      approx. +/-5.8 Mio. years)
    !   second [integer]  (seconds of day, 0,...,86399)
    PRIVATE
    LOGICAL :: init   = .FALSE.
    INTEGER :: day    = 0
    INTEGER :: second = 0
  END TYPE time_days

  ! --- CONSTANTS
  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsbalone_forcing'

!>>>>> TODO: module for general constants / subroutine initialising general constants?
  REAL(wp), PARAMETER :: gamma = 6.0E-3_wp  ! copied from jsbach3 - mo_jsbalone_forcing
                                            ! fixed temperature gradient of atmosphere [K/m] (see Knorr, p. 31). Note that
                                            ! for the standard ICAO-atmosphere this value is different, namely 6.5E-3 !!!
  REAL(wp), PARAMETER :: solar_const =  1361.371_wp      ! default solar constant (AMIP) [W/m2]
  REAL(wp), PARAMETER :: conv_prec_day = 1._wp/86400._wp ! copied from jsbach3 - mo_jsbalone_forcing
                                                         ! Conversion factor for precipitation in [mm/day] to [kg/m^2/s]:
                                                         ! 1 mm/day is equivalent to 1 kg /86400 s / m^2
  REAL(wp), PARAMETER :: molarMassDryAir_kg = amd * 1.e-3_wp  !  28.970e-3_wp   -  Mass of 1 mol of dry air in kg
  REAL(wp), PARAMETER :: molarMassCO2_kg = amco2 * 1.e-3_wp   ! 44.0095e-3_wp   -  Mass of 1 mol CO2 in kg

  REAL(wp), PARAMETER :: defHeightWind = -1._wp        ! -1: for ECHAM lowest layer values
  REAL(wp), PARAMETER :: defHeightHumidity = -1._wp    ! -1: for ECHAM lowest layer values

  !originally from echam mo_physc2.f90 and in jsb3 initialised with jsbalone_iniphy
  REAL(wp), PARAMETER :: cb = 5._wp
  REAL(wp), PARAMETER :: cc = 5._wp
  REAL(wp), PARAMETER :: ckap = 0.4_wp

  CHARACTER(len=*), PARAMETER :: name_of_unit_attribute = 'units'

  ! expected varnames in forcing files
  CHARACTER(len=*), PARAMETER :: max_temp_varname = 'tmax'
  CHARACTER(len=*), PARAMETER :: min_temp_varname = 'tmin'
  CHARACTER(len=*), PARAMETER :: air_temp_varname = 'air_temp'
  CHARACTER(len=*), PARAMETER :: precipitation_varname = 'precip'
  CHARACTER(len=*), PARAMETER :: rel_humidity_varname = 'rel_humidity'
  CHARACTER(len=*), PARAMETER :: qair_varname = 'qair'
  CHARACTER(len=*), PARAMETER :: shortwave_varname = 'shortwave'
  CHARACTER(len=*), PARAMETER :: longwave_varname = 'longwave'
  CHARACTER(len=*), PARAMETER :: CO2_varname = 'CO2'
  CHARACTER(len=*), PARAMETER :: wspeed_varname = 'wspeed'
  ! ... and deposition for QUINCY
  CHARACTER(len=*), PARAMETER :: nhxdep_varname = 'NHx_deposition'
  CHARACTER(len=*), PARAMETER :: noydep_varname = 'NOy_deposition'
  CHARACTER(len=*), PARAMETER :: pdep_varname = 'pdep'

  ! expected units in forcing files (units of precipitation and co2 are specified in the namelist parameters)
  CHARACTER(len=*), PARAMETER :: exp_unit_temperature = 'degC'
  CHARACTER(len=*), PARAMETER :: exp_unit_radiation = 'W/m2'
  CHARACTER(len=*), PARAMETER :: exp_unit_qair = 'g/g'
  CHARACTER(len=*), PARAMETER :: exp_unit_wspeed = 'm/s'
  CHARACTER(len=*), PARAMETER :: unit_kg_m2_s = 'kg/m2/s'
  CHARACTER(len=*), PARAMETER :: exp_unit_deposition = unit_kg_m2_s
  CHARACTER(len=*), PARAMETER :: unit_precip_mm_day = 'mm/day'
  CHARACTER(len=*), PARAMETER :: unspecified_exp_unit = ''

  ! options controlling types of forcing data
  INTEGER,PARAMETER  :: DAILY_    = 1
  INTEGER,PARAMETER  :: SUBDAILY_ = 2
  INTEGER,PARAMETER  :: CONST_    = 3
  INTEGER,PARAMETER  :: TIMESTEP_ = 4
  ! additional option for CO2
  INTEGER,PARAMETER  :: GHG_SCENARIO_  = 5

  ! options controlling atmospheric humidity input data
  INTEGER,PARAMETER :: NONE_ = 0         ! calculate specific humidity from temperature
  INTEGER,PARAMETER :: RH_   = 1         ! relative humidity as forcing data
  INTEGER,PARAMETER :: QAIR_ = 2         ! specific humidity as forcing data


  ! ======================================================================================================= !
  ! TYPEs for handling input data
  ! ======================================================================================================= !


  ! ======================================================================================================= !
  !>
  !> Forcing-related namelist options
  !>
  TYPE forcing_options_type
    INTEGER  :: cyclic_nyears                  !! Number of years for cyclic forcing; default: 0 = no cyclic forcing
    INTEGER  :: cyclic_start_year              !! First year for cyclic forcing; default: 9999, i.e. use experiment start year
    INTEGER  :: forcing_synchron_factor
    INTEGER  :: forcing_steps_per_day          !! The number of time steps included by the used forcing file(s) per day, e.g. 72
    LOGICAL  :: forcing_set_ocean_to_constants !! Specifies if external forcing is only defined for land cells (true)
    LOGICAL  :: forcing_set_miss_to_constants  !! Set missing value in forcing over land cells to constant values (= false)
    LOGICAL  :: air_temperature_as_timestep    !! Specifies if the air temperatures are given per timestep
    INTEGER  :: type_of_qair_forcing           !! Possible values: see 'options controlling atmospheric humidity input data'
    LOGICAL  :: precip_in_mm_per_day           !! Specifies units of precipitation input data: true: mm/day; false: kg/(m^2*s)
    REAL(wp) :: conv_CO2_2_MassRatio           !! Factor for the conversion of CO2 input data to kg(CO2)/kg(dry air)
    REAL(wp) :: heightWind                     !! Reference height for wind speed [m]; -1 for ECHAM lowest layer
    REAL(wp) :: heightHumidity                 !! Reference height for specific humidity [m]; -1 for ECHAM lowest layer
  END TYPE forcing_options_type
  TYPE(forcing_options_type), ALLOCATABLE, SAVE :: forcing_options(:)

  ! ======================================================================================================= !
  !>
  !> this TYPE is NOT used with QUINCY in ICON-Land
  !> implemented here for compatibility with QS (the quincy standalone model)
  !>
  !  QUINCY model
  TYPE t_q_forcing_options
    ! namelist
    ! ..
    ! calculated
    ! ..
    ! read from forcing files
    INTEGER  :: forcing_year = 0
    INTEGER  :: forcing_doy  = 0
    REAL(wp) :: forcing_hour = 0.0_wp
  END TYPE t_q_forcing_options
  TYPE(t_q_forcing_options), SAVE :: forcing_options_quincy

  ! ======================================================================================================= !
  !>
  !> information for single forcing variable and the file it is read from
  !>
  TYPE input_variable_type
    CHARACTER(NF90_MAX_NAME) :: variable_name=""       !! Name of the variable (name in netcdf-data file)
    CHARACTER(NF90_MAX_NAME) :: file_name_prefix = ""  !! Prefix of files in which the variable is located
    CHARACTER(NF90_MAX_NAME) :: expectedUnit = ""      !! Unit for the variable expected in the file
    INTEGER                  :: frequency = 0          !! Input data frequency (see 'options controlling types of forcing data')
    REAL(wp)                 :: ocean_value = 0._wp    !! Value used for ocean cells
    REAL(wp)                 :: constantValue = 0._wp  !! Value in case of frequency 'constant'
    REAL(wp)                 :: missval = 1.0e+20_wp   !! Missing value
    INTEGER                  :: currentDayStart = 0    !! First time step in file corresponding to the current day
    INTEGER                  :: currentDayEnd = 0      !! Last time step in file corresponding to the current day
    INTEGER                  :: length_time_series     !! The number of time steps found in the data (only if my_process_is_stdio())
    REAL(wp), POINTER        :: data_read(:,:,:) => NULL()  !! Pointer to read data
    REAL(wp), POINTER        :: timevalues(:) => NULL()     !! The values of the time variable (only if my_process_is_stdio())
    TYPE(t_input_file)       :: input_file
  END TYPE input_variable_type

  ! ======================================================================================================= !
  !>
  !> All theoretically possible input data, i.e. some are not used. (Input data is controlled by jsb_forcing_nml)
  !>
  TYPE forcing_input_type
     TYPE(input_variable_type) :: air_temp_min  !! daily minimum of air temperature at surface [Celsius]
     TYPE(input_variable_type) :: air_temp_max  !! daily maximum of air temperature at surface [Celsius]
     TYPE(input_variable_type) :: air_temp      !! air temperature [Celsius]
     TYPE(input_variable_type) :: precipitation !! precipitation (including snow); units may be
                                                !! [mm/day] or [kg/m^2/s] (see option  "precip_in_mm_per_day")
     TYPE(input_variable_type) :: shortwave     !! incident shortwave radiation [W/m^2]
     TYPE(input_variable_type) :: longwave      !! incident longwave radiation [W/m^2]
     TYPE(input_variable_type) :: CO2_concentr  !! CO2-concentration [mol(CO2)/mol(air)], [kg(CO2)/kg(dry Air)] or [ppmv]
     TYPE(input_variable_type) :: wind_speed    !! wind speed [m/s]
     TYPE(input_variable_type) :: rel_humidity  !! relative humidity [%]
     TYPE(input_variable_type) :: qair          !! specific humidity [kg/kg]
     ! ... required for Quincy
     TYPE(input_variable_type) :: nhx_dep       !! NHx deposition (dry+wet) [kg/m2/s]
     TYPE(input_variable_type) :: noy_dep       !! NOy deposition (dry+wet) [kg/m2/s]
     TYPE(input_variable_type) :: p_dep         !! P deposition [kg/m2/s]
  END TYPE forcing_input_type
  TYPE(forcing_input_type), ALLOCATABLE, SAVE :: forcing_input(:)

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Reads or updates the external forcing (depending on timestep and temporal resolution of the file)
  !> and derives interface variables that drive the standalone model
  !>
  !> called by 'mo_jsbach_model:jsbach_model:run_one_timestep'
  !>
  SUBROUTINE get_standalone_driver(                                                        &
    & model_id,                                                                            &  ! in
    & model_scheme,                                                                        &
    & startblk, endblk, startidx, endidx,                                                  &
    & current_datetime, next_datetime,                                                     &
    & elevation, sinlat, coslat, lon, evapo_act2pot_proc,                                  &
    & cos_zenith_angle,                                                                    &  ! in
    & CO2_air,                                                                             &  ! out
    & t_air, q_air, rain, snow, wind_air, wind_10m, lw_srf_down, swvis_srf_down,           &
    & swnir_srf_down, swpar_srf_down, fract_par_diffuse, press_srf,                        &
    & nhx_deposition, noy_deposition, nhx_n15_deposition, noy_n15_deposition,              &
    & p_deposition                                                                         &
    & )

    INTEGER,  INTENT(IN)  :: model_id
    INTEGER,  INTENT(IN)  :: model_scheme
    INTEGER,  INTENT(IN)  :: startblk
    INTEGER,  INTENT(IN)  :: endblk
    INTEGER,  INTENT(IN)  :: startidx(:)
    INTEGER,  INTENT(IN)  :: endidx(:)
    TYPE(t_datetime),  INTENT(IN), POINTER :: current_datetime
    TYPE(t_datetime),  INTENT(IN), POINTER :: next_datetime
    REAL(wp), INTENT(IN)  :: elevation(:,:)
    REAL(wp), INTENT(IN)  :: sinlat(:,:)
    REAL(wp), INTENT(IN)  :: coslat(:,:)
    REAL(wp), INTENT(IN)  :: lon(:,:)
    REAL(wp), INTENT(IN)  :: evapo_act2pot_proc(:,:)
    REAL(wp), INTENT(IN)  :: cos_zenith_angle(:,:)
    REAL(wp), INTENT(OUT) :: CO2_air(:,:)
    REAL(wp), INTENT(OUT) :: t_air(:,:)
    REAL(wp), INTENT(OUT) :: q_air(:,:)
    REAL(wp), INTENT(OUT) :: rain(:,:)
    REAL(wp), INTENT(OUT) :: snow(:,:)
    REAL(wp), INTENT(OUT) :: wind_air(:,:)
    REAL(wp), INTENT(OUT) :: wind_10m(:,:)
    REAL(wp), INTENT(OUT) :: lw_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: swvis_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: swnir_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: swpar_srf_down(:,:)
    REAL(wp), INTENT(OUT) :: fract_par_diffuse(:,:)
    REAL(wp), INTENT(OUT) :: press_srf(:,:)
    REAL(wp), INTENT(OUT), OPTIONAL :: nhx_deposition(:,:)
    REAL(wp), INTENT(OUT), OPTIONAL :: noy_deposition(:,:)
    REAL(wp), INTENT(OUT), OPTIONAL :: nhx_n15_deposition(:,:)
    REAL(wp), INTENT(OUT), OPTIONAL :: noy_n15_deposition(:,:)
    REAL(wp), INTENT(OUT), OPTIONAL :: p_deposition(:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)              :: hlp_r
    INTEGER               :: nproma, iblk, nblks, ics, ice, ic, it
    INTEGER               :: day_of_year
    REAL(wp), ALLOCATABLE :: t_air_C(:,:), &
      &                      tmax_data(:,:), &
      &                      tmin_data(:,:), &
      &                      rad_uv_down(:,:), &
      &                      sat_spec_humidity(:,:), &
      &                      rel_humidity(:,:)
    REAL(wp), ALLOCATABLE :: temporary_data1(:,:)
    ! TODO: help1 is rad_sw_down, help2 is rad_sw_down_pot, but these are (currently?) not required - remove?
    REAL(wp), ALLOCATABLE :: help1(:,:)
    REAL(wp), ALLOCATABLE :: help2(:,:)
    REAL(wp), ALLOCATABLE :: tmp(:,:)
    REAL(wp) :: argMin, argMax
    LOGICAL  :: use_quincy

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_standalone_driver'
    ! ----------------------------------------------------------------------------------------------------- !

    lstart  = is_time_experiment_start(current_datetime)
    lresume = is_time_restart(current_datetime)
    use_quincy = (model_scheme == MODEL_QUINCY)

    nproma = SIZE(sinlat,1)
    nblks  = SIZE(sinlat,2)
    ! ----------------------------------------------------------------------------------------------------- !

    ALLOCATE(help1(nproma, nblks))
    ALLOCATE(help2(nproma, nblks))
    ALLOCATE(temporary_data1(nproma, nblks))
    !$ACC ENTER DATA CREATE(help1, help2, temporary_data1)

    ! ----------------------------------------------------------------------------------------------------- !

    ! Initialization needed for grid cells beyond icz
!$OMP PARALLEL DO PRIVATE(iblk, ic)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO iblk = 1, nblks
      DO ic = 1, nproma
        CO2_air(ic, iblk)           = 0._wp
        t_air(ic, iblk)             = 0._wp
        q_air(ic, iblk)             = 0._wp
        rain(ic, iblk)              = 0._wp
        snow(ic, iblk)              = 0._wp
        wind_air(ic, iblk)          = 0._wp
        wind_10m(ic, iblk)          = 0._wp
        lw_srf_down(ic, iblk)       = 0._wp
        swvis_srf_down(ic, iblk)    = 0._wp
        swnir_srf_down(ic, iblk)    = 0._wp
        swpar_srf_down(ic, iblk)    = 0._wp
        fract_par_diffuse(ic, iblk) = 0._wp
        press_srf(ic, iblk)         = 0._wp
        IF (use_quincy) THEN
          nhx_deposition(ic, iblk)     = 0._wp
          noy_deposition(ic, iblk)     = 0._wp
          nhx_n15_deposition(ic, iblk) = 0._wp
          noy_n15_deposition(ic, iblk) = 0._wp
          p_deposition(ic, iblk)       = 0._wp
        END IF
      ENDDO
    ENDDO
    !$ACC END PARALLEL
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> 1.0 Temperature
    !>
    !>   Depending on the forcing frequency different variables are used: Tmin and Tmax need daily forcing,
    !>   while air_temp needs timestep-wise forcing.
    !>

    ALLOCATE(t_air_C(nproma, nblks))
    !$ACC ENTER DATA CREATE(t_air_C)

    IF(forcing_options(model_id)%air_temperature_as_timestep) THEN
#ifdef _OPENACC
      CALL finish(routine, 'Forcing option %air_temperature_as_timestep not ported to GPU, yet. Stop.')
#endif
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%air_temp, t_air_C)
      !$ACC UPDATE DEVICE(t_air_C) ASYNC(1)
    ELSE
      ALLOCATE(tmin_data(nproma, nblks))
      ALLOCATE(tmax_data(nproma, nblks))
      !$ACC ENTER DATA CREATE(tmin_data, tmax_data)

      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%air_temp_min, tmin_data)
      !$ACC UPDATE DEVICE(tmin_data) ASYNC(1)

      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%air_temp_max, tmax_data)
      !$ACC UPDATE DEVICE(tmax_data) ASYNC(1)

        ! estimate temperature at time step from daily temperature maximum and temperature minimum
      CALL instantly_from_daily_temp_2(startblk, endblk, startidx, endidx, &
        & lon, sinlat, coslat, tmin_data, tmax_data, next_datetime, t_air_C)
    ENDIF

    ! convert degC to Kelvin
!$OMP PARALLEL DO PRIVATE(iblk, ic)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO iblk = 1, nblks
      DO ic = 1, nproma
        t_air(ic, iblk) = t_air_C(ic, iblk) + tmelt
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> 2.0 Precipitation
    !>
    CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
      &                            forcing_input(model_id)%precipitation, temporary_data1)
    !$ACC UPDATE DEVICE(temporary_data1) ASYNC(1)

    ! convert if required
    IF(forcing_options(model_id)%precip_in_mm_per_day) THEN
#ifdef _OPENACC
      CALL finish(routine, 'Forcing option %precip_in_mm_per_day not ported to GPU, yet. Stop.')
#endif
!$OMP PARALLEL DO PRIVATE(iblk, ic)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO iblk = 1, nblks
        DO ic = 1, nproma
          temporary_data1(ic, iblk) = temporary_data1(ic, iblk) * conv_prec_day
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
    END IF

    ! compute rain and snow
    IF(forcing_options(model_id)%air_temperature_as_timestep) THEN
#ifdef _OPENACC
      CALL finish(routine, 'Forcing option %air_temperature_as_timestep not ported to GPU, yet. Stop.')
#endif
      IF(forcing_input(model_id)%precipitation%frequency == DAILY_) THEN
        ! In case of daily precipitation the minimum and maximum read for air_temp are used
!$OMP PARALLEL DO PRIVATE(iblk, ic, it, argMin, argMax)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(it, argMin, argMax)
        DO iblk = 1, nblks
          DO ic = 1, nproma
              argMin = forcing_input(model_id)%air_temp%data_read(ic, iblk,1)
              argMax = forcing_input(model_id)%air_temp%data_read(ic, iblk,1)

              DO it = 1, SIZE(forcing_input(model_id)%air_temp%data_read,3)
                  IF(argMin > forcing_input(model_id)%air_temp%data_read(ic, iblk, it)) THEN
                    argMin = forcing_input(model_id)%air_temp%data_read(ic, iblk, it)
                  ENDIF
                  IF(argMax < forcing_input(model_id)%air_temp%data_read(ic, iblk, it)) THEN
                    argMax = forcing_input(model_id)%air_temp%data_read(ic, iblk, it)
                  ENDIF
              ENDDO

              CALL snow_and_rain_from_precip(rain(ic, iblk), snow(ic, iblk), temporary_data1(ic, iblk), &
                &                             0.5_wp*(argMin + argMax))
          ENDDO
        ENDDO
        !$ACC END PARALLEL
!$OMP END PARALLEL DO

      ELSE
        ! else currentvalue
!$OMP PARALLEL DO PRIVATE(iblk, ic)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            CALL snow_and_rain_from_precip(rain(ic, iblk), snow(ic, iblk), temporary_data1(ic, iblk), &
            & t_air_C(ic, iblk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

      ENDIF
    ELSE

!$OMP PARALLEL DO PRIVATE(iblk, ic)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO iblk = 1, nblks
        DO ic = 1, nproma
          CALL snow_and_rain_from_precip(rain(ic, iblk), snow(ic, iblk), temporary_data1(ic, iblk), &
            &                            0.5_wp*(tmax_data(ic, iblk) + tmin_data(ic, iblk)))
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP

    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 3.0 Compute surface pressure
    !>
!$OMP PARALLEL DO PRIVATE(iblk, ics, ice)
    DO iblk = startblk, endblk
      ics = startidx(iblk)
      ice = endidx  (iblk)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = ics, ice
        press_srf(ic,iblk) = bottom_pressure(elevation(ic,iblk), t_air_C(ic,iblk))
      END DO
      !$ACC END PARALLEL LOOP
    ENDDO
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> 4.0 Compute specific humidity
    !>
    ALLOCATE(sat_spec_humidity(nproma, nblks), tmp(nproma, nblks))
    !$ACC ENTER DATA CREATE(sat_spec_humidity, tmp)

!$OMP PARALLEL DO PRIVATE(iblk, ic)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO iblk = 1, nblks
      DO ic = 1, nproma
        sat_spec_humidity(ic, iblk) = 0._wp
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

    SELECT CASE(forcing_options(model_id)%type_of_qair_forcing)
    CASE(NONE_)
#ifdef _OPENACC
      CALL finish(routine, 'Forcing option %type_of_qair_forcing==NONE_ not ported to GPU, yet. Stop.')
#endif
      ! ratio of actual to potential evapotranspiration from the previous day
      IF (forcing_options(model_id)%air_temperature_as_timestep) THEN
        !$ACC KERNELS DEFAULT(PRESENT)
        tmp(:,:) = MINVAL(forcing_input(model_id)%air_temp%data_read(:,:,:), DIM=3)
        !$ACC END KERNELS
        CALL vapour_pressure_from_evapor(startblk, endblk, startidx, endidx, evapo_act2pot_proc(:,:), &
                                         tmp(:,:), t_air_C(:,:), temporary_data1(:,:))
      ELSE
        CALL vapour_pressure_from_evapor(startblk, endblk, startidx, endidx, evapo_act2pot_proc(:,:),  &
          &                              tmin_data(:,:), t_air_C(:,:), temporary_data1(:,:))
      ENDIF
      CALL fun_specific_humidity_2d(temporary_data1(:,:), press_srf(:,:), q_air(:,:))

    CASE(RH_)
#ifdef _OPENACC
      CALL finish(routine, 'Forcing option %type_of_qair_forcing==RH_ not ported to GPU, yet. Stop.')
#endif
      ALLOCATE(rel_humidity(nproma, nblks))
      !$ACC ENTER DATA CREATE(rel_humidity)
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        &                            forcing_input(model_id)%rel_humidity, rel_humidity)
      !$ACC UPDATE DEVICE(rel_humidity) ASYNC(1)
      ! check data
      hlp_r = MAXVAL(MAXVAL(rel_humidity,DIM=2))
      IF(hlp_r < 1._wp) THEN
        WRITE (message_text,*) 'Maximum value of relative humidity is ', hlp_r
        CALL message(TRIM(routine), TRIM(message_text))
        CALL message(TRIM(routine),'relative humidity data probably represent fractions instead of percentages.')
      ENDIF
      ! Convert percentage data to fractional data:
!$OMP PARALLEL DO PRIVATE(ic, iblk)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO iblk=1,nblks
        DO ic=1,nproma
          rel_humidity(ic,iblk) = rel_humidity(ic,iblk)/100._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
      ! Get saturation spec. humidity from air temperature
      DO iblk = startblk, endblk
        ics = startidx(iblk)
        ice = endidx  (iblk)
        CALL sat_specific_humidity(t_air(ics:ice,iblk), press_srf(ics:ice,iblk), sat_spec_humidity(ics:ice,iblk))
      END DO
      ! Get specific humidity
      IF (     forcing_input(model_id)%rel_humidity%frequency == TIMESTEP_ &
        & .OR. forcing_input(model_id)%rel_humidity%frequency == SUBDAILY_) THEN
!$OMP PARALLEL DO PRIVATE(ic, iblk)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
          DO iblk=1,nblks
            DO ic=1,nproma
              q_air(ic,iblk) = rel_humidity(ic,iblk) * sat_spec_humidity(ic,iblk) !TODO: Why *?
            END DO
          END DO
          !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
      ELSE
        IF (forcing_options(model_id)%air_temperature_as_timestep) THEN
          DO iblk = startblk, endblk
            ics = startidx(iblk)
            ice = endidx  (iblk)
            CALL spec_humidity_from_rel_humidity(ice-ics+1, sinlat(ics:ice,iblk), coslat(ics:ice,iblk), &
              & rel_humidity(ics:ice,iblk), MINVAL(forcing_input(model_id)%air_temp%data_read(ics:ice,iblk,:),DIM=2), &
              & MAXVAL(forcing_input(model_id)%air_temp%data_read(ics:ice,iblk,:),DIM=2), &
              & elevation(ics:ice,iblk), q_air(ics:ice,iblk))
          ENDDO
        ELSE
          DO iblk = startblk, endblk
            ics = startidx(iblk)
            ice = endidx  (iblk)
            CALL spec_humidity_from_rel_humidity( &
              & ice-ics+1, sinlat(ics:ice,iblk), coslat(ics:ice,iblk), rel_humidity(ics:ice,iblk), &
              & tmin_data(ics:ice,iblk), tmax_data(ics:ice,iblk), elevation(ics:ice,iblk), q_air(ics:ice,iblk))
          ENDDO
        ENDIF
        ! constrain daily specific humidity to be lower than or equal to the saturation spec. humidity
!$OMP PARALLEL DO PRIVATE(ic, iblk)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk=1,nblks
          DO ic=1,nproma
            q_air(ic,iblk) = MIN(q_air(ic,iblk), sat_spec_humidity(ic,iblk))
          END DO
        END DO
        !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
      ENDIF
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(rel_humidity)
      DEALLOCATE(rel_humidity)

    ! @TODO add docu

    CASE(QAIR_)
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%qair, q_air(:,:))
      !$ACC UPDATE DEVICE(q_air) ASYNC(1)

      ! constrain specific humidity to be lower than or equal to the saturation spec. humidity
      DO iblk = startblk, endblk
        ics = startidx(iblk)
        ice = endidx  (iblk)
        CALL sat_specific_humidity(t_air(ics:ice,iblk), press_srf(ics:ice,iblk), sat_spec_humidity(ics:ice,iblk))
      ENDDO
!$OMP PARALLEL DO PRIVATE(ic, iblk)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO iblk=1,nblks
        DO ic=1,nproma
          q_air(ic,iblk) = MIN(q_air(ic,iblk),sat_spec_humidity(ic,iblk))
        END DO
      END DO
      !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
    END SELECT

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(sat_spec_humidity, tmp)
    DEALLOCATE(sat_spec_humidity, tmp)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 5.0 Compute shortwave radiation
    !>
    day_of_year=day_in_year(next_datetime)

    ALLOCATE(rad_uv_down(nproma, nblks))
    !$ACC ENTER DATA CREATE(rad_uv_down)

    CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
      & forcing_input(model_id)%shortwave, temporary_data1)
    !$ACC UPDATE DEVICE(temporary_data1) ASYNC(1)

!$OMP PARALLEL
    CALL shortwave_from_direct_shortwave( &
      ! in
      & model_id,               &
      & day_of_year,            &
      & cos_zenith_angle (:,:), &
      & coslat           (:,:), &
      & sinlat           (:,:), &
      & press_srf        (:,:), &
      & temporary_data1  (:,:), &
      ! out
      & rad_uv_down      (:,:), & ! rad_UV
      & swpar_srf_down   (:,:), & ! rad_PAR
      & swnir_srf_down   (:,:), & ! rad_NIR
      & fract_par_diffuse(:,:), & ! fract_PAR_diffuse
      & help1            (:,:), & ! rad_sw
      & help2            (:,:))   ! rad_sw_pot
!$OMP END PARALLEL

!$OMP PARALLEL DO PRIVATE(ic, iblk)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO iblk=1,nblks
      DO ic=1,nproma
        swvis_srf_down(ic,iblk) = rad_uv_down(ic,iblk) + swpar_srf_down(ic,iblk)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(rad_uv_down)
    DEALLOCATE(rad_uv_down)

    ! ----------------------------------------------------------------------------------------------------- !
    !> 6.0 Compute longwave radiation
    !>
    IF (forcing_input(model_id)%longwave%frequency == TIMESTEP_ .OR. forcing_input(model_id)%longwave%frequency == SUBDAILY_) THEN
#ifdef _OPENACC
      CALL finish(routine, 'Forcing input %longwave%frequency==TIMESTEP_ or SUBDAILY_  not ported to GPU, yet. Stop.')
#endif
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%longwave, lw_srf_down)
      !$ACC UPDATE DEVICE(lw_srf_down) ASYNC(1)
    ELSE
      ! jsbach3 driver: diurnal longwave downward radiation is generated such that the average radiation matches the
      ! observations but the diurnal cycle follows temperature to help calculation of a correct
      ! net radiation flux. Note that reduces biases resulting from assuming a constant flux, but
      ! this introduces uncertainty due to biases in the assumed dial course of temperature
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        &                            forcing_input(model_id)%longwave, temporary_data1)
      !$ACC UPDATE DEVICE(temporary_data1) ASYNC(1)

      IF (forcing_options(model_id)%air_temperature_as_timestep) THEN
#ifdef _OPENACC
        CALL finish(routine, 'Forcing option %air_temperature_as_timestep for longwave radiation not ported to GPU, yet. Stop.')
#endif
!$OMP PARALLEL
        !$ACC KERNELS DEFAULT(PRESENT)
        lw_srf_down(:,:) = longwave_from_daily_longwave(temporary_data1(:,:), &
          & MINVAL(forcing_input(model_id)%air_temp%data_read(:,:,:),DIM=3), &
          & MAXVAL(forcing_input(model_id)%air_temp%data_read(:,:,:),DIM=3), t_air_C(:,:))
        !$ACC END KERNELS
!$OMP END PARALLEL
      ELSE
!$OMP PARALLEL DO PRIVATE(ic, iblk)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO iblk = 1, nblks
          DO ic = 1, nproma
            lw_srf_down(ic, iblk) = longwave_from_daily_longwave(temporary_data1(ic, iblk), tmin_data(ic, iblk), &
              &                                                  tmax_data(ic, iblk), t_air_C(ic, iblk))
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO
      ENDIF
    ENDIF

    ! ----------------------------------------------------------------------------------------------------- !
    !> 7.0 Wind
    !>
    CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
      & forcing_input(model_id)%wind_speed, wind_10m)
    !$ACC UPDATE DEVICE(wind_10m) ASYNC(1)

!$OMP PARALLEL DO PRIVATE(ic, iblk)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO iblk = 1, nblks
      DO ic = 1, nproma
        wind_air(ic, iblk) = wind_10m(ic, iblk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> 8.0 CO2 concentration
    !>
    CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
      & forcing_input(model_id)%CO2_concentr, CO2_air)
    !$ACC UPDATE DEVICE(CO2_air) ASYNC(1)

    ! and convert to correct unit
!$OMP PARALLEL DO PRIVATE(ic, iblk)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) COPYIN(forcing_options(model_id)) ASYNC(1)
    DO iblk = 1, nblks
      DO ic = 1, nproma
        CO2_air(ic, iblk) = CO2_air(ic, iblk) * forcing_options(model_id)%conv_CO2_2_MassRatio
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------------------------------------- !
    !> 9.0 N & P deposition
    !>
    SELECT CASE (model_scheme)
    CASE (MODEL_QUINCY)
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%nhx_dep, nhx_deposition)
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%noy_dep, noy_deposition)
      CALL get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
        & forcing_input(model_id)%p_dep, p_deposition)
      !$ACC UPDATE DEVICE(nhx_deposition, noy_deposition, p_deposition)

      ! Convert from kg/m2/s -> mumol/m2/s
      ! -- This unit should have been tested above if specified correctly using expectedUnit
      IF(.NOT. (forcing_input(model_id)%nhx_dep%expectedUnit == unit_kg_m2_s)         &
          & .OR. .NOT. (forcing_input(model_id)%noy_dep%expectedUnit == unit_kg_m2_s) &
          & .OR. .NOT. (forcing_input(model_id)%p_dep%expectedUnit == unit_kg_m2_s)) THEN
        WRITE (message_text,*) 'Deposition should be read using the unit "' //TRIM(unit_kg_m2_s)                    &
          & //'" - a specified expected unit does not match ("'//TRIM(forcing_input(model_id)%nhx_dep%expectedUnit) &
          & //'","'//TRIM(forcing_input(model_id)%noy_dep%expectedUnit)                                             &
          & //'","'//TRIM(forcing_input(model_id)%p_dep%expectedUnit)//'").'
        CALL finish(TRIM(routine), message_text)
      ENDIF

!$OMP PARALLEL DO PRIVATE(ic, iblk)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) COPYIN(forcing_options(model_id)) ASYNC(1)
      DO iblk = 1, nblks
        DO ic = 1, nproma
          ! -- 1 kg/m2/s = 1 * 1000 / molar_mass_N / 1e-06 [mu mol/m2/s]
          nhx_deposition(ic, iblk) = nhx_deposition(ic, iblk) * 1000._wp / molar_mass_N / 1.e-6_wp
          noy_deposition(ic, iblk) = noy_deposition(ic, iblk) * 1000._wp / molar_mass_N / 1.e-6_wp
          p_deposition  (ic, iblk) = p_deposition  (ic, iblk) * 1000._wp / molar_mass_P / 1.e-6_wp

          ! SZ: NHX_deposition_N15 = NOY_deposition_N15 = NOY_deposition f(delta_n15=0.0)
          nhx_n15_deposition(ic, iblk) = nhx_deposition(ic, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          noy_n15_deposition(ic, iblk) = noy_deposition(ic, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        ENDDO
      ENDDO
      !$ACC END PARALLEL LOOP
!$OMP END PARALLEL DO

    END SELECT

    ! ----------------------------------------------------------------------------------------------------- !
    !> 10.0 deallocate local variables
    !>
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(t_air_C)
    DEALLOCATE(t_air_C)

    IF(.NOT. forcing_options(model_id)%air_temperature_as_timestep) THEN
      !$ACC EXIT DATA DELETE(tmin_data, tmax_data)
      DEALLOCATE(tmin_data)
      DEALLOCATE(tmax_data)
    ENDIF

    !$ACC EXIT DATA DELETE(temporary_data1, help1, help2)
    DEALLOCATE(temporary_data1)
    DEALLOCATE(help1, help2)

  END SUBROUTINE get_standalone_driver

  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  SUBROUTINE setup_forcing(no_of_models)

    INTEGER, INTENT(in) :: no_of_models

    CHARACTER(len=*), PARAMETER :: routine = modname//':setup_forcing'

    ALLOCATE(forcing_input(no_of_models))
    ALLOCATE(forcing_options(no_of_models))

    CALL message(routine, 'Allocated forcing_input and forcing_options ('//TRIM(int2string(no_of_models))//')')

  END SUBROUTINE setup_forcing

  ! ======================================================================================================= !
  !>
  !> Initialize and configure reading the forcing data (read namelist, set options...)
  !>
  SUBROUTINE init_forcing(model_id, namelist_filename, dtime, steplen)

    USE mo_jsb_namelist_iface, ONLY: POSITIONED, open_nml, position_nml
    USE mo_jsb_io_netcdf_iface, ONLY: t_input_file, netcdf_open_input

    INTEGER,  INTENT(in) :: model_id
    CHARACTER(len=*), INTENT(in) :: namelist_filename
    REAL(wp), INTENT(in) :: dtime, steplen

    TYPE(t_input_file) :: input_file

    ! --- Namelist parameters
    INTEGER  :: cyclic_nyears           ! Number of years to use for cyclic forcing; default: 0
    INTEGER  :: cyclic_start_year       ! Start year of cyclic forcing; default: 9999
    INTEGER  :: forcing_synchron_factor ! Factor defines how many time steps the JS4 driver performs for
                                        ! each atmospheric forcing time step that is read in
    INTEGER  :: forcing_steps_per_day   ! The number of time steps included by the used forcing file(s) per day, e.g. 72
    LOGICAL                :: forcing_set_ocean_to_constants   !! Some external forcing is only defined for land cells (= true)
    !! false: forcing is expected for all variables and cells (e.g. echam)
    LOGICAL                :: forcing_set_miss_to_constants    !! Set missing value in forcing over land cells to constant values (= false)
    REAL(wp)               :: forcing_height_wind              !! jsb3: Defines lowest layer height-> where measurements are taken // -1: for ECHAM lowest layer values
    REAL(wp)               :: forcing_height_humidity          !! jsb3: Defines lowest layer height-> where measurements are taken // -1: for ECHAM lowest layer values

    !! ++ as in jsb3: "HeightHumidity and HeightTemperature have to be equal"
    CHARACTER(NF90_MAX_NAME) :: forcing_temp_file_prefix       !! Prefix for files with temperature forcing data
    CHARACTER(len=128)     :: forcing_temp_frequ               !! Frequency of temperature forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_temp_const_tmin          !! Constant value for minimum daily temperature
    REAL(wp)               :: forcing_temp_const_tmax          !! Constant value for maximum daily temperature
    REAL(wp)               :: forcing_temp_ocean               !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_temp_unit                !! Expected unit of temp input

    CHARACTER(NF90_MAX_NAME) :: forcing_precip_file_prefix     !! Prefix for files with precipitation forcing data (unit: mm/day
    !! or kg/m^2/s, depending on 'forcing_precip_in_mm_per_day')
    CHARACTER(len=128)     :: forcing_precip_frequ             !! Frequency of precipitation forcing (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_precip_const_precip      !! Constant value for precipitation
    REAL(wp)               :: forcing_precip_ocean             !! Constant value for ocean cells
    LOGICAL                :: forcing_precip_in_mm_per_day     !! Precipitation input in mm/day (true) or kg/m2/s (false)
    CHARACTER(len=128)     :: forcing_precip_unit              !! Expected unit of precip input

    CHARACTER(NF90_MAX_NAME) :: forcing_sw_file_prefix         !! Prefix for files with shortwave radiation data
    CHARACTER(len=128)     :: forcing_sw_frequ                 !! Frequency of shortwave forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_sw_const_shortwave       !! Constant value for downward shortwave radiation
    REAL(wp)               :: forcing_sw_ocean                 !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_sw_unit                  !! Expected unit of shortwave input (default: [W/m^2])

    CHARACTER(NF90_MAX_NAME) :: forcing_lw_file_prefix         !! Prefix for files with longwave radiation data
    CHARACTER(len=128)     :: forcing_lw_frequ                 !! Frequency of longwave forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_lw_const_longwave        !! Constant value for downward longwave radiation
    REAL(wp)               :: forcing_lw_ocean                 !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_lw_unit                  !! Expected unit of longwave input (default: [W/m^2])

    CHARACTER(NF90_MAX_NAME) :: forcing_co2_file_prefix        !! Prefix for files with CO2 data (unit: mol(CO2)/mol(air)
    !! or kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit'
    CHARACTER(len=128)     :: forcing_co2_frequ                !! Frequency of CO2 forcing data (DAILY/CONST/TIMESTEP/GHG_SCENARIO)
    REAL(wp)               :: forcing_co2_const_co2            !! Constant value for CO2-concentration (unit: mol(CO2)/mol(air) or
    !! kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit')
    REAL(wp)               :: forcing_co2_ocean                !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_co2_unit                 !! Unit of CO2 input (PPMV/MOL_PER_MOL/KG_PER_KG/'')

    CHARACTER(NF90_MAX_NAME) :: forcing_wind_file_prefix       !! Prefix for files with windspeed data
    CHARACTER(len=128)     :: forcing_wind_frequ               !! Frequency of windspeed forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_wind_const_wspeed        !! Constant value for windspeed [m/s]
    REAL(wp)               :: forcing_wind_ocean               !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_wind_unit                !! Expected unit of wind input (default: [m/s])

    CHARACTER(NF90_MAX_NAME) :: forcing_qair_file_prefix       !! Prefix for files with qair data
    CHARACTER(len=128)     :: forcing_qair_frequ               !! Frequency of qair forcing data (DAILY/CONST)
    CHARACTER(len=128)     :: forcing_qair_type                !! Type of input data used for forcing by atm. humidity
    REAL(wp)               :: forcing_qair_const_rh            !! Constant value for relative humidity
    REAL(wp)               :: forcing_qair_ocean               !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_qair_unit                !! Expected unit of qair input

    CHARACTER(NF90_MAX_NAME) :: forcing_dep_file_prefix        !! Prefix for files with NHx, NOy and P deposition data
    CHARACTER(len=128)     :: forcing_dep_frequ                !! Frequency of deposition data (DAILY/CONST)
    REAL(wp)               :: forcing_nhx_dep_const            !! Constant value for NHx deposition
    REAL(wp)               :: forcing_noy_dep_const            !! Constant value for NOy deposition
    REAL(wp)               :: forcing_p_dep_const              !! Constant value for P deposition
    CHARACTER(len=128)     :: forcing_dep_unit                 !! Expected unit of deposition input

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_forcing'

    NAMELIST /jsb_forcing_nml/ cyclic_nyears, cyclic_start_year, forcing_synchron_factor, forcing_steps_per_day,              &
      & forcing_set_ocean_to_constants, forcing_set_miss_to_constants, forcing_height_wind, forcing_height_humidity,          &
      & forcing_temp_file_prefix, forcing_temp_frequ, forcing_temp_const_tmin, forcing_temp_const_tmax, forcing_temp_ocean,   &
      & forcing_temp_unit,                                                                                                    &
      & forcing_precip_file_prefix, forcing_precip_frequ, forcing_precip_const_precip,                                        &
      & forcing_precip_in_mm_per_day, forcing_precip_ocean, forcing_precip_unit,                                              &
      & forcing_sw_file_prefix, forcing_sw_frequ, forcing_sw_ocean, forcing_sw_unit,                                          &
      & forcing_lw_file_prefix, forcing_lw_frequ, forcing_lw_ocean, forcing_lw_unit,                                          &
      & forcing_co2_file_prefix, forcing_co2_frequ, forcing_co2_const_co2, forcing_co2_ocean, forcing_co2_unit,               &
      & forcing_wind_file_prefix, forcing_wind_frequ, forcing_wind_const_wspeed, forcing_wind_ocean, forcing_wind_unit,       &
      & forcing_qair_file_prefix, forcing_qair_frequ, forcing_qair_type, forcing_qair_const_rh, forcing_qair_ocean,           &
      & forcing_qair_unit, forcing_dep_file_prefix, forcing_dep_frequ, forcing_nhx_dep_const,                                 &
      & forcing_noy_dep_const, forcing_p_dep_const, forcing_dep_unit

    TYPE(t_jsb_model), POINTER :: model
    INTEGER:: nml_handler, read_status, f_unit

    !------------------------------------------------------------------------------------
    model => Get_model(model_id)
    delta_time    = dtime
    time_step_len = steplen
    time_steps_per_day = 86400/INT(delta_time)
    ! Assertion: delta_time needs to be smaller or equal 86400
    IF ((delta_time > 86400._wp)) CALL finish(TRIM(routine), &
      & 'Violation of assertion: Reading forcing data not implemented for delta_time > day.')

    IF (my_process_is_stdio()) THEN
      ! ----- define default values for jsb_forcing_nml namelist parameters
      cyclic_nyears = 0             ! No cyclic forcing
      cyclic_start_year = 9999      ! Use year of experiment start for first cyclic year
      forcing_synchron_factor  = 1  ! Factor defines how many time steps the JS4 driver performs for
                                    ! each atmospheric forcing time step that is read in
      forcing_steps_per_day    = 1  ! Number of forcing time steps per day, e.g. 1 for daily forcing

      forcing_set_ocean_to_constants = .TRUE.
      forcing_set_miss_to_constants = .FALSE.
      forcing_height_wind = defHeightWind
      forcing_height_humidity = defHeightHumidity

      forcing_temp_file_prefix = 'climate_'
      forcing_temp_frequ = 'DAILY'
      forcing_temp_const_tmin = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_temp_const_tmax = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_temp_ocean = 0._wp
      forcing_temp_unit = exp_unit_temperature

      forcing_precip_file_prefix = 'climate_'
      forcing_precip_frequ = 'DAILY'
      forcing_precip_const_precip = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_precip_in_mm_per_day = .TRUE.
      forcing_precip_ocean = 0._wp
      forcing_precip_unit = unit_precip_mm_day

      forcing_sw_file_prefix = 'climate_'
      forcing_sw_frequ = 'DAILY'
      forcing_sw_ocean = 0._wp
      forcing_sw_unit = exp_unit_radiation

      forcing_lw_file_prefix = 'climate_'
      forcing_lw_frequ = 'DAILY'
      forcing_lw_ocean = 100._wp
      forcing_lw_unit = exp_unit_radiation

      forcing_co2_file_prefix = 'climate_'
      forcing_co2_frequ = 'DAILY'
      forcing_co2_const_co2 = 3.67e-4_wp
      forcing_co2_unit = unspecified_exp_unit
      forcing_co2_ocean = 3.67e-4_wp

      forcing_wind_file_prefix = 'climate_'
      forcing_wind_frequ = 'DAILY'
      forcing_wind_const_wspeed = 0.0_wp !jsb3: HUGE(0.0_wp)
      forcing_wind_ocean = 0._wp
      forcing_wind_unit = exp_unit_wspeed

      forcing_qair_file_prefix = 'climate_'
      forcing_qair_frequ = 'DAILY'
      forcing_qair_type = 'QAIR'
      forcing_qair_const_rh = 100._wp
      forcing_qair_ocean = 1.0e-05
      forcing_qair_unit = exp_unit_qair

      ! N and P deposition required for QUINCY
      forcing_dep_file_prefix = 'deposition_'
      forcing_dep_frequ       = 'DAILY' ! 'CONST'
      ! SZ: "background N deposition is 2 kg N / ha / yr = 6.3419e-12 kg / m2 / s (2/365/24/3600/10000)
      !      the background N deposition is the sum of NHx and NOy, therefore "/ 2._wp"
      forcing_nhx_dep_const   = 6.3419e-12_wp / 2._wp
      forcing_noy_dep_const   = 6.3419e-12_wp / 2._wp
      ! SZ: "Default P deposition is assumed to be stoichmetrically balanced (which would be higher than typical background values.)
      !      Thus, total N deposition to total P deposition = 14 g N/ g P -> 14 / molar_mass_N * molar_mass_P"
      !      (Sterner & Elser, Ecological Stoichmetry, 2002, Princeton University Press)
      forcing_p_dep_const     = 6.3419e-12 / (14.0_wp / molar_mass_N * molar_mass_P)
      forcing_dep_unit        = exp_unit_deposition

      ! ----- read the namelist

      ! nml_handler = open_nml ('NAMELIST_JS4-standalone_master')
      nml_handler = open_nml (TRIM(namelist_filename))
      f_unit = position_nml ('jsb_forcing_nml', nml_handler, status=read_status)
      IF (read_status .EQ. POSITIONED) READ(f_unit, jsb_forcing_nml)

      forcing_options(model_id)%cyclic_nyears           = cyclic_nyears
      IF (cyclic_start_year == 9999) cyclic_start_year = get_year_at_experiment_start()
      forcing_options(model_id)%cyclic_start_year       = cyclic_start_year
      forcing_options(model_id)%forcing_synchron_factor = forcing_synchron_factor
      forcing_options(model_id)%forcing_steps_per_day   = forcing_steps_per_day

      IF ((forcing_steps_per_day > time_steps_per_day)) CALL finish(TRIM(routine), &
        & 'Violation of assertion: Forcing time resolution higher than model time step.')
      IF ((forcing_steps_per_day < 1)) CALL finish(TRIM(routine), &
        & 'Violation of assertion: Forcing time step longer than 1 day is not yet implemented.')

      IF (cyclic_nyears > 0) THEN
        CALL message(routine, 'Using cyclic forcing from ' // TRIM(int2string(cyclic_nyears)) // ' years ' // &
          &                   'starting in ' // TRIM(int2string(cyclic_start_year)))
      END IF

        ! ----- derive options - forcing heights: as in jsb3 it is expected that "HeightHumidity and HeightTemperature have to be equal"
      forcing_options(model_id)%heightWind = forcing_height_wind
      forcing_options(model_id)%heightHumidity = forcing_height_humidity
      ! --- Temperature timestep (if daily two variables - max and min - need to be read)
      forcing_options(model_id)%air_temperature_as_timestep &
        & = (frequency_key_to_constant_val('temperature', forcing_temp_frequ) == TIMESTEP_ .OR. &
        &    frequency_key_to_constant_val('temperature', forcing_temp_frequ) == SUBDAILY_)
      ! --- Boolean indicating if forcing is only specified for land cells
      forcing_options(model_id)%forcing_set_ocean_to_constants = forcing_set_ocean_to_constants
      IF(forcing_options(model_id)%forcing_set_ocean_to_constants) THEN
        CALL message(TRIM(routine),'forcing_set_ocean_to_constants = true: ocean cells will be set to the specified constants. '&
          & //'On ERROR: make sure that your jsbach4 land sea mask matches the land sea mask of your forcing!')
      ELSE
        CALL message(TRIM(routine),'forcing_set_ocean_to_constants = false: input is expected for all cells for all variables. '&
          & //'On ERROR: make sure that your forcing is actually defined for all cells for all variables!')
      ENDIF
      forcing_options(model_id)%forcing_set_miss_to_constants = forcing_set_miss_to_constants
      IF(forcing_options(model_id)%forcing_set_miss_to_constants) THEN
        CALL message(TRIM(routine),&
          & 'forcing_set_miss_to_constants == true: Missing values over land will be replace with ocean constants')
      ENDIF
      ! --- Precipitation unit
      forcing_options(model_id)%precip_in_mm_per_day = forcing_precip_in_mm_per_day
      IF(forcing_options(model_id)%precip_in_mm_per_day) THEN
        CALL message(TRIM(routine),'Precipitation forcing in mm/day')
      ELSE
        CALL message(TRIM(routine),'Precipitation forcing in kg/m2/s')
      ENDIF
      ! --- QAIR type
      SELECT CASE(TRIM(tolower(forcing_qair_type)))
      CASE("rh")
        forcing_options(model_id)%type_of_qair_forcing = RH_
      CASE("qair")
        forcing_options(model_id)%type_of_qair_forcing = QAIR_
      CASE("none")
        forcing_options(model_id)%type_of_qair_forcing = NONE_
        CALL message(TRIM(routine),'WARNING: QAIR will be diagnosed from Air temperature (type_of_qair_forcing = none)')
      CASE default
        CALL finish(TRIM(routine),'Please specifiy either rh, qair or none for "FORCING_QAIR_TYPE" in jsb_forcing_nml')
      END SELECT
      CALL message(TRIM(routine),"Type of air humidity forcing: "//TRIM(forcing_qair_type))
      ! --- CO2 unit
      SELECT CASE(TRIM(tolower(forcing_co2_unit)))
      CASE("ppmv","1.e-6")
        CALL message(TRIM(routine),'CO2 forcing in ppmv')
        forcing_options(model_id)%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg * 0.000001_wp
      CASE("mol_per_mol", "mol/mol")
        CALL message(TRIM(routine),'CO2 forcing in mol(CO2)/mol(Dry Air)')
        forcing_options(model_id)%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg
      CASE("kg_per_kg", "kg/kg")
        CALL message(TRIM(routine),'CO2 forcing in kg(CO2)/kg(Dry Air)')
        forcing_options(model_id)%conv_CO2_2_MassRatio = 1.0_wp
      CASE default
        CALL message(TRIM(routine),'WARNING: "FORCING_CO2_UNIT" missing in namelist jsb_forcing_nml.')
        CALL message(TRIM(routine),'Assumption: unit is mol(CO2)/mol(DryAir). Verify your settings!')
        forcing_options(model_id)%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg
      END SELECT
      ! --- get dimensions of CO2 data
      IF (frequency_key_to_constant_val(CO2_varname, forcing_co2_frequ) == GHG_SCENARIO_) THEN
        input_file =  netcdf_open_input(TRIM(forcing_co2_file_prefix))
        IF (input_file%get_dimlen('lon') == 1 .AND. input_file%get_dimlen('lat') == 1) THEN
          CALL message(TRIM(routine),'CO2 forcing set to a globally constant value')
        ELSE
          CALL finish(TRIM(routine),'Dimensions in GHG file differ from expectations')
        ENDIF
        CALL input_file%Close()
      ENDIF
    ENDIF

    ! ----- broadcast information
    IF (my_process_is_mpi_parallel()) THEN
      ! forcing options
      CALL p_bcast(forcing_options(model_id)%cyclic_nyears, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%cyclic_start_year, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%forcing_synchron_factor, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%forcing_steps_per_day, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%forcing_set_ocean_to_constants, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%forcing_set_miss_to_constants, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%air_temperature_as_timestep, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%precip_in_mm_per_day, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%type_of_qair_forcing, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%conv_CO2_2_MassRatio, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%heightWind, p_io, mpi_comm)
      CALL p_bcast(forcing_options(model_id)%heightHumidity, p_io, mpi_comm)

      ! other namelist parameters
      CALL p_bcast(forcing_temp_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_temp_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_temp_const_tmin, p_io, mpi_comm)
      CALL p_bcast(forcing_temp_const_tmax, p_io, mpi_comm)
      CALL p_bcast(forcing_temp_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_temp_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_precip_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_precip_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_precip_const_precip, p_io, mpi_comm)
      CALL p_bcast(forcing_precip_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_precip_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_sw_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_sw_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_sw_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_sw_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_lw_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_lw_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_lw_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_lw_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_co2_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_co2_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_co2_const_co2, p_io, mpi_comm)
      CALL p_bcast(forcing_co2_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_co2_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_wind_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_wind_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_wind_const_wspeed, p_io, mpi_comm)
      CALL p_bcast(forcing_wind_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_wind_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_qair_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_qair_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_qair_const_rh, p_io, mpi_comm)
      CALL p_bcast(forcing_qair_ocean, p_io, mpi_comm)
      CALL p_bcast(forcing_qair_unit, p_io, mpi_comm)

      CALL p_bcast(forcing_dep_file_prefix, p_io, mpi_comm)
      CALL p_bcast(forcing_dep_frequ, p_io, mpi_comm)
      CALL p_bcast(forcing_nhx_dep_const, p_io, mpi_comm)
      CALL p_bcast(forcing_noy_dep_const, p_io, mpi_comm)
      CALL p_bcast(forcing_p_dep_const, p_io, mpi_comm)
      CALL p_bcast(forcing_dep_unit, p_io, mpi_comm)
    ENDIF

    ! ----- init forcing_input
    ! --- temperature
    IF(forcing_options(model_id)%air_temperature_as_timestep) THEN
      CALL init_forcing_input_variable(air_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, forcing_temp_unit, forcing_temp_const_tmax, forcing_input(model_id)%air_temp)
    ELSE
      CALL init_forcing_input_variable(min_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, forcing_temp_unit, forcing_temp_const_tmin, forcing_input(model_id)%air_temp_min)
      CALL init_forcing_input_variable(max_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, forcing_temp_unit, forcing_temp_const_tmax, forcing_input(model_id)%air_temp_max)
    ENDIF
    ! --- precipitation
    CALL init_forcing_input_variable(precipitation_varname, forcing_precip_frequ, forcing_precip_ocean, &
      & forcing_precip_file_prefix, forcing_precip_unit, forcing_precip_const_precip, forcing_input(model_id)%precipitation)
    ! --- atmospheric humidity
    SELECT CASE(forcing_options(model_id)%type_of_qair_forcing)
    CASE(RH_)
      CALL init_forcing_input_variable(rel_humidity_varname, forcing_qair_frequ, forcing_qair_ocean, &
        forcing_qair_file_prefix, forcing_qair_unit ,forcing_qair_const_rh, forcing_input(model_id)%rel_humidity)
    CASE(QAIR_)
      CALL init_forcing_input_variable(qair_varname, forcing_qair_frequ, forcing_qair_ocean, &
        forcing_qair_file_prefix, forcing_qair_unit ,forcing_qair_const_rh, forcing_input(model_id)%qair)
    CASE(NONE_)
      ! no extra forcing needs to be initialised
    CASE default
      CALL finish(TRIM(routine),'Please specifiy either qair, rh or none for "forcing_qair_type" in jsb_forcing_nml')
    END SELECT
    ! --- shortwave radiation
    CALL init_forcing_input_variable(shortwave_varname, forcing_sw_frequ, forcing_sw_ocean, &
      & forcing_sw_file_prefix, forcing_sw_unit, forcing_sw_const_shortwave, forcing_input(model_id)%shortwave)
    ! --- longwave radiation
    CALL init_forcing_input_variable(longwave_varname, forcing_lw_frequ, forcing_lw_ocean, &
      & forcing_lw_file_prefix, forcing_lw_unit, forcing_lw_const_longwave, forcing_input(model_id)%longwave)
    ! --- CO2
    CALL init_forcing_input_variable(CO2_varname, forcing_co2_frequ, forcing_co2_ocean, &
      & forcing_co2_file_prefix, forcing_co2_unit, forcing_co2_const_co2, forcing_input(model_id)%CO2_concentr)
    ! --- Wind
    CALL init_forcing_input_variable(wspeed_varname, forcing_wind_frequ, forcing_wind_ocean, &
      & forcing_wind_file_prefix, forcing_wind_unit, forcing_wind_const_wspeed, forcing_input(model_id)%wind_speed)
    ! --- deposition variables required for quincy
    SELECT CASE (model%config%model_scheme)
    CASE (MODEL_QUINCY)
      CALL init_forcing_input_variable(nhxdep_varname, forcing_dep_frequ, forcing_nhx_dep_const, &
        & forcing_dep_file_prefix, forcing_dep_unit, forcing_nhx_dep_const, forcing_input(model_id)%nhx_dep)
      CALL init_forcing_input_variable(noydep_varname, forcing_dep_frequ, forcing_noy_dep_const, &
        & forcing_dep_file_prefix, forcing_dep_unit, forcing_noy_dep_const, forcing_input(model_id)%noy_dep)
      CALL init_forcing_input_variable(pdep_varname, forcing_dep_frequ, forcing_p_dep_const, &
        & forcing_dep_file_prefix, forcing_dep_unit, forcing_p_dep_const, forcing_input(model_id)%p_dep)
    END SELECT

  END SUBROUTINE init_forcing

  ! ======================================================================================================= !
  !> Provides forcing data for const, timestep or daily forcing (reads data from netcdf if required)
  !>
  SUBROUTINE get_data_and_update_field(startblk, endblk, startidx, endidx, current_datetime, model_id, &
    &                                  forcing_input_variable, current_forcing_array)

    INTEGER,                   INTENT(IN)    :: startblk, endblk, startidx(:), endidx(:), model_id
    TYPE(t_datetime), POINTER, INTENT(in)    :: current_datetime
    TYPE(input_variable_type), INTENT(INOUT) :: forcing_input_variable
    REAL(wp),                  INTENT(OUT)   :: current_forcing_array(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_data_and_update_field'

    LOGICAL ::  new_year_forcing     ! jsb3: defined comparing current_date with next_date, because with stepwise forcing
    LOGICAL ::  new_day_forcing      ! the dates in the forcing file shall correspond to the dates of stepwise output.

    INTEGER           :: year, month, day, hour, minute, second, ydate, cyclic_year
    INTEGER           :: frequency, currentTimeStep, thisDayStart, thisDayEnd, ith_timestep, index_time, nc_var_id
    LOGICAL           :: do_advance
    CHARACTER(len=8)  :: year_label
    CHARACTER(len=NF90_MAX_NAME)  :: readUnit

    INTEGER :: nproma, nblks
    REAL(wp), POINTER :: return_pointer(:,:,:) !temporary pointer
    REAL(wp) :: hlp_r(1)

    TYPE(t_jsb_model), POINTER :: model

    nproma = SIZE(current_forcing_array,1)
    nblks  = SIZE(current_forcing_array,2)

    current_forcing_array(:,:) = 0._wp

    model => Get_model(model_id)

    !---------------------------
    frequency = forcing_input_variable%frequency

    ! ----- in case of constant forcing, only set correct constants and retun
    IF (frequency == CONST_) THEN
      current_forcing_array(:,:) = forcing_input_variable%constantValue
      CALL set_ocean_cells(startblk, endblk, startidx, endidx, model_id, &
        & forcing_input_variable%ocean_value, current_forcing_array)
      CALL set_miss_cells(startblk, endblk, startidx, endidx, model_id, &
        & forcing_input_variable%missval, forcing_input_variable%ocean_value, current_forcing_array)

      RETURN
    ENDIF

    ! ----- times
    new_year_forcing = is_newyear(current_datetime, delta_time)
    new_day_forcing  = is_newday (current_datetime, delta_time)

    CALL get_date_components(current_datetime, year, month, day, hour, minute, second)

    ! Find year of forcing file and date for new datetime
    do_advance = .TRUE.
    IF (forcing_options(model_id)%cyclic_nyears > 0) THEN
      cyclic_year = forcing_options(model_id)%cyclic_start_year &
        &           + MOD(year - get_year_at_experiment_start(), forcing_options(model_id)%cyclic_nyears)
      ! If current year is a leap year but year from forcing file is not, read in Feb 28 twice
      IF (get_year_length(year) == 366 .AND. get_year_length(cyclic_year) == 365) THEN
        IF (month == 2 .AND. day == 29) THEN
          day = 28
          do_advance = .FALSE.
        END IF
      END IF
      year = cyclic_year
    END IF

    WRITE(year_label,'(i8.4)') year
    ydate = ISIGN(1,year)*(IABS(year)*10000+month*100+day)
    IF (forcing_input_variable%frequency == TIMESTEP_) THEN
      ith_timestep = IMerge_HMS2Sec(hour, minute, second)/int(delta_time) + 1 !+1 because 00:00/X => 0; and 1 based vector
    ELSE IF (forcing_input_variable%frequency == SUBDAILY_) THEN
      ith_timestep = IMerge_HMS2Sec(hour, minute, second)/ &
        & int(86400._wp/forcing_options(model_id)%forcing_steps_per_day) + 1 !+1 because 00:00/X => 0; and 1 based vector
    END IF

    ! ----- CO2 ghg_scenario exception: in case of ghg_scenario a single file is used, which contains one CO2 value per year
    IF ((TRIM(forcing_input_variable%variable_name) == CO2_varname) .AND. (frequency == GHG_SCENARIO_)) THEN
      IF (lstart .OR. lresume) THEN
        IF (my_process_is_stdio()) THEN
        ! in case of ghg_scenario the forcing_input_variable%file_name_prefix is expected to contain the whole filename
          IF (.NOT. forcing_input_variable%input_file%is_open) forcing_input_variable%input_file = &
            &  jsb_netcdf_open_input(TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)))
          CALL get_time_dimension(forcing_input_variable%input_file, forcing_input_variable)
          ! Assert: test that the netcdf file unit for this variable is the expected unit or at least equals unspecified_exp_unit
!          CALL test_unit(forcing_input_variable%input_file%file_id, TRIM(forcing_input_variable%variable_name), &
!            & TRIM(forcing_input_variable%expectedUnit), readUnit)
          CALL forcing_input_variable%input_file%Close()
        END IF

        IF(ASSOCIATED(forcing_input_variable%data_read)) DEALLOCATE(forcing_input_variable%data_read)
        ALLOCATE(forcing_input_variable%data_read(nproma, nblks, forcing_options(model_id)%forcing_steps_per_day))
      ENDIF

      IF (new_year_forcing .OR. lstart .OR. lresume) THEN
        IF (my_process_is_stdio()) THEN
          ! find current year in the timevalues
          currentTimeStep = 1
          DO WHILE (currentTimeStep <= forcing_input_variable%length_time_series)
            IF ( year .EQ. INT(forcing_input_variable%timevalues(currentTimeStep)) ) THEN
              EXIT
            ENDIF
            currentTimeStep = currentTimeStep + 1
          ENDDO
          ! Assert that year has been found
          IF ( year .NE. INT(forcing_input_variable%timevalues(currentTimeStep)) ) THEN
            CALL finish(TRIM(routine), 'The ghg_scenario file "'//TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)) &
              &//'" does not contain the required input year '//TRIM(ADJUSTL(year_label)))
          ENDIF

          forcing_input_variable%input_file = &
            & jsb_netcdf_open_input(TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)))
          ! Read globally constant CO2 value
          CALL nf(nf90_inq_varid(forcing_input_variable%input_file%file_id, CO2_varname, nc_var_id), routine)
          CALL nf(nf90_get_var(forcing_input_variable%input_file%file_id, nc_var_id, hlp_r, &
            & (/1,1,currentTimeStep/), (/1,1,1/)), routine)
          CALL forcing_input_variable%input_file%Close()
        ENDIF

        IF(my_process_is_mpi_parallel()) THEN
          CALL p_bcast(hlp_r, p_io, mpi_comm)
        ENDIF
        forcing_input_variable%data_read(:,:,1) = hlp_r(1)
        WRITE (message_text,*) 'Global GHG scenario CO2 value: ', &
          & forcing_input_variable%data_read(1,1,1)
        CALL message('----- ', message_text)

      ENDIF

      current_forcing_array(:,:) = forcing_input_variable%data_read(:,:,1)
      RETURN
    ENDIF ! forcing_input_variable%variable_name == CO2_varname .AND. frequency == GHG_SCENARIO_

    ! ----- For all other variables and frequencies:
    !Assertion: Rest of the routine only implemented for daily, subdaily or timestep frequency
    IF ( (frequency /= TIMESTEP_) .AND. (frequency /= DAILY_)  .AND. (frequency /= SUBDAILY_)) &
      &  CALL finish(TRIM(routine), 'Violation of assertion: routine only implemented for ' //&
      &  'daily, subdaily or timestep frequency, please check program flow')

    ! TODO: test in case of lresume

    ! --- if new current run or new file - get the timeseries of the forcing file and the current position in the timeseries
    IF (new_year_forcing .OR. lstart .OR. lresume) THEN

      IF (forcing_input_variable%input_file%is_open) CALL forcing_input_variable%input_file%Close()
      IF (my_process_is_stdio()) THEN
        forcing_input_variable%input_file = jsb_netcdf_open_input( &
          & TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)) // TRIM(ADJUSTL(year_label)) // '.nc')

        CALL get_time_dimension(forcing_input_variable%input_file, forcing_input_variable)

        ! Assert: test that the netcdf file unit for this variable is the expected unit or at least equals unspecified_exp_unit
        CALL test_unit(forcing_input_variable%input_file%file_id, TRIM(forcing_input_variable%variable_name), &
          & TRIM(forcing_input_variable%expectedUnit), readUnit)
        ! Read missing value for variable
        IF (forcing_options(model_id)%forcing_set_miss_to_constants) THEN
          CALL get_missval(forcing_input_variable%input_file%file_id, TRIM(forcing_input_variable%variable_name), &
            & forcing_input_variable%missval)
        ENDIF
        ! TODO: further assertions? e.g. that lon lat are correct?

        currentTimeStep = 1
        thisDayStart = 0
        thisDayEnd = 0
        ! Determine current date in the file
        DO WHILE (currentTimeStep <= forcing_input_variable%length_time_series)
          IF ( ydate == INT(forcing_input_variable%timevalues(currentTimeStep)) ) THEN
            IF (thisDayStart .EQ. 0) thisDayStart = currentTimeStep
            thisDayEnd = currentTimeStep
            currentTimeStep = currentTimeStep + 1
          ELSEIF ( ydate < INT(forcing_input_variable%timevalues(currentTimeStep)) ) THEN
            currentTimeStep = currentTimeStep + 1
            EXIT
          ELSE
            currentTimeStep = currentTimeStep + 1
          ENDIF
        ENDDO
        currentTimeStep = currentTimeStep - 1
        ! Assert that current date has been found
        IF (( currentTimeStep > forcing_input_variable%length_time_series) &
          & .OR. (( currentTimeStep == 1) .AND. ( ydate /= INT(forcing_input_variable%timevalues(1)))))  THEN
          WRITE (message_text,*) 'The forcing file "'&
            &//trim(TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)) // TRIM(ADJUSTL(year_label)) // '.nc') &
            &//'" does not contain the required input date ',ydate, ' for the variable "' &
            &//TRIM(forcing_input_variable%variable_name)//'".'
          CALL finish(TRIM(routine), message_text)
        ENDIF

        forcing_input_variable%currentDayStart = thisDayStart
        forcing_input_variable%currentDayEnd = thisDayEnd

      ENDIF

      IF(my_process_is_mpi_parallel()) THEN
        CALL p_bcast(forcing_input_variable%currentDayStart, p_io, mpi_comm)
        CALL p_bcast(forcing_input_variable%currentDayEnd, p_io, mpi_comm)
      ENDIF

      IF(ASSOCIATED(forcing_input_variable%data_read)) DEALLOCATE(forcing_input_variable%data_read)
      ALLOCATE(forcing_input_variable%data_read(nproma, nblks, forcing_options(model_id)%forcing_steps_per_day))

      IF (forcing_input_variable%input_file%is_open) CALL forcing_input_variable%input_file%Close()
      forcing_input_variable%input_file = jsb_netcdf_open_input( &
        & TRIM(ADJUSTL(forcing_input_variable%file_name_prefix)) // TRIM(ADJUSTL(year_label)) // '.nc', model%grid_id)

    ENDIF

    thisDayStart = forcing_input_variable%currentDayStart
    thisDayEnd = forcing_input_variable%currentDayEnd

    ! --- if its a new day or new run - read the data
    IF (new_day_forcing .OR. lstart .OR. lresume) THEN
      ! check timeseries
      IF (my_process_is_stdio()) THEN
        IF (TRIM(forcing_input_variable%variable_name) == 'precip') THEN
          WRITE(message_text,*) ' Reading forcing for ', ydate
          CALL message('----- ', message_text)
        END IF
        !Assertion: currentTimeStep needs to be contained in the file
        IF(thisDayStart > forcing_input_variable%length_time_series) CALL finish(TRIM(routine), &
          & 'Day not found: ' //TRIM(int2string(ydate)))

        !Assertion: check for required number of timesteps
        IF( frequency == TIMESTEP_ .OR. frequency == SUBDAILY_) THEN
          IF(thisDayEnd-thisDayStart+1 /= forcing_options(model_id)%forcing_steps_per_day) CALL finish(TRIM(routine), &
             & 'Forcing file does not contain the demanded number of time steps per day: '  &
             & //TRIM(int2string(forcing_options(model_id)%forcing_steps_per_day)) //' -- but: '                      &
             & //TRIM(int2string(thisDayEnd-thisDayStart+1)))
        ELSEIF( frequency == DAILY_ ) THEN
           IF(thisDayStart /= thisDayEnd) CALL finish(TRIM(routine), &
             & 'Forcing file contains more than one time step per day: ' //TRIM(int2string(thisDayEnd-thisDayStart+1)) &
             & //' -- daily reading only implemented for files with one time step per day.')
        ENDIF
      ENDIF
      ! get data
      return_pointer => forcing_input_variable%input_file%Read_2d_time(forcing_input_variable%variable_name, &
        & start_time_step=thisDayStart, end_time_step=thisDayEnd, fill_array=forcing_input_variable%data_read)

      ! Set the ocean cells of the read data to the constant ocean value and increase the start and end for the next day to be read
      IF (forcing_input_variable%frequency == TIMESTEP_ .OR. forcing_input_variable%frequency == SUBDAILY_) THEN
        DO index_time = 1, forcing_options(model_id)%forcing_steps_per_day
          CALL set_ocean_cells(startblk, endblk, startidx, endidx, model_id, &
            & forcing_input_variable%ocean_value, forcing_input_variable%data_read(:,:,index_time))
          CALL set_miss_cells(startblk, endblk, startidx, endidx, model_id, &
            & forcing_input_variable%missval, forcing_input_variable%ocean_value, &
            & forcing_input_variable%data_read(:,:,index_time))
        ENDDO
        IF (do_advance) THEN
          forcing_input_variable%currentDayStart=forcing_input_variable%currentDayEnd+1
          forcing_input_variable%currentDayEnd=forcing_input_variable%currentDayEnd+forcing_options(model_id)%forcing_steps_per_day
        END IF
      ELSE IF (forcing_input_variable%frequency == DAILY_) THEN
        CALL set_ocean_cells(startblk, endblk, startidx, endidx, model_id, &
          & forcing_input_variable%ocean_value, forcing_input_variable%data_read(:,:,1))
        CALL set_miss_cells(startblk, endblk, startidx, endidx, model_id, &
          & forcing_input_variable%missval, forcing_input_variable%ocean_value, forcing_input_variable%data_read(:,:,1))
        IF (do_advance) THEN
          forcing_input_variable%currentDayStart=forcing_input_variable%currentDayStart+1
          forcing_input_variable%currentDayEnd=forcing_input_variable%currentDayEnd+1
        END IF
      ENDIF
    ENDIF

    IF(forcing_input_variable%frequency == TIMESTEP_ .OR. forcing_input_variable%frequency == SUBDAILY_) THEN
      current_forcing_array = forcing_input_variable%data_read(:,:,ith_timestep)
    ELSEIF(forcing_input_variable%frequency == DAILY_) THEN
      current_forcing_array = forcing_input_variable%data_read(:,:,1)
    ENDIF

  END SUBROUTINE get_data_and_update_field

  !>
  !! Sets given values for input_variable_type fields in forcing_input
  !!
  SUBROUTINE init_forcing_input_variable(variable_name, forcing_frequ_keyword, forcing_ocean_value, &
    & forcing_file_prefix, expected_unit, constantValue, forcing_input_variable)

    CHARACTER(len=*),INTENT(IN) :: variable_name
    CHARACTER(len=*),INTENT(IN) :: forcing_frequ_keyword
    REAL(wp),        INTENT(IN) :: forcing_ocean_value
    CHARACTER(len=*),INTENT(IN) :: forcing_file_prefix
    CHARACTER(len=*),INTENT(IN) :: expected_unit
    REAL(wp),        INTENT(IN) :: constantValue
    TYPE(input_variable_type),INTENT(INOUT) :: forcing_input_variable ! structure filled in this subroutine

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_forcing_input_variable'

    forcing_input_variable%frequency = frequency_key_to_constant_val(variable_name, forcing_frequ_keyword)
    forcing_input_variable%variable_name = variable_name
    forcing_input_variable%file_name_prefix = forcing_file_prefix
    forcing_input_variable%ocean_value = forcing_ocean_value
    forcing_input_variable%expectedUnit = expected_unit
    forcing_input_variable%constantValue = constantValue
    IF (ASSOCIATED(forcing_input_variable%timevalues)) DEALLOCATE(forcing_input_variable%timevalues)
    NULLIFY(forcing_input_variable%timevalues)

  END SUBROUTINE init_forcing_input_variable

  ! ======================================================================================================= !
  !>
  !> Assigns forcing frequency keyword to integer constant
  !>
  INTEGER FUNCTION frequency_key_to_constant_val(variable_name, forcing_frequ_keyword)

    CHARACTER(len=*),INTENT(IN) :: variable_name
    CHARACTER(len=*),INTENT(IN) :: forcing_frequ_keyword

    CHARACTER(len=*), PARAMETER :: routine = modname//':frequency_key_to_constant_val'

    SELECT CASE(TRIM(tolower(forcing_frequ_keyword)))
    CASE("monthly")
      CALL finish(TRIM(routine),'Reading files with monthly frequency is currently not supported, please implement.')
    CASE("daily")
      frequency_key_to_constant_val = DAILY_
    CASE("subdaily")
      frequency_key_to_constant_val = SUBDAILY_
    CASE("const")
      frequency_key_to_constant_val = CONST_
    CASE("timestep")
      frequency_key_to_constant_val = TIMESTEP_
    CASE default
      ! Test for ghg scenario exception for CO2
      IF ((variable_name == CO2_varname) .AND. (TRIM(tolower(forcing_frequ_keyword)) == "ghg_scenario")) THEN
        frequency_key_to_constant_val = GHG_SCENARIO_
      ELSE
        CALL finish(TRIM(routine), 'Please specify a valid reading frequency for the variable '//TRIM(variable_name))
      ENDIF
    END SELECT
    CALL message(TRIM(routine), 'Frequency for forcing variable "'//TRIM(variable_name)//'" is: '//TRIM(forcing_frequ_keyword))
  END FUNCTION frequency_key_to_constant_val

  ! ======================================================================================================= !
  !>
  !> Gets length and content of the time dimension
  !>
  SUBROUTINE get_time_dimension(input_file, forcing_input_variable)

    TYPE(t_input_file),        INTENT(IN)    :: input_file
    TYPE(input_variable_type), INTENT(INOUT) :: forcing_input_variable

    INTEGER                                  :: nc_dim_id, nc_var_id

    CHARACTER(LEN=*), PARAMETER              :: routine = modname//"::get_time_dimension"

    CALL nf(nf90_inq_dimid(input_file%file_id, 'time', nc_dim_id), routine)
    CALL nf(nf90_inquire_dimension(input_file%file_id, nc_dim_id, len = forcing_input_variable%length_time_series), routine)
    IF(ASSOCIATED(forcing_input_variable%timevalues)) DEALLOCATE(forcing_input_variable%timevalues)
    ALLOCATE(forcing_input_variable%timevalues(forcing_input_variable%length_time_series))

    ! We need the timevalues for get_data_and_update_field
    CALL nf(nf90_inq_varid(input_file%file_id, 'time', nc_var_id), routine)
    CALL nf(nf90_get_var(input_file%file_id, nc_var_id, forcing_input_variable%timevalues, &
      & (/ 1 /), (/ forcing_input_variable%length_time_series /)), routine)

  END SUBROUTINE get_time_dimension

  ! ======================================================================================================= !
  !>
  !> Tests if the read unit equals the expected unit (with some exception tests)
  !>
  !> TODO: via mo_jsb_io_netcdf_iface / mo_jsb_io_netcdf functions / mo_input?
  !>
  SUBROUTINE test_unit(file_id, variable_name, expectedUnit, readUnit)

    INTEGER, INTENT(IN) :: file_id
    CHARACTER(len=*),INTENT(IN) :: variable_name
    CHARACTER(len=*),INTENT(IN) :: expectedUnit
    CHARACTER(len=*),INTENT(OUT) :: readUnit

    CHARACTER(len=*), PARAMETER :: routine = modname//':test_unit'

    INTEGER :: nc_var_id, status

    IF (TRIM(expectedUnit) .NE. TRIM(unspecified_exp_unit)) THEN
      CALL nf(nf90_inq_varid(file_id, TRIM(variable_name), nc_var_id), routine)
      status = nf90_get_att(file_id, nc_var_id, name_of_unit_attribute, readUnit)
      IF (status /= NF90_NOERR) THEN
        CALL message(TRIM(routine), &
             & 'WARNING: netcdf error when accessing units attribute of variable '//TRIM(variable_name) &
             & //'(expected '//TRIM(expectedUnit)//')')
      ELSE
        IF(TRIM(readUnit) .NE. TRIM(expectedUnit)) THEN
          IF(TRIM(readUnit) .EQ. TRIM(unspecified_exp_unit)) THEN
            CALL message(TRIM(routine), &
              & 'WARNING: Forcing file does not contain the expected unit "' //TRIM(expectedUnit) &
              & //'" but: "'//TRIM(unspecified_exp_unit)//'" for the variable '//TRIM(variable_name))
          ELSEIF (SCAN(readUnit, ACHAR(0)) > 0) THEN
            CALL message(TRIM(routine), 'WARNING: Forcing file contains unit with NULL character for variable '&
              & //TRIM(variable_name)//' the expected unit "' //TRIM(expectedUnit) &
              & //'" has therefore not been tested against the unit read from the file  "'//TRIM(readUnit)//'".')
          ELSE
            CALL finish(TRIM(routine), &
              & 'Forcing file does not contain the expected unit "' //TRIM(expectedUnit) &
              & //'" but: "'//TRIM(readUnit)//'" for the variable '//TRIM(variable_name))
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE test_unit

  ! ======================================================================================================= !
  !>
  !> Reads and returns missing value for variable fields
  !>
  SUBROUTINE get_missval(file_id, variable_name, missval)

    INTEGER,         INTENT(IN)  :: file_id
    CHARACTER(len=*),INTENT(IN)  :: variable_name
    REAL(wp),        INTENT(OUT) :: missval

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_missval'

    INTEGER :: nc_var_id, status

    CALL nf(nf90_inq_varid(file_id, TRIM(variable_name), nc_var_id), routine)
    status = nf90_get_att(file_id, nc_var_id, 'missing_value', missval)

    IF (status == NF90_NOERR) THEN
      CALL message(TRIM(routine), &
        & 'Found missing value '//TRIM(real2string(missval))//' for variable '//TRIM(variable_name))
    ELSE
      status = nf90_get_att(file_id, nc_var_id, '_FillValue', missval)
      IF (status == NF90_NOERR) THEN
        CALL message(TRIM(routine), &
          & 'Found missing value '//TRIM(real2string(missval))//' for variable '//TRIM(variable_name))
      ELSE
        CALL finish(TRIM(routine), &
          & 'No missing value found for variable '//TRIM(variable_name))
      ENDIF
    ENDIF

  END SUBROUTINE get_missval

  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  SUBROUTINE finalize_external_forcing(no_of_models)

    INTEGER, INTENT(in) :: no_of_models

    INTEGER :: jg

    DO jg=1,no_of_models
      IF(forcing_options(jg)%air_temperature_as_timestep) THEN
        IF (forcing_input(jg)%air_temp%input_file%is_open) CALL forcing_input(jg)%air_temp%input_file%Close()
      ELSE
        IF (forcing_input(jg)%air_temp_min%input_file%is_open) CALL forcing_input(jg)%air_temp_min%input_file%Close()
        IF (forcing_input(jg)%air_temp_max%input_file%is_open) CALL forcing_input(jg)%air_temp_max%input_file%Close()
      END IF
      IF (forcing_input(jg)%precipitation%input_file%is_open) CALL forcing_input(jg)%precipitation%input_file%Close()

      SELECT CASE(forcing_options(jg)%type_of_qair_forcing)
      CASE(RH_)
        IF (forcing_input(jg)%rel_humidity%input_file%is_open) CALL forcing_input(jg)%rel_humidity%input_file%Close()
      CASE(QAIR_)
        IF (forcing_input(jg)%qair%input_file%is_open) CALL forcing_input(jg)%qair%input_file%Close()
      END SELECT

      IF (forcing_input(jg)%shortwave%input_file%is_open) CALL forcing_input(jg)%shortwave%input_file%Close()

      IF (forcing_input(jg)%longwave%input_file%is_open) CALL forcing_input(jg)%longwave%input_file%Close()

      IF (forcing_input(jg)%wind_speed%input_file%is_open) CALL forcing_input(jg)%wind_speed%input_file%Close()

      IF (forcing_input(jg)%CO2_concentr%input_file%is_open) CALL forcing_input(jg)%CO2_concentr%input_file%Close()

      IF (forcing_input(jg)%nhx_dep%input_file%is_open) CALL forcing_input(jg)%nhx_dep%input_file%Close()

      IF (forcing_input(jg)%noy_dep%input_file%is_open) CALL forcing_input(jg)%noy_dep%input_file%Close()

      IF (forcing_input(jg)%p_dep%input_file%is_open) CALL forcing_input(jg)%p_dep%input_file%Close()

    END DO

    DEALLOCATE(forcing_input, forcing_options)

  END SUBROUTINE finalize_external_forcing

  ! ======================================================================================================= !
  !>
  !> Converts daily relative humidity to specific humidity considering the daily course of temperature
  !>
  !> @NOTE Using this routine only makes sense together with the generated diurnal cycle of temperature as in
  !>       'SUBROUTINE instantly_from_daily_temp_2()' otherwise the results will be faulty!
  !>
  SUBROUTINE spec_humidity_from_rel_humidity(ice, sinlat, coslat, rel_humidity, Tmin, Tmax, elevation, qair)

    INTEGER, INTENT(IN)  :: ice
    REAL(wp),INTENT(IN)  :: sinlat(1:ice)
    REAL(wp),INTENT(IN)  :: coslat(1:ice)
    REAL(wp),INTENT(IN)  :: rel_humidity(1:ice)  !! mean daily relative humidity
    REAL(wp),INTENT(IN)  :: Tmin(1:ice)          !! minimum temperature (= temperature at sunrise)
    REAL(wp),INTENT(IN)  :: Tmax(1:ice)          !! maximum temperature (= at 2pm)
    REAL(wp),INTENT(IN)  :: elevation(1:ice)     !! height above sea leavel
    REAL(wp),INTENT(OUT) :: qair(1:ice)          !! the specific humidity of the day

    ! constants
    REAL(wp), PARAMETER :: pio180 = pi/180._wp     ! conversion factor from degrees to radians ("pi over 180")
    REAL(wp), PARAMETER :: rat1=14._wp/24._wp      ! day fraction at which maximum temperature shall occur
    REAL(wp), PARAMETER :: rat2=4._wp/24._wp       ! minimum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat3=20._wp/24._wp      ! maximum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat4=2._wp/24._wp       ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                                   ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2
    ! Temperature calculations
    REAL(wp) :: declination       ! solar declination
    REAL(wp) :: cpds, spds        ! cosine and sine contributions to zenith angle
    REAL(wp) :: fract_daylight     ! fraction of day between sunrise and sunset
    REAL(wp) :: fract_since_sunset ! time between current time step and sunset expressed as fraction of day
    REAL(wp) :: fract_sunrise      ! time of sunrise expressed as fraction of day since midnight
    REAL(wp) :: fract_sunset       ! time of sunset expressed as fraction of day since midnight
    REAL(wp) :: temp_sunset       ! temperature at sunset
    REAL(wp) :: dtmp              ! mean day temperature
    REAL(wp) :: dtran             ! day temperature range: maximum - minimum
    REAL(wp) :: sd                ! modulation factor for temperature
    REAL(wp) :: arg,nstep
    REAL(wp),ALLOCATABLE :: tmp(:)  ! diurnal temperature

    ! Qair calculations
    REAL(wp) :: rh_test,rh_last,qsat_last ! temporary RH estimates
    REAL(wp),ALLOCATABLE :: qsat(:)      ! diurnal saturated specific humidity
    REAL(wp),ALLOCATABLE :: inv_qsat(:)  ! inverse diurnal saturated specific humidity
    REAL(wp),ALLOCATABLE :: qair_test(:) ! diurnal actual specific humidity, temporary
    REAL(wp),ALLOCATABLE :: pressure(:)  ! diurnal air pressure
    REAL(wp) :: tstep        ! fractional time step of the day

    ! helpers
    INTEGER    :: i, j, jj
    INTEGER    :: offset(1)

#ifdef _OPENACC
    CALL finish('spec_humidity_from_rel_humidity', 'Not ported to GPU, yet. Stop.')
#endif

    CALL inquire_declination(declination)

    nstep = (86400.0_wp/delta_time)

    IF(.NOT.(ALLOCATED(tmp))) ALLOCATE(tmp(time_steps_per_day))
    IF(.NOT.(ALLOCATED(pressure))) ALLOCATE(pressure(time_steps_per_day))
    IF(.NOT.(ALLOCATED(qsat))) ALLOCATE(qsat(time_steps_per_day))
    IF(.NOT.(ALLOCATED(inv_qsat))) ALLOCATE(inv_qsat(time_steps_per_day))
    IF(.NOT.(ALLOCATED(qair_test))) ALLOCATE(qair_test(time_steps_per_day))

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) PRIVATE(tmp, pressure, qsat, inv_qsat, qair_test)
    DO i = 1,ice
      spds = sinlat(i)*SIN(declination)  ! sine contribution to solar zenith angle
      cpds = coslat(i)*COS(declination)  ! cosine contribution to solar zenith angle

      arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
      IF (arg > 1._wp) THEN               ! polar night:
        fract_daylight = 0._wp            ! ... day length is zero
      ELSE IF (arg < -1._wp) THEN         ! polar day:
        fract_daylight = 1._wp            ! ... day length is the whole day
      ELSE                                ! normal day / night:
        fract_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
      END IF
      dtmp = 0.5_wp*(Tmax(i)+Tmin(i)) ! The mean day temperature
      dtran = Tmax(i)- Tmin(i)        ! temperature range
      DO j = 1,time_steps_per_day
        tstep=j/nstep
        ! -------- normal day ---------------------------------------------
        IF (fract_daylight>=rat2 .AND. fract_daylight<=rat3) THEN
          fract_sunrise = 0.5_wp - fract_daylight/2._wp             ! time of sunrise expressed as fraction of day since midnight
          fract_sunset = fract_sunrise+fract_daylight                ! time of sunset expressed as fraction of day since midnight

          ! For time steps at daylight ...
          IF (tstep>fract_sunrise .AND. tstep<fract_sunset) THEN
            sd = COS(pi*(tstep-rat1)/(fract_daylight/2._wp+rat4)) ! modulation factor for day-temperature such that
                                                                 !   temperature is at maximum at rat1 (2 pm) and at average at
                                                                 !   sunrise-rat4/2 and sunset+rat4/2
            tmp(j) = dtmp + sd*dtran/2._wp                       ! temperature modulated by temperature range

          ! For time steps at night ...
          ELSE
            sd = COS(pi*(fract_sunset-rat1)/(fract_daylight/2._wp+rat4)) ! modulation factor for temperature at sunset
            temp_sunset = dtmp + dtran/2._wp*sd                        ! temperature at sunset
            fract_since_sunset = MOD(tstep-fract_sunset+1._wp,1._wp)     ! day fraction since sunset
            tmp(j) = Tmin(i)+(temp_sunset-Tmin(i))*(1._wp-fract_since_sunset/(1._wp-fract_daylight))! linear interpolation
                                                                     ! ... of temperature such that temperature is equal to the
                                                                     ! ... temperature at sunset and drops until the minimum ...
                                                                     ! ... temperature is met at sunrise.
          END IF
        ! --- long day (close to polar day) ------------------------------------------------------
        ELSEIF (fract_daylight>rat3) THEN
          sd = COS(pi*(tstep-rat1)/(fract_daylight/2._wp+rat4)) ! modulation factor for day-temperature such that ...
                                                               ! ... temperature is at maximum at rat1 (2 pm) and at average
                                                               ! ... at sunrise-rat4/2 and sunset+rat4/2
          tmp(j) = dtmp + dtran/2._wp*sd                       ! temperature modulated by temperature range
        ! ------------------------------- short day (close to polar night) ---------------------------------------------------
        ELSE
          tmp(j) = dtmp                                        ! temperature is set to mean temperture at all time steps
        END IF
        tmp(j) = tmp(j) + tmelt

        ! from function bottom pressure; in jsb3: p0sl_bg = p_sealevel; grav = Gravity; rd = GasConstantDryAir
        pressure(j) = &
          & p0sl_bg * ( 1._wp - elevation(i)*gamma /(tmp(j) + elevation(i)*gamma))**(grav/rd/gamma)

        CALL sat_specific_humidity_scalar(tmp(j), pressure(j), qsat(j)) ! present saturated qair

      ENDDO

      ! implied RH from minimum saturated specific humidity
      rh_test=SUM(MINVAL(qsat(:))/qsat(:))/nstep
      IF( rel_humidity(i) <= rh_test ) THEN
        ! Simply case, relative humitidy that assuming a constant qair there is never
        ! saturation during the day. In this case, qair can be inferred by inverting
        ! RH = 1/n * sum ( QAIR / QSAT ) = QAIR / n * sum ( 1 / QSAT )
        ! -> QAIR = n * RH / sum ( 1 / QSAT )
        inv_qsat(:) = 1._wp / qsat(:)
        qair(i) = nstep * rel_humidity(i) / SUM( inv_qsat(:) )
      ELSE
        ! The case where RH implies saturation at night. Under these conditions
        ! search for the qair value which is closed to giving the prescribed RH
        ! taking account of the implied period during the day where saturation occurs.
        ! Surely there's a better way to do this, but I cannot be bothered now.
        offset(:)=MAXLOC(qsat(:))
        IF( rel_humidity(i) >= 1._wp ) THEN
           qair(i)=qsat(offset(1))
        ELSE
          j=0
          jj=offset(1)
          rh_test=1._wp
          DO WHILE ( rh_test > rel_humidity(i) .AND. j < time_steps_per_day )
            ! save value from last time
            qsat_last=qsat(jj)
            rh_last=rh_test

            ! find tome step after hottest time of the day
            IF(j+offset(1)>time_steps_per_day) THEN
              jj=j+offset(1)-time_steps_per_day
            ELSE
              jj=j+offset(1)
            ENDIF

            ! create qair values constrained by qsat at low temperatures
            WHERE(qsat(:) <= qsat(jj))
              qair_test(:)=qsat(:)
            ELSEWHERE
              qair_test(:)=qsat(jj)
            END WHERE

            ! test whether target RH is reached
            rh_test = SUM(qair_test(:) / qsat(:)) /nstep
            j=j+1
          ENDDO
          ! The above solution is inprecise. Attempt to reconcile this by linear approximation
          qair(i)=qsat(jj)+(qsat_last-qsat(jj))*(rel_humidity(i)-rh_test)/(rh_last-rh_test)
        ENDIF
      ENDIF
    ENDDO
    !$ACC END PARALLEL LOOP

    IF((ALLOCATED(tmp))) DEALLOCATE(tmp)
    IF((ALLOCATED(pressure))) DEALLOCATE(pressure)
    IF((ALLOCATED(qsat))) DEALLOCATE(qsat)
    IF((ALLOCATED(inv_qsat))) DEALLOCATE(inv_qsat)
    IF((ALLOCATED(qair_test))) DEALLOCATE(qair_test)

  END SUBROUTINE spec_humidity_from_rel_humidity

  ! ======================================================================================================= !
  !>
  !> Returns saturation specific humidity for given temperature and pressure
  !> (at some atmospheric level or at surface)
  !>
  !>   Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case,
  !>   but the saturation vapor pressure of water over a liquid surface resp. ice
  !>   (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)
  !>
  SUBROUTINE sat_specific_humidity(temp, pressure, qsat)

    USE mo_phy_schemes, ONLY: qsat_water

    REAL(wp), INTENT(IN) :: temp(:)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(:)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sat_specific_humidity'

    INTEGER :: ic, nc

    nc = SIZE(temp)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
#ifndef _OPENACC
      IF (temp(ic) < 50._wp .OR. temp(ic) > 400._wp) THEN
        WRITE (message_text,*) &
          & 'Temperature out of bound: ', temp(ic), 'K. One thing to check: consistency of land-sea mask and forcing.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSE
#endif
        qsat(ic) = qsat_water(temp(ic), pressure(ic), use_convect_tables=.TRUE.)
#ifndef _OPENACC
      ENDIF
#endif
    ENDDO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE sat_specific_humidity

  SUBROUTINE sat_specific_humidity_scalar(temp, pressure, qsat)

    !$ACC ROUTINE SEQ

    USE mo_phy_schemes, ONLY: qsat_water

    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)

    REAL(wp), INTENT(IN) :: temp      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat

    CHARACTER(len=*), PARAMETER :: routine = modname//':sat_specific_humidity_scalar'

    IF (temp < 50._wp .OR. temp > 400._wp) THEN
!       WRITE (message_text,*) &
!          & 'Temperature out of bound: ', temp(i), 'K. One thing to check: consistency of land-sea mask and forcing.'
!        CALL finish(TRIM(routine), TRIM(message_text))
    ELSE
      qsat = qsat_water(temp, pressure, use_convect_tables=.TRUE.)
    ENDIF

  END SUBROUTINE sat_specific_humidity_scalar

  ! ======================================================================================================= !
  !>
  !> Returns saturation specific humidity for given air and vapour pressure
  !>
  !
  ! docu: fun_specific_humidity()
  !
  ! Specific humidity is computed here by assuming that each component of a mixture of dry air and vapor can be considered as an
  ! ideal gas. So, let m_d and m_v the masses of the two components in a fixed volume V. Then the specific humidity "q" is
  ! defined as
  ! (1) q = m_v/(m_v+m_d)
  ! Closely related is the moisture mixing ratio "r"
  ! (2) r = m_v/m_d
  ! Hence,
  ! (3) q = r/(1+r)
  ! Next it is assumed that both, dray air and vapor behave like ideal gases:
  ! (4) p_d*V = m_d*R_d*T and p_v*V = m_v*R_v*T,
  ! where p_d is the partial pressure of the dry air and p_v the vapor pressure, R_d and R_v the gas constants and T the
  ! temperature. Noting that air pressure "p" is the sum of the two partial pressures (p=p_d+p_v) the division of the two gas
  ! equation by each other gives by using (2):
  ! (5) r = eps*p_v/(p-p_v),
  ! where eps=R_d/R_v is the ratio of gas constants. Entering ths into (3) one finally obtains for the specific humidity:
  !                 p_v
  ! (6) q = eps -------------
  !             p-(1-eps)*p_v
  !
  PURE ELEMENTAL FUNCTION fun_specific_humidity(vapour_pressure, air_pressure)
    REAL(wp),INTENT(in) :: vapour_pressure        !< [N/m^2]
    REAL(wp),INTENT(in) :: air_pressure           !< air pressure at bottom [N/m^2]
    REAL(wp)            :: fun_specific_humidity  !< (from [0,1])

    fun_specific_humidity = rdv * vapour_pressure / (air_pressure - (1._wp - rdv) * vapour_pressure)
  END FUNCTION fun_specific_humidity

  SUBROUTINE fun_specific_humidity_2d(vapour_pressure, air_pressure, fun_specific_humidity)
    REAL(wp), INTENT(IN)  :: vapour_pressure(:,:) !! [N/m^2]
    REAL(wp), INTENT(IN)  :: air_pressure(:,:) !! air pressure at bottom [N/m^2]
    REAL(wp), INTENT(OUT) :: fun_specific_humidity(:,:) !! (from [0,1])
    INTEGER               :: nproma, nblks, ic, iblk

#ifdef _OPENACC
    CALL finish('fun_specific_humidity_2d', 'Not ported to GPU, yet. Stop.')
#endif

    nproma = SIZE(fun_specific_humidity,1)
    nblks  = SIZE(fun_specific_humidity,2)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO iblk = 1, nblks
      DO ic = 1, nproma
        fun_specific_humidity(ic,iblk) = rdv * vapour_pressure(ic,iblk) / &
        & (air_pressure(ic,iblk) - (1._wp - rdv) * vapour_pressure(ic,iblk))
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE fun_specific_humidity_2d

  ! ======================================================================================================= !
  !>
  !> Computes vapor pressure from actual and potential evapotranspiration
  !>
  !
  ! docu: vapour_pressure_from_evapor()
  !
  ! Vapor pressure is computed according to a method described in [1] from actual and potential evapotranspiration.
  ! Typically the mean vapor pressure at a day is approximately equal to the saturation vapor pressure at sunrise
  ! "e_sat(T_min)", where temperature is at its minimum. Only in arid regions this approximation gets bad; in these regions
  ! the saturation vapor pressure is not reached at sunrise. The following model accounts for this observation by using the
  ! ratio of actual to potential evapotranspiration as a measure of aridity. Let "E_act(d)" and "E_pot(d)" denote the total
  ! actual and total potential evapotranspiration at day "d" and "fe(d)" the ratio of these values (Eq. (22) in [1]):
  ! (1) fe(d) = E_act(d)/E_pot(d).
  ! For draught situations  E_pot(d)=infty so that fe(d)=0. In the other extreme, fe(d)=1, E_act(d)=E_pot(d), i.e. actual transpi-
  ! ration is only limited by potential evapotranspiration and this means the climate is extremely wet (sufficient soil moisture
  ! and finite capacity of the atmosphere to absorb water vapor). Let the vapor pressure at sunrise be (Eq. (21) in [1])
  ! (2) e0_vap(d) = (h0 + (1-h0)*fe(d-1))*e_sat(T_min),
  ! where "e_sat(T_min)" is the saturation vapor pressure at sunrise temperature "T_min", "fe(d-1)" the yesterdays
  ! evapotranspiration ratio and "h0" a tuning parameter that has the meaning of the humidity at sunrise for total draught
  ! (fe=0). Here the yesterdays value of the evapotranspiration ratio "fe" is taken. The vapor pressure "e_vap(t)" in the actual
  ! timestep "t" is then estimated from the current temperature "T" by (Eq. (20) in [1])
  ! (3) e_vap(t) = e0_vap(d) + h1*fe(d-1)*(e_sat(T)-e0_vap(d)),
  ! where "h1" is another tuning constant that has the meaning of a daily amplitude of vapor pressure under extreme moist
  ! conditions (fe=1). According to [1] one finds h0=0.96 and h1=0.49 (see pp. 65/66).
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  !
  ! Above values for h0 and h1 lead to extremely moist conditions, especially in desert areas.
  ! A comparison of qair from echam/jsbach output with qair calculated in this routine using different h0 values shows best
  ! agreement for h0=0.6. The same h0 is found when comparing CRUNCEP qair with qair calculated from CRUNCEP temperatures in
  ! this module. For desert studies it might be good to reduce h0 even further. As the specific humidity does not have a
  ! significant diurnal cycle, h1 is set to 0.
  !
  SUBROUTINE vapour_pressure_from_evapor(startblk, endblk, startidx, endidx, evapo_ratio_mean, Tmin, Tair, vpress)

    ! USE mo_phy_schemes, ONLY: sat_pres_water
    USE mo_phy_schemes, ONLY: tlucua, jptlucu1, jptlucu2

    INTEGER,  INTENT(IN) :: startblk, endblk, startidx(:), endidx(:)
    REAL(wp), INTENT(IN) :: evapo_ratio_mean(:,:)  !! fe(d-1) of Eq. (1)
    REAL(wp), INTENT(IN) :: Tmin(:,:)              !! minimum temperature (= temperature at sunrise)
    REAL(wp), INTENT(IN) :: Tair(:,:)              !! temperature at current time step
    REAL(wp), INTENT(OUT) :: vpress(:,:)

    REAL(wp), PARAMETER :: h0 = 0.6_wp  !! (0.96_wp) tuning constant for computing vapor pressure
    REAL(wp), PARAMETER :: h1 = 0.0_wp  !! (0.49_wp) tuning constant for computing vapor pressure

    REAL(wp) :: evapo_sat        ! saturation pressure at actual temperature
    REAL(wp) :: evapo_sat_min    ! saturation pressure at minimum temperature of day (= at sunrise)
    REAL(wp) :: e0_vap           ! vapor pressure at sunrise

    INTEGER :: ib, ic, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':vapour_pressure_from_evapor'

    DO ib = startblk, endblk
      ics = startidx(ib)
      ice = endidx  (ib)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic = ics, ice

        !Assertion: the expected temperature range is indirectly defined by the tlucua bounds -> consistency check
#ifndef _OPENACC
        IF (Tmin(ic,ib) + tmelt < 50._wp .OR. Tmin(ic,ib) + tmelt > 400._wp) THEN
          WRITE (message_text,*) &
            & 'min temperature out of range: ', Tmin(ic,ib), ' degC. One thing to check: consistency of land-sea mask and forcing.'
          CALL finish(TRIM(routine), TRIM(message_text))
        ENDIF
        IF (Tair(ic,ib) + tmelt < 50._wp .OR. Tair(ic,ib) + tmelt > 400._wp) THEN
          WRITE (message_text,*) &
            & 'air temperature out of bounds: ', Tair(ic,ib), ' degC. One thing to check: consistency of land-sea mask and forcing.'
          CALL finish(TRIM(routine), TRIM(message_text))
        ENDIF
#endif

        !jsb3: rdv = eps
        evapo_sat_min = tlucua(INT(1000._wp * (Tmin(ic,ib) + tmelt)))/rdv ! saturation pressure at minimum temperature
        evapo_sat     = tlucua(INT(1000._wp * (Tair(ic,ib) + tmelt)))/rdv ! saturation pressure at actual temperature
        ! evapo_sat_min = sat_pres_water(Tmin(i,j) + tmelt) ! saturation pressure at minimum temperature
        ! evapo_sat     = sat_pres_water(Tair(i,j) + tmelt) ! saturation pressure at actual temperature
        e0_vap = evapo_sat_min * (h0 + (1._wp - h0) * evapo_ratio_mean(ic,ib))                               ! Eq. (2) from above
        vpress(ic,ib) = e0_vap + h1 * evapo_ratio_mean(ic,ib) * (evapo_sat - evapo_sat_min) ! Eq. (3) from above

      ENDDO
      !$ACC END PARALLEL
    ENDDO

  END SUBROUTINE  vapour_pressure_from_evapor

  ! ======================================================================================================= !
  !>
  !> Sets all not land (i.e. not notsea ;)) cells to the given constant
  !>
  !
  ! docu:
  !   Some external forcing is only defined for land cells.
  !   If the external forcing is only defined for land cells is specified by the namelist option
  !   "forcing_options%forcing_set_ocean_to_constants". If true: ocean cells need to be replaced
  !   by constants (reasonable values for the different variables) to ensure numeric stability.
  !   If false: forcing is expected for all variables for all cells (e.g. forcing from echam)
  !
  SUBROUTINE set_ocean_cells(startblk, endblk, startidx, endidx, model_id, constant, array)

    USE mo_jsb_model_class, ONLY: t_jsb_model
    USE mo_jsb_class,       ONLY: Get_model

    INTEGER,  INTENT(IN)    :: startblk, endblk, startidx(:), endidx(:)
    INTEGER,  INTENT(IN)    :: model_id
    REAL(wp), INTENT(IN)    :: constant
    REAL(wp), INTENT(INOUT) :: array(:,:)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    INTEGER :: ib, ic, ics, ice

    IF (.NOT. forcing_options(model_id)%forcing_set_ocean_to_constants) RETURN

    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)

!$OMP PARALLEL DO PRIVATE(ib, ic, ics, ice)
    DO ib = startblk, endblk
      ics = startidx(ib)
      ice = endidx  (ib)
      DO ic = ics, ice
        IF (tile%fract(ic,ib) <= 0._wp) array(ic,ib) = constant
      END DO
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE set_ocean_cells

  ! ======================================================================================================= !
  !>
  !> Replace missing values with the given constant
  !>
  !
  ! docu:
  !    Land sea masks of model grids and offline forcing data don't alway match exactly. If
  !    the forcing dataset has missing values on land cell, the model will most likely
  !    crash with temperature bound violations due to the missing values. This routine
  !    replaces the missing value with a constant. Currently, the "ocean" value is
  !    used for this.
  SUBROUTINE set_miss_cells(startblk, endblk, startidx, endidx, model_id, missval, constant, array)

    INTEGER, INTENT(IN)     :: startblk, endblk, startidx(:), endidx(:)
    INTEGER,  INTENT(IN)    :: model_id
    REAL(wp), INTENT(IN)    :: missval, constant
    REAL(wp), INTENT(INOUT) :: array(:,:)

    INTEGER :: ib, ic, ics, ice

    IF (.NOT. forcing_options(model_id)%forcing_set_miss_to_constants) RETURN

    ! Use REAL values for comparison in case the netCDF missing_value
    ! attribute lacks double precision.

!$OMP PARALLEL DO PRIVATE(ib, ic, ics, ice)
    DO ib = startblk, endblk
      ics = startidx(ib)
      ice = endidx  (ib)
      DO ic = ics, ice
        IF (ABS(array(ic,ib) - missval) <= EPSILON(1.0_wp)) array(ic,ib) = constant
      END DO
    END DO
!$OMP END PARALLEL DO

  END SUBROUTINE set_miss_cells

  ! ======================================================================================================= !
  !>
  !> Estimate timestep temperature from daily min- and max temperature
  !>
  !
  ! docu: instantly_from_daily_temp_2()
  !
  ! The same as subroutine instantly_from_daily_temp_1(), but with temperature minimum and maximum as input instead of
  ! mean temperature and temperature range.
  ! In contrast to Marko's subroutine daytemp() (from which this code is derived) the routine
  ! correctly accounts for local time, i.e. temperatures are correct at all longitudes!
  !
  ! More precisely the temperature development is modelled as follows:
  ! Let: T_min = T_mean - 0.5*T_range,  and: T_max = T_mean + 0.5*T_range. Then for a
  !     (i) normal day (more than 4, but less than 20 hours of daylight):
  !         --> during daylight: temperature behaves like a positive sine-wave, such that T_max is reached at 14:00 and
  !             T_min at sunrise.
  !         --> during night: starting at sunset the temperature drops linearly until T_min at sunrise.
  !    (ii) short day (less than 4 hours of daylight):
  !         temperature is set to T_mean at all time steps, i.e. temperature is assumed to be constant.
  !   (iii) long day (more than 20 hours of daylight):
  !         temperature is modelled as in the daylight-case of (i), but for the whole day, even during night hours.
  !
  ! REMARKS:
  !      -- The time computations are based on the Julian/Gregorian calender. Hence there is a slight mismatch between the
  !         correct begin of a year, as defined by the orbit around the sun, and the begin of the year according to the calender.
  !         The mismatch is about 3 days in 10000 years. This affects the computation of daylength in the current routine, and
  !         thus the "phase" of the temperature curve, i.e. that e.g. the time when the temperature starts rising in the morning,
  !         may be somewhat shifted.
  !
  SUBROUTINE instantly_from_daily_temp_2(startblk, endblk, startidx, endidx, lon, sinlat, coslat, dtmp_min, dtmp_max, date, tmp)

    REAL(wp), PARAMETER :: pio180 = pi/180._wp ! conversion factor from degrees to radians ("pi over 180")
    REAL(wp), PARAMETER :: rat1=14._wp/24._wp  ! day fraction at which maximum temperature shall occur
    REAL(wp), PARAMETER :: rat2=4._wp/24._wp   ! minimum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat3=20._wp/24._wp  ! maximum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat4=2._wp/24._wp   ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                               ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2

    INTEGER,  INTENT(IN) :: startblk, endblk, startidx(:), endidx(:)
    REAL(wp), INTENT (IN) :: lon(:,:)       ! minimal temperature at the considered day
    REAL(wp), INTENT (IN) :: sinlat(:,:)    ! maximal range at that day
    REAL(wp), INTENT (IN) :: coslat(:,:)    ! maximal range at that day
    REAL(wp), INTENT (IN) :: dtmp_min(:,:)  ! minimal temperature at the considered day
    REAL(wp), INTENT (IN) :: dtmp_max(:,:)  ! maximal range at that day
    TYPE (t_datetime), POINTER, INTENT (IN)  :: date
    REAL(wp), INTENT (OUT) :: tmp(:,:)      ! estimated temperature


    INTEGER    :: yearday           ! the number of the day of the current timestep in the current year
    REAL(wp)   :: year_instant      ! the time of the current time step in yearday-fraction-of-day format
    REAL(wp)   :: fract_timestep     ! the current time step expressed as the fraction of the current day
    REAL(wp)   :: fract_localtime    ! local time at a particular grid point expressed as the fraction of the current ...
                                    ! ... day since midnight (mod 1
    REAL(wp)   :: fract_UTC_offset   ! offset between UTC and local time (depends on longitude)
    REAL(wp)   :: declination       ! solar declination
    REAL(wp)   :: cpds, spds        ! cosine and sine contributions to zenith angle
    REAL(wp)   :: fract_daylight     ! fraction of day between sunrise and sunset
    REAL(wp)   :: fract_since_sunset ! time between current time step and sunset expressed as fraction of day
    REAL(wp)   :: fract_sunrise      ! time of sunrise expressed as fraction of day since midnight
    REAL(wp)   :: fract_sunset       ! time of sunset expressed as fraction of day since midnight
    REAL(wp)   :: temp_sunset       ! temperature at sunset
    REAL(wp)   :: dtmp              ! mean day temperature
    REAL(wp)   :: dtran             ! day temperature range: maximum - minimum

    REAL(wp)   :: sd                ! modulation factor for temperature
    REAL(wp)   :: arg
    INTEGER    :: ic, ib, ics, ice

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO ib = 1,SIZE(tmp,2)
      DO ic = 1, SIZE(tmp, 1)
        tmp(ic,ib) = 0._wp
      ENDDO
    ENDDO
    !$ACC END PARALLEL LOOP

    CALL inquire_declination(declination)

    year_instant=get_year_day(date)            ! conversion of current time step from day-seconds to yearday-dayfraction format
    yearday = INT(year_instant)                ! number of day of the current time step since the beginning of the current year
    fract_timestep=year_instant-REAL(yearday,wp)   ! fractional part of the day of the current time step since midnight in UTC

    DO ib = startblk, endblk
      ics = startidx(ib)
      ice = endidx(ib)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic = ics, ice
        fract_UTC_offset = lon(ic,ib)/360._wp   ! difference between local time and UTC as fraction of a day
        fract_localtime = MOD(fract_timestep+fract_UTC_offset+2._wp,1._wp) ! local time as fraction of the local day (modulo 1)
        spds = sinlat(ic,ib)*SIN(declination)  ! sine contribution to solar zenith angle
        cpds = coslat(ic,ib)*COS(declination)  ! cosine contribution to solar zenith angle
        arg = -spds/cpds                        ! faulty?: for cpds=0 arg gets infinite
        IF (arg > 1._wp) THEN               ! polar night:
          fract_daylight = 0._wp            ! ... day length is zero
        ELSEIF (arg < -1._wp) THEN         ! polar day:
          fract_daylight = 1._wp            ! ... day length is the whole day
        ELSE                                ! normal day / night:
          fract_daylight = ACOS(arg)/pi    ! ... day length expressed as fraction of the whole day
        ENDIF

        !! BLARPP empirical correction to account for the fact that the scheme does actually not
        !! consered dtmp=0.5*(max+min). The night-day turn temperature is adjusted such that the mean
        !! temperature is within 0.05 K of the 0.5*(tmax+tmin) temperature
        IF (fract_daylight>=rat2 ) THEN ! -------- normal day ---------------------------------------------
          IF(fract_daylight < 0.36_wp)THEN
            arg=-8._wp * (fract_daylight-0.36_wp)**2 + 0.685_wp
          ELSE
            arg=-0.8_wp * (fract_daylight-0.36_wp)**2 + 0.685_wp
          ENDIF
        ELSE
          arg = 0.5_wp
        ENDIF

        dtmp = arg * dtmp_max(ic,ib) + (1._wp-arg) * dtmp_min(ic,ib) ! The dusk day temperature
        !dtmp = 0.5_wp*(dtmp_max(ic,ib)+dtmp_min(ic,ib)) ! The mean day temperature
        dtran = dtmp_max(ic,ib)- dtmp_min(ic,ib)        ! temperature range
        IF (fract_daylight>=rat2 .AND. fract_daylight<=rat3) THEN ! -------- normal day ---------------------------------------------
          fract_sunrise = 0.5_wp - fract_daylight/2._wp             ! time of sunrise expressed as fraction of day since midnight
          fract_sunset = fract_sunrise+fract_daylight                ! time of sunset expressed as fraction of day since midnight
          IF (fract_localtime > fract_sunrise .AND. fract_localtime <= fract_sunset) THEN  ! For time steps at daylight ...
            sd = COS(pi*(fract_localtime-rat1)/(fract_daylight/2._wp+rat4)) ! ... modulation factor for day-temperature such that
                                                                    ! .. temperature is at maximum at rat1 (2 pm) and at average at
                                                                    ! .. sunrise-rat4/2 and sunset+rat4/2
            IF (sd >= 0) THEN
              tmp(ic,ib) = dtmp + sd*(dtmp_max(ic,ib)-dtmp)                       ! temperature modulated by temperature range
            ELSE
              tmp(ic,ib) = dtmp + sd*(dtmp-dtmp_min(ic,ib))                       ! temperature modulated by temperature range
            ENDIF
            !tmp(ic,ib) = dtmp + sd*dtran/2._wp                       ! temperature modulated by temperature range
          ELSE                                                 ! For time steps at night ...
            sd = COS(pi*(fract_sunset-rat1)/(fract_daylight/2._wp+rat4)) ! modulation factor for temperature at sunset
            IF (sd >= 0) THEN
              temp_sunset = dtmp + sd*(dtmp_max(ic,ib)-dtmp)                       ! temperature modulated by temperature range
            ELSE
              temp_sunset = dtmp + sd*(dtmp-dtmp_min(ic,ib))                       ! temperature modulated by temperature range
            ENDIF
            !temp_sunset = dtmp + dtran/2._wp*sd                     ! temperature at sunset
            fract_since_sunset = MOD(fract_localtime-fract_sunset+1._wp,1._wp) ! day fraction since sunset
            tmp(ic,ib) = dtmp_min(ic,ib)+(temp_sunset-dtmp_min(ic,ib)) &
              & *(1._wp-fract_since_sunset/(1._wp-fract_daylight))! linear interpolation
                                                                ! ... of temperature such that temperature is equal to the
                                                                ! ... temperature at sunset and drops until the minimum ...
                                                                ! ... temperature is met at sunrise.
          ENDIF
        ELSEIF (fract_daylight>rat3) THEN ! --- long day (close to polar day) ------------------------------------------------------
          ! case is now obsolete
          sd = COS(pi*(fract_localtime-rat1)/(fract_daylight/2._wp+rat4)) ! modulation factor for day-temperature such that ...
                                                                       ! ... temperature is at maximum at rat1 (2 pm) and at average
                                                                       ! ... at sunrise-rat4/2 and sunset+rat4/2
          tmp(ic,ib) = dtmp + dtran/2._wp*sd                             ! temperature modulated by temperature range
        ELSE ! ------------------------------- short day (close to polar night) ---------------------------------------------------
          tmp(ic,ib) = dtmp                                        ! temperature is set to mean temperture at all time steps
        ENDIF
      ENDDO
      !$ACC END PARALLEL
    ENDDO

    RETURN
  END SUBROUTINE INSTANTLY_FROM_DAILY_TEMP_2

  ! ======================================================================================================= !
  !>
  !> Calculate bottom pressure from air temperature and elevation
  !>
  !
  ! docu: bottom pressure()
  !
  ! To estimate bottom pressure a so called "polytrope atmosphere" is assumed for the troposhere. Here a fixed temperature
  ! gradient "gamma" is prescribed, i.e. air temperature at height z is set to
  !     (1)     T(z)=T0-gamma*z,
  ! where T0 is the temperature at the bottom. Assuming further hydrostatic conditions the vertical pressure gradient is given by
  !              wp
  !     (2)     ---- = -g*rho,
  !              dz
  ! where rho is the air density and g the acceleration due to gravity. Pressure and Temperature are related by the state
  ! equation of an ideal gas
  !     (3)      p = rho*R*T,
  ! where R is the specific gas constant and T the temperature. By eliminating rho from (3) and (2) and entering (1) for T one
  ! obtains
  !              wp        p*g
  !     (4)     ---- = - --------------
  !              dz      R(T0-gamma*z)
  ! By integration one obtains the final result
  !                              gamma                              g
  !     (5)      p = p0 * (1 - z*-----)^k , with the exponent k= ------- ,
  !                               T0                             R*gamma
  ! where p0 and T0 are the bottom pressure and temperature. It should be noted that the expression in brackets can get negative
  ! only for extreme values of T0 and z. Assuming e.g. T0 > -80 Celsius (approx. 200 Kelvin), the bracket stays positive for
  ! z<30km, which is true even at the Mount Everest. (Eq. (5) coincides with Eq. (31) from Knorr for z << 1!)
  !
  ! If temperature T(z) at height z is known instead of T0, using Eq. (1) Eq. (5) translates to:
  !                                  gamma
  !     (6)      p = p0 * (1 - z*------------)^k
  !                              T(z)+z*gamma
  !
  FUNCTION bottom_pressure(elevation,air_temperature)

  !$ACC ROUTINE SEQ

    REAL(wp), INTENT(IN)  :: elevation
    REAL(wp), INTENT(IN)  :: air_temperature
    REAL(wp)              :: bottom_pressure

    bottom_pressure = &
      & p0sl_bg * (1._wp - elevation * gamma / (tmelt + air_temperature + elevation * gamma)) ** (grav / rd / gamma)

  END FUNCTION bottom_pressure

  ! ======================================================================================================= !
  !>
  !> Split precipitation in snow and rain depending on mean air temperature
  !>
  !
  ! docu: snow_and_rain_from_precip()
  !
  ! Precipitation is divided into rain and snow according to a method by M.S. Wigmosta, L. Vail and D.P. Lettenmaier, "A
  ! distributed hydrology-vegetation model for complex terrain", Water Resources Research 30 (1994) 1665-1679: depending on
  ! average day temperature "T", the snow part of precipitation "P_sn" is computed from the total precipitation "P" as:
  !               / P               for T < -1.1 Celsius
  ! (1)   P_sn = |  P*(3.3 - T)/4.4 for -1.1 Celsius <= T <= 3.3 Celsius
  !               \ 0               for T > 3.3 Celsius
  ! Hence the rain part ist then
  ! (2)   P_rain = P-P_sn.
  !
  ELEMENTAL SUBROUTINE snow_and_rain_from_precip(rain,snow,precipitation,air_temp_daily)

    !$ACC ROUTINE SEQ

    REAL(wp),INTENT(IN)  :: precipitation   !< total precipitation (rain+snow) [kg/m^2/s]
    REAL(wp),INTENT(IN)  :: air_temp_daily  !< daily mean air_temperature at surface [Celsius]
    REAL(wp),INTENT(OUT) :: rain            !< rain [kg/m^2/s]
    REAL(wp),INTENT(OUT) :: snow            !< snow [kg/m^2/s]

    IF(air_temp_daily < -1.1_wp) THEN
       snow = precipitation
    ELSEIF(air_temp_daily > 3.3_wp) THEN
       snow = 0.0_wp
    ELSE
       snow = precipitation*(3.3_wp-air_temp_daily)/4.4_wp
    ENDIF
    ! set the rain part according to Eq. (2):
    rain = precipitation - snow

  END SUBROUTINE snow_and_rain_from_precip

  ! ======================================================================================================= !
  !>
  !> computes shortwave parts (PAR,NIR,fraction of diffuse PAR)
  !>
  !
  ! For timestep data the observed radiation is used directly;
  ! in daily version, the instantaneous flux at noon required to calculate irradiation at any
  ! time of the day is approximated by solving the integral over the positive sun angle, following Prentice et al. Ecol. Mod. 1993
  !
  ! ATTENTION: TR: there maybe considerable deviations in the daily sum of insolation  if daily data is used and
  ! if the orbit is different then todays orbit (i.e. ATTENTION: this may be incorrect for paleo situations)
  !
  ! Rest of the routine follows [1]. It computes PAR and NIR, together with the diffuse part of PAR, from "fpar", which is the
  ! quotient of the actual PAR flux "Rpar_act" (i.e. at cloudy sky) to the potential PAR flux "Rpar_pot" (i.e. at clear sky), both
  ! measured at the ground:
  !                       Rpar_act
  ! (1)           fpot := -------- ,
  !                       Rpar_pot
  ! Following [1] the various components of the total shortwave radiation are computed as follows: The total flux in the PAR band
  ! incident to the outer atmosphere "Rpar_top" is computed from
  ! (2)           Rpar_top = S0*cos(theta)*E,
  ! where S0=0.44*1360 W/m^2 =~= 600 W/m^2 is the PAR part of the solar constant (0.4-0.7 nm), "theta" is the zenith angle
  ! and
  ! (3)           E = (r0/r)^2 = 1.000110+0.034221 cos(psi)+0.001280*sin(psi)+0.000719*cos(2*psi)+0.000077*sin(2*psi)
  ! the inverse squared distance earth-sun "r" in units of the average distance "r0", where "psi" is the day-angle, i.e. the
  ! day of the year expressed as an angle (see e.g. [4], Eq. (3.3), [1] Eq. (24)).
  ! The direct part of potential PAR reaching ground "Rpar_pot_dir" is according to [2] approximately given by
  ! (4)          Rpar_pot_dir = Rpar_top*exp(-0.185*(p/p0)/cos(theta)),
  ! where "p" is the current and "p0" the sealevel pressure (see [1] Eqs. (27)-(29) and [2] Eqs. (1)-(3)), and the exponential
  ! term accounts for direct radiation reduced by Rayleigh scattering ( p/p0/cos(theta) is a measure for the mass of air along
  ! the beam ). According to [2] the diffuse part of potential PAR ("Rpar_pot_diff") reaching ground is approximately 40% of the
  ! difference between Rpar_top and the direct radiation (Equation (3) n [2] contains a printing error?):
  ! (5)           Rpar_pot_diff = 0.4*(Rpar_top - Rpar_pot_dir)
  ! so that potential PAR is
  ! (6)           Rpar_pot = Rpar_pot_dir + Rpar_pot_diff = Rpar_top*[0.4 + 0.6*exp(-0.185*(p/p0)/cos(theta))].
  ! By (1) the actual PAR then follows as
  ! (7)           Rpar_act = fpar * Rpar_pot.
  ! Next the near infrared (NIR) part of shortwave radiation has to be accounted for. Following [3] this is done by introducing
  ! a correction factor "F", such that the total downwards shortwave radiation "Rtot" reaching bottom is obtained from Rpar_act by
  ! (8)           Rtot = Rpar_act/F,
  ! where, using mfpar=1-fpar, the inverse factor is given by
  ! (9)           1/F = 1 + (1.185-0.437*mfpar-0.494*mfpar^2)* exp[(0.0305-0.208*mfpar+0.495*mfpar^2)/cos(theta)]
  ! Therefore the downward radiation reaching the surface in the near infrared band "Rnir" is
  ! (10)          Rnir = Rtot - Rpar_act
  ! Finally the fraction of direct radiation in PAR at the actual PAR
  !                           Rpar_act_dir
  ! (11)          fdirPar =   ------------
  !                             Rpar_act
  ! is estimated from the direct part in potential PAR
  ! (12)          Rpar_pot_dir = S0 * exp(-0.185*(p/p0)/cos(theta))
  ! total potential PAR "Rpar_pot" and "fpar" by (see Knorr Eq. (30))
  !                            /  0 for fpar < 0.2
  !                           /                   2/3
  !                          /      / 0.9 - fpar \      Rpar_pot_dir
  ! (13)          fdirPar = (  (1- (--------------)   ) ------------       for 0.2 <= fpar <= 0.9
  !                          \      \    0.7     /         Rpar_pot
  !                           \
  !                            \  1 for fpar > 0.9
  ! where the quotient showing up in (12) is
  !
  !               Rpar_pot_dir         exp(-0.185*(p/p0)/cos(theta))
  ! (14)          ------------ =   ---------------------------------------
  !                 Rpar_pot       0.4 + 0.6*exp(-0.185*(p/p0)/cos(theta))
  !
  ! It follows that the fraction of diffuse radiation "fdiffPar" in PAR is
  !
  ! (15)         fdiffPar = 1 - fdirPar.
  !
  ! REMARK: Eq. (13) is taken from [2] but contains an additional approximation: In [2] fdirPar does not depend on fpar, but on
  ! the transmission ratio for TOTAL SHORTWAVE flux, whereas fpar is defined as the transmission ratio for PAR only (see Eq. (1)).
  ! Hence, (13) includes the additional approximation of the shortwave transmission ratio by the PAR transmission ratio (fpar).
  !
  ! Literature:
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] A. Weiss & J.M. Norman, "Partitioning solar radiation into direct and diffuse visible and near-infrared components",
  ! Agricultural and Forest Meteorology, 34 (1985) 205-213.
  ! [3] R.T. Pinker & I. Laszlo, "Global distribution of photosynthetically active radiation as observed from satellites",
  ! J. Climate 5 (1992) 56-65.
  ! [4] G.W. Paltridge and C.M.R. Platt, "Radiative Processes in Meteorology and Climatology" (Elsevier, Amsterdam, 1976).
  !
  SUBROUTINE shortwave_from_direct_shortwave( &
    & model_id, &
    & day_of_year, &
    & cos_zenith, &
    & cos_lat, &
    & sin_lat, &
    & air_pressure, &
    & rad_sw_in, &
    & rad_UV, &
    & rad_PAR, &
    & rad_NIR, &
    & fract_PAR_diffuse, &
    & rad_sw, &
    & rad_sw_pot )

    INTEGER,  INTENT(in)  :: model_id           !< model ID
    INTEGER,  INTENT(in)  :: day_of_year        !< day in year (from [1,365])
    REAL(wp), INTENT(in), DIMENSION(:,:) :: &
      & cos_zenith,        & !< cosine of zenith angle
      & cos_lat,           & !< cosine of latitude
      & sin_lat,           & !< sine of latitude
      & air_pressure,      & !< air pressure at bottom (depending on local elevation) [N/m^2]
      & rad_sw_in            !< "Rtot": either timestep or daily average total solar shortwave radiation :
                             ! ... direct+diffuse, 280-3000 nm, (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(wp), INTENT(out), DIMENSION(:,:) :: &
      & rad_UV,            & !< "RUV": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                             ! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
      & rad_PAR,           & !< "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                             ! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
      & rad_NIR,           & !< "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                             ! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
      & fract_PAR_diffuse, & !< "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                             ! .. (values in [0,1])
      & rad_sw,            & !< ..
      & rad_sw_pot           !< ..

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: fsw           !! "fsw": fraction of potential radiation solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                              !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(wp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(wp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r
    REAL(wp) :: delta,K

    REAL(wp), PARAMETER  ::  fract_UV = 0.1_wp ! UV fraction of (UV + PAR)

    INTEGER :: i, j, nc, nblks

    REAL(wp), ALLOCATABLE :: u(:,:), v(:,:), hh(:,:), sinehh(:,:), rad_sw_in_corr(:,:)

    nc = SIZE(cos_zenith,1)
    nblks = SIZE(cos_zenith,2)

    ALLOCATE(u(nc, nblks), v(nc, nblks))
    ALLOCATE(hh(nc, nblks), sinehh(nc, nblks))
    ALLOCATE(rad_sw_in_corr(nc, nblks))
    !$ACC ENTER DATA CREATE(u, v, hh, sinehh, rad_sw_in_corr)

    day_angle = 2._wp*pi*(day_of_year+10._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations

    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)
    K=13750.98708_wp

    IF (      forcing_input(model_id)%shortwave%frequency /= TIMESTEP_ &
      & .AND. forcing_input(model_id)%shortwave%frequency /= SUBDAILY_) THEN
      ! today's sum of top of the atmosphere radiation
      delta=-23.4_wp*pi/180._wp*COS(day_angle)
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO j=1,nblks
        DO i=1,nc
          u(i,j) = sin_lat(i,j) * SIN(delta)
          v(i,j) = cos_lat(i,j) * COS(delta)
          IF (u(i,j).GE.v(i,j)) THEN
            hh(i,j) = pi
          ELSE IF (u(i,j).LE.(0._wp-v(i,j))) THEN
            hh(i,j) = 0.0_wp
          ELSE
            hh(i,j) = ACOS(-u(i,j) / v(i,j))
          END IF
          sinehh(i,j) = SIN(hh(i,j))

          !estimate of noon radiation giving rise to the recorded radiation sum
          rad_sw_in_corr(i,j) = 0._wp
          IF((u(i,j) * hh(i,j) + v(i,j) * sinehh(i,j)) .GT. 1e-5_wp)THEN
            rad_sw_in_corr(i,j) = rad_sw_in(i,j) / (2._wp * (u(i,j) * hh(i,j) + v(i,j) * sinehh(i,j)) * K) * 86400._wp
          END IF
          rad_sw(i,j) = rad_sw_in_corr(i,j) * cos_zenith(i,j)
        END DO
      END DO
      !$ACC END PARALLEL

    ELSE
      ! subdaily or time step forcing, use radiation directly
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO j=1,nblks
        DO i=1,nc
          rad_sw(i,j) = rad_sw_in(i,j)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO j=1,nblks
      DO i=1,nc
        IF(cos_zenith(i,j) > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
          Rpar_top   = solar_const*cos_zenith(i,j)*r_sun_earth_2       ! Eq. (2) from above
          hlp3_r = EXP(-0.185_wp * (air_pressure(i,j) / p0sl_bg) / cos_zenith(i,j)) ! in jsb3: p0sl_bg = p_sealevel
          hlp4_r = 0.4_wp + 0.6_wp * hlp3_r
          rad_sw_pot(i,j) = Rpar_top * hlp4_r
          fsw = MAX(MIN(rad_sw(i,j) / rad_sw_pot(i,j), 1._wp), 0._wp)
          hlp_r = 1._wp - fsw
          hlp2_r = hlp_r * hlp_r
          F_invers=1._wp + (1.185_wp - 0.437_wp * hlp_r - 0.494_wp * hlp2_r) * &
                 & EXP(MIN(2.534_wp, (0.0305_wp - 0.208_wp * hlp_r + 0.495_wp * hlp2_r) / cos_zenith(i,j))) ! Eq. (9) from above
          rad_UV(i,j) = rad_sw(i,j) * fract_UV / F_invers
          rad_PAR(i,j) = rad_sw(i,j) * (1._wp - fract_UV) / F_invers   ! This is "Rtot"; Eq. (8) from above
          rad_NIR(i,j) = rad_sw(i,j) - rad_PAR(i,j) - rad_UV(i,j) ! This is "Rnir"; Eq. (10) from above
          IF(fsw <= 0.2_wp) THEN
            fract_PAR_diffuse(i,j) = 1._wp                                                 ! Eq. (12)-(15) from above
          ELSE IF(fsw >= 0.9_wp) THEN
            fract_PAR_diffuse(i,j) = 0.08_wp                                               ! Eq. (12)-(15) from above
          ELSE
            fract_PAR_diffuse(i,j) = 1._wp - &
                                  & (1._wp- ((0.9_wp -fsw)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
          END IF
        ELSE ! zenith angle larger than 89 degrees - changed Aug. 2017 following T. Raddatz, Jan. 2015
          rad_UV(i,j) = 0._wp
          rad_PAR(i,j) = rad_sw(i,j) * 0.5_wp
          rad_NIR(i,j) = rad_sw(i,j) * 0.5_wp
          fract_PAR_diffuse(i,j) = 1._wp ! (could also be another value in this case)
          rad_sw_pot(i,j)  = rad_sw(i,j)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(u, v, hh, sinehh, rad_sw_in_corr)
    DEALLOCATE(u, v, hh, sinehh, rad_sw_in_corr)

  END SUBROUTINE shortwave_from_direct_shortwave

  ! ======================================================================================================= !
  !>
  !> Computes longwave from observed longwave and diural course associated with air temperature
  !>
  !
  ! docu longwave_from_daily_longwave()
  !
  ! Downward long wave radiation flux "R_d" [W/m^2] is according to [1],[2] computed by
  ! (1) R_d = r_cloud * epsA * sigma * T^4,
  ! where "T" is the air temperature [K] (in the bottom layer?), sigma = 5.6703e-8 W/(m^2 K^4) is the Stefan-Boltzmann constant,
  ! The mean daily longwave is known from observations, thus we are concerned only with the diurnal course associated with T
  ! [1] W. Knorr, "Satellite remote sensing and modelling of the global CO2 exchange of land vegetation", Examensarbeit 49,
  ! (Max Planck Institut fuer Meteorologie, Hamburg, 1998).
  ! [2] W. Brutsaert, "Evaporation into the Atmosphere", (Reidel, Dordrecht, 1982), pp. 299.
  !
  ELEMENTAL FUNCTION longwave_from_daily_longwave(longwave_down_daily, Tmin, Tmax, air_temp)

    !$ACC ROUTINE SEQ

    REAL(wp)            :: longwave_from_daily_longwave !! hourly longwave downward radiation [W/m^2]
    REAL(wp),INTENT(in) :: longwave_down_daily          !! daily average downward longwave radiation [W/m^2]
    REAL(wp),INTENT(in) :: Tmin                         !! daily mimimum temperature [Celsius]
    REAL(wp),INTENT(in) :: Tmax                         !! daily maximum temperature [Celsius]
    REAL(wp),INTENT(in) :: air_temp                     !! air temperature at bottom [Celsius]

    REAL(wp) :: dtmp    ! daily average temperature
    REAL(wp) :: hlp_r   ! helper

    hlp_r = air_temp + tmelt                                ! Celsius --> Kelvin
    dtmp = 0.5_wp * (Tmax + Tmin) + tmelt
    longwave_from_daily_longwave = longwave_down_daily * (hlp_r**4) / (dtmp**4)      ! lw-radiation flux Eq. (1)

  END FUNCTION longwave_from_daily_longwave

  ! ======================================================================================================= !
  ! Procedures that exist in ECHAM but not in ICON have to be defined here
  ! ======================================================================================================= !


  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  FUNCTION IMerge_HMS2Sec (kh, km, ks) RESULT (ix)
    !+
    !
    ! IMerge_HMS2Sec [function, integer]
    !    merge hour, minute and second to internal time
    !    (
    !    hour   [integer] input (hour of day, 0...23)
    !    minute [integer] input (minute of hour, 0...59)
    !    second [integer] input (second of minute, 0...59)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: kh, km, ks
    INTEGER             :: ix

    ix = kh*3600 + km*60 + ks

  END FUNCTION IMerge_HMS2Sec

  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  !> Do not put this in mo_jsb_time_iface. This function is only necessary in mo_jsb4_forcing.
  !>
  INTEGER FUNCTION day_in_year(date)

    TYPE(t_datetime), POINTER, INTENT(in)  :: date

    day_in_year = INT(get_year_day(date))

  END FUNCTION day_in_year

#endif
END MODULE mo_jsb4_forcing
