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
MODULE mo_jsb4_forcing_echam
#ifndef __NO_JSBACH__

  USE mo_time_control,    ONLY: lresume, lstart, get_year_day
  USE mo_time_conversion, ONLY: IMerge_HMS2Sec, day_in_year, time_days
  USE mo_jsb_time,        ONLY: get_date_components, t_datetime, is_newyear, is_newday, &
    &                           get_year_length, get_year_at_experiment_start

  ! Apparently in echam and icon
  USE mo_math_constants,         ONLY: pi
  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util_string,            ONLY: tolower, int2string, real2string
  USE mo_physical_constants,     ONLY: p0sl_bg, vtmpc2, rdv
  USE mo_jsb_physical_constants, ONLY: grav, & ! at the end grav comes from mo_physical_constants, however like this it is consistent with JSBACH
    & tmelt, rd, amco2, amd, cpd, rvd1, cpvd1, von_karman, rd_o_cpd
  USE mo_jsb_convect_tables_iface, ONLY: tlucua, jptlucu1, jptlucu2
  USE mo_orbit_solar,            ONLY: inquire_declination

  ! jsb4
  USE mo_jsb_class,              ONLY: jsbach
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
  USE mo_jsb_parallel,           ONLY: my_process_is_stdio, my_process_is_mpi_parallel, p_io, mpi_comm, p_bcast
  USE mo_jsb_io_netcdf,          ONLY: t_input_file, jsb_netcdf_open_input
  USE mo_jsb_io_netcdf_iface,    ONLY:                 &
    & nf, nf_get_vara_double,                          & ! former nf_check
    & nf_inq_dimid, nf_inq_dimlen, nf_inq_varid,       & ! from netcdf.inc
    & NF_MAX_NAME, nf_get_att_text, nf_get_att_double, & ! from netcdf.inc
    & NF_NOERR                                           ! from netcdf.inc

  USE mo_jsb_class,              ONLY: get_model

  dsl4jsb_Use_processes HYDRO_, A2L_
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(A2L_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: get_interface_variables_from_external_forcing, init_forcing, finalize_external_forcing, &
    &       sat_specific_humidity, forcing_options

  REAL(wp), SAVE :: delta_time, time_step_len
  INTEGER,  SAVE :: time_steps_per_day

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

  ! expected units in forcing files (units of precipitation and co2 are specified in the namelist parameters)
  CHARACTER(len=*), PARAMETER :: exp_unit_temperature = 'degC'
  CHARACTER(len=*), PARAMETER :: exp_unit_radiation = 'W/m2'
  CHARACTER(len=*), PARAMETER :: exp_unit_qair = 'g/g'
  CHARACTER(len=*), PARAMETER :: exp_unit_wspeed = 'm/s'
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

  ! --- Structures for handling input data
  ! options relating forcing read from namelist
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
    ! quincy @TODO  modify the use of the date-information to how this is done with jsbach4
    INTEGER  :: forcing_year = 0
    INTEGER  :: forcing_doy  = 0
    REAL(wp) :: forcing_hour = 0.0_wp
    ! quincy - end
  END TYPE forcing_options_type
  TYPE(forcing_options_type), SAVE :: forcing_options

  ! information for single forcing variable and the file it is read from
  TYPE input_variable_type
    CHARACTER(NF_MAX_NAME)   :: variable_name=""       !! Name of the variable (name in netcdf-data file)
    CHARACTER(NF_MAX_NAME)   :: file_name_prefix = ""  !! Prefix of files in which the variable is located
    CHARACTER(NF_MAX_NAME)   :: expectedUnit = ""      !! Unit for the variable expected in the file
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

  ! All theoretically possible input data, i.e. some are not used. (Input data is controlled by jsb_forcing_nml)
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
  END TYPE forcing_input_type
  TYPE(forcing_input_type),SAVE :: forcing_input

CONTAINS

  !>
  !! @brief provides the variables required for jsb4 interface calls in this timestep based on (namelist defined) external forcing
  !!
  !! @par History
  !! Used parts of jsbach3 mo_jsbalone_forcing: read_forcing / update_forcing and mo_jsbalone: update_driving  by Julia Nabel
  !!
  SUBROUTINE get_interface_variables_from_external_forcing(model_id, ics, ice, icz, nblks, &
    & current_datetime, next_datetime, sinlat, coslat, lon, evapo_act2pot_proc,            &
    & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, &
    & cos_zenith_angle, CO2_air,                                                           &
    & t_air, q_air, rain, snow, wind_air, wind_10m, lw_srf_down, swvis_srf_down,           &
    & swnir_srf_down, swpar_srf_down, fract_par_diffuse, press_srf,                        &
    & drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch,                                   &
    & t_srf_lake, fact_q_air_lake, fact_qsat_srf_lake, rough_h_srf_lake, rough_m_srf_lake, &
    & drag_wtr, t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr,                        &
    & drag_ice, t_acoef_ice, t_bcoef_ice, q_acoef_ice, q_bcoef_ice)

    INTEGER,  INTENT(IN)  :: model_id, ics, ice, icz, nblks
    TYPE(t_datetime),  INTENT(IN), POINTER :: current_datetime
    TYPE(t_datetime),  INTENT(IN), POINTER :: next_datetime
    REAL(wp), INTENT(IN)  :: sinlat(:,:)
    REAL(wp), INTENT(IN)  :: coslat(:,:)
    REAL(wp), INTENT(IN)  :: lon(:,:)
    REAL(wp), INTENT(IN)  :: evapo_act2pot_proc(:,:)
    REAL(wp), INTENT(IN)  :: t_srf_proc(:,:)          ! In setups with lakes refering to the land tile
    REAL(wp), INTENT(IN)  :: fact_q_air_proc(:,:)     !       ""
    REAL(wp), INTENT(IN)  :: fact_qsat_srf_proc(:,:)  !       ""
    REAL(wp), INTENT(IN)  :: rough_h_srf_proc(:,:)    !       ""
    REAL(wp), INTENT(IN)  :: rough_m_srf_proc(:,:)    !       ""
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
    REAL(wp), INTENT(OUT) :: drag_srf(:,:)
    REAL(wp), INTENT(OUT) :: t_acoef(:,:)
    REAL(wp), INTENT(OUT) :: t_bcoef(:,:)
    REAL(wp), INTENT(OUT) :: q_acoef(:,:)
    REAL(wp), INTENT(OUT) :: q_bcoef(:,:)
    REAL(wp), INTENT(OUT) :: pch(:,:)

    ! Optional input for lakes
    REAL(wp), OPTIONAL, INTENT(IN)  :: t_srf_lake(:,:)
    REAL(wp), OPTIONAL, INTENT(IN)  :: fact_q_air_lake(:,:)
    REAL(wp), OPTIONAL, INTENT(IN)  :: fact_qsat_srf_lake(:,:)
    REAL(wp), OPTIONAL, INTENT(IN)  :: rough_h_srf_lake(:,:)
    REAL(wp), OPTIONAL, INTENT(IN)  :: rough_m_srf_lake(:,:)
    ! Optional output for lakes
    REAL(wp), OPTIONAL, INTENT(OUT) :: drag_wtr(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: t_acoef_wtr(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: t_bcoef_wtr(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: q_acoef_wtr(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: q_bcoef_wtr(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: drag_ice(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: t_acoef_ice(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: t_bcoef_ice(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: q_acoef_ice(:,:)
    REAL(wp), OPTIONAL, INTENT(OUT) :: q_bcoef_ice(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_interface_variables_from_external_forcing'

    TYPE(t_jsb_model),  POINTER  :: model
    CLASS(t_jsb_tile_abstract), POINTER :: tile

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Real2D_onDomain :: elevation
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Real2D_onDomain :: DEBUG_VAR

    REAL(wp) :: hlp_r
    INTEGER  :: ic, iblk
    INTEGER  :: day_of_year

    REAL(wp), ALLOCATABLE :: t_air_C(:,:), tmax_data(:,:), tmin_data(:,:), rad_uv_down(:,:), &
      &                      sat_spec_humidity(:,:), rel_humidity(:,:), pch_lake(:,:)
    REAL(wp), ALLOCATABLE :: temporary_data1(:,:)

    ! TODO: help1 is rad_sw_down, help2 is rad_sw_down_pot, but these are (currently?) not required - remove?
    REAL(wp), ALLOCATABLE :: help1(:,:), help2(:,:)

    !---------------------------------------------------------------------------------

    model => get_model(model_id)
    CALL model%Get_top_tile(tile)

    IF (iblk == nblks) THEN
      ic = icz
    ELSE
      ic = ice
    END IF

    ! Initialization needed for grid cells beyond icz
    CO2_air(:,:)           = 0._wp
    t_air(:,:)             = 0._wp
    q_air(:,:)             = 0._wp
    rain(:,:)              = 0._wp
    snow(:,:)              = 0._wp
    wind_air(:,:)          = 0._wp
    wind_10m(:,:)          = 0._wp
    lw_srf_down(:,:)       = 0._wp
    swvis_srf_down(:,:)    = 0._wp
    swnir_srf_down(:,:)    = 0._wp
    swpar_srf_down(:,:)    = 0._wp
    fract_par_diffuse(:,:) = 0._wp
    press_srf(:,:)         = 0._wp

    ALLOCATE(help1(ice, nblks))
    ALLOCATE(help2(ice, nblks))

    ALLOCATE(temporary_data1(ice, nblks))

    ! --- temperature
    ! for temperature the variables used depend on the forcing frequency: tmin&tmax for daily data; t_air_C for timestep
    ALLOCATE(t_air_C(ice, nblks))
    IF(forcing_options%air_temperature_as_timestep) THEN
      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, forcing_input%air_temp, t_air_C)
    ELSE
      ALLOCATE(tmin_data(ice, nblks))
      ALLOCATE(tmax_data(ice, nblks))

      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
        & forcing_input%air_temp_min, tmin_data)
      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
        & forcing_input%air_temp_max, tmax_data)
      ! estimate temperature at time step from daily temperature maximum and temperature minimum
      CALL instantly_from_daily_temp_2(nblks, ice, icz, lon, sinlat, coslat, tmin_data, tmax_data, next_datetime, t_air_C)
    ENDIF
    ! convert degC to Kelvin
    t_air(:,:) = t_air_C(:,:) + tmelt

    ! --- precipitation
    CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
      &                            forcing_input%precipitation, temporary_data1)
    ! convert if required
    IF(forcing_options%precip_in_mm_per_day) temporary_data1 = temporary_data1 * conv_prec_day
    ! compute rain and snow
    IF(forcing_options%air_temperature_as_timestep) THEN
      IF(forcing_input%precipitation%frequency == DAILY_) THEN
        ! In case of daily precipitation the minimum and maximum read for air_temp are used
        CALL snow_and_rain_from_precip(rain(:,:), snow(:,:), temporary_data1(:,:), &
          & 0.5_wp*(MAXVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3) + MINVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3)))
      ELSE
        ! else currentvalue
        CALL snow_and_rain_from_precip(rain(:,:), snow(:,:), temporary_data1(:,:), t_air_C(:,:))
      ENDIF
    ELSE
      CALL snow_and_rain_from_precip(rain(:,:), snow(:,:), temporary_data1(:,:), &
        & 0.5_wp*(tmax_data(:,:) + tmin_data(:,:)))
    ENDIF

    ! get elevation
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onDomain(HYDRO_, elevation)  ! in

    ! --- compute bottom pressure
    DO iblk = 1, nblks
      IF (iblk == nblks) THEN
        ic = icz
      ELSE
        ic = ice
      END IF

      press_srf(1:ic,iblk) = bottom_pressure(elevation(1:ic,iblk), t_air_C(1:ic,iblk))
    ENDDO

    ! --- compute specific humidity
    ALLOCATE(sat_spec_humidity(ice, nblks))
    sat_spec_humidity(:,:) = 0._wp
    SELECT CASE(forcing_options%type_of_qair_forcing)
    CASE(NONE_)
      ! ratio of actual to potential evapotranspiration from the previous day
      IF (forcing_options%air_temperature_as_timestep) THEN
        temporary_data1 = vapour_pressure_from_evapor(nblks, ice, icz, evapo_act2pot_proc(:,:), &
          & MINVAL(forcing_input%air_temp%data_read(:,:,:), DIM=3), t_air_C(:,:))
      ELSE
        temporary_data1 = vapour_pressure_from_evapor(nblks, ice, icz, evapo_act2pot_proc(:,:), tmin_data(:,:), t_air_C(:,:))
      ENDIF
      q_air(:,:) = fun_specific_humidity(temporary_data1(:,:), press_srf(:,:))

    CASE(RH_)
      ALLOCATE(rel_humidity(ice, nblks))
      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
        &                            forcing_input%rel_humidity, rel_humidity)
      ! check data
      hlp_r = MAXVAL(MAXVAL(rel_humidity,DIM=2))
      IF(hlp_r < 1._wp) THEN
        WRITE (message_text,*) 'Maximum value of relative humidity is ', hlp_r
        CALL message(TRIM(routine), TRIM(message_text))
        CALL message(TRIM(routine),'relative humidity data probably represent fractions instead of percentages.')
      ENDIF
      ! Convert percentage data to fractional data:
      rel_humidity = rel_humidity/100._wp
      ! Get saturation spec. humidity from air temperature
      DO iblk = 1, nblks
        IF (iblk == nblks) THEN
          ic = icz
        ELSE
          ic = ice
        END IF
        CALL sat_specific_humidity(ice, t_air(1:ic,iblk), press_srf(1:ic,iblk), sat_spec_humidity(1:ic,iblk))
      ENDDO
      ! Get specific humidity
      IF (forcing_input%rel_humidity%frequency == TIMESTEP_ .OR. forcing_input%rel_humidity%frequency == SUBDAILY_) THEN
        q_air(:,:) = rel_humidity(:,:) * sat_spec_humidity(:,:) !TODO: Why *?
      ELSE
        IF (forcing_options%air_temperature_as_timestep) THEN
          DO iblk = 1, nblks
            IF (iblk == nblks) THEN
              ic = icz
            ELSE
              ic = ice
            END IF
            CALL spec_humidity_from_rel_humidity(ic, sinlat(1:ic,iblk), coslat(1:ic,iblk), &
              & rel_humidity(1:ic,iblk), MINVAL(forcing_input%air_temp%data_read(1:ic,iblk,:),DIM=2), &
              & MAXVAL(forcing_input%air_temp%data_read(1:ic,iblk,:),DIM=2), elevation(1:ic,iblk), q_air(1:ic,iblk))
          ENDDO
        ELSE
          DO iblk = 1, nblks
            IF (iblk == nblks) THEN
              ic = icz
            ELSE
              ic = ice
            END IF
            CALL spec_humidity_from_rel_humidity( &
              & ic, sinlat(1:ic,iblk), coslat(1:ic,iblk), rel_humidity(1:ic,iblk), &
              & tmin_data(1:ic,iblk), tmax_data(1:ic,iblk), elevation(1:ic,iblk), q_air(1:ic,iblk))
          ENDDO
        ENDIF
        ! constrain daily specific humidity to be lower than or equal to the saturation spec. humidity
        q_air(:,:) = MIN(q_air, sat_spec_humidity)
      ENDIF
      DEALLOCATE(rel_humidity)
    CASE(QAIR_)
      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, forcing_input%qair, q_air(:,:))
      ! constrain specific humidity to be lower than or equal to the saturation spec. humidity
      DO iblk = 1, nblks
        IF (iblk == nblks) THEN
          ic = icz
        ELSE
          ic = ice
        END IF
        CALL sat_specific_humidity(ic, t_air(1:ic,iblk), press_srf(1:ic,iblk), sat_spec_humidity(1:ic,iblk))
      ENDDO
      q_air(:,:) = MIN(q_air(:,:),sat_spec_humidity(:,:))
    END SELECT
    DEALLOCATE(sat_spec_humidity)

    ! --- compute shortwave radiation
    day_of_year=day_in_year(next_datetime)
    ALLOCATE(rad_uv_down(ice, nblks))
    CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
      & forcing_input%shortwave, temporary_data1)
    CALL shortwave_from_direct_shortwave( &
      ! in
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

    swvis_srf_down(:,:) = rad_uv_down(:,:) + swpar_srf_down(:,:)
    DEALLOCATE(rad_uv_down)

    ! --- compute longwave radiation
    IF (forcing_input%longwave%frequency == TIMESTEP_ .OR. forcing_input%longwave%frequency == SUBDAILY_) THEN
      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
        & forcing_input%longwave, lw_srf_down)
    ELSE
      ! jsbach3 driver: diurnal longwave downward radiation is generated such that the average radiation matches the
      ! observations but the diurnal cycle follows temperature to help calculation of a correct
      ! net radiation flux. Note that reduces biases resulting from assuming a constant flux, but
      ! this introduces uncertainty due to biases in the assumed dial course of temperature
      CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
        &                            forcing_input%longwave, temporary_data1)

      IF (forcing_options%air_temperature_as_timestep) THEN
        lw_srf_down = longwave_from_daily_longwave(temporary_data1(:,:), MINVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3), &
          & MAXVAL(forcing_input%air_temp%data_read(:,:,:),DIM=3), t_air_C(:,:))
      ELSE
        lw_srf_down = longwave_from_daily_longwave(temporary_data1(:,:), tmin_data(:,:), tmax_data(:,:), t_air_C(:,:))
      ENDIF
    ENDIF

    ! --- wind
    CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, forcing_input%wind_speed, wind_10m)
    wind_air = wind_10m

    ! --- CO2 concentration
    CALL get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, forcing_input%CO2_concentr, CO2_air)
    ! and convert to correct unit
    CO2_air = CO2_air * forcing_options%conv_CO2_2_MassRatio

    ! --- compute remaining required input for jsb4 interface
    CALL get_remaining_input_for_jsb4_interface(nblks, ice, icz, t_air, press_srf, q_air, wind_air, &
      & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, &
      & drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)

    IF (model%config%use_lakes) THEN
      ALLOCATE(pch_lake(ice, nblks))
      ! Note: pch_lake is only defined, as pch should not be overwritten. As pch is only used for
      !       calculations on vegetated tiles, pch_lake is not needed.
      CALL get_remaining_input_for_jsb4_interface(nblks, ice, icz, t_air, press_srf, q_air, wind_air, &
        & t_srf_lake, fact_q_air_lake, fact_qsat_srf_lake, rough_h_srf_lake, rough_m_srf_lake,        &
        & drag_wtr, t_acoef_wtr, t_bcoef_wtr, q_acoef_wtr, q_bcoef_wtr, pch_lake)
      ! As input variables for the lake tile refer to the frozen and unfrozen fractions, we cannot
      ! provide separate exchange coefficients for water and ice. Both are set to identical values.
      drag_ice    = drag_wtr
      t_acoef_ice = t_acoef_wtr
      t_bcoef_ice = t_bcoef_wtr
      q_acoef_ice = q_acoef_wtr
      q_bcoef_ice = q_bcoef_wtr
      DEALLOCATE (pch_lake)
    END IF

    ! ---- clean up
    DEALLOCATE(t_air_C)
    IF(.NOT. forcing_options%air_temperature_as_timestep) THEN
      DEALLOCATE(tmin_data)
      DEALLOCATE(tmax_data)
    ENDIF

    DEALLOCATE(temporary_data1)
    DEALLOCATE(help1, help2)

  END SUBROUTINE get_interface_variables_from_external_forcing


  !>
  !! @brief init reading forcing data (read namelist, set options...)
  !!
  !! @par History
  !! Copied first version from jsbach3 mo_jsbalone_forcing     by Julia Nabel
  !! -- after r304: removed all currently not required options
  !!
  SUBROUTINE init_forcing(model_id, dtime, steplen)

    USE mo_jsb_namelist_iface, ONLY: POSITIONED, open_nml, position_nml
    USE mo_jsb_io_netcdf_iface, ONLY: t_input_file, netcdf_open_input

    INTEGER,  INTENT(in) :: model_id
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
    CHARACTER(NF_MAX_NAME) :: forcing_temp_file_prefix         !! Prefix for files with temperature forcing data
    CHARACTER(len=128)     :: forcing_temp_frequ               !! Frequency of temperature forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_temp_const_tmin          !! Constant value for minimum daily temperature
    REAL(wp)               :: forcing_temp_const_tmax          !! Constant value for maximum daily temperature
    REAL(wp)               :: forcing_temp_ocean               !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_temp_unit                !! Expected unit of temp input

    CHARACTER(NF_MAX_NAME) :: forcing_precip_file_prefix       !! Prefix for files with precipitation forcing data (unit: mm/day
    !! or kg/m^2/s, depending on 'forcing_precip_in_mm_per_day')
    CHARACTER(len=128)     :: forcing_precip_frequ             !! Frequency of precipitation forcing (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_precip_const_precip      !! Constant value for precipitation
    REAL(wp)               :: forcing_precip_ocean             !! Constant value for ocean cells
    LOGICAL                :: forcing_precip_in_mm_per_day     !! Precipitation input in mm/day (true) or kg/m2/s (false)
    CHARACTER(len=128)     :: forcing_precip_unit              !! Expected unit of precip input

    CHARACTER(NF_MAX_NAME) :: forcing_sw_file_prefix           !! Prefix for files with shortwave radiation data
    CHARACTER(len=128)     :: forcing_sw_frequ                 !! Frequency of shortwave forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_sw_const_shortwave       !! Constant value for downward shortwave radiation
    REAL(wp)               :: forcing_sw_ocean                 !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_sw_unit                  !! Expected unit of shortwave input (default: [W/m^2])

    CHARACTER(NF_MAX_NAME) :: forcing_lw_file_prefix           !! Prefix for files with longwave radiation data
    CHARACTER(len=128)     :: forcing_lw_frequ                 !! Frequency of longwave forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_lw_const_longwave        !! Constant value for downward longwave radiation
    REAL(wp)               :: forcing_lw_ocean                 !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_lw_unit                  !! Expected unit of longwave input (default: [W/m^2])

    CHARACTER(NF_MAX_NAME) :: forcing_co2_file_prefix          !! Prefix for files with CO2 data (unit: mol(CO2)/mol(air)
    !! or kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit'
    CHARACTER(len=128)     :: forcing_co2_frequ                !! Frequency of CO2 forcing data (DAILY/CONST/TIMESTEP/GHG_SCENARIO)
    REAL(wp)               :: forcing_co2_const_co2            !! Constant value for CO2-concentration (unit: mol(CO2)/mol(air) or
    !! kg(CO2)/kg(air) or ppmv, depending on 'forcing_co2_unit')
    REAL(wp)               :: forcing_co2_ocean                !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_co2_unit                 !! Unit of CO2 input (PPMV/MOL_PER_MOL/KG_PER_KG/'')

    CHARACTER(NF_MAX_NAME) :: forcing_wind_file_prefix         !! Prefix for files with windspeed data
    CHARACTER(len=128)     :: forcing_wind_frequ               !! Frequency of windspeed forcing data (DAILY/CONST/TIMESTEP)
    REAL(wp)               :: forcing_wind_const_wspeed        !! Constant value for windspeed [m/s]
    REAL(wp)               :: forcing_wind_ocean               !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_wind_unit                !! Expected unit of wind input (default: [m/s])

    CHARACTER(NF_MAX_NAME) :: forcing_qair_file_prefix         !! Prefix for files with qair data
    CHARACTER(len=128)     :: forcing_qair_frequ               !! Frequency of qair forcing data (DAILY/CONST)
    CHARACTER(len=128)     :: forcing_qair_type                !! Type of input data used for forcing by atm. humidity
    REAL(wp)               :: forcing_qair_const_rh            !! Constant value for relative humidity
    REAL(wp)               :: forcing_qair_ocean               !! Constant value for ocean cells
    CHARACTER(len=128)     :: forcing_qair_unit                !! Expected unit of qair input

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
      & forcing_qair_unit

    TYPE(t_jsb_model), POINTER :: model
    INTEGER:: nml_handler, read_status, f_unit

    !------------------------------------------------------------------------------------

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

      ! ----- read the namelist
      model => get_model(model_id)

      ! nml_handler = open_nml ('NAMELIST_JS4-standalone_master')
      nml_handler = open_nml (TRIM(model%namelist_filename))
      f_unit = position_nml ('jsb_forcing_nml', nml_handler, status=read_status)
      IF (read_status .EQ. POSITIONED) READ(f_unit, jsb_forcing_nml)

      forcing_options%cyclic_nyears           = cyclic_nyears
      IF (cyclic_start_year == 9999) cyclic_start_year = get_year_at_experiment_start()
      forcing_options%cyclic_start_year       = cyclic_start_year
      forcing_options%forcing_synchron_factor = forcing_synchron_factor
      forcing_options%forcing_steps_per_day   = forcing_steps_per_day

      IF ((forcing_steps_per_day > time_steps_per_day)) CALL finish(TRIM(routine), &
        & 'Violation of assertion: Forcing time resolution higher than model time step.')
      IF ((forcing_steps_per_day < 1)) CALL finish(TRIM(routine), &
        & 'Violation of assertion: Forcing time step longer than 1 day is not yet implemented.')

      IF (cyclic_nyears > 0) THEN
        CALL message(routine, 'Using cyclic forcing from ' // TRIM(int2string(cyclic_nyears)) // ' years ' // &
          &                   'starting in ' // TRIM(int2string(cyclic_start_year)))
      END IF

        ! ----- derive options - forcing heights: as in jsb3 it is expected that "HeightHumidity and HeightTemperature have to be equal"
      forcing_options%heightWind = forcing_height_wind
      forcing_options%heightHumidity = forcing_height_humidity
      ! --- Temperature timestep (if daily two variables - max and min - need to be read)
      forcing_options%air_temperature_as_timestep &
        & = (frequency_key_to_constant_val('temperature', forcing_temp_frequ) == TIMESTEP_ .OR. &
        &    frequency_key_to_constant_val('temperature', forcing_temp_frequ) == SUBDAILY_)
      ! --- Boolean indicating if forcing is only specified for land cells
      forcing_options%forcing_set_ocean_to_constants = forcing_set_ocean_to_constants
      IF(forcing_options%forcing_set_ocean_to_constants) THEN
        CALL message(TRIM(routine),'forcing_set_ocean_to_constants = true: ocean cells will be set to the specified constants. '&
          & //'On ERROR: make sure that your jsbach4 land sea mask matches the land sea mask of your forcing!')
      ELSE
        CALL message(TRIM(routine),'forcing_set_ocean_to_constants = false: input is expected for all cells for all variables. '&
          & //'On ERROR: make sure that your forcing is actually defined for all cells for all variables!')
      ENDIF
      forcing_options%forcing_set_miss_to_constants = forcing_set_miss_to_constants
      IF(forcing_options%forcing_set_miss_to_constants) THEN
        CALL message(TRIM(routine),&
          & 'forcing_set_miss_to_constants == true: Missing values over land will be replace with ocean constants')
      ENDIF
      ! --- Precipitation unit
      forcing_options%precip_in_mm_per_day = forcing_precip_in_mm_per_day
      IF(forcing_options%precip_in_mm_per_day) THEN
        CALL message(TRIM(routine),'Precipitation forcing in mm/day')
      ELSE
        CALL message(TRIM(routine),'Precipitation forcing in kg/m2/s')
      ENDIF
      ! --- QAIR type
      SELECT CASE(TRIM(tolower(forcing_qair_type)))
      CASE("rh")
        forcing_options%type_of_qair_forcing = RH_
      CASE("qair")
        forcing_options%type_of_qair_forcing = QAIR_
      CASE("none")
        forcing_options%type_of_qair_forcing = NONE_
        CALL message(TRIM(routine),'WARNING: QAIR will be diagnosed from Air temperature (type_of_qair_forcing = none)')
      CASE default
        CALL finish(TRIM(routine),'Please specifiy either rh, qair or none for "FORCING_QAIR_TYPE" in jsb_forcing_nml')
      END SELECT
      CALL message(TRIM(routine),"Type of air humidity forcing: "//TRIM(forcing_qair_type))
      ! --- CO2 unit
      SELECT CASE(TRIM(tolower(forcing_co2_unit)))
      CASE("ppmv","1.e-6")
        CALL message(TRIM(routine),'CO2 forcing in ppmv')
        forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg * 0.000001_wp
      CASE("mol_per_mol", "mol/mol")
        CALL message(TRIM(routine),'CO2 forcing in mol(CO2)/mol(Dry Air)')
        forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg
      CASE("kg_per_kg", "kg/kg")
        CALL message(TRIM(routine),'CO2 forcing in kg(CO2)/kg(Dry Air)')
        forcing_options%conv_CO2_2_MassRatio = 1.0_wp
      CASE default
        CALL message(TRIM(routine),'WARNING: "FORCING_CO2_UNIT" missing in namelist jsb_forcing_nml.')
        CALL message(TRIM(routine),'Assumption: unit is mol(CO2)/mol(DryAir). Verify your settings!')
        forcing_options%conv_CO2_2_MassRatio = molarMassCO2_kg / molarMassDryAir_kg
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
      CALL p_bcast(forcing_options%cyclic_nyears, p_io, mpi_comm)
      CALL p_bcast(forcing_options%cyclic_start_year, p_io, mpi_comm)
      CALL p_bcast(forcing_options%forcing_synchron_factor, p_io, mpi_comm)
      CALL p_bcast(forcing_options%forcing_steps_per_day, p_io, mpi_comm)
      CALL p_bcast(forcing_options%forcing_set_ocean_to_constants, p_io, mpi_comm)
      CALL p_bcast(forcing_options%forcing_set_miss_to_constants, p_io, mpi_comm)
      CALL p_bcast(forcing_options%air_temperature_as_timestep, p_io, mpi_comm)
      CALL p_bcast(forcing_options%precip_in_mm_per_day, p_io, mpi_comm)
      CALL p_bcast(forcing_options%type_of_qair_forcing, p_io, mpi_comm)
      CALL p_bcast(forcing_options%conv_CO2_2_MassRatio, p_io, mpi_comm)
      CALL p_bcast(forcing_options%heightWind, p_io, mpi_comm)
      CALL p_bcast(forcing_options%heightHumidity, p_io, mpi_comm)

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
    ENDIF

    ! ----- init forcing_input
    ! --- temperature
    IF(forcing_options%air_temperature_as_timestep) THEN
      CALL init_forcing_input_variable(air_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, forcing_temp_unit, forcing_temp_const_tmax, forcing_input%air_temp)
    ELSE
      CALL init_forcing_input_variable(min_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, forcing_temp_unit, forcing_temp_const_tmin, forcing_input%air_temp_min)
      CALL init_forcing_input_variable(max_temp_varname, forcing_temp_frequ, forcing_temp_ocean, &
        & forcing_temp_file_prefix, forcing_temp_unit, forcing_temp_const_tmax, forcing_input%air_temp_max)
    ENDIF
    ! --- precipitation
    CALL init_forcing_input_variable(precipitation_varname, forcing_precip_frequ, forcing_precip_ocean, &
      & forcing_precip_file_prefix, forcing_precip_unit, forcing_precip_const_precip, forcing_input%precipitation)
    ! --- atmospheric humidity
    SELECT CASE(forcing_options%type_of_qair_forcing)
    CASE(RH_)
      CALL init_forcing_input_variable(rel_humidity_varname, forcing_qair_frequ, forcing_qair_ocean, &
        forcing_qair_file_prefix, forcing_qair_unit ,forcing_qair_const_rh, forcing_input%rel_humidity)
    CASE(QAIR_)
      CALL init_forcing_input_variable(qair_varname, forcing_qair_frequ, forcing_qair_ocean, &
        forcing_qair_file_prefix, forcing_qair_unit ,forcing_qair_const_rh, forcing_input%qair)
    CASE(NONE_)
      ! no extra forcing needs to be initialised
    CASE default
      CALL finish(TRIM(routine),'Please specifiy either qair, rh or none for "forcing_qair_type" in jsb_forcing_nml')
    END SELECT
    ! --- shortwave radiation
    CALL init_forcing_input_variable(shortwave_varname, forcing_sw_frequ, forcing_sw_ocean, &
      & forcing_sw_file_prefix, forcing_sw_unit, forcing_sw_const_shortwave, forcing_input%shortwave)
    ! --- longwave radiation
    CALL init_forcing_input_variable(longwave_varname, forcing_lw_frequ, forcing_lw_ocean, &
      & forcing_lw_file_prefix, forcing_lw_unit, forcing_lw_const_longwave, forcing_input%longwave)
    ! --- CO2
    CALL init_forcing_input_variable(CO2_varname, forcing_co2_frequ, forcing_co2_ocean, &
      & forcing_co2_file_prefix, forcing_co2_unit, forcing_co2_const_co2, forcing_input%CO2_concentr)
    ! --- Wind
    CALL init_forcing_input_variable(wspeed_varname, forcing_wind_frequ, forcing_wind_ocean, &
      & forcing_wind_file_prefix, forcing_wind_unit, forcing_wind_const_wspeed, forcing_input%wind_speed)

  END SUBROUTINE init_forcing

  !>
  !! @brief Provides forcing data for const, timestep or daily forcing (reads data from netcdf if required)
  !!
  !! @par History
  !! First version by Julia Nabel (14.01.16)
  !!
  SUBROUTINE get_data_and_update_field(nblks, ics, ice, current_datetime, model_id, &
    &                                  forcing_input_variable, current_forcing_array)

    INTEGER,                   INTENT(IN)    :: nblks, ics, ice, model_id
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
    CHARACTER(len=NF_MAX_NAME)  :: readUnit

    REAL(wp), POINTER :: return_pointer(:,:,:) !temporary pointer
    REAL(wp) :: hlp_r(1)

    TYPE(t_jsb_model), POINTER :: model

    model => Get_model(model_id)

    !---------------------------
    frequency = forcing_input_variable%frequency

    ! ----- in case of constant forcing, only set correct constants and retun
    IF (frequency == CONST_) THEN
      current_forcing_array(:,:) = forcing_input_variable%constantValue
      CALL set_ocean_cells(nblks, ics, ice, model_id, forcing_input_variable%ocean_value, current_forcing_array)
      CALL set_miss_cells(nblks, forcing_input_variable%missval, forcing_input_variable%ocean_value, current_forcing_array)
      RETURN
    ENDIF

    ! ----- times
    new_year_forcing = is_newyear(current_datetime, delta_time)
    new_day_forcing  = is_newday (current_datetime, delta_time)

    CALL get_date_components(current_datetime, year, month, day, hour, minute, second)

    ! Find year of forcing file and date for new datetime
    do_advance = .TRUE.
    IF (forcing_options%cyclic_nyears > 0) THEN
      cyclic_year = forcing_options%cyclic_start_year &
        &           + MOD(year - get_year_at_experiment_start(), forcing_options%cyclic_nyears)
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
        & int(86400._wp/forcing_options%forcing_steps_per_day) + 1 !+1 because 00:00/X => 0; and 1 based vector
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
!            & TRIM(forcing_input_variable%expectedUnit), TRIM(readUnit))
          CALL forcing_input_variable%input_file%Close()
        END IF

        IF(ASSOCIATED(forcing_input_variable%data_read)) DEALLOCATE(forcing_input_variable%data_read)
        ALLOCATE(forcing_input_variable%data_read(ice, nblks, forcing_options%forcing_steps_per_day))
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
          CALL nf(nf_inq_varid(forcing_input_variable%input_file%file_id, CO2_varname, nc_var_id), routine)
          CALL nf(nf_get_vara_double(forcing_input_variable%input_file%file_id, nc_var_id, &
            & (/1,1,currentTimeStep/), (/1,1,1/), hlp_r), routine)
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
    ENDIF

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
          & TRIM(forcing_input_variable%expectedUnit), TRIM(readUnit))
        ! Read missing value for variable
        IF (forcing_options%forcing_set_miss_to_constants) THEN
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
      ALLOCATE(forcing_input_variable%data_read(ice, nblks, forcing_options%forcing_steps_per_day))

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
          IF(thisDayEnd-thisDayStart+1 /= forcing_options%forcing_steps_per_day) CALL finish(TRIM(routine), &
             & 'Forcing file does not contain the demanded number of time steps per day: '  &
             & //TRIM(int2string(forcing_options%forcing_steps_per_day)) //' -- but: '                      &
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
        DO index_time = 1, forcing_options%forcing_steps_per_day
          CALL set_ocean_cells(nblks, ics, ice, model_id, &
            & forcing_input_variable%ocean_value, forcing_input_variable%data_read(:,:,index_time))
          CALL set_miss_cells(nblks, forcing_input_variable%missval, forcing_input_variable%ocean_value, &
            & forcing_input_variable%data_read(:,:,index_time))
        ENDDO
        IF (do_advance) THEN
          forcing_input_variable%currentDayStart=forcing_input_variable%currentDayEnd+1
          forcing_input_variable%currentDayEnd=forcing_input_variable%currentDayEnd+forcing_options%forcing_steps_per_day
        END IF
      ELSE IF (forcing_input_variable%frequency == DAILY_) THEN
        CALL set_ocean_cells(nblks, ics, ice, model_id, forcing_input_variable%ocean_value, &
          &                  forcing_input_variable%data_read(:,:,1))
        CALL set_miss_cells(nblks, forcing_input_variable%missval, forcing_input_variable%ocean_value, &
          &                 forcing_input_variable%data_read(:,:,1))
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
  !! @brief sets given values for input_variable_type fields in forcing_input
  !!
  !! @par History
  !! Based on jsbach3 mo_jsbalone_forcing init_variable_access    by Julia Nabel
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

  !>
  !! @brief assigns forcing frequency keyword to integer constant
  !!
  !! @par History
  !! Based on jsbach3 mo_jsbalone_forcing init_variable_access by Julia Nabel (2016-01-14)
  !!
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

  !>
  !! @brief gets the length and the content of the time dimension
  !!
  !! @par History
  !! First version by Julia Nabel (2016-02-25)
  !!
  SUBROUTINE get_time_dimension(input_file, forcing_input_variable)

    TYPE(t_input_file),        INTENT(IN)    :: input_file
    TYPE(input_variable_type), INTENT(INOUT) :: forcing_input_variable

    INTEGER                                  :: nc_dim_id, nc_var_id

    CHARACTER(LEN=*), PARAMETER              :: routine = modname//"::get_time_dimension"

    CALL nf(nf_inq_dimid(input_file%file_id, 'time', nc_dim_id), routine)
    CALL nf(nf_inq_dimlen(input_file%file_id, nc_dim_id, forcing_input_variable%length_time_series), routine)
    IF(ASSOCIATED(forcing_input_variable%timevalues)) DEALLOCATE(forcing_input_variable%timevalues)
    ALLOCATE(forcing_input_variable%timevalues(forcing_input_variable%length_time_series))

    ! We need the timevalues for get_data_and_update_field
    CALL nf(nf_inq_varid(input_file%file_id, 'time', nc_var_id), routine)
    CALL nf(nf_get_vara_double(input_file%file_id, nc_var_id, (/ 1 /), &
      & (/ forcing_input_variable%length_time_series /), forcing_input_variable%timevalues), routine)

  END SUBROUTINE get_time_dimension

  !>
  !! @brief tests if the read unit equals the expected unit (with some exception tests)
  !!
  !! @par History
  !! First version by Julia Nabel (2016-02-25)
  !!
  !! TODO: via mo_jsb_io_netcdf_iface / mo_jsb_io_netcdf functions / mo_input?
  !!
  SUBROUTINE test_unit(file_id, variable_name, expectedUnit, readUnit)

    INTEGER, INTENT(IN) :: file_id
    CHARACTER(len=*),INTENT(IN) :: variable_name
    CHARACTER(len=*),INTENT(IN) :: expectedUnit
    CHARACTER(len=*),INTENT(IN) :: readUnit

    CHARACTER(len=*), PARAMETER :: routine = modname//':test_unit'

    INTEGER :: nc_var_id, status

    IF (TRIM(expectedUnit) .NE. TRIM(unspecified_exp_unit)) THEN
      CALL nf(nf_inq_varid(file_id, TRIM(variable_name), nc_var_id), routine)
      status = nf_get_att_text(file_id, nc_var_id, name_of_unit_attribute, readUnit)
      IF (status /= NF_NOERR) THEN
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

  !>
  !! @brief reads and returns missing value for variable fields
  !!
  !! @par History
  !! First version by Tobias Stacke (2022-12-23)
  !!
  SUBROUTINE get_missval(file_id, variable_name, missval)

    INTEGER,         INTENT(IN)  :: file_id
    CHARACTER(len=*),INTENT(IN)  :: variable_name
    REAL(wp),        INTENT(OUT) :: missval

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_missval'

    INTEGER :: nc_var_id, status

    CALL nf(nf_inq_varid(file_id, TRIM(variable_name), nc_var_id), routine)
    status = nf_get_att_double(file_id, nc_var_id, 'missing_value', missval)

    IF (status == NF_NOERR) THEN
      CALL message(TRIM(routine), &
        & 'Found missing value '//TRIM(real2string(missval))//' for variable '//TRIM(variable_name))
    ELSE
      status = nf_get_att_double(file_id, nc_var_id, '_FillValue', missval)
      IF (status == NF_NOERR) THEN
        CALL message(TRIM(routine), &
          & 'Found missing value '//TRIM(real2string(missval))//' for variable '//TRIM(variable_name))
      ELSE
        CALL finish(TRIM(routine), &
          & 'No missing value found for variable '//TRIM(variable_name))
      ENDIF
    ENDIF

  END SUBROUTINE get_missval

  SUBROUTINE finalize_external_forcing

    IF(forcing_options%air_temperature_as_timestep) THEN
      IF (forcing_input%air_temp%input_file%is_open) CALL forcing_input%air_temp%input_file%Close()
    ELSE
      IF (forcing_input%air_temp_min%input_file%is_open) CALL forcing_input%air_temp_min%input_file%Close()
      IF (forcing_input%air_temp_max%input_file%is_open) CALL forcing_input%air_temp_max%input_file%Close()
    END IF

    IF (forcing_input%precipitation%input_file%is_open) CALL forcing_input%precipitation%input_file%Close()

    SELECT CASE(forcing_options%type_of_qair_forcing)
    CASE(RH_)
      IF (forcing_input%rel_humidity%input_file%is_open) CALL forcing_input%rel_humidity%input_file%Close()
    CASE(QAIR_)
      IF (forcing_input%qair%input_file%is_open) CALL forcing_input%qair%input_file%Close()
    END SELECT

    IF (forcing_input%shortwave%input_file%is_open) CALL forcing_input%shortwave%input_file%Close()

    IF (forcing_input%longwave%input_file%is_open) CALL forcing_input%longwave%input_file%Close()

    IF (forcing_input%wind_speed%input_file%is_open) CALL forcing_input%wind_speed%input_file%Close()

    IF (forcing_input%CO2_concentr%input_file%is_open) CALL forcing_input%CO2_concentr%input_file%Close()

  END SUBROUTINE finalize_external_forcing

  !>
  !! @brief calculates the remaining input required for jsbach4
  !!
  !! @par History
  !! Copied required parts from jsbach3: mo_jsbalone: update_driving and adapted for jsb4   by Julia Nabel
  !!
  SUBROUTINE get_remaining_input_for_jsb4_interface(nblks, ice, icz, temp_air, pressure, qair, wind, &
      & t_srf_proc, fact_q_air_proc, fact_qsat_srf_proc, rough_h_srf_proc, rough_m_srf_proc, &
      & drag_srf, t_acoef, t_bcoef, q_acoef, q_bcoef, pch)

    INTEGER,  INTENT(IN)  :: nblks, ice, icz
    REAL(wp), INTENT(IN)  :: temp_air(:,:)
    REAL(wp), INTENT(IN)  :: pressure(:,:)
    REAL(wp), INTENT(IN)  :: qair(:,:)
    REAL(wp), INTENT(IN)  :: wind(:,:)
    REAL(wp), INTENT(IN)  :: t_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: fact_q_air_proc(:,:)
    REAL(wp), INTENT(IN)  :: fact_qsat_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: rough_h_srf_proc(:,:)
    REAL(wp), INTENT(IN)  :: rough_m_srf_proc(:,:)

    REAL(wp), INTENT(OUT) :: drag_srf(ice,nblks)
    REAL(wp), INTENT(OUT) :: t_acoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: t_bcoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: q_acoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: q_bcoef(ice,nblks)
    REAL(wp), INTENT(OUT) :: pch(ice,nblks)

    INTEGER  :: ic, iblk
    REAL(wp) :: zcons9, zcons11, zcons12, zsigma
    REAL(wp), DIMENSION(ice):: air_pressure, zdu2, ztvd, ztvir, qsat_surf
    REAL(wp), DIMENSION(ice):: ztvs, zg, zgh, zril, zcons, zchnl, zcfnchl, zcfhl, zscfl, zucfhl

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_remaining_input_for_jsb4_interface'

    ! in jsb3: rd = GasConstantDryAir; cpd = SpecificHeatDryAirConstPressure; grav = Gravity
    zcons9 = 3._wp * cb
    zcons11 = 3._wp * cb * cc
    zcons12 = time_step_len * grav / rd

    ! Initializations needed for grid cells beyond icz
    drag_srf(:,:) = 0._wp
    t_acoef(:,:)  = 0._wp
    t_bcoef(:,:)  = 0._wp
    q_acoef(:,:)  = 0._wp
    q_bcoef(:,:)  = 0._wp
    pch(:,:)      = 0._wp

    DO iblk = 1,nblks

      IF (iblk == nblks) THEN
        ic = icz
      ELSE
        ic = ice
      END IF

      !------------------------------------------------------------------------------------
      ! Approximation of cdrag
      !------------------------------------------------------------------------------------
      ! squared wind shear ! minimum wind speed square from echam/mo_surface_land.f90:precalc_land
      zdu2(1:ic) = MAX(wind(1:ic,iblk)**2, 1._wp)

      ! virtual potential air temperature (see mo_surface_boundary.f90)
      ! according to Saucier, WJ Principles of Meteoroligical Analyses
      ! tv = t * (1 + 0.61 * q )    ! virtual temperature
      ! td = t * ( 1000 / p_mb ) ^ R/cdp  ! potential temperature
      ! tvd = tair * (100000/p_pa)^rd_o_cpd * (1 + rvd1 * qair) ! virtual potential temperature
      zsigma    = 0.99615_wp        ! corresponds to lowest level of 47 layer echam
      air_pressure(1:ic) = zsigma * pressure(1:ic, iblk)
      ztvd(1:ic) = ( temp_air(1:ic, iblk) * ( 100000._wp/air_pressure(1:ic))**rd_o_cpd ) * ( 1._wp + rvd1 * qair(1:ic, iblk))
      ztvir(1:ic) = temp_air(1:ic, iblk) * ( 1._wp + rvd1 * qair(1:ic, iblk) )

      ! virtual potential surface temperature
      CALL sat_specific_humidity(ic, t_srf_proc(1:ic,iblk), pressure(1:ic,iblk), qsat_surf(1:ic))
      ztvs(1:ic) = t_srf_proc(1:ic,iblk) * ( 100000._wp/pressure(1:ic,iblk))**rd_o_cpd * &
        & ( 1._wp + rvd1 * ( fact_qsat_srf_proc(1:ic,iblk) * qsat_surf(1:ic) + ( 1._wp - fact_q_air_proc(1:ic,iblk) ) &
        & * qair(1:ic,iblk)))

      ! geopotential of the surface layer (see echam's auxhybc.f90 & geopot.f90)
      ! adapted according to jsb3 offline modifications by Marvin, Philipp et al. Jan. 24, 2017 (r8893)
      IF(forcing_options%heightWind > 0._wp .AND. forcing_options%heightHumidity > 0._wp) THEN
        zg(1:ic) = forcing_options%heightWind * grav
        zgh(1:ic) = forcing_options%heightHumidity * grav ! jsb3 "HeightHumidity and HeightTemperature have to be equal"
      ELSE
        zg(1:ic) = ztvir(1:ic) * rd * LOG(1._wp / zsigma)
        zgh(1:ic) = zg(1:ic)
      ENDIF

      ! Richardson number (dry, Brinkop & Roeckner 1995, Tellus)
      ! ztvd, ztvs are now virtual potential temperatures, changed by Thomas Raddatz 07.2014
      zril(1:ic) = zg(1:ic) * ( ztvd(1:ic) - ztvs(1:ic) ) / ( zdu2(1:ic) * (ztvd(1:ic) + ztvs(1:ic))/2._wp )

      ! Neutral drag coefficient for momentum and heat
      zchnl(1:ic) = von_karman**2 / (LOG(1._wp + zg(1:ic) / (grav * rough_m_srf_proc(1:ic,iblk) )) &
                * LOG( ( grav * rough_m_srf_proc(1:ic,iblk) + zgh(1:ic) ) / (grav * rough_h_srf_proc(1:ic,iblk) )))

      ! account for stable/unstable case: helper variables
      zscfl(1:ic) = SQRT ( 1._wp + 5._wp * ABS(zril(1:ic)))
      zucfhl(1:ic) = 1._wp / (1._wp + zcons11 * zchnl(1:ic) * SQRT(ABS(zril(1:ic)) * (1._wp &
        & + zgh(1:ic) / (grav * rough_h_srf_proc(1:ic,iblk)))))

      ! ignoring cloud water correction (see mo_surface_land.f90)
      zcons(1:ic) = zcons12 * pressure(1:ic,iblk) / ( temp_air(1:ic,iblk) * (1._wp + rvd1 * qair(1:ic,iblk)))
      zcfnchl(1:ic)  = zcons(1:ic) * SQRT(zdu2(1:ic)) * zchnl(1:ic)

      ! Stable / Unstable case
      WHERE ( zril(1:ic) .GT. 0._wp )
        zcfhl(1:ic) = zcfnchl(1:ic) / (1._wp + zcons9 * zril(1:ic) * zscfl(1:ic))
      ELSEWHERE
        zcfhl(1:ic) = zcfnchl(1:ic) * (1._wp - zcons9 * zril(1:ic) * zucfhl(1:ic))
      END WHERE

      drag_srf(1:ic,iblk) = zcfhl(1:ic)
      pch(1:ic,iblk) = zcfhl(1:ic) / zcfnchl(1:ic) * zchnl(1:ic)

      !---------------------------------------------------------------------------------------------------------------
      ! Computation of Richtmeyr-Morton Coefficients
      ! This follows now Jan Polcher's explicit solution, i.e. atmospheric conditions at t+1 are assumed to be valid
      !---------------------------------------------------------------------------------------------------------------
      t_acoef(1:ic,iblk) = 0.0_wp
      t_bcoef(1:ic,iblk) = cpd * (1 + cpvd1 * qair(1:ic,iblk)) * temp_air(1:ic,iblk) + zgh(1:ic)
      q_acoef(1:ic,iblk) = 0.0_wp
      q_bcoef(1:ic,iblk) = qair(1:ic,iblk)

    ENDDO

  END SUBROUTINE get_remaining_input_for_jsb4_interface

  !>
  !! @brief converts daily relative humidity to specific humidity considering the daily course of temperature
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: spec_humidity_from_rel_humidity and adapted for jsb4 input   by Julia Nabel
  !!
  SUBROUTINE spec_humidity_from_rel_humidity(ice, sinlat, coslat, rel_humidity, Tmin, Tmax, elevation, qair)
    ! This subroutine converts daily relative humidity to specific humidity considering the daily course of temperature
    ! ops! Using this routine only makes sense together with the generated diurnal cycle of temperature as in
    ! SUBROUTINE instantly_from_daily_temp_2
    ! Otherwise the results will be faulty!

    IMPLICIT NONE

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

    CALL inquire_declination(declination)

    nstep = (86400.0_wp/delta_time)

    IF(.NOT.(ALLOCATED(tmp))) ALLOCATE(tmp(time_steps_per_day))
    IF(.NOT.(ALLOCATED(pressure))) ALLOCATE(pressure(time_steps_per_day))
    IF(.NOT.(ALLOCATED(qsat))) ALLOCATE(qsat(time_steps_per_day))
    IF(.NOT.(ALLOCATED(inv_qsat))) ALLOCATE(inv_qsat(time_steps_per_day))
    IF(.NOT.(ALLOCATED(qair_test))) ALLOCATE(qair_test(time_steps_per_day))

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
      ENDDO

      ! from function bottom pressure; in jsb3: p0sl_bg = p_sealevel; grav = Gravity; rd = GasConstantDryAir
      pressure(:) = &
        & p0sl_bg * ( 1._wp - elevation(i)*gamma /(tmelt+tmp(:) + elevation(i)*gamma))**(grav/rd/gamma)

      CALL sat_specific_humidity(time_steps_per_day, tmp(:) + tmelt, pressure(:), qsat(:)) ! present saturated qair

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

    IF((ALLOCATED(tmp))) DEALLOCATE(tmp)
    IF((ALLOCATED(pressure))) DEALLOCATE(pressure)
    IF((ALLOCATED(qsat))) DEALLOCATE(qsat)
    IF((ALLOCATED(inv_qsat))) DEALLOCATE(inv_qsat)
    IF((ALLOCATED(qair_test))) DEALLOCATE(qair_test)
  END SUBROUTINE spec_humidity_from_rel_humidity

  !>
  !! @brief Returns saturation specific humidity for given temperature and pressure
  !!
  !! @par History
  !! Copied from jsbach3: mo_atmosphere: sat_specific_humidity    by Julia Nabel (2016)
  !!
  SUBROUTINE sat_specific_humidity(temp_dim, temp, pressure, qsat)
    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)

    INTEGER,  INTENT(IN) :: temp_dim
    REAL(wp), INTENT(IN) :: temp(:)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(:)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sat_specific_humidity'

    INTEGER :: i
    INTEGER  :: it(1:temp_dim)
    REAL(wp) :: tluc(1:temp_dim)

    !Assertion: the expected temperature range is indirectly defined by the tlucua bounds -> consistency check
    DO i = 1,temp_dim
      IF ((temp(i) .LT. REAL(LBOUND(tlucua,dim=1)/1000._wp)) .OR. (temp(i) .GT. REAL(UBOUND(tlucua,dim=1)/1000._wp))) THEN
        WRITE (message_text,*) &
          & 'Inconsistency with tlucua bounds (would give a seg. fault or a floating invalid operation exception)', &
          & '- temperature was: ', temp(i), 'K. One thing to check: consistency of land-sea mask and forcing.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ENDIF
    ENDDO

    it = NINT(temp*1000._wp)
    tluc = tlucua(it)

    WHERE (it >= jptlucu1 .AND. it <= jptlucu2)
      qsat = tluc / (pressure - rvd1*tluc)
    ELSEWHERE
      qsat = HUGE(1._wp)
    END WHERE

  END SUBROUTINE sat_specific_humidity

  !>
  !! @brief Returns saturation specific humidity for given air and vapour pressure
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: specific_humidity    by Julia Nabel
  !
  ! === specific_humidity() =======================================================================================================
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
  PURE elemental FUNCTION fun_specific_humidity(vapour_pressure, air_pressure)
    REAL(wp)            :: fun_specific_humidity !! (from [0,1])
    REAL(wp),INTENT(in) :: vapour_pressure   !! [N/m^2]
    REAL(wp),INTENT(in) :: air_pressure      !! air pressure at bottom [N/m^2]

    !in jsbach3 rdv=eps
    fun_specific_humidity = rdv*vapour_pressure/(air_pressure - (1._wp-rdv)*vapour_pressure)
  END FUNCTION fun_specific_humidity

  !>
  !! @brief computes vapor pressure from actual and potential evapotranspiration
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: vapour_pressure_from_evapor    by Julia Nabel
  !
  ! === vapour_pressure_from_evapor() =============================================================================================
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
  !vg>>
  ! Above values for h0 and h1 lead to extremely moist conditions, especially in desert areas.
  ! A comparison of qair from echam/jsbach output with qair calculated in this routine using different h0 values shows best
  ! agreement for h0=0.6. The same h0 is found when comparing CRUNCEP qair with qair calculated from CRUNCEP temperatures in
  ! this module. For desert studies it might be good to reduce h0 even further. As the specific humidity does not have a
  ! significant diurnal cycle, h1 is set to 0.
  !vg<<
  FUNCTION vapour_pressure_from_evapor(nblks, ice, icz, evapo_ratio_mean, Tmin, Tair)
    INTEGER,  INTENT(IN) :: nblks, ice, icz
    REAL(wp), INTENT(IN) :: evapo_ratio_mean(:,:)  !! fe(d-1) of Eq. (1)
    REAL(wp), INTENT(IN) :: Tmin(:,:)              !! minimum temperature (= temperature at sunrise)
    REAL(wp), INTENT(IN) :: Tair(:,:)              !! temperature at current time step

    CHARACTER(len=*), PARAMETER :: routine = modname//':vapour_pressure_from_evapor'

    REAL(wp), PARAMETER :: h0 = 0.6_wp  !! (0.96_wp) tuning constant for computing vapor pressure
    REAL(wp), PARAMETER :: h1 = 0.0_wp  !! (0.49_wp) tuning constant for computing vapor pressure

    REAL(wp) :: vapour_pressure_from_evapor(ice, nblks) !! the vapour pressure
    REAL(wp) :: evapo_sat        ! saturation pressure at actual temperature
    REAL(wp) :: evapo_sat_min    ! saturation pressure at minimum temperature of day (= at sunrise)
    REAL(wp) :: e0_vap           ! vapor pressure at sunrise

    INTEGER :: i,j,ic

    DO j = 1,nblks
      IF (j == nblks) THEN
        ic = icz
      ELSE
        ic = ice
      END IF
      DO i = 1,ic
        !Assertion: the expected temperature range is indirectly defined by the tlucua bounds -> consistency check
        IF ((Tmin(i,j) .LT. REAL(LBOUND(tlucua,dim=1)/1000._wp)-tmelt) &
            & .OR. (Tmin(i,j) .GT. REAL(UBOUND(tlucua,dim=1)/1000._wp)-tmelt)) THEN
          WRITE (message_text,*) &
            & 'Inconsistency with tlucua bounds (would give a seg. fault or a floating invalid operation exception)', &
            & '- tmin temperature was: ', Tmin(i,j), ' degC. One thing to check: consistency of land-sea mask and forcing.'
          CALL finish(TRIM(routine), TRIM(message_text))
        ENDIF
        IF ((Tair(i,j) .LT. REAL(LBOUND(tlucua,dim=1)/1000._wp)-tmelt) &
            & .OR. (Tair(i,j) .GT. REAL(UBOUND(tlucua,dim=1)/1000._wp)-tmelt)) THEN
          WRITE (message_text,*) &
            & 'Inconsistency with tlucua bounds (would give a seg. fault or a floating invalid operation exception)', &
            & '- air temperature was: ', Tmin(i,j), ' degC. One thing to check: consistency of land-sea mask and forcing.'
          CALL finish(TRIM(routine), TRIM(message_text))
        ENDIF

        !jsb3: rdv = eps
        evapo_sat_min = tlucua(INT(1000._wp * (Tmin(i,j) + tmelt)))/rdv ! saturation pressure at minimum temperature
        evapo_sat     = tlucua(INT(1000._wp * (Tair(i,j) + tmelt)))/rdv ! saturation pressure at actual temperature
        e0_vap = evapo_sat_min * (h0 + (1._wp - h0) * evapo_ratio_mean(i,j))                               ! Eq. (2) from above
        vapour_pressure_from_evapor(i,j) = e0_vap + h1 * evapo_ratio_mean(i,j) * (evapo_sat - evapo_sat_min) ! Eq. (3) from above

      ENDDO
    ENDDO
  END FUNCTION vapour_pressure_from_evapor


  !>
  !! @brief sets all not land (i.e. not notsea ;)) cells to the given constant
  !!
  !! @par History
  !! First version    by Julia Nabel (20.01.16)
  !
  !   Some external forcing is only defined for land cells.
  !   If the external forcing is only defined for land cells is specified by the namelist option
  !   "forcing_options%forcing_set_ocean_to_constants". If true: ocean cells need to be replaced
  !   by constants (reasonable values for the different variables) to ensure numeric stability.
  !   If false: forcing is expected for all variables for all cells (e.g. forcing from echam)
  SUBROUTINE set_ocean_cells(nblks, ics, ice, model_id, constant, array)

    USE mo_jsb_model_class, ONLY: t_jsb_model
    USE mo_jsb_class,       ONLY: Get_model

    !dsl4jsb_Def_memory(SRF)

    INTEGER, INTENT(IN) :: nblks, ics, ice, model_id
    REAL(wp), INTENT(IN) :: constant
    REAL(wp), INTENT(INOUT) :: array(:,:)
    INTEGER :: iblk

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    REAL(wp) :: fract(ice-ics+1)

    IF (.NOT. forcing_options%forcing_set_ocean_to_constants) RETURN

    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)

    DO iblk = 1, nblks
      CALL tile%Get_fraction(ics, ice, iblk, fract=fract(:))
      WHERE (fract(:) == 0._wp)
        array(:,iblk) = constant
      END WHERE
    ENDDO
  END SUBROUTINE set_ocean_cells

  !>
  !! @brief replace missing values with the given constant
  !!
  !! @par History
  !! First version by Tobias Stacke (2022-12-13)
  !
  !  Land sea masks of model grids and offline forcing data don't alway match exactly. If
  !    the forcing dataset has missing values on land cell, the model will most likely
  !    crash with temperature bound violations due to the missing values. This routine
  !    replaces the missing value with a constant. Currently, the "ocean" value is
  !    used for this.
  SUBROUTINE set_miss_cells(nblks, missval, constant, array)

    INTEGER, INTENT(IN)     :: nblks
    REAL(wp), INTENT(IN)    :: missval, constant
    REAL(wp), INTENT(INOUT) :: array(:,:)
    INTEGER :: iblk

    IF (.NOT. forcing_options%forcing_set_miss_to_constants) RETURN

    ! Use REAL values for comparison in case the netCDF missing_value
    ! attribute lacks double precision.

    DO iblk = 1, nblks
      WHERE (ABS(REAL(array(:,iblk)) - REAL(missval)) <= EPSILON(1.0))
        array(:,iblk) = constant
      END WHERE
    ENDDO

  END SUBROUTINE set_miss_cells

  !>
  !! @brief estimate timestep temperature from daily min- and max temperature
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: instantly_from_daily_temp_2 and adapted for 2D arrays    by Julia Nabel
  !
  ! --- instantly_from_daily_temp_2() --------------------------------------------------------------------------------------
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
  SUBROUTINE instantly_from_daily_temp_2(nblks, ice, icz, lon, sinlat, coslat, dtmp_min, dtmp_max, date, tmp)

    IMPLICIT NONE

    REAL(wp), PARAMETER :: pio180 = pi/180._wp ! conversion factor from degrees to radians ("pi over 180")
    REAL(wp), PARAMETER :: rat1=14._wp/24._wp  ! day fraction at which maximum temperature shall occur
    REAL(wp), PARAMETER :: rat2=4._wp/24._wp   ! minimum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat3=20._wp/24._wp  ! maximum length of normal day (as fraction of day)
    REAL(wp), PARAMETER :: rat4=2._wp/24._wp   ! day fraction of dawn times (before sunrise + after sunset), used to set ...
                                               ! ... temperature of a normal day to T_mean at sunrise-rat4/2 and sunset+rat4/2

    INTEGER, INTENT(IN) :: nblks, ice, icz
    REAL(wp), INTENT (IN) :: lon(:,:)       ! minimal temperature at the considered day
    REAL(wp), INTENT (IN) :: sinlat(:,:)    ! maximal range at that day
    REAL(wp), INTENT (IN) :: coslat(:,:)    ! maximal range at that day
    REAL(wp), INTENT (IN) :: dtmp_min(:,:)  ! minimal temperature at the considered day
    REAL(wp), INTENT (IN) :: dtmp_max(:,:)  ! maximal range at that day
    TYPE (time_days), INTENT (IN)  :: date  ! the considered simulation instant at that day in time-seconds format
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
    INTEGER    :: i, j, ic

    tmp(:,:) = 0._wp

    CALL inquire_declination(declination)

    year_instant=get_year_day(date)            ! conversion of current time step from day-seconds to yearday-dayfraction format
    yearday = INT(year_instant)                ! number of day of the current time step since the beginning of the current year
    fract_timestep=year_instant-REAL(yearday,wp)   ! fractional part of the day of the current time step since midnight in UTC

    DO j = 1,nblks
      IF (j == nblks) THEN
        ic = icz
      ELSE
        ic = ice
      END IF

      DO i = 1,ic
        fract_UTC_offset = lon(i,j)/360._wp   ! difference between local time and UTC as fraction of a day
        fract_localtime = MOD(fract_timestep+fract_UTC_offset+2._wp,1._wp) ! local time as fraction of the local day (modulo 1)
        spds = sinlat(i,j)*SIN(declination)  ! sine contribution to solar zenith angle
        cpds = coslat(i,j)*COS(declination)  ! cosine contribution to solar zenith angle
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

        dtmp = arg * dtmp_max(i,j) + (1._wp-arg) * dtmp_min(i,j) ! The dusk day temperature
        !dtmp = 0.5_wp*(dtmp_max(i,j)+dtmp_min(i,j)) ! The mean day temperature
        dtran = dtmp_max(i,j)- dtmp_min(i,j)        ! temperature range
        IF (fract_daylight>=rat2 .AND. fract_daylight<=rat3) THEN ! -------- normal day ---------------------------------------------
          fract_sunrise = 0.5_wp - fract_daylight/2._wp             ! time of sunrise expressed as fraction of day since midnight
          fract_sunset = fract_sunrise+fract_daylight                ! time of sunset expressed as fraction of day since midnight
          IF (fract_localtime > fract_sunrise .AND. fract_localtime <= fract_sunset) THEN  ! For time steps at daylight ...
            sd = COS(pi*(fract_localtime-rat1)/(fract_daylight/2._wp+rat4)) ! ... modulation factor for day-temperature such that
                                                                    ! .. temperature is at maximum at rat1 (2 pm) and at average at
                                                                    ! .. sunrise-rat4/2 and sunset+rat4/2
            IF (sd >= 0) THEN
              tmp(i,j) = dtmp + sd*(dtmp_max(i,j)-dtmp)                       ! temperature modulated by temperature range
            ELSE
              tmp(i,j) = dtmp + sd*(dtmp-dtmp_min(i,j))                       ! temperature modulated by temperature range
            ENDIF
            !tmp(i,j) = dtmp + sd*dtran/2._wp                       ! temperature modulated by temperature range
          ELSE                                                 ! For time steps at night ...
            sd = COS(pi*(fract_sunset-rat1)/(fract_daylight/2._wp+rat4)) ! modulation factor for temperature at sunset
            IF (sd >= 0) THEN
              temp_sunset = dtmp + sd*(dtmp_max(i,j)-dtmp)                       ! temperature modulated by temperature range
            ELSE
              temp_sunset = dtmp + sd*(dtmp-dtmp_min(i,j))                       ! temperature modulated by temperature range
            ENDIF
            !temp_sunset = dtmp + dtran/2._wp*sd                     ! temperature at sunset
            fract_since_sunset = MOD(fract_localtime-fract_sunset+1._wp,1._wp) ! day fraction since sunset
            tmp(i,j) = dtmp_min(i,j)+(temp_sunset-dtmp_min(i,j)) &
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
          tmp(i,j) = dtmp + dtran/2._wp*sd                             ! temperature modulated by temperature range
        ELSE ! ------------------------------- short day (close to polar night) ---------------------------------------------------
          tmp(i,j) = dtmp                                        ! temperature is set to mean temperture at all time steps
        ENDIF
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE INSTANTLY_FROM_DAILY_TEMP_2

  !>
  !! @brief calculate bottom pressure from air temperature and elevation
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: bottom pressure and only slightly adapted for jsb4 context    by Julia Nabel
  !
  ! === bottom pressure() =======================================================================================================
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
  ELEMENTAL FUNCTION bottom_pressure(elevation,air_temperature)
    REAL(wp)              :: bottom_pressure
    REAL(wp), INTENT(IN)  :: elevation
    REAL(wp), INTENT(IN)  :: air_temperature

    ! in jsb3: p0sl_bg = p_sealevel; grav = Gravity; rd = GasConstantDryAir
    bottom_pressure = &
      & p0sl_bg * ( 1._wp - elevation*gamma /(tmelt + air_temperature + elevation*gamma))**(grav/rd/gamma)
  END FUNCTION bottom_pressure

  !>
  !! @brief split precipitation in snow and rain depending on mean air temperature
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: snow_and_rain_from_precip    by Julia Nabel
  !
  ! === snow_and_rain_from_precip() ==============================================================================================
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
    REAL(wp),INTENT(OUT) :: rain          !! rain [kg/m^2/s]
    REAL(wp),INTENT(OUT) :: snow          !! snow [kg/m^2/s]
    REAL(wp),INTENT(IN)  :: precipitation !! total precipitation (rain+snow) [kg/m^2/s]
    REAL(wp),INTENT(IN)  :: air_temp_daily!! daily mean air_temperature at surface [Celsius]

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

  !>
  !! @brief computes shortwave parts (PAR,NIR,fraction of diffuse PAR)
  !! ATTENTION: TR: there maybe considerable deviations in the daily sum of insolation if daily data is used and
  !! if the orbit is different then todays orbit (i.e. ATTENTION: this may be incorrect for paleo situations)
  !!
  !! @par History
  !! JN: copied from jsbach3: mo_jsbalone_forcing: shortwave_from_direct_shortwave (+ comment parts from shortwave_from_fpar_orig)
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
  ELEMENTAL SUBROUTINE shortwave_from_direct_shortwave(day_of_year,cos_zenith,cos_lat,sin_lat,air_pressure, &
    & rad_sw_in, rad_UV,rad_PAR,rad_NIR,fract_PAR_diffuse,rad_sw,rad_sw_pot )

    INTEGER,INTENT(in)   :: day_of_year      !! day in year (from [1,365])
    REAL(wp),INTENT(in)  :: cos_zenith       !! cosine of zenith angle
    REAL(wp),INTENT(in)  :: cos_lat          !! cosine of latitude
    REAL(wp),INTENT(in)  :: sin_lat          !! sine of latitude
    REAL(wp),INTENT(in)  :: air_pressure     !! air pressure at bottom (depending on local elevation) [N/m^2]
    REAL(wp),INTENT(in)  :: rad_sw_in        !! "Rtot": either timestep or daily average total solar shortwave radiation :
                                             !! ... direct+diffuse, 280-3000 nm, (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(wp),INTENT(out) :: rad_UV           !! "RUV": solar radiation flux (direct+diffuse) from the UV band 280-400 nm
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: rad_PAR          !! "Rpar_act": solar radiation flux (direct+diffuse) from the visible band 400-700 nm .
                                             !! .. arriving actually at the bottom (includes cloud shading and scattering) [W/m^2]

    REAL(wp),INTENT(out) :: rad_NIR          !! "Rnir": flux (direct+diffuse) from the near infrared band 700-3000 nm ..
                                             !! .. actually arriving at the bottom (includes cloud shading and scattering) [W/m^2]
    REAL(wp),INTENT(out) :: fract_PAR_diffuse !! "fdiffPar": Fraction from  rad_PAR() that comes down as diffuse radiation ..
                                             !! .. (values in [0,1])
    REAL(wp),INTENT(out) :: rad_sw,rad_sw_pot

    REAL(wp) :: day_angle     !! day of year expressed as angle in radians
    REAL(wp) :: r_sun_earth_2 !! inverse squared earth-sun-distance normalized to mean distance
    REAL(wp) :: Rpar_top      !! PAR flux incident to the outer atmosphere at particular zenith angle [W/m^2]
    REAL(wp) :: fsw           !! "fsw": fraction of potential radiation solar shortwave radiation : direct+diffuse, 280-3000 nm, ,..
                              !! ... (i.e. PAR+NIR) that actually arrives at the bottom [W/m^2]
    REAL(wp) :: F_invers      !! factor 1/F to obtain total shortwave from visible radiation (see Eq. (9))
    REAL(wp) :: hlp_r,hlp2_r,hlp3_r,hlp4_r
    REAL(wp) :: rad_sw_in_corr
    REAL(wp) :: delta,v,u,hh,sinehh,K

    REAL(wp), PARAMETER  ::  fract_UV = 0.1_wp ! UV fraction of (UV + PAR)

    day_angle = 2._wp*pi*(day_of_year+10._wp)/365._wp !! ATTENTION: this may be incorrect for paleo situations

    hlp_r=2._wp*day_angle
    r_sun_earth_2 = &         ! Eq. (3) from above
      & 1.000110_wp+3.4221E-2_wp*COS(day_angle)+1.280E-3_wp*SIN(day_angle)+7.19E-4_wp*COS(hlp_r)+7.7E-5_wp*SIN(hlp_r)
    K=13750.98708_wp

    IF (forcing_input%shortwave%frequency /= TIMESTEP_ .AND. forcing_input%shortwave%frequency /= SUBDAILY_) THEN
      ! today's sum of top of the atmosphere radiation
      delta=-23.4_wp*pi/180._wp*COS(day_angle)
      u=sin_lat*SIN(delta)
      v=cos_lat*COS(delta)
      IF (u.GE.v) THEN
         hh=pi
      ELSEIF (u.LE.(0._wp-v)) THEN
         hh=0.0_wp
      ELSE
         hh=ACOS(-u/v)
      ENDIF
      sinehh=SIN(hh)

      !estimate of noon radiation giving rise to the recorded radiation sum
      rad_sw_in_corr = 0._wp
      IF(u*hh+v*sinehh.GT.1e-5_wp)THEN
         rad_sw_in_corr = rad_sw_in / (2._wp*(u*hh+v*sinehh) * K) * 86400._wp
      ENDIF
      rad_sw = rad_sw_in_corr * cos_zenith

    ELSE
      ! subdaily or time step forcing, use radiation directly
      rad_sw = rad_sw_in
    ENDIF

    IF(cos_zenith > 0.017452406_wp) THEN! zenith angle smaller than 89 degrees  (cos(89) = 0.017452406)
      Rpar_top   = solar_const*cos_zenith*r_sun_earth_2       ! Eq. (2) from above
      hlp3_r = EXP(-0.185_wp*(air_pressure/p0sl_bg)/cos_zenith) ! in jsb3: p0sl_bg = p_sealevel
      hlp4_r = 0.4_wp+0.6_wp*hlp3_r
      rad_sw_pot = Rpar_top * hlp4_r
      fsw = MAX(MIN(rad_sw/rad_sw_pot,1._wp),0._wp)
      hlp_r=1._wp-fsw
      hlp2_r=hlp_r*hlp_r
      F_invers=1._wp + (1.185_wp-0.437_wp*hlp_r-0.494_wp*hlp2_r)* &
        & EXP(MIN(2.534_wp,(0.0305_wp-0.208_wp*hlp_r+0.495_wp*hlp2_r)/cos_zenith)) ! Eq. (9) from above
      rad_UV = rad_sw * fract_UV / F_invers
      rad_PAR = rad_sw * (1._wp - fract_UV) / F_invers   ! This is "Rtot"; Eq. (8) from above
      rad_NIR = rad_sw - rad_PAR - rad_UV ! This is "Rnir"; Eq. (10) from above
      IF(fsw <= 0.2_wp) THEN
         fract_PAR_diffuse = 1._wp                                                 ! Eq. (12)-(15) from above
      ELSEIF(fsw >= 0.9_wp) THEN
         fract_PAR_diffuse = 0.08_wp                                               ! Eq. (12)-(15) from above
      ELSE
         fract_PAR_diffuse = 1._wp - &
           & (1._wp- ((0.9_wp -fsw)/0.7_wp)**(2._wp/3._wp))*hlp3_r/hlp4_r ! Eq. (12)-(15) from above
      ENDIF
    ELSE ! zenith angle larger than 89 degrees - changed Aug. 2017 following T. Raddatz, Jan. 2015
      rad_UV = 0._wp
      rad_PAR = rad_sw * 0.5_wp
      rad_NIR = rad_sw * 0.5_wp
      fract_PAR_diffuse = 1._wp ! (could also be another value in this case)
      rad_sw_pot  = rad_sw
    ENDIF

  END SUBROUTINE shortwave_from_direct_shortwave

  !>
  !! @brief computes longwave from observed longwave and diural course associated with air temperature
  !!
  !! @par History
  !! Copied from jsbach3: mo_jsbalone_forcing: longwave_from_daily_longwave    by Julia Nabel
  !
  ! === longwave_from_daily_longwave() ====================================================================================
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

#endif

END MODULE mo_jsb4_forcing_echam
