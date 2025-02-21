!
! This module is the interface between nwp_nh_interface to the radiation schemes
! (ecRad and RRTM).
!
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

MODULE mo_nwp_rad_interface

  USE mo_exception,            ONLY: finish, message_text
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_parallel_config,      ONLY: nproma
  USE mo_impl_constants,       ONLY: MODIS
  USE mo_kind,                 ONLY: wp
  USE mo_nwp_lnd_types,        ONLY: t_lnd_prog, t_wtr_prog, t_lnd_diag
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_radiation_config,     ONLY: albedo_type, albedo_fixed,                            &
    &                                irad_co2, irad_n2o, irad_ch4, irad_cfc11, irad_cfc12, &
    &                                tsi_radt, ssi_radt, isolrad, cosmu0_dark
  USE mo_radiation,            ONLY: pre_radiation_nwp_steps
  USE mo_nwp_rrtm_interface,   ONLY: nwp_rrtm_radiation,             &
    &                                nwp_rrtm_radiation_reduced
  USE mo_albedo,               ONLY: sfc_albedo, sfc_albedo_modis, sfc_albedo_scm
  USE mtime,                   ONLY: datetime, timedelta, max_timedelta_str_len,                  &
    &                                operator(+), operator(-), newTimedelta, deallocateTimedelta, &
    &                                getPTStringFromSeconds, newDatetime, deallocateDatetime
  USE mo_timer,                ONLY: timer_start, timer_stop, timers_level, timer_preradiaton
  USE mo_nwp_gpu_util,         ONLY: gpu_d2h_nh_nwp
  USE mo_bc_greenhouse_gases,  ONLY: bc_greenhouse_gases_time_interpolation


  USE mo_bc_solar_irradiance,  ONLY: read_bc_solar_irradiance, ssi_time_interpolation
  USE mo_bcs_time_interpolation,ONLY: t_time_interpolation_weights,   &
    &                                 calculate_time_interpolation_weights
  USE mo_o3_util,              ONLY: o3_interface
  USE mo_nwp_aerosol,          ONLY: nwp_aerosol_interface, nwp_aerosol_cleanup
  USE mo_fortran_tools,        ONLY: set_acc_host_or_device
  USE mo_run_config,           ONLY: msg_level

  IMPLICIT NONE

  PRIVATE



  PUBLIC :: nwp_radiation
  

 CONTAINS
  
  !---------------------------------------------------------------------------------------
  !>
  !! This subroutine is the interface between nwp_nh_interface to the radiation schemes.
  !! Depending on inwp_radiation, it can call RRTM (1) or ecRad(4).
  !!
  SUBROUTINE nwp_radiation ( lredgrid, p_sim_time, mtime_datetime, pt_patch,pt_par_patch, &
    & ext_data, lnd_diag, pt_prog, pt_diag, prm_diag, lnd_prog, wtr_prog, zf, zh, dz, lacc)

    CHARACTER(len=*), PARAMETER :: &
      &  routine = 'mo_nwp_rad_interface:nwp_radiation'

    LOGICAL,                 INTENT(in)    :: lredgrid        !< use reduced grid for radiation
    LOGICAL, OPTIONAL,       INTENT(in)    :: lacc ! If true, use openacc

    REAL(wp),                INTENT(in)    :: p_sim_time   !< simulation time
    REAL(wp),                INTENT(in)    :: zf(:,:,:)    !< model full layer height
    REAL(wp),                INTENT(in)    :: zh(:,:,:)    !< model half layer height
    REAL(wp),                INTENT(in)    :: dz(:,:,:)    !< Layer thickness

    TYPE(datetime), POINTER, INTENT(in)    :: mtime_datetime
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_patch     !<grid/patch info.
    TYPE(t_patch), TARGET,   INTENT(in)    :: pt_par_patch !<grid/patch info (parent grid)
    TYPE(t_external_data),   INTENT(inout) :: ext_data
    TYPE(t_lnd_diag),        INTENT(in)    :: lnd_diag   !<diag vars for sfc
    TYPE(t_nh_prog), TARGET, INTENT(inout) :: pt_prog    !<the prognostic variables
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: pt_diag    !<the diagnostic variables
    TYPE(t_nwp_phy_diag),    INTENT(inout) :: prm_diag
    TYPE(t_lnd_prog),        INTENT(inout) :: lnd_prog   ! time level new
    TYPE(t_wtr_prog),        INTENT(in)    :: wtr_prog   ! time level new

    TYPE(datetime) , POINTER :: radiation_time => NULL() !< date and time for radiative transfer and time for radiative transfer
    TYPE(timedelta), POINTER :: td_radiation_offset      !< Offset to center of radiation time step
    TYPE(t_time_interpolation_weights) :: radiation_time_interpolation_weights

    TYPE(datetime), POINTER :: prev_radtime
    TYPE(timedelta), POINTER :: td_dt_rad

    REAL(wp) :: &
      & zaeq1(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq2(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq3(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq4(nproma,pt_patch%nlev,pt_patch%nblks_c), &
      & zaeq5(nproma,pt_patch%nlev,pt_patch%nblks_c)

    REAL(wp), ALLOCATABLE :: & ! These fields are currently only used for Kinne-Aerosol
      & od_lw(:,:,:,:)  , & !< LW optical thickness of aerosols
      & od_sw(:,:,:,:)  , & !< SW aerosol optical thickness
      & g_sw (:,:,:,:)  , & !< SW aerosol asymmetry factor
      & ssa_sw(:,:,:,:)     !< SW aerosol single scattering albedo

    CHARACTER(len=max_timedelta_str_len) :: dstring
    INTEGER :: jg
    INTEGER :: nbands_lw, nbands_sw    !< Number of short and long wave bands
    REAL(wp), POINTER :: wavenum1_sw(:), wavenum2_sw(:)
    LOGICAL :: lzacc

    REAL(wp):: zsct                    ! solar constant (at time of year)
    REAL(wp):: dsec                    ! [s] time increment of radiative transfer wrt. datetime

    LOGICAL :: is_new_day
    LOGICAL :: is_new_month

    !-------------------------------------------------

    wavenum1_sw => NULL()
    wavenum2_sw => NULL()
    ! patch ID
    jg = pt_patch%id

    IF (timers_level > 6) CALL timer_start(timer_preradiaton)

    CALL set_acc_host_or_device(lzacc, lacc)

    !-------------------------------------------------------------------------
    !  Update irradiance
    !-------------------------------------------------------------------------
    ! Only run on first radiation time of the day.

    td_dt_rad => newTimedelta('-',0,0,0,0,0, second=NINT(atm_phy_nwp_config(jg)%dt_rad), ms=0)
    prev_radtime => newDatetime(mtime_datetime + td_dt_rad)

    is_new_day = prev_radtime%date%day /= mtime_datetime%date%day
    is_new_month = is_new_day .AND. prev_radtime%date%month /= mtime_datetime%date%month

    IF (is_new_day) THEN
      IF (isolrad == 2) CALL read_bc_solar_irradiance(mtime_datetime%date%year,.TRUE.)
    END IF

    CALL deallocateTimedelta(td_dt_rad)
    CALL deallocateDatetime(prev_radtime)

    !-------------------------------------------------------------------------
    !> Radiation setup
    !-------------------------------------------------------------------------


    ! Aerosol
    !$ACC DATA CREATE(zaeq1, zaeq2, zaeq3, zaeq4, zaeq5) IF(lzacc)
    CALL nwp_aerosol_interface(mtime_datetime, pt_patch, ext_data, pt_diag, prm_diag,     &
      &                        zf(:,:,:), zh(:,:,:), dz(:,:,:),                           &
      &                        atm_phy_nwp_config(jg)%dt_rad,                             &
      &                        atm_phy_nwp_config(jg)%inwp_radiation,                     &
      &                        nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw,            &
      &                        zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                         &
      &                        od_lw, od_sw, ssa_sw, g_sw, lacc=lzacc)

    ! Ozone
    CALL o3_interface(mtime_datetime, p_sim_time, pt_patch, pt_diag, &
      &               ext_data%atm%o3, prm_diag, atm_phy_nwp_config(jg)%dt_rad, lacc=lzacc)

    IF(ANY((/irad_co2,irad_cfc11,irad_cfc12,irad_n2o,irad_ch4/) == 4)) THEN
      ! Interpolate greenhouse gas concentrations to the current date and time,
      !   placing the annual means at the mid points of the current and preceding or following year,
      !   if the current date is in the 1st or 2nd half of the year, respectively.
      ! The data file containing the greenhouse gas concentration is read in the initialisation
      !   of the NWP physics
      CALL bc_greenhouse_gases_time_interpolation( &
          & mtime_datetime, &
          & print_report=( &
          &     (msg_level >= 11 .AND. is_new_day) .OR. &
          &     (msg_level >= 5 .AND. is_new_month)) &
        )
    END IF



    ! Calculation of zenith angle optimal during dt_rad.
    ! (For radheat, actual zenith angle is calculated separately.)
    CALL pre_radiation_nwp_steps (                        &
      & kbdim        = nproma,                            & !in
      & cosmu0_dark  = cosmu0_dark,                       & !in
      & p_inc_rad    = atm_phy_nwp_config(jg)%dt_rad,     & !in
      & p_inc_radheat= atm_phy_nwp_config(jg)%dt_fastphy, & !in
      & p_sim_time   = p_sim_time,                        & !in
      & pt_patch     = pt_patch,                          & !in
      & zsmu0        = prm_diag%cosmu0(:,:),              & !out
      & zsct         = zsct,                              & !out, optional
      & lacc         = lzacc                               ) !in


    ! Compute tile-based and aggregated surface-albedo
    !
    IF ( albedo_type == MODIS ) THEN
      ! MODIS albedo
      CALL sfc_albedo_modis(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag, lzacc)
    ELSE IF ( albedo_type == 3 ) THEN
      ! globally fixed albedo value for SCM and RCEMIP applications
      CALL sfc_albedo_scm(pt_patch, albedo_fixed, prm_diag, lzacc)
    ELSE
      ! albedo based on tabulated bare soil values
      CALL sfc_albedo(pt_patch, ext_data, lnd_prog, wtr_prog, lnd_diag, prm_diag, lzacc)
    ENDIF

    IF (timers_level > 6) CALL timer_stop(timer_preradiaton)
    
    !-------------------------------------------------------------------------
    !> Radiation
    !-------------------------------------------------------------------------
    !

    SELECT CASE (atm_phy_nwp_config(jg)%inwp_radiation)
    CASE (1) ! RRTM

    
      IF ( .NOT. lredgrid ) THEN
          
        CALL nwp_rrtm_radiation ( mtime_datetime, pt_patch, ext_data, &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,        &
          & pt_diag, prm_diag, lnd_prog, lacc=.FALSE. )
       
      ELSE 

        CALL nwp_rrtm_radiation_reduced ( mtime_datetime, pt_patch,pt_par_patch, ext_data, &
          & zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                             &
          & pt_diag, prm_diag, lnd_prog, lacc=.FALSE. )
          
      ENDIF


    CASE (4) ! ecRad
      CALL finish(routine,  &
        &      'atm_phy_nwp_config(jg)%inwp_radiation = 4 needs -D__ECRAD.')

    CASE DEFAULT !Invalid inwp_radiation
      WRITE (message_text, '(a,i2,a)') 'inwp_radiation = ', atm_phy_nwp_config(jg)%inwp_radiation, &
        &                  ' not valid. Valid choices are 0: none, 1:RRTM, 4:ecRad '
      CALL finish(routine,message_text)
    END SELECT ! inwp_radiation

    !$ACC END DATA


    CALL nwp_aerosol_cleanup(od_lw, od_sw, ssa_sw, g_sw)

  END SUBROUTINE nwp_radiation


END MODULE mo_nwp_rad_interface

