! This module provides parameters controlling online emission module.
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

MODULE mo_oem_nml

    USE mo_oem_config, ONLY: config_vertical_profile_nc   => vertical_profile_nc,   &
                           & config_hour_of_day_nc        => hour_of_day_nc,        &
                           & config_day_of_week_nc        => day_of_week_nc,        &
                           & config_month_of_year_nc      => month_of_year_nc,      &
                           & config_hour_of_year_nc       => hour_of_year_nc,       &
                           & config_gridded_emissions_nc  => gridded_emissions_nc,  &
                           & config_ens_reg_nc            => ens_reg_nc,            & 
                           & config_ens_lambda_nc         => ens_lambda_nc,         &
                           & config_vegetation_indices_nc => vegetation_indices_nc, &  
                           & config_vprm_par              => vprm_par,              &
                           & config_vprm_lambda           => vprm_lambda,           &
                           & config_vprm_alpha            => vprm_alpha,            &
                           & config_vprm_beta             => vprm_beta,             &
                           & config_vprm_tmin             => vprm_tmin,             &
                           & config_vprm_tmax             => vprm_tmax,             &
                           & config_vprm_topt             => vprm_topt,             &
                           & config_vprm_tlow             => vprm_tlow,             &
                           & config_lcut_area             => lcut_area,             &  
                           & config_lon_cut_start         => lon_cut_start,         &
                           & config_lon_cut_end           => lon_cut_end,           &
                           & config_lat_cut_start         => lat_cut_start,         &
                           & config_lat_cut_end           => lat_cut_end

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: MAX_CHAR_LENGTH
  USE mo_mpi,                ONLY: my_process_is_stdio
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_io_units,           ONLY: nnml, nnml_output
  USE mo_master_control,     ONLY: use_restart_namelists
  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings
  USE mo_io_units,           ONLY: filename_max

  IMPLICIT NONE
  PRIVATE
  PUBLIC:: read_oemctrl_namelist

  !-----------------------------------
  ! namelist variables and parameters
  !-----------------------------------
  !
  CHARACTER(LEN=filename_max) :: vertical_profile_nc,   & !< name of the oae vertical profile
    &                            hour_of_day_nc,        & !< name of the oae hour of day file
    &                            day_of_week_nc,        & !< name of the oae day of week file
    &                            month_of_year_nc,      & !< name of the oae month of year file
    &                            hour_of_year_nc,       & !< name of the oae hour of year file
    &                            gridded_emissions_nc,  & !< name of the oae gridded emission file
    &                            ens_reg_nc,            & !< name of file with ensemble-regions
    &                            ens_lambda_nc,         & !< name of file with ensemble-lambdas
    &                            vegetation_indices_nc    !< name of file with MODIS reflectances 
  REAL(wp), DIMENSION(8) ::      vprm_par,              & !< VPRM parameter values for PAR_0
    &                            vprm_lambda,           & !< VPRM parameter values for lambda
    &                            vprm_alpha,            & !< VPRM parameter values for alpha
    &                            vprm_beta,             & !< VPRM parameter values for beta
    &                            vprm_tmin,             & !< VPRM parameter values for T_min
    &                            vprm_tmax,             & !< VPRM parameter values for T_max
    &                            vprm_topt,             & !< VPRM parameter values for T_opt
    &                            vprm_tlow                !< VPRM parameter values for T_low
  LOGICAL ::                     lcut_area                !< Switch to turn on/off to select an
                                                          !< area where no fluxes are applied
  REAL(wp) ::                    lon_cut_start,         & !< longitude start coordinate 
    &                            lon_cut_end,           & !< longitude end coordinate
    &                            lat_cut_start,         & !< latitude start coordinate
    &                            lat_cut_end              !< latitude end coordinate

  REAL(wp) ::                    restart_init_time

  !
  NAMELIST /oemctrl_nml/ vertical_profile_nc,   &
    &                    hour_of_day_nc,        &
    &                    day_of_week_nc,        &
    &                    month_of_year_nc,      &
    &                    hour_of_year_nc,       &
    &                    gridded_emissions_nc,  &
    &                    ens_reg_nc,            &
    &                    ens_lambda_nc,         &
    &                    vegetation_indices_nc, &
    &                    vprm_par,              &
    &                    vprm_lambda,           &
    &                    vprm_alpha,            &
    &                    vprm_beta,             &
    &                    vprm_tmin,             &
    &                    vprm_tmax,             &
    &                    vprm_topt,             &
    &                    vprm_tlow,             & 
    &                    lcut_area,             &
    &                    lon_cut_start,         & 
    &                    lon_cut_end,           & 
    &                    lat_cut_start,         & 
    &                    lat_cut_end              

!==============================================================================
! Module procedure in "mo_oem_nml"
!==============================================================================

CONTAINS


  SUBROUTINE read_oemctrl_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, iz_err
    INTEGER :: iunit

    iz_err = 0

    !0!CHARACTER(len=*), PARAMETER ::  &
    !0!  &  routine = 'mo_oem_nml:read_oemctrl_namelist'

    !-----------------------
    ! 1. default settings   
    !-----------------------

    vertical_profile_nc   = ''
    hour_of_day_nc        = ''
    day_of_week_nc        = ''
    month_of_year_nc      = ''
    hour_of_year_nc       = ''
    gridded_emissions_nc  = ''
    ens_reg_nc            = ''
    ens_lambda_nc         = ''
    vegetation_indices_nc = ''
    vprm_par              = (/  3.139508E+02,  3.132859E+02, &
                          &     5.149856E+02,  1.009878E+02, &
                          &     6.820000E+02,  1.322951E+03, &
                          &     5.794354E+02,  0.000000E+00  &
                          & /)
    vprm_lambda           = (/ -1.968102E-01, -1.808074E-01, &
                          &    -1.428753E-01, -2.019100E-01, &
                          &    -1.141000E-01, -8.061275E-02, &
                          &    -1.704107E-01,  0.000000E+00  &
                          & /)
    vprm_alpha            = (/  2.237246E-01,  1.273649E-01, &
                          &     1.714644E-01,  5.186480E-02, &
                          &     4.900000E-03,  6.764283E-02, &
                          &     8.623752E-02,  0.000000E+00  &
                          & /)
    vprm_beta             = (/ -6.411373E-01,  1.140385E+00, &
                          &     1.072367E-02, -1.675250E-01, &
                          &     0.000000E+00,  5.772793E-01, &
                          &     3.629423E-01,  0.000000E+00  &
                          & /)
    vprm_tmin             = (/  0.0,  0.0,  0.0,  2.0,  2.0,  5.0,  2.0, 0.0 /)
    vprm_tmax             = (/ 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 0.0 /)
    vprm_topt             = (/ 20.0, 20.0, 20.0, 20.0, 20.0, 22.0, 18.0, 0.0 /)
    vprm_tlow             = (/  4.0,  0.0,  2.0,  4.0,  0.0,  0.0,  0.0, 0.0 /)
    lcut_area             = .FALSE.
    lon_cut_start         = 0.0
    lon_cut_end           = 0.0
    lat_cut_start         = 0.0
    lat_cut_end           = 0.0

    !--------------------------------------------------------------------
    ! 2. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('oemctrl_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, oemctrl_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, oemctrl_nml)          ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, oemctrl_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 3. Fill the configuration state
    !----------------------------------------------------

    config_vertical_profile_nc   = vertical_profile_nc
    config_hour_of_day_nc        = hour_of_day_nc
    config_day_of_week_nc        = day_of_week_nc
    config_month_of_year_nc      = month_of_year_nc
    config_hour_of_year_nc       = hour_of_year_nc
    config_gridded_emissions_nc  = gridded_emissions_nc
    config_ens_reg_nc            = ens_reg_nc
    config_ens_lambda_nc         = ens_lambda_nc
    config_vegetation_indices_nc = vegetation_indices_nc
    config_vprm_par              = vprm_par
    config_vprm_lambda           = vprm_lambda
    config_vprm_alpha            = vprm_alpha
    config_vprm_beta             = vprm_beta
    config_vprm_tmin             = vprm_tmin
    config_vprm_tmax             = vprm_tmax
    config_vprm_topt             = vprm_topt
    config_vprm_tlow             = vprm_tlow
    config_lcut_area             = lcut_area
    config_lon_cut_start         = lon_cut_start
    config_lon_cut_end           = lon_cut_end
    config_lat_cut_start         = lat_cut_start
    config_lat_cut_end           = lat_cut_end
    
    !-----------------------------------------------------
    ! 4. write the contents of the namelist to an ASCII file
    !-----------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=oemctrl_nml)


  END SUBROUTINE read_oemctrl_namelist


!------------------------------------------------------------------------------
! End of module mo_oem_nml
!------------------------------------------------------------------------------

END MODULE mo_oem_nml


