!
! Module to calculate source terms for different aerosol species.
!
! This module contains routines to calculate source terms for different aerosol species.
! Implemented up to now:
! - Mineral dust source term based Kok et al. (2014). This aerosol emission is then
!   converted to an AOD source term for the 2D-aerosol scheme of ICON
!
! Literature references:
! Fecan et al. (1998)   - Fecan, F., Marticorena, B., & Bergametti, G. (1998, December).
!                         Parametrization of the increase of the aeolian erosion threshold wind
!                         friction velocity due to soil moisture for arid and semi-arid areas.
!                         In Annales Geophysicae (Vol. 17, No. 1, pp. 149-157). Springer-Verlag.
! Grythe et al. (2014)  - Grythe, H., Stroem, J., Krejci, R., Quinn, P., & Stohl, A., 2014.
!                         A review of sea-spray aerosol source functions using a large global set
!                         of sea salt aerosol concentration measurements.
!                         Atmos. Chem. Phys., 14(3), 1277-1297.
! Jaegle et al. (2011)  - Jaegle, L., Quinn, P. K., Bates, T. S., Alexander, B., and Lin, J.-T., 2011:
!                         Global distribution of sea salt aerosols: new constraints from in situ
!                         and remote sensing observations
!                         Atmos. Chem. Phys., 11, 3137-3157.
! Kok et al. (2012)     - Kok, J. F., E. J. Parteli, T. I. Michaels, and D. B. Karam, 2012:
!                         The physics of wind-blown sand and dust. Rep. prog. Phys., 75(10), 106901.
! Kok et al. (2014)     - Kok, J., N. Mahowald, G. Fratini, J. Gillies, M. Ishizuka, J. Leys, M. Mikami,
!                         M.-S. Park, S.-U. Park, R. Van Pelt, et al., 2014: An improved dust emission
!                         model - part 1: Model description and comparison against measurements.
!                         Atmos. Chem. Phys., 14(23), 13023-13041.
! Kok et al. (2014b)    - Kok, J., S. Albani, N.M. Mahowald, and D.S. Ward, 2014: An improved dust
!                         emission model - Part 2: Evaluation in the Community Earth System Model,
!                         with implications for the use of dust source functions
!                         Atmos. Chem. Phys., 14, 13043-13061
! Raupach, M. R. (1993) - Dry deposition of gases and particles to vegetation.
!                         Clean Air: Journal of the Clean Air Society of Australia and New Zealand,
!                         27(4), 200.
! Rieger et al. (2017)  - Rieger D., Steiner A., Bachmann V., Gasch P., Foerstner J., Deetz K.,
!                         Vogel B., and Vogel H., 2017: Impact of the 4 April 2014 Saharan dust
!                         outbreak on the photovoltaic power generation in Germany.
!                         Atmos. Chem. and Phys. 17 (21), 13391
! Shao, Y., & Lu, H. (2000) - A simple expression for wind erosion threshold friction velocity.
!                         Journal of Geophysical Research: Atmospheres, 105(D17), 22437-22443.
! Zender et al. (2003)  - Zender, C. S., Bian, H., & Newman, D. (2003).
!                         Mineral Dust Entrainment and Deposition (DEAD) model:
!                         Description and 1990s dust climatology.
!                         Journal of Geophysical Research: Atmospheres, 108(D14).
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
MODULE mo_aerosol_sources

  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish, message, message_text
  USE mo_util_phys,                     ONLY: calc_ustar
  USE mo_util_string,                   ONLY: t_keyword_list, associate_keyword, with_keywords, int2string
  USE mo_physical_constants,            ONLY: grav
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_aerosol_sources_types,         ONLY: t_dust_source_const, t_fire_source_info
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_io_config,                     ONLY: default_read_method
  USE mo_io_units,                      ONLY: filename_max
  USE mo_run_config,                    ONLY: msg_level
  USE mo_read_interface,                ONLY: openInputFile, closeFile, on_cells, t_stream_id, read_2D, read_2D_1time
  USE mtime,                            ONLY: datetime

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_aerosol_sources'

  PUBLIC :: aerosol_dust_aod_source, aerosol_ssa_aod_source
  PUBLIC :: aerosol2d_read_data, calc_anthro_aod, calc_so4_nucleation
  PUBLIC :: inquire_fire2d_data

CONTAINS

  !>
  !! SUBROUTINE aerosol_dust_aod_source
  !!
  !! Calculates the source term for mineral dust optical depth based
  !! on external data
  !!
  SUBROUTINE aerosol_dust_aod_source (this_source, dzsoil, w_so, h_snow, w_so_ice, &
    &                                 soiltyp, plcov, lc_class, rho_a, tcm, u, v, aod_flux, dust_flux_out)

    TYPE(t_dust_source_const), INTENT(in) :: &
      &  this_source           !< Constant information for dust emission
    REAL(wp), INTENT(in)  :: &
      &  dzsoil,             & !< Thickness of uppermost soil layer (m)
      &  w_so,               & !< Soil water content (m H2O) (liquid+ice)
      &  w_so_ice,           & !< Soil ice content (m H2O)
      &  h_snow,             & !< Snow height (m)
      &  plcov,              & !< Plant cover fraction (-)
      &  rho_a,              & !< Air density in the lowermost layer (kg m-3)
      &  tcm,                & !< Transfer coefficient for momentum
      &  u,                  & !< Meridional wind speed (m s-1)
      &  v                     !< Zonal wind speed (m s-1)
    INTEGER, INTENT(in)   :: &
      &  soiltyp,            & !< Soiltype index (0-9)
      &  lc_class              !< Land use class
    REAL(wp), INTENT(out) :: &
      &  aod_flux              !< Dust optical depth tendency (s-1)
    REAL(wp), OPTIONAL, INTENT(out) :: &
      &  dust_flux_out         !< Mineral dust emission flux (kg m-2 s-1)
    ! Local Variables
    REAL(wp)              :: &
      &  dust_flux,          & !< Mineral dust emission flux (kg m-2 s-1)
      &  f_bare,             & !< Bare soil fraction, deducted from several land use classes (-)
      &  f_clay,             & !< Clay fraction (-)
      &  h_snow_fac,         & !< Snow factor, no emission with snow height > 5cm
      &  wgrav,              & !< Gravimetric soil moisture content (%)
      &  f_z0,               & !< Roughness correction factor (-)
      &  f_eta,              & !< Soil moisture correction factor (-)
      &  ustar,              & !< Friction wind speed (m s-1)
      &  ustart                !< Threshold friction wind speed (m s-1)

    ! Initializations
    dust_flux   = 0._wp

    ! Calculate gravimetric soil water content from soil moisture
    wgrav       = calc_wgrav(w_so, w_so_ice , dzsoil)
    ! Use f_clay and f_bare from look up tables
    f_clay      = this_source%f_clay(soiltyp)
    f_bare      = this_source%f_bare(lc_class)
    ! Calculate factors considering the impact of frozen soil and snow height
    h_snow_fac  = calc_hsnow_fac(h_snow)
    ! Calculate soil moisture correction
    f_eta       = calc_aerosol_dust_correta_fecan1998(wgrav, f_clay)
    ! Calculate roughness correction
    f_z0        = calc_aerosol_dust_corrz0_raupach1993(plcov)
    ! Calculate friction wind speed
    ustar       = calc_ustar(tcm, u, v)
    ! Calculate threshold friction wind speed
    ustart      = calc_aerosol_dust_ustart_shaolu2000 (f_eta, f_z0, rho_a)

    ! Calculate the emission flux
    dust_flux   = calc_aerosol_dust_mflux_kok2014 (ustar, ustart, f_bare, f_clay, rho_a)
    dust_flux   = dust_flux * h_snow_fac
    aod_flux    = calc_dust_aod (dust_flux)
    
    IF(PRESENT(dust_flux_out)) &
      &  dust_flux_out = dust_flux

  END SUBROUTINE aerosol_dust_aod_source

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_correta_fecan1998
  !!
  !! Calculates a soil moisture correction factor for the threshold friction
  !! wind speed according to Fecan et al. (1998).
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_correta_fecan1998 (wgrav, f_clay) RESULT (f_eta)

    REAL(wp), INTENT(in)    :: &
      &  wgrav,                & !< Input: Gravimetric soil moisture content       (kg kg-1)
      &  f_clay                  !< Input: Mass fraction of clay particles in soil (-)
    REAL(wp)                :: &
      &  f_eta                   !< Output: Roughness correction factor (-)
    ! Local Variables
    REAL(wp)                :: &
      &  wgravt,               & !< Gravimetric soil moisture content threshold
      &  atune                   !< Ad-hoc factor introduced by Zender et al. (2003)

    ! Range for atune highly model dependant, Kok et al. (2014b) use atune=1
    ! Zender et al. (2003) and ICON-ART (Rieger et al. 2017) use atune = 5
    atune  = 1._wp
    wgravt = atune * (14._wp * f_clay* f_clay + 17._wp * f_clay)

    ! Note that wgrav contains a soil ice content penalty
    f_eta  = SQRT( 1._wp + 1.21_wp * (MAX(0._wp,wgrav - wgravt))**0.68_wp)

  END FUNCTION calc_aerosol_dust_correta_fecan1998

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_corrz0_raupach1993
  !!
  !! Calculates a roughness correction factor for the threshold friction
  !! wind speed according to Raupach 1993.
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_corrz0_raupach1993 (plcov) RESULT(f_z0)

    REAL(wp), INTENT(in)    :: &
      &  plcov                   !< Input: Fraction of plant cover (-)
    REAL(wp)                :: &
      &  f_z0                    !< Output: Roughness correction factor
    ! Local Variables
    REAL(wp)                :: &
      &  lambda_p                !< plant cover dependant factor

    ! This function requires a check for plcov being smaller than 1 (<0.995)
    lambda_p = - 0.35_wp* LOG(1._wp - plcov)
    f_z0     = SQRT( (1._wp - 0.5_wp * lambda_p) * (1._wp + 45._wp * lambda_p) )

  END FUNCTION calc_aerosol_dust_corrz0_raupach1993

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_ustart_shaolu2000
  !!
  !! Calculates the threshold friction wind speed according to Shao and Lu (2000):
  !! Only for friction wind speeds above this threshold, saltation occurs
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_ustart_shaolu2000 (f_eta, f_z0, rho_a) RESULT(ustart)

    REAL(wp), INTENT(in)    :: &
      &  f_eta,                & !< Input: Soil moisture correction factor
      &  f_z0,                 & !< Input: Roughness correction factor
      &  rho_a                   !< Input: Air density (kg m-3)
    REAL(wp)                :: &
      &  ustart                  !< Output: Threshold friction wind speed (m s-1)
    ! Local Variables
    REAL(wp), PARAMETER     :: &
     &  D_0   = 75.e-6_wp,    & !< Saltation occurs for u_* > u_*t(D_0) (m)
     &  A_N   = 0.0123_wp,    & !< Coefficient (-)       (Shao and Lu, 2000)
     &  gamma = 3.e-4_wp,     & !< Coefficient (kg s-2)  (Shao and Lu, 2000)
     &  rho_s = 2650_wp,      & !< Bulk soil density (kg m-3)
     &  constval = A_N * ( rho_s * grav * D_0 + gamma / D_0)

    ! Use global minimum of threshold friction velocity d_min as in Rieger et al. (2017)
    !d_min = SQRT(gamma / rho_s / grav)
    ! Calculate threshold friction velocity after Shao and Lu (2000), eq. 24
    !ustart = f_eta * f_z0 * SQRT( A_N / rho_a * ( rho_s * grav * d_min + gamma / d_min) )
    ustart = f_eta * f_z0 * SQRT( constval / rho_a )

  END FUNCTION calc_aerosol_dust_ustart_shaolu2000

!-------------------------------------------------------------------------------------------------

  !>
  !! FUNCTION calc_aerosol_dust_mflux_kok2014
  !!
  !! Calculates the mineral dust emission flux according to Kok et al. (2014). 
  !!
  ELEMENTAL FUNCTION calc_aerosol_dust_mflux_kok2014 (ustar, ustart, f_bare, f_clay, rho_a) RESULT(dust_flux)

    REAL(wp), INTENT(in)  :: &
      &  ustar,              & !< Input: Friction wind speed (m s-1)
      &  ustart,             & !< Input: Threshold friction wind speed (m s-1) 
      &  f_bare,             & !< Input: Bare soil fraction (-)
      &  f_clay,             & !< Input: Clay fraction (-)
      &  rho_a                 !< Input: Air density (kg m-3)
    REAL(wp)              :: &
      &  dust_flux             !< Output: Mineral dust emission flux (kg m-2 s-1)
    ! Local Variables
    REAL(wp), PARAMETER     :: &
      &  rho_a0   = 1.225_wp,  & !< Standard atmospheric density at sea level (kg m-3)
      &  C_d0     = 3.9e-5_wp, & !< Dimensionless coeff. (4.4 +/- 0.5)*10^(-5) (Kok et al. 2014)
      &  C_e      = 2.3_wp,    & !< Dimensionless coeff.  2.0 +/- 0.3          (Kok et al. 2014)
      &  C_a      = 1.7_wp,    & !< Dimensionless coeff.  2.7 +/- 1.0          (Kok et al. 2014)
      &  ustarst0 = 0.16_wp      !< From measurements (m s-1)                  (Kok et al. 2012)
    REAL(wp)              :: &
      &  C_d,                & !< Dimensionless coefficient
      &  frac_ust_ust0,      & !< (ustarst - ustarst0) / ustarst0
      &  ustarst               !< Standardized threshold friction wind speed (m s-1)

    ! Calculate standardized threshold friction wind speed as defined in Kok et al. (2014), eq. 6.
    ustarst = ustart * SQRT(rho_a / rho_a0)

    ! Calculate the following fraction only once
    frac_ust_ust0 = (ustarst-ustarst0) / (ustarst0)

    ! Calculate C_d according to Kok et al. (2014), eq. 18b
    C_d = C_d0 * exp( -C_e * frac_ust_ust0 )

    ! Calculate dust emission flux according to Kok et al. (2014), eq. 18a
    dust_flux = C_d * f_bare * f_clay * rho_a * (ustar**2 - ustart**2) / ustarst  &
      &        * ( ustar / ustart )**(C_a * frac_ust_ust0)

    dust_flux = MAX(0._wp, dust_flux)

  END FUNCTION calc_aerosol_dust_mflux_kok2014

!-------------------------------------------------------------------------------------------------

  !>
  !! SUBROUTINE calc_dust_aod
  !!
  !! Calculates the mineral dust optical depth.
  !! There is no information on the size distribution required at this point as the size distribution
  !! of emitted dust is assumed to be constant. Hence, the weighted contribution of differently sized
  !! particles can be considered already during the preprocessing step when deriving the factor k_etot.
  !!
  !! Here, a complex refractive index of m=1.55+0.0025j was used.
  !!
  ELEMENTAL FUNCTION calc_dust_aod (dust_flux) RESULT(aod_flux)

    REAL(wp), INTENT(in)  :: &
      &  dust_flux             !< Mineral dust emission flux (kg m-2 s-1)
    REAL(wp)              :: &
      &  aod_flux              !< Output: Dust optical depth tendency (s-1)
    REAL(wp)              :: &
      &  k_etot                !< Sum of size-distribution weighted mass-specific extinction coefficient (m2 kg-1)

    ! Value calculated for a refractive index of m=1.55+0.0025j
    ! Meaningful range:
    ! For size distribution up to 2.5 microns k_etot=0.0689_wp
    ! For size distribution up to 20  microns k_etot=0.2472_wp
    k_etot = 225._wp
    aod_flux = k_etot * dust_flux

  END FUNCTION calc_dust_aod

  !>
  !! FUNCTION calc_wgrav
  !!
  !! Calculates the gravimetric soil moisture content in percent
  !!
  FUNCTION calc_wgrav(w_so, w_so_ice, dzsoil) RESULT (wgrav)

    REAL(wp), INTENT(in)  :: &
      &  w_so,               & !< Soil water+ice content (m H2O)
      &  w_so_ice,           & !< Soil ice content (m H2O)
      &  dzsoil                !< Thickness of uppermost soil layer (m)
    REAL(wp)              :: &
      &  wgrav                 !< Gravimetric soil moisture content (%)
!    REAL(wp)                 :: &
!      &  rho_s_bulk = 1.5e3_wp, & !< Bulk soil density (kg m-3)
!      &  rho_w= 1000._wp          !< Water density (kg m-3)

    !wgrav=100._wp*(w_so/dzsoil)*rho_s_bulk/rho_w
    ! Give w_so_ice a factor 10 penalty to efficiently inhibit
    ! mineral dust emission on frozen soil
    wgrav=150_wp*( (w_so + 9._wp * w_so_ice )/dzsoil)

  END FUNCTION calc_wgrav

  !>
  !! FUNCTION calc_hsnow_fac
  !!
  !! Calculate a factor to consider the impact of snow height. For a snow height 
  !! above 5 cm, the factor is 0.
  !! Use a simple linear fit in between f(h_snow=0) = 1 and f(h_snow=0.05) = 0
  !!
  ELEMENTAL FUNCTION calc_hsnow_fac(h_snow) RESULT (h_snow_fac)

    REAL(wp), INTENT(in)  :: &
      &  h_snow                !< Snow height (m)
    REAL(wp)              :: &
      &  h_snow_fac            !< Snow height factor

    h_snow_fac = -20._wp * h_snow + 1._wp

    h_snow_fac = MAX(0._wp,h_snow_fac)

  END FUNCTION calc_hsnow_fac

  !>
  !! SUBROUTINE aerosol_ssa_aod_source
  !!
  !! Calculates the source term for sea salt aerosol optical depth
  !!
  SUBROUTINE aerosol_ssa_aod_source (sst, sp_10m, aod_flux)
    REAL(wp), INTENT(in)  :: &
      &  sst,                & !< Sea surface temperature in K
      &  sp_10m                !< Wind speed in 10m height (m s-1)
    REAL(wp), INTENT(out) :: &
      &  aod_flux              !< Sea salt optical depth tendency (s-1)
    !REAL(wp), OPTIONAL, INTENT(out) :: &
    !  &  ssa_flux, ssa_flux_t  !< Sea salt emission flux (mug m-2 s-1)
    ! Local variables
    REAL(wp) ::              &
      &  sst_degc              !< Sea surface temperature in degree celcius

    sst_degc   = sst - 273.15_wp
    ! ssa_flux and ssa_flux_t calculate the sea spray aerosol mass flux. 
    ! For the 2D-aerosol implementation, only the change in AOD (aod_flux) is needed.
    ! For future checking or debugging, we leave the calls commented here.
    ! ssa_flux   = calc_ssa_mflux_grythe2014(sp_10m)
    ! ssa_flux_t = calc_ssa_sst_weighting(sst_degc) * ssa_flux
    aod_flux   = calc_ssa_sst_weighting(sst_degc) * calc_ssa_aod(sp_10m)

  END SUBROUTINE aerosol_ssa_aod_source

  !>
  !! FUNCTION calc_ssa_sst_weighting
  !!
  !! Calculates a sea surface temperature based weighting of
  !! sea salt aerosol emissions based on Jaegle et al. (2011).
  !! This temperature dependency is important to explain for
  !! the relatively high SSA concentrations found in the tropics
  !! (Grythe et al., 2014).
  !!
  ELEMENTAL FUNCTION calc_ssa_sst_weighting(sst) RESULT (temp_wgt)

    REAL(wp), INTENT(in)  :: &
      &  sst                   !< Sea surface temperature (deg C)
    REAL(wp)              :: &
      &  temp_wgt              !< Temperature weighting

    temp_wgt = 0.3_wp              &
      &      + 0.1_wp     * sst    &
      &      - 0.0076_wp  * sst**2 &
      &      + 0.00021_wp * sst**3

  END FUNCTION calc_ssa_sst_weighting

  !>
  !! SUBROUTINE calc_ssa_aod
  !!
  !! Calculates the sea salt aerosol optical depth.
  !! In order to derive mass fluxes instead of aod,
  !! the following factors can be used. 
  !! (assuming sphericity and a density of 2200 kg m-3)
  !! fac1 = 0.00023865_wp
  !! fac2 = 3.6086e-06_wp
  !!
  ELEMENTAL FUNCTION calc_ssa_aod (u10m) RESULT(aod_flux)

    REAL(wp), INTENT(in)  :: &
      &  u10m                  !< Input: Wind speed in 10m (m s-1)
    REAL(wp)              :: &
      &  aod_flux              !< Output: Sea salt optical depth tendency (s-1)
    ! Local Variables
    REAL(wp)              :: &
      &  fac1, fac2            !< Factors derived from integrating Grythe et al. 2014 source
                               !< function without the wind speed dependency times the mass
                               !< extinction coefficient over all diameters (Eq. 7 with Tw=1)
                               !< and converting mass to aod using the assumptions below

    ! Values derived assuming sphericity, a density of 2200 kg m-3 and m=1.5+1e-09j
    fac1 = 1.0438704354309787e-10_wp
    fac2 = 6.225149737200363e-13_wp

    aod_flux = fac1 * u10m**3.5_wp + fac2 * u10m**3

  END FUNCTION calc_ssa_aod

  !>
  !! SUBROUTINE inquire_fire2d_data
  !!
  !! Reads precursor data for wildfire emissions
  !!
  SUBROUTINE inquire_fire2d_data(p_patch, nroot, fire2d_filename, fire_data, mtime_datetime, lrequired)
    TYPE(t_patch), INTENT(in)    :: &
      &  p_patch
    INTEGER, INTENT(in)          :: &
      &  nroot
    CHARACTER(LEN=*), INTENT(in) :: &
      &  fire2d_filename              !< Name of file to read wildfire emission data from
    TYPE(t_fire_source_info), INTENT(inout) :: &
      &  fire_data
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime               !< Current datetime
    LOGICAL, INTENT(in)          :: &
      &  lrequired                    !< Is reading of the dataset mandatory?
    ! Local variables
    INTEGER                      :: &
      &  js
    CHARACTER(LEN=filename_max)  :: &
      &  filename                     !< Filename generated from fire2d_filename
    LOGICAL                      :: &
      &  lexist,                    & !< Does the file exist?
      &  lnewfile                     !< Was the file read before?

    DO js = 1, fire_data%nspecies
      filename = generate_filename_aerosol_data(fire2d_filename, fire_data%species(js)%varname,    &
        &                                       TRIM(p_patch%grid_filename), nroot, p_patch%level, &
        &                                       p_patch%id, mtime_datetime)
      INQUIRE (FILE=filename, EXIST=lexist)

      IF (.NOT.lexist .AND. lrequired) THEN
        message_text = 'Required wildfire data file is not found: '//TRIM(filename)
        CALL finish(TRIM(modname)//':inquire_fire2d_data',message_text)
      ENDIF

      lnewfile = (filename /= fire_data%species(js)%current_filename)
      IF (lnewfile) fire_data%species(js)%current_filename =  filename

      IF (lexist .AND. lnewfile) THEN
        IF (msg_level >= 7) THEN
          message_text = 'Updating wildfire data from file: '//TRIM(fire_data%species(js)%current_filename)
          CALL message(TRIM(modname)//':inquire_fire2d_data', message_text)
        ENDIF
        CALL aerosol2d_read_data(p_patch, fire_data%species(js)%current_filename,  &
          &                      fire_data%species(js)%varname, fire_data%species(js)%var, .TRUE.)
      ENDIF
    ENDDO

  END SUBROUTINE inquire_fire2d_data

  !>
  !! SUBROUTINE aerosol2d_read_data
  !!
  !! Reads precursor data for anthropogenic emissions
  !!
  SUBROUTINE aerosol2d_read_data(p_patch, filename, varname, var, has_time_dim)
    TYPE(t_patch), INTENT(in)    :: &
      &  p_patch
    CHARACTER(LEN=*), INTENT(in) :: &
      &  filename,                  & !< Name of file to read emission data from
      &  varname                      !< Variable name to be read
    REAL(wp), INTENT(inout)      :: &
      &  var(:,:)
    LOGICAL, INTENT(in)          :: &
      &  has_time_dim                 !< Does the variable in the dataset have a time dimension (to be ignored)
    TYPE(t_stream_id) :: stream_id

    CALL openinputfile(stream_id, filename, p_patch, default_read_method)
    
    IF (has_time_dim) THEN
      CALL read_2D_1time(stream_id, on_cells, varname, var)
    ELSE
      CALL read_2D(stream_id, on_cells, varname, var)
    ENDIF

    CALL closeFile    (stream_id)
    
  END SUBROUTINE aerosol2d_read_data

  !>
  !! Function generate_filename_aerosol_data
  !! Generates the filename for the aerosol input files based on keywords
  !!
  FUNCTION generate_filename_aerosol_data(aerosol_filename_in, species, grid_filename, nroot, jlev, idom, mtime_datetime) &
    & RESULT(result_str)

    CHARACTER(MAX_CHAR_LENGTH)          :: result_str
    CHARACTER(LEN=*), INTENT(in)        :: aerosol_filename_in
    CHARACTER(LEN=*), INTENT(in)        :: species
    CHARACTER(LEN=*), INTENT(in)        :: grid_filename
    INTEGER,          INTENT(in)        :: nroot, jlev, idom
    TYPE(datetime), POINTER, INTENT(in) :: mtime_datetime !< Current datetime
    ! Local variables
    CHARACTER(LEN=8)                    :: date_string
    TYPE (t_keyword_list), POINTER      :: keywords => NULL()

    WRITE (date_string,'(i4.4,2(i2.2))')  &
      &  mtime_datetime%date%year, mtime_datetime%date%month, mtime_datetime%date%day

    CALL associate_keyword("<species>",  TRIM(species),                    keywords)
    CALL associate_keyword("<gridfile>", TRIM(grid_filename),              keywords)
    CALL associate_keyword("<nroot>",    TRIM(int2string(nroot,"(i0)")),   keywords)
    CALL associate_keyword("<nroot0>",   TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",     TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",     TRIM(int2string(idom, "(i2.2)")), keywords)
    CALL associate_keyword("<yyyymmdd>", TRIM(date_string),                keywords)

    result_str = TRIM(with_keywords(keywords, TRIM(aerosol_filename_in)))
  END FUNCTION generate_filename_aerosol_data

  !>
  !! FUNCTION calc_anthro_aod
  !!
  !! Calc anthropogenic AOD emissions based on emission data and a scaling factor
  !!
  FUNCTION calc_anthro_aod(emi, scale_fac) RESULT(aod_flux)
    REAL(wp), INTENT(in)  :: &
      &  emi                   !< Input: Precursor for anthropogenic emissions
    REAL(wp)              :: &
      &  scale_fac             !< Species-dependent scaling factor
    REAL(wp)              :: &
      &  aod_flux              !< Output: Anthrop. emission optical depth tendency (s-1)

    aod_flux       = emi * scale_fac

  END FUNCTION calc_anthro_aod

  !>
  !! FUNCTION calc_so4_nucleation
  !!
  !! Calc pseudo-nucleation of so4 particles in remote regions. (Method is described below)
  !!
  SUBROUTINE calc_so4_nucleation(istart, iend, kstart, kend, temp, relhum, cosmu0, aod_so4)
    INTEGER, INTENT(in)   :: &
      &  istart, iend,       & !< Input: Column loop
      &  kstart, kend          !< Input: Vertical area for averaging
    REAL(wp), INTENT(in)  :: &
      &  temp(:,:),          & !< Temperature
      &  relhum(:,:),        & !< Relative humidity
      &  cosmu0(:)             !< Cosine of solar zenith angle
    REAL(wp), INTENT(inout) :: &
      &  aod_so4(:)
    ! Local variables
    REAL(wp)              :: &
      &  aod_so4_target,     & !< maximum value for AOD due to nucleation
      &  aod_fac,            & !< factor accounting for previously existing AOD
      &  ccrit(istart:iend), & !< critical so2 concentration (vertical average) (mug m-3)
      &  ccritref,           & !< Reference value for ccrit at 255K and RH=1i (mug m-3)
      &  ccrit_fac             !< factor to account for ccrit exceeding ccritref (mug m-3)
    INTEGER               :: &
      &  jc, jk

    ccritref       = calc_ccrit(255._wp,1._wp)
    aod_so4_target = 1._wp
    ccrit(:)       = 0._wp

    DO jk = kstart, kend
      DO jc = istart, iend
        ccrit(jc) = ccrit(jc) + calc_ccrit( temp(jc,jk), relhum(jc,jk) )
      ENDDO
    ENDDO

    DO jc = istart, iend
      ccrit(jc) = ccrit(jc) / real( (kend-kstart+1), wp )
      ! 2nd order polynomial with 1 at aod_so4=0, 0 at aod_so4=0.5 and minimum at aod_so4=0.5
      aod_fac   = MIN( 1._wp , MAX( 0._wp , (4._wp*aod_so4(jc)*(aod_so4(jc)-1._wp)+1._wp) ) )
      ! At 255K and RH=1, ccrit is rather low. Compute scaling factor in relation to that low value
      ! (=1 for ccritref, =>0 for large T and low RH)
      ccrit_fac = MIN( 1._wp , MAX( 0._wp , ccritref / ccrit(jc) ) )
      ! Update aod_so4, pull towards target value aod_so4_target (if higher than current concentration):
      ! As the formation of sulfuric acid from sulfur dioxide includes photocatalytic reactions
      ! Multiply with cosine**4 of solar zenith angle to reduce nucleation close to poles
      aod_so4(jc) = MAX( aod_so4(jc) , (aod_so4_target * aod_fac * ccrit_fac * cosmu0(jc)**4) )
    ENDDO

  END SUBROUTINE calc_so4_nucleation

  !>
  !! FUNCTION calc_ccrit
  !!
  !! Calc critical so2 concentration based on Kerminen & Wexler (1994)
  !!
  FUNCTION calc_ccrit(temp,relhum) RESULT(ccrit)
    REAL(wp), INTENT(in)  :: &
      &  temp,               & !< Temperature
      &  relhum                !< Relative humidity
    REAL(wp)              :: &
      &  ccrit
    ccrit = 0.16_wp * EXP(0.1_wp*temp-3.5_wp*relhum-27.7_wp)
  END FUNCTION calc_ccrit

END MODULE mo_aerosol_sources
