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

MODULE mo_radiation_solar_data
  USE mo_kind, ONLY: wp
  USE mo_impl_constants, ONLY: max_dom
  USE mo_radiation_general, ONLY: nbndsw
  PUBLIC

  REAL(wp) :: psctm(max_dom)                !< orbit and time dependent solar constant for radiation time step
  REAL(wp) :: ssi_factor(nbndsw)            !< fraction of TSI in the 14 RRTM SW bands

  REAL(wp), PARAMETER :: &
       !
       ! SRTM default solar flux
       ! sum of 14 bands is: 1.3682223735968237E+03
       ssi_default(nbndsw) =  (/                               &
       & 12.1095682699999987_wp   , 20.3650825429849398_wp   , &
       & 23.7297328613475429_wp   , 22.4276934179066849_wp   , &
       & 55.6266126199999960_wp   , 1.0293153082385436E+02_wp, &
       & 24.2936128100154392_wp   , 3.4574251380000004E+02_wp, &
       & 2.1818712729860400E+02_wp, 3.4719231470000005E+02_wp, &
       & 1.2949501812000000E+02_wp, 50.1522503011060650_wp   , &
       & 3.0799387010047838_wp    , 12.8893773299999985_wp/) , &
       !
       ! Average of 1979-1988 of transient CMIP5 solar flux
       ! sum of 14 bands is: 1361.371
       ssi_amip(nbndsw) =  (/                                             &
       & 11.95053_wp, 20.14766_wp, 23.40394_wp, 22.09458_wp, 55.41401_wp, &
       & 102.5134_wp, 24.69814_wp, 347.5362_wp, 217.2925_wp, 343.4221_wp, &
       & 129.403_wp,  47.14264_wp, 3.172126_wp, 13.18075_wp /),           &
       !
       ! Average 1844-1856 of transient CMIP5 solar flux
       ! sum of 14 bands is: 1360.875
       ssi_cmip5_picontrol(nbndsw) =  (/                                  &
       & 11.95005_wp, 20.14612_wp, 23.40302_wp, 22.09443_wp, 55.41679_wp, &
       & 102.512_wp , 24.69536_wp, 347.4719_wp, 217.2217_wp, 343.2816_wp, &
       & 129.3001_wp, 47.07624_wp, 3.130199_wp, 13.17521_wp /),           &
       !
       ! Average 1850-1873 of transient CMIP6 solar flux
       ! sum of 14 bands is: 1360.744
       ssi_cmip6_picontrol(nbndsw) =  (/                                  &
       & 12.02503_wp, 20.24537_wp, 23.69633_wp, 22.42093_wp, 55.91312_wp, &
       & 103.5685_wp, 24.46918_wp, 346.3545_wp, 217.1642_wp, 344.9984_wp, &
       & 127.7391_wp, 45.95287_wp, 2.957935_wp, 13.2384_wp /),            &
       !
       ! Solar flux for RCE simulations with diurnal cycle
       ! sum of 14 bands is: 1069.315 for diurnal cycle on.
       ssi_RCEdiurnOn(nbndsw) =  (/                                       &
       & 9.386766_wp, 15.82535_wp, 18.38306_wp,  17.3546_wp, 43.52597_wp, &
       & 80.52106_wp, 19.39961_wp, 272.9788_wp, 170.6764_wp, 269.7473_wp, &
       & 101.642_wp,  37.02906_wp, 2.491606_wp, 10.35307_wp/),            &
       !
       ! Solar flux for RCE simulations without diurnal cycle
       ! global mean insolation = 340.3 W/m2 rescaled from ssi_amip above
       ! with constant factor of app. 0.3183092
       ! sum of 14 bands is: 433.3371
       ssi_RCEdiurnOFF(nbndsw) = (/ & 
       & 3.803964_wp, 6.413186_wp, 7.44969_wp,  7.032908_wp, 17.63879_wp, &
       & 32.63096_wp, 7.861645_wp, 110.624_wp,  69.16621_wp, 109.3144_wp, &
       & 41.19017_wp, 15.00594_wp,1.009717_wp,  4.195554_wp  /),          &
       !
       ! Solar flux for RCEmip_analytical simulations without diurnal cycle
       ! global mean insolation = 409.6 W/m2 rescaled from ssi_default above
       ! with constant factor of app. 0.4031362
       ! sum of 14 bands is: 551.58
       ssi_RCEmip_analytical(nbndsw) = (/ &
       & 4.8818059_wp, 8.2099018_wp, 9.5663147_wp, 9.0414162_wp, 22.4251022_wp, &
       & 41.4954300_wp, 9.7936354_wp, 139.38133_wp, 87.9591370_wp, 139.96581_wp, &
       & 52.2041359_wp, 20.2181892_wp, 1.2416348_wp, 5.1961751_wp /)

  ! Spectral and total solar irradiance for use in integration
  REAL(wp) :: ssi_radt(nbndsw), tsi_radt, tsi

END MODULE mo_radiation_solar_data
