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

MODULE mo_radiation_config

  USE mo_kind,           ONLY: wp
  USE mo_io_units,       ONLY: filename_max
  USE mo_impl_constants, ONLY: MAX_CHAR_LENGTH, max_dom

  IMPLICIT NONE
  PUBLIC
  PUBLIC :: cams_aero_filename

  !--------------------------------------------------------------------------
  ! Basic configuration setup for radiation
  !--------------------------------------------------------------------------

  !TYPE t_radiation_config
    !
    !
    LOGICAL :: lradforcing(2) = (/.FALSE.,.FALSE./) !< diagnostic of instantaneous
    !                                               !< aerosol solar (lradforcing(1)) and
    !                                               !< thermal (lradforcing(2)) radiation forcing 
    INTEGER :: isolrad     !< mode of solar constant calculation
    !
    INTEGER :: albedo_type ! 1: albedo based on surface-type specific set of constants
                           !    (see )
                           ! 2: Modis albedo
                           ! 3: fixed albedo with value albedo_fixed

    REAL(wp) :: albedo_fixed   ! value of fixed albedo for albedo_type=3

    INTEGER :: direct_albedo   ! 1: SZA dependence according to Ritter and Geleyn (1992)
                               ! 2: limitation to diffuse albedo according to Zaengl 
                               !    applied to all land points
                               !    Ritter-Geleyn for ice 
                               ! 3: Parameterization after Yang et al (2008) for snow-free land points
                               !    limitation after Zaengl for snow-coverer points
                               !    Ritter-Geleyn implementation for ice
                               ! 4: Parameterization after Briegleb and Ramanathan (1992) for snow-free land points
                               !    limitation after Zaengl for snow-coverer points
  
    INTEGER :: direct_albedo_water ! 1: Ritter and Geleyn (1992)
                                   ! 2: Yang et al (2008)
                                   ! 3: Taylor et al (1996) as in IFS

    INTEGER :: albedo_whitecap ! 0: no whitecap albedo
                               ! 1: Seferian et al (2018) whitecaps albedo from breaking ocean waves

    INTEGER :: icld_overlap    ! method for cloud overlap calculation in shortwave part of RRTM
                               ! 1: maximum-random overlap
                               ! 2: generalized overlap (Hogan, Illingworth, 2000)
                               ! 3: maximum overlap
                               ! 4: random overlap

    INTEGER :: islope_rad(max_dom) ! slope correction for surface radiation
                                   ! 0: none
                                   ! 1: slope correction for solar radiation without shading effects
                                   ! 2: is for slope-dependent radiation with shading and skyview
                                   ! 3: slope-dependent radiation with shading without skyview

    !$ACC DECLARE CREATE(islope_rad)

    ! --- Switches for radiative agents
    !     irad_x=0 : radiation uses tracer x = 0
    !     irad_x=1 : radiation uses tracer x from a tracer variable
    !     irad_x>1 : radiation uses tracer x following external specifications of various kinds:
    !                - globally constant  or spatially varying
    !                - constant in time, constant annual cycle, or transient
    !
    INTEGER  :: irad_h2o    !< water vapor, clouds and ice for radiation
    INTEGER  :: irad_co2    !< CO2
    INTEGER  :: irad_ch4    !< CH4
    INTEGER  :: irad_n2o    !< N2O
    INTEGER  :: irad_o3     !< O3
    INTEGER  :: irad_o2     !< O2
    INTEGER  :: irad_cfc11  !< CFC 11
    INTEGER  :: irad_cfc12  !< CFC 12
    INTEGER  :: irad_aero   !< aerosols
    LOGICAL  :: lrad_aero_diag  !< diagnose aerosols
    ENUM, BIND(C)
        ENUMERATOR :: iRadAeroNone=0,         iRadAeroConst=2,      iRadAeroTegen=6,        iRadAeroCAMSclim=7, &
          &           iRadAeroCAMStd=8,       iRadAeroART=9,        iRadAeroConstKinne=12,  iRadAeroKinne=13,   &
          &           iRadAeroVolc=14,        iRadAeroKinneVolc=15, iRadAeroKinneVolcSP=18, iRadAeroKinneSP=19
    END ENUM
    !
    ! --- Name of the file that contains  dynamic greenhouse values
    !
    CHARACTER(LEN=filename_max)  :: ghg_filename
    !
    !> NetCDF file with CAMS 3D aerosols
    CHARACTER(LEN=filename_max) :: cams_aero_filename
    !
    ! --- Default gas mixing ratios - 1990 values (CMIP5)
    !
    REAL(wp) :: vmr_co2  , mmr_co2   !< CO2
    REAL(wp) :: vmr_n2o  , mmr_n2o   !< N20
    REAL(wp) :: vmr_o2   , mmr_o2    !< O2
    REAL(wp) :: vmr_ch4  , mmr_ch4   !< CH4
    REAL(wp) :: vmr_cfc11, mmr_cfc11 !< CFC 11
    REAL(wp) :: vmr_cfc12, mmr_cfc12 !< CFC 12
    !
    ! --- Different specifications of the zenith angle
    INTEGER  :: izenith           ! circular orbit, no seasonal cycle but with diurnal cycle 
    REAL(wp) :: cos_zenith_fixed  ! fixed cosine of zenith angle for izenith=6
    !
    ! --- Set minimum (pole) and maximum (equator) overlap
    !     decorrelation length scale in m for latitude-dependen function.
    REAL(wp) :: decorr_pole
    REAL(wp) :: decorr_equator
    !$ACC DECLARE CREATE(decorr_pole, decorr_equator)

    !
    ! ecRad specific configuration
    LOGICAL  :: ecrad_llw_cloud_scat    !< Do long wave cloud scattering?

    LOGICAL  :: ecrad_use_general_cloud_optics            ! Use generalized hydrometeors with different optical tables
    ! The next parameters can take different values depending on ecrad_use_general_cloud_optics (e_gen_cop)
    INTEGER  :: ecrad_iliquid_scat      !< Optical properties for liquid cloud scattering
                                        !< 0: SOCRATES (e_gen_cop =F) / Mie Droplet (e_gen_cop=T)
                                        !< 1: Slingo (1989) (only e_gen_cop=F)
    INTEGER  :: ecrad_iice_scat         !< Optical properties for ice cloud scattering
                                        !< 0: Fu et al. (both but extended radii for e_gen_cop=T)
                                        !< 1: Baran et al. (2016) (only e_gen_cop=F)
                                        !< 2: Yi et al. (2013) (only e_gen_cop=F)
                                        !< 10: Rough Fu (only only e_gen_cop=T)
                                        !< 11: Baum (only only e_gen_cop=T)     
    INTEGER  :: ecrad_isnow_scat        !< Optical properties for snow scattering (only e_gen_cop=T)
                                        !< -1: Snow is not considered independently in radiation calculation
                                        !<  0: Fu et al. 
                                        !< 10: Rough Fu 
    INTEGER  :: ecrad_igraupel_scat     !< Optical properties for graupel scattering (only e_gen_cop=T)
                                        !< -1: Graupel is not considered independently in radiation calculation
                                        !<  0: Fu et al. 
                                        !< 10: Rough Fu 
    INTEGER  :: ecrad_irain_scat        !< Optical properties for rain scattering (only e_gen_cop=T)
                                        !< -1: Rain is not considered independently in radiation calculation
                                        !<  0: Mie Rain

    INTEGER  :: ecrad_isolver           !< Radiation solver
                                        !< 0: McICA (Pincus et al. 2003)
                                        !< 1: Tripleclouds (Shonk and Hogan 2008)
                                        !< 2: McICA for OpenACC
                                        !< 3: SPARTACUS (Hogan et al. 2016)
    INTEGER  :: ecrad_igas_model        !< Gas model and spectral bands
                                        !< 0: RRTMG (Iacono et al. 2008)
                                        !< 1: ecckd (Hogan and Matricardi 2020)

    CHARACTER(len=MAX_CHAR_LENGTH) :: &
      &  ecrad_data_path                !< Folder containing optical properties
  
    ! 2.0 Non NAMELIST global variables and parameters
    ! --------------------------------

    REAL(wp) :: rad_csalbw(10) ! slope of solar albedo with respect to soil water content
                               ! as a function of depth of upper soil layer
  
    
    ! vertical profile parameters (vpp) of CH4 and N2O
    REAL(wp), PARAMETER :: vpp_ch4(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp/)
    REAL(wp), PARAMETER :: vpp_n2o(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)
    !
    !
    ! --- solar activity for radiation time step (this time step can be
    !     different from the actual integration time step, in general it
    !     is in the future relative to actual time step)
    !
    REAL(wp) :: ssi_radt(14)  !< spectrally resolved solar irradiance (SSI) 
    !                         !< [W/m2] at 1 AU distance from the sun
  
    REAL(wp) :: tsi_radt !< total solar irradiance (TSI) [W/m2]
    !                    !< at 1 AU distance from the sun
    !                    !< = SUM(ssi_radt(:))
    !
    ! --- solar activity for actual time step (for calculation of heating
    !     rates)
    !
    REAL(wp) :: tsi
    !
    ! Radiative transfer routine skips all points with cosmu0<=0. 
    ! That's why points to be skipped need to be marked with a value <=0
    REAL(wp), PARAMETER :: &
      &  cosmu0_dark = -1.e-9_wp       ! minimum cosmu0, for smaller values no shortwave calculations

    !$ACC DECLARE COPYIN(vpp_ch4, vpp_n2o)
  
  !END TYPE t_radiation_config
  !>
  !!
  !TYPE(t_radiation_config) :: radiation_config(max_dom)

END MODULE mo_radiation_config
