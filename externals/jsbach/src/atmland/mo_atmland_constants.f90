!> atmosphere-land constants, orbital parameters, conversion factors, plus init routine 
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
MODULE mo_atmland_constants
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  ! orbital parameters (not calculated - do I need to include this dynamically later?)
  REAL(wp),SAVE :: obliquity                   !< obliquity of the Earth's orbit
  REAL(wp),SAVE :: perihel                     !< perihel of the Earth's orbit
  REAL(wp),SAVE :: eccentricity                !< eccentricity of the Earth's orbit
  REAL(wp),SAVE :: solstice                    !< day of the year at which solstice occurs 
  REAL(wp),SAVE :: declination                 !< solar declination (degrees, angle sun at equator)
                                               !! not declared here, but calculated in mo_atmland_util:update_orbital_parameters
                                               !! which is called in mo_atmland_forcing:update_forcing
  
  ! solar angle calculation
  LOGICAL,SAVE  :: angle_at_GMT                !< logical to simulate solar angle according to local time or GMT 
  REAL(wp),SAVE :: reference_longitude_GMT     !< reference longitude for GMT correction 
  REAL(wp),SAVE :: min_cos_zenith_angle        !< minimum zenith angle for which PAR is calculated, see cos_zenith_angle
  
  ! solar constant
  REAL(wp),SAVE :: solar_constant              !< solar shortwave radiation (W m-2) 
  REAL(wp),SAVE :: frac_vis_swtot_srf          !< fraction of solar shortwave radiation in the VIS band (here: 320-700nm for albedo) 
  REAL(wp),SAVE :: frac_par_swvis_srf          !< fraction of visible solar shortwave radiation in the PAR band (400-700nm)
  
  ! diffuse fraction calculation
  REAL(wp),SAVE :: k1_fdiff                    !< constants for calculation of diffuse radiation
  REAL(wp),SAVE :: k2_fdiff                    !< calculation from transmission as in 
  REAL(wp),SAVE :: k3_fdiff                    !< calc_diffuse_radiation_fraction
  REAL(wp),SAVE :: k4_fdiff                    !<  
  REAL(wp),SAVE :: min_sw_fdiffuse             !< minimum diffuse fraction in shortwave downward flux
  REAL(wp),SAVE :: max_transmission_fdiff      !< atmospheric transmission rate above which the diffuse fraction becomes minimal
  
  ! diurnal temperature calculation
  REAL(wp),SAVE :: min_daylength               !< ?
  REAL(wp),SAVE :: dusk_time                   !< ?
  REAL(wp),SAVE :: tmp_peak_time               !< ?
  
  ! wind speed
  REAL(wp),SAVE :: min_wind                    !< minimum wind speed allowed to prevent total decoupling of the surface (m/s)
  REAL(wp),SAVE :: max_wind                    !< maximum wind speed allowed to prevent excess cooling of the surface (m/s)

  ! conversion factors
  REAL(wp),SAVE :: standard_press_srf          !< standard pressure in Pa
  REAL(wp),SAVE :: eps_vpd                     !< eps_vpd used for unit conversion Pa -> g g-1

  ! default values of model forcing for model testing:
  REAL(wp),SAVE :: def_t_air                   !< temperature in deg K
  REAL(wp),SAVE :: def_vpd                     !< vpd in Pa
  REAL(wp),SAVE :: def_rh                      !< relative humidity in fractions
  REAL(wp),SAVE :: def_swdown                  !< PAR in W m-2
  REAL(wp),SAVE :: def_fdiffuse                !< diffuse fraction
  REAL(wp),SAVE :: def_cos_angle               !< cosine of zenith angle
  REAL(wp),SAVE :: def_co2_mixing_ratio        !< atmospheric CO2 in ppm
  REAL(wp),SAVE :: def_co2_mixing_ratio_C13    !< atmospheric 13CO2 in ppm
  REAL(wp),SAVE :: def_co2_mixing_ratio_C14    !< atmospheric 14CO2 in pp?
  REAL(wp),SAVE :: def_co2_deltaC13            !< molar mixing ratio of 13C / 12C of atmospheric CO2
  REAL(wp),SAVE :: def_co2_deltaC14            !< molar mixing ratio of 14C / 12C of atmospheric CO2 
  REAL(wp),SAVE :: def_wind                    !< wind in m/s
  REAL(wp),SAVE :: def_n_deposition            !< default deposition of N in mg/m2/day (converted to mumol/m2/sec in mo_atmland_forcing)
  REAL(wp),SAVE :: def_p_deposition            !< default deposition of P in mg/m2/day (converted to mumol/m2/sec in mo_atmland_forcing)

  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_atmland_constants'
 
CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialise atmland constants
  !! called by jsbach_setup_models() !
  !! 
  !! initialise Earth's orbit; solar angle calculation, solar constant, diffuse fraction calculation
  !! diurnal temperature calculation, conversion factors, default values
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_atmland_constants

    USE mo_jsb_math_constants,    ONLY: hourday
    USE mo_isotope_util,          ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_atmland_constants'


    !> initialise Earth's orbit
    !!  
    obliquity                 = 23.446_wp
    perihel                   = 102.04_wp
    eccentricity              = 0.016724_wp
    solstice                  = 174.0_wp   ! 23rd of June

    !> initialise solar angle calculation
    !! 
    angle_at_GMT              = .FALSE.
    reference_longitude_GMT   = 0.0_wp
    min_cos_zenith_angle      = 0.02_wp   ! in JSBACH they use: mo_bethy_constants.f90:  REAL(dp), PARAMETER :: ZenithMinPar = 1.E-3_dp

    !> initialise solar constant
    !! 
    solar_constant            = 1370.0_wp
    frac_vis_swtot_srf        = 0.54_wp    ! UV-A+PAR / UV-A+PAR+NIR; Monteith & Unsworth 1995
    frac_par_swvis_srf        = 0.9_wp     ! PAR/PAR+UV-A

    !> initalise diffuse fraction calculation
    !! 
    k1_fdiff                  = 1.0045_wp
    k2_fdiff                  = 0.0435_wp
    k3_fdiff                  = -3.5227_wp
    k4_fdiff                  = 2.6313_wp
    min_sw_fdiffuse           = 0.1667_wp 
    max_transmission_fdiff    = 0.75_wp

    !> initialise diurnal temperature calculation
    !! 
    min_daylength             = 4._wp / hourday
    tmp_peak_time             = 14.0_wp / hourday
    dusk_time                 = 1.0_wp / hourday

    !> wind speed
    !! 
    min_wind                  = 0.3_wp
    max_wind                  = 5.0_wp

    !> conversion factors
    !! 
    standard_press_srf        = 101300.0_wp  
    eps_vpd                   = 0.622_wp  

    !> default values
    !! 
    def_t_air                 = 298.15_wp
    def_vpd                   = 1000._wp   
    def_rh                    = 0.8_wp      
    def_swdown                = 652.173913_wp
    def_fdiffuse              = 0.2_wp     
    def_cos_angle             = 1.0_wp    
    def_co2_mixing_ratio      = 380.0_wp
    def_co2_deltaC13          = -7.5_wp 
    def_co2_mixing_ratio_C13  = def_co2_mixing_ratio / (1._wp + 1._wp/calc_mixing_ratio_C13C12(def_co2_deltaC13))
    def_co2_deltaC14          = -14.9_wp
    def_co2_mixing_ratio_C14  = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14) * def_co2_mixing_ratio 
    def_wind                  = 1.0_wp         
    def_p_deposition          = 0.8_wp
    def_p_deposition          = 0.015_wp

  END SUBROUTINE init_atmland_constants

#endif
END MODULE mo_atmland_constants
