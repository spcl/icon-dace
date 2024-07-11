!> Contains constants for the radiation processes
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
MODULE mo_rad_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                   &
   & Epar                   = 2.2E5_wp,    & ! Energy content of PAR [J / mol(photons)]=(4.6 mol/MJ PAR)**-1.  R: for update_par
   & SoilReflectivityParMin = 0.0_wp,      & ! Minimum soil reflectivity in PAR region.  R: for update_par
  !  & LaiMax = 8._wp,                       & ! Maximum LAI (used for nitrogen scaling)
   & LaiMin = 1.E-9_wp,                    & ! Min Lai in PAR computation for calc_fapar_lai_cl
   & LaiLimit = 3._wp,                     & ! Min LAI: used for nitro scaling and estimation of fract cover for calc_fapar_lai_cl
   & FcMax = 1.0_wp,                       & ! Max fractional vegetation cover. R: for calc_fapar_lai_cl
   & FcMin = 1.E-3_wp,                     & ! Min fractional vegetation cover. R: for calc_fapar_lai_cl
   & ZenithMinPar = 1.E-3_wp,              & ! Min cos of zenith angle for which PAR is calculated. R: for calc_fapar_lai_cl
   !     ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
   !     SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region
   !
   !     & minStomaConductance = 0.0_wp,        & ! Minimum stomatal conductance [mol H2O /(m^2 s) ??]
   ! Parameters for spatially constant albedo of the soil surface (UseAlbedoSoilConst)
   & AlbedoSoilConstVis = 0.2_wp,          & ! Albedo constant alternatively usable for alb_vis_mineralsoil
   & AlbedoSoilConstNir = 0.2_wp,          & ! See above, for the nir range
                                !
   & SkyViewFactor      = 0.5_wp,          & ! Constant in calculating the sky view fractions (only to scale the e-function)
   & ZenithAngleFactor  = 2._wp,           & ! Factor in solar zenith angle dependence of snow albedo
                                !                                 ! (the increase of snow albedo is the higher this factor
                                ! Parameters for albedo
   & AlbedoInitial    = 0.2_wp,            & ! initial value for albedo in the whole solar range
   & AlbedoVisInitial = 0.1_wp,            & ! initial value for albedo in the visible range
   & AlbedoNirInitial = 0.3_wp,            & ! initial value for albedo in the near infrared
   & MinRadiation     = 1.0e-09_wp,        & ! Minimum amount of radiation for albedo calculations [W/m2]
                                ! Parameters of snow age scheme albedo
   & AlbedoSnowVisMax_age   = 0.90_wp,     & ! Maximum albedo of fresh snow in the visible range
   & AlbedoSnowNirMax_age   = 0.60_wp,     & ! Maximum albedo of fresh snow in the NIR range
   & AlbedoSnowVisAge       = 0.15_wp,     & ! Maximal rel. reduction of snow albedo by aging in the visible range
   & AlbedoSnowNirAge       = 0.5_wp,      & ! Maximal rel. reduction of snow albedo by aging in the NIR range
   & AlbedoCanopySnow_age   = 0.20_wp,     & ! Albedo of snow covered canopy for snow age scheme (BATS)
   & AlbedoSnowAngle        = 0.4_wp,      & ! Maximal rel. reduction of snow absorption by large solar zenith angle
                                !
                                ! Parameters of snow temperature scheme albedo
   & AlbedoSnowVisMax_temp   = 0.90_wp,    & ! Maximum albedo of fresh snow in the visible range
   & AlbedoSnowNirMax_temp   = 0.65_wp,    & ! Maximum albedo of fresh snow in the NIR range
   & AlbedoSnowVisMin_temp   = 0.52_wp,    & ! Minimum albedo of fresh snow in the visible range
   & AlbedoSnowNirMin_temp   = 0.3_wp,     & ! Minimum albedo of fresh snow in the NIR range
   & TempAlbedoSnowMax       = 5.0_wp,     & ! Maximum snow albedo at this temperature below melting point of H2O
   & AlbedoCanopySnow_temp   = 0.20_wp,    & ! Albedo of snow covered canopy for temperature scheme
    !
    ! Parameters for lake albedo
   & AlbedoLakeWater         = 0.07_wp,    & ! Albedo of open lake water
    ! @todo: These should be set dependent on resolution (see mo_surface_ice:init_albedo_ice in ECHAM6.2)
    ! The values below correspond to the T63 settings
   & AlbedoLakeIceMin        = 0.60_wp,    &
   & AlbedoLakeIceMax        = 0.75_wp,    &
   & AlbedoLakeSnowMin       = 0.65_wp,    &
   & AlbedoLakeSnowMax       = 0.80_wp,    &
    !
    ! Parameters of glacier scheme albedo
   & AlbedoGlacierVisMin  = 0.76_wp,       & ! Albedo of glacier in the visible range at the melting point
   & AlbedoGlacierVisMax  = 0.89_wp,       & ! Albedo of glacier in the visible range at hard frost
   & TempAlbedoGlacierMax = 5.0_wp,        & ! Maximum glacier albedo at this temperature below melting point of H2O
   & AlbedoGlacierNirMin  = 0.34_wp,       & ! Albedo of glacier in the NIR range at at the melting point
   & AlbedoGlacierNirMax  = 0.72_wp,       & ! Albedo of glacier in the NIR range at hard frost
    ! Parameters for function vegetation_albedo_simple echam5
   & AlbedoCanopySnow_simple = 0.2_wp,     & ! Albedo of snow on canopy
   & AlbedoSnowMin_simple    = 0.4_wp,     & ! Minimum possible snow albedo
   & AlbedoSnowMax_simple    = 0.8_wp        ! Maximum possible snow albedo

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_rad_constants'

#endif
END MODULE mo_rad_constants
