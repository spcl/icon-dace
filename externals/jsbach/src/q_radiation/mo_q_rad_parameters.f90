!> QUINCY radiation parameters
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and init radiation parameters
!>
MODULE mo_q_rad_parameters
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  ! Parameters required for albedo and fapar calculation
  REAL(wp), SAVE :: &
    & kbl0_vis = def_parameters, &               !< extinction coefficient over black leaves
    & kdf0_vis = def_parameters, &               !< extinction coefficient for diffuse PAR (Spitters eq. 7)
    & kbl0_nir = def_parameters, &               !< extinction coefficient over black leaves in NIR range
    & kdf0_nir = def_parameters, &               !< extinction coefficient for diffuse NIR (Spitters eq. 7)
    & rho2sbeta = def_parameters, &              !< scaling factor of solar angle in reflection calculation
    & k_csf = def_parameters, &                  !< scaling factor in the crown-shape effect on clumping (Campbell & Norman, 1998, eq 15.35)
    & min_cos_zenith_angle_rad = def_parameters  !< minimum coszine zenith angle for radiation calculation

  REAL(wp), SAVE :: &
    & def_alb_vis_soil = def_parameters, &  !< default VIS albedo for soil
    & def_alb_nir_soil = def_parameters, &  !< default NIR albedo for soil
    & def_alb_vis_snow = def_parameters, &  !< default VIS albedo for snow
    & def_alb_nir_snow = def_parameters     !< default NIR albedo for snow

  ! conversion factors
  REAL(wp), SAVE :: &
    & rad2ppfd = def_parameters             !< micro-mol /J

  ! red to far red ratio calculations
  REAL(wp), SAVE :: &
    & rfr_ratio_toc = def_parameters, &     !< top of the canopy red:farred ratio  (Kull & Kruij 1998)
    & k_r2fr_chl = def_parameters           !< extinction of red:farred ratio on leaf chlorophyll (1/mmol N);
                                            !! (Kull & Krujit 1998 for the far red value, tuned to get a R:FR ratio
                                            !! of 0.15 at full canopy closure (a common value)

  REAL(wp), SAVE :: &
    & MinRadiation = def_parameters         !< Minimum amount of radiation for albedo calculations [W/m2]

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_rad_parameters'

CONTAINS

  ! ======================================================================================================= !
  !> initialize parameters for the process: q_radiation
  !>
  SUBROUTINE init_q_rad_parameters
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_q_rad_parameters'

    ! Parameters required for albedo and fapar calculation
    kbl0_vis                  = 0.5_wp        !< Spitters 1986
    kdf0_vis                  = 0.8_wp        !< Spitters 1986
    kbl0_nir                  = 0.5_wp        !< Spitters 1986
    kdf0_nir                  = 0.8_wp        !< Spitters 1986
    rho2sbeta                 = 1.6_wp        !< Spitters 1986
    k_csf                     = 2.2_wp        !< Campbell & Norman 1998, eq 15.35
    min_cos_zenith_angle_rad  = 0.1736482_wp  !< corresponding to 10 deg solar elevation, below which the
                                              !! radiation equations are no longer valid
    def_alb_vis_soil          = 0.15_wp       !< default from Bonan 2008
    def_alb_nir_soil          = 0.30_wp       !< default from Bonan 2008
    def_alb_vis_snow          = 0.6_wp        !< Literature research: from 0.2(old snow) to 0.9(new snow) | @TODO: could vary with snow age
    def_alb_nir_snow          = 0.6_wp        !< Literature research: from 0.2(old snow) to 0.9(new snow) | @TODO: could vary with snow age
    ! conversion factors
    rad2ppfd                  = 4.6_wp        !< Monteith & Unsworth 1995
    ! red to far red ratio calculations
    rfr_ratio_toc             = 1.2_wp                !< Kull & Kruijt 1998
    k_r2fr_chl                = 0.0012_wp - 0.012_wp  !< Kull & Kruijt 1998
    ! parameters also used in jsbach
    MinRadiation              = 1.0e-09_wp
  END SUBROUTINE init_q_rad_parameters

#endif
END MODULE mo_q_rad_parameters
