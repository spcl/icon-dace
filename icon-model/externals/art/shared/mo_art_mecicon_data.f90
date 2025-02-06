!
! mo_art_mecicon_data
! This module provides data key_value_store structures and constants.
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

MODULE mo_art_mecicon_data
! ICON
  USE mo_kind,                          ONLY: wp, dp, sp, i4, i8

  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
! ART
  USE messy_cmn_photol_mem,             ONLY: IP_MAX
  USE mo_physical_constants,            ONLY: argas, stbo, avo, grav
  USE mo_math_constants,                ONLY: pi
  
  
  IMPLICIT NONE

  PUBLIC 
    INTEGER , PARAMETER  :: NKPPCTRL     = 20
    INTEGER , PARAMETER  :: NMAXFIXSTEPS = 50 

    REAL(dp), PARAMETER :: TINY_DP = TINY(0._dp) 
    REAL(dp), PARAMETER :: HUGE_DP = HUGE(0._dp) 
    REAL(dp), PARAMETER :: BIG_DP = 1.0e+40_dp

    CHARACTER(LEN=*), PARAMETER :: HLINE1 = &
    '*************************************'// &
    '*************************************'
    CHARACTER(LEN=*), PARAMETER :: HLINE2 = &
    '-------------------------------------'// &
    '-------------------------------------'
    CHARACTER(LEN=*), PARAMETER :: HLINE3 = &
    '.....................................'// &
    '.....................................'
    ! mz_rs_20070904-
     
    ! STRING LENGTHs
    INTEGER, PARAMETER :: STRLEN_SHORT  = 8
    INTEGER, PARAMETER :: STRLEN_MEDIUM = 24
    INTEGER, PARAMETER :: STRLEN_LONG   = 64
    INTEGER, PARAMETER :: STRLEN_VLONG  = 80
    INTEGER, PARAMETER :: STRLEN_ULONG  = 256
    INTEGER, PARAMETER :: STRLEN_KPPSPECIES =  15

    ! PHYSICAL CONSTANTS
    REAL(dp), PARAMETER :: R_gas = argas
    REAL(dp), PARAMETER :: N_A = avo
    REAL(dp), PARAMETER :: g = grav
    REAL(dp), PARAMETER :: T0      = 298.15_dp ! standard temperature [K]
    REAL(dp), PARAMETER :: T0_INV  = 1._dp / T0      ! 1/T0 [1/K]
    REAL(dp), PARAMETER :: atm2Pa  = 101325._dp      ! conversion from [atm] to [Pa]
    REAL(dp), PARAMETER :: cal2J   = 4.1868_dp       ! conversion from [cal] to [J]
    REAL(dp), PARAMETER :: k_B     = 1.380662e-23_dp ! Boltzmann constant [J/K]
    REAL(dp), PARAMETER :: c_vKar  = 0.4_dp          !  Karman constant [?]

    ! standard atmosphere vertical gradient of the temperature in the troposphere 
    REAL(dp), PARAMETER :: stDTDZ = 0.0065_dp  ! [K/m] ! um_ak_20120103

    ! MXXX = molar mass of element XXX [g/mol]
    REAL(dp), PARAMETER :: MH  =   1.01_dp
    REAL(dp), PARAMETER :: MC  =  12.01_dp
    REAL(dp), PARAMETER :: MN  =  14.01_dp
    REAL(dp), PARAMETER :: MF  =  19.00_dp
    REAL(dp), PARAMETER :: MNa =  22.99_dp
    REAL(dp), PARAMETER :: MO  =  16.00_dp
    REAL(dp), PARAMETER :: MS  =  32.07_dp
    REAL(dp), PARAMETER :: MCl =  35.45_dp
    REAL(dp), PARAMETER :: MBr =  79.90_dp
    REAL(dp), PARAMETER :: MI  = 126.90_dp
    REAL(dp), PARAMETER :: MHg = 200.59_dp
    ! M_XXX = molar mass of compounds [g/mol]
    REAL(dp), PARAMETER :: M_O3  = MO*3._dp      ! molar mass of ozone [g/mol]
    REAL(dp), PARAMETER :: M_O2  = MO*2._dp      ! molar mass of oxygen [g/mol]
    REAL(dp), PARAMETER :: M_H2O = MH*2._dp + MO ! molar mass of H2O [g/mol]
    REAL(dp), PARAMETER :: M_N2  = MN*2._dp      ! molar mass of N2 [g/mol]

    ! DRY AIR AND WATER VAPOUR THERMODYNAMIC CONSTANTS
    REAL(dp), PARAMETER :: tmelt   = 273.15_dp    ! melting temp. of ice/snow [K]
    REAL(dp), PARAMETER :: ttrip   = 273.16_dp    ! triple point of water [K]
    REAL(dp), PARAMETER :: rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
    ! fb_mk_20101021+
    REAL(dp), PARAMETER :: rho_sea = 1025._dp     ! density of sea water in [kg/m3]
    ! fb_mk_20101021-
    REAL(dp), PARAMETER :: M_air   = 28.970_dp    ! molar mass of dry air [g/mol]
    REAL(dp), PARAMETER :: cp_air  = 1005.46_dp   ! specific heat of dry air at
                                                  ! constant pressure [J/K/kg]
    ! mz_ap_20090519+
    REAL(dp), PARAMETER :: alv   = 2.5008e6_dp    ! latent heat for vaporisation 
    !                                             ! [J/kg]
    REAL(dp), PARAMETER :: als   = 2.8345e6_dp    ! latent heat for sublimation
    !                                             ! [J/kg]
    REAL(dp), PARAMETER :: alf   = als-alv        ! latent heat for fusion [J/kg]
  
    ! mz_ap_20090519-

    ! gas constant for dry air [J/K/kg]
    REAL(dp), PARAMETER :: rd      = 1000._dp * R_gas/M_air ! 287.05_dp
    ! gas constant for water vapour
    REAL(dp), PARAMETER :: rv      = 1000._dp * R_gas/M_H2O ! 461.51_dp
    ! specific heat of water vapour at constant pressure [J/K/kg]
    REAL(dp), PARAMETER :: cpv     = 1869.46_dp
    ! dimensionless auxiliary constants
    ! op_re_20130718+
    !!$  REAL(dp), PARAMETER :: vtmpc1  = rv/rd-1.0_dp
    REAL(dp), PARAMETER :: vtmpc1  = M_air / M_H2O - 1.0_dp
    ! op_re_20130718-
    REAL(dp), PARAMETER :: vtmpc2  = cpv/cp_air-1.0_dp
    REAL(dp), PARAMETER :: MM_eps  = M_H2O/M_air ! mz_hr_20070323
 
    ! cloud and radiation
    REAL(dp), SAVE      :: ceffmin = 10.0_dp    ! min eff.radius for ice cloud
    REAL(dp), PARAMETER :: ceffmax = 150.0_dp   ! max eff.radius for ice cloud
    REAL(dp), SAVE      :: ccwmin  = 1.0e-7_dp  ! cloud water limit for cover>0
    REAL(dp), PARAMETER :: cemiss  = 0.996_dp   ! LW emissivity 

    ! PLANETARY PARAMETERS
    REAL(dp), PARAMETER :: radius_earth = 6371000.0_dp ! radius of the Earth [m]
    REAL(dp), PARAMETER :: OneDay       = 86400.0_dp   ! one day [s]
    ! fu_kk_20061002+
    REAL(dp), PARAMETER :: solc  = 1365.0_dp           ! solar constant [W/m2]
    !REAL(dp), PARAMETER :: solc  = 1365.41_dp          ! solar constant [W/m2]
    ! fu_kk_20061002-
    ! *ratio: atmospheric height/radius of the earth.
    REAL(dp), PARAMETER :: crae = 0.1277e-02_dp

    ! fb_mk_20100212+
    !!$  REAL(dp), PARAMETER :: alv   = 2.5008e6_dp ! latent heat for vaporisation 
    !!$  !                                          ! in J/kg
    !!$  REAL(dp), PARAMETER :: als   = 2.8345e6_dp ! latent heat for sublimation
    !!$  !                                          ! in J/kg
    !!$  REAL(dp), PARAMETER :: alf   = als-alv     ! latent heat for fusion in J/kg
    REAL(dp), PARAMETER :: clw   = 4186.84_dp  ! specific heat for liquid water
    !                                          ! J/K/kg
    REAL(dp), PARAMETER :: csw   = 3994._dp    ! specific heat for sea water
    !                                          ! J/K/kg
    REAL(dp), PARAMETER :: ctfreez = 271.38_dp ! temperature at which sea
                                               ! starts freezing/melting
    ! fb_mk_20100212-

    ! mz_ab_20090525+
    REAL(dp), PARAMETER:: AM    = 1.673e-27_dp  ! Atomic mass unit
    REAL(dp), PARAMETER:: ELCH  = 1.602e-19_dp  ! Electron charge

    REAL(dp), PARAMETER:: TWOPI = pi*2._dp      ! Pi*2.
    REAL(dp), PARAMETER:: PI_2  = pi*0.5_dp     ! Pi/2.
    REAL(dp), PARAMETER:: DTR   = pi/180._dp    ! Degrees to radians
    REAL(dp), PARAMETER:: RTD   = 180._dp/pi    ! Radians to degrees

    TYPE t_art_kpp_control
            ! FOR FIXED TIME STEP CONTROL
    ! ... max. number of fixed time steps (sum must be 1)
    ! ... switch for fixed time stepping
    LOGICAL            :: l_fixed_step = .FALSE.
    INTEGER            :: nfsteps = 1
    ! ... number of kpp control parameters
    !
    INTEGER,  DIMENSION(NKPPCTRL)     :: icntrl  = 0
    REAL(dp), DIMENSION(NKPPCTRL)     :: rcntrl  = 0.0_dp
    REAL(dp), DIMENSION(NMAXFIXSTEPS) :: t_steps = 0.0_dp

    INTEGER, DIMENSiON(20) :: istatus
    INTEGER                :: ierr, xNacc, xNrej
    REAL(dp)               :: time_step

  END TYPE t_art_kpp_control

  TYPE t_art_mecicon_photo
    REAL(dp), DIMENSION(IP_MAX):: jx_fjx
  END TYPE t_art_mecicon_photo

  TYPE t_art_mecicon_gen_control
    LOGICAL  :: l_aero       = .FALSE. ! switch for aero chemistry
    LOGICAL  :: l_tag        = .FALSE. ! switch for tagging
    LOGICAL  :: l_dbl        = .FALSE. ! switch for doubling
    LOGICAL  :: l_skipkpp    = .FALSE. ! skip KPP calculations?
    REAL(dp) :: zmix         = 25._dp  ! ocean mixing height [m] 
                                       ! http://en.wikipedia.org/wiki/Mixed_layer
    LOGICAL  :: l_ff         = .FALSE. ! frostflower model run?
    
  END TYPE t_art_mecicon_gen_control

  TYPE t_art_mecicon_utils
    REAL(dp), POINTER :: reac_rates(:,:,:,:)  !< this is for the diagnostic 
                                              !  output of the reaction rates (s-1)
    INTEGER, ALLOCATABLE :: mapping_indices_kpp(:)

    REAL(wp), POINTER :: &
      &  o3_column(:,:,:)          !< ozone column in DU
  END TYPE t_art_mecicon_utils

  TYPE t_art_mecicon
    TYPE(t_art_kpp_control)          :: kpp_control
    TYPE(t_art_mecicon_photo)        :: photo
    TYPE(t_art_mecicon_gen_control)  :: gen_control
    TYPE(t_art_mecicon_utils)        :: utils
  END TYPE t_art_mecicon

  END MODULE
