!> Contains the routines for the carbon processes
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
MODULE mo_carbon_process
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_carbon_constants, ONLY: fract_green_aboveGround, fract_wood_aboveGround, &
                               & days_per_year, sec_per_year, sec_per_day,        &
                               & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, &
                               & i_lctlib_nonsoluble, i_lctlib_humus, molarMassCO2_kg

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_Cpools, relocate_carbon_fire, relocate_carbon_damage, &
          & yasso, add_litter_to_yasso_pool, get_per_tile

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_process'

CONTAINS

  ! Calculations to allocate NPP to the plant C pools.
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE calc_Cpools( &
    & LAI,                        & ! in
    & LAI_previous,               & ! in
    & max_LAI,                    & ! in
    & NPP_pot_yDayMean,           & ! in
    !
    & fract_npp_2_woodPool,       & ! in
    & fract_npp_2_reservePool,    & ! in
    & fract_npp_2_exudates,       & ! in
    & fract_green_2_herbivory,    & ! in
    !
    & tau_c_woods,                & ! in
    !
    & LAI_shed_constant,          & ! in
    & Max_C_content_woods,        & ! in
    & specific_leaf_area_C,       & ! in
    & reserveC2leafC,             & ! in
    & pheno_type,                 & ! in
    !
    & c_green,                    & ! inout
    & c_woods,                    & ! inout
    & c_reserve,                  & ! inout
    & c_crop_harvest,             & ! inout
    !
    & soilResp_rate,              & ! out
    & NPP_flux_correction,        & ! out
    & excess_NPP,                 & ! out
    & root_exudates,              & ! out
    & cflux_c_greenwood_2_litter, & ! out
    & cflux_herbivory,            & ! out
    & cflux_herb_2_littergreen,   & ! out
    & cflux_herb_2_atm,           & ! out
    & NPP_act,                    & ! out
    !
    ! & fract_litter_wood_new,      & ! inout, optional; R: only for spitfire
    ! Yasso variables
    & temp2_30d,                  & ! in
    & precip_30d,                 & ! in
    !
    & c_acid_ag1,                 & ! inout
    & c_water_ag1,                & ! inout
    & c_ethanol_ag1,              & ! inout
    & c_nonsoluble_ag1,           & ! inout
    & c_acid_bg1,                 & ! inout
    & c_water_bg1,                & ! inout
    & c_ethanol_bg1,              & ! inout
    & c_nonsoluble_bg1,           & ! inout
    & c_humus_1,                  & ! inout
    & c_acid_ag2,                 & ! inout
    & c_water_ag2,                & ! inout
    & c_ethanol_ag2,              & ! inout
    & c_nonsoluble_ag2,           & ! inout
    & c_acid_bg2,                 & ! inout
    & c_water_bg2,                & ! inout
    & c_ethanol_bg2,              & ! inout
    & c_nonsoluble_bg2,           & ! inout
    & c_humus_2,                  & ! inout
    !
    & LeafLit_coef_acid,          & ! in
    & LeafLit_coef_water,         & ! in
    & LeafLit_coef_ethanol,       & ! in
    & LeafLit_coef_nonsoluble,    & ! in
    & LeafLit_coef_humus,         & ! in
    & WoodLit_coef_acid,          & ! in
    & WoodLit_coef_water,         & ! in
    & WoodLit_coef_ethanol,       & ! in
    & WoodLit_coef_nonsoluble,    & ! in
    & WoodLit_coef_humus,         & ! in
    & WoodLitterSize,             & ! in
    !
    & c_decomp_humus_1_sum,       & ! out, optional, necessary with yasso
    & c_decomp_humus_2_sum,       & ! out, optional, necessary with yasso
    & c_into_humus_1_sum,         & ! out, optional, necessary with yasso
    & c_into_humus_2_sum          & ! out, optional, necessary with yasso
    & )

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY : tmelt     ! tmelt =273.15_wp

    ! Arguments
    REAL(wp), INTENT(in)    :: LAI                    ! Yesterdays mean LAI
    REAL(wp), INTENT(in)    :: LAI_previous           ! The day before yesterdays mean LAI
    REAL(wp), INTENT(in)    :: max_LAI                ! Maximum value of LAI
    REAL(wp), INTENT(in)    :: NPP_pot_yDayMean       ! Yesterdays mean NPP rate [mol(C)/m^2 s]
    REAL(wp), INTENT(in)    :: fract_npp_2_woodPool    ! Fraction of NPP to be put maximally into the green pool
    REAL(wp), INTENT(in)    :: fract_npp_2_reservePool ! Fraction of NPP to be put into the optimally into the reserve pool
    REAL(wp), INTENT(in)    :: fract_npp_2_exudates    ! Fraction of NPP to be put into the optimally into the root exudates
    REAL(wp), INTENT(in)    :: fract_green_2_herbivory
    REAL(wp), INTENT(in)    :: tau_c_woods        ! Time constant by which woods Pool is depreciated  !! [days] !!
                                                 ! Note, in the lct lib it is given in [years] but when it is read in
                                                 ! by JS4 it is converted to [days].
    REAL(wp), INTENT(in)    :: LAI_shed_constant      ! Leaf shedding at a constant rate for evergreens etc. [1/day]
    REAL(wp), INTENT(in)    :: Max_C_content_woods    ! Maximum carbon content of wood pool (from lctLib)
    REAL(wp), INTENT(in)    :: specific_leaf_area_C   ! Specific leaf area (from lctLib) [m^2(leaf)/mol(Carbon)]
    REAL(wp), INTENT(in)    :: reserveC2leafC         ! Ratio of max. carbon in reserve pool to max. carbon of leaves
    INTEGER,  INTENT(in)    :: pheno_type             ! crop =5
    !
    REAL(wp), INTENT(inout) :: c_green            ! Green carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods            ! Wood carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve          ! Reserve carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_crop_harvest     ! Crop harvest carbon pool: on input last value; updated on output
                                                     !    [mol(C)/m^2(canopy)]
    !
    REAL(wp), INTENT(out)   :: soilResp_rate          ! Soil (=heterotrophic) respiration rate  [mol(C)/m^2 s] Note: this is a loss,
                                                     !    therefore negative!
    !REAL(wp), INTENT(out)   :: soilResp_rate_pot      ! soil respiration without N limitation
    REAL(wp), INTENT(out)   :: NPP_flux_correction    ! Amount by which the NPP rate entering the routine has to be corrected. This
                                                     !    correction arises either because otherwise the reserve pool would get
                                                     !    negative (positive correction), or the wood pool would exceed its
                                                     !    maximum value (negative correction). [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: excess_NPP             ! Part of NPP that could not be stored in one of the plant carbon pools
                                                     !    (green, wood, reserve), but had to be  thrown away (a posteriori
                                                     !    reduction of NPP) [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: root_exudates          ! Value of NPP_2_rootExudates had to be dropped into the litter_green_bg pool.
                                                     !    [mol(C)/m^2 s]
    REAL(wp), INTENT(out)   :: cflux_c_greenwood_2_litter       ! Total carbon flux from vegetation to litter (green+woody, ag+bg)
    REAL(wp), INTENT(out)   :: cflux_herbivory
    REAL(wp), INTENT(out)   :: cflux_herb_2_littergreen
    REAL(wp), INTENT(out)   :: cflux_herb_2_atm
    REAL(wp), INTENT(out)   :: NPP_act               ! Actual NPP after N-limitation and excess carbon drop [mol(C)/m^2 s]
    !
    ! added for thonicke fire algorithm (spitfire) needed to correctly split the carbon pools into the different fuel classes
    !real(wp),optional,intent(inout) :: fract_litter_wood_new   ! new fraction in above ground wood litter pool

    ! Meteorology for Yasso
    REAL(wp), INTENT(in)    :: temp2_30d              ! 30 day mean temperature
    REAL(wp), INTENT(in)    :: precip_30d             ! 30 day mean precipitation

    ! Yasso pools
    ! Size class 1; green litter
    !   Aboveground
    REAL(wp), INTENT(inout) :: c_acid_ag1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground
    REAL(wp), INTENT(inout) :: c_acid_bg1        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2; woody litter
    !   Aboveground
    REAL(wp), INTENT(inout) :: c_acid_ag2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   Belowground
    REAL(wp), INTENT(inout) :: c_acid_bg2        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2     ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Parameters for Yasso
       ! Vegetation dependent coefficients to separate leaf litter into classes of chemical composition
    REAL(wp), INTENT(in)    :: LeafLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp), INTENT(in)    :: LeafLit_coef_humus         ! fraction going to the humus pool
       ! Vegetation dependent coefficients to separate wood litter into classes of chemical composition
    REAL(wp), INTENT(in)    :: WoodLit_coef_acid          ! fraction going to the acid soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_water         ! fraction going to the water soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_ethanol       ! fraction going to the ethanol soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_nonsoluble    ! fraction going to the non soluble pools
    REAL(wp), INTENT(in)    :: WoodLit_coef_humus         ! fraction going to the humus pool
    REAL(wp), INTENT(in)    :: WoodLitterSize             ! size of coarse debris

    !
    REAL(wp), INTENT(inout) :: c_decomp_humus_1_sum   ! annual cflux sum of humus decomposed (leaf) [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_decomp_humus_2_sum   ! annual cflux sum of humus decomposed (wood) [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_into_humus_1_sum     ! annual cflux sum into humus (leaf) [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_into_humus_2_sum     ! annual cflux sum into humus (wood) [mol(C)/m^2(canopy)]

    ! locals
    LOGICAL  :: is_crop                 ! logical crop mask

    REAL(wp) :: c_green_max,c_reserve_optimal
    REAL(wp) :: NPP_2_greenPool, NPP_2_woodPool, NPP_2_reservePool, NPP_2_rootExudates, Green_2_herbivory
    REAL(wp) :: C_2_litter_greenPool, C_2_litter_woodPool
    REAL(wp) :: excess_carbon
    REAL(wp) :: leaf_shedding

    REAL(wp) :: c_green_pot, c_woods_pot         !! values of potential carbon pools i.e. before accounting for N-limitation
    REAL(wp) :: c_reserve_pot

    REAL(wp) :: Cflx_NPP_2_green_pot,Cflx_NPP_2_wood_pot,Cflx_NPP_2_reserve_pot
    REAL(wp) :: Cflx_2_crop_harvest, Cflx_crop_harvest_2_atmos
    REAL(wp) :: Cflx_faeces_2_LG, Cflx_faeces_2_atm
    ! R: only for spitfire
    !REAL(wp) :: litter_wood_new, litter_wood_total

    ! local variables for Yasso
    REAL(wp), DIMENSION(18) :: Yasso_io_pools
    REAL(wp), DIMENSION(2)  :: Weather
    REAL(wp), DIMENSION(18) :: Yasso_out
    REAL(wp), DIMENSION(5)  :: LeafLit_coefV
    REAL(wp), DIMENSION(5)  :: WoodLit_coefV
    REAL(wp)                :: leafLitter      ! amount of leafLitter
    REAL(wp)                :: sizeZero

    !REAL(wp) :: Nflux_NPP_2_green_pot,Nflux_NPP_2_wood_pot
    !REAL(wp) :: Nflux_litterGreen_ag_2_slow_pot, Nflux_litterGreen_bg_2_slow_pot
    !REAL(wp) :: Nflux_litter_wood_ag_2_slow_pot,Nflux_litter_wood_bg_2_slow_pot
    !REAL(wp) :: minNflux_litter_green_ag_pot,minNflux_litter_green_bg_pot
    !REAL(wp) :: minNflux_litter_wood_ag_pot,minNflux_litter_wood_bg_pot
    !REAL(wp) :: N2O_emissions_litter_pot

    !REAL(wp) :: minN_plant_demand,minN_soil_demand,total_minN_demand
    REAL(wp) :: NPP_2_green_via_mobileN, NPP_2_wood_via_mobileN
    !REAL(wp) :: Nmobile_2_green, Nmobile_2_wood, Nmobile_2_green_pot, Nmobile_2_wood_pot
    !REAL(wp) :: redFact_Nmobile
    REAL(wp) :: redFact_Nlimitation
    !REAL(wp) :: Nflx_faeces_2_SMINN, minNflux_litter_active
    !REAL(wp) :: Nflx_crop_harvest_2_SMINN

    !LOGICAL  :: with_Nitrogen        !! .FALSE.: only Carbon pools are updated; .TRUE.: in addition also Nitrogen pools are updated
    !LOGICAL  :: with_yasso           !! .FALSE.: cbalance model for litter and soil carbon is used; .TRUE.: yasso model is used
    !REAL(wp), PARAMETER :: N2O_ef_grazing = 0._wp  !! N2O grazing switched off; reasonable value for the emission factor: 0.0125_wp


    ! R: Normally parameters should be in mo_carbon_constants, but the following parameters are used only in this subroutine!
    !    Therefore I put them here. Is there any not structural advantage to put them into mo_carbon_constants?
    ! Parameter for Subroutine calc_Cpools
    REAL(wp), PARAMETER :: Q10                    = 1.8_wp          ! Empirical parameter for temperature dependence of
                                                                    !   heterotrophic respiration rate
    REAL(wp), PARAMETER :: kappa                  = 1.0_wp          ! Empirical parameter for soil moisture dependence of
                                                                    !   heterotrophic respiration rate
    REAL(wp), PARAMETER :: referenceTemp_Q10      = Tmelt           ! Reference temperature in the Q10 formula [Kelvin]
    REAL(wp), PARAMETER :: tau_c_reserve          = 1.0_wp *365._wp ! time constant by which reserve pool is depreciated [days]

    REAL(wp), PARAMETER :: tau_c_crop_harvest     = 1._wp*365._wp   ! time constant by which the crop harvest pool decays [days]
    REAL(wp), PARAMETER :: greenC2leafC           = 4.0_wp          ! ratio of carbon green pool (leaves, fine roots, starches,
                                                                    !    sugars) to leaf carbon
    REAL(wp), PARAMETER :: fract_C_faeces2_LG     = 0.3_wp          ! fraction of C from faces put into the litter green pool
    REAL(wp), PARAMETER :: fract_C_crop_harvest   = 0.5_wp          ! fraction of C from ag litter flux put to the crop harvest pool
    REAL(wp), PARAMETER :: N2O_rate_nitrification = 0.001_wp        ! JSBACH assumes lower bound of estimation in the literature as:
                                                                    !   (0.1%) due to nitrification
    REAL(wp), PARAMETER :: N2O_rate_denitrification = 0.00125_wp    ! 0.0025_wp  ! & (0.2%) due to denitrification
                                                                    ! actually 0.2% of denitrified N, not of N inputs
                                                                    ! assumption: 50% loss of NO3 by leaching
    REAL(wp), PARAMETER :: N2O_rate               = 0.005_wp        ! JSBACH assumes 0.5-1%
                                                                    !    N2O emissions rate per day (<0.1%-0.2% nitrification );
                                                                    !    0.2-4.7% due to denit in Xu Ri, 2008 * therein
    REAL(wp), PARAMETER :: sminn_NH4_fraction     = 0.4_wp          ! soil mineral NH4 is 40% of SMINN pool and rest in the form of
                                                                    !  NO3. Ref: Xu Ri, 2008
    ! Note: All fields with intent(out) are non-defined on return, even if the calling routine has set them.

    soilResp_rate            = 0.0_wp
    !soilResp_rate_pot        = 0.0_wp
    NPP_act                  = 0.0_wp
    ! redFact_Nlimitation      = 1.0_wp

    ! R: only for spitfire:
    !IF (PRESENT(fract_litter_wood_new)) THEN
    !   litter_wood_new = 0.0_wp
    !END IF

    ! Preparations
    !-------------------------------------------

    IF (pheno_type == 5) THEN
       is_crop = .TRUE.
    ELSE
       is_crop = .FALSE.
    END IF

    ! Initializations
    NPP_flux_correction        = 0.0_wp
    c_green_max                = 0.0_wp
    c_reserve_optimal          = 0.0_wp
    excess_NPP                 = 0.0_wp
    root_exudates              = 0.0_wp
    cflux_c_greenwood_2_litter = 0.0_wp
    cflux_herbivory            = 0.0_wp
    cflux_herb_2_littergreen   = 0.0_wp
    cflux_herb_2_atm           = 0.0_wp
    Cflx_2_crop_harvest        = 0.0_wp

    ! with yasso:
    Yasso_io_pools(1)  = c_acid_ag1
    Yasso_io_pools(2)  = c_water_ag1
    Yasso_io_pools(3)  = c_ethanol_ag1
    Yasso_io_pools(4)  = c_nonsoluble_ag1
    Yasso_io_pools(5)  = c_acid_bg1
    Yasso_io_pools(6)  = c_water_bg1
    Yasso_io_pools(7)  = c_ethanol_bg1
    Yasso_io_pools(8)  = c_nonsoluble_bg1
    Yasso_io_pools(9)  = c_humus_1
    Yasso_io_pools(10) = c_acid_ag2
    Yasso_io_pools(11) = c_water_ag2
    Yasso_io_pools(12) = c_ethanol_ag2
    Yasso_io_pools(13) = c_nonsoluble_ag2
    Yasso_io_pools(14) = c_acid_bg2
    Yasso_io_pools(15) = c_water_bg2
    Yasso_io_pools(16) = c_ethanol_bg2
    Yasso_io_pools(17) = c_nonsoluble_bg2
    Yasso_io_pools(18) = c_humus_2

    Weather(1) = temp2_30d
    Weather(2) = precip_30d

    LeafLit_coefV(i_lctlib_acid)       = LeafLit_coef_acid
    LeafLit_coefV(i_lctlib_water)      = LeafLit_coef_water
    LeafLit_coefV(i_lctlib_ethanol)    = LeafLit_coef_ethanol
    LeafLit_coefV(i_lctlib_nonsoluble) = LeafLit_coef_nonsoluble
    LeafLit_coefV(i_lctlib_humus)      = LeafLit_coef_humus
    WoodLit_coefV(i_lctlib_acid)       = WoodLit_coef_acid
    WoodLit_coefV(i_lctlib_water)      = WoodLit_coef_water
    WoodLit_coefV(i_lctlib_ethanol)    = WoodLit_coef_ethanol
    WoodLit_coefV(i_lctlib_nonsoluble) = WoodLit_coef_nonsoluble
    WoodLit_coefV(i_lctlib_humus)      = WoodLit_coef_humus

    Yasso_out(1:18) = 0.0_wp           ! array to store the Yasso pools/fluxes
    leafLitter      = 0.0_wp           ! amount of leafLitter
    sizeZero        = 0.0_wp

    ! on all veg tiles:

    c_green_max       = greenC2leafC * LAI / specific_leaf_area_C
    c_reserve_optimal = reserveC2leafC * max_LAI / specific_leaf_area_C

    !! ============== NON-N-LIMITED CARBON AND NITROGEN ALLOCATION =========================================================
    !!
    !! Perform those C- and N-allocation steps that are not restricted by N-availability
    !! 1. Leaf shedding
    !! 2. Wood shedding
    !! 3. Depletion of the reserve pool
    !! 4. Decomposition of slow soil pool (whether this is really independent of N-availability is debatable!!)
    !!
    !!
    !! ----------------------------------------------------------------------------
    !! 1. Leaf shedding: Transfer of C from green pool to leaf litter litter pool

!!$ tr
!!$       IF (LAI >= LAI_previous) THEN                              !! If LAI is increasing or constant the leaf shedding is given
!!$          leaf_shedding = LAI * LAI_shed_constant                 !!    by a constant loss of leaves for evergreens, raingreens
!!$                                                                  !!    and grasses.
!!$       ELSE                                                       !! Otherwise
!!$          leaf_shedding = LAI_previous - LAI                      !!    the leaf shedding is given by the decrease in LAI.
!!$       END IF
!!$ tr

    leaf_shedding = MAX(LAI * LAI_shed_constant,LAI_previous - LAI)  !! Leaf shedding is assured at a minimum constant rate
                                                                     !! which maybe exceeded by the decrease in LAI.

    C_2_litter_greenPool =  &                                     !! Maximally the whole carbon from the green pool is transferred
       MIN(greenC2leafC * leaf_shedding / specific_leaf_area_C, & !!    by leaf and root shedding to the green litter pool
           c_green)
    c_green = c_green - C_2_litter_greenPool              !! This removes carbon from the green pool

    excess_carbon = max(0.0_wp,c_green - c_green_max)               !! If by the reduction of LAI the maximum value of the green
                                                                    !!    pool is still smaller than the current value, there is
    IF (excess_carbon > 0._wp) THEN
      C_2_litter_greenPool = C_2_litter_greenPool + excess_carbon   !!    excess C that has also to be shedded to the green litter
      c_green = c_green - excess_carbon                             !!    pool and also subtracted from the green pool
    END IF

    IF (is_crop) THEN
      !!  Move shedded Carbon pro rata to crop harvest flux
      Cflx_2_crop_harvest = fract_green_aboveGround * fract_C_crop_harvest * C_2_litter_greenPool
    END IF

    cflux_c_greenwood_2_litter = C_2_litter_greenPool

    !! ---------------------------------------------------------------------------------
    !! Herbivory loss from c_green & Npool_green w.r.t PFTs in lctlib

    Green_2_herbivory  = c_green * fract_green_2_herbivory             !! Grazing flux indirect comput.from NPP losses
    c_green        = c_green - Green_2_herbivory                  !! Loss due to grazing

    Cflx_faeces_2_LG  = Green_2_herbivory * fract_C_faeces2_LG             !! Cflux_faeces_2_LG to be stored in litter green pool
                                                                           !!      (not considered in cflux_c_greenwood_2_litter)
                                                                           !!      (no life times (tau) enter Cflx_faeces_2_LG!)
    Cflx_faeces_2_atm = Green_2_herbivory * (1.0_wp - fract_C_faeces2_LG)  !! Cflux_faeces_2_atm to be lost to atmosphere
                                                                           !!      (dk: not added to CO2?)
    !! ---------------------------------------------------------------------------------
    !! 2. Wood shedding: Transfer of C from wood pool to wood litter pools
    !!
                                                            !! Assuming that forests continuously die at the inverse lifetime of trees
    C_2_litter_woodPool = c_woods / tau_c_woods             ! The bigger tau gets, the weaker the C pools will be depreciated:
                                                            ! tau_c_woods: [days] -> [timestep]
                                                            !! .. the shedded wood (MAX FUNCTION FOR NONWOODY PFTS)
    ! R: only for spitfire
    ! if (present(fract_litter_wood_new))     litter_wood_new = fract_wood_aboveGround * C_2_litter_woodPool
    c_woods = c_woods -  C_2_litter_woodPool                 !! .. and then all is subtracted from the wood pool.

    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter + C_2_litter_woodPool

    !! ----------------------------------------------------------------------------
    !! 3. Depletion of reserve pool: Transfer of C from reserve pool to green litter pool
    !!    Some organic carbon of the reserve pool is always lost with mortality of plants.
    !!    Note that the reserve pool contains no Nitrogen (starches and sugar are free of N).
    !!    Note also that the green litter pool has no fixed C/N-ratio (therefore no compensation
    !!    flux is needed)

    C_2_litter_greenPool = c_reserve / tau_c_reserve
    IF (is_crop) THEN
      Cflx_2_crop_harvest = Cflx_2_crop_harvest + &
          fract_green_aboveGround * fract_C_crop_harvest * C_2_litter_greenPool
    END IF
    c_reserve = c_reserve - C_2_litter_greenPool
    c_crop_harvest = c_crop_harvest + Cflx_2_crop_harvest

    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter + C_2_litter_greenPool

    !!
    !! ============== START OF POTENTIAL PLANT-CARBON ALLOCATION TO ESTIMATE PLANT-N-DEMAND =================================
    !!                       (exception: case NPP<0 handles actual allocation)

    c_green_pot   = c_green
    c_woods_pot   = c_woods
    c_reserve_pot = c_reserve

    !! Preliminary determination of distribution of NPP to the various pools (may be corrected afterwards according to
    !!    avaialable Nitrogen)

    IF (NPP_pot_yDayMean <= 0.0_wp) THEN                    !! In case of negative or zero NPP

      NPP_2_reservePool =  NPP_pot_yDayMean * sec_per_day   !! .. the plant looses C and
      c_reserve     = c_reserve + NPP_2_reservePool         !! .. we try to take it from the reserve pool
      IF (c_reserve < 0.0_wp) THEN                          !! But if thereby the reserve pool gets negative
        NPP_2_reservePool = NPP_2_reservePool - c_reserve   !! .. we absorb negative NPP only up to zero reserve pool and
        NPP_flux_correction = -c_reserve / sec_per_day      !! .. give the deficit as a flux correction back to the
        c_reserve = 0.0_wp                                  !! .. calling routine. Hence all Carbon from the reserve pool
      END IF                                                !! .. is used up.

      !! For correct handling of nitrogen limitation below (which handles only the case of positive NPP) set NPP to zero:
      NPP_2_reservePool = 0.0_wp
      NPP_2_greenPool   = 0.0_wp   !! All other transfer rates are zero
      NPP_2_woodPool    = 0.0_wp
      NPP_2_rootExudates= 0.0_wp
      c_reserve_pot = c_reserve    !! remember modified reserve pool

    ELSE  !! NPP is positive

      NPP_2_woodPool    = fract_npp_2_woodPool * NPP_pot_yDayMean * sec_per_day  !! NPP it is distributed according to predefined
                                                                                 !!    relative fraction
      NPP_2_reservePool = fract_npp_2_reservePool * NPP_pot_yDayMean * sec_per_day
      NPP_2_rootExudates= fract_npp_2_exudates * NPP_pot_yDayMean * sec_per_day
      NPP_2_greenPool   = (1.0_wp - fract_npp_2_woodPool - fract_npp_2_reservePool                 &
           - fract_npp_2_exudates) * NPP_pot_yDayMean * sec_per_day

      !! Growth of Wood Pool (structural carbon of living plants)

      c_woods_pot = c_woods_pot + NPP_2_woodPool                  !! Then it is attempted to put all NPP share into it
      excess_carbon = max(0.0_wp,c_woods_pot-Max_C_content_woods) !! .. but thereby the pool may get too large
      IF (excess_carbon > 0._wp) THEN
        c_woods_pot = c_woods_pot - excess_carbon                 !! .. so that the excess carbon has once more to be
                                                                  !! .. subtracted
        NPP_2_greenPool = NPP_2_greenPool + excess_carbon         !! .. and instead made available to the green pool
        NPP_2_woodPool = NPP_2_woodPool - excess_carbon           !! .. and the actual amount of NPP put to the wood pool
      END IF
                                                                  !! .. is less
      !! Growth of Reserve Pool (sugars and starches) and determination how much will enter the Green Pool
      IF (c_reserve_pot < c_reserve_optimal) THEN                 !! .. If the reserve pool is smaller than optimal,
        c_reserve_pot = c_reserve_pot +  NPP_2_reservePool        !! .... it is filled by the available NPP
        excess_carbon =                                         & !! .... Thereby it may happen that it gets larger than
             max(0.0_wp,c_reserve_pot - c_reserve_optimal)        !! .... optimal so that there is excess carbon
        IF (excess_carbon > 0._wp) THEN
          c_reserve_pot = c_reserve_pot - excess_carbon           !! .... that needs not be taken up
          NPP_2_greenPool = NPP_2_greenPool + excess_carbon       !! .... but can better be used to increase the green
          NPP_2_reservePool = NPP_2_reservePool - excess_carbon   !! .... pool and the actual amount of NPP put to the
        END IF
                                                                  !! .... reserve pool is less.
      ELSE                                                        !! .... Otherwise (reserve pool is larger than optimal)
        NPP_2_greenPool = NPP_2_greenPool + NPP_2_reservePool     !! .... all NPP is left for the green pool
        NPP_2_reservePool = 0.0_wp                                !! .... so that nothing is stored in the reserve pool.
      END IF

      !! Growth of Green Pool (leaves and fine roots): (in case of too much NPP, try to put it into the reserve pool).

      c_green_pot = c_green_pot + NPP_2_greenPool                 !! Green pool is filled by the available NPP.
      excess_carbon = max(0.0_wp,c_green_pot - c_green_max)       !! .. Thereby it may get larger than appropriate for
                                                                  !!    current LAI.
      IF (excess_carbon > 0._wp) THEN
        NPP_2_greenPool = NPP_2_greenPool - excess_carbon         !! .. Hence the actual amount of NPP put to the green
                                                                  !!    pool is less.
        c_green_pot = c_green_pot - excess_carbon                 !! .. This excess carbon needs not be taken up but
      END IF
      IF (c_reserve_pot < c_reserve_optimal) THEN                 !! .... if the reserve pool is smaller than optimal
        IF (excess_carbon > 0._wp) THEN
          c_reserve_pot = c_reserve_pot + excess_carbon           !! .... it is tried to put the carbon there,
          NPP_2_reservePool = NPP_2_reservePool + excess_carbon   !! .... which means that additional NPP is put to the
                                                                  !!      reserve pool.
        END IF
        excess_carbon =                                        &  !! .... Thereby it may happen that the reserve pool
             max(0.0_wp,c_reserve_pot - c_reserve_optimal)        !!      increases beyond the optimal value.
        IF (excess_carbon > 0._wp) THEN
          c_reserve_pot = c_reserve_pot - excess_carbon           !! .... In that case the excess carbon is once more
                                                                  !!      removed from the reserve pool
          NPP_2_reservePool = NPP_2_reservePool - excess_carbon   !! .... so that the actual amount of NPP put into the
                                                                  !!      reserve pool is less.
        END IF
      END IF

      excess_NPP = excess_carbon / sec_per_day

      root_exudates = NPP_2_rootExudates / sec_per_day

    END IF !! NPP > 0 end

    !! ============== START OF POTENTIAL LITTER-CARBON ALLOCATION ==============================================

    !! ============== use nitrogen from N mobile pool to satisfy (at least partly) for potential carbon fluxes ===
    ! Remarks: Nmobile_2_green + Nmobile_2_wood <= Npool_mobile
    ! Npool_mobile accounts minN_plant_demand by trasfering N into green and wood pool

    Cflx_NPP_2_green_pot    = NPP_2_greenPool     !! potential growth of green carbon
    Cflx_NPP_2_wood_pot     = NPP_2_woodPool      !! potential growth of wood carbon

    NPP_2_green_via_mobileN = 0.0_wp
    NPP_2_wood_via_mobileN  = 0.0_wp

    Cflx_NPP_2_green_pot = Cflx_NPP_2_green_pot - NPP_2_green_via_mobileN
    Cflx_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot  - NPP_2_wood_via_mobileN

    ! R: These lines were calculated also with_Nitrogen=.FALSE. in JS3. However, it was not used further when with_Nitrogen=.FALSE.
    !    In JS3 it should stand within the brackets.
    !Nflux_NPP_2_green_pot = Cflx_NPP_2_green_pot/cn_green
    !Nflux_NPP_2_wood_pot  = Cflx_NPP_2_wood_pot/cn_woods
    !minN_plant_demand = Nflux_NPP_2_green_pot + Nflux_NPP_2_wood_pot

    !! ============== COLLECT POTENTIAL C-FLUXES ===============================================================

    Cflx_NPP_2_reserve_pot        = NPP_2_reservePool                                     !! potential growth of reserve carbon

    redFact_Nlimitation = 1.0_wp   !! This assures absence of nitrogen limitation in this mixed C-N code

    !! ============== UPDATE ALL C-POOLS =======================================================================

    c_green       = c_green       + redFact_Nlimitation * Cflx_NPP_2_green_pot
    c_woods       = c_woods       + redFact_Nlimitation * Cflx_NPP_2_wood_pot
    c_reserve     = c_reserve     +  Cflx_NPP_2_reserve_pot        !!dk removed Nlimit, as no dependent on N

    !! ============== COMPUTE ACTUAL NPP (is needed outside the routine e.g. to recompute transpiration through stomata) =======

    IF (NPP_pot_yDayMean > 0.0_wp) THEN
      NPP_act = (NPP_2_green_via_mobileN + NPP_2_wood_via_mobileN + NPP_2_rootExudates +               &
        & redFact_Nlimitation * (Cflx_NPP_2_green_pot + Cflx_NPP_2_wood_pot) + Cflx_NPP_2_reserve_pot  &
        & ) / sec_per_day
    ELSE
      NPP_act = NPP_pot_yDayMean + NPP_flux_correction  ! NPP_pot_yDayMean is negative: NPP_flux_correction is zero
                                                                 ! except c_reserve was < 0.0_wp.
    END IF


    !! ============== COMPUTE SOIL RESPIRATION ===================================================================
    leafLitter    = cflux_c_greenwood_2_litter - C_2_litter_woodPool + Cflx_faeces_2_LG - Cflx_2_crop_harvest
    soilResp_rate = 0._wp

    !! call yasso for leaf litter; this time with nutrient effects
    CALL yasso(Yasso_io_pools(1:9), Weather, leafLitter, LeafLit_coefV,               &
      &        sizeZero, Yasso_out, fract_green_aboveGround, NPP_2_rootExudates,redFact_Nlimitation)

    c_acid_ag1        = Yasso_out(1)
    c_water_ag1       = Yasso_out(2)
    c_ethanol_ag1     = Yasso_out(3)
    c_nonsoluble_ag1  = Yasso_out(4)
    c_acid_bg1        = Yasso_out(5)
    c_water_bg1       = Yasso_out(6)
    c_ethanol_bg1     = Yasso_out(7)
    c_nonsoluble_bg1  = Yasso_out(8)
    c_humus_1         = Yasso_out(9)

    c_decomp_humus_1_sum = c_decomp_humus_1_sum + Yasso_out(17)
    c_into_humus_1_sum   = c_into_humus_1_sum + Yasso_out(18)

    ! rescale respiration (from leaf litter decomposition) from "per day" to "per second"
    soilResp_rate = Yasso_out(10) / sec_per_day

    !! call yasso for woody litter; this time with nutrient effects
    CALL yasso(Yasso_io_pools(10:18), Weather, C_2_litter_woodPool, WoodLit_coefV,       &
      & WoodLitterSize, Yasso_out, fract_wood_aboveGround, 0.0_wp, redFact_Nlimitation) ! no exudates

    c_acid_ag2        = Yasso_out(1)
    c_water_ag2       = Yasso_out(2)
    c_ethanol_ag2     = Yasso_out(3)
    c_nonsoluble_ag2  = Yasso_out(4)
    c_acid_bg2        = Yasso_out(5)
    c_water_bg2       = Yasso_out(6)
    c_ethanol_bg2     = Yasso_out(7)
    c_nonsoluble_bg2  = Yasso_out(8)
    c_humus_2         = Yasso_out(9)

    !! Yasso total respiration
    ! add respiration from wood litter decomposition and rescale respiration from "per day" to "per second"
    soilResp_rate = soilResp_rate  + Yasso_out(10) / sec_per_day

    c_decomp_humus_2_sum = c_decomp_humus_2_sum + Yasso_out(17)
    c_into_humus_2_sum   = c_into_humus_2_sum + Yasso_out(18)


    ! always:
    !! ============== Decay of crop harvest =============

    Cflx_crop_harvest_2_atmos = c_crop_harvest / tau_c_crop_harvest
    c_crop_harvest = c_crop_harvest - Cflx_crop_harvest_2_atmos
    soilResp_rate = soilResp_rate - Cflx_crop_harvest_2_atmos / sec_per_day
    !soilResp_rate_pot = soilResp_rate_pot - Cflx_crop_harvest_2_atmos / sec_per_day

    !! ============== COMPUTE Herbivory RESPIRATION =====

    cflux_herb_2_atm = -(Cflx_faeces_2_atm) /sec_per_day

    !! ============== Litter and herbivory diagnostics (flux from day-1 to s-1)
    cflux_c_greenwood_2_litter = cflux_c_greenwood_2_litter /sec_per_day
    cflux_herbivory    = Green_2_herbivory / sec_per_day
    cflux_herb_2_littergreen = Cflx_faeces_2_LG / sec_per_day   ! DSG: goes into yasso



    ! R: only for spitfire:
    ! fraction of new wood litter, assume old and new litter are respired in the same way, respiration does not change the fraction
    !        if (present(fract_litter_wood_new)) then
    !           if (with_yasso) then
    ! litter_wood_total = c_water_ag2 + c_acid_ag2 + c_ethanol_ag2 + c_nonsoluble_ag2
    !           else
    !             litter_wood_total = c_litter_wood_ag
    !           endif
    !           if ((litter_wood_total >  EPSILON(1._wp)) .AND. (litter_wood_new < litter_wood_total)) then
    !              fract_litter_wood_new = litter_wood_new / litter_wood_total
    !           else
    !              fract_litter_wood_new = 1._wp    ! fract_litter_wood_new=1 means fuel fractions will be reset to original values
    !           endif
    !        endif

  END SUBROUTINE calc_Cpools

  ! --- relocate_carbon_fire() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of fire for the carbon pools. More precisely:
  ! It is assumed that for the burned area the carbon from the above ground litter pools (c_litter_green_ag,
  ! c_litter_wood_ag) is completely released to the atmosphere and the carbon from the living plant pools
  ! (c_green_st, c_reserve_st, c_woods_st) is partly released to the atmosphere and partly relocated
  ! into the litter pools.
  !
  SUBROUTINE relocate_carbon_fire( &
    & burned_fract,                      & ! in from disturbed_frac
    & fire_fract_wood_2_atmos,           & ! in
    & LeafLit_coef,                      & ! in
    & WoodLit_coef,                      & ! in
    & c_green,                           & ! inout
    & c_woods,                           & ! inout
    & c_reserve,                         & ! inout
    & cflux_dist_green_2_soil,           & ! inout
    & cflux_dist_woods_2_soil,           & ! inout
    & c_acid_ag1,                        & ! inout
    & c_water_ag1,                       & ! inout
    & c_ethanol_ag1,                     & ! inout
    & c_nonsoluble_ag1,                  & ! inout
    & c_acid_ag2,                        & ! inout
    & c_water_ag2,                       & ! inout
    & c_ethanol_ag2,                     & ! inout
    & c_nonsoluble_ag2,                  & ! inout
    & c_acid_bg1,                        & ! inout
    & c_water_bg1,                       & ! inout
    & c_ethanol_bg1,                     & ! inout
    & c_nonsoluble_bg1,                  & ! inout
    & c_acid_bg2,                        & ! inout
    & c_water_bg2,                       & ! inout
    & c_ethanol_bg2,                     & ! inout
    & c_nonsoluble_bg2,                  & ! inout
    & c_humus_1,                         & ! inout
    & c_humus_2,                         & ! inout
    & co2flux_fire_all_2_atm             & ! out
    & )

    ! Arguments
    REAL(wp), INTENT(in)    :: burned_fract(:)              ! Fraction of the vegetated area of each tile burned till the last
                                                            ! call of this routine
    REAL(wp), INTENT(in)    :: fire_fract_wood_2_atmos      ! Fraction of above ground wood immediately emitted to atmosphere
                                                            ! by fire
    REAL(wp), INTENT(in)    :: LeafLit_coef(:)              ! Factor to spread non woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(in)    :: WoodLit_coef(:)              ! Factor to spread woody litterfall into yasso pools [ ]
    REAL(wp), INTENT(inout) :: c_green(:)                   ! Value of green carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_woods(:)                   ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_reserve(:)                 ! Value of reserve carbon pool [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(out)   :: co2flux_fire_all_2_atm(:)    ! CO2 immediately released by fire per tile [kg(CO2)m-2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_green_2_soil(:)   ! Amount of carbon relocated by wind and fire damage
                                                            ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]
    REAL(wp), INTENT(inout) :: cflux_dist_woods_2_soil(:)   ! Amount of carbon relocated by wind damage and fire
                                                            ! .. to the wood litter pools [mol(C)/m^2(canopy)s-1]
    !   YASSO
    ! SIZE CLASS 1 ; GREEN
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg1(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg1(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg1(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg1(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_1(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! SIZE CLASS 2 ; WOOD
    !   above ground C pools
    REAL(wp), INTENT(inout) :: c_acid_ag2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_ag2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_ag2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_ag2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    REAL(wp), INTENT(inout) :: c_acid_bg2(:)          ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_water_bg2(:)         ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_ethanol_bg2(:)       ! Ethanol soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_nonsoluble_bg2(:)    ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    REAL(wp), INTENT(inout) :: c_humus_2(:)           ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp)  :: non_woody_litter, fraction ! Helper
    INTEGER   :: ic, nc

    nc = SIZE(burned_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc
      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Amount of carbon released from the green and reserve pools entering the green litter pools
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + &
        & (c_green(ic) + c_reserve(ic)) * (1._wp - fract_green_aboveGround) * burned_fract(ic) / sec_per_day

      ! Amount of carbon released from the wood pool entering the woody litter pools
      cflux_dist_woods_2_soil(ic)  = cflux_dist_woods_2_soil(ic) + &
        & c_woods(ic) * (1._wp - fire_fract_wood_2_atmos * fract_wood_aboveGround) * burned_fract(ic) / sec_per_day

      ! Determin the amount of carbon released from the wood pool, the living tissue pools, and the
      ! litter pools going to the atmosphere
      co2flux_fire_all_2_atm(ic) = ( &
          ! above ground fraction of living biomass
        & non_woody_litter * fract_green_aboveGround                          &
          ! add the above ground YASSO pools
        & + c_water_ag1(ic) + c_acid_ag1(ic) + c_ethanol_ag1(ic)              &
        & + c_nonsoluble_ag1(ic)                                              &
        & + c_water_ag2(ic) + c_acid_ag2(ic) + c_ethanol_ag2(ic)              &
        & + c_nonsoluble_ag2(ic)                                              &
          ! add the fire affected fraction of the above ground wood pools
        & + c_woods(ic) * fire_fract_wood_2_atmos * fract_wood_aboveGround    &
          ! scale with the burned fraction
        & ) * burned_fract(ic) * molarMassCO2_kg / sec_per_day

      ! Reduce the above ground yasso pools accordingly
      fraction = 1._wp - burned_fract(ic)
      c_acid_ag1      (ic) = fraction * c_acid_ag1       (ic)
      c_water_ag1     (ic) = fraction * c_water_ag1      (ic)
      c_ethanol_ag1   (ic) = fraction * c_ethanol_ag1    (ic)
      c_nonsoluble_ag1(ic) = fraction * c_nonsoluble_ag1 (ic)

      c_acid_ag2      (ic) = fraction * c_acid_ag2       (ic)
      c_water_ag2     (ic) = fraction * c_water_ag2      (ic)
      c_ethanol_ag2   (ic) = fraction * c_ethanol_ag2    (ic)
      c_nonsoluble_ag2(ic) = fraction * c_nonsoluble_ag2 (ic)

      ! Add below ground fraction of the burned vegetation to the below ground green litter pools
      fraction = (1._wp - fract_green_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Distribute the remaings of burned wood to the yasso litter pools
      !  Above ground
      fraction = (1._wp - fire_fract_wood_2_atmos) * fract_wood_aboveGround * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._wp - fract_wood_aboveGround) * burned_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._wp - burned_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_fire

  ! --- relocate_carbon_damage() ---------------------------------------------------------------------------------------------------
  !
  ! This routine models the consequences of damages to the vegetation (e.g. wind break) for the carbon pools. More precisely:
  ! It is assumed that for the damaged area the carbon from the living plant pools (c_green_st, c_reserve_st and
  ! c_woods_st) is partly relocated into the litter pools.
  !
  SUBROUTINE relocate_carbon_damage(  &
    & damaged_fract,                  & ! in
    & LeafLit_coef,                   & ! in
    & WoodLit_coef,                   & ! in
    & c_green,                        & ! inout
    & c_woods,                        & ! inout
    & c_reserve,                      & ! inout
    & cflux_dist_green_2_soil,        & ! inout
    & cflux_dist_woods_2_soil,        & ! inout
    & c_acid_ag1,                     & ! inout
    & c_water_ag1,                    & ! inout
    & c_ethanol_ag1,                  & ! inout
    & c_nonsoluble_ag1,               & ! inout
    & c_acid_ag2,                     & ! inout
    & c_water_ag2,                    & ! inout
    & c_ethanol_ag2,                  & ! inout
    & c_nonsoluble_ag2,               & ! inout
    & c_acid_bg1,                     & ! inout
    & c_water_bg1,                    & ! inout
    & c_ethanol_bg1,                  & ! inout
    & c_nonsoluble_bg1,               & ! inout
    & c_acid_bg2,                     & ! inout
    & c_water_bg2,                    & ! inout
    & c_ethanol_bg2,                  & ! inout
    & c_nonsoluble_bg2,               & ! inout
    & c_humus_1,                      & ! inout
    & c_humus_2                       & ! inout
    & )

    ! Arguments
    real(wp),intent(in)    :: damaged_fract(:)         ! Fraction of the vegetated area of each tile damaged since the last
                                                       ! call of this routine
    real(wp),intent(in)    :: LeafLit_coef(:)          ! Factor to spread non woody litter into yasso pools [ ]
    real(wp),intent(in)    :: WoodLit_coef(:)          ! Factor to spread woody litterfall into yasso pools [ ]
    real(wp),intent(inout) :: c_green(:)               ! Value of green carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_woods(:)               ! Value of wood carbon pool [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_reserve(:)             ! Value of reserve carbon pool [mol(C)/m^(canopy)2]
    real(wp),intent(inout) :: cflux_dist_green_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
    real(wp),intent(inout) :: cflux_dist_woods_2_soil(:)          ! Amount of carbon relocated by wind damage
                                                                  ! (and later also by fire)
                                                                  ! .. to the green litter pools [mol(C)/m^2(canopy)s-1]

    !   YASSO
    ! Size class 1: green
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg1(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg1(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg1(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg1(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_1(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]
    ! Size class 2: wood
    !   above ground C pools
    real(wp),intent(inout) :: c_acid_ag2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_ag2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_ag2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_ag2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    !   below ground C pools
    real(wp),intent(inout) :: c_acid_bg2(:)        ! Acid soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_water_bg2(:)       ! Water soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_ethanol_bg2(:)     ! Ethanol soluble litter pool for Yasso [mol(C)/m2(canopy)]
    real(wp),intent(inout) :: c_nonsoluble_bg2(:)  ! Non-soluble litter pool for Yasso [mol(C)/m^2(canopy)]
    real(wp),intent(inout) :: c_humus_2(:)         ! Humus litter pool for Yasso [mol(C)/m^2(canopy)]

    ! Local variables
    REAL(wp)  :: non_woody_litter, fraction  ! Helper
    INTEGER   :: ic, nc

    ! Initializations
    nc = SIZE(damaged_fract)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(non_woody_litter, fraction)
    DO ic = 1, nc

      non_woody_litter = c_green(ic) + c_reserve(ic)

      ! Amount of carbon transfered from damaged wood pool to the wood litter pools
      cflux_dist_woods_2_soil(ic) = cflux_dist_woods_2_soil(ic) + c_woods(ic) * damaged_fract(ic) / sec_per_day

      ! Amount of carbon transfered from the damaged non-woody vegetation pools to the green litter pools
      cflux_dist_green_2_soil(ic) = cflux_dist_green_2_soil(ic) + non_woody_litter * damaged_fract(ic) / sec_per_day

      ! Add the carbon of damaged non-woody vegetation pools to the yasso soil pools
      !  Above ground
      fraction = fract_green_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag1(ic), c_water_ag1(ic), c_ethanol_ag1(ic), c_nonsoluble_ag1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)
      !  Below ground
      fraction = (1._wp - fract_green_aboveGround) * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg1(ic), c_water_bg1(ic), c_ethanol_bg1(ic), c_nonsoluble_bg1(ic), c_humus_1(ic), &
        &                          non_woody_litter, fraction, LeafLit_coef)

      ! Add the carbon of damaged woody vegetation pool to the yasso soil pools
      !  Above ground
      fraction = fract_wood_aboveGround * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_ag2(ic), c_water_ag2(ic), c_ethanol_ag2(ic), c_nonsoluble_ag2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)
      !  Below ground
      fraction = (1._wp - fract_wood_aboveGround)  * damaged_fract(ic)
      CALL distribute_yasso_litter(c_acid_bg2(ic), c_water_bg2(ic), c_ethanol_bg2(ic), c_nonsoluble_bg2(ic), c_humus_2(ic), &
        &                          c_woods(ic), fraction, WoodLit_coef)

      ! Reduce the living plant carbon pools accordingly
      fraction = 1._wp - damaged_fract(ic)
      c_green(ic)   = c_green(ic)   * fraction
      c_woods(ic)   = c_woods(ic)   * fraction
      c_reserve(ic) = c_reserve(ic) * fraction

    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE relocate_carbon_damage

#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE yasso (Yasso_io_pools, Weather, litter, Lit_coefV,WoodLitterSize, Yasso_out, fract_aboveground, &
                         NPP_2_rootExudates, redFact_Nlimit)
    !! DSG: 11.03.2013
    !! modified Yasso soil C model using an Euler scheme and readable coding style
    !! AWEN (Acid, Water, Ethanol, Nonsoluble) pools are the yasso pools without the humus pools.
    !! JSBACH specific modification of model structure: AWEN pools are separated into aboveground and belowground pools
    !!
    !! Herbivory flux is treated like normal litter, because the herbivory fluxes in JSBACH are not evaluated, yet.
    !! If the represenation of herbivory in JSBACH is evaluated I recommend to treat herbivory flux in YASSO different from general
    !! litter.
    !!
    !! The routine needs as input 1. Yasso pools
    !!                            2. Weather: precipitation and Temperature
    !!                            3. daily litter flux and NPP_2_rootExudates
    !!                            4. PFT specific parameters
    !!
    !! Even when this subroutine is called every day, the dimensions of the variables refere to years!
    !! output pools: AWEN + humus
    !! The routine must be called twice, first for non-woody litter, then for woody litter
    !!-------------------------------------------------------------------------------------------------------------------------
    !! Yasso - Soil carbon model (Jari Liski)
    !!
    !! The code for the yasso model has been developed and is provided by the Finnish Environment Institute SYKE. It is
    !! distributed as part of the ICON Earth System Model under the license conditions as stated in the header of this module.
    !!
    !! model to calculate the amount of soil organic carbon, changes in the amount of soil organic carbon and
    !! heterotrophic soil respiration.
    !! For documention of yasso see Liski et al 2005, Tuomi et al 2008,2009,2011
    !!
    !! implementation: Tea Thum (FMI), Petri Raisanen (FMI), Daniel Goll (MPI-M)
    !!-------------------------------------------------------------------------------------------------------------------------

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: tmelt
    IMPLICIT NONE

    REAL(wp), DIMENSION(9),  INTENT(IN)  :: Yasso_io_pools     ! Yasso pools IN [mol(c)/m2]
    REAL(wp),                INTENT(IN)  :: WoodLitterSize     ! Size of woody litter;
                                                               !   is zero when the subroutine is called for non-woody litter
    REAL(wp), DIMENSION(2),  INTENT(IN)  :: Weather            ! climatic drivers: air temperature and precipitation, 15-d means
                                                               ! Weather(1) in [K] ; Weather(2) in [mm/s]
    REAL(wp),                INTENT(IN)  :: litter             ! fresh litter inputs  (above and belowground) [mol(C)/m2/day]
    REAL(wp), DIMENSION(5),  INTENT(IN)  :: Lit_coefV          ! fractions to seperate incoming litter to the yasso pool    [ ]
    REAL(wp),                INTENT(IN)  :: fract_aboveground  ! parameter to seperate litter influx into above- and belowground
                                                               !     part                                                   [ ]
    REAL(wp),                INTENT(IN)  :: NPP_2_rootExudates ! root exudates [daily]  [mol(C)/m2/day]

    REAL(wp), DIMENSION(18), INTENT(OUT) :: Yasso_out          ! updated Y pools(1:9) [mol(c)/m2] & respiration(10) [mol(c)/m2/d]
                                                               ! & Cflx_2_humus (11) [mol(c)/m2/d] ... see update_Cpools for more
                                                               ! information

    REAL(wp), OPTIONAL, INTENT(IN)       ::  redFact_Nlimit

    !local
    REAL(wp)     :: Cflx_litter_2_humus            ! Litter influx to humus soluble pool

    ! Yasso pools   [mol(c)/m2]
    REAL(wp)     :: c_acid_ag    ! Aboveground
    REAL(wp)     :: c_water_ag
    REAL(wp)     :: c_ethanol_ag
    REAL(wp)     :: c_nonsoluble_ag
    REAL(wp)     :: c_acid_bg    ! Belowground
    REAL(wp)     :: c_water_bg
    REAL(wp)     :: c_ethanol_bg
    REAL(wp)     :: c_nonsoluble_bg

    REAL(wp)     :: c_humus

    ! Aboveground
    REAL(wp)     :: Cflx_from_acid_ag              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_water_ag             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_ethanol_ag           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_nonsoluble_ag        ! Loss flux of non soluble pool         [mol(c)/m2/a]

    REAL(wp)     :: Cflx_2_acid_ag                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_water_ag                ! Gain flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_ethanol_ag              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_nonsoluble_ag           ! Gain flux of non soluble pool         [mol(c)/m2/a]

    ! Belowground
    REAL(wp)     :: Cflx_from_acid_bg              ! Loss flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_water_bg             ! Loss flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_ethanol_bg           ! Loss flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_nonsoluble_bg        ! Loss flux of non soluble pool         [mol(c)/m2/a]
    REAL(wp)     :: Cflx_from_humus                ! Loss flux of humus pool               [mol(c)/m2/a]

    REAL(wp)     :: Cflx_2_acid_bg                 ! Gain flux of acid soluble pool        [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_water_bg                ! Gain flux of water soluble pool       [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_ethanol_bg              ! Gain flux of ethanol soluble pool     [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_nonsoluble_bg           ! Gain flux of non soluble pool         [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humus                   ! Gain flux of humus pool               [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humusAG                 ! Gain flux of humus pool; from ag only [mol(c)/m2/a]
    REAL(wp)     :: Cflx_2_humusBG                 ! Gain flux of humus pool; from bg only [mol(c)/m2/a]

    ! Respiration is output of Yasso therefore a daily flux
    REAL(wp)     :: soilResp_rateYasso             ! Flux to the atmosphere                [mol(c)/m2/d]
    REAL(wp)     :: soilResp_rateLitterAG          ! Flux to the atmosphere; from AWEN ag only  [mol(c)/m2/d]
    REAL(wp)     :: soilResp_rateLitterBG          ! Flux to the atmosphere; from AWEN bg only  [mol(c)/m2/d]

    ! Decomposition rates [1/a]
    REAL(wp)     :: d_acid
    REAL(wp)     :: d_water
    REAL(wp)     :: d_ethanol
    REAL(wp)     :: d_nonsoluble
    REAL(wp)     :: d_humus

    ! Dcomposition rate for N litter green (needed for N cycle)
    REAL(wp)     :: d_litter_green              ! the rate with which the sum of AWEN pools decays
    REAL(wp)     :: Pseudo_litter_green         ! Sum of all AWEN pools

    ! Yasso parameters (Tuomi '11)
    ! reference decompositation rates for T = 0C and unlimited water availability  in [1/yr]...
    REAL(wp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_wp    !  ... acid-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_water      = -5.9_wp     !  ... water-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_wp    !  ... ethanol-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_wp   !  ... nonsoluble-soluble fraction
    REAL(wp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_wp  !  ... humus
    ! parameters for temperature dependence
    REAL(wp), PARAMETER   :: temp_p1                    = 0.095_wp    ! first  order temperature dependence [1/C]
    REAL(wp), PARAMETER   :: temp_p2                    = -0.0014_wp  ! second order temperature dependence [1/C]
    ! parameter for precipitation dependence
    REAL(wp), PARAMETER   :: precip_p1                  = -1.21_wp    ! first  order precipitation dependence [a/m]
    ! parameters for size dependence
    REAL(wp), PARAMETER   :: size_p1                    = -1.71_wp    ! first  order size dependence [1/cm]
    REAL(wp), PARAMETER   :: size_p2                    = 0.86_wp     ! second order size dependence [1/cm]
    REAL(wp), PARAMETER   :: size_p3                    = -0.306_wp   ! size dependence power [ ]

    ! Fractions to spit the decomposition flux of each pool into fluxes to the other AWEN pools and atmosphere
    ! REMARK: there is no a priori (to the calibration procedure of YASSO) assumption about possible fluxes between pools,
    ! therefore all fluxes are theoretically possible, however some of the fraction are zero which means the corresponding flux
    ! does not exists. However, after a recalibration of the model this fraction may change.
    !
    ! fractions of decomposition flux of acid soluble pool going to the other AWEN pools [ ]
    ! (8,11,14)
    REAL(wp), PARAMETER   :: A_2_W                      =  0.99_wp
    REAL(wp), PARAMETER   :: A_2_E                      =  0.00_wp
    REAL(wp), PARAMETER   :: A_2_N                      =  0.00_wp
    ! fractions of decomposition flux of water soluble pool going to the other AWEN pools [ ]
    ! (5,12,15)
    REAL(wp), PARAMETER   :: W_2_A                      = 0.48_wp
    REAL(wp), PARAMETER   :: W_2_E                      = 0.00_wp
    REAL(wp), PARAMETER   :: W_2_N                      = 0.015_wp
    ! fractions of decomposition flux of ethanol soluble pool going to the other AWEN pools [ ]
    ! (6,9,16)
    REAL(wp), PARAMETER   :: E_2_A                      = 0.01_wp
    REAL(wp), PARAMETER   :: E_2_W                      = 0.00_wp
    REAL(wp), PARAMETER   :: E_2_N                      = 0.95_wp
    ! fractions of decomposition flux of non soluble pool going to the other AWEN pools [ ]
    ! (7,10,13)
    REAL(wp), PARAMETER   :: N_2_A                      = 0.83_wp
    REAL(wp), PARAMETER   :: N_2_W                      = 0.01_wp
    REAL(wp), PARAMETER   :: N_2_E                      = 0.02_wp
    ! fraction of decomposition fluxes of AWEN pools which enters the humus pool [ ]
    REAL(wp), PARAMETER   :: AWEN_2_H                   = 0.0045_wp

    ! climatic drivers
    REAL(wp)     :: precip                  ! 15 runnning mean of precipitation      [m/a]
    REAL(wp)     :: temp                    ! 15 runnning mean of 2m air temperature [C]

    REAL(wp)     :: d_temp                  ! term which accounts for the temperature influence on litter decomposition
    REAL(wp)     :: d_precip                ! term which accounts for the precipitation influence on litter decomposition
    REAL(wp)     :: d_size                  ! term which accounts for the litter size influence on litter decomposition

    REAL(wp)     :: redFactor  ! dummy for redFact_Nlimit

    ! time stepping; the parameterization of yasso is done on a annual time
    ! step, thus the fluxes have to be scaled up to annual rates.
    REAL(wp)     :: dt
    dt = 1.0_wp/days_per_year

    IF(PRESENT(redFact_Nlimit)) THEN ! nitrogen cycle active & we know already the Nlimitation factor
        redFactor  = redFact_Nlimit
    ELSE
        redFactor  = 1.0_wp
    ENDIF

   ! 0. Preparations
    ! initialize yasso pools (AWEN + humus)
    c_acid_ag        = Yasso_io_pools(1)
    c_water_ag       = Yasso_io_pools(2)
    c_ethanol_ag     = Yasso_io_pools(3)
    c_nonsoluble_ag  = Yasso_io_pools(4)

    c_acid_bg        = Yasso_io_pools(5)
    c_water_bg       = Yasso_io_pools(6)
    c_ethanol_bg     = Yasso_io_pools(7)
    c_nonsoluble_bg  = Yasso_io_pools(8)

    c_humus          = Yasso_io_pools(9)

    ! get the sum of the AWEN pools before decomposition
    Pseudo_litter_green = c_acid_ag          &
                        + c_water_ag         &
                        + c_ethanol_ag       &
                        + c_nonsoluble_ag    &
                        + c_acid_bg          &
                        + c_water_bg         &
                        + c_ethanol_bg       &
                        + c_nonsoluble_bg

    ! Change units of climatic forcing variables
    precip = Weather(2)*sec_per_year/1000._wp  ! mm/s -> m/a
    temp   = Weather(1) - Tmelt                ! K -> C

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!! START CALCULATION OF SOIL CARBON !!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! 1. Calculate decomposition rates

    ! Temperature dependence of decomposition
    d_temp   = EXP(temp_p1*temp + temp_p2*temp**2.0_wp)
    ! Precipitation dependence of decomposition
    d_precip = 1.0_wp - EXP(precip_p1*precip)
    ! Litter size dependence of decomposition -- no effect if WoodlitterSize = 0.0
    d_size   = MIN(1.0_wp,(1.0_wp + size_p1 * WoodLitterSize + size_p2 * WoodLitterSize**2.0_wp)**size_p3)

     ! decomposition rates accounting for temperature, precipitation, litter size, and nutrient limitation
     d_acid       =  redFactor * ref_decomp_rate_acid       * d_temp * d_precip * d_size
     d_water      =  redFactor * ref_decomp_rate_water      * d_temp * d_precip * d_size
     d_ethanol    =  redFactor * ref_decomp_rate_ethanol    * d_temp * d_precip * d_size
     d_nonsoluble =  redFactor * ref_decomp_rate_nonsoluble * d_temp * d_precip * d_size

     ! the decomposition of humus is not limited by nutrients; thus no redFactor
     d_humus      =              ref_decomp_rate_humus      * d_temp * d_precip        ! no size effect on humus

    ! 2. Calculate fluxes

    ! loss fluxes (negative values):
    Cflx_from_acid_ag         = c_acid_ag       * d_acid
    Cflx_from_water_ag        = c_water_ag      * d_water
    Cflx_from_ethanol_ag      = c_ethanol_ag    * d_ethanol
    Cflx_from_nonsoluble_ag   = c_nonsoluble_ag * d_nonsoluble

    Cflx_from_acid_bg         = c_acid_bg       * d_acid
    Cflx_from_water_bg        = c_water_bg      * d_water
    Cflx_from_ethanol_bg      = c_ethanol_bg    * d_ethanol
    Cflx_from_nonsoluble_bg   = c_nonsoluble_bg * d_nonsoluble

    Cflx_from_humus           = c_humus         * d_humus

    ! gain fluxes (positive values):
    ! fixed fractions of each loss flux enters another pool; REMARK: the fraction can be zero (see above why)
    Cflx_2_acid_ag           = ABS(Cflx_from_water_ag      * W_2_A  &    ! returns positive fluxes
                                 + Cflx_from_ethanol_ag    * E_2_A  &
                                 + Cflx_from_nonsoluble_ag * N_2_A)

    Cflx_2_water_ag          = ABS(Cflx_from_acid_ag       * A_2_W  &
                                 + Cflx_from_ethanol_ag    * E_2_W  &
                                 + Cflx_from_nonsoluble_ag * N_2_W)

    Cflx_2_ethanol_ag        = ABS(Cflx_from_acid_ag       * A_2_E  &
                                 + Cflx_from_water_ag      * W_2_E  &
                                 + Cflx_from_nonsoluble_ag * N_2_E)

    Cflx_2_nonsoluble_ag     = ABS(Cflx_from_acid_ag       * A_2_N  &
                                 + Cflx_from_water_ag      * W_2_N  &
                                 + Cflx_from_ethanol_ag    * E_2_N)

    Cflx_2_acid_bg           = ABS(Cflx_from_water_bg      * W_2_A  &
                                 + Cflx_from_ethanol_bg    * E_2_A  &
                                 + Cflx_from_nonsoluble_bg * N_2_A)

    Cflx_2_water_bg          = ABS(Cflx_from_acid_bg       * A_2_W  &
                                 + Cflx_from_ethanol_bg    * E_2_W  &
                                 + Cflx_from_nonsoluble_bg * N_2_W)

    Cflx_2_ethanol_bg        = ABS(Cflx_from_acid_bg       * A_2_E  &
                                 + Cflx_from_water_bg      * W_2_E  &
                                 + Cflx_from_nonsoluble_bg * N_2_E)

    Cflx_2_nonsoluble_bg     = ABS(Cflx_from_acid_bg       * A_2_N  &
                                 + Cflx_from_water_bg      * W_2_N  &
                                 + Cflx_from_ethanol_bg    * E_2_N)

    Cflx_2_humusAG        = ABS(Cflx_from_acid_ag                &
                              + Cflx_from_water_ag               &
                              + Cflx_from_ethanol_ag             &
                              + Cflx_from_nonsoluble_ag          &
                              ) * AWEN_2_H

    Cflx_2_humusBG        = ABS(Cflx_from_acid_bg                &
                              + Cflx_from_water_bg               &
                              + Cflx_from_ethanol_bg             &
                              + Cflx_from_nonsoluble_bg          &
                              ) * AWEN_2_H

    Cflx_2_humus          = Cflx_2_humusAG + Cflx_2_humusBG

    ! the remaining fractions of the loss fluxes enter the atmosphere as respiration
    soilResp_rateYasso    =  (Cflx_from_acid_ag + Cflx_from_acid_bg)               &
                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                             + (Cflx_from_water_ag + Cflx_from_water_bg)           &
                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                             + (Cflx_from_ethanol_ag + Cflx_from_ethanol_bg)       &
                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                             + (Cflx_from_nonsoluble_ag + Cflx_from_nonsoluble_bg) &
                             * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)         &
                             + Cflx_from_humus

    ! litter ag & bg respiration (needed for N cycle)
    soilResp_rateLitterAG   =  (Cflx_from_acid_ag )                                &
                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                             + (Cflx_from_water_ag )                               &
                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                             + (Cflx_from_ethanol_ag )                             &
                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                             + (Cflx_from_nonsoluble_ag )                          &
                             * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)

    soilResp_rateLitterBG   =  (Cflx_from_acid_bg )                                &
                             * (1.0_wp - A_2_W - A_2_E - A_2_N - AWEN_2_H)         &
                             + (Cflx_from_water_bg )                               &
                             * (1.0_wp - W_2_A - W_2_E - W_2_N - AWEN_2_H)         &
                             + (Cflx_from_ethanol_bg )                             &
                             * (1.0_wp - E_2_A - E_2_W - E_2_N - AWEN_2_H)         &
                             + (Cflx_from_nonsoluble_bg )                          &
                             * (1.0_wp - N_2_A - N_2_W - N_2_E - AWEN_2_H)

    ! 3. update Yasso pools

    ! 3a. redistribution among yasso pools
    c_acid_ag           = MAX(0.0_wp,                                         &
                             c_acid_ag                                        & ! old pool
                             + (Cflx_from_acid_ag                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_ag) *dt)
    c_water_ag          = MAX(0.0_wp,                                         & ! and so on .....
                             c_water_ag                                       &
                             + (Cflx_from_water_ag                            &
                               + Cflx_2_water_ag) *dt)
    c_ethanol_ag        = MAX(0.0_wp,                                         &
                             c_ethanol_ag                                     &
                             + (Cflx_from_ethanol_ag                          &
                               + Cflx_2_ethanol_ag) *dt)
    c_nonsoluble_ag     = MAX(0.0_wp,                                         &
                             c_nonsoluble_ag                                  &
                             + (Cflx_from_nonsoluble_ag                       &
                               + Cflx_2_nonsoluble_ag) *dt)

    c_acid_bg           = MAX(0.0_wp,                                         &
                             c_acid_bg                                        & ! old pool
                             + (Cflx_from_acid_bg                             & ! incoming flux from AWEN pools
                               + Cflx_2_acid_bg) *dt)
    c_water_bg          = MAX(0.0_wp,                                         & ! ...
                             c_water_bg                                       &
                             + (Cflx_from_water_bg                            &
                               + Cflx_2_water_bg) *dt                         &
                             + NPP_2_rootExudates)                               ! exudates are carbohydrates only and belowground
    c_ethanol_bg        = MAX(0.0_wp,                                         &
                             c_ethanol_bg                                     &
                             + (Cflx_from_ethanol_bg                          &
                               + Cflx_2_ethanol_bg) *dt)
    c_nonsoluble_bg     = MAX(0.0_wp,                                         &
                             c_nonsoluble_bg                                  &
                             + (Cflx_from_nonsoluble_bg                       &
                               + Cflx_2_nonsoluble_bg) *dt)

    c_humus          = MAX(0.0_wp,                            &
                             c_humus                          &
                             + (Cflx_from_humus               &
                               + Cflx_2_humus) *dt)

    ! 3b. distribute litter
    ! ag
    CALL distribute_yasso_litter(c_acid_ag, c_water_ag, c_ethanol_ag, c_nonsoluble_ag, c_humus, &
      &                          litter, fract_aboveground, Lit_coefV)
    ! bg
    CALL distribute_yasso_litter(c_acid_bg, c_water_bg, c_ethanol_bg, c_nonsoluble_bg, c_humus, &
      &                          litter, 1.0_wp - fract_aboveground, Lit_coefV)

    ! Cflx_litter_2_humus is not only required for distribution, but possibly also as output for the eq script
    Cflx_litter_2_humus = Lit_coefV(i_lctlib_humus) * litter

    ! compute d_litter_green; this is given by the change in AWEN pools
    !          Flx                 Pool(t)                   Pool(t+1)
    !    d_P = ----     &   Flx = ----------   --->   d_P = -----------  - 1
    !          Pool(t)             Pool(t+1)                 Pool(t)

    IF (Pseudo_litter_green .GT. 0.0_wp) THEN
       d_litter_green = (c_acid_ag            &
                         + c_water_ag         &
                         + c_ethanol_ag       &
                         + c_nonsoluble_ag    &
                         + c_acid_bg          &
                         + c_water_bg         &
                         + c_ethanol_bg       &
                         + c_nonsoluble_bg)   &
                         / Pseudo_litter_green - 1._wp
    ELSE
      d_litter_green = 0.0_wp
    ENDIF

    ! 4.1 Write pools into the output array
    Yasso_out(1)           = c_acid_ag
    Yasso_out(2)           = c_water_ag
    Yasso_out(3)           = c_ethanol_ag
    Yasso_out(4)           = c_nonsoluble_ag
    Yasso_out(5)           = c_acid_bg
    Yasso_out(6)           = c_water_bg
    Yasso_out(7)           = c_ethanol_bg
    Yasso_out(8)           = c_nonsoluble_bg
    Yasso_out(9)           = c_humus
    ! 4.2 Write respiration into the output array
    Yasso_out(10)           = soilResp_rateYasso * dt ! convert back to daily

    IF (.NOT.PRESENT(redFact_Nlimit)) THEN
      ! the routine might have been called to diagnose the nutrient demand, so we give out
      ! 1. flux for immobilisation: this is the sum of all carbon coming from the non-humus pools to humus
      Yasso_out(11)           = Cflx_2_humusAG * dt ! daily ! the flux is positive
      Yasso_out(12)           = Cflx_2_humusBG * dt ! daily ! the flux is positive
      ! 2. flux for mineralisation:  this is the flux from humus to atmosphere (Cflx_slow_2_atmos)
      Yasso_out(13)           =  - Cflx_from_humus * dt ! daily ! sign CHECKED!
      ! ... as well as all fluxes from AWEN pools to atmosphere
      Yasso_out(14)           = -soilResp_rateLitterAG * dt !daily ! the flux is negative
      Yasso_out(15)           = -soilResp_rateLitterBG * dt !daily ! the flux is negative
      ! 3. decomposition rate for Npool litter green
      Yasso_out(16)           = d_litter_green
    ENDIF

    ! 4.3 Write diagnostics to determine equilibrium of humus pools (according to jsb3 svn rev. )
    Yasso_out(17)           = -1._wp * Cflx_from_humus * dt
    Yasso_out(18)           = (Cflx_2_humus * dt) + Cflx_litter_2_humus

  END SUBROUTINE yasso


  ! ====================================================================================================== !
  !
  !> Distributes incomming litter to yasso pools (i.e. proportional to litter coefficients and passed factors)
  !
#ifndef _OPENACC
  PURE &
#endif
  SUBROUTINE distribute_yasso_litter(c_acid, c_water, c_ethanol, c_nonsoluble, c_humus, &
    &                          litter, fraction, coefficients)

    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(INOUT)             :: c_acid        !< yasso acid pool
    REAL(wp), INTENT(INOUT)             :: c_water       !< yasso water pool
    REAL(wp), INTENT(INOUT)             :: c_ethanol     !< yasso ethanol pool
    REAL(wp), INTENT(INOUT)             :: c_nonsoluble  !< yasso nonsoluble pool
    REAL(wp), INTENT(INOUT)             :: c_humus       !< yasso humus pool
    REAL(wp), INTENT(IN)                :: litter        !< litter input  (above and belowground) [mol(C)/m2/day]
    REAL(wp), INTENT(IN)                :: fraction      !< fraction of litter put into yasso pools []
    REAL(wp), DIMENSION(5), INTENT(IN)  :: coefficients  !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !

    CALL add_litter_to_yasso_pool(c_acid,       litter, fraction, coefficients(i_lctlib_acid))
    CALL add_litter_to_yasso_pool(c_water,      litter, fraction, coefficients(i_lctlib_water))
    CALL add_litter_to_yasso_pool(c_ethanol,    litter, fraction, coefficients(i_lctlib_ethanol))
    CALL add_litter_to_yasso_pool(c_nonsoluble, litter, fraction, coefficients(i_lctlib_nonsoluble))
    CALL add_litter_to_yasso_pool(c_humus,      litter, fraction, coefficients(i_lctlib_humus))

  END SUBROUTINE distribute_yasso_litter

  ! ====================================================================================================== !
  !
  !> Adds incomming litter share (i.e. proportional to litter coefficient and factor)
  !
#ifndef _OPENACC
  ELEMENTAL &
#endif
  SUBROUTINE add_litter_to_yasso_pool(this_pool, litter, fraction, coefficient)

!$ACC ROUTINE SEQ

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    REAL(wp), INTENT(INOUT) :: this_pool     !< target yasso pool
    REAL(wp), INTENT(IN)    :: litter        !< litter input  (above and belowground) [mol(C)/m2/day]
    REAL(wp), INTENT(IN)    :: fraction      !< fraction of litter put into this pool (BG vs AG) [ ]
    REAL(wp), INTENT(IN)    :: coefficient   !< litter coefficient (acid vs water vs ...) [ ]
    ! -------------------------------------------------------------------------------------------------- !
    this_pool           = MAX(0.0_wp, this_pool + (fraction * coefficient * litter))

  END SUBROUTINE add_litter_to_yasso_pool

  ! --- get_per_tile() ---------------------------------------------------------------------------------------------------
  !
  ! calculate per tile density from density per canopy m2
  !
  SUBROUTINE get_per_tile(  &
    & Cpool_per_tile,                 & ! out
    & Cpool_per_canopy,               & ! in
    & veg_fract_correction,           & ! in
    & fract_fpc_max)                    ! in

    REAL(wp), INTENT(IN), DIMENSION(:) :: &
      & Cpool_per_canopy, &
      & fract_fpc_max,    &
      & veg_fract_correction
    REAL(wp), INTENT(OUT), DIMENSION(:) :: &
      & Cpool_per_tile

    INTEGER :: ic

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, SIZE(fract_fpc_max)
      Cpool_per_tile(ic) = Cpool_per_canopy(ic) * veg_fract_correction(ic) * fract_fpc_max(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE get_per_tile

#endif
END MODULE mo_carbon_process
