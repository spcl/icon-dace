!> QUINCY biophysical rate modifiers
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
!>#### routines for the task update_sb_rate_modifier, i.e., biophysical rate modifiers
!>
MODULE mo_q_sb_rate_modifier
#ifndef __NO_QUINCY__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, message_text, finish
  USE mo_jsb_control,     ONLY: debug_on

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_rate_modifier

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_sb_rate_modifier'

CONTAINS

  ! ======================================================================================================= !
  !>biophysical rate modifier
  !>
  SUBROUTINE update_sb_rate_modifier(tile, options)

    USE mo_jsb_class,              ONLY: Get_model
    USE mo_jsb_process_class,      ONLY: SB_, SPQ_
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,         ONLY: t_jsb_task_options
    USE mo_jsb_model_class,        ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,       ONLY: t_lctlib_element
    USE mo_jsb_physical_constants, ONLY: Tzero
    USE mo_sb_constants,           ONLY: ea_depolymerisation, ea_mic_uptake, ea_sorption, ea_desorption, ea_hsc, &
      &                                  kmr_depolymerisation, kmr_mic_uptake, kmr_hsc, kmr_sorption, &
      &                                  ea_nitrification, ed_nitrification, ea_denitrification, ea_gasdiffusion, &
      &                                  t_opt_nitrification, lambda_anvf, k_anvf, &
      &                                  ea_decomposition, ed_decomposition,t_ref_decomposition, kmr_decomposition, &
      &                                  kmr_asymb_bnf
    USE mo_veg_constants,          ONLY: ea_symb_bnf, ed_symb_bnf, t_opt_symb_bnf
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model              !< the model
    TYPE(t_lctlib_element), POINTER   :: lctlib             !< land-cover-type library - parameter across pft's
    INTEGER                           :: iblk, ics, ice, nc !< dimensions
    INTEGER                           :: is, ic             !< dimensions looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_rate_modifier'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! SB_
    dsl4jsb_Real3D_onChunk :: rtm_decomposition
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation
    dsl4jsb_Real3D_onChunk :: rtm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rtm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rtm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: rmm_decomposition
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk :: num_sl_above_bedrock
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    dsl4jsb_Real3D_onChunk :: w_soil_fc_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sat_sl
    dsl4jsb_Real3D_onChunk :: w_soil_pot_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_decomposition)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_mic_uptake)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_plant_uptake)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)                     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_asymb_bnf)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_decomposition)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_plant_uptake)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)                     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_asymb_bnf)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)! out
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_,     num_sl_above_bedrock)        ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_,     t_soil_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     w_soil_sl)                   ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     w_soil_fc_sl)                ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     w_soil_sat_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     w_soil_pot_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SPQ_,     soil_depth_sl)               ! in
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc
      DO is = 1,INT(num_sl_above_bedrock(ic))
        !>1.0 soil organic rate calc modifiers
        !>

        !>  1.1 simple model
        !>
        rtm_decomposition(ic, is)    = calc_peaked_temperature_rate_modifier(t_soil_sl(ic, is), &
          &                                                          ea_decomposition, &
          &                                                          ed_decomposition, &
          &                                                          t_ref_decomposition)
        rmm_decomposition(ic, is)    = MAX(1._wp - w_soil_pot_sl(ic, is) / kmr_decomposition, 0.0_wp)
        !>  1.2 JSM
        !>
        rtm_depolymerisation(ic, is) = calc_temperature_rate_modifier(t_soil_sl(ic, is),ea_depolymerisation)
        rtm_mic_uptake(ic, is)       = calc_temperature_rate_modifier(t_soil_sl(ic, is),ea_mic_uptake)
        rtm_sorption(ic, is)         = calc_temperature_rate_modifier(t_soil_sl(ic, is),ea_sorption)
        rtm_desorption(ic, is)       = calc_temperature_rate_modifier(t_soil_sl(ic, is),ea_desorption)
        rtm_hsc(ic, is)              = calc_temperature_rate_modifier(t_soil_sl(ic, is),ea_hsc)

        rmm_depolymerisation(ic, is) = calc_oxygen_rate_modifier(w_soil_sl(ic, is),w_soil_fc_sl(ic, is),kmr_depolymerisation)
        rmm_mic_uptake(ic, is)       = calc_oxygen_rate_modifier(w_soil_sl(ic, is),w_soil_fc_sl(ic, is),kmr_mic_uptake)
        rmm_sorption(ic, is)         = calc_convert_bulk_to_liquid(w_soil_sl(ic, is),soil_depth_sl(ic, is),kmr_sorption)
        rmm_desorption(ic, is)       = calc_convert_bulk_to_liquid(w_soil_sl(ic, is),soil_depth_sl(ic, is),kmr_sorption)
        rmm_hsc(ic, is)              = calc_substrate_rate_modifier(w_soil_sl(ic, is),w_soil_fc_sl(ic, is),kmr_hsc)

        !>2.0 nitrification-denitrification model
        !>
        rtm_nitrification(ic, is)             = calc_peaked_temperature_rate_modifier(t_soil_sl(ic, is), ea_nitrification, &
          &                                     ed_nitrification, t_opt_nitrification)
        rtm_denitrification(ic, is)           = calc_temperature_rate_modifier(t_soil_sl(ic, is), ea_denitrification)
        rtm_gasdiffusion(ic, is)              = calc_temperature_rate_modifier(t_soil_sl(ic, is), ea_gasdiffusion)
        rmm_nitrification(ic, is)             = calc_moisture_modifier(w_soil_sl(ic, is), w_soil_fc_sl(ic, is))
        rmm_gasdiffusion(ic, is)              = 1.0_wp - calc_moisture_modifier(w_soil_sl(ic, is), w_soil_fc_sl(ic, is))
        anaerobic_volume_fraction_sl(ic, is)  = calc_anaerobic_volume_fraction(w_soil_sl(ic, is), w_soil_sat_sl(ic, is), &
          &                                     lambda_anvf, k_anvf)

        !>  2.1 N fixation model
        !>   Using the symbiotic T-function as no other data available
        !>
        rtm_asymb_bnf(ic, is) = calc_peaked_temperature_rate_modifier(t_soil_sl(ic, is), ea_symb_bnf, &
          &                     ed_symb_bnf, t_opt_symb_bnf)
        rmm_asymb_bnf(ic, is) = MAX(1._wp - w_soil_pot_sl(ic, is) / kmr_asymb_bnf, 0.0_wp)

        !>3.0 root model
        !>
        rtm_plant_uptake(ic, is) = calc_peaked_temperature_rate_modifier(t_soil_sl(ic, is), ea_decomposition, &
          &                        ed_decomposition, t_ref_decomposition)
        rmm_plant_uptake(ic, is) = MAX(1._wp - w_soil_pot_sl(ic, is) / lctlib%phi_leaf_min, 0.0_wp)
      ENDDO
    ENDDO
  END SUBROUTINE update_sb_rate_modifier


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates temperature response factor for various soil biogeochemical reactions
  !!  according to the Arrhenius equation as defined in
  !!  Wang, G., Post, W.M., Mayes, M.A., Frerichs, J.T., Sindhu, J., 2012.
  !!  Parameter estimation for models of ligninolytic and cellulolytic enzyme kinetics.
  !!  Soil Biology and Biochemistry 48, 28-38.
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Soil temperatures (t_soil_sl) lower eps8 are caught by using eps8 instead of t_soil_sl for calculations,
  !!  that is, t_soil_sl is NOT modified in case: t_soil_sl < eps8
  !!
  !! Input: soil temperature, activiation energy
  !!
  !! Output: rate modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_temperature_rate_modifier(t_soil_sl,Ea) RESULT(rate_modifier)

    USE mo_jsb_physical_constants, ONLY: r_gas
    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_sb_constants,           ONLY: t_ref_tresponse

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: t_soil_sl     !< soil temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp           ! either t_soil_sl or eps8
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_temperature_rate_modifier'

    ! temperature must not be negative in the below function
    IF(t_soil_sl < eps8) THEN
      hlp = eps8
    ELSE
      hlp = t_soil_sl
    END IF

    IF(Ea > eps8)THEN
      rate_modifier = exp(-Ea/r_gas*(1._wp/hlp - 1._wp/t_ref_tresponse))
    ELSE
      rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_temperature_rate_modifier

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates temperature response factor for various soil biogeochemical reactions
  !!  according to the peaked Arrhenius equation as defined in
  !!  Medlyn et al. 2002, PCE, eq. 18
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil temperature, (de-)activiation energy, optimum temperature
  !!
  !! Output: rate modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_peaked_temperature_rate_modifier(t_soil_sl,Ea,Ed,t_opt) RESULT(rate_modifier)

    USE mo_jsb_physical_constants, ONLY: r_gas
    USE mo_jsb_math_constants,     ONLY: eps8
    USE mo_sb_constants,           ONLY: t_ref_tresponse

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: t_soil_sl     !< soil temperature [K]
    REAL(wp), INTENT(in)  :: Ea            !< activation energy [J mol-1]
    REAL(wp), INTENT(in)  :: Ed            !< deactivation energy [J mol-1]
    REAL(wp), INTENT(in)  :: t_opt         !< optimum temperature [K]
    REAL(wp)              :: rate_modifier !< temperature rate modifier [unitless]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp1, hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_temperature_rate_modifier'

    IF(Ea>eps8)THEN
       hlp1 = t_soil_sl - t_opt
       hlp2 = t_soil_sl * t_opt * r_gas
       rate_modifier = (Ed * exp(Ea * hlp1 / hlp2)) / (Ed - Ea * (1._wp - exp(Ed * hlp1 / hlp2)))
    ELSE
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_peaked_temperature_rate_modifier


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM
  !
  !-----------------------------------------------------------------------------------------------------
  !> Oxygen limitation of decomposition and uptake due to soil moisture response function
  !!  according to Davidson et al's DAMM model (2012)
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter
  !!
  !! Output: rate modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_oxygen_rate_modifier(w_soil_sl,w_soil_fc_sl,param) RESULT(rate_modifier)

    USE mo_sb_constants,        ONLY: k1_afps
    USE mo_jsb_math_constants,  ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl     !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl  !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param         !<
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    REAL(wp)              :: hlp
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_oxygen_rate_modifier'

    IF(param>eps8) THEN
       IF(w_soil_fc_sl > eps8) THEN
          air_filled_pore_space = (w_soil_fc_sl - w_soil_sl) / w_soil_fc_sl
          IF(air_filled_pore_space < 0.005_wp) air_filled_pore_space = 0.005_wp
          rate_modifier = (air_filled_pore_space**k1_afps) / (param + (air_filled_pore_space**k1_afps))
       ELSE
          rate_modifier = 1._wp
       ENDIF
    ELSE
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_oxygen_rate_modifier


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM
  !
  !-----------------------------------------------------------------------------------------------------
  !> Converts bulk density to liquid concentration
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter
  !!
  !! Output: conversion modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_convert_bulk_to_liquid(w_soil_sl, soil_depth_sl, param) RESULT(factor)

    USE mo_jsb_math_constants,  ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    !0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: soil_depth_sl       !< soil depth [m]
    REAL(wp), INTENT(in)  :: param               !<
    REAL(wp)              :: factor              !< conversion factor [unitless]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_convert_bulk_to_liquid'
    REAL(wp)              :: hlp1

    IF(param>eps8.AND.w_soil_sl > eps8) THEN
      hlp1 = (w_soil_sl / soil_depth_sl)
      IF (hlp1 < param) THEN
        factor = param
      ELSE
        factor = hlp1
      ENDIF
    ELSE
       factor = 1._wp
    ENDIF

  END FUNCTION calc_convert_bulk_to_liquid


  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module, derived from JSM
  !
  !-----------------------------------------------------------------------------------------------------
  !> Substrate limitation of decomposition and uptake due to soil moisture response function
  !!  according to Davidson et al's DAMM model (2012)
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameter
  !!
  !! Output: rate modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  FUNCTION calc_substrate_rate_modifier(w_soil_sl,w_soil_fc_sl,param) RESULT(rate_modifier)

    USE mo_sb_constants,        ONLY: k2_afps
    USE mo_jsb_math_constants,  ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl        !< soil moisture at field capacity [m]
    REAL(wp), INTENT(in)  :: param               !<
    REAL(wp)              :: rate_modifier !< oxygen availability rate modifier [unitless]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_substrate_rate_modifier'

    IF(param>eps8.AND.w_soil_fc_sl > eps8) THEN
       hlp1 = (w_soil_sl / w_soil_fc_sl) ** k2_afps
       IF(hlp1 > param) THEN
          rate_modifier = 1._wp / hlp1
       ELSE
          rate_modifier = 1._wp / param
       ENDIF
    ELSE
      rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_substrate_rate_modifier

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module
  !
  !-----------------------------------------------------------------------------------------------------
  !> Anaerobic fraction of the soil layer
  !!  according to Zaehle et al. 2011
  !!  NB if sensitivity of zero is specified, 1 will be returned for all layers
  !!
  !! Input: soil moisture, sensitivity parameters
  !!
  !! Output: rate modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_anaerobic_volume_fraction(w_soil_sl,w_soil_sat_sl,lambda,k) RESULT(rate_modifier)

    USE mo_jsb_math_constants,    ONLY: eps4

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_sat_sl       !< soil moisture at saturation [m]
    REAL(wp), INTENT(in)  :: lambda              !< lambda in Weibull function
    REAL(wp), INTENT(in)  :: k                   !< k in Weibull function
    REAL(wp)              :: rate_modifier       !< oxygen availability rate modifier [unitless]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_anaerobic_volume_fraction'

    air_filled_pore_space = (w_soil_sat_sl - w_soil_sl ) / w_soil_sat_sl
    IF(air_filled_pore_space > eps4) THEN
       rate_modifier = 1.0_wp - exp(-(lambda * (1.0_wp - air_filled_pore_space))**k)
    ELSE
       rate_modifier = 1._wp - exp(-(lambda * (1.0_wp - eps4))**k)
    ENDIF

  END FUNCTION calc_anaerobic_volume_fraction

  !-----------------------------------------------------------------------------------------------------
  ! Functionality used by the SB module
  !
  !-----------------------------------------------------------------------------------------------------
  !> soil moisture response function for the nitrification model ((FC-soil_moisture)/FC)
  !!
  !! Input: soil moisture, sensitivity parameter
  !!
  !! Output: rate modifier (unitless)
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_moisture_modifier(w_soil_sl,w_soil_fc_sl) RESULT(rate_modifier)

    USE mo_jsb_math_constants,  ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)  :: w_soil_sl           !< soil moisture [m]
    REAL(wp), INTENT(in)  :: w_soil_fc_sl        !< soil moisture at field capacity [m]
    REAL(wp)              :: rate_modifier       !< rate modifier [unitless]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)              :: air_filled_pore_space
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_moisture_modifier'

    IF(w_soil_fc_sl > eps8) THEN
       air_filled_pore_space = (w_soil_fc_sl - w_soil_sl) / w_soil_fc_sl
       rate_modifier = 1.0_wp - air_filled_pore_space
    ELSE
       rate_modifier = 1._wp
    ENDIF

  END FUNCTION calc_moisture_modifier

#endif
END MODULE mo_q_sb_rate_modifier
