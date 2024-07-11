!> Contains the routines for the phenology processes
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
MODULE mo_pheno_process
#ifndef __NO_JSBACH__

   USE mo_kind, ONLY: wp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: calc_summergreen_phenology , calc_evergreen_phenology, &
             calc_crop_phenology, calc_raingreen_phenology, calc_grass_phenology, & !, get_letItDie, get_letItGrow, get_shedRate_RG
             get_foliage_projected_cover, derive_maxLai_according_to_allometric_relationships, track_max_green_pool

   CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_process'

 CONTAINS

  ! --- update_growth_phase() -----------------------------------------------------------------------------------------------------
  !
  ! This subroutines contain the two models for determining dates of budburst ("spring event") and leaf fall ("autumn event") for
  ! the summergreen and evergreen phenologies. For operating the dynamical phenology it needs only be known in what "phase" it
  ! currently is, i.e. whether the current time step falls into the time from spring to autumn event ("growth"), or into the time
  ! from autumn to spring event ("rest") --- the budburst and leaf fall dates itself are not needed. Hence, this routine does not
  ! provide these dates. Instead, it updates the field "growth_phase": A positive entry means that at the particular grid
  ! cell the evergreens and summergreens are growing, whereas a negative entry means that they are at rest.
  !
  ! For the spring event the "alternating model" of Murray et al. [1] is used:
  ! Let "S(d)" denote the value of the heat sum at day "d":
  ! (1) S(d) = SUM(d'=d0,d) MAX(T(d')-T_alt,0),
  ! where "T(d)" is the mean day temperature at day "d", "T_alt" is the "alternating temperature" (which has the function of a
  ! cutoff temperature in the heat sum), and "d0" is the starting date for temperture summation. This starting date is determined
  ! by a temperature criterion (see below) --- once  more we need not know the date, but only whether summation has started:
  ! this information is kept also in the field "heat_sum", which is set to a value smaller than -1 during times where no
  ! heat summation takes place.
  ! Another key quantity of the alternating model the number of chill days "C(d)": this is the number of days with a mean day
  ! temperatu. below the alternating temperature "T_alt", where counting is started here at the day "d_a" of the last autumn event:
  ! (2) C(d) = SUM(d'=d_a,d) STEP(T(d)-T_alt);
  ! here STEP() is the Heaviside step function. From C(d) a critical heatsum "S_crit(d)" is computed:
  ! (3) S_crit(d) = S_crit_min + S_crit_range * exp(-C(d)/C_decay),
  ! where "S_crit_min", "S_crit_range" and "C_decay" are parameters: "S_crit_min" and "S_crit_range" define minimum value and
  ! maximum range of the critical heatsum, whereas "C_decay" determines how fast "S_crit(d)" decreases with increasing number of
  ! chill days. Finally, the spring event happens when first
  ! (4)  S(d) >= S_crit(d).
  ! Technically, at this date for the considered grid point the associated entry of the array growth_phase() is set to +1 and
  ! "springEvent_flag" is set to "true" NO!!. Moreover the chill days count is reset to zero and heat_sum set to -99 to indicate
  ! that summation has stopped.
  !
  ! The "autumn event" is calculated from a pseudo soil temperature: It happens when during the growth phase first the pseudo soil
  ! temperature Ts(d) falls below the critical soil temperature "Ts_crit" (variable "autumn_event_temp"); to prevent that this
  ! event is detected in spring the condition is added that the mean day air temperature T(d) is smaller than the soil temperature,
  ! i.e. the autum event happens, when first
  ! (5) T(d) < Ts(d) < Ts_crit
  ! At this event for the considered grid point the associated entries of the array growth_phase() are set to "-1".
  !  --- An additional daylength-criterion has not implemented yet!
  !
  ! It remains to determine the date "d0" for the start of heat summation: The idea is that heat summation starts when
  ! first the pseudo soil temperature during the rest phase gets larger than a critical temperature "heat_summation_start_temp"
  ! (a tuning parameter). Technically this means to initialize "heat_sum" at the particular grid point by zero; since this
  ! value is larger than -1 it indicates that heat summation has started.
  !
  ! This routine requires that the soil temperature is updated before calling it.
  !
  !
  ! [1] M.B. murray, M.G.R. Cannell and R.I. Smith, "Date of budburst of fifteen tree species in Britain following climatic
  !     warming", J. Appl. Ecology 26 (1989) 693-700). See also: A. Botta, N. Viovy, P. Ciais, P. Friedlingstein and P. Monfray,
  !     "A global prognostic scheme of leaf onset using satellite data", Global Change Biol. 6 (2000) 709-725.
  !




   SUBROUTINE calc_summergreen_phenology(day_of_year,            & ! <- in
                                         is_newday,              & ! <- in
                                         timeStep_in_days,       & ! <- in
                                         lat,                    & ! <- in
                                         previous_day_temp_mean, & ! <- in
                                         pseudo_soil_temp,       & ! <- in
                                         LaiMax_corrected,       & ! <- in
                                         days_since_growth_begin_SG, & ! <- inout
                                         chill_days_SG,          & ! <- inout
                                         heat_sum_SG,            & ! <- inout
                                         growth_phase_SG,        & ! <- inout
                                         lai )                     ! <- inout

    !$ACC ROUTINE SEQ

     !-----------------------------------------------------------------------
     !  DECLARATIONS
     USE mo_pheno_parameters,  ONLY: pheno_param_jsbach


     !-----------------------------------------------------------------------
     !  ARGUMENTS
     INTEGER,  intent(in)    :: day_of_year

     LOGICAL,  intent(in)    :: is_newday

     REAL(wp), intent(in)    :: timeStep_in_days

     REAL(wp), intent(in)    :: lat,                &
       &                        previous_day_temp_mean, &
       &                        LaiMax_corrected,   &
       &                        pseudo_soil_temp

     REAL(wp), intent(inout) :: days_since_growth_begin_SG, &
       &                        chill_days_SG,      &
       &                        heat_sum_SG,        &
       &                        growth_phase_SG,    &
       &                        lai                   ! Total LAI.


     !-----------------------------------------------------------------------
     !  LOCAL VARIABLES
     REAL(wp)          :: xtmp,   &  ! argument to calls of get_letItDie
                                     ! (necessary to allow for inlining of these functions)
       &                  grRate     ! argument to calls of get_letItGrow

     !-----------------------------------------------------------------------
     ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
     ! IF (one_of('calc_phenology', tile%components) < 1) RETURN

     !-----------------------------------------------------------------------
     ! CONTENT

     ! Determine if heat summation for summergreen has to start with zero. Note, this is exactly the same as for evergreens.
     IF  (growth_phase_SG < -0.5_wp                   & ! we are during REST PHASE ...
        .AND. heat_sum_SG < -1._wp                    & ! and heat summation has not yet started...
        .AND.                                         & ! and ...
         (                                            & ! ..
            (lat >= 0._wp .AND. day_of_year == 1)     & ! we are at northern hemisphere at January 1 ...
            .OR.                                      & !   or ...
            (lat <  0._wp .AND. day_of_year == 183)   & ! we are at southern hemisphere at July 2 (July 1 in leap years) ...
         )                                            &
         ) THEN                                         ! then  =>
           heat_sum_SG = 0._wp                          !          heat summation has to start, so initialize it by zero
     END IF

     ! To do only at first time step of a day:
     ! Update output variables: days_since_growth_begin_SG, chill_days_EG, heat_sum_SG, growth_phase_SG
     IF (is_newday) THEN

        ! Eventually update counter for days since last begin of growth phase. Note, this is the same as for evergreens.
        IF(days_since_growth_begin_SG > 0.0_wp .AND. days_since_growth_begin_SG < 365.0_wp) THEN
           days_since_growth_begin_SG = days_since_growth_begin_SG + 1.0_wp
        END IF

        ! If heat summation for summergreen has started count chill_days_EG or increase heat_sum_SG.
        ! Note, this is the same as for evergreens.
        IF(heat_sum_SG > -1._wp) THEN                                            ! heat summation has started ..
           IF(pseudo_soil_temp < pheno_param_jsbach%SG%alternation_temp) THEN ! the day is a chill day ..
              chill_days_SG = chill_days_SG + 1._wp                     ! then => increment chill days counter
                                                                        !         and limit the number of chill days to prevent
                                                                        !         values beyond any limit in polar regions
              chill_days_SG = MIN(chill_days_SG, pheno_param_jsbach%EG_SG%max_chill_days)
           ELSE                                                       ! else
                                                                        ! then => increase heat sum:
                                                                        !         add mean day temperature excess above
                                                                        !         alternating temperature
              heat_sum_SG = heat_sum_SG + pseudo_soil_temp - pheno_param_jsbach%SG%alternation_temp
           END IF
        END IF


        ! Update growth_phase_SG:
        ! Check begin and end of growth and vegetative phase
        ! For summergreen growth_phase_SG values are:
        !    1.0  .. in vegetative phase (i.e. from end of growth phase to autum event)
        !    0.0  .. during growth phase (i.e. some weeks after spring event)
        !   -1.0  .. during rest phase (i.e. from autumn event to spring event)
        IF (growth_phase_SG < -0.5_wp ) THEN                                  ! IF we are in the REST PHASE
           IF (heat_sum_SG > -1._wp                                            & ! if heat summation has already started ...
                   .AND.                                                       & ! and ...
               heat_sum_SG - pheno_param_jsbach%SG%heat_sum_min >                      & ! criticality condition is fulfilled ..
               pheno_param_jsbach%SG%heat_sum_range * EXP(-chill_days_SG / pheno_param_jsbach%SG%chill_decay_const)) THEN
                 growth_phase_SG = 0._wp                                          ! --> growth phase starts
                 days_since_growth_begin_SG = 1.0_wp                              ! start counting days since begin of growth phase
                 heat_sum_SG = -99._wp                                            ! --> set heat_sum to a value smaller than -1 to
                                                                                  ! indicate that heat summation is no more needed
                 chill_days_SG = 0._wp                                            ! --> reset chill days count
           END IF

        ELSE IF (growth_Phase_SG > 0.5_wp) THEN                               ! IF we are in the VEGETATIVE PHASE
           IF (previous_day_temp_mean < pseudo_soil_temp                       & ! Check end of vegetative phase
                                                                                 ! mean day temperature is smaller than soil temp.
                      .AND.                                                    & ! and ..
              pseudo_soil_temp < pheno_param_jsbach%SG%autumn_event_temp               & ! soil temperature falls below critical ..
                      .AND.                                                    & ! and ..
              days_since_growth_begin_SG > pheno_param_jsbach%SG%minLength_gvPhase)    & ! growth + vegetative phase have lasted the
                                           THEN                                  ! minimum time ..
                    growth_phase_SG = -1._wp                                     ! --> rest phase starts
           ELSE
              IF (days_since_growth_begin_SG > pheno_param_jsbach%SG%maxLength_gvPhase) & ! if growth + vegetative phase has reached
                                                                                  ! maximum length
                                       THEN                                       ! then =>
                    growth_phase_SG = -1._wp                                      ! --> rest phase starts
              END IF
           END IF

        ELSE IF (growth_Phase_SG > -0.5_wp .AND. growth_Phase_SG < 0.5_wp ) THEN   ! IF we are in the GROWTH PHASE
           IF(days_since_growth_begin_SG > pheno_param_jsbach%SG%growthPhaseLength) THEN   ! if growth phase of SGs has ended
                    growth_phase_SG = 1.0_wp                                       ! --> vegetative phase starts
           END IF
        END IF
     END IF ! from is_newday


     !!!!!!!! FINALLY, what it is all for:
     ! Update the output variable lai
     IF (growth_Phase_SG < -0.5_wp) THEN                                      ! if we are in the REST PHASE
        xtmp = pheno_param_jsbach%SG%shedRate_rest
        lai =  get_letItDie(timeStep_in_days,lai,xtmp)
     ELSE IF(growth_Phase_SG > 0.5_wp) THEN                                   ! if we are in the VEGETATIVE PHASE
        xtmp = pheno_param_jsbach%SG%shedRate_veget
        lai  = get_letItDie(timeStep_in_days,lai,xtmp)
     ELSE IF (growth_Phase_SG > -0.5_wp .AND. growth_Phase_SG < 0.5_wp ) THEN ! if we are in the GROWTH PHASE
        IF (lai < pheno_param_jsbach%all%laiSeed) THEN  ! if lai is smaller as pheno_param_jsbach%all%laiSeed then
                                                ! at least set it to pheno_param_jsbach%all%laiSeed but stay below laiMax
           lai = MIN(LaiMax_corrected-EPSILON(1.0_wp),pheno_param_jsbach%all%laiSeed)
        ENDIF
        IF(lai > pheno_param_jsbach%all%LAI_negligible) THEN ! if lai is not negligible then let it grow (otherwise there is nothing to do)
           grRate = pheno_param_jsbach%SG%growthRate
           lai = get_letItGrow(timeStep_in_days,lai, grRate, 0.0_wp, LaiMax_corrected)
        END IF
     END IF

   END SUBROUTINE calc_summergreen_phenology


  SUBROUTINE calc_evergreen_phenology(day_of_year,            & ! <- in
                                      is_newday,              & ! <- in
                                      timeStep_in_days,       & ! <- in
                                      lat,                    & ! <- in
                                      pseudo_soil_temp,       & ! <- in
                                      LaiMax_corrected,       & ! <- in
                                      days_since_growth_begin_EG, & ! <- inout
                                      chill_days_EG,          & ! <- inout
                                      heat_sum_EG,            & ! <- inout
                                      growth_phase_EG,        & ! <- inout
                                      lai )                     ! <- inout

    !$ACC ROUTINE SEQ

     !-----------------------------------------------------------------------
     !  DECLARATIONS
     USE mo_pheno_parameters,  ONLY: pheno_param_jsbach

     !-----------------------------------------------------------------------
     !  ARGUMENTS
     INTEGER,  intent(in)    :: day_of_year

     LOGICAL,  intent(in)    :: is_newday
     REAL(wp), intent(in)    :: timeStep_in_days
     REAL(wp), intent(in)    ::                     &
                                lat,                &
                                LaiMax_corrected,   &
                                pseudo_soil_temp

     REAL(wp), intent(inout) ::                     &
                                days_since_growth_begin_EG, &
                                chill_days_EG,      &
                                heat_sum_EG,        &
                                growth_phase_EG,    &
                                lai                   ! Total LAI.

     !-----------------------------------------------------------------------
     !  LOCAL VARIABLES
     REAL(wp)          :: xtmp,   &  ! argument to calls of get_letItDie
                                     ! (necessary to allow for inlining of these functions)
                          grRate     ! argument to calls of get_letItGrow

     !-----------------------------------------------------------------------
     ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
     ! IF (one_of('calc_phenology', tile%components) < 1) RETURN

     !-----------------------------------------------------------------------
     ! CONTENT

     ! Determine if HEAT SUMMATION FOR EVERGREEN HAS TO START from zero.  Note, this is nearly the same as for summergreens.
     IF ( growth_phase_EG > 0.5_wp                   & ! we are during vegetative phase (Note: there is no rest phase for EG)
         .AND.                                       & ! and ..
         heat_sum_EG < -1._wp                        & ! heat summation has not yet started
         .AND.                                       & ! and ..
         (                                           & ! ..
           (lat >= 0._wp .AND. day_of_year == 1)     & ! we are at northern hemisphere at January 1 ..
           .OR.                                      & ! or ..
           (lat <  0._wp .AND. day_of_year == 183)   & ! we are at southern hemisphere at July 2 (July 1 in leap years) ..
         )                                           & ! ..
        ) THEN                                          ! --> heat summation has to start,
         heat_sum_EG = 0._wp                            !     so initialize it by zero
     END IF

     ! To do only at first time step of a day:
     ! Update output variables: days_since_growth_begin_EG, chill_days_EG, heat_sum_EG, growth_phase_EG
     IF (is_newday) THEN
         ! Eventually update counter for days since last begin of growth phase (stop counting after 365 days).
         ! Note, this is the same as for summergreens.
         IF(days_since_growth_begin_EG > 0.0_wp .AND. days_since_growth_begin_EG < 365.0_wp) THEN
            days_since_growth_begin_EG = days_since_growth_begin_EG + 1.0_wp
         END IF

         ! If heat summation for evergreen has started count chill_days_EG or increase heat_sum_EG.
         ! Note, this is the same as for summergreens.
         IF(heat_sum_EG > -1.0_wp ) THEN                                           ! heat summation has started
            IF(pseudo_soil_temp < pheno_param_jsbach%EG%alternation_temp ) THEN    ! the day is a chill day ..
               chill_days_EG = chill_days_EG + 1._wp                                                  ! --> increment chill days counter
               chill_days_EG = &                                                                      !     limit the number of chill days to prevent
                    MIN(chill_days_EG, pheno_param_jsbach%EG_SG%max_chill_days)                       !     values beyond any limit in polar regions
            ELSE                                                                                      ! increase heat sum:
               heat_sum_EG = heat_sum_EG +  pseudo_soil_temp - pheno_param_jsbach%EG%alternation_temp ! --> add mean day temperature excess
                                                                                                      !     above alternating temperature
            END IF
         END IF

         ! Check begin and end of growth and vegetative phase
         ! For evergreen growth_phase_EG values are:
         !    1.0  .. in vegetative phase (i.e. from ent of growth phase to spring event)
         !    0.0  .. during growth phase (i.e. some weeks after spring event)
         ! DOES NOT EXIST FOR EVERGREEN:   -1.0  .. during rest phase (i.e. from autumn event to spring event)
         !  -99 .. Heat summation has not yet started (at all or because growth phase started)

         IF(growth_phase_EG < 0.5_wp ) THEN                                               ! evergreens are during GROWTH
            IF(days_since_growth_begin_EG > pheno_param_jsbach%EG%growthPhaseLength) THEN ! growth phase of EGs has ended
               growth_phase_EG = 1.0_wp                                                   ! --> vegetative phase starts
            END IF
         ELSE                                                                     ! everergreen are in VEGETATIVE PHASE
            IF(heat_sum_EG > -1._wp                                             & ! if heat summation has started ..
                 .AND.                                                          & ! and ..
               heat_sum_EG - pheno_param_jsbach%EG%heat_sum_min >               & ! criticality condition is fulfilled ..
               pheno_param_jsbach%EG%heat_sum_range * EXP(-chill_days_EG /      & ! ..
               pheno_param_jsbach%EG%chill_decay_const)) THEN                     ! ..
               growth_phase_EG = 0._wp                                            ! --> growth phase starts
               days_since_growth_begin_EG = 1._wp                                 ! --> reinitialize counter for days since
                                                                                  !     growth phase begin
               heat_sum_EG = -99._wp                                              ! --> set heat_sum to a value smaller than
                                                                                  !     -1 to indicate that heat summation is
                                                                                  !     no more needed
               chill_days_EG = 0._wp                                              ! --> reset chill days count
            END IF
         END IF
     END IF ! from is_newday


     !!!!!!!! FINALLY, what it is all for:
     ! Update the output variable lai
     IF(ABS(growth_phase_EG) < EPSILON(1.0_wp)) THEN    ! we are during the GROWTH PHASE
        IF (lai < pheno_param_jsbach%all%laiSeed) THEN
           ! start growth at least with seed value smaller than laiMax
           lai = MIN(LaiMax_corrected-EPSILON(1.0_wp),pheno_param_jsbach%all%laiSeed)
        END IF
        IF(lai > pheno_param_jsbach%all%LAI_negligible) THEN
           ! otherwise nothing to do
           grRate = pheno_param_jsbach%EG%growthRate
           lai = get_letItGrow(timeStep_in_days,lai, grRate, 0.0_wp, LaiMax_corrected)
        END IF
     ELSE                                   ! we are during the VEGETATIVE PHASE
        xtmp = pheno_param_jsbach%EG%shedRate_vegetative
        lai  = get_letItDie(timeStep_in_days,lai,xtmp)
     END IF


   END SUBROUTINE calc_evergreen_phenology


   SUBROUTINE calc_crop_phenology(day_of_year,                  & ! <- in
                                  is_newday,                    & ! <- in
                                  timeStep_in_days,             & ! <- in
                                  NPP_pot_rate_ca,              & ! <- in
                                  lat,                          & ! <- in
                                  pseudo_soil_temp,             & ! <- in
                                  LaiMax_corrected,             & ! <- in
                                  previous_day_NPP_pot_rate_ca, & ! <- in
                                  w_soil_filling,               & ! <- in
                                  previous_day_temp_min,        & ! <- in
                                  previous_day_temp_max,        & ! <- in
                                  lctlib_specificLeafArea_C,    & ! <- in
                                  heat_sum_CRP,                 & ! <- inout
                                  heat_sum_winter,              & ! <- inout
                                  growth_phase_CRP,             & ! <- inout
                                  lai )                           ! <- inout

    !$ACC ROUTINE SEQ

     !-----------------------------------------------------------------------
     !  DECLARATIONS
     USE mo_pheno_parameters,  ONLY: pheno_param_jsbach

     !-----------------------------------------------------------------------
     !  ARGUMENTS
     INTEGER,  intent(in)    :: day_of_year

     LOGICAL,  intent(in)    :: is_newday
     REAL(wp), intent(in)    :: timeStep_in_days ! express the time step as a fraction of a day. = dtime / 86400

     REAL(wp), intent(in)    :: lctlib_specificLeafArea_C

     REAL(wp), intent(in)    ::                     &
                                NPP_pot_rate_ca,    &
                                lat,                &
                                LaiMax_corrected,   &
                                previous_day_NPP_pot_rate_ca, &
                                pseudo_soil_temp,      &
                                w_soil_filling,        & ! Fractional filling of the soil water buckets.
                                                         ! R: in JSBACH3 wurde fuer Crops immer der soil layer 1,
                                                         !    d.h. der oberste Layer verwendet!
                                previous_day_temp_min, &
                                previous_day_temp_max

     REAL(wp), intent(inout) ::                   &
                                heat_sum_CRP,     &
                                heat_sum_winter,  &
                                growth_phase_CRP, &
                                lai                 ! Total LAI.
     !-----------------------------------------------------------------------
     !  LOCAL VARIABLES
     REAL(wp)          :: xtmp,         &  ! argument to calls of get_letItDie
                                           ! (necessary to allow for inlining of these functions)
                          grRate,       &  ! argument to calls of get_letItGrow
                          temp_min_max

     LOGICAL           :: NPP_positive     ! flag to indicate the existence of a grass or crop type with a
                                           ! positive NPP

     !-----------------------------------------------------------------------
     ! CONTENT

     ! Process?
     ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
     ! IF (one_of('calc_phenology', tile%components) < 1) RETURN

     ! NPP mask needed for crop growth phase and harvest
     ! R: In JSBACH3 we check grass or crop:
     !    IF (pheno_type_of_tile(...) == 4 .OR. pheno_type_of_tile(...) == 5) THEN ! R: grass and crop
     !       IF (previous_day_NPP_pot_rate_ca(...) > 0.0_wp) positive_NPP(...) = .true.
     !
     !    NPP_positive is only used to decide if the mean day temperature shall be added to the
     !    heat sum for crops (heat_sum_CRP and heat_sum_winter), which are only used to decide
     !    if the crops are in the rest phase (lai gets low and C of leafes is put into the litter C
     !    pool, which represents "harvest") or if they are still in the growth phase.
     !    In JSBACH4 we check only crop because within the JSBACH4 free tile structure it would be
     !    complicated to check if any grass tile somewhere has positive NPP.
     !    However, this parametrization could be made in a new procedure that represents the conditions
     !    needed by crops to honor temperatures...


     IF (previous_day_NPP_pot_rate_ca > - 1.E-10_wp) THEN
       NPP_positive = .TRUE.  ! Do not set to zero, as crops could never grow then...
     ELSE
       NPP_positive = .FALSE.
     END IF

     ! To do only at first time step of a day:
     ! Update output variables: days_since_growth_begin_EG, chill_days_EG, heat_sum_EG, growth_phase_EG
     IF (is_newday) THEN

        ! Crops growth_phase_CRP includes 2 Dimensions: the kind of cropping (fS, sS, W) and the phase (growth, rest).
        ! For extra-tropical crops growth_phase_CRP values are:
        !      2.0  .. during growth phase of winter crops (i.e. from autumn to spring)
        !      1.0  .. during growth phase of second summer crops (i.e. from summer to autumn)
        !      0.0  .. during growth phase of first  summer crops (i.e. from spring to autumn)
        !     -1.0  .. during rest   phase of all    summer crops (i.e. from autumn to spring)
        !     -2.0  .. during rest   phase of winter crops (i.e. from spring to autumn)
        ! Note, these phases exclude each other because there are either winter or summer crops for a grid cell (but never both).
        !
        !             |SPRING      |SUMMER          |AUTUMN          WINTER       |SPRING     |SUMMER    ...
        !------------------------------------------------------------------------------------------
        ! Summer      |-1  0       |0 or 1          |                             |           |
        ! Crops       |     x-can switch to -1 if heatsum is reached-x            |           |
        !------------------------------------------------------------------------------------------
        ! Winter      |2/-2        |-2              |-1 (goto SC) or -2           |           |
        ! Crops       |            |                |x-can switch to 2 if heatsum is reached-x|


        ! At the beginning of spring:
        IF (lat >= 0._wp .AND. day_of_year == 70      & ! we are at northern hemisphere at March 11 (March 10 in leap years)
           .OR.                                       & ! or ..
           lat <  0._wp .AND. day_of_year == 252      & ! we are at southern hemisphere at September 9 (September 8 in leap years) ..
           ) THEN
             heat_sum_CRP = 0._wp                         ! => heat summation has to start, so initialize it by zero
             IF (ABS(growth_phase_CRP) < 1.5_wp) THEN     ! if there are no winter crops at all (neither in rest nor in growth) ..
                growth_phase_CRP = -1._wp                 !     => crop phase is set to rest phase of summer crops
             END IF
        END IF

        ! At the beginning of summer:
        IF (lat >= 0._wp .AND. day_of_year == 172     & ! we are at northern hemisphere at June 21 (June 20 in leap years)
           .OR.                                       & ! or ..
           lat <  0._wp .AND. day_of_year == 354      & ! we are at southern hemisphere at December 21 (December 20 in leap years)
           ) THEN

             IF (ABS(growth_phase_CRP) < 1.5_wp .AND.           &         ! if there are no winter crops at all (neither in rest nor in
                 heat_sum_CRP > pheno_param_jsbach%CRP%heat_sum_harvest & ! growth) and the heat sum is enough for harvesting
                 ) THEN
                 heat_sum_CRP = 0._wp                               ! then heat summation has to start again, so initialize it by
                 growth_phase_CRP = 1._wp                           ! zero and crop phase is set to growth phase of double cropping
                                                                    ! (crops already have been harvested try double cropping.)
             ELSE IF (growth_phase_CRP > 1.5_wp) THEN             ! if there are winter crops still in the growing phase ..
                growth_phase_CRP = -2._wp                         ! => then switch to the rest phase of winter crops
             END IF
             IF (heat_sum_winter > 0._wp)  THEN              ! Heat sum for winter is multiplied by -1 to indicate ..
                heat_sum_winter = -1._wp * heat_sum_winter   ! => that there is no heat summation during summer
             END IF
        END IF

        ! In autumn:
        IF (lat >= 0._wp .AND. day_of_year == 289    & ! we are at northern hemisphere at October 16 (October 15 in leap
           .OR.                                      & ! years) or ..
           lat <  0._wp .AND. day_of_year == 105     & ! we are at southern hemisphere at April 15 (April 14 in leap years)
           ) THEN

            ! If there are winter crops and they were not harvested last winter or spring and summer crops would have been
            ! harvested already  --> crop phase is SET TO REST PHASE OF SUMMER CROPS
            IF (ABS(growth_phase_CRP) > 1.5_wp .AND.      &                             ! winter crops in growth phase..
                ABS(heat_sum_winter) < pheno_param_jsbach%CRP%heat_sum_harvest .AND. &  ! and they are not ready for harvest..
                heat_sum_CRP > pheno_param_jsbach%CRP%heat_sum_harvest) THEN            ! and summer crops are ready for harvest..
                   growth_phase_CRP = -1._wp                                            ! => rest phase of summer crops
            END IF

            ! If there are winter crops and summer crops would have been harvested already twice
            ! --> crop phase is SET TO REST PHASE OF SUMMER CROPS
            IF (ABS(growth_phase_CRP) > 1.5_wp .AND. &                                 ! winter crops in growth phase..
               heat_sum_CRP > 2.0_wp * pheno_param_jsbach%CRP%heat_sum_harvest) THEN   ! and summer crops would have been harvested
                                                                                       ! already twice..
                  growth_phase_CRP = -1._wp                                            ! => rest phase of summer crops
            END IF
            ! If there is no double cropping and summer crops are still not harvested and winter crops would have been harvested
            ! last spring and it is still warm then try winter crops --> crop phase is SET TO REST PHASE OF WINTER CROPS
            ! R: Note, here was an error in the "IF-Tree" in the JSBACH3 version! As I corrected this error there is a difference
            !    between JSBACH3 and 4 here! However, the error only inhibited winter crops, so it wasn't so bad...
            IF (ABS(growth_phase_CRP - 1._wp) > 0.5_wp .AND. &                         ! growth_phase_CRP is not 1 (growth phase sec. crop)
                heat_sum_CRP < pheno_param_jsbach%CRP%heat_sum_harvest .AND. &         ! and summer crops are not ready for harvest..
                ABS(heat_sum_winter) > pheno_param_jsbach%CRP%heat_sum_harvest .AND. & ! and winter crops are ready for harvest..
                pseudo_soil_temp > pheno_param_jsbach%CRP%crit_temp) THEN              ! it is warm enough..
                    growth_phase_CRP = -2._wp                                          ! => growth phase of winter crops
            END IF
            heat_sum_winter = EPSILON(1._wp) ! heat summation for winter crops has to start anyway
        END IF ! is_newday

       ! Each time step:
       ! increase heat sum
       temp_min_max = (previous_day_temp_min + previous_day_temp_max) / 2.0_wp
       IF(temp_min_max > pheno_param_jsbach%CRP%gdd_temp .AND.      & ! it is warm enough to count gdd ..
          NPP_positive) THEN                                          ! and yesterdays productivity of a crop or grass type is
                                                                      ! positive
             heat_sum_CRP = heat_sum_CRP + temp_min_max - pheno_param_jsbach%CRP%gdd_temp ! => add mean day temp. excess above critical
                                                                                          !    temperature to the heat sum for summer crops
             IF (heat_sum_winter > 0._wp) THEN                           ! additionally, if we are in the growing period of
                                                                         ! winter crops
                heat_sum_winter = heat_sum_winter + temp_min_max - pheno_param_jsbach%CRP%gdd_temp  ! => add mean day temperature excess
                                                                                                    !    above critical temperature to
                                                                                                    !    the heat sum for winter crops
             END IF
       END IF

       ! set growth phase for crops
       IF(pseudo_soil_temp > pheno_param_jsbach%CRP%crit_temp) THEN             ! it is warm enough for crop growth ..
          IF (ABS(growth_phase_CRP) < 1.5_wp) THEN                              ! .. and if there are no winter crops ..
             IF (heat_sum_CRP > pheno_param_jsbach%CRP%heat_sum_harvest) THEN   ! .. and the heat sum for summer crops is high enough
                growth_phase_CRP = -1._wp                                       ! => crops are harvested,  i.e. rest phase starts,
             ELSE
                IF (growth_phase_CRP < 0.5_wp) growth_phase_CRP = 0._wp   ! => crops are still in the growth phase ..
             END IF
          ELSE                                                          ! .. but if there are winter crops and
             IF (heat_sum_winter > pheno_param_jsbach%CRP%heat_sum_harvest .OR. &   ! .. the heat sum for winter crops is high enough
                 heat_sum_winter < 0._wp) THEN                                      ! .. or their rest period (summer) started
                    growth_phase_CRP = -2._wp                                       ! => crops are harvested, i.e. rest phase starts,
             ELSE
                growth_phase_CRP = 2._wp                                    ! => crops are still in the growth phase
             END IF
          END IF
       END IF
     END IF

     !!!!!!!! FINALLY, what it is all for:
     ! Update the output variable lai
     IF (ABS(lat) > 30._wp) THEN     ! We are outside the TROPICS

        ! If crops are in the growth phase and there is sufficient water in upper soil layer, then there are good growth
        ! conditions. (Note: addition of 0.001 for numerical reasons)
        IF(growth_phase_CRP > -0.5_wp                           & ! there is ANY growth phase of any crop..
                  .AND.                                         &
           w_soil_filling > pheno_param_jsbach%all%wilt_point + 1E-3_wp & ! soil water is above wilt point
           ) THEN

                ! If also temperature is fine growth should start at least with seed value ( < laiMax) if upper soil layer
                ! has at least a little bit water
                IF(pseudo_soil_temp > pheno_param_jsbach%CRP%crit_temp) THEN ! it is warm enough
                   IF(lai < pheno_param_jsbach%all%laiSeed .AND.           & ! the lai is smaller than lai of seed..
                      w_soil_filling > pheno_param_jsbach%CRP%sproud) THEN   ! and a critical fract. of soil water bucket level is reached
                      lai = MIN(LaiMax_corrected-EPSILON(1.0_wp),pheno_param_jsbach%all%laiSeed) ! => at least with seed value ( < laiMax)
                   ENDIF

                   ! If also yesterdays net production is positive let it grow
                   IF (previous_day_NPP_pot_rate_ca > 0._wp) THEN
                      grRate = pheno_param_jsbach%CRP%leafAlloc_fract * lctlib_specificLeafArea_C * & ! => grRate is high depending on NPP
                                 MAX(0.0_wp,NPP_pot_rate_ca) * 86400.0_wp / MAX(lai,1.E-6_wp)         !    and is in units of 1/days
                      lai = get_letItGrow(timeStep_in_days, lai, grRate, 0.0_wp, LaiMax_corrected)
                   ELSE                                                                ! else ..
                      xtmp = pheno_param_jsbach%CRP%shedRate_growth
                      lai = get_letItDie(timeStep_in_days,lai,xtmp)                           ! do a tiny bit of leaf shedding
                   END IF
                ELSE                                                 ! it is not warm enough
                   xtmp = pheno_param_jsbach%CRP%shedRate_growth
                   lai = get_letItDie(timeStep_in_days,lai,xtmp)                       ! do a tiny bit of leaf shedding
                END IF
        ELSE                                                      ! there is NOT ANY growth phase of any crop
                                                                  ! OR soil water is below wilt point
           xtmp = pheno_param_jsbach%CRP%shedRate_rest
           lai = get_letItDie(timeStep_in_days,lai,xtmp)             ! shed leaves fast
        END IF

     ELSE                              ! We are in the TROPICS
        ! If temperature is fine and there is sufficient water in upper soil layer then there are good growth conditions
        !  (addition of 0.001 for numerical reasons)
        IF (pseudo_soil_temp > pheno_param_jsbach%CRP%crit_temp         & ! it is warm enough..
                  .AND.                                         &
           w_soil_filling > pheno_param_jsbach%all%wilt_point + 1E-3_wp & ! and soil water is above wilt point
          ) THEN
               ! growth should start at least with with seed value ( < laiMax) if upper soil layer has at least a little bit
               ! water
               IF(lai < pheno_param_jsbach%all%laiSeed .AND.             & ! the lai is smaller than lai of seed..
                  w_soil_filling > pheno_param_jsbach%CRP%sproud) THEN     ! and a critical fraction of soil water bucket level is reached
                     lai = MIN(LaiMax_corrected-EPSILON(1.0_wp),pheno_param_jsbach%all%laiSeed)
               ENDIF

               ! If yesterdays net production is positive let it grow
               IF(previous_day_NPP_pot_rate_ca > 0._wp ) THEN
                  grRate = pheno_param_jsbach%CRP%leafAlloc_fract * lctlib_specificLeafArea_C *  &  ! => growth rate
                                MAX(0.0_wp,NPP_pot_rate_ca) * 86400.0_wp / MAX(lai,1.E-6_wp)        !    .. is in units of 1/days
                  lai = get_letItGrow(timeStep_in_days,lai, grRate, 0.0_wp, LaiMax_corrected)
               ELSE
                  xtmp = pheno_param_jsbach%CRP%shedRate_growth
                  lai = get_letItDie(timeStep_in_days,lai,xtmp)     ! .. do a tiny bit of leaf shedding
               END IF
        ELSE                                                     ! it is NOT warm enough..
           xtmp = pheno_param_jsbach%CRP%shedRate_rest
           lai  = get_letItDie(timeStep_in_days,lai,xtmp)      ! .. shed leaves fast
        END IF

     END IF ! outside the TROPICS

   END SUBROUTINE calc_crop_phenology



  SUBROUTINE calc_raingreen_phenology(timeStep_in_days,      & ! <- in
                                                w_soil_filling,        & ! <- in
                                                LaiMax_corrected,      & ! <- in
                                                previous_day_NPP_pot_rate_ca, & ! <- in
                                                lai )                    ! <- inout

    !$ACC ROUTINE SEQ

     !-----------------------------------------------------------------------
     !  DECLARATIONS
     USE mo_pheno_parameters,  ONLY: pheno_param_jsbach

     !-----------------------------------------------------------------------
     !  ARGUMENTS

     REAL(wp), intent(in)    :: timeStep_in_days

     REAL(wp), intent(in)    ::                   &
                                w_soil_filling,   &
                                LaiMax_corrected, &
                                previous_day_NPP_pot_rate_ca

     REAL(wp), intent(inout) :: lai  ! Total LAI.

     !-----------------------------------------------------------------------
     !  LOCAL VARIABLES
     REAL(wp)          :: xtmp,   &       ! argument to calls of get_letItDie
                                          ! (necessary to allow for inlining of these functions)
                          grRate, &       ! argument to calls of get_letItGrow
                          average_filling ! average filling of the soil buckets

     !-----------------------------------------------------------------------
     ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
     ! IF (one_of('calc_phenology', tile%components) < 1) RETURN

     !-----------------------------------------------------------------------
     ! CONTENT

     ! R: Note, as in JSBACH3 all soil layers had the same water content we skiped the soil layers for w_soil_filling in JSBACH4
     !    (see update_phenology). Therefore, the following line is unnecessary.
     !average_filling = SUM(w_soil_filling(j,1:nsoil,i)) / REAL(nsoil,wp)
     average_filling = w_soil_filling

     IF (average_filling > pheno_param_jsbach%all%wilt_point + 1E-3_wp) THEN  ! IF soil water is above wilt point
                                                                              ! (Note: addition of 1E-3_dp for numerical reasons)
        IF (lai < pheno_param_jsbach%all%laiSeed &                              ! IF LAI dropped below seed value..
                   .AND.                 &
            average_filling > pheno_param_jsbach%RG%bucketFill_leafout) THEN     ! and there is sufficient water
               lai = MIN(LaiMax_corrected-EPSILON(1.0),pheno_param_jsbach%all%laiSeed) ! => then start growth at least with seed value
                                                                                       ! but smaller than LAIMax
        END IF
        IF(previous_day_NPP_pot_rate_ca > 0._wp) THEN   ! If yesterdays net production is positive
           grRate = pheno_param_jsbach%RG%growthRate
           xtmp   = get_shedRate_RG(average_filling)
           lai    = get_letItGrow(timeStep_in_days, lai, grRate, xtmp, LaiMax_corrected) ! start growing.
        END IF
     ELSE                                                             ! soil water below wilt point
        xtmp = pheno_param_jsbach%RG%shedRate_drySeason
        lai = get_letItDie(timeStep_in_days,lai, xtmp)       ! shed leaves as fast as possible.
     END IF

   END SUBROUTINE calc_raingreen_phenology


  SUBROUTINE calc_grass_phenology(t_air_in_Celcius,      & ! <- in
                                            timeStep_in_days,      & ! <- in
                                            w_soil_filling,        & ! <- in
                                            LaiMax_corrected,      & ! <- in
                                            previous_day_NPP_pot_rate_ca, & ! <- in
                                            lai )                    ! <- inout

    !$ACC ROUTINE SEQ

     !-----------------------------------------------------------------------
     !  DECLARATIONS
     USE mo_pheno_parameters,  ONLY: pheno_param_jsbach

     !-----------------------------------------------------------------------
     !  ARGUMENTS

     REAL(wp), intent(in)    :: timeStep_in_days
     REAL(wp), intent(in)    ::                   &
                                t_air_in_Celcius, &
                                w_soil_filling,   &
                                LaiMax_corrected, &
                                previous_day_NPP_pot_rate_ca

     REAL(wp), intent(inout) :: lai  ! Total LAI.
     !-----------------------------------------------------------------------
     !  LOCAL VARIABLES
     REAL(wp)          :: xtmp       ! argument to calls of get_letItDie
                                     ! (necessary to allow for inlining of these functions)

     !-----------------------------------------------------------------------
     ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
     ! IF (one_of('calc_phenology', tile%components) < 1) RETURN

     !-----------------------------------------------------------------------
     ! CONTENT

     ! Surprisingly the current concept for grasses is more restrictive than that for raingreens!

     IF (t_air_in_Celcius > pheno_param_jsbach%GRS%crit_temp   & ! IF temperature is fine..
                      .AND.                                    & ! and..
         w_soil_filling > pheno_param_jsbach%all%wilt_point  + 1E-3_wp & ! there is sufficient moisture in upper soil layer
        ) THEN                                                           ! (Note: addition of 1E-3_dp for numerical reasons)

          IF(lai < pheno_param_jsbach%all%laiSeed) THEN
             lai = MIN(LaiMax_corrected-EPSILON(1.0),pheno_param_jsbach%all%laiSeed) ! start growth at least with seed, value, but smaller
                                                                                     ! than laiMax.
          ENDIF
          IF (previous_day_NPP_pot_rate_ca > 0._wp) THEN                    ! If yesterdays net production is positive
             xtmp = pheno_param_jsbach%GRS%growthRate
             lai  = get_letItGrow(timeStep_in_days, lai, xtmp, 0.0_wp, LaiMax_corrected)    ! let it grow
          ELSE
             xtmp = pheno_param_jsbach%GRS%shedRate_growth
             lai  = get_letItGrow(timeStep_in_days, lai, 0.0_wp, xtmp, LaiMax_corrected)    ! no growth only a bit of leaf shedding
          END IF
     ELSE                                                        ! IF temperature is not fine...
        xtmp = pheno_param_jsbach%GRS%shedRate_drySeason
        lai  = get_letItDie(timeStep_in_days, lai,xtmp)             ! let it die as fast as possible
     END IF

   END SUBROUTINE calc_grass_phenology


  ! --- derive_maxLai_according_to_allometric_relationships() ---------------------------------------------------------------------
  !
  ! calculate the maximum LAI for current vegetation biomass and given pft dependant parameters
  !
  SUBROUTINE derive_maxLai_according_to_allometric_relationships( &
                                          & lctlib_alpha_nr_ind,            & ! <- in
                                          & lctlib_beta_nr_ind,             & ! <- in
                                          & lctlib_alpha_leaf,              & ! <- in
                                          & lctlib_beta_leaf,               & ! <- in
                                          & lctlib_clumpinessFactor,        & ! <- in
                                          & lctlib_specificLeafArea_C,      & ! <- in
                                          & c_woods,                        & ! <- in
                                          & current_max_green,              & ! <- InOut
                                          & veg_carbon_at_max_green,        & ! <- InOut
                                          & number_of_individuals,          & ! <- InOut
                                          & biomass_per_individual,         & ! <- InOut
                                          & veg_fract_correction,           & ! <- InOut
                                          & maxLAI_allom )                    ! <- inout

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_pheno_parameters,  ONLY: pheno_param_jsbach
    USE mo_carbon_constants, ONLY: molarMassC_kg

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    REAL(wp), intent(in)    ::  lctlib_alpha_nr_ind,            &
                              & lctlib_beta_nr_ind,             &
                              & lctlib_alpha_leaf,              &
                              & lctlib_beta_leaf,               &
                              & lctlib_clumpinessFactor,        &
                              & lctlib_specificLeafArea_C

    REAL(wp), intent(in)    ::  c_woods
    REAL(wp), intent(inout) ::  current_max_green,              &
                              & veg_carbon_at_max_green,        &
                              & number_of_individuals,          &
                              & biomass_per_individual,         &
                              & veg_fract_correction,           &
                              & maxLAI_allom

    !-----------------------------------------------------------------------
    !  LOCAL VARIABLES     TODO_JN -> KN select those from above which are only required internally
    REAL(wp) :: kgPerHa, &
          & leafC_per_individual,       & ! Leaf C per individual according to allometric relationships [kg/tree]
          & leafC                         ! Leaf C according to allometric relationships [mol(C) m-2(canopy)]


    !-----------------------------------------------------------------------
    ! CONTENT
    kgPerHa = pheno_param_jsbach%FR%kgCtokg * molarMassC_kg * pheno_param_jsbach%FR%m2toha

    !For numerical reasons funct_allometry%FOM_alloc_vegetation_carbon_at_max_green_pool cannot be zero.
    !We also want to set max_LAI to the minimum value after harvest,
    !    ie. when cbalance%c_woods(kidx0:kidx1,itile)=zero.
    IF ((veg_carbon_at_max_green .GT. pheno_param_jsbach%FR%min_c) .AND. (c_woods .GT. pheno_param_jsbach%FR%min_c)) THEN

      !calculate number of individuals, assuming that forest is at self-thinning
      number_of_individuals = exp(MIN((log(veg_carbon_at_max_green * kgPerHa) - lctlib_alpha_nr_ind) &
          & / lctlib_beta_nr_ind, pheno_param_jsbach%FR%log_maxind))

      !calculate biomass per individual
      biomass_per_individual = MAX((veg_carbon_at_max_green * kgPerHa) / number_of_individuals, pheno_param_jsbach%FR%min_zero)

      !calculate Cleaf per individual based on allometric relationship
      leafC_per_individual = exp(lctlib_alpha_leaf + lctlib_beta_leaf * log(biomass_per_individual))

      !calculate C leaf at stand level
      leafC = (leafC_per_individual * number_of_individuals) / kgPerHa

      !derive maximally supported LAI
      maxLAI_allom = MAX(leafC * lctlib_specificLeafArea_C,  pheno_param_jsbach%FR%min_maxLAI)

    ELSE
      maxLAI_allom = pheno_param_jsbach%FR%min_maxLAI
    ENDIF

    ! recalculate correction factor for the canopy cover
    veg_fract_correction = 1.0_wp - exp(-maxLAI_allom / lctlib_clumpinessFactor )

    ! reset tracking variables
    veg_carbon_at_max_green = 0.0_wp
    current_max_green = 0.0_wp

  END SUBROUTINE derive_maxLai_according_to_allometric_relationships

  ! --- track_max_green_pool() ---------------------------------------------------------------------------------------------------
  !
  ! track the maximum of the green pool and the according total vegetation carbon
  !
  SUBROUTINE track_max_green_pool( &
      & c_green,                    & ! <- in
      & c_woods,                    & ! <- in
      & c_reserve,                  & ! <- in
      & current_max_green,          & ! <- InOut
      & veg_carbon_at_max_green )     ! <- inout

    !$ACC ROUTINE SEQ

     !-----------------------------------------------------------------------
     !  ARGUMENTS
     REAL(wp), intent(in)    :: c_green, c_woods, c_reserve
     REAL(wp), intent(inout) :: current_max_green, veg_carbon_at_max_green

     IF (current_max_green + EPSILON(1.0_wp) .LT. c_green) THEN
       current_max_green = c_green
       veg_carbon_at_max_green = c_green + c_woods + c_reserve
     ENDIF

  END SUBROUTINE track_max_green_pool



  ! --- get_letItDie() ----------------------------------------------------------------------------------------------------------------
  !
  ! This is the function get_letItGrow() for growth rate k=0, i.e. it describes only leaf shedding:
  !      x(t+tau) = x(t)*exp(-p*tau),
  ! where "tau" is the time step and "p" the leaf shedding rate. This function is introduced only to save computation time.
  REAL(wp) FUNCTION get_letItDie(timeStep_in_days,x,p)  !! returns x(t+tau)

    REAL(wp),INTENT(in) ::                  &
                          timeStep_in_days, & ! the time step expressed in days (set in "initPhenology")
                          x,                & ! lai at time t
                          p                   ! leaf shedding rate (units: 1/days)

    !$ACC ROUTINE SEQ

    get_letItDie = EXP(-p * timeStep_in_days) * x

  END FUNCTION get_letItDie


  ! --- get_letItGrow() ---------------------------------------------------------------------------------------------------------------
  !
  ! Performs one time step for the growth of the phenology state.
  ! The growth model for the phenology state (x=LAI/LaiMax_dyn) involves two elements: logistic growth and exponential leaf
  ! shedding. Logistic growth guarantees that growth stops at the carrying capacity (here equal 1, because this is the maximum
  ! value of the pheno state). Together with a term accounting for leaf shedding the growth equation is
  !            dx
  ! (1)       ---- = k*(1-x)*x -p*x,
  !            dt
  ! where "k'" is the growth rate scaled by the maximum LAI
  !
  !                                             x(t)
  ! (2)      x(t+tau) = (k-p) ----------------------------------------.
  !                            k*x(t) + exp(-(k-p)*tau)*(k-p-k*x(t))
  !
  REAL(wp) FUNCTION get_letItGrow(timeStep_in_days,x,k,p,z)  !! returns result of integration x(t+tau)

    REAL(wp),INTENT(in) ::                   &
                           timeStep_in_days, & ! the time step length (delta t) expressed in days
                           x,                & ! LAI at time t
                           k,                & ! growth rate          (units: 1/days)
                           p,                & ! leaf shedding rate   (units: 1/days)
                           z                   ! maximum LAI

    REAL(wp)            :: hlp1,hlp2,numerator,denominator

    !$ACC ROUTINE SEQ

    hlp1 = k-p
    hlp2 = k * x / z   ! x/z is the phenoState
    ! Note: make sure that x <= z on input. Call to finish in case of x>z  has been removed from here to allow for inlining.

    numerator = hlp1
    denominator = hlp2 + EXP(-hlp1*timeStep_in_days) * (hlp1 - hlp2)
    IF(ABS(ABS(numerator) - ABS(denominator)) <= 100._wp * EPSILON(1.0_wp)) THEN !! This prevents zero/zero divisions
       get_letItGrow = x
    ELSE
       get_letItGrow = x * numerator / denominator
    ENDIF
  END FUNCTION get_letItGrow


  ! --- shedRate_RG ---------------------------------------------------------------------------------------------------------------
  !
  ! Computes the shedding rate for the raingreens during wet season as a function of the soil water gauge
  ! The idea is, that for high water availability (bucket filling > bucketFill_critical) leafs are shedded only because of their
  ! aging. Below this critical value the shedding rate increases because the plants adapt to lower water availability. At the
  ! wilting point the shedding is so large, that it is larger than the assumed fixed growth rate.
  !
  REAL(wp) FUNCTION get_shedRate_RG(bucket_filling)

    USE mo_pheno_parameters,  ONLY: pheno_param_jsbach

    REAL(wp), INTENT(in) :: bucket_filling     ! average filling of the soil water buckets

    REAL(wp)             :: shedRate, slope

    !$ACC ROUTINE SEQ

    slope    = pheno_param_jsbach%RG%growthRate / &
      &        (pheno_param_jsbach%RG%bucketFill_critical - pheno_param_jsbach%all%wilt_point)
    shedRate = pheno_param_jsbach%RG%growthRate + pheno_param_jsbach%RG%shedRate_aging - &
      &        slope * (bucket_filling - pheno_param_jsbach%all%wilt_point)
    get_shedRate_RG = MIN( MAX( pheno_param_jsbach%RG%shedRate_aging , shedRate) , &
      &                    pheno_param_jsbach%RG%growthRate + pheno_param_jsbach%RG%shedRate_aging)

  END FUNCTION get_shedRate_RG

  REAL(wp) FUNCTION get_foliage_projected_cover(fpc_max, lai, clumpiness)

    REAL(wp), INTENT(in) :: &
      & fpc_max,            &
      & lai,                &
      & clumpiness

    !$ACC ROUTINE SEQ

    get_foliage_projected_cover = fpc_max * (1._wp - EXP(-lai/clumpiness))

  END FUNCTION get_foliage_projected_cover

#endif
END MODULE mo_pheno_process
