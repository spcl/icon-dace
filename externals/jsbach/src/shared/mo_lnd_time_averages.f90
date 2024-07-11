!> defines routine for averaging variables
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
!>#### defines a generic routine providing various ways for the averaging of variables over time for all model processes
!>
MODULE mo_lnd_time_averages
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp


  IMPLICIT NONE
  PUBLIC

  !!----------------------------------------------------------------------------------------------------
  !! parameter definition: moving average period lengths
  !!   you can add new parameters here, please sort by increasing numbers
  !!   the default unit of the averaging period is 'day'
  !!   but can also be 'week' or 'year'
  !!
  ! generic memory time-scales (in days):
    REAL(wp),SAVE :: mavg_period_weekly            = 7.0_wp               !< 7 days; constant across the model
    REAL(wp),SAVE :: mavg_period_monthly           = 30.0_wp              !< 30 days; constant across the model
    REAL(wp),SAVE :: mavg_period_yearly            = 365.0_wp             !< 365 days

  ! specific memory time-scales (in days):
    REAL(wp),SAVE :: mavg_period_tsoa              = 2.0_wp               !< memory time_scale for state of acclimation
    REAL(wp),SAVE :: mavg_period_tuptake           = 3.0_wp               !< nutrient uptake processes (demand, BNF)
    REAL(wp),SAVE :: mavg_period_tlabile           = 7.0_wp               !< labile pool
    REAL(wp),SAVE :: mavg_period_tphen             = 7.0_wp               !< phenological processes
    REAL(wp),SAVE :: mavg_period_tenzyme           = 7.0_wp               !< memory time-scale for enzyme allocation (SB_)
    REAL(wp),SAVE :: mavg_period_tfrac             = 10.0_wp              !< leaf N fractions
    REAL(wp),SAVE :: mavg_period_tcnl              = 20.0_wp              !< leaf stoichiometry
    REAL(wp),SAVE :: mavg_period_talloc            = 30.0_wp              !< biomass allocation
    REAL(wp),SAVE :: mavg_period_tacclim           = 30.0_wp              !< respiration acclimation
    REAL(wp),SAVE :: mavg_period_tmic              = 30.0_wp              !< microbial community acclimation
    REAL(wp),SAVE :: mavg_period_tgrowth           = 365.0_wp             !< vegetation growth
    REAL(wp),SAVE :: mavg_period_tvegdyn           = 365.0_wp             !< vegetation dynamics processes
    REAL(wp),SAVE :: mavg_period_talloc_dynamic    = 365.0_wp * 30.0_wp   !< biomass allocation for empirical dynamic allocation (30 years)
  !!----------------------------------------------------------------------------------------------------


  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_time_averages'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> calculate and return the moving average of a single variable or an array
  !!
  !! this function is used by multiple tasks from several processes
  !!
  !! the default avg_period_length has the unit 'day';
  !!  see predefined averaging-periods 'mavg_period_*' which give the period in days by default \n
  !! the avg_period_unit can be either day (default), week or year
  !!
  !! this function calculates the average across a certain period preceding to the current time step; \n
  !! actually, the average gets rather estimated than exactly calculated using the weighted arithmetic mean, \n
  !!
  !! if the optional parameter 'do_calc' is present: the moving average is calculated only if do_calc = TRUE !
  !!
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_time_mavg(current_avg, new_value, avg_period_length, &   ! default  input arguments
                                         do_calc, avg_period_unit) &                    ! optional input arguments
                                         RESULT(new_avg)

    USE mo_kind,                ONLY: wp
    USE mo_jsb_math_constants,  ONLY: dtime, one_day, one_year

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in)                      :: current_avg        !< current average value
    REAL(wp), INTENT(in)                      :: new_value          !< new value
    REAL(wp), INTENT(in)                      :: avg_period_length  !< length of the period over which the variable is averaged
    LOGICAL,          OPTIONAL, INTENT(in)    :: do_calc            !< optional, do calculate new_avg from the new_value T/F
    CHARACTER(LEN=*), OPTIONAL, INTENT(in)    :: avg_period_unit    !< optional, for different units than 'day' for avg_period length
    REAL(wp)                                  :: new_avg            !< function output: the new average
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                                  :: new_value_weight   ! weight of the new value
    REAL(wp)                                  :: avg_period2seconds ! for converting the period length (day, week, year) into seconds
    LOGICAL                                   :: l_do_calc          ! True/ False of the Optional variable do_calc present
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//':calc_time_mavg'


    !> 1.0 deal with the optional arguments
    !!
    !> 1.1 do_calc: set TRUE per default, i.e. can only be FALSE if 'do_calc' is is present
    !!
    l_do_calc = .TRUE.
    ! check whether the optional statement is provided
    IF(PRESENT(do_calc))  l_do_calc = do_calc

    !> 1.2 avg_period_unit: calc period length in seconds
    !!
    avg_period2seconds = one_day  ! (one_day == 86400._wp seconds)
    IF(l_do_calc .AND. PRESENT(avg_period_unit)) THEN
      SELECT CASE(TRIM(avg_period_unit))
      !CASE('day')
      ! nothing to do, as 'avg_period2seconds = one_day' is set by default anyway
      CASE('week')
        avg_period2seconds = avg_period2seconds * 7.0_wp
      CASE('year')
        avg_period2seconds = avg_period2seconds * one_year ! (one_year == 365.0_wp days)
      END SELECT  ! select SLM model !
    END IF

    !> 2.0 calc weight of the new_value and calc the new_avg
    !! calc the new average as weighted mean of the current_avg and the new_value weighted by new_value_weight
    IF(l_do_calc) THEN
      new_value_weight = 1.0_wp / (avg_period_length * avg_period2seconds / dtime )  ! dtime: length of one time step in seconds
      new_avg = (current_avg * (1.0_wp - new_value_weight)) + (new_value * new_value_weight)
    ELSE
      new_avg = current_avg
    END IF
  END FUNCTION calc_time_mavg

#endif
END MODULE mo_lnd_time_averages
