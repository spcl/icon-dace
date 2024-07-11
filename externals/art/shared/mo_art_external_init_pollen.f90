!
! mo_art_external_init_pollen
! This module provides lookup tables for the modal parameters of pollen
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

MODULE mo_art_external_init_pollen
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_exception,                     ONLY: message
  USE mo_parallel_config,               ONLY: nproma
  USE mo_run_config,                    ONLY: dtime
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH
! ART
  USE mo_art_external_types,            ONLY: t_art_pollen_table, &
    &                                         t_art_pollen_properties
  USE mo_art_emission_pollen_atab,      ONLY: art_pollen_read_atab

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_external_init_pollen'

  PUBLIC :: art_extinit_pollen
  PUBLIC :: art_extinit_pollen_fordiag

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_pollen_fordiag(p_patch,nblks, p_pollen_prop, dict_tracer, &
  &                                   cart_input_folder)
!<
! SUBROUTINE art_extinit_pollen_fordiag
! This SR sets the minimal information needed to construct the pollen
! diagnostics. Necessary since art_extinit_pollen can (currently) only be called
! after diagnostics are created.
! Based on: A. Pauling (2011)
! Part of Module: mo_art_external_init_pollen
! Author: Sven Werchner, KIT
! Initial Release: 2020-10-19
! Modifications:
! DD-MM-YYYY: <author_name>,<institute>
!>
  TYPE(t_patch),INTENT(in)                     :: &
    &  p_patch
  INTEGER, INTENT(in)                          :: &
    &  nblks                  !< number of blocks
  TYPE(t_art_pollen_properties), INTENT(inout) :: &
    &  p_pollen_prop          !< pollen properties
  TYPE(t_key_value_store), INTENT(in)          :: &
    &  dict_tracer            !< Tracer index dictionary
  CHARACTER(LEN=*), INTENT(in)                 :: &
    &  cart_input_folder
  !Local variables
  CHARACTER(len=MAX_CHAR_LENGTH),ALLOCATABLE :: &
    &  vname(:)
  INTEGER         :: &
    &  n,itr,iv,     &        !< counter
    &  ierr                   !< error code
  TYPE(t_art_pollen_table),POINTER :: &
    &  this_pollen_table      !< pointer to pollen table

  ALLOCATE( p_pollen_prop%pollen_type(p_pollen_prop%npollen_types) )
  !$ACC ENTER DATA CREATE(p_pollen_prop%pollen_type)
  ALLOCATE( vname(p_pollen_prop%npollen_types) )

  ! Create a storage container: Table index (integer) --> PollenType(String)
  CALL p_pollen_prop%dict_pollen%init(.FALSE.)

  ! implemented pollen types
  vname = (/'pollbetu','pollpoac','pollambr','pollalnu', 'pollcory'/)

  n = 0 ! counter for pollen_table; storage index

  DO iv = 1,SIZE(vname)

    CALL dict_tracer%get(TRIM(vname(iv)),itr,ierr)
    IF (ierr == SUCCESS) THEN

      n = n + 1
      this_pollen_table => p_pollen_prop%pollen_type(n)

      this_pollen_table%sname = TRIM(vname(iv))
      ALLOCATE( this_pollen_table%fr_cov(nproma, nblks) ) ! used as init-flag in diag-create

      !$ACC ENTER DATA CREATE(this_pollen_table%fr_cov)
      SELECT CASE(TRIM(vname(iv)))
        CASE('pollbetu')
          this_pollen_table%shortname        = 'BETU'
        CASE('pollambr')
          this_pollen_table%shortname        = 'AMBR'
        CASE('pollpoac')
          this_pollen_table%shortname        = 'POAC'
        CASE('pollalnu')
          this_pollen_table%shortname        = 'ALNU'
        CASE('pollcory')
          this_pollen_table%shortname        = 'CORY'
        CASE DEFAULT
          ! cannot be reached
      END SELECT
    ELSE
      CALL message(TRIM(routine)//':art_extinit_pollen','FAILED: Init datatype for: '//TRIM(vname(iv)))
    ENDIF
  END DO
  p_pollen_prop%npollen_used = n

END SUBROUTINE art_extinit_pollen_fordiag
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_pollen(p_patch,nblks, p_pollen_prop, dict_tracer, cart_input_folder)
!<
! SUBROUTINE art_extinit_pollen
! This SR sets the values for the modal properties of different pollen types.
! Based on: A. Pauling (2011)
! Part of Module: mo_art_external_init_pollen
! Author: Jonas Straub, KIT
! Initial Release: 2016-12-21
! Modifications:
! 08-02-2017: Jonas Straub. KIT
! - tuned parameter list: only p_patch%nblks instead of p_patch
! - moved read ext.data to mo_art_external_state
!>
  TYPE(t_patch),INTENT(in)                     :: &
    &  p_patch
  INTEGER, INTENT(in)                          :: &
    &  nblks                  !< number of blocks
  TYPE(t_art_pollen_properties), INTENT(inout) :: &
    &  p_pollen_prop          !< pollen properties
  TYPE(t_key_value_store), INTENT(in)          :: &
    &  dict_tracer            !< Tracer index dictionary
  CHARACTER(LEN=*), INTENT(in)                 :: &
    &  cart_input_folder

  !local
  INTEGER         :: &
    &  itr,ip,       &        !< counter
    &  ierr                   !< error code
  TYPE(t_art_pollen_table),POINTER :: &
    &  this_pollen_table      !< pointer to pollen table
  CHARACTER(len=MAX_CHAR_LENGTH),ALLOCATABLE :: &
    &  vname(:)

  IF (p_pollen_prop%npollen_used==0) THEN
    CALL message('mo_art_external_init_pollen:art_extinit_pollen',  &
      &          'no pollen tracer declared in XML, but iart_pollen > 0')
  END IF

  DO ip = 1, p_pollen_prop%npollen_used

    this_pollen_table => p_pollen_prop%pollen_type(ip)

    !!!! ALLOCATION  !!!!!
    ! fraction of plant depending surface coverage
    ALLOCATE( this_pollen_table%f_q_alt(     nproma, nblks) )
    ! for blooming model
    ALLOCATE( this_pollen_table%tthrs_red(      nproma, nblks ) )
    ! tuning factor derived quantities
    ALLOCATE( this_pollen_table%no_max_day(     nproma, nblks ) )
    ALLOCATE( this_pollen_table%no_max_timestep(nproma, nblks ) )
    !$ACC ENTER DATA CREATE(this_pollen_table%f_q_alt, this_pollen_table%no_max_day, this_pollen_table%no_max_timestep)
    IF(TRIM(this_pollen_table%sname) == 'pollambr') THEN
      ALLOCATE(this_pollen_table%sobs_sum(  nproma, nblks ))
      ALLOCATE(this_pollen_table%rh_sum(    nproma, nblks ))
      !$ACC ENTER DATA CREATE(this_pollen_table%rh_sum, this_pollen_table%sobs_sum)
    END IF


    !!!!  INITIALIZATION  !!!!
    CALL message(TRIM(routine)//':art_extinit_pollen','Init datatype for: '//TRIM(this_pollen_table%sname))

    SELECT CASE(TRIM(this_pollen_table%sname))

    ! BETU (betula), birch
    CASE('pollbetu')
      ! Maximum number of pollen that can be produced on one m2 during one day.
      this_pollen_table%no_max_day(:,:)       = 1.0166E7_wp * this_pollen_table%tune(:,:)
      ! Loss of pollen from the reservoir due to random processes (Animals, ...) per timestep:
      ! Half-life of 43200 seconds (= 12 hours) when only random processes exist.
      this_pollen_table%psi_random       = EXP( LOG(0.5_wp) * dtime / 43200._wp )
      ! Suppression of emission in the aftermath of precipitation:
      ! Evaporation set to 0 when rel. hum. is 100%
      ! Drying needs x hours when rel. hum. is y%
      ! Leads to the formula: coeff = x * (1 - y)
      ! Here: coeff = 3 * (1 - 0.3) = 2.1
      this_pollen_table%xi_r_precip      = 1._wp !if set to 0, no suppression of emission after
                                                 !             precipitation
      this_pollen_table%frac_xi_evap     = this_pollen_table%xi_r_precip * dtime  &
        &                                     / (2.1_wp * 3600._wp)
      ! Maximum number of pollen is reached after 16 h under ideal conditions.
      this_pollen_table%no_max_timestep(:,:)  = this_pollen_table%no_max_day(:,:) * dtime  &
        &                                    / (16._wp*3600._wp)

      ! jul_days_excl defines the start day, from which the model starts to calculate
      ! the temperature sum for birch (temperature sum threshold [tthrs] is derived from pollen data,
      ! by defining the start of the pollen season by number of pollen > 30)
      ! CAUTION: Day 1 is 1.Dec!
      this_pollen_table%jul_days_excl    = 40
      ! pollen-specific base temperature for temperature sum (in deg. C) only days with a
      ! temperature > t_base are considered in the calculation of the temperature sum
      this_pollen_table%t_base           = 9._wp
      this_pollen_table%linit            = .TRUE. ! successful initialized

    ! AMBR (ambrosia), ragweed
    CASE('pollambr')
      ! Tuning based on daily measurements 2012 and 2013 in France, close to ragweed plants.
      ! Maximum number of pollen that can be produced on one m2 during one day.
      this_pollen_table%no_max_day(:,:)       = 0.07535E4_wp * this_pollen_table%tune(:,:)
      ! Loss of pollen from the reservoir due to random processes (Animals, ...) per timestep:
      this_pollen_table%psi_random       = EXP( LOG(0.5_wp) * dtime / 811._wp )
      ! Suppression of emission in the aftermath of precipitation:
      ! Evaporation set to 0 when rel. hum. is 100%
      ! Drying needs x hours when rel. hum. is y%
      ! Leads to the formula: coeff = x * (1 - y)
      ! Here: coeff = 2 * (1 - 0.3) = 1.4
      this_pollen_table%xi_r_precip      = 1._wp
      this_pollen_table%frac_xi_evap     = this_pollen_table%xi_r_precip * dtime  &
        &                                    / (1.4_wp * 3600._wp)
      ! Maximum number of pollen is reached after 9h under ideal conditions.
      this_pollen_table%no_max_timestep(:,:)  = this_pollen_table%no_max_day(:,:) * dtime  &
        &                                    / (9._wp * 3600._wp)

      ! no temperature sum for ambr, therefore not jul_days_excl not needed
      ! the season parametrization (sdes) is based on a climatology.
      this_pollen_table%jul_days_excl    = 1000  !dummy value

      ! sum of radiation since midnight (ambrosia emissions)
      this_pollen_table%sobs_sum(:,:) = 0._wp
      this_pollen_table%rh_sum(:,:) = 0._wp

      this_pollen_table%linit            = .TRUE. ! successful initialized

    ! POAC (poaceae), grasses
    CASE('pollpoac')
      ! Maximum number of pollen that can be produced on one m2 during one day.
      this_pollen_table%no_max_day(:,:)       = 1.09E5_wp * this_pollen_table%tune(:,:)
      ! Maximum number of pollen is reached after 16h under ideal conditions.
      this_pollen_table%no_max_timestep(:,:)  = this_pollen_table%no_max_day(:,:) * dtime  &
        &                                    / (16._wp * 3600._wp)
      ! Loss of pollen from the reservoir due to random processes (Animals, ...) per timestep:
      ! Half-life of 811 seconds (Estimation: 5% per timestep = 60 Sec) when only random
      ! processes exist. 
      this_pollen_table%psi_random       = EXP( LOG(0.5_wp) * dtime / 811._wp )
      ! Suppression of emission in the aftermath of precipitation:
      ! Evaporation set to 0 when rel. hum. is 100%
      ! Drying needs x hours when rel. hum. is y%
      ! Leads to the formula: coeff = x * (1 - y)
      ! Here: coeff = 2 * (1 - 0.3) = 1.4. This is the parameter to tune. Values
      this_pollen_table%xi_r_precip      = 1._wp
      this_pollen_table%frac_xi_evap     = this_pollen_table%xi_r_precip * dtime  &
        &                                    / (1.4_wp * 3600._wp)

      ! jul_days_excl defines the start day, from which the model starts to calculate
      ! the temperature sum for poac (temperature sum threshold [tthrs] is derived from pollen data,
      ! by defining the start of the pollen season by number of pollen >= 20)
      ! CAUTION: Day 1 is 1.Dec!
      this_pollen_table%jul_days_excl    = 46

      ! pollen-specific base temperature for temperature sum (in deg. C) only days with a
      ! temperature > t_base are considered in the calculation of the temperature sum
      this_pollen_table%t_base           = 3._wp
      this_pollen_table%linit            = .TRUE. ! successful initialized

    ! ALNU (alnus), alder
    CASE('pollalnu')
      ! overall tuning factor for Zink emission formula. this is the tuning factor for ALNU.
      ! Maximum number of pollen that can be produced on one m2 during one day.
      this_pollen_table%no_max_day(:,:)       = 9.8E5_wp * this_pollen_table%tune(:,:)
      ! Maximum number of pollen is reached after 16 h under ideal conditions.
      this_pollen_table%no_max_timestep(:,:)  = this_pollen_table%no_max_day(:,:) * dtime  &
        &                                    / (16._wp*3600._wp)
      ! Loss of pollen from the reservoir due to random processes (Animals, ...)
      ! per timestep:
      ! Half-life of 43200 seconds (= 12 hours) when only random processes
      ! exist.
      this_pollen_table%psi_random       = EXP( LOG(0.5_wp) * dtime / 43200._wp )

      ! Suppression of emission in the aftermath of precipitation:
      ! Evaporation set to 0 when rel. hum. is 100%
      ! Drying needs x hours when rel. hum. is y%
      ! Leads to the formula: coeff = x * (1 - y)
      ! Here: coeff = 3 * (1 - 0.3) = 2.1
      this_pollen_table%xi_r_precip      = 1._wp !if set to 0, no suppression of emission after precipitation
      this_pollen_table%frac_xi_evap     = this_pollen_table%xi_r_precip * dtime  &
        &                                    / (2.1_wp * 3600._wp)

      ! jul_days_excl defines the start day, from which the model starts to calculate
      ! the temperature sum for alnu (temperature sum threshold [tthrs] is derived from pollen data,
      ! by defining the start of the pollen season by number of pollen >= 5)
      ! CAUTION: Day 1 is 1.Dec!
      this_pollen_table%jul_days_excl    = 14

      ! pollen-specific base temperature for temperature sum (in deg. C) only days with a
      ! temperature > t_base are considered in the calculation of the temperature sum
      this_pollen_table%t_base           = 5.2_wp
      this_pollen_table%linit            = .TRUE. ! successful initialized

    ! CORY (Corylus), hazel
    CASE('pollcory')
      ! Maximum number of pollen that can be produced on one m2 during one day.
      this_pollen_table%no_max_day(:,:)       = 9.8E5_wp * this_pollen_table%tune(:,:)
      ! Maximum number of pollen is reached after 16 h under ideal conditions.
      this_pollen_table%no_max_timestep(:,:)  = this_pollen_table%no_max_day(:,:) * dtime  &
        &                                    / (16._wp*3600._wp)
      ! Loss of pollen from the reservoir due to random processes (Animals, ...)
      ! per timestep:
      ! Half-life of 43200 seconds (= 12 hours) when only random processes
      ! exist.
      ! was: Psi_others(lalnu) = log(2._wp) * dtime / (12._wp * 3600._wp)
      this_pollen_table%psi_random       = EXP( LOG(0.5_wp) * dtime / 43200._wp )

      ! Suppression of emission in the aftermath of precipitation:
      ! Evaporation set to 0 when rel. hum. is 100%
      ! Drying needs x hours when rel. hum. is y%
      ! Leads to the formula: coeff = x * (1 - y)
      ! Here: coeff = 3 * (1 - 0.3) = 2.1
      this_pollen_table%xi_r_precip      = 1._wp !if set to 0, no suppression of emission after precipitation
      this_pollen_table%frac_xi_evap     = this_pollen_table%xi_r_precip * dtime  &
        &                                    / (2.1_wp * 3600._wp)

      ! jul_days_excl defines the start day, from which the model starts to calculate
      ! the temperature sum for cory (temperature sum threshold [tthrs] is derived from pollen data,
      ! by defining the start of the pollen season by number of pollen >= 5)
      ! CAUTION: Day 1 is 1.Dec!
      this_pollen_table%jul_days_excl    = 3

      ! pollen-specific base temperature for temperature sum (in deg. C) only days with a
      ! temperature > t_base are considered in the calculation of the temperature sum
      this_pollen_table%t_base           = 4.3_wp
      this_pollen_table%linit            = .TRUE. ! successful initialized

    CASE DEFAULT
      ! cannot be reached
    END SELECT

    !$ACC UPDATE DEVICE(this_pollen_table%frac_xi_evap) &
    !$ACC   DEVICE(this_pollen_table%psi_random, this_pollen_table%xi_r_precip)

    this_pollen_table%tthrs_red(:,:)   = 0._wp

    ! fill pollen dictionary. Link in both ways: pollen name to table index
    CALL p_pollen_prop%dict_pollen%put(TRIM(this_pollen_table%sname),ip)
    !p_pollen_prop%dict_pollen%put(n,this_pollen_table%sname)

    ! fill meta_data from atab-file (needed for calc_saisl)
    CALL art_pollen_read_atab(TRIM(this_pollen_table%sname),cart_input_folder,this_pollen_table%pol_atab,p_patch)

  ENDDO

END SUBROUTINE art_extinit_pollen
!!
END MODULE mo_art_external_init_pollen
