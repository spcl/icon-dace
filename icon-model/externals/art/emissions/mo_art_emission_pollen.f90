!
! mo_art_emission_pollen
! This module provides the emission routine of pollen grains
! Based on Zink et al. - EMPOL 1.0: a new parameterization of pollen emission
!                        in numerical weather prediction models.
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

MODULE mo_art_emission_pollen
! ICON
  USE mo_kind,                          ONLY: wp
  USE mtime,                            ONLY: datetime, getDayOfYearFromDateTime
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_impl_constants,                ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_exception,                     ONLY: finish
  USE mo_math_constants,                ONLY: deg2rad
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_util_mtime,                    ONLY: getElapsedSimTimeInSeconds
  USE mo_fortran_tools,                 ONLY: assert_acc_device_only

  ! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_external_types,            ONLY: t_art_pollen_properties,t_art_pollen_table
  USE mo_art_emission_pollen_atab,      ONLY: t_art_pollen_atab_header_info, &
                                          &   t_art_all_stns,t_art_pol_atab

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_emiss_pollen, art_pollen_get_nstns, art_prepare_tsum, art_prepare_saisl, art_prepare_sdes

  CHARACTER(LEN=MAX_CHAR_LENGTH) :: thisroutine = "mo_art_emission_pollen"

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emiss_pollen( p_dtime, current_date, p_rho,                     &
    &                        mode_name, p_pollen_prop, p_tracer, dict_tracer,  &
    &                        p_temp, p_tke,                                    &
    &                        p_rain_gsp_rate, p_rain_con_rate, p_rh_2m,        &
    &                        p_swflxsfc,                                       &
    &                        p_dz_nlev,                                        &
    &                        p_llsm_atm_c,                                     &
    &                        jb, istart, iend, lacc )

!  USE mo_nonhydro_state,       ONLY: p_nh_state

!  TYPE(t_art_pollen_table), INTENT(INOUT) ::  &
!    &   pollen_data                                   !< All data regarding pollen type
!                                                      !  emissions
  TYPE(t_art_pollen_properties), INTENT(INOUT) :: &
    &  p_pollen_prop                                   !< pollen properties

  TYPE(datetime), INTENT(IN)  :: &
    &  current_date                                   !< Date and time information

  TYPE(t_key_value_store), INTENT(IN) :: &
    &  dict_tracer

  REAL(wp),INTENT(IN)::                       &
    &                  p_dtime,               &       !< time step
    &                  p_rho(:),              &       !< density of air   (kg/m^3)
    &                  p_temp(:),             &       !< air temperature  (K)
    &                  p_rain_gsp_rate(:),    &       !< grid-scale rain rate     (kg/(m^2 s))
    &                  p_rain_con_rate(:),    &       !< convective rain rate     (kg/(m^2 s))
    &                  p_rh_2m(:),            &       !< relative humidity at surfcae (%)
    &                  p_dz_nlev(:),          &       !< geometric height of gridbox nlev
    &                  p_tke(:),              &       !< turbulent kinetic energy  (m^2/s^2)
    &                  p_swflxsfc(:)                  !< short-wave radiation flux at the surface (W/m^2)

  REAL(wp),INTENT(INOUT)::                &
    &                  p_tracer(:,:)                  !< concentration of the pollen species

  INTEGER,INTENT(IN)          ::       &
    &     jb,                    &                    !< block loop index
    &     istart, iend                                !< start and end of nproma

  LOGICAL,INTENT(IN)          ::       &
    &     p_llsm_atm_c(:)                             !< TRUE, if cell-center is landpoint

  CHARACTER(LEN=IART_VARNAMELEN),INTENT(IN) :: &
    &  mode_name                                      !< Mode name,

  LOGICAL, OPTIONAL, INTENT(IN) :: lacc

! Local Variables
  INTEGER            :: &
    &  jc,              & !< counter for nproma loop
    &  ipoll_tr,        & !< Index of pollen in tracer container
    &  ipoll,           & !< Index of pollen in pollen table
    &  ierror,          & !< Error return value of tracer dictionary
    &  elapsed_seconds    !< Elapsed seconds since simulation start
  REAL(wp)           :: &
    &  phi_biol,        & !< biological influences
    &  phi_met,         & !< meteorological influences (step 1)
    &  phi_plant,       & !< incorporates the fraction of plant
    &  psi_precip,      & !< rain emptying the reservoir
    &  psi_emis,        & !< pollen grains emitted from the reservoir (m-3)
    &  psi_wet,         & !< Suppression of emission in the aftermath of precipitation
    &  delta_r_precip,  & !< amount of water evaporating from r_precip during one ti
    &  c_max,           & !< highest concentration possible (1/m3)
    &  c_quell,         & !< quelled pollen concentration   (1/m3)
    &  c_buff,          & !< buffer for unit transformation (kg/m3)
    &  f_e_tke,         & !< influence of wind speed
    &  max_emiss_day,   & !< maximum of daily emission per m2
    &  precip,          & !< absolute rain rate at surface (kg/(m2 s))
    &  time_sec,        & !< current time in seconds, for checking the start of a new day
    &  hour,            & !< current hour
    &  minute,          & !< current minute
    &  second,          & !< current second
    &  nn11,            & !< neural network layer 1 node 1 (Ambrosia emission)
    &  nn12,            & !< neural network layer 1 node 2 (Ambrosia emission)
    &  nn13,            &
    &  nn14,            &
    &  nn15,            &
    &  nn16,            &
    &  nn17,            &
    &  nn18,            &
    &  nn19,            &
    &  nn110,           &
    &  nn21,            & ! neural network layer 2 node 1
    &  nn22,            &
    &  nn23,            &
    &  nn24,            &
    &  nn25,            &
    &  nn26,            &
    &  nn27,            &
    &  nn28,            &
    &  nn29,            &
    &  nn210

  REAL(wp) :: &
    &  rh(istart:iend),       & !< relative humidity  (1)
    &  f_r_rh(istart:iend),   & !< influence of relative humidity
    &  f_r_t(istart:iend),    & !< influence of temperature
    &  res_new(istart:iend),  & !< Number of pollen released into the reservoir  (1)
    &  f_e_rh(istart:iend)      !< influence of moisture on the foliage

  TYPE(t_art_pollen_table),POINTER  ::  &
    &   pollen_data       !< All data regarding one specific pollen type

  CALL assert_acc_device_only('art_emiss_pollen', lacc)


  !-----------------------------------------------------------------------------------------
  !--   Start Routine ----------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  CALL dict_tracer%get(TRIM(mode_name), ipoll_tr, ierror)
    IF(ierror /= SUCCESS) CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
                            &          'tracer "'//TRIM(mode_name)//'" not found in tracer dictionary.')
  CALL p_pollen_prop%dict_pollen%get(TRIM(mode_name), ipoll, ierror)
    IF(ierror /= SUCCESS) CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
                            &          '"'//TRIM(mode_name)//'" not found in pollen table dictionary.')
  pollen_data => p_pollen_prop%pollen_type(ipoll)

  hour   = REAL(current_date%time%hour,wp)
  minute = REAL(current_date%time%minute,wp)
  second = REAL(current_date%time%second,wp)

  ! ----------------------------------
  ! --- Calculate and add the emission flux
  ! ----------------------------------

  !$ACC DATA CREATE(f_e_rh, f_r_rh, f_r_t, res_new, rh) &
  !$ACC   PRESENT(pollen_data, pollen_data%fe_plant, pollen_data%fr_cov, pollen_data%f_q_alt, pollen_data%f_q_seas) &
  !$ACC   PRESENT(pollen_data%res_new_sum, pollen_data%res_old, pollen_data%rh_sum, pollen_data%r_precip) &
  !$ACC   PRESENT(pollen_data%sobs_sum, p_dz_nlev, p_llsm_atm_c, p_rain_con_rate, p_rain_gsp_rate, p_rho, p_rh_2m) &
  !$ACC   PRESENT(p_swflxsfc, p_temp, p_tke, p_tracer)

  SELECT CASE (TRIM(mode_name))

    CASE ('pollbetu')  ! Formulas for birch

!NEC$ ivdep
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = istart, iend  !<Loop over all horizontal grid points
        IF( p_llsm_atm_c(jc) ) THEN
          !---------------------------------------------------------------------!
          ! Preliminary calculation: rh                                         !
          !---------------------------------------------------------------------!
          rh(jc) = p_rh_2m(jc) * 0.01_wp

          !...........................................................
          ! Meteorological influences:
          !...........
          f_r_rh(jc) = 1._wp / (1._wp + EXP(21.00_wp * rh(jc) - 15._wp))

          f_r_t(jc) = 1.04_wp * (1._wp / (1._wp + EXP(-0.27_wp * p_temp(jc) +  76._wp))) *  &
            &                   (1._wp / (1._wp + EXP( 0.45_wp * p_temp(jc) - 137._wp)))


        ENDIF ! p_llsm_atm_c
      ENDDO ! jc
      !$ACC END PARALLEL

    CASE ('pollambr')  ! Formulas for ragweed

      ! get elapsed simulation tim in seconds
      elapsed_seconds = getElapsedSimTimeInSeconds(current_date)           

!NEC$ ivdep
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR &
      !$ACC   PRIVATE(nn11, nn12, nn13, nn14, nn15, nn16, nn17, nn18, nn19, nn110) &
      !$ACC   PRIVATE(nn21, nn22, nn23, nn24, nn25, nn26, nn27, nn28, nn29, nn210)
      DO jc = istart, iend  !<Loop over all horizontal grid points
        IF( p_llsm_atm_c(jc) ) THEN
          !---------------------------------------------------------------------!
          ! Preliminary calculation: rh                                         !
          !---------------------------------------------------------------------!
          rh(jc) = p_rh_2m(jc) * 0.01_wp

          !...........................................................
          ! Meteorological influences:
          !...........
          !f_r_rh(jc) = 1._wp / (1._wp + EXP(20._wp *  rh(jc) - 12._wp)) !orig
          !f_r_t(jc)  = 1._wp / (1._wp + EXP(- p_temp(jc) * 0.2_wp + 60._wp)) !orig

          !ANN model 5 of Burki et al. (2019): predictors RH_sum, sobs_sum
          !The neural network was trained on 2014, 2015 and 2016 POLEMIC data (Serbia)

          ! sobs_sum and rh_sum are incremented at or near secs 
          !"0000","0750","1500","2250","3000","3750","4500","5250"
          IF (MODULO(elapsed_seconds, 450) < INT(p_dtime)) THEN
            ! FIXME: why is 102000._wp needed ????
            pollen_data%sobs_sum(jc,jb) = pollen_data%sobs_sum(jc,jb) + p_swflxsfc(jc)/102000._wp
            pollen_data%rh_sum(jc,jb)   = pollen_data%rh_sum(jc,jb)   + (1._wp - rh(jc))/100._wp
          ENDIF

          nn11=TANH(-0.78277_wp-1.3121_wp*pollen_data%rh_sum(jc,jb)+2.41133_wp*pollen_data%sobs_sum(jc,jb))
          nn12=TANH(0.03143_wp-0.82722_wp*pollen_data%rh_sum(jc,jb)+0.8307_wp*pollen_data%sobs_sum(jc,jb))
          nn13=TANH(-0.59948_wp+0.39919_wp*pollen_data%rh_sum(jc,jb)+2.19072_wp*pollen_data%sobs_sum(jc,jb))
          nn14=TANH(-1.74111_wp-0.35065_wp*pollen_data%rh_sum(jc,jb)-0.01616_wp*pollen_data%sobs_sum(jc,jb))
          nn15=TANH(-1.30666_wp+0.86516_wp*pollen_data%rh_sum(jc,jb)-1.08275_wp*pollen_data%sobs_sum(jc,jb))
          nn16=TANH(3.90802_wp+211.08407_wp*pollen_data%rh_sum(jc,jb)-98.20912_wp*pollen_data%sobs_sum(jc,jb))
          nn17=TANH(0.32123_wp+0.75077_wp*pollen_data%rh_sum(jc,jb)+1.49074_wp*pollen_data%sobs_sum(jc,jb))
          nn18=TANH(-0.54378_wp+0.79411_wp*pollen_data%rh_sum(jc,jb)-1.12405_wp*pollen_data%sobs_sum(jc,jb))
          nn19=TANH(-0.1206_wp-0.56818_wp*pollen_data%rh_sum(jc,jb)+1.54418_wp*pollen_data%sobs_sum(jc,jb))
          nn110=TANH(0.29761_wp-0.4393_wp*pollen_data%rh_sum(jc,jb)+0.10877_wp*pollen_data%sobs_sum(jc,jb))

          nn21=TANH(0.5584_wp-0.15288_wp*nn11+0.4372_wp*nn12+1.47799_wp*nn13 &
            &   -1.04862_wp*nn14+1.0169_wp*nn15-0.50665_wp*nn16-1.01584_wp*nn17 &
            &   +1.14501_wp*nn18+0.34865_wp*nn19+1.3187_wp*nn110)

          nn22=TANH(0.58017_wp+1.49962_wp*nn11-1.12488_wp*nn12+1.05889_wp*nn13 &
            &    +0.43026_wp*nn14-0.763_wp*nn15+1.51859_wp*nn16+0.0287_wp*nn17 &
            &    +0.83981_wp*nn18-0.79095_wp*nn19-0.47969_wp*nn110)

          nn23=TANH(-0.54507_wp-1.78196_wp*nn11+1.66104_wp*nn12-1.51524_wp*nn13 &
            &   +0.10023_wp*nn14+0.65049_wp*nn15+1.53689_wp*nn16-0.66237_wp*nn17 &
            &   -0.43884_wp*nn18+0.87595_wp*nn19-0.3894_wp*nn110)

          nn24=TANH(-0.14949_wp-0.66903_wp*nn11+53.1201_wp*nn12-3.56556_wp*nn13 &
            &    +1.20451_wp*nn14+0.13177_wp*nn15+2.91022_wp*nn16-23.59427_wp*nn17 &
            &    +1.31307_wp*nn18-29.14264_wp*nn19+0.91035_wp*nn110)

          nn25=TANH(-0.60736_wp+0.82775_wp*nn11-1.88571_wp*nn12-2.24339_wp*nn13 &
            &    +1.2093_wp*nn14+0.26426_wp*nn15-1.41972_wp*nn16+2.19654_wp*nn17 &
            &    -0.35508_wp*nn18-0.89001_wp*nn19-1.26121_wp*nn110)

          nn26=TANH(2.12732_wp-0.4517_wp*nn11+2.30073_wp*nn12-0.51188_wp*nn13 &
            &    -0.38535_wp*nn14-0.22083_wp*nn15-1.18406_wp*nn16+1.14319_wp*nn17 &
            &    -0.16739_wp*nn18+0.08313_wp*nn19-0.34594_wp*nn110)

          nn27=TANH(-0.44247_wp-0.63948_wp*nn11-0.80694_wp*nn12+0.24845_wp*nn13 &
            &    +0.44157_wp*nn14-1.04759_wp*nn15-0.41582_wp*nn16-0.77383_wp*nn17 &
            &    -0.27629_wp*nn18-0.52017_wp*nn19+1.73734_wp*nn110)

          nn28=TANH(0.18233_wp+0.30421_wp*nn11-2.08994_wp*nn12-0.55899_wp*nn13 &
            &    +0.16582_wp*nn14+0.15934_wp*nn15+2.19585_wp*nn16+0.87302_wp*nn17 &
            &    -0.6635_wp*nn18+0.89907_wp*nn19-1.65624_wp*nn110)

          nn29=TANH(0.85865_wp-0.40122_wp*nn11+0.75878_wp*nn12+1.15085_wp*nn13 &
            &    +0.8702_wp*nn14-1.00696_wp*nn15-0.15103_wp*nn16-1.21695_wp*nn17 &
            &    +1.0322_wp*nn18-0.51579_wp*nn19-0.0329_wp*nn110)

          nn210=TANH(0.02923_wp+1.73947_wp*nn11+2.5155_wp*nn12+0.50039_wp*nn13 &
            &    -0.41072_wp*nn14+0.38135_wp*nn15-1.11433_wp*nn16-0.04792_wp*nn17 &
            &    +2.17579_wp*nn18+1.15646_wp*nn19-15.12584_wp*nn110)


          f_r_rh(jc) = (-0.02991_wp+1.70659_wp*nn21+0.28384_wp*nn22+0.29133_wp*nn23 &
            &      -0.04257_wp*nn24-2.30272_wp*nn25-0.27945_wp*nn26+1.00396_wp*nn27 &
            &      -0.6744_wp*nn28-1.12556_wp*nn29+1.05196_wp*nn210)


          ! the neural network can produce values below zero. Set them to zero
          IF (f_r_rh(jc) < 0._wp) THEN
            f_r_rh(jc) = 0._wp
          ENDIF

          f_r_t(jc) = 1._wp !set to 1 to make dummy factor of it 

        ENDIF ! p_llsm_atm_c
      ENDDO ! jc
      !$ACC END PARALLEL

    CASE ('pollpoac')  ! Formulas for grasses

!NEC$ ivdep
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = istart, iend  !<Loop over all horizontal grid points
        IF( p_llsm_atm_c(jc) ) THEN
          !---------------------------------------------------------------------!
          ! Preliminary calculation: rh                                         !
          !---------------------------------------------------------------------!
          rh(jc) = p_rh_2m(jc) * 0.01_wp

          !...........................................................
          ! Meteorological influences:
          !...........
          f_r_rh(jc) =   2._wp *  (1._wp / (1._wp + EXP(-rh(jc) * 30._wp + 9._wp))) *   &
            &         EXP(-((2._wp * rh(jc)-1._wp)**2)) *                               &
            &         (0.5_wp * (-1._wp) * EXP(-((2._wp * rh(jc)-1._wp)**2)) + 1._wp) * &
            &         (1._wp / (1._wp + EXP(45._wp * rh(jc) - 40._wp)))

          f_r_t(jc) = 1._wp / (1._wp + EXP(-p_temp(jc) * 0.272_wp + 78._wp))

        ENDIF ! p_llsm_atm_c
      ENDDO ! jc
      !$ACC END PARALLEL

    CASE ('pollalnu')  ! Formulas for alder

!NEC$ ivdep
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO jc = istart, iend  !<Loop over all horizontal grid points
        IF( p_llsm_atm_c(jc) ) THEN
          !---------------------------------------------------------------------!
          ! Preliminary calculation: rh                                         !
          !---------------------------------------------------------------------!
          rh(jc) = p_rh_2m(jc) * 0.01_wp

          !...........................................................
          ! Meteorological influences:
          !...........
          f_r_rh(jc) = 1._wp / (1._wp + EXP(21.00_wp * rh(jc) -  15._wp))

          f_r_t(jc) = 1._wp * (1._wp / (1._wp + EXP(2.2_wp*(-0.27_wp * p_temp(jc) + 75.2_wp))))

        ENDIF ! p_llsm_atm_c
      ENDDO ! jc
      !$ACC END PARALLEL

  CASE ('pollcory')  ! Formulas for cory

!NEC$ ivdep
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO jc = istart, iend  !<Loop over all horizontal grid points
      IF( p_llsm_atm_c(jc) ) THEN

        !---------------------------------------------------------------------!
        ! Preliminary calculation: rh                                         !
        !---------------------------------------------------------------------!
        rh(jc) = p_rh_2m(jc) * 0.01_wp

        !...........................................................
        ! Meteorological influences:
        !...........
        ! Formulas for hazel from Christina Endler (DWD)
        f_r_rh(jc) = 1._wp / (1._wp + EXP(20.00_wp * rh(jc) -  16._wp))
        f_r_t(jc) = 1._wp *                  & 
                (1._wp / (1._wp + EXP(2.4_wp*(-0.27_wp * p_temp(jc) +  74.9_wp)))) 


      ENDIF ! p_llsm_atm_c
    ENDDO ! jc
    !$ACC END PARALLEL

  CASE DEFAULT
    CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
      &          'No explicit emission treatment for selected mode_name '//TRIM(mode_name)//'.')
  END SELECT

!NEC$ ivdep
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
  !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(max_emiss_day, phi_biol, phi_met, phi_plant)
  DO jc = istart, iend  !<Loop over all horizontal grid points

    IF( p_llsm_atm_c(jc) ) THEN

      !---------------------------------------------------------------------!
      ! The Reservoir                                                       !
      !---------------------------------------------------------------------!

      !...........................................................
      ! Biological influence:
      !...........
      ! Determines the maximum amount of pollen that can be produced per day.
      ! The switch 'phi_biol' is 1 if this amount hasn't been reached, and
      ! turns to 0 as soon as the maximum possible amount has been released
      ! from the flowers into the reservoir on a given day.
      ! The value of this daily maximum depends on the time of year. Thus, the
      ! calculation of 'phi_biol' requires a description of the pollen season
      ! which is called 'SDES'. 'SDES' is between 0 and 1, the
      ! height of the maximum is determined by 'no_max_day'.
      !...........
      max_emiss_day = pollen_data%no_max_day(jc,jb) * pollen_data%fr_cov(jc,jb)  &
        &              * pollen_data%f_q_alt(jc,jb)                       &
        &              * pollen_data%f_q_seas(jc,jb)

      IF (pollen_data%res_new_sum(jc,jb) < max_emiss_day) THEN
        phi_biol = 1._wp
      ELSE
        phi_biol = 0._wp
      ENDIF

      !...........................................................
      ! Amount of pollen available:
      !...........

      phi_plant = pollen_data%f_q_seas(jc,jb) * pollen_data%no_max_timestep(jc,jb)  &
        &       * pollen_data%fr_cov(jc,jb) * pollen_data%f_q_alt(jc,jb)

      !...........................................................
      ! Meteorological influences:
      !...........

      phi_met = f_r_t(jc) * f_r_rh(jc)

      !...........................................................
      ! Release of the pollen from the flowers,
      ! filling of the reservoir:
      !...........
      res_new(jc) = phi_biol * phi_met * phi_plant

    ENDIF ! p_llsm_atm_c
  ENDDO ! jc

!NEC$ ivdep
  !$ACC LOOP GANG(STATIC: 1) VECTOR &
  !$ACC   PRIVATE(c_buff, c_max, c_quell, delta_r_precip, f_e_tke, precip, psi_emis, psi_precip, psi_wet, time_sec)
  DO jc = istart, iend  !<Loop over all horizontal grid points

    IF( p_llsm_atm_c(jc) ) THEN

      !---------------------------------------------------------------------!
      ! Preliminary calculations: precip                                    !
      !---------------------------------------------------------------------!
      !rain rate at surface
      precip = p_rain_gsp_rate(jc) + p_rain_con_rate(jc)

      !...........................................................
      ! Emptying of the reservoir:
      !...........
      ! loss of pollen from the reservoir due to precipitation
      IF (precip >= 0.0005_wp) THEN
        psi_precip = 0._wp
      ELSEIF (precip < 0._wp) THEN
        psi_precip = 1._wp
      ELSE
        psi_precip = -2000._wp * precip + 1._wp
      END IF

      !...........................................................
      ! New content of the reservoir:
      !...........

      pollen_data%res_old(jc,jb) = (pollen_data%psi_random * pollen_data%res_old(jc,jb) + res_new(jc)) &
        &                          * psi_precip

      ! no negative content
      IF (pollen_data%res_old(jc,jb)<0._wp) THEN
        pollen_data%res_old(jc,jb) = 0._wp
      ENDIF

      !...........................................................
      ! Reduction of emission under moist conditions:
      !...........

      IF (rh(jc) > 0.95_wp) THEN
        f_e_rh(jc) = 0._wp
      ELSE IF (rh(jc) <= 0.95_wp .AND. rh(jc) > 0.9_wp) THEN
        f_e_rh(jc) = 0.5_wp
      ELSE IF (rh(jc) <= 0.9_wp) THEN
        f_e_rh(jc) = 1._wp
      ENDIF

      !...........................................................
      ! Suppression of emission after precipitation:
      ! JS: psi_wet is an additional term to EMPOL 1.0
      !...........
      pollen_data%r_precip(jc,jb) = pollen_data%r_precip(jc,jb) + precip * p_dtime
      ! ISNAN-Check
      IF (.NOT. (pollen_data%r_precip(jc,jb)==pollen_data%r_precip(jc,jb))) THEN
        pollen_data%r_precip(jc,jb) = 0._wp
      ENDIF

      IF (pollen_data%r_precip(jc,jb) >= pollen_data%xi_r_precip) THEN
        pollen_data%r_precip(jc,jb) = pollen_data%xi_r_precip
      ENDIF

      IF (precip < 0.0001_wp) THEN
        delta_r_precip = (1._wp - rh(jc)) * pollen_data%frac_xi_evap
        pollen_data%r_precip(jc,jb) = pollen_data%r_precip(jc,jb) - delta_r_precip
        IF (pollen_data%r_precip(jc,jb) < 0._wp) THEN
          pollen_data%r_precip(jc,jb) = 0._wp
        ENDIF
      ENDIF

      IF (pollen_data%r_precip(jc,jb) > EPSILON(0._wp)) THEN
        psi_wet = 0._wp
      ELSE
        psi_wet = 1._wp
      ENDIF

      !...........................................................
      ! Influence of TKE:
      !...........

      f_e_tke = 1._wp / (1._wp + EXP(-2.1_wp * p_tke(jc) + 4._wp)) - 0.017_wp
      IF (f_e_tke < 0._wp) THEN
        f_e_tke = 0._wp
      ENDIF

      !...........................................................
      ! Calculation of the emission flux:
      !...........
      c_max    = pollen_data%res_old(jc,jb) / p_dz_nlev(jc)                  !(1/m3)
      c_quell  = c_max * f_e_tke * f_e_rh(jc) * psi_wet                      !(1/m3)
      psi_emis = pollen_data%res_old(jc,jb) * f_e_tke * f_e_rh(jc) * psi_wet !(1/m2)
      pollen_data%fe_plant(jc,jb) = psi_emis / p_dtime                       !(1/s m2) emission flux

      ! pick out emitted pollen from reservoir
      pollen_data%res_old(jc,jb) = pollen_data%res_old(jc,jb) - psi_emis


      ! previous concentration + current concentration
      c_buff = (p_tracer(jc,ipoll_tr) * p_rho(jc)) + c_quell    !(1/m3)

      !...........................................................
      ! Update pollen tracer concentration with current emission
      !...........
      p_tracer(jc,ipoll_tr) = c_buff / p_rho(jc)        ! (1/kg)

      !...........................................................
      ! Sum up ruptured pollen per day
      ! OR: If it's a new day --> clean up
      !...........
      time_sec = hour * 3600._wp  &
        &      + minute * 60._wp  &
        &      + second

      IF( time_sec - p_dtime >= 0) THEN
        pollen_data%res_new_sum(jc,jb) = pollen_data%res_new_sum(jc,jb) + res_new(jc)
      ELSE
        pollen_data%res_new_sum(jc,jb) = res_new(jc)  !first timestep for the day
      ENDIF

    ENDIF ! p_llsm_atm_c
  ENDDO ! jc
  !$ACC END PARALLEL

  !$ACC WAIT
  !$ACC END DATA

END SUBROUTINE art_emiss_pollen
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_pollen_get_nstns( p_pollen_prop,  mode_name, n_stns )
!<
! SUBROUTINE art_pollen_get_nstns
! retrieves the n_stns information dependent on mode_name for the use in the
! season-length calculation
!
! Part of Module: mo_art_emission_pollen
! Author: Sven Werchner, KIT
! Initial Release: 2021-03-16
! Modifications:
!<

  TYPE(t_art_pollen_properties),INTENT(INOUT) :: &
    & p_pollen_prop               !< All data related to specific pollen type emissions

  CHARACTER(LEN=IART_VARNAMELEN),INTENT(IN) :: &
    &  mode_name                  !< Mode name,
  INTEGER,INTENT(out)       :: &
    &  n_stns


! Local variables

  TYPE(t_art_pollen_table),POINTER :: &
    &  pollen_data                 !< All data related to specific pollen type emissions


  INTEGER                    :: &
    &  ipoll,                   & !< Index of pollen in pollen table
    &  ierr                       !< Error return value


  CALL p_pollen_prop%dict_pollen%get(mode_name, ipoll, ierr)
  IF(ierr /= SUCCESS) CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
                          &          'ipoll not found in pollen table dictionary. &
                          &           mode_name:'//TRIM(mode_name) )

  pollen_data => p_pollen_prop%pollen_type(ipoll)
  IF ( TRIM(mode_name) /= 'pollpoac' .AND.  &
    &  TRIM(mode_name) /= 'pollambr') THEN
    n_stns = pollen_data%pol_atab%all_stns%n_stns
  ELSE
    n_stns = 1
  END IF

END SUBROUTINE art_pollen_get_nstns
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_tsum( current_date, mode_name,                &
  &                          p_pollen_prop, t_2m, jb, istart, iend )
!<
! SUBROUTINE art_prepare_emission_pollen
! This module calculates the phenological state of the plants (Variable f_q_seas)
! that
! is used in the pollen emission calculation (module pol_emissions). f_q_seas is
! zero
! before and after the pollen season. During the pollen season it ranges from
! zero and  (almost) one. The higher f_q_seas the more plants are flowering.
!
! The current implementation includes birch, alder and grasses. For birch and
! alder, a temperature
! sum model for the start and the end of the pollen season is provided. This
! model is optimized
! for Swiss pollen data. For grasses, the implemented approach includes a
! temperature sum model
! for the start of the pollen season. The end of the pollen season is calculated
! via the
! climatological length of the grass pollen season. Most of the subroutines can
! handle further species some have to be adapted though. For one further species
! (Ambrosia or any other) the structures are already implemented.
!
! Required input per species: cumulative temperature threshold fields for the
! start
! (birch, alder and grasses) and the end (birch and alder) of the pollen season.
! If the length
! of the pollen season is calculated (birch and alder) a set of stations with
! daily climatological
! t2m temperature (12 UTC) is needed in ATAB format. For grasses, the
! climatological
! length of the pollen season has to be provided as well. The calculation
! of the start of the pollen season requires daily 14h UTC t2m temperature fields
! in ATAB format.
! If used in the "operational" mode (calc_sdes_t2m = .FALSE.), the variables
! saisn and ctsum have to be provided instead of the (external) daily t2m
! temperature fields.
!
! Based on Pauling et al. 2014 - Toward optimized temperature sum parameterizations
!                      for forecasting the start of the pollen season, Aerobiologia
!
! Part of Module: mo_art_emission_pollen
! Author: Andreas Pauling, MeteoSwiss
! Initial Release: 2015-01-20
! Modifications:
! 2017-12-12: Jonas Straub, KIT
! - Migration from COSMO-ART to ICON-ART
! 2021-03-16: Sven Werchner, KIT
! - moving determination of n_stns out to separate SR
!<

  TYPE(t_art_pollen_properties),INTENT(INOUT) :: &
    & p_pollen_prop               !< All data related to specific pollen type emissions

  TYPE(datetime),INTENT(IN)  :: &
    &  current_date               !< Date & time in mtime format

  INTEGER,INTENT(IN)         :: &
    &  jb,                      & !< block loop index
    &  istart, iend               !< start and end of nproma

  REAL(wp),INTENT(IN)        :: &
    &  t_2m(:)                    !< temperature in 2m  (K)

  CHARACTER(LEN=IART_VARNAMELEN),INTENT(IN) :: &
    &  mode_name                  !< Mode name,


! Local variables

  TYPE(t_art_pollen_table),POINTER :: &
    &  pollen_data                 !< All data related to specific pollen type emissions


  INTEGER                    :: &
    &  dayinyear,               & !< day in year, with 1st January (= first day)
    &  ipoll,                   & !< Index of pollen in pollen table
    &  ierr,                    & !< Error return value
    &  doy_dec1                   !< days since 1st December (= first day)


  CALL p_pollen_prop%dict_pollen%get(mode_name, ipoll, ierr)
    IF(ierr /= SUCCESS) CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
                            &          'ipoll not found in pollen table dictionary. &
                            &           mode_name:'//TRIM(mode_name) )

  pollen_data => p_pollen_prop%pollen_type(ipoll)

  !-----------------------------------------------------------------------------------------
  !--   Run calculations once a day: at 12 UTC ---------------------------------------------
  !-----------------------------------------------------------------------------------------
    IF (current_date%date%month == 12) THEN
      doy_dec1 = current_date%date%day
    ELSE
      dayinyear = getDayOfYearFromDateTime(current_date, ierr)  ! similar to nzjulianday in COSMO-ART
      doy_dec1 = dayinyear + 31
    ENDIF

    !-----------------------------------------------------------------------------------------
    !--   Start Routine ----------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------	

    IF (TRIM(mode_name) /= 'pollambr' .AND. doy_dec1 > pollen_data%jul_days_excl) THEN
      CALL art_calc_t2m_sum_oper(istart, iend,                 &
        &                        t_2m(:),                      &
        &                        pollen_data%ctsum(:,jb),      &
        &                        pollen_data%saisn(:,jb),      &
        &                        pollen_data%saisa(:,jb),      &
        &                        pollen_data%saisl(:,jb),      &
        &                        pollen_data%tthrs_red(:,jb),  &
        &                        pollen_data%tthrs(:,jb),      &
        &                        pollen_data%tthre(:,jb),      &
        &                        pollen_data%t_base,           &
        &                        mode_name,                    &
        &                        doy_dec1 - pollen_data%jul_days_excl)
    ENDIF

! JF:   END IF

END SUBROUTINE art_prepare_tsum

SUBROUTINE art_prepare_saisl ( p_pollen_prop, current_date, mode_name, saisl_stns )

  TYPE(t_art_pollen_properties),INTENT(INOUT) :: &
    & p_pollen_prop               !< All data related to specific pollen type emissions

  TYPE(datetime),INTENT(IN)  :: &
    &  current_date               !< Date & time in mtime format

  REAL(wp),INTENT(out)        :: &
    &  saisl_stns(:)

  CHARACTER(LEN=IART_VARNAMELEN),INTENT(IN) :: &
    &  mode_name                  !< Mode name,

! Local variables

  TYPE(t_art_pollen_table),POINTER :: &
    &  pollen_data                 !< All data related to specific pollen type emissions


  INTEGER                    :: &
    &  dayinyear,               & !< day in year, with 1st January (= first day)
    &  ipoll,                   & !< Index of pollen in pollen table
    &  ierr,                    & !< Error return value
    &  doy_dec1                   !< days since 1st December (= first day)


  CALL p_pollen_prop%dict_pollen%get(mode_name, ipoll, ierr)
  IF(ierr /= SUCCESS) CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
                          &          'ipoll not found in pollen table dictionary. &
                          &           mode_name:'//TRIM(mode_name) )

  pollen_data => p_pollen_prop%pollen_type(ipoll)

  IF (current_date%date%month == 12) THEN
    doy_dec1 = current_date%date%day
  ELSE
    dayinyear = getDayOfYearFromDateTime(current_date, ierr)  ! similar to nzjulianday in COSMO-ART
    doy_dec1 = dayinyear + 31
  ENDIF

  CALL art_calc_saisl(dayinyear,                            &
    &                 pollen_data%t_base,                   &
    &                 pollen_data%jul_days_excl,            &
    &                 saisl_stns(:),                        &
    &                 pollen_data%saisa(:,:),               &
    &                 pollen_data%ctsum(:,:),               &
    &                 pollen_data%tthrs_red(:,:),           &
    &                 pollen_data%tthre(:,:),               &
    &                 doy_dec1 - pollen_data%jul_days_excl, &
    &                 pollen_data%pol_atab)

END SUBROUTINE art_prepare_saisl


SUBROUTINE art_prepare_sdes (p_pollen_prop, p_patch, &
    &                        jb, istart, iend, mode_name, saisl_stns)

  TYPE(t_art_pollen_properties),INTENT(INOUT) :: &
    & p_pollen_prop               !< All data related to specific pollen type emissions

  TYPE(t_patch), TARGET ::  &
    &  p_patch                 !< Patch on which computation is performed

  INTEGER,INTENT(IN)         :: &
    &  jb,                      & !< block loop index
    &  istart, iend               !< start and end of nproma

  REAL(wp),INTENT(IN)        :: &
    &  saisl_stns(:)

  CHARACTER(LEN=IART_VARNAMELEN),INTENT(IN) :: &
    &  mode_name                  !< Mode name,

! Local variables

  TYPE(t_art_pollen_table),POINTER   :: &
    &  pollen_data                 !< All data related to specific pollen type emissions

  INTEGER                    :: &
    &  ipoll,                   & !< Index of pollen in pollen table
    &  ierr                       !< Error return value


  CALL p_pollen_prop%dict_pollen%get(mode_name, ipoll, ierr)
  IF(ierr /= SUCCESS) CALL finish (TRIM(thisroutine)//':art_emiss_pollen', &
                            &          'ipoll not found in pollen table dictionary. &
                            &           mode_name:'//TRIM(mode_name) )

  pollen_data => p_pollen_prop%pollen_type(ipoll)

  CALL art_calc_sdes(istart, iend, mode_name, jb,  &
    &            pollen_data%f_q_seas(:,jb),       &
    &            pollen_data%fr_cov(:,jb),         &
    &            pollen_data%saisn(:,jb),          &
    &            pollen_data%saisl(:,jb),          &
    &            saisl_stns,                       &
    &            p_patch,                          &
    &            pollen_data%pol_atab%header_info)

END SUBROUTINE art_prepare_sdes
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_calc_t2m_sum_oper( istart, iend, t_2m, ctsum, &
  &                               saisn, saisa, saisl,       &
  &                               tthrs_red, tthrs, tthre,   &
  &                               t_base, mode_name,         &
  &                               days_ctsum )
!<
! SUBROUTINE art_calc_t2m_sum_oper
! Building the 2m-temp sum by taking value from the day before
!
! Based on Pauling et al. 2014 - Toward optimized temperature sum parameterizations
!                      for forecasting the start of the pollen season, Aerobiologia
!
! Part of Module: mo_art_emission_pollen
! Author: Andreas Pauling, MeteoSwiss
! Initial Release: 2015-01-20
! Modifications:
! 2018-01-10: Jonas Straub, KIT
! - Migration from COSMO-ART to ICON-ART
!<

  INTEGER, INTENT(IN)       :: &
    &  istart, iend,           & !< loop indices
    &  days_ctsum                !< days since start of temperature cummulation

  REAL(wp), INTENT(INOUT)   :: &
    &  ctsum(:),               &  !< cumulated temperature sum (adegreeC)
    &  tthrs_red(:),           &  !< reduction of the threshold to account for the fact that
                                  !  there are already pollen emitted before Pollen>30		
    &  saisn(:),               &  !< number of days since the start of the pollen
                                  !  season if present day is during the season; zero outside the season (1)
    &  saisa(:)                   !< as saisn, but contains length of season if _a_fter the
                                  !  season and not zeros (contains _a_ll days of the season) (1)

  REAL(wp),INTENT(IN)       :: &
    &  t_2m(:),                &   !< temperature in 2 m height (K)
    &  tthrs(:),               &   !< cumulated temperature sum threshold for the start of the pollen season
    &  tthre(:),               &   !< cumulated temperature sum threshold for the end   of the pollen season
    &  saisl(:)                    !< length of pollen seasons

  REAL(wp),INTENT(IN)       :: &
    &  t_base                      !< optimized start of flowering model (adegreeC resp. K)	

  CHARACTER(LEN=IART_VARNAMELEN), INTENT(IN) :: &
    &  mode_name                   !< Mode name,

! local variables
  REAL(wp)                  :: &
    &  t2m_cels,               & !< t2m temperature in Celsius (adegreeC)
    &  tthrs_red_offset          !< offset value for temperature threshold reduction

  INTEGER                   :: &
    &  jc                        !< loop index for nproma


  !-----------------------------------------------------------------------------------------
  !--   START  routine   -------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SELECT CASE (TRIM(mode_name))
      ! The temperature sums for the start (end) of the season were derived from pollen measurement data. 
      ! To obtain statistically stable results for the temperature sums we defined the beginning of the 
      ! pollen season as the day when pollen concentration > XX Pollen /m**3, while the value of XX
      ! depends on the pollen species (e.g. 20 for grass pollen).
      ! To account for the fact that pollen is already in the atmosphere, we let the season start earlier 
      ! in the model by reducing the temperature sum by a certain value (tthrs_red_offset). 
      ! Since this threshold reduction applied to the temperature sum, we have to include a weighting as
      ! is is done in the calculation of the temperature sum (weighting by the time of the year).
    CASE('pollbetu')
      ! reduction of the threshold to account for the fact that there are already pollen
      ! emitted before Pollen>30 (984=1day*12deg.*82weight)
      tthrs_red_offset = 984._wp

    CASE('pollpoac')
      ! reduction of the threshold to account for the fact that there are already pollen
      ! pollen emitted before Pollen>=20 (10272=6days*16deg.*107weight)
      tthrs_red_offset = 10272._wp

    CASE('pollalnu')
      ! reduction of the threshold to account for the fact that there are already pollen
      ! emitted before Pollen>=5
      ! reducing tthrs not necessary as Pollen>=5 is the start (and not Pollen>30)
      tthrs_red_offset = 0.0_wp

    CASE('pollcory')
      ! reduction of the threshold to account for the fact that there are already pollen
      ! emitted before Pollen>=5
      ! For hazel, this correction is already included in tthrs! 
      tthrs_red_offset = 0.0_wp

    ! for pollambr tthrs_red_offset is not needed, since temperature sums are not calculated.
    ! the season parametrization (sdes) is based on a climatology.

    CASE DEFAULT
      CALL finish(TRIM(thisroutine)//':art_calc_t2m_sum_oper', &
        &         'T2M sum calculation for '//TRIM(mode_name)//' during runtime not implemented.')
  END SELECT

  IF ( TRIM(mode_name) == 'pollpoac') THEN

!NEC$ ivdep
    DO jc = istart, iend  !<Loop over all horizontal grid points

      ! set t2m_cels to zero if base temperature is not reached
      t2m_cels = t_2m(jc) - 273.15_wp
      IF (t2m_cels > t_base) THEN
        ctsum(jc) = ctsum(jc) + t2m_cels * days_ctsum
      ENDIF

      ! reduction of the threshold to account for the fact that there are already pollen
      tthrs_red(jc) = MAX(tthrs(jc) - tthrs_red_offset, 0.0_wp)

      ! set 2D-Array sais_now to 1 where tthrs_red reached and > 0 (no season), and
      ! saisl not yet reached
      IF (ctsum(jc)     > tthrs_red(jc) .AND. &
        tthrs_red(jc) > 0._wp           .AND. &
        saisa(jc)     < saisl(jc)) THEN
        saisn(jc) = saisn(jc) + 1.0_wp
        saisa(jc) = saisa(jc) + 1.0_wp  ! save saisn that keeps length of season even after the season
      ENDIF

      IF (saisa(jc) > saisl(jc)) saisn(jc) = 0._wp

    END DO !jc

  ELSE

!NEC$ ivdep
    DO jc = istart, iend  !<Loop over all horizontal grid points

      ! set t2m_cels to zero if base temperature is not reached
      t2m_cels = t_2m(jc) - 273.15_wp
      IF (t2m_cels > t_base) THEN
        ctsum(jc) = ctsum(jc) + t2m_cels * days_ctsum
      ENDIF

      ! reduction of the threshold to account for the fact that there are already pollen
      tthrs_red(jc) = MAX(tthrs(jc) - tthrs_red_offset, 0.0_wp)

      ! set 2D-Array sais_now to 1 where tthrs_red reached and > 0 (no season), and
      ! tthre not yet reached
      IF (ctsum(jc)     > tthrs_red(jc) .AND. &
        ctsum(jc)     < tthre(jc)       .AND. &
        t2m_cels      >  t_base         .AND. &
        tthrs_red(jc) > 0._wp) THEN
        saisn(jc) = saisn(jc) + 1.0_wp
        saisa(jc) = saisa(jc) + 1.0_wp   ! save saisn that keeps length of season even after the season
      ENDIF

      IF (ctsum(jc) > tthre(jc)) saisn(jc) = 0._wp

    END DO !jc

  ENDIF

END SUBROUTINE art_calc_t2m_sum_oper
!!!
!!!-------------------------------------------------------------------------
!!!
SUBROUTINE art_calc_saisl(dayinyear,t_base,jul_days_excl,saisl_stns,saisa,ctsum,tthrs_red,tthre,   &
  &                       days_ctsum,pol_atab)
!<
! SUBROUTINE art_calc_saisl
! Calculates the length of the season by using saisn, ctsum and climatologies of temp-measuring
! stations
!
! Based on Pauling et al. 2014 - Toward optimized temperature sum parameterizations
!                      for forecasting the start of the pollen season, Aerobiologia
!
! Part of Module: mo_art_emission_pollen
! Author: Andreas Pauling, MeteoSwiss
! Initial Release: 2015-01-20
! Modifications:
! 2018-01-17: Jonas Straub, KIT
! - Migration from COSMO-ART to ICON-ART
!<


  TYPE(t_art_pol_atab),INTENT(IN)       :: &
    &  pol_atab

  INTEGER,INTENT(IN)         :: &
    &  jul_days_excl,           & !< pollen-specific start day for temperature sum minus 1
    &  dayinyear,               & !< day in year, since 1st January
    &  days_ctsum                 !< days since start of temperature cummulation


  REAL(wp),INTENT(INOUT)     :: &
    &  saisl_stns(:),           &  !< length of pollen seasons at gridpoint

    &  t_base,                  &  !< optimized start of flowering model (adegreeC resp. K)	
    &  tthrs_red(:,:),          &  !< reduction of the threshold to account for the fact that
    &  tthre(:,:),              &  !< cumulated temperature sum threshold for the end of the
                                   !  there are already pollen emitted before Pollen>30		
    &  ctsum(:,:),              &  !< cumulated temperature sum (adegreeC)
    &  saisa(:,:)                  !< as saisn, but contains length of season if _a_fter the
                                   !  season and not zeros (contains _a_ll days of the season) (1)
  REAL(wp),ALLOCATABLE       :: &
    &  array(:,:)                  !<

! local variables
  TYPE(t_art_all_stns)       :: &
    &  all_stns

  INTEGER                    ::   &
    &  tomorrow_end_clim,         & !< =tomorrow till last day of climatological period (e.g.31 May)
    &  day,                       & !< loop index for tomorrow_end_clim
    &  stn,                       & !< loop index for stations
    &  tri_iidx_loc,tri_iblk_loc, & !< indices of station locations
    &  nrow,                      & !< number of rows (days) in input file
    &  n_stns,                    & !< number of stations in input file (=ncol)
    &  doy_dec1                     !< doy_dec1 are the days of the year with start Dec 1

  INTEGER,ALLOCATABLE        ::   &
    &  iarray(:)

  REAL(wp),ALLOCATABLE               :: &
    &  t2m_norm_stns(:,:),          & !< array with all temp means stored in input file
    &  t2m_clim_stns(:,:),          & !< climatology (tomorrow -> 31 May) of stations;
                                      !  extract of t2m_norm_stns(:,:)
    &  t2m_clim_stns_w(:,:),        & !< t2m_clim_stns weighted with doi
    &  t2m_clim_intg_distsum(:),    & !< distance weighted sum
    &  t2m_clim_intg_distw(:),      & !< inverse dist. weighted climate data at current gridpoint
    &  t2m_clim_intg_stns(:,:),     & !< integrated temp over all days at all stations
    &  t2m_clim_csum_off(:),        & !<
    &  t2m_clim_csum_off_stns(:,:), & !<
    &  t2m_clim_csum_stns(:,:),     & !<
    &  sais_clim(:),                & !< logical for all days: during season=1., else=0.
    &  sais_clim_stns(:,:),         & !< logical for all days: during season=1., else=0.
    &  doy_dec1_stns(:,:),          & !<
    &  thr_s_stns(:,:),             & !<
    &  thr_e_stns(:,:),             & !<
    &  ctsum_stns(:),               & !<
    &  saisa_stns(:),               & !<
    &  thr_s_stns_tmp1(:),          & !<
    &  thr_e_stns_tmp1(:)             !<

  !-----------------------------------------------------------------------------------------
  !--   PREPARE external station climatologies for season calculation   --------------------
  !-----------------------------------------------------------------------------------------		

  !-----------------------------------------------------------------------------------------
  !--   ALLOCATIONS of local fields           ----------------------------------------------
  !-----------------------------------------------------------------------------------------

  n_stns = pol_atab%all_stns%n_stns
  nrow = pol_atab%n_days
  ALLOCATE(iarray(nrow))
  iarray = pol_atab%idays

  ! define timespan (in days since Jan 1) used for allocation. tomorrow_end_clim refers to DOI
  !                                                            in atab file
  tomorrow_end_clim     = iarray(size(iarray)) - dayinyear

  ALLOCATE( t2m_clim_intg_distsum (        tomorrow_end_clim) )
  ALLOCATE( t2m_clim_intg_distw   (        tomorrow_end_clim) )
  ALLOCATE( t2m_clim_csum_off     (        tomorrow_end_clim) )
  ALLOCATE( sais_clim             (        tomorrow_end_clim) )

  ALLOCATE( t2m_clim_intg_stns    (n_stns, tomorrow_end_clim) )
  ALLOCATE( t2m_clim_stns         (n_stns, tomorrow_end_clim) )
  ALLOCATE( t2m_clim_stns_w       (n_stns, tomorrow_end_clim) )
  ALLOCATE( t2m_clim_csum_stns    (n_stns, tomorrow_end_clim) )
  ALLOCATE( t2m_clim_csum_off_stns(n_stns, tomorrow_end_clim) )
  ALLOCATE( sais_clim_stns        (n_stns, tomorrow_end_clim) )
  ALLOCATE( doy_dec1_stns         (n_stns, tomorrow_end_clim) )
  ALLOCATE( thr_s_stns            (n_stns, tomorrow_end_clim) )
  ALLOCATE( thr_e_stns            (n_stns, tomorrow_end_clim) )

  ALLOCATE( ctsum_stns            (n_stns                   ) )
  ALLOCATE( saisa_stns            (n_stns                   ) )
  ALLOCATE( thr_s_stns_tmp1       (n_stns                   ) )
  ALLOCATE( thr_e_stns_tmp1       (n_stns                   ) )

  ALLOCATE( array                 (n_stns, nrow             ) )
  ALLOCATE( t2m_norm_stns         (n_stns, nrow             ) )

  ! Initializations
  t2m_clim_intg_distsum(:)     = 0._wp
  t2m_clim_intg_distw(:)       = 0._wp
  t2m_clim_intg_stns(:,:)      = 0._wp
  t2m_clim_csum_off_stns(:,:)  = 0._wp
  sais_clim(:)                 = 0._wp
  doy_dec1_stns(:,:)           = 0._wp
! t2m_norm_stns(:,:)           = REAL(array(:,:),wp)
  t2m_clim_stns(:,:)           = 0._wp
  t2m_clim_stns_w(:,:)         = 0._wp
  t2m_clim_csum_stns(:,:)      = 0._wp
  sais_clim_stns(:,:)          = 0._wp
  thr_s_stns(:,:)              = 0._wp
  thr_e_stns(:,:)              = 0._wp
  ctsum_stns(:)                = 0._wp
  saisa_stns(:)                = 0._wp
  thr_s_stns_tmp1(:)           = 0._wp
  thr_e_stns_tmp1(:)           = 0._wp
  saisl_stns(:)                = 0._wp

  all_stns = pol_atab%all_stns
  array = pol_atab%temp_array
  t2m_norm_stns(:,:)           = REAL(array(:,:),wp)


  DO stn=1,n_stns
    IF(all_stns%p(stn)%ithis_nlocal_pts >0) THEN
      tri_iidx_loc = all_stns%p(stn)%tri_iidx_loc
      tri_iblk_loc = all_stns%p(stn)%tri_iblk_loc
      ! get thresholds for ctsum at the stations from the whole field
      ctsum_stns(stn)      = ctsum(tri_iidx_loc,tri_iblk_loc)
      saisa_stns(stn)      = saisa(tri_iidx_loc,tri_iblk_loc)
      thr_s_stns_tmp1(stn) = tthrs_red(tri_iidx_loc,tri_iblk_loc)
      thr_e_stns_tmp1(stn) = tthre(tri_iidx_loc,tri_iblk_loc)
    ENDIF
  END DO

  thr_s_stns = SPREAD(thr_s_stns_tmp1, 2, tomorrow_end_clim)
  thr_e_stns = SPREAD(thr_e_stns_tmp1, 2, tomorrow_end_clim)


  ! Select t_2m climatology tomorrow to 31 May
  t2m_clim_stns = t2m_norm_stns(:,dayinyear - iarray(1) + 2 : SIZE(iarray))

  ! Set t_2m < t_base(isp) to zero and calculate weighted temperature sum
  WHERE (t2m_clim_stns < t_base) t2m_clim_stns = 0._wp

  ! Create 2D-array conformable to t2m_clim_stns
  ! subset of doy_dec1: +31-jul_days_excl(isp), because weights have to refer to the Dec. 1
  doy_dec1_stns = &
    REAL(SPREAD((/(doy_dec1, doy_dec1=days_ctsum+1,iarray(SIZE(iarray))+31-jul_days_excl)/), &
    &            1, n_stns),wp)

  t2m_clim_stns_w = t2m_clim_stns * doy_dec1_stns

  DO day=1,tomorrow_end_clim
    t2m_clim_csum_stns(:,day) = SUM(t2m_clim_stns_w(:,1:day), DIM=2)
  END DO

  ! add temperature sum until today to correct the offset due to past temperatures
  t2m_clim_csum_off_stns = SPREAD(ctsum_stns, 2, tomorrow_end_clim) &
  &                      + t2m_clim_csum_stns


  ! Set 2D-array sais_clim_stns to 1, where during season
  WHERE (t2m_clim_csum_off_stns > thr_s_stns .AND. &
         t2m_clim_csum_off_stns < thr_e_stns) sais_clim_stns = 1._wp

  ! length of the season = past days during season and coming days until end of  season
  saisl_stns = saisa_stns + SUM(sais_clim_stns, DIM=2)

  DEALLOCATE( iarray                 )

  DEALLOCATE( t2m_clim_intg_distsum  )
  DEALLOCATE( t2m_clim_intg_distw    )
  DEALLOCATE( t2m_clim_csum_off      )
  DEALLOCATE( sais_clim              )

  DEALLOCATE( t2m_clim_intg_stns     )
  DEALLOCATE( t2m_clim_stns          )
  DEALLOCATE( t2m_clim_stns_w        )
  DEALLOCATE( t2m_clim_csum_stns     )
  DEALLOCATE( t2m_clim_csum_off_stns )
  DEALLOCATE( sais_clim_stns         )
  DEALLOCATE( doy_dec1_stns          )
  DEALLOCATE( thr_s_stns             )
  DEALLOCATE( thr_e_stns             )

  DEALLOCATE( ctsum_stns             )
  DEALLOCATE( saisa_stns             )
  DEALLOCATE( thr_s_stns_tmp1        )
  DEALLOCATE( thr_e_stns_tmp1        )

  DEALLOCATE( array                  )
  DEALLOCATE( t2m_norm_stns          )

END SUBROUTINE art_calc_saisl
!!!
!!!-------------------------------------------------------------------------
!!!
SUBROUTINE art_calc_sdes(istart, iend, mode_name, jb, f_q_seas, fr_cov, &
  &                      saisn, saisl, saisl_stns,  &
  &                      p_patch, header_info)
!<
! SUBROUTINE art_calc_sdes
! Calculation of f_q_seas (season parameter in art_emiss_pollen()) using saisl und saisn
!
! Based on Pauling et al. 2014 - Toward optimized temperature sum parameterizations
!                      for forecasting the start of the pollen season, Aerobiologia
!
! Part of Module: mo_art_emission_pollen
! Author: Andreas Pauling, MeteoSwiss
! Initial Release: 2015-01-20
! Modifications:
! 2018-01-17: Jonas Straub, KIT
! - Migration from COSMO-ART to ICON-ART
!<
  REAL(wp),INTENT(INOUT)      :: &
    &  saisl(:),                 &  !< length of pollen seasons
    &  f_q_seas(:)                  !< state of pollen season; parameter for emission routine
  REAL(wp),INTENT(IN)         :: &
    &  fr_cov(:),                &  !< fraction of plant coverage in grid box
    &  saisn(:),                 &  !< number of days since the start of the pollen
                                    !  season if present day is during the season; zero outside
                                    !  the season (1)
    &  saisl_stns(:)

  CHARACTER(LEN=IART_VARNAMELEN),INTENT(IN) :: &
    &  mode_name                   !< Mode name
  TYPE(t_patch), INTENT(IN) :: &
    &  p_patch
  TYPE(t_art_pollen_atab_header_info),INTENT(IN) :: &
    &  header_info           !< extract header information
  INTEGER,INTENT(IN)        :: &
    &  istart, iend,           & !< loop indices
    &  jb

 ! local variables
  INTEGER           :: &
    &  jc, n_stns
  REAL(wp)          :: &
    &  tune_flower, scale, lon, lat
  REAL(wp), ALLOCATABLE :: &
    &  diff_lon(:),           &
    &  diff_lat(:),           &
    &  dist(:),               &
    &  stns_lon_deg(:),       &
    &  stns_lat_deg(:),       &
    &  stns_lon_rad(:),       &
    &  stns_lat_rad(:)

  IF (TRIM(mode_name) /= 'pollpoac') THEN

    n_stns = header_info%data_cols
    ALLOCATE(diff_lon(n_stns),diff_lat(n_stns),dist(n_stns))
    ALLOCATE(stns_lon_deg(n_stns),stns_lat_deg(n_stns))
    ALLOCATE(stns_lon_rad(n_stns),stns_lat_rad(n_stns))

    dist(:)                      = 0._wp
    diff_lon(:)                  = 0._wp
    diff_lat(:)                  = 0._wp
    stns_lat_rad(:)              = 0._wp
    stns_lon_rad(:)              = 0._wp
    stns_lat_deg(:)              = 0._wp
    stns_lon_deg(:)              = 0._wp

    ! Initialize Local fields with input data
    ! in atab: lon values are -180>lon>180, maybe icon: 0>lon>360 ?
    stns_lon_deg(:) = REAL(header_info%longitude,wp) + 360 !* deg2rad
    stns_lon_deg(:) = MOD(stns_lon_deg,360._wp) !* deg2rad
    stns_lat_deg(:) = REAL(header_info%latitude,wp)  !* deg2rad

    stns_lon_rad = stns_lon_deg * deg2rad
    stns_lat_rad = stns_lat_deg * deg2rad

    DO jc = istart, iend  !<Loop over all horizontal grid points

      IF (fr_cov(jc) <= 0._wp) CYCLE

      ! coordinates of current model gridpoint
      lon = p_patch%cells%center(jc,jb)%lon
      lat = p_patch%cells%center(jc,jb)%lat

      ! calculate distance between all stations and current gridpoint
      diff_lon = (lon - stns_lon_rad) * COS(lon)
      diff_lat = lat  - stns_lat_rad
      dist = SQRT(diff_lon * diff_lon + diff_lat * diff_lat)

      ! inverse distance weighting
      saisl(jc) =  SUM((1._wp/dist) * saisl_stns(:)) / SUM(1._wp/dist)

    END DO !jc

  ENDIF

  ! Calculation of f_q_seas using saisl und saisn
  SELECT CASE (TRIM(mode_name))

  CASE('pollbetu')

!NEC$ ivdep
    DO jc = istart, iend  !<Loop over all horizontal grid points

      IF (fr_cov(jc) <= 0._wp) CYCLE

      f_q_seas(jc) = 0._wp ! set f_q_seas to zero before updating

      ! Calculation of f_q_seas using saisl und saisn
      ! These curves have been updated by Simon Adamov (MeteoSwiss).
      ! The curves were fit to the pollen measurements in Switzerland from 2000-2020 using
      ! Maximum Likelihood and Gradient Descent
      IF (saisn(jc) > 0._wp) THEN ! Ensures the sdes curve is zero outside season
         IF (saisl(jc) > 0._wp) THEN ! ensures that there is no division by zero
            f_q_seas(jc) = &
               ! (maximum / (1 + exp(-slope1 * (x - midpoint1)))) *
               ! (maximum / (1 + exp( slope2 * (x - midpoint2))))
               (1.0_wp/(1._wp + EXP((-0.713_wp*32.266_wp/saisl(jc))* &
                                    (saisn(jc) - 3.626_wp*saisl(jc)/32.266_wp))))* &
               (1.0_wp/(1._wp + EXP((0.412_wp*32.266_wp/saisl(jc))* &
                                    (saisn(jc) - 22.856_wp*saisl(jc)/32.266_wp))))
         END IF
      END IF

      IF (f_q_seas(jc) < 0._wp) f_q_seas(jc) = 0._wp ! for rounding reasons

    END DO !jc

  CASE('pollpoac')

!NEC$ ivdep
      DO jc = istart, iend  !<Loop over all horizontal grid points

        IF (fr_cov(jc) <= 0._wp) CYCLE

        f_q_seas(jc) = 0._wp ! set f_q_seas to zero before updating

        ! Calculation of f_q_seas using saisl und saisn
        ! These curves have been updated by Simon Adamov (MeteoSwiss).
        ! The curves were fit to the pollen measurements in Switzerland from 2000-2020 using
        ! Maximum Likelihood and Gradient Descent
        IF (saisn(jc) > 0._wp) THEN ! Ensures the sdes curve is zero outside season
           IF (saisl(jc) > 0._wp) THEN ! ensures that there is no division by zero
              f_q_seas(jc) = &
                 ! (maximum / (1 + exp(-slope1 * (x - midpoint1)))) *
                 ! (maximum / (1 + exp( slope2 * (x - midpoint2))))
                 (1.01_wp/(1._wp + EXP((-0.275_wp*86.538_wp/saisl(jc))* &
                                       (saisn(jc) - 16.878_wp*saisl(jc)/86.538_wp))))* &
                 (1.01_wp/(1._wp + EXP((0.156_wp*86.538_wp/saisl(jc))* &
                                     (saisn(jc) - 60.440_wp*saisl(jc)/86.538_wp))))
           END IF
        END IF

        IF (f_q_seas(jc) < 0._wp) f_q_seas(jc) = 0._wp ! for rounding reasons

      END DO !jc

  CASE('pollalnu')

!NEC$ ivdep
      DO jc = istart, iend  !<Loop over all horizontal grid points

        IF (fr_cov(jc) <= 0._wp) CYCLE

        f_q_seas(jc) = 0._wp ! set f_q_seas to zero before updating

        ! Calculation of f_q_seas using saisl und saisn
        ! These curves have been updated by Simon Adamov (MeteoSwiss).
        ! The curves were fit to the pollen measurements in Switzerland from 2000-2020 using
        ! Maximum Likelihood and Gradient Descent
        IF (saisn(jc) > 0._wp) THEN ! Ensures the sdes curve is zero outside season
           IF (saisl(jc) > 0._wp) THEN ! ensures that there is no division by zero
              f_q_seas(jc) = &
                 ! (maximum / (1 + exp(-slope1 * (x - midpoint1)))) *
                 ! (maximum / (1 + exp( slope2 * (x - midpoint2))))
                 (1.62_wp/(1._wp + EXP((-0.272_wp*42.788_wp/saisl(jc))* &
                                       (saisn(jc) - 19.450_wp*saisl(jc)/42.788_wp))))* &
                 (1.62_wp/(1._wp + EXP((0.213_wp*42.788_wp/saisl(jc))* &
                                       (saisn(jc) - 23.435_wp*saisl(jc)/42.788_wp))))
           END IF
        END IF

        IF (f_q_seas(jc) < 0._wp) f_q_seas(jc) = 0._wp ! for rounding reasons

    END DO !jc

  CASE('pollcory')

!NEC$ ivdep
    DO jc = istart, iend  !<Loop over all horizontal grid points

      IF (fr_cov(jc) <= 0._wp) CYCLE

        f_q_seas(jc) = 0._wp ! set f_q_seas to zero before updating

        ! Calculation of f_q_seas using saisl und saisn
        ! These curves have been updated by Simon Adamov (MeteoSwiss).
        ! The curves were fit to the pollen measurements in Switzerland from 2000-2020 using
        ! Maximum Likelihood and Gradient Descent

        IF (saisn(jc) > 0._wp) THEN ! Ensures the sdes curve is zero outside season
           ! ensures that there is no division by zero
           IF (saisl(jc) > 0._wp) THEN
              f_q_seas(jc)      = &
                 ! (maximum / (1 + exp(-slope1 * (x - midpoint1)))) *
                 ! (maximum / (1 + exp( slope2 * (x - midpoint2))))
                 (1.31_wp/(1._wp + EXP((-0.477_wp*36.695_wp/saisl(jc))* &
                                       (saisn(jc) - 12.500_wp*saisl(jc)/36.695_wp))))* &
                 (1.31_wp/(1._wp + EXP((0.222_wp*36.695_wp/saisl(jc))* &
                                       (saisn(jc) - 19.711_wp*saisl(jc)/36.695_wp))))
           END IF
        END IF

        IF (f_q_seas(jc) < 0._wp) f_q_seas(jc) = 0._wp ! for rounding reasons

      END DO !jc

    CASE DEFAULT
      CALL finish('mo_art_emission_pollen:art_calc_sdes', &
        &         'ART: Pollen emissions for '//TRIM(mode_name)//' not implemented')
  END SELECT

  !$ACC UPDATE DEVICE(f_q_seas) ASYNC(1)

  IF (TRIM(mode_name) /= 'pollpoac') THEN
    DEALLOCATE( dist         )
    DEALLOCATE( diff_lon     )
    DEALLOCATE( diff_lat     )
    DEALLOCATE( stns_lon_rad )
    DEALLOCATE( stns_lat_rad )
    DEALLOCATE( stns_lon_deg )
    DEALLOCATE( stns_lat_deg )
  ENDIF

END SUBROUTINE art_calc_sdes

END MODULE mo_art_emission_pollen
