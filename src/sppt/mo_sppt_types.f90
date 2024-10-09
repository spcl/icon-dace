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

! Type definition for SPPT (Stochastic Perturbation of Physics Tendencies)

MODULE mo_sppt_types

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_sppt


  ! SPPT state vector - Note there are no local variables except for subroutine locals.
  TYPE t_sppt


    ! a) fiels used to save the current state of prognostic variables

    REAL(wp), POINTER         &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
    , CONTIGUOUS              &
#endif
    & ::                      &
    & temp_now(:,:,:),        &        ! current temperature  (nproma,nlev,nblks_c)    [K]

    & qv_now(:,:,:),          &        ! current value of tracer - water vapor         [kg/kg]
    & qi_now(:,:,:),          &        ! current value of tracer - cloud ice           [kg/kg]
    & qr_now(:,:,:),          &        ! current value of tracer - rain                [kg/kg]
    & qs_now(:,:,:),          &        ! current value of tracer - snow                [kg/kg]
    & qc_now(:,:,:),          &        ! current value of tracer - cloud water         [kg/kg]
    & qg_now(:,:,:),          &        ! current value of tracer - graupel             [kg/kg]

    ! b) fields use for random number patterns
    & rn_3d(:,:,:),           &        ! final 3D array of random numbers              [-]
    & rn_2d_now(:,:),         &        ! utility array for random numbers              [-]
    & rn_2d_new(:,:),         &        ! utility array for random numbers              [-]

    ! c) fields used for calculations of tendencies 
    & ddt_temp_fast(:,:,:), &          ! fast physics tendency for air temperature     [K/s]
    & ddt_u_fast(:,:,:),    &          ! fast physics tendency for wind component u    [m/s^2]
    & ddt_v_fast(:,:,:),    &          ! fast physics tendency for wind component v    [m/s^2]

    & ddt_qv_fast(:,:,:),   &
    & ddt_qi_fast(:,:,:),   &
    & ddt_qr_fast(:,:,:),   &
    & ddt_qs_fast(:,:,:),   &
    & ddt_qc_fast(:,:,:),   &
    & ddt_qg_fast(:,:,:),   &
    & ddt_qv(:,:,:),        &          ! tendencies fast and slow physics combined
    & ddt_qi(:,:,:),        &
    & ddt_qr(:,:,:),        &
    & ddt_qs(:,:,:),        &
    & ddt_qc(:,:,:),        &
    & ddt_qg(:,:,:)         &

    & => NULL()

  END TYPE t_sppt


END MODULE mo_sppt_types

