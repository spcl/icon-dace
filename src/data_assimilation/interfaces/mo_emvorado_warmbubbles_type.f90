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

MODULE mo_emvorado_warmbubbles_type

  USE mo_kind,                  ONLY: wp
  USE mo_impl_constants,        ONLY: max_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_warmbubble, t_bubblecontainer, autobubs_list

  TYPE t_warmbubble
    CHARACTER(len=12) :: ctype_tempdist  ! Type of perturbation 'cos-hrd', 'cos-instant'
    LOGICAL  ::   ltempdist              ! Switch to set temperature disturbance to ACTIVE
    LOGICAL  ::   ladd_bubblenoise_t     ! Switch to overlay random noise on the disturbance (not yet implemented)
    LOGICAL  ::   lbub_rhconst           ! Switch to activate a moisture increment such that rel. hum. stays constant during heating 
    REAL(wp) ::   htempdist              ! Time for beginning of temperature disturbance since model start time [s]
    REAL(wp) ::   centlon                ! Center (lon) of temperature disturbance [deg]
    REAL(wp) ::   centlat                ! Center (lat) of temperature disturbance [deg]
    REAL(wp) ::   centz                  ! Center (Z) of temperature disturbance [m]
    REAL(wp) ::   timespan               ! Total duration for release of temperature disturbance [s]
    REAL(wp) ::   timecounter            ! Actual duration during which the bubble was already active [s]
    REAL(wp) ::   radx                   ! Length scale / radius (X) of temperature disturbance [m]
    REAL(wp) ::   rady                   ! Length scale / radius (Y) of temperature disturbance [m]
    REAL(wp) ::   radz                   ! Length scale / radius (Z) of temperature disturbance [m]
    REAL(wp) ::   rotangle               ! Rotation angle of main axes of temperature disturbance [degrees]
    REAL(wp) ::   heatingrate            ! Constant heating rate [K/s] for 'cos-hrd' bubble
    REAL(wp) ::   dT                     ! Temperature increment of 'cos-instant' bubble [K]
    REAL(wp) ::   dT_bubblenoise         ! In case of ladd_bubblenoise_t=.true., relative noise level, such that
                                         !   dT          = dT          * (1 + dT_bubblenoise * random_noise[-1,1])   ('cos-instant')
                                         !   heatingrate = heatingrate * (1 + dT_bubblenoise * random_noise[-1,1])   ('cos-hrd')
  END TYPE t_warmbubble

  TYPE t_bubblecontainer
    INTEGER :: num_bubs = 0                   ! Number of elements of the list "bubs"
    TYPE(t_warmbubble), POINTER :: bubs(:) => NULL()
  END TYPE t_bubblecontainer
  
  TYPE(t_bubblecontainer) :: autobubs_list(max_dom)
  
END MODULE mo_emvorado_warmbubbles_type
