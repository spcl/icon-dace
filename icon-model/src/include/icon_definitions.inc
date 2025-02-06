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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!
#if defined(__NO_TIMERS__)
#define start_timer(timer,level) ! no timer
#define stop_timer(timer,level)  ! no timer
#else
#define start_timer(timer,level) IF (timers_level >= level) CALL timer_start(timer)
#define stop_timer(timer,level) IF (timers_level >= level) CALL timer_stop(timer)
#endif

#if defined(__NO_TIMERS__)
#define start_detail_timer(timer,level) ! no timer
#define stop_detail_timer(timer,level)  ! no timer
#else
#define start_detail_timer(timer,level) IF (timers_level >= level) CALL timer_start(timer)
#define stop_detail_timer(timer,level)  IF (timers_level >= level) CALL timer_stop(timer)
#endif


#if defined(__NO_SYNC_TIMERS__) || defined(__NO_TIMERS__)
#define start_sync_timer(timer) ! no timer
#define stop_sync_timer(timer)  ! no timer
#else
#define start_sync_timer(timer) IF (activate_sync_timers) CALL timer_start(timer)
#define stop_sync_timer(timer)  IF (activate_sync_timers) CALL timer_stop(timer)
#endif
