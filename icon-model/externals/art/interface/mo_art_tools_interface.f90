!
! Provides interface to the ART-routine for using tools
!
! This module provides an interface to the ART-routine art_ini_tracer.
! The interface is written in such a way, that ICON will compile and run
! properly, even if the ART-routines are not available at compile time.
!
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

MODULE mo_art_tools_interface
  USE mo_kind,                          ONLY: wp
  USE mo_run_config,                    ONLY: lart
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art, timer_art_toolInt
  USE mo_art_unit_conversion,           ONLY: art_massmix2density

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC  :: art_tools_interface 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_tools_interface(defcase, prog_list, tracer_now, tracer_new, rho)
  !>
  !! Interface for ART tools
  !!
  !! @par Revision History
  !! Initial revision by Daniel Rieger, KIT (2013-12-16)
  
  CHARACTER(len=*),INTENT(in)        :: & 
    &  defcase                            !< definition of case 
  TYPE(t_var_list_ptr),TARGET,INTENT(in) :: &
    &  prog_list                          !< prognostic state list
  REAL(wp),INTENT(in)                :: &
    &  tracer_now(:,:,:,:),             & !< tracer concentrations at timelevel nnow
    &  rho(:,:,:)                         !< Density
  REAL(wp),INTENT(inout)             :: &
    &  tracer_new(:,:,:,:)                !< tracer concentrations at timelevel nnew
  
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art)
    IF (timers_level > 3) CALL timer_start(timer_art_toolInt)

    IF (TRIM(defcase) .EQ. 'unit_conversion') THEN
      CALL art_massmix2density(prog_list, tracer_now, tracer_new, rho)
    ENDIF

    IF (timers_level > 3) CALL timer_stop(timer_art_toolInt)
    IF (timers_level > 3) CALL timer_stop(timer_art)
  ENDIF
    
END SUBROUTINE art_tools_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tools_interface
