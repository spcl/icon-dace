!
! Provides interface to the ART-routine for initializing ART-tracer compounds to the tracer field.
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

MODULE mo_art_tracer_interface
  USE mo_kind,                          ONLY: wp, i8
  USE mo_exception,                     ONLY: message,finish
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_fortran_tools,                 ONLY: t_ptr_2d3d
  USE mo_advection_config,              ONLY: t_advection_config
  USE mo_nwp_phy_types,                 ONLY: t_nwp_phy_tend
  USE mo_nonhydro_types,                ONLY: t_nh_prog
  USE mo_run_config,                    ONLY: lart, iforcing
  USE mo_impl_constants,                ONLY: iaes
  USE mo_radiation_config,              ONLY: irad_o3
  USE mo_aes_rad_config,              ONLY: aes_rad_config
  USE mtime,                            ONLY: MAX_TIMEDELTA_STR_LEN,                            &
                                          &   timedelta, newTimedelta, deallocateTimedelta,     &
                                          &   getTotalMilliSecondsTimeDelta, getPTStringFromMS
  USE mo_timer,                         ONLY: timers_level, timer_start, timer_stop,   &
                                          &   timer_art_tracInt
  USE mo_art_tracer,                    ONLY: art_tracer
  USE mo_art_init,                      ONLY: art_init

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_tracer_interface 

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_tracer_interface(defcase,jg,nblks_c,this_list,vname_prefix, &
   &                            ptr_arr,advconf,phy_tend, p_prog,          &
   &                            timelev,ldims)
!! Interface for ART-routine art_ini_tracer 
!!
!! This interface calls the ART-routine art_ini_tracer, if ICON has been 
!! built including the ART-package. Otherwise, this is simply a dummy 
!! routine.
!!
!! @par Revision History
!! Initial revision by Kristina Lundgren, KIT (2012-04-03)
  INTEGER,INTENT(in)             :: &
    &   jg,                         & !< patch id
    &   nblks_c                       !< patch block
  TYPE(t_var_list_ptr),INTENT(INOUT) :: &
    &   this_list                     !< current list: prognostic or phys. tend.  
  TYPE(t_ptr_2d3d),INTENT(inout),OPTIONAL         :: &
    &   ptr_arr(:)                    !< pointer to each element in list
  TYPE(t_advection_config),INTENT(inout),OPTIONAL ::  &
    &   advconf                       !< advection config
      
  TYPE(t_nwp_phy_tend),INTENT(inout),OPTIONAL     :: &
    &   phy_tend                      !< physical tendencies
      
  TYPE(t_nh_prog),INTENT(inout),OPTIONAL :: &
    &   p_prog                        !< prognostic variables
      
  INTEGER,INTENT(in), OPTIONAL   :: &
    &   timelev,                    & !< drieg : why is timelevel optional?
    &   ldims(3)                      !< local dimensions, for checking

  CHARACTER(len=*), INTENT(IN)   :: & 
    &   vname_prefix                  !< list name
  CHARACTER(len=*), INTENT(IN)   :: & 
    &   defcase                       !< definition of case 

  !-----------------------------------------------------------------------
 
  IF (lart) THEN
    IF (timers_level > 3) CALL timer_start(timer_art_tracInt)

    IF (TRIM(defcase) .NE. 'diag') THEN

      IF (.NOT. PRESENT(advconf)) THEN
              CALL finish('mo_art_tracer_interface:art_tracer_interface',  &
         &      'advconf not present at defcase'//TRIM(defcase))
      ENDIF
      IF (.NOT. PRESENT(ptr_arr)) THEN
              CALL finish('mo_art_tracer_interface:art_tracer_interface',  &
         &      'ptr_arr not present at defcase'//TRIM(defcase))
      ENDIF
      
      CALL message('','ART: Definition of tracers for defcase: '//TRIM(defcase))
        
      IF (TRIM(defcase) .EQ. 'prog') THEN 
        CALL art_tracer(defcase, jg, nblks_c, this_list, vname_prefix, ptr_arr, advconf,  & 
          &             p_prog=p_prog,timelev=timelev, ldims=ldims)
      ELSE
        CALL art_tracer(defcase, jg, nblks_c, this_list, vname_prefix, ptr_arr, advconf,  &
          &             phy_tend=phy_tend, ldims=ldims)
      ENDIF
    
      IF (TRIM(defcase) .EQ. 'prog') THEN
        IF (timelev .EQ. 1) THEN 
          IF ( iforcing == iaes) THEN
            irad_o3 = aes_rad_config(jg)%irad_o3
          END IF
        END IF
      END IF

    ENDIF

    IF (timers_level > 3) CALL timer_stop(timer_art_tracInt)
  ENDIF ! lart

END SUBROUTINE art_tracer_interface 
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tracer_interface

