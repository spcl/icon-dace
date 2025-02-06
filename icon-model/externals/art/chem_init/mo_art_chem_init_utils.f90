!
! mo_art_chem_init_utils
! This module provides the utilities for external data set interpolation
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_art_chem_init_utils
  ! ICON
  USE mo_impl_constants,      ONLY: MAX_CHAR_LENGTH
  USE mo_exception,           ONLY: finish

  ! ART
  USE mo_art_data,            ONLY: p_art_data
  USE mo_art_atmo_data,       ONLY: t_art_atmo
  USE mo_art_chem_init_types, ONLY: t_chem_init_state



  IMPLICIT NONE


  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_art_chem_init_utils'


  PUBLIC :: allocate_chem_init_atm, allocate_chem_init_chem
  PUBLIC :: deallocate_chem_init_atm, deallocate_chem_init_chem

CONTAINS



SUBROUTINE allocate_chem_init_atm (chem_init,jg, nlev_in_chem_init)
!<
! SUBROUTINE allocate_chem_init_atm                   
! This subroutine allocates arrays for external dataset interpolation
! Part of Module: mo_art_chem_init_utils
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-11-15                
! Modifications:
!>


  INTEGER, INTENT(IN) :: jg                   !< patch id, number of blocks centre and edges
  INTEGER, INTENT(IN) :: nlev_in_chem_init    !< number of model levels external data

  TYPE(t_chem_init_state),INTENT(inout)    :: &
              &  chem_init                    !< ART chem_init fields

  ! local variables
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
        &   routine = modname//':allocate_chem_init_atm'
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo


  IF (nlev_in_chem_init == 0) THEN
    CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
  END IF
  
  
  
  ! Allocate atmospheric input data
  ALLOCATE(chem_init%chem_init_in%temp     (art_atmo%nproma,nlev_in_chem_init,art_atmo%nblks))
  ALLOCATE(chem_init%chem_init_in%q        (art_atmo%nproma,nlev_in_chem_init,art_atmo%nblks))
  ALLOCATE(chem_init%chem_init_in%temp_v   (art_atmo%nproma,nlev_in_chem_init,art_atmo%nblks))
  ALLOCATE(chem_init%chem_init_in%pres     (art_atmo%nproma,nlev_in_chem_init,art_atmo%nblks))
  ALLOCATE(chem_init%chem_init_in%z3d      (art_atmo%nproma,nlev_in_chem_init,art_atmo%nblks))
  ALLOCATE(chem_init%chem_init_in%phi_sfc  (art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(chem_init%chem_init_in%ps       (art_atmo%nproma,art_atmo%nblks))


  chem_init%chem_init_in%linitialized = .TRUE.

END SUBROUTINE allocate_chem_init_atm


SUBROUTINE allocate_chem_init_chem (chem_init,jg,   &
                       &            nlev_in_chem_init)
!<
! SUBROUTINE allocate_chem_init_chem                   
! This subroutine allocates arrays for external dataset interpolation
! Part of Module: mo_art_chem_init_utils
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-11-15                
! Modifications:
!>


  INTEGER, INTENT(IN)  :: jg                    !< patch id
  INTEGER, INTENT(IN)  :: nlev_in_chem_init     !< number of model levels in external data

  TYPE(t_chem_init_state), INTENT(inout)    :: &
    &  chem_init                    !< ART chem_init fields

  ! Local variables
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
       &    routine = modname//':allocate_chem_init_chem'
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo   !< Pointer to ART atmo fields

  art_atmo => p_art_data(jg)%atmo
    
  IF (nlev_in_chem_init == 0) THEN
    CALL finish(routine, "Number of input levels <nlev_in> not yet initialized.")
  END IF

  ALLOCATE(chem_init%chem_init_chem%spec(art_atmo%nproma,                        &
                         &               nlev_in_chem_init,                      &
                         &               art_atmo%nblks,                         &
                         &               chem_init%chem_init_chem%n_spec_readin))

  ALLOCATE(chem_init%chem_init_chem%spec_interp(art_atmo%nproma,                      &
                         &                      art_atmo%nlev,                        &
                         &                      art_atmo%nblks,                       &
                         &                      chem_init%chem_init_chem%n_spec_readin))

  chem_init%chem_init_chem%linitialized = .TRUE.

END SUBROUTINE allocate_chem_init_chem


SUBROUTINE deallocate_chem_init_atm (chem_init)
!<
! SUBROUTINE deallocate_chem_init_atm                   
! This subroutine deallocates arrays for external dataset interpolation
! Part of Module: mo_art_chem_init_utils
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-15                
! Modifications:
!>
  TYPE(t_chem_init_state),INTENT(inout)    :: &
              &  chem_init                    !< ART chem_init fields

  ! local variables
  
  
  ! Deallocate atmospheric input data
  IF (ALLOCATED(chem_init%chem_init_in%temp)) DEALLOCATE(chem_init%chem_init_in%temp)
  IF (ALLOCATED(chem_init%chem_init_in%q)) DEALLOCATE(chem_init%chem_init_in%q)
  IF (ALLOCATED(chem_init%chem_init_in%temp_v)) DEALLOCATE(chem_init%chem_init_in%temp_v)
  IF (ALLOCATED(chem_init%chem_init_in%pres)) DEALLOCATE(chem_init%chem_init_in%pres)
  IF (ALLOCATED(chem_init%chem_init_in%z3d)) DEALLOCATE(chem_init%chem_init_in%z3d)
  IF (ALLOCATED(chem_init%chem_init_in%phi_sfc)) DEALLOCATE(chem_init%chem_init_in%phi_sfc)
  IF (ALLOCATED(chem_init%chem_init_in%ps)) DEALLOCATE(chem_init%chem_init_in%ps)


  chem_init%chem_init_in%linitialized = .FALSE.

END SUBROUTINE deallocate_chem_init_atm


SUBROUTINE deallocate_chem_init_chem (chem_init)
!<
! SUBROUTINE deallocate_chem_init_chem                   
! This subroutine deallocates arrays for external dataset interpolation
! Part of Module: mo_art_chem_init_utils
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-15                
! Modifications:
!>


  TYPE(t_chem_init_state), INTENT(inout)    :: &
    &  chem_init                    !< ART chem_init fields

  ! Local variables
  DEALLOCATE(chem_init%chem_init_chem%spec)

  DEALLOCATE(chem_init%chem_init_chem%spec_interp)

  chem_init%chem_init_chem%linitialized = .FALSE.

END SUBROUTINE deallocate_chem_init_chem
END MODULE mo_art_chem_init_utils

