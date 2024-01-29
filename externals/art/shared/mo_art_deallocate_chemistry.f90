!
! mo_art_deallocate_chemistry
! This module provides a routine to clean up chemistry structures
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

MODULE mo_art_deallocate_chemistry
! ICON Routines
  USE mo_art_config,                    ONLY: art_config
! ART Routines
  USE mo_art_chem_data,                 ONLY: t_art_chem
  USE mo_art_chem_init_meta,            ONLY: art_delete_OH_chem_meta
  USE mo_art_read_linoz,                ONLY: art_linoz_deallocate
  USE mo_art_photo_init,                ONLY: art_photo_deallocate
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_deallocate_chemistry

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_deallocate_chemistry(jg, chem)
!<
! SUBROUTINE art_deallocate_chemistry
! This subroutine does the clean up for ART chemistry structures
! Based on: -
! Part of Module: mo_art_deallocate_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-15
! Modifications:
!>
  INTEGER, INTENT(in) :: &
    &  jg
  TYPE(t_art_chem),INTENT(inout) :: &
    &  chem              !< number of model domains

  IF (chem%param%OH_chem_meta%is_init) THEN
    CALL art_delete_OH_chem_meta(chem%param%OH_chem_meta)
  END IF

  IF (chem%param%linoz%is_init) THEN
    CALL art_linoz_deallocate(chem%param%linoz)
  END IF

  IF (chem%param%OH_chem_meta%is_init .OR. &
    &  art_config(jg)%lart_mecca) THEN

    CALL art_photo_deallocate(chem%photo)
  END IF

  IF (art_config(jg)%lart_chemtracer) THEN
    DEALLOCATE(chem%param%o2_column)
  END IF

  IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
    DEALLOCATE(chem%mecicon%utils%mapping_indices_kpp)
#endif
  END IF

  DEALLOCATE(chem%water_tracers)
  DEALLOCATE(chem%CO2_mmr_depos)

  DEALLOCATE(chem%vmr2Nconc)
  
 
    
END SUBROUTINE art_deallocate_chemistry
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_deallocate_chemistry
