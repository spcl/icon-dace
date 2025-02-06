!
! mo_art_fplume_write
! This module stores variables for output in external types structure.
!
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

MODULE mo_art_fplume_write
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_kind,                ONLY: wp
  USE mo_art_external_types,  ONLY: t_art_volc_fplume
  USE mo_parallel_config,     ONLY: nproma
  USE mo_model_domain,        ONLY: p_patch

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_write_fplume_values

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_write_fplume_values(fplume,jg,height,MER,iidx,iblks,fplume_on)

  TYPE(t_art_volc_fplume),INTENT(INOUT):: &
    &  fplume
  INTEGER, INTENT(IN)           :: jg,iidx,iblks
  REAL(wp), INTENT(IN)          :: height,MER
  LOGICAL, INTENT(IN)           :: fplume_on


    IF (fplume_on) THEN
      fplume%plume_H(iidx,iblks)         = height
      fplume%plume_MER(iidx,iblks)       = MER
    ELSE
      fplume%plume_H(iidx,iblks)                     = 0.0_wp
      fplume%plume_MER(iidx,iblks)                   = 0.0_wp
    ENDIF

  RETURN
END SUBROUTINE art_write_fplume_values
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_fplume_write
