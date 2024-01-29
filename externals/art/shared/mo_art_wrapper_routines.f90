!
! mo_art_wrapper_routines
! This module collects the routines that use ICON information and therefore
! have to be available for another host model or have to be redone for another
! host model (or using ART as bix model)
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

MODULE mo_art_wrapper_routines
  USE mo_loopindices,                   ONLY: get_indices_c
  USE mo_impl_constants,                ONLY: min_rlcell_int          
  USE mo_impl_constants_grf,            ONLY: grf_bdywidth_c          
  USE mo_model_domain,                  ONLY: p_patch

  !ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_get_indices_c

CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_get_indices_c(jg, jb, i_startidx, i_endidx)
  INTEGER, INTENT(in)     ::  &
    &  jg, jb
  INTEGER, INTENT(inout)  ::  &
    &  i_startidx, i_endidx

  ! local variables
  INTEGER ::       &
    &  i_rlstart,  &
    &  i_rlend
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  i_rlstart   =  grf_bdywidth_c+1
  i_rlend     =  min_rlcell_int
  CALL get_indices_c(p_patch(jg), jb, art_atmo%i_startblk, art_atmo%i_endblk,   &
            &        i_startidx,  i_endidx, i_rlstart, i_rlend)
END SUBROUTINE art_get_indices_c

END MODULE mo_art_wrapper_routines
