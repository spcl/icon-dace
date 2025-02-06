!
! mo_art_emission_pntSrc
! This module provides an flexible emission routine for point sources.
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

MODULE mo_art_emission_pntSrc
! ICON
  USE mo_kind,                          ONLY: wp
  USE mtime,                            ONLY: datetime
! ART
  USE mo_art_pntSrc_types,              ONLY: t_art_all_pntSrc
  USE mo_art_impl_constants,            ONLY: UNDEF_INT_ART

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_emission_pntSrc'

  PUBLIC :: art_emission_pntSrc

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_emission_pntSrc(all_pntSrc, current_date, dtime, rho, cell_area, dz, tracer)
!<
! SUBROUTINE art_emission_pntSrc
! This subroutine provides an flexible emission routine for point sources.
! Based on: Werchner (2016) - Bachelorthesis, KIT
! Part of Module: mo_art_emission_pntSrc
! Author: Sven Werchner, Daniel Rieger, KIT
! Initial Release: 2017-01-26
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  TYPE(t_art_all_pntSrc), INTENT(inout) :: &
    &  all_pntSrc                         !< Container with information about all point sources
  TYPE(datetime), INTENT(in),POINTER :: &
    &  current_date                       !< mtime object containing current date
  REAL(wp), INTENT(in)               :: &
    &  dtime,                           & !< Model time step
    &  rho(:,:,:),                      & !< Density of air
    &  dz(:,:,:),                       & !< Layer height
    &  cell_area(:,:)                     !< Cell area
  REAL(wp), INTENT(inout)            :: &
    &  tracer(:,:,:,:)                    !< Mass mixing ratio of selected tracers
! Local variables
  REAL(wp)              :: &
    &  emiss_fct             !< dummy variable for emission for loop optimization
  INTEGER               :: &
    &  jsource,            & !< Counter for point source scenarios
    &  jlocsource,         & !< Counter for sources on PE domain
    &  jc, jk, jb,         & !< Location of source
    &  jk_top, jk_bot,     & !< Top and bottom vertical index of source
    &  itr,                & !< Index of tracer
    &  itr0                  !< Index for number concentration for aerosol tracer

  DO jsource = 1, all_pntSrc%nsources
    IF (all_pntSrc%p(jsource)%isActive(current_date)) THEN
      ! Abbreviations
      jc           = all_pntSrc%p(jsource)%tri_iidx_loc
      jb           = all_pntSrc%p(jsource)%tri_iblk_loc
      jk_top       = all_pntSrc%p(jsource)%k_index
      jk_bot       = all_pntSrc%p(jsource)%k_index_bot
      itr          = all_pntSrc%p(jsource)%itr
      itr0         = all_pntSrc%p(jsource)%itr0

      ! Only do emissions, if current PE domain includes this point source
      DO jlocsource = 1, all_pntSrc%p(jsource)%ithis_nlocal_pts

        emiss_fct = all_pntSrc%p(jsource)%source_strength * dtime / cell_area(jc,jb)
        ! for emission in single level jk_top == jk_bot
        DO jk = jk_top, jk_bot
          tracer(jc,jk,jb,itr) = tracer(jc,jk,jb,itr)                                             &   
            &                  + all_pntSrc%p(jsource)%height_factor(jk)                          &   
            &                  * emiss_fct / ( rho(jc,jk,jb) * dz(jc,jk,jb) )
        ENDDO

        IF (itr0 /= UNDEF_INT_ART) THEN
          emiss_fct = all_pntSrc%p(jsource)%emiss_rate0 * dtime / cell_area(jc,jb)
          DO jk = jk_top, jk_bot
            tracer(jc,jk,jb,itr0) = tracer(jc,jk,jb,itr0)                                         &   
              &                   + all_pntSrc%p(jsource)%height_factor(jk)                       &   
              &                   * emiss_fct / ( rho(jc,jk,jb) * dz(jc,jk,jb) )
          ENDDO
        ENDIF

      ENDDO !jlocsource
    ENDIF !pntSrc%isActive
  ENDDO !jsource
END SUBROUTINE art_emission_pntSrc
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_pntSrc
