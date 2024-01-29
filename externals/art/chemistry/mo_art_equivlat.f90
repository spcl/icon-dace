!
! mo_art_equivlat
! This module provides the subroutines to calculate the
! equivalent latitude
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

MODULE mo_art_equivlat
  ! ICON
  USE mo_kind,                   ONLY: wp
  USE mo_math_constants,         ONLY: pi
  USE mo_physical_constants,     ONLY: earth_radius
  USE mo_model_domain,           ONLY: t_patch
  USE mo_sync,                   ONLY: global_min, global_max, global_sum_array

  ! ART
  USE mo_art_data,               ONLY: p_art_data
  USE mo_art_atmo_data,          ONLY: t_art_atmo
  USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
      
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif
      
  IMPLICIT NONE


  PRIVATE
  
  PUBLIC   :: art_equivalent_latitude       
         
CONTAINS

SUBROUTINE art_equivalent_latitude(jg, pv, equivLat_gridpoints, equivLat_isolines, n_isolines, &
  &                            select_hemisphere)
!<
! SUBROUTINE art_equivalent_latitude             
! This subroutine calculates the equivalent latitude
! Part of Module: mo_equivlat
! Author: Christian Stassen 
! Initial Release:  2014-10-27                
! Modifications:
!>

  INTEGER, INTENT(in)   ::  &
    &  jg                  !< patch id
  REAL(wp), INTENT(IN)  ::  &
    &  pv(:,:,:)                   !< Potential vorticity
  REAL(wp), INTENT(INOUT) ::  &
    &  equivLat_gridpoints(:,:,:)  !< equivalent Latitude at each gridpoint
  REAL(wp), INTENT(INOUT) :: &
    &  equivLat_isolines(:,:)      !< Isolines of equivalent latitude
  INTEGER, INTENT(IN)     :: &
    &  n_isolines                  !< Number of isolines
  CHARACTER(LEN=2), INTENT(IN)  ::   &
    &  select_hemisphere           !< Hemisphere on which equiv. lat. is computed

  !Local variables
  INTEGER ::                     &
     &    jc, jk, jb, i,         &                   !< loop indizes
     &    i_startidx, i_endidx

  REAL(wp),ALLOCATABLE :: minPV(:), &
                        & maxPV(:), &
                        & enclosed_area(:,:)

  REAL(wp),ALLOCATABLE :: pv_isolines(:,:)

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  ! ----------------------------------
  ! --- Start routine
  ! ----------------------------------

  art_atmo => p_art_data(jg)%atmo


  ! ----------------------------------
  ! --- Allocation
  ! ----------------------------------

  IF (.NOT. ALLOCATED(pv_isolines) )   ALLOCATE( pv_isolines(n_isolines,art_atmo%nlev)  )
  IF (.NOT. ALLOCATED(maxPV) )         ALLOCATE( maxPV(art_atmo%nlev)                   )
  IF (.NOT. ALLOCATED(minPV) )         ALLOCATE( minPV(art_atmo%nlev)                   )
  IF (.NOT. ALLOCATED(enclosed_area) ) ALLOCATE( enclosed_area(n_isolines,art_atmo%nlev))


  ! ----------------------------------
  ! --- Get Equivalent Latitude
  ! ----------------------------------

  ! ----------------------------------
  ! --- Find Max/Min PV at 
  ! --- hemisphere of interest
  ! ----------------------------------
  SELECT CASE(select_hemisphere)
    ! ----------------------------------
    ! --- First north hemisphere
    ! ----------------------------------
    CASE('NH')

    !Init values
    equivLat_isolines(:,:) = 0.
    maxPV(:) = -1.e20_wp
    minPV(:) = 1.e-20_wp
    enclosed_area(:,:) = 0.

     DO jb = art_atmo%i_startblk, art_atmo%i_endblk
       CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

       DO jk = 1,art_atmo%nlev
         DO jc = i_startidx,i_endidx
           IF (art_atmo%lat(jc,jb) > 0.) THEN
             IF (maxPV(jk) < pv(jc,jk,jb)) maxPV(jk) = pv(jc,jk,jb)
             IF (minPV(jk) > pv(jc,jk,jb)) minPV(jk) = pv(jc,jk,jb)
           END IF
         ENDDO
       ENDDO
     ENDDO

    maxPV(:) = global_max(maxPV(:))
    minPV(:) = global_min(minPV(:))

    ! ----------------------------------
    ! --- Calculate n_isoline values
    ! ----------------------------------

    DO jk = 1,art_atmo%nlev
      DO i = 1, n_isolines
        pv_isolines(i,jk) = minPV(jk) + REAL(i)*(maxPV(jk)-minPV(jk)) / (REAL(n_isolines))
      END DO
    END DO

    ! ----------------------------------
    ! --- Calculate geographic area
    ! --- enclosed by each pv_isoline
    ! ----------------------------------

     DO jb = art_atmo%i_startblk, art_atmo%i_endblk
       CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

       DO jk = 1,art_atmo%nlev
         DO i = 1, n_isolines
           DO jc = i_startidx,i_endidx
             ! <= is right for south hemisphere; should be >= for nh
             IF ( pv(jc,jk,jb) >= pv_isolines(i,jk) )        &
               &  enclosed_area(i,jk) = enclosed_area(i,jk)  &
               &  + art_atmo%cell_area(jc,jb) / (earth_radius * earth_radius)
           END DO
         ENDDO
       ENDDO
     ENDDO


     DO jk = 1,art_atmo%nlev
       DO i = 1, n_isolines
         enclosed_area(i,jk) = global_sum_array(enclosed_area(i,jk))
       ENDDO
     ENDDO

    ! ----------------------------------
    ! --- Calculate equivalent latitude 
    ! --- of isolines
    ! ----------------------------------

     DO jk = 1,art_atmo%nlev
       DO i = 1, n_isolines
         equivLat_isolines(i,jk) = 180./pi*asin(1.-(enclosed_area(i,jk)/(2.*pi)))
       ENDDO
     ENDDO

    ! ----------------------------------
    ! --- Get equivalent latitude 
    ! --- for each gridpoint 
    ! ----------------------------------

     DO jb = art_atmo%i_startblk, art_atmo%i_endblk
       CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

       DO jk=1,art_atmo%nlev
         DO i = 1, n_isolines
           DO jc=i_startidx,i_endidx
             IF   ( art_atmo%lat(jc,jb) > 0. .AND.  &
               &    ABS(pv(jc,jk,jb) - pv_isolines(i,jk)) <=    &
               &    (maxPV(jk)-minPV(jk))/(2*(n_isolines-1)) )  &
               &  equivLat_gridpoints(jc,jk,jb) = equivLat_isolines(i,jk)                
           END DO
         ENDDO
       ENDDO
     ENDDO

    ! ----------------------------------
    ! --- Second south hemisphere
    ! ----------------------------------
    CASE('SH')

    !Init values
    equivLat_isolines(:,:) = 0.
    maxPV(:) = -1.e20_wp
    minPV(:) = 1.e-20_wp
    enclosed_area(:,:) = 0.

     DO jb = art_atmo%i_startblk, art_atmo%i_endblk
       CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

       DO jk = 1,art_atmo%nlev
         DO jc = i_startidx,i_endidx
           IF (art_atmo%lat(jc,jb) < 0.) THEN
             IF(maxPV(jk) < pv(jc,jk,jb)) maxPV(jk) = pv(jc,jk,jb)
             IF(minPV(jk) > pv(jc,jk,jb)) minPV(jk) = pv(jc,jk,jb)
           END IF
         ENDDO
       ENDDO
     ENDDO

    maxPV(:) = global_max(maxPV(:))
    minPV(:) = global_min(minPV(:))

    ! ----------------------------------
    ! --- Calculate n_isoline values
    ! ----------------------------------

    DO jk = 1,art_atmo%nlev
      DO i = 1, n_isolines
        pv_isolines(i,jk) = minPV(jk) + REAL(i)*(maxPV(jk)-minPV(jk)) / (REAL(n_isolines))
      END DO
    END DO

    ! ----------------------------------
    ! --- Calculate geographic area
    ! --- enclosed by each pv_isoline
    ! ----------------------------------

     DO jb = art_atmo%i_startblk, art_atmo%i_endblk
       CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

       DO jk = 1,art_atmo%nlev
         DO i = 1, n_isolines
           DO jc = i_startidx,i_endidx
             ! <= is right for south hemisphere; should be >= for nh
             IF ( pv(jc,jk,jb) <= pv_isolines(i,jk) )        &
               &  enclosed_area(i,jk) = enclosed_area(i,jk)  &
               &  + art_atmo%cell_area(jc,jb) / (earth_radius * earth_radius)
           END DO
         ENDDO
       ENDDO
     ENDDO


     DO jk = 1,art_atmo%nlev
       DO i = 1, n_isolines
         enclosed_area(i,jk) = global_sum_array(enclosed_area(i,jk))
       ENDDO
     ENDDO

    ! ----------------------------------
    ! --- Calculate equivalent latitude 
    ! --- of isolines
    ! ----------------------------------

     DO jk = 1,art_atmo%nlev
       DO i = 1, n_isolines
         equivLat_isolines(i,jk) = -180./pi*asin(1.-(enclosed_area(i,jk)/(2.*pi)))
         IF( equivLat_isolines(i,jk) > 0.0_wp ) equivLat_isolines(i,jk) = 0.0_wp
       ENDDO
     ENDDO

    ! ----------------------------------
    ! --- Get equivalent latitude 
    ! --- for each gridpoint 
    ! ----------------------------------

     DO jb = art_atmo%i_startblk, art_atmo%i_endblk
       CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

       DO jk=1,art_atmo%nlev
         DO i = 1, n_isolines
           DO jc=i_startidx,i_endidx
             IF   ( art_atmo%lat(jc,jb) < 0. .AND.  &
               &    ABS(pv(jc,jk,jb) - pv_isolines(i,jk)) <=    &
               &    (maxPV(jk)-minPV(jk))/(2*(n_isolines-1)) )  &
               &  equivLat_gridpoints(jc,jk,jb) = equivLat_isolines(i,jk)        
           END DO
         ENDDO
       ENDDO
     ENDDO

  END SELECT

END SUBROUTINE art_equivalent_latitude
END MODULE mo_art_equivlat
