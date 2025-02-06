
!
! mo_art_cell_loop
! This module provides the interface for looping over cells.
!
! These routines are based on atm_phy_aes/mo_omp_loop
! Orig author:Marco Giorgetta, MPI-M
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

MODULE mo_art_cell_loop

  USE mo_kind,                 ONLY: wp

  USE mtime,                   ONLY: datetime
  
  USE mo_art_chem_types,       ONLY: t_chem_meta_passive
  USE mo_art_chem_types_param, ONLY: t_chem_meta_lt, t_chem_meta_linoz, &
                                 &   t_chem_meta_simnoy
  USE mo_art_data,             ONLY: p_art_data
  USE mo_art_atmo_data,        ONLY: t_art_atmo
  USE mo_art_wrapper_routines, ONLY: art_get_indices_c


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: art_loop_cell_tracer, art_loop_cell_tracer_td
  PUBLIC :: art_loop_cell_array,  art_loop_cell_array_td



  INTERFACE art_loop_cell_array
    MODULE PROCEDURE art_loop_cell_array_3d
  END INTERFACE art_loop_cell_array

  INTERFACE art_loop_cell_array_td
    MODULE PROCEDURE art_loop_cell_td
  END INTERFACE art_loop_cell_array_td

  INTERFACE art_loop_cell_tracer
   MODULE PROCEDURE art_loop_cell_tracer_lt
   MODULE PROCEDURE art_loop_cell_tracer_linoz
   MODULE PROCEDURE art_loop_cell_tracer_simnoy
  END INTERFACE art_loop_cell_tracer


  INTERFACE art_loop_cell_tracer_td
   MODULE PROCEDURE art_loop_cell_tracer_td_lt
   MODULE PROCEDURE art_loop_cell_tracer_td_passive
   MODULE PROCEDURE art_loop_cell_tracer_td_linoz
   MODULE PROCEDURE art_loop_cell_tracer_td_simnoy
  END INTERFACE art_loop_cell_tracer_td

CONTAINS

  SUBROUTINE art_loop_cell_array_3d(jg ,array,&
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)     :: jg
    REAL(wp), INTENT(inout) :: array(:,:,:)
    
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             nproma,nlev,array)
         USE mo_kind,      ONLY: wp
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         INTEGER        ,INTENT(in) :: nproma,nlev
         REAL(wp), INTENT(inout) :: array(:,:,:)
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo



!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)
       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       art_atmo%nproma,art_atmo%nlev,array)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_array_3d


  SUBROUTINE art_loop_cell_tracer_lt(jg ,tracer,&
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_lt), INTENT(inout) :: tracer
    
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             nproma,nlev,tracer)
          USE mo_art_chem_types_param,   ONLY: t_chem_meta_lt

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         INTEGER        ,INTENT(in) :: nproma,nlev
         TYPE(t_chem_meta_lt), INTENT(inout) :: tracer
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo



!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)
       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       art_atmo%nproma,art_atmo%nlev,tracer)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_lt

 SUBROUTINE art_loop_cell_tracer_linoz(jg  ,tracer,&
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_linoz), INTENT(inout) :: tracer
    
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             nproma,nlev,tracer)
          USE mo_art_chem_types_param,   ONLY: t_chem_meta_linoz

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         INTEGER        ,INTENT(in) :: nproma,nlev
         TYPE(t_chem_meta_linoz), INTENT(inout) :: tracer
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo


!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       art_atmo%nproma,art_atmo%nlev,tracer)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_linoz


 SUBROUTINE art_loop_cell_tracer_simnoy(jg  ,tracer,&
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_simnoy), INTENT(inout) :: tracer
    
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             nproma,nlev,tracer)
          USE mo_art_chem_types_param,   ONLY: t_chem_meta_simnoy

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         INTEGER        ,INTENT(in) :: nproma,nlev
         TYPE(t_chem_meta_simnoy), INTENT(inout) :: tracer
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo


!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       art_atmo%nproma,art_atmo%nlev,tracer)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_simnoy

 SUBROUTINE art_loop_cell_tracer_td_linoz(jg  ,tracer,&
       &                   current_date, &
       &                   p_dtime, &
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_linoz), INTENT(inout) :: tracer
    
    TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
    REAL(wp), INTENT(IN) :: p_dtime
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             tracer,current_date, p_dtime)
          USE mo_art_chem_types_param,   ONLY: t_chem_meta_linoz
          USE mtime,                     ONLY: datetime
          USE mo_kind,                   ONLY: wp
          
          

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         TYPE(t_chem_meta_linoz), INTENT(inout) :: tracer
         TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
         REAL(wp), INTENT(IN) :: p_dtime
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       tracer,current_date, p_dtime)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_td_linoz


 SUBROUTINE art_loop_cell_tracer_td_simnoy(jg  ,tracer,&
       &                   current_date, &
       &                   p_dtime, &
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_simnoy), INTENT(inout) :: tracer
    
    TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
    REAL(wp), INTENT(IN) :: p_dtime
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             tracer,current_date, p_dtime)
          USE mo_art_chem_types_param,   ONLY: t_chem_meta_simnoy
          USE mtime,                     ONLY: datetime
          USE mo_kind,                   ONLY: wp
          
          

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         TYPE(t_chem_meta_simnoy), INTENT(inout) :: tracer
         TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
         REAL(wp), INTENT(IN) :: p_dtime
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       tracer,current_date, p_dtime)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_td_simnoy

 SUBROUTINE art_loop_cell_tracer_td_passive(jg  ,tracer,&
       &                   current_date, &
       &                   p_dtime, &
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_passive), INTENT(inout) :: tracer
    
    TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
    REAL(wp), INTENT(IN) :: p_dtime
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             tracer,current_date, p_dtime)
          USE mo_art_chem_types,         ONLY: t_chem_meta_passive
          USE mtime,                     ONLY: datetime
          USE mo_kind,                   ONLY: wp
          
          

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         TYPE(t_chem_meta_passive), INTENT(inout) :: tracer
         TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
         REAL(wp), INTENT(IN) :: p_dtime
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       tracer,current_date, p_dtime)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_td_passive

 SUBROUTINE art_loop_cell_tracer_td_lt(jg  ,tracer,&
       &                   current_date, &
       &                   p_dtime, &
       &                   routine)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(t_chem_meta_lt), INTENT(inout) :: tracer
    
    TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
    REAL(wp), INTENT(IN) :: p_dtime
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             tracer,current_date, p_dtime)
          USE mo_art_chem_types_param,   ONLY: t_chem_meta_lt
          USE mtime,                     ONLY: datetime
          USE mo_kind,                   ONLY: wp
          
          

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         TYPE(t_chem_meta_lt), INTENT(inout) :: tracer
         TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
         REAL(wp), INTENT(IN) :: p_dtime
         
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo

!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       tracer,current_date, p_dtime)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_tracer_td_lt


  SUBROUTINE art_loop_cell_td(jg,current_date, &
       &                   routine, array_in)

    ! Arguments
    !
    INTEGER, INTENT(in)                 :: jg
    TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
    REAL(wp), OPTIONAL, INTENT(INOUT) :: array_in(:,:,:)
    
    !
    INTERFACE

    
       !
       SUBROUTINE routine(jg,jb,jcs,jce ,&
            &             current_date, array_in)

          USE mtime,                   ONLY: datetime
          USE mo_kind                 ,ONLY: wp

        
         !
         INTEGER        ,INTENT(in) :: jg,jb,jcs,jce
         
         TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< actual date
         REAL(wp), OPTIONAL, INTENT(INOUT) :: array_in(:,:,:)
         

         !
       END SUBROUTINE routine
       !
    END INTERFACE

    ! Local variables
    !
    INTEGER  :: jb          !< index of block loop
    INTEGER  :: jcs, jce    !< start and end indices of columns loop
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo

    art_atmo => p_art_data(jg)%atmo


!$OMP PARALLEL DO PRIVATE(jb,jcs,jce)
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
       !
       CALL art_get_indices_c(jg,jb,jcs,jce)

       IF (jcs>jce) CYCLE
       !
       CALL routine(jg,jb,jcs,jce ,&
            &       current_date, array_in)
       !
    END DO
!$OMP END PARALLEL DO 

  END SUBROUTINE art_loop_cell_td




END MODULE mo_art_cell_loop


