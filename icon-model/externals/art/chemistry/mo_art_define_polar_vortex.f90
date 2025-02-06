!
! mo_define_polar_vortex
! This module provides the subroutines to define the Polar Vortex
! using Nash Criteria 1996
!
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

MODULE mo_art_define_polar_vortex
  USE mo_kind,                        ONLY: wp
  USE mo_model_domain,                ONLY: t_patch
  USE mo_util_sort,                   ONLY: quicksort
  USE mo_sync,                        ONLY: global_sum_array
  USE mo_mpi,                         ONLY: get_mpi_all_workroot_id, get_my_mpi_work_id,      &
                                        &   my_process_is_mpi_workroot, p_comm_work,          &
                                        &   my_process_is_work, get_my_mpi_work_comm_size,    &
                                        &   p_send, p_recv, p_scatterv, p_gatherv, p_bcast

  USE mtime,                          ONLY: datetime, newDatetime

  ! ART
  USE mo_art_data,                    ONLY: p_art_data
  USE mo_art_atmo_data,               ONLY: t_art_atmo
  USE mo_art_wrapper_routines,        ONLY: art_get_indices_c
  USE mo_art_potential_vorticity,     ONLY: art_get_potential_vorticity
  USE mo_art_equivlat,                ONLY: art_equivalent_latitude
  USE mo_art_config,                  ONLY: art_config
      
  IMPLICIT NONE

  INTEGER, SAVE ::  &
    &  no_of_vortex_inits      !< Number of vortex init dates
  CHARACTER(LEN=2), ALLOCATABLE, SAVE ::  &
    &  hemisphere(:)           !< hemisphere on which equivLat is computed
  
  TYPE t_vortex_init_date
    TYPE(datetime), POINTER :: p
  END TYPE t_vortex_init_date
  
  TYPE(t_vortex_init_date),ALLOCATABLE ::  &
    &  vortex_init_date(:)    !< Date to intialize the vortextracer

  PRIVATE
  
  PUBLIC   :: no_of_vortex_inits,               &
      &       hemisphere,                       &
      &       vortex_init_date,                 &
      &       art_init_polar_vortex_boundary,   &
      &       art_define_polar_vortex_boundary   
         
CONTAINS

! --- Routine to get the intial date for the vortex
! --- tracer. This routine is only called once 
! --- at model start
! --- Author: Christian Stassen 24.03.2015
! --- Changelog: -
SUBROUTINE art_init_polar_vortex_boundary(jg)
  !<
  ! SUBROUTINE art_init_polar_vortex_boundary
  ! Routine to get the initial date for the vortex
  ! tracer. This routine is only called once 
  ! --- at model start
  ! Part of Module: mo_define_polar_vortex
  ! Author: Christian Stassen, KIT
  ! Initial Release: 2015-03-24
  !>

  CHARACTER(LEN=120)                :: vortex_init_dateandhemis_read_in
  CHARACTER(LEN=120), ALLOCATABLE   :: vortex_init_dateandhemis_string(:)
  CHARACTER(LEN=120), ALLOCATABLE   :: vortex_init_date_string(:)

  !Local variables
  INTEGER :: jg                  !< loop index
  INTEGER :: i, j, k, ivortex
  CHARACTER(LEN=4)  :: cyear
  CHARACTER(LEN=16) :: nh_default_date, sh_default_date


  !vortex_init_dateandhemis_read_in = TRIM(adjustl(art_config(jg)%cart_vortex_init_date))
  !Until used in namelist, only hardcoded vortex_init
  vortex_init_dateandhemis_read_in = 'NH-2013-11-05T00:00:00Z'


  IF ( TRIM(vortex_init_dateandhemis_read_in) .EQ. '' ) THEN
    no_of_vortex_inits = 0

  ELSE
    ! ----------------------------------
    ! --- Read in the string specifier for vortex
    ! ----------------------------------

    no_of_vortex_inits = 1

    DO i=1,LEN_TRIM(vortex_init_dateandhemis_read_in)
      IF ( TRIM(vortex_init_dateandhemis_read_in(i:i)) .EQ. ' ' ) THEN
        no_of_vortex_inits = no_of_vortex_inits + 1
      END IF
    END DO

    ! ----------------------------------
    ! --- Allocate vortex_init_date with the 
    ! --- number of dates found
    ! ----------------------------------

    IF ( .NOT. ALLOCATED(hemisphere) )       ALLOCATE(hemisphere(no_of_vortex_inits))
    IF ( .NOT. ALLOCATED(vortex_init_date) ) ALLOCATE(vortex_init_date(no_of_vortex_inits))


    ALLOCATE(vortex_init_date_string(no_of_vortex_inits))
    ALLOCATE(vortex_init_dateandhemis_string(no_of_vortex_inits))

    ! ----------------------------------
    ! --- Select hermisphere and date
    ! ----------------------------------
    j = 1
    k = 1
    DO i=1,LEN_TRIM(vortex_init_dateandhemis_read_in)
      IF ( trim(adjustl(vortex_init_dateandhemis_read_in(i:i))) == ' ' ) THEN
        j = j + 1
        k = 1
      ELSE
        vortex_init_dateandhemis_string(j)(k:k) = TRIM(vortex_init_dateandhemis_read_in(i:i))
        k = k + 1
      END IF
    END DO

    hemisphere(:) = vortex_init_dateandhemis_string(:)(1:2)
!    vortex_init_date_string(:) = vortex_init_dateandhemis_string(:)(4:k-1)
    vortex_init_date_string(:) = vortex_init_dateandhemis_string(:)(4:k-1)

    ! ----------------------------------
    ! --- Convert the string to type datetime
    ! ----------------------------------
    DO j=1, no_of_vortex_inits
      vortex_init_date(j)%p => newDatetime(trim(adjustl(vortex_init_date_string(j))))
    END DO

  END IF

END SUBROUTINE art_init_polar_vortex_boundary


! --- Routine to set the intial values for the vortex
! --- tracer. This routine is called every time
! --- act. date == vortex_date
! --- Author: Christian Stassen 24.03.2015
! --- Changelog: -
SUBROUTINE art_define_polar_vortex_boundary(current_date,ptr_vortex_log,  &
                          &                 select_hemisphere,p_patch, opt_lev_bot,opt_lev_top)
  !<
  ! SUBROUTINE define_polar_vortex_boundary
  ! Routine to get the initial date for the vortex
  ! tracer. This routine is called every time
  ! act. date == vortex_date
  ! Part of Module: mo_define_polar_vortex
  ! Author: Christian Stassen, KIT
  ! Initial Release: 2015-03-24
  !>
  TYPE(t_patch), TARGET, INTENT(IN) ::  &
    &  p_patch                !< patch on which computation is performed
  TYPE(datetime),POINTER,INTENT(IN) :: &
    &  current_date
  INTEGER, INTENT(IN), OPTIONAL  :: &
    &  opt_lev_bot, opt_lev_top
  INTEGER ::    &
    &  lev_bot, &
    &  lev_top
  LOGICAL, ALLOCATABLE, INTENT(INOUT) ::  &
    &  ptr_vortex_log(:,:,:)    !< Logical: True in vortex

  !Local variables
  INTEGER  ::  &
    &  jg                  !< loop index
  INTEGER ::                     &
     &    jc, jk, jb, i, j,      &                   !< loop indizes
     &    i_startidx, i_endidx                       !< index of start and end block 
  REAL(wp),ALLOCATABLE ::  &
    &  pv(:,:,:)                   !< Potential vorticity
  REAL(wp), ALLOCATABLE ::  &
    &  nash_equivLat(:)           !< Equivalent latitde at Nash-Criteria
  REAL(wp), ALLOCATABLE ::  &
    &  nash_crit(:,:,:)           !< Nash-Criteria
  REAL(wp),ALLOCATABLE ::       &
    &  equivLat_isolines(:,:),  &
    &  equivLat_gridpoints(:,:,:)  !< equivalent latitude & isolines
  REAL(wp), ALLOCATABLE :: &
    &  windspeed(:,:,:)           !< abs. horiz. windspeed
  REAL(wp) ::              &
    &  min_equivLat  = 60.   !< min. equivLat to start search for Nash (sign doesn't matter)
  REAL(wp) ::                  &
    &  breakup_windspeed  = 20.   !< min. windspeed for vortex
  INTEGER  ::                    &
    &  n_isolines          = 200   !< number of isolines
  CHARACTER(LEN=2) ::     &
    &  select_hemisphere           !< hemisphere on which equivLat is computed
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo           !< ART atmo fields

  ! ----------------------------------
  ! --- Start routine
  ! ----------------------------------

  jg        = p_patch%id
  art_atmo  => p_art_data(jg)%atmo

  ! ----------------------------------
  ! --- Get the loop indices
  ! ----------------------------------
    

  ! ----------------------------------
  ! --- Check for optional arguments
  ! ----------------------------------

  IF (PRESENT(opt_lev_bot)) THEN
    lev_bot = opt_lev_bot
  ELSE
    lev_bot = art_atmo%nlev
  END IF

  IF(PRESENT(opt_lev_top)) THEN
    lev_top = opt_lev_top
  ELSE
    lev_top = 1
  END IF

  ! ----------------------------------
  ! --- Allocation
  ! ----------------------------------
  IF ( .NOT. ALLOCATED(ptr_vortex_log) ) THEN
    ALLOCATE( ptr_vortex_log(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks) )
  END IF

  ptr_vortex_log(:,:,:)       = .FALSE.
  select_hemisphere           = 'NA'

  ! ----------------------------------
  ! --- Check if act. date = vortex init date
  ! --- If True => calculate the vortex
  ! --- If False => do nothing (vortex is false everywhere)
  ! ----------------------------------

  IF ( no_of_vortex_inits .NE. 0 ) THEN
    DO j = 1, no_of_vortex_inits !no_of_vortex_inits
      IF ( (vortex_init_date(j)%p%date%year  == current_date%date%year            &
         &       .OR.  vortex_init_date(j)%p%date%year  == 0)                     &
         &  .AND. (vortex_init_date(j)%p%date%month == current_date%date%month    &
         &       .OR.  vortex_init_date(j)%p%date%month == 0)                     &
         &  .AND. (vortex_init_date(j)%p%date%day   == current_date%date%day      &
         &       .OR.  vortex_init_date(j)%p%date%day   == 0)                     &
         &  .AND. (vortex_init_date(j)%p%time%hour  == current_date%time%hour) )  &
      &  THEN !vortex_init_date==act.date

      ! ----------------------------------
      ! --- Allocation
      ! ----------------------------------

      IF (.NOT. ALLOCATED(pv))   &
        &  ALLOCATE( pv(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
      IF (.NOT. ALLOCATED(nash_equivLat))  &
        &  ALLOCATE( nash_equivLat(art_atmo%nlev))
      IF (.NOT. ALLOCATED(nash_crit)) &
        &  ALLOCATE( nash_crit(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
      IF (.NOT. ALLOCATED(windspeed)) &
        &  ALLOCATE( windspeed(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
      IF (.NOT. ALLOCATED(equivLat_isolines)) &
        &  ALLOCATE( equivLat_isolines(n_isolines,art_atmo%nlev))
      IF (.NOT. ALLOCATED(equivLat_gridpoints))   &
        &  ALLOCATE( equivLat_gridpoints(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))

      pv(:,:,:)                   = 0.
      nash_equivLat(:)            = 0.
      nash_crit(:,:,:)            = 0.
      windspeed(:,:,:)            = 0.
      equivLat_isolines(:,:)      = 0. 
      equivLat_gridpoints(:,:,:)  = 0.

      ! ----------------------------------
      ! --- Set hemisphere
      ! ----------------------------------

      select_hemisphere = hemisphere(j)

      ! ----------------------------------
      ! --- Calculate Potential Vorticity (PV)
      ! ----------------------------------

      CALL art_get_potential_vorticity(p_patch, art_atmo%vor, art_atmo%dpres_mc, pv)

      ! ----------------------------------
      ! --- Get Equivalent Latitude
      ! ----------------------------------

      CALL art_equivalent_latitude(jg, pv, equivLat_gridpoints, equivLat_isolines, n_isolines, &
        &                      select_hemisphere)

      ! ----------------------------------
      ! --- Get abs. horiz. windspeed 
      ! --- on NH/SH
      ! ----------------------------------

       DO jb = art_atmo%i_startblk, art_atmo%i_endblk
         CALL art_get_indices_c(jg, jb,  i_startidx, i_endidx)

         DO jk=lev_top, lev_bot
           DO jc=i_startidx,i_endidx
             IF( select_hemisphere == 'SH' .AND. p_patch%cells%center(jc,jb)%lat < 0.) THEN
                 windspeed(jc,jk,jb) = sqrt( art_atmo%u(jc,jk,jb)*art_atmo%u(jc,jk,jb)   &
                   &                       + art_atmo%v(jc,jk,jb)*art_atmo%v(jc,jk,jb) )
             ELSE IF( select_hemisphere == 'NH' .AND. p_patch%cells%center(jc,jb)%lat > 0.) THEN
                 windspeed(jc,jk,jb) = sqrt( art_atmo%u(jc,jk,jb)*art_atmo%u(jc,jk,jb)   &
                   &                       + art_atmo%v(jc,jk,jb)*art_atmo%v(jc,jk,jb) )
             ELSE
                 windspeed(jc,jk,jb) = 0.
             END IF
           ENDDO
         ENDDO
       ENDDO 

      ! ----------------------------------
      ! --- Get Nash-Criteria 
      ! ----------------------------------

      DO jk=lev_top, lev_bot
        CALL get_nash_criteria( p_patch, equivLat_gridpoints(:,jk,:), pv(:,jk,:),                &
                              & windspeed(:,jk,:), equivLat_isolines(:,jk), nash_crit(:,jk,:),   &
                              & nash_equivLat(jk), min_equivLat, breakup_windspeed)
      END DO

      ! ----------------------------------
      ! --- Get Vortex
      ! ----------------------------------

      SELECT CASE ( select_hemisphere )
      CASE ( 'SH' )
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb,  i_startidx, i_endidx)
          DO jk=lev_top, lev_bot
            DO jc=i_startidx,i_endidx
              IF( equivLat_gridpoints(jc,jk,jb) <= nash_equivLat(jk) &
                & .AND. nash_equivLat(jk) .NE. 0.0_wp) THEN

                  ptr_vortex_log(jc,jk,jb) = .TRUE.
              END IF
            END DO
          END DO
        END DO

      CASE ( 'NH' )
        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb,  i_startidx, i_endidx)
          DO jk=lev_top, lev_bot
            DO jc=i_startidx,i_endidx
              IF( equivLat_gridpoints(jc,jk,jb) >= nash_equivLat(jk) &
                  & .AND. nash_equivLat(jk) .NE. 0.0_wp) THEN

                 ptr_vortex_log(jc,jk,jb) = .TRUE.
              END IF
            ENDDO
          ENDDO                                                                            
        ENDDO

      CASE DEFAULT
          ptr_vortex_log(:,:,:) = .FALSE.
      END SELECT

      ! ----------------------------------
      ! --- Deallocation
      ! ----------------------------------

      IF(ALLOCATED(pv) )                  DEALLOCATE( pv )
      IF(ALLOCATED(nash_equivLat) )       DEALLOCATE( nash_equivLat )
      IF(ALLOCATED(nash_crit) )           DEALLOCATE( nash_crit )
      IF(ALLOCATED(windspeed) )           DEALLOCATE( windspeed )
      IF(ALLOCATED(equivLat_isolines) )   DEALLOCATE( equivLat_isolines )
      IF(ALLOCATED(equivLat_gridpoints) ) DEALLOCATE( equivLat_gridpoints )

      END IF !vortex_init_date==act.date
    END DO !no_of_vortex_inits
  END IF !no_of_vortex_inits


END SUBROUTINE art_define_polar_vortex_boundary

SUBROUTINE get_nash_criteria( p_patch, equivLat_in, pv_in, windspeed_in, equivLat_isolines, &
                            & local_nash, equivLat_nash, min_equivLat, breakup_windspeed)
!<
! SUBROUTINE get_nash_criteria
! Routine to get points where Nash-Criteria holds
! Part of Module: mo_define_polar_vortex
! Author: Christian Stassen, KIT
! Initial Release: 2015-03-24
!>
  TYPE(t_patch), TARGET, INTENT(IN)  ::  &
    &  p_patch              !< patch on which computation is performed
  REAL(wp), INTENT(IN)               ::  &
    &  equivLat_in(:,:)     !< Equivalent Latitude (nproma, nblks)   
  REAL(wp), INTENT(IN)               ::  &
    &  pv_in(:,:)           !< Potential Vorticty (nproma, nblks)   
  REAL(wp), INTENT(IN)               ::  &
    &  windspeed_in(:,:)    !< Total windspeed (nproma, nblks)     
  REAL(wp), INTENT(INOUT)            ::  &
    &  local_nash(:,:)      !< Nash-Criteria (nproma, nblks)#
  REAL(wp), INTENT(INOUT)            ::  &
    &  equivLat_nash        !< Equivalent Latitude at max. Nash-Criteria
  REAL(wp), INTENT(IN)               ::  &
    &  equivLat_isolines(:) !< Isolines of equivalent latitude
  REAL(wp), INTENT(IN)               ::  &
    &  min_equivLat         !< Minimum of equivalent latitude at which the vortex
                            !  is located
  REAL(wp), INTENT(IN)               ::  &
    &  breakup_windspeed    !< minimum windspeed for the vortex

  ! --- Global variables
  REAL(wp), ALLOCATABLE              :: global_equivLat(:)
  REAL(wp), ALLOCATABLE              :: global_pv(:)
  REAL(wp), ALLOCATABLE              :: mean_global_windspeed(:) 
  REAL(wp), ALLOCATABLE              :: global_nash(:) 


  !Local variables
  INTEGER ::                     &
     &    jc, jb, i, j, jg,      &                   !< loop indizes
     &    i_startidx, i_endidx    

  ! -- Variables needed for communication
  INTEGER              :: work_comm_size             !< total number of working procs
  INTEGER              :: sendcount, &               !< number of items to send
                        & recvcount                  !< number of items to recieve
  INTEGER, ALLOCATABLE :: sendcounts(:), &           !< number of items to send to proc (i)
                        & recvcounts(:), &           !< number of items to recieve from proc (i)
                        & displ(:)                   !< number of displacements of proc (i)

  REAL(wp) :: totNoCells                             !< total number of cells

  INTEGER, ALLOCATABLE  :: global_field_index(:)     !< global index
  INTEGER, ALLOCATABLE  :: global_orderd_index(:)    !< global orderd index

  INTEGER               :: n_isolines
  INTEGER, ALLOCATABLE  :: size_global_windspeed(:)

  REAL(wp), ALLOCATABLE :: mean_local_windspeed(:)
  REAL(wp), ALLOCATABLE :: size_local_windspeed(:)

  REAL(wp)              :: mean_nash, size_nash, & 
                         & mean_equivLat, size_equivLat

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  ! ----------------------------------
  ! --- Start routine
  ! ----------------------------------

  ! ----------------------------------
  ! --- Get the loop indices
  ! ----------------------------------

  jg = p_patch%id
  art_atmo => p_art_data(jg)%atmo

  ! ----------------------------------
  ! --- Get field sizes
  ! ----------------------------------

  sendcount = 0
  totNoCells = 0.
  
   DO jb = art_atmo%i_startblk, art_atmo%i_endblk
     CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

     DO jc=i_startidx,i_endidx
       sendcount = sendcount + 1
     ENDDO
   ENDDO

  totNoCells = p_patch%n_patch_cells_g
  work_comm_size = get_my_mpi_work_comm_size()
  n_isolines = SIZE(equivLat_isolines(:), DIM=1)

  ! ----------------------------------
  ! --- Allocation
  ! ----------------------------------

  IF (my_process_is_mpi_workroot()) THEN
    ALLOCATE( global_equivLat(INT(totNoCells)),             &
            & global_pv(INT(totNoCells)),                   &
            & global_nash(INT(totNoCells)),                 &
            & global_field_index(INT(totNoCells)),          &
            & global_orderd_index(INT(totNoCells))          )

  END IF

  ALLOCATE( mean_local_windspeed(n_isolines), &
          & size_local_windspeed(n_isolines)  )

  ALLOCATE( mean_global_windspeed(n_isolines), &
          & size_global_windspeed(n_isolines)  )

  ALLOCATE( sendcounts(0:work_comm_size-1) , &
          & recvcounts(0:work_comm_size-1) , &
          & displ(0:work_comm_size-1)        ) 

  ! ----------------------------------
  ! --- Get mean windspeed along
  ! --- the equivalent latitudes
  ! ----------------------------------

  mean_local_windspeed(:) = 0.
  size_local_windspeed(:) = 0

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
    DO jc = i_startidx,i_endidx
      DO i = 1,n_isolines      !in der i Schleife tritt beim 2. Aufrauf ein Fehler auf!
        IF( equivLat_in(jc,jb) .EQ. equivLat_isolines(i) ) THEN 
          mean_local_windspeed(i) = mean_local_windspeed(i) + windspeed_in(jc,jb)
          size_local_windspeed(i) = size_local_windspeed(i) + 1.0_wp
        END IF
      ENDDO
    ENDDO
  ENDDO

  DO i = 1,n_isolines
    mean_global_windspeed(i) = global_sum_array(mean_local_windspeed(i))
    size_global_windspeed(i) = global_sum_array(size_local_windspeed(i))

    IF (size_global_windspeed(i) .NE. 0.0_wp) THEN
      mean_global_windspeed(i) = mean_global_windspeed(i) / size_global_windspeed(i)
    ELSE
      mean_global_windspeed(i) = 0.0_wp
    END IF
  END DO        

  ! ----------------------------------
  ! --- Gather field
  ! ----------------------------------

  IF (my_process_is_work()) THEN
    CALL p_send( sendcount, get_mpi_all_workroot_id(), &
               & get_my_mpi_work_id(), 1, p_comm_work)
    recvcount = sendcount   !Send and recieve are equal 
  END IF

  IF (my_process_is_mpi_workroot()) THEN
    DO i=0,work_comm_size-1
      CALL p_recv(recvcounts(i), i, i, 1, p_comm_work)

    END DO
    sendcounts(:) = recvcounts(:)   !Send and recieve are equal 
  END IF

  IF (my_process_is_mpi_workroot()) THEN
    DO i=0,work_comm_size-1
        displ(i) = i*recvcounts(i)
    END DO
  END IF

  CALL p_bcast(displ(:), get_mpi_all_workroot_id(), p_comm_work)
  CALL p_bcast(recvcounts(:), get_mpi_all_workroot_id(), p_comm_work)

  CALL p_gatherv( equivLat_in(:,:), sendcount,           &
                & global_equivLat(:), recvcounts(:),     &
                & displ(:), get_mpi_all_workroot_id(), p_comm_work )

  CALL p_gatherv( pv_in(:,:), sendcount,            &
                & global_pv(:), recvcounts(:),      &
                & displ(:), get_mpi_all_workroot_id(), p_comm_work )

  ! ----------------------------------
  ! --- Get "global" index to be sorted
  ! ----------------------------------

  IF (my_process_is_mpi_workroot()) THEN
    DO i = 1,INT(totNoCells)
      global_field_index(i) = i
    END DO
  END IF

  ! ----------------------------------
  ! --- Sort global field
  ! ----------------------------------

  IF (my_process_is_mpi_workroot()) THEN
    global_orderd_index(:) = global_field_index(:)
    CALL quicksort(global_equivLat(:), global_orderd_index(:))
  END IF
  

  ! ----------------------------------
  ! --- Compute Nash Criteria
  ! ----------------------------------

  IF (my_process_is_mpi_workroot()) THEN
    global_nash(:) = 0.
    DO i = 2,INT(totNoCells)-1
      DO j = 1,n_isolines
        IF ( mean_global_windspeed(j) >= breakup_windspeed       &
          & .AND. global_equivLat(i) .EQ. equivLat_isolines(j)   &
          & .AND. ABS(global_equivLat(i)) >= ABS(min_equivLat)   &
          & .AND. (  global_equivLat(global_orderd_index(i-1))   &
          &        - global_equivLat(global_orderd_index(i+1))) .NE. 0.0_wp ) THEN        

          global_nash(global_orderd_index(i)) =                 &
            &  ABS((global_pv(global_orderd_index(i-1))         &
            &       - global_pv(global_orderd_index(i+1)))      &
            & / (global_equivLat(global_orderd_index(i-1))      &
            &    - global_equivLat(global_orderd_index(i+1)))   &
            & * (mean_global_windspeed(j)))

        END IF
      END DO
    END DO
  END IF


  ! ----------------------------------
  ! --- Send back
  ! ----------------------------------

  CALL p_scatterv( global_nash(:), sendcounts(:),             &
                 & displ(:), local_nash(:,:), recvcount,      &
                 & get_mpi_all_workroot_id(), p_comm_work     )


  ! ----------------------------------
  ! --- Get mean nash
  ! ----------------------------------

  mean_nash = 0.
  size_nash = 0.
  mean_equivLat = 0.
  size_equivLat = 0.

   DO jb = art_atmo%i_startblk, art_atmo%i_endblk
     CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

     DO jc = i_startidx,i_endidx
       IF (local_nash(jc,jb) > 0.) THEN
         mean_nash = mean_nash + local_nash(jc,jb)
         size_nash = size_nash + 1.0_wp
       END IF
     ENDDO
   ENDDO

  IF (global_sum_array(size_nash) > 0.0_wp) THEN
    mean_nash = global_sum_array(mean_nash) / global_sum_array(size_nash)
  END IF

  ! ----------------------------------
  ! --- Get mean of equivLat
  ! ----------------------------------

  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

    DO jc = i_startidx,i_endidx
      IF (local_nash(jc,jb) > mean_nash) THEN
        mean_equivLat = mean_equivLat + equivLat_in(jc,jb)
        size_equivLat = size_equivLat + 1
      END IF
    END DO
  END DO

  IF (global_sum_array(size_equivLat) > 0.) THEN
    equivLat_nash = global_sum_array(mean_equivLat) / global_sum_array(size_equivLat)
  ELSE
    equivLat_nash = 0.
  END IF

  ! ----------------------------------
  ! --- DEALLOCATION
  ! ----------------------------------

  IF (my_process_is_mpi_workroot()) THEN

    IF ( ALLOCATED(global_equivLat) ) DEALLOCATE( global_equivLat )
    IF ( ALLOCATED(global_pv) ) DEALLOCATE( global_pv ) 
    IF ( ALLOCATED(global_nash) ) DEALLOCATE( global_nash) 
    IF ( ALLOCATED(global_field_index) ) DEALLOCATE( global_field_index )
    IF ( ALLOCATED(global_orderd_index) ) DEALLOCATE( global_orderd_index )

  END IF

END SUBROUTINE get_nash_criteria
END MODULE mo_art_define_polar_vortex

