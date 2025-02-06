!
! mo_art_pntSrc_types
! This module provides data structures required by the pntSrc component
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

MODULE mo_art_pntSrc_types
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_math_constants,                ONLY: deg2rad
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_expression,                    ONLY: expression
  USE mo_fortran_tools,                 ONLY: init
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_parallel_config,               ONLY: nproma, p_test_run
  USE mo_gnat_gridsearch,               ONLY: gnat_init_grid, gnat_destroy, t_gnat_tree, &
                                          &   gnat_query_containing_triangles,           &
                                          &   gnat_merge_distributed_queries, gk
  USE mo_grid_config,                   ONLY: grid_sphere_radius
  USE mtime,                            ONLY: datetime, newDatetime, datetimeToString,   &
                                          &   max_datetime_str_len, newEvent, timedelta, &
                                          &   newTimedelta, ASSIGNMENT(=), OPERATOR(+),  &
                                          &   deallocateTimedelta, deallocateDatetime,   &
                                          &   event, isCurrentEventActive
! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN, UNDEF_REAL_ART, UNDEF_INT_ART

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_pntSrc_types'

  TYPE t_art_pntSrc
    CHARACTER(LEN=IART_VARNAMELEN) :: &
      &  id                       !< Name of source
    REAL(wp)                 :: &
      &  lon,                   & !< longitude of the emission source
      &  lat,                   & !< latitude of the emission source
      &  height,                & !< geometric height of the emission source in [m]
      &  height_bot,            & !< geometric bottom height of the emission source in [m]
      &  source_strength,       & !< Source strength of tracer per s-1, already converted into correct (tracer type specific) units
      &  emiss_rate0              !< for aerosol tracer (modal distribution) emission -> emission rate of number concentration
     REAL(wp),POINTER         :: &
      &  height_factor(:)         !< in case of emission profile, this height factor distributes 
                                  !  the normalized emission to model levels
    INTEGER                  :: &
      &  itr,                   & !< Index of corresponding tracer in tracer container
      &  itr0,                  & !< in case of aerosol tracer -> index of corresponding nmb_conc
      &  ithis_nlocal_pts,      & !< local number of sources on 'subdomain' in parallel mode  
      &  tri_iidx_loc,          & !< location idx of source
      &  tri_iblk_loc,          & !< location blk of source
      &  k_index,               & !< (Top)  vertical level of source
      &  k_index_bot              !< Bottom vertical level of source (uniform profile if defined)
    TYPE(datetime), POINTER, PRIVATE :: &
      &  start_time,            & !< Start time of the emission event
      &  end_time                 !< End time of the emission event
    TYPE(event), POINTER, PRIVATE    :: &
      &  emissEvent               !< Emission event
    CONTAINS
      PROCEDURE, PUBLIC :: init     => init_pntSrc
      PROCEDURE, PUBLIC :: print    => print_pntSrc_summary
      PROCEDURE, PUBLIC :: isActive => is_pntSrc_active
  END TYPE t_art_pntSrc

  TYPE t_art_all_pntSrc
    TYPE(t_art_pntSrc), POINTER :: &
      &  p(:)                        !< List of all point sources
    INTEGER                     :: &
      &  nsources                    !< number of active source scenarios
  END TYPE t_art_all_pntSrc

  PUBLIC :: t_art_all_pntSrc, t_art_pntSrc

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE init_pntSrc(this_pntSrc, tc_dt_model, tc_exp_refdate, p_patch, z_ifc, id, lon, lat,    &
  &                    height, itr, source_strength, startTime, endTime, lexclude_end,            &
  &                    emiss_profile, height_bot, itr0, emiss_rate0 )
!<
! SUBROUTINE init_pntSrc
! This routine initializes one point source based on a given patch,
! longitude, latitude, height, start, endtime and optional bottom height
! Based on: Werchner (2016) - Bachelorthesis, KIT
! Part of Module: mo_art_pntSrc_types
! Author: Sven Werchner, Daniel Rieger, KIT
! Initial Release: 2017-01-25
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
! 2019-08-07: Lukas Muser, KIT
! - expression for emission profile can be handed in which is evaluated during init routine
!
!>
  CLASS(t_art_pntSrc), INTENT(inout)    :: &
    &  this_pntSrc                           !< Point source to be initialized
  TYPE(timedelta), POINTER, INTENT(in)  :: &
    &  tc_dt_model                           !< Model timestep
  TYPE(datetime), POINTER, INTENT(in)   :: &
    &  tc_exp_refdate                        !< Experiment reference date
  TYPE(t_patch), INTENT(in)             :: &
    &  p_patch
  REAL(wp), INTENT(in)                  :: &
    &  z_ifc(:,:,:),                       & !< Height of half levels
    &  lon,                                & !< longitude of the emission source
    &  lat,                                & !< latitude of the emission source
    &  height,                             & !< height of the emission source above ground (m)
    &  source_strength                       !< Source strength of tracer per s-1
  INTEGER, INTENT(in)                   :: &
    &  itr                                   !< Index of tracer to apply this source scenario
  CHARACTER(LEN=*), INTENT(in)          :: &
    &  id                                    !< Name of source
  CHARACTER(LEN=*), INTENT(in)          :: &
    &  startTime, endTime                    !< Start and end time of source scenario, format: 'YYYY-MM-DDTHH:MM:SS'
  LOGICAL, INTENT(in), OPTIONAL         :: &
    &  lexclude_end                          !< Exclude endTime from active time interval of point source
  CHARACTER(LEN=*),INTENT(in), OPTIONAL :: &
    &  emiss_profile                         !< arithmetic expression of emission profile
  REAL(wp), INTENT(in), OPTIONAL        :: &
    &  height_bot,                         & !< bottom height of emission source with uniform profile above ground (m)
    &  emiss_rate0                           !< 
  INTEGER, INTENT(in), OPTIONAL         :: &
    &  itr0                                  !< index of corresponding nmb_conc for aero tracer
! Local variables
  TYPE(t_gnat_tree)         :: &
    &  gnat                      !< Grid search object
  TYPE(expression)          :: &
    &  emiss_expr                !< expression of emission profile
  TYPE(timedelta), POINTER  :: &
    &  td_minus_second           !< time delta of minus one second
  INTEGER                   :: &
    &  nblks, npromz,          & !< calculated for call to GNAT
    &  nlocsrc,                & !< counter for number of local sources
    &  gnat_jc,gnat_jb,        & !< coordinates for GNAT
    &  jc, jb, jk,             & !< Actual horizontal/vertical indices
    &  ierror
  REAL(wp), POINTER         :: &
    &  val_2d(:,:),            & !< value of evaluated arithmetic expression
    &  z_star(:,:)               !< z_ifc in 2D, z_star(nlev+1,1)
  INTEGER, ALLOCATABLE      :: &
    &  tri_idx(:,:,:),         &
    &  owner(:)                  !< rank of sender PE for each station (GNAT)
  REAL(gk), ALLOCATABLE     :: &
    &  in_points(:,:,:),       & !< geographical locations (GNAT)
    &  min_dist(:,:),          & !< minimal distance (GNAT)
    &  lat_idx(:,:),           &
    &  lon_idx(:,:)

  ierror = SUCCESS
  
  ! ----------------------------------
  ! --- Get the source position
  ! ----------------------------------

  nblks=1/nproma+1
  npromz=1-nproma*(nblks-1)
  
  ! --- Allocate strcutures demanded by GNAT
  ALLOCATE(tri_idx(2,nproma,nblks))
  ALLOCATE(owner(1))
  ALLOCATE(in_points(nproma,nblks,2))
  ALLOCATE(min_dist(nproma,nblks))
  ALLOCATE(lat_idx(nproma,nblks))
  ALLOCATE(lon_idx(nproma,nblks))
  
  ! --- Save lon/lat into structures demanded by GNAT
  gnat_jb=1
  gnat_jc=1

  lat_idx(gnat_jc,gnat_jb) = lat
  lon_idx(gnat_jc,gnat_jb) = lon

  in_points(gnat_jc,gnat_jb,1) = lon_idx(gnat_jc,gnat_jb) * deg2rad
  in_points(gnat_jc,gnat_jb,2) = lat_idx(gnat_jc,gnat_jb) * deg2rad

  ! --- Build GNAT data structure
  CALL gnat_init_grid(gnat, p_patch)

  ! --- Perform proximity query
  CALL gnat_query_containing_triangles(gnat, p_patch, in_points(:,:,:),              &
    &                                  nproma, nblks, npromz,                        &
    &                                  grid_sphere_radius,p_test_run,tri_idx(:,:,:), &
    &                                  min_dist(:,:))

  CALL gnat_merge_distributed_queries(p_patch, 1, nproma, nblks, min_dist,  &
    &                                 tri_idx(:,:,:), in_points(:,:,:),               &
    &                                 owner(:), this_pntSrc%ithis_nlocal_pts)

  ! --- Cleanup GNAT
  CALL gnat_destroy(gnat)

  ! --- Allocate metadata fields and dummy variables
  ALLOCATE(this_pntSrc%height_factor(p_patch%nlev))
  ALLOCATE(val_2d(p_patch%nlev+1,1))
  ALLOCATE(z_star(p_patch%nlev+1,1))
  CALL init(this_pntSrc%height_factor, lacc=.FALSE.)

  ! --- Save locations into this_pntSrc
  gnat_jc=0
  gnat_jb=1

  DO nlocsrc = 1, this_pntSrc%ithis_nlocal_pts
    gnat_jc=gnat_jc+1
    IF(gnat_jc>nproma) THEN
      gnat_jc=1
      gnat_jb=gnat_jb+1
    ENDIF
    
    this_pntSrc%tri_iidx_loc=tri_idx(1,gnat_jc,gnat_jb)
    this_pntSrc%tri_iblk_loc=tri_idx(2,gnat_jc,gnat_jb)

    jc = this_pntSrc%tri_iidx_loc
    jb = this_pntSrc%tri_iblk_loc

 
    ! There are 3 different ways of setting the emission height in the pntSrc.xml
    !
    ! I)   Set (only) a positive value for height
    !      -> emission in a single level at the height above ground specified in the xml or in
    !         combination with emiss_profile this determines the max height of the profile
    ! II)  Set only a negative value for height
    !      -> emission between ABS(height) and ground level, either uniform profile or profile is
    !         given as emiss_profile
    ! III) Set positive value for height and also height_bot
    !      -> emission between height and height_bot, either uniform profile or profile is given as
    !         emiss_profile

    ! ----------------------------------
    ! --- Emission height handling
    ! ----------------------------------

    this_pntSrc%height_bot  = UNDEF_REAL_ART
    this_pntSrc%k_index_bot = UNDEF_INT_ART

    IF (height >= 0._wp) THEN
      this_pntSrc%height      = height + z_ifc(jc,p_patch%nlev+1,jb)      ! Adding surface height
      IF (PRESENT(height_bot)) THEN
        IF (height_bot /= UNDEF_REAL_ART) THEN

          ! emission with (uniform) profile from height_bot to height above surface
          this_pntSrc%height_bot = height_bot + z_ifc(jc,p_patch%nlev+1,jb) ! Adding surface height

          ! Find index of layer that contains the release height_bot
          DO jk= 1, p_patch%nlev
            IF( this_pntSrc%height_bot >= z_ifc(jc, jk+1 ,jb) .AND. &
              & this_pntSrc%height_bot <  z_ifc(jc, jk ,  jb)) THEN
              this_pntSrc%k_index_bot = jk
              EXIT
            ELSE
              ! value given for height_bot might be out of range, emission is set to single level
              ! emission
            ENDIF
          ENDDO

        ENDIF
      ELSE ! no height_bot given
        ! emission in a single model level
        ! this_pntSrc%height_bot and this_pntSrc%k_index_bot stay UNDEF
      ENDIF

    ELSE ! Negative height means (uniform) profile from surface up to given ABS(height)
      this_pntSrc%height      = -1._wp * height + z_ifc(jc,p_patch%nlev+1,jb) ! Adding surface height
      this_pntSrc%height_bot  = z_ifc(jc,p_patch%nlev+1,jb)                   ! Start at surface height
      this_pntSrc%k_index_bot = p_patch%nlev
    ENDIF ! height >= 0

    ! Find index of layer that contains the release height
    DO jk= 1, p_patch%nlev
      IF( this_pntSrc%height >= z_ifc(jc, jk+1 ,jb) .AND. &
        & this_pntSrc%height <  z_ifc(jc, jk ,  jb)) THEN
        this_pntSrc%k_index = jk
        EXIT
      ENDIF
    ENDDO

    ! For emission in single level set k_index_bot = k_index_top and set height factor 1.0 in
    ! respective level
    IF (this_pntSrc%k_index_bot == UNDEF_INT_ART) THEN
      this_pntSrc%k_index_bot = this_pntSrc%k_index
      this_pntSrc%height_factor(this_pntSrc%k_index) = 1.0_wp
    ENDIF
    
    ! ----------------------------------
    ! --- Emission profile handling
    ! ----------------------------------

    ! For a given emission profile, evaluate profile and determine height_factor
    SELECT CASE(emiss_profile)
      CASE('')
        ! IF height_bot OR height < 0 is given, but no emiss_profile specified, then assume a 
        ! constant emission profile between height_bot and height with 
        ! height_factor(jk) = dz(jk) / column_height
        IF (this_pntSrc%height_bot /= UNDEF_REAL_ART) THEN
          DO jk= this_pntSrc%k_index, this_pntSrc%k_index_bot
            this_pntSrc%height_factor(jk) = ( z_ifc(jc,jk,jb) - z_ifc(jc,jk+1,jb)  )              &
              &                           / ( z_ifc(jc,this_pntSrc%k_index,jb)                    &
              &                           -   z_ifc(jc,this_pntSrc%k_index_bot+1,jb) )
          ENDDO
        ENDIF
        ! For emission in single level, nothing to do

      CASE DEFAULT
        ! Emission profile is given in pntSrc.xml and will be evaluated. The result is stored as
        ! the height factor in the pntSrc meta data
        
        ! Assume that the emission profile reaches from ground to pntSrc height if it is not
        ! specified differently
        IF (this_pntSrc%height_bot == UNDEF_REAL_ART) THEN
          PRINT *,'ART: WARNING: emission profile is used without any bottom boundary, '//        &
            &     'surface level is used as default!'
          this_pntSrc%height_bot  = z_ifc(jc,p_patch%nlev+1,jb)
          this_pntSrc%k_index_bot = p_patch%nlev
        ENDIF

        emiss_expr = expression(emiss_profile)
        IF (emiss_expr%err_no /= SUCCESS) CALL finish(TRIM(routine)//':init_pntSrc',              &
                                            &         'Invalid expression for emiss_profile: '    &
                                            &         //emiss_profile//'.')
        ! link variables to expression
        z_star(:,1) = ( z_ifc(jc,:,jb) - this_pntSrc%height_bot )                                 &
          &         / ( this_pntSrc%height - this_pntSrc%height_bot )
        z_star(this_pntSrc%k_index_bot+1,1) = 0._wp
        z_star(this_pntSrc%k_index,1)       = 1._wp
        CALL emiss_expr%link("z_star", z_star)
        
        ! evaluate expression
        CALL emiss_expr%evaluate(val_2d)
        IF (emiss_expr%err_no /= SUCCESS) CALL finish(TRIM(routine)//':init_pntSrc', 'Evaluation  &
                                            &         of emiss_profile: '//emiss_profile//        &
                                            &         ' failed.')
        
        DO jk= this_pntSrc%k_index, this_pntSrc%k_index_bot
          ! height factor for level jk is the emission profile evaluated at the upper boundary
          ! minus the lower boundary of the model level
          this_pntSrc%height_factor(jk) = val_2d(jk,1) - val_2d(jk+1,1)
        ENDDO
        ! Print warning if SUM(height_factor) /= 1
        IF (ABS(SUM(this_pntSrc%height_factor) - 1._wp) > 1.0e-10_wp) THEN
          PRINT *,'ART: WARNING: '//TRIM(routine)//':init_pntSrc: sum of height factor is not '// &
            &     'equal to 1. Hence, the source strength is scaled by the factor = ',            &
            &      SUM(this_pntSrc%height_factor) 
        ENDIF
        
        ! Clean up expression
        CALL emiss_expr%finalize()
    END SELECT ! emiss_profile

  ENDDO ! nlocsrc

  this_pntSrc%id              = TRIM(ADJUSTL(id))
  this_pntSrc%lon             = lon
  this_pntSrc%lat             = lat
  this_pntSrc%itr             = itr
  this_pntSrc%source_strength = source_strength

  IF (PRESENT(itr0)) THEN
    this_pntSrc%itr0 = itr0
  ELSE
    this_pntSrc%itr0 = UNDEF_INT_ART
  ENDIF

  IF (PRESENT(emiss_rate0)) THEN
    this_pntSrc%emiss_rate0 = emiss_rate0
  ELSE
    this_pntSrc%emiss_rate0 = UNDEF_REAL_ART
  ENDIF

  ! ----------------------------------
  ! --- Time handling
  ! ----------------------------------
  this_pntSrc%start_time => newDatetime(startTime, errno=ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':init_pntSrc', &
                            &         'Could not create datetime object from string '//TRIM(startTime)//'.')
  this_pntSrc%end_time => newDatetime(endTime, errno=ierror)
    IF(ierror /= SUCCESS) CALL finish(TRIM(routine)//':init_pntSrc', &
                            &         'Could not create datetime object from string '//TRIM(endTime)//'.')
  IF (PRESENT(lexclude_end)) THEN
    IF (lexclude_end) THEN
      td_minus_second => newTimedelta('-',0,0,0,0,0,second=1,ms=0)
      this_pntSrc%end_time = this_pntSrc%end_time + td_minus_second
      CALL deallocateTimedelta(td_minus_second)
    ENDIF
  ENDIF

  this_pntSrc%emissEvent => newEvent('Emission',               &
    &                                tc_exp_refdate,           &
    &                                this_pntSrc%start_time,   &
    &                                this_pntSrc%end_time,     &
    &                                tc_dt_model)

  ! ----------------------------------
  ! --- Clean up
  ! ----------------------------------
  DEALLOCATE(val_2d)
  DEALLOCATE(z_star)
  DEALLOCATE(tri_idx)
  DEALLOCATE(owner)
  DEALLOCATE(in_points)
  DEALLOCATE(min_dist)
  DEALLOCATE(lat_idx)
  DEALLOCATE(lon_idx)

END SUBROUTINE init_pntSrc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE print_pntSrc_summary(this_pntSrc)
!<
! SUBROUTINE print_pntSrc_summary
! Print all information of current source to standard output
! Based on: -
! Part of Module: mo_art_pntSrc_types
! Author: Daniel Rieger, KIT
! Initial Release: 2017-01-26
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CLASS(t_art_pntSrc),INTENT(inout)   :: &
    &  this_pntSrc                         !< Point source
! Local variables
  CHARACTER(LEN=max_datetime_str_len) :: &
    &  dstring                             !< Datetime as string

  WRITE (message_text,*) '==ART: PRINTOUT INFO PNTSRC '//TRIM(ADJUSTL(this_pntSrc%id))//'=='
  CALL message ('', message_text)
  WRITE (message_text,*) 'NAME            DATA'
  CALL message ('', message_text)
  WRITE (message_text,'(A16,E13.6)') 'lon:            ',this_pntSrc%lon
  CALL message ('', message_text)
  WRITE (message_text,'(A16,E13.6)') 'lat:            ',this_pntSrc%lat
  CALL message ('', message_text)
  WRITE (message_text,'(A16,E13.6)') 'source strength:',this_pntSrc%source_strength
  CALL message ('', message_text)
  WRITE (message_text,'(A16,I3)')    'itr:            ',this_pntSrc%itr
  CALL message ('', message_text)

  ! mtime dates
  CALL datetimeToString(this_pntSrc%start_time, dstring)
  WRITE (message_text,*)    'Start time:     '//TRIM(ADJUSTL(dstring))
  CALL message ('', message_text)
  CALL datetimeToString(this_pntSrc%end_time, dstring)
  WRITE (message_text,*)    'End time:       '//TRIM(ADJUSTL(dstring))
  CALL message ('', message_text)

END SUBROUTINE print_pntSrc_summary
!!
!!-------------------------------------------------------------------------
!!
LOGICAL FUNCTION is_pntSrc_active(this_pntSrc, current_date)
!<
! FUNCTION is_pntSrc_active
! Returns if the point source is currently active
! Based on: -
! Part of Module: mo_art_pntSrc_types
! Author: Daniel Rieger, KIT
! Initial Release: 2017-02-20
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CLASS(t_art_pntSrc),INTENT(inout)   :: &
    &  this_pntSrc                         !< Point source
  TYPE(datetime), POINTER, INTENT(in) :: &
    &  current_date                        !< Current date
! Local variables
  
  IF (isCurrentEventActive(this_pntSrc%emissEvent, current_date)) THEN
    is_pntSrc_active = .TRUE.
  ELSE
    is_pntSrc_active = .FALSE.
  ENDIF
END FUNCTION is_pntSrc_active
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_pntSrc_types
