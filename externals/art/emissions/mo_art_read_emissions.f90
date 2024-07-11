
!
! mo_art_read_emissions
! This module provides routines for automatically find the emission files with names
! according to the ICON-ART file name convention (see mo_art_create_filenames for details)
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

MODULE mo_art_read_emissions
! ICON
  USE mo_kind,                    ONLY: wp, i8
  USE mo_physical_constants,      ONLY: argas, amd
  USE mo_run_config,              ONLY: ntracer
  ! This is needed because of the call of openInputFile. If this is in a wrapper
  ! routine (eventually), this could be deleted as well
  USE mo_model_domain,            ONLY: t_patch
  USE mo_exception,               ONLY: finish, message

  USE mo_art_config,              ONLY: IART_PATH_LEN, art_config
  USE mo_read_interface,          ONLY: t_stream_id, openInputFile, closeFile,  &
    &                               read_2d_1time, read_2d_1lev_1time,          &
    &                               read_3d_1time, on_cells
  USE mo_io_config,               ONLY: default_read_method
! EXTERNALS
  USE mtime,                      ONLY: OPERATOR(>=), datetime,timedelta,      &
                                    &   OPERATOR(+),  OPERATOR(-),             &
                                    &   OPERATOR(==),                          &
                                    &   newDatetime, OPERATOR(>),OPERATOR(<),  &
                                    &   max_datetime_str_len, newTimeDelta,    &
                                    &   datetimeToPosixString,                 &
                                    &   divisionquotienttimespan,              &
                                    &   divideTwoDatetimeDiffsInSeconds,       &
                                    &   deallocateTimeDelta,                   &
                                    &   deallocateDatetime,                    &
                                    &   getNoOfDaysInYearDateTime
! ART
  USE mo_art_emiss_types,         ONLY: t_art_emiss_storage,          &
                                   &    t_art_emiss_type_container
  USE mo_art_prescribed_types,    ONLY: t_art_emiss_prescribed, &
                                   &    t_art_chem_prescribed,  &
                                   &    t_art_prescribed
  USE mo_art_create_filenames,    ONLY: art_create_filenames
  USE mo_art_bvoc,                ONLY: bvoc_guenther2012, sel_bccnum

  USE mo_art_data,                ONLY: p_art_data
  USE mo_art_atmo_data,           ONLY: t_art_atmo
  USE mo_art_wrapper_routines,    ONLY: art_get_indices_c
  USE mo_art_vinterp,             ONLY: art_prepare_vinterp,     &
                                    &   art_prepare_vinterp_pres
  USE mo_art_vinterp_chem_init,   ONLY: art_vinterp_chem_init
  
  


  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: art_read_emissions, art_init_emission_struct
  PUBLIC :: art_add_emission_to_tracers
  PUBLIC :: art_convert_emission_to_mmr

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init_emission_struct(p_patch, mtime_current,  emiss_struct)
!<
! SUBROUTINE art_init_emission_struct                   
! This subroutine sets the boundary dates of the respective emission dataset,
! and searches and reads the emission files closest to the simulation time.
! Finally, the current emission is calculated by linear interpolation
! Part of Module: mo_art_read_emission
! Author: Michael Weimer, KIT
! Initial Release: 2016-03-08              
! Modifications:
!>
  TYPE(t_patch), TARGET, INTENT(IN)   :: &
    &  p_patch                             !< patch on which computation is performed
  
  TYPE(datetime), POINTER :: &
    &  mtime_current                       !< current simulation time

  CLASS(t_art_prescribed), INTENT(INOUT), TARGET :: &
    &  emiss_struct                        !< structure including all information about
                                           !  the current emission dataset which is read

  ! Local variables
  TYPE(datetime), POINTER ::             &
    &  mtime_current_0minsec,            & !< mtime_current with minute and second set to 0
    &  datetime_before, datetime_after,  & !< pointer to datetime before, after and
    &  first_date, last_date                              !< pointer to first and last_date in 
                                                          !  emiss_struct (read from xml file)
  
  REAL(wp), POINTER ::         &
    &  this_emiss_2d(:,:),     &    !< pointer to emission data at current,
    &  emiss_2d_before(:,:),   &    !                           before
    &  emiss_2d_after(:,:),    &    !                           and after emission dataset
    &  this_emiss_3d(:,:,:,:),   &    !< pointer to emission data at current,
    &  emiss_3d_before(:,:,:,:), &    !                           before
    &  emiss_3d_after(:,:,:,:)        !                           and after emission dataset

  TYPE(divisionquotienttimespan) ::  &
    &  quotient                     !< quotient and remainder of division of 
                                    !  two datetime differences
  INTEGER(i8)  ::    &
    &  denominator                  !< denominator of the above fraction
  INTEGER ::      &
    &  ierr,      &                 !< error handler of creating datetime
    &  jg
  REAL(wp) ::  &
    &  x                            !< time fraction (t - t1)/(t2 - t1) of linear interpolation
  TYPE(timedelta), POINTER :: &
    &  td_one_hour => NULL()        !< timedelta of one hour (or minus hour)
  INTEGER(i8) :: &
    &  before_orig_year,  &
    &  after_orig_year

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  jg       =  p_patch%id
  art_atmo => p_art_data(jg)%atmo

  ! --------------------------------------------------------------------------------
  ! 1. Allocations and definitions: Initialization of pointers
  ! --------------------------------------------------------------------------------

  ! The usage of the original mtime_current with minutes and seconds leads to
  ! problems since emissions can be only hourly in the maximum. So set minute
  ! and second to 0, but only in initialisation. This essentially means that in
  ! the first time step, if the current time is between one hour around an
  ! emission dataset, the linear interpolation coefficient is set to 0 or 1.
  mtime_current_0minsec => newDatetime(mtime_current)
  mtime_current_0minsec%time%minute = 0
  mtime_current_0minsec%time%second = 0
  mtime_current_0minsec%time%ms = 0

  first_date        => emiss_struct%first_date
  last_date         => emiss_struct%last_date
  SELECT TYPE (presc => emiss_struct)
    TYPE IS (t_art_emiss_prescribed)
      ALLOCATE(presc%emiss_2d_read(art_atmo%nproma,art_atmo%nblks))
      ALLOCATE(presc%emiss_2d(art_atmo%nproma,art_atmo%nblks,3))
      presc%emiss_2d(:,:,:) = 0._wp

      emiss_2d_before   => presc%emiss_2d(:,:,1)
      this_emiss_2d     => presc%emiss_2d(:,:,2)
      emiss_2d_after    => presc%emiss_2d(:,:,3)

    TYPE IS (t_art_chem_prescribed)
      ALLOCATE(presc%vinterp_3d(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks,presc%num_vars,3))

      IF ((TRIM(presc%unit_vertical) == '') .AND.   &
         &     (TRIM(presc%vert_coord_filename) == '')) THEN

        presc%prescribed_state%chem_init_in%nlev_chem_init = 1

        ALLOCATE(presc%prescribed_state%chem_init_chem%spec(art_atmo%nproma,1,  &
                             &                        art_atmo%nblks,presc%num_vars))
        presc%prescribed_state%chem_init_in%linitialized = .TRUE.
      END IF

      presc%vinterp_3d(:,:,:,:,:) = 0._wp

      emiss_3d_before   => presc%vinterp_3d(:,:,:,:,1)
      this_emiss_3d     => presc%vinterp_3d(:,:,:,:,2)
      emiss_3d_after    => presc%vinterp_3d(:,:,:,:,3)

  END SELECT
    
  emiss_struct%datetime(1)%ptr => newDatetime(mtime_current_0minsec, ierr)
  emiss_struct%datetime(3)%ptr => newDatetime(mtime_current_0minsec, ierr)


  datetime_before   => emiss_struct%datetime(1)%ptr
  datetime_after    => emiss_struct%datetime(3)%ptr

  ! --------------------------------------------------------------------------------
  ! 2. On initialisation (in mo_art_init):
  !   -  set datetimes of external emission dataset correctly
  !       --> read first_date and last_date of dataset from file first_and_last_date.txt
  !       --> calculate datetimes of emission files that exist the nearest around 
  !           the current datetime (--> datetime_before and datetime_after) (see below for details)
  !   -  read the datasets which exist at datetime_before and datetime_after
  ! --------------------------------------------------------------------------------


  !! to ensure that art_emission_chemtracer reads the data at the first timestep, too:

  IF (last_date < first_date) THEN
    CALL finish('mo_art_read_emissions:art_init_emission_struct',     &
      &         'last date is earlier than first date.')
  END IF

  ! datetime_before and datetime_after outside reading have the same year as the current datetime,
  ! but if datetime is not between first and last date, years of datetime_before and datetime_after
  ! have to be set to a value that is between first and last date (only for reading the file).
  ! After reading the file, the years are reset to values around current datetime.

  before_orig_year = datetime_before%date%year
  after_orig_year = datetime_after%date%year

  IF (mtime_current_0minsec > last_date) THEN
    datetime_before%date%year = last_date%date%year
    datetime_after%date%year  = last_date%date%year

    IF (datetime_after > last_date) THEN
      datetime_after%date%year = datetime_after%date%year - 1
    END IF

    CALL art_determine_datetimes_before_or_after(p_patch, emiss_struct,1,  &
                               &                 before_orig_year - datetime_before%date%year)

  ELSE IF (first_date > mtime_current_0minsec) THEN

    datetime_before%date%year = first_date%date%year
    datetime_after%date%year  = first_date%date%year

    IF (datetime_before < first_date) THEN
      datetime_before%date%year = datetime_before%date%year + 1
    END IF

    CALL art_determine_datetimes_before_or_after(p_patch,emiss_struct,3,  &
                                 &               after_orig_year - datetime_after%date%year)
  END IF


  ! if file exists at current datetime (= start date or after call of 
  ! art_determine_datetimes_before_after)
  ! in art_search_read_file_around_datetime, first search is done at one hour 
  ! before or after the committed datetime --> one hour is added to datetime_before
  
  td_one_hour => newTimeDelta('+PT1H')
  datetime_before = datetime_before + td_one_hour
  CALL deallocateTimeDelta(td_one_hour)

  td_one_hour => newTimeDelta('-PT1H')
  datetime_after = datetime_after + td_one_hour
  CALL deallocateTimeDelta(td_one_hour)


  CALL art_search_read_file_around_datetime(p_patch,emiss_struct,-1,1,  &
                             &              before_orig_year - datetime_before%date%year)

  CALL art_search_read_file_around_datetime(p_patch,emiss_struct,+1,3,  &
                             &              after_orig_year - datetime_after%date%year )

  ! datetime_before and datetime_after are reset to original years 
  ! or one year around if year changed while searching
  datetime_before%date%year = mtime_current_0minsec%date%year 
  datetime_after%date%year = mtime_current_0minsec%date%year

  IF (datetime_before > mtime_current_0minsec) THEN
    datetime_before%date%year = datetime_before%date%year - 1
  END IF

  IF (mtime_current_0minsec > datetime_after) THEN
    datetime_after%date%year = datetime_after%date%year + 1
  END IF     

  IF (mtime_current_0minsec == datetime_before) THEN
    x = 0._wp
  ELSE IF (mtime_current_0minsec == datetime_after) THEN
    x = 1._wp
  ELSE
    ! calculation of time fraction for linear interpolation
    ! Please note: for calculation the original current time is used
    CALL divideTwoDatetimeDiffsInSeconds(mtime_current,datetime_before, &
                          &              datetime_after,datetime_before,        &
                          &              denominator,quotient)

    x = REAL(quotient%remainder_in_ms,wp) / 1000._wp / REAL(denominator,wp)
  END IF

  SELECT TYPE (presc => emiss_struct)
    TYPE IS (t_art_emiss_prescribed)
      this_emiss_2d = emiss_2d_before + x * (emiss_2d_after - emiss_2d_before)
    TYPE IS (t_art_chem_prescribed)
      this_emiss_3d = emiss_3d_before + x * (emiss_3d_after - emiss_3d_before)
  END SELECT

  emiss_struct%type_is_init = .TRUE.

  CALL deallocateDatetime(mtime_current_0minsec)

END SUBROUTINE art_init_emission_struct
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_emissions(p_patch, mtime_current,  emiss_struct)
!<
! SUBROUTINE art_read_emissions
! This subroutine searches and reads the emission files closest  to the simulation time 
! if this is necessary. Finally, the current emission is calculated by linear interpolation
! Part of Module: mo_art_read_emission
! Author: Michael Weimer, KIT
! Initial Release: about 2015-04-01
! Modifications:
! 2016-03-08: Michael Weimer, KIT
! - separated initial call from operational subroutine
!>
  TYPE(t_patch), TARGET, INTENT(IN)   :: &
    &  p_patch                             !< patch on which computation is performed
  
  TYPE(datetime), POINTER :: &
    &  mtime_current                       !< current simulation time

  CLASS(t_art_prescribed), INTENT(INOUT), TARGET :: &
    &  emiss_struct                        !< structure including all information about the current
                                           !  emission dataset which is read

  ! Local variables
  TYPE(datetime), POINTER ::                          &
    &  datetime_before, datetime_after, & !< pointer to datetime before and after emission dataset
    &  first_date, last_date              !< pointer to first and last_date in emiss_struct 
                                          !  (read from 'first_and_last_date.txt')
  
  REAL(wp), POINTER ::         &
    &  this_emiss_2d(:,:),     &    !< pointer to emission data at current,
    &  emiss_2d_before(:,:),   &    !                           before
    &  emiss_2d_after(:,:),    &    !                           and after emission dataset
    &  this_emiss_3d(:,:,:,:),   &    !< pointer to emission data at current,
    &  emiss_3d_before(:,:,:,:), &    !                           before
    &  emiss_3d_after(:,:,:,:)        !                           and after emission dataset

  INTEGER ::                  &
    &  i, jk,                 &     !< loop indices
    &  jg
  INTEGER ::  jb, istart, iend
  TYPE(divisionquotienttimespan) :: &
    &  quotient                     !< quotient and remainder of division
                                    !  of two datetime differences
  INTEGER(i8)  ::    &
    &  denominator                  !< denominator of the fraction above
  REAL(wp) ::  &
    &  x                            !< time fraction (t - t1)/(t2 - t1) of linear interpolation
  TYPE(timedelta), POINTER :: &
     & td_one_hour => NULL()        !< timedelta of one hour (or minus one hour)
  INTEGER(i8) :: &
    &  after_orig_year
  TYPE(t_art_atmo), POINTER :: &
    & art_atmo


  jg = p_patch%id
  art_atmo => p_art_data(jg)%atmo

  ! --------------------------------------------------------------------------------
  ! 1. Definitions: Initialization of pointers
  ! --------------------------------------------------------------------------------

  SELECT TYPE (presc => emiss_struct)
    TYPE IS (t_art_emiss_prescribed)
      emiss_2d_before   => presc%emiss_2d(:,:,1)
      this_emiss_2d     => presc%emiss_2d(:,:,2)
      emiss_2d_after    => presc%emiss_2d(:,:,3)
    TYPE IS (t_art_chem_prescribed)
      emiss_3d_before   => presc%vinterp_3d(:,:,:,:,1)
      this_emiss_3d     => presc%vinterp_3d(:,:,:,:,2)
      emiss_3d_after    => presc%vinterp_3d(:,:,:,:,3)
  END SELECT

  datetime_before   => emiss_struct%datetime(1)%ptr
  datetime_after    => emiss_struct%datetime(3)%ptr

  first_date        => emiss_struct%first_date
  last_date         => emiss_struct%last_date

  ! -----------------------------------------------------------------------------------------------
  ! 2. tasks during operation:
  !   - when current datetime exceeds datetime_after:
  !     --> datetime_before and emission can be copied from datetime and emission after the dataset
  !     --> search a new file after current datetime (with the correct year, similar to first call)
  !   - linear interpolation is calculated and saved into this_emission
  ! -----------------------------------------------------------------------------------------------

  after_orig_year = datetime_after%date%year

  IF  (mtime_current > datetime_after)  THEN
    datetime_before = datetime_after

    SELECT TYPE (presc => emiss_struct)
      TYPE IS (t_art_emiss_prescribed)
        emiss_2d_before = emiss_2d_after
      TYPE IS (t_art_chem_prescribed)
        emiss_3d_before = emiss_3d_after
    END SELECT

    IF (datetime_after >= last_date) THEN
      datetime_after%date%year = last_date%date%year

      IF (datetime_after >= last_date) THEN
        datetime_after%date%year = datetime_after%date%year - 1
      END IF

    ELSE IF (first_date > datetime_after) THEN
      datetime_after%date%year = first_date%date%year
    
      IF (datetime_after == last_date) THEN
        datetime_after%date%year = datetime_after%date%year - 1
      ELSE
        CALL art_determine_datetimes_before_or_after(p_patch,emiss_struct,3, &
                                  &                  after_orig_year - datetime_after%date%year)

        td_one_hour => newTimeDelta('-PT1H')
        datetime_after = datetime_after + td_one_hour
        CALL deallocateTimeDelta(td_one_hour)
      END IF

    END IF

    CALL art_search_read_file_around_datetime(p_patch,emiss_struct,1,3,  &
                       &                      after_orig_year - datetime_after%date%year)

    datetime_after%date%year = mtime_current%date%year
  
    IF (mtime_current >= datetime_after) THEN
      datetime_after%date%year = datetime_after%date%year + 1
    END IF     
  END IF ! datetime > datetime_after

  ! calculation of time fraction for linear interpolation
  CALL divideTwoDatetimeDiffsInSeconds(mtime_current,datetime_before,datetime_after,  &
                             &         datetime_before,denominator,quotient)

  IF (quotient%remainder_in_ms == 0) THEN
    x = 1._wp
  ELSE
    x = REAL(quotient%remainder_in_ms,wp) / 1000._wp / REAL(denominator,wp)
  END IF
    

  SELECT TYPE (presc => emiss_struct)
    TYPE IS (t_art_emiss_prescribed)
      !2D data
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, istart, iend)

        DO i = istart,iend      
          this_emiss_2d(i,jb) = emiss_2d_before(i,jb) & 
           &   + x * (emiss_2d_after(i,jb) - emiss_2d_before(i,jb))
        END DO
      END DO
    TYPE IS (t_art_chem_prescribed)
      !3D data
      DO jb = art_atmo%i_startblk, art_atmo%i_endblk
        CALL art_get_indices_c(jg, jb, istart, iend)

        DO jk = 1,art_atmo%nlev
          DO i = istart, iend
            this_emiss_3d(i,jk,jb,:) = emiss_3d_before(i,jk,jb,:)  &
             &  + x * (emiss_3d_after(i,jk,jb,:) - emiss_3d_before(i,jk,jb,:))
          END DO
        END DO
      END DO

  END SELECT

END SUBROUTINE art_read_emissions
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_determine_datetimes_before_or_after(p_patch,emiss_struct,index_emiss,diff_orig_year)
!<
! SUBROUTINE art_determine_datetimes_before_or_after
! This subtroutine calculates datetime_before and datetime_after
! depending on whether datetime_before is before first_date or datetime_after is after last_date.
! Part of Module: mo_art_read_emission
! Author: Michael Weimer, KIT
! Initial Release: about 2015-04-01
! Modifications:
! 2017-06-16: Michael Weimer, KIT
! - completely renewed this subroutine so that it also works with data sets
!   with length less than a year
!>
  IMPLICIT NONE
  TYPE(t_patch), INTENT(IN) :: p_patch      !< patch on which calculation is performed
  TYPE(datetime), POINTER :: date         !< datetime_before or datetime_after
  TYPE(datetime), POINTER :: ref1, ref2   !< reference datetimes (either first_date or last_date)
  CLASS(t_art_prescribed), INTENT(INOUT), TARGET ::  &
      &  emiss_struct                   !< same as in art_read_emissions
  INTEGER, INTENT(IN) :: index_emiss    !< necessary parameter for 
                                        !  art_search_read_file_around_datetime: index of datetime
                                        !  and emission in emiss_struct
  INTEGER(i8), INTENT(IN) ::  &
    & diff_orig_year

  INTEGER :: &
    &  year_hour_inc            !< either -1 or 1, depends on direction in which searching has to 
                                !  be done
  TYPE(datetime), POINTER :: &
    &  datetime_tmp => NULL()   !< temporarily saved calday + caltime of date1

  TYPE(timedelta), POINTER  ::&
    &  difference1 => NULL(), &
    &  difference2 => NULL() !< temporarily saved time difference between date next time with
                             !  emission file

  LOGICAL ::          &
    &  condition,     & !< conditions comparing date with reference dates
    &  condition2,    &
    &  not_condition2, &
    &  condition3

  ! Here, we give an example for index_emiss == 1, so selecting the correct file
  ! before simulation time if the simulation time exceeds last date and is
  ! beyond the data set

  ! date is the datetime_before
  date => emiss_struct%datetime(index_emiss)%ptr

  SELECT CASE (index_emiss)
    CASE(1)
      ! we want to determine the file at the end of the dataset which is closest
      ! to the simulation time: Is it last_date or is it another file within the
      ! dataset?

      ! Since we want to search for datetime_before the searching direction is
      ! negative in time (looking for the file before date)
      year_hour_inc = -1
      ref1  => emiss_struct%first_date
      ref2  => emiss_struct%last_date
      condition = (ref1 < date)
      condition2 = (date < ref2)
      not_condition2 = (date > ref2)
    CASE(3)
      year_hour_inc = +1
      ref1 => emiss_struct%last_date
      ref2 => emiss_struct%first_date
      condition = (date < ref1)
      not_condition2 = (ref2 > date)
      condition2 = (ref2 < date)
    CASE DEFAULT
      CALL finish('mo_art_read_emissions:art_determine_datetimes_before_or_after', &
              & 'parameter index_emiss has to be 1 or 3.')
  END SELECT

  ! if first date is earlier than datetime_before everything is okay. If not,
  ! just select last date as datetime_before ("else" part)
  IF (condition) THEN
    difference1 => newTimeDelta('PT0H')
    difference2 => newTimeDelta('PT0H')

    IF (condition2) THEN
      ! if datetime_before is earlier than last_date it is within the boundaries
      ! of the dataset, so temporarily save dateime_before, search for the file before 
      ! datetime_before and calculate the difference between them

      ! Additionally, calculate the difference of datetime_before with year
      ! increased by 1 to last_date
      datetime_tmp => newDatetime(date)

      CALL art_search_read_file_around_datetime(p_patch,emiss_struct,year_hour_inc,index_emiss, &
        &                                       diff_orig_year,.FALSE.)
      difference2 = datetime_tmp - date

      datetime_tmp%date%year = datetime_tmp%date%year - year_hour_inc

      difference1 = datetime_tmp - ref2

      CALL deallocateDatetime(datetime_tmp)
    ELSE IF (not_condition2) THEN
      ! Otherwise, do it the other way around: first calculate the difference to
      ! last date then calculate the difference to the closest file within the
      ! data set
      difference1 = date - ref2
      
      date%date%year = date%date%year + year_hour_inc

      SELECT CASE (index_emiss)
        CASE(1)
          condition3 = (date > ref1)
        CASE(3)
          condition3 = (ref1 > date)
      END SELECT

      IF (condition3) THEN
        datetime_tmp => newDatetime(date)

        CALL art_search_read_file_around_datetime(p_patch,emiss_struct,year_hour_inc,index_emiss, &
          &                                       diff_orig_year + year_hour_inc,.FALSE.)
        difference2 = datetime_tmp - date

        CALL deallocateDatetime(datetime_tmp)
      ELSE
        ! if datetime_before is beyond first date, data set is shorter than a
        ! year, so just select last date as datetime_before
        difference2 = difference1
      END IF
    END IF

    ! create absolute values of the differences
    difference1%sign = '+'
    difference2%sign = '+'

    ! if difference to last date is equal or lower than the other difference
    ! choose last date as datetime_before
    IF (difference2 >= difference1) THEN
      date = ref2
    END IF

    CALL deallocateTimeDelta(difference1)
    CALL deallocateTimeDelta(difference2)
  ELSE
    date = ref2
  END IF
    
END SUBROUTINE art_determine_datetimes_before_or_after
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_search_read_file_around_datetime(p_patch,emiss_struct,hour_inc, index_emiss, &
  &                                             diff_orig_year, with_read)
!<
! SUBROUTINE art_search_read_file_around_datetime
! This subroutine searches the next file before or after (depending on hour_inc) 
! emiss_struct%datetime(index_emiss) and reads the file if that is wanted
! Part of Module: mo_art_read_emission
! Author: Michael Weimer, KIT
! Initial Release: about 2015-04-01
! Modifications:
! 2016-06-19: Michael Weimer, KIT
! - Included correct treatment of leap years
!>
  IMPLICIT NONE
  TYPE(t_patch), INTENT(IN) :: &
    &  p_patch          !< patch on which calculation is performed
  CLASS(t_art_prescribed), INTENT(INOUT) ::   &
    &  emiss_struct     !< same as in art_read_emissions
  INTEGER, INTENT(IN) :: &
     &  hour_inc,    &  !< direction of search by choosing the increment either -1, 0 or +1 for
                        !  searching before, here or after datetime, respectively
     &  index_emiss     !< index of emission for which datetime has to be calculated (either 1 or 3
                        !  for emission_before or emission_after)
  TYPE(timedelta), POINTER :: &
    &  td_one_hour => NULL()  !< timedelta of one hour
  INTEGER(i8), INTENT(IN) ::  &
    &  diff_orig_year   !< difference of the years between current simulation time 
                        !  and the actual searched file
  LOGICAL, OPTIONAL, INTENT(IN) ::  &
    &  with_read        !< should the emission be read actually?

! local variables
  INTEGER :: n          !< loop index
  LOGICAL :: l_exist,  l_with_read !< logicals if file exists or emission should be read
  CHARACTER(LEN = max_datetime_str_len)  ::  &
    &  datetime_str     !< string of a datetime
  CHARACTER(LEN = IART_PATH_LEN) ::  &
    &  emissiondataset  !< path and name of the searched file
  INTEGER(i8) ::  &
    & year              !< year of the date of the searched emission file
  INTEGER :: &
    & jg,    &
    & day, month,  &      !< day and month of the searched emission file
    & idx_date, len_date  !< index of the date and length of the date in the
                          !  ICON-ART name convention
  TYPE(timedelta), POINTER :: &
    &  td_diff_years => NULL()  !< time delta with diff_orig_year as years
  CHARACTER (LEN = 7) ::  &
    &  td_str

  ! --------------------------------------
  ! 1. Initializations
  ! --------------------------------------

  jg = p_patch%id

  IF (diff_orig_year >= 0) THEN
    WRITE(td_str,'(A2,I4.4,A1)') '+P', ABS(diff_orig_year), 'Y'
  ELSE
    WRITE(td_str,'(A2,I4.4,A1)') '-P', ABS(diff_orig_year), 'Y'
  END IF

  td_diff_years => newTimeDelta(TRIM(td_str))

  IF (PRESENT(with_read)) THEN
    l_with_read = with_read
  ELSE
    l_with_read = .TRUE.
  END IF

  n       = 0
  l_exist = .FALSE.

  SELECT CASE (hour_inc)
    CASE(-1)
      td_one_hour => newTimeDelta('-PT1H')
    CASE(0)
      td_one_hour => newTimeDelta('PT0H')
    CASE(1)
      td_one_hour => newTimeDelta('PT1H')
    CASE DEFAULT
      CALL finish('mo_art_read_emissions: art_search_read_file_around_datetime', &
           & 'hour_inc neither -1, 0 nor 1.')
  END SELECT

  CALL datetimeToPosixString(emiss_struct%datetime(index_emiss)%ptr,datetime_str,"%Y-%m-%d-%H")


  
  ! in art_create_filenames, many things are checked which is inefficient if it
  ! is done within the while loop below. Therefore, it is done once and the
  ! position and length of the date within the file name is determined
  CALL art_create_filenames(jg,TRIM(emiss_struct%path),emiss_struct%iType_data, &
    &                       emissiondataset,.TRUE.,TRIM(datetime_str)  )

  len_date = LEN_TRIM(datetime_str)
  idx_date = INDEX(emissiondataset,TRIM(datetime_str))

  ! --------------------------------------
  ! 2. Execute searching and reading:
  !    - while file cannot be found add or subtract one hour of datetime_after or datetime_before,
  !      respectively
  !    - if file is found read it (if exception with leap years is okay)
  ! --------------------------------------
  
  ! Stop if a file cannot be found about 11 years before or after datetime
  DO WHILE ((n < 100000) .AND. (.NOT. l_exist))
    n = n+1

    emiss_struct%datetime(index_emiss)%ptr =  &
      &  emiss_struct%datetime(index_emiss)%ptr + td_one_hour


    CALL datetimeToPosixString(emiss_struct%datetime(index_emiss)%ptr,datetime_str,"%Y-%m-%d-%H")
    
    emissiondataset(idx_date:idx_date + len_date-1) = TRIM(datetime_str)
    
    INQUIRE(FILE = TRIM(emissiondataset), EXIST=l_exist)
    

    ! There is an exception with resepct to leap years:
    ! if e.g. last year is a leap year and the simulation time is in no leap
    ! year after last date it can happen that here a data set of 29 February is
    ! read which is however a non-existing date in the simulation time's year.
    ! So, data sets of 29 February are excluded in this case

    ! The following if statements are ordered according to their probable
    ! occurrence in the model (beginning with the most unprobable one).
    day = emiss_struct%datetime(index_emiss)%ptr%date%day
    month = emiss_struct%datetime(index_emiss)%ptr%date%month
    year = emiss_struct%datetime(index_emiss)%ptr%date%year

    ! if date of found emission file is 29 February
    IF ((month == 2) .AND. (day == 29) .AND. (l_exist)) THEN

      ! if one of the boundary dates is a leap year
      IF ((emiss_struct%first_date_is_leap_year)  &
           .OR. (emiss_struct%last_date_is_leap_year)) THEN
       
        ! if year of found file is beyond boundary date (so if year has to be
        ! changed after the reading)
        IF ((year + diff_orig_year > emiss_struct%last_date%date%year)  &
           .OR. (year + diff_orig_year < emiss_struct%first_date%date%year)) THEN

          ! if the actual date of the emission (before or after simulation time)
          ! is no leap year
          IF (getNoOfDaysInYearDateTime(  &
             &  emiss_struct%datetime(index_emiss)%ptr + td_diff_years) /= 366)  &
             &  THEN
            l_exist = .FALSE.
          END IF
        END IF
      END IF
    END IF

  END DO


  IF (n >= 100000) THEN
     CALL finish('mo_art_read_emissions: art_search_read_file_around_datetime', &
           & 'no file found around datetime. Last searched file:'//TRIM(emissiondataset))
  ELSE
    IF (l_with_read) THEN
      CALL art_open_read_close(p_patch,emissiondataset,emiss_struct,index_emiss)
    END IF
  END IF

  CALL deallocateTimeDelta(td_one_hour)
  CALL deallocateTimeDelta(td_diff_years)

END SUBROUTINE art_search_read_file_around_datetime
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_open_read_close(p_patch,file_namepath,emiss_struct,index_emiss)
!<
! SUBROUTINE art_open_read_close
! This subroutine opens, reads and closes a committed dataset and writes it into 
! emiss_struct%emission(index_emiss).
! In addition, it copies the emission into user given lowest levels in case of
! two-dimensional emission data.
! Part of Module: mo_art_read_emission
! Author: Michael Weimer, KIT
! Initial Release: about 2015-04-01
! Modifications:
! 2017-08-11: Michael Weimer, KIT
! - included vertical interpolation of 3D datasets (also for pressure levels)
!>
  IMPLICIT NONE
  TYPE(t_patch), INTENT(IN)      :: &
    &  p_patch              !< patch on which computation is performed
  CHARACTER(LEN = *), INTENT(IN) ::  &
    &  file_namepath        !< path and file name to be read
  CLASS(t_art_prescribed), INTENT(INOUT) :: &
    &  emiss_struct         !< same as in art_read_emissions
  INTEGER, OPTIONAL, INTENT(IN)  ::  &
    &  index_emiss          !< index of the emissions

  INTEGER ::                 &
    &  jg,                   &
    &  varidx,               & !< loop indices
    &  jc, jk, jb,           &
    &  i_startidx,           &
    &  i_endidx
  TYPE(t_stream_id) :: stream_id

  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  jg       =  p_patch%id
  art_atmo => p_art_data(jg)%atmo

  CALL openInputFile( stream_id, TRIM(file_namepath), p_patch, &
    &                 default_read_method)

  SELECT TYPE(presc => emiss_struct)
    !-------
    TYPE IS (t_art_emiss_prescribed)
    !-------

      presc%emiss_2d(:,:,index_emiss) = 0.0_wp

      DO varidx = 1,presc%num_vars
        presc%emiss_2d_read(:,:) = 0.0_wp
        CALL message('mo_art_read_emissions:art_open_read_close',  &
               &     'Reading '//TRIM(presc%vname(varidx))//' from '//TRIM(file_namepath)//'.')
        CALL read_2d_1time(stream_id, on_cells,TRIM(presc%vname(varidx)),   &
                  &        fill_array= presc%emiss_2d_read(:,:))
        presc%emiss_2d(:,:,index_emiss) =  &
          &   presc%emiss_2d(:,:,index_emiss) + presc%emiss_2d_read(:,:)
  
      END DO

    !-------
    TYPE IS (t_art_chem_prescribed)
    !-------
      ! This should be SI units but for testing use this one
      IF (TRIM(presc%unit_vertical) == 'hPa') THEN
        CALL art_prepare_vinterp_pres(jg, presc%prescribed_state,         &
                 &                    TRIM(file_namepath))

      ELSE
        IF (.NOT. presc%prescribed_state%chem_init_in%linitialized) THEN
          CALL art_prepare_vinterp(presc%prescribed_state,p_patch,                         &
                  &                TRIM(presc%path)//'/'//TRIM(presc%vert_coord_filename), &
                  &                TRIM(file_namepath))
        END IF
      END IF

      IF (presc%prescribed_state%chem_init_in%nlev_chem_init == 1) THEN

        presc%vinterp_3d(:,:,:,:,index_emiss) = 0._wp

        DO varidx = 1,presc%num_vars
          CALL message('mo_art_read_emissions:art_open_read_close',  &
                 &     'Reading '//TRIM(presc%vname(varidx))//' from '//TRIM(file_namepath)//'.')

           CALL read_2d_1lev_1time(stream_id, on_cells,TRIM(presc%vname(varidx)),           &
                      &     fill_array=presc%prescribed_state%chem_init_chem%spec(:,1,:,varidx))

          DO jk = presc%var_dep(varidx)%upp_idx,presc%var_dep(varidx)%bot_idx
            presc%vinterp_3d(:,jk,:,varidx,index_emiss) =                       &
                 &   presc%prescribed_state%chem_init_chem%spec(:,1,:,varidx)   &
                 &   * presc%var_dep(varidx)%mol_weight / amd * 1000._wp
          END DO

        END DO

      ELSE
        DO varidx = 1,presc%num_vars
          CALL message('mo_art_read_emissions:art_open_read_close',  &
                 &     'Reading '//TRIM(presc%vname(varidx))//' from '//TRIM(file_namepath)//'.')
          CALL read_3d_1time(stream_id, on_cells,TRIM(presc%vname(varidx)), &
                 &           fill_array=presc%prescribed_state%chem_init_chem%spec(:,:,:,varidx))

          DO jb = art_atmo%i_startblk, art_atmo%i_endblk
            CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
 
            DO jk = 1, presc%prescribed_state%chem_init_in%nlev_chem_init
              DO jc = i_startidx, i_endidx
                presc%prescribed_state%chem_init_chem%spec(jc,jk,jb,varidx)          &
                  &   = presc%prescribed_state%chem_init_chem%spec(jc,jk,jb,varidx)  &
                  &     * presc%var_dep(varidx)%mol_weight / amd * 1000._wp
              END DO
            END DO
          END DO

        END DO

        CALL art_vinterp_chem_init(presc%prescribed_state,presc%num_vars,jg)

        presc%vinterp_3d(:,:,:,:,index_emiss)                                      &
               &    = presc%prescribed_state%chem_init_chem%spec_interp(:,:,:,:)

      END IF

  END SELECT

  CALL closeFile(stream_id)

END SUBROUTINE art_open_read_close
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_add_emission_to_tracers(p_tracer_now,emiss,p_patch,dtime,mtime_current)
!<
! SUBROUTINE art_add_emission_to_tracer
! This subroutine calculates the current emission or computes online emissions
! and adds them to the tracer mixing ratio
! Part of Module: mo_art_read_emissions
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-09
! Modifications:
!>
  IMPLICIT NONE

  TYPE(t_patch), INTENT(IN) :: &
   &  p_patch                  !< patch on which computation is performed

  REAL(wp), INTENT(INOUT) :: &
   &  p_tracer_now(:,:,:,:)     !< tracer mixing ratio

  REAL(wp), INTENT(IN) :: &
   &  dtime                     !< time step    

  TYPE(datetime), POINTER ::  &
   &  mtime_current             !< current datetime

  TYPE(t_art_emiss_storage), INTENT(INOUT) :: &
    &  emiss                    !< emission storage

  !local variables
  REAL(wp), ALLOCATABLE ::  &
   &  cbio(:,:,:),          &  !< biogenic emission of MEGAN species
   &  emiss_mmr(:)             !< emission value converted to mass mixing ratio

  INTEGER :: i,j,itr,        & !< loop index
    &    jk, idx_cbio          !< loop index and index of biogenic emission

  REAL(wp), ALLOCATABLE :: dummy_wp(:)   !< for online biogenic emissions some of the 
                                         !    parameters are not needed
  REAL(wp), POINTER :: &
    &  pft(:,:,:)              !< plant functional types (read from external file)

  TYPE(t_art_emiss_type_container)  :: &
    &  tracer_emission         !< temporal save of the actual emission container

  INTEGER    :: iget_error  !< error handler for getting an emission element

  INTEGER :: jb, istart, iend , jg                 !< loop indices

  TYPE(t_art_atmo),POINTER    :: &
      &  art_atmo                     !< Pointer to ART atmo fields

  !!------------------------------------------------------------------------------
  !! --- Time dependent emissions ---
  !!------------------------------------------------------------------------------

  jg   = p_patch%id

  art_atmo => p_art_data(jg)%atmo

  ALLOCATE(dummy_wp(art_atmo%nproma))
  
  !-------------------------------------------------------------------------------
  ! --- Initialisation
  !-------------------------------------------------------------------------------
  
  ALLOCATE(cbio(art_atmo%nproma,art_atmo%nblks,sel_bccnum))
  ALLOCATE(emiss_mmr(art_atmo%nlev))
  

  cbio(:,:,:) = 0._wp
  dummy_wp(:) = 1._wp
  
  !-------------------------------------------------------------------------------
  ! --- read emission data when datetime has changed and convert the values to vmr
  !-------------------------------------------------------------------------------
  
  DO itr = 1,ntracer
  
    CALL emiss%get(itr,tracer_emission,iget_error)
  
    IF (iget_error == 0) THEN
      IF (tracer_emission%bioonl%idx > 0) THEN
        IF (MAXVAL(cbio) <= 0.0_wp) THEN
          pft => tracer_emission%bioonl%pft

          DO jb = art_atmo%i_startblk, art_atmo%i_endblk
            CALL art_get_indices_c(jg, jb, istart, iend)

            CALL bvoc_guenther2012(art_atmo%temp(:,art_atmo%nlev,jb),     &
                               &   art_atmo%swflx_par_sfc(:,jb),          &
!                               &   dummy_wp,                              &
!                               &   dummy_wp,                              &
!                               &   dummy_wp,                              &
                               &   pft(:,jb,:),                           &
                               &   art_atmo%lai(:,jb),                    &
!                               &   art_atmo%sza_deg(:,jb),                &
                               &   istart,                                &
                               &   iend,                                  &
                               &   cbio(:,jb,:),                          &
                               &   dummy_wp)

            DO i = istart,iend
              DO idx_cbio = 1,sel_bccnum
                ! cbio is given in units of kg m-2 h-1 and has to be divided by 3600 
                ! to convert it to SI units
                cbio(i,jb,idx_cbio) = cbio(i,jb,idx_cbio) / 3600._wp
              END DO

            END DO
          END DO

          NULLIFY(pft)
        END IF ! cbio <= 0

        DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          CALL art_get_indices_c(jg, jb, istart, iend)
          DO i = istart,iend
            tracer_emission%bioonl%cbio(i,jb) = cbio(i,jb,tracer_emission%bioonl%idx)
          END DO
        END DO

        IF (.NOT. art_config(jg)%lart_emiss_turbdiff) THEN

          DO jb = art_atmo%i_startblk, art_atmo%i_endblk
            CALL art_get_indices_c(jg, jb, istart, iend)
            DO i = istart,iend
              CALL art_convert_emission_to_mmr(emiss_mmr(:),                        &
                      &                        tracer_emission%bioonl%cbio(i,jb),   &
                      &                        art_atmo%rho(i,:,jb),                &
                      &                        art_atmo%dz(i,:,jb),                 &
                      &                        dtime,                               &
                      &                        tracer_emission%bioonl%num_emiss_lev,&
                      &                        art_atmo%nlev)
    
              DO jk = art_atmo%nlev,art_atmo%nlev-tracer_emission%bioonl%num_emiss_lev+1,-1
                p_tracer_now(i,jk,jb,itr) = p_tracer_now(i,jk,jb,itr) &
                      &                     + emiss_mmr(jk)*tracer_emission%bioonl%scaling_factor
              END DO
            END DO
          END DO
        END IF ! lart_emiss_turbdiff
  
      END IF ! bioonl > 0
  
      IF (tracer_emission%num_types_prescribed == 0) THEN
        ! if standard_value of emission is zero, nothing has to be done
        IF (tracer_emission%std%val > 0.0_wp) THEN
          SELECT CASE (tracer_emission%std%mode)
          CASE(1)
            IF (.NOT. art_config(jg)%lart_emiss_turbdiff) THEN
              DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                CALL art_get_indices_c(jg, jb, istart, iend)
                DO i = istart, iend
                  CALL art_convert_emission_to_mmr(emiss_mmr(:),                     &
                                 &                 tracer_emission%std%val,          &
                                 &                 art_atmo%rho(i,:,jb),             &
                                 &                 art_atmo%dz(i,:,jb),              &
                                 &                 dtime,                            &
                                 &                 tracer_emission%std%num_emiss_lev,&
                                 &                 art_atmo%nlev)
                  p_tracer_now(i,:,jb,itr) = p_tracer_now(i,:,jb,itr) + emiss_mmr(:)
                END DO
              END DO
            END IF
          CASE(2)
            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
              CALL art_get_indices_c(jg, jb, istart, iend)
              DO i = istart, iend
                p_tracer_now(i,art_atmo%nlev,jb,itr) = p_tracer_now(i,art_atmo%nlev,jb,itr)  &
                    &     + tracer_emission%std%val * dtime
              END DO
            END DO
          END SELECT
           
        END IF ! standard value > 0
      
      ELSE ! num_types_prescribed == 0
  
        DO j = 1,tracer_emission%num_types_prescribed
          IF (.NOT. tracer_emission%types(j)%type_is_init) THEN
            CALL art_init_emission_struct(p_patch,mtime_current,tracer_emission%types(j))
          ELSE
            CALL art_read_emissions(p_patch,mtime_current,tracer_emission%types(j))
          END IF

          ! 2D emissions
          IF (.NOT. art_config(jg)%lart_emiss_turbdiff) THEN
            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
              CALL art_get_indices_c(jg, jb, istart, iend)

              DO i = istart, iend
                CALL art_convert_emission_to_mmr(emiss_mmr(:),                               &
                                &                tracer_emission%types(j)%emiss_2d(i,jb,2),  &
                                &                art_atmo%rho(i,:,jb),                       &
                                &                art_atmo%dz(i,:,jb),                        &
                                &                dtime,                                      &
                                &                tracer_emission%types(j)%num_emiss_lev,     &
                                &                art_atmo%nlev)

                DO jk = art_atmo%nlev,art_atmo%nlev-tracer_emission%types(j)%num_emiss_lev+1,-1
                  p_tracer_now(i,jk,jb,itr) = p_tracer_now(i,jk,jb,itr) + emiss_mmr(jk)  &
                        &                     * tracer_emission%types(j)%scaling_factor
                END DO
              END DO 
            END DO
          END IF ! lart_emiss_turbfdiff
        END DO ! j

      END IF  ! num_types_prescribed == 0

      CALL emiss%update(itr,tracer_emission,iget_error)
    END IF ! iget_error == 0
  
  END DO  ! iTR
  !!------------------------------------------------------------------------------
  !! --- Deallocation
  !!------------------------------------------------------------------------------
    

  DEALLOCATE(cbio)
  DEALLOCATE(emiss_mmr)
  
  NULLIFY(art_atmo)
END SUBROUTINE art_add_emission_to_tracers
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_convert_emission_to_mmr(emiss_mmr,emiss_kgm2s1,rho,height_layers, &
  &                                    dtime,num_emiss_lev, nlev,emiss_3d)
!<
! SUBROUTINE art_convert emission to mmr
! This subroutine converts the emission as mass flux density in kg m-2 s-1 to
! mass mixing ratio
! Part of Module: mo_art_read_emissions
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-09
! Modifications:
!>
  IMPLICIT NONE
  REAL(wp), INTENT(INOUT) :: emiss_mmr(:)   !< emission in units of kg kg-1
  REAL(wp), INTENT(IN)  ::  &
    & emiss_kgm2s1,         &             !< emission in units of kg m-2 s-1
    & rho(:),               &             !< density
    & height_layers(:),     &             !< height of the model layers
    & dtime                               !< time step
  REAL(wp), INTENT(IN), OPTIONAL :: &
    &  emiss_3d(:)                        !< 3D emission

  INTEGER, INTENT(IN) ::    &
    & num_emiss_lev, nlev                 !< number of model layers in which emission should be
                                          !  included  / total number of model levels
  ! local variables
  INTEGER  :: jk                          !< vertical loop index
  REAL(wp) :: kg_air(nlev)                !< mass of the air relative to grid box base area

  emiss_mmr(:) = 0.0_wp
  kg_air    = 0.0_wp

  DO jk = 1,nlev
    kg_air(jk) = rho(jk) * height_layers(jk)

    IF (PRESENT(emiss_3d)) THEN
      emiss_mmr(jk) = emiss_3d(jk) * dtime  / kg_air(jk)
    END IF
  END DO

  DO jk = nlev,nlev-num_emiss_lev+1,-1
    emiss_mmr(jk) = emiss_kgm2s1 * dtime  &
                 / SUM(kg_air(nlev-num_emiss_lev+1:nlev))
  END DO

END SUBROUTINE art_convert_emission_to_mmr
!
!--------------------------------------------------------------
!

END MODULE mo_art_read_emissions



