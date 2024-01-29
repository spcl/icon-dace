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

MODULE mo_art_oem_emission

!------------------------------------------------------------------------------
!
! Description:
!   This module contains subroutines for computatoin of the emissions by 
!   scaling the gridded emissions with temporal and vertical profiles and
!   adds them to the tracers.
!==============================================================================
  ! ICON
  USE mo_kind,                   ONLY: wp, i8
  USE mo_exception,              ONLY: message, message_text
  USE mo_model_domain,           ONLY: t_patch, p_patch
  USE mo_nonhydro_state,         ONLY: p_nh_state
  USE mo_parallel_config,        ONLY: nproma, idx_1d, blk_no,     &
                                   &   idx_no
  USE mo_grid_config,            ONLY: n_dom
  USE mo_dynamics_config,        ONLY: nnow_rcf
  USE mo_loopindices,            ONLY: get_indices_c
  USE mo_impl_constants,         ONLY: min_rlcell
  USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c

  USE mo_time_config,            ONLY: time_config !configure_time
  USE mtime,                     ONLY: MAX_DATETIME_STR_LEN, &
                                   &   datetimeToString,     &
                                   &   julianday,            &
                                   &   newJulianday,         &
                                   &   getJulianDayFromDatetime,  &
                                   &   getDatetimeFromJulianDay,  &
                                   &   getnoofdaysinyeardatetime, &
                                   &   datetime,             &
                                   &   no_of_ms_in_a_day,    &
                                   &   timedelta,            &
                                   &   newTimedelta,         &
                                   &   newDatetime,          &
                                   &   OPERATOR(+)

  ! new
  USE mo_async_latbc_types,      ONLY: t_latbc_data
  USE mo_initicon_types,         ONLY: t_pi_atm
  USE mo_util_mtime,             ONLY: getElapsedSimTimeInSeconds

  ! ART
  USE mo_art_atmo_data,          ONLY: t_art_atmo
  USE mo_art_data,               ONLY: p_art_data
  USE mo_art_wrapper_routines,   ONLY: art_get_indices_c

  ! OEM
  USE mo_art_oem_types,          ONLY: p_art_oem_data,          &
                                   &   t_art_oem_data,          &
                                   &   t_art_oem_config,        &
                                   &   t_art_oem_ensemble


!---------------------------------------------------------------------------

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: art_oem_compute_emissions,       &
    &       art_oem_extract_time_information

  ! Constant variable
  INTEGER,  PARAMETER :: ncat_max = 200
  INTEGER,  PARAMETER :: tp_param_hourofday = 24
  INTEGER,  PARAMETER :: tp_param_dayofweek = 7
  INTEGER,  PARAMETER :: tp_param_monthofyear = 12
  INTEGER(KIND=2), PARAMETER :: tp_param_hour = 8784

  TYPE(julianday), POINTER :: jd, jdref

  CHARACTER(LEN=3), DIMENSION(7) :: day_of_week = &
     & (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)

  TYPE(t_art_oem_data),     POINTER :: oem_data         !< OEM data structure -> data
  TYPE(t_art_oem_config),   POINTER :: oem_config
  TYPE(t_art_oem_ensemble), POINTER :: oem_ensemble


!==============================================================================
! Module procedures
!==============================================================================


CONTAINS


!==============================================================================
!==============================================================================
!+ Initialize and read gridded emissions from NetCDF file
!------------------------------------------------------------------------------

  SUBROUTINE art_oem_compute_emissions(p_tracer_now,p_patch,dtime,mtime_current,ierror,yerrmsg)

!-----------------------------------------------------------------------------
! Description: This subroutine uses the temporal and vertical profiles to 
! compute the gridded emissions and adds them to the OEM-tracer
!-----------------------------------------------------------------------------

  IMPLICIT NONE

  REAL(wp), INTENT(inout) :: &
   &  p_tracer_now(:,:,:,:)     !< tracer mixing ratio

  TYPE(t_patch), INTENT(in) :: &
   &  p_patch                  !< patch on which computation is performed

  REAL(wp), INTENT(in) :: &
   &  dtime                     !< time step    

  TYPE(datetime), POINTER ::  &
   &  mtime_current             !< current datetime

  INTEGER, INTENT(OUT)             :: ierror
  CHARACTER(LEN= *), INTENT(OUT)   :: yerrmsg


  !---------------------------------------------------------------------------
  ! Local variables
  INTEGER :: grd_index, hod, dow, moy, hoy, cat, &
    &        hod_n, dow_n, moy_n, hoy_n, minutes, &
    &        tp_cat_idx, vp_cat_idx, tp_country_idx, &
    &        gridded_idx, nc, nt, is, ie, jg, &
    &        nlev, k, nr, nens, ens_int_idx, trcr_idx, &
    &        jb, jc, nblks_c, i_startblk, i_endblk, &
    &        table_nr, ens_count, nt_old

  REAL(KIND=wp) :: newemis, temp_scaling_fact, lambda, z_mass, sim_time, &
    &              n1, n2, temp_scaling_fact_now, temp_scaling_fact_next

  CHARACTER(LEN=2) :: numstring

  TYPE(datetime) :: datetime_next

  TYPE(timedelta), POINTER :: mtime_td

  TYPE(t_art_atmo), POINTER :: art_atmo

  CHARACTER(*), PARAMETER :: routine = "art_oem_compute_emissions"

  ierror = 0
  yerrmsg = ''

!- End of header
!==============================================================================

    jg   = p_patch%id
    art_atmo => p_art_data(jg)%atmo

    oem_data => p_art_oem_data%data_fields
    oem_config => p_art_oem_data%configure
    oem_ensemble => p_art_oem_data%ensemble

    nblks_c = art_atmo%nblks
    nlev = art_atmo%nlev

!------------------------------------------------------------------------------
! Section 1: Calculate emission field
!------------------------------------------------------------------------------

    IF (oem_config%emis_tracer>0) THEN

      ! read start- and end-block for this PE:
      i_startblk = art_atmo%i_startblk
      i_endblk   = art_atmo%i_endblk

      ! Extract the different time information for this timestep
      CALL art_oem_extract_time_information(time_config%tc_current_date,hod,dow,moy,hoy,minutes)
      ! Extract the time information for one hour later
      mtime_td => newTimedelta("PT01H")
      datetime_next = time_config%tc_current_date + mtime_td
      CALL art_oem_extract_time_information(datetime_next,hod_n,dow_n,moy_n,hoy_n,minutes)
      n2 = minutes/60._wp
      n1 = 1._wp - n2

      nt_old = 0
      ens_count = 0
      DO nt = 1,oem_config%emis_tracer
        trcr_idx = oem_config%emis_idx(nt)
        DO nc = 1, ncat_max
          IF(oem_config%ycatl_l(nt,nc) /= "") THEN
            ! get index of this category
            grd_index = get_gridded_emissions_idx(oem_config%ycatl_l(nt,nc))
            vp_cat_idx = get_category_idx(oem_config%yvpl_l(nt,nc), 'vp')
            tp_cat_idx = get_category_idx(oem_config%ytpl_l(nt,nc), 'tp')

            DO jb = i_startblk, i_endblk
              ! read indices within this block:
              CALL art_get_indices_c(jg, jb, is, ie)
              DO jc = is, ie
                newemis = oem_data%gridded_emissions(jc,jb,grd_index)
                tp_country_idx = get_country_idx(oem_data%country_ids(jc,jb))

                ! get temporal scaling factor for this and the next hour and interpolate linearly
                temp_scaling_fact = 0._wp

                SELECT CASE (oem_config%itype_tscale_l(nt))
                  CASE(0)
                     temp_scaling_fact = 1._wp
                  CASE(1)
                    temp_scaling_fact_now = temporal_scaling_factor(tp_cat_idx,  &
                       & hod, dow, moy, tp_country_idx)
                    temp_scaling_fact_next = temporal_scaling_factor(tp_cat_idx,  &
                       & hod_n, dow_n, moy_n, tp_country_idx)
                    temp_scaling_fact = n1*temp_scaling_fact_now + n2*temp_scaling_fact_next
                  CASE(2)
                    temp_scaling_fact_now = oem_data%tp_hourofyear(hoy, tp_cat_idx, tp_country_idx)
                    temp_scaling_fact_next = oem_data%tp_hourofyear(hoy_n, tp_cat_idx, tp_country_idx)
                    temp_scaling_fact = n1*temp_scaling_fact_now + n2*temp_scaling_fact_next
                END SELECT

                lambda = 1._wp
                !------------------------------------------------------------------------------
                ! Section 1.1: Emission fields for ensemble members
                !------------------------------------------------------------------------------
                IF ( ANY( oem_ensemble%ens_name==oem_config%emis_name(nt) ) ) THEN
                  IF (nt_old/=nt) THEN
                    ens_count = ens_count+1
                    nt_old = nt
                  ENDIF
                  nr = oem_data%reg_map(jc,jb)
                  DO nens = 1,SIZE(oem_data%lambda_mat, dim=4)
                    DO table_nr=1,200
                      IF (oem_ensemble%ens_table(1,table_nr)==nens .AND. oem_ensemble%ens_name(table_nr)==oem_config%emis_name(nt)) THEN
                        lambda = oem_data%lambda_mat(ens_count,nc,nr,nens)
                        !read index of ensemble tracer
                        ens_int_idx = oem_ensemble%ens_table(2,table_nr)
                        DO k = 1,nlev
                          ! calculate air mass
                          z_mass = p_nh_state(jg)%diag%airmass_now(jc,k,jb)
                          p_tracer_now(jc,k,jb,ens_int_idx) = p_tracer_now(jc,k,jb,ens_int_idx) + dtime*newemis &
                            & * oem_data%vert_scaling_fact(jc,k,jb,vp_cat_idx) * temp_scaling_fact * lambda &
                            & / z_mass
                        ENDDO !k=1,nlev
                      ENDIF ! ens_table(1,table_nr)==nens
                    ENDDO ! table_nr=1,200
                  ENDDO ! nens = 1,SIZE(lambda_mat, dim=4)
                ENDIF !lens(nt)==.TRUE.
                !------------------------------------------------------------------------------
                ! Section 1.2: Emission fields for non ensemble members
                !------------------------------------------------------------------------------
                DO k = 1,nlev
                  ! calculate air mass
                  z_mass = p_nh_state(jg)%diag%airmass_now(jc,k,jb)
                  p_tracer_now(jc,k,jb,trcr_idx) = p_tracer_now(jc,k,jb,trcr_idx) + dtime * newemis &
                    & * oem_data%vert_scaling_fact(jc,k,jb,vp_cat_idx) * temp_scaling_fact &
                    & / z_mass
                ENDDO
              ENDDO ! jc = is, ie
            ENDDO ! jb = i_startblk, i_endblk
          ENDIF ! (ycatl_trcr(nt,nc) /= "")
        ENDDO ! nc = 1,20
      ENDDO ! nt = 1,oem_config%emis_tracer
    ENDIF


!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

  END SUBROUTINE art_oem_compute_emissions



!==============================================================================
!==============================================================================
!+ Find the index of a gridded emission field based on its name
!------------------------------------------------------------------------------

  FUNCTION get_gridded_emissions_idx(name)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: get_gridded_emissions_idx
    INTEGER :: i

    get_gridded_emissions_idx = 0
    DO i = 1, SIZE(oem_config%gridded_emissions_idx)
      IF(name == oem_config%gridded_emissions_idx(i)) THEN
        get_gridded_emissions_idx = i
        EXIT
      END IF
    END DO

  END FUNCTION get_gridded_emissions_idx


!==============================================================================
!==============================================================================
!+ Extract time information for temporal scaling factor for the current UTC date
!------------------------------------------------------------------------------

  SUBROUTINE art_oem_extract_time_information(date, hour_of_day, day_of_week, month_of_year, &
    &                                         hour_of_year, minutes)

    IMPLICIT NONE

    ! Parameters
    TYPE(datetime), INTENT(IN) :: date
    INTEGER, INTENT(OUT) :: hour_of_day
    INTEGER, INTENT(OUT) :: day_of_week
    INTEGER, INTENT(OUT) :: month_of_year
    INTEGER, INTENT(OUT) :: hour_of_year
    INTEGER, INTENT(OUT) :: minutes

    ! Local variables
    CHARACTER(LEN=14) :: yactdate1 ! yyyymmddhhmmss
    CHARACTER(LEN=28) :: yactdate2 ! wd   dd.mm.yy  hh mm ss UTC
    INTEGER :: nactday ! day of the year
    REAL(KIND=wp) :: y1, y2, acthour ! actual hour of the day
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: time_string
    INTEGER :: errno, ndays
    TYPE(datetime) :: refdate


    CALL datetimeToString(time_config%tc_current_date, time_string)
    !format: 2017-01-12T05:57:00.000

    jd => newJulianday(0_i8, 0_i8)
    CALL getJulianDayFromDatetime(time_config%tc_current_date,jd,errno) !(dt, jd, errno)

    ndays = getNoOfDaysInYearDateTime(time_config%tc_current_date,errno) !(dt, errno)

    ! ref-date (01 Jan T00 of this year):
    refdate = time_config%tc_current_date
    refdate%date%month = 1
    refdate%date%day = 1
    refdate%time%hour = 0
    ! Julian day of refdate:
    jdref => newJulianday(0_i8, 0_i8)
    CALL getJulianDayFromDatetime(refdate,jdref,errno) !(dt, jd, errno)

    y1 = REAL(jd%day,wp) + REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)
    y2 = REAL(jdref%day,wp) + REAL(jdref%ms,wp)/REAL(no_of_ms_in_a_day,wp)

    hour_of_day = extract_hour(time_string)
    day_of_week = extract_day_of_week(jd)
    month_of_year = extract_month(time_string)
    hour_of_year = int(24._wp*(y1-y2)) !24 * (nactday-1) + hour_of_day ! 1 Jan, 00 UTC -> 1
    minutes = extract_minutes(time_string)

  END SUBROUTINE art_oem_extract_time_information


!==============================================================================
!==============================================================================
!+ Calculate the temporal scaling factor
!------------------------------------------------------------------------------

  FUNCTION temporal_scaling_factor(category_idx, hour_of_day, day_of_week, &
    &                              month_of_year, country_id)

    IMPLICIT NONE

    ! Result variable
    REAL(KIND=wp) :: temporal_scaling_factor

    ! Parameters
    INTEGER, INTENT(IN) :: category_idx
    INTEGER, INTENT(IN) :: hour_of_day     ! UTC
    INTEGER, INTENT(IN) :: day_of_week
    INTEGER, INTENT(IN) :: month_of_year
    INTEGER, INTENT(IN) :: country_id

    ! SC_TODO: hour_of_day should be local time and thus should be computed
    ! based on country_id and date (i.e. daylight saving time)

    temporal_scaling_factor = oem_data%tp_hourofday(hour_of_day, category_idx, country_id) * &
       &                      oem_data%tp_dayofweek(day_of_week, category_idx, country_id) * &
       &                      oem_data%tp_monthofyear(month_of_year, category_idx, country_id)

  END FUNCTION temporal_scaling_factor


!==============================================================================
!==============================================================================
!+ Get the category index for vertical or temporal profile
!------------------------------------------------------------------------------

  FUNCTION get_category_idx(category, list)
    IMPLICIT NONE

    INTEGER :: get_category_idx
    INTEGER :: i

    CHARACTER(LEN=*), INTENT(IN) :: category, list

    get_category_idx = 0

    IF (TRIM(list) == 'tp') THEN
      DO i = 1, SIZE(oem_config%tp_category)
        IF(category == oem_config%tp_category(i)) THEN
          get_category_idx = i
        END IF
      END DO
    ELSE
!    IF (TRIM(list) == 'vp') THEN
      DO i = 1, SIZE(oem_config%vp_category)
        IF(category == oem_config%vp_category(i)) THEN
          get_category_idx = i
        END IF
      END DO
    ENDIF

  END FUNCTION get_category_idx


!==============================================================================
!==============================================================================
!+ Get the country index for the temporal profile
!------------------------------------------------------------------------------

  FUNCTION get_country_idx(emissions_id)
    IMPLICIT NONE

    INTEGER :: get_country_idx
    INTEGER :: i, emissions_id

    get_country_idx = 0

      DO i = 1, oem_config%tp_ncountry
        IF(emissions_id == oem_data%tp_countryid(i)) THEN
          get_country_idx = i
        END IF
      END DO

  END FUNCTION get_country_idx

!==============================================================================
!==============================================================================
!+ Get the hour from the date
!------------------------------------------------------------------------------

  FUNCTION extract_hour(date)
    IMPLICIT NONE

    INTEGER :: extract_hour
    CHARACTER(LEN=MAX_DATETIME_STR_LEN), INTENT(IN) :: date
    !format: 2017-01-12T05:57:00.000

    read(date(12:13),'(I2)') extract_hour
    extract_hour = extract_hour + 1

  END FUNCTION extract_hour

!==============================================================================
!==============================================================================
!+ Get the hour from the date
!------------------------------------------------------------------------------

  FUNCTION extract_minutes(date)
    IMPLICIT NONE

    INTEGER :: extract_minutes
    CHARACTER(LEN=MAX_DATETIME_STR_LEN), INTENT(IN) :: date
    !format: 2017-01-12T05:57:00.000

    read(date(15:16),'(I2)') extract_minutes

  END FUNCTION extract_minutes
  
!==============================================================================
!==============================================================================
!+ Extract week day from a date
!------------------------------------------------------------------------------

  FUNCTION extract_day_of_week(jd2)

!-----------------------------------------------------------------------------
! Description: Extract week day from a date formatted as "wd   dd.mm.yy  hh mm ss UTC"
! Return value is an integer between 1 and 7
!-----------------------------------------------------------------------------

  IMPLICIT NONE

   INTEGER :: extract_day_of_week
   INTEGER :: i
   REAL(KIND=wp) :: y
   TYPE(julianday), POINTER :: jd2 ! julidan day: [day ms]

   y = REAL(jd2%day,wp) + (REAL(jd2%ms,wp)/REAL(no_of_ms_in_a_day,wp))

   extract_day_of_week = int(MOD(y+0.5_wp,7._wp)+1._wp)

  END FUNCTION extract_day_of_week

!==============================================================================
!==============================================================================
!+ Extract month from the date
!------------------------------------------------------------------------------

   FUNCTION extract_month(date)

!------------------------------------------------------------------------------
! Description: Extract month from a date formatted as "yyyymmddhhmmss"
! Return value is a integer between 1 and 12
!------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: extract_month
    CHARACTER(LEN=MAX_DATETIME_STR_LEN), INTENT(IN) :: date
    !format: 2017-01-12T05:57:00.000

    read(date(6:7),'(I2)') extract_month

  END FUNCTION extract_month


!==============================================================================


END MODULE mo_art_oem_emission

