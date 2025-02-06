!
! Source Module for OEM: Online Anthropogenic Emissions
!--------------------------------------------------------------------------------------
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

MODULE mo_art_oem_emission

  !------------------------------------------------------------------------------------
  ! Description:
  !   This module contains subroutines for computatoin of the emissions by 
  !   scaling the gridded emissions with temporal and vertical profiles and
  !   adds them to the tracers.
  
  !! Modifications: 
  !! 2024: Arash Hamzehloo, Empa
  !! - OEM was substantially refactored & ported to GPUs.
  !!
  !====================================================================================
    ! ICON 
    USE mo_parallel_config,        ONLY: nproma
    USE mo_kind,                   ONLY: wp, i8
    USE mo_model_domain,           ONLY: t_patch
    USE mo_nonhydro_state,         ONLY: p_nh_state
  
    USE mo_time_config,            ONLY: time_config
    USE mtime,                     ONLY: MAX_DATETIME_STR_LEN,      &
                                     &   datetimeToString,          &
                                     &   julianday,                 &
                                     &   newJulianday,              &
                                     &   getJulianDayFromDatetime,  &
                                     &   getnoofdaysinyeardatetime, &
                                     &   datetime,                  &
                                     &   timedelta,                 &
                                     &   newTimedelta,              &
                                     &   newDatetime,               &
                                     &   no_of_ms_in_a_day,         &
                                     &   OPERATOR(+)
  
    ! ART
    USE mo_art_atmo_data,          ONLY: t_art_atmo
    USE mo_art_data,               ONLY: p_art_data
    USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
  
    ! OEM
    USE mo_art_oem_types,          ONLY: p_art_oem_data,     &
                                     &   t_art_oem_data,     &
                                     &   t_art_oem_config
  
  !---------------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    PRIVATE
    PUBLIC ::                       &
      &  art_oem_compute_emissions, &
      &  art_oem_extract_time_information
  
    ! Constant variable
    INTEGER,  PARAMETER :: tp_param_hourofday = 24
    INTEGER,  PARAMETER :: tp_param_dayofweek = 7
    INTEGER,  PARAMETER :: tp_param_monthofyear = 12
    INTEGER(KIND=2), PARAMETER :: tp_param_hour = 8784
  
    TYPE(julianday), POINTER :: jd, jdref
  
    CHARACTER(LEN=3), DIMENSION(7) :: day_of_week = &
      & (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)
  
    INTEGER :: global_iteration, hour_iteration, hour_min

    INTEGER, allocatable :: hod(:), dow(:), moy(:), hoy(:)

    INTEGER, allocatable :: grd_index(:,:), vp_cat_idx(:,:), tp_cat_idx(:,:), tp_country_idx(:,:)
 
    REAL(KIND=wp), allocatable :: tsf_now(:,:,:,:), tsf_next(:,:,:,:), total_weight(:,:,:,:,:)

    INTEGER :: ncat_max = 200

  !====================================================================================
  ! Module procedures
  !====================================================================================
   
  CONTAINS
  
    SUBROUTINE art_oem_compute_emissions(p_tracer_now,p_patch,dtime,mtime_current,ierror,yerrmsg)
  
  !-----------------------------------------------------------------------------------
  ! Description: This subroutine uses the temporal and vertical profiles to compute
  ! the gridded emissions and add them to the OEM-tracers.
  !-----------------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    REAL(wp), INTENT(INOUT) :: &
     &  p_tracer_now(:,:,:,:)     !< tracer mixing ratio
  
    TYPE(t_patch), INTENT(IN) :: &
     &  p_patch                  !< patch on which computation is performed
  
    REAL(wp), INTENT(IN) :: &
     &  dtime                     !< time step    
  
    TYPE(datetime), POINTER ::  &
     &  mtime_current             !< current datetime
  
    INTEGER, INTENT(INOUT)             :: ierror
    CHARACTER(LEN= *), INTENT(INOUT)   :: yerrmsg
  
    !---------------------------------------------------------------------------------
    ! Local variables
    INTEGER ::                 &
      &  nc, nt, is, ie, jg,   &
      &  nlev, k, trcr_idx, jb, jc,          &
      &  nblks_c, i_startblk, i_endblk, idx, &
      &  ii, jj, kk, ll, &
      &  min, f_sec, &
      &  nr, nens, ens_int_idx, table_nr, ens_count, nt_old
  
    REAL(KIND=wp) ::                         &
      &  n1, n2, temp_scaling_fact_now,      &
      &  temp_scaling_fact_next, minutes, temp_scaling_fact, lambda
  
       
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo
  
    TYPE(datetime) :: datetime_next
    TYPE(timedelta), POINTER :: mtime_td

    CHARACTER(*), PARAMETER :: routine = "art_oem_compute_emissions"

  ! End of header
  !==================================================================================== 
    ierror = 0
    yerrmsg = '   '
    jg   = p_patch%id
    nblks_c = p_art_data(jg)%atmo%nblks
    nlev = p_art_data(jg)%atmo%nlev
    
  !------------------------------------------------------------------------------------
  
    IF (p_art_oem_data%configure%emis_tracer>0) THEN

      ! Read start- and end-block for this PE:
      i_startblk = p_art_data(jg)%atmo%i_startblk
      i_endblk   = p_art_data(jg)%atmo%i_endblk

      global_iteration = global_iteration + 1

      if (global_iteration==1) then
        f_sec=INT(dtime)
        allocate(hod(2), dow(2), moy(2), hoy(2))
      else
        f_sec=0
      endif 

      IF (mtime_current%time%minute==0 .AND. mtime_current%time%second == f_sec ) THEN

        ! Extract the different time information for this timestep
        CALL art_oem_extract_time_information(time_config%tc_current_date,hod(1),dow(1),moy(1),hoy(1),minutes)
        ! Extract the time information for one hour later
        mtime_td => newTimedelta("PT01H")
        datetime_next = time_config%tc_current_date + mtime_td
        CALL art_oem_extract_time_information(datetime_next,hod(2),dow(2),moy(2),hoy(2),minutes)

        hour_iteration = 0

      END IF

      hour_iteration = hour_iteration + 1
      hour_min = hour_iteration*INT(dtime/60._wp)

!------------------------------------------------------------------------------------
! Section 1:
!     
! The following section is executed only once at the beginning of the 
! simulation. Here, indices of the grid emission, vertical profile, temporal 
! profile, and countries are updated. 
!
! After that, the total "weight" due to the "gridded" emission and vertical 
! scaling are evaluated. 
!------------------------------------------------------------------------------------

      IF (global_iteration == 1) THEN
        allocate(grd_index(p_art_oem_data%configure%emis_tracer, ncat_max))
        allocate(vp_cat_idx(p_art_oem_data%configure%emis_tracer, ncat_max))
        allocate(tp_cat_idx(p_art_oem_data%configure%emis_tracer, ncat_max))

        DO nt = 1, p_art_oem_data%configure%emis_tracer 
          DO nc = 1, ncat_max

            IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN

              grd_index(nt,nc) = 0

              DO ii = 1, SIZE(p_art_oem_data%configure%gridded_emissions_idx)
                IF(p_art_oem_data%configure%ycatl_l(nt,nc) == p_art_oem_data%configure%gridded_emissions_idx(ii)) THEN
                  grd_index(nt,nc) = ii
                  EXIT
                END IF
              END DO

              vp_cat_idx(nt,nc) = 0

              DO jj = 1, SIZE(p_art_oem_data%configure%vp_category)
                IF(p_art_oem_data%configure%yvpl_l(nt,nc) == p_art_oem_data%configure%vp_category(jj)) THEN
                  vp_cat_idx(nt,nc) = jj
                  EXIT
                END IF
              END DO

              tp_cat_idx(nt,nc) = 0

              DO kk = 1, SIZE(p_art_oem_data%configure%tp_category)
                IF(p_art_oem_data%configure%ytpl_l(nt,nc) == p_art_oem_data%configure%tp_category(kk)) THEN
                  tp_cat_idx(nt,nc) = kk
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO

        DO jb = i_startblk, i_endblk
          is = 1

          IF (jb == i_startblk) THEN

            allocate(tp_country_idx(nproma, i_endblk))
            allocate(total_weight(nproma, nlev, p_art_oem_data%configure%emis_tracer,ncat_max, i_endblk))
            allocate(tsf_now(nproma, p_art_oem_data%configure%emis_tracer,ncat_max, i_endblk))
            allocate(tsf_next(nproma, p_art_oem_data%configure%emis_tracer,ncat_max, i_endblk))

            total_weight(:,:,:,:,:) = 0._wp  
            tsf_now(:,:,:,:) = 0._wp  
            tsf_next(:,:,:,:) = 0._wp
            lambda =  0._wp  

          ENDIF

          DO nt = 1, p_art_oem_data%configure%emis_tracer 
            DO nc = 1, ncat_max

              IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN

                DO k=1, nlev 
                  DO jc = is, nproma

                    total_weight(jc,k,nt,nc,jb) = dtime * p_art_oem_data%data_fields%gridded_emissions(jc,jb,grd_index(nt,nc)) &
                    & * p_art_oem_data%data_fields%vert_scaling_fact(jc,k,jb,vp_cat_idx(nt,nc))

                    tp_country_idx(jc,jb) = 0
                    DO ll = 1, p_art_oem_data%configure%tp_ncountry
                      IF(p_art_oem_data%data_fields%country_ids(jc,jb) == p_art_oem_data%data_fields%tp_countryid(ll)) THEN
                        tp_country_idx(jc,jb) = ll
                      END IF
                    END DO
                    
                  END DO

                END DO

              END IF

            END DO
          END DO

        END DO

        !$ACC ENTER DATA COPYIN(p_nh_state, p_art_oem_data, p_patch, p_art_data(jg), &  
        !$ACC           & p_art_oem_data%configure%ycatl_l, p_art_oem_data%configure%gridded_emissions_idx, &
        !$ACC           & p_art_oem_data%configure%vp_category, p_art_oem_data%configure%yvpl_l, & 
        !$ACC           & p_art_oem_data%configure%ytpl_l,  p_art_oem_data%configure%tp_category, p_art_oem_data%configure%itype_tscale_l, &
        !$ACC           & p_art_oem_data%data_fields%tp_hourofday, p_art_oem_data%data_fields%tp_dayofweek, p_art_oem_data%data_fields%tp_monthofyear, &
        !$ACC           & p_art_oem_data%configure%tp_ncountry, p_art_oem_data%data_fields%country_ids, p_art_oem_data%data_fields%tp_countryid, & 
        !$ACC           & p_art_oem_data%data_fields%tp_hourofyear, p_art_oem_data%data_fields%tp_hourofyear, &
        !$ACC           & p_art_oem_data%data_fields%gridded_emissions, p_art_oem_data%configure%emis_tracer, &
        !$ACC           & p_art_oem_data%configure%emis_idx, p_art_oem_data%data_fields%vert_scaling_fact, p_art_data(jg)%atmo%i_endblk, &
        !$ACC           & p_art_oem_data%ensemble%ens_name, p_art_oem_data%configure%emis_name, p_art_oem_data%data_fields%reg_map, &
        !$ACC           & p_art_oem_data%ensemble%ens_table, p_art_oem_data%data_fields%lambda_mat)


        !$ACC ENTER DATA COPYIN(tp_country_idx, tp_cat_idx, hod, dow, moy, hoy, &
        !$ACC           & p_tracer_now, total_weight, tsf_now, tsf_next, &
        !$ACC           & ncat_max)

      ENDIF ! global_iteration == 1

!------------------------------------------------------------------------------------
! End of Section 1
!------------------------------------------------------------------------------------
      !$ACC DATA COPYIN(hour_min, mtime_current)
!------------------------------------------------------------------------------------
! Section 2:
!     
! The following section is executed at one hour (physical time) intervals. 
! Here, the temporal scaling factors for every minute of the next hour are 
! evaluated.
!------------------------------------------------------------------------------------

      !$ACC PARALLEL DEFAULT(PRESENT)
      !
      ! Note that on Piz Daint using VECTOR_LENGTH(256) may enhance the performance.  

      IF (mtime_current%time%minute==0 .AND. mtime_current%time%second == f_sec ) THEN
        !$ACC LOOP SEQ
        DO jb = i_startblk, i_endblk
          !$ACC LOOP SEQ
          DO nc = 1, ncat_max
            !$ACC LOOP SEQ
            DO nt = 1, p_art_oem_data%configure%emis_tracer
              IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN

                is = 1

                !$ACC LOOP GANG VECTOR
                DO jc = is, nproma

                  SELECT CASE (p_art_oem_data%configure%itype_tscale_l(nt))
  
                  CASE(0)

                    tsf_now(jc,nt,nc,jb) =  1._wp
                    tsf_next(jc,nt,nc,jb) =  1._wp

                  CASE(1)

                    tsf_now(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofday(hod(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                    &                       p_art_oem_data%data_fields%tp_dayofweek(dow(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                    &                       p_art_oem_data%data_fields%tp_monthofyear(moy(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))

                    tsf_next(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofday(hod(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                    &                        p_art_oem_data%data_fields%tp_dayofweek(dow(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                    &                        p_art_oem_data%data_fields%tp_monthofyear(moy(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))

                  CASE(2)

                    tsf_now(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofyear(hoy(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))
                    tsf_next(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofyear(hoy(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))

                  END SELECT

                END DO ! jc = is, nproma
              END IF ! p_art_oem_data%configure%ycatl_l(nt,nc) /= ""
            END DO ! nt = 1, p_art_oem_data%configure%emis_tracer
          END DO ! nc = 1, ncat_max
        END DO ! jb = i_startblk, i_endblk 
        
      END IF ! mtime_current%time%minute==0 .AND. mtime_current%time%second == f_sec 
         
      !$ACC END PARALLEL
!------------------------------------------------------------------------------------
! End of Section 2
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Section 3:
!
! This section is exectuted at every time-step and includes the main loops over 
! the tracers. The final temporal scaling facatore is calculated using a linear 
! interpolation to the values obtained in section 2.
! 
! Please note that in the code provided, the locations of the 'gang' and 'vector',
! (and 'worker', if needed) directives are case-sensitive. The GPU porting approach 
! outlined below is optimised for scenarios involving a limited number of tracers  
! and categories, with a high number of ensemble members (typically exceeding 100).      
!------------------------------------------------------------------------------------
      !$ACC PARALLEL DEFAULT(PRESENT)
      !
      ! Note that on Piz Daint using VECTOR_LENGTH(256) may enhance the performance.

      nt_old = 0
      ens_count = 0
      n2 = (hour_min)/60._wp
      n1 = 1._wp - n2
      
      !$ACC LOOP SEQ
      DO jb = i_startblk, i_endblk   
        !$ACC LOOP SEQ      
        DO  nt = 1, p_art_oem_data%configure%emis_tracer  
          !$ACC LOOP SEQ
          DO nc = 1, ncat_max

            IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN
              trcr_idx = p_art_oem_data%configure%emis_idx(nt)

              is = 1
              
              ! Evaluate the ensemble tracers 
              IF ( ANY( p_art_oem_data%ensemble%ens_name==p_art_oem_data%configure%emis_name(nt) ) ) THEN
                IF (nt_old/=nt) THEN
                  ens_count = ens_count+1
                  nt_old = nt
                ENDIF
 
                !$ACC LOOP GANG COLLAPSE(2) PRIVATE(ens_int_idx) 
                DO nens = 1,SIZE(p_art_oem_data%data_fields%lambda_mat, dim=4)  
                  DO table_nr=1,200
                    IF (p_art_oem_data%ensemble%ens_table(1,table_nr)==nens .AND. p_art_oem_data%ensemble%ens_name(table_nr)==p_art_oem_data%configure%emis_name(nt)) THEN
                      
                      ens_int_idx = p_art_oem_data%ensemble%ens_table(2,table_nr)

                      !$ACC LOOP VECTOR COLLAPSE(2) PRIVATE(nr, lambda) 
                      DO k = 1,nlev
                        DO jc = is, nproma
                          nr = p_art_oem_data%data_fields%reg_map(jc,jb)
                          lambda = p_art_oem_data%data_fields%lambda_mat(ens_count,nc,nr,nens)
                          
                          p_tracer_now(jc,k,jb,ens_int_idx) = p_tracer_now(jc,k,jb,ens_int_idx) + total_weight(jc,k,nt,nc,jb) & 
                          &                               * (n1*tsf_now(jc,nt,nc,jb)+n2*tsf_next(jc,nt,nc,jb)) &
                          &                               * lambda &  
                          &                              / p_nh_state(jg)%diag%airmass_now(jc,k,jb)

                        

                        ENDDO ! jc = is, nproma 
                      ENDDO ! k = 1,nlev
                    ENDIF ! ens_table
                  ENDDO ! table_nr=1,200
                ENDDO ! nens = 1,SIZE(lambda_mat, dim=4)
              ENDIF ! lens(nt)==.TRUE.
 
              ! Evaluate the emission tracers 
              !$ACC LOOP GANG VECTOR COLLAPSE(2)
              DO k = 1,nlev 
                DO jc = is, nproma

                  p_tracer_now(jc,k,jb,trcr_idx) = p_tracer_now(jc,k,jb,trcr_idx) + total_weight(jc,k,nt,nc,jb) &
                  &                               * (n1*tsf_now(jc,nt,nc,jb)+n2*tsf_next(jc,nt,nc,jb)) &  
                  &                              / p_nh_state(jg)%diag%airmass_now(jc,k,jb)

                ENDDO 
              ENDDO

            ENDIF ! p_art_oem_data%configure%ycatl_l(nt,nc) /= ""

          ENDDO ! nc 
        ENDDO ! nt
      ENDDO ! jb = i_startblk, i_endblk

      !$ACC END PARALLEL
!------------------------------------------------------------------------------------
! End of Section 3
!------------------------------------------------------------------------------------

      !$ACC END DATA
    ENDIF ! p_art_oem_data%configure%emis_tracer>0

  !-----------------------------------------------------------------------------------
  ! End of the Subroutine
  !------------------------------------------------------------------------------------
  
    END SUBROUTINE art_oem_compute_emissions
  
  !====================================================================================
  
    SUBROUTINE art_oem_extract_time_information(date, hour_of_day, day_of_week, month_of_year, &
      &                                         hour_of_year, minutes)
  
      IMPLICIT NONE
  
      ! Parameters
      TYPE(datetime), INTENT(IN) :: date
      INTEGER, INTENT(INOUT) :: hour_of_day
      INTEGER, INTENT(INOUT) :: day_of_week
      INTEGER, INTENT(INOUT) :: month_of_year
      INTEGER, INTENT(INOUT) :: hour_of_year
      REAL(KIND=wp), INTENT(INOUT) :: minutes
  
      ! Local variables
      CHARACTER(LEN=14) :: yactdate1 ! yyyymmddhhmmss
      CHARACTER(LEN=28) :: yactdate2 ! wd   dd.mm.yy  hh mm ss UTC
      INTEGER :: nactday ! day of the year
      REAL(KIND=wp) :: y1, y2, acthour ! actual hour of the day
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: time_string
      INTEGER :: errno, ndays
      TYPE(datetime) :: refdate
  
      CALL datetimeToString(date, time_string)
      !format: 2017-01-12T05:57:00.000
  
      jd => newJulianday(0_i8, 0_i8)
    
      ! ref-date (01 Jan T00 of this year):
      refdate = date
      refdate%date%month = 1
      refdate%date%day = 1
      refdate%time%hour = 0
      ! Julian day of refdate:
      jdref => newJulianday(0_i8, 0_i8) 
  
      y1 = REAL(jd%day,wp) + REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)
      y2 = REAL(jdref%day,wp) + REAL(jdref%ms,wp)/REAL(no_of_ms_in_a_day,wp)
  
      read(time_string(12:13),'(I2)') hour_of_day
      hour_of_day = hour_of_day + 1
  
      day_of_week = int(MOD((REAL(jd%day,wp) + (REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)))+0.5_wp,7._wp)+1._wp)
  
      read(time_string(6:7),'(I2)') month_of_year
  
      hour_of_year = int(24._wp*(y1-y2)) !24 * (nactday-1) + hour_of_day ! 1 Jan, 00 UTC -> 1
  
      read(time_string(15:16),'(F2.0)') minutes
  
    END SUBROUTINE art_oem_extract_time_information
  
  !====================================================================================
  
  END MODULE mo_art_oem_emission
