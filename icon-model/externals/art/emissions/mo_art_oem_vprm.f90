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

MODULE mo_art_oem_vprm

  !------------------------------------------------------------------------------
  !
  ! Description:
  !   This module contains subroutines for the reading in of gridded vegetation 
  !   indices (LSWI and EVI) and computation of the biospheric fluxes for the 
  !   Online Emission Module (OEM).
  !==============================================================================
    ! ICON
    USE mo_kind,                   ONLY: wp, i8
    USE mo_exception,              ONLY: message, message_text
    USE mo_var_list,               ONLY: t_var_list_ptr
    USE mo_model_domain,           ONLY: t_patch, p_patch
    USE mo_nonhydro_state,         ONLY: p_nh_state
    USE mo_parallel_config,        ONLY: nproma, idx_1d, blk_no,     &
                                     &   idx_no
    USE mo_math_constants,         ONLY: rad2deg
  
    USE mo_time_config,            ONLY: time_config !configure_time
    USE mtime,                     ONLY: MAX_DATETIME_STR_LEN, &
                                     &   datetimeToString,     &
                                     &   julianday,            &
                                     &   newJulianday,         &
                                     &   getJulianDayFromDatetime,  &
                                     &   getDatetimeFromJulianDay,  &
                                     &   datetime,             &
                                     &   no_of_ms_in_a_day,    &
                                     &   timedelta,            &
                                     &   newTimedelta,         &
                                     &   OPERATOR(+)
  
    USE mo_nwp_phy_state,          ONLY: prm_diag
  
    ! ART
    USE mo_art_atmo_data,          ONLY: t_art_atmo
    USE mo_art_data,               ONLY: p_art_data
    USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
  
    ! OEM
    USE mo_art_oem_types,          ONLY: p_art_oem_data,          &
                                     &   t_art_oem_data,          &
                                     &   t_art_oem_config
  
    USE mo_oem_config,             ONLY: vprm_par,          & 
                                     &   vprm_lambda,       & 
                                     &   vprm_alpha,        & 
                                     &   vprm_beta,         & 
                                     &   vprm_tmin,         & 
                                     &   vprm_tmax,         & 
                                     &   vprm_topt,         & 
                                     &   vprm_tlow,         & 
                                     &   lcut_area,         & 
                                     &   lon_cut_start,     &
                                     &   lon_cut_end,       &
                                     &   lat_cut_start,     & 
                                     &   lat_cut_end          
  
  !---------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    PRIVATE
    PUBLIC :: art_oem_compute_biosphere_fluxes, &
      &       art_oem_extract_dos
  
    TYPE(t_art_atmo),         POINTER :: art_atmo     !< ART atmo fields
    TYPE(t_art_oem_data),     POINTER :: oem_data     !< OEM data structure -> data
    TYPE(t_art_oem_config),   POINTER :: oem_config   !< OEM data structure -> config
  
    ! Constant variable
    INTEGER,  PARAMETER :: tp_param_hourofday = 24
    INTEGER,  PARAMETER :: tp_param_dayofweek = 7
    INTEGER,  PARAMETER :: tp_param_monthofyear = 12
    INTEGER(KIND=2), PARAMETER :: tp_param_hour = 8784
  
    TYPE(julianday), POINTER :: jd, jdref
  
    CHARACTER(LEN=3), DIMENSION(7) :: day_of_week = &
      & (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)
  
  
  !==============================================================================
  ! Module procedures
  !==============================================================================
  
  
  CONTAINS
  
  
  SUBROUTINE art_oem_compute_biosphere_fluxes(p_tracer_now,p_patch,dtime,mtime_current,ierror,yerrmsg)
  
  !-----------------------------------------------------------------------------
  ! Description: This subroutine computes the VPRM flux field and adds them to the OEM-tracer
  !-----------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    REAL(wp), INTENT(inout) :: &
     &  p_tracer_now(:,:,:,:)     !< tracer mixing ratio
  
    TYPE(t_patch), INTENT(IN) :: &
     &  p_patch                   !< patch on which computation is performed
  
    REAL(wp), INTENT(IN) :: &
     &  dtime                     !< time step    
  
    TYPE(datetime), POINTER ::  &
     &  mtime_current             !< current datetime
  
    INTEGER, INTENT(OUT)             :: ierror
    CHARACTER(LEN= *), INTENT(OUT)   :: yerrmsg
  
  
    !---------------------------------------------------------------------------
    ! Local variables
    INTEGER :: dos, is, ie, i_startblk, i_endblk, &
      &        jb, jc, jg, k, l, nr, nt, nblks_c, &
      &        table_nr, trcr_idx
  
    REAL(KIND=wp) :: z_mass, t_degc, newflux, veg_frac, a1, a2, a3, &
      &              t_scale, w_scale, p_scale, gee, evi_thresh, pabs
  
    REAL(KIND=wp), PARAMETER :: eps_div = 1.e-12_wp
  
    CHARACTER(LEN=2) :: numstring
  
    TYPE(datetime) :: datetime_next
  
    TYPE(timedelta), POINTER :: mtime_td
  
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo
  
    CHARACTER(*), PARAMETER :: routine = "art_oem_compute_biosphere_fluxes"
  
  
  
  !- End of header
  !==============================================================================
  
      jg   = p_patch%id
      art_atmo => p_art_data(jg)%atmo
  
      oem_data => p_art_oem_data%data_fields
      oem_config => p_art_oem_data%configure
  
      nblks_c = art_atmo%nblks
  
  !------------------------------------------------------------------------------
  ! Section 1: Calculate VPRM flux field
  !------------------------------------------------------------------------------
  
      IF (oem_config%vprm_tracer>0) THEN
  
        ! read start- and end-block for this PE:
        i_startblk = art_atmo%i_startblk
        i_endblk   = art_atmo%i_endblk
  
        ! Extract the day of simulation for this timestep
        CALL art_oem_extract_dos(time_config%tc_current_date,dos)
  
        DO nt = 1, oem_config%vprm_tracer
          trcr_idx = oem_config%vprm_idx(nt)
          DO jb = i_startblk, i_endblk
            ! read indices within this block:
            CALL art_get_indices_c(jg, jb, is, ie)
            DO jc = is, ie
              ! Initialize new flux for VPRM tracer
              newflux = 0.0_wp
              ! Get 2-meter temperature in degC
              t_degc = art_atmo%t_2m(jc,jb) - 273.15_wp
              ! Compute respiration 
              IF (oem_config%vprm_flux_type(nt) == 'resp' .OR.         &
                &  (oem_config%vprm_flux_type(nt) == 'resp_x' .AND. (  &
                &  art_atmo%lon(jc,jb)*rad2deg < lon_cut_start .OR.    &
                &  art_atmo%lon(jc,jb)*rad2deg > lon_cut_end   .OR.    &
                &  art_atmo%lat(jc,jb)*rad2deg < lat_cut_start .OR.    &
                &  art_atmo%lat(jc,jb)*rad2deg > lat_cut_end   )))     &
              & THEN
                DO l = 1, 8 
                  veg_frac = oem_data%vprm_lu_class_fraction(jc,jb,l)
                  IF (veg_frac < 1.E-8_wp) CYCLE ! Fluxes are zero
                  ! Respiration
                  newflux = newflux + (vprm_alpha(l) * t_degc + vprm_beta(l)) * veg_frac
                END DO
                ! Avoid negative values
                newflux = MAX(newflux, 0.0_wp)
                ! Convert units from umol m-2 s-1 to kg m-2 s-1
                newflux = newflux * 44.01_wp * 1.e-9_wp
            ! ! Compute photosynthetic uptake (GEE)
              END IF ! oem_config%vprm_flux_type(nt) == 'resp'
              IF (oem_config%vprm_flux_type(nt) == 'gpp' .OR.         &
                &  (oem_config%vprm_flux_type(nt) == 'gpp_x' .AND. (  &
                &  art_atmo%lon(jc,jb)*rad2deg < lon_cut_start .OR.   &
                &  art_atmo%lon(jc,jb)*rad2deg > lon_cut_end   .OR.   &
                &  art_atmo%lat(jc,jb)*rad2deg < lat_cut_start .OR.   &
                &  art_atmo%lat(jc,jb)*rad2deg > lat_cut_end   )))    &
              & THEN
                ! Initialize flux component gee for GPP tracer
                gee = 0.0_wp
                ! Radiation
                DO l = 1, 8 
                  a1 = t_degc - vprm_tmin(l)
                  a2 = t_degc - vprm_tmax(l)
                  a3 = t_degc - vprm_topt(l)
                  IF (a1 < 0._wp .OR. a2 > 0._wp) CYCLE ! No photosynthesis
    
                  ! Temperature sensitivity on photosynthesis
                  t_scale = (a1 * a2 / (a1 * a2 - a3**2 + eps_div))
                  IF (t_scale < 0._wp) CYCLE ! No photosynthesis
    
                  ! Effect of water stress on GEE
                  ! Modification due to different dependency on ground water
                  IF (l == oem_data%i_vprm_lc_shrub .OR. & ! Grassland and shrubland
                      l == oem_data%i_vprm_lc_grass)     & ! are xeric systems
                  THEN  
                    IF (oem_data%lswi_max(jc,jb) < 1.E-7_wp) THEN 
                      w_scale = 0._wp
                    ELSE
                      w_scale = (oem_data%lswi(jc,jb,dos) - oem_data%lswi_min(jc,jb)) / &
                        &       (oem_data%lswi_max(jc,jb) - oem_data%lswi_min(jc,jb) + eps_div)
                    END IF
                  ELSE
                    w_scale = (1._wp + oem_data%lswi(jc,jb,dos)) / (1._wp + oem_data%lswi_max(jc,jb))
                  END IF
       
                  ! Effect of leaf phenology
                  IF (l == oem_data%i_vprm_lc_evergreen) THEN          ! Evergreen
                    p_scale = 1._wp
                  ELSE IF (l == oem_data%i_vprm_lc_savanna .OR. &
                    &      l == oem_data%i_vprm_lc_grass) THEN        ! Savanna or grassland
                    p_scale = (1._wp + oem_data%lswi(jc,jb,dos)) / 2._wp
                  ELSE                                                ! Other vegetation types
                    evi_thresh = oem_data%evi_min(jc,jb) + &
                      &          0.55_wp * (oem_data%evi_max(jc,jb) - oem_data%evi_min(jc,jb))
                    IF (oem_data%evi(jc,jb,dos) >= evi_thresh) THEN  ! Full canopy period
                      p_scale = 1._wp
                    ELSE
                      p_scale = (1._wp + oem_data%lswi(jc,jb,dos)) / 2._wp  ! Bad-burst to full canopy period                  
                    END IF
                  END IF
                  ! VPRM vegetation class fraction
                  veg_frac = oem_data%vprm_lu_class_fraction(jc,jb,l)
                  ! Shortwave downward photosynthetically active flux at the surface [W/m2]
                  pabs = prm_diag(jg)%swflxsfc(jc,jb) + prm_diag(jg)%swflx_up_sfc(jc,jb)
                  gee = gee + (vprm_lambda(l) * t_scale * p_scale * w_scale *  &
                    &          oem_data%evi(jc,jb,dos) * 1._wp / (1._wp + pabs/vprm_par(l)) * pabs) &
                    &          * veg_frac
                END DO
                ! Switch the sign and avoid negative values
                gee = -1.0_wp * gee
                gee = MAX(gee, 0.0_wp)
                ! Convert units from umol m-2 s-1 to kg m-2 s-1
                newflux = newflux + gee * 44.01_wp * 1.e-9_wp
              END IF ! oem_config%vprm_flux_type(nt) == 'gpp'
              ! Calculate air mass
              z_mass = p_nh_state(jg)%diag%airmass_now(jc,art_atmo%nlev,jb)
              ! Add flux to tracer
              p_tracer_now(jc,art_atmo%nlev,jb,trcr_idx) = p_tracer_now(jc,art_atmo%nlev,jb,trcr_idx) + &
                 &                             dtime * newflux /z_mass
            ENDDO ! jc = is, ie
          ENDDO ! jb = i_startblk, i_endblk
        ENDDO ! nt = 1, oem_config%vprm_tracer
      ENDIF ! oem_config%vprm_tracer>0
  
  
  !------------------------------------------------------------------------------
  ! End of the Subroutine
  !------------------------------------------------------------------------------
  
    END SUBROUTINE art_oem_compute_biosphere_fluxes
  
  !==============================================================================
  !==============================================================================
  !+ Extract which day of the simulation we are at [1, 2, 3, ...]
  !------------------------------------------------------------------------------
  
    SUBROUTINE art_oem_extract_dos(date, day_of_year)
  
      IMPLICIT NONE
  
      ! Parameters
      TYPE(datetime), INTENT(IN) :: date
      INTEGER, INTENT(OUT) :: day_of_year
  
      ! Local variables
      REAL(KIND=wp) :: y1, y2
      CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: time_string
      INTEGER :: errno, hour_of_year
      TYPE(datetime) :: refdate
  
  
      CALL datetimeToString(date, time_string)
      !format: 2017-01-12T05:57:00.000
  
      ! Get current date
      jd => newJulianday(0_i8, 0_i8)
      CALL getJulianDayFromDatetime(date,jd,errno) !(dt, jd, errno)
    
      refdate = time_config%tc_startdate
      ! Julian day of refdate:
      jdref => newJulianday(0_i8, 0_i8)
      CALL getJulianDayFromDatetime(refdate,jdref,errno) !(dt, jd, errno)
  
      y1 = REAL(jd%day,wp) + REAL(jd%ms,wp)/REAL(no_of_ms_in_a_day,wp)
      y2 = REAL(jdref%day,wp) + REAL(jdref%ms,wp)/REAL(no_of_ms_in_a_day,wp)
      hour_of_year = int(24._wp*(y1-y2)) !24 * (nactday-1) + hour_of_day ! 1 Jan, 00 UTC -> 1
      day_of_year = hour_of_year/24_i8 + 1
  
    END SUBROUTINE art_oem_extract_dos
  !==============================================================================
  !==============================================================================
  
  END MODULE mo_art_oem_vprm
  
  
