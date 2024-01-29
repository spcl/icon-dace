!
! mo_art_fplume_types
! This module initializes the data structure required by the FPLUME
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

MODULE mo_art_fplume_types
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_exception,                     ONLY: message, message_text,finish
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_parallel_config,               ONLY: nproma, p_test_run
  USE mo_physical_constants,            ONLY: amd
  USE mo_math_constants,                ONLY: pi
  USE mo_grid_config,                   ONLY: grid_sphere_radius
  USE mtime,                            ONLY: datetime,newDatetime,datetimeToString,   &
                                          &   max_datetime_str_len,newEvent,timedelta, &
                                          &   event, isCurrentEventActive
  USE mo_tracer_metadata_types,         ONLY: t_aero_meta, t_chem_meta
  USE mo_var_list,                      ONLY: get_tracer_info_dyn_by_idx
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic
  USE mo_gnat_gridsearch,               ONLY: gnat_init_grid, gnat_destroy,t_gnat_tree, &
                                          &   gnat_query_containing_triangles,gk,       &
                                          &   gnat_merge_distributed_queries
  USE mo_fortran_tools,                 ONLY: init
! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN, IART_XMLTAGLEN,                    &
                                          &   UNDEF_REAL_ART,UNDEF_INT_ART
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_art_config,                    ONLY: IART_PATH_LEN,art_config
  USE mo_art_fplume_utilities,          ONLY: get_input_rea


  IMPLICIT NONE

  PRIVATE

  TYPE t_fplume_init
    REAL(wp)                 :: &
      &  plume_Hdt,               & !< plume height
      &  plume_Mdt,               & !< MFR
      &  plume_udt,               & !< exit velocity
      &  plume_Tdt,               & !< exit temperature
      &  plume_Tvdt,              &
      &  plume_Tldt,              &
      &  plume_Tsdt,              &
      &  plume_wvdt=0.0_wp,       &
      &  plume_wldt=0.0_wp,       &
      &  plume_wsdt=0.0_wp,       &
      &  plume_fine_ash_fraction, &
      &  MER_SO2
    TYPE(datetime), POINTER :: &
      &  start_time,            & !< Start time of the emission event
      &  end_time                 !< End time of the emission event
    TYPE(event), POINTER    :: &
      &  emissEvent               !< Emission event
    CHARACTER(LEN=max_datetime_str_len) :: &
      & phase_beg,phase_end
    CONTAINS
      PROCEDURE, PUBLIC:: init => init_fplume_phase
      PROCEDURE, PUBLIC :: isActive => is_phase_active
  END TYPE t_fplume_init

  TYPE t_fplume_phases
    ! from .inp file
    TYPE(t_fplume_init), POINTER :: &
      &  phase(:)                        !< List of all phase
    INTEGER                     :: &
      &  nphases,                  & !< number of phase
      &  plume_modv,               &
      &  tri_iidx_loc,             & !< location idx of source
      &  tri_iblk_loc,             &  !< location blk of source
      &  ithis_nlocal_pts            !< Local (PE) volcanoes
    LOGICAL :: plume_moist_air      = .TRUE.
    LOGICAL :: plume_wind_coupling  = .TRUE.
    LOGICAL :: plume_reentrainment  = .TRUE.
    LOGICAL :: plume_latent_heat    = .TRUE.
    CHARACTER(LEN=32)         :: plume_solve_for
    REAL(wp)                 :: &
      &  fplume_min_height = 0.0_wp,& !< Minimum height above which FPlume is used (below Mastin)
      &  plume_n_MFR(2),        & !< MFR search range
      &  lon,                   & !< longitude of the emission source
      &  lat,                   & !< latitude of the emission source
      &  plume_zv,              & !< vent height
      &  Cw   = 2000.0_wp,      & !< water (generic) used if latent_heat=.false.
      &  Cp   = 1600.0_wp,      & !< solid particles (pyroclasts)
      &  Ca   = 1000.0_wp,      & !< air
      &  plume_xi  = 0.23_wp,   & !< Factor (Bursik 2001).
      &  plume_zmin_wind  = 0.0_wp,& !< Ignore wind entrainment below this zvalue (low jet region)
      &  plume_c_umbrella = 1.32_wp,&   !< Thickness of umbrella relative to Hb (>1)
      &  plume_a_s_jet    = 0.075_wp,&  !< Default (constant) value in jet region
      &  plume_a_s_plume  = 0.12_wp,&   !< Default (constant) value in plume region
      &  plume_a_v        = 0.3_wp,&    !< Default (constant) value
      &  cfactor_asha = 0.33_wp, &
      &  cfactor_ashb = 0.34_wp, &
      &  cfactor_ashc = 0.33_wp
    CHARACTER(LEN=32) ::        &
      &  plume_type_as,         &
      &  plume_type_av
    ! from .tgsd file
    INTEGER               :: nc
    REAL(wp), ALLOCATABLE :: fc    (:)                ! fc  (nc)
    REAL(wp), ALLOCATABLE :: rhop  (:)                ! rhop(nc)
    REAL(wp), ALLOCATABLE :: diam  (:)                ! diam(nc)
    REAL(wp), ALLOCATABLE :: sphe  (:)                ! sphe(nc)
    REAL(wp), ALLOCATABLE :: psi   (:)                ! psi (nc)
    REAL(wp)              :: vset_aggr
    REAL(wp)              :: diam_aggr
    REAL(wp)              :: Dfo
    LOGICAL :: plume_aggregation    = .FALSE.
    CHARACTER(LEN=32) :: plume_type_aggr
    CONTAINS
      PROCEDURE, PUBLIC :: init     => art_init_fplume_locs
  END TYPE t_fplume_phases

  TYPE t_art_all_volc_fplume
    TYPE(t_fplume_phases), POINTER :: &
      &  p(:)                        !< List of all point sources
  END TYPE t_art_all_volc_fplume

  PUBLIC :: t_fplume_init, t_fplume_phases, t_art_all_volc_fplume

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init_fplume_locs(this_volc,p_patch,lon,lat)
  CLASS(t_fplume_phases),INTENT(inout)     :: &
    &  this_volc                           !< Point source to be initialized
  TYPE(t_patch), INTENT(in)                :: &
    &  p_patch                        !< Current domain
  REAL(wp), INTENT(in)                    :: lon, lat
! Local variables
  TYPE(t_gnat_tree)         :: &
    &  gnat                      !< Grid search object
  INTEGER                   :: &
    &  nblks, npromz,          & !< calculated for call to GNAT
    &  gnat_jc,gnat_jb,        & !< coordinates for GNAT
    &  nlocsrc,istat
  INTEGER, ALLOCATABLE      :: &
    &  tri_idx(:,:,:),         &
    &  owner(:)                  !< rank of sender PE for each station (GNAT)
  REAL(gk), ALLOCATABLE     :: &
    &  in_points(:,:,:),       & !< geographical locations (GNAT)
    &  min_dist(:,:),          & !< minimal distance (GNAT)
    &  lat_idx(:,:),           &
    &  lon_idx(:,:)
  REAL(wp)                       :: &
    &  rvoid(2)
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: mess
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: thisroutine

  thisroutine = 'art_init_fplume'

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

!  lat = p_art_data%fplume_init_all%p(i_volc)%lat
!  lon = p_art_data%fplume_init_all%p(i_volc)%lon
  lat_idx(gnat_jc,gnat_jb) = lat
  lon_idx(gnat_jc,gnat_jb) = lon

  in_points(gnat_jc,gnat_jb,1) = lon_idx(gnat_jc,gnat_jb) * pi/180._wp
  in_points(gnat_jc,gnat_jb,2) = lat_idx(gnat_jc,gnat_jb) * pi/180._wp

  ! --- Build GNAT data structure
  CALL gnat_init_grid(gnat, p_patch)
  ! --- Perform proximity query
  CALL gnat_query_containing_triangles(gnat, p_patch, in_points(:,:,:),nproma, nblks, npromz,    &
    &                                  grid_sphere_radius,p_test_run,tri_idx(:,:,:),min_dist(:,:))
  CALL gnat_merge_distributed_queries(p_patch, 1, nproma, nblks, min_dist, tri_idx(:,:,:),       &
    &                                 in_points(:,:,:),owner(:),                                 &
    &                                 this_volc%ithis_nlocal_pts)

  ! --- Cleanup GNAT
  CALL gnat_destroy(gnat)

  gnat_jc=0
  gnat_jb=1


  DO nlocsrc = 1, this_volc%ithis_nlocal_pts
    gnat_jc=gnat_jc+1
    IF(gnat_jc>nproma) THEN
      gnat_jc=1
      gnat_jb=gnat_jb+1
    ENDIF
    this_volc%tri_iidx_loc=tri_idx(1,gnat_jc,gnat_jb)
    this_volc%tri_iblk_loc=tri_idx(2,gnat_jc,gnat_jb)
  ENDDO

  DEALLOCATE(tri_idx)
  DEALLOCATE(owner)
  DEALLOCATE(in_points)
  DEALLOCATE(min_dist)
  DEALLOCATE(lat_idx)
  DEALLOCATE(lon_idx)

  RETURN
END SUBROUTINE art_init_fplume_locs

SUBROUTINE init_fplume_phase(this_phase,phase_b,phase_e,phase_SO2,phase_Hdt,phase_Mdt,         &
                          &  phase_udt,phase_Tdt,phase_Tvdt,phase_Tldt,phase_Tsdt,phase_wvdt,  &
                          &  phase_wldt,phase_wsdt,tc_exp_refdate,tc_dt_model,                 &
                          &  phase_fine_ash_fraction)

  CLASS(t_fplume_init),INTENT(INOUT) :: &
    & this_phase
  CHARACTER(LEN=max_datetime_str_len):: phase_b, phase_e
  REAL(wp), INTENT(IN)    :: &
    &  phase_Mdt,              &          ! MER (kg/s) of each eruptionphase
    &  phase_Hdt,              &          ! Height (m agl) of each eruption phase
    &  phase_udt,              &          ! Exit velocity
    &  phase_Tdt,              &          ! Exit temperature
    &  phase_Tvdt,             &          ! Exit vapour temperature
    &  phase_Tldt,             &          ! Exit liquid water temperature
    &  phase_Tsdt,             &          ! Exit solid water temperature
    &  phase_wvdt,             &          ! Exit water vapour fraction
    &  phase_wldt,             &          ! Exit Water liquid fraction
    &  phase_wsdt,             &          ! Exit water solid  fraction
    &  phase_fine_ash_fraction,&          ! fine ash fraction
    &  phase_SO2
  TYPE(timedelta), POINTER, INTENT(IN) :: &
    &  tc_dt_model                    !< Model timestep
  TYPE(datetime), POINTER, INTENT(IN)  :: &
    &  tc_exp_refdate                 !< Experiment reference date configuration
  INTEGER:: ierror
  CHARACTER(LEN=MAX_CHAR_LENGTH)::thisroutine

  thisroutine = 'init_fplume_phase'

  this_phase%phase_beg = TRIM(phase_b)
  this_phase%phase_end = TRIM(phase_e)
  this_phase%plume_Hdt = phase_Hdt
  this_phase%plume_Mdt = phase_Mdt
  this_phase%plume_udt = phase_udt
  this_phase%plume_Tdt = phase_Tdt
  this_phase%plume_Tvdt = phase_Tvdt
  this_phase%plume_Tldt = phase_Tldt
  this_phase%plume_Tsdt = phase_Tsdt
  this_phase%plume_wvdt = phase_wvdt
  this_phase%plume_wldt = phase_wldt
  this_phase%plume_wsdt = phase_wsdt
  this_phase%plume_fine_ash_fraction = phase_fine_ash_fraction
  this_phase%MER_SO2  = phase_SO2

  ! ----------------------------------
  ! --- Time handling
  ! ----------------------------------

  ALLOCATE(this_phase%start_time)
  ALLOCATE(this_phase%end_time)
  ALLOCATE(this_phase%emissEvent)

  this_phase%start_time => newDatetime(this_phase%phase_beg,errno=ierror)
  IF(ierror /= SUCCESS) CALL finish(TRIM(thisroutine)//':init_fplume1',        &
                  & 'Could not create datetime object from string '           &
                  & //TRIM(this_phase%phase_beg)//'.')
  this_phase%end_time   => newDatetime(this_phase%phase_end,errno=ierror)
  IF(ierror /= SUCCESS) CALL finish(TRIM(thisroutine)//':init_fplume2',        &
                  & 'Could not create datetime object from string '           &
                  & //TRIM(this_phase%phase_end)//'.')

  this_phase%emissEvent => newEvent('FPLUME',tc_exp_refdate,                    &
                  & this_phase%start_time, this_phase%end_time, tc_dt_model)


RETURN
END SUBROUTINE init_fplume_phase

LOGICAL FUNCTION is_phase_active(this_phase, current_date)
!<
! FUNCTION is_phase_active
! Returns if the phase is currently active
! Based on: mo_art_pntSrc_types
! Part of Module: mo_art_fplume_types
! Author: Julia Bruckert, KIT
! Initial Release: 2020-02-14
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  CLASS(t_fplume_init),INTENT(INOUT)   :: &
    &  this_phase                          !< Fplume phase
  TYPE(datetime), POINTER, INTENT(IN) :: &
    &  current_date                        !< Current date
! Local variables


  IF (isCurrentEventActive(this_phase%emissEvent, current_date)) THEN
    is_phase_active = .TRUE.
  ELSE
    is_phase_active = .FALSE.
  ENDIF
END FUNCTION is_phase_active

END MODULE mo_art_fplume_types
