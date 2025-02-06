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

! Configuration setup for SPPT (Stochastic Pertubation of Physics Tendencies)

MODULE mo_sppt_config

  USE mo_kind,                    ONLY: wp,i8
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_impl_constants,          ONLY: max_dom, SUCCESS, min_rlcell_int
  USE mo_time_config,             ONLY: time_config
  USE mo_run_config,              ONLY: msg_level
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config
  USE mtime,                      ONLY: datetime, timedelta, event, getPTStringFromMS, &
    &                                   MAX_TIMEDELTA_STR_LEN, newTimedelta, newEvent
  USE mo_util_mtime,              ONLY: mtime_timedelta_from_fseconds
  USE mo_physical_constants,      ONLY: p0ref, rd, grav
  USE mo_vertical_coord_table,    ONLY: vct_a
  USE mo_model_domain,            ONLY: t_patch
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_sync,                    ONLY: global_max, global_min
  USE mo_grid_config,             ONLY: n_dom, l_limited_area
#ifndef __NO_ICON_LES__
  USE mo_ls_forcing_nml,          ONLY: is_ls_forcing
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_sppt_config'

  ! FUNCTIONS/SUBROUTINES
  PUBLIC :: configure_sppt
  PUBLIC :: crosscheck_sppt

  ! Types:
  PUBLIC :: t_sppt_config

  ! VARIABLES
  PUBLIC :: sppt_config


  !---------------------------------------------------------------------------
  !
  ! Basic configuration setup for SPPT
  !
  !---------------------------------------------------------------------------

  TYPE:: t_sppt_config

    ! namelist variables

    LOGICAL :: lsppt      ! forecast with SPPT

    REAL(wp) :: hinc_rn   ! time increment (s) for drawing a new field of random numbers
    REAL(wp) :: dlat_rn   ! random number coarse grid point distance in meridional direction (deg)
    REAL(wp) :: dlon_rn   ! random number coarse grid point distance in zonal direction (deg)
    REAL(wp) :: range_rn  ! max magnitude of random numbers
    REAL(wp) :: stdv_rn   ! standard deviation of the gaussian distribution of random numbers

    REAL(wp), ALLOCATABLE :: taper(:)              ! < factor for vertical tapering

    TYPE(event), POINTER  :: read_rapa_Event   => NULL()

    ! derived variables
    !
    !< patch bounding box in geographical cordinates
    REAL(wp) :: bbmin_lon
    REAL(wp) :: bbmax_lon
    REAL(wp) :: bbmin_lat
    REAL(wp) :: bbmax_lat

    !< data for the coarse random field
    INTEGER :: coarse_nlat
    INTEGER :: coarse_nlon

    TYPE(timedelta):: mtime_hinc_rn    !  hinc_rn converted into type timedelta 

    TYPE(datetime) :: validity_date_rn_2d_new ! date at which the event read_rapa_event was active the last time

  END TYPE t_sppt_config

  TYPE(t_sppt_config) :: sppt_config(max_dom)


  CONTAINS

  !------------------------------------------------------------------------------------
  !
  ! Further configurations of SPPT (Stochastic Perturbation of Physics Tendencies)
  !
  !------------------------------------------------------------------------------------

  SUBROUTINE configure_sppt(n_dom, p_patch, mtime_current)

    ! Subroutine arguments
    INTEGER,        INTENT(IN)          :: n_dom          ! number of domains
    TYPE(t_patch),  INTENT(IN)          :: p_patch(:)
    TYPE(datetime), POINTER, INTENT(IN) :: mtime_current  !< current datetime (mtime)

    ! Local variables
    REAL(wp), ALLOCATABLE :: sigm_coord(:)
    REAL(wp) :: t00
    REAL(wp) :: p50, p100, zpno, zpnu   ! utility variables - naming convention from COSMO

    INTEGER  :: nlev
    INTEGER  :: jk, jg
    INTEGER  :: ist
    INTEGER  :: jb, jc
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end

    INTEGER  ::    k50        , & ! indexes for vertical tapering
                   k100       , & !
                   k870

    REAL(wp), PARAMETER                        ::  &
      &                   h_scal  = 10000.0_wp    , &
      &                   delta_t = 75.0_wp       , &
      &                   t0sl    = 288.150_wp

    CHARACTER(len=*), PARAMETER ::  &
      &      routine = modname//":configure_sppt"

    !$ACC ENTER DATA CREATE(sppt_config)

    DO jg = 1, n_dom

      IF (.NOT. sppt_config(jg)%lsppt) CYCLE

      nlev = p_patch(jg)%nlev

      !---------------------------------
      ! Create mtime event for SPPT
      !--------------------------------

      CALL setup_sppt_events(sppt_config(jg)%read_rapa_Event, sppt_config(jg)%hinc_rn)

      !---------------------------------
      ! Convert hinc_rn into mtime timedelta
      !--------------------------------

      CALL mtime_timedelta_from_fseconds(sppt_config(jg)%hinc_rn, mtime_current, sppt_config(jg)%mtime_hinc_rn)

      !---------------------------------
      ! Create factors for tapering
      !--------------------------------

      ! Allocate array for vertical tapering
      ALLOCATE(sppt_config(jg)%taper(nlev), stat=ist)
      IF(ist /= SUCCESS) CALL finish(routine, "memory allocation failure for array taper")

      ! Allocate local coordinate array
      ALLOCATE(sigm_coord(nlev+1), stat=ist)
      IF(ist /= SUCCESS) CALL finish(routine, "memory allocation failure for array sigm_coord")

      ! Initiate tapering factor
      sppt_config(jg)%taper(:) = 1.0_wp
      zpno = 0.0_wp


      ! Calculate sigm_coord
      t00   = t0sl - delta_t
      !
      DO jk=1,nlev+1
        sigm_coord(jk) = EXP( - grav/rd*h_scal/t00 *         &
                         LOG( (EXP(vct_a(jk)/h_scal)*t00 + delta_t)/(t00 + delta_t)) )
      ENDDO ! end of jk


      ! Find indices of vertical levels of interest
      DO jk = 1, nlev

        zpno = p0ref*sigm_coord(jk  )
        zpnu = p0ref*sigm_coord(jk+1)

        IF ( (zpno <= 870.0E2_wp) .AND. (870.0E2_wp < zpnu) ) THEN
          k870 = jk
        ENDIF

        IF ( (zpno <= 100.0E2_wp) .AND. (100.0E2_wp < zpnu) ) THEN
          k100 = jk
          p100 = zpno
        ENDIF

        IF ( (zpno <= 50.0E2_wp) .AND. (50.0E2_wp < zpnu) )  THEN
          k50 = jk
          p50 = zpno
        ENDIF

      ENDDO ! end of jk


      ! Calculate tapering factor
      DO jk = 1,nlev

        IF(jk <= k50) THEN
          sppt_config(jg)%taper(jk) = 0.0_wp
        ENDIF

        IF(jk > k50 .AND. jk < k100) THEN
          sppt_config(jg)%taper(jk) = (p0ref*0.5_wp*(sigm_coord(jk) + sigm_coord(jk+1))-p50)/(p100-p50)
       ENDIF

       IF(jk >= k100 .AND. jk <= k870) THEN
          sppt_config(jg)%taper(jk) = 1.0_wp
        ENDIF

        IF(jk > k870) THEN  ! Note one could also adapt the privious if-statement
          sppt_config(jg)%taper(jk) = 1.0_wp
        ENDIF

      ENDDO ! end of jk

      DEALLOCATE(sigm_coord, stat=ist)
      IF(ist /= SUCCESS) CALL finish(routine, "memory deallocation failure for array sigm_coord")

      !---------------------------------------
      ! Compute bounding box of given patch
      !---------------------------------------

      sppt_config(jg)%bbmin_lon = HUGE(0._wp)
      sppt_config(jg)%bbmin_lat = HUGE(0._wp)
      sppt_config(jg)%bbmax_lon = -1._wp*HUGE(0._wp)
      sppt_config(jg)%bbmax_lat = -1._wp*HUGE(0._wp)

      rl_start = 1
      rl_end   = min_rlcell_int

      i_startblk = p_patch(jg)%cells%start_block(rl_start)
      i_endblk   = p_patch(jg)%cells%end_block(rl_end)

      DO jb = i_startblk, i_endblk

        CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk,  &
                           i_startidx, i_endidx, rl_start, rl_end)

        DO jc = i_startidx, i_endidx
          sppt_config(jg)%bbmin_lon = MIN(sppt_config(jg)%bbmin_lon,p_patch(jg)%cells%center(jc,jb)%lon)
          sppt_config(jg)%bbmin_lat = MIN(sppt_config(jg)%bbmin_lat,p_patch(jg)%cells%center(jc,jb)%lat)
          sppt_config(jg)%bbmax_lon = MAX(sppt_config(jg)%bbmax_lon,p_patch(jg)%cells%center(jc,jb)%lon)
          sppt_config(jg)%bbmax_lat = MAX(sppt_config(jg)%bbmax_lat,p_patch(jg)%cells%center(jc,jb)%lat)
        ENDDO

      ENDDO  ! jb

      sppt_config(jg)%bbmin_lon = global_min(sppt_config(jg)%bbmin_lon)
      sppt_config(jg)%bbmin_lat = global_min(sppt_config(jg)%bbmin_lat)
      sppt_config(jg)%bbmax_lon = global_max(sppt_config(jg)%bbmax_lon)
      sppt_config(jg)%bbmax_lat = global_max(sppt_config(jg)%bbmax_lat)

      IF (msg_level >= 10) THEN
        WRITE(message_text,'(a,i2,a,2e16.8,a,2e16.8,a)') 'DOM',jg,' bounding box MAX(lat,lon) MIN(lat,lon): (', &
          &   sppt_config(jg)%bbmax_lat,sppt_config(jg)%bbmax_lon, ') (', &
          &   sppt_config(jg)%bbmin_lat,sppt_config(jg)%bbmin_lon, ')'
        CALL message(routine, message_text)
      ENDIF

      ! get dimension of coarse random data field
      ! there are ceil(range/rn) boxes, but ceil(range/rn)+1 vertices!
      sppt_config(jg)%coarse_nlon = ceiling((sppt_config(jg)%bbmax_lon - sppt_config(jg)%bbmin_lon) &
        &                           /sppt_config(jg)%dlon_rn)+1
      sppt_config(jg)%coarse_nlat = ceiling((sppt_config(jg)%bbmax_lat - sppt_config(jg)%bbmin_lat) &
        &                           /sppt_config(jg)%dlat_rn)+1

      !$ACC UPDATE DEVICE(sppt_config(jg:jg)) ASYNC(1)
      !$ACC ENTER DATA COPYIN(sppt_config(jg)%taper) ASYNC(1)

    ENDDO  ! jg

  END SUBROUTINE configure_sppt


  !------------------------------------------------------------------------------------
  !
  ! Setup mtime events for SPPT utilizing newEvent() subroutine
  !
  !   - read_rapa (read pattern/field of random numbers) 
  !
  !------------------------------------------------------------------------------------
  !
  SUBROUTINE setup_sppt_events(read_rapa_Event, hinc_rn)

    ! Input arguments
    TYPE(event), POINTER, INTENT(INOUT) :: read_rapa_Event ! and all consecutive events

    REAL(wp),             INTENT(IN)    :: hinc_rn         ! increment in seconds to read/create 
                                                           ! new random numbers (namelist parameter)

    ! Local
    TYPE(timedelta), POINTER               :: eventInterval    => NULL()
    CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN)   :: td_string

    !---------------------------------
    ! Create event read_rapa_Event for reading/creating random numbers for perturbation
    !--------------------------------

    CALL getPTStringFromMS(INT(hinc_rn*1000.0_wp,i8), td_string)

    eventInterval => newTimedelta(td_string)

    read_rapa_Event => newEvent( 'read_rapa', time_config%tc_exp_startdate,  &    ! "anchor date"
      &                           time_config%tc_exp_startdate,              &    ! start date
      &                           time_config%tc_exp_stopdate,               &    ! stop date
      &                           eventInterval)


  END SUBROUTINE setup_sppt_events


  !------------------------------------------------------------------------------------
  !
  ! Crosscheck for SPPT (Stochastic Perturbation of Physics Tendencies)
  !
  !------------------------------------------------------------------------------------
  !
  SUBROUTINE crosscheck_sppt()

    INTEGER  :: jg
    CHARACTER(len=*), PARAMETER :: routine =  modname//'::crosscheck_sppt'

    !--------------------------------
    ! Global cross checks - exit if ...
    !--------------------------------

    ! ... l_limited_area=.FALSE., i.e. global run
    IF(.NOT. l_limited_area) THEN
      CALL finish(routine, "Global SPPT runs are currently not supported.")
    ENDIF

#ifndef __NO_ICON_LES__
    ! ... large scale forcing is switched on
    IF(is_ls_forcing) THEN
        CALL finish(routine, "SPPT and large scale forcing not supported.")
      ENDIF
#endif

    !------------------------------
    ! Cross checks for all domains - exit if ...
    !------------------------------

    DO jg=1,n_dom

      ! Exit if higher order microphysic schemes are used
      IF(atm_phy_nwp_config(jg)%inwp_gscp > 2) THEN
        CALL finish(routine, "SPPT in combination with higher order microphysic schemes not supported/tested.")
      ENDIF

    END DO

  END SUBROUTINE crosscheck_sppt

END MODULE
