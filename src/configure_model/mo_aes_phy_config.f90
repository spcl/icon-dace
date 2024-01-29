! Configuration of the AES physics package.
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

MODULE mo_aes_phy_config

  USE mo_exception     ,ONLY: message, message_text, print_value, finish
  USE mo_kind          ,ONLY: wp
  USE mo_impl_constants,ONLY: max_dom

  USE mtime            ,ONLY: OPERATOR(<), OPERATOR(>), OPERATOR(==),                                                              &
       &                      t_datetime =>datetime , newDatetime , datetimeToString , max_datetime_str_len ,                      &
       &                      t_timedelta=>timedelta, newTimedelta, timedeltaToString, max_timedelta_str_len, deallocateTimeDelta, &
       &                      t_event    =>event    , newEvent    , eventGroup       , addEventToEventGroup ,                      &
       &                      getTotalMilliSecondsTimeDelta
  USE mo_event_manager ,ONLY: addEventGroup, getEventGroup, printEventGroup

  USE mo_master_config ,ONLY: experimentStartDate, experimentStopDate

  USE mo_vertical_coord_table ,ONLY: vct_a

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                    name    !< name for this unit

  PUBLIC ::         aes_phy_config     !< user specified configuration parameters
  PUBLIC ::         aes_phy_tc         !< derived time control (tc) parameters
  PUBLIC ::    init_aes_phy_config     !< allocate and initialize aes_phy_config
  PUBLIC ::    eval_aes_phy_config     !< evaluate aes_phy_config
  PUBLIC ::    eval_aes_phy_tc         !< evaluate aes_phy_tc
  PUBLIC ::   print_aes_phy_config     !< print out

  PUBLIC ::                  dt_zero   !< a zero (0 sec) time interval 'PT0S'

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_phy'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the AES physics package
  !!
  TYPE t_aes_phy_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! dynamics physics coupling
     !                      !  If negative tracer mass fractions are found
     !                      !  in the dynamics to physics interface, then ...
     INTEGER  :: iqneg_d2p  !  1: ... they are reported
     !                      !  2: ... they are set to zero
     !                      !  3: ... they are reported and set to zero
     !
     !                      !  If negative tracer mass fractions are found
     !                      !  in the physics to dynamics interface, then ...
     INTEGER  :: iqneg_p2d  !  1: ... they are reported
     !                      !  2: ... they are set to zero
     !                      !  3: ... they are reported and set to zero
     !
     ! time control of processes
     ! - sd = start date
     ! - ed = end date
     ! - dt = time interval
     !
     ! forcing control of processes
     ! - fc
     !   fc = 0: diagnostic, tendencies of the param. are not used in integration
     !   fc = 1: prognostic, tendencies of the param. are used to update the model state
     !
     ! atmospheric physics
     CHARACTER(len=max_timedelta_str_len) :: dt_rad  !< time  step of LW radiation
     CHARACTER(len=max_datetime_str_len ) :: sd_rad  !< start time of LW radiation
     CHARACTER(len=max_datetime_str_len ) :: ed_rad  !< end   time of LW radiation
     INTEGER                              :: fc_rht
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_vdf  !< time  step of vertical diffusion
     CHARACTER(len=max_datetime_str_len ) :: sd_vdf  !< start time of vertical diffusion
     CHARACTER(len=max_datetime_str_len ) :: ed_vdf  !< end   time of vertical diffusion
     INTEGER                              :: fc_vdf
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_mig  !< time  step of cloud microphysics (graupel)
     CHARACTER(len=max_datetime_str_len ) :: sd_mig  !< start time of cloud microphysics (graupel)
     CHARACTER(len=max_datetime_str_len ) :: ed_mig  !< end   time of cloud microphysics (graupel)
     INTEGER                              :: fc_mig
     INTEGER                              :: if_mig
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_two  !< time  step of cloud microphysics (twomom)
     CHARACTER(len=max_datetime_str_len ) :: sd_two  !< start time of cloud microphysics (twomom)
     CHARACTER(len=max_datetime_str_len ) :: ed_two  !< end   time of cloud microphysics (twomom)
     INTEGER                              :: fc_two
     !
     ! atmospheric chemistry
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_car  !< time  step of lin. Cariolle ozone chemistry
     CHARACTER(len=max_datetime_str_len ) :: sd_car  !< start time of lin. Cariolle ozone chemistry
     CHARACTER(len=max_datetime_str_len ) :: ed_car  !< end   time of lin. Cariolle ozone chemistry
     INTEGER                              :: fc_car
     !
     CHARACTER(len=max_timedelta_str_len) :: dt_art  !< time  step of ART chemistry
     CHARACTER(len=max_datetime_str_len ) :: sd_art  !< start time of ART chemistry
     CHARACTER(len=max_datetime_str_len ) :: ed_art  !< end   time of ART chemistry
     INTEGER                              :: fc_art
     !
     ! surface
     LOGICAL                              :: lsstice !< .true. for inst. 6hourly sst and ice (prelim)
     LOGICAL                              :: l2moment !< .true. for 2-moment microphysics scheme
     LOGICAL                              :: lmlo    !< .true. for mixed layer ocean
     LOGICAL                              :: lice    !< .true. for sea-ice temperature calculation
     LOGICAL                              :: ljsb    !< .true. for calculating the JSBACH land surface
     LOGICAL                              :: llake   !< .true. for using lakes in JSBACH
     LOGICAL                              :: lamip   !< .true. for AMIP simulations
     LOGICAL                              :: use_shflx_adjustment
     LOGICAL                              :: suppress_shflx_adjustment_over_ice
        !
     ! vertical range parameters
     REAL(wp)                             :: zmaxcloudy !< maximum height (m)   for cloud related computations
     INTEGER                              :: jks_cloudy !< vertical start index for cloud related computations
     !                                                  !  to be diagnosed from zmaxcloudy and vct_a(:)
     !
  END TYPE t_aes_phy_config


  TYPE t_aes_phy_tc
     !
     ! mtime datetime and time delta, and events
     ! -----------------------------------------
     !
     ! - sd_prc = start date
     ! - ed_prc = end date
     ! - dt_prc = time interval
     ! - ev_prd = event
     !
     ! - dt_prc_sec = time interval in seconds
     !
     ! atmospheric physics
     TYPE(t_datetime ), POINTER :: datetime
     TYPE(t_timedelta), POINTER :: dt_phy
     REAL(wp)                   :: dt_phy_sec
     !
     ! atmospheric processes
     !
     TYPE(t_timedelta), POINTER :: dt_rad
     TYPE(t_datetime ), POINTER :: sd_rad
     TYPE(t_datetime ), POINTER :: ed_rad
     TYPE(t_event    ), POINTER :: ev_rad
     REAL(wp)                   :: dt_rad_sec
     LOGICAL                    :: is_in_sd_ed_interval_rad
     LOGICAL                    :: is_active_rad
     !
     TYPE(t_timedelta), POINTER :: dt_vdf
     TYPE(t_datetime ), POINTER :: sd_vdf
     TYPE(t_datetime ), POINTER :: ed_vdf
     TYPE(t_event    ), POINTER :: ev_vdf
     REAL(wp)                   :: dt_vdf_sec
     LOGICAL                    :: is_in_sd_ed_interval_vdf
     LOGICAL                    :: is_active_vdf
     !
     TYPE(t_timedelta), POINTER :: dt_mig
     TYPE(t_datetime ), POINTER :: sd_mig
     TYPE(t_datetime ), POINTER :: ed_mig
     TYPE(t_event    ), POINTER :: ev_mig
     REAL(wp)                   :: dt_mig_sec
     LOGICAL                    :: is_in_sd_ed_interval_mig
     LOGICAL                    :: is_active_mig
     !
     TYPE(t_timedelta), POINTER :: dt_two
     TYPE(t_datetime ), POINTER :: sd_two
     TYPE(t_datetime ), POINTER :: ed_two
     TYPE(t_event    ), POINTER :: ev_two
     REAL(wp)                   :: dt_two_sec
     LOGICAL                    :: is_in_sd_ed_interval_two
     LOGICAL                    :: is_active_two
     !
     ! atmospheric chemistry
     !
     TYPE(t_timedelta), POINTER :: dt_car
     TYPE(t_datetime ), POINTER :: sd_car
     TYPE(t_datetime ), POINTER :: ed_car
     TYPE(t_event    ), POINTER :: ev_car
     REAL(wp)                   :: dt_car_sec
     LOGICAL                    :: is_in_sd_ed_interval_car
     LOGICAL                    :: is_active_car
     !
     TYPE(t_timedelta), POINTER :: dt_art
     TYPE(t_datetime ), POINTER :: sd_art
     TYPE(t_datetime ), POINTER :: ed_art
     TYPE(t_event    ), POINTER :: ev_art
     REAL(wp)                   :: dt_art_sec
     LOGICAL                    :: is_in_sd_ed_interval_art
     LOGICAL                    :: is_active_art
     !
  END TYPE t_aes_phy_tc

  !>
  !! Configuration/logicals/timecontrol state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_phy_config), TARGET :: aes_phy_config (max_dom)
  TYPE(t_aes_phy_tc)    , TARGET :: aes_phy_tc     (max_dom)
  
  !>
  !! Events and event group
  !!
  INTEGER                   :: aes_phy_events
  TYPE(eventGroup), POINTER :: aes_phy_event_group

  !>
  !! For convenience
  !!
  TYPE(t_timedelta) , POINTER :: dt_zero

CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_phy_config
    !
    dt_zero =>  newTimedelta ('PT0S')
    !
    ! AES physics configuration
    ! ---------------------------
    !
    ! dynamics physics coupling
    aes_phy_config(:)%iqneg_d2p  = 0
    aes_phy_config(:)%iqneg_p2d  = 0
    !
    ! time control parameters
    aes_phy_config(:)% dt_rad = ''
    aes_phy_config(:)% sd_rad = ''
    aes_phy_config(:)% ed_rad = ''
    aes_phy_config(:)% fc_rht = 1
    !
    aes_phy_config(:)% dt_vdf = ''
    aes_phy_config(:)% sd_vdf = ''
    aes_phy_config(:)% ed_vdf = ''
    aes_phy_config(:)% fc_vdf = 1
    !
    aes_phy_config(:)% dt_mig = ''
    aes_phy_config(:)% sd_mig = ''
    aes_phy_config(:)% ed_mig = ''
    aes_phy_config(:)% fc_mig = 1
    aes_phy_config(:)% if_mig = 1
    !
    aes_phy_config(:)% dt_two = ''
    aes_phy_config(:)% sd_two = ''
    aes_phy_config(:)% ed_two = ''
    aes_phy_config(:)% fc_two = 1
    !
    aes_phy_config(:)% dt_car = ''
    aes_phy_config(:)% sd_car = ''
    aes_phy_config(:)% ed_car = ''
    aes_phy_config(:)% fc_car = 1
    !
    aes_phy_config(:)% dt_art = ''
    aes_phy_config(:)% sd_art = ''
    aes_phy_config(:)% ed_art = ''
    aes_phy_config(:)% fc_art = 1
    !
    ! logical switches
    aes_phy_config(:)% ljsb  = .FALSE.
    aes_phy_config(:)% llake = .FALSE.
    aes_phy_config(:)% lamip = .FALSE.
    aes_phy_config(:)% l2moment  = .FALSE.
    aes_phy_config(:)% lmlo  = .FALSE.
    aes_phy_config(:)% lice  = .FALSE.
    !
    aes_phy_config(:)% lsstice          = .FALSE.
    !
    aes_phy_config(:)% use_shflx_adjustment = .FALSE.
    aes_phy_config(:)% suppress_shflx_adjustment_over_ice = .FALSE.
    !
    ! vertical range parameters
    aes_phy_config(:)% zmaxcloudy = 33000.0_wp
    !
  END SUBROUTINE init_aes_phy_config

  !----

  !>
  !! Check the aes_phy_config state
  !!
  SUBROUTINE eval_aes_phy_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg, jk
    CHARACTER(LEN=2)    :: cg
    !
    DO jg = 1,ng
       !
       ! time control of parameterizations
       !
       WRITE(cg,'(i0)') jg
       !
       CALL eval_aes_phy_config_details(TRIM(cg),                'rad' ,&
            &                             aes_phy_config (jg)% dt_rad  ,&
            &                             aes_phy_config (jg)% sd_rad  ,&
            &                             aes_phy_config (jg)% ed_rad  ,&
            &                             aes_phy_config (jg)% fc_rht  )
       !
       CALL eval_aes_phy_config_details(TRIM(cg),                'vdf' ,&
            &                             aes_phy_config (jg)% dt_vdf  ,&
            &                             aes_phy_config (jg)% sd_vdf  ,&
            &                             aes_phy_config (jg)% ed_vdf  ,&
            &                             aes_phy_config (jg)% fc_vdf  )
       !
       CALL eval_aes_phy_config_details(TRIM(cg),                'mig' ,&
            &                             aes_phy_config (jg)% dt_mig  ,&
            &                             aes_phy_config (jg)% sd_mig  ,&
            &                             aes_phy_config (jg)% ed_mig  ,&
            &                             aes_phy_config (jg)% fc_mig  )
       !
       CALL eval_aes_phy_config_details(TRIM(cg),                'two' ,&
            &                             aes_phy_config (jg)% dt_two  ,&
            &                             aes_phy_config (jg)% sd_two  ,&
            &                             aes_phy_config (jg)% ed_two  ,&
            &                             aes_phy_config (jg)% fc_two  )
       !
       CALL eval_aes_phy_config_details(TRIM(cg),                'car' ,&
            &                             aes_phy_config (jg)% dt_car  ,&
            &                             aes_phy_config (jg)% sd_car  ,&
            &                             aes_phy_config (jg)% ed_car  ,&
            &                             aes_phy_config (jg)% fc_car  )
       !
       CALL eval_aes_phy_config_details(TRIM(cg),                'art' ,&
            &                             aes_phy_config (jg)% dt_art  ,&
            &                             aes_phy_config (jg)% sd_art  ,&
            &                             aes_phy_config (jg)% ed_art  ,&
            &                             aes_phy_config (jg)% fc_art  )
       !
       ! vertical range for cloud related computations
       !
       aes_phy_config(jg)% jks_cloudy = 1
       DO jk = 1,SIZE(vct_a)-1
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > aes_phy_config(jg)% zmaxcloudy) THEN
             aes_phy_config(jg)% jks_cloudy = aes_phy_config(jg)% jks_cloudy + 1
          ELSE
             EXIT
          END IF
       END DO
       !
    END DO
    !
  CONTAINS
    !
    SUBROUTINE eval_aes_phy_config_details  (cg,process ,&
         &                                   config_dt  ,&
         &                                   config_sd  ,&
         &                                   config_ed  ,&
         &                                   config_fc  )
      !
      CHARACTER(LEN=*),PARAMETER  :: method_name ='eval_aes_phy_config_details'
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process
      !
      ! sd, ed and tc arguments are empty strings or 'P...'  strings
      CHARACTER(len=max_timedelta_str_len), INTENT(inout) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(inout) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(inout) :: config_ed
      !
      ! forcing control
      INTEGER                             , INTENT(in)    :: config_fc
      !
      ! mtime time control (TC) variables
      TYPE(t_timedelta), POINTER :: tc_dt
      !
      ! 1. if dt='' or dt contains only blanks, then use dt='PT0S',
      !    because MTIME cannot digest empty strings
      !
      IF (TRIM(config_dt)=='') THEN
         config_dt='PT0S'
      END IF
      !
      !
      ! 2. if dt<0 then stop
      !
      tc_dt => newTimedelta (config_dt)
      IF (tc_dt < dt_zero) THEN
         CALL finish(method_name,'negative aes_phy_config('//TRIM(cg)//')% dt_'//TRIM(process)//' is not allowed')
      END IF
      !
      !
      ! 3. if dt is zero in any format, then set 'PT0S'
      !
      IF (tc_dt == dt_zero) THEN
         config_dt = 'PT0S'
      END IF
      !
      !
      ! 4. if dt>0 check start and end dates
      !
      IF (tc_dt > dt_zero) THEN
         !
         ! if start and end dates are empty strings or contain only blanks
         ! then use the start and stop dates of the experiment
         !
         IF (TRIM(config_sd) == '') config_sd = experimentStartDate
         IF (TRIM(config_ed) == '') config_ed = experimentStopDate
         !
      END IF
      !
      !
      ! 5. fc must be 0 or 1
      !
      SELECT CASE(config_fc)
      CASE(0,1)
         ! OK
      CASE DEFAULT
         ! not allowed
         CALL finish(method_name,'aes_phy_config('//TRIM(cg)//')% fc_'//TRIM(process)//' must be 0 or 1')
      END SELECT
      !
      CALL deallocateTimeDelta(tc_dt)
      !
    END SUBROUTINE eval_aes_phy_config_details
    !
  END SUBROUTINE eval_aes_phy_config

  !----

  !>
  !! Evaluate the configuration state
  !!
  SUBROUTINE eval_aes_phy_tc(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    ! AES physics timecontrol
    ! -------------------------
    !
    ! mtime events
    !
    aes_phy_events        =  addEventGroup("aes_phy_events_group")
    aes_phy_event_group   => getEventGroup( aes_phy_events )
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL eval_aes_phy_tc_details(cg,                     'rad'    ,&
            &                         aes_phy_config(jg)% dt_rad     ,&
            &                         aes_phy_config(jg)% sd_rad     ,&
            &                         aes_phy_config(jg)% ed_rad     ,&
            &                         aes_phy_tc    (jg)% dt_rad     ,&
            &                         aes_phy_tc    (jg)% sd_rad     ,&
            &                         aes_phy_tc    (jg)% ed_rad     ,&
            &                         aes_phy_tc    (jg)% ev_rad     ,&
            &                         aes_phy_tc    (jg)% dt_rad_sec )
       !
       CALL eval_aes_phy_tc_details(cg,                     'vdf'    ,&
            &                         aes_phy_config(jg)% dt_vdf     ,&
            &                         aes_phy_config(jg)% sd_vdf     ,&
            &                         aes_phy_config(jg)% ed_vdf     ,&
            &                         aes_phy_tc    (jg)% dt_vdf     ,&
            &                         aes_phy_tc    (jg)% sd_vdf     ,&
            &                         aes_phy_tc    (jg)% ed_vdf     ,&
            &                         aes_phy_tc    (jg)% ev_vdf     ,&
            &                         aes_phy_tc    (jg)% dt_vdf_sec )
       !
       CALL eval_aes_phy_tc_details(cg,                     'mig'    ,&
            &                         aes_phy_config(jg)% dt_mig     ,&
            &                         aes_phy_config(jg)% sd_mig     ,&
            &                         aes_phy_config(jg)% ed_mig     ,&
            &                         aes_phy_tc    (jg)% dt_mig     ,&
            &                         aes_phy_tc    (jg)% sd_mig     ,&
            &                         aes_phy_tc    (jg)% ed_mig     ,&
            &                         aes_phy_tc    (jg)% ev_mig     ,&
            &                         aes_phy_tc    (jg)% dt_mig_sec )
       !
       CALL eval_aes_phy_tc_details(cg,                     'two'    ,&
            &                         aes_phy_config(jg)% dt_two     ,&
            &                         aes_phy_config(jg)% sd_two     ,&
            &                         aes_phy_config(jg)% ed_two     ,&
            &                         aes_phy_tc    (jg)% dt_two     ,&
            &                         aes_phy_tc    (jg)% sd_two     ,&
            &                         aes_phy_tc    (jg)% ed_two     ,&
            &                         aes_phy_tc    (jg)% ev_two     ,&
            &                         aes_phy_tc    (jg)% dt_two_sec )
       !
       CALL eval_aes_phy_tc_details(cg,                     'car'    ,&
            &                         aes_phy_config(jg)% dt_car     ,&
            &                         aes_phy_config(jg)% sd_car     ,&
            &                         aes_phy_config(jg)% ed_car     ,&
            &                         aes_phy_tc    (jg)% dt_car     ,&
            &                         aes_phy_tc    (jg)% sd_car     ,&
            &                         aes_phy_tc    (jg)% ed_car     ,&
            &                         aes_phy_tc    (jg)% ev_car     ,&
            &                         aes_phy_tc    (jg)% dt_car_sec )
       !
       CALL eval_aes_phy_tc_details(cg,                     'art'    ,&
            &                         aes_phy_config(jg)% dt_art     ,&
            &                         aes_phy_config(jg)% sd_art     ,&
            &                         aes_phy_config(jg)% ed_art     ,&
            &                         aes_phy_tc    (jg)% dt_art     ,&
            &                         aes_phy_tc    (jg)% sd_art     ,&
            &                         aes_phy_tc    (jg)% ed_art     ,&
            &                         aes_phy_tc    (jg)% ev_art     ,&
            &                         aes_phy_tc    (jg)% dt_art_sec )
       !
    END DO
    !
  CONTAINS
    !
    SUBROUTINE eval_aes_phy_tc_details  (cg, process ,&
         &                               config_dt   ,&
         &                               config_sd   ,&
         &                               config_ed   ,&
         &                               tc_dt       ,&
         &                               tc_sd       ,&
         &                               tc_ed       ,&
         &                               tc_ev       ,&
         &                               dt_sec      )
      !
      CHARACTER(LEN=*),PARAMETER  :: method_name ='eval_aes_phy_tc_details'
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in)    :: cg
      CHARACTER(len=*)                    , INTENT(in)    :: process
      !
      ! configuration strings
      CHARACTER(len=max_timedelta_str_len), INTENT(in) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_ed
      !
      ! mtime time control (TC) variables
      TYPE(t_timedelta), POINTER :: tc_dt
      TYPE(t_datetime ), POINTER :: tc_sd
      TYPE(t_datetime ), POINTER :: tc_ed
      TYPE(t_event    ), POINTER :: tc_ev
      !
      REAL(wp), INTENT(out) :: dt_sec

      LOGICAL :: lret
      !
      tc_dt => newTimedelta (config_dt)
      IF (tc_dt > dt_zero) THEN
         tc_sd => newDatetime  (config_sd)
         tc_ed => newDatetime  (config_ed)
         tc_ev => newEvent(process//'_d'//cg, &
              &            tc_sd,             & ! <- start date as reference date!
              &            tc_sd,             &
              &            tc_ed,             &
              &            tc_dt)
         lret   = addEventToEventGroup(tc_ev, aes_phy_event_group)
         IF (.NOT.lret) THEN
            CALL finish(method_name,'addEventToEventGroup returned .FALSE. for event '//process//'_d'//cg)
         END IF
         dt_sec = REAL(getTotalMilliSecondsTimeDelta(tc_dt,tc_sd),wp)/1000._wp
      END IF
      !
    END SUBROUTINE eval_aes_phy_tc_details
    !
  END SUBROUTINE eval_aes_phy_tc

  !----

  !>
  !! Print out the user controlled configuration state and the derived logicals and time controls
  !!
  SUBROUTINE print_aes_phy_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','AES physics configuration')
    CALL message    ('','=========================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL message    ('','User controlled parameters')
       CALL message    ('','..........................')
       CALL message    ('','')
       CALL message    ('','dynamics physics coupling')
       CALL message    ('','treatment of negative tracer mass fractions:')
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% iqneg_d2p  ',aes_phy_config(jg)% iqneg_d2p    )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% iqneg_p2d  ',aes_phy_config(jg)% iqneg_p2d    )
       CALL message    ('','')
       CALL message    ('','time control parameters')
       CALL message    ('','')
       CALL print_aes_phy_config_details(cg,                     'rad' ,&
            &                              aes_phy_config(jg)% dt_rad  ,&
            &                              aes_phy_config(jg)% sd_rad  ,&
            &                              aes_phy_config(jg)% ed_rad  ,&
            &                              aes_phy_config(jg)% fc_rht  )
       !
       CALL print_aes_phy_config_details(cg,                     'vdf' ,&
            &                              aes_phy_config(jg)% dt_vdf  ,&
            &                              aes_phy_config(jg)% sd_vdf  ,&
            &                              aes_phy_config(jg)% ed_vdf  ,&
            &                              aes_phy_config(jg)% fc_vdf  )
       !
       CALL print_aes_phy_config_details(cg,                     'mig' ,&
            &                              aes_phy_config(jg)% dt_mig  ,&
            &                              aes_phy_config(jg)% sd_mig  ,&
            &                              aes_phy_config(jg)% ed_mig  ,&
            &                              aes_phy_config(jg)% fc_mig  )
       !
       CALL print_aes_phy_config_details(cg,                     'two' ,&
            &                              aes_phy_config(jg)% dt_two  ,&
            &                              aes_phy_config(jg)% sd_two  ,&
            &                              aes_phy_config(jg)% ed_two  ,&
            &                              aes_phy_config(jg)% fc_two  )
       !
       CALL print_aes_phy_config_details(cg,                     'car' ,&
            &                              aes_phy_config(jg)% dt_car  ,&
            &                              aes_phy_config(jg)% sd_car  ,&
            &                              aes_phy_config(jg)% ed_car  ,&
            &                              aes_phy_config(jg)% fc_car  )
       !
       CALL print_aes_phy_config_details(cg,                     'art' ,&
            &                              aes_phy_config(jg)% dt_art  ,&
            &                              aes_phy_config(jg)% sd_art  ,&
            &                              aes_phy_config(jg)% ed_art  ,&
            &                              aes_phy_config(jg)% fc_art  )
       !
       CALL message    ('','logical switches')
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% lmlo ',    aes_phy_config(jg)% lmlo  )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% l2moment ',aes_phy_config(jg)% l2moment  )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% lice ',    aes_phy_config(jg)% lice  )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% ljsb ',    aes_phy_config(jg)% ljsb  )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% llake',    aes_phy_config(jg)% llake )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% lamip',    aes_phy_config(jg)% lamip )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% lsstice ', aes_phy_config(jg)% lsstice  )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% use_shflx_adjustment ', &
           &                                                            aes_phy_config(jg)% use_shflx_adjustment)
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% suppress_shflx_adjustment_over_ice ', &
           &                                                            aes_phy_config(jg)% suppress_shflx_adjustment_over_ice)
       CALL message    ('','')
       !
       CALL message    ('','vertical ranges')
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% zmaxcloudy ',aes_phy_config(jg)% zmaxcloudy )
       CALL print_value('    aes_phy_config('//TRIM(cg)//')% jks_cloudy ',aes_phy_config(jg)% jks_cloudy )
       CALL message    ('','')
       !
       CALL message    ('','Derived time control')
       CALL message    ('','....................')
       CALL message    ('','')
       !
       CALL print_aes_phy_tc_details(cg,                 'rad'    ,&
            &                          aes_phy_tc(jg)% dt_rad     ,&
            &                          aes_phy_tc(jg)% sd_rad     ,&
            &                          aes_phy_tc(jg)% ed_rad     ,&
            &                          aes_phy_tc(jg)% dt_rad_sec )
       !
       CALL print_aes_phy_tc_details(cg,                 'vdf'    ,&
            &                          aes_phy_tc(jg)% dt_vdf     ,&
            &                          aes_phy_tc(jg)% sd_vdf     ,&
            &                          aes_phy_tc(jg)% ed_vdf     ,&
            &                          aes_phy_tc(jg)% dt_vdf_sec )
       !
       CALL print_aes_phy_tc_details(cg,                 'mig'    ,&
            &                          aes_phy_tc(jg)% dt_mig     ,&
            &                          aes_phy_tc(jg)% sd_mig     ,&
            &                          aes_phy_tc(jg)% ed_mig     ,&
            &                          aes_phy_tc(jg)% dt_mig_sec )
       !
       CALL print_aes_phy_tc_details(cg,                 'two'    ,&
            &                          aes_phy_tc(jg)% dt_two     ,&
            &                          aes_phy_tc(jg)% sd_two     ,&
            &                          aes_phy_tc(jg)% ed_two     ,&
            &                          aes_phy_tc(jg)% dt_two_sec )
       !
       CALL print_aes_phy_tc_details(cg,                 'car'    ,&
            &                          aes_phy_tc(jg)% dt_car     ,&
            &                          aes_phy_tc(jg)% sd_car     ,&
            &                          aes_phy_tc(jg)% ed_car     ,&
            &                          aes_phy_tc(jg)% dt_car_sec )
       !
       CALL print_aes_phy_tc_details(cg,                 'art'    ,&
            &                          aes_phy_tc(jg)% dt_art     ,&
            &                          aes_phy_tc(jg)% sd_art     ,&
            &                          aes_phy_tc(jg)% ed_art     ,&
            &                          aes_phy_tc(jg)% dt_art_sec )
       !
       CALL message    ('','')
       CALL message    ('','------------------------------------------------------------------------')
       CALL message    ('','')
       !
    END DO
    !
    CALL message    ('','Events on all domains')
    CALL message    ('','.....................')
    CALL message    ('','')
    CALL printEventGroup(aes_phy_events)
    CALL message    ('','')
    CALL message    ('','========================================================================')
    !
  CONTAINS
    !
    SUBROUTINE print_aes_phy_config_details  (cg, process ,&
         &                                    config_dt   ,&
         &                                    config_sd   ,&
         &                                    config_ed   ,&
         &                                    config_fc   )
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*)                    , INTENT(in) :: cg
      CHARACTER(len=*)                    , INTENT(in) :: process
      !
      ! configuration strings
      CHARACTER(len=max_timedelta_str_len), INTENT(in) :: config_dt
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_sd
      CHARACTER(len=max_datetime_str_len ), INTENT(in) :: config_ed
      !
      ! forcing control
      INTEGER                             , INTENT(in) :: config_fc
      !
      CALL message       ('    aes_phy_config('//cg//')% dt_'//process,config_dt )
      IF (config_dt /= 'PT0S') THEN
         CALL message    ('    aes_phy_config('//cg//')% sd_'//process,config_sd )
         CALL message    ('    aes_phy_config('//cg//')% ed_'//process,config_ed )
         CALL print_value('    aes_phy_config('//cg//')% fc_'//process,config_fc )
         SELECT CASE(config_fc)
         CASE(0)
            CALL message ('',process//' tendencies are diagnostic')
         CASE(1)
            CALL message ('',process//' tendencies are used to update the model state')
         END SELECT
      END IF
      CALL message   ('','')
      !
    END SUBROUTINE print_aes_phy_config_details
    !
    !
    SUBROUTINE print_aes_phy_tc_details  (cg, process ,&
         &                                tc_dt       ,&
         &                                tc_sd       ,&
         &                                tc_ed       ,&
         &                                dt_sec      )
      !
      ! grid and name of evaluated configuration
      CHARACTER(len=*) , INTENT(in) :: cg
      CHARACTER(len=*) , INTENT(in) :: process
      !
      ! mtime time control (TC) variables
      TYPE(t_timedelta), POINTER    :: tc_dt
      TYPE(t_datetime ), POINTER    :: tc_sd
      TYPE(t_datetime ), POINTER    :: tc_ed
      !
      REAL(wp)         , INTENT(in) :: dt_sec
      !
      CHARACTER(LEN=MAX_TIMEDELTA_STR_LEN) :: td_string
      CHARACTER(LEN=MAX_DATETIME_STR_LEN ) :: dt_string
      !
      CALL timedeltaToString  (tc_dt, td_string)
      CALL message   ('    aes_phy_tc('//cg//')% dt_'//process,td_string)
      IF (tc_dt > dt_zero) THEN
         CALL datetimeToString(tc_sd, dt_string)
         CALL message('    aes_phy_tc('//cg//')% sd_'//process,dt_string)
         CALL datetimeToString(tc_ed, dt_string)
         CALL message('    aes_phy_tc('//cg//')% ed_'//process,dt_string)
         WRITE (message_text,'(f8.3,a)') dt_sec,' sec'
         CALL message('    aes_phy_tc('//cg//')% dt_'//process//'_sec',TRIM(message_text))
      END IF
      CALL message   ('','')
      !
    END SUBROUTINE print_aes_phy_tc_details
    !
  END SUBROUTINE print_aes_phy_config

  !----

END MODULE mo_aes_phy_config
