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

! Namelist variables for assimilation schemes

MODULE mo_assimilation_nml

  USE mo_kind,                ONLY: wp,i4
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_namelist,            ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio 

  USE mo_impl_constants,      ONLY: max_dom
  
  USE mo_master_config,       ONLY: isRestart
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,   &
                            &       open_and_restore_namelist, close_tmpfile
  USE mo_assimilation_config, ONLY: assimilation_config !,n_noobs
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings
  USE mo_run_config,          ONLY: ldass_lhn, dtime

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_assimilation_namelist

  !---------------------------------------------------------------
  ! Namelist variables 
  !---------------------------------------------------------------
  
  LOGICAL                          ::           &
    llhn(max_dom)             ,& ! on/off switch for latent heat nudging (lhn)
    llhnverif(max_dom)        ,& ! on/off switch for verification against radar
    lhn_artif        ,& ! apply artificial profile
    lhn_filt         ,& ! vertical filtering of lhn t-increments
    lhn_relax        ,& ! horizontal filtering of lhn t-increments
    lhn_limit        ,& ! impose absolute limit on lhn t-increm. (abs_lhn_lim)
    lhn_limitp       ,& ! impose absolute limit on lhn t-increm. (abs_lhn_lim), restrict increments to Gaussian profile
    lhn_no_ttend     ,& ! do not apply the t increments (except for artif points), but only the q increments
    lhn_hum_adj      ,& ! apply a humidity adjustment along with t-increments
    lhn_incloud      ,& ! apply the LHN-scaling in cloudy layers only
    lhn_spqual       ,& ! switch for use of a spatial quality function
    lhn_black        ,& ! use blacklist for radar data
    lhn_qrs          ,& ! calculate the integrated precipitation flux
    lhn_logscale     ,& ! apply logarithmic scaling factors
    lhn_wweight      ,& ! apply a weighting with respect to the mean horizontal wind
    lhn_bright       ,& ! apply bright band detection
    lhn_diag         ,& ! produce more detailed diagnostic output during lhn
    lhn_artif_only   ,& ! apply only artificial temperature profile instead of applying modelled tt_lheat profile
    lhn_refbias         ! apply bias corretion of refernece precipitation

  INTEGER ::  &
    nlhn_start       ,& ! start of latent heat nudging period in timesteps
    nlhn_end         ,& ! end of latent heat nudging period in timesteps
    nlhnverif_start  ,& ! start of latent heat nudging period in timesteps
    nlhnverif_end    ,& ! end of latent heat nudging period in timesteps
    nlhn_relax       ,& ! number of interations of horizontal filtering
    nradar(max_dom)  ,& ! max. number of radar stations within input data
    lhn_updt_rule(max_dom)  ! Rule for updates of temperature/humidity/hydrometeors

  REAL (KIND=wp)                   ::           &
    lhn_coef          ,& ! factor for reduction of lhn t-increments
    lhn_dt_obs        ,& ! time step of input data in minutes
    abs_lhn_lim       ,& ! absolute limit for lhn t-increments (used if lhn_limit or lhn_limitp)
    fac_lhn_artif     ,& ! factor when artificial profile will applied
    fac_lhn_artif_tune ,& ! tuning factor of increments when artificial profile will applied
    fac_lhn_up        ,& ! limiting factor for upscaling of model heating profile
    fac_lhn_down      ,& ! limiting factor for downscaling model heating profile
    thres_lhn         ,& ! threshold of rain rates to be consinderd within lhn approach
    rqrsgmax          ,& ! ratio of maximum of qrsgflux, needed for reference precipitation
    tt_artif_max      ,& ! maximum of latent heat used as artificial profile
    zlev_artif_max    ,& ! altidude of maximum of artificial profile
    std_artif_max     ,& ! parameter to define vertical width of artifical temperature profile
    start_fadeout     ,& ! time relative to lhn_end, when the lhn coefficient in decreased toward 0
    ref_bias0         ,& ! starting ratio of model precipitation to reference precipitation used in LHN
    dtrefbias         ,& ! adaptation time scale for bias correction
    rttend            ,& ! ratio of temperature increment to be applied
    bbthres           ,& ! threshold of precipitation rate used in bright band detection
    hzerolim             ! limitation of hzerocl used in bright band detection

  CHARACTER (LEN=100)              ::           &
    radar_in             ,& ! directory for reading radar-files
    radardata_file(max_dom)       ,& ! filename of radar data
    blacklist_file(max_dom)       ,& ! filename of blacklist for radar data
    height_file(max_dom)             ! filename of radar beam heights

  LOGICAL :: dace_coupling
  INTEGER :: dace_time_ctrl(3)
  INTEGER :: dace_debug       ! Debugging level for DACE interface
  CHARACTER(LEN=255) :: &
    dace_output_file         ! filename for stdout redirection
  CHARACTER(LEN=255) :: &
    dace_namelist_file       ! filename of the file containing the dace namelist

! CHARACTER (LEN=12)               ::           &
!    noobs_date (n_noobs)    ! array of missing observations

  ! NAMELIST/assimilation_nml/  llhn         ,llhnverif                  ,           &
  NAMELIST/assimilation_nml/  nlhn_start   ,nlhn_end                   ,           &
                              nlhnverif_start ,nlhnverif_end           ,           &
                              lhn_coef, fac_lhn_up  ,fac_lhn_down      ,           &
                              thres_lhn    ,                                       &  ! noobs_date
                              rqrsgmax     , rttend                    ,           &
                              radar_in     , lhn_updt_rule             ,           &
                              lhn_black    ,blacklist_file             ,           &
                              lhn_artif    ,fac_lhn_artif, fac_lhn_artif_tune ,    &
                              lhn_filt     ,lhn_hum_adj, lhn_no_ttend  ,           &
                              lhn_limit    ,lhn_limitp  ,abs_lhn_lim   ,           &
                              lhn_relax    ,nlhn_relax                 ,           &
                              lhn_incloud  ,lhn_diag, lhn_qrs          ,           &
                              lhn_logscale ,lhn_wweight                ,           &
                              lhn_bright   ,bbthres ,hzerolim          ,           &
                              lhn_artif_only   ,                                   &
                              height_file  ,lhn_spqual                 ,           &
                              lhn_dt_obs   ,nradar, radardata_file     ,           &
                              tt_artif_max ,zlev_artif_max, std_artif_max,         &
                              start_fadeout,                                       &
                              lhn_refbias  ,ref_bias0 , dtrefbias      ,           &
                              dace_coupling ,dace_time_ctrl, dace_debug,           &
                              dace_output_file, dace_namelist_file
CONTAINS
  !>
  !!
  SUBROUTINE read_assimilation_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: iunit
    INTEGER :: jg
    CHARACTER(LEN=*),PARAMETER :: routine='mo_assimilation_nml:read_assimilation_namelist'

    !------------------------------------------------------------
    ! Set up the default values
    !------------------------------------------------------------
    dace_coupling      = .false.
    dace_time_ctrl     = 0
    dace_debug         = 0
    dace_output_file   = ""
    dace_namelist_file = 'namelist'

    llhn(:)               = ldass_lhn
    llhnverif(:)          = .TRUE.
    lhn_artif             = .TRUE.
    lhn_filt           = .TRUE.
    lhn_relax          = .FALSE.
    lhn_limit          = .TRUE.
    lhn_limitp         = .FALSE.
    lhn_hum_adj        = .TRUE.
    lhn_no_ttend       = .FALSE.
    lhn_spqual         = .FALSE.
    lhn_black          = .FALSE.
    lhn_incloud        = .TRUE.
    lhn_qrs            = .TRUE.
    lhn_logscale       = .TRUE.
    lhn_wweight        = .FALSE.
    lhn_diag           = .FALSE.
    lhn_artif_only     = .FALSE.
    lhn_bright         = .FALSE.
    lhn_refbias        = .FALSE.
    nlhn_start         = -9999
    nlhn_end           = -9999
    nlhnverif_start    = -9999
    nlhnverif_end      = -9999
    nradar(:)          = 200
    lhn_updt_rule(:)   = 0
    nlhn_relax         = 2_i4
    lhn_dt_obs         = 300.0_wp
    lhn_coef           = 1.0_wp
    abs_lhn_lim        = 50._wp / 3600._wp    ! max. change in heating: 4K/h
    fac_lhn_artif      = 5._wp
    fac_lhn_artif_tune = 1._wp
    fac_lhn_up         = 2.0_wp
    fac_lhn_down       = 1.0_wp / 2.0_wp
    thres_lhn          = 0.1_wp / 3600._wp
    radar_in           = './'
    radardata_file(:)     = 'radardata.nc'
    blacklist_file(:)     = 'radarblacklist.nc'
    height_file(:)        = 'radarheight.nc'
!    noobs_date  (:)    = '            '
    rqrsgmax           = 1.0_wp
    rttend             = 1.0_wp
    ref_bias0          = 1.0_wp
    dtrefbias          = dtime
    tt_artif_max       = 0.0015
    zlev_artif_max     = 1000.
    std_artif_max      = 4.
    start_fadeout      = 1.0
    bbthres            = 2.5
    hzerolim           = 2500.

    !------------------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above by 
    ! values in the restart file
    !------------------------------------------------------------------------
    IF (isRestart()) THEN
      funit = open_and_restore_namelist('assimilation_nml')
      READ(funit,NML=assimilation_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------------------
    ! Read user's (new) specifications. (Done so far by all MPI processes)
    !------------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('assimilation_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, assimilation_nml)  ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, assimilation_nml)                                      ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, assimilation_nml)  ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !-----------------------------------------------------
    ! Sanity check
    !-----------------------------------------------------

    !-----------------------------------------------------
    ! 4. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=assimilation_nml)
      CALL store_and_close_namelist(funit, 'assimilation_nml')
    ENDIF
    
    ! write the contents of the namelist to an ASCII file
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=assimilation_nml)

    !-----------------------------------------------------
    ! 5. Fill configuration state
    !-----------------------------------------------------

    DO jg= 1,max_dom
        assimilation_config(jg)%dace_coupling   = dace_coupling
        assimilation_config(jg)%dace_time_ctrl  = dace_time_ctrl
        assimilation_config(jg)%dace_debug      = dace_debug
        assimilation_config(jg)%dace_output_file= dace_output_file
        assimilation_config(jg)%dace_namelist_file = dace_namelist_file
        assimilation_config(jg)%llhn            = llhn(jg)
        assimilation_config(jg)%llhnverif       = llhnverif(jg)
        assimilation_config(jg)%lhn_artif       = lhn_artif
        assimilation_config(jg)%lhn_filt        = lhn_filt
        assimilation_config(jg)%lhn_relax       = lhn_relax
        assimilation_config(jg)%lhn_limit       = lhn_limit
        assimilation_config(jg)%lhn_limitp      = lhn_limitp
        assimilation_config(jg)%lhn_hum_adj     = lhn_hum_adj
        assimilation_config(jg)%lhn_no_ttend    = lhn_no_ttend
        assimilation_config(jg)%lhn_diag        = lhn_diag
        assimilation_config(jg)%lhn_qrs         = lhn_qrs
        assimilation_config(jg)%lhn_logscale    = lhn_logscale
        assimilation_config(jg)%lhn_wweight     = lhn_wweight
        assimilation_config(jg)%lhn_spqual      = lhn_spqual
        assimilation_config(jg)%lhn_black       = lhn_black
        assimilation_config(jg)%lhn_incloud     = lhn_incloud
        assimilation_config(jg)%lhn_artif_only  = lhn_artif_only
        assimilation_config(jg)%lhn_bright      = lhn_bright
        assimilation_config(jg)%lhn_refbias     = lhn_refbias
        assimilation_config(jg)%nlhn_start      = nlhn_start
        assimilation_config(jg)%nlhn_end        = nlhn_end
        assimilation_config(jg)%nlhnverif_start = nlhnverif_start
        assimilation_config(jg)%nlhnverif_end   = nlhnverif_end
        assimilation_config(jg)%nradar          = nradar(jg)
        assimilation_config(jg)%lhn_updt_rule   = lhn_updt_rule(jg)
        assimilation_config(jg)%nlhn_relax      = nlhn_relax
        assimilation_config(jg)%lhn_dt_obs      = lhn_dt_obs/60.
        assimilation_config(jg)%lhn_coef        = lhn_coef
        assimilation_config(jg)%abs_lhn_lim     = abs_lhn_lim
        assimilation_config(jg)%fac_lhn_artif   = fac_lhn_artif
        assimilation_config(jg)%fac_lhn_artif_tune = fac_lhn_artif_tune
        assimilation_config(jg)%fac_lhn_up      = fac_lhn_up
        assimilation_config(jg)%fac_lhn_down    = fac_lhn_down
        assimilation_config(jg)%thres_lhn       = thres_lhn
        assimilation_config(jg)%rqrsgmax        = rqrsgmax
        assimilation_config(jg)%rttend          = rttend  
        assimilation_config(jg)%ref_bias0       = ref_bias0
        assimilation_config(jg)%dtrefbias       = dtrefbias
        assimilation_config(jg)%tt_artif_max    = tt_artif_max
        assimilation_config(jg)%zlev_artif_max  = zlev_artif_max
        assimilation_config(jg)%std_artif_max   = std_artif_max
        assimilation_config(jg)%start_fadeout   = start_fadeout
        assimilation_config(jg)%bbthres         = bbthres
        assimilation_config(jg)%hzerolim        = hzerolim
        assimilation_config(jg)%radar_in        = radar_in
        assimilation_config(jg)%radardata_file  = radardata_file(jg)
        assimilation_config(jg)%blacklist_file  = blacklist_file(jg)
        assimilation_config(jg)%height_file     = height_file(jg)
    ENDDO 

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=assimilation_nml)                    
      CALL store_and_close_namelist(funit, 'assimilation_nml')             
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=assimilation_nml)


  END SUBROUTINE read_assimilation_namelist
  !-------------

END MODULE mo_assimilation_nml
