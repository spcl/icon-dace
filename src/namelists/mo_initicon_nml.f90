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

MODULE mo_initicon_nml
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
  USE mo_kind,               ONLY: wp, i8
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_impl_constants,     ONLY: max_dom, vname_len,                       &
    &                              max_var_ml, MODE_IFSANA, MODE_DWDANA,     &
    &                              MODE_IAU, MODE_IAU_OLD, MODE_COMBINED,    &
    &                              MODE_COSMO, MODE_ICONVREMAP, ivexpol
  USE mo_io_units,           ONLY: nnml, nnml_output, filename_max
  USE mo_namelist,           ONLY: position_nml, positioned, open_nml, close_nml
  USE mo_mpi,                ONLY: my_process_is_stdio 
  USE mo_initicon_config,    ONLY: initicon_config,      &
    & config_init_mode           => init_mode,           &
    & config_nlevsoil_in         => nlevsoil_in,         &
    & config_zpbl1               => zpbl1,               &
    & config_zpbl2               => zpbl2,               &
    & config_lread_ana           => lread_ana,           &
    & config_qcana_mode          => qcana_mode,          &
    & config_qiana_mode          => qiana_mode,          &
    & config_qrsgana_mode        => qrsgana_mode,        &
    & config_qnxana_2mom_mode    => qnxana_2mom_mode,    &
    & config_use_lakeiceana      => use_lakeiceana,      &
    & config_lconsistency_checks => lconsistency_checks, &
    & config_ifs2icon_filename   => ifs2icon_filename,   &
    & config_dwdfg_filename      => dwdfg_filename,      &
    & config_dwdana_filename     => dwdana_filename,     &
    & config_l_coarse2fine_mode  => l_coarse2fine_mode,  &
    & config_lp2cintp_incr       => lp2cintp_incr,       &
    & config_lp2cintp_sfcana     => lp2cintp_sfcana,     &
    & config_lvert_remap_fg      => lvert_remap_fg,      &
    & config_ltile_coldstart     => ltile_coldstart,     &
    & config_ltile_init          => ltile_init,          &
    & config_icpl_da_sfcevap     => icpl_da_sfcevap,     &
    & config_smi_relax_timescale => smi_relax_timescale, &
    & config_icpl_da_skinc       => icpl_da_skinc,       &
    & config_icpl_da_snowalb     => icpl_da_snowalb,     &
    & config_icpl_da_sfcfric     => icpl_da_sfcfric,     &
    & config_scalfac_da_sfcfric  => scalfac_da_sfcfric,  &
    & config_icpl_da_tkhmin      => icpl_da_tkhmin,      &
    & config_icpl_da_seaice      => icpl_da_seaice,      &
    & config_dt_ana              => dt_ana,              &
    & config_adjust_tso_tsnow    => adjust_tso_tsnow,    &
    & config_filetype            => filetype,            &
    & config_dt_iau              => dt_iau,              &
    & config_iterate_iau         => iterate_iau,         &
    & config_timeshift           => timeshift,           &
    & config_type_iau_wgt        => type_iau_wgt,        &
    & config_niter_divdamp       => niter_divdamp,       &
    & config_niter_diffu         => niter_diffu,         &
    & config_itype_vert_expol    => itype_vert_expol,    &
    & config_ana_varnames_map_file => ana_varnames_map_file, &
    & config_pinit_seed          => pinit_seed,          &
    & config_pinit_amplitude     => pinit_amplitude,     &
    & config_lcouple_ocean_coldstart => lcouple_ocean_coldstart, &
    & config_fire2d_filename     => fire2d_filename

  USE mo_nml_annotate,       ONLY: temp_defaults, temp_settings


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_initicon_namelist

CONTAINS

!-------------------------------------------------------------------------
!
!
 !!  Initialization of the initicon coordinate namelist

 SUBROUTINE read_initicon_namelist( filename )
    
  CHARACTER(LEN=*), INTENT(IN) :: filename

  !local variable
  INTEGER :: i_status
  INTEGER :: z_go_init(7)   ! for consistency check
  INTEGER :: iunit
  INTEGER :: jg

  CHARACTER(len=*), PARAMETER ::  &
    &  routine = 'mo_initicon_nml: read_initicon_namelist'

    ! Variables of that type contain a list of all mandatory input fields
  ! This list can include a subset or the entire set of mandatory fields.
  TYPE t_check_input
    CHARACTER(LEN=vname_len) :: list(max_var_ml)
  END TYPE t_check_input


  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the init_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  INTEGER  :: init_mode     ! initialization mode
  INTEGER  :: nlevsoil_in   ! number of soil levels of input data

  REAL(wp) :: zpbl1, zpbl2  ! AGL heights used for vertical gradient computation
  LOGICAL  :: lread_ana     ! If .TRUE., read analysis fields are read from analysis file
                            ! dwdana_filename. If .FALSE., ICON is soleyly started 
                            ! from first guess fields.   
  LOGICAL  :: lconsistency_checks    ! check validity of input fields (FG and ANA)

  LOGICAL  :: l_coarse2fine_mode(max_dom)  ! If true, apply special corrections for interpolation from coarse
                                           ! to fine resolutions over mountainous terrain
  LOGICAL  :: lp2cintp_incr(max_dom) ! If true, perform parent-to-child interpolation of atmospheric data
                                     ! assimilation increments
  LOGICAL  :: lp2cintp_sfcana(max_dom) ! If true, perform parent-to-child interpolation of
                                       ! surface analysis data
  LOGICAL  :: lcouple_ocean_coldstart  ! If true, initialize newly defined land points from ICON-O 
                                       ! with default T and Q profiles 

  LOGICAL  :: ltile_coldstart  ! If true, initialize tile-based surface fields from first guess with tile-averaged fields

  LOGICAL  :: ltile_init       ! If true, initialize tile-based surface fields from first guess without tiles

  INTEGER  :: icpl_da_sfcevap  ! Type of coupling between data assimilation and model parameters 
                               ! affecting surface evaporation (plants + bare soil)

  REAL(wp) :: smi_relax_timescale ! Time scale (days) for ICON-internal soil moisture relaxation

  INTEGER  :: icpl_da_skinc    ! Coupling between data assimilation and skin conductivity

  INTEGER  :: icpl_da_snowalb  ! Coupling between data assimilation and snow albedo

  INTEGER  :: icpl_da_sfcfric  ! Coupling between data assimilation and surface friction (roughness length and SSO blocking)

  REAL(wp) :: scalfac_da_sfcfric ! Scaling factor for adaptive surface friction

  INTEGER  :: icpl_da_tkhmin   ! Coupling between data assimilation and near-surface profiles of minimum vertical diffusion

  INTEGER  :: icpl_da_seaice   ! Coupling between data assimilation and sea ice

  REAL(wp) :: dt_ana           ! Time interval of assimilation cycle [s] (relevant for icpl_da_sfcevap >= 2)

  LOGICAL  :: adjust_tso_tsnow ! Apply T increments for lowest model level also to snow and upper soil layers

  LOGICAL  :: use_lakeiceana   ! If true, use ice fraction analysis data also over lakes (otherwise sea points only)

  LOGICAL  :: lvert_remap_fg   ! If true, vertical remappting of first guess input is performed

  INTEGER  :: qcana_mode, qiana_mode, qrsgana_mode, qnxana_2mom_mode  ! mode of processing QC/QI/QR/QS/QG/QH/QNX increments

  INTEGER  :: filetype      ! One of CDI's FILETYPE\_XXX constants. Possible values: 2 (=FILETYPE\_GRB2), 4 (=FILETYPE\_NC2)

  REAL(wp) :: dt_iau        ! Time interval during which incremental analysis update (IAU) is performed [s]. 
                            ! Only required for init_mode=MODE_IAU, MODE_IAU_OLD

  !> Allows IAU runs to start earlier than the nominal simulation start
  ! date without showing up in the output metadata:
  REAL(wp) :: dt_shift

  INTEGER  :: type_iau_wgt  ! Type of weighting function for IAU.
                            ! 1: Top-hat
                            ! 2: SIN2
                            ! Only required for init_mode=MODE_IAU, MODE_IAU_OLD
  LOGICAL  :: iterate_iau   ! if .TRUE., iterate IAU phase with halved dt_iau in first iteration

  INTEGER  :: niter_divdamp ! number of divergence damping iterations on wind increment from DA
  INTEGER  :: niter_diffu   ! number of diffusion iterations on wind increment from DA

  INTEGER :: itype_vert_expol ! Type of vertical extrapolation of initial data. 
                              ! 1: Linear extrapolation (standard setting) 
                              ! 2: Blending with climatology 
                              ! (intended for simulations with the upper-atmosphere configuration)


  TYPE(t_check_input) :: check_ana(max_dom)  ! patch-specific list of mandatory analysis fields.
                                             ! This list can include a subset or the 
                                             ! entire set of default analysis fields.

  TYPE(t_check_input) :: check_fg(max_dom)  ! A small subset of first guess input fields are 
                                            ! declared 'optional' in ICON. By adding them to this list, 
                                            ! they will become mandatory, meaning that the model 
                                            ! aborts if any of these fields is missing. This list can 
                                            ! include a subset of the optional first guess fields, 
                                            ! or even the entire set of first guess fields (however, 
                                            ! it only effects the optional ones). On default this list 
                                            ! is empty, meaning that optional first guess fields 
                                            ! experience a cold-start initialization if they are missing. 
                                            ! The model does not abort.

   ! patch-specific list of mandatory first guess fields.
                                             ! This list can include a subset or the 
                                             ! entire set of default first guess fields.

  ! IFS2ICON input filename, may contain keywords, by default
  ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: ifs2icon_filename

  ! DWD-FG input filename, may contain keywords, by default
  ! dwdfg_filename = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdfg_filename

  ! DWD-ANA input filename, may contain keywords, by default
  ! dwdana_filename = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  CHARACTER(LEN=filename_max) :: dwdana_filename

  ! analysis file: dictionary which maps internal variable names onto
  ! GRIB2 shortnames or NetCDF var names.
  CHARACTER(LEN=filename_max) :: ana_varnames_map_file

  ! perturb initial conditions. perturbation is only applied for pinit_seed /= 0
  INTEGER(i8) :: pinit_seed = 0_i8
  REAL(wp) :: pinit_amplitude = 0._wp

  CHARACTER(LEN=filename_max) :: & !< Filename that contains wildfire precursor emissions (2d-aerosol, iprog_aero=3)
    &  fire2d_filename             !< Allowed keywords: <species>, <gridfile>, <nroot>, <nroot0>, <jlev>, <idom>, <yyyymmdd>

  NAMELIST /initicon_nml/ init_mode, zpbl1, zpbl2, l_coarse2fine_mode,      &
                          nlevsoil_in, lread_ana,                           &
                          lconsistency_checks,                              &
                          ifs2icon_filename, dwdfg_filename,                &
                          dwdana_filename, filetype, dt_iau, dt_shift,      &
                          type_iau_wgt, check_ana, check_fg,                &
                          ana_varnames_map_file, lp2cintp_incr,             &
                          lp2cintp_sfcana, use_lakeiceana,                  &
                          ltile_coldstart, ltile_init,                      &
                          lvert_remap_fg, iterate_iau, niter_divdamp,       &
                          niter_diffu, qcana_mode, qiana_mode, qrsgana_mode,&
                          qnxana_2mom_mode, itype_vert_expol, pinit_seed,   &
                          pinit_amplitude, icpl_da_sfcevap, dt_ana,         &
                          icpl_da_skinc, icpl_da_snowalb, adjust_tso_tsnow, &
                          icpl_da_sfcfric, lcouple_ocean_coldstart,         &
                          icpl_da_tkhmin, icpl_da_seaice, fire2d_filename,  &
                          scalfac_da_sfcfric, smi_relax_timescale

  !------------------------------------------------------------
  ! 2.0 set up the default values for initicon
  !------------------------------------------------------------
  !
  !
  init_mode   = MODE_IFSANA    ! Start from IFS analysis
  nlevsoil_in = 4              ! number of soil levels of input data
  zpbl1       = 500._wp        ! AGL heights used for computing vertical 
  zpbl2       = 1000._wp       ! gradients
  lread_ana   = .TRUE.         ! true: read analysis fields from file dwdana_filename
                               ! false: start ICON from first guess file (no analysis)
  qcana_mode  = 0              ! 1: add QC increments on QV increments (0: ignore them)
                               ! 2: add QC increments on QV increments in case of subsaturation and to QC otherwise
  qiana_mode  = 0              ! 0/1: use/ignore QI increments
  qrsgana_mode= 0              ! 0/1: use/ignore QR/QS/QG/QH increments
  qnxana_2mom_mode= 0          ! 0/1: use/ignore QNX increments
  use_lakeiceana = .FALSE.     ! do not use ice fraction analysis data over freshwater lakes
  lconsistency_checks = .TRUE. ! check validity of input fields  
  filetype    = -1             ! "-1": undefined
  dt_iau      = 10800._wp      ! 3-hour interval for IAU
  iterate_iau = .FALSE.        ! no iteration of IAU
  dt_shift    = 0._wp          ! do not shift actual simulation start backward
  niter_diffu = 10             ! number of diffusion iterations on wind increment from DA
  niter_divdamp = 25           ! number of divergence damping iterations on wind increment from DA
  type_iau_wgt= 1              ! Top-hat weighting function
  itype_vert_expol = ivexpol%lin ! linear vertical extrapolation of initial data

  DO jg=1,SIZE(check_ana)
    check_ana(jg)%list(:) = '' ! list of mandatory analysis fields. This list can include a subset 
  ENDDO                        ! or the entire set of default analysis fields. If any of these fields
                               ! is missing in the analysis file, the model aborts. On default 
                               ! this list is empty, meaning that fields which are missing in the 
                               ! analysis file (when compared to the default set), are simply 
                               ! taken from the first guess.

  DO jg=1,SIZE(check_fg)
    check_fg(jg)%list(:) = '' ! A small subset of first guess input fields are declared 'optional' in ICON. 
  ENDDO                       ! By adding them to this list, they will become mandatory, meaning that the model 
                              ! aborts if any of these fields is missing. This list can include a subset 
                              ! of the optional first guess fields, or even the entire set of first guess fields 
                              ! (however, it only effects the optional first guess fields). On default 
                              ! this list is empty, meaning that optional first guess fields experience 
                              ! a cold-start initialization if they are missing. The model does not abort.

  ana_varnames_map_file = " "

  ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdfg_filename    = "<path>dwdFG_R<nroot>B<jlev>_DOM<idom>.nc"
  dwdana_filename   = "<path>dwdana_R<nroot>B<jlev>_DOM<idom>.nc"
  l_coarse2fine_mode(:) = .FALSE. ! true: apply corrections for coarse-to-fine-mesh interpolation
  lp2cintp_incr(:)      = .FALSE. ! true: perform parent-to-child interpolation of atmospheric data assimilation increments
  lp2cintp_sfcana(:)    = .FALSE. ! true: perform parent-to-child interpolation of surface analysis data
  lcouple_ocean_coldstart= .TRUE. ! true: initialize newly defined land points from ICON-O with default T and Q profiles
  ltile_coldstart       = .FALSE. ! true: initialize tile-based surface fields from first guess with tile-averaged fields
  ltile_init            = .FALSE. ! true: initialize tile-based surface fields from first guess without tiles
  lvert_remap_fg        = .FALSE. ! true: perform vertical remapping of first-guess input

  icpl_da_sfcevap = 0   ! Coupling between data assimilation and parameters affecting surface evaporation
                        ! 0: none
                        ! 1: use filtered T2M bias 
                        ! 2: use filtered T2M bias and filtered RH increment at lowest model level
                        ! 3: use filtered T and RH increments at lowest model level
                        ! 4: as 3, but uses cr_bsmin instead of c_soil for adapting bare-soil evaporation

  smi_relax_timescale = 20._wp ! Time scale (days) for ICON-internal soil moisture relaxation

  icpl_da_skinc = 0     ! Coupling between data assimilation and skin conductivity
                        ! 0: off, 1: on, 2: as 1, plus soil heat conductivity and capacity

  icpl_da_snowalb = 0   ! Coupling between data assimilation and snow albedo
                        ! 0: off, 1: on, 2: as 1, plus sea-ice albedo

  icpl_da_sfcfric = 0   ! Coupling between data assimilation and surface friction (roughness length and SSO blocking)
                        ! 0: off, 1:on

  scalfac_da_sfcfric = 2.5_wp  ! scaling factor for adaptive surface friction

  icpl_da_tkhmin   = 0  ! Coupling between data assimilation and near-surface profile of minimum vertical diffusion for heat
                        ! 0: off, 1:on

  icpl_da_seaice   = 0  ! Coupling between data assimilation and sea ice
                        ! 0: off, 1:adjustment of t_seaice to filtered DA increment, 2: 1+ adaptive tuning of bottom heat flux

  adjust_tso_tsnow = .FALSE. ! If .TRUE., apply T increments for lowest model level also to snow and upper soil layers

  dt_ana  = 10800._wp   ! Time interval of assimilation cycle (relevant for icpl_da_sfcevap >= 2; set 3600 s for ICON-D2

  pinit_seed        = 0_i8        ! =0: do not perturb initial data. >0: perturb initial data with this as seed
  pinit_amplitude   = 0._wp       ! amplitude of the initial perturbation for numerical tolerance test

  fire2d_filename = 'gfas2d_emi_<species>_<gridfile>_<yyyymmdd>.nc'

  !------------------------------------------------------------
  ! 3.0 Read the initicon namelist.
  !------------------------------------------------------------
  ! (done so far by all MPI processes)
  !
  CALL open_nml(TRIM(filename))
  CALL position_nml ('initicon_nml', status=i_status)
  IF (my_process_is_stdio()) THEN
    iunit = temp_defaults()
    WRITE(iunit, initicon_nml)  ! write defaults to temporary text file
  END IF
  SELECT CASE (i_status)
  CASE (positioned)
    READ (nnml, initicon_nml)                                      ! overwrite default settings
    IF (my_process_is_stdio()) THEN
      iunit = temp_settings()
      WRITE(iunit, initicon_nml)  ! write settings to temporary text file
    END IF
  END SELECT
  CALL close_nml


  !------------------------------------------------------------
  ! 4.0 check the consistency of the parameters
  !------------------------------------------------------------
  !
  z_go_init = (/MODE_IFSANA,MODE_DWDANA,MODE_IAU,MODE_IAU_OLD,MODE_COMBINED,MODE_COSMO,MODE_ICONVREMAP/)
  IF (ALL(z_go_init /= init_mode)) THEN
    CALL finish( TRIM(routine),                         &
      &  'Invalid initialization mode. init_mode must be between 1 and 7')
  ENDIF

  ! Check whether init_mode and lread_ana are consistent
  IF (ANY((/MODE_COMBINED,MODE_COSMO,MODE_ICONVREMAP/)==init_mode) .AND. lread_ana) THEN
    lread_ana = .FALSE.
    WRITE(message_text,'(a,i2,a)') 'init_mode=', init_mode, &
      '. no analysis required => lread_ana re-set to .FALSE.'
    CALL message(TRIM(routine),message_text)
  ENDIF

  ! Check whether an analysis file is provided, if lread_ana=.TRUE.
  IF (lread_ana) THEN
    IF (dwdana_filename ==' ') THEN
    CALL finish( TRIM(routine),                         &
      &  'dwdana_filename required, but missing.')
    ENDIF
  ENDIF


  ! Setting the first entry of lp2cintp_incr / lp2cintp_sfcana to true activates parent-to-child interpolation
  ! of DA increments / surface analysis for all domains
  IF (lp2cintp_incr(1)) THEN
    lp2cintp_incr(2:max_dom) = .TRUE.
  ENDIF
  IF (lp2cintp_sfcana(1)) THEN
    lp2cintp_sfcana(2:max_dom) = .TRUE.
  ENDIF
  ! To simplify runtime flow control, set the switches for the global domain to false
  lp2cintp_incr(1)   = .FALSE.
  lp2cintp_sfcana(1) = .FALSE.

  ! make sure that dt_shift is negative or 0.
  IF ( dt_shift > 0._wp ) THEN
    WRITE(message_text,'(a,f8.2,a)') 'dt_shift=', dt_shift, &
      ' not allowed. Must be NEGATIVE or 0.'
    CALL finish(TRIM(routine),message_text)
  ENDIF

  ! this is needed because the I/O of the filtered T increment is controlled via icpl_da_sfcevap >= 3
  IF (icpl_da_snowalb >= 1 .AND. icpl_da_sfcevap < 3) THEN
    WRITE(message_text,'(a)') 'icpl_da_snowalb >= 1 must be combined with icpl_da_sfcevap >= 3'
    CALL finish(TRIM(routine),message_text)
  ENDIF
  IF (icpl_da_seaice >= 1 .AND. icpl_da_sfcevap < 3) THEN
    WRITE(message_text,'(a)') 'icpl_da_seaice >= 1 must be combined with icpl_da_sfcevap >= 3'
    CALL finish(TRIM(routine),message_text)
  ENDIF

  IF (icpl_da_tkhmin >= 1 .AND. (icpl_da_skinc == 0 .OR. icpl_da_sfcevap <= 2) ) THEN
    WRITE(message_text,'(a)') 'icpl_da_tkhmin = 1 must be combined with icpl_da_sfcevap > 2 and icpl_da_skinc > 0'
    CALL finish(TRIM(routine),message_text)
  ENDIF

  ! IAU iteration is meaningless if the model starts without backward time shift
  IF (dt_shift == 0._wp) THEN
    iterate_iau = .FALSE.
  END IF

  ! Check setting for vertical extrapolation
  SELECT CASE(itype_vert_expol)
  CASE(ivexpol%lin, ivexpol%upatmo) 
    ! Ok
  CASE DEFAULT
    CALL finish( TRIM(routine),'Invalid value for itype_vert_expol.' )
  END SELECT


  !------------------------------------------------------------
  ! 5.0 Fill the configuration state
  !------------------------------------------------------------

  config_init_mode           = init_mode
  config_nlevsoil_in         = nlevsoil_in
  config_zpbl1               = zpbl1
  config_zpbl2               = zpbl2
  config_lread_ana           = lread_ana
  config_qcana_mode          = qcana_mode
  config_qiana_mode          = qiana_mode
  config_qrsgana_mode        = qrsgana_mode
  config_qnxana_2mom_mode    = qnxana_2mom_mode
  config_use_lakeiceana      = use_lakeiceana
  config_lconsistency_checks = lconsistency_checks
  config_ifs2icon_filename   = ifs2icon_filename
  config_dwdfg_filename      = dwdfg_filename
  config_dwdana_filename     = dwdana_filename
  config_l_coarse2fine_mode  = l_coarse2fine_mode
  config_lp2cintp_incr       = lp2cintp_incr
  config_lp2cintp_sfcana     = lp2cintp_sfcana
  config_ltile_coldstart     = ltile_coldstart
  config_ltile_init          = ltile_init
  config_icpl_da_sfcevap     = icpl_da_sfcevap
  config_smi_relax_timescale = smi_relax_timescale
  config_icpl_da_skinc       = icpl_da_skinc
  config_icpl_da_snowalb     = icpl_da_snowalb
  config_icpl_da_sfcfric     = icpl_da_sfcfric
  config_scalfac_da_sfcfric  = scalfac_da_sfcfric
  config_icpl_da_tkhmin      = icpl_da_tkhmin
  config_icpl_da_seaice      = icpl_da_seaice
  config_dt_ana              = dt_ana
  config_adjust_tso_tsnow    = adjust_tso_tsnow
  config_lvert_remap_fg      = lvert_remap_fg
  config_filetype            = filetype
  config_dt_iau              = dt_iau
  config_iterate_iau         = iterate_iau
  config_timeshift%dt_shift  = dt_shift
  config_type_iau_wgt        = type_iau_wgt
  config_ana_varnames_map_file = ana_varnames_map_file
  config_niter_divdamp         = niter_divdamp
  config_niter_diffu           = niter_diffu
  config_itype_vert_expol      = itype_vert_expol
  config_pinit_seed            = pinit_seed
  config_pinit_amplitude       = pinit_amplitude
  config_lcouple_ocean_coldstart = lcouple_ocean_coldstart
  config_fire2d_filename       = TRIM(fire2d_filename)

  DO jg=1,max_dom
    initicon_config(jg)%ana_checklist = check_ana(jg)%list
    initicon_config(jg)%fg_checklist  = check_fg(jg)%list
  ENDDO

  ! write the contents of the namelist to an ASCII file

  IF(my_process_is_stdio()) WRITE(nnml_output,nml=initicon_nml)

END SUBROUTINE read_initicon_namelist

END MODULE mo_initicon_nml
