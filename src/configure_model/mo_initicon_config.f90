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

MODULE mo_initicon_config

  USE mo_kind,               ONLY: wp, i8
  USE mo_exception,          ONLY: finish
  USE mo_util_string,        ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                              int2string
  USE mo_impl_constants,     ONLY: max_dom, vname_len, max_var_ml, MAX_CHAR_LENGTH,  &
    &                              MODE_IFSANA, MODE_COMBINED, MODE_COSMO,           &
    &                              MODE_IAU, MODE_IAU_OLD, MODE_ICONVREMAP, nclass_aero, &
    &                              ivexpol
  USE mo_io_units,           ONLY: filename_max
  USE mo_io_util,            ONLY: get_filetype
  USE mo_model_domain,       ONLY: t_patch
  USE mo_grid_config,        ONLY: l_limited_area, nroot
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mtime,                 ONLY: timedelta
  USE mo_upatmo_config,      ONLY: upatmo_config

  IMPLICIT NONE

  PRIVATE


  ! Types
  PUBLIC :: t_initicon_config

  ! Variables
  PUBLIC :: init_mode, nlevsoil_in, zpbl1, zpbl2
  PUBLIC :: dt_iau
  PUBLIC :: type_iau_wgt
  PUBLIC :: iterate_iau
  PUBLIC :: l_sst_in
  PUBLIC :: use_lakeiceana
  PUBLIC :: lread_ana
  PUBLIC :: lread_vn
  PUBLIC :: lread_tke
  PUBLIC :: lconsistency_checks
  PUBLIC :: l_coarse2fine_mode
  PUBLIC :: lp2cintp_incr, lp2cintp_sfcana
  PUBLIC :: qcana_mode, qiana_mode, qrsgana_mode, qnxana_2mom_mode
  PUBLIC :: lcouple_ocean_coldstart
  PUBLIC :: ltile_coldstart
  PUBLIC :: ltile_init
  PUBLIC :: icpl_da_sfcevap, icpl_da_skinc, icpl_da_snowalb, icpl_da_sfcfric, icpl_da_tkhmin, dt_ana
  PUBLIC :: adjust_tso_tsnow, icpl_da_seaice
  PUBLIC :: lvert_remap_fg
  PUBLIC :: ifs2icon_filename
  PUBLIC :: dwdfg_filename
  PUBLIC :: dwdana_filename
  PUBLIC :: filetype
  PUBLIC :: ana_varnames_map_file
  PUBLIC :: init_mode_soil
  PUBLIC :: is_iau_active
  PUBLIC :: iau_wgt_dyn, iau_wgt_adv
  PUBLIC :: niter_divdamp, niter_diffu
  PUBLIC :: t_timeshift
  PUBLIC :: timeshift
  PUBLIC :: initicon_config
  PUBLIC :: aerosol_fg_present
  PUBLIC :: lanaread_tseasfc
  PUBLIC :: itype_vert_expol
  PUBLIC :: pinit_seed
  PUBLIC :: pinit_amplitude
  PUBLIC :: fire2d_filename

  ! Subroutines
  PUBLIC :: configure_initicon

  ! Functions
  PUBLIC :: generate_filename
  PUBLIC :: fgFilename
  PUBLIC :: fgFiletype
  PUBLIC :: anaFilename
  PUBLIC :: anaFiletype

  TYPE t_timeshift
    REAL(wp)                 :: dt_shift
    TYPE(timedelta), POINTER :: mtime_shift
    TYPE(timedelta), POINTER :: mtime_absshift   ! absolute value
  END TYPE t_timeshift

  ! ----------------------------------------------------------------------------
  ! 1.0 Namelist variables for the init_icon preprocessing program
  ! ----------------------------------------------------------------------------
  !
  TYPE :: t_initicon_config
    CHARACTER(LEN=vname_len) :: ana_checklist(max_var_ml) ! list of mandatory analysis fields. 
                                                        ! This list can include a subset or the 
                                                        ! entire set of default analysis fields.
    CHARACTER(LEN=vname_len) :: fg_checklist(max_var_ml) ! list of mandatory first guess fields. 
                                                        ! This list can include a subset or the 
                                                        ! entire set of default first guess fields.
  END TYPE t_initicon_config
  !
  ! probably those which are domain-dependent could be included into aboves type lateron
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

  INTEGER  :: icpl_da_sfcevap  ! Type of coupling between data assimilation and model parameters affecting surface evaporation (plants + bare soil)

  INTEGER  :: icpl_da_skinc    ! Coupling between data assimilation and skin conductivity

  INTEGER  :: icpl_da_snowalb  ! Coupling between data assimilation and snow albedo

  INTEGER  :: icpl_da_sfcfric  ! Coupling between data assimilation and surface friction (roughness length and SSO blocking)

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

  TYPE(t_timeshift) :: &    ! Allows IAU runs to start earlier than the nominal simulation start date 
    &  timeshift            ! without showing up in the output metadata

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

   
  ! ----------------------------------------------------------------------------
  ! Derived variables / variables based on input file contents
  ! ----------------------------------------------------------------------------

  !> control variable that specifies if u/v or vn are read as wind
  !  field input
  LOGICAL :: lread_vn  = .FALSE.

  !> control variable that specifies if TKE has been found in the
  !  input (used for MODE_ICONVREMAP only)
  LOGICAL :: lread_tke = .FALSE.

  !> logical switch, if sea surface temperature is provided as input
  LOGICAL :: l_sst_in  = .TRUE.

  !> initialization mode of soil model (coldstart, warmstart,
  !  warmstart+IAU)
  INTEGER :: init_mode_soil

  !> determines whether IAU is active at current time
  LOGICAL :: is_iau_active = .FALSE.

  REAL(wp):: iau_wgt_dyn = 0._wp    !< IAU weight for dynamics fields 
  REAL(wp):: iau_wgt_adv = 0._wp    !< IAU weight for tracer fields

  !> registers if aerosol fields have been read from the first-guess
  !  data
  LOGICAL :: aerosol_fg_present(max_dom,nclass_aero) = .FALSE.

  !> registers if SST and sea ice fraction data have been read from
  !  analysis
  LOGICAL :: lanaread_tseasfc(max_dom) = .FALSE.

  TYPE(t_initicon_config), TARGET :: initicon_config(0:max_dom)

  ! perturb initial conditions. perturbation is only applied for pinit_seed /= 0
  INTEGER(i8) :: pinit_seed = 0_i8
  REAL(wp) :: pinit_amplitude = 0._wp

  CHARACTER(LEN=filename_max) :: & !< Filename that contains wildfire precursor emissions (2d-aerosol, iprog_aero=3)
    &  fire2d_filename             !< Allowed keywords: <species>, <gridfile>, <nroot>, <nroot0>, <jlev>, <idom>, <yyyymmdd>

CONTAINS

  !! setup additional initicon control variables
  !!
  !! Setup of additional initicon control variables depending on the 
  !! initicon-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  SUBROUTINE configure_initicon()
    !
    CHARACTER(len=*), PARAMETER :: routine = 'mo_initicon_config:configure_initicon'

    !-----------------------------------------------------------------------
    !


    IF ( ANY((/MODE_IFSANA,MODE_COMBINED,MODE_COSMO/) == init_mode) ) THEN
       init_mode_soil = 1   ! full coldstart is executed
       ! i.e. w_so_ice and h_snow are re-diagnosed
    ELSE IF (l_limited_area .AND. init_mode == MODE_ICONVREMAP .AND. .NOT. lread_ana) THEN
       init_mode_soil = 1   ! same initialization for limited-area cold start
    ELSE IF ( ANY((/MODE_IAU, MODE_IAU_OLD/) == init_mode) ) THEN
       init_mode_soil = 3  ! warmstart (within assimilation cycle) with analysis increments for h_snow
    ELSE
       init_mode_soil = 2  ! warmstart with full fields for h_snow from snow analysis
    ENDIF

    !
    ! set switch(es) for vertical extrapolation of initial data
    !
    ! just to make sure
    IF (.NOT. ALLOCATED(upatmo_config)) THEN 
      CALL finish('mo_initicon_config:configure_initicon', &
        &         'upatmo_config is not allocated')
    ENDIF
    SELECT CASE(itype_vert_expol)
    CASE(ivexpol%lin)
      ! linear extrapolation is the standard case 
      upatmo_config(:)%exp%l_expol = .FALSE.
    CASE(ivexpol%upatmo) 
      ! this case is intended for (but not necessarily limited to) 
      ! upper-atmosphere simulations, further specifiers for 
      ! this extrapolation can be set in 'upatmo_nml'
      upatmo_config(:)%exp%l_expol = .TRUE.
    END SELECT
    upatmo_config(:)%exp%l_initicon_config  = .TRUE.

  END SUBROUTINE configure_initicon



  FUNCTION generate_filename(input_filename, model_base_dir, &
    &                        nroot, jlev, idom)  RESULT(result_str)
    CHARACTER(len=*), INTENT(IN)   :: input_filename, &
      &                               model_base_dir
    INTEGER,          INTENT(IN)   :: nroot, jlev, idom
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str
    TYPE (t_keyword_list), POINTER :: keywords => NULL()

    CALL associate_keyword("<path>",   TRIM(model_base_dir),             keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i0)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
    ! replace keywords in "input_filename", which is by default
    ! ifs2icon_filename = "<path>ifs2icon_R<nroot>B<jlev>_DOM<idom>.nc"
    result_str = TRIM(with_keywords(keywords, TRIM(input_filename)))

  END FUNCTION generate_filename


  FUNCTION fgFilename(p_patch) RESULT(resultVar)
    CHARACTER(LEN = filename_max) :: resultVar
    TYPE(t_patch), INTENT(IN) :: p_patch

    resultVar = generate_filename(dwdfg_filename, getModelBaseDir(), nroot, p_patch%level, p_patch%id)
  END FUNCTION fgFilename

  FUNCTION anaFilename(p_patch) RESULT(resultVar)
    CHARACTER(LEN = filename_max) :: resultVar
    TYPE(t_patch), INTENT(IN) :: p_patch

    resultVar = generate_filename(dwdana_filename, getModelBaseDir(), nroot, p_patch%level, p_patch%id)
  END FUNCTION anaFilename

  INTEGER FUNCTION fgFiletype() RESULT(resultVar)
    IF(filetype == -1) THEN
        ! get_filetype() ONLY uses the suffix, which IS already a part of the template IN dwdfg_filename.
        ! This IS why it suffices to USE the dwdfg_filename directly here without expanding it first via generate_filename().
        resultVar = get_filetype(TRIM(dwdfg_filename))
    ELSE
        resultVar = filetype
    END IF
  END FUNCTION fgFiletype

  INTEGER FUNCTION anaFiletype() RESULT(resultVar)
    IF(filetype == -1) THEN
        ! get_filetype() ONLY uses the suffix, which IS already a part of the template IN dwdana_filename.
        ! This IS why it suffices to USE the dwdana_filename directly here without expanding it first via generate_filename().
        resultVar = get_filetype(TRIM(dwdana_filename))
    ELSE
        resultVar = filetype
    END IF
  END FUNCTION anaFiletype

END MODULE mo_initicon_config
