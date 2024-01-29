!
! mo_art_diag_state
! This module creates the diagnostic state for ART variables
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

MODULE mo_art_diag_state
! ICON
  USE mo_kind,                          ONLY: wp, i8
  USE mo_model_domain,                  ONLY: p_patch
  USE mo_parallel_config,               ONLY: nproma
  USE mo_run_config,                    ONLY: ntracer
  USE mo_exception,                     ONLY: finish, message, message_text
  USE mo_cdi_constants,                 ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_zaxis_type,                    ONLY: ZA_HYBRID, ZA_SURFACE, ZA_ATMOSPHERE,    &
    &                                         ZA_PRES_FL_BOT_TOP
  USE mo_cdi,                           ONLY: DATATYPE_PACK16,                         &
    &                                         DATATYPE_PACK24, CDI_DATATYPE_PACK32,    &
    &                                         TSTEP_ACCUM, TSTEP_AVG, TSTEP_MAX,       &
    &                                         GRID_UNSTRUCTURED
  USE mo_cf_convention,                 ONLY: t_cf_var
  USE mo_grib2,                         ONLY: t_grib2_var, grib2_var, t_grib2_int_key, &
    &                                         OPERATOR(+)
  USE mo_var_list,                      ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: CLASS_DEFAULT, CLASS_CHEM,               &
    &                                         CLASS_CHEM_STAT, CLASS_CHEM_OPTP,        &
    &                                         CLASS_DISTR, CLASS_DISTR_STAT
  USE mo_var_metadata,                  ONLY: create_vert_interp_metadata,             &
    &                                         create_hor_interp_metadata,              &
    &                                         vintp_types, get_var_name
  USE mo_impl_constants,                ONLY: VINTP_METHOD_LIN,                        &
    &                                         HINTP_TYPE_LONLAT_NNB,                   &
    &                                         HINTP_TYPE_LONLAT_BCTR,                  &
    &                                         SUCCESS, MAX_CHAR_LENGTH
  USE mo_var_groups,                    ONLY: groups, MAX_GROUPS
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_util_string,                   ONLY: tolower, str_replace, split_string
  USE mo_action,                        ONLY: ACTION_RESET, new_action, actions
  USE mtime,                            ONLY: max_timedelta_str_len, getPTStringFromMS

  USE mo_art_general_interface,         ONLY: getNetcdfPrecision
  ! ART
  USE mo_art_config,                    ONLY: art_config, npreslay
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_data,                 ONLY: t_art_chem, nphot
  USE mo_art_diag_types,                ONLY: t_art_diag, art_diag_tracer_index
  USE mo_art_modes_linked_list,         ONLY: p_mode_state, t_mode
  USE mo_art_modes,                     ONLY: t_fields_radio, t_fields_2mom
  USE mo_art_read_xml,                  ONLY: t_xml_file, art_open_xml_file, art_close_xml_file,  &
                                          &   art_get_childnumber_xml, art_read_elements_xml
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN, IART_MAX_RADIOACT,                 &
                                          &   IART_MAX_DIAGCONT, IART_ACC_DRYDEPO, IART_ACC_SEDIM,&
                                          &   IART_ACC_EMISS, IART_EMISS, IART_ACC_WETDEPO_GSCP,  &
                                          &   IART_ACC_WETDEPO_CON, IART_ACC_WETDEPO_RRSFC
  USE mo_art_io_constants,              ONLY: IART_FILENAMELEN
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
#ifdef __ART_GPL
  USE messy_mecca_kpp_global,           ONLY: NREACT, ihs_MAX
  USE mo_art_psc_types,                 ONLY: t_art_psc
#endif
  USE mo_art_external_types,            ONLY: t_art_volc_fplume
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string
  USE mo_art_external_types,            ONLY: t_art_pollen_table, t_art_external
  USE mo_art_external_init_pollen,      ONLY: art_extinit_pollen_fordiag

#include "add_var_acc_macro.inc"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_diag_state'

  INTEGER                               :: &
    &  ndiag_xml                             !< number of diagn. vars in XML
  TYPE(t_key_value_store), ALLOCATABLE  :: &
    &  meta_storage(:)                       !< Key-Value storage to be filled with data from 
                                             !  child elements

  PUBLIC :: art_create_diagnostics

CONTAINS

SUBROUTINE art_create_diagnostics(jg, p_diag_list, p_prog_list)
!<
! SUBROUTINE art_create_diagnostics
! This subroutine creates the diagnostic state for ART variables
! Part of Module: mo_art_diag_state
! Author: Daniel Rieger, KIT
! Initial Release: 2014-06-18
!>

  INTEGER,INTENT(in)          :: &
    &  jg                           !< patch id
  TYPE(t_var_list_ptr),INTENT(inout) :: &
    &  p_diag_list                  !< current list: diagnostic
  TYPE(t_var_list_ptr),INTENT(in)    :: &
    &  p_prog_list                  !< current list: prognostic
! Local Variables
#ifdef __ART_GPL
  TYPE(t_art_psc),POINTER    :: &
    &  PSC                          !< Pointer to ART PSC fields (polar stratospheric clouds)
#endif
  TYPE(t_art_volc_fplume),POINTER :: & !< Pointer to FPLUME data
    &  fplume
  TYPE(t_art_diag),POINTER    :: &
    &  art_diag                     !< Pointer to ART diagnostic fields
  TYPE(t_cf_var)              :: &
    &  cf_desc                      !< NETCDF variable description
  TYPE(t_grib2_var)           :: &
    &  grib2_desc,               &  !< GRIB2 variable description
    &  grib2_tracer                 !< GRIB2 variable description of tracer
  TYPE(t_mode), POINTER       :: &
    &  this_mode                    !< Pointer to current mode
  TYPE(t_xml_file)            :: &
    &  diag_xmlfile                 !< Diagnostics XML file with GRIB2 meta data
  LOGICAL                     :: &
    &  lexist_diagxml,           &  !< flag if diagnostics XML file exists
    &  l_diag                       !< flag indicating if definition for variable exists in XML file
  LOGICAL                     :: &
    &  art_groups(MAX_GROUPS)       !< group flags for ART diagnostic variables
  INTEGER                     :: &
    &  jt, jp,                   &  !< Photolysis / tracer container index, Pressure layer index
    &  shape2d_c(2),             &  !< Shape for 2D diagnostic variables
    &  shape3d_c(3),             &  !< Shape for 3D diagnostic variables
    &  shape2d_t(3),             &  !< Shape for 2D diagnostic variables, tracer name in 
                                    !    3rd dimension
    &  shape4d_photo(4),         &  !< Shape for photo diagnostic variables
    &  shape4d_reac_rates(4),    &  !< Shape for reac_rate diagnostic variables
    &  shape4d_ihs(4),           &  !< Shape for cgaml diagnostic variables
    &  shape3d_NSB(3),           &  !< Shape for PSC diagnostic variables
    &  shape4d_NSB(4),           &  !< Shape for PSC diagnostic variables
    &  num_radioact                 !< number of radionuclides
  CHARACTER(LEN=3)            :: &
    &  ctrc                         !< Character array to write photolysis container index
  REAL(wp)                    :: &
    &  pres_bot(npreslay),       &  !< bottom of pressure layers [Pa]
    &  pres_top(npreslay),       &  !< top    of pressure layers [Pa]
    &  volc_pres_bot(npreslay),  &  !< bottom of pressure layers for volcanic ash [Pa]
    &  volc_pres_top(npreslay),  &  !< top    of pressure layers for volcanic ash [Pa]
    &  volc_ash_cld_thr(2)          !< array that holds thresholds [ug/m^3] for volanic
                                    !  ash cloud definition
  CHARACTER(LEN=3)            :: &
    &  fl_bot(npreslay),         &  !< associated bottom flight level
    &  fl_top(npreslay),         &  !< associated top    flight level
    &  volc_fl_bot(npreslay),    &  !< associated bottom flight level for volcanic ash
    &  volc_fl_top(npreslay)        !< associated top    flight level for volcanic ash
  CHARACTER(LEN=15)           :: &  !< names of Aeronet wavelengths
    &  tau_dust_shortnames(9),   &
    &  tau_seas_shortnames(9),   &
    &  tau_volc_shortnames(9),   &
    &  tau_soot_shortnames(9)
  CHARACTER(LEN=16)           :: &  !< names of Ceilometer wavelengths
    &  bsc_dust_shortnames(3),   &
    &  bsc_seas_shortnames(3),   &
    &  bsc_volc_shortnames(3),   &
    &  bsc_soot_shortnames(3),   &
    &  ceil_dust_shortnames(3),  &
    &  ceil_seas_shortnames(3),  &
    &  ceil_volc_shortnames(3),  &
    &  ceil_soot_shortnames(3),  &
    &  sat_dust_shortnames(3),   &
    &  sat_seas_shortnames(3),   &
    &  sat_volc_shortnames(3),   &
    &  sat_soot_shortnames(3)

  CHARACTER(LEN=20) :: &
    &  radioact_shortnames(IART_MAX_RADIOACT)  !< Shortnames of radioactive species 
                                               !    e.g. Cs-137 for Caesium137
  CHARACTER(LEN=20) :: &
    &  radio_shortname
  CHARACTER(LEN=1) :: &
    &  radio_search, radio_replace
  CHARACTER(LEN=128)          :: &
    &  var_shortname, var_description
  CHARACTER(LEN=max_timedelta_str_len) :: &
    &  radioact_maxtint
  INTEGER                     :: &
    &  var_class                    !< Variable class, used to set the correct PDT
  TYPE(t_art_chem), POINTER   :: &
    &  art_chem                     !< Pointer to ART diagnostic fields
  TYPE(t_art_external), POINTER   :: &
    &  art_ext                      !< Pointer to ART diagnostic fields
  TYPE(t_art_pollen_table), POINTER         :: &
    &  this_pollen_table            !< pointer to pollen table
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: &
    &  thisroutine = "mo_art_diag_state:art_create_diagnostics"
  INTEGER                         :: &
    &  idx_diag_xml,                 &  !< loop index for diagn. vars
    &  ntrac_diag                       !< number of tracers, for which a certain diagnostics 
                                        !  shall be computed
  CHARACTER(LEN=IART_FILENAMELEN) :: &
    &  x_path_meta
  CHARACTER(LEN=3)                :: &
    &  idx_diag_str
  CHARACTER(LEN=30), ALLOCATABLE  :: &
    & shortnames_diag(:)                !< shortnames of all tracers, for which a certain 
                                        !  diagnostics shall be compted
  INTEGER, ALLOCATABLE            :: &
    & parameterNumbers_diag(:)          !< parameterNumbers of all tracers, for which a certain 
                                        !  diagnostics shall be compted
  CHARACTER(LEN=30)  :: &
    & tracer_name                       !< tracer_name of one tracer, for which a diagnostics 
                                        !  shall be computed
  INTEGER, POINTER                :: &
    & jsp_idx                           !< index of tracer in tracer container
  INTEGER                         :: &
    & plume_ns                          !< number of vertical levels in FPLUME

  INTEGER                         :: &  !< precision of float to be used in netcdf
    & datatype_flt

  TYPE(t_art_atmo),POINTER    :: &
    &  art_atmo           !< Pointer to ART atmospheric fields


  ! Associate Pointers for short references
  art_diag => p_art_data(jg)%diag
  art_chem => p_art_data(jg)%chem
  art_atmo => p_art_data(jg)%atmo

  fplume   => p_art_data(jg)%ext%volc_fplume
#ifdef __ART_GPL
  PSC      => p_art_data(jg)%chem%PSC_meta
#endif
  art_ext  => p_art_data(jg)%ext

  plume_ns = 300

  shape2d_c                 = (/ nproma,                      p_patch(jg)%nblks_c      /)
  shape3d_c                 = (/ nproma, p_patch(jg)%nlev   , p_patch(jg)%nblks_c      /)
  shape4d_photo             = (/ nproma, p_patch(jg)%nlev   , p_patch(jg)%nblks_c , 72 /)
#ifdef __ART_GPL
  shape4d_reac_rates        = (/ nproma, p_patch(jg)%nlev   , p_patch(jg)%nblks_c , NREACT /)
  shape4d_ihs               = (/ nproma, p_patch(jg)%nlev   , p_patch(jg)%nblks_c , ihs_MAX /)
  shape3d_NSB               = (/ nproma,                      p_patch(jg)%nblks_c , PSC%NSB /)
  shape4d_NSB               = (/ nproma, p_patch(jg)%nlev   , p_patch(jg)%nblks_c , PSC%NSB /)
#endif

  ! metadata for aerosol optical depths at specific Aeronet wavelengths
  tau_dust_shortnames   = (/ 'aod_dust_340nm ', 'aod_dust_380nm ', 'aod_dust_440nm ',  &
    &                        'aod_dust_500nm ', 'aod_dust_550nm ', 'aod_dust_675nm ',  &
    &                        'aod_dust_870nm ', 'aod_dust_1020nm', 'aod_dust_1064nm'/)
  tau_seas_shortnames   = (/ 'aod_seas_340nm ', 'aod_seas_380nm ', 'aod_seas_440nm ',  &
    &                        'aod_seas_500nm ', 'aod_seas_550nm ', 'aod_seas_675nm ',  &
    &                        'aod_seas_870nm ', 'aod_seas_1020nm', 'aod_seas_1064nm'/)
  tau_volc_shortnames   = (/ 'aod_ash_340nm  ', 'aod_ash_380nm  ', 'aod_ash_440nm  ',  &
    &                        'aod_ash_500nm  ', 'aod_ash_550nm  ', 'aod_ash_675nm  ',  &
    &                        'aod_ash_870nm  ', 'aod_ash_1020nm ', 'aod_ash_1064nm '/)
  tau_soot_shortnames   = (/ 'aod_soot_340nm ', 'aod_soot_380nm ', 'aod_soot_440nm ',  &
    &                        'aod_soot_500nm ', 'aod_soot_550nm ', 'aod_soot_675nm ',  &
    &                        'aod_soot_870nm ', 'aod_soot_1020nm', 'aod_soot_1064nm'/)

  ! metadata for aerosol backscatter at specific Ceilometer wavelengths
  bsc_dust_shortnames   = (/ 'bsc_dust_355nm  ', 'bsc_dust_532nm  ', 'bsc_dust_1064nm '/)
  bsc_seas_shortnames   = (/ 'bsc_seas_355nm  ', 'bsc_seas_532nm  ', 'bsc_seas_1064nm '/)
  bsc_volc_shortnames   = (/ 'bsc_ash_355nm   ', 'bsc_ash_532nm   ', 'bsc_ash_1064nm  '/)
  bsc_soot_shortnames   = (/ 'bsc_soot_355nm  ', 'bsc_soot_532nm  ', 'bsc_soot_1064nm '/)

  ceil_dust_shortnames  = (/ 'ceil_dust_355nm ', 'ceil_dust_532nm ', 'ceil_dust_1064nm'/)
  ceil_seas_shortnames  = (/ 'ceil_seas_355nm ', 'ceil_seas_532nm ', 'ceil_seas_1064nm'/)
  ceil_volc_shortnames  = (/ 'ceil_ash_355nm  ', 'ceil_ash_532nm  ', 'ceil_ash_1064nm '/)
  ceil_soot_shortnames  = (/ 'ceil_soot_355nm ', 'ceil_soot_532nm ', 'ceil_soot_1064nm'/)

  sat_dust_shortnames   = (/ 'sat_dust_355nm  ', 'sat_dust_532nm  ', 'sat_dust_1064nm '/)
  sat_seas_shortnames   = (/ 'sat_seas_355nm  ', 'sat_seas_532nm  ', 'sat_seas_1064nm '/)
  sat_volc_shortnames   = (/ 'sat_ash_355nm   ', 'sat_ash_532nm   ', 'sat_ash_1064nm  '/)
  sat_soot_shortnames   = (/ 'sat_soot_355nm  ', 'sat_soot_532nm  ', 'sat_soot_1064nm '/)

  ! metadata related to pressure / flight layers for e.g. mineral dust and radionuclides
  pres_bot = (/    -1.0_wp, 84307.0_wp, 69682.0_wp, 59524.0_wp, 50600.0_wp, 37601.0_wp, 23842.0_wp /)
  pres_top = (/ 84307.0_wp, 69682.0_wp, 59524.0_wp, 50600.0_wp, 37601.0_wp, 23842.0_wp, 14748.0_wp /)
  fl_bot   = (/ 'SFC', '050', '100', '140', '180', '250', '350' /)
  fl_top   = (/ '050', '100', '140', '180', '250', '350', '450' /)
  ! metadata related to pressure / flight layers for volcanic ash
  volc_pres_bot   = (/    -1.0_wp, 46500.0_wp, 24000.0_wp,    -1.0_wp, 70000.0_wp, 38500.0_wp, 20000.0_wp /)
  volc_pres_top   = (/ 46500.0_wp, 24000.0_wp,  9100.0_wp, 70000.0_wp, 38500.0_wp, 20000.0_wp, 10000.0_wp /)
  volc_fl_bot     = (/ 'SFC', '200', '350', 'SFC', '100', '245', '390' /)
  volc_fl_top     = (/ '200', '350', '550', '100', '245', '390', '530' /)
  volc_ash_cld_thr = (/2000.0_wp, 4000.0_wp/)

  var_class = CLASS_DEFAULT

  datatype_flt = getNetcdfPrecision()


  ! Inquire existence of XML file with GRIB2 meta data for diagnostic variables
  IF (TRIM(art_config(jg)%cart_diagnostics_xml) /= '') THEN
    INQUIRE(FILE=TRIM(art_config(jg)%cart_diagnostics_xml),EXIST = lexist_diagxml)
  ELSE
    lexist_diagxml = .FALSE.
  END IF

  IF (lexist_diagxml) THEN

    ! Open XML file with GRIB2 meta data for diagnostic variables
    CALL art_open_xml_file(TRIM(art_config(jg)%cart_diagnostics_xml), diag_xmlfile)

    CALL art_get_childnumber_xml(diag_xmlfile, '/diagnostics/', ndiag_xml)

    ALLOCATE(meta_storage(ndiag_xml))

    DO idx_diag_xml = 1, ndiag_xml

      ! Set x_path to current diagnostic variable in XML file
      WRITE(idx_diag_str,'(I3)') idx_diag_xml
      x_path_meta = '/diagnostics/aerosol['//TRIM(ADJUSTL(idx_diag_str))//']/'

      ! Create a new storage container
      CALL meta_storage(idx_diag_xml)%init(.FALSE.)

      ! Read meta data of current diagnostic variable in XML file
      CALL art_read_elements_xml(diag_xmlfile, TRIM(x_path_meta), meta_storage(idx_diag_xml))

    END DO

    ! Close XML file with GRIB2 meta data for diagnostic variables
    CALL art_close_xml_file(diag_xmlfile)

  END IF

  IF (art_config(jg)%lart_aerosol) THEN

    IF (art_config(jg)%iart_radioact == 1) THEN

      ! metadata related to radioactive tracers
      num_radioact = 0

      radio_search = '_'
      radio_replace = '-'

      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          TYPE IS (t_fields_radio)
            num_radioact = num_radioact+1
            IF (num_radioact <= IART_MAX_RADIOACT) THEN
              radio_shortname = TRIM(fields%name)
              radioact_shortnames(num_radioact) = str_replace(radio_shortname,radio_search, &
                &                                             radio_replace)
            ELSE
              CALL finish(thisroutine,'Number of defined radionuclides exceeds IART_MAX_RADIOACT.')
            ENDIF
        END SELECT
        this_mode => this_mode%next_mode
      END DO

    ENDIF

  ENDIF

  !===============================================================================================
  !create new diagnostics for investigating the soiling of pv panels
  !diagnostics may also be used for other art-tracers than dust
  IF (jg == 1) THEN
    ALLOCATE(art_diag_tracer_index(IART_MAX_DIAGCONT, ntracer))
    art_diag_tracer_index(:,:) = 0
  END IF

  ! ----------------------------------
  ! --- 1.0 Aerosol
  ! ----------------------------------


  IF (art_config(jg)%lart_diag_out) THEN
    IF (art_config(jg)%lart_aerosol) THEN

      ALLOCATE (art_diag%so4_sol_aeronet(1))
      cf_desc    = t_cf_var('AOD_550_so4_sol','-', &
        &                     'SO4 sol Optical Depth',datatype_flt)
      grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('AOD_550_so4_sol', grib2_desc, var_class, l_diag)
      CALL add_var( p_diag_list, 'AOD_550_so4_sol', art_diag%so4_sol_aeronet(1)%tau,    &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
        &           ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS"),        &
        &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
        &           fallback_type=HINTP_TYPE_LONLAT_NNB) )

      ALLOCATE (art_diag%ash_insol_aeronet(1))
      cf_desc    = t_cf_var('AOD_550_ash_insol','-', &
          &                     'Ash insol Optical Depth',datatype_flt)
      grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED,GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('AOD_550_ash_insol', grib2_desc, var_class, l_diag)
      CALL add_var( p_diag_list, 'AOD_550_ash_insol',art_diag%ash_insol_aeronet(1)%tau,    &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,var_class=var_class, &
        &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"),        &
        &           hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
        &           fallback_type=HINTP_TYPE_LONLAT_NNB) )

      ALLOCATE (art_diag%ash_mixed_aeronet(1))
      cf_desc    = t_cf_var('AOD_550_ash_mixed','-', &
        &                     'Ash mixed Optical Depth',datatype_flt)
      grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16,GRID_UNSTRUCTURED,GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('AOD_550_ash_mixed', grib2_desc, var_class, l_diag)
      CALL add_var( p_diag_list,'AOD_550_ash_mixed',art_diag%ash_mixed_aeronet(1)%tau,    &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,var_class=var_class, &
        &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"),        &
        &           hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
        &           fallback_type=HINTP_TYPE_LONLAT_NNB) )

      ALLOCATE (art_diag%ash_giant_aeronet(1))
      cf_desc    = t_cf_var('AOD_550_ash_giant','-', &
        &                     'Ash giant Optical Depth',datatype_flt)
      grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED,GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('AOD_550_ash_giant', grib2_desc, var_class, l_diag)
      CALL add_var( p_diag_list, 'AOD_550_ash_giant',art_diag%ash_giant_aeronet(1)%tau,    &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,var_class=var_class, &
        &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"),        &
        &           hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
        &           fallback_type=HINTP_TYPE_LONLAT_NNB) )


      ALLOCATE (art_diag%dust_aeronet(9))
      DO jt=1,9
        ! Mineral Dust Optical Depth
        cf_desc    = t_cf_var(TRIM(tau_dust_shortnames(jt)),'-', &
          &                     'Mineral Dust Optical Depth',datatype_flt)
        grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(tau_dust_shortnames(jt), grib2_desc, var_class, l_diag)
        IF (TRIM(tolower(tau_dust_shortnames(jt))) == 'aod_dust_550nm' .AND. l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, TRIM(tau_dust_shortnames(jt)), art_diag%dust_aeronet(jt)%tau, &
          &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc,grib2_desc,var_class=var_class, &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                    &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,&
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB) )
      ENDDO     ! jt
      DO jt=1,9
        ! Total column of Mineral Dust Optical Depth
        WRITE(var_shortname,  '(A1,A)')  'T', TRIM(tau_dust_shortnames(jt))
        cf_desc    = t_cf_var(TRIM(var_shortname),'-', &
          &                   'Total column of Mineral Dust Optical Depth',datatype_flt)
        grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (TRIM(tolower(var_shortname)) == 'taod_dust_550nm' .AND. l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var( p_diag_list, TRIM(var_shortname), art_diag%dust_aeronet(jt)%tau_vi,              &
          &           GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE, cf_desc, grib2_desc, var_class=var_class, &
          &           ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                          &
          &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,     &
          &                                                  fallback_type=HINTP_TYPE_LONLAT_NNB))
      ENDDO     ! jt
      ALLOCATE (art_diag%dust_ceilo(3))
      DO jt=1,3
        ! Mineral Dust Backscatter
        cf_desc    = t_cf_var(TRIM(bsc_dust_shortnames(jt)),'m-1 sr-1',  &
          &                     'Mineral Dust Backscatter',datatype_flt)
        grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(bsc_dust_shortnames(jt), grib2_desc, var_class, l_diag)
        CALL add_var(p_diag_list, TRIM(bsc_dust_shortnames(jt)), art_diag%dust_ceilo(jt)%bsc,   &
          &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class,  &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS") )
      ENDDO     ! jt
      ALLOCATE (art_diag%dust_att(3))
      DO jt=1,3
        ! Mineral Dust attenuated backscatter (ceilometer)
        cf_desc    = t_cf_var(TRIM(ceil_dust_shortnames(jt)),'m-1 sr-1', &
          &                     'Mineral Dust Attenuated Backscatter Ceilometer',datatype_flt)
        grib2_desc = grib2_var(0, 20, 105, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(ceil_dust_shortnames(jt), grib2_desc, var_class, l_diag)
        IF (TRIM(tolower(ceil_dust_shortnames(jt))) == 'ceil_dust_1064nm' .AND. l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list,TRIM(ceil_dust_shortnames(jt)),art_diag%dust_att(jt)%ceil_bsc, &
          &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class,  &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups )
      ENDDO     ! jt
      ALLOCATE (art_diag%dust_sat(3))
      DO jt=1,3
        ! Mineral Dust attenuated backscatter (satellite)
        cf_desc    = t_cf_var(TRIM(sat_dust_shortnames(jt)),'m-1 sr-1', &
          &                     'Mineral Dust Attenuated Backscatter Satellite',datatype_flt)
        grib2_desc = grib2_var(0, 20, 106, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(sat_dust_shortnames(jt), grib2_desc, var_class, l_diag)
        CALL add_var(p_diag_list,TRIM(sat_dust_shortnames(jt)),art_diag%dust_sat(jt)%sat_bsc,   &
          &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class,  &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS") )
      ENDDO     ! jt
      ! Threshold friction velocity (see Vogel et.al 2006, eq.3.4)
      cf_desc    = t_cf_var('ustar_threshold', 'm s-1', &
           &                'Threshold friction velocity for dust emission', datatype_flt)
      grib2_desc = grib2_var(0, 2, 203, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('ustar_thres', grib2_desc, var_class, l_diag)
      IF (l_diag) THEN
        art_groups = groups("ART_ROUTINE_DIAG")
      ELSE
        art_groups = groups()
      ENDIF
      CALL add_var(p_diag_list, 'ustar_thres', art_diag%ustar_threshold, GRID_UNSTRUCTURED_CELL,&
        &          ZA_SURFACE,  cf_desc, grib2_desc, var_class=var_class,                       &
        &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups )
      ! Friction velocity (JF: not really ART specific -- move to another place?!)
      cf_desc    = t_cf_var('ustar', 'm s-1', &
           &                'Friction velocity', datatype_flt)
      grib2_desc = grib2_var(0, 2, 30, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('ustar', grib2_desc, var_class, l_diag)
      IF (l_diag) THEN
        art_groups = groups("ART_ROUTINE_DIAG")
      ELSE
        art_groups = groups()
      ENDIF
      CALL add_var(p_diag_list, 'ustar', art_diag%ustar, GRID_UNSTRUCTURED_CELL,                &
        &          ZA_SURFACE, cf_desc, grib2_desc, var_class=var_class,                        &
        &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups )
      ! Total mineral dust mass concentration
      cf_desc    = t_cf_var('dust_total_mc', 'kg m-3', &
        &                   'Total mineral dust mass concentration', datatype_flt)
      grib2_desc = grib2_var(0, 20, 0, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('dust_total_mc', grib2_desc, var_class, l_diag)
      IF (l_diag) THEN
        art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
      ELSE
        art_groups = groups("ART_DIAGNOSTICS")
      ENDIF
      CALL add_var( p_diag_list, 'dust_total_mc', art_diag%dust_total_mc,                         &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
        &           ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                       &
        &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
        &                                                  fallback_type=HINTP_TYPE_LONLAT_NNB ), &
        &           vert_interp=create_vert_interp_metadata(vert_intp_type =                      &
        &                                                              vintp_types("P","Z","I"),  &
        &                                                   vert_intp_method=VINTP_METHOD_LIN ) )
      ! Maximum total mineral dust mass concentration between given pressure levels
      ALLOCATE (art_diag%dust_max_total_mc(npreslay))
      DO jp=1, npreslay
        art_diag%dust_max_total_mc(jp)%pres_bot = pres_bot(jp)
        art_diag%dust_max_total_mc(jp)%pres_top = pres_top(jp)
        WRITE(var_shortname,  '(A,A3,A1,A3)')    'dust_max_total_mc_', fl_bot(jp), '_', fl_top(jp)
        WRITE(var_description,'(A,A3,A4,A3,A1)')  &
          &  'Maximum total mineral dust mass concentration between given pressure levels   &
          &  (Flight level ', fl_bot(jp), ' to ', fl_top(jp), ')'
        cf_desc    = t_cf_var(TRIM(var_shortname), 'kg m-3', TRIM(var_description),         &
          &                   datatype_flt)
        grib2_desc = grib2_var(0, 20, 61, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%dust_max_total_mc(jp)%maximum,   &
          &          GRID_UNSTRUCTURED_CELL, ZA_PRES_FL_BOT_TOP, cf_desc, grib2_desc,            &
          &          var_class=var_class, ldims=shape2d_c, lrestart=.FALSE.,                     &
          &          in_group=art_groups,                                                        &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))
      END DO
      ! Total column of mineral dust mass concentration
      cf_desc    = t_cf_var('dust_total_mc_vi', 'kg m-2', &
        &                   'Total column of mineral dust mass concentration', datatype_flt)
      grib2_desc = grib2_var(0, 20, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      IF (lexist_diagxml) &
        &  CALL art_get_GRIB2_diag('dust_total_mc_vi', grib2_desc, var_class, l_diag)
      IF (l_diag) THEN
        art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
      ELSE
        art_groups = groups("ART_DIAGNOSTICS")
      ENDIF
      CALL add_var(p_diag_list, 'dust_total_mc_vi', art_diag%dust_total_mc_vi,                  &
        &          GRID_UNSTRUCTURED_CELL,ZA_ATMOSPHERE,cf_desc,grib2_desc,var_class=var_class, &
        &          ldims=shape2d_c,lrestart=.FALSE.,in_group=art_groups,                        &
        &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
        &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))


      !==============================
      ! ACCUMULATED DRY DEPOSITION: [kg m-2], [# m-2],...
      !==============================
      !----------------------------------
      ! 1) First create 3d variable, which contains all tracers in the rightmost dimension
      !----------------------------------

      ! basic netcdf/grib2 settings:
      cf_desc    = t_cf_var('acc_drydepo', 'tracer-unit m-2',                                     &
           &                'accumulated dry deposition of tracer', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      ! get additional grib2 info from diagnostics.xml file for the data container
      IF (lexist_diagxml) THEN

         CALL art_get_GRIB2_diag('acc_drydepo', grib2_desc, var_class, l_diag,                    &
              & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,            &
              & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN
            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)

            ! add new list element to p_diag_list
            CALL add_var( p_diag_list, 'acc_drydepo', art_diag%acc_drydepo,                       &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                       &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.,            &
                 & lopenacc=.TRUE.)
            __acc_attach(art_diag%acc_drydepo)

            !----------------------------------
            ! 2) Add references for all tracers  art_diag%acc_drydepo(:,:,1/2/...)
            !----------------------------------
            ALLOCATE(art_diag%acc_drydepo_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag
              WRITE(tracer_name,'(A)')shortnames_diag(jt)

              ! reset grib2_desc for each tracer, get tracer specific settings
              grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
              CALL art_get_GRIB2_diag('acc_drydepo', grib2_desc, var_class, l_diag,               &
                &                     itrac = jt, p_prog_list = p_prog_list,                      &
                &                     shortnames_active_tracer = shortnames_diag,                 &
                &                     parameterNumber_list = parameterNumbers_diag,               &
                &                     ntrac_active_tracer = ntrac_diag)

              ! Get index of tracer in tracer container as well as grib2 metadata from tracer
              ! Note: p_prog list could be replaced by p_nh_state_lists(jg)%prog_list(1), 
              !       but this is not working...
              CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

              IF ( ASSOCIATED(jsp_idx) ) THEN
                WRITE(message_text,'(A,A,A,I3,A,I3)')                                  &
                  &    'ART: a reference is added for diagnostic : acc_drydepo_',      &
                  &    TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                CALL message (TRIM(thisroutine), message_text)
                IF (jg == 1) art_diag_tracer_index(IART_ACC_DRYDEPO, jsp_idx) = jt

                ! Take grib2info 'modeNumber' from tracer and write it to diagnostic variable 
                ! (if available)
                CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,   &
                  &                            grib2_desc)
                ! Set grib2 'parameterNumber' for current tracer
                grib2_desc%number = parameterNumbers_diag(jt)

                ! netcdf-settings
                cf_desc = t_cf_var('acc_drydepo_'//TRIM(tracer_name), 'tracer-unit m-2',     &
                  &                'accumulated dry deposition of '//TRIM(tracer_name),      &
                  &                datatype_flt)
                ! For each specified tracer, add list reference to basic data container
                ! art_diag%acc_drydepo; e.g. for tracer 1 (dusta): write data to 
                ! art_diag%acc_drydepo(:,:,1) and specify it in output-namelist as 
                ! acc_drydepo_dusta (c.f. grib2-shortName)

                CALL add_ref( p_diag_list, 'acc_drydepo', 'acc_drydepo_'//TRIM(tracer_name),  &
                  &           art_diag%acc_drydepo_ptr(jt)%p_2d, GRID_UNSTRUCTURED_CELL,      &
                  &           ZA_SURFACE, cf_desc, grib2_desc,                                &
                  &           hor_interp=create_hor_interp_metadata(hor_intp_type =           &
                  &                                                 HINTP_TYPE_LONLAT_BCTR),  &
                  &           ref_idx = jt, var_class=var_class, ldims=shape2d_c,             &
                  &           lrestart=.FALSE.,loutput=.TRUE., isteptype=TSTEP_ACCUM,         &
                  &           in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),          &
                  &           idx_tracer = jsp_idx, idx_diag = jt)

              ELSE
                WRITE(message_text,'(a,a,a)')                                            &
                  &   'ART: WARNING: diagnostic acc_drydepo was requested for tracer: ', &
                  &   tracer_name, ', but such a tracer does not exis.'
                CALL message (TRIM(thisroutine), message_text)

              END IF
            ENDDO
          ELSE
            WRITE(message_text,'(A)') 'ART: WARNING: diagnostic acc_drydepo could not '// &
              &                       'be found in diagnostics.xml file.'
            CALL message (TRIM(thisroutine), message_text)

         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml


!NOTE: check for ncdf output:
!should work with:
!  <aerosol id="acc_drydepo">
!    <tracername_list type="char">dusta,dustb,dustc,dusta0,dustb0,dustc0</tracername_list>
!  </aerosol>

      !==============================
      ! ACCUMULATED SEDIMENTATION: [kg m-2], [# m-2],...
      !==============================
      cf_desc    = t_cf_var('acc_sedim', 'tracer-unit m-2',                                       &
           &                'accumulated sedimentation of tracer', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      IF (lexist_diagxml) THEN
         CALL art_get_GRIB2_diag('acc_sedim', grib2_desc, var_class, l_diag,                      &
              & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,            &
              & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN
            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)

            CALL add_var( p_diag_list, 'acc_sedim', art_diag%acc_sedim,                           &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                       &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.)

            ALLOCATE(art_diag%acc_sedim_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag
               WRITE(tracer_name,'(A)')shortnames_diag(jt)

               ! reset grib2_desc for each tracer, get tracer specific settings
               grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
               CALL art_get_GRIB2_diag('acc_sedim', grib2_desc,var_class,l_diag,itrac = jt,          &
                 &    p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,         &
                 &    parameterNumber_list=parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

               CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

               IF ( ASSOCIATED(jsp_idx) ) THEN
                  WRITE(message_text,'(A,A,A,I3,A,I3)')  & 
                    &     'ART: a reference is added for diagnostic : acc_sedim_',&
                    &     TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                  CALL message (TRIM(thisroutine), message_text)
                  IF (jg == 1) art_diag_tracer_index(IART_ACC_SEDIM, jsp_idx) = jt

                  CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,     &
                    &                            grib2_desc)
                  grib2_desc%number = parameterNumbers_diag(jt)
                  cf_desc = t_cf_var('acc_sedim_'//TRIM(tracer_name), 'tracer-unit m-2',         &
                    &                'accumulated sedimentation of '//TRIM(tracer_name),         &
                    &                datatype_flt)
                  ! For each specified tracer, add list reference to basic data container
                  CALL add_ref(p_diag_list, 'acc_sedim', 'acc_sedim_'//TRIM(tracer_name),       &
                    &          art_diag%acc_sedim_ptr(jt)%p_2d, GRID_UNSTRUCTURED_CELL,         &
                    &          ZA_SURFACE, cf_desc, grib2_desc,                                 &
                    &          hor_interp=create_hor_interp_metadata(hor_intp_type =            &
                    &                                                HINTP_TYPE_LONLAT_BCTR),   &
                    &          ref_idx = jt, var_class=var_class, ldims=shape2d_c,              &
                    &          lrestart=.FALSE.,loutput=.TRUE., isteptype=TSTEP_ACCUM,          &
                    &          in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),           &
                    &          idx_tracer = jsp_idx, idx_diag = jt)
               ELSE
                  WRITE(message_text,'(a,a,a)')                                          &
                    &   'ART: WARNING: diagnostic acc_sedim was requested for tracer: ', &
                    &   tracer_name, ', but such a tracer does not exis.'
                  CALL message (TRIM(thisroutine), message_text)
               END IF
            ENDDO
         ELSE
            WRITE(message_text,'(A)') 'ART: WARNING: diagnostic acc_sedim could not be '// &
              &                       'found in diagnostics.xml file.'
            CALL message (TRIM(thisroutine), message_text)
         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml


      !==============================
      ! ACCUMULATED WET DEPOSITION BY GRID SCALE PRECIPITATION: [kg m-2], [# m-2],...
      !==============================
      cf_desc = t_cf_var('acc_wetdepo_gscp', 'tracer-unit m-2',                                &
        &                'accumulated wet deposition by grid scale precipitation of tracer',   &
        &                datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      IF (lexist_diagxml) THEN
         CALL art_get_GRIB2_diag('acc_wetdepo_gscp', grib2_desc, var_class, l_diag,               &
              & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,            &
              & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN
            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)

            CALL add_var( p_diag_list, 'acc_wetdepo_gscp', art_diag%acc_wetdepo_gscp,             &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                       &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.)

            ALLOCATE(art_diag%acc_wetdepo_gscp_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag
               WRITE(tracer_name,'(A)')shortnames_diag(jt)

               ! reset grib2_desc for each tracer, get tracer specific settings
               grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
               CALL art_get_GRIB2_diag('acc_wetdepo_gscp', grib2_desc, var_class, l_diag,            &
                 &                     itrac = jt, p_prog_list = p_prog_list,                        &
                 &                     shortnames_active_tracer = shortnames_diag,                   &
                 &                     parameterNumber_list = parameterNumbers_diag,                 &
                 &                     ntrac_active_tracer = ntrac_diag)

               CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

               IF ( ASSOCIATED(jsp_idx) ) THEN
                  WRITE(message_text,'(A,A,A,I3,A,I3)')                                 &
                    &   'ART: a reference is added for diagnostic : acc_wetdepo_gscp_', &
                    &   TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                  CALL message (TRIM(thisroutine), message_text)
                  IF (jg == 1) art_diag_tracer_index(IART_ACC_WETDEPO_GSCP, jsp_idx) = jt

                  CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,   &
                    &                            grib2_desc)
                  grib2_desc%number = parameterNumbers_diag(jt)
                  cf_desc = t_cf_var('acc_wetdepo_gscp_'//TRIM(tracer_name),'tracer-unit m-2', &
                    &         'accumulated wet deposition by grid scale precipitation of '//   &
                    &         TRIM(tracer_name), datatype_flt)
                  ! For each specified tracer, add list reference to basic data container
                  CALL add_ref(p_diag_list, 'acc_wetdepo_gscp',                                 &
                    &          'acc_wetdepo_gscp_'//TRIM(tracer_name),                          &
                    &          art_diag%acc_wetdepo_gscp_ptr(jt)%p_2d,                          &
                    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
                    &          hor_interp=create_hor_interp_metadata(hor_intp_type =            &
                    &                                                HINTP_TYPE_LONLAT_BCTR),   &
                    &          ref_idx = jt, var_class=var_class, ldims=shape2d_c,              &
                    &          lrestart=.FALSE.,loutput=.TRUE., isteptype=TSTEP_ACCUM,          &
                    &          in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),           &
                    &          idx_tracer = jsp_idx, idx_diag = jt)
               ELSE
                  WRITE(message_text,'(a,a,a)')  &
                    &   'ART: WARNING: diagnostic acc_wetdepo_gscp was requested for tracer: ', &
                    &   tracer_name, ', but such a tracer does not exist.'
                  CALL message (TRIM(thisroutine), message_text)
               END IF
            ENDDO
         ELSE
           WRITE(message_text,'(A)') 'ART: WARNING: diagnostic acc_wetdepo_gscp could not be '//&
             &                       'found in diagnostics.xml file.'
            CALL message (TRIM(thisroutine), message_text)
         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml


      !==============================
      ! ACCUMULATED WET DEPOSITION BY CONVECTIVE PRECIPITATION: [kg m-2], [# m-2],...
      !==============================
      cf_desc = t_cf_var('acc_wetdepo_con', 'tracer-unit m-2',                                &
        &                'accumulated wet deposition by convective precipitation of tracer',  &
        &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      IF (lexist_diagxml) THEN
         CALL art_get_GRIB2_diag('acc_wetdepo_con', grib2_desc, var_class, l_diag,            &
           & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,           &
           & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN
            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)

            CALL add_var( p_diag_list, 'acc_wetdepo_con', art_diag%acc_wetdepo_con,           &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                   &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.)

            ALLOCATE(art_diag%acc_wetdepo_con_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag
               WRITE(tracer_name,'(A)')shortnames_diag(jt)

               ! reset grib2_desc for each tracer, get tracer specific settings
               grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
               CALL art_get_GRIB2_diag('acc_wetdepo_con', grib2_desc, var_class, l_diag,             &
                 &                     itrac = jt, p_prog_list = p_prog_list,                        &
                 &                     shortnames_active_tracer = shortnames_diag,                   &
                 &                     parameterNumber_list = parameterNumbers_diag,                 &
                 &                     ntrac_active_tracer = ntrac_diag)

               CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

               IF ( ASSOCIATED(jsp_idx) ) THEN
                  WRITE(message_text,'(A,A,A,I3,A,I3)')  &
                    &     'ART: a reference is added for diagnostic : acc_wetdepo_con_',&
                    &     TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                  CALL message (TRIM(thisroutine), message_text)
                  IF (jg == 1) art_diag_tracer_index(IART_ACC_WETDEPO_CON, jsp_idx) = jt

                  CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,    &
                    &                            grib2_desc)
                  grib2_desc%number = parameterNumbers_diag(jt)
                  cf_desc = t_cf_var('acc_wetdepo_con_'//TRIM(tracer_name), 'tracer-unit m-2',    &
                    &       'accumulated wet deposition by convective precipitation of '//        &
                    &       TRIM(tracer_name), datatype_flt)
                  ! For each specified tracer, add list reference to basic data container
                  CALL add_ref(p_diag_list, 'acc_wetdepo_con',                                    &
                    &          'acc_wetdepo_con_'//TRIM(tracer_name),                             &
                    &          art_diag%acc_wetdepo_con_ptr(jt)%p_2d,                             &
                    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                    &          hor_interp=create_hor_interp_metadata(hor_intp_type =              &
                    &                                                HINTP_TYPE_LONLAT_BCTR),     &
                    &          ref_idx = jt, var_class=var_class, ldims=shape2d_c,                &
                    &          lrestart=.FALSE.,loutput=.TRUE., isteptype=TSTEP_ACCUM,            &
                    &          in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),             &
                    &          idx_tracer = jsp_idx, idx_diag = jt)
               ELSE
                  WRITE(message_text,'(a,a,a)')  & 
                    &      'ART: WARNING: diagnostic acc_wetdepo_con was requested for tracer: ', &
                    &      tracer_name, ', but such a tracer does not exis.'
                  CALL message (TRIM(thisroutine), message_text)
               END IF
            ENDDO
         ELSE
           WRITE(message_text,'(A)') 'ART: WARNING: diagnostic acc_wetdepo_con could not be '// &
             &                       'found in diagnostics.xml file.'
            CALL message (TRIM(thisroutine), message_text)
         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml


      !==============================
      ! ACCUMULATED WET DEPOSITION IF PRECIPITATION REACHES SURFACE: [kg m-2], [# m-2],...
      !==============================
      cf_desc = t_cf_var('acc_wetdepo_rrsfc', 'tracer-unit m-2',                                  &
        &                'accumulated wet deposition of tracer if precipitation reaches surface', &
        &                datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      IF (lexist_diagxml) THEN
         CALL art_get_GRIB2_diag('acc_wetdepo_rrsfc', grib2_desc, var_class, l_diag,              &
              & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,            &
              & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN
            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)

            CALL add_var( p_diag_list, 'acc_wetdepo_rrsfc', art_diag%acc_wetdepo_rrsfc,           &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                       &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.)

            ALLOCATE(art_diag%acc_wetdepo_rrsfc_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag

               WRITE(tracer_name,'(A)')shortnames_diag(jt)

               ! reset grib2_desc for each tracer, get tracer specific settings
               grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
               CALL art_get_GRIB2_diag('acc_wetdepo_rrsfc', grib2_desc,                              &
                 &                     var_class, l_diag, itrac = jt, p_prog_list = p_prog_list,     &
                 &                     shortnames_active_tracer = shortnames_diag,                   &
                 &                     parameterNumber_list = parameterNumbers_diag,                 &
                 &                     ntrac_active_tracer = ntrac_diag)

               CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

               IF ( ASSOCIATED(jsp_idx) ) THEN
                  WRITE(message_text,'(A,A,A,I3,A,I3)')                                   &
                    &    'ART: a reference is added for diagnostic : acc_wetdepo_rrsfc_', &
                    &    TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                  CALL message (TRIM(thisroutine), message_text)
                  IF (jg == 1) art_diag_tracer_index(IART_ACC_WETDEPO_RRSFC, jsp_idx) = jt

                  !Take grib2 info from tracer (if available)
                  CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,    &
                    &                            grib2_desc)
                  grib2_desc%number = parameterNumbers_diag(jt)
                  cf_desc = t_cf_var('acc_wetdepo_rrsfc_'//TRIM(tracer_name), 'tracer-unit m-2',  &
                    &                'accumulated wet deposition of '//TRIM(tracer_name)//        &
                    &                ' if precipitation reaches surface', datatype_flt)
                  ! For each specified tracer, add list reference to basic data container
                  CALL add_ref(p_diag_list, 'acc_wetdepo_rrsfc',                                 &
                    &          'acc_wetdepo_rrsfc_'//TRIM(tracer_name),                          &
                    &          art_diag%acc_wetdepo_rrsfc_ptr(jt)%p_2d, GRID_UNSTRUCTURED_CELL,  &
                    &          ZA_SURFACE, cf_desc, grib2_desc,                                  &
                    &          hor_interp=create_hor_interp_metadata(hor_intp_type =             &
                    &                                                HINTP_TYPE_LONLAT_BCTR),    &
                    &          ref_idx = jt, var_class=var_class, ldims=shape2d_c,               &
                    &          lrestart=.FALSE.,loutput=.TRUE., isteptype=TSTEP_ACCUM,           &
                    &          in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),            &
                    &          idx_tracer = jsp_idx, idx_diag = jt)
               ELSE
                 WRITE(message_text,'(a,a,a)') &
                   &   'ART: WARNING: diagnostic acc_wetdepo_rrsfc was requested for tracer: ', &
                   &   tracer_name, ', but such a tracer does not exis.'
                  CALL message (TRIM(thisroutine), message_text)
               END IF
            ENDDO
         ELSE
           WRITE(message_text,'(A)') 'ART: WARNING: diagnostic acc_wetdepo_rrsfc could not be '// &
             &                       'found in diagnostics.xml file.'
            CALL message (TRIM(thisroutine), message_text)
         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml


      !==============================
      ! TRACER EMISSION  [kg m-2 s-1], [# m-2 s-1],...
      !==============================
      cf_desc    = t_cf_var('emiss', 'tracer-unit m-2 s-1',                                      &
           &                'emission of tracer', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      IF (lexist_diagxml) THEN
         CALL art_get_GRIB2_diag('emiss', grib2_desc, var_class, l_diag,                         &
              & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,           &
              & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN
            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)

            CALL add_var( p_diag_list, 'emiss', art_diag%emiss,                                  &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                      &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.)

            ALLOCATE(art_diag%emiss_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag
               WRITE(tracer_name,'(A)')shortnames_diag(jt)

               ! reset grib2_desc for each tracer, get tracer specific settings
               grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
               CALL art_get_GRIB2_diag('emiss', grib2_desc, var_class, l_diag, itrac = jt,           &
                    & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,         &
                    & parameterNumber_list=parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

               CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

               IF ( ASSOCIATED(jsp_idx) ) THEN
                  WRITE(message_text,'(A,A,A,I3,A,I3)')  &
                    &   'ART: a reference is added for diagnostic : emiss_',&
                    &   TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                  CALL message (TRIM(thisroutine), message_text)
                  IF (jg == 1) art_diag_tracer_index(IART_EMISS, jsp_idx) = jt

                  CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,   &
                    &                            grib2_desc)
                  grib2_desc%number = parameterNumbers_diag(jt)
                  cf_desc    = t_cf_var('emiss_'//TRIM(tracer_name), 'tracer-unit m-2 s-1',      &
                    &                   'emission of '//TRIM(tracer_name), datatype_flt)
                  ! For each specified tracer, add list reference to basic data container
                  CALL add_ref(p_diag_list, 'emiss', 'emiss_'//TRIM(tracer_name),                &
                    &          art_diag%emiss_ptr(jt)%p_2d, GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  &
                    &          cf_desc, grib2_desc,                                              &
                    &          hor_interp=create_hor_interp_metadata(hor_intp_type =             &
                    &                                                HINTP_TYPE_LONLAT_BCTR),    &
                    &          ref_idx = jt, var_class=var_class, ldims=shape2d_c,               &
                    &          lrestart=.FALSE.,loutput=.TRUE.,                                  &
                    &          in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),            &
                    &          idx_tracer = jsp_idx, idx_diag = jt)
               ELSE
                  WRITE(message_text,'(a,a,a)')                                      &
                    &   'ART: WARNING: diagnostic emiss was requested for tracer: ', &
                    &   tracer_name, ', but such a tracer does not exis.'
                  CALL message (TRIM(thisroutine), message_text)
               END IF
            ENDDO
         ELSE
           WRITE(message_text,'(A)') 'ART: WARNING: diagnostic emiss could not be '// &
             &                       'found in diagnostics.xml file.'
           CALL message (TRIM(thisroutine), message_text)
         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml



      !==============================
      ! ACCUMULATED TRACER EMISSION  [kg m-2], [# m-2],...
      !==============================
      cf_desc    = t_cf_var('acc_emiss', 'tracer-unit m-2',                                      &
           &                'accumulated emission of tracer', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)

      IF (lexist_diagxml) THEN
         CALL art_get_GRIB2_diag('acc_emiss', grib2_desc, var_class, l_diag,                     &
              & p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,           &
              & parameterNumber_list = parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

         IF (ntrac_diag > 0) THEN

            shape2d_t = (/ nproma, p_patch(jg)%nblks_c , ntrac_diag /)
            CALL add_var( p_diag_list, 'acc_emiss', art_diag%acc_emiss,                          &
                 & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                      &
                 & ldims=shape2d_t,lrestart=.FALSE.,loutput=.FALSE.,lcontainer=.TRUE.)

            ALLOCATE(art_diag%acc_emiss_ptr(ntrac_diag))
            DO jt = 1,ntrac_diag
               WRITE(tracer_name,'(A)')shortnames_diag(jt)

               ! reset grib2_desc for each tracer, get tracer specific settings
               grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
               CALL art_get_GRIB2_diag('acc_emiss', grib2_desc,var_class,l_diag,itrac = jt,          &
                 &  p_prog_list = p_prog_list, shortnames_active_tracer = shortnames_diag,           &
                 &  parameterNumber_list=parameterNumbers_diag, ntrac_active_tracer = ntrac_diag)

               CALL get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2_tracer)

               IF ( ASSOCIATED(jsp_idx) ) THEN
                 WRITE(message_text,'(A,A,A,I3,A,I3)')                                 & 
                   &   'ART: a reference is added for diagnostic : acc_emiss_',        &
                   &   TRIM(tracer_name),', idx_diag = ',jt, ', idx_tracer = ',jsp_idx
                  CALL message (TRIM(thisroutine), message_text)
                  IF (jg == 1) art_diag_tracer_index(IART_ACC_EMISS, jsp_idx) = jt

                  CALL get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('discipline', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('parameterCategory', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('constituentType', grib2_tracer, grib2_desc)
                  CALL get_and_set_int_grib2_key('numberOfModeOfDistribution', grib2_tracer,   &
                    &                            grib2_desc)
                  grib2_desc%number = parameterNumbers_diag(jt)
                  cf_desc = t_cf_var('acc_emiss_'//TRIM(tracer_name), 'tracer-unit m-2 s-1',   &
                    &                'accumulated emission of '//TRIM(tracer_name), datatype_flt)
                  ! For each specified tracer, add list reference to basic data container
                  CALL add_ref(p_diag_list, 'acc_emiss',                                          &
                    &          'acc_emiss_'//TRIM(tracer_name), art_diag%acc_emiss_ptr(jt)%p_2d,  &
                    &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                    &          hor_interp=create_hor_interp_metadata(hor_intp_type =              &
                    &                                                HINTP_TYPE_LONLAT_BCTR),     &
                    &          ref_idx = jt, var_class=var_class, ldims=shape2d_c,                &
                    &          lrestart=.FALSE.,loutput=.TRUE., isteptype=TSTEP_ACCUM,            &
                    &          in_group=groups("ART_DIAGNOSTICS","ART_ROUTINE_DIAG"),             &
                    &          idx_tracer = jsp_idx, idx_diag = jt)
               ELSE
                  WRITE(message_text,'(a,a,a)')                                          &
                    &   'ART: WARNING: diagnostic acc_emiss was requested for tracer: ', &
                    &    tracer_name, ', but such a tracer does not exis.'
                  CALL message (TRIM(thisroutine), message_text)
               END IF
            ENDDO
         ELSE
            WRITE(message_text,'(A)') 'ART: WARNING: diagnostic acc_emiss could not be '// &
              &                       'found in diagnostics.xml file.'
            CALL message (TRIM(thisroutine), message_text)
         END IF ! ntrac_diag > 0
      END IF ! lexist_diagxml



      IF (art_config(jg)%iart_seasalt > 0) THEN
        ALLOCATE (art_diag%seas_aeronet(9))
        DO jt=1,9
          ! Sea salt Optical Depth
          cf_desc    = t_cf_var(TRIM(tau_seas_shortnames(jt)),'-', &
            &                   'Sea Salt Optical Depth',datatype_flt)
          grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(tau_seas_shortnames(jt), grib2_desc, var_class, l_diag)
          IF (TRIM(tolower(tau_seas_shortnames(jt))) == 'aod_seas_550nm' .AND. l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var(p_diag_list, TRIM(tau_seas_shortnames(jt)), art_diag%seas_aeronet(jt)%tau, &
            &          GRID_UNSTRUCTURED_CELL,ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class,   &
            &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,&
            &                                                fallback_type=HINTP_TYPE_LONLAT_NNB) )
        ENDDO     ! jt
        DO jt=1,9
          ! Total column of Sea salt Optical Depth
          WRITE(var_shortname,  '(A1,A)')  'T', TRIM(tau_seas_shortnames(jt))
          cf_desc    = t_cf_var(TRIM(var_shortname),'-', &
            &                   'Total column of Sea Salt Optical Depth',datatype_flt)
          grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
          IF (TRIM(tolower(var_shortname)) == 'taod_seas_550nm' .AND. l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var( p_diag_list, TRIM(var_shortname), art_diag%seas_aeronet(jt)%tau_vi,              &
            &           GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE, cf_desc, grib2_desc, var_class=var_class, &
            &           ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                          &
            &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,     &
            &                                                  fallback_type=HINTP_TYPE_LONLAT_NNB) )
        ENDDO     ! jt
        ALLOCATE (art_diag%seas_ceilo(3))
        DO jt=1,3
          ! Sea salt Backscatter
          cf_desc    = t_cf_var(TRIM(bsc_seas_shortnames(jt)),'m-1 sr-1', &
            &                   'Sea Salt Backscatter',datatype_flt)
          grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            & CALL art_get_GRIB2_diag(bsc_seas_shortnames(jt), grib2_desc, var_class, l_diag)
          CALL add_var(p_diag_list, TRIM(bsc_seas_shortnames(jt)), art_diag%seas_ceilo(jt)%bsc,  &
            &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape3d_c, lrestart=.FALSE.)! in_group=groups("ART_DIAGNOSTICS"))
        ENDDO     ! jt
        ALLOCATE (art_diag%seas_att(3))
        DO jt=1,3
          ! Sea salt attenuated backscatter (ceilometer)
          cf_desc    = t_cf_var(TRIM(ceil_seas_shortnames(jt)),'m-1 sr-1', &
            &                     'Sea Salt Attenuated Backscatter Ceilometer',datatype_flt)
          grib2_desc = grib2_var(0, 20, 105, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(ceil_seas_shortnames(jt), grib2_desc, var_class, l_diag)
          IF (TRIM(tolower(ceil_seas_shortnames(jt))) == 'ceil_seas_1064nm' .AND. l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var( p_diag_list, TRIM(ceil_seas_shortnames(jt)), art_diag%seas_att(jt)%ceil_bsc,  &
            &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
            &           ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups )
        ENDDO     ! jt
        ALLOCATE (art_diag%seas_sat(3))
        DO jt=1,3
          ! Sea salt attenuated backscatter (satellite)
          cf_desc    = t_cf_var(TRIM(sat_seas_shortnames(jt)),'m-1 sr-1', &
            &                     'Sea Salt Attenuated Backscatter Satellite',datatype_flt)
          grib2_desc = grib2_var(0, 20, 106, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            & CALL art_get_GRIB2_diag(sat_seas_shortnames(jt), grib2_desc, var_class, l_diag)
          CALL add_var(p_diag_list,TRIM(sat_seas_shortnames(jt)),art_diag%seas_sat(jt)%sat_bsc,  &
            &          GRID_UNSTRUCTURED_CELL,ZA_HYBRID,cf_desc,grib2_desc, var_class=var_class, &
            &          ldims=shape3d_c, lrestart=.FALSE. )! in_group=groups("ART_DIAGNOSTICS"))
        ENDDO     ! jt
        ! Total sea salt mass concentration
        cf_desc    = t_cf_var('seas_total_mc', 'kg m-3', &
          &                   'Total sea salt mass concentration', datatype_flt)
        grib2_desc = grib2_var(0, 20, 0, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag('seas_total_mc', grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var( p_diag_list, 'seas_total_mc', art_diag%seas_total_mc,                           &
          &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class,   &
          &           ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                         &
          &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,    &
          &                                                  fallback_type=HINTP_TYPE_LONLAT_NNB ),   &
          &           vert_interp=create_vert_interp_metadata(vert_intp_type =                        &
          &                                                              vintp_types("P","Z","I"),    &
          &                                                   vert_intp_method=VINTP_METHOD_LIN ) )
      ENDIF !iart_seasalt

      ALLOCATE (art_diag%volc_aeronet(9))
      DO jt=1,9
        ! Volcanic ash Optical Depth
        cf_desc    = t_cf_var(TRIM(tau_volc_shortnames(jt)),'-', &
          &                     'Volcanic Ash Optical Depth',datatype_flt)
        grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (TRIM(tolower(tau_volc_shortnames(jt))) == 'aod_ash_550nm' .AND. l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(tau_volc_shortnames(jt), grib2_desc, var_class, l_diag)
        CALL add_var(p_diag_list, TRIM(tau_volc_shortnames(jt)), art_diag%volc_aeronet(jt)%tau, &
          &          GRID_UNSTRUCTURED_CELL,ZA_HYBRID,cf_desc,grib2_desc, var_class=var_class,  &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                    &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,&
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))
      ENDDO     ! jt
      DO jt=1,9
        ! Total column of Volcanic ash Optical Depth
        WRITE(var_shortname,  '(A1,A)')  'T', TRIM(tau_volc_shortnames(jt))
        cf_desc    = t_cf_var(TRIM(var_shortname),'-', &
          &                   'Total column of Volcanic Ash Optical Depth',datatype_flt)
        grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (TRIM(tolower(var_shortname)) == 'taod_ash_550nm' .AND. l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var( p_diag_list, TRIM(var_shortname), art_diag%volc_aeronet(jt)%tau_vi,              &
          &           GRID_UNSTRUCTURED_CELL, ZA_ATMOSPHERE, cf_desc, grib2_desc, var_class=var_class, &
          &           ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                          &
          &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,     &
          &                                                  fallback_type=HINTP_TYPE_LONLAT_NNB))
      ENDDO     ! jt
      ALLOCATE (art_diag%volc_ceilo(3))
      DO jt=1,3
        ! Volcanic ash backscatter
        cf_desc    = t_cf_var(TRIM(bsc_volc_shortnames(jt)),'m-1 sr-1', &
          &                     'Volcanic Ash Backscatter',datatype_flt)
        grib2_desc = grib2_var(0, 20, 103, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(bsc_volc_shortnames(jt), grib2_desc, var_class, l_diag)
        CALL add_var(p_diag_list, TRIM(bsc_volc_shortnames(jt)), art_diag%volc_ceilo(jt)%bsc,   &
          &          GRID_UNSTRUCTURED_CELL,ZA_HYBRID,cf_desc, grib2_desc, var_class=var_class, &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS") )
      ENDDO     ! jt
      ALLOCATE (art_diag%volc_att(3))
      DO jt=1,3
        ! Volcanic ash attenuated backscatter (ceilometer)
        cf_desc    = t_cf_var(TRIM(ceil_volc_shortnames(jt)),'m-1 sr-1', &
          &                     'Volcanic Ash Attenuated Backscatter Ceilometer',datatype_flt)
        grib2_desc = grib2_var(0, 20, 105, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(ceil_volc_shortnames(jt), grib2_desc, var_class, l_diag)
        IF (TRIM(tolower(ceil_volc_shortnames(jt))) == 'ceil_ash_1064nm' .AND. l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list,TRIM(ceil_volc_shortnames(jt)),art_diag%volc_att(jt)%ceil_bsc, &
          &          GRID_UNSTRUCTURED_CELL,ZA_HYBRID,cf_desc, grib2_desc, var_class=var_class, &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups )
      ENDDO     ! jt
      ALLOCATE (art_diag%volc_sat(3))
      DO jt=1,3
        ! Volcanic ash attenuated backscatter (satellite)
        cf_desc    = t_cf_var(TRIM(sat_volc_shortnames(jt)),'m-1 sr-1', &
          &                     'Volcanic Ash Attenuated Backscatter Satellite',datatype_flt)
        grib2_desc = grib2_var(0, 20, 106, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          & CALL art_get_GRIB2_diag(sat_volc_shortnames(jt), grib2_desc, var_class, l_diag)
        CALL add_var(p_diag_list, TRIM(sat_volc_shortnames(jt)), art_diag%volc_sat(jt)%sat_bsc, &
          &          GRID_UNSTRUCTURED_CELL,ZA_HYBRID,cf_desc, grib2_desc, var_class=var_class, &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS") )
      ENDDO     ! jt
  
      IF (art_config(jg)%iart_fire > 0) THEN
        ALLOCATE (art_diag%soot_aeronet(9))
        DO jt=1,9
          ALLOCATE (art_diag%soot_aeronet(jt)%tau(nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c))
          ! Soot Optical Depth
          cf_desc    = t_cf_var(TRIM(tau_soot_shortnames(jt)),'-', &
            &                     'Soot Optical Depth',datatype_flt)
          grib2_desc = grib2_var(0, 20, 102, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (TRIM(tolower(tau_soot_shortnames(jt))) == 'aod_soot_550nm' .AND. l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(tau_soot_shortnames(jt), grib2_desc, var_class, l_diag)
          CALL add_var( p_diag_list, TRIM(tau_soot_shortnames(jt)), art_diag%soot_aeronet(jt)%tau,    &
            &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
            &           ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                       &
            &           hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
            &                                                  fallback_type=HINTP_TYPE_LONLAT_NNB))
        ENDDO     ! jt
        ALLOCATE (art_diag%soot_ceilo(3))
        DO jt=1,3
          ALLOCATE (art_diag%soot_ceilo(jt)%bsc(nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c))
          ! Soot ash backscatter
          cf_desc    = t_cf_var(TRIM(bsc_soot_shortnames(jt)),'m-1 sr-1', &
            &                     'Soot Backscatter',datatype_flt)
          grib2_desc = grib2_var(0, 20, 103, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(bsc_soot_shortnames(jt), grib2_desc, var_class, l_diag)
          CALL add_var( p_diag_list, TRIM(bsc_soot_shortnames(jt)), art_diag%soot_ceilo(jt)%bsc,      &
            &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
            &           ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS") )
        ENDDO     ! jt
        ALLOCATE (art_diag%soot_att(3))
        DO jt=1,3
          ALLOCATE (art_diag%soot_att(jt)%ceil_bsc(nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c))
          ! Soot attenuated backscatter (ceilometer)
          cf_desc    = t_cf_var(TRIM(ceil_soot_shortnames(jt)),'m-1 sr-1', &
            &                     'Soot Attenuated Backscatter Ceilometer',datatype_flt)
          grib2_desc = grib2_var(0, 20, 105, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(ceil_soot_shortnames(jt), grib2_desc, var_class, l_diag)
          IF (TRIM(tolower(ceil_soot_shortnames(jt))) == 'ceil_soot_1064nm' .AND. l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var( p_diag_list, TRIM(ceil_soot_shortnames(jt)), art_diag%soot_att(jt)%ceil_bsc,  &
            &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
            &           ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups )
        ENDDO     ! jt
        ALLOCATE (art_diag%soot_sat(3))
        DO jt=1,3
          ALLOCATE (art_diag%soot_sat(jt)%sat_bsc(nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c))
          ! Soot attenuated backscatter (satellite)
          cf_desc    = t_cf_var(TRIM(sat_soot_shortnames(jt)),'m-1 sr-1', &
            &                     'Soot Attenuated Backscatter Satellite',datatype_flt)
          grib2_desc = grib2_var(0, 20, 106, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(sat_soot_shortnames(jt), grib2_desc, var_class, l_diag)
          CALL add_var( p_diag_list, TRIM(sat_soot_shortnames(jt)), art_diag%soot_sat(jt)%sat_bsc,    &
            &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
            &           ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS") )
        ENDDO     ! jt
      ENDIF ! iart_fire
 

      this_mode => p_mode_state(jg)%p_mode_list%p%first_mode

      DO WHILE(ASSOCIATED(this_mode))
        ! Select type of mode
        SELECT TYPE (fields=>this_mode%fields)
          CLASS is (t_fields_2mom)
            cf_desc    = t_cf_var('diam_'//TRIM(this_mode%fields%name), 'm', &
              &                   'Diameter of mode '//TRIM(this_mode%fields%name), datatype_flt)
            grib2_desc = grib2_var(0, 254, 201, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
            CALL add_var(p_diag_list, 'diam_'//TRIM(this_mode%fields%name), fields%diameter,           &
              &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, &
              &          ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS"))
          CLASS DEFAULT
            ! Nothing to do here
        END SELECT
        this_mode => this_mode%next_mode
      ENDDO

      IF (art_config(jg)%iart_pollen > 0) THEN
        CALL art_extinit_pollen_fordiag(p_patch(jg),p_patch(jg)%nblks_c, art_ext%pollen_prop, &
          &                     p_art_data(jg)%dict_tracer,TRIM(art_config(jg)%cart_input_folder))

        DO jt=1, art_ext%pollen_prop%npollen_types
          this_pollen_table => art_ext%pollen_prop%pollen_type(jt)
          ! check if pollen type is active -> leave loop
          IF (.NOT.ASSOCIATED (this_pollen_table%fr_cov) .OR.  &
            & TRIM(this_pollen_table%shortname) == '') EXIT
  
          ! precipitation reservoir (liquid water on the flowers, preventing them from flowering)
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'rprec'
          WRITE(var_description,'(A,A)') 'Precipitation reservoir of ', &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%r_precip(nproma, p_patch(jg)%nblks_c))
          cf_desc   = t_cf_var(TRIM(var_shortname), 'm-2', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(2, 0, 13, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%r_precip,              &
            &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%r_precip)
          ! number of pollen in the reservoir (previous timestep)
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'reso'
          WRITE(var_description,'(A,A)') 'Pollen reservoir (previous timestep) of ',              &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%res_old(nproma, p_patch(jg)%nblks_c))
          cf_desc   = t_cf_var(TRIM(var_shortname), 'm-2', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 192, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%res_old,               &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%res_old)
          ! sum of released pollen into the reservoir (daily sum)
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'ress'
          WRITE(var_description,'(A,A)') 'Pollen reservoir (daily sum) of ',                      &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%res_new_sum(nproma, p_patch(jg)%nblks_c))
          cf_desc   = t_cf_var(TRIM(var_shortname), 'm-2', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 193, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%res_new_sum,           &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          isteptype=TSTEP_ACCUM,                                                     &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%res_new_sum)
          ! state of pollen season (eq. zero before and after season, the higher, 
          !                         the more plants are flowering)
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'sdes'
          WRITE(var_description,'(A,A)') 'State of pollen season of ',            &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%f_q_seas(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), '-', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 194, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%f_q_seas,              &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%f_q_seas)
          ! cumulated weighted temperature sum of daily (12 UTC) values
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'ctsum'
          WRITE(var_description,'(A,A)') 'Cumulated weighted 2m temperature sum of ',             &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%ctsum(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), 'K', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 196, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%ctsum,                 &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups )
          ! number of days since the start of pollen season
          ! (if present day is in the season: zero outside season)
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'saisn'
          WRITE(var_description,'(A,A,A)') 'Number of days since start of pollen season of ',     &
            &                             TRIM(this_pollen_table%shortname),                      &
            &                             ' (if present day is in the season: zero outside season)'
          ALLOCATE (this_pollen_table%saisn(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), 'd', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 200, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%saisn,                &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups )
          ! length of pollen season
          IF (TRIM(this_pollen_table%shortname) /= 'POAC') THEN
            WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'saisl'
            WRITE(var_description,'(A,A)') 'Length of pollen season of ',  &
              &                            TRIM(this_pollen_table%shortname)
            ALLOCATE (this_pollen_table%saisl(nproma, p_patch(jg)%nblks_c))
            cf_desc   = t_cf_var(TRIM(var_shortname), 'd', TRIM(var_description), datatype_flt)
            grib2_desc = grib2_var(0, 20, 201, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
            IF (lexist_diagxml) &
              &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
            IF (l_diag) THEN
              art_groups = groups("ART_ROUTINE_DIAG")
            ELSE
              art_groups = groups()
            ENDIF
            CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%saisl,               &
              &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,&
              &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups )
          ENDIF
          ! pollen number emission flux
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'fe'
          WRITE(var_description,'(A,A)') 'Emission flux of ', TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%fe_plant(nproma, p_patch(jg)%nblks_c))
          cf_desc = t_cf_var(TRIM(var_shortname),'m-2 s-1',TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 202, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%fe_plant,             &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%fe_plant)
          ! number of days since the start of pollen season
          ! (if present day is outside the season: length of current season)
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'saisa'
          WRITE(var_description,'(A,A,A)') 'Number of days since the start of pollen season of ', &
            &                              TRIM(this_pollen_table%shortname),                     &
            &                              ' (if present day is outside the season: '//           &
            &                              'length of current season)'
          ALLOCATE (this_pollen_table%saisa(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), 'd', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 200, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%saisa,                &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups )
          ! field of pollen tuning factors
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'tune'
          WRITE(var_description,'(A,A)') 'Pollen tuning factor of ',                              &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%tune(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), '-', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 211, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%tune,                  &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%tune)
          ! Cumulated 2m temperature sum threshold for the start of the pollen season
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'tthrs'
          WRITE(var_description,'(A,A)') 'Cumulated 2m temperature sum threshold ' //             &
                                         'for the start of the pollen season of ',                &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%tthrs(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), '-', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 197, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%tthrs,                  &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%tthrs)
          ! Cumulated 2m temperature sum threshold for the end of the pollen season
          WRITE(var_shortname,  '(A,A)') TRIM(this_pollen_table%shortname), 'tthre'
          WRITE(var_description,'(A,A)') 'Cumulated 2m temperature sum threshold ' //             &
                                         'for the end of the pollen season of ',                &
            &                            TRIM(this_pollen_table%shortname)
          ALLOCATE (this_pollen_table%tthre(nproma, p_patch(jg)%nblks_c))
          cf_desc    = t_cf_var(TRIM(var_shortname), '-', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 198, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(TRIM(var_shortname), grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups()
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), this_pollen_table%tthre,                  &
            &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,  &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          lopenacc=.TRUE.)
          __acc_attach(this_pollen_table%tthre)

        ENDDO ! jt

      END IF ! lart_aerosol (pollen)
    END IF ! lart_aerosol


    ! ----------------------------------
    ! --- 2.0 Chemistry
    ! ----------------------------------

    IF (art_config(jg)%lart_chem) THEN
      cf_desc    = t_cf_var('art_o3', 'kg/kg', 'art_o3', datatype_flt)
      grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
 
      CALL add_var( p_diag_list, 'art_o3', art_atmo%o3_ext,                  &
             &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
             &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
             &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
             &             ldims=shape3d_c,                                          &
             &             lcontainer=.False., lrestart=.FALSE., loutput=.TRUE. )


      cf_desc    = t_cf_var('OH_Nconc', '# / cm3', 'OH_Nconc', datatype_flt)
      grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
   
      CALL add_var( p_diag_list, 'OH_Nconc', art_chem%param%OH_chem_meta%OH_Nconc,       &
             &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
             &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
             &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
             &             ldims=shape3d_c,                                          &
             &             lcontainer=.False., lrestart=.FALSE., loutput=.TRUE. )


      IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
        cf_desc    = t_cf_var('reaction_rates', '1/s', 'MECCA reaction rates', datatype_flt)
        grib2_desc = grib2_var(0, 254, 218, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_diag_list, 'reac_rates', art_chem%mecicon%utils%reac_rates,  &
          &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,      &
          &             vert_interp=create_vert_interp_metadata(vert_intp_type =     &
          &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN), &
          &             ldims=shape4d_reac_rates,                                    &
          &             lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )

        ALLOCATE(art_diag%art_reac_rates_ptr(NREACT))
        DO jt =1,NREACT
          WRITE(ctrc,'(I3.3)')jt
          CALL add_ref( p_diag_list, 'reac_rates', 'reac_rate'//TRIM(ADJUSTL(ctrc)),        &
            &             art_diag%art_reac_rates_ptr(jt)%p_3d,                             &
            &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                                &
            &             cf_desc, grib2_desc,ref_idx=jt, ldims=shape3d_c,lrestart=.FALSE., &
            &             vert_interp=create_vert_interp_metadata(                          &
            &             vert_intp_type=vintp_types("P","Z","I"),                          &
            &             vert_intp_method=VINTP_METHOD_LIN),                               &
            &             in_group=groups("ART_DIAGNOSTICS"))
        ENDDO
#endif
      END IF
    END IF

    ! ----------------------------------
    ! --- 3.0 Cloud Microphysics and Aerosol-Cloud-Interactions
    ! ----------------------------------

    IF (art_config(jg)%iart_aci_warm > 0) THEN
      ! Number of calls of activation routine (accumulated)
      cf_desc    = t_cf_var('ncalls_warm', '#', &
        &                   'Accumulated number of calls of activation routine', datatype_flt)
      grib2_desc = grib2_var(0, 254, 210, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ncalls_warm', art_diag%ncalls_warm,             &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,      &
        &           ldims=shape3d_c, lrestart=.FALSE.)
      ! Number of activated particles (accumulated, warm-phase)
      cf_desc    = t_cf_var('aci_nnuctot_warm', '# m-3', &
        &                   'Accumulated number of activated particles', datatype_flt)
      grib2_desc = grib2_var(0, 254, 211, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'aci_nnuctot_warm', art_diag%aci_nnuctot_warm,             &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,                &
        &           ldims=shape3d_c, lrestart=.FALSE.)
      ! Number of activated particles according to FHH theory (accumulated, warm-phase)
      cf_desc    = t_cf_var('aci_nnucfhh_warm', '# m-3', &
        &                   'Accumulated number of activated particles FHH theory', datatype_flt)
      grib2_desc = grib2_var(0, 254, 212, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'aci_nnucfhh_warm', art_diag%aci_nnucfhh_warm,             &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,                &
        &           ldims=shape3d_c, lrestart=.FALSE.)
      ! Maximum supersaturation over liquid water
      cf_desc    = t_cf_var('smax_water', 'nounit', &
        &                   'Maximum supersaturation over liquid water', datatype_flt)
      grib2_desc = grib2_var(0, 254, 213, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'smax_water', art_diag%smax_water,                         &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,                &
        &           ldims=shape3d_c, lrestart=.FALSE.)
    END IF !iart_aci_warm > 0

    IF (art_config(jg)%iart_aci_cold > 0) THEN
      ! Number of calls of ice nucleation routine (accumulated)
      cf_desc    = t_cf_var('ncalls_cold', '#', &
        &                   'Accumulated number of calls of ice nucleation routine',datatype_flt)
      grib2_desc = grib2_var(0, 254, 214, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ncalls_cold', art_diag%ncalls_cold,             &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,      &
        &           ldims=shape3d_c, lrestart=.FALSE.)
      ! Number of total nucleated (hom. freez. + het. nuc.) particles (accumulated, cold-phase)
      cf_desc    = t_cf_var('aci_nnuctot_cold', '# m-3', &
        &                   'Accumulated number of total nucleated particles', datatype_flt)
      grib2_desc = grib2_var(0, 254, 215, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'aci_nnuctot_cold', art_diag%aci_nnuctot_cold,         &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,            &
        &           ldims=shape3d_c, lrestart=.FALSE.)
      ! Number of het. nucleated particles (accumulated, cold-phase)
      cf_desc    = t_cf_var('aci_nnuchet_cold', '# m-3', &
        &                   'Accumulated number of het. nucleated particles', datatype_flt)
      grib2_desc = grib2_var(0, 254, 216, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'aci_nnuchet_cold', art_diag%aci_nnuchet_cold,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,          &
        &           ldims=shape3d_c, lrestart=.FALSE.)
      ! Maximum supersaturation over ice
      cf_desc    = t_cf_var('smax_ice', 'nounit', &
        &                   'Maximum supersaturation over ice', datatype_flt)
      grib2_desc = grib2_var(0, 254, 217, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'smax_ice', art_diag%smax_ice,                 &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,    &
        &           ldims=shape3d_c, lrestart=.FALSE.)
    ENDIF !iart_aci_cold > 0

    ! ----------------------------------
    ! --- 4.0 Debug
    ! ----------------------------------

    ! ART DEBUG 3D
    cf_desc    = t_cf_var('art_3d_dbg01', 'nounit', 'Art 3D Debug Variable 01', datatype_flt)
    grib2_desc = grib2_var(0, 254, 245, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ART_3D_DBG01', art_diag%art_3d_dbg01,           &
      &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,      &
      &           ldims=shape3d_c, lrestart=.FALSE.)

    ! ART DEBUG 2D
    cf_desc    = t_cf_var('art_2d_dbg01', 'nounit', 'Art 2D Debug Variable 01', datatype_flt)
    grib2_desc = grib2_var(0, 254, 246, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ART_2D_DBG01', art_diag%art_2d_dbg01,           &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
      &           ldims=shape2d_c, lrestart=.FALSE.)

    ! ART DEBUG 3D
    ! GV: This call to add_var resulted in the error "double free or corruption"
    !     when compiling icon with art enabled. Since the variables here are not
    !     used anywhere else, the code is commented out.
    ! cf_desc    = t_cf_var('art_3d_dbg_init', 'nounit', 'Art 3D Debug INIT', datatype_flt)
    ! grib2_desc = grib2_var(0, 254, 245, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
    ! CALL add_var( p_diag_list, 'ART_3D_DBG_INIT', art_diag%art_3d_dbg_init,           &
    !   &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,      &
    !   &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
    !     &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
    !     &             ldims=shape3d_c,                                          &
    !     &             lrestart=.FALSE., loutput=.TRUE.)


    ! ----------------------------------
    ! --- 5.0 Tropopause diagnosticcs
    ! ----------------------------------

    ! MS: Note that the following two are also allocated in mo_art_collect_atmo_state.f90:art_collect_atmo_state, so it 
    !     should be decided which code to keep.

    cf_desc    = t_cf_var('art_ptropo', 'Pa', 'tropopause height', datatype_flt)
    grib2_desc = grib2_var(0, 254, 246, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'art_ptropo', art_atmo%ptropo,           &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
      &           ldims=shape2d_c, lrestart=.FALSE., lopenacc=.TRUE.)
    __acc_attach(art_atmo%ptropo)


    cf_desc    = t_cf_var('art_ktrpwmop', 'nounit', 'level of tropopause', datatype_flt)
    grib2_desc = grib2_var(0, 254, 246, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'art_ktrpwmop', art_atmo%ktrpwmop1_real,         &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
      &           ldims=shape2d_c, lrestart=.FALSE.) !, lopenacc=.TRUE.)  ! This variable is already allocated on GPU in
    !__acc_attach(art_atmo%ktrpwmop1_real)                                ! shared/mo_art_collect_atmo_state.f90:art_collect_atmo_state

    ! ----------------------------------
    ! --- 6.0 Volcanic ash diagnostics
    ! ----------------------------------
    IF (art_config(jg)%lart_aerosol) THEN
      IF (art_config(jg)%iart_volcano > 0) THEN
        ! Total volcanic ash mass concentration
        cf_desc    = t_cf_var('ash_total_mc', 'kg m-3', &
          &                   'Total volcanic ash mass concentration', datatype_flt)
        grib2_desc = grib2_var(0, 20, 0, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag('ash_total_mc', grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, 'ash_total_mc', art_diag%ash_total_mc, GRID_UNSTRUCTURED_CELL,  &
          &          ZA_HYBRID,  cf_desc, grib2_desc, var_class=var_class, ldims=shape3d_c,       &
          &          lrestart=.FALSE., in_group=art_groups,                                       &
          &          hor_interp=create_hor_interp_metadata( hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
          &                                                 fallback_type=HINTP_TYPE_LONLAT_NNB ))
        ! Maximum total volcanic ash mass concentration between given pressure levels
        ALLOCATE (art_diag%ash_max_total_mc(npreslay))
        DO jp=1, npreslay
          art_diag%ash_max_total_mc(jp)%pres_bot = volc_pres_bot(jp)
          art_diag%ash_max_total_mc(jp)%pres_top = volc_pres_top(jp)
          WRITE(var_shortname,  '(A,A3,A1,A3)')    'ash_max_total_mc_', volc_fl_bot(jp), '_', volc_fl_top(jp)
          WRITE(var_description,'(A,A3,A4,A3,A1)')  &
            &  'Maximum total volcanic ash mass concentration between given pressure levels   &
            &  (Flight level ', volc_fl_bot(jp), ' to ', volc_fl_top(jp), ')'
          cf_desc    = t_cf_var(TRIM(var_shortname), 'kg m-3', TRIM(var_description),         &
            &                   datatype_flt)
          grib2_desc = grib2_var(0, 20, 61, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%ash_max_total_mc(jp)%maximum,    &
            &          GRID_UNSTRUCTURED_CELL, ZA_PRES_FL_BOT_TOP, cf_desc,grib2_desc,             &
            &          var_class=var_class, ldims=shape2d_c, lrestart=.FALSE.,                     &
            &          in_group=art_groups,                                                        &
            &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
            &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))
        END DO
        ! Total column of volcanic ash mass concentration
        cf_desc    = t_cf_var('ash_total_mc_vi', 'kg m-2', &
          &                   'Total column of volcanic ash mass concentration', datatype_flt)
        grib2_desc = grib2_var(0, 20, 1, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag('ash_total_mc_vi', grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, 'ash_total_mc_vi', art_diag%ash_total_mc_vi,                    &
          &          GRID_UNSTRUCTURED_CELL,ZA_ATMOSPHERE,cf_desc,grib2_desc,var_class=var_class, &
          &          ldims=shape2d_c,lrestart=.FALSE.,in_group=art_groups,                        &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))
        ! Height of maximal total volcanic ash mass concentration
        cf_desc    = t_cf_var('ash_hml_max', 'm', &
          &                   'Height of maximal total volcanic ash mass concentration', &
          &                   datatype_flt)
        grib2_desc = grib2_var(0, 20, 62, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag('ash_hml_max', grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, 'ash_hml_max', art_diag%ash_hml_max,                            &
          &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,   &
          &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                      &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))
        ! Height of volvanic ash cloud base/top for different threshold concentrations
        ALLOCATE (art_diag%ash_cloud(2))
        DO jt=1, 2
          art_diag%ash_cloud(jt)%threshold = volc_ash_cld_thr(jt)
          WRITE(var_shortname, '(A,I4)')  'ash_cloud_base_',INT(art_diag%ash_cloud(jt)%threshold)
          WRITE(var_description, '(A,I4,A7)') &
            & 'lowest height in meter where the mass concentration &
            &of ash is higher than ',INT(art_diag%ash_cloud(jt)%threshold),'ug/m3'
          cf_desc    = t_cf_var(TRIM(var_shortname), 'm', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 63, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%ash_cloud(jt)%base,             &
            &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,&
            &                                                fallback_type=HINTP_TYPE_LONLAT_NNB) )
          WRITE(var_shortname, '(A,I4)')  'ash_cloud_top_',INT(art_diag%ash_cloud(jt)%threshold)
          WRITE(var_description, '(A,I4,A7)') &
            & 'greatest height in meter where the mass concentration &
            &of ash is higher than ',INT(art_diag%ash_cloud(jt)%threshold),'ug/m3'
          cf_desc    = t_cf_var(TRIM(var_shortname), 'm', TRIM(var_description), datatype_flt)
          grib2_desc = grib2_var(0, 20, 64, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%ash_cloud(jt)%top,              &
            &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class, &
            &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                    &
            &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,&
            &                                                fallback_type=HINTP_TYPE_LONLAT_NNB) )
        END DO
      ENDIF
    ENDIF !lart_aerosol

  ENDIF !lart_diag_out

  ! ----------------------------------
  ! --- 7.0 Radioactive tracers diagnostics
  ! ----------------------------------
  IF (art_config(jg)%lart_aerosol) THEN
    IF (art_config(jg)%iart_radioact == 1) THEN
      CALL getPTStringFromMS(NINT(1000*art_config(jg)%radioact_maxtint, i8), radioact_maxtint)
      ! Wet and dry deposition values will be calculated independent 
      ! from art_config(jg)%lart_diag_out
      ALLOCATE (art_diag%radioact(num_radioact))
      DO jt=1, num_radioact
        ! Accumulated wet deposition of radionuclides
        WRITE(var_shortname,  '(A,A1)')  TRIM(radioact_shortnames(jt)), 'w'
        WRITE(var_description,'(A,A)')   TRIM(radioact_shortnames(jt)), ' - wet deposition'
        cf_desc    = t_cf_var(TRIM(var_shortname), 'Bq m-2', TRIM(var_description), datatype_flt)
        grib2_desc = grib2_var(0, 18, 11, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%radioact(jt)%wetdepo,             &
          &          GRID_UNSTRUCTURED_CELL, ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,   &
          &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                      &
          &          isteptype=TSTEP_ACCUM,                                                       &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB ))
        ! Accumulated dry deposition of radionuclides
        WRITE(var_shortname,  '(A,A1)')  TRIM(radioact_shortnames(jt)), 'd'
        WRITE(var_description,'(A,A)')   TRIM(radioact_shortnames(jt)), ' - dry deposition'
        cf_desc    = t_cf_var(TRIM(var_shortname), 'Bq m-2', TRIM(var_description), datatype_flt)
        grib2_desc = grib2_var(0, 18, 12, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%radioact(jt)%drydepo,             &
          &          GRID_UNSTRUCTURED_CELL,ZA_SURFACE,cf_desc,grib2_desc,var_class=var_class,    &
          &          ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,                      &
          &          isteptype=TSTEP_ACCUM,                                                       &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB ))
        ! Averaged air concentration [Bq m-3] of radionuclides
        WRITE(var_shortname,  '(A1,A)')  'A', TRIM(radioact_shortnames(jt))
        WRITE(var_description,'(A,A)')   'Averaged air concentration of ', &
          &                              TRIM(radioact_shortnames(jt))
        cf_desc    = t_cf_var(TRIM(var_shortname), 'Bq m-3', TRIM(var_description), datatype_flt)
        grib2_desc = grib2_var(0, 18, 10, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%radioact(jt)%avg,                 &
          &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class,    &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                      &
          &          isteptype=TSTEP_AVG,                                                         &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB ), &
          &          vert_interp=create_vert_interp_metadata(vert_intp_type =                     &
          &                                                             vintp_types("P","Z","I"), &
          &                                                  vert_intp_method=VINTP_METHOD_LIN ) )
        ! Maximum air concentration [Bq m-3] of radionuclides in given time interval
        WRITE(var_shortname,  '(A,A4)')  TRIM(radioact_shortnames(jt)), '_MAX'
        WRITE(var_description,'(5A)')    'Maximum air concentration of ', TRIM(radioact_shortnames(jt)), &
          &                              ' since end of previous full ', TRIM(radioact_maxtint(3:)),     &
          &                              ' since model start'
        cf_desc    = t_cf_var(TRIM(var_shortname), 'Bq m-3', TRIM(var_description), datatype_flt)
        grib2_desc = grib2_var(0, 18, 10, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
        IF (lexist_diagxml) &
          &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
        IF (l_diag) THEN
          art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
        ELSE
          art_groups = groups("ART_DIAGNOSTICS")
        ENDIF
        CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%radioact(jt)%maxtint,             &
          &          GRID_UNSTRUCTURED_CELL, ZA_HYBRID,cf_desc,grib2_desc,var_class=var_class,    &
          &          ldims=shape3d_c, lrestart=.FALSE., in_group=art_groups,                      &
          &          isteptype=TSTEP_MAX, initval=0._wp, resetval=0._wp,                          &
          &          action_list=actions(new_action(ACTION_RESET, TRIM(radioact_maxtint))),       &
          &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,  &
          &                                                fallback_type=HINTP_TYPE_LONLAT_NNB ), &
          &          vert_interp=create_vert_interp_metadata(vert_intp_type =                     &
          &                                                             vintp_types("P","Z","I"), &
          &                                                  vert_intp_method=VINTP_METHOD_LIN ) )

        ! Maximal maximum air concentration of radionuclides in time interval between given pressure levels
        ALLOCATE (art_diag%radioact(jt)%maxtint_layer(npreslay))
        DO jp=1, npreslay
          art_diag%radioact(jt)%maxtint_layer(jp)%pres_bot = pres_bot(jp)
          art_diag%radioact(jt)%maxtint_layer(jp)%pres_top = pres_top(jp)
          WRITE(var_shortname,  '(A,A5,A3,A1,A3)')  TRIM(radioact_shortnames(jt)), '_MAX_', fl_bot(jp), '_', fl_top(jp)
          WRITE(var_description,'(A,A3,A4,A3,A1)')  &
            &  'Maximal maximum air concentration of radionuclides in time interval between given pressure levels   &
            &  (Flight level ', fl_bot(jp), ' to ', fl_top(jp), ')'
          cf_desc    = t_cf_var(TRIM(var_shortname), 'Bq m-3', TRIM(var_description),         &
            &                   datatype_flt)
          grib2_desc = grib2_var(0, 18, 15, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
          IF (lexist_diagxml) &
            &  CALL art_get_GRIB2_diag(var_shortname, grib2_desc, var_class, l_diag)
          IF (l_diag) THEN
            art_groups = groups("ART_DIAGNOSTICS", "ART_ROUTINE_DIAG")
          ELSE
            art_groups = groups("ART_DIAGNOSTICS")
          ENDIF
          CALL add_var(p_diag_list, TRIM(var_shortname), art_diag%radioact(jt)%maxtint_layer(jp)%maximum, &
            &          GRID_UNSTRUCTURED_CELL, ZA_PRES_FL_BOT_TOP, cf_desc, grib2_desc,                   &
            &          var_class=var_class, ldims=shape2d_c, lrestart=.FALSE., in_group=art_groups,       &
            &          isteptype=TSTEP_MAX, initval=0._wp, resetval=0._wp,                                &
            &          action_list=actions(new_action(ACTION_RESET, TRIM(radioact_maxtint))),             &
            &          hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR,        &
            &                                                fallback_type=HINTP_TYPE_LONLAT_NNB))
        END DO
      END DO
    ENDIF
  ENDIF

  ! ----------------------------------
  ! --- 8.0 Chemistry
  ! ----------------------------------

  IF (art_config(jg)%lart_chem) THEN
    ! ----------------------------------
    ! --- 8.1 Photolysis
    ! ----------------------------------

    IF (p_art_data(jg)%chem%param%OH_chem_meta%is_init   &
      &  .OR. art_config(jg)%lart_mecca) THEN
      cf_desc    = t_cf_var('photolysis_fields', '1/s', 'photolysis rates', datatype_flt)
      grib2_desc = grib2_var(0, 254, 218, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'photo', art_chem%photo%rate,                      &
        &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
        &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
        &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
        &             ldims=shape4d_photo,                                          &
        &             lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )

      ALLOCATE(art_diag%art_photolysis_ptr(nphot))
      DO jt =1,nphot
        WRITE(ctrc,'(I2)')jt
        CALL add_ref( p_diag_list, 'photo', 'photo'//TRIM(ADJUSTL(ctrc)),       &
          &             art_diag%art_photolysis_ptr(jt)%p_3d,                   &
          &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                      &
          &             cf_desc, grib2_desc, ref_idx=jt, ldims=shape3d_c, lrestart=.FALSE., &
          &             vert_interp=create_vert_interp_metadata(                &
          &             vert_intp_type=vintp_types("P","Z","I"),                &
          &             vert_intp_method=VINTP_METHOD_LIN))
      ENDDO
    ENDIF

    ! ----------------------------------
    ! --- 8.2 O3 Column (full chemistry)
    ! ----------------------------------

    IF (art_config(jg)%lart_mecca) THEN
#ifdef __ART_GPL
       cf_desc    = t_cf_var('art_full_chemistry_o3_col', 'DU', 'Ozone column', datatype_flt)
       grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( p_diag_list, 'art_full_chemistry_o3_col',                       &
         &           art_chem%mecicon%utils%o3_column,                               &
         &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
         &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
         &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
         &           ldims=shape3d_c, lrestart=.FALSE.)
#endif
    ENDIF

    ! ----------------------------------
    ! --- 8.3 PSC structures
    ! ----------------------------------

    IF (art_config(jg)%lart_psc) THEN
#ifdef __ART_GPL
      cf_desc    = t_cf_var('sts_liqsur', 'cm2 cm-3', 'liquid area density of STS', datatype_flt)
      grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'sts_liqsur', PSC%liqsur,                          &
        &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
        &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
        &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
        &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))


      cf_desc    = t_cf_var('cgaml', '-', 'STS uptake coefficient of the reaction', datatype_flt)
      grib2_desc = grib2_var(0, 254, 218, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'cgaml', PSC%cgaml,                                    &
        &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,           &
        &             vert_interp=create_vert_interp_metadata(vert_intp_type =          &
        &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),      &
        &             ldims=shape4d_ihs,                                                &
        &             lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )

      ALLOCATE(art_diag%art_cgaml_ptr(ihs_MAX))
      DO jt =1,ihs_MAX
        WRITE(ctrc,'(I2.2)')jt
        CALL add_ref( p_diag_list, 'cgaml', 'cgaml'//TRIM(ADJUSTL(ctrc)),       &
          &             art_diag%art_cgaml_ptr(jt)%p_3d,                        &
          &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                      &
          &             cf_desc, grib2_desc,ref_idx=jt, ldims=shape3d_c, lrestart=.FALSE., &
          &             vert_interp=create_vert_interp_metadata(                &
          &             vert_intp_type=vintp_types("P","Z","I"),                &
          &             vert_intp_method=VINTP_METHOD_LIN),                     &
          &             in_group=groups("ART_DIAGNOSTICS"))
     ENDDO

     cf_desc    = t_cf_var('dens_ice', 'm-3', 'number density of ice particles', datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'dens_ice', PSC%dens_ice,                          &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('radius_ice', 'm', 'radius of ice particles', datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'radius_ice', PSC%radius_ice,                      &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('radius_STS', 'm', 'radius of STS particles', datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'radius_STS', PSC%radius_STS,                      &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('dens_NAT', 'm-3', 'number density of NAT particles', datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'dens_NAT', PSC%dens_NAT,                          &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape4d_NSB, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('radius_NAT', 'm', 'radius of NAT particles', datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'radius_NAT', PSC%radius_NAT,                      &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape4d_NSB, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('HNO3_Nconc_s', 'cm-3', 'number concentration of HNO3 in NAT',  &
                    &      datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'HNO3_Nconc_s', PSC%HNO3_Nconc_s,                  &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('HNO3_Nconc_l', 'cm-3', 'number concentration of HNO3 in STS',  &
                     &     datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'HNO3_Nconc_l', PSC%HNO3_Nconc_l,                  &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))

     cf_desc    = t_cf_var('ice_vmr_Marti', 'mol mol-1',                           &
                &   'volume mixing ratio of solid water by Marti and Mauersberger', datatype_flt)
     grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'ice_vmr_Marti', PSC%ice_vmr_Marti,                &
       &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,        &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &           ldims=shape3d_c, lrestart=.FALSE.,in_group=groups("ART_DIAGNOSTICS"))


     cf_desc = t_cf_var('NAT_sedi_rel_difference', '-',                                          &
       &       'relative difference of NAT mass bef and aft sedi (aft - bef) * 2 / (aft + bef)', &
       &       datatype_flt)
     grib2_desc = grib2_var(0, 254, 218, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'NAT_sedi_rel_difference', PSC%NAT_sedi_rel_diff,  &
       &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &             ldims=shape3d_NSB,                                            &
       &             lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )
     ALLOCATE(art_diag%art_NAT_rel_diff_ptr(PSC%NSB))
     DO jt =1,PSC%NSB
       WRITE(ctrc,'(I2.2)')jt

       CALL add_ref( p_diag_list, 'NAT_sedi_rel_difference',                   &
         &             'NAT_sedi_rel_diff'//TRIM(ADJUSTL(ctrc)),               &
         &             art_diag%art_NAT_rel_diff_ptr(jt)%p_2d,                 &
         &             GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                     &
         &             cf_desc, grib2_desc, ref_idx=jt, ldims=shape2d_c,       &
         &             lrestart=.FALSE.,                                       &
         &             in_group=groups("ART_DIAGNOSTICS"))
     ENDDO


     cf_desc  = t_cf_var('NAT_sedi_vel', 'm s-1',                                   &
              &   'sedimentation velocity of NAT particles',datatype_flt)
     grib2_desc = grib2_var(0, 254, 218, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
     CALL add_var( p_diag_list, 'NAT_sedi_vel', PSC%v_sed_NAT_out,                 &
       &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID, cf_desc, grib2_desc,       &
       &             vert_interp=create_vert_interp_metadata(vert_intp_type =      &
       &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),  &
       &             ldims=shape4d_NSB,                                            &
       &             lcontainer=.TRUE., lrestart=.FALSE., loutput=.TRUE. )

     ALLOCATE(art_diag%art_NAT_sedi_vel_ptr(PSC%NSB))
     DO jt =1,PSC%NSB
       WRITE(ctrc,'(I2.2)')jt
       CALL add_ref( p_diag_list, 'NAT_sedi_vel',                              &
         &             'NAT_sedi_vel'//TRIM(ADJUSTL(ctrc)),                    &
         &             art_diag%art_NAT_sedi_vel_ptr(jt)%p_3d,                 &
         &             GRID_UNSTRUCTURED_CELL, ZA_HYBRID,                      &
         &             cf_desc, grib2_desc, ref_idx=jt,ldims=shape3d_c, lrestart=.FALSE., &
         &             vert_interp=create_vert_interp_metadata(                &
         &             vert_intp_type=vintp_types("P","Z","I"),                &
         &             vert_intp_method=VINTP_METHOD_LIN),in_group=groups("ART_DIAGNOSTICS"))
     ENDDO
#endif
    END IF !lart_psc

    ! ----------------------------------
    ! --- 8.4 SO2 Column
    ! ----------------------------------
    IF (art_config(jg)%lart_chemtracer) THEN

       cf_desc    = t_cf_var('art_so2_col', 'DU', 'SO2 column', datatype_flt)
       grib2_desc = grib2_var(0, 254, 253, DATATYPE_PACK16, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( p_diag_list, 'art_so2_col', art_chem%so2_column,                             &
         &           GRID_UNSTRUCTURED_CELL, ZA_HYBRID,  cf_desc, grib2_desc,                     &
         &             vert_interp=create_vert_interp_metadata(vert_intp_type =                   &
         &             vintp_types("P","Z","I"),vert_intp_method=VINTP_METHOD_LIN),               &
         &           ldims=shape3d_c, lrestart=.FALSE., in_group=groups("ART_DIAGNOSTICS"))

    ENDIF 

  END IF !lart_chem

  ! ----------------------------------
  ! --- 9.0 Clean up
  ! ----------------------------------

  IF (lexist_diagxml) THEN
    DO idx_diag_xml = 1, ndiag_xml
      ! Set meta data storage free again
      CALL meta_storage(idx_diag_xml)%destruct
    END DO
    DEALLOCATE(meta_storage)
  END IF

  ! --------------------------------
  ! --- 10.0 Write FPLUME Output
  ! --------------------------------
  
  IF (art_config(jg)%iart_fplume/=0) THEN
    CALL add_var(p_diag_list,'plume_height', fplume%plume_H,grid_unstructured_cell,za_surface,&
        & t_cf_var('plume_height', 'm', 'plume height',datatype_flt),                         &
        & grib2_var(255, 255, 255, datatype_pack16,GRID_UNSTRUCTURED,grid_cell),              &
        & ldims=shape2d_c, lrestart=.FALSE.,in_group=groups("ART_FPLUME"))

    CALL add_var(p_diag_list,'plume_MFR', fplume%plume_MER,grid_unstructured_cell,za_surface,&
        & t_cf_var('plume_MFR', 'kg/s', 'plume MFR',datatype_flt),                           &
        & grib2_var(255, 255, 255, datatype_pack16,GRID_UNSTRUCTURED,grid_cell),             &
        & ldims=shape2d_c, lrestart=.FALSE.,in_group=groups("ART_FPLUME"))

    CALL add_var(p_diag_list,'MER_transport',fplume%MER_transport,grid_unstructured_cell,za_surface,&
        & t_cf_var('MER_transport', 'kg/s', 'Amount of very fine ash for transport',datatype_flt),  &
        & grib2_var(255, 255, 255, datatype_pack16,GRID_UNSTRUCTURED,grid_cell),&
        & ldims=shape2d_c, lrestart=.FALSE.,in_group=groups("ART_FPLUME"))

  ENDIF

END SUBROUTINE art_create_diagnostics
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_get_GRIB2_diag(shortname, grib2, var_class, l_diag,           &
  &                           itrac, p_prog_list, shortnames_active_tracer,  &
  &                           parameterNumber_list, ntrac_active_tracer)

  CHARACTER(LEN=*), INTENT(in)             :: &
    &  shortname                                !< shortName of the variable
  TYPE(t_grib2_var), INTENT(inout)         :: &
    &  grib2                                    !< GRIB2 variable description
  INTEGER, INTENT(out)                     :: &
    &  var_class                                !< Variable class, used to set the correct PDT
  LOGICAL, INTENT(out)                     :: &
    &  l_diag                                   !< logical indicating if definition was found in XML file
  INTEGER,  OPTIONAL, INTENT(in)           :: &
    & itrac                                     !< index of tracer in tracer loop of current diagnostic
  TYPE(t_var_list_ptr), OPTIONAL, INTENT(in)   :: &
    & p_prog_list                               !< current prognostic list, to read names of active tracers
  INTEGER,  OPTIONAL, INTENT(out)          :: &
    & ntrac_active_tracer                       !< number of active tracers, for which a certain 
                                                !  diagnostics is requested
  CHARACTER(len=30), ALLOCATABLE, OPTIONAL, INTENT(out) :: &
    & shortnames_active_tracer(:)               !< shortnames of active tracers, for which a 
                                                !  certain diagnostics shall be compted
  INTEGER, ALLOCATABLE, OPTIONAL, INTENT(out) :: &
    & parameterNumber_list(:)                   !< List of parameterNumbers for each tracer of a 
                                                !  certain diagnostics
  ! Local
  CHARACTER(:), ALLOCATABLE                :: &
    &  name_diag_xml                            !< Name of diagnostic variable in XML file
  INTEGER                                  :: &
    &  igrib2_val,                            & !< Integer value of GRIB2 key read from XML file
    &  jkey, jkey_end,                        & !< loop index and loop end for GRIB2 keys
    &  idx_diag_xml,                          & !< loop index for diagn. vars
    &  ierror
  CHARACTER(:), ALLOCATABLE                :: &
    &  cdiag_val                                !< character array with tracer names for which a 
                                                !  diagnostics shall be performed, read from XML 
                                                !  file
  CHARACTER(LEN=MAX_CHAR_LENGTH)           :: &
    &  character_list                           !< ... to save cdiag_val to character with 
                                                !      explicit lenght
  CHARACTER(LEN=IART_VARNAMELEN)           :: &
    &  cgrib2_keys(22)                          !< Array with optional GRIB2 keys
  INTEGER                                  :: &
    &  pos(100), ilength(100)                   !< DUMMY dimension (allows for 100 tracers to be 
                                                !  specified)
  INTEGER                                  :: &
    &  i,j                                      !< Loop index, counter
  INTEGER                                  :: &
    & nnum_diag                                 !< number of parameterNumbers specified for all 
                                                !  tracers of a diagnostics block
  CHARACTER(LEN=IART_VARNAMELEN), ALLOCATABLE :: &
    & shortnames(:)                             !< shortnames of all tracers, for which a certain 
                                                !  diagnostics is requested
  INTEGER, ALLOCATABLE                     :: &
    & ind_active_tracer(:)                      !< index of shortname in _list block, if shortname 
                                                !  is not an active tracer, ind_active_tracer = 0
  INTEGER, ALLOCATABLE                     :: &
    & productDefinitionTemplate_list(:)         !< List of productDefinitionTemplate Numbers for 
                                                !  each tracer of a certain diagnostics
  CHARACTER(LEN=IART_VARNAMELEN)           :: &
    & tracer_name                               !< tracer_name of one tracer, for which a diagnostics shall be computed
  CHARACTER(LEN=IART_VARNAMELEN)           :: &
    & current_tracer_name                       !< name of current tracer as CharArray
  INTEGER                                  :: &
    & ntrac,                                  & !< number of tracers, for which a certain diagnostics is requested
    & iv                                        !< loop index


  var_class = CLASS_DEFAULT

  l_diag = .FALSE.

  ! Optional aerosol GRIB2 metadata
  cgrib2_keys = (/ 'productDefinitionTemplate                         ',  &
    &              'constituentType                                   ',  &
    &              'aerosolType                                       ',  &
    &              'typeOfFirstFixedSurface                           ',  &
    &              'typeOfSecondFixedSurface                          ',  &
    &              'scaledValueOfFirstFixedSurface                    ',  &
    &              'scaledValueOfSecondFixedSurface                   ',  &
    &              'scaleFactorOfFirstFixedSurface                    ',  &
    &              'scaleFactorOfSecondFixedSurface                   ',  &
    &              'typeOfSizeInterval                                ',  &
    &              'typeOfWavelengthInterval                          ',  &
    &              'scaledValueOfFirstWavelength                      ',  &
    &              'scaleFactorOfFirstWavelength                      ',  &
    &              'localInformationNumber                            ',  &
    &              'typeOfDistributionFunction                        ',  &
    &              'numberOfModeOfDistribution                        ',  &
    &              'modeNumber                                        ',  &
    &              'numberOfDistributionFunctionParameters            ',  &
    &              'scaledValueOfDistributionFunctionParameter_1      ',  &
    &              'scaleFactorOfDistributionFunctionParameter_1      ',  &
    &              'scaledValueOfDistributionFunctionParameter_2      ',  &
    &              'scaleFactorOfDistributionFunctionParameter_2      ' /)


  DO idx_diag_xml = 1, ndiag_xml

    !set default number of tracers per diagnostics-block
    IF (PRESENT(ntrac_active_tracer)) ntrac_active_tracer = 0

    CALL key_value_storage_as_string(meta_storage(idx_diag_xml), 'name', name_diag_xml, ierror)
    IF (ierror /= SUCCESS) name_diag_xml = ''

    IF ( TRIM(ADJUSTL(tolower(shortname))) == TRIM(ADJUSTL(tolower(name_diag_xml))) ) THEN

      l_diag = .TRUE.

      ! Set basic GRIB2 meta data
      CALL meta_storage(idx_diag_xml)%get('discipline',igrib2_val,ierror)
      IF (ierror == SUCCESS) grib2%discipline = igrib2_val
      CALL meta_storage(idx_diag_xml)%get('parameterCategory',igrib2_val,ierror)
      IF (ierror == SUCCESS) grib2%category = igrib2_val
      CALL meta_storage(idx_diag_xml)%get('parameterNumber',igrib2_val,ierror)
      IF (ierror == SUCCESS) grib2%number = igrib2_val
      CALL meta_storage(idx_diag_xml)%get('bitsPerValue',igrib2_val,ierror)
      IF (ierror == SUCCESS) THEN
        SELECT CASE(igrib2_val)
        CASE(24)
          grib2%bits = DATATYPE_PACK24
        CASE(32)
          grib2%bits = CDI_DATATYPE_PACK32
        END SELECT
      ENDIF

      ! set variable class according to productDefinitionTemplate
      CALL meta_storage(idx_diag_xml)%get(TRIM(cgrib2_keys(1)),igrib2_val,ierror)
      IF (ierror == SUCCESS) THEN
        SELECT CASE (igrib2_val)
        CASE (40)
          var_class = CLASS_CHEM
        CASE (42)
          var_class = CLASS_CHEM_STAT
        CASE (48)
          var_class = CLASS_CHEM_OPTP
        CASE (57)
          var_class = CLASS_DISTR
        CASE (67)
          var_class = CLASS_DISTR_STAT
        CASE default
          var_class = CLASS_DEFAULT
        END SELECT
      END IF


      !--- read '*_list'-tags from the diagnostic.xml file if available
      ! get optional tracer names for which a diagnostic-block shall be computed
      CALL key_value_storage_as_string(meta_storage(idx_diag_xml), 'tracername_list', &
        &                              cdiag_val, ierror)
      IF (ierror == SUCCESS) THEN
         character_list = TRIM(cdiag_val)
         ! get number of tracer names indicated in xml-block
         CALL split_string(character_list, ntrac, pos, ilength)

         ALLOCATE(shortnames(ntrac))
         ! get tracer shortnames
         DO i=1,ntrac
            IF (i < ntrac) THEN
               shortnames(i) = character_list(pos(i):pos(i)+ilength(i)-1)
            ELSE
               shortnames(i) = character_list(pos(i):pos(i)+ilength(i))
            END IF
         END DO

         ! check if shortname is an actual active tracer name
         ntrac_active_tracer = 0
         ALLOCATE(ind_active_tracer(ntrac)) !if shortname is not an active tracer, 
                                            !  ind_active_tracer = 0
         DO i=1,ntrac
            ind_active_tracer(i) = 0
            WRITE(tracer_name,'(A)')shortnames(i)
            DO iv = 1, p_prog_list%p%nvars
               current_tracer_name = TRIM(p_prog_list%p%vl(iv)%p%info%name(1: &
                    & (LEN(TRIM(p_prog_list%p%vl(iv)%p%info%name))-4))) !e.g. dusta from dusta.TL1
               IF ( current_tracer_name == TRIM(tracer_name) ) THEN
                  ntrac_active_tracer = ntrac_active_tracer + 1
                  ind_active_tracer(i) = i
                  EXIT
               END IF
            END DO
         END DO

         !store only shortnames of active tracers
         ALLOCATE(shortnames_active_tracer(ntrac_active_tracer))
         j = 0
         DO i=1,ntrac
            IF ( ind_active_tracer(i) == i ) THEN
               j = j + 1
               shortnames_active_tracer(j) = shortnames(i)
            END IF
         END DO

         ! get parameterNumber for each tracer in this diagnostic block,
         ! but only if it is an active tracer
         CALL key_value_storage_as_string(meta_storage(idx_diag_xml), 'parameterNumber_list', &
           &                              cdiag_val, ierror)
         IF (ierror == SUCCESS) THEN
            character_list = TRIM(cdiag_val)
            CALL split_string(character_list, nnum_diag, pos, ilength)
            IF (nnum_diag == ntrac) THEN
               ALLOCATE(parameterNumber_list(ntrac_active_tracer))
               j = 0
               DO i=1,nnum_diag
                  IF ( ind_active_tracer(i) == i ) THEN
                     j = j + 1
                     IF (i < nnum_diag) THEN
                        !convert string (e.g. '192,192,192,193,193,193') to number
                        READ(character_list(pos(i):pos(i)+ilength(i)-1), *) parameterNumber_list(j)
                     ELSE
                        READ(character_list(pos(i):pos(i)+ilength(i)), *) parameterNumber_list(j)
                     END IF
                  END IF !shortname is an active tracer
               END DO
            END IF !(nnum_diag == ntrac)
         END IF !(ierror == SUCCESS)'parameterNumber_list'

         ! get productDefinitionTemplate for each tracer in this diagnostic block,
         ! but only if it is an active tracer
         CALL key_value_storage_as_string(meta_storage(idx_diag_xml), 'productDefinitionTemplate_list', &
           &                              cdiag_val, ierror)
         IF (ierror == SUCCESS) THEN
            character_list = TRIM(cdiag_val)
            CALL split_string(character_list, nnum_diag, pos, ilength)
            IF (nnum_diag == ntrac) THEN
               ALLOCATE(productDefinitionTemplate_list(ntrac_active_tracer))
               j = 0
               DO i=1,nnum_diag
                  IF ( ind_active_tracer(i) == i ) THEN
                     j = j + 1
                     IF (i < nnum_diag) THEN
                        !convert string (e.g. '67,67,67,67,67,42') to number
                        READ(character_list(pos(i):pos(i)+ilength(i)-1), *)  &
                          &                                       productDefinitionTemplate_list(j)
                     ELSE
                        READ(character_list(pos(i):pos(i)+ilength(i)), *)    &
                          &                                       productDefinitionTemplate_list(j)
                     END IF
                  END IF !shortname is an active tracer
               END DO

               ! set variable class for current tracer number itrac
               IF ( PRESENT(itrac) ) THEN
                  SELECT CASE (productDefinitionTemplate_list(itrac))
                  CASE (40)
                     var_class = CLASS_CHEM
                  CASE (42)
                     var_class = CLASS_CHEM_STAT
                  CASE (48)
                     var_class = CLASS_CHEM_OPTP
                  CASE (57)
                     var_class = CLASS_DISTR
                  CASE (67)
                     var_class = CLASS_DISTR_STAT
                  CASE default
                     var_class = CLASS_DEFAULT
                  END SELECT
               END IF
            ELSE
               WRITE(message_text,'(A)') 'ART:WARNING: number of entries in '// &
                 &                       'productDefinitionTemplate_list /= tracername_list! '
               CALL message (TRIM(routine)//':ART:', message_text)
               WRITE(message_text,'(A,A)') 'An error will occure when writing diagnostic: ', &
                 &                         shortname
               CALL message (TRIM(routine)//':ART:', message_text)
            END IF !(nnum_diag == ntrac)
         END IF !(ierror == SUCCESS)'productDefinitionTemplate_list'
      END IF !(ierror == SUCCESS)'tracername_list'
      !---

      IF ( ANY ( (/CLASS_DISTR_STAT, CLASS_DISTR/) == var_class ) ) THEN
         ! set also: typeOfDistributionFunction, numberOfModeOfDistribution,
         !           modeNumber, numberOfDistributionFunctionParameters
         jkey_end = 18
      ELSE
         jkey_end = 14
      END IF

      ! set (scalar) integer valued additional GRIB2 keys
      DO jkey = 2, jkey_end ! 14 or 18
        CALL meta_storage(idx_diag_xml)%get(TRIM(cgrib2_keys(jkey)),igrib2_val,ierror)
        IF (ierror == SUCCESS) THEN
          grib2 = grib2 + t_grib2_int_key(TRIM(cgrib2_keys(jkey)), igrib2_val)
        END IF
      END DO

      RETURN

    END IF

  END DO

END SUBROUTINE art_get_GRIB2_diag



SUBROUTINE get_tracer_index_and_grib_meta(p_prog_list, tracer_name, jsp_idx, grib2)
  !<
  ! SUBROUTINE get_tracer_index_and_grib_meta
  ! get index of a given tracer in the tracer
  ! container as well as grib2 metadata
  ! given the tracer name, e.g. 'dusta' (as
  ! defined in the tracer.xml file)
  ! Based on: -
  ! Part of Module: mo_art_diag_state
  ! Author: Andrea Steiner, DWD
  ! Initial Release: 2017-11-24
  ! Modifications:
  ! 2017-11-24: Andrea Steiner, DWD
  ! -
  !>

  TYPE(t_var_list_ptr),INTENT(in)  :: p_prog_list         !< current prognostic list
  CHARACTER(LEN=*),    INTENT(in)  :: tracer_name         !< name of tracer, for which to return the index in tracer container
  INTEGER, POINTER,    INTENT(out) :: jsp_idx             !< index of tracer in tracer container
  TYPE(t_grib2_var),   INTENT(out) :: grib2               !< grib2 metadata

  !Local variables
  CHARACTER(LEN=50)                :: current_tracer_name !< name of current tracer
  INTEGER                          :: iv                  !< loop index                 


  ! Loop through the p_prog_list
  DO iv = 1, p_prog_list%p%nvars

     ! Plain variable name (i.e. without TIMELEVEL_SUFFIX)
     current_tracer_name = tolower(get_var_name(p_prog_list%p%vl(iv)%p%info)) !e.g. dusta from dusta.TL1

     ! If the requestet tracer name is reached get index in tracer container 
     ! as well as grib2 metadata
     IF ( TRIM(current_tracer_name) == TRIM(tolower(tracer_name)) ) THEN
        jsp_idx => p_prog_list%p%vl(iv)%p%info%ncontained
        grib2 = p_prog_list%p%vl(iv)%p%info%grib2
        EXIT
     ELSE
        NULLIFY(jsp_idx)
        grib2 = grib2_var(-1, -1, -1, -1, -1, -1)
     END IF

  ENDDO !Loop through p_prog_list

END SUBROUTINE get_tracer_index_and_grib_meta


SUBROUTINE get_and_set_int_grib2_key(grib2key, grib2_in, grib2_inout)
  !<
  ! SUBROUTINE get_and_set_int_grib2_key
  ! search grib2 derived data type 'grib2_in' for additional grib key
  ! 'grib2key' and transfers grib2 key to 'grib2_inout'.
  !
  ! Attention: the value of grib2 key must be an integer
  !
  ! Example:
  ! get ModeNumber from grib2_tracer and write it into grib2_desc of diagnostic variable:
  ! get_and_set_int_grib2_key('modeNumber', grib2_tracer, grib2_desc)
  !
  ! Based on: -
  ! Part of Module: mo_art_diag_state
  ! Author: Andrea Steiner, DWD
  ! Initial Release: 2017-12-19
  ! Modifications:
  ! 2017-11-24: Andrea Steiner, DWD
  ! -
  !>

  CHARACTER(LEN=*), INTENT(in)     :: grib2key    !< grib2key e.g. 'modeNumber'
  TYPE(t_grib2_var), INTENT(in)    :: grib2_in    !< GRIB2 variable description e.g. of tracer
  TYPE(t_grib2_var), INTENT(inout) :: grib2_inout !< GRIB2 variable description e.g. of diagnostic 
                                                  !  variable
  !Local variables
  INTEGER                          :: i                !< Loop index
  INTEGER                          :: grib2val         !< Value of grib2key in grib2_in

  ! Set general grib2 settings if requested
  IF ( TRIM(grib2key) == 'discipline' ) THEN
     grib2_inout%discipline = grib2_in%discipline
  ELSE IF ( TRIM(grib2key) == 'parameterCategory' ) THEN
     grib2_inout%category = grib2_in%category
  ELSE IF ( TRIM(grib2key) == 'parameterNumber' ) THEN
     grib2_inout%number = grib2_in%number
  ELSE

     ! Loop over all additional grib2 keys in grib2_in
     DO i=1,grib2_in%additional_keys%nint_keys

        ! Search grib2_in for grib2key, get value
        IF ( TRIM(grib2_in%additional_keys%int_key(i)%key) == TRIM(grib2key) ) THEN

           grib2val = grib2_in%additional_keys%int_key(i)%val
            ! append grib2 info
           grib2_inout = grib2_inout + t_grib2_int_key(TRIM(grib2key), grib2val)

           !WRITE(message_text,'(A,A,A,I6,A)') 'ART: INFO: additional grib-key ', &
           !  &   TRIM(grib2key),' = ',grib2val,' was transferred from tracer to diagnostic.'
           !CALL message (TRIM(routine)//':get_and_set_int_grib2_key', message_text)

        END IF
     END DO
  END IF

END SUBROUTINE get_and_set_int_grib2_key

!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_diag_state

