!
! mo_art_tracer_def_wrapper
! This module provides a wrapper between the reading of the XML lists of tracers and
! the tracer definitions via add_ref
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

MODULE mo_art_tracer_def_wrapper
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mo_impl_constants,                ONLY: VINTP_METHOD_LIN, SUCCESS
  USE mo_physical_constants,            ONLY: amd

  USE mo_var_metadata_types,            ONLY: t_var_metadata, t_var_metadata_dynamic,  &
                                          &   POST_OP_NONE,                            &
                                          &   POST_OP_SCALE,                           &
                                          &   t_vert_interp_meta,                      &
                                          &   CLASS_DEFAULT, CLASS_CHEM, CLASS_DISTR
  USE mo_var_metadata,                  ONLY: create_vert_interp_metadata,             &
                                          &   vintp_types, post_op, get_timelevel_string
  USE mo_var_groups,                    ONLY: groups, MAX_GROUPS, var_groups_dyn
  USE mo_fortran_tools,                 ONLY: t_ptr_2d3d
  USE mo_advection_config,              ONLY: t_advection_config
  USE mo_cf_convention,                 ONLY: t_cf_var
  USE mo_grib2,                         ONLY: grib2_var, t_grib2_var, t_grib2_int_key, &
    &                                         OPERATOR(+)
  USE mo_cdi,                           ONLY: DATATYPE_PACK16,                         &
                                          &   DATATYPE_PACK24, CDI_DATATYPE_PACK32,    &
                                          &   GRID_UNSTRUCTURED
  USE mo_cdi_constants,                 ONLY: GRID_CELL
  USE mo_art_config,                    ONLY: t_art_config
  USE mo_var_list,                      ONLY: add_ref, find_list_element, t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_tracer_metadata,               ONLY: create_tracer_metadata_aero

  USE mo_tracer_metadata_types,         ONLY: t_tracer_meta, t_aero_meta, t_chem_meta
  USE mo_key_value_store,               ONLY: t_key_value_store
  USE mo_util_string,                   ONLY: tolower
  USE mo_grid_config,                   ONLY: l_limited_area
  USE mo_limarea_config,                ONLY: latbc_config

  USE mo_art_general_interface,         ONLY: getNetcdfPrecision
! ART
  USE mo_art_impl_constants,            ONLY: IART_AERO_TR, IART_CHEM_TR,    &
                                          &   UNDEF_REAL_ART, UNDEF_INT_ART, &
                                          &   IART_VARNAMELEN
  USE mo_art_read_xml,                  ONLY: t_xml_file, art_read_elements_xml

  USE mo_art_chem_types_param,          ONLY: t_chem_meta_lt, t_chem_meta_OH,       &
                                          &   t_chem_meta_linoz,                    &
                                          &   t_chem_meta_simnoy, t_chem_meta_cold
  USE mo_art_chem_types,                ONLY: t_chem_meta_param,   &
                                          &   t_chem_meta_passive, &
                                          &   t_chem_meta_mecca
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string

  IMPLICIT NONE
  
  PRIVATE
  
  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_tracer_def_wrapper'
  
  PUBLIC  :: art_tracer_def_wrapper
  
  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_tracer_def_wrapper(IART_TRACER_TYPE, defcase, art_config, advconf, &
  &                               this_list, vname_prefix, ptr_arr, meta_storage, &
  &                               xmlfile, xpath,                                 &
  &                               tracer_idx, tracer_name_in, tracer_mode,        &
  &                               tracer_sub, modeNumber_in, timelev, ldims)
!<
! SUBROUTINE art_tracer_def_wrapper
! The subroutine art_tracer_def_wrapper gets information of the current
! tracer that shall be created. The subroutine add_ref is called in order 
! to create the tracer. Arguments are set automatically depending on the 
! defcase (whether a prognostic or a tendency field is to be created).
! The index of the current field in the complete tracer vector is returned.
! Part of Module: mo_art_tracer_def_wrapper
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-02
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  INTEGER, INTENT(in)                      :: &
    &  IART_TRACER_TYPE                         !< Tracer type (IART_AERO_TR, IART_CHEM_TR, 
                                                !               IART_PASS_TR)
  CHARACTER(len=*), INTENT(in)             :: &
    &  defcase,                               & !< cases:  prog, conv, turb
    &  vname_prefix                             !< prefix for variable names (usually none)
  CHARACTER(len=100), INTENT(in)           :: &
    &  tracer_name_in
  CHARACTER(len=IART_VARNAMELEN),OPTIONAL  :: &
    &  tracer_mode                              !< Name of mode the aerosol tracer is contained in
  CHARACTER(len=*), INTENT(in), OPTIONAL   :: & ! (currently) only used for IART_AERO_TR
    &  tracer_sub
  TYPE(t_art_config), INTENT(inout)        :: &
    &  art_config                               !< ART configuration state
  TYPE(t_advection_config), INTENT(inout)  :: &
    &  advconf                                  !< Advection configuration state
  TYPE(t_var_list_ptr), INTENT(inout)          :: &
    &  this_list                                !< current list (prog, diag, tend)
  TYPE(t_ptr_2d3d),INTENT(inout)           :: &
    &  ptr_arr(:)                               !< Pointer to mixing ratio field
  TYPE(t_key_value_store), INTENT(in)      :: &
    &  meta_storage                             !< Meta data storage container
  TYPE(t_xml_file),INTENT(in)              :: &
    &  xmlfile                                  !< XML file with tracer metadata
  CHARACTER(LEN=*),INTENT(in)              :: &
    &  xpath                                    !< X-Path to current element (tracer) in XML file
  INTEGER, INTENT(inout)                   :: &
    &  tracer_idx                               !< index of tracer in container
  CHARACTER(len=*), INTENT(in), OPTIONAL   :: & ! (currently) only used for IART_AERO_TR
    &  modeNumber_in                            !  modeNumber for GRIB2
  INTEGER,INTENT(in),OPTIONAL              :: &
    &  timelev,                               & !< necessary for prognostic tracer list. 
                                                !  Not for physical tendencies 
    &  ldims(3)                                 !< local dimensions, for checking
! Local variables
  TYPE(t_var),POINTER             :: &
    &  target_element                           !< Pointer to element in this_list
  TYPE(t_var_metadata),POINTER             :: &
    &  target_info                              !< Pointer to metadata of target_element
  TYPE(t_var_metadata_dynamic),POINTER     :: &
    &  target_info_dyn                          !< Pointer to dynamic metadata of target_element 
                                                !  (i.e. tracer meta)
  TYPE(t_cf_var)                           :: &
    &  cf                                       !< NETCDF variable metadata
  TYPE(t_grib2_var)                        :: &
    &  grib2                                    !< grib2 variable metadata
  TYPE(t_vert_interp_meta)                 :: &
    &  vert_interp                              !< Vertical interpolation metadata
  CLASS(t_tracer_meta), ALLOCATABLE        :: &
    &  tracer_meta                              !< Tracer metadata
  REAL(wp)                                 :: &
    &  mol_weight,                            & !< Molar mass [g mol-1]
    &  rho,                                   & !< Density [kg m-3]
    &  sol                                      !< Solubility
  ! for post operations
  INTEGER                                  :: &
    &  post_op_id                               !< identifier for the post operation
  REAL(wp)                                 :: &
    &  scale_factor                             !< scaling factor for post operation
  INTEGER                                  :: &
    &  idx,                                   & !< Index of tracer in container
    &  ierror,                                & !< Error return code
    &  iconv, iturb,                          & !< integer value of lturb/lconv 0 = false, 1 = true
    &  ivadv_tracer,                          & !< Vertical advection scheme
    &  ihadv_tracer,                          & !< Horizontal advection scheme
    &  ised_tracer,                           & !< Sedimentation scheme
    &  iwash_tracer,                          & !< Washout scheme
    &  moment,                                & !< Moment of aerosol distribution (e.g. 0=number, 
                                                !  3=proportional to mass)
    &  var_class,                             & !< Variable class, used to set the correct PDT
    &  igrib2_val,                            & !< Integer value of GRIB2 key read from XML file
    &  jkey                                     !< loop index for GRIB2 keys

  CHARACTER(:), ALLOCATABLE                :: &
    &  target_name,                           & !< Kind of field (prog/tend,...)
    &  tracer_name,                           & !< Name of tracer including pre- and suffix
    &  units,                                 & !< Unit of tracer (derived from meta_storage)
    &  c_initc,                               & !< character value of initc tag ("file", ...)
    &  c_latbc,                               & !< character value of latbc tag ("file", ...)
    &  c_meteogram,                           & !< character value of meteogram tag ("file",...)
    &  transport_template,                    & !< Name of transport template to be applied 
                                                !  on tracer
    &  c_solve

  CHARACTER(len=IART_VARNAMELEN)           :: &
    &  cgrib2_keys(11)                          !< Array with optional GRIB2 keys
  CHARACTER(len=4)                         :: &
    &  suffix                                   !< Suffix of the tracer name
  LOGICAL                                  :: &
    &  in_group(MAX_GROUPS),                  & !< Output groups
    &  lrestart,                              & !< Variable needed for restart?
    &  lis_tracer,                            & !< Is this variable a tracer?
    &  lturb_tracer,                          & !< Participation in turbulence scheme?
    &  lconv_tracer,                          & !< Participation in convection scheme?
    &  ldep_tracer,                           & !< Participation in dry deposition?
    &  linsol,                                & !< TRUE for insoluble tracer
    &  lfeedback,                             & !< Child -> Parent feedback?
    &  l_initc,                               & !< logical =.TRUE., if initc tag = "file"
    &  l_latbc,                               & !< logical =.TRUE., if latbc tag = "file"
    &  l_modeNumber,                          & !< logical =.TRUE., if modeNumber tag is given
    &  l_meteogram                              !< logical =.TRUE., if meteogram tag = "true"

  INTEGER :: lfeedback_in
  INTEGER :: idx_group_mtgrm


  lconv_tracer = .FALSE.
  lturb_tracer = .FALSE.
!  Definitions
  ierror      = SUCCESS
  in_group(:) = .FALSE.
  
  mol_weight = UNDEF_REAL_ART
  ivadv_tracer = 0
  ihadv_tracer = 0

  ! post operation: start with empty values
  post_op_id   = POST_OP_NONE
  scale_factor = 1._wp

  var_class  = CLASS_DEFAULT

! Always required metadata
  CALL key_value_storage_as_string(meta_storage,'unit',units,ierror)
    IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
                             &  'Required metadata unit not present: '//TRIM(tracer_name_in))

! Defcase dependent settings
  SELECT CASE(TRIM(defcase))
    CASE('prog')
      target_name    = TRIM(vname_prefix)//'tracer'//get_timelevel_string(timelev)
      target_element => find_list_element (this_list, target_name)
      target_info    => target_element%info
      idx            =  target_info%ncontained+1  ! index in 4D tracer container
      WRITE (message_text,'(a,i3,a,a)') 'ART: Tracer index ',idx,  &
                &                       ' assigned to ',TRIM(tracer_name_in)

      CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)

      IF (PRESENT(timelev)) THEN               !< For prognostic tracers, multiple timelevels exist
        WRITE(suffix,'(".TL",i1)') timelev     !< These are distinguished by the suffix of the name
      ELSE
        CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
          &         'timelev is not present for defcase prog')
      ENDIF
      tracer_name  = TRIM(vname_prefix)//TRIM(tracer_name_in)//suffix

      cf%standard_name = TRIM(vname_prefix)//TRIM(tracer_name_in)
      cf%units         = TRIM(units)
      cf%long_name     = TRIM(tracer_name_in)//'_mixing_ratio'
      cf%datatype      = getNetcdfPrecision()

      lrestart   = .TRUE.
      lis_tracer = .TRUE.

      ! Note that sedimentation, deposition and washout are currently handled
      ! by the ART mode structure so these flags do not change any result
      CALL meta_storage%get('ised',ised_tracer,ierror)
      IF (ierror /= SUCCESS) ised_tracer  = 1      !< Default value
      CALL meta_storage%get('iwash',iwash_tracer,ierror)
      IF (ierror /= SUCCESS) iwash_tracer = 1      !< Default value
      CALL meta_storage%get('ldep',ldep_tracer,ierror)
      IF (ierror /= SUCCESS) ldep_tracer  = .TRUE. !< Default value

      ! Check if tracers tendencies due to turbulence and convection are considered
      CALL meta_storage%get('iturb',iturb,ierror)
      IF (ierror /= SUCCESS) THEN
        lturb_tracer = .TRUE. !< Default value
        WRITE (message_text,*) 'ART: WARNING: metadata iturb for '//TRIM(tracer_name_in)  &
                  &          //' not found, using default'
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
      ELSE
        IF (iturb == 1) THEN
          lturb_tracer = .TRUE.
        ELSE
          lturb_tracer = .FALSE.
        ENDIF
      ENDIF
      CALL meta_storage%get('iconv',iconv,ierror)
      IF (ierror /= SUCCESS) THEN
        lconv_tracer = .TRUE. !< Default value
        WRITE (message_text,*) 'ART: WARNING: metadata iconv for '//TRIM(tracer_name_in)  &
                   &         //' not found, using default'
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
      ELSE
        IF (iconv == 1) THEN
          lconv_tracer = .TRUE.
        ELSE
          lconv_tracer = .FALSE.
        ENDIF
      ENDIF

      ! Get and set transport template
      CALL key_value_storage_as_string(meta_storage,'transport',transport_template,ierror)
      IF (ierror /= SUCCESS) THEN
        CALL set_transport(TRIM(tracer_name_in), 'default',                      &
          &                advconf%ivadv_tracer(idx), advconf%itype_vlimit(idx), &
          &                advconf%ihadv_tracer(idx), advconf%itype_hlimit(idx))
      ELSE
        CALL set_transport(TRIM(tracer_name_in), TRIM(transport_template),       &
          &                advconf%ivadv_tracer(idx), advconf%itype_vlimit(idx), &
          &                advconf%ihadv_tracer(idx), advconf%itype_hlimit(idx))
        IF (TRIM(tolower(transport_template)) == 'off') THEN
          IF (lconv_tracer) THEN
            WRITE (message_text,*) 'ART: WARNING: transport for '//TRIM(tracer_name_in)  &
                         &       //' is off but iconv is .true.'
            CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
          ENDIF
          IF (lturb_tracer) THEN
            WRITE (message_text,*) 'ART: WARNING: transport for '//TRIM(tracer_name_in)  &
                        &        //' is off but iturb is .true.'
            CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
          ENDIF
        ENDIF
      ENDIF
      ! Overwrite tracer metadata settings (although those are currently not used)
      ivadv_tracer = advconf%ivadv_tracer(idx)
      ihadv_tracer = advconf%ihadv_tracer(idx)

    CASE('conv')
      CALL meta_storage%get('iconv',iconv,ierror)
      IF (ierror /= SUCCESS) THEN
        lconv_tracer = .TRUE. !< Default value
        WRITE (message_text,*) 'ART: WARNING: metadata iconv for '//TRIM(tracer_name_in)  &
                     &       //' not found, using default'
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
      ELSE
        IF (iconv == 1) THEN
          lconv_tracer = .TRUE.
        ELSE
          lconv_tracer = .FALSE.
        ENDIF
      ENDIF

      IF(lconv_tracer) THEN
        target_name    =  'ddt_tracer_pconv'
        target_element => find_list_element (this_list, target_name)
        target_info    => target_element%info
        idx            =  target_info%ncontained+1  ! index in 4D tracer container
        WRITE (message_text,'(a,i3,a,a)') 'ART: Tracer index ',idx,' assigned to ',  &
                    &                     TRIM(tracer_name_in)
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)

        tracer_name  = TRIM(vname_prefix)//TRIM(tracer_name_in)//'_conv'

        cf%standard_name = tracer_name
        cf%units         = TRIM(units)//' s-1'
        cf%long_name     = 'tendency_of_'//TRIM(tracer_name_in)//'_mixing_ratio_due_to_convection'
        cf%datatype      = getNetcdfPrecision()

        lrestart   = .TRUE.
        lis_tracer = .FALSE.

        ised_tracer  = 0
        iwash_tracer = 0
        lturb_tracer = .FALSE.
        ldep_tracer  = .FALSE.
      ENDIF

    CASE('turb')
      CALL meta_storage%get('iturb',iturb,ierror)
      IF (ierror /= SUCCESS) THEN
        lturb_tracer = .TRUE. !< Default value
        WRITE (message_text,*) 'ART: WARNING: metadata iturb for '//TRIM(tracer_name_in)  &
                    &        //' not found, using default'
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
      ELSE
        IF (iturb == 1) THEN
          lturb_tracer = .TRUE.
        ELSE
          lturb_tracer = .FALSE.
        ENDIF
      ENDIF

      IF(lturb_tracer) THEN
        target_name    =  'ddt_tracer_turb'
        target_element => find_list_element(this_list, target_name)
        target_info    => target_element%info
        idx            =  target_info%ncontained+1  ! index in 4D tracer container
        WRITE (message_text,'(a,i3,a,a)') 'ART: Tracer index ',idx,' assigned to ',  &
                   &                      TRIM(tracer_name_in)
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)

        tracer_name  = TRIM(vname_prefix)//TRIM(tracer_name_in)//'_turb'
      
        cf%standard_name = tracer_name
        cf%units         = TRIM(units)//' s-1'
        cf%long_name     = 'tendency_of_'//TRIM(tracer_name_in)//'_mixing_ratio_due_to_turbulence'
        cf%datatype      = getNetcdfPrecision()

        lrestart   = .FALSE.
        lis_tracer = .FALSE.

        ised_tracer  = 0
        iwash_tracer = 0
        lconv_tracer = .FALSE.
        ldep_tracer  = .FALSE.
      ENDIF
    

    CASE DEFAULT
      CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
        &         'defcase unknown')
  END SELECT
! General settings for all defcases and tracer types

  grib2             = grib2_var(-1, -1, -1, -1, -1, -1) ! initalize GRIB2 structure

  grib2%bits = DATATYPE_PACK16
  CALL meta_storage%get('bitsPerValue',igrib2_val,ierror)
  IF (ierror == SUCCESS) THEN
    SELECT CASE(igrib2_val)
    CASE(24)
      grib2%bits = DATATYPE_PACK24
    CASE(32)
      grib2%bits = CDI_DATATYPE_PACK32
    END SELECT
  ENDIF

  grib2%gridtype    = GRID_UNSTRUCTURED
  grib2%subgridtype = GRID_CELL

  CALL meta_storage%get('discipline',grib2%discipline,ierror)
  IF (ierror /= SUCCESS) grib2%discipline = 0  !< Dummy value in case no grib meta data exists
  CALL meta_storage%get('parameterCategory',grib2%category,ierror)
  IF (ierror /= SUCCESS) grib2%category = 254  !< Dummy value in case no grib meta data exists
  CALL meta_storage%get('parameterNumber',grib2%number,ierror)
  IF (ierror /= SUCCESS) grib2%number = 1      !< Dummy value in case no grib meta data exists

  vert_interp = create_vert_interp_metadata(       &
    & vert_intp_type   = vintp_types('P','Z','I'), &
    & vert_intp_method = VINTP_METHOD_LIN,         &
    & l_loglin         = .FALSE.,                  &
    & l_extrapol       = .FALSE.,                  &
    & l_pd_limit       = .FALSE.,                  &
    & lower_limit      = 0.0_wp)

! Tracer type dependent settings and creation of tracer metadata container
  SELECT CASE (IART_TRACER_TYPE)
    CASE(IART_AERO_TR)
      ALLOCATE(t_aero_meta :: tracer_meta)

      ! Set-up post operation for the aerosol tracers
      CALL meta_storage%get('scale_factor', scale_factor, ierror)
      IF (ierror == SUCCESS) THEN
        post_op_id   = POST_OP_SCALE
      ELSE
        ! safety meassure - probably not necessary
        post_op_id   = POST_OP_NONE
        scale_factor = 1._wp
      END IF

      IF(TRIM(defcase) == 'prog') THEN
        l_initc = .FALSE.
        l_latbc = .FALSE.
        l_meteogram = .FALSE.
        in_group = groups('ART_AEROSOL')
        CALL key_value_storage_as_string(meta_storage,'initc', c_initc, ierror)
        IF (ierror == SUCCESS) THEN
          IF (TRIM(c_initc) == 'file') l_initc = .TRUE.
        ENDIF
        CALL key_value_storage_as_string(meta_storage,'latbc', c_latbc, ierror)
        IF (ierror == SUCCESS) THEN
          IF (TRIM(c_latbc) == 'file') l_latbc = .TRUE.
        ENDIF
        CALL key_value_storage_as_string(meta_storage,'meteogram', c_meteogram, ierror)
        IF (ierror == SUCCESS) THEN
          IF (TRIM(c_meteogram) == 'true') l_meteogram = .TRUE.
        ENDIF
        IF ( l_initc .AND. l_latbc ) THEN
          in_group = groups('ART_AEROSOL','tracer_fg_in','LATBC_PREFETCH_VARS')
        ELSE
          IF ( l_initc ) in_group = groups('ART_AEROSOL','tracer_fg_in')
          IF ( l_latbc ) THEN
            in_group = groups('ART_AEROSOL','LATBC_PREFETCH_VARS')
            IF ( l_limited_area .AND. latbc_config%init_latbc_from_fg .AND. &
              &  art_config%iart_init_aero == 0 ) THEN
              CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
                &         'Do not use init_latbc_from_fg = .TRUE. without initial values!' )
            ENDIF
          ENDIF
        ENDIF
        IF ( l_meteogram ) THEN
            idx_group_mtgrm           = var_groups_dyn%group_id("METEOGRAM")
            in_group(idx_group_mtgrm) = .TRUE.
        END IF
      ENDIF

      ! Required aerosol metadata
      CALL meta_storage%get('moment',moment,ierror)
      IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper',                &
                               &         'Required aerosol metadata moment not present: '//       &
                               &         TRIM(tracer_name_in))
      IF(moment == 3) THEN ! density only needed for mass species, not for number of mixed aerosol
        CALL meta_storage%get('rho',rho,ierror)
        IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper',              &
                                 &         'Required aerosol metadata rho not present: '//        &
                                 &         TRIM(tracer_name_in))
        CALL meta_storage%get('mol_weight',mol_weight,ierror)
        IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper',              &
                                 &         'Required aerosol metadata mol_weight not present: '// &
                                 &         TRIM(tracer_name_in))
        CALL meta_storage%get('sol',sol,ierror)
        IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper',              &
                                 &         'Required aerosol metadata sol not present: '//        &
                                 &         TRIM(tracer_name_in))
        IF (ABS(sol) < 1.0e-20_wp) THEN
          linsol = .TRUE.
        ELSE
          linsol = .FALSE.
        ENDIF
      ELSE ! tracer is specific number
        rho        = UNDEF_REAL_ART
        mol_weight = UNDEF_REAL_ART
        linsol     = .FALSE.
      ENDIF

      CALL meta_storage%get('lfeedback', lfeedback_in, ierror)
      IF (ierror /= SUCCESS) THEN
        lfeedback = .FALSE.
      ELSE
        IF (lfeedback_in == 1) THEN
          lfeedback = .TRUE.
        ELSE
          lfeedback = .FALSE.
        ENDIF
      ENDIF


      SELECT TYPE(meta=>tracer_meta)
        TYPE IS(t_aero_meta)
          meta = create_tracer_metadata_aero(            &
            &              lis_tracer   = lis_tracer,            &
            &              name         = tracer_name,           &
            &              ihadv_tracer = ihadv_tracer,          &
            &              ivadv_tracer = ivadv_tracer,          &
            &              lturb_tracer = lturb_tracer,          &
            &              lconv_tracer = lconv_tracer,          &
            &              ised_tracer  = ised_tracer,           &
            &              ldep_tracer  = ldep_tracer,           &
            &              iwash_tracer = iwash_tracer,          &
            &              moment       = moment,                &
            &              mode         = TRIM(tracer_mode),     &
            &              substance    = TRIM(tracer_sub),      &
            &              rho          = rho,                   &
            &              mol_weight   = mol_weight,            &
            &              linsol       = linsol,                &
            &              lfeedback    = lfeedback)
        CLASS DEFAULT
          CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
            &      'tracer_meta has unknown type')
      END SELECT

    CASE(IART_CHEM_TR)

      ! Set-up post operation for the chemical tracer variables
      CALL meta_storage%get('mol_weight', mol_weight, ierror)
      IF (ierror == SUCCESS) THEN
        post_op_id   = POST_OP_SCALE  
        scale_factor = amd / 1000._wp / mol_weight
      END IF

      CALL key_value_storage_as_string(meta_storage,'c_solve', c_solve, ierror)

      IF (ierror == SUCCESS) THEN
        SELECT CASE (ADJUSTL(TRIM(c_solve)))
            CASE('param')

              ALLOCATE(t_chem_meta_param :: tracer_meta)

            CASE('lt')

              ALLOCATE(t_chem_meta_lt :: tracer_meta)

            CASE('cold')

              ALLOCATE(t_chem_meta_cold :: tracer_meta)

            CASE('OH')

              ALLOCATE(t_chem_meta_OH :: tracer_meta)

            CASE('linoz')

              ALLOCATE(t_chem_meta_linoz :: tracer_meta)

            CASE('simnoy')

              ALLOCATE(t_chem_meta_simnoy :: tracer_meta)

            CASE('mecca')

              ALLOCATE(t_chem_meta_mecca :: tracer_meta)

            CASE('passive')

              ALLOCATE(t_chem_meta_passive :: tracer_meta)
              post_op_id = POST_OP_NONE

            CASE default

              ALLOCATE(t_chem_meta_passive :: tracer_meta)
              post_op_id = POST_OP_NONE
        END SELECT
      ELSE
        ALLOCATE(t_chem_meta_passive :: tracer_meta)
        post_op_id = POST_OP_NONE
      END IF


      IF(TRIM(defcase) == 'prog') THEN
        l_initc = .FALSE.
        l_latbc = .FALSE.
        in_group = groups('ART_CHEMISTRY')
        CALL key_value_storage_as_string(meta_storage,'initc', c_initc, ierror)
        IF (ierror == SUCCESS) THEN
          IF (TRIM(c_initc) == 'file') l_initc = .TRUE.
        ENDIF
        CALL key_value_storage_as_string(meta_storage,'latbc', c_latbc, ierror)
        IF (ierror == SUCCESS) THEN
          IF (TRIM(c_latbc) == 'file') l_latbc = .TRUE.
        ENDIF
        IF ( l_initc .AND. l_latbc ) THEN
          in_group = groups('ART_CHEMISTRY','tracer_fg_in','LATBC_PREFETCH_VARS')
        ELSE
          IF ( l_initc ) in_group = groups('ART_CHEMISTRY','tracer_fg_in')
          IF ( l_latbc ) THEN
            in_group = groups('ART_CHEMISTRY','LATBC_PREFETCH_VARS')
            IF ( l_limited_area .AND. latbc_config%init_latbc_from_fg .AND. &
              &  art_config%iart_init_gas == 0 ) THEN
              CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
                &         'Do not use init_latbc_from_fg = .TRUE. without initial values!' )
            ENDIF
          ENDIF
        ENDIF
      ENDIF

      ! Required chemical metadata
      CALL meta_storage%get('lfeedback', lfeedback_in, ierror)
      IF (ierror /= SUCCESS) THEN
        lfeedback = .FALSE.
      ELSE
        IF (lfeedback_in == 1) THEN
          lfeedback = .TRUE.
        ELSE
          lfeedback = .FALSE.
        ENDIF
      ENDIF



      SELECT TYPE(meta=>tracer_meta)
        CLASS IS(t_chem_meta)
          CALL meta%set_tracer_meta(               &
            &              lis_tracer      = lis_tracer,    &
            &              name            = tracer_name_in,&
            &              ihadv_tracer    = ihadv_tracer,  &
            &              ivadv_tracer    = ivadv_tracer,  &
            &              lturb_tracer    = lturb_tracer,  &
            &              lconv_tracer    = lconv_tracer,  &
            &              ised_tracer     = ised_tracer,   &
            &              ldep_tracer     = ldep_tracer,   &
            &              iwash_tracer    = iwash_tracer,  &
            &              lfeedback       = lfeedback)
        CLASS DEFAULT
          CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
            &      'tracer_meta has unknown type')
      END SELECT

    CASE DEFAULT
      CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
        &      'IART_TRACER_TYPE unknown')
  END SELECT


  ! Optional GRIB2 metadata for ART
  cgrib2_keys = (/ 'productDefinitionTemplate                         ',  &
    &              'constituentType                                   ',  &
    &              'typeOfDistributionFunction                        ',  &
    &              'numberOfModeOfDistribution                        ',  &
    &              'modeNumber                                        ',  &
    &              'numberOfDistributionFunctionParameters            ',  &
    &              'localInformationNumber                            ',  &
    &              'scaledValueOfDistributionFunctionParameter_1      ',  &
    &              'scaleFactorOfDistributionFunctionParameter_1      ',  &
    &              'scaledValueOfDistributionFunctionParameter_2      ',  &
    &              'scaleFactorOfDistributionFunctionParameter_2      ' /)

  ! set variable class according to productDefinitionTemplate
  CALL meta_storage%get(TRIM(cgrib2_keys(1)),igrib2_val,ierror)
  IF (ierror == SUCCESS) THEN
    SELECT CASE (igrib2_val)
    CASE (40)
      var_class = CLASS_CHEM
    CASE (57)
      var_class = CLASS_DISTR
    CASE default
      var_class = CLASS_DEFAULT
    END SELECT
  END IF

  ! set (scalar) integer valued additional GRIB2 keys
  l_modeNumber = .FALSE.
  DO jkey = 2, 7
    CALL meta_storage%get(TRIM(cgrib2_keys(jkey)),igrib2_val,ierror)
    IF (ierror == SUCCESS) THEN
      grib2 = grib2 + t_grib2_int_key(TRIM(cgrib2_keys(jkey)), igrib2_val)
      IF ( TRIM(cgrib2_keys(jkey)) == 'modeNumber' ) l_modeNumber = .TRUE.
    END IF
  END DO
  ! treatment of modeNumber specified via tag modeNumber_list
  IF ( PRESENT(modeNumber_in) .AND. .NOT.l_modeNumber ) THEN
    IF (LEN(TRIM(modeNumber_in)) > 0) THEN
      READ(modeNumber_in, *, iostat=ierror) igrib2_val
      IF (ierror == SUCCESS) THEN
        grib2 = grib2 + t_grib2_int_key('modeNumber', igrib2_val)
      END IF
    END IF
  END IF

! Add tracer to reference container
  IF (TRIM(defcase) == 'prog' .OR. lconv_tracer .OR. lturb_tracer) THEN
    ! No add_ref call for tracer tendencies if process is deactivated
    CALL add_ref(this_list,                                &
      &          target_name,                              &
      &          tracer_name,                              &
      &          ptr_arr(idx)%p_3d,                        &
      &          target_info%hgrid,                        &
      &          target_info%vgrid,                        &
      &          cf,                                       &
      &          grib2,                                    &
      &          ref_idx          = idx,                   &
      &          ldims            = ldims,                 &
      &          loutput          = .TRUE.,                &
      &          lrestart         = lrestart,              &
      &          tlev_source      = 1,                     &
      &          vert_interp      = vert_interp,           &
      &          tracer_info      = tracer_meta,           &
      &          var_class        = var_class,             &
      &          post_op          = post_op(post_op_id,    &
      &                             arg1 = scale_factor),  &
      &          in_group         = in_group)
  ENDIF

  IF(TRIM(defcase) =='prog') THEN
    ! Get all other metadata
    NULLIFY(target_element)
    target_element  => find_list_element (this_list, tracer_name)
    IF (ASSOCIATED(target_element)) THEN
      WRITE (message_text,*) 'ART: Setting optional metadata for '//TRIM(tracer_name_in)//'.'
        CALL message (TRIM(routine)//':art_tracer_def_wrapper', message_text)
    ELSE 
      CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
        &         'Element '//TRIM(tracer_name_in)//' not found in container')
    ENDIF

    target_info_dyn => target_element%info_dyn
    CALL art_read_elements_xml(xmlfile,TRIM(xpath),target_info_dyn%tracer%opt_meta)
    ! Dumping opt_meta-storage
!    This has to be reimplemented in mo_key_value_storage
!    CALL target_info_dyn%tracer%opt_meta%dump("optional metadata for "//TRIM(tracer_name_in))
    ! Small check if required metadata is available (i.e. if reading of meta worked)
    CALL key_value_storage_as_string(target_info_dyn%tracer%opt_meta,'unit',units, ierror)
      IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_tracer_def_wrapper', &
                               &         'Reading optional metadata for '//TRIM(tracer_name_in)  &
                               &       //' failed.')
  ENDIF

  ! Get the number of convection tracers
  IF (lconv_tracer) THEN
    art_config%nconv_tracer = art_config%nconv_tracer + 1
  ENDIF

  ! Get the number of turbulence tracers
  IF (lturb_tracer) THEN
    art_config%nturb_tracer = art_config%nturb_tracer + 1
  END IF

  ! The integers which describe the position of a single variable in the
  ! complete tracer vector needs to be saved only once at the beginning
  IF (TRIM(defcase) == 'prog') then  
    tracer_idx = idx                 
  ENDIF

END SUBROUTINE art_tracer_def_wrapper
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_transport(tracer_name_in, transport_template,                     &
  &                      ivadv_tracer, itype_vlimit, ihadv_tracer, itype_hlimit)
!<
! SUBROUTINE set_transport
! This subroutine sets the transport flags ivadv_tracer, itype_vlimit, 
! ihadv_tracer, itype_hlimit according to the chosen transport template
! Part of Module: mo_art_tracer_def_wrapper
! Author: Daniel Rieger, KIT
! Initial Release: 2016-12-20
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  CHARACTER(LEN=*), INTENT(in) :: &
    &  tracer_name_in,            & !< Name of tracer
    &  transport_template           !< Transport template to be applied
  INTEGER, INTENT(out)         :: &
    &  ivadv_tracer,              & !< Vertical advection scheme
    &  itype_vlimit,              & !< Vertical flux limiter
    &  ihadv_tracer,              & !< Horizontal advection scheme
    &  itype_hlimit                 !< Horizontal flux limiter

  ivadv_tracer = UNDEF_INT_ART
  itype_vlimit = UNDEF_INT_ART
  ihadv_tracer = UNDEF_INT_ART
  itype_hlimit = UNDEF_INT_ART

  IF (TRIM(tolower(transport_template))=='default') THEN
    WRITE (message_text,*) 'ART: WARNING: Tracer '//TRIM(tracer_name_in)  &
                &        //' transport template set to default.'
    CALL message (TRIM(routine)//':set_transport', message_text)
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 2  ! monotonic reconstruction limiter
    ihadv_tracer = 22 ! combination of miura and miura with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF

  IF (TRIM(tolower(transport_template))=='on') THEN ! interim template
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 2  ! monotonic reconstruction limiter
    ihadv_tracer = 22 ! combination of miura and miura with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF

  IF (TRIM(tolower(transport_template))=='off') THEN ! interim template
    ivadv_tracer = 0  ! turned off
    itype_vlimit = 0  ! turned off
    ihadv_tracer = 0  ! turned off
    itype_hlimit = 0  ! turned off
  ENDIF

  IF (TRIM(tolower(transport_template))=='stdchem') THEN ! interim template
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 3  ! positive definte flux limiter
    ihadv_tracer = 22 ! combination of miura and miura with subcycling
    itype_hlimit = 4  ! positive definte flux limiter
  ENDIF

  IF (TRIM(tolower(transport_template))=='stdchem_amip') THEN 
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 1  ! semi-monotonic reconstruction filter
    ihadv_tracer = 22 ! combination of miura and miura with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF
  IF (TRIM(tolower(transport_template))=='stdchem_amip2') THEN 
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 3  ! positive definte flux limiter
    ihadv_tracer = 52 ! combination of hybrid FFSL/Miura3 with subcycling
    itype_hlimit = 4  ! positive definte flux limiter
  ENDIF
   IF (TRIM(tolower(transport_template))=='qv') THEN ! interim template
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 1  ! semi-monotonic reconstruction filter
    ihadv_tracer = 52 ! combination of hybrid FFSL/Miura3 with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF
  IF (TRIM(tolower(transport_template))=='ccmi') THEN ! interim template
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 1  ! semi-monotonic reconstruction filter
    ihadv_tracer = 22 ! combination of miura and miura with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF

  IF (TRIM(tolower(transport_template))=='stdaero') THEN ! interim template
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 2  ! monotonic reconstruction limiter
    ihadv_tracer = 22 ! combination of miura and miura with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF

  IF (TRIM(tolower(transport_template))=='hadv52aero') THEN ! interim template
    ivadv_tracer = 3  ! 3rd order PPM, handles CFL>1
    itype_vlimit = 2  ! monotonic reconstruction limiter
    ihadv_tracer = 52 ! combination of hybrid FFSL/Miura3 with subcycling
    itype_hlimit = 3  ! monotonic flux limiter
  ENDIF

  IF (ivadv_tracer == UNDEF_INT_ART .OR. itype_vlimit == UNDEF_INT_ART .OR. &
    & ihadv_tracer == UNDEF_INT_ART .OR. itype_hlimit == UNDEF_INT_ART) THEN
    CALL finish(TRIM(routine)//':set_transport', &
      &         'Transport template '//TRIM(transport_template)//' of tracer ' &
      &        //TRIM(tracer_name_in)//' unknown.')
  ENDIF

END SUBROUTINE set_transport
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_tracer_def_wrapper

