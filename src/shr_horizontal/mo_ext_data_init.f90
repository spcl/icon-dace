! Initialization/reading reading of external datasets
!
! This module contains read and initialization routines for the external data state.
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ext_data_init

  USE mo_kind,               ONLY: wp
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: inwp, iaes, io3_clim, io3_ape,                                   &
    &                              max_char_length, min_rlcell_int, min_rlcell,                     &
    &                              MODIS, GLOBCOVER2009, GLC2000, SUCCESS, SSTICE_ANA_CLINC,        &
    &                              SSTICE_CLIM, LSS_TERRA
  USE mo_math_constants,     ONLY: dbl_eps, rad2deg
  USE mo_physical_constants, ONLY: ppmv2gg, zemiss_def, tmelt
  USE mo_run_config,         ONLY: msg_level, iforcing, check_uuid_gracefully
  USE mo_impl_constants_grf, ONLY: grf_bdywidth_c
  USE mo_lnd_nwp_config,     ONLY: ntiles_total, ntiles_lnd, ntiles_water, lsnowtile, frlnd_thrhld, &
                                   frlndtile_thrhld, frlake_thrhld, frsea_thrhld, isub_water,       &
                                   isub_seaice, isub_lake, sstice_mode, sst_td_filename,            &
                                   ci_td_filename, itype_lndtbl, c_soil, c_soil_urb, cskinc,        &
                                   lterra_urb, itype_eisa, cr_bsmin, itype_evsl
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config, iprog_aero
  USE mo_extpar_config,      ONLY: itopo, itype_lwemiss, extpar_filename, generate_filename,    &
    &                              generate_td_filename, extpar_varnames_map_file,              &
    &                              n_iter_smooth_topo, i_lctype, nclass_lu, nhori, nmonths_ext, &
    &                              itype_vegetation_cycle, read_nc_via_cdi, pp_sso
  USE mo_initicon_config,    ONLY: icpl_da_sfcevap, dt_ana, icpl_da_skinc, icpl_da_seaice, icpl_da_snowalb
  USE mo_radiation_config,   ONLY: irad_o3, albedo_type, islope_rad,    &
    &                              irad_aero, iRadAeroTegen, iRadAeroART
  USE mo_process_topo,       ONLY: smooth_topo_real_data, postproc_sso, smooth_frland
  USE mo_model_domain,       ONLY: t_patch
  USE mo_exception,          ONLY: message, message_text, finish
  USE mo_grid_config,        ONLY: n_dom, nroot
  USE mo_intp_data_strc,     ONLY: t_int_state
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_mpi,                ONLY: my_process_is_stdio, p_io, p_bcast, &
    &                              p_comm_work_test, p_comm_work, my_process_is_mpi_workroot
  USE mo_sync,               ONLY: global_sum_array
  USE mo_parallel_config,    ONLY: p_test_run, nproma
  USE mo_ext_data_types,     ONLY: t_external_data
  USE mo_nonhydro_types,     ONLY: t_nh_diag
  USE mo_ext_data_state,     ONLY: construct_ext_data, levelname, cellname, o3name, o3unit, &
    &                              nlev_o3, nmonths
  USE mo_master_config,      ONLY: getModelBaseDir
  USE mo_time_config,        ONLY: time_config
  USE mo_io_config,          ONLY: default_read_method
  USE mo_read_interface,     ONLY: openInputFile, closeFile, on_cells, t_stream_id, &
    &                              read_2D, read_2D_int, read_3D_extdim, read_2D_extdim
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf
  USE turb_data,             ONLY: c_lnd, c_sea
  USE mo_util_cdi,           ONLY: get_cdi_varID, test_cdi_varID, read_cdi_2d,     &
    &                              read_cdi_3d, t_inputParameters,                 &
    &                              makeInputParameters, deleteInputParameters,     &
    &                              has_filetype_netcdf
  USE mo_util_uuid_types,    ONLY: t_uuid, uuid_string_length
  USE mo_util_uuid,          ONLY: OPERATOR(==), uuid_unparse
  USE mo_dictionary,         ONLY: t_dictionary
  USE mo_nwp_tuning_config,  ONLY: itune_albedo
  USE mo_cdi,                ONLY: FILETYPE_GRB2, streamOpenRead, streamInqFileType, &
    &                              streamInqVlist, vlistInqVarZaxis, zaxisInqSize,   &
    &                              vlistNtsteps, vlistInqVarGrid, cdiInqAttTxt,    &
    &                              vlistInqVarIntKey, CDI_GLOBAL, gridInqUUID, &
    &                              streamClose, cdiStringError, cdi_undefid
  USE mo_math_gradients,     ONLY: grad_fe_cell
  USE mo_fortran_tools,      ONLY: var_scale
  USE mtime,                 ONLY: datetime, newDatetime, deallocateDatetime,        &
    &                              MAX_DATETIME_STR_LEN, datetimetostring,           &
    &                              OPERATOR(+)
  USE mo_util_mtime,         ONLY: assumePrevMidnight
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,         &
    &                                  calculate_time_interpolation_weights
  USE mo_coupling_config,    ONLY: is_coupled_to_ocean
  USE mo_grid_config,        ONLY: l_scm_mode
  USE mo_scm_nml,            ONLY: i_scm_netcdf
  USE mo_nh_torus_exp,       ONLY: read_ext_scm_nc

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_ext_data_init'

  ! Number of landcover classes provided by external parameter data
  ! Needs to be changed into a variable if landcover classifications
  ! with a different number of classes become available
  INTEGER, PARAMETER :: num_lcc = 23, n_param_lcc = 8

  LOGICAL, ALLOCATABLE :: is_frglac_in(:) !< checks whether the extpar file contains fr_glac
  LOGICAL :: read_netcdf_parallel         !< control variable if NetCDF extpar data are read via parallel NetCDF or cdilib


  PUBLIC :: init_ext_data
  PUBLIC :: init_index_lists
  PUBLIC :: interpol_monthly_mean
  PUBLIC :: diagnose_ext_aggr
  PUBLIC :: vege_clim

!-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !! Init external data for atmosphere
  !!
  !! 1. Build data structure, including field lists and
  !!    memory allocation.
  !! 2. External data are read in from netCDF file or set analytically
  !!
  SUBROUTINE init_ext_data (p_patch, p_int_state, ext_data)

    TYPE(t_patch), INTENT(INOUT)         :: p_patch(:)
    TYPE(t_int_state), INTENT(IN)        :: p_int_state(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)


    ! external parameters for SCM run
    INTEGER  :: soiltyp_scm                       ! soil type
    REAL(wp) :: fr_land_scm                       ! land fraction
    REAL(wp) :: plcov_mx_scm                      ! maximum plant cover
    REAL(wp) :: lai_mx_scm                        ! maximum leaf area index
    REAL(wp) :: rootdp_scm                        ! root depth
    REAL(wp) :: rsmin_scm                         ! minimum stomata resistance
    REAL(wp) :: z0_scm                            ! roughness length
    REAL(wp) :: topo_scm                          ! topographic height
    REAL(wp) :: emis_rad_scm                      ! emissivity
    REAL(wp) :: lu_class_fr_scm(num_lcc)          ! lu_class_fraction
    CHARACTER(len=max_char_length) :: lctype_scm  ! type of data source for land use
    
    INTEGER :: jg, ilcc
    INTEGER :: cdi_extpar_id(n_dom)  !< CDI stream ID (for each domain)
    INTEGER :: cdi_filetype(n_dom)   !< CDI filetype (for each domain)
    ! dictionary which maps internal variable names onto
    ! GRIB2 shortnames or NetCDF var names.
    TYPE (t_dictionary) :: extpar_varnames_dict

    TYPE(datetime), POINTER :: this_datetime
    CHARACTER(len=*), PARAMETER :: routine = modname//':init_ext_data'
    LOGICAL :: is_mpi_workroot

    is_mpi_workroot = my_process_is_mpi_workroot()

    !-------------------------------------------------------------------------
    CALL message(routine, 'Start')

    !-------------------------------------------------------------------------
    !  1.  inquire external files for their data structure
    !-------------------------------------------------------------------------

    ALLOCATE(is_frglac_in(n_dom))
    ! Set default value for is_frglac_in. Will be overwritten, if external data
    ! contain fr_glac
    is_frglac_in(1:n_dom) = .FALSE.

    ! initialize stream-IDs as "uninitialized"
    cdi_extpar_id(:) = cdi_undefid
    ! open CDI stream (files):
    IF (iforcing == inwp) CALL inquire_external_files(p_patch, cdi_extpar_id, cdi_filetype)

    ! read the map file (internal -> GRIB2) into dictionary data structure:
    CALL extpar_varnames_dict%init(.FALSE.)
    IF (ANY(cdi_filetype(:) == FILETYPE_GRB2)) THEN
      IF (extpar_varnames_map_file /= ' ') &
        & CALL extpar_varnames_dict%loadfile(TRIM(extpar_varnames_map_file))
      read_netcdf_parallel = .FALSE. ! GRIB2 can only be read using cdi library
    ELSE IF (read_nc_via_cdi) THEN
      read_netcdf_parallel = .FALSE.
    ELSE
      read_netcdf_parallel = .TRUE.
    END IF

    !------------------------------------------------------------------
    !  2.  construct external fields for the model
    !------------------------------------------------------------------

    ! top-level procedure for building data structures for
    ! external data.
    
    IF (i_scm_netcdf > 0) THEN
      nclass_lu = num_lcc ! 3rd dim of lu_class_fraction, has to agree with num_lcc
    ENDIF
   
    CALL construct_ext_data(p_patch, ext_data)

    !-------------------------------------------------------------------------
    !  3.  read the data into the fields
    !-------------------------------------------------------------------------

    ! Check, whether external data should be read from file

    SELECT CASE(itopo)

    CASE(0) ! itopo, do not read external data
      !
      CALL message(routine,'Running with analytical topography' )
      !
      ! initalize external data with meaningful data, in the case that they
      ! are not read in from file.
      IF ( iforcing == inwp ) THEN

        ! SCM netcdf input or surface parameters
        IF ( l_scm_mode .AND. (i_scm_netcdf==1) ) THEN

          !read external parameters from netCDF file
          ! TODO: read external parameters for unified SCM formal 'uf'
          CALL read_ext_scm_nc(num_lcc,soiltyp_scm,fr_land_scm,plcov_mx_scm,lai_mx_scm,rootdp_scm, & 
            &                  rsmin_scm,z0_scm,topo_scm,emis_rad_scm,lu_class_fr_scm,lctype_scm)
          DO jg = 1, n_dom
            !set external parameters
            ext_data(jg)%atm%fr_land(:,:)     = fr_land_scm  ! land fraction
            IF (fr_land_scm >= 0.5_wp ) THEN
              ext_data(jg)%atm%llsm_atm_c(:,:)= .TRUE.       ! land-sea mask
            ELSE
              ext_data(jg)%atm%llsm_atm_c(:,:)= .FALSE.
            ENDIF
            ext_data(jg)%atm%llake_c(:,:)     = .FALSE.      ! lake mask
            ext_data(jg)%atm%plcov_mx(:,:)    = plcov_mx_scm ! plant cover
            ext_data(jg)%atm%lai_mx(:,:)      = lai_mx_scm   ! max Leaf area index
            ext_data(jg)%atm%rootdp(:,:)      = rootdp_scm   ! root depth
            ext_data(jg)%atm%rsmin(:,:)       = rsmin_scm    ! minimal stomata resistence
            ext_data(jg)%atm%soiltyp(:,:)     = soiltyp_scm  ! soil type
            ext_data(jg)%atm%z0(:,:)          = z0_scm       ! roughness length
            ext_data(jg)%atm%topography_c(:,:)= topo_scm     ! topographic height
            ext_data(jg)%atm%emis_rad(:,:)    = emis_rad_scm ! emissivity
            !
            DO ilcc=1, num_lcc
              ext_data(jg)%atm%lu_class_fraction(:,:,ilcc) = lu_class_fr_scm(ilcc) ! fraction of LU class
            ENDDO
            IF (TRIM(lctype_scm) .EQ. "GLC2000") THEN
              i_lctype(jg) = GLC2000
            ELSE IF (TRIM(lctype_scm) .EQ. "GLOBCOVER2009" ) THEN
              i_lctype(jg) = GLOBCOVER2009
            ELSE
              CALL finish(routine,'Unknown landcover data source')
            ENDIF
            !ext_data(jg)%atm%i_lc_water        = 21
  
            !Special setup for tiles
            ext_data(jg)%atm%soiltyp_t(:,:,:) = soiltyp_scm ! soil type
            ext_data(jg)%atm%frac_t(:,:,:)    = 0._wp       ! set all tiles to 0
            ext_data(jg)%atm%frac_t(:,:,isub_water) = 1._wp ! set only ocean to 1
            ext_data(jg)%atm%lc_class_t(:,:,:) = 1          ! land cover class
            
          END DO

        ELSE
          DO jg = 1, n_dom
            ext_data(jg)%atm%fr_land(:,:)     = 0._wp       ! land fraction
            ext_data(jg)%atm%llsm_atm_c(:,:)  = .FALSE.     ! land-sea mask
            ext_data(jg)%atm%llake_c(:,:)     = .FALSE.     ! lake mask
!
            ext_data(jg)%atm%urb_isa(:,:)     = 0._wp       ! impervious surface area fraction of the urban canopy
            IF (lterra_urb) THEN
              ext_data(jg)%atm%fr_paved(:,:)    = 0._wp       ! impervious surface area (ISA) fraction
              ext_data(jg)%atm%urb_ai(:,:)      = 2._wp       ! surface area index of the urban canopy
              ext_data(jg)%atm%urb_alb_red(:,:) = 0.9_wp      ! albedo reduction factor for the urban canopy
              ext_data(jg)%atm%urb_fr_bld(:,:)  = 0.667_wp    ! building area fraction with respect to urban tile
              ext_data(jg)%atm%urb_h2w(:,:)     = 1.5_wp      ! street canyon H/W ratio
              ext_data(jg)%atm%urb_h_bld(:,:)   = 15._wp      ! building height
              ext_data(jg)%atm%urb_alb_th(:,:)  = 0.14_wp     ! thermal albedo of urban material
              ext_data(jg)%atm%urb_alb_so(:,:)  = 0.101_wp    ! solar albedo of urban material
              ext_data(jg)%atm%urb_hcap(:,:)    = 1250000._wp ! volumetric heat capacity of urban material
              ext_data(jg)%atm%urb_hcon(:,:)    = 0.767_wp    ! thermal conductivity of urban material
              ext_data(jg)%atm%ahf(:,:)         = 0._wp       ! anthropogenic heat flux
            ENDIF
!
            ext_data(jg)%atm%plcov_mx(:,:)    = 0.5_wp      ! plant cover
            ext_data(jg)%atm%lai_mx(:,:)      = 3._wp       ! max Leaf area index
            ext_data(jg)%atm%rootdp(:,:)      = 1._wp       ! root depth
            ext_data(jg)%atm%skinc(:,:)       = 30._wp      ! skin conductivity
            ext_data(jg)%atm%rsmin(:,:)       = 150._wp     ! minimal stomata resistence
            ext_data(jg)%atm%soiltyp(:,:)     = 8           ! soil type
            ext_data(jg)%atm%z0(:,:)          = 0.001_wp    ! roughness length
            ext_data(jg)%atm%topography_c(:,:)= 0.0_wp      ! topographic height
            i_lctype(jg) = GLOBCOVER2009
  
            !Special setup for tiles
            ext_data(jg)%atm%soiltyp_t(:,:,:) = 8           ! soil type
            ext_data(jg)%atm%frac_t(:,:,:)    = 0._wp       ! set all tiles to 0
            ext_data(jg)%atm%frac_t(:,:,isub_water) = 1._wp ! set only ocean to 1
            ext_data(jg)%atm%lc_class_t(:,:,:) = 1          ! land cover class
  
          END DO

        ENDIF
     
        IF (l_scm_mode) THEN
          ! necessary call to initialize some extpar variables without reading the extpar file
          CALL read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
            &                     extpar_varnames_dict)
          CALL message(routine,'read_ext_data_atm completed' )
        ENDIF

        DO jg = 1, n_dom
          ext_data(jg)%atm%emis_rad(:,:)    = zemiss_def ! longwave surface emissivity
        END DO

        ! call read_ext_data_atm to read O3
        IF ( irad_o3 == io3_clim .OR. irad_o3 == io3_ape .OR. sstice_mode == SSTICE_CLIM) THEN
          CALL read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
            &                     extpar_varnames_dict)
          CALL message(routine,'read_ext_data_atm completed' )
        END IF

      END IF


    CASE(1) ! itopo, read external data from file

      CALL message(routine,'Start reading external data from file' )

      CALL read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
        &                     extpar_varnames_dict)

      CALL message(routine,'Finished reading external data' )


      DO jg = 1, n_dom
         IF ( iforcing == inwp ) THEN
            IF (pp_sso > 0) THEN
               CALL postproc_sso ( p_patch(jg)                  ,&
                &                 p_int_state(jg)               ,&
                &                 ext_data(jg)%atm%fr_glac      ,&
                &                 ext_data(jg)%atm%topography_c ,&
                &                 ext_data(jg)%atm%sso_stdh     ,&
                &                 ext_data(jg)%atm%sso_sigma )
            ENDIF
         ENDIF

         ! topography smoothing
         IF ( n_iter_smooth_topo(jg) > 0 ) THEN
            IF ( iforcing == inwp ) THEN
               CALL smooth_topo_real_data ( p_patch(jg)                          ,&
                 &                          p_int_state(jg)                      ,&
                 &                          ext_data(jg)%atm%fr_land             ,&
                 &                          ext_data(jg)%atm%topography_c        ,&
                 &                          fr_lake = ext_data(jg)%atm%fr_lake   ,&
                 &                          sso_stdh = ext_data(jg)%atm%sso_stdh )
            ELSE IF (iforcing == iaes) THEN
               CALL smooth_topo_real_data ( p_patch(jg)                          ,&
                 &                          p_int_state(jg)                      ,&
                 &                          ext_data(jg)%atm%fr_land             ,&
                 &                          ext_data(jg)%atm%topography_c )
            ENDIF
         ENDIF ! n_iter_smooth_topo(jg) > 0

         IF ( iforcing == inwp ) THEN
           ! smooth land fraction for adaptive tuning of sea ice bottom heat flux and sea ice albedo
           IF (icpl_da_seaice >= 2 .OR. icpl_da_snowalb >= 2) THEN
             CALL smooth_frland (p_patch(jg),                &
               &                 p_int_state(jg),            &
               &                 ext_data(jg)%atm%fr_land,   &
               &                 ext_data(jg)%atm%fr_land_smt)
           ENDIF
         ENDIF

         ! calculate gradient of orography for resolved surface drag
         !
         call grad_fe_cell  ( ext_data(jg)%atm%topography_c, &
           &                  p_patch(jg),                   &
           &                  p_int_state(jg),               &
           &                  ext_data(jg)%atm%grad_topo )
      END DO


      IF ( iforcing == inwp ) THEN

        ! Get interpolated ndviratio, alb_dif, albuv_dif and albni_dif. Interpolation
        ! is done in time, based on ini_datetime (midnight). Fields are updated on a
        ! daily basis.

        ! When initializing the model we set the target hour to 0 (midnight) as well.
        ! When restarting, the target interpolation time must be set to cur_datetime
        ! midnight.
        !

        this_datetime => newDatetime(time_config%tc_current_date)

        ! always assume midnight
        DO jg = 1, n_dom
          CALL interpol_monthly_mean(p_patch(jg),                        &! in
            &                        assumePrevMidnight(this_datetime),  &! in
            &                        ext_data(jg)%atm_td%ndvi_mrat,      &! in
            &                        ext_data(jg)%atm%ndviratio          )! out
          IF (itype_vegetation_cycle > 1) THEN
            CALL interpol_monthly_mean(p_patch(jg),                      &! in
              &                        assumePrevMidnight(this_datetime),&! in
              &                        ext_data(jg)%atm_td%t2m_m,        &! in
              &                        ext_data(jg)%atm%t2m_clim,        &! out
              &                        ext_data(jg)%atm%t2m_climgrad     )! optional out
          ENDIF
        ENDDO

        IF ( albedo_type == MODIS) THEN
          DO jg = 1, n_dom
            CALL interpol_monthly_mean(p_patch(jg),                      &! in
              &                        assumePrevMidnight(this_datetime),&! in
              &                        ext_data(jg)%atm_td%alb_dif,      &! in
              &                        ext_data(jg)%atm%alb_dif          )! out

            CALL interpol_monthly_mean(p_patch(jg),                      &! in
              &                        assumePrevMidnight(this_datetime),&! in
              &                        ext_data(jg)%atm_td%albuv_dif,    &! in
              &                        ext_data(jg)%atm%albuv_dif        )! out

            CALL interpol_monthly_mean(p_patch(jg),                      &! in
              &                        assumePrevMidnight(this_datetime),&! in
              &                        ext_data(jg)%atm_td%albni_dif,    &! in
              &                        ext_data(jg)%atm%albni_dif        )! out
          ENDDO
        ENDIF  ! albedo_type

        IF ( itype_lwemiss == 2) THEN
          DO jg = 1, n_dom
            CALL interpol_monthly_mean(p_patch(jg),                      &! in
              &                        assumePrevMidnight(this_datetime),&! in
              &                        ext_data(jg)%atm_td%lw_emiss,     &! in
              &                        ext_data(jg)%atm%emis_rad         )! out
          ENDDO
        ENDIF  ! lwemiss

        ! clean up
        CALL deallocateDatetime(this_datetime)

      END IF ! inwp

    CASE DEFAULT ! itopo

      CALL finish(routine, 'topography selection not supported' )

    END SELECT ! itopo

    ! close CDI stream (file):
    IF (is_mpi_workroot) THEN
      DO jg=1,n_dom
        IF (cdi_extpar_id(jg) /= cdi_undefid) &
             CALL streamClose(cdi_extpar_id(jg))
      END DO
    END IF

    ! destroy variable name dictionary:
    CALL extpar_varnames_dict%finalize()

  END SUBROUTINE init_ext_data


  !-------------------------------------------------------------------------
  ! Open ExtPar file and investigate the data structure of the
  ! external parameters.
  !
  ! Note: This subroutine opens the file and returns a CDI file ID.
  !-------------------------------------------------------------------------
  SUBROUTINE inquire_extpar_file(p_patch, jg, cdi_extpar_id, cdi_filetype, &
    &                            is_frglac_in)
    TYPE(t_patch), INTENT(IN)      :: p_patch(:)
    INTEGER,       INTENT(IN)      :: jg
    INTEGER,       INTENT(INOUT)   :: cdi_extpar_id     !< CDI stream ID
    INTEGER,       INTENT(INOUT)   :: cdi_filetype      !< CDI filetype
    LOGICAL,       INTENT(OUT)     :: is_frglac_in      !< check for fr_glac in Extpar file

    ! local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//'::inquire_extpar_file'
    INTEGER                 :: mpi_comm, vlist_id, lu_class_fraction_id, zaxis_id, var_id
    INTEGER                 :: horizon_id
    LOGICAL                 :: l_exist
    CHARACTER(filename_max) :: extpar_file !< file name for reading in
    INTEGER :: extpar_file_namelen

    TYPE(t_uuid)            :: extpar_uuidOfHGrid             ! uuidOfHGrid contained in the
                                                              ! extpar file

    CHARACTER(len=uuid_string_length) :: grid_uuid_unparsed   ! unparsed grid uuid (human readable)
    CHARACTER(len=uuid_string_length) :: extpar_uuid_unparsed ! same for extpar-file uuid

    LOGICAL                 :: lmatch                         ! for comparing UUIDs

    INTEGER                 :: cdiGridID

    INTEGER :: lu_var_id, localInformationNumber
    INTEGER :: ret
    CHARACTER(len=max_char_length) :: rawdata_attr

    !---------------------------------------------!
    ! Check validity of external parameter file   !
    !---------------------------------------------!
    IF (my_process_is_mpi_workroot()) THEN
      ! generate file name
      extpar_file = generate_filename(extpar_filename,                   &
        &                             getModelBaseDir(),                 &
        &                             TRIM(p_patch(jg)%grid_filename),   &
        &                             nroot,                             &
        &                             p_patch(jg)%level, p_patch(jg)%id)
      extpar_file_namelen = LEN_TRIM(extpar_file)
      CALL message(routine, "extpar_file = "//extpar_file(1:extpar_file_namelen))

      INQUIRE (FILE=extpar_file, EXIST=l_exist)
      IF (.NOT.l_exist)  CALL finish(routine,'external data file is not found.')

      ! open file
      cdi_extpar_id = streamOpenRead(extpar_file(1:extpar_file_namelen))
      IF (cdi_extpar_id < 0) THEN
        WRITE (message_text, '(133a)') "Cannot open external parameter file ", &
             cdiStringError(cdi_extpar_id)
        CALL finish(routine, TRIM(message_text))
      END IF
      cdi_filetype  = streamInqFileType(cdi_extpar_id)

      IF (islope_rad(jg) >= 2) THEN
      ! get the number of horizon sectors
        horizon_id = get_cdi_varID(cdi_extpar_id, "HORIZON")
        vlist_id   = streamInqVlist(cdi_extpar_id)
        zaxis_id   = vlistInqVarZaxis(vlist_id, horizon_id)
        nhori = zaxisInqSize(zaxis_id)
        !$ACC UPDATE DEVICE(nhori) ASYNC(1)
        WRITE(message_text,'(A,I4)')  &
          & 'Number of horizon sectors in external data file = ', nhori
        CALL message(routine, TRIM(message_text))
      ENDIF

      ! get the number of landuse classes
      lu_class_fraction_id = get_cdi_varID(cdi_extpar_id, "LU_CLASS_FRACTION")
      vlist_id             = streamInqVlist(cdi_extpar_id)
      zaxis_id             = vlistInqVarZaxis(vlist_id, lu_class_fraction_id)
      IF (l_scm_mode) THEN
        nclass_lu(jg)      = num_lcc
      ELSE
        nclass_lu(jg)      = zaxisInqSize(zaxis_id)
      ENDIF

      ! get time dimension from external data file
      nmonths_ext(jg)      = vlistNtsteps(vlist_id)

      ! make sure that num_lcc is equal to nclass_lu. If not, then the internal
      ! land-use lookup tables and the external land-use class field are inconsistent.

      IF (nclass_lu(jg) /= num_lcc) THEN
        WRITE(message_text,'(A,I3,A,I3)')  &
          & 'Number of land-use classes in external file ', nclass_lu(jg), &
          & ' does not match ICON-internal value num_lcc ', num_lcc
        CALL finish(routine,TRIM(message_text))
      ENDIF

      IF ( msg_level>10 ) THEN
        WRITE(message_text,'(A,I4)')  &
          & 'Number of land_use classes in external data file = ', nclass_lu(jg)
        CALL message(routine,message_text)

        WRITE(message_text,'(A,I4)')  &
          & 'Number of months in external data file = ', nmonths_ext(jg)
        CALL message(routine,message_text)
      ENDIF


      ! Compare UUID of external parameter file with UUID of grid.
      !
      ! get horizontal grid UUID contained in extpar file
      ! use lu_class_fraction as sample field
      cdiGridID = vlistInqVarGrid(vlist_id, lu_class_fraction_id)
      CALL gridInqUUID(cdiGridID, extpar_uuidOfHGrid%DATA)
      !
      ! --- compare UUID of horizontal grid file with UUID from extpar file
      lmatch = (p_patch(jg)%grid_uuid == extpar_uuidOfHGrid)

      IF (.NOT. lmatch) THEN
        CALL uuid_unparse(p_patch(jg)%grid_uuid, grid_uuid_unparsed)
        CALL uuid_unparse(extpar_uuidOfHGrid   , extpar_uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from gridfile: ', TRIM(grid_uuid_unparsed)
        CALL message(routine,message_text)
        WRITE(message_text,'(a,a)') 'uuidOfHgrid from extpar file: ', TRIM(extpar_uuid_unparsed)
        CALL message(routine,message_text)

        WRITE(message_text,'(a)') 'Extpar file and horizontal grid file do not match!'
        IF (check_uuid_gracefully) THEN
          CALL message(routine, TRIM(message_text))
        ELSE
          CALL finish(routine, TRIM(message_text))
        END IF
      ENDIF



      ! Determine which data source has been used to generate the
      ! external perameters: For NetCDF format, we check the
      ! global attribute "rawdata". For GRIB2 format we check the
      ! key "localInformationNumber".
      IF (has_filetype_netcdf(cdi_filetype)) THEN
        ret      = cdiInqAttTxt(vlist_id, CDI_GLOBAL, 'rawdata', max_char_length, rawdata_attr)
        IF (INDEX(rawdata_attr,'GLC2000') /= 0) THEN
          i_lctype(jg) = GLC2000
        ELSE IF (INDEX(rawdata_attr,'GLOBCOVER2009') /= 0) THEN
          i_lctype(jg) = GLOBCOVER2009
        ELSE
          CALL finish(routine,'Unknown landcover data source')
        ENDIF
      ELSE IF (cdi_filetype == FILETYPE_GRB2) THEN
        lu_var_id              = get_cdi_varID(cdi_extpar_id, 'LU_CLASS_FRACTION')
        localInformationNumber = vlistInqVarIntKey(vlist_id, lu_var_id, "localInformationNumber")
        SELECT CASE (localInformationNumber)
        CASE (2)  ! 2 = GLC2000
          i_lctype(jg) = GLC2000
        CASE (1)  ! 1 = ESA GLOBCOVER
          i_lctype(jg) = GLOBCOVER2009
        CASE DEFAULT
          CALL finish(routine,'Unknown landcover data source')
        END SELECT
      END IF

      ! Check whether external parameter file contains MODIS albedo-data
      IF ( albedo_type == MODIS ) THEN
        IF ( (test_cdi_varID(cdi_extpar_id, 'ALB')   == -1) .OR.    &
          &  (test_cdi_varID(cdi_extpar_id, 'ALNID') == -1) .OR.    &
          &  (test_cdi_varID(cdi_extpar_id, 'ALUVD') == -1) ) THEN
          CALL finish(routine,'MODIS albedo fields missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Check whether external parameter file contains monthly longwave surface emissivity data
      IF ( itype_lwemiss == 2 ) THEN
        IF (test_cdi_varID(cdi_extpar_id, 'EMISS')   == -1) THEN
          CALL finish(routine,'Monthly longwave surface emissvity data missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Check whether external parameter file contains SST climatology
      IF ( sstice_mode == SSTICE_ANA_CLINC ) THEN
        IF ( test_cdi_varID(cdi_extpar_id, 'T_SEA')  == -1 ) THEN
          CALL finish(routine,'SST climatology missing in '//TRIM(extpar_filename))
        ENDIF
      ENDIF

      ! Search for glacier fraction in Extpar file
      !
      IF (has_filetype_netcdf(cdi_filetype)) THEN
        var_id = test_cdi_varID(cdi_extpar_id,'ICE')
      ELSE IF (cdi_filetype == FILETYPE_GRB2) THEN
        var_id = test_cdi_varID(cdi_extpar_id,'FR_ICE')
      ENDIF
      IF (var_id == -1) THEN
        is_frglac_in = .FALSE.
      ELSE
        is_frglac_in = .TRUE.
      ENDIF

    ENDIF ! my_process_is_mpi_workroot()

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF
    ! broadcast nhori from I-Pe to WORK Pes
    CALL p_bcast(nhori, p_io, mpi_comm)
    !$ACC UPDATE DEVICE(nhori) ASYNC(1)
    ! broadcast nclass_lu from I-Pe to WORK Pes
    CALL p_bcast(nclass_lu(jg), p_io, mpi_comm)
    ! broadcast nmonths from I-Pe to WORK Pes
    CALL p_bcast(nmonths_ext(jg), p_io, mpi_comm)
    ! broadcast is_frglac_in from I-Pe to WORK Pes
    CALL p_bcast(is_frglac_in, p_io, mpi_comm)
    ! broadcast i_lctype from I-Pe to WORK Pes
    CALL p_bcast(i_lctype(jg), p_io, mpi_comm)
    ! broadcast cdi filetype
    CALL p_bcast(cdi_filetype, p_io, mpi_comm)

  END SUBROUTINE inquire_extpar_file

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE inquire_external_files(p_patch, cdi_extpar_id, cdi_filetype)

    !-------------------------------------------------------
    !
    ! open netcdf files and investigate the data structure
    ! of the external parameters
    !
    !-------------------------------------------------------

    TYPE(t_patch), INTENT(IN)      :: p_patch(:)
    INTEGER,       INTENT(INOUT)   :: cdi_extpar_id(:)  !< CDI stream ID
    INTEGER,       INTENT(INOUT)   :: cdi_filetype(:)   !< CDI filetype

    INTEGER :: jg, mpi_comm
    INTEGER :: no_cells
    INTEGER :: ncid, dimid

    LOGICAL :: l_exist

    CHARACTER(len=*), PARAMETER :: &
      routine = modname//':inquire_external_files'

    CHARACTER(filename_max) :: ozone_file  !< file name for reading in

!--------------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    DO jg= 1,n_dom

      !------------------------------------------------!
      ! 1. Check validity of external parameter file   !
      !------------------------------------------------!

      IF ( itopo == 1 ) THEN
        CALL inquire_extpar_file(p_patch, jg, cdi_extpar_id(jg), cdi_filetype(jg), &
          &                      is_frglac_in(jg))
      END IF

      !------------------------------------------------!
      ! 2. Check validity of ozone file                !
      !------------------------------------------------!

      ! default values for nlev_o3 and nmonths
      nlev_o3 = 1
      nmonths   = 1

      O3 : IF ( irad_o3 == io3_clim .OR. irad_o3 == io3_ape ) THEN

        IF(iforcing == inwp .AND. irad_o3 == io3_ape) THEN
          levelname = 'level'
          cellname  = 'ncells'
          o3name    = 'O3'
          o3unit    = 'g/g'
        ELSE ! o3_clim
          levelname = 'plev'
          cellname  = 'ncells'
          o3name    = 'O3'
          o3unit    = 'g/g' !this unit ozon will have after being read out and converted from ppmv
        ENDIF

        IF_IO : IF(my_process_is_stdio()) THEN

          WRITE(ozone_file,'(a,i2.2,a)') 'o3_icon_DOM',jg,'.nc'

          ! Note resolution assignment is done per script by symbolic links

          INQUIRE (FILE=ozone_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            WRITE(0,*) 'DOMAIN=',jg
            CALL finish(routine,'ozone file of domain is not found.')
          ENDIF

          !
          ! open file
          !
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid), routine)
          WRITE(0,*)'open ozone file'

          !
          ! get number of cells
          !
          CALL nf(nf_inq_dimid (ncid, TRIM(cellname), dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, no_cells), routine)
          WRITE(0,*)'number of cells are', no_cells

          !
          ! check the number of cells and verts
          !
          IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
            CALL finish(routine,&
              & 'Number of patch cells and cells in ozone file do not match.')
          ENDIF

          !
          ! check the time structure
          !
          CALL nf(nf_inq_dimid (ncid, 'time', dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nmonths), routine)
          WRITE(message_text,'(A,I4)')  &
            & 'Number of months in ozone file = ', nmonths
          CALL message(routine,message_text)

          !
          ! check the vertical structure
          !
          CALL nf(nf_inq_dimid (ncid,TRIM(levelname), dimid), routine)
          CALL nf(nf_inq_dimlen(ncid, dimid, nlev_o3), routine)

          WRITE(message_text,'(A,I4)')  &
            & 'Number of pressure levels in ozone file = ', nlev_o3
          CALL message(routine,message_text)

          !
          ! close file
          !
          CALL nf(nf_close(ncid), routine)

        END IF IF_IO ! pe

        CALL p_bcast(nlev_o3, p_io, mpi_comm)
        CALL p_bcast(nmonths,   p_io, mpi_comm)

      END IF O3 !o3

    ENDDO ! ndom

  END SUBROUTINE inquire_external_files

  !-------------------------------------------------------------------------
  !! Read atmospheric external data from netcdf
  !!
  SUBROUTINE read_ext_data_atm (p_patch, ext_data, nlev_o3, cdi_extpar_id, &
    &                           extpar_varnames_dict)

    TYPE(t_patch), TARGET, INTENT(IN)    :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)
    INTEGER,               INTENT(IN)    :: nlev_o3

    INTEGER,               INTENT(IN)    :: cdi_extpar_id(:)      !< CDI stream ID
    TYPE (t_dictionary),   INTENT(IN)    :: extpar_varnames_dict  !< variable names dictionary (for GRIB2)

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_ext_data_atm'
    ! input file for topography_c for mpi-physics
    CHARACTER(len=max_char_length) :: land_sso_fn, land_frac_fn

    CHARACTER(filename_max) :: ozone_file  !< file name for reading in
    CHARACTER(filename_max) :: sst_td_file !< file name for reading in
    CHARACTER(filename_max) :: ci_td_file  !< file name for reading in

    INTEGER :: jg, jc, jb, i, mpi_comm, ilu,im
    INTEGER :: jk
    INTEGER :: ncid, varid
    TYPE(t_stream_id) :: stream_id

    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk   !> blocks
    INTEGER :: i_startidx, i_endidx   !< slices
    INTEGER :: i_nchdom               !< domain index

    REAL(wp):: zdummy_o3lev(nlev_o3) ! will be used for pressure and height levels
    REAL(wp):: albfac, albthresh             ! for MODIS albedo tuning

    REAL(wp), DIMENSION(num_lcc*n_param_lcc)         :: lu_glc2000   ! < lookup table landuse class GLC2000
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009   ! < lookup table landuse class GlobCover2009
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009_v2 ! < modified lookup table landuse class GlobCover2009
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009_v3 ! < even less evaporating lookup table landuse class GlobCover2009
    REAL(wp), DIMENSION(num_lcc*n_param_lcc), TARGET :: lu_gcv2009_v4 ! < retuned lookup table landuse class GlobCover2009
    REAL(wp), POINTER :: lu_gcv(:)

    LOGICAL :: l_exist
    CHARACTER(filename_max) :: extpar_file

    TYPE(t_inputParameters) :: parameters
    LOGICAL :: is_mpi_workroot

    is_mpi_workroot = my_process_is_mpi_workroot()

!                    z0         pcmx      laimx rd      rsmin      snowalb snowtile skinc
!
 DATA lu_glc2000 /   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp, 100._wp, & ! evergreen broadleaf forest
                 &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 150.0_wp,  0.31_wp, 1._wp, 100._wp, & ! deciduous broadleaf closed forest
                 &   0.15_wp,  0.8_wp,  4.0_wp, 2.0_wp, 150.0_wp,  0.31_wp, 1._wp, 100._wp, & ! deciduous broadleaf open   forest
                 &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.27_wp, 1._wp, 100._wp, & ! evergreen needleleaf forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.33_wp, 1._wp, 100._wp, & ! deciduous needleleaf forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 150.0_wp,  0.29_wp, 1._wp, 100._wp, & ! mixed leaf trees
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! fresh water flooded trees
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! saline water flooded trees
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! mosaic tree / natural vegetation
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp, 150.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! burnt tree cover
                 &   0.20_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! evergreen shrubs closed-open
                 &   0.15_wp,  0.8_wp,  1.5_wp, 2.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! decidous shrubs closed-open
                 &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp,  40.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! herbaceous vegetation closed-open
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  40.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! sparse herbaceous or grass
                 &   0.05_wp,  0.8_wp,  2.0_wp, 0.4_wp,  40.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! flooded shrubs or herbaceous
                 &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! cultivated & managed areas
                 &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! mosaic crop / tree / natural vegetation
                 &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 100.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! mosaic crop / shrub / grass
                 &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! bare areas
                 &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, 100._wp, & ! water
                 &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! snow & ice
                 &   1.00_wp,  0.2_wp,  1.0_wp, 0.6_wp, 120.0_wp,  -1.0_wp, 1._wp, 100._wp, & ! artificial surface
                 &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp,  40.0_wp,  -1.0_wp,-1._wp, 100._wp  / ! undefined

 DATA lu_gcv2009 /   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp,  30._wp, & ! irrigated croplands
                 &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp,  30._wp, & ! rainfed croplands
                 &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 120.0_wp,  0.55_wp, 1._wp,  10._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                 &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 100.0_wp,  0.72_wp, 1._wp,  30._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp,  50._wp, & ! closed broadleaved evergreen forest
                 &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 150.0_wp,  0.31_wp, 1._wp,  50._wp, & ! closed broadleaved deciduous forest
                 &   0.15_wp,  0.8_wp,  4.0_wp, 2.0_wp, 150.0_wp,  0.31_wp, 1._wp,  30._wp, & ! open broadleaved deciduous forest
                 &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.27_wp, 1._wp,  50._wp, & ! closed needleleaved evergreen forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 150.0_wp,  0.33_wp, 1._wp,  50._wp, & ! open needleleaved deciduous forest
                 &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 150.0_wp,  0.29_wp, 1._wp,  50._wp, & ! mixed broadleaved and needleleaved forest
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.60_wp, 1._wp,  30._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                 &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.65_wp, 1._wp,  10._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                 &   0.15_wp,  0.8_wp,  2.5_wp, 1.5_wp, 120.0_wp,  0.65_wp, 1._wp,  50._wp, & ! closed to open shrubland
                 &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp,  40.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open herbaceous vegetation
                 &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  40.0_wp,  0.76_wp, 1._wp,  10._wp, & ! sparse vegetation
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed to open forest regulary flooded
                 &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed forest or shrubland permanently flooded
                 &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  40.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open grassland regularly flooded
                 &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 120.0_wp,  0.50_wp, 1._wp, 200._wp, & ! artificial surfaces
                 &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 120.0_wp,  0.76_wp, 1._wp, 200._wp, & ! bare areas
                 &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp,-1._wp, 200._wp, & ! water bodies
                 &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 200._wp, & ! permanent snow and ice
                 &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp, 200._wp  / ! undefined

! Tuned version of gcv2009 based on IFS values (Juergen Helmert und Martin Koehler)
 DATA lu_gcv2009_v2 /  0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 180.0_wp,  0.72_wp, 1._wp,  30._wp, & ! irrigated croplands
                   &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 140.0_wp,  0.72_wp, 1._wp,  30._wp, & ! rainfed croplands
                   &   0.25_wp,  0.8_wp,  3.0_wp, 1.0_wp, 130.0_wp,  0.55_wp, 1._wp,  10._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.07_wp,  0.9_wp,  3.5_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp,  30._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp,  50._wp, & ! closed broadleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 175.0_wp,  0.31_wp, 1._wp,  50._wp, & ! closed broadleaved deciduous forest
                   &   0.15_wp,  0.8_wp,  4.0_wp, 1.5_wp, 175.0_wp,  0.31_wp, 1._wp,  30._wp, & ! open broadleaved deciduous forest
                   &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 250.0_wp,  0.27_wp, 1._wp,  50._wp, & ! closed needleleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 250.0_wp,  0.33_wp, 1._wp,  50._wp, & ! open needleleaved deciduous forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 210.0_wp,  0.29_wp, 1._wp,  50._wp, & ! mixed broadleaved and needleleaved forest
                   &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.60_wp, 1._wp,  30._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                   &   0.20_wp,  0.8_wp,  2.5_wp, 1.0_wp, 150.0_wp,  0.65_wp, 1._wp,  10._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 1.5_wp, 225.0_wp,  0.65_wp, 1._wp,  50._wp, & ! closed to open shrubland
                   &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp, 100.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open herbaceous vegetation
                   &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp,  80.0_wp,  0.76_wp, 1._wp,  10._wp, & ! sparse vegetation
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed to open forest regulary flooded
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 150.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed forest or shrubland permanently flooded
                   &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  80.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open grassland regularly flooded
                   &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 180.0_wp,  0.50_wp, 1._wp, 200._wp, & ! artificial surfaces
                   &   0.05_wp,  0.05_wp, 0.6_wp, 0.3_wp, 200.0_wp,  0.76_wp, 1._wp, 200._wp, & ! bare areas
                   &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 150.0_wp,  -1.0_wp,-1._wp, 200._wp, & ! water bodies
                   &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 200._wp, & ! permanent snow and ice
                   &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp, 200._wp  / ! undefined

! Even more tuned version of gcv2009 by Guenther Zaengl (appears to produce the smallest temperature biases)
 DATA lu_gcv2009_v3 /  0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 190.0_wp,  0.72_wp, 1._wp,  30._wp, & ! irrigated croplands
                   &   0.07_wp,  0.9_wp,  3.3_wp, 1.0_wp, 170.0_wp,  0.72_wp, 1._wp,  30._wp, & ! rainfed croplands
                   &   0.25_wp,  0.8_wp,  3.0_wp, 0.5_wp, 160.0_wp,  0.55_wp, 1._wp,  10._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.07_wp,  0.9_wp,  3.5_wp, 0.7_wp, 150.0_wp,  0.72_wp, 1._wp,  30._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 280.0_wp,  0.38_wp, 1._wp,  50._wp, & ! closed broadleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  6.0_wp, 1.0_wp, 225.0_wp,  0.31_wp, 1._wp,  50._wp, & ! closed broadleaved deciduous forest
                   &   0.15_wp,  0.8_wp,  4.0_wp, 1.5_wp, 225.0_wp,  0.31_wp, 1._wp,  30._wp, & ! open broadleaved deciduous forest
                   &   1.00_wp,  0.8_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.27_wp, 1._wp,  50._wp, & ! closed needleleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.33_wp, 1._wp,  50._wp, & ! open needleleaved deciduous forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.8_wp, 270.0_wp,  0.29_wp, 1._wp,  50._wp, & ! mixed broadleaved and needleleaved forest
                   &   0.20_wp,  0.8_wp,  2.5_wp, 0.8_wp, 200.0_wp,  0.60_wp, 1._wp,  30._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                   &   0.20_wp,  0.8_wp,  2.5_wp, 0.6_wp, 200.0_wp,  0.65_wp, 1._wp,  10._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 0.9_wp, 265.0_wp,  0.65_wp, 1._wp,  50._wp, & ! closed to open shrubland
                   &   0.03_wp,  0.9_wp,  3.1_wp, 0.4_wp, 140.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open herbaceous vegetation
                   &   0.05_wp,  0.5_wp,  0.6_wp, 0.2_wp, 120.0_wp,  0.76_wp, 1._wp,  10._wp, & ! sparse vegetation
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed to open forest regulary flooded
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed forest or shrubland permanently flooded
                   &   0.05_wp,  0.8_wp,  2.0_wp, 0.7_wp, 120.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open grassland regularly flooded
                   &   1.00_wp,  0.2_wp,  1.6_wp, 0.2_wp, 300.0_wp,  0.50_wp, 1._wp, 200._wp, & ! artificial surfaces
                   &   0.05_wp,  0.05_wp, 0.6_wp,0.05_wp, 300.0_wp,  0.76_wp, 1._wp, 200._wp, & ! bare areas
                   &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 150.0_wp,  -1.0_wp,-1._wp, 200._wp, & ! water bodies
                   &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 200._wp, & ! permanent snow and ice
                   &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp, 200._wp  / ! undefined

! Yet another tuned version by Guenther Zaengl (adjusted to resistance-based bare soil evaporation scheme)
 DATA lu_gcv2009_v4 /  0.25_wp,  0.9_wp,  3.3_wp, 1.0_wp, 225.0_wp,  0.72_wp, 1._wp, 100._wp, & ! irrigated croplands
                   &   0.10_wp,  0.9_wp,  3.3_wp, 1.0_wp, 140.0_wp,  0.72_wp, 1._wp,  50._wp, & ! rainfed croplands
                   &   0.30_wp,  0.8_wp,  3.0_wp, 1.0_wp, 130.0_wp,  0.55_wp, 1._wp,  30._wp, & ! mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.10_wp,  0.9_wp,  3.5_wp, 1.0_wp, 120.0_wp,  0.72_wp, 1._wp,  40._wp, & ! mosaic vegetation (50-70%) - cropland (20-50%)
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 250.0_wp,  0.38_wp, 1._wp,  40._wp, & ! closed broadleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp,1.25_wp, 300.0_wp,  0.31_wp, 1._wp,  30._wp, & ! closed broadleaved deciduous forest
                   &   0.50_wp,  0.8_wp,  4.0_wp, 1.5_wp, 225.0_wp,  0.31_wp, 1._wp,  50._wp, & ! open broadleaved deciduous forest
                   &   1.00_wp,  0.8_wp,  5.0_wp,0.75_wp, 300.0_wp,  0.27_wp, 1._wp,  50._wp, & ! closed needleleaved evergreen forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 0.6_wp, 300.0_wp,  0.33_wp, 1._wp,  10._wp, & ! open needleleaved deciduous forest
                   &   1.00_wp,  0.9_wp,  5.0_wp, 1.0_wp, 270.0_wp,  0.29_wp, 1._wp,  15._wp, & ! mixed broadleaved and needleleaved forest
                   &   0.15_wp,  0.8_wp,  2.5_wp, 1.1_wp, 170.0_wp,  0.60_wp, 1._wp,  30._wp, & ! mosaic shrubland (50-70%) - grassland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 0.9_wp, 170.0_wp,  0.65_wp, 1._wp,  30._wp, & ! mosaic grassland (50-70%) - shrubland (20-50%)
                   &   0.15_wp,  0.8_wp,  2.5_wp, 1.5_wp, 180.0_wp,  0.65_wp, 1._wp,  75._wp, & ! closed to open shrubland
                   &   0.03_wp,  0.9_wp,  3.1_wp, 0.6_wp, 100.0_wp,  0.76_wp, 1._wp,  70._wp, & ! closed to open herbaceous vegetation
                   &   0.05_wp,  0.5_wp,  0.6_wp, 0.3_wp, 140.0_wp,  0.76_wp, 1._wp,  15._wp, & ! sparse vegetation
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp,  50._wp, & ! closed to open forest regulary flooded
                   &   1.00_wp,  0.8_wp,  5.0_wp, 1.0_wp, 190.0_wp,  0.30_wp, 1._wp,  80._wp, & ! closed forest or shrubland permanently flooded
                   &   0.05_wp,  0.8_wp,  2.0_wp, 1.0_wp,  80.0_wp,  0.76_wp, 1._wp,  30._wp, & ! closed to open grassland regularly flooded
                   &   1.00_wp,  0.2_wp,  1.6_wp, 0.6_wp, 300.0_wp,  0.50_wp, 1._wp, 200._wp, & ! artificial surfaces
                   &   0.02_wp,  0.01_wp, 0.2_wp, 0.3_wp, 300.0_wp,  0.76_wp, 1._wp, 200._wp, & ! bare areas
                   &   0.0002_wp,0.0_wp,  0.0_wp, 0.0_wp, 150.0_wp,  -1.0_wp,-1._wp, 200._wp, & ! water bodies
                   &   0.01_wp,  0.0_wp,  0.0_wp, 0.0_wp, 120.0_wp,  -1.0_wp, 1._wp, 200._wp, & ! permanent snow and ice
                   &   0.00_wp,  0.0_wp,  0.0_wp, 0.0_wp, 250.0_wp,  -1.0_wp,-1._wp, 200._wp  / ! undefined

    !----------------------------------------------------------------------

    IF(p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    ENDIF

    IF ( itopo == 1 .AND. iforcing /= inwp ) THEN

      DO jg = 1,n_dom

        ! Read topography and land-sea mask

        IF (n_dom > 1) THEN
          WRITE(land_sso_fn , '(a,i2.2,a)') 'bc_land_sso_DOM' , jg, '.nc'
          WRITE(land_frac_fn , '(a,i2.2,a)') 'bc_land_frac_DOM' , jg, '.nc'
        ELSE
          land_sso_fn  = 'bc_land_sso.nc'
          land_frac_fn = 'bc_land_frac.nc'
        ENDIF

        CALL openInputFile(stream_id, land_sso_fn, p_patch(jg), default_read_method)
        CALL read_2D(stream_id, on_cells, 'oromea', &
          &          ext_data(jg)%atm%topography_c)
        CALL closeFile(stream_id)

        CALL openInputFile(stream_id, land_frac_fn, p_patch(jg), default_read_method)
        CALL read_2D(stream_id, on_cells, 'notsea', &
          &          ext_data(jg)%atm%fr_land)
        CALL closeFile(stream_id)

      END DO

    END IF

    ! If ocean coupling is used, then read the land sea masks
    IF ( is_coupled_to_ocean() ) THEN
      DO jg = 1,n_dom
        IF ( iforcing == inwp ) THEN

          ! --- option NWP grids: see below in read_extdata

        ELSE

          ! --- option MPI grids for coupling: Read fraction of land (land-sea mask) from
          ! interpolated ocean grid (ocean: integer 0/1 lsm). lsm_ctr_c is the fraction of land.
          ! 0.0 is ocean, 1.0 is land, fractions on coasts (lakes are land).
          ! Used for lsmnolake in src/atm_coupling/mo_atmo_coupling_frame.f90.

          CALL openInputFile(stream_id, p_patch(jg)%grid_filename, p_patch(jg), default_read_method)

          CALL read_2D (stream_id, on_cells, 'cell_sea_land_mask', ext_data(jg)%atm%lsm_ctr_c)

          CALL closeFile(stream_id)

        END IF
      END DO
    END IF

    !------------------------------------------------!
    ! Read data from ExtPar file                     !
    !------------------------------------------------!

    IF ( (itopo == 1 .OR. l_scm_mode) .AND. iforcing == inwp ) THEN
      DO jg = 1,n_dom

        ! Preset parameter fields with the correct table values
        ilu = 0
        IF (i_lctype(jg) == GLC2000) THEN
          ext_data(jg)%atm%i_lc_snow_ice = 21
          ext_data(jg)%atm%i_lc_water    = 20
          ext_data(jg)%atm%i_lc_urban    = 22
          ext_data(jg)%atm%i_lc_shrub_eg = 11
          ext_data(jg)%atm%i_lc_shrub    = 12
          ext_data(jg)%atm%i_lc_grass    = 13
          ext_data(jg)%atm%i_lc_bare_soil= 19
          ext_data(jg)%atm%i_lc_sparse   = 14
          DO i = 1, num_lcc*n_param_lcc, n_param_lcc
            ilu=ilu+1
            ext_data(jg)%atm%z0_lcc(ilu)          = lu_glc2000(i  )  ! Land-cover related roughness length
            ext_data(jg)%atm%plcovmax_lcc(ilu)    = lu_glc2000(i+1)  ! Maximum plant cover fraction for each land-cover class
            ext_data(jg)%atm%laimax_lcc(ilu)      = lu_glc2000(i+2)  ! Maximum leaf area index for each land-cover class
            ext_data(jg)%atm%rootdmax_lcc(ilu)    = lu_glc2000(i+3)  ! Maximum root depth for each land-cover class
            ext_data(jg)%atm%stomresmin_lcc(ilu)  = lu_glc2000(i+4)  ! Minimum stomata resistance for each land-cover class
            ext_data(jg)%atm%snowalb_lcc(ilu)     = lu_glc2000(i+5)  ! Albedo in case of snow cover for each land-cover class
            ext_data(jg)%atm%snowtile_lcc(ilu)    = &
              &          MERGE(.TRUE.,.FALSE.,lu_glc2000(i+6)>0._wp) ! Existence of snow tiles for land-cover class
          ENDDO
        ELSE IF (i_lctype(jg) == GLOBCOVER2009) THEN
          SELECT CASE (itype_lndtbl)
          CASE (1)
            lu_gcv => lu_gcv2009
          CASE (2)
            lu_gcv => lu_gcv2009_v2
          CASE (3)
            lu_gcv => lu_gcv2009_v3
          CASE (4)
            lu_gcv => lu_gcv2009_v4
          END SELECT

          ext_data(jg)%atm%i_lc_snow_ice    = 22
          ext_data(jg)%atm%i_lc_water       = 21
          ext_data(jg)%atm%i_lc_urban       = 19
          ext_data(jg)%atm%i_lc_shrub_eg    = 12
          ext_data(jg)%atm%i_lc_shrub       = 13
          ext_data(jg)%atm%i_lc_grass       = 14
          ext_data(jg)%atm%i_lc_bare_soil   = 20
          ext_data(jg)%atm%i_lc_sparse      = 15
          ext_data(jg)%atm%i_lc_crop_irrig  = 1
          ext_data(jg)%atm%i_lc_crop_rain   = 2
          ext_data(jg)%atm%i_lc_crop_mos    = 3
          ext_data(jg)%atm%i_lc_veg_mos     = 4 
          ext_data(jg)%atm%i_lc_forest_b_eg = 5
          ext_data(jg)%atm%i_lc_forest_b_d  = 6
          ext_data(jg)%atm%i_lc_woodland    = 7
          ext_data(jg)%atm%i_lc_forest_n_eg = 8
          ext_data(jg)%atm%i_lc_forest_n_d  = 9
          ext_data(jg)%atm%i_lc_forest_bn   = 10
          ext_data(jg)%atm%i_lc_shrub_mos   = 11
          ext_data(jg)%atm%i_lc_forest_rf   = 16
          ext_data(jg)%atm%i_lc_forest_pf   = 17
          ext_data(jg)%atm%i_lc_grass_rf    = 18
          
          DO i = 1, num_lcc*n_param_lcc, n_param_lcc
            ilu=ilu+1
            ext_data(jg)%atm%z0_lcc(ilu)          = lu_gcv(i  )  ! Land-cover related roughness length
            ext_data(jg)%atm%plcovmax_lcc(ilu)    = lu_gcv(i+1)  ! Maximum plant cover fraction for each land-cover class
            ext_data(jg)%atm%laimax_lcc(ilu)      = lu_gcv(i+2)  ! Maximum leaf area index for each land-cover class
            ext_data(jg)%atm%rootdmax_lcc(ilu)    = lu_gcv(i+3)  ! Maximum root depth for each land-cover class
            ext_data(jg)%atm%stomresmin_lcc(ilu)  = lu_gcv(i+4)  ! Minimum stomata resistance for each land-cover class
            ext_data(jg)%atm%snowalb_lcc(ilu)     = lu_gcv(i+5)  ! Albedo in case of snow cover for each land-cover class
            ext_data(jg)%atm%snowtile_lcc(ilu)    = &
              &          MERGE(.TRUE.,.FALSE.,lu_gcv(i+6)>0._wp) ! Existence of snow tiles for land-cover class
            IF (cskinc <= 0._wp) THEN
              ext_data(jg)%atm%skinc_lcc(ilu)     = lu_gcv(i+7)  ! Skin conductivity for each land use class
            ELSE
              ext_data(jg)%atm%skinc_lcc(ilu)     = cskinc       ! Constant value specified in namelist
            ENDIF
          ENDDO
        ENDIF

        !$ACC UPDATE &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_snow_ice) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_water) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_urban) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_shrub_eg) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_shrub) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_grass) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_bare_soil) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_sparse) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_crop_irrig) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_crop_rain) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_crop_mos) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_veg_mos) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_b_eg) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_b_d) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_woodland) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_n_eg) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_n_d) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_bn) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_shrub_mos) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_rf) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_forest_pf) &
        !$ACC   DEVICE(ext_data(jg)%atm%i_lc_grass_rf) &
        !$ACC   ASYNC(1)

        ! Urban canopy parameters
        DO ilu = 1, num_lcc
          IF (ilu == ext_data(jg)%atm%i_lc_urban) THEN
            ext_data(jg)%atm%fr_paved_lcc(ilu) = 1._wp        ! Impervious surface area (ISA) fraction for urban land use class
            ext_data(jg)%atm%ahf_lcc(ilu)      = 15._wp       ! Anthropogenic heat flux for urban land use class
          ELSE
            ext_data(jg)%atm%fr_paved_lcc(ilu) = 0._wp
            ext_data(jg)%atm%ahf_lcc(ilu)      = 0._wp
          ENDIF
        ENDDO

        ! Derived parameter: minimum allowed land-cover related roughness length in the
        ! presence of low ndvi and/or snow cover
        DO ilu = 1, num_lcc
          IF (ilu == ext_data(jg)%atm%i_lc_urban .OR. ilu == ext_data(jg)%atm%i_lc_water) THEN
            ext_data(jg)%atm%z0_lcc_min(ilu) = ext_data(jg)%atm%z0_lcc(ilu) ! no reduction in urban regions and over water
          ELSE IF (pp_sso == 2 .AND. ext_data(jg)%atm%z0_lcc(ilu) >= 0.5_wp) THEN   ! if MERIT/REMA orography is used:
            ext_data(jg)%atm%z0_lcc_min(ilu) = 0.75_wp*ext_data(jg)%atm%z0_lcc(ilu) ! 75% for nominal roughness lengths >= 50 cm
          ELSE IF (ext_data(jg)%atm%z0_lcc(ilu) > 0.1_wp) THEN
            ext_data(jg)%atm%z0_lcc_min(ilu) = 0.3_wp*ext_data(jg)%atm%z0_lcc(ilu) ! 30% for nominal roughness lengths > 10 cm
          ELSE
            ext_data(jg)%atm%z0_lcc_min(ilu) = 0.1_wp*ext_data(jg)%atm%z0_lcc(ilu) ! 10% otherwise
          ENDIF
        ENDDO

      ENDDO ! jg
    END IF

    IF (itopo == 1 .AND. iforcing == inwp) THEN
      DO jg = 1,n_dom
        ! Start reading external parameter data
        ! The cdi-based read routines are used for GRIB2 input data only due to performance problems
        IF (read_netcdf_parallel) THEN
          extpar_file = generate_filename(extpar_filename, getModelBaseDir(), &
            &                             TRIM(p_patch(jg)%grid_filename),    &
            &                              nroot,                             &
            &                             p_patch(jg)%level, p_patch(jg)%id)
          CALL openinputfile(stream_id, extpar_file, p_patch(jg), default_read_method)
        ELSE
          parameters = makeInputParameters(cdi_extpar_id(jg), p_patch(jg)%n_patch_cells_g, p_patch(jg)%comm_pat_scatter_c, &
          &                                opt_dict=extpar_varnames_dict)
        ENDIF

        !--------------------------------------------------------------------
        !
        ! Read topography for triangle centers (triangular grid)
        !
        !--------------------------------------------------------------------
        CALL read_extdata('topography_c', ext_data(jg)%atm%topography_c)

        ! If ocean coupling and TERRA is used, then read the land sea masks

        IF ( is_coupled_to_ocean() .AND. atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN

          ! --- option NWP grids for coupling: Read fraction of land (land-sea mask) from
          ! interpolated ocean grid (ocean: integer 0/1 lsm). lsm_ctr_c is the fraction of land.
          ! 0.0 is ocean, 1.0 is land, fractions on coasts (lakes are land).
          ! Used in routine lsm_ocean_atmo.

          CALL read_extdata('cell_sea_land_mask', ext_data(jg)%atm%lsm_ctr_c)

        ENDIF

        !
        ! other external parameters on triangular grid
        !
        CALL read_extdata('FR_LAND', ext_data(jg)%atm%fr_land)

        SELECT CASE ( iforcing )
        CASE ( inwp )
          CALL read_extdata('NDVI_MAX',  ext_data(jg)%atm%ndvi_max)
          CALL read_extdata('SOILTYP',   arr2di=ext_data(jg)%atm%soiltyp)
          CALL read_extdata('T_CL',      ext_data(jg)%atm%t_cl)
          CALL read_extdata('SSO_STDH',  ext_data(jg)%atm%sso_stdh)
          CALL read_extdata('SSO_THETA', ext_data(jg)%atm%sso_theta)
          CALL read_extdata('SSO_GAMMA', ext_data(jg)%atm%sso_gamma)
          CALL read_extdata('SSO_SIGMA', ext_data(jg)%atm%sso_sigma)
          CALL read_extdata('FR_LAKE',   ext_data(jg)%atm%fr_lake)
          CALL read_extdata('DEPTH_LK',  ext_data(jg)%atm%depth_lk)

          IF (islope_rad(jg) >= 2) THEN
            CALL read_extdata('HORIZON', arr3d=ext_data(jg)%atm%horizon,ltime=.FALSE.)
            CALL read_extdata('SKYVIEW', ext_data(jg)%atm%skyview)
  
            rl_start = 1
            rl_end   = min_rlcell
            i_startblk = p_patch(jg)%cells%start_block(rl_start)
            i_endblk   = p_patch(jg)%cells%end_block(rl_end)
  
            ! Test consistency of horizon
!$OMP DO PRIVATE(jb,jc,im,i_startidx,i_endidx)
            DO jb = i_startblk, i_endblk
              CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
              DO im = 1, nhori
                DO jc = i_startidx,i_endidx
                  IF (ext_data(jg)%atm%horizon(jc,jb,im) > 90.0_wp .OR. &
                      ext_data(jg)%atm%horizon(jc,jb,im) < 0.0_wp) THEN
                    WRITE (message_text,'(A,F8.3)') 'ext_data(jg)%atm%horizon(jc,jb,im) = ',&
                                                     ext_data(jg)%atm%horizon(jc,jb,im)
                    CALL message(routine, message_text)
                    WRITE (message_text,'(A,3I8)') '   DEBUG point: ', jc, jb, im
                    CALL message(routine, message_text)
                    CALL finish(routine,'HORIZON is out of bounds!')
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
!$OMP END DO
          ENDIF
  
          CALL read_extdata('LU_CLASS_FRACTION', arr3d=ext_data(jg)%atm%lu_class_fraction,ltime=.FALSE.) 

          ! The following fields are only required without surface tiles
          IF (ntiles_lnd == 1) THEN
            CALL read_extdata('PLCOV_MX', ext_data(jg)%atm%plcov_mx)
            CALL read_extdata('LAI_MX',   ext_data(jg)%atm%lai_mx)
            CALL read_extdata('ROOTDP',   ext_data(jg)%atm%rootdp)
            CALL read_extdata('RSMIN',    ext_data(jg)%atm%rsmin)
            CALL read_extdata('FOR_D',    ext_data(jg)%atm%for_d)
            CALL read_extdata('FOR_E',    ext_data(jg)%atm%for_e)
          ENDIF

          IF (atm_phy_nwp_config(jg)%itype_z0 == 1) THEN
            ! only read, if contribution from sub-scale orography should be included in z0
            CALL read_extdata('Z0', ext_data(jg)%atm%z0)
          ENDIF

          IF (is_frglac_in(jg)) THEN
            ! for backward compatibility with extpar files generated prior to 2014-01-31
             CALL read_extdata('ICE', ext_data(jg)%atm%fr_glac)
          ELSE
            ! for new extpar files (generated after 2014-01-31)
            ! take it from lu_class_fraction
            ext_data(jg)%atm%fr_glac(:,:) = ext_data(jg)%atm%lu_class_fraction(:,:,ext_data(jg)%atm%i_lc_snow_ice)
          ENDIF

          IF (itype_lwemiss == 2) THEN
            CALL read_extdata('EMISS',   arr3d=ext_data(jg)%atm_td%lw_emiss)
          ELSE IF (itype_lwemiss == 1) THEN
            CALL read_extdata('EMIS_RAD', ext_data(jg)%atm%emis_rad)
          ELSE
            ext_data(jg)%atm%emis_rad(:,:)= zemiss_def
          ENDIF

          ! Copy sso_stdh to sso_stdh_raw before applying correction for orography filtering
          ext_data(jg)%atm%sso_stdh_raw(:,:) = ext_data(jg)%atm%sso_stdh(:,:)


          IF ( iprog_aero > 1) THEN
            CALL read_extdata('emi_bc',  arr2d=ext_data(jg)%atm%emi_bc )
            CALL read_extdata('emi_oc',  arr2d=ext_data(jg)%atm%emi_oc )
            CALL read_extdata('emi_so2', arr2d=ext_data(jg)%atm%emi_so2)
          ENDIF
          ! Read time dependent data
          IF ( irad_aero == iRadAeroTegen .OR. irad_aero == iRadAeroART) THEN
            CALL read_extdata('AER_SS',   arr3d=ext_data(jg)%atm_td%aer_ss)
            CALL read_extdata('AER_DUST', arr3d=ext_data(jg)%atm_td%aer_dust)
            CALL read_extdata('AER_ORG',  arr3d=ext_data(jg)%atm_td%aer_org)
            CALL read_extdata('AER_SO4',  arr3d=ext_data(jg)%atm_td%aer_so4)
            CALL read_extdata('AER_BC',   arr3d=ext_data(jg)%atm_td%aer_bc)
          ENDIF  ! irad_aero
          CALL read_extdata('NDVI_MRAT', arr3d=ext_data(jg)%atm_td%ndvi_mrat)

          IF (sstice_mode == SSTICE_ANA_CLINC) THEN
            CALL read_extdata('T_SEA', arr3d=ext_data(jg)%atm_td%sst_m)
          ENDIF

          IF (itype_vegetation_cycle > 1) THEN
            CALL read_extdata('T_2M_CLIM', arr3d=ext_data(jg)%atm_td%t2m_m)
            CALL read_extdata('TOPO_CLIM',   ext_data(jg)%atm%topo_t2mclim)
          ENDIF

          !--------------------------------
          ! If MODIS albedo is used
          !--------------------------------
          IF ( albedo_type == MODIS) THEN
            CALL read_extdata('ALB',   arr3d=ext_data(jg)%atm_td%alb_dif)
            CALL read_extdata('ALUVD', arr3d=ext_data(jg)%atm_td%albuv_dif)
            CALL read_extdata('ALNID', arr3d=ext_data(jg)%atm_td%albni_dif)

            rl_start = 1
            rl_end   = min_rlcell

            i_startblk = p_patch(jg)%cells%start_block(rl_start)
            i_endblk   = p_patch(jg)%cells%end_block(rl_end)

            albthresh = 0.3_wp ! threshold value for albedo modification

!$OMP PARALLEL
            ! Scale from [%] to [1]
            CALL var_scale(ext_data(jg)%atm_td%alb_dif(:,:,:), 1._wp/100._wp)
            CALL var_scale(ext_data(jg)%atm_td%albuv_dif(:,:,:), 1._wp/100._wp)
            CALL var_scale(ext_data(jg)%atm_td%albni_dif(:,:,:), 1._wp/100._wp)
!$OMP BARRIER


            IF (itune_albedo >= 1) THEN
              ! Test: reduce albedo over land where modis albedo is higher than 0.3 (variable albthresh)
!$OMP DO PRIVATE(jb,jc,im,i_startidx,i_endidx,albfac)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

                DO im = 1, 12
                  DO jc = i_startidx,i_endidx
                    IF (ext_data(jg)%atm%soiltyp(jc,jb) >= 2 .AND. ext_data(jg)%atm%soiltyp(jc,jb) <= 8) THEN
                      IF (ext_data(jg)%atm_td%alb_dif(jc,jb,im) > albthresh) THEN
                        albfac = (albthresh+2._wp*ext_data(jg)%atm_td%alb_dif(jc,jb,im))/ &
                          (3._wp*ext_data(jg)%atm_td%alb_dif(jc,jb,im))
                        ext_data(jg)%atm_td%alb_dif(jc,jb,im)   = albfac*ext_data(jg)%atm_td%alb_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albuv_dif(jc,jb,im) = albfac*ext_data(jg)%atm_td%albuv_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albni_dif(jc,jb,im) = albfac*ext_data(jg)%atm_td%albni_dif(jc,jb,im)
                      ENDIF
                    ENDIF
                  ENDDO
               ENDDO
              ENDDO
!$OMP END DO
            ENDIF  ! Sahara albedo tuning

            IF (itune_albedo >= 2) THEN
              ! Increase albedo over the Antarctic plateau by 5% (from 70% to 75%) in order to get rid of summertime warm bias
!$OMP DO PRIVATE(jb,jc,im,i_startidx,i_endidx,albfac)
              DO jb = i_startblk, i_endblk
                CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

                DO im = 1, 12
                  DO jc = i_startidx,i_endidx
                    IF (ext_data(jg)%atm%soiltyp(jc,jb) == 1 .AND. p_patch(jg)%cells%center(jc,jb)%lat*rad2deg < -65._wp ) THEN
                      IF (ext_data(jg)%atm%topography_c(jc,jb) > 1000._wp) THEN
                        albfac = MIN(1._wp,1.e-3_wp*(ext_data(jg)%atm%topography_c(jc,jb)-1000._wp))
                        ext_data(jg)%atm_td%alb_dif(jc,jb,im)   = 0.05_wp*albfac + ext_data(jg)%atm_td%alb_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albuv_dif(jc,jb,im) = 0.05_wp*albfac + ext_data(jg)%atm_td%albuv_dif(jc,jb,im)
                        ext_data(jg)%atm_td%albni_dif(jc,jb,im) = 0.05_wp*albfac + ext_data(jg)%atm_td%albni_dif(jc,jb,im)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
!$OMP END DO
            ENDIF  ! Antarctic albedo tuning
!$OMP END PARALLEL

          END IF  !  albedo_type

        END SELECT ! iforcing

        IF (read_netcdf_parallel) THEN
          CALL closeFile(stream_id)
        ELSE
          CALL deleteInputParameters(parameters)
        ENDIF

        !
        ! derived external parameter fields
        !

        ! land sea mask at cell centers (LOGICAL)
        !

        ! adjust atmo LSM to ocean LSM for coupled simulation and initialize new land points (TERRA only)

        IF ( is_coupled_to_ocean() .AND. atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN
          CALL lsm_ocean_atmo ( p_patch(jg), ext_data(jg) )
        ENDIF

        i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

        rl_start = 1
        rl_end   = min_rlcell

        i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
        i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

        DO jb = i_startblk, i_endblk
          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            &                i_startidx, i_endidx, rl_start, rl_end)

          ! Loop starts with 1 instead of i_startidx
          ! because the start index is missing in RRTM
          DO jc = 1,i_endidx
            IF (ext_data(jg)%atm%fr_land(jc,jb) > 0.5_wp) THEN
              ext_data(jg)%atm%llsm_atm_c(jc,jb) = .TRUE.  ! land point
            ELSE
              ext_data(jg)%atm%llsm_atm_c(jc,jb) = .FALSE. ! water point
            ENDIF
            IF (ext_data(jg)%atm%fr_lake(jc,jb) >= 0.5_wp) THEN
              ext_data(jg)%atm%llake_c(jc,jb) = .TRUE.     ! lake point
            ELSE
              ext_data(jg)%atm%llake_c(jc,jb) = .FALSE.    ! no lake point
            ENDIF
          ENDDO
        ENDDO


      ENDDO  ! jg

    ENDIF ! (itopo == 1 and inwp)

    !-------------------------------------------------------
    ! Read ozone
    !-------------------------------------------------------

    IF ( iforcing == inwp .AND. (irad_o3 == io3_clim .OR. irad_o3 == io3_ape) ) THEN

      DO jg = 1,n_dom

        WRITE(ozone_file,'(a,I2.2,a)') 'o3_icon_DOM',jg,'.nc'

        IF(my_process_is_stdio()) THEN
          ! open file
          !
          CALL nf(nf_open(TRIM(ozone_file), NF_NOWRITE, ncid), routine)
          WRITE(0,*)'read ozone levels'
          CALL nf(nf_inq_varid(ncid, TRIM(levelname), varid), routine)
          CALL nf(nf_get_var_double(ncid, varid, zdummy_o3lev(:)), routine)
          CALL nf(nf_close(ncid), routine)
          !
        ENDIF ! pe

        CALL p_bcast(zdummy_o3lev(:), p_io, mpi_comm)

        DO jk=1,nlev_o3
          ext_data(jg)%atm_td%pfoz(jk)=zdummy_o3lev(jk)
        ENDDO

        ! define half levels of ozone pressure grid
        ! upper boundary: ph =      0.Pa -> extrapolation of uppermost value
        ! lower boundary: ph = 125000.Pa -> extrapolation of lowermost value
        ext_data(jg)%atm_td%phoz(1)           = 0._wp
        ext_data(jg)%atm_td%phoz(2:nlev_o3) = (ext_data(jg)%atm_td%pfoz(1:nlev_o3-1) &
          &                                   +  ext_data(jg)%atm_td%pfoz(2:nlev_o3))*.5_wp
        ext_data(jg)%atm_td%phoz(nlev_o3+1) = 125000._wp

        DO i=1,nlev_o3

          WRITE(message_text,'(a,i4,f12.4,f12.4)')'full/half level press ozone ', &
                              i, ext_data(jg)%atm_td%pfoz(i), ext_data(jg)%atm_td%phoz(i+1)
          CALL message(routine, TRIM(message_text))
        ENDDO

        CALL openinputfile(stream_id, ozone_file, p_patch(jg), default_read_method)

        CALL read_3D_extdim(stream_id, on_cells, TRIM(o3name), &
          &                 ext_data(jg)%atm_td%O3)

        ! convert from ppmv to g/g only in case of APE ozone
        ! whether o3mr2gg or ppmv2gg is used to convert O3 to gg depends on the units of
        ! the incoming ozone file.  Often, the incoming units are not ppmv.
        !
        IF(iforcing == inwp .AND. irad_o3 == io3_ape) THEN
           ! ozone input expected in units of ppmv
           WRITE(message_text,'(a,f12.4,f12.4)')'MAX/MIN o3 ppmv', &
                & MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)), MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))
           CALL message(routine, TRIM(message_text))
           ext_data(jg)%atm_td%O3(:,:,:,:)= ext_data(jg)%atm_td%O3(:,:,:,:)*ppmv2gg
        END IF

        WRITE(message_text,'(a,e12.4,e12.4)')'MAX/MIN o3 g/g', &
           MAXVAL(ext_data(jg)%atm_td%O3(:,:,:,:)), MINVAL(ext_data(jg)%atm_td%O3(:,:,:,:))
        CALL message(routine, TRIM(message_text))

        ! close file
        CALL closeFile(stream_id)

      ENDDO ! ndom
    END IF ! irad_o3 and inwp


    !------------------------------------------
    ! Read time dependent SST and ICE Fraction
    !------------------------------------------
    IF (sstice_mode == SSTICE_CLIM .AND. iforcing == inwp) THEN

      DO jg = 1,n_dom
       !Read the climatological values for SST and ice cover

        DO im=1,12

         sst_td_file= generate_td_filename(sst_td_filename,                &
           &                             getModelBaseDir(),                &
           &                             TRIM(p_patch(jg)%grid_filename),  &
           &                             im,clim=.TRUE.                   )

         IF(is_mpi_workroot) THEN

          CALL message  (routine, TRIM(sst_td_file))

          INQUIRE (FILE=sst_td_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(routine,'td sst external data file is not found.')
          ENDIF

         ENDIF

         CALL openinputfile(stream_id, sst_td_file, p_patch(jg), &
          &                        default_read_method)
         CALL read_2D (stream_id, on_cells, 'SST', &
          &            ext_data(jg)%atm_td%sst_m(:,:,im))
         CALL closeFile(stream_id)

         ci_td_file= generate_td_filename(ci_td_filename,                  &
           &                             getModelBaseDir(),                &
           &                             TRIM(p_patch(jg)%grid_filename),  &
           &                             im,clim=.TRUE.                   )

         IF(my_process_is_stdio()) THEN

          CALL message  (routine, TRIM(ci_td_file))

          INQUIRE (FILE=ci_td_file, EXIST=l_exist)
          IF (.NOT.l_exist) THEN
            CALL finish(routine,'td ci external data file is not found.')
          ENDIF

         ENDIF

         CALL openinputfile(stream_id, ci_td_file, p_patch(jg), default_read_method)
         CALL read_2D(stream_id, on_cells, 'CI', &
          &           ext_data(jg)%atm_td%fr_ice_m(:,:,im))
         CALL closeFile(stream_id)

        END DO

      END DO ! ndom

   END IF ! sstice_mode and inwp

   CONTAINS

     ! Wrapper routine for reading input data via cdilib for GRIB2 or via optimized netcdf routines
     !
     SUBROUTINE read_extdata(varname,arr2d,arr2di,arr3d,ltime)
       CHARACTER(LEN=*), INTENT(IN) :: varname  ! name of input variable
       REAL(wp), OPTIONAL, INTENT(INOUT) :: arr2d(:,:), arr3d(:,:,:) ! alternative I/O arrays
       INTEGER, OPTIONAL, INTENT(INOUT) :: arr2di(:,:)
       LOGICAL, OPTIONAL, INTENT(IN) :: ltime ! .true. if third dimension is time

       LOGICAL :: dim3_is_time

       IF (PRESENT(arr3d) .AND. PRESENT(ltime)) THEN
         dim3_is_time = ltime
       ELSE
         dim3_is_time = .TRUE.
       ENDIF

       IF (PRESENT(arr2d)) THEN
         IF (read_netcdf_parallel) THEN
           CALL read_2D(stream_id, on_cells, TRIM(varname), arr2d)
         ELSE
           CALL read_cdi_2d(parameters, TRIM(varname), arr2d)
         ENDIF
       ELSE IF (PRESENT(arr2di)) THEN
         IF (read_netcdf_parallel) THEN
           CALL read_2D_int(stream_id, on_cells, TRIM(varname), arr2di)
         ELSE
           CALL read_cdi_2d(parameters, TRIM(varname), arr2di)
         ENDIF
       ELSE IF (PRESENT(arr3d)) THEN
         IF (read_netcdf_parallel) THEN
           CALL read_2D_extdim(stream_id, on_cells, TRIM(varname), arr3d)
         ELSE IF (dim3_is_time) THEN
           CALL read_cdi_2d(parameters, SIZE(arr3d,3), TRIM(varname), arr3d)
         ELSE
           CALL read_cdi_3d(parameters, TRIM(varname), SIZE(arr3d,3), arr3d, opt_lev_dim=3 )
         ENDIF
       ENDIF

     END SUBROUTINE read_extdata

  END SUBROUTINE read_ext_data_atm
  !-------------------------------------------------------------------------



  SUBROUTINE init_index_lists (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch(:)
    TYPE(t_external_data), INTENT(INOUT) :: ext_data(:)

    INTEGER :: i_lu, jb,jc, jg, i_count, i_count_flk, ic, jt, jt_in
    INTEGER :: i_count_sea             ! number of sea points
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !< slices
    INTEGER :: i_nchdom                !< domain index
    LOGICAL  :: tile_mask(num_lcc), lhave_urban
    REAL(wp) :: tile_frac(num_lcc), sum_frac, dtdz_clim, t2mclim_hc, lat
    INTEGER  :: lu_subs, it_count(ntiles_total)
    INTEGER  :: npoints, npoints_sea, npoints_lake
    INTEGER  :: i_lc_water
    INTEGER, ALLOCATABLE :: icount_falseglac(:)

    REAL(wp), POINTER  ::  &  !< pointer to proportion of actual value/maximum
      &  ptr_ndviratio(:,:)   !< NDVI (for starting time of model integration)

    REAL(wp) :: scalfac       ! scaling factor for inflating dominating land fractions
                              ! to fr_land (or fr_land + fr_lake)
    REAL(wp) :: zfr_land      ! fr_land derived from land tile fractions

!!$    CHARACTER(len=max_char_length), PARAMETER :: &
!!$      routine = modname//':init_index_lists'

    !-------------------------------------------------------------------------

    WRITE(message_text,'(a,i4)') &
      &  'Index list generation - number of land tiles: ', ntiles_lnd
    CALL message('', TRIM(message_text))
    WRITE(message_text,'(a,i4)')  'Total number of tiles: ', ntiles_total
    CALL message('', TRIM(message_text))

    ! climatological temperature gradient used for height correction of T2M climatology
    dtdz_clim = -5.e-3_wp  ! -5 K/km

    DO jg = 1, n_dom

       ptr_ndviratio => ext_data(jg)%atm%ndviratio(:,:)

       i_nchdom  = MAX(1,p_patch(jg)%n_childdom)

       i_lc_water = ext_data(jg)%atm%i_lc_water

       ! Initialization of index list counts - moved here in order to avoid uninitialized elements
       ! along nest boundaries

       ext_data(jg)%atm%gp_count_t(:,:) = 0
       ext_data(jg)%atm%lp_count_t(:,:) = 0

       ALLOCATE(icount_falseglac(p_patch(jg)%nblks_c))
       icount_falseglac(:) = 0

!$OMP PARALLEL PRIVATE(rl_start,rl_end,i_startblk,i_endblk)
       !
       ! exclude the boundary interpolation zone of nested domains
       rl_start = grf_bdywidth_c+1
       rl_end   = min_rlcell_int

       i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
       i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)
#ifdef __SX__
! turn off OpenMP on the NEC until MAXLOC bug (not threadsafe) is fixed
!$OMP SINGLE
#else
!$OMP DO PRIVATE(jb,jc,i_lu,i_startidx,i_endidx,i_count,i_count_sea,i_count_flk,tile_frac,lhave_urban,&
!$OMP            tile_mask,lu_subs,sum_frac,scalfac,zfr_land,it_count,ic,jt,jt_in,t2mclim_hc,lat ) ICON_OMP_DEFAULT_SCHEDULE
#endif
       DO jb=i_startblk, i_endblk

         CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

         i_count                       = 0   ! counter for land points
         i_count_sea                   = 0   ! counter for sea points
         i_count_flk                   = 0   ! counter for lake points

         it_count(:)                   = 0 ! counter for tiles

         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%lc_class_t(jc,jb,:) = -1    ! dummy value for undefined points

           IF (ext_data(jg)%atm%fr_land(jc,jb)> frlnd_thrhld) THEN ! searching for land-points
             i_count=i_count+1
             ext_data(jg)%atm%list_land%idx(i_count,jb) = jc  ! write index of land-points

             tile_frac(:)= ext_data(jg)%atm%lu_class_fraction(jc,jb,:)
             tile_mask(:)=.true.
             tile_mask(i_lc_water)=.false. ! exclude water points

             lhave_urban = .false.

             ext_data(jg)%atm%list_land%ncount(jb) = i_count

             IF (ntiles_lnd == 1) THEN

               ! i_lu=1 contains grid-box mean values from EXTPAR!
               !
               ext_data(jg)%atm%lc_frac_t(jc,jb,1)  = 1._wp
               ext_data(jg)%atm%lc_class_t(jc,jb,1) = MAXLOC(tile_frac,1,tile_mask)
               lhave_urban = ext_data(jg)%atm%lc_class_t(jc,jb,1) == ext_data(jg)%atm%i_lc_urban
               !
               ! Urban Canopy Parameters (UCPs)
               !
               ! impervious surface area fraction of the urban canopy
               ext_data(jg)%atm%urb_isa_t(jc,jb,1)       = ext_data(jg)%atm%fr_paved_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))

               IF (lterra_urb) THEN
                 ! impervious surface area (ISA) fraction
                 ext_data(jg)%atm%fr_paved_t(jc,jb,1)    = ext_data(jg)%atm%fr_paved_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))

                 ! building area fraction with respect to urban tile
                 ext_data(jg)%atm%urb_fr_bld_t(jc,jb,1)  = 0.667_wp

                 ! street canyon H/W ratio
                 ext_data(jg)%atm%urb_h2w_t(jc,jb,1)     = 1.5_wp

                 ! surface area index of the urban canopy
                 ext_data(jg)%atm%urb_ai_t(jc,jb,1)      = (1.0_wp + 2.0_wp * ext_data(jg)%atm%urb_h2w_t(jc,jb,1))  &
                   &                                     * (1.0_wp - ext_data(jg)%atm%urb_fr_bld_t(jc,jb,1))        &
                   &                                     + ext_data(jg)%atm%urb_fr_bld_t(jc,jb,1)

                 ! albedo reduction factor for the urban canopy
                 ! Reduce the effective albedo according to the building density,
                 ! the reduction factor is based on Monte-Carlo simulations.
                 ext_data(jg)%atm%urb_alb_red_t(jc,jb,1) = EXP(-0.6_wp * ext_data(jg)%atm%urb_h2w_t(jc,jb,1))       &
                   &                                     * (1.0_wp - ext_data(jg)%atm%urb_fr_bld_t(jc,jb,1))        &
                   &                                     + ext_data(jg)%atm%urb_fr_bld_t(jc,jb,1)

                 ! building height
                 ext_data(jg)%atm%urb_h_bld_t(jc,jb,1)   = 15._wp

                 ! thermal albedo of urban material
                 ext_data(jg)%atm%urb_alb_th_t(jc,jb,1)  = 0.14_wp

                 ! solar albedo of urban material, times albedo reduction factor for the urban canopy
                 ext_data(jg)%atm%urb_alb_so_t(jc,jb,1)  = 0.101_wp * ext_data(jg)%atm%urb_alb_red_t(jc,jb,1)

                 ! volumetric heat capacity of urban material
                 ext_data(jg)%atm%urb_hcap_t(jc,jb,1)    = 1250000._wp

                 ! thermal conductivity of urban material
                 ext_data(jg)%atm%urb_hcon_t(jc,jb,1)    = 0.767_wp

                 ! anthropogenic heat flux
                 ext_data(jg)%atm%ahf_t(jc,jb,1)         = ext_data(jg)%atm%ahf_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))
               ENDIF
               !
               ! Vegetation Parameters etc.
               !
               ! root depth
               ext_data(jg)%atm%rootdp_t (jc,jb,1)  = ext_data(jg)%atm%rootdp(jc,jb)

               ! plant cover
               ext_data(jg)%atm%plcov_t  (jc,jb,1)  = ptr_ndviratio(jc,jb)                                                   &
                 &     * MIN(ext_data(jg)%atm%ndvi_max(jc,jb),ext_data(jg)%atm%plcov_mx(jc,jb))

               ! transpiration area index
               ext_data(jg)%atm%tai_t    (jc,jb,1)  = ext_data(jg)%atm%plcov_t(jc,jb,1)                                      &
                 &                                  * ext_data(jg)%atm%lai_mx(jc,jb)

               ! surface area index
               IF (lterra_urb) THEN
                 ext_data(jg)%atm%sai_t  (jc,jb,1)  = c_lnd * (1.0_wp - ext_data(jg)%atm%urb_isa_t(jc,jb,1))                 &
                                                    + ext_data(jg)%atm%urb_ai_t(jc,jb,1)*ext_data(jg)%atm%urb_isa_t(jc,jb,1) &
                                                    + ext_data(jg)%atm%tai_t(jc,jb,1)
               ELSE
                 ext_data(jg)%atm%sai_t  (jc,jb,1)  = c_lnd + ext_data(jg)%atm%tai_t(jc,jb,1)
               END IF

               ! evaporative soil area index
               IF (icpl_da_sfcevap >= 4 .OR. itype_evsl == 5) THEN
                 ext_data(jg)%atm%eai_t(jc,jb,1)    = MERGE(0.75_wp,                                                         &
                   2.0_wp-1.25_wp*tile_frac(ext_data(jg)%atm%i_lc_urban),lhave_urban)
                 ext_data(jg)%atm%r_bsmin(jc,jb)    = cr_bsmin
               ELSE
                 ext_data(jg)%atm%eai_t(jc,jb,1)    = MERGE(c_soil_urb,c_soil,lhave_urban)
                 ext_data(jg)%atm%r_bsmin(jc,jb)    = 50._wp ! previously hard-coded in TERRA
               ENDIF

               IF (lterra_urb .AND. ((itype_eisa == 2) .OR. (itype_eisa == 3))) THEN
                 ext_data(jg)%atm%eai_t(jc,jb,1)    = ext_data(jg)%atm%eai_t(jc,jb,1)                                        &
                                                    * (1.0_wp - ext_data(jg)%atm%urb_isa_t(jc,jb,1))
               END IF

               ! skin conductivity
               ext_data(jg)%atm%skinc_t(jc,jb,1)    = ext_data(jg)%atm%skinc_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))

               ! minimum stomatal resistance
               ext_data(jg)%atm%rsmin2d_t(jc,jb,1)  = ext_data(jg)%atm%rsmin(jc,jb)

               ! soil type
               ext_data(jg)%atm%soiltyp_t(jc,jb,1)  = ext_data(jg)%atm%soiltyp(jc,jb)


               ! Workaround for GLC2000 hole below 60 deg S
               ! (only necesary for old extpar files generated prior to 2014-01-31)
               IF (is_frglac_in(jg)) THEN
                 IF (tile_frac(ext_data(jg)%atm%lc_class_t(jc,jb,1))<=0._wp) &
                   ext_data(jg)%atm%lc_class_t(jc,jb,1) = ext_data(jg)%atm%i_lc_snow_ice
               ENDIF

               ! static index list and corresponding counter
               ext_data(jg)%atm%idx_lst_lp_t(i_count,jb,1)  = jc
               ext_data(jg)%atm%lp_count_t(jb,1)            = i_count

               ! initialize dynamic index list (in case of lsnowtile=true) with the same values
               ext_data(jg)%atm%idx_lst_t(i_count,jb,1) = jc
               ext_data(jg)%atm%gp_count_t(jb,1)        = i_count

               ! initialize snowtile flag with 1 if the tile is eligible for separate treatment of
               ! a snow-covered and a snow-free part, otherwise with -1
               IF (ext_data(jg)%atm%snowtile_lcc(ext_data(jg)%atm%lc_class_t(jc,jb,1))) THEN
                 ext_data(jg)%atm%snowtile_flag_t(jc,jb,1)  = 1
               ELSE
                 ext_data(jg)%atm%snowtile_flag_t(jc,jb,1)  = -1
               ENDIF

             ELSE
               ext_data(jg)%atm%lc_frac_t(jc,jb,:)  = 0._wp ! to be really safe

               DO i_lu = 1, ntiles_lnd
                 lu_subs = MAXLOC(tile_frac,1,tile_mask)
                 ! Note that we have to take into account points with fr_land > frlnd_thrhld but
                 ! maximum_tile_frac <frlndtile_thrhld (for tile=1).
                 ! This e.g. may be the case at non-dominant land points with very inhomogeneous land class coverage.
                 ! That's why checking for (tile_frac(lu_subs) >= frlndtile_thrhld) in the next If-statement is not enough.
                 ! We accept class fractions for tile=1 even if tile_frac << frlndtile_thrhld.
                 !
                 ! The additional check tile_frac(lu_subs) >= 1.e-03_wp is only added for backward compatibility and is
                 ! required by all extpar-files generated prior to 2014-01-30. In these files it is not assured that
                 ! SUM(tile_frac(:))=1. I.e. glacier below 60degS are missing, such that SUM(tile_frac(:))=0 in these cases.
                 !
                 IF ( (i_lu==1 .AND. tile_frac(lu_subs) >= 1.e-03_wp) .OR. (tile_frac(lu_subs) >= frlndtile_thrhld) ) THEN
                   it_count(i_lu)    = it_count(i_lu) + 1
                   tile_mask(lu_subs)= .FALSE.

                   ! static index list and corresponding counter
                   ext_data(jg)%atm%idx_lst_lp_t(it_count(i_lu),jb,i_lu) = jc
                   ext_data(jg)%atm%lp_count_t(jb,i_lu)                  = it_count(i_lu)

                   ! initialize dynamic index list (in case of lsnowtile=true) with the same values
                   ext_data(jg)%atm%idx_lst_t(it_count(i_lu),jb,i_lu) = jc
                   ext_data(jg)%atm%gp_count_t(jb,i_lu)               = it_count(i_lu)

                   ! initialize snowtile flag with 1 if the tile is eligible for separate treatment of
                   ! a snow-covered and a snow-free part, otherwise with -1
                   IF (ext_data(jg)%atm%snowtile_lcc(lu_subs)) THEN
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)  = 1
                   ELSE
                     ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)  = -1
                   ENDIF

                   ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = tile_frac(lu_subs)
                   ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = lu_subs
                   IF (lu_subs == ext_data(jg)%atm%i_lc_urban) lhave_urban = .true.
                 ELSE
                   EXIT ! no more land cover classes exceeding the threshold
                 ENDIF

               END DO

               ! fix for non-dominant land points
               !!! only active for 'old' extpar datasets (20131009 and earlier) !!!
               IF (ext_data(jg)%atm%fr_land(jc,jb) < 0.5_wp) THEN
                 IF (ext_data(jg)%atm%soiltyp(jc,jb) == 9) THEN  ! sea water
                   ! reset soil type to sandy loam ...
                   ext_data(jg)%atm%soiltyp(jc,jb) = 4
                 ENDIF
                 IF (ptr_ndviratio(jc,jb) <= 0.0_wp) THEN  ! here: 0=extpar_missval
                   ! ... and reset ndviratio
                   ptr_ndviratio(jc,jb) = 0.5_wp
                 ENDIF
                 IF (ext_data(jg)%atm%ndvi_max(jc,jb) <= 0.0_wp ) THEN  ! here: 0=extpar_missval
                   ! ... and reset ndvi_max to meaningful value (needed for plant cover)
                   ext_data(jg)%atm%ndvi_max(jc,jb) = 0.8_wp
                 ENDIF
               ENDIF

               sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))


               DO i_lu = 1, ntiles_lnd

                 !  Workaround for GLC2000 hole below 60 deg S
                 ! (only necesary for old extpar files generated prior to 2014-01-31)
                 IF (is_frglac_in(jg)) THEN
                   IF ( sum_frac < 1.e-10_wp) THEN
                     IF (i_lu == 1) THEN
                       it_count(i_lu)    = it_count(i_lu) + 1
                       ! static index list and corresponding counter
                       ext_data(jg)%atm%idx_lst_lp_t(it_count(i_lu),jb,i_lu) = jc
                       ext_data(jg)%atm%lp_count_t(jb,i_lu)               = it_count(i_lu)

                       ! initialize dynamic index list (in case of lsnowtile=true) with the same values
                       ext_data(jg)%atm%idx_lst_t(it_count(i_lu),jb,i_lu) = jc
                       ext_data(jg)%atm%gp_count_t(jb,i_lu)               = it_count(i_lu)

                       ! the snowtile flag is initialized with 1 here because the snow/ice class is
                       ! supposed to be located on the snowtiles (if activated)
                       ext_data(jg)%atm%snowtile_flag_t(jc,jb,i_lu)         = 1

                       ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = ext_data(jg)%atm%i_lc_snow_ice
                       ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = ext_data(jg)%atm%fr_land(jc,jb)
                     ELSE
                       ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = -1
                       ext_data(jg)%atm%lc_frac_t(jc,jb,i_lu)  = 0._wp
                     ENDIF
                   END IF  ! sum_frac < 1.e-10_wp
                 ENDIF  ! is_frglac_in(jg)

                 ! consistency corrections for glaciered points
                 !
                 ! a) plausibility check for glacier points based on T2M climatology (if available):
                 !    if the warmest month exceeds 10 deg C, then it is unlikely for glaciers to exist
                 !    This correction requires a monthly T2M climatology, which is available only 
                 !    if itype_vegetation_cycle > 1
                 IF (itype_vegetation_cycle > 1) THEN
                   IF (ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) == ext_data(jg)%atm%i_lc_snow_ice) THEN
                     ! Calculate height-corrected annual maximum of T2M climatology, 
                     ! including contribution from SSO standard deviation. 
                     ! This is used below to reset misclassified glacier points 
                     ! (e.g. salt lakes) to bare soil
                     !
                     t2mclim_hc = MAXVAL(ext_data(jg)%atm_td%t2m_m(jc,jb,:)) + dtdz_clim *             &
                       ( ext_data(jg)%atm%topography_c(jc,jb) + 1.5_wp*ext_data(jg)%atm%sso_stdh(jc,jb) - &
                         ext_data(jg)%atm%topo_t2mclim(jc,jb) )

                     IF (t2mclim_hc > (tmelt + 10._wp)) THEN
                       ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) = ext_data(jg)%atm%i_lc_bare_soil
                       ext_data(jg)%atm%fr_glac(jc,jb)     = 0._wp
                       icount_falseglac(jb) = icount_falseglac(jb) + 1
                     ENDIF
                   ENDIF
                 ENDIF

                 !
                 lu_subs = ext_data(jg)%atm%lc_class_t(jc,jb,i_lu)
                 IF (lu_subs < 0) CYCLE
                 !
                 ! Urban Canopy Parameters (UCPs)
                 !
                 ! impervious surface area fraction of the urban canopy
                 ext_data(jg)%atm%urb_isa_t(jc,jb,i_lu)       = ext_data(jg)%atm%fr_paved_lcc(lu_subs)

                 IF (lterra_urb) THEN
                   ! impervious surface area (ISA) fraction
                   ext_data(jg)%atm%fr_paved_t(jc,jb,i_lu)    = ext_data(jg)%atm%fr_paved_lcc(lu_subs)

                   ! building area fraction with respect to urban tile
                   ext_data(jg)%atm%urb_fr_bld_t(jc,jb,i_lu)  = 0.667_wp

                   ! street canyon H/W ratio
                   ext_data(jg)%atm%urb_h2w_t(jc,jb,i_lu)     = 1.5_wp

                   ! surface area index of the urban canopy
                   ext_data(jg)%atm%urb_ai_t(jc,jb,i_lu)      = (1.0_wp + 2.0_wp * ext_data(jg)%atm%urb_h2w_t(jc,jb,i_lu)) &
                     &                                        * (1.0_wp - ext_data(jg)%atm%urb_fr_bld_t(jc,jb,i_lu))       &
                     &                                        + ext_data(jg)%atm%urb_fr_bld_t(jc,jb,i_lu)

                   ! albedo reduction factor for the urban canopy
                   ! Reduce the effective albedo according to the building density,
                   ! the reduction factor is based on Monte-Carlo simulations.
                   ext_data(jg)%atm%urb_alb_red_t(jc,jb,i_lu) = EXP(-0.6_wp * ext_data(jg)%atm%urb_h2w_t(jc,jb,i_lu))      &
                     &                                        * (1.0_wp - ext_data(jg)%atm%urb_fr_bld_t(jc,jb,i_lu))       &
                     &                                        + ext_data(jg)%atm%urb_fr_bld_t(jc,jb,i_lu)

                   ! building height
                   ext_data(jg)%atm%urb_h_bld_t(jc,jb,i_lu)   = 15._wp

                   ! thermal albedo of urban material
                   ext_data(jg)%atm%urb_alb_th_t(jc,jb,i_lu)  = 0.14_wp

                   ! solar albedo of urban material, times albedo reduction factor for the urban canopy
                   ext_data(jg)%atm%urb_alb_so_t(jc,jb,i_lu)  = 0.101_wp * ext_data(jg)%atm%urb_alb_red_t(jc,jb,i_lu)

                   ! volumetric heat capacity of urban material
                   ext_data(jg)%atm%urb_hcap_t(jc,jb,i_lu)    = 1250000._wp

                   ! thermal conductivity of urban material
                   ext_data(jg)%atm%urb_hcon_t(jc,jb,i_lu)    = 0.767_wp

                   ! anthropogenic heat flux
                   ext_data(jg)%atm%ahf_t(jc,jb,i_lu)         = ext_data(jg)%atm%ahf_lcc(lu_subs)
                 ENDIF
                 !
                 ! Vegetation Parameters etc.
                 !
                 ! root depth
                 ext_data(jg)%atm%rootdp_t (jc,jb,i_lu)  = ext_data(jg)%atm%rootdmax_lcc(lu_subs)

                 ! plant cover
                 ext_data(jg)%atm%plcov_t  (jc,jb,i_lu)  = ptr_ndviratio(jc,jb)                                            &
                   & * MIN(ext_data(jg)%atm%ndvi_max(jc,jb),ext_data(jg)%atm%plcovmax_lcc(lu_subs))

                 ! transpiration area index
                 ext_data(jg)%atm%tai_t    (jc,jb,i_lu)  = ext_data(jg)%atm%plcov_t(jc,jb,i_lu)                            &
                   & * ext_data(jg)%atm%laimax_lcc(lu_subs)

                 ! surface area index
                 IF (lterra_urb) THEN
                   ext_data(jg)%atm%sai_t  (jc,jb,i_lu)  = c_lnd * (1.0_wp - ext_data(jg)%atm%urb_isa_t(jc,jb,i_lu))       &
                                          + ext_data(jg)%atm%urb_ai_t(jc,jb,i_lu) * ext_data(jg)%atm%urb_isa_t(jc,jb,i_lu) &
                                          + ext_data(jg)%atm%tai_t(jc,jb,i_lu)
                 ELSE
                   ext_data(jg)%atm%sai_t  (jc,jb,i_lu)  = c_lnd + ext_data(jg)%atm%tai_t(jc,jb,i_lu)
                 END IF

                 ! evaporative soil area index
                 IF (icpl_da_sfcevap >= 4 .OR. itype_evsl == 5) THEN
                   ext_data(jg)%atm%eai_t (jc,jb,i_lu)   = MERGE(0.75_wp,2.0_wp,lu_subs == ext_data(jg)%atm%i_lc_urban)
                   ext_data(jg)%atm%r_bsmin(jc,jb)       = cr_bsmin
                   ! on non-urban tiles, the eai is reduced only if no urban tile is present on the grid point
                   IF (.NOT. lhave_urban) ext_data(jg)%atm%eai_t(jc,jb,i_lu) =                                             &
                                          2.0_wp - 1.25_wp*tile_frac(ext_data(jg)%atm%i_lc_urban)
                 ELSE
                   ext_data(jg)%atm%eai_t (jc,jb,i_lu)   = MERGE(c_soil_urb,c_soil,lu_subs == ext_data(jg)%atm%i_lc_urban)
                   ext_data(jg)%atm%r_bsmin(jc,jb)       = 50._wp ! previously hard-coded in TERRA
                 ENDIF

                 IF (lterra_urb .AND. ((itype_eisa == 2) .OR. (itype_eisa == 3))) THEN
                   ext_data(jg)%atm%eai_t (jc,jb,i_lu)   = ext_data(jg)%atm%eai_t(jc,jb,i_lu)                              &
                                                         * (1.0_wp - ext_data(jg)%atm%urb_isa_t(jc,jb,i_lu))
                 END IF

                 ! skin conductivity
                 lat = p_patch(jg)%cells%center(jc,jb)%lat*rad2deg
                 IF (itype_lndtbl == 4 .AND. lat > -10._wp .AND. lat < 42.5_wp) THEN
                   ext_data(jg)%atm%skinc_t(jc,jb,i_lu)  = MIN(200._wp,ext_data(jg)%atm%skinc_lcc(lu_subs)*                &
                                                          (1._wp+MIN(1._wp,0.4_wp*(42.5_wp-lat),0.4_wp*(lat+10._wp))) )
                 ELSE
                   ext_data(jg)%atm%skinc_t(jc,jb,i_lu)  = ext_data(jg)%atm%skinc_lcc(lu_subs)
                 ENDIF

                 ! minimum stomatal resistance
                 ext_data(jg)%atm%rsmin2d_t(jc,jb,i_lu)  = ext_data(jg)%atm%stomresmin_lcc(lu_subs)

                 ! soil type
                 ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu)  = ext_data(jg)%atm%soiltyp(jc,jb)


                 ! consistency corrections for glaciered points (continued)
                 !
                 ! b) set soiltype to ice if landuse = ice (already done in extpar for dominant glacier points)
                 IF (ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) == ext_data(jg)%atm%i_lc_snow_ice) &
                   & ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) = 1
                 !
                 ! c) set soiltype to rock or sandy loam if landuse /= ice and soiltype = ice
                 IF (ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) == 1 .AND. &
                   & ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) /= ext_data(jg)%atm%i_lc_snow_ice) THEN
                   IF (ext_data(jg)%atm%lc_class_t(jc,jb,i_lu) == ext_data(jg)%atm%i_lc_bare_soil) THEN
                     ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) = 2 ! assume rock in case of bare soil
                   ELSE
                     ext_data(jg)%atm%soiltyp_t(jc,jb,i_lu) = 4 ! otherwise assume sandy loam
                   ENDIF
                 ENDIF

               END DO
             END IF ! ntiles
           ELSE  ! fr_land(jc,jb)<= frlnd_thrhld
             ! measures for land-specific fields that are also defined on non-dominating land points:
             !
             ! for_d, for_e: only accessed via land point index list
             !               -> non-dominant land points do not matter when running without tiles
             ! rootdp, rsmin, lai_mx, plcov_mx, ndvi_max : only accessed via land point index list
             ! ndvi_mrat -> ndviratio : only accessed via land point index list
             ! soiltyp :
             !
             ! glacier fraction
             ext_data(jg)%atm%fr_glac(jc,jb) = 0._wp  ! for frlnd_thrhld=0.5 (i.e. without tiles) this is
                                                      ! identical to what has previously been done within
                                                      ! EXTPAR crosschecks.
           ENDIF



           !
           ! searching for lake-points
           !
           IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN
             i_count_flk=i_count_flk+1
             ext_data(jg)%atm%list_lake%idx(i_count_flk,jb) = jc  ! write index of lake-points
             ext_data(jg)%atm%list_lake%ncount(jb) = i_count_flk

             ! set land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,isub_lake) = ext_data(jg)%atm%i_lc_water
             ! set also area fractions
             ext_data(jg)%atm%lc_frac_t(jc,jb,isub_lake)  = ext_data(jg)%atm%fr_lake(jc,jb)

             ! set surface area index (needed by turbtran)
             ext_data(jg)%atm%sai_t    (jc,jb,isub_lake)  = c_sea
           ENDIF

           !
           ! searching for sea points
           !
           IF (1._wp-ext_data(jg)%atm%fr_land(jc,jb)-ext_data(jg)%atm%fr_lake(jc,jb) &
             &   >= frsea_thrhld) THEN
             i_count_sea=i_count_sea + 1

             ! Ensure that sea and lake tiles do not coexist on any grid point (this is already done
             ! in extpar but might not be fulfilled in data sets generated with other software)
             IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN
               CALL finish('', "Lake and sea tiles must not coexist on any grid point")
             ENDIF

             ext_data(jg)%atm%list_sea%idx(i_count_sea,jb) = jc  ! write index of sea-points
             ext_data(jg)%atm%list_sea%ncount(jb) = i_count_sea
             ! set land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,isub_water) = ext_data(jg)%atm%i_lc_water
             ! set also area fractions
             ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water)  = 1._wp                        &
               &         -ext_data(jg)%atm%fr_land(jc,jb) - ext_data(jg)%atm%fr_lake(jc,jb)
             ! fix potential truncation errors
             IF (ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water) < 1.e-10_wp) &
               ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water) = 0._wp

             ! set surface area index (needed by turbtran)
             ext_data(jg)%atm%sai_t    (jc,jb,isub_water)  = c_sea

             ! set land-cover class for seaice tile
             ! sea-ice and sea have the same land cover class. This is consistent with the
             ! applied GRIB2 tile template, where sea-ice and sea are treated as two
             ! attributes of the same tile. Per definition, different attributes of the
             ! same tile have the same land-cover class
             ext_data(jg)%atm%lc_class_t(jc,jb,isub_seaice) = ext_data(jg)%atm%lc_class_t(jc,jb,isub_water)
           ENDIF

           !
           ! index list for seaice points is generated in mo_nwp_sfc_utils/init_sea_lists
           !
           ! note that in principle, sai_t for seaice should be different from sai_t for
           ! open water points. However, for the time being, sai_t=c_sea is also used
           ! for seaice points.

         END DO ! jc



         ! Inflate dominating land-tile fractions to fr_land or fr_land + fr_lake, depending
         ! on whether a lake tile is present (fr_lake >= frlake_thrhld), or not
         ! (fr_lake < frlake_thrhld).
         IF (ntiles_lnd > 1) THEN
           ! Inflate fractions for land points

           DO ic = 1, ext_data(jg)%atm%list_land%ncount(jb)

             jc = ext_data(jg)%atm%list_land%idx(ic,jb)

             ! sum up fractions of dominating land tiles
             sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))

             IF (ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld) THEN
               ! cell with lake point
               ! inflate dominating land fractions to fr_land
               scalfac = ext_data(jg)%atm%fr_land(jc,jb)/sum_frac
             ELSE
               ! cell without lake point
               ! inflate dominating land fractions to (fr_land + fr_lake)
               scalfac = (ext_data(jg)%atm%fr_land(jc,jb) + ext_data(jg)%atm%fr_lake(jc,jb))/sum_frac
             ENDIF

             ! inflate land fractions
             DO jt = 1, ntiles_total
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt) = ext_data(jg)%atm%lc_frac_t(jc,jb,jt) * scalfac
             ENDDO
           ENDDO  ! ic


           ! Inflate fractions for
           ! - sea-water only points
           ! - lake only points.
           ! As a side effect, fractions for land-only points (with 0<fr_sea<frsea_thrhld)
           ! are also corrected.
           ! Note that, for simplicity, we loop over all points. At mixed land/water points this
           ! should do no harm, since these have already been inflated in the loop above.
           DO jc = i_startidx, i_endidx
             sum_frac = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd)) + &
                        SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water:isub_lake))

             DO jt = 1, ntiles_total + MIN(2,ntiles_water)
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt) = ext_data(jg)%atm%lc_frac_t(jc,jb,jt) / sum_frac
             ENDDO
           ENDDO  ! jc

           ! Ensure consistency between fr_land and the adjusted sum of the tile fractions
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%fr_land(jc,jb)     = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))
           ENDDO  ! jc

         ELSE ! overwrite fractional settings over water points if tile approach is turned off
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%lc_frac_t(jc,jb,1) = 1._wp
           ENDDO
         ENDIF


         ! Compute inverse of fr_land based on land tile fractions.
         ! Required for proper aggregation of land-only variables
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%inv_frland_from_tiles(jc,jb) = 0._wp
           zfr_land = SUM(ext_data(jg)%atm%lc_frac_t(jc,jb,1:ntiles_lnd))

           IF (zfr_land > 0._wp) THEN
             ext_data(jg)%atm%inv_frland_from_tiles(jc,jb) = 1._wp/zfr_land
           ENDIF

         ENDDO  ! jc


         IF (lsnowtile) THEN ! copy external data fields to snow tile grid points
           DO jt = ntiles_lnd+1, ntiles_total

             jt_in = jt - ntiles_lnd
             ext_data(jg)%atm%lp_count_t(jb,jt)     = ext_data(jg)%atm%lp_count_t(jb,jt_in)
             ext_data(jg)%atm%idx_lst_lp_t(:,jb,jt) = ext_data(jg)%atm%idx_lst_lp_t(:,jb,jt_in)
             !
             ! the following two fields are reset in init_snowtile_lists, but presetting them here
             ! avoids complications in initicon
             ext_data(jg)%atm%gp_count_t(jb,jt)     = ext_data(jg)%atm%gp_count_t(jb,jt_in)
             ext_data(jg)%atm%idx_lst_t(:,jb,jt)    = ext_data(jg)%atm%idx_lst_t(:,jb,jt_in)

!CDIR NODEP
             DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
               jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
!
               ext_data(jg)%atm%urb_isa_t(jc,jb,jt)     = ext_data(jg)%atm%urb_isa_t(jc,jb,jt_in)
               IF (lterra_urb) THEN
                 ext_data(jg)%atm%fr_paved_t(jc,jb,jt)    = ext_data(jg)%atm%fr_paved_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_ai_t(jc,jb,jt)      = ext_data(jg)%atm%urb_ai_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_alb_red_t(jc,jb,jt) = ext_data(jg)%atm%urb_alb_red_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_fr_bld_t(jc,jb,jt)  = ext_data(jg)%atm%urb_fr_bld_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_h2w_t(jc,jb,jt)     = ext_data(jg)%atm%urb_h2w_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_h_bld_t(jc,jb,jt)   = ext_data(jg)%atm%urb_h_bld_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_alb_th_t(jc,jb,jt)  = ext_data(jg)%atm%urb_alb_th_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_alb_so_t(jc,jb,jt)  = ext_data(jg)%atm%urb_alb_so_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_hcap_t(jc,jb,jt)    = ext_data(jg)%atm%urb_hcap_t(jc,jb,jt_in)
                 ext_data(jg)%atm%urb_hcon_t(jc,jb,jt)    = ext_data(jg)%atm%urb_hcon_t(jc,jb,jt_in)
                 ext_data(jg)%atm%ahf_t(jc,jb,jt)         = ext_data(jg)%atm%ahf_t(jc,jb,jt_in)
               ENDIF
!
               ext_data(jg)%atm%rootdp_t(jc,jb,jt)      = ext_data(jg)%atm%rootdp_t(jc,jb,jt_in)
               ext_data(jg)%atm%plcov_t(jc,jb,jt)       = ext_data(jg)%atm%plcov_t(jc,jb,jt_in)
               ext_data(jg)%atm%tai_t(jc,jb,jt)         = ext_data(jg)%atm%tai_t(jc,jb,jt_in)
               ext_data(jg)%atm%sai_t(jc,jb,jt)         = ext_data(jg)%atm%sai_t(jc,jb,jt_in)
               ext_data(jg)%atm%eai_t(jc,jb,jt)         = ext_data(jg)%atm%eai_t(jc,jb,jt_in)
               ext_data(jg)%atm%skinc_t(jc,jb,jt)       = ext_data(jg)%atm%skinc_t(jc,jb,jt_in)
               ext_data(jg)%atm%rsmin2d_t(jc,jb,jt)     = ext_data(jg)%atm%rsmin2d_t(jc,jb,jt_in)
               ext_data(jg)%atm%soiltyp_t(jc,jb,jt)     = ext_data(jg)%atm%soiltyp_t(jc,jb,jt_in)
               ext_data(jg)%atm%lc_class_t(jc,jb,jt)    = ext_data(jg)%atm%lc_class_t(jc,jb,jt_in)
               ext_data(jg)%atm%lc_frac_t(jc,jb,jt)     = ext_data(jg)%atm%lc_frac_t(jc,jb,jt_in)
             ENDDO

           ENDDO
         ENDIF


         ! Initialize frac_t with lc_frac_t on all static grid points
         ! Recall: frac_t differs from lc_frac_t in the presence of time-dependent sub-lists
         !         (snow tiles or sea ice)
         ! In this case, frac_t refers to the time-dependent sub-tiles.
         ! ** Aggregation operations always must use frac_t **
         DO jt = 1, ntiles_lnd
           DO jc = i_startidx, i_endidx
             ext_data(jg)%atm%frac_t(jc,jb,jt)  = ext_data(jg)%atm%lc_frac_t(jc,jb,jt)
           ENDDO
         ENDDO
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%frac_t(jc,jb,isub_water) = ext_data(jg)%atm%lc_frac_t(jc,jb,isub_water)
         ENDDO
         DO jc = i_startidx, i_endidx
           ext_data(jg)%atm%frac_t(jc,jb,isub_lake)  = ext_data(jg)%atm%lc_frac_t(jc,jb,isub_lake)
           !
           ! If tiles are active, ensure consistency between fr_lake and the rescaled tile fraction
           IF (ntiles_lnd > 1) ext_data(jg)%atm%fr_lake(jc,jb) = ext_data(jg)%atm%frac_t(jc,jb,isub_lake)
         ENDDO
         ! frac_t(jc,jb,isub_seaice) is set in init_sea_lists

       END DO !jb
#ifndef __SX__
!$OMP END DO

!$OMP SINGLE
#endif
       ! Some useful diagnostics
       npoints      = ext_data(jg)%atm%list_land %get_sum_global(i_startblk,i_endblk)
       npoints_sea  = ext_data(jg)%atm%list_sea  %get_sum_global(i_startblk,i_endblk)
       npoints_lake = ext_data(jg)%atm%list_lake %get_sum_global(i_startblk,i_endblk)
       !
       WRITE(message_text,'(a,i3,a,i10)') 'Number of land points in domain',jg,':', npoints
       CALL message('', TRIM(message_text))
       WRITE(message_text,'(a,i3,a,i10)') 'Number of sea points in domain',jg,':', npoints_sea
       CALL message('', TRIM(message_text))
       WRITE(message_text,'(a,i3,a,i10)') 'Number of lake points in domain',jg,':', npoints_lake
       CALL message('', TRIM(message_text))
       !
       !
       DO i_lu = 1, ntiles_lnd
         npoints = SUM(ext_data(jg)%atm%gp_count_t(i_startblk:i_endblk,i_lu))
         npoints = global_sum_array(npoints)
         WRITE(message_text,'(a,i2,a,i10)') 'Number of points in tile',i_lu,':',npoints
         CALL message('', TRIM(message_text))
       ENDDO

       npoints = SUM(icount_falseglac(i_startblk:i_endblk))
       npoints = global_sum_array(npoints)
       WRITE(message_text,'(a,i3,a,i10)') 'Number of corrected false glacier points in domain',jg,':', npoints
       CALL message('', TRIM(message_text))
!$OMP END SINGLE NOWAIT


       !
       ! For consistency: remove depth_lk information, where fr_lake < frlake_thrhld.
       ! Boundary interpolation zone of nested domains is explicitly included.
       !
       ! In case of ntiles_lnd > 1, fr_lake ranges from 0<=fr_lake<=1 at nest 
       ! boundaries, whereas at prognostic points fr_lake ranges from 
       ! frlake_thrhld<=fr_lake<=1. I.e. fr_lake consistency adjustment 
       ! was not performed for nest boundary.
       rl_start = 1
       rl_end   = min_rlcell_int

       i_startblk = p_patch(jg)%cells%start_blk(rl_start,1)
       i_endblk   = p_patch(jg)%cells%end_blk(rl_end,i_nchdom)

!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
       DO jb=i_startblk, i_endblk

         CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
            & i_startidx, i_endidx, rl_start, rl_end)

         DO jc = i_startidx, i_endidx
           !
           ! For consistency: remove depth_lk information, where fr_lake=0
           ext_data(jg)%atm%depth_lk(jc,jb) = MERGE(ext_data(jg)%atm%depth_lk(jc,jb), &
             &                                      -1._wp,                           &
             &                                      ext_data(jg)%atm%fr_lake(jc,jb) >= frlake_thrhld)
         ENDDO

       ENDDO  ! jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      DEALLOCATE(icount_falseglac)

      !$ACC UPDATE DEVICE(ext_data(jg)%atm%list_sea%idx, ext_data(jg)%atm%list_sea%ncount) ASYNC(1)
    END DO  !jg

  END SUBROUTINE init_index_lists



  !-------------------------------------------------------------------------
  !! Diagnose aggregated external fields
  !!
  !! Aggregated external fields are diagnosed based on tile based external
  !! fields. In addition, fr_land, fr_lake and depth_lk are re-diagnosed,
  !! in order to be consistent with tile-information. Note that the latter 
  !! re-diagnosis has been moved to init_index_lists in order not to 
  !! compromise restart reproducibility.
  !!
  SUBROUTINE diagnose_ext_aggr (p_patch, ext_data)

    TYPE(t_patch)        , INTENT(IN)    :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data

    INTEGER  :: jb,jt,ic,jc
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk    !> blocks
    INTEGER  :: i_startidx, i_endidx
    INTEGER  :: i_count
    REAL(wp) :: area_frac

!!$    CHARACTER(len=max_char_length), PARAMETER :: &
!!$      routine = modname//':diagnose_ext_aggr'

    !-------------------------------------------------------------------------


    ! exclude the boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! Fill nest boundary points of sai with c_sea because the initial call of turbtran
    ! may produce invalid operations otherwise
    ext_data%atm%sai(:,1:i_startblk) = c_sea

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,ic,i_startidx,i_endidx,i_count,jc,area_frac)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

!
      ext_data%atm%urb_isa    (i_startidx:i_endidx,jb) = 0._wp
      IF (lterra_urb) THEN
        ext_data%atm%fr_paved   (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_ai     (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_alb_red(i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_fr_bld (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_h2w    (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_h_bld  (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_alb_th (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_alb_so (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_hcap   (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%urb_hcon   (i_startidx:i_endidx,jb) = 0._wp
        ext_data%atm%ahf        (i_startidx:i_endidx,jb) = 0._wp
      ENDIF
!
      ext_data%atm%plcov      (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%rootdp     (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%lai        (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%skinc      (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%rsmin      (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%tai        (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%eai        (i_startidx:i_endidx,jb) = 0._wp
      ext_data%atm%sai        (i_startidx:i_endidx,jb) = 0._wp


      DO jt = 1, ntiles_total
        i_count = ext_data%atm%gp_count_t(jb,jt)
        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_t(ic,jb,jt)

          ! note that frac_t must be re-scaled such that sum(frac_t(1:ntiles_lnd)) = 1
          ! therefore we multiply by inv_frland_from_tiles
          area_frac = ext_data%atm%frac_t(jc,jb,jt)           &
            &       * ext_data%atm%inv_frland_from_tiles(jc,jb)
!
          ! impervious surface area fraction of the urban canopy (aggregated)
          ext_data%atm%urb_isa(jc,jb) = ext_data%atm%urb_isa(jc,jb)           &
              &              + ext_data%atm%urb_isa_t(jc,jb,jt) * area_frac

          IF (lterra_urb) THEN
            ! impervious surface area (ISA) fraction (aggregated)
            ext_data%atm%fr_paved(jc,jb) = ext_data%atm%fr_paved(jc,jb)       &
              &              + ext_data%atm%fr_paved_t(jc,jb,jt) * area_frac

            ! surface area index of the urban canopy (aggregated)
            ext_data%atm%urb_ai(jc,jb) = ext_data%atm%urb_ai(jc,jb)           &
              &              + ext_data%atm%urb_ai_t(jc,jb,jt) * area_frac

            ! albedo reduction factor for the urban canopy (aggregated)
            ext_data%atm%urb_alb_red(jc,jb) = ext_data%atm%urb_alb_red(jc,jb) &
              &              + ext_data%atm%urb_alb_red_t(jc,jb,jt) * area_frac

            ! building area fraction with respect to urban tile (aggregated)
            ext_data%atm%urb_fr_bld(jc,jb) = ext_data%atm%urb_fr_bld(jc,jb)   &
              &              + ext_data%atm%urb_fr_bld_t(jc,jb,jt) * area_frac

            ! street canyon H/W ratio (aggregated)
            ext_data%atm%urb_h2w(jc,jb) = ext_data%atm%urb_h2w(jc,jb)         &
              &              + ext_data%atm%urb_h2w_t(jc,jb,jt) * area_frac

            ! building height (aggregated)
            ext_data%atm%urb_h_bld(jc,jb) = ext_data%atm%urb_h_bld(jc,jb)     &
              &              + ext_data%atm%urb_h_bld_t(jc,jb,jt) * area_frac

            ! thermal albedo of urban material (aggregated)
            ext_data%atm%urb_alb_th(jc,jb) = ext_data%atm%urb_alb_th(jc,jb)   &
              &              + ext_data%atm%urb_alb_th_t(jc,jb,jt) * area_frac

            ! solar albedo of urban material (aggregated)
            ext_data%atm%urb_alb_so(jc,jb) = ext_data%atm%urb_alb_so(jc,jb)   &
              &              + ext_data%atm%urb_alb_so_t(jc,jb,jt) * area_frac

            ! volumetric heat capacity of urban material (aggregated)
            ext_data%atm%urb_hcap(jc,jb) = ext_data%atm%urb_hcap(jc,jb)       &
              &              + ext_data%atm%urb_hcap_t(jc,jb,jt) * area_frac

            ! thermal conductivity of urban material (aggregated)
            ext_data%atm%urb_hcon(jc,jb) = ext_data%atm%urb_hcon(jc,jb)       &
              &              + ext_data%atm%urb_hcon_t(jc,jb,jt) * area_frac

            ! anthropogenic heat flux (aggregated)
            ext_data%atm%ahf(jc,jb) = ext_data%atm%ahf(jc,jb)                 &
              &              + ext_data%atm%ahf_t(jc,jb,jt) * area_frac
          ENDIF
!
          ! plant cover (aggregated)
          ext_data%atm%plcov(jc,jb) = ext_data%atm%plcov(jc,jb)       &
            &              + ext_data%atm%plcov_t(jc,jb,jt) * area_frac

          ! root depth (aggregated)
          ext_data%atm%rootdp(jc,jb) = ext_data%atm%rootdp(jc,jb)     &
            &              + ext_data%atm%rootdp_t(jc,jb,jt) * area_frac

          ! surface area index (aggregated)
          ext_data%atm%lai(jc,jb) = ext_data%atm%lai(jc,jb)           &
            &              + ( ext_data%atm%tai_t(jc,jb,jt)           &
            &              /(ext_data%atm%plcov_t(jc,jb,jt)+dbl_eps) * area_frac )

          ! evaporative soil area index (aggregated)
          ext_data%atm%eai(jc,jb) = ext_data%atm%eai(jc,jb)           &
            &              +  ext_data%atm%eai_t(jc,jb,jt) * area_frac

          ! transpiration area index (aggregated)
          ext_data%atm%tai(jc,jb) = ext_data%atm%tai(jc,jb)           &
            &              +  ext_data%atm%tai_t(jc,jb,jt) * area_frac

          ! skin conductivity (aggregated)
          ext_data%atm%skinc(jc,jb) = ext_data%atm%skinc(jc,jb)       &
            &              + ext_data%atm%skinc_t(jc,jb,jt) * area_frac

          ! minimal stomata resistance (aggregated)
          ext_data%atm%rsmin(jc,jb) = ext_data%atm%rsmin(jc,jb)       &
            &              + ext_data%atm%rsmin2d_t(jc,jb,jt) * area_frac

        ENDDO  !ic
      ENDDO  !jt


      ! aggregate fields with water tiles
      DO jt = 1,ntiles_total + ntiles_water
        DO jc = i_startidx, i_endidx

          area_frac = ext_data%atm%frac_t(jc,jb,jt)

          ! surface area index (aggregated)
          ext_data%atm%sai(jc,jb) = ext_data%atm%sai(jc,jb)           &
            &             +  ext_data%atm%sai_t(jc,jb,jt) * area_frac
        ENDDO  ! jc
      ENDDO  !jt

    ENDDO  !jb
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE diagnose_ext_aggr



  !-------------------------------------------------------------------------
  !! Get interpolated field from monthly mean climatology
  !!
  !! Get interpolated field from monthly mean climatology. A linear interpolation
  !! in time between successive months is performed, assuming that the monthly field
  !! applies to the 15th of the month.
  !!
  SUBROUTINE interpol_monthly_mean(p_patch, mtime_date, monthly_means, out_field, out_diff)

    TYPE(datetime)   , INTENT(IN)  :: mtime_date
    TYPE(t_patch)    , INTENT(IN)  :: p_patch
    REAL(wp),          INTENT(IN)  :: monthly_means(:,:,:)  ! monthly mean climatology
    REAL(wp),          INTENT(OUT) :: out_field(:,:)        ! interpolated output field
    REAL(wp), OPTIONAL,INTENT(OUT) :: out_diff(:,:)         ! difference between adjacent monthly means

    INTEGER                             :: jc, jb               !< loop index
    INTEGER                             :: i_startblk, i_endblk
    INTEGER                             :: rl_start, rl_end
    INTEGER                             :: i_startidx, i_endidx
    INTEGER                             :: mo1, mo2             !< nearest months
    REAL(wp)                            :: zw1, zw2
    TYPE(datetime), POINTER             :: mtime_hour

    TYPE(t_time_interpolation_weights)  :: current_time_interpolation_weights
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dtime_string
    
    CHARACTER(len=*), PARAMETER :: routine = modname//': interpol_monthly_mean'

    !---------------------------------------------------------------
    ! Find the 2 nearest months mo1, mo2 and the weights zw1, zw2
    ! to the actual date and time

    mtime_hour => newDatetime(mtime_date)
    mtime_hour%time%minute = 0
    mtime_hour%time%second = 0
    mtime_hour%time%ms     = 0     
    current_time_interpolation_weights = calculate_time_interpolation_weights(mtime_hour)
    call deallocateDatetime(mtime_hour)    
    mo1 = current_time_interpolation_weights%month1
    mo2 = current_time_interpolation_weights%month2
    zw1 = current_time_interpolation_weights%weight1
    zw2 = current_time_interpolation_weights%weight2

    ! consistency check
    IF ((MIN(mo1,mo2) < 1) .OR. (MAX(mo1,mo2) > SIZE(monthly_means,3))) THEN
      CALL datetimeToString(mtime_date, dtime_string)
      CALL message('','Result of call to calculate_time_interpolation_weights:')
      WRITE (message_text,'(a,a)')      '   mtime_date = ', TRIM(dtime_string)
      CALL message('', message_text)
      WRITE (message_text,'(a,i2.2)')   '   mo1        = ', mo1
      CALL message('', message_text)
      WRITE (message_text,'(a,i2.2)')   '   mo2        = ', mo2
      CALL message('', message_text)
      WRITE (message_text,'(a,f25.15)') '   weight1    = ', zw1
      CALL message('', message_text)
      WRITE (message_text,'(a,f25.15)') '   weight2    = ', zw2
      CALL message('', message_text)      
      CALL finish(routine, "Error!")
    END IF

    ! include boundary interpolation zone and halo points 
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    ! Get interpolated field
    !
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx)
    DO jb=i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
         & i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx, i_endidx
        out_field(jc,jb) = zw1*monthly_means(jc,jb,mo1) &
          &              + zw2*monthly_means(jc,jb,mo2)
      ENDDO
      IF (PRESENT(out_diff)) THEN
        DO jc = i_startidx, i_endidx
          out_diff(jc,jb) = monthly_means(jc,jb,mo2) - monthly_means(jc,jb,mo1)
        ENDDO
      ENDIF
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE interpol_monthly_mean


  !-------------------------------------------------------------------------
  !! Improves specifiation of vegetation climatology based on monthly climatology
  !! of 2m-temperature. In particular, a distinction between deciduous and evergreen
  !! vegetation classes is made.
  !!
  SUBROUTINE vege_clim (p_patch, ext_data, nh_diag)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data
    TYPE(t_nh_diag),  INTENT(IN)         :: nh_diag

    INTEGER  :: jb,jt,ic,jc,i
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk,i_startidx, i_endidx
    INTEGER  :: i_count,ilu

    REAL(wp) :: t2mclim_hc(nproma),t_asyfac(nproma),tdiff_norm,wfac,dtdz_clim,trans_width,trh_bias, &
                skinc_fac,lat,scal,dtfac_skinc
    REAL(wp), DIMENSION(num_lcc) :: laimin,threshold_temp,temp_asymmetry,rd_fac

    INTEGER, PARAMETER :: nparam = 4  ! Number of parameters used in lookup table 

    REAL(wp), DIMENSION(num_lcc*nparam), TARGET :: vege_table ! < lookup table with control parameter specifications

    !-------------------------------------------------------------------------

    !                 lai_min   T_thresh   T_asy    rd_fac
    DATA vege_table /  0.5_wp,  283.0_wp,  2.0_wp,  5.0_wp,   & ! 1 - irrigated croplands
                   &   0.5_wp,  283.0_wp,  2.0_wp,  5.0_wp,   & ! 2 - rainfed croplands
                   &   0.7_wp,  282.0_wp,  1.0_wp,  4.0_wp,   & ! 3 - mosaic cropland (50-70%) - vegetation (20-50%)
                   &   0.7_wp,  281.0_wp,  1.0_wp,  3.0_wp,   & ! 4 - mosaic vegetation (50-70%) - cropland (20-50%)
                   &   5.0_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 5 - closed broadleaved evergreen forest
                   &   0.5_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 6 - closed broadleaved deciduous forest
                   &   0.5_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 7 - open broadleaved deciduous forest
                   &   3.0_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 8 - closed needleleaved evergreen forest
                   &   1.5_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 9 - open needleleaved evergreen or deciduous forest
                   &   1.5_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 10- mixed broadleaved and needleleaved forest
                   &   1.0_wp,  280.0_wp,  0.0_wp,  1.0_wp,   & ! 11- mosaic shrubland (50-70%) - grassland (20-50%)
                   &   1.0_wp,  280.0_wp,  0.0_wp,  1.0_wp,   & ! 12- mosaic grassland (50-70%) - shrubland (20-50%)
                   &   1.5_wp,  280.0_wp,  0.0_wp,  1.0_wp,   & ! 13- closed to open shrubland
                   &   1.0_wp,  279.0_wp,  0.0_wp,  1.0_wp,   & ! 14- closed to open herbaceous vegetation
                   &   0.3_wp,  278.0_wp,  0.0_wp,  1.0_wp,   & ! 15- sparse vegetation
                   &   5.0_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 16- closed to open forest regulary flooded
                   &   5.0_wp,  282.0_wp,  0.0_wp,  1.0_wp,   & ! 17- closed forest or shrubland permanently flooded
                   &   1.0_wp,  279.0_wp,  0.0_wp,  1.0_wp,   & ! 18- closed to open grassland regularly flooded
                   &   1.0_wp,  281.0_wp,  0.0_wp,  1.0_wp,   & ! 19- artificial surfaces
                   &   0.2_wp,  280.0_wp,  0.0_wp,  1.0_wp,   & ! 20- bare areas
                   &   0.0_wp,    0.0_wp,  0.0_wp,  1.0_wp,   & ! 21- water bodies
                   &   0.0_wp,    0.0_wp,  0.0_wp,  1.0_wp,   & ! 22- permanent snow and ice
                   &   0.0_wp,    0.0_wp,  0.0_wp,  1.0_wp    / ! 23- undefined

    !-------------------------------------------------------------------------

    CALL message('','Modify NDVI-based plant cover properties using T2M climatology')

    ! store table values according to landuse class
    ilu = 0
    DO i = 1, num_lcc*nparam, nparam
      ilu=ilu+1
      laimin(ilu)         = vege_table(i)
      threshold_temp(ilu) = vege_table(i+1)
      temp_asymmetry(ilu) = vege_table(i+2)
      rd_fac(ilu)         = vege_table(i+3)
    ENDDO

    ! climatological temperature gradient used for height correction of coarse input data
    dtdz_clim = -5.e-3_wp  ! -5 K/km

    ! transition width for temperature-dependent tai limitation
    trans_width = 2.0_wp

    ! adaptation factor to analysis interval for adaptive skin conductivity
    dtfac_skinc = (10800._wp/dt_ana)**(2._wp/3._wp)

    ! exclude the boundary interpolation zone of nested domains
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jt,ic,i_startidx,i_endidx,i_count,jc,ilu,t2mclim_hc,t_asyfac,tdiff_norm,wfac,trh_bias,skinc_fac,lat,scal)
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        & i_startidx, i_endidx, rl_start, rl_end)

      ! height-corrected climatological 2m-temperature
      DO jc = i_startidx, i_endidx
        ext_data%atm%t2m_clim_hc(jc,jb) = ext_data%atm%t2m_clim(jc,jb) + dtdz_clim * &
          (ext_data%atm%topography_c(jc,jb) - ext_data%atm%topo_t2mclim(jc,jb))

        t2mclim_hc(jc) = ext_data%atm%t2m_clim_hc(jc,jb) ! local copy needed for DA coupling (icpl_da_sfcevap > 0)
        t_asyfac(jc)   = SIGN(MIN(1._wp,ABS(ext_data%atm%t2m_climgrad(jc,jb))/2.5_wp), &
                                 -1._wp*ext_data%atm%t2m_climgrad(jc,jb))
      ENDDO

      IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
        DO jc = i_startidx, i_endidx
          t2mclim_hc(jc) = t2mclim_hc(jc) + 1.5_wp*SIGN(MIN(2.5_wp,ABS(nh_diag%t2m_bias(jc,jb))),nh_diag%t2m_bias(jc,jb))
        ENDDO
      ELSE IF (icpl_da_sfcevap >= 3) THEN
        DO jc = i_startidx, i_endidx
          t2mclim_hc(jc) = t2mclim_hc(jc) -                                                               &
            6._wp*SIGN(MIN(0.625_wp,ABS(10800._wp/dt_ana*nh_diag%t_avginc(jc,jb))),nh_diag%t_avginc(jc,jb))
        ENDDO
      ENDIF

      DO jt = 1, ntiles_total
        i_count = ext_data%atm%lp_count_t(jb,jt)
        IF (i_count == 0) CYCLE ! skip loop if the index list for the given tile is empty

        DO ic = 1, i_count
          jc = ext_data%atm%idx_lst_lp_t(ic,jb,jt)

          ilu = ext_data%atm%lc_class_t(jc,jb,jt)

          ! modification of tai/sai
          tdiff_norm = (t2mclim_hc(jc) - (threshold_temp(ilu)+t_asyfac(jc)*temp_asymmetry(ilu)) ) / trans_width
          tdiff_norm = MIN(1._wp,MAX(-1._wp,tdiff_norm))

          wfac = 0.5_wp*(1._wp + tdiff_norm) ! weighting factor of original ndvi-based tai specification

          ext_data%atm%tai_t(jc,jb,jt) = MIN(ext_data%atm%tai_t(jc,jb,jt), &
            wfac*ext_data%atm%tai_t(jc,jb,jt) + (1._wp-wfac)*laimin(ilu) )

          ext_data%atm%laifac_t(jc,jb,jt) = &
            (wfac*ext_data%atm%laimax_lcc(ilu) + (1._wp-wfac)*laimin(ilu))/MAX(0.01_wp,ext_data%atm%laimax_lcc(ilu))

          IF (lterra_urb) THEN
            ext_data%atm%sai_t(jc,jb,jt) = c_lnd * (1.0_wp - ext_data%atm%urb_isa_t(jc,jb,jt))                &
                                         + ext_data%atm%urb_ai_t(jc,jb,jt) * ext_data%atm%urb_isa_t(jc,jb,jt) &
                                         + ext_data%atm%tai_t(jc,jb,jt)
          ELSE
            ext_data%atm%sai_t(jc,jb,jt) = c_lnd + ext_data%atm%tai_t(jc,jb,jt)
          END IF


          ! modification of root depth
          ext_data%atm%rootdp_t(jc,jb,jt) = ext_data%atm%rootdmax_lcc(ilu) ! reset to table-based value

          IF (t2mclim_hc(jc) <= threshold_temp(ilu)-temp_asymmetry(ilu)) THEN
            ext_data%atm%rootdp_t(jc,jb,jt) = ext_data%atm%rootdp_t(jc,jb,jt)/rd_fac(ilu)
          ELSE IF (t2mclim_hc(jc) <= threshold_temp(ilu)+temp_asymmetry(ilu)+trans_width .AND. &
                   ext_data%atm%t2m_climgrad(jc,jb) > 0._wp) THEN
            wfac = (t2mclim_hc(jc)-(threshold_temp(ilu)-temp_asymmetry(ilu)))/(2._wp*temp_asymmetry(ilu)+trans_width)
            ext_data%atm%rootdp_t(jc,jb,jt) = ext_data%atm%rootdp_t(jc,jb,jt)*(wfac + (1._wp-wfac)/rd_fac(ilu))
          ELSE IF (t2mclim_hc(jc) <= threshold_temp(ilu)+trans_width .AND. &
                   ext_data%atm%t2m_climgrad(jc,jb) <= 0._wp) THEN
            wfac = (t2mclim_hc(jc)-(threshold_temp(ilu)-temp_asymmetry(ilu)))/(temp_asymmetry(ilu)+trans_width)
            ext_data%atm%rootdp_t(jc,jb,jt) = ext_data%atm%rootdp_t(jc,jb,jt)*(wfac + (1._wp-wfac)/rd_fac(ilu))
          ENDIF

          IF (icpl_da_sfcevap >= 3 .AND. icpl_da_skinc >= 1) THEN
            trh_bias = 10800._wp/dt_ana*(125._wp*nh_diag%rh_avginc(jc,jb) - &
              4._wp*(nh_diag%t_avginc(jc,jb)-0.5_wp*nh_diag%t_wgt_avginc(jc,jb)))
          ELSE IF (icpl_da_sfcevap >= 3) THEN
            trh_bias = 10800._wp/dt_ana*(100._wp*nh_diag%rh_avginc(jc,jb)-4._wp*nh_diag%t_avginc(jc,jb))
          ELSE IF (icpl_da_sfcevap >= 2) THEN
            trh_bias = nh_diag%t2m_bias(jc,jb) + 100._wp*10800._wp/dt_ana*nh_diag%rh_avginc(jc,jb)
          ELSE IF (icpl_da_sfcevap == 1) THEN
            trh_bias = nh_diag%t2m_bias(jc,jb)
          ELSE
            trh_bias = 0._wp
          ENDIF

          IF (icpl_da_sfcevap >= 4) THEN
            IF (trh_bias < 0._wp) THEN
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)*(1._wp-0.75_wp*trh_bias)
              ext_data%atm%r_bsmin(jc,jb)      = cr_bsmin*(1._wp-trh_bias)
            ELSE
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)/(1._wp+0.75_wp*trh_bias)
              ext_data%atm%r_bsmin(jc,jb)      = cr_bsmin/(1._wp+trh_bias)
            ENDIF
          ELSE IF (icpl_da_sfcevap >= 1) THEN
            IF (trh_bias < 0._wp) THEN
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)*(1._wp-0.5_wp*trh_bias)
              ext_data%atm%eai_t(jc,jb,jt)     = MERGE(c_soil_urb,c_soil,ilu == ext_data%atm%i_lc_urban) /     &
                                                 (1._wp-0.25_wp*trh_bias)
            ELSE
              ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)/(1._wp+0.5_wp*trh_bias)
              ext_data%atm%eai_t(jc,jb,jt)     = MIN(MERGE(c_soil_urb,c_soil,ilu == ext_data%atm%i_lc_urban) * &
                                                 (1._wp+0.25_wp*trh_bias), 2._wp)
            ENDIF

            IF (lterra_urb .AND. ((itype_eisa == 2) .OR. (itype_eisa == 3))) THEN
              ext_data%atm%eai_t(jc,jb,jt)     = ext_data%atm%eai_t(jc,jb,jt)                                  &
                                               * (1.0_wp - ext_data%atm%urb_isa_t(jc,jb,jt))
            END IF

          ELSE
            ext_data%atm%rsmin2d_t(jc,jb,jt) = ext_data%atm%stomresmin_lcc(ilu)
          ENDIF

          ! Tuning factor for skin conductivity
          IF (icpl_da_skinc >= 1) THEN
            scal = MERGE(4._wp, 2.5_wp, icpl_da_skinc == 1)
            IF (nh_diag%t_wgt_avginc(jc,jb) < 0._wp) THEN
              skinc_fac = MAX(0.1_wp,1._wp+dtfac_skinc*scal*nh_diag%t_wgt_avginc(jc,jb))
            ELSE
              skinc_fac = 1._wp/MAX(0.1_wp,1._wp-dtfac_skinc*scal*nh_diag%t_wgt_avginc(jc,jb))
            ENDIF

            lat = p_patch%cells%center(jc,jb)%lat*rad2deg
            IF (itype_lndtbl == 4 .AND. lat > -10._wp .AND. lat < 42.5_wp) THEN
              ext_data%atm%skinc_t(jc,jb,jt) = skinc_fac*MIN(200._wp,ext_data%atm%skinc_lcc(ilu)*          &
                                               (1._wp+MIN(1._wp,0.4_wp*(42.5_wp-lat),0.4_wp*(lat+10._wp))) )
            ELSE
              ext_data%atm%skinc_t(jc,jb,jt) = skinc_fac*ext_data%atm%skinc_lcc(ilu)
            ENDIF
          ENDIF

        ENDDO
      ENDDO

    ENDDO
!$OMP END DO
!$OMP END PARALLEL


  END SUBROUTINE vege_clim


  !-------------------------------------------------------------------------
  !! adjust atmo LSM to ocean LSM for coupled simulation and initialize new land points
  !!
  !! coupled A-O: ocean LSM dominates over atmospheric LSM: sea (0), land (1)
  !!   some fixes are done later in init_index_lists (also in mo_ext_data_init)
  !!
  !! ICON-O and ICON-A decision tree for LSM:
  !!
  !!                ICON-A (fr_land + fr_lake)
  !!                      0    0.6    1
  !!            ___________________________________
  !!             0    |   0   !0     !0       ICON-O deletes all ICON-A land
  !! ICON-O      0.4  |   0   !0.4   !0.4
  !! lsm_ctr_c   1    |  !1   !1      1
  !!
  !! note: lsm_ctr_c is real variable (float)
  !!       identical ocean/atmo grids: lsm_ctr_c is 0. or 1. only
  !!       different ocean/atmo grids: lsm_ctr_c is float (fraction)
  !! ICON-seamless prototype 2 - uncommon grids (A/O) support fractional lsm_ctr_c at ocean coast
  !! This routine doesn't support the case ntiles=1.  Additional surface parameters
  !! would have to be implemented
  !!
  !-------------------------------------------------------------------------

  SUBROUTINE lsm_ocean_atmo (p_patch, ext_data)

    TYPE(t_patch), INTENT(IN)            :: p_patch
    TYPE(t_external_data), INTENT(INOUT) :: ext_data

    INTEGER  :: jb,jc
    INTEGER  :: rl_start, rl_end
    INTEGER  :: i_startblk, i_endblk,i_startidx, i_endidx
    REAL(wp) :: lat, lon
    REAL(wp), PARAMETER :: eps = 1.e-10_wp

    !-----------------------------------------------------------------------

    CALL message('','Adjust atmo LSM to ocean LSM for coupled simulation.')

    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,lat,lon)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)

      DO jc = i_startidx,i_endidx

        ! land-sea-mask switched by ocean, 0: no change, 1: new land, 2: new ocean
        ext_data%atm%lsm_switch     (jc,jb)   = 0
 
        ! ICON-O is ocean but ICON-A is fractional land:
        ! convert ICON-A land and lake to pure ocean (flooded grid-point)

        IF ( ext_data%atm%lsm_ctr_c(jc,jb) == 0._wp ) THEN

          IF (( ext_data%atm%fr_land(jc,jb) + ext_data%atm%fr_lake(jc,jb)) > 0._wp ) THEN
            ext_data%atm%lsm_switch (jc,jb)   = 1
            ext_data%atm%fr_land    (jc,jb)   = 0._wp
            ext_data%atm%fr_lake    (jc,jb)   = 0._wp
            ext_data%atm%topography_c(jc,jb)  = 0._wp      ! topography to zero
          ENDIF

        ! ICON-O is land but ICON-A is fractional ocean:  make artificial ICON-A pure land
        ! convert ICON-A ocean to artificial pure land (rising grid-point): (alternative: artificial lake)

        ELSE IF ( ext_data%atm%lsm_ctr_c(jc,jb) > (1._wp-eps) ) THEN

          IF (( ext_data%atm%fr_land(jc,jb) + ext_data%atm%fr_lake(jc,jb)) < 1._wp ) THEN
            ! eps necessary because lsm_ctr_c has input values of 0.999999 instead of 1.0)

            ! switch land lakes from ocean to flake (Caspian and Dead Seas)

            lat = p_patch%cells%center(jc,jb)%lat*rad2deg
            lon = p_patch%cells%center(jc,jb)%lon*rad2deg

               ! Dead Sea    (< -390m)
            IF ( ( (lon > 35._wp).AND.(lon < 36._wp) .AND. (lat > 31._wp).AND.(lat < 32._wp) ) .OR. &
               ! Caspian Sea (< -25m)
               & ( (lon > 46._wp).AND.(lon < 55._wp) .AND. (lat > 36._wp).AND.(lat < 48._wp) ) ) THEN

              ext_data%atm%lsm_switch (jc,jb)   = 6
              ext_data%atm%fr_lake    (jc,jb)   = 1._wp - ext_data%atm%fr_land(jc,jb)
              ext_data%atm%depth_lk   (jc,jb)   = 50._wp

            ELSE

              ext_data%atm%lsm_switch (jc,jb)   = 2
              ! set fraction to land
              ext_data%atm%fr_land    (jc,jb)   = 1._wp
              ext_data%atm%fr_lake    (jc,jb)   = 0._wp
              ! topography values
              ext_data%atm%topography_c(jc,jb)  = 0._wp       ! alternative: add little hight (1m?)
              ext_data%atm%soiltyp    (jc,jb)   = 4           ! soil type to sandy loam (sfc_terra_data)
              ext_data%atm%lu_class_fraction(jc,jb,:)                       = 0._wp
              ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_grass) = 1._wp 
                                                              ! grass-land (read_ext_data_atm)
              ! land use variables will then automatically be set in init_index_lists:
              !   frac_t, lc_frac_t, fr_glac, lc_class_t, llsm_atm_c, llake_c,
              !   plcov_mx, lai_mx, rootdp, skinc, rsmin, z0
              ! initialization of soil moisture (w_so) and temperature (t_so) in new_land_from_ocean

            ENDIF

          ENDIF

        ! ICON-O is fractional ocean (coast)
        ! attention: fr_lake is not allowed on the coast

        ELSE IF ( ext_data%atm%lsm_ctr_c(jc,jb) > 0._wp ) THEN

          ! ICON-A is pure land: reduce land & delete lake accordingly
          IF (( ext_data%atm%fr_land  (jc,jb) + ext_data%atm%fr_lake  (jc,jb)) == 1._wp ) THEN

            ext_data%atm%lsm_switch (jc,jb)  = 3
           !ext_data%atm%fr_land    (jc,jb)  = ext_data%atm%fr_land(jc,jb) * ext_data%atm%lsm_ctr_c(jc,jb)
           !ext_data%atm%fr_lake    (jc,jb)  = ext_data%atm%fr_lake(jc,jb) * ext_data%atm%lsm_ctr_c(jc,jb)
            ext_data%atm%fr_land    (jc,jb)  = ext_data%atm%lsm_ctr_c(jc,jb)
            ext_data%atm%fr_lake    (jc,jb)  = 0._wp

          ! ICON-A is fraction: adjust ICON-A lsm to ICON-O lsm, no lake fraction
          !   attention: prevent small fr_land being deleted with soiltyp=9
          ELSE IF ( (ext_data%atm%fr_land(jc,jb) > frlnd_thrhld)    .OR. &
            &       (ext_data%atm%fr_lake(jc,jb) > frlake_thrhld) ) THEN

            ext_data%atm%lsm_switch (jc,jb)  = 4
           !lsm_atm                          = ext_data%atm%fr_land(jc,jb) + ext_data%atm%fr_lake  (jc,jb)
           !ext_data%atm%fr_land    (jc,jb)  = ext_data%atm%fr_land(jc,jb) * ext_data%atm%lsm_ctr_c(jc,jb) / lsm_atm
           !ext_data%atm%fr_lake    (jc,jb)  = ext_data%atm%fr_lake(jc,jb) * ext_data%atm%lsm_ctr_c(jc,jb) / lsm_atm
            ext_data%atm%fr_land    (jc,jb)  = ext_data%atm%lsm_ctr_c(jc,jb)
            ext_data%atm%fr_lake    (jc,jb)  = 0._wp

          ! ICON-A is pure ocean: add land in ICON-A according to ICON-O lsm, no lake fraction (rising coastal area)
          !   attention: small fr_land is deleted elsewhere
          ELSE IF ( (ext_data%atm%fr_land(jc,jb)) <= frlnd_thrhld    .AND. &
            &       (ext_data%atm%fr_lake(jc,jb)) <= frlake_thrhld ) THEN

            ext_data%atm%lsm_switch (jc,jb)  = 5
            ! set fraction to land
            ext_data%atm%fr_land    (jc,jb)  = ext_data%atm%lsm_ctr_c(jc,jb)
            ext_data%atm%fr_lake    (jc,jb)  = 0._wp
            ! topography values
            ext_data%atm%topography_c(jc,jb) = 0._wp       ! alternative: add little hight (1m?)
            ext_data%atm%soiltyp    (jc,jb)  = 4            ! soil type to sandy loam (sfc_terra_data)
            ext_data%atm%lu_class_fraction(jc,jb,:)                       = 0._wp
            ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_grass) = 1._wp  ! grass-land (read_ext_data_atm)
            ! land use variables will then automatically be set in init_index_lists:
            !   frac_t, lc_frac_t, fr_glac, lc_class_t, llsm_atm_c, llake_c,
            !   plcov_mx, lai_mx, rootdp, skinc, rsmin, z0
            ! initialization of soil moisture (w_so) and temperature (t_so) in new_land_from_ocean

          ENDIF

        ! ICON-O is < 0.0 or > 1.0; not allowed (ICON-O input data problems)
        ELSE   
          ext_data%atm%lsm_switch   (jc,jb)  = 10
        ENDIF

      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

  END SUBROUTINE lsm_ocean_atmo


END MODULE mo_ext_data_init

