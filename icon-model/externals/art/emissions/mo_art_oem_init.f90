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

MODULE mo_art_oem_init

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_var_list,                ONLY: t_var_list_ptr
  USE mo_tracer_metadata_types,   ONLY: t_chem_meta
  USE mo_model_domain,            ONLY: p_patch
  USE mo_parallel_config,         ONLY: idx_1d

  USE mo_var_metadata_types,      ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_physical_constants,      ONLY: grav, rd, p0sl_bg
  USE mo_nh_vert_interp,          ONLY: prepare_extrap, prepare_lin_intp, lin_intp
  USE mo_nh_vert_interp_ipz,      ONLY: z_at_plevels
  USE mo_impl_constants,          ONLY: min_rlcell
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_ext_data_state,          ONLY: ext_data
  USE mo_ext_data_types,          ONLY: t_external_data

  USE mo_time_config,             ONLY: time_config
  USE mtime,                      ONLY: MAX_DATETIME_STR_LEN, &
                                    &   datetimeToString,     &
                                    &   getnoofdaysinyeardatetime, &
                                    &   datetime,             &
                                    &   newDatetime,          &
                                    &   OPERATOR(==)

  !ART
  USE mo_art_data,                ONLY: p_art_data
  USE mo_art_atmo_data,           ONLY: t_art_atmo
  USE mo_art_wrapper_routines,    ONLY: art_get_indices_c

  !OEM
  USE mo_art_oem_ensemble,        ONLY: art_oem_ens_init
  USE mo_oem_config,              ONLY: gridded_emissions_nc,              &
                                    &   vertical_profile_nc,               &
                                    &   hour_of_day_nc,                    &
                                    &   day_of_week_nc,                    &
                                    &   month_of_year_nc,                  &
                                    &   hour_of_year_nc,                   &
                                    &   vegetation_indices_nc
  USE mo_art_oem_types,           ONLY: p_art_oem_data,                    &
                                    &   t_art_oem_data,                    &
                                    &   t_art_oem_config,                  &
                                    &   t_art_oem_ensemble

  USE netcdf,                     ONLY: nf90_noerr, nf90_open,             &
                                    &   nf90_nowrite, nf90_strerror,       & 
                                    &   nf90_inq_dimid, nf90_inq_varid,    &
                                    &   nf90_inquire_dimension,            &
                                    &   nf90_get_var, nf90_close


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_oem_init, art_oem_init_and_read_gridded_emissions, &
    &       art_oem_init_vertical_profile_fields, &
    &       art_oem_init_temporal_profile_fields, &
    &       art_oem_read_temporal_profile_from_file, &
    &       art_oem_interpolate_vertical_profiles, &
    &       art_oem_read_vertical_profile_from_file, &
    &       art_oem_cleanup_vertical_profile_fields, &
    &       art_oem_init_and_read_vegetation_indices

  !---------------------------------------------------------------
  ! Pointer to OEM data structures

  TYPE(t_art_atmo),         POINTER :: art_atmo     !< ART atmo fields
  TYPE(t_art_oem_data),     POINTER :: oem_data     !< OEM data structure -> data
  TYPE(t_art_oem_config),   POINTER :: oem_config   !< OEM data structure -> config
  TYPE(t_art_oem_ensemble), POINTER :: oem_ensemble !< OEM data structure -> ensemble

  !---------------------------------------------------------------
  ! Vertical profile variables

  INTEGER :: vp_ncategory, & ! Number of categories for vertical profiles
    &        vp_nlevel       ! Number of levels for vertical profiles

  REAL(KIND=wp),  DIMENSION(:), ALLOCATABLE :: &
    & vp_layer_bot,        & ! bottom of layer above ground
    & vp_layer_top           ! top of layer above ground

  REAL(KIND=wp),  DIMENSION(:,:), ALLOCATABLE :: &
    & vp_factor                  ! vertical scale factor 

  !------------------------------------------------------------------
  ! Temporal profile variables

  INTEGER :: tp_ncategory ! Number of categories for temporal profiles

  ! Constant variable
  INTEGER,  PARAMETER :: ncat_max = 200
  INTEGER,  PARAMETER :: tp_param_hourofday = 24
  INTEGER,  PARAMETER :: tp_param_dayofweek = 7
  INTEGER,  PARAMETER :: tp_param_monthofyear = 12
  INTEGER,  PARAMETER :: tp_param_hour = 8784

  CHARACTER(LEN=3), DIMENSION(7) :: day_of_week = &
    & (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)

  REAL(KIND=wp), PARAMETER :: h_scal_bg = 10000._wp  ! [m]    scale height
  REAL(KIND=wp), PARAMETER :: t0sl_bg   = 288.15_wp  ! [K]    sea level temperature
  REAL(KIND=wp), PARAMETER :: del_t_bg  = 75._wp     ! [K]    difference between sea level
  !                                                           temperature and asymptotic
  !                                                           stratospheric temperature


!------------------------------------------------------------------------------------

CONTAINS


!
!------------------------------------------------------------------------------------
!

SUBROUTINE art_oem_init(jg, p_prog_list, ierror, yerrmsg)

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg           !< patch on which computation is performed
  TYPE(t_var_list_ptr), INTENT(IN)   :: p_prog_list  !< list of prognostic variables
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER(LEN= *), INTENT(INOUT)   :: yerrmsg

  !--------------------------------------------------------------------
  ! local variables

  CHARACTER(:), ALLOCATABLE :: oem_type_string, oem_cat_string, oem_vp_string, &
     &                         oem_tp_string, oem_init_string, oem_ftype_string

  CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE :: ycatl_temp, yvpl_temp, ytpl_temp, ens_name

  CHARACTER(LEN=20) :: tracer_name

  INTEGER :: n_cat, iemis, ivprm, ic, n, ens_nr, ierr, iv, irs, oem_itscale

  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element

  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element



  CHARACTER(*), PARAMETER :: routine = "art_oem_init"


!------------------------------------------------------------------------------
! Section 1: Initialize config-vars and get number and types of OEM-tracers
!------------------------------------------------------------------------------

  ierr = 0
  ierror = 0

  art_atmo => p_art_data(jg)%atmo

  oem_data => p_art_oem_data%data_fields
  oem_config => p_art_oem_data%configure
  oem_ensemble => p_art_oem_data%ensemble

  ! Auxillary variable to read ensemble names
  ALLOCATE(ens_name(2)); ens_name=''

  ! Set initial value to config-variables
  oem_config%emis_tracer = 0
  oem_config%ens_tracer = 0
  oem_config%vprm_tracer = 0
  oem_ensemble%ens_name = ''

  DO iv = 1, p_prog_list%p%nvars
    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info
    IF (info_dyn%tracer%lis_tracer) THEN
      SELECT TYPE(meta => info_dyn%tracer)
        CLASS IS(t_chem_meta)

          ! Check for 'oem_type', which labels a tracer as "OEM-tracer"
          CALL meta%opt_meta%get("oem_type",oem_type_string,ierr)
          IF (ierr/=0) THEN
            CYCLE
          ENDIF

          SELECT CASE(oem_type_string)
            CASE("emis") 
              ! Count OEM-tracers with emissions
              oem_config%emis_tracer = oem_config%emis_tracer + 1
            CASE("vprm")
              ! Count OEM-tracers with VPRM fluxes
              oem_config%vprm_tracer = oem_config%vprm_tracer + 1
            CASE("ens")
              ! Count ensemble tracers
              oem_config%ens_tracer = oem_config%ens_tracer + 1
              ! Read the tracer-name and number of ensemble-member
              ! (e.g. 'CH4_A-01' --> 'CH4_A' + '1')
              tracer_name=TRIM(meta%name)
              DO ic=1,len(tracer_name)
                IF (tracer_name(ic:ic) == "-") tracer_name(ic:ic) = " "
              ENDDO
              READ(tracer_name,fmt=*) ens_name(:)
              READ(ens_name(2),*) ens_nr
              ! The ens_table contains:
              !   1) the number of the ensemble member (e.g. '1')
              !   2) the tracer-index in the ICON data structure
              ! Additionally, the name of the reference-tracer (e.g. 'CH4_A')
              ! is stored in ens_name.
              ! With this information, the emissions can be added to the correct tracer
              DO n=1,ncat_max
                IF (oem_ensemble%ens_name(n)=='') THEN
                  oem_ensemble%ens_name(n) = ens_name(1)
                  oem_ensemble%ens_table(1,n) = ens_nr
                  oem_ensemble%ens_table(2,n) = info%ncontained
                  EXIT
                ENDIF
              ENDDO !n=1,ncat_max
          END SELECT !oem_type_string=="emis"/"vprm"/"ens"
      END SELECT !SELECT TYPE(meta => info_dyn%tracer)
    ENDIF !info_dyn%tracer%lis_tracer
  END DO


!------------------------------------------------------------------------------
! Section 2: Now the number of the different OEM-tracers are known, allocate the variables
!------------------------------------------------------------------------------

  IF (oem_config%emis_tracer>0) THEN
    ALLOCATE(oem_config%ycatl_l(oem_config%emis_tracer,ncat_max)); oem_config%ycatl_l=''
    ALLOCATE(oem_config%yvpl_l(oem_config%emis_tracer,ncat_max)); oem_config%yvpl_l=''
    ALLOCATE(oem_config%ytpl_l(oem_config%emis_tracer,ncat_max)); oem_config%ytpl_l=''
    ALLOCATE(oem_config%emis_name(oem_config%emis_tracer)); oem_config%emis_name=''
    ALLOCATE(oem_config%itype_tscale_l(oem_config%emis_tracer)); oem_config%itype_tscale_l = 0
    ALLOCATE(oem_config%emis_idx(oem_config%emis_tracer)); oem_config%emis_idx = 0
    ALLOCATE(oem_config%vprm_idx(oem_config%vprm_tracer)); oem_config%vprm_idx = 0
  ENDIF

  IF (oem_config%vprm_tracer>0) THEN
    ALLOCATE(oem_config%vprm_name(oem_config%vprm_tracer)); oem_config%vprm_name=''
    ALLOCATE(oem_config%vprm_flux_type(oem_config%vprm_tracer)); oem_config%vprm_flux_type=''
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Read categories and properties of emission, background and restart tracer
!------------------------------------------------------------------------------

  iemis = 1
  ivprm = 1
  irs = 1

  DO iv = 1, p_prog_list%p%nvars
    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info
    IF (info_dyn%tracer%lis_tracer) THEN
      SELECT TYPE(meta => info_dyn%tracer)
        CLASS IS(t_chem_meta)

          ! Get the type of OEM-tracer
          CALL meta%opt_meta%get("oem_type",oem_type_string,ierr)
          IF (ierr==0) THEN

            ! For emission-tracers, all the emission categories and vertical-
            ! and temporal-profiles need to be read:
            IF (oem_type_string=="emis") THEN
              oem_config%emis_name(iemis) = TRIM(meta%name)
              oem_config%emis_idx(iemis) = info%ncontained  ! index in 4D tracer container
              ! Tracer catepories
              CALL meta%opt_meta%get("oem_cat",oem_cat_string,ierr)
              IF (ierr/=0) THEN
                CALL finish(routine,'oem_cat is missing for the OEM tracer')
              ENDIF
              ! Vertical profiles
              CALL meta%opt_meta%get("oem_vp",oem_vp_string,ierr)
              IF (ierr/=0) THEN
                CALL finish(routine,'oem_vp is missing for the OEM tracer')
              ENDIF
              ! Temporal profiles
              CALL meta%opt_meta%get("oem_tp",oem_tp_string,ierr)
              IF (ierr/=0) THEN
                CALL finish(routine,'oem_tp is missing for the OEM tracer')
              ENDIF
              ! Type of time profile (hod or hoy)
              CALL meta%opt_meta%get("oem_tscale",oem_itscale,ierr)
              IF (ierr/=0) THEN
                CALL finish(routine,'oem_tscale is missing for the OEM tracer')
              ENDIF
              ! set itype_tscale_l
              IF (oem_itscale==1) THEN
                oem_config%itype_tscale_l(iemis)=1
              ENDIF
              IF (oem_itscale==2) THEN
                oem_config%itype_tscale_l(iemis)=2
              ENDIF
              ! Read emission categories
              ! Determine the number of categories per tracer
              n_cat = 1
              DO ic=1,len(oem_cat_string)
                if (oem_cat_string(ic:ic) == ',') n_cat = n_cat + 1
              ENDDO
              ! Read OEM categories into temporary variable
              ALLOCATE(ycatl_temp(n_cat))
              READ(UNIT=oem_cat_string,fmt=*) ycatl_temp(:)
              ! Write categories into config-variable
              DO n=1,n_cat
                oem_config%ycatl_l(iemis,n)=ycatl_temp(n)
              ENDDO
              DEALLOCATE(ycatl_temp)
              ! Read vertical categories
              ! Determine the number of categories per tracer
              n_cat = 1
              DO ic=1,len(oem_vp_string)
                if (oem_vp_string(ic:ic) == ',') n_cat = n_cat + 1
              ENDDO
              ! Read OEM categories into temporary variable
              ALLOCATE(yvpl_temp(n_cat))
              READ(UNIT=oem_vp_string,fmt=*) yvpl_temp(:)
              ! Write categories into config-variable
              DO n=1,n_cat
                oem_config%yvpl_l(iemis,n)=yvpl_temp(n)
              ENDDO
              DEALLOCATE(yvpl_temp)
              ! Read temporal categories
              ! Determine the number of categories per tracer
              n_cat = 1
              DO ic=1,len(oem_tp_string)
                if (oem_tp_string(ic:ic) == ',') n_cat = n_cat + 1
              ENDDO
              ! Read OEM categories into temporary variable
              ALLOCATE(ytpl_temp(n_cat))
              READ(UNIT=oem_tp_string,fmt=*) ytpl_temp(:)
              ! Write categories into config-variable
              DO n=1,n_cat
                oem_config%ytpl_l(iemis,n)=ytpl_temp(n)
              ENDDO
              DEALLOCATE(ytpl_temp)
              iemis = iemis + 1
            ELSEIF (oem_type_string=="vprm") THEN
              oem_config%vprm_name(ivprm) = TRIM(meta%name)
              oem_config%vprm_idx(ivprm) = info%ncontained  ! index in 4D tracer container
              ! Tracer flux type
              CALL meta%opt_meta%get("oem_ftype", oem_ftype_string, ierr)
              IF (ierr/=0) THEN
                CALL finish(routine, &
                   &        'oem_ftype ("resp" or "gpp" is missing for the OEM tracer')
              ENDIF
              ! Write the flux type into config-variable
              oem_config%vprm_flux_type(ivprm)=TRIM(oem_ftype_string)
              ivprm = ivprm + 1
            ENDIF !oem_type_string=="vprm"
          ENDIF !ierr==0
      END SELECT !SELECT TYPE(meta => info_dyn%tracer)
    ENDIF !info_dyn%tracer%lis_tracer
  END DO


!------------------------------------------------------------------------------
! Section 4: Calling the further initialitation routines
!------------------------------------------------------------------------------

  IF (oem_config%emis_tracer>0) THEN
    CALL art_oem_init_vertical_profile_fields(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Vertical profiles initialized.')

    CALL art_oem_init_temporal_profile_fields(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Temporal profiles initialized.')

    CALL art_oem_init_and_read_gridded_emissions(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Gridded emissions initialized and read.')

    CALL art_oem_read_temporal_profile_from_file(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Temporal profiles read.')

    CALL art_oem_read_vertical_profile_from_file(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Vertical profiles read.')

    CALL art_oem_interpolate_vertical_profiles(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Vertical profiles interpolated.')
  ENDIF

  IF (oem_config%ens_tracer>0) THEN
    CALL art_oem_ens_init(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Ensemble scaling factors initialized')
  ENDIF

  IF (oem_config%vprm_tracer>0) THEN
    CALL art_oem_init_and_read_vegetation_indices(jg, ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Vegetation indices initialized and read.')

    CALL art_oem_init_vprm_lu_class_fraction(jg, ext_data(jg), ierror, yerrmsg)
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'VPRM land-use class fractions initialized.')
  ENDIF


!------------------------------------------------------------------------------
! Section 5: Calling cleanup routines
!------------------------------------------------------------------------------

  IF (oem_config%emis_tracer>0) THEN
    CALL art_oem_cleanup_vertical_profile_fields()
    IF (ierror /= 0) THEN
      CALL finish (TRIM(routine), yerrmsg)
    ENDIF
    CALL message(routine,'Vertical profile fields cleaned up')
  ENDIF


!------------------------------------------------------------------------------
! End of Subroutine art_oem_init
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_init



!==============================================================================
!==============================================================================
!+ Initialize and read gridded emissions from NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE art_oem_init_and_read_gridded_emissions(jg, ierror, yerrmsg)

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg         !< patch on which computation is performed
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER(LEN= *), INTENT(INOUT)   :: yerrmsg

  !---------------------------------------------------------------------
  ! Local variables

  REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
    & gridded_emissions_tot       ! Annual mean "2D" emissions fields (ncells, ntracers)

  INTEGER, DIMENSION(:), ALLOCATABLE :: &
    & country_ids_tot             ! EMEP country code on the whole domain (ncells)

  CHARACTER(LEN=100), DIMENSION(600) :: unique_cat

  LOGICAL :: file_exists

  INTEGER :: i, n, m, ncid, dimid, n_glb_cells, varid, jc, nc, jb, &
    &        i_startblk, i_endblk, glbidx, is, ie, locidx, ge_ncategory

  CHARACTER(*), PARAMETER :: routine = "art_oem_init_and_read_gridded_emissions"


!------------------------------------------------------------------------------
! Section 1: Get the unique category names
!------------------------------------------------------------------------------

  yerrmsg = '   '
  ierror = nf90_noerr

  ! Check if the file exists
  INQUIRE(FILE=gridded_emissions_nc,EXIST=file_exists)
  IF (.NOT.(file_exists)) THEN
    ierror=1
    yerrmsg = "OAE gridded emissions file does not exist! "
    RETURN
  ENDIF

  ! Find the unique category names
  ge_ncategory = 1
  unique_cat = ''
  DO i = 1, oem_config%emis_tracer
    inner: DO n = 1, ncat_max
      IF (oem_config%ycatl_l(i,n) /= '') THEN
        DO m = 1,ge_ncategory
          IF (unique_cat(1) == '') THEN
            unique_cat(1) = oem_config%ycatl_l(i,n)
            CYCLE inner
          ENDIF
          IF(unique_cat(m) == oem_config%ycatl_l(i,n)) THEN
            CYCLE inner
          ENDIF
        ENDDO
        ge_ncategory = ge_ncategory + 1
        unique_cat(ge_ncategory) = oem_config%ycatl_l(i,n)
      ENDIF
    ENDDO inner
  ENDDO

!------------------------------------------------------------------------------
! Section 2: Read-in the gridded emissions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 2.1: Get the dimensions
!------------------------------------------------------------------------------

  ! Open the NetCDF file
  ierror = nf90_open(gridded_emissions_nc, nf90_nowrite, ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Get number of global cells
  ierror = nf90_inq_dimid(ncid, 'cell', dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_inquire_dimension(ncid, dimid, len=n_glb_cells)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 2.2: Allocate variables
!------------------------------------------------------------------------------

  ALLOCATE(gridded_emissions_tot(ge_ncategory,n_glb_cells)); gridded_emissions_tot = 0.0_wp
  ALLOCATE(country_ids_tot(n_glb_cells)); country_ids_tot = 0
  ALLOCATE(oem_config%gridded_emissions_idx(ge_ncategory)); oem_config%gridded_emissions_idx = ''
  ALLOCATE(oem_data%country_ids(art_atmo%nproma,art_atmo%nblks)); oem_data%country_ids = 0
  ALLOCATE(oem_data%gridded_emissions(art_atmo%nproma,art_atmo%nblks,ge_ncategory)); oem_data%gridded_emissions = 0.0_wp

  DO i = 1, ge_ncategory
    oem_config%gridded_emissions_idx(i) = ''
  ENDDO
  DO i = 1, ge_ncategory
    oem_config%gridded_emissions_idx(i) = unique_cat(i)
  ENDDO

!------------------------------------------------------------------------------
! Section 2.3: Read emissions
!------------------------------------------------------------------------------

  DO i = 1, ge_ncategory
    ierror = nf90_inq_varid(ncid, oem_config%gridded_emissions_idx(i), varid )
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_get_var(ncid, varid, &
       &     gridded_emissions_tot(i,:), start = (/1/), count=(/n_glb_cells/))
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  ENDDO

  ! Read the country_ids
  ierror = nf90_inq_varid(ncid, "country_ids", varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_get_var(ncid, varid, country_ids_tot, start = (/1/), count=(/n_glb_cells/))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Close the NetCDF file
  ierror = nf90_close(ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 3: Distribute the values
!------------------------------------------------------------------------------

  i_startblk = art_atmo%i_startblk
  i_endblk = art_atmo%i_endblk

  ! distribute emission values
  DO nc = 1, ge_ncategory
    DO jb = i_startblk, i_endblk
      ! read indices within this block:
      CALL art_get_indices_c(jg, jb, is, ie)
      DO jc = is, ie
        ! get local index from jc and jb:
        locidx = idx_1d(jc,jb)
        ! get global index from local index:
        glbidx = p_patch(jg)%cells%decomp_info%glb_index(locidx)
        oem_data%gridded_emissions(jc,jb,nc) = gridded_emissions_tot(nc,glbidx)
      ENDDO
    ENDDO
  ENDDO

  ! distribute coundry indices
  DO jb = i_startblk, i_endblk
    CALL art_get_indices_c(jg, jb, is, ie)
    DO jc = is, ie
      locidx = idx_1d(jc,jb)
      glbidx = p_patch(jg)%cells%decomp_info%glb_index(locidx)
      oem_data%country_ids(jc,jb) = country_ids_tot(glbidx)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_init_and_read_gridded_emissions



!==============================================================================
!==============================================================================
!+ Allocate the data structure necessary for the vertical profile
!------------------------------------------------------------------------------

SUBROUTINE art_oem_init_vertical_profile_fields(jg, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine reads in the parameters of the vertical scaling profiles
!   for OEM. It also allocates and initializes the vertical scaling variables.  
!
!------------------------------------------------------------------------------
   IMPLICIT NONE

   INTEGER, INTENT(IN)                :: jg         !< patch on which computation is performed
   INTEGER, INTENT(INOUT)             :: ierror
   CHARACTER (LEN= *),INTENT(INOUT)   :: yerrmsg

   !---------------------------------------------
   ! Local variables:

   INTEGER :: m, n, ncid, dimid, i
   CHARACTER (LEN=100), DIMENSION(600) :: unique_cat
   LOGICAL :: file_exists

   CHARACTER(*), PARAMETER :: routine = "art_oem_init_vertical_profile_fields"

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Open the NetCDF file and read in the parameters
!------------------------------------------------------------------------------

  ierror = nf90_noerr
  ! Check if the file exists
  INQUIRE(FILE=vertical_profile_nc,EXIST=file_exists)
  IF (.NOT.(file_exists)) THEN
    ierror = 1
    yerrmsg = "OAE vertical profile file does not exist! "
    RETURN
  ENDIF

  ! Open the vertical profile NetCDF file
  ierror = nf90_open(vertical_profile_nc, nf90_nowrite, ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Get level dimension information
  ierror = nf90_inq_dimid(ncid, "level", dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, dimid, len = vp_nlevel)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Find the unique category names
  vp_ncategory = 1
  unique_cat = ''
  DO i= 1,oem_config%emis_tracer
    inner: DO n = 1, ncat_max
      IF (oem_config%yvpl_l(i,n) /= '') THEN
        DO m = 1,vp_ncategory
          IF (unique_cat(1) == '') THEN
            unique_cat(1) = oem_config%yvpl_l(i,n)
            CYCLE inner
          ENDIF
          IF(unique_cat(m) == oem_config%yvpl_l(i,n)) THEN
            CYCLE inner
          ENDIF
        ENDDO
        vp_ncategory = vp_ncategory + 1
        unique_cat(vp_ncategory) = oem_config%yvpl_l(i,n)
      ENDIF
    ENDDO inner
  ENDDO

  ! close the file
  ierror = nf90_close(ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

!-------------------------------------------------------------------
! Section 2: Allocate and initialize the vertical profile variables
!-------------------------------------------------------------------

  ALLOCATE(vp_layer_bot(vp_nlevel)); vp_layer_bot = 0.0_wp
  ALLOCATE(vp_layer_top(vp_nlevel)); vp_layer_top = 0.0_wp
  ALLOCATE(vp_factor(vp_ncategory,vp_nlevel)); vp_factor = 0.0_wp
  ALLOCATE(oem_config%vp_category(vp_ncategory)); oem_config%vp_category = ''
  ALLOCATE(oem_data%vert_scaling_fact(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks,vp_ncategory))
  oem_data%vert_scaling_fact = 0.0_wp

  ! Save the vertical profile category names
  DO n = 1,vp_ncategory
    oem_config%vp_category(n) = unique_cat(n)
  ENDDO


!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_init_vertical_profile_fields



!==============================================================================
!==============================================================================
!+ Allocate the data structures necessary for the temporal profiles
!------------------------------------------------------------------------------

SUBROUTINE art_oem_init_temporal_profile_fields(jg, ierror, yerrmsg)

!------------------------------------------------------------------------------
! Description: This subroutine reads in the parameters for the temporal profiles
!              and allocates and initializes the necessary variables.
!-----------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg         !< patch on which computation is performed
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER (LEN= *), INTENT(INOUT)  :: yerrmsg

  !--------------------------------------------------
  ! Local variables

  INTEGER :: ncid, dimid, tp_nvar, i, n, m
  CHARACTER (LEN=100), DIMENSION(600) :: unique_cat
  LOGICAL :: file_exists

  CHARACTER(*), PARAMETER :: routine = "art_oem_init_temporal_profile_fields"

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Open NetCDF file and read in variables
!------------------------------------------------------------------------------

  ierror = nf90_noerr
  yerrmsg='     '

!------------------------------------------------------------------------------
! Section 1.1: Check if all files exist and get dimensions for hour_of_day
!------------------------------------------------------------------------------

  IF ( ANY(oem_config%itype_tscale_l(:)==1) ) THEN
    INQUIRE(FILE=hour_of_day_nc,EXIST=file_exists)
    IF (.NOT.(file_exists)) THEN
      ierror=1
      yerrmsg = "OEM hour_of_day file does not exist! "
      RETURN
    ENDIF

    INQUIRE(FILE=day_of_week_nc,EXIST=file_exists)
    IF (.NOT.(file_exists)) THEN
      ierror=1
      yerrmsg = "OEM day_of_week file does not exist! "
      RETURN
    ENDIF

    INQUIRE(FILE=month_of_year_nc,EXIST=file_exists)
    IF (.NOT.(file_exists)) THEN
      ierror=1
      yerrmsg = "OEM month_of_year file does not exist! "
      RETURN
    ENDIF

    ! Open the NetCDF file
    ierror = nf90_open(month_of_year_nc, nf90_nowrite, ncid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ! Get countryID dimension information
    ierror = nf90_inq_dimid(ncid, "country", dimid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ierror = nf90_inquire_dimension(ncid, dimid, len = oem_config%tp_ncountry)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
    ierror = nf90_close(ncid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF ! hour of day

!------------------------------------------------------------------------------
! Section 1.2: Check if all files exist and get dimensions for hour_of_year
!------------------------------------------------------------------------------

  IF ( ANY(oem_config%itype_tscale_l(:)==2) ) THEN
    INQUIRE(FILE=hour_of_year_nc,EXIST=file_exists)
    IF (.NOT.(file_exists)) THEN
      ierror=1
      yerrmsg = "OEM hour_of_year file does not exist! "
      RETURN
    ENDIF

    ! Open the NetCDF file
    IF ( .NOT. ANY(oem_config%itype_tscale_l(:)==1) ) THEN
      ierror = nf90_open(hour_of_year_nc, nf90_nowrite, ncid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
      ! Get countryID dimension information
      ierror = nf90_inq_dimid(ncid, "country", dimid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_inquire_dimension(ncid, dimid, len = oem_config%tp_ncountry)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_close(ncid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF ! .NOT. hour of day
  ENDIF ! hour of year

!------------------------------------------------------------------------------
! Section 2: Find the unique category names
!------------------------------------------------------------------------------

  tp_ncategory = 1
  unique_cat = ''
  DO i= 1, oem_config%emis_tracer
    inner: DO n = 1, ncat_max
      IF (oem_config%ytpl_l(i,n) /= '') THEN
        DO m = 1,tp_ncategory
          IF (unique_cat(1) == '') THEN
            unique_cat(1) = oem_config%ytpl_l(i,n)
            CYCLE inner
          ENDIF
          IF(unique_cat(m) == oem_config%ytpl_l(i,n)) THEN
            CYCLE inner
          ENDIF
        ENDDO
        tp_ncategory = tp_ncategory + 1
        unique_cat(tp_ncategory) = oem_config%ytpl_l(i,n)
      ENDIF
    ENDDO inner
  ENDDO

!-------------------------------------------------------------------
! Section 3: Allocate and initialize the temporal profile variables
!-------------------------------------------------------------------

  ! Allocate fields to store temporal profile information
  ALLOCATE(oem_config%tp_category(tp_ncategory)); oem_config%tp_category = ''
  ALLOCATE(oem_data%tp_countryid(oem_config%tp_ncountry)); oem_data%tp_countryid = 0

  IF ( ANY(oem_config%itype_tscale_l(:)==1) ) THEN
    ALLOCATE(oem_data%tp_hourofday(tp_param_hourofday, tp_ncategory,oem_config%tp_ncountry))
    oem_data%tp_hourofday = 0.0_wp
    ALLOCATE(oem_data%tp_dayofweek(tp_param_dayofweek, tp_ncategory,oem_config%tp_ncountry))
    oem_data%tp_dayofweek = 0.0_wp
    ALLOCATE(oem_data%tp_monthofyear(tp_param_monthofyear, tp_ncategory,oem_config%tp_ncountry))
    oem_data%tp_monthofyear = 0.0_wp
  ENDIF

  IF ( ANY(oem_config%itype_tscale_l(:)==2) ) THEN
    ALLOCATE(oem_data%tp_hourofyear(tp_param_hour, tp_ncategory, oem_config%tp_ncountry))
    oem_data%tp_hourofyear = 0.0_wp
  ENDIF

  ! Save the category names
  DO n = 1, tp_ncategory
    oem_config%tp_category(n) = unique_cat(n)
  ENDDO

 
!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_init_temporal_profile_fields



!==============================================================================
!==============================================================================
!+ Read in the temporal profile variables
!------------------------------------------------------------------------------

SUBROUTINE art_oem_read_temporal_profile_from_file(jg, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:  This subroutine reads in the temporal profile variables
!   from multiple NetCDF files and distributes them to all processors.
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg         !< patch on which computation is performed
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER (LEN= *),INTENT(INOUT)   :: yerrmsg

  !----------------------------------------------------------------
  ! Local variables

  INTEGER :: i, hdid, dwid, myid, hyid, n, varid, tp_nvar, tp_id
  CHARACTER (LEN=20) :: var_name, tp_cat

  CHARACTER(*), PARAMETER :: routine = "art_oem_read_temporal_profile_from_file"

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Open NetCDF file and read in variables
!------------------------------------------------------------------------------

  ierror = nf90_noerr
  yerrmsg = '   '

!------------------------------------------------------------------------------
! Section 1.1: For hour_of_day case
!------------------------------------------------------------------------------

  IF ( ANY(oem_config%itype_tscale_l(:)==1) ) THEN
    ! Open the NetCDF files
    ierror = nf90_open(hour_of_day_nc, nf90_nowrite, hdid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ierror = nf90_open(day_of_week_nc, nf90_nowrite, dwid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ierror = nf90_open(month_of_year_nc, nf90_nowrite, myid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ! Read hourofday variable
    DO i = 1, tp_ncategory
      ierror = nf90_inq_varid(hdid, oem_config%tp_category(i), varid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_get_var(hdid, varid, oem_data%tp_hourofday(:,i,:), &
         &                  start = (/1,1/), count=(/oem_config%tp_ncountry,tp_param_hourofday/), &
         &                  map = (/ tp_param_hourofday, 1 /))
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
    ENDDO

    ! Read dayofweek variable
    DO i = 1, tp_ncategory
      ierror = nf90_inq_varid(dwid, oem_config%tp_category(i), varid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_get_var(dwid, varid, oem_data%tp_dayofweek(:,i,:), &
        &                   start = (/1,1/), count=(/oem_config%tp_ncountry,tp_param_dayofweek /), &
        &                   map = (/ tp_param_dayofweek, 1 /))
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
    ENDDO

    ! Read monthofyear variable
    DO i = 1, tp_ncategory
      ierror = nf90_inq_varid(myid, oem_config%tp_category(i), varid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_get_var(myid, varid, oem_data%tp_monthofyear(:,i,:), &
         &                  start = (/1,1/), count=(/oem_config%tp_ncountry,tp_param_monthofyear /), &
         &                  map = (/ tp_param_monthofyear, 1 /))
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
    ENDDO

    ! Read countryID variable
    ierror = nf90_inq_varid(myid, "country", varid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ierror = nf90_get_var(myid, varid, oem_data%tp_countryid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ! Close the NetCDF files
    ierror = nf90_close(hdid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ierror = nf90_close(dwid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ierror = nf90_close(myid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF ! hour of day

!------------------------------------------------------------------------------
! Section 1.2: For hour_of_year case
!------------------------------------------------------------------------------

  IF ( ANY(oem_config%itype_tscale_l(:)==2) ) THEN
    ! Open the file
    ierror = nf90_open(hour_of_year_nc, nf90_nowrite, hyid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

    ! Read hourofyear variable
    DO i = 1, tp_ncategory
      ierror = nf90_inq_varid(hyid, oem_config%tp_category(i), varid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_get_var(hyid, varid, oem_data%tp_hourofyear(:,i,:), &
        &                   start = (/1,1/), count=(/oem_config%tp_ncountry,tp_param_hour/), &
        &                   map = (/ tp_param_hour, 1 /))
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
    ENDDO

    IF ( .NOT. ANY(oem_config%itype_tscale_l(:)==1) ) THEN
      ! Read countryID variable if it hasn't been read already
      ierror = nf90_inq_varid(hyid, "country", varid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF

      ierror = nf90_get_var(hyid, varid, oem_data%tp_countryid)
      IF (ierror /= nf90_noerr) THEN
        yerrmsg = TRIM(nf90_strerror(ierror))
        RETURN
      ENDIF
    ENDIF ! .NOT. hourofday

    ! Close the NetCDF file
    ierror = nf90_close(hyid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  ENDIF ! hourofyear


!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_read_temporal_profile_from_file



!==============================================================================
!==============================================================================
!+ Read the vertical profiles from NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE art_oem_read_vertical_profile_from_file(jg, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:  This subroutine reads in the vertical scaling factor, layer_top,
!   and layer_bot variables from NetCDF file and distributes them to all 
!   processors.
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg         !< patch on which computation is performed
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER (LEN= *),INTENT(INOUT)   :: yerrmsg

  !----------------------------------------------------------
  ! Local variables

  INTEGER :: i, ncid, varid

  CHARACTER(*), PARAMETER :: routine = "art_oem_read_vertical_profile_from_file"

!- End of header
!==============================================================================

!------------------------------------------------------------------------------
! Section 1: Open NetCDF file and read in variables
!------------------------------------------------------------------------------

  yerrmsg  = '   '
  ierror   = nf90_noerr

  ! Open the file
  ierror = nf90_open(vertical_profile_nc, nf90_nowrite, ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  DO i = 1, vp_ncategory
    ! Get the id for the vp_factor variable
    ierror = nf90_inq_varid(ncid, oem_config%vp_category(i), varid)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
    ! Read in the vp_factor variable
    ierror = nf90_get_var(ncid, varid, vp_factor(i,:))
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  END DO

  ! Get the id for the layer_bot variable
  ierror = nf90_inq_varid(ncid, "layer_bot", varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Read in the layer_bot variable
  ierror = nf90_get_var(ncid, varid, vp_layer_bot)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Get the id for the layer_top variable
  ierror = nf90_inq_varid(ncid, "layer_top", varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Read in the layer_top variable
  ierror = nf90_get_var(ncid, varid, vp_layer_top)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Close the file
  ierror = nf90_close(ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF


!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_read_vertical_profile_from_file


!==============================================================================
!==============================================================================
!+ Interpolate the vertical profiles to model grid
!------------------------------------------------------------------------------

SUBROUTINE art_oem_interpolate_vertical_profiles(jg, ierror, yerrmsg)

!-----------------------------------------------------------------------------
! Description: This subroutine interpolates the vertical profiles from levels
!              read from the NetCDF file to model vertical levels. 
!
! Method:      This code is directly adapted from the quasi-mass conserving
!              vertical interpolation scheme distvert_3D from the art_interpol.f90
!              file in Int2lm.  
!-----------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)              :: jg         !< patch on which computation is performed
  INTEGER, INTENT(INOUT)           :: ierror
  CHARACTER (LEN= *),INTENT(INOUT) :: yerrmsg


  !--------------------------------------------
  ! Local variables

  INTEGER ::         ncontribLevels,      & ! number of model levels
    &                invertedLevel,       & ! contains the model level 
    &                i, zncat, jb, jc, zke_emiss, zke_lm, zke_contrib, &
    &                nlev, nlev1, is, ie, i_startblk, i_endblk

  REAL(KIND=wp) ::  contributionPart, ub_diff, lb_diff, &
    &               ref_prs_u, ref_prs_l, dpres

  INTEGER, DIMENSION(:), ALLOCATABLE :: &
    &               levelIndizes                     !indices of model levels

  REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: &
    &               levelContribution                ! relative contribution to a certain model level

  REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE ::  &
    &               ub_model,                      & ! upper bounds (geometric height) model levels 
    &               lb_model                         ! lower bounds (geometric height) model levels

  CHARACTER(*), PARAMETER :: routine = "art_oem_interpolate_vertical_profiles"

!- End of header
!==============================================================================

  !----------------------------------------------------------------------------                  
  ! Check for each of the emissions data levels the corresponding model levels.
  ! Each emissions data level represents a volume defined by its upper and lower
  ! bounds lb_emiss and ub_emiss.
  ! Those are distributed to all model levels, which are themselves assumed to
  ! represent volumes ranging from lb_model to ub_model
  ! The distribution is weighted according to the vertical pressure distribution
  ! to reach quasi - mass consistency
  !----------------------------------------------------------------------------

  ierror = 0
  yerrmsg = ''

  i_startblk = art_atmo%i_startblk
  i_endblk = art_atmo%i_endblk

  nlev = art_atmo%nlev
  nlev1 = nlev + 1

  ALLOCATE(levelIndizes(nlev1))
  ALLOCATE(levelContribution(nlev))
  ALLOCATE(ub_model(art_atmo%nproma,art_atmo%nblks,nlev))
  ALLOCATE(lb_model(art_atmo%nproma,art_atmo%nblks,nlev))

  ! Setup arrays containing the upper and lower boundary values (geometric height)
  ! for each model / emissions level
  lb_model(:,:,:)    = 0.0_wp
  lb_model(:,:,1) = art_atmo%z_ifc(:,nlev1,:) - ext_data(jg)%atm%topography_c(:,:)
  ub_model(:,:,:)    = 0.0_wp
  ub_model(:,:,nlev) = art_atmo%z_ifc(:,1,:) - ext_data(jg)%atm%topography_c(:,:)

  ! the model levels are vice versa, we start at the top...
  ! further on, those values include the topography
  ! so we need to subtract this part.
  DO i = 1, nlev-1
    lb_model(:,:,i+1) = art_atmo%z_ifc(:,nlev1-i,:) - ext_data(jg)%atm%topography_c(:,:)
    ub_model(:,:,i  ) = lb_model(:,:,i+1)
  ENDDO

  ! the model levels height depends on topography,
  ! so we need to treat each point separately
  ! btw. there is a lot of optimization potential down here... :)

  DO jb = i_startblk, i_endblk
    CALL art_get_indices_c(jg, jb, is, ie)
    DO jc = is, ie
      DO zke_emiss = 1,vp_nlevel
        ! we loop through each emissions level and look for model levels that fall within
        ! (or model levels where the emissions level falls in)
        ncontribLevels      = 1 !_iintegers
        levelContribution(:)  = 0.0_wp
        levelIndizes(:)     = -1
        ub_diff         = 0.0_wp
        lb_diff = 0.0_wp
        DO zke_lm = 1, nlev
          ! the model levels start from the top
          ! we use ke1, which is ke + 1
          ! because our index (zke_lm) starts a 1
          invertedLevel = nlev1-zke_lm

          ! Possible relationships between model and input level:

          ! chem   model  model
          ! |x+1|  |___|  |___|
          ! |___|  |   |  |   |
          ! |   |  | c |  |   |
          ! |   |  |___|  |   |
          ! | x |  | b |  | d |
          ! |   |  |___|  |   |
          ! |___|  | a |  |   |
          ! |   |  |   |  |___|
          ! |x-1|  |___|  |   |

          ! the emissions data level zke_emiss shares at least partly some height 
          ! with the model level zke_lm
          ! (true for a, b, c and d in the graphic)

          IF (ub_model(jc,jb,zke_lm) >= vp_layer_bot(zke_emiss) .AND. &
            & lb_model(jc,jb,zke_lm) <= vp_layer_top(zke_emiss)) THEN
            ub_diff = ub_model(jc,jb,zke_lm) - vp_layer_top(zke_emiss)
            lb_diff = lb_model(jc,jb,zke_lm) - vp_layer_bot(zke_emiss)
            !WRITE(*,*) "match!!"

            contributionPart = 0.0_wp

            ! now look if we have a full match,
            ! either the total model level fits into the input data level
            ! or vice versa..
            ! (true for b and d in the graphic)

            IF ( (ub_diff >= 0.0_wp) .AND. (lb_diff <= 0.0_wp) .OR. &
              &  (ub_diff <= 0.0_wp) .AND. (lb_diff >= 0.0_wp)) THEN
              ! add the whole model level pressure differences to the contributions...
              contributionPart = 1.0_wp

            ELSE
              ! this level only partly fits.
              ! so the remaining possibilities are
              ! * the model level is at the upper end of the emissions level
              ! (true for c in the graphic)

              IF (ub_diff > 0.0_wp) THEN
                contributionPart = (vp_layer_top(zke_emiss) - lb_model(jc,jb,zke_lm))/ &
                &                  (ub_model(jc,jb,zke_lm) - lb_model(jc,jb,zke_lm))
              ELSE
                ! * or the model level is at the lower end of the input level
                ! (true for a in the graphic)
                contributionPart = (ub_model(jc,jb,zke_lm) - vp_layer_bot(zke_emiss))/ &
                  &                (ub_model(jc,jb,zke_lm) - lb_model(jc,jb,zke_lm))
              ENDIF !ub_diff>0
            ENDIF !(ub_diff >= 0.0_wp) .AND. (lb_diff <= 0.0_wp) .OR. & (ub_diff <= 0.0_wp) .AND. (lb_diff >= 0.0_wp))

            ! the contributions are weighted according to pressure difference within the model level
            ! first, caclulate reference pressere at half-levels mass point (there is no dpres for reference atm)
            ref_prs_u = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg) &
              &  *LOG((EXP(art_atmo%z_ifc(jc,invertedLevel,jb)/h_scal_bg) &
              &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))
            ref_prs_l = p0sl_bg*EXP(-grav/rd*h_scal_bg/(t0sl_bg-del_t_bg)   &
              &  *LOG((EXP(art_atmo%z_ifc(jc,invertedLevel+1,jb)/h_scal_bg) &
              &  *(t0sl_bg-del_t_bg) +del_t_bg)/t0sl_bg))
            dpres = ref_prs_l - ref_prs_u
            levelContribution(ncontribLevels) = dpres * contributionPart
            levelIndizes(ncontribLevels) = invertedLevel
            ncontribLevels = ncontribLevels + 1

          ENDIF !(ub_model(jc,jb,zke_lm) >= vp_layer_bot(zke_emiss) .AND. lb_model(jc,jb,zke_lm) <= vp_layer_top(zke_emiss))
        ENDDO !zke_lm = 1, nlev

        ! now we weigh the contributions to each model level by the pressure differences found
        IF (SUM(levelContribution) > 0._wp) THEN
          ! make _relative_ contributions
          levelContribution = levelContribution / SUM(levelContribution)
        ELSE
          ! this would be really strange, as you would need some emissions levels that do not at
          ! all fall within the model vertical range (i. e. stratospheric :) ).
          ! As we would lose mass we cannot allow this to happen.
          CALL finish(routine,"ERROR: No model level found to put input data in. Check your vertical levels!")
        ENDIF

        ! so now we have two arrays. One filled with the relative contributions
        ! of the model levels in terms of pressure differences;
        ! and one containing the corresponding indizes.
        DO zncat = 1, vp_ncategory
          DO zke_contrib = 1, ncontribLevels
            IF (levelIndizes(zke_contrib) > 0 .AND. levelContribution(zke_contrib) > 0.0_wp) THEN !iintegers wp
              oem_data%vert_scaling_fact(jc,levelIndizes(zke_contrib),jb,zncat) = &
                 & oem_data%vert_scaling_fact(jc,levelIndizes(zke_contrib),jb,zncat) + &
                 & (levelContribution(zke_contrib) * vp_factor(zncat,zke_emiss))
            ENDIF
          ENDDO !zke_contrib = 1, ncontribLevels
        ENDDO !zncat = 1, vp_ncategory
      ENDDO !zke_emiss = 1, vp_nlevel
    ENDDO !jc = is, ie
  ENDDO !jb = i_startblk, i_endblk


!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------ 

END SUBROUTINE art_oem_interpolate_vertical_profiles



!==============================================================================
!==============================================================================
!+ Cleanup the vertical profile data structures 
!------------------------------------------------------------------------------

  SUBROUTINE art_oem_cleanup_vertical_profile_fields()

!------------------------------------------------------------------------------

    IMPLICIT NONE

    IF (ALLOCATED(vp_layer_bot)) DEALLOCATE(vp_layer_bot)
    IF (ALLOCATED(vp_layer_top)) DEALLOCATE(vp_layer_top)
    IF (ALLOCATED(vp_factor)) DEALLOCATE(vp_factor)

  END SUBROUTINE art_oem_cleanup_vertical_profile_fields

!==============================================================================
!==============================================================================
!+ Initialize and read vegetation indices (EVI and LSWI) from NetCDF file
!------------------------------------------------------------------------------

SUBROUTINE art_oem_init_and_read_vegetation_indices(jg, ierror, yerrmsg)

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg           !< patch on which computation is performed
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER(LEN= *), INTENT(INOUT)   :: yerrmsg

  !---------------------------------------------------------------------
  ! Local variables

  LOGICAL :: file_exists

  INTEGER :: i, n, m, ncid, dimid, n_glb_cells, varid, jc, jb,           &
    &        i_startblk, i_endblk, glbidx, is, ie, locidx, dos_start,    &
    &        dos_stop, ndays, trec, nd, ncell, nmodis_days
  REAL, DIMENSION(:,:), ALLOCATABLE :: &
    & lswi_tot, &       ! gridded LSWI on whole domain
    & evi_tot           ! gridded EVI on whole domain

  REAL, DIMENSION(:), ALLOCATABLE :: &
    & lswi_max_tot, &     ! Maximum LSWI over a season on the whole domain
    & lswi_min_tot, &     ! Minimum LSWI over a season on the whole domain
    & evi_max_tot,  &     ! Maximum EVI over a season on the whole domain
    & evi_min_tot         ! Minimum EVI over a season on the whole domain
  CHARACTER(10), DIMENSION(:), ALLOCATABLE :: times ! dtime strings (e.g. 2022-11-01)
    
  CHARACTER(*), PARAMETER :: routine = "art_oem_init_and_read_vegetation_indices"

  TYPE(datetime), POINTER :: modis_date
  

!------------------------------------------------------------------------------
! Section 1: Read in the variables
!------------------------------------------------------------------------------

  yerrmsg = '   '
  ierror = nf90_noerr

  ! Check if the file exists
  INQUIRE(FILE=vegetation_indices_nc, EXIST=file_exists)
  IF (.NOT.(file_exists)) THEN
    ierror=1
    yerrmsg = "Vegetation indices file does not exist! "
    RETURN
  ENDIF

  ! Open the NetCDF file
  ierror = nf90_open(vegetation_indices_nc, nf90_nowrite, ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Get the dimension id for 'time'
  ierror = nf90_inq_dimid(ncid, 'time', dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Get the dimension length for 'time'
  ierror = nf90_inquire_dimension(ncid, dimid, len=ndays)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Get the variable id for 'date_str'
  ierror = nf90_inq_varid(ncid, 'date_str', varid )
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Allocate times array
  ALLOCATE(times(ndays))
  DO i = 1, ndays
    times(i) = ''
  ENDDO

  ! Get the character strings for 'date_str'
    ierror = nf90_get_var(ncid, varid, times)
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF

  ! Map MODIS dates to simulation dates ("dos" = day of simulation)
  dos_start = 0
  dos_stop = 0
  do i = 1, ndays
  ! Get the simulation's start and stop dates
    modis_date => newDatetime(trim(times(i))//trim('T00:00:00'))
    IF ((time_config%tc_startdate%date%year  == modis_date%date%year ) .AND. &
      & (time_config%tc_startdate%date%month == modis_date%date%month) .AND. &
      & (time_config%tc_startdate%date%day   == modis_date%date%day  )) THEN
      dos_start = i
    ENDIF

    IF ((time_config%tc_stopdate%date%year  == modis_date%date%year ) .AND. &
      & (time_config%tc_stopdate%date%month == modis_date%date%month) .AND. &
      & (time_config%tc_stopdate%date%day   == modis_date%date%day  )) THEN
      dos_stop = i
    ENDIF
  enddo

  IF (dos_start < 1) THEN
    ierror = 1
    yerrmsg = 'Variable dos_start must be positive (i.e., MODIS data starts too late?).'
    RETURN
  ENDIF

  IF (dos_stop < 1) THEN
    ierror = 1
    yerrmsg = 'Variable dos_stop must be positive (i.e., MODIS data ends too early?).'
    RETURN
  ENDIF

  ! Get number of global cells
  ierror = nf90_inq_dimid(ncid, 'cell', dimid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ierror = nf90_inquire_dimension(ncid, dimid, len=n_glb_cells)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Deallocate times array
  IF (ALLOCATED(times)) DEALLOCATE(times)

  ! Set number of MODIS days relevant to keep in memory
  nmodis_days = dos_stop - dos_start + 1

  ! Allocate and read in LSWI field from netCDF file
  ALLOCATE(lswi_tot(nmodis_days, n_glb_cells)); lswi_tot = 0.0
  ierror = nf90_inq_varid(ncid, 'lswi', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  DO trec = 1, nmodis_days
    ierror = nf90_get_var(ncid, varid, lswi_tot(trec,:), &
      &                   start = (/1, trec+dos_start/), &
      &                   count=(/n_glb_cells, 1/) & 
      &                   )
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  ENDDO

  ! Allocate and read in LSWI field from netCDF file
  ALLOCATE(evi_tot(nmodis_days, n_glb_cells)); evi_tot = 0.0
  ierror = nf90_inq_varid(ncid, 'evi', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  DO trec = 1, nmodis_days
    ierror = nf90_get_var(ncid, varid, evi_tot(trec,:), &
      &                   start = (/1, trec+dos_start/), &
      &                   count=(/n_glb_cells, 1/) & 
      &                   )
    IF (ierror /= nf90_noerr) THEN
      yerrmsg = TRIM(nf90_strerror(ierror))
      RETURN
    ENDIF
  ENDDO

  ! Allocate and read in LSWI_min field from netCDF file
  ALLOCATE(lswi_min_tot(n_glb_cells)); lswi_min_tot = 0.0
  ierror = nf90_inq_varid(ncid, 'lswi_min', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_get_var(ncid, varid, lswi_min_tot(:))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Allocate and read in EVI_min field from netCDF file
  ALLOCATE(evi_min_tot(n_glb_cells)); evi_min_tot = 0.0
  ierror = nf90_inq_varid(ncid, 'evi_min', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_get_var(ncid, varid, evi_min_tot(:))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Allocate and read in LSWI_max field from netCDF file
  ALLOCATE(lswi_max_tot(n_glb_cells)); lswi_max_tot = 0.0
  ierror = nf90_inq_varid(ncid, 'lswi_max', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_get_var(ncid, varid, lswi_max_tot(:))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Allocate and read in EVI_max field from netCDF file
  ALLOCATE(evi_max_tot(n_glb_cells)); evi_max_tot = 0.0
  ierror = nf90_inq_varid(ncid, 'evi_max', varid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF
  ierror = nf90_get_var(ncid, varid, evi_max_tot(:))
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

  ! Close the NetCDF file
  ierror = nf90_close(ncid)
  IF (ierror /= nf90_noerr) THEN
    yerrmsg = TRIM(nf90_strerror(ierror))
    RETURN
  ENDIF

!------------------------------------------------------------------------------
! Section 2.1: Allocate variables
!------------------------------------------------------------------------------

  ALLOCATE(oem_data%lswi(art_atmo%nproma,art_atmo%nblks,nmodis_days)); oem_data%lswi = 0
  ALLOCATE(oem_data%evi(art_atmo%nproma,art_atmo%nblks,nmodis_days)); oem_data%evi = 0
  ALLOCATE(oem_data%lswi_max(art_atmo%nproma,art_atmo%nblks)); oem_data%lswi_max = 0
  ALLOCATE(oem_data%evi_max(art_atmo%nproma,art_atmo%nblks)); oem_data%evi_max = 0
  ALLOCATE(oem_data%lswi_min(art_atmo%nproma,art_atmo%nblks)); oem_data%lswi_min = 0
  ALLOCATE(oem_data%evi_min(art_atmo%nproma,art_atmo%nblks)); oem_data%evi_min = 0

!------------------------------------------------------------------------------
! Section 2.2: Distribute the values
!------------------------------------------------------------------------------

  i_startblk = art_atmo%i_startblk
  i_endblk = art_atmo%i_endblk

  ! Distribute reflectances values
  DO nd = 1, nmodis_days
    DO jb = i_startblk, i_endblk
      ! read indices within this block:
      CALL art_get_indices_c(jg, jb, is, ie)
      DO jc = is, ie
        ! get local index from jc and jb:
        locidx = idx_1d(jc,jb)
        ! get global index from local index:
        glbidx = p_patch(jg)%cells%decomp_info%glb_index(locidx)
        oem_data%lswi(jc,jb,nd) = lswi_tot(nd,glbidx)
        oem_data%evi(jc,jb,nd) = evi_tot(nd,glbidx)
        oem_data%lswi_max(jc,jb) = lswi_max_tot(glbidx)
        oem_data%evi_max(jc,jb) = evi_max_tot(glbidx)
        oem_data%lswi_min(jc,jb) = lswi_min_tot(glbidx)
        oem_data%evi_min(jc,jb) = evi_min_tot(glbidx)
      ENDDO
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_init_and_read_vegetation_indices



!==============================================================================
!==============================================================================
!+ Initialize VPRM land-use classes 
!------------------------------------------------------------------------------

SUBROUTINE art_oem_init_vprm_lu_class_fraction(jg, ext_data, ierror, yerrmsg)

  IMPLICIT NONE

  INTEGER, INTENT(IN)                :: jg           !< patch on which computation is performed
  TYPE(t_external_data), INTENT(IN)  :: ext_data     !< external data
  INTEGER, INTENT(INOUT)             :: ierror
  CHARACTER(LEN= *), INTENT(INOUT)   :: yerrmsg

  !---------------------------------------------------------------------
  ! Local variables

  INTEGER :: i, n, m, jc, jb, i_startblk, i_endblk, is, ie, ncell

  INTEGER, PARAMETER :: num_vprm_lc_classes = 8

  CHARACTER(*), PARAMETER :: routine = "art_oem_init_vprm_lu_class_fraction"

!------------------------------------------------------------------------------
! Section 1: Set the indices for the VPRM land-use classes
!------------------------------------------------------------------------------

  oem_data%i_vprm_lc_evergreen = 1
  oem_data%i_vprm_lc_deciduous = 2
  oem_data%i_vprm_lc_mixed = 3
  oem_data%i_vprm_lc_shrub = 4
  oem_data%i_vprm_lc_savanna = 5
  oem_data%i_vprm_lc_crop = 6
  oem_data%i_vprm_lc_grass = 7
  oem_data%i_vprm_lc_urban = 8

!------------------------------------------------------------------------------
! Section 2.1: Allocate variables
!------------------------------------------------------------------------------

  ALLOCATE(oem_data%vprm_lu_class_fraction(art_atmo%nproma,art_atmo%nblks,num_vprm_lc_classes))
  oem_data%vprm_lu_class_fraction = 0

!------------------------------------------------------------------------------
! Section 2.2: Map the ICON to VPRM land-use classes
!------------------------------------------------------------------------------

  i_startblk = art_atmo%i_startblk
  i_endblk = art_atmo%i_endblk

  DO jb = i_startblk, i_endblk
    ! read indices within this block:
    CALL art_get_indices_c(jg, jb, is, ie)
    DO jc = is, ie
      ! Index 1: Evergreen
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_evergreen)    = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_b_eg)  + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_n_eg)
      ! Index 2: Deciduous
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_deciduous)    = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_b_d)   + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_woodland)     + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_n_d)
      ! Index 3: Mixed Forest
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_mixed)        = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_bn)    + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_rf)    + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_forest_pf)
      ! Index 4: Shrubland
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_shrub)        = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_shrub_mos)    + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_shrub_eg)     + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_shrub)        + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_sparse) * 0.5
      ! Index 5: Savanna
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_savanna)      = 0
      ! Index 6: Cropland
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_crop)         = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_crop_irrig)   + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_crop_rain)    + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_crop_mos)
      ! Index 7: Grassland
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_grass)        = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_veg_mos)      + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_grass)        + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_grass_rf)
      ! Index 8: Urban Area
      oem_data%vprm_lu_class_fraction(jc,jb,oem_data%i_vprm_lc_urban)        = &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_sparse) * 0.5 + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_urban)        + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_bare_soil)    + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_water)        + &
        & ext_data%atm%lu_class_fraction(jc,jb,ext_data%atm%i_lc_snow_ice)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE art_oem_init_vprm_lu_class_fraction
!==============================================================================
!==============================================================================

END MODULE mo_art_oem_init

