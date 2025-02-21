!NEC$ options "-O1"
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

MODULE mo_nwp_phy_state

!
! Description:  Contains the data structures
!  to store the physical model state and other auxiliary variables
!  in order to run the ECHAM physics.
!  Constructors and destructors for these data structures as well as
!  initialization of fields are also defined here.
!  This module should be an analogon to 'mo_hydro_state.f90'
!
!  TODO/To think about:
!     - should physics be called before or after dynamics?
!     - allocate fluxes at edges instead at the centers?
!     - horizontal/vertical tracer flux (reconstruct q'v_n' into q'u' and q'v') ?
!     - provide the "virt_inc" with meaning
!     - where to provide the lat/lon info for radiation?
!     - how to implement the echam-modules - rewriting them or "capsulate"?
!     - revision of fields if there are needed or tp be replaced
!     - fill the physics tendency construction/destruction subroutine
!     - later implement already calculated icon gradients for echam physics
!     - think about variables for flexible time steps

USE mo_kind,                ONLY: wp, i8
USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag, t_nwp_phy_tend, t_nwp_phy_stochconv
USE mo_impl_constants,      ONLY: success, &
  &                               VINTP_METHOD_LIN,VINTP_METHOD_QV,   &
  &                               TASK_COMPUTE_RH, TASK_COMPUTE_PV,   &
  &                               TASK_COMPUTE_SDI2,                  &
  &                               TASK_COMPUTE_LPI,                   &
  &                               TASK_COMPUTE_CEILING,               &
  &                               TASK_COMPUTE_VIS,                   &
  &                               TASK_COMPUTE_HBAS_SC,               &
  &                               TASK_COMPUTE_HTOP_SC,               &
  &                               TASK_COMPUTE_TWATER,                &
  &                               TASK_COMPUTE_Q_SEDIM,               &
  &                               TASK_COMPUTE_DBZ850,                &
  &                               TASK_COMPUTE_DBZLMX_LOW,            &
  &                               TASK_COMPUTE_DBZCMAX,               &
  &                               TASK_COMPUTE_WSHEAR_U,              &
  &                               TASK_COMPUTE_WSHEAR_V,              &
  &                               TASK_COMPUTE_LAPSERATE,             &
  &                               TASK_COMPUTE_SRH,                   &
  &                               TASK_COMPUTE_INVERSION,             &
  &                               ivdiff,                             &
  &                               LSS_JSBACH,                         &
  &                               HINTP_TYPE_LONLAT_NNB,              &
  &                               HINTP_TYPE_LONLAT_BCTR,             &
  &                               HINTP_TYPE_LONLAT_RBF,              &
  &                               RTTOV_BT_CL,                        &
  &                               RTTOV_RAD_CL, RTTOV_RAD_CS,         &
  &                               iss, iorg, ibc, iso4,               &
  &                               idu, nclass_aero, vname_len
USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL,             &
  &                               GRID_CELL
USE mo_parallel_config,     ONLY: nproma
USE mo_run_config,          ONLY: nqtendphy, iqv, iqc, iqi, iqr, iqs, iqg, iqh, lart, ldass_lhn
USE mo_exception,           ONLY: message, finish !,message_text
USE mo_model_domain,        ONLY: t_patch, p_patch, p_patch_local_parent
USE mo_grid_config,         ONLY: n_dom, n_dom_start, nexlevs_rrg_vnest
USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config, icpl_aero_conv, iprog_aero
USE turb_data,              ONLY: ltkecon
USE mo_initicon_config,     ONLY: icpl_da_sfcevap, icpl_da_snowalb, icpl_da_skinc, icpl_da_seaice
USE mo_radiation_config,    ONLY: irad_aero, iRadAeroTegen, iRadAeroART
USE mo_lnd_nwp_config,      ONLY: ntiles_total, ntiles_water, nlev_soil
USE mo_nwp_vdiff_interface, ONLY: nwp_vdiff_setup
USE mo_var_list,            ONLY: add_var, add_ref, t_var_list_ptr
USE mo_var_list_register,   ONLY: vlr_add, vlr_del
USE mo_var_groups,          ONLY: groups, MAX_GROUPS
USE mo_var_metadata_types,  ONLY: POST_OP_SCALE, POST_OP_LIN2DBZ, CLASS_SYNSAT, CLASS_CHEM, &
  &                               CLASS_TILE, CLASS_TILE_LAND
USE mo_var_metadata,        ONLY: create_vert_interp_metadata,  &
  &                               create_hor_interp_metadata,   &
  &                               vintp_types, post_op
USE mo_nwp_parameters,      ONLY: t_phy_params
USE mo_cf_convention,       ONLY: t_cf_var
USE mo_grib2,               ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+)
USE mo_cdi,                 ONLY: TSTEP_MIN, TSTEP_MAX, TSTEP_INSTANT, TSTEP_CONSTANT, &
    &                             TSTEP_AVG, TSTEP_ACCUM, DATATYPE_PACK16,             &
    &                             DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED, GRID_LONLAT, DATATYPE_INT32
USE mo_zaxis_type,          ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF,          &
  &                               ZA_SURFACE, ZA_HEIGHT_2M, ZA_HEIGHT_10M,       &
  &                               ZA_HEIGHT_2M_LAYER, ZA_TOA, ZA_DEPTH_BELOW_LAND,   &
  &                               ZA_PRESSURE_0, ZA_PRESSURE_400, ZA_SRH, &
  &                               ZA_PRESSURE_800, ZA_CLOUD_BASE, ZA_CLOUD_TOP,  &
  &                               ZA_ISOTHERM_ZERO, ZA_ECHOTOP, ZA_WSHEAR, ZA_PRESSURE_LAPSERATE
USE mo_physical_constants,  ONLY: grav

USE mo_ls_forcing_nml,      ONLY: is_ls_forcing
USE mo_les_nml,             ONLY: turb_profile_list, turb_tseries_list

USE mo_advection_config,     ONLY: advection_config
USE mo_synsat_config,        ONLY: lsynsat, num_images, get_synsat_name, num_sensors, &
  &                                total_numchans, get_synsat_grib_triple
USE mo_art_config,           ONLY: nart_tendphy



USE mo_comin_config,         ONLY: comin_config
USE mo_action,               ONLY: ACTION_RESET, new_action, actions
USE mo_io_config,            ONLY: lflux_avg, lnetcdf_flt64_output, gust_interval, &
  &                                celltracks_interval, echotop_meta, &
  &                                maxt_interval, precip_interval, t_var_in_output, &
  &                                totprec_d_interval, itype_hzerocl, &
  &                                uh_max_zmin, uh_max_zmax, luh_max_out, uh_max_nlayer, &
  &                                sunshine_interval, n_wshear, n_srh
USE mtime,                   ONLY: max_timedelta_str_len, getPTStringFromMS
USE mo_name_list_output_config, ONLY: is_variable_in_output
USE mo_util_string,          ONLY: real2string
USE mo_sbm_storage,          ONLY: construct_sbm_storage, destruct_sbm_storage
USE mo_coupling_config,      ONLY: is_coupled_to_waves
USE mo_netcdf,               ONLY: NF_MAX_NAME

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




IMPLICIT NONE
PRIVATE

!public interface
!
! subroutines
PUBLIC :: construct_nwp_phy_state
PUBLIC :: destruct_nwp_phy_state
!variables
PUBLIC :: prm_diag 
PUBLIC :: prm_nwp_tend
PUBLIC :: prm_nwp_stochconv
PUBLIC :: phy_params
PUBLIC :: prm_nwp_diag_list  !< variable lists
PUBLIC :: prm_nwp_tend_list  !< variable lists
PUBLIC :: prm_nwp_stochconv_list  !< variable lists

!

!!--------------------------------------------------------------------------
!!                          STATE VARIABLES 
!!--------------------------------------------------------------------------
  TYPE(t_nwp_phy_diag), ALLOCATABLE, TARGET :: prm_diag(:) !< shape: (n_dom)
  TYPE(t_nwp_phy_tend), ALLOCATABLE :: prm_nwp_tend(:)     !< shape: (n_dom)
  TYPE(t_nwp_phy_stochconv), ALLOCATABLE, TARGET :: prm_nwp_stochconv(:) !< shape: (n_dom)  
!-------------------------------------------------------------------------
  
!!--------------------------------------------------------------------------
!!                          VARIABLE LISTS
!!--------------------------------------------------------------------------
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_nwp_diag_list(:)  !< shape: (n_dom)
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_nwp_tend_list(:)  !< shape: (n_dom)
  TYPE(t_var_list_ptr),ALLOCATABLE :: prm_nwp_stochconv_list(:)  !< shape: (n_dom)

!!-------------------------------------------------------------------------
!! Parameters of various physics parameterizations that have to be 
!! domain-dependent (computed during physics initialization phase)
!!-------------------------------------------------------------------------
  TYPE (t_phy_params), ALLOCATABLE :: phy_params(:)  !< shape: (n_dom)


CONTAINS

!-------------------------------------------------------------------------

SUBROUTINE construct_nwp_phy_state( p_patch, var_in_output )
  TYPE(t_patch), TARGET, INTENT(in) :: p_patch(n_dom)
  TYPE(t_var_in_output), INTENT(in) :: var_in_output(n_dom)
  CHARACTER(LEN=21)                 :: listname
  INTEGER :: jg, ist, nblks_c, nlev, nlevp1
  CHARACTER(*), PARAMETER :: routine = 'mo_nwp_phy_state:construct_nwp_state'

!-------------------------------------------------------------------------

  CALL message(routine, 'start to construct 3D state vector')

  ! Allocate pointer arrays prm_diag_nwp and prm_nwp_tend, 
  ! as well as the corresponding list arrays.

  ALLOCATE(prm_diag(n_dom), prm_nwp_diag_list(n_dom),STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of diagnostic physical array and list failed')
  ENDIF

  ALLOCATE(prm_nwp_tend(n_dom), prm_nwp_tend_list(n_dom), STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of tendency physical array and list failed')
  ENDIF

  ALLOCATE(prm_nwp_stochconv(n_dom), prm_nwp_stochconv_list(n_dom), STAT=ist)
  IF(ist/=success)THEN
    CALL finish ('mo_nwp_phy_state:construct_nwp_state', &
      'allocation of stochastic convection physical array and list failed')
  ENDIF 

  !$ACC ENTER DATA COPYIN(prm_diag, prm_nwp_tend)

  DO jg = 1, n_dom

     !determine size of arrays
     nblks_c = p_patch(jg)%nblks_c
     
     ! number of vertical levels
     nlev   = p_patch(jg)%nlev
     nlevp1 = p_patch(jg)%nlevp1

     WRITE(listname,'(a,i2.2)') 'prm_diag_of_domain_',jg
     CALL new_nwp_phy_diag_list( jg, nlev, nlevp1, nblks_c, listname,  &
       &                         prm_nwp_diag_list(jg), prm_diag(jg), var_in_output(jg) )
     WRITE(listname,'(a,i2.2)') 'prm_tend_of_domain_',jg
     CALL new_nwp_phy_tend_list ( jg, nlev, nblks_c, listname, &
       &                          prm_nwp_tend_list(jg), prm_nwp_tend(jg))
     WRITE(listname,'(a,i2.2)') 'prm_stch_of_domain_',jg
     CALL new_nwp_phy_stochconv_list ( jg, nblks_c, listname, &
       &                          prm_nwp_stochconv_list(jg), prm_nwp_stochconv(jg))

     ! create additional SBM-specific storage, if SBM microphysics is selected
     IF (atm_phy_nwp_config(jg)%inwp_gscp == 8) THEN
       CALL construct_sbm_storage (p_patch)
     ENDIF
  ENDDO

  ! Allocate variable of type t_phy_params containing domain-dependent parameters
  ALLOCATE(phy_params(n_dom), STAT=ist)
  !$ACC ENTER DATA CREATE(phy_params)
  IF(ist/=success) CALL finish(routine, 'allocation of phy_params array failed')


  CALL message(routine, 'construction of state vector finished')

END SUBROUTINE construct_nwp_phy_state

!
SUBROUTINE destruct_nwp_phy_state()

  INTEGER :: jg, ist  !< grid level/domain index

  CHARACTER(len=*), PARAMETER :: &
    routine = 'mo_nwp_phy_state:destruct_nwp_phy_state'
  CALL message(routine, 'start to destruct 3D state vector')

  !$ACC WAIT(1)
  DO jg = 1,n_dom
    CALL vlr_del(prm_nwp_diag_list(jg))
    CALL vlr_del(prm_nwp_tend_list(jg))
    CALL vlr_del(prm_nwp_stochconv_list(jg))
        

    IF (ALLOCATED(prm_diag(jg)%synsat_image))   DEALLOCATE(prm_diag(jg)%synsat_image)
    IF (ASSOCIATED(prm_diag(jg)%qrs_flux) .AND. .NOT. ldass_lhn) THEN
      !$ACC EXIT DATA DELETE(prm_diag(jg)%qrs_flux)
      DEALLOCATE(prm_diag(jg)%qrs_flux)
    ENDIF

  END DO

  !$ACC EXIT DATA DELETE(prm_diag, prm_nwp_tend)

  DEALLOCATE(prm_diag, prm_nwp_diag_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish (routine, &
       &  'deallocation of NWP physics diagnostic array and list failed')
  ENDIF
 
  DEALLOCATE(prm_nwp_tend, prm_nwp_tend_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish (routine, &
       &  'deallocation of NWP physics tendencies array and list failed')
  ENDIF

  DEALLOCATE(prm_nwp_stochconv, prm_nwp_stochconv_list, STAT=ist)
  IF(ist/=success)THEN
    CALL finish (routine, &
       &  'deallocation of NWP stochastic convection array and list failed')
  ENDIF 

  !$ACC EXIT DATA DELETE(phy_params)
  DEALLOCATE(phy_params, STAT=ist)
  IF(ist/=success)THEN
    CALL finish (routine, 'deallocation of phy_params array failed')
  ENDIF

  ! destruct SBM-specific storage
  IF (atm_phy_nwp_config(jg)%inwp_gscp == 8) THEN
    CALL destruct_sbm_storage ()
  ENDIF

  CALL message(routine, 'destruction of 3D state vector finished')

END SUBROUTINE destruct_nwp_phy_state

     !
SUBROUTINE new_nwp_phy_diag_list( k_jg, klev, klevp1, kblks,    &
                     &  listname, diag_list, diag, var_in_output)

    INTEGER,INTENT(IN) :: klev, klevp1, kblks, k_jg !< dimension sizes

    CHARACTER(len=*),INTENT(IN)     :: listname
    CHARACTER(len=*), PARAMETER :: avg_pfx = 'avg_', tot_pfx = 'tot_'
    CHARACTER(LEN=1)                :: csfc
    CHARACTER(LEN=2)                :: caer

    TYPE(t_var_list_ptr)    ,INTENT(INOUT) :: diag_list
    TYPE(t_nwp_phy_diag),INTENT(INOUT) :: diag
    TYPE(t_var_in_output), INTENT(IN)  :: var_in_output

    ! Local variables
    INTEGER :: n_updown = 8 !> number of up/downdrafts variables

    TYPE(t_cf_var)    ::    cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d(3), shape3dsubs(3), &
      &        shape3dsubsw(3), shape3d_synsat(3),     &
      &        shape2d_synsat(2), shape3d_aero(3), shape3dechotop(3), shape3dwshear(3),shape3d_hail(3)
    INTEGER :: shape3dkp1(3), shape3dflux(3), shape3d_uh_max(3), shape3dturb(3), shape3dsrh(3)
    INTEGER :: ibits,  kcloud
    INTEGER :: jsfc, ist
    CHARACTER(len=NF_MAX_NAME) :: long_name
    CHARACTER(len=21) :: name
    CHARACTER(len=3)  :: prefix
    CHARACTER(len=8)  :: meaning
    CHARACTER(len=10) :: varunits  ! variable units, depending on "lflux_avg"
    INTEGER :: a_steptype
    LOGICAL :: lrestart, lrestart_flux

    LOGICAL :: lradiance, lcloudy
    INTEGER :: ichan, idiscipline, icategory, inumber, &
      &        wave_no, wave_no_scalfac, iimage, isens, k
    CHARACTER(LEN=vname_len) :: shortname
    CHARACTER(LEN=128)         :: longname, unit
    CHARACTER(len=max_timedelta_str_len) :: gust_int, celltracks_int,   &
      &                                     echotop_int
    ! For lpi_con_max need an hourly reset for the first 48 h,
    ! a 3-hourly reset for day 3 and 4, and a 6 hourly reset thereafter.
    ! lpi_stop 3 is not needed - it is the end of the simulation
    ! We use the same variables for mlpi_con_max and lfd_con_max.
    CHARACTER(len=max_timedelta_str_len) :: lpi_int1, lpi_int2, lpi_int3
    CHARACTER(len=max_timedelta_str_len) :: lpi_start1, lpi_start2, lpi_start3
    CHARACTER(len=max_timedelta_str_len) :: lpi_end1, lpi_end2
    !
    INTEGER :: constituentType                 ! for variable of class 'chem'

    INTEGER :: datatype_flt
    LOGICAL :: in_group(MAX_GROUPS)            ! for adding a variable to one or more groups

    CHARACTER(len=*), PARAMETER :: &
      routine = 'mo_nwp_phy_state:new_nwp_phy_diag_list'

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! bits "entropy" of horizontal slice

    shape2d        = (/nproma,               kblks            /)
    shape3d        = (/nproma, klev,         kblks            /)
    shape3dkp1     = (/nproma, klevp1,       kblks            /)
    shape3dsubs    = (/nproma, kblks,        ntiles_total     /)
    shape3dsubsw   = (/nproma, kblks,        ntiles_total+ntiles_water /)
    shape3d_aero   = (/nproma, nclass_aero,  kblks            /)
    shape3dechotop = (/nproma, echotop_meta(k_jg)%nechotop, kblks/)
    shape3d_uh_max = (/nproma, kblks,        uh_max_nlayer    /)
    shape3dwshear  = (/nproma, n_wshear,     kblks/)
    shape3dsrh     = (/nproma, n_srh,        kblks/)
    shape3d_hail   = (/nproma, 5,            kblks/)

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY( &
      &     diag%aercl_bc, &
      &     diag%aercl_du, &
      &     diag%aercl_or, &
      &     diag%aercl_ss, &
      &     diag%aercl_su, &
      &     diag%aerosol, &
      &     diag%albdif_t, &
      &     diag%albnirdif_t, &
      &     diag%albvisdif_t, &
      &     diag%bruvais, &
      &     diag%buffer_rrg, &
      &     diag%buffer_rttov, &
      &     diag%cloudtop, &
      &     diag%si, &
      &     diag%sli, &
      &     diag%swiss12, &
      &     diag%swiss00, &
      &     diag%cape_3km, &
      &     diag%lfc_ml, &
      &     diag%lcl_ml, &
      &     diag%cape_mu, &
      &     diag%clc_rad, &
      &     diag%ceiling_height, &
      &     diag%cin_mu, &
      &     diag%dbz3d_lin, &
      &     diag%dbz_850, &
      &     diag%dbzlmx_low, &
      &     diag%dbz_cmax, &
      &     diag%dbz_ctmax, &
      &     diag%dhail, &
      &     diag%dhail_av, &
      &     diag%dhail_mx, &
      &     diag%dhail_sd, &
      &     diag%dursun, &
      &     diag%dursun_m, &
      &     diag%dursun_r, &
      &     diag%echotop, &
      &     diag%echotopinm, &
      &     diag%graupel_gsp, &
      &     diag%graupel_gsp_rate, &
      &     diag%hail_gsp, &
      &     diag%hail_gsp_rate, &
      &     diag%hbas_sc, &
      &     diag%heatcap_fac, &
      &     diag%heatcond_fac, &
      &     diag%htop_sc, &
      &     diag%ice_gsp, &
      &     diag%ice_gsp_rate, &
      &     diag%inversion_height, &
      &     diag%lapse_rate, &
      &     diag%lhn_diag, &
      &     diag%liqfl_turb, &
      &     diag%low_ent_zone, &
      &     diag%lpi, &
      &     diag%lwflxsfc_t, &
      &     diag%mech_prod, &
      &     diag%p_cbase, &
      &     diag%p_ctop, &
      &     diag%pv, &
      &     diag%q_sedim, &
      &     diag%qrs_flux, &
      &     diag%qvtend_lhn, &
      &     diag%reff_qc, &
      &     diag%reff_qg, &
      &     diag%reff_qh, &
      &     diag%reff_qi, &
      &     diag%reff_qr, &
      &     diag%reff_qs, &
      &     diag%rh, &
      &     diag%sdi2, &
      &     diag%snowalb_fac, &
      &     diag%srh, &
      &     diag%swflxsfc_t, &
      &     diag%synsat_arr, &
      &     diag%t_cbase, &
      &     diag%t_ctop, &
      &     diag%tcond10_max, &
      &     diag%tcond_max, &
      &     diag%tetfl_turb, &
      &     diag%tkred_sfc, &
      &     diag%tkred_sfc_h, &
      &     diag%tt_lheat, &
      &     diag%ttend_lhn, &
      &     diag%twater, &
      &     diag%uh_max_3d, &
      &     diag%vapfl_turb, &
      &     diag%vorw_ctmax, &
      &     diag%wdur, &
      &     diag%wup_mask, &
      &     diag%w_ctmax, &
      &     diag%wshear_u, &
      &     diag%wshear_v, &
      &     diag%z_pbl)


    ! Register a field list and apply default settings

    CALL vlr_add(diag_list, TRIM(listname), patch_id=k_jg, lrestart=.TRUE.)
   
    !------------------------------
    ! Meteorological quantities
    !------------------------------

    !-------------------
    ! Clouds and precip
    !------------------
    ! 2D and 3D variables

    kcloud= 3


    ! &      diag%tot_prec_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_rate', 'kg m-2 s-1', 'total precipitation rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec_rate', diag%tot_prec_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    

    ! &      diag%prec_gsp_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_gsp_rate', 'kg m-2 s-1',                      &
      &          'gridscale precipitation rate', datatype_flt)
    grib2_desc = grib2_var(0, 1, 54, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_gsp_rate', diag%prec_gsp_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    

    ! &      diag%rain_gsp_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp_rate', 'kg m-2 s-1', 'gridscale rain rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp_rate', diag%rain_gsp_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT,                                   &
                & in_group=groups("iau_init_vars"),                          &
                & lopenacc=.TRUE. )
    


    ! &      diag%snow_gsp_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp_rate', 'kg m-2 s-1', 'gridscale snow rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp_rate', diag%snow_gsp_rate,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & isteptype=TSTEP_INSTANT,                                   &
                & in_group=groups("iau_init_vars"),                          &
                & lopenacc=.TRUE. )
    

    ! For graupel scheme 
    IF (atm_phy_nwp_config(k_jg)%lhave_graupel) THEN
      
      ! &      diag%graupel_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('graupel_gsp_rate', 'kg m-2 s-1', 'gridscale graupel rate', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 75, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'graupel_gsp_rate', diag%graupel_gsp_rate,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d,                                             &
                  & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    
    ENDIF

    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (1,2,3,4,5,6,7,8)

       ! &      diag%ice_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('ice_gsp_rate', 'kg m-2 s-1', 'gridscale ice rate', &
        &                   datatype_flt)
      ! Note (GZ): the current grib encoding pertains to ice pellets and is thus inappropriate
      grib2_desc = grib2_var(0, 1, 68, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'ice_gsp_rate', diag%ice_gsp_rate,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, isteptype=TSTEP_INSTANT,                    &
                  & in_group=groups("iau_init_vars"),                          &
                  & lopenacc=.TRUE.  )
      
    END SELECT

    !For two moment microphysics
    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (4,5,6,7,8)

       ! &      diag%hail_gsp_rate(nproma,nblks_c)
      cf_desc    = t_cf_var('hail_gsp_rate', 'kg m-2 s-1', 'gridscale hail rate', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 73, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hail_gsp_rate', diag%hail_gsp_rate,            &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
      
    END SELECT


    ! &      diag%rain_con_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con_rate', 'kg m-2 s-1', 'convective rain rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con_rate', diag%rain_con_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE.)
    


    ! &      diag%snow_con_rate(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con_rate', 'kg m-2 s-1', 'convective snow rate', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con_rate', diag%snow_con_rate,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE.)
    

    ! &      diag%rain_con_rate_3d(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('rain_con_rate_3d', 'kg m-2 s-1',                &
      &          '3d convective rain rate', datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con_rate_3d', diag%rain_con_rate_3d,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & lrestart = lart,                                            & ! .TRUE. for ART 
                & ldims=shape3dkp1,                                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE.)
    


    ! &      diag%snow_con_rate_3d(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('snow_con_rate_3d', 'kg m-2 s-1',                   &
      &          '3d convective snow rate', datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con_rate_3d', diag%snow_con_rate_3d,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                & ldims=shape3dkp1,                                           &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE.)
    


    ! &      diag%rain_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp ', 'kg m-2 ', 'gridscale rain ', datatype_flt)
    grib2_desc = grib2_var(0, 1, 77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp', diag%rain_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d, in_group=groups("precip_vars"),             &
                & isteptype=TSTEP_ACCUM ,                                    &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                   &
                & initval=0._wp, resetval=0._wp,                             &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE. )
    

    ! &      diag%rain_gsp0(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_gsp0', 'kg m-2 ', 'gridscale rain0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 77, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_gsp0', diag%rain_gsp0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM, lopenacc=.TRUE. )
    


    ! &      diag%snow_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp', 'kg m-2 ', 'gridscale snow', datatype_flt)
    grib2_desc = grib2_var(0, 1, 56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp', diag%snow_gsp,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & in_group=groups("precip_vars"),                            &
                & isteptype=TSTEP_ACCUM ,                                    &
                & hor_interp=create_hor_interp_metadata(                     &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                   &
                & initval=0._wp, resetval=0._wp,                             &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE. )
    

    ! &      diag%snow_gsp0(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_gsp0', 'kg m-2 ', 'gridscale snow0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 56, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_gsp0', diag%snow_gsp0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM, lopenacc=.TRUE. )
    


    !Surface precipitation variables for graupel scheme and two moment microphysics
    IF (atm_phy_nwp_config(k_jg)%lhave_graupel) THEN
       ! &      diag%graupel_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('graupel_gsp', 'kg m-2', 'gridscale graupel',      &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 75, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'graupel_gsp', diag%graupel_gsp,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM,                                     &
                  & hor_interp=create_hor_interp_metadata(                     &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                   &
                  & initval=0._wp, resetval=0._wp,                             &
                  & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                  & lopenacc=.TRUE.)
    
    ENDIF

    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (1,2,3,4,5,6,7,8)

       ! &      diag%ice_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('ice_gsp', 'kg m-2', 'gridscale ice', datatype_flt)
      ! Note (GZ): the current grib encoding pertains to ice pellets and is thus inappropriate
      grib2_desc = grib2_var(0, 1, 68, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'ice_gsp', diag%ice_gsp,                        &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d,                                             &
             !!!  & in_group=groups("precip_vars"), !! GZ: removed from group until correct grib encoding is available
                  & isteptype=TSTEP_ACCUM,                                     &
                  & hor_interp=create_hor_interp_metadata(                     &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                   &
                  & initval=0._wp, resetval=0._wp,                             &
                  & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                  & lopenacc=.TRUE. )
    
    END SELECT

    SELECT CASE (atm_phy_nwp_config(k_jg)%inwp_gscp)
    CASE (4,5,6,7,8)

       ! &      diag%hail_gsp(nproma,nblks_c)
      cf_desc    = t_cf_var('hail_gsp', 'kg m-2', 'gridscale hail', datatype_flt)
      grib2_desc = grib2_var(0, 1, 73, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hail_gsp', diag%hail_gsp,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                  & ldims=shape2d, in_group=groups("precip_vars"),             &
                  & isteptype=TSTEP_ACCUM,                                     &
                  & hor_interp=create_hor_interp_metadata(                     &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                   &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                   &
                  & initval=0._wp, resetval=0._wp,                             &
                  & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                  & lopenacc=.TRUE.)
    

    END SELECT

    ! &      diag%prec_gsp(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_gsp', 'kg m-2', 'gridscale precip', datatype_flt)
    grib2_desc = grib2_var(0, 1, 54, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_gsp', diag%prec_gsp,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE.)
    



    ! &      diag%rain_con(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con', 'kg m-2 ', 'convective rain', datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con', diag%rain_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, in_group=groups("precip_vars"),              &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE. )
    


    ! &      diag%rain_con0(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_con0', 'kg m-2 ', 'convective rain0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 76, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_con0', diag%rain_con0,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,   &
                & ldims=shape2d,                                             &
                & loutput=.false., lrestart=.TRUE.,                          &
                & isteptype=TSTEP_ACCUM, &
                & lopenacc=.TRUE. )
    


    ! &      diag%snow_con(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con', 'kg m-2', 'convective snow', datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con', diag%snow_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE. )
    

    ! &      diag%snow_con0(nproma,nblks_c)
    cf_desc    = t_cf_var('snow_con0', 'kg m-2', 'convective snow0', datatype_flt)
    grib2_desc = grib2_var(0, 1, 55, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'snow_con0', diag%snow_con0,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & loutput=.false., lrestart=.TRUE.,                           &
                & isteptype=TSTEP_ACCUM, &
                & lopenacc=.TRUE. )
    



    ! &      diag%prec_con(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_con', 'kg m-2', 'convective precip', datatype_flt)
    grib2_desc = grib2_var(0, 1, 37, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_con', diag%prec_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE. )
    


    ! &      diag%tot_prec(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec', 'kg m-2', 'total precip', datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec', diag%tot_prec,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
                & in_group=groups("precip_vars"),                             &
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE. )
    

    ! &      diag%tot_prec_d(nproma,nblks_c) = diag%prec_gsp_d(nproma,nblks_c) + diag%prec_con_d(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_d', 'kg m-2', 'total precip since end of previous full '// &
                 TRIM(totprec_d_interval(k_jg)(3:))//' since model start', datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
         &           + t_grib2_int_key("typeOfStatisticalProcessing", 4)
    CALL add_var( diag_list, 'tot_prec_d', diag%tot_prec_d,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
!                & in_group=groups("precip_vars"),                             &  ! is an optional diagnostic
                & isteptype=TSTEP_ACCUM,                                      &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,totprec_d_interval(k_jg))), &
                & lopenacc=.TRUE. )
    

    ! &      diag%prec_gsp_d(nproma,nblks_c) !! have to be in restart
    cf_desc    = t_cf_var('prec_gsp_d', 'kg m-2', 'gridscale precip since end of previous full '// &
                          TRIM(totprec_d_interval(k_jg)(3:))//' since model start', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &  ! no grib definitions available yet
         &           + t_grib2_int_key("typeOfStatisticalProcessing", 4)
    CALL add_var( diag_list, 'prec_gsp_d', diag%prec_gsp_d,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
!                & in_group=groups("precip_vars"),                             &  ! is an optional diagnostic
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,totprec_d_interval(k_jg))), &
                & lopenacc=.TRUE.)
    

    ! &      diag%prec_con_d(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_con_d', 'kg m-2', 'convective precip since end of previous full '// &
                          TRIM(totprec_d_interval(k_jg)(3:))//' since model start', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &  ! no grib definitions available yet
         &           + t_grib2_int_key("typeOfStatisticalProcessing", 4)
    CALL add_var( diag_list, 'prec_con_d', diag%prec_con_d,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,                                              &
!                & in_group=groups("precip_vars"),                             &  ! is an optional diagnostic
                & isteptype=TSTEP_ACCUM ,                                     &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,totprec_d_interval(k_jg))), &
                & lopenacc=.TRUE. )
    


    ! &      diag%prec_con_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_con_rate_avg', 'kg m-2 s-1',                  &
      &          'convective precip rate, time average', datatype_flt)
    grib2_desc = grib2_var(0, 1, 37, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_con_rate_avg', diag%prec_con_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,  lrestart=.FALSE.,                           &
                & in_group=groups("additional_precip_vars"),                  &
                & isteptype=TSTEP_AVG,                                        &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE.                                             )
    

    ! &      diag%prec_gsp_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('prec_gsp_rate_avg', 'kg m-2 s-1',                  &
      &          'gridscale precip rate, time average', datatype_flt)
    grib2_desc = grib2_var(0, 1, 54, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'prec_gsp_rate_avg', diag%prec_gsp_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars"),                  &
                & isteptype=TSTEP_AVG,                                        &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE.                                             )
    

    ! &      diag%tot_prec_rate_avg(nproma,nblks_c)
    cf_desc    = t_cf_var('tot_prec_rate_avg', 'kg m-2 s-1',                  &
      &          'total precip rate, time average', datatype_flt)
    grib2_desc = grib2_var(0, 1, 52, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_prec_rate_avg', diag%tot_prec_rate_avg,     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & isteptype=TSTEP_AVG,                                        &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & initval=0._wp, resetval=0._wp,                              &
                & action_list=actions(new_action(ACTION_RESET,precip_interval(k_jg))), &
                & lopenacc=.TRUE.                                             )
    

    ! &      diag%cape(nproma,nblks_c)
    cf_desc    = t_cf_var('cape', 'J kg-1 ', 'conv avail pot energy', datatype_flt)
    grib2_desc = grib2_var(0, 7, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)      &
    &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
    CALL add_var( diag_list, 'cape', diag%cape,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & in_group=groups("additional_precip_vars"),                  &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                & lopenacc=.TRUE. )
    

    ! &      diag%cape_ml(nproma,nblks_c)
    ! typeOfLevel ZA_SURFACE is changed to 192 in vlistDefVarIntKey
    cf_desc    = t_cf_var('cape_ml', 'J kg-1 ', 'cape of mean surface layer parcel', datatype_flt)
    grib2_desc = grib2_var(0, 7, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
                &           + t_grib2_int_key("typeOfFirstFixedSurface", 192)
    CALL add_var( diag_list, 'cape_ml', diag%cape_ml,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                & lopenacc=.TRUE. )
    
    
    ! &      diag%cin_ml(nproma,nblks_c)
    ! typeOfLevel ZA_SURFACE is changed to 192 in vlistDefVarIntKey
    cf_desc    = t_cf_var('cin_ml', 'J kg-1 ', 'convective inhibition of mean surface layer parcel', datatype_flt)
    grib2_desc = grib2_var(0, 7, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
         &           + t_grib2_int_key("typeOfFirstFixedSurface", 192)
    CALL add_var( diag_list, 'cin_ml', diag%cin_ml,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE.,                            &
!!$                & lmiss=.TRUE., missval=-999.9_wp,                            &
                & hor_interp=create_hor_interp_metadata(                      &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                & lopenacc=.TRUE. )
    

    ! &      diag%gust10(nproma,nblks_c)
    CALL getPTStringFromMS(NINT(1000*gust_interval(k_jg), i8), gust_int)
    cf_desc    = t_cf_var('gust10', 'm s-1 ', 'gust at 10 m since end of previous full '// &
                          TRIM(gust_int(3:))//' since model start', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'gust10', diag%gust10,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.TRUE., in_group=groups("pbl_vars"),   &
                & isteptype=TSTEP_MAX,                                           &
                & initval=0._wp, resetval=0._wp,                                 &
                & action_list=actions(new_action(ACTION_RESET, TRIM(gust_int))), &
                lopenacc=.TRUE.                                                  )
    

    ! &      diag%dyn_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('dyn_gust', 'm s-1 ', 'dynamical gust', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'dyn_gust', diag%dyn_gust,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.TRUE., isteptype=TSTEP_INSTANT,     &
                & hor_interp=create_hor_interp_metadata(                       &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,             &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),             &
                & loutput=.TRUE., lopenacc=.TRUE.                              )
    

    ! &      diag%con_gust(nproma,nblks_c)
    cf_desc    = t_cf_var('con_gust', 'm s-1 ', 'convective contribution to wind gust', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'con_gust', diag%con_gust,                        &
                & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,  &
                & ldims=shape2d, lrestart=.TRUE., isteptype=TSTEP_INSTANT,     &
                & loutput=.TRUE., lopenacc=.TRUE.                              )
    
   
    ! &      diag%rain_upd(nproma,nblks_c)
    cf_desc    = t_cf_var('rain_upd', 'kg m-2 s-1', 'rain in updroughts', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rain_upd', diag%rain_upd,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. )
    

    ! &      diag%con_udd(nproma,nlev,nblks,8)
    cf_desc    = t_cf_var('con_udd', 'unit ', 'convective up/downdraft fields', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'con_udd', diag%con_udd,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=(/nproma,klev,kblks,n_updown/),                       &
                & lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! &      diag%mbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mbas_con', '', 'cloud base level index', datatype_flt)
    grib2_desc = grib2_var(0, 6, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'mbas_con', diag%mbas_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,           &
                & lopenacc=.TRUE.)
    

    ! &      diag%mtop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('mtop_con', '', 'cloud top level index', datatype_flt)
    grib2_desc = grib2_var(0, 6, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'mtop_con', diag%mtop_con,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,           &
                & lopenacc=.TRUE.)
    

    ! &      diag%locum(nproma,nblks_c)
    cf_desc    = t_cf_var('locum', '', 'convective activity indicator', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'locum', diag%locum,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,           &
                & lopenacc=.TRUE.)
    

    ! &      diag%ldshcv(nproma,nblks_c)
    cf_desc    = t_cf_var('ldshcv', '', 'shallow convection indicator', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ldshcv', diag%ldshcv,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,           &
                & lopenacc=.TRUE. )
    

    ! &      diag%ktype(nproma,nblks_c)
    cf_desc    = t_cf_var('ktype', '', 'type of convection', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ktype', diag%ktype,                              &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ),                    &
                & lopenacc=.TRUE.)
    

    ! &      diag%k650(nproma,nblks_c)
    cf_desc    = t_cf_var('k650', '', 'level index corresponding to the HAG of the 650hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k650', diag%k650,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ) )

    ! &      diag%k850(nproma,nblks_c)
    cf_desc    = t_cf_var('k850', '', 'level index corresponding to the HAG of the 850hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k850', diag%k850,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ), lopenacc=.TRUE. )
    

    ! &      diag%k950(nproma,nblks_c)
    cf_desc    = t_cf_var('k950', '', 'level index corresponding to the HAG of the 950hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k950', diag%k950,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ),                    &
                & lopenacc=.TRUE.  )
    

    ! &      diag%k800(nproma,nblks_c)
    cf_desc    = t_cf_var('k800', '', 'level index corresponding to the HAG of the 800hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k800', diag%k800,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ),                    &
                & lopenacc=.TRUE.  )
    

    ! &      diag%k700(nproma,nblks_c)
    cf_desc    = t_cf_var('k700', '', 'level index corresponding to the HAG of the 700hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k700', diag%k700,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ), lopenacc=.TRUE. )
    

    ! &      diag%k400(nproma,nblks_c)
    cf_desc    = t_cf_var('k400', '', 'level index corresponding to the HAG of the 400hPa level', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'k400', diag%k400,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,  &
                & isteptype=TSTEP_CONSTANT,                                    &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ),                    &
                & lopenacc=.TRUE.  )
    


    ! &      diag%ktop_envel(nproma,nblks_c)
    cf_desc    = t_cf_var('ktop_envel', '', 'level index of upper boundary of SSO envelope layer', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ktop_envel', diag%ktop_envel,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,                 &
                & grib2_desc,ldims=shape2d, lrestart=.TRUE., loutput=.FALSE.,  &
                & lopenacc=.TRUE.  )
    

    !        diag%snowlmt(nproma,nblks_c)
    cf_desc    = t_cf_var('snowlmt', 'm', 'Height of snow fall limit above MSL', &
      &                   datatype_flt)
    grib2_desc = grib2_var(0, 1, 204, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'snowlmt', diag%snowlmt,                          &
                & GRID_UNSTRUCTURED_CELL, ZA_ISOTHERM_ZERO, cf_desc,           &
                & grib2_desc,ldims=shape2d, lrestart=.FALSE.,                  &
                & loutput=.TRUE.,                                              &
                & lmiss=.TRUE., missval=-999._wp,                              &
                & hor_interp=create_hor_interp_metadata(                       &
                &    hor_intp_type=HINTP_TYPE_LONLAT_NNB ),                    &
                & lopenacc=.TRUE.  )
    


    ! &      diag%hpbl(nproma,nblks_c)
    ! height of bdl for use in non-LES simulations

    IF ( var_in_output%hpbl ) THEN


      cf_desc    = t_cf_var('hpbl', 'm', 'boundary layer height above sea level', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hpbl', diag%hpbl,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

    ENDIF



    ! &      diag%clc(nproma,nlev,nblks_c)
    cf_desc      = t_cf_var('clc', '',  'cloud cover', datatype_flt)
    new_cf_desc  = t_cf_var('clc', '%', 'cloud cover', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    lrestart     = (atm_phy_nwp_config(k_jg)%inwp_turb == ivdiff .OR. lart)
    CALL add_var( diag_list, 'clc', diag%clc,                                 &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
      & ldims=shape3d, lrestart=lrestart,                                     &
      & in_group=groups("cloud_diag"),                                        &
      & vert_interp=create_vert_interp_metadata(                              &
      &             vert_intp_type=vintp_types("P","Z","I"),                  &
      &             vert_intp_method=VINTP_METHOD_LIN,                        &
      &             l_loglin=.FALSE.,                                         &
      &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
      &             lower_limit=0._wp ),                                      &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc),                                  &
      & lopenacc=.TRUE.  )
    

    ! &      diag%clc_rad(nproma,nlev,nblks_c)
    IF ( atm_phy_nwp_config(k_jg)%luse_clc_rad ) THEN

      ! clc_rad is the modified cloud cover used in radiation schemes and RTTOV,
      ! if reff and qr, qs, qg are active in radiation.
      ! In this case, it differs from the above "normal" clc in that
      ! it is set to 1.0 for gridpoints where qr, qg are present, as required by
      ! the radiative transfer schemes of RTTOV, RRTM and ECRAD.
      ! We need a separate clc_rad for this, because this modification is not
      ! wanted in the "normal" clc above, which is the basis for diagnostics
      ! like ceiling, clct, clch, clcm, clcl. For example,
      ! if it rains, the ceiling should be where the cloud base (qc) is and not
      ! at the ground, where only raindrops are present.
      ! Currently saved with shortname DUMMY_1 in Grib Files.
           
      cf_desc      = t_cf_var('clc_rad', '',  'cloud cover for radiation scheme', datatype_flt)
      new_cf_desc  = t_cf_var('clc_rad', '%', 'cloud cover for radiation scheme', datatype_flt)
      grib2_desc   = grib2_var(0, 254, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'clc_rad', diag%clc_rad,                            &
           & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
           & ldims=shape3d, lrestart=.FALSE.,                                      &
           & vert_interp=create_vert_interp_metadata(                              &
           &             vert_intp_type=vintp_types("P","Z","I"),                  &
           &             vert_intp_method=VINTP_METHOD_LIN,                        &
           &             l_loglin=.FALSE.,                                         &
           &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
           &             lower_limit=0._wp ),                                      &
           & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
           &                 new_cf=new_cf_desc),                                  &
           & lopenacc=.TRUE.  )
      

    END IF


    ! &      diag%clct(nproma,nblks_c)
    cf_desc      = t_cf_var('clct', '',  'total cloud cover', datatype_flt)
    new_cf_desc  = t_cf_var('clct', '%', 'total cloud cover', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clct', diag%clct,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & in_group=groups("additional_precip_vars"),                            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc),                                  &
      & lopenacc=.TRUE.  )
    

    ! &      diag%clct_mod(nproma,nblks_c)
    cf_desc      = t_cf_var('clct_mod', '', 'modified total cloud cover for media', &
      &                     datatype_flt)
    grib2_desc   = grib2_var(0, 6, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clct_mod', diag%clct_mod,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & in_group=groups("additional_precip_vars"),                            &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE.  )
    

    ! &      diag%clch(nproma,nblks_c)
    cf_desc      = t_cf_var('clch', '', 'high level clouds',  datatype_flt)
    new_cf_desc  = t_cf_var('clch', '%', 'high level clouds', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clch', diag%clch,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_0, cf_desc, grib2_desc,           &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc),                                  &
      & lopenacc=.TRUE.  )
    

    ! &      diag%clcm(nproma,nblks_c)
    cf_desc      = t_cf_var('clcm', '',  'mid level clouds', datatype_flt)
    new_cf_desc  = t_cf_var('clcm', '%', 'mid level clouds', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'clcm', diag%clcm,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_400, cf_desc, grib2_desc,         &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc),                                  &
      & lopenacc=.TRUE.  )
    

    ! &      diag%clcl(nproma,nblks_c)
    cf_desc      = t_cf_var('clcl', '',  'low level clouds', datatype_flt)
    new_cf_desc  = t_cf_var('clcl', '%', 'low level clouds', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 1)
    CALL add_var( diag_list, 'clcl', diag%clcl,                               &
      & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_800, cf_desc, grib2_desc,         &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp,                          &
      &                 new_cf=new_cf_desc),                                  &
      & lopenacc=.TRUE.  )
    

    ! &      diag%cldepth(nproma,nblks_c)
    cf_desc      = t_cf_var('cldepth', '',  'modified cloud depth for media', datatype_flt)
    grib2_desc   = grib2_var(0, 6, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cldepth', diag%cldepth,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE.  )
    

    ! &      diag%fac_ccqc(nproma,nblks_c)
    cf_desc    = t_cf_var('fac_ccqc', ' ','factor for cloud cover - cloud water relationship', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'fac_ccqc', diag%fac_ccqc,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE., initval=1._wp, lopenacc=.TRUE. )
    

    ! &      diag%hbas_con(nproma,nblks_c)
    cf_desc    = t_cf_var('hbas_con', 'm', 'height of convective cloud base', datatype_flt)
    grib2_desc = grib2_var(0, 6, 26, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'hbas_con', diag%hbas_con,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc,           &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
!!$      & lmiss=.TRUE., missval=-500._wp,                                       &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB), &
      & lopenacc=.TRUE.)
    

    ! &      diag%htop_con(nproma,nblks_c)
    cf_desc    = t_cf_var('htop_con', 'm', 'height of convective cloud top', datatype_flt)
    grib2_desc = grib2_var(0, 6, 27, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'htop_con', diag%htop_con,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
!!$      & lmiss=.TRUE., missval=-500._wp,                                       &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB), &
      & lopenacc=.TRUE.)
    


    ! &      diag%htop_dc(nproma,nblks_c)
    cf_desc    = t_cf_var('htop_dc', 'm', 'height of top of dry convection', datatype_flt)
    grib2_desc = grib2_var(0, 6, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    CALL add_var( diag_list, 'htop_dc', diag%htop_dc,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE.,                                      &
!!$      & lmiss=.TRUE., missval=-999._wp,                                       &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB), &
      & lopenacc=.TRUE.)
    

    ! &      diag%acdnc(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('acdnc', 'm-3', 'cloud droplet number concentration', datatype_flt)
    grib2_desc = grib2_var(0, 6, 30, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'acdnc', diag%acdnc,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,               &
      & ldims=shape3d, lrestart=.FALSE.,                                      &
      & isteptype=TSTEP_CONSTANT, lopenacc=.TRUE. )
    

    IF (ldass_lhn) THEN
      ! &      diag%tt_lheat(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('tt_lheat', 'K s-1',                &
        &          '3d latent heat relaese', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'tt_lheat', diag%tt_lheat,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    
  
      ! &      diag%qrs_flux(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('qrs_flux', 'kg m-2 s-1',                &
        &          '3d precipitation flux', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'qrs_flux', diag%qrs_flux,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    
  
      ! &      diag%lhn_diag(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('lhn_diag', '-',                &
        &          'diagnose of LHN', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'lhn_diag', diag%lhn_diag,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB), &
                  & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
      
  
      ! &      diag%lhn_diag(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('ttend_lhn', 'K s-1',                &
        &          'tempature increment', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'ttend_lhn', diag%ttend_lhn,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
      

      ! &      diag%lhn_diag(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('qvtend_lhn', 'g/g s-1',                &
        &          'qv increment', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 255, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'qvtend_lhn', diag%qvtend_lhn,                       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & lrestart = .FALSE., & ! .TRUE. may be necessary for ART (to be evaluated)
                  & ldims=shape3d,                                              &  
                  & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                  &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
      
    ELSE
      ALLOCATE (diag%qrs_flux(1,1,kblks))
      !$ACC ENTER DATA CREATE(diag%qrs_flux)
    ENDIF


    !      diag%tot_cld(nproma,nlev,nblks_c,3)
    cf_desc    = t_cf_var('tot_cld', ' ','total cloud variables (qv,qc,qi)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld', diag%tot_cld,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=(/nproma,klev,kblks,3/) ,                               &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                & initval=0.0_wp,                                               &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
                &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                & lopenacc=.TRUE.  )
    

    ! Diagnostic variables for any of the stochastic convection schemes
    IF (atm_phy_nwp_config(k_jg)%lstoch_sde .or. atm_phy_nwp_config(k_jg)%lstoch_expl &
         & .or. atm_phy_nwp_config(k_jg)%lstoch_deep) THEN    

      ! &      diag%mf_b(nproma,nblks_c)
      cf_desc    = t_cf_var('mf_b', 'm s-1 ', 'bulk mass flux ', DATATYPE_FLT32)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mf_b', diag%mf_b,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,lrestart=.FALSE.)
      ! &      diag%mf_p(nproma,nblks_c)
      cf_desc    = t_cf_var('mf_p', 'm s-1 ', 'perturbed mass flux ', DATATYPE_FLT32)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mf_p', diag%mf_p,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,lrestart=.FALSE.)
      ! &      diag%mf_num(nproma,nblks_c)
      cf_desc    = t_cf_var('mf_num', '- ', 'number of clouds per grid box ', DATATYPE_FLT32)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mf_num', diag%mf_num,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape2d,lrestart=.FALSE.)
   ELSE
      ALLOCATE(diag%mf_num(0,kblks), diag%mf_b(0,kblks), diag%mf_p(0,kblks))
   ENDIF

    ALLOCATE( diag%tot_ptr(kcloud))

    !QV
    CALL add_ref( diag_list, 'tot_cld',                                            &
                & tot_pfx//'qv_dia', diag%tot_ptr(iqv)%p_3d,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                & t_cf_var(tot_pfx//'qv_dia', 'kg kg-1',                           &
                &          'total specific humidity (diagnostic)', datatype_flt),  &
                & grib2_var(0, 1, 211, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ref_idx=iqv,                                                     &
                & ldims=shape3d,                                                   &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P", "Z", "I"),           &
                &             vert_intp_method=VINTP_METHOD_QV,                    &
                &             l_satlimit=.FALSE.,                                  & 
                &             lower_limit=2.5e-7_wp, l_restore_pbldev=.FALSE. ),   &
                & in_group=groups("cloud_diag") )
    

    !QC
    CALL add_ref( diag_list, 'tot_cld',                                            &
                & tot_pfx//'qc_dia', diag%tot_ptr(iqc)%p_3d,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                & t_cf_var(tot_pfx//'qc_dia', 'kg kg-1',                           &
                & 'total specific cloud water content (diagnostic)', datatype_flt),&
                & grib2_var(0, 1, 212, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ref_idx=iqc,                                                     &
                & ldims=shape3d,                                                   &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P","Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_LIN,                   &
                &             l_loglin=.FALSE.,                                    &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,              &
                &             lower_limit=0._wp ),                                 &
                & in_group=groups("cloud_diag") )
    

    !QI
    CALL add_ref( diag_list, 'tot_cld',                                            &
                & tot_pfx//'qi_dia', diag%tot_ptr(iqi)%p_3d,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                & t_cf_var(tot_pfx//'qi_dia', 'kg kg-1',                           &
                & 'total specific cloud ice content (diagnostic)', datatype_flt),  &
                & grib2_var(0, 1, 213, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
                & ref_idx=iqi,                                                     &
                & ldims=shape3d,                                                   &
                & vert_interp=create_vert_interp_metadata(                         &
                &             vert_intp_type=vintp_types("P","Z","I"),             &
                &             vert_intp_method=VINTP_METHOD_LIN,                   &
                &             l_loglin=.FALSE.,                                    &
                &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,              &
                &             lower_limit=0._wp ),                                 &
                & in_group=groups("cloud_diag") )
    

    !      diag%tot_cld_vi(nproma,nblks_c,3)
    cf_desc     = t_cf_var('tot_cld_vi', 'kg m-2','vertical integr total cloud variables', datatype_flt)
    grib2_desc   = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tot_cld_vi', diag%tot_cld_vi,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                &                                 ldims=(/nproma,kblks,3/),   &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       &
                & lopenacc=.TRUE.  )
    

    ! fill the seperate variables belonging to the container tot_cld_vi
    ALLOCATE( diag%tci_ptr(3))
       
    !TQV_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                        &
      & 'tqv_dia', diag%tci_ptr(iqv)%p_2d,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                       &
      & t_cf_var('tqv_dia', 'kg m**-2',                           &
      & 'total column integrated water vapour (diagnostic)', datatype_flt),   &
      & grib2_var( 0, 1, 214, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
      & ref_idx=iqv,                                              &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !TQC_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & 'tqc_dia', diag%tci_ptr(iqc)%p_2d,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
      & t_cf_var('tqc_dia', 'kg m**-2',                            &
      & 'total column integrated cloud water (diagnostic)', datatype_flt),    &
      & grib2_var( 0, 1, 215, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
      & ref_idx=iqc,                                               &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    !TQI_DIA
    CALL add_ref( diag_list, 'tot_cld_vi',                         &
      & 'tqi_dia', diag%tci_ptr(iqi)%p_2d,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
      & t_cf_var('tqi_dia', 'kg m**-2',                            &
      & 'total column integrated cloud ice (diagnostic)', datatype_flt),      &
      & grib2_var(0, 1, 216, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
      & ref_idx=iqi,                                               &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("additional_precip_vars"))

    ! &      diag%qc_sgs(nproma,nlev,nblks_c)
    cf_desc      = t_cf_var('qc_sgs', 'kg kg-1',  'subgrid-scale cloud water', datatype_flt)
    grib2_desc   = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qc_sgs', diag%qc_sgs,                          &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,           &
      & ldims=shape3d, lrestart=.TRUE.,                                      &
      & in_group=groups("cloud_diag"),                                       &
      & lopenacc=.TRUE.  )
     


    !------------------
    ! Radiation
    !------------------
    ! 2D variables

    !        diag%albdif    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albdif', '', 'Shortwave albedo for diffuse radiation', &
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albdif', '%','Shortwave albedo for diffuse radiation', &
      &                    datatype_flt)
    grib2_desc  = grib2_var(0, 19, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albdif', diag%albdif,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, in_group=groups("rad_vars"),                         &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),   &
      & lopenacc=.TRUE.)
      

    !        diag%albvisdif    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albvisdif', '', 'UV visible albedo for diffuse radiation', &
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albvisdif', '%','UV visible albedo for diffuse radiation', &
      &                    datatype_flt)
    grib2_desc  = grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albvisdif', diag%albvisdif,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, in_group=groups("rad_vars"),                         &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),   &
      & lopenacc=.TRUE.)
      

    !        diag%albvisdir    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albvisdir', '', 'UV visible albedo for direct radiation', &
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albvisdir', '%','UV visible albedo for direct radiation', &
      &                    datatype_flt)
    grib2_desc  = grib2_var(192, 128, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albvisdir', diag%albvisdir,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),   &
      & lopenacc=.TRUE.)
      

    !        diag%albnirdif    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albnirdif', '',  'Near IR albedo for diffuse radiation',&
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albnirdif', '%', 'Near IR albedo for diffuse radiation',&
      &                    datatype_flt)
    grib2_desc  = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albnirdif', diag%albnirdif,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, in_group=groups("rad_vars"),                         &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),   &
      & lopenacc=.TRUE.)
      

    !        diag%albnirdir    (nproma,       nblks),          &
    cf_desc     = t_cf_var('albnirdir', '',  'Near IR albedo for direct radiation',&
      &                    datatype_flt)
    new_cf_desc = t_cf_var('albnirdir', '%', 'Near IR albedo for direct radiation',&
      &                    datatype_flt)
    grib2_desc  = grib2_var(192, 128, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'albnirdir', diag%albnirdir,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),   &
      & lopenacc=.TRUE.)
      


    ! longwave surface emissivity
    !
    ! lw_emiss     diag%lw_emiss(nproma,nblks_c)
    cf_desc    = t_cf_var('lw_emiss', '-', 'longwave surface emissivity', datatype_flt)
    grib2_desc = grib2_var( 2, 3, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lw_emiss', diag%lw_emiss,         &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
      &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
      &           lrestart=.TRUE., lopenacc=.TRUE.)
      


    IF (icpl_da_snowalb >= 1) THEN
      ! Factor for adaptive snow albedo tuning
      !
      ! snowalb_fac     diag%snowalb_fac(nproma,nblks_c)
      cf_desc    = t_cf_var('snowalb_fac', '-', 'tuning factor for snow albedo', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'snowalb_fac', diag%snowalb_fac,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
        &           initval=1.0_wp, lrestart=.TRUE., lopenacc=.TRUE.)
        
    ENDIF

    IF (icpl_da_seaice >= 2) THEN
      ! Factor for adaptive bottom heat flux tuning
      !
      ! hflux_si_fac     diag%hflux_si_fac(nproma,nblks_c)
      cf_desc    = t_cf_var('hflux_si_fac', '-', 'tuning factor for seaice bottom heat flux', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hflux_si_fac', diag%hflux_si_fac, &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
        &           initval=0.0_wp, lrestart=.TRUE., lopenacc=.TRUE.)
        
    ENDIF

    IF (icpl_da_skinc >= 2) THEN
      ! Factors for adaptive tuning of soil heat capacity and conductivity
      !
      ! heatcap_fac     diag%heatcap_fac(nproma,nblks_c)
      cf_desc    = t_cf_var('heatcap_fac', '-', 'tuning factor for soil heat capacity', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'heatcap_fac', diag%heatcap_fac,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
        &           initval=1.0_wp, lrestart=.TRUE., lopenacc=.TRUE.)
        

      ! heatcond_fac     diag%heatcond_fac(nproma,nblks_c)
      cf_desc    = t_cf_var('heatcond_fac', '-', 'tuning factor for soil heat conductivity', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'heatcond_fac', diag%heatcond_fac,   &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
        &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
        &           initval=1.0_wp, lrestart=.TRUE., lopenacc=.TRUE.)
        

    ENDIF

    ! Factor for adaptive surface friction tuning
    !
    ! sfcfric_fac     diag%sfcfric_fac(nproma,nblks_c)
    cf_desc    = t_cf_var('sfcfric_fac', '-', 'tuning factor for surface friction', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sfcfric_fac', diag%sfcfric_fac,   &
      &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  &
      &           grib2_desc, ldims=shape2d, loutput=.TRUE.,    &
      &           initval=1.0_wp, lrestart=.TRUE., lopenacc=.TRUE.)
      

    ! These variables only make sense if the land-surface scheme is switched on.
    IF ( atm_phy_nwp_config(k_jg)%inwp_surface > 0 ) THEN

      !        diag%albdif_t (nproma, nblks, ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('albdif_t', '', &
        &                   'tile-based shortwave albedo for diffusive radiation',&
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 19, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'albdif_t', diag%albdif_t,                &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE.,       &
        & loutput=.FALSE., lopenacc=.TRUE.)
      

      ! fill the seperate variables belonging to the container albdif_t
      ALLOCATE(diag%albdif_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'albdif_t',                               &
           & 'albdif_t_'//TRIM(ADJUSTL(csfc)),                             &
           & diag%albdif_t_ptr(jsfc)%p_2d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('albdif_t_'//csfc, '', '', datatype_flt),    &
           & grib2_var(0, 19, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
           & ref_idx=jsfc, ldims=shape2d,                                  &
           & var_class=CLASS_TILE,                                         &
           & lrestart=.TRUE.                  )
      ENDDO

      !        diag%albvisdif_t (nproma, nblks, ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('albvisdif_t', '', &
        &                   'tile-based UV visible albedo for diffusive radiation',&
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'albvisdif_t', diag%albvisdif_t,               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE., &
        & loutput=.FALSE., lopenacc=.TRUE.)
        

      ! fill the seperate variables belonging to the container albvisdif_t
      ALLOCATE(diag%albvisdif_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'albvisdif_t',                            &
           & 'albvisdif_t_'//TRIM(ADJUSTL(csfc)),                          &
           & diag%albvisdif_t_ptr(jsfc)%p_2d,                              &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('albvisdif_t_'//csfc, '', '', datatype_flt), &
           & grib2_var(0, 19, 222, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
           & ref_idx=jsfc, ldims=shape2d,                                  &
           & var_class=CLASS_TILE,                                         &
           & lrestart=.TRUE.                  )
      ENDDO

      !        diag%albnirdif_t (nproma, nblks, ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('albnirdif_t', '', &
        &                   'tile-based near IR albedo for diffuse radiation',&
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'albnirdif_t', diag%albnirdif_t,               &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE., &
        & loutput=.FALSE., lopenacc=.TRUE.)
        


      ! fill the seperate variables belonging to the container albnirdif_t
      ALLOCATE(diag%albnirdif_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'albnirdif_t',                            &
           & 'albnirdif_t_'//TRIM(ADJUSTL(csfc)),                          &
           & diag%albnirdif_t_ptr(jsfc)%p_2d,                              &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('albnirdif_t_'//csfc, '', '', datatype_flt), &
           & grib2_var(0, 19, 223, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
           & ref_idx=jsfc, ldims=shape2d,                                  &
           & var_class=CLASS_TILE,                                         &
           & lrestart=.TRUE.                  )
      ENDDO


      ! &      diag%swflxsfc_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('sob_s_t', 'W m-2', 'tile-based shortwave net flux at surface', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'sob_s_t', diag%swflxsfc_t,                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE., &
        & loutput=.FALSE., lopenacc=.TRUE.)
      

      ! fill the seperate variables belonging to the container swflxsfc_t
      ALLOCATE(diag%swflxsfc_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'sob_s_t',                                &
           & 'sob_s_t_'//TRIM(ADJUSTL(csfc)),                              &
           & diag%swflxsfc_t_ptr(jsfc)%p_2d,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('swflxsfc_t_'//csfc, '', '', datatype_flt),  &
           & grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
           & ref_idx=jsfc, ldims=shape2d,                                  &
           & var_class=CLASS_TILE,                                         &
           & lrestart=.TRUE., loutput=.TRUE.,                              &
           & in_group=groups("rad_vars"))
      ENDDO

      ! &      diag%lwflxsfc_t(nproma,nblks_c,ntiles_total+ntiles_water)
      cf_desc    = t_cf_var('thb_s_t', 'W m-2', 'tile-based longwave net flux at surface', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'thb_s_t', diag%lwflxsfc_t,                 &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
        & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE.,         &
        & loutput=.FALSE., lopenacc=.TRUE.)
      

      ! fill the seperate variables belonging to the container lwflxsfc_t
      ALLOCATE(diag%lwflxsfc_t_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        WRITE(csfc,'(i1)') jsfc 
        CALL add_ref( diag_list, 'thb_s_t',                                &
           & 'thb_s_t_'//TRIM(ADJUSTL(csfc)),                              &
           & diag%lwflxsfc_t_ptr(jsfc)%p_2d,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('lwflxsfc_t_'//csfc, '', '', datatype_flt),  &
           & grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
           & ref_idx=jsfc, ldims=shape2d,                                  &
           & var_class=CLASS_TILE,                                         &
           & lrestart=.TRUE., loutput=.TRUE.,                              &
           & in_group=groups("rad_vars"))
      ENDDO
    ENDIF

    !        diag%cosmu0    (nproma,nblks)
    cf_desc    = t_cf_var('cosmu0', '-', 'Cosine of solar zenith angle', datatype_flt)
    grib2_desc = grib2_var(192, 214, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cosmu0', diag%cosmu0,                  &
    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
    & ldims=shape2d, in_group=groups("rad_vars"), lopenacc=.TRUE.)
    

    ! &      diag%tsfctrad(nproma,nblks_c)
    cf_desc    = t_cf_var('tsfctrad', 'K', 'surface temperature at trad', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tsfctrad', diag%tsfctrad,                     &
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
            & ldims=shape2d, lopenacc=.TRUE.                                )
    

    ! &       diag% flxdwswtoa(nproma,       nblks),          &
    cf_desc    = t_cf_var('sod_t', 'W m-2', 'downward shortwave flux at TOA', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 4, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sod_t', diag%flxdwswtoa,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
      & ldims=shape2d, lrestart=.FALSE.,in_group=groups("rad_vars"),     &
      & lopenacc=.TRUE. )
    

    ! &      diag%swflxsfc(nproma,nblks_c)
    cf_desc    = t_cf_var('sob_s', 'W m-2', 'shortwave net flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sob_s', diag%swflxsfc,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & in_group=groups("rad_vars"),                                        &
      & lopenacc=.TRUE. )
    

    ! &      diag%lwflxclr_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('thbclr_s', 'W m-2', 'net longwave clear-sky flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thbclr_s', diag%lwflxclr_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
      & ldims=shape2d, lopenacc=.TRUE. )
    

    ! &      diag%swflxclr_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('sobclr_s', 'W m-2', 'net shortwave clear-sky flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sobclr_s', diag%swflxclr_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
      & ldims=shape2d, lopenacc=.TRUE. )
    

    ! &      diag%swflxtoa(nproma,nblks_c)
    cf_desc    = t_cf_var('sob_t', 'W m-2', 'shortwave net flux at TOA', datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sob_t', diag%swflxtoa,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE.,                                    &
      & in_group=groups("rad_vars"), lopenacc=.TRUE. )
    

    ! &      diag%lwflxtoa(nproma,nblks_c)
    cf_desc    = t_cf_var('thb_t', 'W m-2', 'thermal net flux at TOA', datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thb_t', diag%lwflxtoa,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE.,                                    &
      & in_group=groups("rad_vars"), lopenacc=.TRUE.                        )
    

    ! &      diag%trsolclr_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsolclr_sfc', '', 'shortwave clear-sky transmissivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsolclr_sfc', diag%trsolclr_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, loutput=.FALSE., lopenacc=.TRUE.                     )
    

    ! &      diag%trsol_up_toa(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_up_toa', '', 'shortwave upward transmissivity at TOA', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_up_toa', diag%trsol_up_toa,             &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.   )
    

    ! &      diag%trsol_up_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_up_sfc', '', 'shortwave upward transmissivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_up_sfc', diag%trsol_up_sfc,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE., lopenacc=.TRUE.   )
    

    ! &      diag%trsol_nir_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_nir_sfc', '', 'near-infrared downward transmissivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_nir_sfc', diag%trsol_nir_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE., lopenacc=.TRUE.    )
    

    ! &      diag%trsol_vis_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_vis_sfc', '', 'visible downward transmissivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_vis_sfc', diag%trsol_vis_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE., lopenacc=.TRUE.    )
    

    ! &      diag%trsol_par_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_par_sfc', '', 'photosynthetically active downward transmissivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_par_sfc', diag%trsol_par_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE., lopenacc=.TRUE.    )
    

    ! &      diag%trsol_dn_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('trsol_dn_sfc_diff', '', 'shortwave diffuse downward transmissivity at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsol_dn_sfc_diff', diag%trsol_dn_sfc_diff,   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., loutput=.FALSE., lopenacc=.TRUE.    )
    

    ! &      diag%swflx_up_toa(nproma,nblks_c)
    cf_desc    = t_cf_var('sou_t', 'W m-2', 'shortwave upward flux at TOA', datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sou_t', diag%swflx_up_toa,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,                &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars"),       &
      & lopenacc=.TRUE.                                                     )
    

    ! &      diag%swflx_up_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('sou_s', 'W m-2', 'shortwave upward flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sou_s', diag%swflx_up_sfc,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars"),       &
      & lopenacc=.TRUE. )
    

    ! &      diag%swflx_nir_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('swflx_nir_sfc', 'W m-2', 'downward near-infrared flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_nir_sfc', diag%swflx_nir_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE.)
    

    ! &      diag%swflx_vis_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('swflx_vis_sfc', 'W m-2', 'downward visible flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_vis_sfc', diag%swflx_vis_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE.)
    

    ! &      diag%swflx_par_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('swflx_par_sfc', 'W m-2', 'downward photosynthetically active flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_par_sfc', diag%swflx_par_sfc,           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE.)
    

    ! &      diag%fr_nir_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('fr_nir_sfc_diff', '', 'diffuse fraction of downward near-infrared flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'fr_nir_sfc_diff', diag%fr_nir_sfc_diff,       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE.)
    

    ! &      diag%fr_vis_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('fr_vis_sfc_diff', '', 'diffuse fraction of downward visible flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'fr_vis_sfc_diff', diag%fr_vis_sfc_diff,       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE.)
    

    ! &      diag%fr_par_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('fr_par_sfc_diff', '', 'diffuse fraction of downward photosynthetically active flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'fr_par_sfc_diff', diag%fr_par_sfc_diff,       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d, lrestart=.TRUE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE.)
    

    ! &      diag%swflx_dn_sfc_diff(nproma,nblks_c)
    cf_desc    = t_cf_var('sodifd_s', 'W m-2', 'shortwave diffuse downward flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 4, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sodifd_s', diag%swflx_dn_sfc_diff,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("rad_vars"),        &
      & lopenacc=.TRUE. )
    

    ! &      diag%lwflxsfc(nproma,nblks_c)
    cf_desc    = t_cf_var('thb_s', 'W m-2', 'longwave net flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thb_s', diag%lwflxsfc,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,                                                      &
      & in_group=groups("rad_vars"), lopenacc=.TRUE. )
    

    ! This is an auxiliary storage field needed because the final diagnostic quantity (below)
    ! is updated each physics time step following the time evolution of the ground temperature
    ! &      diag%lwflx_up_sfc_rs(nproma,nblks_c)
    cf_desc    = t_cf_var('lwflx_up_sfc_rs', 'W m-2', 'longwave upward flux at surface', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflx_up_sfc_rs', diag%lwflx_up_sfc_rs,       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,  lrestart=.TRUE., loutput=.FALSE., lopenacc=.TRUE.   )
    

    ! &      diag%lwflx_up_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('thu_s', 'W m-2', 'longwave upward flux at surface', datatype_flt)
    grib2_desc = grib2_var(0, 5, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'thu_s', diag%lwflx_up_sfc,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
      & ldims=shape2d,  lrestart=.TRUE.,                                    &
      & in_group=groups("rad_vars"), lopenacc=.TRUE. )
    

    IF (lflux_avg ) THEN
        prefix = "a"
        meaning = "mean"
        varunits= "W/m**2"
        a_steptype= TSTEP_AVG
    ELSE
        prefix = "acc"
        meaning = "acc." 
        varunits= "J/m**2"
        a_steptype= TSTEP_ACCUM
    END IF
    WRITE(name,'(A,A5)') TRIM(prefix),"thb_t"
    WRITE(long_name,'(A26,A4,A18)') "TOA net thermal radiation ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%lwflxtoa_a,                  &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,               &
      & ldims=shape2d, isteptype=a_steptype, in_group=groups("rad_vars"),  &
      & hor_interp=create_hor_interp_metadata(                             &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                           &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                          &
      & lopenacc=.TRUE.)
    


    ! &      diag%swflxtoa_a(nproma,nblks_c)
    WRITE(name,'(A,A5)') TRIM(prefix),"sob_t"
    WRITE(long_name,'(A26,A4,A18)') "TOA net solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name) , diag%swflxtoa_a,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,               &
      & ldims=shape2d,                                                     &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                 &
      & hor_interp=create_hor_interp_metadata(                             &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                           &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                          &
      & lopenacc=.TRUE. )
    
    
    ! &      diag%lwflxsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A5)') TRIM(prefix),"thb_s"
    WRITE(long_name,'(A30,A4,A18)') "surface net thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%lwflxsfc_a,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE. )
    

    ! &      diag%lwflxclrsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A8)') TRIM(prefix),"thbclr_s"
    WRITE(long_name,'(A40,A4,A18)') "clear-sky surface net thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 5, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%lwflxclrsfc_a,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, isteptype=a_steptype,                              &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE. )
    

    ! &      diag%swflxsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A5)') TRIM(prefix),"sob_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface net solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%swflxsfc_a,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE.)
    

    ! &      diag%swflxclrsfc_a(nproma,nblks_c)
    WRITE(name,'(A,A8)') TRIM(prefix),"sobclr_s"
    WRITE(long_name,'(A40,A4,A18)') "Clear-sky surface net solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%swflxclrsfc_a,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, isteptype=a_steptype,                              &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE.)
    

    ! &       diag%asod_t(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"sod_t"
    WRITE(long_name,'(A30,A4,A18)') "Top down solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asod_t,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,              &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("rad_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE. )
    

    ! &       diag%asou_t(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"sou_t"
    WRITE(long_name,'(A30,A4,A18)') "Top up solar radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asou_t,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_TOA, cf_desc, grib2_desc,             &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &       diag%athd_s(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"thd_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 5, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%athd_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &       diag%athu_s(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"thu_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface up thermal radiation ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 5, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%athu_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &       diag%asodird_s(nproma,nblks)
    WRITE(name,'(A,A8)') TRIM(prefix),"sodird_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down solar direct rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asodird_s,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &       diag%asodifd_s(nproma,nblks)
    WRITE(name,'(A,A8)') TRIM(prefix),"sodifd_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down solar diff. rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asodifd_s,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &       diag%asod_s(nproma,nblks)
    WRITE(name,'(A,A5)') TRIM(prefix),"sod_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface down solar rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asod_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d, isteptype=a_steptype,                             &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &       diag%asodifu_s(nproma,nblks)
    WRITE(name,'(A,A8)') TRIM(prefix),"sodifu_s"
    WRITE(long_name,'(A30,A4,A18)') "Surface up solar diff. rad. ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%asodifu_s,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   &
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &      diag%aswflx_par_sfc(nproma,nblks_c)
    WRITE(name,'(A,A13)') TRIM(prefix),"swflx_par_sfc"
    WRITE(long_name,'(A30,A4,A18)') "Downward PAR flux ", meaning, &
                                  &" since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 4, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%aswflx_par_sfc,            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,         &
      & ldims=shape2d,                                                   & 
      & isteptype=a_steptype, in_group=groups("rad_vars"),               &
      & hor_interp=create_hor_interp_metadata(                           &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                         &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                        &
      & lopenacc=.TRUE. )
    

    ! &      diag%vio3(nproma,nblks_c)
    cf_desc    = t_cf_var('vio3', '', 'vertically integrated ozone amount', datatype_flt)
    grib2_desc = grib2_var(0, 14, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'vio3', diag%vio3,                           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. )
    

    ! &      diag%hmo3(nproma,nblks_c)
    cf_desc    = t_cf_var('hmo3', 'Pa', 'height of O3 maximum (Pa)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'hmo3', diag%hmo3,                           &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. )
    

    ! Reference pressure used for vertical distribution of aerosol optical depths
    ! &      diag%pref_aerdis(nproma,nblks_c)
    cf_desc    = t_cf_var('pref_aerdis', '',                                                             &
      &                   'Reference pressure used for vertical distribution of aerosol optical depths', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'pref_aerdis', diag%pref_aerdis,              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. ) 
    


    IF (irad_aero == iRadAeroTegen .OR. irad_aero == iRadAeroART) THEN ! Tegen aerosol climatology, time-interpolated values 
                                                    ! (needed as state fields for coupling with microphysics and convection)
      IF (atm_phy_nwp_config(k_jg)%icpl_aero_gscp > 1 .OR. icpl_aero_conv > 1 .OR. iprog_aero > 0) THEN
        lrestart = .TRUE.
      ELSE
        lrestart = .FALSE.
      ENDIF

      ! &      diag%aercl_ss(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_ss', '', 'sea salt aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_ss', diag%aercl_ss,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
      

      ! &      diag%aercl_or(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_or', '', 'organic aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_or', diag%aercl_or,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
      

      ! &      diag%aercl_bc(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_bc', '', 'black carbon aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_bc', diag%aercl_bc,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
       

      ! &      diag%aercl_su(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_su', '', 'total sulfate aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_su', diag%aercl_su,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
       

      ! &      diag%aercl_du(nproma,nblks_c)
      cf_desc    = t_cf_var('aercl_du', '', 'total soil dust aerosol climatology', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aercl_du', diag%aercl_du,                       &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
       

      ! &      diag%aerosol(nproma,nclass_aero,nblks_c)
      cf_desc    = t_cf_var('aerosol', '', '', DATATYPE_FLT32)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'aerosol', diag%aerosol,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape3d_aero, lrestart=.FALSE., lcontainer=.TRUE.,          &
        & loutput=.FALSE.,                                                  &
        & in_group=groups("iau_restore_vars"),                              &
        & lopenacc=.TRUE. )
      

      ALLOCATE(diag%aerosol_ptr(nclass_aero))
      DO k = 1, nclass_aero
        SELECT CASE (k)
        CASE (iss)
          caer='ss'
          constituentType = 62008
          cf_desc    = t_cf_var('aer_'//caer, '', &
            &                   'sea salt aerosol', DATATYPE_FLT32)
        CASE (iorg)
          caer='or'
          constituentType = 62010
          cf_desc    = t_cf_var('aer_'//caer, '', &
            &                   'organic aerosol', DATATYPE_FLT32)
        CASE (ibc)
          caer='bc'
          constituentType = 62009
          cf_desc    = t_cf_var('aer_'//caer, '', &
            &                   'black carbon aerosol', DATATYPE_FLT32)
        CASE (iso4)
          caer='su'
          constituentType = 62006
          cf_desc    = t_cf_var('aer_'//caer, '', &
            &                   'total sulfate aerosol', DATATYPE_FLT32)
        CASE (idu)
          caer='du'
          constituentType = 62001
          cf_desc    = t_cf_var('aer_'//caer, '', &
            &                   'total soil dust aerosol', DATATYPE_FLT32)
        END SELECT

        grib2_desc = grib2_var(0, 20, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL)   &
          &           + t_grib2_int_key("constituentType", constituentType)
        IF (iprog_aero > 1 .OR. iprog_aero == 1 .AND. k == idu) THEN
          CALL add_ref( diag_list, 'aerosol',                                    &
                & 'aer_'//TRIM(caer), diag%aerosol_ptr(k)%p_2d,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                & cf_desc,                                                       &
                & grib2_desc, ref_idx=k,                                         &
                & ldims=shape2d, lrestart=lrestart, opt_var_ref_pos = 2,         &
                & var_class=CLASS_CHEM,                                          &
                & in_group=groups("dwd_fg_sfc_vars","mode_iau_fg_in",            &
                & "mode_iau_old_fg_in","mode_dwd_fg_in")                         )
        ELSE
          CALL add_ref( diag_list, 'aerosol',                                    &
                & 'aer_'//TRIM(caer), diag%aerosol_ptr(k)%p_2d,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                & cf_desc,                                                       &
                & grib2_desc, ref_idx=k,                                         &
                & ldims=shape2d, lrestart=lrestart, opt_var_ref_pos = 2,         &
                & var_class=CLASS_CHEM                                           )
        ENDIF
      ENDDO

    ENDIF

    IF ( (irad_aero == iRadAeroTegen .OR. irad_aero == iRadAeroART) .AND.  &
      &  (atm_phy_nwp_config(k_jg)%icpl_aero_gscp == 1 .OR. icpl_aero_conv == 1) ) THEN
      lrestart = .TRUE.
    ELSE
      lrestart = .FALSE.
    ENDIF

    ! &      diag%cloud_num(nproma,nblks_c)
    cf_desc    = t_cf_var('cloud_num', 'm-3', 'cloud droplet number concentration', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'cloud_num', diag%cloud_num,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
    

    ! estimated inversion strength
    ! &      diag%conv_eis(nproma,nblks_c)
    cf_desc    = t_cf_var('conv_eis', 'm-3', 'estimated inversion strength', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'conv_eis', diag%conv_eis,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=lrestart, lopenacc=.TRUE. )
    

    !------------------
    !Radiation 3D variables

    ! &      diag%lwflxall(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('lwflxall', 'W m-2 ', 'longwave net flux', datatype_flt)
    grib2_desc = grib2_var(0, 5, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflxall', diag%lwflxall,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
      & ldims=shape3dkp1, lopenacc=.TRUE. )
    

    ! &      diag%trsolall(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('trsolall', '', 'shortwave net tranmissivity', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'trsolall', diag%trsolall,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
      & ldims=shape3dkp1, lopenacc=.TRUE.                                 )
    
      
    ! Set dimensions for 3D radiative flux variables
    IF (atm_phy_nwp_config(k_jg)%l_3d_rad_fluxes) THEN
       shape3dflux = (/ nproma, klevp1, kblks /)
       lrestart_flux = .TRUE.
    ELSE
       shape3dflux = (/ 1, 1, kblks/)
       lrestart_flux = .FALSE.
    END IF

    ! &      diag%lwflx_up(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('lwflx_up', 'W m-2 ', 'longwave upward flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 201, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflx_up', diag%lwflx_up,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE. )
    

    ! &      diag%lwflx_dn(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('lwflx_dn', 'W m-2 ', 'longwave downward flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflx_dn', diag%lwflx_dn,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
    
  
    ! &      diag%swflx_up(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('swflx_up', 'W m-2 ', 'shortwave upward flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_up', diag%swflx_up,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
      
  
    ! &      diag%swflx_dn(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('swflx_dn', 'W m-2 ', 'shortwave downward flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 204, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_dn', diag%swflx_dn,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
    

    ! &      diag%lwflx_up_clr(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('lwflx_up_clr', 'W m-2 ', 'longwave upward clear-sky flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 205, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflx_up_clr', diag%lwflx_up_clr,           &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
    
 
    ! &      diag%lwflx_dn_clr(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('lwflx_dn_clr', 'W m-2 ', 'longwave downward clear-sky flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 206, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lwflx_dn_clr', diag%lwflx_dn_clr,           &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
    
  
    ! &      diag%swflx_up_clr(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('swflx_up_clr', 'W m-2 ', 'shortave upward clear-sky flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 207, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_up_clr', diag%swflx_up_clr,           &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
    
  
    ! &      diag%swflx_dn_clr(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('swflx_dn_clr', 'W m-2 ', 'shortave downward clear-sky flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 208, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'swflx_dn_clr', diag%swflx_dn_clr,           &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, ldims=shape3dflux, lrestart=lrestart_flux, &
      & lopenacc=.TRUE.)
    


    !------------------
    !Turbulence 2D variables
    
    ! &      diag%shfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('shfl_s', 'W m-2 ', 'surface sensible heat flux', datatype_flt)
    grib2_desc = grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'shfl_s', diag%shfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    

    WRITE(name,'(A,A6)') TRIM(prefix),"shfl_s"
    WRITE(long_name,'(A27,A4,A18)') "surface sensible heat flux ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%ashfl_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("pbl_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE.                                                   )
      


    ! &      diag%lhfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('lhfl_s', 'W m-2 ', 'surface latent heat flux', datatype_flt)
    grib2_desc = grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_s', diag%lhfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE.)
    
                
    WRITE(name,'(A,A6)') TRIM(prefix),"lhfl_s"
    WRITE(long_name,'(A27,A4,A18)') "surface latent heat flux ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%alhfl_s,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("pbl_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ), lopenacc=.TRUE. )
    


    ! &      diag%lhfl_bs(nproma,nblks_c)
    cf_desc    = t_cf_var('lhfl_bs', 'W m-2 ', 'latent heat flux from bare soil', &
      &          datatype_flt)
    grib2_desc = grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_bs', diag%lhfl_bs,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE.)
    

    WRITE(name,'(A,A7)') TRIM(prefix),"lhfl_bs"
    WRITE(long_name,'(A27,A4,A18)') "latent heat flux from bare soil", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%alhfl_bs,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & isteptype=a_steptype,                                             &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE.                                                   )
    

    ! &      diag%lhfl_pl(nproma,nlev_soil,nblks_c)
    cf_desc    = t_cf_var('lhfl_pl', 'W m-2 ', 'latent heat flux from plants', &
      &          datatype_flt)
    grib2_desc = grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_pl', diag%lhfl_pl,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
      & ldims=(/nproma,nlev_soil,kblks/), lrestart=.FALSE.,               &
      & lopenacc=.TRUE.)
    
              
    WRITE(name,'(A,A7)') TRIM(prefix),"lhfl_pl"
    WRITE(long_name,'(A27,A4,A18)') "latent heat flux from plants", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%alhfl_pl,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
      & ldims=(/nproma,nlev_soil,kblks/), lrestart=.TRUE.,                &
      & isteptype=a_steptype,                                             &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE.                                                   )
    

    ! &      diag%qhfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('qhfl_s', 'Kg m-2 s-1', 'surface moisture flux', datatype_flt)
    grib2_desc = grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qhfl_s', diag%qhfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.TRUE., lopenacc=.TRUE.,                  &
      & in_group=groups("pbl_vars"))
    

    WRITE(name,'(A,A6)') TRIM(prefix),"qhfl_s"
    WRITE(long_name,'(A23,A4,A18)') "surface moisture flux ", meaning, &
                                  & " since model start"
    IF (lflux_avg ) THEN
      cf_desc = t_cf_var(name, 'Kg m-2 s-1', long_name, datatype_flt)
    ELSE
      cf_desc = t_cf_var(name, 'Kg m-2', long_name, datatype_flt)
    ENDIF
    grib2_desc = grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%aqhfl_s  ,                  &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & isteptype=a_steptype, in_group=groups("pbl_vars"),                &
      & hor_interp=create_hor_interp_metadata(                            &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                          &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                         &
      & lopenacc=.TRUE.                                                   )
      


    cf_desc    = t_cf_var('qcfl_s', 'kg m-2 s-1',                         &
      &          'surface cloud water deposition flux due to diffusion',  &
      &          datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qcfl_s', diag%qcfl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lopenacc=.TRUE.)
    

    cf_desc    = t_cf_var('qifl_s', 'kg m-2 s-1',                         &
      &          'surface cloud ice deposition flux due to diffusion',    &
      &          datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qifl_s', diag%qifl_s,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lopenacc=.TRUE.)
    

    ! Effective radius fields
    ! These are not additive fields and therefore nearest neig. interoplation is used for horizontal.
    ! Linear interpolation is used for vertical because nearest neig is not yet implemented

    IF (atm_phy_nwp_config(k_jg)%icalc_reff > 0) THEN 

      cf_desc      = t_cf_var('reff_qc', 'm',  'effective radius of cloud water', datatype_flt)
      grib2_desc  = grib2_var(0, 1, 129, ibits, GRID_UNSTRUCTURED, GRID_CELL)    ! Corresponds to RECLOUD
      CALL add_var( diag_list, 'reff_qc', diag%reff_qc,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
        & ldims=shape3d, lrestart=.TRUE.,                                       &
        & initval=1.0e-5_wp,                                                    & 
        & vert_interp=create_vert_interp_metadata(                              &
        &             vert_intp_type=vintp_types("P","Z","I"),                  &
        &             vert_intp_method=VINTP_METHOD_LIN,                        &
        &             l_loglin=.FALSE.,                                         &
        &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
        &             lower_limit=0._wp ),                                      &
        & hor_interp=create_hor_interp_metadata(                                &
        &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
        & lopenacc=.TRUE.                                                       )
       


      cf_desc      = t_cf_var('reff_qi', 'm',  'effective radius of cloud ice', datatype_flt)
      grib2_desc   = grib2_var(0, 1, 131, ibits, GRID_UNSTRUCTURED, GRID_CELL)    ! Corresponds to REICE
      CALL add_var( diag_list, 'reff_qi', diag%reff_qi,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
        & ldims=shape3d, lrestart=.TRUE.,                                       &
        & initval=1.5e-5_wp,                                                    &
        & vert_interp=create_vert_interp_metadata(                              &
        &             vert_intp_type=vintp_types("P","Z","I"),                  &
        &             vert_intp_method=VINTP_METHOD_LIN,                        &
        &             l_loglin=.FALSE.,                                         &
        &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
        &             lower_limit=0._wp ),                                      &
        & hor_interp=create_hor_interp_metadata(                                &
        &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
        & lopenacc=.TRUE.                                                       )
       


      cf_desc      = t_cf_var('reff_qr', 'm',  'effective radius of rain droplets', datatype_flt)
      grib2_desc   = grib2_var(0, 1, 130, ibits, GRID_UNSTRUCTURED, GRID_CELL)    ! Corresponds to RERAIN
      CALL add_var( diag_list, 'reff_qr', diag%reff_qr,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
        & ldims=shape3d, lrestart=.TRUE.,                                       &
        & initval=5.0e-4_wp,                                                    &
        & vert_interp=create_vert_interp_metadata(                              &
        &             vert_intp_type=vintp_types("P","Z","I"),                  &
        &             vert_intp_method=VINTP_METHOD_LIN,                        &
        &             l_loglin=.FALSE.,                                         &
        &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
        &             lower_limit=0._wp ),                                      &
        & hor_interp=create_hor_interp_metadata(                                &
        &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
        & lopenacc=.TRUE.                                                       )
      

      IF (iqs > 0) THEN
        cf_desc      = t_cf_var('reff_qs', 'm',  'effective radius of snow', datatype_flt)
        grib2_desc   = grib2_var(0, 1, 132, ibits, GRID_UNSTRUCTURED, GRID_CELL)    ! Corresponds to RESNOW
        CALL add_var( diag_list, 'reff_qs', diag%reff_qs,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
          & ldims=shape3d, lrestart=.TRUE.,                                       &
          & initval=1.0e-4_wp,                                                    &
          & vert_interp=create_vert_interp_metadata(                              &
          &             vert_intp_type=vintp_types("P","Z","I"),                  &
          &             vert_intp_method=VINTP_METHOD_LIN,                        &
          &             l_loglin=.FALSE.,                                         &
          &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
          &             lower_limit=0._wp ),                                      &
          & hor_interp=create_hor_interp_metadata(                                &
          &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
          & lopenacc=.TRUE.                                                       )
        

      ENDIF

      IF (iqg > 0) THEN
        cf_desc      = t_cf_var('reff_qg', 'm',  'effective radius of graupel', datatype_flt)
        grib2_desc   = grib2_var(0, 1, 133, ibits, GRID_UNSTRUCTURED, GRID_CELL)    ! Corresponds to REGRAUPEL
        CALL add_var( diag_list, 'reff_qg', diag%reff_qg,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
          & ldims=shape3d, lrestart=.TRUE.,                                       &
          & initval=3.0e-4_wp,                                                    &
          & vert_interp=create_vert_interp_metadata(                              &
          &             vert_intp_type=vintp_types("P","Z","I"),                  &
          &             vert_intp_method=VINTP_METHOD_LIN,                        &
          &             l_loglin=.FALSE.,                                         &
          &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
          &             lower_limit=0._wp ),                                      &
          & hor_interp=create_hor_interp_metadata(                                &
          &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
          & lopenacc=.TRUE.                                                       )
        
      END IF

      IF (iqh > 0) THEN
        cf_desc      = t_cf_var('reff_qh', 'm',  'effective radius of hail', datatype_flt)
        grib2_desc   = grib2_var(0, 1, 134, ibits, GRID_UNSTRUCTURED, GRID_CELL)    ! Corresponds to REHAIL
        CALL add_var( diag_list, 'reff_qh', diag%reff_qh,                         &
          & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,            &
          & ldims=shape3d, lrestart=.TRUE.,                                       &
          & initval=1.0e-3_wp,                                                    &
          & vert_interp=create_vert_interp_metadata(                              &
          &             vert_intp_type=vintp_types("P","Z","I"),                  &
          &             vert_intp_method=VINTP_METHOD_LIN,                        &
          &             l_loglin=.FALSE.,                                         &
          &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,                   &
          &             lower_limit=0._wp ),                                      &
          & hor_interp=create_hor_interp_metadata(                                &
          &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
          & lopenacc=.TRUE.                                                       )
        
      END IF
    END IF  ! IF effective radius


    ! &      diag%tcm(nproma,nblks_c)
    cf_desc    = t_cf_var('tcm', ' ','turbulent transfer coefficients for momentum', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tcm', diag%tcm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    

    ! &      diag%tch(nproma,nblks_c)
    cf_desc    = t_cf_var('tch', ' ','turbulent transfer coefficients for heat', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tch', diag%tch,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    
    
    ! &      diag%tfm(nproma,nblks_c)
    cf_desc    = t_cf_var('tfm', ' ','factor of laminar transfer of momentum', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfm', diag%tfm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lopenacc=.TRUE. )
    

    ! &      diag%tfh(nproma,nblks_c)
    cf_desc    = t_cf_var('tfh', ' ','factor of laminar transfer of scalars', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfh', diag%tfh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lopenacc=.TRUE. )
    

    ! &      diag%tfv(nproma,nblks_c)
    cf_desc    = t_cf_var('tfv', ' ','laminar reduction factor for evaporation', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfv', diag%tfv,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lopenacc=.TRUE. )
    

    ! &      diag%tvm(nproma,nblks_c)
    cf_desc    = t_cf_var('tvm', 'm s-1','turbulent transfer velocity for momentum', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvm', diag%tvm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    

    ! &      diag%tvh(nproma,nblks_c)
    cf_desc    = t_cf_var('tvh', 'm s-1 ','turbulent transfer velocity for heat', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvh', diag%tvh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    

    ! &      diag%tkr(nproma,nblks_c)
    cf_desc    = t_cf_var('tkr', 'm2 s-1 ','turbulent reference surface diffusion coefficient', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkr', diag%tkr,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    

    ! &      diag%tkred_sfc(nproma,nblks_c)
    cf_desc    = t_cf_var('tkred_sfc', ' ','reduction factor for minimum diffusion coefficients for EPS perturbations', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkred_sfc', diag%tkred_sfc,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE., initval=1._wp, lopenacc=.TRUE. )
    

    ! &      diag%tkred_sfc_h(nproma,nblks_c)
    cf_desc    = t_cf_var('tkred_sfc_h', ' ','reduction factor for minimum diffusion coefficient for heat for model-DA coupling', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkred_sfc_h', diag%tkred_sfc_h,             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.TRUE., initval=1._wp, lopenacc=.TRUE. )
    

    ! &      diag%pat_len(nproma,nblks_c)
    cf_desc    = t_cf_var('pat_len', 'm','length scale of sub-grid scale roughness elements', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'pat_len', diag%pat_len,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. )
    

    !        diag%rlamh_fac_t (nproma, nblks, ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('rlamh_fac_t', '', 'scaling factor for rlam_heat', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rlamh_fac_t', diag%rlamh_fac_t,          &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,       &
      & ldims=shape3dsubsw, lcontainer=.TRUE., lrestart=.FALSE.,       &
      & loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the seperate variables belonging to the container rlamh_fac_t
    ALLOCATE(diag%rlamh_fac_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc 
      CALL add_ref( diag_list, 'rlamh_fac_t',                            &
         & 'rlamh_fac_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%rlamh_fac_ptr(jsfc)%p_2d,                                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
         & t_cf_var('rlamh_fac_'//TRIM(csfc), '', '', datatype_flt),     &
         & grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                                  &
         & var_class=CLASS_TILE,                                         &
         & lrestart=.FALSE.)
    ENDDO



    ! &      diag%gz0(nproma,nblks_c)
    cf_desc     = t_cf_var('gz0', 'm2 s-2 ','roughness length times gravity', datatype_flt)
    new_cf_desc = t_cf_var( 'z0',       'm','roughness length',               datatype_flt)
    grib2_desc = grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'gz0', diag%gz0,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
      & ldims=shape2d,                                                    &
      & post_op=post_op(POST_OP_SCALE, arg1=1._wp/grav,                   &
      &                 new_cf=new_cf_desc),                              &
      & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in",               &
      &                 "mode_iau_fg_in","mode_iau_old_fg_in",            &
      &                 "mode_iniana","iau_restore_vars"),                &
      & initval=0.01_wp, lopenacc=.TRUE. )
    


    IF (is_coupled_to_waves()) THEN
      ! &      diag%z0_waves(nproma,nblks_c)
      cf_desc     = t_cf_var('z0_waves', 'm','wave-dependent roughness length', datatype_flt)
      grib2_desc = grib2_var(10, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'z0_waves', diag%z0_waves,                   &
        & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
        & ldims=shape2d, loutput=.TRUE.,                                    &
        & lopenacc=.FALSE. )
    END IF


    ! &      diag%t_2m(nproma,nblks_c)
    IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
      in_group = groups("pbl_vars","dwd_fg_atm_vars","mode_iau_ana_in")
    ELSE
      in_group = groups("pbl_vars","dwd_fg_atm_vars")
    ENDIF
    cf_desc    = t_cf_var('t_2m', 'K ','temperature in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 't_2m', diag%t_2m,                           &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & in_group=in_group,                                                &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE. )
    


    ! &      diag%t_tilemax_inst_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('t_tilemax_inst_2m', 'K ', &
      &                   'instantaneous temperature in 2m, maximum over tiles', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 't_tilemax_inst_2m', diag%t_tilemax_inst_2m, &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE. )
    


    ! &      diag%t_tilemin_inst_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('t_tilemin_inst_2m', 'K ', &
      &                   'instantaneous temperature in 2m, minimum over tiles', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 't_tilemin_inst_2m', diag%t_tilemin_inst_2m, &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE.,                                  &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE. )
    



    ! &      diag%t_2m_land(nproma,nblks_c)
    cf_desc    = t_cf_var('t_2m_land', 'K ','temperature in 2m over land fraction', &
      &          datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 181)
    CALL add_var( diag_list, 't_2m_land', diag%t_2m_land,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M_LAYER, cf_desc, grib2_desc,  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars"),     &
      & lopenacc=.TRUE. )
    

    ! &      diag%tmax_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('tmax_2m', 'K ','Max 2m temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tmax_2m', diag%tmax_2m,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & isteptype=TSTEP_MAX, initval=-999._wp, resetval=-999._wp,         &
      & action_list=actions(new_action(ACTION_RESET,maxt_interval(k_jg))),&
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE. ) 
    

    ! &      diag%tmin_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('tmin_2m', 'K ','Min 2m temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tmin_2m', diag%tmin_2m,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.TRUE.,                                   &
      & isteptype=TSTEP_MIN, initval=999._wp, resetval=999._wp,           &
      & action_list=actions(new_action(ACTION_RESET,maxt_interval(k_jg))),&
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE. )
    

    ! &      diag%qv_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('qv_2m', 'kg kg-1 ','specific water vapor content in 2m', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qv_2m', diag%qv_2m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars"),     &
      & lopenacc=.TRUE. )
    

    ! &      diag%rh_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('rh_2m', '%','relative humidity in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rh_2m', diag%rh_2m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. )
    

    ! &      diag%rh_2m_land(nproma,nblks_c)
    cf_desc    = t_cf_var('rh_2m_land', '%','relative humidity in 2m over land fraction', &
      &          datatype_flt)
    grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 181)
    CALL add_var( diag_list, 'rh_2m_land', diag%rh_2m_land,               &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M_LAYER, cf_desc, grib2_desc,  &
      & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE. )
    

    ! &      diag%td_2m(nproma,nblks_c)
    cf_desc    = t_cf_var('td_2m', 'K ','dew-point in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'td_2m', diag%td_2m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars","dwd_fg_atm_vars"), &
      & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
      &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
      & lopenacc=.TRUE. )
    

    ! &      diag%td_2m_land(nproma,nblks_c)
    cf_desc    = t_cf_var('td_2m_land', 'K ','dew-point in 2m over land fraction', &
      &           datatype_flt)
    grib2_desc = grib2_var(0, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 181)
    CALL add_var( diag_list, 'td_2m_land', diag%td_2m_land,               &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M_LAYER, cf_desc, grib2_desc,  &
      & ldims=shape2d, lrestart=.FALSE., in_group=groups("pbl_vars"),     &
      & lopenacc=.TRUE. )
    

    ! &      diag%u_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('u_10m', 'm s-1 ','zonal wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'u_10m', diag%u_10m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
      & ldims=shape2d, lrestart=lart, in_group=groups("pbl_vars","dwd_fg_atm_vars"), &
      & lopenacc=.TRUE. )
    

    ! &      diag%v_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('v_10m', 'm s-1 ','meridional wind in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'v_10m', diag%v_10m,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
      & ldims=shape2d, lrestart=lart, in_group=groups("pbl_vars","dwd_fg_atm_vars"), &
      & lopenacc=.TRUE. )
    

    ! &      diag%sp_10m(nproma,nblks_c)
    cf_desc    = t_cf_var('sp_10m', 'm s-1 ','wind speed in 10m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sp_10m', diag%sp_10m,                       &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,       &
      & ldims=shape2d, lrestart=.TRUE., lopenacc=.TRUE. )
    

    !tiled quantities
    ! &      diag%shfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('shfl_s_t', 'W m-2 ', 'tile-based surface sensible heat flux', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'shfl_s_t', diag%shfl_s_t,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                       &
      & in_group=groups("iau_init_vars"),                                           &
      & lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container shfl_s_t
    ALLOCATE(diag%shfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'shfl_s_t',                            &
         & 'shfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%shfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('shfl_s_t_'//csfc, '', '', datatype_flt), &
         & grib2_var(0, 0, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%lhfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('lhfl_s_t', 'W m-2 ', 'tile-based surface latent heat flux', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_s_t', diag%lhfl_s_t,                                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,   &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container lhfl_s_t
    ALLOCATE(diag%lhfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'lhfl_s_t',                            &
         & 'lhfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%lhfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('lhfl_s_t_'//csfc, '', '', datatype_flt),         &
         & grib2_var(0, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%lhfl_bs_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('lhfl_bs_t', 'W m-2 ', 'tile-based latent heat flux from bare soil', &
      &                   datatype_flt)
    grib2_desc = grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_bs_t', diag%lhfl_bs_t,                              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubs,    &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container lhfl_bs_t
    ALLOCATE(diag%lhfl_bs_t_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'lhfl_bs_t',                           &
         & 'lhfl_bs_t_'//TRIM(ADJUSTL(csfc)),                         &
         & diag%lhfl_bs_t_ptr(jsfc)%p_2d,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('lhfl_bs_t_'//csfc, '', '', datatype_flt),        &
         & grib2_var(2, 0, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE_LAND,                                 &
         & lrestart=.FALSE., loutput=.TRUE.)
    ENDDO

    ! &      diag%lhfl_pl_t(nproma,nlev_soil,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('lhfl_pl_t', 'W m-2 ', 'tile-based latent heat flux from plants', &
      &                   datatype_flt)
    grib2_desc = grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'lhfl_pl_t', diag%lhfl_pl_t,                 &
      & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc, &
      & ldims=(/nproma,nlev_soil,kblks,ntiles_total/), lcontainer=.TRUE., &
      & lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container lhfl_pl_t
    ALLOCATE(diag%lhfl_pl_t_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'lhfl_pl_t',                           &
         & 'lhfl_pl_t_'//TRIM(ADJUSTL(csfc)),                         &
         & diag%lhfl_pl_t_ptr(jsfc)%p_3d,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,               &
         & t_cf_var('lhfl_pl_t_'//csfc, '', '', datatype_flt),        &
         & grib2_var(2, 0, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
         & ref_idx=jsfc,                                              &
         & var_class=CLASS_TILE_LAND,                                 &
         & ldims=(/nproma,nlev_soil,kblks/), lrestart=.FALSE.,        &
         & loutput=.TRUE.)
    ENDDO

    ! &      diag%qhfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('qhfl_s_t', 'Kg m-2 s-1','tile based surface moisture flux', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'qhfl_s_t', diag%qhfl_s_t,                                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,   &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                          &
      & in_group=groups("iau_init_vars"),                                              &
      & lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container qhfl_s_t
    ALLOCATE(diag%qhfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'qhfl_s_t',                            &
         & 'qhfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%qhfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('qhfl_s_t_'//csfc, '', '', datatype_flt),         &
         & grib2_var(2, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tcm_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tcm_t', ' ', &
      & 'tile-based turbulent transfer coefficients for momentum', datatype_flt)
    grib2_desc = grib2_var(0, 2, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tcm_t', diag%tcm_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container tcm_t
    ALLOCATE(diag%tcm_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tcm_t',                               &
         & 'tcm_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tcm_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tcm_t_'//csfc, '', '', datatype_flt),            &
         & grib2_var(0, 2, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tch_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tch_t', ' ', &
         &                'tile-based turbulent transfer coefficients for heat', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tch_t', diag%tch_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)                       
    

    ! fill the separate variables belonging to the container tch_t
    ALLOCATE(diag%tch_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tch_t',                               &
         & 'tch_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tch_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tch_t_'//csfc, '', '', datatype_flt),            &
         & grib2_var(0, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tfv_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tfv_t', ' ', &
         &                'tile-based laminar reduction factor for evaporation', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tfv_t', diag%tfv_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container tfv_t
    ALLOCATE(diag%tfv_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tfv_t',                               &
         & 'tfv_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tfv_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tfv_t_'//csfc, '', '', datatype_flt),            &
         & grib2_var(0, 4, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tvm_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tvm_t', 'm s-1', &
      & 'tile-based turbulent transfer velocity for momentum', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvm_t', diag%tvm_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container tvm_t
    ALLOCATE(diag%tvm_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tvm_t',                               &
         & 'tvm_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tvm_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tvm_t_'//csfc, '', '', DATATYPE_FLT32),          &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
!        & lrestart=.TRUE., loutput=.TRUE.)  !nec. for rest. only if 'tvm' substitutes 'tcm'
         & lrestart=.FALSE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tvh_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tvh_t', 'm s-1', &
         &                'tile-based turbulent transfer velocity for heat', &
         &                DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvh_t', diag%tvh_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container tvh_t
    ALLOCATE(diag%tvh_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tvh_t',                               &
         & 'tvh_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tvh_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tvh_t_'//csfc, '', '', DATATYPE_FLT32),          &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
!        & lrestart=.TRUE., loutput=.TRUE.) !nec. for rest. only if 'tvm' substitutes 'tcm'
         & lrestart=.FALSE., loutput=.TRUE.)
    ENDDO

    ! &      diag%tkr_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tkr_t', 'm2 s-1', &
      & 'tile-based turbulent reference surface diffusion coefficient', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkr_t', diag%tkr_t,                                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    

    ! fill the separate variables belonging to the container tkr_t
    ALLOCATE(diag%tkr_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tkr_t',                               &
         & 'tkr_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%tkr_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('tkr_t_'//csfc, '', '', DATATYPE_FLT32),    &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
!        & lrestart=.TRUE., loutput=.TRUE.) !nec. for rest. only if 'imode_trancnf>=4'
         & lrestart=.FALSE., loutput=.TRUE.)
    ENDDO


    ! &      diag%gz0_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('gz0_t', 'm2 s-2 ', 'tile-based roughness length times gravity', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'gz0_t', diag%gz0_t,                                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                        &
      & in_group=groups("iau_restore_vars"),                                         &
      & lopenacc=.TRUE.)
    



    ! fill the separate variables belonging to the container gz0_t
    ALLOCATE(diag%gz0_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'gz0_t',                               &
         & 'gz0_t_'//TRIM(ADJUSTL(csfc)),                             &
         & diag%gz0_t_ptr(jsfc)%p_2d,                                 &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('gz0_t_'//csfc, '', '', datatype_flt),            &
         & grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO


    ! &      diag%tvs_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tvs_s_t', 'm2 s-2 ',                              &
         &                'tile-based turbulence velocity scale at surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tvs_s_t', diag%tvs_s_t,                                &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    


    ! fill the separate variables belonging to the container tvs_s_t
    ALLOCATE(diag%tvs_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tvs_s_t',                                &
         & 'tvs_s_t_'//TRIM(ADJUSTL(csfc)),                              &
         & diag%tvs_s_t_ptr(jsfc)%p_2d,                                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
         & t_cf_var('tvs_s_t_'//csfc, '', '', datatype_flt),             &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
         & ref_idx=jsfc, ldims=shape2d,                                  &
         & var_class=CLASS_TILE,                                         &
         & lrestart=.TRUE., loutput=.FALSE.)
    ENDDO



    ! &      diag%tkvm_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tkvm_s_t', 'm s-2 ',                                       &
         &                'tile-based turbulent diff. coeff for momentum at surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvm_s_t', diag%tkvm_s_t,                              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    


    ! fill the separate variables belonging to the container tkvm_s_t
    ALLOCATE(diag%tkvm_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tkvm_s_t',                                &
         & 'tkvm_s_t_'//TRIM(ADJUSTL(csfc)),                              &
         & diag%tkvm_s_t_ptr(jsfc)%p_2d,                                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
         & t_cf_var('tkvm_s_t_'//csfc, '', '', datatype_flt),             &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
         & ref_idx=jsfc, ldims=shape2d,                                   &
         & var_class=CLASS_TILE,                                          &
         & lrestart=.TRUE., loutput=.FALSE.)
    ENDDO


    ! &      diag%tkvh_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('tkvh_s_t', 'm s-2 ',                                      &
         &                'tile-based turbulent diff. coeff for heat at surface',    &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvh_s_t', diag%tkvh_s_t,                              &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw, &
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    


    ! fill the separate variables belonging to the container tkvm_s_t
    ALLOCATE(diag%tkvh_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'tkvh_s_t',                                &
         & 'tkvh_s_t_'//TRIM(ADJUSTL(csfc)),                              &
         & diag%tkvh_s_t_ptr(jsfc)%p_2d,                                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
         & t_cf_var('tkvh_s_t_'//csfc, '', '', datatype_flt),             &
         & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
         & ref_idx=jsfc, ldims=shape2d,                                   &
         & var_class=CLASS_TILE,                                          &
         & lrestart=.TRUE., loutput=.FALSE.)
    ENDDO



    ! &      diag%u_10m_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('u_10m_t', 'm s-1 ', 'tile-based zonal wind in 2m', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'u_10m_t', diag%u_10m_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    
  
    ! fill the separate variables belonging to the container u_10m_t
    ALLOCATE(diag%u_10m_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'u_10m_t',                             &
         & 'u_10m_t_'//TRIM(ADJUSTL(csfc)),                           &
         & diag%u_10m_t_ptr(jsfc)%p_2d,                               &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                     &
         & t_cf_var('u_10m_t_'//csfc, '', '', datatype_flt),          &
         & grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO


    ! &      diag%v_10m_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('v_10m_t', 'm s-1 ', 'tile-based meridional wind in 2m', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'v_10m_t', diag%v_10m_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    


    ! fill the separate variables belonging to the container v_10m_t
    ALLOCATE(diag%v_10m_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'v_10m_t',                             &
         & 'v_10m_t_'//TRIM(ADJUSTL(csfc)),                           &
         & diag%v_10m_t_ptr(jsfc)%p_2d,                               &
         & GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M,                     &
         & t_cf_var('v_10m_t_'//csfc, '', '', datatype_flt),          &
         & grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.)
    ENDDO



    ! &      diag%umfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('umfl_s_t', 'N m-2 ', 'u-momentum flux at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'umfl_s_t', diag%umfl_s_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    
  
    ! fill the separate variables belonging to the container umfl_s_t
    ALLOCATE(diag%umfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'umfl_s_t',                            &
         & 'umfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%umfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('umfl_s_t_'//csfc, '', '', datatype_flt),         &
         & grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.,                           &
         & isteptype=TSTEP_INSTANT )
    ENDDO


    ! &      diag%vmfl_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('vmfl_s_t', 'N m-2 ', 'v-momentum flux at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'vmfl_s_t', diag%vmfl_s_t,                                  &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape3dsubsw,&
      & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
    
  
    ! fill the separate variables belonging to the container vmfl_s_t
    ALLOCATE(diag%vmfl_s_t_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      WRITE(csfc,'(i1)') jsfc
      CALL add_ref( diag_list, 'vmfl_s_t',                            &
         & 'vmfl_s_t_'//TRIM(ADJUSTL(csfc)),                          &
         & diag%vmfl_s_t_ptr(jsfc)%p_2d,                              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
         & t_cf_var('vmfl_s_t_'//csfc, '', '', datatype_flt),         &
         & grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
         & ref_idx=jsfc, ldims=shape2d,                               &
         & var_class=CLASS_TILE,                                      &
         & lrestart=.TRUE., loutput=.TRUE.,                           &
         & isteptype=TSTEP_INSTANT )
    ENDDO


    ! &      diag%umfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('umfl_s', 'N m-2', 'u-momentum flux at the surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'umfl_s', diag%umfl_s,                            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & lrestart=.FALSE., loutput=.TRUE.,                                      &
      & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    

    ! &      diag%vmfl_s(nproma,nblks_c)
    cf_desc    = t_cf_var('vmfl_s', 'N m-2', 'v-momentum flux at the surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'vmfl_s', diag%vmfl_s,                            &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & lrestart=.FALSE., loutput=.TRUE.,                                      &
      & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    

    IF (lflux_avg ) THEN
        prefix = "a"
        meaning = "mean"
        varunits= "N/m**2"
        a_steptype= TSTEP_AVG 
    ELSE
        prefix = "acc"
        meaning = "acc." 
        varunits= "Ns/m**2"    ! or "kg/(m*s)"
        a_steptype= TSTEP_ACCUM     
    END IF

    WRITE(name,'(A,A6)') TRIM(prefix),"umfl_s"
    WRITE(long_name,'(A26,A4,A18)') "u-momentum flux flux at surface ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 2, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
!   aumfl_s and avmfl_s are needed for the restart only to get the same output values
!   They are not important to obtain bit identical model results
    CALL add_var( diag_list, TRIM(name), diag%aumfl_s,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                              &
      & lopenacc=.TRUE.                                                        )
      

    WRITE(name,'(A,A6)') TRIM(prefix),"vmfl_s"
    WRITE(long_name,'(A26,A4,A18)') "v-momentum flux flux at surface ", meaning, &
                                  & " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(0, 2, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%avmfl_s,                         &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                              &
      & lopenacc=.TRUE.                                                        )
      


    !------------------------------
    ! SSO surface stress
    !------------------------------

    cf_desc    = t_cf_var('str_u_sso', 'N m-2', 'zonal sso surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 128, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'str_u_sso', diag%str_u_sso,                     &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,  &
      & lopenacc=.TRUE.)
    

    cf_desc    = t_cf_var('str_v_sso', 'N m-2', 'meridional sso surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 128, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'str_v_sso', diag%str_v_sso,                        &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d , &
      & lopenacc=.TRUE.)
    

    WRITE(name,'(A,A9)') TRIM(prefix),"str_u_sso"
    WRITE(long_name,'(A25,A4,A18)') "zonal sso surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(192, 128, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%astr_u_sso,                      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                              &
      & lopenacc=.TRUE.                                                        )
    

    WRITE(name,'(A,A9)') TRIM(prefix),"str_v_sso"
    WRITE(long_name,'(A30,A4,A18)') "meridional sso surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(192, 128, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%astr_v_sso,                      &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                              &
      & lopenacc=.TRUE.                                                        )
    


    !------------------------------
    ! grid-scale surface stress
    !------------------------------

    cf_desc    = t_cf_var('drag_u_grid', 'N m-2', 'zonal resolved surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 150, 168, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'drag_u_grid', diag%drag_u_grid,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d, &
      & lopenacc=.TRUE.                                                         )
    

    cf_desc    = t_cf_var('drag_v_grid', 'N m-2', 'meridional resolved surface stress', datatype_flt)
    grib2_desc = grib2_var(192, 150, 169, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'drag_v_grid', diag%drag_v_grid,                   &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d, &
      & lopenacc=.TRUE.                                                         )
    

    WRITE(name,'(A,A11)') TRIM(prefix),"drag_u_grid"
    WRITE(long_name,'(A30,A4,A18)') "zonal resolved surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(192, 150, 168, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%adrag_u_grid,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                              &
      & lopenacc=.TRUE.                                                        )
    

    WRITE(name,'(A,A11)') TRIM(prefix),"drag_v_grid"
    WRITE(long_name,'(A35,A4,A18)') "meridional resolved surface stress ", meaning, " since model start"
    cf_desc    = t_cf_var(name, varunits, long_name, datatype_flt)
    grib2_desc = grib2_var(192, 150, 169, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, TRIM(name), diag%adrag_v_grid,                    &
      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,&
      & isteptype=a_steptype, lrestart=.TRUE., loutput=.TRUE.,                 &
      & hor_interp=create_hor_interp_metadata(                                 &
      &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                               &
      &    fallback_type=HINTP_TYPE_LONLAT_RBF ),                              &
      & lopenacc=.TRUE.                                                        )
    


    !------------------------------
    ! LES
    !------------------------------

    !Anurag Dipankar, MPI (7 Oct 2013)
    !Diagnostics for LES physics

    IF ( atm_phy_nwp_config(k_jg)%is_les_phy ) THEN

      ! &      diag%z_pbl(nproma,nblks_c)
      cf_desc    = t_cf_var('z_pbl', 'm', 'boundary layer height above sea level', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'z_pbl', diag%z_pbl,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%bruvais(nproma,nlev+1,nblks_c)
      cf_desc    = t_cf_var('bruvais', '1/s**2', 'Brunt Vaisala Frequency', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'bruvais', diag%bruvais,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
        & ldims=shape3dkp1, initval=0.0_wp, lrestart=.FALSE.,               &
        & lopenacc=.TRUE.                                                   )
      

      ! &      diag%mech_prod(nproma,nlev+1,nblks_c)
      cf_desc    = t_cf_var('mech_prod', 'm**2/s**3', 'mechanical production term in TKE Eq', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mech_prod', diag%mech_prod,                     &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,       &
        & ldims=shape3dkp1, initval=0.0_wp, lrestart=.FALSE.,                   &
        & lopenacc=.TRUE.                                                       )
      

      !1D and 0D diagnostic variables that can not be part of add_var
      ALLOCATE( diag%turb_diag_1dvar(klevp1,SIZE(turb_profile_list,1)),  &
                diag%turb_diag_0dvar(SIZE(turb_tseries_list,1)), STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:new_nwp_phy_diag_list', &
                    'allocation of 1D and 0D diag var failed!')
      ENDIF
      diag%turb_diag_1dvar = 0._wp
      diag%turb_diag_0dvar = 0._wp

      !  
      !Some diagnostics specific to HDCP2
      !

      ! &      diag%t_cbase(nproma,nblks_c) 
      cf_desc    = t_cf_var('t_cbase', 'K', 'cloud base temperature', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 't_cbase', diag%t_cbase,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%p_cbase(nproma,nblks_c) 
      cf_desc    = t_cf_var('p_cbase', 'Pa', 'cloud base pressure', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'p_cbase', diag%p_cbase,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%t_ctop(nproma,nblks_c) 
      cf_desc    = t_cf_var('t_ctop', 'K', 'cloud top temperature', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 't_ctop', diag%t_ctop,                         &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,          &
        & ldims=shape2d, lrestart=.FALSE. )

      ! &      diag%p_ctop(nproma,nblks_c) 
      cf_desc    = t_cf_var('p_ctop', 'K', 'cloud top pressure', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'p_ctop', diag%p_ctop,                           &
        & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP, cf_desc, grib2_desc,            &
        & ldims=shape2d, lrestart=.FALSE. )

    END IF

    !------------------
    !Turbulence 3D variables
    !------------------

    IF (.NOT.atm_phy_nwp_config(k_jg)%is_les_phy) THEN
    ! &      diag%tkvm(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('tkvm', 'm**2/s', ' turbulent diffusion coefficients for momentum', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, in_group=groups("pbl_vars"), lopenacc=.TRUE. )
    

   ! &      diag%tkvh(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('tkvh', 'm**2/s', ' turbulent diffusion coefficients for heat', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, in_group=groups("pbl_vars"), lopenacc=.TRUE.) 
    
    ELSE
     ! &      diag%tkvm(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('tkvm', 'kg/(ms)', ' mass weighted turbulent viscosity', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 2, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'tkvm', diag%tkvm,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
        & ldims=shape3dkp1, lrestart=.FALSE., in_group=groups("pbl_vars"),    &
        & lopenacc=.TRUE. )
      

     ! &      diag%tkvh(nproma,nlevp1,nblks_c)
      cf_desc    = t_cf_var('tkvh', 'kg/(ms)', ' mass weighted turbulent diffusivity', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 0, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'tkvh', diag%tkvh,                             &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
        & ldims=shape3dkp1, lrestart=.FALSE., in_group=groups("pbl_vars"),    &
        & lopenacc=.TRUE. )
       
    ENDIF

   ! &      diag%rcld(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('rcld', '', 'standard deviation of the saturation deficit', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'rcld', diag%rcld,                             &
      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
      & ldims=shape3dkp1, lopenacc=.TRUE. )
    

! SO FAR UNUSED
!!$   ! &      diag%edr(nproma,nlevp1,nblks_c)
!!$    cf_desc    = t_cf_var('edr', '', 'eddy dissipation rate', datatype_flt)
!!$    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
!!$    CALL add_var( diag_list, 'edr', diag%edr,                               &
!!$      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
!!$      & ldims=shape3dkp1, lrestart=.FALSE. ) 

! turbulent flux profiles
    IF (atm_phy_nwp_config(k_jg)%l_3d_turb_fluxes) THEN
       shape3Dturb = (/ nproma, klevp1, kblks /)

       !      diag%tetfl_turb(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('tetfl_turb', 'K m-2 s-1', 'vertical turbulent theta flux', &
          &                datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( diag_list, 'tetfl_turb', diag%tetfl_turb,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
         & ldims=shape3dturb, lrestart=.FALSE. ,lopenacc = .TRUE.)
       


       !       diag%vapfl_turb(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('vapfl_turb', 'kg m-2 s-1', 'vertical turbulent water vapour flux', &
            &                datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( diag_list, 'vapfl_turb', diag%vapfl_turb,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
         & ldims=shape3dturb, lrestart=.FALSE.,lopenacc = .TRUE.)
       


        !      diag%liqfl_turb(nproma,nlevp1,nblks_c)
       cf_desc    = t_cf_var('liqfl_turb', 'kg  m-2 s-1', 'vertical turbulent liquid water flux', &
            &                datatype_flt)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( diag_list, 'liqfl_turb', diag%liqfl_turb,                             &
         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,        &
         & ldims=shape3dturb, lrestart=.FALSE. ,lopenacc = .TRUE.)
       

    END IF

    !------------------
    ! Optional computation of diagnostic fields
    !------------------

    ! Note: These tasks are registered for the post-processing scheduler
    !       which takes care of the regular update:
    ! 
    ! &     relative humidity
    !    
    IF (var_in_output%rh) THEN
      cf_desc    = t_cf_var('rh', '%', 'relative humidity', datatype_flt)
      grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "rh", diag%rh,                                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3d,                                                 &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.TRUE.,             &
                    &             lower_limit=0._wp ),                               &
                    & l_pp_scheduler_task=TASK_COMPUTE_RH, lrestart=.FALSE.,         &
                    & lopenacc=.TRUE.)
      
    END IF
    
    ! &     potential vorticity
    
    IF (var_in_output%pv) THEN
      cf_desc    = t_cf_var('pv', 'K m2 kg-1 s-1', 'potential vorticity', DATATYPE_FLT32)
      grib2_desc = grib2_var(0, 2, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "pv", diag%pv,                                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3d,                                                 &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE.),                               &
                    & l_pp_scheduler_task=TASK_COMPUTE_PV, lrestart=.FALSE.,         &
                    & lopenacc=.TRUE.)
      
    END IF

    IF (var_in_output%sdi2) THEN
      cf_desc    = t_cf_var('sdi2', 's-1', 'supercell detection index (SDI2)', datatype_flt)
      grib2_desc = grib2_var(0, 7, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "sdi2", diag%sdi2,                                             &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_SDI2, lrestart=.FALSE., lopenacc=.TRUE.)
       
    END IF

    IF (var_in_output%dhail_mx .OR. var_in_output%dhail_av .OR. var_in_output%dhail_sd ) THEN

      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)

      ! Expected hailstone diameter of each sample size. Each sample is the maximum 
      ! over the time periode.
      ! This variable is used for hailcast computation
      ! grib2: no grib definition, not for nwp operational output
      cf_desc    = t_cf_var('dhail', 'mm', 'expected hailsize since end of previous full '//&
                             TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                     &
                    & "dhail", diag%dhail,                                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ldims=shape3d_hail,                                          &
                    & lrestart=.FALSE., loutput=.FALSE., isteptype=TSTEP_MAX,      &
                    & resetval=0.0_wp, initval=0.0_wp,                             &
                    & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                    & lopenacc = .TRUE. )
      

      ! Average expected hailstone diameter over all samples size. Each sample is the maximum over 
      ! the time periode.
      ! GRIB2: no grib definition, not for nwp operational output
      cf_desc    = t_cf_var('dhail_av', 'mm', 'average expected hailsize since end of previous full '//&
                             TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                     &
                    & "dhail_av", diag%dhail_av,                                   &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ldims=shape2d,                                               &
                    & lrestart=.FALSE., loutput=.TRUE., isteptype=TSTEP_MAX,       &
                    & resetval=0.0_wp, initval=0.0_wp,                             &
                    & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                    & lopenacc = .TRUE. )
      

      ! Maximum expected hailstone diameter over all samples size. Each sample is  
      ! the maximum over the time periode
      cf_desc    = t_cf_var('dhail_mx', 'mm', 'maximum expected hailsize since end of previous full ' //&
                             TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var(0, 1, 238, ibits, GRID_UNSTRUCTURED, GRID_CELL)       &
                    + t_grib2_int_key("typeOfStatisticalProcessing", 2)
      CALL add_var( diag_list,                                                     &
                    & "dhail_mx", diag%dhail_mx,                                   &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ldims=shape2d,                                               &
                    & lrestart=.FALSE., loutput=.TRUE., isteptype=TSTEP_MAX,       &
                    & resetval=0.0_wp, initval=0.0_wp,                             &
                    & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                    & lopenacc = .TRUE. )
      

      ! Standard deviation of the expected hailstone diameter over all samples size. 
      ! Each sample is the maximum over the time periode
      ! GRIB2: no grib definition, not for nwp operational output
      cf_desc    = t_cf_var('dhail_sd', 'mm', 'standard deviation of expected hailsize since end of previous full '//&
                          TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                     &
                    & "dhail_sd", diag%dhail_sd,                                   &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ldims=shape2d,                                               &
                    & lrestart=.FALSE., loutput=.TRUE., isteptype=TSTEP_MAX,       &
                    & resetval=0.0_wp, initval=0.0_wp,                             &
                    & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                    & lopenacc = .TRUE. )
      



    ! updraft duration diag%wdur(nproma,nblks_c)
    ! 
      cf_desc    = t_cf_var('wdur', 's', 'updraft duration', datatype_flt)
      grib2_desc = grib2_var(0, 2, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                   &
                    & "wdur", diag%wdur,                                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                        &
                    & cf_desc, grib2_desc,                                       &
                    & ldims=shape2d,                                             &
                    & lrestart=.TRUE.,                                           &
                    & loutput=.TRUE.,                                            &
                    & lopenacc = .TRUE. )
      


    ! mask of maximum updraft velocity diag%wup_mask(nproma,nblks_c)
    ! 
      cf_desc    = t_cf_var('wup_mask', '', 'updraft mask', DATATYPE_INT32)
      grib2_desc = grib2_var(0, 2, 211, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                     &
                    & "wup_mask", diag%wup_mask,                                   &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ldims=shape2d,                                               &
                    & lrestart=.TRUE.,                                             &
                    & loutput=.TRUE.,                                              &
                    & lopenacc = .TRUE. )
      

    END IF


    IF (var_in_output%lpi) THEN
      cf_desc    = t_cf_var('lpi', 'J kg-1', 'lightning potential index (LPI)', datatype_flt)
      grib2_desc = grib2_var(0, 17, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "lpi", diag%lpi,                                               &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_LPI, lrestart=.FALSE.,        &
                    & lopenacc=.TRUE.)
      
    END IF

    IF (var_in_output%lpi_max) THEN
      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)
      cf_desc    = t_cf_var('lpi_max', 'J kg-1',                   &
           &                 'lightning potential index, maximum since end of previous full '// &
                             TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 17, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'lpi_max', diag%lpi_max,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                  & resetval=0.0_wp, initval=0.0_wp,                         &
                  & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                  & lopenacc=.TRUE.)
      
    END IF



    ! (Modified) convective lightning potential index lpi_con_max, mlpi_con_max 
    ! and convection index koi. 
    ! If one of these is requested in the output_nml the diagnosis is run.
    !
    IF (var_in_output%lpi_con .OR. var_in_output%lpi_con_max .OR.             &
      &  var_in_output%mlpi_con .OR. var_in_output%mlpi_con_max .OR.          &
      &  var_in_output%koi                                                    ) THEN

      ! &      diag%lpi_con(nproma,nblks_c)
      cf_desc    = t_cf_var('lpi_con', 'J kg-1 ',                               &
       &           'subgrid-scale lightning potential index', datatype_flt)
      grib2_desc = grib2_var(0,17,5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'lpi_con', diag%lpi_con   ,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & hor_interp=create_hor_interp_metadata(                      &
                  & hor_intp_type=HINTP_TYPE_LONLAT_NNB), lopenacc=.TRUE. )
      

      ! For the time being we use as reset intervals:
      !  0-48 h - hourly
      ! 48-72 h - 3 hourly
      ! 72-   h - 6 hourly
      ! maximisation
      lpi_int1(:) = ' '
      lpi_int2(:) = ' '
      lpi_int3(:) = ' '
      lpi_start1(:) = ' '
      lpi_start2(:) = ' '
      lpi_start3(:) = ' '
      lpi_end1(:) = ' '
      lpi_end2(:) = ' '
      CALL getPTStringFromMS(3600000_i8, lpi_int1)
      CALL getPTStringFromMS(10800000_i8, lpi_int2)
      CALL getPTStringFromMS(21600000_i8, lpi_int3)
      CALL getPTStringFromMS(0_i8, lpi_start1)
      CALL getPTStringFromMS(172800000_i8, lpi_start2)
      CALL getPTStringFromMS(259200000_i8, lpi_start3)
      CALL getPTStringFromMS(172800000_i8, lpi_end1)
      CALL getPTStringFromMS(259200000_i8, lpi_end2)


      ! &      diag%lpi_con_max(nproma,nblks_c)
      cf_desc    = t_cf_var('lpi_con_max', 'J kg-1 ',                          &
          &  'subgrid-scale lightning potential index, maximum since end of previous full ' &
          &//'01H (- +48h), 03H (+48 - +72h) and 06h (+72h -) since model start', datatype_flt)
      grib2_desc = grib2_var(0,17,5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'lpi_con_max', diag%lpi_con_max   ,             &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                    & cf_desc, grib2_desc,                                     &
                    & ldims=shape2d,                                           &
                    & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                    & resetval=0.0_wp, initval=0.0_wp,                         &
                    & action_list=actions(                                     &
                    &    new_action( ACTION_RESET, TRIM(lpi_int1),             &
                    &                TRIM(lpi_start1), TRIM(lpi_end1) ),       &
                    &    new_action( ACTION_RESET, TRIM(lpi_int2),             &
                    &                TRIM(lpi_start2), TRIM(lpi_end2) ),       &
                    &    new_action( ACTION_RESET, TRIM(lpi_int3),             &
                    &                TRIM(lpi_start3)                 )     ), &
        & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
        &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                    & lopenacc=.TRUE.)
      


      ! &      diag%mlpi_con(nproma,nblks_c)
      cf_desc    = t_cf_var('mlpi_con', 'J kg-1 ', 'modified lightning potential index', datatype_flt)
      grib2_desc = grib2_var(0, 17, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mlpi_con', diag%mlpi_con   ,                    &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE.)
      


      ! &      diag%mlpi_con_max(nproma,nblks_c)
      cf_desc    = t_cf_var('mlpi_con_max', 'J kg-1 ',                         &
          &  'modified lightning potential index, maximum since end of previous full '    &
          &//'01H (- +48h), 03H (+48 - +72h) and 06h (+72h -) since model start', datatype_flt)
      grib2_desc = grib2_var(0, 17, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'mlpi_con_max', diag%mlpi_con_max   ,           &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                    & cf_desc, grib2_desc,                                     &
                    & ldims=shape2d,                                           &
                    & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                    & resetval=0.0_wp, initval=0.0_wp,                         &
                    & action_list=actions(                                     &
                    &    new_action( ACTION_RESET, TRIM(lpi_int1),             &
                    &                TRIM(lpi_start1), TRIM(lpi_end1) ),       &
                    &    new_action( ACTION_RESET, TRIM(lpi_int2),             &
                    &                TRIM(lpi_start2), TRIM(lpi_end2) ),       &
                    &    new_action( ACTION_RESET, TRIM(lpi_int3),             &
                    &                TRIM(lpi_start3)                 )     ), &
        & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
        &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                  & lopenacc=.TRUE.)
      


      ! &      diag%koi(nproma,nblks_c)
      cf_desc    = t_cf_var('koi', 'K', 'convection index', datatype_flt)
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'koi', diag%koi   ,                         &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_NNB), lopenacc=.TRUE. )
      

    ELSE ! dummy allocation
      ALLOCATE(diag%lpi_con (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%lpi_con failed')
      !
      ALLOCATE(diag%lpi_con_max (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%lpi_con_max failed')
      !
      ALLOCATE(diag%mlpi_con (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%mlpi_con failed')
      !
      ALLOCATE(diag%mlpi_con_max (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%mlpi_con_max failed')
      !
      ALLOCATE(diag%koi (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%koi failed')

      !$ACC ENTER DATA CREATE(diag%lpi_con, diag%lpi_con_max, diag%mlpi_con, diag%mlpi_con_max, diag%koi)
    ENDIF



    ! Lightning flash density lfd_con, lfd_con_max. 
    ! If one of these fields is requested in the output_nml the diagnosis is run.
    !
    IF ( var_in_output%lfd_con .OR.  var_in_output%lfd_con_max ) THEN

      ! &      diag%lfd_con(nproma,nblks_c)
      cf_desc    = t_cf_var('lfd_con', 'km-2 day-1', 'lightning flash density km-2 day-1', datatype_flt)
      grib2_desc = grib2_var(0,17,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)        &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 8)
      CALL add_var( diag_list, 'lfd_con', diag%lfd_con   ,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & isteptype=TSTEP_INSTANT,                                    &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_NNB), lopenacc=.TRUE. )
      


      ! For the time being we use as reset intervals:
      !  0-48 h - hourly
      ! 48-72 h - 3 hourly
      ! 72-   h - 6 hourly
      ! maximisation
      lpi_int1(:) = ' '
      lpi_int2(:) = ' '
      lpi_int3(:) = ' '
      lpi_start1(:) = ' '
      lpi_start2(:) = ' '
      lpi_start3(:) = ' '
      lpi_end1(:) = ' '
      lpi_end2(:) = ' '
      CALL getPTStringFromMS(3600000_i8, lpi_int1)
      CALL getPTStringFromMS(10800000_i8, lpi_int2)
      CALL getPTStringFromMS(21600000_i8, lpi_int3)
      CALL getPTStringFromMS(0_i8, lpi_start1)
      CALL getPTStringFromMS(172800000_i8, lpi_start2)
      CALL getPTStringFromMS(259200000_i8, lpi_start3)
      CALL getPTStringFromMS(172800000_i8, lpi_end1)
      CALL getPTStringFromMS(259200000_i8, lpi_end2)


      ! &      diag%lfd_con_max(nproma,nblks_c)
      cf_desc    = t_cf_var('lfd_con_max', 'km-2 day-1',                       &
          &  'maximum lightning flash density km-2 day-1 since end of previous full '     &
          &//'01H (- +48h), 03H (+48 - +72h) and 06h (+72h -) since model start', datatype_flt)
      grib2_desc = grib2_var(0,17,4, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)          &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 8)
      CALL add_var( diag_list, 'lfd_con_max', diag%lfd_con_max   ,             &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                    & cf_desc, grib2_desc,                                     &
                    & ldims=shape2d,                                           &
                    & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                    & resetval=0.0_wp, initval=0.0_wp,                         &
                    & action_list=actions(                                     &
                    &    new_action( ACTION_RESET, TRIM(lpi_int1),             &
                    &                TRIM(lpi_start1), TRIM(lpi_end1) ),       &
                    &    new_action( ACTION_RESET, TRIM(lpi_int2),             &
                    &                TRIM(lpi_start2), TRIM(lpi_end2) ),       &
                    &    new_action( ACTION_RESET, TRIM(lpi_int3),             &
                    &                TRIM(lpi_start3)                 )     ), &
        & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_BCTR, &
        &                                       fallback_type=HINTP_TYPE_LONLAT_RBF), &
                    & lopenacc=.TRUE.)
      
    ELSE ! dummy allocation
      ALLOCATE(diag%lfd_con (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%lfd_con failed')
      !
      ALLOCATE(diag%lfd_con_max (0,kblks), STAT=ist)
      IF (ist/=SUCCESS)  CALL finish(TRIM(routine), 'dummy allocation for diag%lfd_con_max failed')

      !$ACC ENTER DATA CREATE(diag%lfd_con, diag%lfd_con_max)
    ENDIF



    IF (var_in_output%ceiling) THEN
      cf_desc    = t_cf_var('ceiling', 'm', 'ceiling height', datatype_flt)
      grib2_desc = grib2_var(0, 6, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)          &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
      CALL add_var( diag_list,                                                       &
                    & "ceiling", diag%ceiling_height,                                &
                    & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,                         &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_CEILING, lrestart=.FALSE.,    &
                    & lopenacc=.TRUE. )
    
    END IF

    IF (var_in_output%vis) THEN
      cf_desc    = t_cf_var('vis', 'm', 'near surface visibility', datatype_flt)
      grib2_desc = grib2_var(0, 19, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)          
      CALL add_var( diag_list,                                                       &
                    & "vis", diag%vis,                                               &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_VIS, lrestart=.FALSE.,        &
                    & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%hbas_sc) THEN
      cf_desc    = t_cf_var('hbas_sc', 'm', 'cloud base above msl, shallow convection', datatype_flt)
      grib2_desc = grib2_var(0, 6, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
        &              + t_grib2_int_key("typeOfSecondFixedSurface", 101)
      CALL add_var( diag_list,                                                       &
                    & "hbas_sc", diag%hbas_sc,                                       &
                    & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,                         &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_HBAS_SC, lrestart=.FALSE.,    &
                    & lopenacc=.TRUE. )
      
    END IF


    IF (var_in_output%inversion_height) THEN

      cf_desc    = t_cf_var('inversion_height', 'm', 'lowest inversion height', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)     
      CALL add_var( diag_list,                                                       &
                    & "inversion_height", diag%inversion_height,                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP,                          &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_INVERSION, lrestart=.FALSE.,  &
                    & lopenacc=.TRUE. )
      

      cf_desc    = t_cf_var('low_ent_zone', 'm', &
      &     'lowest limit of the entrainment zone corresponding to the lowest inversion', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "low_ent_zone", diag%low_ent_zone,                             &
                    & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP,                          &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & lrestart=.FALSE.,                                              &
                    & lopenacc=.TRUE. )
      

    END IF


    IF (var_in_output%htop_sc) THEN
      cf_desc    = t_cf_var('htop_sc', 'm', 'cloud top above msl, shallow convection', datatype_flt)
      grib2_desc = grib2_var(0, 6, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
        &              + t_grib2_int_key("typeOfSecondFixedSurface", 101)
      CALL add_var( diag_list,                                                       &
                    & "htop_sc", diag%htop_sc,                                       &
                    & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_TOP,                          &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_HTOP_SC, lrestart=.FALSE.,    &
                    & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%twater) THEN
      cf_desc    = t_cf_var('twater', 'kg m-2', 'Total column integrated water', datatype_flt)
      grib2_desc = grib2_var(0, 1, 78, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "twater", diag%twater,                                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_TWATER, lrestart=.FALSE.,     &
                    & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%q_sedim) THEN
      cf_desc    = t_cf_var('q_sedim', 'kg kg-1', 'Specific content of precipitation particles', datatype_flt)
      grib2_desc = grib2_var(0, 1, 196, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
        &              + t_grib2_int_key("typeOfSecondFixedSurface", 150)
      CALL add_var( diag_list,                                                       &
                    & "q_sedim", diag%q_sedim,                                       &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_Q_SEDIM, lrestart=.FALSE.,    &
                    & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%tcond_max) THEN
      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)
      cf_desc    = t_cf_var('tcond_max', 'kg m-2',                           &
        &                   'total column-integrated condensate, max. since end of previous full '// &
        &                   TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 1, 81, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &              + t_grib2_int_key("typeOfSecondFixedSurface", 8)
      CALL add_var( diag_list, 'tcond_max', diag%tcond_max,                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                  & resetval=0.0_wp, initval=0.0_wp,                         &
                  & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%tcond10_max) THEN
      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)
      cf_desc    = t_cf_var('tcond10_max', 'kg m-2',                         &
        &                   'total column-integrated condensate above z(T=-10 degC), max. since end of previous full '// &
        &                   TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 1, 81, ibits, GRID_UNSTRUCTURED, GRID_CELL)       &
        &              + t_grib2_int_key("typeOfFirstFixedSurface",           20)  &
        &              + t_grib2_int_key("typeOfSecondFixedSurface",           8)  &
        &              + t_grib2_int_key("scaledValueOfFirstFixedSurface", 26315)  &
        &              + t_grib2_int_key("scaleFactorOfFirstFixedSurface",     2)
      CALL add_var( diag_list, 'tcond10_max', diag%tcond10_max,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                  & resetval=0.0_wp, initval=0.0_wp,                         &
                  & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                  & lopenacc=.TRUE. )
      
    END IF

    luh_max_out(k_jg, :) = (/var_in_output%uh_max_low, var_in_output%uh_max_med, var_in_output%uh_max/)
    uh_max_zmin          = (/                   0._wp,                 2000._wp,             2000._wp/)
    uh_max_zmax          = (/                3000._wp,                 5000._wp,             8000._wp/)

    IF (ANY(luh_max_out(k_jg, :))) THEN

      cf_desc    = t_cf_var('uh_max_3d', 'm2 s-2', 'dummy', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'uh_max_3d', diag%uh_max_3d,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                  &
                  & cf_desc, grib2_desc, ldims=shape3d_uh_max,           &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,&
                  & lopenacc=.TRUE. )
      

      ALLOCATE( diag%uh_max_ptr(uh_max_nlayer) )

      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)

      DO k = 1,uh_max_nlayer

        SELECT CASE (k)
          CASE (1)
            shortname   = 'uh_max_low'
          CASE (2)
            shortname   = 'uh_max_med'
          CASE (3)
            shortname   = 'uh_max'
        END SELECT

        longname = 'updraft helicity '// &
                   TRIM(real2string(uh_max_zmin(k)))//'-'//TRIM(real2string(uh_max_zmax(k)))// &
                   ' m, max. since end of previous full '//TRIM(celltracks_int(3:))//' since model start'

        IF (luh_max_out(k_jg, k)) THEN
          cf_desc    = t_cf_var(TRIM(shortname), 'm2 s-2', TRIM(longname), datatype_flt)
          grib2_desc = grib2_var( 0, 7, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)                  &
            &           + t_grib2_int_key("typeOfFirstFixedSurface",          102)                &
            &           + t_grib2_int_key("typeOfSecondFixedSurface",         102)                &
            &           + t_grib2_int_key("scaledValueOfFirstFixedSurface",  int(uh_max_zmin(k))) &
            &           + t_grib2_int_key("scaledValueOfSecondFixedSurface", int(uh_max_zmax(k)))
          CALL add_ref( diag_list, 'uh_max_3d', TRIM(shortname), diag%uh_max_ptr(k)%p_2d, &
                      & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                               &
                      & cf_desc, grib2_desc, ref_idx=k, ldims=shape2d,                    &
                      & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,             &
                      & resetval=0.0_wp, initval=0.0_wp,                                  &
                      & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ) )
        END IF
      END DO

    END IF ! ANY(luh_max_out)

    IF (var_in_output%vorw_ctmax) THEN
      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)
      cf_desc    = t_cf_var('vorw_ctmax', 's-1',                   &
        &                   'Maximum rotation amplitude since end of previous full '// &
        &                   TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 206, ibits, GRID_UNSTRUCTURED, GRID_CELL)   &
        &           + t_grib2_int_key("typeOfFirstFixedSurface",          102)  &
        &           + t_grib2_int_key("typeOfSecondFixedSurface",         102)  &
        &           + t_grib2_int_key("scaledValueOfFirstFixedSurface",     0)  &
        &           + t_grib2_int_key("scaledValueOfSecondFixedSurface", 3000)
      CALL add_var( diag_list, 'vorw_ctmax', diag%vorw_ctmax,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                  & resetval=0.0_wp, initval=0.0_wp,                         &
                  & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%w_ctmax) THEN
      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)
      cf_desc    = t_cf_var('w_ctmax', ' m s-1',                   &
        &                   'Maximum updraft track since end of previous full '// &
        &                   TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 207, ibits, GRID_UNSTRUCTURED, GRID_CELL)    &
        &           + t_grib2_int_key("typeOfFirstFixedSurface",           102)  &
        &           + t_grib2_int_key("typeOfSecondFixedSurface",          102)  &
        &           + t_grib2_int_key("scaledValueOfFirstFixedSurface",      0)  &
        &           + t_grib2_int_key("scaledValueOfSecondFixedSurface", 10000)
      CALL add_var( diag_list, 'w_ctmax', diag%w_ctmax,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,    &
                  & resetval=0.0_wp, initval=0.0_wp,                         &
                  & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ) , &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%wshear_u) THEN
      cf_desc    = t_cf_var('wshear_u', 'm s-1', 'U-component of vertical wind shear vector', datatype_flt)
      grib2_desc = grib2_var(0, 2, 208, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "wshear_u", diag%wshear_u,                                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_WSHEAR,                             &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3dwshear,                                           &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_WSHEAR_U, lrestart=.FALSE. )
    END IF
    
    IF (var_in_output%wshear_v) THEN
      cf_desc    = t_cf_var('wshear_v', 'm s-1', 'V-component of vertical wind shear vector', datatype_flt)
      grib2_desc = grib2_var(0, 2, 209, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "wshear_v", diag%wshear_v,                                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_WSHEAR,                             &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3dwshear,                                           &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_WSHEAR_V, lrestart=.FALSE. )
    END IF
    
    IF (var_in_output%lapserate) THEN
      cf_desc    = t_cf_var('lapse_rate', 'K/m', 'Temperature lapse rate 500hPa - 850hPa', datatype_flt)
      grib2_desc = grib2_var(0, 0, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "lapse_rate", diag%lapse_rate,                                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_PRESSURE_LAPSERATE,                 &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape2d,                                                 &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_LAPSERATE, lrestart=.FALSE.)
    END IF

    IF (var_in_output%srh) THEN
      cf_desc    = t_cf_var('srh', 'm2 s-2', 'Storm relative helicity', datatype_flt)
      grib2_desc = grib2_var(0, 7, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list,                                                       &
                    & "srh", diag%srh,                                               &
                    & GRID_UNSTRUCTURED_CELL, ZA_SRH,                                &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3dsrh,                                              &
                    & isteptype=TSTEP_INSTANT,                                       &
                    & l_pp_scheduler_task=TASK_COMPUTE_SRH, lrestart=.FALSE.)
    END IF

    IF (var_in_output%cloudtop) THEN
      ! &      diag%cloudtop(nproma,nblks_c)
      cf_desc    = t_cf_var('cloudtop', 'm', 'Cloud top height', datatype_flt)
      grib2_desc = grib2_var(0, 6, 12, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'cloudtop', diag%cloudtop,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & in_group=groups("additional_precip_vars"),                  &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      
    ENDIF

    IF (var_in_output%si .OR. var_in_output%sli .OR. var_in_output%swiss12 .OR. var_in_output%swiss00) THEN
      ! &      diag%si(nproma,nblks_c)
      cf_desc    = t_cf_var('si', 'K', 'Showalter Index', datatype_flt)
      grib2_desc = grib2_var(0, 7, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'si', diag%si,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & in_group=groups("additional_precip_vars"),                  &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      

      ! &      diag%sli(nproma,nblks_c)
      cf_desc    = t_cf_var('sli', 'K', 'Surface Lifted Index', datatype_flt)
      grib2_desc = grib2_var(0, 7, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 10)
      CALL add_var( diag_list, 'sli', diag%sli,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & in_group=groups("additional_precip_vars"),                  &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      

      ! &      diag%swiss12(nproma,nblks_c)
      cf_desc    = t_cf_var('swiss12', 'K', 'SWISS12 Index', datatype_flt)
      grib2_desc = grib2_var(215, 7, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( diag_list, 'swiss12', diag%swiss12,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & in_group=groups("additional_precip_vars"),                  &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      

      ! &      diag%swiss00(nproma,nblks_c)
      cf_desc    = t_cf_var('swiss00', 'K', 'SWISS00 Index', datatype_flt)
      grib2_desc = grib2_var(215, 7, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( diag_list, 'swiss00', diag%swiss00,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
                  & in_group=groups("additional_precip_vars"),                  &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &    fallback_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      
    ENDIF



    IF (var_in_output%cape_3km .OR. var_in_output%lcl_ml .OR. var_in_output%lfc_ml) THEN
      ! &      diag%lcl_ml(nproma,nblks_c)
      ! typeOfLevel ZA_SURFACE is changed to 192 in vlistDefVarIntKey
      cf_desc    = t_cf_var('lcl_ml', 'm ', 'Lifted Condensation Level of Mean Layer parcel (HAG)', datatype_flt)
      grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 5)                  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 192)
      CALL add_var( diag_list, 'lcl_ml', diag%lcl_ml,                           &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
  !!$                & lmiss=.TRUE., missval=-999.9_wp,                            &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      

          ! &      diag%lfc_ml(nproma,nblks_c)
      ! typeOfLevel ZA_SURFACE is changed to 192 in vlistDefVarIntKey
      cf_desc    = t_cf_var('lfc_ml', 'm ', 'Level of Free Convection of Mean Layer parcel (HAG)', datatype_flt)
      grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)         &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 194)                &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 192)
      CALL add_var( diag_list, 'lfc_ml', diag%lfc_ml,                           &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d, lrestart=.FALSE.,                            &
  !!$                & lmiss=.TRUE., missval=-999.9_wp,                            &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                  & lopenacc=.TRUE. )
      
      
      cf_desc    = t_cf_var('cape_3km', 'J kg-1 ',                                 &
            'Mean Layer CAPE, with endpoint at 3000 m', datatype_flt)
      grib2_desc = grib2_var(0, 7, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)        &
            &           + t_grib2_int_key("typeOfFirstFixedSurface", 103)         &
            &           + t_grib2_int_key("typeOfSecondFixedSurface", 192)        &
            &           + t_grib2_int_key("scaledValueOfFirstFixedSurface", 3000) &
            &           + t_grib2_int_key("scaleFactorOfFirstFixedSurface", 0)
      CALL add_var( diag_list, 'cape_3km', diag%cape_3km,                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                         &
                    & cf_desc, grib2_desc,                                        &
                    & ldims=shape2d,                                              &
                    & isteptype=TSTEP_INSTANT,                                    &
                    & hor_interp=create_hor_interp_metadata(                      &
                    &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                    & lrestart=.FALSE., lopenacc=.TRUE.)
      
    ENDIF

    IF (var_in_output%cape_mu .OR. var_in_output%cin_mu) THEN
      cf_desc    = t_cf_var('cape_mu', 'J kg-1 ',                                 &
           'Most unstable CAPE, approximated by parcel with largest Theta_e in 3000 m layer', datatype_flt)
      grib2_desc = grib2_var(0, 7, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)        &
           &           + t_grib2_int_key("typeOfFirstFixedSurface", 193)
      CALL add_var( diag_list, 'cape_mu', diag%cape_mu,                           &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                         &
                    & cf_desc, grib2_desc,                                        &
                    & ldims=shape2d,                                              &
                    & isteptype=TSTEP_INSTANT,                                    &
                    & hor_interp=create_hor_interp_metadata(                      &
                    &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                    & lrestart=.FALSE., lopenacc=.TRUE.)
      

      cf_desc    = t_cf_var('cin_mu', 'J kg-1 ',                                  &
           'Most unstable convective inhibition, approximated by parcel with largest Theta_e in 3000 m layer', datatype_flt)
      grib2_desc = grib2_var(0, 7, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)        &
           &           + t_grib2_int_key("typeOfFirstFixedSurface", 193)
      CALL add_var( diag_list, 'cin_mu', diag%cin_mu,                             &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                         &
                    & cf_desc, grib2_desc,                                        &
                    & ldims=shape2d,                                              &
                    & isteptype=TSTEP_INSTANT,                                    &
                    & hor_interp=create_hor_interp_metadata(                      &
                    &    hor_intp_type=HINTP_TYPE_LONLAT_NNB),                    &
                    & lrestart=.FALSE., lopenacc=.TRUE.)
      
    END IF


    IF (var_in_output%dbz .OR. var_in_output%dbz850 .OR. var_in_output%dbzcmax .OR. var_in_output%dbzctmax .OR. &
         var_in_output%dbzlmx_low .OR. var_in_output%echotop .OR. var_in_output%echotopinm) THEN
      cf_desc    = t_cf_var('dbz', 'dBZ',&
        &                   'radar reflectivity in dBZ', datatype_flt)
      grib2_desc = grib2_var( 0, 15, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      ! internal storage not in dBZ, but linear! dBZ conversion is a post-op before output!
      ! It is explicitly computed in mo_nh_stepping.f90, so no l_pp_scheduler_task necessary
      CALL add_var( diag_list, 'dbz', diag%dbz3d_lin,                          &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
                  & ldims=shape3d, loutput=.TRUE., lrestart=.FALSE.,           &
                  & vert_interp=create_vert_interp_metadata(                   &
                  &             vert_intp_type=vintp_types("P","Z","I"),       &
                  &             vert_intp_method=VINTP_METHOD_LIN,             &
                  &             l_loglin=.FALSE., l_pd_limit=.FALSE.,          &
                  &             l_extrapol=.FALSE., lower_limit=0._wp),        &
                  & post_op=post_op(POST_OP_LIN2DBZ, arg1=1e-15_wp),           &
                  & lopenacc=.TRUE. )
      
    END IF
    
    IF (var_in_output%dbzcmax) THEN
      cf_desc    = t_cf_var('dbz_cmax', 'dBZ',                            &
        &                   'Column maximum reflectivity', datatype_flt)
      grib2_desc = grib2_var( 0, 15, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface",           1) &
        &           + t_grib2_int_key("typeOfSecondFixedSurface",          8)
      CALL add_var( diag_list, 'dbz_cmax', diag%dbz_cmax,                    &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.FALSE., loutput=.TRUE.,                        &
                  & l_pp_scheduler_task=TASK_COMPUTE_DBZCMAX,                &
                  & post_op=post_op(POST_OP_LIN2DBZ, arg1=1e-15_wp),         &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%dbz850) THEN
      cf_desc    = t_cf_var('dbz_850', 'dBZ',                             &
        &                   'Reflectivity in approx. 850 hPa', datatype_flt)
      grib2_desc = grib2_var( 0, 15, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface",           1)
      CALL add_var( diag_list, 'dbz_850', diag%dbz_850,                      &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.FALSE., loutput=.TRUE.,                        &
                  & l_pp_scheduler_task=TASK_COMPUTE_DBZ850,                 &
                  & post_op=post_op(POST_OP_LIN2DBZ, arg1=1e-15_wp),         &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%dbzlmx_low) THEN
      ! NOTE: In eccodes the shortName DBZLMX_LOW is associated with the fixed bounds
      !       500 m and 2500 m. But in practice the bounds 1000 m and 2000 m were
      !       found to be more appropriate to approximate typical radar composites.
      !       However, there might be further adaptions in the near future, and
      !       to avoid several consecutive adaptions of eccodes until consolidation
      !       of the layer bounds, we for now set the "real" bounds to 1000 m and 2000 m
      !       in mo_pp_tasks.f90, but we do not change them here in order not to
      !       break the connection to the eccodes shortName.
      cf_desc    = t_cf_var('dbzlmx_low', 'dBZ',                             &
        &                   'Reflectivity layer maximum 500 - 2500 m AGL', datatype_flt)
      grib2_desc = grib2_var( 0, 15, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface",          103)  &
        &           + t_grib2_int_key("typeOfSecondFixedSurface",         103)  &
        &           + t_grib2_int_key("scaleFactorOfFirstFixedSurface",     0)  &
        &           + t_grib2_int_key("scaledValueOfFirstFixedSurface",  2500)  &
        &           + t_grib2_int_key("scaleFactorOfSecondFixedSurface",    0)  &
        &           + t_grib2_int_key("scaledValueOfSecondFixedSurface", 500)
      CALL add_var( diag_list, 'dbzlmx_low', diag%dbzlmx_low,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                      &
                  & cf_desc, grib2_desc,                                     &
                  & ldims=shape2d,                                           &
                  & lrestart=.FALSE., loutput=.TRUE.,                        &
                  & l_pp_scheduler_task=TASK_COMPUTE_DBZLMX_LOW,             &
                  & post_op=post_op(POST_OP_LIN2DBZ, arg1=1e-15_wp),         &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%dbzctmax) THEN
      celltracks_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*celltracks_interval(k_jg), i8), celltracks_int)
      cf_desc    = t_cf_var('dbz_ctmax', 'dBZ',                   &
        &                   'Column and time maximum reflectivity since end of previous full '// &
        &                   TRIM(celltracks_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 15, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)    &
        &           + t_grib2_int_key("typeOfFirstFixedSurface",           1)   &
        &           + t_grib2_int_key("typeOfSecondFixedSurface",          8)
      CALL add_var( diag_list, 'dbz_ctmax', diag%dbz_ctmax,                     &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                         &
                  & cf_desc, grib2_desc,                                        &
                  & ldims=shape2d,                                              &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,       &
                  & resetval=0.0_wp, initval=0.0_wp,                            &
                  & action_list=actions( new_action( ACTION_RESET, celltracks_int ) ), & 
                  & post_op=post_op(POST_OP_LIN2DBZ, arg1=1e-15_wp),            &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%echotop) THEN
      echotop_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*echotop_meta(k_jg)%time_interval, i8), echotop_int)
      cf_desc    = t_cf_var('echotop', 'Pa',                   &
        &                   'Minimum pressure of exceeding radar reflectivity threshold since end of previous full '// &
        &                   TRIM(echotop_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)     &
           &        + t_grib2_int_key("typeOfFirstFixedSurface",          25)
!!$ This seems not to have any effect for creating a bitmap with missing values:
!!$         &        + t_grib2_dbl_key("missingValue",                 -999d0)   &
!!$         &        + t_grib2_int_key("bitmapPresent",                     1)
      CALL add_var( diag_list, 'echotop', diag%echotop,                         &
                  & GRID_UNSTRUCTURED_CELL, ZA_ECHOTOP,                         &
                  & cf_desc, grib2_desc,                                        &
                  & ldims=shape3dechotop,                                       &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MIN,       &
                  & resetval=-999.0_wp, initval=-999.0_wp,                      &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &            hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
!!$                  & lmiss=.TRUE., missval=-999.0_wp,                            &
                  & action_list=actions( new_action( ACTION_RESET, echotop_int ) ), &
                  & lopenacc=.TRUE. )
      
    END IF

    IF (var_in_output%echotopinm) THEN
      echotop_int(:) = ' '
      CALL getPTStringFromMS(NINT(1000*echotop_meta(k_jg)%time_interval, i8), echotop_int)
      cf_desc    = t_cf_var('echotopinm', 'm',                                  &
        &                   'Maximum height of exceeding radar reflectivity threshold since end of previous full '// &
        &                   TRIM(echotop_int(3:))//' since model start', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)     &
           &        + t_grib2_int_key("typeOfFirstFixedSurface",          25)
!!$ This seems not to have any effect for creating a bitmap with missing values:
!!$         &        + t_grib2_dbl_key("missingValue",                 -999d0)   &
!!$         &        + t_grib2_int_key("bitmapPresent",                     1)
      CALL add_var( diag_list, 'echotopinm', diag%echotopinm,                   &
                  & GRID_UNSTRUCTURED_CELL, ZA_ECHOTOP,                         &
                  & cf_desc, grib2_desc,                                        &
                  & ldims=shape3dechotop,                                       &
                  & lrestart=.TRUE., loutput=.TRUE., isteptype=TSTEP_MAX,       &
                  & resetval=-999.0_wp, initval=-999.0_wp,                      &
                  & hor_interp=create_hor_interp_metadata(                      &
                  &            hor_intp_type=HINTP_TYPE_LONLAT_NNB),            &
!!$                  & lmiss=.TRUE., missval=-999.0_wp,                            &
                  & action_list=actions( new_action( ACTION_RESET, echotop_int ) ), &
                  & lopenacc=.TRUE. )
      
    END IF


    IF (var_in_output%dursun_m .OR. var_in_output%dursun_r) THEN
      ! &      diag%dursun_m(nproma,nblks_c)
      cf_desc    = t_cf_var('maximum_duration_of_sunshine', 's', 'Possible astronomical maximum of sunshine', datatype_flt)
      grib2_desc = grib2_var(0, 6, 205, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
        &        + t_grib2_int_key("typeOfStatisticalProcessing", 11)
      CALL add_var( diag_list, 'dursun_m', diag%dursun_m,                       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
        &           ldims=shape2d,                                              &
        &           lrestart=.TRUE., isteptype=TSTEP_ACCUM ,                    &
        &           initval=0._wp, resetval=0._wp,                              &
        &           action_list=actions(new_action(ACTION_RESET,sunshine_interval(k_jg))), &
        &           lopenacc=.TRUE. )
        
      ! &      diag%dursun_r(nproma,nblks_c)
      cf_desc    = t_cf_var('relative_duration_of_sunshine', '%', 'relative duration of sunshine', datatype_flt)
      grib2_desc = grib2_var(0, 6, 206, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'dursun_r', diag%dursun_r,                       &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
        &           ldims=shape2d,                                              &
        &           lrestart=.TRUE., isteptype=TSTEP_ACCUM ,                    &
        &           initval=0._wp, resetval=0._wp,                              &
        &           action_list=actions(new_action(ACTION_RESET,sunshine_interval(k_jg))), &
        &           lopenacc=.TRUE.)
        
    END IF
    IF (var_in_output%dursun) THEN
      ! &      diag%dursun(nproma,nblks_c)
      cf_desc    = t_cf_var('duration_of_sunshine', 's', 'sunshine duration', datatype_flt)
      grib2_desc = grib2_var(0, 6, 33, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'dursun', diag%dursun,                           &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
        &           ldims=shape2d,                                              &
        &           lrestart=.TRUE., isteptype=TSTEP_ACCUM ,                    &
        &           initval=0._wp, resetval=0._wp,                              &
        &           action_list=actions(new_action(ACTION_RESET,sunshine_interval(k_jg))), &
        &           lopenacc=.TRUE.)
        
    END IF



    !  Height of 0 deg C level
    cf_desc    = t_cf_var('hzerocl', '', 'height of 0 deg C level', datatype_flt)
    grib2_desc = grib2_var(0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)             &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    IF (itype_hzerocl == 2) THEN
      CALL add_var( diag_list, 'hzerocl', diag%hzerocl,                                &
        &           GRID_UNSTRUCTURED_CELL, ZA_ISOTHERM_ZERO, cf_desc, grib2_desc,     &
        &           ldims=shape2d, lrestart=.FALSE.,                                   &
        &           lmiss=.TRUE., missval=-999.0_wp,                                   &
        &           hor_interp=create_hor_interp_metadata(                             &
        &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB), lopenacc=.TRUE. )
    ELSE
      CALL add_var( diag_list, 'hzerocl', diag%hzerocl,                                &
        &           GRID_UNSTRUCTURED_CELL, ZA_ISOTHERM_ZERO, cf_desc, grib2_desc,     &
        &           ldims=shape2d, lrestart=.FALSE.,                                   &
        &           hor_interp=create_hor_interp_metadata(                             &
        &                      hor_intp_type=HINTP_TYPE_LONLAT_NNB), lopenacc=.TRUE. )
    END IF
    


    !  significant weather WW
    cf_desc    = t_cf_var('ww', '', 'significant weather', datatype_flt)
    grib2_desc = grib2_var(0, 19, 25, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'ww', diag%iww,                                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE.,                                   &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB), &
                & lopenacc=.TRUE. )
    
      
    ! mask field to distinguish between tropics and extratropics (for tuning purposes)
    cf_desc    = t_cf_var('tropics_mask', '', 'tropics_mask', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'tropics_mask', diag%tropics_mask,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,                  &
                & lopenacc=.TRUE. )
    

    ! mask field to distinguish between inner tropics and elsewhere (for tuning purposes)
    cf_desc    = t_cf_var('innertropics_mask', '', 'innertropics_mask', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'innertropics_mask', diag%innertropics_mask,            &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,                  &
                & lopenacc=.TRUE. )
    

    !  mask field that is used for latitude-dependent SSO tuning parameters
    cf_desc    = t_cf_var('sso_lat_mask', '', 'sso_lat_mask', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, 'sso_lat_mask', diag%sso_lat_mask,                      &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
                & ldims=shape2d, lrestart=.FALSE., loutput=.FALSE.,                  &
                & lopenacc=.TRUE. )
    

    ! buffer field needed for the combination of vertical nesting with a reduced radiation grid
    IF (k_jg > n_dom_start) THEN
      cf_desc    = t_cf_var('buffer_rrg', '', 'buffer_rrg', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'buffer_rrg', diag%buffer_rrg,            &
                  & GRID_UNSTRUCTURED_CELL, -1, cf_desc, grib2_desc,     & ! use -1 as hackish default value for vgrid
                  & ldims=(/nproma, 3*nexlevs_rrg_vnest, p_patch_local_parent(k_jg)%nblks_c/), &
                  & lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE. )
      
    ENDIF

    ! buffer field needed for the combination of vertical nesting with a reduced radiation grid
    IF (lsynsat(k_jg)) THEN
      IF  ((k_jg > n_dom_start) .AND. (p_patch(k_jg)%nshift > 0)) THEN
        cf_desc    = t_cf_var('buffer_rttov', '', 'buffer_rttov', datatype_flt)
        grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( diag_list, 'buffer_rttov', diag%buffer_rttov,        &
                    & GRID_UNSTRUCTURED_CELL, -1, cf_desc, grib2_desc,     & ! use -1 as hackish default value for vgrid
                    & ldims=(/nproma, 5*p_patch(k_jg)%nshift, p_patch_local_parent(k_jg)%nblks_c/), &
                    & lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
        
      ENDIF
      
      shape3d_synsat = (/nproma, num_images, p_patch(k_jg)%nblks_c /)
      shape2d_synsat = (/nproma,             p_patch(k_jg)%nblks_c /)

      ! introduce container variable for RTTOV synthetic satellite imagery:
      cf_desc    = t_cf_var('rttov_channels', '', '', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'rttov_channels', diag%synsat_arr,                 &
        &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
        &           ldims=shape3d_synsat ,                                        &
        &           lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
        &           var_class=CLASS_SYNSAT, lopenacc=.TRUE.)
      

      ! add reference variables for the different images:
      ALLOCATE(diag%synsat_image(num_images))

      iimage = 0
      sensor_loop: DO isens = 1, num_sensors

        DO ichan = 1,total_numchans(isens)

          DO k=1,4
            ! the following translation can be derived by gazing at the corresponding RTTOV loop
            lradiance = ((MOD(iimage,4)+1) == RTTOV_RAD_CL) .OR. ((MOD(iimage,4)+1) == RTTOV_RAD_CS)
            lcloudy   = ((MOD(iimage,4)+1) == RTTOV_BT_CL)  .OR. ((MOD(iimage,4)+1) == RTTOV_RAD_CL)
            iimage = iimage + 1
          
            IF (lradiance) THEN
              unit = "mW/cm-1/sr/sq.m"
            ELSE
              unit = "K"
            END IF

            CALL get_synsat_name(lradiance, lcloudy, ichan, shortname, longname)
            CALL get_synsat_grib_triple(lradiance, lcloudy, ichan,       &
              &                         idiscipline, icategory, inumber, &
              &                         wave_no, wave_no_scalfac)
            
            cf_desc    = t_cf_var(shortname, unit, longname, datatype_flt)
            grib2_desc = grib2_var(idiscipline, icategory, inumber, ibits, GRID_UNSTRUCTURED, GRID_CELL)   &
              &           + t_grib2_int_key("scaledValueOfCentralWaveNumber", wave_no)                  &
              &           + t_grib2_int_key("scaleFactorOfCentralWaveNumber", wave_no_scalfac)          &
              &           + t_grib2_int_key("satelliteSeries", 333)                                     &
              &           + t_grib2_int_key("satelliteNumber",  72)                                     &
              &           + t_grib2_int_key("instrumentType",  207)
            CALL add_ref( diag_list, 'rttov_channels', TRIM(shortname),                  &
              &           diag%synsat_image(iimage)%p,                                   &
              &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
              &           cf_desc, grib2_desc, ref_idx=iimage, ldims=shape2d_synsat,     &
              &           opt_var_ref_pos = 2, lrestart=.FALSE., loutput=.TRUE.,         &
              &           in_group=groups("RTTOV"), var_class=CLASS_SYNSAT )
          END DO
        END DO
      END DO sensor_loop

    ENDIF

    ! vars for global diagnostics as in src/atm_phy_echam/mo_echam_phy_memory.f90
    IF (var_in_output%tas_gmean) THEN
      cf_desc    = t_cf_var('tas_gmean', 'K', 'global mean temperature at 2m', datatype_flt,'tas_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'tas_gmean', diag%tas_gmean,            &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%rsdt_gmean) THEN
      cf_desc    = t_cf_var('rsdt_gmean', 'W m-2', 'global mean toa incident shortwave radiation', datatype_flt,'rsdt_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'rsdt_gmean', diag%rsdt_gmean,          &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%rsut_gmean) THEN
      cf_desc    = t_cf_var('rsut_gmean', 'W m-2', 'global mean toa outgoing shortwave radiation', datatype_flt,'rsut_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'rsut_gmean', diag%rsut_gmean,          &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%rlut_gmean) THEN
      cf_desc    = t_cf_var('rlut_gmean', 'W m-2', 'global mean toa outgoing longwave radiation', datatype_flt,'rlut_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'rlut_gmean', diag%rlut_gmean,          &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%prec_gmean) THEN
      cf_desc    = t_cf_var('prec_gmean', 'kg m-2 s-1', 'global mean precipitation flux', datatype_flt,'prec_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'prec_gmean', diag%prec_gmean,          &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%evap_gmean) THEN
      cf_desc    = t_cf_var('evap_gmean', 'kg m-2 s-1', 'global mean evaporation flux', datatype_flt,'evap_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'evap_gmean', diag%evap_gmean,          &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%pme_gmean) THEN
      cf_desc    = t_cf_var('pme_gmean', 'kg m-2 s-1', 'global mean P-E', datatype_flt,'pme_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'pme_gmean', diag%pme_gmean,          &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF
    IF (var_in_output%radtop_gmean) THEN
      cf_desc    = t_cf_var('radtop_gmean', 'W m-2', 'global mean toa total radiation', datatype_flt,'radtop_gmean')
      grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_LONLAT)
      CALL add_var( diag_list, 'radtop_gmean', diag%radtop_gmean,      &
                  & GRID_LONLAT, ZA_SURFACE, cf_desc, grib2_desc,                &
                  & lrestart = .FALSE., ldims=(/1/),                             &
                  & lopenacc=.TRUE.)
      
    ENDIF

    ! Initialize JSBACH + VDIFF state.
    IF (atm_phy_nwp_config(k_jg)%inwp_surface == LSS_JSBACH .OR. atm_phy_nwp_config(k_jg)%inwp_turb == ivdiff) THEN
      CALL diag%nwp_vdiff_state%init(nproma, p_patch(k_jg)%nlev, p_patch(k_jg)%nblks_c, diag_list)
      CALL nwp_vdiff_setup( p_patch(k_jg) )
    END IF

    CALL message('mo_nwp_phy_state:construct_nwp_phy_diag', &
                 'construction of NWP physical fields finished')  

END SUBROUTINE new_nwp_phy_diag_list


SUBROUTINE new_nwp_phy_stochconv_list( k_jg, kblks,    &
                     &  listname, stochconv_list, phy_stochconv)

    INTEGER,INTENT(IN) :: kblks, k_jg !< dimension sizes

    CHARACTER(len=*),INTENT(IN)     :: listname

    TYPE(t_var_list_ptr)    ,INTENT(INOUT) :: stochconv_list
    TYPE(t_nwp_phy_stochconv),INTENT(INOUT) :: phy_stochconv

    ! Local variables
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2),shape3dstoch(3)
    INTEGER :: ibits

    INTEGER :: datatype_flt

    CHARACTER(len=*), PARAMETER :: &
      routine = 'mo_nwp_phy_state:new_nwp_phy_stochconv_list'

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! bits "entropy" of horizontal slice

    shape2d        = (/nproma,               kblks            /)
    shape3dstoch   = (/atm_phy_nwp_config(k_jg)%nclds , nproma, kblks/)

    ! Register a field list and apply default settings

    CALL vlr_add(stochconv_list, TRIM(listname), patch_id=k_jg, lrestart=.TRUE.)

    !------------------------------
    ! Meteorological quantities
    !------------------------------

    
    ! Variables for SDE stochastic convection schemes
    IF (atm_phy_nwp_config(k_jg)%lstoch_sde ) THEN
       ! stochastic shallow convection (using stochastic differential equations, SDE)
       ! number of active clouds, normalised by grid point area
       cf_desc    = t_cf_var('clnum_a', 'm-2', 'number of active clouds per unit area', datatype_flt32)
       grib2_desc = grib2_var(0, 6, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL) 
       CALL add_var( stochconv_list,                                           &
            & 'clnum_a', phy_stochconv%clnum_a,  &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,                                              &
            & loutput=.TRUE.,lrestart=.TRUE.,                             &
            & in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","iau_restore_vars"))
        
       ! perturbed cloud-base mass flux for active clouds
       cf_desc    = t_cf_var('clmf_a', 'kg s-1 m-2', 'cloud base mass flux associated with active clouds', datatype_flt32)
       grib2_desc = grib2_var(0, 3, 200, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clmf_a', phy_stochconv%clmf_a,    &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,                                              &
            & loutput=.TRUE.,lrestart=.TRUE.,                             &
            & in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","iau_restore_vars"))

       ! number of passive clouds, normalised by grid point area
       cf_desc    = t_cf_var('clnum_p', 'm-2', 'number of passive clouds per unit area', datatype_flt32)
       grib2_desc = grib2_var(0, 6, 211, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clnum_p', phy_stochconv%clnum_p,  &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,                                              &
            & loutput=.TRUE.,lrestart=.TRUE. ,                            &
            & in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","iau_restore_vars"))
       
       ! perturbed cloud-base mass flux for passive clouds
       cf_desc    = t_cf_var('clmf_p', 'kg s-1 m-2','cloud base mass flux associated with passive clouds', datatype_flt32)
       grib2_desc = grib2_var(0, 3, 199, ibits, GRID_UNSTRUCTURED, GRID_CELL) 
       CALL add_var( stochconv_list,                                          &
            & 'clmf_p', phy_stochconv%clmf_p,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape2d,                                             &
            & loutput=.TRUE.,lrestart=.TRUE.,                            &
            & in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","iau_restore_vars"))

    ELSE IF (atm_phy_nwp_config(k_jg)%lstoch_expl) THEN
       ! stochastic shallow convection (explicit)
       ! use variables only locally in cumaster, do not write for restart
       ! number of active clouds, normalised by grid point area
       cf_desc    = t_cf_var('clnum_a', 'm-2', 'number of active clouds per unit area', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clnum_a', phy_stochconv%clnum_a,  &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,lrestart=.FALSE.)
        
       ! perturbed cloud-base mass flux for active clouds
       cf_desc    = t_cf_var('clmf_a', 'kg s-1 m-2', 'cloud base mass flux associated with active clouds', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clmf_a', phy_stochconv%clmf_a,    &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,lrestart=.FALSE.)

       ! number of passive clouds, normalised by grid point area
       cf_desc    = t_cf_var('clnum_p', 'm-2', 'number of passive clouds per unit area', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clnum_p', phy_stochconv%clnum_p,  &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d, lrestart=.FALSE.)
       
       ! perturbed cloud-base mass flux for passive clouds
       cf_desc    = t_cf_var('clmf_p', 'kg s-1 m-2','cloud base mass flux associated with passive clouds', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'clmf_p', phy_stochconv%clmf_p,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape2d,lrestart=.FALSE.)

       ! explicit cloud ensemble: individual cloud's base mass flux
       cf_desc    = t_cf_var('mf_i', 'kg m-2 -s','cloud base mass flux', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'mf_i', phy_stochconv%mf_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))
       
       ! explicit cloud ensemble: individual cloud's current age
       cf_desc    = t_cf_var('time_i', 's','cloud age', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'time_i', phy_stochconv%time_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))

       ! explicit cloud ensemble: individual cloud's expected lifetime
       cf_desc    = t_cf_var('life_i', 's','cloud lifetime', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'life_i', phy_stochconv%life_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))

       ! explicit cloud ensemble: individual cloud's base updraft area
       cf_desc    = t_cf_var('area_i', 'm2','cloud updraft area', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'area_i', phy_stochconv%area_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))
       
       ! explicit cloud ensemble: individual cloud's type: active or passive
       cf_desc    = t_cf_var('type_i', '-','cloud active or passive', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'type_i', phy_stochconv%type_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))
       
       ! explicit cloud ensemble: individual cloud's convection type: shallow/deep
       cf_desc    = t_cf_var('ktype_i', '-','cloud convection type', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'ktype_i', phy_stochconv%ktype_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))

       ! explicit cloud ensemble: individual cloud's cloud base (model level index)
       cf_desc    = t_cf_var('base_i', '-','cloud base level', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'base_i', phy_stochconv%base_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))

       ! explicit cloud ensemble: individual cloud's maximum depth (in layers)
       cf_desc    = t_cf_var('depth_i', '-','cloud depth', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'depth_i', phy_stochconv%depth_i,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))

       ! explicit cloud ensemble: is this space in the ensemble used or not (logical)
       cf_desc    = t_cf_var('used_cell', '-','cell use tag', datatype_flt32)
       grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                          &
            & 'used_cell', phy_stochconv%used_cell,   &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE,cf_desc, grib2_desc, &
            & ldims=shape3dstoch,lrestart=.FALSE.,                       &
            & in_group=groups("iau_init_vars"))

    ELSE
       ALLOCATE(phy_stochconv%clmf_p(0,kblks),phy_stochconv%clmf_a(0,kblks),phy_stochconv%clnum_p(0,kblks),phy_stochconv%clnum_a(0,kblks))
       ALLOCATE(phy_stochconv%time_i(0,0,kblks),phy_stochconv%mf_i(0,0,kblks),phy_stochconv%life_i(0,0,kblks),phy_stochconv%area_i(0,0,kblks))
       ALLOCATE(phy_stochconv%depth_i(0,0,kblks),phy_stochconv%type_i(0,0,kblks),phy_stochconv%ktype_i(0,0,kblks),phy_stochconv%used_cell(0,0,kblks))
    ENDIF

    IF (atm_phy_nwp_config(k_jg)%lstoch_deep) THEN
       ! stochastic deep convection
       ! number of deep clouds
       cf_desc    = t_cf_var('clnum_d', 'm-2', 'number of deep clouds per unit area', datatype_flt32)
       ! NOTE: "highjacking" grib codes for shallow active clouds here, no grib codes defined yet for deep
       ! should not conflict, as deep/shallow stochastic schemes should not be used together
       grib2_desc = grib2_var(0, 6, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clnum_d', phy_stochconv%clnum_d,  &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,                                              &
            & loutput=.TRUE.,lrestart=.TRUE.,                             &
            & in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","iau_restore_vars"))
        
       ! perturbed cloud-base mass flux for deep clouds
       cf_desc    = t_cf_var('clmf_d', 'kg s-1 m-2', 'cloud base mass flux associated with deep clouds', datatype_flt32)
       ! NOTE: "highjacking" grib codes for shallow active clouds here, no grib codes defined yet for deep
       ! should not conflict, as deep/shallow stochastic schemes should not be used together      
       grib2_desc = grib2_var(0, 3, 200, ibits, GRID_UNSTRUCTURED, GRID_CELL)
       CALL add_var( stochconv_list,                                           &
            & 'clmf_d', phy_stochconv%clmf_d,    &
            & GRID_UNSTRUCTURED_CELL, ZA_CLOUD_BASE, cf_desc, grib2_desc, &
            & ldims=shape2d,                                              &
            & loutput=.TRUE.,lrestart=.TRUE.,                             &
            & in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","iau_restore_vars"))
    ELSE
       ALLOCATE(phy_stochconv%clmf_d(0,kblks),phy_stochconv%clnum_d(0,kblks))
    ENDIF
    

END SUBROUTINE new_nwp_phy_stochconv_list

SUBROUTINE new_nwp_phy_tend_list( k_jg, klev,  kblks,   &
                     & listname, phy_tend_list, phy_tend)

    INTEGER,INTENT(IN) :: k_jg, klev, kblks !< dimension sizes

    CHARACTER(len=*),INTENT(IN) :: listname

    TYPE(t_var_list_ptr)    ,INTENT(INOUT) :: phy_tend_list
    TYPE(t_nwp_phy_tend),INTENT(INOUT) :: phy_tend

    ! Local variables

    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d(3), shape3dkp1(3), shape4d(4), shape4d_conv(4), shape4d_gscp(4)
    INTEGER :: ibits, ktracer, ist, ntr_conv
    LOGICAL :: lrestart
    INTEGER :: datatype_flt
    INTEGER :: ncomin_tendphy_turb, ncomin_tendphy_conv

    ncomin_tendphy_turb = comin_config%comin_icon_domain_config(k_jg)%nturb_tracer
    ncomin_tendphy_conv = comin_config%comin_icon_domain_config(k_jg)%nconv_tracer

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    shape3d    = (/nproma, klev  , kblks            /)
    shape3dkp1 = (/nproma, klev+1, kblks            /)

    IF (lart) THEN
     shape4d    = (/nproma, klev  , kblks, nqtendphy+nart_tendphy+ncomin_tendphy_turb /)
    ELSE
     shape4d    = (/nproma, klev  , kblks, nqtendphy+ncomin_tendphy_turb /)
    ENDIF 
      
    ! dimension of convective tracer field
    ntr_conv = nqtendphy + ncomin_tendphy_conv
    IF (lart)                                        ntr_conv = ntr_conv + nart_tendphy
    IF (atm_phy_nwp_config(k_jg)%ldetrain_conv_prec) ntr_conv = ntr_conv + 2 ! plus qr and qs

    shape4d_conv = (/nproma, klev  , kblks, ntr_conv /)
    shape4d_gscp = (/nproma, klev  , kblks, 6 /)

    NULLIFY( &
      & phy_tend%ddt_temp_drag, &
      & phy_tend%ddt_temp_gscp, &
      & phy_tend%ddt_temp_pconv, &
      & phy_tend%ddt_tke, &
      & phy_tend%ddt_tke_hsh, &
      & phy_tend%ddt_tke_pconv, &
      & phy_tend%ddt_tracer_gscp, &
      & phy_tend%ddt_tracer_pconv, &
      & phy_tend%ddt_u_gwd, &
      & phy_tend%ddt_u_pconv, &
      & phy_tend%ddt_u_sso, &
      & phy_tend%ddt_v_gwd, &
      & phy_tend%ddt_v_pconv, &
      & phy_tend%ddt_v_sso, &
      & phy_tend%ddt_w_turb )

    CALL vlr_add(phy_tend_list, TRIM(listname), patch_id=k_jg ,lrestart=.TRUE.)
    
    !------------------------------
    ! Temperature tendencies
    !------------------------------

   ! &      phy_tend%ddt_temp_radsw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radsw', 'K s-1', &
         &                            'short wave radiative temperature tendency', datatype_flt)
    grib2_desc = grib2_var(0, 4, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radsw', phy_tend%ddt_temp_radsw,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3d,                                                  &
                & in_group=groups("phys_tendencies","iau_init_vars"),             &
                & lopenacc=.TRUE. )
                

   ! &      phy_tend%ddt_temp_radlw(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_radlw', 'K s-1', &
         &                            'long wave radiative temperature tendency', datatype_flt)
    grib2_desc = grib2_var(0, 5, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_radlw', phy_tend%ddt_temp_radlw,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3d,                                                  &
                & in_group=groups("phys_tendencies","iau_init_vars"),             &
                & lopenacc=.TRUE. )
                

    ! &      phy_tend%ddt_temp_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_turb', 'K s-1', &
         &                            'turbulence temperature tendency', datatype_flt)
    grib2_desc = grib2_var(192, 162, 121, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_turb', phy_tend%ddt_temp_turb,         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3d, lrestart=.FALSE.,                                &
                & in_group=groups("phys_tendencies","iau_init_vars"),             &
                & lopenacc=.TRUE. )
    

    ! &      phy_tend%ddt_temp_clcov(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_clcov', 'K s-1', &
         &                            'sgs condensation temperature tendency', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_clcov', phy_tend%ddt_temp_clcov,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & ldims=shape3d, lrestart=.TRUE.,                                 &
                & in_group=groups("phys_tendencies"), lopenacc=.TRUE. )
    

    cf_desc    = t_cf_var('ddt_temp_drag', 'K s-1', &
           &                'sso + gwdrag temperature tendency', datatype_flt)
    grib2_desc = grib2_var(192, 162, 125, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_drag', phy_tend%ddt_temp_drag,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, lrestart=.FALSE.,                              &
                  & in_group=groups("phys_tendencies"), lopenacc=.TRUE. )
                  

    ! &      phy_tend%ddt_temp_pconv(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_temp_pconv', 'K s-1', &
           &                            'convective temperature tendency', datatype_flt)
    grib2_desc = grib2_var(0, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_temp_pconv', phy_tend%ddt_temp_pconv,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, in_group=groups("phys_tendencies"),            &
                  & lopenacc=.TRUE. )
                  

    IF (atm_phy_nwp_config(k_jg)%is_les_phy .OR. &
        is_variable_in_output(var_name="ddt_temp_gscp")) THEN


      ! &      phy_tend%ddt_temp_gscp(nproma,nlev,nblks)
      cf_desc    = t_cf_var('ddt_temp_gscp', 'K s-1', &
           &                            'microphysical temperature tendency', datatype_flt)
      grib2_desc = grib2_var(192, 162, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_temp_gscp', phy_tend%ddt_temp_gscp,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, lrestart=.FALSE.,                              &
                  & in_group=groups("phys_tendencies"), lopenacc=.TRUE. )
      
    END IF

    !------------------------------
    ! Zonal Wind tendencies
    !------------------------------

   ! &      phy_tend%ddt_u_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_turb', 'm s-2', &
         &                            'turbulence tendency of zonal wind', datatype_flt)
    grib2_desc = grib2_var(192, 162, 119, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_turb', phy_tend%ddt_u_turb,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3d, lrestart=.FALSE.,                                &
                & in_group=groups("phys_tendencies"), lopenacc=.TRUE.)
    

    ! &      phy_tend%ddt_u_sso(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_sso', 'm s-2', &
            &                            'sso tendency of zonal wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 194, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_sso', phy_tend%ddt_u_sso,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                   & vert_interp=create_vert_interp_metadata(                     &
                   & vert_intp_type=vintp_types("P","Z","I"),                     &
                   & vert_intp_method=VINTP_METHOD_LIN),                          &
                   & ldims=shape3d, in_group=groups("phys_tendencies"),           &
                   & lopenacc=.TRUE. )
    
   
    ! &      phy_tend%ddt_u_gwd(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_gwd', 'm s-2', &
            &                            'GWD tendency of zonal wind', datatype_flt)
    grib2_desc = grib2_var(192, 128, 220, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_gwd', phy_tend%ddt_u_gwd,              &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                   & vert_interp=create_vert_interp_metadata(                     &
                   & vert_intp_type=vintp_types("P","Z","I"),                     &
                   & vert_intp_method=VINTP_METHOD_LIN),                          &
                   & ldims=shape3d, in_group=groups("phys_tendencies"),           &
                   & lopenacc=.TRUE. )
    
   
    ! &      phy_tend%ddt_u_pconv(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_u_pconv', 'm s-2', &
            &                            'convective tendency of zonal wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_u_pconv', phy_tend%ddt_u_pconv,          &
                   & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
                   & vert_interp=create_vert_interp_metadata(                     &
                   & vert_intp_type=vintp_types("P","Z","I"),                     &
                   & vert_intp_method=VINTP_METHOD_LIN),                          &
                   & ldims=shape3d, in_group=groups("phys_tendencies"),           &
                   & lopenacc=.TRUE. )
       


    !------------------------------
    ! Meridional Wind tendencies
    !------------------------------

    ! &      phy_tend%ddt_v_turb(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_turb', 'm s-2', &
         &                            'turbulence tendency of meridional wind', datatype_flt)
    grib2_desc = grib2_var(192, 162, 120, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_turb', phy_tend%ddt_v_turb,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3d, lrestart=.FALSE.,                                &
                & in_group=groups("phys_tendencies"), lopenacc=.TRUE.)
    

    ! &      phy_tend%ddt_v_sso(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_sso', 'm s-2', &
           &                            'sso tendency of meridional wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_sso', phy_tend%ddt_v_sso,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, in_group=groups("phys_tendencies"),            &
                  & lopenacc=.TRUE. )
    
  
    ! &      phy_tend%ddt_v_gwd(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_gwd', 'm s-2', &
           &                            'GWD tendency of meridional wind', datatype_flt)
    grib2_desc = grib2_var(192, 128, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_gwd', phy_tend%ddt_v_gwd,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, in_group=groups("phys_tendencies"),            &
                  & lopenacc=.TRUE. )
    
  
    ! &      phy_tend%ddt_v_pconv(nproma,nlev,nblks)
    cf_desc    = t_cf_var('ddt_v_pconv', 'm s-2', &
           &                            'convective tendency of meridional wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 193, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_v_pconv', phy_tend%ddt_v_pconv,           &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, in_group=groups("phys_tendencies"),            &
                  & lopenacc=.TRUE. )
    


    !------------------------------
    ! Vertical Wind tendencies
    !------------------------------

    IF ( atm_phy_nwp_config(k_jg)%is_les_phy ) THEN
      ! &      phy_tend%ddt_w_turb(nproma,nlev+1,nblks)
      cf_desc    = t_cf_var('ddt_w_turb', 'm s-2', &
           &                            'turbulence tendency of vertical wind', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_w_turb', phy_tend%ddt_w_turb,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & vert_interp=create_vert_interp_metadata(                      &
                  & vert_intp_type=vintp_types("P","Z","I"),                      &
                  & vert_intp_method=VINTP_METHOD_LIN),                           &
                  & ldims=shape3d, lrestart=.FALSE.,                              &
                  & initval=0._wp, in_group=groups("phys_tendencies"),            &
                  & lopenacc=.TRUE. )
      
    END IF



    !------------------------------
    ! Moist tracer tendencies
    !------------------------------

    ! --- Turbulence moist tracer tendencies

    ! &      phy_tend%ddt_tracer_turb(nproma,nlev,nblks,nqtendphy),          &
    cf_desc    = t_cf_var('ddt_tracer_turb', 's-1', &
         &                            'turbulence tendency of tracers', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_turb', phy_tend%ddt_tracer_turb,        &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape4d,&
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., &
                & in_group=groups("iau_init_vars"),                     &
                & lopenacc=.TRUE.)
    

    IF (lart) THEN
     ktracer=nqtendphy+nart_tendphy+ncomin_tendphy_turb
    ELSE
     ktracer=nqtendphy+ncomin_tendphy_turb
    ENDIF
    ALLOCATE( phy_tend%tracer_turb_ptr(ktracer) )
    !$ACC ENTER DATA CREATE(phy_tend%tracer_turb_ptr)

    !qv
    CALL add_ref( phy_tend_list, 'ddt_tracer_turb',                               &
                & 'ddt_qv_turb', phy_tend%tracer_turb_ptr(iqv)%p_3d,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('ddt_qv_turb', 'kg kg**-1 s**-1',                      &
                & 'turbulence tendency of specific humidity', datatype_flt),      &
                & grib2_var(192, 162, 122, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ref_idx=iqv, ldims=shape3d, lrestart=.FALSE.,                   &
                & in_group=groups("phys_tendencies") )
    !qc
    CALL add_ref( phy_tend_list, 'ddt_tracer_turb',                               &
                & 'ddt_qc_turb', phy_tend%tracer_turb_ptr(iqc)%p_3d,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('ddt_qc_turb', 'kg kg**-1 s**-1',                      &
                & 'turbulence tendency of specific cloud water', datatype_flt),   &
                & grib2_var(192, 162, 201, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ref_idx=iqc, ldims=shape3d, lrestart=.FALSE.,                   &
                & in_group=groups("phys_tendencies") )
    !qi
    CALL add_ref( phy_tend_list, 'ddt_tracer_turb',                               &
                & 'ddt_qi_turb', phy_tend%tracer_turb_ptr(iqi)%p_3d,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                           &
                & t_cf_var('ddt_qi_turb', 'kg kg**-1 s**-1',                      &
                & 'turbulence tendency of specific cloud ice', datatype_flt),     &
                & grib2_var(192, 162, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ref_idx=iqi, ldims=shape3d, lrestart=.FALSE.,                   &
                & in_group=groups("phys_tendencies") )


    ! --- Convection moist tracer tendencies

    cf_desc    = t_cf_var('ddt_tracer_pconv', 'kg m-3 s-1', &
           &                            'convective tendency of tracers', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tracer_pconv', phy_tend%ddt_tracer_pconv,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape4d_conv,&
                  & initval=0._wp, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., &
                  & lopenacc=.TRUE.)
    

    IF (lart) THEN
      ktracer=nqtendphy+nart_tendphy+ncomin_tendphy_conv
    ELSE
      ktracer=nqtendphy+ncomin_tendphy_conv
    ENDIF
    IF (atm_phy_nwp_config(k_jg)%ldetrain_conv_prec) ktracer = ktracer+2

    ALLOCATE( phy_tend%tracer_conv_ptr(ktracer) )
    !$ACC ENTER DATA CREATE(phy_tend%tracer_conv_ptr)

    !qv
    CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                                &
                  & 'ddt_qv_conv', phy_tend%tracer_conv_ptr(iqv)%p_3d,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  & t_cf_var('ddt_qv_conv', 'kg m-3 s-1',                             &
                  & 'convective tendency of absolute humidity', datatype_flt),        &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & vert_interp=create_vert_interp_metadata(                          &
                  & vert_intp_type=vintp_types("P","Z","I"),                          &
                  & vert_intp_method=VINTP_METHOD_LIN),                               &
                  & ref_idx=iqv, ldims=shape3d, in_group=groups("phys_tendencies") )
    !qc
    CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                                &
                  & 'ddt_qc_conv', phy_tend%tracer_conv_ptr(iqc)%p_3d,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  & t_cf_var('ddt_qc_conv', 'kg m-3 s-1',                             &
                  & 'convective tendency of cloud water mass density', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & vert_interp=create_vert_interp_metadata(                          &
                  & vert_intp_type=vintp_types("P","Z","I"),                          &
                  & vert_intp_method=VINTP_METHOD_LIN),                               &
                  & ref_idx=iqc, ldims=shape3d, in_group=groups("phys_tendencies") )
    !qi
    CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                                &
                  & 'ddt_qi_conv', phy_tend%tracer_conv_ptr(iqi)%p_3d,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  & t_cf_var('ddt_qi_conv', 'kg m-3 s-1',                             &
                  & 'convective tendency of cloud ice mass density', datatype_flt),   &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & vert_interp=create_vert_interp_metadata(                          &
                  & vert_intp_type=vintp_types("P","Z","I"),                          &
                  & vert_intp_method=VINTP_METHOD_LIN),                               &
                  & ref_idx=iqi, ldims=shape3d, in_group=groups("phys_tendencies") )

    IF (atm_phy_nwp_config(k_jg)%ldetrain_conv_prec) THEN
      !qr
      CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                              &
                  & 'ddt_qr_conv', phy_tend%tracer_conv_ptr(iqr)%p_3d,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  & t_cf_var('ddt_qr_conv', 'kg m-3 s-1',                             &
                  & 'convective tendency of rain mass density', datatype_flt),        &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & vert_interp=create_vert_interp_metadata(                          &
                  & vert_intp_type=vintp_types("P","Z","I"),                          &
                  & vert_intp_method=VINTP_METHOD_LIN),                               &
                  & ref_idx=iqr, ldims=shape3d, in_group=groups("phys_tendencies") )
      !qs
      CALL add_ref( phy_tend_list, 'ddt_tracer_pconv',                              &
                  & 'ddt_qs_conv', phy_tend%tracer_conv_ptr(iqs)%p_3d,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
                  & t_cf_var('ddt_qs_conv', 'kg m-3 s-1',                             &
                  & 'convective tendency of snow mass density', datatype_flt),        &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
                  & vert_interp=create_vert_interp_metadata(                          &
                  & vert_intp_type=vintp_types("P","Z","I"),                          &
                  & vert_intp_method=VINTP_METHOD_LIN),                               &
                  & ref_idx=iqs, ldims=shape3d, in_group=groups("phys_tendencies") )
    ENDIF


    ! --- Microphysics moist tracer tendencies

    IF ( is_variable_in_output(var_name="ddt_qv_gscp") .OR.   &
      &  is_variable_in_output(var_name="ddt_qc_gscp") .OR.   &
      &  is_variable_in_output(var_name="ddt_qi_gscp") .OR.   &
      &  is_variable_in_output(var_name="ddt_qr_gscp") .OR.   &
      &  is_variable_in_output(var_name="ddt_qs_gscp") .OR.   &
      &  is_variable_in_output(var_name="ddt_qg_gscp") ) THEN
      cf_desc    = t_cf_var('ddt_tracer_gscp', 's-1', &
           &                            'microphysics tendency of tracers', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( phy_tend_list, 'ddt_tracer_gscp', phy_tend%ddt_tracer_gscp,              &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, ldims=shape4d_gscp,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., lopenacc=.TRUE.)
      

      IF (lart) THEN
       ktracer=6 + nart_tendphy 
      ELSE
       ktracer=6
      ENDIF
      ALLOCATE( phy_tend%tracer_gscp_ptr(ktracer) )

      !qv
      IF ( iqv /= 0 ) THEN
        CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                    & 'ddt_qv_gscp', phy_tend%tracer_gscp_ptr(iqv)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var('ddt_qv_gscp', 'kg kg**-1 s**-1',                     &
                    & 'microphysics tendency of specific humidity', datatype_flt),   &
                    & grib2_var(192, 162, 204, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & vert_interp=create_vert_interp_metadata(                       &
                    & vert_intp_type=vintp_types("P","Z","I"),                       &
                    & vert_intp_method=VINTP_METHOD_LIN),                            &
                    & ref_idx=iqv,                                                   &
                    & ldims=shape3d, in_group=groups("phys_tendencies") )
      END IF

      !qc
      IF ( iqc /= 0 ) THEN
        CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                    & 'ddt_qc_gscp', phy_tend%tracer_gscp_ptr(iqc)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var('ddt_qc_gscp', 'kg kg**-1 s**-1',                     &
                    & 'microphysics tendency of specific cloud water', datatype_flt),&
                    & grib2_var(192, 162, 205, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & vert_interp=create_vert_interp_metadata(                       &
                    & vert_intp_type=vintp_types("P","Z","I"),                       &
                    & vert_intp_method=VINTP_METHOD_LIN),                            &
                    & ref_idx=iqc,                                                   &
                    & ldims=shape3d, in_group=groups("phys_tendencies") )
       END IF

      !qi
      IF ( iqi /= 0 ) THEN
        CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                    & 'ddt_qi_gscp', phy_tend%tracer_gscp_ptr(iqi)%p_3d,             &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var('ddt_qi_gscp', 'kg kg**-1 s**-1',                     &
                    & 'microphysics tendency of specific cloud ice', datatype_flt),  &
                    & grib2_var(192, 162, 206, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & vert_interp=create_vert_interp_metadata(                       &
                    & vert_intp_type=vintp_types("P","Z","I"),                       &
                    & vert_intp_method=VINTP_METHOD_LIN),                            &
                    & ref_idx=iqi,                                                   &
                    & ldims=shape3d, in_group=groups("phys_tendencies") )
      END IF

      !qr
      IF ( iqr /= 0 ) THEN
        CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                  & 'ddt_qr_gscp', phy_tend%tracer_gscp_ptr(iqr)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                  & t_cf_var('ddt_qr_gscp', 'kg kg**-1 s**-1',                       &
                  & 'microphysics tendency of rain', datatype_flt),                  &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                  & vert_interp=create_vert_interp_metadata(                         &
                  & vert_intp_type=vintp_types("P","Z","I"),                         &
                  & vert_intp_method=VINTP_METHOD_LIN),                              &
                  & ref_idx=iqr,                                                     &
                  & ldims=shape3d, in_group=groups("phys_tendencies") )
      END IF

      !qs
      IF ( iqs /= 0 ) THEN
        CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                  & 'ddt_qs_gscp', phy_tend%tracer_gscp_ptr(iqs)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                  & t_cf_var('ddt_qs_gscp', 'kg kg**-1 s**-1',                       &
                  & 'microphysics tendency of snow', datatype_flt),                  &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                  & vert_interp=create_vert_interp_metadata(                         &
                  & vert_intp_type=vintp_types("P","Z","I"),                         &
                  & vert_intp_method=VINTP_METHOD_LIN),                              &
                  & ref_idx=iqs,                                                     &
                  & ldims=shape3d, in_group=groups("phys_tendencies") )
      END IF

      !qg
      IF ( iqg /= 0 ) THEN
        CALL add_ref( phy_tend_list, 'ddt_tracer_gscp',                              &
                  & 'ddt_qg_gscp', phy_tend%tracer_gscp_ptr(iqg)%p_3d,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                  & t_cf_var('ddt_qg_gscp', 'kg kg**-1 s**-1',                       &
                  & 'microphysics tendency of graupel', datatype_flt),               &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
                  & vert_interp=create_vert_interp_metadata(                         &
                  & vert_intp_type=vintp_types("P","Z","I"),                         &
                  & vert_intp_method=VINTP_METHOD_LIN),                              &
                  & ref_idx=iqg,                                                     &
                  & ldims=shape3d, in_group=groups("phys_tendencies") )
      END IF

    END IF


    !------------------------------
    ! TKE tendency
    !------------------------------

    !      phy_tend%ddt_tke(nproma,nlevp1,nblks)
    cf_desc    = t_cf_var('ddt_tke', 'm s-2'          , &
           &                'tendency of turbulent velocity scale', datatype_flt)
    grib2_desc = grib2_var(0, 19, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke', phy_tend%ddt_tke,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3dkp1, lopenacc=.TRUE.,                              &
                & in_group=groups("phys_tendencies") )
    
  
    IF (ltkecon) THEN
      lrestart = .TRUE.
    ELSE
      lrestart = .FALSE.
    ENDIF
  
    !      phy_tend%ddt_tke_pconv(nproma,nlevp1,nblks)
    cf_desc    = t_cf_var('ddt_tke_pconv', 'm**2 s**-3'          , &
         &                'TKE tendency due to sub-grid scale convection', datatype_flt)
    grib2_desc = grib2_var(0, 19, 219, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke_pconv', phy_tend%ddt_tke_pconv,       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3dkp1, lrestart=lrestart, lopenacc=.TRUE.,           &
                & in_group=groups("phys_tendencies") )  
    


    !      phy_tend%ddt_tke_hsh(nproma,nlevp1,nblks)
    cf_desc    = t_cf_var('ddt_tke_hsh', 'm**2 s**-3'          , &
         &                'TKE tendency horizonzal shear production', datatype_flt)
    grib2_desc = grib2_var(0, 19, 221, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( phy_tend_list, 'ddt_tke_hsh', phy_tend%ddt_tke_hsh,           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
                & vert_interp=create_vert_interp_metadata(                        &
                & vert_intp_type=vintp_types("P","Z","I"),                        &
                & vert_intp_method=VINTP_METHOD_LIN),                             &
                & ldims=shape3dkp1, lrestart=.FALSE., lopenacc=.TRUE.,            &
                & in_group=groups("phys_tendencies") )  
    

 
   !Anurag Dipankar, MPIM (2013-May-31)
   !Large-scale tendencies for idealized testcases: add_var doesn't work
   !for 1D variables so using ALLOCATE-DEALLOCATE
   !Therefore, these variables can't go into restart/output
   !Initialize them all to 0 
    IF(is_ls_forcing)THEN

      ALLOCATE(phy_tend%ddt_u_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%u_ls failed')
      ENDIF
      phy_tend%ddt_u_ls = 0._wp

      ALLOCATE(phy_tend%ddt_v_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%v_ls failed')
      ENDIF
      phy_tend%ddt_v_ls = 0._wp

      ALLOCATE(phy_tend%ddt_temp_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%temp_ls failed')
      ENDIF
      phy_tend%ddt_temp_ls = 0._wp

      ALLOCATE(phy_tend%ddt_tracer_ls(klev,nqtendphy),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%tracer_ls failed')
      ENDIF
      phy_tend%ddt_tracer_ls = 0._wp

      ! Added by Christopher Moseley: 7 output variables for LS tendencies

      ALLOCATE(phy_tend%ddt_temp_subs_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_subs_ls failed')
      ENDIF
      phy_tend%ddt_temp_subs_ls = 0._wp

      ALLOCATE(phy_tend%ddt_qv_subs_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_qv_subs_ls failed')
      ENDIF
      phy_tend%ddt_qv_subs_ls = 0._wp

      ALLOCATE(phy_tend%ddt_temp_adv_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_adv_ls failed')
      ENDIF
      phy_tend%ddt_temp_adv_ls = 0._wp

      ALLOCATE(phy_tend%ddt_qv_adv_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_qv_adv_ls failed')
      ENDIF
      phy_tend%ddt_qv_adv_ls = 0._wp

      ALLOCATE(phy_tend%ddt_u_adv_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_u_adv_ls failed')
      ENDIF
      phy_tend%ddt_u_adv_ls = 0._wp

      ALLOCATE(phy_tend%ddt_v_adv_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_v_adv_ls failed')
      ENDIF
      phy_tend%ddt_v_adv_ls = 0._wp

      ALLOCATE(phy_tend%ddt_temp_nud_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_nud_ls failed')
      ENDIF
      phy_tend%ddt_temp_nud_ls = 0._wp

      ALLOCATE(phy_tend%ddt_qv_nud_ls(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_qv_nud_ls failed')
      ENDIF
      phy_tend%ddt_qv_nud_ls = 0._wp

      ALLOCATE(phy_tend%wsub(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%wsub failed')
      ENDIF
      phy_tend%wsub = 0._wp

      ALLOCATE(phy_tend%temp_nudge(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%temp_nudge failed')
      ENDIF
      phy_tend%temp_nudge = 0._wp

      ALLOCATE(phy_tend%q_nudge(klev,3),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%q_nudge failed')
      ENDIF
      phy_tend%q_nudge = 0._wp

      ALLOCATE(phy_tend%u_nudge(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%u_nudge failed')
      ENDIF
      phy_tend%u_nudge = 0._wp

      ALLOCATE(phy_tend%v_nudge(klev),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%v_nudge failed')
      ENDIF
      phy_tend%v_nudge = 0._wp

      ALLOCATE(phy_tend%ddt_temp_sim_rad(nproma,klev,kblks),STAT=ist)
      IF (ist/=SUCCESS)THEN
        CALL finish('mo_nwp_phy_state:construct_nwp_phy_tend', &
                    'allocation for phy_tend%ddt_temp_sim_rad failed')
      ENDIF
      phy_tend%ddt_temp_sim_rad = 0._wp

    END IF

    CALL message('mo_nwp_phy_state:construct_nwp_phy_tend', &
      'construction of NWP physical tendency fields finished')

END SUBROUTINE new_nwp_phy_tend_list
!
!-------------------------------------------------------------------------

END MODULE mo_nwp_phy_state
!<
