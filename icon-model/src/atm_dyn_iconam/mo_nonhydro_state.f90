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

! Constructs and destructs the state vector of the nonhydrostatic.
!
! Constructs and destructs the state vector of the nonhydrostatic
! model variables. They are subdivided in several classes: prognostics
! and diagnostics.

!NEC$ options "-O0"

MODULE mo_nonhydro_state

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS, vname_len, vlname_len,                 &
    &                                INWP, iaes, n_camsaermr,                        &
    &                                VINTP_METHOD_VN,                                &
    &                                VINTP_METHOD_QV, VINTP_METHOD_PRES,             &
    &                                VINTP_METHOD_LIN,                               &
    &                                VINTP_METHOD_LIN_NLEVP1,                        &
    &                                TASK_INTP_MSL, HINTP_TYPE_NONE,                 &
    &                                MODE_IAU, MODE_IAU_OLD,                         &
    &                                TASK_COMPUTE_OMEGA, TLEV_NNOW_RCF,              &
    &                                MODE_ICONVREMAP,HINTP_TYPE_LONLAT_RBF,          &
    &                                HINTP_TYPE_LONLAT_BCTR,                         &
    &                                MIN_RLCELL_INT, MIN_RLCELL
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
    &                                GRID_UNSTRUCTURED_VERT, GRID_CELL, GRID_EDGE,   &
    &                                GRID_VERTEX
  USE mo_exception,            ONLY: message, finish
  USE mo_master_control,       ONLY: get_my_process_name
  USE mo_model_domain,         ONLY: t_patch
  USE mo_nonhydro_types,       ONLY: t_nh_state, t_nh_state_lists,       &
                                     t_nh_prog, t_nh_diag,               &
    &                                t_nh_ref, t_nh_metrics
  USE mo_grid_config,          ONLY: n_dom, l_limited_area, ifeedback_type
  USE mo_nonhydrostatic_config,ONLY: itime_scheme, igradp_method, ndyn_substeps_max, &
    &                                lcalc_dpsdt
  USE mo_dynamics_config,      ONLY: nsav1, nsav2
  USE mo_parallel_config,      ONLY: nproma
  USE mo_run_config,           ONLY: iforcing, ntracer, iqm_max, iqt,           &
    &                                iqv, iqc, iqi, iqr, iqs,                   &
    &                                ico2, ich4, in2o, io3,                     &
    &                                iqni, iqg, iqh, iqnr, iqns,                & 
    &                                iqng, iqnh, iqnc, inccn, ininpot, ininact, &
    &                                iqgl, iqhl,                                &
    &                                iqtke, ltestcase, lart,                    &
    &                                iqbin, iqb_i, iqb_e, iqb_s            
  USE mo_coupling_config,      ONLY: is_coupled_to_ocean
  USE mo_radiation_config,     ONLY: irad_aero, iRadAeroCAMSclim, iRadAeroCAMStd
  USE mo_io_config,            ONLY: inextra_2d, inextra_3d, lnetcdf_flt64_output, &
    &                                t_var_in_output
  USE mo_limarea_config,       ONLY: latbc_config
  USE mo_advection_config,     ONLY: t_advection_config, advection_config
  USE mo_turbdiff_config,      ONLY: turbdiff_config
  USE mo_diffusion_config,     ONLY: diffusion_config
  USE mo_initicon_config,      ONLY: init_mode, qcana_mode, qiana_mode, qrsgana_mode, &
    &                                icpl_da_sfcevap, icpl_da_skinc, icpl_da_sfcfric
  USE mo_nudging_config,       ONLY: nudging_config, indg_type
  USE mo_var_list,             ONLY: add_var, find_list_element, add_ref, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_var,                  ONLY: t_var
  USE mo_var_groups,           ONLY: MAX_GROUPS, groups
  USE mo_var_metadata,         ONLY: create_vert_interp_metadata,            &
    &                                create_hor_interp_metadata,             &
    &                                vintp_types, get_timelevel_string
  USE mo_tracer_metadata,      ONLY: create_tracer_metadata,                 &
    &                                create_tracer_metadata_hydro
  USE mo_advection_utils,      ONLY: add_tracer_ref
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+)
  USE mo_gribout_config,       ONLY: gribout_config
#ifdef __ICON_ART
  USE mo_art_tracer_interface, ONLY: art_tracer_interface
  USE mo_art_diagnostics_interface, ONLY: art_diagnostics_interface_init
#endif
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_zaxis_type,           ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF, ZA_HEIGHT_2M, ZA_HEIGHT_10M, &
    &                                ZA_REFERENCE_HALF_HHL, ZA_SURFACE, ZA_MEANSEA
  USE mo_cdi,                  ONLY: DATATYPE_FLT32, DATATYPE_FLT64,                 &
    &                                DATATYPE_PACK16, DATATYPE_PACK24,               &
    &                                DATATYPE_INT, TSTEP_CONSTANT, GRID_UNSTRUCTURED
  USE mo_aes_vdf_config,       ONLY: aes_vdf_config
  USE mo_turb_vdiff_params,    ONLY: VDIFF_TURB_3DSMAGORINSKY
  USE mo_loopindices,          ONLY: get_indices_c

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nonhydro_state'


  PUBLIC :: construct_nh_state    ! Constructor for the nonhydrostatic state
  PUBLIC :: destruct_nh_state     ! Destructor for the nonhydrostatic state
  PUBLIC :: duplicate_prog_state  ! Copy the prognostic state
  PUBLIC :: new_zd_metrics
  PUBLIC :: p_nh_state            ! state vector of nonhydrostatic variables (variable)
  PUBLIC :: p_nh_state_lists      ! lists for state vector of nonhydrostatic variables (variable)
  

  TYPE(t_nh_state),       TARGET, ALLOCATABLE :: p_nh_state(:)
  TYPE(t_nh_state_lists), TARGET, ALLOCATABLE :: p_nh_state_lists(:)

  CONTAINS

!-------------------------------------------------------------------------
!!            SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
!-------------------------------------------------------------------------
!
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! Top-level procedure for building the prognostic and diagnostic states.
  !! It calls constructors to single time level prognostic states, and 
  !! diagnostic states.
  !! Initialization of all components with zero.
  !!
  SUBROUTINE construct_nh_state(p_patch,  n_timelevels, var_in_output)
    TYPE(t_patch),         INTENT(IN) ::      & ! patch
      &  p_patch(n_dom)
    INTEGER,               INTENT(IN) ::      & ! number of timelevels
      &  n_timelevels
    TYPE(t_var_in_output), INTENT(IN) ::      & !< switches for optional diagnostics
      &  var_in_output(:)

    INTEGER  :: ntl,      &! local number of timelevels (with extra timelevels)
                ist,      &! status
                jg,       &! grid level counter
                jt         ! time level counter
    LOGICAL  :: l_extra_timelev

    CHARACTER(len=vlname_len) :: listname
    CHARACTER(LEN=vname_len) :: varname_prefix

    CHARACTER(*), PARAMETER :: routine = modname//'::construct_nh_state'
!-----------------------------------------------------------------------

    CALL message (routine, 'Construction of NH state started')


    ALLOCATE(p_nh_state(n_dom), p_nh_state_lists(n_dom), stat=ist)
    IF (ist /= success) CALL finish(routine, &
      &  'allocation of nonhydrostatic state  array and list failed')

    !$ACC ENTER DATA COPYIN(p_nh_state)

    DO jg = 1, n_dom

      ntl = n_timelevels

      ! As grid nesting is not called at every dynamics time step, an extra time
      ! level is needed for full-field interpolation and boundary-tendency calculation
      IF (n_dom > 1) THEN
        ntl = ntl + 1
        nsav1(jg) = ntl
      ENDIF

      ! In the presence of grid nesting and incremental feedback, another extra time level is needed to
      ! compute the feedback increments
      ! This extra time level is also used to store the driving-model data in the
      ! limited-area mode
      IF (ifeedback_type == 1 .AND. jg > 1 .OR. l_limited_area .AND. jg == 1) ntl = ntl + 1
      nsav2(jg) = ntl

      !
      ! Allocate pointer array p_nh_state(jg)%prog, as well as the
      ! corresponding list array for each grid level.
      !
      ! create state arrays
      ALLOCATE(p_nh_state(jg)%prog(1:ntl), STAT=ist)
      IF (ist/=SUCCESS) CALL finish(routine,                                   &
        'allocation of prognostic state array failed')
      !$ACC ENTER DATA COPYIN(p_nh_state(jg)%prog, p_nh_state(jg)%diag, p_nh_state(jg)%metrics, p_nh_state(jg)%ref)

      ! create state list
      ALLOCATE(p_nh_state_lists(jg)%prog_list(1:ntl), STAT=ist)
      IF (ist/=SUCCESS) CALL finish(routine,                                   &
        'allocation of prognostic state list array failed')

      ! Create tracer list (no extra timelevels).
      ! The list will be filled in `configure_advection`, once the total number
      ! of tracers is known.
      ALLOCATE(p_nh_state_lists(jg)%tracer_list(1:n_timelevels), STAT=ist)
      IF (ist/=SUCCESS) CALL finish(routine,                                   &
        'allocation of prognostic tracer list array failed')


      !
      ! Build lists for every timelevel
      !
      DO jt = 1, ntl

        ! Tracer fields do not need extra time levels because feedback is not incremental
        ! and the nest-call frequency is always synchronized with the advection time step
        l_extra_timelev = jt > n_timelevels

        WRITE(listname,'(a,i2.2,a,i2.2)') 'nh_state_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        ! Build prog state list
        ! includes memory allocation
        !
        ! varname_prefix = 'nh_prog_'
        varname_prefix = ''
        CALL new_nh_state_prog_list(p_patch(jg), p_nh_state(jg)%prog(jt),  &
          &  p_nh_state_lists(jg)%prog_list(jt), listname, TRIM(varname_prefix), &
          &  l_extra_timelev, jt)

      ENDDO ! jt

      !
      ! Build diag state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_diag_of_domain_',jg
      CALL new_nh_state_diag_list(p_patch(jg), p_nh_state(jg)%diag, &
        &  p_nh_state_lists(jg)%diag_list, listname, var_in_output(jg))

#ifdef __ICON_ART
      ! art: add ART diagnostics to diag list
      IF (lart) THEN
        CALL art_diagnostics_interface_init(jg, p_nh_state_lists(jg)%diag_list, &
          & p_prog_list=p_nh_state_lists(jg)%prog_list(1))
      ENDIF
#endif

      !
      ! Build metrics state list
      ! includes memory allocation
      !
      WRITE(listname,'(a,i2.2)') 'nh_state_metrics_of_domain_',jg
      CALL new_nh_metrics_list(p_patch(jg), p_nh_state(jg)%metrics, &
        &  p_nh_state_lists(jg)%metrics_list, listname )

      !
      ! Build ref state list
      ! includes memory allocation
      !
      IF ( ltestcase ) THEN      ! It is much cleaner to test this here rather than in new_nh_state_ref_list
        WRITE(listname,'(a,i2.2)') 'nh_state_ref_of_domain_',jg
        CALL new_nh_state_ref_list(p_patch(jg), p_nh_state(jg)%ref, &
          &  p_nh_state_lists(jg)%ref_list, listname)
      ELSE
        NULLIFY(p_nh_state(jg)%ref%vn_ref, p_nh_state(jg)%ref%w_ref)   ! See new_nh_state_ref_list
      ENDIF

    ENDDO ! jg

    CALL message (routine, 'NH state construction completed')

  END SUBROUTINE construct_nh_state


!-------------------------------------------------------------------------
!
!
  !>
  !! Destructor for prognostic and diagnostic states.
  !!
  !! It calls destructors to
  !! single time level prognostic states, and diagnostic states.
  !!
  SUBROUTINE destruct_nh_state(p_nh_state, p_nh_state_lists)
!
    TYPE(t_nh_state),       ALLOCATABLE, INTENT(INOUT) :: & ! nh state at different grid levels
      &  p_nh_state(:)

    TYPE(t_nh_state_lists), ALLOCATABLE, INTENT(INOUT) :: & ! lists of nh state at different grid levels
      &  p_nh_state_lists(:)

    INTEGER  :: ntl_prog, & ! number of timelevels prog state
                ntl_tra,  & ! number of timelevels 
                ist, &      ! status
                jg,  &      ! grid level counter
                jt          ! time level counter

    CHARACTER(len=*), PARAMETER :: routine = modname//'::destruct_nh_state'
!-----------------------------------------------------------------------

    CALL message(routine, 'Destruction of NH state started')


    !$ACC WAIT(1)
    DO jg = 1, n_dom

      ntl_prog = SIZE(p_nh_state(jg)%prog(:))
      IF(ntl_prog==0)THEN
        CALL finish(routine, 'prognostic array has no timelevels')
      ENDIF

      ntl_tra = SIZE(p_nh_state_lists(jg)%tracer_list(:))
      IF(ntl_tra==0)THEN
        CALL finish(routine, 'tracer list has no timelevels')
      ENDIF

      ! delete reference state list elements
      IF ( ltestcase ) THEN
        CALL vlr_del(p_nh_state_lists(jg)%ref_list)
      ENDIF

      ! delete diagnostic state list elements
      CALL vlr_del(p_nh_state_lists(jg)%diag_list)

      ! delete metrics state list elements
      CALL vlr_del(p_nh_state_lists(jg)%metrics_list)


      ! delete prognostic state list elements
      DO jt = 1, ntl_prog
        CALL vlr_del(p_nh_state_lists(jg)%prog_list(jt))
      ENDDO

      ! delete tracer list list elements
      DO jt = 1, ntl_tra
        CALL vlr_del(p_nh_state_lists(jg)%tracer_list(jt))
      ENDDO

      !$ACC EXIT DATA DELETE(p_nh_state(jg)%prog, p_nh_state(jg)%metrics, p_nh_state(jg)%ref, p_nh_state(jg)%diag)

      ! destruct state lists and arrays
      DEALLOCATE(p_nh_state_lists(jg)%prog_list, STAT=ist )
      IF (ist/=SUCCESS) CALL finish(routine,&
        & 'deallocation of prognostic state list array failed')

      DEALLOCATE(p_nh_state_lists(jg)%tracer_list, STAT=ist )
      IF (ist/=SUCCESS) CALL finish(routine,&
        & 'deallocation of tracer list array failed')

      DEALLOCATE(p_nh_state(jg)%prog )
      IF (ist/=SUCCESS) CALL finish (routine,&
        & 'deallocation of prognostic state array failed')

    ENDDO

    !$ACC EXIT DATA DELETE(p_nh_state)

    DEALLOCATE (p_nh_state, p_nh_state_lists, STAT=ist)
    IF (ist /= SUCCESS) CALL finish(routine,'deallocation of nonhydrostatic state vector failed')

    CALL message(routine, 'NH state destruction completed')

  END SUBROUTINE destruct_nh_state
  !-------------------------------------------------------------------------
  !
  !
  !>
  !!
  !! duplicate prognostic state 
  !!
  SUBROUTINE duplicate_prog_state ( p_prog_i, p_prog_d)

      TYPE(t_nh_prog), INTENT(IN)      :: &  !< prognostic state vector to be copied
     &  p_prog_i
      TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< duplicated prognostic state vector
     &  p_prog_d

    !--------------------------------------------------------------
     p_prog_d%w(:,:,:)              = p_prog_i%w(:,:,:)
     p_prog_d%vn(:,:,:)             = p_prog_i%vn(:,:,:)
     p_prog_d%rho(:,:,:)            = p_prog_i%rho(:,:,:)
     p_prog_d%theta_v(:,:,:)        = p_prog_i%theta_v(:,:,:)
     IF (ASSOCIATED(p_prog_i%exner)) &
       p_prog_d%exner(:,:,:)          = p_prog_i%exner(:,:,:)
     IF (ASSOCIATED(p_prog_i%tracer)) THEN
      p_prog_d%tracer(:,:,:,:)       = p_prog_i%tracer(:,:,:,:)
     END IF
     IF (ASSOCIATED(p_prog_i%tke )) THEN
      p_prog_d%tke(:,:,:)            = p_prog_i%tke(:,:,:)
     END IF

  END SUBROUTINE duplicate_prog_state

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Allocation of components of prognostic state.
  !!
  !! Initialization of components with zero.
  !!
  SUBROUTINE new_nh_state_prog_list ( p_patch, p_prog, p_prog_list,  &
    &                                 listname, vname_prefix, l_extra_timelev, timelev)
!
    TYPE(t_patch), TARGET, INTENT(IN) :: & !< current patch
      &  p_patch

    TYPE(t_nh_prog),  INTENT(INOUT)   :: & !< current prognostic state
      &  p_prog 

    TYPE(t_var_list_ptr), INTENT(INOUT)   :: p_prog_list !< current prognostic state list

    CHARACTER(len=*), INTENT(IN)      :: & !< list name
      &  listname, vname_prefix

    LOGICAL, INTENT(IN) :: l_extra_timelev  !< specifies extra time levels for which 
                                            !< not all variables are allocated

    INTEGER, INTENT(IN) :: timelev


    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e       !< number of edge blocks to allocate

    INTEGER :: nlev, nlevp1

    INTEGER :: shape3d_c(3), shape3d_e(3), shape3d_chalf(3), &
      &        shape4d_c(4)

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for some thermodynamic fields

    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    CHARACTER(len=4) :: suffix
    CHARACTER(len=2) :: passive_tracer_suffix

    TYPE(t_advection_config), POINTER :: advconf

    INTEGER           :: jt
    INTEGER           :: ipassive        ! loop counter
    INTEGER           :: dummy_idx, vntl, tlen

    CHARACTER(LEN=vname_len)             :: tracer_container_name
    CHARACTER(LEN=vname_len+LEN(suffix)) :: tracer_name

    TYPE(t_var), POINTER :: target_element
    INTEGER              :: tracer_idx

    INTEGER                       :: iqb
    LOGICAL :: ingroup(MAX_GROUPS)
    !**
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! pointer to advection_config(jg) to save some paperwork
    advconf => advection_config(p_patch%id)

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! predefined array shapes
    shape3d_e     = (/nproma, nlev,   nblks_e  /)
    shape3d_c     = (/nproma, nlev,   nblks_c  /)
    shape3d_chalf = (/nproma, nlevp1, nblks_c  /)
    shape4d_c     = (/nproma, nlev,   nblks_c, ntracer /)

    ! Suffix (mandatory for time level dependent variables)

    suffix = get_timelevel_string(timelev)

    vntl = LEN_TRIM(vname_prefix)
    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_prog_list, TRIM(listname), patch_id=p_patch%id, &
      &          lrestart=.TRUE., model_type=get_my_process_name())

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(p_prog%w, p_prog%vn, p_prog%rho, p_prog%exner, p_prog%theta_v, p_prog%tracer, p_prog%tke)


    !------------------------------
    ! Meteorological quantities
    !------------------------------

    ! vn           p_prog%vn(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge', datatype_flt)
    grib2_desc = grib2_var(0, 2, 34, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_prog_list, vname_prefix(1:vntl)//'vn'//suffix, p_prog%vn,     &
      &           GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
      &           vert_interp=create_vert_interp_metadata(                      &
      &             vert_intp_method=VINTP_METHOD_VN,                           &
      &             l_hires_intp=.FALSE., l_restore_fricred=.FALSE.),           &
      &           ldims=shape3d_e,                                              &
      &           in_group=groups("nh_prog_vars","dwd_fg_atm_vars",             &
      &                           "mode_dwd_fg_in","mode_iau_fg_in",            &
      &                           "mode_iau_old_fg_in","LATBC_PREFETCH_VARS",   &
      &                           "iau_restore_vars"),                          &
      &           lopenacc = .TRUE. )
    __acc_attach(p_prog%vn)

    ! w            p_prog%w(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'Vertical velocity', datatype_flt)
    grib2_desc = grib2_var(0, 2, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_prog_list, vname_prefix(1:vntl)//'w'//suffix, p_prog%w,      &
      &          GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
      &          ldims=shape3d_chalf,                                          & 
      &          vert_interp=create_vert_interp_metadata(                      &
      &             vert_intp_type=vintp_types("P","Z","I"),                   &
      &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                &
      &          in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars", &
      &                          "dwd_fg_atm_vars","mode_dwd_fg_in",           &
      &                          "mode_iau_fg_in","mode_iau_old_fg_in",        &
      &                          "LATBC_PREFETCH_VARS",                        &
      &                          "mode_iniana","icon_lbc_vars",                &
      &                          "iau_restore_vars"),                          &
      &          lopenacc = .TRUE.)
    __acc_attach(p_prog%w)

    ! rho          p_prog%rho(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('air_density', 'kg m-3', 'density', datatype_flt)
    grib2_desc = grib2_var(0, 3, 10, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_prog_list, vname_prefix(1:vntl)//'rho'//suffix, p_prog%rho,  &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,   &
      &           ldims=shape3d_c,                                             &
      &           vert_interp=create_vert_interp_metadata(                     &
      &              vert_intp_type=vintp_types("P","Z","I"),                  &
      &              vert_intp_method=VINTP_METHOD_LIN ),                      &
      &           in_group=groups("nh_prog_vars","dwd_fg_atm_vars",            &
      &                           "mode_dwd_fg_in","mode_iau_fg_in",           &
      &                           "mode_iau_old_fg_in","LATBC_PREFETCH_VARS",  &
      &                           "iau_restore_vars"),                         &
      &           lopenacc = .TRUE. )
    __acc_attach(p_prog%rho)

    ! theta_v      p_prog%theta_v(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('virtual_potential_temperature', 'K', &
      &                   'virtual potential temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 15, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_prog_list, vname_prefix(1:vntl)//'theta_v'//suffix, p_prog%theta_v, &
      &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
      &           ldims=shape3d_c,                                                    &
      &           in_group=groups("nh_prog_vars","dwd_fg_atm_vars",                   &
      &           "mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in",             &
      &           "LATBC_PREFETCH_VARS",                                              &
      &           "iau_restore_vars"),                                                &
      &           lopenacc = .TRUE. )
    __acc_attach(p_prog%theta_v)


    IF (.NOT. l_extra_timelev) THEN
      ! exner        p_prog%exner(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('exner_pressure', '-', 'exner pressure', datatype_flt)
      grib2_desc = grib2_var(0, 3, 26, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_prog_list, vname_prefix(1:vntl)//'exner'//suffix, p_prog%exner, &
        &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
        &           ldims=shape3d_c,                                                  &
        &           in_group=groups("nh_prog_vars","iau_restore_vars"),               &
        &           lopenacc = .TRUE. )
      __acc_attach(p_prog%exner)

      ! Tracer array for (model) internal use

      ! tracer         p_prog%tracer(nproma,nlev,nblks_c,ntracer)
      IF (ntracer > 0) THEN
        tracer_container_name = vname_prefix(1:vntl)//'tracer'//suffix
        cf_desc    = t_cf_var('tracer', 'kg kg-1', 'tracer', datatype_flt)
        grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_prog_list, tracer_container_name, p_prog%tracer,          &
          &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
          &           ldims=shape4d_c, tlev_source=TLEV_NNOW_RCF,                 &
          &           lcontainer=.TRUE., lrestart=.FALSE.,                        &
          &           in_group=groups("iau_restore_vars"),                        &
          &           loutput=.FALSE.,                                            &
          &           lopenacc = .TRUE. )
        __acc_attach(p_prog%tracer)
      ENDIF

      ALLOCATE( p_prog%tracer_ptr(ntracer) )
      !$ACC ENTER DATA CREATE(p_prog%tracer_ptr)

      IF ( iforcing == inwp .OR. iforcing == iaes ) THEN
        
        ! References to individual tracers, for I/O and setting of additional metadata
        ! ----------------------------------------------------------------------------
        ! The section below creates references to individual tracer fields (qv, qc, qi, ...)
        ! in the "tracer" container. This is done for active tracers, i.e. for tracers
        ! with corresponding index variables (iqv, iqc, iqi, ...) set to a non-zero value
        ! in the model configuration. The setup of the active, i.e. non-zero tracer indices
        ! depends on the selected forcing and is currently handled in these subroutines:
        ! - for iforcing = inwp   : mo_nml_crosscheck:atm_crosscheck
        ! - for iforcing = iaes : mo_aes_phy_init:init_aes_phy_itracer
        !
        ! For additional tracers, which do not have specific index variables, the indices
        ! need to be set via add_tracer_ref. 
        !
        ! Note that we make use of 
        ! - create_tracer_metadata
        ! - create_tracer_metadata_hydroMass
        ! - create_tracer_metadata_hydroNr
        ! for creating tracer-specific metadata. As a side effect, tracers are added to 
        ! distinct tracer groups, depending on the create_tracer_metadata[...] routine 
        ! used. Thus, make sure to use the right one when adding additional tracers.
        ! create_tracer_metadata[...] are described in more detail in mo_tracer_metadata.

        !QV
        IF ( iqv /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqv))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqv)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(iqv)%p_3d,                      &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(vname_prefix(1:vntl)//'specific_humidity',            &
            &                    'kg kg-1','Specific humidity',datatype_flt),          &
            &           grib2_var( 0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=iqv,                                                   &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
            &                       name        = tracer_name,                         &
            &                       lfeedback   = .TRUE.,                              &
            &                       ihadv_tracer=advconf%ihadv_tracer(iqv),            &
            &                       ivadv_tracer=advconf%ivadv_tracer(iqv)),           &
            &           hor_interp=create_hor_interp_metadata(                         &
            &                       hor_intp_type=HINTP_TYPE_LONLAT_BCTR,              &
            &                       fallback_type=HINTP_TYPE_LONLAT_RBF),              &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_QV,                  &
            &                       l_satlimit=.FALSE.,                                &
            &                       lower_limit=2.5e-7_wp, l_restore_pbldev=.FALSE. ), &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
            &                           "dwd_fg_atm_vars","mode_dwd_ana_in",           &
            &                           "LATBC_PREFETCH_VARS","mode_iau_fg_in",        &
            &                           "mode_iau_old_fg_in","mode_iau_ana_in",        &
            &                           "mode_iau_anaatm_in",                          &
            &                           "mode_iau_old_ana_in",                         &
            &                           "mode_iniana","icon_lbc_vars") )
          __acc_attach(p_prog%tracer_ptr(iqv)%p_3d)
        END IF

        IF ( latbc_config%latbc_contains_qcqi ) THEN
          ingroup=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",               &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",          &
            &                           "mode_iau_ana_in", "mode_iau_anaatm_in",     &
            &                           "mode_iau_fg_in","mode_iau_old_fg_in",       &
            &                           "LATBC_PREFETCH_VARS",                       &
            &                           "mode_iniana","icon_lbc_vars")
        ELSE
          ingroup=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",               &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",          &
            &                           "mode_iau_ana_in", "mode_iau_anaatm_in",     &
            &                           "mode_iau_fg_in","mode_iau_old_fg_in",       &
            &                           "mode_iniana","icon_lbc_vars")
        ENDIF

        !QC
        IF ( iqc /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqc))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqc)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                          &
            &         tracer_name, p_prog%tracer_ptr(iqc)%p_3d,                      &
            &         GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &         t_cf_var(tracer_name(1:tlen+vntl),                             &
            &          'kg kg-1', 'specific cloud water content', datatype_flt),     &
            &         grib2_var(0, 1, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &         ref_idx=iqc,                                                   &
            &         ldims=shape3d_c,                                               &
            &         tlev_source=TLEV_NNOW_RCF,                                     &              ! output from nnow_rcf slice
            &         tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,    &
            &                     name        = tracer_name,                         &
            &                     lfeedback   = .TRUE.,                              &
            &                     ihadv_tracer=advconf%ihadv_tracer(iqc),            &
            &                     ivadv_tracer=advconf%ivadv_tracer(iqc)),           &
            &         vert_interp=create_vert_interp_metadata(                       &
            &                     vert_intp_type=vintp_types("P","Z","I"),           &
            &                     vert_intp_method=VINTP_METHOD_LIN,                 &
            &                     l_loglin=.FALSE.,                                  &
            &                     l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                     lower_limit=0._wp  ),                              &
            &         in_group=ingroup)
          __acc_attach(p_prog%tracer_ptr(iqc)%p_3d)
        END IF ! iqc

        !QI
        IF ( iqi /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqi))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqi)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                          &
            &         tracer_name, p_prog%tracer_ptr(iqi)%p_3d,                      &
            &         GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &         t_cf_var(tracer_name(1:vntl+tlen),                             &
            &          'kg kg-1','specific cloud ice content', datatype_flt),        &
            &         grib2_var(0, 1, 82, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &         ref_idx=iqi,                                                   &
            &         ldims=shape3d_c,                                               &
            &         tlev_source=TLEV_NNOW_RCF,                                     &              ! output from nnow_rcf slice
            &         tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,    &
            &                     name        = tracer_name,                         &
            &                     lfeedback   = .TRUE.,                              &
            &                     ihadv_tracer=advconf%ihadv_tracer(iqi),            &
            &                     ivadv_tracer=advconf%ivadv_tracer(iqi)),           &
            &         vert_interp=create_vert_interp_metadata(                       &
            &                     vert_intp_type=vintp_types("P","Z","I"),           &
            &                     vert_intp_method=VINTP_METHOD_LIN,                 &
            &                     l_loglin=.FALSE.,                                  &
            &                     l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                     lower_limit=0._wp  ),                              &
            &         in_group=ingroup )
          __acc_attach(p_prog%tracer_ptr(iqi)%p_3d)
        END IF ! iqi

        !QR
        IF ( iqr /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqr))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqr)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(iqr)%p_3d,                      &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            'kg kg-1','specific rain content', datatype_flt),             &
            &           grib2_var(0, 1, 24, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=iqr,                                                   &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     &              ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,    &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(iqr),            &
            &                       ivadv_tracer=advconf%ivadv_tracer(iqr)),           &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       l_loglin=.FALSE.,                                  &
            &                       l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                       lower_limit=0._wp  ),                              &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",            &
            &                           "mode_iau_ana_in", "mode_iau_anaatm_in",       &
            &                           "mode_iau_fg_in","mode_iau_old_fg_in",         &
            &                           "LATBC_PREFETCH_VARS",                         &
            &                           "mode_iniana","icon_lbc_vars") )
          __acc_attach(p_prog%tracer_ptr(iqr)%p_3d)
        END IF

        !QS
        IF ( iqs /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqs))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqs)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(iqs)%p_3d,                      &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:tlen+vntl),                             &
            &            'kg kg-1','specific snow content', datatype_flt),             &
            &           grib2_var(0, 1, 25, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=iqs,                                                   &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,    &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(iqs),            &
            &                       ivadv_tracer=advconf%ivadv_tracer(iqs)),           &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       l_loglin=.FALSE.,                                  &
            &                       l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                       lower_limit=0._wp  ),                              &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",            &
            &                           "mode_iau_ana_in", "mode_iau_anaatm_in",       &
            &                           "mode_iau_fg_in","mode_iau_old_fg_in",         &
            &                           "LATBC_PREFETCH_VARS",                         &
            &                           "mode_iniana","icon_lbc_vars") )
          __acc_attach(p_prog%tracer_ptr(iqs)%p_3d)
        END IF

        !QG
        IF ( iqg /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqg))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqg)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(iqg)%p_3d,                      &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            'kg kg-1','specific graupel content', datatype_flt),          &
            &           grib2_var(0, 1, 32, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=iqg,                                                   &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,    &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(iqg),            &
            &                       ivadv_tracer=advconf%ivadv_tracer(iqg)),           &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       l_loglin=.FALSE.,                                  &
            &                       l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                       lower_limit=0._wp  ),                              &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",            &
            &                           "mode_iau_ana_in", "mode_iau_anaatm_in",       &
            &                           "mode_iau_fg_in","mode_iau_old_fg_in",         &
            &                           "LATBC_PREFETCH_VARS",                         &
            &                           "mode_iniana","icon_lbc_vars") )
          __acc_attach(p_prog%tracer_ptr(iqg)%p_3d)
        END IF


        !hail
        IF ( iqh /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqh))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqh)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                          &
            &           tracer_name, p_prog%tracer_ptr(iqh)%p_3d,                    &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
            &           t_cf_var(tracer_name(1:vntl+tlen),                           &
            &            'kg kg-1 ','specific hail content', datatype_flt),          &
            &           grib2_var(0, 1, 71, ibits, GRID_UNSTRUCTURED, GRID_CELL),    &
            &           ref_idx=iqh,                                                 &
            &           ldims=shape3d_c,                                             &
            &           tlev_source=TLEV_NNOW_RCF,                                   & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,  &
            &                     name        = tracer_name,                         &
            &                     ihadv_tracer=advconf%ihadv_tracer(iqh),            &
            &                     ivadv_tracer=advconf%ivadv_tracer(iqh)),           &
            &           vert_interp=create_vert_interp_metadata(                     &
            &                     vert_intp_type=vintp_types("P","Z","I"),           &
            &                     vert_intp_method=VINTP_METHOD_LIN,                 &
            &                     l_loglin=.FALSE.,                                  &
            &                     l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                     lower_limit=0._wp  ),                              &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",&
            &                         "dwd_fg_atm_vars","mode_dwd_fg_in",            &
            &                         "mode_iau_ana_in","mode_iau_anaatm_in",        &
            &                         "mode_iau_fg_in",                              &
            &                         "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqh)%p_3d)
        END IF


        ! liquid water (meltwater) on graupel (shortname "QG_LIQ")
        IF ( iqgl /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqgl))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqgl)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                           &
            &           tracer_name, p_prog%tracer_ptr(iqgl)%p_3d,                    &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                         &
            &           t_cf_var(tracer_name(1:vntl+tlen),                            &
            &            'kgkg-1 ','Specific mass of liquid water coating on graupel',&
            &            datatype_flt),                                               &
            &           grib2_var(0, 1, 113, ibits, GRID_UNSTRUCTURED, GRID_CELL),    & ! from shortName.def
            &           ref_idx=iqgl,                                                 &
            &           ldims=shape3d_c,                                              &
            &           tlev_source=TLEV_NNOW_RCF,                                    & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,   &
            &                     name        = tracer_name,                          &
            &                     ihadv_tracer=advconf%ihadv_tracer(iqgl),            &
            &                     ivadv_tracer=advconf%ivadv_tracer(iqgl)),           &
            &           vert_interp=create_vert_interp_metadata(                      &
            &                      vert_intp_type=vintp_types("P","Z","I"),           &
            &                      vert_intp_method=VINTP_METHOD_LIN,                 &
            &                      l_loglin=.FALSE.,                                  &
            &                      l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &                      lower_limit=0._wp  ),                              &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars", &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",           &
            &                           "mode_iau_ana_in","mode_iau_anaatm_in",       &
            &                           "mode_iau_fg_in",                             &
            &                           "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqgl)%p_3d)
        END IF

        ! liquid water (meltwater) on hail  (shortname "QH_LIQ")
        IF ( iqhl /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqhl))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqhl)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                          &
            &           tracer_name, p_prog%tracer_ptr(iqhl)%p_3d,                   &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
            &           t_cf_var(tracer_name(1:tlen+vntl),                           &
            &            'kgkg-1 ','Specific mass of liquid water coating on hail',  &
            &            datatype_flt),                                              &
            &           grib2_var(0, 1, 110, ibits, GRID_UNSTRUCTURED, GRID_CELL),   & ! from shortName.def
            &           ref_idx=iqhl,                                                &
            &           ldims=shape3d_c,                                             &
            &           tlev_source=TLEV_NNOW_RCF,                                   & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,  &
            &                     name        = tracer_name,                         &
            &                     ihadv_tracer=advconf%ihadv_tracer(iqhl),           &
            &                     ivadv_tracer=advconf%ivadv_tracer(iqhl)),          &
            &           vert_interp=create_vert_interp_metadata(                     &
            &                     vert_intp_type=vintp_types("P","Z","I"),           &
            &                     vert_intp_method=VINTP_METHOD_LIN,                 &
            &                     l_loglin=.FALSE.,                                  &
            &                     l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &           lower_limit=0._wp  ),                                        &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",&
            &                         "dwd_fg_atm_vars","mode_dwd_fg_in",            &
            &                         "mode_iau_ana_in","mode_iau_anaatm_in",        &
            &                         "mode_iau_fg_in",                              &
            &                         "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqhl)%p_3d)
        END IF

        !33 drop mass bins
        DO iqb = iqb_i, iqb_s
          IF ( iqbin(iqb) /= 0 ) THEN
              tracer_name = TRIM(vname_prefix)//TRIM(advconf%tracer_names(iqbin(iqb)))
              CALL add_ref( p_prog_list, tracer_container_name,                      &
                    & TRIM(tracer_name)//suffix, p_prog%tracer_ptr(iqbin(iqb))%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(TRIM(tracer_name),                                    &
                    &  'kgkg-1 ','drop mass bins', datatype_flt),                    &
                    & grib2_var(0, 1, 71, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                    & ref_idx=iqbin(iqb),                                            &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata_hydro(lis_tracer=.TRUE.,    &
                    &             name        = TRIM(tracer_name)//suffix,           &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqbin(iqb)),     &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqbin(iqb))),    &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          END IF
        END DO

        !33 ccn number bins
        DO iqb = iqb_s+1, iqb_e
          IF ( iqbin(iqb) /= 0 ) THEN
            tracer_name = TRIM(vname_prefix)//TRIM(advconf%tracer_names(iqbin(iqb)))
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & TRIM(tracer_name)//suffix, p_prog%tracer_ptr(iqbin(iqb))%p_3d, &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(TRIM(tracer_name),                                    &
                    &  ' kg-1 ','ccn number bins', datatype_flt),                    &
                    & grib2_var(0, 6, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                    & ref_idx=iqbin(iqb),                                            &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = TRIM(tracer_name)//suffix,           &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqbin(iqb)),     &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqbin(iqb))),    &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          END IF
        END DO

        !O3
        IF ( iqt <= io3 .AND. io3 <= ntracer) THEN
          tlen = LEN_TRIM(advconf%tracer_names(io3))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(io3)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(io3)%p_3d,                      &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            'kg kg-1','o3_mass_mixing_ratio', datatype_flt),              &
            &           grib2_var(0,14,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=io3,                                                   &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(io3),            &
            &                       ivadv_tracer=advconf%ivadv_tracer(io3)),           &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       lower_limit=0.0_wp )                               )
          __acc_attach(p_prog%tracer_ptr(io3)%p_3d)
        END IF

        !CO2
        IF ( iqt <= ico2 .AND. ico2 <= ntracer) THEN
          tlen = LEN_TRIM(advconf%tracer_names(ico2))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(ico2)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(ico2)%p_3d,                     &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            'kg kg-1','co2_mass_mixing_ratio', datatype_flt),             &
            &           grib2_var(0,14,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=ico2,                                                  &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(ico2),           &
            &                       ivadv_tracer=advconf%ivadv_tracer(ico2)),          &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       lower_limit=0.0_wp )                               )
          __acc_attach(p_prog%tracer_ptr(ico2)%p_3d)
        END IF

        !CH4
        IF ( iqt <= ich4 .AND. ich4 <= ntracer ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(ich4))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(ich4)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(ich4)%p_3d,                     &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            'kg kg-1','ch4_mass_mixing_ratio', datatype_flt),             &
            &           grib2_var(0,14,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=ich4,                                                  &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(ich4),           &
            &                       ivadv_tracer=advconf%ivadv_tracer(ich4)),          &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       lower_limit=0.0_wp )                               )
          __acc_attach(p_prog%tracer_ptr(ich4)%p_3d)
        END IF

        !N2O
        IF ( iqt <= in2o .AND. in2o <= ntracer) THEN
          tlen = LEN_TRIM(advconf%tracer_names(in2o))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(in2o)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(in2o)%p_3d,                     &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            'kg kg-1','n2o_mass_mixing_ratio', datatype_flt),             &
            &           grib2_var(0,14,255, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=in2o,                                                  &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(in2o),           &
            &                       ivadv_tracer=advconf%ivadv_tracer(in2o)),          &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       lower_limit=0.0_wp )                               )
          __acc_attach(p_prog%tracer_ptr(in2o)%p_3d)
        END IF


        !ice number concentration
        IF ( iqni /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqni))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqni)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                            &
            &           tracer_name, p_prog%tracer_ptr(iqni)%p_3d,                     &
            &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
            &           t_cf_var(tracer_name(1:vntl+tlen),                             &
            &            ' kg-1 ','number concentration cloud ice', datatype_flt),     &
            &           grib2_var(0, 6, 29, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            &           ref_idx=iqni,                                                  &
            &           ldims=shape3d_c,                                               &
            &           tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
            &           tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
            &                       name        = tracer_name,                         &
            &                       ihadv_tracer=advconf%ihadv_tracer(iqni),           &
            &                       ivadv_tracer=advconf%ivadv_tracer(iqni)),          &
            &           vert_interp=create_vert_interp_metadata(                       &
            &                       vert_intp_type=vintp_types("P","Z","I"),           &
            &                       vert_intp_method=VINTP_METHOD_LIN,                 &
            &                       l_loglin=.FALSE.,                                  &
            &                       l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
            &             lower_limit=0._wp  ),                                        &
            &           in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
            &                           "dwd_fg_atm_vars","mode_dwd_fg_in",            &
            &                           "mode_iau_ana_in", "mode_iau_anaatm_in",       &
            &                           "mode_iau_fg_in",                              &
            &                           "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqni)%p_3d)
        END IF

        !rain droplet concentration
        IF ( iqnr /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(iqnr))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqnr)(1:tlen)//suffix
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(iqnr)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),                             &
                    &  ' kg-1 ','number_concentration_rain_droplet', datatype_flt),  &
                    & grib2_var(0, 1, 100, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                    & ref_idx=iqnr,                                                  &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqnr),           &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqnr)),          &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqnr)%p_3d)
        END IF


        !snow concentration
        IF ( iqns /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(iqns))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqns)(1:tlen)//suffix
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(iqns)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),                             &
                    &  ' kg-1 ','number_concentration_snow', datatype_flt),          &
                    & grib2_var(0, 1, 101, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                    & ref_idx=iqns,                                                  &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqns),           &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqns)),          &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqns)%p_3d)
        END IF

        !graupel concentration
        IF ( iqng /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(iqng))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqng)(1:tlen)//suffix
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(iqng)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),                             &
                    &  ' kg-1 ','number_concentration_graupel', datatype_flt),       &
                    & grib2_var(0, 1, 102, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                    & ref_idx=iqng,                                                  &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqng),           &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqng)),          &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqng)%p_3d)
        END IF

        !hail concentration
        IF ( iqnh /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(iqnh))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqnh)(1:tlen)//suffix
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(iqnh)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),                             &
                    &  ' kg-1 ','number_concentration_hail', datatype_flt),          &
                    & grib2_var(0, 1, 103, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
                    & ref_idx=iqnh,                                                  &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqnh),           &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqnh)),          &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqnh)%p_3d)
        END IF

        ! cloud droplet concentration
        ! QNC  pdis=0 pcat=6 pnum=28 #DWD: Number of cloud droplets per unit mass of air. paramId=502315
        IF ( iqnc /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(iqnc))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqnc)(1:tlen)//suffix
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(iqnc)%p_3d,                     &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),                             &
                    &  ' kg-1 ','number_concentration_cloud_droplets', datatype_flt),&
                    & grib2_var(0, 6, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                    & ref_idx=iqnc,                                                  &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(iqnc),           &
                    &             ivadv_tracer=advconf%ivadv_tracer(iqnc)),          &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars",  &
                    &                 "dwd_fg_atm_vars","mode_dwd_fg_in",            &
                    &                 "mode_iau_ana_in","mode_iau_anaatm_in",        &
                    &                 "mode_iau_fg_in",                              &
                    &                 "LATBC_PREFETCH_VARS")  )
          __acc_attach(p_prog%tracer_ptr(iqnc)%p_3d)
        END IF


        IF ( ininact /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(ininact))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(ininact)(1:tlen)//suffix
            ! NO OFFICIAL GRIB CODINGS YET! THE "255, 255, 255" HAS TO BE ADAPTED WHEN THESE CODINGS BECOME AVAILABLE!
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(ininact)%p_3d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),' kg-1 ',                    &
                    & 'number_concentration_activated_ice_nuclei', datatype_flt),    &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & ref_idx=ininact,                                               &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(ininact),        &
                    &             ivadv_tracer=advconf%ivadv_tracer(ininact)),       &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars")  )
          __acc_attach(p_prog%tracer_ptr(ininact)%p_3d)
        END IF ! inwp_gscp==4,5,6,7,8

        ! concentration of cloud condensation nuclei
        IF ( inccn /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(inccn))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(inccn)(1:tlen)//suffix
            ! NO OFFICIAL GRIB CODINGS YET! THE "255, 255, 255" HAS TO BE ADAPTED WHEN THESE CODINGS BECOME AVAILABLE!
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(inccn)%p_3d,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),' kg-1 ',                    &
                    & 'number_concentration_cloud_condensation_nuclei',datatype_flt),&
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & ref_idx=inccn,                                                 &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(inccn),          &
                    &             ivadv_tracer=advconf%ivadv_tracer(inccn)),         &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars")  )
          __acc_attach(p_prog%tracer_ptr(inccn)%p_3d)
        END IF

        IF ( ininpot /= 0 ) THEN
            tlen = LEN_TRIM(advconf%tracer_names(ininpot))
            tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(ininpot)(1:tlen)//suffix
            ! NO OFFICIAL GRIB CODINGS YET! THE "255, 255, 255" HAS TO BE ADAPTED WHEN THESE CODINGS BECOME AVAILABLE!
            CALL add_ref( p_prog_list, tracer_container_name,                        &
                    & tracer_name, p_prog%tracer_ptr(ininpot)%p_3d,                  &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & t_cf_var(tracer_name(1:vntl+tlen),' kg-1 ',                    &
                    & 'number_concentration_potential_ice_nuclei', datatype_flt),    &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
                    & ref_idx=ininpot,                                               &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     & ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name        = tracer_name,                         &
                    &             ihadv_tracer=advconf%ihadv_tracer(ininpot),        &
                    &             ivadv_tracer=advconf%ivadv_tracer(ininpot)),       &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,            &
                    &             lower_limit=0._wp  ),                              &
                    & in_group=groups("atmo_ml_vars", "atmo_pl_vars", "atmo_zl_vars")  )
          __acc_attach(p_prog%tracer_ptr(ininpot)%p_3d)
        END IF


        IF ( iqtke /= 0 ) THEN
          tlen = LEN_TRIM(advconf%tracer_names(iqtke))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(iqtke)(1:tlen)//suffix
          cf_desc    = t_cf_var(tracer_name(1:vntl+tlen), 'm s-1',         &
            &          'turbulent velocity scale (at full levels)', datatype_flt)
          grib2_desc = grib2_var(0, 19, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_ref( p_prog_list, tracer_container_name,                          &
                    & tracer_name, p_prog%tracer_ptr(iqtke)%p_3d,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & cf_desc, grib2_desc,                                           &
                    & ref_idx=iqtke,                                                 &
                    & ldims=shape3d_c,                                               &
                    & tlev_source=TLEV_NNOW_RCF,                                     &              ! output from nnow_rcf slice
                    & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,          &
                    &             name = tracer_name),                               &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE.,                                  &
                    &             l_pd_limit=.FALSE.,                                &
                    &             lower_limit=0._wp  )  )
          __acc_attach(p_prog%tracer_ptr(iqtke)%p_3d)
        ENDIF

#ifdef __ICON_ART
        ! art
        IF (lart) THEN
          CALL art_tracer_interface('prog',p_patch%id,p_patch%nblks_c,p_prog_list,vname_prefix,&
            &                       ptr_arr=p_prog%tracer_ptr,advconf=advconf,p_prog=p_prog,   &
            &                       timelev=timelev,ldims=shape3d_c)
        ENDIF
#endif

        ! tke            p_prog%tke(nproma,nlevp1,nblks_c)
        ! for output take field from nnow_rcf slice
        cf_desc    = t_cf_var('specific_kinetic_energy_of_air', 'm2 s-2',         &
          &          'turbulent kinetic energy', datatype_flt)
        grib2_desc = grib2_var(0, 19, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_prog_list, vname_prefix(1:vntl)//'tke'//suffix, p_prog%tke, &
          &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                  &
          &           cf_desc, grib2_desc, ldims=shape3d_chalf,                   &
          &           tlev_source=TLEV_NNOW_RCF,                                  &
          &           vert_interp=create_vert_interp_metadata(                    &
          &                       vert_intp_type=vintp_types("P","Z","I"),        &
          &                       vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),     &
          &           in_group=groups("atmo_ml_vars", "atmo_pl_vars",             &
          &           "atmo_zl_vars", "dwd_fg_atm_vars", "mode_dwd_fg_in",        &
          &           "mode_iau_fg_in","mode_iau_old_fg_in",                      &
          &           "mode_iniana","icon_lbc_vars","iau_restore_vars"),          &
          &           lopenacc = .TRUE.  )
        __acc_attach(p_prog%tke)
        
      ELSE

        ! add refs to tracers with generic name q1 ... qn
        ! (used for example with test cases)

        DO jt = 1, ntracer - advection_config(p_patch%id)%npassive_tracer
          tlen = LEN_TRIM(advconf%tracer_names(jt))
          tracer_name = vname_prefix(1:vntl)//advconf%tracer_names(jt)(1:tlen)//suffix
          CALL add_ref( p_prog_list, tracer_container_name,                     &
            & tracer_name, p_prog%tracer_ptr(jt)%p_3d,                          &
            & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                             &
            & t_cf_var(tracer_name(1:vntl+tlen),                                &
            &          'kg/kg','Tracer mixing ratio '//tracer_name(1:vntl+tlen),&
            & datatype_flt),                                                    &
            & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),          &
            & ref_idx=jt,                                                       &
            & ldims=shape3d_c,                                                  &
            & tlev_source=TLEV_NNOW_RCF,                                        &              ! output from nnow_rcf slice
            & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,             &
            &                       name        = tracer_name,                  &
            &                       lfeedback   = .TRUE.,                       &
            &                       ihadv_tracer=advconf%ihadv_tracer(jt),      &
            &                       ivadv_tracer=advconf%ivadv_tracer(jt)),     &
            & vert_interp=create_vert_interp_metadata(                          &
            &             vert_intp_type=vintp_types("P","Z","I"),              &
            &             vert_intp_method=VINTP_METHOD_LIN,                    &
            &             l_loglin=.FALSE.,                                     &
            &             l_extrapol=.FALSE., l_pd_limit=.FALSE.,               &
            &             lower_limit=0._wp  )  )
          __acc_attach(p_prog%tracer_ptr(jt)%p_3d)
        END DO
      ENDIF


      ! add references to additional passive tracers, if existing
      DO ipassive=1,advection_config(p_patch%id)%npassive_tracer

        ! Determine index of the following tracer within the 4D tracer container.
        ! We need this information, in order to pass some metadata from the configure 
        ! state into the tracer_info metadata storage.
        !
        ! get pointer to target element (in this case 4D tracer container)
        target_element => find_list_element (p_prog_list, 'tracer'//suffix)
        tracer_idx = target_element%info%ncontained+1

        WRITE(passive_tracer_suffix,'(I2)') ipassive
        tracer_name = 'Qpassive_'//passive_tracer_suffix(1+MERGE(1,0,ipassive<=9):)
        tlen = LEN_TRIM(tracer_name)
        tracer_name(tlen+1:) = suffix
        cf_desc    = t_cf_var(tracer_name(1:tlen), 'kg kg-1', 'passive tracer', datatype_flt)
        grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_tracer_ref( p_prog_list, tracer_container_name,                     &
          &                  tracer_name,                                            &
          &                  dummy_idx,                                              &
          &                  p_prog%tracer_ptr(:),                                   &
          &                  cf_desc, grib2_desc,                                    &
          &                  advection_config(p_patch%id),                           &
          &                  ldims=shape3d_c,                                        &
          &                  loutput=.TRUE.,                                         &
          &                  lrestart=.FALSE.,                                       &
          &                  tlev_source=TLEV_NNOW_RCF,                              &  ! output from nnow_rcf slice
          &                  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,   &
          &                              name = tracer_name,                         &
          &                              lfeedback   = .TRUE.,                       &
          &                              ihadv_tracer=advconf%ihadv_tracer(tracer_idx), &
          &                              ivadv_tracer=advconf%ivadv_tracer(tracer_idx)) )
      ENDDO

    ENDIF ! allocation only if not extra_timelev


  END SUBROUTINE new_nh_state_prog_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic state.
  !!
  !! Initialization of components with zero.
  !!
  SUBROUTINE new_nh_state_diag_list ( p_patch, p_diag, p_diag_list,  &
    &                                 listname, var_in_output )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_diag),  INTENT(INOUT)   :: &  !< diagnostic state
      &  p_diag 

    TYPE(t_var_list_ptr), INTENT(INOUT)   :: &  !< diagnostic state list
      &  p_diag_list
    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname
    TYPE(t_var_in_output), INTENT(IN) :: &  !< optional diagnostic switches
      &  var_in_output

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e, &    !< number of edge blocks to allocate
               nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev, nlevp1

    INTEGER :: n_timlevs     !< number of time levels for advection 
                             !< tendency fields

    INTEGER :: shape2d_c(2),  shape3d_c(3), shape3d_e(3),          &
      &        shape3d_e2(3), shape3d_v(3), shape3d_chalf(3),      &
      &        shape3d_ehalf(3), shape4d_chalf(4), shape4d_e(4),   &
      &        shape4d_entl(4), shape4d_chalfntl(4), shape4d_c(4), &
      &        shape2d_extra(3), shape3d_extra(4), shape3d_ubcc(3),&
      &        shape3d_ubcp2(3), shape4d_cams(4)
 
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for some thermodynamic fields

    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    INTEGER :: jt

    LOGICAL :: lrestart

    CHARACTER(LEN=3) :: ctrc
    CHARACTER(len=4) :: suffix
    !--------------------------------------------------------------

    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    IF (itime_scheme >= 2) THEN
     n_timlevs = 2
    ELSE
     n_timlevs = 1
    ENDIF

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! predefined array shapes
    shape2d_c     = (/nproma,          nblks_c    /)
    shape2d_extra = (/nproma, nblks_c, inextra_2d /)
    shape3d_c     = (/nproma, nlev   , nblks_c    /)
    shape3d_e     = (/nproma, nlev   , nblks_e    /)
    shape3d_e2    = (/nproma, 2      , nblks_e    /)
    shape3d_v     = (/nproma, nlev   , nblks_v    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)
    shape3d_ehalf = (/nproma, nlevp1 , nblks_e    /)
    shape3d_ubcp2 = (/nproma, nblks_c, ndyn_substeps_max+2 /)
    shape3d_ubcc  = (/nproma, nblks_c, 2  /)
    shape3d_extra = (/nproma, nlev   , nblks_c, inextra_3d  /)
    shape4d_cams  = (/nproma, nlev   , nblks_c, n_camsaermr /)
    shape4d_c     = (/nproma, nlev   , nblks_c, ntracer     /)
    shape4d_chalf = (/nproma, nlevp1 , nblks_c, ntracer     /)
    shape4d_e     = (/nproma, nlev   , nblks_e, ntracer     /)
    shape4d_entl  = (/nproma, nlev   , nblks_e, n_timlevs   /)
    shape4d_chalfntl = (/nproma, nlevp1, nblks_c, n_timlevs /)

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(p_diag%u, &
    &       p_diag%v, &
    &       p_diag%vt, &
    &       p_diag%omega_z, &
    &       p_diag%vor, &
    &       p_diag%ddt_vn_phy, &
    &       p_diag%ddt_vn_dyn, &
    &       p_diag%ddt_ua_dyn, &
    &       p_diag%ddt_va_dyn, &
    &       p_diag%ddt_vn_dmp, &
    &       p_diag%ddt_ua_dmp, &
    &       p_diag%ddt_va_dmp, &
    &       p_diag%ddt_vn_hdf, &
    &       p_diag%ddt_ua_hdf, &
    &       p_diag%ddt_va_hdf, &
    &       p_diag%ddt_vn_adv, &
    &       p_diag%ddt_ua_adv, &
    &       p_diag%ddt_va_adv, &
    &       p_diag%ddt_vn_cor, &
    &       p_diag%ddt_ua_cor, &
    &       p_diag%ddt_va_cor, &
    &       p_diag%ddt_vn_pgr, &
    &       p_diag%ddt_ua_pgr, &
    &       p_diag%ddt_va_pgr, &
    &       p_diag%ddt_vn_phd, &
    &       p_diag%ddt_ua_phd, &
    &       p_diag%ddt_va_phd, &
    &       p_diag%ddt_vn_iau, &
    &       p_diag%ddt_ua_iau, &
    &       p_diag%ddt_va_iau, &
    &       p_diag%ddt_vn_ray, &
    &       p_diag%ddt_ua_ray, &
    &       p_diag%ddt_va_ray, &
    &       p_diag%ddt_vn_grf, &
    &       p_diag%ddt_ua_grf, &
    &       p_diag%ddt_va_grf, &
    &       p_diag%ddt_exner_phy, &
    &       p_diag%ddt_temp_dyn, &
    &       p_diag%ddt_tracer_adv, &
    &       p_diag%tracer_vi, &
    &       p_diag%exner_pr, &
    &       p_diag%exner_dyn_incr, &
    &       p_diag%temp, &
    &       p_diag%tempv, &
    &       p_diag%temp_ifc, &
    &       p_diag%pres, &
    &       p_diag%pres_ifc, &
    &       p_diag%pres_sfc, &
    &       p_diag%pres_sfc_old, &
    &       p_diag%ddt_pres_sfc, &
    &       p_diag%pres_msl, &
    &       p_diag%dpres_mc, &
    &       p_diag%omega, &
    &       p_diag%hfl_tracer, &
    &       p_diag%vfl_tracer, &
    &       p_diag%div, &
    &       p_diag%div_ic, &
    &       p_diag%hdef_ic, &
    &       p_diag%kh_smag_e, &
    &       p_diag%dwdx, &
    &       p_diag%dwdy, &
    &       p_diag%mass_fl_e, &
    &       p_diag%mass_fl_e_sv, &
    &       p_diag%rho_ic, &
    &       p_diag%theta_v_ic, &
    &       p_diag%w_concorr_c, &
    &       p_diag%vn_ie, &
    &       p_diag%ddt_vn_apc_pc, &
    &       p_diag%ddt_vn_cor_pc, &
    &       p_diag%ddt_w_adv_pc, &
    &       p_diag%airmass_now, &
    &       p_diag%airmass_new, &
    &       p_diag%grf_tend_vn, &
    &       p_diag%grf_tend_w, &
    &       p_diag%grf_tend_rho, &
    &       p_diag%grf_tend_mflx, &
    &       p_diag%grf_bdy_mflx, &
    &       p_diag%grf_tend_thv, &
    &       p_diag%grf_tend_tracer, &
    &       p_diag%vn_ie_int, &
    &       p_diag%vn_ie_ubc, &
    &       p_diag%w_int, &
    &       p_diag%w_ubc, &
    &       p_diag%theta_v_ic_int, &
    &       p_diag%theta_v_ic_ubc, &
    &       p_diag%rho_ic_int, &
    &       p_diag%rho_ic_ubc, &
    &       p_diag%mflx_ic_int, &
    &       p_diag%mflx_ic_ubc, &
    &       p_diag%vn_incr, &
    &       p_diag%exner_incr, &
    &       p_diag%rho_incr, &
    &       p_diag%rhov_incr, &
    &       p_diag%rhoc_incr, &
    &       p_diag%rhoi_incr, &
    &       p_diag%rhor_incr, &
    &       p_diag%rhos_incr, &
    &       p_diag%rhog_incr, &
    &       p_diag%rhoh_incr, &
    &       p_diag%rhonc_incr, &
    &       p_diag%rhoni_incr, &
    &       p_diag%rhonr_incr, &
    &       p_diag%rhons_incr, &
    &       p_diag%rhong_incr, &
    &       p_diag%rhonh_incr, &
    &       p_diag%camsaermr,  &
    &       p_diag%extra_2d,   &
    &       p_diag%extra_3d)



    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_diag_list, TRIM(listname), patch_id=p_patch%id, &
      &          lrestart=.TRUE., model_type=get_my_process_name())

    ! u           p_diag%u(nproma,nlev,nblks_c)
    !
    ! ATTENTION: the vertical interpolation of "u" and "v" results from a vertical
    !            interpolation of the normal velocity "vn" on the edge points, followed
    !            by the edge-to-center opertor that computes (u,v) from vn. Therefore the
    !            the "vert_intp_method" specified for "vn" is used and and definition
    !            of "vert_intp_method" for "u" and "v" is ignored.
    !
    cf_desc    = t_cf_var('eastward_wind', 'm s-1', 'Zonal wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'u', p_diag%u,                                   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I") ),        &
                & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars", &
                &                 "dwd_fg_atm_vars","mode_dwd_ana_in",          &
                &                 "mode_iau_ana_in","mode_iau_old_ana_in",      &
                &                 "mode_iau_anaatm_in","LATBC_PREFETCH_VARS",   &
                &                 "mode_iniana","icon_lbc_vars"),               & 
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%u)

    ! v           p_diag%v(nproma,nlev,nblks_c)
    !
    ! ATTENTION: see comment for "u".
    !
    cf_desc    = t_cf_var('northward_wind', 'm s-1', 'Meridional wind', datatype_flt)
    grib2_desc = grib2_var(0, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'v', p_diag%v,                                   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I") ),        &
                & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars", &
                &                 "dwd_fg_atm_vars","mode_dwd_ana_in",          &
                &                 "mode_iau_ana_in","mode_iau_old_ana_in",      &
                &                 "mode_iau_anaatm_in","LATBC_PREFETCH_VARS",   &
                &                 "mode_iniana","icon_lbc_vars"),               &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%v)

    IF (irad_aero == iRadAeroCAMSclim .OR. irad_aero == iRadAeroCAMStd) THEN
      cf_desc    = t_cf_var('CAMS_aerosols', 'kg kg-1', 'CAMS aerosols mixing ratios', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'camsaermr', p_diag%camsaermr,                   &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & initval=0._wp,                                                &
                  & lcontainer=.TRUE., ldims=shape4d_cams, lrestart=.FALSE.)
    END IF

    ! vt           p_diag%vt(nproma,nlev,nblks_e)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('tangential_wind', 'm s-1', 'tangential-component of wind', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 2, 35, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'vt', p_diag%vt,                                          &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,             &
                & ldims=shape3d_e,                                                       &
                & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_NONE ), &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%vt)

    ! omega_z      p_diag%omega_z(nproma,nlev,nblks_v)
    !
    cf_desc    = t_cf_var('atmospheric_relative_vorticity', 'm s-1', 'vertical vorticity', &
      &          datatype_flt)
    grib2_desc = grib2_var(0, 2, 12, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
    CALL add_var( p_diag_list, 'omega_z', p_diag%omega_z,                       &
                & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_v, lrestart=.FALSE., in_group=groups("atmo_derived_vars"), &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%omega_z)

    ! ddt_vn_phy   p_diag%ddt_vn_phy(nproma,nlev,nblks_e)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('normal_wind_physical_tendency', 'm s-2',             &
      &                   'normal wind physical tendency', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'ddt_vn_phy', p_diag%ddt_vn_phy,                 &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e,                                              &
                & in_group=groups("iau_init_vars"),                             &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%ddt_vn_phy)


    IF (var_in_output%ddt_vn_dyn .OR. var_in_output%ddt_ua_dyn .OR. var_in_output%ddt_va_dyn) THEN
      ! ddt_vn_dyn   p_diag%ddt_vn_dyn(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_dynamics', 'm s-2', &
        &                   'normal wind tendency by dynamics', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 220, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_dyn', p_diag%ddt_vn_dyn,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_dyn)
      p_diag%ddt_vn_dyn_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_dyn_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_dyn .OR. var_in_output%ddt_va_dyn) THEN
      ! ddt_ua_dyn   p_diag%ddt_ua_dyn(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_dynamics', 'm s-2', &
        &                   'zonal wind tendency by dynamics', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 230, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_dyn', p_diag%ddt_ua_dyn,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_dyn)
      p_diag%ddt_ua_dyn_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_dyn_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_dyn .OR. var_in_output%ddt_va_dyn) THEN
      ! ddt_va_dyn   p_diag%ddt_va_dyn(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_dynamics', 'm s-2', &
        &                   'meridional wind tendency by dynamics', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 240, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_dyn', p_diag%ddt_va_dyn,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_dyn)
      p_diag%ddt_va_dyn_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_dyn_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_dmp .OR. var_in_output%ddt_ua_dmp .OR. var_in_output%ddt_va_dmp) THEN
     ! ddt_vn_dmp   p_diag%ddt_vn_dmp(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_divergence_damping', 'm s-2', &
        &                   'normal wind tendency by divergence damping', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 221, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_dmp', p_diag%ddt_vn_dmp,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_dmp)
      p_diag%ddt_vn_dmp_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_dmp_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_dmp .OR. var_in_output%ddt_va_dmp) THEN
      ! ddt_ua_dmp   p_diag%ddt_ua_dmp(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_divergence_damping', 'm s-2', &
        &                   'zonal wind tendency by divergence damping', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 231, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_dmp', p_diag%ddt_ua_dmp,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_dmp)
      p_diag%ddt_ua_dmp_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_dmp_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_dmp .OR. var_in_output%ddt_va_dmp) THEN
      ! ddt_va_dmp   p_diag%ddt_va_dmp(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_divergence_damping', 'm s-2', &
        &                   'meridional wind tendency by divergence damping', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 241, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_dmp', p_diag%ddt_va_dmp,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_dmp)
      p_diag%ddt_va_dmp_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_dmp_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_hdf .OR. var_in_output%ddt_ua_hdf .OR. var_in_output%ddt_va_hdf) THEN
      ! ddt_vn_hdf   p_diag%ddt_vn_hdf(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_horizontal_diffusion', 'm s-2', &
        &                   'normal wind tendency by horizontal diffusion', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 222, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_hdf', p_diag%ddt_vn_hdf,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_hdf)
      p_diag%ddt_vn_hdf_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_hdf_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_hdf .OR. var_in_output%ddt_va_hdf) THEN
      ! ddt_ua_hdf   p_diag%ddt_ua_hdf(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_horizontal_diffusion', 'm s-2', &
        &                   'zonal wind tendency by horizontal diffusion', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 232, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_hdf', p_diag%ddt_ua_hdf,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_hdf)
      p_diag%ddt_ua_hdf_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_hdf_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_hdf .OR. var_in_output%ddt_va_hdf) THEN
      ! ddt_va_hdf   p_diag%ddt_va_hdf(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_horizontal_diffusion', 'm s-2', &
        &                   'meridional wind tendency by horizontal diffusion', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 242, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_hdf', p_diag%ddt_va_hdf,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_hdf)
      p_diag%ddt_va_hdf_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_hdf_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_adv .OR. var_in_output%ddt_ua_adv .OR. var_in_output%ddt_va_adv) THEN
      ! ddt_vn_adv   p_diag%ddt_vn_adv(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_advection', 'm s-2', &
        &                   'normal wind tendency by advection', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 223, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_adv', p_diag%ddt_vn_adv,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_adv)
      p_diag%ddt_vn_adv_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_adv_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_adv .OR. var_in_output%ddt_va_adv) THEN
      ! ddt_ua_adv   p_diag%ddt_ua_adv(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_advection', 'm s-2', &
        &                   'zonal wind tendency by advection', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 233, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_adv', p_diag%ddt_ua_adv,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_adv)
      p_diag%ddt_ua_adv_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_adv_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_adv .OR. var_in_output%ddt_va_adv) THEN
      ! ddt_va_adv   p_diag%ddt_va_adv(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_advection', 'm s-2', &
        &                   'meridional wind tendency by advection', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 243, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_adv', p_diag%ddt_va_adv,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_adv)
      p_diag%ddt_va_adv_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_adv_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_cor .OR. var_in_output%ddt_ua_cor .OR. var_in_output%ddt_va_cor) THEN
      ! ddt_vn_cor   p_diag%ddt_vn_cor(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_coriolis_effect', 'm s-2', &
        &                   'normal wind tendency by coriolis effect', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 224, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_cor', p_diag%ddt_vn_cor,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_cor)
      p_diag%ddt_vn_cor_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_cor_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_cor .OR. var_in_output%ddt_va_cor) THEN
      ! ddt_ua_cor   p_diag%ddt_ua_cor(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_coriolis_effect', 'm s-2', &
        &                   'zonal wind tendency by coriolis effect', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 234, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_cor', p_diag%ddt_ua_cor,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_cor)
      p_diag%ddt_ua_cor_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_cor_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_cor .OR. var_in_output%ddt_va_cor) THEN
      ! ddt_va_cor   p_diag%ddt_va_cor(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_coriolis_effect', 'm s-2', &
        &                   'meridional wind tendency by coriolis effect', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 244, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_cor', p_diag%ddt_va_cor,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_cor)
      p_diag%ddt_va_cor_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_cor_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_pgr .OR. var_in_output%ddt_ua_pgr .OR. var_in_output%ddt_va_pgr) THEN
      ! ddt_vn_pgr   p_diag%ddt_vn_pgr(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_pressure_gradient', 'm s-2', &
        &                   'normal wind tendency by pressure gradient', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 225, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_pgr', p_diag%ddt_vn_pgr,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_pgr)
      p_diag%ddt_vn_pgr_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_pgr_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_pgr .OR. var_in_output%ddt_va_pgr) THEN
      ! ddt_ua_pgr   p_diag%ddt_ua_pgr(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_pressure_gradient', 'm s-2', &
        &                   'zonal wind tendency by pressure gradient', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 235, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_pgr', p_diag%ddt_ua_pgr,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_pgr)
      p_diag%ddt_ua_pgr_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_pgr_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_pgr .OR. var_in_output%ddt_va_pgr) THEN
      ! ddt_va_pgr   p_diag%ddt_va_pgr(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_pressure_gradient', 'm s-2', &
        &                   'meridional wind tendency by pressure gradient', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 245, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_pgr', p_diag%ddt_va_pgr,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_pgr)
      p_diag%ddt_va_pgr_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_pgr_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_phd .OR. var_in_output%ddt_ua_phd .OR. var_in_output%ddt_va_phd) THEN
      ! ddt_vn_phd   p_diag%ddt_vn_phd(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_physics_in_dynamics', 'm s-2', &
        &                   'normal wind tendency by physics in dynamics', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 226, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_phd', p_diag%ddt_vn_phd,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_phd)
      p_diag%ddt_vn_phd_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_phd_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_phd .OR. var_in_output%ddt_va_phd) THEN
      ! ddt_ua_phd   p_diag%ddt_ua_phd(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_physics_in_dynamics', 'm s-2', &
        &                   'zonal wind tendency by physics in dynamics', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 236, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_phd', p_diag%ddt_ua_phd,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_phd)
      p_diag%ddt_ua_phd_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_phd_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_phd .OR. var_in_output%ddt_va_phd) THEN
      ! ddt_va_phd   p_diag%ddt_va_phd(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_physics_in_dynamics', 'm s-2', &
        &                   'meridional wind tendency by physics in dynamics', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 246, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_phd', p_diag%ddt_va_phd,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_phd)
      p_diag%ddt_va_phd_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_phd_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_iau .OR. var_in_output%ddt_ua_iau .OR. var_in_output%ddt_va_iau) THEN
      ! ddt_vn_iau   p_diag%ddt_vn_iau(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_incremental_analysis_update', 'm s-2', &
        &                   'normal wind tendency by incremental analysis update', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 228, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_iau', p_diag%ddt_vn_iau,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_iau)
      p_diag%ddt_vn_iau_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_iau_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_iau .OR. var_in_output%ddt_va_iau) THEN
      ! ddt_ua_iau   p_diag%ddt_ua_iau(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_incremental_analysis_update', 'm s-2', &
        &                   'zonal wind tendency by incremental analysis update', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 238, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_iau', p_diag%ddt_ua_iau,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_iau)
      p_diag%ddt_ua_iau_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_iau_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_iau .OR. var_in_output%ddt_va_iau) THEN
      ! ddt_va_iau   p_diag%ddt_va_iau(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_incremental_analysis_update', 'm s-2', &
        &                   'meridional wind tendency by incremental analysis update', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 248, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_iau', p_diag%ddt_va_iau,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_iau)
      p_diag%ddt_va_iau_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_iau_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_ray .OR. var_in_output%ddt_ua_ray .OR. var_in_output%ddt_va_ray) THEN
      ! ddt_vn_ray   p_diag%ddt_vn_ray(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_rayleigh_damping', 'm s-2', &
         &                  'normal wind tendency by rayleigh damping', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 229, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_ray', p_diag%ddt_vn_ray,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_ray)
      p_diag%ddt_vn_ray_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_ray_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_ray .OR. var_in_output%ddt_va_ray) THEN
      ! ddt_ua_ray   p_diag%ddt_ua_ray(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_rayleigh_damping', 'm s-2', &
        &                   'zonal wind tendency by rayleigh damping', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 239, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_ray', p_diag%ddt_ua_ray,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_ray)
      p_diag%ddt_ua_ray_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_ray_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_ray .OR. var_in_output%ddt_va_ray) THEN
      ! ddt_va_ray   p_diag%ddt_va_ray(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_rayleigh_damping', 'm s-2', &
        &                   'meridional wind tendency by rayleigh damping', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 249, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_ray', p_diag%ddt_va_ray,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_ray)
      p_diag%ddt_va_ray_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_ray_is_associated) ASYNC(1)
    END IF


    IF (var_in_output%ddt_vn_grf .OR. var_in_output%ddt_ua_grf .OR. var_in_output%ddt_va_grf) THEN
      ! ddt_vn_grf   p_diag%ddt_vn_grf(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('normal_wind_tendency_by_grid_refinement', 'm s-2',  &
         &                  'normal wind tendency by grid refinement', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 250, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_grf', p_diag%ddt_vn_grf,                &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_vn_grf)
      p_diag%ddt_vn_grf_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_vn_grf_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_grf .OR. var_in_output%ddt_va_grf) THEN
      ! ddt_ua_grf   p_diag%ddt_ua_grf(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('zonal_wind_tendency_by_grid_refinement', 'm s-2',   &
        &                   'zonal wind tendency by grid refinement', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 251, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_ua_grf', p_diag%ddt_ua_grf,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_ua_grf)
      p_diag%ddt_ua_grf_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_ua_grf_is_associated) ASYNC(1)
    END IF

    IF (var_in_output%ddt_ua_grf .OR. var_in_output%ddt_va_grf) THEN
      ! ddt_va_grf   p_diag%ddt_va_grf(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('meridional_wind_tendency_by_grid_refinement', 'm s-2', &
        &                   'meridional wind tendency by grid refinement', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 252, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_va_grf', p_diag%ddt_va_grf,                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_c ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%ddt_va_grf)
      p_diag%ddt_va_grf_is_associated=.TRUE.
      !$ACC UPDATE DEVICE(p_diag%ddt_va_grf_is_associated) ASYNC(1)
    END IF


    ! ddt_exner_phy  p_diag%ddt_exner_phy(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('exner_pressure_physical_tendency', 's-1',            &
      &                   'exner pressure physical tendency', datatype_flt)
    grib2_desc = grib2_var(0, 3, 197, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_exner_phy', p_diag%ddt_exner_phy,           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%ddt_exner_phy)


    ! ddt_temp_dyn  p_diag%ddt_temp_dyn(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('dynamical_temperature_tendency', 'K s-1',            &
      &                   'dynamical temperature tendency', datatype_flt)
    grib2_desc = grib2_var(0, 3, 197, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_temp_dyn', p_diag%ddt_temp_dyn,             &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & in_group=groups("iau_init_vars"),                             &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%ddt_temp_dyn)

    ! exner_pr    p_diag%exner_pr(nproma,nlev,nblks_c)
    ! *** needs to be saved for restart ***
    cf_desc    = t_cf_var('exner_perturbation_pressure', '-', 'exner perturbation pressure', datatype_flt)
    grib2_desc = grib2_var(0, 3, 26, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_pr', p_diag%exner_pr,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%exner_pr)

    ! exner_dyn_incr    p_diag%exner_dyn_incr(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('exner_dynamics_increment', '-', 'exner dynamics increment', &
         &                datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'exner_dyn_incr', p_diag%exner_dyn_incr,         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & in_group=groups("iau_init_vars"),                             &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%exner_dyn_incr)

    ! pres_sfc     p_diag%pres_sfc(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('surface_air_pressure', 'Pa', 'surface pressure', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_sfc', p_diag%pres_sfc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d_c, lrestart=.FALSE.,                            &
                & in_group=groups("dwd_fg_atm_vars", "LATBC_PREFETCH_VARS" ),   &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%pres_sfc)


    IF (lcalc_dpsdt) THEN
     ! pres_sfc_old     p_diag%pres_sfc_old(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('surface_pressure', 'Pa', 'surface pressure', datatype_flt)
      grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'pres_sfc_old', p_diag%pres_sfc_old,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE., loutput=.FALSE.,           &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%pres_sfc_old)

      ! ddt_pres_sfc     p_diag%ddt_pres_sfc(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('surface_pressure tendency', 'Pa s-1', 'surface pressure tendency', datatype_flt)
      grib2_desc = grib2_var(0, 3, 2, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_pres_sfc', p_diag%ddt_pres_sfc,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                  & ldims=shape2d_c, lrestart=.FALSE.,                            &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%ddt_pres_sfc)
    ENDIF

    ! temp         p_diag%temp(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('air_temperature', 'K', 'Temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'temp', p_diag%temp,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & hor_interp=create_hor_interp_metadata(                        &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,              &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),              &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & in_group=groups("atmo_ml_vars","atmo_pl_vars","atmo_zl_vars", &
                &                 "dwd_fg_atm_vars","mode_dwd_ana_in",          &
                &                 "mode_iau_ana_in","mode_iau_old_ana_in",      &
                &                 "mode_iau_anaatm_in","LATBC_PREFETCH_VARS",   &
                &                 "mode_iniana","icon_lbc_vars"),               &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%temp)

    ! tempv        p_diag%tempv(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_temperature', 'K', 'Virtual temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'tempv', p_diag%tempv,                           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%tempv)

    ! temp_ifc     p_diag%temp_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('air_temperature', 'K', 'temperature at half level', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'temp_ifc', p_diag%temp_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=lart,                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),       &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%temp_ifc)


    ! pres         p_diag%pres(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('air_pressure', 'Pa', 'Pressure', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'pres', p_diag%pres,                             &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE. ,                           &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("Z","I"),              &
                &             vert_intp_method=VINTP_METHOD_PRES ),             &
                & in_group=groups("atmo_ml_vars","atmo_zl_vars",                &
                & "dwd_fg_atm_vars","mode_dwd_ana_in",                          &
                & "mode_iau_ana_in","mode_iau_old_ana_in",                      &
                & "mode_iau_anaatm_in","LATBC_PREFETCH_VARS",                   &
                & "mode_iniana","icon_lbc_vars"),                               &
                &  lopenacc = .TRUE. )
    __acc_attach(p_diag%pres)

    ! pres_ifc     p_diag%pres_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('air_pressure', 'Pa', 'pressure at half level', datatype_flt)
    grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'pres_ifc', p_diag%pres_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE.,                        &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("Z","I"),              &
                &             vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),       &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%pres_ifc)


    ! dpres_mc     p_diag%dpres_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('pressure_thickness', 'Pa', 'pressure thickness', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'dpres_mc', p_diag%dpres_mc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%dpres_mc)


    ! div          p_diag%div(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('divergence_of_wind', 's-1', 'Divergence', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'div', p_diag%div,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & in_group=groups("atmo_derived_vars"),                         &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%div)

    IF (turbdiff_config(p_patch%id)%itype_sher >= 1 .OR. turbdiff_config(p_patch%id)%ltkeshs) THEN
      ! div_ic          p_diag%div_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('divergence at half levels', 's-1', 'divergence at half levels', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'div_ic', p_diag%div_ic,                          &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & ldims=shape3d_chalf, lrestart=.FALSE.,                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%div_ic)


      ! hdef_ic          p_diag%hdef_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('horizontal wind field deformation', 's-2', 'Deformation', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'hdef_ic', p_diag%hdef_ic,                            &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf, lrestart=.FALSE.,                             &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%hdef_ic)

    ELSE ! dummy allocation to satisfy the strict NAG compiler
      ALLOCATE(p_diag%div_ic(1,1,nblks_c), p_diag%hdef_ic(1,1,nblks_c))
    ENDIF

    IF (turbdiff_config(p_patch%id)%itype_sher >= 2) THEN
      ! dwdx          p_diag%dwdx(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Zonal gradient of vertical wind', 's-1', 'Zonal gradient of vertical wind', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'dwdx', p_diag%dwdx,                                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf, lrestart=.FALSE.,                             &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%dwdx)

      ! dwdy          p_diag%dwdy(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Meridional gradient of vertical wind', 's-1', 'Meridional gradient of vertical wind', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'dwdy', p_diag%dwdy,                                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
                  & ldims=shape3d_chalf, lrestart=.FALSE.,                             &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%dwdy)

    ELSE ! dummy allocation to satisfy the strict NAG compiler
      ALLOCATE(p_diag%dwdx(1,1,nblks_c), p_diag%dwdy(1,1,nblks_c))
    ENDIF

    IF (diffusion_config(p_patch%id)%lhdiff_q) THEN
      ! kh_smag_e   p_diag%kh_smag_e(nproma,nlev,nblks_e)
      cf_desc    = t_cf_var('horizontal Smagorinsky diffusion coefficient', 'm^2 s-1', &
        &                   'horizontal Smagorinsky diffusion coefficient', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'kh_smag_e', p_diag%kh_smag_e,                  &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & cf_desc, grib2_desc,                                         &
                  & ldims=shape3d_e ,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                            &
                  & lopenacc = .TRUE.                                            )
      __acc_attach(p_diag%kh_smag_e)
    END IF

    ! vor          p_diag%vor(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('relative_vorticity_on_cells', 's-1', 'Vorticity', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 12, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'vor', p_diag%vor,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, lrestart=.FALSE.,                            &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN ),              &
                & in_group=groups("atmo_derived_vars"),                         &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%vor)


    ! mass_fl_e    p_diag%mass_fl_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('horizontal_mass_flux_at_edges', 'kg m-1 s-1',        &
       &         'horizontal mass flux at edges', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'mass_fl_e', p_diag%mass_fl_e,                   &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, lrestart=.FALSE.,                            &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%mass_fl_e)


    ! mass_fl_e_sv    p_diag%mass_fl_e_sv(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('storage_field_for_horizontal_mass_flux_at_edges', 'kg m-1 s-1',  &
       &         'storage field for horizontal mass flux at edges', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'mass_fl_e_sv', p_diag%mass_fl_e_sv,             &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, lrestart=.FALSE.,                            &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%mass_fl_e_sv)


    ! rho_ic       p_diag%rho_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('density', 'kg m-3', 'density at half level', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'rho_ic', p_diag%rho_ic,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE.,                           &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%rho_ic)


    ! w_concorr_c  p_diag%w_concorr_c(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('contravariant_vertical_correction', 'm s-1',            &
      &                   'contravariant vertical correction', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'w_concorr_c', p_diag%w_concorr_c,                  &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE.,                           &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%w_concorr_c)


    ! theta_v_ic   p_diag%theta_v_ic(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_potential_temperature_at_half_levels', 'K',&
      &                   'virtual potential temperature at half levels', datatype_flt)
    grib2_desc = grib2_var( 0, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'theta_v_ic', p_diag%theta_v_ic,                    &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE.,                           &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%theta_v_ic)


    ! vn_ie        p_diag%vn_ie(nproma,nlevp1,nblks_e)
    !
    cf_desc    = t_cf_var('normal_wind_at_half_level', 'm s-1',                  &
      &                   'normal wind at half level', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 34, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'vn_ie', p_diag%vn_ie,                         &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & ldims=shape3d_ehalf, lrestart=.FALSE.,                         &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%vn_ie)


    ! ddt_vn_apc_pc   p_diag%ddt_vn_apc_pc(nproma,nlev,nblks_e,n_timlevs)
    cf_desc    = t_cf_var('advective+coriolis_normal_wind_tendency__predictor_corrector', 'm s-2',    &
      &                   'advective+coriolis normal wind tendency, predictor/corrector', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 201, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'ddt_vn_apc_pc', p_diag%ddt_vn_apc_pc,            &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=shape4d_entl ,&
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,          &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%ddt_vn_apc_pc)

    ALLOCATE(p_diag%ddt_vn_apc_pc_ptr(n_timlevs))
    DO jt =1,n_timlevs
      suffix = get_timelevel_string(jt)
      CALL add_ref( p_diag_list, 'ddt_vn_apc_pc',                                &
                  & 'ddt_vn_apc_pc'//suffix, p_diag%ddt_vn_apc_pc_ptr(jt)%p_3d,  &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                  & t_cf_var('ddt_vn_apc_pc'//suffix, 'm s-2','', datatype_flt), &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_EDGE),&
                  & ref_idx=jt,                                                  &
                  & ldims=shape3d_e, lrestart=.FALSE. )
    ENDDO


    IF (var_in_output%ddt_vn_adv .OR. var_in_output%ddt_ua_adv .OR. var_in_output%ddt_va_adv .OR. &
      & var_in_output%ddt_vn_cor .OR. var_in_output%ddt_ua_cor .OR. var_in_output%ddt_va_cor) THEN
      ! ddt_vn_cor_pc   p_diag%ddt_vn_cor_pc(nproma,nlev,nblks_e,n_timlevs)
      cf_desc    = t_cf_var('coriolis_normal_wind_tendency__predictor_corrector', 'm s-2',    &
        &                   'coriolis normal wind tendency, predictor/corrector', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 201, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'ddt_vn_cor_pc', p_diag%ddt_vn_cor_pc,            &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,     &
                  & ldims=shape4d_entl ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,          &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%ddt_vn_cor_pc)

      ALLOCATE(p_diag%ddt_vn_cor_pc_ptr(n_timlevs))
      DO jt =1,n_timlevs
        suffix = get_timelevel_string(jt)
        CALL add_ref( p_diag_list, 'ddt_vn_cor_pc',                                &
                    & 'ddt_vn_cor_pc'//suffix, p_diag%ddt_vn_cor_pc_ptr(jt)%p_3d,  &
                    & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                        &
                    & t_cf_var('ddt_vn_cor_pc'//suffix, 'm s-2','', datatype_flt), &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_EDGE),&
                    & ref_idx=jt,                                                  &
                    & ldims=shape3d_e, lrestart=.FALSE. )
      ENDDO
    END IF


    ! ddt_w_adv_pc    p_diag%ddt_w_adv_pc(nproma,nlevp1,nblks_c,n_timlevs)
    ! *** needs to be saved for restart (TL nnow) ***
    cf_desc    = t_cf_var('advective_vertical_wind_tendency__predictor_corrector', 'm s-2',    &
      &                   'advective vertical wind tendency, predictor/corrector', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'ddt_w_adv_pc', p_diag%ddt_w_adv_pc,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & ldims=shape4d_chalfntl,                                        &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,          &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%ddt_w_adv_pc)

    ALLOCATE(p_diag%ddt_w_adv_pc_ptr(n_timlevs))
    DO jt =1,n_timlevs
      suffix = get_timelevel_string(jt)
      CALL add_ref( p_diag_list, 'ddt_w_adv_pc',                                  &
                  & 'ddt_w_adv_pc'//suffix, p_diag%ddt_w_adv_pc_ptr(jt)%p_3d,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                    &
                  & t_cf_var('ddt_w_adv_pc'//suffix, 'm s-2','', datatype_flt),   &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                  & ref_idx=jt,                                                   &
                  & ldims=shape3d_chalf, lrestart=.TRUE. )
    ENDDO



    ! airmass_now   p_diag%airmass_now(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('airmass_now', 'kg m-2',&
      &                   'mass of air in layer at physics time step now', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'airmass_now', p_diag%airmass_now,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_c, loutput=.FALSE., lrestart=.FALSE.,              &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%airmass_now)


    ! airmass_new   p_diag%airmass_new(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('airmass_new', 'kg m-2',&
      &                   'mass of air in layer at physics time step new', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'airmass_new', p_diag%airmass_new,                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_c, loutput=.FALSE., lrestart=.FALSE.,              &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%airmass_new)
      

    ! grf_tend_vn  p_diag%grf_tend_vn(nproma,nlev,nblks_e)
    !
    ! restart needed for boundary nudging in case l_limited_area
    lrestart = p_patch%id == 1 .AND. l_limited_area

    cf_desc    = t_cf_var('normal_wind_tendency', 'm s-2',                       &
      &                   'normal wind tendency (grid refinement)', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'grf_tend_vn', p_diag%grf_tend_vn,             &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=shape3d_e, lrestart=lrestart,                            &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%grf_tend_vn)


    ! grf_tend_mflx  p_diag%grf_tend_mflx(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('normal_mass_flux_tendency', 'kg m-2 s-2',             &
      &                   'normal mass flux tendency (grid refinement)', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'grf_tend_mflx', p_diag%grf_tend_mflx,         &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=shape3d_e, lrestart=.FALSE.,                             &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%grf_tend_mflx)


    ! grf_tend_w  p_diag%grf_tend_w(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('vertical_wind_tendency', 'm s-2',                     &
      &                   'vertical wind tendency (grid refinement)', datatype_flt)
    grib2_desc = grib2_var( 0, 2, 214, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'grf_tend_w', p_diag%grf_tend_w,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                & ldims=shape3d_chalf, lrestart=.FALSE.,                         &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%grf_tend_w)

    ! grf_tend_rho   p_diag%grf_tend_rho(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('density_tendency', 'kg m-3 s-1',                      &
      &                   'density tendency (grid refinement)', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'grf_tend_rho', p_diag%grf_tend_rho,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=shape3d_c, lrestart=.FALSE.,                             &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%grf_tend_rho)

    ! grf_tend_thv   p_diag%grf_tend_thv(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('virtual_potential_temperature_tendency', 'K s-1',     &
      &                   'virtual potential temperature tendency (grid refinement)', &
      &                   datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'grf_tend_thv', p_diag%grf_tend_thv,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                & ldims=shape3d_c, lrestart=.FALSE.,                             &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%grf_tend_thv)


    ! Storage fields for vertical nesting; for vn, the middle index (2) addresses 
    ! the field and its temporal tendency

    ! vn_ie_int   p_diag%vn_ie_int(nproma,2,nblks_e)
    !
    cf_desc    = t_cf_var('vn_ie_int', 'm s-1',                               &
      &                   'normal velocity and tendency at parent interface level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'vn_ie_int', p_diag%vn_ie_int,                 &
                & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e2, lrestart=.FALSE.,                         &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%vn_ie_int)


    ! vn_ie_ubc   p_diag%vn_ie_ubc(nproma,2,nblks_e)
    !
    cf_desc    = t_cf_var('vn_ie_ubc', 'm s-1',                               &
      &                   'normal velocity and tendency at child upper boundary', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_diag_list, 'vn_ie_ubc', p_diag%vn_ie_ubc,                 &
                & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e2, lrestart=.FALSE.,                         &
                & lopenacc = .TRUE.  )
    __acc_attach(p_diag%vn_ie_ubc)


    ! w_int       p_diag%w_int(nproma,nblks_c,ndyn_substeps_max+2)
    !
    cf_desc    = t_cf_var('w_int', 'm s-1',                                   &
      &                   'vertical velocity at parent interface level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'w_int', p_diag%w_int,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcp2, lrestart=.FALSE., loutput=.FALSE.,     &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%w_int)


    ! w_ubc       p_diag%w_ubc(nproma,nblks_c,2)
    !
    cf_desc    = t_cf_var('w_ubc', 'm s-1',                                   &
      &                   'vertical velocity and tendency at child upper boundary', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'w_ubc', p_diag%w_ubc,                         &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcc, lrestart=.FALSE.,                       &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%w_ubc)


    ! theta_v_ic_int    p_diag%theta_v_ic_int(nproma,nblks_c,ndyn_substeps_max+2)
    !
    cf_desc    = t_cf_var('theta_v_ic_int', 'K',                              &
      &                   'potential temperature at parent interface level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'theta_v_ic_int', p_diag%theta_v_ic_int,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcp2, lrestart=.FALSE., loutput=.FALSE.,     &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%theta_v_ic_int)


    ! theta_v_ic_ubc    p_diag%theta_v_ic_ubc(nproma,nblks_c,2)
    !
    cf_desc    = t_cf_var('theta_v_ic_ubc', 'K',                              &
      &                   'potential temperature and tendency at child upper boundary', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'theta_v_ic_ubc', p_diag%theta_v_ic_ubc,       &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcc, lrestart=.FALSE., loutput=.TRUE.,       &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%theta_v_ic_ubc)


    ! rho_ic_int    p_diag%rho_ic_int(nproma,nblks_c,ndyn_substeps_max+2)
    !
    cf_desc    = t_cf_var('rho_ic_int', 'K',                                  &
      &                   'density at parent interface level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'rho_ic_int', p_diag%rho_ic_int,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcp2, lrestart=.FALSE., loutput=.FALSE.,     &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%rho_ic_int)


    ! rho_ic_ubc    p_diag%rho_ic_ubc(nproma,nblks_c,2)
    !
    cf_desc    = t_cf_var('rho_ic_ubc', 'K',                                  &
      &                   'density and tendency at child upper boundary', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'rho_ic_ubc', p_diag%rho_ic_ubc,               &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcc, lrestart=.FALSE., loutput=.TRUE.,       &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%rho_ic_ubc)


    ! mflx_ic_int  p_diag%mflx_ic_int(nproma,nblks_c,ndyn_substeps_max+2)
    !
    cf_desc    = t_cf_var('mflx_ic_int', 'kg m-3',                            &
      &                   'mass flux at parent interface level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'mflx_ic_int', p_diag%mflx_ic_int,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcp2, lrestart=.FALSE., loutput=.FALSE.,     &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%mflx_ic_int)


    ! mflx_ic_ubc  p_diag%mflx_ic_ubc(nproma,nblks_c,2)
    !
    cf_desc    = t_cf_var('mflx_ic_ubc', 'kg m-3',                            &
      &                   'mass flux and tendency at child upper boundary', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_diag_list, 'mflx_ic_ubc', p_diag%mflx_ic_ubc,             &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_ubcc, lrestart=.FALSE., loutput=.TRUE.,       &
                & lopenacc = .TRUE. )
    __acc_attach(p_diag%mflx_ic_ubc)


    !
    ! tracers
    !
    IF ( ntracer >0 ) THEN

      ! grf_tend_tracer   p_diag%grf_tend_tracer(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('tracer_tendency', 'kg kg-1 s-1',                   &
        &                   'tracer tendency for grid refinement', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'grf_tend_tracer', p_diag%grf_tend_tracer,        &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                  & ldims=shape4d_c ,                                              &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,          &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%grf_tend_tracer)

      ALLOCATE(p_diag%ddt_grf_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I3.3)')jt
        CALL add_ref( p_diag_list, 'grf_tend_tracer',                              &
                    & 'ddt_grf_q'//ctrc, p_diag%ddt_grf_trc_ptr(jt)%p_3d,          &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                        &
                    & t_cf_var('ddt_grf_q'//ctrc, 'kg kg-1 s**-1','',datatype_flt),&
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                    & ref_idx=jt,                                                  &
                    & ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO


      ! hfl_tracer   p_diag%hfl_tracer(nproma,nlev,nblks_e,ntracer)
      !
      cf_desc    = t_cf_var('horizontal tracer flux', 'kg m-1 s-1',               &
        &                   'horizontal tracer flux', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'hfl_tracer', p_diag%hfl_tracer,                 &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & ldims=shape4d_e ,&
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%hfl_tracer)

      ALLOCATE(p_diag%hfl_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I3.3)')jt
        CALL add_ref( p_diag_list, 'hfl_tracer',                                    &
                    & 'hfl_q'//ctrc, p_diag%hfl_trc_ptr(jt)%p_3d,                   &
                    & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE,                         &
                    & t_cf_var('hfl_q'//ctrc, 'kg m-1 s-1','', datatype_flt),       &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE),&
                    & ref_idx=jt,                                                   &
                    & ldims=shape3d_e, lrestart=.FALSE. )
      ENDDO


      ! vfl_tracer   p_diag%vfl_tracer(nproma,nlevp1,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('vertical_tracer_flux', 'kg m-1 s-1',                    &
        &                   'vertical tracer flux', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'vfl_tracer', p_diag%vfl_tracer,                    &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                  & ldims=shape4d_chalf ,                                            &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,            &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%vfl_tracer)

      ALLOCATE(p_diag%vfl_trc_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I3.3)')jt
        CALL add_ref( p_diag_list, 'vfl_tracer',                                  &
                    & 'vfl_q'//ctrc, p_diag%vfl_trc_ptr(jt)%p_3d,                 &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF,                  &
                    & t_cf_var('vfl_q'//ctrc, 'kg m-1 s-1','', datatype_flt),     &
                    & grib2_var(255,255, 255, ibits, GRID_UNSTRUCTURED,GRID_CELL),&
                    & ref_idx=jt,                                                 &
                    & ldims=shape3d_chalf, lrestart=.FALSE. )
      ENDDO


      ! ddt_tracer_adv   p_diag%ddt_tracer_adv(nproma,nlev,nblks_c,ntracer)
      !
      cf_desc    = t_cf_var('advective tracer tendency', 'kg kg-1 s-1',         &
        &                   'advective tracer tendency', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'ddt_tracer_adv', p_diag%ddt_tracer_adv,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape4d_c,                                            &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       &
                  & in_group=groups("iau_init_vars"),                           &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%ddt_tracer_adv)

      ALLOCATE(p_diag%ddt_trc_adv_ptr(ntracer))
      DO jt =1,ntracer
        WRITE(ctrc,'(I3.3)')jt
        CALL add_ref( p_diag_list, 'ddt_tracer_adv',                                &
                    & 'ddt_adv_q'//ctrc, p_diag%ddt_trc_adv_ptr(jt)%p_3d,           &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                            &
                    & t_cf_var('ddt_adv_q'//ctrc, 'kg kg-1 s-1','', datatype_flt),  &
                    & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                    & ref_idx=jt,                                                   &
                    & ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO


      ! Vertical integrals of mass-related tracers,  tracer_vi(nproma,nblks_c,ntracer)
      cf_desc    = t_cf_var('tracer_vi', '', 'tracer_vi', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'tracer_vi', p_diag%tracer_vi,                  &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,     &
                  & ldims=(/nproma, nblks_c, ntracer/), lrestart=.FALSE.,        &
                  & loutput=.FALSE., lcontainer=.TRUE.,                          &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%tracer_vi)


      ALLOCATE(p_diag%tracer_vi_ptr(ntracer))

      ! Q1 vertical integral: tqv(nproma,nblks_c)
      IF ( iqv /= 0 ) THEN
        cf_desc    = t_cf_var('tqv', 'kg m-2', 'total column integrated water vapour', &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 64, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqv',                             &
                    & p_diag%tracer_vi_ptr(iqv)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqv,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      END IF

      ! Q2 vertical integral: tqc(nproma,nblks_c)
      IF ( iqc /= 0 ) THEN
        cf_desc    = t_cf_var('tqc', 'kg m-2', 'total column integrated cloud water', &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 69, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqc',                             &
                    & p_diag%tracer_vi_ptr(iqc)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqc,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      END IF

      ! Q3 vertical integral: tqi(nproma,nblks_c)
      IF ( iqi /= 0 ) THEN
        cf_desc    = t_cf_var('tqi', 'kg m-2', 'total column integrated cloud ice', &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 70, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqi',                             &
                    & p_diag%tracer_vi_ptr(iqi)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqi,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      END IF

      IF ( iqr /= 0 ) THEN
      ! IF ( iqm_max >= 4 ) THEN
        ! Q4 vertical integral: tqr(nproma,nblks_c)
        cf_desc    = t_cf_var('tqr', 'kg m-2', 'total column integrated rain',     &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 45, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqr',                             &
                    & p_diag%tracer_vi_ptr(iqr)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqr,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      ENDIF ! iqr /= 0

      IF ( iqs /= 0 ) THEN
      ! IF ( iqm_max >= 5 ) THEN
        ! Q5 vertical integral: tqs(nproma,nblks_c)
        cf_desc    = t_cf_var('tqs', 'kg m-2', 'total column integrated snow',     &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 46, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqs',                             &
                    & p_diag%tracer_vi_ptr(iqs)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqs,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      ENDIF  ! iqs /= 0


      IF ( iqg /= 0 ) THEN
      ! IF ( ANY((/2,4,5,6/) == atm_phy_nwp_config(p_patch%id)%inwp_gscp ) ) THEN
        !
        ! Q6 vertical integral: tqg(nproma,nblks_c)
        cf_desc    = t_cf_var('tqg', 'kg m-2', 'total column integrated graupel',  &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 74, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqg',                             &
                    & p_diag%tracer_vi_ptr(iqg)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqg,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      ENDIF

      ! Note that hail is only taken into account by schemes 4, 5, 6 and 7
      !
      IF ( iqh /= 0 ) THEN
      !   IF ( ANY((/4,5,6/) == atm_phy_nwp_config(p_patch%id)%inwp_gscp ) ) THEN
        ! Q7 vertical integral: tqh(nproma,nblks_c)
        cf_desc    = t_cf_var('tqh', 'kg m-2', 'total column integrated hail',     &
          &          datatype_flt)
        grib2_desc = grib2_var( 0, 1, 72, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_ref( p_diag_list, 'tracer_vi', 'tqh',                             &
                    & p_diag%tracer_vi_ptr(iqh)%p_2d,                              &
                    & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
                    & cf_desc, grib2_desc,                                         &
                    & ref_idx=iqh,                                                 &
                    & ldims=shape2d_c, lrestart=.FALSE.)
      ENDIF

    ENDIF  !  ntracer >0


    !
    ! (Transformed) Analysis increments
    !
    IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
      ! vn_incr   p_diag%vn_incr(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('vn_incr', ' ',                   &
        &                   'vn increment from DA', datatype_flt)
      grib2_desc = grib2_var( 0, 2, 34, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_diag_list, 'vn_incr', p_diag%vn_incr,                     &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e,                                            &
                  & lrestart=.FALSE., loutput=.FALSE.,                          &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%vn_incr)


      ! exner_incr   p_diag%exner_incr(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('exner_incr', ' ',                   &
        &                   'exner increment from DA', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 26, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'exner_incr', p_diag%exner_incr,               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                           &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%exner_incr)


      ! rho_incr   p_diag%rho_incr(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('rho_incr', 'kg m-3',                   &
        &                   'density increment from DA', datatype_flt)
      grib2_desc = grib2_var( 0, 3, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'rho_incr', p_diag%rho_incr,                   &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                           &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%rho_incr)


      ! rhov_incr  p_diag%rhov_incr(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('rhov_incr', 'kg m-3',                   &
        &                   'partial density of water vapour increment from DA', datatype_flt)
      grib2_desc = grib2_var(  0, 1, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'rhov_incr', p_diag%rhov_incr,                 &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & lrestart=.FALSE., loutput=.TRUE.,                           &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%rhov_incr)

      IF (qcana_mode > 0) THEN
        ! rhoc_incr  p_diag%rhoc_incr(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('rhoc_incr', 'kg m-3',                   &
          &                   'partial density of cloud water increment from DA', datatype_flt)
        grib2_desc = grib2_var( 0, 6, 38, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_diag_list, 'rhoc_incr', p_diag%rhoc_incr,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                    & ldims=shape3d_c,                                               &
                    & lrestart=.FALSE., loutput=.TRUE.,                              &
                    & lopenacc = .TRUE. )
        __acc_attach(p_diag%rhoc_incr)
      ENDIF

      IF (qiana_mode > 0) THEN
        ! rhoi_incr  p_diag%rhoi_incr(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('rhoi_incr', ' ',                   &
          &                   'partial density of cloud ice increment from DA', datatype_flt)
        grib2_desc = grib2_var( 0, 6, 39, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_diag_list, 'rhoi_incr', p_diag%rhoi_incr,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                    & ldims=shape3d_c, &
                    & lrestart=.FALSE., loutput=.TRUE.,                              &
                    & lopenacc = .TRUE. )
        __acc_attach(p_diag%rhoi_incr)
      END IF
        
      IF (qrsgana_mode > 0) THEN
        ! rhor_incr  p_diag%rhor_incr(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('rhor_incr', ' ',                   &
          &                   'partial density of rain increment from DA', datatype_flt)
        grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_diag_list, 'rhor_incr', p_diag%rhor_incr,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                    & ldims=shape3d_c, &
                    & lrestart=.FALSE., loutput=.TRUE.,                              &
                    & lopenacc = .TRUE. )
        __acc_attach(p_diag%rhor_incr)

        ! rhos_incr  p_diag%rhos_incr(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('rhos_incr', ' ',                   &
          &                   'partial density of snow increment from DA', datatype_flt)
        grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_diag_list, 'rhos_incr', p_diag%rhos_incr,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                    & ldims=shape3d_c, &
                    & lrestart=.FALSE., loutput=.TRUE.,                              &
                    & lopenacc = .TRUE. )
        __acc_attach(p_diag%rhos_incr)


        ! rhog_incr  p_diag%rhog_incr(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('rhog_incr', ' ',                   &
          &                   'partial density of graupel increment from DA', datatype_flt)
        grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_diag_list, 'rhog_incr', p_diag%rhog_incr,                    &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                    & ldims=shape3d_c, &
                    & lrestart=.FALSE., loutput=.TRUE.,                              &
                    & lopenacc = .TRUE. )
        __acc_attach(p_diag%rhog_incr)

        IF ( atm_phy_nwp_config(p_patch%id)%l2moment ) THEN

          ! rhoh_incr  p_diag%rhoh_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhoh_incr', ' ',                   &
            &                   'partial density of hail increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhoh_incr', p_diag%rhoh_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhoh_incr)
        END IF
      END IF

      IF (atm_phy_nwp_config(p_patch%id)%l2moment) THEN
        IF (qcana_mode > 0) THEN
          ! rhonc_incr  p_diag%rhonc_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhonc_incr', ' ',                   &
               &                'number density of cloud water increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhonc_incr', p_diag%rhonc_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhonc_incr)
        END IF
        IF (qiana_mode > 0) THEN
          ! rhoni_incr  p_diag%rhoni_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhoni_incr', ' ',                   &
               &                'number density of cloud ice increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhoni_incr', p_diag%rhoni_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhoni_incr)
        END IF
        IF (qrsgana_mode > 0) THEN
          ! rhonr_incr  p_diag%rhonr_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhonr_incr', ' ',                   &
               &                'number density of rain increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhonr_incr', p_diag%rhonr_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhonr_incr)

          ! rhons_incr  p_diag%rhons_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhons_incr', ' ',                   &
               &                'number density of snow increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhons_incr', p_diag%rhons_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhons_incr)
          
          ! rhong_incr  p_diag%rhong_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhong_incr', ' ',                   &
               &                'number density of graupel increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhong_incr', p_diag%rhong_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhong_incr)

          ! rhonh_incr  p_diag%rhonh_incr(nproma,nlev,nblks_c)
          !
          cf_desc    = t_cf_var('rhonh_incr', ' ',                   &
               &                'number density of hail increment from DA', datatype_flt)
          grib2_desc = grib2_var( 255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
          CALL add_var( p_diag_list, 'rhonh_incr', p_diag%rhonh_incr,                    &
                      & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,     &
                      & ldims=shape3d_c, &
                      & lrestart=.FALSE., loutput=.TRUE.,                              &
                      & lopenacc = .TRUE. )
          __acc_attach(p_diag%rhonh_incr)

        END IF
      END IF
    ENDIF  ! init_mode = MODE_IAU, MODE_IAU_OLD

    ! From a logical point of view, the following two fields make sense only in combination with an IAU cycle,
    ! but allocating them anyway allows using the analysis interpolation mode without changing the namelists  
    !
    IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
      !  Filtered T2M bias
      cf_desc    = t_cf_var('t2m_bias', 'K', 'Filtered T2M bias', datatype_flt)
      grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
                 + t_grib2_int_key("typeOfGeneratingProcess", 206)
      CALL add_var( p_diag_list, 't2m_bias', p_diag%t2m_bias,                         &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_2M, cf_desc, grib2_desc,        &
        &           ldims=shape2d_c, lrestart=.true.,                                 &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in") )
    ENDIF

    IF (icpl_da_sfcevap >= 2) THEN
      !  Time-filtered near-surface level RH increment from data assimilation
      cf_desc    = t_cf_var('rh_avginc', '1', 'Filtered RH increment', datatype_flt)
      grib2_desc = grib2_var(0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
                 + t_grib2_int_key("typeOfGeneratingProcess", 207)     &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)      &
                 + t_grib2_int_key("scaledValueOfFirstFixedSurface", 20)
      CALL add_var( p_diag_list, 'rh_avginc', p_diag%rh_avginc,                            &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,            &
        &           ldims=shape2d_c, lrestart=.true.,                                      &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in"), & 
        &           lopenacc = .TRUE. )
      __acc_attach(p_diag%rh_avginc)
    ENDIF

    IF (icpl_da_sfcevap >= 3) THEN
      !  Time-filtered near-surface level T increment from data assimilation
      cf_desc    = t_cf_var('t_avginc', 'K', 'Filtered T increment', datatype_flt)
      grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
                 + t_grib2_int_key("typeOfGeneratingProcess", 207)     &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)      &
                 + t_grib2_int_key("scaledValueOfFirstFixedSurface", 20)
      CALL add_var( p_diag_list, 't_avginc', p_diag%t_avginc,                       &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,     &
        &           ldims=shape2d_c, lrestart=.true.,                               &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in") )
    ENDIF

    IF (icpl_da_sfcfric >= 1) THEN
      !  Time-filtered near-surface level wind speed increment from data assimilation
      cf_desc    = t_cf_var('vabs_avginc', 'm/s', 'Filtered wind speed increment', datatype_flt)
      grib2_desc = grib2_var(0, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
                 + t_grib2_int_key("typeOfGeneratingProcess", 207)     &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)      &
                 + t_grib2_int_key("scaledValueOfFirstFixedSurface", 20)
      CALL add_var( p_diag_list, 'vabs_avginc', p_diag%vabs_avginc,                 &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,     &
        &           ldims=shape2d_c, lrestart=.true.,                               &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in") )
    ENDIF

    IF (icpl_da_skinc >= 1) THEN
      !  Time-filtered near-surface level T increment from data assimilation, weighted with COS(local time)
      cf_desc    = t_cf_var('t_wgt_avginc', 'K', 'Weighted filtered T increment', datatype_flt)
      grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
                 + t_grib2_int_key("typeOfGeneratingProcess", 208)     &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)      &
                 + t_grib2_int_key("scaledValueOfFirstFixedSurface", 20)
      CALL add_var( p_diag_list, 't_wgt_avginc', p_diag%t_wgt_avginc,               &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,     &
        &           ldims=shape2d_c, lrestart=.true.,                               &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in") )

    ENDIF

    IF (icpl_da_sfcevap >= 5) THEN
      !  Time-filtered near-surface level T increment from data assimilation, weighted with peak in afternoon
      cf_desc    = t_cf_var('t_daywgt_avginc', 'K', 'Weighted filtered T increment', datatype_flt)
      grib2_desc = grib2_var(0, 210, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)&
                 + t_grib2_int_key("typeOfGeneratingProcess", 207)         &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)          &
                 + t_grib2_int_key("typeOfFirstFixedSurface", 103)
      CALL add_var( p_diag_list, 't_daywgt_avginc', p_diag%t_daywgt_avginc,         &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,     &
        &           ldims=shape2d_c, lrestart=.true.,                               &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in") )

      !  Time-filtered near-surface level RH increment from data assimilation, weighted with peak in afternoon
      cf_desc    = t_cf_var('rh_daywgt_avginc', '1', 'Weighted filtered RH increment', datatype_flt)
      grib2_desc = grib2_var(0, 210, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)&
                 + t_grib2_int_key("typeOfGeneratingProcess", 207)         &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)          &
                 + t_grib2_int_key("typeOfFirstFixedSurface", 103)
      CALL add_var( p_diag_list, 'rh_daywgt_avginc', p_diag%rh_daywgt_avginc,       &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,     &
        &           ldims=shape2d_c, lrestart=.true.,                               &
        &           in_group=groups("mode_iau_fg_in","mode_dwd_fg_in","mode_combined_in") )
    ENDIF

    IF (latbc_config%fac_latbc_presbiascor > 0._wp) THEN
      !  Time-filtered near-surface level pressure increment from data assimilation
      cf_desc    = t_cf_var('p_avginc', 'Pa', 'Filtered P increment', datatype_flt)
      grib2_desc = grib2_var(0, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL) &
                 + t_grib2_int_key("typeOfGeneratingProcess", 207)         &
                 + t_grib2_int_key("typeOfSecondFixedSurface", 1)          &
                 + t_grib2_int_key("scaledValueOfFirstFixedSurface", 20)
      CALL add_var( p_diag_list, 'p_avginc', p_diag%p_avginc,                       &
        &           GRID_UNSTRUCTURED_CELL, ZA_HEIGHT_10M, cf_desc, grib2_desc,     &
        &           ldims=shape2d_c, lrestart=.true.,                               &
        &           in_group=groups("mode_iau_fg_in") )
    ENDIF


    !----------------------
    ! optional diagnostics
    !----------------------

    ! pres_msl           p_diag%pres_msl(nproma,nblks_c)
    !
    ! Note: This task is registered for the post-processing scheduler
    !        which takes care of the regular update.
    !
    IF (var_in_output%pres_msl .OR. is_coupled_to_ocean()) THEN
      cf_desc    = t_cf_var('mean sea level pressure', 'Pa', &
        &                   'mean sea level pressure', datatype_flt)
      grib2_desc = grib2_var(0, 3, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'pres_msl', p_diag%pres_msl,                     &
        &           GRID_UNSTRUCTURED_CELL, ZA_MEANSEA, cf_desc, grib2_desc,      &
        &           ldims=shape2d_c, lrestart=.FALSE.,                            &
        &           l_pp_scheduler_task=TASK_INTP_MSL,                            &
        &           lopenacc = .TRUE. )
      __acc_attach(p_diag%pres_msl)
    END IF

    ! vertical velocity ( omega=dp/dt ) 
    !
    ! Note: This task is registered for the post-processing scheduler
    !       which takes care of the regular update:
    ! 
    IF (var_in_output%omega) THEN
      cf_desc    = t_cf_var('omega', 'Pa s-1', 'vertical velocity', datatype_flt)
      grib2_desc = grib2_var(0, 2, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list,                                                     &
                    & "omega", p_diag%omega,                                         &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                          &
                    & cf_desc, grib2_desc,                                           &
                    & ldims=shape3d_c,                                               &
                    & vert_interp=create_vert_interp_metadata(                       &
                    &             vert_intp_type=vintp_types("P","Z","I"),           &
                    &             vert_intp_method=VINTP_METHOD_LIN,                 &
                    &             l_loglin=.FALSE., l_extrapol=.FALSE.),             &
                    & in_group=groups("atmo_derived_vars"),                          &
                    & l_pp_scheduler_task=TASK_COMPUTE_OMEGA, lrestart=.FALSE.,      &
                    & lopenacc = .TRUE. )
      __acc_attach(p_diag%omega)
    END IF

    !---------------------- End of optional diagnostics ----------------------


    IF(inextra_2d > 0) THEN


      ! extra_2d   p_diag%extra_2d(nproma,nblks_c,inextra_2d)
      !
      cf_desc    = t_cf_var('extra_field_2D', '-', 'extra field 2D', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'extra_2d', p_diag%extra_2d,                   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & initval=0._wp,                                              &
                  & lcontainer=.TRUE., ldims=shape2d_extra, lrestart=.FALSE.,   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%extra_2d)

      ALLOCATE(p_diag%extra_2d_ptr(inextra_2d))
      DO jt =1,inextra_2d
        ! GRIB2: extra fields are encoded as "DUMMY_x" (where x=jt)
        grib2_desc = grib2_var( 0, 254, jt, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        WRITE(ctrc,'(I2)')jt
        CALL add_ref( p_diag_list, 'extra_2d', 'extra_2d'//TRIM(ADJUSTL(ctrc)), &
          &           p_diag%extra_2d_ptr(jt)%p_2d,                             &
          &           GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                       &
          &           cf_desc, grib2_desc,                                      &
          &           ref_idx=jt,                                               &
          &           ldims=shape2d_c, lrestart=.FALSE. )
      ENDDO
    ENDIF


    IF(inextra_3d > 0) THEN

      ! extra_3d   p_diag%extra_3d(nproma,nlev,nblks_c,inextra_3d)
      !
      cf_desc    = t_cf_var('extra_fields_3D', '-', 'extra fields 3D', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_diag_list, 'extra_3d', p_diag%extra_3d,                   &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & initval=0._wp, ldims=shape3d_extra,                         &
                  & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       & 
                  & lopenacc = .TRUE. )
      __acc_attach(p_diag%extra_3d)

      ALLOCATE(p_diag%extra_3d_ptr(inextra_3d))
      DO jt =1,inextra_3d
        ! GRIB2: extra fields are encoded as "DUMMY_x" (where x=jt)
        grib2_desc = grib2_var( 0, 254, jt, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        WRITE(ctrc,'(I2)')jt
        CALL add_ref( p_diag_list, 'extra_3d', 'extra_3d'//TRIM(ADJUSTL(ctrc)), &
          &           p_diag%extra_3d_ptr(jt)%p_3d,                             &
          &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                     &
          &           cf_desc, grib2_desc,                                      &
          &           ref_idx=jt,                                               &
          &           ldims=shape3d_c, lrestart=.FALSE. )
      ENDDO
    ENDIF
    
  END SUBROUTINE new_nh_state_diag_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of reference state.
  !!
  !! Initialization of components with zero.
  !!
  SUBROUTINE new_nh_state_ref_list ( p_patch, p_ref, p_ref_list,  &
    &                                 listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_ref),  INTENT(INOUT)   :: &  !< reference state
      &  p_ref 

    TYPE(t_var_list_ptr), INTENT(INOUT)   :: &  !< reference state list
      &  p_ref_list
    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, nblks_e    !< number of cell/edge blocks to allocate

    INTEGER :: nlev, nlevp1  !< number of vertical full/half levels

    INTEGER :: shape3d_e(3), shape3d_chalf(3)
 
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    !--------------------------------------------------------------

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(p_ref%vn_ref, p_ref%w_ref)


    ! We ONLY need the variables IN p_ref for test runs.
    IF(.NOT. ltestcase ) RETURN     ! WS: this is bad coding!  The test needs to be in front of the subroutine call.


    !determine size of arrays
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ! predefined array shapes
    shape3d_e     = (/nproma, nlev   , nblks_e    /)
    shape3d_chalf = (/nproma, nlevp1 , nblks_c    /)


    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_ref_list, TRIM(listname), patch_id=p_patch%id, &
      &          lrestart=.FALSE., model_type=get_my_process_name())

    ! vn_ref     p_ref%vn_ref(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('normal_velocity', 'm s-1', 'velocity normal to edge', datatype_flt)
    grib2_desc = grib2_var(0, 2, 34, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_ref_list, 'vn_ref', p_ref%vn_ref,                              &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,       &   
                & ldims=shape3d_e, lrestart=.FALSE.,  loutput=.FALSE.,             &
                & isteptype=TSTEP_CONSTANT,                                        &
                & lopenacc = .TRUE. )
    __acc_attach(p_ref%vn_ref)


    ! w_ref      p_ref%w_ref(nproma,nlev+1,nblks_c)
    !
    cf_desc    = t_cf_var('upward_air_velocity', 'm s-1', 'Vertical velocity', datatype_flt)
    grib2_desc = grib2_var(0, 2, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_ref_list, 'w_ref', p_ref%w_ref,                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, lrestart=.FALSE., loutput=.FALSE.,          &
                & isteptype=TSTEP_CONSTANT,                                        &
                & lopenacc = .TRUE. )
    __acc_attach(p_ref%w_ref)

  END SUBROUTINE new_nh_state_ref_list


  !---------------------------------------------------------------------------
  !>
  !! Allocates all metric coefficients defined in type metrics_3d of the patch.
  !!
  SUBROUTINE new_nh_metrics_list ( p_patch, p_metrics, p_metrics_list,  &
    &                              listname )
!
    TYPE(t_patch), TARGET, INTENT(IN) :: &  !< current patch
      &  p_patch

    TYPE(t_nh_metrics),  INTENT(INOUT):: &  !< diagnostic state
      &  p_metrics 

    TYPE(t_var_list_ptr), INTENT(INOUT) :: p_metrics_list   !< diagnostic state list

    CHARACTER(len=*), INTENT(IN)      :: &  !< list name
      &  listname

    ! local variables
    CHARACTER(*), PARAMETER    :: routine = modname//"::new_nh_metrics_list"
    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e, &    !< number of edge blocks to allocate
               nblks_v       !< number of vertex blocks to allocate

    INTEGER :: nlev, nlevp1, jg

    INTEGER :: shape1d_c(1), shape1d_chalf(1), shape2d_c(2), shape2d_e(2),  &
      &        shape3d_c(3), shape3d_e(3), shape3d_v(3),                 &
      &        shape3d_chalf(3), shape3d_ehalf(3), shape2d_ccubed(3),    &
      &        shape2d_ecubed(3), shape3d_vhalf(3), shape2d_esquared(3), &
      &        shape3d_esquared(4), shape3d_e8(4)
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for selected fields
    INTEGER :: datatype_flt       !< floating point accuracy in NetCDF output
    LOGICAL :: group(MAX_GROUPS)
    INTEGER :: ic, jb, jc, i_startblk, i_endblk, i_startidx, i_endidx

    !--------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nblks_v = p_patch%nblks_v

    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    jg     = p_patch%id

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF


    ! Moved here from set_nh_metrics because we need bdy_halo_c_dim/blk/idx for add_var
    i_startblk = p_patch%cells%start_block(min_rlcell_int-1)
    i_endblk   = p_patch%cells%end_block(min_rlcell)
    ic = 0
    DO jb = i_startblk, i_endblk

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                         i_startidx, i_endidx, min_rlcell_int-1, min_rlcell)

      DO jc = i_startidx, i_endidx
        IF (p_patch%cells%refin_ctrl(jc,jb)>=1 .AND. &
          & p_patch%cells%refin_ctrl(jc,jb)<=4) THEN
          ic = ic+1
        ENDIF
      ENDDO
    ENDDO
    p_metrics%bdy_halo_c_dim = ic
    !$ACC UPDATE DEVICE(p_metrics%bdy_halo_c_dim) ASYNC(1)


    ! predefined array shapes
    shape1d_c        = (/nlev                        /)
    shape1d_chalf    = (/nlevp1                      /)
    shape2d_c        = (/nproma,          nblks_c    /)
    shape2d_e        = (/nproma,          nblks_e    /)
    shape2d_esquared = (/nproma, 2      , nblks_e    /)    
    shape2d_ccubed   = (/nproma, 3      , nblks_c    /)     
    shape2d_ecubed   = (/nproma, 3      , nblks_e    /)     
    shape3d_c        = (/nproma, nlev   , nblks_c    /)     
    shape3d_chalf    = (/nproma, nlevp1 , nblks_c    /)      
    shape3d_e        = (/nproma, nlev   , nblks_e    /)     
    shape3d_ehalf    = (/nproma, nlevp1 , nblks_e    /)     
    shape3d_esquared = (/2     , nproma , nlev   , nblks_e /)
    shape3d_e8       = (/8     , nproma , nlev   , nblks_e /)
    shape3d_v        = (/nproma, nlev   , nblks_v    /)     
    shape3d_vhalf    = (/nproma, nlevp1 , nblks_v    /)


    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(p_metrics%z_ifc, &
    &       p_metrics%z_mc, &
    &       p_metrics%ddqz_z_full, &
    &       p_metrics%geopot, &
    &       p_metrics%geopot_agl, &
    &       p_metrics%geopot_agl_ifc, &
    &       p_metrics%dgeopot_mc, &
    &       p_metrics%rayleigh_w, &
    &       p_metrics%rayleigh_vn, &
    &       p_metrics%enhfac_diffu, &
    &       p_metrics%scalfac_dd3d, &
    &       p_metrics%vwind_expl_wgt, &
    &       p_metrics%vwind_impl_wgt, &
    &       p_metrics%theta_ref_mc, &
    &       p_metrics%theta_ref_me, &
    &       p_metrics%theta_ref_ic, &
    &       p_metrics%tsfc_ref, &
    &       p_metrics%exner_ref_mc, &
    &       p_metrics%rho_ref_mc, &
    &       p_metrics%rho_ref_me, &
    &       p_metrics%zd_intcoef, &
    &       p_metrics%zd_geofac, &
    &       p_metrics%zd_e2cell, &
    &       p_metrics%zd_diffcoef, &
    &       p_metrics%inv_ddqz_z_full_e, &
    &       p_metrics%inv_ddqz_z_full_v, &
    &       p_metrics%inv_ddqz_z_half, &
    &       p_metrics%inv_ddqz_z_half_e, &
    &       p_metrics%inv_ddqz_z_half_v, &
    &       p_metrics%wgtfac_v, &
    &       p_metrics%mixing_length_sq, &
    &       p_metrics%rho_ref_corr, &
    &       p_metrics%ddxn_z_full, &
    &       p_metrics%ddxn_z_full_c, &
    &       p_metrics%ddxn_z_full_v, &
    &       p_metrics%ddxn_z_half_e, &
    &       p_metrics%ddxn_z_half_c, &
    &       p_metrics%ddxt_z_full, &
    &       p_metrics%ddxt_z_full_c, &
    &       p_metrics%ddxt_z_full_v, &
    &       p_metrics%ddxt_z_half_e, &
    &       p_metrics%ddxt_z_half_c, &
    &       p_metrics%ddxt_z_half_v, &
    &       p_metrics%ddqz_z_full_e, &
    &       p_metrics%ddqz_z_half, &
    &       p_metrics%inv_ddqz_z_full, &
    &       p_metrics%wgtfac_c, &
    &       p_metrics%wgtfac_e, &
    &       p_metrics%wgtfacq_c, &
    &       p_metrics%wgtfacq_e, &
    &       p_metrics%wgtfacq1_c, &
    &       p_metrics%wgtfacq1_e, &
    &       p_metrics%coeff_gradekin, &
    &       p_metrics%coeff1_dwdz, &
    &       p_metrics%coeff2_dwdz, &
    &       p_metrics%zdiff_gradp, &
    &       p_metrics%coeff_gradp, &
    &       p_metrics%exner_exfac, &
    &       p_metrics%d_exner_dz_ref_ic, &
    &       p_metrics%d2dexdz2_fac1_mc, &
    &       p_metrics%d2dexdz2_fac2_mc, &
    &       p_metrics%pg_exdist, &
    &       p_metrics%vertidx_gradp, &
    &       p_metrics%zd_indlist, &
    &       p_metrics%zd_blklist, &
    &       p_metrics%zd_edgeidx, &
    &       p_metrics%zd_edgeblk, &
    &       p_metrics%zd_vertidx, &
    &       p_metrics%pg_edgeidx, &
    &       p_metrics%pg_edgeblk, &
    &       p_metrics%pg_vertidx, &
    &       p_metrics%nudge_c_idx, &
    &       p_metrics%nudge_e_idx, &
    &       p_metrics%nudge_c_blk, &
    &       p_metrics%nudge_e_blk, &
    &       p_metrics%bdy_halo_c_idx, &
    &       p_metrics%bdy_halo_c_blk, &
    &       p_metrics%ovlp_halo_c_dim, &
    &       p_metrics%ovlp_halo_c_idx, &
    &       p_metrics%ovlp_halo_c_blk, &
    &       p_metrics%bdy_mflx_e_idx, &
    &       p_metrics%bdy_mflx_e_blk, &
    &       p_metrics%nudgecoeff_vert, &
    &       p_metrics%mask_prog_halo_c, &
    &       p_metrics%mask_mtnpoints, &
    &       p_metrics%mask_mtnpoints_g, &
    &       p_metrics%deepatmo_gradh_mc, &
    &       p_metrics%deepatmo_divh_mc, &
    &       p_metrics%deepatmo_vol_mc, &
    &       p_metrics%deepatmo_invr_mc, &
    &       p_metrics%deepatmo_divzU_mc, &
    &       p_metrics%deepatmo_divzL_mc, &
    &       p_metrics%deepatmo_gradh_ifc, &
    &       p_metrics%deepatmo_invr_ifc)


    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_metrics_list, TRIM(listname), patch_id=p_patch%id, &
      &          lrestart=.FALSE., model_type=get_my_process_name())

    ! geometric height at the vertical interface of cells
    ! z_ifc        p_metrics%z_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_half_level_center', 'm',         &
      &                   'geometric height at half level center', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 6, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 101)
    IF (init_mode == MODE_ICONVREMAP) THEN
      group=groups("dwd_fg_atm_vars", "LATBC_PREFETCH_VARS", "mode_dwd_fg_in",  &
        &          "mode_iniana","icon_lbc_vars")
    ELSE
      group=groups("dwd_fg_atm_vars", "LATBC_PREFETCH_VARS",                    &
        &          "mode_iniana","icon_lbc_vars")
    ENDIF
    CALL add_var( p_metrics_list, 'z_ifc', p_metrics%z_ifc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF_HHL, cf_desc, grib2_desc, &
                & ldims=shape3d_chalf,                                          & 
                & vert_interp=create_vert_interp_metadata(                      &
                &   vert_intp_type=vintp_types("P","Z","I"),                    &
                &   vert_intp_method=VINTP_METHOD_LIN_NLEVP1 ),                 &
                & in_group=group, isteptype=TSTEP_CONSTANT,                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%z_ifc)

    ! geometric height at full levels
    ! z_mc         p_metrics%z_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geometric_height_at_full_level_center', 'm',         &
      &                   'geometric height at full level center', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'z_mc', p_metrics%z_mc,                       &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &    vert_intp_type=vintp_types("P","Z","I"),                   &
                &    vert_intp_method=VINTP_METHOD_LIN ),                       &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%z_mc)

    ! slope of the terrain in normal direction (full level)
    ! ddxn_z_full  p_metrics%ddxn_z_full(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',             &
      &                   'terrain slope in normal direction', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxn_z_full', p_metrics%ddxn_z_full,         &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, loutput=.TRUE.,                              &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%ddxn_z_full)


    ! slope of the terrain in tangential direction (full level)
    ! ddxt_z_full  p_metrics%ddxt_z_full(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',         &
      &                   'terrain slope in tangential direction', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddxt_z_full', p_metrics%ddxt_z_full,         &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, loutput=.TRUE.,                              &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%ddxt_z_full)

#ifndef __NO_ICON_LES__
    IF (atm_phy_nwp_config(jg)%is_les_phy &
       & .OR. aes_vdf_config(jg)%turb==VDIFF_TURB_3DSMAGORINSKY) THEN
#else
    IF (aes_vdf_config(jg)%turb==VDIFF_TURB_3DSMAGORINSKY) THEN
#endif
      ! slope of the terrain in normal direction (half level)
      ! ddxn_z_half_e  p_metrics%ddxn_z_full(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',           &
           &                   'terrain slope in normal direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddxn_z_half_e', p_metrics%ddxn_z_half_e,   &
           & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
           & ldims=shape3d_ehalf, loutput=.TRUE.,                               &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxn_z_half_e)


      ! slope of the terrain in normal direction (half level)
      ! ddxn_z_half_c  p_metrics%ddxn_z_full(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',           &
           &                   'terrain slope in normal direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'ddxn_z_half_c', p_metrics%ddxn_z_half_c,   &
           & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
           & ldims=shape3d_chalf, loutput=.TRUE.,                               &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxn_z_half_c)


      ! slope of the terrain in normal direction (full level)
      ! ddxn_z_full_c  p_metrics%ddxn_z_full(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',           &
           &                   'terrain slope in normal direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'ddxn_z_full_c', p_metrics%ddxn_z_full_c,   &
           & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,         &
           & ldims=shape3d_c, loutput=.TRUE.,                                   &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxn_z_full_c)


      ! slope of the terrain in normal direction (full level)
      ! ddxn_z_full_v  p_metrics%ddxn_z_full(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_normal_direction', '-',           &
           &                   'terrain slope in normal direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddxn_z_full_v', p_metrics%ddxn_z_full_v,   &
           & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE, cf_desc, grib2_desc,         &
           & ldims=shape3d_v, loutput=.TRUE.,                                   &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxn_z_full_v)


      ! slope of the terrain in tangential direction (full level)
      ! ddxt_z_half_e  p_metrics%ddxt_z_full(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',       &
           &                   'terrain slope in tangential direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'ddxt_z_half_e', p_metrics%ddxt_z_half_e,   &
           & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
           & ldims=shape3d_ehalf, loutput=.TRUE.,                               &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxt_z_half_e)


      ! slope of the terrain in tangential direction (full level)
      ! ddxt_z_half_c  p_metrics%ddxt_z_full(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',       &
           &                   'terrain slope in tangential direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'ddxt_z_half_c', p_metrics%ddxt_z_half_c,   &
           & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
           & ldims=shape3d_chalf, loutput=.TRUE.,                               &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxt_z_half_c)


      ! slope of the terrain in tangential direction (full level)
      ! ddxt_z_half_v  p_metrics%ddxt_z_full(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',       &
           &                   'terrain slope in tangential direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddxt_z_half_v', p_metrics%ddxt_z_half_v,   &
           & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE_HALF, cf_desc, grib2_desc,    &
           & ldims=shape3d_vhalf, loutput=.TRUE.,                               &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxt_z_half_v)


      ! slope of the terrain in tangential direction (full level)
      ! ddxt_z_full_c  p_metrics%ddxt_z_full_c(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',       &
           &                   'terrain slope in tangential direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'ddxt_z_full_c', p_metrics%ddxt_z_full_c,   &
           & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,         &
           & ldims=shape3d_c, loutput=.TRUE.,                                   &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxt_z_full_c)


      ! slope of the terrain in tangential direction (full level)
      ! ddxt_z_full_v  p_metrics%ddxt_z_full(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('terrain_slope_in_tangential_direction', '-',       &
           &                   'terrain slope in tangential direction', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'ddxt_z_full_v', p_metrics%ddxt_z_full_v,   &
           & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE, cf_desc, grib2_desc,         &
           & ldims=shape3d_v, loutput=.TRUE.,                                   &
           & isteptype=TSTEP_CONSTANT,                                          &
           & lopenacc = .TRUE. )
      __acc_attach(p_metrics%ddxt_z_full_v)

    ENDIF  !is_les_phy

    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full_e  p_metrics%ddqz_z_full_e(nproma,nlev,nblks_e)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant (edge)', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 192, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_metrics_list, 'ddqz_z_full_e', p_metrics%ddqz_z_full_e,     &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_e, loutput=.TRUE.,                              &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%ddqz_z_full_e)


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_half  p_metrics%ddqz_z_half(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_half', p_metrics%ddqz_z_half,         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, loutput=.TRUE.,                          &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%ddqz_z_half)


    ! functional determinant of the metrics [sqrt(gamma)]
    ! ddqz_z_full  p_metrics%ddqz_z_full(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('metrics_functional_determinant', '-',                &
      &                   'metrics functional determinant', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'ddqz_z_full', p_metrics%ddqz_z_full,         &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c, loutput=.TRUE.,                              &
                & isteptype=TSTEP_CONSTANT,                                     &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%ddqz_z_full)


    ! geopotential at full level cell center
    ! geopot       p_metrics%geopot(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
       &                  'geopotential at full level cell centre', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot', p_metrics%geopot,                   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & vert_interp=create_vert_interp_metadata(                      &
                &             vert_intp_type=vintp_types("P","Z","I"),          &
                &             vert_intp_method=VINTP_METHOD_LIN,                &
                &             l_extrapol=.TRUE., l_pd_limit=.FALSE.),           &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%geopot)


    ! geopotential above groundlevel at cell center
    ! geopot_agl   p_metrics%geopot_agl(nproma,nlev  ,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl', p_metrics%geopot_agl,           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%geopot_agl)


    ! geopotential above groundlevel at cell center
    ! geopot_agl_ifc  p_metrics%geopot_agl_ifc(nproma,nlevp1,nblks_c)
    !
    cf_desc    = t_cf_var('geopotential', 'm2 s-2',                             &
      &                   'geopotential above groundlevel at cell center', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'geopot_agl_ifc', p_metrics%geopot_agl_ifc,   &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf,                                          &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%geopot_agl_ifc)


    ! difference of geopotential between the half levels
    ! dgeopot_mc   p_metrics%dgeopot_mc(nproma,nlev,nblks_c)
    !
    cf_desc    = t_cf_var('dgeopot_mc', 'm2 s-2',                               &
      &                   'geopotential difference between half levels', datatype_flt)
    grib2_desc = grib2_var( 0, 3, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'dgeopot_mc', p_metrics%dgeopot_mc,           &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                & ldims=shape3d_c,                                              &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%dgeopot_mc)

    
    ! Rayleigh damping coefficient for w
    cf_desc = t_cf_var('rayleigh_w', '-',                                      &
     &                'Rayleigh damping coefficient for w', datatype_flt)
    grib2_desc = grib2_var( 0, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'rayleigh_w', p_metrics%rayleigh_w,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc, &
        & ldims = shape1d_chalf , &
        & lopenacc = .TRUE. )
    __acc_attach(p_metrics%rayleigh_w)

    ! Rayleigh damping coefficient for vn
    cf_desc = t_cf_var('rayleigh_vn', '-',                                      &
     &                'Rayleigh damping coefficient for vn', datatype_flt)
    grib2_desc = grib2_var( 0, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'rayleigh_vn', p_metrics%rayleigh_vn,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims = shape1d_c , &
        & lopenacc = .TRUE. )
   __acc_attach(p_metrics%rayleigh_vn)

    ! Background nabla2 diffusion coefficient for upper sponge layer
   cf_desc = t_cf_var('enhfac_diffu', '-',                                      &
     &                'Background nabla2 diffusion coefficient for upper sponge layer', datatype_flt)
   grib2_desc = grib2_var( 0, 1, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
   CALL add_var( p_metrics_list, 'enhfac_diffu', p_metrics%enhfac_diffu,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims = shape1d_c , &
        & lopenacc = .TRUE. )
   __acc_attach(p_metrics%enhfac_diffu)

    ! Scaling factor for 3D divergence damping terms
   cf_desc = t_cf_var('scalfac_dd3d', '-',                                      &
     &                'Scaling factor for 3D divergence damping terms', datatype_flt)
   grib2_desc = grib2_var( 0, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
   CALL add_var( p_metrics_list, 'scalfac_dd3d', p_metrics%scalfac_dd3d,           &
        & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc, &
        & ldims = shape1d_c , &
        & lopenacc = .TRUE. )
   __acc_attach(p_metrics%scalfac_dd3d)

    ! Horizontal mask field for 3D divergence damping terms
    ! hmask_dd3d   p_metrics%hmask_dd3d(nproma,nblks_e)
    !
    cf_desc    = t_cf_var('Mask field for 3D divergence damping term', '-',       &
      &                   'Mask field for 3D divergence damping term', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( p_metrics_list, 'hmask_dd3d', p_metrics%hmask_dd3d,      &
                & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
                & ldims=shape2d_e, loutput=.FALSE.,                        &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%hmask_dd3d)

    !----------------------------------------------------------------------------

    ! Explicit weight in vertical wind solver
    ! vwind_expl_wgt   p_metrics%vwind_expl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Explicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Explicit weight in vertical wind solver', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_expl_wgt', p_metrics%vwind_expl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d_c,                                              &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%vwind_expl_wgt)


    ! Implicit weight in vertical wind solver
    ! vwind_impl_wgt  p_metrics%vwind_impl_wgt(nproma,nblks_c)
    !
    cf_desc    = t_cf_var('Implicit_weight_in_vertical_wind_solver', '-',       &
      &                   'Implicit weight in vertical wind solver', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'vwind_impl_wgt', p_metrics%vwind_impl_wgt,   &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,      &
                & ldims=shape2d_c,                                              &
                & lopenacc = .TRUE. )
    __acc_attach(p_metrics%vwind_impl_wgt)


      ! weighting factor for interpolation from full to half levels
      ! wgtfac_c     p_metrics%wgtfac_c(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for interpolation from full to half levels', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfac_c', p_metrics%wgtfac_c,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & ldims=shape3d_chalf, loutput=.FALSE.,                       &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfac_c)


      ! weighting factor for interpolation from full to half levels
      ! wgtfac_e     p_metrics%wgtfac_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for interpolation from full to half levels', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfac_e', p_metrics%wgtfac_e,             &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & ldims=shape3d_ehalf, loutput=.FALSE.,                       &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfac_e)


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_c    p_metrics%wgtfacq_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to surface', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq_c', p_metrics%wgtfacq_c,           &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ccubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfacq_c)


      ! weighting factor for quadratic interpolation to surface
      ! wgtfacq_e    p_metrics%wgtfacq_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to surface', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq_e', p_metrics%wgtfacq_e,           &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ecubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfacq_e)


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_c    p_metrics%wgtfacq1_c(nproma,3,nblks_c)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to model top', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'wgtfacq1_c', p_metrics%wgtfacq1_c,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ccubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfacq1_c)


      ! weighting factor for quadratic interpolation to model top
      ! wgtfacq1_e   p_metrics%wgtfacq1_e(nproma,3,nblks_e)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                            &
      &                     'weighting factor for quadratic interpolation to model top', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'wgtfacq1_e', p_metrics%wgtfacq1_e,         &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape2d_ecubed, loutput=.FALSE.,                      &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfacq1_e)

      ! coefficients for more accurate discretization of grad(E_kin)
      ! coeff_gradekin   p_metrics%coeff_gradekin(nproma,2,nblks_e)
      !
      cf_desc    = t_cf_var('coefficients', '-',                        &
      &                     'coefficients for kinetic energy gradient', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'coeff_gradekin', p_metrics%coeff_gradekin, &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape2d_esquared, loutput=.FALSE.,                    &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%coeff_gradekin)

      ! Inverse layer thickness of full levels
      ! inv_ddqz_z_full   p_metrics%inv_ddqz_z_full(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Inverse_layer_thickness', 'm-1',                   &
      &                     'Inverse layer thickness of full levels', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_full', p_metrics%inv_ddqz_z_full, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,    &
                  & ldims=shape3d_c, loutput=.FALSE.,                             &
                  & isteptype=TSTEP_CONSTANT,                                     &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%inv_ddqz_z_full)


      ! Coefficients for second-order accurate dw/dz term
      ! coeff1_dwdz  p_metrics%coeff1_dwdz(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Coefficient', '',      &
      &                     'Coefficient for second-order accurate dw/dz term', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'coeff1_dwdz', p_metrics%coeff1_dwdz,           &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                  & ldims=shape3d_c, loutput=.FALSE.,                               &
                  & isteptype=TSTEP_CONSTANT,                                       &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%coeff1_dwdz)
      !
      ! coeff2_dwdz  p_metrics%coeff2_dwdz(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Coefficient', '',      &
      &                     'Coefficient for second-order accurate dw/dz term', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'coeff2_dwdz', p_metrics%coeff2_dwdz,           &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                  & ldims=shape3d_c, loutput=.FALSE.,                               &
                  & isteptype=TSTEP_CONSTANT,                                       &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%coeff2_dwdz)

      ! Reference atmosphere field exner
      ! exner_ref_mc  p_metrics%exner_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',                &
      &                     'Reference atmosphere field exner', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_ref_mc', p_metrics%exner_ref_mc,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                  & ldims=shape3d_c,                                                &
                  & isteptype=TSTEP_CONSTANT,                                       &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%exner_ref_mc)


      ! Vertical index of neighbor points needed for Taylor-expansion-based pressure gradient
      ! vertidx_gradp  p_metrics%vertidx_gradp(2,nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Vertical_index', '-',                              &
      &                     'Vertical index', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'vertidx_gradp', p_metrics%vertidx_gradp,   &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_esquared, loutput=.FALSE.,                    &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%vertidx_gradp)

      IF (igradp_method <= 3) THEN
        ! Height differences between local edge point and neighbor cell points used for
        ! pressure gradient computation
        ! zdiff_gradp  p_metrics%zdiff_gradp(2,nproma,nlev,nblks_e)
        !
        cf_desc    = t_cf_var('Height_differences', 'm',                          &
        &                     'Height differences', datatype_flt)
        grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
        CALL add_var( p_metrics_list, 'zdiff_gradp', p_metrics%zdiff_gradp,       &
                    & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                    & ldims=shape3d_esquared, loutput=.FALSE.,                    &
                    & isteptype=TSTEP_CONSTANT,                                   &
                    & lopenacc = .TRUE. )
        __acc_attach(p_metrics%zdiff_gradp)
      ELSE
        ! Coefficients for cubic interpolation of Exner pressure
        ! coeff_gradp  p_metrics%coeff_gradp(8,nproma,nlev,nblks_e)
        !
        cf_desc    = t_cf_var('Interpolation_coefficients', '-',                  &
        &                     'Interpolation coefficients', datatype_flt)
        grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
        CALL add_var( p_metrics_list, 'coeff_gradp', p_metrics%coeff_gradp,       &
                    & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                    & ldims=shape3d_e8, loutput=.FALSE.,                          &
                    & isteptype=TSTEP_CONSTANT,                                   &
                    & lopenacc = .TRUE. )
        __acc_attach(p_metrics%coeff_gradp)
      ENDIF


      ! Extrapolation factor for Exner pressure
      ! exner_exfac  p_metrics%exner_exfac(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Extrapolation_factor_for_Exner_pressure', '-',     &
      &                     'Extrapolation factor for Exner pressure', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'exner_exfac', p_metrics%exner_exfac,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c, loutput=.FALSE.,                           &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%exner_exfac)

      ! Reference atmosphere field theta
      ! theta_ref_mc  p_metrics%theta_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_mc', p_metrics%theta_ref_mc,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%theta_ref_mc)


      ! Reference atmosphere field theta (edges)
      ! theta_ref_me  p_metrics%theta_ref_me(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'theta_ref_me', p_metrics%theta_ref_me,     &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%theta_ref_me)


      ! Reference atmosphere field theta
      ! theta_ref_ic  p_metrics%theta_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_theta', 'K',            &
      &                     'Reference atmosphere field theta', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'theta_ref_ic', p_metrics%theta_ref_ic,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,&
                  & ldims=shape3d_chalf,                                        &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%theta_ref_ic)


      ! Reference surface temperature
      ! tsfc_ref  p_metrics%tsfc_ref(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_surface_temperature', 'K',               &
      &                     'Reference surface temperature', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'tsfc_ref', p_metrics%tsfc_ref,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%tsfc_ref)
      


      ! Reference atmosphere field density
      ! rho_ref_mc  p_metrics%rho_ref_mc(nproma,nlev,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_density', '-',          &
      &                     'Reference atmosphere field density', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'rho_ref_mc', p_metrics%rho_ref_mc,         &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%rho_ref_mc)


      ! Reference atmosphere field density (edges)
      ! rho_ref_me  p_metrics%rho_ref_me(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_density', '-',          &
      &                     'Reference atmosphere field density', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'rho_ref_me', p_metrics%rho_ref_me,         &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,  &
                  & ldims=shape3d_e,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%rho_ref_me)


      ! Reference atmosphere field exner
      ! d_exner_dz_ref_ic  p_metrics%d_exner_dz_ref_ic(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',                  &
      &                     'Reference atmosphere field exner', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'd_exner_dz_ref_ic', p_metrics%d_exner_dz_ref_ic, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                  & ldims=shape3d_chalf,  loutput=.FALSE.,                            &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%d_exner_dz_ref_ic)


      IF (igradp_method <= 3) THEN
        ! Reference atmosphere field exner
        ! d2dexdz2_fac1_mc  p_metrics%d2dexdz2_fac1_mc(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',                &
        &                     'Reference atmosphere field exner', datatype_flt)
        grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_metrics_list, 'd2dexdz2_fac1_mc', p_metrics%d2dexdz2_fac1_mc, &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                    & ldims=shape3d_c,  loutput=.FALSE.,                              &
                    & isteptype=TSTEP_CONSTANT,                                       &
                    & lopenacc = .TRUE. )
        __acc_attach(p_metrics%d2dexdz2_fac1_mc)


        ! Reference atmosphere field exner
        ! d2dexdz2_fac2_mc  p_metrics%d2dexdz2_fac2_mc(nproma,nlev,nblks_c)
        !
        cf_desc    = t_cf_var('Reference_atmosphere_field_exner', '-',                &
        &                     'Reference atmosphere field exner', datatype_flt)
        grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_metrics_list, 'd2dexdz2_fac2_mc', p_metrics%d2dexdz2_fac2_mc, &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,      &
                    & ldims=shape3d_c,  loutput=.FALSE.,                              &
                    & isteptype=TSTEP_CONSTANT,                                       &
                    & lopenacc = .TRUE. )
        __acc_attach(p_metrics%d2dexdz2_fac2_mc)
      ENDIF

      ! index lists for halo points belonging to the nest boundary region
      ! p_metrics%bdy_halo_c_idx
      ! Note: if bdy_halo_c_dim == 0, the array is still allocated with size 1.
      !       ldims has to reflect that, because otherwise it causes serialbox to crash
      !
      cf_desc    = t_cf_var('bdy_halo_c_idx', '-',                                      &
      &                     'index lists for halo points belonging to the nest boundary region', &
      &                     DATATYPE_INT)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'bdy_halo_c_idx', p_metrics%bdy_halo_c_idx, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                  & ldims=(/MAX(p_metrics%bdy_halo_c_dim,1)/), lopenacc = .TRUE. )
      __acc_attach(p_metrics%bdy_halo_c_idx)

      ! block lists for halo points belonging to the nest boundary region
      ! bdy_halo_c_blk
      ! Note: if bdy_halo_c_dim == 0, the array is still allocated with size 1.
      !       ldims has to reflect that, because otherwise it causes serialbox to crash
      !
      cf_desc    = t_cf_var('bdy_halo_c_blk', '-',                                      &
      &                     'block lists for halo points belonging to the nest boundary region', &
      &                     DATATYPE_INT)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'bdy_halo_c_blk', p_metrics%bdy_halo_c_blk, &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,        &
                  & ldims=(/MAX(p_metrics%bdy_halo_c_dim,1)/), lopenacc = .TRUE. )
      __acc_attach(p_metrics%bdy_halo_c_blk)
      
      ! mask field that excludes boundary halo points
      ! mask_prog_halo_c  p_metrics%mask_prog_halo_c(nproma,nblks_c)
      ! Note: Here "loutput" is set to .FALSE. since the output
      !       scheme operates on REAL model variables only and
      !       throws an error on this.
      !
      cf_desc    = t_cf_var('mask_field', '-',                                      &
      &                     'mask field that excludes boundary halo points', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'mask_prog_halo_c', p_metrics%mask_prog_halo_c, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
                  & ldims=shape2d_c, loutput=.FALSE.,                               &
                  & isteptype=TSTEP_CONSTANT,                                       &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%mask_prog_halo_c)

      ! Mask field for mountain or upper slope points
      ! p_metrics%mask_mtnpoints(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('Mask field for mountain points', '-',              &
      &                     'Mask field for mountain points', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'mask_mtnpoints', p_metrics%mask_mtnpoints, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%mask_mtnpoints)

      ! Mask field for mountain or upper slope points (used for gust parameterization)
      ! p_metrics%mask_mtnpoints_g(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('Mask field for mountain points', '-',              &
      &                     'Mask field for mountain points', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'mask_mtnpoints_g', p_metrics%mask_mtnpoints_g, &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%mask_mtnpoints_g)

      ! Slope angle
      ! p_metrics%slope_angle(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('Slope angle', 'rad',                               &
      &                     'Slpe angle', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'slope_angle', p_metrics%slope_angle,       &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%slope_angle)

      ! Slope azimuth
      ! p_metrics%slope_azimuth(nproma,nblks_c)
      !
      cf_desc    = t_cf_var('Slope azimuth', 'rad',                             &
      &                     'Slpe azimuth', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'slope_azimuth', p_metrics%slope_azimuth,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                  & ldims=shape2d_c,                                            &
                  & isteptype=TSTEP_CONSTANT,                                   &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%slope_azimuth)


      IF (nudging_config(jg)%nudge_type==indg_type%ubn .OR. &
        & nudging_config(jg)%nudge_type==indg_type%globn) THEN
        !
        ! vertically varying nudging coefficient
        ! p_metrics%nudgecoeff_vert(nlev)
        !
        cf_desc    = t_cf_var('nudgecoeff_vert', '-',              &
        &                     'vertically varying nudging coefficient', datatype_flt)
        grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( p_metrics_list, 'nudgecoeff_vert', p_metrics%nudgecoeff_vert, &
                    & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,  &
                    & ldims=shape1d_c,                                            &
                    & isteptype=TSTEP_CONSTANT,                                   &
                    & lopenacc = .TRUE. )
        __acc_attach(p_metrics%nudgecoeff_vert)
      ENDIF


    !Add LES related variables : Anurag Dipankar MPIM (2013-04)
#ifndef __NO_ICON_LES__
    IF(atm_phy_nwp_config(jg)%is_les_phy &
      & .OR. aes_vdf_config(jg)%turb==VDIFF_TURB_3DSMAGORINSKY)THEN
#else
    IF(aes_vdf_config(jg)%turb==VDIFF_TURB_3DSMAGORINSKY)THEN
#endif
      ! inv_ddqz_z_half_e  p_metrics%inv_ddqz_z_half_e(nproma,nlevp1,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',                    &
        &                   'metrics functional determinant', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_half_e', p_metrics%inv_ddqz_z_half_e, &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                  & ldims=shape3d_ehalf,                                              &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%inv_ddqz_z_half_e)

      ! inv_ddqz_z_full_e  p_metrics%inv_ddqz_z_full_e(nproma,nlev,nblks_e)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',                    &
        &                   'metrics functional determinant (edge)', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_full_e', p_metrics%inv_ddqz_z_full_e, &
                  & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,        &
                  & ldims=shape3d_e,                                                  &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%inv_ddqz_z_full_e)

      ! inv_ddqz_z_full_v  p_metrics%inv_ddqz_z_full_v(nproma,nlev,nblks_v)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',                    &
        &                   'metrics functional determinant', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_full_v', p_metrics%inv_ddqz_z_full_v, &
                  & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE, cf_desc, grib2_desc,        &
                  & ldims=shape3d_v, loutput=.TRUE.,                                  &
                  & isteptype=TSTEP_CONSTANT ,                                        &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%inv_ddqz_z_full_v)


      ! inv_ddqz_z_half  p_metrics%inv_ddqz_z_half(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',                    &
        &                   'metrics functional determinant', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_half', p_metrics%inv_ddqz_z_half,     &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                  & ldims=shape3d_chalf,                                              &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%inv_ddqz_z_half)


      ! inv_ddqz_z_half_v   p_metrics%inv_ddqz_z_half_v(nproma,nlevp1,nblks_v)
      !
      cf_desc    = t_cf_var('metrics_functional_determinant', '-',                    &
      &                     'metrics functional determinant', datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'inv_ddqz_z_half_v', p_metrics%inv_ddqz_z_half_v, &
                  & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                  & ldims=shape3d_vhalf,                                              &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%inv_ddqz_z_half_v)


      ! mixing_length_sq  p_metrics%mixing_length_sq(nproma,nlevp1,nblks_c)
      !
      cf_desc    = t_cf_var('mixing_length_sq', 'm2','square of mixing length for Smagorinsky model', &
                             datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_metrics_list, 'mixing_length_sq', p_metrics%mixing_length_sq,   &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                  & ldims=shape3d_chalf,                                              &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%mixing_length_sq)


      ! weighting factor for interpolation from full to half levels
      ! wgtfac_v     p_metrics%wgtfac_v(nproma,nlevp1,nblks_v)
      !
      cf_desc    = t_cf_var('weighting_factor', '-',                                  &
      &                     'weighting factor for interpolation from full to half levels(verts)', &
      &                     datatype_flt)
      grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_VERTEX)
      CALL add_var( p_metrics_list, 'wgtfac_v', p_metrics%wgtfac_v,                   &
                  & GRID_UNSTRUCTURED_VERT, ZA_REFERENCE_HALF, cf_desc, grib2_desc,   &
                  & ldims=shape3d_vhalf,                                              &
                  & isteptype=TSTEP_CONSTANT,                                         &
                  & lopenacc = .TRUE. )
      __acc_attach(p_metrics%wgtfac_v)

    END IF !if is_les_phy 

    !----------------------------------------------------------------------------

    !--------------------------------------------------------
    ! metrical modification factors for the deep-atmosphere:
    !--------------------------------------------------------

    ! p_metrics%deepatmo_gradh_mc(nlev):
    ! metrical modification factor for horizontal gradient at full levels
    cf_desc    = t_cf_var('deepatmo_gradh_mc', '-',                                                       &
      &                   'deep-atmosphere metrical modification factor for horizontal gradient', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_gradh_mc', ptr=p_metrics%deepatmo_gradh_mc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE, cf=cf_desc, grib2=grib2_desc,         &
      &           ldims=shape1d_c, initval=1._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_gradh_mc)

    ! p_metrics%deepatmo_divh_mc(nlev):
    ! metrical modification factor for horizontal part of divergence at full levels
    cf_desc    = t_cf_var('deepatmo_divh_mc', '-',                                                      &
      &                   'deep-atmosphere metrical modification factor for horizontal part of divergence', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_divh_mc', ptr=p_metrics%deepatmo_divh_mc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE, cf=cf_desc, grib2=grib2_desc,       &
      &           ldims=shape1d_c, initval=1._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_divh_mc)

    ! p_metrics%deepatmo_vol_mc(nlev):
    ! metrical modification factor for cell volume at full levels
    cf_desc    = t_cf_var('deepatmo_vol_mc', '-',                                                     &
      &                   'deep-atmosphere metrical modification factor for cell volume', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_vol_mc', ptr=p_metrics%deepatmo_vol_mc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE, cf=cf_desc, grib2=grib2_desc,     &
      &           ldims=shape1d_c, initval=1._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_vol_mc)

    ! p_metrics%deepatmo_invr_mc(nlev):
    ! metrical modification factor: inverse of radial distance of full levels from center of Earth
    cf_desc    = t_cf_var('deepatmo_invr_mc', 'm-1',                                                    &
      &                   'deep-atmosphere metrical modification factor -- inverse radial distance from Earth center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_invr_mc', ptr=p_metrics%deepatmo_invr_mc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE, cf=cf_desc, grib2=grib2_desc,       &
      &           ldims=shape1d_c, initval=0._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_invr_mc)

    ! p_metrics%deepatmo_divzU_mc(nlev):
    ! metrical modification factor for vertical part of divergence at full levels
    cf_desc    = t_cf_var('deepatmo_divzU_mc', '-',                                                       &
      &                   'deep-atmosphere metrical modification factor for vertical part of divergence', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_divzU_mc', ptr=p_metrics%deepatmo_divzU_mc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE, cf=cf_desc, grib2=grib2_desc,         &
      &           ldims=shape1d_c, initval=1._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_divzU_mc)

    ! p_metrics%deepatmo_divzL_mc(nlev):
    ! metrical modification factor for vertical part of divergence at full levels
    cf_desc    = t_cf_var('deepatmo_divzL_mc', '-',                                                       &
      &                   'deep-atmosphere metrical modification factor for vertical part of divergence', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_divzL_mc', ptr=p_metrics%deepatmo_divzL_mc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE, cf=cf_desc, grib2=grib2_desc,         &
      &           ldims=shape1d_c, initval=1._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_divzL_mc)

    ! p_metrics%deepatmo_gradh_ifc(nlevp1):
    ! metrical modification factor for horizontal gradient at half levels
    cf_desc    = t_cf_var('deepatmo_gradh_ifc', '-',                                                        &
      &                   'deep-atmosphere metrical modification factor for horizontal gradient', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_gradh_ifc', ptr=p_metrics%deepatmo_gradh_ifc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE_HALF, cf=cf_desc, grib2=grib2_desc,      &
      &           ldims=shape1d_chalf, initval=1._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_gradh_ifc)

    ! p_metrics%deepatmo_invr_ifc(nlevp1):
    ! metrical modification factor: inverse radius from center of Earth of half levels
    cf_desc    = t_cf_var('deepatmo_invr_ifc', 'm-1',                                                     &
      &                   'deep-atmosphere metrical modification factor -- inverse radial distance from Earth center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( this_list=p_metrics_list, varname='deepatmo_invr_ifc', ptr=p_metrics%deepatmo_invr_ifc, &
      &           hgrid=GRID_UNSTRUCTURED_CELL, vgrid=ZA_REFERENCE_HALF, cf=cf_desc, grib2=grib2_desc,    &
      &           ldims=shape1d_chalf, initval=0._wp, loutput=.FALSE., lopenacc = .TRUE.)
    __acc_attach(p_metrics%deepatmo_invr_ifc)

  END SUBROUTINE new_nh_metrics_list

  
  SUBROUTINE new_zd_metrics(p_metrics, p_metrics_list, numpoints)

    TYPE(t_nh_metrics),  INTENT(INOUT):: &  !< diagnostic state
         &  p_metrics

    TYPE(t_var_list_ptr), INTENT(INOUT) :: p_metrics_list   !< diagnostic state list
    
    INTEGER, INTENT(INOUT) :: numpoints

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits
    INTEGER :: datatype_flt

    ibits = DATATYPE_PACK16
    IF ( lnetcdf_flt64_output ) THEN
       datatype_flt = DATATYPE_FLT64
    ELSE
       datatype_flt = DATATYPE_FLT32
    ENDIF

    
    ! Missing description
    ! p_metrics%zd_indlist
    !
    cf_desc    = t_cf_var('zd_indlist', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_indlist', p_metrics%zd_indlist,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/4, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_indlist)

    ! Missing description
    ! p_metrics%zd_blklist
    !
    cf_desc    = t_cf_var('zd_blklist', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_blklist', p_metrics%zd_blklist,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/4, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_blklist)

    ! Missing description
    ! p_metrics%zd_vertidx
    !
    cf_desc    = t_cf_var('zd_vertidx', '-',                                       &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_vertidx', p_metrics%zd_vertidx,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/4, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_vertidx)


    ! Missing description
    ! p_metrics%zd_edgeidx
    !
    cf_desc    = t_cf_var('zd_edgeidx', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_edgeidx', p_metrics%zd_edgeidx,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/3, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_edgeidx)

    ! Missing description
    ! p_metrics%zd_edgeblk
    !
    cf_desc    = t_cf_var('zd_edgeblk', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_edgeblk', p_metrics%zd_edgeblk,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/3, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_edgeblk)

    ! Missing description
    ! p_metrics%zd_geofac
    !
    cf_desc    = t_cf_var('zd_geofac', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_geofac', p_metrics%zd_geofac,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/4, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_geofac)

    ! Missing description
    ! p_metrics%zd_e2cell
    !
    cf_desc    = t_cf_var('zd_e2cell', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_e2cell', p_metrics%zd_e2cell,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/3, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_e2cell)

    ! Missing description
    ! p_metrics%zd_intcoef
    !
    cf_desc    = t_cf_var('zd_intcoef', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_intcoef', p_metrics%zd_intcoef,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/3, numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_intcoef)

    ! Missing description
    ! p_metrics%zd_diffcoef
    !
    cf_desc    = t_cf_var('zd_diffcoef', '-',                                             &
         &                'Missing description' , datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_metrics_list, 'zd_diffcoef', p_metrics%zd_diffcoef,                &
         &           GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,              &
         &           ldims = (/numpoints/), lopenacc = .TRUE. )
    __acc_attach(p_metrics%zd_diffcoef)
    
  END SUBROUTINE new_zd_metrics
  
END MODULE mo_nonhydro_state
