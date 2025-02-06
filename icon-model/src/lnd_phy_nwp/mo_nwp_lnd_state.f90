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

! Contains the data structures to store the physical model state and
! other auxiliary variables in order to run the NWP land physics.
! Constructors and destructors for these data structures as well as
! initialization of fields are also defined here.  This module should
! be an analogon to 'mo_hydro_state.f90'
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

MODULE mo_nwp_lnd_state

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS, MAX_CHAR_LENGTH, HINTP_TYPE_LONLAT_NNB, &
    &                                HINTP_TYPE_LONLAT_BCTR, TLEV_NNOW_RCF,           &
    &                                ALB_SI_MISSVAL, TASK_COMPUTE_SMI
  USE mo_cdi_constants,        ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL
  USE mo_physical_constants,   ONLY: tmelt
  USE mo_master_control,       ONLY: get_my_process_name
  USE mo_parallel_config,      ONLY: nproma
  USE mo_nwp_lnd_types,        ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_exception,            ONLY: message, finish
  USE mo_model_domain,         ONLY: t_patch
  USE mo_grid_config,          ONLY: n_dom
  USE mo_io_config,            ONLY: var_in_output
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config
  USE mo_coupling_config,      ONLY: is_coupled_to_ocean
  USE mo_lnd_nwp_config,       ONLY: nlev_soil, nlev_snow, ntiles_total, &
    &                                lmulti_snow, ntiles_water, lseaice, llake, &
    &                                itype_interception, l2lay_rho_snow, itype_trvg, &
    &                                itype_snowevap, groups_smi, zml_soil
  USE mo_io_config,            ONLY: lnetcdf_flt64_output, runoff_interval, &
    &                                melt_interval
  USE mo_gribout_config,       ONLY: gribout_config
  USE mo_var_list,             ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,    ONLY: vlr_add, vlr_del
  USE mo_var_groups,           ONLY: groups
  USE mo_var_metadata_types,   ONLY: POST_OP_SCALE, CLASS_TILE, CLASS_TILE_LAND
  USE mo_var_metadata,         ONLY: create_hor_interp_metadata, &
    &                                post_op, get_timelevel_string
  USE mo_cf_convention,        ONLY: t_cf_var
  USE mo_grib2,                ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+)
  USE mo_cdi,                  ONLY: DATATYPE_PACK16, DATATYPE_PACK24, DATATYPE_FLT32, &
    &                                TSTEP_ACCUM, GRID_UNSTRUCTURED, DATATYPE_FLT64
  USE mo_zaxis_type,           ONLY: ZA_SURFACE, ZA_SNOW,                        &
    &                                ZA_SNOW_HALF, ZA_DEPTH_BELOW_LAND,          &
    &                                ZA_DEPTH_BELOW_LAND_P1,                     &
    &                                ZA_DEPTH_RUNOFF_S, ZA_DEPTH_RUNOFF_G,       &
    &                                ZA_SEDIMENT_BOTTOM_TW_HALF, ZA_LAKE_BOTTOM, &
    &                                ZA_LAKE_BOTTOM_HALF, ZA_MIX_LAYER
  USE sfc_terra_data,          ONLY: zzhls, zdzhs, zdzms
  USE mo_action_types,         ONLY: ACTION_RESET, actions, new_action


#include "add_var_acc_macro.inc"


  IMPLICIT NONE
  PRIVATE


  !public interface
  !
  ! subroutines
  PUBLIC :: construct_nwp_lnd_state
  PUBLIC :: destruct_nwp_lnd_state

  !
  !variables
  PUBLIC :: p_lnd_state  !> state vector (variable) for land scheme


! complete state vector
!
  TYPE(t_lnd_state), TARGET, ALLOCATABLE :: p_lnd_state(:) 

  CONTAINS


!-------------------------------------------------------------------------
!!            SUBROUTINES FOR BUILDING AND DELETING VARIABLE LISTS 
!-------------------------------------------------------------------------
!
  !>
  !! Constructor for prognostic and diagnostic states.
  !!
  !! It calls constructors to
  !! single time level prognostic states, and diagnostic states.
  !! Initialization of all components with zero.
  !!
  SUBROUTINE construct_nwp_lnd_state(p_patch, p_lnd_state, l_smi, n_timelevels)
!
    TYPE(t_patch), TARGET, INTENT(IN)   :: p_patch(n_dom) !< patch
    LOGICAL,               INTENT(IN)   :: l_smi(n_dom)   !< Flag. TRUE if computation of 
                                                          !< soil moisture index desired
    INTEGER, OPTIONAL, INTENT(IN)       :: n_timelevels   !< number of timelevels

    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(n_dom)
                                           ! nh state at different grid levels
    INTEGER :: ntl, &! local number of timelevels
      &        ist, &! status
      &        jg,  &! grid level counter
      &        jt,  &! time level counter
      &        kso, &! soil level counter
      &      nblks_c ! number of blocks of cells

    CHARACTER(len=MAX_CHAR_LENGTH) :: listname, varname_prefix

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:construct_nwp_lnd_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'land state construction started')

    !$ACC ENTER DATA COPYIN(p_lnd_state)
    ! Allocate and compute fields with information about soil levels and layers
    ALLOCATE (zzhls          (nlev_soil)       , & ! depth of the half level soil layers in m
              zdzhs          (nlev_soil)       , & ! layer thickness between half levels
              zdzms          (nlev_soil)  )        ! distance between main levels

    zzhls(1) = 2.0_wp*zml_soil(1) ! depth of first half level
    zdzhs(1) = zzhls(1)           ! layer thickness betw. half levels of uppermost layer
    zdzms(1) = zml_soil(1)        ! layer thickness between soil surface and main level
                                  !   of uppermost layer

    DO kso = 2,nlev_soil
      zzhls(kso) = zzhls(kso-1) + 2.0_wp*(zml_soil(kso) -zzhls(kso-1))
      zdzhs(kso) = zzhls(kso) - zzhls(kso-1)       ! layer thickness betw. half levels
      zdzms(kso) = zml_soil(kso) - zml_soil(kso-1) ! layer thickness betw. main levels
    ENDDO
    !$ACC ENTER DATA COPYIN(zzhls, zdzhs, zdzms)
    DO jg = 1, n_dom

      IF(PRESENT(n_timelevels))THEN
        ntl = n_timelevels
      ELSE
        ntl = 1
      ENDIF
   
      !
      !determine size of arrays
      nblks_c = p_patch(jg)%nblks_c


      !
      !create state arrays and lists
      !
      ALLOCATE(p_lnd_state(jg)%prog_lnd(1:ntl), &
           &   p_lnd_state(jg)%lnd_prog_nwp_list(1:ntl),STAT=ist)
      !$ACC ENTER DATA COPYIN(p_lnd_state(jg)%prog_lnd)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine), &
             'allocation of land prognostic state array failed')
      ENDIF

      ALLOCATE(p_lnd_state(jg)%prog_wtr(1:ntl), &
               p_lnd_state(jg)%wtr_prog_nwp_list(1:ntl),STAT=ist)
      !$ACC ENTER DATA COPYIN(p_lnd_state(jg)%prog_wtr)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine), &
           'allocation of water prognostic state array failed')
      ENDIF

      !
      !construct prognostic state
      !
      DO jt = 1, ntl

        WRITE(listname,'(a,i2.2,a,i2.2)') 'lnd_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        !varname_prefix = 'lnd_prog_'
        varname_prefix = ''
        CALL  new_nwp_lnd_prog_list(jg, nblks_c, TRIM(listname),             &
          &     TRIM(varname_prefix), p_lnd_state(jg)%lnd_prog_nwp_list(jt), &
          &     p_lnd_state(jg)%prog_lnd(jt), jt)


        WRITE(listname,'(a,i2.2,a,i2.2)') 'wtr_prog_of_domain_',jg, &
          &                               '_and_timelev_',jt

        varname_prefix = ''
        CALL new_nwp_wtr_prog_list(jg, nblks_c, TRIM(listname),             &
          &     TRIM(varname_prefix), p_lnd_state(jg)%wtr_prog_nwp_list(jt), &
          &     p_lnd_state(jg)%prog_wtr(jt), jt)

      ENDDO

      !
      !construct diagnostic state
      !
      WRITE(listname,'(a,i2.2)') 'lnd_diag_of_domain_',jg

      !varname_prefix = 'lnd_diag_'
      varname_prefix = ''
      CALL new_nwp_lnd_diag_list(jg, nblks_c, TRIM(listname),        &
        &   TRIM(varname_prefix), p_lnd_state(jg)%lnd_diag_nwp_list, &
        &   p_lnd_state(jg)%diag_lnd, l_smi(jg))

    ENDDO !ndom

    CALL message (TRIM(routine), 'land state construction completed')

  END SUBROUTINE construct_nwp_lnd_state

!-------------------------------------------------------------------------
!
!
  !>
  !! Destructor for prognostic and diagnostic states.
  !!
  !! It calls destructors to
  !! single time level prognostic states, and diagnostic states.
  !!
  SUBROUTINE destruct_nwp_lnd_state()
    !
    INTEGER :: ntl, &! local number of timelevels
      &        ist, &! status
      &         jg, &! grid level counter
      &         jt   ! time level counter

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = 'mo_nwp_lnd_state:destruct_nwp_lnd_state'
!-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'Destruction of nwp land state started')

    ! Deallocate fields with information about soil levels and layers
    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(zzhls, zdzhs, zdzms)
    DEALLOCATE (zzhls, zdzhs, zdzms)
    
    DO jg = 1, n_dom

      ntl = SIZE(p_lnd_state(jg)%prog_lnd(:))
      IF(ntl==0)THEN
        CALL finish(TRIM(routine), 'prognostic array has no timelevels')
      ENDIF


      DO jt = 1, ntl
        ! delete prognostic state list elements
        CALL vlr_del(p_lnd_state(jg)%lnd_prog_nwp_list(jt))
      ENDDO

      DO jt = 1, ntl
        ! delete prognostic state list elements
        CALL vlr_del(p_lnd_state(jg)%wtr_prog_nwp_list(jt))
      ENDDO

      ! delete diagnostic state list elements
      CALL vlr_del(p_lnd_state(jg)%lnd_diag_nwp_list)

      !$ACC EXIT DATA DELETE(p_lnd_state(jg)%prog_lnd)
      ! destruct state lists and arrays
      DEALLOCATE(p_lnd_state(jg)%prog_lnd,&
           &     p_lnd_state(jg)%lnd_prog_nwp_list, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine),  &
             'deallocation of land prognostic state array failed')
      ENDIF

      IF (lseaice .OR. llake) THEN
        ! destruct state lists and arrays
        !$ACC EXIT DATA DELETE(p_lnd_state(jg)%prog_wtr)
        DEALLOCATE(p_lnd_state(jg)%prog_wtr,&
             &     p_lnd_state(jg)%wtr_prog_nwp_list, STAT=ist)
        IF(ist/=SUCCESS)THEN
          CALL finish (TRIM(routine),  &
               'deallocation of water prognostic state array failed')
        ENDIF
      ENDIF

    ENDDO

    !$ACC EXIT DATA DELETE(p_lnd_state)
    DEALLOCATE(p_lnd_state, STAT=ist)
    IF(ist/=success)THEN
      CALL finish (TRIM(routine), 'deallocation of p_lnd_state failed')
    ENDIF

    CALL message (TRIM(routine), 'Destruction of nwp land state completed')

  END SUBROUTINE destruct_nwp_lnd_state

!-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of prognostic land state.
  !!
  !! Initialization of components with zero.
  !!
  SUBROUTINE new_nwp_lnd_prog_list( p_jg, kblks, listname, vname_prefix, &
    &                               prog_list, p_prog_lnd, timelev )

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg  !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(LEN=2)            :: csfc

    TYPE(t_var_list_ptr),INTENT(INOUT) :: prog_list
    TYPE(t_lnd_prog),INTENT(INOUT) :: p_prog_lnd

    INTEGER, INTENT(IN) :: timelev

    ! Local variables
    !
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape3d_subsw(3)
    INTEGER :: shape4d_snow_subs(4)
    INTEGER :: ibits
    INTEGER :: jsfc               !< tile counter
    INTEGER :: DATATYPE_PACK_VAR  !< variable packing accuracy for t_so
    INTEGER :: datatype_flt

    CHARACTER(len=4) suffix

!-----------------------------------------------------------------------

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    ibits = DATATYPE_PACK16 ! packing accuracy of horizontal slice

    IF (gribout_config(p_jg)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for t_so
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for all variables
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    ! predefined array shapes
    shape2d              = (/nproma,            kblks            /)
    shape3d_subs         = (/nproma,            kblks, ntiles_total /)
    shape3d_subsw        = (/nproma,            kblks, ntiles_total+ntiles_water /)
    shape4d_snow_subs    = (/nproma, nlev_snow, kblks, ntiles_total /)


    ! Suffix (mandatory for time level dependent variables)

    suffix = get_timelevel_string(timelev)


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_prog_lnd%t_s_t, &
    &       p_prog_lnd%t_sk_t, &
    &       p_prog_lnd%t_g, &
    &       p_prog_lnd%t_g_t, &
    &       p_prog_lnd%w_i_t, &
    &       p_prog_lnd%w_p_t, &
    &       p_prog_lnd%w_s_t, &
    &       p_prog_lnd%t_so_t, &
    &       p_prog_lnd%w_so_t, &
    &       p_prog_lnd%w_so_ice_t, &
    &       p_prog_lnd%t_snow_t, &
    &       p_prog_lnd%w_snow_t, &
    &       p_prog_lnd%rho_snow_t, &
    &       p_prog_lnd%t_snow_mult_t, &
    &       p_prog_lnd%wtot_snow_t, &
    &       p_prog_lnd%wliq_snow_t, &
    &       p_prog_lnd%rho_snow_mult_t, &
    &       p_prog_lnd%dzh_snow_t)


    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(prog_list, TRIM(listname), patch_id=p_jg, &
      &          lrestart=.TRUE., model_type=get_my_process_name())

    !------------------------------

    ! & p_prog_lnd%t_g(nproma,nblks_c), STAT = ist)
    cf_desc    = t_cf_var('t_g', 'K', 'weighted surface temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_g'//suffix, p_prog_lnd%t_g,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape2d,                                                      &
         & tlev_source=TLEV_NNOW_RCF,                      &! for output take field from nnow_rcf slice
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",     &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in",              &
         &                 "mode_combined_in","mode_cosmo_in","mode_iniana"),  &
         & lopenacc=.TRUE. ) 
    __acc_attach(p_prog_lnd%t_g)


    ! & p_prog_lnd%t_g_t(nproma,nblks_c,ntiles_total+ntiles_water), STAT = ist)
    cf_desc    = t_cf_var('t_g_t', 'K', 'weighted surface temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_g_t'//suffix, p_prog_lnd%t_g_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,           &
         & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE.,           &
         & loutput=.FALSE.,                                                    &
         & tlev_source=TLEV_NNOW_RCF,                                          &
         & in_group=groups("iau_restore_vars"),                                &
         & lopenacc=.TRUE. )  
    __acc_attach(p_prog_lnd%t_g_t)


    ! fill the separate variables belonging to the container t_gt
    ALLOCATE(p_prog_lnd%t_gt_ptr(ntiles_total+ntiles_water))
      DO jsfc = 1,ntiles_total+ntiles_water
        NULLIFY(p_prog_lnd%t_gt_ptr(jsfc)%p_2d, p_prog_lnd%t_gt_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'t_g_t'//suffix,                &
               & vname_prefix//'t_g_t_'//TRIM(ADJUSTL(csfc))//suffix,          &
               & p_prog_lnd%t_gt_ptr(jsfc)%p_2d,                               &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('t_g_t_'//TRIM(csfc), '', '', datatype_flt),       &
               & grib2_var(0, 0, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
               & ref_idx=jsfc,                                                 &
               & ldims=shape2d,                                                &
               & var_class=CLASS_TILE,                                         &
               & tlev_source=TLEV_NNOW_RCF,                                    & ! for output take field from nnow_rcf slice
               & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
      ENDDO

    IF ( atm_phy_nwp_config(p_jg)%inwp_surface > 0 ) THEN

    ! & p_prog_lnd%t_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('t_s_t', 'K', 'temperature of ground surface', datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_s_t'//suffix, p_prog_lnd%t_s_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,&
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%t_s_t)

    ! fill the separate variables belonging to the container t_s
    ALLOCATE(p_prog_lnd%t_s_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      NULLIFY(p_prog_lnd%t_s_ptr(jsfc)%p_2d, p_prog_lnd%t_s_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_s_t'//suffix,              &
           & vname_prefix//'t_s_t_'//TRIM(ADJUSTL(csfc))//suffix,          &
           & p_prog_lnd%t_s_ptr(jsfc)%p_2d,                                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('t_s_t_'//csfc, '', '', datatype_flt),               &
           & grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
           & ref_idx=jsfc,                                                 &
           & ldims=shape2d,                                                &
           & var_class=CLASS_TILE,                                         &
           & tlev_source=TLEV_NNOW_RCF, in_group=groups("land_tile_vars") ) ! for output take field from nnow_rcf slice
    ENDDO


    ! & p_prog_lnd%t_sk_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('t_sk_t', 'K', 'skin temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_sk_t'//suffix, p_prog_lnd%t_sk_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
         & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE.,       &
         & loutput=.FALSE.,                                                &
         & tlev_source=TLEV_NNOW_RCF,                                      &
         & in_group=groups("iau_restore_vars"),                            &
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%t_sk_t)

    ! fill the separate variables belonging to the container t_sk
    ALLOCATE(p_prog_lnd%t_sk_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      NULLIFY(p_prog_lnd%t_sk_ptr(jsfc)%p_2d, p_prog_lnd%t_sk_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_sk_t'//suffix,             &
           & vname_prefix//'t_sk_t_'//TRIM(ADJUSTL(csfc))//suffix,         &
           & p_prog_lnd%t_sk_ptr(jsfc)%p_2d,                               &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('t_sk_t_'//csfc, '', '', datatype_flt),              &
           & grib2_var(0, 0, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
           & ref_idx=jsfc,                                                 &
           & ldims=shape2d,                                                &
           & var_class=CLASS_TILE, tlev_source=TLEV_NNOW_RCF,              & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t")         )
    ENDDO




    ! & p_prog_lnd%w_i_t(nproma,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_i_t', 'm H2O', 'water content of interception water', datatype_flt)
    new_cf_desc = t_cf_var('w_i_t', 'kg/m**2', 'weighted water content of interception water', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_i_t'//suffix, p_prog_lnd%w_i_t,     &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,               &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,               &
         & loutput=.FALSE.,                                                       &
         & tlev_source=TLEV_NNOW_RCF,                                             &
         & in_group=groups("iau_restore_vars"),                                   &
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%w_i_t)

    ! fill the separate variables belonging to the container w_i
    ALLOCATE(p_prog_lnd%w_i_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_i_ptr(jsfc)%p_2d, p_prog_lnd%w_i_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_i_t'//suffix,                &
           & vname_prefix//'w_i_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
           & p_prog_lnd%w_i_ptr(jsfc)%p_2d,                                  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
           & t_cf_var('w_i_t_'//csfc, '', '', datatype_flt),               &
           & grib2_var(2, 0, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ref_idx=jsfc,                                                   &
           & ldims=shape2d,                                                  &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),          &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO


    IF (itype_interception == 2) THEN

      ! & p_prog_lnd%w_p_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('w_p_t', 'm H2O', 'water content of interception water', datatype_flt)
      grib2_desc = grib2_var(2, 0, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'w_p_t'//suffix, p_prog_lnd%w_p_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,&
           & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%w_p_t)

      ! fill the separate variables belonging to the container w_p
      ALLOCATE(p_prog_lnd%w_p_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%w_p_ptr(jsfc)%p_2d, p_prog_lnd%w_p_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'w_p_t'//suffix,                &
             & vname_prefix//'w_p_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
             & p_prog_lnd%w_p_ptr(jsfc)%p_2d,                                  &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
             & t_cf_var('w_p_t_'//csfc, '', '', datatype_flt),               &
             & grib2_var(2, 0, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
             & ref_idx=jsfc,                                                   &
             & ldims=shape2d,                                                  &
             & var_class=CLASS_TILE_LAND,                                      &
             & tlev_source=TLEV_NNOW_RCF,                                      &
             & in_group=groups("land_tile_vars")) ! for output take field from nnow_rcf slice
      ENDDO

      ! & p_prog_lnd%w_s_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('w_s_t', 'm H2O', 'water content of interception water', datatype_flt)
      grib2_desc = grib2_var(2, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'w_s_t'//suffix, p_prog_lnd%w_s_t,     &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,grib2_desc,                 &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,&
           & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%w_s_t)

      ! fill the separate variables belonging to the container w_s
      ALLOCATE(p_prog_lnd%w_s_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%w_s_ptr(jsfc)%p_2d, p_prog_lnd%w_s_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'w_s_t'//suffix,                &
             & vname_prefix//'w_s_t_'//TRIM(ADJUSTL(csfc))//suffix,            &
             & p_prog_lnd%w_s_ptr(jsfc)%p_2d,                                  &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
             & t_cf_var('w_s_t_'//csfc, '', '', datatype_flt),               &
             & grib2_var(2, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
             & ref_idx=jsfc,                                                   &
             & ldims=shape2d,                                                  &
             & var_class=CLASS_TILE_LAND,                                      &
             & tlev_source=TLEV_NNOW_RCF,                                      &
             & in_group=groups("land_tile_vars")) ! for output take field from nnow_rcf slice
      ENDDO
    END IF  ! itype_interception == 2


    ! & p_prog_lnd%t_so_t(nproma,nlev_soil+1,nblks_c,ntiles_total) 
    cf_desc    = t_cf_var('t_so_t', 'K', 'soil temperature (main level)', datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_so_t'//suffix, p_prog_lnd%t_so_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND_P1, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+1,kblks,ntiles_total/),                      &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                 &
         & tlev_source=TLEV_NNOW_RCF,                                            &
         & in_group=groups("iau_restore_vars"),                                  &
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%t_so_t)

    ! fill the separate variables belonging to the container t_so
    !
    ALLOCATE(p_prog_lnd%t_so_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%t_so_ptr(jsfc)%p_2d, p_prog_lnd%t_so_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_so_t'//suffix,               &
           & vname_prefix//'t_so_t_'//TRIM(ADJUSTL(csfc))//suffix,           &
           & p_prog_lnd%t_so_ptr(jsfc)%p_3d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND_P1,                 &
           & t_cf_var('t_so_t_'//csfc, '', '', datatype_flt),                &
           & grib2_var(2, 3, 18, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL),  &
           & ref_idx=jsfc,                                                   &
           & ldims=(/nproma,nlev_soil+1,kblks/),                             &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    ENDDO



   ! & p_prog_lnd%w_so_t(nproma,nlev_soil,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_so_t', 'm H20', 'total water content (ice + liquid water)', &
         &                datatype_flt)

    new_cf_desc = t_cf_var('w_so_t', 'kg/m**2', 'total water content (ice + liquid water)', datatype_flt)
    grib2_desc  = grib2_var(2, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_so_t'//suffix, p_prog_lnd%w_so_t,  &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc,     &
         & ldims=(/nproma,nlev_soil,kblks,ntiles_total/),                        &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                 &
         & tlev_source=TLEV_NNOW_RCF,                                            &
         & in_group=groups("iau_restore_vars"),                                  &
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%w_so_t)

    ! fill the separate variables belonging to the container w_so
    ALLOCATE(p_prog_lnd%w_so_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_so_ptr(jsfc)%p_2d, p_prog_lnd%w_so_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_so_t'//suffix,               &
           & vname_prefix//'w_so_t_'//TRIM(ADJUSTL(csfc))//suffix,           &
           & p_prog_lnd%w_so_ptr(jsfc)%p_3d,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                    &
           & t_cf_var('w_so_t_'//csfc, '', '', datatype_flt),              &
           & grib2_var(2, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ref_idx=jsfc,                                                   &
           & ldims=(/nproma,nlev_soil,kblks/),                               &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), & 
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),          &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO



    ! & p_prog_lnd%w_so_ice_t(nproma,nlev_soil,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_so_ice_t', 'm H20', 'ice content', datatype_flt)
    new_cf_desc = t_cf_var('w_so_ice_t', 'kg/m**2', 'ice content', datatype_flt)
    grib2_desc = grib2_var(2, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_so_ice_t'//suffix,                 &
         & p_prog_lnd%w_so_ice_t, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,   &
         & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil,kblks,ntiles_total/),   &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                 &
         & tlev_source=TLEV_NNOW_RCF,                                            &
         & in_group=groups("iau_restore_vars"),                                  &
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%w_so_ice_t)

    ! fill the separate variables belonging to the container w_so_ice
    ALLOCATE(p_prog_lnd%w_so_ice_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_so_ice_ptr(jsfc)%p_2d, p_prog_lnd%w_so_ice_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_so_ice_t'//suffix,           &
           & vname_prefix//'w_so_ice_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%w_so_ice_ptr(jsfc)%p_3d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                    &
           & t_cf_var('w_so_ice_t_'//csfc, '', '', datatype_flt),          &
           & grib2_var(2, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ref_idx=jsfc,                                                   &
           & ldims=(/nproma,nlev_soil,kblks/),                               &
           & var_class=CLASS_TILE_LAND,                                      &
           & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),          &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO




    ! & p_prog_lnd%t_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('t_snow_t', 'K', 'temperature of the snow-surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_snow_t'//suffix, p_prog_lnd%t_snow_t,&
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,               &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,                &
         & loutput=.FALSE.,                                                        &
         & tlev_source=TLEV_NNOW_RCF,                                              &
         & in_group=groups("iau_restore_vars"),                                    &
         & lopenacc=.TRUE. ) 
    __acc_attach(p_prog_lnd%t_snow_t)

    ! fill the separate variables belonging to the container t_snow
    ALLOCATE(p_prog_lnd%t_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%t_snow_ptr(jsfc)%p_2d, p_prog_lnd%t_snow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'t_snow_t'//suffix,             &
             & vname_prefix//'t_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
             & p_prog_lnd%t_snow_ptr(jsfc)%p_2d,                             &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
             & t_cf_var('t_snow_t_'//csfc, '', '', datatype_flt),          &
             & grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
             & ref_idx=jsfc,                                                 &
             & ldims=shape2d,                                                &
             & var_class=CLASS_TILE_LAND,                                    &
             & tlev_source=TLEV_NNOW_RCF,                                    & ! for output take field from nnow_rcf slice
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    ENDDO


    ! & p_prog_lnd%w_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc     = t_cf_var('w_snow_t', 'm H2O', 'water equivalent of snow', datatype_flt)
    new_cf_desc = t_cf_var('w_snow_t', 'kg m-2', 'water equivalent of snow', datatype_flt)
    grib2_desc = grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'w_snow_t'//suffix, p_prog_lnd%w_snow_t,&
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,               &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,                &
         & loutput=.FALSE.,                                                        &
         & tlev_source=TLEV_NNOW_RCF,                                              &
         & in_group=groups("iau_restore_vars"),                                    &
         & lopenacc=.TRUE.)
    __acc_attach(p_prog_lnd%w_snow_t)

    ! fill the separate variables belonging to the container w_snow
    ALLOCATE(p_prog_lnd%w_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%w_snow_ptr(jsfc)%p_2d, p_prog_lnd%w_snow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc  
      CALL add_ref( prog_list, vname_prefix//'w_snow_t'//suffix,           &
           & vname_prefix//'w_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%w_snow_ptr(jsfc)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
           & t_cf_var('w_snow_t_'//csfc, '', '', datatype_flt),          &
           & grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
           & ref_idx=jsfc,                                                 &
           & ldims=shape2d,                                                &
           & var_class=CLASS_TILE_LAND,                                    &
           & tlev_source=TLEV_NNOW_RCF,                                    & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t"),        &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc) )
    ENDDO


    ! & p_prog_lnd%rho_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('rho_snow_t', 'kg/m**3', 'snow density', datatype_flt)
    grib2_desc = grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'rho_snow_t'//suffix, p_prog_lnd%rho_snow_t, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,                    &
         & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,                     &
         & loutput=.FALSE.,                                                             &
         & tlev_source=TLEV_NNOW_RCF,                                                   &
         & in_group=groups("iau_restore_vars"),                                         &
         & lopenacc=.TRUE.)
        __acc_attach(p_prog_lnd%rho_snow_t)

    ! fill the separate variables belonging to the container rho_snow
    ALLOCATE(p_prog_lnd%rho_snow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_prog_lnd%rho_snow_ptr(jsfc)%p_2d, p_prog_lnd%rho_snow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( prog_list, vname_prefix//'rho_snow_t'//suffix,           &
           & vname_prefix//'rho_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
           & p_prog_lnd%rho_snow_ptr(jsfc)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
           & t_cf_var('rho_snow_t_'//csfc, '', '', datatype_flt),          &
           & grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
           & ref_idx=jsfc,                                                   &
           & ldims=shape2d,                                                  &
           & var_class=CLASS_TILE_LAND,                                      &
           & tlev_source=TLEV_NNOW_RCF,                                      & ! for output take field from nnow_rcf slice
           & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    END DO


    IF (lmulti_snow .OR. l2lay_rho_snow) THEN

      ! & p_prog_lnd%rho_snow_mult_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('rho_snow_mult_t', 'kg/m**3', 'snow density', datatype_flt)
      grib2_desc = grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'rho_snow_mult_t'//suffix,              &
           & p_prog_lnd%rho_snow_mult_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,            &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                           &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                   &
           & tlev_source=TLEV_NNOW_RCF,                                              &
           & in_group=groups("iau_restore_vars"),                                    &
           & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%rho_snow_mult_t)

      ! fill the separate variables belonging to the container rho_snow_mult
      !
      ALLOCATE(p_prog_lnd%rho_snow_mult_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_2d, p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'rho_snow_mult_t'//suffix,       &
             & vname_prefix//'rho_snow_mult_t_'//TRIM(ADJUSTL(csfc))//suffix,   &
             & p_prog_lnd%rho_snow_mult_ptr(jsfc)%p_3d,                         &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                 &
             & t_cf_var('rho_snow_mult_t_'//csfc, '', '', datatype_flt),      &
             & grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
             & ref_idx=jsfc,                                                    &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,               &
             & var_class=CLASS_TILE_LAND,                                       &
             & tlev_source=TLEV_NNOW_RCF,                                       &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice 
      ENDDO

    ENDIF

    IF (lmulti_snow) THEN

      ! & p_prog_lnd%t_snow_mult_t(nproma,nlev_snow+1,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('t_snow_mult_t', 'K', 'temperature of the snow-surface', &
           & datatype_flt)
      grib2_desc = grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_snow_mult_t'//suffix, p_prog_lnd%t_snow_mult_t, &
       & GRID_UNSTRUCTURED_CELL, ZA_SNOW_HALF, cf_desc, grib2_desc,               &
       & ldims=(/nproma,nlev_snow+1,kblks,ntiles_total/),                         &
       & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                    &
       & tlev_source=TLEV_NNOW_RCF,                                               &
       & in_group=groups("iau_restore_vars"),                                     &
       & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%t_snow_mult_t)

      ! fill the separate variables belonging to the container t_snow_mult
      !
      ALLOCATE(p_prog_lnd%t_snow_mult_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%t_snow_mult_ptr(jsfc)%p_2d, p_prog_lnd%t_snow_mult_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'t_snow_mult_t'//suffix,      &
             & vname_prefix//'t_snow_mult_t_'//TRIM(ADJUSTL(csfc))//suffix,  &
             & p_prog_lnd%t_snow_mult_ptr(jsfc)%p_3d,                        &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW_HALF,                         &
             & t_cf_var('t_snow_mult_t_'//csfc, '', '', datatype_flt),     &
             & grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
             & ref_idx=jsfc,                                                 &
             & ldims=(/nproma,nlev_snow+1,kblks/),                           &
             & var_class=CLASS_TILE_LAND,                                    &
             & lrestart=.TRUE.,                                              &
             & tlev_source=TLEV_NNOW_RCF,                                    &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice 
      ENDDO



      ! & p_prog_lnd%wtot_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('wtot_snow_t', 'm H2O', 'total water content in snow', datatype_flt)
      grib2_desc = grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'wtot_snow_t'//suffix,                &
           & p_prog_lnd%wtot_snow_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,              &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                 &
           & tlev_source=TLEV_NNOW_RCF,                                            &
           & in_group=groups("iau_restore_vars"),                                  &
           & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%wtot_snow_t)

      ! fill the separate variables belonging to the container wtot_snow
      !
      ALLOCATE(p_prog_lnd%wtot_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%wtot_snow_ptr(jsfc)%p_2d, p_prog_lnd%wtot_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'wtot_snow_t'//suffix,          &
             & vname_prefix//'wtot_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,      &
             & p_prog_lnd%wtot_snow_ptr(jsfc)%p_3d,                            &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                &
             & t_cf_var('wtot_snow_t_'//csfc, '', '', datatype_flt),           &
             & grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
             & ref_idx=jsfc,                                                   &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
             & var_class=CLASS_TILE_LAND, tlev_source=TLEV_NNOW_RCF,           &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice
      ENDDO



      ! & p_prog_lnd%wliq_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('wliq_snow_t', 'm H2O', 'liquid water content in snow', datatype_flt)
      grib2_desc = grib2_var(0, 1, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'wliq_snow_t'//suffix,                &
           & p_prog_lnd%wliq_snow_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,              &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                 &
           & tlev_source=TLEV_NNOW_RCF,                                            &
           & in_group=groups("iau_restore_vars"),                                  &
           & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%wliq_snow_t)

      ! fill the separate variables belonging to the container wliq_snow
      !
      ALLOCATE(p_prog_lnd%wliq_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%wliq_snow_ptr(jsfc)%p_2d, p_prog_lnd%wliq_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( prog_list, vname_prefix//'wliq_snow_t'//suffix,          &
             & vname_prefix//'wliq_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,      &
             & p_prog_lnd%wliq_snow_ptr(jsfc)%p_3d,                            &
             & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                &
             & t_cf_var('wliq_snow_t_'//csfc, '', '', datatype_flt),           &
             & grib2_var(0, 1, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
             & ref_idx=jsfc,                                                   &
             & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
             & var_class=CLASS_TILE_LAND, tlev_source=TLEV_NNOW_RCF,           &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice
      ENDDO


      ! & p_prog_lnd%dzh_snow_t(nproma,nlev_snow,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('dzh_snow_t', 'm', 'layer thickness between half levels in snow', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'dzh_snow_t'//suffix,                 &
           & p_prog_lnd%dzh_snow_t, GRID_UNSTRUCTURED_CELL, ZA_SNOW,               &
           & cf_desc, grib2_desc, ldims=shape4d_snow_subs,                         &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,                 &
           & tlev_source=TLEV_NNOW_RCF,                                            &
           & in_group=groups("iau_restore_vars"),                                  &
           & lopenacc=.TRUE.)
      __acc_attach(p_prog_lnd%dzh_snow_t)

      ! fill the separate variables belonging to the container dzh_snow
      !
      ALLOCATE(p_prog_lnd%dzh_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_prog_lnd%dzh_snow_ptr(jsfc)%p_2d, p_prog_lnd%dzh_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc  
        CALL add_ref( prog_list, vname_prefix//'dzh_snow_t'//suffix,             &
               & vname_prefix//'dzh_snow_t_'//TRIM(ADJUSTL(csfc))//suffix,       &
               & p_prog_lnd%dzh_snow_ptr(jsfc)%p_3d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_SNOW,                                &
               & t_cf_var('dzh_snow_t_'//csfc, '', '', datatype_flt),            &
               & grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL),       &
               & ref_idx=jsfc,                                                   &
               & ldims=(/nproma,nlev_snow,kblks/), lrestart=.TRUE.,              &
               & var_class=CLASS_TILE_LAND, tlev_source=TLEV_NNOW_RCF,           &
               & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") ) ! for output take field from nnow_rcf slice
      ENDDO

    ENDIF  ! lmulti_snow

    END IF !inwp_surface > 0

  END SUBROUTINE new_nwp_lnd_prog_list



  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of prognostic lake/seaice state.
  !!
  !! Initialization of components with zero.
  !!
  !! Modifications by Dmitrii Mironov, DWD (2016-08-04)
  !! - Prognostic sea-ice albedo is added.
  !!
  SUBROUTINE new_nwp_wtr_prog_list( p_jg, kblks, listname, vname_prefix, &
    &                               prog_list, p_prog_wtr, timelev )

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg  !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix

    TYPE(t_var_list_ptr),INTENT(INOUT) :: prog_list
    TYPE(t_wtr_prog),INTENT(INOUT) :: p_prog_wtr

    INTEGER, INTENT(IN) :: timelev

    ! Local variables
    !
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2)
    INTEGER :: ibits

    CHARACTER(len=4) :: suffix
    INTEGER          :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

!-----------------------------------------------------------------------

    ibits = DATATYPE_PACK16 ! "entropy" of horizontal slice

    ! predefined array shapes
    shape2d              = (/nproma, kblks/)

    ! Suffix (mandatory for time level dependent variables)

    suffix = get_timelevel_string(timelev)


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_prog_wtr%t_ice, &
    &       p_prog_wtr%h_ice, &
    &       p_prog_wtr%t_snow_si, &
    &       p_prog_wtr%h_snow_si, &
    &       p_prog_wtr%alb_si,    &
    &       p_prog_wtr%t_snow_lk,  &
    &       p_prog_wtr%h_snow_lk,  &
    &       p_prog_wtr%t_mnw_lk,   &
    &       p_prog_wtr%t_wml_lk,   &
    &       p_prog_wtr%h_ml_lk,    &
    &       p_prog_wtr%t_bot_lk,   &
    &       p_prog_wtr%c_t_lk,     &
    &       p_prog_wtr%t_b1_lk,    &
    &       p_prog_wtr%h_b1_lk)


    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(prog_list, TRIM(listname), patch_id=p_jg, &
      &          lrestart=.TRUE., model_type=get_my_process_name())

    !------------------------------

    ! h_ice-field is also needed by turbdiff irrespective of whether 
    ! the seaice model is used or not.

    !
    ! sea-ice and lake-model
    !
    ! since it is currently not envisaged to have mixed sea-lake gridpoints, t_ice 
    ! is used for both sea- and lake-ice temperatures. This is in accordance with 
    ! the COSMO implementation 
    ! & p_prog_wtr%t_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('t_ice', 'K', 'sea/lake-ice temperature', datatype_flt)
    grib2_desc = grib2_var(10, 2, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_ice'//suffix, p_prog_wtr%t_ice,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                           &
         & initval = tmelt,                                                    &
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in", &
         &                 "mode_iau_old_fg_in","mode_combined_in","mode_cosmo_in", &
         &                 "mode_iniana","iau_restore_vars"), lopenacc=.TRUE. )
    __acc_attach(p_prog_wtr%t_ice)


    ! since it is currently not envisaged to have mixed sea-lake gridpoints, h_ice 
    ! is used for both sea- and lake-ice thickness. This is in accordance with 
    ! the COSMO implementation 
    ! & p_prog_wtr%h_ice(nproma,nblks_c)
    cf_desc    = t_cf_var('h_ice', 'm', 'sea/lake-ice depth', datatype_flt)
    grib2_desc = grib2_var(10, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'h_ice'//suffix, p_prog_wtr%h_ice,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                           &
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in", &
         &                 "mode_iau_old_fg_in","mode_combined_in",            &
         &                 "mode_cosmo_in","mode_iniana","iau_restore_vars"),  &
         & lopenacc=.TRUE. )   
    __acc_attach(p_prog_wtr%h_ice)


    !
    ! sea-ice model specific
    !

    ! & p_prog_wtr%t_snow_si(nproma,nblks_c)
    cf_desc    = t_cf_var('t_snow_si', 'K', 'temperature of snow on sea ice', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'t_snow_si'//suffix, p_prog_wtr%t_snow_si,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
         & tlev_source=TLEV_NNOW_RCF, lopenacc=.TRUE.)
    __acc_attach(p_prog_wtr%t_snow_si)


    ! & p_prog_wtr%h_snow_si(nproma,nblks_c)
    cf_desc    = t_cf_var('h_snow_si', 'm', 'depth of snow on sea ice', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'h_snow_si'//suffix, p_prog_wtr%h_snow_si,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
         & tlev_source=TLEV_NNOW_RCF, lopenacc=.TRUE.)
    __acc_attach(p_prog_wtr%h_snow_si)

    ! & p_prog_wtr%alb_si(nproma,nblks_c)
    cf_desc     = t_cf_var('sea_ice_albedo', '-', 'sea ice albedo (diffuse)', datatype_flt)
    new_cf_desc = t_cf_var('sea_ice_albedo', '%', 'sea ice albedo (diffuse)', datatype_flt)
    grib2_desc = grib2_var(0, 19, 234, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prog_list, vname_prefix//'alb_si'//suffix, p_prog_wtr%alb_si,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                             &
         & initval = ALB_SI_MISSVAL,                                             &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), & 
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in", "mode_iau_fg_in", &
         &  "mode_iau_old_fg_in","mode_cosmo_in","mode_iniana",                  &
         &  "mode_combined_in","iau_restore_vars"),                              &
         & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),     &
         & lopenacc=.TRUE. )   
    __acc_attach(p_prog_wtr%alb_si)


    !
    !FLAKE
    !

    IF (llake) THEN   ! lake model switched on

      ! p_prog_wtr%t_snow_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_snow_lk', 'K', 'temperature of snow on lake ice', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_snow_lk'//suffix, p_prog_wtr%t_snow_lk,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
           & tlev_source=TLEV_NNOW_RCF, lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%t_snow_lk)


      ! p_prog_wtr%h_snow_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('h_snow_lk', 'm', 'depth of snow on lake ice', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'h_snow_lk'//suffix, p_prog_wtr%h_snow_lk,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, ldims=shape2d,     &
           & tlev_source=TLEV_NNOW_RCF, lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%h_snow_lk)


      ! p_prog_wtr%t_mnw_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_mnw_lk', 'K', 'mean temperature of the water column', datatype_flt)
      grib2_desc = grib2_var(1, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( prog_list, vname_prefix//'t_mnw_lk'//suffix, p_prog_wtr%t_mnw_lk,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF,   &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in","mode_iniana",  &
           & "iau_restore_vars"),                                                                                     &
           & lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%t_mnw_lk)


      ! p_prog_wtr%t_wml_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_wml_lk', 'K', 'mixed-layer temperature', datatype_flt)
      grib2_desc = grib2_var(1, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( prog_list, vname_prefix//'t_wml_lk'//suffix, p_prog_wtr%t_wml_lk,                                 &
           & GRID_UNSTRUCTURED_CELL, ZA_MIX_LAYER, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF,     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in","mode_iniana",  &
           & "iau_restore_vars"),                                                                                     &
           & lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%t_wml_lk)


      ! p_prog_wtr%h_ml_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('h_ml_lk', 'm', 'mixed-layer thickness', datatype_flt)
      grib2_desc = grib2_var(1, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfFirstFixedSurface", 1)
      CALL add_var( prog_list, vname_prefix//'h_ml_lk'//suffix, p_prog_wtr%h_ml_lk,                                   &
           & GRID_UNSTRUCTURED_CELL, ZA_MIX_LAYER, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF,     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in","mode_iniana",  &
           & "iau_restore_vars"),                                                                                     &
           & lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%h_ml_lk)


      ! p_prog_wtr%t_bot_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_bot_lk', 'K', 'temperature at the water-bottom sediment interface', &
        &          datatype_flt)
      grib2_desc = grib2_var(1, 2, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_bot_lk'//suffix, p_prog_wtr%t_bot_lk,                                    &
           & GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM_HALF, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF, &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in","mode_iniana",     &
           & "iau_restore_vars"),                                                                                        &
           & lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%t_bot_lk)


      ! p_prog_wtr%c_t_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('c_t_lk', '-', 'shape factor (temp. profile in lake thermocline)', &
        &          datatype_flt)
      grib2_desc = grib2_var(1, 2, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
        &           + t_grib2_int_key("typeOfSecondFixedSurface", 162)
      CALL add_var( prog_list, vname_prefix//'c_t_lk'//suffix, p_prog_wtr%c_t_lk,                                     &
           & GRID_UNSTRUCTURED_CELL, ZA_MIX_LAYER, cf_desc, grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF,     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in","mode_iniana",  &
           & "iau_restore_vars"),                                                                                     &
           & lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%c_t_lk)


      ! p_prog_wtr%t_b1_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('t_b1_lk', 'K',                                         &
        &          'temperature at the bottom of the upper layer of the sediments', &
        &          datatype_flt)
      grib2_desc = grib2_var(1, 2, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( prog_list, vname_prefix//'t_b1_lk'//suffix, p_prog_wtr%t_b1_lk, &
           & GRID_UNSTRUCTURED_CELL, ZA_SEDIMENT_BOTTOM_TW_HALF, cf_desc,           &
           & grib2_desc, ldims=shape2d, tlev_source=TLEV_NNOW_RCF,                  &
           & in_group=groups("iau_restore_vars"),                                   &
           & lopenacc=.TRUE. )   
      __acc_attach(p_prog_wtr%t_b1_lk)


      ! p_prog_wtr%h_b1_lk(nproma,nblks_c)
      cf_desc    = t_cf_var('h_b1_lk', 'm',                                         &
        &          'thickness of the upper layer of the sediments', datatype_flt)
      grib2_desc = grib2_var(1, 2, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      &           + t_grib2_int_key("typeOfSecondFixedSurface", 165)
      CALL add_var( prog_list, vname_prefix//'h_b1_lk'//suffix, p_prog_wtr%h_b1_lk,      &
           & GRID_UNSTRUCTURED_CELL, ZA_LAKE_BOTTOM, cf_desc, grib2_desc, ldims=shape2d, &
           & tlev_source=TLEV_NNOW_RCF,                                                  &
           & in_group=groups("iau_restore_vars"),                                        &
           & lopenacc=.TRUE. )
      __acc_attach(p_prog_wtr%h_b1_lk)

    ENDIF  ! llake

  END SUBROUTINE new_nwp_wtr_prog_list


  !-------------------------------------------------------------------------
  !>
  !! Allocation of components of diagnostic land state.
  !!
  !! Initialization of components with zero.
  !!
  SUBROUTINE new_nwp_lnd_diag_list( p_jg, kblks, listname, vname_prefix, &
    &                               diag_list, p_diag_lnd, l_smi)

    INTEGER,INTENT(IN) ::  kblks !< dimension sizes
    INTEGER,INTENT(IN) ::  p_jg !< patch id

    CHARACTER(len=*),INTENT(IN) :: listname, vname_prefix
    CHARACTER(LEN=2)            :: csfc

    TYPE(t_var_list_ptr),INTENT(INOUT) :: diag_list
    TYPE(t_lnd_diag),INTENT(INOUT) :: p_diag_lnd
    LOGICAL,         INTENT(IN)    :: l_smi   !< Flag. TRUE if computation 
                                              !< of soil moisture index desired

    ! Local variables
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape2d(2), shape3d_subs(3), shape3d_subsw(3)
    INTEGER :: jsfc          !< tile counter
    INTEGER :: ibits
    INTEGER :: datatype_flt

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

!--------------------------------------------------------------

    ibits = DATATYPE_PACK16 ! packing accuracy of horizontal slice

    ! predefined array shapes
    shape2d       = (/nproma, kblks            /)
    shape3d_subs  = (/nproma, kblks, ntiles_total /)
    shape3d_subsw = (/nproma, kblks, ntiles_total+ntiles_water /)


    !------------------------------
    ! Ensure that all pointers have a defined association status.
    !------------------------------
    NULLIFY(p_diag_lnd%qv_s, &
    &       p_diag_lnd%t_s, &
    &       p_diag_lnd%t_sk, &
    &       p_diag_lnd%t_seasfc, &
    &       p_diag_lnd%w_i, &
    &       p_diag_lnd%w_p, &
    &       p_diag_lnd%w_s, &
    &       p_diag_lnd%t_so, &
    &       p_diag_lnd%w_so, &
    &       p_diag_lnd%w_so_ice, &
    &       p_diag_lnd%smi, &
    &       p_diag_lnd%runoff_s, &
    &       p_diag_lnd%runoff_g, &
    &       p_diag_lnd%fr_seaice, &
    &       p_diag_lnd%qv_s_t, &
    &       p_diag_lnd%runoff_s_t, &
    &       p_diag_lnd%runoff_g_t, &
    &       p_diag_lnd%rstom, &
    &       p_diag_lnd%rstom_t, &
    &       p_diag_lnd%t_snow, &
    &       p_diag_lnd%rho_snow, &
    &       p_diag_lnd%w_snow, &
    &       p_diag_lnd%h_snow, &
    &       p_diag_lnd%h_snow_t, &
    &       p_diag_lnd%freshsnow, &
    &       p_diag_lnd%freshsnow_t, &
    &       p_diag_lnd%snowfrac, &
    &       p_diag_lnd%snowfrac_t, &
    &       p_diag_lnd%snowfrac_lc_t, &
    &       p_diag_lnd%snowfrac_lcu_t, &
    &       p_diag_lnd%hsnow_max, &
    &       p_diag_lnd%snow_age, &
    &       p_diag_lnd%qi_snowdrift_flx, &
    &       p_diag_lnd%t_snow_mult, &
    &       p_diag_lnd%rho_snow_mult, &
    &       p_diag_lnd%wliq_snow, &
    &       p_diag_lnd%wtot_snow, &
    &       p_diag_lnd%dzh_snow,  &
    &       p_diag_lnd%resid_wso, &
    &       p_diag_lnd%resid_wso_t, &
    &       p_diag_lnd%resid_wso_inst_t, &
    &       p_diag_lnd%condhf_ice, &
    &       p_diag_lnd%meltpot_ice)


    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(diag_list, TRIM(listname), patch_id=p_jg, &
      &          lrestart=.TRUE., model_type=get_my_process_name())

    !------------------------------

    ! & p_diag_lnd%qv_s(nproma,nblks_c)
    cf_desc    = t_cf_var('qv_s', 'kg kg-1', 'specific humidity at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)  &
      ! Quick hack: shortName.def should be revised, instead
      &           + t_grib2_int_key("scaleFactorOfFirstFixedSurface", 0)
    CALL add_var( diag_list, vname_prefix//'qv_s', p_diag_lnd%qv_s,          &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape2d,                                                  &
           & initval=0.001_wp,                                               &
           & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in", &
           &                 "mode_iau_fg_in","mode_iau_old_fg_in",          &
           &                 "mode_combined_in","mode_cosmo_in",             &
           &                 "mode_iniana"), lopenacc=.TRUE. )      
    __acc_attach(p_diag_lnd%qv_s)


    ! & p_diag_lnd%fr_seaice(nproma,nblks_c)
    cf_desc    = t_cf_var('fr_seaice', '1', 'fraction of sea ice', datatype_flt)
    grib2_desc = grib2_var(10, 2, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'fr_seaice', p_diag_lnd%fr_seaice,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                   &
           & initval=0._wp,                                                    &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_ana_in",&
           &                 "mode_iau_old_ana_in","mode_combined_in",         &
           &                 "mode_iniana","iau_restore_vars"), lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%fr_seaice)

    IF (is_coupled_to_ocean()) THEN
      ! & p_diag_lnd%condhf_ice(nproma,nblks_c)
      cf_desc    = t_cf_var('condhf_ice', 'W m-2', 'conductive heat flux at sea-ice bottom', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'condhf_ice', p_diag_lnd%condhf_ice,&
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=.TRUE., lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%condhf_ice)

      ! & p_diag_lnd%meltpot_ice(nproma,nblks_c)
      cf_desc    = t_cf_var('meltpot_ice', 'W m-2', 'melt potential at sea-ice top', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'meltpot_ice', p_diag_lnd%meltpot_ice,&
            & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
            & ldims=shape2d, lrestart=.FALSE., lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%meltpot_ice)
    END IF

    ! & p_diag_lnd%qv_s_t(nproma,nblks_c,ntiles_total+ntiles_water)
    cf_desc    = t_cf_var('qv_s_t', 'kg kg-1', 'specific humidity at the surface', datatype_flt)
    grib2_desc = grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'qv_s_t', p_diag_lnd%qv_s_t,        &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE.,         &
           & loutput=.FALSE.,                                                  &
           & in_group=groups("iau_restore_vars"),                              &
           & initval=0.001_wp, lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%qv_s_t)

    ! fill the separate variables belonging to the container qv_s_t
    ALLOCATE(p_diag_lnd%qv_st_ptr(ntiles_total+ntiles_water))
    DO jsfc = 1,ntiles_total+ntiles_water
      NULLIFY(p_diag_lnd%qv_st_ptr(jsfc)%p_2d, p_diag_lnd%qv_st_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'qv_s_t',                       &
             & vname_prefix//'qv_s_t_'//ADJUSTL(TRIM(csfc)),                 &
             & p_diag_lnd%qv_st_ptr(jsfc)%p_2d,                              &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
             & t_cf_var('qv_s_t_'//csfc, '', '', datatype_flt),            &
             & grib2_var(0, 1, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
             & ref_idx=jsfc,                                                 &
             & ldims=shape2d,                                                &
             & var_class=CLASS_TILE,                                         &
             & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
    ENDDO

    IF ( atm_phy_nwp_config(p_jg)%inwp_surface > 0) THEN

    ! & p_diag_lnd%t_s(nproma,nblks_c)
    cf_desc    = t_cf_var('t_s', 'K', 'weighted temperature of ground surface', datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_s', p_diag_lnd%t_s,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE., lopenacc=.TRUE.      )
    __acc_attach(p_diag_lnd%t_s)


    ! & p_diag_lnd%t_sk(nproma,nblks_c)
    cf_desc    = t_cf_var('t_sk', 'K', 'skin temperature', datatype_flt)
    grib2_desc = grib2_var(0, 0, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_sk', p_diag_lnd%t_sk,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE., in_group=              &
         & groups("dwd_fg_sfc_vars", "mode_iau_fg_in", "mode_dwd_fg_in", "mode_combined_in"), &
         & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%t_sk)


    ! & p_diag_lnd%t_seasfc(nproma,nblks_c)
    cf_desc    = t_cf_var('t_seasfc', 'K', 'sea surface temperature', datatype_flt)
    grib2_desc = grib2_var(10, 3, 0, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_seasfc', p_diag_lnd%t_seasfc,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                       &
         & lmiss=.TRUE., missval=0.0_wp,                                         &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), & 
         & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in",                  &
         &                 "mode_iau_ana_in","mode_iau_old_ana_in"),             &
         & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%t_seasfc)


    ! & p_diag_lnd%w_i(nproma,nblks_c)
    cf_desc     = t_cf_var('w_i', 'm H2O', 'weighted water content of interception water', &
         &                datatype_flt)
    new_cf_desc = t_cf_var('w_i', 'kg m-2', 'weighted water content of interception water', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 0, 13, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_i', p_diag_lnd%w_i,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",       &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in","mode_combined_in", &
         &                 "mode_cosmo_in","mode_iniana"),                       &
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
         & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%w_i)


    IF (itype_interception == 2) THEN

      ! & p_diag_lnd%w_p(nproma,nblks_c)
      cf_desc     = t_cf_var('w_p', 'm H2O', 'water content of pond interception water', &
           &                datatype_flt)
      new_cf_desc = t_cf_var('w_p', 'kg m-2', 'water content of pond interception water', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'w_p', p_diag_lnd%w_p,                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc,  grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
           & in_group=groups("land_vars"),                                         &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
           & lopenacc=.TRUE. )
      __acc_attach(p_diag_lnd%w_p)

      ! & p_diag_lnd%w_s(nproma,nblks_c)
      cf_desc     = t_cf_var('w_s', 'm H2O', 'water content of interception snow', &
           &                datatype_flt)
      new_cf_desc = t_cf_var('w_s', 'kg m-2', 'water content of interception snow', &
           &                datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'w_s', p_diag_lnd%w_s,                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,              &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
           & in_group=groups("land_vars"),                                         &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
           & lopenacc=.TRUE. )
      __acc_attach(p_diag_lnd%w_s)

    END IF  ! itype_interception == 2

    IF (itype_snowevap == 3) THEN
      ! maximum snow depth reached within current snow-cover period
      cf_desc    = t_cf_var('hsnow_max', 'm', 'maximum snow depth', datatype_flt)
      grib2_desc = grib2_var(0, 1, 235, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'hsnow_max', p_diag_lnd%hsnow_max,                &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
             & ldims=shape2d, lrestart=.TRUE.,                                   &
             & in_group=groups("dwd_fg_sfc_vars","mode_iau_fg_in",               &
             &                "mode_dwd_fg_in","mode_iniana","mode_combined_in", &
             &                "iau_restore_vars"),                               &
             & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%hsnow_max)

      ! duration of current snow-cover peiod
      cf_desc    = t_cf_var('snow_age', 'd', 'duration of snow cover', datatype_flt)
      grib2_desc = grib2_var(0, 1, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'snow_age', p_diag_lnd%snow_age,                  &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
             & ldims=shape2d, lrestart=.TRUE.,                                   &
             & in_group=groups("dwd_fg_sfc_vars","mode_iau_fg_in",               &
             &                "mode_dwd_fg_in","mode_iniana","mode_combined_in", &
             &                "iau_restore_vars"),                               &
             & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%snow_age)

      ! cloud-ice flux from drifting snow
      cf_desc    = t_cf_var('qi_snowdrift_flx', 'kg m-2 s', &
          & 'upward cloud-ice surface flux from drifting snow', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, 'qi_snowdrift_flx', p_diag_lnd%qi_snowdrift_flx,              &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
             & ldims=shape2d, lrestart=.FALSE.,                                  &
             & in_group=groups('land_vars'),                                     &
             & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%qi_snowdrift_flx)
    ENDIF

    ! & p_diag_lnd%t_so(nproma,nlev_soil+1,nblks_c)
    cf_desc    = t_cf_var('t_so', 'K', 'weighted soil temperature (main level)', &
         &                datatype_flt)
    grib2_desc = grib2_var(2, 3, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_so', p_diag_lnd%t_so,              &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND_P1, cf_desc, grib2_desc,  &
         & ldims=(/nproma,nlev_soil+1,kblks/),                                   &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars","dwd_fg_sfc_vars",                        &
         &                 "mode_dwd_fg_in","mode_iau_fg_in","mode_iau_old_fg_in",&
         &                 "mode_dwd_ana_in","mode_iau_ana_in",                  &
         &                 "mode_iau_old_ana_in","mode_combined_in",             &
         &                 "mode_cosmo_in","mode_iniana"), lopenacc=.TRUE.       )
      __acc_attach(p_diag_lnd%t_so)


    ! & p_diag_lnd%w_so(nproma,nlev_soil,nblks_c)
    cf_desc     = t_cf_var('w_so', 'm H20', 'total water content (ice + liquid water)', datatype_flt)
    new_cf_desc = t_cf_var('w_so', 'kg m-2', 'total water content (ice + liquid water)', datatype_flt)
    grib2_desc  = grib2_var(2, 3, 20, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_so', p_diag_lnd%w_so,              &
         & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND, cf_desc, grib2_desc,     &
         & ldims=(/nproma,nlev_soil,kblks/),                                     &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_ana_in",      &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in","mode_iau_ana_in", &
         &                 "mode_iau_old_ana_in","mode_combined_in",             &
         &                 "mode_cosmo_in"),                                     &
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ),&
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
         & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%w_so)

    ! & p_diag_lnd%w_so_ice(nproma,nlev_soil,nblks_c)
    cf_desc      = t_cf_var('w_so_ice', 'm H20',   'ice content', datatype_flt)
    new_cf_desc  = t_cf_var('w_so_ice', 'kg m-2', 'ice content', datatype_flt)
    grib2_desc   = grib2_var(2, 3, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)      
    CALL add_var( diag_list, vname_prefix//'w_so_ice',                           &
         & p_diag_lnd%w_so_ice, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,     &
         & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil,kblks/),                &
         & lrestart=.FALSE., loutput=.TRUE.,                                     &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",       &
         &                 "mode_iau_fg_in","mode_iau_old_fg_in","mode_iniana"),&
         & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ),&
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
         & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%w_so_ice)

    IF (l_smi) THEN
      ! & p_diag_lnd%smi(nproma,nlev_soil,nblks_c)
      cf_desc      = t_cf_var('smi', '--',   'soil moisture index', datatype_flt)
      grib2_desc   = grib2_var(2, 3, 200, ibits, GRID_UNSTRUCTURED, GRID_CELL)      
      CALL add_var( diag_list, vname_prefix//'smi',                                &
           & p_diag_lnd%smi, GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,          &
           & cf_desc, grib2_desc, ldims=(/nproma,nlev_soil,kblks/),                &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & in_group=groups(groups_smi),                                      &
           & hor_interp=create_hor_interp_metadata(hor_intp_type=HINTP_TYPE_LONLAT_NNB ), &
           & l_pp_scheduler_task=TASK_COMPUTE_SMI,                                 &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%smi)
    ENDIF  ! l_smi


    ! & p_diag_lnd%runoff_s(nproma,nblks_c)
    cf_desc    = t_cf_var('runoff_s', 'kg m-2', &
         &                'weighted surface water runoff', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s', p_diag_lnd%runoff_s,        &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S, cf_desc, grib2_desc,       &
           & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                      &
           & isteptype=TSTEP_ACCUM,                                                &
           & hor_interp=create_hor_interp_metadata(                                &
           &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                              &
           &    fallback_type=HINTP_TYPE_LONLAT_NNB),                                &
           & initval=0._wp, resetval=0._wp,                              &
           & action_list=actions(new_action(ACTION_RESET,runoff_interval(p_jg))),  &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%runoff_s)


    ! & p_diag_lnd%runoff_g(nproma,nblks_c)
    cf_desc    = t_cf_var('runoff_g', 'kg m-2', &
         &                'weighted soil/ground water runoff', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g', p_diag_lnd%runoff_g,        &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
           & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                      &
           & isteptype=TSTEP_ACCUM,                                                &
           & hor_interp=create_hor_interp_metadata(                                &
           &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                              &
           &    fallback_type=HINTP_TYPE_LONLAT_NNB),                              &
           & initval=0._wp, resetval=0._wp,                              &
           & action_list=actions(new_action(ACTION_RESET,runoff_interval(p_jg))),  &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%runoff_g)

    IF (var_in_output(p_jg)%res_soilwatb) THEN
      ! & p_diag_lnd%resid_wso(nproma,nblks_c)
      ! resid_wso is a budget considering lateral fluxes and sources / sinks located between
      ! Earth surface and last hydrological soil layer -> best suited: vertical axis of w_so
      cf_desc    = t_cf_var('resid_wso', 'kg m-2', &
          &                'residuum of the mass content of water in soil', datatype_flt)
      grib2_desc = grib2_var(2, 3, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'resid_wso', p_diag_lnd%resid_wso,     &
            & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
            & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
            & isteptype=TSTEP_ACCUM,                                                &
            & hor_interp=create_hor_interp_metadata(                                &
            &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                              &
            &    fallback_type=HINTP_TYPE_LONLAT_NNB),                              &
            & initval=0._wp, resetval=0._wp,                                        &
            & action_list=actions(new_action(ACTION_RESET,runoff_interval(p_jg))),  &
            & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%resid_wso)
    ENDIF

    ! & p_diag_lnd%runoff_s_inst_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_s_inst_t', 'kg m-2', &
         &                'surface water runoff; instantaneous value', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s_inst_t', p_diag_lnd%runoff_s_inst_t,  &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S, cf_desc, grib2_desc,       &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,    &
           & in_group=groups("iau_init_vars"),                                             &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%runoff_s_inst_t)

    ! & p_diag_lnd%runoff_g_inst_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_g_inst_t', 'kg m-2', &
         &                'soil/ground water runoff; instantaneous value', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g_inst_t', p_diag_lnd%runoff_g_inst_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,      &
           & in_group=groups("iau_init_vars"),                                               &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%runoff_g_inst_t) 

    IF (var_in_output(p_jg)%res_soilwatb) THEN    
      ! & p_diag_lnd%resid_wso_inst_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('resid_wso_inst_t', 'kg m-2', &
          &                'residuum of the mass content of water in soil; instanteaneous value', datatype_flt)
      grib2_desc = grib2_var(2, 3, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'resid_wso_inst_t', p_diag_lnd%resid_wso_inst_t, &
            & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,                 &
            & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       &
            & in_group=groups("iau_init_vars"),                                               &
            & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%resid_wso_inst_t) 
    ENDIF

    ! & p_diag_lnd%runoff_s_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_s_t', 'kg m-2', &
         &                'surface water runoff', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_s_t', p_diag_lnd%runoff_s_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S, cf_desc, grib2_desc,       &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., &
           & initval=0._wp, resetval=0._wp,                              &
           & action_list=actions(new_action(ACTION_RESET,runoff_interval(p_jg))), &
           & in_group=groups("iau_init_vars"),                                    &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%runoff_s_t)

    ! & p_diag_lnd%runoff_g_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('runoff_g_t', 'kg m-2', &
         &                'soil/ground water runoff', datatype_flt)
    grib2_desc = grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'runoff_g_t', p_diag_lnd%runoff_g_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
           & ldims=shape3d_subsw, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE., &
           & initval=0._wp, resetval=0._wp,                              &
           & action_list=actions(new_action(ACTION_RESET,runoff_interval(p_jg))), &
           & in_group=groups("iau_init_vars"),                                    &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%runoff_g_t)

    IF (var_in_output(p_jg)%res_soilwatb) THEN
      ! & p_diag_lnd%resid_wso_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('resid_wso_t', 'kg m-2', &
      &                'residuum of the mass content of water in soil', datatype_flt)
      grib2_desc = grib2_var(2, 3, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'resid_wso_t', p_diag_lnd%resid_wso_t, &
            & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G, cf_desc, grib2_desc,       &
            & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,&
            & initval=0._wp, resetval=0._wp,                              &
            & action_list=actions(new_action(ACTION_RESET,runoff_interval(p_jg))), &
            & in_group=groups("iau_init_vars"),                                    &
            & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%resid_wso_t)
    ENDIF

    ! fill the separate variables belonging to the container runoff_s
    ALLOCATE(p_diag_lnd%runoff_s_ptr(ntiles_total + ntiles_water))
    DO jsfc = 1,ntiles_total + ntiles_water
      NULLIFY(p_diag_lnd%runoff_s_ptr(jsfc)%p_2d, p_diag_lnd%runoff_s_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'runoff_s_t',                       &
               & vname_prefix//'runoff_s_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%runoff_s_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_S,                      &
               & t_cf_var('runoff_s_t_'//csfc, '', '', datatype_flt),            &
               & grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
               & ref_idx=jsfc,                                                   &
               & ldims=shape2d,                                                  &
!               & isteptype=TSTEP_ACCUM,                                          &
               & var_class=CLASS_TILE_LAND )
    END DO

    ! fill the separate variables belonging to the container runoff_g
    ALLOCATE(p_diag_lnd%runoff_g_ptr(ntiles_total + ntiles_water))
    DO jsfc = 1,ntiles_total + ntiles_water
      NULLIFY(p_diag_lnd%runoff_g_ptr(jsfc)%p_2d, p_diag_lnd%runoff_g_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'runoff_g_t',                       &
               & vname_prefix//'runoff_g_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%runoff_g_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G,                      &
               & t_cf_var('runoff_g_t_'//csfc, '', '', datatype_flt),            &
               & grib2_var(2, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL),        &
               & ref_idx=jsfc,                                                   &
               & ldims=shape2d,                                                  &
!               & isteptype=TSTEP_ACCUM,                                          &
               & var_class=CLASS_TILE_LAND )
    END DO


    IF (var_in_output(p_jg)%res_soilwatb) THEN
      ! fill the separate variables belonging to the container resid_wso
      ALLOCATE(p_diag_lnd%resid_wso_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_diag_lnd%resid_wso_ptr(jsfc)%p_2d, p_diag_lnd%resid_wso_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'resid_wso_t',                     &
                & vname_prefix//'resid_wso_t_'//ADJUSTL(TRIM(csfc)),              &
                & p_diag_lnd%resid_wso_ptr(jsfc)%p_2d,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_RUNOFF_G,                      &
                & t_cf_var('resid_wso_t_'//csfc, '', '', datatype_flt),           &
                & grib2_var(2, 3, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
                & ref_idx=jsfc,                                                   &
                & ldims=shape2d,                                                  &
  !               & isteptype=TSTEP_ACCUM,                                          &
                & var_class=CLASS_TILE_LAND )
       END DO
    ENDIF
    
    IF (var_in_output(p_jg)%snow_melt) THEN
        NULLIFY(p_diag_lnd%snow_melt,p_diag_lnd%snow_melt_flux_t)
        
        ! & lnd_diag%snow_melt(nproma,nblks_c)
        cf_desc     = t_cf_var('surface_snow_melt_amount', 'kg m-2', &
          &           'snow melt amount', datatype_flt)
        grib2_desc  = grib2_var(0, 1, 16, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( diag_list, vname_prefix//'snow_melt', p_diag_lnd%snow_melt, &
             & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,          & 
             & ldims=shape2d, lrestart=.TRUE., loutput=.TRUE.,                    &
             & isteptype=TSTEP_ACCUM,                                             &
             & initval=0._wp, resetval=0._wp, lopenacc=.TRUE.,                    &
             & action_list=actions(new_action(ACTION_RESET,melt_interval(p_jg))) )
        __acc_attach(p_diag_lnd%snow_melt)

        ! & lnd_diag%snow_melt_flux_t(nproma,nblks_c,ntiles_total)
        cf_desc    = t_cf_var('surface_snow_melt_flux_t', 'kg m-2 s-1', &
         &                'snow melt flux tile', datatype_flt)
        grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
        CALL add_var( diag_list, vname_prefix//'snow_melt_flux_t',                     &
           & p_diag_lnd%snow_melt_flux_t, GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, &
           & grib2_desc, ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,      &
           & loutput=.FALSE., initval=0._wp, resetval=0._wp, lopenacc=.TRUE.)
        __acc_attach(p_diag_lnd%snow_melt_flux_t)

        ! fill the separate variables belonging to the container snow_melt
        ALLOCATE(p_diag_lnd%snow_melt_ptr(ntiles_total))
        DO jsfc = 1,ntiles_total
          NULLIFY(p_diag_lnd%snow_melt_ptr(jsfc)%p_2d,                            &
          &       p_diag_lnd%snow_melt_ptr(jsfc)%p_3d)
          WRITE(csfc,'(i2)') jsfc 
          CALL add_ref( diag_list, vname_prefix//'snow_melt_flux_t',              &
                & vname_prefix//'snow_melt_flux_t_'//ADJUSTL(TRIM(csfc)),         &
                & p_diag_lnd%snow_melt_ptr(jsfc)%p_2d,                            &
                & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_LAND,                    &
                & t_cf_var('snow_melt_flux_t_'//csfc, '', '', datatype_flt),      &
                & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),  &
                & ref_idx=jsfc,                                                   &
                & ldims=shape2d,                                                  &
                & var_class=CLASS_TILE_LAND )
        END DO
    ENDIF

    ! & p_diag_lnd%rstom(nproma,nblks_c)
    cf_desc    = t_cf_var('rstom', 's m-1','stomatal resistance', datatype_flt)
    grib2_desc = grib2_var(2, 0, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rstom', p_diag_lnd%rstom,             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE.,                                     &
           & loutput=.TRUE., lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%rstom)


    ! & p_diag_lnd%rstom_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('rstom_t', 's m-1','tile based stomatal resistance', &
      &          datatype_flt)
    grib2_desc = grib2_var(2, 0, 195, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rstom_t', p_diag_lnd%rstom_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,        &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,        &
           & loutput=.FALSE., lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%rstom_t) 

    IF (itype_trvg == 3) THEN ! extended plant evaporation scheme
      ! & p_diag_lnd%plantevap(nproma,nblks_c)
      cf_desc    = t_cf_var('plantevap', 'kg m-2', &
         &       'function of time-integrated plant evaporation', datatype_flt)
      grib2_desc = grib2_var(2, 0, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'plantevap', p_diag_lnd%plantevap,   &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d,in_group=groups("dwd_fg_sfc_vars", "mode_iau_fg_in",   &
           & "mode_dwd_fg_in","mode_iniana","mode_combined_in"), lrestart=.FALSE.,&
           & loutput=.TRUE., lopenacc=.TRUE. )
      __acc_attach(p_diag_lnd%plantevap)

      ! & p_diag_lnd%plantevap_t(nproma,nblks_c,ntiles_total)
      cf_desc    = t_cf_var('plantevap_t', 'kg m-2', &
         &       'function of time-integrated plant evaporation', datatype_flt)
      grib2_desc = grib2_var(2, 0, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'plantevap_t', p_diag_lnd%plantevap_t,  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,                &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,                &
           & in_group=groups("iau_restore_vars"),                                    &
           & loutput=.FALSE.,                                                        &
           & lopenacc=.TRUE.)
      __acc_attach(p_diag_lnd%plantevap_t)

      ! fill the separate variables belonging to the container plantevap
      ALLOCATE(p_diag_lnd%plantevap_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_diag_lnd%plantevap_ptr(jsfc)%p_2d, p_diag_lnd%plantevap_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'plantevap_t',                    &
               & vname_prefix//'plantevap_t_'//ADJUSTL(TRIM(csfc)),              &
               & p_diag_lnd%plantevap_ptr(jsfc)%p_2d,                            &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                             &
               & t_cf_var('plantevap_t_'//csfc, '', '', datatype_flt),           &
               & grib2_var(2, 0, 198, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
               & ref_idx=jsfc,                                                   &
               & ldims=shape2d,                                                  &
               & var_class=CLASS_TILE_LAND,in_group=groups("land_tile_vars",     &
               & "dwd_fg_sfc_vars_t") )
      END DO
    ENDIF

    ! & p_diag_lnd%t_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('t_snow', 'K', 'weighted temperature of the snow-surface', &
         &                datatype_flt)
    grib2_desc = grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'t_snow', p_diag_lnd%t_snow,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                    &
         & in_group=groups("land_vars", "snow_vars","dwd_fg_sfc_vars",         &
         &                 "mode_dwd_ana_in","mode_iau_fg_in",                 &
         &                 "mode_iau_old_fg_in","mode_combined_in",            &
         &                 "mode_cosmo_in","mode_iniana"),                     &
         & lopenacc=.TRUE.                                                     )
    __acc_attach(p_diag_lnd%t_snow)


    ! & p_diag_lnd%w_snow(nproma,nblks_c)
    cf_desc     = t_cf_var('w_snow', 'm H2O',   'weighted water eqivalent of snow', datatype_flt)
    new_cf_desc = t_cf_var('w_snow', 'kg m-2', 'weighted water equivalent of snow', datatype_flt)
    grib2_desc  = grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'w_snow', p_diag_lnd%w_snow,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,           &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                    &
         & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_dwd_fg_in",     &
         &                 "mode_iau_old_ana_in",                              &
         &                 "mode_combined_in","mode_cosmo_in","mode_iniana"),  &
         & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),  &
         & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%w_snow)


    ! & p_diag_lnd%rho_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('rho_snow', 'kg m-3', 'weighted snow density', datatype_flt)
    grib2_desc = grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'rho_snow', p_diag_lnd%rho_snow,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,  cf_desc, grib2_desc,             &
         & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                      &
         & in_group=groups("land_vars", "snow_vars","dwd_fg_sfc_vars",           &
         &                 "mode_dwd_fg_in","mode_iau_fg_in",                    &
         &                 "mode_iau_old_ana_in","mode_combined_in",             &
         &                 "mode_cosmo_in","mode_iniana"), lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%rho_snow)


    ! & p_diag_lnd%h_snow(nproma,nblks_c)
    cf_desc    = t_cf_var('h_snow', 'm', 'weighted snow depth', datatype_flt)
    grib2_desc = grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'h_snow', p_diag_lnd%h_snow,        &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                  &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in", &
           &                 "mode_iau_ana_in","mode_iau_old_ana_in",              &
           &                 "mode_combined_in","mode_iniana","iau_restore_vars"), &
           &                 lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%h_snow)

    ! & p_diag_lnd%h_snow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('h_snow_t', 'm', 'snow height', datatype_flt)
    grib2_desc = grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'h_snow_t', p_diag_lnd%h_snow_t,    &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,          &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,          &
           & in_group=groups("iau_restore_vars"),                              &
           & loutput=.FALSE.,                                                  &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%h_snow_t)

    ! fill the separate variables belonging to the container h_snow
    ALLOCATE(p_diag_lnd%h_snow_ptr(ntiles_total))
      DO jsfc = 1,ntiles_total
        NULLIFY(p_diag_lnd%h_snow_ptr(jsfc)%p_2d, p_diag_lnd%h_snow_ptr(jsfc)%p_3d)
        WRITE(csfc,'(i2)') jsfc 
        CALL add_ref( diag_list, vname_prefix//'h_snow_t',                     &
               & vname_prefix//'h_snow_t_'//ADJUSTL(TRIM(csfc)),               &
               & p_diag_lnd%h_snow_ptr(jsfc)%p_2d,                             &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('h_snow_t_'//csfc, '', '', datatype_flt),          &
               & grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
               & ref_idx=jsfc,                                                 &
               & ldims=shape2d,                                                &
               & var_class=CLASS_TILE_LAND,                                    &
               & in_group=groups("land_tile_vars","dwd_fg_sfc_vars_t") )
      ENDDO


    ! & p_diag_lnd%freshsnow(nproma,nblks_c)
    cf_desc    = t_cf_var('freshsnow', '1', &
           & 'weighted indicator for age of snow in top of snow layer', datatype_flt)
    grib2_desc = grib2_var(0, 1, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'freshsnow', p_diag_lnd%freshsnow,     &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_ana_in","mode_iau_fg_in",&
           &                 "mode_iau_ana_in","mode_iau_old_ana_in",             &
           &                 "mode_combined_in","mode_cosmo_in","mode_iniana"),   &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%freshsnow)


    ! & p_diag_lnd%freshsnow_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('freshsnow_t', '1', &
         &                'indicator for age of snow in top of snow layer', datatype_flt)
    grib2_desc = grib2_var(0, 1, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'freshsnow_t', p_diag_lnd%freshsnow_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,             &
           & in_group=groups("iau_restore_vars"),                                 &
           & loutput=.FALSE.,                                                     &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%freshsnow_t)

    ! fill the separate variables belonging to the container freshsnow
    ALLOCATE(p_diag_lnd%freshsnow_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%freshsnow_ptr(jsfc)%p_2d, p_diag_lnd%freshsnow_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'freshsnow_t',                   &
               & vname_prefix//'freshsnow_t_'//ADJUSTL(TRIM(csfc)),           &
               & p_diag_lnd%freshsnow_ptr(jsfc)%p_2d,                         &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                          &
               & t_cf_var('freshsnow_t_'//csfc, '', '', datatype_flt),      &
               & grib2_var(0, 1, 203, ibits, GRID_UNSTRUCTURED, GRID_CELL),   &
               & ref_idx=jsfc,                                                &
               & ldims=shape2d,                                               &
               & var_class=CLASS_TILE_LAND,                                   &
               & in_group=groups("dwd_fg_sfc_vars_t") )
    END DO


    ! & p_diag_lnd%snowfrac(nproma,nblks_c)
    cf_desc    = t_cf_var('snowfrac', '1' , 'snow-cover fraction', datatype_flt)
    new_cf_desc= t_cf_var('snowfrac', '% ', 'snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac', p_diag_lnd%snowfrac,       &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("land_vars"),                                        &
           & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),    &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%snowfrac)


    ! & p_diag_lnd%snowfrac_lc(nproma,nblks_c)
    cf_desc    = t_cf_var('snowfrac_lc', '1 ', 'snow-cover fraction', DATATYPE_FLT32)
    new_cf_desc= t_cf_var('snowfrac_lc', '% ', 'snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(0, 1, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_lc', p_diag_lnd%snowfrac_lc, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,             &
           & ldims=shape2d, lrestart=.FALSE., loutput=.TRUE.,                     &
           & in_group=groups("land_vars","dwd_fg_sfc_vars","mode_iau_fg_in"),     &
           & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc),    &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%snowfrac_lc)


    ! local snow-cover referring to a tile; zero on snow-free tiles, one (except for very low snow depths) on snow-tile
    !
    ! & p_diag_lnd%snowfrac_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_t', '1 ', 'local tile-based snow-cover fraction', datatype_flt)
    new_cf_desc= t_cf_var('snowfrac_t', '% ', 'local tile-based snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_t', p_diag_lnd%snowfrac_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,           &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,           &
           & in_group=groups("iau_restore_vars"),                               &
           & loutput=.FALSE.,                                                   &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%snowfrac_t)

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%snowfrac_ptr(jsfc)%p_2d, p_diag_lnd%snowfrac_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'snowfrac_t',                     &
               & vname_prefix//'snowfrac_t_'//ADJUSTL(TRIM(csfc)),             &
               & p_diag_lnd%snowfrac_ptr(jsfc)%p_2d,                           &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('snowfrac_t_'//csfc, '', '', datatype_flt),        &
               & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
               & ref_idx=jsfc,                                                 &
               & ldims=shape2d,                                                &
               & var_class=CLASS_TILE_LAND,                                    &
               & in_group=groups("land_tile_vars"),                            &
               & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )
    END DO


    ! snow-cover fraction referring to a landuse tile of a grid point; has the same value on a snow tile
    ! and the corresponding snow-free tile
    !
    ! & p_diag_lnd%snowfrac_lc_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_lc_t', '1 ', 'tile-based snow-cover fraction', datatype_flt)
    new_cf_desc= t_cf_var('snowfrac_lc_t', '% ', 'tile-based snow-cover fraction', DATATYPE_FLT32)
    grib2_desc = grib2_var(0, 1, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_lc_t', p_diag_lnd%snowfrac_lc_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE.,            &
           & in_group=groups("iau_restore_vars"),                                &
           & loutput=.FALSE.,                                                    &
           & lopenacc=.TRUE. )
    __acc_attach(p_diag_lnd%snowfrac_lc_t)

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_lc_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_2d, p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'snowfrac_lc_t',                  &
               & vname_prefix//'snowfrac_lc_t_'//ADJUSTL(TRIM(csfc)),          &
               & p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_2d,                        &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('snowfrac_lc_t_'//csfc, '', '', datatype_flt),     &
               & grib2_var(0, 1, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL),     &
               & ref_idx=jsfc,                                                 &
               & ldims=shape2d,                                                &
               & var_class=CLASS_TILE_LAND,                                    &
               & in_group=groups("land_tile_vars", "dwd_fg_sfc_vars_t"),       &
               & post_op=post_op(POST_OP_SCALE, arg1=100._wp, new_cf=new_cf_desc) )
    END DO


    ! same as snowfrac_lc_t if no melting occurs on a given grid point; in the case of melting snow,
    ! this quantity carries the unmodifed snow-cover fraction diagnosed from the snow depth, whereas
    ! snowfrac_lc_t carries an artificially reduced snow-cover fraction
    !
    ! & p_diag_lnd%snowfrac_lcu_t(nproma,nblks_c,ntiles_total)
    cf_desc    = t_cf_var('snowfrac_lcu_t', '1 ', 'tile-based snow-cover fraction', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( diag_list, vname_prefix//'snowfrac_lcu_t', p_diag_lnd%snowfrac_lcu_t, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
           & ldims=shape3d_subs, lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,  &
           & lopenacc=.TRUE.)
    __acc_attach(p_diag_lnd%snowfrac_lcu_t)

    ! fill the separate variables belonging to the container snowfrac
    ALLOCATE(p_diag_lnd%snowfrac_lcu_ptr(ntiles_total))
    DO jsfc = 1,ntiles_total
      NULLIFY(p_diag_lnd%snowfrac_lcu_ptr(jsfc)%p_2d, p_diag_lnd%snowfrac_lcu_ptr(jsfc)%p_3d)
      WRITE(csfc,'(i2)') jsfc 
      CALL add_ref( diag_list, vname_prefix//'snowfrac_lcu_t',                 &
               & vname_prefix//'snowfrac_lcu_t_'//ADJUSTL(TRIM(csfc)),         &
               & p_diag_lnd%snowfrac_lc_ptr(jsfc)%p_2d,                        &
               & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
               & t_cf_var('snowfrac_lcu_t_'//csfc, '', '', datatype_flt),      &
               & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
               & ref_idx=jsfc,                                                 &
               & ldims=shape2d, lrestart=.TRUE.,                               &
               & var_class=CLASS_TILE_LAND )
    END DO


    IF (lmulti_snow .OR. l2lay_rho_snow) THEN

      ! & p_diag_lnd%rho_snow_mult(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('rho_snow_mult', 'kg m-3', 'weighted snow density', datatype_flt)
      grib2_desc = grib2_var(0, 1, 61, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'rho_snow_mult',                      &
           & p_diag_lnd%rho_snow_mult, GRID_UNSTRUCTURED_CELL, ZA_SNOW,            &
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars"),               &
           & lopenacc=.TRUE. )
      __acc_attach(p_diag_lnd%rho_snow_mult)

    ENDIF

    IF (lmulti_snow) THEN

      ! & p_diag_lnd%t_snow_mult(nproma,nlev_snow+1,nblks_c)
      cf_desc    = t_cf_var('t_snow_mult', 'K', 'weighted temperature of the snow', datatype_flt)
      grib2_desc = grib2_var(0, 0, 18, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'t_snow_mult', p_diag_lnd%t_snow_mult, &
       & GRID_UNSTRUCTURED_CELL, ZA_SNOW_HALF, cf_desc, grib2_desc,                 &
       & ldims=(/nproma,nlev_snow+1,kblks/),                                        &
       & lrestart=.FALSE., loutput=.TRUE.,                                          &
       & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",       &
       &                 "mode_iau_old_fg_in","multisnow_vars") )


      ! & p_diag_lnd%wliq_snow(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('wliq_snow', 'm H2O', 'weighted liquid water content in snow', &
        &                   datatype_flt)
      new_cf_desc= t_cf_var('wliq_snow', 'kg m-2', 'weighted liquid water content in snow', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 210, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'wliq_snow',                          &
           & p_diag_lnd%wliq_snow, GRID_UNSTRUCTURED_CELL, ZA_SNOW,                & 
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars", "snow_vars") )


      ! & p_diag_lnd%wtot_snow(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('wtot_snow', 'm H2O', 'weighted total water content in snow', &
           &                datatype_flt)
      new_cf_desc= t_cf_var('wtot_snow', 'kg m-2', 'weighted total water content in snow', &
        &                   datatype_flt)
      grib2_desc = grib2_var(0, 1, 60, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'wtot_snow',                          &
           & p_diag_lnd%wtot_snow, GRID_UNSTRUCTURED_CELL, ZA_SNOW,                &
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & post_op=post_op(POST_OP_SCALE, arg1=1000._wp, new_cf=new_cf_desc),    & 
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars","snow_vars") )


      ! & p_diag_lnd%dzh_snow(nproma,nlev_snow,nblks_c)
      cf_desc    = t_cf_var('dzh_snow', 'm', &
           &                'weighted layer thickness between half levels in snow', &
           &                datatype_flt)
      grib2_desc = grib2_var(0, 1, 11, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( diag_list, vname_prefix//'dzh_snow',                           &
           & p_diag_lnd%dzh_snow, GRID_UNSTRUCTURED_CELL, ZA_SNOW,                 &
           & cf_desc, grib2_desc, ldims=(/nproma, nlev_snow, kblks/),              &
           & lrestart=.FALSE., loutput=.TRUE.,                                     &
           & in_group=groups("dwd_fg_sfc_vars","mode_dwd_fg_in","mode_iau_fg_in",  &
           &                 "mode_iau_old_fg_in","multisnow_vars","snow_vars"))

    ENDIF  ! lmulti_snow

    ENDIF  ! inwp_surface > 0

    

  END SUBROUTINE  new_nwp_lnd_diag_list


END MODULE mo_nwp_lnd_state
