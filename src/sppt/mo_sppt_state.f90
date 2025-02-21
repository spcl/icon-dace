!
! Construction of a data object which is used to store fields used for SPPT
! (Stochastic Perturbation of Physics Tendencies)
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

MODULE mo_sppt_state

  USE mo_impl_constants,          ONLY: SUCCESS, MAX_CHAR_LENGTH, HINTP_TYPE_LONLAT_BCTR, &
   &                                    HINTP_TYPE_LONLAT_RBF,VINTP_METHOD_LIN
  USE mo_exception,               ONLY: message, finish
  USE mo_sppt_types,              ONLY: t_sppt
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_config,             ONLY: n_dom
  USE mo_gribout_config,          ONLY: gribout_config
  USE mo_model_domain,            ONLY: t_patch
  USE mo_var_list,                ONLY: add_var, t_var_list_ptr
  USE mo_var_list_register,       ONLY: vlr_add, vlr_del
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_CELL,GRID_CELL
  USE mo_zaxis_type,              ONLY: ZA_REFERENCE
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                     ONLY: GRID_UNSTRUCTURED, DATATYPE_FLT32,      &
   &                                    DATATYPE_PACK16, DATATYPE_PACK24
  USE mo_var_metadata,            ONLY: create_vert_interp_metadata,            &
   &                                    create_hor_interp_metadata,             &
   &                                    vintp_types
  USE mo_run_config,              ONLY: iqg
  USE mo_sppt_config,             ONLY: sppt_config


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

  PUBLIC :: sppt

  PUBLIC :: construct_sppt_state
  PUBLIC :: destruct_sppt_state

  ! state variable
  TYPE(t_sppt),          ALLOCATABLE :: sppt(:)   ! n_dom

  ! variable list
  TYPE (t_var_list_ptr), ALLOCATABLE :: sppt_list(:)  ! n_dom


  CONTAINS


  !---------------------------------------------------------------------
  ! Constructor for sppt state
  !---------------------------------------------------------------------
  SUBROUTINE construct_sppt_state (p_patch)

    TYPE(t_patch),       INTENT(IN) :: p_patch(:)

    ! local variables
    INTEGER :: jg
    INTEGER :: ist
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname
    CHARACTER(*), PARAMETER :: routine = 'mo_sppt_state:construct_sppt_state'

    ! Allocate pointer arrays sppt, as well as the corresponding list arrays.
    ALLOCATE(sppt(n_dom), sppt_list(n_dom),STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'allocation of sppt array and list failed')
    ENDIF

    !$ACC ENTER DATA COPYIN(sppt)

    DO jg = 1, n_dom

      WRITE(listname,'(a,i2.2)') 'sppt_of_domain_',jg

      CALL new_sppt_list(p_patch(jg), listname, sppt_list(jg), sppt(jg))

    ENDDO
    CALL message(routine, 'construction of sppt state finished')

  END SUBROUTINE construct_sppt_state


  !---------------------------------------------------------------------
  ! Destructor for sppt state
  !---------------------------------------------------------------------

  SUBROUTINE destruct_sppt_state ()

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(*), PARAMETER :: routine = 'mo_sppt_state:destruct_sppt_state'

    !--------------------------------------------------------------

    ! delete sppt varlist
    DO jg = 1, n_dom
      CALL vlr_del(sppt_list(jg))
    ENDDO

    !$ACC WAIT(1)
    DO jg = 1, n_dom
      !$ACC EXIT DATA DELETE(sppt_config(jg)%taper)
      DEALLOCATE(sppt_config(jg)%taper, STAT=ist)
      IF(ist/=SUCCESS)THEN
        CALL finish (TRIM(routine), 'deallocation of sppt_config(:)%taper failed')
      ENDIF
    ENDDO

    !$ACC EXIT DATA DELETE(sppt, sppt_config)

    DEALLOCATE(sppt, sppt_list, STAT=ist)

    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'deallocation of sppt array and list failed')
    ENDIF

    CALL message(routine, 'destruction of sppt state finished')

  END SUBROUTINE destruct_sppt_state


  !---------------------------------------------------------------------
  ! Constructor for sppt state
  !---------------------------------------------------------------------

  SUBROUTINE new_sppt_list (p_patch, listname, sppt_list, sppt)

    TYPE(t_patch)         , INTENT(IN   ) :: p_patch
    CHARACTER(len=*)      , INTENT(IN   ) :: listname
    TYPE(t_var_list_ptr)  , INTENT(INOUT) :: sppt_list
    TYPE(t_sppt)          , INTENT(INOUT) :: sppt

    ! local variables
    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d_c(3)
    INTEGER :: shape2d_c(2)

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt
    INTEGER :: DATATYPE_PACK_VAR  !< variable "entropy" for some thermodynamic fields

    INTEGER :: nlev
    INTEGER :: nblks_c

    CHARACTER(*), PARAMETER :: routine = 'mo_sppt_state:new_sppt_list'

    !---------------------------------------------------------------------

    nlev    = p_patch%nlev
    nblks_c = p_patch%nblks_c

    ibits        = DATATYPE_PACK16   ! "entropy" of horizontal slice
    datatype_flt = DATATYPE_FLT32

    IF (gribout_config(p_patch%id)%lgribout_24bit) THEN  ! analysis
      ! higher accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK24
    ELSE
      ! standard accuracy for atmospheric thermodynamic fields
      DATATYPE_PACK_VAR = DATATYPE_PACK16
    ENDIF

    ! predefined array shapes
    shape3d_c     = (/nproma, nlev,   nblks_c /)
    shape2d_c     = (/nproma,         nblks_c /)

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------

    NULLIFY(sppt%temp_now,      &
     &      sppt%qv_now,        &
     &      sppt%qi_now,        &
     &      sppt%qr_now,        &
     &      sppt%qs_now,        &
     &      sppt%qc_now,        &
     &      sppt%qg_now,        &
     &      sppt%rn_3d,         &
     &      sppt%rn_2d_now,     &
     &      sppt%rn_2d_new,     &
     &      sppt%ddt_temp_fast, &
     &      sppt%ddt_u_fast,    &
     &      sppt%ddt_v_fast,    &
     &      sppt%ddt_qv_fast,   &
     &      sppt%ddt_qi_fast,   &
     &      sppt%ddt_qr_fast,   &
     &      sppt%ddt_qs_fast,   &
     &      sppt%ddt_qc_fast,   &
     &      sppt%ddt_qg_fast,   &
     &      sppt%ddt_qv,        &        ! tendencies fast and slow physics combined
     &      sppt%ddt_qi,        &
     &      sppt%ddt_qr,        &
     &      sppt%ddt_qs,        &
     &      sppt%ddt_qc,        &
     &      sppt%ddt_qg         )

    ! Register a field list and apply default settings
    CALL vlr_add(sppt_list, TRIM(listname), patch_id=p_patch%id, lrestart=.FALSE.)

    ! a) Fields to save current values

    !temp_now         sppt%temp_mow(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('temp_now', 'K', 'current value of air temperature sppt-only)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'temp_now', sppt%temp_now,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    

    ! qv_now         sppt%qv_now(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qv_now', 'kg kg-1', 'current value of tracer water vapour - sppt)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'qv_now', sppt%qv_now,                                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,                       &
                & ldims=shape3d_c, lrestart=.FALSE.,                                               &
                & hor_interp=create_hor_interp_metadata(                                           &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                                 &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                                 &
                & vert_interp=create_vert_interp_metadata(                                         &
                &             vert_intp_type=vintp_types("P","Z","I"),                             &
                &             vert_intp_method=VINTP_METHOD_LIN ),                                 &
                & lopenacc=.TRUE. )
      


    ! qi_now         sppt%qi_now(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qi_now', 'kg kg-1', 'current value of tracer ice - sppt)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'qi_now', sppt%qi_now,                                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,                       &
                & ldims=shape3d_c, lrestart=.FALSE.,                                               &
                & hor_interp=create_hor_interp_metadata(                                           &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                                 &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                                 &
                & vert_interp=create_vert_interp_metadata(                                         &
                &             vert_intp_type=vintp_types("P","Z","I"),                             &
                &             vert_intp_method=VINTP_METHOD_LIN ),                                 &
                & lopenacc=.TRUE. )
    


    ! qr_now         sppt%qr_now(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qr_now', 'kg kg-1', 'current value of tracer rain - sppt)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'qr_now', sppt%qr_now,                                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,                       &
                & ldims=shape3d_c, lrestart=.FALSE.,                                               &
                & hor_interp=create_hor_interp_metadata(                                           &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                                 &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                                 &
                & vert_interp=create_vert_interp_metadata(                                         &
                &             vert_intp_type=vintp_types("P","Z","I"),                             &
                &             vert_intp_method=VINTP_METHOD_LIN ),                                 &
                & lopenacc=.TRUE. )
    


    ! qs_now         sppt%qs_now(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qs_now', 'kg kg-1', 'current value of tracer snow - sppt)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'qs_now', sppt%qs_now,                                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,                       &
                & ldims=shape3d_c, lrestart=.FALSE.,                                               &
                & hor_interp=create_hor_interp_metadata(                                           &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                                 &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                                 &
                & vert_interp=create_vert_interp_metadata(                                         &
                &             vert_intp_type=vintp_types("P","Z","I"),                             &
                &             vert_intp_method=VINTP_METHOD_LIN ),                                 &
                & lopenacc=.TRUE. )
    


    ! qc_now         sppt%qc_now(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('qc_now', 'kg kg-1', 'current value of tracer cloud ice - sppt)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'qc_now', sppt%qc_now,                                                &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,                       &
                & ldims=shape3d_c, lrestart=.FALSE.,                                               &
                & hor_interp=create_hor_interp_metadata(                                           &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                                 &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                                 &
                & vert_interp=create_vert_interp_metadata(                                         &
                &             vert_intp_type=vintp_types("P","Z","I"),                             &
                &             vert_intp_method=VINTP_METHOD_LIN ),                                 &
                & lopenacc=.TRUE. )
    

    IF ( iqg /= 0 ) THEN
      ! qg_now         sppt%qg_now(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('qg_now', 'kg kg-1', 'current value of water tracer graupel - sppt)', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( sppt_list, 'qg_now', sppt%qg_now,                                                &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,                       &
                  & ldims=shape3d_c, lrestart=.FALSE.,                                               &
                  & hor_interp=create_hor_interp_metadata(                                           &
                  &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                                 &
                  &            fallback_type=HINTP_TYPE_LONLAT_RBF),                                 &
                  & vert_interp=create_vert_interp_metadata(                                         &
                  &             vert_intp_type=vintp_types("P","Z","I"),                             &
                  &             vert_intp_method=VINTP_METHOD_LIN ),                                 &
                  & lopenacc=.TRUE. )
      
    ENDIF

    ! b) fields for random number generation

    ! rn_3d         sppt%rn_3d(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('rn_3d', '-', 'field holding random numbers (sppt-only)', datatype_flt)
    grib2_desc = grib2_var(0, 19, 235, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'rn_3d', sppt%rn_3d,                                    &
               & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
               & ldims=shape3d_c, lrestart=.FALSE., loutput=.TRUE.,                  &
               & hor_interp=create_hor_interp_metadata(                              &
               &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
               &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
               & vert_interp=create_vert_interp_metadata(                            &
               &             vert_intp_type=vintp_types("P","Z","I"),                &
               &             vert_intp_method=VINTP_METHOD_LIN ),                    &
               & lopenacc=.TRUE. )
   


    ! rn_2d_now     sppt%rn_2d_now(nproma,nblks_c)
    cf_desc    = t_cf_var('rn_2d_now', '-', 'utility field for gaussian random number for sppt ', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'rn_2d_now', sppt%rn_2d_now,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
         & ldims=shape2d_c, lrestart=.FALSE.,                                  &
         & lopenacc=.TRUE. )



    ! rn_2d_new     sppt%rn_2d_new(nproma,nblks_c)
    cf_desc    = t_cf_var('rn_2d_new', '-', 'utility field for gaussian random number for sppt ', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'rn_2d_new', sppt%rn_2d_new,                       &
         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,           &
         & ldims=shape2d_c, lrestart=.FALSE.,                                   &
         & lopenacc=.TRUE. )



    ! c) additional fields for SPPT

    ! ddt_temp_fast         sppt%ddt_temp_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_temp_fast', 'K s-1', 'fast physics tendencies for air temperature)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_temp_fast', sppt%ddt_temp_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    

    ! ddt_u_fast         sppt%ddt_u_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_u_fast', 'm s-2', 'fast physics tendencies for u component)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_u_fast', sppt%ddt_u_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_v_fast         sppt%ddt_v_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_v_fast', 'm s-2', 'fast physics tendencies for v component)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_v_fast', sppt%ddt_v_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qv_fast         sppt%ddt_qv_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qv_fast', 'kg kg-1 s-1', 'fast physics tendencies for qv tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qv_fast', sppt%ddt_qv_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qi_fast         sppt%ddt_qi_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qi_fast', 'kg kg-1 s-1', 'fast physics tendencies for qi tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qi_fast', sppt%ddt_qi_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qr_fast         sppt%ddt_qr_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qr_fast', 'kg kg-1 s-1', 'fast physics tendencies for qr tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qr_fast', sppt%ddt_qr_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qs_fast         sppt%ddt_qs_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qs_fast', 'kg kg-1 s-1', 'fast physics tendencies for qs tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qs_fast', sppt%ddt_qs_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qc_fast         sppt%ddt_qc_fast(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qc_fast', 'kg kg-1 s-1', 'fast physics tendencies for qc tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qc_fast', sppt%ddt_qc_fast,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    

    IF ( iqg /= 0 ) THEN
      ! ddt_qg_fast         sppt%ddt_qg_fast(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('ddt_qg_fast', 'kg kg-1 s-1', 'fast physics tendencies for qg tracer)', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( sppt_list, 'ddt_qg_fast', sppt%ddt_qg_fast,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                  & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                  & hor_interp=create_hor_interp_metadata(                              &
                  &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                  & vert_interp=create_vert_interp_metadata(                            &
                  &             vert_intp_type=vintp_types("P","Z","I"),                &
                  &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                  & lopenacc=.TRUE. )
      
    ENDIF

    ! ddt_qv         sppt%ddt_qv(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qv', 'kg kg-1 s-1', 'total physics tendencies for qv tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qv', sppt%ddt_qv,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qi         sppt%ddt_qi(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qi', 'kg kg-1 s-1', 'total physics tendencies for qi tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qi', sppt%ddt_qi,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qr         sppt%ddt_qr(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qr', 'kg kg-1 s-1', 'total physics tendencies for qr tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qr', sppt%ddt_qr,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qs         sppt%ddt_qs(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qs', 'kg kg-1 s-1', 'total physics tendencies for qs tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qs', sppt%ddt_qs,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    


    ! ddt_qc         sppt%ddt_qc(nproma,nlev,nblks_c)
    cf_desc    = t_cf_var('ddt_qc', 'kg kg-1 s-1', 'total physics tendencies for qc tracer)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( sppt_list, 'ddt_qc', sppt%ddt_qc,                               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                & hor_interp=create_hor_interp_metadata(                              &
                &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                & vert_interp=create_vert_interp_metadata(                            &
                &             vert_intp_type=vintp_types("P","Z","I"),                &
                &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                & lopenacc=.TRUE. )
    

    IF ( iqg /= 0 ) THEN
      ! ddt_qg         sppt%ddt_qg(nproma,nlev,nblks_c)
      cf_desc    = t_cf_var('ddt_qg', 'kg kg-1 s-1', 'total physics tendencies for qg tracer)', datatype_flt)
      grib2_desc = grib2_var(255, 255, 255, DATATYPE_PACK_VAR, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( sppt_list, 'ddt_qg', sppt%ddt_qg,                               &
                  & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,          &
                  & ldims=shape3d_c, lrestart=.FALSE.,                                  &
                  & hor_interp=create_hor_interp_metadata(                              &
                  &            hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                    &
                  &            fallback_type=HINTP_TYPE_LONLAT_RBF),                    &
                  & vert_interp=create_vert_interp_metadata(                            &
                  &             vert_intp_type=vintp_types("P","Z","I"),                &
                  &             vert_intp_method=VINTP_METHOD_LIN ),                    &
                  & lopenacc=.TRUE. )
      
    ENDIF

  END SUBROUTINE new_sppt_list

END MODULE mo_sppt_state
