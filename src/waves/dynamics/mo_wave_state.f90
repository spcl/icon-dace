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

! Contains the variables to set up the wave  model.

MODULE mo_wave_state

  USE mo_master_control,            ONLY: get_my_process_name
  USE mo_exception,                 ONLY: message, finish
  USE mo_parallel_config,           ONLY: nproma
  USE mo_model_domain,              ONLY: t_patch
  USE mo_grid_config,               ONLY: n_dom, l_limited_area, ifeedback_type
  USE mo_impl_constants,            ONLY: success, max_char_length, VNAME_LEN, TLEV_NNOW_RCF, &
    &                                     HINTP_TYPE_LONLAT_NNB, HINTP_TYPE_LONLAT_BCTR
  USE mo_var_list,                  ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_var_list_register,         ONLY: vlr_add, vlr_del
  USE mo_var_groups,                ONLY: groups
  USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL, &
       &                                  GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_cdi,                       ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED, &
       &                                  DATATYPE_PACK16, DATATYPE_INT
  USE mo_zaxis_type,                ONLY: ZA_SURFACE, ZA_REFERENCE
  USE mo_cf_convention,             ONLY: t_cf_var
  USE mo_grib2,                     ONLY: t_grib2_var, grib2_var
  USE mo_io_config,                 ONLY: lnetcdf_flt64_output
  USE mo_run_config,                ONLY: ntracer
  USE mo_var_metadata,              ONLY: get_timelevel_string, create_hor_interp_metadata
  USE mo_tracer_metadata,           ONLY: create_tracer_metadata

  USE mo_wave_types,                ONLY: t_wave_prog, t_wave_source, t_wave_diag, &
       &                                  t_wave_state, t_wave_state_lists
  USE mo_wave_config,               ONLY: t_wave_config, wave_config
  USE mo_energy_propagation_config, ONLY: t_energy_propagation_config, energy_propagation_config


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_state'

  PUBLIC :: construct_wave_state    ! Constructor for the wave state
  PUBLIC :: destruct_wave_state     ! Destructor

  PUBLIC :: p_wave_state            ! state vector of wave variables
  PUBLIC :: p_wave_state_lists      ! lists for state vector of wave variables

  TYPE(t_wave_state),       TARGET, ALLOCATABLE :: p_wave_state(:)
  TYPE(t_wave_state_lists), TARGET, ALLOCATABLE :: p_wave_state_lists(:)

CONTAINS

  SUBROUTINE construct_wave_state(p_patch, n_timelevels)

    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    INTEGER,       INTENT(IN) :: n_timelevels

    CHARACTER(len=max_char_length) :: listname
    CHARACTER(len=*), PARAMETER :: routine = modname//'::construct_wave_state'

    INTEGER :: ntl,  &! local number of timelevels
         ist,        &! status
         jg,         &! grid level counter
         jt           ! time level counter

    CALL message (routine, 'Construction of wave state started')

    ALLOCATE (p_wave_state(n_dom),p_wave_state_lists(n_dom), stat=ist)
    IF (ist /= success) THEN
       CALL finish(TRIM(routine),'allocation for wave state failed')
    END IF

    DO jg = 1, n_dom

       ntl = n_timelevels

       ! As grid nesting is not called at every dynamics time step, an extra time
       ! level is needed for full-field interpolation and boundary-tendency calculation
       IF (n_dom > 1) THEN
          ntl = ntl + 1
       END IF

       IF (ifeedback_type == 1 .AND. jg > 1 .OR. l_limited_area .AND. jg == 1) ntl = ntl + 1

       ALLOCATE(p_wave_state(jg)%prog(1:ntl), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic state array failed')

       ! create state list
       ALLOCATE(p_wave_state_lists(jg)%prog_list(1:ntl), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic state list array failed')

       !
       ! Build lists for every timelevel
       !
       DO jt = 1, ntl

          WRITE(listname,'(a,i2.2,a,i2.2)') 'wave_state_prog_of_domain_',jg, &
               &                            '_and_timelev_',jt

          ! Build prog state list
          ! includes memory allocation
          CALL new_wave_state_prog_list(p_patch(jg), p_wave_state(jg)%prog(jt), &
               & p_wave_state_lists(jg)%prog_list(jt), &
               & listname, jt)

       END DO

       ! Build source state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_source_of_domain_',jg
       CALL new_wave_state_source_list(&
            p_patch(jg), &
            p_wave_state(jg)%source, &
            p_wave_state_lists(jg)%source_list, &
            listname)

       ! Build diag state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_diag_of_domain_',jg
       CALL new_wave_state_diag_list(&
            p_patch(jg), &
            p_wave_state(jg)%diag, &
            p_wave_state_lists(jg)%diag_list, &
            listname)

    END DO

    CALL message (routine, 'wave state construction completed')

  END SUBROUTINE construct_wave_state


  SUBROUTINE new_wave_state_prog_list(p_patch, p_prog, p_prog_list, listname, timelev)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_prog),     INTENT(INOUT) :: p_prog
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_prog_list !< current prognostic state list
    CHARACTER(len=*),      INTENT(IN)    :: listname
    INTEGER,               INTENT(IN)    :: timelev

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c    !< number of cell blocks to allocate

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output

    CHARACTER(len=4)         :: suffix
    CHARACTER(len=VNAME_LEN) :: freq_ind_str, dir_ind_str
    CHARACTER(LEN=VNAME_LEN) :: tracer_container_name
    CHARACTER(len=VNAME_LEN) :: tracer_name

    TYPE(t_energy_propagation_config), POINTER :: enprop_conf
    TYPE(t_wave_config),               POINTER :: wc


    INTEGER :: shape3d_c(3), shape4d_c(4)
    INTEGER :: jt, nlev

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! number of vertical levels
    nlev    = p_patch%nlev

    ! pointer to energy_propagation_config(jg) to save some paperwork
    enprop_conf => energy_propagation_config(p_patch%id)
    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    shape4d_c    = (/nproma, nlev, nblks_c, ntracer/)
    shape3d_c    = (/nproma, nlev, nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF (lnetcdf_flt64_output) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    END IF

    ! Suffix (mandatory for time level dependent variables)
    suffix = get_timelevel_string(timelev)

    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_prog_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &          model_type=get_my_process_name())


    tracer_container_name = 'tracer'//suffix
    cf_desc    = t_cf_var('tracer', '', 'spectral bin of wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( p_prog_list, tracer_container_name, p_prog%tracer, &
         & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE, cf_desc, grib2_desc,&
         & ldims=shape4d_c ,                                         &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

    ALLOCATE(p_prog%tracer_ptr(ntracer))

    DO jt = 1, ntracer
       write(freq_ind_str,'(I0.3)') wc%freq_ind(jt)
       write(dir_ind_str,'(I0.3)') wc%dir_ind(jt)

       tracer_name = TRIM(enprop_conf%tracer_names(jt))//'_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)//suffix

       CALL add_ref( p_prog_list, tracer_container_name,                          &
            & TRIM(tracer_name), p_prog%tracer_ptr(jt)%p_3d,                      &
            & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE,                               &
            & t_cf_var(TRIM(tracer_name), '-','spectral bin of wave energy',      &
            & datatype_flt),                                                      &
            & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),      &
            & ref_idx=jt,                                                         &
            & ldims=shape3d_c,                                                    &
            & loutput=.TRUE., lrestart=.TRUE.,                                   &
            & tlev_source=TLEV_NNOW_RCF,                                          &
            & tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,               &
            &                       name        = TRIM(tracer_name),              &
            &                       lfeedback   = .TRUE.,                         &
            &                       ihadv_tracer= 2,                              &
            &                       ivadv_tracer= 0 ),                            &
            & in_group=groups("wave_spectrum"))

    END DO

  END SUBROUTINE new_wave_state_prog_list



  SUBROUTINE new_wave_state_source_list(p_patch, p_source, p_source_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_source),   INTENT(INOUT) :: p_source
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_source_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_source_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c
    INTEGER :: jt
    INTEGER :: ist
    INTEGER :: shape3d_tr_c(3), shape2d_c(2)
    CHARACTER(len=VNAME_LEN) :: sl_name, fl_name, llws_name
    CHARACTER(len=3) :: freq_ind_str, dir_ind_str

    TYPE(t_wave_config), POINTER :: wc

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    nblks_c = p_patch%nblks_c

    shape3d_tr_c = (/nproma, nblks_c, ntracer/)
    shape2d_c    = (/nproma, nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(p_source_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &           model_type=get_my_process_name())


    ! fl          p_source%fl(nproma,nblks_c,ntracer)
    cf_desc    = t_cf_var('fl', '-', 'DIAG. MTRX OF FUNC. DERIVATIVE', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_source_list, 'fl', p_source%fl,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_tr_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! sl          p_source%sl(nproma,nblks_c,ntracer)
    cf_desc    = t_cf_var('sl', '-', 'TOTAL SOURCE FUNCTION', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_source_list, 'sl', p_source%sl,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_tr_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! llws        p_source%llws(nproma,nblks_c,ntracer)
    cf_desc    = t_cf_var('llws', '-', '1 where sinput is positive', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_source_list, 'llws', p_source%llws,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_tr_c,                                      &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ALLOCATE(p_source%sl_ptr(ntracer),           &
      &      p_source%fl_ptr(ntracer),           &
      &      p_source%llws_ptr(ntracer), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine,               &
          'allocation of sl_ptr, fl_ptr, llws_ptr failed')

    DO jt = 1, ntracer
      write(freq_ind_str,'(I0.3)') wc%freq_ind(jt)
      write(dir_ind_str,'(I0.3)') wc%dir_ind(jt)

      sl_name = 'sl_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)
      CALL add_ref(p_source_list, 'sl',                                     &
           & sl_name, p_source%sl_ptr(jt)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
           & t_cf_var(sl_name, '-',sl_name, datatype_flt),                  &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ref_idx=jt, ldims=shape2d_c, lrestart=.FALSE., loutput=.TRUE.)

      fl_name = 'fl_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)
      CALL add_ref(p_source_list, 'fl',                                     &
           & fl_name, p_source%fl_ptr(jt)%p_2d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
           & t_cf_var(fl_name, '-',fl_name, datatype_flt),                  &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ref_idx=jt, ldims=shape2d_c, lrestart=.FALSE., loutput=.TRUE.)

      llws_name = 'llws_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)
      CALL add_ref(p_source_list, 'llws',                                   &
           & llws_name, p_source%llws_ptr(jt)%p,                            &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
           & t_cf_var(llws_name, '-',llws_name, datatype_int),              &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ref_idx=jt, ldims=shape2d_c, lrestart=.TRUE., loutput=.TRUE.,  &
           & in_group=groups("wave_debug"))
    END DO

  END SUBROUTINE new_wave_state_source_list



  SUBROUTINE new_wave_state_diag_list(p_patch, p_diag, p_diag_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_diag),     INTENT(INOUT) :: p_diag
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_diag_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_diag_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c, nblks_e, nlev
    INTEGER :: nfreqs, ndirs
    INTEGER :: jg, jf, jt
    INTEGER :: ist
    INTEGER :: shape2d_c(2), shape2d_e(2)
    INTEGER :: shape3d_freq_c(3), shape3d_freq_e(3)
    INTEGER :: shape3d_tr_c(3), shape4d_tr_e(4)
    INTEGER :: shape1d_freq_p4(1), shape1d_dir_2(2)
    INTEGER :: shape4d_freq_p4_2_dir_18(4)

    CHARACTER(len=3) :: freq_ind_str, dir_ind_str
    CHARACTER(len=VNAME_LEN) :: out_name

    TYPE(t_wave_config),      POINTER :: wc

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    jg      = p_patch%id
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e
    nlev    = p_patch%nlev


    nfreqs =  wave_config(jg)%nfreqs
    ndirs = wave_config(jg)%ndirs

    shape1d_freq_p4   = (/nfreqs+4/)
    shape1d_dir_2     = (/ndirs, 2/)
    shape2d_c         = (/nproma, nblks_c/)
    shape2d_e         = (/nproma, nblks_e/)
    shape3d_freq_c    = (/nproma, nblks_c, nfreqs/)
    shape3d_freq_e    = (/nproma, nblks_e, nfreqs/)
    shape3d_tr_c      = (/nproma, nblks_c, ntracer/)
    shape4d_tr_e      = (/nproma, nlev, nblks_e, ntracer/)

    shape4d_freq_p4_2_dir_18 = (/18,nfreqs+4,2,ndirs/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(p_diag_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &         model_type=get_my_process_name())


    !wave group velocity
    cf_desc    = t_cf_var('gv_c', 'm s-1', 'group velocity at cells', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'gv_c_freq', p_diag%gv_c,             &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_c,                                    &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    cf_desc    = t_cf_var('gv_e', 'm s-1', 'group velocity at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(p_diag_list, 'gv_e_freq', p_diag%gv_e,             &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_e,                                    &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ALLOCATE(p_diag%gv_e_freq_ptr(nfreqs), p_diag%gv_c_freq_ptr(nfreqs), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine,                            &
          'allocation of gv_e_freq_ptr, gv_c_freq_ptr failed')

    DO jf = 1, nfreqs
      write(freq_ind_str,'(I0.3)') jf
      out_name = 'gv_c_'//TRIM(freq_ind_str)
      CALL add_ref(p_diag_list, 'gv_c_freq',                                &
           & TRIM(out_name), p_diag%gv_c_freq_ptr(jf)%p_2d,                 &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
           & t_cf_var(TRIM(out_name), 'm s-1','group velocity at cells',    &
           & datatype_flt),                                                 &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ref_idx=jf, ldims=shape2d_c, lrestart=.FALSE., loutput=.TRUE., &
           & in_group=groups("wave_phy_ext"))

      out_name = 'gv_e_'//TRIM(freq_ind_str)
      CALL add_ref(p_diag_list, 'gv_e_freq',                                &
           & TRIM(out_name), p_diag%gv_e_freq_ptr(jf)%p_2d,                 &
           & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE,                            &
           & t_cf_var(TRIM(out_name), 'm s-1','group velocity at edges',    &
           & datatype_flt),                                                 &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE), &
           & ref_idx=jf, ldims=shape2d_e, lrestart=.FALSE., loutput=.TRUE., &
           & in_group=groups("wave_phy_ext"))

    END DO

    cf_desc    = t_cf_var('normal_group_velocity', 'm s-1', 'group velocity normal to edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(p_diag_list, 'gvn_e', p_diag%gvn_e,                &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape4d_tr_e, lrestart=.FALSE., loutput=.TRUE.)

    cf_desc    = t_cf_var('tangential_group_velocity', 'm s-1', 'group velocity tangential to edge', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(p_diag_list, 'gvt_e', p_diag%gvt_e,                &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape4d_tr_e, lrestart=.FALSE., loutput=.TRUE.)

    !Wave physics group
    cf_desc    = t_cf_var('emean', 'm^2', 'total wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'emean', p_diag%emean,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('emeanws', 'm^2', 'wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'emeanws', p_diag%emeanws,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('femean', 'm^2', 'mean frequency wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'femean', p_diag%femean,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('hrms_frac', '-', 'square ratio (Hrms / Hmax)**2', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hrmc_frac', p_diag%hrms_frac,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape2d_c)

    cf_desc    = t_cf_var('wbr_frac', '-', 'fraction of breaking waves', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'wbr_frac', p_diag%wbr_frac,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('wave_num_c', '1/m', 'wave number at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'wave_num_c', p_diag%wave_num_c,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_c,                                    &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ALLOCATE(p_diag%wave_num_c_ptr(nfreqs))
    !
    DO jf = 1, nfreqs
      write(freq_ind_str,'(I0.3)') jf
      out_name = 'wave_num_c_'//TRIM(freq_ind_str)
      CALL add_ref(p_diag_list, 'wave_num_c',                               &
           & TRIM(out_name), p_diag%wave_num_c_ptr(jf)%p_2d,                &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
           & t_cf_var(TRIM(out_name), 'm-1','wave number at cell center',   &
           & datatype_flt),                                                 &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ref_idx=jf, ldims=shape2d_c, lrestart=.FALSE., loutput=.TRUE.)
    END DO


    cf_desc    = t_cf_var('wave_num_e', 'm-1', 'wave number at edge midpoint', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(p_diag_list, 'wave_num_e', p_diag%wave_num_e,      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_e,                                    &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ALLOCATE(p_diag%wave_num_e_ptr(nfreqs))
    !
    DO jf = 1, nfreqs
      write(freq_ind_str,'(I0.3)') jf
      out_name = 'wave_num_e_'//TRIM(freq_ind_str)
      CALL add_ref(p_diag_list, 'wave_num_e',                               &
           & TRIM(out_name), p_diag%wave_num_e_ptr(jf)%p_2d,                &
           & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE,                            &
           & t_cf_var(TRIM(out_name), 'm-1','wave number at edge midpoint', &
           & datatype_flt),                                                 &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE), &
           & ref_idx=jf, ldims=shape2d_e, lrestart=.FALSE., loutput=.TRUE.)
    END DO


    cf_desc    = t_cf_var('f1mean', 'm^2', 'mean frequency wave energy based on F-moment', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'f1mean', p_diag%f1mean,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('femeanws', 'm^2', 'mean frequency wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'femeanws', p_diag%femeanws,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('akmean', '', 'Mean wavenumber based on SQRT(1/K)-moment, wm1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'akmean', p_diag%akmean,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('xkmean', '', 'Mean wavenumber based on SQRT(K)-moment, wm2', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'xkmean', p_diag%xkmean,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('swell_mask', '-', 'swell mask', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'swell_mask', p_diag%swell_mask,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    cf_desc    = t_cf_var('swell_mask_tr', '-', 'swell mask for tracers', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'swell_mask_tr', p_diag%swell_mask_tr,&
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,    &
         & ldims=shape3d_tr_c)

    ALLOCATE(p_diag%swmask_ptr(ntracer), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine,               &
          'allocation of swmask_ptr failed')

    DO jt = 1, ntracer
      write(freq_ind_str,'(I0.3)') wc%freq_ind(jt)
      write(dir_ind_str,'(I0.3)') wc%dir_ind(jt)

      out_name = 'swmask_'//TRIM(freq_ind_str)//'_'//TRIM(dir_ind_str)
      CALL add_ref(p_diag_list, 'swell_mask_tr',                            &
           & out_name, p_diag%swmask_ptr(jt)%p,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                            &
           & t_cf_var(out_name, '-',out_name, datatype_int),                &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ref_idx=jt, ldims=shape2d_c, lrestart=.TRUE., loutput=.TRUE.,  &
           & in_group=groups("wave_debug"))
    END DO


    cf_desc    = t_cf_var('last_prog_freq_ind', '-', 'last frequency index of the prognostic range', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'last_prog_freq_ind', p_diag%last_prog_freq_ind, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('ALPHAJ', '-', 'JONSWAP alpha', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'ALPHAJ', p_diag%ALPHAJ,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.TRUE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('FP', 'Hz', 'JONSWAP peak frequency', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'FP', p_diag%FP,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.TRUE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('ET', '-', 'JONSWAP spectra', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'ET', p_diag%ET,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape3d_freq_c)

    cf_desc    = t_cf_var('FLMINFR', '-', 'THE MINIMUM VALUE IN SPECTRAL', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'FLMINFR', p_diag%FLMINFR,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape3d_freq_c)

    cf_desc    = t_cf_var('friction_velocity', 'm s-1', 'friction velocity', datatype_flt)
    grib2_desc = grib2_var(10, 0, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'ustar', p_diag%ustar,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('roughness length', 'm', 'Surface roughness length', datatype_flt)
    grib2_desc = grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'z0', p_diag%z0,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('wave_stress', '(m/s)**2', 'wave stress', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tauw', p_diag%tauw,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.TRUE., loutput=.TRUE.,                         &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('integrated_energy_flux', '-', 'integrated energy flux', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'phiaw', p_diag%phiaw,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('tauhf1', '-', 'high-fequency stress', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tauhf1', p_diag%tauhf1,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('phihf1', '-', 'high-frequency energy flux into ocean', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'phihf1', p_diag%phihf1,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('tauhf', '-', 'high-fequency stress', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tauhf', p_diag%tauhf,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('phihf', '-', 'high-frequency energy flux into ocean', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'phihf', p_diag%phihf,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('xlevtail', '-', 'tail level', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'xlevtail', p_diag%xlevtail,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! nonlinear transfer function coefficients for shallow water
    cf_desc    = t_cf_var('enh', '-', 'nonlinear transfer function coefficients', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'enh', p_diag%enh,                    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    ! for discrete approximation of nonlinear transfer
    cf_desc    = t_cf_var('IKP', '-', 'IKP', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'IKP', p_diag%IKP,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('IKP1', '-', 'IKP1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'IKP1', p_diag%IKP1,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('IKM', '-', 'IKM', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'IKM', p_diag%IKM,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('IKM1', '-', 'IKM1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'IKM1', p_diag%IKM1,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('K1W', '-', 'K1W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'K1W', p_diag%K1W,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('K2W', '-', 'K2W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'K2W', p_diag%K2W,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('K11W', '-', 'K11W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'K11W', p_diag%K11W,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('K21W', '-', 'K21W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'K21W', p_diag%K21W,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('non_lin_tr_ind', '-', 'non_lin_tr_ind', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'non_lin_tr_ind', p_diag%non_lin_tr_ind, &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc,    &
         & lrestart=.FALSE., loutput=.TRUE.,                           &
         & ldims=shape4d_freq_p4_2_dir_18)

    cf_desc    = t_cf_var('JA1', '-', 'JA1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'JA1', p_diag%JA1,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('JA2', '-', 'JA2', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'JA2', p_diag%JA2,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_dir_2)

    cf_desc    = t_cf_var('AF11', '-', 'AF11', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'AF11', p_diag%AF11,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('FKLAP', '-', 'FKLAP', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'FKLAP', p_diag%FKLAP,                &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('FKLAP1', '-', 'FKLAP1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'FKLAP1', p_diag%FKLAP1,              &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('FKLAM', '-', 'FKLAM', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'FKLAM', p_diag%FKLAM,                &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    cf_desc    = t_cf_var('FKLAM1', '-', 'FKLAM1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'FKLAM1', p_diag%FKLAM1,              &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape1d_freq_p4)

    ! wave group
    ! total
    cf_desc    = t_cf_var('hs', 'm', 'Total significant wave height', datatype_flt)
    grib2_desc = grib2_var(10, 0, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hs', p_diag%hs,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & hor_interp=create_hor_interp_metadata(                   &
         &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                 &
         &    fallback_type=HINTP_TYPE_LONLAT_NNB),                 &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('hs_dir', 'deg', 'Total mean wave direction', datatype_flt)
    grib2_desc = grib2_var(10, 0, 10, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hs_dir', p_diag%hs_dir,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('tpp', 's', 'Total wave peak period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 34, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tpp', p_diag%tpp,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('tmp', 's', 'Total wave mean period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tmp', p_diag%tmp,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('tm1', 's', 'Total m1 wave period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tm1', p_diag%tm1,                    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('tm2', 's', 'Total m2 wave period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tm2', p_diag%tm2,                    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('ds', 'deg', 'Total directional wave spread', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'ds', p_diag%ds,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! wind sea
    cf_desc    = t_cf_var('emean_sea', 'm^2', 'Wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'emean_sea', p_diag%emean_sea,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('femean_sea', 'm^2', 'Wind sea mean frequency wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'femean_sea', p_diag%femean_sea,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('f1mean_sea', 'm^2', 'Wind sea mean frequency', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'f1mean_sea', p_diag%f1mean_sea,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('hs_sea', 'm', 'Sea significant wave height', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hs_sea', p_diag%hs_sea,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('hs_sea_dir', 'deg', 'Sea mean wave direction', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hs_sea_dir', p_diag%hs_sea_dir,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('pp_sea', 's', 'Sea wave peak period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'pp_sea', p_diag%pp_sea,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('mp_sea', 's', 'Sea wave mean period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'mp_sea', p_diag%mp_sea,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('m1_sea', 's', 'Sea m1 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'm1_sea', p_diag%m1_sea,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('m2_sea', 's', 'Sea m2 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'm2_sea', p_diag%m2_sea,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('ds_sea', 'deg', 'Sea directional wave spread', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'ds_sea', p_diag%ds_sea,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

     ! swell
    cf_desc    = t_cf_var('emean_swell', 'm^2', 'Swell wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'emean_swell', p_diag%emean_swell,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('femean_swell', 'm^2', 'Swell mean frequency wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'femean_swell', p_diag%femean_swell,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('f1mean_swell', 'm^2', 'Swell mean frequency', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'f1mean_swell', p_diag%f1mean_swell,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('hs_swell', 'm', 'Swell significant wave height', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hs_swell', p_diag%hs_swell,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('hs_swell_dir', 'deg', 'Swell mean wave direction', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'hs_swell_dir', p_diag%hs_swell_dir,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('pp_swell', 's', 'Swell wave peak period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'pp_swell', p_diag%pp_swell,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('mp_swell', 's', 'Swell wave mean period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'mp_swell', p_diag%mp_swell,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('m1_swell', 's', 'Swell m1 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'm1_swell', p_diag%m1_swell,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('m2_swell', 's', 'Swell m2 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'm2_swell', p_diag%m2_swell,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('ds_swell', 'deg', 'Swell directional wave spread', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'ds_swell', p_diag%ds_swell,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('drag', '-', 'Drag coefficient', datatype_flt)
    grib2_desc = grib2_var(10, 0, 16, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'drag', p_diag%drag,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    cf_desc    = t_cf_var('tauwn', '-', 'Normalized wave stress', datatype_flt)
    grib2_desc = grib2_var(10, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'tauwn', p_diag%tauwn,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    cf_desc    = t_cf_var('beta', '-', 'Charnock parameter', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'beta', p_diag%beta,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    cf_desc    = t_cf_var('u_stokes', 'ms-1', 'U-component surface Stokes drift', datatype_flt)
    grib2_desc = grib2_var(10, 0, 21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'u_stokes', p_diag%u_stokes,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c)

    cf_desc    = t_cf_var('v_stokes', 'ms-1', 'V-component surface Stokes drift', datatype_flt)
    grib2_desc = grib2_var(10, 0, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_diag_list, 'v_stokes', p_diag%v_stokes,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c)


  END SUBROUTINE new_wave_state_diag_list

  !>
  !! Destruction of wave-specific variable lists and memory deallocation
  !!
  SUBROUTINE destruct_wave_state ()

    INTEGER :: ist
    INTEGER :: jg, jt
    CHARACTER(len=*), PARAMETER :: routine = modname//'::destruct_wave_state'

    DO jg = 1, n_dom
      ! delete prognostic state list elements
      DO jt = 1, SIZE(p_wave_state_lists(jg)%prog_list(:))
        CALL vlr_del(p_wave_state_lists(jg)%prog_list(jt))
      ENDDO

      ! delete source state list elements
      CALL vlr_del(p_wave_state_lists(jg)%source_list)

      ! delete diagnostics state list elements
      CALL vlr_del(p_wave_state_lists(jg)%diag_list)

      ! deallocate state lists and arrays
      DEALLOCATE(p_wave_state_lists(jg)%prog_list, stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'deallocation for prog_list array failed')
      END IF

      DEALLOCATE(p_wave_state(jg)%prog, stat=ist)
      IF (ist /= success) THEN
        CALL finish(TRIM(routine),'deallocation of prognostic state array failed')
      END IF
    ENDDO

    ! deallocate states
    DEALLOCATE(p_wave_state, p_wave_state_lists, stat=ist)
    IF (ist /= success) THEN
      CALL finish(TRIM(routine),'deallocation for p_wave_state failed')
    END IF

  END SUBROUTINE destruct_wave_state


END MODULE mo_wave_state
