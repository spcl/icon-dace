! Construction of a data object which is used to store mass fluxes
! and velocities that are used to drive the tracer transport schemes.
! In contrast to the mass fluxes in the dynamical core, these mass fluxes
! are averaged over the dynamics substeps.
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

MODULE mo_prepadv_state

  USE mo_impl_constants,          ONLY: SUCCESS, MAX_CHAR_LENGTH, vname_len
  USE mo_exception,               ONLY: message, finish
  USE mo_master_control,          ONLY: get_my_process_name
  USE mo_model_domain,            ONLY: t_patch
  USE mo_prepadv_types,           ONLY: t_prepare_adv, t_step_adv
  USE mo_parallel_config,         ONLY: nproma
  USE mo_grid_config,             ONLY: n_dom
  USE mo_run_config,              ONLY: ntracer
  USE mo_var_list_register,       ONLY: vlr_add, vlr_del
  USE mo_var_list,                ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_zaxis_type,              ONLY: ZA_REFERENCE, ZA_REFERENCE_HALF, ZA_SURFACE
  USE mo_cdi,                     ONLY: DATATYPE_PACK16, GRID_UNSTRUCTURED, TSTEP_INSTANT, &
    &                                   DATATYPE_FLT32
  USE mo_cdi_constants,           ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, GRID_CELL, GRID_EDGE
  USE mo_cf_convention,           ONLY: t_cf_var
  USE mo_grib2,                   ONLY: t_grib2_var, grib2_var

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: prep_adv
  PUBLIC :: prep_adv_list
  PUBLIC :: jstep_adv

  PUBLIC :: construct_prepadv_state
  PUBLIC :: destruct_prepadv_state


  ! state variable
  TYPE(t_prepare_adv), ALLOCATABLE :: prep_adv(:)  ! n_dom
  TYPE(t_step_adv),    ALLOCATABLE :: jstep_adv(:) ! n_dom

  ! variable list
  TYPE (t_var_list_ptr), ALLOCATABLE :: prep_adv_list(:)  ! n_dom

CONTAINS

  !>
  !! Constructor for prepadv state
  SUBROUTINE construct_prepadv_state (p_patch)

    TYPE(t_patch), INTENT(in) :: p_patch(:)

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(len=MAX_CHAR_LENGTH) :: listname

    CHARACTER(*), PARAMETER :: routine = 'mo_prepadv_state:construct_prepadv_state'


    ! Allocate pointer arrays prep_adv, as well as the corresponding list arrays.
    !   
    ALLOCATE(prep_adv(n_dom), prep_adv_list(n_dom),STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'allocation of prep_adv array and list failed')
    ENDIF

    !$ACC ENTER DATA COPYIN(prep_adv)

    DO jg = 1, n_dom
      WRITE(listname,'(a,i2.2)') 'prepadv_of_domain_',jg
      CALL new_prep_adv_list( p_patch(jg), listname, prep_adv_list(jg), prep_adv(jg))
    ENDDO

    CALL message(routine, 'construction of prep_adv state finished')

  END SUBROUTINE construct_prepadv_state



  !>
  !! Destructor for prepadv state
  !!
  SUBROUTINE destruct_prepadv_state ()

    ! local variables
    INTEGER :: jg
    INTEGER :: ist                             !< error status
    CHARACTER(*), PARAMETER :: routine = 'mo_prepadv_state:destruct_prepadv_state'

    !--------------------------------------------------------------

    ! delete prep_adv varlist
    DO jg = 1, n_dom
      CALL vlr_del(prep_adv_list(jg))
    ENDDO

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(prep_adv)

    DEALLOCATE(prep_adv, prep_adv_list, STAT=ist)
    IF(ist/=SUCCESS)THEN
      CALL finish (TRIM(routine), 'deallocation of prep_adv array and list failed')
    ENDIF

    CALL message(routine, 'destruction of prep_adv state finished')

  END SUBROUTINE destruct_prepadv_state



  !>
  !! Constructor for prepadv state
  !!
  SUBROUTINE new_prep_adv_list (p_patch, listname, prep_adv_list, prep_adv)

    TYPE(t_patch)         , INTENT(IN   ) :: p_patch
    CHARACTER(len=*)      , INTENT(IN   ) :: listname
    TYPE(t_var_list_ptr)  , INTENT(INOUT) :: prep_adv_list
    TYPE(t_prepare_adv)   , INTENT(INOUT) :: prep_adv

    ! local variables
    INTEGER :: nblks_c, &    !< number of cell blocks to allocate
               nblks_e       !< number of edge blocks to allocate

    INTEGER :: nlev, nlevp1

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: shape3d_e(3), shape3d_chalf(3), shape3d_tracer(3), &
      &        shape2d_c(2)

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt

    CHARACTER(*), PARAMETER :: routine = 'mo_prepadv_state:new_prep_adv_list'

    CHARACTER(LEN=2) :: ctrnam
    CHARACTER(LEN=vname_len+LEN(ctrnam)) :: tracer_name
    INTEGER :: tlen
    INTEGER :: jt

    !--------------------------------------------------------------

    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e


    ! number of vertical levels
    nlev   = p_patch%nlev
    nlevp1 = p_patch%nlevp1

    ibits        = DATATYPE_PACK16   ! "entropy" of horizontal slice
    datatype_flt = DATATYPE_FLT32 

    shape3d_e     = (/nproma, nlev          , nblks_e /)
    shape3d_chalf = (/nproma, nlevp1        , nblks_c /)
    shape3d_tracer= (/nproma, MAX(1,ntracer), nblks_c /)
    shape2d_c     = (/nproma,                 nblks_c /)

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(prep_adv%mass_flx_me, &
      &     prep_adv%mass_flx_ic, &
      &     prep_adv%vol_flx_ic,  &
      &     prep_adv%vn_traj,     &
      &     prep_adv%q_int,       &
      &     prep_adv%q_ubc        )

    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(prep_adv_list, TRIM(listname), patch_id=p_patch%id, &
      &          lrestart=.FALSE., model_type=get_my_process_name())


    ! mass_flx_me      prep_adv%mass_flx_me(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('mass_flx_me', 'kg m-2 s-1', &
      &                   'horizontal mass flux (averaged over dynamics substeps)', &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( prep_adv_list, 'mass_flx_me', prep_adv%mass_flx_me,              &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_e, loutput=.FALSE.,                                &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%mass_flx_me)


    ! mass_flx_ic      prep_adv%mass_flx_ic(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('mass_flx_ic', 'kg m-2 s-1', &
      &                   'vertical mass flux (averaged over dynamics substeps)',  &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prep_adv_list, 'mass_flx_ic', prep_adv%mass_flx_ic,              &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, loutput=.FALSE.,                            &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%mass_flx_ic)

    ! vol_flx_ic      prep_adv%vol_flx_ic(nproma,nlevp1,nblks_c)
    cf_desc    = t_cf_var('vol_flx_ic', 'm s-1', &
      &                   'vertical volume flux (averaged over dynamics substeps)',  &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prep_adv_list, 'vol_flx_ic', prep_adv%vol_flx_ic,               &
                & GRID_UNSTRUCTURED_CELL, ZA_REFERENCE_HALF, cf_desc, grib2_desc,  &
                & ldims=shape3d_chalf, loutput=.FALSE.,                            &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%vol_flx_ic)


    ! vn_traj          prep_adv%vn_traj(nproma,nlev,nblks_e)
    cf_desc    = t_cf_var('vn_traj', 'm s-1', &
      &                   'velocity normal to edge (averaged over dynamics substeps)',  &
      &                   datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var( prep_adv_list, 'vn_traj', prep_adv%vn_traj,                      &
                & GRID_UNSTRUCTURED_EDGE, ZA_REFERENCE, cf_desc, grib2_desc,       &
                & ldims=shape3d_e, loutput=.FALSE.,                                &
                & isteptype=TSTEP_INSTANT, lopenacc=.TRUE. )
    __acc_attach(prep_adv%vn_traj)


    ! q_int        prep_adv%q_int(nproma,ntracer,nblks_c)
    !
    cf_desc    = t_cf_var('q_int', 'kg kg-1',                        &
      &                   'q at parent interface level', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prep_adv_list, 'q_int', prep_adv%q_int,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_tracer,                                       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       &
                & lopenacc = .TRUE. )
    __acc_attach(prep_adv%q_int)

    ALLOCATE(prep_adv%q_int_ptr(ntracer))
    DO jt =1,ntracer
      WRITE(ctrnam,'(I2)') jt
      tracer_name = 'q_int'//ctrnam(1+MERGE(1,0,jt<=9):)
      tlen = LEN_TRIM(tracer_name)
      CALL add_ref( prep_adv_list, 'q_int',                                       &
                  & tracer_name(1:tlen), prep_adv%q_int_ptr(jt)%p_2d,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                  & t_cf_var(tracer_name(1:tlen), 'kg kg-1','', datatype_flt),    &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                  & ref_idx=jt, opt_var_ref_pos=2,                                &
                  & ldims=shape2d_c, lrestart=.FALSE. )
    ENDDO


    ! q_ubc        prep_adv%q_ubc(nproma,ntracer,nblks_c)
    !
    cf_desc    = t_cf_var('q_ubc', 'kg kg-1',                      &
      &                   'q at child upper boundary', datatype_flt)
    grib2_desc = grib2_var( 255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( prep_adv_list, 'q_ubc', prep_adv%q_ubc,                     &
                & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
                & ldims=shape3d_tracer,                                       &
                & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,       &
                & lopenacc = .TRUE. )
    __acc_attach(prep_adv%q_ubc)

    ALLOCATE(prep_adv%q_ubc_ptr(ntracer))
    DO jt =1,ntracer
      WRITE(ctrnam,'(I2)') jt
      tracer_name = 'q_ubc'//ctrnam(1+MERGE(1,0,jt<=9):)
      tlen = LEN_TRIM(tracer_name)
      CALL add_ref( prep_adv_list, 'q_ubc',                                       &
                  & tracer_name(1:tlen), prep_adv%q_ubc_ptr(jt)%p_2d,             &
                  & GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                           &
                  & t_cf_var(tracer_name(1:tlen), 'kg kg-1','', datatype_flt),    &
                  & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL),&
                  & ref_idx=jt, opt_var_ref_pos=2,                                &
                  & ldims=shape2d_c, lrestart=.FALSE. )
    ENDDO

  END SUBROUTINE new_prep_adv_list


END MODULE mo_prepadv_state

