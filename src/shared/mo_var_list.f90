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

#include "omp_definitions.inc"

MODULE mo_var_list

#if defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR)
#ifdef VARLIST_INITIZIALIZE_WITH_NAN
  USE, INTRINSIC :: ieee_arithmetic
#endif
#endif
  USE mo_kind,             ONLY: sp, dp, i8
  USE mo_cf_convention,    ONLY: t_cf_var
  USE mo_grib2,            ONLY: t_grib2_var, grib2_var
  USE mo_var_groups,       ONLY: var_groups_dyn, groups
  USE mo_var_metadata_types,ONLY: t_var_metadata, t_vert_interp_meta, &
    & t_union_vals, CLASS_TILE, t_hor_interp_meta, t_post_op_meta,    &
    & CLASS_TILE_LAND, t_var_metadata_dynamic
  USE mo_var_metadata,     ONLY: create_vert_interp_metadata,       &
    & create_hor_interp_metadata, get_var_timelevel, get_var_name,  &
    & set_var_metadata, set_var_metadata_dyn
  USE mo_tracer_metadata_types, ONLY: t_tracer_meta
  USE mo_var,              ONLY: t_var, t_var_ptr, level_type_ml
  USE mo_exception,        ONLY: message, finish, message_text
  USE mo_util_texthash,    ONLY: text_hash_c
  USE mo_util_string,      ONLY: tolower
  USE mo_impl_constants,   ONLY: REAL_T, SINGLE_T, BOOL_T, INT_T, &
    & vlname_len, vname_len, TIMELEVEL_SUFFIX
  USE mo_fortran_tools,    ONLY: init_contiguous_dp, init_contiguous_sp, &
    &                            init_contiguous_i4, init_contiguous_l
  USE mo_action_types,     ONLY: t_var_action

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: add_var, add_ref, find_list_element
  PUBLIC :: t_var_list_ptr
  PUBLIC :: get_tracer_info_dyn_by_idx, find_tracer_by_index

  TYPE :: t_var_list
    CHARACTER(len=256) :: filename = ''
    CHARACTER(len=vlname_len) :: vlname = ''
    INTEGER(i8) :: memory_used = 0_i8
    CHARACTER(len=8) :: post_suf = '', rest_suf = '', init_suf = '', &
      & model_type = ''
    INTEGER :: patch_id = -1, nvars = 0, vlevel_type = level_type_ml, &
      & output_type = -1, restart_type = -1, compression_type = -1
    LOGICAL :: loutput = .TRUE., lrestart = .FALSE., linitial = .FALSE., &
      & restart_opened = .FALSE., output_opened = .FALSE., lmiss = .FALSE., &
      & lmask_boundary = .TRUE. , first = .FALSE.
    INTEGER, ALLOCATABLE :: tl(:), hgrid(:), key(:), key_notl(:)
    LOGICAL, ALLOCATABLE :: lout(:)
    TYPE(t_var_ptr), ALLOCATABLE :: vl(:)
  END TYPE t_var_list

  TYPE :: t_var_list_ptr
    TYPE(t_var_list), POINTER :: p => NULL()
  CONTAINS
    PROCEDURE :: register => register_list_element
    PROCEDURE :: delete => delete_list
    PROCEDURE :: print => print_var_list
  END TYPE t_var_list_ptr

 INTERFACE add_var  ! create a new list entry
    MODULE PROCEDURE add_var_list_element_5d
    MODULE PROCEDURE add_var_list_element_r4d
    MODULE PROCEDURE add_var_list_element_r3d
    MODULE PROCEDURE add_var_list_element_r2d
    MODULE PROCEDURE add_var_list_element_r1d
    MODULE PROCEDURE add_var_list_element_s4d
    MODULE PROCEDURE add_var_list_element_s3d
    MODULE PROCEDURE add_var_list_element_s2d
    MODULE PROCEDURE add_var_list_element_s1d
    MODULE PROCEDURE add_var_list_element_i4d
    MODULE PROCEDURE add_var_list_element_i3d
    MODULE PROCEDURE add_var_list_element_i2d
    MODULE PROCEDURE add_var_list_element_i1d
    MODULE PROCEDURE add_var_list_element_l4d
    MODULE PROCEDURE add_var_list_element_l3d
    MODULE PROCEDURE add_var_list_element_l2d
    MODULE PROCEDURE add_var_list_element_l1d
  END INTERFACE add_var

  INTERFACE add_ref
    MODULE PROCEDURE add_var_list_reference_r4d
    MODULE PROCEDURE add_var_list_reference_r3d
    MODULE PROCEDURE add_var_list_reference_r2d
    MODULE PROCEDURE add_var_list_reference_s4d
    MODULE PROCEDURE add_var_list_reference_s3d
    MODULE PROCEDURE add_var_list_reference_s2d
    MODULE PROCEDURE add_var_list_reference_i2d
  END INTERFACE add_ref

  CHARACTER(*), PARAMETER :: modname = "mo_var_list"

CONTAINS

  ! remove all elements of a linked list
  SUBROUTINE delete_list(this)
    CLASS(t_var_list_ptr), INTENT(INOUT) :: this
    INTEGER :: i, n

    IF (ASSOCIATED(this%p)) THEN
      IF (ALLOCATED(this%p%vl)) THEN
        n = SIZE(this%p%vl)
        !$ACC WAIT(1)
        DO i = 1, n
          IF (ASSOCIATED(this%p%vl(i)%p)) THEN
            IF (this%p%vl(i)%p%info%allocated) THEN
              SELECT CASE(this%p%vl(i)%p%info%data_type)
              CASE(REAL_T)
                !$ACC EXIT DATA DELETE(this%p%vl(i)%p%r_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%r_ptr)
              CASE(SINGLE_T)
                !$ACC EXIT DATA DELETE(this%p%vl(i)%p%s_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%s_ptr)
              CASE(INT_T)
                !$ACC EXIT DATA DELETE(this%p%vl(i)%p%i_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%i_ptr)
              CASE(BOOL_T)
                !$ACC EXIT DATA DELETE(this%p%vl(i)%p%l_ptr) IF(this%p%vl(i)%p%info%lopenacc)
                DEALLOCATE(this%p%vl(i)%p%l_ptr)
              END SELECT
            END IF
            DEALLOCATE(this%p%vl(i)%p)
          END IF
        END DO
      END IF
    END IF
  END SUBROUTINE delete_list

  !------------------------------------------------------------------------------------------------
  !
  ! Get a copy of the dynamic metadata concerning a var_list element by index of the element
  !
  SUBROUTINE get_tracer_info_dyn_by_idx (this_list, ncontained, info_dyn)
    !
    TYPE(t_var_list_ptr),         INTENT(in)  :: this_list    ! list
    INTEGER,                      INTENT(in)  :: ncontained   ! index of variable in container
    TYPE(t_var_metadata_dynamic), INTENT(out) :: info_dyn     ! dynamic variable meta data
    !
    TYPE(t_var), POINTER :: element
    !
    element => find_tracer_by_index (this_list, ncontained)
    IF (ASSOCIATED (element)) THEN
      info_dyn = element%info_dyn
    ENDIF
    !
  END SUBROUTINE get_tracer_info_dyn_by_idx


  !-----------------------------------------------------------------------------
  !
  ! Overloaded to search for a tracer by its index (ncontained)
  !
  FUNCTION find_tracer_by_index (this_list, ncontained) RESULT(ret_list_elem)
    !
    TYPE(t_var_list_ptr), INTENT(in) :: this_list
    INTEGER,              INTENT(in) :: ncontained
    !
    INTEGER :: iv
    TYPE(t_var),POINTER :: ret_list_elem
    !
    DO iv=1, this_list%p%nvars
      ret_list_elem => this_list%p%vl(iv)%p
      IF (ret_list_elem%info_dyn%tracer%lis_tracer) THEN
        IF(ncontained == ret_list_elem%info%ncontained) THEN
          RETURN
        ENDIF
      ENDIF
    ENDDO
    !
    NULLIFY (ret_list_elem)
    !
  END FUNCTION find_tracer_by_index

  !-----------------------------------------------------------------------------
  ! add a list element to the linked list
  SUBROUTINE register_list_element(this, varp)
    CLASS(t_var_list_ptr), INTENT(INOUT) :: this
    TYPE(t_var), INTENT(IN), POINTER :: varp
    INTEGER :: iv, na , nv
    TYPE(t_var_ptr), ALLOCATABLE :: vtmp(:)
    CHARACTER(*), PARAMETER :: routine = modname//":register_list_element"
    INTEGER, ALLOCATABLE :: itmp1(:), itmp2(:), itmp3(:), itmp4(:)
    LOGICAL, ALLOCATABLE :: ltmp(:)

    IF (.NOT.ASSOCIATED(this%p)) CALL finish(routine, "not a valid var_list")
    IF (.NOT.ASSOCIATED(varp)) CALL finish(routine, "not a valid var")
    na = 0
    nv = this%p%nvars
    IF (nv .EQ. 0) THEN
      na = 16
    ELSE IF (SIZE(this%p%vl) .EQ. nv) THEN
      na = nv + MAX(8, nv / 8)
    END IF
    IF (na .GT. 0) THEN
      ALLOCATE(itmp1(na), itmp2(na), itmp3(na), itmp4(na), ltmp(na), vtmp(na))
      IF (nv .GT. 0) THEN
        itmp1(1:nv) = this%p%tl(1:nv)
        itmp2(1:nv) = this%p%hgrid(1:nv)
        itmp3(1:nv) = this%p%key(1:nv)
        itmp4(1:nv) = this%p%key_notl(1:nv)
        ltmp(1:nv) = this%p%lout(1:nv)
        DO iv = 1, nv
          vtmp(iv)%p => this%p%vl(iv)%p
        END DO
      END IF
      CALL MOVE_ALLOC(itmp1, this%p%tl)
      CALL MOVE_ALLOC(itmp2, this%p%hgrid)
      CALL MOVE_ALLOC(itmp3, this%p%key)
      CALL MOVE_ALLOC(itmp4, this%p%key_notl)
      CALL MOVE_ALLOC(ltmp, this%p%lout)
      CALL MOVE_ALLOC(vtmp, this%p%vl)
    END IF
    nv = nv + 1
    this%p%vl(nv)%p => varp
    this%p%tl(nv) = get_var_timelevel(varp%info%name)
    this%p%hgrid(nv) = varp%info%hgrid
    this%p%key(nv) = text_hash_c(TRIM(varp%info%name))
    this%p%key_notl(nv) = text_hash_c(tolower(get_var_name(varp%info)))
    this%p%lout(nv) = varp%info%loutput
    this%p%nvars = nv
  END SUBROUTINE register_list_element

  SUBROUTINE inherit_var_list_metadata(this, info)
    CLASS(t_var_list_ptr), INTENT(IN) :: this
    TYPE(t_var_metadata), INTENT(OUT) :: info

    info%grib2               = grib2_var(-1, -1, -1, -1, -1, -1)
    info%lrestart            = this%p%lrestart
    info%lmiss               = this%p%lmiss
    info%lmask_boundary      = this%p%lmask_boundary
    info%vert_interp         = create_vert_interp_metadata()
    info%hor_interp          = create_hor_interp_metadata()
    info%in_group(:)         = groups()
  END SUBROUTINE inherit_var_list_metadata

  !------------------------------------------------------------------------------------------------
  ! Create a list new entry
  SUBROUTINE add_var_list_element_5d(data_type, list, varname, hgrid, vgrid, &
    & cf, grib2, ldims, new_elem, loutput, lcontainer, lrestart,             &
    & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,       &
    & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,       &
    & tracer_info, p5_r, p5_s, p5_i, p5_l, initval_r, initval_s, initval_i,  &
    & initval_l, resetval_r, resetval_s, resetval_i, resetval_l, new_element,&
    & missval_r, missval_s, missval_i, missval_l, var_class, lopenacc)
    INTEGER, INTENT(IN) :: data_type, hgrid, vgrid, ldims(:)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: list
    CHARACTER(*), INTENT(IN) :: varname
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    TYPE(t_var), POINTER, INTENT(OUT) :: new_elem
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, lrestart_cont, &
      & lmiss, in_group(:), initval_l, resetval_l, missval_l, lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, l_pp_scheduler_task, &
      & initval_i, resetval_i, missval_i, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS, TARGET, OPTIONAL :: p5_r(:,:,:,:,:)
    REAL(sp), CONTIGUOUS, TARGET, OPTIONAL :: p5_s(:,:,:,:,:)
    INTEGER, CONTIGUOUS, TARGET, OPTIONAL :: p5_i(:,:,:,:,:)
    LOGICAL, CONTIGUOUS, TARGET, OPTIONAL :: p5_l(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    REAL(dp), INTENT(IN), OPTIONAL :: initval_r, resetval_r, missval_r
    REAL(sp), INTENT(IN), OPTIONAL :: initval_s, resetval_s, missval_s
    TYPE(t_var), POINTER, INTENT(OUT), OPTIONAL :: new_element
    TYPE(t_union_vals) :: missval, initval, resetval, ivals
    INTEGER :: d(5), istat, ndims
    LOGICAL :: referenced, is_restart_var
    CHARACTER(*), PARAMETER :: routine = modname//":add_var_list_element_5d"

    ndims = SIZE(ldims)
    ! Check for a variable of the same name in this list
    ! This consistency check only makes sense inside individual lists.
    ! For single-domain setups and/or when using internal post-processing 
    ! (e.g. lon-lat or vertically interpolated output)  
    ! duplicate names may exist in different lists
    IF (PRESENT(hor_interp)) THEN
      IF (ASSOCIATED(find_list_element(list, varname,                        &
        &                              opt_lonlat_id=hor_interp%lonlat_id))) &
        & CALL finish(routine, "duplicate entry ("//TRIM(varname)//") in var_list ("//TRIM(list%p%vlname)//")")
    ELSE
      IF (ASSOCIATED(find_list_element(list, varname))) &
        & CALL finish(routine, "duplicate entry ("//TRIM(varname)//") in var_list ("//TRIM(list%p%vlname)//")")
    END IF
    is_restart_var = list%p%lrestart
    IF (PRESENT(lrestart)) THEN
      is_restart_var = lrestart
      IF (.NOT.list%p%lrestart .AND. lrestart) &
        & CALL finish(routine, 'for list '//TRIM(list%p%vlname)//' restarting not enabled, '// &
                           & 'but restart of '//TRIM(varname)//' requested.')
    ENDIF
    IF (is_restart_var .AND. (.NOT. ANY(data_type == (/REAL_T, SINGLE_T, INT_T/)))) &
      & CALL finish(routine, 'unsupported data_type for "'//TRIM(varname)//'": '// &
        & 'data_type of restart variables must be floating-point or integer type.')
    ALLOCATE(new_elem)
    CALL inherit_var_list_metadata(list, new_elem%info)
    ! init local fields
    missval = new_elem%info%missval
    initval = new_elem%info%initval
    resetval= new_elem%info%resetval
    ! and set meta data
    referenced = ANY([PRESENT(p5_r), PRESENT(p5_s), PRESENT(p5_i), PRESENT(p5_l)])
    IF (PRESENT(missval_r))  missval%rval  = missval_r
    IF (PRESENT(missval_s))  missval%sval  = missval_s
    IF (PRESENT(missval_i))  missval%ival  = missval_i
    IF (PRESENT(missval_l))  missval%lval  = missval_l
    IF (PRESENT(initval_r))  initval%rval  = initval_r
    IF (PRESENT(initval_s))  initval%sval  = initval_s
    IF (PRESENT(initval_i))  initval%ival  = initval_i
    IF (PRESENT(initval_l))  initval%lval  = initval_l
    IF (PRESENT(resetval_r)) resetval%rval = resetval_r
    IF (PRESENT(resetval_s)) resetval%sval = resetval_s
    IF (PRESENT(resetval_i)) resetval%ival = resetval_i
    IF (PRESENT(resetval_l)) resetval%lval = resetval_l
    CALL set_var_metadata(new_elem%info, ldims, name=varname,       &
      & hgrid=hgrid, vgrid=vgrid, cf=cf, grib2=grib2, loutput=loutput,       &
      & lcontainer=lcontainer, lrestart=lrestart, missval=missval,           &
      & lrestart_cont=lrestart_cont, initval=initval, isteptype=isteptype,   &
      & resetval=resetval, tlev_source=tlev_source, vert_interp=vert_interp, &
      & hor_interp=hor_interp, l_pp_scheduler_task=l_pp_scheduler_task,      &
      & post_op=post_op, action_list=action_list, var_class=var_class,       &
      & data_type=data_type, lopenacc=lopenacc, lmiss=lmiss, in_group=in_group)
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn (new_elem%info_dyn, tracer_info=tracer_info)
    new_elem%info%ndims = ndims
    new_elem%info%used_dimensions(1:ndims) = ldims(1:ndims)
    new_elem%info%dom = list%p%patch_id
    IF(PRESENT(info)) info => new_elem%info
    NULLIFY(new_elem%r_ptr, new_elem%s_ptr, new_elem%i_ptr, new_elem%l_ptr)
    d(1:ndims)    = new_elem%info%used_dimensions(1:ndims)
    d((ndims+1):) = 1
#if    defined (VARLIST_INITIZIALIZE_WITH_NAN) \
    && (defined (__INTEL_COMPILER) || defined (__PGI) || defined (NAGFOR))
    ivals%rval = ieee_value(ptr, ieee_signaling_nan)
    ivals%sval = ieee_value(ptr, ieee_signaling_nan)
#endif
    IF (ANY([PRESENT(initval_r), PRESENT(initval_s), PRESENT(initval_i), PRESENT(initval_l)])) THEN
      ivals = initval
    ELSE IF (PRESENT(lmiss)) THEN
      ivals = missval
    END IF
    SELECT CASE(data_type)
    CASE(REAL_T)
      IF (referenced) THEN
        new_elem%r_ptr => p5_r
      ELSE
        new_elem%var_base_size = 8
        ALLOCATE(new_elem%r_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%r_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_dp(new_elem%r_ptr, PRODUCT(d(1:5)), ivals%rval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%r_ptr) ASYNC(1) IF(new_elem%info%lopenacc)
    CASE(SINGLE_T)
      IF (referenced) THEN
        new_elem%s_ptr => p5_s
      ELSE
        new_elem%var_base_size = 4
        ALLOCATE(new_elem%s_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%s_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_sp(new_elem%s_ptr, PRODUCT(d(1:5)), ivals%sval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%s_ptr) ASYNC(1) IF(new_elem%info%lopenacc)
    CASE(INT_T)
      IF (referenced) THEN
        new_elem%i_ptr => p5_i
      ELSE
        new_elem%var_base_size = 4
        ALLOCATE(new_elem%i_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of arrayb'//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%i_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_i4(new_elem%i_ptr, PRODUCT(d(1:5)), ivals%ival)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%i_ptr) ASYNC(1) IF(new_elem%info%lopenacc)
    CASE(BOOL_T)
      IF (referenced) THEN
        new_elem%l_ptr => p5_l
      ELSE
        new_elem%var_base_size = 4
        ALLOCATE(new_elem%l_ptr(d(1), d(2), d(3), d(4), d(5)), STAT=istat)
        IF (istat /= 0) CALL finish(routine, 'allocation of array '//TRIM(varname)//' failed')
        !$ACC ENTER DATA CREATE(new_elem%l_ptr) IF(new_elem%info%lopenacc)
      END IF
      !ICON_OMP PARALLEL
      CALL init_contiguous_l(new_elem%l_ptr, PRODUCT(d(1:5)), ivals%lval)
      !ICON_OMP END PARALLEL
      !$ACC UPDATE DEVICE(new_elem%l_ptr) ASYNC(1) IF(new_elem%info%lopenacc)
    END SELECT
    CALL register_list_element(list, new_elem)
    IF (.NOT.referenced) list%p%memory_used = list%p%memory_used + &
      & INT(new_elem%var_base_size, i8) * INT(PRODUCT(d(1:5)),i8)
    new_elem%info%allocated = .TRUE.
    IF (PRESENT(new_element)) new_element => new_elem
  END SUBROUTINE add_var_list_element_5d

  SUBROUTINE add_var_list_element_r4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%r_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_r4d

  SUBROUTINE add_var_list_element_r3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta), INTENT(in), OPTIONAL :: tracer_info
    REAL(dp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%r_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_r3d

  SUBROUTINE add_var_list_element_r2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    REAL(dp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%r_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_r2d

  SUBROUTINE add_var_list_element_r1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(dp), POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(dp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(REAL_T, this_list, varname, hgrid, vgrid, &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_r=p5, initval_r=initval, resetval_r=resetval, missval_r=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%r_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_r1d

  SUBROUTINE add_var_list_element_s4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(sp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,&
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%s_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_s4d

  SUBROUTINE add_var_list_element_s3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(:)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info
    REAL(sp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,&
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%s_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_s3d

  SUBROUTINE add_var_list_element_s2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element, tracer_info, &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    CLASS(t_tracer_meta),    INTENT(in), OPTIONAL :: tracer_info
    REAL(sp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,&
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element,   &
      & tracer_info=tracer_info)
    ptr => element%s_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_s2d

  SUBROUTINE add_var_list_element_s1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    REAL(sp), POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    REAL(sp), CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(SINGLE_T, this_list, varname, hgrid, vgrid,& 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_s=p5, initval_s=initval, resetval_s=resetval, missval_s=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%s_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_s1d

  SUBROUTINE add_var_list_element_i4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_i4d

  SUBROUTINE add_var_list_element_i3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_i3d

  SUBROUTINE add_var_list_element_i2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_i2d

  SUBROUTINE add_var_list_element_i1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, lmiss, in_group(:), lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class, initval, resetval, missval
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    INTEGER, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(INT_T, this_list, varname, hgrid, vgrid,  &
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_i=p5, initval_i=initval, resetval_i=resetval, missval_i=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%i_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_i1d

  SUBROUTINE add_var_list_element_l4d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,:,:,:,1)
  END SUBROUTINE add_var_list_element_l4d

  SUBROUTINE add_var_list_element_l3d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,:,:,1,1)
  END SUBROUTINE add_var_list_element_l3d

  SUBROUTINE add_var_list_element_l2d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,:,1,1,1)
  END SUBROUTINE add_var_list_element_l2d

  SUBROUTINE add_var_list_element_l1d(this_list, varname, ptr, hgrid, vgrid, &
    & cf, grib2, ldims, loutput, lcontainer, lrestart, lrestart_cont,     &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, info,    &
    & p5, vert_interp, hor_interp, in_group, new_element,        &
    & l_pp_scheduler_task, post_op, action_list, var_class, lopenacc)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: varname
    LOGICAL, POINTER, INTENT(OUT) :: ptr(:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ldims(1)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lcontainer, lrestart, &
      & lrestart_cont, initval, resetval, lmiss, missval, in_group(:), &
      & lopenacc
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, l_pp_scheduler_task, &
      & tlev_source, var_class
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    LOGICAL, CONTIGUOUS, TARGET, OPTIONAL :: p5(:,:,:,:,:)
    TYPE(t_vert_interp_meta), INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var), POINTER :: element

    CALL add_var_list_element_5d(BOOL_T, this_list, varname, hgrid, vgrid, & 
      & cf, grib2, ldims, element, loutput, lcontainer, lrestart,          &
      & lrestart_cont, isteptype, lmiss, tlev_source, info, vert_interp,   &
      & hor_interp, in_group, l_pp_scheduler_task, post_op, action_list,   &
      & p5_l=p5, initval_l=initval, resetval_l=resetval, missval_l=missval,&
      & var_class=var_class, lopenacc=lopenacc, new_element=new_element)
    ptr => element%l_ptr(:,1,1,1,1)
  END SUBROUTINE add_var_list_element_l1d

  SUBROUTINE add_var_list_reference_util(target_element, new_list_element,      &
    & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, ldims, &
    & dtype, icontainer, vrp, loutput, lrestart, lrestart_cont, isteptype,      &
    & lmiss, tlev_source, tracer_info, info, vert_interp, hor_interp, in_group, &
    & new_element, l_pp_scheduler_task, post_op, action_list, idx_diag,         &
    & var_class, opt_var_ref_pos, initval_r, initval_s, initval_i, missval_r,   &
    & missval_s, missval_i, resetval_r, resetval_s, resetval_i, idx_tracer)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    TYPE(t_var), INTENT(OUT), POINTER :: target_element, new_list_element
    CHARACTER(*), INTENT(IN) :: target_name, refname
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(:), dtype
    INTEGER, INTENT(OUT) :: icontainer, vrp
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, idx_tracer, idx_diag
    REAL(dp), INTENT(IN), OPTIONAL :: initval_r, resetval_r, missval_r
    REAL(sp), INTENT(IN), OPTIONAL :: initval_s, resetval_s, missval_s
    INTEGER, INTENT(IN), OPTIONAL :: initval_i, resetval_i, missval_i
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    TYPE(t_var_metadata), POINTER :: target_info, ref_info
    TYPE(t_union_vals) :: missvalt, initvalt, resetvalt
    CHARACTER(*), PARAMETER :: routine = modname//":add_var_list_reference_util"
    INTEGER :: var_ref_pos, ndims, di(5), di3, max_ref, ts_pos

    ndims = SIZE(ldims)
    target_element => find_list_element(this_list, target_name)
    IF (.NOT. ASSOCIATED(target_element)) THEN
      CALL finish(routine, "target element "//TRIM(target_name)//" not found!")
    END IF
    target_info => target_element%info
    IF (PRESENT(opt_var_ref_pos)) THEN
      var_ref_pos = opt_var_ref_pos
      IF (.NOT. target_info%lcontainer) &
        &  CALL finish(routine, "invalid container index: Target is not a container variable!")
      IF ((target_info%var_ref_pos /= var_ref_pos) .AND. &
        & (target_info%var_ref_pos /= -1)) THEN
        CALL finish(routine, "Container index does not match the previously set value!")
      END IF
      target_info%var_ref_pos = var_ref_pos
    ELSE
      var_ref_pos = ndims + 1
    END IF
    di3 = MERGE(4, 0, ndims.EQ.3)
   IF (.NOT. ANY(NDIMS .EQ. (/ 2, 3, 4 /))) CALL finish(routine, "Internal error 1!")
    SELECT CASE(var_ref_pos)
    CASE(1)
      di = (/ 2, 3, di3, 0, 0 /)
    CASE(2)
      di = (/ 1, 3, di3, 0, 0 /)
    CASE(3)
      di = (/ 1, 2, di3, 0, 0 /)
    CASE(4)
      IF (NDIMS.EQ.2) CALL finish(routine, "Internal error 2!")
      di = (/ 1, 2, 3, 0, 0 /)
    CASE(5)
      IF (NDIMS.EQ.2 .OR. NDIMS.EQ.3) CALL finish(routine, "Internal error 3!")
      di = (/ 1, 2, 3, 4, 0 /)
    CASE DEFAULT
      CALL finish(routine, "Internal error 4!")
    END SELECT
    IF (target_info%lcontainer) THEN
      max_ref = 0
      IF (ASSOCIATED(target_element%r_ptr)) THEN
        max_ref = SIZE(target_element%r_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%s_ptr)) THEN
        max_ref = SIZE(target_element%s_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%i_ptr)) THEN
        max_ref = SIZE(target_element%i_ptr, var_ref_pos)
      ELSE IF (ASSOCIATED(target_element%l_ptr)) THEN
        max_ref = SIZE(target_element%l_ptr, var_ref_pos)
      END IF
      ! Counting the number of existing references is deactivated, 
      ! if the slice index to be referenced is given explicitly.
        target_info%ncontained = target_info%ncontained+1
        ! only check validity of given slice index
        IF ( (ref_idx > max_ref) .OR. (ref_idx < 1)) THEN
          WRITE (message_text, "(2(a,i3),a)") 'Slice idx ', ref_idx, ' for ' // &
            & TRIM(refname) // ' out of allowable range [1,',max_ref,']'
          CALL finish(routine, message_text)
        ENDIF
      IF (ANY(ldims(1:ndims) /= target_info%used_dimensions(di(1:ndims)))) &
        & CALL finish(routine, TRIM(refname)//' dimensions requested and available differ.')
    ENDIF
    ! add list entry
    ALLOCATE(new_list_element)
    IF (PRESENT(new_element)) new_element=>new_list_element
    new_list_element%ref_to => target_element
    ref_info => new_list_element%info
    CALL inherit_var_list_metadata(this_list, ref_info)
    ! init local fields
    missvalt  = ref_info%missval
    initvalt  = ref_info%initval
    resetvalt = ref_info%resetval
    IF (PRESENT(missval_r))  missvalt%rval  = missval_r
    IF (PRESENT(missval_s))  missvalt%sval  = missval_s
    IF (PRESENT(missval_i))  missvalt%ival  = missval_i
    IF (PRESENT(initval_r))  initvalt%rval  = initval_r
    IF (PRESENT(initval_s))  initvalt%sval  = initval_s
    IF (PRESENT(initval_i))  initvalt%ival  = initval_i
    IF (PRESENT(resetval_r)) resetvalt%rval = resetval_r
    IF (PRESENT(resetval_s)) resetvalt%sval = resetval_s
    IF (PRESENT(resetval_i)) resetvalt%ival = resetval_i
    CALL set_var_metadata(ref_info, ldims, name=refname, hgrid=hgrid, &
      & vgrid=vgrid, cf=cf, grib2=grib2, loutput=loutput, lrestart=lrestart,   &
      & missval=missvalt, lrestart_cont=lrestart_cont, initval=initvalt,       &
      & isteptype=isteptype, resetval=resetvalt, tlev_source=tlev_source,      &
      & vert_interp=vert_interp, hor_interp=hor_interp, post_op=post_op,       &
      & l_pp_scheduler_task=l_pp_scheduler_task, action_list=action_list,      &
      & var_class=var_class, data_type=dtype, lmiss=lmiss, in_group=in_group,  &
      & idx_tracer=idx_tracer, idx_diag=idx_diag, lopenacc=target_info%lopenacc)
    ! set dynamic metadata, i.e. polymorphic tracer metadata
    CALL set_var_metadata_dyn(new_list_element%info_dyn, tracer_info=tracer_info)
    ref_info%ndims = ndims
    ref_info%used_dimensions(:) = 0
    ref_info%used_dimensions(1:ndims) = target_info%used_dimensions(di(1:ndims))
    ref_info%dom = target_info%dom
    IF (PRESENT(var_class)) THEN
      IF (ANY((/CLASS_TILE, CLASS_TILE_LAND/) == var_class)) THEN
        ! automatically add tile to its variable specific tile-group
        ts_pos = INDEX(target_info%name, TIMELEVEL_SUFFIX)
        ts_pos = MERGE(ts_pos-1, vname_len, ts_pos .GT. 0)
        IF (PRESENT(in_group)) then
          ref_info%in_group = groups(target_info%name(1:ts_pos), groups_in=in_group)
        ELSE
          ref_info%in_group = groups(target_info%name(1:ts_pos))
        END IF
      END IF
    END IF
    IF (target_info%lcontainer) THEN
      ref_info%lcontained                   = .TRUE.
      ref_info%used_dimensions(ndims+1)     = 1
      ref_info%var_ref_pos = var_ref_pos
      ref_info%maxcontained = max_ref
      ref_info%ncontained = ref_idx
    ENDIF
    icontainer = MERGE(ref_info%ncontained, 1, target_info%lcontainer)
    vrp = var_ref_pos
    IF(PRESENT(info)) info => ref_info
    CALL register_list_element(this_list, new_list_element)
  END SUBROUTINE add_var_list_reference_util

  SUBROUTINE add_var_list_reference_r4d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r4d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, vrp, loutput=loutput, lrestart=lrestart,&
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval)
    IF (.NOT. ASSOCIATED(target_element%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%r_ptr(icontainer,:,:,:,:)
    CASE(2)
      ptr => target_element%r_ptr(:,icontainer,:,:,:)
    CASE(3)
      ptr => target_element%r_ptr(:,:,icontainer,:,:)
    CASE(4)
      ptr => target_element%r_ptr(:,:,:,icontainer,:)
    CASE(5)
      ptr => target_element%r_ptr(:,:,:,:,icontainer)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%r_ptr => target_element%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%rval
  END SUBROUTINE add_var_list_reference_r4d

  SUBROUTINE add_var_list_reference_r3d(this_list, target_name, refname, ptr,       &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r3d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, vrp, loutput=loutput, lrestart=lrestart,&
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval)
    IF (.NOT. ASSOCIATED(target_element%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%r_ptr(icontainer,:,:,:,1)
    CASE(2)
      ptr => target_element%r_ptr(:,icontainer,:,:,1)
    CASE(3)
      ptr => target_element%r_ptr(:,:,icontainer,:,1)
    CASE(4)
      ptr => target_element%r_ptr(:,:,:,icontainer,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%r_ptr => target_element%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%rval
  END SUBROUTINE add_var_list_reference_r3d

  SUBROUTINE add_var_list_reference_r2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class,     &
    & idx_tracer, idx_diag)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(dp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf                  
    TYPE(t_grib2_var), INTENT(IN) :: grib2               
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, idx_tracer, idx_diag
    REAL(dp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info         
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info                
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp         
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp          
    TYPE(t_var), POINTER, OPTIONAL :: new_element         
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op            
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list         
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_r2d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, REAL_T, icontainer, vrp, loutput=loutput, lrestart=lrestart,&
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_r=initval,                &
      & missval_r=missval, resetval_r=resetval, idx_tracer=idx_tracer,     &
      & idx_diag=idx_diag)
    IF (.NOT. ASSOCIATED(target_element%r_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%r_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%r_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%r_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%r_ptr => target_element%r_ptr
    IF (.NOT. ASSOCIATED(new_list_element%r_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%rval
  END SUBROUTINE add_var_list_reference_r2d

  SUBROUTINE add_var_list_reference_s4d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(4)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s4d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, vrp, loutput=loutput, lrestart=lrestart,&
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%s_ptr(icontainer,:,:,:,:)
    CASE(2)
      ptr => target_element%s_ptr(:,icontainer,:,:,:)
    CASE(3)
      ptr => target_element%s_ptr(:,:,icontainer,:,:)
    CASE(4)
      ptr => target_element%s_ptr(:,:,:,icontainer,:)
    CASE(5)
      ptr => target_element%s_ptr(:,:,:,:,icontainer)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%s_ptr => target_element%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%sval
  END SUBROUTINE add_var_list_reference_s4d

  SUBROUTINE add_var_list_reference_s3d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(3)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s3d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, vrp, loutput=loutput, lrestart=lrestart, &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%s_ptr(icontainer,:,:,:,1)
    CASE(2)
      ptr => target_element%s_ptr(:,icontainer,:,:,1)
    CASE(3)
      ptr => target_element%s_ptr(:,:,icontainer,:,1)
    CASE(4)
      ptr => target_element%s_ptr(:,:,:,icontainer,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%s_ptr => target_element%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%sval
  END SUBROUTINE add_var_list_reference_s3d

  SUBROUTINE add_var_list_reference_s2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    REAL(sp), POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos
    REAL(sp), INTENT(IN), OPTIONAL :: initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_s2d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, SINGLE_T, icontainer, vrp, loutput=loutput, lrestart=lrestart, &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_s=initval,                &
      & missval_s=missval, resetval_s=resetval)
    IF (.NOT. ASSOCIATED(target_element%s_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%s_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%s_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%s_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%s_ptr => target_element%s_ptr
    IF (.NOT. ASSOCIATED(new_list_element%s_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%sval
  END SUBROUTINE add_var_list_reference_s2d

  SUBROUTINE add_var_list_reference_i2d(this_list, target_name, refname, ptr,    &
    & hgrid, vgrid, cf, grib2, ref_idx, ldims, loutput, lrestart, lrestart_cont, &
    & initval, isteptype, resetval, lmiss, missval, tlev_source, tracer_info,    &
    & info, vert_interp, hor_interp, in_group, new_element,             &
    & l_pp_scheduler_task, post_op, action_list, opt_var_ref_pos, var_class)
    TYPE(t_var_list_ptr), INTENT(INOUT) :: this_list
    CHARACTER(*), INTENT(IN) :: target_name, refname
    INTEGER, POINTER :: ptr(:,:)
    INTEGER, INTENT(IN) :: hgrid, vgrid, ref_idx, ldims(2)
    TYPE(t_cf_var), INTENT(IN) :: cf
    TYPE(t_grib2_var), INTENT(IN) :: grib2
    LOGICAL, INTENT(IN), OPTIONAL :: loutput, lrestart, lrestart_cont, &
      & lmiss, in_group(:)
    INTEGER, INTENT(IN), OPTIONAL :: isteptype, tlev_source, var_class, &
      & l_pp_scheduler_task, opt_var_ref_pos, initval, resetval, missval
    CLASS(t_tracer_meta), INTENT(IN), OPTIONAL :: tracer_info
    TYPE(t_var_metadata), POINTER, OPTIONAL :: info
    TYPE(t_vert_interp_meta),INTENT(IN), OPTIONAL :: vert_interp
    TYPE(t_hor_interp_meta), INTENT(IN), OPTIONAL :: hor_interp
    TYPE(t_var), POINTER, OPTIONAL :: new_element
    TYPE(t_post_op_meta), INTENT(IN), OPTIONAL :: post_op
    TYPE(t_var_action), INTENT(IN), OPTIONAL :: action_list
    CHARACTER(*), PARAMETER :: routine = modname//"::add_var_list_reference_i2d"
    TYPE(t_var), POINTER :: target_element, new_list_element
    INTEGER :: icontainer, vrp

    CALL add_var_list_reference_util(target_element, new_list_element,     &
      & this_list, target_name, refname, hgrid, vgrid, cf, grib2, ref_idx, &
      & ldims, INT_T, icontainer, vrp, loutput=loutput, lrestart=lrestart, &
      & lrestart_cont=lrestart_cont, isteptype=isteptype, lmiss=lmiss,     &
      & tlev_source=tlev_source, tracer_info=tracer_info, info=info,       &
      & vert_interp=vert_interp, hor_interp=hor_interp, in_group=in_group, &
      & new_element=new_element, l_pp_scheduler_task=l_pp_scheduler_task,  &
      & post_op=post_op, action_list=action_list, var_class=var_class,     &
      & opt_var_ref_pos=opt_var_ref_pos, initval_i=initval,                &
      & missval_i=missval, resetval_i=resetval)
    IF (.NOT. ASSOCIATED(target_element%i_ptr)) &
      & CALL finish(routine, TRIM(refname)//' not created.')
    SELECT CASE(vrp)
    CASE(1)
      ptr => target_element%i_ptr(icontainer,:,:,1,1)
    CASE(2)
      ptr => target_element%i_ptr(:,icontainer,:,1,1)
    CASE(3)
      ptr => target_element%i_ptr(:,:,icontainer,1,1)
    CASE default
      CALL finish(routine, "internal error!")
    END SELECT
    new_list_element%i_ptr => target_element%i_ptr
    IF (.NOT. ASSOCIATED(new_list_element%i_ptr)) &
      & WRITE (0,*) 'problem with association of ptr for '//TRIM(refname)
    IF (PRESENT(lmiss)) ptr = new_list_element%info%missval%ival
  END SUBROUTINE add_var_list_reference_i2d

  SUBROUTINE print_var_list(this, lshort)
    CLASS(t_var_list_ptr), INTENT(IN) :: this
    LOGICAL, OPTIONAL, INTENT(IN) :: lshort
    TYPE(t_var), POINTER :: le
    INTEGER :: j
    LOGICAL :: short

    short = .FALSE.
    IF (PRESENT(lshort)) short = lshort
    CALL message('','')
    CALL message('','')
    CALL message('','Status of variable list '//TRIM(this%p%vlname)//':')
    CALL message('','')
    DO j = 1, this%p%nvars
      le => this%p%vl(j)%p
      IF (le%info%lcontainer) CYCLE
      IF (le%info%name == '') CYCLE
      IF (short) THEN
        CALL le%print_short()
      ELSE
        CALL le%print_rigorous()
      END IF
    END DO
  END SUBROUTINE print_var_list

  FUNCTION find_list_element(this, vname, opt_hgrid, opt_with_tl, opt_output, &
    &                        opt_lonlat_id)  RESULT(element)
    TYPE(t_var_list_ptr), INTENT(IN) :: this
    CHARACTER(*), INTENT(IN) :: vname
    INTEGER, OPTIONAL, INTENT(IN) :: opt_hgrid
    LOGICAL, OPTIONAL, INTENT(IN) :: opt_with_tl, opt_output
    !> optional: distinguish variables of different lon-lat output grids
    INTEGER, OPTIONAL, INTENT(IN) :: opt_lonlat_id
    TYPE(t_var), POINTER :: element
    INTEGER :: key, hgrid, time_lev, iv
    LOGICAL :: with_tl, omit_output, with_output

    NULLIFY(element)
    with_tl = .TRUE.
    IF (PRESENT(opt_with_tl)) with_tl = opt_with_tl
    hgrid = -1
    IF (PRESENT(opt_hgrid)) hgrid = opt_hgrid
    IF (with_tl) THEN
      key = text_hash_c(TRIM(vname))
    ELSE
      key = text_hash_c(tolower(vname))
    END IF
    with_output = .TRUE.
    omit_output = .NOT.PRESENT(opt_output)
    IF (.NOT.omit_output) with_output = opt_output
    time_lev = get_var_timelevel(vname)
    DO iv = 1, this%p%nvars
      IF (-1 .NE. hgrid .AND. this%p%hgrid(iv) .NE. hgrid) CYCLE
      IF (MERGE(.FALSE., with_output .NEQV. this%p%lout(iv), omit_output)) CYCLE
      IF (time_lev .NE. this%p%tl(iv)) CYCLE
      IF (key .NE. MERGE(this%p%key(iv), this%p%key_notl(iv), with_tl)) CYCLE
      IF (PRESENT(opt_lonlat_id)) THEN
         IF (this%p%vl(iv)%p%info%hor_interp%lonlat_id /= opt_lonlat_id) CYCLE
      END IF
      element => this%p%vl(iv)%p
    ENDDO
  END FUNCTION find_list_element
  
END MODULE mo_var_list
