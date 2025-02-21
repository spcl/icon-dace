! Routines for handling proxy variables e.g. accumulation buffers
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

MODULE mo_derived_variable_handling

  USE mo_kind,                ONLY: wp, sp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_dynamics_config,     ONLY: nnow, nnew, nold
  USE mo_impl_constants,      ONLY: vname_len, REAL_T, TIMELEVEL_SUFFIX, max_dom
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
                              & GRID_ZONAL, GRID_UNSTRUCTURED_VERT
  USE mo_name_list_output_types, ONLY: t_output_name_list
  USE mo_name_list_output_metadata, ONLY: metainfo_get_timelevel
  USE mo_var_list,            ONLY: add_var, t_var_list_ptr
  USE mo_var,                 ONLY: t_var, t_var_ptr
  USE mo_var_metadata,        ONLY: get_var_name
  USE mo_var_metadata_types,  ONLY: t_var_metadata
  USE mo_var_list_register_utils, ONLY: vlr_find
  USE mo_exception,           ONLY: finish
  USE mtime,                  ONLY: newEvent, event, isCurrentEventActive, datetime
  USE mo_time_config,         ONLY: time_config
  USE mo_cdi,                 ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_LONLAT, TSTEP_CONSTANT
  USE mo_util_texthash,       ONLY: text_hash_c
  USE mo_mpi,                 ONLY: p_bcast, i_am_accel_node
  USE ISO_C_BINDING,          ONLY: C_INT

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

  PUBLIC :: init_statistics, update_statistics, statistics_active_on_dom

  TYPE :: t_derivate_var
    TYPE(t_var_ptr) :: dst, src(3)
    INTEGER :: tls(3) = -999, counter = 0
  END TYPE t_derivate_var

  TYPE :: t_derivate_var_alloctble
    TYPE(t_derivate_var), ALLOCATABLE :: a
  END TYPE t_derivate_var_alloctble

  TYPE :: t_derivate_event
    TYPE(t_derivate_var_alloctble), ALLOCATABLE :: vars(:)
    TYPE(event), POINTER :: mtime_event
    CHARACTER(LEN=vname_len) :: eString
    INTEGER :: eKey = 0
  END TYPE t_derivate_event

  TYPE :: t_derivate_event_alloctble
    TYPE(t_derivate_event), ALLOCATABLE :: a
  END TYPE t_derivate_event_alloctble

  TYPE :: t_derivate_op
    TYPE(t_derivate_event_alloctble), ALLOCATABLE :: events(:)
  END TYPE t_derivate_op

  INTEGER, PARAMETER :: nops = 5
  TYPE(t_derivate_op), TARGET :: ops(nops)
  CHARACTER(*), PARAMETER :: dlim = '|'
  ENUM, BIND(C)
    ENUMERATOR :: O_MEAN = 1, O_MAX, O_MIN, O_MEAN_SQ, O_ACC
    ENUMERATOR :: F_ACC = 100, F_MAX, F_MIN, F_ACC_SQ, F_ASS, F_ASS_SQ, F_WGT, F_MISS
  END ENUM
  INTEGER(C_INT), PARAMETER :: opcodes(nops) = [O_MEAN, O_MAX, O_MIN, O_MEAN_SQ, O_ACC]
  CHARACTER(*), PARAMETER :: opnames(nops)   = ["mean  ", "max   ", "min   ", "square", "acc   "]
  INTEGER, PARAMETER :: oplen(nops)          = [4, 3, 3, 6, 3]
  ! internal processing rules: which function to use for which operator
  INTEGER(C_INT), PARAMETER :: opfuncs(nops) = [F_ACC, F_MAX, F_MIN, F_ACC_SQ, F_ACC]
  CHARACTER(*), PARAMETER :: modname = 'mo_derived_variable_handling'

  !> Flag indicating that statistics accumulation is active on a domain.
  LOGICAL :: statistics_active_on_dom(max_dom) = .FALSE.

CONTAINS

  ! wire up namelist mvstream associations
  SUBROUTINE init_statistics(p_onl, need_bc, root, comm, patch_2d)
    TYPE(t_output_name_list), POINTER, INTENT(IN) :: p_onl
    TYPE(t_patch), INTENT(IN), OPTIONAL :: patch_2d
    LOGICAL, INTENT(IN) :: need_bc
    INTEGER, INTENT(IN) :: root, comm
    TYPE t_vl_arr
      CHARACTER(LEN=vname_len), POINTER :: p(:)
    END TYPE t_vl_arr
    TYPE(t_vl_arr) :: vls(4)
    INTEGER :: iop, jop, loplen, ivl

    iop = 0
    loplen = LEN_TRIM(p_onl%operation)
    DO jop = 1, nops
      IF (loplen .NE. oplen(jop)) CYCLE
      IF (p_onl%operation(1:loplen) /= opnames(jop)(1:oplen(jop))) CYCLE
      iop = jop
      EXIT
    END DO
    vls(1)%p => p_onl%ml_varlist
    vls(2)%p => p_onl%pl_varlist
    vls(3)%p => p_onl%hl_varlist
    vls(4)%p => p_onl%il_varlist
    IF (iop .NE. 0) THEN
      statistics_active_on_dom(p_onl%dom) = .TRUE.
      DO ivl = 1, 4
        IF (vls(ivl)%p(1)(1:1) /= ' ') THEN
          IF (PRESENT(patch_2d)) &
            & CALL init_op(iop, p_onl, vls(ivl)%p, patch_2d)
          IF (need_bc) CALL p_bcast(vls(ivl)%p, root, comm)
        END IF
      END DO
    END IF
  END SUBROUTINE init_statistics

  SUBROUTINE init_op(iop, p_onl, vlist, patch_2d)
    INTEGER, INTENT(IN) :: iop
    TYPE(t_output_name_list), TARGET :: p_onl
    CHARACTER(LEN=vname_len), INTENT(INOUT) :: vlist(:)
    TYPE(t_patch), INTENT(IN), TARGET :: patch_2d
    INTEGER :: eKey, ie, iv, nv_scan, nv_new, nv_old, ne, dns_len, it_len, st_len
    TYPE(t_derivate_event), POINTER :: ederiv
    TYPE(t_derivate_event_alloctble), ALLOCATABLE :: etmp(:)
    TYPE(t_derivate_var_alloctble), ALLOCATABLE :: vderiv(:), vtmp(:)
    TYPE(t_var), POINTER :: vl_elem
    CHARACTER(LEN=vname_len) :: dname, eString, dname_suffix
    TYPE(t_var_list_ptr) :: src_list
    CHARACTER(*), PARAMETER :: routine = modname//":init_op"

    IF (ANY(1 < [p_onl%stream_partitions_ml, p_onl%stream_partitions_pl,  &
      &          p_onl%stream_partitions_hl, p_onl%stream_partitions_il])) &
      & CALL finish(routine, "only supported on global domain 1 " // &
      &                      "and without stream partitioning!")
    nv_scan = 1
    DO WHILE(.NOT.(vlist(nv_scan + 1)(1:1) == ' '))
      nv_scan = nv_scan + 1
    END DO
    ALLOCATE(vderiv(nv_scan))
    ! uniq identifier for an event based on output start/end/interval
    st_len = LEN_TRIM(p_onl%output_start(1))
    it_len = LEN_TRIM(p_onl%output_interval(1))
    WRITE(eString, "(5a)") p_onl%output_start(1)(:st_len), '_', &
      & TRIM(p_onl%output_end(1)), '_', p_onl%output_interval(1)(:it_len)
    WRITE(dname_suffix, "(8a,i0)") dlim, opnames(iop)(1:oplen(iop)), dlim, &
        & p_onl%output_interval(1)(:it_len), dlim, &
        & p_onl%output_start(1)(:st_len), dlim, 'DOM', p_onl%dom
    dns_len = LEN_TRIM(dname_suffix)
    ! this has the advantage that we can compute a uniq id without creating the event itself
    ! fill main dictionary of variables for different event
    eKey = text_hash_c(eString)
    NULLIFY(ederiv)
    ne = 0
    IF (ALLOCATED(ops(iop)%events)) ne = SIZE(ops(iop)%events)
    DO ie = 2, ne
      IF (ops(iop)%events(ie)%a%eKey .EQ. eKey) THEN
        IF (ops(iop)%events(ie)%a%eString /= eString) THEN
          ederiv => ops(iop)%events(ie)%a
          EXIT
        END IF
      END IF
    END DO
    nv_old = 0
    IF (.NOT.ASSOCIATED(ederiv)) THEN
      IF (ne .GT. 0) CALL MOVE_ALLOC(ops(iop)%events, etmp)
      ALLOCATE(ops(iop)%events(ne + 1))
      DO ie = 1, ne
        CALL MOVE_ALLOC(etmp(ie)%a, ops(iop)%events(ie)%a)
      END DO
      ne = ne + 1
      ALLOCATE(ops(iop)%events(ne)%a)
      ederiv => ops(iop)%events(ne)%a
      ederiv%eString = eString
      ederiv%eKey = eKey
      ederiv%mtime_event => newEvent(eString, p_onl%output_start(1), &
        & p_onl%output_start(1), p_onl%output_end(1), p_onl%output_interval(1))
    ELSE
      IF (ALLOCATED(ederiv%vars)) THEN
        CALL MOVE_ALLOC(ederiv%vars, vtmp)
        nv_old = SIZE(vtmp)
      END IF
    END IF
    nv_new = 0
    DO iv = 1, nv_scan
      ! collect data variables only
      IF (INDEX(vlist(iv),':') > 0) CYCLE ! to avoid e.g. "grid:clon" stuff
      ! check for already created meanStream variable (maybe from another output_nml with the same output_interval)
      ! names consist of original spot-value names PLUS event information (start + interval of output)
      WRITE(dname, "(2a)") TRIM(vlist(iv)), dname_suffix(:dns_len)
      vl_elem => vlr_find(dname, opt_patch_id=p_onl%dom)
      IF (.NOT.ASSOCIATED(vl_elem)) THEN !not found -->> create a new one
        ALLOCATE(vderiv(nv_new+1)%a) ! staging a new var entry
        CALL find_src_element(TRIM(vlist(iv)), vderiv(nv_new+1)%a)
        IF (TSTEP_CONSTANT .EQ. vderiv(nv_new+1)%a%src(1)%p%info%isteptype) THEN
          DEALLOCATE(vderiv(nv_new+1)%a) ! discard staged entry
          dname = vlist(iv) ! no aggregation needed, since constant
        ELSE
          nv_new = nv_new + 1 ! new entry is valid, so keep
          CALL copy_var_to_list(vderiv(nv_new)%a) ! add_var to store accumulation
        END IF
      END IF
      vlist(iv) = dname
    END DO
    ALLOCATE(ederiv%vars(nv_new + nv_old))
    DO iv = 1, nv_old
      CALL MOVE_ALLOC(vtmp(iv)%a, ederiv%vars(iv)%a)
    END DO
    DO iv = 1, nv_new
      CALL MOVE_ALLOC(vderiv(iv)%a, ederiv%vars(iv+nv_old)%a)
    END DO
  CONTAINS

  SUBROUTINE find_src_element(vname, deriv)
    CHARACTER(*), INTENT(IN) :: vname
    TYPE(t_derivate_var), INTENT(INOUT) :: deriv
    INTEGER :: k, l, tls(3)
    CHARACTER(4) :: tl_suff
    INTEGER, PARAMETER :: grids(5) = [GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_EDGE, &
      & GRID_UNSTRUCTURED_VERT, GRID_LONLAT, GRID_ZONAL]

    DO k = 1, 5 ! scan for simple (instant) variable
      vl_elem => vlr_find(vname, opt_patch_id=p_onl%dom, &
        & opt_hgrid=grids(k), opt_list=src_list, opt_cs=.FALSE., opt_output=.TRUE.)
      IF (ASSOCIATED(vl_elem)) THEN
        deriv%tls(1) = -1
        deriv%src(1)%p => vl_elem
        RETURN
      END IF
    END DO
    l = 0
    tls = [nold(1), nnow(1), nnew(1)]
    DO k = 1, 3 ! scan for time-levels of variable
      WRITE(tl_suff, '(a3,i1)') TIMELEVEL_SUFFIX, tls(k)
      vl_elem => vlr_find(vname//tl_suff, &
        & opt_patch_id=p_onl%dom, opt_list=src_list, opt_output=.TRUE.)
      IF (ASSOCIATED(vl_elem)) THEN
        l = l + 1
        deriv%src(l)%p => vl_elem
        deriv%tls(l) = tls(k)
      END IF
    END DO
    IF (l .EQ. 0) CALL finish(routine, 'no source variable: '//TRIM(vname))
  END SUBROUTINE find_src_element

  SUBROUTINE copy_var_to_list(deriv)
    TYPE(t_derivate_var), INTENT(INOUT) :: deriv
    TYPE(t_var_metadata), POINTER :: info

    info => deriv%src(1)%p%info
    CALL add_var(REAL_T, src_list, dname, info%hgrid, info%vgrid, info%cf, &
      & info%grib2, info%used_dimensions(1:info%ndims), vl_elem, &
      & tlev_source=info%tlev_source, isteptype=info%isteptype, &
      & post_op=info%post_op, initval_r=info%initval%rval, &
      & resetval_r=info%resetval%rval, lmiss=info%lmiss, &
      & missval_r=info%missval%rval, &
      & vert_interp=info%vert_interp, hor_interp=info%hor_interp, &
      & in_group=info%in_group, &
      & loutput=.TRUE., lrestart=.FALSE., var_class=info%var_class, &
      & lopenacc = info%lopenacc )
    SELECT CASE(info%hgrid)
    CASE(GRID_UNSTRUCTURED_CELL)
      vl_elem%info%subset = patch_2d%cells%owned
    CASE(GRID_UNSTRUCTURED_EDGE)
      vl_elem%info%subset = patch_2d%edges%owned
    CASE(GRID_UNSTRUCTURED_VERT)
      vl_elem%info%subset = patch_2d%verts%owned
    END SELECT
    vl_elem%info%cf%datatype = &
      & MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    IF ("" == info%cf%short_name) &
      & vl_elem%info%cf%short_name = get_var_name(deriv%src(1)%p%info)
    deriv%dst%p => vl_elem
  END SUBROUTINE copy_var_to_list

  END SUBROUTINE init_op

  SUBROUTINE perform_op(src, dest, funccode, weight, miss, miss_s)
    TYPE(t_var), POINTER, INTENT(IN) :: src
    TYPE(t_var), POINTER, INTENT(INOUT) :: dest
    INTEGER(C_INT), INTENT(IN) :: funccode
    REAL(wp), INTENT(IN), OPTIONAL :: weight, miss
    REAL(sp), INTENT(IN), OPTIONAL :: miss_s
    REAL(wp) :: miss_src, miss_dst, weight_dst
    INTEGER :: ic,si,sj,sk,sl,sm,ei,ej,ek,el,em,sbi,ebi,bk,sbk,ebk
    REAL(wp), POINTER :: sd5d(:,:,:,:,:)
    REAL(sp), POINTER :: ss5d(:,:,:,:,:)
    LOGICAL :: lacc

    IF (PRESENT(weight)) weight_dst = weight
    IF (PRESENT(miss)) miss_dst = miss
    IF (PRESENT(miss)) miss_src = MERGE(miss, REAL(miss_s, wp), ASSOCIATED(src%r_ptr))
    sm = 1; sl = 1; sk = 1; sj = 1; si = 1
    em = SIZE(dest%r_ptr, 5); el = SIZE(dest%r_ptr, 4); ek = SIZE(dest%r_ptr, 3)
    ej = SIZE(dest%r_ptr, 2); ei = SIZE(dest%r_ptr, 1)
    bk = 0; sbk = -1; ebk = HUGE(ebk); sbi = -1; ebi = HUGE(ebi)
    NULLIFY(sd5d, ss5d)
    IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,:,:,:)
    IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,:,:,:)
    IF (funccode .NE. F_MISS .AND. src%info%lcontained) THEN
      ic = src%info%ncontained
      SELECT CASE(src%info%var_ref_pos)
      CASE(1)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(ic:ic,:,:,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(ic:ic,:,:,:,:)
      CASE(2)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,ic:ic,:,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,ic:ic,:,:,:)
      CASE(3)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,ic:ic,:,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,ic:ic,:,:)
      CASE(4)
        IF (ASSOCIATED(src%r_ptr)) sd5d => src%r_ptr(:,:,:,ic:ic,:)
        IF (ASSOCIATED(src%s_ptr)) ss5d => src%s_ptr(:,:,:,ic:ic,:)
      END SELECT
    ELSE IF (funccode .NE. F_MISS) THEN
      IF (dest%info%ndims .EQ. 3) THEN
        bk = 3
        sbk = dest%info%subset%start_block
        ebk = dest%info%subset%end_block
      ELSE IF (dest%info%ndims .EQ. 2 .AND. &
        & GRID_ZONAL .NE. dest%info%hgrid) THEN
        bk = 2
        sbk = dest%info%subset%start_block
        ebk = dest%info%subset%end_block
      END IF
      sbi = dest%info%subset%start_index
      ebi = dest%info%subset%end_index
    END IF
    lacc = dest%info%lopenacc .AND. i_am_accel_node
    CALL perform_op_5d()
  CONTAINS

  SUBROUTINE perform_op_5d()
    INTEGER :: i,j,k,l,m,lblk
    REAL(wp), POINTER :: tmp1(:,:,:,:,:), tmp2(:,:,:,:,:)

    IF (ASSOCIATED(sd5d)) THEN
      tmp1 => sd5d
      
    ELSE IF (funccode .NE. F_MISS .OR. funccode .NE. F_WGT) THEN
        ALLOCATE(tmp1(si:ei,sj:ej,sk:ek,sl:el,sm:em))
        !$ACC ENTER DATA CREATE(tmp1) IF(lacc)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(ss5d, tmp1) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)
      DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        tmp1(i,j,k,l,m) = REAL(ss5d(i,j,k,l,m), wp)
      END DO;      END DO;    END DO;  END DO;END DO
    END IF
    tmp2 => dest%r_ptr
    
    SELECT CASE(funccode)
    CASE(F_ACC)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)
       DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
       tmp2(i,j,k,l,m) = tmp2(i,j,k,l,m) + tmp1(i,j,k,l,m)
       END DO;      END DO;    END DO;  END DO;END DO
    CASE(F_MAX)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        IF (tmp1(i,j,k,l,m) .GT. tmp2(i,j,k,l,m)) tmp2(i,j,k,l,m) = tmp1(i,j,k,l,m)
        END DO;      END DO;    END DO;  END DO;END DO
    CASE(F_MIN)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        IF (tmp1(i,j,k,l,m) .LT. tmp2(i,j,k,l,m)) tmp2(i,j,k,l,m) = tmp1(i,j,k,l,m)
        END DO;      END DO;    END DO;  END DO;END DO
    CASE(F_ACC_SQ)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        tmp2(i,j,k,l,m) = tmp2(i,j,k,l,m) + tmp1(i,j,k,l,m) * tmp1(i,j,k,l,m)
        END DO;      END DO;    END DO;  END DO;END DO
    CASE(F_ASS)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        tmp2(i,j,k,l,m) = tmp1(i,j,k,l,m)
        END DO;      END DO;    END DO;  END DO;END DO               
    CASE(F_ASS_SQ)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        tmp2(i,j,k,l,m) = tmp1(i,j,k,l,m) * tmp1(i,j,k,l,m) 
        END DO;      END DO;    END DO;  END DO;END DO          
    CASE(F_WGT)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        tmp2(i,j,k,l,m) = tmp2(i,j,k,l,m) * weight_dst
        END DO;      END DO;    END DO;  END DO;END DO
    CASE(F_MISS)
!ICON_OMP PARALLEL DO PRIVATE(lblk) COLLAPSE(4)
!$ACC PARALLEL LOOP PRESENT(tmp1, tmp2) COLLAPSE(5) GANG VECTOR ASYNC(1) IF(lacc)       
        DO m = sm, em;  DO l = sl, el;    DO k = sk, ek;      DO j = sj, ej;        DO i = si, ei;          lblk = MERGE(0, MERGE(k, j, bk.EQ.3), bk.EQ.0);          IF (lblk.LT.sbk .OR. lblk.GT.ebk) CYCLE;          IF (lblk.EQ.sbk .AND. i.LT.sbi) CYCLE;          IF (lblk.EQ.ebk .AND. i.GT.ebi) CYCLE
        IF (tmp1(i,j,k,l,m) .EQ. miss_src) tmp2(i,j,k,l,m) = miss_dst
        END DO;      END DO;    END DO;  END DO;END DO
    END SELECT
    IF (ASSOCIATED(ss5d) .AND. ASSOCIATED(tmp1)) THEN
      DEALLOCATE(tmp1)
      !$ACC EXIT DATA DELETE(tmp1) IF(lacc)
    END IF
  END SUBROUTINE perform_op_5d

  END SUBROUTINE perform_op

  !! Execute the accumulation forall internal variables and compute mean values
  !! if the corresponding event is active
  SUBROUTINE update_statistics()
    INTEGER :: iop

    DO iop = 1, nops
      IF (ALLOCATED(ops(iop)%events)) &
        CALL update_op()
    END DO
  CONTAINS

  SUBROUTINE update_op()
    INTEGER :: tl, iv, ie, it, ne, nv
    INTEGER, POINTER :: ct
    TYPE(t_var), POINTER :: src, dst
    TYPE(t_derivate_event), POINTER :: ederiv
    TYPE(datetime) :: mtime_date
    LOGICAL :: isactive

    ne = SIZE(ops(iop)%events)
    mtime_date = time_config%tc_current_date
    DO ie = 1, ne
      ederiv => ops(iop)%events(ie)%a
      isactive = LOGICAL(isCurrentEventActive(ederiv%mtime_event, mtime_date))
      nv = 0
      IF (ALLOCATED(ederiv%vars)) nv = SIZE(ederiv%vars)
      DO iv = 1, nv
        dst => ederiv%vars(iv)%a%dst%p
        ct => ederiv%vars(iv)%a%counter
        src => ederiv%vars(iv)%a%src(1)%p
        IF (ederiv%vars(iv)%a%tls(1) .NE. -1) THEN
          tl = metainfo_get_timelevel(dst%info, dst%info%dom)
          it = MAXLOC(MERGE(1, 0, tl .EQ. ederiv%vars(iv)%a%tls(:)), 1)
          src => ederiv%vars(iv)%a%src(it)%p
        END IF
        IF (ct .EQ. 0) THEN ! initial assignment
          CALL perform_op(src, dst, MERGE(F_ASS_SQ, F_ASS, opcodes(iop) .EQ. O_MEAN_SQ))
        ELSE ! actual update
          CALL perform_op(src, dst, opfuncs(iop))
        END IF
        ct = ct + 1
        IF (isactive) THEN ! output step, so weighting is applied this time for time mean operators
          IF ((O_MEAN .EQ. opcodes(iop) .OR. O_MEAN_SQ .EQ. opcodes(iop)) .AND. ct .GT. 1) &
            & CALL perform_op(src, dst, F_WGT, weight=(1._wp / REAL(ct, wp)))
          IF (dst%info%lmiss) & ! (re)set missval where applicable
            & CALL perform_op(src, dst, F_MISS, miss=dst%info%missval%rval, &
                &             miss_s=src%info%missval%sval)
          ct = 0
        END IF
      END DO
    END DO
  END SUBROUTINE update_op

  END SUBROUTINE update_statistics

END MODULE mo_derived_variable_handling
