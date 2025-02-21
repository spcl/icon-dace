! contains lhs-type for use in solver-backends
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








MODULE mo_ocean_solve_lhs

  USE mo_exception, ONLY: finish
  USE mo_kind, ONLY: sp, wp, i8
  USE mo_mpi, ONLY: p_comm_work, p_max, p_pe_work
  USE mo_timer, ONLY: timer_start, timer_stop, new_timer
  USE mo_ocean_solve_aux, ONLY: t_ocean_solve_parm
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_trivial_transfer, ONLY: t_trivial_transfer
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_nml, ONLY: l_lhs_direct
  USE mo_run_config, ONLY: ltimer
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lhs

  CHARACTER(LEN=*), PARAMETER :: module_name = "mo_ocean_solve_lhs"

  TYPE :: t_lhs
    PRIVATE
    INTEGER :: nblk_loc, nblk_a_loc, nidx_loc, nidx_e_loc, &
      & nnzero_loc, nnzero_cal, nextra, nnzs, timer_upd, &
      & nindep_grp, timer_init, grp_nelem(64), grp_nelem_max
    INTEGER, PUBLIC :: timer
    LOGICAL :: is_const, have_sp
    CLASS(t_lhs_agen), POINTER :: agen => NULL()
    CLASS(t_transfer), POINTER :: trans => NULL()
    LOGICAL, PUBLIC :: is_init = .false.
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: x_t, ax_t, dcoef_c_wp
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:,:) :: coef_c_sp
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) :: coef_l_wp, coef_c_wp
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:) :: dcoef_c_sp
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: blk_loc, idx_loc, blk_cal, &
      & idx_cal, grp_smap_blk, grp_smap_idx, grp_smap_inz
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: grp_map_idx, grp_map_blk, &
      & inz_t, dnz_cal







  CONTAINS
! for actual left-hand-sides
    PROCEDURE, PRIVATE :: apply_wp => lhs_apply_wp
    PROCEDURE, PRIVATE :: apply_sp => lhs_apply_sp
    GENERIC, PUBLIC :: apply => apply_wp, apply_sp
    PROCEDURE, PRIVATE :: apply_noaii_wp => lhs_apply_noaii_wp
    PROCEDURE, PRIVATE :: apply_noaii_sp => lhs_apply_noaii_sp
    GENERIC, PUBLIC :: apply_noaii => apply_noaii_wp, apply_noaii_sp
    PROCEDURE, PRIVATE :: get_invaii_wp => lhs_get_invaii_wp
    PROCEDURE, PRIVATE :: get_invaii_sp => lhs_get_invaii_sp
    GENERIC, PUBLIC :: get_invaii => get_invaii_wp, get_invaii_sp
! for transfer from local to calc fields (maybe on different task)
    PROCEDURE, PUBLIC :: update => lhs_update
! service routines
    PROCEDURE, PUBLIC :: construct => lhs_construct
    PROCEDURE, PUBLIC :: dump_matrix => lhs_dump_matrix
    PROCEDURE, PUBLIC :: dump_ax => lhs_dump_ax
    PROCEDURE, PRIVATE :: find_nnzero => lhs_find_nnzero
    PROCEDURE, PRIVATE :: create_matrix_init => lhs_create_matrix_init
    PROCEDURE, PRIVATE :: create_matrix_redo => lhs_create_matrix_redo
    PROCEDURE, PRIVATE :: create_matrix => lhs_create_matrix
    PROCEDURE, PRIVATE :: doit_wp => lhs_doit_wp
    PROCEDURE, PRIVATE :: doit_sp => lhs_doit_sp
    GENERIC, PRIVATE :: doit => doit_wp, doit_sp
    PROCEDURE, PRIVATE :: noaii_doit_wp => lhs_noaii_doit_wp
    PROCEDURE, PRIVATE :: noaii_doit_sp => lhs_noaii_doit_sp
    GENERIC, PRIVATE :: noaii_doit => noaii_doit_wp, noaii_doit_sp
  END TYPE t_lhs

CONTAINS

! interface routine to compute D^(-1) of lhs-matrix
  SUBROUTINE lhs_get_invaii_wp(this, inv_aii, from_sp)
    CLASS(t_lhs), INTENT(INOUT), TARGET :: this
    REAL(KIND=wp), DIMENSION(:,:), INTENT(OUT), POINTER :: inv_aii
    INTEGER, INTENT(IN), OPTIONAL :: from_sp
    INTEGER :: nidx, inz, iblk, iidx, jidx, jblk
    LOGICAL :: gotsome(this%trans%nidx)
    CHARACTER(LEN=*), PARAMETER :: routine = module_name // &
      & 'lhs_get_invaii_wp'

    IF (l_lhs_direct) &
       CALL finish(routine, "Jacobi-Preconditioner not implemented!")
!if first call
    IF (.NOT.ALLOCATED(this%dcoef_c_wp)) THEN
! alloc vector for diagonal elements
      ALLOCATE(this%dcoef_c_wp(this%trans%nidx, this%trans%nblk), &
        & this%dnz_cal(this%trans%nidx, this%trans%nblk))
! find diagonal elements in lhs-matrix
!ICON_OMP PARALLEL DO PRIVATE(iidx, nidx, gotsome, inz)
      DO iblk = 1, this%trans%nblk
        this%dcoef_c_wp(:, iblk) = 0._wp
        this%dnz_cal(:, iblk) = 1
        nidx = MERGE(this%trans%nidx, this%trans%nidx_e, &
          & iblk .NE. this%trans%nblk)
        gotsome(:) = .false.
        DO inz = 1, this%nnzero_cal
          IF (.NOT.ANY(.NOT.gotsome(:))) CYCLE
          DO iidx = 1, nidx
            IF (gotsome(iidx)) CYCLE
            IF (this%blk_cal(iidx, iblk, inz) .EQ. iblk) THEN
              IF (this%idx_cal(iidx, iblk, inz) .EQ. iidx) THEN
! check, if diagonal element is non-zero, (if not heal it)
                IF (this%coef_c_wp(iidx, iblk, inz) .NE. 0.0_wp) THEN
                  this%dcoef_c_wp(iidx, iblk) = 1._wp / &
                    & this%coef_c_wp(iidx, iblk, inz)
                ELSE
                  this%dcoef_c_wp(iidx, iblk) = 1._wp
                END IF
                this%dnz_cal(iidx, iblk) = inz
                gotsome(iidx) = .true.
              END IF
            END IF
          END DO
        END DO
      END DO
!ICON_OMP END PARALLEL DO
    ELSE IF ((this%have_sp .AND. PRESENT(from_sp)) .OR. .NOT.this%have_sp) THEN
! else : indices of diagonal elements are already known
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iidx, inz, jidx, jblk)
      DO iblk = 1, this%trans%nblk
        DO iidx = 1, this%trans%nidx
          inz = this%dnz_cal(iidx, iblk)
          jidx = this%idx_cal(iidx, iblk, inz)
          jblk = this%blk_cal(iidx, iblk, inz)
! check, if diagonal element is non-zero, (if not heal it)
          IF(this%coef_c_wp(jidx, jblk, inz) .NE. 0._wp) THEN
            this%dcoef_c_wp(iidx, iblk) = 1._wp / &
              & this%coef_c_wp(jidx, jblk, inz)
          ELSE
            this%dcoef_c_wp(iidx, iblk) = 1._wp
          END IF
        END DO
      END DO
!ICON_OMP END PARALLEL DO
    END IF
    inv_aii => this%dcoef_c_wp
  END SUBROUTINE lhs_get_invaii_wp

! as lhs_get_invaii_wp, but sp-variant
  SUBROUTINE lhs_get_invaii_sp(this, inv_aii)
    CLASS(t_lhs), INTENT(INOUT), TARGET :: this
    REAL(KIND=sp), DIMENSION(:,:), INTENT(OUT), POINTER :: inv_aii
    REAL(KIND=wp), DIMENSION(:,:), POINTER :: dummy
    INTEGER :: iblk
    CHARACTER(LEN=*), PARAMETER :: routine = module_name // &
      & 'lhs_get_invaii_sp'

    IF (l_lhs_direct) &
       CALL finish(routine, "Jacobi-Preconditioner not implemented!")
    IF (.NOT.ALLOCATED(this%dcoef_c_sp)) &
      & ALLOCATE(this%dcoef_c_sp(this%trans%nidx, this%trans%nblk), &
        & this%dnz_cal(this%trans%nidx, this%trans%nblk))
    CALL this%get_invaii(dummy, 1)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, this%trans%nblk
      this%dcoef_c_sp(:, iblk) = REAL(dummy(:, iblk), sp)
    END DO
!ICON_OMP END PARALLEL DO
    inv_aii => this%dcoef_c_sp
  END SUBROUTINE lhs_get_invaii_sp

! dump matrix to text-file
  SUBROUTINE lhs_dump_matrix(this, id, prefix, lprecon)
    CLASS(t_lhs), INTENT(IN), TARGET :: this
    INTEGER, INTENT(IN) :: id
    LOGICAL, INTENT(IN) :: lprecon
    CHARACTER(LEN=*), INTENT(IN) :: prefix
    CHARACTER(LEN=128) :: fileName
    INTEGER :: inz, iblk, iidx
    INTEGER, PARAMETER :: fileNo = 501
    CHARACTER(LEN=*), PARAMETER :: routine = module_name // &
      & "::lhs_dump_matrix()"

    IF (lprecon)  & !THEN
      & CALL finish(routine, &
        & "cannot dump preconditioner matrix if no precon present")
    WRITE(fileName, "(A,I0.4,A,i0.4,a)") &
      & TRIM(prefix)//"_",id,"_",p_pe_work,".txt"
    OPEN(fileNo, FILE=TRIM(fileName), STATUS='new')
    DO inz = 1, SIZE(this%coef_l_wp, 3)
      DO iblk = 1, this%nblk_loc
        DO iidx = 1, this%nidx_loc
          WRITE(fileNo, "(2(a,2(i8.8,a)),es12.5)") &
            & "(" , iidx, ":" ,iblk,"), ", &
            & "(" , this%idx_loc(iidx, iblk, inz), &
            & ":" , this%blk_loc(iidx, iblk, inz), ")", &
            & this%coef_l_wp(iidx, iblk, inz)
        END DO
      END DO
    END DO
    CLOSE(fileNo)
  END SUBROUTINE lhs_dump_matrix

! dump result of vector -- lhs-matrix multiplication to text-file
  SUBROUTINE lhs_dump_ax(this, prefix, ax)
    CLASS(t_lhs), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: prefix
    REAL(KIND=wp), INTENT(IN) :: ax(:,:)
    CHARACTER(LEN=128) :: fileName
    INTEGER :: iblk, iidx
    INTEGER, PARAMETER :: fileNo = 501

    CALL this%update()
    write(fileName, "(A,i0.4,a)") &
      & TRIM(prefix)//"_",p_pe_work,".txt"
    open (fileNo, FILE=TRIM(fileName), STATUS='new')
    DO iblk = 1, this%nblk_a_loc
      DO iidx = 1, this%nidx_loc
        WRITE(fileNo, "(a,2(i8.8,a),es12.5)") &
          & "(" , iidx, ":" ,iblk,"), ", &
          & ax(iidx, iblk)
      END DO
    END DO
    CLOSE(fileNo)
  END SUBROUTINE lhs_dump_ax

! interface routine to update lhs-matrix
  SUBROUTINE lhs_update(this)
    CLASS(t_lhs), INTENT(INOUT) :: this
    CHARACTER(LEN=*),PARAMETER :: routine = module_name// &
      & "::lhs_update()"

    IF (ltimer) CALL timer_start(this%timer_upd)
    IF (.NOT.this%is_init) &
      &CALL finish(routine, "t_lhs was not initiaized-...!")
    IF (.NOT.this%is_const .AND. .NOT.l_lhs_direct) THEN
! re-compute coeffs
      CALL this%create_matrix()
! transfer from worker- to solve-PEs
      CALL this%trans%into(this%coef_l_wp, this%coef_c_wp, 2)
      IF (this%have_sp .AND. this%trans%is_solver_pe) &
        & this%coef_c_sp(:,:,:) = REAL(this%coef_c_wp(:,:,:), sp)
    END IF
    IF (ltimer) CALL timer_stop(this%timer_upd)
  END SUBROUTINE lhs_update

! interface to init lhs-object using the provided t_agen - object (i.e. matrix generator)
  SUBROUTINE lhs_construct(this, have_sp, par, agen, trans)
    CLASS(t_lhs), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN) :: have_sp
    TYPE(t_ocean_solve_parm), INTENT(IN) :: par
    CLASS(t_lhs_agen), TARGET, INTENT(IN) :: agen
    CLASS(t_transfer), TARGET, INTENT(IN) :: trans
    CHARACTER(LEN=*), PARAMETER :: routine = module_name//&
      & "::lhs_construct()"

    IF (ltimer) THEN
      this%timer_init = new_timer("lhs init")
      CALL timer_start(this%timer_init)
    END IF
    IF (agen%is_init) THEN
      this%agen => agen
    ELSE
      CALL finish(routine, &
        & "matrix generating type/func has to be initialized...")
    END IF
    IF (trans%is_init) THEN
      this%trans => trans
      IF (l_lhs_direct) THEN
        SELECT TYPE (trans)
        CLASS IS (t_trivial_transfer)
        CLASS DEFAULT
          CALL finish(routine, &
            & "direct use of t_agen in lhs only for t_trivial_transfer type!")
        END SELECT
      END IF
    ELSE
      CALL finish(routine, &
        & "transfer type has to be initialized...")
    END IF
    IF (par%nidx .GT. 0) THEN
      this%nidx_loc = par%nidx
      this%nidx_e_loc = par%nidx_e
    ELSE
      CALL finish(routine, &
        & "no valid nidx a.k.a. nproma was provided")
    END IF
    IF (par%nblk .GT. 0) THEN
      this%nblk_loc = par%nblk
      this%nblk_a_loc = par%nblk_a
    ELSE
      this%nblk_loc = 0
      this%nblk_a_loc = 1
    END IF
    this%have_sp = have_sp
    this%is_const = this%agen%is_const
    IF (.NOT.l_lhs_direct) THEN
      IF (.NOT.agen%use_shortcut) THEN
        CALL this%find_nnzero()
        ALLOCATE(this%coef_l_wp(this%nidx_loc, this%nblk_loc, this%nnzero_cal))
        ALLOCATE(this%idx_loc(this%nidx_loc, this%nblk_loc, this%nnzero_cal))
        ALLOCATE(this%blk_loc(this%nidx_loc, this%nblk_loc, this%nnzero_cal))
      END IF
! create initial lhs-matrix
      CALL this%create_matrix(.true.)
      CALL trans%into_once(this%coef_l_wp, this%coef_c_wp, 4)
      IF (have_sp .AND. this%trans%is_solver_pe) THEN
        ALLOCATE(this%coef_c_sp(trans%nidx, trans%nblk, this%nnzero_cal))
        this%coef_c_sp(:,:,:) = REAL(this%coef_c_wp(:,:,:), sp)
      END IF
      CALL trans%into_once(this%idx_loc, this%blk_loc, &
        & this%idx_cal, this%blk_cal, 4)
      IF (this%is_const) &
        & DEALLOCATE(this%x_t, this%ax_t, this%coef_l_wp, this%idx_loc, this%blk_loc)
    END IF
    IF (ltimer) THEN
      CALL timer_stop(this%timer_init)
      this%timer_upd = new_timer("lhs update")
      this%timer = new_timer("lhs apply")
    END IF
    this%is_init = .true.
  END SUBROUTINE lhs_construct

! backend routine to update or create the lhs-matrix
  SUBROUTINE lhs_create_matrix(this, on_init_in)
    CLASS(t_lhs), INTENT(INOUT) :: this
    LOGICAL, INTENT(IN), OPTIONAL :: on_init_in
    LOGICAL :: on_init

    on_init = .false.
    IF (PRESENT(on_init_in)) on_init = on_init_in
    IF (on_init) THEN
      CALL this%create_matrix_init()
    ELSE
      CALL this%create_matrix_redo()
    END IF
  END SUBROUTINE lhs_create_matrix

! backend routine to update the lhs matrix from a provided "t_agen" (the generator for matrix 'A')
  SUBROUTINE lhs_create_matrix_redo(this)
    CLASS(t_lhs), INTENT(INOUT) :: this
    INTEGER :: inz, igrp, ielem

    IF (this%agen%use_shortcut) THEN
      CALL this%agen%matrix_shortcut(this%idx_loc, this%blk_loc, this%coef_l_wp)
      RETURN
    END IF
! iterate over groups (see lhs_create_matrix_init)
    DO igrp = 1, this%nindep_grp
      DO ielem = 1, this%grp_nelem(igrp)
        this%x_t(this%grp_map_idx(ielem, igrp), &
          & this%grp_map_blk(ielem, igrp)) = 1._wp
      END DO
      CALL this%agen%apply(this%x_t, this%ax_t)
! update coeffs from this group
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(inz)
      DO ielem = 1, this%grp_nelem(igrp)
        DO inz = 1, this%nnzero_loc
          IF (this%grp_smap_idx(inz, ielem, igrp) .GT. 0) THEN
            this%coef_l_wp(this%grp_smap_idx(inz, ielem, igrp), &
              & this%grp_smap_blk(inz, ielem, igrp), &
              & this%grp_smap_inz(inz, ielem, igrp)) = &
              & this%ax_t(this%grp_smap_idx(inz, ielem, igrp), &
              & this%grp_smap_blk(inz, ielem, igrp))
            this%ax_t(this%grp_smap_idx(inz, ielem, igrp), &
              & this%grp_smap_blk(inz, ielem, igrp)) = 0._wp
          END IF
        END DO
        this%x_t(this%grp_map_idx(ielem, igrp), &
          & this%grp_map_blk(ielem, igrp)) = 0._wp
      END DO
!ICON_OMP END PARALLEL DO
    END DO
! in case another PE has more groups
    DO igrp = 1, this%nextra
      CALL this%agen%apply(this%x_t, this%ax_t)
    END DO
  END SUBROUTINE lhs_create_matrix_redo

! backend routine for creating the lhs matrix from a provided "t_agen" (the generator for matrix 'A')
  SUBROUTINE lhs_create_matrix_init(this)
    CLASS(t_lhs), INTENT(INOUT) :: this
    INTEGER :: iidx, iblk, inz, nidx, jidx, jblk, igrp, jnz, ielem
    INTEGER(KIND=i8) :: grp_codom(this%nidx_loc, this%nblk_loc)
    INTEGER :: grp_map(this%nidx_loc, this%nblk_a_loc)
    INTEGER :: image_idx(64), image_blk(64), elem_ct(64)
    INTEGER(KIND=i8) :: grp_id_image
    LOGICAL :: next
    INTEGER :: ngid, jgid, igid, sgid, sgid_blk, sgid_idx
    INTEGER, ALLOCATABLE, DIMENSION(:) :: gid_blk, gid_idx, gid_list
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: idx_s_loc, blk_s_loc, inz_s_loc

    IF (this%agen%use_shortcut) THEN
      ALLOCATE(this%coef_l_wp(this%nidx_loc, this%nblk_loc, 1))
      CALL this%agen%matrix_shortcut(this%idx_loc, this%blk_loc, this%coef_l_wp)
      this%nnzero_cal = SIZE(this%idx_loc, 3)
      RETURN
    END IF
    ALLOCATE(idx_s_loc(this%nnzero_loc, this%nidx_loc, this%nblk_a_loc), &
      & blk_s_loc(this%nnzero_loc, this%nidx_loc, this%nblk_a_loc), &
      & inz_s_loc(this%nnzero_loc, this%nidx_loc, this%nblk_a_loc))
! allocate and init necessary arrays
    IF (.NOT.ALLOCATED(this%x_t)) THEN
      ALLOCATE(this%x_t(this%nidx_loc, this%nblk_a_loc), &
        this%ax_t(this%nidx_loc, this%nblk_a_loc))
!ICON_OMP PARALLEL WORKSHARE
      this%x_t(:,:) = 0._wp
      this%ax_t(:,:) = 0._wp
!ICON_OMP END PARALLEL WORKSHARE
    END IF
    IF (.NOT.ALLOCATED(this%inz_t)) &
      ALLOCATE(this%inz_t(this%nidx_loc, this%nblk_loc))
!ICON_OMP PARALLEL WORKSHARE
    this%idx_loc(:,:,:) = 1
    this%blk_loc(:,:,:) = 1
    this%inz_t(:,:) = 0
    this%coef_l_wp(:,:,:) = 0._wp
    image_idx(:) = 0
    image_blk(:) = 0
    grp_codom(:,:) = 0
    this%grp_nelem(:) = 0
    grp_map(:,:) = 127
!ICON_OMP END PARALLEL WORKSHARE
! set a valid dummy
    IF (this%nblk_loc > 0) THEN
      this%idx_loc(1,1,:) = 2
      this%blk_loc(1,1,:) = 2
    END IF
    this%nindep_grp = 0
    ngid = this%trans%ngid_a_l
! ordering by global indices is necessary to get bit-identical solutions, even
! if varying the decomposition for whatever reason
! get the global indices of the elements (from trans) -  and block them
    ALLOCATE(gid_blk(ngid), gid_idx(ngid), gid_list(ngid))
    DO jgid = 1, ngid
      gid_blk(jgid) = (jgid + this%nidx_loc - 1) / this%nidx_loc
      gid_idx(jgid) = jgid - (gid_blk(jgid)-1) * this%nidx_loc
      gid_list(jgid) = this%trans%globalID_loc(gid_idx(jgid), gid_blk(jgid))
    END DO
! sort by global IDs
    DO jgid = 2, ngid
      sgid =  gid_list(jgid)
      IF (sgid .LT. 0) EXIT
      igid = jgid
      sgid_blk = gid_blk(jgid)
      sgid_idx = gid_idx(jgid)
      DO WHILE (igid .GT. 1)
        IF (gid_list(igid-1) .GT. sgid) THEN
          gid_list(igid) = gid_list(igid-1)
          gid_blk(igid) = gid_blk(igid-1)
          gid_idx(igid) = gid_idx(igid-1)
          igid = igid - 1
        ELSE
          EXIT
        END IF
      END DO
      gid_list(igid) = sgid
      gid_blk(igid) = sgid_blk
      gid_idx(igid) = sgid_idx
    END DO
! scan x-vector in order of ascending global IDs
    ngid = jgid - 1
    this%nextra = p_max(ngid, comm=p_comm_work) - ngid
    DEALLOCATE(gid_list)
    DO igid = 1, ngid
      iidx = gid_idx(igid)
      iblk = gid_blk(igid)
      idx_s_loc(:, iidx, iblk) = 0
      blk_s_loc(:, iidx, iblk) = 0
      inz_s_loc(:, iidx, iblk) = 0
      this%x_t(iidx, iblk) = 1._wp
      CALL this%agen%apply(this%x_t, this%ax_t)
      inz = 0
! find non-zero elements in ax = Ax and store their location and value (twice)
      DO jblk = 1, this%nblk_loc
        nidx = MERGE(this%nidx_loc, this%nidx_e_loc, jblk .NE. this%nblk_loc)
        DO jidx = 1, nidx
          IF (this%ax_t(jidx, jblk) .NE. 0._wp) THEN
            inz = inz + 1
            IF (inz .GT. 63) THEN
              CALL finish("lhs_create_matrix_init", &
                & "stencils with more than 63 entries are not supported")
            END IF
            this%inz_t(jidx, jblk) = this%inz_t(jidx, jblk) + 1
            inz_s_loc(inz, iidx, iblk) = this%inz_t(jidx, jblk)
            idx_s_loc(inz, iidx, iblk) = jidx
            blk_s_loc(inz, iidx, iblk) = jblk
            this%idx_loc(jidx, jblk, this%inz_t(jidx, jblk)) = iidx
            this%blk_loc(jidx, jblk, this%inz_t(jidx, jblk)) = iblk
            this%coef_l_wp(jidx, jblk, this%inz_t(jidx, jblk)) = &
              & this%ax_t(jidx, jblk)
            this%ax_t(jidx, jblk) = 0._wp
            image_idx(inz) = jidx
            image_blk(inz) = jblk
          END IF
        END DO
      END DO
      this%x_t(iidx, iblk) = 0._wp
      grp_id_image = 1_i8
      next = .true.
! grouping speeds up the matrix update done before every solve (if lhs-generator is not flagged constant)
! because agen%apply has to be called for each group only once, instead of for every single element
! find group of elements of x, such that the image under the lhs-matrix does not interfere with that of the other group members
      DO igrp = 1, this%nindep_grp
        next = .false.
        DO jnz = 1, inz
          IF (IAND(grp_id_image, &
            & grp_codom(image_idx(jnz), image_blk(jnz))) &
            & .NE. 0) THEN
            next = .true.
            EXIT
          END IF
        END DO
        IF (.NOT.next) EXIT
        grp_id_image = grp_id_image * 2_i8
      END DO
! if no group was found, open up a new group
      IF (next) THEN
        this%nindep_grp = this%nindep_grp + 1
        IF (this%nindep_grp > SIZE(this%grp_nelem)) &
          CALL finish("lhs_create_matrix_init", &
            "Too many independent groups. Increase t_lhs::grp_nelem")
        igrp = this%nindep_grp
      END IF
! store assignment of x-element to group
      grp_map(iidx, iblk) = igrp
      this%grp_nelem(igrp) = this%grp_nelem(igrp) + 1
      DO jnz = 1, inz
        grp_codom(image_idx(jnz), image_blk(jnz)) = &
          & grp_codom(image_idx(jnz), image_blk(jnz)) &
          & + grp_id_image
      END DO
      image_idx(1:inz) = 0
      image_blk(1:inz) = 0
    END DO
    DEALLOCATE(this%inz_t)
    DEALLOCATE(gid_blk, gid_idx)
    this%grp_nelem_max = MAXVAL(this%grp_nelem(:))
    ALLOCATE(this%grp_map_idx(this%grp_nelem_max, this%nindep_grp), &
      & this%grp_map_blk(this%grp_nelem_max, this%nindep_grp))
    elem_ct(:) = 0
! reorder indices info to grouped ordering
    DO iblk = 1, this%nblk_a_loc
      DO iidx = 1, this%nidx_loc
        IF (grp_map(iidx, iblk) .EQ. 127) CYCLE
        elem_ct(grp_map(iidx, iblk)) = &
          & elem_ct(grp_map(iidx, iblk)) + 1
        this%grp_map_idx(elem_ct(grp_map(iidx, iblk)), &
          & grp_map(iidx, iblk)) = iidx
        this%grp_map_blk(elem_ct(grp_map(iidx, iblk)), &
          & grp_map(iidx, iblk)) = iblk
      END DO
    END DO
! allocate arrays to hold grouped matrix info
    ALLOCATE( this%grp_smap_idx(this%nnzero_cal, &
          & this%grp_nelem_max, this%nindep_grp), &
      & this%grp_smap_blk(this%nnzero_cal, &
          & this%grp_nelem_max, this%nindep_grp), &
      & this%grp_smap_inz(this%nnzero_cal, &
          & this%grp_nelem_max, this%nindep_grp))
! store group ordered matrix info (to be used, when updating the matrix coeffs)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ielem, inz, iidx, iblk)
    DO igrp = 1, this%nindep_grp
      DO ielem = 1, this%grp_nelem(igrp)
        DO inz = 1, this%nnzero_loc
          iidx = this%grp_map_idx(ielem, igrp)
          iblk = this%grp_map_blk(ielem, igrp)
          this%grp_smap_idx(inz, ielem, igrp) = &
            & idx_s_loc(inz, iidx, iblk)
          this%grp_smap_blk(inz, ielem, igrp) = &
            & blk_s_loc(inz, iidx, iblk)
          this%grp_smap_inz(inz, ielem, igrp) = &
            & inz_s_loc(inz, iidx, iblk)
        END DO
      END DO
    END DO
!ICON_OMP END PARALLEL DO
    DEALLOCATE(idx_s_loc, blk_s_loc, inz_s_loc)
! in case another PE own more active elements
    DO iidx = 1, this%nextra
      CALL this%agen%apply(this%x_t, this%ax_t)
    END DO
    this%nextra = p_max(this%nindep_grp, comm=p_comm_work) - this%nindep_grp
  END SUBROUTINE lhs_create_matrix_init

! backend routine applying lhs-matrix

  PURE &

  SUBROUTINE lhs_doit_wp(this, x, ax, a , b, i, lacc)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: a
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: i, b
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    REAL(KIND=wp) :: x_t(this%trans%nidx)
    INTEGER :: iidx, iblk, inz
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(this, this%trans, x, a, b, i) &
    !$ACC   COPY(ax) &
    !$ACC   CREATE(x_t) IF(lzacc)

!ICON_OMP PARALLEL
! apply ax(i) = sum(A(i,j)*x(j))
!ICON_OMP DO PRIVATE(inz, iidx, x_t)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO iblk = 1, this%trans%nblk
      ax(:, iblk) = 0._wp
      DO inz = 1, SIZE(a, 3)
        FORALL(iidx = 1:this%trans%nidx) x_t(iidx) = &
          & x(i(iidx, iblk, inz), b(iidx, iblk, inz))
        ax(:, iblk) = ax(:, iblk) + x_t(:) * a(:, iblk, inz)
      END DO
    END DO
    !$ACC END KERNELS
!ICON_OMP END DO NOWAIT
! zero all non-active elements
!ICON_OMP DO SCHEDULE(DYNAMIC)
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO iblk = this%trans%nblk + 1, SIZE(ax, 2)
      ax(:, iblk) = 0.0_wp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
    !$ACC END DATA
  END SUBROUTINE lhs_doit_wp

! backend routine applying lhs-matrix, but omitting diagonal elements
  PURE SUBROUTINE lhs_noaii_doit_wp(this, x, ax, a , b, i)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: a
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: i, b
    INTEGER :: iidx, iblk, inz

!ICON_OMP PARALLEL
! apply ax(i) = sum(A(i,j)*x(j))
!ICON_OMP DO PRIVATE(inz, iidx)
    DO iblk = 1, this%trans%nblk
      ax(:, iblk) = 0._wp
      DO inz = 1, SIZE(a, 3)
        DO iidx = 1, this%trans%nidx
! exclude diagonal elements
          IF (i(iidx, iblk, inz) .NE. iidx) THEN
            IF (b(iidx, iblk, inz) .NE. iblk) THEN
              ax(iidx, iblk) = ax(iidx, iblk) + &
                & x(i(iidx, iblk, inz), b(iidx, iblk, inz)) &
                & * a(iidx, iblk, inz)
            END IF
          END IF
        END DO
      END DO
    END DO
!ICON_OMP END DO NOWAIT
! zero all non-active elements
!ICON_OMP DO SCHEDULE(DYNAMIC, 1)
    DO iblk = this%trans%nblk + 1, SIZE(ax, 2)
      ax(:, iblk) = 0.0_wp
    END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
  END SUBROUTINE lhs_noaii_doit_wp

! interface for solvers, applying lhs-matrix
  SUBROUTINE lhs_apply_wp(this, x, ax, opt_direct, lacc)
    CLASS(t_lhs), INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    LOGICAL, INTENT(IN), OPTIONAL :: opt_direct
    LOGICAL, INTENT(in), OPTIONAL :: lacc

    LOGICAL :: lzacc
    LOGICAL :: l_direct
    CHARACTER(LEN=*),PARAMETER :: routine = module_name//":lhs_apply_wp()"

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (.NOT.this%is_init) CALL finish(routine, "t_lhs was not initiaized-...!")
    l_direct = l_lhs_direct
    IF (PRESENT(opt_direct)) l_direct = opt_direct
    IF (ltimer) CALL timer_start(this%timer)
    IF (.NOT.l_direct) THEN
      CALL this%doit(x, ax, this%coef_c_wp, &
       & this%blk_cal, this%idx_cal, lacc=lzacc)
    ELSE
      CALL this%agen%apply(x, ax, lacc=lzacc)
    END IF
    IF (ltimer) CALL timer_stop(this%timer)
  END SUBROUTINE lhs_apply_wp

! interface for solvers, applying lhs-matrix, but omitting diagonal elements
  SUBROUTINE lhs_apply_noaii_wp(this, x, ax)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    CHARACTER(LEN=*),PARAMETER :: routine = module_name//":lhs_apply_noaii_wp()"

    IF (.NOT.this%is_init) CALL finish(routine, "t_lhs was not initiaized-...!")
    IF (ltimer) CALL timer_start(this%timer)
    CALL this%noaii_doit(x, ax, this%coef_c_wp, &
      & this%blk_cal, this%idx_cal)
    IF (ltimer) CALL timer_stop(this%timer)
  END SUBROUTINE lhs_apply_noaii_wp

! sp-variant of lhs_doit_wp

  PURE &

  SUBROUTINE lhs_doit_sp(this, x, ax, a, b, i, lacc)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=sp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    REAL(KIND=sp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: a
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: i, b
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    REAL(KIND=sp) :: x_t(this%trans%nidx)
    INTEGER :: iidx, iblk, inz
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(this, this%trans, x, a, b, i) &
    !$ACC   COPY(ax) &
    !$ACC   CREATE(x_t) IF(lzacc)

!ICON_OMP PARALLEL
!ICON_OMP DO PRIVATE(inz, iidx, x_t)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO iblk = 1, this%trans%nblk
      ax(:, iblk) = 0._wp
      DO inz = 1, SIZE(a, 3)
        FORALL(iidx = 1:this%trans%nidx) x_t(iidx) = &
          & x(i(iidx, iblk, inz), b(iidx, iblk, inz))
        ax(:, iblk) = ax(:, iblk) + x_t(:) * a(:, iblk, inz)
      END DO
    END DO
    !$ACC END KERNELS
!ICON_OMP END DO NOWAIT
!ICON_OMP DO SCHEDULE(DYNAMIC, 1)
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO iblk = this%trans%nblk + 1, SIZE(ax, 2)
      ax(:, iblk) = 0.0_sp
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
    !$ACC END DATA
  END SUBROUTINE lhs_doit_sp

! sp-variant of lhs_apply_wp
  SUBROUTINE lhs_apply_sp(this, x, ax, opt_direct)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=sp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    LOGICAL, INTENT(IN), OPTIONAL :: opt_direct
    LOGICAL :: l_direct
    CHARACTER(LEN=*),PARAMETER :: routine = module_name//":lhs_apply_sp()"

    IF (.NOT.this%is_init) CALL finish(routine, "t_lhs was not initiaized-...!")
    l_direct = l_lhs_direct
    IF (PRESENT(opt_direct)) l_direct = opt_direct
    IF (ltimer) CALL timer_start(this%timer)
    IF (.NOT.l_direct) THEN
      CALL this%doit(x, ax, this%coef_c_sp, &
        & this%blk_cal, this%idx_cal)
    ELSE
       CALL finish(routine, "l_lhs_direct mode not possible with sp-solve")
    END IF
    IF (ltimer) CALL timer_stop(this%timer)
  END SUBROUTINE lhs_apply_sp

! sp-variant of lhs_noaii_doit_wp
  SUBROUTINE lhs_noaii_doit_sp(this, x, ax, a , b, i)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=sp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    REAL(KIND=sp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: a
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: i, b
    INTEGER :: iidx, iblk, inz

!ICON_OMP PARALLEL
!ICON_OMP DO PRIVATE(inz, iidx)
    DO iblk = 1, this%trans%nblk
      ax(:, iblk) = 0._sp
      DO inz = 1, SIZE(a, 3)
        DO iidx = 1, this%trans%nidx
          IF (i(iidx, iblk, inz) .NE. iidx) THEN
            IF (b(iidx, iblk, inz) .NE. iblk) THEN
              ax(iidx, iblk) = ax(iidx, iblk) + &
                & x(i(iidx, iblk, inz), b(iidx, iblk, inz)) &
                & * a(iidx, iblk, inz)
            END IF
          END IF
        END DO
      END DO
    END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP DO SCHEDULE(DYNAMIC, 1)
    DO iblk = this%trans%nblk + 1, SIZE(ax, 2)
      ax(:, iblk) = 0.0_sp
    END DO
!ICON_OMP END DO NOWAIT
!ICON_OMP END PARALLEL
  END SUBROUTINE lhs_noaii_doit_sp

! sp-variant of lhs_apply_noaii_wp
  SUBROUTINE lhs_apply_noaii_sp(this, x, ax)
    CLASS(t_lhs), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: x
    REAL(KIND=sp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: ax
    CHARACTER(LEN=*),PARAMETER :: routine = module_name//":lhs_apply_noaii_sp()"

    IF (.NOT.this%is_init) CALL finish(routine, "t_lhs was not initiaized-...!")
    IF (ltimer) CALL timer_start(this%timer)
    CALL this%noaii_doit(x, ax, this%coef_c_sp, this%blk_cal, this%idx_cal)
    IF (ltimer) CALL timer_stop(this%timer)
  END SUBROUTINE lhs_apply_noaii_sp

! find the max number of non-zero elements in the lhs-matrix
  SUBROUTINE lhs_find_nnzero(this)
    CLASS(t_lhs), INTENT(INOUT) :: this
    INTEGER :: iidx, iblk, inz, jidx, jblk, inz_max, &
      & nelems, nmax_elems, nidx

! alloc temporary arrays, if not done, yet
    IF (.NOT.ALLOCATED(this%x_t)) &
      & ALLOCATE(this%x_t(this%nidx_loc, this%nblk_a_loc), &
        & this%ax_t(this%nidx_loc, this%nblk_a_loc), &
        & this%inz_t(this%nidx_loc, this%nblk_loc))
! init temporary arrays
!ICON_OMP PARALLEL WORKSHARE
    this%x_t(:,:) = 0._wp
    this%ax_t(:,:) = 0._wp
    this%inz_t(:,:) = 0
!ICON_OMP END PARALLEL WORKSHARE
    nelems = this%nblk_a_loc * this%nidx_loc
    nmax_elems = p_max(nelems, comm=p_comm_work)
    this%nextra = nmax_elems - nelems
    inz_max = 0
! scan x-vector
    DO iblk = 1, this%nblk_a_loc
      DO iidx = 1, this%nidx_loc
        this%x_t(iidx, iblk) = 1._wp
        CALL this%agen%apply(this%x_t, this%ax_t)
        inz = 0
! find non-zero elements of Ax and count then
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(jidx, nidx)
        DO jblk = 1, this%nblk_loc
          nidx = MERGE(this%nidx_loc, this%nidx_e_loc, jblk .NE. this%nblk_loc)
          DO jidx = 1, nidx
            IF (this%ax_t(jidx, jblk) .NE. 0._wp) THEN
!ICON_OMP ATOMIC
              inz = inz + 1
!ICON_OMP ATOMIC
              this%inz_t(jidx,jblk) = this%inz_t(jidx,jblk) + 1
              this%ax_t(jidx, jblk) = 0._wp
            END IF
          END DO
        END DO
!ICON_OMP END PARALLEL DO
        inz_max = MAX(inz, inz_max)
        this%x_t(iidx, iblk) = 0._wp
      END DO
    END DO
! in case we have fewer active elements than other solver-PEs
    DO iidx = 1, this%nextra
      CALL this%agen%apply(this%x_t, this%ax_t)
    END DO
    this%nnzero_loc = inz_max
    inz_max = MAXVAL(this%inz_t(:,:))
    this%nnzero_cal = p_max(inz_max, comm=p_comm_work)
  END SUBROUTINE lhs_find_nnzero

END MODULE mo_ocean_solve_lhs
