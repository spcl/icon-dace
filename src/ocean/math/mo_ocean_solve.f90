! contains general interface to the actual solver backends (init, solve, destruct)
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!





MODULE mo_ocean_solve
  !-------------------------------------------------------------------------
  USE mo_exception, ONLY: finish
  USE mo_kind, ONLY: wp
  USE mo_ocean_solve_lhs_type, ONLY: t_lhs_agen
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_trivial_transfer, ONLY: t_trivial_transfer
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_ocean_solve_gmres, ONLY: t_ocean_solve_gmres
  USE mo_ocean_solve_cg, ONLY: t_ocean_solve_cg
  USE mo_ocean_solve_cgj, ONLY: t_ocean_solve_cgj
  USE mo_ocean_solve_bicgStab, ONLY: t_ocean_solve_bicgStab
  USE mo_ocean_solve_legacy_gmres, ONLY: t_ocean_solve_legacy_gmres
  USE mo_ocean_solve_mres, ONLY: t_ocean_solve_mres
  USE mo_ocean_solve_aux, ONLY: solve_gmres, solve_cg, solve_mres, &
    & solve_precon_jac, solve_bcgs, solve_legacy_gmres, &
    & t_ocean_solve_parm
  USE mo_run_config, ONLY: ltimer
  USE mo_timer, ONLY: new_timer, timer_start, timer_stop
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_ocean_solve
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve'

! general solver object type
  TYPE :: t_ocean_solve
    PRIVATE
! abstract backend-type
    CLASS(t_ocean_solve_backend), ALLOCATABLE :: act
! holds local portion of guess->solution
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: x_loc_wp
! holds local portion of rhs
    REAL(KIND=wp), POINTER, DIMENSION(:,:), PUBLIC :: b_loc_wp
! holds residuals of each solver iteration
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:), PUBLIC :: res_loc_wp
    INTEGER :: sol_type, timer, timer_init
! name of actual backend chosen
    CHARACTER(LEN=64), PUBLIC :: sol_type_name
    LOGICAL, PUBLIC :: is_init = .false.
  CONTAINS
! interfaces
    PROCEDURE :: dump_matrix => ocean_solve_dump_matrix
    PROCEDURE :: construct => ocean_solve_construct
    PROCEDURE :: solve => ocean_solve_solve
  END TYPE t_ocean_solve

CONTAINS

! dumps lhs-matrix to text-file
  SUBROUTINE ocean_solve_dump_matrix(this, id, lprecon_in)
    CLASS(t_ocean_solve), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: id
    LOGICAL, INTENT(IN), OPTIONAL :: lprecon_in
    LOGICAL :: lprecon
    CHARACTER(LEN=*), PARAMETER :: routine = this_mod_name // &
      & "::t_ocean_solve::ocean_solve_dump_matrix()"

    lprecon = .false.
    IF (PRESENT(lprecon_in)) lprecon = lprecon_in
    CALL this%act%dump_matrix(id, lprecon)
  END SUBROUTINE ocean_solve_dump_matrix

! init solver object (allocate backend and initialize it)
  SUBROUTINE ocean_solve_construct(this, st, par, par_sp, lhs_agen, trans, lacc)
    CLASS(t_ocean_solve), TARGET, INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: st
    TYPE(t_ocean_solve_parm), INTENT(IN) :: par, par_sp
    CLASS(t_lhs_agen), TARGET, INTENT(IN) :: lhs_agen
    CLASS(t_transfer), TARGET, INTENT(IN) :: trans
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = this_mod_name// &
      & "::t_ocean_solve::ocean_solve_construct()"

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (ALLOCATED(this%act)) CALL finish(routine, "already initialized!")
    IF (ltimer) THEN
      this%timer_init = new_timer("solver init")
      CALL timer_start(this%timer_init)
    END IF
! decide which backend-solver to use
    SELECT CASE(st)
    CASE(solve_gmres)
      this%timer = new_timer("gmres_solve")
      WRITE(this%sol_type_name, "(a)") "GMRES"
      ALLOCATE(t_ocean_solve_gmres :: this%act)
    CASE(solve_cg)
      SELECT CASE(par%pt)
      CASE(solve_precon_jac)
        this%timer = new_timer("cgj_solve")
        WRITE(this%sol_type_name, "(a)") "CG+JAC"
        ALLOCATE(t_ocean_solve_cgj :: this%act)
      CASE DEFAULT
        this%timer = new_timer("cg_solve")
        WRITE(this%sol_type_name, "(a)") "CG"
        ALLOCATE(t_ocean_solve_cg :: this%act)
      END SELECT
    CASE(solve_bcgs)
      this%timer = new_timer("bcgs_solve")
      WRITE(this%sol_type_name, "(a)") "BCGS"
      ALLOCATE(t_ocean_solve_bicgStab :: this%act)
    CASE(solve_legacy_gmres)
      this%timer = new_timer("leg_gmres_solve")
      SELECT TYPE (trans)
      CLASS IS (t_trivial_transfer)
      CLASS DEFAULT
        CALL finish(routine, &
          & "legacy_gmres must use trivial transfer infrastructure")
      END SELECT
      WRITE(this%sol_type_name, "(a)") "GMRES(legacy)"
      ALLOCATE(t_ocean_solve_legacy_gmres :: this%act)
    CASE(solve_mres)
      this%timer = new_timer("mres_solve")
      WRITE(this%sol_type_name, "(a)") "MINRES"
      ALLOCATE(t_ocean_solve_mres :: this%act)
    CASE DEFAULT
      CALL finish(routine, "unrecognized solver_type")
    END SELECT
    this%sol_type = st
! init backend
    CALL this%act%construct(par, par_sp, lhs_agen, trans, lacc=lzacc)
! init arrays / pointers
    NULLIFY(this%b_loc_wp)
    ALLOCATE(this%x_loc_wp(par%nidx, par%nblk_a), this%res_loc_wp(2))
    this%x_loc_wp(:, par%nblk_a) = 0._wp
    this%act%res_loc_wp => this%res_loc_wp
    this%act%x_loc_wp => this%x_loc_wp
    IF (ltimer) THEN
      CALL timer_stop(this%timer_init)
      IF (.NOT.trans%is_solver_pe) THEN
        CALL timer_start(this%timer)
        CALL timer_start(this%act%lhs%timer)
        CALL timer_stop(this%act%lhs%timer)
        CALL timer_start(trans%timer_glob_sum)
        CALL timer_stop(trans%timer_glob_sum)
        CALL timer_start(trans%timer_sync)
        CALL timer_stop(trans%timer_sync)
        CALL timer_stop(this%timer)
      END IF
    END IF
    this%is_init = .true.
  END SUBROUTINE ocean_solve_construct

! general interface for solve
  SUBROUTINE ocean_solve_solve(this, niter, niter_sp, lacc)
    CLASS(t_ocean_solve), INTENT(INOUT) :: this
    INTEGER, INTENT(OUT) :: niter, niter_sp
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = this_mod_name// &
      & '::t_ocean_solve::ocean_solve'

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (.NOT.ALLOCATED(this%act)) &
      & CALL finish(routine, "solve needs to be initialized")
    IF (ltimer) CALL timer_start(this%timer)
! update rhs-pointer
    this%act%b_loc_wp => this%b_loc_wp
! call backend
    CALL this%act%solve(niter, niter_sp, MERGE(1, 0, this%sol_type .NE. solve_legacy_gmres), lacc=lzacc)
    IF (ltimer) CALL timer_stop(this%timer)

  END SUBROUTINE ocean_solve_solve

END MODULE mo_ocean_solve
