!>
!! @file test_perf_stripes.f90
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#include "fc_feature_defs.inc"
PROGRAM test_perf_stripes
  USE mpi
  USE yaxt, ONLY: xt_idxlist, xt_idxlist_delete, xt_idxvec_new, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, xt_redist, &
       xt_redist_p2p_new, xt_redist_delete, xt_redist_p2p_off_new, &
       xt_redist_p2p_blocks_off_new, xt_redist_s_exchange, &
       xt_stripe, xt_idxstripes_new, xt_int_kind, xt_initialize, xt_finalize

  USE ftest_common, ONLY: finish_mpi, init_mpi, timer, treset, tstart, &
       tstop, treport, id_map, test_abort, factorize, regular_deco, &
       cmp_arrays

  ! PGI compilers up to at least version 15 do not handle generic
  ! interfaces correctly
#if defined __PGI
  USE xt_redist_int_i2, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange
#endif
  IMPLICIT NONE
  ! global extents including halos:

  INTEGER, PARAMETER :: nlev = 20
  INTEGER, PARAMETER :: undef_int = - 1
  INTEGER(xt_int_kind), PARAMETER :: undef_index = -1
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size

  INTEGER, PARAMETER :: grid_kind_test = 1
  INTEGER, PARAMETER :: grid_kind_toy  = 2
  INTEGER, PARAMETER :: grid_kind_tp10 = 3
  INTEGER, PARAMETER :: grid_kind_tp04 = 4
  INTEGER, PARAMETER :: grid_kind_tp6M = 5
  INTEGER :: grid_kind = grid_kind_test

  CHARACTER(len=10) :: grid_label
  INTEGER :: g_ie, g_je ! global domain extents
  INTEGER :: ie, je, ke ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  ! process rank, process coords within (0:, 0:) process space
  INTEGER :: mype, mypx, mypy
  LOGICAL :: lroot ! true only for proc 0

  INTEGER, ALLOCATABLE :: g_id(:,:) ! global id
  ! global "tripolar-like" toy bounds exchange
  INTEGER, ALLOCATABLE :: g_tpex(:, :)
  TYPE(xt_xmap) :: xmap_tpex_2d, xmap_tpex_3d, xmap_tpex_3d_ws
  TYPE(xt_redist) :: redist_tpex_2d, redist_surf_tpex_2d, redist_tpex_3d, &
       redist_tpex_3d_ws, redist_tpex_3d_wb
  TYPE(Xt_idxlist) :: loc_id_3d_ws, loc_tpex_3d_ws

  INTEGER(xt_int_kind), ALLOCATABLE :: loc_id_2d(:,:), loc_tpex_2d(:,:)
  INTEGER(xt_int_kind), ALLOCATABLE :: loc_id_3d(:,:,:), loc_tpex_3d(:,:,:)
  INTEGER, ALLOCATABLE :: fval_2d(:,:), gval_2d(:,:)
  INTEGER, ALLOCATABLE :: fval_3d(:,:,:), gval_3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)
  LOGICAL, PARAMETER :: full_test = .TRUE.
  LOGICAL :: verbose
  CHARACTER(len=*), PARAMETER :: filename = 'test_perf_stripes.f90'

  TYPE(timer) :: t_all, t_surf_redist, t_exch_surf
  TYPE(timer) :: t_xmap_2d, t_redist_2d, t_exch_2d
  TYPE(timer) :: t_xmap_3d, t_redist_3d, t_exch_3d
  TYPE(timer) :: t_xmap_3d_ws, t_redist_3d_ws, t_exch_3d_ws, t_exch_3d_wb
  TYPE(timer) :: t_redist_3d_wb

  !WRITE(0,*) '(debug) test_perf_stripes: verbose=', verbose

  CALL treset(t_all, 'all')
  CALL treset(t_surf_redist, 'surf_redist')
  CALL treset(t_exch_surf, 'exch_surf')
  CALL treset(t_xmap_2d, 'xmap_2d')
  CALL treset(t_redist_2d, 'redist_2d')
  CALL treset(t_exch_2d, 'exch_2d')
  CALL treset(t_xmap_3d, 'xmap_3d')
  CALL treset(t_redist_3d, 'redist_3d')
  CALL treset(t_exch_3d, 'exch_3d')

  CALL treset(t_xmap_3d_ws, 'xmap_3d_ws')
  CALL treset(t_redist_3d_ws, 'redist_3d_ws')
  CALL treset(t_exch_3d_ws, 'exch_3d_ws')

  CALL treset(t_redist_3d_wb, 'redist_3d_wb')
  CALL treset(t_exch_3d_wb, 'exch_3d_wb')

  CALL init_mpi

  CALL tstart(t_all)

  ! mpi & decomposition & allocate mem:
  CALL init_all

  ALLOCATE(fval_3d(nlev,ie,je), gval_3d(nlev,ie,je))

  ! full global index space:
  CALL id_map(g_id)

  ! local window of global index space:
  CALL get_window(g_id, loc_id_2d)

  ! define bounds exchange for full global index space
  CALL def_exchange(g_id, g_tpex)
  !g_tpex = g_id
  DEALLOCATE(g_id)
  ! local window of global bounds exchange:
  CALL get_window(g_tpex, loc_tpex_2d)
  DEALLOCATE(g_tpex)

  IF (full_test) THEN
  ! xmap: loc_id_2d -> loc_tpex_2d
  CALL tstart(t_xmap_2d)
  CALL gen_xmap_2d(loc_id_2d, loc_tpex_2d, xmap_tpex_2d)
  CALL tstop(t_xmap_2d)

  ! transposition: loc_id_2d:data -> loc_tpex_2d:data
  CALL tstart(t_redist_2d)
  CALL gen_redist(xmap_tpex_2d, MPI_INTEGER, MPI_INTEGER, redist_tpex_2d)
  CALL tstop(t_redist_2d)

  ! test 2d-to-2d transposition:
  fval_2d = INT(loc_id_2d)
  CALL tstart(t_exch_2d)
  CALL xt_redist_s_exchange(redist_tpex_2d, fval_2d, gval_2d)
  CALL tstop(t_exch_2d)
  IF (cmp_arrays(gval_2d, loc_tpex_2d)) &
       CALL test_abort('array eqivalence test failed', filename, __LINE__)

  ! define positions of surface elements within (k,i,j) array
  CALL id_map(id_pos)
  CALL id_map(pos3d_surf)
  CALL gen_pos3d_surf(pos3d_surf)

  ! generate surface transposition:
  CALL tstart(t_surf_redist)
  CALL gen_off_redist(xmap_tpex_2d, MPI_INTEGER, id_pos(:,:)-1, &
       MPI_INTEGER, pos3d_surf(:,:)-1, redist_surf_tpex_2d)
  CALL tstop(t_surf_redist)
  DEALLOCATE(id_pos, pos3d_surf)

  ! 2d to surface boundsexchange:
  gval_3d = -1
  CALL tstart(t_exch_surf)
  CALL xt_redist_s_exchange(redist_surf_tpex_2d, &
       RESHAPE(fval_2d, (/ 1, ie, je /)), gval_3d)
  CALL tstop(t_exch_surf)

  ! check surface:
  IF (cmp_arrays(gval_3d(1, :, :), loc_tpex_2d)) &
       CALL test_abort('surface check failed', filename, __LINE__)

  IF (nlev>1) THEN
    ! check sub surface:
    IF (ANY(gval_3d(2, :, :) /= -1)) &
         CALL test_abort('surface check failed', filename, __LINE__)
  ENDIF
  endif

  ! inflate (i,j) -> (k,i,j)
  CALL inflate_idx(1, loc_id_2d, loc_id_3d)
  DEALLOCATE(loc_id_2d)
  CALL inflate_idx(1, loc_tpex_2d, loc_tpex_3d)
  DEALLOCATE(loc_tpex_2d)

  IF (full_test) THEN

  ! xmap: loc_id_3d -> loc_tpex_3d
  CALL tstart(t_xmap_3d)
  CALL gen_xmap_3d(loc_id_3d, loc_tpex_3d, xmap_tpex_3d)
  CALL tstop(t_xmap_3d)

  ! transposition: loc_id_3d:data -> loc_tpex_3d:data
  CALL tstart(t_redist_3d)
  CALL gen_redist(xmap_tpex_3d, MPI_INTEGER, MPI_INTEGER, redist_tpex_3d)
  CALL tstop(t_redist_3d)

  CALL xt_xmap_delete(xmap_tpex_3d)

  ! test 3d-to-3d transposition:
  fval_3d = INT(loc_id_3d)
  gval_3d = -1
  CALL tstart(t_exch_3d)
  CALL xt_redist_s_exchange(redist_tpex_3d, fval_3d, gval_3d)
  CALL tstop(t_exch_3d)

  CALL xt_redist_delete(redist_tpex_3d)

  ! check 3d exchange:
  IF (cmp_arrays(gval_3d, loc_tpex_3d)) &
       CALL test_abort('3D array eqivalence test failed', filename, __LINE__)
  endif


  ! gen stripes, xmap, redist:
  CALL gen_stripes(loc_id_3d, loc_id_3d_ws)
  CALL gen_stripes(loc_tpex_3d, loc_tpex_3d_ws)

  xmap_tpex_3d_ws = xt_xmap_all2all_new(loc_id_3d_ws, loc_tpex_3d_ws, &
       MPI_COMM_WORLD)

  CALL tstart(t_redist_3d_ws)
  redist_tpex_3d_ws  = xt_redist_p2p_new(xmap_tpex_3d_ws, MPI_INTEGER)
  CALL tstop(t_redist_3d_ws)

  ! test redist_tpex_3d_ws:
  fval_3d = INT(loc_id_3d)
  gval_3d = -1

  CALL tstart(t_exch_3d_ws)
  CALL xt_redist_s_exchange(redist_tpex_3d_ws, fval_3d, gval_3d)
  CALL tstop(t_exch_3d_ws)
  if (full_test) then
  ! check 3d exchange:
  IF (cmp_arrays(gval_3d, loc_tpex_3d)) &
       CALL test_abort('3D array eqivalence test (using stripes) failed', &
       filename, __LINE__)
  endif

  CALL tstart(t_redist_3d_wb)
  CALL gen_redist_3d_wb(xmap_tpex_2d, MPI_INTEGER, redist_tpex_3d_wb)
  CALL tstop(t_redist_3d_wb)

  fval_3d = INT(loc_id_3d)
  gval_3d = -1
  CALL tstart(t_exch_3d_wb)
  CALL xt_redist_s_exchange(redist_tpex_3d_wb, fval_3d, gval_3d)
  CALL tstop(t_exch_3d_wb)

  DEALLOCATE(fval_3d)
  CALL xt_redist_delete(redist_tpex_3d_wb)

  IF (cmp_arrays(gval_3d, loc_tpex_3d)) &
       CALL test_abort('3D array eqivalence test after redist failed', &
       filename, __LINE__)

  ! cleanup:
  IF (full_test) THEN
    CALL xt_redist_delete(redist_tpex_3d_ws)
    CALL xt_xmap_delete(xmap_tpex_3d_ws)
    CALL xt_xmap_delete(xmap_tpex_2d)
    CALL xt_redist_delete(redist_tpex_2d)
    CALL xt_redist_delete(redist_surf_tpex_2d)
  ENDIF

  CALL tstop(t_all)

  IF (verbose) WRITE(0,*) 'timer report for nprocs=',nprocs

  CALL treport(t_all, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_surf_redist, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_surf, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_xmap_2d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_redist_2d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_2d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_xmap_3d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_redist_3d, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_3d, TRIM(grid_label), mpi_comm_world)

  CALL treport(t_xmap_3d_ws, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_redist_3d_ws, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_3d_ws, TRIM(grid_label), mpi_comm_world)

  CALL treport(t_redist_3d_wb, TRIM(grid_label), mpi_comm_world)
  CALL treport(t_exch_3d_wb, TRIM(grid_label), mpi_comm_world)

  DEALLOCATE(loc_tpex_3d, loc_id_3d, fval_2d, gval_2d, gval_3d)
  CALL xt_finalize()
  CALL finish_mpi

CONTAINS

  SUBROUTINE gen_redist_3d_wb(xmap_2d, dt, redist_3d)
    TYPE(xt_xmap), INTENT(in) :: xmap_2d
    INTEGER, INTENT(in) :: dt
    TYPE(xt_redist), INTENT(out) :: redist_3d

    INTEGER :: block_disp(ie,je), block_size(ie,je)
    INTEGER :: i, j
    ! data(k,i,j)
    DO j = 1, je
      DO i = 1, ie
        block_disp(i,j) = ( (j-1) * ie + i - 1 ) * nlev
        block_size(i,j) =  nlev
      ENDDO
    ENDDO
    !WRITE(0,*) '(gen_redist_3d_wb) call redist with field sizes =',ie*je
    redist_3d = xt_redist_p2p_blocks_off_new(xmap_2d, block_disp, block_size, &
         SIZE(block_size), block_disp, block_size, SIZE(block_size),dt)

  END SUBROUTINE gen_redist_3d_wb

  SUBROUTINE inflate_idx(inflate_pos, idx_2d, idx_3d)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::inflate_idx: '
    INTEGER, INTENT(in) :: inflate_pos
    INTEGER(xt_int_kind), INTENT(in) :: idx_2d(:,:)
    INTEGER(xt_int_kind), ALLOCATABLE, INTENT(out) :: idx_3d(:,:,:)

    INTEGER :: i, j, k

    SELECT CASE(inflate_pos)
    CASE(1)
      ALLOCATE(idx_3d(ke, ie, je))
      DO j=1,je
        DO i=1,ie
          DO k=1,ke
            idx_3d(k,i,j) = INT(k + (idx_2d(i,j)-1) * ke, xt_int_kind)
          ENDDO
        ENDDO
      ENDDO
    CASE(3)
      ALLOCATE(idx_3d(ie, je, ke))
      DO k=1,ke
        DO j=1,je
          DO i=1,ie
            idx_3d(i,j,k) = INT(idx_2d(i,j) + (k-1) * g_ie * g_je, xt_int_kind)
          ENDDO
        ENDDO
      ENDDO
    CASE DEFAULT
      CALL test_abort(context//' unsupported inflate position', &
           filename, __LINE__)
    END SELECT

  END SUBROUTINE inflate_idx

  SUBROUTINE gen_pos3d_surf(pos)
    INTEGER, INTENT(inout) :: pos(:,:)

    ! positions for zero based arrays ([k,i,j] dim order):
    ! old pos = i + j*ie
    ! new pos = k + (i + j*ie)*nlev

    INTEGER :: ii,jj, i,j,k, p,q

    k = 0 ! surface
    DO jj=1,je
      DO ii=1,ie
        p = pos(ii,jj) - 1 ! shift to 0-based index
        j = p/ie
        i = MOD(p,ie)
        q = k +  (i + j*ie)*nlev
        pos(ii,jj) = q + 1 ! shift to 1-based index
      ENDDO
    ENDDO

  END SUBROUTINE gen_pos3d_surf

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror
    CHARACTER(len=20) :: grid_str

    CALL xt_initialize(MPI_COMM_WORLD)

    CALL get_environment_variable('YAXT_TEST_PERF_GRID', grid_str)

    verbose = .TRUE.

    SELECT CASE (TRIM(ADJUSTL(grid_str)))
    CASE('TOY')
      grid_kind = grid_kind_toy
      grid_label = 'TOY'
      g_ie = 66
      g_je = 36
    CASE('TP10')
      grid_kind = grid_kind_tp10
      grid_label = 'TP10'
      g_ie = 362
      g_je = 192
    CASE('TP04')
      grid_kind = grid_kind_tp04
      grid_label = 'TP04'
      g_ie = 802
      g_je = 404
    CASE('TP6M')
      grid_kind = grid_kind_tp6m
      grid_label = 'TP6M'
      g_ie = 3602
      g_je = 2394
    CASE default
      grid_kind = grid_kind_test
      grid_label = 'TEST'
      g_ie = 32
      g_je = 12
      verbose = .FALSE.
    END SELECT

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'MPI_COMM_SIZE failed', filename, __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) &
         CALL test_abort(context//'MPI_COMM_RANK failed', filename, __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
      verbose = .FALSE.
    ENDIF

    CALL factorize(nprocs, nprocx, nprocy)
    IF (lroot .AND. verbose) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    IF (lroot .AND. verbose) WRITE(0,*) 'g_ie, g_je=',g_ie, g_je
    mypy = mype / nprocx
    mypx = MOD(mype, nprocx)

    CALL deco
    ke = nlev

    ALLOCATE(g_id(g_ie, g_je), g_tpex(g_ie, g_je))

    ALLOCATE(fval_2d(ie,je), gval_2d(ie,je))
    ALLOCATE(loc_id_2d(ie,je), loc_tpex_2d(ie,je))
    ALLOCATE(id_pos(ie,je), pos3d_surf(ie,je))

    fval_2d = undef_int
    gval_2d = undef_int
    loc_id_2d = INT(undef_int, xt_int_kind)
    loc_tpex_2d = INT(undef_int, xt_int_kind)
    id_pos = undef_int
    pos3d_surf = undef_int

  END SUBROUTINE init_all

  SUBROUTINE gen_redist(xmap, send_dt, recv_dt, redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: send_dt, recv_dt
    TYPE(xt_redist),INTENT(out) :: redist

    INTEGER :: dt

    IF (send_dt /= recv_dt) &
         CALL test_abort('gen_redist: (send_dt /= recv_dt) unsupported', &
         filename, __LINE__)
    dt = send_dt
    redist = xt_redist_p2p_new(xmap, dt)

  END SUBROUTINE gen_redist

  SUBROUTINE gen_off_redist(xmap, send_dt, send_off, recv_dt, recv_off, redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER,INTENT(in) :: send_dt, recv_dt
    INTEGER,INTENT(in) :: send_off(:,:), recv_off(:,:)
    TYPE(xt_redist),INTENT(out) :: redist

    INTEGER :: dt

    IF (send_dt /= recv_dt) &
         CALL test_abort('gen_off_redist: (send_dt /= recv_dt) unsupported', &
         filename, __LINE__)
    dt = send_dt

    redist = xt_redist_p2p_off_new(xmap, send_off, recv_off, dt)
  END SUBROUTINE gen_off_redist

  SUBROUTINE get_window(gval, win)
    INTEGER, INTENT(in) :: gval(:,:)
    INTEGER(xt_int_kind), INTENT(out) :: win(:,:)

    INTEGER :: i, j, ig, jg

    DO j = 1, je
      jg = p_joff + j
      DO i = 1, ie
        ig = p_ioff + i
        win(i,j) = INT(gval(ig,jg), xt_int_kind)
      ENDDO
    ENDDO

  END SUBROUTINE get_window

  SUBROUTINE gen_stripes(local_idx, local_stripes)
    CHARACTER(len=*), PARAMETER :: context = 'gen_stripes: '

    INTEGER(xt_int_kind), INTENT(in) :: local_idx(:,:,:)
    TYPE(Xt_idxlist), INTENT(out) :: local_stripes

    TYPE(xt_stripe), ALLOCATABLE :: stripes(:,:)
    INTEGER :: i, j, k, ni, nj, nk

    ! FIXME: assert nk matches xt_int_kind representable values
    nk = SIZE(local_idx,1)
    ni = SIZE(local_idx,2)
    nj = SIZE(local_idx,3)

    ALLOCATE(stripes(ni,nj))

    DO j = 1, nj
      DO i = 1, ni
        ! start, nstrides, stride
        stripes(i,j) = xt_stripe(local_idx(1,i,j), 1_xt_int_kind, nk)
        DO k = 1, nk
          IF (local_idx(1,i,j)-1+k /= local_idx(k,i,j)) &
               CALL test_abort(context//'stripe condition violated', &
               filename, __LINE__)
        ENDDO
      ENDDO
    ENDDO

    local_stripes = xt_idxstripes_new(stripes, SIZE(stripes))

  END SUBROUTINE gen_stripes

  SUBROUTINE gen_xmap_2d(local_src_idx, local_dst_idx, xmap)
    INTEGER(xt_int_kind), INTENT(in) :: local_src_idx(:,:)
    INTEGER(xt_int_kind), INTENT(in) :: local_dst_idx(:,:)
    TYPE(xt_xmap), INTENT(out) :: xmap

    TYPE(Xt_idxlist) :: src_idxlist, dst_idxlist

    src_idxlist = xt_idxvec_new(local_src_idx, SIZE(local_src_idx))
    dst_idxlist = xt_idxvec_new(local_dst_idx, SIZE(local_dst_idx))
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist,  MPI_COMM_WORLD)

    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)

  END SUBROUTINE gen_xmap_2d

  SUBROUTINE gen_xmap_3d(local_src_idx, local_dst_idx, xmap)
    INTEGER(xt_int_kind), INTENT(in) :: local_src_idx(:,:,:)
    INTEGER(xt_int_kind), INTENT(in) :: local_dst_idx(:,:,:)
    TYPE(xt_xmap), INTENT(out) :: xmap

    TYPE(Xt_idxlist) :: src_idxlist, dst_idxlist

    src_idxlist = xt_idxvec_new(local_src_idx, SIZE(local_src_idx))
    dst_idxlist = xt_idxvec_new(local_dst_idx, SIZE(local_dst_idx))
    xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist,  MPI_COMM_WORLD)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(dst_idxlist)


  END SUBROUTINE gen_xmap_3d


  SUBROUTINE def_exchange(g_id, g_tpex)
    INTEGER, INTENT(in) :: g_id(:, :)
    INTEGER, INTENT(out) :: g_tpex(:, :)
    LOGICAL, PARAMETER :: increased_north_halo = .FALSE.
    LOGICAL, PARAMETER :: with_north_halo = .true.
    INTEGER :: i, j
    INTEGER :: g_core_is, g_core_ie, g_core_js, g_core_je
    INTEGER :: north_halo

    ! global core domain:
    g_core_is = nhalo + 1
    g_core_ie = g_ie-nhalo
    g_core_js = nhalo + 1
    g_core_je = g_je-nhalo

    ! global tripolar boundsexchange:
    g_tpex = undef_index
    g_tpex(g_core_is:g_core_ie, g_core_js:g_core_je) &
         = g_id(g_core_is:g_core_ie, g_core_js:g_core_je)

    IF (with_north_halo) THEN

      ! north inversion, (maybe with increased north halo)
      IF (increased_north_halo) THEN
        north_halo = nhalo+1
      ELSE
        north_halo = nhalo
      ENDIF

      IF (2*north_halo > g_core_je) &
           CALL test_abort('def_exchange: grid too small (or halo too large&
           &) for tripolar north exchange', &
           filename, __LINE__)
      DO j = 1, north_halo
        DO i = g_core_is, g_core_ie
          g_tpex(i,j) = g_tpex(g_core_ie + (g_core_is-i), 2*north_halo + (1-j))
        ENDDO
      ENDDO

    ELSE

      DO j = 1, nhalo
        DO i = nhalo+1, g_ie-nhalo
          g_tpex(i,j) = g_id(i,j)
        ENDDO
      ENDDO

    ENDIF

    ! south: no change
    DO j = g_core_je+1, g_je
      DO i = nhalo+1, g_ie-nhalo
        g_tpex(i,j) = g_id(i,j)
      ENDDO
    ENDDO

    ! PBC
    DO j = 1, g_je
      DO i = 1, nhalo
        g_tpex(g_core_is-i,j) = g_tpex(g_core_ie+(1-i),j)
      ENDDO
      DO i = 1, nhalo
        g_tpex(g_core_ie+i,j) = g_tpex(nhalo+i,j)
      ENDDO
    ENDDO

    CALL check_g_idx (g_tpex)

  END SUBROUTINE def_exchange

  SUBROUTINE check_g_idx(gidx)
    INTEGER,INTENT(in) :: gidx(:,:)

    IF (ANY(gidx == undef_index)) THEN
      CALL test_abort('check_g_idx: check failed', filename, __LINE__)
    ENDIF
  END SUBROUTINE check_g_idx

  SUBROUTINE deco
    INTEGER :: cx0(0:nprocx-1), cxn(0:nprocx-1)
    INTEGER :: cy0(0:nprocy-1), cyn(0:nprocy-1)

    cx0 = 0
    cxn = 0
    CALL regular_deco(g_ie-2*nhalo, cx0, cxn)

    cy0 = 0
    cyn = 0
    CALL regular_deco(g_je-2*nhalo, cy0, cyn)

    ! process local deco variables:
    ie = cxn(mypx) + 2*nhalo
    je = cyn(mypy) + 2*nhalo
    p_ioff = cx0(mypx)
    p_joff = cy0(mypy)

  END SUBROUTINE deco

END PROGRAM test_perf_stripes
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
