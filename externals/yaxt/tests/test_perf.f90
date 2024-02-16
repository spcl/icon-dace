!>
!! @file test_perf.f90
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
PROGRAM test_perf
  USE mpi
  USE yaxt, ONLY: xt_finalize, xt_idxlist, xt_idxvec_new, xt_xmap, &
       xt_xmap_all2all_new, xt_redist, xt_redist_p2p_new, &
       xt_redist_p2p_off_new, xt_redist_s_exchange, xt_idxlist_delete, &
       xt_xmap_delete, xt_redist_delete, xt_int_kind, xt_initialize, &
       xi => xt_int_kind
  USE ftest_common, ONLY: finish_mpi, init_mpi, treset, tstart, tstop, &
       treport, timer, id_map, factorize, regular_deco, set_verbose, icmp, &
       test_abort
  USE test_redist_common, only: xt_redist_s_exchange
  USE yaxt, ONLY: xt_finalize
#ifdef __PGI
  ! PGI up to at least 15.4 has a bug that prevents proper import of
  ! multiply extended generics. This is a separate bug from the one exhibited
  ! in 12.7 and older (see test_xmap_intersection_parallel_f.f90 for that)
  USE xt_redist_int_i2, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange
#endif
  IMPLICIT NONE
  ! global extents including halos:

  INTEGER, PARAMETER :: nlev = 30
  INTEGER, PARAMETER :: undef_int = (HUGE(undef_int)-1)/2 - 1
  INTEGER(xt_int_kind), PARAMETER :: undef_index = -1_xi
  INTEGER, PARAMETER :: nhalo = 1 ! 1dim. halo border size

  CHARACTER(len=8) :: grid_label
  INTEGER :: g_ie, g_je ! global domain extents
  INTEGER :: ie, je, ke ! local extents, including halos
  INTEGER :: p_ioff, p_joff ! offsets within global domain
  INTEGER :: nprocx, nprocy ! process space extents
  INTEGER :: nprocs ! == nprocx*nprocy
  ! process rank, process coords within (0:, 0:) process space
  INTEGER :: mype, mypx, mypy
  LOGICAL :: lroot ! true only for proc 0
  CHARACTER(len=*), PARAMETER :: filename = 'test_perf.f90'

  INTEGER(xt_int_kind), ALLOCATABLE :: g_id(:,:) ! global id
  ! global "tripolar-like" toy bounds exchange
  INTEGER(xt_int_kind), ALLOCATABLE :: g_tpex(:, :)

  INTEGER(xt_int_kind), ALLOCATABLE :: loc_id_2d(:,:), loc_tpex_2d(:,:)
  INTEGER(xt_int_kind), ALLOCATABLE :: loc_id_3d(:,:,:), loc_tpex_3d(:,:,:)
  INTEGER, ALLOCATABLE :: fval_2d(:,:), gval_2d(:,:)
  INTEGER, ALLOCATABLE :: fval_3d(:,:,:), gval_3d(:,:,:)
  INTEGER, ALLOCATABLE :: id_pos(:,:), pos3d_surf(:,:)
  INTEGER, ALLOCATABLE :: ref_m1(:,:)
  LOGICAL :: verbose

  TYPE(timer) :: t_all, t_surf_redist, t_exch_surf
  TYPE(timer) :: t_idxlist_2d, t_xmap_2d, t_redist_2d, t_exch_2d
  TYPE(timer) :: t_idxlist_3d, t_xmap_3d, t_redist_3d, t_exch_3d

  ! yaxt data:
  TYPE(xt_idxlist) :: id_2d_idxlist, tpex_2d_idxlist
  TYPE(xt_idxlist) :: id_3d_idxlist, tpex_3d_idxlist
  TYPE(xt_xmap) :: tpex_2d_xmap, tpex_3d_xmap
  TYPE(xt_redist) :: tpex_2d_redist, surf_tpex_redist
  TYPE(xt_redist) :: tpex_3d_redist

  CALL treset(t_all, 'all')
  CALL treset(t_surf_redist, 'surf_redist')
  CALL treset(t_exch_surf, 'exch_surf')
  CALL treset(t_idxlist_2d, 'idxlist_2d')
  CALL treset(t_xmap_2d, 'xmap_2d')
  CALL treset(t_redist_2d, 'redist_2d')
  CALL treset(t_exch_2d, 'exch_2d')
  CALL treset(t_idxlist_3d, 'idxlist_3d')
  CALL treset(t_xmap_3d, 'xmap_3d')
  CALL treset(t_redist_3d, 'redist_3d')
  CALL treset(t_exch_3d, 'exch_3d')

  CALL init_mpi

  CALL tstart(t_all)

  ! mpi & decomposition & allocate mem:
  CALL init_all

  ! full global index space:
  CALL id_map(g_id)

  ! local window of global index space:
  CALL get_window(g_id, loc_id_2d)


  ! define bounds exchange for full global index space
  CALL def_exchange(g_id, g_tpex)
  DEALLOCATE(g_id)
  ! local window of global bounds exchange:
  CALL get_window(g_tpex, loc_tpex_2d)
  DEALLOCATE(g_tpex)

  ! gen xmap: loc_id_2d -> loc_tpex_2d
  CALL tstart(t_idxlist_2d)
  id_2d_idxlist = xt_idxvec_new(loc_id_2d)
  tpex_2d_idxlist = xt_idxvec_new(loc_tpex_2d)
  CALL tstop(t_idxlist_2d)
  CALL tstart(t_xmap_2d)
  tpex_2d_xmap = xt_xmap_all2all_new(id_2d_idxlist, tpex_2d_idxlist, MPI_COMM_WORLD)
  CALL tstop(t_xmap_2d)

  ! gen redist: loc_id_2d:data -> loc_tpex_2d:data
  CALL tstart(t_redist_2d)
  tpex_2d_redist = xt_redist_p2p_new(tpex_2d_xmap,  MPI_INTEGER)
  CALL tstop(t_redist_2d)

  ! test 2d-to-2d redist
  CALL copy_xi(loc_id_2d, fval_2d)
  CALL tstart(t_exch_2d)
  CALL xt_redist_s_exchange(tpex_2d_redist, fval_2d, gval_2d)
  CALL tstop(t_exch_2d)
  CALL icmp('2d to 2d check', gval_2d, loc_tpex_2d, mype)

  ! define positions of surface elements within (i,k,j) array
  CALL id_map(id_pos, 0)
  CALL id_map(pos3d_surf, 0)
  CALL gen_pos3d_surf(pos3d_surf)

  ! gen surface redist:
  CALL tstart(t_surf_redist)
  surf_tpex_redist = xt_redist_p2p_off_new(tpex_2d_xmap, id_pos, pos3d_surf, MPI_INTEGER)
  CALL tstop(t_surf_redist)
  DEALLOCATE(id_pos, pos3d_surf)

  ! 2d to surface boundsexchange:
  ALLOCATE(gval_3d(ie,nlev,je))
  gval_3d = -1
  CALL tstart(t_exch_surf)
  CALL xt_redist_s_exchange(surf_tpex_redist, fval_2d, gval_3d)
  CALL tstop(t_exch_surf)

  CALL icmp('surface check', gval_3d(:,1,:), loc_tpex_2d, mype)

  ! check sub surface:
  ALLOCATE(ref_m1(SIZE(loc_tpex_2d, 1), SIZE(loc_tpex_2d, 2)))
  ref_m1 = -1
  CALL icmp('sub surface check', gval_3d(:,2,:), ref_m1, mype)

  DEALLOCATE(ref_m1)
  ! cleanup 2d:
  CALL xt_idxlist_delete(id_2d_idxlist)
  CALL xt_idxlist_delete(tpex_2d_idxlist)
  CALL xt_xmap_delete(tpex_2d_xmap)
  CALL xt_redist_delete(surf_tpex_redist)
  CALL xt_redist_delete(tpex_2d_redist)

  ! inflate 2d -> 3d
  CALL inflate_idx(3, loc_id_2d, loc_id_3d)
  DEALLOCATE(loc_id_2d)
  CALL inflate_idx(3, loc_tpex_2d, loc_tpex_3d)
  DEALLOCATE(loc_tpex_2d)

  ! gen xmap: loc_id_3d -> loc_tpex_3d
  CALL tstart(t_idxlist_3d)
  id_3d_idxlist = xt_idxvec_new(loc_id_3d)
  tpex_3d_idxlist = xt_idxvec_new(loc_tpex_3d)
  CALL tstop(t_idxlist_3d)
  CALL tstart(t_xmap_3d)
  tpex_3d_xmap = xt_xmap_all2all_new(id_3d_idxlist, tpex_3d_idxlist, MPI_COMM_WORLD)
  CALL tstop(t_xmap_3d)

  ! gen redist: loc_id_3d:data -> loc_tpex_3d:data
  CALL tstart(t_redist_3d)
  tpex_3d_redist = xt_redist_p2p_new(tpex_3d_xmap,  MPI_INTEGER)
  CALL tstop(t_redist_3d)

  ! test 3d-to-3d redist:
  DEALLOCATE(gval_3d)
  ALLOCATE(fval_3d(ie,je,nlev), gval_3d(ie,je,nlev))
  CALL copy_xi_3d(loc_id_3d, fval_3d)
  CALL tstart(t_exch_3d)
  CALL xt_redist_s_exchange(tpex_3d_redist, fval_3d, gval_3d)
  CALL tstop(t_exch_3d)
  DEALLOCATE(fval_3d)
  CALL icmp('3d to 3d check', gval_3d, loc_tpex_3d, mype)
  DEALLOCATE(gval_3d)

  ! cleanup 3d:
  CALL xt_idxlist_delete(id_3d_idxlist)
  CALL xt_idxlist_delete(tpex_3d_idxlist)
  CALL xt_xmap_delete(tpex_3d_xmap)
  CALL xt_redist_delete(tpex_3d_redist)

  CALL tstop(t_all)

  IF (verbose) WRITE(0,*) 'timer report for nprocs=',nprocs

  CALL tvec_report((/t_all, t_surf_redist, t_exch_surf, t_idxlist_2d, t_xmap_2d, &
       t_redist_2d, t_exch_2d, t_idxlist_3d, t_xmap_3d,t_redist_3d, t_exch_3d/), &
       TRIM(grid_label), mpi_comm_world)

  DEALLOCATE(loc_tpex_3d, loc_id_3d, fval_2d, gval_2d)

  ! finalize:
  CALL xt_finalize
  CALL finish_mpi

CONTAINS

  SUBROUTINE copy_xi(s, d)
    INTEGER(xt_int_kind), INTENT(in) :: s(:,:)
    INTEGER, INTENT(out) :: d(:,:)
    INTEGER :: i, j, m, n
    m = SIZE(s, 1)
    n = SIZE(s, 2)
    DO j = 1, n
      DO i = 1, m
        d(i, j) = INT(s(i, j))
      END DO
    END DO
  END SUBROUTINE copy_xi

  SUBROUTINE copy_xi_3d(s, d)
    INTEGER(xt_int_kind), INTENT(in) :: s(:,:,:)
    INTEGER, INTENT(out) :: d(:,:,:)
    INTEGER :: i, j, k, m, n, o
    m = SIZE(s, 1)
    n = SIZE(s, 2)
    o = SIZE(s, 3)
    DO k = 1, o
      DO j = 1, n
        DO i = 1, m
          d(i, j, k) = INT(s(i, j, k))
        END DO
      END DO
    END DO
  END SUBROUTINE copy_xi_3d

  SUBROUTINE tvec_report(tvec, extra_label, comm)
    TYPE(timer), INTENT(in) :: tvec(:)
    CHARACTER(len=*), INTENT(in) :: extra_label
    INTEGER, INTENT(in) :: comm
    INTEGER :: it
    DO it = 1, SIZE(tvec)
      CALL treport(tvec(it), extra_label, comm)
    ENDDO
  END SUBROUTINE tvec_report

  SUBROUTINE inflate_idx(inflate_pos, idx_2d, idx_3d)
    CHARACTER(len=*), PARAMETER :: context = 'test_perf::inflate_idx: '
    INTEGER, INTENT(in) :: inflate_pos
    INTEGER(xt_int_kind), INTENT(in) :: idx_2d(:,:)
    INTEGER(xt_int_kind), ALLOCATABLE, INTENT(out) :: idx_3d(:,:,:)

    INTEGER :: i, j, k

    IF (inflate_pos == 3) THEN
      ALLOCATE(idx_3d(ie, je, ke))
      DO k=1,ke
        DO j=1,je
          DO i=1,ie
            idx_3d(i,j,k) = idx_2d(i,j) + INT(k-1,xi) * INT(g_ie, xi)&
                 & * INT(g_je, xi)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      CALL test_abort(context//' unsupported inflate position', &
           filename, __LINE__)
    ENDIF

  END SUBROUTINE inflate_idx

  SUBROUTINE gen_pos3d_surf(pos)
    INTEGER, INTENT(inout) :: pos(:,:)
    ! positions for zero based arrays (ECHAM grid point dim order)
    ! old pos = i + j*ie
    ! new pos = i + k*ie + j*ie*nlev
    INTEGER :: ii,jj, i,j,k, p

    k = 0 ! surface
    DO jj=1,je
      DO ii=1,ie
        p = pos(ii,jj)
        j = p/ie
        i = MOD(p,ie)
        pos(ii,jj) = i + k*ie + j*ie*nlev
      ENDDO
    ENDDO

  END SUBROUTINE gen_pos3d_surf

  SUBROUTINE init_all
    CHARACTER(len=*), PARAMETER :: context = 'init_all: '
    INTEGER :: ierror
    CHARACTER(len=20) :: grid_str

    CALL get_environment_variable('YAXT_TEST_PERF_GRID', grid_str)

    verbose = .TRUE.

    SELECT CASE (ADJUSTL(grid_str))
    CASE('TOY')
      grid_label = 'TOY'
      g_ie = 66
      g_je = 36
    CASE('TP10')
      grid_label = 'TP10'
      g_ie = 362
      g_je = 192
    CASE('TP04')
      grid_label = 'TP04'
      g_ie = 802
      g_je = 404
    CASE('TP6M')
      grid_label = 'TP6M'
      g_ie = 3602
      g_je = 2394
    CASE default
      grid_label = 'TEST'
      g_ie = 32
      g_je = 12
      verbose = .FALSE.
    END SELECT

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_COMM_SIZE failed', &
         filename, __LINE__)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
    IF (ierror /= MPI_SUCCESS) CALL test_abort(context//'MPI_COMM_RANK failed', &
         filename, __LINE__)
    IF (mype==0) THEN
      lroot = .true.
    ELSE
      lroot = .FALSE.
      verbose = .FALSE.
    ENDIF
    CALL set_verbose(verbose)
    CALL factorize(nprocs, nprocx, nprocy)
    IF (verbose) WRITE(0,*) 'nprocx, nprocy=',nprocx, nprocy
    IF (verbose) WRITE(0,*) 'g_ie, g_je=',g_ie, g_je
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
    loc_id_2d = undef_index
    loc_tpex_2d = undef_index
    id_pos = undef_int
    pos3d_surf = undef_int

    CALL xt_initialize(MPI_COMM_WORLD)

  END SUBROUTINE init_all

  SUBROUTINE get_window(gval, win)
    INTEGER(xt_int_kind), INTENT(in) :: gval(:,:)
    INTEGER(xt_int_kind), INTENT(out) :: win(:,:)
    INTEGER :: i, j, ig, jg
    DO j = 1, je
      jg = p_joff + j
      DO i = 1, ie
        ig = p_ioff + i
        win(i,j) =  gval(ig,jg)
      ENDDO
    ENDDO
  END SUBROUTINE get_window

  SUBROUTINE def_exchange(g_id, g_tpex)
    INTEGER(xt_int_kind), INTENT(in) :: g_id(:, :)
    INTEGER(xt_int_kind), INTENT(out) :: g_tpex(:, :)
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
           CALL test_abort('def_exchange: grid too small (or halo too large) &
           &for tripolar north exchange', &
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
    INTEGER(xt_int_kind),INTENT(in) :: gidx(:,:)

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

END PROGRAM test_perf
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
