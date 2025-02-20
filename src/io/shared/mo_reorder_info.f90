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

! the Intel compiler has an interesting behaviour: if optimization is off,
! actual arguments to dummy arguments with attribute CONTIGUOUS will always
! be copied
#if ((! defined __INTEL_COMPILER) || defined __OPTIMIZE__)
#define USE_CONTIGUOUS 1
#else
#undef USE_CONTIGUOUS
#endif
#undef USE_CONTIGUOUS
MODULE mo_reorder_info
  USE mo_kind, ONLY: i4, i8, dp, sp
  USE mo_mpi, ONLY: p_bcast, p_comm_remote_size, p_allgather, p_allgatherv, &
       p_comm_size, p_int_i8
  USE mo_exception, ONLY: finish, message_text
  USE mo_communication,             ONLY: idx_no, blk_no
#ifndef NOMPI
  USE mpi
#endif
  USE mo_util_sort, ONLY: quicksort
#ifdef HAVE_CDI_PIO
  USE yaxt, ONLY: xt_idxlist, xt_is_null, xt_idxvec_new, &
       xt_idxlist_delete, xt_idxstripes_from_idxlist_new, &
       xt_int_kind
#endif
  IMPLICIT NONE
  PRIVATE
  !------------------------------------------------------------------------------------------------
  ! TYPE t_reorder_info describes how local cells/edges/verts
  ! have to be reordered to get the global array.
  ! Below, "points" refers to either cells, edges or verts.
  !
  ! TODO[FP] Note that the "reorder_info" contains fields of *global*
  !          size (reorder_index). On the compute PEs these fields
  !          could be deallocated after the call to
  !          "transfer_reorder_info" in the setup phase!

  TYPE t_reorder_info
    INTEGER                    :: n_glb  !< Global number of points per physical patch, per layer
    INTEGER                    :: n_own  ! Number of own points (without halo, only belonging to phyiscal patch)
    !> Only set on compute PEs, set to 0 on IO PEs
    INTEGER, ALLOCATABLE       :: own_idx(:), own_blk(:)
    !> dest idx and blk for own points, only set on sequential/test PEs
    INTEGER, ALLOCATABLE       :: pe_own(:)
    !> n_own, gathered for all compute PEs (set on all PEs)
    INTEGER, ALLOCATABLE       :: pe_off(:)
    !> this ranks reordering indices, i.e. place each packed element goes into
    !! in final array
    INTEGER, ALLOCATABLE       :: reorder_index_own(:)
    !> offset of contributions of all ranks concatenated
    !! (should ONLY be set on I/O servers)
    INTEGER, ALLOCATABLE       :: reorder_index(:)

#ifdef HAVE_CDI_PIO
    ! describe which parts of a global array this MPI rank provides for an
    ! array of nlev in reorder_idxlst_xt(nlev)
    TYPE(xt_idxlist), ALLOCATABLE :: reorder_idxlst_xt(:)
#endif

    ! Index how to reorder the contributions of all compute PEs
    ! into the global array (set on all PEs)
  END TYPE t_reorder_info
  INTERFACE ri_cpy_part2whole
    MODULE PROCEDURE ri_part2whole_1d_dp_dp, ri_part2whole_1d_sp_sp, &
         ri_part2whole_1d_sp_dp, ri_part2whole_1d_dp_sp, ri_part2whole_1d_i4_i4
    MODULE PROCEDURE ri_part2whole_2d_dp_dp, ri_part2whole_2d_sp_sp, &
         ri_part2whole_2d_sp_dp, ri_part2whole_2d_dp_sp, ri_part2whole_2d_i4_i4
  END INTERFACE ri_cpy_part2whole
  INTERFACE ri_cpy_blk2part
    MODULE PROCEDURE ri_blk2part_2d_dp_dp, ri_blk2part_2d_sp_sp, &
         ri_blk2part_2d_sp_dp, ri_blk2part_2d_i4_i4, ri_blk2part_2d_i4_dp
    MODULE PROCEDURE ri_blk2part_3d_dp_dp, ri_blk2part_3d_sp_sp, &
         ri_blk2part_3d_sp_dp, ri_blk2part_3d_i4_i4, ri_blk2part_3d_i4_dp
  END INTERFACE ri_cpy_blk2part
  PUBLIC :: t_reorder_info
  PUBLIC :: transfer_reorder_info
  PUBLIC :: mask2reorder_info
  PUBLIC :: ri_cpy_part2whole
  PUBLIC :: ri_cpy_blk2part
  PUBLIC :: release_reorder_info
  CHARACTER(len=*), PARAMETER :: modname = 'mo_reorder_info'
CONTAINS

  SUBROUTINE release_reorder_info(ri)
    TYPE(t_reorder_info), INTENT(INOUT) :: ri

#ifdef HAVE_CDI_PIO
    INTEGER :: i, n
#endif
    IF (ALLOCATED(ri%reorder_index_own)) DEALLOCATE(ri%reorder_index_own)
    IF (ALLOCATED(ri%reorder_index)) DEALLOCATE(ri%reorder_index)
    IF (ALLOCATED(ri%own_idx)) DEALLOCATE(ri%own_idx)
    IF (ALLOCATED(ri%own_blk)) DEALLOCATE(ri%own_blk)
    IF (ALLOCATED(ri%pe_own)) DEALLOCATE(ri%pe_own)
    IF (ALLOCATED(ri%pe_off)) DEALLOCATE(ri%pe_off)
#ifdef HAVE_CDI_PIO
    IF (ALLOCATED(ri%reorder_idxlst_xt)) THEN
      n = SIZE(ri%reorder_idxlst_xt)
      DO i = 1, n
        IF (.NOT. xt_is_null(ri%reorder_idxlst_xt(i))) &
             CALL xt_idxlist_delete(ri%reorder_idxlst_xt(i))
      END DO
    END IF
#endif

  END SUBROUTINE release_reorder_info

  SUBROUTINE mask2reorder_info(ri, mask, n_points_g, glb_index, group_comm, &
       retained_occupation_mask, create_idxlist)
    TYPE(t_reorder_info), INTENT(inout) :: ri
    LOGICAL, INTENT(in) :: mask(:)
    INTEGER, INTENT(in) :: n_points_g, glb_index(:), group_comm
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: glb_index, mask
#endif
    INTEGER(i8), ALLOCATABLE, OPTIONAL, INTENT(out) :: &
         retained_occupation_mask(:)
    LOGICAL, OPTIONAL, INTENT(in) :: create_idxlist

    INTEGER :: n_points, i, il, n, group_comm_size
    INTEGER :: ierror
    INTEGER, ALLOCATABLE :: glbidx_own(:), permutation(:), occ_pfxsum(:), &
         buf(:)
    INTEGER(i8), ALLOCATABLE :: occupation_mask(:)
    INTEGER(i8), PARAMETER :: nbits_i8 = BIT_SIZE(occupation_mask)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_reorder_data"
#ifdef HAVE_CDI_PIO
    TYPE(xt_idxlist) :: idxvec
#endif
    LOGICAL :: create_idxlist_

    IF (PRESENT(create_idxlist)) THEN
      create_idxlist_ = create_idxlist
    ELSE
      create_idxlist_ = .FALSE.
    END IF
    n_points = SIZE(mask)
    n = COUNT(mask)
    ! Get number of owned cells/edges/verts (without halos, physical patch only)
    ri%n_own = n

    ! Set index arrays to own cells/edges/verts
    ! Global index of my own points
    ALLOCATE(ri%own_idx(n), ri%own_blk(n), ri%reorder_index_own(n))
    group_comm_size = p_comm_size(group_comm)
    ALLOCATE(ri%pe_own(0:group_comm_size-1), ri%pe_off(0:group_comm_size-1))
    ALLOCATE(glbidx_own(n), permutation(n), buf(n))
    n = 0
    DO i = 1, n_points
      IF (mask(i)) THEN
        n = n+1
        ri%own_idx(n) = idx_no(i)
        ri%own_blk(n) = blk_no(i)
        glbidx_own(n) = glb_index(i)
      ENDIF
    ENDDO

    ! sort used global indices so that target-side re-ordering has lower overhead
    DO i = 1, n
      permutation(i) = i
    END DO
    CALL quicksort(glbidx_own, permutation)
    CALL permute(ri%own_idx, permutation, buf)
    CALL permute(ri%own_blk, permutation, buf)
    DEALLOCATE(buf, permutation)
    ! Gather the number of own points for every PE into ri%pe_own
    CALL p_allgather(ri%n_own, ri%pe_own, comm=group_comm)

    ! Get offset within result array
    il = 0
    DO i = 0, group_comm_size-1
      ri%pe_off(i) = il
      il = il + ri%pe_own(i)
    ENDDO

    ! Get global number of points for current (physical!) patch
    ri%n_glb = il

    ! Compute the global ordering of the local data when it is gathered
    ! on PE 0 or I/O servers exactly as it is retrieved later during I/O
    ! 1. set bit for every global index used
    ! todo: occupation_mask is still a global size array but only uses one bit
    ! per grid cell, i.e. decomposition makes sense but takes additional coding
    ! effort
    CALL indices2occmask(occupation_mask, n_points_g, ri%n_own, glbidx_own)
#ifndef NOMPI
    ! 2. reduce occupation over all ranks in group_comm
    CALL mpi_allreduce(mpi_in_place, occupation_mask, SIZE(occupation_mask), &
      &                p_int_i8, mpi_bor, group_comm, ierror)
    IF (ierror /= MPI_SUCCESS) CALL finish(routine, 'mpi_allreduce failed')
#endif

    ! 3. now compute number of bits set in preceding entries of occupation_mask.
    ! After the following loop, occ_pfxsum(i) equals the number of set bits
    ! (i.e. number of global indices used in output) contained in
    ! occupation_mask(0:i-1)
    CALL occmask_pfxsum(occ_pfxsum, occupation_mask, ri%n_glb)
    ! 4. given the above two arrays, one can now compute for each
    ! global index its position in the output array in O(1) note: this
    ! routine produces 0-based indices because that's what's needed
    ! for CDI-PIO
    CALL glb_idx2reorder_idx(ri%reorder_index_own, glbidx_own(1:ri%n_own), &
         occ_pfxsum, occupation_mask)
#ifdef HAVE_CDI_PIO
    IF (create_idxlist_) THEN
      ALLOCATE(ri%reorder_idxlst_xt(1))
      idxvec = xt_idxvec_new(INT(ri%reorder_index_own, xt_int_kind))
      ri%reorder_idxlst_xt(1) = xt_idxstripes_from_idxlist_new(idxvec)
      CALL xt_idxlist_delete(idxvec)
    END IF
#endif
    DO i = 1, ri%n_own
      ri%reorder_index_own(i) = ri%reorder_index_own(i) + 1
    END DO
    DEALLOCATE(glbidx_own)
    IF (PRESENT(retained_occupation_mask)) THEN
      CALL MOVE_ALLOC(occupation_mask, retained_occupation_mask)
    ELSE
      DEALLOCATE(occupation_mask)
    END IF
  END SUBROUTINE mask2reorder_info

  !> create bmask(0:nocc-1) such that bits 0..n_points_g-1 can be
  !! represented and given an array of indices, for each pos in indices,
  !! set the corresponding bit pos-1 in bmask while all other bits are unset
  SUBROUTINE indices2occmask(bmask, n_points_g, nidx, indices)
    INTEGER(i8), ALLOCATABLE, INTENT(out) :: bmask(:)
    INTEGER, INTENT(in) :: n_points_g, nidx, indices(nidx)
    INTEGER(i8) :: pos, apos, bpos
    INTEGER(i8), PARAMETER :: nbits_i8 = BIT_SIZE(pos)
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::indices2occmask"
    INTEGER :: i, n, nocc
    n = nidx
    nocc = INT((INT(n_points_g, i8) + nbits_i8 - 1_i8) / nbits_i8)
    ALLOCATE(bmask(0:nocc-1))
    bmask = 0_i8
    DO i = 1, n
      pos = INT(indices(i), i8)
      apos = (pos - 1_i8)/nbits_i8
      bpos = MOD(pos - 1_i8, nbits_i8)
      bmask(apos) = IBSET(bmask(apos), bpos)
    END DO
  END SUBROUTINE indices2occmask

  !> given a bit-mask bmask(0:nocc-1), compute for each entry i the sum over the
  !! number of bits set in bmask(0:i-1) and store in pfxsum(i)
  SUBROUTINE occmask_pfxsum(pfxsum, bmask, n_glb)
    INTEGER, ALLOCATABLE, INTENT(out) :: pfxsum(:)
    INTEGER(i8), INTENT(in) :: bmask(0:)
    INTEGER, INTENT(in) :: n_glb ! number of bits set in total
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: bmask
#endif
    INTEGER(i8) :: occ_temp, occ_accum
    INTEGER :: i, nocc, ub
    CHARACTER(len=*), PARAMETER :: routine = modname//'::occmask_pfxsum'

    nocc = SIZE(bmask)
    ub = nocc-1
    ALLOCATE(pfxsum(0:ub))
    occ_accum = 0_i8
    DO i = 0, ub
      pfxsum(i) = INT(occ_accum)
      occ_temp = POPCNT(bmask(i))
      occ_accum = occ_accum + occ_temp
    END DO
    IF (occ_accum /= INT(n_glb, i8)) THEN
      WRITE (message_text, '(2(a,i0))') 'Bit-counting failed: n=', occ_accum, &
           & ' /= n_glb=', n_glb
      CALL finish(routine,TRIM(message_text))
    ENDIF
  END SUBROUTINE occmask_pfxsum

  SUBROUTINE glb_idx2reorder_idx(reorder_index, glbidx, pfxsum, bmask)
    INTEGER, INTENT(out) :: reorder_index(:)
    INTEGER, INTENT(in) :: glbidx(:), pfxsum(0:)
    INTEGER(i8), INTENT(in) :: bmask(0:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: reorder_index, glbidx, pfxsum, bmask
#endif
    INTEGER(i8) :: pos, apos, bpos, occ_accum
    INTEGER(i8), PARAMETER :: nbits_i8 = BIT_SIZE(pos)
    INTEGER :: i, n
    n = SIZE(glbidx)
    DO i = 1, n
      pos = INT(glbidx(i), i8)
      apos = (pos - 1_i8)/nbits_i8
      bpos = MOD(pos - 1_i8, nbits_i8)
      occ_accum = IAND(ISHFT(1_i8, bpos) - 1_i8, bmask(apos))
      reorder_index(i) = pfxsum(apos) + POPCNT(occ_accum)
    END DO
  END SUBROUTINE glb_idx2reorder_idx

  SUBROUTINE permute(a, permutation, a_perm)
    INTEGER, ALLOCATABLE, INTENT(inout) :: a(:)
    INTEGER, ALLOCATABLE, INTENT(inout) :: a_perm(:)
    INTEGER, INTENT(in) :: permutation(:)
    CHARACTER(len=*), PARAMETER :: routine=modname//"::permute"
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: permutation
#endif
    INTEGER, ALLOCATABLE :: temp(:)
    INTEGER :: i, n
    n = SIZE(a)
    IF (SIZE(permutation) /= n) THEN
      CALL finish(routine, "invalid permutation")
    END IF
    IF (ALLOCATED(a_perm)) THEN
      IF (SIZE(a_perm) /= n) THEN
        DEALLOCATE(a_perm)
        ALLOCATE(a_perm(n))
      END IF
    ELSE
      ALLOCATE(a_perm(n))
    END IF
    DO i = 1, n
      a_perm(i) = a(permutation(i))
    END DO
    CALL MOVE_ALLOC(a, temp)
    CALL MOVE_ALLOC(a_perm, a)
    CALL MOVE_ALLOC(temp, a_perm)
  END SUBROUTINE permute
  !------------------------------------------------------------------------------------------------
  !> Transfers reorder_info from clients to servers for IO/async latbc
  !
  SUBROUTINE transfer_reorder_info(ri, is_server, bcast_root, &
       client_server_intercomm)
    TYPE(t_reorder_info), INTENT(INOUT) :: ri ! Result: reorder info
    LOGICAL, INTENT(in) :: is_server
    INTEGER, INTENT(in) :: bcast_root, client_server_intercomm

    INTEGER :: other_group_size, dummy(1)
    INTEGER, ALLOCATABLE :: zcounts(:)
    ! Transfer the global number of points, this is not yet known on IO PEs
    CALL p_bcast(ri%n_glb,  bcast_root, client_server_intercomm)

    other_group_size = p_comm_remote_size(client_server_intercomm)
    IF (is_server) THEN

      ! On IO PEs: n_own = 0, own_idx and own_blk are not allocated
      ri%n_own = 0
      ! pe_own/pe_off must be allocated for client_group_size, not for p_n_work
      ALLOCATE(ri%pe_own(0:other_group_size-1))
      ALLOCATE(ri%pe_off(0:other_group_size-1))

      ALLOCATE(ri%reorder_index(ri%n_glb))
    ENDIF

    CALL p_bcast(ri%pe_own, bcast_root, client_server_intercomm)
    CALL p_bcast(ri%pe_off, bcast_root, client_server_intercomm)

    IF (is_server) THEN
      CALL p_allgatherv(dummy(1:0), ri%reorder_index, &
        &               ri%pe_own, ri%pe_off, client_server_intercomm)
    ELSE
      ALLOCATE(zcounts(0:other_group_size-1))
      zcounts = 0
      CALL p_allgatherv(ri%reorder_index_own, dummy(1:0), &
        &               zcounts, zcounts, client_server_intercomm)
    END IF

  END SUBROUTINE transfer_reorder_info

  SUBROUTINE ri_part2whole_1d_dp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:)
    REAL(dp), INTENT(inout) :: whole_data(:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_dp_dp

  SUBROUTINE ri_part2whole_1d_sp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:)
    REAL(sp), INTENT(inout) :: whole_data(:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_sp_sp

  SUBROUTINE ri_part2whole_1d_sp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:)
    REAL(dp), INTENT(inout) :: whole_data(:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_sp_dp

  SUBROUTINE ri_part2whole_1d_dp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:)
    REAL(sp), INTENT(inout) :: whole_data(:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = REAL(part_data(i), sp)
    END DO
  END SUBROUTINE ri_part2whole_1d_dp_sp

  SUBROUTINE ri_part2whole_1d_i4_i4(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    INTEGER(i4), INTENT(in) :: part_data(:)
    INTEGER(i4), INTENT(inout) :: whole_data(:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, n, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    DO i = 1, n
      whole_data(ri%reorder_index(ofs + i)) = part_data(i)
    END DO
  END SUBROUTINE ri_part2whole_1d_i4_i4

  SUBROUTINE ri_part2whole_2d_dp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:,:)
    REAL(dp), INTENT(inout) :: whole_data(:,:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
!$omp do
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_dp_dp

  SUBROUTINE ri_part2whole_2d_sp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:,:)
    REAL(sp), INTENT(inout) :: whole_data(:,:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
!$omp do
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_sp_sp

  SUBROUTINE ri_part2whole_2d_sp_dp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(sp), INTENT(in) :: part_data(:,:)
    REAL(dp), INTENT(inout) :: whole_data(:,:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
!$omp do
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_sp_dp

  SUBROUTINE ri_part2whole_2d_dp_sp(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    REAL(dp), INTENT(in) :: part_data(:,:)
    REAL(sp), INTENT(inout) :: whole_data(:,:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
!$omp do
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = REAL(part_data(i, k), sp)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_dp_sp

  SUBROUTINE ri_part2whole_2d_i4_i4(ri, part_idx, part_data, whole_data)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER, INTENT(in) :: part_idx
    INTEGER(i4), INTENT(in) :: part_data(:,:)
    INTEGER(i4), INTENT(inout) :: whole_data(:,:)
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, whole_data
#endif
    INTEGER :: i, k, n, nlev, ofs

    ofs = ri%pe_off(part_idx)
    n = ri%pe_own(part_idx)
    nlev = SIZE(part_data, 2)
!$omp do
    DO k = 1, nlev
      DO i = 1, n
        whole_data(ri%reorder_index(ofs + i), k) = part_data(i, k)
      END DO
    END DO
  END SUBROUTINE ri_part2whole_2d_i4_i4

  SUBROUTINE ri_blk2part_2d_dp_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(dp), INTENT(in) :: blk_data(:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
!$omp do
    DO i = 1, n
      part_data(offset + i) = blk_data(ri%own_idx(i), ri%own_blk(i))
    END DO
!$omp end do nowait
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_dp_dp

  SUBROUTINE ri_blk2part_2d_sp_sp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(sp), INTENT(in) :: blk_data(:,:)
    REAL(sp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
!$omp do
    DO i = 1, n
      part_data(offset + i) = blk_data(ri%own_idx(i), ri%own_blk(i))
    END DO
!$omp end do nowait
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_sp_sp

  SUBROUTINE ri_blk2part_2d_sp_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(sp), INTENT(in) :: blk_data(:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n

    n = ri%n_own
!$omp do
    DO i = 1, n
      part_data(offset + i) = REAL(blk_data(ri%own_idx(i), ri%own_blk(i)), dp)
    END DO
!$omp end do nowait
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_sp_dp

  SUBROUTINE ri_blk2part_2d_i4_i4(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER(i4), INTENT(in) :: blk_data(:,:)
    INTEGER(i4), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
!$omp do
    DO i = 1, n
      part_data(offset + i) = blk_data(ri%own_idx(i), ri%own_blk(i))
    END DO
!$omp end do nowait
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_i4_i4

  SUBROUTINE ri_blk2part_2d_i4_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER(i4), INTENT(in) :: blk_data(:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, n
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
!$omp do
    DO i = 1, n
      part_data(offset + i) = REAL(blk_data(ri%own_idx(i), ri%own_blk(i)), dp)
    END DO
!$omp end do nowait
    offset = offset + n
  END SUBROUTINE ri_blk2part_2d_i4_dp

  SUBROUTINE ri_blk2part_3d_dp_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(dp), INTENT(in) :: blk_data(:,:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, j, n, nlev, ofs
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
    nlev = SIZE(blk_data, 2)
!$omp do
    DO j = 1, nlev
      ofs = offset + (j - 1) * n
      DO i = 1, n
        part_data(ofs + i) = blk_data(ri%own_idx(i), j, ri%own_blk(i))
      END DO
    END DO
!$omp end do nowait
    offset = offset + nlev * n
  END SUBROUTINE ri_blk2part_3d_dp_dp

  SUBROUTINE ri_blk2part_3d_sp_sp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(sp), INTENT(in) :: blk_data(:,:,:)
    REAL(sp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, j, n, nlev, ofs
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
    nlev = SIZE(blk_data, 2)
!$omp do
    DO j = 1, nlev
      ofs = offset + (j - 1) * n
      DO i = 1, n
        part_data(ofs + i) = blk_data(ri%own_idx(i), j, ri%own_blk(i))
      END DO
    END DO
!$omp end do nowait
    offset = offset + n * nlev
  END SUBROUTINE ri_blk2part_3d_sp_sp

  SUBROUTINE ri_blk2part_3d_sp_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    REAL(sp), INTENT(in) :: blk_data(:,:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, j, n, nlev, ofs
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
    nlev = SIZE(blk_data, 2)
!$omp do
    DO j = 1, nlev
      ofs = offset + (j - 1) * n
      DO i = 1, n
        part_data(ofs + i) = REAL(blk_data(ri%own_idx(i), j, ri%own_blk(i)), dp)
      END DO
    END DO
!$omp end do nowait
    offset = offset + n * nlev
  END SUBROUTINE ri_blk2part_3d_sp_dp

  SUBROUTINE ri_blk2part_3d_i4_i4(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER(i4), INTENT(in) :: blk_data(:,:,:)
    INTEGER(i4), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, j, n, nlev, ofs
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
    nlev = SIZE(blk_data, 2)
!$omp do
    DO j = 1, nlev
      ofs = offset + (j - 1) * n
      DO i = 1, n
        part_data(ofs + i) = blk_data(ri%own_idx(i), j, ri%own_blk(i))
      END DO
    END DO
!$omp end do nowait
    offset = offset + n * nlev
  END SUBROUTINE ri_blk2part_3d_i4_i4

  SUBROUTINE ri_blk2part_3d_i4_dp(ri, blk_data, part_data, offset)
    TYPE(t_reorder_info), INTENT(IN) :: ri
    INTEGER(i4), INTENT(in) :: blk_data(:,:,:)
    REAL(dp), INTENT(inout) :: part_data(:)
    INTEGER, INTENT(inout) :: offset
    INTEGER :: i, j, n, nlev, ofs
#ifdef USE_CONTIGUOUS
    CONTIGUOUS :: part_data, blk_data
#endif

    n = ri%n_own
    nlev = SIZE(blk_data, 2)
!$omp do
    DO j = 1, nlev
      ofs = offset + (j - 1) * n
      DO i = 1, n
        part_data(ofs + i) = REAL(blk_data(ri%own_idx(i), j, ri%own_blk(i)), dp)
      END DO
    END DO
!$omp end do nowait
    offset = offset + n * nlev
  END SUBROUTINE ri_blk2part_3d_i4_dp

END MODULE mo_reorder_info
