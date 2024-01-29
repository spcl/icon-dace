! Module for reading multifile restart files
!
! Note: The single file implementation of the restart input can be
!       found in the module "mo_load_singlefile_restart"
!
!
! All work processes are used to read the data from the
! multifile. Each processor handles all domains of its horizontal
! "chunk", since the domain loop is located outside of this routine,
! in the calling "src/drivers" routine.
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
#include "omp_definitions.inc"

MODULE mo_load_multifile_restart

  USE ISO_C_BINDING,             ONLY: C_PTR, C_LOC, C_F_POINTER
  USE mo_c_restart_util,         ONLY: checkMultifileDir
  USE mo_impl_constants,         ONLY: SUCCESS, SINGLE_T, REAL_T, INT_T
  USE mo_communication,          ONLY: t_comm_pattern, t_p_comm_pattern, exchange_data
  USE mo_communication_factory,  ONLY: setup_comm_pattern
  USE mo_decomposition_tools,    ONLY: t_glb2loc_index_lookup, init_glb2loc_index_lookup, set_inner_glb_index, &
    &                                  deallocate_glb2loc_index_lookup
  USE mo_dynamics_config,        ONLY: nnow, nnow_rcf
  USE mo_exception,              ONLY: finish, warning
  USE mo_kind,                   ONLY: sp, dp
  USE mo_model_domain,           ONLY: t_patch
  USE mo_mpi,                    ONLY: p_comm_work, p_n_work, p_pe_work, my_process_is_work, p_bcast, &
    &                                  p_sum, p_alltoall, p_alltoallv
  USE mo_multifile_restart_util, ONLY: multifilePayloadPath, rBuddy, rGroup, vNames_glbIdx
  USE mo_parallel_config,        ONLY: nproma, idx_no, blk_no, restart_load_scale_max
  USE mo_restart_nml_and_att,    ONLY: getAttributesForRestarting, ocean_initFromRestart_OVERRIDE
  USE mo_key_value_store,        ONLY: t_key_value_store
  USE mo_restart_var_data,       ONLY: get_var_3d_ptr, has_valid_time_level
  USE mo_var,                    ONLY: t_var_ptr
  USE mo_timer,                  ONLY: timer_start, timer_stop, timer_load_restart_io, timers_level, &
    &                                  timer_load_restart_comm_setup, timer_load_restart_communication, &
    &                                  timer_load_restart_get_var_id
  USE mo_netcdf_errhandler,      ONLY: nf
  USE mo_netcdf
  USE mo_cdi_constants,      ONLY: GRID_UNSTRUCTURED_CELL, GRID_UNSTRUCTURED_VERT, GRID_UNSTRUCTURED_EDGE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: multifileCheckRestartFiles, multifileReadPatch

  TYPE t_PayloadFile
    INTEGER :: ncid, iCnts(3), iVarIds(3)
  END TYPE t_PayloadFile

  CHARACTER(*), PARAMETER :: modname = "mo_load_multifile_restart"

CONTAINS

  SUBROUTINE multifileCheckRestartFiles(filename)
    CHARACTER(*), INTENT(IN) :: filename
    TYPE(t_key_value_store), POINTER :: restartAttributes
    INTEGER :: n_dom, multifile_file_count
    CHARACTER(*), PARAMETER :: routine = modname//":multifileCheckRestartFiles"

    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('multifile_n_dom_active', n_dom)
    CALL restartAttributes%get('multifile_file_count', multifile_file_count)
    IF(checkMultifileDir(filename, n_dom, multifile_file_count) /= SUCCESS) &
      CALL finish(routine, "'"//filename//"' is not a valid restart multifile")
  END SUBROUTINE multifileCheckRestartFiles

  ! Opens all the payload files that should be opened by this
  ! process, creating t_PayloadFile objects for them, AND
  ! calculating the total numbers of cells/edges/vertices contained
  ! IN these files.
  !
  ! This FUNCTION performs the distribution of the payload files
  ! among the restart processes, subsequent code just has to iterate
  ! over the array returned by this FUNCTION.
  SUBROUTINE openPayloadFiles(mfPath, mfCnt, dom, tCnt, files, at_once)
    CHARACTER(*), INTENT(IN) :: mfPath
    INTEGER, INTENT(IN) :: mfCnt, dom
    INTEGER, INTENT(OUT) :: tCnt(3)
    TYPE(t_PayloadFile), INTENT(INOUT), ALLOCATABLE :: files(:)
    LOGICAL, INTENT(OUT) :: at_once
    INTEGER :: myFCnt, myFF, i
    CHARACTER(*), PARAMETER :: routine = modname//":openPayloadFiles"

    IF(.NOT.my_process_is_work()) CALL finish(routine, "this is not a work process")
! number of reader ranks buffering one full 3d field of at max restart_load_scale_max worker ranks
! default : restart_load_scale_max = 1
! if restart_load_scale_max == 0 - - - - - - - -> read data in a LEVEL-based loop
! if restart_load_scale_max <  #ranks / #files -> read data in a LEVEL-based loop
! if restart_load_scale_max >= #ranks / #files -> read data in a VARIABLE-based loop
! if restart_load_scale_max < 0  - - - - - - - -> read data in a VARIABLE-based loop
    at_once = p_n_work .LE. mfCnt * &
      & MERGE(restart_load_scale_max, p_n_work, restart_load_scale_max .GE. 0)
    myFF = -1
    IF (p_n_work .GE. mfCnt) THEN
      myFCnt = MERGE(1, 0, rBuddy(nr_in=mfCnt) .EQ. p_pe_work)
      IF (myFCnt .GT. 0) myFF = rGroup(nr_in=mfCnt)
    ELSE
      myFCnt = COUNT( (/ (rGroup(nr_in=p_n_work, nw_in=mfCnt, pe_in=i) .EQ. p_pe_work, &
        &                    i = p_pe_work, mfCnt-1) /) )
      i = p_pe_work
      DO WHILE(myFF .LT. 0)
        IF (rGroup(nr_in=p_n_work, nw_in=mfCnt, pe_in=i) .EQ. p_pe_work) myFF = i
        i = i + 1
      END DO
    END IF
#ifdef DEBUG
    WRITE(0, "(3(a,i4),a)") "restart: rank ", p_pe_work, " will start reading from file ", &
      &                      myFF, ", ", myFCnt, " files in total"
#endif
    ALLOCATE(files(myFCnt))
    tCnt(:) = 0
    DO i = 1, myFCnt
      CALL payloadfile_open(files(i), myFF - 1 + i)
      tCnt(:) = tCnt(:) + files(i)%iCnts(:)
    END DO
  CONTAINS

    SUBROUTINE payloadFile_open(me, partId)
      TYPE(t_PayloadFile), INTENT(INOUT) :: me
      INTEGER, INTENT(in) :: partId
      INTEGER :: iG, ndim, dimid(4)
      CHARACTER(*), PARAMETER :: routine = modname//":payloadFile_init"
      CHARACTER(:), ALLOCATABLE :: pathname
  
      IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
      CALL multifilePayloadPath(mfPath, dom, partId, pathname)
      CALL nf(nf_open(pathname, NF_NOWRITE, me%ncid), routine)
      DO iG = 1, 3
        CALL nf(nf_inq_varid(me%ncid, vNames_glbIdx(iG), me%iVarIds(iG)), routine)
        CALL nf(nf_inq_varndims(me%ncid, me%iVarIds(iG), ndim), routine)
        CALL nf(nf_inq_vardimid(me%ncid, me%iVarIds(iG), dimid(1:ndim)), routine)
        CALL nf(nf_inq_dimlen(me%ncid, dimid(1), me%iCnts(iG)), routine)
      END DO
      IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
    END SUBROUTINE payloadFile_open
  END SUBROUTINE openPayloadFiles

  ! This asserts that the given pattern redistributes a given array
  ! exactly as described by the global indices arrays.
  SUBROUTINE checkRedistributionPattern(pat, provGlbIdces, reqdGlbIdces)
    CLASS(t_comm_pattern), INTENT(IN), POINTER :: pat
    INTEGER, INTENT(IN) :: provGlbIdces(:), reqdGlbIdces(:)
    INTEGER :: iSize, oSize, i, error 
    CHARACTER(*), PARAMETER :: routine = modname//":checkRedistributionPattern"
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: input, output

    iSize = SIZE(provGlbIdces)
    oSize = SIZE(reqdGlbIdces)
    ALLOCATE(input(nproma, blk_no(iSize)), output(nproma, blk_no(oSize)))
!ICON_OMP PARALLEL PRIVATE(i)
!ICON_OMP DO SCHEDULE(STATIC)
    DO i = 1, blk_no(oSize)
      output(:,i) = REAL(-1, dp)
    END DO
!ICON_OMP DO SCHEDULE(STATIC)
    DO i = 1, iSize
      input(idx_no(i), blk_no(i)) = REAL(provGlbIdces(i), dp)
    END DO
!ICON_OMP END PARALLEL
    CALL exchange_data(pat, output, input)
    error = 0
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, oSize
      IF(INT(output(idx_no(i), blk_no(i))) /= reqdGlbIdces(i)) THEN
!ICON_OMP CRITICAL
        error = error + 1
        PRINT "(6(i6,1x))", i, blk_no(i), idx_no(i), oSize, &
          & INT(output(idx_no(i), blk_no(i))), reqdGlbIdces(i)
!ICON_OMP END CRITICAL
      END IF
    END DO
    IF (error .GT. 0) THEN
      PRINT "(2(a,i8))", "check comm pattern: found deviations: ", &
        & error, " out of ", oSize
      CALL finish(routine, "comm pattern does not work as expected")
    END IF
  END SUBROUTINE checkRedistributionPattern

  ! Create a t_comm_pattern that IS able to redistribute the DATA as
  ! it IS READ from the multifile to the domain decomposition of the
  ! work processes.  The tricky part about this IS, that
  ! setup_comm_pattern() requires an owner(:) array as input, which
  ! we have to create first.
  !
  ! The resulting t_comm_pattern will also initialize the halo points.
  FUNCTION makeRedistributionPattern(pIds, rIds) RESULT(pat)
    CLASS(t_comm_pattern), POINTER :: pat
    INTEGER, INTENT(IN) :: pIds(:), rIds(:)
    INTEGER :: glbSize, locSize, cuml, ii
    INTEGER, ALLOCATABLE :: own_src(:), own_dst(:), bd(:)
    TYPE(t_glb2loc_index_lookup) :: lupTbl
    CHARACTER(*), PARAMETER :: routine = modname//":makeRedistributionPattern"

    IF(timers_level >= 7) CALL timer_start(timer_load_restart_comm_setup)
    locSize = SIZE(pIds)
    glbSize = p_sum(locSize, comm=p_comm_work)
    CALL get_owners()
    CALL init_glb2loc_index_lookup(lupTbl, glbSize)
    CALL set_inner_glb_index(lupTbl, pIds, [(ii, ii = 1, locSize)])
    CALL setup_comm_pattern(SIZE(rIds), own_dst, rIds, lupTbl, locSize, &
      & own_src, pIds, pat, inplace=.FALSE., comm=p_comm_work)
    CALL deallocate_glb2loc_index_lookup(lupTbl)
#ifdef DEBUG
    CALL checkRedistributionPattern(pat, pGlbIds, rGlbIds)
#endif
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_comm_setup)
  CONTAINS

    SUBROUTINE get_owners()
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dI_s, dI_c, dI_i, dI_o, o_m, &
        & iPk_s, iPk_c, iPk_i, iPk_o, nI_s, nI_c, nI_i, nI_o, oPk_i, oPk_o, oPk_c
      INTEGER :: i, j

#ifdef DEBUG 
      ! consistency checks
      IF(1 > MINVAL(pIds,1)) &
        & CALL finish(routine, "global index out of RANGE (<1)")
      IF(glbSize < MAXVAL(pIds,1)) &
        &  CALL finish(routine, "global index out of range (>globalSize)")
#endif
      ALLOCATE(bd(0:p_n_work), nI_i(p_n_work), nI_o(p_n_work))
      bd(:) = [(INT(REAL(i,dp)*(REAL(glbSize,dp)/REAL(p_n_work,dp))), i=0,p_n_work)]
      bd(p_n_work) = glbSize
      CALL get_disps(dI_s, nI_s, iPk_s, pIds)
      CALL get_disps(dI_c, nI_c, iPk_c, rIds)
      CALL p_alltoall(nI_s, nI_i, p_comm_work)
      CALL p_alltoall(nI_c, nI_o, p_comm_work)
      ALLOCATE(dI_o(p_n_work+1), dI_i(p_n_work+1))
      dI_o(1) = 0
      dI_i(1) = 0
      DO i = 2, p_n_work+1
        dI_o(i) = dI_o(i-1) + nI_o(i-1)
        dI_i(i) = dI_i(i-1) + nI_i(i-1)
      END DO
      ALLOCATE(oPk_i(dI_i(p_n_work+1)), o_m(bd(p_pe_work)+1:bd(p_pe_work+1)), &
        & oPk_o(dI_o(p_n_work+1)), iPk_o(dI_o(p_n_work+1)), &
        & oPk_c(dI_c(p_n_work+1)), iPk_i(dI_i(p_n_work+1)))
      CALL p_alltoallv(iPk_s, nI_s, dI_s, iPk_i, nI_i, dI_i, p_comm_work)
      CALL p_alltoallv(iPk_c, nI_c, dI_c, iPk_o, nI_o, dI_o, p_comm_work)
      DEALLOCATE(dI_s, nI_s, nI_i, iPk_s, iPk_c)
!ICON_OMP PARALLEL PRIVATE(i)
!ICON_OMP DO SCHEDULE(GUIDED)
      DO i = 1, p_n_work
        IF (dI_i(i)+1 .LE. dI_i(i+1)) &
          & oPk_i(dI_i(i)+1:dI_i(i+1)) = i - 1
      END DO
!ICON_OMP DO SCHEDULE(GUIDED)
      DO i = 1, dI_i(p_n_work+1)
        o_m(iPk_i(i)) = oPk_i(i)
      END DO
!ICON_OMP DO SCHEDULE(GUIDED)
      DO i = 1, dI_o(p_n_work+1)
        oPk_o(i) = o_m(iPk_o(i))
      END DO
!ICON_OMP END PARALLEL
      CALL p_alltoallv(oPk_o, nI_o, dI_o, oPk_c, nI_c, dI_c, p_comm_work)
      DEALLOCATE(o_m, oPk_o, iPk_o, iPk_i, dI_i, nI_o, dI_o)
      ALLOCATE(own_dst(SIZE(rIds)), own_src(locSize))
!ICON_OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(cuml,j)
      DO i = 1, p_n_work
        IF (dI_c(i+1)-dI_c(i) .GT. 0) THEN
          cuml = dI_c(i)
          DO j = 1, SIZE(rIds)
            IF (rIds(j) .GT. bd(i-1) .AND. rIds(j) .LE. bd(i)) THEN
              cuml = cuml + 1
              own_dst(j) = oPk_c(cuml)
            END IF
          END DO
        END IF
      END DO
      own_src(:) = p_pe_work
    END SUBROUTINE get_owners

    SUBROUTINE get_disps(dId, nId, idPk, id)
      INTEGER, INTENT(INOUT), ALLOCATABLE :: dId(:), nId(:), idPk(:)
      INTEGER, INTENT(IN), CONTIGUOUS :: id(:)
      INTEGER :: i, j

      ALLOCATE(dId(p_n_work+1), idPk(SIZE(id)), nId(p_n_work))
      cuml = 0
      DO i = 1, p_n_work
        dId(i) = cuml
        DO j = 1, SIZE(id)
          IF (id(j) .GT. bd(i-1) .AND. id(j) .LE. bd(i)) THEN
            cuml = cuml + 1
            idPk(cuml) = id(j)
          END IF
        END DO
        nId(i) = cuml - dId(i)
      END DO
      dId(p_n_work+1) = cuml
    END SUBROUTINE get_disps
  END FUNCTION makeRedistributionPattern

  SUBROUTINE multifileReadPatch(vDat, ptc, multifilePath)
    TYPE(t_var_ptr), INTENT(IN) :: vDat(:)
    TYPE(t_patch), TARGET, INTENT(IN) :: ptc
    CHARACTER(*), INTENT(IN) :: multifilePath
    TYPE(t_p_comm_pattern) :: cpat(3)
    LOGICAL :: load_var_at_once, int_is_int
    INTEGER :: hi, hmap(MAX(GRID_UNSTRUCTURED_CELL,GRID_UNSTRUCTURED_VERT,GRID_UNSTRUCTURED_EDGE))
    TYPE(t_PayloadFile), ALLOCATABLE :: files(:)

    hmap(GRID_UNSTRUCTURED_CELL) = 1
    hmap(GRID_UNSTRUCTURED_VERT) = 2
    hmap(GRID_UNSTRUCTURED_EDGE) = 3
    CALL construct()
    CALL readData()
    DO hi = 1, 3
      CALL cpat(hi)%p%delete()
    END DO
    IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
    IF(ALLOCATED(files)) THEN
      DO hi = 1, SIZE(files)
        CALL nf(nf_close(files(hi)%ncid), modname//":multifileReadPatch")
      END DO
    END IF
    IF(timers_level >= 7) CALL timer_stop(timer_load_restart_io)
  CONTAINS

  SUBROUTINE construct()
    INTEGER :: mfileCnt, ierr, cFId, tCnt(3), cOff(3), iG, i, n
    INTEGER, ALLOCATABLE :: glbidx_read(:)
    INTEGER, POINTER :: glb_index(:)
    REAL(dp), ALLOCATABLE :: buffer(:)
    TYPE(t_key_value_store), POINTER :: restartAttributes
    CHARACTER(*), PARAMETER :: routine = modname//":multifileReadPatch:construct"

    CALL getAttributesForRestarting(restartAttributes)
    CALL restartAttributes%get('multifile_file_count', mfileCnt)
    IF (p_pe_work .EQ. 0) WRITE(0, *) "reading from ", mfileCnt, " files/patch."
    CALL restartAttributes%get('int_is_int', int_is_int, opt_err=ierr)
    IF (ierr .NE. 0) int_is_int = .FALSE.
    CALL openPayloadFiles(multifilePath, mfileCnt, ptc%id, tCnt, files, load_var_at_once)
    cOff(:) = 0
    ALLOCATE(glbidx_read(MAXVAL(tCnt(:),1)))
    DO iG = 1, 3
      DO cFId = 1, SIZE(files)
        n = files(cFId)%iCnts(iG)
        IF (n .LE. 0) CYCLE
        IF (timers_level >= 7) CALL timer_start(timer_load_restart_io)
        IF (int_is_int) THEN
          CALL nf(nf_get_vara_int(files(cFId)%ncid, files(cFId)%iVarIds(iG), [1,1], [n,1], &
            &                     glbidx_read(cOff(iG)+1:cOff(iG)+n)), routine)
        ELSE
          ALLOCATE(buffer(n))
          CALL nf(nf_get_vara_double(files(cFId)%ncid, files(cFId)%iVarIds(iG), [1,1], [n,1], &
            &                        buffer), routine)
          !ICON_OMP PARALLEL DO SCHEDULE(STATIC)
          DO i = 1, n
            glbidx_read(cOff(iG) + i) = INT(buffer(i))
          END DO
          DEALLOCATE(buffer)
        END IF
        IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
        cOff(iG) = cOff(iG) + n
      END DO
      SELECT CASE (iG)
      CASE (GRID_UNSTRUCTURED_CELL)
        glb_index => ptc%cells%decomp_info%glb_index
      CASE (GRID_UNSTRUCTURED_VERT)
        glb_index => ptc%verts%decomp_info%glb_index
      CASE (GRID_UNSTRUCTURED_EDGE)
        glb_index => ptc%edges%decomp_info%glb_index
      END SELECT
      cpat(ig)%p => makeRedistributionPattern(glbidx_read(:cOff(iG)), glb_index)
    END DO
  END SUBROUTINE construct

  SUBROUTINE readData()
    INTEGER :: iV, lId, llId, lCnt, fId, vIDs(SIZE(files), SIZE(vDat)), st(3), ct(3)
    INTEGER :: ofs, i, hgrid, nblk, nds(SIZE(vDat)), pCt(SIZE(files)), max_r, max_e
    REAL(dp), ALLOCATABLE, TARGET :: buf_e(:), buf_r(:)
    REAL(dp), POINTER :: ptr_3d_d(:,:,:)
    REAL(dp), CONTIGUOUS, POINTER :: buf_3d_d(:,:,:), buf_d(:)
    REAL(sp), POINTER :: ptr_3d_s(:,:,:)
    REAL(sp), CONTIGUOUS, POINTER :: buf_3d_s(:,:,:), buf_s(:)
    INTEGER, POINTER :: ptr_3d_i(:,:,:)
    INTEGER, CONTIGUOUS, POINTER :: buf_3d_i(:,:,:), buf_i(:)
    CHARACTER(*), PARAMETER :: routine = modname//":multifileReadPatch:readData"
    LOGICAL :: en_bloc
    TYPE(C_PTR) :: cptr_r, cptr_e

    nds(:) = -1; max_r = 1; max_e = 1
    IF (timers_level >= 7) CALL timer_start(timer_load_restart_get_var_id)
    IF (ALLOCATED(files)) THEN
      IF (SIZE(files) .GT. 0) THEN
        DO iV = 1, SIZE(vDat)
          IF (.NOT.has_valid_time_level(vDat(iV)%p%info, ptc%id, nnow(ptc%id), nnow_rcf(ptc%id))) CYCLE
          i = nf_inq_varid(files(1)%ncid, vDat(iV)%p%info%name, vIDs(1,iV))
          IF (i .EQ. NF_NOERR) THEN
            CALL nf(nf_inq_varndims(files(1)%ncid, vIDs(1, iV), nds(iV)), routine)
            DO fId = 2, SIZE(files)
              CALL nf(nf_inq_varid(files(fId)%ncid, vDat(iV)%p%info%name, vIDs(fId,iV)), routine)
            END DO
            lCnt = MERGE(vDat(iV)%p%info%used_dimensions(2), 1, vDat(iV)%p%info%ndims .GT. 2)
            en_bloc = load_var_at_once .AND. lCnt .GT. 1
            hgrid = hmap(vDat(iV)%p%info%hgrid)
            nblk = blk_no(SUM(files(:)%iCnts(hgrid)))*nproma
            max_r = MAX(max_r, MERGE(MAXVAL(files(:)%iCnts(hgrid),1)*lCnt, nblk, en_bloc))
            max_e = MAX(max_e, MERGE(nblk*lCnt, 1, en_bloc))
          ELSE
            IF (ocean_initFromRestart_OVERRIDE .OR. vDat(iV)%p%info%lrestart_cont) THEN
              CALL warning(routine, "variable not found: "//TRIM(vDat(iV)%p%info%NAME))
            ELSE
              CALL finish(routine, "variable not found: "//TRIM(vDat(iV)%p%info%NAME))
            END IF
          END IF
        END DO
      END IF
    END IF
    IF (timers_level >= 7) CALL timer_stop(timer_load_restart_get_var_id)
    CALL p_bcast(nds, 0, comm=p_comm_work)
    ALLOCATE(buf_r(max_r), buf_e(max_e))
    cptr_r = C_LOC(buf_r(1))
    cptr_e = C_LOC(buf_e(1))
    DO iV = 1, SIZE(vDat)
      IF (nds(iV) .EQ. -1) CYCLE
      hgrid = hmap(vDat(iV)%p%info%hgrid)
      pCt = files(:)%iCnts(hgrid)
      nblk = blk_no(SUM(pCt(:)))
      lCnt = MERGE(vDat(iV)%p%info%used_dimensions(2), 1, vDat(iV)%p%info%ndims .GT. 2)
      en_bloc = load_var_at_once .AND. lCnt .GT. 1
      SELECT CASE(vDat(iV)%p%info%data_type)
      CASE(REAL_T)
        CALL get_var_3d_ptr(vDat(iV)%p, ptr_3d_d)
        CALL C_F_POINTER(cptr_r, buf_d, [MERGE(MAXVAL(pCt,1)*lCnt, nblk*nproma, en_bloc)])
        IF (en_bloc) CALL C_F_POINTER(cptr_e, buf_3d_d, [nproma,lcnt,nblk])
        DO llId = 1, MERGE(1, lCnt, en_bloc)
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
          ofs = 0
          st(:) = [1, llId, 1]
          DO fId = 1, SIZE(files)
            IF (pCt(fId) .LT. 1) CYCLE
            ct(:) = [pCt(fId), MERGE(lCnt, 1, en_bloc), 1]
            CALL nf(nf_get_vara_double(files(fId)%ncid, vIDs(fId, iV), st(:nds(iV)), ct(:nds(iV)), &
              & buf_d(MERGE(1, ofs+1, en_bloc):MERGE(ct(1)*lcnt, ofs+ct(1), en_bloc))), routine)
            IF (en_bloc) THEN
              !ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
              DO lId = 1, lCnt
                DO i = 1, pCt(fId)
                  buf_3d_d(idx_no(i+ofs), lId, blk_no(i+ofs)) = buf_d(i+(lId-1)*pCt(fId))
                END DO
              END DO
            END IF
            ofs = ofs + ct(1)
          END DO
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
          IF (.NOT.en_bloc) CALL C_F_POINTER(cptr_r, buf_3d_d, [nproma,1,nblk])
          CALL exchange_data(cpat(hgrid)%p, ptr_3d_d(:,llId:MERGE(lCnt, llId, en_bloc),:), buf_3d_d)
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
        END DO
      CASE(SINGLE_T)
        CALL get_var_3d_ptr(vDat(iV)%p, ptr_3d_s)
        CALL C_F_POINTER(cptr_r, buf_s, [MERGE(MAXVAL(pCt,1)*lCnt, nblk*nproma, en_bloc)])
        IF (en_bloc) CALL C_F_POINTER(cptr_e, buf_3d_s, [nproma,lcnt,nblk])
        DO llId = 1, MERGE(1, lCnt, en_bloc)
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
          ofs = 0
          st(:) = [1, llId, 1]
          DO fId = 1, SIZE(files)
            IF (pCt(fId) .LT. 1) CYCLE
            ct(:) = [pCt(fId), MERGE(lCnt, 1, en_bloc), 1]
            CALL nf(nf_get_vara_real(files(fId)%ncid, vIDs(fId, iV), st(:nds(iV)), ct(:nds(iV)), &
              & buf_s(MERGE(1, ofs+1, en_bloc):MERGE(ct(1)*lcnt, ofs+ct(1), en_bloc))), routine)
            IF (en_bloc) THEN
              !ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
              DO lId = 1, lCnt
                DO i = 1, pCt(fId)
                  buf_3d_s(idx_no(i+ofs), lId, blk_no(i+ofs)) = buf_s(i+(lId-1)*pCt(fId))
                END DO
              END DO
            END IF
            ofs = ofs + ct(1)
          END DO
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
          IF (.NOT.en_bloc) CALL C_F_POINTER(cptr_r, buf_3d_s, [nproma,1,nblk])
          CALL exchange_data(cpat(hgrid)%p, ptr_3d_s(:,llId:MERGE(lCnt, llId, en_bloc),:), buf_3d_s)
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
        END DO
      CASE(INT_T)
        CALL get_var_3d_ptr(vDat(iV)%p, ptr_3d_i)
        IF (int_is_int) THEN
          CALL C_F_POINTER(cptr_r, buf_i, [MERGE(MAXVAL(pCt,1)*lCnt, nblk*nproma, en_bloc)])
          IF (en_bloc) CALL C_F_POINTER(cptr_e, buf_3d_i, [nproma,MERGE(lcnt,1,en_bloc),nblk])
        ELSE
          CALL C_F_POINTER(cptr_r, buf_d, [MERGE(MAXVAL(pCt,1)*lCnt, nblk*nproma, en_bloc)])
          CALL C_F_POINTER(cptr_e, buf_3d_i, [nproma,MERGE(lcnt,1,en_bloc),nblk])
        END IF
        DO llId = 1, MERGE(1, lCnt, en_bloc)
          ofs = 0 
          st(:) = [1, llId, 1]
          IF(timers_level >= 7) CALL timer_start(timer_load_restart_io)
          DO fId = 1, SIZE(files)
            IF (pCt(fId) .LT. 1) CYCLE
            ct(:) = [pCt(fId), MERGE(lCnt, 1, en_bloc), 1]
            IF (int_is_int) THEN
              CALL nf(nf_get_vara_int(files(fId)%ncid, vIDs(fId, iV), st(:nds(iV)), ct(:nds(iV)), &
                & buf_i(MERGE(1, ofs+1, en_bloc):MERGE(ct(1)*lcnt, ofs+ct(1), en_bloc))), routine)
              IF (en_bloc) THEN
                !ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
                DO lId = 1, lCnt
                  DO i = 1, pCt(fId)
                    buf_3d_i(idx_no(i+ofs), lId, blk_no(i+ofs)) = buf_i(i+(lId-1)*pCt(fId))
                  END DO
                END DO
              END IF
            ELSE
              CALL nf(nf_get_vara_double(files(fId)%ncid, vIDs(fId, iV), st(:nds(iV)), ct(:nds(iV)), &
                & buf_d(1:ct(1)*MERGE(lcnt, 1, en_bloc))), routine)
              !ICON_OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2)
              DO lId = 1, MERGE(lCnt, 1, en_bloc)
                DO i = 1, pCt(fId)
                  buf_3d_i(idx_no(i+ofs), lId, blk_no(i+ofs)) = INT(buf_d(i+(lId-1)*pCt(fId)))
                END DO
              END DO
            END IF
            ofs = ofs + ct(1)
          END DO
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_io)
          IF (timers_level >= 7) CALL timer_start(timer_load_restart_communication)
          IF (.NOT.en_bloc .AND. int_is_int) CALL C_F_POINTER(cptr_r, buf_3d_i, [nproma,1,nblk])
          CALL exchange_data(cpat(hgrid)%p, ptr_3d_i(:,llId:MERGE(lCnt, llId, en_bloc),:), buf_3d_i)
          IF (timers_level >= 7) CALL timer_stop(timer_load_restart_communication)
        END DO
      END SELECT
    END DO
  END SUBROUTINE readData

  END SUBROUTINE multifileReadPatch

END MODULE mo_load_multifile_restart
