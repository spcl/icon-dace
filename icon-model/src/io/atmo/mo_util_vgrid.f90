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

! Module handling the reading / on-the-fly generation and the writing of the vertical grid.

MODULE mo_util_vgrid

  USE ISO_C_BINDING,                        ONLY: C_DOUBLE
  USE mo_cdi,                               ONLY: streamDefTimestep, streamOpenWrite, gridCreate, institutInq, vlistCreate, &
                                                & vlistInqVarZaxis, vlistInqVarGrid, streamInqVlist, streamOpenRead, &
                                                & ZAXIS_REFERENCE, zaxisCreate, TIME_CONSTANT, vlistDefVar, FILETYPE_NC2, &
                                                & DATATYPE_INT32, DATATYPE_FLT64, CDI_UNDEFID, CDI_GLOBAL, cdiDefAttInt, &
                                                & cdiInqAttInt, zaxisDestroy, gridDestroy, vlistDestroy, streamClose, &
                                                & streamWriteVarSlice, streamWriteVar, streamDefVlist, vlistDefVarTsteptype, &
                                                & vlistDefVarDatatype, vlistDefVarName, zaxisDefNumber, zaxisDefUUID, &
                                                & gridDefPosition, gridInqUUID, gridDefNumber, gridDefUUID, zaxisDefLevels, &
                                                & gridDefNvertex, vlistDefInstitut, zaxisInqUUID, streamReadVar, &
                                                & GRID_UNSTRUCTURED, TSTEP_CONSTANT
  USE mo_kind,                              ONLY: wp, dp
  USE mo_exception,                         ONLY: finish, message, message_text, warning
  !
  USE mo_grid_config,                       ONLY: n_dom, vertical_grid_filename, vct_filename, create_vgrid
  USE mo_sleve_config,                      ONLY: lread_smt
  USE mo_nonhydrostatic_config,             ONLY: ivctype
  USE mo_parallel_config,                   ONLY: nproma
  USE mo_gribout_config,                    ONLY: gribout_config
  USE mo_run_config,                        ONLY: number_of_grid_used, msg_level, check_uuid_gracefully
  !
  USE mo_impl_constants,                    ONLY: SUCCESS, MAX_CHAR_LENGTH
  USE mo_model_domain,                      ONLY: t_patch
  USE mo_ext_data_types,                    ONLY: t_external_data
  USE mo_intp_data_strc,                    ONLY: t_int_state
  USE mo_nh_testcases_nml,                  ONLY: layer_thickness, n_flat_level
  USE mo_init_vgrid,                        ONLY: init_hybrid_coord, init_sleve_coord,                  &
    &                                             prepare_vcoord, init_vert_coord
  USE mo_process_topo,                      ONLY: compute_smooth_topo
  USE mo_nh_init_nest_utils,                ONLY: topo_blending_and_fbk
  USE mo_communication,                     ONLY: exchange_data, idx_no, blk_no
  USE mo_util_string,                       ONLY: int2string
  USE mo_util_uuid_types,                   ONLY: t_uuid, UUID_STRING_LENGTH
  USE mo_util_uuid,                         ONLY: uuid_generate, OPERATOR(==), uuid_unparse
  USE mo_util_cdi,                          ONLY: get_cdi_varID, read_cdi_3d, cdiGetStringError, &
    &                                             t_inputParameters, makeInputParameters, deleteInputParameters
  USE mo_mpi,                               ONLY: my_process_is_mpi_workroot,                           &
    &                                             p_comm_work, p_bcast, p_io
  USE mo_util_vgrid_types,                  ONLY: vgrid_buffer

  IMPLICIT NONE

  PRIVATE

  ! variables and data types
  ! PUBLIC :: ...
  ! subroutines
  PUBLIC :: construct_vertical_grid

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_vgrid'


CONTAINS

  !----------------------------------------------------------------------------------------------------
  !> Import vertical grid/ define vertical coordinate
  !
  !  Note: As in the initial implementation of this code snippet, we
  !        use only the domain 1 here, and we implicitly assume that
  !        all other domains (with smaller patch%nlevs) share their
  !        levels with domain 1.
  !
  SUBROUTINE construct_vertical_grid(p_patch, p_int_state, ext_data, &
    &                                vct_a, vct_b, vct, nflatlev)
    TYPE(t_patch),          INTENT(INOUT) :: p_patch(:)
    TYPE(t_int_state),      INTENT(IN)    :: p_int_state(:)
    TYPE(t_external_data),  INTENT(INOUT) :: ext_data(:)          ! (1,..., n_dom)
    REAL(wp),               INTENT(INOUT) :: vct_a(:)             ! param. A of the vertical coordinate
    REAL(wp),               INTENT(INOUT) :: vct_b(:)             ! param. B of the vertical coordinate
    REAL(wp),               INTENT(INOUT) :: vct  (:)             ! param. A and B of the vertical coordinate
    INTEGER,                INTENT(INOUT) :: nflatlev(:)

    ! local variables
    CHARACTER(*), PARAMETER           :: routine = modname//"::construct_vertical_grid"
    INTEGER                           :: nlevp1, nblks_c, jg, error_status, j, iidx, &
      &                                  jc_c, jb_c, nlevels
    REAL(wp), ALLOCATABLE             :: topography_smt(:,:)
    REAL(C_DOUBLE), ALLOCATABLE       :: r_in(:,:)
    INTEGER, ALLOCATABLE              :: glbidx(:)
    CHARACTER(len=UUID_STRING_LENGTH) :: uuid_unparsed
    LOGICAL                           :: lvct_from_file  ! TRUE: vct_a is read from file

    !--- Initialize vertical coordinate table vct_a, vct_b (grid stretching)


    IF (TRIM(vertical_grid_filename(1)) == "") THEN

      ! skip the following paragraph if we read the 3D vertical grid information from
      ! file:
      !
      ! initialize 1D vertical coordinate tables
      IF (ivctype == 1) THEN
        CALL init_hybrid_coord(p_patch(1)%nlev, vct_filename, vct_a, vct_b, layer_thickness, n_flat_level)
      ELSE IF (ivctype == 2) THEN
        CALL init_sleve_coord(p_patch(1)%nlev, vct_filename, vct_a, vct_b)
      ENDIF

      lvct_from_file = (TRIM(vct_filename) /= "")
      CALL prepare_vcoord(lvct_from_file, ivctype, p_patch(1)%nlev, vct_a, vct_b, vct, nflatlev)

    END IF


    !--- Allocate 3D half level coordinate arrays

    ALLOCATE(vgrid_buffer(n_dom), STAT=error_status)
    IF (error_status /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    DO jg = 1,n_dom
      nlevp1   = p_patch(jg)%nlev + 1
      nblks_c  = p_patch(jg)%nblks_c

      ! Note: We allocate only the temporary field for the
      !       vertex-half-level coordinates "z_ifv". The
      !       corresponding cell-half-level coordinates "z_ifc" are
      !       a proper model variable that is created in
      !       "mo_nonhydro_state"
      !

      ! Note: This array is deallocated after being used in
      !       mo_nonhydro_state::new_nh_metrics_list
      ALLOCATE(vgrid_buffer(jg)%z_ifc(nproma,nlevp1,nblks_c), STAT=error_status)
      IF (error_status /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    END DO

    !--- initialize 3D half level coordinate

    IF (vertical_grid_filename(1) == "") THEN
      ! skip the following paragraph if we read vertical grid from
      ! file:
      !

      ! Perform topography smoothing and feedback
      IF (n_dom > 1) CALL topo_blending_and_fbk(jg=1, ext_data=ext_data(:))

      DO jg = 1,n_dom
        nlevp1   = p_patch(jg)%nlev + 1
        ALLOCATE(topography_smt(nproma,p_patch(jg)%nblks_c))

        ! Compute smooth topography when SLEVE coordinate is used
        IF ( ivctype == 2 .AND. .NOT. lread_smt ) THEN
          CALL compute_smooth_topo(p_patch    = p_patch(jg),                   & !in
            &                      p_int      = p_int_state(jg),               & !in
            &                      topo_c     = ext_data(jg)%atm%topography_c, & !in
            &                      niter      = 25,                            & !in
            &                      topo_smt_c = topography_smt)                  !out
        ENDIF

        ! total shift of model top with respect to global domain
        IF (jg > 1) THEN
          nflatlev(jg) = nflatlev(1) - p_patch(jg)%nshift_total
        ENDIF

        IF (jg > 1 .AND. p_patch(jg)%nshift_total > 0 .AND. nflatlev(jg) <= 1) THEN
          CALL finish(routine, 'flat_height too close to the top of the innermost nested domain')
        ENDIF

        ! Initialize vertical coordinate for cell points
        CALL init_vert_coord(vct_a, vct_b, ext_data(jg)%atm%topography_c, topography_smt, &
          &                  vgrid_buffer(jg)%z_ifc, p_patch(jg)%nlev,                    &
          &                  p_patch(jg)%nblks_c, p_patch(jg)%npromz_c,                   &
          &                  p_patch(jg)%nshift_total, nflatlev(jg) )
        DEALLOCATE(topography_smt)
      END DO

      ! Parallel generation of UUID for distributed 3D field (to
      ! this end, we create a local copy). No global array involved.
      DO jg = 1,n_dom

        nlevels = SIZE(vgrid_buffer(jg)%z_ifc,2)
        ALLOCATE(r_in(nlevels, p_patch(jg)%n_patch_cells), &
          &      glbidx(p_patch(jg)%n_patch_cells))
        iidx = 0
        DO j = 1, p_patch(jg)%n_patch_cells
          jc_c = idx_no(j)
          jb_c = blk_no(j)
          IF (p_patch(jg)%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE
          iidx = iidx + 1
          r_in(:,iidx) = vgrid_buffer(jg)%z_ifc(jc_c,1:nlevels,jb_c)
          glbidx(iidx) = p_patch(jg)%cells%decomp_info%glb_index(j)
        ENDDO
        CALL uuid_generate(p_comm_work, r_in(:,1:iidx), glbidx(1:iidx), &
          &                p_patch(jg)%n_patch_cells_g, vgrid_buffer(jg)%uuid)
        IF (my_process_is_mpi_workroot()) THEN
          CALL uuid_unparse(vgrid_buffer(jg)%uuid, uuid_unparsed)
          WRITE (0,*) "parallel calculation of vgrid UUID. Generated UUID: ", uuid_unparsed
        END IF
        DEALLOCATE(r_in,glbidx)

        ! broadcast the generated UUID st. it is available on all
        ! workers:
        CALL p_bcast(vgrid_buffer(jg)%uuid%data, SIZE(vgrid_buffer(jg)%uuid%data, 1), 0, p_comm_work)

      ENDDO ! jg

      !--- If the user did not provide an external vertical grid
      !--- file, this file will be created:
      IF (create_vgrid) THEN
        DO jg = 1,n_dom
          CALL write_vgrid_file(p_patch(jg), vct_a, vct_b, nflatlev(jg), &
            &                   "vgrid_DOM"//TRIM(int2string(jg, "(i2.2)"))//".nc")
        END DO
      ENDIF  ! create_vgrid

    ELSE

      !--- The user has provided an external vertical grid file; we read its
      !--- contents here:
      DO jg = 1,n_dom
        CALL read_vgrid_file(p_patch(jg), vct_a, vct_b, nflatlev(jg), TRIM(vertical_grid_filename(jg)))

        ! Reset topography_c to z_ifc(nlevp1) in order to ensure identity
        ! in case of nesting (topography blending/feedback)
        nlevp1 = p_patch(jg)%nlev + 1
        ext_data(jg)%atm%topography_c(:,:) = vgrid_buffer(jg)%z_ifc(:,nlevp1,:)

      END DO

    END IF
    !$ACC UPDATE DEVICE(vct_a) ASYNC(1)

  END SUBROUTINE construct_vertical_grid


  !----------------------------------------------------------------------------------------------------
  !> Creates a file for vertical grids.
  !!
  !! Initial implementation by L. Kornblueh (MPI-M), F. Prill (DWD) : 2014-01-29
  !!
  SUBROUTINE write_vgrid_file(p_patch, vct_a, vct_b, nflat, filename)
    TYPE(t_patch),          INTENT(IN)    :: p_patch
    REAL(wp),               INTENT(IN)    :: vct_a(:)             ! param. A of the vertical coordinate
    REAL(wp),               INTENT(IN)    :: vct_b(:)             ! param. B of the vertical coordinate
    INTEGER,                INTENT(INOUT) :: nflat
    CHARACTER(LEN=*),       INTENT(IN)    :: filename
    ! local variables
    CHARACTER(*), PARAMETER    :: routine = modname//"::write_vgrid_file"
    INTEGER                    :: error_status, jk, output_type, nlevp1,       &
      &                           gridtype, iret,                              &
      &                           cdiFileID, cdiVarID_c, cdiVlistID,           &
      &                           cdiInstID, cdiCellGridID, cdiZaxisID,        &
      &                           cdiVarID_vct_a, cdiVarID_vct_b, cdiColumnGridID
    REAL(wp), ALLOCATABLE      :: r1d(:)                                     ! field gathered on I/O processor
    REAL(dp), ALLOCATABLE      :: levels(:)
    CHARACTER(len=132)         :: message_text
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText
    INTEGER :: oneInt(1)

    CALL message(routine, "create vertical grid description file.")

    nlevp1 = p_patch%nlevp1
    IF (my_process_is_mpi_workroot()) THEN
      gridtype    = GRID_UNSTRUCTURED
      output_type = FILETYPE_NC2

      !--- create meta-data
      cdiVlistID = vlistCreate()
      ! define institute
      cdiInstID = institutInq(gribout_config(p_patch%id)%generatingCenter,          &
        &                     gribout_config(p_patch%id)%generatingSubcenter, '', '')
      CALL vlistDefInstitut(cdiVlistID, cdiInstID)

      !--- create grids
      cdiCellGridID   = gridCreate(gridtype, p_patch%n_patch_cells_g)
      CALL gridDefNvertex(cdiCellGridID, p_patch%cells%max_connectivity)
      cdiColumnGridID = gridCreate(gridtype, 1)

      !--- add vertical grid descriptions
      cdiZaxisID = zaxisCreate(ZAXIS_REFERENCE, nlevp1)
      ALLOCATE(levels(nlevp1), STAT=error_status)
      IF (error_status /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
      DO jk = 1, nlevp1
        levels(jk) = REAL(jk,dp)
      END DO
      CALL zaxisDefLevels(cdiZaxisID, levels)
      DEALLOCATE(levels, STAT=error_status)
      IF (error_status /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

      !--- copy UUID for horizontal grid
      !
      ! cells
      CALL gridDefUUID(cdiCellGridID, p_patch%grid_uuid%data)
      CALL gridDefNumber(cdiCellGridID, number_of_grid_used(p_patch%id))
      CALL gridDefPosition(cdiCellGridID, 1)

      CALL zaxisDefNumber(cdiZaxisID, ivctype)

      !--- add variables
      cdiVarID_vct_a    = vlistDefVar(cdiVlistID, cdiColumnGridID, cdiZaxisID, TIME_CONSTANT)
      CALL vlistDefVarName(cdiVlistID, cdiVarID_vct_a, "vct_a")
      CALL vlistDefVarDatatype(cdiVlistID, cdiVarID_vct_a, DATATYPE_FLT64)
      CALL vlistDefVarTsteptype(cdiVlistID, cdiVarID_vct_a, TSTEP_CONSTANT)
      cdiVarID_vct_b    = vlistDefVar(cdiVlistID, cdiColumnGridID, cdiZaxisID, TIME_CONSTANT)
      CALL vlistDefVarName(cdiVlistID, cdiVarID_vct_b, "vct_b")
      CALL vlistDefVarDatatype(cdiVlistID, cdiVarID_vct_b, DATATYPE_FLT64)
      CALL vlistDefVarTsteptype(cdiVlistID, cdiVarID_vct_b, TSTEP_CONSTANT)
      cdiVarID_c        = vlistDefVar(cdiVlistID, cdiCellGridID,   cdiZaxisID, TIME_CONSTANT)
      CALL vlistDefVarName(cdiVlistID, cdiVarID_c, "z_ifc")
      CALL vlistDefVarDatatype(cdiVlistID, cdiVarID_c, DATATYPE_FLT64)
      CALL vlistDefVarTsteptype(cdiVlistID, cdiVarID_c, TSTEP_CONSTANT)
      !--- add "nflat"
      oneInt(1) = nflat
      iret = cdiDefAttInt(cdiVlistID, CDI_GLOBAL, "nflat", DATATYPE_INT32,  1, oneInt)

      !--- open file via CDI
      cdiFileID   = streamOpenWrite(TRIM(filename), output_type)
      ! check if the file could be opened
      IF (cdiFileID < 0) THEN
        CALL cdiGetStringError(cdiFileID, cdiErrorText)
        WRITE(message_text,'(4a)') 'File ', TRIM(filename), &
             ' cannot be opened: ', TRIM(cdiErrorText)
        CALL finish(routine, TRIM(message_text))
      ENDIF

      ! assign the vlist (which must have ben set before)
      CALL streamDefVlist(cdiFileID, cdiVlistID)
      ! streamDefTimestep is required, even without time axis!
      iret = streamDefTimestep(cdiFileID, 0)
      !--- write 1D coordinate arrays:
      CALL streamWriteVar(cdiFileID, cdiVarID_vct_a, vct_a, 0)
      CALL streamWriteVar(cdiFileID, cdiVarID_vct_b, vct_b, 0)
    END IF


    !--- collect CELL coordinate field from compute PEs and write them
    !--- via CDI:

    ! allocate global 2D slice (only on work root PE)
    ALLOCATE(r1d(MERGE(p_patch%n_patch_cells_g, 0, my_process_is_mpi_workroot())), STAT=error_status)
    IF (error_status /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    DO jk=1,nlevp1
      CALL exchange_data(in_array=vgrid_buffer(p_patch%id)%z_ifc(:,jk,:), &
        &                out_array=r1d, gather_pattern=p_patch%comm_pat_gather_c)
      IF (my_process_is_mpi_workroot()) THEN
        CALL streamWriteVarSlice(cdiFileID, cdiVarID_c, jk-1, r1d, 0)

        !--- set UUID for vertical grid
        CALL zaxisDefUUID(cdiZaxisID, vgrid_buffer(p_patch%id)%uuid%data)
      END IF
    END DO
    DEALLOCATE(r1d, STAT=error_status)
    IF (error_status /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    !--- close file
    IF (my_process_is_mpi_workroot()) THEN
      CALL streamClose(cdiFileID)
      IF (cdiCellGridID   /= CDI_UNDEFID) CALL gridDestroy(cdiCellGridID)
      IF (cdiZaxisID      /= CDI_UNDEFID) CALL zaxisDestroy(cdiZaxisID)
      CALL vlistDestroy(cdiVlistID)
    END IF

  END SUBROUTINE write_vgrid_file


  !----------------------------------------------------------------------------------------------------
  !> Opens and reads an existing file with a vertical grid definition.
  !!
  !! Initial implementation by L. Kornblueh (MPI-M), F. Prill (DWD) : 2014-01-29
  !!
  !! Output of this subroutine are:
  !!   vct_a, vct_b, vgrid_buffer%z_ifc
  !!
  !!
  SUBROUTINE read_vgrid_file(p_patch, vct_a, vct_b, nflat, filename)
    TYPE(t_patch),          INTENT(IN)    :: p_patch
    REAL(wp),               INTENT(INOUT) :: vct_a(:)             ! param. A of the vertical coordinate
    REAL(wp),               INTENT(INOUT) :: vct_b(:)             ! param. B of the vertical coordinate
    INTEGER,                INTENT(INOUT) :: nflat
    CHARACTER(LEN=*),       INTENT(IN)    :: filename
    ! local variables
    CHARACTER(*), PARAMETER  :: routine = modname//"::read_vgrid_file"
    LOGICAL                  :: lexists
    INTEGER                  :: cdiFileID, cdiVarID_vct_a, cdiVarID_vct_b,       &
      &                         nlevp1, nmiss, iret, cdiVlistID, cdiGridID,      &
      &                         cdiVarID_z_ifc, cdiZaxisID, oneInt(1)
    TYPE(t_uuid)             :: vfile_uuidOfHGrid            ! uuidOfHGrid contained in the
                                                             ! vertical grid file
    LOGICAL                  :: lmatch                       ! check matching of uuid's

    CHARACTER(len=uuid_string_length) :: uuid_unparsed
    TYPE(t_inputParameters)  :: parameters   !TODO: Move to caller?
    CHARACTER(LEN=MAX_CHAR_LENGTH) :: cdiErrorText

    CALL message(routine, "read vertical grid description file.")

    IF (my_process_is_mpi_workroot()) THEN
      INQUIRE(file=TRIM(filename), exist=lexists)
      IF (.NOT. lexists)  CALL finish(routine, "File "//TRIM(filename)//" not found!")

      !--- open file
      cdiFileID = streamOpenRead(TRIM(filename))
      ! check if the file could be opened
      IF (cdiFileID < 0) THEN
        CALL cdiGetStringError(cdiFileID, cdiErrorText)
        WRITE(message_text,'(4a)') 'File ', TRIM(filename), &
             ' cannot be opened: ', TRIM(cdiErrorText)
        CALL finish(routine, TRIM(message_text))
      ENDIF

      !--- read vct_a, vct_b
      cdiVarID_vct_a     = get_cdi_varID(cdiFileID, "vct_a")
      CALL streamReadVar(cdiFileID, cdiVarID_vct_a, vct_a, nmiss)
      cdiVarID_vct_b     = get_cdi_varID(cdiFileID, "vct_b")
      CALL streamReadVar(cdiFileID, cdiVarID_vct_b, vct_b, nmiss)
      cdiVlistID         = streamInqVlist(cdiFileID)
      iret = cdiInqAttInt(cdiVlistID, CDI_GLOBAL, "nflat", 1, oneInt)
      nflat = oneInt(1)

      !--- get UUID for horizontal grid contained in vertical grid file
      cdiVarID_z_ifc     = get_cdi_varID(cdiFileID, "z_ifc")
      cdiGridID          = vlistInqVarGrid(cdiVlistID, cdiVarID_z_ifc)
      CALL gridInqUUID(cdiGridID, vfile_uuidOfHGrid%data)

      ! compare horizontal UUID contained in the vertical grid file with the
      ! UUID of the horizontal grid file
      lmatch = (p_patch%grid_uuid == vfile_uuidOfHGrid)
      IF (.NOT. lmatch) THEN
        WRITE(message_text,'(a)') 'uuidOfHGrid: Horizontal and vertical grid ', &
          &                       'file do not match!'
        IF (check_uuid_gracefully) THEN
          CALL warning(routine, TRIM(message_text))
        ELSE
          CALL finish(routine, TRIM(message_text))
        END IF
      ENDIF

      !--- get UUID for vertical grid
      cdiZaxisID = vlistInqVarZaxis(cdiVlistID, cdiVarID_z_ifc)
      CALL zaxisInqUUID(cdiZaxisID, vgrid_buffer(p_patch%id)%uuid%data)



      IF (msg_level > 10) THEN
        ! vertical UUID
        CALL uuid_unparse(vgrid_buffer(p_patch%id)%uuid, uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfVGrid from vert. grid file:', TRIM(uuid_unparsed)
        CALL message(TRIM(routine), TRIM(message_text))

        ! vertical grid horizontal UUID
        CALL uuid_unparse(vfile_uuidOfHGrid, uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfHGrid from vert. grid file:', TRIM(uuid_unparsed)
        CALL message(TRIM(routine), TRIM(message_text))

        ! horizontal UUID
        CALL uuid_unparse(p_patch%grid_uuid, uuid_unparsed)
        WRITE(message_text,'(a,a)') 'uuidOfHGrid from hor. grid file:', TRIM(uuid_unparsed)
        CALL message(TRIM(routine), TRIM(message_text))
      ENDIF

    END IF

    ! broadcast data:
    CALL p_bcast(vct_a,                         p_io, p_comm_work)
    CALL p_bcast(vct_b,                         p_io, p_comm_work)
    CALL p_bcast(nflat,                         p_io, p_comm_work)
    CALL p_bcast(vgrid_buffer(p_patch%id)%uuid%data, SIZE(vgrid_buffer(p_patch%id)%uuid%data, 1), p_io, p_comm_work)

    !--- read 3D fields
    nlevp1 = p_patch%nlevp1
    parameters = makeInputParameters(cdiFileID, p_patch%n_patch_cells_g, p_patch%comm_pat_scatter_c)
    CALL read_cdi_3d(parameters, 'z_ifc', nlevp1, vgrid_buffer(p_patch%id)%z_ifc )
    CALL deleteInputParameters(parameters)

    !--- close file
    IF (my_process_is_mpi_workroot()) THEN
      CALL streamClose(cdiFileID)
    END IF

  END SUBROUTINE read_vgrid_file

END MODULE mo_util_vgrid
