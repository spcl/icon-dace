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

! The base CLASS for the PUBLIC INTERFACE for restart writing.

MODULE mo_restart_descriptor
  USE mtime,                        ONLY: datetime, datetimeToString, MAX_DATETIME_STR_LEN
  USE mo_exception,                 ONLY: finish
  USE mo_cf_convention, ONLY: cf_global_info
  USE mo_grid_config,               ONLY: n_dom
  USE mo_impl_constants,            ONLY: SUCCESS
  USE mo_kind,                      ONLY: dp, wp
  USE mo_model_domain,              ONLY: t_patch
  USE mo_mpi, ONLY: my_process_is_work, p_bcast, p_comm_work_2_restart, &
    & p_pe_work, my_process_is_mpi_test, process_mpi_restart_size
  USE mo_restart_nml_and_att,       ONLY: restartAttributeList_make, bcastNamelistStore, &
    & restartAttributeList_write_to_ncdf
  USE mo_key_value_store,           ONLY: t_key_value_store
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_restart_util, ONLY: t_restart_args, create_restart_file_link, &
    & getRestartFilename, restartBcastRoot, t_rfids
  USE mo_restart_var_data,          ONLY: has_valid_time_level
  USE mo_var_list_register_utils,   ONLY: vlr_replicate
#ifndef __NO_ICON_UPATMO__
  USE mo_upatmo_flowevent_utils,    ONLY: t_upatmoRestartAttributes, upatmoRestartAttributesSet
#endif
  USE mo_cdi,                       ONLY: FILETYPE_NC2, FILETYPE_NC4
  USE mo_restart_patch_data, ONLY: t_restartPatchData
#ifndef NOMPI
  USE mo_async_restart_patch_data, ONLY: t_asyncPatchData
#endif
  USE mo_sync_restart_patch_data, ONLY: t_syncPatchData
  USE mo_multifile_restart_patch_data, ONLY: t_multifilePatchData
  USE mo_var_metadata_types, ONLY: t_var_metadata
  USE mo_netcdf_errhandler, ONLY: nf

  USE mo_netcdf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_RestartDescriptor

  ! This IS the actual INTERFACE to the restart writing code (apart from the restart_main_proc PROCEDURE). Its USE IS as follows:
  !
  ! First, AND ONLY once during a run, a t_RestartDescriptor IS created.
  ! This IS done via the factory FUNCTION createRestartDescriptor(), which IS found IN mo_restart.
  !
  ! Then, for each restart that IS to be written, the updatePatch()
  ! method IS used to set the current time dependend information for
  ! each patch.
  ! Once all patches are updated, a single CALL to writeRestart() triggers the actual restart writing.
  ! The updatePatch() - writeRestart() sequence can be repeated ANY number of times.
  !
  ! Finally, destruct() must be called for cleanup. This IS especially important IN the CASE of asynchronous restart writing,
  ! because the destruct() CALL will signal the restart PEs to finish their work, AND wait for them to stop.
  TYPE, ABSTRACT :: t_RestartDescriptor
    CLASS(t_RestartPatchData), POINTER :: patchData(:)
#ifndef NOMPI
    TYPE(t_asyncPatchData), ALLOCATABLE :: aPatchData(:)
#endif
    TYPE(t_syncPatchData), ALLOCATABLE :: sPatchData(:)
    TYPE(t_multifilePatchData), ALLOCATABLE :: mPatchData(:)
    CHARACTER(:), ALLOCATABLE :: modelType
  CONTAINS
    PROCEDURE :: updatePatch => restartDescriptor_updatePatch
    PROCEDURE :: writeFiles => restartDescriptor_writeFiles
    PROCEDURE :: transferGlobalParameters => restartDescriptor_transferGlobalParameters
    PROCEDURE :: defineRestartAttributes => restartDescriptor_defineRestartAttributes
    PROCEDURE(i_construct), DEFERRED :: construct
    PROCEDURE(i_writeRestart), DEFERRED :: writeRestart
    PROCEDURE(i_destruct), DEFERRED :: destruct
  END TYPE t_RestartDescriptor

  ABSTRACT INTERFACE
    SUBROUTINE i_construct(me, modelType)
      IMPORT t_RestartDescriptor
      CLASS(t_RestartDescriptor), INTENT(INOUT), TARGET :: me
      CHARACTER(LEN = *), INTENT(IN) :: modelType
    END SUBROUTINE i_construct

    SUBROUTINE i_writeRestart(me, this_datetime, jstep, opt_output_jfile)
      IMPORT t_RestartDescriptor, datetime
      CLASS(t_RestartDescriptor), INTENT(INOUT), TARGET :: me
      TYPE(datetime), POINTER, INTENT(IN) :: this_datetime
      INTEGER, INTENT(IN) :: jstep
      INTEGER, INTENT(IN), OPTIONAL :: opt_output_jfile(:)
    END SUBROUTINE i_writeRestart

    SUBROUTINE i_destruct(me)
      IMPORT t_RestartDescriptor
      CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
    END SUBROUTINE i_destruct
  END INTERFACE

  CHARACTER(*), PARAMETER :: modname = "mo_restart_descriptor"

CONTAINS

  ! Transfers the modelType, n_dom, AND the namelistArchive to the dedicated restart processes.
  SUBROUTINE restartDescriptor_transferGlobalParameters(me)
    CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
#ifndef NOMPI
    INTEGER :: error, length, bcast_root
    CHARACTER(*), PARAMETER :: routine = modname//":restartDescriptor_transferGlobalParameters"

    bcast_root = restartBcastRoot()
    length = LEN(me%modelType)
    CALL p_bcast(length, bcast_root, p_comm_work_2_restart)
    IF (length /= LEN(me%modelType)) THEN
      DEALLOCATE(me%modelType)
      ALLOCATE(CHARACTER(LEN = length) :: me%modelType, STAT = error)
      IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
    END IF
    CALL p_bcast(me%modelType, bcast_root, p_comm_work_2_restart)
    CALL p_bcast(n_dom, bcast_root, p_comm_work_2_restart)
    CALL vlr_replicate(bcast_root, p_comm_work_2_restart)
    CALL bcastNamelistStore(bcast_root, p_comm_work_2_restart)
#endif
  END SUBROUTINE restartDescriptor_transferGlobalParameters

  ! Update the internal description of the given patch. This should
  ! be called once for every patch before every CALL to
  ! writeRestart().
  SUBROUTINE restartDescriptor_updatePatch(me, patch, opt_pvct, opt_t_elapsed_phy, &
                                        &opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd, &
                                        &opt_nlev_snow, opt_nice_class, opt_ndom,  &
#ifndef __NO_ICON_UPATMO__
                                        &opt_upatmo_restart_atts, &
#endif
                                        &opt_ocean_zlevels, &
                                        &opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces )
    CLASS(t_RestartDescriptor), INTENT(INOUT) :: me
    TYPE(t_patch), INTENT(IN) :: patch
    INTEGER, INTENT(IN), OPTIONAL :: opt_depth_lnd, opt_ndyn_substeps, opt_jstep_adv_marchuk_order, &
                                   & opt_nlev_snow, opt_nice_class, opt_ndom, opt_ocean_zlevels
    REAL(wp), INTENT(IN), OPTIONAL :: opt_pvct(:), opt_t_elapsed_phy(:), opt_ocean_zheight_cellMiddle(:), &
         & opt_ocean_zheight_cellInterfaces(:)
#ifndef __NO_ICON_UPATMO__
    TYPE(t_upatmoRestartAttributes), INTENT(IN), OPTIONAL :: opt_upatmo_restart_atts
#endif
    INTEGER :: jg
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":restartDescriptor_updatePatch"

    IF (my_process_is_mpi_test()) RETURN
    IF(.NOT.my_process_is_work()) CALL finish(routine, "assertion failed")
    DO jg = 1, SIZE(me%patchData)
      IF (patch%id == me%patchData(jg)%description%id) THEN
        CALL me%patchData(jg)%description%update(patch, opt_pvct, opt_t_elapsed_phy,                &
          &      opt_ndyn_substeps, opt_jstep_adv_marchuk_order, opt_depth_lnd,                     &
          &      opt_nlev_snow, opt_nice_class, opt_ndom,                                           &
#ifndef __NO_ICON_UPATMO__
          &      opt_upatmo_restart_atts,                                                           &
#endif
          &      opt_ocean_zlevels, opt_ocean_zheight_cellMiddle, opt_ocean_zheight_cellInterfaces)
      END IF
    END DO
  END SUBROUTINE restartDescriptor_updatePatch

  SUBROUTINE restartDescriptor_defineRestartAttributes(me, rAttribs, rArgs)
    CLASS(t_RestartDescriptor), INTENT(IN) :: me
    TYPE(t_key_value_store), ALLOCATABLE, INTENT(OUT) :: rAttribs
    TYPE(t_restart_args), INTENT(IN) :: rArgs
    INTEGER :: jg, j
    CHARACTER(LEN=15) :: prefix
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: dstring

    ! first the attributes that are independent of the domain
    CALL restartAttributelist_make(rAttribs)
    ! put CF-Convention required restart attributes
    CALL rAttribs%put('bool_int_is_int', 1)
    CALL rAttribs%put('title',       cf_global_info%title)
    CALL rAttribs%put('institution', cf_global_info%institution)
    CALL rAttribs%put('source',      cf_global_info%source)
    CALL rAttribs%put('history',     cf_global_info%history)
    CALL rAttribs%put('references',  cf_global_info%references)
    CALL rAttribs%put('comment',     cf_global_info%comment)
    CALL datetimeToString(rArgs%restart_datetime, dstring)
    CALL rAttribs%put('tc_startdate', dstring)
    CALL rAttribs%put('n_dom', me%patchData(1)%description%opt_ndom)
    CALL rAttribs%put('jstep', rArgs%jstep )
    IF (ALLOCATED(rArgs%output_jfile)) THEN
      DO j = 1, SIZE(rArgs%output_jfile)
        WRITE(prefix, "(a,i2.2)") 'output_jfile_', j
        CALL rAttribs%put(prefix(1:15), rArgs%output_jfile(j) )
      END DO
    END IF
    ! now the stuff that depends on the domain
    DO jg = 1, me%patchData(1)%description%opt_ndom
      CALL put_dom_rstrt_attr(rattribs, me%patchData(jg)%description, jg)
    END DO
  END SUBROUTINE restartDescriptor_defineRestartAttributes

  SUBROUTINE put_dom_rstrt_attr(rattribs, desc, jg)
    TYPE(t_key_value_store), INTENT(INOUT) :: rAttribs
    CLASS(t_restart_patch_description), INTENT(in) :: desc
    INTEGER, INTENT(in) :: jg
    CHARACTER(LEN=2) :: domStr
    CHARACTER(LEN=25) :: prefix
    INTEGER :: j
    IF (desc%id .NE. jg) &
      CALL finish(modname//":defineRestartAttributes", "mismatch of DOM-ID")
    WRITE(domStr, '(i2.2)') desc%id
    CALL rAttribs%put('nold_DOM'//domStr, desc%nold)
    CALL rAttribs%put('nnow_DOM'//domStr, desc%nnow)
    CALL rAttribs%put('nnew_DOM'//domStr, desc%nnew)
    CALL rAttribs%put('nnow_rcf_DOM'//domStr, desc%nnow_rcf)
    CALL rAttribs%put('nnew_rcf_DOM'//domStr, desc%nnew_rcf)
    IF (ALLOCATED(desc%opt_t_elapsed_phy)) THEN
      DO j = 1, SIZE(desc%opt_t_elapsed_phy)
        WRITE (prefix, '(2(a,i2.2))') 't_elapsed_phy_DOM', jg, '_PHY', j
        CALL rAttribs%put(prefix, desc%opt_t_elapsed_phy(j))
      END DO
    END IF
    IF (desc%opt_ndyn_substeps%present) &
      & CALL rAttribs%put('ndyn_substeps_DOM'//domStr, desc%opt_ndyn_substeps%v)
    IF (desc%opt_jstep_adv_marchuk_order%present) &
      & CALL rAttribs%put('jstep_adv_marchuk_order_DOM'//domStr, desc%opt_jstep_adv_marchuk_order%v)
#ifndef __NO_ICON_UPATMO__
    CALL upatmoRestartAttributesSet(desc%id, desc%opt_upatmo_restart_atts, rAttribs)
#endif
  END SUBROUTINE put_dom_rstrt_attr

  SUBROUTINE restartDescriptor_writeFiles(me, rArgs, isSync)
    CLASS(t_restartDescriptor), INTENT(INOUT), TARGET :: me
    TYPE(t_restart_args), INTENT(IN) :: rArgs
    LOGICAL, INTENT(in) :: isSync
    CHARACTER(*), PARAMETER :: routine = modname//":restartDescriptor_writeFiles"
    CHARACTER(:), ALLOCATABLE :: fname
    TYPE(t_rfids) :: rfids
    CLASS(t_RestartPatchData), POINTER :: pData
    TYPE(t_restart_patch_description), POINTER :: desc
    LOGICAL :: lIsWriteProcess
    INTEGER :: jg
    TYPE(t_key_value_store), ALLOCATABLE :: rAttribs

    DO jg = 1, SIZE(me%patchData, 1)
      pData => me%patchData(jg)
      desc => pData%description
      IF (SIZE(pData%varData) .LE. 0) CYCLE
      IF (isSync) THEN
        lIsWriteProcess = p_pe_work .EQ. 0
      ELSE
        lIsWriteProcess = .true.
        IF (.NOT.desc%l_dom_active) CYCLE
        IF (p_pe_work .NE. MOD(desc%id-1, process_mpi_restart_size)) CYCLE
      END IF
      IF (lIsWriteProcess) THEN
        IF (.NOT.ALLOCATED(rAttribs)) &
          & CALL me%defineRestartAttributes(rAttribs, rArgs)
        CALL restartfile_open()
      END IF
      IF (ALLOCATED(me%sPatchData)) THEN
        CALL me%sPatchData(jg)%writeData(rfids%ncid)
#ifndef NOMPI
      ELSE IF (ALLOCATED(me%aPatchData)) THEN
        CALL me%aPatchData(jg)%writeData(rfids%ncid)
#endif
      ELSE
        CALL finish(routine, "multifile does not use this routine!")
      END IF
      IF(lIsWriteProcess) THEN
        CALL create_restart_file_link(fname, TRIM(rArgs%modelType), &
          & desc%id, opt_ndom=desc%opt_ndom)
        CALL nf(nf90_close(rfids%ncid), routine)
      END IF
    END DO
    IF (ALLOCATED(rAttribs)) CALL rAttribs%destruct()
  CONTAINS

    SUBROUTINE restartfile_open()
      CHARACTER(len=MAX_DATETIME_STR_LEN) :: datetimeString
      INTEGER :: i, ncid, tvid
      REAL(dp) :: date_dayas
      TYPE(t_var_metadata), POINTER :: ci
#ifdef DEBUG
      WRITE (nerr,'(a,i6)') routine//' p_pe=',p_pe
#endif
      ! assume all restart variables uses the same file format
      CALL datetimeToString(rArgs%restart_datetime, datetimeString)
      CALL getRestartFilename(desc%base_filename, desc%id, rArgs, fname, date_dayas)
      SELECT CASE(pData%restartType)
      CASE(FILETYPE_NC2)
        WRITE(0,*) "Write netCDF2 restart for: "//TRIM(datetimeString)
        CALL nf(nf90_create(fname, NF90_64BIT_OFFSET, ncid), routine)
      CASE(FILETYPE_NC4)
        WRITE(0,*) "Write netCDF4 restart for: "//TRIM(datetimeString)
        CALL nf(nf90_create(fname, NF90_NETCDF4, ncid), routine)
      CASE default
        CALL finish(routine, "file format for restart variables must be NetCDF")
      END SELECT
      CALL rfids%init(ncid, desc%n_patch_elem_g(1:3), tvid)
      ! set global attributes
      CALL restartAttributeList_write_to_ncdf(rAttribs, ncid)
#ifdef DEBUG
      WRITE (nerr, '(2(a,i6))' )routine//' p_pe=',p_pe,' open netCDF file with ID=',ncid
#endif
      ! go over the all restart variables in the associated array AND define those
      ! that have a valid time level
      DO i = 1, SIZE(pData%varData)
        ci => pData%varData(i)%p%info
        IF(has_valid_time_level(ci, desc%id, desc%nnow, desc%nnow_rcf)) &
          & CALL rfids%def_ncdfvar(ci, desc%hmap(ci%hgrid))
      ENDDO
      CALL nf(nf90_set_fill(ncid, NF90_NOFILL, i), routine)
      CALL nf(nf90_enddef(ncid), routine)
      CALL nf(nf90_put_var(ncid, tvid, REAL(date_dayas), [1]), routine)
    END SUBROUTINE restartfile_open
  END SUBROUTINE restartDescriptor_writeFiles

END MODULE mo_restart_descriptor
