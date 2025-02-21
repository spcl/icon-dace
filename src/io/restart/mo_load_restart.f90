! Module for reading restart files
!
! Note: The asynchronous implementation of the restart output can be
!       found in the module "mo_async_restart"
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

MODULE mo_load_restart
  USE mo_exception,          ONLY: message, finish, warning
  USE mo_load_multifile_restart, ONLY: multifileReadPatch, multifileCheckRestartFiles
  USE mo_load_singlefile_restart, ONLY: singlefileReadPatch, singlefileCheckRestartFiles
  USE mo_model_domain,       ONLY: t_patch
  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_workroot, my_process_is_stdio
  USE mo_multifile_restart_util, ONLY: multifileRestartLinkName
  USE mo_restart_nml_and_att,ONLY: restartAttributeList_read
  USE mo_restart_util,       ONLY: restartSymlinkName
  USE mo_var,                ONLY: t_var_ptr
  USE mo_timer,              ONLY: timer_start, timer_stop, timer_load_restart, timer_load_restart_io, &
    & timer_load_restart_comm_setup, timer_load_restart_communication, &
    & timer_load_restart_get_var_id, timers_level
  USE mo_var_list_register_utils, ONLY: vlr_select_restart_vars, vlr_collect_modelTypes
  USE mo_master_control,     ONLY: get_my_process_name
  USE mo_netcdf_errhandler,  ONLY: nf
  USE mo_netcdf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_restart_files, read_restart_header

  CHARACTER(*), PARAMETER :: modname = "mo_load_restart"

CONTAINS

  !Checks the existence of either a single- OR multifile restart
  !symlink, throwing a warning IF both exist.  Returns the path to the
  !symlink that IS to be used, out_lIsMultifile IS TRUE IF that
  !symlink points to a multifile.
  !
  !The RESULT of the first CALL IS cached, so that this FUNCTION may
  !be called several times without negative consequences.
  SUBROUTINE findRestartFile(modelType, out_lIsMultifile, linkname)
    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: linkname
    CHARACTER(*), INTENT(IN) :: modelType
    LOGICAL, INTENT(OUT) :: out_lIsMultifile

    LOGICAL, SAVE :: cache_isValid = .FALSE., cache_isMultifile
    CHARACTER(:), ALLOCATABLE :: singlefileLinkName, multifileLinkName
    LOGICAL :: haveSinglefileLink, haveMultifileLink
    CHARACTER(*), PARAMETER :: routine = modname//":findRestartFile"

    IF(modelType == '') CALL finish(routine, "assertion failed: invalid modelType")

    IF(.NOT.cache_isValid) THEN
      !get the two possible paths
      CALL restartSymlinkName(modelType, 1, singlefileLinkName, -1)
      CALL multifileRestartLinkName(modelType, multifileLinkName)
      !check whether the respective files exist
      INQUIRE(file = singlefileLinkName, exist = haveSinglefileLink)



      INQUIRE(file = multifileLinkName, exist = haveMultifileLink)

      !determine which path to USE
      IF(haveMultifileLink) THEN
        cache_isMultifile = .TRUE.
        IF (haveSinglefileLink) CALL warning(routine, &
          & "both kinds, single- and multi-file restart-links, present, choosing multi-file")
      ELSE
        IF(haveSinglefileLink) THEN
          cache_isMultifile = .FALSE.
        ELSE
          WRITE (0,*) "singlefileLinkName = ", singlefileLinkName
          WRITE (0,*) "multifileLinkName = ", multifileLinkName
          CALL finish(routine, "no functional symlink found")
        END IF
      END IF
      cache_isValid = .TRUE.
    END IF
    out_lIsMultifile = cache_isMultifile
    IF(cache_isMultifile) THEN
      CALL multifileRestartLinkName(modelType, linkname)
    ELSE
      CALL restartSymlinkName(modelType, 1, linkname, -1)
    END IF
  END SUBROUTINE findRestartFile

  ! Read all namelists used in the previous run
  ! and store them in a buffer. These values will overwrite the
  ! model default, and will later be overwritten if the user has
  ! specified something different for this integration.
  !
  ! Collective across p_comm_work.
  !
  ! This IS used for both single- AND multifile restart files.
  SUBROUTINE readRestartAttributeFile(attributeFile)
    CHARACTER(*), INTENT(IN) :: attributeFile
    INTEGER :: ncid
    CHARACTER(*), PARAMETER :: routine = modname//":readRestartAttributeFile"

    IF(my_process_is_mpi_workroot()) CALL nf(nf_open(attributeFile, NF_NOWRITE, ncid), routine)
    CALL restartAttributeList_read(0, p_comm_work, ncid=ncid)
    IF(my_process_is_mpi_workroot()) CALL nf(nf_close(ncid), routine)
    CALL message(routine, "read namelists and attributes from restart file")
  END SUBROUTINE readRestartAttributeFile

  ! Reads attributes and namelists for all available domains from restart file.
  SUBROUTINE read_restart_header(modelType)
    CHARACTER(*), INTENT(IN) :: modelType
    CHARACTER(:), ALLOCATABLE :: filename, mfaname
    LOGICAL :: lIsMultifile

    ! Must NOT USE a timer here as the timers are NOT initialized yet,
    CALL findRestartFile(modelType, lIsMultifile, filename)
    IF(lIsMultifile) THEN
      mfaname = filename//"/attributes.nc"
      CALL readRestartAttributeFile(mfaname)
      IF(my_process_is_stdio()) CALL multifileCheckRestartFiles(filename)
    ELSE
      CALL readRestartAttributeFile(filename)
      IF(my_process_is_stdio()) CALL singlefileCheckRestartFiles(modelType)
    END IF
  END SUBROUTINE read_restart_header

  SUBROUTINE read_restart_files(p_patch, opt_ndom)
    TYPE(t_patch), INTENT(in) :: p_patch
    INTEGER, OPTIONAL, INTENT(in) :: opt_ndom
    CHARACTER(:), ALLOCATABLE :: restartPath
    LOGICAL :: is_MF, is_timer_init
    INTEGER :: ndom_deopt, imodel, mt_len
    TYPE(t_var_ptr), ALLOCATABLE :: varData(:)
    CHARACTER(LEN=8), ALLOCATABLE :: mTypes(:)

    IF(timers_level >= 5) CALL timer_start(timer_load_restart)

    ! Make sure that all the subcounters are recognized as subcounters on all work processes.
    IF(timers_level >= 7) THEN
      CALL timer_start(timer_load_restart_io)
      CALL timer_stop(timer_load_restart_io)
      CALL timer_start(timer_load_restart_communication)
      CALL timer_stop(timer_load_restart_communication)
    END IF
    is_timer_init = .FALSE.
    IF (my_process_is_mpi_workroot()) &
       WRITE(0,'(a,i0)') "restart: reading restart data for patch ", p_patch%id
    CALL vlr_collect_modelTypes(p_patch%id, mTypes)
    DO imodel = 1, SIZE(mTypes)
      mt_len = LEN_TRIM(mTypes(imodel))
      CALL vlr_select_restart_vars(varData, p_patch%id, mTypes(imodel)(1:mt_len))
      !determine whether we have a multifile to READ
      CALL findRestartFile(mTypes(imodel)(1:mt_len), is_MF, restartPath)
      IF(is_MF) THEN
        !multifile loading also uses these two timers
        IF(timers_level >= 7 .AND. .NOT.is_timer_init) THEN
          CALL timer_start(timer_load_restart_comm_setup)
          CALL timer_stop(timer_load_restart_comm_setup)
          CALL timer_start(timer_load_restart_get_var_id)
          CALL timer_stop(timer_load_restart_get_var_id)
          is_timer_init = .TRUE.
        END IF
        !load the multifile
        CALL multifileReadPatch(varData, p_patch, restartPath)
      ELSE
        !can't USE restartPath here because the path IS patch dependent IN the single file CASE
        ndom_deopt = -1
        IF (PRESENT(opt_ndom)) ndom_deopt = opt_ndom
        CALL singlefileReadPatch(varData, mTypes(imodel)(1:mt_len), p_patch, ndom_deopt)
      END IF
      DEALLOCATE(varData)
    END DO
    IF(timers_level >= 5) CALL timer_stop(timer_load_restart)
    CALL message('','')
    CALL message('','----------------------------------------------')
    CALL message('','')
  END SUBROUTINE read_restart_files

END MODULE mo_load_restart
