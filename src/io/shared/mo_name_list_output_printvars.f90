! Auxiliary module: Print a list of all output variables.
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

MODULE mo_name_list_output_printvars

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_char
#ifdef GRIBAPI
  USE mo_cdi,                               ONLY: streamOpenWrite, FILETYPE_GRB2, gridCreate,                &
    &                                             GRID_UNSTRUCTURED, TAXIS_ABSOLUTE,                         &
    &                                             vlistDestroy, streamClose, streamDefVlist,                 &
    &                                             streamWriteVarSlice, gridDestroy,            &
    &                                             streamOpenRead, streamInqVlist, vlistinqvarname,           &
    &                                             taxisDestroy, vlistCreate, taxisCreate, institutInq,       &
    &                                             vlistDefTaxis, vlistDefInstitut
  USE mo_mpi,                               ONLY: get_my_global_mpi_id
  USE mo_util_file,                         ONLY: util_unlink
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxisList, t_verticalAxis
  USE mo_level_selection_types,             ONLY: t_level_selection
  USE mo_name_list_output_zaxes,            ONLY: setup_ml_axes_atmo, setup_zaxes_oce
  USE mo_master_control,                    ONLY: my_process_is_atmo, my_process_is_ocean, my_process_is_jsbach, &
    &                                             my_process_is_waves
#ifndef __NO_ICON_WAVES__
  USE mo_waves_vertical_axes,               ONLY: setup_zaxes_waves
#endif
#ifndef __NO_JSBACH__
  USE mo_atm_phy_nwp_config,                ONLY: atm_phy_nwp_config
  USE mo_grid_config,                       ONLY: n_dom
  USE mo_impl_constants,                    ONLY: LSS_JSBACH
  USE mo_jsb_vertical_axes,                 ONLY: setup_zaxes_jsbach
#endif
  USE mo_util_cdi,                          ONLY: create_cdi_variable
#endif
  USE mo_cdi,                               ONLY: cdi_max_name
  USE mo_gribout_config,                    ONLY: t_gribout_config
  USE mo_kind,                              ONLY: wp, dp
  USE mo_impl_constants,                    ONLY: SUCCESS, vname_len
  USE mo_cf_convention,                     ONLY: t_cf_var
  USE mo_exception,                         ONLY: finish, message_text
  USE mo_var_metadata_types,                ONLY: t_var_metadata
  USE mo_var_list_register,                 ONLY: t_vl_register_iter
  USE mo_var,                               ONLY: t_var, level_type_ml
  USE mo_dictionary,                        ONLY: t_dictionary
  USE mo_util_sort,                         ONLY: quicksort
  USE mo_util_string,                       ONLY: remove_duplicates, toupper, tolower, int2string

  IMPLICIT NONE
  PRIVATE

  ! subroutines
  PUBLIC :: print_var_list

  !------------------------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_printvars'


CONTAINS

  !------------------------------------------------------------------------------------------------
  !> @return variable name without time level or tile suffix
  !
  CHARACTER(LEN=vname_len) FUNCTION get_var_basename(info)
    TYPE(t_var_metadata) :: info
    INTEGER :: endidx
    CHARACTER(LEN=1) :: suffix_str

    ! first cut off the time level suffix:
    endidx = INDEX(info%name,'.TL')
    IF (endidx .EQ. 0) THEN
      endidx = LEN_TRIM(info%name)
    ELSE
      endidx = endidx - 1
    END IF
    suffix_str = " "
    ! condense three-digit suffices "001, 002, 003, ..." into "*"
    ! for container variables
    IF (endidx > 3) THEN
      IF (info%lcontained .AND. &
          .NOT. is_number(info%name(endidx-3:endidx-3)) .AND. &
        &       is_number(info%name(endidx-2:endidx-2)) .AND. &
        &       is_number(info%name(endidx-1:endidx-1)) .AND. &
        &       is_number(info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 3
        suffix_str = "*"
      END IF
    END IF
    ! condense two-digit suffices "_01, _02, _03, ..." into "*"
    IF (endidx > 2) THEN
      IF ((info%name(endidx-2:endidx-2) == "_")    .AND. &
        &  is_number(info%name(endidx-1:endidx-1)) .AND. &
        &  is_number(info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 2
        suffix_str = "*"
      END IF
    END IF
    ! condense one-digit suffices "_1, _2, _3, ..." into "_*"
    IF ((endidx > 1) .AND. (suffix_str == " ")) THEN
      IF ((info%name(endidx-1:endidx-1) == "_") .AND. &
        &  is_number(info%name(endidx  :endidx  ))) THEN
        endidx = endidx - 1
        suffix_str = "*"
      END IF
    END IF
    get_var_basename = info%name(1:endidx)//suffix_str

  CONTAINS
    LOGICAL FUNCTION is_number(chr)
      CHARACTER, INTENT(IN) :: chr
      is_number = (IACHAR(chr) - IACHAR('0')) <= 9
    END FUNCTION is_number
  END FUNCTION get_var_basename

#ifdef GRIBAPI
  !------------------------------------------------------------------------------------------------
  !> @return GRIB2 short name of a given variable.
  !
  !  In order to make sure that we really get the true GRIB2 shortName
  !  (which depends on the definition files used), we choose a rather
  !  awkward approach: A temporary file is created and read in
  !  immediately afterwards. This is slow!
  !
  SUBROUTINE identify_grb2_shortname(info, vname, verticalAxisList, gribout_config, i_lctype, &
    &                                out_varnames_dict)
    TYPE (t_var_metadata),                      INTENT(IN)    :: info
    CHARACTER(kind=c_char, LEN=cdi_max_name+1), INTENT(OUT)   :: vname
    TYPE(t_verticalAxisList),                   INTENT(INOUT) :: verticalAxisList
    TYPE(t_gribout_config),                     INTENT(IN)    :: gribout_config
    INTEGER,                                    INTENT(IN)    :: i_lctype
    TYPE(t_dictionary),                         INTENT(IN)    :: out_varnames_dict
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::identify_grb2_shortname"
    CHARACTER(*), PARAMETER :: tmp_filename_base = "tmpfile.grb"
    TYPE(t_verticalAxis), POINTER :: zaxis
    INTEGER :: tmp_vlistID, tmp_gridID, tmp_zaxisID, tmp_varID, tmp_streamID, &
      & nmiss, tmp_taxisID, tmp_cdiInstID,  ierrstat
    REAL(wp) :: tmp_var1(1)
    CHARACTER(LEN=LEN(tmp_filename_base)+32) :: tmp_filename

    ! retrieve vertical axis object
    zaxis => verticalAxisList%getEntry(icon_zaxis_type=info%vgrid)
    IF (.NOT. ASSOCIATED(zaxis)) THEN
      WRITE (message_text,*) "variable '", TRIM(info%name), "' :Zaxis no. ", info%vgrid, " undefined."
      CALL finish(routine, message_text)
    END IF
    tmp_zaxisID = zaxis%cdi_id
    ! create dummy variable list
    tmp_vlistID = vlistCreate()
    tmp_cdiInstID = institutInq(gribout_config%generatingCenter,          &
      &                         gribout_config%generatingSubcenter, '', '')
    ! define Institute
    CALL vlistDefInstitut(tmp_vlistID, tmp_cdiInstID)
    ! create a dummy time axis 
    tmp_taxisID = taxisCreate(TAXIS_ABSOLUTE) 
    ! assign the time axis to the variable list 
    CALL vlistDefTaxis(tmp_vlistID, tmp_taxisID) 
    ! create dummy grid (size=1):
    tmp_gridID  = gridCreate(GRID_UNSTRUCTURED, 1)
    ! add the variable in question to the vlist:
    tmp_varID = create_cdi_variable(tmp_vlistID, tmp_gridID, tmp_zaxisID,        &
      &                             info, 0._dp, FILETYPE_GRB2,                  &
      &                             gribout_config, i_lctype, out_varnames_dict)
    ! open temporary GRIB2 file, write variable:
    tmp_filename = tmp_filename_base//"."//int2string(get_my_global_mpi_id(),'(i0)')
    tmp_streamID = streamOpenWrite(TRIM(tmp_filename), FILETYPE_GRB2) 
    CALL streamDefVlist(tmp_streamID, tmp_vlistID) 
    tmp_var1(:) = 0._wp
    nmiss       = 0
    CALL streamWriteVarSlice(tmp_streamID, tmp_varID, 0, tmp_var1, nmiss) 
    CALL streamClose(tmp_streamID) 
    CALL vlistDestroy(tmp_vlistID)
    CALL taxisDestroy(tmp_taxisID)
    CALL gridDestroy(tmp_gridID) 
    ! Re-open file:
    tmp_streamID = streamOpenRead(TRIM(tmp_filename)) 
    ! Get the variable list of the dataset 
    tmp_vlistID = streamInqVlist(tmp_streamID) 
    CALL vlistInqVarName(tmp_vlistID, tmp_varID, vname)   
    ! Close the input stream 
    CALL streamClose(tmp_streamID) 
    ierrstat = util_unlink(TRIM(tmp_filename))
  END SUBROUTINE identify_grb2_shortname
#endif

  !------------------------------------------------------------------------------------------------
  !> Print list of all output variables (LaTeX table formatting).
  !
  SUBROUTINE print_var_list(out_varnames_dict, print_patch_id, gribout_config, &
    &                       i_lctype)
    TYPE(t_dictionary),     INTENT(IN) :: out_varnames_dict
    TYPE(t_gribout_config), INTENT(IN) :: gribout_config
    INTEGER,                INTENT(IN) :: i_lctype, print_patch_id

    CHARACTER(*), PARAMETER :: routine = modname//"::print_var_list"
    INTEGER,      PARAMETER :: max_str_len = cdi_max_name + 1 + 128 + vname_len + 99
    CHARACTER(*), PARAMETER :: varprefix = "\varname{", CR = " \\[0.5em]"
    CHARACTER(*), PARAMETER :: descrprefix = "\vardescr{"
    INTEGER,      PARAMETER :: PREF = LEN_TRIM(varprefix) + 1
#ifdef GRIBAPI
    TYPE(t_level_selection),  POINTER              :: tmp_level_selection => NULL()
    TYPE(t_verticalAxisList), TARGET               :: tmp_verticalAxisList
    TYPE(t_verticalAxisList), POINTER              :: it
#endif
    TYPE (t_var_metadata),    POINTER              :: info
    TYPE(t_var),     POINTER              :: elem
    TYPE(t_cf_var),           POINTER              :: this_cf
    INTEGER                                        :: i, iv, nout_vars, iout_var, ierrstat
    CHARACTER(kind=c_char, LEN = cdi_max_name + 1) :: vname
    CHARACTER(LEN=max_str_len), ALLOCATABLE        :: out_vars(:)
    CHARACTER(len=128)                             :: descr_string
    TYPE(t_vl_register_iter) :: vl_iter
    ! ---------------------------------------------------------------------------

#ifdef GRIBAPI
    IF (my_process_is_ocean()) THEN
      CALL setup_zaxes_oce(tmp_verticalAxisList, tmp_level_selection)

    ELSE IF (my_process_is_atmo()) THEN
      CALL setup_ml_axes_atmo(tmp_verticalAxisList, tmp_level_selection, print_patch_id)
#ifndef __NO_JSBACH__
      IF (ANY(atm_phy_nwp_config(1:n_dom)%inwp_surface == LSS_JSBACH)) &
          & CALL setup_zaxes_jsbach(tmp_verticalAxisList)
#endif
    ELSE IF (my_process_is_jsbach()) THEN
#ifndef __NO_JSBACH__
      CALL setup_zaxes_jsbach(tmp_verticalAxisList)
#endif
    ELSE IF (my_process_is_waves()) THEN
#ifndef __NO_ICON_WAVES__
      CALL setup_zaxes_waves(tmp_verticalAxisList, tmp_level_selection, print_patch_id)
#endif
    ENDIF
    it => tmp_verticalAxisList
    DO
      IF (.NOT. ASSOCIATED(it%axis)) CALL finish(routine, "Internal error!")
      CALL it%axis%cdiZaxisCreate()
      IF (.NOT. ASSOCIATED(it%next)) EXIT
      it => it%next
    END DO
#endif

    ! count the no. of output variables:
    nout_vars = 0
    DO WHILE(vl_iter%next())
      IF (vl_iter%cur%p%patch_id /= print_patch_id) CYCLE
      IF(.NOT.vl_iter%cur%p%loutput) CYCLE
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE
      DO iv = 1, vl_iter%cur%p%nvars
        IF (vl_iter%cur%p%vl(iv)%P%info%loutput) nout_vars = nout_vars + 1
      END DO
    END DO
    ! allocate sufficient space
    ALLOCATE(out_vars(nout_vars), STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
    iout_var = 0
    DO WHILE(vl_iter%next())
      IF (vl_iter%cur%p%patch_id /= print_patch_id) CYCLE
      ! Inspect only model level variables
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE
      ! Do not inspect element if output is disabled
      IF(.NOT.vl_iter%cur%p%loutput) CYCLE
      DO iv = 1, vl_iter%cur%p%nvars
        elem => vl_iter%cur%p%vl(iv)%p
        info => elem%info
        ! Do not inspect element if output is disabled
        IF (.NOT.info%loutput) CYCLE
        this_cf => info%cf
        IF (info%post_op%lnew_cf) this_cf => info%post_op%new_cf
        ! if no short is available and if the variable is a
        ! "reference" into another variable, then search for this
        ! source variable:
        IF ((LEN_TRIM(this_cf%long_name) == 0) .AND. info%lcontained) THEN
          this_cf => info%cf
          IF (ASSOCIATED(elem%ref_to)) THEN
            IF (elem%ref_to%info%ncontained > 0) this_cf => elem%ref_to%info%cf
          END IF
        END IF
#ifdef GRIBAPI
        CALL identify_grb2_shortname(info, vname, tmp_verticalAxisList,     &
          &                          gribout_config, i_lctype,       &
          &                          out_varnames_dict)
#else
        vname = info%name
#endif
        vname = tolower(vname)
        IF ((vname(1:3) == "var") .OR. (vname(1:5) == "param")) THEN
          vname = ""
        ELSE
          vname = varprefix//TRIM(vname)//"}"
        END IF

        descr_string = this_cf%long_name
        ! upcase first letter of description string:
        descr_string(1:1) = toupper(descr_string(1:1))

        iout_var = iout_var + 1
        WRITE (out_vars(iout_var),'(8a)') varprefix, &
          & tolower(get_var_basename(info)), "} & ", &
          & TRIM(vname), ' & ', descrprefix, TRIM(descr_string), "}"
      ENDDO
    ENDDO
#ifdef GRIBAPI
    CALL tmp_verticalAxisList%finalize()
#endif
    ! sort and remove duplicates
    CALL quicksort(out_vars(1:nout_vars))
    CALL remove_duplicates(out_vars(1:nout_vars), nout_vars)

    ! print table, but add a gap when new alphabetical letter starts:
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) "List of output variables"
    WRITE (0,*) "----------------------------------------------------------------------"
    WRITE (0,*) " "
    DO i = 1, nout_vars-1
      ! check for the initial character of the variable name:
      WRITE(0,*) TRIM(out_vars(i)) // CR(1:MERGE(3,LEN(CR), &
        &  out_vars(i)(PREF:PREF) == out_vars(i+1)(PREF:PREF)))
    END DO
    WRITE (0,*) TRIM(out_vars(nout_vars)) // CR(1:3)
    WRITE (0,*) " "

    DEALLOCATE(out_vars, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE print_var_list

END MODULE mo_name_list_output_printvars
