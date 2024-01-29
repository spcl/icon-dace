!
! This module provides an input instruction list that is used to
! determine, which variables may be read from which input file.
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

MODULE mo_input_instructions

    USE mo_exception,          ONLY: message, finish, message_text
    USE mo_impl_constants,     ONLY: SUCCESS, MODE_DWDANA, MODE_ICONVREMAP, MODE_IAU, MODE_IAU_OLD, &
      &                              MODE_COMBINED, MODE_COSMO, max_ntracer, vname_len
    USE mo_initicon_config,    ONLY: initicon_config, lread_ana, ltile_coldstart, lp2cintp_incr,    &
      &                              lp2cintp_sfcana, lvert_remap_fg
    USE mo_initicon_types,     ONLY: ana_varnames_dict
    USE mo_input_request_list, ONLY: t_InputRequestList
    USE mo_lnd_nwp_config,     ONLY: lsnowtile
    USE mo_model_domain,       ONLY: t_patch
    USE mo_util_string,        ONLY: difference, add_to_list, one_of
    USE mo_util_table,         ONLY: t_table, initialize_table, add_table_column, set_table_entry,  &
      &                              print_table, finalize_table
    USE mo_var_list_register_utils, ONLY: vlr_group

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: t_readInstructionList, readInstructionList_make, t_readInstructionListPtr
    PUBLIC :: kInputSourceNone, kInputSourceFg, kInputSourceAna, kInputSourceBoth, kInputSourceCold
    PUBLIC :: kStateNoFetch, kStateFailedFetch, kStateRead,  kInputSourceAnaI, kInputSourceFgAnaI

    ! The possible RETURN values of readInstructionList_sourceOfVar().
    ! kInputSourceBoth  : First guess and analysis increment read from file
    ! kInputSourceAnaI  : Analysis interpolated from parent grid
    ! kInputSourceFgAnaI: First guess read from file, analysis increment interpolated from parent grid
    ENUM, BIND(C)
        ENUMERATOR :: kInputSourceUnset = 1, kInputSourceNone, kInputSourceFg, kInputSourceAna, &
          &           kInputSourceBoth, kInputSourceCold, kInputSourceAnaI, kInputSourceFgAnaI
    END ENUM

    ! The possible values for statusFg AND statusAna:
    ENUM, BIND(C)
        ENUMERATOR :: kStateNoFetch = 1, kStateFailedFetch, kStateRead, kStateFailedOptFetch
    END ENUM

    ! The readInstructionList is used to tell the fetch_dwd*()
    ! routines which fields may be read from first guess and/or
    ! analysis, and to signal whether we have found data in the first
    ! guess file, so that we can fail correctly.  I. e. for each
    ! variable, we expect a pair of calls to
    !
    !     wantVarXXX()
    !     handleErrorXXX()
    !
    ! The wantVarXXX() CALL tells the caller whether they should try
    ! to READ the variable, the handleErrorXXX() CALL tells the
    ! ReadInstructionList whether that READ was successfull, possibly
    ! panicking with a `finish()` CALL IN CASE there are no options
    ! left to READ the DATA.  If a variable is declared as optional,
    ! handleErrorXXX() won't panick if no data is available.  Instead
    ! sourceOfVar will be set to kInputSourceCold, indicating that the
    ! variable should experience some sort of coldstart
    ! initialization.
    ! If there exists a fallback variable for a particular variable, 
    ! optionalReadResultXXX() has to be used instead of handleErrorXXX(). 
    ! The former won't panick if no data is available.
    TYPE :: t_readInstructionList
        TYPE(t_readInstruction), POINTER :: list(:)
        INTEGER :: nInstructions

    CONTAINS
      ! tell an InputRequestList what variables are required:
        PROCEDURE :: fileRequests => readInstructionList_fileRequests

        ! inquire whether an attempt should be made to read a variable
        PROCEDURE :: wantVar => readInstructionList_wantVar
        PROCEDURE :: wantVarFg => readInstructionList_wantVarFg
        PROCEDURE :: wantVarAna => readInstructionList_wantVarAna

        ! inform the ReadInstructionList about the result of a read,
        ! possibly triggering a `finish()` call if there is no
        ! alternative left to read the variable
        PROCEDURE :: handleError => readInstructionList_handleError
        PROCEDURE :: handleErrorFg => readInstructionList_handleErrorFg
        PROCEDURE :: handleErrorAna => readInstructionList_handleErrorAna

        ! inform the ReadInstructionList about the result of a read (nofail variant)
        PROCEDURE :: optionalReadResult => readInstructionList_optionalReadResult
        PROCEDURE :: optionalReadResultFg => readInstructionList_optionalReadResultFg
        PROCEDURE :: optionalReadResultAna => readInstructionList_optionalReadResultAna

        ! returns kInputSourceNone, kInputSourceFg, kInputSourceAna,
        ! OR kInputSourceBoth:
        PROCEDURE :: sourceOfVar => readInstructionList_sourceOfVar

        ! overrides the automatic calculation of the input source,
        ! argument must be 'kInputSource\(None\|Fg\|Ana\|Both\)':
        PROCEDURE :: setSource => readInstructionList_setSource

        ! returns the result of a read attempt for a particular field:
        PROCEDURE :: fetchStatus => readInstructionList_fetchStatus
        
        ! print a table with the information obtained via
        ! `handleErrorXXX()`, and `setSource()`; THIS DEPENDS ON THE
        ! ACTUAL READ ATTEMPTS AND THEIR RESULTS, NOT ON `lReadFg` or
        ! `lReadAna`.
        PROCEDURE :: printSummary => readInstructionList_printSummary

        PROCEDURE :: destruct => readInstructionList_destruct

        PROCEDURE, PRIVATE :: construct => readInstructionList_construct
        PROCEDURE, PRIVATE :: resize => readInstructionList_resize
        PROCEDURE, PRIVATE :: findInstruction => readInstructionList_findInstruction
    END TYPE t_readInstructionList

    TYPE :: t_readInstructionListPtr
        TYPE(t_readInstructionList), POINTER :: ptr
    END TYPE t_readInstructionListPtr

    TYPE :: t_readInstruction
        CHARACTER(LEN = vname_len) :: varName
        ! These reflect the result of interpreting the variable
        ! groups, they are not used in the table output:
        LOGICAL :: lReadFg, lReadAna

        ! TRUE: variable is read if it is present in the file, but it
        ! is not necessary to start the run. (re-set to
        ! .FALSE. according to the `fg_checklist` namelist parameter):
        LOGICAL :: lOptionalFg
                              
        ! Panick if reading from analysis fails (set according to the
        ! `ana_checklist` namelist parameter):
        LOGICAL :: lRequireAna

        ! These reflect what read attempts really have been made, and
        ! what their results were. This is the basis for the table
        ! output.
        INTEGER :: statusFg, statusAna

        ! If this is not kInputSourceUnset, it overrides the automatic
        ! calculation of the input source.
        INTEGER :: sourceOverride
    CONTAINS
        PROCEDURE :: source => readInstruction_source
    END TYPE t_readInstruction

    CHARACTER(LEN = *), PARAMETER :: modname = 'mo_input_instructions'

CONTAINS

    SUBROUTINE collectGroupFg(outGroup, outGroupSize, init_mode)
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize
        INTEGER, INTENT(IN) :: init_mode

        ! local tracerGroup to add all vars in group tracer_fg_in to first guess
        CHARACTER(LEN = vname_len) :: tracerGroup(max_ntracer)
        INTEGER :: tracerGroupSize

        CHARACTER(len=vname_len) :: jsbGroup(1000)
        INTEGER :: jsbGroupSize

        SELECT CASE(init_mode)
            CASE(MODE_DWDANA, MODE_ICONVREMAP)
                CALL vlr_group('mode_dwd_fg_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
                CALL vlr_group('tracer_fg_in', tracerGroup, tracerGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
                CALL vlr_group('jsb_init_vars', jsbGroup, jsbGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
                ! fgGroup += tracerGroup
                CALL add_to_list(outGroup, outGroupSize, tracerGroup, tracerGroupSize)
                CALL add_to_list(outGroup, outGroupSize, jsbGroup, jsbGroupSize)
            CASE(MODE_IAU)
                CALL vlr_group('mode_iau_fg_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
                CALL vlr_group('tracer_fg_in', tracerGroup, tracerGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
                CALL vlr_group('jsb_init_vars', jsbGroup, jsbGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
                ! fgGroup += tracerGroup
                CALL add_to_list(outGroup, outGroupSize, tracerGroup, tracerGroupSize)
                CALL add_to_list(outGroup, outGroupSize, jsbGroup, jsbGroupSize)
            CASE(MODE_IAU_OLD)
                CALL vlr_group('mode_iau_old_fg_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_COMBINED)
                CALL vlr_group('mode_combined_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_COSMO)
                CALL vlr_group('mode_cosmo_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE DEFAULT
                outGroupSize = 0
        END SELECT
    END SUBROUTINE collectGroupFg

    ! Sub-list of optional first guess fields
    ! Fields in this list are read from the first guess field if they are present, 
    ! but they are not needed to start the model.
    !
    ! ToDo: 
    ! So far, this list has to be created manually. In the future this 
    ! should be done automatically via (add_var) metadata flags. 
    !
    ! The JSBACH fields are marked optional to allow for initialization from files
    ! without JSBACH fields. In that case JSBACH uses its own initialization files
    ! -- either from climatologies or from IFS, depending on JSBACH configuration.
    ! (Actually, the optional fields overwrite JSBACH's initialization.) If JSBACH
    ! fields are expected, one of them can be marked mandatory in the initicon_nml:
    ! check_fg list.
    SUBROUTINE collectGroupFgOpt(outGroup, outGroupSize)
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize

        CHARACTER(len=vname_len) :: jsbGroup(1000)
        INTEGER :: jsbGroupSize

        outGroup(1:31) = (/'alb_si       ','rho_snow_mult','aer_ss       ','aer_or       ', &
          &                'aer_bc       ','aer_su       ','aer_du       ','plantevap    ', &
          &                't_sk         ','t2m_bias     ','hsnow_max    ','snow_age     ', &
          &                'qg           ','qh           ','qnc          ','qni          ', &
          &                'qnr          ','qns          ','qng          ','qnh          ', &
          &                'rh_avginc    ','t_avginc     ','t_wgt_avginc ','p_avginc     ', &
          &                'clmf_a       ','clmf_p       ','clmf_d       ','clnum_a      ', &
          &                'clnum_p      ','clnum_d      ','vabs_avginc  '/)
        outGroupSize  = 31

        CALL vlr_group('jsb_init_vars', jsbGroup, jsbGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
        CALL add_to_list(outGroup, outGroupSize, jsbGroup, jsbGroupSize)

    END SUBROUTINE collectGroupFgOpt

    SUBROUTINE collectGroupAna(outGroup, outGroupSize, init_mode)
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize
        INTEGER, INTENT(IN) :: init_mode

        SELECT CASE(init_mode)
            CASE(MODE_DWDANA, MODE_ICONVREMAP)
                CALL vlr_group('mode_dwd_ana_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_IAU)
                CALL vlr_group('mode_iau_ana_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE(MODE_IAU_OLD)
                CALL vlr_group('mode_iau_old_ana_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE DEFAULT
                outGroupSize = 0
        END SELECT
    END SUBROUTINE collectGroupAna

    SUBROUTINE collectGroupAnaAtm(outGroup, outGroupSize, init_mode)
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize
        INTEGER, INTENT(IN) :: init_mode

        SELECT CASE(init_mode)
            CASE(MODE_IAU, MODE_IAU_OLD)
                CALL vlr_group('mode_iau_anaatm_in', outGroup, outGroupSize, loutputvars_only=.FALSE., lremap_lonlat=.FALSE.)
            CASE DEFAULT
                outGroupSize = 0
        END SELECT
    END SUBROUTINE collectGroupAnaAtm

    SUBROUTINE copyGroup(inGroup, inGroupSize, outGroup, outGroupSize)
        CHARACTER(LEN = vname_len), INTENT(IN) :: inGroup(:)
        INTEGER, INTENT(IN) :: inGroupSize
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: outGroup(:)
        INTEGER, INTENT(OUT) :: outGroupSize

        outGroup(1:inGroupSize) = inGroup(1:inGroupSize)
        outGroupSize = inGroupSize
    END SUBROUTINE copyGroup

    SUBROUTINE mergeAnaIntoFg(anaGroup, anaGroupSize, fgGroup, fgGroupSize, init_mode)
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: anaGroup(:), fgGroup(:)
        INTEGER, INTENT(INOUT) :: anaGroupSize, fgGroupSize
        INTEGER, INTENT(IN)    :: init_mode

        ! fgGroup += anaGroup
        CALL add_to_list(fgGroup, fgGroupSize, anaGroup, anaGroupSize)

        ! Remove fields 'u', 'v', 'temp', 'pres' except in VREMAP mode, where the diagnostic variable set
        ! is allowed as an alternative to the corresponding prognostic variable set (vn, theta_v, rho)
        IF (init_mode /= MODE_ICONVREMAP) CALL difference(fgGroup, fgGroupSize, (/'u   ','v   ','temp','pres'/), 4)

        ! anaGroup = --
        anaGroupSize = 0
    END SUBROUTINE mergeAnaIntoFg

    SUBROUTINE collectGroups(p_patch, init_mode, fgGroup, fgGroupSize, &
      &                      fgOptGroup, fgOptGroupSize, anaGroup, anaGroupSize)
        TYPE(t_patch), INTENT(IN) :: p_patch
        INTEGER, INTENT(in) :: init_mode
        CHARACTER(LEN = vname_len), INTENT(INOUT) :: anaGroup(:), fgGroup(:), fgOptGroup(:)
        INTEGER, INTENT(OUT) :: anaGroupSize, fgGroupSize, fgOptGroupSize

        CHARACTER(LEN = *), PARAMETER :: routine = modname//':collectGroups'
        CHARACTER(LEN=vname_len), DIMENSION(200) :: anaAtmGroup
        INTEGER :: anaAtmGroupSize, jg
        LOGICAL :: lRemoveSnowfrac

        ! get the raw DATA
        CALL collectGroupFg(fgGroup, fgGroupSize, init_mode)
        CALL collectGroupFgOpt(fgOptGroup, fgOptGroupSize)
        CALL collectGroupAna(anaGroup, anaGroupSize, init_mode)
        CALL collectGroupAnaAtm(anaAtmGroup, anaAtmGroupSize, init_mode)
        jg = p_patch%id


        ! integrate the information of these three groups, the init_mode, AND some flags to produce the effective fgGroup AND anaGroup
        SELECT CASE(init_mode)
            CASE(MODE_DWDANA, MODE_ICONVREMAP)
                IF(.NOT.lread_ana) THEN
                    ! lump together fgGroup and anaGroup
                    CALL mergeAnaIntoFg(anaGroup, anaGroupSize, fgGroup, fgGroupSize, init_mode)
                ENDIF
                IF (init_mode == MODE_ICONVREMAP) CALL add_to_list(fgGroup, fgGroupSize, 'smi')

            CASE(MODE_IAU, MODE_IAU_OLD)
                ! in case of tile coldstart, we can omit snowfrac_lc
                ! Remove field 'snowfrac_lc' from FG list
                lRemoveSnowfrac = ltile_coldstart
                IF(init_mode == MODE_IAU .AND. .NOT. lsnowtile) lRemoveSnowfrac = .TRUE.
                IF(lRemoveSnowfrac) CALL difference(fgGroup, fgGroupSize, (/'snowfrac_lc'/), 1)

                IF (.NOT. lp2cintp_incr(jg) .AND. .NOT. lp2cintp_sfcana(jg) ) THEN
                    ! full ANA read

                ELSE IF (lp2cintp_incr(jg) .AND. .NOT. lp2cintp_sfcana(jg)) THEN
                    ! SFC-ANA read
                    ! atmospheric analysis fieds are interpolated from parent domain, 
                    ! however surface analysis fields are read from file

                    ! Remove fields atmospheric analysis fields from anaGroup
                    CALL difference(anaGroup, anaGroupSize, anaAtmGroup, anaAtmGroupSize)

                ELSE IF (lp2cintp_incr(jg) .AND. lp2cintp_sfcana(jg) ) THEN
                    ! no ANA-read
                    ! lump together fgGroup and anaGroup
                    CALL mergeAnaIntoFg(anaGroup, anaGroupSize, fgGroup, fgGroupSize, init_mode)

                ELSE
                    WRITE(message_text,'(a,l1,a,l1,a)') 'Combination lp2cintp_incr=',lp2cintp_incr(jg), &
                    &                       ' and lp2cintp_sfcana=',lp2cintp_sfcana(jg),' not allowed'
                    CALL finish(routine, message_text)
                ENDIF

                ! when vertical remapping of the FG-fields is applied, z_ifc is required 
                ! as FG input field
                IF (lvert_remap_fg) THEN
                  CALL add_to_list(fgGroup, fgGroupSize, 'z_ifc')
                ENDIF

            CASE(MODE_COMBINED,MODE_COSMO)
                ! add SMI to the default list
                ! I.e. ICON tries to read SMI, with W_SO being the fallback-option
                CALL add_to_list(fgGroup, fgGroupSize, 'smi')

                ! no analysis group
                anaGroupSize = 0

            CASE DEFAULT
                fgGroupSize = 0
                anaGroupSize = 0

        END SELECT


    END SUBROUTINE collectGroups


    !-------------
    !>
    !! SUBROUTINE readInstructionList_make
    !! Generates a mapping from icon variable names to a pair of logicals that define
    !! from which file each variable IS allowed to be READ.
    !!
    !! This IS based on:
    !!  1. the group specifications IN the add_var() calls
    !!  2. user input via initicon_config(:)%ana_checklist,
    !!     ANY of these variables must be READ from the analysis file
    !!  3. user input via initicon_config(:)%fg_checklist,
    !!     ANY of these first guess fields must be read from the first guess file.
    !!     I.e. optional first guess fields are turned into mandatory ones
    !!  3. some other interesting rules, like the substitution of 'smi' for 'w_so'
    !!     IN the CASE of MODE_COSMO AND MODE_COMBINED.
    !!
    !! These are the variable groups which are used to determine from which file a variable IS READ:
    !!     MODE_DWDANA    : mode_dwd_fg_in + mode_dwd_ana_in
    !!     MODE_ICONVREMAP: mode_dwd_fg_in + mode_dwd_ana_in
    !!     MODE_IAU       : mode_iau_fg_in + mode_iau_ana_in - mode_iau_anaatm_in
    !!     MODE_IAU_OLD   : mode_iau_old_fg_in + mode_iau_old_ana_in - mode_iau_anaatm_in
    !!     MODE_COMBINED  : mode_combined_in
    !!     MODE_COSMO     : mode_cosmo_in
    !!     MODE_IFSANA    : <NONE>
    !!
    !! In contrast to the old create_input_groups() SUBROUTINE, this does NOT check the contents of the files itself.
    !! It ONLY provides concise instructions of which fields to READ from which file.
    !!
    !! Returns a new instruction list object.
    FUNCTION readInstructionList_make(p_patch, init_mode) RESULT(resultVar)
        TYPE(t_patch), INTENT(IN) :: p_patch
        INTEGER, INTENT(IN) :: init_mode
        ! the resulting list of variable names to be READ together with flags defining which input file may be used
        TYPE(t_readInstructionList), POINTER :: resultVar

        ! local variables
        CHARACTER(LEN = *), PARAMETER :: routine = modname//':readInstructionList_make'
        INTEGER :: ivar, error
        TYPE(t_readInstruction), pointer :: curInstruction

        ! lists of variable names
        CHARACTER(LEN=vname_len), DIMENSION(1000) :: grp_vars_fg, grp_vars_optfg, grp_vars_ana

        ! the corresponding sizes of the lists above
        INTEGER :: ngrp_vars_fg, ngrp_vars_optfg, ngrp_vars_ana
        !
        !-------------------------------------------------------------------------

        ! create a list of instructions according to the current init_mode AND configuration flags
        CALL collectGroups(p_patch, init_mode, grp_vars_fg   , ngrp_vars_fg   , &
          &                                    grp_vars_optfg, ngrp_vars_optfg, &
          &                                    grp_vars_ana  , ngrp_vars_ana    )
        ALLOCATE(resultVar, STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure");
        CALL resultVar%construct()
        DO ivar = 1, ngrp_vars_fg
            curInstruction => resultVar%findInstruction(grp_vars_fg(ivar))
            curInstruction%lReadFg = .TRUE.
            !
            ! mark optional first guess fields in the instruction list. 
            IF (one_of(TRIM(grp_vars_fg(ivar)), grp_vars_optfg)/=-1) THEN
                WRITE(message_text,'(a,a,a,i2)') 'Declare ',TRIM(grp_vars_fg(ivar)),' as OPTIONAL for DOM ', p_patch%id
                CALL message(routine, message_text)
                curInstruction%lOptionalFg=.TRUE.
            ENDIF
        END DO
        DO ivar = 1, ngrp_vars_ana
            curInstruction => resultVar%findInstruction(grp_vars_ana(ivar))
            curInstruction%lReadFg = .TRUE.
            curInstruction%lReadAna = .TRUE.
        END DO

        ! Allow the user to override the DEFAULT settings for optional fields via fg_checklist
        ! I.e. the user can change an optional field into a mandatory one and force a model abort, 
        ! if the field is not available as input.
        ! 
        ! translate GRIB2 varname to internal netcdf varname
        ! If requested GRIB2 varname is not found in the dictionary
        ! (i.e. due to typos) -> Model abort
        DO ivar=1,SIZE(initicon_config(p_patch%id)%fg_checklist)
            IF (initicon_config(p_patch%id)%fg_checklist(ivar) == ' ') EXIT

            curInstruction => resultVar%findInstruction(TRIM(ana_varnames_dict%get( &
                 &                                                  initicon_config(p_patch%id)%fg_checklist(ivar), &
                 &                                                  linverse=.TRUE.)), opt_expand=.FALSE.)
            ! Note that depending on the Namelist settings, not every field listed in 
            ! fg_checklist is part of the instruction list. Therefore, curInstruction 
            ! may be non-associated. 
            IF (ASSOCIATED(curInstruction) .AND. curInstruction%lOptionalFg) THEN
                curInstruction%lOptionalFg = .FALSE.

                WRITE(message_text,'(a,a,a,i2)') 'Transform ',TRIM(ana_varnames_dict%get( &
                 &                                initicon_config(p_patch%id)%fg_checklist(ivar), linverse=.TRUE.)), &
                 &                               ' into a mandatory first guess field for DOM', p_patch%id
                CALL message(routine, message_text)
            ENDIF
        ENDDO

        ! Allow the user to override the DEFAULT settings via ana_checklist. These variables must be READ from analysis.
        IF( lread_ana ) THEN
            ! translate GRIB2 varname to internal netcdf varname
            ! If requested GRIB2 varname is not found in the dictionary
            ! (i.e. due to typos) -> Model abort
            DO ivar=1,SIZE(initicon_config(p_patch%id)%ana_checklist)
                IF (initicon_config(p_patch%id)%ana_checklist(ivar) == ' ') EXIT

                curInstruction => resultVar%findInstruction(TRIM(ana_varnames_dict%get( &
                &                                                      initicon_config(p_patch%id)%ana_checklist(ivar), &
                &                                                      linverse=.TRUE.)))
                curInstruction%lReadAna = .TRUE.
                curInstruction%lRequireAna = .TRUE.
            ENDDO
        END IF

    END FUNCTION readInstructionList_make

    SUBROUTINE readInstructionList_construct(me)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        INTEGER :: error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_construct"

        ALLOCATE(me%list(8), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")
        me%nInstructions = 0
    END SUBROUTINE readInstructionList_construct

    FUNCTION readInstructionList_findInstruction(me, varName, opt_expand) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        LOGICAL, INTENT(IN), OPTIONAL  :: opt_expand  ! TRUE: expand list, if variable is not found 
        TYPE(t_readInstruction), POINTER :: resultVar

        INTEGER :: i
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_findInstruction"
        LOGICAL :: expand   ! TRUE/FALSE: expand/do not expand list, if variable is not found  

        IF (PRESENT(opt_expand)) THEN
            expand = opt_expand
        ELSE
            expand = .TRUE. 
        ENDIF
        ! initialize
        resultVar => NULL() 

        ! try to find it IN the current list
        DO i = 1, me%nInstructions
            resultVar => me%list(i)
            IF(resultVar%varName == varName) RETURN
        END DO

        IF (expand) THEN
            ! the variable IS NOT found IN the list, expand the list
            ! first increase the allocation IF necessary
            IF(me%nInstructions == SIZE(me%list, 1)) CALL me%resize(2*me%nInstructions)
            ! THEN add a new entry at the END
            me%nInstructions = me%nInstructions + 1
            resultVar => me%list(me%nInstructions)
            resultVar%varName = varName
            resultVar%lReadFg = .FALSE.
            resultVar%lOptionalFg = .FALSE.
            resultVar%lReadAna = .FALSE.
            resultVar%lRequireAna = .FALSE.
            resultVar%statusFg = kStateNoFetch
            resultVar%statusAna = kStateNoFetch
            resultVar%sourceOverride = kInputSourceUnset
        ENDIF
    END FUNCTION readInstructionList_findInstruction

    SUBROUTINE readInstructionList_resize(me, newSize)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: newSize

        TYPE(t_readInstruction), POINTER :: tempList(:)
        INTEGER :: i, error
        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_resize"

        ! just IN CASE we are actually shrinking the list, avoid overrunning the new buffer
        IF(newSize < me%nInstructions) me%nInstructions = newSize

        ! get a new allocation
        ALLOCATE(tempList(newSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

        ! copy the list contents over
        DO i = 1, me%nInstructions
            tempList(i)%varName = me%list(i)%varName
            tempList(i)%lReadFg = me%list(i)%lReadFg
            tempList(i)%lOptionalFg = me%list(i)%lOptionalFg
            tempList(i)%lReadAna = me%list(i)%lReadAna
            tempList(i)%lRequireAna = me%list(i)%lRequireAna
            tempList(i)%statusFg = me%list(i)%statusFg
            tempList(i)%statusAna = me%list(i)%statusAna
            tempList(i)%sourceOverride = me%list(i)%sourceOverride
        END DO

        ! replace the old list with the new one
        DEALLOCATE(me%list)
        me%list => tempList
        tempList => NULL()
    END SUBROUTINE readInstructionList_resize

    SUBROUTINE readInstructionList_fileRequests(me, requestList, lIsFg)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        TYPE(t_inputRequestList), INTENT(INOUT) :: requestList
        LOGICAL, INTENT(IN) :: lIsFg

        INTEGER :: i

        DO i = 1, me%nInstructions
            IF(lIsFg) THEN
                IF(me%list(i)%lReadFg) CALL requestList%request(TRIM(me%list(i)%varName))
            ELSE
                IF(me%list(i)%lReadAna) CALL requestList%request(TRIM(me%list(i)%varName))
            END IF
        END DO
    END SUBROUTINE readInstructionList_fileRequests

    LOGICAL FUNCTION readInstructionList_wantVar(me, varName, lIsFg) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        LOGICAL, INTENT(IN) :: lIsFg

        IF(lIsFg) THEN
            resultVar = me%wantVarFg(varName)
        ELSE
            resultVar = me%wantVarAna(varName)
        END IF
    END FUNCTION readInstructionList_wantVar

    LOGICAL FUNCTION readInstructionList_wantVarFg(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName

        INTEGER :: i

        resultVar = .FALSE.
        DO i = 1, me%nInstructions
            IF(me%list(i)%varName == varName) THEN
                resultVar = me%list(i)%lReadFg
                RETURN
            END IF
        END DO
    END FUNCTION readInstructionList_wantVarFg

    LOGICAL FUNCTION readInstructionList_wantVarAna(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(IN) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName

        INTEGER :: i

        resultVar = .FALSE.
        DO i = 1, me%nInstructions
            IF(me%list(i)%varName == varName) THEN
                resultVar = me%list(i)%lReadAna
                RETURN
            END IF
        END DO
    END FUNCTION readInstructionList_wantVarAna

    SUBROUTINE readInstructionList_handleError(me, lSuccess, varName, caller, lIsFg)
        LOGICAL, INTENT(IN) :: lSuccess, lIsFg
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        IF(lIsFg) THEN
            CALL me%handleErrorFg(lSuccess, varName, caller)
        ELSE
            CALL me%handleErrorAna(lSuccess, varName, caller)
        END IF
    END SUBROUTINE readInstructionList_handleError

    SUBROUTINE readInstructionList_handleErrorFg(me, lSuccess, varName, caller)
        LOGICAL, INTENT(IN) :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_handleErrorFg"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadFg) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        IF(lSuccess) THEN
            instruction%statusFg = kStateRead
        ELSE
            IF (instruction%lOptionalFg) THEN
                ! field is not necessary for starting the model
                ! Therefore we set the status to kStateFailedOptFetch 
                ! instead of kStateFailedFetch
                instruction%statusFg = kStateFailedOptFetch
            ELSE
                instruction%statusFg = kStateFailedFetch

                ! check whether we can READ this variable from the analysis file as well
                IF(instruction%lReadAna) THEN
                    ! now this variable must be READ from the analysis file
                    instruction%lReadFg = .FALSE.
                ELSE
                    CALL finish(caller, "failed to read variable '"//varName//"' from the first guess file, &
                    &and reading from analysis file is not allowed")
                END IF
            ENDIF 
        END IF
    END SUBROUTINE readInstructionList_handleErrorFg

    SUBROUTINE readInstructionList_handleErrorAna(me, lSuccess, varName, caller)
        LOGICAL, INTENT(IN) :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_handleErrorAna"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadAna) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        IF(lSuccess) THEN
            instruction%statusAna = kStateRead
        ELSE
            instruction%statusAna = kStateFailedFetch

            ! check whether this variable was already READ from the first guess file; IF NOT, we have to error OUT
            IF(.NOT.instruction%lReadFg) THEN
                CALL finish(caller, "failed to read variable '"//varName//"' from the analysis file, &
                    &and reading from first guess is not allowed or failed")
            END IF

            ! check whether we are forced by the user to READ analysis DATA for this variable (ana_checklist)
            IF(instruction%lRequireAna) THEN
                CALL finish(caller, "failed to read variable '"//varName//"' from the analysis file (required by ana_checklist)")
            END IF
        END IF
    END SUBROUTINE readInstructionList_handleErrorAna


    SUBROUTINE readInstructionList_optionalReadResult(me, lSuccess, varName, caller, lIsFg)
        LOGICAL, INTENT(IN) :: lSuccess, lIsFg
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        IF(lIsFg) THEN
            CALL me%optionalReadResultFg(lSuccess, varName, caller)
        ELSE
            CALL me%optionalReadResultAna(lSuccess, varName, caller)
        END IF
    END SUBROUTINE readInstructionList_optionalReadResult

    SUBROUTINE readInstructionList_optionalReadResultFg(me, lSuccess, varName, caller)
        LOGICAL, INTENT(IN) :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_optionalReadResultFg"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadFg) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        instruction%statusFg = MERGE(kStateRead, kStateFailedFetch, lsuccess)
    END SUBROUTINE readInstructionList_optionalReadResultFg

    SUBROUTINE readInstructionList_optionalReadResultAna(me, lSuccess, varName, caller)
        LOGICAL, INTENT(IN) :: lSuccess
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName, caller

        CHARACTER(LEN = *), PARAMETER :: routine = modname//"readInstructionList_optionalReadResultAna"
        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)

        ! sanity check
        IF(.NOT.instruction%lReadAna) CALL finish(routine, "internal error: variable '"//varName//"' was read even though there &
            &are no read instructions for it")

        instruction%statusAna = MERGE(kStateRead, kStateFailedFetch, lSuccess)
    END SUBROUTINE readInstructionList_optionalReadResultAna


    INTEGER FUNCTION readInstruction_source(me) RESULT(resultVar)
        CLASS(t_readInstruction), INTENT(IN) :: me

        IF(me%sourceOverride /= kInputSourceUnset) THEN
            resultVar = me%sourceOverride
        ELSE IF (me%statusAna == kStateRead) THEN
          resultVar = kInputSourceAna
        ELSE IF (me%statusFg == kStateRead) THEN
          resultVar = kInputSourceFg
        ELSE IF (me%statusFg == kStateFailedOptFetch) THEN
          resultVar = kInputSourceCold
        ELSE
          resultVar = kInputSourceNone
        END IF
    END FUNCTION readInstruction_source

    INTEGER FUNCTION readInstructionList_sourceOfVar(me, varName) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName

        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)
        resultVar = instruction%source()
    END FUNCTION readInstructionList_sourceOfVar

    SUBROUTINE readInstructionList_setSource(me, varName, source)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: source

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_setSource"
        TYPE(t_readInstruction), POINTER :: instruction

        SELECT CASE(source)
            CASE(kInputSourceNone, kInputSourceFg, kInputSourceAna, kInputSourceBoth, &
              &  kInputSourceCold, kInputSourceAnaI, kInputSourceFgAnaI)

                instruction => me%findInstruction(varName)
                instruction%sourceOverride = source
            CASE DEFAULT
                CALL finish(routine, "assertion failed: illegal source argument")
        END SELECT
    END SUBROUTINE readInstructionList_setSource


    INTEGER FUNCTION readInstructionList_fetchStatus(me, varName, lIsFg) RESULT(resultVar)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: varName
        LOGICAL, INTENT(IN) :: lIsFg

        TYPE(t_readInstruction), POINTER :: instruction

        instruction => me%findInstruction(varName)
        IF(lIsFg) THEN
            resultVar = instruction%statusFg
        ELSE
            resultVar = instruction%statusAna
        ENDIF
    END FUNCTION readInstructionList_fetchStatus




    ! The table that is printed by this function deliberately depends on the actual read attempts and their results, 
    ! not on `lReadFg` or `lReadAna`. This is due to the fact that there are existing discrepancies between the 
    ! input groups and the actual read attempts made by the `fetch...()` routines in `mo_initicon_io`: 
    ! The table is supposed to show the reality of which data was read from where, and which inputs were used,
    ! not some hypothetical this-is-what-should-have-been-done info.
    SUBROUTINE readInstructionList_printSummary(me, jg)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me
        INTEGER, INTENT(IN) :: jg

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":readInstructionList_printSummary"
        INTEGER :: i
        TYPE(t_readInstruction), POINTER :: curInstruction
        TYPE(t_table) :: table
        CHARACTER(LEN = *), PARAMETER :: variableCol = "variable", &
                                       & fgAttemptCol = "FG read attempt", &
                                       & fgSuccessCol = "FG data found", &
                                       & anaAttemptCol = "ANA read attempt", &
                                       & anaSuccessCol = "ANA data found", &
                                       & anaRequiredCol = "ANA DATA required (ana_checklist)", &
                                       & useCol = "data used from"

        IF(me%nInstructions == 0) THEN
            ! Don't print empty tables, just give a message that there are no input instructions.
            WRITE(0, '(a,i0)') "no input results available for domain ", jg
            RETURN
        END IF

        ! Print the title of the list.
        WRITE(0,*) ""
        WRITE(0, '(a,i0,a)') "input results for domain ", jg, ":"

        CALL initialize_table(table)
        CALL add_table_column(table, variableCol)
        CALL add_table_column(table, fgAttemptCol)
        CALL add_table_column(table, fgSuccessCol)
        CALL add_table_column(table, anaAttemptCol)
        CALL add_table_column(table, anaSuccessCol)
        CALL add_table_column(table, anaRequiredCol)
        CALL add_table_column(table, useCol)

        DO i = 1, me%nInstructions
            curInstruction => me%list(i)

            CALL set_table_entry(table, i, variableCol, TRIM(curInstruction%varName))

            ! Display whether an actual read attempt from FG has been made, and its result.
            ! read attempt == there was a call from a `fetch...()` routine in `mo_initicon_io`
            SELECT CASE(curInstruction%statusFg)
                CASE(kStateNoFetch)
                    CALL set_table_entry(table, i, fgAttemptCol, "no")
                CASE(kStateFailedFetch, kStateFailedOptFetch)
                    CALL set_table_entry(table, i, fgAttemptCol, "yes")
                    CALL set_table_entry(table, i, fgSuccessCol, "no")
                CASE(kStateRead)
                    CALL set_table_entry(table, i, fgAttemptCol, "yes")
                    CALL set_table_entry(table, i, fgSuccessCol, "yes")
                CASE DEFAULT
                    CALL finish(routine, "unexpected value in statusFg")
            END SELECT

            ! As above, but for ANA.
            SELECT CASE(curInstruction%statusAna)
                CASE(kStateNoFetch)
                    CALL set_table_entry(table, i, anaAttemptCol, "no")
                CASE(kStateFailedFetch, kStateFailedOptFetch)
                    CALL set_table_entry(table, i, anaAttemptCol, "yes")
                    CALL set_table_entry(table, i, anaSuccessCol, "no")
                CASE(kStateRead)
                    CALL set_table_entry(table, i, anaAttemptCol, "yes")
                    CALL set_table_entry(table, i, anaSuccessCol, "yes")
                CASE DEFAULT
                    CALL finish(routine, "unexpected value in statusAna")
            END SELECT

            IF(curInstruction%lRequireAna) THEN
                CALL set_table_entry(table, i, anaRequiredCol, "required")
            END IF

            ! Display which data is actually used. If there was a call to `setSource()`, that is what is printed,
            ! otherwise this column is computed from the two `statusXXX` variables printed above.
            SELECT CASE(curInstruction%source())
                CASE(kInputSourceNone)
                    CALL set_table_entry(table, i, useCol, "none")
                CASE(kInputSourceCold)
                    CALL set_table_entry(table, i, useCol, "cold!")
                CASE(kInputSourceFg)
                    CALL set_table_entry(table, i, useCol, "fg")
                CASE(kInputSourceAna)
                    CALL set_table_entry(table, i, useCol, "ana")
                CASE(kInputSourceBoth)
                    CALL set_table_entry(table, i, useCol, "both")
                CASE(kInputSourceAnaI)
                    CALL set_table_entry(table, i, useCol, "ana(intp)")
                CASE(kInputSourceFgAnaI)
                    CALL set_table_entry(table, i, useCol, "fg,ana(intp)")
                CASE DEFAULT
                    CALL finish(routine, "unexpected RESULT from curInstruction%source()")
            END SELECT
        END DO

        CALL print_table(table, opt_delimiter = " | ")

        CALL finalize_table(table)

        WRITE(0,*) ""   ! Separate the table from the following messages.
    END SUBROUTINE readInstructionList_printSummary

    SUBROUTINE readInstructionList_destruct(me)
        CLASS(t_readInstructionList), INTENT(INOUT) :: me

        DEALLOCATE(me%list)
    END SUBROUTINE readInstructionList_destruct

END MODULE mo_input_instructions
