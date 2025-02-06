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

! This module contains the I/O routines for initicon

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_initicon_io

  USE mo_kind,                ONLY: wp, dp
  USE mo_io_units,            ONLY: filename_max
  USE mo_parallel_config,     ONLY: nproma
  USE mo_run_config,          ONLY: msg_level, iqv, iqc, iqi, iqr, iqs, iqg, &
                                    iqh, iqnc, iqni, iqnr, iqns, iqng, iqnh
  USE mo_dynamics_config,     ONLY: nnow, nnow_rcf
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_state, t_nh_prog, t_nh_diag
  USE mo_nwp_phy_types,       ONLY: t_nwp_phy_diag,t_nwp_phy_stochconv
  USE mo_nwp_lnd_types,       ONLY: t_lnd_state, t_lnd_prog, t_lnd_diag, t_wtr_prog
  USE mo_initicon_types,      ONLY: t_initicon_state, t_pi_atm_in, t_pi_atm, &
    &                               alb_snow_var, geop_ml_var
  USE mo_input_instructions,  ONLY: t_readInstructionListPtr, kInputSourceFg, &
    &                               kInputSourceAna, kInputSourceBoth, kStateFailedFetch, &
    &                               kInputSourceCold
  USE mo_initicon_config,     ONLY: init_mode, l_sst_in, generate_filename,             &
    &                               ifs2icon_filename, lread_vn, lread_tke,             &
    &                               lp2cintp_incr, lp2cintp_sfcana, ltile_coldstart,    &
    &                               lvert_remap_fg, aerosol_fg_present, nlevsoil_in, qcana_mode, qiana_mode, qrsgana_mode, &
    &                               qnxana_2mom_mode, icpl_da_sfcevap, icpl_da_skinc, icpl_da_snowalb, icpl_da_sfcfric
  USE mo_nh_init_nest_utils,  ONLY: interpolate_scal_increments, interpolate_sfcana
  USE mo_nh_init_utils,       ONLY: convert_omega2w, compute_input_pressure_and_height
  USE mo_impl_constants,      ONLY: max_dom, MODE_ICONVREMAP,          &
    &                               MODE_IAU, MODE_IAU_OLD, MODE_IFSANA, MODE_COMBINED, &
    &                               MODE_COSMO, iss, iorg, ibc, iso4, idu, SUCCESS,     &
    &                               vname_len
  USE mo_exception,           ONLY: message, finish, message_text, warning
  USE mo_grid_config,         ONLY: n_dom, nroot, l_limited_area
  USE mo_mpi,                 ONLY: p_io, p_bcast, p_comm_work,    &
    &                               my_process_is_mpi_workroot,    &
    &                               my_process_is_stdio
  USE mo_io_config,           ONLY: default_read_method
  USE mo_limarea_config,      ONLY: latbc_config
  USE mo_read_interface,      ONLY: t_stream_id, openInputFile, closeFile, on_cells, on_edges, &
    &                               read_2d_1time, read_2d_1lev_1time, read_3d_1time
  USE mo_netcdf_errhandler,   ONLY: nf
  USE mo_netcdf
  USE mo_nwp_sfc_tiles,       ONLY: t_tileinfo_icon, trivial_tile_att
  USE mo_lnd_nwp_config,      ONLY: ntiles_total,  l2lay_rho_snow, &
    &                               ntiles_water, lmulti_snow, lsnowtile, &
    &                               isub_lake, isub_water, llake, lprog_albsi, itype_trvg, &
    &                               itype_snowevap, itype_canopy, nlev_soil
  USE mo_master_config,       ONLY: getModelBaseDir
  USE mo_nwp_sfc_interp,      ONLY: smi_to_wsoil
  USE mo_initicon_utils,      ONLY: allocate_extana_atm, allocate_extana_sfc, &
       &                            init_qnx_from_qx_twomom, init_qnxinc_from_qxinc_twomom
  USE mo_physical_constants,  ONLY: cpd, rd, cvd_o_rd, p0ref, vtmpc1, tmelt
  USE mo_fortran_tools,       ONLY: init
  USE mo_input_request_list,  ONLY: t_InputRequestList
  USE mo_util_string,         ONLY: int2string
  USE mo_atm_phy_nwp_config,  ONLY: iprog_aero, atm_phy_nwp_config
  USE mo_var_metadata_types,  ONLY: t_var_metadata
  USE mo_var_list_register,   ONLY: t_vl_register_iter
  USE mo_var,                 ONLY: level_type_ml
  USE mo_var_groups,          ONLY: var_groups_dyn
  USE mo_var_metadata,        ONLY: get_var_name

  ! High level overview of how `mo_initicon` reads input data
  ! =========================================================
  !
  ! Knobs that control what is read
  ! -------------------------------
  !
  !  1. variable groups
  !
  !  2. init_mode
  !
  !  3. further flags and the `ana_checklist` and 'fg_checklist' fields from the namelists
  !
  !
  ! Steps of reading
  ! ----------------
  !
  !  1. A list with input instructions is generated.
  !     The code for this is found in `mo_input_instructions`.
  !
  !     This is where the variable groups enter the process, all further processing depends on the variable groups only indirectly via `t_ReadInstructions`.
  !     This is also the place where the `ana_checklist` and `fg_checklist` namelist parameter is evaluated.
  !
  !     The input instructions serve a dual purpose (which is a source of confusion, unfortunately):
  !     They encode from which input file a read attempt for a variable should be made,
  !     but they also collect whether such an attempt has been made, and whether that attempt was successfull.
  !
  !  2. The input instructions are used to generate a `t_InputRequestList`.
  !     This is done by the `fileRequests()` method of `t_ReadInstructionList`, and is quite straight-forward:
  !     All variables that can be read from a file are requested from that file.
  !
  !  3. The files are read, and a file inventory output is produced in the process (mo_initicon: read_dwdfg() and read_dwdana()).
  !     Note that the file inventory contains only information about variables that have been requested!
  !
  !  4. The data is fetched from the `t_InputRequestList`.
  !     This happens in `fetch_dwd...()` routines in `mo_initicon_io`.
  !
  !     This is a complicated process that takes much more conditions into account then the code generating the `t_ReadInstructions` does.
  !     Especially, this is the place where all the different namelist flags are honored, and it's the place where optional reading is handled
  !     (variables that are read if they are present in the file, but which are not necessary to start the run).
  !
  !     The `fetch_dwd...()` routines generally first ask the `t_ReadInstructions` whether they should attempt to read a variable,
  !     then they try to fetch the data from the `t_InputRequestList`, and finally inform the `t_ReadInstructions` about the result of this operation.
  !     It is then the duty of `t_ReadInstructions` to check whether a failure to read is fatal or whether it can be compensated.
  !
  !     This ask-fetch-inform triple is usually encapsulated within the wrapper functions `fetch2d()` etc. within this file.
  !     Exceptions exist, like the all-or-nothing behavior that's implemented for the aerosol variables, but they are not frequent.
  !
  !  5. The `t_ReadInstructions` is asked to output a table describing what read attempts have been made
  !     and which data is actually used for the model run.
  !
  !     Note that this is about what read attempts have *really* been made (= there was an ask-fetch-inform triple for this variable),
  !     not what read attempts *should* have been made (= `t_ReadInstructions::wantVar()` would have returned `.true.` if it had been called).
  !
  !
  ! How to add variables
  ! --------------------
  !
  ! There are two places that need to be changed to add support for a new input variable:
  !
  !  1. The variable has to be added to the relevant variable groups, so that input instructions are generated for it, and so that it's requested from the right files.
  !
  !  2. If the new variable is an optional first guess field, it has to be added manually in the SUBROUTINE collectGroupFgOpt 
  !     in mo_input_instructions
  !
  !  2. Code has to be added to one of the `fetch_dwd...()` routines in this file to retrieve the data from the file(s)
  !     and to inform the `t_ReadInstructions` about the succes/failure to do so.
  !
  !     In most cases, this is simply done by calling the right wrapper routine (`fetch2d()` and friends).
  !     If the appropriate wrapper routine is missing, just copy-paste-modify one of the existing wrappers;
  !     I have only generated the ones that were actually needed to implement the current functionality, so the list of wrappers is incomplete.

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_initicon_io'


  PUBLIC :: read_extana_atm
  PUBLIC :: read_extana_sfc

  PUBLIC :: fetch_dwdfg_atm
  PUBLIC :: fetch_dwdfg_atm_ii
  PUBLIC :: fetch_dwdfg_sfc
  PUBLIC :: fetch_dwdana_atm
  PUBLIC :: fetch_dwdana_sfc

! PUBLIC :: process_input_dwdfg_atm     !This is currently not needed.
! PUBLIC :: process_input_dwdfg_atm_ii
  PUBLIC :: process_input_dwdfg_sfc
  PUBLIC :: process_input_dwdana_atm
  PUBLIC :: process_input_dwdana_sfc

  PUBLIC :: height_or_lev


  TYPE :: t_fetchParams
    TYPE(t_readInstructionListPtr), ALLOCATABLE :: inputInstructions(:)
    CLASS(t_InputRequestList), POINTER :: requestList
    CHARACTER(LEN = :), ALLOCATABLE :: routine
    LOGICAL :: isFg
  END TYPE t_fetchParams

  CONTAINS


  ! Hack to determine the dimension name for phase2 simulations.
  ! We now also determine the number of height levels
  SUBROUTINE height_or_lev (ncid, dimid, nlev)
    INTEGER, INTENT(IN   ) :: ncid
    INTEGER, INTENT(  OUT) :: dimid
    INTEGER, INTENT(  OUT) :: nlev

    INTEGER, PARAMETER :: max_n_height = 9
    INTEGER :: nlevs(max_n_height)
    INTEGER :: dimids(max_n_height)
    INTEGER :: retval
    INTEGER :: i
    CHARACTER(len=23) :: dimstring

    CHARACTER(len=*), PARAMETER :: routine = modname//"height_or_lev"

    retval = nf90_inq_dimid(ncid, 'lev', dimid)

    ! the "lev" branch
    IF (retval == nf90_noerr) THEN
      CALL nf(nf90_inquire_dimension(ncid, dimid, len = nlev), routine)
    ! the "hate*" branch
    ELSE
      dimstring = "height"
      DO i = 1, max_n_height
        retval = nf90_inq_dimid(ncid, TRIM(dimstring), dimids(i))
        IF (retval == nf90_noerr) THEN
          CALL nf(nf90_inquire_dimension(ncid, dimids(i), len = nlevs(i)), routine)
        ELSE
          nlevs(i) = HUGE(1)
        ENDIF

        WRITE(dimstring, "(A,I0)") "height_",i+1
      ENDDO
      i     = MINLOC(nlevs,1)
      nlev  = nlevs(i)
      dimid = dimids(i)
    ENDIF

  END SUBROUTINE height_or_lev




  !>
  !! Read horizontally interpolated external analysis (atmosphere only)
  !!
  !! Reads horizontally interpolated external analysis atmosphere data
  !! (currently IFS or COSMO) and reads in vertical coordinate table.
  !!
  SUBROUTINE read_extana_atm (p_patch, initicon)

    TYPE(t_patch), TARGET,  INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    INTEGER :: jg, jlev, jc, jk, jb, i_endidx, jg1
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels, nlev_in, nhyi
    INTEGER :: ncid, dimid, varid, ierrstat
    TYPE(t_stream_id) :: stream_id
    INTEGER :: psvar_ndims, geopvar_ndims, itemp(7)

    REAL(wp), ALLOCATABLE               :: psfc(:,:), phi_sfc(:,:), z_ifc_in(:,:,:), &
      &                                    w_ifc(:,:,:), omega(:,:,:)

    CHARACTER(LEN=10) :: psvar

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_extana_atm'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)
    LOGICAL :: lread_process
    LOGICAL :: lread_qr, lread_qs ! are qr, qs provided as input?
    !-------------------------------------------------------------------------

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()
    nlev_in = 0

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level
      ifs2icon_file(jg) = " "

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      !
      ! generate file name
      !
      ifs2icon_file(jg) = generate_filename(ifs2icon_filename, getModelBaseDir(), &
        &                                   nroot, jlev, jg)

      IF (my_process_is_stdio()) THEN
        ! consistency check: check for duplicate file names which may
        ! occur, for example, if the keyword pattern (namelist
        ! parameter) has been defined ambiguously by the user.
        DO jg1 = 1,(jg-1)
          IF (.NOT. p_patch(jg1)%ldom_active) CYCLE
          IF (ifs2icon_file(jg1) == ifs2icon_file(jg)) THEN
            CALL finish(routine, "Error! Namelist parameter ifs2icon_filename="//TRIM(ifs2icon_filename)//&
              &"has been defined ambiguously for domains "//TRIM(int2string(jg1, '(i0)'))//" and "//&
              &TRIM(int2string(jg, '(i0)'))//"!")
          END IF
        END DO
      END IF

    END DO


    DO jg = 1, n_dom

      jlev = p_patch(jg)%level

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! Read in data from IFS2ICON
      !
      IF( lread_process ) THEN

        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(routine,'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf90_open(TRIM(ifs2icon_file(jg)), NF90_NOWRITE, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf90_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf90_inquire_dimension(ncid, dimid, len = no_cells), routine)

        !
        ! get number of vertical levels
        !
        CALL height_or_lev(ncid, dimid, no_levels)

        !
        ! check the number of cells
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(routine,&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        !
        ! DEFINE the number of vertical levels
        !
        IF(nlev_in /= 0)  THEN
          IF (nlev_in /= no_levels) THEN
            CALL finish(routine, 'nlev_in has already been defined differently!')
          END IF
        ELSE
          nlev_in = no_levels
        END IF


        !
        ! Check if rain water (QR) is provided as input
        !
        IF (nf90_inq_varid(ncid, 'QR', varid) == nf90_noerr) THEN
          lread_qr = .true.
        ELSE
          lread_qr = .false.
          CALL message(routine, 'Rain water (QR) not available in input data')
        ENDIF

        !
        ! Check if snow water (QS) is provided as input
        !
        IF (nf90_inq_varid(ncid, 'QS', varid) == nf90_noerr) THEN
          lread_qs = .true.
        ELSE
          lread_qs = .false.
          CALL message(routine, 'Snow water (QS) not available in input data')
        ENDIF

        !
        ! Check if normal velocity component (VN) is provided as input
        !
        lread_vn = nf90_inq_varid(ncid, 'VN', varid) == nf90_noerr


        IF (init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED) THEN
          !
          ! get number of hybrid coefficients
          !
          CALL nf(nf90_inq_dimid(ncid, 'nhyi', dimid), routine)
          CALL nf(nf90_inquire_dimension(ncid, dimid, len = nhyi), routine)
          !
          ! Check if surface pressure (PS) or its logarithm (LNPS) is provided as input
          !
          IF (nf90_inq_varid(ncid, 'PS', varid) == nf90_noerr) THEN
            psvar = 'PS'
          ELSE IF (nf90_inq_varid(ncid, 'LNPS', varid) == nf90_noerr) THEN
            psvar = 'LNPS'
          ENDIF

          ! Find out the dimension of psvar for reading purpose
          IF (nf90_inq_varid(ncid, TRIM(psvar), varid) == nf90_noerr) THEN
            CALL nf(nf90_inquire_variable(ncid, varid, ndims = psvar_ndims), routine)
          ELSE
            CALL finish(routine, 'surface pressure var '//TRIM(psvar)//' is missing')
          ENDIF

          !
          ! Check if model-level surface Geopotential is provided as GEOSP or GEOP_ML
          !
          IF (nf90_inq_varid(ncid, 'GEOSP', varid) == nf90_noerr) THEN
            geop_ml_var = 'GEOSP'
          ELSE IF (nf90_inq_varid(ncid, 'GEOP_ML', varid) == nf90_noerr) THEN
            geop_ml_var = 'GEOP_ML'
          ELSE
            CALL finish(routine,'Could not find model-level sfc geopotential')
          ENDIF

          ! Find out the dimension of geop_ml_var for reading purpose
          IF (nf90_inq_varid(ncid, TRIM(geop_ml_var), varid) == nf90_noerr) THEN
            CALL nf(nf90_inquire_variable(ncid, varid, ndims = geopvar_ndims), routine)
          ELSE
            CALL finish(routine,'surface geopotential var '//TRIM(geop_ml_var)//' is missing')
          ENDIF

        ENDIF

      ENDIF ! pe_io

      !
      ! open file (NetCDF file!)
      !
      CALL openInputFile(stream_id, ifs2icon_file(jg), p_patch(jg), &
        &                       default_read_method)

      itemp(1) = nlev_in
      itemp(2) = MERGE(1, 0, lread_qs)
      itemp(3) = MERGE(1, 0, lread_qr)
      itemp(4) = MERGE(1, 0, lread_vn)
      IF (init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED) THEN
        itemp(5) = nhyi
        itemp(6) = psvar_ndims
        itemp(7) = geopvar_ndims
      END IF

      CALL p_bcast(itemp, p_io, p_comm_work)

      nlev_in  = itemp(1)
      lread_qs = itemp(2) /= 0
      lread_qr = itemp(3) /= 0
      lread_vn = itemp(4) /= 0
      IF (init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED) THEN
        nhyi = itemp(5)
        psvar_ndims = itemp(6)
        geopvar_ndims = itemp(7)
      END IF

      IF (msg_level >= 10) THEN
        IF (init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED) THEN
          WRITE(message_text,'(2a)') 'surface pressure variable: ', TRIM(psvar)
          CALL message(routine, message_text)
          WRITE(message_text,'(2a)') 'Model-level surface geopotential: ', &
               TRIM(geop_ml_var)
          CALL message(routine, message_text)
        ENDIF
        IF (.NOT. lread_vn) THEN
          WRITE(message_text,'(a)') 'No direct input of vn! vn derived from (u,v).'
          CALL message(routine, message_text)
        ELSE
          WRITE(message_text,'(a)') 'Direct input of vn!'
          CALL message(routine, message_text)
        ENDIF
      ENDIF

      ! allocate data structure
      CALL allocate_extana_atm(nblks_c  = p_patch(jg)%nblks_c, &
        &                      nblks_e  = p_patch(jg)%nblks_e, &
        &                      nlev_in  = nlev_in,             &
        &                      atm_in   = initicon(jg)%atm_in, &
        &                      const    = initicon(jg)%const  ) !inout

      ! allocate local temporary arrays:
      ALLOCATE(psfc(nproma, p_patch(jg)%nblks_c), phi_sfc(nproma, p_patch(jg)%nblks_c))

      ! start reading atmospheric fields
      !
      CALL read_3d_1time(stream_id, on_cells, 'T', fill_array=initicon(jg)%atm_in%temp)

      IF (lread_vn) THEN
        CALL read_3d_1time(stream_id, on_edges, 'VN', fill_array=initicon(jg)%atm_in%vn)
      ELSE
        CALL read_3d_1time(stream_id, on_cells, 'U', fill_array=initicon(jg)%atm_in%u)

        CALL read_3d_1time(stream_id, on_cells, 'V', fill_array=initicon(jg)%atm_in%v)
      ENDIF


      IF (init_mode == MODE_COSMO) THEN
        ! allocate temporary array:
        ALLOCATE(w_ifc(nproma,nlev_in+1, p_patch(jg)%nblks_c))

        CALL read_3d_1time(stream_id, on_cells, 'W', fill_array=w_ifc)
      ELSE
        ! Note: in this case, input vertical velocity is in fact omega (Pa/s)

        ! allocate temporary array:
        ALLOCATE(omega(nproma,nlev_in,p_patch(jg)%nblks_c))

        CALL read_3d_1time(stream_id, on_cells, 'W', fill_array=omega)
      ENDIF

      CALL read_3d_1time(stream_id, on_cells, 'QV', fill_array=initicon(jg)%atm_in%qv)
      CALL read_3d_1time(stream_id, on_cells, 'QC', fill_array=initicon(jg)%atm_in%qc)
      CALL read_3d_1time(stream_id, on_cells, 'QI', fill_array=initicon(jg)%atm_in%qi)

      IF (lread_qr) THEN
        CALL read_3d_1time(stream_id, on_cells, 'QR', fill_array=initicon(jg)%atm_in%qr)
      ELSE
        initicon(jg)%atm_in%qr(:,:,:)=0._wp
      ENDIF

      IF (lread_qs) THEN
        CALL read_3d_1time(stream_id, on_cells, 'QS', fill_array=initicon(jg)%atm_in%qs)
      ELSE
        initicon(jg)%atm_in%qs(:,:,:)=0._wp
      ENDIF

      ! Allocate and read in vertical coordinate tables
      !
      ! Note that here the IFS input vertical grid is set up. This has
      ! to be distinguished from vct_a, vct_b, vct for the ICON
      ! vertical grid.
      !
      IF (init_mode == MODE_IFSANA .OR. init_mode == MODE_COMBINED) THEN

        ! read surface presure and (surface) geopotential
        IF (psvar_ndims==2)THEN
          CALL read_2d_1time(stream_id, on_cells, TRIM(psvar), &
            &                     fill_array=psfc)
        ELSEIF(psvar_ndims==3)THEN
          CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(psvar), &
            &                     fill_array=psfc)
        ELSE
          CALL finish(routine,'surface pressure var '//TRIM(psvar)//' dimension mismatch')
        END IF
        
        IF (geopvar_ndims==2)THEN
          CALL read_2d_1time(stream_id, on_cells, TRIM(geop_ml_var), &
            &                     fill_array=phi_sfc)
        ELSEIF(geopvar_ndims==3)THEN
          CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(geop_ml_var), &
            &                     fill_array=phi_sfc)
        ELSE
          CALL finish(routine,'surface geopotential var '//TRIM(geop_ml_var)//' dimension mismatch')
        END IF

        CALL initicon(jg)%const%vct%construct(ncid, p_io, p_comm_work)

      ELSE IF (init_mode == MODE_COSMO) THEN ! in case of COSMO-DE initial data
        
        ! allocate temporary array:
        ALLOCATE(z_ifc_in(nproma,nlev_in+1, p_patch(jg)%nblks_c), STAT=ierrstat)
        IF (ierrstat /= SUCCESS) CALL finish(routine, "ALLOCATE failed!")

        CALL read_3d_1time(stream_id, on_cells, 'HHL', fill_array=z_ifc_in)
        CALL read_3d_1time(stream_id, on_cells, 'P', fill_array=initicon(jg)%atm_in%pres)

        ! Interpolate input 'z3d' and 'w' from interface levels to main levels
!$OMP PARALLEL
!$OMP DO PRIVATE (jk,jc,jb,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1,p_patch(jg)%nblks_c

          i_endidx = MERGE(nproma, p_patch(jg)%npromz_c, &
               jb /= p_patch(jg)%nblks_c)

          DO jk = 1, nlev_in
            DO jc = 1, i_endidx

              initicon(jg)%const%z_mc_in(jc,jk,jb) = (z_ifc_in(jc,jk,jb) + z_ifc_in(jc,jk+1,jb)) * 0.5_wp
              initicon(jg)%atm_in%w(jc,jk,jb)      = (w_ifc(jc,jk,jb)    + w_ifc(jc,jk+1,jb)) * 0.5_wp
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL

      ELSE

        CALL finish(routine,'Incorrect init_mode')

      ENDIF ! init_mode = MODE_COSMO

      ! close file
      !
      IF(lread_process) CALL nf(nf90_close(ncid), routine)
      CALL closeFile(stream_id)

      IF (ANY(init_mode == (/ MODE_COMBINED, MODE_IFSANA /) )) THEN
        ! compute pressure and height of input data, using the IFS routines
        CALL compute_input_pressure_and_height(p_patch(jg), psfc, phi_sfc, initicon(jg))

        CALL convert_omega2w(omega, initicon(jg)%atm_in%w,                        &
          &                  initicon(jg)%atm_in%pres,  initicon(jg)%atm_in%temp, &
          &                  p_patch(jg)%nblks_c, p_patch(jg)%npromz_c,           &
          &                  initicon(jg)%atm_in%nlev )
      END IF

      ! cleanup
      IF (ALLOCATED(omega))    DEALLOCATE(omega)
      IF (ALLOCATED(psfc))     DEALLOCATE(psfc)
      IF (ALLOCATED(phi_sfc))  DEALLOCATE(phi_sfc)
      IF (ALLOCATED(z_ifc_in)) DEALLOCATE(z_ifc_in)
      IF (ALLOCATED(w_ifc))    DEALLOCATE(w_ifc)

    ENDDO ! loop over model domains

  END SUBROUTINE read_extana_atm





  !>
  !! Read horizontally interpolated external analysis (surface only)
  !!
  !! Reads horizontally interpolated external analysis surface data
  !! Currently, only IFS data are processed here
  !!
  SUBROUTINE read_extana_sfc (p_patch, initicon)

    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT) :: initicon(:)

    INTEGER :: jg, jlev, jg1
    LOGICAL :: l_exist

    INTEGER :: no_cells, no_levels
    INTEGER :: ncid, dimid, varid
    INTEGER :: geop_sfc_var_ndims     ! dimension of geop_sfc_var
    TYPE(t_stream_id) :: stream_id

    CHARACTER(LEN=10) :: geop_sfc_var ! surface-level surface geopotential

    CHARACTER(len=*), PARAMETER :: &
      routine = 'mo_nh_initicon:read_extana_sfc'

    CHARACTER(LEN=filename_max) :: ifs2icon_file(max_dom)
    LOGICAL :: lread_process

    !-------------------------------------------------------------------------

    ! flag. if true, then this PE reads data from file and broadcasts
    lread_process = my_process_is_mpi_workroot()

    DO jg = 1, n_dom

      jlev = p_patch(jg)%level
      ifs2icon_file(jg) = " "

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      !
      ! generate file name
      !
      ifs2icon_file(jg) = generate_filename(ifs2icon_filename, getModelBaseDir(), &
        &                                   nroot, jlev, jg)

      IF (my_process_is_stdio()) THEN
        ! consistency check: check for duplicate file names which may
        ! occur, for example, if the keyword pattern (namelist
        ! parameter) has been defined ambiguously by the user.
        DO jg1 = 1,(jg-1)
          IF (.NOT. p_patch(jg1)%ldom_active) CYCLE
          IF (ifs2icon_file(jg1) == ifs2icon_file(jg)) THEN
            CALL finish(routine, "Error! Namelist parameter ifs2icon_filename="//TRIM(ifs2icon_filename)//&
              &"has been defined ambiguously for domains "//TRIM(int2string(jg1, '(i0)'))//" and "//&
              &TRIM(int2string(jg, '(i0)'))//"!")
          END IF
        END DO
      END IF

    END DO


    DO jg = 1, n_dom

      jlev = p_patch(jg)%level

      ! Skip reading the atmospheric input data if a model domain
      ! is not active at initial time
      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      ! Read in data from IFS2ICON
      !
      IF(lread_process) THEN
        INQUIRE (FILE=ifs2icon_file(jg), EXIST=l_exist)
        IF (.NOT.l_exist) THEN
          CALL finish(routine,'IFS2ICON file is not found: '//TRIM(ifs2icon_file(jg)))
        ENDIF

        !
        ! open file
        !
        CALL nf(nf90_open(TRIM(ifs2icon_file(jg)), NF90_NOWRITE, ncid), routine)

        !
        ! get number of cells
        !
        CALL nf(nf90_inq_dimid(ncid, 'ncells', dimid), routine)
        CALL nf(nf90_inquire_dimension(ncid, dimid, len = no_cells), routine)

        !
        ! get number of vertical levels
        !
        CALL height_or_lev(ncid, dimid, no_levels)

        !
        ! check the number of cells and vertical levels
        !
        IF(p_patch(jg)%n_patch_cells_g /= no_cells) THEN
          CALL finish(routine,&
          & 'Number of patch cells and cells in IFS2ICON file do not match.')
        ENDIF

        ! Check, if the surface-level surface geopotential (GEOP_SFC) is available.
        ! If GEOP_SFC is missing, a warning will be issued and the model-level surface
        ! geopotential (GEOSP or GEOP_ML) will be used instead.
        IF (nf90_inq_varid(ncid, 'GEOP_SFC', varid) == nf90_noerr) THEN
          geop_sfc_var = 'GEOP_SFC'
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'surface-level surface geopotential is missing. ', &
            &  'use model-level surface geopotential, instead.'
          CALL message(routine, message_text)

          ! use model level geopotential instead
          geop_sfc_var = geop_ml_var
        ENDIF
        ! inquire dimension of geop_sfc_var which can be either 2 or 3
        IF (nf90_inq_varid(ncid, TRIM(geop_sfc_var), varid) == nf90_noerr) THEN
          CALL nf(nf90_inquire_variable(ncid, varid, ndims = geop_sfc_var_ndims), routine)
        ELSE
          CALL finish(routine,'surface geopotential variable '//TRIM(geop_sfc_var)//' is missing')
        ENDIF

        ! Check, if the snow albedo ALB_SNOW is available.
        ! If ALB_SNOW is missing, a warning will be issued and RHO_SNOW
        ! will be used instead to determine FRESHSNOW.
        IF (nf90_inq_varid(ncid, 'ALB_SNOW', varid) == nf90_noerr) THEN
          WRITE (message_text,'(a,a)')                            &
            &  'snow albedo available, ', &
            &  'used to determine freshsnow.'
          alb_snow_var = 'ALB_SNOW'
        ELSE

          WRITE (message_text,'(a,a)')                            &
            &  'snow albedo is missing. ', &
            &  'use snow density value, instead.'
          CALL message(routine, message_text)

          alb_snow_var = 'RHO_SNOW'
        ENDIF


        ! Check, if sea surface temperature field is provided as input
        IF (nf90_inq_varid(ncid, 'SST', varid) == nf90_noerr) THEN
          l_sst_in = .TRUE.
        ELSE
          WRITE (message_text,'(a,a)')                            &
            &  'sea surface temperature not available. ', &
            &  'initialize with skin temperature, instead.'
          CALL message(routine, message_text)
          l_sst_in = .FALSE.
        ENDIF

      ENDIF  ! p_io

      !
      ! open file
      !
      CALL openInputFile(stream_id, ifs2icon_file(jg), p_patch(jg), &
        &                       default_read_method)

      CALL p_bcast(l_sst_in, p_io, p_comm_work)

      CALL p_bcast(alb_snow_var, p_io, p_comm_work)

      CALL p_bcast(geop_sfc_var_ndims, p_io, p_comm_work)

      ! allocate data structure
      CALL allocate_extana_sfc(nblks_c     = p_patch(jg)%nblks_c,  &
        &                      nlevsoil_in = nlevsoil_in,          &
        &                      sfc_in      = initicon(jg)%sfc_in   ) !inout


      ! start reading surface fields
      !
      IF (geop_sfc_var_ndims == 3) THEN
        CALL read_2d_1lev_1time(stream_id, on_cells, TRIM(geop_sfc_var), &
          &                     fill_array=initicon(jg)%sfc_in%phi)
      ELSE IF (geop_sfc_var_ndims == 2) THEN
        CALL read_2d_1time(stream_id, on_cells, TRIM(geop_sfc_var), &
          &                     fill_array=initicon(jg)%sfc_in%phi)
      ELSE
        CALL finish(routine,"geop_sfc_var: Dimension mismatch")
      ENDIF

      CALL read_2d_1time(stream_id, on_cells, 'SKT', &
        &                fill_array=initicon(jg)%sfc_in%tskin)

      IF ( l_sst_in) THEN
        CALL read_2d_1time(stream_id, on_cells, 'SST', &
          &                fill_array=initicon(jg)%sfc_in%sst)
      ELSE
       initicon(jg)%sfc_in%sst(:,:)=0.0_wp
      END IF

      CALL read_2d_1time(stream_id, on_cells, 'T_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%tsnow)
      CALL read_2d_1time(stream_id, on_cells, TRIM(alb_snow_var), &
        &                fill_array=initicon(jg)%sfc_in%snowalb)
      CALL read_2d_1time(stream_id, on_cells, 'W_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%snowweq)
      CALL read_2d_1time(stream_id, on_cells, 'RHO_SNOW', &
        &                fill_array=initicon(jg)%sfc_in%snowdens)
      CALL read_2d_1time(stream_id, on_cells, 'W_I', &
        &                fill_array=initicon(jg)%sfc_in%skinres)
      CALL read_2d_1time(stream_id, on_cells, 'LSM', &
        &                fill_array=initicon(jg)%sfc_in%ls_mask)
      CALL read_2d_1time(stream_id, on_cells, 'CI', &
        &                fill_array=initicon(jg)%sfc_in%seaice)
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL1', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,1))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL2', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,2))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL3', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,3))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'STL4', &
        &                     fill_array=initicon(jg)%sfc_in%tsoil(:,:,4))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL1', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,1))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL2', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,2))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL3', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,3))
      CALL read_2d_1lev_1time(stream_id, on_cells, 'SMIL4', &
        &                     fill_array=initicon(jg)%sfc_in%wsoil(:,:,4))

      ! Checks
      IF(ANY(initicon(jg)%sfc_in%snowweq < 0._wp)) THEN
        CALL warning(routine,'W_SNOW contains negative values: Setting to Zero')
        initicon(jg)%sfc_in%snowweq =  &
          & MERGE(initicon(jg)%sfc_in%snowweq, 0._wp, (initicon(jg)%sfc_in%snowweq >= 0._wp))
      ENDIF

      ! close file
      !
      IF(lread_process) CALL nf(nf90_close(ncid), routine)
      CALL closeFile(stream_id)


    ENDDO ! loop over model domains

  END SUBROUTINE read_extana_sfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  !! Some wrapper routines to fetch DATA from an InputRequestList that first ask the inputInstructions whether a READ attempt should be made,
  !! AND proceed to tell the inputInstructions whether that attempt was successful.
  !! In the CASE of tiled input, these also take care of checking the ltile_coldstart flag, falling back to reading AND copying untiled input DATA IF it's set.

  SUBROUTINE fetch2d(params, varName, level, jg, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    REAL(dp), INTENT(in) :: level
    INTEGER, INTENT(in) :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:)
    TYPE(t_tileinfo_icon) :: tileinfo_icon
    LOGICAL :: fetchResult

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
      tileinfo_icon = trivial_tile_att%getTileinfo_icon()
      fetchResult = params%requestList%fetch2d(varName, level, tileinfo_icon%idx, jg, field)
      CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetch2d

  SUBROUTINE fetch3d(params, varName, jg, field, found)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, INTENT(in) :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:,:)
    LOGICAL, INTENT(OUT), OPTIONAL :: found ! allows to do the error handling in the calling routine

    LOGICAL :: fetchResult
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        tileinfo_icon = trivial_tile_att%getTileinfo_icon()
        fetchResult = params%requestList%fetch3d(varName, tileinfo_icon%idx, jg, field)
        IF (PRESENT(found)) THEN
          found = fetchResult
        ELSE
          CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
        ENDIF
    END IF
  END SUBROUTINE fetch3d

  ! Wrapper around fetch3D which allows to have a correct summary table entry
  !  for optional input fields. It continues in case of missing field in input
  !  data, sets correct "FG read attempt" etc. entries in the summary table
  !  and returns a logical flag if success or not:
  SUBROUTINE fetch3d_with_status (caller, context, params, varName, jg, field3d, lfound)
      
    CHARACTER (len=*), INTENT(in)      :: caller   ! string to identify the caller in the diagn. message below
    CHARACTER (len=*), INTENT(in)      :: context  ! string to identify the calling context in the diagn. message below
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER (len=*), INTENT(in)      :: varName  ! name of vari
    INTEGER, INTENT(in)                :: jg
    REAL(wp), INTENT(inout)            :: field3d(:,:,:)
    LOGICAL, INTENT (out)              :: lfound

    CALL fetch3d(params, TRIM(varName), jg, field3d, lfound)
    ! Error handler and summary table entry
    CALL params%inputInstructions(jg)%ptr%handleError(lfound, TRIM(varName), params%routine, params%isFg)
    IF (.NOT. lfound) THEN
      CALL message(caller, 'Variable '//TRIM(varName)//' not found in '//TRIM(context))
    END IF
    
  END SUBROUTINE fetch3d_with_status

  SUBROUTINE fetchSurface(params, varName, jg, field, found)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, INTENT(in) :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:)
    LOGICAL, INTENT(OUT), OPTIONAL :: found ! allows to do the error handling in the calling routine

    LOGICAL :: fetchResult
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        tileinfo_icon = trivial_tile_att%getTileinfo_icon()
        fetchResult = params%requestList%fetchSurface(varName, tileinfo_icon%idx, jg, field)
        IF (PRESENT(found)) THEN
          found = fetchResult
        ELSE
          CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
        ENDIF
    END IF
  END SUBROUTINE fetchSurface

  ! Wrapper for requestList%fetchTiledSurface() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  SUBROUTINE fetchTiledSurface(params, varName, jg, tileCount, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, INTENT(in) :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    INTEGER :: jt
    LOGICAL :: fetchResult
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            fetchResult = .TRUE.
            DO jt = 1, tileCount
                tileinfo_icon = trivial_tile_att%getTileinfo_icon()
                fetchResult = fetchResult.AND.params%requestList%fetchSurface(varName, tileinfo_icon%idx, jg, field(:,:,jt))
            END DO
        ELSE
            !True tiled input.
            fetchResult = params%requestList%fetchTiledSurface(varName, jg, field)
        END IF
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetchTiledSurface

  ! Wrapper for requestList%fetchTiled3d() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  SUBROUTINE fetchTiled3d(params, varName, jg, tileCount, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, INTENT(in) :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:,:)

    INTEGER :: jt
    LOGICAL :: fetchResult
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            fetchResult = .TRUE.
            DO jt = 1, tileCount
                tileinfo_icon = trivial_tile_att%getTileinfo_icon()
                fetchResult = fetchResult.AND.params%requestList%fetch3d(varName, tileinfo_icon%idx, jg, field(:,:,:,jt))
            END DO
        ELSE
            !True tiled input.
            fetchResult = params%requestList%fetchTiled3d(varName, jg, field)
        END IF
        CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
    END IF
  END SUBROUTINE fetchTiled3d


  ! Wrapper for requestList%fetchTiled3d() 
  ! - that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  ! - has a fallback option, if the primary input field is not found.
  SUBROUTINE fetchTiled3dWithFallback(params, varName, varNameFallback, jg, tileCount, field, opt_field_fallback)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName, varNameFallback
    INTEGER, INTENT(in) :: jg, tileCount
    REAL(wp),           TARGET, INTENT(INOUT) :: field(:,:,:,:)
    REAL(wp), OPTIONAL, TARGET, INTENT(INOUT) :: opt_field_fallback(:,:,:,:)  ! optional target field for 
                                                                              ! fallback input
    ! local
    INTEGER :: jt
    LOGICAL :: fetchResult
    REAL(wp), POINTER :: ptr_field_fb(:,:,:,:)
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    ! set pointer to fallback target field
    IF (PRESENT(opt_field_fallback)) THEN
      ptr_field_fb => opt_field_fallback
    ELSE
      ptr_field_fb => field
    ENDIF

    tileinfo_icon = trivial_tile_att%getTileinfo_icon()
    IF(params%inputInstructions(jg)%ptr%wantVar(varName, params%isFg)) THEN
        IF(ltile_coldstart) THEN
            !Fake tiled input by copying the input field to all tiles.
            fetchResult = .TRUE.
            DO jt = 1, tileCount
                fetchResult = fetchResult.AND.  &
                  &           params%requestList%fetch3d(varName, tileinfo_icon%idx, jg, field(:,:,:,jt))
            END DO
            IF (.NOT.fetchResult) THEN
                CALL params%inputInstructions(jg)%ptr%optionalReadResult(fetchResult, varName, params%routine, params%isFg)
                ! try fallback
                fetchResult = .TRUE.
                DO jt = 1, tileCount
                   fetchResult = fetchResult.AND. &
                     &           params%requestList%fetch3d(varNameFallback, tileinfo_icon%idx, jg, ptr_field_fb(:,:,:,jt))
                END DO
                CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varNameFallback, params%routine, params%isFg)
            ELSE
                CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
            ENDIF
        ELSE
            !True tiled input.
            fetchResult = params%requestList%fetchTiled3d(varName, jg, field)
            IF (.NOT.fetchResult) THEN
                CALL params%inputInstructions(jg)%ptr%optionalReadResult(fetchResult, varName, params%routine, params%isFg)
                ! try fallback
                fetchResult = params%requestList%fetchTiled3d(varNameFallback, jg, ptr_field_fb)
                CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varNameFallback, params%routine, params%isFg)
            ELSE
                CALL params%inputInstructions(jg)%ptr%handleError(fetchResult, varName, params%routine, params%isFg)
            ENDIF
        END IF
    END IF
  END SUBROUTINE fetchTiled3dWithFallback


  SUBROUTINE fetchRequired3d(params, varName, jg, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, INTENT(in) :: jg
    REAL(wp), INTENT(INOUT) :: field(:,:,:)
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    tileinfo_icon = trivial_tile_att%getTileinfo_icon()
    CALL params%requestList%fetchRequired3d(varName, tileinfo_icon%idx, jg, field)
    CALL params%inputInstructions(jg)%ptr%handleError(.TRUE., varName, params%routine, params%isFg)
  END SUBROUTINE fetchRequired3d

  ! Wrapper for requestList%fetchRequiredTiledSurface() that falls back to reading copies of untiled input IF ltile_coldstart IS set.
  SUBROUTINE fetchRequiredTiledSurface(params, varName, jg, tileCount, field)
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    CHARACTER(LEN = *), INTENT(IN) :: varName
    INTEGER, INTENT(IN) :: jg, tileCount
    REAL(wp), INTENT(INOUT) :: field(:,:,:)

    INTEGER :: jt
    TYPE(t_tileinfo_icon) :: tileinfo_icon

    IF(ltile_coldstart) THEN
        !Fake tiled input by copying the input field to all tiles.
        tileinfo_icon = trivial_tile_att%getTileinfo_icon()
        DO jt = 1, tileCount
            CALL params%requestList%fetchRequiredSurface(varName, tileinfo_icon%idx, jg, field(:,:,jt))
        END DO
    ELSE
        !True tiled input.
        CALL params%requestList%fetchRequiredTiledSurface(varName, jg, field)
    END IF
    CALL params%inputInstructions(jg)%ptr%handleError(.TRUE., varName, params%routine, params%isFg)
  END SUBROUTINE fetchRequiredTiledSurface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  !! Fetch the DWD first guess from the request list (atmosphere only)
  !! First guess (FG) is read for theta_v, rho, vn, w, tke,
  !! whereas DA output is read for T, p, u, v,
  !! qv, qc, qi, qr, qs, qg.
  SUBROUTINE fetch_dwdfg_atm(requestList, p_patch, p_nh_state, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(INOUT) :: p_patch(:)
    TYPE(t_nh_state), INTENT(INOUT), TARGET :: p_nh_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(*), PARAMETER :: routine = modname//':fetch_dwdfg_atm'
    INTEGER                  :: jg, nlev_in, jb,jk,jc,nlen
    TYPE(t_nh_prog), POINTER :: prognosticFields
    REAL(wp), POINTER        :: my_ptr3d(:,:,:)
    TYPE(t_fetchParams)      :: params
    REAL(wp), ALLOCATABLE    :: z_ifc_in(:,:,:)

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            ! save some paperwork
            prognosticFields => p_nh_state(jg)%prog(nnow(jg))

            ! request the first guess fields (atmosphere only)
            CALL fetchRequired3d(params, 'theta_v', jg, prognosticFields%theta_v)
            CALL fetchRequired3d(params, 'rho', jg, prognosticFields%rho)
            CALL fetchRequired3d(params, 'w', jg, prognosticFields%w)
            CALL fetchRequired3d(params, 'tke', jg, prognosticFields%tke)

            ! Only needed for FG-only runs; usually read from ANA
            my_ptr3d => prognosticFields%tracer(:,:,:,iqv)
            CALL fetch3d(params, 'qv', jg, my_ptr3d)

            my_ptr3d => prognosticFields%tracer(:,:,:,iqc)
            CALL fetch3d_with_status(routine, 'dwdfg file', params, 'qc', jg, my_ptr3d, &
                 atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqc))

            my_ptr3d => prognosticFields%tracer(:,:,:,iqi)
            CALL fetch3d_with_status(routine, 'dwdfg file', params, 'qi', jg, my_ptr3d, &
                 atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqi))

            IF ( iqr /= 0 ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqr)
              CALL fetch3d_with_status(routine, 'dwdfg file', params, 'qr', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqr))
            END IF

            IF ( iqs /= 0 ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqs)
              CALL fetch3d_with_status(routine, 'dwdfg file', params, 'qs', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqs))
            END IF

            IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqg)
              CALL fetch3d_with_status(routine, 'dwdfg file', params, 'qg', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqg))
            END IF

            IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
              my_ptr3d => prognosticFields%tracer(:,:,:,iqh)
              CALL fetch3d_with_status(routine, 'dwdfg file', params, 'qh', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqh))
              my_ptr3d => prognosticFields%tracer(:,:,:,iqnc)
              CALL fetch3d_with_status (routine, 'dwdfg file', params, 'qnc', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqnc))
              my_ptr3d => prognosticFields%tracer(:,:,:,iqni)
              CALL fetch3d_with_status (routine, 'dwdfg file', params, 'qni', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqni))
              my_ptr3d => prognosticFields%tracer(:,:,:,iqnr)
              CALL fetch3d_with_status (routine, 'dwdfg file', params, 'qnr', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqnr))
              my_ptr3d => prognosticFields%tracer(:,:,:,iqns)
              CALL fetch3d_with_status (routine, 'dwdfg file', params, 'qns', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqns))
              my_ptr3d => prognosticFields%tracer(:,:,:,iqng)
              CALL fetch3d_with_status (routine, 'dwdfg file', params, 'qng', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqng))
              my_ptr3d => prognosticFields%tracer(:,:,:,iqnh)
              CALL fetch3d_with_status (routine, 'dwdfg file', params, 'qnh', jg, my_ptr3d, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_fg(iqnh))
              ! If any of the number conc. qnx could not be read from fg, this will diagnose them from qx:
              CALL init_qnx_from_qx_twomom (routine, p_patch(jg), prognosticFields, &
                                            .NOT. atm_phy_nwp_config(jg)%lhydrom_read_from_fg(:) )
            END IF

            IF (lvert_remap_fg) THEN

                ! the number of input and output levels must be the same for this mode
                CALL allocate_extana_atm(nblks_c  = p_patch(jg)%nblks_c, & 
                  &                      nblks_e  = p_patch(jg)%nblks_e, & 
                  &                      nlev_in  = p_patch(jg)%nlev,    & 
                  &                      atm_in   = initicon(jg)%atm_in, & 
                  &                      const    = initicon(jg)%const  ) !inout

                ! allocate temporary array:
                nlev_in = SIZE(initicon(jg)%const%z_mc_in,2)
                ALLOCATE(z_ifc_in(nproma,nlev_in+1, p_patch(jg)%nblks_c))

                CALL fetchRequired3d(params, 'z_ifc', jg, z_ifc_in)
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,nlen) ICON_OMP_DEFAULT_SCHEDULE
                DO jb = 1, p_patch(jg)%nblks_c
                  
                  IF (jb /= p_patch(jg)%nblks_c) THEN
                    nlen = nproma
                  ELSE
                    nlen = p_patch(jg)%npromz_c
                  ENDIF

                  ! interpolate half-level variables to full levels and diagnose pressure and temperature 
                  DO jk = 1, nlev_in
                    DO jc = 1, nlen
                      
                      initicon(jg)%const%z_mc_in(jc,jk,jb) = (z_ifc_in(jc,jk,jb) + &
                        &   z_ifc_in(jc,jk+1,jb)) * 0.5_wp
                      
                    END DO
                  END DO
                END DO
!$OMP END DO
!$OMP END PARALLEL

                ! clean-up
                DEALLOCATE(z_ifc_in)

                ! fetch additional tracers in first guess
                IF ( init_mode /= MODE_IAU_OLD ) &
                  &  CALL fetch_tracer_fg('tracer_fg_in', params, jg,     &
                  &                        atm_in  = initicon(jg)%atm_in, &
                  &                        nblks_c = p_patch(jg)%nblks_c, &
                  &                        nlev_in = p_patch(jg)%nlev)

            ELSE

                ! fetch additional tracers in first guess
                IF ( init_mode /= MODE_IAU_OLD ) &
                  &  CALL fetch_tracer_fg('tracer_fg_in', params, jg,        &
                  &                        tracer = prognosticFields%tracer)

            END IF

            CALL fetchRequired3d(params, 'vn', jg, prognosticFields%vn)
        END IF
    END DO

  END SUBROUTINE fetch_dwdfg_atm

  !>
  !! Fetch DWD first guess from the request list (atmosphere only) and store to initicon input state
  !! for subsequent vertical remapping to the current ICON grid
  !! First guess (FG) is read for z_ifc, theta_v, rho, vn, w, tke,
  !! whereas DA output is read for T, p, u, v,
  !! qv, qc, qi, qr, qs, qg.
  !!
  SUBROUTINE fetch_dwdfg_atm_ii(requestList, p_patch, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_atm_ii'
    INTEGER :: jg, jb, jk, jc, i_endidx, nlev_in
    REAL(wp) :: tempv, exner
    REAL(dp), POINTER :: levelValues(:)
    TYPE(t_fetchParams) :: params
    REAL(wp), ALLOCATABLE :: z_ifc_in(:,:,:), w_ifc(:,:,:), tke_ifc(:,:,:)
    LOGICAL :: lfound_thv, lfound_rho, lfound_vn, lfound_hus, lfound_clw, lfound_cli, lfound_qr, lfound_qs, lfound_qg

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN  ! Skip reading the atmospheric input data if a model domain is not active at initial time

            ! determine number of HALF LEVELS of generalized Z-AXIS
            levelValues => requestList%getLevels('z_ifc', jg)
            IF(.NOT.ASSOCIATED(levelValues)) CALL finish(routine, "no DATA found for domain "//TRIM(int2string(jg))//" of &
                                                                  &required variable 'z_ifc'")

            CALL allocate_extana_atm(nblks_c  = p_patch(jg)%nblks_c,    & 
              &                      nblks_e  = p_patch(jg)%nblks_e,    & 
              &                      nlev_in  = SIZE(levelValues,1)-1,  & 
              &                      atm_in   = initicon(jg)%atm_in  ,  & 
              &                      const    = initicon(jg)%const    ) !inout

            ! allocate temporary array:
            nlev_in = SIZE(initicon(jg)%const%z_mc_in,2)
            ALLOCATE(z_ifc_in(nproma,nlev_in+1, p_patch(jg)%nblks_c), &
              &      w_ifc   (nproma,nlev_in+1, p_patch(jg)%nblks_c), &
              &      tke_ifc (nproma,nlev_in+1, p_patch(jg)%nblks_c) )

            ! start reading first guess (atmosphere only)
            CALL fetchRequired3d(params, 'z_ifc', jg, z_ifc_in)
            CALL fetch3d(params, 'theta_v', jg, initicon(jg)%atm_in%theta_v, lfound_thv)
            CALL fetch3d(params, 'rho', jg, initicon(jg)%atm_in%rho, lfound_rho)
            IF (.NOT. (lfound_thv .AND. lfound_rho) ) THEN ! try fetching the diagnostic thermodynamic variables
              CALL fetchRequired3d(params, 'temp', jg, initicon(jg)%atm_in%temp)
              CALL fetchRequired3d(params, 'pres', jg, initicon(jg)%atm_in%pres)
            ENDIF
            CALL fetchRequired3d(params, 'w',   jg, w_ifc)
            ! If the TKE field is not in the input data, a cold-start initialization is executed in init_nwp_phy
            CALL fetch3d(params, 'tke', jg, tke_ifc, lread_tke)

            ! Check for both DWD/Sapphire variables due to naming mismatch
            CALL fetch3d(params, 'hus', jg, initicon(jg)%atm_in%qv, lfound_hus)
            CALL fetch3d(params, 'clw', jg, initicon(jg)%atm_in%qc, lfound_clw)
            CALL fetch3d(params, 'cli', jg, initicon(jg)%atm_in%qi, lfound_cli)
!$OMP PARALLEL
            IF (.NOT. lfound_hus) CALL fetchRequired3d(params, 'qv', jg, initicon(jg)%atm_in%qv)
            IF (.NOT. lfound_clw) CALL fetchRequired3d(params, 'qc', jg, initicon(jg)%atm_in%qc)
            IF (.NOT. lfound_cli) CALL fetchRequired3d(params, 'qi', jg, initicon(jg)%atm_in%qi)
!$OMP END PARALLEL

            CALL fetch3d(params, 'qr', jg, initicon(jg)%atm_in%qr, lfound_qr)
            CALL fetch3d(params, 'qs', jg, initicon(jg)%atm_in%qs, lfound_qs)
!$OMP PARALLEL
            IF (.NOT. lfound_qr) CALL init(initicon(jg)%atm_in%qr(:,:,:), lacc=.FALSE.)
            IF (.NOT. lfound_qs) CALL init(initicon(jg)%atm_in%qs(:,:,:), lacc=.FALSE.)
!$OMP END PARALLEL
            
            ! Set up flags for later calls
            IF (lfound_hus) latbc_config%latbc_contains_hus = .TRUE.
            
            ! fetch additional tracers in first guess
            CALL fetch_tracer_fg('tracer_fg_in', params, jg,       &
              &                   atm_in  = initicon(jg)%atm_in,   &
              &                   nblks_c = p_patch(jg)%nblks_c,   &
              &                   nlev_in = SIZE(levelValues,1)-1)

            IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
              CALL fetch3d(params, 'qg', jg, initicon(jg)%atm_in%qg, lfound_qg)
              IF (.NOT. lfound_qg) THEN
                CALL message(routine,'Graupel (QG) not available in input data')
                initicon(jg)%atm_in%qg(:,:,:) = 0._wp
              END IF
            ELSE
              initicon(jg)%atm_in%qg(:,:,:) = 0._wp
            END IF

            CALL fetch3d(params, 'vn', jg, initicon(jg)%atm_in%vn, lfound_vn)
            IF (lfound_vn) THEN
              lread_vn = .TRUE.  ! Tell the vertical interpolation routine that vn needs to be processed
            ELSE ! try fetching U and V components
              CALL fetchRequired3d(params, 'u', jg, initicon(jg)%atm_in%u)
              CALL fetchRequired3d(params, 'v', jg, initicon(jg)%atm_in%v)
              lread_vn = .FALSE.
            ENDIF

            ! Interpolate half level variables from interface levels to main levels, and convert thermodynamic variables
            ! into temperature and pressure as expected by the vertical interpolation routine  
!$OMP PARALLEL
!$OMP DO PRIVATE (jk,jc,jb,i_endidx,tempv,exner) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = 1,p_patch(jg)%nblks_c

                i_endidx = MERGE(nproma, p_patch(jg)%npromz_c, jb /= p_patch(jg)%nblks_c)

                DO jk = 1, initicon(jg)%atm_in%nlev

                  DO jc = 1, i_endidx
                    initicon(jg)%const%z_mc_in(jc,jk,jb) = (z_ifc_in(jc,jk,jb) + z_ifc_in(jc,jk+1,jb)) * 0.5_wp
                    initicon(jg)%atm_in%w(jc,jk,jb) = (w_ifc(jc,jk,jb) + w_ifc(jc,jk+1,jb)) * 0.5_wp
                  END DO

                  IF (lread_tke) THEN
                    DO jc = 1, i_endidx
                      initicon(jg)%atm_in%tke(jc,jk,jb) = (tke_ifc(jc,jk,jb) + tke_ifc(jc,jk+1,jb)) * 0.5_wp
                    END DO
                  ELSE
                    initicon(jg)%atm_in%tke(:,:,jb) = 0._wp
                  ENDIF

                  IF (lfound_thv .AND. lfound_rho) THEN
                    DO jc = 1, i_endidx

                      exner = (initicon(jg)%atm_in%rho(jc,jk,jb)*initicon(jg)%atm_in%theta_v(jc,jk,jb)*rd/p0ref)**(1._wp/cvd_o_rd)
                      tempv = initicon(jg)%atm_in%theta_v(jc,jk,jb)*exner

                      initicon(jg)%atm_in%pres(jc,jk,jb) = exner**(cpd/rd)*p0ref
                      initicon(jg)%atm_in%temp(jc,jk,jb) = tempv / (1._wp + vtmpc1*initicon(jg)%atm_in%qv(jc,jk,jb) - &
                        (initicon(jg)%atm_in%qc(jc,jk,jb) + initicon(jg)%atm_in%qi(jc,jk,jb) +                        &
                         initicon(jg)%atm_in%qr(jc,jk,jb) + initicon(jg)%atm_in%qs(jc,jk,jb) +                        &
                         initicon(jg)%atm_in%qg(jc,jk,jb)                                      ) )

                    END DO
                  ENDIF
                END DO
              END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

              ! cleanup
              IF (ALLOCATED(z_ifc_in)) DEALLOCATE(z_ifc_in)
              IF (ALLOCATED(w_ifc))    DEALLOCATE(w_ifc)
              IF (ALLOCATED(tke_ifc))  DEALLOCATE(tke_ifc)

        END IF !ldom_active
    ENDDO

  END SUBROUTINE fetch_dwdfg_atm_ii


  !>
  !! Fetch DA-analysis DATA from the request list (atmosphere only)
  !!
  !! Depending on the initialization mode, either full fields or increments
  !! are read (atmosphere only). The following full fields are read, if available:
  !!     u, v, t, p, qv, qi, qc, qr, qs, qg
  !!
  SUBROUTINE fetch_dwdana_atm(requestList, p_patch, p_nh_state, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch), INTENT(IN) :: p_patch(:)
    TYPE(t_nh_state), INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_atm'
    INTEGER :: jg
    TYPE(t_pi_atm), POINTER :: my_ptr
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN  ! Skip reading the atmospheric input data if a model domain is not active at initial time
            ! Depending on the initialization mode chosen (incremental vs. non-incremental)
            ! input fields are stored in different locations.
            IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
                ! Skip this domain if its interpolated from its parent in process_input_dwdana_atm()
                IF (lp2cintp_incr(jg)) CYCLE
                my_ptr => initicon(jg)%atm_inc
            ELSE
                my_ptr => initicon(jg)%atm
            ENDIF

            ! start reading DA output (atmosphere only)
            ! The dynamical variables temp, pres, u and v, which need further processing,
            ! are either stored in initicon(jg)%atm or initicon(jg)%atm_inc, depending on whether
            ! IAU is used or not. The moisture variables, which can be taken over directly from
            ! the Analysis, are written to the NH prognostic state
            CALL fetch3d(params, 'temp', jg, my_ptr%temp)
            CALL fetch3d(params, 'pres', jg, my_ptr%pres)
            CALL fetch3d(params, 'u', jg, my_ptr%u)
            CALL fetch3d(params, 'v', jg, my_ptr%v)

            IF ( ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qv') == kInputSourceFg
                CALL fetch3d(params, 'qv', jg, my_ptr%qv)
                ! check whether we are using DATA from both FG and ANA input, so that it's correctly listed in the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qv') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('qv', kInputSourceBoth)
                END IF
            ELSE
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqv)
                CALL fetch3d(params, 'qv', jg, my_ptr3d)
            ENDIF

            IF (init_mode == MODE_IAU .AND. qcana_mode > 0) THEN
              lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qc') == kInputSourceFg
              CALL fetch3d_with_status (routine, 'dwdana file', params, 'qc', jg, my_ptr%qc, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqc))
              IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qc') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qc', kInputSourceBoth)
              END IF
            ENDIF

            IF (init_mode == MODE_IAU .AND. qiana_mode > 0) THEN
              lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qi') == kInputSourceFg
              CALL fetch3d_with_status (routine, 'dwdana file', params, 'qi', jg, my_ptr%qi, &
                   atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqi))
              IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qi') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qi', kInputSourceBoth)
              END IF
            ENDIF

            IF (init_mode == MODE_IAU .AND. qrsgana_mode > 0) THEN
              IF ( iqr /= 0 ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qr') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qr', jg, my_ptr%qr, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqr))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qr') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qr', kInputSourceBoth)
                END IF
              END IF
              IF ( iqs /= 0 ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qs') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qs', jg, my_ptr%qs, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqs))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qs') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qs', kInputSourceBoth)
                END IF
              END IF
              IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qg') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qg', jg, my_ptr%qg, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqg))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qg') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qg', kInputSourceBoth)
                END IF
              END IF

              IF ( atm_phy_nwp_config(jg)%l2moment ) THEN

                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qh') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qh', jg, my_ptr%qh, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqh))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qh') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qh', kInputSourceBoth)
                END IF

              END IF
            END IF

            IF (init_mode == MODE_IAU .AND. atm_phy_nwp_config(jg)%l2moment) THEN

              IF (qnxana_2mom_mode > 0) THEN

                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qnc') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qnc', jg, my_ptr%qnc, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqnc))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qnc') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qnc', kInputSourceBoth)
                END IF
                
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qni') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qni', jg, my_ptr%qni, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqni))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qni') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qni', kInputSourceBoth)
                END IF
                
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qnr') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qnr', jg, my_ptr%qnr, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqnr))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qnr') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qnr', kInputSourceBoth)
                END IF
                
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qns') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qns', jg, my_ptr%qns, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqns))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qns') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qns', kInputSourceBoth)
                END IF
                
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qng') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qng', jg, my_ptr%qng, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqng))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qng') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qng', kInputSourceBoth)
                END IF
                
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('qnh') == kInputSourceFg
                CALL fetch3d_with_status (routine, 'dwdana file', params, 'qnh', jg, my_ptr%qnh, &
                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(iqnh))
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('qnh') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('qnh', kInputSourceBoth)
                END IF

                ! QNX increments which were not available from dwdana file are diagnosed based on the QX from fg and ana increments:
                CALL init_qnxinc_from_qxinc_twomom ( routine, p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                                     initicon(jg), &
                                                     atm_phy_nwp_config(jg)%lhydrom_read_from_fg(:), &
                                                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(:), &
                                                     .NOT. atm_phy_nwp_config(jg)%lhydrom_read_from_ana(:) )
              ELSE

                ! QNX increments not read from dwdana file, but diagnosed based on the QX from fg and ana increments:
                CALL init_qnxinc_from_qxinc_twomom ( routine, p_patch(jg), p_nh_state(jg)%prog(nnow(jg)), &
                                                     initicon(jg), &
                                                     atm_phy_nwp_config(jg)%lhydrom_read_from_fg(:), &
                                                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(:), &
                                                     atm_phy_nwp_config(jg)%lhydrom_read_from_ana(:) .OR. .TRUE. )

              END IF

            END IF

            ! In initialization modes other than IAU, analyses are usually read only in the special operation
            ! mode for interpolating uninitialized analyses. These currently do not contain cloud and precip variables,
            ! but processing them would be enabled here
            IF ( .NOT. ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode) ) THEN
              my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqc)
              CALL fetch3d(params, 'qc', jg, my_ptr3d)
              my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqi)
              CALL fetch3d(params, 'qi', jg, my_ptr3d)
              IF ( iqr /= 0 ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqr)
                CALL fetch3d(params, 'qr', jg, my_ptr3d)
              END IF
              IF ( iqs /= 0 ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqs)
                CALL fetch3d(params, 'qs', jg, my_ptr3d)
              END IF
              IF ( atm_phy_nwp_config(jg)%lhave_graupel ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqg)
                CALL fetch3d(params, 'qg', jg, my_ptr3d)
              END IF

              IF ( atm_phy_nwp_config(jg)%l2moment ) THEN
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqh)
                CALL fetch3d(params, 'qh', jg, my_ptr3d)
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqnc)
                CALL fetch3d(params, 'qnc', jg, my_ptr3d)
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqni)
                CALL fetch3d(params, 'qni', jg, my_ptr3d)
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqnr)
                CALL fetch3d(params, 'qnr', jg, my_ptr3d)
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqns)
                CALL fetch3d(params, 'qns', jg, my_ptr3d)
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqng)
                CALL fetch3d(params, 'qng', jg, my_ptr3d)
                my_ptr3d => p_nh_state(jg)%prog(nnow(jg))%tracer(:,:,:,iqnh)
                CALL fetch3d(params, 'qnh', jg, my_ptr3d)
              END IF
            ENDIF

        END IF
    ENDDO ! loop over model domains

  END SUBROUTINE fetch_dwdana_atm

  !>
  !! Fetch additional tracers in the first guess from the request list
  !! Lopp over all ICON vars and collect info similar to collect_group
  SUBROUTINE fetch_tracer_fg(grp_name, params, jg, tracer, atm_in, nblks_c, nlev_in)
    CHARACTER(LEN=*),    INTENT(IN)    :: grp_name
    TYPE(t_fetchParams), INTENT(INOUT) :: params
    INTEGER,             INTENT(IN)    :: jg
    REAL(wp), POINTER,   INTENT(INOUT), OPTIONAL :: tracer(:,:,:,:)
    TYPE(t_pi_atm_in),   INTENT(INOUT), OPTIONAL :: atm_in
    INTEGER,             INTENT(IN), OPTIONAL    :: nlev_in, nblks_c
    INTEGER                           :: i, grp_id, idx
    TYPE(t_var_metadata), POINTER     :: info
    CHARACTER(:), ALLOCATABLE :: vname
    REAL(wp), POINTER                 :: my_ptr3d(:,:,:)
    TYPE(t_vl_register_iter) :: vl_iter

    idx = 0
    grp_id = var_groups_dyn%group_id(grp_name)

    ! loop over all variable lists and variables
    DO WHILE(vl_iter%next())
      IF (vl_iter%cur%p%vlevel_type /= level_type_ml) CYCLE
      ! do not inspect variable list if its domain does not match:
      IF (vl_iter%cur%p%patch_id /= jg)  CYCLE
      LOOPVAR : DO i = 1, vl_iter%cur%p%nvars
        info => vl_iter%cur%p%vl(i)%p%info
        ! Do not inspect element if it is a container
        IF (.NOT. info%lcontainer .AND. info%in_group(grp_id)) THEN
          idx = idx + 1
          vname = get_var_name(info)
          IF ( PRESENT(atm_in) ) THEN
            atm_in%tracer(idx)%var_element => vl_iter%cur%p%vl(i)%p
            ! allocate source array for vertical interpolation
            ALLOCATE(atm_in%tracer(idx)%field(nproma,nlev_in,nblks_c))
!$OMP PARALLEL
            CALL init(atm_in%tracer(idx)%field(:,:,:), lacc=.FALSE.)      !_jf: necessary?
!$OMP END PARALLEL
            ! request the first guess fields
            my_ptr3d => atm_in%tracer(idx)%field(:,:,:)
            CALL fetch3d(params, vname, jg, my_ptr3d)
          ENDIF
          IF ( PRESENT(tracer) ) THEN
            ! request the first guess fields
            my_ptr3d => tracer(:,:,:,info%ncontained)
            CALL fetch3d(params, vname, jg, my_ptr3d)
          END IF
        END IF
      ENDDO LOOPVAR ! loop over vlist "i"
    ENDDO ! i = 1, nvar_lists
  END SUBROUTINE fetch_tracer_fg

  !XXX: Cannot be moved into fetch_dwdana_atm() since it needs input from fetch_dwdana_sfc()
  SUBROUTINE process_input_dwdana_atm (p_patch, initicon)
    TYPE(t_patch),          INTENT(IN)    :: p_patch(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)

    INTEGER :: jg

    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active) THEN
        IF (ANY((/MODE_IAU,MODE_IAU_OLD/) == init_mode)) THEN
          IF (lp2cintp_incr(jg)) THEN
            ! Perform parent-to-child interpolation of atmospheric DA increments
            CALL interpolate_scal_increments(initicon, p_patch(jg)%parent_id, jg)
          END IF
        END IF
      END IF
    END DO
  END SUBROUTINE process_input_dwdana_atm


  !>
  !! Fetch DWD first guess DATA from request list (land/surface only)
  SUBROUTINE fetch_dwdfg_sfc(requestList, p_patch, prm_diag, prm_nwp_stochconv, p_nh_state, p_lnd_state, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_nwp_phy_diag),      INTENT(INOUT) :: prm_diag(:)
    TYPE(t_nwp_phy_stochconv), INTENT(INOUT) :: prm_nwp_stochconv(:)
    TYPE(t_nh_state),  TARGET, INTENT(INOUT) :: p_nh_state(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    INTEGER :: jg, error
    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(dp), POINTER :: levels(:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_nh_diag),  POINTER :: nh_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    TYPE(t_fetchParams) :: params

    CHARACTER(len=*), PARAMETER :: routine = modname//':fetch_dwdfg_sfc'

    ! Workaround for the intel compiler botching implicit allocation.
    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")

    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN ! Skip reading the atmospheric input data if a model domain is not active at initial time
            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            nh_diag  => p_nh_state(jg)%diag
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

            ! COSMO-DE does not provide sea ice field. In that case set fr_seaice to 0
            IF (init_mode /= MODE_COSMO) THEN
                CALL fetchSurface(params, 'fr_seaice', jg, lnd_diag%fr_seaice)
            ELSE
!$OMP PARALLEL
                CALL init(lnd_diag%fr_seaice(:,:), lacc=.FALSE.)
!$OMP END PARALLEL
            ENDIF ! init_mode /= MODE_COSMO

            ! sea-ice related fields
            CALL fetchSurface(params, 't_ice', jg, wtr_prog%t_ice)
            CALL fetchSurface(params, 'h_ice', jg, wtr_prog%h_ice)
            IF ( lprog_albsi ) THEN  ! prognostic sea-ice albedo
              CALL fetchSurface(params, 'alb_si', jg, wtr_prog%alb_si)
            ENDIF

            !These two fields are required for the processing step below, AND they are NOT initialized before this SUBROUTINE IS called, so they are fetched as required.
            !This diverges from the code that I found which READ them conditionally.
            CALL fetchRequiredTiledSurface(params, 't_g', jg, ntiles_total + ntiles_water, lnd_prog%t_g_t)
            ! Limit SST in case of interpolated data to physically reasonable values
            IF (ltile_coldstart) lnd_prog%t_g_t(:,:,isub_water) = MIN(305._wp,lnd_prog%t_g_t(:,:,isub_water))
            CALL fetchRequiredTiledSurface(params, 'qv_s', jg, ntiles_total + ntiles_water, lnd_diag%qv_s_t)

            CALL fetchTiledSurface(params, 'freshsnow', jg, ntiles_total, lnd_diag%freshsnow_t)
            IF (lsnowtile .AND. .NOT. ltile_coldstart) THEN
                CALL fetchTiledSurface(params, 'snowfrac_lc', jg, ntiles_total, lnd_diag%snowfrac_lc_t)
            ENDIF
            CALL fetchTiledSurface(params, 'w_snow', jg, ntiles_total, lnd_prog%w_snow_t)
            CALL fetchTiledSurface(params, 'w_i', jg, ntiles_total, lnd_prog%w_i_t)
            CALL fetchTiledSurface(params, 'h_snow', jg, ntiles_total, lnd_diag%h_snow_t)
            CALL fetchTiledSurface(params, 't_snow', jg, ntiles_total, lnd_prog%t_snow_t)
            CALL fetchTiledSurface(params, 'rho_snow', jg, ntiles_total, lnd_prog%rho_snow_t)

            ! t_sk is always fetched if present in the input data because using itype_canopy=1 with input data containing
            ! t_sk causes an error otherwise
            CALL fetchTiledSurface(params, 't_sk', jg, ntiles_total + ntiles_water, lnd_prog%t_sk_t)

            IF (itype_trvg == 3) THEN
              CALL fetchTiledSurface(params, 'plantevap', jg, ntiles_total, lnd_diag%plantevap_t)
            ENDIF
            IF (icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) THEN
              CALL fetchSurface(params, 't2m_bias', jg, nh_diag%t2m_bias)
            ENDIF
            IF (icpl_da_sfcevap >= 2) THEN
              CALL fetchSurface(params, 'rh_avginc', jg, nh_diag%rh_avginc)
            ENDIF
            IF (icpl_da_sfcevap >= 3 .OR. icpl_da_snowalb >= 1) THEN
              CALL fetchSurface(params, 't_avginc', jg, nh_diag%t_avginc)
            ENDIF
            IF (icpl_da_sfcfric >= 1) THEN
              CALL fetchSurface(params, 'vabs_avginc', jg, nh_diag%vabs_avginc)
            ENDIF
            IF (icpl_da_skinc >= 1) THEN
              CALL fetchSurface(params, 't_wgt_avginc', jg, nh_diag%t_wgt_avginc)
            ENDIF
            IF (icpl_da_sfcevap >= 5) THEN
              CALL fetchSurface(params, 't_daywgt_avginc', jg, nh_diag%t_daywgt_avginc)
              CALL fetchSurface(params, 'rh_daywgt_avginc', jg, nh_diag%rh_daywgt_avginc)
            ENDIF
            IF (latbc_config%fac_latbc_presbiascor > 0._wp) THEN
              CALL fetchSurface(params, 'p_avginc', jg, nh_diag%p_avginc)
            ENDIF

            IF (itype_snowevap == 3) THEN
              CALL fetchSurface(params, 'hsnow_max', jg, lnd_diag%hsnow_max)
              CALL fetchSurface(params, 'snow_age',  jg, lnd_diag%snow_age)
            ENDIF

            IF (lmulti_snow) THEN
                ! multi layer snow fields
                CALL fetchTiled3d(params, 't_snow_mult', jg, ntiles_total, lnd_prog%t_snow_mult_t)
                CALL fetchTiled3d(params, 'rho_snow_mult', jg, ntiles_total, lnd_prog%rho_snow_mult_t)
                CALL fetchTiled3d(params, 'wtot_snow', jg, ntiles_total, lnd_prog%wtot_snow_t)
                CALL fetchTiled3d(params, 'wliq_snow', jg, ntiles_total, lnd_prog%wliq_snow_t)
                CALL fetchTiled3d(params, 'dzh_snow', jg, ntiles_total, lnd_prog%dzh_snow_t)
            ELSE IF (l2lay_rho_snow) THEN
                CALL fetchTiled3d(params, 'rho_snow_mult', jg, ntiles_total, lnd_prog%rho_snow_mult_t)
                IF (inputInstructions(jg)%ptr%sourceOfVar('rho_snow_mult') == kInputSourceCold) THEN
                    ! initialize top-layer snow density with average density if no input field is available
                    lnd_prog%rho_snow_mult_t(:,1,:,:) = lnd_prog%rho_snow_t(:,:,:)
                ENDIF
            END IF ! lmulti_snow

            ! Consistency check for number of soil levels, account for the additional level of t_so compared to smi
            IF (ASSOCIATED(params%requestList%findIconName('t_so'))) THEN
                      levels => params%requestList%getLevels('t_so',jg)
            ELSE
              CALL finish(routine, "t_so not found in input data")
            ENDIF
            IF ( (SIZE(levels)-1) /= nlev_soil ) THEN
              WRITE(message_text,'(a,i3,a,i3,a)') 'Numbers of soil levels in initial data (=',SIZE(levels)-1, &
                &                                  ') and model (=',nlev_soil,') do not match!'
              CALL finish(routine, message_text)
            ENDIF

            ! multi layer fields
            !
            ! Note that either w_so OR smi is written to w_so_t. Which one is required depends
            ! on the initialization mode. Checking grp_vars_fg takes care of this. In case
            ! that smi is read, it is lateron converted to w_so (see smi_to_wsoil)
            SELECT CASE(init_mode)
                CASE(MODE_COMBINED,MODE_COSMO,MODE_ICONVREMAP)
                    CALL fetchTiled3dWithFallback(params          = params,         &
                      &                           varName         = 'smi' ,         &
                      &                           varNameFallback = 'w_so',         &
                      &                           jg              = jg,             &
                      &                           tileCount       = ntiles_total,   &
                      &                           field           = lnd_prog%w_so_t )
                CASE DEFAULT
                    CALL fetchTiled3d(params, 'w_so', jg, ntiles_total, lnd_prog%w_so_t)
            END SELECT
            SELECT CASE(init_mode)
              CASE(MODE_COMBINED,MODE_COSMO)
                ! no w_so_ice available
              CASE DEFAULT
                CALL fetchTiled3d(params, 'w_so_ice', jg, ntiles_total, lnd_prog%w_so_ice_t) ! w_so_ice is re-diagnosed in terra_multlay_init
            END SELECT


            CALL fetchTiled3d(params, 't_so', jg, ntiles_total, lnd_prog%t_so_t)

            ! Skipped in MODE_COMBINED and in MODE_COSMO (i.e. when starting from ICON/COSMO soil)
            ! Instead z0 is re-initialized (see mo_nwp_phy_init)
            CALL fetchSurface(params, 'gz0', jg, prm_diag(jg)%gz0)

            ! first guess for fresh water lake fields
            IF(ASSOCIATED(wtr_prog%t_mnw_lk)) CALL fetchSurface(params, 't_mnw_lk', jg, wtr_prog%t_mnw_lk)
            IF(ASSOCIATED(wtr_prog%t_wml_lk)) CALL fetchSurface(params, 't_wml_lk', jg, wtr_prog%t_wml_lk)
            IF(ASSOCIATED(wtr_prog%h_ml_lk)) CALL fetchSurface(params, 'h_ml_lk', jg, wtr_prog%h_ml_lk)
            IF(ASSOCIATED(wtr_prog%t_bot_lk)) CALL fetchSurface(params, 't_bot_lk', jg, wtr_prog%t_bot_lk)
            IF(ASSOCIATED(wtr_prog%c_t_lk)) CALL fetchSurface(params, 'c_t_lk', jg, wtr_prog%c_t_lk)
            ! The following two variables are actually constants and would not need to be read at all
            ! currently removed from input variable groups, which are defined via add_var
            IF(ASSOCIATED(wtr_prog%t_b1_lk)) CALL fetchSurface(params, 't_b1_lk', jg, wtr_prog%t_b1_lk)
            IF(ASSOCIATED(wtr_prog%h_b1_lk)) CALL fetchSurface(params, 'h_b1_lk', jg, wtr_prog%h_b1_lk)

            IF(iprog_aero >= 1) THEN
                IF (iprog_aero == 1) THEN
                  aerosol_fg_present(jg,1:4) = .FALSE.
                  aerosol_fg_present(jg,5)   = .TRUE.
                ELSE
                  aerosol_fg_present(jg,:)   = .TRUE.
                ENDIF

                my_ptr2d => prm_diag(jg)%aerosol(:,iss,:)
                CALL fetchSurface(params, 'aer_ss', jg, my_ptr2d)
                IF (inputInstructions(jg)%ptr%sourceOfVar('aer_ss') == kInputSourceCold) aerosol_fg_present(jg,iss) = .FALSE.
                !
                my_ptr2d => prm_diag(jg)%aerosol(:,iorg,:)
                CALL fetchSurface(params, 'aer_or', jg, my_ptr2d)
                IF (inputInstructions(jg)%ptr%sourceOfVar('aer_or') == kInputSourceCold) aerosol_fg_present(jg,iorg) = .FALSE.
                !
                my_ptr2d => prm_diag(jg)%aerosol(:,ibc,:)
                CALL fetchSurface(params, 'aer_bc', jg, my_ptr2d)
                IF (inputInstructions(jg)%ptr%sourceOfVar('aer_bc') == kInputSourceCold) aerosol_fg_present(jg,ibc) = .FALSE.
                !
                my_ptr2d => prm_diag(jg)%aerosol(:,iso4,:)
                CALL fetchSurface(params, 'aer_su', jg, my_ptr2d)
                IF (inputInstructions(jg)%ptr%sourceOfVar('aer_su') == kInputSourceCold) aerosol_fg_present(jg,iso4) = .FALSE.
                !
                my_ptr2d => prm_diag(jg)%aerosol(:,idu,:)
                CALL fetchSurface(params, 'aer_du', jg, my_ptr2d)
                IF (inputInstructions(jg)%ptr%sourceOfVar('aer_du') == kInputSourceCold) aerosol_fg_present(jg,idu) = .FALSE.

                ! Reading of each of these variables is purely OPTIONAL, and inform the InputInstructions about this.
                IF (.NOT.aerosol_fg_present(jg,iss))  CALL inputInstructions(jg)%ptr%setSource('aer_ss', kInputSourceCold)
                IF (.NOT.aerosol_fg_present(jg,iorg)) CALL inputInstructions(jg)%ptr%setSource('aer_or', kInputSourceCold)
                IF (.NOT.aerosol_fg_present(jg,ibc))  CALL inputInstructions(jg)%ptr%setSource('aer_bc', kInputSourceCold)
                IF (.NOT.aerosol_fg_present(jg,iso4)) CALL inputInstructions(jg)%ptr%setSource('aer_su', kInputSourceCold)
                IF (.NOT.aerosol_fg_present(jg,idu))  CALL inputInstructions(jg)%ptr%setSource('aer_du', kInputSourceCold)
            ELSE
                aerosol_fg_present(jg,:) = .FALSE.
            END IF

            IF (atm_phy_nwp_config(jg)%lstoch_sde) THEN
               CALL fetchSurface(params, 'clmf_a' , jg, prm_nwp_stochconv(jg)%clmf_a)
               CALL fetchSurface(params, 'clmf_p' , jg, prm_nwp_stochconv(jg)%clmf_p)
               CALL fetchSurface(params, 'clnum_a', jg, prm_nwp_stochconv(jg)%clnum_a)
               CALL fetchSurface(params, 'clnum_p', jg, prm_nwp_stochconv(jg)%clnum_p)
            END IF
            IF (atm_phy_nwp_config(jg)%lstoch_deep) THEN
               CALL fetchSurface(params, 'clmf_d' , jg, prm_nwp_stochconv(jg)%clmf_d)
               CALL fetchSurface(params, 'clnum_d', jg, prm_nwp_stochconv(jg)%clnum_d)
            END IF
        END IF
    END DO

    CALL fetch_dwdfg_jsb(requestList, inputInstructions)

  END SUBROUTINE fetch_dwdfg_sfc


  !! processing of DWD first-guess DATA
  SUBROUTINE process_input_dwdfg_sfc(p_patch, inputInstructions, p_lnd_state, ext_data)
    TYPE(t_patch),             INTENT(IN)      :: p_patch(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT)   :: p_lnd_state(:)
    TYPE(t_external_data)    , INTENT(INOUT)   :: ext_data(:)

    INTEGER :: jg, jt, jb, jc, ic, i_endidx
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_wtr_prog), POINTER :: wtr_prog_now
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    CHARACTER(LEN = *), PARAMETER :: routine = modname//":process_input_dwdfg_sfc"


    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            lnd_prog     => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            wtr_prog_now => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))
            lnd_diag     => p_lnd_state(jg)%diag_lnd
            DO jt=1, ntiles_total + ntiles_water
                ! Copy t_g_t and qv_s_t to t_g and qv_s, respectively, on limited-area domains.
                ! This is needed to avoid invalid operations in the initialization of the turbulence scheme
                IF (jt == 1 .AND. (jg > 1 .OR. l_limited_area) ) THEN
!$OMP PARALLEL DO PRIVATE(jb,jc,i_endidx)
                    ! include boundary interpolation zone of nested domains and halo points
                    DO jb = 1, p_patch(jg)%nblks_c
                        IF (jb == p_patch(jg)%nblks_c) THEN
                            i_endidx = p_patch(jg)%npromz_c
                        ELSE
                            i_endidx = nproma
                        END IF

                        DO jc = 1, i_endidx
                            lnd_prog%t_g(jc,jb)  = lnd_prog%t_g_t(jc,jb,1)
                            lnd_diag%qv_s(jc,jb) = lnd_diag%qv_s_t(jc,jb,1)
                        END DO  ! jc
                    END DO  ! jb
!$OMP END PARALLEL DO
                END IF
            END DO


            ! When starting from ICON or COSMO soil (i.e. MODE_COMBINED or MODE_COSMODE).
            ! SMI is read if available, with W_SO being the fallback option. If SMI is 
            ! read, it is directly stored in w_so_t. Here, it is converted into w_so
            IF (ANY((/MODE_COMBINED,MODE_COSMO,MODE_ICONVREMAP/) == init_mode)) THEN
                IF (inputInstructions(jg)%ptr%sourceOfVar('smi')==kInputSourceFg) THEN
                    DO jt=1, ntiles_total
                        CALL smi_to_wsoil(p_patch(jg), ext_data(jg)%atm%list_land, &
                          &               ext_data(jg)%atm%soiltyp, lnd_prog%w_so_t(:,:,:,jt))
                    END DO
                ENDIF
            END IF


            ! Initialize t_s_t, which is not read in
            !
!$OMP PARALLEL DO PRIVATE(jb,jt,ic,jc)
            DO jb = 1, p_patch(jg)%nblks_c

              ! set t_s_t to t_so_t(1) for land tiles
              DO jt = 1, ntiles_total
                DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
                  jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
                  lnd_prog%t_s_t(jc,jb,jt) = lnd_prog%t_so_t(jc,1,jb,jt)
                ENDDO
              ENDDO  ! ntiles

              IF (llake) THEN
                ! take lake surface temperature from t_wml_lk
                DO ic = 1, ext_data(jg)%atm%list_lake%ncount(jb)
                  jc = ext_data(jg)%atm%list_lake%idx(ic,jb)
                  lnd_prog%t_s_t(jc,jb,isub_lake) = wtr_prog_now%t_wml_lk(jc,jb)
                ENDDO
              ELSE
                ! without lake model, take t_g as an estimate of the lake temperature
                DO ic = 1, ext_data(jg)%atm%list_lake%ncount(jb)
                  jc = ext_data(jg)%atm%list_lake%idx(ic,jb)
                  lnd_prog%t_s_t(jc,jb,isub_lake) = MAX(tmelt, lnd_prog%t_g_t(jc,jb,isub_lake))
                ENDDO
              ENDIF

              IF (itype_canopy == 1 .OR. inputInstructions(jg)%ptr%sourceOfVar('t_sk') == kInputSourceCold) THEN
                DO jt = 1, ntiles_total + ntiles_water
                  IF (jb == p_patch(jg)%nblks_c) THEN
                    i_endidx = p_patch(jg)%npromz_c
                  ELSE
                    i_endidx = nproma
                  END IF
                  DO jc = 1, i_endidx
                    lnd_prog%t_sk_t(jc,jb,jt) = lnd_prog%t_s_t(jc,jb,jt)
                  ENDDO
                ENDDO
              ELSE
                ! Ensure that t_sk is filled with meaningful values on all existing land points (may not be the case after changing extpar data)
                DO jt = 1, ntiles_total
                  DO ic = 1, ext_data(jg)%atm%lp_count_t(jb,jt)
                    jc = ext_data(jg)%atm%idx_lst_lp_t(ic,jb,jt)
                    IF (lnd_prog%t_sk_t(jc,jb,jt) < 100._wp) lnd_prog%t_sk_t(jc,jb,jt) = lnd_prog%t_s_t(jc,jb,jt)
                  ENDDO
                ENDDO  ! ntiles
              ENDIF

              ! NOTE: Initialization of sea-water and sea-ice tiles 
              ! is done later in mo_nwp_sfc_utils:nwp_surface_init for 
              ! two reasons:
              ! I)  index lists for sea-ice and sea-water are 
              !     not yet available at this point.
              ! II) If an SST-analysis is read in, t_s_t(isub_water) 
              !     must be updated another time, anyway.

            ENDDO  ! jb
!$OMP END PARALLEL DO

        END IF

    END DO

    CALL process_input_dwdfg_jsb
  END SUBROUTINE process_input_dwdfg_sfc


  !>
  !! Fetch DWD analysis DATA from the request list (land/surface only)
  !!
  !! Analysis is read for:
  !! fr_seaice, t_ice, h_ice, t_so, h_snow, w_snow, w_i, t_snow, rho_snow, freshsnow, w_so
  !!
  !!
  SUBROUTINE fetch_dwdana_sfc(requestList, p_patch, p_lnd_state, initicon, inputInstructions)
    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_patch),             INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state), TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_initicon_state), INTENT(INOUT), TARGET :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    INTEGER :: jg
    INTEGER, PARAMETER :: jt = 1
    REAL(wp), POINTER :: my_ptr2d(:,:)
    REAL(wp), POINTER :: my_ptr3d(:,:,:)
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag
    TYPE(t_wtr_prog), POINTER :: wtr_prog
    TYPE(t_fetchParams) :: params
    LOGICAL :: lHaveFg

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdana_sfc'

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .FALSE.

    DO jg = 1, n_dom
        IF(p_patch(jg)%ldom_active) THEN
            IF (ANY((/MODE_IAU, MODE_IAU_OLD /) == init_mode) .AND. lp2cintp_sfcana(jg)) CYCLE  ! Skip this domain if it will be interpolated from its parent domain in process_input_dwdana_sfc()

            lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
            lnd_diag => p_lnd_state(jg)%diag_lnd
            wtr_prog => p_lnd_state(jg)%prog_wtr(nnow_rcf(jg))

            ! sea-ice fraction
            CALL fetchSurface(params, 'fr_seaice', jg, lnd_diag%fr_seaice)
            ! sea-ice temperature
            CALL fetchSurface(params, 't_ice', jg, wtr_prog%t_ice)
            ! sea-ice height
            CALL fetchSurface(params, 'h_ice', jg, wtr_prog%h_ice)

            ! sea-surface temperature: fetch T_SEA or, alternatively, T_SO(0)
            my_ptr2d => initicon(jg)%sfc%sst(:,:)
            CALL fetch2d(params, 't_seasfc', 0.0_wp, jg, my_ptr2d)
            !
            IF (inputInstructions(jg)%ptr%fetchStatus('t_seasfc', lIsFg=.FALSE.) == kStateFailedFetch) THEN
              ! T_SO(0). Note that the file may contain a 3D field, of which we ONLY fetch the level at 0.0.
              lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('t_so') == kInputSourceFg
              CALL fetch2d(params, 't_so', 0.0_wp, jg, my_ptr2d)
              ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
              IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('t_so') == kInputSourceAna) THEN
                  CALL inputInstructions(jg)%ptr%setSource('t_so', kInputSourceBoth)
              END IF
            ENDIF


            ! h_snow
            IF ( init_mode == MODE_IAU ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('h_snow') == kInputSourceFg
                my_ptr2d => initicon(jg)%sfc_inc%h_snow(:,:)
                CALL fetchSurface(params, 'h_snow', jg, my_ptr2d)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('h_snow') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('h_snow', kInputSourceBoth)
                END IF
            ELSE
                my_ptr2d => lnd_diag%h_snow_t(:,:,jt)
                CALL fetchSurface(params, 'h_snow', jg, my_ptr2d)
            ENDIF

            ! w_snow
            my_ptr2d => lnd_prog%w_snow_t(:,:,jt)
            CALL fetchSurface(params, 'w_snow', jg, my_ptr2d)

            ! w_i
            my_ptr2d => lnd_prog%w_i_t(:,:,jt)
            CALL fetchSurface(params, 'w_i', jg, my_ptr2d)

            ! t_snow
            my_ptr2d => lnd_prog%t_snow_t(:,:,jt)
            CALL fetchSurface(params, 't_snow', jg, my_ptr2d)

            ! rho_snow
            my_ptr2d => lnd_prog%rho_snow_t(:,:,jt)
            CALL fetchSurface(params, 'rho_snow', jg, my_ptr2d)

            ! freshsnow
            IF ( init_mode == MODE_IAU ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('freshsnow') == kInputSourceFg
                my_ptr2d => initicon(jg)%sfc_inc%freshsnow(:,:)
                CALL fetchSurface(params, 'freshsnow', jg, my_ptr2d)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('freshsnow') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('freshsnow', kInputSourceBoth)
                END IF
            ELSE
                my_ptr2d => lnd_diag%freshsnow_t(:,:,jt)
                CALL fetchSurface(params, 'freshsnow', jg, my_ptr2d)
            ENDIF

            ! w_so
            IF ( (init_mode == MODE_IAU) .OR. (init_mode == MODE_IAU_OLD) ) THEN
                lHaveFg = inputInstructions(jg)%ptr%sourceOfVar('w_so') == kInputSourceFg
                my_ptr3d => initicon(jg)%sfc_inc%w_so(:,:,:)
                CALL fetch3d(params, 'w_so', jg, my_ptr3d)
                ! check whether we are using DATA from both FG AND ANA input, so that it's correctly listed IN the input source table
                IF(lHaveFg.AND.inputInstructions(jg)%ptr%sourceOfVar('w_so') == kInputSourceAna) THEN
                    CALL inputInstructions(jg)%ptr%setSource('w_so', kInputSourceBoth)
                END IF
            ELSE
                my_ptr3d => lnd_prog%w_so_t(:,:,:,jt)
                CALL fetch3d(params, 'w_so', jg, my_ptr3d)
            END IF

            ! t_2m bias
            IF ((icpl_da_sfcevap == 1 .OR. icpl_da_sfcevap == 2) .AND. init_mode == MODE_IAU ) THEN
               my_ptr2d => initicon(jg)%sfc_inc%t_2m(:,:)
               CALL fetchSurface(params, 't_2m', jg, my_ptr2d)
            ENDIF

        END IF
    END DO ! loop over model domains

  END SUBROUTINE fetch_dwdana_sfc


  ! Fill remaining tiles for snow variables if tile approach is used
  ! Only fields that are actually read from the snow analysis are copied; note that MODE_IAU is mandatory when using tiles
  SUBROUTINE process_input_dwdana_sfc (p_patch, p_lnd_state, initicon, inputInstructions)
    TYPE(t_patch),                  INTENT(IN)    :: p_patch(:)
    TYPE(t_lnd_state),      TARGET, INTENT(INOUT) :: p_lnd_state(:)
    TYPE(t_initicon_state),         INTENT(INOUT) :: initicon(:)
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(n_dom)

    INTEGER :: jg, jt, jb, jc, i_endidx
    TYPE(t_lnd_prog), POINTER :: lnd_prog
    TYPE(t_lnd_diag), POINTER :: lnd_diag

    jt = 1
    DO jg = 1, n_dom
      IF(p_patch(jg)%ldom_active .AND. ANY((/MODE_IAU, MODE_IAU_OLD /) == init_mode)) THEN
        IF (lp2cintp_sfcana(jg)) THEN
          ! Perform parent-to-child interpolation of surface fields read from the analysis
          CALL interpolate_sfcana(initicon, inputInstructions, p_patch(jg)%parent_id, jg, p_lnd_state(:))

        ELSE IF (ntiles_total>1 .AND. init_mode == MODE_IAU_OLD) THEN
          ! MODE_IAU_OLD: H_SNOW, FRESHSNOW, W_SNOW and RHO_SNOW are read from analysis (full fields)
          ! Since only the first tile index is filled (see above), we fill (here) the remaining tiles
          ! if tile approach is used. Note that for MODE_IAU this copy is skipped on purpose,
          ! since for ltile_coldstart=.FALSE. tile information would be overwritten.
          lnd_prog => p_lnd_state(jg)%prog_lnd(nnow_rcf(jg))
          lnd_diag => p_lnd_state(jg)%diag_lnd

!$OMP PARALLEL DO PRIVATE(jb,jc,jt,i_endidx)
          DO jb = 1, p_patch(jg)%nblks_c
            i_endidx = MERGE(nproma, p_patch(jg)%npromz_c, &
              &              jb /= p_patch(jg)%nblks_c)
            DO jt = 2, ntiles_total
              DO jc = 1, i_endidx
                lnd_diag%freshsnow_t(jc,jb,jt) = lnd_diag%freshsnow_t(jc,jb,1)
                lnd_diag%h_snow_t(jc,jb,jt)    = lnd_diag%h_snow_t(jc,jb,1)
                lnd_prog%w_snow_t(jc,jb,jt)    = lnd_prog%w_snow_t(jc,jb,1)
                lnd_prog%rho_snow_t(jc,jb,jt)  = lnd_prog%rho_snow_t(jc,jb,1)
              END DO
            END DO
          END DO
!$OMP END PARALLEL DO
        END IF
      END IF
    END DO ! loop over model domains
  END SUBROUTINE process_input_dwdana_sfc


  SUBROUTINE fetch_dwdfg_jsb(requestList, inputInstructions)

    USE mo_impl_constants, ONLY: REAL_T, LSS_JSBACH
    USE mo_var_groups, ONLY: var_groups_dyn
    USE mo_var_list_register, ONLY: t_vl_register_iter
    USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config

    CLASS(t_InputRequestList), POINTER, INTENT(INOUT) :: requestList
    TYPE(t_readInstructionListPtr), INTENT(INOUT) :: inputInstructions(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':fetch_dwdfg_jsb'

    TYPE(t_fetchParams) :: params
    TYPE(t_vl_register_iter) :: iter
    INTEGER :: grp
    INTEGER :: i
    INTEGER :: jg

    INTEGER :: nvars, nfound

    LOGICAL :: found

    IF (.NOT. ANY(atm_phy_nwp_config(:)%inwp_surface == LSS_JSBACH)) RETURN

    grp = var_groups_dyn%group_id('jsb_init_vars')

    ALLOCATE(params%inputInstructions(SIZE(inputInstructions, 1)))
    params%inputInstructions = inputInstructions
    params%requestList => requestList
    params%routine = routine
    params%isFg = .TRUE.

    nvars = 0
    nfound = 0

    DO WHILE (iter%next())
      jg = iter%cur%p%patch_id

      IF (atm_phy_nwp_config(jg)%inwp_surface /= LSS_JSBACH) CYCLE

      ASSOCIATE (vl => iter%cur%p%vl)
        DO i = 1, iter%cur%p%nvars
          IF (.NOT. vl(i)%p%info%in_group(grp)) CYCLE

          found = .FALSE.
          nvars = nvars + 1

          IF (vl(i)%p%info%data_type /= REAL_T) THEN
            WRITE(message_text,'(a,a,a,i2)') 'Variable "', TRIM(vl(i)%p%info%name), '" has unhandled type: ', vl(i)%p%info%data_type
            CALL warning(routine, message_text)
            CYCLE
          END IF

          SELECT CASE (vl(i)%p%info%ndims)
          CASE (2)
            CALL fetchSurface(params, TRIM(vl(i)%p%info%name), jg, vl(i)%p%r_ptr(:,:,1,1,1), found=found)
          CASE (3)
            CALL fetch3d(params, TRIM(vl(i)%p%info%name), jg, vl(i)%p%r_ptr(:,:,:,1,1), found=found)
          CASE DEFAULT
            WRITE(message_text,'(a,a,a,i2)') 'Variable "', TRIM(vl(i)%p%info%name), '" has unhandled number of dimensions: ', vl(i)%p%info%ndims
            CALL warning(routine, message_text)
          END SELECT

          CALL params%inputInstructions(jg)%ptr%handleError(found, TRIM(vl(i)%p%info%name), params%routine, params%isFg)

          IF (found) THEN
            nfound = nfound + 1
          END IF
        END DO
      END ASSOCIATE
    END DO

    WRITE(message_text,'(a,i4,a,i4,a)') 'Fetching JSBACH fields... Found ', nfound, ' of ', nvars, ' vars.'
    CALL message(routine, message_text)

  END SUBROUTINE fetch_dwdfg_jsb

  SUBROUTINE process_input_dwdfg_jsb
#ifndef __NO_JSBACH__
    USE mo_impl_constants, ONLY: LSS_JSBACH
    USE mo_jsb_model_init, ONLY: jsbach_init_after_restart
    USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':process_input_dwdfg_jsb'

    INTEGER :: jg

    IF (.NOT. ANY(atm_phy_nwp_config(:)%inwp_surface == LSS_JSBACH)) RETURN

    CALL message(routine, 'Sanitizing JSBACH fields.')

    DO jg = 1, n_dom
      IF (atm_phy_nwp_config(jg)%inwp_surface /= LSS_JSBACH) CYCLE

      ! From JSBACH's point of view, we are doing a restart because we replace its initial state
      ! with the state read in from the FG file.
      ! The main purpose of this routine here is processing some of the fields just read in. If the
      ! FG file has less than double precision, rounding errors could lead to constraint violations
      ! and model aborts. The routine checks the constraints and corrects the fields.
      CALL jsbach_init_after_restart(jg)
    END DO

#endif
  END SUBROUTINE process_input_dwdfg_jsb

END MODULE mo_initicon_io

