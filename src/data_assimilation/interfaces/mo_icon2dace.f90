! Interface DACE observation operators from ICON
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

MODULE mo_icon2dace

!-----------------------------------------------------------------------------
!
! Description:
!   Interface DACE observation operators from ICON
!
!-----------------------------------------------------------------------------
  !-------------
  ! Modules used
  !-------------

  !-------------
  ! ICON modules
  !-------------
  USE mo_kind,           ONLY: wp                    ! working precision kind
  USE mo_exception,      ONLY: finish, message
  USE mo_namelist,       ONLY: open_nml,            &! open namelist file
                               close_nml,           &! close namelist file
                               position_nml,        &! position to nml group
                               nnml,                &! Fortran namelist unit
                               POSITIONED            ! position_nml: OK return flag
  USE mo_util_vcs,       ONLY: util_repository_url, &!
                               util_branch_name,    &!
                               util_revision_key     !
  USE mo_grid_config,    ONLY: nroot, start_lev
  USE mo_gribout_config, ONLY: gribout_config
  USE mo_parallel_config,ONLY: nproma, idx_1d
  USE mo_model_domain,   ONLY: p_patch,             &!
               t_patch_icon => t_patch
  USE mo_communication,  ONLY: t_comm_pattern, t_comm_gather_pattern
  USE mo_loopindices,    ONLY: get_indices_c, get_indices_v
  USE mo_impl_constants, ONLY: min_rlcell, min_rlvert, min_rledge
  USE mo_sync,           ONLY: SYNC_C, SYNC_V, sync_patch_array

  USE mo_time_config,    ONLY: time_config
  USE mtime,             ONLY: &
                         MD => MAX_DATETIME_STR_LEN
  USE mo_util_mtime,     ONLY: getElapsedSimTimeInSeconds

  USE mo_ext_data_state, ONLY: ext_data       ! fr_land etc.
  USE mo_nonhydro_state, ONLY: p_nh_state
  USE mo_nonhydro_types, ONLY: t_nh_prog, t_nh_diag
  USE mo_nwp_phy_state,  ONLY: prm_diag
  USE mo_nwp_phy_types,  ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_state,  ONLY: p_lnd_state
  USE mo_nwp_lnd_types,  ONLY: t_lnd_diag, t_lnd_prog
  USE mo_run_config,     ONLY: iqv, iqc, iqi, iqg, iqr, iqs
  USE mo_dynamics_config,ONLY: nnow, nnow_rcf
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info
  USE mo_assimilation_config,ONLY: assimilation_config
  USE mo_nh_diagnose_pres_temp, ONLY: diagnose_pres_temp
  USE mo_intp_rbf,              ONLY: rbf_vec_interpol_cell
  USE mo_intp_data_strc,        ONLY: p_int_state
  USE mo_atm_phy_nwp_config,    ONLY: atm_phy_nwp_config 
  !-------------------
  ! ICON event control
  !-------------------
  USE mtime,             ONLY: datetime,  newDatetime,  deallocateDatetime, &
       &                       timedelta, newTimedelta, deallocateTimedelta,&
       &                       datetimetostring, timedeltaToString,         &
       &                       MAX_DATETIME_STR_LEN, MAX_TIMEDELTA_STR_LEN, &
       &                       MAX_MTIME_ERROR_STR_LEN, mtime_strerror,     &
       &                       OPERATOR(+), OPERATOR(-), OPERATOR(>=),      &
       &                       ASSIGNMENT(=), event, eventGroup, newEvent,  &
       &                       addEventToEventGroup, isCurrentEventActive
  USE mo_event_manager,  ONLY: addEventGroup, getEventGroup, printEventGroup


  IMPLICIT NONE
  !----------------
  ! Public entities
  !----------------
  PRIVATE
  PUBLIC :: init_dace      ! initialise DACE as a subsystem
  PUBLIC :: init_dace_op   ! initialise DACE observation operators
  PUBLIC :: run_dace_op    ! run the DACE observation operators for a time step
  PUBLIC :: finish_dace    ! write fof-Files, clean up DACE
  PUBLIC :: mec_event      ! event handle for invoking the MEC
  PUBLIC :: dace_op_init   ! DACE operators initialized?

  !-----------------
  ! Module variables
  !-----------------
  logical, save        :: dace_op_init = .false. ! DACE operators initialized?
  PROTECTED            :: dace_op_init
  !----------------------
  ! Event control (timer)
  !----------------------
  TYPE(event),    POINTER :: mec_Event   => NULL()
  TYPE(datetime), POINTER :: mec_RefDate => NULL()
  !----------------------------------------------------------------------------
  !============================================================================
contains
  !============================================================================
  subroutine init_dace (comm, p_io, ldetached)
    !----------------------------------------
    ! set the DACE MPI communicator from ICON
    ! initialize domain decomposition, grid
    !----------------------------------------
    integer ,intent(in) :: comm  ! communicator to use
    integer ,intent(in) :: p_io  ! PE to use for I/O
    logical, intent(in) :: ldetached ! indicates that current PE is a detached IO-PE
    !----------------
    ! Local variables
    !----------------
    integer       :: n
    integer       :: ios     ! I/O status
    integer       :: iu      ! temp. unit number
    character(MD) :: cdate   ! yyyy-mm-ddThh:mm:ss.000
    character(14) :: adate   ! yyyymmddhhmmss
    character(14) :: refdate ! yyyymmddhhmmss

    call message ("","")
    call message ("icon2dace","initializing DACE coupling")

    CALL finish ("init_dace","DACE coupling requested but not compiled in")

  end subroutine init_dace
  !============================================================================
  !============================================================================
  subroutine init_dace_op ()
    !-----------------------------------------------------------------------
    ! initialise the DACE observation operators from ICON
    ! 1) set grid according to the specifications by ICON
    ! 2) read relevant namelists in DACE
    ! 3) read observation CDFIN files, set up operators
    ! 4) distribute observations over PEs, set up interpolation coefficients
    ! 5) return information on required time steps to ICON
    !-----------------------------------------------------------------------
    CALL finish ("init_dace_op","DACE coupling requested but not compiled in")

  end subroutine init_dace_op
  !============================================================================
  !============================================================================
  subroutine run_dace_op (mtime_current)
    TYPE(datetime), POINTER :: mtime_current    !< current datetime (mtime)
    !---------------------------------------------------------
    ! run the dace observation operators for a given time step
    !---------------------------------------------------------
    CALL finish ("run_dace_op","DACE coupling requested but not compiled in")

  end subroutine run_dace_op
  !============================================================================
  subroutine finish_dace ()
    !-------------------
    ! 1) write fof-Files
    ! 2) clean up DACE
    !-------------------
    CALL finish ("finish_dace","DACE coupling requested but not compiled in")

  end subroutine finish_dace
  !============================================================================
  !============================================================================

end module mo_icon2dace
