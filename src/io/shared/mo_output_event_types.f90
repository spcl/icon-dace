! Type definitions for handling of regular output steps and ready
! file events on multiple I/O PEs.
!
! See "mo_output_event_handler" for a detailed description.
!
! -----------------------------------------------------------------------------------!
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

MODULE mo_output_event_types
  USE mtime,                 ONLY: datetime
  USE mo_kind,               ONLY: wp
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_sim_step_info
  PUBLIC :: t_event_data
  PUBLIC :: t_event_step_data
  PUBLIC :: t_event_step
  PUBLIC :: t_output_event
  PUBLIC :: t_par_output_event
  PUBLIC :: MAX_FILENAME_STR_LEN
  PUBLIC :: MAX_EVENT_NAME_STR_LEN
  PUBLIC :: DEFAULT_EVENT_NAME


  !> max. length of filename
  INTEGER, PARAMETER :: MAX_FILENAME_STR_LEN   =  256

  !> max. length of output event name
  INTEGER, PARAMETER :: MAX_EVENT_NAME_STR_LEN = 1024

  !>  default event name (ie. events that DO NOT write ready files
  CHARACTER (LEN=*), PARAMETER :: DEFAULT_EVENT_NAME = "default"
  

  !---------------------------------------------------------------
  ! types

  !---------------------------------------------------------------
  ! DEFINITIONS FOR CONVERSION "time stamp -> simulation step"

  !> Data structure containing all necessary data for mapping a time
  !! stamp onto a corresponding simulation step index.
  !!
  !! These data members correspond to the master_time_control_nml in
  !! the following way:
  !!
  !!   ***************************  EXPERIMENT RUN   *************
  !!   *                                                         *
  !!   *                    ***  JOB RUN   ***                   *
  !!   *                    *                *                   *
  !! --|--------------------[----------------]-------------------|---------------> (time axis)
  !!   ^                    ^                ^                   ^
  !!   tc_exp_startdate     tc_startdate     tc_stopdate         tc_exp_stopdate   (master_time_control_nml notation)
  !!
  !!   sim_start            run_start        restart_time        sim_end           (t_sim_step_info)
  !!
  !!
  !!
  TYPE t_sim_step_info
    !> simulation start/end time stamp
    TYPE(datetime) :: sim_start, sim_end
    !> start of this run (-> restart)
    TYPE(datetime) :: run_start
    !> end of this run (-> restart)
    TYPE(datetime) :: restart_time
    !> model domain start/end time
    TYPE(datetime) :: dom_start_time, dom_end_time
    !> [s] length of a time step
    REAL(wp)       :: dtime
    !> initial time loop counter (important for restart)
    INTEGER        :: jstep0
  END TYPE t_sim_step_info


  !---------------------------------------------------------------
  ! EVENT-SPECIFIC DATA
  !
  !  We collect these members in derived data types in order to
  !  separate the trigger mechanism from the event (the handling of
  !  output steps) itself.

  !> Container for event-specific data, for all steps of an event.
  !
  TYPE t_event_data
    CHARACTER(LEN=MAX_EVENT_NAME_STR_LEN) :: name                             !< output event name
    !> simulation start
    TYPE(datetime)                        :: sim_start
  END TYPE t_event_data


  !> Container for event-specific data, for a single event step.
  !
  !  @note Note that it is possible to have different date-time stamps
  !        for a single model step: This happens when output events
  !        that were chosen at different intervals are triggered at
  !        the same mode step which is the earliest dynamics/advection
  !        step available.
  !
  TYPE t_event_step_data
    !> rank of participating PE
    INTEGER                               :: i_pe
    !> file counter, "part of file" counter
    INTEGER                               :: jfile, jpart
    !> Flag. .TRUE. if file is to be opened in this step
    LOGICAL                               :: l_open_file

    !> tag, e.g. for MPI isend/irecv messages, unique at least on sender side
    INTEGER                               :: i_tag
    !> mtime time stamp
    TYPE(datetime) :: datetime
    !> output file name
    CHARACTER(LEN=MAX_FILENAME_STR_LEN)   :: filename_string
  END TYPE t_event_step_data


  !---------------------------------------------------------------
  ! DEFINITION OF EVENT STEPS
  
  !> Single step of an event.
  !
  !  Events are triggered by the simulation step (INTEGER) only to
  !  avoid ambiguities.  The corresponding event time stamp string
  !  must not necessarily match this simulation step exactly.
  !
  TYPE t_event_step
    INTEGER                               :: i_sim_step                       !< simulation step that triggers event
    INTEGER                               :: n_pes                            !< no. PEs participating in this step
    !> exact (model step) time stamp
    TYPE(datetime) :: exact_date
    !
    TYPE(t_event_step_data), ALLOCATABLE  :: event_step_data(:)               !< event data for each PE
  END TYPE t_event_step

  
  !> List of steps for an event.
  !
  TYPE t_output_event
    TYPE(t_event_data)                    :: event_data                       !< event data
    INTEGER                               :: n_event_steps                    !< total no. of event steps
    INTEGER                               :: i_event_step                     !< current event step
    TYPE(t_event_step), ALLOCATABLE       :: event_step(:)                    !< event steps (1,...,n_event_steps)
  END TYPE t_output_event


  !---------------------------------------------------------------
  ! MPI-PARALLEL EVENTS
 
  !> List of steps of an event that is performed on several PEs in
  !  parallel.
  !
  TYPE t_par_output_event
    TYPE(t_output_event), ALLOCATABLE     :: output_event                     !< event data structure

    ! --- MPI related fields:
    INTEGER                               :: icomm                            !< MPI communicator
    INTEGER                               :: iroot                            !< MPI rank of root PE
    INTEGER                               :: isend_req                        !< MPI isend request token
    INTEGER                               :: isend_buf                        !< MPI isend data buffer
    INTEGER                               :: irecv_nreq                       !< no. of MPI irecv requests
    INTEGER, ALLOCATABLE                  :: irecv_req(:)                     !< MPI irecv request tokens
    INTEGER, ALLOCATABLE                  :: irecv_buf(:)                     !< MPI irecv data buffer

    TYPE(t_par_output_event), POINTER     :: next                             !< neighbor in linked list
  END TYPE t_par_output_event

END MODULE mo_output_event_types
