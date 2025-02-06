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

! Contains code for age tacer dynamics

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_ocean_age_tracer
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_ocean_nml,           ONLY: green_start_date, green_stop_date,        &
                                  & green_duration,                           &
                                  & i_sea_ice,                                &
                                  & age_tracer_inv_relax_time,                &
                                  & l_relaxage_ice,                           &
                                  & age_idx, green_idx, diagnose_age,         &
                                  & diagnose_green
  USE mo_run_config,          ONLY: dtime
  USE mo_model_domain,        ONLY: t_patch, t_patch_3d
  USE mo_dynamics_config,     ONLY: nnew
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_sea_ice_types,       ONLY: t_sea_ice
  USE mo_timer,               ONLY: ltimer, timer_start, timer_stop, &
                                  & timer_extra10, timer_extra11
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device
  USE mo_time_config,         ONLY: time_config
  USE mtime,                  ONLY: datetime, datetimeToString,               &
                                  & deallocateDatetime, MAX_DATETIME_STR_LEN, &
                                  & OPERATOR(-), OPERATOR(+), OPERATOR(>),    &
                                  & OPERATOR(*), ASSIGNMENT(=), OPERATOR(==), &
                                  & OPERATOR(>=), OPERATOR(<=), OPERATOR(/=), &
                                  & OPERATOR(<)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_age_tracer

CONTAINS

  SUBROUTINE calc_age_tracer(patch_3d, ocean_state, p_ice, lacc)
    TYPE(t_patch_3d ), TARGET, INTENT(in)         :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET, INTENT(in) :: ocean_state
    TYPE (t_sea_ice),              INTENT(IN)     :: p_ice
    LOGICAL, INTENT(in), OPTIONAL                 :: lacc
    REAL(wp), DIMENSION(:,:,:), POINTER           :: age_tracer
    REAL(wp), DIMENSION(:,:,:), POINTER           :: green_tracer
    TYPE(datetime), POINTER                       :: ocean_current_time
    REAL(wp)                                      :: relax_strength
    REAL(wp)                                      :: green_target
    CHARACTER(LEN=MAX_DATETIME_STR_LEN)           :: current_date_string
    CHARACTER(LEN=30)                             :: target_string
    INTEGER                                       :: jc, jb, jk, elev
    INTEGER                                       :: i_startidx_c, i_endidx_c
    TYPE(t_subset_range), POINTER                 :: all_cells
    TYPE(t_patch), POINTER                        :: p_patch
    LOGICAL                                       :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    ! --- calendar
    ocean_current_time  => time_config%tc_current_date

    ! --- age_tracer
    IF (diagnose_age) THEN
        age_tracer => ocean_state%p_prog(nnew(1))%tracer(:,:,:,age_idx)
        !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        age_tracer(:,:,:) = age_tracer(:,:,:) + patch_3d%wet_c(:,:,:) * dtime    
        !$ACC END KERNELS
    END IF

    ! --- green_tracer  
    IF (diagnose_green) THEN
      green_tracer => ocean_state%p_prog(nnew(1))%tracer(:,:,:,green_idx)
      green_target = 0
      IF (ocean_current_time >= green_start_date) THEN
        IF (ocean_current_time < green_stop_date) THEN
          ! We are in the impulse phase: time to Greeeennnnn
          green_target = 1._wp / green_duration
        END IF
      END IF
    END IF
    
    ! --- Surface relaxation phase
    p_patch   => patch_3D%p_patch_2D(1)
    all_cells => p_patch%cells%all
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = i_startidx_c, i_endidx_c

        ! --- Calculate the relaxation strength
        IF (l_relaxage_ice .AND. i_sea_ice >= 1) THEN
          relax_strength = (1.0_wp - p_ice%concsum(jc,jb)) * age_tracer_inv_relax_time  * dtime
        ELSE
          relax_strength = age_tracer_inv_relax_time  * dtime
        END IF
        
        ! --- relax the age to zero
        IF (diagnose_age) THEN
            age_tracer(jc,1,jb) = (1._wp - relax_strength) * age_tracer(jc,1,jb)
        END IF

        ! --- relax the Green's tracer to the target
        IF (diagnose_green) THEN
          green_tracer(jc,1,jb) = green_tracer(jc,1, jb) - relax_strength * (green_tracer(jc,1,jb) - green_target)
        END IF
      
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
    
    ! --- Check everything is still positive
    DO jb = all_cells%start_block, all_cells%end_block
      elev = p_patch%nlev
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jk = 1, elev
        DO jc = i_startidx_c, i_endidx_c
          
          ! --- age_tracer
          IF (diagnose_age) THEN
            IF (age_tracer(jc,jk,jb) < 0) THEN
                age_tracer(jc,jk,jb) = 0
            END IF
          END IF

          ! --- green_tracer
          IF (diagnose_green) THEN
            IF (green_tracer(jc,jk,jb) < 0)  THEN
              green_tracer(jc,jk,jb) = 0
            END IF
          END IF

        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END DO
    !$ACC WAIT(1)
    
    ! --- Clean up the pointers
    NULLIFY(age_tracer)
    NULLIFY(green_tracer)
    NULLIFY(ocean_current_time)
    NULLIFY(all_cells)
    NULLIFY(p_patch)
    
  END SUBROUTINE calc_age_tracer

END MODULE mo_ocean_age_tracer
