! Provide an implementation of the ocean tracer nuding functionality.
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

!----------------------------
#include "omp_definitions.inc"
#include "icon_definitions.inc"
!----------------------------
MODULE mo_ocean_nudging
!-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  USE mo_kind,                      ONLY: wp
  USE mo_impl_constants,            ONLY: sea_boundary, sea, min_dolic
  USE mo_math_constants,            ONLY: pi
  USE mo_ocean_nml,                 ONLY: n_zlev, no_tracer,              &
    & threshold_min_t, threshold_max_t, threshold_min_s, threshold_max_s, &
    & type_3dimrelax_temp, para_3dimrelax_temp,                           &
    & type_3dimrelax_salt, para_3dimrelax_salt!,                           &
  USE mo_util_dbg_prnt,             ONLY: dbg_print
  USE mo_parallel_config,           ONLY: nproma
  USE mo_dynamics_config,           ONLY: nold, nnew
  USE mo_run_config,                ONLY: dtime, ltimer, debug_check_level
  USE mo_ocean_types,               ONLY: t_hydro_ocean_state
  USE mo_ocean_nudging_types,       ONLY: t_ocean_nudge
  USE mo_model_domain,              ONLY: t_patch, t_patch_3d
  USE mo_exception,                 ONLY: finish !, message_text, message
  USE mo_operator_ocean_coeff_3d,   ONLY: t_operator_coeff
  USE mo_grid_subset,               ONLY: t_subset_range, get_index_range
  USE mo_fortran_tools,             ONLY: set_acc_host_or_device


  IMPLICIT NONE

  TYPE(t_ocean_nudge)  :: ocean_nudge

  PRIVATE

  CHARACTER(LEN=12)           :: str_module = 'oceTracer   '  ! Output of module for 1 line debug
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug

  ! Public interface
  PUBLIC :: nudge_ocean_tracers, ocean_nudge


CONTAINS

  !-------------------------------------------------------------------------
  !
  !>
  !! !  SUBROUTINE nudge saLINITY AND TEMERATURE
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE nudge_ocean_tracers(patch_3d, p_os, lacc)
    TYPE(t_patch_3d ),TARGET, INTENT(inout)      :: patch_3d
    TYPE(t_hydro_ocean_state), TARGET :: p_os
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    !Local variables

    REAL(wp) :: z_relax
    REAL(wp) :: z_c(nproma,n_zlev,patch_3d%p_patch_2d(1)%alloc_cell_blocks)
    INTEGER :: jc, blockNo, start_index, end_index, level
    TYPE(t_subset_range), POINTER :: cells_in_domain
    TYPE(t_patch), POINTER :: patch_2D
    INTEGER, POINTER :: dolic_c(:,:)
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    patch_2d  => patch_3d%p_patch_2d(1)
    cells_in_domain => patch_2d%cells%ALL
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c

    ! Final step: 3-dim temperature relaxation
    !  - strict time constant, i.e. independent of layer thickness
    !  - additional forcing Term F_T = -1/tau(T-T*) [ K/s ]
    !    when using the sign convention
    !      dT/dt = Operators + F_T
    !    i.e. F_T <0 for  T-T* >0 (i.e. decreasing temperature if it is warmer than relaxation data)
    !  - discretized:
    !    tracer = tracer - 1/(para_3dimRelax_Temp[months]) * (tracer(1)-data_3dimRelax_Temp)
    IF (no_tracer>=1 .AND. type_3dimrelax_temp >0) THEN

      ! calculate relaxation term
      z_relax = 1.0_wp/(para_3dimrelax_temp*2.592e6_wp)
!      ocean_nudge%forc_3dimrelax_temp(:,:,:) = -z_relax * ocean_nudge%relax_3dim_coefficient(:,:,:) &
!        & * ( p_os%p_prog(nnew(1))%tracer(:,:,:,1) - ocean_nudge%data_3dimrelax_temp(:,:,:))

    DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
      CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      DO jc = start_index, end_index
        DO level = 1, dolic_c(jc,blockNo)
          ocean_nudge%forc_3dimrelax_temp(jc,level,blockNo) = z_relax * dtime &
            * (ocean_nudge%data_3dimrelax_temp(jc,level,blockNo) - p_os%p_prog(nnew(1))%tracer(jc,level,blockNo,1) )

          ! add relaxation term to new temperature
          p_os%p_prog(nnew(1))%tracer(jc,level,blockNo,1) = p_os%p_prog(nnew(1))%tracer(jc,level,blockNo,1) + &
                                                ocean_nudge%forc_3dimRelax_temp(jc,level,blockNo)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('3d_relax: tracer T forc', ocean_nudge%forc_3dimRelax_Temp, str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('3d_relax: tracer T data', ocean_nudge%data_3dimRelax_Temp, str_module,idt_src, in_subset=cells_in_domain)      
      z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,1)
      CALL dbg_print('3d_relax: tracer T trac', z_c, str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF

    ! Final step: 3-dim salinity relaxation
    !  - additional forcing Term F_S = -1/tau(S-S*) [ psu/s ]
    !    when using the sign convention
    !      dS/dt = Operators + F_S
    !    i.e. F_S <0 for  S-S* >0 (i.e. decreasing salinity if it is larger than relaxation data)
    !  - discretized:
    !    tracer = tracer - 1/(para_3dimRelax_salt[months]) * (tracer(1)-data_3dimRelax_salt)
    IF (no_tracer>=2 .AND. type_3dimrelax_salt >0) THEN

      ! calculate relaxation term
      z_relax = 1.0_wp/(para_3dimrelax_salt*2.592e6_wp)
!      ocean_nudge%forc_3dimrelax_salt(:,:,:) = -z_relax * ocean_nudge%relax_3dim_coefficient(:,:,:) &
!        & ( p_os%p_prog(nnew(1))%tracer(:,:,:,2) -       &
!        & ocean_nudge%forc_3dimrelax_salt(:,:,:))

      DO blockNo = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, blockNo, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO jc = start_index, end_index
          DO level = 1, dolic_c(jc,blockNo)
            ocean_nudge%forc_3dimrelax_salt(jc,level,blockNo) = z_relax * dtime &
              * (ocean_nudge%data_3dimrelax_salt(jc,level,blockNo) - p_os%p_prog(nnew(1))%tracer(jc,level,blockNo,2) )

            ! add relaxation term to new temperature
            p_os%p_prog(nnew(1))%tracer(jc,level,blockNo,2) = p_os%p_prog(nnew(1))%tracer(jc,level,blockNo,2) + &
                                                  ocean_nudge%forc_3dimRelax_salt(jc,level,blockNo)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
      !$ACC WAIT(1)

      ! add relaxation term to new salinity

      !---------DEBUG DIAGNOSTICS-------------------------------------------
      idt_src=1  ! output print level (1-5, fix)
      CALL dbg_print('3d_relax: tracer S forc'  ,ocean_nudge%forc_3dimrelax_salt,str_module,idt_src, in_subset=cells_in_domain)
      CALL dbg_print('3d_relax: tracer S data'  ,ocean_nudge%data_3dimrelax_salt,str_module,idt_src, in_subset=cells_in_domain)
      z_c(:,:,:) =  p_os%p_prog(nnew(1))%tracer(:,:,:,2)
      CALL dbg_print('3d_relax: tracer S trac'  ,z_c                           ,str_module,idt_src, in_subset=cells_in_domain)
      !---------------------------------------------------------------------

    END IF

  END SUBROUTINE nudge_ocean_tracers
  !-------------------------------------------------------------------------


END MODULE mo_ocean_nudging


