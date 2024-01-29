! Interface between NWP physics and the hydrological discharge model, through a coupler.
! Based on the ocean coupling interface
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
!
!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_hydrodisc_coupling

  USE mo_kind                ,ONLY: wp
  USE mo_model_domain        ,ONLY: t_patch
  USE mo_nwp_lnd_types       ,ONLY: t_lnd_diag
  USE mo_nwp_phy_types       ,ONLY: t_nwp_phy_diag
  USE mo_lnd_nwp_config      ,ONLY: ntiles_total, isub_lake
  USE mo_ext_data_types      ,ONLY: t_external_data
  USE mo_fortran_tools       ,ONLY: init
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_atm_phy_nwp_config  ,ONLY: atm_phy_nwp_config
  USE mo_impl_constants      ,ONLY: min_rlcell, LSS_TERRA, SUCCESS
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_run_config          ,ONLY: ltimer, dtime
  USE mo_timer               ,ONLY: timer_start, timer_stop, timer_coupling_put

  USE mo_coupling_config     ,ONLY: is_coupled_to_hydrodisc
#ifdef YAC_coupling
  USE mo_atmo_coupling_frame ,ONLY: field_id, CPF_RUNOFFS, CPF_RUNOFFG
  USE mo_yac_finterface      ,ONLY: yac_fput, yac_dble_ptr, YAC_ACTION_COUPLING, YAC_ACTION_OUT_OF_BOUND
#endif

  USE mo_exception           ,ONLY: warning, message, finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nwp_couple_hydrodisc

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_nwp_hydrodisc_coupling' ! Output of module for debug

CONTAINS


  !>
  !! SUBROUTINE nwp_couple_hydrodisc -- the interface between
  !! NWP physics and the hydrodisc, through a coupler
  !!
  !! This subroutine is called from nwp_nh_interface.

  SUBROUTINE nwp_couple_hydrodisc( p_patch, lnd_diag, prm_diag, ext_data )

    ! Arguments

    TYPE(t_patch),   TARGET, INTENT(INOUT)  :: p_patch
    TYPE(t_lnd_diag),        INTENT(INOUT)  :: lnd_diag
    TYPE(t_nwp_phy_diag),    INTENT(INOUT)  :: prm_diag
    TYPE(t_external_data),   INTENT(INOUT)  :: ext_data

    ! Local variables

    LOGICAL               :: write_coupler_restart
    INTEGER               :: nbr_hor_cells         ! inner points
    INTEGER               :: jg                    ! patch ID
    INTEGER               :: jb                    ! block loop count
    INTEGER               :: jc                    ! nproma loop count
    INTEGER               :: error                 
    INTEGER               :: info, ierror          ! return values from cpl_put/get calls
    INTEGER               :: rl_start, rl_end
    INTEGER               :: i_startblk, i_endblk  ! blocks
    INTEGER               :: i_startidx, i_endidx  ! slices
    INTEGER               :: isubs                 ! tile index
    REAL(wp), TARGET, ALLOCATABLE :: buffer(:,:)   ! buffer transferred to YAC coupler
    CHARACTER(LEN=*), PARAMETER   :: routine = 'nwp_couple_hydrodisc'

#ifdef YAC_coupling
    TYPE(yac_dble_ptr)    :: ptr(1,1)
#endif

#ifndef YAC_coupling
    CALL finish('nwp_couple_hydrodisc: unintentionally called. Check your source code and configure.')
#else

    ALLOCATE(buffer(nproma, p_patch%nblks_c), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    jg          = p_patch%id

    ! include boundary interpolation zone of nested domains and halo points
    rl_start = 1
    rl_end   = min_rlcell

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    nbr_hor_cells = p_patch%n_patch_cells

    !-------------------------------------------------------------------------
    ! If running in atm-hydrological discharge coupled mode, exchange information 
    !-------------------------------------------------------------------------
    !
    ! Possible fields that contain information to be sent to the hydrological
    ! discharge model include
    !
    ! 1. lnd_diag%runoff_s_inst_t(:,:,:)    averaged surface water runoff   [kg/m2/s] 
    ! 2. lnd_diag%runoff_g_inst_t(:,:,:)    averaged soil water runoff      [kg/m2/s] 
    !
    !-------------------------------------------------------------------------

    !------------------------------------------------
    !  Send surface water runoff
    !    field_id(CPF_RUNOFFS) represents "surface water runoff"
    !------------------------------------------------

    !$OMP PARALLEL
    CALL init(buffer(:,:))
    !$OMP END PARALLEL

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, isubs) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)
      IF ( atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN

        ! aggregate over tiles

        DO isubs = 1, ntiles_total
          IF ( isubs == isub_lake ) THEN
            DO jc = i_startidx, i_endidx
              ! Take P-E over the lake as runoff
              buffer(jc,jb) = buffer(jc,jb) + ( prm_diag%tot_prec_rate(jc,jb) + prm_diag%qhfl_s_t(jc,jb,isubs) )  &
                &           * ext_data%atm%frac_t(jc,jb,isubs)
            ENDDO
          ELSE
            DO jc = i_startidx, i_endidx
              buffer(jc,jb) = buffer(jc,jb) + lnd_diag%runoff_s_inst_t(jc,jb,isubs) / dtime                       &
                &           * ext_data%atm%frac_t(jc,jb,isubs)
            ENDDO
          ENDIF
        ENDDO  ! isubs
        
      ELSE ! JSBACH

        DO jc = i_startidx, i_endidx
          buffer(jc,jb) = lnd_diag%runoff_s_inst_t(jc,jb,1) / dtime
        ENDDO

      ENDIF
    ENDDO ! jb

!ICON_OMP_END_PARALLEL_DO

    ptr(1,1)%p(1:nbr_hor_cells) => buffer

    IF (ltimer) CALL timer_start(timer_coupling_put)
    CALL yac_fput ( field_id(CPF_RUNOFFS), SIZE(ptr(:,1:1),1), 1, ptr(:,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_hydrodisc', 'YAC says fput called after end of run - id=1, surface water runoff')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
      CALL message('nwp_couple_hydrodisc', 'YAC says it is put for restart - ids 1, surface water runoff')
    ENDIF

    !------------------------------------------------
    !  Send soil water runoff
    !    field_id(CPF_RUNOFFG) represents "soil water runoff"
    !------------------------------------------------

    !$OMP PARALLEL
    CALL init(buffer(:,:))
    !$OMP END PARALLEL

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, isubs) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, rl_start, rl_end)

      IF ( atm_phy_nwp_config(jg)%inwp_surface == LSS_TERRA ) THEN

        ! aggregate over tiles
        DO isubs = 1, ntiles_total
          DO jc = i_startidx, i_endidx
            buffer (jc,jb) = buffer(jc,jb) + lnd_diag%runoff_g_inst_t(jc,jb,isubs) / dtime  &
                &            * ext_data%atm%frac_t(jc,jb,isubs)
          ENDDO
        ENDDO  ! isubs
        
      ELSE ! JSBACH

        DO jc = i_startidx, i_endidx
          buffer(jc,jb) = lnd_diag%runoff_g_inst_t(jc,jb,1) / dtime
        ENDDO

      ENDIF
    ENDDO ! jb
!ICON_OMP_END_PARALLEL_DO

    ptr(1,1)%p(1:nbr_hor_cells) => buffer

    IF (ltimer) CALL timer_start(timer_coupling_put)
    CALL yac_fput ( field_id(CPF_RUNOFFG), SIZE(ptr(:,1:1),1), 1, ptr(:,1:1), info, ierror )
    IF ( info > YAC_ACTION_COUPLING .AND. info < YAC_ACTION_OUT_OF_BOUND ) THEN
      write_coupler_restart = .TRUE.
    ELSE
      write_coupler_restart = .FALSE.
    ENDIF

    IF ( info == YAC_ACTION_OUT_OF_BOUND ) &
         & CALL warning('nwp_couple_hydrodisc', 'YAC says fput called after end of run - id=2, ground water runoff')

    IF (ltimer) CALL timer_stop(timer_coupling_put)

    IF ( write_coupler_restart ) THEN
      CALL message('nwp_couple_hydrodisc', 'YAC says it is put for restart - ids 2, ground water runoff')
    ENDIF

    DEALLOCATE(buffer, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

!YAC_coupling
#endif

  END SUBROUTINE nwp_couple_hydrodisc

END MODULE mo_nwp_hydrodisc_coupling
