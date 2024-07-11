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
  USE mo_ext_data_types      ,ONLY: t_external_data
  USE mo_lnd_nwp_config      ,ONLY: ntiles_total, ntiles_water
  USE mo_fortran_tools       ,ONLY: init, assert_acc_host_only
  USE mo_parallel_config     ,ONLY: nproma
  USE mo_atm_phy_nwp_config  ,ONLY: atm_phy_nwp_config
  USE mo_impl_constants      ,ONLY: start_prog_cells, end_prog_cells, LSS_TERRA, SUCCESS
  USE mo_physical_constants  ,ONLY: rhoh2o
  USE mo_run_config          ,ONLY: dtime, msg_level
  USE mo_loopindices         ,ONLY: get_indices_c

  USE mo_coupling_utils      ,ONLY: cpl_def_field, cpl_put_field

  USE mo_exception           ,ONLY: finish, message, message_text
  USE mo_sync                ,ONLY: global_sum_array

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_nwp_hydrodisc_coupling, nwp_couple_hydrodisc

  CHARACTER(len=*), PARAMETER :: str_module = 'mo_nwp_hydrodisc_coupling' ! Output of module for debug

  INTEGER :: field_id_runoffs, field_id_runoffg

CONTAINS

  !>
  !! Registers fields required for the coupling between NWP physics and
  !! the hydrological discharge model
  !!
  !! This subroutine is called from construct_atmo_coupling.
  !!
  SUBROUTINE construct_nwp_hydrodisc_coupling( &
    p_patch, ext_data, comp_id, grid_id, cell_point_id, timestepstring)

    TYPE(t_patch), TARGET, INTENT(IN) :: p_patch(:)
    TYPE(t_external_data), INTENT(IN) :: ext_data(:)
    INTEGER, INTENT(IN) :: comp_id
    INTEGER, INTENT(IN) :: grid_id
    INTEGER, INTENT(IN) :: cell_point_id
    CHARACTER(LEN=*), INTENT(IN) :: timestepstring

    TYPE(t_patch), POINTER :: patch_horz

    INTEGER :: jg

    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':construct_nwp_hydrodisc_coupling'

    jg = 1
    patch_horz => p_patch(jg)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "surface_water_runoff", 1, field_id_runoffs)

    CALL cpl_def_field( &
      comp_id, cell_point_id, timestepstring, &
      "soil_water_runoff", 1, field_id_runoffg)

  END SUBROUTINE construct_nwp_hydrodisc_coupling

  !>
  !! SUBROUTINE nwp_couple_hydrodisc -- the interface between
  !! NWP physics and the hydrodisc, through a coupler
  !!
  !! This subroutine is called from nwp_nh_interface.

  SUBROUTINE nwp_couple_hydrodisc( p_patch, lnd_diag, prm_diag, ext_data, lacc )

    ! Arguments

    TYPE(t_patch),   TARGET, INTENT(INOUT)  :: p_patch
    TYPE(t_lnd_diag),        INTENT(INOUT)  :: lnd_diag
    TYPE(t_nwp_phy_diag),    INTENT(INOUT)  :: prm_diag
    TYPE(t_external_data),   INTENT(INOUT)  :: ext_data
    LOGICAL, OPTIONAL,       INTENT(IN)     :: lacc

    ! Local variables

    INTEGER               :: jg                    ! patch ID
    INTEGER               :: jb                    ! block loop count
    INTEGER               :: jc                    ! nproma loop count
    INTEGER               :: error
    INTEGER               :: i_startblk, i_endblk  ! blocks
    INTEGER               :: i_startidx, i_endidx  ! slices
    INTEGER               :: isubs                 ! tile index
    REAL(wp), TARGET, ALLOCATABLE :: buffer(:,:,:) ! buffer transferred to YAC coupler
    CHARACTER(LEN=*), PARAMETER   :: routine = str_module // ':nwp_couple_hydrodisc'
    REAL(wp)              :: diag_tmp

    ! This routine hasn't been ported yet.
    CALL assert_acc_host_only(routine, lacc)

    ALLOCATE(buffer(nproma, p_patch%nblks_c, 2), STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "memory allocation failure")

    jg         = p_patch%id

    i_startblk = p_patch%cells%start_block(start_prog_cells)
    i_endblk   = p_patch%cells%end_block(end_prog_cells)

    !-------------------------------------------------------------------------
    ! If running in atm-hydrological discharge coupled mode, exchange information
    !-------------------------------------------------------------------------
    !
    ! Possible fields that contain information to be sent to the hydrological
    ! discharge model include
    !
    ! 1. lnd_diag%runoff_s_inst_t(:,:,:)  tile-specific instantaneous surface water runoff [kg/m2]
    ! 2. lnd_diag%runoff_g_inst_t(:,:,:)  tile-specific instantaneous ground water runoff  [kg/m2]
    !
    !-------------------------------------------------------------------------

    !$OMP PARALLEL
    CALL init(buffer(:,:,:))
    !$OMP END PARALLEL

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, i_startidx, i_endidx, isubs) ICON_OMP_RUNTIME_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                       & i_startidx, i_endidx, start_prog_cells, end_prog_cells)

      ! aggregate over tiles
      DO isubs = 1, ntiles_total + ntiles_water
        DO jc = i_startidx, i_endidx
          buffer(jc,jb,1) = buffer(jc,jb,1) + lnd_diag%runoff_s_inst_t(jc,jb,isubs) / dtime &
            &           * ext_data%atm%frac_t(jc,jb,isubs)
          buffer(jc,jb,2) = buffer(jc,jb,2) + lnd_diag%runoff_g_inst_t(jc,jb,isubs) / dtime &
            &           * ext_data%atm%frac_t(jc,jb,isubs)
        ENDDO
      ENDDO  ! isubs
    ENDDO ! jb

!ICON_OMP_END_PARALLEL_DO

    ! Online diagnose for global sum of runoff (m3/s) before sending to YAC:
    IF (msg_level >= 10) THEN
      diag_tmp = global_sum_array(buffer(:,:,1) * p_patch%cells%area(:,:) / rhoh2o)
      WRITE(message_text,'(a,f15.3)') ' NWP-HD: global total surface runoff (m3/s) :' , diag_tmp
      CALL message(routine, message_text)

      diag_tmp = global_sum_array(buffer(:,:,2) * p_patch%cells%area(:,:) / rhoh2o)
      WRITE(message_text,'(a,f15.3)') ' NWP-HD: global total ground runoff (m3/s) :' , diag_tmp
      CALL message(routine, message_text)
    ENDIF

    !------------------------------------------------
    !  Send surface water runoff
    !    "surface water runoff"
    !------------------------------------------------
    CALL cpl_put_field( &
      routine, field_id_runoffs, 'surface water runoff', &
      p_patch%n_patch_cells, buffer(:,:,1))

    !------------------------------------------------
    !  Send soil water runoff
    !    "soil water runoff"
    !------------------------------------------------
    CALL cpl_put_field( &
      routine, field_id_runoffg, 'ground water runoff', &
      p_patch%n_patch_cells, buffer(:,:,2))

    DEALLOCATE(buffer, STAT = error)
    IF(error /= SUCCESS) CALL finish(routine, "DEALLOCATE failed!")

  END SUBROUTINE nwp_couple_hydrodisc

END MODULE mo_nwp_hydrodisc_coupling
