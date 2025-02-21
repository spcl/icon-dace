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


    CALL finish('nwp_couple_hydrodisc: unintentionally called. Check your source code and configure.')

  END SUBROUTINE nwp_couple_hydrodisc

END MODULE mo_nwp_hydrodisc_coupling
