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

MODULE mo_emvorado_init







  IMPLICIT NONE

  PUBLIC ::  init_emvorado_mpi, prep_emvorado_domains

CONTAINS

  SUBROUTINE prep_emvorado_domains (n_dom_model, radar_flag_doms_model)
     
    INTEGER, INTENT(in) :: n_dom_model                            ! Number of model domains
    LOGICAL, INTENT(in) :: radar_flag_doms_model (1:n_dom_model)  ! Switch for running EMVORADO for each domain
    





  END SUBROUTINE prep_emvorado_domains
  
  SUBROUTINE init_emvorado_mpi (luse_radarfwo,                                              & ! INPUT
                                comm_world_icon, my_world_id_icon, nproc_icon,              & ! INPUT
                                comm_work_icon, my_work_id_icon, num_work_icon,             & ! INPUT
                                lwork_pe_icon,                                              & ! INPUT
                                nprocio_radar_icon, radar_master_icon, radario_master_icon, & ! INPUT
                                ierror, errmsg                                              & ! OUTPUT
                                )
  
    ! INPUT parameters:
    !------------------
    INTEGER, INTENT(in) :: comm_work_icon , my_work_id_icon , num_work_icon
    INTEGER, INTENT(in) :: comm_world_icon, my_world_id_icon, &
                           nproc_icon

    INTEGER, INTENT(in) :: radar_master_icon,   & ! Start-PE of comm_radar in comm_world_icon
                           radario_master_icon, & ! Start-PE of comm_radario in comm_world_icon
                           nprocio_radar_icon     ! Number of asynchronous radar IO PEs provided by ICON

    LOGICAL, INTENT(in) :: luse_radarfwo(:), lwork_pe_icon

    INTEGER,          INTENT(out)   :: ierror
    CHARACTER(len=*), INTENT(inout) :: errmsg

    ierror = 0

    
  END SUBROUTINE init_emvorado_mpi


END MODULE mo_emvorado_init
