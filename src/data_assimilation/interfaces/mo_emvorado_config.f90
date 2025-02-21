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

MODULE mo_emvorado_config
  
  USE mo_run_config,      ONLY: msg_level
  USE mo_exception,       ONLY: message






  IMPLICIT NONE

  PUBLIC :: config_emvorado

  LOGICAL :: radar_is_initialized = .FALSE.


CONTAINS


  SUBROUTINE config_emvorado (n_dom_model, radar_flag_doms_model)

    INTEGER, INTENT(in) :: n_dom_model                            ! Number of model domains
    LOGICAL, INTENT(in) :: radar_flag_doms_model (1:n_dom_model)  ! Switch for running EMVORADO for each domain

    INTEGER             :: jg
    character(len=3)    :: cjg

!!$ timer_start(timer_radar_ini) is done in organize_radar('init')


!!$ timer_stop(timer_radar_ini) is done in organize_radar('init')
    
  END SUBROUTINE config_emvorado
  

END MODULE mo_emvorado_config
