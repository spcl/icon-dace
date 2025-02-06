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

#ifdef HAVE_RADARFWO
  USE radar_data_namelist, ONLY: crosscheck_domains_radar_nml
  USE radar_organize,      ONLY: organize_radar
#endif

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

#ifdef HAVE_RADARFWO
    IF ( .NOT. radar_is_initialized ) THEN

     ! First call of radar operator on this PE, therefore a call to radar station config ("initialization") is necessary:
      DO jg=1, n_dom_model
        IF (radar_flag_doms_model(jg)) THEN
          IF (msg_level >= 11) THEN
            cjg(:) = ' '
            WRITE (cjg, '(i3)') jg
            CALL message('config_emvorado(): ', &
                 'initializing EMVORADO on domain '//TRIM(ADJUSTL(cjg))//' (namelist RADARSIM_PARAMS + general setup)')
          END IF
          CALL organize_radar('init', 1, 1, jg)   ! the time levels for dyn and qx are not relevant here, "1" are only dummies
        END IF
      END DO
      CALL crosscheck_domains_radar_nml ()
      radar_is_initialized = .TRUE.
    END IF
#endif

!!$ timer_stop(timer_radar_ini) is done in organize_radar('init')
    
  END SUBROUTINE config_emvorado
  

END MODULE mo_emvorado_config
