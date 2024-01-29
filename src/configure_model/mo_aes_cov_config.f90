! Configuration of the parameterization for cloud cover,
! that is used in the AES physics package.
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

MODULE mo_aes_cov_config

  USE mo_exception            ,ONLY: finish, message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_cov_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_cov_config   !< allocate and initialize aes_cov_config
  PUBLIC ::    eval_aes_cov_config   !< evaluate aes_cov_config
  PUBLIC ::   print_aes_cov_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_cov'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the AES microphysics
  !!
  TYPE t_aes_cov_config
     !
     ! configuration parameters
     ! ------------------------
     !
     REAL(wp) :: cqx      !          critical mass fraction of cloud water + cloud ice in air, in kg/kg
     !
  END TYPE t_aes_cov_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_cov_config), TARGET :: aes_cov_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_cov_config
    !
    ! AES cloud cover configuration
    ! -------------------------------
    !
    aes_cov_config(:)% cqx      = 1.0e-8_wp
    !
  END SUBROUTINE init_aes_cov_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_aes_cov_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    DO jg = 1,ng
       !
       IF (aes_cov_config(jg)% cqx < 0.0_wp .OR. 1.0_wp < aes_cov_config(jg)% cqx) THEN
          WRITE(cg,'(i0)') jg
          CALL finish('eval_aes_cov_config', &
               &      'aes_cov_config('//TRIM(cg)//')% cqx <0 or >1 is not allowed')
       END IF
       !
    END DO
    !
  END SUBROUTINE eval_aes_cov_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_cov_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','AES cloud cover configuration')
    CALL message    ('','=============================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL message    ('','---      0/1 cloud cover based on mass fraction of cloud water + cloud ice in air')
       CALL print_value('    aes_cov_config('//TRIM(cg)//')% cqx      ',aes_cov_config(jg)% cqx     )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_aes_cov_config

  !----

END MODULE mo_aes_cov_config
