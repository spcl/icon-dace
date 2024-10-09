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

! Configuration of the upper and lower levels for WMO_tropopause calculation
! that is used in the AES physics package.

MODULE mo_aes_wmo_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom

  USE mo_vertical_coord_table ,ONLY: vct_a

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_wmo_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_wmo_config   !< allocate and initialize aes_wmo_config
  PUBLIC ::    eval_aes_wmo_config   !< evaluate aes_wmo_config
  PUBLIC ::   print_aes_wmo_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_wmo'
  
  !>
  !! Configuration type containing parameters and switches for the configuration of the AES microphysics
  !!
  TYPE t_aes_wmo_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! vertical loop limits
     REAL(wp) :: zmaxwmo  !          maximum height (m) for tropopause calculation
     REAL(wp) :: zminwmo  !          minimum height (m) for tropopause calculation
     INTEGER  :: jkswmo   !          vertical start index for tropopause calculation
     INTEGER  :: jkewmo   !          vertical end   index for tropopause calculation
     !                               diagnosed in eval_aes_wmo_config
  END TYPE t_aes_wmo_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_wmo_config), TARGET :: aes_wmo_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_wmo_config
    !
    ! WMO_tropopause configuration
    ! ----------------------------
    !
    aes_wmo_config(:)% zmaxwmo  = 38000.
    aes_wmo_config(:)% zminwmo  =  5000.
    !
  END SUBROUTINE init_aes_wmo_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_aes_wmo_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg, jk, klev
    !
    klev = SIZE(vct_a)-1
    !
    DO jg = 1,ng
       !
       ! diagnose jkswmo = start index
       !
       aes_wmo_config(jg)% jkswmo = 1
       DO jk = 1,klev
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > aes_wmo_config(jg)% zmaxwmo) THEN
             aes_wmo_config(jg)% jkswmo = jk + 1
          ELSE
             EXIT
          END IF
       END DO
       !
       ! diagnose jkewmo = end index
       !
       aes_wmo_config(jg)% jkewmo = aes_wmo_config(jg)% jkswmo
       DO jk = aes_wmo_config(jg)% jkswmo,klev
          IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp > aes_wmo_config(jg)% zminwmo) THEN
             aes_wmo_config(jg)% jkewmo = jk + 1
          ELSE
             EXIT
          END IF
       END DO
       !
    END DO
    !
  END SUBROUTINE eval_aes_wmo_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_wmo_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','WMO tropopause configuration')
    CALL message    ('','============================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    aes_wmo_config('//TRIM(cg)//')% zmaxwmo ',aes_wmo_config(jg)% zmaxwmo )
       CALL print_value('    aes_wmo_config('//TRIM(cg)//')% zminwmo ',aes_wmo_config(jg)% zminwmo )
       CALL print_value('    aes_wmo_config('//TRIM(cg)//')% jkswmo  ',aes_wmo_config(jg)% jkswmo  )
       CALL print_value('    aes_wmo_config('//TRIM(cg)//')% jkewmo  ',aes_wmo_config(jg)% jkewmo    )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_aes_wmo_config

  !----

END MODULE mo_aes_wmo_config
