! Configuration of the parameterization for cloud optical properties,
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

MODULE mo_aes_cop_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_physical_constants   ,ONLY: tmelt
  USE mtime,                   ONLY: OPERATOR(>)

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_cop_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_cop_config   !< allocate and initialize aes_cop_config
!!$  PUBLIC ::    eval_aes_cop_config   !< evaluate aes_cop_config
  PUBLIC ::   print_aes_cop_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_cop'
  
  !>
  !! Configuration type containing parameters for the configuration of the cloud optical properties
  !!  and parameters of aes cloud microphysics, still used in mo_aes_convect_tables
  !!
  TYPE t_aes_cop_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! cloud droplet number concentration
     REAL(wp) :: cn1lnd   ! [1e6/m3] over land, p <= 100 hPa
     REAL(wp) :: cn2lnd   ! [1e6/m3] over land, p >= 800 hPa
     REAL(wp) :: cn1sea   ! [1e6/m3] over sea , p <= 100 hPa
     REAL(wp) :: cn2sea   ! [1e6/m3] over sea , p >= 800 hPa
     !
     ! cloud optics
     REAL(wp) :: cinhomi  !          ice    cloud inhomogeneity factor
     REAL(wp) :: cinhoms  !          snow   cloud inhomogeneity factor
     REAL(wp) :: cinhoml  !          liquid cloud inhomogeneity factor
     !
     ! freezing/deposition/sublimation for mo_aes_convect_tables:
     REAL(wp) :: cthomi   ! [K]      maximum temperature for homogeneous freezing
     REAL(wp) :: csecfrl  ! [kg/kg]  minimum in-cloud water mass mixing ratio in mixed phase clouds
     !
  END TYPE t_aes_cop_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_cop_config), TARGET :: aes_cop_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_cop_config

    !
    ! aes cloud microphysics active
    !
    ! cloud optical properties configuration
    ! --------------------------------------
    !
    ! cloud droplet number concentration
    aes_cop_config(:)% cn1lnd   =  20._wp
    aes_cop_config(:)% cn2lnd   = 180._wp
    aes_cop_config(:)% cn1sea   =  20._wp
    aes_cop_config(:)% cn2sea   =  80._wp
    !
    ! cloud optics
    aes_cop_config(:)% cinhomi  = 0.80_wp
    aes_cop_config(:)% cinhoms  = 0.80_wp
    aes_cop_config(:)% cinhoml  = 0.80_wp
    !
    ! freezing/deposition/sublimation
    aes_cop_config(:)% cthomi   = tmelt-35.0_wp
    aes_cop_config(:)% csecfrl  = 1.5e-5_wp
    !
  END SUBROUTINE init_aes_cop_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_aes_cop_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_aes_cop_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_cop_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','cloud optical properties configuration')
    CALL message    ('','======================================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn1lnd   ',aes_cop_config(jg)% cn1lnd  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn2lnd   ',aes_cop_config(jg)% cn2lnd  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn1sea   ',aes_cop_config(jg)% cn1sea  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn2sea   ',aes_cop_config(jg)% cn2sea  )
       CALL message    ('','')
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhomi  ',aes_cop_config(jg)% cinhomi )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoms  ',aes_cop_config(jg)% cinhoms )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml  ',aes_cop_config(jg)% cinhoml )
       CALL message    ('','')
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cthomi   ',aes_cop_config(jg)% cthomi  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% csecfrl  ',aes_cop_config(jg)% csecfrl )
       CALL message    ('','')

       !
    END DO
    !
  END SUBROUTINE print_aes_cop_config

  !----

END MODULE mo_aes_cop_config
