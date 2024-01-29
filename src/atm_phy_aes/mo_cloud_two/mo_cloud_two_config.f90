!
! Configuration of the parameterization for two-moment bulk microphysics by
! Seifert and Beheng (2006) with prognostic cloud droplet number from DWD
! (inwp_gscp=4) to be used in the sapphire physics package.
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

MODULE mo_cloud_two_config

  USE mo_exception            ,ONLY: message, print_value
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom

  USE mo_cloud_two_types      ,ONLY: t_cloud_two_config

  IMPLICIT NONE
  PRIVATE

  ! configuration
  PUBLIC ::         cloud_two_config   !< user specified configuration parameters
  PUBLIC ::    init_cloud_two_config   !< allocate and initialize cloud_two_config
!!$  PUBLIC ::    eval_cloud_two_config   !< evaluate cloud_two_config
  PUBLIC ::   print_cloud_two_config   !< print out

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_cloud_two_config), TARGET :: cloud_two_config(max_dom)
  
CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_cloud_two_config
    !
    ! two-moment bulk microphysics configuration
    ! ------------------------------------------
    !
    ! no parameter settings yet available...
    !
    ! general thresholds
    !
    ! grid scale microphysics
    !
  END SUBROUTINE init_cloud_two_config

  !----

!!$  !>
!!$  !! Evaluate additional derived parameters
!!$  !!
!!$  SUBROUTINE eval_cloud_two_config
!!$    !
!!$    ...
!!$    !
!!$  END SUBROUTINE eval_cloud_two_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_cloud_two_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('',' Two-moment bulk microphysics configuration')
    CALL message    ('','=============================0=============')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       !
       ! kept as example... nothing yet available...
       !
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% qi0            ',cloud_two_config(jg)% qi0            )
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% qc0            ',cloud_two_config(jg)% qc0            )
       !CALL message    ('','')
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% zceff_min      ',cloud_two_config(jg)% zceff_min      )
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% v0snow         ',cloud_two_config(jg)% v0snow         )
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% zvz0i          ',cloud_two_config(jg)% zvz0i          )
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% icesedi_exp    ',cloud_two_config(jg)% icesedi_exp    )
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% mu_rain        ',cloud_two_config(jg)% mu_rain        )
       !CALL print_value('    cloud_two_config('//TRIM(cg)//')% rain_n0_factor ',cloud_two_config(jg)% rain_n0_factor )
       CALL message    ('','')
       !
    END DO
    !
  END SUBROUTINE print_cloud_two_config

  !----

END MODULE mo_cloud_two_config
