! Configuration of the Carbon cycle settings
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

MODULE mo_ccycle_config

  USE mo_exception     ,ONLY: message, print_value, finish
  USE mo_kind          ,ONLY: wp
  USE mo_impl_constants,ONLY: max_dom
  USE mo_grid_config   ,ONLY: n_dom

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                  name       !< name for this unit

  ! configuration
  PUBLIC ::       t_ccycle_config
  PUBLIC ::         ccycle_config
  PUBLIC ::    init_ccycle_config       !< initialize ccycle_config
  PUBLIC ::   print_ccycle_config       !< print out

  ! Named constants
  PUBLIC :: CCYCLE_MODE_NONE
  PUBLIC :: CCYCLE_MODE_INTERACTIVE
  PUBLIC :: CCYCLE_MODE_PRESCRIBED

  PUBLIC :: CCYCLE_CO2CONC_CONST
  PUBLIC :: CCYCLE_CO2CONC_FROMFILE

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'ccycle'

  !> No carbon cycle.
  INTEGER, PARAMETER :: CCYCLE_MODE_NONE = 0
  !> Atmospheric CO2 concentration fully interacts with land-surface scheme.
  INTEGER, PARAMETER :: CCYCLE_MODE_INTERACTIVE = 1
  !> Atmospheric CO2 concentration is prescribed.
  INTEGER, PARAMETER :: CCYCLE_MODE_PRESCRIBED = 2

  !> Constant atmospheric CO2 concentration (for CCYCLE_MODE_PRESCRIBED).
  INTEGER, PARAMETER :: CCYCLE_CO2CONC_CONST = 2
  !> CO2 volume mixing ratio is read from scenario file `bc_greenhouse_gases.nc` (for
  !! CCYCLE_MODE_PRESCRIBED).
  INTEGER, PARAMETER :: CCYCLE_CO2CONC_FROMFILE = 4

  !>
  !! Configuration type containing switches for the configuration of the carbon cycle
  !!
  TYPE t_ccycle_config
     !
     ! configuration parameters
     ! ------------------------
     !
     INTEGER  :: iccycle   !< c-cycle mode
     INTEGER  :: ico2conc  !< co2 concentration provided to land and ocean
     !
     REAL(wp) :: vmr_co2   !< co2 volume mixing ratio for c-cycle
     !
  END TYPE t_ccycle_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_ccycle_config) :: ccycle_config(max_dom)

CONTAINS

  !----
  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_ccycle_config
    !
    ! Carbon cycle configuration
    ! --------------------------
    !
    ccycle_config(:)% iccycle  = CCYCLE_MODE_NONE
    !                                         ! 0: no c-cycle
    !                                           1: c-cycle with interactive atm. co2 concentration
    !                                           2: c-cycle with prescribed  atm. co2 concentration
    !
    ! For iccycle = 2:
    ccycle_config(:)% ico2conc = CCYCLE_CO2CONC_CONST
    !                                         ! 2: constant  co2 concentration vmr_co2
    !                                           4: transient co2 concentration scenario from file
    !
    ! For ico2conc = 2:
    ccycle_config(:)% vmr_co2  = 284.3e-06_wp ! co2 volume mixing ratio of 1850 (CMIP6)
    !
  END SUBROUTINE init_ccycle_config

  !----

  !>
  !! Print out the user controlled configuration
  !!
  SUBROUTINE print_ccycle_config
    !
    INTEGER           :: jg
    CHARACTER(LEN=2)  :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','Carbon cycle configuration')
    CALL message    ('','==========================')
    CALL message    ('','')
    !
    DO jg = 1,n_dom
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    ccycle_config('//TRIM(cg)//')% iccycle  ',ccycle_config(jg)% iccycle )
       CALL print_value('    ccycle_config('//TRIM(cg)//')% ico2conc ',ccycle_config(jg)% ico2conc)
       CALL print_value('    ccycle_config('//TRIM(cg)//')% vmr_co2  ',ccycle_config(jg)% vmr_co2 )
       CALL message    ('','')
       !
       SELECT CASE(ccycle_config(jg)% iccycle)
       CASE (CCYCLE_MODE_NONE)
          CALL message ('','C-cycle is switched off')
       CASE (CCYCLE_MODE_INTERACTIVE)
          CALL message ('','C-cycle is used with interactive atmospheric CO2 concentration')
       CASE (CCYCLE_MODE_PRESCRIBED)
          CALL message ('','C-cycle is used with prescribed atmospheric CO2 concentration')
          !
          SELECT CASE(ccycle_config(jg)% ico2conc)
          CASE (CCYCLE_CO2CONC_CONST)
             CALL print_value('    CO2 volume mixing ratio is constant (ppv)',ccycle_config(jg)% vmr_co2)
          CASE (CCYCLE_CO2CONC_FROMFILE)
             CALL message('','CO2 volume mixing ratio is read from scenario file bc_greenhouse_gases.nc')
          CASE default
             CALL finish('print_ccycle_config','ccycle_config(jg)% ico2conc invalid, must be 2 or 4')
          END SELECT
          !
       CASE default
          CALL finish('print_ccycle_config','ccycle_config(jg)% iccycle invalid, must be 0, 1 or 2')
       END SELECT
       !
    END DO
    !
  END SUBROUTINE print_ccycle_config

  !----

END MODULE mo_ccycle_config
