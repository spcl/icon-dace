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

!> @brief
!! Module defining data types for atmosphere energy diagnostics.
!!
!! Author: Marco Giorgetta, MPI-M, 2024

MODULE mo_atm_energy_types

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_atm_energy_config, t_atm_energy

  !-----------------------------------------------------------------------------

  TYPE t_atm_energy_config
     !
     ! configuration parameters
     ! ------------------------
     !
     LOGICAL :: l_atm_energy  !< main switch: .TRUE. = on, .FALSE. = off
     !
  END TYPE t_atm_energy_config

  !-----------------------------------------------------------------------------

  TYPE t_atm_energy
     !
     ! atmosphere energy state
     ! =======================
     !
     ! - kind of energy:
     !
     !   ein : moist internal
     !   ekh : horizontal kinetic
     !   ekv : vertical kinetic
     !   egp : geopotential
     !   eto : total = ein + ekh + ekv + egp
     !
     ! - at different stages of the time step calculations:
     !
     !   ..1.. : 1st stage, for difference ..2.. - ..1..,
     !           used also for output of energies diagnosed after the physics state is updated from the last parameterization
     !           -> output variables 'exphy..'
     !
     !   ..2.. : 2nd stage, for difference ..2.. - ..1..,
     !           used also for output of energies diagnosed after the dynamics state is updated from physics, at the end of the time step
     !           -> output variables 'ex..'
     !
     !   ..3.. : 3rd stage, for difference ..2.. - ..3..,
     !           used also for output of energies diagnosed after the dynamics state is updated from dynamics
     !           -> output variables 'exdyn..'
     !
     ! - spatial context:
     !
     !   ..   : [J/m3] density
     !   ..vi : [J/m2] content   = vertical integral of density
     !   ..hi : [J/m]  profile   = layerwise horizontal integral of density
     !   ..ti : [J]    total     = total integral of density
     !
     REAL(wp),POINTER ::                                                                            &
          !
          ! - moist internal energy
          !
          & ein1(:,:,:) => NULL(), ein1vi(:,:) => NULL(), ein1hi(:) => NULL(), ein1ti(:) => NULL(), &
          & ein2(:,:,:) => NULL(), ein2vi(:,:) => NULL(), ein2hi(:) => NULL(), ein2ti(:) => NULL(), &
          & ein3(:,:,:) => NULL(), ein3vi(:,:) => NULL(), ein3hi(:) => NULL(), ein3ti(:) => NULL(), &
          !
          ! - horizontal kinetic energy
          !
          & ekh1(:,:,:) => NULL(), ekh1vi(:,:) => NULL(), ekh1hi(:) => NULL(), ekh1ti(:) => NULL(), &
          & ekh2(:,:,:) => NULL(), ekh2vi(:,:) => NULL(), ekh2hi(:) => NULL(), ekh2ti(:) => NULL(), &
          & ekh3(:,:,:) => NULL(), ekh3vi(:,:) => NULL(), ekh3hi(:) => NULL(), ekh3ti(:) => NULL(), &
          !
          ! - vertical kinetic energy
          !
          & ekv1(:,:,:) => NULL(), ekv1vi(:,:) => NULL(), ekv1hi(:) => NULL(), ekv1ti(:) => NULL(), &
          & ekv2(:,:,:) => NULL(), ekv2vi(:,:) => NULL(), ekv2hi(:) => NULL(), ekv2ti(:) => NULL(), &
          & ekv3(:,:,:) => NULL(), ekv3vi(:,:) => NULL(), ekv3hi(:) => NULL(), ekv3ti(:) => NULL(), &
          !
          ! - potential energy
          !
          & egp1(:,:,:) => NULL(), egp1vi(:,:) => NULL(), egp1hi(:) => NULL(), egp1ti(:) => NULL(), &
          & egp2(:,:,:) => NULL(), egp2vi(:,:) => NULL(), egp2hi(:) => NULL(), egp2ti(:) => NULL(), &
          & egp3(:,:,:) => NULL(), egp3vi(:,:) => NULL(), egp3hi(:) => NULL(), egp3ti(:) => NULL(), &
          !
          ! - total energy
          !
          & eto1(:,:,:) => NULL(), eto1vi(:,:) => NULL(), eto1hi(:) => NULL(), eto1ti(:) => NULL(), &
          & eto2(:,:,:) => NULL(), eto2vi(:,:) => NULL(), eto2hi(:) => NULL(), eto2ti(:) => NULL(), &
          & eto3(:,:,:) => NULL(), eto3vi(:,:) => NULL(), eto3hi(:) => NULL(), eto3ti(:) => NULL()
     !
     !
     ! atmosphere energy tendency
     ! ==========================
     !
     !   ..dyn.. : due to dynamics            from difference ..2.. - ..1..
     !   ..phy.. : due to physics             from difference ..2.. - ..3..
     !   ..cld.. : due to cloud microphysics  from difference ..2.. - ..1..
     !   ..rad.. : due to radiation           from difference ..2.. - ..1..
     !   ..tmx.. : due to turbulent mixing    from difference ..2.. - ..1..
     !
     REAL(wp), POINTER ::                                                                                     &
          !
          !   dynamics
          !
          &   eindyn(:,:,:) => NULL(), eindynvi(:,:) => NULL(), eindynhi(:) => NULL(), eindynti(:) => NULL(), &
          &   ekhdyn(:,:,:) => NULL(), ekhdynvi(:,:) => NULL(), ekhdynhi(:) => NULL(), ekhdynti(:) => NULL(), &
          &   ekvdyn(:,:,:) => NULL(), ekvdynvi(:,:) => NULL(), ekvdynhi(:) => NULL(), ekvdynti(:) => NULL(), &
          &   egpdyn(:,:,:) => NULL(), egpdynvi(:,:) => NULL(), egpdynhi(:) => NULL(), egpdynti(:) => NULL(), &
          &   etodyn(:,:,:) => NULL(), etodynvi(:,:) => NULL(), etodynhi(:) => NULL(), etodynti(:) => NULL(), &
          !
          !   microphysics
          !
          &   eincld(:,:,:) => NULL(), eincldvi(:,:) => NULL(), eincldhi(:) => NULL(), eincldti(:) => NULL(), &
          !
          !   radiation
          !
          &   einrad(:,:,:) => NULL(), einradvi(:,:) => NULL(), einradhi(:) => NULL(), einradti(:) => NULL(), &
          !
          !   turbulent mixing
          !
          &   eintmx(:,:,:) => NULL(), eintmxvi(:,:) => NULL(), eintmxhi(:) => NULL(), eintmxti(:) => NULL(), &
          &   ekhtmx(:,:,:) => NULL(), ekhtmxvi(:,:) => NULL(), ekhtmxhi(:) => NULL(), ekhtmxti(:) => NULL(), &
          &   ekvtmx(:,:,:) => NULL(), ekvtmxvi(:,:) => NULL(), ekvtmxhi(:) => NULL(), ekvtmxti(:) => NULL(), &
          !
          !   physics
          !
          &   einphy(:,:,:) => NULL(), einphyvi(:,:) => NULL(), einphyhi(:) => NULL(), einphyti(:) => NULL(), &
          &   ekhphy(:,:,:) => NULL(), ekhphyvi(:,:) => NULL(), ekhphyhi(:) => NULL(), ekhphyti(:) => NULL(), & 
          &   ekvphy(:,:,:) => NULL(), ekvphyvi(:,:) => NULL(), ekvphyhi(:) => NULL(), ekvphyti(:) => NULL(), & 
          &   egpphy(:,:,:) => NULL(), egpphyvi(:,:) => NULL(), egpphyhi(:) => NULL(), egpphyti(:) => NULL(), &
          &   etophy(:,:,:) => NULL(), etophyvi(:,:) => NULL(), etophyhi(:) => NULL(), etophyti(:) => NULL()

  END TYPE t_atm_energy

  !-----------------------------------------------------------------------------

END MODULE mo_atm_energy_types
