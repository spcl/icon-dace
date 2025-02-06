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

! Contains surface type (tile) indices used by the turbulent mixing parameterization.
!
! Contains subroutines for initializing the AES physics package

MODULE mo_aes_sfc_indices

  USE mo_exception, ONLY: message, message_text

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nsfc_type, iwtr, iice, ilnd, igbm   !< index variables
  PUBLIC :: csfc                                !< sfc names
  PUBLIC :: init_sfc_indices                    !< subroutine

  INTEGER :: nsfc_type   !< total number of surface types
  INTEGER :: iwtr = 1    !< index for water-covered surface
  INTEGER :: iice = 2    !< index for ice-covered   surface
  INTEGER :: ilnd = 3    !< index for land          surface
  INTEGER :: igbm        !< index for grid-box mean

  CHARACTER(LEN=3) :: csfc(3) = (/'wtr','ice','lnd'/)

CONTAINS
  !>
  !! Set surface indices according to the simulation setup
  !! (e.g., dynamical core test, aqua-planet, or real
  !! climate simulation).
  !!
  SUBROUTINE init_sfc_indices( ctest_name )

    CHARACTER(len=*),INTENT(IN) :: ctest_name

    SELECT CASE(TRIM(ctest_name))
    CASE('APE','APE_aes','RCE','RCE_glb','RCE_Tconst','RCE_Tprescr','aes_bubble','RCEhydro','CBL_flxconst','RCEMIP_analytical', &
      &  'dcmip_tc_52')
      ! Aqua-planet simulation, no land, no ice;
      ! No needed to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

      iwtr      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      ilnd      = 999

    CASE('aes_bubble_land')

      ! Terra-planet bubble simulation, no ocean, no ice, no lakes, no glaciers;
      ! No need to distinguish the aggregated grid-box mean and the value on different types of surface

      ilnd      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      iwtr      = 999
      csfc(1:3) = (/'lnd', 'wtr', 'ice'/)

    CASE('APEi','APEc','APEc_nh')
      ! Aqua-planet simulation with ice, but no land;

      iwtr      = 1
      iice      = 2
      nsfc_type = 2
      igbm      = 0
      ilnd      = 999

    CASE('JWw-Moist','LDF-Moist','jabw_m')
      ! Baroclinic wave test, no land, no ice.

      iwtr      = 1
      nsfc_type = 1
      igbm      = 0
      iice      = 999
      ilnd      = 999

      ! maybe worth trying later:
      ! iwtr      = 1
      ! iice      = 2
      ! nsfc_type = 2
      ! igbm      = 0
      ! ilnd      = 999

      ! A wild idea: a dry or completely frozen planet
      ! iice      = 1
      ! ilnd      = 2
      ! nsfc_type = 2
      ! igbm      = 0
      ! iwtr      = 999

    CASE('TPEo','TPEc')
      ! Terra-planet simulation
      ! no ocean, no ice but lakes and ice on lakes ... therefore have to use iice and iwtr
      ! No need to distinguish the aggregated grid-box mean
      ! and the value on different types of surface

      nsfc_type = 3
      iwtr      = 1
      iice      = 2
      ilnd      = 3
      igbm      = 0

    CASE DEFAULT
      ! Standard setup for real-world climate simulation.
      ! Three surface types are considered.

      iwtr      = 1
      iice      = 2
      ilnd      = 3
      nsfc_type = 3
      igbm      = 0

    END SELECT


    WRITE(message_text,*) " "
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(i3,a)')    &
      & nsfc_type, " surface type(s) activated."
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,'(a,4i4,a)') &
      & "Indices for water, ice, land, and grid-box mean are ", &
      & iwtr, iice, ilnd, igbm, ", respectively."
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

    WRITE(message_text,*) " "
    CALL message("mo_aes_sfc_indices/init_sfc_indices",TRIM(message_text))

  END SUBROUTINE init_sfc_indices
  !-------------

END MODULE mo_aes_sfc_indices
