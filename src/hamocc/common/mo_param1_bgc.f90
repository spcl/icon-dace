! @brief bgc tracer parameters
!
! Definition of parameters
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

MODULE mo_param1_bgc
  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC


  ! advected tracers
  INTEGER, PARAMETER :: i_base_adv = 18,              &
       &                isco212    = 1,               &
       &                ialkali    = 2,               &
       &                iphosph    = 3,               &
       &                iano3      = 4,               &
       &                igasnit    = 5,               &
       &                iphy       = 6,               &
       &                izoo       = 7,               &
       &                icya       = 8,               &
       &                ioxygen    = 9,               &
       &                isilica    = 10,              &
       &                idoc       = 11,              &
       &                ian2o      = 12,              &
       &                idet       = 13,              &         
       &                iiron      = 14,              &               
       &                idms       = 15,              &
       &                ih2s       = 16,              &
       &                idust      = 17,              &
       &                iagesc     = 18

  ! extended N-cycle parameters
  INTEGER :: iammo, iano2, i_amm_adv

  ! total number of advected tracers
  INTEGER :: ntraad

  ! non-advected (fast sinking) tracers
  INTEGER ::             icalc,         &
       &                 iopal,         &
       &                 i_base  


  INTEGER :: n_bgctra



  ! atmosphere
  INTEGER, PARAMETER :: iatmco2    = 1,               &
       &                iatmo2     = 2,               &
       &                iatmn2     = 3,               &
       &                i_base_atm = 3

  INTEGER, PARAMETER :: natm = i_base_atm 

  ! sediment
  INTEGER, PARAMETER :: issso12   = 1,                &
       &                isssc12   = 2,                &
       &                issssil   = 3,                &
       &                issster   = 4,                &
       &                nsss_base = 4
  INTEGER, PARAMETER :: nsedtra = nsss_base 

  ! pore water tracers, index must be the same as for ocetra otherwise problems in dipowa.f90!
  INTEGER, PARAMETER :: ipowaic    = 1,               &
       &                ipowaal    = 2,               &
       &                ipowaph    = 3,               &
       &                ipowaox    = 4,               &
       &                ipown2     = 5,               &
       &                ipowno3    = 6,               &
       &                ipowasi    = 7,               &
       &                ipowafe    = 8,               &
       &                ipowh2s    = 9,               &
       &                npowa_base = 9

  ! Extended N-cycle parameters
  INTEGER :: ipownh4, ipowno2, npowa_ammo
  
  INTEGER :: npowtra

 ! Diagnostics

  INTEGER, PARAMETER :: kphosy    = 1,               &
      &                 ksred     = 2,               &
      &                 kdenit    = 3,               &
      &                 kremin    = 4,               &
      &                 kgraz     = 5,               &
      &                 knfix     = 6,               &
      &                 kbacfra   = 7,               &
      &                 kpho_cya  = 8,               &
      &                 kcyaloss  = 9,               &
      &                 kn2b      = 10,               &
      &                 kh2ob     = 11,               &
      &                 kdelsil   = 12,               &
      &                 kdelcar   = 13,               &
      &                 kbacfrac  = 14,               &
      &                 kdmsprod  = 15,               &
      &                 kdmsbac   = 16,               &
      &                 kdmsuv    = 17,               &
      &                 keuexp    = 18,               &
      &                 kplim     = 19,               &
      &                 kflim     = 20,               &
      &                 knlim     = 21,               &
      &                 kgraton   = 22,               &
      &                 kexudp    = 23,               &
      &                 kexudz    = 24,               &
      &                 kzdy      = 25,               &
      &                 kpdy      = 26,               &
      &                 kaou      = 27,               &
      &                 kcLlim    = 28,               &
      &                 kcTlim    = 29,               &
      &                 kcPlim    = 30,               &
      &                 kcFlim    = 31,               &
      &                 kh2sprod  = 32,               &
      &                 kh2sloss  = 33,               &
      &                 kwpoc     = 34,               &
      &                 kwopal    = 35,               &
      &                 kwcal    = 36,               &
      &                 kwdust    = 37,               &
      &                 nbgctend_base  = 37 
 
  INTEGER, PARAMETER :: kcflux     = 1,               &
      &                 koflux     = 2,               &
      &                 knflux     = 3,               &
      &                 knfixd     = 4,               &
      &                 kprcaca    = 5,               &
      &                 kprorca    = 6,               &
      &                 ksilpro    = 7,               &
      &                 kcoex90    = 8,               &
      &                 kcalex90   = 9,               &
      &                 kopex90    = 10,               &
      &                 kn2oflux   = 11,               &
      &                 korginp    = 12,              &
      &                 ksilinp    = 13,              &
      &                 kcalinp    = 14,              &
      &                 kprodus    = 15,               &
      &                 kdmsflux   = 16,               &
      &                 kcoex1000  = 17,               &
      &                 kcalex1000 = 18,               &
      &                 kopex1000  = 19,               &
      &                 kcoex2000  = 20,               &
      &                 kcalex2000 = 21,               &
      &                 kopex2000  = 22,               &
      &                 kpco2      = 23,               &
      &                 klysocl    = 24,               &
      &                 knitinp    = 25,               &
      &                 kcflux_cpl = 26,               &
      &                 nbgcflux_base   = 26  

  INTEGER :: nbgcflux, nbgctend

  ! for extended N-cycle
  INTEGER :: nbgcflux_ammo, nbgctend_ammo
  INTEGER :: knh3flux
  INTEGER :: kgppnh, kcyapro, kdnrn, kdnra, kradnrn, kradnra
  INTEGER :: kanam, kraanam, kraanamox, kranitox, kammox, kraammox
  INTEGER :: knitox, kradeni

  ! Sediment
  INTEGER, PARAMETER::  &
      &                 nsed_diag_base =3,&
      &                 isremino =1,&
      &                 isreminn =2,&
      &                 isremins =3   

  INTEGER :: nsed_diag, nsed_diag_ammo
  INTEGER :: ksammox, ksnitox, ksanam, ksdnrn, ksdnra, ksnrn2

CONTAINS

  SUBROUTINE set_tracer_indices
    USE mo_hamocc_nml, ONLY : l_N_cycle,i_settling

      ! water column tracers (advected)
      iammo      = MERGE(i_base_adv+1,0,l_N_cycle)
      iano2      = MERGE(i_base_adv+2,0,l_N_cycle)
      i_amm_adv  = MERGE(2,0,l_N_cycle)
      ntraad = i_base_adv + i_amm_adv

      ! water column tracers (non-advected) -> need to be the last tracers
      ! within the 4D tracer array
      icalc    = ntraad+1
      iopal    = ntraad+2
      i_base   = 2
      n_bgctra = ntraad + i_base

      IF (i_settling==2) THEN
          ntraad=n_bgctra
      END IF
      
      ! porewater tracers
      ipownh4  = MERGE(npowa_base+1,0,l_N_cycle)
      ipowno2  = MERGE(npowa_base+2,0,l_N_cycle)
      npowa_ammo = MERGE(2,0,l_N_cycle)
      npowtra = npowa_base + npowa_ammo

      ! 3D diagnostic indices
      kgppnh = MERGE(nbgctend_base+1,0,l_N_cycle)      ! pp on nh4
      kcyapro = MERGE(nbgctend_base+2,0,l_N_cycle)     ! cyano pp on nh4 and no2
      kdnrn = MERGE(nbgctend_base+3,0,l_N_cycle)       ! dnrn of no3 to no2
      kdnra = MERGE(nbgctend_base+4,0,l_N_cycle)       ! dnrn of no3 to nh4
      kanam = MERGE(nbgctend_base+5,0,l_N_cycle)       ! anammox of nh4 and no2
      kammox = MERGE(nbgctend_base+6,0,l_N_cycle)     ! nitrification of nh4
      knitox = MERGE(nbgctend_base+7,0,l_N_cycle)     ! nitrification of no2
      ! kradnrn = MERGE(nbgctend_base+8,0,l_N_cycle)     ! LR: not include day-1 variables for now
      ! kradnra = MERGE(nbgctend_base+9,0,l_N_cycle)     ! LR: not include day-1 variables for now
      ! kraanam = MERGE(nbgctend_base+10,0,l_N_cycle)     ! LR: not include day-1 variables for now
      ! kranitox = MERGE(nbgctend_base+11,0,l_N_cycle)   ! LR: not include day-1 variables for now
      ! kraammox = MERGE(nbgctend_base+12,0,l_N_cycle)   ! LR: not include day-1 variables for now
      ! kradeni = MERGE(nbgctend_base+13,0,l_N_cycle)    ! LR: not include day-1 variables for now

      nbgctend_ammo = MERGE(7,0,l_N_cycle)
!      nbgctend_ammo = MERGE(13,0,l_N_cycle)

      nbgctend = nbgctend_base + nbgctend_ammo

      ! 2D diagnostoc indices
      knh3flux = MERGE(nbgcflux_base+1,0,l_N_cycle)
      nbgcflux_ammo = MERGE(1,0,l_N_cycle)
      nbgcflux = nbgcflux_base + nbgcflux_ammo

      ! Sediment diagnostic indices
      ksammox = MERGE(nsed_diag_base+1,0,l_N_cycle) ! sed nitrification of nh4
      ksnitox = MERGE(nsed_diag_base+2,0,l_N_cycle) ! sed nitrification of no2
      ksanam = MERGE(nsed_diag_base+3,0,l_N_cycle)  ! sed anammox of nh4 and no2
      ksdnrn = MERGE(nsed_diag_base+4,0,l_N_cycle)  ! sed dnrn of no3 to nh4
      ksdnra = MERGE(nsed_diag_base+5,0,l_N_cycle)  ! sed dnrn of no3 to no2
      ksnrn2 = MERGE(nsed_diag_base+6,0,l_N_cycle)  ! sed denitrification of no2 (in wc this uses the normal denitr. tend)
      nsed_diag_ammo = MERGE(6,0,l_N_cycle)
      nsed_diag = nsed_diag_base + nsed_diag_ammo

  END SUBROUTINE set_tracer_indices
END MODULE mo_param1_bgc


