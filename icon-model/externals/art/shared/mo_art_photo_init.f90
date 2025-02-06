!
! mo_art_photo_init
! This module provides initialization routines for CLOUDJ
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

MODULE mo_art_photo_init
  USE FJX_INIT_MOD,      ONLY: INIT_FJX
  USE FJX_CMN_MOD,       ONLY: JXL1_, AN_, JVN_, S_
  USE mo_art_chem_data,  ONLY: t_art_photolysis,  &
    &                          nphot
  USE mo_exception,                     ONLY: finish
    
  IMPLICIT NONE
  
  PRIVATE

  CHARACTER*6, DIMENSION(JVN_)   :: TITLJXX
  INTEGER                        :: NJXX 

  PUBLIC ::art_photo_init, NJXX
  PUBLIC ::TITLJXX
  PUBLIC ::art_photo_deallocate

CONTAINS

SUBROUTINE art_photo_init(art_photo, nproma, nlev, nblks)

!<
! SUBROUTINE art_photo init
! This subroutine provides an interface to initialization routines in
! CLOUDJ. 
! Part of Module: mo_art_photo_init
! Author: Jennifer Schroeter, KIT
! Initial Release: 2016-09-22
! Modifications:
!>
  TYPE(t_art_photolysis), INTENT(inout) :: &
    &  art_photo                   !< Pointer to ART photolysis fields
  INTEGER, INTENT(in) :: &
    &  nproma, nlev, nblks
  ! local variables
  INTEGER :: &
    &  L_, L1_

  call INIT_FJX (TITLJXX,JVN_,NJXX)

  ALLOCATE(art_photo%heating_rates(nproma,nlev,nblks,S_+2))

  ALLOCATE(art_photo%input_fastjx_clf (nproma,nlev,nblks))
  ALLOCATE(art_photo%input_fastjx_rh (nproma,nlev,nblks))
  ALLOCATE(art_photo%input_fastjx_lwp (nproma,nlev,nblks))
  ALLOCATE(art_photo%input_fastjx_iwp (nproma,nlev,nblks))
  ALLOCATE(art_photo%output_fastjx_photo_new (nproma,nlev,nblks,nphot))
  ALLOCATE(art_photo%input_fastjx_reice (nproma,nlev,nblks))
  ALLOCATE(art_photo%input_fastjx_reliq (nproma,nlev,nblks))
  ALLOCATE(art_photo%fjx_zkap (nproma,nblks))
  ALLOCATE(art_photo%fjx_zland (nproma,nblks))
  ALLOCATE(art_photo%fjx_zglac (nproma,nblks))
  ALLOCATE(art_photo%input_fastjx_lwc (nproma,nlev,nblks))
  ALLOCATE(art_photo%input_fastjx_iwc (nproma,nlev,nblks))
  ALLOCATE(art_photo%fjx_cdnc(nproma,nlev,nblks))

  L_     = nlev-1
  L1_    = L_+1

  ALLOCATE(art_photo%PPP(L1_+1))
  ALLOCATE(art_photo%ZZZ(L1_+1))
  ALLOCATE(art_photo%DDD(L1_)) 
  ALLOCATE(art_photo%TTT(L1_)) 
  ALLOCATE(art_photo%RRR(L1_)) 
  ALLOCATE(art_photo%OOO(L1_)) 
  ALLOCATE(art_photo%LWP(L1_)) 
  ALLOCATE(art_photo%IWP(L1_)) 
  ALLOCATE(art_photo%REFFL(L1_))
  ALLOCATE(art_photo%REFFI(L1_))
  ALLOCATE(art_photo%CLF(L1_)) 
  ALLOCATE(art_photo%CLDIW(L1_))
  ALLOCATE(art_photo%AERSP(L1_,AN_))
  ALLOCATE(art_photo%NDXAER(L1_,AN_))
  ALLOCATE(art_photo%VALJXX(L_,JVN_))
  ALLOCATE(art_photo%SKPERD(S_+2,JXL1_))

END SUBROUTINE art_photo_init

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE art_photo_deallocate(art_photo)

!<
! SUBROUTINE art_photo_deallocate
! This subroutine deallocates the arrays for photolysis
! Part of Module: mo_art_photo_init
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-12
! Modifications:
!>
  TYPE(t_art_photolysis), INTENT(inout)    :: &
    &  art_photo                   !< ART photolysis fields

  DEALLOCATE(art_photo%heating_rates)
  DEALLOCATE(art_photo%input_fastjx_clf )
  DEALLOCATE(art_photo%input_fastjx_rh )
  DEALLOCATE(art_photo%input_fastjx_lwp )
  DEALLOCATE(art_photo%input_fastjx_iwp )
  DEALLOCATE(art_photo%output_fastjx_photo_new )
  DEALLOCATE(art_photo%input_fastjx_reice )
  DEALLOCATE(art_photo%input_fastjx_reliq )
  DEALLOCATE(art_photo%fjx_zkap )
  DEALLOCATE(art_photo%fjx_zland )
  DEALLOCATE(art_photo%fjx_zglac )
  DEALLOCATE(art_photo%input_fastjx_lwc )
  DEALLOCATE(art_photo%input_fastjx_iwc )
  DEALLOCATE(art_photo%fjx_cdnc)
  DEALLOCATE(art_photo%PPP)
  DEALLOCATE(art_photo%ZZZ)
  DEALLOCATE(art_photo%DDD) 
  DEALLOCATE(art_photo%TTT) 
  DEALLOCATE(art_photo%RRR) 
  DEALLOCATE(art_photo%OOO) 
  DEALLOCATE(art_photo%LWP) 
  DEALLOCATE(art_photo%IWP) 
  DEALLOCATE(art_photo%REFFL)
  DEALLOCATE(art_photo%REFFI)
  DEALLOCATE(art_photo%CLF) 
  DEALLOCATE(art_photo%CLDIW)
  DEALLOCATE(art_photo%AERSP)
  DEALLOCATE(art_photo%NDXAER)
  DEALLOCATE(art_photo%VALJXX)
  DEALLOCATE(art_photo%SKPERD)
END SUBROUTINE art_photo_deallocate

END MODULE
