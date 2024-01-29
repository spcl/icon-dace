!---------------------------------------------------------
!---Read LINOZ tables from 'Linoz2004Br.dat' for linearized Ozonechemistry
!---------------------------------------------------------
!
! mo_art_read_linoz
! This module provides a simple linear ozone chemistry
! first introduced by McLinden 2000
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

MODULE mo_art_read_linoz
  USE mo_art_chem_data,                 ONLY: t_art_linoz

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: art_linoz_read
  PUBLIC  :: art_linoz_deallocate

CONTAINS

SUBROUTINE art_linoz_read(art_linoz, nproma, nlev)

!<
! SUBROUTINE art_linoz_read
! Routine reads in values needed for linearized ansatz
! based on McLinden et al.
! File 'Linoz2004Br.dat' has to be linked in the output directory.
! Part of Module: mo_art_read_linoz
! Author: Christian Stassen, KIT
! Initial Release: 2015-03-24
!>

  IMPLICIT NONE 

  INTEGER, INTENT(in) :: &
    &  nproma, nlev              !< dimensions of arrays


  TYPE(t_art_linoz), INTENT(inout)    :: &
    &  art_linoz                  !< Pointer to ART chem fields

  character*50 TITEL
  character*120 linozinifile
  integer n,m,k,l


  IF (.NOT. ALLOCATED (art_linoz%tparm)) ALLOCATE(art_linoz%tparm(25,18,12,7))
  linozinifile = 'Linoz2004Br.dat'

  open(61,file=linozinifile, STATUS = 'old')
  read(61,901) TITEL
  do n=1,7
      read(61,901) TITEL
      do m=1,12
          do k=1,18
              read(61,902) (art_linoz%tparm(l,k,m,n),l=1,25)
          enddo
      enddo
  enddo
 
  901 format(A50)
  902 format(20x,6e10.3/(8e10.3))
  close(61)


  IF (.NOT. ALLOCATED(art_linoz%linoz_tab1)) ALLOCATE(art_linoz%linoz_tab1(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linoz%linoz_tab2)) ALLOCATE(art_linoz%linoz_tab2(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linoz%linoz_tab3)) ALLOCATE(art_linoz%linoz_tab3(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linoz%linoz_tab4)) ALLOCATE(art_linoz%linoz_tab4(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linoz%linoz_tab5)) ALLOCATE(art_linoz%linoz_tab5(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linoz%linoz_tab6)) ALLOCATE(art_linoz%linoz_tab6(nproma,nlev))
  IF (.NOT. ALLOCATED(art_linoz%linoz_tab7)) ALLOCATE(art_linoz%linoz_tab7(nproma,nlev))

  art_linoz%is_init = .TRUE.



END SUBROUTINE art_linoz_read

SUBROUTINE art_linoz_deallocate(art_linoz)
!<
! SUBROUTINE linoz_read
! This routines deallocates all allocated arrays of LINOZ
! Part of Module: mo_art_read_linoz
! Author: Michael Weimer, KIT
! Initial Release: 2020-06-12
!>
  TYPE(t_art_linoz), INTENT(inout)    :: &
    &  art_linoz                    !< Pointer to ART chem fields

  DEALLOCATE(art_linoz%linoz_tab1)
  DEALLOCATE(art_linoz%linoz_tab2)
  DEALLOCATE(art_linoz%linoz_tab3)
  DEALLOCATE(art_linoz%linoz_tab4)
  DEALLOCATE(art_linoz%linoz_tab5)
  DEALLOCATE(art_linoz%linoz_tab6)
  DEALLOCATE(art_linoz%linoz_tab7)
  DEALLOCATE(art_linoz%tparm)

  art_linoz%is_init = .FALSE.

END SUBROUTINE art_linoz_deallocate
END MODULE mo_art_read_linoz
