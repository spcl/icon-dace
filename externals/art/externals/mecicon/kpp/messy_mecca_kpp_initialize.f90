! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Initialization File
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
! SPDX-License-Identifier: GPL-3.0-only  
! ---------------------------------------------------------------

MODULE messy_mecca_kpp_initialize 

 USE messy_mecca_kpp_parameters, ONLY: NVAR, NFIX 
 USE mo_kind,                    ONLY: dp 

 IMPLICIT NONE 
 PUBLIC

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialize - function to initialize concentrations
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE initialize ( )


  USE messy_mecca_kpp_global

  INTEGER :: i
  REAL(kind=dp) :: x

  CFACTOR = 1.000000e+00_dp

  x = (0.)*CFACTOR
  DO i = 1, NVAR
    VAR(i) = x
  END DO

  x = (0.)*CFACTOR
  DO i = 1, NFIX
    FIX(i) = x
  END DO

! constant rate coefficients
! END constant rate coefficients

! INLINED initializations

! End INLINED initializations

      
END SUBROUTINE initialize

! End of Initialize function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 END MODULE messy_mecca_kpp_initialize 

