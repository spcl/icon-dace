! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The ODE Function of Chemical Model File
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

MODULE messy_mecca_kpp_function 
  USE mo_kind,                 ONLY: dp

  USE messy_mecca_kpp_parameters
  IMPLICIT NONE 
  PUBLIC :: Fun
  PUBLIC :: A

! A - Rate for each equation
  REAL(kind=dp) :: A(NREACT)

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Fun - time derivatives of variables - Agregate form
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      Vdot      - Time derivative of variable species concentrations
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Fun ( V, F, RCT, Vdot )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! Vdot - Time derivative of variable species concentrations
  REAL(kind=dp) :: Vdot(NVAR)


! Computation of equation rates
  A(1) = RCT(1)*V(12)*F(1)
  A(2) = RCT(2)*V(8)*F(1)
  A(3) = RCT(3)*V(8)*V(11)
  A(4) = RCT(4)*V(10)*V(11)
  A(5) = RCT(5)*V(3)*V(11)
  A(6) = RCT(6)*V(3)*V(10)
  A(7) = RCT(7)*V(4)*V(12)
  A(8) = RCT(8)*V(12)*F(2)
  A(9) = RCT(9)*V(1)*V(12)
  A(10) = RCT(10)*V(5)*V(11)
  A(11) = RCT(11)*V(8)*V(9)
  A(12) = RCT(12)*V(9)*V(11)
  A(13) = RCT(13)*V(6)*V(9)
  A(14) = RCT(14)*V(9)*V(10)
  A(15) = RCT(15)*V(7)*V(10)
  A(16) = RCT(16)*F(1)
  A(17) = RCT(17)*V(11)
  A(18) = RCT(18)*V(11)
  A(19) = RCT(19)*V(9)
  A(20) = RCT(20)*V(6)
  A(21) = RCT(21)*V(2)
  A(22) = RCT(22)*V(7)

! Aggregate function
  Vdot(1) = -A(9)
  Vdot(2) = A(13)-A(21)
  Vdot(3) = A(4)-A(5)-A(6)
  Vdot(4) = A(6)-A(7)+A(15)
  Vdot(5) = 2*A(9)-A(10)+A(11)+A(19)
  Vdot(6) = A(12)-A(13)+A(15)-A(20)+A(21)
  Vdot(7) = A(14)-A(15)-A(22)
  Vdot(8) = A(1)-A(2)-A(3)+A(8)-A(11)+2*A(16)+A(18)+A(19)+A(20)
  Vdot(9) = A(10)-A(11)-A(12)-A(13)-A(14)-A(19)+A(20)+A(21)+A(22)
  Vdot(10) = -A(4)+A(5)-A(6)+2*A(7)-A(14)-A(15)+A(22)
  Vdot(11) = A(2)-A(3)-A(4)-A(5)-A(10)-A(12)-A(17)-A(18)
  Vdot(12) = -A(1)-A(7)-A(8)-A(9)+A(17)
      
END SUBROUTINE Fun

! End of Fun function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 END MODULE messy_mecca_kpp_function 

