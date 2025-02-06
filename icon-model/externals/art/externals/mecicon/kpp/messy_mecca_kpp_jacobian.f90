! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! The ODE Jacobian of Chemical Model File
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

MODULE messy_mecca_kpp_Jacobian 
  USE mo_kind,                 ONLY: dp

  USE messy_mecca_kpp_parameters
  USE messy_mecca_kpp_JacobianSP 
  USE mo_kind, ONLY:dp

  IMPLICIT NONE
  PUBLIC :: Jac_SP

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP - the Jacobian of Variables in sparse matrix representation
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      JVS       - sparse Jacobian of variables
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP ( V, F, RCT, JVS )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)


! Local variables
! B - Temporary array
  REAL(kind=dp) :: B(37)

! B(1) = dA(1)/dV(12)
  B(1) = RCT(1)*F(1)
! B(3) = dA(2)/dV(8)
  B(3) = RCT(2)*F(1)
! B(5) = dA(3)/dV(8)
  B(5) = RCT(3)*V(11)
! B(6) = dA(3)/dV(11)
  B(6) = RCT(3)*V(8)
! B(7) = dA(4)/dV(10)
  B(7) = RCT(4)*V(11)
! B(8) = dA(4)/dV(11)
  B(8) = RCT(4)*V(10)
! B(9) = dA(5)/dV(3)
  B(9) = RCT(5)*V(11)
! B(10) = dA(5)/dV(11)
  B(10) = RCT(5)*V(3)
! B(11) = dA(6)/dV(3)
  B(11) = RCT(6)*V(10)
! B(12) = dA(6)/dV(10)
  B(12) = RCT(6)*V(3)
! B(13) = dA(7)/dV(4)
  B(13) = RCT(7)*V(12)
! B(14) = dA(7)/dV(12)
  B(14) = RCT(7)*V(4)
! B(15) = dA(8)/dV(12)
  B(15) = RCT(8)*F(2)
! B(17) = dA(9)/dV(1)
  B(17) = RCT(9)*V(12)
! B(18) = dA(9)/dV(12)
  B(18) = RCT(9)*V(1)
! B(19) = dA(10)/dV(5)
  B(19) = RCT(10)*V(11)
! B(20) = dA(10)/dV(11)
  B(20) = RCT(10)*V(5)
! B(21) = dA(11)/dV(8)
  B(21) = RCT(11)*V(9)
! B(22) = dA(11)/dV(9)
  B(22) = RCT(11)*V(8)
! B(23) = dA(12)/dV(9)
  B(23) = RCT(12)*V(11)
! B(24) = dA(12)/dV(11)
  B(24) = RCT(12)*V(9)
! B(25) = dA(13)/dV(6)
  B(25) = RCT(13)*V(9)
! B(26) = dA(13)/dV(9)
  B(26) = RCT(13)*V(6)
! B(27) = dA(14)/dV(9)
  B(27) = RCT(14)*V(10)
! B(28) = dA(14)/dV(10)
  B(28) = RCT(14)*V(9)
! B(29) = dA(15)/dV(7)
  B(29) = RCT(15)*V(10)
! B(30) = dA(15)/dV(10)
  B(30) = RCT(15)*V(7)
! B(32) = dA(17)/dV(11)
  B(32) = RCT(17)
! B(33) = dA(18)/dV(11)
  B(33) = RCT(18)
! B(34) = dA(19)/dV(9)
  B(34) = RCT(19)
! B(35) = dA(20)/dV(6)
  B(35) = RCT(20)
! B(36) = dA(21)/dV(2)
  B(36) = RCT(21)
! B(37) = dA(22)/dV(7)
  B(37) = RCT(22)

! Construct the Jacobian terms from B's
! JVS(1) = Jac_FULL(1,1)
  JVS(1) = -B(17)
! JVS(2) = Jac_FULL(1,12)
  JVS(2) = -B(18)
! JVS(3) = Jac_FULL(2,2)
  JVS(3) = -B(36)
! JVS(4) = Jac_FULL(2,6)
  JVS(4) = B(25)
! JVS(5) = Jac_FULL(2,9)
  JVS(5) = B(26)
! JVS(6) = Jac_FULL(3,3)
  JVS(6) = -B(9)-B(11)
! JVS(7) = Jac_FULL(3,10)
  JVS(7) = B(7)-B(12)
! JVS(8) = Jac_FULL(3,11)
  JVS(8) = B(8)-B(10)
! JVS(9) = Jac_FULL(4,3)
  JVS(9) = B(11)
! JVS(10) = Jac_FULL(4,4)
  JVS(10) = -B(13)
! JVS(11) = Jac_FULL(4,7)
  JVS(11) = B(29)
! JVS(12) = Jac_FULL(4,10)
  JVS(12) = B(12)+B(30)
! JVS(13) = Jac_FULL(4,11)
  JVS(13) = 0
! JVS(14) = Jac_FULL(4,12)
  JVS(14) = -B(14)
! JVS(15) = Jac_FULL(5,1)
  JVS(15) = 2*B(17)
! JVS(16) = Jac_FULL(5,5)
  JVS(16) = -B(19)
! JVS(17) = Jac_FULL(5,8)
  JVS(17) = B(21)
! JVS(18) = Jac_FULL(5,9)
  JVS(18) = B(22)+B(34)
! JVS(19) = Jac_FULL(5,11)
  JVS(19) = -B(20)
! JVS(20) = Jac_FULL(5,12)
  JVS(20) = 2*B(18)
! JVS(21) = Jac_FULL(6,2)
  JVS(21) = B(36)
! JVS(22) = Jac_FULL(6,6)
  JVS(22) = -B(25)-B(35)
! JVS(23) = Jac_FULL(6,7)
  JVS(23) = B(29)
! JVS(24) = Jac_FULL(6,9)
  JVS(24) = B(23)-B(26)
! JVS(25) = Jac_FULL(6,10)
  JVS(25) = B(30)
! JVS(26) = Jac_FULL(6,11)
  JVS(26) = B(24)
! JVS(27) = Jac_FULL(7,7)
  JVS(27) = -B(29)-B(37)
! JVS(28) = Jac_FULL(7,9)
  JVS(28) = B(27)
! JVS(29) = Jac_FULL(7,10)
  JVS(29) = B(28)-B(30)
! JVS(30) = Jac_FULL(8,6)
  JVS(30) = B(35)
! JVS(31) = Jac_FULL(8,7)
  JVS(31) = 0
! JVS(32) = Jac_FULL(8,8)
  JVS(32) = -B(3)-B(5)-B(21)
! JVS(33) = Jac_FULL(8,9)
  JVS(33) = -B(22)+B(34)
! JVS(34) = Jac_FULL(8,10)
  JVS(34) = 0
! JVS(35) = Jac_FULL(8,11)
  JVS(35) = -B(6)+B(33)
! JVS(36) = Jac_FULL(8,12)
  JVS(36) = B(1)+B(15)
! JVS(37) = Jac_FULL(9,2)
  JVS(37) = B(36)
! JVS(38) = Jac_FULL(9,5)
  JVS(38) = B(19)
! JVS(39) = Jac_FULL(9,6)
  JVS(39) = -B(25)+B(35)
! JVS(40) = Jac_FULL(9,7)
  JVS(40) = B(37)
! JVS(41) = Jac_FULL(9,8)
  JVS(41) = -B(21)
! JVS(42) = Jac_FULL(9,9)
  JVS(42) = -B(22)-B(23)-B(26)-B(27)-B(34)
! JVS(43) = Jac_FULL(9,10)
  JVS(43) = -B(28)
! JVS(44) = Jac_FULL(9,11)
  JVS(44) = B(20)-B(24)
! JVS(45) = Jac_FULL(9,12)
  JVS(45) = 0
! JVS(46) = Jac_FULL(10,3)
  JVS(46) = B(9)-B(11)
! JVS(47) = Jac_FULL(10,4)
  JVS(47) = 2*B(13)
! JVS(48) = Jac_FULL(10,7)
  JVS(48) = -B(29)+B(37)
! JVS(49) = Jac_FULL(10,9)
  JVS(49) = -B(27)
! JVS(50) = Jac_FULL(10,10)
  JVS(50) = -B(7)-B(12)-B(28)-B(30)
! JVS(51) = Jac_FULL(10,11)
  JVS(51) = -B(8)+B(10)
! JVS(52) = Jac_FULL(10,12)
  JVS(52) = 2*B(14)
! JVS(53) = Jac_FULL(11,3)
  JVS(53) = -B(9)
! JVS(54) = Jac_FULL(11,5)
  JVS(54) = -B(19)
! JVS(55) = Jac_FULL(11,8)
  JVS(55) = B(3)-B(5)
! JVS(56) = Jac_FULL(11,9)
  JVS(56) = -B(23)
! JVS(57) = Jac_FULL(11,10)
  JVS(57) = -B(7)
! JVS(58) = Jac_FULL(11,11)
  JVS(58) = -B(6)-B(8)-B(10)-B(20)-B(24)-B(32)-B(33)
! JVS(59) = Jac_FULL(11,12)
  JVS(59) = 0
! JVS(60) = Jac_FULL(12,1)
  JVS(60) = -B(17)
! JVS(61) = Jac_FULL(12,4)
  JVS(61) = -B(13)
! JVS(62) = Jac_FULL(12,7)
  JVS(62) = 0
! JVS(63) = Jac_FULL(12,9)
  JVS(63) = 0
! JVS(64) = Jac_FULL(12,10)
  JVS(64) = 0
! JVS(65) = Jac_FULL(12,11)
  JVS(65) = B(32)
! JVS(66) = Jac_FULL(12,12)
  JVS(66) = -B(1)-B(14)-B(15)-B(18)
      
END SUBROUTINE Jac_SP

! End of Jac_SP function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Jac_SP_Vec - function for sparse multiplication: sparse Jacobian times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JUV       - Jacobian times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Jac_SP_Vec ( JVS, UV, JUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JUV - Jacobian times user vector
  REAL(kind=dp) :: JUV(NVAR)

  JUV(1) = JVS(1)*UV(1)+JVS(2)*UV(12)
  JUV(2) = JVS(3)*UV(2)+JVS(4)*UV(6)+JVS(5)*UV(9)
  JUV(3) = JVS(6)*UV(3)+JVS(7)*UV(10)+JVS(8)*UV(11)
  JUV(4) = JVS(9)*UV(3)+JVS(10)*UV(4)+JVS(11)*UV(7)+JVS(12)*UV(10)+JVS(14)*UV(12)
  JUV(5) = JVS(15)*UV(1)+JVS(16)*UV(5)+JVS(17)*UV(8)+JVS(18)*UV(9)+JVS(19)*UV(11)+JVS(20)*UV(12)
  JUV(6) = JVS(21)*UV(2)+JVS(22)*UV(6)+JVS(23)*UV(7)+JVS(24)*UV(9)+JVS(25)*UV(10)+JVS(26)*UV(11)
  JUV(7) = JVS(27)*UV(7)+JVS(28)*UV(9)+JVS(29)*UV(10)
  JUV(8) = JVS(30)*UV(6)+JVS(32)*UV(8)+JVS(33)*UV(9)+JVS(35)*UV(11)+JVS(36)*UV(12)
  JUV(9) = JVS(37)*UV(2)+JVS(38)*UV(5)+JVS(39)*UV(6)+JVS(40)*UV(7)+JVS(41)*UV(8)+JVS(42)*UV(9)+JVS(43)*UV(10)+JVS(44)&
             &*UV(11)
  JUV(10) = JVS(46)*UV(3)+JVS(47)*UV(4)+JVS(48)*UV(7)+JVS(49)*UV(9)+JVS(50)*UV(10)+JVS(51)*UV(11)+JVS(52)*UV(12)
  JUV(11) = JVS(53)*UV(3)+JVS(54)*UV(5)+JVS(55)*UV(8)+JVS(56)*UV(9)+JVS(57)*UV(10)+JVS(58)*UV(11)
  JUV(12) = JVS(60)*UV(1)+JVS(61)*UV(4)+JVS(65)*UV(11)+JVS(66)*UV(12)
      
END SUBROUTINE Jac_SP_Vec

! End of Jac_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacTR_SP_Vec - sparse multiplication: sparse Jacobian transposed times vector
!   Arguments :
!      JVS       - sparse Jacobian of variables
!      UV        - User vector for variables
!      JTUV      - Jacobian transposed times user vector
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacTR_SP_Vec ( JVS, UV, JTUV )

! JVS - sparse Jacobian of variables
  REAL(kind=dp) :: JVS(LU_NONZERO)
! UV - User vector for variables
  REAL(kind=dp) :: UV(NVAR)
! JTUV - Jacobian transposed times user vector
  REAL(kind=dp) :: JTUV(NVAR)

  JTUV(1) = JVS(1)*UV(1)+JVS(15)*UV(5)+JVS(60)*UV(12)
  JTUV(2) = JVS(3)*UV(2)+JVS(21)*UV(6)+JVS(37)*UV(9)
  JTUV(3) = JVS(6)*UV(3)+JVS(9)*UV(4)+JVS(46)*UV(10)+JVS(53)*UV(11)
  JTUV(4) = JVS(10)*UV(4)+JVS(47)*UV(10)+JVS(61)*UV(12)
  JTUV(5) = JVS(16)*UV(5)+JVS(38)*UV(9)+JVS(54)*UV(11)
  JTUV(6) = JVS(4)*UV(2)+JVS(22)*UV(6)+JVS(30)*UV(8)+JVS(39)*UV(9)
  JTUV(7) = JVS(11)*UV(4)+JVS(23)*UV(6)+JVS(27)*UV(7)+JVS(40)*UV(9)+JVS(48)*UV(10)
  JTUV(8) = JVS(17)*UV(5)+JVS(32)*UV(8)+JVS(41)*UV(9)+JVS(55)*UV(11)
  JTUV(9) = JVS(5)*UV(2)+JVS(18)*UV(5)+JVS(24)*UV(6)+JVS(28)*UV(7)+JVS(33)*UV(8)+JVS(42)*UV(9)+JVS(49)*UV(10)+JVS(56)&
              &*UV(11)
  JTUV(10) = JVS(7)*UV(3)+JVS(12)*UV(4)+JVS(25)*UV(6)+JVS(29)*UV(7)+JVS(43)*UV(9)+JVS(50)*UV(10)+JVS(57)*UV(11)
  JTUV(11) = JVS(8)*UV(3)+JVS(19)*UV(5)+JVS(26)*UV(6)+JVS(35)*UV(8)+JVS(44)*UV(9)+JVS(51)*UV(10)+JVS(58)*UV(11)+JVS(65)&
               &*UV(12)
  JTUV(12) = JVS(2)*UV(1)+JVS(14)*UV(4)+JVS(20)*UV(5)+JVS(36)*UV(8)+JVS(52)*UV(10)+JVS(66)*UV(12)
      
END SUBROUTINE JacTR_SP_Vec

! End of JacTR_SP_Vec function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 END MODULE messy_mecca_kpp_Jacobian 

