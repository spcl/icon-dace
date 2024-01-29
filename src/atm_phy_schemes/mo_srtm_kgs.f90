!----------------------------------------------------------------------------
! Copyright (c) 2002-2020, Atmospheric & Environmental Research, Inc. (AER)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC.,
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.
!                        (http://www.rtweb.aer.com/)
!----------------------------------------------------------------------------

! SPDX-License-Identifier: BSD-3-Clause
! Code has been modified for the use in ICON (and formerly in ECHAM)

MODULE mo_yoesrta16

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA16* - SRTM COEFFICIENTS FOR INTERVAL 16
  !     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg=16, ng16 = 16, ngs15 = 0

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(3,jpg)
  REAL(wp) :: sfluxref(jpg)
  REAL(wp) :: rayl            ,strrat1
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng16),absa(585,ng16)
  REAL(wp) :: kbc(5,13:59,ng16),absb(235,ng16)
  REAL(wp) :: selfrefc(10,ng16),forrefc(3,ng16)
  REAL(wp) :: sfluxrefc(ng16)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT1 : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta16

MODULE mo_yoesrta17

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA17* - SRTM COEFFICIENTS FOR INTERVAL 17
  !     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER ::  jpg = 16, ng17 = 16, ngs16=16

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg)  ,forref(4,jpg)
  REAL(wp) :: sfluxref(jpg,5)
  REAL(wp) :: rayl              ,strrat
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng17),absa(585,ng17)
  REAL(wp) :: kbc(5,5,13:59,ng17),absb(1175,ng17)
  REAL(wp) :: selfrefc(10,ng17),forrefc(4,ng17)
  REAL(wp) :: sfluxrefc(ng17,5)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,1,13,1),absb(1,1))
  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta17

MODULE mo_yoesrta18

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA18* - SRTM COEFFICIENTS FOR INTERVAL 16
  !     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng18 = 16, ngs17=32

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(3,jpg)
  REAL(wp) :: sfluxref(jpg,9)
  REAL(wp) :: rayl            ,strrat

  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng18) ,absa(585,ng18)
  REAL(wp) :: kbc(5,13:59,ng18),absb(235,ng18)
  REAL(wp) :: selfrefc(10,ng18),forrefc(3,ng18)
  REAL(wp) :: sfluxrefc(ng18,9)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta18

MODULE mo_yoesrta19

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA19* - SRTM COEFFICIENTS FOR INTERVAL 19
  !     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng19 = 16

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(3,jpg)
  REAL(wp) :: sfluxref(jpg,9)
  REAL(wp) :: rayl            ,strrat
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng19) ,absa(585,ng19)
  REAL(wp) :: kbc(5,13:59,ng19),absb(235,ng19)
  REAL(wp) :: selfrefc(10,ng19),forrefc(3,ng19)
  REAL(wp) :: sfluxrefc(ng19,9)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))
  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta19

MODULE mo_yoesrta20

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA20* - SRTM COEFFICIENTS FOR INTERVAL 20
  !     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng20 = 16

  REAL(wp) :: ka(5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(4,jpg)
  REAL(wp) :: sfluxref(jpg)  ,absch4(jpg)
  REAL(wp) :: rayl
  INTEGER  :: layreffr

  REAL(wp) :: kac(5,13,ng20)   ,absa(65,ng20)
  REAL(wp) :: kbc(5,13:59,ng20),absb(235,ng20)
  REAL(wp) :: selfrefc(10,ng20),forrefc(4,ng20)
  REAL(wp) :: sfluxrefc(ng20)  ,absch4c(ng20)

  !EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! ABSCH4  : REAL
  ! RAYL    : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! ABSCH4C : REAL     Reduced g-point array for ABSCH4
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta20

MODULE mo_yoesrta21

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA21* - SRTM COEFFICIENTS FOR INTERVAL 21
  !     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng21 = 16

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg)  ,forref(4,jpg)
  REAL(wp) :: sfluxref(jpg,9)
  REAL(wp) :: rayl              ,strrat
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng21)   ,absa(585,ng21)
  REAL(wp) :: kbc(5,5,13:59,ng21),absb(1175,ng21)
  REAL(wp) :: selfrefc(10,ng21)  ,forrefc(4,ng21)
  REAL(wp) :: sfluxrefc(ng21,9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))
  !EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13,1),ABSB(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta21

MODULE mo_yoesrta22

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA22* - SRTM COEFFICIENTS FOR INTERVAL 22
  !     BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng22 = 16

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(3,jpg)
  REAL(wp) :: sfluxref(jpg,9)
  REAL(wp) :: rayl            ,strrat
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng22) ,absa(585,ng22)
  REAL(wp) :: kbc(5,13:59,ng22),absb(235,ng22)
  REAL(wp) :: selfrefc(10,ng22),forrefc(3,ng22)
  REAL(wp) :: sfluxrefc(ng22,9)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta22

MODULE mo_yoesrta23

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA23* - SRTM COEFFICIENTS FOR INTERVAL 23
  !     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng23 = 16

  REAL(wp) :: ka(5,13,jpg)
  REAL(wp) :: selfref(10,jpg),forref(3,jpg)
  REAL(wp) :: sfluxref(jpg)  ,rayl(jpg)
  REAL(wp) :: givfac
  INTEGER  :: layreffr

  REAL(wp) :: kac(5,13,ng23)   ,absa(65,ng23)
  REAL(wp) :: selfrefc(10,ng23),forrefc(3,ng23)
  REAL(wp) :: sfluxrefc(ng23)  ,raylc(ng23)

  !EQUIVALENCE (KA(1,1,1),ABSA(1,1))
  EQUIVALENCE (kac(1,1,1),absa(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! GIVFAC  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! RAYLC   : REAL     Reduced g-point array for RAYL
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta23

MODULE mo_yoesrta24

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA24* - SRTM COEFFICIENTS FOR INTERVAL 24
  !     BAND 24: 12850-16000 cm-1 (low - H2O,O2; high - O2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng24 = 16

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(3,jpg)
  REAL(wp) :: sfluxref(jpg,9)
  REAL(wp) :: abso3a(jpg), abso3b(jpg), rayla(jpg,9), raylb(jpg)
  REAL(wp) :: strrat
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng24) ,absa(585,ng24)
  REAL(wp) :: kbc(5,13:59,ng24),absb(235,ng24)
  REAL(wp) :: selfrefc(10,ng24),forrefc(3,ng24)
  REAL(wp) :: sfluxrefc(ng24,9)
  REAL(wp) :: abso3ac(ng24), abso3bc(ng24), raylac(ng24,9), raylbc(ng24)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! ABSO3A  : REAL
  ! ABSO3B  : REAL
  ! RAYLA   : REAL
  ! RAYLB   : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! ABSO3AC : REAL     Reduced g-point array for ABSO3A
  ! ABSO3BC : REAL     Reduced g-point array for ABSO3B
  ! RAYLAC  : REAL     Reduced g-point array for RAYLA
  ! RAYLBC  : REAL     Reduced g-point array for RAYLB
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta24

MODULE mo_yoesrta25

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA25* - SRTM COEFFICIENTS FOR INTERVAL 25
  !     BAND 25: 16000-22650 cm-1 (low - H2O; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng25 = 16

  REAL(wp) :: ka(5,13,jpg)
  REAL(wp) :: sfluxref(jpg)
  REAL(wp) :: rayl(jpg), abso3a(jpg), abso3b(jpg)
  INTEGER  :: layreffr

  REAL(wp) :: kac(5,13,ng25) ,absa(65,ng25)
  REAL(wp) :: sfluxrefc(ng25)
  REAL(wp) :: raylc(ng25), abso3ac(ng25), abso3bc(ng25)

  !EQUIVALENCE (KA(1,1,1),ABSA(1,1))
  EQUIVALENCE (kac(1,1,1),absa(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! ABSO3A  : REAL
  ! ABSO3B  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! RAYLC   : REAL     Reduced g-point array for RAYL
  ! ABSO3AC : REAL     Reduced g-point array for ABSO3A
  ! ABSO3BC : REAL     Reduced g-point array for ABSO3B
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta25

MODULE mo_yoesrta26

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA26* - SRTM COEFFICIENTS FOR INTERVAL 26
  !     BAND 26: 22650-29000 cm-1 (low - nothing; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng26 = 16

  REAL(wp) :: sfluxref(jpg), rayl(jpg)

  REAL(wp) :: sfluxrefc(ng26), raylc(ng26)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! RAYLC   : REAL     Reduced g-point array for RAYL
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta26

MODULE mo_yoesrta27

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA27* - SRTM COEFFICIENTS FOR INTERVAL 27
  !     BAND 27: 29000-38000 cm-1 (low - O3; high - O3)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng27 = 16

  REAL(wp) :: ka(5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: sfluxref(jpg)  ,rayl(jpg)
  REAL(wp) :: scalekur
  INTEGER  :: layreffr

  REAL(wp) :: kac(5,13,ng27)   ,absa(65,ng27)
  REAL(wp) :: kbc(5,13:59,ng27),absb(235,ng27)
  REAL(wp) :: sfluxrefc(ng27)  ,raylc(ng27)

  !EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! FRACREFA: REAL
  ! KA      : REAL
  ! KB      : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! SCALEKUR: REAL
  ! LAYREFFR:INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! RAYLC   : REAL     Reduced g-point array for RAYL
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta27

MODULE mo_yoesrta28

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA28* - SRTM COEFFICIENTS FOR INTERVAL 28
  !     BAND 28: 38000-50000 cm-1 (low - O3, O2; high - O3, O2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng28 = 16

  REAL(wp) :: ka(9,5,13,jpg)
  REAL(wp) :: kb(5,5,13:59,jpg)
  REAL(wp) :: sfluxref(jpg,5)
  REAL(wp) :: rayl              ,strrat
  INTEGER  :: layreffr

  REAL(wp) :: kac(9,5,13,ng28)   ,absa(585,ng28)
  REAL(wp) :: kbc(5,5,13:59,ng28),absb(1175,ng28)
  REAL(wp) :: sfluxrefc(ng28,5)

  !EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1,1),absa(1,1)), (kbc(1,1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SFLUXREF: REAL
  ! RAYL    : REAL
  ! STRRAT  : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta28

MODULE mo_yoesrta29

  USE mo_kind, ONLY : wp

  IMPLICIT NONE

  PUBLIC

  SAVE

  !     -----------------------------------------------------------------
  !*    ** *YOESRTA29* - SRTM COEFFICIENTS FOR INTERVAL 29
  !     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: jpg = 16, ng29 = 16

  REAL(wp) :: ka(5,13,jpg)
  REAL(wp) :: kb(5,13:59,jpg)
  REAL(wp) :: selfref(10,jpg),forref(4,jpg)
  REAL(wp) :: sfluxref(jpg)  ,absh2o(jpg)  , absco2(jpg)
  REAL(wp) :: rayl
  INTEGER  :: layreffr

  REAL(wp) :: kac(5,13,ng29)   ,absa(65,ng29)
  REAL(wp) :: kbc(5,13:59,ng29),absb(235,ng29)
  REAL(wp) :: selfrefc(10,ng29),forrefc(4,ng29)
  REAL(wp) :: sfluxrefc(ng29)  ,absh2oc(ng29)  , absco2c(ng29)

  !EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
  EQUIVALENCE (kac(1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

  !     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
  !     M. J. IACONO          AER             12/09/03

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! KA      : REAL
  ! KB      : REAL
  ! SELFREF : REAL
  ! FORREF  : REAL
  ! SFLUXREF: REAL
  ! ABSH2O  : REAL
  ! ABSCO2  : REAL
  ! RAYL    : REAL
  ! LAYREFFR: INTEGER
  ! KAC     : REAL     Reduced g-point array for KA
  ! KBC     : REAL     Reduced g-point array for KB
  ! SELFREFC: REAL     Reduced g-point array for SELFREF
  ! FORREFC : REAL     Reduced g-point array for FORREF
  !SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
  ! ABSH2OC : REAL     Reduced g-point array for ABSH2O
  ! ABSCO2C : REAL     Reduced g-point array for ABSCO2
  !     -----------------------------------------------------------------
END MODULE mo_yoesrta29

