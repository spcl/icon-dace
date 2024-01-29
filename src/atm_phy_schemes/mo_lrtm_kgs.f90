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


MODULE mo_rrlw_kg01

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 1
  ! band 1:  10-250 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mn2 : real
  ! kbo_mn2 : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no1  = 16

  REAL(wp) :: fracrefao(no1)  , fracrefbo(no1)
  REAL(wp) :: kao(5,13,no1)
  REAL(wp) :: kbo(5,13:59,no1)
  REAL(wp) :: kao_mn2(19,no1) , kbo_mn2(19,no1)
  REAL(wp) :: selfrefo(10,no1), forrefo(4,no1)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 1
  ! band 1:  10-250 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! absa    : real
  ! absb    : real
  ! ka_mn2  : real
  ! kb_mn2  : real
  ! selfref : real
  ! forref  : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng1  = 10

  REAL(wp) :: fracrefa(ng1)  , fracrefb(ng1)
  REAL(wp) :: ka(5,13,ng1)   , absa(65,ng1)
  REAL(wp) :: kb(5,13:59,ng1), absb(235,ng1)
  REAL(wp) :: ka_mn2(19,ng1) , kb_mn2(19,ng1)
  REAL(wp) :: selfref(10,ng1), forref(4,ng1)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg01

MODULE mo_rrlw_kg02

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 2
  ! band 2:  250-500 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no2  = 16

  REAL(wp) :: fracrefao(no2)   , fracrefbo(no2)
  REAL(wp) :: kao(5,13,no2)
  REAL(wp) :: kbo(5,13:59,no2)
  REAL(wp) :: selfrefo(10,no2) , forrefo(4,no2)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 2
  ! band 2:  250-500 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! absa    : real
  ! absb    : real
  ! selfref : real
  ! forref  : real
  !
  ! refparam: real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng2  = 12

  REAL(wp) :: fracrefa(ng2)  , fracrefb(ng2)
  REAL(wp) :: ka(5,13,ng2)   , absa(65,ng2)
  REAL(wp) :: kb(5,13:59,ng2), absb(235,ng2)
  REAL(wp) :: selfref(10,ng2), forref(4,ng2)

  REAL(wp) :: refparam(13)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg02

MODULE mo_rrlw_kg03

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 3
  ! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mn2o: real
  ! kbo_mn2o: real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no3  = 16

  REAL(wp) :: fracrefao(no3,9) ,fracrefbo(no3,5)
  REAL(wp) :: kao(9,5,13,no3)
  REAL(wp) :: kbo(5,5,13:59,no3)
  REAL(wp) :: kao_mn2o(9,19,no3), kbo_mn2o(5,19,no3)
  REAL(wp) :: selfrefo(10,no3)
  REAL(wp) :: forrefo(4,no3)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 3
  ! band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! ka_mn2o : real
  ! kb_mn2o : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng3  = 16

  REAL(wp) :: fracrefa(ng3,9) ,fracrefb(ng3,5)
  REAL(wp) :: ka(9,5,13,ng3)  ,absa(585,ng3)
  REAL(wp) :: kb(5,5,13:59,ng3),absb(1175,ng3)
  REAL(wp) :: ka_mn2o(9,19,ng3), kb_mn2o(5,19,ng3)
  REAL(wp) :: selfref(10,ng3)
  REAL(wp) :: forref(4,ng3)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

END MODULE mo_rrlw_kg03

MODULE mo_rrlw_kg04

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 4
  ! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no4  = 16

  REAL(wp) :: fracrefao(no4,9)  ,fracrefbo(no4,5)
  REAL(wp) :: kao(9,5,13,no4)
  REAL(wp) :: kbo(5,5,13:59,no4)
  REAL(wp) :: selfrefo(10,no4)  ,forrefo(4,no4)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 4
  ! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  ! absa    : real
  ! absb    : real
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! selfref : real
  ! forref  : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng4  = 14

  REAL(wp) :: fracrefa(ng4,9)  ,fracrefb(ng4,5)
  REAL(wp) :: ka(9,5,13,ng4)   ,absa(585,ng4)
  REAL(wp) :: kb(5,5,13:59,ng4),absb(1175,ng4)
  REAL(wp) :: selfref(10,ng4)  ,forref(4,ng4)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

END MODULE mo_rrlw_kg04

MODULE mo_rrlw_kg05

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 5
  ! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mo3 : real
  ! selfrefo: real
  ! forrefo : real
  ! ccl4o   : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no5  = 16

  REAL(wp) :: fracrefao(no5,9) ,fracrefbo(no5,5)
  REAL(wp) :: kao(9,5,13,no5)
  REAL(wp) :: kbo(5,5,13:59,no5)
  REAL(wp) :: kao_mo3(9,19,no5)
  REAL(wp) :: selfrefo(10,no5)
  REAL(wp) :: forrefo(4,no5)
  REAL(wp) :: ccl4o(no5)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 5
  ! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! ka_mo3  : real
  ! selfref : real
  ! forref  : real
  ! ccl4    : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng5  = 16

  REAL(wp) :: fracrefa(ng5,9) ,fracrefb(ng5,5)
  REAL(wp) :: ka(9,5,13,ng5)   ,absa(585,ng5)
  REAL(wp) :: kb(5,5,13:59,ng5),absb(1175,ng5)
  REAL(wp) :: ka_mo3(9,19,ng5)
  REAL(wp) :: selfref(10,ng5)
  REAL(wp) :: forref(4,ng5)
  REAL(wp) :: ccl4(ng5)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

END MODULE mo_rrlw_kg05

MODULE mo_rrlw_kg06

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 6
  ! band 6:  820-980 cm-1 (low - h2o; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  ! kao     : real
  ! kao_mco2: real
  ! selfrefo: real
  ! forrefo : real
  !cfc11adjo: real
  ! cfc12o  : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no6  = 16

  REAL(wp) , DIMENSION(no6) :: fracrefao
  REAL(wp) :: kao(5,13,no6)
  REAL(wp) :: kao_mco2(19,no6)
  REAL(wp) :: selfrefo(10,no6)
  REAL(wp) :: forrefo(4,no6)

  REAL(wp) , DIMENSION(no6) :: cfc11adjo
  REAL(wp) , DIMENSION(no6) :: cfc12o

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 6
  ! band 6:  820-980 cm-1 (low - h2o; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  ! ka      : real
  ! ka_mco2 : real
  ! selfref : real
  ! forref  : real
  !cfc11adj : real
  ! cfc12   : real
  !
  ! absa    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng6  = 8

  REAL(wp) , DIMENSION(ng6) :: fracrefa
  REAL(wp) :: ka(5,13,ng6),absa(65,ng6)
  REAL(wp) :: ka_mco2(19,ng6)
  REAL(wp) :: selfref(10,ng6)
  REAL(wp) :: forref(4,ng6)

  REAL(wp) , DIMENSION(ng6) :: cfc11adj
  REAL(wp) , DIMENSION(ng6) :: cfc12

  EQUIVALENCE (ka(1,1,1),absa(1,1))

END MODULE mo_rrlw_kg06

MODULE mo_rrlw_kg07

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 7
  ! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mco2: real
  ! kbo_mco2: real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no7  = 16

  REAL(wp) , DIMENSION(no7) :: fracrefbo
  REAL(wp) :: fracrefao(no7,9)
  REAL(wp) :: kao(9,5,13,no7)
  REAL(wp) :: kbo(5,13:59,no7)
  REAL(wp) :: kao_mco2(9,19,no7)
  REAL(wp) :: kbo_mco2(19,no7)
  REAL(wp) :: selfrefo(10,no7)
  REAL(wp) :: forrefo(4,no7)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 7
  ! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! ka_mco2 : real
  ! kb_mco2 : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng7  = 12

  REAL(wp) , DIMENSION(ng7) :: fracrefb
  REAL(wp) :: fracrefa(ng7,9)
  REAL(wp) :: ka(9,5,13,ng7) ,absa(585,ng7)
  REAL(wp) :: kb(5,13:59,ng7),absb(235,ng7)
  REAL(wp) :: ka_mco2(9,19,ng7)
  REAL(wp) :: kb_mco2(19,ng7)
  REAL(wp) :: selfref(10,ng7)
  REAL(wp) :: forref(4,ng7)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg07

MODULE mo_rrlw_kg08

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 8
  ! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mco2: real
  ! kbo_mco2: real
  ! kao_mn2o: real
  ! kbo_mn2o: real
  ! kao_mo3 : real
  ! selfrefo: real
  ! forrefo : real
  ! cfc12o  : real
  !cfc22adjo: real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no8  = 16

  REAL(wp) , DIMENSION(no8) :: fracrefao
  REAL(wp) , DIMENSION(no8) :: fracrefbo
  REAL(wp) , DIMENSION(no8) :: cfc12o
  REAL(wp) , DIMENSION(no8) :: cfc22adjo

  REAL(wp) :: kao(5,13,no8)
  REAL(wp) :: kao_mco2(19,no8)
  REAL(wp) :: kao_mn2o(19,no8)
  REAL(wp) :: kao_mo3(19,no8)
  REAL(wp) :: kbo(5,13:59,no8)
  REAL(wp) :: kbo_mco2(19,no8)
  REAL(wp) :: kbo_mn2o(19,no8)
  REAL(wp) :: selfrefo(10,no8)
  REAL(wp) :: forrefo(4,no8)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 8
  ! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! ka_mco2 : real
  ! kb_mco2 : real
  ! ka_mn2o : real
  ! kb_mn2o : real
  ! ka_mo3  : real
  ! selfref : real
  ! forref  : real
  ! cfc12   : real
  ! cfc22adj: real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng8  = 8

  REAL(wp) , DIMENSION(ng8) :: fracrefa
  REAL(wp) , DIMENSION(ng8) :: fracrefb
  REAL(wp) , DIMENSION(ng8) :: cfc12
  REAL(wp) , DIMENSION(ng8) :: cfc22adj

  REAL(wp) :: ka(5,13,ng8)    ,absa(65,ng8)
  REAL(wp) :: kb(5,13:59,ng8) ,absb(235,ng8)
  REAL(wp) :: ka_mco2(19,ng8)
  REAL(wp) :: ka_mn2o(19,ng8)
  REAL(wp) :: ka_mo3(19,ng8)
  REAL(wp) :: kb_mco2(19,ng8)
  REAL(wp) :: kb_mn2o(19,ng8)
  REAL(wp) :: selfref(10,ng8)
  REAL(wp) :: forref(4,ng8)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg08

MODULE mo_rrlw_kg09

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 9
  ! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mn2o: real
  ! kbo_mn2o: real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no9  = 16

  REAL(wp) , DIMENSION(no9) :: fracrefbo

  REAL(wp) :: fracrefao(no9,9)
  REAL(wp) :: kao(9,5,13,no9)
  REAL(wp) :: kbo(5,13:59,no9)
  REAL(wp) :: kao_mn2o(9,19,no9)
  REAL(wp) :: kbo_mn2o(19,no9)
  REAL(wp) :: selfrefo(10,no9)
  REAL(wp) :: forrefo(4,no9)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 9
  ! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! ka_mn2o : real
  ! kb_mn2o : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng9  = 12

  REAL(wp) , DIMENSION(ng9) :: fracrefb
  REAL(wp) :: fracrefa(ng9,9)
  REAL(wp) :: ka(9,5,13,ng9) ,absa(585,ng9)
  REAL(wp) :: kb(5,13:59,ng9) ,absb(235,ng9)
  REAL(wp) :: ka_mn2o(9,19,ng9)
  REAL(wp) :: kb_mn2o(19,ng9)
  REAL(wp) :: selfref(10,ng9)
  REAL(wp) :: forref(4,ng9)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg09

MODULE mo_rrlw_kg10

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 10
  ! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no10 = 16

  REAL(wp) , DIMENSION(no10) :: fracrefao
  REAL(wp) , DIMENSION(no10) :: fracrefbo

  REAL(wp) :: kao(5,13,no10)
  REAL(wp) :: kbo(5,13:59,no10)
  REAL(wp) :: selfrefo(10,no10)
  REAL(wp) :: forrefo(4,no10)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 10
  ! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng10 = 6

  REAL(wp) , DIMENSION(ng10) :: fracrefa
  REAL(wp) , DIMENSION(ng10) :: fracrefb

  REAL(wp) :: ka(5,13,ng10)   , absa(65,ng10)
  REAL(wp) :: kb(5,13:59,ng10), absb(235,ng10)
  REAL(wp) :: selfref(10,ng10)
  REAL(wp) :: forref(4,ng10)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg10

MODULE mo_rrlw_kg11

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 11
  ! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! kao_mo2 : real
  ! kbo_mo2 : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no11 = 16

  REAL(wp) , DIMENSION(no11) :: fracrefao
  REAL(wp) , DIMENSION(no11) :: fracrefbo

  REAL(wp) :: kao(5,13,no11)
  REAL(wp) :: kbo(5,13:59,no11)
  REAL(wp) :: kao_mo2(19,no11)
  REAL(wp) :: kbo_mo2(19,no11)
  REAL(wp) :: selfrefo(10,no11)
  REAL(wp) :: forrefo(4,no11)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 11
  ! band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! ka_mo2  : real
  ! kb_mo2  : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng11 = 8

  REAL(wp) , DIMENSION(ng11) :: fracrefa
  REAL(wp) , DIMENSION(ng11) :: fracrefb

  REAL(wp) :: ka(5,13,ng11)   , absa(65,ng11)
  REAL(wp) :: kb(5,13:59,ng11), absb(235,ng11)
  REAL(wp) :: ka_mo2(19,ng11)
  REAL(wp) :: kb_mo2(19,ng11)
  REAL(wp) :: selfref(10,ng11)
  REAL(wp) :: forref(4,ng11)

  EQUIVALENCE (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg11

MODULE mo_rrlw_kg12

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 12
  ! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  ! kao     : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no12 = 16

  REAL(wp) :: fracrefao(no12,9)
  REAL(wp) :: kao(9,5,13,no12)
  REAL(wp) :: selfrefo(10,no12)
  REAL(wp) :: forrefo(4,no12)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 12
  ! band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  ! ka      : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng12 = 8

  REAL(wp) :: fracrefa(ng12,9)
  REAL(wp) :: ka(9,5,13,ng12) ,absa(585,ng12)
  REAL(wp) :: selfref(10,ng12)
  REAL(wp) :: forref(4,ng12)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1))

END MODULE mo_rrlw_kg12

MODULE mo_rrlw_kg13

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 13
  ! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  ! kao     : real
  ! kao_mco2: real
  ! kao_mco : real
  ! kbo_mo3 : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no13 = 16

  REAL(wp) , DIMENSION(no13) :: fracrefbo

  REAL(wp) :: fracrefao(no13,9)
  REAL(wp) :: kao(9,5,13,no13)
  REAL(wp) :: kao_mco2(9,19,no13)
  REAL(wp) :: kao_mco(9,19,no13)
  REAL(wp) :: kbo_mo3(19,no13)
  REAL(wp) :: selfrefo(10,no13)
  REAL(wp) :: forrefo(4,no13)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 13
  ! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  ! ka      : real
  ! ka_mco2 : real
  ! ka_mco  : real
  ! kb_mo3  : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng13 = 4

  REAL(wp) , DIMENSION(ng13) :: fracrefb

  REAL(wp) :: fracrefa(ng13,9)
  REAL(wp) :: ka(9,5,13,ng13) ,absa(585,ng13)
  REAL(wp) :: ka_mco2(9,19,ng13)
  REAL(wp) :: ka_mco(9,19,ng13)
  REAL(wp) :: kb_mo3(19,ng13)
  REAL(wp) :: selfref(10,ng13)
  REAL(wp) :: forref(4,ng13)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1))

END MODULE mo_rrlw_kg13

MODULE mo_rrlw_kg14

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 14
  ! band 14:  2250-2380 cm-1 (low - co2; high - co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  !fracrefbo: real
  ! kao     : real
  ! kbo     : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no14 = 16

  REAL(wp) , DIMENSION(no14) :: fracrefao
  REAL(wp) , DIMENSION(no14) :: fracrefbo

  REAL(wp) :: kao(5,13,no14)
  REAL(wp) :: kbo(5,13:59,no14)
  REAL(wp) :: selfrefo(10,no14)
  REAL(wp) :: forrefo(4,no14)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 14
  ! band 14:  2250-2380 cm-1 (low - co2; high - co2)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  !fracrefb : real
  ! ka      : real
  ! kb      : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng14 = 2

  REAL(wp) , DIMENSION(ng14) :: fracrefa
  REAL(wp) , DIMENSION(ng14) :: fracrefb

  REAL(wp) :: ka(5,13,ng14)   ,absa(65,ng14)
  REAL(wp) :: kb(5,13:59,ng14),absb(235,ng14)
  REAL(wp) :: selfref(10,ng14)
  REAL(wp) :: forref(4,ng14)

  EQUIVALENCE (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg14

MODULE mo_rrlw_kg15

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 15
  ! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  ! kao     : real
  ! kao_mn2 : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no15 = 16

  REAL(wp) :: fracrefao(no15,9)
  REAL(wp) :: kao(9,5,13,no15)
  REAL(wp) :: kao_mn2(9,19,no15)
  REAL(wp) :: selfrefo(10,no15)
  REAL(wp) :: forrefo(4,no15)


  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 15
  ! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  ! ka      : real
  ! ka_mn2  : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng15 = 2

  REAL(wp) :: fracrefa(ng15,9)
  REAL(wp) :: ka(9,5,13,ng15) ,absa(585,ng15)
  REAL(wp) :: ka_mn2(9,19,ng15)
  REAL(wp) :: selfref(10,ng15)
  REAL(wp) :: forref(4,ng15)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1))

END MODULE mo_rrlw_kg15
MODULE mo_rrlw_kg16

  USE mo_kind ,ONLY : wp

  IMPLICIT NONE

  PUBLIC

  !-----------------------------------------------------------------
  ! rrtmg_lw ORIGINAL abs. coefficients for interval 16
  ! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefao: real
  ! kao     : real
  ! kbo     : real
  ! selfrefo: real
  ! forrefo : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: no16 = 16

  REAL(wp) , DIMENSION(no16) :: fracrefbo

  REAL(wp) :: fracrefao(no16,9)
  REAL(wp) :: kao(9,5,13,no16)
  REAL(wp) :: kbo(5,13:59,no16)
  REAL(wp) :: selfrefo(10,no16)
  REAL(wp) :: forrefo(4,no16)

  !-----------------------------------------------------------------
  ! rrtmg_lw COMBINED abs. coefficients for interval 16
  ! band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
  !
  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !-----------------------------------------------------------------
  !
  !  name     type     purpose
  !  ----   : ----   : ---------------------------------------------
  !fracrefa : real
  ! ka      : real
  ! kb      : real
  ! selfref : real
  ! forref  : real
  !
  ! absa    : real
  ! absb    : real
  !-----------------------------------------------------------------

  INTEGER, PARAMETER :: ng16 = 2

  REAL(wp) , DIMENSION(ng16) :: fracrefb

  REAL(wp) :: fracrefa(ng16,9)
  REAL(wp) :: ka(9,5,13,ng16) ,absa(585,ng16)
  REAL(wp) :: kb(5,13:59,ng16), absb(235,ng16)
  REAL(wp) :: selfref(10,ng16)
  REAL(wp) :: forref(4,ng16)

  EQUIVALENCE (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

END MODULE mo_rrlw_kg16
