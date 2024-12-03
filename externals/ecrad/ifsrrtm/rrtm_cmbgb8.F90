! # 1 "ifsrrtm/rrtm_cmbgb8.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_cmbgb8.f90"
! this file has been modified for the use in icon

!***************************************************************************
subroutine rrtm_cmbgb8
!***************************************************************************

!     band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!     abozzo 201306 updated to rrtmg v4.85
!***************************************************************************

! parameters
use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook, jphook

use yoerrto8 , only : kao     ,kbo     ,selfrefo   ,forrefo, fracrefao  ,&
 & fracrefbo, kao_mco2, kao_mn2o ,kao_mo3, kbo_mco2, kbo_mn2o, &
 & cfc12o   , cfc22adjo  
use yoerrta8 , only : ka      ,kb      ,selfref    ,forref, fracrefa   ,&
 & fracrefb , ka_mco2, ka_mn2o ,ka_mo3, kb_mco2, kb_mn2o,&
 & cfc12    , cfc22adj  
use yoerrtrwt, only : rwgt
use yoerrtftr, only : ngc      ,ngs      ,ngn      

implicit none

integer(kind=jpim) :: igc, ipr, iprsm, jp, jt

real(kind=jprb) :: z_sumf1, z_sumf2, z_sumk, z_sumk1, z_sumk2, z_sumk3, z_sumk4, z_sumk5
real(kind=jphook) :: zhook_handle

if (lhook) call dr_hook('rrtm_cmbgb8',0,zhook_handle)
do jt = 1,5
  do jp = 1,13
    iprsm = 0
    do igc = 1,ngc(8)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(7)+igc)
        iprsm = iprsm + 1
        z_sumk = z_sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
      enddo
      ka(jt,jp,igc) = z_sumk
    enddo
  enddo
enddo
do jt = 1,5
  do jp = 13,59
    iprsm = 0
    do igc = 1,ngc(8)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(7)+igc)
        iprsm = iprsm + 1
        z_sumk = z_sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+112)
      enddo
      kb(jt,jp,igc) = z_sumk
    enddo
  enddo
enddo

do jt = 1,10
  iprsm = 0
  do igc = 1,ngc(8)
    z_sumk = 0.0_jprb
    do ipr = 1, ngn(ngs(7)+igc)
      iprsm = iprsm + 1
      z_sumk = z_sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
    enddo
    selfref(jt,igc) = z_sumk
  enddo
enddo
do jt = 1,4
   iprsm = 0
   do igc = 1,ngc(8)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(7)+igc)
         iprsm = iprsm + 1
         z_sumk = z_sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
      enddo
      forref(jt,igc) = z_sumk
   enddo
enddo

do jt = 1,19
iprsm = 0
do igc = 1,ngc(8)
  z_sumk1= 0.0_jprb
  z_sumk2= 0.0_jprb
  z_sumk3= 0.0_jprb
  z_sumk4= 0.0_jprb
  z_sumk5= 0.0_jprb
  do ipr = 1, ngn(ngs(7)+igc)
    iprsm = iprsm + 1
    z_sumk1= z_sumk1+ kao_mco2(jt,iprsm)*rwgt(iprsm+112)
    z_sumk2= z_sumk2+ kbo_mco2(jt,iprsm)*rwgt(iprsm+112)
    z_sumk3= z_sumk3+ kao_mo3(jt,iprsm)*rwgt(iprsm+112)
    z_sumk4= z_sumk4+ kao_mn2o(jt,iprsm)*rwgt(iprsm+112)
    z_sumk5= z_sumk5+ kbo_mn2o(jt,iprsm)*rwgt(iprsm+112)
  enddo
  ka_mco2(jt,igc) = z_sumk1
  kb_mco2(jt,igc) = z_sumk2
  ka_mo3(jt,igc) = z_sumk3
  ka_mn2o(jt,igc) = z_sumk4
  kb_mn2o(jt,igc) = z_sumk5
enddo
enddo



iprsm = 0
do igc = 1,ngc(8)
  z_sumf1= 0.0_jprb
  z_sumf2= 0.0_jprb
  z_sumk1= 0.0_jprb
  z_sumk2= 0.0_jprb
  do ipr = 1, ngn(ngs(7)+igc)
    iprsm = iprsm + 1
    z_sumf1= z_sumf1+ fracrefao(iprsm)
    z_sumf2= z_sumf2+ fracrefbo(iprsm)
    z_sumk1= z_sumk1+ cfc12o(iprsm)*rwgt(iprsm+112)
    z_sumk2= z_sumk2+ cfc22adjo(iprsm)*rwgt(iprsm+112)
  enddo
  fracrefa(igc) = z_sumf1
  fracrefb(igc) = z_sumf2
  cfc12(igc) = z_sumk1
  cfc22adj(igc) = z_sumk2
enddo


!$acc update device(fracrefa, fracrefb, cfc12, cfc22adj, ka, kb, ka_mco2, &
!$acc               ka_mn2o, ka_mo3, kb_mco2, kb_mn2o, selfref, forref)

if (lhook) call dr_hook('rrtm_cmbgb8',1,zhook_handle)
end subroutine rrtm_cmbgb8
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define __sizeof_int__ 4
! #define __sizeof_pointer__ 8
! #define __gfortran__ 1
! #define __gfc_real_16__ 1
! #define __stdc_hosted__ 0
! #define __no_math_errno__ 1
! #define __sizeof_float__ 4
! #define __pic__ 2
! #define _language_fortran 1
! #define __sizeof_long__ 8
! #define __gfc_int_8__ 1
! #define __dynamic__ 1
! #define __sizeof_short__ 2
! #define __gnuc__ 13
! #define __sizeof_long_double__ 16
! #define __biggest_alignment__ 16
! #define __atomic_relaxed 0
! #define _lp64 1
! #define __ecrad_little_endian 1
! #define __gfc_int_1__ 1
! #define __order_big_endian__ 4321
! #define __byte_order__ __order_little_endian__
! #define __sizeof_size_t__ 8
! #define __pic__ 2
! #define __sizeof_double__ 8
! #define __atomic_consume 1
! #define __gnuc_minor__ 3
! #define __gfc_int_16__ 1
! #define __lp64__ 1
! #define __atomic_seq_cst 5
! #define __sizeof_long_long__ 8
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

