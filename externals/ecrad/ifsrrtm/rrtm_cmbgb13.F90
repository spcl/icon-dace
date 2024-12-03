! # 1 "ifsrrtm/rrtm_cmbgb13.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_cmbgb13.f90"
! this file has been modified for the use in icon

!***************************************************************************
subroutine rrtm_cmbgb13
!***************************************************************************

!     band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!     abozzo 201306 updated to rrtmg v4.85
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!***************************************************************************

! parameters
use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook, jphook

use yoerrto13, only : kao     ,selfrefo, forrefo   ,fracrefao, fracrefbo, &
                     & kao_mco2, kao_mco, kbo_mo3
use yoerrta13, only : ka      ,selfref, forref    ,fracrefa, fracrefb, &
                     & ka_mco2, ka_mco, kb_mo3
use yoerrtrwt, only : rwgt
use yoerrtftr, only : ngc      ,ngs      ,ngn      

implicit none

integer(kind=jpim) :: igc, ipr, iprsm, jn, jp, jt

real(kind=jprb) :: z_sumf, z_sumk, z_sumk1, z_sumk2
real(kind=jphook) :: zhook_handle


if (lhook) call dr_hook('rrtm_cmbgb13',0,zhook_handle)
do jn = 1,9
  do jt = 1,5
    do jp = 1,13
      iprsm = 0
      do igc = 1,ngc(13)
        z_sumk = 0.0_jprb
        do ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1

          z_sumk = z_sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
        enddo

        ka(jn,jt,jp,igc) = z_sumk
      enddo
    enddo
  enddo
enddo

do jn = 1,9
   do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(13)
        z_sumk1 = 0.0_jprb
        z_sumk2 = 0.0_jprb
         do ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            z_sumk1 = z_sumk1 + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+192)
            z_sumk2 = z_sumk2 + kao_mco(jn,jt,iprsm)*rwgt(iprsm+192)
         enddo
         ka_mco2(jn,jt,igc) = z_sumk1
         ka_mco(jn,jt,igc) = z_sumk2
      enddo
   enddo
enddo

do jt = 1,19
   iprsm = 0
   do igc = 1,ngc(13)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(12)+igc)
         iprsm = iprsm + 1
         z_sumk = z_sumk + kbo_mo3(jt,iprsm)*rwgt(iprsm+192)
      enddo
      kb_mo3(jt,igc) = z_sumk
   enddo
enddo


do jt = 1,10
  iprsm = 0
  do igc = 1,ngc(13)
    z_sumk = 0.0_jprb
    do ipr = 1, ngn(ngs(12)+igc)
      iprsm = iprsm + 1

      z_sumk = z_sumk + selfrefo(jt,iprsm)*rwgt(iprsm+192)
    enddo

    selfref(jt,igc) = z_sumk
  enddo
enddo

do jt = 1,4
   iprsm = 0
   do igc = 1,ngc(13)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(12)+igc)
         iprsm = iprsm + 1
         z_sumk = z_sumk + forrefo(jt,iprsm)*rwgt(iprsm+192)
      enddo
      forref(jt,igc) = z_sumk
   enddo
enddo

iprsm = 0
do igc = 1,ngc(13)
   z_sumf = 0.0_jprb
   do ipr = 1, ngn(ngs(12)+igc)
      iprsm = iprsm + 1
      z_sumf = z_sumf + fracrefbo(iprsm)
   enddo
   fracrefb(igc) = z_sumf
enddo

do jp = 1,9
  iprsm = 0
  do igc = 1,ngc(13)
    z_sumf = 0.0_jprb
    do ipr = 1, ngn(ngs(12)+igc)
      iprsm = iprsm + 1

      z_sumf = z_sumf + fracrefao(iprsm,jp)
    enddo

    fracrefa(igc,jp) = z_sumf
  enddo
enddo


!$acc update device(fracrefa, fracrefb, ka, selfref, forref, ka_mco2, ka_mco, &
!$acc               kb_mo3)

if (lhook) call dr_hook('rrtm_cmbgb13',1,zhook_handle)
end subroutine rrtm_cmbgb13
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

