! # 1 "ifsrrtm/rrtm_cmbgb10.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_cmbgb10.f90"
! this file has been modified for the use in icon

!***************************************************************************
subroutine rrtm_cmbgb10
!***************************************************************************

!     band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!     abozzo updated to rrtmg v4.85
!***************************************************************************

! parameters
use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook, jphook

use yoerrto10, only : kao     ,kbo      ,fracrefao   ,fracrefbo, selfrefo,forrefo
use yoerrta10, only : ka      ,kb       ,fracrefa    ,fracrefb, selfref,forref
use yoerrtrwt, only : rwgt
use yoerrtftr, only : ngc      ,ngs      ,ngn      

implicit none

integer(kind=jpim) :: igc, ipr, iprsm, jp, jt

real(kind=jprb) :: z_sumf1, z_sumf2, z_sumk
real(kind=jphook) :: zhook_handle

if (lhook) call dr_hook('rrtm_cmbgb10',0,zhook_handle)
do jt = 1,5
  do jp = 1,13
    iprsm = 0
    do igc = 1,ngc(10)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(9)+igc)
        iprsm = iprsm + 1

        z_sumk = z_sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
      enddo

      ka(jt,jp,igc) = z_sumk
    enddo
  enddo
enddo
do jt = 1,5
  do jp = 13,59
    iprsm = 0
    do igc = 1,ngc(10)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(ngs(9)+igc)
        iprsm = iprsm + 1

        z_sumk = z_sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+144)
      enddo

      kb(jt,jp,igc) = z_sumk
    enddo
  enddo
enddo

  do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(10)
            z_sumk = 0.0_jprb
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               z_sumk = z_sumk + selfrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            selfref(jt,igc) = z_sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(10)
            z_sumk = 0.0_jprb
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               z_sumk = z_sumk + forrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            forref(jt,igc) = z_sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(10)
         z_sumf1= 0.0_jprb
         z_sumf2= 0.0_jprb
         do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            z_sumf1= z_sumf1+ fracrefao(iprsm)
            z_sumf2= z_sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = z_sumf1
         fracrefb(igc) = z_sumf2
      enddo


!$acc update device(fracrefa, fracrefb, ka, kb, selfref, forref)

if (lhook) call dr_hook('rrtm_cmbgb10',1,zhook_handle)
end subroutine rrtm_cmbgb10
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

