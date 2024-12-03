! # 1 "ifsrrtm/rrtm_cmbgb1.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_cmbgb1.f90"
! this file has been modified for the use in icon

!***************************************************************************
subroutine rrtm_cmbgb1
!***************************************************************************

!  the subroutines cmbgb1->cmbgb16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. the data are combined with appropriate weighting following the
!  g-point mapping arrays specified in rrtminit.  plank fraction data
!  in arrays fracrefa and fracrefb are combined without weighting.  all
!  g-point reduced data are put into new arrays for use in rrtm.

!  band 1:  10-250 cm-1 (low - h2o; high - h2o)

! abozzo may 2013 update to the last version of rrtmg

!band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                       (high key - h2o; high minor - n2)
!***************************************************************************

! parameters
use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook, jphook

use yoerrto1 , only : kao, kbo, selfrefo, forrefo, fracrefao,fracrefbo,kao_mn2, kbo_mn2
use yoerrta1 , only : ka , kb , selfref , forref , fracrefa ,fracrefb,  ka_mn2, kb_mn2
use yoerrtrwt, only : rwgt
use yoerrtftr, only : ngc      ,ngn      

implicit none

integer(kind=jpim) :: igc, ipr, iprsm, jp, jt

real(kind=jprb) :: z_sumf1, z_sumf2, z_sumk, z_sumk1, z_sumk2
real(kind=jphook) :: zhook_handle

if (lhook) call dr_hook('rrtm_cmbgb1',0,zhook_handle)
do jt = 1,5
  do jp = 1,13
    iprsm = 0
    do igc = 1,ngc(1)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(igc)
        iprsm = iprsm + 1

        z_sumk = z_sumk + kao(jt,jp,iprsm)*rwgt(iprsm)
      enddo

      ka(jt,jp,igc) = z_sumk
    enddo
  enddo
  do jp = 13,59
    iprsm = 0
    do igc = 1,ngc(1)
      z_sumk = 0.0_jprb
      do ipr = 1, ngn(igc)
        iprsm = iprsm + 1

        z_sumk = z_sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
      enddo

      kb(jt,jp,igc) = z_sumk
    enddo
  enddo
enddo

do jt = 1,10
  iprsm = 0
  do igc = 1,ngc(1)
    z_sumk = 0.0_jprb
    do ipr = 1, ngn(igc)
      iprsm = iprsm + 1

      z_sumk = z_sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
    enddo

    selfref(jt,igc) = z_sumk
  enddo
enddo

do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(1)
            z_sumk = 0.0_jprb
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               z_sumk = z_sumk + forrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            forref(jt,igc) = z_sumk
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(1)
            z_sumk1 = 0.0_jprb
            z_sumk2 = 0.0_jprb
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               z_sumk1 = z_sumk1 + kao_mn2(jt,iprsm)*rwgt(iprsm)
               z_sumk2 = z_sumk2 + kbo_mn2(jt,iprsm)*rwgt(iprsm)
            enddo
            ka_mn2(jt,igc) = z_sumk1
            kb_mn2(jt,igc) = z_sumk2
         enddo
      enddo

iprsm = 0
do igc = 1,ngc(1)
  z_sumf1 = 0.0_jprb
  z_sumf2 = 0.0_jprb
  do ipr = 1, ngn(igc)
    iprsm = iprsm + 1

    z_sumf1= z_sumf1+ fracrefao(iprsm)
    z_sumf2= z_sumf2+ fracrefbo(iprsm)
  enddo

  fracrefa(igc) = z_sumf1
  fracrefb(igc) = z_sumf2
enddo

!$acc update device(fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, &
!$acc               selfref, forref)

if (lhook) call dr_hook('rrtm_cmbgb1',1,zhook_handle)
end subroutine rrtm_cmbgb1
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

