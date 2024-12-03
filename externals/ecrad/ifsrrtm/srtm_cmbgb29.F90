! # 1 "ifsrrtm/srtm_cmbgb29.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_cmbgb29.f90"
! this file has been modified for the use in icon

subroutine srtm_cmbgb29

!     band 29:  820-2600 cm-1 (low - h2o; high - co2)
!-----------------------------------------------------------------------

use parkind1  ,only : jpim , jprb
use ecradhook   ,only : lhook, dr_hook, jphook

use yoesrtm  , only : ngn
use yoesrtwn , only : ngc, ngs, rwgt
!use yoesrtwn , only : ngc, ngs, ngn, rwgt
use yoesrta29, only : ka, kb, selfref, forref, sfluxref, &
                    & absh2o, absco2, &
                    & kac, kbc, selfrefc, forrefc, sfluxrefc, &
                    & absh2oc, absco2c

implicit none

! local variables
integer(kind=jpim) :: jt, jp, igc, ipr, iprsm
real(kind=jprb)    :: zsumk, zsumf1, zsumf2, zsumf3

real(kind=jphook) :: zhook_handle
!     ------------------------------------------------------------------
if (lhook) call dr_hook('srtm_cmbgb29',0,zhook_handle)

do jt = 1,5
  do jp = 1,13
    iprsm = 0
    do igc = 1,ngc(14)
      zsumk = 0.
      do ipr = 1, ngn(ngs(13)+igc)
        iprsm = iprsm + 1
        zsumk = zsumk + ka(jt,jp,iprsm)*rwgt(iprsm+208)
      enddo
      kac(jt,jp,igc) = zsumk
    enddo
  enddo

  do jp = 13,59
    iprsm = 0
    do igc = 1,ngc(14)
      zsumk = 0.
      do ipr = 1, ngn(ngs(13)+igc)
        iprsm = iprsm + 1
        zsumk = zsumk + kb(jt,jp,iprsm)*rwgt(iprsm+208)
      enddo
      kbc(jt,jp,igc) = zsumk
    enddo
  enddo
enddo

do jt = 1,10
  iprsm = 0
  do igc = 1,ngc(14)
    zsumk = 0.
    do ipr = 1, ngn(ngs(13)+igc)
      iprsm = iprsm + 1
      zsumk = zsumk + selfref(jt,iprsm)*rwgt(iprsm+208)
    enddo
    selfrefc(jt,igc) = zsumk
  enddo
enddo

do jt = 1,4
  iprsm = 0
  do igc = 1,ngc(14)
    zsumk = 0.
    do ipr = 1, ngn(ngs(13)+igc)
      iprsm = iprsm + 1
      zsumk = zsumk + forref(jt,iprsm)*rwgt(iprsm+208)
    enddo
    forrefc(jt,igc) = zsumk
  enddo
enddo

iprsm = 0
do igc = 1,ngc(14)
  zsumf1 = 0.
  zsumf2 = 0.
  zsumf3 = 0.
  do ipr = 1, ngn(ngs(13)+igc)
    iprsm = iprsm + 1
    zsumf1 = zsumf1 + sfluxref(iprsm)
    zsumf2 = zsumf2 + absco2(iprsm)*rwgt(iprsm+208)
    zsumf3 = zsumf3 + absh2o(iprsm)*rwgt(iprsm+208)
  enddo
  sfluxrefc(igc) = zsumf1
  absco2c(igc) = zsumf2
  absh2oc(igc) = zsumf3
enddo

!$acc update device(kac, kbc, selfrefc, forrefc, sfluxrefc, absh2oc, absco2c)

!     -----------------------------------------------------------------
if (lhook) call dr_hook('srtm_cmbgb29',1,zhook_handle)
end subroutine srtm_cmbgb29

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

