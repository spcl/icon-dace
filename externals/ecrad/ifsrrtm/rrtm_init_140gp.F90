! # 1 "ifsrrtm/rrtm_init_140gp.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_init_140gp.f90"
!***************************************************************************
subroutine rrtm_init_140gp(cdirectory)
!***************************************************************************
!     reformatted for f90 by jjmorcrette, ecmwf, 980714

!     jjmorcrette 20110613 flexible number of g-points

! parameters
use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook, dr_hook, jphook

use parrrtm  , only : jpband   ,jpg
use yoerrtm  , only : jpgpt
use yoerrtwn , only : ng       
use yoerrtftr, only : ngc      ,ngn      ,ngm     , wt
! output
use yoerrtbg2, only : corr1    ,corr2
use yoerrtrwt, only : frefa    ,frefb    ,frefadf  ,frefbdf   ,rwgt
!use yomlun   , only : nulout

implicit none

character(len=*), intent(in) :: cdirectory

real(kind=jprb) :: zwtsm(jpg)

integer(kind=jpim) :: i, ibnd, ig, igc, igcsm, ind, ipr, iprsm, ipt

real(kind=jprb) :: zfp, zrtfp, zwtsum
real(kind=jphook) :: zhook_handle

!#include "surrtmcf.intfb.h"

! # 1 "./include/surrtftr.intfb.h" 1
interface
subroutine surrtftr
end subroutine surrtftr
end interface
! # 34 "ifsrrtm/rrtm_init_140gp.f90" 2


! # 1 "./include/rrtm_kgb1.intfb.h" 1
interface
subroutine rrtm_kgb1(directory)
character(len=*), intent(in) :: directory
end subroutine rrtm_kgb1
end interface
! # 36 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb10.intfb.h" 1
interface
subroutine rrtm_kgb10
end subroutine rrtm_kgb10
end interface
! # 37 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb11.intfb.h" 1
interface
subroutine rrtm_kgb11
end subroutine rrtm_kgb11
end interface
! # 38 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb12.intfb.h" 1
interface
subroutine rrtm_kgb12
end subroutine rrtm_kgb12
end interface
! # 39 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb13.intfb.h" 1
interface
subroutine rrtm_kgb13
end subroutine rrtm_kgb13
end interface
! # 40 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb14.intfb.h" 1
interface
subroutine rrtm_kgb14
end subroutine rrtm_kgb14
end interface
! # 41 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb15.intfb.h" 1
interface
subroutine rrtm_kgb15
end subroutine rrtm_kgb15
end interface
! # 42 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb16.intfb.h" 1
interface
subroutine rrtm_kgb16
end subroutine rrtm_kgb16
end interface
! # 43 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb2.intfb.h" 1
interface
subroutine rrtm_kgb2
end subroutine rrtm_kgb2
end interface
! # 44 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb3.intfb.h" 1
interface
subroutine rrtm_kgb3
end subroutine rrtm_kgb3
end interface
! # 45 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb4.intfb.h" 1
interface
subroutine rrtm_kgb4
end subroutine rrtm_kgb4
end interface
! # 46 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb5.intfb.h" 1
interface
subroutine rrtm_kgb5
end subroutine rrtm_kgb5
end interface
! # 47 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb6.intfb.h" 1
interface
subroutine rrtm_kgb6
end subroutine rrtm_kgb6
end interface
! # 48 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb7.intfb.h" 1
interface
subroutine rrtm_kgb7
end subroutine rrtm_kgb7
end interface
! # 49 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb8.intfb.h" 1
interface
subroutine rrtm_kgb8
end subroutine rrtm_kgb8
end interface
! # 50 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_kgb9.intfb.h" 1
interface
subroutine rrtm_kgb9
end subroutine rrtm_kgb9
end interface
! # 51 "ifsrrtm/rrtm_init_140gp.f90" 2


! # 1 "./include/rrtm_cmbgb1.intfb.h" 1
interface
subroutine rrtm_cmbgb1
end subroutine rrtm_cmbgb1
end interface
! # 53 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb10.intfb.h" 1
interface
subroutine rrtm_cmbgb10
end subroutine rrtm_cmbgb10
end interface
! # 54 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb11.intfb.h" 1
interface
subroutine rrtm_cmbgb11
end subroutine rrtm_cmbgb11
end interface
! # 55 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb12.intfb.h" 1
interface
subroutine rrtm_cmbgb12
end subroutine rrtm_cmbgb12
end interface
! # 56 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb13.intfb.h" 1
interface
subroutine rrtm_cmbgb13
end subroutine rrtm_cmbgb13
end interface
! # 57 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb14.intfb.h" 1
interface
subroutine rrtm_cmbgb14
end subroutine rrtm_cmbgb14
end interface
! # 58 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb15.intfb.h" 1
interface
subroutine rrtm_cmbgb15
end subroutine rrtm_cmbgb15
end interface
! # 59 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb16.intfb.h" 1
interface
subroutine rrtm_cmbgb16
end subroutine rrtm_cmbgb16
end interface
! # 60 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb2.intfb.h" 1
interface
subroutine rrtm_cmbgb2
end subroutine rrtm_cmbgb2
end interface
! # 61 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb3.intfb.h" 1
interface
subroutine rrtm_cmbgb3
end subroutine rrtm_cmbgb3
end interface
! # 62 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb4.intfb.h" 1
interface
subroutine rrtm_cmbgb4
end subroutine rrtm_cmbgb4
end interface
! # 63 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb5.intfb.h" 1
interface
subroutine rrtm_cmbgb5
end subroutine rrtm_cmbgb5
end interface
! # 64 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb6.intfb.h" 1
interface
subroutine rrtm_cmbgb6
end subroutine rrtm_cmbgb6
end interface
! # 65 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb7.intfb.h" 1
interface
subroutine rrtm_cmbgb7
end subroutine rrtm_cmbgb7
end interface
! # 66 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb8.intfb.h" 1
interface
subroutine rrtm_cmbgb8
end subroutine rrtm_cmbgb8
end interface
! # 67 "ifsrrtm/rrtm_init_140gp.f90" 2

! # 1 "./include/rrtm_cmbgb9.intfb.h" 1
interface
subroutine rrtm_cmbgb9
end subroutine rrtm_cmbgb9
end interface
! # 68 "ifsrrtm/rrtm_init_140gp.f90" 2

if (lhook) call dr_hook('rrtm_init_140gp',0,zhook_handle)

!call surrtmcf
call surrtftr

! read the absorption-related coefficients over the 16 x 16 g-points

call rrtm_kgb1(cdirectory)
call rrtm_kgb2
call rrtm_kgb3
call rrtm_kgb4
call rrtm_kgb5
call rrtm_kgb6
call rrtm_kgb7
call rrtm_kgb8
call rrtm_kgb9
call rrtm_kgb10
call rrtm_kgb11
call rrtm_kgb12
call rrtm_kgb13
call rrtm_kgb14
call rrtm_kgb15
call rrtm_kgb16

!  calculate lookup tables for functions needed in routine taumol (taugb2)

corr1(0) = 1.0_jprb
corr1(200) = 1.0_jprb
corr2(0) = 1.0_jprb
corr2(200) = 1.0_jprb
do i = 1,199
  zfp = 0.005_jprb*real(i)
  zrtfp = sqrt(zfp)
  corr1(i) = zrtfp/zfp
  corr2(i) = (1.0_jprb-zrtfp)/(1.0_jprb-zfp)
enddo

!  perform g-point reduction from 16 per band (256 total points) to
!  a band dependant number (140 total points) for all absorption
!  coefficient input data and planck fraction input data.
!  compute relative weighting for new g-point combinations.

igcsm = 0
do ibnd = 1,jpband
  iprsm = 0
  if (ngc(ibnd) < 16) then
    do igc = 1,ngc(ibnd)
      igcsm = igcsm + 1
      zwtsum = 0.0_jprb
      do ipr = 1, ngn(igcsm)
        iprsm = iprsm + 1
        zwtsum = zwtsum + wt(iprsm)
      enddo
      zwtsm(igc) = zwtsum
    enddo

    do ig = 1,ng(ibnd)
      ind = (ibnd-1)*16 + ig
      rwgt(ind) = wt(ig)/zwtsm(ngm(ind))
    enddo
  else
    do ig = 1,ng(ibnd)
      igcsm = igcsm + 1
      ind = (ibnd-1)*16 + ig
      rwgt(ind) = 1.0_jprb
    enddo
  endif
enddo

!  initialize arrays for combined planck fraction data.

do ipt = 1,13
  do ipr = 1, jpgpt
    frefa(ipr,ipt) = 0.0_jprb
    frefadf(ipr,ipt) = 0.0_jprb
  enddo
enddo
do ipt = 1,6
  do ipr = 1, jpgpt
    frefb(ipr,ipt) = 0.0_jprb
    frefbdf(ipr,ipt) = 0.0_jprb
  enddo
enddo

!  reduce g-points for relevant data in each lw spectral band.

call rrtm_cmbgb1
call rrtm_cmbgb2
call rrtm_cmbgb3
call rrtm_cmbgb4
call rrtm_cmbgb5
call rrtm_cmbgb6
call rrtm_cmbgb7
call rrtm_cmbgb8
call rrtm_cmbgb9
call rrtm_cmbgb10
call rrtm_cmbgb11
call rrtm_cmbgb12
call rrtm_cmbgb13
call rrtm_cmbgb14
call rrtm_cmbgb15
call rrtm_cmbgb16

if (lhook) call dr_hook('rrtm_init_140gp',1,zhook_handle)
end subroutine rrtm_init_140gp
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

