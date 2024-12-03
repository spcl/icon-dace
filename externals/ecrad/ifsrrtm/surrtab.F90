! # 1 "ifsrrtm/surrtab.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/surrtab.f90"
subroutine surrtab

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package ** aer's rrtm lw radiation **

!     j.-j. morcrette       e.c.m.w.f.      98/07/14

!     -----------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook, dr_hook, jphook

use yoerrtab , only : trans, bpade

implicit none

integer(kind=jpim) :: itr

real(kind=jprb) :: ztau, ztfn
real(kind=jphook) :: zhook_handle

if (lhook) call dr_hook('surrtab',0,zhook_handle)
bpade=1.0_jprb/0.278_jprb
trans(0)   =1.0_jprb
trans(5000)=0.0_jprb
do itr=1,4999
  ztfn=real(itr)/5000._jprb
  ztau=bpade*ztfn/(1.0_jprb-ztfn)
  trans(itr)=exp(-ztau)
enddo

!     -----------------------------------------------------------------

if (lhook) call dr_hook('surrtab',1,zhook_handle)
end subroutine surrtab
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

