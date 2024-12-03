! # 1 "ifsrrtm/yoerrtrf.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerrtrf.f90"
! this file has been modified for the use in icon

module yoerrtrf

use parkind1  ,only : jprb

implicit none

public

save

!     -----------------------------------------------------------------
!*    ** *yoerrtrf* - rrtm reference atmosphere
!     -----------------------------------------------------------------

real(kind=jprb) , dimension(59) :: pref
real(kind=jprb) , dimension(59) :: preflog
real(kind=jprb) , dimension(59) :: tref
real(kind=jprb)  :: chi_mls(7,59)

!$acc declare create(pref, preflog, tref, chi_mls)

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package ** rrtm lw radiation **

!     j.-j. morcrette       e.c.m.w.f.      98/01/15

!  name     type     purpose
!  ----  :  ----   : ---------------------------------------------------
! pref   :  real    
! preflog: real
! tref   : real
!     -----------------------------------------------------------------
end module yoerrtrf
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

