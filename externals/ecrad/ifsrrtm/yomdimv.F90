! # 1 "ifsrrtm/yomdimv.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yomdimv.f90"
module yomdimv

use parkind1  ,only : jpim

implicit none

save

!     ------------------------------------------------------------------

type :: tdimv

!*    dimensions of model working arrays

! === vertical resolution =====================================================

! nflevg : number of levels in grid point space
! nflevl : number of levels in fourier and legendre space
! nflevlmx : maximum nflevl among all pes
! nflsur : over dimensioning of nflevl for technical reasons, always odd
! nflsul : number of additional levels for semi-lagrangian
! nflsa  = 1    -nflsul
! nflen  = nflevg+nflsul
! niolevg : number of levels in the whole atmosphere (used for i/os and definitions) ; 
!           nflevg can be a truncation of niolevg

integer(kind=jpim) :: nflevg
integer(kind=jpim) :: nflevl
integer(kind=jpim) :: nflevlmx
integer(kind=jpim) :: nflsur
integer(kind=jpim) :: nflsul
integer(kind=jpim) :: nflsa
integer(kind=jpim) :: nflen
integer(kind=jpim) :: nflevsf
integer(kind=jpim) :: niolevg

end type tdimv

type(tdimv), pointer :: yrdimv => null()

!     ------------------------------------------------------------------

end module yomdimv
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

