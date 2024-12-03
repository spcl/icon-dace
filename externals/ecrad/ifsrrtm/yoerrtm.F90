! # 1 "ifsrrtm/yoerrtm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerrtm.f90"
module yoerrtm

use parkind1  ,only : jpim
use parrrtm   ,only : jpgmax

implicit none

public

save

!     ------------------------------------------------------------------
!     parameters relevant to aer's rrtm-lw radiation scheme: part 2

!     20110613 jjmorcrette

!     modified to allow possibilities of different g-point numbers.  
!     ------------------------------------------------------------------

!integer(kind=jpim) :: jpgpt
!integer(kind=jpim) :: jpglw

!integer(kind=jpim) :: ng1
!integer(kind=jpim) :: ng2
!integer(kind=jpim) :: ng3
!integer(kind=jpim) :: ng4
!integer(kind=jpim) :: ng5
!integer(kind=jpim) :: ng6
!integer(kind=jpim) :: ng7
!integer(kind=jpim) :: ng8
!integer(kind=jpim) :: ng9
!integer(kind=jpim) :: ng10
!integer(kind=jpim) :: ng11
!integer(kind=jpim) :: ng12
!integer(kind=jpim) :: ng13
!integer(kind=jpim) :: ng14
!integer(kind=jpim) :: ng15
!integer(kind=jpim) :: ng16

!integer(kind=jpim) :: ngs1
!integer(kind=jpim) :: ngs2
!integer(kind=jpim) :: ngs3
!integer(kind=jpim) :: ngs4
!integer(kind=jpim) :: ngs5
!integer(kind=jpim) :: ngs6
!integer(kind=jpim) :: ngs7
!integer(kind=jpim) :: ngs8
!integer(kind=jpim) :: ngs9
!integer(kind=jpim) :: ngs10
!integer(kind=jpim) :: ngs11
!integer(kind=jpim) :: ngs12
!integer(kind=jpim) :: ngs13
!integer(kind=jpim) :: ngs14
!integer(kind=jpim) :: ngs15
!integer(kind=jpim) :: ngs16

integer(kind=jpim), parameter :: jpgpt  = 140
integer(kind=jpim), parameter :: jpglw  = 140

!-- ngnn : number of g-points in each longwave spectral band
integer(kind=jpim), parameter :: ng1  = 10
integer(kind=jpim), parameter :: ng2  = 12
integer(kind=jpim), parameter :: ng3  = 16
integer(kind=jpim), parameter :: ng4  = 14
integer(kind=jpim), parameter :: ng5  = 16
integer(kind=jpim), parameter :: ng6  = 8
integer(kind=jpim), parameter :: ng7  = 12
integer(kind=jpim), parameter :: ng8  = 8
integer(kind=jpim), parameter :: ng9  = 12
integer(kind=jpim), parameter :: ng10 = 6
integer(kind=jpim), parameter :: ng11 = 8
integer(kind=jpim), parameter :: ng12 = 8
integer(kind=jpim), parameter :: ng13 = 4
integer(kind=jpim), parameter :: ng14 = 2
integer(kind=jpim), parameter :: ng15 = 2
integer(kind=jpim), parameter :: ng16 = 2
!-- ngsnn: accumulated number of g-points at the beginning of spectral band nn+1
integer(kind=jpim), parameter :: ngs1  = 10
integer(kind=jpim), parameter :: ngs2  = 22
integer(kind=jpim), parameter :: ngs3  = 38
integer(kind=jpim), parameter :: ngs4  = 52
integer(kind=jpim), parameter :: ngs5  = 68
integer(kind=jpim), parameter :: ngs6  = 76
integer(kind=jpim), parameter :: ngs7  = 88
integer(kind=jpim), parameter :: ngs8  = 96
integer(kind=jpim), parameter :: ngs9  = 108
integer(kind=jpim), parameter :: ngs10 = 114
integer(kind=jpim), parameter :: ngs11 = 122
integer(kind=jpim), parameter :: ngs12 = 130
integer(kind=jpim), parameter :: ngs13 = 134
integer(kind=jpim), parameter :: ngs14 = 136
integer(kind=jpim), parameter :: ngs15 = 138


integer(kind=jpim) :: ngn(jpgmax), ngblw(jpgmax)

!     ------------------------------------------------------------------
end module yoerrtm

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

