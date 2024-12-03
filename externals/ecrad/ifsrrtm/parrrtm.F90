! # 1 "ifsrrtm/parrrtm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/parrrtm.f90"
module parrrtm

use parkind1  ,only : jpim

implicit none

public

save

!     ------------------------------------------------------------------
!     parameters relevant to aer's rrtm-lw radiation scheme

!     19980714  jjmorcrette
!     20110322  jjmorcrette : additional comments
!     20110603  jjmorcrette reduced number of g-points
!     ------------------------------------------------------------------

!-- basic spectral information unrelated to number of g-points
! jpg    : maximum possible number of g-points in each band of rrtm_lw
! jpband : number of longwave spectral bands
! jpxsec : number of cross-sections for active trace gases
! jpinpx : maximum dimension of the array of active trace gases 
! jpgpt  : total number of g-points in the (operational) spectrally-reduced rrtm_lw

integer(kind=jpim), parameter :: jpg    = 16
integer(kind=jpim), parameter :: jpband = 16
integer(kind=jpim), parameter :: jpxsec = 4
integer(kind=jpim), parameter :: jpinpx = 35
integer(kind=jpim), parameter :: jpgmax = 256

!-- configuration for eps with 70 g-points

!integer(kind=jpim), parameter :: jpgpt  = 70

!integer(kind=jpim), parameter :: ng1  = 4
!integer(kind=jpim), parameter :: ng2  = 7
!integer(kind=jpim), parameter :: ng3  = 8
!integer(kind=jpim), parameter :: ng4  = 7
!integer(kind=jpim), parameter :: ng5  = 8
!integer(kind=jpim), parameter :: ng6  = 4
!integer(kind=jpim), parameter :: ng7  = 6
!integer(kind=jpim), parameter :: ng8  = 4
!integer(kind=jpim), parameter :: ng9  = 6
!integer(kind=jpim), parameter :: ng10 = 3
!integer(kind=jpim), parameter :: ng11 = 4
!integer(kind=jpim), parameter :: ng12 = 4
!integer(kind=jpim), parameter :: ng13 = 2
!integer(kind=jpim), parameter :: ng14 = 1
!integer(kind=jpim), parameter :: ng15 = 1
!integer(kind=jpim), parameter :: ng16 = 1

!integer(kind=jpim), parameter :: ngs1  = 4
!integer(kind=jpim), parameter :: ngs2  = 11
!integer(kind=jpim), parameter :: ngs3  = 19
!integer(kind=jpim), parameter :: ngs4  = 26
!integer(kind=jpim), parameter :: ngs5  = 34
!integer(kind=jpim), parameter :: ngs6  = 38
!integer(kind=jpim), parameter :: ngs7  = 44
!integer(kind=jpim), parameter :: ngs8  = 48
!integer(kind=jpim), parameter :: ngs9  = 54
!integer(kind=jpim), parameter :: ngs10 = 57
!integer(kind=jpim), parameter :: ngs11 = 61
!integer(kind=jpim), parameter :: ngs12 = 65
!integer(kind=jpim), parameter :: ngs13 = 67
!integer(kind=jpim), parameter :: ngs14 = 68
!integer(kind=jpim), parameter :: ngs15 = 69


!-- configuration with 140 g-points

!integer(kind=jpim), parameter :: jpgptf = 140
!integer(kind=jpim), parameter :: jpgptr = 140
!integer(kind=jpim), parameter :: jpgpt  = 140

!integer(kind=jpim), parameter :: ng1  = 8
!integer(kind=jpim), parameter :: ng2  = 14
!integer(kind=jpim), parameter :: ng3  = 16
!integer(kind=jpim), parameter :: ng4  = 14
!integer(kind=jpim), parameter :: ng5  = 16
!integer(kind=jpim), parameter :: ng6  = 8
!integer(kind=jpim), parameter :: ng7  = 12
!integer(kind=jpim), parameter :: ng8  = 8
!integer(kind=jpim), parameter :: ng9  = 12
!integer(kind=jpim), parameter :: ng10 = 6
!integer(kind=jpim), parameter :: ng11 = 8
!integer(kind=jpim), parameter :: ng12 = 8
!integer(kind=jpim), parameter :: ng13 = 4
!integer(kind=jpim), parameter :: ng14 = 2
!integer(kind=jpim), parameter :: ng15 = 2
!integer(kind=jpim), parameter :: ng16 = 2

!integer(kind=jpim), parameter :: ngs1  = 8
!integer(kind=jpim), parameter :: ngs2  = 22
!integer(kind=jpim), parameter :: ngs3  = 38
!integer(kind=jpim), parameter :: ngs4  = 52
!integer(kind=jpim), parameter :: ngs5  = 68
!integer(kind=jpim), parameter :: ngs6  = 76
!integer(kind=jpim), parameter :: ngs7  = 88
!integer(kind=jpim), parameter :: ngs8  = 96
!integer(kind=jpim), parameter :: ngs9  = 108
!integer(kind=jpim), parameter :: ngs10 = 114
!integer(kind=jpim), parameter :: ngs11 = 122
!integer(kind=jpim), parameter :: ngs12 = 130
!integer(kind=jpim), parameter :: ngs13 = 134
!integer(kind=jpim), parameter :: ngs14 = 136
!integer(kind=jpim), parameter :: ngs15 = 138

!     ------------------------------------------------------------------
end module parrrtm
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

