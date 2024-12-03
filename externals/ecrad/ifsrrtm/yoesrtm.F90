! # 1 "ifsrrtm/yoesrtm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoesrtm.f90"
module yoesrtm

use parkind1  ,only : jpim
use parsrtm   ,only : jpgmax

implicit none

public

save

!     ------------------------------------------------------------------
!     parameters relevant to aer's rrtm-sw radiation scheme: part 2

!     20110610 jjmorcrette

!     modified to allow possibilities of different g-point numbers.  
!     ------------------------------------------------------------------

!integer(kind=jpim) :: jpgpt
!integer(kind=jpim) :: jpgsw

!integer(kind=jpim) :: ng16
!integer(kind=jpim) :: ng17
!integer(kind=jpim) :: ng18
!integer(kind=jpim) :: ng19
!integer(kind=jpim) :: ng20
!integer(kind=jpim) :: ng21
!integer(kind=jpim) :: ng22
!integer(kind=jpim) :: ng23
!integer(kind=jpim) :: ng24
!integer(kind=jpim) :: ng25
!integer(kind=jpim) :: ng26
!integer(kind=jpim) :: ng27
!integer(kind=jpim) :: ng28
!integer(kind=jpim) :: ng29

!-- ngnn : number of g-points in bands nn=16 to 29
!- as used for the (operational) 112 g-points version of rrtm_sw
integer(kind=jpim), parameter :: jpgpt  = 112
integer(kind=jpim), parameter :: jpgsw  = 112

integer(kind=jpim), parameter :: ng16 = 6
integer(kind=jpim), parameter :: ng17 = 12
integer(kind=jpim), parameter :: ng18 = 8
integer(kind=jpim), parameter :: ng19 = 8
integer(kind=jpim), parameter :: ng20 = 10
integer(kind=jpim), parameter :: ng21 = 10
integer(kind=jpim), parameter :: ng22 = 2
integer(kind=jpim), parameter :: ng23 = 10
integer(kind=jpim), parameter :: ng24 = 8
integer(kind=jpim), parameter :: ng25 = 6
integer(kind=jpim), parameter :: ng26 = 6
integer(kind=jpim), parameter :: ng27 = 8
integer(kind=jpim), parameter :: ng28 = 6
integer(kind=jpim), parameter :: ng29 = 12


integer(kind=jpim) :: ngn(jpgmax), ngbsw(jpgmax)

!     ------------------------------------------------------------------
end module yoesrtm

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

