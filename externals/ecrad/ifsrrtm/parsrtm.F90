! # 1 "ifsrrtm/parsrtm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/parsrtm.f90"
module parsrtm

use parkind1  ,only : jpim

implicit none

public

save

!     ------------------------------------------------------------------
!     parameters relevant to aer's rrtm-sw radiation scheme

!     030224  jjmorcrette

!     modified for g-point reduction from 224 to 112.  
!     swap code below to restore 224 g-point set. 
!     mar2004 mjiacono, aer
!     20110322 jjmorcrette : additional comments
!     20110603 jjmorcrette reduced number of g-points
!     ------------------------------------------------------------------

!-- basic spectral information unrelated to number of g-points
! jpg     : integer : maximum number of g-points in a given spectral band
! jpband  : integer : total number of spectral bands 
! jpsw    : integer : total number of shortwave spectral bands
! jpb1    : integer : starting index of shortwave spectrum
! jpb2    : integer : end index of shortwave spectrum

integer(kind=jpim), parameter :: jpg    = 16
integer(kind=jpim), parameter :: jpband = 29
integer(kind=jpim), parameter :: jpsw   = 14
integer(kind=jpim), parameter :: jpb1   = 16
integer(kind=jpim), parameter :: jpb2   = 29
integer(kind=jpim), parameter :: jpgmax = 224

!-- other information that could be relevant for rrtm_sw
!-- nb: the following parameters are unused within the ecmwf ifs. 
!       they relate to the description of the optical properties 
!       in the original cloud model embedded in rrtm_sw
!integer(kind=jpim), parameter :: jmcmu  = 32
!integer(kind=jpim), parameter :: jmumu  = 32
!integer(kind=jpim), parameter :: jmphi  = 3
!integer(kind=jpim), parameter :: jmxang = 4
!integer(kind=jpim), parameter :: jmxstr = 16

!-- original spectral grid before spectral averaging
!-- original from aer, inc with 224 g-points
integer(kind=jpim), parameter :: ngs16 = 0
integer(kind=jpim), parameter :: ngs17 = 16
integer(kind=jpim), parameter :: ngs18 = 32
integer(kind=jpim), parameter :: ngs19 = 48
integer(kind=jpim), parameter :: ngs20 = 64
integer(kind=jpim), parameter :: ngs21 = 80
integer(kind=jpim), parameter :: ngs22 = 96
integer(kind=jpim), parameter :: ngs23 = 112
integer(kind=jpim), parameter :: ngs24 = 128
integer(kind=jpim), parameter :: ngs25 = 144
integer(kind=jpim), parameter :: ngs26 = 160
integer(kind=jpim), parameter :: ngs27 = 176
integer(kind=jpim), parameter :: ngs28 = 192
integer(kind=jpim), parameter :: ngs29 = 208

!-------------------------------------------------------------------------------
!-- ngnn : number of g-points in bands nn=16 to 29
!- as used in the ng g-points version of rrtm_sw
!-------------------------------------------------------------------------------
!-- configuration with 14 spectral intervals
!   and a total of 56 g-points (14xvariable number)

!integer(kind=jpim), parameter :: jpgpt  = 56
!
!integer(kind=jpim), parameter :: ng16 = 3
!integer(kind=jpim), parameter :: ng17 = 6
!integer(kind=jpim), parameter :: ng18 = 4
!integer(kind=jpim), parameter :: ng19 = 4
!integer(kind=jpim), parameter :: ng20 = 5
!integer(kind=jpim), parameter :: ng21 = 5
!integer(kind=jpim), parameter :: ng22 = 1
!integer(kind=jpim), parameter :: ng23 = 5
!integer(kind=jpim), parameter :: ng24 = 4
!integer(kind=jpim), parameter :: ng25 = 3
!integer(kind=jpim), parameter :: ng26 = 3
!integer(kind=jpim), parameter :: ng27 = 4
!integer(kind=jpim), parameter :: ng28 = 3
!integer(kind=jpim), parameter :: ng29 = 6
!-------------------------------------------------------------------------------
!-- configuration with 14 spectral intervals
!   and a total of 112 g-points (14xvariable number)
!
!integer(kind=jpim), parameter :: jpgpt  = 112
!
!integer(kind=jpim), parameter :: ng16 = 6
!integer(kind=jpim), parameter :: ng17 = 12
!integer(kind=jpim), parameter :: ng18 = 8
!integer(kind=jpim), parameter :: ng19 = 8
!integer(kind=jpim), parameter :: ng20 = 10
!integer(kind=jpim), parameter :: ng21 = 10
!integer(kind=jpim), parameter :: ng22 = 2
!integer(kind=jpim), parameter :: ng23 = 10
!integer(kind=jpim), parameter :: ng24 = 8
!integer(kind=jpim), parameter :: ng25 = 6
!integer(kind=jpim), parameter :: ng26 = 6
!integer(kind=jpim), parameter :: ng27 = 8
!integer(kind=jpim), parameter :: ng28 = 6
!integer(kind=jpim), parameter :: ng29 = 12

!-------------------------------------------------------------------------------
!-- configuration with 14 spectral intervals 
!   and a total of 224 g-points (14x16)
! 
!integer(kind=jpim), parameter :: jpgpt  = 224

!integer(kind=jpim), parameter :: ng16 = 16
!integer(kind=jpim), parameter :: ng17 = 16
!integer(kind=jpim), parameter :: ng18 = 16
!integer(kind=jpim), parameter :: ng19 = 16
!integer(kind=jpim), parameter :: ng20 = 16
!integer(kind=jpim), parameter :: ng21 = 16
!integer(kind=jpim), parameter :: ng22 = 16
!integer(kind=jpim), parameter :: ng23 = 16
!integer(kind=jpim), parameter :: ng24 = 16
!integer(kind=jpim), parameter :: ng25 = 16
!integer(kind=jpim), parameter :: ng26 = 16
!integer(kind=jpim), parameter :: ng27 = 16
!integer(kind=jpim), parameter :: ng28 = 16
!integer(kind=jpim), parameter :: ng29 = 16

!     ------------------------------------------------------------------
end module parsrtm

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

