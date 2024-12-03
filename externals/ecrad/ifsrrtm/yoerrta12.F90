! # 1 "ifsrrtm/yoerrta12.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerrta12.f90"
! this file has been modified for the use in icon

module yoerrta12

use parkind1  ,only : jpim     ,jprb

implicit none

public

save

!     -----------------------------------------------------------------
!*    ** *yoerrta12* - rrtm coefficients for interval 12
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!     abozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

integer(kind=jpim), parameter :: ng12 = 8

real(kind=jprb) :: fracrefa(ng12,9)
real(kind=jprb) :: ka(9,5,13,ng12) ,absa(585,ng12)
real(kind=jprb) :: selfref(10,ng12)
real(kind=jprb) :: forref(4,ng12)

real(kind=jprb) :: strrat

equivalence (ka(1,1,1,1),absa(1,1))

!$acc declare create(fracrefa, ka, absa, selfref, forref)

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package ** rrtm lw radiation **

!     j.-j. morcrette       e.c.m.w.f.      98/07/14

!  name     type     purpose
!  ----   : ----   : ---------------------------------------------------
! absa    : real     absorption coefficient of major absorber for m reference tropospheric 
!                    pressures and n reference tropospheric temperatures 
! fracrefa: real     distance from r and t reference tabulated points (troposphere)
! ka      : real     absorption coefficient of major absorber (equiv. to absa)   
! selfref : real     self broadening coefficient for water vapour
! strrat  : real     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
end module yoerrta12
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

