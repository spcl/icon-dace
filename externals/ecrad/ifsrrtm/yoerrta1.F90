! # 1 "ifsrrtm/yoerrta1.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerrta1.f90"
! this file has been modified for the use in icon

module yoerrta1

use parkind1  ,only : jpim     ,jprb

implicit none

public

save

!     -----------------------------------------------------------------
!*    ** *yoerrta1* - rrtm coefficients for interval 1
!     band 1:  10-250 cm-1 (low - h2o; high - h2o)
!     abozzo may 2013 update to last version of rrtmg
!     band 1:  10-350 cm-1
!     -----------------------------------------------------------------

integer(kind=jpim), parameter :: ng1  = 10

real(kind=jprb) :: fracrefa(ng1)  , fracrefb(ng1)
real(kind=jprb) :: ka(5,13,ng1)   , absa(65,ng1)
real(kind=jprb) :: kb(5,13:59,ng1), absb(235,ng1)
real(kind=jprb) :: ka_mn2(19,ng1) , kb_mn2(19,ng1)
real(kind=jprb) :: selfref(10,ng1), forref(4,ng1)

equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

!$acc declare create(fracrefa, fracrefb, ka, absa, kb, absb, ka_mn2, kb_mn2, &
!$acc                selfref, forref)

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package ** rrtm lw radiation **

!     j.-j. morcrette       e.c.m.w.f.      98/07/14

!  name     type     purpose
!  ----   : ----   : ---------------------------------------------------
! absa    : real     absorption coefficient of major absorber for m reference tropospheric 
!                    pressures and n reference tropospheric temperatures
! absb    : real     absorption coefficient of secondary absorber for m reference stratospheric
!                    pressures and n reference stratospheric temperatures 
! fracrefa: real     distance from r and t reference tabulated points (troposphere)
! fracrefb: real     distance from r and t reference tabulated points (stratosphere)
! forref  : real     foreign broadening coefficient for water vapour
! ka      : real     absorption coefficient of major absorber    
! kb      : real     absorption coefficient of secondary absorber    
! selfref : real     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
end module yoerrta1
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

