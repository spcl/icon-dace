! # 1 "ifsrrtm/yoerrta8.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerrta8.f90"
! this file has been modified for the use in icon

module yoerrta8

use parkind1  ,only : jpim     ,jprb

implicit none

public

save

!     -----------------------------------------------------------------
!*    ** *yoerrta8* - rrtm coefficients for interval 8
!     band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!     abozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

integer(kind=jpim), parameter :: ng8  = 8

real(kind=jprb) , dimension(ng8) :: fracrefa
real(kind=jprb) , dimension(ng8) :: fracrefb
real(kind=jprb) , dimension(ng8) :: cfc12
real(kind=jprb) , dimension(ng8) :: cfc22adj

real(kind=jprb) :: ka(5,13,ng8)    ,absa(65,ng8)
real(kind=jprb) :: kb(5,13:59,ng8) ,absb(235,ng8)
real(kind=jprb) :: ka_mco2(19,ng8)
real(kind=jprb) :: ka_mn2o(19,ng8)
real(kind=jprb) :: ka_mo3(19,ng8)
real(kind=jprb) :: kb_mco2(19,ng8)
real(kind=jprb) :: kb_mn2o(19,ng8)
real(kind=jprb) :: selfref(10,ng8)
real(kind=jprb) :: forref(4,ng8)


equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

!$acc declare create(fracrefa, fracrefb, cfc12, cfc22adj, ka, absa, kb, absb, &
!$acc                ka_mco2, ka_mn2o, ka_mo3, kb_mco2, kb_mn2o, selfref, forref)

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package *

!     j.-j. morcrette       e.c.m.w.f.      98/07/14

!  name     type     purpose
!  ----   : ----   : ---------------------------------------------------
! absa    : real     absorption coefficient of major absorber for m reference tropospheric 
!                    pressures and n reference tropospheric temperatures 
! absb    : real     absorption coefficient of secondary absorber for m reference stratospheric
!                    pressures and n reference stratospheric temperatures 
! absco2a : real     absorption coefficient for co2 for m reference tropospheric 
!                    pressures and n reference tropospheric temperatures 
! absco2b : real     absorption coefficient for co2 for m reference stratospheric 
!                    pressures and n reference stratospheric temperatures     
! absn2oa : real     absorption coefficient for n2o for m reference tropospheric 
!                    pressures and n reference tropospheric temperatures     
! absn2ob : real     absorption coefficient for n2o for m reference stratospheric 
!                    pressures and n reference stratospheric temperatures 
! cfc12   : real     absorption coefficient for cfc-12
! cfc22adj: real     absorption coefficient for cfc-22 (adjusted)
! fracrefa: real     distance from r and t reference tabulated points (troposphere)
! fracrefb: real     distance from r and t reference tabulated points (stratosphere)
! h2oref  : real     reference profile for h2o
! ka      : real     absorption coefficient of major absorber (equiv. to absa)   
! kb      : real     absorption coefficient of secondary absorber (equiv. to absb)   
! n2oref  : real     reference profile for n2o
! o3ref   : real     reference profile for o3
! selfref : real     self broadening coefficient for water vapour
!     -----------------------------------------------------------------
end module yoerrta8
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

