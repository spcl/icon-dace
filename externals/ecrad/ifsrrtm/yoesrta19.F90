! # 1 "ifsrrtm/yoesrta19.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoesrta19.f90"
! this file has been modified for the use in icon

module yoesrta19

use parkind1  ,only : jpim     ,jprb,jprd

implicit none

public

save

!     -----------------------------------------------------------------
!*    ** *yoesrta19* - srtm coefficients for interval 19
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)
!     -----------------------------------------------------------------

integer(kind=jpim), parameter :: jpg = 16, ng19 = 16

real(kind=jprb) :: ka(9,5,13,jpg) 
real(kind=jprb) :: kb(5,13:59,jpg)
real(kind=jprd) :: ka_d(9,5,13,jpg) 
real(kind=jprd) :: kb_d(5,13:59,jpg)
real(kind=jprb) :: selfref(10,jpg),forref(3,jpg)
real(kind=jprb) :: sfluxref(jpg,9)
real(kind=jprb) :: rayl            ,strrat
integer(kind=jpim) :: layreffr

real(kind=jprb) :: kac(9,5,13,ng19) ,absa(585,ng19)
real(kind=jprb) :: kbc(5,13:59,ng19),absb(235,ng19)
real(kind=jprb) :: selfrefc(10,ng19),forrefc(3,ng19)
real(kind=jprb) :: sfluxrefc(ng19,9)

!equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!$acc declare create(kac, absa, kbc, absb, selfrefc, forrefc, sfluxrefc)

!     -----------------------------------------------------------------
!        * e.c.m.w.f. physics package ** rrtm sw radiation **

!     j.-j. morcrette       e.c.m.w.f.      02/10/29
!     m. j. iacono          aer             12/09/03

!  name     type     purpose
!  ----   : ----   : ---------------------------------------------------
! ka      : real     absorption coefficient of major absorber
! kb      : real     absorption coefficient of secondary absorber
! selfref : real     self brodening coefficient for water vapour
! forref  : real     foreign broadening coefficient for water vapour
! sfluxref: real     incident solar radiation in the spectral interval
! rayl    : real     rayleigh scattering parameter
! strrat  : real     weighting factor for the transition between tropospheric 
!                    and stratospheric computations
! layreffr: integer  reference level for the transition
! kac     : real     reduced g-point array for ka
! kbc     : real     reduced g-point array for kb
! selfrefc: real     reduced g-point array for selfref
! forrefc : real     reduced g-point array for forref
!sfluxrefc: real     reduced g-point array for sfluxref
!     -----------------------------------------------------------------
end module yoesrta19

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

