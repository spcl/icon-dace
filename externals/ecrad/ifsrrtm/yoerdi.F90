! # 1 "ifsrrtm/yoerdi.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoerdi.f90"
module yoerdi

use parkind1  ,only : jprb

implicit none

public

save

!     -----------------------------------------------------------------
!*    ** *yoerdi* - coefficients within radiation interface
!     -----------------------------------------------------------------

real(kind=jprb) :: rrae
real(kind=jprb) :: rsundur
real(kind=jprb) :: rcardi
real(kind=jprb) :: rch4
real(kind=jprb) :: rn2o
real(kind=jprb) :: rno2
real(kind=jprb) :: ro3
real(kind=jprb) :: rccl4
real(kind=jprb) :: rcfc11
real(kind=jprb) :: rcfc12
real(kind=jprb) :: rcfc22
real(kind=jprb) :: repclc
real(kind=jprb) :: reph2o
real(kind=jprb) :: rcco2, rcch4, rcn2o, rcno2, rccfc11, rccfc12, rccfc22, rcccl4
real(kind=jprb) :: rsolinc

!        * e.c.m.w.f. physics package *

!     original  j.-j. morcrette       e.c.m.w.f.      89/07/14
!     modified  p. viterbo    99/03/26    surface tiling
!     modified  p. viterbo    24/05/2004  surf library
!     modified jjmorcrette    2005/01/19  ghg and solar constant variability

!  name     type     purpose
!  ----  :  ----   : ---------------------------------------------------
! rrae   : effect of earth's curvature on cosine solar zenith angle
! rsundur: minimum direct solar for computing solar duration
! rcardi : specific atmospheric content in co2
! rch4, rn2o, rno2, ro3, rcfc11, rcfc12 mass mixing ratio of various trace gases
! rcch4, rcn2o, ... mass mixing ratio of various trace gases in climate mode
! repclc : security to avoid zero or one cloud covers
! reph2o : security to avoid water vapour content in a layer
!          to be more than the respective value at saturation.
!     -----------------------------------------------------------------
end module yoerdi
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

