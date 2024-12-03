! # 1 "ifsrrtm/yoesw.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/yoesw.f90"
module yoesw

use parkind1  ,only : jpim     ,jprb

implicit none

save

!     ------------------------------------------------------------------
!*    ** *yoesw* - coefficients for shortwave radiation transfer
!     ------------------------------------------------------------------

real(kind=jprb) :: apad(6,3,7)
real(kind=jprb) :: bpad(6,3,7)
real(kind=jprb) :: rray(6,6)
real(kind=jprb), allocatable :: rsun(:)
real(kind=jprb) :: rpdh1
real(kind=jprb) :: rpdu1
real(kind=jprb) :: rpnh
real(kind=jprb) :: rpnu
real(kind=jprb) :: rswce(6)
real(kind=jprb) :: rswcp(6)
real(kind=jprb) :: rtdh2o
real(kind=jprb) :: rtdumg
real(kind=jprb) :: rth2o
real(kind=jprb) :: rtumg
real(kind=jprb) :: d(6,3)
real(kind=jprb) :: rexpo3(6,2,7)
integer(kind=jpim) :: nexpo3(6)

real(kind=jprb) :: ryfwca(6)
real(kind=jprb) :: ryfwcb(6)
real(kind=jprb) :: ryfwcc(6)
real(kind=jprb) :: ryfwcd(6)
real(kind=jprb) :: ryfwce(6)
real(kind=jprb) :: ryfwcf(6)

real(kind=jprb) :: rebcua(6)
real(kind=jprb) :: rebcub(6)
real(kind=jprb) :: rebcuc(6)
real(kind=jprb) :: rebcud(6)
real(kind=jprb) :: rebcue(6)
real(kind=jprb) :: rebcuf(6)
real(kind=jprb) :: rebcug(16)
real(kind=jprb) :: rebcuh(16)
real(kind=jprb) :: rebcui(6)
real(kind=jprb) :: rebcuj(6)

real(kind=jprb) :: raswca(6)
real(kind=jprb) :: raswcb(6)
real(kind=jprb) :: raswcc(6)
real(kind=jprb) :: raswcd(6)
real(kind=jprb) :: raswce(6)
real(kind=jprb) :: raswcf(6)

real(kind=jprb) :: rfueta(16,3),rfuetb(16,4), rfuetc(16,4)
real(kind=jprb) :: rfulio(16,3)
real(kind=jprb) :: rhsavi(16,3)
real(kind=jprb) :: rlilia(16,5),rlilib(16,4)

real(kind=jprb) :: rflaa0(6)
real(kind=jprb) :: rflaa1(6)
real(kind=jprb) :: rflbb0(6)
real(kind=jprb) :: rflbb1(6)
real(kind=jprb) :: rflbb2(6)
real(kind=jprb) :: rflbb3(6)
real(kind=jprb) :: rflcc0(6)
real(kind=jprb) :: rflcc1(6)
real(kind=jprb) :: rflcc2(6)
real(kind=jprb) :: rflcc3(6)

real(kind=jprb) :: rfuaa0(6)
real(kind=jprb) :: rfuaa1(6)
real(kind=jprb) :: rfubb0(6)
real(kind=jprb) :: rfubb1(6)
real(kind=jprb) :: rfubb2(6)
real(kind=jprb) :: rfubb3(6)
real(kind=jprb) :: rfucc0(6)
real(kind=jprb) :: rfucc1(6)
real(kind=jprb) :: rfucc2(6)
real(kind=jprb) :: rfucc3(6)
real(kind=jprb) :: rfldd0(6)
real(kind=jprb) :: rfldd1(6)
real(kind=jprb) :: rfldd2(6)
real(kind=jprb) :: rfldd3(6)

real(kind=jprb) :: reffia

real(kind=jprb) :: rtaua(6,6)
real(kind=jprb) :: rpiza(6,6)
real(kind=jprb) :: rcga(6,6)
real(kind=jprb) :: raer(6,6)

integer(kind=jpim) :: nmpsrtm(14), ntyps

real(kind=jprb) :: radjust

!        * e.c.m.w.f. physics package *

!     j.-j. morcrette       e.c.m.w.f.      89/07/14

!  name     type     purpose
!  ----  :  ----   : ---------------------------------------------------
!  apad  :  real     pade approximants numerator
!  bpad  :  real     pade approximants denominator
!  d     :  real     transmission limit for infinite absorber amount
!  rray  :  real     rayleigh scattering coefficients
!  rsun  :  real     solar fraction in spectral intervals
!  rpdh1 :  1 + exponent pressure dependence h2o
!  rpdu1 :  1 + exponent pressure dependence uniformly mixed gases
!  rpnh  :  reference pressure factor for h2o
!  rpnu  :  reference pressure factor for uniformly mixed gases
!  rswce :  e-type, h2o continuum absorption coefficient 
!  rswcp :  p-type, h2o continuum absorption coefficient 
!  rtdh2o:  exponent temperature dependence h2o
!  rtdumg:  exponent temperature dependence uniformly mixed gases
!  rth2o :  reference temperature h2o
!  rtumg :  reference temperature uniformly mixed gases
!     -----------------------------------------------------------------

!        * e.c.m.w.f. physics package *

!     j.-j. morcrette       e.c.m.w.f.      89/07/14

!  name     type     purpose
!  ----  :  ----   : ---------------------------------------------------
!*    fouquart (1987) water cloud optical properties

! ryfwca :  real   : c1 in optical thickness formula
! ryfwcb :  real   : c2 in optical thickness formula
! ryfwcc :  real   : single scattering albedo parameter
! ryfwcd :  real   : single scattering albedo parameter
! ryfwce :  real   : single scattering albedo parameter
! ryfwcf :  real   : assymetry factor

!*    slingo (1989) water cloud optical properties

! raswca :  real   : c1 in optical thickness formula
! raswcb :  real   : c2 in optical thickness formula
! raswcc :  real   : single scattering albedo parameter
! raswcd :  real   : single scattering albedo parameter
! raswce :  real   : single scattering albedo parameter
! raswcf :  real   : assymetry factor

!*   lindner,li (2000) water cloud optical properties (rrtm)

! rlilia : real    : mass absorption coefficients (polynomial developm)
! rlilib : real    : 1-ssa coefficients  (polynomial developm)

!*    ice cloud optical properties derived from ebert-curry (1992)

! rebcua :  real   : c1 in optical thickness formula
! rebcub :  real   : c2 in optical thickness formula
! rebcuc :  real   : 1-c3  in single scattering albedo formula
! rebcud :  real   : c4 in single scattering albedo formula
! rebcue :  real   : c5 in assymetry factor formula
! rebcuf :  real   : c6 in assymetry factor formula
! rebcug :  real   : c7 in mass absorption coefficient formula
! rebcuh :  real   : c8 in mass absorption coefficient formula
! rebcui :  real   : c7 in mass absorption coefficient spectral formula
! rebcuj :  real   : c8 in mass absorption coefficient spectral formula

!*    ice cloud optical properties derived from sun-shine (1995)

! rshsue :  real   : e in single scattering albedo formula
! rshsuf :  real   : f in single scattering albedo formula
! rshsuh :  real   : h in assymetry factor formula
! rshsuk :  real   : k in assymetry factor formula
! rshsua :  real   : alpha in ssa correction factor formula
! rshsug :  real   : gamma in assymetry correction factor formula
! rshsufa:  real   : coefficients in temperature correction factor

! reffia :  real   : c9  in effective radius formula

!*    ice cloud optical properties derived from fu-liou (1993)

! rfulio :  real   : coefficients in expression for lw extinction coeff.
! rflaa  :  real   : coefficients in expression for sw extinction coeff.
! rflbb  :  real   : coefficients in expression for sw single scatt.alb.
! rflcc  :  real   : coefficients in expression for sw assymetry factor
! rfldd  :  real   : coefficients in expression for sw assymetry factor

!*    ice cloud optical properties derived from fu (1996) & fu et al. (1998)

! rfueta :  real   : coefficients in expression for lw extinction coeff.
! rfuaa  :  real   : coefficients in expression for sw extinction coeff.
! rfubb  :  real   : coefficients in expression for sw single scatt.alb.
! rfucc  :  real   : coefficients in expression for sw assymetry factor

!     -----------------------------------------------------------------

!        * e.c.m.w.f. physics package *

!     j.-j. morcrette       e.c.m.w.f.      89/07/14

!  name     type     purpose
!  ----  :  ----   : -------
!  rtaua :  real     s.w. normalized optical thickness at 0.55 micron
!  rpiza :  real     s.w. single scattering albedo
!  rcga  :  real     s.w. assymetry factor
!  raer  :  real     l.w. absorption coefficients
!     -----------------------------------------------------------------

!        * e.c.m.w.f. physics package *

!     j.-j. morcrette       e.c.m.w.f.      89/07/14

!  name     type     purpose
!  ----  :  ----   : -------
!rtweight:  real     s.w. integrated weight 
! nmpsrtm: integer  : indices for mapping sw[1:6] albedo into srtm[1:14]  
!     -----------------------------------------------------------------
end module yoesw
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

