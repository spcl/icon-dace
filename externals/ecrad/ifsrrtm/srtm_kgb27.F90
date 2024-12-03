! # 1 "ifsrrtm/srtm_kgb27.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb27.f90"
subroutine srtm_kgb27

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 16: 29000-38000 cm-1 (low - o3; high - o3)
!     reformatted for f90 by jjmorcrette, ecmwf
!     g.mozdzynski march 2011 read constants from files
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  , only : jprb
use ecradhook   , only : lhook, dr_hook, jphook
use yomlun    , only : nulrad
use yommp0    , only : nproc, myproc
use mpl_module, only : mpl_broadcast
use yomtag    , only : mtagrad
use yoesrta27 , only : ka, kb, sfluxref, rayl, scalekur, layreffr, &
  &  ka_d, kb_d

!     ------------------------------------------------------------------

implicit none

! kurucz
!     the following values were obtained using the "low resolution"
!     version of the kurucz solar source function.  for unknown reasons,
!     the total irradiance in this band differs from the corresponding
!     total in the "high-resolution" version of the kurucz function.
!     therefore, below these values are scaled by the factor scalekur.
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 34 "ifsrrtm/srtm_kgb27.f90" 2

if (lhook) call dr_hook('srtm_kgb27',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb27:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb27:')
endif

sfluxref = (/ &
 & 14.0526_jprb    , 11.4794_jprb    , 8.72590_jprb    , 5.56966_jprb    , &
 & 3.80927_jprb    , 1.57690_jprb    , 1.15099_jprb    , 1.10012_jprb    , &
 & 0.658212_jprb   , 5.86859e-02_jprb, 5.56186e-02_jprb, 4.68040e-02_jprb, &
 & 3.64897e-02_jprb, 3.58053e-02_jprb, 1.38130e-02_jprb, 1.90193e-03_jprb /)  

!     rayleigh extinction coefficient at v = 2925 cm-1.
rayl = (/ &
 & 3.44534e-06_jprb,4.14480e-06_jprb,4.95069e-06_jprb,5.81204e-06_jprb, &
 & 6.69748e-06_jprb,7.56488e-06_jprb,8.36344e-06_jprb,9.04135e-06_jprb, &
 & 9.58324e-06_jprb,9.81542e-06_jprb,9.75119e-06_jprb,9.74533e-06_jprb, &
 & 9.74139e-06_jprb,9.73525e-06_jprb,9.73577e-06_jprb,9.73618e-06_jprb /)  

scalekur = 50.15_jprb/48.37_jprb

layreffr = 32

!     ------------------------------------------------------------------

!     the array ka contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  the first 
!     index in the array, js, runs from 1 to 9, and corresponds to 
!     different values of the binary species parameter.  for instance, 
!     js=1 refers to dry air, js = 2 corresponds to the paramter value 1/8, 
!     js = 3 corresponds to the parameter value 2/8, etc.  the second index
!     in the array, jt, which runs from 1 to 5, corresponds to different
!     temperatures.  more specifically, jt = 3 means that the data are for
!     the reference temperature tref for this  pressure level, jt = 2 refers
!     to tref-15, jt = 1 is for tref-30, jt = 4 is for tref+15, and jt = 5
!     is for tref+30.  the third index, jp, runs from 1 to 13 and refers
!     to the jpth reference pressure level (see taumol.f for these levels
!     in mb).  the fourth index, ig, goes from 1 to 16, and indicates
!     which g-interval the absorption coefficients are for.
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
!     the array kb contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels < ~100mb and temperatures. the first 
!     index in the array, jt, which runs from 1 to 5, corresponds to 
!     different temperatures.  more specifically, jt = 3 means that the 
!     data are for the reference temperature tref for this pressure 
!     level, jt = 2 refers to the temperature tref-15, jt = 1 is for
!     tref-30, jt = 4 is for tref+15, and jt = 5 is for tref+30.  
!     the second index, jp, runs from 13 to 59 and refers to the jpth
!     reference pressure level (see taumol.f for the value of these
!     pressure levels in mb).  the third index, ig, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.
!     -----------------------------------------------------------------
  
     
!     -----------------------------------------------------------------
if (lhook) call dr_hook('srtm_kgb27',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb27:error reading file radsrtm")

end subroutine srtm_kgb27
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

