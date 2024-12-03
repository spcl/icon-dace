! # 1 "ifsrrtm/srtm_kgb25.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb25.f90"
subroutine srtm_kgb25

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 25: 16000-22650 cm-1 (low - h2o; high - nothing)
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
use yoesrta25 , only : ka, sfluxref, rayl, abso3a, abso3b, layreffr, ka_d

!     ------------------------------------------------------------------

implicit none

! kurucz
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 28 "ifsrrtm/srtm_kgb25.f90" 2

if (lhook) call dr_hook('srtm_kgb25',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d
  ka = real(ka_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb25:')
endif

sfluxref = (/ &
 & 42.6858_jprb , 45.7720_jprb, 44.9872_jprb, 45.9662_jprb    , &
 & 46.5458_jprb , 41.6926_jprb, 32.2893_jprb, 24.0928_jprb    , &
 & 16.7686_jprb , 1.86048_jprb, 1.54057_jprb, 1.23503_jprb    , &
 & 0.915085_jprb,0.590099_jprb,0.218622_jprb, 3.21287e-02_jprb /)  

!     rayleigh extinction coefficient at v = 2925 cm-1.
rayl = (/ &
 & 9.81132e-07_jprb,8.25605e-07_jprb,6.71302e-07_jprb,5.53556e-07_jprb,  &
 & 3.97383e-07_jprb,3.68206e-07_jprb,4.42379e-07_jprb,4.57799e-07_jprb, &
 & 4.22683e-07_jprb,3.87113e-07_jprb,3.79810e-07_jprb,3.63192e-07_jprb, &
 & 3.51921e-07_jprb,3.34231e-07_jprb,3.34294e-07_jprb,3.32673e-07_jprb /)  
     
abso3a = (/ &
 & 2.32664e-02_jprb,5.76154e-02_jprb,0.125389_jprb,0.250158_jprb, &
 & 0.378756_jprb   ,0.402196_jprb   ,0.352026_jprb,0.352036_jprb, &
 & 0.386253_jprb   ,0.414598_jprb   ,0.420079_jprb,0.435471_jprb, &
 & 0.445487_jprb   ,0.459549_jprb   ,0.452920_jprb,0.456838_jprb /)  

abso3b = (/      &
 & 1.76917e-02_jprb,4.64185e-02_jprb,1.03640e-01_jprb,0.189469_jprb, &
 & 0.303858_jprb   ,0.400248_jprb   ,0.447357_jprb   ,0.470009_jprb, &
 & 0.498673_jprb   ,0.515696_jprb   ,0.517053_jprb   ,0.517930_jprb, &
 & 0.518345_jprb   ,0.524952_jprb   ,0.508244_jprb   ,0.468981_jprb /)  

layreffr = 2

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
if (lhook) call dr_hook('srtm_kgb25',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb25:error reading file radsrtm")

end subroutine srtm_kgb25
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

