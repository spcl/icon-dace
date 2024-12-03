! # 1 "ifsrrtm/srtm_kgb28.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb28.f90"
subroutine srtm_kgb28

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 28: 38000-50000 cm-1 (low - o3,o2; high - o3,o2)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
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
use yoesrta28 , only : ka, kb, sfluxref, rayl, strrat, layreffr, ka_d, kb_d

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
! # 29 "ifsrrtm/srtm_kgb28.f90" 2

if (lhook) call dr_hook('srtm_kgb28',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb28:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb28:')
endif

sfluxref(:,1) = (/ &
 & 1.06156_jprb    , 0.599910_jprb   , 0.422462_jprb   , 0.400077_jprb   , &
 & 0.282221_jprb   , 0.187893_jprb   , 6.77357e-02_jprb, 3.04572e-02_jprb, &
 & 2.00442e-02_jprb, 2.30786e-03_jprb, 2.08824e-03_jprb, 1.42604e-03_jprb, &
 & 9.67384e-04_jprb, 6.35362e-04_jprb, 1.47727e-04_jprb, 6.87639e-06_jprb /)  
sfluxref(:,2) = (/ &
 & 1.07598_jprb    , 0.585099_jprb   , 0.422852_jprb   , 0.400077_jprb   , &
 & 0.282221_jprb   , 0.187893_jprb   , 6.69686e-02_jprb, 3.09070e-02_jprb, &
 & 2.02400e-02_jprb, 2.47760e-03_jprb, 1.89411e-03_jprb, 1.41122e-03_jprb, &
 & 1.12449e-03_jprb, 5.73505e-04_jprb, 2.04160e-04_jprb, 1.58371e-05_jprb /)  
sfluxref(:,3) = (/ &
 & 0.461647_jprb   , 0.406113_jprb   , 0.332506_jprb   , 0.307508_jprb   , &
 & 0.211167_jprb   , 0.235457_jprb   , 0.495886_jprb   , 0.363921_jprb   , &
 & 0.192700_jprb   , 2.04678e-02_jprb, 1.55407e-02_jprb, 1.03882e-02_jprb, &
 & 1.10778e-02_jprb, 1.00504e-02_jprb, 4.93497e-03_jprb, 5.73410e-04_jprb /)  
sfluxref(:,4) = (/ &
 & 0.132669_jprb   , 0.175058_jprb   , 0.359263_jprb   , 0.388142_jprb   , &
 & 0.350359_jprb   , 0.475892_jprb   , 0.489593_jprb   , 0.408437_jprb   , &
 & 0.221049_jprb   , 1.94514e-02_jprb, 1.54848e-02_jprb, 1.44999e-02_jprb, &
 & 1.44568e-02_jprb, 1.00527e-02_jprb, 4.95897e-03_jprb, 5.73327e-04_jprb /)  
sfluxref(:,5) = (/ &
 & 7.54800e-02_jprb, 0.232246_jprb   , 0.359263_jprb   , 0.388142_jprb   , &
 & 0.350359_jprb   , 0.426317_jprb   , 0.493485_jprb   , 0.432016_jprb   , &
 & 0.239203_jprb   , 1.74951e-02_jprb, 1.74477e-02_jprb, 1.83566e-02_jprb, &
 & 1.44818e-02_jprb, 1.01048e-02_jprb, 4.97487e-03_jprb, 5.66831e-04_jprb /)  

!     rayleigh extinction coefficient at v = ????? cm-1.
rayl = 2.02e-05_jprb

strrat = 6.67029e-07_jprb

layreffr = 58

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

if (lhook) call dr_hook('srtm_kgb28',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb28:error reading file radsrtm")

end subroutine srtm_kgb28
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

