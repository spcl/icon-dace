! # 1 "ifsrrtm/srtm_kgb22.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb22.f90"
subroutine srtm_kgb22

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 16:  7700-8050 cm-1 (low - h2o,o2; high - o2)
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
use yoesrta22 , only : ka, kb, selfref, forref, sfluxref, rayl, strrat, layreffr  ,&
  &  ka_d, kb_d

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
! # 30 "ifsrrtm/srtm_kgb22.f90" 2

if (lhook) call dr_hook('srtm_kgb22',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb22:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb22:')
endif

sfluxref(:, 1) = (/ &
 & 3.71641_jprb    ,3.63190_jprb    ,3.44795_jprb    ,3.17936_jprb    , &
 & 2.86071_jprb    ,2.48490_jprb    ,2.02471_jprb    ,1.52475_jprb    , &
 & 1.03811_jprb    ,0.113272_jprb   ,9.37115e-02_jprb,7.38969e-02_jprb, &
 & 5.44713e-02_jprb,3.45905e-02_jprb,1.30293e-02_jprb,1.84198e-03_jprb /)  
sfluxref(:, 2) = (/ &
 & 3.73933_jprb    ,3.60360_jprb    ,3.43370_jprb    ,3.19749_jprb    ,  &
 & 2.87747_jprb    ,2.47926_jprb    ,2.02175_jprb    ,1.52010_jprb    , &
 & 1.03612_jprb    ,0.113265_jprb   ,9.37145e-02_jprb,7.38951e-02_jprb, &
 & 5.44714e-02_jprb,3.45906e-02_jprb,1.30293e-02_jprb,1.84198e-03_jprb /)  
sfluxref(:, 3) = (/ &
 & 3.73889_jprb    ,3.60279_jprb    ,3.43404_jprb    ,3.20560_jprb    , &
 & 2.87367_jprb    ,2.47515_jprb    ,2.02412_jprb    ,1.52315_jprb    , &
 & 1.03146_jprb    ,0.113272_jprb   ,9.36707e-02_jprb,7.39080e-02_jprb, &
 & 5.44598e-02_jprb,3.45906e-02_jprb,1.30293e-02_jprb,1.84198e-03_jprb /)  
sfluxref(:, 4) = (/ &
 & 3.73801_jprb    ,3.60530_jprb    ,3.43659_jprb    ,3.20640_jprb    , &
 & 2.87039_jprb    ,2.47330_jprb    ,2.02428_jprb    ,1.52509_jprb    , &
 & 1.03037_jprb    ,0.112553_jprb   ,9.35352e-02_jprb,7.39675e-02_jprb, &
 & 5.43951e-02_jprb,3.45669e-02_jprb,1.30292e-02_jprb,1.84198e-03_jprb /)  
sfluxref(:, 5) = (/ &
 & 3.73809_jprb    ,3.60996_jprb    ,3.43602_jprb    ,3.20364_jprb    , &
 & 2.87005_jprb    ,2.47343_jprb    ,2.02353_jprb    ,1.52617_jprb    , &
 & 1.03138_jprb    ,0.111172_jprb   ,9.29885e-02_jprb,7.35034e-02_jprb, &
 & 5.42427e-02_jprb,3.45732e-02_jprb,1.30169e-02_jprb,1.84550e-03_jprb /)  
sfluxref(:, 6) = (/ &
 & 3.73872_jprb    ,3.62054_jprb    ,3.42934_jprb    ,3.20110_jprb    , &
 & 2.86886_jprb    ,2.47379_jprb    ,2.02237_jprb    ,1.52754_jprb    ,  &
 & 1.03228_jprb    ,0.111597_jprb   ,9.12252e-02_jprb,7.33115e-02_jprb, &
 & 5.35600e-02_jprb,3.45187e-02_jprb,1.30184e-02_jprb,1.84551e-03_jprb /)  
sfluxref(:, 7) = (/ &
 & 3.73969_jprb    ,3.65461_jprb    ,3.40646_jprb    ,3.19082_jprb    , &
 & 2.86919_jprb    ,2.47289_jprb    ,2.02312_jprb    ,1.52629_jprb    , &
 & 1.03329_jprb    ,0.111611_jprb   ,9.16275e-02_jprb,7.14731e-02_jprb, &
 & 5.31771e-02_jprb,3.44980e-02_jprb,1.30190e-02_jprb,1.84551e-03_jprb /)  
sfluxref(:, 8) = (/ &
 & 3.73995_jprb    ,3.65348_jprb    ,3.43707_jprb    ,3.16351_jprb    , &
 & 2.87003_jprb    ,2.47392_jprb    ,2.02114_jprb    ,1.52548_jprb    ,  &
 & 1.03306_jprb    ,0.111088_jprb   ,9.12422e-02_jprb,7.11146e-02_jprb, &
 & 5.31333e-02_jprb,3.45302e-02_jprb,1.30209e-02_jprb,1.84554e-03_jprb /)  
sfluxref(:, 9) = (/ &
 & 3.73788_jprb    ,3.65004_jprb    ,3.46938_jprb    ,3.15236_jprb    , &
 & 2.86381_jprb    ,2.47393_jprb    ,2.01715_jprb    ,1.52134_jprb    , &
 & 1.03163_jprb    ,0.111259_jprb   ,9.12948e-02_jprb,7.09999e-02_jprb, &
 & 5.31792e-02_jprb,3.44955e-02_jprb,1.30189e-02_jprb,1.84551e-03_jprb /)  

!     rayleigh extinction coefficient at v = 8000 cm-1.
rayl = 1.54e-08_jprb

strrat = 0.022708_jprb

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


forref(:, 1) = (/ 0.351362e-07_jprb, 0.341136e-07_jprb, 0.181317e-06_jprb /)
forref(:, 2) = (/ 0.109648e-06_jprb, 0.344240e-06_jprb, 0.139709e-05_jprb /)
forref(:, 3) = (/ 0.374823e-06_jprb, 0.103424e-05_jprb, 0.188717e-05_jprb /)
forref(:, 4) = (/ 0.580041e-06_jprb, 0.116876e-05_jprb, 0.121183e-05_jprb /)
forref(:, 5) = (/ 0.115608e-05_jprb, 0.148110e-05_jprb, 0.836083e-06_jprb /)
forref(:, 6) = (/ 0.181460e-05_jprb, 0.133313e-05_jprb, 0.500167e-06_jprb /)
forref(:, 7) = (/ 0.199096e-05_jprb, 0.115276e-05_jprb, 0.432994e-06_jprb /)
forref(:, 8) = (/ 0.183730e-05_jprb, 0.122260e-05_jprb, 0.433248e-06_jprb /)
forref(:, 9) = (/ 0.198386e-05_jprb, 0.100130e-05_jprb, 0.269712e-06_jprb /)
forref(:,10) = (/ 0.276382e-05_jprb, 0.749215e-06_jprb, 0.236919e-06_jprb /)
forref(:,11) = (/ 0.298202e-05_jprb, 0.629688e-06_jprb, 0.228388e-06_jprb /)
forref(:,12) = (/ 0.364604e-05_jprb, 0.455336e-06_jprb, 0.206130e-06_jprb /)
forref(:,13) = (/ 0.373339e-05_jprb, 0.245210e-06_jprb, 0.201987e-06_jprb /)
forref(:,14) = (/ 0.480378e-05_jprb, 0.177591e-06_jprb, 0.171458e-06_jprb /)
forref(:,15) = (/ 0.521700e-05_jprb, 0.203358e-06_jprb, 0.189559e-06_jprb /)
forref(:,16) = (/ 0.542717e-05_jprb, 0.219022e-06_jprb, 0.218271e-06_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.538526e-04_jprb, 0.464603e-04_jprb, 0.400828e-04_jprb, 0.345807e-04_jprb, 0.298339e-04_jprb, &
 & 0.257386e-04_jprb, 0.222055e-04_jprb, 0.191574e-04_jprb, 0.165277e-04_jprb, 0.142590e-04_jprb /)  
selfref(:, 2) = (/ &
 & 0.162409e-03_jprb, 0.128347e-03_jprb, 0.101430e-03_jprb, 0.801571e-04_jprb, 0.633460e-04_jprb, &
 & 0.500607e-04_jprb, 0.395616e-04_jprb, 0.312645e-04_jprb, 0.247075e-04_jprb, 0.195257e-04_jprb /)  
selfref(:, 3) = (/ &
 & 0.262882e-03_jprb, 0.212793e-03_jprb, 0.172247e-03_jprb, 0.139427e-03_jprb, 0.112860e-03_jprb, &
 & 0.913557e-04_jprb, 0.739487e-04_jprb, 0.598584e-04_jprb, 0.484529e-04_jprb, 0.392206e-04_jprb /)  
selfref(:, 4) = (/ &
 & 0.242873e-03_jprb, 0.204225e-03_jprb, 0.171726e-03_jprb, 0.144399e-03_jprb, 0.121421e-03_jprb, &
 & 0.102099e-03_jprb, 0.858516e-04_jprb, 0.721899e-04_jprb, 0.607022e-04_jprb, 0.510426e-04_jprb /)  
selfref(:, 5) = (/ &
 & 0.235614e-03_jprb, 0.207814e-03_jprb, 0.183293e-03_jprb, 0.161666e-03_jprb, 0.142591e-03_jprb, &
 & 0.125766e-03_jprb, 0.110927e-03_jprb, 0.978381e-04_jprb, 0.862939e-04_jprb, 0.761119e-04_jprb /)  
selfref(:, 6) = (/ &
 & 0.205508e-03_jprb, 0.190174e-03_jprb, 0.175985e-03_jprb, 0.162854e-03_jprb, 0.150702e-03_jprb, &
 & 0.139458e-03_jprb, 0.129052e-03_jprb, 0.119423e-03_jprb, 0.110513e-03_jprb, 0.102267e-03_jprb /)  
selfref(:, 7) = (/ &
 & 0.185027e-03_jprb, 0.175148e-03_jprb, 0.165796e-03_jprb, 0.156944e-03_jprb, 0.148565e-03_jprb, &
 & 0.140633e-03_jprb, 0.133124e-03_jprb, 0.126016e-03_jprb, 0.119288e-03_jprb, 0.112919e-03_jprb /)  
selfref(:, 8) = (/ &
 & 0.192634e-03_jprb, 0.180192e-03_jprb, 0.168554e-03_jprb, 0.157668e-03_jprb, 0.147484e-03_jprb, &
 & 0.137959e-03_jprb, 0.129048e-03_jprb, 0.120713e-03_jprb, 0.112917e-03_jprb, 0.105624e-03_jprb /)  
selfref(:, 9) = (/ &
 & 0.161632e-03_jprb, 0.155919e-03_jprb, 0.150408e-03_jprb, 0.145092e-03_jprb, 0.139963e-03_jprb, &
 & 0.135016e-03_jprb, 0.130244e-03_jprb, 0.125640e-03_jprb, 0.121199e-03_jprb, 0.116915e-03_jprb /)  
selfref(:,10) = (/ &
 & 0.120880e-03_jprb, 0.125265e-03_jprb, 0.129810e-03_jprb, 0.134520e-03_jprb, 0.139400e-03_jprb, &
 & 0.144458e-03_jprb, 0.149699e-03_jprb, 0.155130e-03_jprb, 0.160758e-03_jprb, 0.166591e-03_jprb /)  
selfref(:,11) = (/ &
 & 0.104705e-03_jprb, 0.111761e-03_jprb, 0.119291e-03_jprb, 0.127330e-03_jprb, 0.135910e-03_jprb, &
 & 0.145068e-03_jprb, 0.154843e-03_jprb, 0.165277e-03_jprb, 0.176414e-03_jprb, 0.188302e-03_jprb /)  
selfref(:,12) = (/ &
 & 0.846335e-04_jprb, 0.951236e-04_jprb, 0.106914e-03_jprb, 0.120166e-03_jprb, 0.135060e-03_jprb, &
 & 0.151800e-03_jprb, 0.170616e-03_jprb, 0.191763e-03_jprb, 0.215532e-03_jprb, 0.242246e-03_jprb /)  
selfref(:,13) = (/ &
 & 0.669754e-04_jprb, 0.781902e-04_jprb, 0.912829e-04_jprb, 0.106568e-03_jprb, 0.124413e-03_jprb, &
 & 0.145245e-03_jprb, 0.169566e-03_jprb, 0.197959e-03_jprb, 0.231107e-03_jprb, 0.269805e-03_jprb /)  
selfref(:,14) = (/ &
 & 0.597091e-04_jprb, 0.722265e-04_jprb, 0.873679e-04_jprb, 0.105684e-03_jprb, 0.127839e-03_jprb, &
 & 0.154639e-03_jprb, 0.187057e-03_jprb, 0.226272e-03_jprb, 0.273707e-03_jprb, 0.331087e-03_jprb /)  
selfref(:,15) = (/ &
 & 0.640410e-04_jprb, 0.771879e-04_jprb, 0.930338e-04_jprb, 0.112133e-03_jprb, 0.135152e-03_jprb, &
 & 0.162897e-03_jprb, 0.196338e-03_jprb, 0.236644e-03_jprb, 0.285225e-03_jprb, 0.343778e-03_jprb /)  
selfref(:,16) = (/ &
 & 0.666420e-04_jprb, 0.801056e-04_jprb, 0.962892e-04_jprb, 0.115742e-03_jprb, 0.139126e-03_jprb, &
 & 0.167233e-03_jprb, 0.201019e-03_jprb, 0.241630e-03_jprb, 0.290446e-03_jprb, 0.349125e-03_jprb /)  

if (lhook) call dr_hook('srtm_kgb22',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb22:error reading file radsrtm")

end subroutine srtm_kgb22
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

