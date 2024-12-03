! # 1 "ifsrrtm/srtm_kgb23.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb23.f90"
subroutine srtm_kgb23

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 16:  8050-12850 cm-1 (low - h2o; high - nothing)
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
use yoesrta23 , only : ka, selfref, forref, sfluxref, rayl, givfac, layreffr  ,&
   &   ka_d

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
! # 29 "ifsrrtm/srtm_kgb23.f90" 2

if (lhook) call dr_hook('srtm_kgb23',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d
  ka = real(ka_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb23:')
endif

sfluxref = (/ &
 & 53.2101_jprb , 51.4143_jprb, 49.3348_jprb, 45.4612_jprb    , &
 & 40.8294_jprb , 35.1801_jprb, 28.6947_jprb, 21.5751_jprb    , &
 & 14.6388_jprb , 1.59111_jprb, 1.31860_jprb, 1.04018_jprb    , &
 & 0.762140_jprb,0.484214_jprb,0.182275_jprb, 2.54948e-02_jprb /)  

!     rayleigh extinction coefficient at all v 
rayl = (/ &
 & 5.94837e-08_jprb,5.70593e-08_jprb,6.27845e-08_jprb,5.56602e-08_jprb, &
 & 5.25571e-08_jprb,4.73388e-08_jprb,4.17466e-08_jprb,3.98097e-08_jprb, &
 & 4.00786e-08_jprb,3.67478e-08_jprb,3.45186e-08_jprb,3.46156e-08_jprb, &
 & 3.32155e-08_jprb,3.23642e-08_jprb,2.72590e-08_jprb,2.96813e-08_jprb /)  

!     average giver et al. correction factor for this band.
givfac = 1.029_jprb

layreffr = 6

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

forref(:, 1) = (/ 0.315770e-07_jprb, 0.671978e-07_jprb, 0.440649e-06_jprb /)
forref(:, 2) = (/ 0.313674e-06_jprb, 0.285252e-06_jprb, 0.421024e-05_jprb /)
forref(:, 3) = (/ 0.135818e-05_jprb, 0.145071e-05_jprb, 0.611285e-05_jprb /)
forref(:, 4) = (/ 0.534065e-05_jprb, 0.586268e-05_jprb, 0.933970e-05_jprb /)
forref(:, 5) = (/ 0.964007e-05_jprb, 0.107110e-04_jprb, 0.104486e-04_jprb /)
forref(:, 6) = (/ 0.302775e-04_jprb, 0.357530e-04_jprb, 0.340724e-04_jprb /)
forref(:, 7) = (/ 0.102437e-03_jprb, 0.108475e-03_jprb, 0.105245e-03_jprb /)
forref(:, 8) = (/ 0.146054e-03_jprb, 0.141490e-03_jprb, 0.133071e-03_jprb /)
forref(:, 9) = (/ 0.163978e-03_jprb, 0.150208e-03_jprb, 0.142864e-03_jprb /)
forref(:,10) = (/ 0.220412e-03_jprb, 0.182943e-03_jprb, 0.150941e-03_jprb /)
forref(:,11) = (/ 0.228877e-03_jprb, 0.197679e-03_jprb, 0.163220e-03_jprb /)
forref(:,12) = (/ 0.234177e-03_jprb, 0.217734e-03_jprb, 0.185038e-03_jprb /)
forref(:,13) = (/ 0.257187e-03_jprb, 0.241570e-03_jprb, 0.221178e-03_jprb /)
forref(:,14) = (/ 0.272455e-03_jprb, 0.270637e-03_jprb, 0.256269e-03_jprb /)
forref(:,15) = (/ 0.339445e-03_jprb, 0.300268e-03_jprb, 0.286574e-03_jprb /)
forref(:,16) = (/ 0.338841e-03_jprb, 0.355428e-03_jprb, 0.353794e-03_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.100945e-04_jprb, 0.801113e-05_jprb, 0.635771e-05_jprb, 0.504554e-05_jprb, 0.400419e-05_jprb, &
 & 0.317777e-05_jprb, 0.252191e-05_jprb, 0.200141e-05_jprb, 0.158834e-05_jprb, 0.126052e-05_jprb /)  
selfref(:, 2) = (/ &
 & 0.107573e-04_jprb, 0.999809e-05_jprb, 0.929245e-05_jprb, 0.863661e-05_jprb, 0.802706e-05_jprb, &
 & 0.746053e-05_jprb, 0.693399e-05_jprb, 0.644460e-05_jprb, 0.598976e-05_jprb, 0.556702e-05_jprb /)  
selfref(:, 3) = (/ &
 & 0.350389e-04_jprb, 0.319234e-04_jprb, 0.290850e-04_jprb, 0.264989e-04_jprb, 0.241428e-04_jprb, &
 & 0.219962e-04_jprb, 0.200404e-04_jprb, 0.182586e-04_jprb, 0.166351e-04_jprb, 0.151560e-04_jprb /)  
selfref(:, 4) = (/ &
 & 0.122993e-03_jprb, 0.110885e-03_jprb, 0.999691e-04_jprb, 0.901277e-04_jprb, 0.812551e-04_jprb, &
 & 0.732559e-04_jprb, 0.660443e-04_jprb, 0.595426e-04_jprb, 0.536809e-04_jprb, 0.483963e-04_jprb /)  
selfref(:, 5) = (/ &
 & 0.206434e-03_jprb, 0.187435e-03_jprb, 0.170185e-03_jprb, 0.154522e-03_jprb, 0.140301e-03_jprb, &
 & 0.127388e-03_jprb, 0.115664e-03_jprb, 0.105019e-03_jprb, 0.953540e-04_jprb, 0.865783e-04_jprb /)  
selfref(:, 6) = (/ &
 & 0.590645e-03_jprb, 0.533109e-03_jprb, 0.481177e-03_jprb, 0.434305e-03_jprb, 0.391998e-03_jprb, &
 & 0.353812e-03_jprb, 0.319346e-03_jprb, 0.288238e-03_jprb, 0.260160e-03_jprb, 0.234817e-03_jprb /)  
selfref(:, 7) = (/ &
 & 0.163029e-02_jprb, 0.148773e-02_jprb, 0.135763e-02_jprb, 0.123891e-02_jprb, 0.113057e-02_jprb, &
 & 0.103170e-02_jprb, 0.941483e-03_jprb, 0.859153e-03_jprb, 0.784023e-03_jprb, 0.715462e-03_jprb /)  
selfref(:, 8) = (/ &
 & 0.204528e-02_jprb, 0.189258e-02_jprb, 0.175128e-02_jprb, 0.162053e-02_jprb, 0.149954e-02_jprb, &
 & 0.138758e-02_jprb, 0.128398e-02_jprb, 0.118812e-02_jprb, 0.109941e-02_jprb, 0.101733e-02_jprb /)  
selfref(:, 9) = (/ &
 & 0.210589e-02_jprb, 0.197078e-02_jprb, 0.184434e-02_jprb, 0.172601e-02_jprb, 0.161528e-02_jprb, &
 & 0.151164e-02_jprb, 0.141466e-02_jprb, 0.132390e-02_jprb, 0.123896e-02_jprb, 0.115947e-02_jprb /)  
selfref(:,10) = (/ &
 & 0.245098e-02_jprb, 0.233745e-02_jprb, 0.222918e-02_jprb, 0.212592e-02_jprb, 0.202745e-02_jprb, &
 & 0.193353e-02_jprb, 0.184397e-02_jprb, 0.175856e-02_jprb, 0.167710e-02_jprb, 0.159941e-02_jprb /)  
selfref(:,11) = (/ &
 & 0.267460e-02_jprb, 0.253325e-02_jprb, 0.239936e-02_jprb, 0.227255e-02_jprb, 0.215244e-02_jprb, &
 & 0.203868e-02_jprb, 0.193093e-02_jprb, 0.182888e-02_jprb, 0.173222e-02_jprb, 0.164067e-02_jprb /)  
selfref(:,12) = (/ &
 & 0.304510e-02_jprb, 0.283919e-02_jprb, 0.264720e-02_jprb, 0.246820e-02_jprb, 0.230130e-02_jprb, &
 & 0.214568e-02_jprb, 0.200059e-02_jprb, 0.186531e-02_jprb, 0.173918e-02_jprb, 0.162157e-02_jprb /)  
selfref(:,13) = (/ &
 & 0.338445e-02_jprb, 0.314719e-02_jprb, 0.292655e-02_jprb, 0.272139e-02_jprb, 0.253060e-02_jprb, &
 & 0.235319e-02_jprb, 0.218822e-02_jprb, 0.203482e-02_jprb, 0.189217e-02_jprb, 0.175952e-02_jprb /)  
selfref(:,14) = (/ &
 & 0.388649e-02_jprb, 0.357018e-02_jprb, 0.327961e-02_jprb, 0.301269e-02_jprb, 0.276750e-02_jprb, &
 & 0.254226e-02_jprb, 0.233535e-02_jprb, 0.214528e-02_jprb, 0.197068e-02_jprb, 0.181029e-02_jprb /)  
selfref(:,15) = (/ &
 & 0.412547e-02_jprb, 0.387413e-02_jprb, 0.363810e-02_jprb, 0.341646e-02_jprb, 0.320831e-02_jprb, &
 & 0.301285e-02_jprb, 0.282930e-02_jprb, 0.265693e-02_jprb, 0.249506e-02_jprb, 0.234305e-02_jprb /)  
selfref(:,16) = (/ &
 & 0.534327e-02_jprb, 0.482967e-02_jprb, 0.436544e-02_jprb, 0.394583e-02_jprb, 0.356655e-02_jprb, &
 & 0.322373e-02_jprb, 0.291387e-02_jprb, 0.263378e-02_jprb, 0.238062e-02_jprb, 0.215179e-02_jprb /)  
     
!     -----------------------------------------------------------------
if (lhook) call dr_hook('srtm_kgb23',1,zhook_handle)
return

1001 continue
call abor1("srtm_kgb23:error reading file radsrtm")

end subroutine srtm_kgb23
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

