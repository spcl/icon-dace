! # 1 "ifsrrtm/srtm_kgb29.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb29.f90"
subroutine srtm_kgb29

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 29:   820-2600 cm-1 (low - h2o; high - co2)
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
use yoesrta29 , only : ka, kb, selfref, forref, sfluxref, rayl, &
 & absh2o, absco2, layreffr  , ka_d, kb_d

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
! # 29 "ifsrrtm/srtm_kgb29.f90" 2

if (lhook) call dr_hook('srtm_kgb29',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) ka_d,kb_d
  close(nulrad,err=1000)
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb29:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb29:')
endif

sfluxref = (/ &
 & 1.32880_jprb    , 2.14018_jprb    , 1.97612_jprb    , 1.79000_jprb    , &
 & 1.51242_jprb    , 1.22977_jprb    , 1.06052_jprb    , 0.800996_jprb   , &
 & 0.748053_jprb   , 8.64369e-02_jprb, 7.10675e-02_jprb, 5.62425e-02_jprb, &
 & 4.46988e-02_jprb, 3.07441e-02_jprb, 1.16728e-02_jprb, 1.65573e-03_jprb /)  

absco2 = (/ &
 & 2.90073e-06_jprb, 2.12382e-05_jprb, 1.03032e-04_jprb, 1.86481e-04_jprb, &
 & 4.31997e-04_jprb, 6.08238e-04_jprb, 2.17603e-03_jprb, 4.64479e-02_jprb, &
 & 2.96956_jprb    , 14.9569_jprb    , 28.4831_jprb    , 61.3998_jprb    , &
 & 164.129_jprb    , 832.282_jprb    , 4995.02_jprb    , 12678.1_jprb     /)  
     
absh2o = (/ &
 & 2.99508e-04_jprb, 3.95012e-03_jprb, 1.49316e-02_jprb, 3.24384e-02_jprb, &
 & 6.92879e-02_jprb, 0.123523_jprb   , 0.360985_jprb   , 1.86434_jprb    , &
 & 10.38157_jprb   , 0.214129_jprb   , 0.213914_jprb   , 0.212781_jprb   , &
 & 0.215562_jprb   , 0.218087_jprb   , 0.220918_jprb   , 0.218546_jprb    /)  
     
!     rayleigh extinction coefficient at v = 2200 cm-1.
rayl = 9.30e-11_jprb

layreffr = 49

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


forref(:, 1) = (/ 0.299818e-05_jprb, 0.209282e-05_jprb, 0.988353e-04_jprb, 0.632178e-03_jprb /)
forref(:, 2) = (/ 0.633648e-05_jprb, 0.509214e-04_jprb, 0.650535e-03_jprb, 0.264019e-02_jprb /)
forref(:, 3) = (/ 0.636782e-04_jprb, 0.136577e-03_jprb, 0.166500e-02_jprb, 0.750821e-02_jprb /)
forref(:, 4) = (/ 0.472314e-03_jprb, 0.988296e-03_jprb, 0.585751e-02_jprb, 0.187352e-01_jprb /)
forref(:, 5) = (/ 0.558635e-02_jprb, 0.856489e-02_jprb, 0.157438e-01_jprb, 0.181471e-01_jprb /)
forref(:, 6) = (/ 0.217395e-01_jprb, 0.229156e-01_jprb, 0.230125e-01_jprb, 0.143821e-01_jprb /)
forref(:, 7) = (/ 0.277222e-01_jprb, 0.299252e-01_jprb, 0.208929e-01_jprb, 0.826748e-02_jprb /)
forref(:, 8) = (/ 0.252119e-01_jprb, 0.262911e-01_jprb, 0.187663e-01_jprb, 0.417110e-02_jprb /)
forref(:, 9) = (/ 0.304941e-01_jprb, 0.175545e-01_jprb, 0.971224e-02_jprb, 0.142023e-02_jprb /)
forref(:,10) = (/ 0.327200e-01_jprb, 0.215788e-01_jprb, 0.346831e-02_jprb, 0.157989e-02_jprb /)
forref(:,11) = (/ 0.324955e-01_jprb, 0.228571e-01_jprb, 0.171749e-02_jprb, 0.226853e-02_jprb /)
forref(:,12) = (/ 0.326588e-01_jprb, 0.198544e-01_jprb, 0.532339e-06_jprb, 0.279086e-02_jprb /)
forref(:,13) = (/ 0.345157e-01_jprb, 0.168679e-01_jprb, 0.505361e-06_jprb, 0.276647e-02_jprb /)
forref(:,14) = (/ 0.448765e-01_jprb, 0.123791e-02_jprb, 0.488367e-06_jprb, 0.122245e-02_jprb /)
forref(:,15) = (/ 0.486925e-01_jprb, 0.464371e-06_jprb, 0.464241e-06_jprb, 0.753846e-06_jprb /)
forref(:,16) = (/ 0.530511e-01_jprb, 0.376234e-06_jprb, 0.409824e-06_jprb, 0.470650e-06_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.118069e+00_jprb, 0.713523e-01_jprb, 0.431199e-01_jprb, 0.260584e-01_jprb, 0.157477e-01_jprb, &
 & 0.951675e-02_jprb, 0.575121e-02_jprb, 0.347560e-02_jprb, 0.210039e-02_jprb, 0.126932e-02_jprb /)  
selfref(:, 2) = (/ &
 & 0.137081e-01_jprb, 0.139046e-01_jprb, 0.141040e-01_jprb, 0.143061e-01_jprb, 0.145112e-01_jprb, &
 & 0.147193e-01_jprb, 0.149303e-01_jprb, 0.151443e-01_jprb, 0.153614e-01_jprb, 0.155816e-01_jprb /)  
selfref(:, 3) = (/ &
 & 0.166575e-01_jprb, 0.164916e-01_jprb, 0.163273e-01_jprb, 0.161647e-01_jprb, 0.160037e-01_jprb, &
 & 0.158443e-01_jprb, 0.156864e-01_jprb, 0.155302e-01_jprb, 0.153755e-01_jprb, 0.152224e-01_jprb /)  
selfref(:, 4) = (/ &
 & 0.597379e-01_jprb, 0.509517e-01_jprb, 0.434579e-01_jprb, 0.370662e-01_jprb, 0.316145e-01_jprb, &
 & 0.269647e-01_jprb, 0.229988e-01_jprb, 0.196162e-01_jprb, 0.167311e-01_jprb, 0.142703e-01_jprb /)  
selfref(:, 5) = (/ &
 & 0.227517e+00_jprb, 0.198401e+00_jprb, 0.173011e+00_jprb, 0.150870e+00_jprb, 0.131563e+00_jprb, &
 & 0.114726e+00_jprb, 0.100044e+00_jprb, 0.872415e-01_jprb, 0.760769e-01_jprb, 0.663411e-01_jprb /)  
selfref(:, 6) = (/ &
 & 0.453235e+00_jprb, 0.414848e+00_jprb, 0.379712e+00_jprb, 0.347552e+00_jprb, 0.318116e+00_jprb, &
 & 0.291173e+00_jprb, 0.266512e+00_jprb, 0.243940e+00_jprb, 0.223279e+00_jprb, 0.204368e+00_jprb /)  
selfref(:, 7) = (/ &
 & 0.569263e+00_jprb, 0.516415e+00_jprb, 0.468473e+00_jprb, 0.424982e+00_jprb, 0.385528e+00_jprb, &
 & 0.349737e+00_jprb, 0.317269e+00_jprb, 0.287815e+00_jprb, 0.261095e+00_jprb, 0.236856e+00_jprb /)  
selfref(:, 8) = (/ &
 & 0.490314e+00_jprb, 0.448042e+00_jprb, 0.409413e+00_jprb, 0.374116e+00_jprb, 0.341861e+00_jprb, &
 & 0.312387e+00_jprb, 0.285455e+00_jprb, 0.260844e+00_jprb, 0.238355e+00_jprb, 0.217805e+00_jprb /)  
selfref(:, 9) = (/ &
 & 0.258162e+00_jprb, 0.265085e+00_jprb, 0.272193e+00_jprb, 0.279493e+00_jprb, 0.286988e+00_jprb, &
 & 0.294684e+00_jprb, 0.302586e+00_jprb, 0.310701e+00_jprb, 0.319033e+00_jprb, 0.327588e+00_jprb /)  
selfref(:,10) = (/ &
 & 0.332019e+00_jprb, 0.331902e+00_jprb, 0.331784e+00_jprb, 0.331666e+00_jprb, 0.331549e+00_jprb, &
 & 0.331431e+00_jprb, 0.331314e+00_jprb, 0.331197e+00_jprb, 0.331079e+00_jprb, 0.330962e+00_jprb /)  
selfref(:,11) = (/ &
 & 0.357523e+00_jprb, 0.353154e+00_jprb, 0.348839e+00_jprb, 0.344576e+00_jprb, 0.340366e+00_jprb, &
 & 0.336207e+00_jprb, 0.332099e+00_jprb, 0.328041e+00_jprb, 0.324032e+00_jprb, 0.320073e+00_jprb /)  
selfref(:,12) = (/ &
 & 0.294662e+00_jprb, 0.299043e+00_jprb, 0.303488e+00_jprb, 0.308000e+00_jprb, 0.312579e+00_jprb, &
 & 0.317226e+00_jprb, 0.321941e+00_jprb, 0.326727e+00_jprb, 0.331585e+00_jprb, 0.336514e+00_jprb /)  
selfref(:,13) = (/ &
 & 0.227445e+00_jprb, 0.241545e+00_jprb, 0.256519e+00_jprb, 0.272422e+00_jprb, 0.289311e+00_jprb, &
 & 0.307247e+00_jprb, 0.326294e+00_jprb, 0.346523e+00_jprb, 0.368005e+00_jprb, 0.390820e+00_jprb /)  
selfref(:,14) = (/ &
 & 0.616203e-02_jprb, 0.113523e-01_jprb, 0.209144e-01_jprb, 0.385307e-01_jprb, 0.709852e-01_jprb, &
 & 0.130776e+00_jprb, 0.240929e+00_jprb, 0.443865e+00_jprb, 0.817733e+00_jprb, 0.150651e+01_jprb /)  
selfref(:,15) = (/ &
 & 0.279552e-03_jprb, 0.808472e-03_jprb, 0.233812e-02_jprb, 0.676192e-02_jprb, 0.195557e-01_jprb, &
 & 0.565555e-01_jprb, 0.163560e+00_jprb, 0.473020e+00_jprb, 0.136799e+01_jprb, 0.395626e+01_jprb /)  
selfref(:,16) = (/ &
 & 0.261006e-03_jprb, 0.771043e-03_jprb, 0.227776e-02_jprb, 0.672879e-02_jprb, 0.198777e-01_jprb, &
 & 0.587212e-01_jprb, 0.173470e+00_jprb, 0.512452e+00_jprb, 0.151385e+01_jprb, 0.447209e+01_jprb /)  
     
!     -----------------------------------------------------------------
if (lhook) call dr_hook('srtm_kgb29',1,zhook_handle)
return

1000 continue
call abor1("srtm_kgb29:error closing file radsrtm")
1001 continue
call abor1("srtm_kgb29:error reading file radsrtm")

end subroutine srtm_kgb29
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

