! # 1 "ifsrrtm/rrtm_kgb5.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb5.f90"
subroutine rrtm_kgb5

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation
!     g.mozdzynski march 2011 read constants from files
!     abozzo 201306 updated to rrtmg v4.85
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto5 , only : kao     ,kbo     ,selfrefo   ,forrefo, fracrefao  ,&
 & fracrefbo, ccl4o  , kao_mo3, kao_d, kbo_d

use yommp0    , only : nproc, myproc

!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 30 "ifsrrtm/rrtm_kgb5.f90" 2

if (lhook) call dr_hook('rrtm_kgb5',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb5:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb5:')
endif


! planck fraction mapping level : p = 473.42 mb, t = 259.83
      fracrefao(:, 1) = (/ &
       & 1.4111e-01_jprb,1.4222e-01_jprb,1.3802e-01_jprb,1.3101e-01_jprb,1.2244e-01_jprb,1.0691e-01_jprb, &
       & 8.8703e-02_jprb,6.7130e-02_jprb,4.5509e-02_jprb,4.9866e-03_jprb,4.1214e-03_jprb,3.2557e-03_jprb, &
       & 2.3805e-03_jprb,1.5450e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 2) = (/ &
       & 1.4152e-01_jprb,1.4271e-01_jprb,1.3784e-01_jprb,1.3075e-01_jprb,1.2215e-01_jprb,1.0674e-01_jprb, &
       & 8.8686e-02_jprb,6.7135e-02_jprb,4.5508e-02_jprb,4.9866e-03_jprb,4.1214e-03_jprb,3.2558e-03_jprb, &
       & 2.3805e-03_jprb,1.5450e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 3) = (/ &
       & 1.4159e-01_jprb,1.4300e-01_jprb,1.3781e-01_jprb,1.3094e-01_jprb,1.2192e-01_jprb,1.0661e-01_jprb, &
       & 8.8529e-02_jprb,6.7127e-02_jprb,4.5511e-02_jprb,4.9877e-03_jprb,4.1214e-03_jprb,3.2558e-03_jprb, &
       & 2.3805e-03_jprb,1.5450e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 4) = (/ &
       & 1.4162e-01_jprb,1.4337e-01_jprb,1.3774e-01_jprb,1.3122e-01_jprb,1.2172e-01_jprb,1.0641e-01_jprb, &
       & 8.8384e-02_jprb,6.7056e-02_jprb,4.5514e-02_jprb,4.9880e-03_jprb,4.1214e-03_jprb,3.2557e-03_jprb, &
       & 2.3805e-03_jprb,1.5450e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 5) = (/ &
       & 1.4161e-01_jprb,1.4370e-01_jprb,1.3770e-01_jprb,1.3143e-01_jprb,1.2173e-01_jprb,1.0613e-01_jprb, &
       & 8.8357e-02_jprb,6.6874e-02_jprb,4.5509e-02_jprb,4.9883e-03_jprb,4.1214e-03_jprb,3.2558e-03_jprb, &
       & 2.3804e-03_jprb,1.5450e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 6) = (/ &
       & 1.4154e-01_jprb,1.4405e-01_jprb,1.3771e-01_jprb,1.3169e-01_jprb,1.2166e-01_jprb,1.0603e-01_jprb, &
       & 8.8193e-02_jprb,6.6705e-02_jprb,4.5469e-02_jprb,4.9902e-03_jprb,4.1214e-03_jprb,3.2558e-03_jprb, &
       & 2.3804e-03_jprb,1.5450e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 7) = (/ &
       & 1.4126e-01_jprb,1.4440e-01_jprb,1.3790e-01_jprb,1.3214e-01_jprb,1.2153e-01_jprb,1.0603e-01_jprb, &
       & 8.7908e-02_jprb,6.6612e-02_jprb,4.5269e-02_jprb,4.9900e-03_jprb,4.1256e-03_jprb,3.2558e-03_jprb, &
       & 2.3804e-03_jprb,1.5451e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 8) = (/ &
       & 1.4076e-01_jprb,1.4415e-01_jprb,1.3885e-01_jprb,1.3286e-01_jprb,1.2147e-01_jprb,1.0612e-01_jprb, &
       & 8.7579e-02_jprb,6.6280e-02_jprb,4.4977e-02_jprb,4.9782e-03_jprb,4.1200e-03_jprb,3.2620e-03_jprb, &
       & 2.3820e-03_jprb,1.5452e-03_jprb,5.8423e-04_jprb,8.2275e-05_jprb/)
      fracrefao(:, 9) = (/ &
       & 1.4205e-01_jprb,1.4496e-01_jprb,1.4337e-01_jprb,1.3504e-01_jprb,1.2260e-01_jprb,1.0428e-01_jprb, &
       & 8.4946e-02_jprb,6.3625e-02_jprb,4.2951e-02_jprb,4.7313e-03_jprb,3.9157e-03_jprb,3.0879e-03_jprb, &
       & 2.2666e-03_jprb,1.5193e-03_jprb,5.7469e-04_jprb,8.1674e-05_jprb/)

! planck fraction mapping level : p = 0.2369280 mbar, t = 253.60 k
      fracrefbo(:, 1) = (/ &
       & 1.4075e-01_jprb,1.4196e-01_jprb,1.3833e-01_jprb,1.3345e-01_jprb,1.2234e-01_jprb,1.0718e-01_jprb, &
       & 8.8004e-02_jprb,6.6308e-02_jprb,4.5028e-02_jprb,4.9029e-03_jprb,4.0377e-03_jprb,3.1870e-03_jprb, &
       & 2.3503e-03_jprb,1.5146e-03_jprb,5.7165e-04_jprb,8.2371e-05_jprb/)
      fracrefbo(:, 2) = (/ &
       & 1.4081e-01_jprb,1.4225e-01_jprb,1.3890e-01_jprb,1.3410e-01_jprb,1.2254e-01_jprb,1.0680e-01_jprb, &
       & 8.7391e-02_jprb,6.5819e-02_jprb,4.4725e-02_jprb,4.9121e-03_jprb,4.0420e-03_jprb,3.1869e-03_jprb, &
       & 2.3504e-03_jprb,1.5146e-03_jprb,5.7165e-04_jprb,8.2371e-05_jprb/)
      fracrefbo(:, 3) = (/ &
       & 1.4087e-01_jprb,1.4227e-01_jprb,1.3920e-01_jprb,1.3395e-01_jprb,1.2270e-01_jprb,1.0694e-01_jprb, &
       & 8.7229e-02_jprb,6.5653e-02_jprb,4.4554e-02_jprb,4.8797e-03_jprb,4.0460e-03_jprb,3.1939e-03_jprb, &
       & 2.3505e-03_jprb,1.5146e-03_jprb,5.7165e-04_jprb,8.1910e-05_jprb/)
      fracrefbo(:, 4) = (/ &
       & 1.4089e-01_jprb,1.4238e-01_jprb,1.3956e-01_jprb,1.3379e-01_jprb,1.2284e-01_jprb,1.0688e-01_jprb, &
       & 8.7192e-02_jprb,6.5490e-02_jprb,4.4390e-02_jprb,4.8395e-03_jprb,4.0173e-03_jprb,3.2070e-03_jprb, &
       & 2.3559e-03_jprb,1.5146e-03_jprb,5.7165e-04_jprb,8.2371e-05_jprb/)
      fracrefbo(:, 5) = (/ &
       & 1.4091e-01_jprb,1.4417e-01_jprb,1.4194e-01_jprb,1.3457e-01_jprb,1.2167e-01_jprb,1.0551e-01_jprb, &
       & 8.6450e-02_jprb,6.4889e-02_jprb,4.3584e-02_jprb,4.7551e-03_jprb,3.9509e-03_jprb,3.1374e-03_jprb, &
       & 2.3226e-03_jprb,1.4942e-03_jprb,5.7545e-04_jprb,8.0887e-05_jprb/)



ccl4o( :) = (/&
 & 26.1407_jprb,  53.9776_jprb,  63.8085_jprb,  36.1701_jprb,&
 & 15.4099_jprb, 10.23116_jprb,  4.82948_jprb,  5.03836_jprb,&
 & 1.75558_jprb,  0.0_jprb     ,  0.0_jprb     ,  0.0_jprb     ,&
 & 0.0_jprb     ,  0.0_jprb     ,  0.0_jprb     ,  0.0_jprb      /)  


!     ------------------------------------------------------------------

!     the array kao contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels > ~100mb, temperatures, and ratios
!     of water vapor to co2.  the first index in the array, js, runs
!     from 1 to 9 and corresponds to different water vapor to co2 ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = h2o/(h20 + (rat) * co2), where rat is the ratio of the integrated
!     line strength in the band of co2 to that of h2o.  for instance,
!     js=1 refers to dry air (eta = 0), js = 9 corresponds to eta = 1.0.
!     the 2nd index in the array, jt, which runs from 1 to 5, corresponds 
!     to different temperatures.  more specifically, jt = 3 means that the 
!     data are for the reference temperature tref for this  pressure 
!     level, jt = 2 refers to the temperature tref-15, 
!     jt = 1 is for tref-30, jt = 4 is for tref+15, and jt = 5
!     is for tref+30.  the third index, jp, runs from 1 to 13 and refers
!     to the reference pressure level (e.g. jp = 1 is for a
!     pressure of 1053.63 mb).  the fourth index, ig, goes from 1 to 16,
!     and tells us which g-interval the absorption coefficients are for.



!     the array kb contains absorption coefs for each of the 16 g-intervals
!     for a range of pressure levels  < ~100mb, temperatures, and ratios
!     of o3 to co2.  the first index in the array, js, runs from 1 to 5, 
!     and corresponds to different o3 to co2 ratios, as expressed through 
!     the binary species parameter eta, defined as eta = o3/(o3+rat*co2), 
!     where rat is the ratio of the integrated line strength in the band 
!     of co2 to that of o3.  for instance, js=1 refers to no o3 (eta = 0) 
!     and js = 5 corresponds to eta = 1.0.  the second index, jt, which
!     runs from 1 to 5, corresponds to different temperatures.  more 
!     specifically, jt = 3 means that the data are for the corresponding 
!     reference temperature tref for this  pressure level, jt = 2 refers 
!     to the tref-15, jt = 1 is for tref-30, jt = 4 is for tref+15, and
!     jt = 5 is for tref+30.  the third index, jp, runs from 13 to 59 and
!     refers to the corresponding pressure level in pref (e.g. jp = 13 is
!     for a pressure of 95.5835 mb).  the fourth index, ig, goes from 1 to
!     16, and tells us which g-interval the absorption coefficients are for.



!     the array kao_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   the first index in the array, js, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference mls column amount value of gas 1 
!     to that of gas2.  the second index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the third index 
!     runs over the g-channel (1 to 16).

      kao_mo3( 1, :, 1) = (/ &
     & 9.31040e-03_jprb, 1.01286e-02_jprb, 1.10186e-02_jprb, 1.19869e-02_jprb, 1.30403e-02_jprb, &
     & 1.41862e-02_jprb, 1.54328e-02_jprb, 1.67890e-02_jprb, 1.82644e-02_jprb, 1.98694e-02_jprb, &
     & 2.16154e-02_jprb, 2.35149e-02_jprb, 2.55813e-02_jprb, 2.78293e-02_jprb, 3.02749e-02_jprb, &
     & 3.29353e-02_jprb, 3.58295e-02_jprb, 3.89781e-02_jprb, 4.24034e-02_jprb/)
      kao_mo3( 2, :, 1) = (/ &
     & 1.11200e-02_jprb, 1.20461e-02_jprb, 1.30493e-02_jprb, 1.41360e-02_jprb, 1.53133e-02_jprb, &
     & 1.65886e-02_jprb, 1.79701e-02_jprb, 1.94666e-02_jprb, 2.10878e-02_jprb, 2.28440e-02_jprb, &
     & 2.47465e-02_jprb, 2.68074e-02_jprb, 2.90399e-02_jprb, 3.14583e-02_jprb, 3.40782e-02_jprb, &
     & 3.69162e-02_jprb, 3.99907e-02_jprb, 4.33211e-02_jprb, 4.69289e-02_jprb/)
      kao_mo3( 3, :, 1) = (/ &
     & 1.21630e-02_jprb, 1.31401e-02_jprb, 1.41956e-02_jprb, 1.53359e-02_jprb, 1.65679e-02_jprb, &
     & 1.78988e-02_jprb, 1.93366e-02_jprb, 2.08899e-02_jprb, 2.25680e-02_jprb, 2.43808e-02_jprb, &
     & 2.63394e-02_jprb, 2.84552e-02_jprb, 3.07410e-02_jprb, 3.32104e-02_jprb, 3.58782e-02_jprb, &
     & 3.87603e-02_jprb, 4.18739e-02_jprb, 4.52377e-02_jprb, 4.88716e-02_jprb/)
      kao_mo3( 4, :, 1) = (/ &
     & 1.26231e-02_jprb, 1.36243e-02_jprb, 1.47049e-02_jprb, 1.58713e-02_jprb, 1.71301e-02_jprb, &
     & 1.84888e-02_jprb, 1.99553e-02_jprb, 2.15380e-02_jprb, 2.32463e-02_jprb, 2.50901e-02_jprb, &
     & 2.70801e-02_jprb, 2.92280e-02_jprb, 3.15463e-02_jprb, 3.40484e-02_jprb, 3.67489e-02_jprb, &
     & 3.96637e-02_jprb, 4.28097e-02_jprb, 4.62051e-02_jprb, 4.98699e-02_jprb/)
      kao_mo3( 5, :, 1) = (/ &
     & 1.33345e-02_jprb, 1.43736e-02_jprb, 1.54938e-02_jprb, 1.67012e-02_jprb, 1.80027e-02_jprb, &
     & 1.94057e-02_jprb, 2.09180e-02_jprb, 2.25481e-02_jprb, 2.43053e-02_jprb, 2.61994e-02_jprb, &
     & 2.82411e-02_jprb, 3.04419e-02_jprb, 3.28142e-02_jprb, 3.53714e-02_jprb, 3.81279e-02_jprb, &
     & 4.10992e-02_jprb, 4.43021e-02_jprb, 4.77545e-02_jprb, 5.14760e-02_jprb/)
      kao_mo3( 6, :, 1) = (/ &
     & 1.43294e-02_jprb, 1.54133e-02_jprb, 1.65791e-02_jprb, 1.78331e-02_jprb, 1.91819e-02_jprb, &
     & 2.06328e-02_jprb, 2.21935e-02_jprb, 2.38721e-02_jprb, 2.56778e-02_jprb, 2.76200e-02_jprb, &
     & 2.97091e-02_jprb, 3.19562e-02_jprb, 3.43733e-02_jprb, 3.69732e-02_jprb, 3.97698e-02_jprb, &
     & 4.27779e-02_jprb, 4.60136e-02_jprb, 4.94939e-02_jprb, 5.32375e-02_jprb/)
      kao_mo3( 7, :, 1) = (/ &
     & 1.48298e-02_jprb, 1.59503e-02_jprb, 1.71554e-02_jprb, 1.84517e-02_jprb, 1.98458e-02_jprb, &
     & 2.13453e-02_jprb, 2.29581e-02_jprb, 2.46928e-02_jprb, 2.65585e-02_jprb, 2.85652e-02_jprb, &
     & 3.07235e-02_jprb, 3.30449e-02_jprb, 3.55417e-02_jprb, 3.82272e-02_jprb, 4.11155e-02_jprb, &
     & 4.42221e-02_jprb, 4.75634e-02_jprb, 5.11572e-02_jprb, 5.50225e-02_jprb/)
      kao_mo3( 8, :, 1) = (/ &
     & 1.41792e-02_jprb, 1.53141e-02_jprb, 1.65398e-02_jprb, 1.78637e-02_jprb, 1.92935e-02_jprb, &
     & 2.08378e-02_jprb, 2.25057e-02_jprb, 2.43071e-02_jprb, 2.62526e-02_jprb, 2.83539e-02_jprb, &
     & 3.06234e-02_jprb, 3.30745e-02_jprb, 3.57218e-02_jprb, 3.85810e-02_jprb, 4.16690e-02_jprb, &
     & 4.50042e-02_jprb, 4.86064e-02_jprb, 5.24969e-02_jprb, 5.66988e-02_jprb/)
      kao_mo3( 9, :, 1) = (/ &
     & 8.82784e-03_jprb, 9.48321e-03_jprb, 1.01872e-02_jprb, 1.09435e-02_jprb, 1.17560e-02_jprb, &
     & 1.26287e-02_jprb, 1.35662e-02_jprb, 1.45734e-02_jprb, 1.56553e-02_jprb, 1.68175e-02_jprb, &
     & 1.80660e-02_jprb, 1.94072e-02_jprb, 2.08480e-02_jprb, 2.23958e-02_jprb, 2.40584e-02_jprb, &
     & 2.58445e-02_jprb, 2.77631e-02_jprb, 2.98242e-02_jprb, 3.20383e-02_jprb/)
      kao_mo3( 1, :, 2) = (/ &
     & 4.28238e-02_jprb, 4.51015e-02_jprb, 4.75003e-02_jprb, 5.00266e-02_jprb, 5.26873e-02_jprb, &
     & 5.54896e-02_jprb, 5.84409e-02_jprb, 6.15491e-02_jprb, 6.48227e-02_jprb, 6.82704e-02_jprb, &
     & 7.19014e-02_jprb, 7.57256e-02_jprb, 7.97532e-02_jprb, 8.39949e-02_jprb, 8.84623e-02_jprb, &
     & 9.31673e-02_jprb, 9.81225e-02_jprb, 1.03341e-01_jprb, 1.08838e-01_jprb/)
      kao_mo3( 2, :, 2) = (/ &
     & 4.83672e-02_jprb, 5.07219e-02_jprb, 5.31911e-02_jprb, 5.57806e-02_jprb, 5.84962e-02_jprb, &
     & 6.13440e-02_jprb, 6.43303e-02_jprb, 6.74621e-02_jprb, 7.07464e-02_jprb, 7.41905e-02_jprb, &
     & 7.78023e-02_jprb, 8.15899e-02_jprb, 8.55619e-02_jprb, 8.97273e-02_jprb, 9.40955e-02_jprb, &
     & 9.86763e-02_jprb, 1.03480e-01_jprb, 1.08518e-01_jprb, 1.13801e-01_jprb/)
      kao_mo3( 3, :, 2) = (/ &
     & 5.24315e-02_jprb, 5.48650e-02_jprb, 5.74115e-02_jprb, 6.00762e-02_jprb, 6.28645e-02_jprb, &
     & 6.57822e-02_jprb, 6.88354e-02_jprb, 7.20302e-02_jprb, 7.53734e-02_jprb, 7.88717e-02_jprb, &
     & 8.25324e-02_jprb, 8.63630e-02_jprb, 9.03714e-02_jprb, 9.45658e-02_jprb, 9.89549e-02_jprb, &
     & 1.03548e-01_jprb, 1.08354e-01_jprb, 1.13383e-01_jprb, 1.18645e-01_jprb/)
      kao_mo3( 4, :, 2) = (/ &
     & 5.65191e-02_jprb, 5.90383e-02_jprb, 6.16699e-02_jprb, 6.44187e-02_jprb, 6.72901e-02_jprb, &
     & 7.02894e-02_jprb, 7.34224e-02_jprb, 7.66951e-02_jprb, 8.01137e-02_jprb, 8.36846e-02_jprb, &
     & 8.74147e-02_jprb, 9.13111e-02_jprb, 9.53812e-02_jprb, 9.96326e-02_jprb, 1.04074e-01_jprb, &
     & 1.08712e-01_jprb, 1.13558e-01_jprb, 1.18620e-01_jprb, 1.23907e-01_jprb/)
      kao_mo3( 5, :, 2) = (/ &
     & 6.03171e-02_jprb, 6.29114e-02_jprb, 6.56172e-02_jprb, 6.84394e-02_jprb, 7.13830e-02_jprb, &
     & 7.44532e-02_jprb, 7.76555e-02_jprb, 8.09955e-02_jprb, 8.44791e-02_jprb, 8.81125e-02_jprb, &
     & 9.19023e-02_jprb, 9.58550e-02_jprb, 9.99778e-02_jprb, 1.04278e-01_jprb, 1.08763e-01_jprb, &
     & 1.13441e-01_jprb, 1.18320e-01_jprb, 1.23409e-01_jprb, 1.28717e-01_jprb/)
      kao_mo3( 6, :, 2) = (/ &
     & 6.51092e-02_jprb, 6.77827e-02_jprb, 7.05660e-02_jprb, 7.34635e-02_jprb, 7.64801e-02_jprb, &
     & 7.96204e-02_jprb, 8.28898e-02_jprb, 8.62934e-02_jprb, 8.98367e-02_jprb, 9.35255e-02_jprb, &
     & 9.73658e-02_jprb, 1.01364e-01_jprb, 1.05526e-01_jprb, 1.09859e-01_jprb, 1.14370e-01_jprb, &
     & 1.19066e-01_jprb, 1.23955e-01_jprb, 1.29045e-01_jprb, 1.34344e-01_jprb/)
      kao_mo3( 7, :, 2) = (/ &
     & 7.09653e-02_jprb, 7.37378e-02_jprb, 7.66187e-02_jprb, 7.96121e-02_jprb, 8.27225e-02_jprb, &
     & 8.59543e-02_jprb, 8.93125e-02_jprb, 9.28018e-02_jprb, 9.64275e-02_jprb, 1.00195e-01_jprb, &
     & 1.04109e-01_jprb, 1.08177e-01_jprb, 1.12403e-01_jprb, 1.16795e-01_jprb, 1.21358e-01_jprb, &
     & 1.26099e-01_jprb, 1.31026e-01_jprb, 1.36145e-01_jprb, 1.41464e-01_jprb/)
      kao_mo3( 8, :, 2) = (/ &
     & 7.69193e-02_jprb, 7.97926e-02_jprb, 8.27733e-02_jprb, 8.58653e-02_jprb, 8.90728e-02_jprb, &
     & 9.24002e-02_jprb, 9.58518e-02_jprb, 9.94324e-02_jprb, 1.03147e-01_jprb, 1.07000e-01_jprb, &
     & 1.10997e-01_jprb, 1.15143e-01_jprb, 1.19444e-01_jprb, 1.23906e-01_jprb, 1.28535e-01_jprb, &
     & 1.33336e-01_jprb, 1.38317e-01_jprb, 1.43484e-01_jprb, 1.48844e-01_jprb/)
      kao_mo3( 9, :, 2) = (/ &
     & 4.57962e-02_jprb, 4.76027e-02_jprb, 4.94805e-02_jprb, 5.14323e-02_jprb, 5.34611e-02_jprb, &
     & 5.55700e-02_jprb, 5.77620e-02_jprb, 6.00405e-02_jprb, 6.24089e-02_jprb, 6.48707e-02_jprb, &
     & 6.74296e-02_jprb, 7.00895e-02_jprb, 7.28542e-02_jprb, 7.57281e-02_jprb, 7.87153e-02_jprb, &
     & 8.18203e-02_jprb, 8.50478e-02_jprb, 8.84027e-02_jprb, 9.18898e-02_jprb/)
      kao_mo3( 1, :, 3) = (/ &
     & 1.12607e-01_jprb, 1.16047e-01_jprb, 1.19591e-01_jprb, 1.23244e-01_jprb, 1.27009e-01_jprb, &
     & 1.30888e-01_jprb, 1.34886e-01_jprb, 1.39006e-01_jprb, 1.43252e-01_jprb, 1.47628e-01_jprb, &
     & 1.52137e-01_jprb, 1.56785e-01_jprb, 1.61574e-01_jprb, 1.66509e-01_jprb, 1.71595e-01_jprb, &
     & 1.76836e-01_jprb, 1.82238e-01_jprb, 1.87804e-01_jprb, 1.93541e-01_jprb/)
      kao_mo3( 2, :, 3) = (/ &
     & 1.14531e-01_jprb, 1.17850e-01_jprb, 1.21266e-01_jprb, 1.24781e-01_jprb, 1.28397e-01_jprb, &
     & 1.32119e-01_jprb, 1.35948e-01_jprb, 1.39888e-01_jprb, 1.43943e-01_jprb, 1.48115e-01_jprb, &
     & 1.52407e-01_jprb, 1.56825e-01_jprb, 1.61370e-01_jprb, 1.66047e-01_jprb, 1.70860e-01_jprb, &
     & 1.75812e-01_jprb, 1.80907e-01_jprb, 1.86150e-01_jprb, 1.91546e-01_jprb/)
      kao_mo3( 3, :, 3) = (/ &
     & 1.13986e-01_jprb, 1.17222e-01_jprb, 1.20551e-01_jprb, 1.23974e-01_jprb, 1.27494e-01_jprb, &
     & 1.31114e-01_jprb, 1.34837e-01_jprb, 1.38666e-01_jprb, 1.42604e-01_jprb, 1.46653e-01_jprb, &
     & 1.50817e-01_jprb, 1.55099e-01_jprb, 1.59503e-01_jprb, 1.64032e-01_jprb, 1.68690e-01_jprb, &
     & 1.73480e-01_jprb, 1.78406e-01_jprb, 1.83472e-01_jprb, 1.88682e-01_jprb/)
      kao_mo3( 4, :, 3) = (/ &
     & 1.13713e-01_jprb, 1.16892e-01_jprb, 1.20160e-01_jprb, 1.23519e-01_jprb, 1.26972e-01_jprb, &
     & 1.30522e-01_jprb, 1.34171e-01_jprb, 1.37922e-01_jprb, 1.41778e-01_jprb, 1.45742e-01_jprb, &
     & 1.49817e-01_jprb, 1.54005e-01_jprb, 1.58311e-01_jprb, 1.62737e-01_jprb, 1.67287e-01_jprb, &
     & 1.71964e-01_jprb, 1.76771e-01_jprb, 1.81714e-01_jprb, 1.86794e-01_jprb/)
      kao_mo3( 5, :, 3) = (/ &
     & 1.12321e-01_jprb, 1.15413e-01_jprb, 1.18591e-01_jprb, 1.21856e-01_jprb, 1.25211e-01_jprb, &
     & 1.28658e-01_jprb, 1.32200e-01_jprb, 1.35840e-01_jprb, 1.39580e-01_jprb, 1.43423e-01_jprb, &
     & 1.47372e-01_jprb, 1.51429e-01_jprb, 1.55599e-01_jprb, 1.59883e-01_jprb, 1.64284e-01_jprb, &
     & 1.68808e-01_jprb, 1.73455e-01_jprb, 1.78231e-01_jprb, 1.83138e-01_jprb/)
      kao_mo3( 6, :, 3) = (/ &
     & 1.14158e-01_jprb, 1.17218e-01_jprb, 1.20360e-01_jprb, 1.23586e-01_jprb, 1.26899e-01_jprb, &
     & 1.30300e-01_jprb, 1.33793e-01_jprb, 1.37379e-01_jprb, 1.41061e-01_jprb, 1.44842e-01_jprb, &
     & 1.48724e-01_jprb, 1.52711e-01_jprb, 1.56804e-01_jprb, 1.61007e-01_jprb, 1.65322e-01_jprb, &
     & 1.69754e-01_jprb, 1.74304e-01_jprb, 1.78976e-01_jprb, 1.83773e-01_jprb/)
      kao_mo3( 7, :, 3) = (/ &
     & 1.21015e-01_jprb, 1.23989e-01_jprb, 1.27036e-01_jprb, 1.30157e-01_jprb, 1.33355e-01_jprb, &
     & 1.36632e-01_jprb, 1.39990e-01_jprb, 1.43429e-01_jprb, 1.46954e-01_jprb, 1.50565e-01_jprb, &
     & 1.54264e-01_jprb, 1.58055e-01_jprb, 1.61939e-01_jprb, 1.65918e-01_jprb, 1.69995e-01_jprb, &
     & 1.74172e-01_jprb, 1.78452e-01_jprb, 1.82836e-01_jprb, 1.87329e-01_jprb/)
      kao_mo3( 8, :, 3) = (/ &
     & 1.33952e-01_jprb, 1.36939e-01_jprb, 1.39992e-01_jprb, 1.43114e-01_jprb, 1.46305e-01_jprb, &
     & 1.49567e-01_jprb, 1.52902e-01_jprb, 1.56311e-01_jprb, 1.59797e-01_jprb, 1.63360e-01_jprb, &
     & 1.67002e-01_jprb, 1.70726e-01_jprb, 1.74533e-01_jprb, 1.78424e-01_jprb, 1.82403e-01_jprb, &
     & 1.86470e-01_jprb, 1.90627e-01_jprb, 1.94878e-01_jprb, 1.99223e-01_jprb/)
      kao_mo3( 9, :, 3) = (/ &
     & 1.01003e-01_jprb, 1.03713e-01_jprb, 1.06495e-01_jprb, 1.09352e-01_jprb, 1.12285e-01_jprb, &
     & 1.15297e-01_jprb, 1.18390e-01_jprb, 1.21566e-01_jprb, 1.24827e-01_jprb, 1.28176e-01_jprb, &
     & 1.31614e-01_jprb, 1.35145e-01_jprb, 1.38770e-01_jprb, 1.42493e-01_jprb, 1.46315e-01_jprb, &
     & 1.50240e-01_jprb, 1.54271e-01_jprb, 1.58409e-01_jprb, 1.62659e-01_jprb/)
      kao_mo3( 1, :, 4) = (/ &
     & 2.35597e-01_jprb, 2.37975e-01_jprb, 2.40376e-01_jprb, 2.42802e-01_jprb, 2.45253e-01_jprb, &
     & 2.47728e-01_jprb, 2.50228e-01_jprb, 2.52753e-01_jprb, 2.55304e-01_jprb, 2.57881e-01_jprb, &
     & 2.60483e-01_jprb, 2.63112e-01_jprb, 2.65767e-01_jprb, 2.68450e-01_jprb, 2.71159e-01_jprb, &
     & 2.73895e-01_jprb, 2.76660e-01_jprb, 2.79452e-01_jprb, 2.82272e-01_jprb/)
      kao_mo3( 2, :, 4) = (/ &
     & 2.27965e-01_jprb, 2.30334e-01_jprb, 2.32728e-01_jprb, 2.35146e-01_jprb, 2.37590e-01_jprb, &
     & 2.40059e-01_jprb, 2.42554e-01_jprb, 2.45075e-01_jprb, 2.47621e-01_jprb, 2.50195e-01_jprb, &
     & 2.52795e-01_jprb, 2.55422e-01_jprb, 2.58077e-01_jprb, 2.60759e-01_jprb, 2.63468e-01_jprb, &
     & 2.66206e-01_jprb, 2.68973e-01_jprb, 2.71768e-01_jprb, 2.74593e-01_jprb/)
      kao_mo3( 3, :, 4) = (/ &
     & 2.25956e-01_jprb, 2.28277e-01_jprb, 2.30622e-01_jprb, 2.32991e-01_jprb, 2.35384e-01_jprb, &
     & 2.37802e-01_jprb, 2.40244e-01_jprb, 2.42712e-01_jprb, 2.45205e-01_jprb, 2.47724e-01_jprb, &
     & 2.50268e-01_jprb, 2.52839e-01_jprb, 2.55436e-01_jprb, 2.58060e-01_jprb, 2.60711e-01_jprb, &
     & 2.63389e-01_jprb, 2.66094e-01_jprb, 2.68827e-01_jprb, 2.71589e-01_jprb/)
      kao_mo3( 4, :, 4) = (/ &
     & 2.28371e-01_jprb, 2.30595e-01_jprb, 2.32840e-01_jprb, 2.35107e-01_jprb, 2.37397e-01_jprb, &
     & 2.39708e-01_jprb, 2.42042e-01_jprb, 2.44399e-01_jprb, 2.46779e-01_jprb, 2.49182e-01_jprb, &
     & 2.51608e-01_jprb, 2.54058e-01_jprb, 2.56532e-01_jprb, 2.59030e-01_jprb, 2.61552e-01_jprb, &
     & 2.64099e-01_jprb, 2.66671e-01_jprb, 2.69267e-01_jprb, 2.71889e-01_jprb/)
      kao_mo3( 5, :, 4) = (/ &
     & 2.42563e-01_jprb, 2.44620e-01_jprb, 2.46695e-01_jprb, 2.48787e-01_jprb, 2.50897e-01_jprb, &
     & 2.53024e-01_jprb, 2.55170e-01_jprb, 2.57334e-01_jprb, 2.59516e-01_jprb, 2.61717e-01_jprb, &
     & 2.63936e-01_jprb, 2.66174e-01_jprb, 2.68431e-01_jprb, 2.70708e-01_jprb, 2.73003e-01_jprb, &
     & 2.75318e-01_jprb, 2.77653e-01_jprb, 2.80008e-01_jprb, 2.82382e-01_jprb/)
      kao_mo3( 6, :, 4) = (/ &
     & 2.54052e-01_jprb, 2.56017e-01_jprb, 2.57997e-01_jprb, 2.59992e-01_jprb, 2.62003e-01_jprb, &
     & 2.64029e-01_jprb, 2.66071e-01_jprb, 2.68129e-01_jprb, 2.70203e-01_jprb, 2.72293e-01_jprb, &
     & 2.74398e-01_jprb, 2.76521e-01_jprb, 2.78659e-01_jprb, 2.80814e-01_jprb, 2.82986e-01_jprb, &
     & 2.85175e-01_jprb, 2.87380e-01_jprb, 2.89603e-01_jprb, 2.91842e-01_jprb/)
      kao_mo3( 7, :, 4) = (/ &
     & 2.54061e-01_jprb, 2.55982e-01_jprb, 2.57917e-01_jprb, 2.59867e-01_jprb, 2.61832e-01_jprb, &
     & 2.63811e-01_jprb, 2.65806e-01_jprb, 2.67815e-01_jprb, 2.69840e-01_jprb, 2.71880e-01_jprb, &
     & 2.73936e-01_jprb, 2.76007e-01_jprb, 2.78093e-01_jprb, 2.80196e-01_jprb, 2.82314e-01_jprb, &
     & 2.84449e-01_jprb, 2.86599e-01_jprb, 2.88766e-01_jprb, 2.90949e-01_jprb/)
      kao_mo3( 8, :, 4) = (/ &
     & 2.72482e-01_jprb, 2.73916e-01_jprb, 2.75358e-01_jprb, 2.76807e-01_jprb, 2.78264e-01_jprb, &
     & 2.79729e-01_jprb, 2.81201e-01_jprb, 2.82681e-01_jprb, 2.84169e-01_jprb, 2.85665e-01_jprb, &
     & 2.87168e-01_jprb, 2.88680e-01_jprb, 2.90199e-01_jprb, 2.91726e-01_jprb, 2.93262e-01_jprb, &
     & 2.94805e-01_jprb, 2.96357e-01_jprb, 2.97917e-01_jprb, 2.99485e-01_jprb/)
      kao_mo3( 9, :, 4) = (/ &
     & 1.93414e-01_jprb, 1.95498e-01_jprb, 1.97605e-01_jprb, 1.99734e-01_jprb, 2.01886e-01_jprb, &
     & 2.04062e-01_jprb, 2.06261e-01_jprb, 2.08483e-01_jprb, 2.10730e-01_jprb, 2.13001e-01_jprb, &
     & 2.15296e-01_jprb, 2.17616e-01_jprb, 2.19961e-01_jprb, 2.22331e-01_jprb, 2.24727e-01_jprb, &
     & 2.27148e-01_jprb, 2.29596e-01_jprb, 2.32070e-01_jprb, 2.34571e-01_jprb/)
      kao_mo3( 1, :, 5) = (/ &
     & 5.30785e-01_jprb, 5.30477e-01_jprb, 5.30169e-01_jprb, 5.29861e-01_jprb, 5.29553e-01_jprb, &
     & 5.29246e-01_jprb, 5.28938e-01_jprb, 5.28631e-01_jprb, 5.28324e-01_jprb, 5.28017e-01_jprb, &
     & 5.27711e-01_jprb, 5.27404e-01_jprb, 5.27098e-01_jprb, 5.26792e-01_jprb, 5.26486e-01_jprb, &
     & 5.26180e-01_jprb, 5.25875e-01_jprb, 5.25569e-01_jprb, 5.25264e-01_jprb/)
      kao_mo3( 2, :, 5) = (/ &
     & 5.33406e-01_jprb, 5.32997e-01_jprb, 5.32587e-01_jprb, 5.32178e-01_jprb, 5.31769e-01_jprb, &
     & 5.31360e-01_jprb, 5.30952e-01_jprb, 5.30544e-01_jprb, 5.30137e-01_jprb, 5.29729e-01_jprb, &
     & 5.29322e-01_jprb, 5.28916e-01_jprb, 5.28509e-01_jprb, 5.28103e-01_jprb, 5.27697e-01_jprb, &
     & 5.27292e-01_jprb, 5.26887e-01_jprb, 5.26482e-01_jprb, 5.26077e-01_jprb/)
      kao_mo3( 3, :, 5) = (/ &
     & 5.39814e-01_jprb, 5.39234e-01_jprb, 5.38655e-01_jprb, 5.38077e-01_jprb, 5.37499e-01_jprb, &
     & 5.36922e-01_jprb, 5.36345e-01_jprb, 5.35769e-01_jprb, 5.35194e-01_jprb, 5.34620e-01_jprb, &
     & 5.34045e-01_jprb, 5.33472e-01_jprb, 5.32899e-01_jprb, 5.32327e-01_jprb, 5.31756e-01_jprb, &
     & 5.31185e-01_jprb, 5.30614e-01_jprb, 5.30045e-01_jprb, 5.29475e-01_jprb/)
      kao_mo3( 4, :, 5) = (/ &
     & 5.39054e-01_jprb, 5.38348e-01_jprb, 5.37643e-01_jprb, 5.36938e-01_jprb, 5.36235e-01_jprb, &
     & 5.35532e-01_jprb, 5.34831e-01_jprb, 5.34130e-01_jprb, 5.33431e-01_jprb, 5.32732e-01_jprb, &
     & 5.32034e-01_jprb, 5.31337e-01_jprb, 5.30641e-01_jprb, 5.29946e-01_jprb, 5.29252e-01_jprb, &
     & 5.28559e-01_jprb, 5.27866e-01_jprb, 5.27175e-01_jprb, 5.26484e-01_jprb/)
      kao_mo3( 5, :, 5) = (/ &
     & 5.29240e-01_jprb, 5.28475e-01_jprb, 5.27711e-01_jprb, 5.26949e-01_jprb, 5.26187e-01_jprb, &
     & 5.25427e-01_jprb, 5.24668e-01_jprb, 5.23909e-01_jprb, 5.23152e-01_jprb, 5.22396e-01_jprb, &
     & 5.21641e-01_jprb, 5.20888e-01_jprb, 5.20135e-01_jprb, 5.19383e-01_jprb, 5.18633e-01_jprb, &
     & 5.17883e-01_jprb, 5.17135e-01_jprb, 5.16388e-01_jprb, 5.15642e-01_jprb/)
      kao_mo3( 6, :, 5) = (/ &
     & 5.21746e-01_jprb, 5.20815e-01_jprb, 5.19886e-01_jprb, 5.18958e-01_jprb, 5.18032e-01_jprb, &
     & 5.17107e-01_jprb, 5.16184e-01_jprb, 5.15263e-01_jprb, 5.14343e-01_jprb, 5.13425e-01_jprb, &
     & 5.12509e-01_jprb, 5.11594e-01_jprb, 5.10681e-01_jprb, 5.09770e-01_jprb, 5.08860e-01_jprb, &
     & 5.07952e-01_jprb, 5.07045e-01_jprb, 5.06140e-01_jprb, 5.05237e-01_jprb/)
      kao_mo3( 7, :, 5) = (/ &
     & 5.26752e-01_jprb, 5.25550e-01_jprb, 5.24352e-01_jprb, 5.23156e-01_jprb, 5.21963e-01_jprb, &
     & 5.20772e-01_jprb, 5.19584e-01_jprb, 5.18399e-01_jprb, 5.17217e-01_jprb, 5.16038e-01_jprb, &
     & 5.14861e-01_jprb, 5.13686e-01_jprb, 5.12515e-01_jprb, 5.11346e-01_jprb, 5.10180e-01_jprb, &
     & 5.09016e-01_jprb, 5.07855e-01_jprb, 5.06697e-01_jprb, 5.05541e-01_jprb/)
      kao_mo3( 8, :, 5) = (/ &
     & 5.23581e-01_jprb, 5.22513e-01_jprb, 5.21446e-01_jprb, 5.20382e-01_jprb, 5.19320e-01_jprb, &
     & 5.18260e-01_jprb, 5.17203e-01_jprb, 5.16147e-01_jprb, 5.15094e-01_jprb, 5.14042e-01_jprb, &
     & 5.12993e-01_jprb, 5.11946e-01_jprb, 5.10901e-01_jprb, 5.09859e-01_jprb, 5.08818e-01_jprb, &
     & 5.07780e-01_jprb, 5.06743e-01_jprb, 5.05709e-01_jprb, 5.04677e-01_jprb/)
      kao_mo3( 9, :, 5) = (/ &
     & 3.80393e-01_jprb, 3.80680e-01_jprb, 3.80967e-01_jprb, 3.81254e-01_jprb, 3.81542e-01_jprb, &
     & 3.81829e-01_jprb, 3.82117e-01_jprb, 3.82405e-01_jprb, 3.82693e-01_jprb, 3.82982e-01_jprb, &
     & 3.83271e-01_jprb, 3.83559e-01_jprb, 3.83849e-01_jprb, 3.84138e-01_jprb, 3.84428e-01_jprb, &
     & 3.84717e-01_jprb, 3.85007e-01_jprb, 3.85298e-01_jprb, 3.85588e-01_jprb/)
      kao_mo3( 1, :, 6) = (/ &
     & 6.14818e-01_jprb, 6.10664e-01_jprb, 6.06539e-01_jprb, 6.02441e-01_jprb, 5.98372e-01_jprb, &
     & 5.94330e-01_jprb, 5.90315e-01_jprb, 5.86327e-01_jprb, 5.82366e-01_jprb, 5.78432e-01_jprb, &
     & 5.74524e-01_jprb, 5.70643e-01_jprb, 5.66788e-01_jprb, 5.62959e-01_jprb, 5.59156e-01_jprb, &
     & 5.55379e-01_jprb, 5.51627e-01_jprb, 5.47901e-01_jprb, 5.44199e-01_jprb/)
      kao_mo3( 2, :, 6) = (/ &
     & 6.10199e-01_jprb, 6.06143e-01_jprb, 6.02114e-01_jprb, 5.98112e-01_jprb, 5.94136e-01_jprb, &
     & 5.90187e-01_jprb, 5.86264e-01_jprb, 5.82367e-01_jprb, 5.78496e-01_jprb, 5.74651e-01_jprb, &
     & 5.70831e-01_jprb, 5.67037e-01_jprb, 5.63268e-01_jprb, 5.59524e-01_jprb, 5.55805e-01_jprb, &
     & 5.52110e-01_jprb, 5.48440e-01_jprb, 5.44795e-01_jprb, 5.41174e-01_jprb/)
      kao_mo3( 3, :, 6) = (/ &
     & 6.02949e-01_jprb, 5.99057e-01_jprb, 5.95190e-01_jprb, 5.91348e-01_jprb, 5.87531e-01_jprb, &
     & 5.83739e-01_jprb, 5.79971e-01_jprb, 5.76227e-01_jprb, 5.72508e-01_jprb, 5.68812e-01_jprb, &
     & 5.65140e-01_jprb, 5.61493e-01_jprb, 5.57868e-01_jprb, 5.54267e-01_jprb, 5.50690e-01_jprb, &
     & 5.47135e-01_jprb, 5.43603e-01_jprb, 5.40094e-01_jprb, 5.36608e-01_jprb/)
      kao_mo3( 4, :, 6) = (/ &
     & 6.05047e-01_jprb, 6.01155e-01_jprb, 5.97289e-01_jprb, 5.93448e-01_jprb, 5.89631e-01_jprb, &
     & 5.85838e-01_jprb, 5.82071e-01_jprb, 5.78327e-01_jprb, 5.74607e-01_jprb, 5.70912e-01_jprb, &
     & 5.67240e-01_jprb, 5.63592e-01_jprb, 5.59967e-01_jprb, 5.56365e-01_jprb, 5.52787e-01_jprb, &
     & 5.49232e-01_jprb, 5.45699e-01_jprb, 5.42190e-01_jprb, 5.38703e-01_jprb/)
      kao_mo3( 5, :, 6) = (/ &
     & 6.03593e-01_jprb, 5.99867e-01_jprb, 5.96164e-01_jprb, 5.92483e-01_jprb, 5.88825e-01_jprb, &
     & 5.85190e-01_jprb, 5.81577e-01_jprb, 5.77987e-01_jprb, 5.74419e-01_jprb, 5.70872e-01_jprb, &
     & 5.67348e-01_jprb, 5.63846e-01_jprb, 5.60365e-01_jprb, 5.56905e-01_jprb, 5.53467e-01_jprb, &
     & 5.50050e-01_jprb, 5.46654e-01_jprb, 5.43279e-01_jprb, 5.39926e-01_jprb/)
      kao_mo3( 6, :, 6) = (/ &
     & 6.03940e-01_jprb, 6.00224e-01_jprb, 5.96531e-01_jprb, 5.92861e-01_jprb, 5.89213e-01_jprb, &
     & 5.85588e-01_jprb, 5.81985e-01_jprb, 5.78404e-01_jprb, 5.74845e-01_jprb, 5.71308e-01_jprb, &
     & 5.67793e-01_jprb, 5.64299e-01_jprb, 5.60827e-01_jprb, 5.57377e-01_jprb, 5.53947e-01_jprb, &
     & 5.50539e-01_jprb, 5.47151e-01_jprb, 5.43785e-01_jprb, 5.40439e-01_jprb/)
      kao_mo3( 7, :, 6) = (/ &
     & 6.06242e-01_jprb, 6.02257e-01_jprb, 5.98299e-01_jprb, 5.94367e-01_jprb, 5.90461e-01_jprb, &
     & 5.86580e-01_jprb, 5.82725e-01_jprb, 5.78895e-01_jprb, 5.75090e-01_jprb, 5.71311e-01_jprb, &
     & 5.67556e-01_jprb, 5.63826e-01_jprb, 5.60120e-01_jprb, 5.56439e-01_jprb, 5.52782e-01_jprb, &
     & 5.49149e-01_jprb, 5.45540e-01_jprb, 5.41954e-01_jprb, 5.38393e-01_jprb/)
      kao_mo3( 8, :, 6) = (/ &
     & 6.11929e-01_jprb, 6.07173e-01_jprb, 6.02454e-01_jprb, 5.97773e-01_jprb, 5.93127e-01_jprb, &
     & 5.88518e-01_jprb, 5.83944e-01_jprb, 5.79406e-01_jprb, 5.74903e-01_jprb, 5.70436e-01_jprb, &
     & 5.66002e-01_jprb, 5.61604e-01_jprb, 5.57239e-01_jprb, 5.52909e-01_jprb, 5.48612e-01_jprb, &
     & 5.44349e-01_jprb, 5.40118e-01_jprb, 5.35921e-01_jprb, 5.31756e-01_jprb/)
      kao_mo3( 9, :, 6) = (/ &
     & 6.21189e-01_jprb, 6.17338e-01_jprb, 6.13511e-01_jprb, 6.09707e-01_jprb, 6.05927e-01_jprb, &
     & 6.02170e-01_jprb, 5.98437e-01_jprb, 5.94726e-01_jprb, 5.91039e-01_jprb, 5.87375e-01_jprb, &
     & 5.83733e-01_jprb, 5.80114e-01_jprb, 5.76517e-01_jprb, 5.72943e-01_jprb, 5.69390e-01_jprb, &
     & 5.65860e-01_jprb, 5.62352e-01_jprb, 5.58865e-01_jprb, 5.55400e-01_jprb/)
      kao_mo3( 1, :, 7) = (/ &
     & 7.41310e-01_jprb, 7.30108e-01_jprb, 7.19075e-01_jprb, 7.08209e-01_jprb, 6.97507e-01_jprb, &
     & 6.86967e-01_jprb, 6.76586e-01_jprb, 6.66362e-01_jprb, 6.56292e-01_jprb, 6.46374e-01_jprb, &
     & 6.36607e-01_jprb, 6.26987e-01_jprb, 6.17512e-01_jprb, 6.08181e-01_jprb, 5.98990e-01_jprb, &
     & 5.89939e-01_jprb, 5.81024e-01_jprb, 5.72244e-01_jprb, 5.63597e-01_jprb/)
      kao_mo3( 2, :, 7) = (/ &
     & 7.38780e-01_jprb, 7.27631e-01_jprb, 7.16651e-01_jprb, 7.05836e-01_jprb, 6.95185e-01_jprb, &
     & 6.84695e-01_jprb, 6.74362e-01_jprb, 6.64186e-01_jprb, 6.54163e-01_jprb, 6.44292e-01_jprb, &
     & 6.34569e-01_jprb, 6.24993e-01_jprb, 6.15562e-01_jprb, 6.06273e-01_jprb, 5.97124e-01_jprb, &
     & 5.88113e-01_jprb, 5.79238e-01_jprb, 5.70498e-01_jprb, 5.61889e-01_jprb/)
      kao_mo3( 3, :, 7) = (/ &
     & 7.33846e-01_jprb, 7.22799e-01_jprb, 7.11919e-01_jprb, 7.01203e-01_jprb, 6.90648e-01_jprb, &
     & 6.80252e-01_jprb, 6.70012e-01_jprb, 6.59927e-01_jprb, 6.49993e-01_jprb, 6.40209e-01_jprb, &
     & 6.30572e-01_jprb, 6.21080e-01_jprb, 6.11731e-01_jprb, 6.02523e-01_jprb, 5.93453e-01_jprb, &
     & 5.84520e-01_jprb, 5.75721e-01_jprb, 5.67055e-01_jprb, 5.58519e-01_jprb/)
      kao_mo3( 4, :, 7) = (/ &
     & 7.21218e-01_jprb, 7.10492e-01_jprb, 6.99926e-01_jprb, 6.89517e-01_jprb, 6.79262e-01_jprb, &
     & 6.69160e-01_jprb, 6.59209e-01_jprb, 6.49405e-01_jprb, 6.39747e-01_jprb, 6.30233e-01_jprb, &
     & 6.20860e-01_jprb, 6.11627e-01_jprb, 6.02531e-01_jprb, 5.93570e-01_jprb, 5.84743e-01_jprb, &
     & 5.76047e-01_jprb, 5.67480e-01_jprb, 5.59040e-01_jprb, 5.50726e-01_jprb/)
      kao_mo3( 5, :, 7) = (/ &
     & 7.10588e-01_jprb, 7.00014e-01_jprb, 6.89596e-01_jprb, 6.79334e-01_jprb, 6.69225e-01_jprb, &
     & 6.59266e-01_jprb, 6.49455e-01_jprb, 6.39790e-01_jprb, 6.30269e-01_jprb, 6.20889e-01_jprb, &
     & 6.11650e-01_jprb, 6.02547e-01_jprb, 5.93581e-01_jprb, 5.84747e-01_jprb, 5.76045e-01_jprb, &
     & 5.67473e-01_jprb, 5.59028e-01_jprb, 5.50709e-01_jprb, 5.42513e-01_jprb/)
      kao_mo3( 6, :, 7) = (/ &
     & 6.98166e-01_jprb, 6.87706e-01_jprb, 6.77402e-01_jprb, 6.67253e-01_jprb, 6.57256e-01_jprb, &
     & 6.47408e-01_jprb, 6.37708e-01_jprb, 6.28154e-01_jprb, 6.18742e-01_jprb, 6.09472e-01_jprb, &
     & 6.00340e-01_jprb, 5.91346e-01_jprb, 5.82486e-01_jprb, 5.73758e-01_jprb, 5.65162e-01_jprb, &
     & 5.56694e-01_jprb, 5.48353e-01_jprb, 5.40138e-01_jprb, 5.32045e-01_jprb/)
      kao_mo3( 7, :, 7) = (/ &
     & 6.76974e-01_jprb, 6.67034e-01_jprb, 6.57240e-01_jprb, 6.47590e-01_jprb, 6.38081e-01_jprb, &
     & 6.28712e-01_jprb, 6.19481e-01_jprb, 6.10385e-01_jprb, 6.01422e-01_jprb, 5.92592e-01_jprb, &
     & 5.83891e-01_jprb, 5.75317e-01_jprb, 5.66870e-01_jprb, 5.58547e-01_jprb, 5.50345e-01_jprb, &
     & 5.42265e-01_jprb, 5.34303e-01_jprb, 5.26457e-01_jprb, 5.18727e-01_jprb/)
      kao_mo3( 8, :, 7) = (/ &
     & 6.30061e-01_jprb, 6.21017e-01_jprb, 6.12102e-01_jprb, 6.03316e-01_jprb, 5.94656e-01_jprb, &
     & 5.86120e-01_jprb, 5.77706e-01_jprb, 5.69414e-01_jprb, 5.61240e-01_jprb, 5.53184e-01_jprb, &
     & 5.45243e-01_jprb, 5.37416e-01_jprb, 5.29702e-01_jprb, 5.22098e-01_jprb, 5.14604e-01_jprb, &
     & 5.07217e-01_jprb, 4.99936e-01_jprb, 4.92760e-01_jprb, 4.85687e-01_jprb/)
      kao_mo3( 9, :, 7) = (/ &
     & 8.97633e-01_jprb, 8.87307e-01_jprb, 8.77100e-01_jprb, 8.67010e-01_jprb, 8.57036e-01_jprb, &
     & 8.47176e-01_jprb, 8.37431e-01_jprb, 8.27797e-01_jprb, 8.18274e-01_jprb, 8.08861e-01_jprb, &
     & 7.99555e-01_jprb, 7.90357e-01_jprb, 7.81265e-01_jprb, 7.72278e-01_jprb, 7.63393e-01_jprb, &
     & 7.54611e-01_jprb, 7.45930e-01_jprb, 7.37349e-01_jprb, 7.28867e-01_jprb/)
      kao_mo3( 1, :, 8) = (/ &
     & 4.87356e-01_jprb, 4.80743e-01_jprb, 4.74220e-01_jprb, 4.67785e-01_jprb, 4.61437e-01_jprb, &
     & 4.55176e-01_jprb, 4.49000e-01_jprb, 4.42907e-01_jprb, 4.36897e-01_jprb, 4.30969e-01_jprb, &
     & 4.25121e-01_jprb, 4.19353e-01_jprb, 4.13663e-01_jprb, 4.08049e-01_jprb, 4.02513e-01_jprb, &
     & 3.97051e-01_jprb, 3.91663e-01_jprb, 3.86349e-01_jprb, 3.81106e-01_jprb/)
      kao_mo3( 2, :, 8) = (/ &
     & 4.86776e-01_jprb, 4.80157e-01_jprb, 4.73627e-01_jprb, 4.67187e-01_jprb, 4.60834e-01_jprb, &
     & 4.54567e-01_jprb, 4.48386e-01_jprb, 4.42289e-01_jprb, 4.36274e-01_jprb, 4.30342e-01_jprb, &
     & 4.24490e-01_jprb, 4.18718e-01_jprb, 4.13024e-01_jprb, 4.07407e-01_jprb, 4.01867e-01_jprb, &
     & 3.96403e-01_jprb, 3.91012e-01_jprb, 3.85695e-01_jprb, 3.80450e-01_jprb/)
      kao_mo3( 3, :, 8) = (/ &
     & 4.86111e-01_jprb, 4.79496e-01_jprb, 4.72972e-01_jprb, 4.66536e-01_jprb, 4.60188e-01_jprb, &
     & 4.53926e-01_jprb, 4.47750e-01_jprb, 4.41657e-01_jprb, 4.35648e-01_jprb, 4.29720e-01_jprb, &
     & 4.23873e-01_jprb, 4.18105e-01_jprb, 4.12416e-01_jprb, 4.06804e-01_jprb, 4.01269e-01_jprb, &
     & 3.95809e-01_jprb, 3.90423e-01_jprb, 3.85111e-01_jprb, 3.79871e-01_jprb/)
      kao_mo3( 4, :, 8) = (/ &
     & 4.85501e-01_jprb, 4.78880e-01_jprb, 4.72350e-01_jprb, 4.65908e-01_jprb, 4.59554e-01_jprb, &
     & 4.53288e-01_jprb, 4.47106e-01_jprb, 4.41009e-01_jprb, 4.34995e-01_jprb, 4.29063e-01_jprb, &
     & 4.23211e-01_jprb, 4.17440e-01_jprb, 4.11747e-01_jprb, 4.06132e-01_jprb, 4.00594e-01_jprb, &
     & 3.95131e-01_jprb, 3.89743e-01_jprb, 3.84428e-01_jprb, 3.79185e-01_jprb/)
      kao_mo3( 5, :, 8) = (/ &
     & 4.83679e-01_jprb, 4.77140e-01_jprb, 4.70691e-01_jprb, 4.64328e-01_jprb, 4.58051e-01_jprb, &
     & 4.51859e-01_jprb, 4.45751e-01_jprb, 4.39726e-01_jprb, 4.33781e-01_jprb, 4.27918e-01_jprb, &
     & 4.22133e-01_jprb, 4.16427e-01_jprb, 4.10798e-01_jprb, 4.05245e-01_jprb, 3.99767e-01_jprb, &
     & 3.94363e-01_jprb, 3.89032e-01_jprb, 3.83773e-01_jprb, 3.78585e-01_jprb/)
      kao_mo3( 6, :, 8) = (/ &
     & 4.72120e-01_jprb, 4.65834e-01_jprb, 4.59630e-01_jprb, 4.53510e-01_jprb, 4.47471e-01_jprb, &
     & 4.41513e-01_jprb, 4.35633e-01_jprb, 4.29833e-01_jprb, 4.24109e-01_jprb, 4.18461e-01_jprb, &
     & 4.12889e-01_jprb, 4.07391e-01_jprb, 4.01966e-01_jprb, 3.96614e-01_jprb, 3.91332e-01_jprb, &
     & 3.86122e-01_jprb, 3.80980e-01_jprb, 3.75907e-01_jprb, 3.70901e-01_jprb/)
      kao_mo3( 7, :, 8) = (/ &
     & 4.58683e-01_jprb, 4.52758e-01_jprb, 4.46909e-01_jprb, 4.41135e-01_jprb, 4.35437e-01_jprb, &
     & 4.29812e-01_jprb, 4.24259e-01_jprb, 4.18779e-01_jprb, 4.13369e-01_jprb, 4.08029e-01_jprb, &
     & 4.02758e-01_jprb, 3.97555e-01_jprb, 3.92419e-01_jprb, 3.87350e-01_jprb, 3.82346e-01_jprb, &
     & 3.77406e-01_jprb, 3.72531e-01_jprb, 3.67719e-01_jprb, 3.62968e-01_jprb/)
      kao_mo3( 8, :, 8) = (/ &
     & 4.56091e-01_jprb, 4.50481e-01_jprb, 4.44940e-01_jprb, 4.39467e-01_jprb, 4.34062e-01_jprb, &
     & 4.28722e-01_jprb, 4.23449e-01_jprb, 4.18240e-01_jprb, 4.13096e-01_jprb, 4.08015e-01_jprb, &
     & 4.02996e-01_jprb, 3.98039e-01_jprb, 3.93143e-01_jprb, 3.88307e-01_jprb, 3.83531e-01_jprb, &
     & 3.78813e-01_jprb, 3.74154e-01_jprb, 3.69552e-01_jprb, 3.65006e-01_jprb/)
      kao_mo3( 9, :, 8) = (/ &
     & 9.11213e-01_jprb, 9.03270e-01_jprb, 8.95396e-01_jprb, 8.87591e-01_jprb, 8.79855e-01_jprb, &
     & 8.72185e-01_jprb, 8.64583e-01_jprb, 8.57046e-01_jprb, 8.49576e-01_jprb, 8.42170e-01_jprb, &
     & 8.34829e-01_jprb, 8.27552e-01_jprb, 8.20339e-01_jprb, 8.13188e-01_jprb, 8.06100e-01_jprb, &
     & 7.99073e-01_jprb, 7.92108e-01_jprb, 7.85204e-01_jprb, 7.78359e-01_jprb/)
      kao_mo3( 1, :, 9) = (/ &
     & 5.56194e-01_jprb, 5.48595e-01_jprb, 5.41100e-01_jprb, 5.33707e-01_jprb, 5.26415e-01_jprb, &
     & 5.19223e-01_jprb, 5.12129e-01_jprb, 5.05132e-01_jprb, 4.98231e-01_jprb, 4.91424e-01_jprb, &
     & 4.84710e-01_jprb, 4.78087e-01_jprb, 4.71556e-01_jprb, 4.65113e-01_jprb, 4.58758e-01_jprb, &
     & 4.52491e-01_jprb, 4.46309e-01_jprb, 4.40211e-01_jprb, 4.34197e-01_jprb/)
      kao_mo3( 2, :, 9) = (/ &
     & 5.56174e-01_jprb, 5.48575e-01_jprb, 5.41079e-01_jprb, 5.33687e-01_jprb, 5.26395e-01_jprb, &
     & 5.19203e-01_jprb, 5.12109e-01_jprb, 5.05112e-01_jprb, 4.98211e-01_jprb, 4.91404e-01_jprb, &
     & 4.84690e-01_jprb, 4.78068e-01_jprb, 4.71536e-01_jprb, 4.65093e-01_jprb, 4.58739e-01_jprb, &
     & 4.52471e-01_jprb, 4.46289e-01_jprb, 4.40191e-01_jprb, 4.34177e-01_jprb/)
      kao_mo3( 3, :, 9) = (/ &
     & 5.55996e-01_jprb, 5.48403e-01_jprb, 5.40913e-01_jprb, 5.33526e-01_jprb, 5.26239e-01_jprb, &
     & 5.19052e-01_jprb, 5.11963e-01_jprb, 5.04971e-01_jprb, 4.98074e-01_jprb, 4.91272e-01_jprb, &
     & 4.84562e-01_jprb, 4.77944e-01_jprb, 4.71417e-01_jprb, 4.64978e-01_jprb, 4.58628e-01_jprb, &
     & 4.52364e-01_jprb, 4.46186e-01_jprb, 4.40092e-01_jprb, 4.34081e-01_jprb/)
      kao_mo3( 4, :, 9) = (/ &
     & 5.55859e-01_jprb, 5.48271e-01_jprb, 5.40786e-01_jprb, 5.33404e-01_jprb, 5.26123e-01_jprb, &
     & 5.18941e-01_jprb, 5.11856e-01_jprb, 5.04869e-01_jprb, 4.97977e-01_jprb, 4.91179e-01_jprb, &
     & 4.84474e-01_jprb, 4.77861e-01_jprb, 4.71337e-01_jprb, 4.64903e-01_jprb, 4.58557e-01_jprb, &
     & 4.52297e-01_jprb, 4.46123e-01_jprb, 4.40033e-01_jprb, 4.34026e-01_jprb/)
      kao_mo3( 5, :, 9) = (/ &
     & 5.54550e-01_jprb, 5.46921e-01_jprb, 5.39397e-01_jprb, 5.31976e-01_jprb, 5.24657e-01_jprb, &
     & 5.17439e-01_jprb, 5.10320e-01_jprb, 5.03300e-01_jprb, 4.96376e-01_jprb, 4.89547e-01_jprb, &
     & 4.82812e-01_jprb, 4.76170e-01_jprb, 4.69619e-01_jprb, 4.63158e-01_jprb, 4.56786e-01_jprb, &
     & 4.50502e-01_jprb, 4.44304e-01_jprb, 4.38192e-01_jprb, 4.32163e-01_jprb/)
      kao_mo3( 6, :, 9) = (/ &
     & 5.53514e-01_jprb, 5.45883e-01_jprb, 5.38358e-01_jprb, 5.30937e-01_jprb, 5.23618e-01_jprb, &
     & 5.16399e-01_jprb, 5.09280e-01_jprb, 5.02260e-01_jprb, 4.95336e-01_jprb, 4.88507e-01_jprb, &
     & 4.81773e-01_jprb, 4.75132e-01_jprb, 4.68582e-01_jprb, 4.62122e-01_jprb, 4.55752e-01_jprb, &
     & 4.49469e-01_jprb, 4.43273e-01_jprb, 4.37162e-01_jprb, 4.31136e-01_jprb/)
      kao_mo3( 7, :, 9) = (/ &
     & 5.49865e-01_jprb, 5.42303e-01_jprb, 5.34846e-01_jprb, 5.27491e-01_jprb, 5.20237e-01_jprb, &
     & 5.13084e-01_jprb, 5.06028e-01_jprb, 4.99070e-01_jprb, 4.92207e-01_jprb, 4.85438e-01_jprb, &
     & 4.78763e-01_jprb, 4.72179e-01_jprb, 4.65686e-01_jprb, 4.59282e-01_jprb, 4.52967e-01_jprb, &
     & 4.46738e-01_jprb, 4.40595e-01_jprb, 4.34536e-01_jprb, 4.28561e-01_jprb/)
      kao_mo3( 8, :, 9) = (/ &
     & 5.25435e-01_jprb, 5.18437e-01_jprb, 5.11533e-01_jprb, 5.04721e-01_jprb, 4.97999e-01_jprb, &
     & 4.91367e-01_jprb, 4.84823e-01_jprb, 4.78366e-01_jprb, 4.71996e-01_jprb, 4.65710e-01_jprb, &
     & 4.59508e-01_jprb, 4.53388e-01_jprb, 4.47350e-01_jprb, 4.41393e-01_jprb, 4.35515e-01_jprb, &
     & 4.29715e-01_jprb, 4.23992e-01_jprb, 4.18345e-01_jprb, 4.12774e-01_jprb/)
      kao_mo3( 9, :, 9) = (/ &
     & 3.48228e-01_jprb, 3.45949e-01_jprb, 3.43686e-01_jprb, 3.41437e-01_jprb, 3.39203e-01_jprb, &
     & 3.36983e-01_jprb, 3.34778e-01_jprb, 3.32588e-01_jprb, 3.30412e-01_jprb, 3.28250e-01_jprb, &
     & 3.26102e-01_jprb, 3.23968e-01_jprb, 3.21848e-01_jprb, 3.19742e-01_jprb, 3.17650e-01_jprb, &
     & 3.15572e-01_jprb, 3.13507e-01_jprb, 3.11456e-01_jprb, 3.09418e-01_jprb/)
      kao_mo3( 1, :,10) = (/ &
     & 8.34107e-01_jprb, 8.27276e-01_jprb, 8.20501e-01_jprb, 8.13781e-01_jprb, 8.07117e-01_jprb, &
     & 8.00507e-01_jprb, 7.93951e-01_jprb, 7.87449e-01_jprb, 7.81000e-01_jprb, 7.74604e-01_jprb, &
     & 7.68260e-01_jprb, 7.61968e-01_jprb, 7.55728e-01_jprb, 7.49539e-01_jprb, 7.43400e-01_jprb, &
     & 7.37312e-01_jprb, 7.31274e-01_jprb, 7.25285e-01_jprb, 7.19345e-01_jprb/)
      kao_mo3( 2, :,10) = (/ &
     & 8.32838e-01_jprb, 8.26022e-01_jprb, 8.19263e-01_jprb, 8.12558e-01_jprb, 8.05908e-01_jprb, &
     & 7.99313e-01_jprb, 7.92772e-01_jprb, 7.86284e-01_jprb, 7.79849e-01_jprb, 7.73467e-01_jprb, &
     & 7.67137e-01_jprb, 7.60859e-01_jprb, 7.54633e-01_jprb, 7.48457e-01_jprb, 7.42332e-01_jprb, &
     & 7.36257e-01_jprb, 7.30232e-01_jprb, 7.24256e-01_jprb, 7.18329e-01_jprb/)
      kao_mo3( 3, :,10) = (/ &
     & 8.31167e-01_jprb, 8.24361e-01_jprb, 8.17611e-01_jprb, 8.10916e-01_jprb, 8.04276e-01_jprb, &
     & 7.97691e-01_jprb, 7.91159e-01_jprb, 7.84681e-01_jprb, 7.78256e-01_jprb, 7.71883e-01_jprb, &
     & 7.65563e-01_jprb, 7.59294e-01_jprb, 7.53077e-01_jprb, 7.46910e-01_jprb, 7.40795e-01_jprb, &
     & 7.34729e-01_jprb, 7.28713e-01_jprb, 7.22746e-01_jprb, 7.16828e-01_jprb/)
      kao_mo3( 4, :,10) = (/ &
     & 8.29026e-01_jprb, 8.22246e-01_jprb, 8.15521e-01_jprb, 8.08851e-01_jprb, 8.02236e-01_jprb, &
     & 7.95675e-01_jprb, 7.89167e-01_jprb, 7.82713e-01_jprb, 7.76312e-01_jprb, 7.69962e-01_jprb, &
     & 7.63665e-01_jprb, 7.57419e-01_jprb, 7.51225e-01_jprb, 7.45081e-01_jprb, 7.38987e-01_jprb, &
     & 7.32943e-01_jprb, 7.26949e-01_jprb, 7.21003e-01_jprb, 7.15107e-01_jprb/)
      kao_mo3( 5, :,10) = (/ &
     & 8.26226e-01_jprb, 8.19471e-01_jprb, 8.12771e-01_jprb, 8.06126e-01_jprb, 7.99536e-01_jprb, &
     & 7.92999e-01_jprb, 7.86515e-01_jprb, 7.80085e-01_jprb, 7.73707e-01_jprb, 7.67382e-01_jprb, &
     & 7.61108e-01_jprb, 7.54885e-01_jprb, 7.48714e-01_jprb, 7.42592e-01_jprb, 7.36521e-01_jprb, &
     & 7.30500e-01_jprb, 7.24527e-01_jprb, 7.18604e-01_jprb, 7.12729e-01_jprb/)
      kao_mo3( 6, :,10) = (/ &
     & 8.33246e-01_jprb, 8.26510e-01_jprb, 8.19828e-01_jprb, 8.13200e-01_jprb, 8.06626e-01_jprb, &
     & 8.00105e-01_jprb, 7.93637e-01_jprb, 7.87221e-01_jprb, 7.80856e-01_jprb, 7.74544e-01_jprb, &
     & 7.68282e-01_jprb, 7.62071e-01_jprb, 7.55910e-01_jprb, 7.49799e-01_jprb, 7.43737e-01_jprb, &
     & 7.37725e-01_jprb, 7.31760e-01_jprb, 7.25845e-01_jprb, 7.19977e-01_jprb/)
      kao_mo3( 7, :,10) = (/ &
     & 8.45693e-01_jprb, 8.38967e-01_jprb, 8.32295e-01_jprb, 8.25675e-01_jprb, 8.19108e-01_jprb, &
     & 8.12594e-01_jprb, 8.06131e-01_jprb, 7.99719e-01_jprb, 7.93359e-01_jprb, 7.87049e-01_jprb, &
     & 7.80789e-01_jprb, 7.74579e-01_jprb, 7.68419e-01_jprb, 7.62307e-01_jprb, 7.56244e-01_jprb, &
     & 7.50230e-01_jprb, 7.44263e-01_jprb, 7.38343e-01_jprb, 7.32471e-01_jprb/)
      kao_mo3( 8, :,10) = (/ &
     & 8.32139e-01_jprb, 8.25565e-01_jprb, 8.19044e-01_jprb, 8.12574e-01_jprb, 8.06156e-01_jprb, &
     & 7.99788e-01_jprb, 7.93470e-01_jprb, 7.87202e-01_jprb, 7.80984e-01_jprb, 7.74815e-01_jprb, &
     & 7.68694e-01_jprb, 7.62622e-01_jprb, 7.56598e-01_jprb, 7.50622e-01_jprb, 7.44692e-01_jprb, &
     & 7.38810e-01_jprb, 7.32974e-01_jprb, 7.27184e-01_jprb, 7.21440e-01_jprb/)
      kao_mo3( 9, :,10) = (/ &
     & 2.34258e-01_jprb, 2.35247e-01_jprb, 2.36239e-01_jprb, 2.37236e-01_jprb, 2.38237e-01_jprb, &
     & 2.39242e-01_jprb, 2.40252e-01_jprb, 2.41265e-01_jprb, 2.42283e-01_jprb, 2.43306e-01_jprb, &
     & 2.44332e-01_jprb, 2.45363e-01_jprb, 2.46398e-01_jprb, 2.47438e-01_jprb, 2.48482e-01_jprb, &
     & 2.49531e-01_jprb, 2.50583e-01_jprb, 2.51641e-01_jprb, 2.52702e-01_jprb/)
      kao_mo3( 1, :,11) = (/ &
     & 8.31308e-01_jprb, 8.22153e-01_jprb, 8.13098e-01_jprb, 8.04143e-01_jprb, 7.95287e-01_jprb, &
     & 7.86528e-01_jprb, 7.77866e-01_jprb, 7.69299e-01_jprb, 7.60827e-01_jprb, 7.52448e-01_jprb, &
     & 7.44161e-01_jprb, 7.35965e-01_jprb, 7.27860e-01_jprb, 7.19844e-01_jprb, 7.11916e-01_jprb, &
     & 7.04075e-01_jprb, 6.96321e-01_jprb, 6.88652e-01_jprb, 6.81068e-01_jprb/)
      kao_mo3( 2, :,11) = (/ &
     & 8.31577e-01_jprb, 8.22400e-01_jprb, 8.13324e-01_jprb, 8.04349e-01_jprb, 7.95472e-01_jprb, &
     & 7.86693e-01_jprb, 7.78011e-01_jprb, 7.69425e-01_jprb, 7.60934e-01_jprb, 7.52537e-01_jprb, &
     & 7.44232e-01_jprb, 7.36019e-01_jprb, 7.27896e-01_jprb, 7.19863e-01_jprb, 7.11919e-01_jprb, &
     & 7.04062e-01_jprb, 6.96292e-01_jprb, 6.88608e-01_jprb, 6.81009e-01_jprb/)
      kao_mo3( 3, :,11) = (/ &
     & 8.31578e-01_jprb, 8.22422e-01_jprb, 8.13368e-01_jprb, 8.04413e-01_jprb, 7.95557e-01_jprb, &
     & 7.86798e-01_jprb, 7.78136e-01_jprb, 7.69569e-01_jprb, 7.61097e-01_jprb, 7.52717e-01_jprb, &
     & 7.44430e-01_jprb, 7.36235e-01_jprb, 7.28129e-01_jprb, 7.20113e-01_jprb, 7.12185e-01_jprb, &
     & 7.04344e-01_jprb, 6.96589e-01_jprb, 6.88920e-01_jprb, 6.81336e-01_jprb/)
      kao_mo3( 4, :,11) = (/ &
     & 8.31261e-01_jprb, 8.22111e-01_jprb, 8.13062e-01_jprb, 8.04112e-01_jprb, 7.95261e-01_jprb, &
     & 7.86507e-01_jprb, 7.77850e-01_jprb, 7.69288e-01_jprb, 7.60820e-01_jprb, 7.52445e-01_jprb, &
     & 7.44163e-01_jprb, 7.35971e-01_jprb, 7.27870e-01_jprb, 7.19858e-01_jprb, 7.11935e-01_jprb, &
     & 7.04098e-01_jprb, 6.96348e-01_jprb, 6.88683e-01_jprb, 6.81102e-01_jprb/)
      kao_mo3( 5, :,11) = (/ &
     & 8.31565e-01_jprb, 8.22404e-01_jprb, 8.13344e-01_jprb, 8.04384e-01_jprb, 7.95523e-01_jprb, &
     & 7.86760e-01_jprb, 7.78092e-01_jprb, 7.69521e-01_jprb, 7.61044e-01_jprb, 7.52660e-01_jprb, &
     & 7.44368e-01_jprb, 7.36168e-01_jprb, 7.28058e-01_jprb, 7.20038e-01_jprb, 7.12106e-01_jprb, &
     & 7.04261e-01_jprb, 6.96503e-01_jprb, 6.88830e-01_jprb, 6.81242e-01_jprb/)
      kao_mo3( 6, :,11) = (/ &
     & 8.17636e-01_jprb, 8.08497e-01_jprb, 7.99461e-01_jprb, 7.90525e-01_jprb, 7.81690e-01_jprb, &
     & 7.72953e-01_jprb, 7.64314e-01_jprb, 7.55771e-01_jprb, 7.47324e-01_jprb, 7.38971e-01_jprb, &
     & 7.30712e-01_jprb, 7.22545e-01_jprb, 7.14469e-01_jprb, 7.06483e-01_jprb, 6.98587e-01_jprb, &
     & 6.90779e-01_jprb, 6.83058e-01_jprb, 6.75424e-01_jprb, 6.67875e-01_jprb/)
      kao_mo3( 7, :,11) = (/ &
     & 7.95247e-01_jprb, 7.86140e-01_jprb, 7.77137e-01_jprb, 7.68238e-01_jprb, 7.59440e-01_jprb, &
     & 7.50743e-01_jprb, 7.42145e-01_jprb, 7.33646e-01_jprb, 7.25245e-01_jprb, 7.16939e-01_jprb, &
     & 7.08729e-01_jprb, 7.00612e-01_jprb, 6.92589e-01_jprb, 6.84658e-01_jprb, 6.76817e-01_jprb, &
     & 6.69066e-01_jprb, 6.61404e-01_jprb, 6.53830e-01_jprb, 6.46342e-01_jprb/)
      kao_mo3( 8, :,11) = (/ &
     & 7.63069e-01_jprb, 7.54006e-01_jprb, 7.45051e-01_jprb, 7.36202e-01_jprb, 7.27458e-01_jprb, &
     & 7.18818e-01_jprb, 7.10281e-01_jprb, 7.01845e-01_jprb, 6.93509e-01_jprb, 6.85272e-01_jprb, &
     & 6.77133e-01_jprb, 6.69091e-01_jprb, 6.61144e-01_jprb, 6.53292e-01_jprb, 6.45533e-01_jprb, &
     & 6.37866e-01_jprb, 6.30290e-01_jprb, 6.22804e-01_jprb, 6.15407e-01_jprb/)
      kao_mo3( 9, :,11) = (/ &
     & 2.03255e-01_jprb, 2.03004e-01_jprb, 2.02753e-01_jprb, 2.02502e-01_jprb, 2.02252e-01_jprb, &
     & 2.02001e-01_jprb, 2.01752e-01_jprb, 2.01502e-01_jprb, 2.01253e-01_jprb, 2.01004e-01_jprb, &
     & 2.00755e-01_jprb, 2.00507e-01_jprb, 2.00259e-01_jprb, 2.00011e-01_jprb, 1.99764e-01_jprb, &
     & 1.99517e-01_jprb, 1.99270e-01_jprb, 1.99024e-01_jprb, 1.98777e-01_jprb/)
      kao_mo3( 1, :,12) = (/ &
     & 4.13201e-01_jprb, 4.05258e-01_jprb, 3.97468e-01_jprb, 3.89828e-01_jprb, 3.82334e-01_jprb, &
     & 3.74985e-01_jprb, 3.67777e-01_jprb, 3.60707e-01_jprb, 3.53774e-01_jprb, 3.46973e-01_jprb, &
     & 3.40303e-01_jprb, 3.33762e-01_jprb, 3.27346e-01_jprb, 3.21054e-01_jprb, 3.14882e-01_jprb, &
     & 3.08829e-01_jprb, 3.02893e-01_jprb, 2.97071e-01_jprb, 2.91360e-01_jprb/)
      kao_mo3( 2, :,12) = (/ &
     & 4.12835e-01_jprb, 4.04897e-01_jprb, 3.97112e-01_jprb, 3.89477e-01_jprb, 3.81988e-01_jprb, &
     & 3.74644e-01_jprb, 3.67440e-01_jprb, 3.60376e-01_jprb, 3.53447e-01_jprb, 3.46651e-01_jprb, &
     & 3.39986e-01_jprb, 3.33449e-01_jprb, 3.27038e-01_jprb, 3.20750e-01_jprb, 3.14582e-01_jprb, &
     & 3.08534e-01_jprb, 3.02602e-01_jprb, 2.96784e-01_jprb, 2.91077e-01_jprb/)
      kao_mo3( 3, :,12) = (/ &
     & 4.13023e-01_jprb, 4.05079e-01_jprb, 3.97289e-01_jprb, 3.89648e-01_jprb, 3.82155e-01_jprb, &
     & 3.74805e-01_jprb, 3.67597e-01_jprb, 3.60527e-01_jprb, 3.53594e-01_jprb, 3.46793e-01_jprb, &
     & 3.40124e-01_jprb, 3.33583e-01_jprb, 3.27167e-01_jprb, 3.20875e-01_jprb, 3.14704e-01_jprb, &
     & 3.08652e-01_jprb, 3.02716e-01_jprb, 2.96894e-01_jprb, 2.91184e-01_jprb/)
      kao_mo3( 4, :,12) = (/ &
     & 4.13397e-01_jprb, 4.05437e-01_jprb, 3.97630e-01_jprb, 3.89973e-01_jprb, 3.82463e-01_jprb, &
     & 3.75099e-01_jprb, 3.67876e-01_jprb, 3.60792e-01_jprb, 3.53844e-01_jprb, 3.47031e-01_jprb, &
     & 3.40348e-01_jprb, 3.33794e-01_jprb, 3.27367e-01_jprb, 3.21063e-01_jprb, 3.14880e-01_jprb, &
     & 3.08817e-01_jprb, 3.02870e-01_jprb, 2.97038e-01_jprb, 2.91318e-01_jprb/)
      kao_mo3( 5, :,12) = (/ &
     & 4.13043e-01_jprb, 4.05106e-01_jprb, 3.97321e-01_jprb, 3.89686e-01_jprb, 3.82198e-01_jprb, &
     & 3.74854e-01_jprb, 3.67651e-01_jprb, 3.60586e-01_jprb, 3.53657e-01_jprb, 3.46861e-01_jprb, &
     & 3.40195e-01_jprb, 3.33658e-01_jprb, 3.27246e-01_jprb, 3.20958e-01_jprb, 3.14790e-01_jprb, &
     & 3.08741e-01_jprb, 3.02808e-01_jprb, 2.96990e-01_jprb, 2.91283e-01_jprb/)
      kao_mo3( 6, :,12) = (/ &
     & 4.13151e-01_jprb, 4.05202e-01_jprb, 3.97406e-01_jprb, 3.89760e-01_jprb, 3.82261e-01_jprb, &
     & 3.74906e-01_jprb, 3.67693e-01_jprb, 3.60619e-01_jprb, 3.53680e-01_jprb, 3.46876e-01_jprb, &
     & 3.40202e-01_jprb, 3.33656e-01_jprb, 3.27237e-01_jprb, 3.20941e-01_jprb, 3.14766e-01_jprb, &
     & 3.08710e-01_jprb, 3.02770e-01_jprb, 2.96945e-01_jprb, 2.91232e-01_jprb/)
      kao_mo3( 7, :,12) = (/ &
     & 4.13052e-01_jprb, 4.05109e-01_jprb, 3.97319e-01_jprb, 3.89678e-01_jprb, 3.82185e-01_jprb, &
     & 3.74835e-01_jprb, 3.67627e-01_jprb, 3.60557e-01_jprb, 3.53624e-01_jprb, 3.46823e-01_jprb, &
     & 3.40154e-01_jprb, 3.33612e-01_jprb, 3.27197e-01_jprb, 3.20905e-01_jprb, 3.14734e-01_jprb, &
     & 3.08681e-01_jprb, 3.02745e-01_jprb, 2.96923e-01_jprb, 2.91213e-01_jprb/)
      kao_mo3( 8, :,12) = (/ &
     & 4.13152e-01_jprb, 4.05209e-01_jprb, 3.97418e-01_jprb, 3.89778e-01_jprb, 3.82284e-01_jprb, &
     & 3.74935e-01_jprb, 3.67727e-01_jprb, 3.60657e-01_jprb, 3.53723e-01_jprb, 3.46923e-01_jprb, &
     & 3.40253e-01_jprb, 3.33712e-01_jprb, 3.27296e-01_jprb, 3.21004e-01_jprb, 3.14833e-01_jprb, &
     & 3.08780e-01_jprb, 3.02844e-01_jprb, 2.97021e-01_jprb, 2.91311e-01_jprb/)
      kao_mo3( 9, :,12) = (/ &
     & 1.31008e-01_jprb, 1.30607e-01_jprb, 1.30208e-01_jprb, 1.29810e-01_jprb, 1.29413e-01_jprb, &
     & 1.29017e-01_jprb, 1.28623e-01_jprb, 1.28229e-01_jprb, 1.27837e-01_jprb, 1.27446e-01_jprb, &
     & 1.27056e-01_jprb, 1.26668e-01_jprb, 1.26280e-01_jprb, 1.25894e-01_jprb, 1.25509e-01_jprb, &
     & 1.25125e-01_jprb, 1.24743e-01_jprb, 1.24361e-01_jprb, 1.23981e-01_jprb/)
      kao_mo3( 1, :,13) = (/ &
     & 4.66826e-01_jprb, 4.71437e-01_jprb, 4.76094e-01_jprb, 4.80798e-01_jprb, 4.85547e-01_jprb, &
     & 4.90344e-01_jprb, 4.95187e-01_jprb, 5.00079e-01_jprb, 5.05019e-01_jprb, 5.10008e-01_jprb, &
     & 5.15046e-01_jprb, 5.20134e-01_jprb, 5.25272e-01_jprb, 5.30461e-01_jprb, 5.35701e-01_jprb, &
     & 5.40993e-01_jprb, 5.46338e-01_jprb, 5.51735e-01_jprb, 5.57185e-01_jprb/)
      kao_mo3( 2, :,13) = (/ &
     & 4.66579e-01_jprb, 4.71199e-01_jprb, 4.75865e-01_jprb, 4.80577e-01_jprb, 4.85336e-01_jprb, &
     & 4.90141e-01_jprb, 4.94995e-01_jprb, 4.99896e-01_jprb, 5.04846e-01_jprb, 5.09845e-01_jprb, &
     & 5.14893e-01_jprb, 5.19992e-01_jprb, 5.25141e-01_jprb, 5.30340e-01_jprb, 5.35592e-01_jprb, &
     & 5.40895e-01_jprb, 5.46251e-01_jprb, 5.51660e-01_jprb, 5.57122e-01_jprb/)
      kao_mo3( 3, :,13) = (/ &
     & 4.66956e-01_jprb, 4.71567e-01_jprb, 4.76224e-01_jprb, 4.80927e-01_jprb, 4.85677e-01_jprb, &
     & 4.90474e-01_jprb, 4.95318e-01_jprb, 5.00209e-01_jprb, 5.05149e-01_jprb, 5.10138e-01_jprb, &
     & 5.15176e-01_jprb, 5.20264e-01_jprb, 5.25402e-01_jprb, 5.30591e-01_jprb, 5.35831e-01_jprb, &
     & 5.41123e-01_jprb, 5.46467e-01_jprb, 5.51864e-01_jprb, 5.57314e-01_jprb/)
      kao_mo3( 4, :,13) = (/ &
     & 4.66456e-01_jprb, 4.71080e-01_jprb, 4.75750e-01_jprb, 4.80467e-01_jprb, 4.85230e-01_jprb, &
     & 4.90040e-01_jprb, 4.94898e-01_jprb, 4.99804e-01_jprb, 5.04759e-01_jprb, 5.09763e-01_jprb, &
     & 5.14817e-01_jprb, 5.19920e-01_jprb, 5.25075e-01_jprb, 5.30280e-01_jprb, 5.35537e-01_jprb, &
     & 5.40846e-01_jprb, 5.46208e-01_jprb, 5.51622e-01_jprb, 5.57091e-01_jprb/)
      kao_mo3( 5, :,13) = (/ &
     & 4.66853e-01_jprb, 4.71456e-01_jprb, 4.76104e-01_jprb, 4.80798e-01_jprb, 4.85539e-01_jprb, &
     & 4.90326e-01_jprb, 4.95160e-01_jprb, 5.00042e-01_jprb, 5.04973e-01_jprb, 5.09952e-01_jprb, &
     & 5.14979e-01_jprb, 5.20057e-01_jprb, 5.25185e-01_jprb, 5.30363e-01_jprb, 5.35592e-01_jprb, &
     & 5.40873e-01_jprb, 5.46205e-01_jprb, 5.51591e-01_jprb, 5.57029e-01_jprb/)
      kao_mo3( 6, :,13) = (/ &
     & 4.66832e-01_jprb, 4.71448e-01_jprb, 4.76110e-01_jprb, 4.80817e-01_jprb, 4.85571e-01_jprb, &
     & 4.90372e-01_jprb, 4.95221e-01_jprb, 5.00118e-01_jprb, 5.05063e-01_jprb, 5.10056e-01_jprb, &
     & 5.15100e-01_jprb, 5.20193e-01_jprb, 5.25336e-01_jprb, 5.30531e-01_jprb, 5.35776e-01_jprb, &
     & 5.41074e-01_jprb, 5.46424e-01_jprb, 5.51826e-01_jprb, 5.57283e-01_jprb/)
      kao_mo3( 7, :,13) = (/ &
     & 4.66679e-01_jprb, 4.71299e-01_jprb, 4.75965e-01_jprb, 4.80677e-01_jprb, 4.85436e-01_jprb, &
     & 4.90241e-01_jprb, 4.95095e-01_jprb, 4.99996e-01_jprb, 5.04946e-01_jprb, 5.09945e-01_jprb, &
     & 5.14993e-01_jprb, 5.20092e-01_jprb, 5.25240e-01_jprb, 5.30440e-01_jprb, 5.35692e-01_jprb, &
     & 5.40995e-01_jprb, 5.46351e-01_jprb, 5.51759e-01_jprb, 5.57222e-01_jprb/)
      kao_mo3( 8, :,13) = (/ &
     & 4.66982e-01_jprb, 4.71598e-01_jprb, 4.76260e-01_jprb, 4.80967e-01_jprb, 4.85721e-01_jprb, &
     & 4.90522e-01_jprb, 4.95371e-01_jprb, 5.00268e-01_jprb, 5.05213e-01_jprb, 5.10206e-01_jprb, &
     & 5.15250e-01_jprb, 5.20343e-01_jprb, 5.25486e-01_jprb, 5.30680e-01_jprb, 5.35926e-01_jprb, &
     & 5.41223e-01_jprb, 5.46573e-01_jprb, 5.51976e-01_jprb, 5.57432e-01_jprb/)
      kao_mo3( 9, :,13) = (/ &
     & 1.13709e-01_jprb, 1.13141e-01_jprb, 1.12576e-01_jprb, 1.12013e-01_jprb, 1.11453e-01_jprb, &
     & 1.10897e-01_jprb, 1.10342e-01_jprb, 1.09791e-01_jprb, 1.09242e-01_jprb, 1.08696e-01_jprb, &
     & 1.08153e-01_jprb, 1.07613e-01_jprb, 1.07075e-01_jprb, 1.06540e-01_jprb, 1.06007e-01_jprb, &
     & 1.05478e-01_jprb, 1.04951e-01_jprb, 1.04426e-01_jprb, 1.03904e-01_jprb/)
      kao_mo3( 1, :,14) = (/ &
     & 5.67608e-01_jprb, 5.55796e-01_jprb, 5.44230e-01_jprb, 5.32904e-01_jprb, 5.21814e-01_jprb, &
     & 5.10955e-01_jprb, 5.00322e-01_jprb, 4.89910e-01_jprb, 4.79714e-01_jprb, 4.69731e-01_jprb, &
     & 4.59956e-01_jprb, 4.50384e-01_jprb, 4.41011e-01_jprb, 4.31834e-01_jprb, 4.22847e-01_jprb, &
     & 4.14048e-01_jprb, 4.05431e-01_jprb, 3.96994e-01_jprb, 3.88732e-01_jprb/)
      kao_mo3( 2, :,14) = (/ &
     & 5.67766e-01_jprb, 5.55948e-01_jprb, 5.44376e-01_jprb, 5.33045e-01_jprb, 5.21950e-01_jprb, &
     & 5.11086e-01_jprb, 5.00448e-01_jprb, 4.90031e-01_jprb, 4.79831e-01_jprb, 4.69844e-01_jprb, &
     & 4.60064e-01_jprb, 4.50488e-01_jprb, 4.41111e-01_jprb, 4.31930e-01_jprb, 4.22939e-01_jprb, &
     & 4.14136e-01_jprb, 4.05516e-01_jprb, 3.97075e-01_jprb, 3.88810e-01_jprb/)
      kao_mo3( 3, :,14) = (/ &
     & 5.67460e-01_jprb, 5.55647e-01_jprb, 5.44080e-01_jprb, 5.32754e-01_jprb, 5.21664e-01_jprb, &
     & 5.10805e-01_jprb, 5.00172e-01_jprb, 4.89760e-01_jprb, 4.79564e-01_jprb, 4.69582e-01_jprb, &
     & 4.59806e-01_jprb, 4.50235e-01_jprb, 4.40862e-01_jprb, 4.31685e-01_jprb, 4.22699e-01_jprb, &
     & 4.13900e-01_jprb, 4.05284e-01_jprb, 3.96847e-01_jprb, 3.88586e-01_jprb/)
      kao_mo3( 4, :,14) = (/ &
     & 5.67925e-01_jprb, 5.56107e-01_jprb, 5.44536e-01_jprb, 5.33205e-01_jprb, 5.22110e-01_jprb, &
     & 5.11246e-01_jprb, 5.00608e-01_jprb, 4.90191e-01_jprb, 4.79991e-01_jprb, 4.70004e-01_jprb, &
     & 4.60224e-01_jprb, 4.50647e-01_jprb, 4.41270e-01_jprb, 4.32088e-01_jprb, 4.23097e-01_jprb, &
     & 4.14293e-01_jprb, 4.05673e-01_jprb, 3.97231e-01_jprb, 3.88966e-01_jprb/)
      kao_mo3( 5, :,14) = (/ &
     & 5.67520e-01_jprb, 5.55733e-01_jprb, 5.44190e-01_jprb, 5.32887e-01_jprb, 5.21818e-01_jprb, &
     & 5.10980e-01_jprb, 5.00366e-01_jprb, 4.89974e-01_jprb, 4.79797e-01_jprb, 4.69831e-01_jprb, &
     & 4.60072e-01_jprb, 4.50516e-01_jprb, 4.41159e-01_jprb, 4.31996e-01_jprb, 4.23023e-01_jprb, &
     & 4.14236e-01_jprb, 4.05633e-01_jprb, 3.97207e-01_jprb, 3.88957e-01_jprb/)
      kao_mo3( 6, :,14) = (/ &
     & 5.67549e-01_jprb, 5.55749e-01_jprb, 5.44195e-01_jprb, 5.32880e-01_jprb, 5.21801e-01_jprb, &
     & 5.10952e-01_jprb, 5.00329e-01_jprb, 4.89927e-01_jprb, 4.79740e-01_jprb, 4.69766e-01_jprb, &
     & 4.59999e-01_jprb, 4.50435e-01_jprb, 4.41070e-01_jprb, 4.31900e-01_jprb, 4.22920e-01_jprb, &
     & 4.14127e-01_jprb, 4.05517e-01_jprb, 3.97086e-01_jprb, 3.88830e-01_jprb/)
      kao_mo3( 7, :,14) = (/ &
     & 5.67727e-01_jprb, 5.55909e-01_jprb, 5.44337e-01_jprb, 5.33005e-01_jprb, 5.21910e-01_jprb, &
     & 5.11046e-01_jprb, 5.00408e-01_jprb, 4.89991e-01_jprb, 4.79791e-01_jprb, 4.69804e-01_jprb, &
     & 4.60024e-01_jprb, 4.50448e-01_jprb, 4.41071e-01_jprb, 4.31890e-01_jprb, 4.22900e-01_jprb, &
     & 4.14096e-01_jprb, 4.05476e-01_jprb, 3.97036e-01_jprb, 3.88771e-01_jprb/)
      kao_mo3( 8, :,14) = (/ &
     & 5.67795e-01_jprb, 5.55965e-01_jprb, 5.44381e-01_jprb, 5.33039e-01_jprb, 5.21933e-01_jprb, &
     & 5.11058e-01_jprb, 5.00410e-01_jprb, 4.89984e-01_jprb, 4.79775e-01_jprb, 4.69779e-01_jprb, &
     & 4.59991e-01_jprb, 4.50407e-01_jprb, 4.41023e-01_jprb, 4.31834e-01_jprb, 4.22836e-01_jprb, &
     & 4.14026e-01_jprb, 4.05400e-01_jprb, 3.96953e-01_jprb, 3.88683e-01_jprb/)
      kao_mo3( 9, :,14) = (/ &
     & 1.32957e-01_jprb, 1.31737e-01_jprb, 1.30528e-01_jprb, 1.29330e-01_jprb, 1.28143e-01_jprb, &
     & 1.26967e-01_jprb, 1.25802e-01_jprb, 1.24648e-01_jprb, 1.23504e-01_jprb, 1.22370e-01_jprb, &
     & 1.21247e-01_jprb, 1.20135e-01_jprb, 1.19032e-01_jprb, 1.17940e-01_jprb, 1.16857e-01_jprb, &
     & 1.15785e-01_jprb, 1.14722e-01_jprb, 1.13669e-01_jprb, 1.12626e-01_jprb/)
      kao_mo3( 1, :,15) = (/ &
     & 1.51281e-01_jprb, 1.53439e-01_jprb, 1.55628e-01_jprb, 1.57848e-01_jprb, 1.60100e-01_jprb, &
     & 1.62384e-01_jprb, 1.64700e-01_jprb, 1.67049e-01_jprb, 1.69432e-01_jprb, 1.71849e-01_jprb, &
     & 1.74301e-01_jprb, 1.76787e-01_jprb, 1.79309e-01_jprb, 1.81866e-01_jprb, 1.84461e-01_jprb, &
     & 1.87092e-01_jprb, 1.89761e-01_jprb, 1.92468e-01_jprb, 1.95213e-01_jprb/)
      kao_mo3( 2, :,15) = (/ &
     & 1.51431e-01_jprb, 1.53589e-01_jprb, 1.55778e-01_jprb, 1.57998e-01_jprb, 1.60250e-01_jprb, &
     & 1.62534e-01_jprb, 1.64850e-01_jprb, 1.67199e-01_jprb, 1.69582e-01_jprb, 1.71999e-01_jprb, &
     & 1.74450e-01_jprb, 1.76937e-01_jprb, 1.79458e-01_jprb, 1.82016e-01_jprb, 1.84610e-01_jprb, &
     & 1.87241e-01_jprb, 1.89909e-01_jprb, 1.92616e-01_jprb, 1.95361e-01_jprb/)
      kao_mo3( 3, :,15) = (/ &
     & 1.51299e-01_jprb, 1.53461e-01_jprb, 1.55654e-01_jprb, 1.57878e-01_jprb, 1.60134e-01_jprb, &
     & 1.62422e-01_jprb, 1.64744e-01_jprb, 1.67098e-01_jprb, 1.69486e-01_jprb, 1.71908e-01_jprb, &
     & 1.74364e-01_jprb, 1.76856e-01_jprb, 1.79383e-01_jprb, 1.81947e-01_jprb, 1.84547e-01_jprb, &
     & 1.87184e-01_jprb, 1.89859e-01_jprb, 1.92572e-01_jprb, 1.95324e-01_jprb/)
      kao_mo3( 4, :,15) = (/ &
     & 1.51281e-01_jprb, 1.53439e-01_jprb, 1.55628e-01_jprb, 1.57848e-01_jprb, 1.60100e-01_jprb, &
     & 1.62384e-01_jprb, 1.64700e-01_jprb, 1.67049e-01_jprb, 1.69432e-01_jprb, 1.71849e-01_jprb, &
     & 1.74301e-01_jprb, 1.76787e-01_jprb, 1.79309e-01_jprb, 1.81866e-01_jprb, 1.84461e-01_jprb, &
     & 1.87092e-01_jprb, 1.89761e-01_jprb, 1.92468e-01_jprb, 1.95213e-01_jprb/)
      kao_mo3( 5, :,15) = (/ &
     & 1.51281e-01_jprb, 1.53439e-01_jprb, 1.55628e-01_jprb, 1.57848e-01_jprb, 1.60100e-01_jprb, &
     & 1.62384e-01_jprb, 1.64700e-01_jprb, 1.67049e-01_jprb, 1.69432e-01_jprb, 1.71849e-01_jprb, &
     & 1.74301e-01_jprb, 1.76787e-01_jprb, 1.79309e-01_jprb, 1.81866e-01_jprb, 1.84461e-01_jprb, &
     & 1.87092e-01_jprb, 1.89761e-01_jprb, 1.92468e-01_jprb, 1.95213e-01_jprb/)
      kao_mo3( 6, :,15) = (/ &
     & 1.51299e-01_jprb, 1.53461e-01_jprb, 1.55654e-01_jprb, 1.57878e-01_jprb, 1.60134e-01_jprb, &
     & 1.62422e-01_jprb, 1.64744e-01_jprb, 1.67098e-01_jprb, 1.69486e-01_jprb, 1.71908e-01_jprb, &
     & 1.74364e-01_jprb, 1.76856e-01_jprb, 1.79383e-01_jprb, 1.81947e-01_jprb, 1.84547e-01_jprb, &
     & 1.87184e-01_jprb, 1.89859e-01_jprb, 1.92572e-01_jprb, 1.95324e-01_jprb/)
      kao_mo3( 7, :,15) = (/ &
     & 1.51299e-01_jprb, 1.53461e-01_jprb, 1.55654e-01_jprb, 1.57878e-01_jprb, 1.60134e-01_jprb, &
     & 1.62422e-01_jprb, 1.64744e-01_jprb, 1.67098e-01_jprb, 1.69486e-01_jprb, 1.71908e-01_jprb, &
     & 1.74364e-01_jprb, 1.76856e-01_jprb, 1.79383e-01_jprb, 1.81947e-01_jprb, 1.84547e-01_jprb, &
     & 1.87184e-01_jprb, 1.89859e-01_jprb, 1.92572e-01_jprb, 1.95324e-01_jprb/)
      kao_mo3( 8, :,15) = (/ &
     & 1.51281e-01_jprb, 1.53439e-01_jprb, 1.55628e-01_jprb, 1.57848e-01_jprb, 1.60100e-01_jprb, &
     & 1.62384e-01_jprb, 1.64700e-01_jprb, 1.67049e-01_jprb, 1.69432e-01_jprb, 1.71849e-01_jprb, &
     & 1.74301e-01_jprb, 1.76787e-01_jprb, 1.79309e-01_jprb, 1.81866e-01_jprb, 1.84461e-01_jprb, &
     & 1.87092e-01_jprb, 1.89761e-01_jprb, 1.92468e-01_jprb, 1.95213e-01_jprb/)
      kao_mo3( 9, :,15) = (/ &
     & 2.44180e-01_jprb, 2.35686e-01_jprb, 2.27487e-01_jprb, 2.19574e-01_jprb, 2.11935e-01_jprb, &
     & 2.04563e-01_jprb, 1.97447e-01_jprb, 1.90578e-01_jprb, 1.83949e-01_jprb, 1.77550e-01_jprb, &
     & 1.71373e-01_jprb, 1.65412e-01_jprb, 1.59658e-01_jprb, 1.54104e-01_jprb, 1.48743e-01_jprb, &
     & 1.43569e-01_jprb, 1.38574e-01_jprb, 1.33754e-01_jprb, 1.29101e-01_jprb/)
      kao_mo3( 1, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 2, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 3, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 4, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 5, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 6, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 7, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 8, :,16) = (/ &
     & 1.02934e-01_jprb, 1.04369e-01_jprb, 1.05825e-01_jprb, 1.07300e-01_jprb, 1.08797e-01_jprb, &
     & 1.10314e-01_jprb, 1.11852e-01_jprb, 1.13412e-01_jprb, 1.14993e-01_jprb, 1.16597e-01_jprb, &
     & 1.18223e-01_jprb, 1.19871e-01_jprb, 1.21543e-01_jprb, 1.23238e-01_jprb, 1.24956e-01_jprb, &
     & 1.26699e-01_jprb, 1.28466e-01_jprb, 1.30257e-01_jprb, 1.32073e-01_jprb/)
      kao_mo3( 9, :,16) = (/ &
     & 3.91531e-01_jprb, 3.78978e-01_jprb, 3.66827e-01_jprb, 3.55067e-01_jprb, 3.43683e-01_jprb, &
     & 3.32664e-01_jprb, 3.21999e-01_jprb, 3.11675e-01_jprb, 3.01683e-01_jprb, 2.92011e-01_jprb, &
     & 2.82648e-01_jprb, 2.73586e-01_jprb, 2.64815e-01_jprb, 2.56325e-01_jprb, 2.48107e-01_jprb, &
     & 2.40152e-01_jprb, 2.32453e-01_jprb, 2.25000e-01_jprb, 2.17787e-01_jprb/)


!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.0689e-05_jprb,1.6987e-05_jprb,1.8993e-05_jprb,3.4470e-05_jprb,4.0873e-05_jprb,4.8275e-05_jprb, &
     &6.1178e-05_jprb,6.4035e-05_jprb,6.6253e-05_jprb,7.8914e-05_jprb,8.1640e-05_jprb,7.9738e-05_jprb, &
     &7.8492e-05_jprb,9.1565e-05_jprb,1.0262e-04_jprb,1.0368e-04_jprb/)
      forrefo(2,:) = (/ &
     &1.1194e-05_jprb,1.6128e-05_jprb,1.7213e-05_jprb,2.6845e-05_jprb,4.1361e-05_jprb,5.1508e-05_jprb, &
     &6.8245e-05_jprb,7.4063e-05_jprb,7.6273e-05_jprb,8.4061e-05_jprb,8.2492e-05_jprb,8.1720e-05_jprb, &
     &7.7626e-05_jprb,1.0096e-04_jprb,1.0519e-04_jprb,1.0631e-04_jprb/)
      forrefo(3,:) = (/ &
     &1.0891e-05_jprb,1.4933e-05_jprb,1.7964e-05_jprb,2.2577e-05_jprb,4.4290e-05_jprb,5.4675e-05_jprb, &
     &7.2494e-05_jprb,7.8410e-05_jprb,7.6948e-05_jprb,7.5742e-05_jprb,7.7654e-05_jprb,8.2760e-05_jprb, &
     &7.8443e-05_jprb,9.8384e-05_jprb,1.0634e-04_jprb,1.0838e-04_jprb/)
      forrefo(4,:) = (/ &
     &1.1316e-05_jprb,1.5470e-05_jprb,2.1246e-05_jprb,3.3349e-05_jprb,4.8704e-05_jprb,5.6424e-05_jprb, &
     &5.8569e-05_jprb,5.8780e-05_jprb,6.0358e-05_jprb,6.1586e-05_jprb,6.4281e-05_jprb,6.9333e-05_jprb, &
     &7.2763e-05_jprb,7.2675e-05_jprb,7.3754e-05_jprb,1.0131e-04_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).
      selfrefo(:, 1) = (/ &
     & 1.27686e-01_jprb, 1.09347e-01_jprb, 9.36410e-02_jprb, 8.01912e-02_jprb, 6.86732e-02_jprb, &
     & 5.88096e-02_jprb, 5.03627e-02_jprb, 4.31290e-02_jprb, 3.69343e-02_jprb, 3.16294e-02_jprb/)
      selfrefo(:, 2) = (/ &
     & 1.40051e-01_jprb, 1.20785e-01_jprb, 1.04170e-01_jprb, 8.98402e-02_jprb, 7.74816e-02_jprb, &
     & 6.68231e-02_jprb, 5.76308e-02_jprb, 4.97030e-02_jprb, 4.28658e-02_jprb, 3.69691e-02_jprb/)
      selfrefo(:, 3) = (/ &
     & 1.42322e-01_jprb, 1.22872e-01_jprb, 1.06080e-01_jprb, 9.15829e-02_jprb, 7.90671e-02_jprb, &
     & 6.82616e-02_jprb, 5.89329e-02_jprb, 5.08790e-02_jprb, 4.39258e-02_jprb, 3.79228e-02_jprb/)
      selfrefo(:, 4) = (/ &
     & 1.53244e-01_jprb, 1.33057e-01_jprb, 1.15530e-01_jprb, 1.00311e-01_jprb, 8.70977e-02_jprb, &
     & 7.56244e-02_jprb, 6.56626e-02_jprb, 5.70130e-02_jprb, 4.95028e-02_jprb, 4.29819e-02_jprb/)
      selfrefo(:, 5) = (/ &
     & 1.71011e-01_jprb, 1.46680e-01_jprb, 1.25810e-01_jprb, 1.07910e-01_jprb, 9.25563e-02_jprb, &
     & 7.93874e-02_jprb, 6.80922e-02_jprb, 5.84040e-02_jprb, 5.00943e-02_jprb, 4.29669e-02_jprb/)
      selfrefo(:, 6) = (/ &
     & 1.76012e-01_jprb, 1.51010e-01_jprb, 1.29560e-01_jprb, 1.11157e-01_jprb, 9.53672e-02_jprb, &
     & 8.18207e-02_jprb, 7.01984e-02_jprb, 6.02270e-02_jprb, 5.16720e-02_jprb, 4.43322e-02_jprb/)
      selfrefo(:, 7) = (/ &
     & 1.85600e-01_jprb, 1.59051e-01_jprb, 1.36300e-01_jprb, 1.16803e-01_jprb, 1.00095e-01_jprb, &
     & 8.57776e-02_jprb, 7.35077e-02_jprb, 6.29930e-02_jprb, 5.39823e-02_jprb, 4.62606e-02_jprb/)
      selfrefo(:, 8) = (/ &
     & 1.88931e-01_jprb, 1.61727e-01_jprb, 1.38440e-01_jprb, 1.18506e-01_jprb, 1.01442e-01_jprb, &
     & 8.68356e-02_jprb, 7.43321e-02_jprb, 6.36290e-02_jprb, 5.44670e-02_jprb, 4.66243e-02_jprb/)
      selfrefo(:, 9) = (/ &
     & 1.91122e-01_jprb, 1.63407e-01_jprb, 1.39710e-01_jprb, 1.19450e-01_jprb, 1.02128e-01_jprb, &
     & 8.73176e-02_jprb, 7.46552e-02_jprb, 6.38290e-02_jprb, 5.45728e-02_jprb, 4.66589e-02_jprb/)
      selfrefo(:,10) = (/ &
     & 1.91334e-01_jprb, 1.64872e-01_jprb, 1.42070e-01_jprb, 1.22421e-01_jprb, 1.05490e-01_jprb, &
     & 9.09008e-02_jprb, 7.83291e-02_jprb, 6.74960e-02_jprb, 5.81612e-02_jprb, 5.01174e-02_jprb/)
      selfrefo(:,11) = (/ &
     & 1.89858e-01_jprb, 1.63934e-01_jprb, 1.41550e-01_jprb, 1.22222e-01_jprb, 1.05534e-01_jprb, &
     & 9.11237e-02_jprb, 7.86814e-02_jprb, 6.79380e-02_jprb, 5.86615e-02_jprb, 5.06517e-02_jprb/)
      selfrefo(:,12) = (/ &
     & 1.89783e-01_jprb, 1.63757e-01_jprb, 1.41300e-01_jprb, 1.21923e-01_jprb, 1.05203e-01_jprb, &
     & 9.07760e-02_jprb, 7.83274e-02_jprb, 6.75860e-02_jprb, 5.83176e-02_jprb, 5.03202e-02_jprb/)
      selfrefo(:,13) = (/ &
     & 1.87534e-01_jprb, 1.62016e-01_jprb, 1.39970e-01_jprb, 1.20924e-01_jprb, 1.04470e-01_jprb, &
     & 9.02541e-02_jprb, 7.79730e-02_jprb, 6.73630e-02_jprb, 5.81967e-02_jprb, 5.02778e-02_jprb/)
      selfrefo(:,14) = (/ &
     & 1.99128e-01_jprb, 1.71410e-01_jprb, 1.47550e-01_jprb, 1.27011e-01_jprb, 1.09332e-01_jprb, &
     & 9.41131e-02_jprb, 8.10128e-02_jprb, 6.97360e-02_jprb, 6.00289e-02_jprb, 5.16731e-02_jprb/)
      selfrefo(:,15) = (/ &
     & 1.99460e-01_jprb, 1.72342e-01_jprb, 1.48910e-01_jprb, 1.28664e-01_jprb, 1.11171e-01_jprb, &
     & 9.60560e-02_jprb, 8.29962e-02_jprb, 7.17120e-02_jprb, 6.19620e-02_jprb, 5.35376e-02_jprb/)
      selfrefo(:,16) = (/ &
     & 1.99906e-01_jprb, 1.72737e-01_jprb, 1.49260e-01_jprb, 1.28974e-01_jprb, 1.11445e-01_jprb, &
     & 9.62982e-02_jprb, 8.32102e-02_jprb, 7.19010e-02_jprb, 6.21288e-02_jprb, 5.36848e-02_jprb/)

if (lhook) call dr_hook('rrtm_kgb5',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb5:error reading file radrrtm")

end subroutine rrtm_kgb5
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

