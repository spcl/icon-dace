! # 1 "ifsrrtm/rrtm_kgb8.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb8.f90"
subroutine rrtm_kgb8

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo updated to rrtmg v4.85
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad
use yommp0    , only : nproc, myproc

use yoerrto8 , only : kao     ,kbo       ,selfrefo ,forrefo,fracrefao ,&
 & fracrefbo, cfc12o  ,cfc22adjo ,kao_mco2,kbo_mco2,&
 & kao_mn2o,kbo_mn2o,kao_mo3  , kao_d, kbo_d


!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 31 "ifsrrtm/rrtm_kgb8.f90" 2

if (lhook) call dr_hook('rrtm_kgb8',0,zhook_handle)

if( myproc==1 )then
  read(nulrad) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb8:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb8:')
endif


! planck fraction mapping level : p=473.4280 mb, t = 259.83 k
      fracrefao(:) = (/ &
      &  1.6004e-01_jprb,1.5437e-01_jprb,1.4502e-01_jprb,1.3084e-01_jprb,1.1523e-01_jprb,9.7743e-02_jprb, &
      &  8.0376e-02_jprb,6.0261e-02_jprb,4.1111e-02_jprb,4.4772e-03_jprb,3.6511e-03_jprb,2.9154e-03_jprb, &
      &  2.1184e-03_jprb,1.3048e-03_jprb,4.6637e-04_jprb,6.5624e-05_jprb/)

! planck fraction mapping level : p=95.5835 mb, t= 215.7 k
      fracrefbo(:) = (/ &
      &  1.4987e-01_jprb,1.4665e-01_jprb,1.4154e-01_jprb,1.3200e-01_jprb,1.1902e-01_jprb,1.0352e-01_jprb, &
      &  8.4939e-02_jprb,6.4105e-02_jprb,4.3190e-02_jprb,4.5129e-03_jprb,3.7656e-03_jprb,2.8733e-03_jprb, &
      &  2.0947e-03_jprb,1.3201e-03_jprb,5.1832e-04_jprb,7.7473e-05_jprb/)

! minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

cfc12o( :) = (/&
 & 85.4027_jprb, 89.4696_jprb, 74.0959_jprb, 67.7480_jprb,&
 & 61.2444_jprb, 59.9073_jprb, 60.8296_jprb, 63.0998_jprb,&
 & 59.6110_jprb, 64.0735_jprb, 57.2622_jprb, 58.9721_jprb,&
 & 43.5505_jprb, 26.1192_jprb, 32.7023_jprb, 32.8667_jprb/)  

!     original cfc22 is multiplied by 1.485 to account for the 780-850 cm-1 
!     and 1290-1335 cm-1 bands.
cfc22adjo( :) = (/&
 & 135.335_jprb, 89.6642_jprb, 76.2375_jprb, 65.9748_jprb,&
 & 63.1164_jprb, 60.2935_jprb, 64.0299_jprb, 75.4264_jprb,&
 & 51.3018_jprb, 7.07911_jprb, 5.86928_jprb, 0.398693_jprb,&
 & 2.82885_jprb, 9.12751_jprb, 6.28271_jprb, 0.0_jprb/)  

!     ------------------------------------------------------------------

!     the array kao contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  the first
!     index in the array, jt, which runs from 1 to 5, corresponds to 
!     different temperatures.  more specifically, jt = 3 means that the 
!     data are for the corresponding tref for this  pressure level, 
!     jt = 2 refers to the temperaturetref-15, jt = 1 is for tref-30, 
!     jt = 4 is for tref+15, and jt = 5 is for tref+30.  the second 
!     index, jp, runs from 1 to 13 and refers to the corresponding 
!     pressure level in pref (e.g. jp = 1 is for a pressure of 1053.63 mb).  
!     the third index, ig, goes from 1 to 16, and tells us which 
!     g-interval the absorption coefficients are for.
!     the array ka contains absorption coef5s at the 16 chosen g-values 
!     for a range of pressure levels > ~100mb and temperatures.  the first
!     index in the array, jt, which runs from 1 to 5, corresponds to 
!     different temperatures.  more specifically, jt = 3 means that the 
!     data are for the cooresponding tref for this  pressure level, 
!     jt = 2 refers to the temperature
!     tref-15, jt = 1 is for tref-30, jt = 4 is for tref+15, and jt = 5
!     is for tref+30.  the second index, jp, runs from 1 to 13 and refers
!     to the corresponding pressure level in pref (e.g. jp = 1 is for a
!     pressure of 1053.63 mb).  the third index, ig, goes from 1 to 16,
!     and tells us which "g-channel" the absorption coefficients are for.



!     the array kbo contains absorption coefs at the 16 chosen g-values 
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



!     the array kao_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kao_mco2(:, 1) = (/ &
     & 8.88964e-07_jprb, 1.13087e-06_jprb, 1.43861e-06_jprb, 1.83010e-06_jprb, 2.32811e-06_jprb, &
     & 2.96165e-06_jprb, 3.76760e-06_jprb, 4.79286e-06_jprb, 6.09712e-06_jprb, 7.75630e-06_jprb, &
     & 9.86699e-06_jprb, 1.25521e-05_jprb, 1.59678e-05_jprb, 2.03130e-05_jprb, 2.58407e-05_jprb, &
     & 3.28727e-05_jprb, 4.18182e-05_jprb, 5.31980e-05_jprb, 6.76745e-05_jprb/)
      kao_mco2(:, 2) = (/ &
     & 1.10492e-05_jprb, 1.35003e-05_jprb, 1.64952e-05_jprb, 2.01545e-05_jprb, 2.46256e-05_jprb, &
     & 3.00885e-05_jprb, 3.67632e-05_jprb, 4.49188e-05_jprb, 5.48835e-05_jprb, 6.70588e-05_jprb, &
     & 8.19351e-05_jprb, 1.00111e-04_jprb, 1.22320e-04_jprb, 1.49455e-04_jprb, 1.82610e-04_jprb, &
     & 2.23121e-04_jprb, 2.72618e-04_jprb, 3.33095e-04_jprb, 4.06988e-04_jprb/)
      kao_mco2(:, 3) = (/ &
     & 1.51034e-05_jprb, 1.81249e-05_jprb, 2.17508e-05_jprb, 2.61020e-05_jprb, 3.13238e-05_jprb, &
     & 3.75901e-05_jprb, 4.51101e-05_jprb, 5.41344e-05_jprb, 6.49640e-05_jprb, 7.79601e-05_jprb, &
     & 9.35562e-05_jprb, 1.12272e-04_jprb, 1.34732e-04_jprb, 1.61686e-04_jprb, 1.94031e-04_jprb, &
     & 2.32847e-04_jprb, 2.79429e-04_jprb, 3.35329e-04_jprb, 4.02411e-04_jprb/)
      kao_mco2(:, 4) = (/ &
     & 1.57088e-05_jprb, 1.89537e-05_jprb, 2.28688e-05_jprb, 2.75928e-05_jprb, 3.32924e-05_jprb, &
     & 4.01695e-05_jprb, 4.84671e-05_jprb, 5.84787e-05_jprb, 7.05584e-05_jprb, 8.51332e-05_jprb, &
     & 1.02719e-04_jprb, 1.23937e-04_jprb, 1.49538e-04_jprb, 1.80427e-04_jprb, 2.17697e-04_jprb, &
     & 2.62666e-04_jprb, 3.16923e-04_jprb, 3.82388e-04_jprb, 4.61376e-04_jprb/)
      kao_mco2(:, 5) = (/ &
     & 3.09299e-05_jprb, 3.73196e-05_jprb, 4.50294e-05_jprb, 5.43320e-05_jprb, 6.55563e-05_jprb, &
     & 7.90995e-05_jprb, 9.54405e-05_jprb, 1.15157e-04_jprb, 1.38948e-04_jprb, 1.67652e-04_jprb, &
     & 2.02288e-04_jprb, 2.44078e-04_jprb, 2.94501e-04_jprb, 3.55342e-04_jprb, 4.28751e-04_jprb, &
     & 5.17327e-04_jprb, 6.24200e-04_jprb, 7.53153e-04_jprb, 9.08745e-04_jprb/)
      kao_mco2(:, 6) = (/ &
     & 1.98653e-05_jprb, 2.38878e-05_jprb, 2.87248e-05_jprb, 3.45413e-05_jprb, 4.15355e-05_jprb, &
     & 4.99459e-05_jprb, 6.00593e-05_jprb, 7.22206e-05_jprb, 8.68445e-05_jprb, 1.04429e-04_jprb, &
     & 1.25575e-04_jprb, 1.51003e-04_jprb, 1.81579e-04_jprb, 2.18346e-04_jprb, 2.62559e-04_jprb, &
     & 3.15724e-04_jprb, 3.79654e-04_jprb, 4.56529e-04_jprb, 5.48971e-04_jprb/)
      kao_mco2(:, 7) = (/ &
     & 1.54276e-06_jprb, 1.90144e-06_jprb, 2.34351e-06_jprb, 2.88836e-06_jprb, 3.55989e-06_jprb, &
     & 4.38754e-06_jprb, 5.40761e-06_jprb, 6.66485e-06_jprb, 8.21439e-06_jprb, 1.01242e-05_jprb, &
     & 1.24780e-05_jprb, 1.53790e-05_jprb, 1.89546e-05_jprb, 2.33614e-05_jprb, 2.87928e-05_jprb, &
     & 3.54869e-05_jprb, 4.37374e-05_jprb, 5.39060e-05_jprb, 6.64388e-05_jprb/)
      kao_mco2(:, 8) = (/ &
     & 1.66907e-06_jprb, 2.11106e-06_jprb, 2.67008e-06_jprb, 3.37714e-06_jprb, 4.27143e-06_jprb, &
     & 5.40254e-06_jprb, 6.83318e-06_jprb, 8.64266e-06_jprb, 1.09313e-05_jprb, 1.38260e-05_jprb, &
     & 1.74872e-05_jprb, 2.21180e-05_jprb, 2.79750e-05_jprb, 3.53830e-05_jprb, 4.47527e-05_jprb, &
     & 5.66036e-05_jprb, 7.15927e-05_jprb, 9.05509e-05_jprb, 1.14529e-04_jprb/)
      kao_mco2(:, 9) = (/ &
     & 1.22817e-06_jprb, 1.56416e-06_jprb, 1.99206e-06_jprb, 2.53703e-06_jprb, 3.23108e-06_jprb, &
     & 4.11501e-06_jprb, 5.24074e-06_jprb, 6.67445e-06_jprb, 8.50037e-06_jprb, 1.08258e-05_jprb, &
     & 1.37874e-05_jprb, 1.75592e-05_jprb, 2.23629e-05_jprb, 2.84807e-05_jprb, 3.62721e-05_jprb, &
     & 4.61950e-05_jprb, 5.88325e-05_jprb, 7.49272e-05_jprb, 9.54249e-05_jprb/)
      kao_mco2(:,10) = (/ &
     & 3.45943e-08_jprb, 3.84726e-08_jprb, 4.27856e-08_jprb, 4.75821e-08_jprb, 5.29164e-08_jprb, &
     & 5.88487e-08_jprb, 6.54460e-08_jprb, 7.27829e-08_jprb, 8.09423e-08_jprb, 9.00164e-08_jprb, &
     & 1.00108e-07_jprb, 1.11331e-07_jprb, 1.23811e-07_jprb, 1.37691e-07_jprb, 1.53128e-07_jprb, &
     & 1.70294e-07_jprb, 1.89385e-07_jprb, 2.10616e-07_jprb, 2.34228e-07_jprb/)
      kao_mco2(:,11) = (/ &
     & 2.89971e-08_jprb, 3.35110e-08_jprb, 3.87275e-08_jprb, 4.47561e-08_jprb, 5.17230e-08_jprb, &
     & 5.97745e-08_jprb, 6.90794e-08_jprb, 7.98327e-08_jprb, 9.22599e-08_jprb, 1.06622e-07_jprb, &
     & 1.23219e-07_jprb, 1.42400e-07_jprb, 1.64567e-07_jprb, 1.90184e-07_jprb, 2.19789e-07_jprb, &
     & 2.54003e-07_jprb, 2.93542e-07_jprb, 3.39237e-07_jprb, 3.92044e-07_jprb/)
      kao_mco2(:,12) = (/ &
     & 2.51330e-08_jprb, 2.96783e-08_jprb, 3.50457e-08_jprb, 4.13837e-08_jprb, 4.88679e-08_jprb, &
     & 5.77056e-08_jprb, 6.81416e-08_jprb, 8.04650e-08_jprb, 9.50171e-08_jprb, 1.12201e-07_jprb, &
     & 1.32492e-07_jprb, 1.56454e-07_jprb, 1.84748e-07_jprb, 2.18160e-07_jprb, 2.57614e-07_jprb, &
     & 3.04203e-07_jprb, 3.59218e-07_jprb, 4.24182e-07_jprb, 5.00895e-07_jprb/)
      kao_mco2(:,13) = (/ &
     & 1.16966e-07_jprb, 1.13960e-07_jprb, 1.11032e-07_jprb, 1.08179e-07_jprb, 1.05400e-07_jprb, &
     & 1.02691e-07_jprb, 1.00053e-07_jprb, 9.74820e-08_jprb, 9.49772e-08_jprb, 9.25368e-08_jprb, &
     & 9.01591e-08_jprb, 8.78425e-08_jprb, 8.55854e-08_jprb, 8.33863e-08_jprb, 8.12437e-08_jprb, &
     & 7.91562e-08_jprb, 7.71223e-08_jprb, 7.51407e-08_jprb, 7.32100e-08_jprb/)
      kao_mco2(:,14) = (/ &
     & 9.17853e-08_jprb, 8.94322e-08_jprb, 8.71395e-08_jprb, 8.49055e-08_jprb, 8.27289e-08_jprb, &
     & 8.06080e-08_jprb, 7.85415e-08_jprb, 7.65279e-08_jprb, 7.45660e-08_jprb, 7.26544e-08_jprb, &
     & 7.07918e-08_jprb, 6.89770e-08_jprb, 6.72086e-08_jprb, 6.54856e-08_jprb, 6.38068e-08_jprb, &
     & 6.21710e-08_jprb, 6.05772e-08_jprb, 5.90242e-08_jprb, 5.75110e-08_jprb/)
      kao_mco2(:,15) = (/ &
     & 8.34607e-08_jprb, 8.13236e-08_jprb, 7.92413e-08_jprb, 7.72122e-08_jprb, 7.52351e-08_jprb, &
     & 7.33087e-08_jprb, 7.14315e-08_jprb, 6.96025e-08_jprb, 6.78202e-08_jprb, 6.60837e-08_jprb, &
     & 6.43915e-08_jprb, 6.27427e-08_jprb, 6.11361e-08_jprb, 5.95707e-08_jprb, 5.80453e-08_jprb, &
     & 5.65590e-08_jprb, 5.51108e-08_jprb, 5.36996e-08_jprb, 5.23246e-08_jprb/)
      kao_mco2(:,16) = (/ &
     & 8.34607e-08_jprb, 8.13236e-08_jprb, 7.92413e-08_jprb, 7.72122e-08_jprb, 7.52351e-08_jprb, &
     & 7.33087e-08_jprb, 7.14315e-08_jprb, 6.96025e-08_jprb, 6.78202e-08_jprb, 6.60837e-08_jprb, &
     & 6.43915e-08_jprb, 6.27427e-08_jprb, 6.11361e-08_jprb, 5.95707e-08_jprb, 5.80453e-08_jprb, &
     & 5.65590e-08_jprb, 5.51108e-08_jprb, 5.36996e-08_jprb, 5.23246e-08_jprb/)

      kao_mo3(:, 1) = (/ &
     & 1.18276e-01_jprb, 1.18009e-01_jprb, 1.17742e-01_jprb, 1.17476e-01_jprb, 1.17210e-01_jprb, &
     & 1.16945e-01_jprb, 1.16681e-01_jprb, 1.16417e-01_jprb, 1.16153e-01_jprb, 1.15891e-01_jprb, &
     & 1.15629e-01_jprb, 1.15367e-01_jprb, 1.15106e-01_jprb, 1.14846e-01_jprb, 1.14586e-01_jprb, &
     & 1.14327e-01_jprb, 1.14069e-01_jprb, 1.13811e-01_jprb, 1.13553e-01_jprb/)
      kao_mo3(:, 2) = (/ &
     & 1.83777e-01_jprb, 1.84268e-01_jprb, 1.84761e-01_jprb, 1.85255e-01_jprb, 1.85751e-01_jprb, &
     & 1.86248e-01_jprb, 1.86746e-01_jprb, 1.87245e-01_jprb, 1.87746e-01_jprb, 1.88248e-01_jprb, &
     & 1.88752e-01_jprb, 1.89257e-01_jprb, 1.89763e-01_jprb, 1.90270e-01_jprb, 1.90779e-01_jprb, &
     & 1.91290e-01_jprb, 1.91801e-01_jprb, 1.92314e-01_jprb, 1.92829e-01_jprb/)
      kao_mo3(:, 3) = (/ &
     & 2.33414e-01_jprb, 2.34511e-01_jprb, 2.35614e-01_jprb, 2.36722e-01_jprb, 2.37836e-01_jprb, &
     & 2.38954e-01_jprb, 2.40078e-01_jprb, 2.41207e-01_jprb, 2.42342e-01_jprb, 2.43481e-01_jprb, &
     & 2.44626e-01_jprb, 2.45777e-01_jprb, 2.46933e-01_jprb, 2.48094e-01_jprb, 2.49261e-01_jprb, &
     & 2.50433e-01_jprb, 2.51611e-01_jprb, 2.52794e-01_jprb, 2.53983e-01_jprb/)
      kao_mo3(:, 4) = (/ &
     & 2.84906e-01_jprb, 2.87358e-01_jprb, 2.89832e-01_jprb, 2.92328e-01_jprb, 2.94844e-01_jprb, &
     & 2.97383e-01_jprb, 2.99943e-01_jprb, 3.02525e-01_jprb, 3.05130e-01_jprb, 3.07757e-01_jprb, &
     & 3.10406e-01_jprb, 3.13078e-01_jprb, 3.15774e-01_jprb, 3.18492e-01_jprb, 3.21234e-01_jprb, &
     & 3.24000e-01_jprb, 3.26789e-01_jprb, 3.29603e-01_jprb, 3.32440e-01_jprb/)
      kao_mo3(:, 5) = (/ &
     & 3.40508e-01_jprb, 3.44095e-01_jprb, 3.47720e-01_jprb, 3.51383e-01_jprb, 3.55084e-01_jprb, &
     & 3.58824e-01_jprb, 3.62604e-01_jprb, 3.66424e-01_jprb, 3.70284e-01_jprb, 3.74184e-01_jprb, &
     & 3.78126e-01_jprb, 3.82109e-01_jprb, 3.86134e-01_jprb, 3.90202e-01_jprb, 3.94312e-01_jprb, &
     & 3.98466e-01_jprb, 4.02663e-01_jprb, 4.06905e-01_jprb, 4.11191e-01_jprb/)
      kao_mo3(:, 6) = (/ &
     & 3.78368e-01_jprb, 3.83690e-01_jprb, 3.89086e-01_jprb, 3.94558e-01_jprb, 4.00107e-01_jprb, &
     & 4.05735e-01_jprb, 4.11441e-01_jprb, 4.17227e-01_jprb, 4.23095e-01_jprb, 4.29046e-01_jprb, &
     & 4.35080e-01_jprb, 4.41199e-01_jprb, 4.47404e-01_jprb, 4.53697e-01_jprb, 4.60078e-01_jprb, &
     & 4.66548e-01_jprb, 4.73110e-01_jprb, 4.79764e-01_jprb, 4.86511e-01_jprb/)
      kao_mo3(:, 7) = (/ &
     & 4.51965e-01_jprb, 4.58461e-01_jprb, 4.65051e-01_jprb, 4.71735e-01_jprb, 4.78516e-01_jprb, &
     & 4.85394e-01_jprb, 4.92371e-01_jprb, 4.99448e-01_jprb, 5.06627e-01_jprb, 5.13909e-01_jprb, &
     & 5.21296e-01_jprb, 5.28789e-01_jprb, 5.36390e-01_jprb, 5.44100e-01_jprb, 5.51920e-01_jprb, &
     & 5.59854e-01_jprb, 5.67901e-01_jprb, 5.76064e-01_jprb, 5.84344e-01_jprb/)
      kao_mo3(:, 8) = (/ &
     & 3.00557e-01_jprb, 3.03974e-01_jprb, 3.07430e-01_jprb, 3.10925e-01_jprb, 3.14460e-01_jprb, &
     & 3.18035e-01_jprb, 3.21651e-01_jprb, 3.25307e-01_jprb, 3.29006e-01_jprb, 3.32746e-01_jprb, &
     & 3.36529e-01_jprb, 3.40355e-01_jprb, 3.44224e-01_jprb, 3.48137e-01_jprb, 3.52095e-01_jprb, &
     & 3.56098e-01_jprb, 3.60146e-01_jprb, 3.64241e-01_jprb, 3.68381e-01_jprb/)
      kao_mo3(:, 9) = (/ &
     & 2.10042e-01_jprb, 2.12905e-01_jprb, 2.15806e-01_jprb, 2.18748e-01_jprb, 2.21729e-01_jprb, &
     & 2.24751e-01_jprb, 2.27814e-01_jprb, 2.30919e-01_jprb, 2.34066e-01_jprb, 2.37256e-01_jprb, &
     & 2.40489e-01_jprb, 2.43767e-01_jprb, 2.47089e-01_jprb, 2.50457e-01_jprb, 2.53870e-01_jprb, &
     & 2.57330e-01_jprb, 2.60837e-01_jprb, 2.64392e-01_jprb, 2.67996e-01_jprb/)
      kao_mo3(:,10) = (/ &
     & 2.09288e-01_jprb, 2.11759e-01_jprb, 2.14259e-01_jprb, 2.16789e-01_jprb, 2.19349e-01_jprb, &
     & 2.21939e-01_jprb, 2.24559e-01_jprb, 2.27210e-01_jprb, 2.29893e-01_jprb, 2.32607e-01_jprb, &
     & 2.35354e-01_jprb, 2.38133e-01_jprb, 2.40944e-01_jprb, 2.43789e-01_jprb, 2.46667e-01_jprb, &
     & 2.49580e-01_jprb, 2.52527e-01_jprb, 2.55508e-01_jprb, 2.58525e-01_jprb/)
      kao_mo3(:,11) = (/ &
     & 2.28947e-01_jprb, 2.30609e-01_jprb, 2.32283e-01_jprb, 2.33969e-01_jprb, 2.35667e-01_jprb, &
     & 2.37378e-01_jprb, 2.39101e-01_jprb, 2.40836e-01_jprb, 2.42584e-01_jprb, 2.44345e-01_jprb, &
     & 2.46118e-01_jprb, 2.47905e-01_jprb, 2.49704e-01_jprb, 2.51516e-01_jprb, 2.53342e-01_jprb, &
     & 2.55181e-01_jprb, 2.57033e-01_jprb, 2.58899e-01_jprb, 2.60778e-01_jprb/)
      kao_mo3(:,12) = (/ &
     & 2.57263e-01_jprb, 2.58272e-01_jprb, 2.59285e-01_jprb, 2.60302e-01_jprb, 2.61323e-01_jprb, &
     & 2.62347e-01_jprb, 2.63376e-01_jprb, 2.64409e-01_jprb, 2.65446e-01_jprb, 2.66487e-01_jprb, &
     & 2.67532e-01_jprb, 2.68581e-01_jprb, 2.69635e-01_jprb, 2.70692e-01_jprb, 2.71753e-01_jprb, &
     & 2.72819e-01_jprb, 2.73889e-01_jprb, 2.74963e-01_jprb, 2.76042e-01_jprb/)
      kao_mo3(:,13) = (/ &
     & 2.43322e-01_jprb, 2.45918e-01_jprb, 2.48541e-01_jprb, 2.51192e-01_jprb, 2.53872e-01_jprb, &
     & 2.56580e-01_jprb, 2.59317e-01_jprb, 2.62083e-01_jprb, 2.64879e-01_jprb, 2.67704e-01_jprb, &
     & 2.70560e-01_jprb, 2.73446e-01_jprb, 2.76363e-01_jprb, 2.79311e-01_jprb, 2.82290e-01_jprb, &
     & 2.85302e-01_jprb, 2.88345e-01_jprb, 2.91421e-01_jprb, 2.94529e-01_jprb/)
      kao_mo3(:,14) = (/ &
     & 2.10568e-01_jprb, 2.16529e-01_jprb, 2.22660e-01_jprb, 2.28964e-01_jprb, 2.35446e-01_jprb, &
     & 2.42113e-01_jprb, 2.48967e-01_jprb, 2.56016e-01_jprb, 2.63265e-01_jprb, 2.70719e-01_jprb, &
     & 2.78383e-01_jprb, 2.86265e-01_jprb, 2.94370e-01_jprb, 3.02704e-01_jprb, 3.11275e-01_jprb, &
     & 3.20088e-01_jprb, 3.29150e-01_jprb, 3.38470e-01_jprb, 3.48052e-01_jprb/)
      kao_mo3(:,15) = (/ &
     & 2.60406e-02_jprb, 2.78779e-02_jprb, 2.98448e-02_jprb, 3.19505e-02_jprb, 3.42048e-02_jprb, &
     & 3.66181e-02_jprb, 3.92017e-02_jprb, 4.19675e-02_jprb, 4.49285e-02_jprb, 4.80985e-02_jprb, &
     & 5.14920e-02_jprb, 5.51250e-02_jprb, 5.90143e-02_jprb, 6.31781e-02_jprb, 6.76356e-02_jprb, &
     & 7.24076e-02_jprb, 7.75163e-02_jprb, 8.29854e-02_jprb, 8.88404e-02_jprb/)
      kao_mo3(:,16) = (/ &
     & 2.31483e-02_jprb, 2.46840e-02_jprb, 2.63217e-02_jprb, 2.80681e-02_jprb, 2.99302e-02_jprb, &
     & 3.19160e-02_jprb, 3.40335e-02_jprb, 3.62914e-02_jprb, 3.86992e-02_jprb, 4.12668e-02_jprb, &
     & 4.40046e-02_jprb, 4.69242e-02_jprb, 5.00374e-02_jprb, 5.33571e-02_jprb, 5.68971e-02_jprb, &
     & 6.06720e-02_jprb, 6.46974e-02_jprb, 6.89897e-02_jprb, 7.35669e-02_jprb/)

      kao_mn2o(:, 1) = (/ &
     & 3.02276e-02_jprb, 3.10321e-02_jprb, 3.18580e-02_jprb, 3.27059e-02_jprb, 3.35764e-02_jprb, &
     & 3.44700e-02_jprb, 3.53875e-02_jprb, 3.63293e-02_jprb, 3.72962e-02_jprb, 3.82889e-02_jprb, &
     & 3.93079e-02_jprb, 4.03541e-02_jprb, 4.14281e-02_jprb, 4.25307e-02_jprb, 4.36627e-02_jprb, &
     & 4.48248e-02_jprb, 4.60178e-02_jprb, 4.72425e-02_jprb, 4.84999e-02_jprb/)
      kao_mn2o(:, 2) = (/ &
     & 6.10132e-02_jprb, 6.17435e-02_jprb, 6.24825e-02_jprb, 6.32304e-02_jprb, 6.39872e-02_jprb, &
     & 6.47531e-02_jprb, 6.55281e-02_jprb, 6.63124e-02_jprb, 6.71061e-02_jprb, 6.79093e-02_jprb, &
     & 6.87221e-02_jprb, 6.95446e-02_jprb, 7.03770e-02_jprb, 7.12194e-02_jprb, 7.20718e-02_jprb, &
     & 7.29344e-02_jprb, 7.38074e-02_jprb, 7.46908e-02_jprb, 7.55848e-02_jprb/)
      kao_mn2o(:, 3) = (/ &
     & 1.04479e-01_jprb, 1.05566e-01_jprb, 1.06664e-01_jprb, 1.07774e-01_jprb, 1.08895e-01_jprb, &
     & 1.10028e-01_jprb, 1.11173e-01_jprb, 1.12329e-01_jprb, 1.13498e-01_jprb, 1.14679e-01_jprb, &
     & 1.15872e-01_jprb, 1.17077e-01_jprb, 1.18295e-01_jprb, 1.19526e-01_jprb, 1.20770e-01_jprb, &
     & 1.22026e-01_jprb, 1.23296e-01_jprb, 1.24578e-01_jprb, 1.25875e-01_jprb/)
      kao_mn2o(:, 4) = (/ &
     & 2.07260e-01_jprb, 2.08126e-01_jprb, 2.08996e-01_jprb, 2.09869e-01_jprb, 2.10746e-01_jprb, &
     & 2.11627e-01_jprb, 2.12511e-01_jprb, 2.13399e-01_jprb, 2.14291e-01_jprb, 2.15187e-01_jprb, &
     & 2.16086e-01_jprb, 2.16989e-01_jprb, 2.17896e-01_jprb, 2.18807e-01_jprb, 2.19721e-01_jprb, &
     & 2.20640e-01_jprb, 2.21562e-01_jprb, 2.22488e-01_jprb, 2.23418e-01_jprb/)
      kao_mn2o(:, 5) = (/ &
     & 3.71566e-01_jprb, 3.71353e-01_jprb, 3.71141e-01_jprb, 3.70928e-01_jprb, 3.70716e-01_jprb, &
     & 3.70504e-01_jprb, 3.70292e-01_jprb, 3.70080e-01_jprb, 3.69869e-01_jprb, 3.69657e-01_jprb, &
     & 3.69446e-01_jprb, 3.69234e-01_jprb, 3.69023e-01_jprb, 3.68812e-01_jprb, 3.68601e-01_jprb, &
     & 3.68390e-01_jprb, 3.68179e-01_jprb, 3.67969e-01_jprb, 3.67758e-01_jprb/)
      kao_mn2o(:, 6) = (/ &
     & 5.28092e-01_jprb, 5.27262e-01_jprb, 5.26433e-01_jprb, 5.25605e-01_jprb, 5.24779e-01_jprb, &
     & 5.23954e-01_jprb, 5.23130e-01_jprb, 5.22307e-01_jprb, 5.21486e-01_jprb, 5.20666e-01_jprb, &
     & 5.19847e-01_jprb, 5.19030e-01_jprb, 5.18214e-01_jprb, 5.17399e-01_jprb, 5.16586e-01_jprb, &
     & 5.15773e-01_jprb, 5.14962e-01_jprb, 5.14153e-01_jprb, 5.13344e-01_jprb/)
      kao_mn2o(:, 7) = (/ &
     & 3.88140e-01_jprb, 3.87956e-01_jprb, 3.87773e-01_jprb, 3.87590e-01_jprb, 3.87407e-01_jprb, &
     & 3.87224e-01_jprb, 3.87041e-01_jprb, 3.86858e-01_jprb, 3.86675e-01_jprb, 3.86492e-01_jprb, &
     & 3.86310e-01_jprb, 3.86127e-01_jprb, 3.85945e-01_jprb, 3.85763e-01_jprb, 3.85580e-01_jprb, &
     & 3.85398e-01_jprb, 3.85216e-01_jprb, 3.85034e-01_jprb, 3.84852e-01_jprb/)
      kao_mn2o(:, 8) = (/ &
     & 3.12991e-01_jprb, 3.12246e-01_jprb, 3.11504e-01_jprb, 3.10763e-01_jprb, 3.10024e-01_jprb, &
     & 3.09287e-01_jprb, 3.08552e-01_jprb, 3.07818e-01_jprb, 3.07086e-01_jprb, 3.06356e-01_jprb, &
     & 3.05628e-01_jprb, 3.04901e-01_jprb, 3.04176e-01_jprb, 3.03453e-01_jprb, 3.02732e-01_jprb, &
     & 3.02012e-01_jprb, 3.01294e-01_jprb, 3.00577e-01_jprb, 2.99863e-01_jprb/)
      kao_mn2o(:, 9) = (/ &
     & 4.11761e-01_jprb, 4.11309e-01_jprb, 4.10858e-01_jprb, 4.10407e-01_jprb, 4.09957e-01_jprb, &
     & 4.09507e-01_jprb, 4.09057e-01_jprb, 4.08608e-01_jprb, 4.08160e-01_jprb, 4.07712e-01_jprb, &
     & 4.07265e-01_jprb, 4.06818e-01_jprb, 4.06371e-01_jprb, 4.05925e-01_jprb, 4.05480e-01_jprb, &
     & 4.05035e-01_jprb, 4.04590e-01_jprb, 4.04146e-01_jprb, 4.03703e-01_jprb/)
      kao_mn2o(:,10) = (/ &
     & 2.84648e-01_jprb, 2.87025e-01_jprb, 2.89421e-01_jprb, 2.91838e-01_jprb, 2.94275e-01_jprb, &
     & 2.96732e-01_jprb, 2.99210e-01_jprb, 3.01708e-01_jprb, 3.04227e-01_jprb, 3.06768e-01_jprb, &
     & 3.09329e-01_jprb, 3.11912e-01_jprb, 3.14517e-01_jprb, 3.17143e-01_jprb, 3.19791e-01_jprb, &
     & 3.22461e-01_jprb, 3.25154e-01_jprb, 3.27869e-01_jprb, 3.30606e-01_jprb/)
      kao_mn2o(:,11) = (/ &
     & 2.75090e-01_jprb, 2.79370e-01_jprb, 2.83716e-01_jprb, 2.88129e-01_jprb, 2.92611e-01_jprb, &
     & 2.97163e-01_jprb, 3.01786e-01_jprb, 3.06480e-01_jprb, 3.11248e-01_jprb, 3.16090e-01_jprb, &
     & 3.21007e-01_jprb, 3.26001e-01_jprb, 3.31072e-01_jprb, 3.36222e-01_jprb, 3.41452e-01_jprb, &
     & 3.46764e-01_jprb, 3.52158e-01_jprb, 3.57636e-01_jprb, 3.63200e-01_jprb/)
      kao_mn2o(:,12) = (/ &
     & 1.67753e-01_jprb, 1.71386e-01_jprb, 1.75098e-01_jprb, 1.78890e-01_jprb, 1.82765e-01_jprb, &
     & 1.86723e-01_jprb, 1.90767e-01_jprb, 1.94899e-01_jprb, 1.99120e-01_jprb, 2.03433e-01_jprb, &
     & 2.07839e-01_jprb, 2.12340e-01_jprb, 2.16939e-01_jprb, 2.21638e-01_jprb, 2.26438e-01_jprb, &
     & 2.31342e-01_jprb, 2.36353e-01_jprb, 2.41472e-01_jprb, 2.46701e-01_jprb/)
      kao_mn2o(:,13) = (/ &
     & 1.40543e-01_jprb, 1.42049e-01_jprb, 1.43571e-01_jprb, 1.45109e-01_jprb, 1.46663e-01_jprb, &
     & 1.48234e-01_jprb, 1.49822e-01_jprb, 1.51427e-01_jprb, 1.53049e-01_jprb, 1.54689e-01_jprb, &
     & 1.56346e-01_jprb, 1.58021e-01_jprb, 1.59713e-01_jprb, 1.61424e-01_jprb, 1.63153e-01_jprb, &
     & 1.64901e-01_jprb, 1.66668e-01_jprb, 1.68453e-01_jprb, 1.70258e-01_jprb/)
      kao_mn2o(:,14) = (/ &
     & 1.51530e-01_jprb, 1.50944e-01_jprb, 1.50360e-01_jprb, 1.49779e-01_jprb, 1.49199e-01_jprb, &
     & 1.48622e-01_jprb, 1.48047e-01_jprb, 1.47474e-01_jprb, 1.46903e-01_jprb, 1.46335e-01_jprb, &
     & 1.45769e-01_jprb, 1.45205e-01_jprb, 1.44643e-01_jprb, 1.44083e-01_jprb, 1.43526e-01_jprb, &
     & 1.42971e-01_jprb, 1.42418e-01_jprb, 1.41867e-01_jprb, 1.41318e-01_jprb/)
      kao_mn2o(:,15) = (/ &
     & 2.20492e-01_jprb, 2.16479e-01_jprb, 2.12539e-01_jprb, 2.08671e-01_jprb, 2.04873e-01_jprb, &
     & 2.01145e-01_jprb, 1.97484e-01_jprb, 1.93890e-01_jprb, 1.90361e-01_jprb, 1.86897e-01_jprb, &
     & 1.83495e-01_jprb, 1.80156e-01_jprb, 1.76877e-01_jprb, 1.73658e-01_jprb, 1.70497e-01_jprb, &
     & 1.67394e-01_jprb, 1.64348e-01_jprb, 1.61356e-01_jprb, 1.58420e-01_jprb/)
      kao_mn2o(:,16) = (/ &
     & 2.19848e-01_jprb, 2.15847e-01_jprb, 2.11919e-01_jprb, 2.08062e-01_jprb, 2.04275e-01_jprb, &
     & 2.00558e-01_jprb, 1.96908e-01_jprb, 1.93324e-01_jprb, 1.89806e-01_jprb, 1.86351e-01_jprb, &
     & 1.82960e-01_jprb, 1.79630e-01_jprb, 1.76361e-01_jprb, 1.73151e-01_jprb, 1.70000e-01_jprb, &
     & 1.66906e-01_jprb, 1.63868e-01_jprb, 1.60886e-01_jprb, 1.57958e-01_jprb/)

!     the array kbo_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kbo_mco2(:, 1) = (/ &
     & 4.74280e-08_jprb, 6.62724e-08_jprb, 9.26042e-08_jprb, 1.29398e-07_jprb, 1.80812e-07_jprb, &
     & 2.52653e-07_jprb, 3.53039e-07_jprb, 4.93310e-07_jprb, 6.89316e-07_jprb, 9.63198e-07_jprb, &
     & 1.34590e-06_jprb, 1.88067e-06_jprb, 2.62790e-06_jprb, 3.67204e-06_jprb, 5.13104e-06_jprb, &
     & 7.16974e-06_jprb, 1.00185e-05_jprb, 1.39991e-05_jprb, 1.95613e-05_jprb/)
      kbo_mco2(:, 2) = (/ &
     & 1.14872e-07_jprb, 1.63356e-07_jprb, 2.32304e-07_jprb, 3.30352e-07_jprb, 4.69783e-07_jprb, &
     & 6.68064e-07_jprb, 9.50033e-07_jprb, 1.35101e-06_jprb, 1.92123e-06_jprb, 2.73213e-06_jprb, &
     & 3.88527e-06_jprb, 5.52513e-06_jprb, 7.85711e-06_jprb, 1.11734e-05_jprb, 1.58893e-05_jprb, &
     & 2.25957e-05_jprb, 3.21326e-05_jprb, 4.56948e-05_jprb, 6.49811e-05_jprb/)
      kbo_mco2(:, 3) = (/ &
     & 3.30676e-07_jprb, 4.76313e-07_jprb, 6.86094e-07_jprb, 9.88267e-07_jprb, 1.42353e-06_jprb, &
     & 2.05048e-06_jprb, 2.95356e-06_jprb, 4.25439e-06_jprb, 6.12813e-06_jprb, 8.82711e-06_jprb, &
     & 1.27148e-05_jprb, 1.83147e-05_jprb, 2.63810e-05_jprb, 3.79998e-05_jprb, 5.47359e-05_jprb, &
     & 7.88430e-05_jprb, 1.13568e-04_jprb, 1.63585e-04_jprb, 2.35632e-04_jprb/)
      kbo_mco2(:, 4) = (/ &
     & 6.58642e-07_jprb, 9.52761e-07_jprb, 1.37822e-06_jprb, 1.99368e-06_jprb, 2.88396e-06_jprb, &
     & 4.17181e-06_jprb, 6.03475e-06_jprb, 8.72960e-06_jprb, 1.26279e-05_jprb, 1.82669e-05_jprb, &
     & 2.64241e-05_jprb, 3.82239e-05_jprb, 5.52929e-05_jprb, 7.99844e-05_jprb, 1.15702e-04_jprb, &
     & 1.67369e-04_jprb, 2.42109e-04_jprb, 3.50223e-04_jprb, 5.06617e-04_jprb/)
      kbo_mco2(:, 5) = (/ &
     & 1.26418e-06_jprb, 1.82095e-06_jprb, 2.62292e-06_jprb, 3.77810e-06_jprb, 5.44204e-06_jprb, &
     & 7.83881e-06_jprb, 1.12911e-05_jprb, 1.62640e-05_jprb, 2.34269e-05_jprb, 3.37445e-05_jprb, &
     & 4.86061e-05_jprb, 7.00131e-05_jprb, 1.00848e-04_jprb, 1.45263e-04_jprb, 2.09239e-04_jprb, &
     & 3.01392e-04_jprb, 4.34131e-04_jprb, 6.25329e-04_jprb, 9.00733e-04_jprb/)
      kbo_mco2(:, 6) = (/ &
     & 2.38529e-06_jprb, 3.43110e-06_jprb, 4.93545e-06_jprb, 7.09937e-06_jprb, 1.02120e-05_jprb, &
     & 1.46895e-05_jprb, 2.11300e-05_jprb, 3.03943e-05_jprb, 4.37205e-05_jprb, 6.28894e-05_jprb, &
     & 9.04630e-05_jprb, 1.30126e-04_jprb, 1.87179e-04_jprb, 2.69247e-04_jprb, 3.87296e-04_jprb, &
     & 5.57104e-04_jprb, 8.01364e-04_jprb, 1.15272e-03_jprb, 1.65812e-03_jprb/)
      kbo_mco2(:, 7) = (/ &
     & 5.41398e-06_jprb, 7.54295e-06_jprb, 1.05091e-05_jprb, 1.46417e-05_jprb, 2.03993e-05_jprb, &
     & 2.84211e-05_jprb, 3.95973e-05_jprb, 5.51683e-05_jprb, 7.68626e-05_jprb, 1.07088e-04_jprb, &
     & 1.49199e-04_jprb, 2.07869e-04_jprb, 2.89610e-04_jprb, 4.03496e-04_jprb, 5.62165e-04_jprb, &
     & 7.83229e-04_jprb, 1.09122e-03_jprb, 1.52033e-03_jprb, 2.11818e-03_jprb/)
      kbo_mco2(:, 8) = (/ &
     & 1.09995e-05_jprb, 1.54018e-05_jprb, 2.15660e-05_jprb, 3.01973e-05_jprb, 4.22831e-05_jprb, &
     & 5.92059e-05_jprb, 8.29017e-05_jprb, 1.16081e-04_jprb, 1.62540e-04_jprb, 2.27592e-04_jprb, &
     & 3.18681e-04_jprb, 4.46226e-04_jprb, 6.24817e-04_jprb, 8.74886e-04_jprb, 1.22504e-03_jprb, &
     & 1.71533e-03_jprb, 2.40185e-03_jprb, 3.36313e-03_jprb, 4.70915e-03_jprb/)
      kbo_mco2(:, 9) = (/ &
     & 3.29051e-05_jprb, 4.59996e-05_jprb, 6.43050e-05_jprb, 8.98950e-05_jprb, 1.25668e-04_jprb, &
     & 1.75678e-04_jprb, 2.45588e-04_jprb, 3.43319e-04_jprb, 4.79942e-04_jprb, 6.70933e-04_jprb, &
     & 9.37930e-04_jprb, 1.31118e-03_jprb, 1.83295e-03_jprb, 2.56237e-03_jprb, 3.58206e-03_jprb, &
     & 5.00753e-03_jprb, 7.00027e-03_jprb, 9.78599e-03_jprb, 1.36803e-02_jprb/)
      kbo_mco2(:,10) = (/ &
     & 1.95126e-05_jprb, 2.65944e-05_jprb, 3.62463e-05_jprb, 4.94013e-05_jprb, 6.73305e-05_jprb, &
     & 9.17669e-05_jprb, 1.25072e-04_jprb, 1.70465e-04_jprb, 2.32332e-04_jprb, 3.16652e-04_jprb, &
     & 4.31575e-04_jprb, 5.88208e-04_jprb, 8.01687e-04_jprb, 1.09264e-03_jprb, 1.48920e-03_jprb, &
     & 2.02968e-03_jprb, 2.76631e-03_jprb, 3.77029e-03_jprb, 5.13865e-03_jprb/)
      kbo_mco2(:,11) = (/ &
     & 8.67271e-05_jprb, 1.19228e-04_jprb, 1.63908e-04_jprb, 2.25332e-04_jprb, 3.09774e-04_jprb, &
     & 4.25860e-04_jprb, 5.85450e-04_jprb, 8.04845e-04_jprb, 1.10646e-03_jprb, 1.52110e-03_jprb, &
     & 2.09112e-03_jprb, 2.87476e-03_jprb, 3.95207e-03_jprb, 5.43309e-03_jprb, 7.46911e-03_jprb, &
     & 1.02681e-02_jprb, 1.41161e-02_jprb, 1.94060e-02_jprb, 2.66783e-02_jprb/)
      kbo_mco2(:,12) = (/ &
     & 3.79194e-07_jprb, 5.51419e-07_jprb, 8.01866e-07_jprb, 1.16606e-06_jprb, 1.69567e-06_jprb, &
     & 2.46582e-06_jprb, 3.58577e-06_jprb, 5.21438e-06_jprb, 7.58268e-06_jprb, 1.10266e-05_jprb, &
     & 1.60348e-05_jprb, 2.33176e-05_jprb, 3.39081e-05_jprb, 4.93087e-05_jprb, 7.17040e-05_jprb, &
     & 1.04271e-04_jprb, 1.51630e-04_jprb, 2.20498e-04_jprb, 3.20644e-04_jprb/)
      kbo_mco2(:,13) = (/ &
     & 1.72555e-07_jprb, 2.29952e-07_jprb, 3.06441e-07_jprb, 4.08373e-07_jprb, 5.44209e-07_jprb, &
     & 7.25229e-07_jprb, 9.66461e-07_jprb, 1.28793e-06_jprb, 1.71634e-06_jprb, 2.28724e-06_jprb, &
     & 3.04805e-06_jprb, 4.06192e-06_jprb, 5.41303e-06_jprb, 7.21356e-06_jprb, 9.61299e-06_jprb, &
     & 1.28106e-05_jprb, 1.70717e-05_jprb, 2.27503e-05_jprb, 3.03177e-05_jprb/)
      kbo_mco2(:,14) = (/ &
     & 7.42245e-09_jprb, 7.17780e-09_jprb, 6.94122e-09_jprb, 6.71243e-09_jprb, 6.49118e-09_jprb, &
     & 6.27723e-09_jprb, 6.07032e-09_jprb, 5.87024e-09_jprb, 5.67675e-09_jprb, 5.48964e-09_jprb, &
     & 5.30870e-09_jprb, 5.13372e-09_jprb, 4.96451e-09_jprb, 4.80087e-09_jprb, 4.64263e-09_jprb, &
     & 4.48961e-09_jprb, 4.34163e-09_jprb, 4.19852e-09_jprb, 4.06014e-09_jprb/)
      kbo_mco2(:,15) = (/ &
     & 7.41847e-09_jprb, 7.17332e-09_jprb, 6.93627e-09_jprb, 6.70705e-09_jprb, 6.48541e-09_jprb, &
     & 6.27109e-09_jprb, 6.06386e-09_jprb, 5.86347e-09_jprb, 5.66970e-09_jprb, 5.48234e-09_jprb, &
     & 5.30117e-09_jprb, 5.12599e-09_jprb, 4.95659e-09_jprb, 4.79280e-09_jprb, 4.63441e-09_jprb, &
     & 4.48126e-09_jprb, 4.33317e-09_jprb, 4.18998e-09_jprb, 4.05152e-09_jprb/)
      kbo_mco2(:,16) = (/ &
     & 7.42855e-09_jprb, 7.18278e-09_jprb, 6.94513e-09_jprb, 6.71535e-09_jprb, 6.49317e-09_jprb, &
     & 6.27834e-09_jprb, 6.07062e-09_jprb, 5.86977e-09_jprb, 5.67557e-09_jprb, 5.48779e-09_jprb, &
     & 5.30622e-09_jprb, 5.13066e-09_jprb, 4.96091e-09_jprb, 4.79678e-09_jprb, 4.63808e-09_jprb, &
     & 4.48462e-09_jprb, 4.33625e-09_jprb, 4.19278e-09_jprb, 4.05406e-09_jprb/)

      kbo_mn2o(:, 1) = (/ &
     & 2.49055e-04_jprb, 2.53574e-04_jprb, 2.58175e-04_jprb, 2.62860e-04_jprb, 2.67629e-04_jprb, &
     & 2.72485e-04_jprb, 2.77429e-04_jprb, 2.82463e-04_jprb, 2.87588e-04_jprb, 2.92806e-04_jprb, &
     & 2.98119e-04_jprb, 3.03528e-04_jprb, 3.09036e-04_jprb, 3.14643e-04_jprb, 3.20352e-04_jprb, &
     & 3.26165e-04_jprb, 3.32083e-04_jprb, 3.38109e-04_jprb, 3.44243e-04_jprb/)
      kbo_mn2o(:, 2) = (/ &
     & 3.79251e-04_jprb, 4.04353e-04_jprb, 4.31117e-04_jprb, 4.59652e-04_jprb, 4.90075e-04_jprb, &
     & 5.22513e-04_jprb, 5.57097e-04_jprb, 5.93970e-04_jprb, 6.33284e-04_jprb, 6.75200e-04_jprb, &
     & 7.19890e-04_jprb, 7.67539e-04_jprb, 8.18340e-04_jprb, 8.72505e-04_jprb, 9.30255e-04_jprb, &
     & 9.91827e-04_jprb, 1.05747e-03_jprb, 1.12747e-03_jprb, 1.20209e-03_jprb/)
      kbo_mn2o(:, 3) = (/ &
     & 7.61140e-04_jprb, 8.36483e-04_jprb, 9.19284e-04_jprb, 1.01028e-03_jprb, 1.11029e-03_jprb, &
     & 1.22019e-03_jprb, 1.34098e-03_jprb, 1.47372e-03_jprb, 1.61959e-03_jprb, 1.77991e-03_jprb, &
     & 1.95610e-03_jprb, 2.14973e-03_jprb, 2.36253e-03_jprb, 2.59639e-03_jprb, 2.85340e-03_jprb, &
     & 3.13585e-03_jprb, 3.44626e-03_jprb, 3.78740e-03_jprb, 4.16230e-03_jprb/)
      kbo_mn2o(:, 4) = (/ &
     & 2.01074e-03_jprb, 2.26915e-03_jprb, 2.56077e-03_jprb, 2.88987e-03_jprb, 3.26126e-03_jprb, &
     & 3.68038e-03_jprb, 4.15337e-03_jprb, 4.68714e-03_jprb, 5.28951e-03_jprb, 5.96929e-03_jprb, &
     & 6.73643e-03_jprb, 7.60217e-03_jprb, 8.57916e-03_jprb, 9.68172e-03_jprb, 1.09260e-02_jprb, &
     & 1.23301e-02_jprb, 1.39147e-02_jprb, 1.57030e-02_jprb, 1.77211e-02_jprb/)
      kbo_mn2o(:, 5) = (/ &
     & 7.43302e-03_jprb, 8.32582e-03_jprb, 9.32585e-03_jprb, 1.04460e-02_jprb, 1.17007e-02_jprb, &
     & 1.31061e-02_jprb, 1.46803e-02_jprb, 1.64436e-02_jprb, 1.84186e-02_jprb, 2.06309e-02_jprb, &
     & 2.31090e-02_jprb, 2.58846e-02_jprb, 2.89937e-02_jprb, 3.24762e-02_jprb, 3.63769e-02_jprb, &
     & 4.07463e-02_jprb, 4.56404e-02_jprb, 5.11223e-02_jprb, 5.72627e-02_jprb/)
      kbo_mn2o(:, 6) = (/ &
     & 2.71911e-02_jprb, 2.94258e-02_jprb, 3.18441e-02_jprb, 3.44612e-02_jprb, 3.72933e-02_jprb, &
     & 4.03582e-02_jprb, 4.36750e-02_jprb, 4.72644e-02_jprb, 5.11487e-02_jprb, 5.53523e-02_jprb, &
     & 5.99014e-02_jprb, 6.48243e-02_jprb, 7.01518e-02_jprb, 7.59172e-02_jprb, 8.21563e-02_jprb, &
     & 8.89082e-02_jprb, 9.62150e-02_jprb, 1.04122e-01_jprb, 1.12679e-01_jprb/)
      kbo_mn2o(:, 7) = (/ &
     & 1.63331e-01_jprb, 1.80469e-01_jprb, 1.99406e-01_jprb, 2.20330e-01_jprb, 2.43449e-01_jprb, &
     & 2.68995e-01_jprb, 2.97221e-01_jprb, 3.28408e-01_jprb, 3.62869e-01_jprb, 4.00945e-01_jprb, &
     & 4.43017e-01_jprb, 4.89503e-01_jprb, 5.40867e-01_jprb, 5.97621e-01_jprb, 6.60330e-01_jprb, &
     & 7.29619e-01_jprb, 8.06179e-01_jprb, 8.90772e-01_jprb, 9.84242e-01_jprb/)
      kbo_mn2o(:, 8) = (/ &
     & 1.32648e+00_jprb, 1.33515e+00_jprb, 1.34387e+00_jprb, 1.35265e+00_jprb, 1.36149e+00_jprb, &
     & 1.37038e+00_jprb, 1.37933e+00_jprb, 1.38835e+00_jprb, 1.39742e+00_jprb, 1.40655e+00_jprb, &
     & 1.41574e+00_jprb, 1.42499e+00_jprb, 1.43429e+00_jprb, 1.44367e+00_jprb, 1.45310e+00_jprb, &
     & 1.46259e+00_jprb, 1.47215e+00_jprb, 1.48176e+00_jprb, 1.49144e+00_jprb/)
      kbo_mn2o(:, 9) = (/ &
     & 3.12620e+00_jprb, 3.03118e+00_jprb, 2.93905e+00_jprb, 2.84972e+00_jprb, 2.76310e+00_jprb, &
     & 2.67911e+00_jprb, 2.59768e+00_jprb, 2.51873e+00_jprb, 2.44217e+00_jprb, 2.36794e+00_jprb, &
     & 2.29596e+00_jprb, 2.22618e+00_jprb, 2.15851e+00_jprb, 2.09290e+00_jprb, 2.02929e+00_jprb, &
     & 1.96761e+00_jprb, 1.90780e+00_jprb, 1.84982e+00_jprb, 1.79359e+00_jprb/)
      kbo_mn2o(:,10) = (/ &
     & 1.60677e-02_jprb, 1.82485e-02_jprb, 2.07254e-02_jprb, 2.35384e-02_jprb, 2.67332e-02_jprb, &
     & 3.03617e-02_jprb, 3.44827e-02_jprb, 3.91629e-02_jprb, 4.44785e-02_jprb, 5.05154e-02_jprb, &
     & 5.73718e-02_jprb, 6.51589e-02_jprb, 7.40027e-02_jprb, 8.40470e-02_jprb, 9.54546e-02_jprb, &
     & 1.08411e-01_jprb, 1.23125e-01_jprb, 1.39836e-01_jprb, 1.58816e-01_jprb/)
      kbo_mn2o(:,11) = (/ &
     & 1.55287e-02_jprb, 1.78265e-02_jprb, 2.04642e-02_jprb, 2.34922e-02_jprb, 2.69683e-02_jprb, &
     & 3.09588e-02_jprb, 3.55397e-02_jprb, 4.07984e-02_jprb, 4.68352e-02_jprb, 5.37653e-02_jprb, &
     & 6.17208e-02_jprb, 7.08535e-02_jprb, 8.13375e-02_jprb, 9.33728e-02_jprb, 1.07189e-01_jprb, &
     & 1.23049e-01_jprb, 1.41257e-01_jprb, 1.62158e-01_jprb, 1.86152e-01_jprb/)
      kbo_mn2o(:,12) = (/ &
     & 7.13719e-03_jprb, 8.18879e-03_jprb, 9.39535e-03_jprb, 1.07797e-02_jprb, 1.23680e-02_jprb, &
     & 1.41903e-02_jprb, 1.62811e-02_jprb, 1.86800e-02_jprb, 2.14324e-02_jprb, 2.45902e-02_jprb, &
     & 2.82134e-02_jprb, 3.23704e-02_jprb, 3.71400e-02_jprb, 4.26122e-02_jprb, 4.88908e-02_jprb, &
     & 5.60945e-02_jprb, 6.43596e-02_jprb, 7.38424e-02_jprb, 8.47224e-02_jprb/)
      kbo_mn2o(:,13) = (/ &
     & 9.28813e-03_jprb, 1.06108e-02_jprb, 1.21218e-02_jprb, 1.38480e-02_jprb, 1.58199e-02_jprb, &
     & 1.80727e-02_jprb, 2.06463e-02_jprb, 2.35864e-02_jprb, 2.69452e-02_jprb, 3.07822e-02_jprb, &
     & 3.51657e-02_jprb, 4.01734e-02_jprb, 4.58941e-02_jprb, 5.24296e-02_jprb, 5.98956e-02_jprb, &
     & 6.84249e-02_jprb, 7.81688e-02_jprb, 8.93002e-02_jprb, 1.02017e-01_jprb/)
      kbo_mn2o(:,14) = (/ &
     & 2.17205e-02_jprb, 2.51661e-02_jprb, 2.91581e-02_jprb, 3.37835e-02_jprb, 3.91425e-02_jprb, &
     & 4.53517e-02_jprb, 5.25458e-02_jprb, 6.08811e-02_jprb, 7.05387e-02_jprb, 8.17282e-02_jprb, &
     & 9.46927e-02_jprb, 1.09714e-01_jprb, 1.27118e-01_jprb, 1.47282e-01_jprb, 1.70645e-01_jprb, &
     & 1.97715e-01_jprb, 2.29078e-01_jprb, 2.65417e-01_jprb, 3.07520e-01_jprb/)
      kbo_mn2o(:,15) = (/ &
     & 4.89156e-02_jprb, 5.70504e-02_jprb, 6.65379e-02_jprb, 7.76033e-02_jprb, 9.05089e-02_jprb, &
     & 1.05561e-01_jprb, 1.23116e-01_jprb, 1.43590e-01_jprb, 1.67469e-01_jprb, 1.95320e-01_jprb, &
     & 2.27802e-01_jprb, 2.65686e-01_jprb, 3.09869e-01_jprb, 3.61401e-01_jprb, 4.21503e-01_jprb, &
     & 4.91600e-01_jprb, 5.73354e-01_jprb, 6.68703e-01_jprb, 7.79910e-01_jprb/)
      kbo_mn2o(:,16) = (/ &
     & 1.13156e-02_jprb, 1.46199e-02_jprb, 1.88891e-02_jprb, 2.44050e-02_jprb, 3.15316e-02_jprb, &
     & 4.07393e-02_jprb, 5.26358e-02_jprb, 6.80061e-02_jprb, 8.78649e-02_jprb, 1.13523e-01_jprb, &
     & 1.46673e-01_jprb, 1.89504e-01_jprb, 2.44841e-01_jprb, 3.16338e-01_jprb, 4.08713e-01_jprb, &
     & 5.28064e-01_jprb, 6.82266e-01_jprb, 8.81496e-01_jprb, 1.13891e+00_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &4.8166e-07_jprb,3.7500e-07_jprb,4.8978e-07_jprb,5.9624e-07_jprb,6.3742e-07_jprb,7.5551e-07_jprb, &
     &7.7706e-07_jprb,6.8681e-07_jprb,7.5212e-07_jprb,8.0956e-07_jprb,7.8117e-07_jprb,7.4835e-07_jprb, &
     &9.4118e-07_jprb,1.2585e-06_jprb,1.4976e-06_jprb,1.4976e-06_jprb/)
      forrefo(2,:) = (/ &
     &3.1320e-07_jprb,4.0764e-07_jprb,4.7468e-07_jprb,5.9976e-07_jprb,7.3324e-07_jprb,8.1488e-07_jprb, &
     &7.6442e-07_jprb,8.2007e-07_jprb,7.7721e-07_jprb,7.6377e-07_jprb,8.0327e-07_jprb,7.1881e-07_jprb, &
     &8.2148e-07_jprb,1.0203e-06_jprb,1.5033e-06_jprb,1.5032e-06_jprb/)
      forrefo(3,:) = (/ &
     &4.1831e-07_jprb,5.5043e-07_jprb,5.7783e-07_jprb,6.1294e-07_jprb,6.3396e-07_jprb,6.2292e-07_jprb, &
     &6.1719e-07_jprb,6.4183e-07_jprb,7.6180e-07_jprb,9.5477e-07_jprb,9.5901e-07_jprb,1.0207e-06_jprb, &
     &1.0387e-06_jprb,1.1305e-06_jprb,1.3602e-06_jprb,1.5063e-06_jprb/)
      forrefo(4,:) = (/ &
     &8.5878e-07_jprb,6.0921e-07_jprb,5.5773e-07_jprb,5.3374e-07_jprb,5.0495e-07_jprb,4.9844e-07_jprb, &
     &5.1536e-07_jprb,5.2908e-07_jprb,4.7977e-07_jprb,5.3177e-07_jprb,4.9266e-07_jprb,4.5403e-07_jprb, &
     &3.9695e-07_jprb,3.4792e-07_jprb,3.4912e-07_jprb,3.4102e-07_jprb/)

!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 3.16029e-02_jprb, 2.74633e-02_jprb, 2.38660e-02_jprb, 2.07399e-02_jprb, 1.80232e-02_jprb, &
     & 1.56624e-02_jprb, 1.36108e-02_jprb, 1.18280e-02_jprb, 1.02787e-02_jprb, 8.93231e-03_jprb/)
      selfrefo(:, 2) = (/ &
     & 3.10422e-02_jprb, 2.71312e-02_jprb, 2.37130e-02_jprb, 2.07254e-02_jprb, 1.81142e-02_jprb, &
     & 1.58320e-02_jprb, 1.38374e-02_jprb, 1.20940e-02_jprb, 1.05703e-02_jprb, 9.23854e-03_jprb/)
      selfrefo(:, 3) = (/ &
     & 3.08657e-02_jprb, 2.69431e-02_jprb, 2.35190e-02_jprb, 2.05301e-02_jprb, 1.79210e-02_jprb, &
     & 1.56435e-02_jprb, 1.36554e-02_jprb, 1.19200e-02_jprb, 1.04051e-02_jprb, 9.08279e-03_jprb/)
      selfrefo(:, 4) = (/ &
     & 3.02668e-02_jprb, 2.64686e-02_jprb, 2.31470e-02_jprb, 2.02422e-02_jprb, 1.77020e-02_jprb, &
     & 1.54806e-02_jprb, 1.35379e-02_jprb, 1.18390e-02_jprb, 1.03533e-02_jprb, 9.05406e-03_jprb/)
      selfrefo(:, 5) = (/ &
     & 2.98317e-02_jprb, 2.61491e-02_jprb, 2.29210e-02_jprb, 2.00914e-02_jprb, 1.76112e-02_jprb, &
     & 1.54371e-02_jprb, 1.35314e-02_jprb, 1.18610e-02_jprb, 1.03968e-02_jprb, 9.11332e-03_jprb/)
      selfrefo(:, 6) = (/ &
     & 2.95545e-02_jprb, 2.59083e-02_jprb, 2.27120e-02_jprb, 1.99100e-02_jprb, 1.74537e-02_jprb, &
     & 1.53004e-02_jprb, 1.34128e-02_jprb, 1.17580e-02_jprb, 1.03074e-02_jprb, 9.03576e-03_jprb/)
      selfrefo(:, 7) = (/ &
     & 2.97352e-02_jprb, 2.60320e-02_jprb, 2.27900e-02_jprb, 1.99517e-02_jprb, 1.74670e-02_jprb, &
     & 1.52916e-02_jprb, 1.33872e-02_jprb, 1.17200e-02_jprb, 1.02604e-02_jprb, 8.98258e-03_jprb/)
      selfrefo(:, 8) = (/ &
     & 2.96543e-02_jprb, 2.59760e-02_jprb, 2.27540e-02_jprb, 1.99316e-02_jprb, 1.74593e-02_jprb, &
     & 1.52937e-02_jprb, 1.33967e-02_jprb, 1.17350e-02_jprb, 1.02794e-02_jprb, 9.00437e-03_jprb/)
      selfrefo(:, 9) = (/ &
     & 2.97998e-02_jprb, 2.60786e-02_jprb, 2.28220e-02_jprb, 1.99721e-02_jprb, 1.74781e-02_jprb, &
     & 1.52955e-02_jprb, 1.33855e-02_jprb, 1.17140e-02_jprb, 1.02512e-02_jprb, 8.97110e-03_jprb/)
      selfrefo(:,10) = (/ &
     & 2.98826e-02_jprb, 2.61096e-02_jprb, 2.28130e-02_jprb, 1.99326e-02_jprb, 1.74159e-02_jprb, &
     & 1.52170e-02_jprb, 1.32957e-02_jprb, 1.16170e-02_jprb, 1.01502e-02_jprb, 8.86867e-03_jprb/)
      selfrefo(:,11) = (/ &
     & 2.94710e-02_jprb, 2.58147e-02_jprb, 2.26120e-02_jprb, 1.98066e-02_jprb, 1.73493e-02_jprb, &
     & 1.51969e-02_jprb, 1.33115e-02_jprb, 1.16600e-02_jprb, 1.02134e-02_jprb, 8.94628e-03_jprb/)
      selfrefo(:,12) = (/ &
     & 2.96297e-02_jprb, 2.59544e-02_jprb, 2.27350e-02_jprb, 1.99149e-02_jprb, 1.74446e-02_jprb, &
     & 1.52808e-02_jprb, 1.33853e-02_jprb, 1.17250e-02_jprb, 1.02706e-02_jprb, 8.99663e-03_jprb/)
      selfrefo(:,13) = (/ &
     & 2.96272e-02_jprb, 2.59013e-02_jprb, 2.26440e-02_jprb, 1.97963e-02_jprb, 1.73067e-02_jprb, &
     & 1.51302e-02_jprb, 1.32275e-02_jprb, 1.15640e-02_jprb, 1.01097e-02_jprb, 8.83833e-03_jprb/)
      selfrefo(:,14) = (/ &
     & 2.89906e-02_jprb, 2.53971e-02_jprb, 2.22490e-02_jprb, 1.94911e-02_jprb, 1.70751e-02_jprb, &
     & 1.49585e-02_jprb, 1.31044e-02_jprb, 1.14800e-02_jprb, 1.00570e-02_jprb, 8.81038e-03_jprb/)
      selfrefo(:,15) = (/ &
     & 2.80884e-02_jprb, 2.46987e-02_jprb, 2.17180e-02_jprb, 1.90970e-02_jprb, 1.67924e-02_jprb, &
     & 1.47659e-02_jprb, 1.29839e-02_jprb, 1.14170e-02_jprb, 1.00392e-02_jprb, 8.82765e-03_jprb/)
      selfrefo(:,16) = (/ &
     & 2.80884e-02_jprb, 2.46987e-02_jprb, 2.17180e-02_jprb, 1.90970e-02_jprb, 1.67924e-02_jprb, &
     & 1.47659e-02_jprb, 1.29839e-02_jprb, 1.14170e-02_jprb, 1.00392e-02_jprb, 8.82765e-03_jprb/)


if (lhook) call dr_hook('rrtm_kgb8',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb8:error reading file radrrtm")

end subroutine rrtm_kgb8
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

