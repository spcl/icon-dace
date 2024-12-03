! # 1 "ifsrrtm/rrtm_kgb1.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb1.f90"
! this file has been modified for the use in icon

subroutine rrtm_kgb1(directory)

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 1:  10-250 cm-1 (low - h2o; high - h2o)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo may 2013 update to rrtmg v4.85
!     band 1:  10-350 cm-1
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook
use yomlun    ,only : nulrad
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad
use yommp0    , only : nproc, myproc

use yoerrto1 , only : kao     ,kbo     ,kao_d,kbo_d,selfrefo   ,fracrefao ,&
 & fracrefbo  ,forrefo, kao_mn2, kbo_mn2

!     ------------------------------------------------------------------

implicit none

character(len=*), intent(in) :: directory

!character(len = 80) :: clzzz
character(len = 255) :: clf1
real(kind=jprb) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 36 "ifsrrtm/rrtm_kgb1.f90" 2

if (lhook) call dr_hook('rrtm_kgb1',0,zhook_handle)

if( myproc==1 )then
  !call getenv("data",clzzz)
  !if(clzzz /= " ") then
  !  clf1=trim(clzzz) // "/radrrtm"
  clf1 = directory // "/radrrtm"

  open(nulrad,file=trim(clf1),convert="big_endian",form="unformatted",action="read",err=1000)




  !else
  !  open(nulrad,file='radrrtm',form="unformatted",action="read",err=1000)
  !endif
  read(nulrad,err=1001) kao_d,kbo_d
 ! convert the data into model actual precision.
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb1:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb1:')
endif

! planck fraction mapping level: p = 212.7250 mbar, t = 223.06 k
fracrefao(:) = (/ &
 & 2.1227e-01_jprb,1.8897e-01_jprb,1.3934e-01_jprb,1.1557e-01_jprb,9.5282e-02_jprb,8.3359e-02_jprb, &
 & 6.5333e-02_jprb,5.2016e-02_jprb,3.4272e-02_jprb,4.0257e-03_jprb,3.1857e-03_jprb,2.6014e-03_jprb, &
 & 1.9141e-03_jprb,1.2612e-03_jprb,5.3169e-04_jprb,7.6476e-05_jprb/)

! planck fraction mapping level: p = 212.7250 mbar, t = 223.06 k
! these planck fractions were calculated using lower atmosphere
! parameters.
fracrefbo(:) = (/ &
 & 2.1227e-01_jprb,1.8897e-01_jprb,1.3934e-01_jprb,1.1557e-01_jprb,9.5282e-02_jprb,8.3359e-02_jprb, &
 & 6.5333e-02_jprb,5.2016e-02_jprb,3.4272e-02_jprb,4.0257e-03_jprb,3.1857e-03_jprb,2.6014e-03_jprb, &
 & 1.9141e-03_jprb,1.2612e-03_jprb,5.3169e-04_jprb,7.6476e-05_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     & 3.6742e-02_jprb,1.0664e-01_jprb,2.6132e-01_jprb,2.7906e-01_jprb,2.8151e-01_jprb,2.7465e-01_jprb, &
     & 2.8530e-01_jprb,2.9123e-01_jprb,3.0697e-01_jprb,3.1801e-01_jprb,3.2444e-01_jprb,2.7746e-01_jprb, &
     & 3.1994e-01_jprb,2.9750e-01_jprb,2.1226e-01_jprb,1.2847e-01_jprb/)
      forrefo(2,:) = (/ &
     & 4.0450e-02_jprb,1.1085e-01_jprb,2.9205e-01_jprb,3.1934e-01_jprb,3.1739e-01_jprb,3.1450e-01_jprb, &
     & 3.2797e-01_jprb,3.2223e-01_jprb,3.3099e-01_jprb,3.4800e-01_jprb,3.4046e-01_jprb,3.5700e-01_jprb, &
     & 3.8264e-01_jprb,3.6679e-01_jprb,3.3481e-01_jprb,3.2113e-01_jprb/)
      forrefo(3,:) = (/ &
     & 4.6952e-02_jprb,1.1999e-01_jprb,3.1473e-01_jprb,3.7015e-01_jprb,3.6913e-01_jprb,3.6352e-01_jprb, &
     & 3.7754e-01_jprb,3.7402e-01_jprb,3.7113e-01_jprb,3.7720e-01_jprb,3.8365e-01_jprb,4.0876e-01_jprb, &
     & 4.2968e-01_jprb,4.4186e-01_jprb,4.3468e-01_jprb,4.7083e-01_jprb/)
      forrefo(4,:) = (/ &
     & 7.0645e-02_jprb,1.6618e-01_jprb,2.8516e-01_jprb,3.1819e-01_jprb,3.0131e-01_jprb,2.9552e-01_jprb, &
     & 2.8972e-01_jprb,2.9348e-01_jprb,2.8668e-01_jprb,2.8483e-01_jprb,2.8130e-01_jprb,2.7757e-01_jprb, &
     & 2.9735e-01_jprb,3.1684e-01_jprb,3.0681e-01_jprb,3.6778e-01_jprb/)


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



     kao_mn2(:, 1) = (/ &
     & 5.12042e-08_jprb, 5.51239e-08_jprb, 5.93436e-08_jprb, 6.38863e-08_jprb, 6.87767e-08_jprb, &
     & 7.40415e-08_jprb, 7.97093e-08_jprb, 8.58110e-08_jprb, 9.23797e-08_jprb, 9.94513e-08_jprb, &
     & 1.07064e-07_jprb, 1.15260e-07_jprb, 1.24083e-07_jprb, 1.33581e-07_jprb, 1.43807e-07_jprb, &
     & 1.54815e-07_jprb, 1.66666e-07_jprb, 1.79424e-07_jprb, 1.93159e-07_jprb/)
      kao_mn2(:, 2) = (/ &
     & 2.30938e-07_jprb, 2.41696e-07_jprb, 2.52955e-07_jprb, 2.64738e-07_jprb, 2.77071e-07_jprb, &
     & 2.89978e-07_jprb, 3.03486e-07_jprb, 3.17623e-07_jprb, 3.32419e-07_jprb, 3.47904e-07_jprb, &
     & 3.64111e-07_jprb, 3.81072e-07_jprb, 3.98824e-07_jprb, 4.17402e-07_jprb, 4.36846e-07_jprb, &
     & 4.57196e-07_jprb, 4.78494e-07_jprb, 5.00784e-07_jprb, 5.24112e-07_jprb/)
      kao_mn2(:, 3) = (/ &
     & 6.70458e-07_jprb, 7.04274e-07_jprb, 7.39795e-07_jprb, 7.77109e-07_jprb, 8.16304e-07_jprb, &
     & 8.57476e-07_jprb, 9.00724e-07_jprb, 9.46154e-07_jprb, 9.93876e-07_jprb, 1.04400e-06_jprb, &
     & 1.09666e-06_jprb, 1.15197e-06_jprb, 1.21008e-06_jprb, 1.27111e-06_jprb, 1.33522e-06_jprb, &
     & 1.40256e-06_jprb, 1.47331e-06_jprb, 1.54761e-06_jprb, 1.62567e-06_jprb/)
      kao_mn2(:, 4) = (/ &
     & 1.84182e-06_jprb, 1.89203e-06_jprb, 1.94360e-06_jprb, 1.99658e-06_jprb, 2.05101e-06_jprb, &
     & 2.10692e-06_jprb, 2.16435e-06_jprb, 2.22335e-06_jprb, 2.28396e-06_jprb, 2.34622e-06_jprb, &
     & 2.41017e-06_jprb, 2.47587e-06_jprb, 2.54337e-06_jprb, 2.61270e-06_jprb, 2.68392e-06_jprb, &
     & 2.75708e-06_jprb, 2.83224e-06_jprb, 2.90944e-06_jprb, 2.98875e-06_jprb/)
      kao_mn2(:, 5) = (/ &
     & 3.41996e-06_jprb, 3.32758e-06_jprb, 3.23770e-06_jprb, 3.15024e-06_jprb, 3.06515e-06_jprb, &
     & 2.98235e-06_jprb, 2.90180e-06_jprb, 2.82341e-06_jprb, 2.74715e-06_jprb, 2.67294e-06_jprb, &
     & 2.60074e-06_jprb, 2.53049e-06_jprb, 2.46214e-06_jprb, 2.39563e-06_jprb, 2.33092e-06_jprb, &
     & 2.26796e-06_jprb, 2.20670e-06_jprb, 2.14709e-06_jprb, 2.08910e-06_jprb/)
      kao_mn2(:, 6) = (/ &
     & 3.38746e-06_jprb, 3.25966e-06_jprb, 3.13669e-06_jprb, 3.01836e-06_jprb, 2.90449e-06_jprb, &
     & 2.79491e-06_jprb, 2.68947e-06_jprb, 2.58801e-06_jprb, 2.49037e-06_jprb, 2.39642e-06_jprb, &
     & 2.30601e-06_jprb, 2.21902e-06_jprb, 2.13530e-06_jprb, 2.05475e-06_jprb, 1.97723e-06_jprb, &
     & 1.90264e-06_jprb, 1.83086e-06_jprb, 1.76179e-06_jprb, 1.69532e-06_jprb/)
      kao_mn2(:, 7) = (/ &
     & 3.17530e-06_jprb, 3.07196e-06_jprb, 2.97199e-06_jprb, 2.87527e-06_jprb, 2.78170e-06_jprb, &
     & 2.69118e-06_jprb, 2.60360e-06_jprb, 2.51887e-06_jprb, 2.43690e-06_jprb, 2.35759e-06_jprb, &
     & 2.28087e-06_jprb, 2.20664e-06_jprb, 2.13483e-06_jprb, 2.06536e-06_jprb, 1.99814e-06_jprb, &
     & 1.93312e-06_jprb, 1.87021e-06_jprb, 1.80934e-06_jprb, 1.75046e-06_jprb/)
      kao_mn2(:, 8) = (/ &
     & 2.84701e-06_jprb, 2.77007e-06_jprb, 2.69521e-06_jprb, 2.62237e-06_jprb, 2.55150e-06_jprb, &
     & 2.48254e-06_jprb, 2.41545e-06_jprb, 2.35017e-06_jprb, 2.28666e-06_jprb, 2.22486e-06_jprb, &
     & 2.16473e-06_jprb, 2.10623e-06_jprb, 2.04930e-06_jprb, 1.99392e-06_jprb, 1.94003e-06_jprb, &
     & 1.88760e-06_jprb, 1.83659e-06_jprb, 1.78695e-06_jprb, 1.73866e-06_jprb/)
      kao_mn2(:, 9) = (/ &
     & 2.79917e-06_jprb, 2.73207e-06_jprb, 2.66658e-06_jprb, 2.60266e-06_jprb, 2.54027e-06_jprb, &
     & 2.47937e-06_jprb, 2.41994e-06_jprb, 2.36192e-06_jprb, 2.30530e-06_jprb, 2.25004e-06_jprb, &
     & 2.19610e-06_jprb, 2.14346e-06_jprb, 2.09208e-06_jprb, 2.04193e-06_jprb, 1.99298e-06_jprb, &
     & 1.94520e-06_jprb, 1.89857e-06_jprb, 1.85306e-06_jprb, 1.80864e-06_jprb/)
      kao_mn2(:,10) = (/ &
     & 2.74910e-06_jprb, 2.64462e-06_jprb, 2.54412e-06_jprb, 2.44743e-06_jprb, 2.35442e-06_jprb, &
     & 2.26495e-06_jprb, 2.17887e-06_jprb, 2.09606e-06_jprb, 2.01641e-06_jprb, 1.93978e-06_jprb, &
     & 1.86606e-06_jprb, 1.79514e-06_jprb, 1.72692e-06_jprb, 1.66129e-06_jprb, 1.59815e-06_jprb, &
     & 1.53742e-06_jprb, 1.47899e-06_jprb, 1.42278e-06_jprb, 1.36871e-06_jprb/)
      kao_mn2(:,11) = (/ &
     & 2.63952e-06_jprb, 2.60263e-06_jprb, 2.56626e-06_jprb, 2.53039e-06_jprb, 2.49503e-06_jprb, &
     & 2.46016e-06_jprb, 2.42578e-06_jprb, 2.39188e-06_jprb, 2.35845e-06_jprb, 2.32549e-06_jprb, &
     & 2.29299e-06_jprb, 2.26094e-06_jprb, 2.22934e-06_jprb, 2.19819e-06_jprb, 2.16747e-06_jprb, &
     & 2.13717e-06_jprb, 2.10731e-06_jprb, 2.07786e-06_jprb, 2.04882e-06_jprb/)
      kao_mn2(:,12) = (/ &
     & 2.94106e-06_jprb, 2.82819e-06_jprb, 2.71966e-06_jprb, 2.61528e-06_jprb, 2.51492e-06_jprb, &
     & 2.41841e-06_jprb, 2.32560e-06_jprb, 2.23635e-06_jprb, 2.15053e-06_jprb, 2.06800e-06_jprb, &
     & 1.98863e-06_jprb, 1.91232e-06_jprb, 1.83893e-06_jprb, 1.76836e-06_jprb, 1.70049e-06_jprb, &
     & 1.63524e-06_jprb, 1.57248e-06_jprb, 1.51214e-06_jprb, 1.45411e-06_jprb/)
      kao_mn2(:,13) = (/ &
     & 2.94607e-06_jprb, 2.87369e-06_jprb, 2.80309e-06_jprb, 2.73422e-06_jprb, 2.66705e-06_jprb, &
     & 2.60152e-06_jprb, 2.53760e-06_jprb, 2.47526e-06_jprb, 2.41445e-06_jprb, 2.35513e-06_jprb, &
     & 2.29726e-06_jprb, 2.24082e-06_jprb, 2.18577e-06_jprb, 2.13207e-06_jprb, 2.07969e-06_jprb, &
     & 2.02859e-06_jprb, 1.97875e-06_jprb, 1.93014e-06_jprb, 1.88272e-06_jprb/)
      kao_mn2(:,14) = (/ &
     & 2.58051e-06_jprb, 2.48749e-06_jprb, 2.39782e-06_jprb, 2.31139e-06_jprb, 2.22807e-06_jprb, &
     & 2.14775e-06_jprb, 2.07033e-06_jprb, 1.99570e-06_jprb, 1.92376e-06_jprb, 1.85441e-06_jprb, &
     & 1.78756e-06_jprb, 1.72313e-06_jprb, 1.66101e-06_jprb, 1.60114e-06_jprb, 1.54342e-06_jprb, &
     & 1.48778e-06_jprb, 1.43415e-06_jprb, 1.38245e-06_jprb, 1.33262e-06_jprb/)
      kao_mn2(:,15) = (/ &
     & 3.03447e-06_jprb, 2.88559e-06_jprb, 2.74401e-06_jprb, 2.60938e-06_jprb, 2.48135e-06_jprb, &
     & 2.35961e-06_jprb, 2.24384e-06_jprb, 2.13375e-06_jprb, 2.02906e-06_jprb, 1.92951e-06_jprb, &
     & 1.83484e-06_jprb, 1.74481e-06_jprb, 1.65921e-06_jprb, 1.57780e-06_jprb, 1.50039e-06_jprb, &
     & 1.42677e-06_jprb, 1.35677e-06_jprb, 1.29020e-06_jprb, 1.22690e-06_jprb/)
      kao_mn2(:,16) = (/ &
     & 1.48655e-06_jprb, 1.48283e-06_jprb, 1.47913e-06_jprb, 1.47543e-06_jprb, 1.47174e-06_jprb, &
     & 1.46806e-06_jprb, 1.46439e-06_jprb, 1.46072e-06_jprb, 1.45707e-06_jprb, 1.45343e-06_jprb, &
     & 1.44979e-06_jprb, 1.44617e-06_jprb, 1.44255e-06_jprb, 1.43894e-06_jprb, 1.43534e-06_jprb, &
     & 1.43176e-06_jprb, 1.42817e-06_jprb, 1.42460e-06_jprb, 1.42104e-06_jprb/)
      kbo_mn2(:, 1) = (/ &
     & 5.12042e-08_jprb, 5.51239e-08_jprb, 5.93436e-08_jprb, 6.38863e-08_jprb, 6.87767e-08_jprb, &
     & 7.40415e-08_jprb, 7.97093e-08_jprb, 8.58110e-08_jprb, 9.23797e-08_jprb, 9.94513e-08_jprb, &
     & 1.07064e-07_jprb, 1.15260e-07_jprb, 1.24083e-07_jprb, 1.33581e-07_jprb, 1.43807e-07_jprb, &
     & 1.54815e-07_jprb, 1.66666e-07_jprb, 1.79424e-07_jprb, 1.93159e-07_jprb/)
      kbo_mn2(:, 2) = (/ &
     & 2.30938e-07_jprb, 2.41696e-07_jprb, 2.52955e-07_jprb, 2.64738e-07_jprb, 2.77071e-07_jprb, &
     & 2.89978e-07_jprb, 3.03486e-07_jprb, 3.17623e-07_jprb, 3.32419e-07_jprb, 3.47904e-07_jprb, &
     & 3.64111e-07_jprb, 3.81072e-07_jprb, 3.98824e-07_jprb, 4.17402e-07_jprb, 4.36846e-07_jprb, &
     & 4.57196e-07_jprb, 4.78494e-07_jprb, 5.00784e-07_jprb, 5.24112e-07_jprb/)
      kbo_mn2(:, 3) = (/ &
     & 6.70458e-07_jprb, 7.04274e-07_jprb, 7.39795e-07_jprb, 7.77109e-07_jprb, 8.16304e-07_jprb, &
     & 8.57476e-07_jprb, 9.00724e-07_jprb, 9.46154e-07_jprb, 9.93876e-07_jprb, 1.04400e-06_jprb, &
     & 1.09666e-06_jprb, 1.15197e-06_jprb, 1.21008e-06_jprb, 1.27111e-06_jprb, 1.33522e-06_jprb, &
     & 1.40256e-06_jprb, 1.47331e-06_jprb, 1.54761e-06_jprb, 1.62567e-06_jprb/)
      kbo_mn2(:, 4) = (/ &
     & 1.84182e-06_jprb, 1.89203e-06_jprb, 1.94360e-06_jprb, 1.99658e-06_jprb, 2.05101e-06_jprb, &
     & 2.10692e-06_jprb, 2.16435e-06_jprb, 2.22335e-06_jprb, 2.28396e-06_jprb, 2.34622e-06_jprb, &
     & 2.41017e-06_jprb, 2.47587e-06_jprb, 2.54337e-06_jprb, 2.61270e-06_jprb, 2.68392e-06_jprb, &
     & 2.75708e-06_jprb, 2.83224e-06_jprb, 2.90944e-06_jprb, 2.98875e-06_jprb/)
      kbo_mn2(:, 5) = (/ &
     & 3.41996e-06_jprb, 3.32758e-06_jprb, 3.23770e-06_jprb, 3.15024e-06_jprb, 3.06515e-06_jprb, &
     & 2.98235e-06_jprb, 2.90180e-06_jprb, 2.82341e-06_jprb, 2.74715e-06_jprb, 2.67294e-06_jprb, &
     & 2.60074e-06_jprb, 2.53049e-06_jprb, 2.46214e-06_jprb, 2.39563e-06_jprb, 2.33092e-06_jprb, &
     & 2.26796e-06_jprb, 2.20670e-06_jprb, 2.14709e-06_jprb, 2.08910e-06_jprb/)
      kbo_mn2(:, 6) = (/ &
     & 3.38746e-06_jprb, 3.25966e-06_jprb, 3.13669e-06_jprb, 3.01836e-06_jprb, 2.90449e-06_jprb, &
     & 2.79491e-06_jprb, 2.68947e-06_jprb, 2.58801e-06_jprb, 2.49037e-06_jprb, 2.39642e-06_jprb, &
     & 2.30601e-06_jprb, 2.21902e-06_jprb, 2.13530e-06_jprb, 2.05475e-06_jprb, 1.97723e-06_jprb, &
     & 1.90264e-06_jprb, 1.83086e-06_jprb, 1.76179e-06_jprb, 1.69532e-06_jprb/)
      kbo_mn2(:, 7) = (/ &
     & 3.17530e-06_jprb, 3.07196e-06_jprb, 2.97199e-06_jprb, 2.87527e-06_jprb, 2.78170e-06_jprb, &
     & 2.69118e-06_jprb, 2.60360e-06_jprb, 2.51887e-06_jprb, 2.43690e-06_jprb, 2.35759e-06_jprb, &
     & 2.28087e-06_jprb, 2.20664e-06_jprb, 2.13483e-06_jprb, 2.06536e-06_jprb, 1.99814e-06_jprb, &
     & 1.93312e-06_jprb, 1.87021e-06_jprb, 1.80934e-06_jprb, 1.75046e-06_jprb/)
      kbo_mn2(:, 8) = (/ &
     & 2.84701e-06_jprb, 2.77007e-06_jprb, 2.69521e-06_jprb, 2.62237e-06_jprb, 2.55150e-06_jprb, &
     & 2.48254e-06_jprb, 2.41545e-06_jprb, 2.35017e-06_jprb, 2.28666e-06_jprb, 2.22486e-06_jprb, &
     & 2.16473e-06_jprb, 2.10623e-06_jprb, 2.04930e-06_jprb, 1.99392e-06_jprb, 1.94003e-06_jprb, &
     & 1.88760e-06_jprb, 1.83659e-06_jprb, 1.78695e-06_jprb, 1.73866e-06_jprb/)
      kbo_mn2(:, 9) = (/ &
     & 2.79917e-06_jprb, 2.73207e-06_jprb, 2.66658e-06_jprb, 2.60266e-06_jprb, 2.54027e-06_jprb, &
     & 2.47937e-06_jprb, 2.41994e-06_jprb, 2.36192e-06_jprb, 2.30530e-06_jprb, 2.25004e-06_jprb, &
     & 2.19610e-06_jprb, 2.14346e-06_jprb, 2.09208e-06_jprb, 2.04193e-06_jprb, 1.99298e-06_jprb, &
     & 1.94520e-06_jprb, 1.89857e-06_jprb, 1.85306e-06_jprb, 1.80864e-06_jprb/)
      kbo_mn2(:,10) = (/ &
     & 2.74910e-06_jprb, 2.64462e-06_jprb, 2.54412e-06_jprb, 2.44743e-06_jprb, 2.35442e-06_jprb, &
     & 2.26495e-06_jprb, 2.17887e-06_jprb, 2.09606e-06_jprb, 2.01641e-06_jprb, 1.93978e-06_jprb, &
     & 1.86606e-06_jprb, 1.79514e-06_jprb, 1.72692e-06_jprb, 1.66129e-06_jprb, 1.59815e-06_jprb, &
     & 1.53742e-06_jprb, 1.47899e-06_jprb, 1.42278e-06_jprb, 1.36871e-06_jprb/)
      kbo_mn2(:,11) = (/ &
     & 2.63952e-06_jprb, 2.60263e-06_jprb, 2.56626e-06_jprb, 2.53039e-06_jprb, 2.49503e-06_jprb, &
     & 2.46016e-06_jprb, 2.42578e-06_jprb, 2.39188e-06_jprb, 2.35845e-06_jprb, 2.32549e-06_jprb, &
     & 2.29299e-06_jprb, 2.26094e-06_jprb, 2.22934e-06_jprb, 2.19819e-06_jprb, 2.16747e-06_jprb, &
     & 2.13717e-06_jprb, 2.10731e-06_jprb, 2.07786e-06_jprb, 2.04882e-06_jprb/)
      kbo_mn2(:,12) = (/ &
     & 2.94106e-06_jprb, 2.82819e-06_jprb, 2.71966e-06_jprb, 2.61528e-06_jprb, 2.51492e-06_jprb, &
     & 2.41841e-06_jprb, 2.32560e-06_jprb, 2.23635e-06_jprb, 2.15053e-06_jprb, 2.06800e-06_jprb, &
     & 1.98863e-06_jprb, 1.91232e-06_jprb, 1.83893e-06_jprb, 1.76836e-06_jprb, 1.70049e-06_jprb, &
     & 1.63524e-06_jprb, 1.57248e-06_jprb, 1.51214e-06_jprb, 1.45411e-06_jprb/)
      kbo_mn2(:,13) = (/ &
     & 2.94607e-06_jprb, 2.87369e-06_jprb, 2.80309e-06_jprb, 2.73422e-06_jprb, 2.66705e-06_jprb, &
     & 2.60152e-06_jprb, 2.53760e-06_jprb, 2.47526e-06_jprb, 2.41445e-06_jprb, 2.35513e-06_jprb, &
     & 2.29726e-06_jprb, 2.24082e-06_jprb, 2.18577e-06_jprb, 2.13207e-06_jprb, 2.07969e-06_jprb, &
     & 2.02859e-06_jprb, 1.97875e-06_jprb, 1.93014e-06_jprb, 1.88272e-06_jprb/)
      kbo_mn2(:,14) = (/ &
     & 2.58051e-06_jprb, 2.48749e-06_jprb, 2.39782e-06_jprb, 2.31139e-06_jprb, 2.22807e-06_jprb, &
     & 2.14775e-06_jprb, 2.07033e-06_jprb, 1.99570e-06_jprb, 1.92376e-06_jprb, 1.85441e-06_jprb, &
     & 1.78756e-06_jprb, 1.72313e-06_jprb, 1.66101e-06_jprb, 1.60114e-06_jprb, 1.54342e-06_jprb, &
     & 1.48778e-06_jprb, 1.43415e-06_jprb, 1.38245e-06_jprb, 1.33262e-06_jprb/)
      kbo_mn2(:,15) = (/ &
     & 3.03447e-06_jprb, 2.88559e-06_jprb, 2.74401e-06_jprb, 2.60938e-06_jprb, 2.48135e-06_jprb, &
     & 2.35961e-06_jprb, 2.24384e-06_jprb, 2.13375e-06_jprb, 2.02906e-06_jprb, 1.92951e-06_jprb, &
     & 1.83484e-06_jprb, 1.74481e-06_jprb, 1.65921e-06_jprb, 1.57780e-06_jprb, 1.50039e-06_jprb, &
     & 1.42677e-06_jprb, 1.35677e-06_jprb, 1.29020e-06_jprb, 1.22690e-06_jprb/)
      kbo_mn2(:,16) = (/ &
     & 1.48655e-06_jprb, 1.48283e-06_jprb, 1.47913e-06_jprb, 1.47543e-06_jprb, 1.47174e-06_jprb, &
     & 1.46806e-06_jprb, 1.46439e-06_jprb, 1.46072e-06_jprb, 1.45707e-06_jprb, 1.45343e-06_jprb, &
     & 1.44979e-06_jprb, 1.44617e-06_jprb, 1.44255e-06_jprb, 1.43894e-06_jprb, 1.43534e-06_jprb, &
     & 1.43176e-06_jprb, 1.42817e-06_jprb, 1.42460e-06_jprb, 1.42104e-06_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 2.16803e+00_jprb, 1.98236e+00_jprb, 1.81260e+00_jprb, 1.65737e+00_jprb, 1.51544e+00_jprb, &
     & 1.38567e+00_jprb, 1.26700e+00_jprb, 1.15850e+00_jprb, 1.05929e+00_jprb, 9.68576e-01_jprb/)
      selfrefo(:, 2) = (/ &
     & 3.70149e+00_jprb, 3.43145e+00_jprb, 3.18110e+00_jprb, 2.94902e+00_jprb, 2.73387e+00_jprb, &
     & 2.53441e+00_jprb, 2.34951e+00_jprb, 2.17810e+00_jprb, 2.01919e+00_jprb, 1.87188e+00_jprb/)
      selfrefo(:, 3) = (/ &
     & 6.17433e+00_jprb, 5.62207e+00_jprb, 5.11920e+00_jprb, 4.66131e+00_jprb, 4.24438e+00_jprb, &
     & 3.86474e+00_jprb, 3.51906e+00_jprb, 3.20430e+00_jprb, 2.91769e+00_jprb, 2.65672e+00_jprb/)
      selfrefo(:, 4) = (/ &
     & 6.56459e+00_jprb, 5.94787e+00_jprb, 5.38910e+00_jprb, 4.88282e+00_jprb, 4.42410e+00_jprb, &
     & 4.00848e+00_jprb, 3.63190e+00_jprb, 3.29070e+00_jprb, 2.98155e+00_jprb, 2.70145e+00_jprb/)
      selfrefo(:, 5) = (/ &
     & 6.49581e+00_jprb, 5.91114e+00_jprb, 5.37910e+00_jprb, 4.89494e+00_jprb, 4.45436e+00_jprb, &
     & 4.05344e+00_jprb, 3.68860e+00_jprb, 3.35660e+00_jprb, 3.05448e+00_jprb, 2.77956e+00_jprb/)
      selfrefo(:, 6) = (/ &
     & 6.50189e+00_jprb, 5.89381e+00_jprb, 5.34260e+00_jprb, 4.84294e+00_jprb, 4.39001e+00_jprb, &
     & 3.97944e+00_jprb, 3.60727e+00_jprb, 3.26990e+00_jprb, 2.96409e+00_jprb, 2.68687e+00_jprb/)
      selfrefo(:, 7) = (/ &
     & 6.64768e+00_jprb, 6.01719e+00_jprb, 5.44650e+00_jprb, 4.92993e+00_jprb, 4.46236e+00_jprb, &
     & 4.03914e+00_jprb, 3.65605e+00_jprb, 3.30930e+00_jprb, 2.99543e+00_jprb, 2.71134e+00_jprb/)
      selfrefo(:, 8) = (/ &
     & 6.43744e+00_jprb, 5.87166e+00_jprb, 5.35560e+00_jprb, 4.88490e+00_jprb, 4.45557e+00_jprb, &
     & 4.06397e+00_jprb, 3.70679e+00_jprb, 3.38100e+00_jprb, 3.08384e+00_jprb, 2.81281e+00_jprb/)
      selfrefo(:, 9) = (/ &
     & 6.55466e+00_jprb, 5.99777e+00_jprb, 5.48820e+00_jprb, 5.02192e+00_jprb, 4.59525e+00_jprb, &
     & 4.20484e+00_jprb, 3.84759e+00_jprb, 3.52070e+00_jprb, 3.22158e+00_jprb, 2.94787e+00_jprb/)
      selfrefo(:,10) = (/ &
     & 6.84510e+00_jprb, 6.26933e+00_jprb, 5.74200e+00_jprb, 5.25902e+00_jprb, 4.81667e+00_jprb, &
     & 4.41152e+00_jprb, 4.04046e+00_jprb, 3.70060e+00_jprb, 3.38933e+00_jprb, 3.10424e+00_jprb/)
      selfrefo(:,11) = (/ &
     & 6.83128e+00_jprb, 6.25536e+00_jprb, 5.72800e+00_jprb, 5.24510e+00_jprb, 4.80291e+00_jprb, &
     & 4.39799e+00_jprb, 4.02722e+00_jprb, 3.68770e+00_jprb, 3.37681e+00_jprb, 3.09212e+00_jprb/)
      selfrefo(:,12) = (/ &
     & 7.35969e+00_jprb, 6.61719e+00_jprb, 5.94960e+00_jprb, 5.34936e+00_jprb, 4.80968e+00_jprb, &
     & 4.32445e+00_jprb, 3.88817e+00_jprb, 3.49590e+00_jprb, 3.14321e+00_jprb, 2.82610e+00_jprb/)
      selfrefo(:,13) = (/ &
     & 7.50064e+00_jprb, 6.80749e+00_jprb, 6.17840e+00_jprb, 5.60744e+00_jprb, 5.08925e+00_jprb, &
     & 4.61894e+00_jprb, 4.19210e+00_jprb, 3.80470e+00_jprb, 3.45310e+00_jprb, 3.13399e+00_jprb/)
      selfrefo(:,14) = (/ &
     & 7.40801e+00_jprb, 6.71328e+00_jprb, 6.08370e+00_jprb, 5.51316e+00_jprb, 4.99613e+00_jprb, &
     & 4.52759e+00_jprb, 4.10298e+00_jprb, 3.71820e+00_jprb, 3.36950e+00_jprb, 3.05351e+00_jprb/)
      selfrefo(:,15) = (/ &
     & 7.51895e+00_jprb, 6.68846e+00_jprb, 5.94970e+00_jprb, 5.29254e+00_jprb, 4.70796e+00_jprb, &
     & 4.18795e+00_jprb, 3.72538e+00_jprb, 3.31390e+00_jprb, 2.94787e+00_jprb, 2.62227e+00_jprb/)
      selfrefo(:,16) = (/ &
     & 7.84774e+00_jprb, 6.80673e+00_jprb, 5.90380e+00_jprb, 5.12065e+00_jprb, 4.44138e+00_jprb, &
     & 3.85223e+00_jprb, 3.34122e+00_jprb, 2.89800e+00_jprb, 2.51357e+00_jprb, 2.18014e+00_jprb/)





if (lhook) call dr_hook('rrtm_kgb1',1,zhook_handle)
return

1000 continue
call abor1("rrtm_kgb1:error opening file radrrtm")
1001 continue
call abor1("rrtm_kgb1:error reading file radrrtm")

!     -----------------------------------------------------------------
end subroutine rrtm_kgb1
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

