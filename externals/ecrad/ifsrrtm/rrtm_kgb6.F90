! # 1 "ifsrrtm/rrtm_kgb6.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb6.f90"
subroutine rrtm_kgb6

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 6:  820-980 cm-1 (low - h2o; high - nothing)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
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

use yoerrto6 , only : kao     ,kao_mco2, selfrefo, forrefo ,fracrefao ,cfc11adjo ,&
 & cfc12o, kao_d 
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
! # 29 "ifsrrtm/rrtm_kgb6.f90" 2

if (lhook) call dr_hook('rrtm_kgb6',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d
  kao = real(kao_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb6:')
endif


! planck fraction mapping level : p = 473.4280 mb, t = 259.83 k
fracrefao(:) = (/ &
  &   1.4353e-01_jprb,1.4774e-01_jprb,1.4467e-01_jprb,1.3785e-01_jprb,1.2376e-01_jprb,1.0214e-01_jprb, &
  &   8.1984e-02_jprb,6.1152e-02_jprb,4.0987e-02_jprb,4.5067e-03_jprb,4.0020e-03_jprb,3.1772e-03_jprb, &
  &   2.3458e-03_jprb,1.5025e-03_jprb,5.7415e-04_jprb,8.2970e-05_jprb/)

! minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12


!     cfc11 is multiplied by 1.385 to account for the 1060-1107 cm-1 band.
cfc11adjo( :) = (/&
 & 0.0_jprb,      0.0_jprb, 36.7627_jprb, 150.757_jprb,    &
 & 81.4109_jprb, 74.9112_jprb, 56.9325_jprb, 49.3226_jprb,  &
 & 57.1074_jprb, 66.1202_jprb, 109.557_jprb, 89.0562_jprb,  &
 & 149.865_jprb, 196.140_jprb, 258.393_jprb, 80.9923_jprb/)  

cfc12o( :) = (/&
 & 62.8368_jprb, 43.2626_jprb, 26.7549_jprb, 22.2487_jprb,&
 & 23.5029_jprb, 34.8323_jprb, 26.2335_jprb, 23.2306_jprb,&
 & 18.4062_jprb, 13.9534_jprb, 22.6268_jprb, 24.2604_jprb,&
 & 30.0088_jprb, 26.3634_jprb, 15.8237_jprb, 57.5050_jprb/)  


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


!     the array kao_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kao_mco2(:, 1) = (/ &
     & 1.45661e-05_jprb, 1.73337e-05_jprb, 2.06273e-05_jprb, 2.45466e-05_jprb, 2.92105e-05_jprb, &
     & 3.47607e-05_jprb, 4.13654e-05_jprb, 4.92251e-05_jprb, 5.85781e-05_jprb, 6.97083e-05_jprb, &
     & 8.29533e-05_jprb, 9.87149e-05_jprb, 1.17471e-04_jprb, 1.39792e-04_jprb, 1.66353e-04_jprb, &
     & 1.97961e-04_jprb, 2.35574e-04_jprb, 2.80335e-04_jprb, 3.33600e-04_jprb/)
      kao_mco2(:, 2) = (/ &
     & 9.96332e-06_jprb, 1.21229e-05_jprb, 1.47506e-05_jprb, 1.79478e-05_jprb, 2.18381e-05_jprb, &
     & 2.65716e-05_jprb, 3.23310e-05_jprb, 3.93389e-05_jprb, 4.78658e-05_jprb, 5.82408e-05_jprb, &
     & 7.08647e-05_jprb, 8.62250e-05_jprb, 1.04914e-04_jprb, 1.27655e-04_jprb, 1.55325e-04_jprb, &
     & 1.88992e-04_jprb, 2.29957e-04_jprb, 2.79801e-04_jprb, 3.40448e-04_jprb/)
      kao_mco2(:, 3) = (/ &
     & 1.14968e-05_jprb, 1.39890e-05_jprb, 1.70215e-05_jprb, 2.07115e-05_jprb, 2.52013e-05_jprb, &
     & 3.06644e-05_jprb, 3.73118e-05_jprb, 4.54002e-05_jprb, 5.52420e-05_jprb, 6.72173e-05_jprb, &
     & 8.17887e-05_jprb, 9.95188e-05_jprb, 1.21092e-04_jprb, 1.47343e-04_jprb, 1.79283e-04_jprb, &
     & 2.18148e-04_jprb, 2.65438e-04_jprb, 3.22980e-04_jprb, 3.92995e-04_jprb/)
      kao_mco2(:, 4) = (/ &
     & 1.02186e-05_jprb, 1.23232e-05_jprb, 1.48613e-05_jprb, 1.79222e-05_jprb, 2.16134e-05_jprb, &
     & 2.60649e-05_jprb, 3.14332e-05_jprb, 3.79071e-05_jprb, 4.57145e-05_jprb, 5.51297e-05_jprb, &
     & 6.64843e-05_jprb, 8.01773e-05_jprb, 9.66905e-05_jprb, 1.16605e-04_jprb, 1.40621e-04_jprb, &
     & 1.69583e-04_jprb, 2.04510e-04_jprb, 2.46631e-04_jprb, 2.97426e-04_jprb/)
      kao_mco2(:, 5) = (/ &
     & 1.03469e-05_jprb, 1.24680e-05_jprb, 1.50239e-05_jprb, 1.81037e-05_jprb, 2.18149e-05_jprb, &
     & 2.62869e-05_jprb, 3.16756e-05_jprb, 3.81690e-05_jprb, 4.59935e-05_jprb, 5.54220e-05_jprb, &
     & 6.67833e-05_jprb, 8.04737e-05_jprb, 9.69704e-05_jprb, 1.16849e-04_jprb, 1.40803e-04_jprb, &
     & 1.69667e-04_jprb, 2.04448e-04_jprb, 2.46359e-04_jprb, 2.96861e-04_jprb/)
      kao_mco2(:, 6) = (/ &
     & 1.71660e-05_jprb, 2.07334e-05_jprb, 2.50420e-05_jprb, 3.02461e-05_jprb, 3.65317e-05_jprb, &
     & 4.41235e-05_jprb, 5.32930e-05_jprb, 6.43680e-05_jprb, 7.77446e-05_jprb, 9.39010e-05_jprb, &
     & 1.13415e-04_jprb, 1.36984e-04_jprb, 1.65451e-04_jprb, 1.99835e-04_jprb, 2.41363e-04_jprb, &
     & 2.91522e-04_jprb, 3.52104e-04_jprb, 4.25276e-04_jprb, 5.13654e-04_jprb/)
      kao_mco2(:, 7) = (/ &
     & 4.78803e-05_jprb, 5.79395e-05_jprb, 7.01119e-05_jprb, 8.48418e-05_jprb, 1.02666e-04_jprb, &
     & 1.24235e-04_jprb, 1.50336e-04_jprb, 1.81920e-04_jprb, 2.20139e-04_jprb, 2.66388e-04_jprb, &
     & 3.22354e-04_jprb, 3.90077e-04_jprb, 4.72028e-04_jprb, 5.71197e-04_jprb, 6.91199e-04_jprb, &
     & 8.36413e-04_jprb, 1.01214e-03_jprb, 1.22477e-03_jprb, 1.48209e-03_jprb/)
      kao_mco2(:, 8) = (/ &
     & 1.27954e-04_jprb, 1.55281e-04_jprb, 1.88445e-04_jprb, 2.28692e-04_jprb, 2.77534e-04_jprb, &
     & 3.36808e-04_jprb, 4.08741e-04_jprb, 4.96037e-04_jprb, 6.01977e-04_jprb, 7.30542e-04_jprb, &
     & 8.86566e-04_jprb, 1.07591e-03_jprb, 1.30570e-03_jprb, 1.58456e-03_jprb, 1.92298e-03_jprb, &
     & 2.33367e-03_jprb, 2.83208e-03_jprb, 3.43694e-03_jprb, 4.17097e-03_jprb/)
      kao_mco2(:, 9) = (/ &
     & 2.93792e-05_jprb, 3.55109e-05_jprb, 4.29223e-05_jprb, 5.18805e-05_jprb, 6.27083e-05_jprb, &
     & 7.57960e-05_jprb, 9.16151e-05_jprb, 1.10736e-04_jprb, 1.33847e-04_jprb, 1.61782e-04_jprb, &
     & 1.95547e-04_jprb, 2.36359e-04_jprb, 2.85689e-04_jprb, 3.45315e-04_jprb, 4.17384e-04_jprb, &
     & 5.04495e-04_jprb, 6.09787e-04_jprb, 7.37054e-04_jprb, 8.90882e-04_jprb/)
      kao_mco2(:,10) = (/ &
     & 5.08569e-05_jprb, 6.24700e-05_jprb, 7.67350e-05_jprb, 9.42574e-05_jprb, 1.15781e-04_jprb, &
     & 1.42220e-04_jprb, 1.74695e-04_jprb, 2.14587e-04_jprb, 2.63588e-04_jprb, 3.23778e-04_jprb, &
     & 3.97712e-04_jprb, 4.88530e-04_jprb, 6.00085e-04_jprb, 7.37114e-04_jprb, 9.05433e-04_jprb, &
     & 1.11219e-03_jprb, 1.36616e-03_jprb, 1.67812e-03_jprb, 2.06131e-03_jprb/)
      kao_mco2(:,11) = (/ &
     & 4.82546e-06_jprb, 6.21462e-06_jprb, 8.00369e-06_jprb, 1.03078e-05_jprb, 1.32752e-05_jprb, &
     & 1.70969e-05_jprb, 2.20188e-05_jprb, 2.83575e-05_jprb, 3.65211e-05_jprb, 4.70348e-05_jprb, &
     & 6.05753e-05_jprb, 7.80138e-05_jprb, 1.00472e-04_jprb, 1.29397e-04_jprb, 1.66647e-04_jprb, &
     & 2.14622e-04_jprb, 2.76407e-04_jprb, 3.55980e-04_jprb, 4.58459e-04_jprb/)
      kao_mco2(:,12) = (/ &
     & 2.41346e-06_jprb, 2.96282e-06_jprb, 3.63723e-06_jprb, 4.46516e-06_jprb, 5.48153e-06_jprb, &
     & 6.72926e-06_jprb, 8.26100e-06_jprb, 1.01414e-05_jprb, 1.24498e-05_jprb, 1.52837e-05_jprb, &
     & 1.87627e-05_jprb, 2.30335e-05_jprb, 2.82765e-05_jprb, 3.47129e-05_jprb, 4.26144e-05_jprb, &
     & 5.23144e-05_jprb, 6.42225e-05_jprb, 7.88410e-05_jprb, 9.67871e-05_jprb/)
      kao_mco2(:,13) = (/ &
     & 2.76412e-06_jprb, 3.46195e-06_jprb, 4.33596e-06_jprb, 5.43062e-06_jprb, 6.80164e-06_jprb, &
     & 8.51879e-06_jprb, 1.06695e-05_jprb, 1.33631e-05_jprb, 1.67367e-05_jprb, 2.09621e-05_jprb, &
     & 2.62542e-05_jprb, 3.28824e-05_jprb, 4.11839e-05_jprb, 5.15813e-05_jprb, 6.46035e-05_jprb, &
     & 8.09134e-05_jprb, 1.01341e-04_jprb, 1.26925e-04_jprb, 1.58969e-04_jprb/)
      kao_mco2(:,14) = (/ &
     & 1.25126e-06_jprb, 1.54971e-06_jprb, 1.91935e-06_jprb, 2.37715e-06_jprb, 2.94416e-06_jprb, &
     & 3.64640e-06_jprb, 4.51615e-06_jprb, 5.59335e-06_jprb, 6.92749e-06_jprb, 8.57985e-06_jprb, &
     & 1.06263e-05_jprb, 1.31610e-05_jprb, 1.63001e-05_jprb, 2.01881e-05_jprb, 2.50034e-05_jprb, &
     & 3.09672e-05_jprb, 3.83536e-05_jprb, 4.75018e-05_jprb, 5.88319e-05_jprb/)
      kao_mco2(:,15) = (/ &
     & 1.59748e-06_jprb, 2.08378e-06_jprb, 2.71812e-06_jprb, 3.54557e-06_jprb, 4.62491e-06_jprb, &
     & 6.03282e-06_jprb, 7.86932e-06_jprb, 1.02649e-05_jprb, 1.33897e-05_jprb, 1.74658e-05_jprb, &
     & 2.27827e-05_jprb, 2.97182e-05_jprb, 3.87649e-05_jprb, 5.05657e-05_jprb, 6.59589e-05_jprb, &
     & 8.60380e-05_jprb, 1.12230e-04_jprb, 1.46394e-04_jprb, 1.90959e-04_jprb/)
      kao_mco2(:,16) = (/ &
     & 1.68148e-06_jprb, 2.17133e-06_jprb, 2.80388e-06_jprb, 3.62071e-06_jprb, 4.67549e-06_jprb, &
     & 6.03756e-06_jprb, 7.79642e-06_jprb, 1.00677e-05_jprb, 1.30006e-05_jprb, 1.67879e-05_jprb, &
     & 2.16786e-05_jprb, 2.79941e-05_jprb, 3.61493e-05_jprb, 4.66803e-05_jprb, 6.02792e-05_jprb, &
     & 7.78398e-05_jprb, 1.00516e-04_jprb, 1.29799e-04_jprb, 1.67612e-04_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &3.2710e-07_jprb,5.2119e-07_jprb,8.4740e-07_jprb,1.6908e-06_jprb,2.3433e-06_jprb,4.4129e-06_jprb, &
     &3.8930e-06_jprb,2.3338e-06_jprb,2.4115e-06_jprb,2.4271e-06_jprb,2.4836e-06_jprb,2.6470e-06_jprb, &
     &2.9559e-06_jprb,2.3940e-06_jprb,2.9711e-06_jprb,2.9511e-06_jprb/)
      forrefo(2,:) = (/ &
     &6.5125e-07_jprb,1.2128e-06_jprb,1.7249e-06_jprb,2.7126e-06_jprb,3.1780e-06_jprb,2.1444e-06_jprb, &
     &1.8265e-06_jprb,1.7385e-06_jprb,1.4574e-06_jprb,1.6135e-06_jprb,2.4966e-06_jprb,2.8127e-06_jprb, &
     &2.5229e-06_jprb,2.3251e-06_jprb,2.5353e-06_jprb,3.0200e-06_jprb/)
      forrefo(3,:) = (/ &
     &1.4969e-06_jprb,1.8516e-06_jprb,2.5791e-06_jprb,2.7846e-06_jprb,1.9789e-06_jprb,1.6688e-06_jprb, &
     &1.1037e-06_jprb,9.9065e-07_jprb,1.1557e-06_jprb,7.0847e-07_jprb,5.7758e-07_jprb,4.0425e-07_jprb, &
     &3.2427e-07_jprb,3.2267e-07_jprb,3.1444e-07_jprb,2.6046e-07_jprb/)
      forrefo(4,:) = (/ &
     &1.7567e-06_jprb,1.6891e-06_jprb,2.1003e-06_jprb,2.0957e-06_jprb,2.3664e-06_jprb,2.1538e-06_jprb, &
     &1.5275e-06_jprb,1.0487e-06_jprb,8.7390e-07_jprb,7.9360e-07_jprb,7.7778e-07_jprb,8.1445e-07_jprb, &
     &8.2121e-07_jprb,5.4395e-07_jprb,3.1273e-07_jprb,3.1848e-07_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 7.73921e-02_jprb, 6.45225e-02_jprb, 5.37930e-02_jprb, 4.48477e-02_jprb, 3.73900e-02_jprb, &
     & 3.11723e-02_jprb, 2.59887e-02_jprb, 2.16670e-02_jprb, 1.80640e-02_jprb, 1.50601e-02_jprb/)
      selfrefo(:, 2) = (/ &
     & 8.47756e-02_jprb, 7.10616e-02_jprb, 5.95660e-02_jprb, 4.99301e-02_jprb, 4.18529e-02_jprb, &
     & 3.50824e-02_jprb, 2.94072e-02_jprb, 2.46500e-02_jprb, 2.06624e-02_jprb, 1.73199e-02_jprb/)
      selfrefo(:, 3) = (/ &
     & 8.84829e-02_jprb, 7.46093e-02_jprb, 6.29110e-02_jprb, 5.30469e-02_jprb, 4.47295e-02_jprb, &
     & 3.77161e-02_jprb, 3.18025e-02_jprb, 2.68160e-02_jprb, 2.26114e-02_jprb, 1.90661e-02_jprb/)
      selfrefo(:, 4) = (/ &
     & 9.27003e-02_jprb, 7.88864e-02_jprb, 6.71310e-02_jprb, 5.71273e-02_jprb, 4.86144e-02_jprb, &
     & 4.13700e-02_jprb, 3.52052e-02_jprb, 2.99590e-02_jprb, 2.54946e-02_jprb, 2.16955e-02_jprb/)
      selfrefo(:, 5) = (/ &
     & 9.14315e-02_jprb, 7.85661e-02_jprb, 6.75110e-02_jprb, 5.80115e-02_jprb, 4.98487e-02_jprb, &
     & 4.28344e-02_jprb, 3.68072e-02_jprb, 3.16280e-02_jprb, 2.71776e-02_jprb, 2.33534e-02_jprb/)
      selfrefo(:, 6) = (/ &
     & 7.72984e-02_jprb, 6.91044e-02_jprb, 6.17790e-02_jprb, 5.52301e-02_jprb, 4.93755e-02_jprb, &
     & 4.41414e-02_jprb, 3.94622e-02_jprb, 3.52790e-02_jprb, 3.15392e-02_jprb, 2.81959e-02_jprb/)
      selfrefo(:, 7) = (/ &
     & 7.46998e-02_jprb, 6.66597e-02_jprb, 5.94850e-02_jprb, 5.30825e-02_jprb, 4.73691e-02_jprb, &
     & 4.22707e-02_jprb, 3.77210e-02_jprb, 3.36610e-02_jprb, 3.00380e-02_jprb, 2.68049e-02_jprb/)
      selfrefo(:, 8) = (/ &
     & 7.59386e-02_jprb, 6.66263e-02_jprb, 5.84560e-02_jprb, 5.12876e-02_jprb, 4.49982e-02_jprb, &
     & 3.94801e-02_jprb, 3.46387e-02_jprb, 3.03910e-02_jprb, 2.66642e-02_jprb, 2.33944e-02_jprb/)
      selfrefo(:, 9) = (/ &
     & 7.26921e-02_jprb, 6.43261e-02_jprb, 5.69230e-02_jprb, 5.03719e-02_jprb, 4.45747e-02_jprb, &
     & 3.94447e-02_jprb, 3.49051e-02_jprb, 3.08880e-02_jprb, 2.73332e-02_jprb, 2.41875e-02_jprb/)
      selfrefo(:,10) = (/ &
     & 7.43684e-02_jprb, 6.58735e-02_jprb, 5.83490e-02_jprb, 5.16840e-02_jprb, 4.57803e-02_jprb, &
     & 4.05509e-02_jprb, 3.59189e-02_jprb, 3.18160e-02_jprb, 2.81818e-02_jprb, 2.49626e-02_jprb/)
      selfrefo(:,11) = (/ &
     & 8.97599e-02_jprb, 7.73727e-02_jprb, 6.66950e-02_jprb, 5.74908e-02_jprb, 4.95569e-02_jprb, &
     & 4.27179e-02_jprb, 3.68227e-02_jprb, 3.17410e-02_jprb, 2.73606e-02_jprb, 2.35848e-02_jprb/)
      selfrefo(:,12) = (/ &
     & 9.12262e-02_jprb, 7.84848e-02_jprb, 6.75230e-02_jprb, 5.80922e-02_jprb, 4.99786e-02_jprb, &
     & 4.29982e-02_jprb, 3.69927e-02_jprb, 3.18260e-02_jprb, 2.73809e-02_jprb, 2.35567e-02_jprb/)
      selfrefo(:,13) = (/ &
     & 9.03254e-02_jprb, 7.83291e-02_jprb, 6.79260e-02_jprb, 5.89046e-02_jprb, 5.10813e-02_jprb, &
     & 4.42970e-02_jprb, 3.84139e-02_jprb, 3.33120e-02_jprb, 2.88877e-02_jprb, 2.50511e-02_jprb/)
      selfrefo(:,14) = (/ &
     & 9.22803e-02_jprb, 7.94172e-02_jprb, 6.83470e-02_jprb, 5.88199e-02_jprb, 5.06209e-02_jprb, &
     & 4.35647e-02_jprb, 3.74921e-02_jprb, 3.22660e-02_jprb, 2.77684e-02_jprb, 2.38977e-02_jprb/)
      selfrefo(:,15) = (/ &
     & 9.36819e-02_jprb, 8.10810e-02_jprb, 7.01750e-02_jprb, 6.07359e-02_jprb, 5.25665e-02_jprb, &
     & 4.54959e-02_jprb, 3.93764e-02_jprb, 3.40800e-02_jprb, 2.94960e-02_jprb, 2.55286e-02_jprb/)
      selfrefo(:,16) = (/ &
     & 1.00195e-01_jprb, 8.58713e-02_jprb, 7.35950e-02_jprb, 6.30737e-02_jprb, 5.40566e-02_jprb, &
     & 4.63286e-02_jprb, 3.97054e-02_jprb, 3.40290e-02_jprb, 2.91641e-02_jprb, 2.49948e-02_jprb/)

if (lhook) call dr_hook('rrtm_kgb6',1,zhook_handle)

return

1001 continue
call abor1("rrtm_kgb6:error reading file radrrtm")

end subroutine rrtm_kgb6
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

