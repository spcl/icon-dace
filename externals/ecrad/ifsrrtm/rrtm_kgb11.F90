! # 1 "ifsrrtm/rrtm_kgb11.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb11.f90"
subroutine rrtm_kgb11

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 11:  1480-1800 cm-1 (low - h2o; high - h2o)
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
use yommp0    , only : nproc, myproc
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad

use yoerrto11, only : kao     ,kbo     ,selfrefo, forrefo, fracrefao ,fracrefbo, &
                    & kao_mo2,kbo_mo2, kao_d, kbo_d

!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 29 "ifsrrtm/rrtm_kgb11.f90" 2

if (lhook) call dr_hook('rrtm_kgb11',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb11:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb11:')
endif

! planck fraction mapping level : p=1053.63 mb, t= 294.2 k
      fracrefao(:) = (/ &
     &  1.4601e-01_jprb,1.3824e-01_jprb,1.4240e-01_jprb,1.3463e-01_jprb,1.1948e-01_jprb,1.0440e-01_jprb, &
     &  8.8667e-02_jprb,6.5792e-02_jprb,4.3893e-02_jprb,4.7941e-03_jprb,4.0760e-03_jprb,3.3207e-03_jprb, &
     &  2.4087e-03_jprb,1.3912e-03_jprb,4.3482e-04_jprb,6.0932e-05_jprb/)

! planck fraction mapping level : p=0.353 mb, t = 262.11 k
      fracrefbo(:) = (/ &
     &  7.2928e-02_jprb,1.4900e-01_jprb,1.6156e-01_jprb,1.5603e-01_jprb,1.3934e-01_jprb,1.1394e-01_jprb, &
     &  8.8783e-02_jprb,6.2411e-02_jprb,4.0191e-02_jprb,4.4587e-03_jprb,3.9533e-03_jprb,3.0847e-03_jprb, &
     &  2.2317e-03_jprb,1.4410e-03_jprb,5.6722e-04_jprb,7.7933e-05_jprb/)

!     ------------------------------------------------------------------

!     the array ka contains absorption coefs at the 16 chosen g-values 
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


!     the array kao_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level below 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kao_mo2(:, 1) = (/ &
     & 2.31723e-06_jprb, 2.28697e-06_jprb, 2.25710e-06_jprb, 2.22762e-06_jprb, 2.19852e-06_jprb, &
     & 2.16981e-06_jprb, 2.14147e-06_jprb, 2.11350e-06_jprb, 2.08590e-06_jprb, 2.05865e-06_jprb, &
     & 2.03176e-06_jprb, 2.00523e-06_jprb, 1.97904e-06_jprb, 1.95319e-06_jprb, 1.92768e-06_jprb, &
     & 1.90250e-06_jprb, 1.87765e-06_jprb, 1.85313e-06_jprb, 1.82893e-06_jprb/)
      kao_mo2(:, 2) = (/ &
     & 1.81980e-06_jprb, 1.81352e-06_jprb, 1.80726e-06_jprb, 1.80101e-06_jprb, 1.79479e-06_jprb, &
     & 1.78860e-06_jprb, 1.78242e-06_jprb, 1.77626e-06_jprb, 1.77013e-06_jprb, 1.76402e-06_jprb, &
     & 1.75792e-06_jprb, 1.75185e-06_jprb, 1.74580e-06_jprb, 1.73977e-06_jprb, 1.73377e-06_jprb, &
     & 1.72778e-06_jprb, 1.72181e-06_jprb, 1.71587e-06_jprb, 1.70994e-06_jprb/)
      kao_mo2(:, 3) = (/ &
     & 2.26922e-06_jprb, 2.25413e-06_jprb, 2.23914e-06_jprb, 2.22425e-06_jprb, 2.20945e-06_jprb, &
     & 2.19476e-06_jprb, 2.18016e-06_jprb, 2.16566e-06_jprb, 2.15126e-06_jprb, 2.13695e-06_jprb, &
     & 2.12274e-06_jprb, 2.10862e-06_jprb, 2.09459e-06_jprb, 2.08066e-06_jprb, 2.06683e-06_jprb, &
     & 2.05308e-06_jprb, 2.03942e-06_jprb, 2.02586e-06_jprb, 2.01239e-06_jprb/)
      kao_mo2(:, 4) = (/ &
     & 2.15555e-06_jprb, 2.14539e-06_jprb, 2.13527e-06_jprb, 2.12520e-06_jprb, 2.11517e-06_jprb, &
     & 2.10520e-06_jprb, 2.09527e-06_jprb, 2.08538e-06_jprb, 2.07555e-06_jprb, 2.06576e-06_jprb, &
     & 2.05601e-06_jprb, 2.04631e-06_jprb, 2.03666e-06_jprb, 2.02706e-06_jprb, 2.01749e-06_jprb, &
     & 2.00798e-06_jprb, 1.99851e-06_jprb, 1.98908e-06_jprb, 1.97970e-06_jprb/)
      kao_mo2(:, 5) = (/ &
     & 2.05821e-06_jprb, 2.04914e-06_jprb, 2.04011e-06_jprb, 2.03111e-06_jprb, 2.02216e-06_jprb, &
     & 2.01324e-06_jprb, 2.00437e-06_jprb, 1.99553e-06_jprb, 1.98673e-06_jprb, 1.97798e-06_jprb, &
     & 1.96926e-06_jprb, 1.96057e-06_jprb, 1.95193e-06_jprb, 1.94333e-06_jprb, 1.93476e-06_jprb, &
     & 1.92623e-06_jprb, 1.91774e-06_jprb, 1.90928e-06_jprb, 1.90087e-06_jprb/)
      kao_mo2(:, 6) = (/ &
     & 2.20148e-06_jprb, 2.18998e-06_jprb, 2.17854e-06_jprb, 2.16717e-06_jprb, 2.15585e-06_jprb, &
     & 2.14459e-06_jprb, 2.13339e-06_jprb, 2.12225e-06_jprb, 2.11117e-06_jprb, 2.10014e-06_jprb, &
     & 2.08918e-06_jprb, 2.07827e-06_jprb, 2.06741e-06_jprb, 2.05662e-06_jprb, 2.04588e-06_jprb, &
     & 2.03519e-06_jprb, 2.02457e-06_jprb, 2.01399e-06_jprb, 2.00348e-06_jprb/)
      kao_mo2(:, 7) = (/ &
     & 2.28960e-06_jprb, 2.27651e-06_jprb, 2.26349e-06_jprb, 2.25054e-06_jprb, 2.23767e-06_jprb, &
     & 2.22487e-06_jprb, 2.21215e-06_jprb, 2.19950e-06_jprb, 2.18692e-06_jprb, 2.17441e-06_jprb, &
     & 2.16198e-06_jprb, 2.14961e-06_jprb, 2.13732e-06_jprb, 2.12509e-06_jprb, 2.11294e-06_jprb, &
     & 2.10085e-06_jprb, 2.08884e-06_jprb, 2.07689e-06_jprb, 2.06501e-06_jprb/)
      kao_mo2(:, 8) = (/ &
     & 2.28564e-06_jprb, 2.27363e-06_jprb, 2.26168e-06_jprb, 2.24980e-06_jprb, 2.23798e-06_jprb, &
     & 2.22622e-06_jprb, 2.21452e-06_jprb, 2.20288e-06_jprb, 2.19131e-06_jprb, 2.17980e-06_jprb, &
     & 2.16834e-06_jprb, 2.15695e-06_jprb, 2.14562e-06_jprb, 2.13434e-06_jprb, 2.12313e-06_jprb, &
     & 2.11197e-06_jprb, 2.10087e-06_jprb, 2.08984e-06_jprb, 2.07886e-06_jprb/)
      kao_mo2(:, 9) = (/ &
     & 2.28505e-06_jprb, 2.27395e-06_jprb, 2.26291e-06_jprb, 2.25192e-06_jprb, 2.24099e-06_jprb, &
     & 2.23011e-06_jprb, 2.21928e-06_jprb, 2.20850e-06_jprb, 2.19778e-06_jprb, 2.18711e-06_jprb, &
     & 2.17649e-06_jprb, 2.16592e-06_jprb, 2.15540e-06_jprb, 2.14494e-06_jprb, 2.13452e-06_jprb, &
     & 2.12416e-06_jprb, 2.11385e-06_jprb, 2.10358e-06_jprb, 2.09337e-06_jprb/)
      kao_mo2(:,10) = (/ &
     & 2.25915e-06_jprb, 2.24938e-06_jprb, 2.23965e-06_jprb, 2.22997e-06_jprb, 2.22032e-06_jprb, &
     & 2.21072e-06_jprb, 2.20116e-06_jprb, 2.19164e-06_jprb, 2.18216e-06_jprb, 2.17272e-06_jprb, &
     & 2.16333e-06_jprb, 2.15397e-06_jprb, 2.14465e-06_jprb, 2.13538e-06_jprb, 2.12614e-06_jprb, &
     & 2.11695e-06_jprb, 2.10779e-06_jprb, 2.09868e-06_jprb, 2.08960e-06_jprb/)
      kao_mo2(:,11) = (/ &
     & 2.52025e-06_jprb, 2.50423e-06_jprb, 2.48831e-06_jprb, 2.47249e-06_jprb, 2.45677e-06_jprb, &
     & 2.44115e-06_jprb, 2.42563e-06_jprb, 2.41021e-06_jprb, 2.39489e-06_jprb, 2.37967e-06_jprb, &
     & 2.36454e-06_jprb, 2.34951e-06_jprb, 2.33457e-06_jprb, 2.31973e-06_jprb, 2.30498e-06_jprb, &
     & 2.29033e-06_jprb, 2.27577e-06_jprb, 2.26130e-06_jprb, 2.24692e-06_jprb/)
      kao_mo2(:,12) = (/ &
     & 2.52634e-06_jprb, 2.51180e-06_jprb, 2.49735e-06_jprb, 2.48299e-06_jprb, 2.46871e-06_jprb, &
     & 2.45451e-06_jprb, 2.44039e-06_jprb, 2.42635e-06_jprb, 2.41239e-06_jprb, 2.39851e-06_jprb, &
     & 2.38472e-06_jprb, 2.37100e-06_jprb, 2.35736e-06_jprb, 2.34380e-06_jprb, 2.33032e-06_jprb, &
     & 2.31691e-06_jprb, 2.30358e-06_jprb, 2.29033e-06_jprb, 2.27716e-06_jprb/)
      kao_mo2(:,13) = (/ &
     & 2.66614e-06_jprb, 2.64897e-06_jprb, 2.63191e-06_jprb, 2.61496e-06_jprb, 2.59812e-06_jprb, &
     & 2.58138e-06_jprb, 2.56476e-06_jprb, 2.54824e-06_jprb, 2.53183e-06_jprb, 2.51552e-06_jprb, &
     & 2.49932e-06_jprb, 2.48322e-06_jprb, 2.46723e-06_jprb, 2.45134e-06_jprb, 2.43555e-06_jprb, &
     & 2.41987e-06_jprb, 2.40428e-06_jprb, 2.38880e-06_jprb, 2.37341e-06_jprb/)
      kao_mo2(:,14) = (/ &
     & 2.96755e-06_jprb, 2.94803e-06_jprb, 2.92864e-06_jprb, 2.90937e-06_jprb, 2.89023e-06_jprb, &
     & 2.87122e-06_jprb, 2.85233e-06_jprb, 2.83357e-06_jprb, 2.81493e-06_jprb, 2.79641e-06_jprb, &
     & 2.77802e-06_jprb, 2.75974e-06_jprb, 2.74159e-06_jprb, 2.72355e-06_jprb, 2.70563e-06_jprb, &
     & 2.68784e-06_jprb, 2.67015e-06_jprb, 2.65259e-06_jprb, 2.63514e-06_jprb/)
      kao_mo2(:,15) = (/ &
     & 1.30668e-06_jprb, 1.31378e-06_jprb, 1.32091e-06_jprb, 1.32808e-06_jprb, 1.33530e-06_jprb, &
     & 1.34255e-06_jprb, 1.34984e-06_jprb, 1.35717e-06_jprb, 1.36454e-06_jprb, 1.37195e-06_jprb, &
     & 1.37941e-06_jprb, 1.38690e-06_jprb, 1.39443e-06_jprb, 1.40200e-06_jprb, 1.40962e-06_jprb, &
     & 1.41727e-06_jprb, 1.42497e-06_jprb, 1.43271e-06_jprb, 1.44049e-06_jprb/)
      kao_mo2(:,16) = (/ &
     & 5.99001e-07_jprb, 6.16844e-07_jprb, 6.35219e-07_jprb, 6.54141e-07_jprb, 6.73626e-07_jprb, &
     & 6.93692e-07_jprb, 7.14356e-07_jprb, 7.35635e-07_jprb, 7.57548e-07_jprb, 7.80114e-07_jprb, &
     & 8.03352e-07_jprb, 8.27282e-07_jprb, 8.51925e-07_jprb, 8.77302e-07_jprb, 9.03435e-07_jprb, &
     & 9.30347e-07_jprb, 9.58060e-07_jprb, 9.86599e-07_jprb, 1.01599e-06_jprb/)

!     the array kbo_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kbo_mo2(:, 1) = (/ &
     & 4.97626e-07_jprb, 5.05955e-07_jprb, 5.14424e-07_jprb, 5.23034e-07_jprb, 5.31789e-07_jprb, &
     & 5.40690e-07_jprb, 5.49739e-07_jprb, 5.58941e-07_jprb, 5.68296e-07_jprb, 5.77808e-07_jprb, &
     & 5.87479e-07_jprb, 5.97312e-07_jprb, 6.07310e-07_jprb, 6.17475e-07_jprb, 6.27810e-07_jprb, &
     & 6.38318e-07_jprb, 6.49002e-07_jprb, 6.59865e-07_jprb, 6.70910e-07_jprb/)
      kbo_mo2(:, 2) = (/ &
     & 3.10232e-06_jprb, 3.06339e-06_jprb, 3.02496e-06_jprb, 2.98700e-06_jprb, 2.94952e-06_jprb, &
     & 2.91252e-06_jprb, 2.87597e-06_jprb, 2.83989e-06_jprb, 2.80426e-06_jprb, 2.76907e-06_jprb, &
     & 2.73433e-06_jprb, 2.70002e-06_jprb, 2.66614e-06_jprb, 2.63269e-06_jprb, 2.59966e-06_jprb, &
     & 2.56704e-06_jprb, 2.53483e-06_jprb, 2.50303e-06_jprb, 2.47162e-06_jprb/)
      kbo_mo2(:, 3) = (/ &
     & 2.91635e-06_jprb, 2.88637e-06_jprb, 2.85669e-06_jprb, 2.82733e-06_jprb, 2.79826e-06_jprb, &
     & 2.76949e-06_jprb, 2.74102e-06_jprb, 2.71284e-06_jprb, 2.68495e-06_jprb, 2.65735e-06_jprb, &
     & 2.63003e-06_jprb, 2.60299e-06_jprb, 2.57623e-06_jprb, 2.54975e-06_jprb, 2.52353e-06_jprb, &
     & 2.49759e-06_jprb, 2.47191e-06_jprb, 2.44650e-06_jprb, 2.42135e-06_jprb/)
      kbo_mo2(:, 4) = (/ &
     & 3.15584e-06_jprb, 3.11986e-06_jprb, 3.08430e-06_jprb, 3.04914e-06_jprb, 3.01438e-06_jprb, &
     & 2.98002e-06_jprb, 2.94605e-06_jprb, 2.91247e-06_jprb, 2.87927e-06_jprb, 2.84645e-06_jprb, &
     & 2.81400e-06_jprb, 2.78192e-06_jprb, 2.75021e-06_jprb, 2.71886e-06_jprb, 2.68787e-06_jprb, &
     & 2.65723e-06_jprb, 2.62694e-06_jprb, 2.59699e-06_jprb, 2.56739e-06_jprb/)
      kbo_mo2(:, 5) = (/ &
     & 2.52067e-06_jprb, 2.50127e-06_jprb, 2.48202e-06_jprb, 2.46291e-06_jprb, 2.44396e-06_jprb, &
     & 2.42515e-06_jprb, 2.40648e-06_jprb, 2.38796e-06_jprb, 2.36958e-06_jprb, 2.35134e-06_jprb, &
     & 2.33324e-06_jprb, 2.31529e-06_jprb, 2.29747e-06_jprb, 2.27978e-06_jprb, 2.26224e-06_jprb, &
     & 2.24482e-06_jprb, 2.22755e-06_jprb, 2.21040e-06_jprb, 2.19339e-06_jprb/)
      kbo_mo2(:, 6) = (/ &
     & 2.37304e-06_jprb, 2.36340e-06_jprb, 2.35380e-06_jprb, 2.34423e-06_jprb, 2.33471e-06_jprb, &
     & 2.32522e-06_jprb, 2.31578e-06_jprb, 2.30637e-06_jprb, 2.29700e-06_jprb, 2.28766e-06_jprb, &
     & 2.27837e-06_jprb, 2.26911e-06_jprb, 2.25989e-06_jprb, 2.25071e-06_jprb, 2.24157e-06_jprb, &
     & 2.23246e-06_jprb, 2.22339e-06_jprb, 2.21436e-06_jprb, 2.20536e-06_jprb/)
      kbo_mo2(:, 7) = (/ &
     & 2.56366e-06_jprb, 2.56395e-06_jprb, 2.56424e-06_jprb, 2.56453e-06_jprb, 2.56482e-06_jprb, &
     & 2.56510e-06_jprb, 2.56539e-06_jprb, 2.56568e-06_jprb, 2.56597e-06_jprb, 2.56625e-06_jprb, &
     & 2.56654e-06_jprb, 2.56683e-06_jprb, 2.56712e-06_jprb, 2.56741e-06_jprb, 2.56769e-06_jprb, &
     & 2.56798e-06_jprb, 2.56827e-06_jprb, 2.56856e-06_jprb, 2.56885e-06_jprb/)
      kbo_mo2(:, 8) = (/ &
     & 2.54502e-06_jprb, 2.55393e-06_jprb, 2.56287e-06_jprb, 2.57185e-06_jprb, 2.58085e-06_jprb, &
     & 2.58989e-06_jprb, 2.59896e-06_jprb, 2.60806e-06_jprb, 2.61719e-06_jprb, 2.62636e-06_jprb, &
     & 2.63555e-06_jprb, 2.64478e-06_jprb, 2.65404e-06_jprb, 2.66334e-06_jprb, 2.67266e-06_jprb, &
     & 2.68202e-06_jprb, 2.69141e-06_jprb, 2.70084e-06_jprb, 2.71030e-06_jprb/)
      kbo_mo2(:, 9) = (/ &
     & 1.84106e-06_jprb, 1.83922e-06_jprb, 1.83737e-06_jprb, 1.83553e-06_jprb, 1.83369e-06_jprb, &
     & 1.83186e-06_jprb, 1.83002e-06_jprb, 1.82819e-06_jprb, 1.82636e-06_jprb, 1.82453e-06_jprb, &
     & 1.82270e-06_jprb, 1.82087e-06_jprb, 1.81905e-06_jprb, 1.81723e-06_jprb, 1.81541e-06_jprb, &
     & 1.81359e-06_jprb, 1.81177e-06_jprb, 1.80996e-06_jprb, 1.80814e-06_jprb/)
      kbo_mo2(:,10) = (/ &
     & 1.83886e-06_jprb, 1.83632e-06_jprb, 1.83379e-06_jprb, 1.83126e-06_jprb, 1.82874e-06_jprb, &
     & 1.82622e-06_jprb, 1.82370e-06_jprb, 1.82119e-06_jprb, 1.81868e-06_jprb, 1.81617e-06_jprb, &
     & 1.81367e-06_jprb, 1.81117e-06_jprb, 1.80867e-06_jprb, 1.80618e-06_jprb, 1.80369e-06_jprb, &
     & 1.80120e-06_jprb, 1.79872e-06_jprb, 1.79624e-06_jprb, 1.79377e-06_jprb/)
      kbo_mo2(:,11) = (/ &
     & 2.30390e-06_jprb, 2.30269e-06_jprb, 2.30148e-06_jprb, 2.30028e-06_jprb, 2.29907e-06_jprb, &
     & 2.29787e-06_jprb, 2.29667e-06_jprb, 2.29546e-06_jprb, 2.29426e-06_jprb, 2.29306e-06_jprb, &
     & 2.29186e-06_jprb, 2.29066e-06_jprb, 2.28946e-06_jprb, 2.28826e-06_jprb, 2.28706e-06_jprb, &
     & 2.28586e-06_jprb, 2.28466e-06_jprb, 2.28347e-06_jprb, 2.28227e-06_jprb/)
      kbo_mo2(:,12) = (/ &
     & 2.38201e-06_jprb, 2.36536e-06_jprb, 2.34882e-06_jprb, 2.33240e-06_jprb, 2.31609e-06_jprb, &
     & 2.29990e-06_jprb, 2.28382e-06_jprb, 2.26785e-06_jprb, 2.25199e-06_jprb, 2.23625e-06_jprb, &
     & 2.22061e-06_jprb, 2.20508e-06_jprb, 2.18967e-06_jprb, 2.17436e-06_jprb, 2.15915e-06_jprb, &
     & 2.14406e-06_jprb, 2.12907e-06_jprb, 2.11418e-06_jprb, 2.09940e-06_jprb/)
      kbo_mo2(:,13) = (/ &
     & 2.33326e-06_jprb, 2.32549e-06_jprb, 2.31775e-06_jprb, 2.31003e-06_jprb, 2.30234e-06_jprb, &
     & 2.29467e-06_jprb, 2.28703e-06_jprb, 2.27941e-06_jprb, 2.27182e-06_jprb, 2.26426e-06_jprb, &
     & 2.25672e-06_jprb, 2.24920e-06_jprb, 2.24171e-06_jprb, 2.23424e-06_jprb, 2.22680e-06_jprb, &
     & 2.21939e-06_jprb, 2.21200e-06_jprb, 2.20463e-06_jprb, 2.19729e-06_jprb/)
      kbo_mo2(:,14) = (/ &
     & 2.75292e-06_jprb, 2.75210e-06_jprb, 2.75129e-06_jprb, 2.75047e-06_jprb, 2.74965e-06_jprb, &
     & 2.74883e-06_jprb, 2.74801e-06_jprb, 2.74720e-06_jprb, 2.74638e-06_jprb, 2.74556e-06_jprb, &
     & 2.74475e-06_jprb, 2.74393e-06_jprb, 2.74311e-06_jprb, 2.74230e-06_jprb, 2.74148e-06_jprb, &
     & 2.74067e-06_jprb, 2.73985e-06_jprb, 2.73904e-06_jprb, 2.73822e-06_jprb/)
      kbo_mo2(:,15) = (/ &
     & 2.55262e-06_jprb, 2.53364e-06_jprb, 2.51480e-06_jprb, 2.49611e-06_jprb, 2.47755e-06_jprb, &
     & 2.45913e-06_jprb, 2.44084e-06_jprb, 2.42269e-06_jprb, 2.40468e-06_jprb, 2.38680e-06_jprb, &
     & 2.36906e-06_jprb, 2.35144e-06_jprb, 2.33396e-06_jprb, 2.31660e-06_jprb, 2.29938e-06_jprb, &
     & 2.28228e-06_jprb, 2.26531e-06_jprb, 2.24847e-06_jprb, 2.23175e-06_jprb/)
      kbo_mo2(:,16) = (/ &
     & 3.11382e-06_jprb, 3.08751e-06_jprb, 3.06141e-06_jprb, 3.03554e-06_jprb, 3.00989e-06_jprb, &
     & 2.98445e-06_jprb, 2.95923e-06_jprb, 2.93422e-06_jprb, 2.90942e-06_jprb, 2.88483e-06_jprb, &
     & 2.86045e-06_jprb, 2.83628e-06_jprb, 2.81231e-06_jprb, 2.78854e-06_jprb, 2.76498e-06_jprb, &
     & 2.74161e-06_jprb, 2.71844e-06_jprb, 2.69547e-06_jprb, 2.67269e-06_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &2.8858e-02_jprb,3.6879e-02_jprb,4.0746e-02_jprb,4.2561e-02_jprb,4.2740e-02_jprb,4.2707e-02_jprb, &
     &4.4109e-02_jprb,4.4540e-02_jprb,4.5206e-02_jprb,4.4679e-02_jprb,4.5034e-02_jprb,4.5364e-02_jprb, &
     &4.6790e-02_jprb,4.7857e-02_jprb,4.8328e-02_jprb,4.8084e-02_jprb/)
      forrefo(2,:) = (/ &
     &2.7887e-02_jprb,3.7376e-02_jprb,4.0980e-02_jprb,4.2986e-02_jprb,4.3054e-02_jprb,4.2975e-02_jprb, &
     &4.3754e-02_jprb,4.4352e-02_jprb,4.4723e-02_jprb,4.6236e-02_jprb,4.5273e-02_jprb,4.5360e-02_jprb, &
     &4.5332e-02_jprb,4.7587e-02_jprb,4.7035e-02_jprb,5.0267e-02_jprb/)
      forrefo(3,:) = (/ &
     &2.5846e-02_jprb,3.6753e-02_jprb,4.2334e-02_jprb,4.3806e-02_jprb,4.3848e-02_jprb,4.3215e-02_jprb, &
     &4.3838e-02_jprb,4.4278e-02_jprb,4.4658e-02_jprb,4.5403e-02_jprb,4.5255e-02_jprb,4.6347e-02_jprb, &
     &4.4722e-02_jprb,4.6612e-02_jprb,4.6836e-02_jprb,4.8720e-02_jprb/)
      forrefo(4,:) = (/ &
     &2.8955e-02_jprb,3.7608e-02_jprb,4.1989e-02_jprb,4.4919e-02_jprb,4.2803e-02_jprb,4.2842e-02_jprb, &
     &4.2632e-02_jprb,4.1056e-02_jprb,4.0086e-02_jprb,4.1401e-02_jprb,4.2746e-02_jprb,4.2142e-02_jprb, &
     &4.1871e-02_jprb,4.3917e-02_jprb,4.5462e-02_jprb,4.8359e-02_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 5.96496e-01_jprb, 5.49171e-01_jprb, 5.05600e-01_jprb, 4.65486e-01_jprb, 4.28555e-01_jprb, &
     & 3.94554e-01_jprb, 3.63250e-01_jprb, 3.34430e-01_jprb, 3.07897e-01_jprb, 2.83468e-01_jprb/)
      selfrefo(:, 2) = (/ &
     & 7.46455e-01_jprb, 6.82459e-01_jprb, 6.23950e-01_jprb, 5.70457e-01_jprb, 5.21550e-01_jprb, &
     & 4.76836e-01_jprb, 4.35956e-01_jprb, 3.98580e-01_jprb, 3.64409e-01_jprb, 3.33167e-01_jprb/)
      selfrefo(:, 3) = (/ &
     & 7.86805e-01_jprb, 7.21186e-01_jprb, 6.61040e-01_jprb, 6.05910e-01_jprb, 5.55378e-01_jprb, &
     & 5.09059e-01_jprb, 4.66605e-01_jprb, 4.27690e-01_jprb, 3.92021e-01_jprb, 3.59327e-01_jprb/)
      selfrefo(:, 4) = (/ &
     & 8.11740e-01_jprb, 7.44359e-01_jprb, 6.82570e-01_jprb, 6.25910e-01_jprb, 5.73954e-01_jprb, &
     & 5.26311e-01_jprb, 4.82622e-01_jprb, 4.42560e-01_jprb, 4.05823e-01_jprb, 3.72136e-01_jprb/)
      selfrefo(:, 5) = (/ &
     & 8.14870e-01_jprb, 7.47200e-01_jprb, 6.85150e-01_jprb, 6.28253e-01_jprb, 5.76081e-01_jprb, &
     & 5.28241e-01_jprb, 4.84374e-01_jprb, 4.44150e-01_jprb, 4.07266e-01_jprb, 3.73446e-01_jprb/)
      selfrefo(:, 6) = (/ &
     & 8.10104e-01_jprb, 7.43259e-01_jprb, 6.81930e-01_jprb, 6.25661e-01_jprb, 5.74035e-01_jprb, &
     & 5.26669e-01_jprb, 4.83212e-01_jprb, 4.43340e-01_jprb, 4.06758e-01_jprb, 3.73195e-01_jprb/)
      selfrefo(:, 7) = (/ &
     & 8.13119e-01_jprb, 7.48127e-01_jprb, 6.88330e-01_jprb, 6.33312e-01_jprb, 5.82692e-01_jprb, &
     & 5.36118e-01_jprb, 4.93267e-01_jprb, 4.53840e-01_jprb, 4.17565e-01_jprb, 3.84189e-01_jprb/)
      selfrefo(:, 8) = (/ &
     & 8.26137e-01_jprb, 7.58984e-01_jprb, 6.97290e-01_jprb, 6.40611e-01_jprb, 5.88539e-01_jprb, &
     & 5.40699e-01_jprb, 4.96748e-01_jprb, 4.56370e-01_jprb, 4.19274e-01_jprb, 3.85193e-01_jprb/)
      selfrefo(:, 9) = (/ &
     & 8.30566e-01_jprb, 7.63984e-01_jprb, 7.02740e-01_jprb, 6.46405e-01_jprb, 5.94587e-01_jprb, &
     & 5.46922e-01_jprb, 5.03079e-01_jprb, 4.62750e-01_jprb, 4.25654e-01_jprb, 3.91532e-01_jprb/)
      selfrefo(:,10) = (/ &
     & 8.67471e-01_jprb, 7.91575e-01_jprb, 7.22320e-01_jprb, 6.59124e-01_jprb, 6.01457e-01_jprb, &
     & 5.48835e-01_jprb, 5.00817e-01_jprb, 4.57000e-01_jprb, 4.17017e-01_jprb, 3.80532e-01_jprb/)
      selfrefo(:,11) = (/ &
     & 8.51029e-01_jprb, 7.79373e-01_jprb, 7.13750e-01_jprb, 6.53652e-01_jprb, 5.98615e-01_jprb, &
     & 5.48212e-01_jprb, 5.02053e-01_jprb, 4.59780e-01_jprb, 4.21067e-01_jprb, 3.85613e-01_jprb/)
      selfrefo(:,12) = (/ &
     & 8.36772e-01_jprb, 7.68751e-01_jprb, 7.06260e-01_jprb, 6.48848e-01_jprb, 5.96104e-01_jprb, &
     & 5.47647e-01_jprb, 5.03129e-01_jprb, 4.62230e-01_jprb, 4.24655e-01_jprb, 3.90136e-01_jprb/)
      selfrefo(:,13) = (/ &
     & 8.36551e-01_jprb, 7.71089e-01_jprb, 7.10750e-01_jprb, 6.55133e-01_jprb, 6.03867e-01_jprb, &
     & 5.56614e-01_jprb, 5.13058e-01_jprb, 4.72910e-01_jprb, 4.35904e-01_jprb, 4.01794e-01_jprb/)
      selfrefo(:,14) = (/ &
     & 8.84307e-01_jprb, 8.11175e-01_jprb, 7.44090e-01_jprb, 6.82553e-01_jprb, 6.26106e-01_jprb, &
     & 5.74326e-01_jprb, 5.26829e-01_jprb, 4.83260e-01_jprb, 4.43294e-01_jprb, 4.06633e-01_jprb/)
      selfrefo(:,15) = (/ &
     & 8.90356e-01_jprb, 8.19830e-01_jprb, 7.54890e-01_jprb, 6.95094e-01_jprb, 6.40035e-01_jprb, &
     & 5.89337e-01_jprb, 5.42655e-01_jprb, 4.99670e-01_jprb, 4.60090e-01_jprb, 4.23646e-01_jprb/)
      selfrefo(:,16) = (/ &
     & 9.67549e-01_jprb, 8.79393e-01_jprb, 7.99270e-01_jprb, 7.26447e-01_jprb, 6.60259e-01_jprb, &
     & 6.00101e-01_jprb, 5.45425e-01_jprb, 4.95730e-01_jprb, 4.50563e-01_jprb, 4.09511e-01_jprb/)

if (lhook) call dr_hook('rrtm_kgb11',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb11:error reading file radrrtm")

end subroutine rrtm_kgb11
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

