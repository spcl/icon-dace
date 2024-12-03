! # 1 "ifsrrtm/rrtm_kgb15.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb15.f90"
subroutine rrtm_kgb15

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     abozzo 2001306 updated to rrtmg v4.85
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!      f. vana  05-mar-2015  support for single precision
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook,   dr_hook, jphook
use yomlun    ,only : nulrad
use mpl_module,only : mpl_broadcast
use yomtag    ,only : mtagrad
use yommp0    , only : nproc, myproc

use yoerrto15, only : kao,kao_mn2  ,selfrefo,forrefo   ,fracrefao, kao_d


!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 31 "ifsrrtm/rrtm_kgb15.f90" 2

if (lhook) call dr_hook('rrtm_kgb15',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d
  kao = real(kao_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb15:')
endif

! planck fraction mapping level : p = 1053. mb, t = 294.2 k
      fracrefao(:, 1) = (/ &
     &  1.0689e-01_jprb,1.1563e-01_jprb,1.2447e-01_jprb,1.2921e-01_jprb,1.2840e-01_jprb,1.2113e-01_jprb, &
     &  1.0643e-01_jprb,8.4987e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 2) = (/ &
     &  1.0782e-01_jprb,1.1637e-01_jprb,1.2290e-01_jprb,1.2911e-01_jprb,1.2841e-01_jprb,1.2113e-01_jprb, &
     &  1.0643e-01_jprb,8.4987e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 3) = (/ &
     &  1.0858e-01_jprb,1.1860e-01_jprb,1.2237e-01_jprb,1.2665e-01_jprb,1.2841e-01_jprb,1.2111e-01_jprb, &
     &  1.0642e-01_jprb,8.4987e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 4) = (/ &
     &  1.1022e-01_jprb,1.1965e-01_jprb,1.2334e-01_jprb,1.2383e-01_jprb,1.2761e-01_jprb,1.2109e-01_jprb, &
     &  1.0642e-01_jprb,8.4987e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 5) = (/ &
     &  1.1342e-01_jprb,1.2069e-01_jprb,1.2360e-01_jprb,1.2447e-01_jprb,1.2340e-01_jprb,1.2020e-01_jprb, &
     &  1.0639e-01_jprb,8.4987e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 6) = (/ &
     &  1.1771e-01_jprb,1.2280e-01_jprb,1.2177e-01_jprb,1.2672e-01_jprb,1.2398e-01_jprb,1.1787e-01_jprb, &
     &  1.0131e-01_jprb,8.4987e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 7) = (/ &
     &  1.2320e-01_jprb,1.2491e-01_jprb,1.2001e-01_jprb,1.2936e-01_jprb,1.2653e-01_jprb,1.1929e-01_jprb, &
     &  9.8955e-02_jprb,7.4887e-02_jprb,6.0142e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 8) = (/ &
     &  1.3105e-01_jprb,1.2563e-01_jprb,1.3055e-01_jprb,1.2854e-01_jprb,1.3402e-01_jprb,1.1571e-01_jprb, &
     &  9.4876e-02_jprb,6.0459e-02_jprb,5.6457e-02_jprb,6.6798e-03_jprb,5.5293e-03_jprb,4.3700e-03_jprb, &
     &  3.2061e-03_jprb,2.0476e-03_jprb,7.7366e-04_jprb,1.0897e-04_jprb/)
      fracrefao(:, 9) = (/ &
     &  1.1375e-01_jprb,1.2090e-01_jprb,1.2348e-01_jprb,1.2458e-01_jprb,1.2406e-01_jprb,1.1921e-01_jprb, &
     &  1.0802e-01_jprb,8.6613e-02_jprb,5.8125e-02_jprb,6.2984e-03_jprb,5.2359e-03_jprb,4.0641e-03_jprb, &
     &  2.9379e-03_jprb,1.9001e-03_jprb,7.2646e-04_jprb,1.0553e-04_jprb/)


!     ------------------------------------------------------------------

!     the array kao contains absorption coefs at the 16 chosen g-values 
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


!     the array ka_mxx contains the absorption coefficient for 
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

       kao_mn2( 1, :, 1) = (/ &
     & 3.24352e-08_jprb, 3.39625e-08_jprb, 3.55618e-08_jprb, 3.72364e-08_jprb, 3.89899e-08_jprb, &
     & 4.08259e-08_jprb, 4.27484e-08_jprb, 4.47614e-08_jprb, 4.68692e-08_jprb, 4.90763e-08_jprb, &
     & 5.13873e-08_jprb, 5.38071e-08_jprb, 5.63409e-08_jprb, 5.89940e-08_jprb, 6.17720e-08_jprb, &
     & 6.46808e-08_jprb, 6.77266e-08_jprb, 7.09158e-08_jprb, 7.42553e-08_jprb/)
       kao_mn2( 2, :, 1) = (/ &
     & 3.44203e-08_jprb, 3.60254e-08_jprb, 3.77053e-08_jprb, 3.94636e-08_jprb, 4.13038e-08_jprb, &
     & 4.32299e-08_jprb, 4.52458e-08_jprb, 4.73557e-08_jprb, 4.95640e-08_jprb, 5.18753e-08_jprb, &
     & 5.42943e-08_jprb, 5.68262e-08_jprb, 5.94761e-08_jprb, 6.22496e-08_jprb, 6.51524e-08_jprb, &
     & 6.81906e-08_jprb, 7.13704e-08_jprb, 7.46986e-08_jprb, 7.81819e-08_jprb/)
       kao_mn2( 3, :, 1) = (/ &
     & 3.44344e-08_jprb, 3.61485e-08_jprb, 3.79480e-08_jprb, 3.98370e-08_jprb, 4.18201e-08_jprb, &
     & 4.39019e-08_jprb, 4.60873e-08_jprb, 4.83815e-08_jprb, 5.07899e-08_jprb, 5.33182e-08_jprb, &
     & 5.59723e-08_jprb, 5.87586e-08_jprb, 6.16836e-08_jprb, 6.47541e-08_jprb, 6.79776e-08_jprb, &
     & 7.13614e-08_jprb, 7.49138e-08_jprb, 7.86429e-08_jprb, 8.25577e-08_jprb/)
       kao_mn2( 4, :, 1) = (/ &
     & 4.21102e-08_jprb, 4.38921e-08_jprb, 4.57493e-08_jprb, 4.76852e-08_jprb, 4.97029e-08_jprb, &
     & 5.18061e-08_jprb, 5.39982e-08_jprb, 5.62831e-08_jprb, 5.86647e-08_jprb, 6.11470e-08_jprb, &
     & 6.37344e-08_jprb, 6.64313e-08_jprb, 6.92422e-08_jprb, 7.21722e-08_jprb, 7.52261e-08_jprb, &
     & 7.84092e-08_jprb, 8.17270e-08_jprb, 8.51852e-08_jprb, 8.87897e-08_jprb/)
       kao_mn2( 5, :, 1) = (/ &
     & 4.78813e-08_jprb, 5.01015e-08_jprb, 5.24246e-08_jprb, 5.48554e-08_jprb, 5.73989e-08_jprb, &
     & 6.00603e-08_jprb, 6.28452e-08_jprb, 6.57592e-08_jprb, 6.88083e-08_jprb, 7.19987e-08_jprb, &
     & 7.53371e-08_jprb, 7.88304e-08_jprb, 8.24855e-08_jprb, 8.63102e-08_jprb, 9.03122e-08_jprb, &
     & 9.44997e-08_jprb, 9.88815e-08_jprb, 1.03466e-07_jprb, 1.08264e-07_jprb/)
       kao_mn2( 6, :, 1) = (/ &
     & 7.03115e-08_jprb, 7.27877e-08_jprb, 7.53511e-08_jprb, 7.80048e-08_jprb, 8.07519e-08_jprb, &
     & 8.35958e-08_jprb, 8.65398e-08_jprb, 8.95875e-08_jprb, 9.27426e-08_jprb, 9.60087e-08_jprb, &
     & 9.93899e-08_jprb, 1.02890e-07_jprb, 1.06514e-07_jprb, 1.10265e-07_jprb, 1.14148e-07_jprb, &
     & 1.18168e-07_jprb, 1.22330e-07_jprb, 1.26638e-07_jprb, 1.31098e-07_jprb/)
       kao_mn2( 7, :, 1) = (/ &
     & 8.86454e-08_jprb, 9.20065e-08_jprb, 9.54951e-08_jprb, 9.91159e-08_jprb, 1.02874e-07_jprb, &
     & 1.06775e-07_jprb, 1.10823e-07_jprb, 1.15025e-07_jprb, 1.19387e-07_jprb, 1.23913e-07_jprb, &
     & 1.28612e-07_jprb, 1.33488e-07_jprb, 1.38550e-07_jprb, 1.43803e-07_jprb, 1.49255e-07_jprb, &
     & 1.54915e-07_jprb, 1.60788e-07_jprb, 1.66885e-07_jprb, 1.73213e-07_jprb/)
       kao_mn2( 8, :, 1) = (/ &
     & 1.34118e-07_jprb, 1.38267e-07_jprb, 1.42545e-07_jprb, 1.46955e-07_jprb, 1.51502e-07_jprb, &
     & 1.56189e-07_jprb, 1.61022e-07_jprb, 1.66004e-07_jprb, 1.71140e-07_jprb, 1.76435e-07_jprb, &
     & 1.81893e-07_jprb, 1.87521e-07_jprb, 1.93323e-07_jprb, 1.99304e-07_jprb, 2.05470e-07_jprb, &
     & 2.11827e-07_jprb, 2.18381e-07_jprb, 2.25138e-07_jprb, 2.32103e-07_jprb/)
       kao_mn2( 9, :, 1) = (/ &
     & 5.08256e-08_jprb, 5.30384e-08_jprb, 5.53476e-08_jprb, 5.77573e-08_jprb, 6.02718e-08_jprb, &
     & 6.28959e-08_jprb, 6.56342e-08_jprb, 6.84917e-08_jprb, 7.14737e-08_jprb, 7.45854e-08_jprb, &
     & 7.78327e-08_jprb, 8.12213e-08_jprb, 8.47574e-08_jprb, 8.84475e-08_jprb, 9.22983e-08_jprb, &
     & 9.63167e-08_jprb, 1.00510e-07_jprb, 1.04886e-07_jprb, 1.09452e-07_jprb/)
       kao_mn2( 1, :, 2) = (/ &
     & 8.23958e-08_jprb, 8.39092e-08_jprb, 8.54504e-08_jprb, 8.70200e-08_jprb, 8.86183e-08_jprb, &
     & 9.02460e-08_jprb, 9.19036e-08_jprb, 9.35917e-08_jprb, 9.53107e-08_jprb, 9.70614e-08_jprb, &
     & 9.88442e-08_jprb, 1.00660e-07_jprb, 1.02509e-07_jprb, 1.04391e-07_jprb, 1.06309e-07_jprb, &
     & 1.08261e-07_jprb, 1.10250e-07_jprb, 1.12275e-07_jprb, 1.14337e-07_jprb/)
       kao_mn2( 2, :, 2) = (/ &
     & 8.52335e-08_jprb, 8.69254e-08_jprb, 8.86509e-08_jprb, 9.04107e-08_jprb, 9.22054e-08_jprb, &
     & 9.40357e-08_jprb, 9.59024e-08_jprb, 9.78061e-08_jprb, 9.97476e-08_jprb, 1.01728e-07_jprb, &
     & 1.03747e-07_jprb, 1.05806e-07_jprb, 1.07907e-07_jprb, 1.10049e-07_jprb, 1.12233e-07_jprb, &
     & 1.14461e-07_jprb, 1.16733e-07_jprb, 1.19050e-07_jprb, 1.21414e-07_jprb/)
       kao_mn2( 3, :, 2) = (/ &
     & 1.04608e-07_jprb, 1.06067e-07_jprb, 1.07546e-07_jprb, 1.09046e-07_jprb, 1.10567e-07_jprb, &
     & 1.12110e-07_jprb, 1.13673e-07_jprb, 1.15259e-07_jprb, 1.16866e-07_jprb, 1.18496e-07_jprb, &
     & 1.20149e-07_jprb, 1.21825e-07_jprb, 1.23524e-07_jprb, 1.25247e-07_jprb, 1.26994e-07_jprb, &
     & 1.28765e-07_jprb, 1.30561e-07_jprb, 1.32382e-07_jprb, 1.34229e-07_jprb/)
       kao_mn2( 4, :, 2) = (/ &
     & 1.17504e-07_jprb, 1.18763e-07_jprb, 1.20036e-07_jprb, 1.21322e-07_jprb, 1.22622e-07_jprb, &
     & 1.23936e-07_jprb, 1.25265e-07_jprb, 1.26607e-07_jprb, 1.27964e-07_jprb, 1.29335e-07_jprb, &
     & 1.30721e-07_jprb, 1.32122e-07_jprb, 1.33538e-07_jprb, 1.34969e-07_jprb, 1.36415e-07_jprb, &
     & 1.37877e-07_jprb, 1.39354e-07_jprb, 1.40848e-07_jprb, 1.42357e-07_jprb/)
       kao_mn2( 5, :, 2) = (/ &
     & 1.23552e-07_jprb, 1.25200e-07_jprb, 1.26870e-07_jprb, 1.28562e-07_jprb, 1.30277e-07_jprb, &
     & 1.32015e-07_jprb, 1.33776e-07_jprb, 1.35560e-07_jprb, 1.37368e-07_jprb, 1.39200e-07_jprb, &
     & 1.41057e-07_jprb, 1.42938e-07_jprb, 1.44845e-07_jprb, 1.46777e-07_jprb, 1.48735e-07_jprb, &
     & 1.50718e-07_jprb, 1.52729e-07_jprb, 1.54766e-07_jprb, 1.56830e-07_jprb/)
       kao_mn2( 6, :, 2) = (/ &
     & 1.29682e-07_jprb, 1.32226e-07_jprb, 1.34820e-07_jprb, 1.37464e-07_jprb, 1.40161e-07_jprb, &
     & 1.42910e-07_jprb, 1.45713e-07_jprb, 1.48571e-07_jprb, 1.51486e-07_jprb, 1.54457e-07_jprb, &
     & 1.57487e-07_jprb, 1.60576e-07_jprb, 1.63726e-07_jprb, 1.66937e-07_jprb, 1.70212e-07_jprb, &
     & 1.73551e-07_jprb, 1.76955e-07_jprb, 1.80426e-07_jprb, 1.83965e-07_jprb/)
       kao_mn2( 7, :, 2) = (/ &
     & 1.77416e-07_jprb, 1.78627e-07_jprb, 1.79846e-07_jprb, 1.81073e-07_jprb, 1.82309e-07_jprb, &
     & 1.83554e-07_jprb, 1.84806e-07_jprb, 1.86068e-07_jprb, 1.87338e-07_jprb, 1.88616e-07_jprb, &
     & 1.89904e-07_jprb, 1.91200e-07_jprb, 1.92505e-07_jprb, 1.93819e-07_jprb, 1.95142e-07_jprb, &
     & 1.96474e-07_jprb, 1.97815e-07_jprb, 1.99165e-07_jprb, 2.00524e-07_jprb/)
       kao_mn2( 8, :, 2) = (/ &
     & 2.20695e-07_jprb, 2.20451e-07_jprb, 2.20208e-07_jprb, 2.19965e-07_jprb, 2.19722e-07_jprb, &
     & 2.19480e-07_jprb, 2.19238e-07_jprb, 2.18996e-07_jprb, 2.18754e-07_jprb, 2.18513e-07_jprb, &
     & 2.18272e-07_jprb, 2.18031e-07_jprb, 2.17790e-07_jprb, 2.17550e-07_jprb, 2.17310e-07_jprb, &
     & 2.17070e-07_jprb, 2.16831e-07_jprb, 2.16591e-07_jprb, 2.16352e-07_jprb/)
       kao_mn2( 9, :, 2) = (/ &
     & 1.23015e-07_jprb, 1.24808e-07_jprb, 1.26626e-07_jprb, 1.28471e-07_jprb, 1.30343e-07_jprb, &
     & 1.32242e-07_jprb, 1.34168e-07_jprb, 1.36123e-07_jprb, 1.38106e-07_jprb, 1.40118e-07_jprb, &
     & 1.42160e-07_jprb, 1.44231e-07_jprb, 1.46332e-07_jprb, 1.48464e-07_jprb, 1.50627e-07_jprb, &
     & 1.52822e-07_jprb, 1.55048e-07_jprb, 1.57307e-07_jprb, 1.59599e-07_jprb/)
       kao_mn2( 1, :, 3) = (/ &
     & 1.87585e-07_jprb, 1.89503e-07_jprb, 1.91440e-07_jprb, 1.93398e-07_jprb, 1.95375e-07_jprb, &
     & 1.97372e-07_jprb, 1.99390e-07_jprb, 2.01429e-07_jprb, 2.03488e-07_jprb, 2.05568e-07_jprb, &
     & 2.07670e-07_jprb, 2.09793e-07_jprb, 2.11938e-07_jprb, 2.14105e-07_jprb, 2.16294e-07_jprb, &
     & 2.18505e-07_jprb, 2.20739e-07_jprb, 2.22996e-07_jprb, 2.25275e-07_jprb/)
       kao_mn2( 2, :, 3) = (/ &
     & 1.82585e-07_jprb, 1.84249e-07_jprb, 1.85929e-07_jprb, 1.87624e-07_jprb, 1.89335e-07_jprb, &
     & 1.91061e-07_jprb, 1.92803e-07_jprb, 1.94561e-07_jprb, 1.96335e-07_jprb, 1.98125e-07_jprb, &
     & 1.99932e-07_jprb, 2.01755e-07_jprb, 2.03594e-07_jprb, 2.05451e-07_jprb, 2.07324e-07_jprb, &
     & 2.09214e-07_jprb, 2.11122e-07_jprb, 2.13047e-07_jprb, 2.14989e-07_jprb/)
       kao_mn2( 3, :, 3) = (/ &
     & 1.64711e-07_jprb, 1.67539e-07_jprb, 1.70417e-07_jprb, 1.73343e-07_jprb, 1.76321e-07_jprb, &
     & 1.79349e-07_jprb, 1.82429e-07_jprb, 1.85562e-07_jprb, 1.88749e-07_jprb, 1.91990e-07_jprb, &
     & 1.95288e-07_jprb, 1.98642e-07_jprb, 2.02053e-07_jprb, 2.05523e-07_jprb, 2.09053e-07_jprb, &
     & 2.12643e-07_jprb, 2.16295e-07_jprb, 2.20010e-07_jprb, 2.23788e-07_jprb/)
       kao_mn2( 4, :, 3) = (/ &
     & 1.67494e-07_jprb, 1.71011e-07_jprb, 1.74601e-07_jprb, 1.78267e-07_jprb, 1.82009e-07_jprb, &
     & 1.85831e-07_jprb, 1.89732e-07_jprb, 1.93715e-07_jprb, 1.97782e-07_jprb, 2.01935e-07_jprb, &
     & 2.06174e-07_jprb, 2.10503e-07_jprb, 2.14922e-07_jprb, 2.19434e-07_jprb, 2.24041e-07_jprb, &
     & 2.28745e-07_jprb, 2.33548e-07_jprb, 2.38451e-07_jprb, 2.43457e-07_jprb/)
       kao_mn2( 5, :, 3) = (/ &
     & 1.97399e-07_jprb, 2.00092e-07_jprb, 2.02821e-07_jprb, 2.05588e-07_jprb, 2.08393e-07_jprb, &
     & 2.11236e-07_jprb, 2.14118e-07_jprb, 2.17039e-07_jprb, 2.20000e-07_jprb, 2.23001e-07_jprb, &
     & 2.26043e-07_jprb, 2.29127e-07_jprb, 2.32252e-07_jprb, 2.35421e-07_jprb, 2.38633e-07_jprb, &
     & 2.41888e-07_jprb, 2.45188e-07_jprb, 2.48533e-07_jprb, 2.51923e-07_jprb/)
       kao_mn2( 6, :, 3) = (/ &
     & 2.24021e-07_jprb, 2.24970e-07_jprb, 2.25923e-07_jprb, 2.26880e-07_jprb, 2.27840e-07_jprb, &
     & 2.28805e-07_jprb, 2.29774e-07_jprb, 2.30747e-07_jprb, 2.31725e-07_jprb, 2.32706e-07_jprb, &
     & 2.33692e-07_jprb, 2.34681e-07_jprb, 2.35675e-07_jprb, 2.36673e-07_jprb, 2.37675e-07_jprb, &
     & 2.38682e-07_jprb, 2.39693e-07_jprb, 2.40708e-07_jprb, 2.41727e-07_jprb/)
       kao_mn2( 7, :, 3) = (/ &
     & 1.98178e-07_jprb, 2.00676e-07_jprb, 2.03205e-07_jprb, 2.05766e-07_jprb, 2.08359e-07_jprb, &
     & 2.10986e-07_jprb, 2.13645e-07_jprb, 2.16337e-07_jprb, 2.19064e-07_jprb, 2.21825e-07_jprb, &
     & 2.24621e-07_jprb, 2.27452e-07_jprb, 2.30319e-07_jprb, 2.33222e-07_jprb, 2.36161e-07_jprb, &
     & 2.39138e-07_jprb, 2.42152e-07_jprb, 2.45204e-07_jprb, 2.48294e-07_jprb/)
       kao_mn2( 8, :, 3) = (/ &
     & 2.83042e-07_jprb, 2.89941e-07_jprb, 2.97009e-07_jprb, 3.04250e-07_jprb, 3.11666e-07_jprb, &
     & 3.19264e-07_jprb, 3.27047e-07_jprb, 3.35019e-07_jprb, 3.43186e-07_jprb, 3.51552e-07_jprb, &
     & 3.60122e-07_jprb, 3.68901e-07_jprb, 3.77893e-07_jprb, 3.87105e-07_jprb, 3.96542e-07_jprb, &
     & 4.06208e-07_jprb, 4.16111e-07_jprb, 4.26254e-07_jprb, 4.36645e-07_jprb/)
       kao_mn2( 9, :, 3) = (/ &
     & 1.98963e-07_jprb, 2.01576e-07_jprb, 2.04224e-07_jprb, 2.06907e-07_jprb, 2.09626e-07_jprb, &
     & 2.12379e-07_jprb, 2.15169e-07_jprb, 2.17996e-07_jprb, 2.20860e-07_jprb, 2.23761e-07_jprb, &
     & 2.26701e-07_jprb, 2.29679e-07_jprb, 2.32696e-07_jprb, 2.35753e-07_jprb, 2.38851e-07_jprb, &
     & 2.41988e-07_jprb, 2.45167e-07_jprb, 2.48388e-07_jprb, 2.51651e-07_jprb/)
       kao_mn2( 1, :, 4) = (/ &
     & 3.75434e-07_jprb, 3.79581e-07_jprb, 3.83775e-07_jprb, 3.88014e-07_jprb, 3.92301e-07_jprb, &
     & 3.96634e-07_jprb, 4.01016e-07_jprb, 4.05446e-07_jprb, 4.09925e-07_jprb, 4.14453e-07_jprb, &
     & 4.19032e-07_jprb, 4.23661e-07_jprb, 4.28341e-07_jprb, 4.33073e-07_jprb, 4.37857e-07_jprb, &
     & 4.42694e-07_jprb, 4.47585e-07_jprb, 4.52529e-07_jprb, 4.57528e-07_jprb/)
       kao_mn2( 2, :, 4) = (/ &
     & 3.76756e-07_jprb, 3.80760e-07_jprb, 3.84805e-07_jprb, 3.88894e-07_jprb, 3.93027e-07_jprb, &
     & 3.97203e-07_jprb, 4.01423e-07_jprb, 4.05689e-07_jprb, 4.10000e-07_jprb, 4.14356e-07_jprb, &
     & 4.18759e-07_jprb, 4.23209e-07_jprb, 4.27706e-07_jprb, 4.32250e-07_jprb, 4.36843e-07_jprb, &
     & 4.41485e-07_jprb, 4.46176e-07_jprb, 4.50917e-07_jprb, 4.55708e-07_jprb/)
       kao_mn2( 3, :, 4) = (/ &
     & 3.76258e-07_jprb, 3.78929e-07_jprb, 3.81619e-07_jprb, 3.84329e-07_jprb, 3.87057e-07_jprb, &
     & 3.89805e-07_jprb, 3.92572e-07_jprb, 3.95359e-07_jprb, 3.98166e-07_jprb, 4.00993e-07_jprb, &
     & 4.03839e-07_jprb, 4.06706e-07_jprb, 4.09594e-07_jprb, 4.12502e-07_jprb, 4.15430e-07_jprb, &
     & 4.18379e-07_jprb, 4.21349e-07_jprb, 4.24341e-07_jprb, 4.27353e-07_jprb/)
       kao_mn2( 4, :, 4) = (/ &
     & 3.17796e-07_jprb, 3.22447e-07_jprb, 3.27166e-07_jprb, 3.31954e-07_jprb, 3.36812e-07_jprb, &
     & 3.41742e-07_jprb, 3.46743e-07_jprb, 3.51818e-07_jprb, 3.56967e-07_jprb, 3.62191e-07_jprb, &
     & 3.67492e-07_jprb, 3.72870e-07_jprb, 3.78328e-07_jprb, 3.83865e-07_jprb, 3.89483e-07_jprb, &
     & 3.95183e-07_jprb, 4.00967e-07_jprb, 4.06835e-07_jprb, 4.12789e-07_jprb/)
       kao_mn2( 5, :, 4) = (/ &
     & 3.33793e-07_jprb, 3.38941e-07_jprb, 3.44169e-07_jprb, 3.49478e-07_jprb, 3.54868e-07_jprb, &
     & 3.60342e-07_jprb, 3.65900e-07_jprb, 3.71544e-07_jprb, 3.77275e-07_jprb, 3.83094e-07_jprb, &
     & 3.89003e-07_jprb, 3.95003e-07_jprb, 4.01096e-07_jprb, 4.07283e-07_jprb, 4.13565e-07_jprb, &
     & 4.19944e-07_jprb, 4.26421e-07_jprb, 4.32999e-07_jprb, 4.39677e-07_jprb/)
       kao_mn2( 6, :, 4) = (/ &
     & 3.60052e-07_jprb, 3.66686e-07_jprb, 3.73442e-07_jprb, 3.80323e-07_jprb, 3.87330e-07_jprb, &
     & 3.94466e-07_jprb, 4.01734e-07_jprb, 4.09136e-07_jprb, 4.16674e-07_jprb, 4.24351e-07_jprb, &
     & 4.32169e-07_jprb, 4.40132e-07_jprb, 4.48241e-07_jprb, 4.56500e-07_jprb, 4.64910e-07_jprb, &
     & 4.73476e-07_jprb, 4.82200e-07_jprb, 4.91084e-07_jprb, 5.00132e-07_jprb/)
       kao_mn2( 7, :, 4) = (/ &
     & 4.14713e-07_jprb, 4.21885e-07_jprb, 4.29181e-07_jprb, 4.36603e-07_jprb, 4.44153e-07_jprb, &
     & 4.51834e-07_jprb, 4.59648e-07_jprb, 4.67598e-07_jprb, 4.75684e-07_jprb, 4.83910e-07_jprb, &
     & 4.92279e-07_jprb, 5.00793e-07_jprb, 5.09453e-07_jprb, 5.18264e-07_jprb, 5.27226e-07_jprb, &
     & 5.36344e-07_jprb, 5.45620e-07_jprb, 5.55055e-07_jprb, 5.64654e-07_jprb/)
       kao_mn2( 8, :, 4) = (/ &
     & 4.15352e-07_jprb, 4.24386e-07_jprb, 4.33617e-07_jprb, 4.43049e-07_jprb, 4.52685e-07_jprb, &
     & 4.62532e-07_jprb, 4.72592e-07_jprb, 4.82872e-07_jprb, 4.93374e-07_jprb, 5.04106e-07_jprb, &
     & 5.15071e-07_jprb, 5.26274e-07_jprb, 5.37721e-07_jprb, 5.49417e-07_jprb, 5.61367e-07_jprb, &
     & 5.73577e-07_jprb, 5.86053e-07_jprb, 5.98800e-07_jprb, 6.11825e-07_jprb/)
       kao_mn2( 9, :, 4) = (/ &
     & 3.33820e-07_jprb, 3.39144e-07_jprb, 3.44553e-07_jprb, 3.50048e-07_jprb, 3.55631e-07_jprb, &
     & 3.61302e-07_jprb, 3.67065e-07_jprb, 3.72919e-07_jprb, 3.78866e-07_jprb, 3.84908e-07_jprb, &
     & 3.91047e-07_jprb, 3.97284e-07_jprb, 4.03620e-07_jprb, 4.10057e-07_jprb, 4.16597e-07_jprb, &
     & 4.23241e-07_jprb, 4.29991e-07_jprb, 4.36849e-07_jprb, 4.43816e-07_jprb/)
       kao_mn2( 1, :, 5) = (/ &
     & 6.99819e-07_jprb, 7.04629e-07_jprb, 7.09472e-07_jprb, 7.14349e-07_jprb, 7.19258e-07_jprb, &
     & 7.24202e-07_jprb, 7.29180e-07_jprb, 7.34192e-07_jprb, 7.39238e-07_jprb, 7.44319e-07_jprb, &
     & 7.49435e-07_jprb, 7.54586e-07_jprb, 7.59773e-07_jprb, 7.64995e-07_jprb, 7.70253e-07_jprb, &
     & 7.75547e-07_jprb, 7.80877e-07_jprb, 7.86245e-07_jprb, 7.91649e-07_jprb/)
       kao_mn2( 2, :, 5) = (/ &
     & 6.98257e-07_jprb, 7.03182e-07_jprb, 7.08143e-07_jprb, 7.13138e-07_jprb, 7.18169e-07_jprb, &
     & 7.23235e-07_jprb, 7.28336e-07_jprb, 7.33474e-07_jprb, 7.38648e-07_jprb, 7.43858e-07_jprb, &
     & 7.49106e-07_jprb, 7.54390e-07_jprb, 7.59711e-07_jprb, 7.65071e-07_jprb, 7.70467e-07_jprb, &
     & 7.75902e-07_jprb, 7.81376e-07_jprb, 7.86887e-07_jprb, 7.92438e-07_jprb/)
       kao_mn2( 3, :, 5) = (/ &
     & 6.98531e-07_jprb, 7.03429e-07_jprb, 7.08361e-07_jprb, 7.13328e-07_jprb, 7.18329e-07_jprb, &
     & 7.23365e-07_jprb, 7.28437e-07_jprb, 7.33545e-07_jprb, 7.38688e-07_jprb, 7.43867e-07_jprb, &
     & 7.49082e-07_jprb, 7.54335e-07_jprb, 7.59623e-07_jprb, 7.64950e-07_jprb, 7.70313e-07_jprb, &
     & 7.75714e-07_jprb, 7.81153e-07_jprb, 7.86630e-07_jprb, 7.92145e-07_jprb/)
       kao_mn2( 4, :, 5) = (/ &
     & 7.37210e-07_jprb, 7.38869e-07_jprb, 7.40532e-07_jprb, 7.42198e-07_jprb, 7.43868e-07_jprb, &
     & 7.45542e-07_jprb, 7.47219e-07_jprb, 7.48901e-07_jprb, 7.50586e-07_jprb, 7.52275e-07_jprb, &
     & 7.53967e-07_jprb, 7.55664e-07_jprb, 7.57364e-07_jprb, 7.59068e-07_jprb, 7.60777e-07_jprb, &
     & 7.62488e-07_jprb, 7.64204e-07_jprb, 7.65924e-07_jprb, 7.67647e-07_jprb/)
       kao_mn2( 5, :, 5) = (/ &
     & 6.07063e-07_jprb, 6.12893e-07_jprb, 6.18779e-07_jprb, 6.24722e-07_jprb, 6.30721e-07_jprb, &
     & 6.36778e-07_jprb, 6.42893e-07_jprb, 6.49067e-07_jprb, 6.55301e-07_jprb, 6.61594e-07_jprb, &
     & 6.67947e-07_jprb, 6.74362e-07_jprb, 6.80838e-07_jprb, 6.87376e-07_jprb, 6.93978e-07_jprb, &
     & 7.00642e-07_jprb, 7.07371e-07_jprb, 7.14164e-07_jprb, 7.21022e-07_jprb/)
       kao_mn2( 6, :, 5) = (/ &
     & 6.13354e-07_jprb, 6.20147e-07_jprb, 6.27016e-07_jprb, 6.33961e-07_jprb, 6.40983e-07_jprb, &
     & 6.48082e-07_jprb, 6.55260e-07_jprb, 6.62518e-07_jprb, 6.69856e-07_jprb, 6.77276e-07_jprb, &
     & 6.84777e-07_jprb, 6.92362e-07_jprb, 7.00030e-07_jprb, 7.07784e-07_jprb, 7.15624e-07_jprb, &
     & 7.23550e-07_jprb, 7.31564e-07_jprb, 7.39667e-07_jprb, 7.47859e-07_jprb/)
       kao_mn2( 7, :, 5) = (/ &
     & 6.86666e-07_jprb, 6.92902e-07_jprb, 6.99195e-07_jprb, 7.05545e-07_jprb, 7.11952e-07_jprb, &
     & 7.18418e-07_jprb, 7.24943e-07_jprb, 7.31526e-07_jprb, 7.38170e-07_jprb, 7.44874e-07_jprb, &
     & 7.51639e-07_jprb, 7.58465e-07_jprb, 7.65353e-07_jprb, 7.72304e-07_jprb, 7.79318e-07_jprb, &
     & 7.86395e-07_jprb, 7.93537e-07_jprb, 8.00744e-07_jprb, 8.08016e-07_jprb/)
       kao_mn2( 8, :, 5) = (/ &
     & 9.39664e-07_jprb, 9.39765e-07_jprb, 9.39866e-07_jprb, 9.39966e-07_jprb, 9.40067e-07_jprb, &
     & 9.40168e-07_jprb, 9.40268e-07_jprb, 9.40369e-07_jprb, 9.40470e-07_jprb, 9.40571e-07_jprb, &
     & 9.40671e-07_jprb, 9.40772e-07_jprb, 9.40873e-07_jprb, 9.40974e-07_jprb, 9.41074e-07_jprb, &
     & 9.41175e-07_jprb, 9.41276e-07_jprb, 9.41377e-07_jprb, 9.41478e-07_jprb/)
       kao_mn2( 9, :, 5) = (/ &
     & 6.02847e-07_jprb, 6.09726e-07_jprb, 6.16684e-07_jprb, 6.23722e-07_jprb, 6.30839e-07_jprb, &
     & 6.38038e-07_jprb, 6.45320e-07_jprb, 6.52684e-07_jprb, 6.60132e-07_jprb, 6.67665e-07_jprb, &
     & 6.75284e-07_jprb, 6.82991e-07_jprb, 6.90785e-07_jprb, 6.98668e-07_jprb, 7.06641e-07_jprb, &
     & 7.14705e-07_jprb, 7.22861e-07_jprb, 7.31110e-07_jprb, 7.39453e-07_jprb/)
       kao_mn2( 1, :, 6) = (/ &
     & 1.13692e-06_jprb, 1.13231e-06_jprb, 1.12772e-06_jprb, 1.12315e-06_jprb, 1.11859e-06_jprb, &
     & 1.11406e-06_jprb, 1.10954e-06_jprb, 1.10505e-06_jprb, 1.10057e-06_jprb, 1.09610e-06_jprb, &
     & 1.09166e-06_jprb, 1.08724e-06_jprb, 1.08283e-06_jprb, 1.07844e-06_jprb, 1.07407e-06_jprb, &
     & 1.06971e-06_jprb, 1.06538e-06_jprb, 1.06106e-06_jprb, 1.05676e-06_jprb/)
       kao_mn2( 2, :, 6) = (/ &
     & 1.13682e-06_jprb, 1.13221e-06_jprb, 1.12762e-06_jprb, 1.12305e-06_jprb, 1.11849e-06_jprb, &
     & 1.11396e-06_jprb, 1.10944e-06_jprb, 1.10495e-06_jprb, 1.10047e-06_jprb, 1.09600e-06_jprb, &
     & 1.09156e-06_jprb, 1.08714e-06_jprb, 1.08273e-06_jprb, 1.07834e-06_jprb, 1.07397e-06_jprb, &
     & 1.06961e-06_jprb, 1.06528e-06_jprb, 1.06096e-06_jprb, 1.05666e-06_jprb/)
       kao_mn2( 3, :, 6) = (/ &
     & 1.13642e-06_jprb, 1.13181e-06_jprb, 1.12722e-06_jprb, 1.12265e-06_jprb, 1.11809e-06_jprb, &
     & 1.11356e-06_jprb, 1.10904e-06_jprb, 1.10455e-06_jprb, 1.10007e-06_jprb, 1.09560e-06_jprb, &
     & 1.09116e-06_jprb, 1.08674e-06_jprb, 1.08233e-06_jprb, 1.07794e-06_jprb, 1.07357e-06_jprb, &
     & 1.06921e-06_jprb, 1.06488e-06_jprb, 1.06056e-06_jprb, 1.05626e-06_jprb/)
       kao_mn2( 4, :, 6) = (/ &
     & 1.13626e-06_jprb, 1.13160e-06_jprb, 1.12696e-06_jprb, 1.12233e-06_jprb, 1.11773e-06_jprb, &
     & 1.11314e-06_jprb, 1.10858e-06_jprb, 1.10403e-06_jprb, 1.09950e-06_jprb, 1.09498e-06_jprb, &
     & 1.09049e-06_jprb, 1.08602e-06_jprb, 1.08156e-06_jprb, 1.07712e-06_jprb, 1.07270e-06_jprb, &
     & 1.06830e-06_jprb, 1.06392e-06_jprb, 1.05955e-06_jprb, 1.05520e-06_jprb/)
       kao_mn2( 5, :, 6) = (/ &
     & 1.22429e-06_jprb, 1.21163e-06_jprb, 1.19909e-06_jprb, 1.18669e-06_jprb, 1.17441e-06_jprb, &
     & 1.16226e-06_jprb, 1.15024e-06_jprb, 1.13834e-06_jprb, 1.12656e-06_jprb, 1.11491e-06_jprb, &
     & 1.10338e-06_jprb, 1.09196e-06_jprb, 1.08067e-06_jprb, 1.06949e-06_jprb, 1.05842e-06_jprb, &
     & 1.04747e-06_jprb, 1.03664e-06_jprb, 1.02591e-06_jprb, 1.01530e-06_jprb/)
       kao_mn2( 6, :, 6) = (/ &
     & 1.02400e-06_jprb, 1.02238e-06_jprb, 1.02077e-06_jprb, 1.01916e-06_jprb, 1.01755e-06_jprb, &
     & 1.01594e-06_jprb, 1.01433e-06_jprb, 1.01273e-06_jprb, 1.01113e-06_jprb, 1.00953e-06_jprb, &
     & 1.00794e-06_jprb, 1.00635e-06_jprb, 1.00476e-06_jprb, 1.00317e-06_jprb, 1.00159e-06_jprb, &
     & 1.00000e-06_jprb, 9.98425e-07_jprb, 9.96848e-07_jprb, 9.95273e-07_jprb/)
       kao_mn2( 7, :, 6) = (/ &
     & 1.08594e-06_jprb, 1.08185e-06_jprb, 1.07778e-06_jprb, 1.07373e-06_jprb, 1.06969e-06_jprb, &
     & 1.06566e-06_jprb, 1.06165e-06_jprb, 1.05766e-06_jprb, 1.05368e-06_jprb, 1.04971e-06_jprb, &
     & 1.04576e-06_jprb, 1.04183e-06_jprb, 1.03791e-06_jprb, 1.03400e-06_jprb, 1.03011e-06_jprb, &
     & 1.02623e-06_jprb, 1.02237e-06_jprb, 1.01852e-06_jprb, 1.01469e-06_jprb/)
       kao_mn2( 8, :, 6) = (/ &
     & 1.25029e-06_jprb, 1.22508e-06_jprb, 1.20038e-06_jprb, 1.17618e-06_jprb, 1.15247e-06_jprb, &
     & 1.12924e-06_jprb, 1.10647e-06_jprb, 1.08416e-06_jprb, 1.06231e-06_jprb, 1.04089e-06_jprb, &
     & 1.01990e-06_jprb, 9.99343e-07_jprb, 9.79196e-07_jprb, 9.59454e-07_jprb, 9.40111e-07_jprb, &
     & 9.21158e-07_jprb, 9.02587e-07_jprb, 8.84390e-07_jprb, 8.66560e-07_jprb/)
       kao_mn2( 9, :, 6) = (/ &
     & 1.21299e-06_jprb, 1.19953e-06_jprb, 1.18622e-06_jprb, 1.17305e-06_jprb, 1.16003e-06_jprb, &
     & 1.14716e-06_jprb, 1.13443e-06_jprb, 1.12184e-06_jprb, 1.10939e-06_jprb, 1.09708e-06_jprb, &
     & 1.08491e-06_jprb, 1.07287e-06_jprb, 1.06096e-06_jprb, 1.04919e-06_jprb, 1.03755e-06_jprb, &
     & 1.02603e-06_jprb, 1.01465e-06_jprb, 1.00339e-06_jprb, 9.92253e-07_jprb/)
       kao_mn2( 1, :, 7) = (/ &
     & 1.53893e-06_jprb, 1.51743e-06_jprb, 1.49623e-06_jprb, 1.47532e-06_jprb, 1.45471e-06_jprb, &
     & 1.43438e-06_jprb, 1.41434e-06_jprb, 1.39458e-06_jprb, 1.37509e-06_jprb, 1.35588e-06_jprb, &
     & 1.33694e-06_jprb, 1.31826e-06_jprb, 1.29984e-06_jprb, 1.28167e-06_jprb, 1.26377e-06_jprb, &
     & 1.24611e-06_jprb, 1.22870e-06_jprb, 1.21153e-06_jprb, 1.19460e-06_jprb/)
       kao_mn2( 2, :, 7) = (/ &
     & 1.53809e-06_jprb, 1.51665e-06_jprb, 1.49552e-06_jprb, 1.47467e-06_jprb, 1.45412e-06_jprb, &
     & 1.43386e-06_jprb, 1.41388e-06_jprb, 1.39418e-06_jprb, 1.37475e-06_jprb, 1.35559e-06_jprb, &
     & 1.33670e-06_jprb, 1.31807e-06_jprb, 1.29970e-06_jprb, 1.28159e-06_jprb, 1.26373e-06_jprb, &
     & 1.24612e-06_jprb, 1.22875e-06_jprb, 1.21163e-06_jprb, 1.19475e-06_jprb/)
       kao_mn2( 3, :, 7) = (/ &
     & 1.53883e-06_jprb, 1.51733e-06_jprb, 1.49613e-06_jprb, 1.47522e-06_jprb, 1.45461e-06_jprb, &
     & 1.43428e-06_jprb, 1.41424e-06_jprb, 1.39448e-06_jprb, 1.37499e-06_jprb, 1.35578e-06_jprb, &
     & 1.33684e-06_jprb, 1.31816e-06_jprb, 1.29974e-06_jprb, 1.28157e-06_jprb, 1.26367e-06_jprb, &
     & 1.24601e-06_jprb, 1.22860e-06_jprb, 1.21143e-06_jprb, 1.19450e-06_jprb/)
       kao_mn2( 4, :, 7) = (/ &
     & 1.53789e-06_jprb, 1.51645e-06_jprb, 1.49532e-06_jprb, 1.47448e-06_jprb, 1.45393e-06_jprb, &
     & 1.43366e-06_jprb, 1.41368e-06_jprb, 1.39398e-06_jprb, 1.37455e-06_jprb, 1.35539e-06_jprb, &
     & 1.33650e-06_jprb, 1.31787e-06_jprb, 1.29950e-06_jprb, 1.28139e-06_jprb, 1.26353e-06_jprb, &
     & 1.24592e-06_jprb, 1.22856e-06_jprb, 1.21143e-06_jprb, 1.19455e-06_jprb/)
       kao_mn2( 5, :, 7) = (/ &
     & 1.54059e-06_jprb, 1.51888e-06_jprb, 1.49747e-06_jprb, 1.47637e-06_jprb, 1.45557e-06_jprb, &
     & 1.43505e-06_jprb, 1.41483e-06_jprb, 1.39489e-06_jprb, 1.37523e-06_jprb, 1.35585e-06_jprb, &
     & 1.33675e-06_jprb, 1.31791e-06_jprb, 1.29934e-06_jprb, 1.28103e-06_jprb, 1.26297e-06_jprb, &
     & 1.24517e-06_jprb, 1.22763e-06_jprb, 1.21033e-06_jprb, 1.19327e-06_jprb/)
       kao_mn2( 6, :, 7) = (/ &
     & 1.70605e-06_jprb, 1.65759e-06_jprb, 1.61052e-06_jprb, 1.56478e-06_jprb, 1.52034e-06_jprb, &
     & 1.47716e-06_jprb, 1.43521e-06_jprb, 1.39445e-06_jprb, 1.35485e-06_jprb, 1.31637e-06_jprb, &
     & 1.27898e-06_jprb, 1.24266e-06_jprb, 1.20737e-06_jprb, 1.17308e-06_jprb, 1.13976e-06_jprb, &
     & 1.10739e-06_jprb, 1.07594e-06_jprb, 1.04539e-06_jprb, 1.01570e-06_jprb/)
       kao_mn2( 7, :, 7) = (/ &
     & 1.39128e-06_jprb, 1.36388e-06_jprb, 1.33702e-06_jprb, 1.31068e-06_jprb, 1.28487e-06_jprb, &
     & 1.25956e-06_jprb, 1.23475e-06_jprb, 1.21044e-06_jprb, 1.18659e-06_jprb, 1.16322e-06_jprb, &
     & 1.14031e-06_jprb, 1.11785e-06_jprb, 1.09584e-06_jprb, 1.07425e-06_jprb, 1.05309e-06_jprb, &
     & 1.03235e-06_jprb, 1.01202e-06_jprb, 9.92088e-07_jprb, 9.72548e-07_jprb/)
       kao_mn2( 8, :, 7) = (/ &
     & 1.15676e-06_jprb, 1.13709e-06_jprb, 1.11775e-06_jprb, 1.09874e-06_jprb, 1.08005e-06_jprb, &
     & 1.06168e-06_jprb, 1.04362e-06_jprb, 1.02587e-06_jprb, 1.00842e-06_jprb, 9.91271e-07_jprb, &
     & 9.74411e-07_jprb, 9.57838e-07_jprb, 9.41547e-07_jprb, 9.25532e-07_jprb, 9.09791e-07_jprb, &
     & 8.94316e-07_jprb, 8.79105e-07_jprb, 8.64153e-07_jprb, 8.49455e-07_jprb/)
       kao_mn2( 9, :, 7) = (/ &
     & 1.53483e-06_jprb, 1.51352e-06_jprb, 1.49252e-06_jprb, 1.47180e-06_jprb, 1.45138e-06_jprb, &
     & 1.43123e-06_jprb, 1.41137e-06_jprb, 1.39178e-06_jprb, 1.37246e-06_jprb, 1.35341e-06_jprb, &
     & 1.33463e-06_jprb, 1.31610e-06_jprb, 1.29784e-06_jprb, 1.27982e-06_jprb, 1.26206e-06_jprb, &
     & 1.24454e-06_jprb, 1.22727e-06_jprb, 1.21024e-06_jprb, 1.19344e-06_jprb/)
       kao_mn2( 1, :, 8) = (/ &
     & 1.70380e-06_jprb, 1.67470e-06_jprb, 1.64609e-06_jprb, 1.61796e-06_jprb, 1.59032e-06_jprb, &
     & 1.56315e-06_jprb, 1.53645e-06_jprb, 1.51020e-06_jprb, 1.48440e-06_jprb, 1.45904e-06_jprb, &
     & 1.43411e-06_jprb, 1.40961e-06_jprb, 1.38553e-06_jprb, 1.36186e-06_jprb, 1.33859e-06_jprb, &
     & 1.31572e-06_jprb, 1.29324e-06_jprb, 1.27115e-06_jprb, 1.24943e-06_jprb/)
       kao_mn2( 2, :, 8) = (/ &
     & 1.70380e-06_jprb, 1.67470e-06_jprb, 1.64609e-06_jprb, 1.61796e-06_jprb, 1.59032e-06_jprb, &
     & 1.56315e-06_jprb, 1.53645e-06_jprb, 1.51020e-06_jprb, 1.48440e-06_jprb, 1.45904e-06_jprb, &
     & 1.43411e-06_jprb, 1.40961e-06_jprb, 1.38553e-06_jprb, 1.36186e-06_jprb, 1.33859e-06_jprb, &
     & 1.31572e-06_jprb, 1.29324e-06_jprb, 1.27115e-06_jprb, 1.24943e-06_jprb/)
       kao_mn2( 3, :, 8) = (/ &
     & 1.70380e-06_jprb, 1.67470e-06_jprb, 1.64609e-06_jprb, 1.61796e-06_jprb, 1.59032e-06_jprb, &
     & 1.56315e-06_jprb, 1.53645e-06_jprb, 1.51020e-06_jprb, 1.48440e-06_jprb, 1.45904e-06_jprb, &
     & 1.43411e-06_jprb, 1.40961e-06_jprb, 1.38553e-06_jprb, 1.36186e-06_jprb, 1.33859e-06_jprb, &
     & 1.31572e-06_jprb, 1.29324e-06_jprb, 1.27115e-06_jprb, 1.24943e-06_jprb/)
       kao_mn2( 4, :, 8) = (/ &
     & 1.70380e-06_jprb, 1.67470e-06_jprb, 1.64609e-06_jprb, 1.61796e-06_jprb, 1.59032e-06_jprb, &
     & 1.56315e-06_jprb, 1.53645e-06_jprb, 1.51020e-06_jprb, 1.48440e-06_jprb, 1.45904e-06_jprb, &
     & 1.43411e-06_jprb, 1.40961e-06_jprb, 1.38553e-06_jprb, 1.36186e-06_jprb, 1.33859e-06_jprb, &
     & 1.31572e-06_jprb, 1.29324e-06_jprb, 1.27115e-06_jprb, 1.24943e-06_jprb/)
       kao_mn2( 5, :, 8) = (/ &
     & 1.70380e-06_jprb, 1.67470e-06_jprb, 1.64609e-06_jprb, 1.61796e-06_jprb, 1.59032e-06_jprb, &
     & 1.56315e-06_jprb, 1.53645e-06_jprb, 1.51020e-06_jprb, 1.48440e-06_jprb, 1.45904e-06_jprb, &
     & 1.43411e-06_jprb, 1.40961e-06_jprb, 1.38553e-06_jprb, 1.36186e-06_jprb, 1.33859e-06_jprb, &
     & 1.31572e-06_jprb, 1.29324e-06_jprb, 1.27115e-06_jprb, 1.24943e-06_jprb/)
       kao_mn2( 6, :, 8) = (/ &
     & 1.70380e-06_jprb, 1.67470e-06_jprb, 1.64609e-06_jprb, 1.61796e-06_jprb, 1.59032e-06_jprb, &
     & 1.56315e-06_jprb, 1.53645e-06_jprb, 1.51020e-06_jprb, 1.48440e-06_jprb, 1.45904e-06_jprb, &
     & 1.43411e-06_jprb, 1.40961e-06_jprb, 1.38553e-06_jprb, 1.36186e-06_jprb, 1.33859e-06_jprb, &
     & 1.31572e-06_jprb, 1.29324e-06_jprb, 1.27115e-06_jprb, 1.24943e-06_jprb/)
       kao_mn2( 7, :, 8) = (/ &
     & 1.71827e-06_jprb, 1.65481e-06_jprb, 1.59370e-06_jprb, 1.53484e-06_jprb, 1.47816e-06_jprb, &
     & 1.42357e-06_jprb, 1.37099e-06_jprb, 1.32036e-06_jprb, 1.27160e-06_jprb, 1.22464e-06_jprb, &
     & 1.17941e-06_jprb, 1.13585e-06_jprb, 1.09390e-06_jprb, 1.05350e-06_jprb, 1.01459e-06_jprb, &
     & 9.77124e-07_jprb, 9.41037e-07_jprb, 9.06284e-07_jprb, 8.72813e-07_jprb/)
       kao_mn2( 8, :, 8) = (/ &
     & 1.77169e-06_jprb, 1.62858e-06_jprb, 1.49703e-06_jprb, 1.37610e-06_jprb, 1.26494e-06_jprb, &
     & 1.16276e-06_jprb, 1.06883e-06_jprb, 9.82495e-07_jprb, 9.03131e-07_jprb, 8.30177e-07_jprb, &
     & 7.63117e-07_jprb, 7.01473e-07_jprb, 6.44810e-07_jprb, 5.92723e-07_jprb, 5.44844e-07_jprb, &
     & 5.00832e-07_jprb, 4.60376e-07_jprb, 4.23187e-07_jprb, 3.89003e-07_jprb/)
       kao_mn2( 9, :, 8) = (/ &
     & 1.70025e-06_jprb, 1.67042e-06_jprb, 1.64110e-06_jprb, 1.61231e-06_jprb, 1.58401e-06_jprb, &
     & 1.55622e-06_jprb, 1.52891e-06_jprb, 1.50208e-06_jprb, 1.47572e-06_jprb, 1.44982e-06_jprb, &
     & 1.42438e-06_jprb, 1.39939e-06_jprb, 1.37483e-06_jprb, 1.35071e-06_jprb, 1.32700e-06_jprb, &
     & 1.30372e-06_jprb, 1.28084e-06_jprb, 1.25836e-06_jprb, 1.23628e-06_jprb/)
       kao_mn2( 1, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 2, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 3, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 4, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 5, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 6, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 7, :, 9) = (/ &
     & 1.74004e-06_jprb, 1.70661e-06_jprb, 1.67383e-06_jprb, 1.64167e-06_jprb, 1.61014e-06_jprb, &
     & 1.57921e-06_jprb, 1.54887e-06_jprb, 1.51912e-06_jprb, 1.48994e-06_jprb, 1.46132e-06_jprb, &
     & 1.43325e-06_jprb, 1.40572e-06_jprb, 1.37871e-06_jprb, 1.35223e-06_jprb, 1.32625e-06_jprb, &
     & 1.30078e-06_jprb, 1.27579e-06_jprb, 1.25128e-06_jprb, 1.22725e-06_jprb/)
       kao_mn2( 8, :, 9) = (/ &
     & 1.08654e-06_jprb, 1.09039e-06_jprb, 1.09425e-06_jprb, 1.09812e-06_jprb, 1.10201e-06_jprb, &
     & 1.10592e-06_jprb, 1.10983e-06_jprb, 1.11376e-06_jprb, 1.11771e-06_jprb, 1.12167e-06_jprb, &
     & 1.12564e-06_jprb, 1.12962e-06_jprb, 1.13363e-06_jprb, 1.13764e-06_jprb, 1.14167e-06_jprb, &
     & 1.14571e-06_jprb, 1.14977e-06_jprb, 1.15384e-06_jprb, 1.15793e-06_jprb/)
       kao_mn2( 9, :, 9) = (/ &
     & 1.74382e-06_jprb, 1.71092e-06_jprb, 1.67864e-06_jprb, 1.64697e-06_jprb, 1.61589e-06_jprb, &
     & 1.58541e-06_jprb, 1.55549e-06_jprb, 1.52615e-06_jprb, 1.49735e-06_jprb, 1.46910e-06_jprb, &
     & 1.44138e-06_jprb, 1.41419e-06_jprb, 1.38751e-06_jprb, 1.36133e-06_jprb, 1.33564e-06_jprb, &
     & 1.31045e-06_jprb, 1.28572e-06_jprb, 1.26146e-06_jprb, 1.23766e-06_jprb/)
       kao_mn2( 1, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 2, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 3, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 4, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 5, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 6, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 7, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 8, :,10) = (/ &
     & 1.73703e-06_jprb, 1.70249e-06_jprb, 1.66863e-06_jprb, 1.63544e-06_jprb, 1.60292e-06_jprb, &
     & 1.57104e-06_jprb, 1.53980e-06_jprb, 1.50917e-06_jprb, 1.47916e-06_jprb, 1.44974e-06_jprb, &
     & 1.42091e-06_jprb, 1.39265e-06_jprb, 1.36496e-06_jprb, 1.33781e-06_jprb, 1.31121e-06_jprb, &
     & 1.28513e-06_jprb, 1.25957e-06_jprb, 1.23452e-06_jprb, 1.20997e-06_jprb/)
       kao_mn2( 9, :,10) = (/ &
     & 1.82903e-06_jprb, 1.78673e-06_jprb, 1.74541e-06_jprb, 1.70505e-06_jprb, 1.66562e-06_jprb, &
     & 1.62710e-06_jprb, 1.58947e-06_jprb, 1.55271e-06_jprb, 1.51680e-06_jprb, 1.48172e-06_jprb, &
     & 1.44745e-06_jprb, 1.41398e-06_jprb, 1.38128e-06_jprb, 1.34933e-06_jprb, 1.31813e-06_jprb, &
     & 1.28765e-06_jprb, 1.25787e-06_jprb, 1.22878e-06_jprb, 1.20036e-06_jprb/)
       kao_mn2( 1, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 2, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 3, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 4, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 5, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 6, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 7, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 8, :,11) = (/ &
     & 1.73118e-06_jprb, 1.69710e-06_jprb, 1.66370e-06_jprb, 1.63095e-06_jprb, 1.59885e-06_jprb, &
     & 1.56737e-06_jprb, 1.53652e-06_jprb, 1.50628e-06_jprb, 1.47663e-06_jprb, 1.44756e-06_jprb, &
     & 1.41907e-06_jprb, 1.39114e-06_jprb, 1.36376e-06_jprb, 1.33691e-06_jprb, 1.31060e-06_jprb, &
     & 1.28480e-06_jprb, 1.25951e-06_jprb, 1.23472e-06_jprb, 1.21041e-06_jprb/)
       kao_mn2( 9, :,11) = (/ &
     & 1.81037e-06_jprb, 1.76948e-06_jprb, 1.72952e-06_jprb, 1.69045e-06_jprb, 1.65228e-06_jprb, &
     & 1.61496e-06_jprb, 1.57848e-06_jprb, 1.54283e-06_jprb, 1.50799e-06_jprb, 1.47393e-06_jprb, &
     & 1.44064e-06_jprb, 1.40810e-06_jprb, 1.37630e-06_jprb, 1.34522e-06_jprb, 1.31484e-06_jprb, &
     & 1.28514e-06_jprb, 1.25611e-06_jprb, 1.22774e-06_jprb, 1.20002e-06_jprb/)
       kao_mn2( 1, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 2, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 3, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 4, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 5, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 6, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 7, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 8, :,12) = (/ &
     & 1.73338e-06_jprb, 1.69915e-06_jprb, 1.66560e-06_jprb, 1.63271e-06_jprb, 1.60046e-06_jprb, &
     & 1.56886e-06_jprb, 1.53788e-06_jprb, 1.50751e-06_jprb, 1.47774e-06_jprb, 1.44856e-06_jprb, &
     & 1.41995e-06_jprb, 1.39191e-06_jprb, 1.36442e-06_jprb, 1.33748e-06_jprb, 1.31107e-06_jprb, &
     & 1.28518e-06_jprb, 1.25980e-06_jprb, 1.23492e-06_jprb, 1.21053e-06_jprb/)
       kao_mn2( 9, :,12) = (/ &
     & 2.04857e-06_jprb, 1.98353e-06_jprb, 1.92055e-06_jprb, 1.85957e-06_jprb, 1.80053e-06_jprb, &
     & 1.74336e-06_jprb, 1.68800e-06_jprb, 1.63441e-06_jprb, 1.58251e-06_jprb, 1.53227e-06_jprb, &
     & 1.48362e-06_jprb, 1.43651e-06_jprb, 1.39090e-06_jprb, 1.34674e-06_jprb, 1.30398e-06_jprb, &
     & 1.26257e-06_jprb, 1.22249e-06_jprb, 1.18367e-06_jprb, 1.14609e-06_jprb/)
       kao_mn2( 1, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 2, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 3, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 4, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 5, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 6, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 7, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 8, :,13) = (/ &
     & 1.73511e-06_jprb, 1.70072e-06_jprb, 1.66702e-06_jprb, 1.63398e-06_jprb, 1.60159e-06_jprb, &
     & 1.56985e-06_jprb, 1.53874e-06_jprb, 1.50824e-06_jprb, 1.47835e-06_jprb, 1.44905e-06_jprb, &
     & 1.42033e-06_jprb, 1.39218e-06_jprb, 1.36459e-06_jprb, 1.33755e-06_jprb, 1.31104e-06_jprb, &
     & 1.28505e-06_jprb, 1.25958e-06_jprb, 1.23462e-06_jprb, 1.21015e-06_jprb/)
       kao_mn2( 9, :,13) = (/ &
     & 2.13403e-06_jprb, 2.05906e-06_jprb, 1.98673e-06_jprb, 1.91694e-06_jprb, 1.84961e-06_jprb, &
     & 1.78463e-06_jprb, 1.72194e-06_jprb, 1.66145e-06_jprb, 1.60309e-06_jprb, 1.54678e-06_jprb, &
     & 1.49244e-06_jprb, 1.44002e-06_jprb, 1.38943e-06_jprb, 1.34062e-06_jprb, 1.29353e-06_jprb, &
     & 1.24809e-06_jprb, 1.20425e-06_jprb, 1.16195e-06_jprb, 1.12113e-06_jprb/)
       kao_mn2( 1, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 2, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 3, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 4, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 5, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 6, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 7, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 8, :,14) = (/ &
     & 1.73398e-06_jprb, 1.69941e-06_jprb, 1.66553e-06_jprb, 1.63233e-06_jprb, 1.59979e-06_jprb, &
     & 1.56790e-06_jprb, 1.53664e-06_jprb, 1.50601e-06_jprb, 1.47598e-06_jprb, 1.44656e-06_jprb, &
     & 1.41772e-06_jprb, 1.38946e-06_jprb, 1.36176e-06_jprb, 1.33461e-06_jprb, 1.30801e-06_jprb, &
     & 1.28193e-06_jprb, 1.25637e-06_jprb, 1.23133e-06_jprb, 1.20678e-06_jprb/)
       kao_mn2( 9, :,14) = (/ &
     & 1.83423e-06_jprb, 1.79123e-06_jprb, 1.74923e-06_jprb, 1.70821e-06_jprb, 1.66816e-06_jprb, &
     & 1.62904e-06_jprb, 1.59085e-06_jprb, 1.55354e-06_jprb, 1.51712e-06_jprb, 1.48154e-06_jprb, &
     & 1.44681e-06_jprb, 1.41288e-06_jprb, 1.37975e-06_jprb, 1.34740e-06_jprb, 1.31581e-06_jprb, &
     & 1.28496e-06_jprb, 1.25483e-06_jprb, 1.22540e-06_jprb, 1.19667e-06_jprb/)
       kao_mn2( 1, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 2, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 3, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 4, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 5, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 6, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 7, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 8, :,15) = (/ &
     & 1.73231e-06_jprb, 1.69765e-06_jprb, 1.66368e-06_jprb, 1.63039e-06_jprb, 1.59776e-06_jprb, &
     & 1.56579e-06_jprb, 1.53445e-06_jprb, 1.50375e-06_jprb, 1.47366e-06_jprb, 1.44417e-06_jprb, &
     & 1.41527e-06_jprb, 1.38695e-06_jprb, 1.35919e-06_jprb, 1.33199e-06_jprb, 1.30534e-06_jprb, &
     & 1.27922e-06_jprb, 1.25362e-06_jprb, 1.22853e-06_jprb, 1.20395e-06_jprb/)
       kao_mn2( 9, :,15) = (/ &
     & 1.71602e-06_jprb, 1.68499e-06_jprb, 1.65452e-06_jprb, 1.62461e-06_jprb, 1.59523e-06_jprb, &
     & 1.56639e-06_jprb, 1.53807e-06_jprb, 1.51026e-06_jprb, 1.48295e-06_jprb, 1.45614e-06_jprb, &
     & 1.42981e-06_jprb, 1.40395e-06_jprb, 1.37857e-06_jprb, 1.35364e-06_jprb, 1.32917e-06_jprb, &
     & 1.30513e-06_jprb, 1.28153e-06_jprb, 1.25836e-06_jprb, 1.23561e-06_jprb/)
       kao_mn2( 1, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 2, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 3, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 4, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 5, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 6, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 7, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 8, :,16) = (/ &
     & 1.73310e-06_jprb, 1.69826e-06_jprb, 1.66413e-06_jprb, 1.63069e-06_jprb, 1.59791e-06_jprb, &
     & 1.56580e-06_jprb, 1.53433e-06_jprb, 1.50349e-06_jprb, 1.47328e-06_jprb, 1.44367e-06_jprb, &
     & 1.41465e-06_jprb, 1.38622e-06_jprb, 1.35836e-06_jprb, 1.33106e-06_jprb, 1.30431e-06_jprb, &
     & 1.27810e-06_jprb, 1.25241e-06_jprb, 1.22724e-06_jprb, 1.20257e-06_jprb/)
       kao_mn2( 9, :,16) = (/ &
     & 1.79375e-06_jprb, 1.75599e-06_jprb, 1.71903e-06_jprb, 1.68284e-06_jprb, 1.64741e-06_jprb, &
     & 1.61273e-06_jprb, 1.57878e-06_jprb, 1.54554e-06_jprb, 1.51301e-06_jprb, 1.48116e-06_jprb, &
     & 1.44998e-06_jprb, 1.41945e-06_jprb, 1.38957e-06_jprb, 1.36032e-06_jprb, 1.33168e-06_jprb, &
     & 1.30365e-06_jprb, 1.27620e-06_jprb, 1.24934e-06_jprb, 1.22304e-06_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.1755e-06_jprb,6.5398e-07_jprb,4.3915e-07_jprb,3.0753e-07_jprb,1.9677e-07_jprb,1.4362e-07_jprb, &
     &9.4598e-08_jprb,1.1848e-07_jprb,1.4280e-07_jprb,1.5821e-07_jprb,1.5816e-07_jprb,1.5769e-07_jprb, &
     &1.5844e-07_jprb,1.6016e-07_jprb,1.6232e-07_jprb,1.6320e-07_jprb/)
      forrefo(2,:) = (/ &
     &1.0703e-06_jprb,6.2783e-07_jprb,4.7122e-07_jprb,2.6300e-07_jprb,1.8538e-07_jprb,1.5076e-07_jprb, &
     &1.9474e-07_jprb,2.9543e-07_jprb,2.0093e-07_jprb,1.5819e-07_jprb,1.5826e-07_jprb,1.5737e-07_jprb, &
     &1.5751e-07_jprb,1.5910e-07_jprb,1.6181e-07_jprb,1.6320e-07_jprb/)
      forrefo(3,:) = (/ &
     &1.0470e-06_jprb,5.8184e-07_jprb,4.8218e-07_jprb,2.7771e-07_jprb,1.9036e-07_jprb,1.5737e-07_jprb, &
     &1.8633e-07_jprb,2.5754e-07_jprb,4.0647e-07_jprb,1.5839e-07_jprb,1.5914e-07_jprb,1.5788e-07_jprb, &
     &1.5731e-07_jprb,1.5836e-07_jprb,1.6103e-07_jprb,1.6320e-07_jprb/)
      forrefo(4,:) = (/ &
     &1.3891e-06_jprb,5.4901e-07_jprb,2.8850e-07_jprb,1.9176e-07_jprb,1.4549e-07_jprb,1.3603e-07_jprb, &
     &1.7472e-07_jprb,2.9796e-07_jprb,3.2452e-07_jprb,2.5231e-07_jprb,2.8195e-07_jprb,1.5527e-07_jprb, &
     &1.5507e-07_jprb,1.5442e-07_jprb,1.5275e-07_jprb,1.6057e-07_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

     selfrefo(:, 1) = (/ &
     & 1.73980e-03_jprb, 1.41928e-03_jprb, 1.15780e-03_jprb, 9.44496e-04_jprb, 7.70490e-04_jprb, &
     & 6.28541e-04_jprb, 5.12744e-04_jprb, 4.18280e-04_jprb, 3.41219e-04_jprb, 2.78356e-04_jprb/)
      selfrefo(:, 2) = (/ &
     & 1.84082e-03_jprb, 1.50228e-03_jprb, 1.22600e-03_jprb, 1.00053e-03_jprb, 8.16525e-04_jprb, &
     & 6.66359e-04_jprb, 5.43811e-04_jprb, 4.43800e-04_jprb, 3.62182e-04_jprb, 2.95574e-04_jprb/)
      selfrefo(:, 3) = (/ &
     & 1.92957e-03_jprb, 1.57727e-03_jprb, 1.28930e-03_jprb, 1.05390e-03_jprb, 8.61484e-04_jprb, &
     & 7.04197e-04_jprb, 5.75627e-04_jprb, 4.70530e-04_jprb, 3.84622e-04_jprb, 3.14399e-04_jprb/)
      selfrefo(:, 4) = (/ &
     & 2.12958e-03_jprb, 1.73572e-03_jprb, 1.41470e-03_jprb, 1.15305e-03_jprb, 9.39798e-04_jprb, &
     & 7.65984e-04_jprb, 6.24317e-04_jprb, 5.08850e-04_jprb, 4.14739e-04_jprb, 3.38034e-04_jprb/)
      selfrefo(:, 5) = (/ &
     & 2.30636e-03_jprb, 1.88401e-03_jprb, 1.53900e-03_jprb, 1.25717e-03_jprb, 1.02695e-03_jprb, &
     & 8.38891e-04_jprb, 6.85270e-04_jprb, 5.59780e-04_jprb, 4.57270e-04_jprb, 3.73533e-04_jprb/)
      selfrefo(:, 6) = (/ &
     & 2.47824e-03_jprb, 2.03278e-03_jprb, 1.66740e-03_jprb, 1.36769e-03_jprb, 1.12185e-03_jprb, &
     & 9.20206e-04_jprb, 7.54803e-04_jprb, 6.19130e-04_jprb, 5.07844e-04_jprb, 4.16561e-04_jprb/)
      selfrefo(:, 7) = (/ &
     & 2.54196e-03_jprb, 2.10768e-03_jprb, 1.74760e-03_jprb, 1.44904e-03_jprb, 1.20148e-03_jprb, &
     & 9.96215e-04_jprb, 8.26019e-04_jprb, 6.84900e-04_jprb, 5.67890e-04_jprb, 4.70870e-04_jprb/)
      selfrefo(:, 8) = (/ &
     & 2.52650e-03_jprb, 2.11773e-03_jprb, 1.77510e-03_jprb, 1.48790e-03_jprb, 1.24717e-03_jprb, &
     & 1.04539e-03_jprb, 8.76251e-04_jprb, 7.34480e-04_jprb, 6.15646e-04_jprb, 5.16039e-04_jprb/)
      selfrefo(:, 9) = (/ &
     & 2.82351e-03_jprb, 2.34652e-03_jprb, 1.95010e-03_jprb, 1.62065e-03_jprb, 1.34686e-03_jprb, &
     & 1.11933e-03_jprb, 9.30232e-04_jprb, 7.73080e-04_jprb, 6.42477e-04_jprb, 5.33939e-04_jprb/)
      selfrefo(:,10) = (/ &
     & 2.98189e-03_jprb, 2.46741e-03_jprb, 2.04170e-03_jprb, 1.68944e-03_jprb, 1.39795e-03_jprb, &
     & 1.15676e-03_jprb, 9.57176e-04_jprb, 7.92030e-04_jprb, 6.55377e-04_jprb, 5.42302e-04_jprb/)
      selfrefo(:,11) = (/ &
     & 2.98239e-03_jprb, 2.46774e-03_jprb, 2.04190e-03_jprb, 1.68954e-03_jprb, 1.39799e-03_jprb, &
     & 1.15675e-03_jprb, 9.57137e-04_jprb, 7.91970e-04_jprb, 6.55305e-04_jprb, 5.42224e-04_jprb/)
      selfrefo(:,12) = (/ &
     & 2.97833e-03_jprb, 2.46461e-03_jprb, 2.03950e-03_jprb, 1.68772e-03_jprb, 1.39661e-03_jprb, &
     & 1.15571e-03_jprb, 9.56370e-04_jprb, 7.91410e-04_jprb, 6.54903e-04_jprb, 5.41942e-04_jprb/)
      selfrefo(:,13) = (/ &
     & 2.97779e-03_jprb, 2.46463e-03_jprb, 2.03990e-03_jprb, 1.68836e-03_jprb, 1.39741e-03_jprb, &
     & 1.15659e-03_jprb, 9.57278e-04_jprb, 7.92310e-04_jprb, 6.55771e-04_jprb, 5.42762e-04_jprb/)
      selfrefo(:,14) = (/ &
     & 2.98326e-03_jprb, 2.46943e-03_jprb, 2.04410e-03_jprb, 1.69203e-03_jprb, 1.40060e-03_jprb, &
     & 1.15936e-03_jprb, 9.59673e-04_jprb, 7.94380e-04_jprb, 6.57557e-04_jprb, 5.44301e-04_jprb/)
      selfrefo(:,15) = (/ &
     & 2.99407e-03_jprb, 2.47825e-03_jprb, 2.05130e-03_jprb, 1.69790e-03_jprb, 1.40539e-03_jprb, &
     & 1.16327e-03_jprb, 9.62862e-04_jprb, 7.96980e-04_jprb, 6.59676e-04_jprb, 5.46028e-04_jprb/)
      selfrefo(:,16) = (/ &
     & 3.00005e-03_jprb, 2.48296e-03_jprb, 2.05500e-03_jprb, 1.70080e-03_jprb, 1.40765e-03_jprb, &
     & 1.16503e-03_jprb, 9.64224e-04_jprb, 7.98030e-04_jprb, 6.60481e-04_jprb, 5.46641e-04_jprb/)

if (lhook) call dr_hook('rrtm_kgb15',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb15:error reading file radrrtm")

end subroutine rrtm_kgb15
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

