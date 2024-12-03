! # 1 "ifsrrtm/rrtm_kgb13.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb13.f90"
subroutine rrtm_kgb13

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
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

use yoerrto13, only : kao, kao_d, selfrefo, forrefo ,fracrefao, fracrefbo, kao_mco2, kao_mco,kbo_mo3


!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 29 "ifsrrtm/rrtm_kgb13.f90" 2

if (lhook) call dr_hook('rrtm_kgb13',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d
  kao = real(kao_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb13:')
endif

! planck fraction mapping level : p=473.4280 mb, t = 259.83 k      
      fracrefao(:, 1) = (/ &
     &  1.7534e-01_jprb,1.7394e-01_jprb,1.6089e-01_jprb,1.3782e-01_jprb,1.0696e-01_jprb,8.5853e-02_jprb, &
     &  6.6548e-02_jprb,4.9053e-02_jprb,3.2064e-02_jprb,3.4820e-03_jprb,2.8763e-03_jprb,2.2204e-03_jprb, &
     &  1.5612e-03_jprb,9.8572e-04_jprb,3.6853e-04_jprb,5.1612e-05_jprb/)
      fracrefao(:, 2) = (/ &
     &  1.7489e-01_jprb,1.7309e-01_jprb,1.5981e-01_jprb,1.3782e-01_jprb,1.0797e-01_jprb,8.6367e-02_jprb, &
     &  6.7042e-02_jprb,4.9257e-02_jprb,3.2207e-02_jprb,3.4820e-03_jprb,2.8767e-03_jprb,2.2203e-03_jprb, &
     &  1.5613e-03_jprb,9.8571e-04_jprb,3.6853e-04_jprb,5.1612e-05_jprb/)
      fracrefao(:, 3) = (/ &
     &  1.7459e-01_jprb,1.7259e-01_jprb,1.5948e-01_jprb,1.3694e-01_jprb,1.0815e-01_jprb,8.7376e-02_jprb, &
     &  6.7339e-02_jprb,4.9541e-02_jprb,3.2333e-02_jprb,3.5019e-03_jprb,2.8958e-03_jprb,2.2527e-03_jprb, &
     &  1.6099e-03_jprb,9.8574e-04_jprb,3.6853e-04_jprb,5.1612e-05_jprb/)
      fracrefao(:, 4) = (/ &
     &  1.7391e-01_jprb,1.7244e-01_jprb,1.5921e-01_jprb,1.3644e-01_jprb,1.0787e-01_jprb,8.7776e-02_jprb, &
     &  6.8361e-02_jprb,4.9628e-02_jprb,3.2578e-02_jprb,3.5117e-03_jprb,2.9064e-03_jprb,2.2571e-03_jprb, &
     &  1.6887e-03_jprb,1.0045e-03_jprb,3.6853e-04_jprb,5.1612e-05_jprb/)
      fracrefao(:, 5) = (/ &
     &  1.7338e-01_jprb,1.7157e-01_jprb,1.5957e-01_jprb,1.3571e-01_jprb,1.0773e-01_jprb,8.7966e-02_jprb, &
     &  6.9000e-02_jprb,5.0300e-02_jprb,3.2813e-02_jprb,3.5470e-03_jprb,2.9425e-03_jprb,2.2552e-03_jprb, &
     &  1.7038e-03_jprb,1.1025e-03_jprb,3.6853e-04_jprb,5.1612e-05_jprb/)
      fracrefao(:, 6) = (/ &
     &  1.7230e-01_jprb,1.7082e-01_jprb,1.5917e-01_jprb,1.3562e-01_jprb,1.0806e-01_jprb,8.7635e-02_jprb, &
     &  6.9815e-02_jprb,5.1155e-02_jprb,3.3139e-02_jprb,3.6264e-03_jprb,2.9436e-03_jprb,2.3417e-03_jprb, &
     &  1.7731e-03_jprb,1.1156e-03_jprb,4.4533e-04_jprb,5.1612e-05_jprb/)
      fracrefao(:, 7) = (/ &
     &  1.7073e-01_jprb,1.6961e-01_jprb,1.5844e-01_jprb,1.3594e-01_jprb,1.0821e-01_jprb,8.7791e-02_jprb, &
     &  7.0502e-02_jprb,5.1904e-02_jprb,3.4107e-02_jprb,3.5888e-03_jprb,2.9574e-03_jprb,2.5851e-03_jprb, &
     &  1.9127e-03_jprb,1.1537e-03_jprb,4.7789e-04_jprb,1.0016e-04_jprb/)
      fracrefao(:, 8) = (/ &
     &  1.6700e-01_jprb,1.6848e-01_jprb,1.5628e-01_jprb,1.3448e-01_jprb,1.1011e-01_jprb,8.9016e-02_jprb, &
     &  7.1973e-02_jprb,5.2798e-02_jprb,3.5650e-02_jprb,3.8534e-03_jprb,3.4142e-03_jprb,2.7799e-03_jprb, &
     &  2.1288e-03_jprb,1.3043e-03_jprb,6.2858e-04_jprb,1.0016e-04_jprb/)
      fracrefao(:, 9) = (/ &
     &  1.6338e-01_jprb,1.5565e-01_jprb,1.4470e-01_jprb,1.3500e-01_jprb,1.1909e-01_jprb,9.8312e-02_jprb, &
     &  7.9023e-02_jprb,5.5728e-02_jprb,3.6831e-02_jprb,3.6569e-03_jprb,3.0552e-03_jprb,2.3431e-03_jprb, &
     &  1.7088e-03_jprb,1.1082e-03_jprb,3.6829e-04_jprb,5.1612e-05_jprb/)

! planck fraction mapping level : p=4.758820 mb, t = 250.85 k
      fracrefbo(:) = (/ &
     &  1.5411e-01_jprb,1.3573e-01_jprb,1.2527e-01_jprb,1.2698e-01_jprb,1.2394e-01_jprb,1.0876e-01_jprb, &
     &  8.9906e-02_jprb,6.9551e-02_jprb,4.8240e-02_jprb,5.2434e-03_jprb,4.3630e-03_jprb,3.4262e-03_jprb, &
     &  2.5124e-03_jprb,1.5479e-03_jprb,3.7294e-04_jprb,5.1050e-05_jprb/)


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

      kao_mco2( 1, :, 1) = (/ &
     & 1.09539e-04_jprb, 1.17067e-04_jprb, 1.25113e-04_jprb, 1.33712e-04_jprb, 1.42902e-04_jprb, &
     & 1.52724e-04_jprb, 1.63221e-04_jprb, 1.74439e-04_jprb, 1.86428e-04_jprb, 1.99241e-04_jprb, &
     & 2.12934e-04_jprb, 2.27569e-04_jprb, 2.43210e-04_jprb, 2.59926e-04_jprb, 2.77790e-04_jprb, &
     & 2.96883e-04_jprb, 3.17287e-04_jprb, 3.39094e-04_jprb, 3.62400e-04_jprb/)
      kao_mco2( 2, :, 1) = (/ &
     & 1.25202e-04_jprb, 1.34718e-04_jprb, 1.44957e-04_jprb, 1.55974e-04_jprb, 1.67829e-04_jprb, &
     & 1.80585e-04_jprb, 1.94311e-04_jprb, 2.09079e-04_jprb, 2.24971e-04_jprb, 2.42069e-04_jprb, &
     & 2.60468e-04_jprb, 2.80265e-04_jprb, 3.01567e-04_jprb, 3.24488e-04_jprb, 3.49150e-04_jprb, &
     & 3.75688e-04_jprb, 4.04242e-04_jprb, 4.34966e-04_jprb, 4.68026e-04_jprb/)
      kao_mco2( 3, :, 1) = (/ &
     & 1.12112e-04_jprb, 1.21090e-04_jprb, 1.30786e-04_jprb, 1.41259e-04_jprb, 1.52571e-04_jprb, &
     & 1.64788e-04_jprb, 1.77984e-04_jprb, 1.92237e-04_jprb, 2.07631e-04_jprb, 2.24257e-04_jprb, &
     & 2.42215e-04_jprb, 2.61611e-04_jprb, 2.82560e-04_jprb, 3.05187e-04_jprb, 3.29625e-04_jprb, &
     & 3.56021e-04_jprb, 3.84530e-04_jprb, 4.15322e-04_jprb, 4.48580e-04_jprb/)
      kao_mco2( 4, :, 1) = (/ &
     & 9.74130e-05_jprb, 1.05372e-04_jprb, 1.13982e-04_jprb, 1.23295e-04_jprb, 1.33369e-04_jprb, &
     & 1.44265e-04_jprb, 1.56053e-04_jprb, 1.68803e-04_jprb, 1.82595e-04_jprb, 1.97514e-04_jprb, &
     & 2.13652e-04_jprb, 2.31109e-04_jprb, 2.49992e-04_jprb, 2.70418e-04_jprb, 2.92512e-04_jprb, &
     & 3.16412e-04_jprb, 3.42265e-04_jprb, 3.70230e-04_jprb, 4.00479e-04_jprb/)
      kao_mco2( 5, :, 1) = (/ &
     & 8.71018e-05_jprb, 9.40759e-05_jprb, 1.01608e-04_jprb, 1.09744e-04_jprb, 1.18531e-04_jprb, &
     & 1.28022e-04_jprb, 1.38272e-04_jprb, 1.49343e-04_jprb, 1.61301e-04_jprb, 1.74216e-04_jprb, &
     & 1.88166e-04_jprb, 2.03232e-04_jprb, 2.19504e-04_jprb, 2.37079e-04_jprb, 2.56062e-04_jprb, &
     & 2.76565e-04_jprb, 2.98709e-04_jprb, 3.22626e-04_jprb, 3.48458e-04_jprb/)
      kao_mco2( 6, :, 1) = (/ &
     & 7.55256e-05_jprb, 8.17414e-05_jprb, 8.84688e-05_jprb, 9.57500e-05_jprb, 1.03630e-04_jprb, &
     & 1.12159e-04_jprb, 1.21390e-04_jprb, 1.31381e-04_jprb, 1.42193e-04_jprb, 1.53896e-04_jprb, &
     & 1.66562e-04_jprb, 1.80270e-04_jprb, 1.95107e-04_jprb, 2.11164e-04_jprb, 2.28543e-04_jprb, &
     & 2.47353e-04_jprb, 2.67710e-04_jprb, 2.89743e-04_jprb, 3.13589e-04_jprb/)
      kao_mco2( 7, :, 1) = (/ &
     & 5.31515e-05_jprb, 6.06869e-05_jprb, 6.92907e-05_jprb, 7.91143e-05_jprb, 9.03306e-05_jprb, &
     & 1.03137e-04_jprb, 1.17759e-04_jprb, 1.34454e-04_jprb, 1.53516e-04_jprb, 1.75281e-04_jprb, &
     & 2.00131e-04_jprb, 2.28504e-04_jprb, 2.60900e-04_jprb, 2.97888e-04_jprb, 3.40121e-04_jprb, &
     & 3.88341e-04_jprb, 4.43397e-04_jprb, 5.06259e-04_jprb, 5.78033e-04_jprb/)
      kao_mco2( 8, :, 1) = (/ &
     & 2.52471e-04_jprb, 2.96005e-04_jprb, 3.47045e-04_jprb, 4.06886e-04_jprb, 4.77045e-04_jprb, &
     & 5.59302e-04_jprb, 6.55742e-04_jprb, 7.68811e-04_jprb, 9.01377e-04_jprb, 1.05680e-03_jprb, &
     & 1.23902e-03_jprb, 1.45267e-03_jprb, 1.70315e-03_jprb, 1.99683e-03_jprb, 2.34114e-03_jprb, &
     & 2.74482e-03_jprb, 3.21811e-03_jprb, 3.77300e-03_jprb, 4.42358e-03_jprb/)
      kao_mco2( 9, :, 1) = (/ &
     & 4.06711e-05_jprb, 4.53161e-05_jprb, 5.04917e-05_jprb, 5.62583e-05_jprb, 6.26836e-05_jprb, &
     & 6.98427e-05_jprb, 7.78194e-05_jprb, 8.67071e-05_jprb, 9.66100e-05_jprb, 1.07644e-04_jprb, &
     & 1.19938e-04_jprb, 1.33636e-04_jprb, 1.48898e-04_jprb, 1.65904e-04_jprb, 1.84852e-04_jprb, &
     & 2.05964e-04_jprb, 2.29487e-04_jprb, 2.55697e-04_jprb, 2.84900e-04_jprb/)
      kao_mco2( 1, :, 2) = (/ &
     & 2.01759e-04_jprb, 2.15641e-04_jprb, 2.30478e-04_jprb, 2.46336e-04_jprb, 2.63285e-04_jprb, &
     & 2.81400e-04_jprb, 3.00761e-04_jprb, 3.21455e-04_jprb, 3.43573e-04_jprb, 3.67212e-04_jprb, &
     & 3.92477e-04_jprb, 4.19482e-04_jprb, 4.48344e-04_jprb, 4.79192e-04_jprb, 5.12162e-04_jprb, &
     & 5.47401e-04_jprb, 5.85064e-04_jprb, 6.25319e-04_jprb, 6.68343e-04_jprb/)
      kao_mco2( 2, :, 2) = (/ &
     & 2.53461e-04_jprb, 2.70916e-04_jprb, 2.89574e-04_jprb, 3.09516e-04_jprb, 3.30832e-04_jprb, &
     & 3.53616e-04_jprb, 3.77969e-04_jprb, 4.03999e-04_jprb, 4.31822e-04_jprb, 4.61561e-04_jprb, &
     & 4.93348e-04_jprb, 5.27324e-04_jprb, 5.63640e-04_jprb, 6.02457e-04_jprb, 6.43948e-04_jprb, &
     & 6.88295e-04_jprb, 7.35697e-04_jprb, 7.86364e-04_jprb, 8.40519e-04_jprb/)
      kao_mco2( 3, :, 2) = (/ &
     & 2.58821e-04_jprb, 2.76943e-04_jprb, 2.96334e-04_jprb, 3.17082e-04_jprb, 3.39283e-04_jprb, &
     & 3.63038e-04_jprb, 3.88457e-04_jprb, 4.15655e-04_jprb, 4.44758e-04_jprb, 4.75899e-04_jprb, &
     & 5.09220e-04_jprb, 5.44874e-04_jprb, 5.83024e-04_jprb, 6.23845e-04_jprb, 6.67525e-04_jprb, &
     & 7.14263e-04_jprb, 7.64273e-04_jprb, 8.17785e-04_jprb, 8.75043e-04_jprb/)
      kao_mco2( 4, :, 2) = (/ &
     & 2.46588e-04_jprb, 2.64630e-04_jprb, 2.83993e-04_jprb, 3.04771e-04_jprb, 3.27071e-04_jprb, &
     & 3.51001e-04_jprb, 3.76683e-04_jprb, 4.04244e-04_jprb, 4.33821e-04_jprb, 4.65563e-04_jprb, &
     & 4.99627e-04_jprb, 5.36183e-04_jprb, 5.75414e-04_jprb, 6.17515e-04_jprb, 6.62697e-04_jprb, &
     & 7.11185e-04_jprb, 7.63220e-04_jprb, 8.19063e-04_jprb, 8.78991e-04_jprb/)
      kao_mco2( 5, :, 2) = (/ &
     & 2.19140e-04_jprb, 2.36464e-04_jprb, 2.55158e-04_jprb, 2.75330e-04_jprb, 2.97097e-04_jprb, &
     & 3.20585e-04_jprb, 3.45929e-04_jprb, 3.73277e-04_jprb, 4.02787e-04_jprb, 4.34630e-04_jprb, &
     & 4.68991e-04_jprb, 5.06068e-04_jprb, 5.46076e-04_jprb, 5.89247e-04_jprb, 6.35831e-04_jprb, &
     & 6.86097e-04_jprb, 7.40338e-04_jprb, 7.98867e-04_jprb, 8.62022e-04_jprb/)
      kao_mco2( 6, :, 2) = (/ &
     & 1.74073e-04_jprb, 1.92221e-04_jprb, 2.12260e-04_jprb, 2.34388e-04_jprb, 2.58824e-04_jprb, &
     & 2.85807e-04_jprb, 3.15603e-04_jprb, 3.48505e-04_jprb, 3.84837e-04_jprb, 4.24957e-04_jprb, &
     & 4.69260e-04_jprb, 5.18181e-04_jprb, 5.72202e-04_jprb, 6.31855e-04_jprb, 6.97727e-04_jprb, &
     & 7.70466e-04_jprb, 8.50789e-04_jprb, 9.39485e-04_jprb, 1.03743e-03_jprb/)
      kao_mco2( 7, :, 2) = (/ &
     & 1.74359e-04_jprb, 1.99276e-04_jprb, 2.27753e-04_jprb, 2.60299e-04_jprb, 2.97497e-04_jprb, &
     & 3.40010e-04_jprb, 3.88599e-04_jprb, 4.44130e-04_jprb, 5.07598e-04_jprb, 5.80135e-04_jprb, &
     & 6.63039e-04_jprb, 7.57789e-04_jprb, 8.66079e-04_jprb, 9.89845e-04_jprb, 1.13130e-03_jprb, &
     & 1.29296e-03_jprb, 1.47773e-03_jprb, 1.68890e-03_jprb, 1.93025e-03_jprb/)
      kao_mco2( 8, :, 2) = (/ &
     & 1.08215e-03_jprb, 1.20760e-03_jprb, 1.34759e-03_jprb, 1.50382e-03_jprb, 1.67815e-03_jprb, &
     & 1.87270e-03_jprb, 2.08980e-03_jprb, 2.33206e-03_jprb, 2.60242e-03_jprb, 2.90411e-03_jprb, &
     & 3.24078e-03_jprb, 3.61648e-03_jprb, 4.03573e-03_jprb, 4.50359e-03_jprb, 5.02568e-03_jprb, &
     & 5.60830e-03_jprb, 6.25846e-03_jprb, 6.98399e-03_jprb, 7.79363e-03_jprb/)
      kao_mco2( 9, :, 2) = (/ &
     & 1.04969e-04_jprb, 1.20766e-04_jprb, 1.38939e-04_jprb, 1.59848e-04_jprb, 1.83903e-04_jprb, &
     & 2.11578e-04_jprb, 2.43418e-04_jprb, 2.80049e-04_jprb, 3.22193e-04_jprb, 3.70678e-04_jprb, &
     & 4.26461e-04_jprb, 4.90638e-04_jprb, 5.64472e-04_jprb, 6.49418e-04_jprb, 7.47147e-04_jprb, &
     & 8.59583e-04_jprb, 9.88940e-04_jprb, 1.13776e-03_jprb, 1.30898e-03_jprb/)
      kao_mco2( 1, :, 3) = (/ &
     & 3.72106e-04_jprb, 3.96252e-04_jprb, 4.21966e-04_jprb, 4.49347e-04_jprb, 4.78506e-04_jprb, &
     & 5.09557e-04_jprb, 5.42623e-04_jprb, 5.77834e-04_jprb, 6.15330e-04_jprb, 6.55260e-04_jprb, &
     & 6.97781e-04_jprb, 7.43060e-04_jprb, 7.91278e-04_jprb, 8.42626e-04_jprb, 8.97304e-04_jprb, &
     & 9.55532e-04_jprb, 1.01754e-03_jprb, 1.08357e-03_jprb, 1.15388e-03_jprb/)
      kao_mco2( 2, :, 3) = (/ &
     & 4.20563e-04_jprb, 4.46162e-04_jprb, 4.73319e-04_jprb, 5.02130e-04_jprb, 5.32693e-04_jprb, &
     & 5.65118e-04_jprb, 5.99516e-04_jprb, 6.36007e-04_jprb, 6.74720e-04_jprb, 7.15789e-04_jprb, &
     & 7.59358e-04_jprb, 8.05579e-04_jprb, 8.54613e-04_jprb, 9.06632e-04_jprb, 9.61817e-04_jprb, &
     & 1.02036e-03_jprb, 1.08247e-03_jprb, 1.14836e-03_jprb, 1.21826e-03_jprb/)
      kao_mco2( 3, :, 3) = (/ &
     & 4.89664e-04_jprb, 5.18321e-04_jprb, 5.48654e-04_jprb, 5.80764e-04_jprb, 6.14752e-04_jprb, &
     & 6.50729e-04_jprb, 6.88812e-04_jprb, 7.29124e-04_jprb, 7.71795e-04_jprb, 8.16963e-04_jprb, &
     & 8.64774e-04_jprb, 9.15384e-04_jprb, 9.68955e-04_jprb, 1.02566e-03_jprb, 1.08569e-03_jprb, &
     & 1.14922e-03_jprb, 1.21648e-03_jprb, 1.28767e-03_jprb, 1.36303e-03_jprb/)
      kao_mco2( 4, :, 3) = (/ &
     & 4.61143e-04_jprb, 4.92198e-04_jprb, 5.25343e-04_jprb, 5.60720e-04_jprb, 5.98480e-04_jprb, &
     & 6.38783e-04_jprb, 6.81799e-04_jprb, 7.27713e-04_jprb, 7.76718e-04_jprb, 8.29023e-04_jprb, &
     & 8.84851e-04_jprb, 9.44438e-04_jprb, 1.00804e-03_jprb, 1.07592e-03_jprb, 1.14837e-03_jprb, &
     & 1.22571e-03_jprb, 1.30825e-03_jprb, 1.39635e-03_jprb, 1.49038e-03_jprb/)
      kao_mco2( 5, :, 3) = (/ &
     & 4.01988e-04_jprb, 4.36672e-04_jprb, 4.74349e-04_jprb, 5.15278e-04_jprb, 5.59737e-04_jprb, &
     & 6.08032e-04_jprb, 6.60495e-04_jprb, 7.17484e-04_jprb, 7.79390e-04_jprb, 8.46638e-04_jprb, &
     & 9.19688e-04_jprb, 9.99041e-04_jprb, 1.08524e-03_jprb, 1.17888e-03_jprb, 1.28059e-03_jprb, &
     & 1.39109e-03_jprb, 1.51111e-03_jprb, 1.64149e-03_jprb, 1.78313e-03_jprb/)
      kao_mco2( 6, :, 3) = (/ &
     & 3.35536e-04_jprb, 3.74371e-04_jprb, 4.17700e-04_jprb, 4.66045e-04_jprb, 5.19985e-04_jprb, &
     & 5.80169e-04_jprb, 6.47318e-04_jprb, 7.22238e-04_jprb, 8.05831e-04_jprb, 8.99098e-04_jprb, &
     & 1.00316e-03_jprb, 1.11927e-03_jprb, 1.24881e-03_jprb, 1.39335e-03_jprb, 1.55461e-03_jprb, &
     & 1.73455e-03_jprb, 1.93530e-03_jprb, 2.15930e-03_jprb, 2.40921e-03_jprb/)
      kao_mco2( 7, :, 3) = (/ &
     & 3.24677e-04_jprb, 3.75160e-04_jprb, 4.33491e-04_jprb, 5.00893e-04_jprb, 5.78774e-04_jprb, &
     & 6.68765e-04_jprb, 7.72749e-04_jprb, 8.92900e-04_jprb, 1.03173e-03_jprb, 1.19215e-03_jprb, &
     & 1.37751e-03_jprb, 1.59170e-03_jprb, 1.83918e-03_jprb, 2.12515e-03_jprb, 2.45558e-03_jprb, &
     & 2.83738e-03_jprb, 3.27856e-03_jprb, 3.78832e-03_jprb, 4.37735e-03_jprb/)
      kao_mco2( 8, :, 3) = (/ &
     & 2.24656e-03_jprb, 2.45550e-03_jprb, 2.68386e-03_jprb, 2.93347e-03_jprb, 3.20629e-03_jprb, &
     & 3.50448e-03_jprb, 3.83041e-03_jprb, 4.18665e-03_jprb, 4.57602e-03_jprb, 5.00160e-03_jprb, &
     & 5.46677e-03_jprb, 5.97519e-03_jprb, 6.53090e-03_jprb, 7.13829e-03_jprb, 7.80217e-03_jprb, &
     & 8.52780e-03_jprb, 9.32091e-03_jprb, 1.01878e-02_jprb, 1.11353e-02_jprb/)
      kao_mco2( 9, :, 3) = (/ &
     & 2.07746e-04_jprb, 2.38909e-04_jprb, 2.74746e-04_jprb, 3.15959e-04_jprb, 3.63355e-04_jprb, &
     & 4.17860e-04_jprb, 4.80541e-04_jprb, 5.52625e-04_jprb, 6.35521e-04_jprb, 7.30852e-04_jprb, &
     & 8.40484e-04_jprb, 9.66561e-04_jprb, 1.11155e-03_jprb, 1.27829e-03_jprb, 1.47004e-03_jprb, &
     & 1.69055e-03_jprb, 1.94414e-03_jprb, 2.23577e-03_jprb, 2.57115e-03_jprb/)
      kao_mco2( 1, :, 4) = (/ &
     & 7.26052e-04_jprb, 7.62476e-04_jprb, 8.00726e-04_jprb, 8.40896e-04_jprb, 8.83081e-04_jprb, &
     & 9.27382e-04_jprb, 9.73905e-04_jprb, 1.02276e-03_jprb, 1.07407e-03_jprb, 1.12795e-03_jprb, &
     & 1.18454e-03_jprb, 1.24396e-03_jprb, 1.30637e-03_jprb, 1.37190e-03_jprb, 1.44073e-03_jprb, &
     & 1.51300e-03_jprb, 1.58890e-03_jprb, 1.66861e-03_jprb, 1.75232e-03_jprb/)
      kao_mco2( 2, :, 4) = (/ &
     & 4.65815e-04_jprb, 5.01167e-04_jprb, 5.39203e-04_jprb, 5.80126e-04_jprb, 6.24154e-04_jprb, &
     & 6.71524e-04_jprb, 7.22489e-04_jprb, 7.77322e-04_jprb, 8.36316e-04_jprb, 8.99788e-04_jprb, &
     & 9.68077e-04_jprb, 1.04155e-03_jprb, 1.12060e-03_jprb, 1.20564e-03_jprb, 1.29714e-03_jprb, &
     & 1.39559e-03_jprb, 1.50151e-03_jprb, 1.61546e-03_jprb, 1.73807e-03_jprb/)
      kao_mco2( 3, :, 4) = (/ &
     & 3.56225e-04_jprb, 3.93073e-04_jprb, 4.33732e-04_jprb, 4.78598e-04_jprb, 5.28105e-04_jprb, &
     & 5.82732e-04_jprb, 6.43010e-04_jprb, 7.09524e-04_jprb, 7.82918e-04_jprb, 8.63903e-04_jprb, &
     & 9.53266e-04_jprb, 1.05187e-03_jprb, 1.16068e-03_jprb, 1.28074e-03_jprb, 1.41322e-03_jprb, &
     & 1.55941e-03_jprb, 1.72071e-03_jprb, 1.89870e-03_jprb, 2.09511e-03_jprb/)
      kao_mco2( 4, :, 4) = (/ &
     & 3.37845e-04_jprb, 3.79675e-04_jprb, 4.26684e-04_jprb, 4.79514e-04_jprb, 5.38884e-04_jprb, &
     & 6.05606e-04_jprb, 6.80589e-04_jprb, 7.64855e-04_jprb, 8.59555e-04_jprb, 9.65980e-04_jprb, &
     & 1.08558e-03_jprb, 1.21999e-03_jprb, 1.37105e-03_jprb, 1.54080e-03_jprb, 1.73157e-03_jprb, &
     & 1.94597e-03_jprb, 2.18691e-03_jprb, 2.45767e-03_jprb, 2.76197e-03_jprb/)
      kao_mco2( 5, :, 4) = (/ &
     & 3.52456e-04_jprb, 4.02782e-04_jprb, 4.60294e-04_jprb, 5.26017e-04_jprb, 6.01126e-04_jprb, &
     & 6.86958e-04_jprb, 7.85046e-04_jprb, 8.97140e-04_jprb, 1.02524e-03_jprb, 1.17163e-03_jprb, &
     & 1.33892e-03_jprb, 1.53010e-03_jprb, 1.74858e-03_jprb, 1.99825e-03_jprb, 2.28358e-03_jprb, &
     & 2.60964e-03_jprb, 2.98226e-03_jprb, 3.40809e-03_jprb, 3.89471e-03_jprb/)
      kao_mco2( 6, :, 4) = (/ &
     & 4.42884e-04_jprb, 5.08187e-04_jprb, 5.83119e-04_jprb, 6.69100e-04_jprb, 7.67758e-04_jprb, &
     & 8.80963e-04_jprb, 1.01086e-03_jprb, 1.15991e-03_jprb, 1.33094e-03_jprb, 1.52718e-03_jprb, &
     & 1.75237e-03_jprb, 2.01075e-03_jprb, 2.30724e-03_jprb, 2.64744e-03_jprb, 3.03780e-03_jprb, &
     & 3.48572e-03_jprb, 3.99969e-03_jprb, 4.58944e-03_jprb, 5.26614e-03_jprb/)
      kao_mco2( 7, :, 4) = (/ &
     & 8.09850e-04_jprb, 9.09940e-04_jprb, 1.02240e-03_jprb, 1.14876e-03_jprb, 1.29074e-03_jprb, &
     & 1.45026e-03_jprb, 1.62950e-03_jprb, 1.83089e-03_jprb, 2.05718e-03_jprb, 2.31143e-03_jprb, &
     & 2.59710e-03_jprb, 2.91808e-03_jprb, 3.27873e-03_jprb, 3.68395e-03_jprb, 4.13926e-03_jprb, &
     & 4.65083e-03_jprb, 5.22564e-03_jprb, 5.87148e-03_jprb, 6.59715e-03_jprb/)
      kao_mco2( 8, :, 4) = (/ &
     & 3.13265e-03_jprb, 3.42454e-03_jprb, 3.74362e-03_jprb, 4.09243e-03_jprb, 4.47375e-03_jprb, &
     & 4.89059e-03_jprb, 5.34627e-03_jprb, 5.84441e-03_jprb, 6.38897e-03_jprb, 6.98426e-03_jprb, &
     & 7.63502e-03_jprb, 8.34642e-03_jprb, 9.12409e-03_jprb, 9.97423e-03_jprb, 1.09036e-02_jprb, &
     & 1.19195e-02_jprb, 1.30301e-02_jprb, 1.42442e-02_jprb, 1.55714e-02_jprb/)
      kao_mco2( 9, :, 4) = (/ &
     & 5.71287e-04_jprb, 6.51252e-04_jprb, 7.42411e-04_jprb, 8.46330e-04_jprb, 9.64794e-04_jprb, &
     & 1.09984e-03_jprb, 1.25379e-03_jprb, 1.42929e-03_jprb, 1.62935e-03_jprb, 1.85742e-03_jprb, &
     & 2.11741e-03_jprb, 2.41380e-03_jprb, 2.75167e-03_jprb, 3.13683e-03_jprb, 3.57591e-03_jprb, &
     & 4.07645e-03_jprb, 4.64705e-03_jprb, 5.29751e-03_jprb, 6.03903e-03_jprb/)
      kao_mco2( 1, :, 5) = (/ &
     & 2.92395e-04_jprb, 3.32719e-04_jprb, 3.78604e-04_jprb, 4.30818e-04_jprb, 4.90232e-04_jprb, &
     & 5.57839e-04_jprb, 6.34771e-04_jprb, 7.22312e-04_jprb, 8.21927e-04_jprb, 9.35278e-04_jprb, &
     & 1.06426e-03_jprb, 1.21104e-03_jprb, 1.37805e-03_jprb, 1.56810e-03_jprb, 1.78435e-03_jprb, &
     & 2.03043e-03_jprb, 2.31045e-03_jprb, 2.62908e-03_jprb, 2.99166e-03_jprb/)
      kao_mco2( 2, :, 5) = (/ &
     & 3.13069e-04_jprb, 3.61343e-04_jprb, 4.17061e-04_jprb, 4.81371e-04_jprb, 5.55597e-04_jprb, &
     & 6.41269e-04_jprb, 7.40151e-04_jprb, 8.54280e-04_jprb, 9.86008e-04_jprb, 1.13805e-03_jprb, &
     & 1.31353e-03_jprb, 1.51608e-03_jprb, 1.74985e-03_jprb, 2.01967e-03_jprb, 2.33110e-03_jprb, &
     & 2.69055e-03_jprb, 3.10543e-03_jprb, 3.58427e-03_jprb, 4.13696e-03_jprb/)
      kao_mco2( 3, :, 5) = (/ &
     & 3.06937e-04_jprb, 3.57841e-04_jprb, 4.17187e-04_jprb, 4.86375e-04_jprb, 5.67038e-04_jprb, &
     & 6.61078e-04_jprb, 7.70714e-04_jprb, 8.98532e-04_jprb, 1.04755e-03_jprb, 1.22128e-03_jprb, &
     & 1.42382e-03_jprb, 1.65996e-03_jprb, 1.93525e-03_jprb, 2.25620e-03_jprb, 2.63038e-03_jprb, &
     & 3.06661e-03_jprb, 3.57519e-03_jprb, 4.16812e-03_jprb, 4.85937e-03_jprb/)
      kao_mco2( 4, :, 5) = (/ &
     & 4.06428e-04_jprb, 4.72379e-04_jprb, 5.49033e-04_jprb, 6.38125e-04_jprb, 7.41674e-04_jprb, &
     & 8.62026e-04_jprb, 1.00191e-03_jprb, 1.16449e-03_jprb, 1.35345e-03_jprb, 1.57308e-03_jprb, &
     & 1.82834e-03_jprb, 2.12503e-03_jprb, 2.46986e-03_jprb, 2.87064e-03_jprb, 3.33647e-03_jprb, &
     & 3.87788e-03_jprb, 4.50715e-03_jprb, 5.23852e-03_jprb, 6.08858e-03_jprb/)
      kao_mco2( 5, :, 5) = (/ &
     & 6.01967e-04_jprb, 6.90414e-04_jprb, 7.91856e-04_jprb, 9.08204e-04_jprb, 1.04165e-03_jprb, &
     & 1.19470e-03_jprb, 1.37023e-03_jprb, 1.57156e-03_jprb, 1.80247e-03_jprb, 2.06731e-03_jprb, &
     & 2.37106e-03_jprb, 2.71944e-03_jprb, 3.11901e-03_jprb, 3.57729e-03_jprb, 4.10290e-03_jprb, &
     & 4.70574e-03_jprb, 5.39716e-03_jprb, 6.19017e-03_jprb, 7.09969e-03_jprb/)
      kao_mco2( 6, :, 5) = (/ &
     & 1.11622e-03_jprb, 1.25799e-03_jprb, 1.41776e-03_jprb, 1.59783e-03_jprb, 1.80077e-03_jprb, &
     & 2.02947e-03_jprb, 2.28723e-03_jprb, 2.57773e-03_jprb, 2.90512e-03_jprb, 3.27408e-03_jprb, &
     & 3.68992e-03_jprb, 4.15856e-03_jprb, 4.68673e-03_jprb, 5.28197e-03_jprb, 5.95282e-03_jprb, &
     & 6.70887e-03_jprb, 7.56094e-03_jprb, 8.52123e-03_jprb, 9.60348e-03_jprb/)
      kao_mco2( 7, :, 5) = (/ &
     & 3.63860e-03_jprb, 3.96164e-03_jprb, 4.31337e-03_jprb, 4.69632e-03_jprb, 5.11327e-03_jprb, &
     & 5.56724e-03_jprb, 6.06151e-03_jprb, 6.59967e-03_jprb, 7.18561e-03_jprb, 7.82356e-03_jprb, &
     & 8.51816e-03_jprb, 9.27443e-03_jprb, 1.00978e-02_jprb, 1.09943e-02_jprb, 1.19705e-02_jprb, &
     & 1.30332e-02_jprb, 1.41904e-02_jprb, 1.54502e-02_jprb, 1.68219e-02_jprb/)
      kao_mco2( 8, :, 5) = (/ &
     & 5.96957e-03_jprb, 6.53049e-03_jprb, 7.14412e-03_jprb, 7.81541e-03_jprb, 8.54977e-03_jprb, &
     & 9.35314e-03_jprb, 1.02320e-02_jprb, 1.11934e-02_jprb, 1.22452e-02_jprb, 1.33958e-02_jprb, &
     & 1.46545e-02_jprb, 1.60315e-02_jprb, 1.75379e-02_jprb, 1.91858e-02_jprb, 2.09886e-02_jprb, &
     & 2.29608e-02_jprb, 2.51182e-02_jprb, 2.74784e-02_jprb, 3.00604e-02_jprb/)
      kao_mco2( 9, :, 5) = (/ &
     & 1.19381e-03_jprb, 1.33882e-03_jprb, 1.50143e-03_jprb, 1.68379e-03_jprb, 1.88831e-03_jprb, &
     & 2.11767e-03_jprb, 2.37488e-03_jprb, 2.66333e-03_jprb, 2.98683e-03_jprb, 3.34961e-03_jprb, &
     & 3.75646e-03_jprb, 4.21272e-03_jprb, 4.72440e-03_jprb, 5.29823e-03_jprb, 5.94176e-03_jprb, &
     & 6.66345e-03_jprb, 7.47280e-03_jprb, 8.38045e-03_jprb, 9.39835e-03_jprb/)
      kao_mco2( 1, :, 6) = (/ &
     & 4.12429e-04_jprb, 4.84830e-04_jprb, 5.69942e-04_jprb, 6.69995e-04_jprb, 7.87613e-04_jprb, &
     & 9.25878e-04_jprb, 1.08842e-03_jprb, 1.27949e-03_jprb, 1.50410e-03_jprb, 1.76814e-03_jprb, &
     & 2.07854e-03_jprb, 2.44343e-03_jprb, 2.87237e-03_jprb, 3.37662e-03_jprb, 3.96938e-03_jprb, &
     & 4.66621e-03_jprb, 5.48536e-03_jprb, 6.44831e-03_jprb, 7.58031e-03_jprb/)
      kao_mco2( 2, :, 6) = (/ &
     & 6.43498e-04_jprb, 7.46132e-04_jprb, 8.65134e-04_jprb, 1.00312e-03_jprb, 1.16311e-03_jprb, &
     & 1.34861e-03_jprb, 1.56371e-03_jprb, 1.81310e-03_jprb, 2.10228e-03_jprb, 2.43758e-03_jprb, &
     & 2.82635e-03_jprb, 3.27714e-03_jprb, 3.79981e-03_jprb, 4.40586e-03_jprb, 5.10855e-03_jprb, &
     & 5.92333e-03_jprb, 6.86806e-03_jprb, 7.96346e-03_jprb, 9.23357e-03_jprb/)
      kao_mco2( 3, :, 6) = (/ &
     & 1.11336e-03_jprb, 1.26910e-03_jprb, 1.44662e-03_jprb, 1.64897e-03_jprb, 1.87962e-03_jprb, &
     & 2.14254e-03_jprb, 2.44224e-03_jprb, 2.78385e-03_jprb, 3.17325e-03_jprb, 3.61712e-03_jprb, &
     & 4.12308e-03_jprb, 4.69981e-03_jprb, 5.35720e-03_jprb, 6.10656e-03_jprb, 6.96073e-03_jprb, &
     & 7.93439e-03_jprb, 9.04424e-03_jprb, 1.03093e-02_jprb, 1.17514e-02_jprb/)
      kao_mco2( 4, :, 6) = (/ &
     & 1.87991e-03_jprb, 2.10276e-03_jprb, 2.35202e-03_jprb, 2.63082e-03_jprb, 2.94268e-03_jprb, &
     & 3.29150e-03_jprb, 3.68168e-03_jprb, 4.11810e-03_jprb, 4.60626e-03_jprb, 5.15228e-03_jprb, &
     & 5.76303e-03_jprb, 6.44617e-03_jprb, 7.21030e-03_jprb, 8.06500e-03_jprb, 9.02102e-03_jprb, &
     & 1.00904e-02_jprb, 1.12865e-02_jprb, 1.26244e-02_jprb, 1.41208e-02_jprb/)
      kao_mco2( 5, :, 6) = (/ &
     & 3.65848e-03_jprb, 4.01372e-03_jprb, 4.40346e-03_jprb, 4.83104e-03_jprb, 5.30015e-03_jprb, &
     & 5.81480e-03_jprb, 6.37943e-03_jprb, 6.99888e-03_jprb, 7.67849e-03_jprb, 8.42408e-03_jprb, &
     & 9.24208e-03_jprb, 1.01395e-02_jprb, 1.11241e-02_jprb, 1.22042e-02_jprb, 1.33893e-02_jprb, &
     & 1.46894e-02_jprb, 1.61158e-02_jprb, 1.76806e-02_jprb, 1.93975e-02_jprb/)
      kao_mco2( 6, :, 6) = (/ &
     & 5.38476e-03_jprb, 5.85088e-03_jprb, 6.35735e-03_jprb, 6.90765e-03_jprb, 7.50560e-03_jprb, &
     & 8.15530e-03_jprb, 8.86124e-03_jprb, 9.62829e-03_jprb, 1.04617e-02_jprb, 1.13673e-02_jprb, &
     & 1.23513e-02_jprb, 1.34205e-02_jprb, 1.45822e-02_jprb, 1.58445e-02_jprb, 1.72160e-02_jprb, &
     & 1.87062e-02_jprb, 2.03255e-02_jprb, 2.20849e-02_jprb, 2.39966e-02_jprb/)
      kao_mco2( 7, :, 6) = (/ &
     & 6.27017e-03_jprb, 6.84772e-03_jprb, 7.47846e-03_jprb, 8.16731e-03_jprb, 8.91960e-03_jprb, &
     & 9.74118e-03_jprb, 1.06384e-02_jprb, 1.16183e-02_jprb, 1.26885e-02_jprb, 1.38573e-02_jprb, &
     & 1.51336e-02_jprb, 1.65276e-02_jprb, 1.80500e-02_jprb, 1.97125e-02_jprb, 2.15283e-02_jprb, &
     & 2.35112e-02_jprb, 2.56769e-02_jprb, 2.80420e-02_jprb, 3.06249e-02_jprb/)
      kao_mco2( 8, :, 6) = (/ &
     & 9.61932e-03_jprb, 1.04802e-02_jprb, 1.14182e-02_jprb, 1.24401e-02_jprb, 1.35534e-02_jprb, &
     & 1.47664e-02_jprb, 1.60880e-02_jprb, 1.75278e-02_jprb, 1.90965e-02_jprb, 2.08056e-02_jprb, &
     & 2.26677e-02_jprb, 2.46964e-02_jprb, 2.69066e-02_jprb, 2.93147e-02_jprb, 3.19383e-02_jprb, &
     & 3.47967e-02_jprb, 3.79110e-02_jprb, 4.13039e-02_jprb, 4.50005e-02_jprb/)
      kao_mco2( 9, :, 6) = (/ &
     & 2.37921e-03_jprb, 2.64556e-03_jprb, 2.94173e-03_jprb, 3.27105e-03_jprb, 3.63724e-03_jprb, &
     & 4.04442e-03_jprb, 4.49718e-03_jprb, 5.00064e-03_jprb, 5.56045e-03_jprb, 6.18293e-03_jprb, &
     & 6.87510e-03_jprb, 7.64475e-03_jprb, 8.50057e-03_jprb, 9.45219e-03_jprb, 1.05103e-02_jprb, &
     & 1.16870e-02_jprb, 1.29953e-02_jprb, 1.44501e-02_jprb, 1.60677e-02_jprb/)
      kao_mco2( 1, :, 7) = (/ &
     & 4.64970e-03_jprb, 5.13188e-03_jprb, 5.66406e-03_jprb, 6.25144e-03_jprb, 6.89972e-03_jprb, &
     & 7.61523e-03_jprb, 8.40493e-03_jprb, 9.27654e-03_jprb, 1.02385e-02_jprb, 1.13003e-02_jprb, &
     & 1.24721e-02_jprb, 1.37655e-02_jprb, 1.51930e-02_jprb, 1.67685e-02_jprb, 1.85075e-02_jprb, &
     & 2.04267e-02_jprb, 2.25450e-02_jprb, 2.48829e-02_jprb, 2.74633e-02_jprb/)
      kao_mco2( 2, :, 7) = (/ &
     & 6.37148e-03_jprb, 6.96805e-03_jprb, 7.62046e-03_jprb, 8.33397e-03_jprb, 9.11428e-03_jprb, &
     & 9.96765e-03_jprb, 1.09009e-02_jprb, 1.19216e-02_jprb, 1.30378e-02_jprb, 1.42585e-02_jprb, &
     & 1.55935e-02_jprb, 1.70536e-02_jprb, 1.86503e-02_jprb, 2.03965e-02_jprb, 2.23062e-02_jprb, &
     & 2.43948e-02_jprb, 2.66789e-02_jprb, 2.91768e-02_jprb, 3.19086e-02_jprb/)
      kao_mco2( 3, :, 7) = (/ &
     & 7.79364e-03_jprb, 8.48097e-03_jprb, 9.22892e-03_jprb, 1.00428e-02_jprb, 1.09285e-02_jprb, &
     & 1.18923e-02_jprb, 1.29411e-02_jprb, 1.40825e-02_jprb, 1.53244e-02_jprb, 1.66759e-02_jprb, &
     & 1.81466e-02_jprb, 1.97470e-02_jprb, 2.14885e-02_jprb, 2.33836e-02_jprb, 2.54458e-02_jprb, &
     & 2.76899e-02_jprb, 3.01320e-02_jprb, 3.27893e-02_jprb, 3.56811e-02_jprb/)
      kao_mco2( 4, :, 7) = (/ &
     & 8.70586e-03_jprb, 9.48737e-03_jprb, 1.03390e-02_jprb, 1.12672e-02_jprb, 1.22786e-02_jprb, &
     & 1.33808e-02_jprb, 1.45820e-02_jprb, 1.58910e-02_jprb, 1.73175e-02_jprb, 1.88721e-02_jprb, &
     & 2.05662e-02_jprb, 2.24124e-02_jprb, 2.44243e-02_jprb, 2.66169e-02_jprb, 2.90062e-02_jprb, &
     & 3.16101e-02_jprb, 3.44477e-02_jprb, 3.75400e-02_jprb, 4.09099e-02_jprb/)
      kao_mco2( 5, :, 7) = (/ &
     & 9.24510e-03_jprb, 1.00865e-02_jprb, 1.10045e-02_jprb, 1.20061e-02_jprb, 1.30988e-02_jprb, &
     & 1.42910e-02_jprb, 1.55916e-02_jprb, 1.70106e-02_jprb, 1.85588e-02_jprb, 2.02479e-02_jprb, &
     & 2.20908e-02_jprb, 2.41013e-02_jprb, 2.62948e-02_jprb, 2.86880e-02_jprb, 3.12990e-02_jprb, &
     & 3.41476e-02_jprb, 3.72555e-02_jprb, 4.06462e-02_jprb, 4.43455e-02_jprb/)
      kao_mco2( 6, :, 7) = (/ &
     & 1.09559e-02_jprb, 1.19933e-02_jprb, 1.31290e-02_jprb, 1.43722e-02_jprb, 1.57331e-02_jprb, &
     & 1.72229e-02_jprb, 1.88537e-02_jprb, 2.06390e-02_jprb, 2.25933e-02_jprb, 2.47327e-02_jprb, &
     & 2.70747e-02_jprb, 2.96384e-02_jprb, 3.24449e-02_jprb, 3.55171e-02_jprb, 3.88802e-02_jprb, &
     & 4.25619e-02_jprb, 4.65921e-02_jprb, 5.10039e-02_jprb, 5.58335e-02_jprb/)
      kao_mco2( 7, :, 7) = (/ &
     & 1.36116e-02_jprb, 1.48659e-02_jprb, 1.62357e-02_jprb, 1.77318e-02_jprb, 1.93657e-02_jprb, &
     & 2.11502e-02_jprb, 2.30991e-02_jprb, 2.52276e-02_jprb, 2.75522e-02_jprb, 3.00910e-02_jprb, &
     & 3.28638e-02_jprb, 3.58921e-02_jprb, 3.91995e-02_jprb, 4.28116e-02_jprb, 4.67565e-02_jprb, &
     & 5.10650e-02_jprb, 5.57704e-02_jprb, 6.09095e-02_jprb, 6.65221e-02_jprb/)
      kao_mco2( 8, :, 7) = (/ &
     & 1.51783e-02_jprb, 1.64551e-02_jprb, 1.78392e-02_jprb, 1.93399e-02_jprb, 2.09667e-02_jprb, &
     & 2.27304e-02_jprb, 2.46424e-02_jprb, 2.67153e-02_jprb, 2.89626e-02_jprb, 3.13988e-02_jprb, &
     & 3.40401e-02_jprb, 3.69035e-02_jprb, 4.00077e-02_jprb, 4.33731e-02_jprb, 4.70216e-02_jprb, &
     & 5.09770e-02_jprb, 5.52651e-02_jprb, 5.99139e-02_jprb, 6.49538e-02_jprb/)
      kao_mco2( 9, :, 7) = (/ &
     & 1.00072e-02_jprb, 1.08638e-02_jprb, 1.17937e-02_jprb, 1.28032e-02_jprb, 1.38991e-02_jprb, &
     & 1.50888e-02_jprb, 1.63803e-02_jprb, 1.77824e-02_jprb, 1.93045e-02_jprb, 2.09568e-02_jprb, &
     & 2.27507e-02_jprb, 2.46980e-02_jprb, 2.68120e-02_jprb, 2.91070e-02_jprb, 3.15984e-02_jprb, &
     & 3.43031e-02_jprb, 3.72393e-02_jprb, 4.04268e-02_jprb, 4.38872e-02_jprb/)
      kao_mco2( 1, :, 8) = (/ &
     & 1.59610e-02_jprb, 1.74387e-02_jprb, 1.90532e-02_jprb, 2.08171e-02_jprb, 2.27444e-02_jprb, &
     & 2.48501e-02_jprb, 2.71508e-02_jprb, 2.96645e-02_jprb, 3.24109e-02_jprb, 3.54115e-02_jprb, &
     & 3.86900e-02_jprb, 4.22720e-02_jprb, 4.61856e-02_jprb, 5.04616e-02_jprb, 5.51334e-02_jprb, &
     & 6.02378e-02_jprb, 6.58147e-02_jprb, 7.19079e-02_jprb, 7.85653e-02_jprb/)
      kao_mco2( 2, :, 8) = (/ &
     & 1.61961e-02_jprb, 1.76986e-02_jprb, 1.93405e-02_jprb, 2.11348e-02_jprb, 2.30955e-02_jprb, &
     & 2.52381e-02_jprb, 2.75794e-02_jprb, 3.01380e-02_jprb, 3.29340e-02_jprb, 3.59893e-02_jprb, &
     & 3.93280e-02_jprb, 4.29766e-02_jprb, 4.69636e-02_jprb, 5.13204e-02_jprb, 5.60815e-02_jprb, &
     & 6.12843e-02_jprb, 6.69697e-02_jprb, 7.31826e-02_jprb, 7.99718e-02_jprb/)
      kao_mco2( 3, :, 8) = (/ &
     & 1.72034e-02_jprb, 1.88241e-02_jprb, 2.05974e-02_jprb, 2.25377e-02_jprb, 2.46609e-02_jprb, &
     & 2.69841e-02_jprb, 2.95261e-02_jprb, 3.23076e-02_jprb, 3.53511e-02_jprb, 3.86813e-02_jprb, &
     & 4.23253e-02_jprb, 4.63126e-02_jprb, 5.06754e-02_jprb, 5.54493e-02_jprb, 6.06728e-02_jprb, &
     & 6.63885e-02_jprb, 7.26426e-02_jprb, 7.94859e-02_jprb, 8.69738e-02_jprb/)
      kao_mco2( 4, :, 8) = (/ &
     & 1.79777e-02_jprb, 1.96517e-02_jprb, 2.14815e-02_jprb, 2.34817e-02_jprb, 2.56682e-02_jprb, &
     & 2.80583e-02_jprb, 3.06709e-02_jprb, 3.35268e-02_jprb, 3.66486e-02_jprb, 4.00611e-02_jprb, &
     & 4.37914e-02_jprb, 4.78690e-02_jprb, 5.23262e-02_jprb, 5.71985e-02_jprb, 6.25245e-02_jprb, &
     & 6.83464e-02_jprb, 7.47104e-02_jprb, 8.16670e-02_jprb, 8.92713e-02_jprb/)
      kao_mco2( 5, :, 8) = (/ &
     & 2.02540e-02_jprb, 2.21214e-02_jprb, 2.41610e-02_jprb, 2.63887e-02_jprb, 2.88218e-02_jprb, &
     & 3.14792e-02_jprb, 3.43816e-02_jprb, 3.75516e-02_jprb, 4.10139e-02_jprb, 4.47954e-02_jprb, &
     & 4.89256e-02_jprb, 5.34366e-02_jprb, 5.83635e-02_jprb, 6.37447e-02_jprb, 6.96220e-02_jprb, &
     & 7.60413e-02_jprb, 8.30523e-02_jprb, 9.07098e-02_jprb, 9.90734e-02_jprb/)
      kao_mco2( 6, :, 8) = (/ &
     & 2.19009e-02_jprb, 2.38517e-02_jprb, 2.59762e-02_jprb, 2.82899e-02_jprb, 3.08097e-02_jprb, &
     & 3.35540e-02_jprb, 3.65427e-02_jprb, 3.97976e-02_jprb, 4.33424e-02_jprb, 4.72030e-02_jprb, &
     & 5.14074e-02_jprb, 5.59863e-02_jprb, 6.09731e-02_jprb, 6.64040e-02_jprb, 7.23187e-02_jprb, &
     & 7.87603e-02_jprb, 8.57755e-02_jprb, 9.34157e-02_jprb, 1.01736e-01_jprb/)
      kao_mco2( 7, :, 8) = (/ &
     & 2.52383e-02_jprb, 2.73978e-02_jprb, 2.97421e-02_jprb, 3.22869e-02_jprb, 3.50496e-02_jprb, &
     & 3.80486e-02_jprb, 4.13042e-02_jprb, 4.48383e-02_jprb, 4.86749e-02_jprb, 5.28397e-02_jprb, &
     & 5.73610e-02_jprb, 6.22690e-02_jprb, 6.75970e-02_jprb, 7.33810e-02_jprb, 7.96598e-02_jprb, &
     & 8.64758e-02_jprb, 9.38751e-02_jprb, 1.01907e-01_jprb, 1.10627e-01_jprb/)
      kao_mco2( 8, :, 8) = (/ &
     & 3.36506e-02_jprb, 3.59288e-02_jprb, 3.83613e-02_jprb, 4.09584e-02_jprb, 4.37313e-02_jprb, &
     & 4.66920e-02_jprb, 4.98531e-02_jprb, 5.32283e-02_jprb, 5.68319e-02_jprb, 6.06795e-02_jprb, &
     & 6.47876e-02_jprb, 6.91739e-02_jprb, 7.38570e-02_jprb, 7.88573e-02_jprb, 8.41960e-02_jprb, &
     & 8.98962e-02_jprb, 9.59824e-02_jprb, 1.02481e-01_jprb, 1.09419e-01_jprb/)
      kao_mco2( 9, :, 8) = (/ &
     & 2.15151e-02_jprb, 2.34420e-02_jprb, 2.55415e-02_jprb, 2.78291e-02_jprb, 3.03215e-02_jprb, &
     & 3.30372e-02_jprb, 3.59961e-02_jprb, 3.92200e-02_jprb, 4.27326e-02_jprb, 4.65598e-02_jprb, &
     & 5.07299e-02_jprb, 5.52734e-02_jprb, 6.02238e-02_jprb, 6.56176e-02_jprb, 7.14944e-02_jprb, &
     & 7.78977e-02_jprb, 8.48744e-02_jprb, 9.24759e-02_jprb, 1.00758e-01_jprb/)
      kao_mco2( 1, :, 9) = (/ &
     & 3.34296e-02_jprb, 3.64437e-02_jprb, 3.97294e-02_jprb, 4.33114e-02_jprb, 4.72164e-02_jprb, &
     & 5.14734e-02_jprb, 5.61143e-02_jprb, 6.11735e-02_jprb, 6.66890e-02_jprb, 7.27016e-02_jprb, &
     & 7.92564e-02_jprb, 8.64022e-02_jprb, 9.41922e-02_jprb, 1.02685e-01_jprb, 1.11943e-01_jprb, &
     & 1.22035e-01_jprb, 1.33038e-01_jprb, 1.45033e-01_jprb, 1.58109e-01_jprb/)
      kao_mco2( 2, :, 9) = (/ &
     & 3.73946e-02_jprb, 4.07543e-02_jprb, 4.44160e-02_jprb, 4.84066e-02_jprb, 5.27558e-02_jprb, &
     & 5.74958e-02_jprb, 6.26616e-02_jprb, 6.82915e-02_jprb, 7.44273e-02_jprb, 8.11144e-02_jprb, &
     & 8.84023e-02_jprb, 9.63449e-02_jprb, 1.05001e-01_jprb, 1.14435e-01_jprb, 1.24717e-01_jprb, &
     & 1.35922e-01_jprb, 1.48135e-01_jprb, 1.61444e-01_jprb, 1.75949e-01_jprb/)
      kao_mco2( 3, :, 9) = (/ &
     & 4.24539e-02_jprb, 4.61192e-02_jprb, 5.01010e-02_jprb, 5.44265e-02_jprb, 5.91255e-02_jprb, &
     & 6.42302e-02_jprb, 6.97756e-02_jprb, 7.57998e-02_jprb, 8.23442e-02_jprb, 8.94535e-02_jprb, &
     & 9.71766e-02_jprb, 1.05566e-01_jprb, 1.14681e-01_jprb, 1.24582e-01_jprb, 1.35338e-01_jprb, &
     & 1.47022e-01_jprb, 1.59716e-01_jprb, 1.73505e-01_jprb, 1.88485e-01_jprb/)
      kao_mco2( 4, :, 9) = (/ &
     & 5.30296e-02_jprb, 5.73416e-02_jprb, 6.20043e-02_jprb, 6.70462e-02_jprb, 7.24980e-02_jprb, &
     & 7.83931e-02_jprb, 8.47676e-02_jprb, 9.16604e-02_jprb, 9.91137e-02_jprb, 1.07173e-01_jprb, &
     & 1.15888e-01_jprb, 1.25311e-01_jprb, 1.35501e-01_jprb, 1.46519e-01_jprb, 1.58433e-01_jprb, &
     & 1.71316e-01_jprb, 1.85246e-01_jprb, 2.00309e-01_jprb, 2.16597e-01_jprb/)
      kao_mco2( 5, :, 9) = (/ &
     & 6.26111e-02_jprb, 6.74018e-02_jprb, 7.25591e-02_jprb, 7.81111e-02_jprb, 8.40878e-02_jprb, &
     & 9.05218e-02_jprb, 9.74482e-02_jprb, 1.04904e-01_jprb, 1.12931e-01_jprb, 1.21572e-01_jprb, &
     & 1.30875e-01_jprb, 1.40889e-01_jprb, 1.51669e-01_jprb, 1.63274e-01_jprb, 1.75767e-01_jprb, &
     & 1.89216e-01_jprb, 2.03694e-01_jprb, 2.19279e-01_jprb, 2.36058e-01_jprb/)
      kao_mco2( 6, :, 9) = (/ &
     & 7.59080e-02_jprb, 8.13446e-02_jprb, 8.71706e-02_jprb, 9.34139e-02_jprb, 1.00104e-01_jprb, &
     & 1.07274e-01_jprb, 1.14957e-01_jprb, 1.23190e-01_jprb, 1.32013e-01_jprb, 1.41468e-01_jprb, &
     & 1.51600e-01_jprb, 1.62458e-01_jprb, 1.74094e-01_jprb, 1.86562e-01_jprb, 1.99924e-01_jprb, &
     & 2.14243e-01_jprb, 2.29587e-01_jprb, 2.46031e-01_jprb, 2.63652e-01_jprb/)
      kao_mco2( 7, :, 9) = (/ &
     & 8.81942e-02_jprb, 9.39942e-02_jprb, 1.00176e-01_jprb, 1.06763e-01_jprb, 1.13784e-01_jprb, &
     & 1.21267e-01_jprb, 1.29242e-01_jprb, 1.37742e-01_jprb, 1.46800e-01_jprb, 1.56454e-01_jprb, &
     & 1.66743e-01_jprb, 1.77708e-01_jprb, 1.89395e-01_jprb, 2.01850e-01_jprb, 2.15124e-01_jprb, &
     & 2.29272e-01_jprb, 2.44349e-01_jprb, 2.60418e-01_jprb, 2.77544e-01_jprb/)
      kao_mco2( 8, :, 9) = (/ &
     & 6.28535e-02_jprb, 6.69314e-02_jprb, 7.12740e-02_jprb, 7.58982e-02_jprb, 8.08225e-02_jprb, &
     & 8.60662e-02_jprb, 9.16502e-02_jprb, 9.75965e-02_jprb, 1.03929e-01_jprb, 1.10671e-01_jprb, &
     & 1.17852e-01_jprb, 1.25498e-01_jprb, 1.33640e-01_jprb, 1.42311e-01_jprb, 1.51544e-01_jprb, &
     & 1.61376e-01_jprb, 1.71846e-01_jprb, 1.82996e-01_jprb, 1.94868e-01_jprb/)
      kao_mco2( 9, :, 9) = (/ &
     & 6.39196e-02_jprb, 6.86702e-02_jprb, 7.37738e-02_jprb, 7.92568e-02_jprb, 8.51473e-02_jprb, &
     & 9.14756e-02_jprb, 9.82742e-02_jprb, 1.05578e-01_jprb, 1.13425e-01_jprb, 1.21855e-01_jprb, &
     & 1.30911e-01_jprb, 1.40641e-01_jprb, 1.51093e-01_jprb, 1.62323e-01_jprb, 1.74387e-01_jprb, &
     & 1.87348e-01_jprb, 2.01272e-01_jprb, 2.16231e-01_jprb, 2.32301e-01_jprb/)
      kao_mco2( 1, :,10) = (/ &
     & 9.44086e-02_jprb, 1.02788e-01_jprb, 1.11911e-01_jprb, 1.21844e-01_jprb, 1.32659e-01_jprb, &
     & 1.44434e-01_jprb, 1.57253e-01_jprb, 1.71211e-01_jprb, 1.86407e-01_jprb, 2.02952e-01_jprb, &
     & 2.20966e-01_jprb, 2.40578e-01_jprb, 2.61932e-01_jprb, 2.85180e-01_jprb, 3.10492e-01_jprb, &
     & 3.38051e-01_jprb, 3.68056e-01_jprb, 4.00723e-01_jprb, 4.36291e-01_jprb/)
      kao_mco2( 2, :,10) = (/ &
     & 1.29528e-01_jprb, 1.39646e-01_jprb, 1.50554e-01_jprb, 1.62315e-01_jprb, 1.74994e-01_jprb, &
     & 1.88664e-01_jprb, 2.03401e-01_jprb, 2.19290e-01_jprb, 2.36419e-01_jprb, 2.54887e-01_jprb, &
     & 2.74798e-01_jprb, 2.96263e-01_jprb, 3.19406e-01_jprb, 3.44356e-01_jprb, 3.71255e-01_jprb, &
     & 4.00256e-01_jprb, 4.31522e-01_jprb, 4.65230e-01_jprb, 5.01571e-01_jprb/)
      kao_mco2( 3, :,10) = (/ &
     & 1.52325e-01_jprb, 1.62991e-01_jprb, 1.74404e-01_jprb, 1.86616e-01_jprb, 1.99684e-01_jprb, &
     & 2.13666e-01_jprb, 2.28628e-01_jprb, 2.44637e-01_jprb, 2.61767e-01_jprb, 2.80096e-01_jprb, &
     & 2.99710e-01_jprb, 3.20696e-01_jprb, 3.43152e-01_jprb, 3.67181e-01_jprb, 3.92892e-01_jprb, &
     & 4.20403e-01_jprb, 4.49841e-01_jprb, 4.81340e-01_jprb, 5.15045e-01_jprb/)
      kao_mco2( 4, :,10) = (/ &
     & 1.59763e-01_jprb, 1.70378e-01_jprb, 1.81698e-01_jprb, 1.93770e-01_jprb, 2.06644e-01_jprb, &
     & 2.20373e-01_jprb, 2.35015e-01_jprb, 2.50629e-01_jprb, 2.67281e-01_jprb, 2.85039e-01_jprb, &
     & 3.03977e-01_jprb, 3.24174e-01_jprb, 3.45712e-01_jprb, 3.68681e-01_jprb, 3.93176e-01_jprb, &
     & 4.19299e-01_jprb, 4.47157e-01_jprb, 4.76866e-01_jprb, 5.08549e-01_jprb/)
      kao_mco2( 5, :,10) = (/ &
     & 1.79202e-01_jprb, 1.91125e-01_jprb, 2.03840e-01_jprb, 2.17402e-01_jprb, 2.31866e-01_jprb, &
     & 2.47292e-01_jprb, 2.63744e-01_jprb, 2.81291e-01_jprb, 3.00005e-01_jprb, 3.19964e-01_jprb, &
     & 3.41251e-01_jprb, 3.63955e-01_jprb, 3.88169e-01_jprb, 4.13994e-01_jprb, 4.41537e-01_jprb, &
     & 4.70912e-01_jprb, 5.02242e-01_jprb, 5.35656e-01_jprb, 5.71293e-01_jprb/)
      kao_mco2( 6, :,10) = (/ &
     & 1.66628e-01_jprb, 1.76984e-01_jprb, 1.87984e-01_jprb, 1.99668e-01_jprb, 2.12078e-01_jprb, &
     & 2.25259e-01_jprb, 2.39259e-01_jprb, 2.54129e-01_jprb, 2.69924e-01_jprb, 2.86700e-01_jprb, &
     & 3.04519e-01_jprb, 3.23446e-01_jprb, 3.43549e-01_jprb, 3.64901e-01_jprb, 3.87580e-01_jprb, &
     & 4.11669e-01_jprb, 4.37255e-01_jprb, 4.64431e-01_jprb, 4.93297e-01_jprb/)
      kao_mco2( 7, :,10) = (/ &
     & 2.03980e-01_jprb, 2.17141e-01_jprb, 2.31152e-01_jprb, 2.46067e-01_jprb, 2.61945e-01_jprb, &
     & 2.78847e-01_jprb, 2.96839e-01_jprb, 3.15993e-01_jprb, 3.36382e-01_jprb, 3.58087e-01_jprb, &
     & 3.81193e-01_jprb, 4.05789e-01_jprb, 4.31972e-01_jprb, 4.59845e-01_jprb, 4.89517e-01_jprb, &
     & 5.21103e-01_jprb, 5.54727e-01_jprb, 5.90520e-01_jprb, 6.28623e-01_jprb/)
      kao_mco2( 8, :,10) = (/ &
     & 1.96161e-04_jprb, 2.07177e-04_jprb, 2.18812e-04_jprb, 2.31101e-04_jprb, 2.44079e-04_jprb, &
     & 2.57787e-04_jprb, 2.72264e-04_jprb, 2.87554e-04_jprb, 3.03703e-04_jprb, 3.20758e-04_jprb, &
     & 3.38772e-04_jprb, 3.57797e-04_jprb, 3.77891e-04_jprb, 3.99113e-04_jprb, 4.21527e-04_jprb, &
     & 4.45200e-04_jprb, 4.70202e-04_jprb, 4.96608e-04_jprb, 5.24498e-04_jprb/)
      kao_mco2( 9, :,10) = (/ &
     & 1.76275e-01_jprb, 1.88091e-01_jprb, 2.00699e-01_jprb, 2.14152e-01_jprb, 2.28507e-01_jprb, &
     & 2.43824e-01_jprb, 2.60168e-01_jprb, 2.77607e-01_jprb, 2.96216e-01_jprb, 3.16071e-01_jprb, &
     & 3.37258e-01_jprb, 3.59865e-01_jprb, 3.83987e-01_jprb, 4.09726e-01_jprb, 4.37190e-01_jprb, &
     & 4.66495e-01_jprb, 4.97765e-01_jprb, 5.31131e-01_jprb, 5.66733e-01_jprb/)
      kao_mco2( 1, :,11) = (/ &
     & 1.99797e-01_jprb, 2.14154e-01_jprb, 2.29543e-01_jprb, 2.46038e-01_jprb, 2.63718e-01_jprb, &
     & 2.82669e-01_jprb, 3.02981e-01_jprb, 3.24753e-01_jprb, 3.48090e-01_jprb, 3.73104e-01_jprb, &
     & 3.99915e-01_jprb, 4.28652e-01_jprb, 4.59455e-01_jprb, 4.92471e-01_jprb, 5.27859e-01_jprb, &
     & 5.65791e-01_jprb, 6.06448e-01_jprb, 6.50027e-01_jprb, 6.96738e-01_jprb/)
      kao_mco2( 2, :,11) = (/ &
     & 2.20638e-01_jprb, 2.35685e-01_jprb, 2.51759e-01_jprb, 2.68929e-01_jprb, 2.87271e-01_jprb, &
     & 3.06863e-01_jprb, 3.27791e-01_jprb, 3.50146e-01_jprb, 3.74026e-01_jprb, 3.99535e-01_jprb, &
     & 4.26784e-01_jprb, 4.55891e-01_jprb, 4.86983e-01_jprb, 5.20195e-01_jprb, 5.55673e-01_jprb, &
     & 5.93570e-01_jprb, 6.34052e-01_jprb, 6.77294e-01_jprb, 7.23486e-01_jprb/)
      kao_mco2( 3, :,11) = (/ &
     & 2.62988e-01_jprb, 2.80924e-01_jprb, 3.00085e-01_jprb, 3.20552e-01_jprb, 3.42414e-01_jprb, &
     & 3.65768e-01_jprb, 3.90715e-01_jprb, 4.17363e-01_jprb, 4.45829e-01_jprb, 4.76237e-01_jprb, &
     & 5.08718e-01_jprb, 5.43414e-01_jprb, 5.80477e-01_jprb, 6.20068e-01_jprb, 6.62359e-01_jprb, &
     & 7.07535e-01_jprb, 7.55791e-01_jprb, 8.07339e-01_jprb, 8.62403e-01_jprb/)
      kao_mco2( 4, :,11) = (/ &
     & 2.43674e-01_jprb, 2.59946e-01_jprb, 2.77304e-01_jprb, 2.95821e-01_jprb, 3.15575e-01_jprb, &
     & 3.36647e-01_jprb, 3.59127e-01_jprb, 3.83108e-01_jprb, 4.08691e-01_jprb, 4.35981e-01_jprb, &
     & 4.65094e-01_jprb, 4.96152e-01_jprb, 5.29282e-01_jprb, 5.64626e-01_jprb, 6.02329e-01_jprb, &
     & 6.42550e-01_jprb, 6.85457e-01_jprb, 7.31229e-01_jprb, 7.80057e-01_jprb/)
      kao_mco2( 5, :,11) = (/ &
     & 2.23323e-01_jprb, 2.37553e-01_jprb, 2.52689e-01_jprb, 2.68791e-01_jprb, 2.85918e-01_jprb, &
     & 3.04136e-01_jprb, 3.23515e-01_jprb, 3.44129e-01_jprb, 3.66057e-01_jprb, 3.89381e-01_jprb, &
     & 4.14192e-01_jprb, 4.40584e-01_jprb, 4.68657e-01_jprb, 4.98520e-01_jprb, 5.30285e-01_jprb, &
     & 5.64074e-01_jprb, 6.00016e-01_jprb, 6.38248e-01_jprb, 6.78917e-01_jprb/)
      kao_mco2( 6, :,11) = (/ &
     & 2.83716e-01_jprb, 3.02622e-01_jprb, 3.22788e-01_jprb, 3.44298e-01_jprb, 3.67241e-01_jprb, &
     & 3.91713e-01_jprb, 4.17816e-01_jprb, 4.45658e-01_jprb, 4.75356e-01_jprb, 5.07033e-01_jprb, &
     & 5.40820e-01_jprb, 5.76859e-01_jprb, 6.15300e-01_jprb, 6.56302e-01_jprb, 7.00037e-01_jprb, &
     & 7.46686e-01_jprb, 7.96443e-01_jprb, 8.49516e-01_jprb, 9.06126e-01_jprb/)
      kao_mco2( 7, :,11) = (/ &
     & 1.00497e-03_jprb, 1.06500e-03_jprb, 1.12863e-03_jprb, 1.19606e-03_jprb, 1.26751e-03_jprb, &
     & 1.34323e-03_jprb, 1.42348e-03_jprb, 1.50852e-03_jprb, 1.59864e-03_jprb, 1.69414e-03_jprb, &
     & 1.79535e-03_jprb, 1.90261e-03_jprb, 2.01628e-03_jprb, 2.13673e-03_jprb, 2.26438e-03_jprb, &
     & 2.39966e-03_jprb, 2.54302e-03_jprb, 2.69494e-03_jprb, 2.85594e-03_jprb/)
      kao_mco2( 8, :,11) = (/ &
     & 3.22623e-04_jprb, 3.39937e-04_jprb, 3.58181e-04_jprb, 3.77404e-04_jprb, 3.97658e-04_jprb, &
     & 4.19000e-04_jprb, 4.41487e-04_jprb, 4.65180e-04_jprb, 4.90146e-04_jprb, 5.16451e-04_jprb, &
     & 5.44167e-04_jprb, 5.73372e-04_jprb, 6.04143e-04_jprb, 6.36567e-04_jprb, 6.70730e-04_jprb, &
     & 7.06726e-04_jprb, 7.44655e-04_jprb, 7.84619e-04_jprb, 8.26727e-04_jprb/)
      kao_mco2( 9, :,11) = (/ &
     & 2.23872e-01_jprb, 2.38360e-01_jprb, 2.53786e-01_jprb, 2.70210e-01_jprb, 2.87697e-01_jprb, &
     & 3.06316e-01_jprb, 3.26140e-01_jprb, 3.47247e-01_jprb, 3.69720e-01_jprb, 3.93647e-01_jprb, &
     & 4.19122e-01_jprb, 4.46246e-01_jprb, 4.75126e-01_jprb, 5.05874e-01_jprb, 5.38613e-01_jprb, &
     & 5.73470e-01_jprb, 6.10583e-01_jprb, 6.50098e-01_jprb, 6.92170e-01_jprb/)
      kao_mco2( 1, :,12) = (/ &
     & 3.52418e-01_jprb, 3.76085e-01_jprb, 4.01341e-01_jprb, 4.28293e-01_jprb, 4.57055e-01_jprb, &
     & 4.87749e-01_jprb, 5.20504e-01_jprb, 5.55458e-01_jprb, 5.92760e-01_jprb, 6.32567e-01_jprb, &
     & 6.75047e-01_jprb, 7.20380e-01_jprb, 7.68757e-01_jprb, 8.20383e-01_jprb, 8.75476e-01_jprb, &
     & 9.34268e-01_jprb, 9.97009e-01_jprb, 1.06396e+00_jprb, 1.13541e+00_jprb/)
      kao_mco2( 2, :,12) = (/ &
     & 3.38812e-01_jprb, 3.61001e-01_jprb, 3.84645e-01_jprb, 4.09836e-01_jprb, 4.36678e-01_jprb, &
     & 4.65278e-01_jprb, 4.95750e-01_jprb, 5.28219e-01_jprb, 5.62814e-01_jprb, 5.99674e-01_jprb, &
     & 6.38949e-01_jprb, 6.80796e-01_jprb, 7.25384e-01_jprb, 7.72892e-01_jprb, 8.23511e-01_jprb, &
     & 8.77446e-01_jprb, 9.34913e-01_jprb, 9.96144e-01_jprb, 1.06138e+00_jprb/)
      kao_mco2( 3, :,12) = (/ &
     & 3.44644e-01_jprb, 3.66671e-01_jprb, 3.90105e-01_jprb, 4.15038e-01_jprb, 4.41564e-01_jprb, &
     & 4.69785e-01_jprb, 4.99810e-01_jprb, 5.31754e-01_jprb, 5.65740e-01_jprb, 6.01897e-01_jprb, &
     & 6.40366e-01_jprb, 6.81293e-01_jprb, 7.24836e-01_jprb, 7.71162e-01_jprb, 8.20448e-01_jprb, &
     & 8.72885e-01_jprb, 9.28673e-01_jprb, 9.88027e-01_jprb, 1.05117e+00_jprb/)
      kao_mco2( 4, :,12) = (/ &
     & 4.20358e-01_jprb, 4.47809e-01_jprb, 4.77053e-01_jprb, 5.08207e-01_jprb, 5.41395e-01_jprb, &
     & 5.76750e-01_jprb, 6.14414e-01_jprb, 6.54538e-01_jprb, 6.97282e-01_jprb, 7.42818e-01_jprb, &
     & 7.91327e-01_jprb, 8.43004e-01_jprb, 8.98056e-01_jprb, 9.56703e-01_jprb, 1.01918e+00_jprb, &
     & 1.08574e+00_jprb, 1.15664e+00_jprb, 1.23217e+00_jprb, 1.31264e+00_jprb/)
      kao_mco2( 5, :,12) = (/ &
     & 4.42756e-01_jprb, 4.72000e-01_jprb, 5.03174e-01_jprb, 5.36408e-01_jprb, 5.71837e-01_jprb, &
     & 6.09606e-01_jprb, 6.49870e-01_jprb, 6.92793e-01_jprb, 7.38551e-01_jprb, 7.87331e-01_jprb, &
     & 8.39333e-01_jprb, 8.94770e-01_jprb, 9.53868e-01_jprb, 1.01687e+00_jprb, 1.08403e+00_jprb, &
     & 1.15563e+00_jprb, 1.23196e+00_jprb, 1.31333e+00_jprb, 1.40007e+00_jprb/)
      kao_mco2( 6, :,12) = (/ &
     & 1.53662e-01_jprb, 1.63104e-01_jprb, 1.73126e-01_jprb, 1.83764e-01_jprb, 1.95055e-01_jprb, &
     & 2.07040e-01_jprb, 2.19762e-01_jprb, 2.33265e-01_jprb, 2.47598e-01_jprb, 2.62811e-01_jprb, &
     & 2.78960e-01_jprb, 2.96100e-01_jprb, 3.14294e-01_jprb, 3.33606e-01_jprb, 3.54104e-01_jprb, &
     & 3.75862e-01_jprb, 3.98956e-01_jprb, 4.23470e-01_jprb, 4.49490e-01_jprb/)
      kao_mco2( 7, :,12) = (/ &
     & 5.41472e-04_jprb, 5.65116e-04_jprb, 5.89793e-04_jprb, 6.15547e-04_jprb, 6.42426e-04_jprb, &
     & 6.70479e-04_jprb, 6.99757e-04_jprb, 7.30313e-04_jprb, 7.62203e-04_jprb, 7.95486e-04_jprb, &
     & 8.30223e-04_jprb, 8.66476e-04_jprb, 9.04312e-04_jprb, 9.43801e-04_jprb, 9.85013e-04_jprb, &
     & 1.02803e-03_jprb, 1.07292e-03_jprb, 1.11977e-03_jprb, 1.16866e-03_jprb/)
      kao_mco2( 8, :,12) = (/ &
     & 5.94251e-04_jprb, 6.17650e-04_jprb, 6.41969e-04_jprb, 6.67246e-04_jprb, 6.93518e-04_jprb, &
     & 7.20824e-04_jprb, 7.49206e-04_jprb, 7.78705e-04_jprb, 8.09366e-04_jprb, 8.41234e-04_jprb, &
     & 8.74356e-04_jprb, 9.08783e-04_jprb, 9.44566e-04_jprb, 9.81757e-04_jprb, 1.02041e-03_jprb, &
     & 1.06059e-03_jprb, 1.10235e-03_jprb, 1.14575e-03_jprb, 1.19087e-03_jprb/)
      kao_mco2( 9, :,12) = (/ &
     & 4.21683e-01_jprb, 4.49025e-01_jprb, 4.78140e-01_jprb, 5.09142e-01_jprb, 5.42155e-01_jprb, &
     & 5.77309e-01_jprb, 6.14741e-01_jprb, 6.54601e-01_jprb, 6.97046e-01_jprb, 7.42242e-01_jprb, &
     & 7.90369e-01_jprb, 8.41617e-01_jprb, 8.96188e-01_jprb, 9.54297e-01_jprb, 1.01617e+00_jprb, &
     & 1.08206e+00_jprb, 1.15222e+00_jprb, 1.22693e+00_jprb, 1.30649e+00_jprb/)
      kao_mco2( 1, :,13) = (/ &
     & 5.61805e-01_jprb, 5.98988e-01_jprb, 6.38631e-01_jprb, 6.80898e-01_jprb, 7.25962e-01_jprb, &
     & 7.74009e-01_jprb, 8.25236e-01_jprb, 8.79853e-01_jprb, 9.38085e-01_jprb, 1.00017e+00_jprb, &
     & 1.06637e+00_jprb, 1.13694e+00_jprb, 1.21219e+00_jprb, 1.29242e+00_jprb, 1.37795e+00_jprb, &
     & 1.46915e+00_jprb, 1.56638e+00_jprb, 1.67005e+00_jprb, 1.78058e+00_jprb/)
      kao_mco2( 2, :,13) = (/ &
     & 5.55938e-01_jprb, 5.91800e-01_jprb, 6.29976e-01_jprb, 6.70615e-01_jprb, 7.13876e-01_jprb, &
     & 7.59927e-01_jprb, 8.08949e-01_jprb, 8.61133e-01_jprb, 9.16683e-01_jprb, 9.75817e-01_jprb, &
     & 1.03877e+00_jprb, 1.10577e+00_jprb, 1.17711e+00_jprb, 1.25304e+00_jprb, 1.33387e+00_jprb, &
     & 1.41992e+00_jprb, 1.51152e+00_jprb, 1.60902e+00_jprb, 1.71282e+00_jprb/)
      kao_mco2( 3, :,13) = (/ &
     & 5.94615e-01_jprb, 6.33277e-01_jprb, 6.74453e-01_jprb, 7.18307e-01_jprb, 7.65012e-01_jprb, &
     & 8.14753e-01_jprb, 8.67729e-01_jprb, 9.24149e-01_jprb, 9.84238e-01_jprb, 1.04823e+00_jprb, &
     & 1.11639e+00_jprb, 1.18898e+00_jprb, 1.26629e+00_jprb, 1.34862e+00_jprb, 1.43631e+00_jprb, &
     & 1.52970e+00_jprb, 1.62916e+00_jprb, 1.73509e+00_jprb, 1.84791e+00_jprb/)
      kao_mco2( 4, :,13) = (/ &
     & 5.48973e-01_jprb, 5.84145e-01_jprb, 6.21570e-01_jprb, 6.61394e-01_jprb, 7.03768e-01_jprb, &
     & 7.48858e-01_jprb, 7.96836e-01_jprb, 8.47889e-01_jprb, 9.02212e-01_jprb, 9.60015e-01_jprb, &
     & 1.02152e+00_jprb, 1.08697e+00_jprb, 1.15661e+00_jprb, 1.23071e+00_jprb, 1.30956e+00_jprb, &
     & 1.39347e+00_jprb, 1.48274e+00_jprb, 1.57774e+00_jprb, 1.67883e+00_jprb/)
      kao_mco2( 5, :,13) = (/ &
     & 1.49742e-01_jprb, 1.59049e-01_jprb, 1.68934e-01_jprb, 1.79434e-01_jprb, 1.90586e-01_jprb, &
     & 2.02432e-01_jprb, 2.15013e-01_jprb, 2.28377e-01_jprb, 2.42571e-01_jprb, 2.57648e-01_jprb, &
     & 2.73661e-01_jprb, 2.90670e-01_jprb, 3.08736e-01_jprb, 3.27925e-01_jprb, 3.48307e-01_jprb, &
     & 3.69955e-01_jprb, 3.92949e-01_jprb, 4.17372e-01_jprb, 4.43312e-01_jprb/)
      kao_mco2( 6, :,13) = (/ &
     & 8.81777e-04_jprb, 9.16690e-04_jprb, 9.52985e-04_jprb, 9.90718e-04_jprb, 1.02994e-03_jprb, &
     & 1.07072e-03_jprb, 1.11312e-03_jprb, 1.15719e-03_jprb, 1.20301e-03_jprb, 1.25064e-03_jprb, &
     & 1.30016e-03_jprb, 1.35163e-03_jprb, 1.40515e-03_jprb, 1.46079e-03_jprb, 1.51862e-03_jprb, &
     & 1.57875e-03_jprb, 1.64126e-03_jprb, 1.70624e-03_jprb, 1.77380e-03_jprb/)
      kao_mco2( 7, :,13) = (/ &
     & 8.84366e-04_jprb, 9.20446e-04_jprb, 9.57999e-04_jprb, 9.97083e-04_jprb, 1.03776e-03_jprb, &
     & 1.08010e-03_jprb, 1.12417e-03_jprb, 1.17003e-03_jprb, 1.21777e-03_jprb, 1.26745e-03_jprb, &
     & 1.31916e-03_jprb, 1.37298e-03_jprb, 1.42899e-03_jprb, 1.48729e-03_jprb, 1.54797e-03_jprb, &
     & 1.61113e-03_jprb, 1.67686e-03_jprb, 1.74527e-03_jprb, 1.81647e-03_jprb/)
      kao_mco2( 8, :,13) = (/ &
     & 8.92597e-04_jprb, 9.33069e-04_jprb, 9.75377e-04_jprb, 1.01960e-03_jprb, 1.06583e-03_jprb, &
     & 1.11416e-03_jprb, 1.16468e-03_jprb, 1.21749e-03_jprb, 1.27269e-03_jprb, 1.33040e-03_jprb, &
     & 1.39073e-03_jprb, 1.45378e-03_jprb, 1.51970e-03_jprb, 1.58861e-03_jprb, 1.66064e-03_jprb, &
     & 1.73594e-03_jprb, 1.81465e-03_jprb, 1.89693e-03_jprb, 1.98294e-03_jprb/)
      kao_mco2( 9, :,13) = (/ &
     & 1.46280e-01_jprb, 1.55378e-01_jprb, 1.65043e-01_jprb, 1.75308e-01_jprb, 1.86212e-01_jprb, &
     & 1.97794e-01_jprb, 2.10097e-01_jprb, 2.23164e-01_jprb, 2.37045e-01_jprb, 2.51788e-01_jprb, &
     & 2.67449e-01_jprb, 2.84084e-01_jprb, 3.01754e-01_jprb, 3.20522e-01_jprb, 3.40458e-01_jprb, &
     & 3.61634e-01_jprb, 3.84127e-01_jprb, 4.08020e-01_jprb, 4.33398e-01_jprb/)
      kao_mco2( 1, :,14) = (/ &
     & 9.20236e-01_jprb, 9.80010e-01_jprb, 1.04367e+00_jprb, 1.11146e+00_jprb, 1.18366e+00_jprb, &
     & 1.26054e+00_jprb, 1.34242e+00_jprb, 1.42962e+00_jprb, 1.52248e+00_jprb, 1.62137e+00_jprb, &
     & 1.72669e+00_jprb, 1.83885e+00_jprb, 1.95829e+00_jprb, 2.08549e+00_jprb, 2.22096e+00_jprb, &
     & 2.36522e+00_jprb, 2.51886e+00_jprb, 2.68247e+00_jprb, 2.85671e+00_jprb/)
      kao_mco2( 2, :,14) = (/ &
     & 8.39823e-01_jprb, 8.95471e-01_jprb, 9.54806e-01_jprb, 1.01807e+00_jprb, 1.08553e+00_jprb, &
     & 1.15746e+00_jprb, 1.23416e+00_jprb, 1.31593e+00_jprb, 1.40313e+00_jprb, 1.49610e+00_jprb, &
     & 1.59523e+00_jprb, 1.70094e+00_jprb, 1.81364e+00_jprb, 1.93382e+00_jprb, 2.06195e+00_jprb, &
     & 2.19858e+00_jprb, 2.34426e+00_jprb, 2.49960e+00_jprb, 2.66522e+00_jprb/)
      kao_mco2( 3, :,14) = (/ &
     & 5.39252e-01_jprb, 5.73971e-01_jprb, 6.10925e-01_jprb, 6.50259e-01_jprb, 6.92125e-01_jprb, &
     & 7.36686e-01_jprb, 7.84117e-01_jprb, 8.34601e-01_jprb, 8.88336e-01_jprb, 9.45530e-01_jprb, &
     & 1.00641e+00_jprb, 1.07120e+00_jprb, 1.14017e+00_jprb, 1.21358e+00_jprb, 1.29171e+00_jprb, &
     & 1.37488e+00_jprb, 1.46340e+00_jprb, 1.55762e+00_jprb, 1.65790e+00_jprb/)
      kao_mco2( 4, :,14) = (/ &
     & 1.14837e-03_jprb, 1.19701e-03_jprb, 1.24770e-03_jprb, 1.30055e-03_jprb, 1.35563e-03_jprb, &
     & 1.41305e-03_jprb, 1.47289e-03_jprb, 1.53528e-03_jprb, 1.60030e-03_jprb, 1.66808e-03_jprb, &
     & 1.73873e-03_jprb, 1.81237e-03_jprb, 1.88913e-03_jprb, 1.96914e-03_jprb, 2.05254e-03_jprb, &
     & 2.13947e-03_jprb, 2.23009e-03_jprb, 2.32454e-03_jprb, 2.42299e-03_jprb/)
      kao_mco2( 5, :,14) = (/ &
     & 1.14611e-03_jprb, 1.19424e-03_jprb, 1.24440e-03_jprb, 1.29666e-03_jprb, 1.35111e-03_jprb, &
     & 1.40786e-03_jprb, 1.46698e-03_jprb, 1.52859e-03_jprb, 1.59279e-03_jprb, 1.65968e-03_jprb, &
     & 1.72938e-03_jprb, 1.80201e-03_jprb, 1.87769e-03_jprb, 1.95655e-03_jprb, 2.03872e-03_jprb, &
     & 2.12434e-03_jprb, 2.21355e-03_jprb, 2.30651e-03_jprb, 2.40338e-03_jprb/)
      kao_mco2( 6, :,14) = (/ &
     & 1.14203e-03_jprb, 1.18930e-03_jprb, 1.23852e-03_jprb, 1.28979e-03_jprb, 1.34317e-03_jprb, &
     & 1.39877e-03_jprb, 1.45666e-03_jprb, 1.51695e-03_jprb, 1.57974e-03_jprb, 1.64513e-03_jprb, &
     & 1.71322e-03_jprb, 1.78413e-03_jprb, 1.85798e-03_jprb, 1.93488e-03_jprb, 2.01497e-03_jprb, &
     & 2.09837e-03_jprb, 2.18522e-03_jprb, 2.27567e-03_jprb, 2.36986e-03_jprb/)
      kao_mco2( 7, :,14) = (/ &
     & 1.11217e-03_jprb, 1.15727e-03_jprb, 1.20421e-03_jprb, 1.25305e-03_jprb, 1.30386e-03_jprb, &
     & 1.35674e-03_jprb, 1.41177e-03_jprb, 1.46902e-03_jprb, 1.52860e-03_jprb, 1.59059e-03_jprb, &
     & 1.65510e-03_jprb, 1.72222e-03_jprb, 1.79207e-03_jprb, 1.86475e-03_jprb, 1.94037e-03_jprb, &
     & 2.01907e-03_jprb, 2.10095e-03_jprb, 2.18616e-03_jprb, 2.27482e-03_jprb/)
      kao_mco2( 8, :,14) = (/ &
     & 1.21596e-03_jprb, 1.25817e-03_jprb, 1.30183e-03_jprb, 1.34702e-03_jprb, 1.39377e-03_jprb, &
     & 1.44214e-03_jprb, 1.49219e-03_jprb, 1.54398e-03_jprb, 1.59757e-03_jprb, 1.65302e-03_jprb, &
     & 1.71039e-03_jprb, 1.76975e-03_jprb, 1.83117e-03_jprb, 1.89473e-03_jprb, 1.96049e-03_jprb, &
     & 2.02853e-03_jprb, 2.09893e-03_jprb, 2.17178e-03_jprb, 2.24716e-03_jprb/)
      kao_mco2( 9, :,14) = (/ &
     & 1.14611e-03_jprb, 1.19424e-03_jprb, 1.24440e-03_jprb, 1.29666e-03_jprb, 1.35111e-03_jprb, &
     & 1.40786e-03_jprb, 1.46698e-03_jprb, 1.52859e-03_jprb, 1.59279e-03_jprb, 1.65968e-03_jprb, &
     & 1.72938e-03_jprb, 1.80201e-03_jprb, 1.87769e-03_jprb, 1.95655e-03_jprb, 2.03872e-03_jprb, &
     & 2.12434e-03_jprb, 2.21355e-03_jprb, 2.30651e-03_jprb, 2.40338e-03_jprb/)
      kao_mco2( 1, :,15) = (/ &
     & 1.29470e+00_jprb, 1.37848e+00_jprb, 1.46768e+00_jprb, 1.56266e+00_jprb, 1.66378e+00_jprb, &
     & 1.77145e+00_jprb, 1.88609e+00_jprb, 2.00814e+00_jprb, 2.13809e+00_jprb, 2.27645e+00_jprb, &
     & 2.42376e+00_jprb, 2.58061e+00_jprb, 2.74761e+00_jprb, 2.92541e+00_jprb, 3.11472e+00_jprb, &
     & 3.31628e+00_jprb, 3.53088e+00_jprb, 3.75938e+00_jprb, 4.00265e+00_jprb/)
      kao_mco2( 2, :,15) = (/ &
     & 7.23701e-01_jprb, 7.68508e-01_jprb, 8.16089e-01_jprb, 8.66616e-01_jprb, 9.20272e-01_jprb, &
     & 9.77250e-01_jprb, 1.03775e+00_jprb, 1.10201e+00_jprb, 1.17024e+00_jprb, 1.24269e+00_jprb, &
     & 1.31963e+00_jprb, 1.40133e+00_jprb, 1.48809e+00_jprb, 1.58023e+00_jprb, 1.67807e+00_jprb, &
     & 1.78196e+00_jprb, 1.89229e+00_jprb, 2.00945e+00_jprb, 2.13386e+00_jprb/)
      kao_mco2( 3, :,15) = (/ &
     & 1.81684e-03_jprb, 1.85424e-03_jprb, 1.89241e-03_jprb, 1.93137e-03_jprb, 1.97114e-03_jprb, &
     & 2.01172e-03_jprb, 2.05313e-03_jprb, 2.09540e-03_jprb, 2.13854e-03_jprb, 2.18257e-03_jprb, &
     & 2.22750e-03_jprb, 2.27336e-03_jprb, 2.32016e-03_jprb, 2.36793e-03_jprb, 2.41668e-03_jprb, &
     & 2.46643e-03_jprb, 2.51721e-03_jprb, 2.56903e-03_jprb, 2.62192e-03_jprb/)
      kao_mco2( 4, :,15) = (/ &
     & 1.84644e-03_jprb, 1.88437e-03_jprb, 1.92309e-03_jprb, 1.96260e-03_jprb, 2.00293e-03_jprb, &
     & 2.04408e-03_jprb, 2.08608e-03_jprb, 2.12894e-03_jprb, 2.17268e-03_jprb, 2.21732e-03_jprb, &
     & 2.26288e-03_jprb, 2.30938e-03_jprb, 2.35683e-03_jprb, 2.40525e-03_jprb, 2.45467e-03_jprb, &
     & 2.50510e-03_jprb, 2.55658e-03_jprb, 2.60910e-03_jprb, 2.66271e-03_jprb/)
      kao_mco2( 5, :,15) = (/ &
     & 1.88579e-03_jprb, 1.92454e-03_jprb, 1.96408e-03_jprb, 2.00443e-03_jprb, 2.04561e-03_jprb, &
     & 2.08764e-03_jprb, 2.13054e-03_jprb, 2.17431e-03_jprb, 2.21898e-03_jprb, 2.26457e-03_jprb, &
     & 2.31110e-03_jprb, 2.35858e-03_jprb, 2.40704e-03_jprb, 2.45650e-03_jprb, 2.50697e-03_jprb, &
     & 2.55848e-03_jprb, 2.61104e-03_jprb, 2.66469e-03_jprb, 2.71943e-03_jprb/)
      kao_mco2( 6, :,15) = (/ &
     & 1.95322e-03_jprb, 1.99316e-03_jprb, 2.03391e-03_jprb, 2.07549e-03_jprb, 2.11793e-03_jprb, &
     & 2.16123e-03_jprb, 2.20542e-03_jprb, 2.25051e-03_jprb, 2.29652e-03_jprb, 2.34347e-03_jprb, &
     & 2.39139e-03_jprb, 2.44028e-03_jprb, 2.49017e-03_jprb, 2.54109e-03_jprb, 2.59304e-03_jprb, &
     & 2.64605e-03_jprb, 2.70015e-03_jprb, 2.75536e-03_jprb, 2.81169e-03_jprb/)
      kao_mco2( 7, :,15) = (/ &
     & 2.13640e-03_jprb, 2.17976e-03_jprb, 2.22400e-03_jprb, 2.26914e-03_jprb, 2.31520e-03_jprb, &
     & 2.36219e-03_jprb, 2.41013e-03_jprb, 2.45905e-03_jprb, 2.50896e-03_jprb, 2.55988e-03_jprb, &
     & 2.61184e-03_jprb, 2.66485e-03_jprb, 2.71893e-03_jprb, 2.77412e-03_jprb, 2.83042e-03_jprb, &
     & 2.88787e-03_jprb, 2.94648e-03_jprb, 3.00629e-03_jprb, 3.06730e-03_jprb/)
      kao_mco2( 8, :,15) = (/ &
     & 2.17014e-03_jprb, 2.21411e-03_jprb, 2.25897e-03_jprb, 2.30474e-03_jprb, 2.35143e-03_jprb, &
     & 2.39908e-03_jprb, 2.44769e-03_jprb, 2.49728e-03_jprb, 2.54788e-03_jprb, 2.59950e-03_jprb, &
     & 2.65217e-03_jprb, 2.70591e-03_jprb, 2.76073e-03_jprb, 2.81667e-03_jprb, 2.87374e-03_jprb, &
     & 2.93197e-03_jprb, 2.99137e-03_jprb, 3.05198e-03_jprb, 3.11382e-03_jprb/)
      kao_mco2( 9, :,15) = (/ &
     & 1.88579e-03_jprb, 1.92454e-03_jprb, 1.96408e-03_jprb, 2.00443e-03_jprb, 2.04561e-03_jprb, &
     & 2.08764e-03_jprb, 2.13054e-03_jprb, 2.17431e-03_jprb, 2.21898e-03_jprb, 2.26457e-03_jprb, &
     & 2.31110e-03_jprb, 2.35858e-03_jprb, 2.40704e-03_jprb, 2.45650e-03_jprb, 2.50697e-03_jprb, &
     & 2.55848e-03_jprb, 2.61104e-03_jprb, 2.66469e-03_jprb, 2.71943e-03_jprb/)
      kao_mco2( 1, :,16) = (/ &
     & 1.48989e+00_jprb, 1.58377e+00_jprb, 1.68356e+00_jprb, 1.78964e+00_jprb, 1.90241e+00_jprb, &
     & 2.02228e+00_jprb, 2.14971e+00_jprb, 2.28516e+00_jprb, 2.42915e+00_jprb, 2.58221e+00_jprb, &
     & 2.74492e+00_jprb, 2.91788e+00_jprb, 3.10174e+00_jprb, 3.29718e+00_jprb, 3.50494e+00_jprb, &
     & 3.72578e+00_jprb, 3.96055e+00_jprb, 4.21010e+00_jprb, 4.47538e+00_jprb/)
      kao_mco2( 2, :,16) = (/ &
     & 2.10609e-03_jprb, 2.14759e-03_jprb, 2.18992e-03_jprb, 2.23307e-03_jprb, 2.27708e-03_jprb, &
     & 2.32196e-03_jprb, 2.36771e-03_jprb, 2.41438e-03_jprb, 2.46196e-03_jprb, 2.51047e-03_jprb, &
     & 2.55995e-03_jprb, 2.61040e-03_jprb, 2.66184e-03_jprb, 2.71430e-03_jprb, 2.76779e-03_jprb, &
     & 2.82234e-03_jprb, 2.87796e-03_jprb, 2.93467e-03_jprb, 2.99251e-03_jprb/)
      kao_mco2( 3, :,16) = (/ &
     & 2.10609e-03_jprb, 2.14759e-03_jprb, 2.18992e-03_jprb, 2.23307e-03_jprb, 2.27708e-03_jprb, &
     & 2.32196e-03_jprb, 2.36771e-03_jprb, 2.41438e-03_jprb, 2.46196e-03_jprb, 2.51047e-03_jprb, &
     & 2.55995e-03_jprb, 2.61040e-03_jprb, 2.66184e-03_jprb, 2.71430e-03_jprb, 2.76779e-03_jprb, &
     & 2.82234e-03_jprb, 2.87796e-03_jprb, 2.93467e-03_jprb, 2.99251e-03_jprb/)
      kao_mco2( 4, :,16) = (/ &
     & 2.10609e-03_jprb, 2.14759e-03_jprb, 2.18992e-03_jprb, 2.23307e-03_jprb, 2.27708e-03_jprb, &
     & 2.32196e-03_jprb, 2.36771e-03_jprb, 2.41438e-03_jprb, 2.46196e-03_jprb, 2.51047e-03_jprb, &
     & 2.55995e-03_jprb, 2.61040e-03_jprb, 2.66184e-03_jprb, 2.71430e-03_jprb, 2.76779e-03_jprb, &
     & 2.82234e-03_jprb, 2.87796e-03_jprb, 2.93467e-03_jprb, 2.99251e-03_jprb/)
      kao_mco2( 5, :,16) = (/ &
     & 2.10609e-03_jprb, 2.14759e-03_jprb, 2.18992e-03_jprb, 2.23307e-03_jprb, 2.27708e-03_jprb, &
     & 2.32196e-03_jprb, 2.36771e-03_jprb, 2.41438e-03_jprb, 2.46196e-03_jprb, 2.51047e-03_jprb, &
     & 2.55995e-03_jprb, 2.61040e-03_jprb, 2.66184e-03_jprb, 2.71430e-03_jprb, 2.76779e-03_jprb, &
     & 2.82234e-03_jprb, 2.87796e-03_jprb, 2.93467e-03_jprb, 2.99251e-03_jprb/)
      kao_mco2( 6, :,16) = (/ &
     & 2.10609e-03_jprb, 2.14759e-03_jprb, 2.18992e-03_jprb, 2.23307e-03_jprb, 2.27708e-03_jprb, &
     & 2.32196e-03_jprb, 2.36771e-03_jprb, 2.41438e-03_jprb, 2.46196e-03_jprb, 2.51047e-03_jprb, &
     & 2.55995e-03_jprb, 2.61040e-03_jprb, 2.66184e-03_jprb, 2.71430e-03_jprb, 2.76779e-03_jprb, &
     & 2.82234e-03_jprb, 2.87796e-03_jprb, 2.93467e-03_jprb, 2.99251e-03_jprb/)
      kao_mco2( 7, :,16) = (/ &
     & 2.09970e-03_jprb, 2.14101e-03_jprb, 2.18312e-03_jprb, 2.22606e-03_jprb, 2.26985e-03_jprb, &
     & 2.31450e-03_jprb, 2.36003e-03_jprb, 2.40645e-03_jprb, 2.45379e-03_jprb, 2.50205e-03_jprb, &
     & 2.55127e-03_jprb, 2.60146e-03_jprb, 2.65263e-03_jprb, 2.70481e-03_jprb, 2.75801e-03_jprb, &
     & 2.81226e-03_jprb, 2.86758e-03_jprb, 2.92399e-03_jprb, 2.98150e-03_jprb/)
      kao_mco2( 8, :,16) = (/ &
     & 2.09970e-03_jprb, 2.14101e-03_jprb, 2.18312e-03_jprb, 2.22606e-03_jprb, 2.26985e-03_jprb, &
     & 2.31450e-03_jprb, 2.36003e-03_jprb, 2.40645e-03_jprb, 2.45379e-03_jprb, 2.50205e-03_jprb, &
     & 2.55127e-03_jprb, 2.60146e-03_jprb, 2.65263e-03_jprb, 2.70481e-03_jprb, 2.75801e-03_jprb, &
     & 2.81226e-03_jprb, 2.86758e-03_jprb, 2.92399e-03_jprb, 2.98150e-03_jprb/)
      kao_mco2( 9, :,16) = (/ &
     & 2.10609e-03_jprb, 2.14759e-03_jprb, 2.18992e-03_jprb, 2.23307e-03_jprb, 2.27708e-03_jprb, &
     & 2.32196e-03_jprb, 2.36771e-03_jprb, 2.41438e-03_jprb, 2.46196e-03_jprb, 2.51047e-03_jprb, &
     & 2.55995e-03_jprb, 2.61040e-03_jprb, 2.66184e-03_jprb, 2.71430e-03_jprb, 2.76779e-03_jprb, &
     & 2.82234e-03_jprb, 2.87796e-03_jprb, 2.93467e-03_jprb, 2.99251e-03_jprb/)

      kao_mco( 1, :, 1) = (/ &
     & 4.58355e-01_jprb, 4.47074e-01_jprb, 4.36070e-01_jprb, 4.25337e-01_jprb, 4.14868e-01_jprb, &
     & 4.04657e-01_jprb, 3.94697e-01_jprb, 3.84982e-01_jprb, 3.75506e-01_jprb, 3.66264e-01_jprb, &
     & 3.57249e-01_jprb, 3.48456e-01_jprb, 3.39879e-01_jprb, 3.31514e-01_jprb, 3.23354e-01_jprb, &
     & 3.15395e-01_jprb, 3.07632e-01_jprb, 3.00061e-01_jprb, 2.92675e-01_jprb/)
      kao_mco( 2, :, 1) = (/ &
     & 7.03080e-01_jprb, 6.84132e-01_jprb, 6.65696e-01_jprb, 6.47756e-01_jprb, 6.30300e-01_jprb, &
     & 6.13314e-01_jprb, 5.96786e-01_jprb, 5.80703e-01_jprb, 5.65053e-01_jprb, 5.49826e-01_jprb, &
     & 5.35009e-01_jprb, 5.20591e-01_jprb, 5.06561e-01_jprb, 4.92910e-01_jprb, 4.79627e-01_jprb, &
     & 4.66701e-01_jprb, 4.54124e-01_jprb, 4.41886e-01_jprb, 4.29978e-01_jprb/)
      kao_mco( 3, :, 1) = (/ &
     & 8.53018e-01_jprb, 8.29537e-01_jprb, 8.06703e-01_jprb, 7.84497e-01_jprb, 7.62903e-01_jprb, &
     & 7.41903e-01_jprb, 7.21481e-01_jprb, 7.01621e-01_jprb, 6.82307e-01_jprb, 6.63526e-01_jprb, &
     & 6.45261e-01_jprb, 6.27499e-01_jprb, 6.10226e-01_jprb, 5.93429e-01_jprb, 5.77094e-01_jprb, &
     & 5.61208e-01_jprb, 5.45760e-01_jprb, 5.30737e-01_jprb, 5.16128e-01_jprb/)
      kao_mco( 4, :, 1) = (/ &
     & 9.58866e-01_jprb, 9.31881e-01_jprb, 9.05654e-01_jprb, 8.80166e-01_jprb, 8.55395e-01_jprb, &
     & 8.31321e-01_jprb, 8.07925e-01_jprb, 7.85187e-01_jprb, 7.63089e-01_jprb, 7.41613e-01_jprb, &
     & 7.20742e-01_jprb, 7.00457e-01_jprb, 6.80744e-01_jprb, 6.61586e-01_jprb, 6.42966e-01_jprb, &
     & 6.24871e-01_jprb, 6.07285e-01_jprb, 5.90194e-01_jprb, 5.73584e-01_jprb/)
      kao_mco( 5, :, 1) = (/ &
     & 1.07140e+00_jprb, 1.04056e+00_jprb, 1.01061e+00_jprb, 9.81521e-01_jprb, 9.53269e-01_jprb, &
     & 9.25829e-01_jprb, 8.99180e-01_jprb, 8.73297e-01_jprb, 8.48160e-01_jprb, 8.23746e-01_jprb, &
     & 8.00035e-01_jprb, 7.77006e-01_jprb, 7.54641e-01_jprb, 7.32919e-01_jprb, 7.11822e-01_jprb, &
     & 6.91333e-01_jprb, 6.71433e-01_jprb, 6.52106e-01_jprb, 6.33336e-01_jprb/)
      kao_mco( 6, :, 1) = (/ &
     & 1.21046e+00_jprb, 1.17478e+00_jprb, 1.14015e+00_jprb, 1.10655e+00_jprb, 1.07393e+00_jprb, &
     & 1.04228e+00_jprb, 1.01156e+00_jprb, 9.81740e-01_jprb, 9.52803e-01_jprb, 9.24720e-01_jprb, &
     & 8.97463e-01_jprb, 8.71011e-01_jprb, 8.45338e-01_jprb, 8.20422e-01_jprb, 7.96240e-01_jprb, &
     & 7.72771e-01_jprb, 7.49993e-01_jprb, 7.27887e-01_jprb, 7.06433e-01_jprb/)
      kao_mco( 7, :, 1) = (/ &
     & 1.57730e+00_jprb, 1.52919e+00_jprb, 1.48255e+00_jprb, 1.43733e+00_jprb, 1.39349e+00_jprb, &
     & 1.35099e+00_jprb, 1.30978e+00_jprb, 1.26983e+00_jprb, 1.23110e+00_jprb, 1.19355e+00_jprb, &
     & 1.15715e+00_jprb, 1.12186e+00_jprb, 1.08764e+00_jprb, 1.05446e+00_jprb, 1.02230e+00_jprb, &
     & 9.91121e-01_jprb, 9.60890e-01_jprb, 9.31583e-01_jprb, 9.03169e-01_jprb/)
      kao_mco( 8, :, 1) = (/ &
     & 2.43678e+00_jprb, 2.36595e+00_jprb, 2.29719e+00_jprb, 2.23042e+00_jprb, 2.16560e+00_jprb, &
     & 2.10266e+00_jprb, 2.04154e+00_jprb, 1.98221e+00_jprb, 1.92460e+00_jprb, 1.86866e+00_jprb, &
     & 1.81435e+00_jprb, 1.76162e+00_jprb, 1.71042e+00_jprb, 1.66070e+00_jprb, 1.61244e+00_jprb, &
     & 1.56557e+00_jprb, 1.52007e+00_jprb, 1.47589e+00_jprb, 1.43300e+00_jprb/)
      kao_mco( 9, :, 1) = (/ &
     & 9.66296e-01_jprb, 9.39903e-01_jprb, 9.14232e-01_jprb, 8.89262e-01_jprb, 8.64973e-01_jprb, &
     & 8.41348e-01_jprb, 8.18369e-01_jprb, 7.96017e-01_jprb, 7.74275e-01_jprb, 7.53128e-01_jprb, &
     & 7.32558e-01_jprb, 7.12549e-01_jprb, 6.93088e-01_jprb, 6.74157e-01_jprb, 6.55744e-01_jprb, &
     & 6.37834e-01_jprb, 6.20413e-01_jprb, 6.03468e-01_jprb, 5.86985e-01_jprb/)
      kao_mco( 1, :, 2) = (/ &
     & 1.15047e+00_jprb, 1.12127e+00_jprb, 1.09281e+00_jprb, 1.06507e+00_jprb, 1.03804e+00_jprb, &
     & 1.01169e+00_jprb, 9.86010e-01_jprb, 9.60983e-01_jprb, 9.36591e-01_jprb, 9.12818e-01_jprb, &
     & 8.89649e-01_jprb, 8.67067e-01_jprb, 8.45059e-01_jprb, 8.23610e-01_jprb, 8.02705e-01_jprb, &
     & 7.82330e-01_jprb, 7.62473e-01_jprb, 7.43119e-01_jprb, 7.24257e-01_jprb/)
      kao_mco( 2, :, 2) = (/ &
     & 1.43243e+00_jprb, 1.39430e+00_jprb, 1.35719e+00_jprb, 1.32106e+00_jprb, 1.28590e+00_jprb, &
     & 1.25167e+00_jprb, 1.21836e+00_jprb, 1.18593e+00_jprb, 1.15436e+00_jprb, 1.12364e+00_jprb, &
     & 1.09373e+00_jprb, 1.06462e+00_jprb, 1.03628e+00_jprb, 1.00870e+00_jprb, 9.81848e-01_jprb, &
     & 9.55714e-01_jprb, 9.30275e-01_jprb, 9.05514e-01_jprb, 8.81412e-01_jprb/)
      kao_mco( 3, :, 2) = (/ &
     & 1.61389e+00_jprb, 1.56911e+00_jprb, 1.52556e+00_jprb, 1.48323e+00_jprb, 1.44207e+00_jprb, &
     & 1.40205e+00_jprb, 1.36314e+00_jprb, 1.32531e+00_jprb, 1.28854e+00_jprb, 1.25278e+00_jprb, &
     & 1.21801e+00_jprb, 1.18421e+00_jprb, 1.15135e+00_jprb, 1.11940e+00_jprb, 1.08834e+00_jprb, &
     & 1.05814e+00_jprb, 1.02877e+00_jprb, 1.00022e+00_jprb, 9.72466e-01_jprb/)
      kao_mco( 4, :, 2) = (/ &
     & 1.78458e+00_jprb, 1.73440e+00_jprb, 1.68564e+00_jprb, 1.63825e+00_jprb, 1.59219e+00_jprb, &
     & 1.54742e+00_jprb, 1.50391e+00_jprb, 1.46163e+00_jprb, 1.42053e+00_jprb, 1.38059e+00_jprb, &
     & 1.34178e+00_jprb, 1.30405e+00_jprb, 1.26739e+00_jprb, 1.23175e+00_jprb, 1.19712e+00_jprb, &
     & 1.16346e+00_jprb, 1.13075e+00_jprb, 1.09896e+00_jprb, 1.06806e+00_jprb/)
      kao_mco( 5, :, 2) = (/ &
     & 1.92622e+00_jprb, 1.87172e+00_jprb, 1.81876e+00_jprb, 1.76730e+00_jprb, 1.71730e+00_jprb, &
     & 1.66871e+00_jprb, 1.62150e+00_jprb, 1.57562e+00_jprb, 1.53104e+00_jprb, 1.48772e+00_jprb, &
     & 1.44563e+00_jprb, 1.40473e+00_jprb, 1.36498e+00_jprb, 1.32636e+00_jprb, 1.28883e+00_jprb, &
     & 1.25237e+00_jprb, 1.21693e+00_jprb, 1.18250e+00_jprb, 1.14905e+00_jprb/)
      kao_mco( 6, :, 2) = (/ &
     & 2.23194e+00_jprb, 2.16782e+00_jprb, 2.10554e+00_jprb, 2.04505e+00_jprb, 1.98630e+00_jprb, &
     & 1.92924e+00_jprb, 1.87381e+00_jprb, 1.81998e+00_jprb, 1.76770e+00_jprb, 1.71691e+00_jprb, &
     & 1.66759e+00_jprb, 1.61968e+00_jprb, 1.57315e+00_jprb, 1.52796e+00_jprb, 1.48406e+00_jprb, &
     & 1.44143e+00_jprb, 1.40002e+00_jprb, 1.35980e+00_jprb, 1.32073e+00_jprb/)
      kao_mco( 7, :, 2) = (/ &
     & 2.64692e+00_jprb, 2.57290e+00_jprb, 2.50096e+00_jprb, 2.43103e+00_jprb, 2.36305e+00_jprb, &
     & 2.29697e+00_jprb, 2.23275e+00_jprb, 2.17031e+00_jprb, 2.10963e+00_jprb, 2.05064e+00_jprb, &
     & 1.99330e+00_jprb, 1.93756e+00_jprb, 1.88338e+00_jprb, 1.83072e+00_jprb, 1.77953e+00_jprb, &
     & 1.72977e+00_jprb, 1.68140e+00_jprb, 1.63438e+00_jprb, 1.58868e+00_jprb/)
      kao_mco( 8, :, 2) = (/ &
     & 2.86812e+00_jprb, 2.80121e+00_jprb, 2.73586e+00_jprb, 2.67204e+00_jprb, 2.60970e+00_jprb, &
     & 2.54882e+00_jprb, 2.48936e+00_jprb, 2.43129e+00_jprb, 2.37457e+00_jprb, 2.31917e+00_jprb, &
     & 2.26507e+00_jprb, 2.21223e+00_jprb, 2.16062e+00_jprb, 2.11022e+00_jprb, 2.06099e+00_jprb, &
     & 2.01291e+00_jprb, 1.96595e+00_jprb, 1.92009e+00_jprb, 1.87529e+00_jprb/)
      kao_mco( 9, :, 2) = (/ &
     & 1.25243e+00_jprb, 1.22790e+00_jprb, 1.20385e+00_jprb, 1.18027e+00_jprb, 1.15716e+00_jprb, &
     & 1.13449e+00_jprb, 1.11227e+00_jprb, 1.09049e+00_jprb, 1.06913e+00_jprb, 1.04819e+00_jprb, &
     & 1.02766e+00_jprb, 1.00754e+00_jprb, 9.87804e-01_jprb, 9.68457e-01_jprb, 9.49490e-01_jprb, &
     & 9.30894e-01_jprb, 9.12662e-01_jprb, 8.94787e-01_jprb, 8.77263e-01_jprb/)
      kao_mco( 1, :, 3) = (/ &
     & 2.55598e+00_jprb, 2.48729e+00_jprb, 2.42045e+00_jprb, 2.35541e+00_jprb, 2.29211e+00_jprb, &
     & 2.23052e+00_jprb, 2.17058e+00_jprb, 2.11225e+00_jprb, 2.05549e+00_jprb, 2.00025e+00_jprb, &
     & 1.94650e+00_jprb, 1.89419e+00_jprb, 1.84329e+00_jprb, 1.79376e+00_jprb, 1.74555e+00_jprb, &
     & 1.69865e+00_jprb, 1.65300e+00_jprb, 1.60858e+00_jprb, 1.56535e+00_jprb/)
      kao_mco( 2, :, 3) = (/ &
     & 2.93113e+00_jprb, 2.85257e+00_jprb, 2.77612e+00_jprb, 2.70172e+00_jprb, 2.62932e+00_jprb, &
     & 2.55885e+00_jprb, 2.49028e+00_jprb, 2.42354e+00_jprb, 2.35859e+00_jprb, 2.29538e+00_jprb, &
     & 2.23386e+00_jprb, 2.17400e+00_jprb, 2.11573e+00_jprb, 2.05903e+00_jprb, 2.00385e+00_jprb, &
     & 1.95015e+00_jprb, 1.89788e+00_jprb, 1.84702e+00_jprb, 1.79752e+00_jprb/)
      kao_mco( 3, :, 3) = (/ &
     & 3.26626e+00_jprb, 3.18025e+00_jprb, 3.09651e+00_jprb, 3.01497e+00_jprb, 2.93558e+00_jprb, &
     & 2.85828e+00_jprb, 2.78302e+00_jprb, 2.70973e+00_jprb, 2.63838e+00_jprb, 2.56891e+00_jprb, &
     & 2.50126e+00_jprb, 2.43540e+00_jprb, 2.37127e+00_jprb, 2.30883e+00_jprb, 2.24803e+00_jprb, &
     & 2.18883e+00_jprb, 2.13120e+00_jprb, 2.07508e+00_jprb, 2.02044e+00_jprb/)
      kao_mco( 4, :, 3) = (/ &
     & 3.65895e+00_jprb, 3.56418e+00_jprb, 3.47187e+00_jprb, 3.38194e+00_jprb, 3.29435e+00_jprb, &
     & 3.20903e+00_jprb, 3.12591e+00_jprb, 3.04495e+00_jprb, 2.96608e+00_jprb, 2.88926e+00_jprb, &
     & 2.81443e+00_jprb, 2.74153e+00_jprb, 2.67053e+00_jprb, 2.60136e+00_jprb, 2.53398e+00_jprb, &
     & 2.46835e+00_jprb, 2.40442e+00_jprb, 2.34214e+00_jprb, 2.28148e+00_jprb/)
      kao_mco( 5, :, 3) = (/ &
     & 4.13692e+00_jprb, 4.03459e+00_jprb, 3.93479e+00_jprb, 3.83746e+00_jprb, 3.74254e+00_jprb, &
     & 3.64997e+00_jprb, 3.55968e+00_jprb, 3.47163e+00_jprb, 3.38576e+00_jprb, 3.30201e+00_jprb, &
     & 3.22034e+00_jprb, 3.14068e+00_jprb, 3.06299e+00_jprb, 2.98723e+00_jprb, 2.91334e+00_jprb, &
     & 2.84128e+00_jprb, 2.77100e+00_jprb, 2.70246e+00_jprb, 2.63561e+00_jprb/)
      kao_mco( 6, :, 3) = (/ &
     & 4.42856e+00_jprb, 4.32480e+00_jprb, 4.22348e+00_jprb, 4.12453e+00_jprb, 4.02790e+00_jprb, &
     & 3.93353e+00_jprb, 3.84137e+00_jprb, 3.75137e+00_jprb, 3.66348e+00_jprb, 3.57765e+00_jprb, &
     & 3.49383e+00_jprb, 3.41198e+00_jprb, 3.33204e+00_jprb, 3.25397e+00_jprb, 3.17774e+00_jprb, &
     & 3.10329e+00_jprb, 3.03058e+00_jprb, 2.95958e+00_jprb, 2.89024e+00_jprb/)
      kao_mco( 7, :, 3) = (/ &
     & 4.31306e+00_jprb, 4.21750e+00_jprb, 4.12406e+00_jprb, 4.03268e+00_jprb, 3.94333e+00_jprb, &
     & 3.85596e+00_jprb, 3.77053e+00_jprb, 3.68699e+00_jprb, 3.60530e+00_jprb, 3.52542e+00_jprb, &
     & 3.44731e+00_jprb, 3.37093e+00_jprb, 3.29624e+00_jprb, 3.22321e+00_jprb, 3.15179e+00_jprb, &
     & 3.08196e+00_jprb, 3.01368e+00_jprb, 2.94691e+00_jprb, 2.88161e+00_jprb/)
      kao_mco( 8, :, 3) = (/ &
     & 4.38922e+00_jprb, 4.32180e+00_jprb, 4.25543e+00_jprb, 4.19007e+00_jprb, 4.12571e+00_jprb, &
     & 4.06234e+00_jprb, 3.99995e+00_jprb, 3.93851e+00_jprb, 3.87802e+00_jprb, 3.81846e+00_jprb, &
     & 3.75981e+00_jprb, 3.70206e+00_jprb, 3.64520e+00_jprb, 3.58922e+00_jprb, 3.53409e+00_jprb, &
     & 3.47981e+00_jprb, 3.42636e+00_jprb, 3.37374e+00_jprb, 3.32192e+00_jprb/)
      kao_mco( 9, :, 3) = (/ &
     & 1.56810e+00_jprb, 1.54211e+00_jprb, 1.51654e+00_jprb, 1.49139e+00_jprb, 1.46667e+00_jprb, &
     & 1.44235e+00_jprb, 1.41844e+00_jprb, 1.39492e+00_jprb, 1.37179e+00_jprb, 1.34905e+00_jprb, &
     & 1.32668e+00_jprb, 1.30469e+00_jprb, 1.28306e+00_jprb, 1.26178e+00_jprb, 1.24086e+00_jprb, &
     & 1.22029e+00_jprb, 1.20006e+00_jprb, 1.18016e+00_jprb, 1.16059e+00_jprb/)
      kao_mco( 1, :, 4) = (/ &
     & 6.58275e+00_jprb, 6.43026e+00_jprb, 6.28130e+00_jprb, 6.13579e+00_jprb, 5.99365e+00_jprb, &
     & 5.85481e+00_jprb, 5.71918e+00_jprb, 5.58669e+00_jprb, 5.45727e+00_jprb, 5.33085e+00_jprb, &
     & 5.20736e+00_jprb, 5.08673e+00_jprb, 4.96889e+00_jprb, 4.85379e+00_jprb, 4.74135e+00_jprb, &
     & 4.63151e+00_jprb, 4.52422e+00_jprb, 4.41942e+00_jprb, 4.31704e+00_jprb/)
      kao_mco( 2, :, 4) = (/ &
     & 6.59883e+00_jprb, 6.45139e+00_jprb, 6.30725e+00_jprb, 6.16633e+00_jprb, 6.02855e+00_jprb, &
     & 5.89386e+00_jprb, 5.76217e+00_jprb, 5.63342e+00_jprb, 5.50756e+00_jprb, 5.38450e+00_jprb, &
     & 5.26419e+00_jprb, 5.14657e+00_jprb, 5.03158e+00_jprb, 4.91916e+00_jprb, 4.80925e+00_jprb, &
     & 4.70180e+00_jprb, 4.59675e+00_jprb, 4.49404e+00_jprb, 4.39363e+00_jprb/)
      kao_mco( 3, :, 4) = (/ &
     & 6.58521e+00_jprb, 6.44452e+00_jprb, 6.30683e+00_jprb, 6.17209e+00_jprb, 6.04023e+00_jprb, &
     & 5.91118e+00_jprb, 5.78489e+00_jprb, 5.66130e+00_jprb, 5.54034e+00_jprb, 5.42198e+00_jprb, &
     & 5.30614e+00_jprb, 5.19277e+00_jprb, 5.08183e+00_jprb, 4.97326e+00_jprb, 4.86701e+00_jprb, &
     & 4.76303e+00_jprb, 4.66127e+00_jprb, 4.56168e+00_jprb, 4.46422e+00_jprb/)
      kao_mco( 4, :, 4) = (/ &
     & 6.33742e+00_jprb, 6.20870e+00_jprb, 6.08260e+00_jprb, 5.95905e+00_jprb, 5.83802e+00_jprb, &
     & 5.71945e+00_jprb, 5.60328e+00_jprb, 5.48947e+00_jprb, 5.37798e+00_jprb, 5.26875e+00_jprb, &
     & 5.16174e+00_jprb, 5.05690e+00_jprb, 4.95419e+00_jprb, 4.85357e+00_jprb, 4.75499e+00_jprb, &
     & 4.65841e+00_jprb, 4.56379e+00_jprb, 4.47110e+00_jprb, 4.38029e+00_jprb/)
      kao_mco( 5, :, 4) = (/ &
     & 5.99732e+00_jprb, 5.88459e+00_jprb, 5.77398e+00_jprb, 5.66545e+00_jprb, 5.55896e+00_jprb, &
     & 5.45447e+00_jprb, 5.35194e+00_jprb, 5.25134e+00_jprb, 5.15263e+00_jprb, 5.05578e+00_jprb, &
     & 4.96075e+00_jprb, 4.86750e+00_jprb, 4.77601e+00_jprb, 4.68623e+00_jprb, 4.59815e+00_jprb, &
     & 4.51171e+00_jprb, 4.42691e+00_jprb, 4.34370e+00_jprb, 4.26205e+00_jprb/)
      kao_mco( 6, :, 4) = (/ &
     & 5.74529e+00_jprb, 5.65249e+00_jprb, 5.56119e+00_jprb, 5.47136e+00_jprb, 5.38299e+00_jprb, &
     & 5.29604e+00_jprb, 5.21049e+00_jprb, 5.12633e+00_jprb, 5.04353e+00_jprb, 4.96206e+00_jprb, &
     & 4.88191e+00_jprb, 4.80306e+00_jprb, 4.72547e+00_jprb, 4.64915e+00_jprb, 4.57405e+00_jprb, &
     & 4.50017e+00_jprb, 4.42748e+00_jprb, 4.35596e+00_jprb, 4.28560e+00_jprb/)
      kao_mco( 7, :, 4) = (/ &
     & 5.87251e+00_jprb, 5.79956e+00_jprb, 5.72753e+00_jprb, 5.65638e+00_jprb, 5.58613e+00_jprb, &
     & 5.51674e+00_jprb, 5.44822e+00_jprb, 5.38054e+00_jprb, 5.31371e+00_jprb, 5.24771e+00_jprb, &
     & 5.18253e+00_jprb, 5.11815e+00_jprb, 5.05458e+00_jprb, 4.99180e+00_jprb, 4.92979e+00_jprb, &
     & 4.86856e+00_jprb, 4.80809e+00_jprb, 4.74837e+00_jprb, 4.68939e+00_jprb/)
      kao_mco( 8, :, 4) = (/ &
     & 5.68503e+00_jprb, 5.62827e+00_jprb, 5.57207e+00_jprb, 5.51644e+00_jprb, 5.46136e+00_jprb, &
     & 5.40684e+00_jprb, 5.35285e+00_jprb, 5.29941e+00_jprb, 5.24650e+00_jprb, 5.19412e+00_jprb, &
     & 5.14226e+00_jprb, 5.09092e+00_jprb, 5.04009e+00_jprb, 4.98977e+00_jprb, 4.93995e+00_jprb, &
     & 4.89063e+00_jprb, 4.84180e+00_jprb, 4.79346e+00_jprb, 4.74560e+00_jprb/)
      kao_mco( 9, :, 4) = (/ &
     & 2.69278e+00_jprb, 2.65058e+00_jprb, 2.60903e+00_jprb, 2.56814e+00_jprb, 2.52789e+00_jprb, &
     & 2.48827e+00_jprb, 2.44927e+00_jprb, 2.41088e+00_jprb, 2.37310e+00_jprb, 2.33590e+00_jprb, &
     & 2.29929e+00_jprb, 2.26325e+00_jprb, 2.22778e+00_jprb, 2.19286e+00_jprb, 2.15849e+00_jprb, &
     & 2.12466e+00_jprb, 2.09136e+00_jprb, 2.05859e+00_jprb, 2.02632e+00_jprb/)
      kao_mco( 1, :, 5) = (/ &
     & 9.12231e+00_jprb, 9.00052e+00_jprb, 8.88036e+00_jprb, 8.76180e+00_jprb, 8.64482e+00_jprb, &
     & 8.52941e+00_jprb, 8.41553e+00_jprb, 8.30318e+00_jprb, 8.19233e+00_jprb, 8.08295e+00_jprb, &
     & 7.97504e+00_jprb, 7.86857e+00_jprb, 7.76352e+00_jprb, 7.65987e+00_jprb, 7.55760e+00_jprb, &
     & 7.45671e+00_jprb, 7.35715e+00_jprb, 7.25893e+00_jprb, 7.16202e+00_jprb/)
      kao_mco( 2, :, 5) = (/ &
     & 8.37315e+00_jprb, 8.27808e+00_jprb, 8.18410e+00_jprb, 8.09118e+00_jprb, 7.99931e+00_jprb, &
     & 7.90849e+00_jprb, 7.81871e+00_jprb, 7.72994e+00_jprb, 7.64217e+00_jprb, 7.55541e+00_jprb, &
     & 7.46963e+00_jprb, 7.38482e+00_jprb, 7.30098e+00_jprb, 7.21809e+00_jprb, 7.13614e+00_jprb, &
     & 7.05512e+00_jprb, 6.97502e+00_jprb, 6.89582e+00_jprb, 6.81753e+00_jprb/)
      kao_mco( 3, :, 5) = (/ &
     & 8.14557e+00_jprb, 8.06533e+00_jprb, 7.98587e+00_jprb, 7.90720e+00_jprb, 7.82930e+00_jprb, &
     & 7.75217e+00_jprb, 7.67580e+00_jprb, 7.60018e+00_jprb, 7.52530e+00_jprb, 7.45117e+00_jprb, &
     & 7.37776e+00_jprb, 7.30508e+00_jprb, 7.23311e+00_jprb, 7.16186e+00_jprb, 7.09130e+00_jprb, &
     & 7.02144e+00_jprb, 6.95227e+00_jprb, 6.88378e+00_jprb, 6.81596e+00_jprb/)
      kao_mco( 4, :, 5) = (/ &
     & 8.18046e+00_jprb, 8.11056e+00_jprb, 8.04126e+00_jprb, 7.97256e+00_jprb, 7.90444e+00_jprb, &
     & 7.83690e+00_jprb, 7.76994e+00_jprb, 7.70355e+00_jprb, 7.63773e+00_jprb, 7.57247e+00_jprb, &
     & 7.50777e+00_jprb, 7.44362e+00_jprb, 7.38002e+00_jprb, 7.31697e+00_jprb, 7.25445e+00_jprb, &
     & 7.19247e+00_jprb, 7.13101e+00_jprb, 7.07008e+00_jprb, 7.00968e+00_jprb/)
      kao_mco( 5, :, 5) = (/ &
     & 8.30092e+00_jprb, 8.23529e+00_jprb, 8.17019e+00_jprb, 8.10559e+00_jprb, 8.04151e+00_jprb, &
     & 7.97794e+00_jprb, 7.91487e+00_jprb, 7.85230e+00_jprb, 7.79022e+00_jprb, 7.72863e+00_jprb, &
     & 7.66753e+00_jprb, 7.60691e+00_jprb, 7.54678e+00_jprb, 7.48711e+00_jprb, 7.42792e+00_jprb, &
     & 7.36920e+00_jprb, 7.31094e+00_jprb, 7.25314e+00_jprb, 7.19580e+00_jprb/)
      kao_mco( 6, :, 5) = (/ &
     & 8.30014e+00_jprb, 8.24466e+00_jprb, 8.18955e+00_jprb, 8.13481e+00_jprb, 8.08044e+00_jprb, &
     & 8.02642e+00_jprb, 7.97277e+00_jprb, 7.91948e+00_jprb, 7.86655e+00_jprb, 7.81396e+00_jprb, &
     & 7.76173e+00_jprb, 7.70985e+00_jprb, 7.65832e+00_jprb, 7.60713e+00_jprb, 7.55628e+00_jprb, &
     & 7.50577e+00_jprb, 7.45560e+00_jprb, 7.40577e+00_jprb, 7.35627e+00_jprb/)
      kao_mco( 7, :, 5) = (/ &
     & 7.95931e+00_jprb, 7.93958e+00_jprb, 7.91989e+00_jprb, 7.90025e+00_jprb, 7.88066e+00_jprb, &
     & 7.86112e+00_jprb, 7.84163e+00_jprb, 7.82219e+00_jprb, 7.80279e+00_jprb, 7.78344e+00_jprb, &
     & 7.76414e+00_jprb, 7.74489e+00_jprb, 7.72568e+00_jprb, 7.70653e+00_jprb, 7.68742e+00_jprb, &
     & 7.66836e+00_jprb, 7.64934e+00_jprb, 7.63038e+00_jprb, 7.61146e+00_jprb/)
      kao_mco( 8, :, 5) = (/ &
     & 9.32576e+00_jprb, 9.31747e+00_jprb, 9.30919e+00_jprb, 9.30092e+00_jprb, 9.29265e+00_jprb, &
     & 9.28439e+00_jprb, 9.27613e+00_jprb, 9.26789e+00_jprb, 9.25965e+00_jprb, 9.25142e+00_jprb, &
     & 9.24320e+00_jprb, 9.23498e+00_jprb, 9.22677e+00_jprb, 9.21857e+00_jprb, 9.21038e+00_jprb, &
     & 9.20219e+00_jprb, 9.19401e+00_jprb, 9.18584e+00_jprb, 9.17767e+00_jprb/)
      kao_mco( 9, :, 5) = (/ &
     & 4.13116e+00_jprb, 4.08426e+00_jprb, 4.03788e+00_jprb, 3.99204e+00_jprb, 3.94671e+00_jprb, &
     & 3.90190e+00_jprb, 3.85760e+00_jprb, 3.81380e+00_jprb, 3.77049e+00_jprb, 3.72768e+00_jprb, &
     & 3.68536e+00_jprb, 3.64351e+00_jprb, 3.60214e+00_jprb, 3.56124e+00_jprb, 3.52081e+00_jprb, &
     & 3.48083e+00_jprb, 3.44131e+00_jprb, 3.40224e+00_jprb, 3.36361e+00_jprb/)
      kao_mco( 1, :, 6) = (/ &
     & 1.21200e+01_jprb, 1.21580e+01_jprb, 1.21961e+01_jprb, 1.22344e+01_jprb, 1.22728e+01_jprb, &
     & 1.23113e+01_jprb, 1.23499e+01_jprb, 1.23886e+01_jprb, 1.24275e+01_jprb, 1.24664e+01_jprb, &
     & 1.25056e+01_jprb, 1.25448e+01_jprb, 1.25841e+01_jprb, 1.26236e+01_jprb, 1.26632e+01_jprb, &
     & 1.27029e+01_jprb, 1.27428e+01_jprb, 1.27827e+01_jprb, 1.28228e+01_jprb/)
      kao_mco( 2, :, 6) = (/ &
     & 1.25231e+01_jprb, 1.25625e+01_jprb, 1.26020e+01_jprb, 1.26417e+01_jprb, 1.26815e+01_jprb, &
     & 1.27214e+01_jprb, 1.27614e+01_jprb, 1.28015e+01_jprb, 1.28418e+01_jprb, 1.28822e+01_jprb, &
     & 1.29228e+01_jprb, 1.29634e+01_jprb, 1.30042e+01_jprb, 1.30451e+01_jprb, 1.30862e+01_jprb, &
     & 1.31274e+01_jprb, 1.31687e+01_jprb, 1.32101e+01_jprb, 1.32517e+01_jprb/)
      kao_mco( 3, :, 6) = (/ &
     & 1.27566e+01_jprb, 1.27983e+01_jprb, 1.28401e+01_jprb, 1.28820e+01_jprb, 1.29241e+01_jprb, &
     & 1.29663e+01_jprb, 1.30087e+01_jprb, 1.30512e+01_jprb, 1.30938e+01_jprb, 1.31366e+01_jprb, &
     & 1.31795e+01_jprb, 1.32225e+01_jprb, 1.32657e+01_jprb, 1.33090e+01_jprb, 1.33525e+01_jprb, &
     & 1.33961e+01_jprb, 1.34399e+01_jprb, 1.34838e+01_jprb, 1.35278e+01_jprb/)
      kao_mco( 4, :, 6) = (/ &
     & 1.27132e+01_jprb, 1.27454e+01_jprb, 1.27777e+01_jprb, 1.28101e+01_jprb, 1.28425e+01_jprb, &
     & 1.28750e+01_jprb, 1.29077e+01_jprb, 1.29403e+01_jprb, 1.29731e+01_jprb, 1.30060e+01_jprb, &
     & 1.30389e+01_jprb, 1.30720e+01_jprb, 1.31051e+01_jprb, 1.31383e+01_jprb, 1.31716e+01_jprb, &
     & 1.32049e+01_jprb, 1.32384e+01_jprb, 1.32719e+01_jprb, 1.33055e+01_jprb/)
      kao_mco( 5, :, 6) = (/ &
     & 1.33151e+01_jprb, 1.33523e+01_jprb, 1.33896e+01_jprb, 1.34271e+01_jprb, 1.34646e+01_jprb, &
     & 1.35022e+01_jprb, 1.35400e+01_jprb, 1.35779e+01_jprb, 1.36158e+01_jprb, 1.36539e+01_jprb, &
     & 1.36921e+01_jprb, 1.37303e+01_jprb, 1.37687e+01_jprb, 1.38072e+01_jprb, 1.38458e+01_jprb, &
     & 1.38846e+01_jprb, 1.39234e+01_jprb, 1.39623e+01_jprb, 1.40013e+01_jprb/)
      kao_mco( 6, :, 6) = (/ &
     & 1.41448e+01_jprb, 1.41902e+01_jprb, 1.42357e+01_jprb, 1.42814e+01_jprb, 1.43272e+01_jprb, &
     & 1.43732e+01_jprb, 1.44194e+01_jprb, 1.44656e+01_jprb, 1.45121e+01_jprb, 1.45586e+01_jprb, &
     & 1.46054e+01_jprb, 1.46522e+01_jprb, 1.46993e+01_jprb, 1.47464e+01_jprb, 1.47938e+01_jprb, &
     & 1.48413e+01_jprb, 1.48889e+01_jprb, 1.49367e+01_jprb, 1.49846e+01_jprb/)
      kao_mco( 7, :, 6) = (/ &
     & 1.56578e+01_jprb, 1.56938e+01_jprb, 1.57299e+01_jprb, 1.57661e+01_jprb, 1.58024e+01_jprb, &
     & 1.58388e+01_jprb, 1.58752e+01_jprb, 1.59117e+01_jprb, 1.59484e+01_jprb, 1.59851e+01_jprb, &
     & 1.60218e+01_jprb, 1.60587e+01_jprb, 1.60957e+01_jprb, 1.61327e+01_jprb, 1.61698e+01_jprb, &
     & 1.62070e+01_jprb, 1.62443e+01_jprb, 1.62817e+01_jprb, 1.63192e+01_jprb/)
      kao_mco( 8, :, 6) = (/ &
     & 1.73627e+01_jprb, 1.74761e+01_jprb, 1.75903e+01_jprb, 1.77052e+01_jprb, 1.78208e+01_jprb, &
     & 1.79373e+01_jprb, 1.80544e+01_jprb, 1.81724e+01_jprb, 1.82911e+01_jprb, 1.84106e+01_jprb, &
     & 1.85309e+01_jprb, 1.86519e+01_jprb, 1.87738e+01_jprb, 1.88964e+01_jprb, 1.90198e+01_jprb, &
     & 1.91441e+01_jprb, 1.92692e+01_jprb, 1.93950e+01_jprb, 1.95217e+01_jprb/)
      kao_mco( 9, :, 6) = (/ &
     & 7.16326e+00_jprb, 7.12921e+00_jprb, 7.09531e+00_jprb, 7.06158e+00_jprb, 7.02800e+00_jprb, &
     & 6.99459e+00_jprb, 6.96133e+00_jprb, 6.92824e+00_jprb, 6.89530e+00_jprb, 6.86252e+00_jprb, &
     & 6.82989e+00_jprb, 6.79742e+00_jprb, 6.76510e+00_jprb, 6.73293e+00_jprb, 6.70092e+00_jprb, &
     & 6.66906e+00_jprb, 6.63736e+00_jprb, 6.60580e+00_jprb, 6.57439e+00_jprb/)
      kao_mco( 1, :, 7) = (/ &
     & 2.09288e+01_jprb, 2.10487e+01_jprb, 2.11692e+01_jprb, 2.12904e+01_jprb, 2.14124e+01_jprb, &
     & 2.15350e+01_jprb, 2.16583e+01_jprb, 2.17823e+01_jprb, 2.19070e+01_jprb, 2.20325e+01_jprb, &
     & 2.21587e+01_jprb, 2.22855e+01_jprb, 2.24132e+01_jprb, 2.25415e+01_jprb, 2.26706e+01_jprb, &
     & 2.28004e+01_jprb, 2.29310e+01_jprb, 2.30623e+01_jprb, 2.31943e+01_jprb/)
      kao_mco( 2, :, 7) = (/ &
     & 2.08159e+01_jprb, 2.09509e+01_jprb, 2.10867e+01_jprb, 2.12234e+01_jprb, 2.13610e+01_jprb, &
     & 2.14994e+01_jprb, 2.16388e+01_jprb, 2.17791e+01_jprb, 2.19202e+01_jprb, 2.20623e+01_jprb, &
     & 2.22053e+01_jprb, 2.23493e+01_jprb, 2.24942e+01_jprb, 2.26400e+01_jprb, 2.27867e+01_jprb, &
     & 2.29345e+01_jprb, 2.30831e+01_jprb, 2.32328e+01_jprb, 2.33834e+01_jprb/)
      kao_mco( 3, :, 7) = (/ &
     & 2.10827e+01_jprb, 2.12409e+01_jprb, 2.14002e+01_jprb, 2.15608e+01_jprb, 2.17225e+01_jprb, &
     & 2.18855e+01_jprb, 2.20497e+01_jprb, 2.22151e+01_jprb, 2.23818e+01_jprb, 2.25497e+01_jprb, &
     & 2.27189e+01_jprb, 2.28893e+01_jprb, 2.30611e+01_jprb, 2.32341e+01_jprb, 2.34084e+01_jprb, &
     & 2.35840e+01_jprb, 2.37609e+01_jprb, 2.39392e+01_jprb, 2.41188e+01_jprb/)
      kao_mco( 4, :, 7) = (/ &
     & 2.13866e+01_jprb, 2.15772e+01_jprb, 2.17694e+01_jprb, 2.19634e+01_jprb, 2.21590e+01_jprb, &
     & 2.23565e+01_jprb, 2.25556e+01_jprb, 2.27566e+01_jprb, 2.29594e+01_jprb, 2.31639e+01_jprb, &
     & 2.33703e+01_jprb, 2.35785e+01_jprb, 2.37886e+01_jprb, 2.40005e+01_jprb, 2.42144e+01_jprb, &
     & 2.44301e+01_jprb, 2.46477e+01_jprb, 2.48673e+01_jprb, 2.50889e+01_jprb/)
      kao_mco( 5, :, 7) = (/ &
     & 1.93714e+01_jprb, 1.95595e+01_jprb, 1.97493e+01_jprb, 1.99410e+01_jprb, 2.01345e+01_jprb, &
     & 2.03300e+01_jprb, 2.05273e+01_jprb, 2.07265e+01_jprb, 2.09277e+01_jprb, 2.11308e+01_jprb, &
     & 2.13359e+01_jprb, 2.15430e+01_jprb, 2.17521e+01_jprb, 2.19632e+01_jprb, 2.21764e+01_jprb, &
     & 2.23917e+01_jprb, 2.26090e+01_jprb, 2.28284e+01_jprb, 2.30500e+01_jprb/)
      kao_mco( 6, :, 7) = (/ &
     & 1.70418e+01_jprb, 1.72109e+01_jprb, 1.73817e+01_jprb, 1.75541e+01_jprb, 1.77283e+01_jprb, &
     & 1.79041e+01_jprb, 1.80818e+01_jprb, 1.82612e+01_jprb, 1.84423e+01_jprb, 1.86253e+01_jprb, &
     & 1.88101e+01_jprb, 1.89967e+01_jprb, 1.91852e+01_jprb, 1.93755e+01_jprb, 1.95678e+01_jprb, &
     & 1.97619e+01_jprb, 1.99580e+01_jprb, 2.01560e+01_jprb, 2.03560e+01_jprb/)
      kao_mco( 7, :, 7) = (/ &
     & 1.31735e+01_jprb, 1.32921e+01_jprb, 1.34118e+01_jprb, 1.35326e+01_jprb, 1.36545e+01_jprb, &
     & 1.37775e+01_jprb, 1.39015e+01_jprb, 1.40267e+01_jprb, 1.41531e+01_jprb, 1.42805e+01_jprb, &
     & 1.44091e+01_jprb, 1.45389e+01_jprb, 1.46698e+01_jprb, 1.48019e+01_jprb, 1.49353e+01_jprb, &
     & 1.50698e+01_jprb, 1.52055e+01_jprb, 1.53424e+01_jprb, 1.54806e+01_jprb/)
      kao_mco( 8, :, 7) = (/ &
     & 4.97361e+00_jprb, 4.96550e+00_jprb, 4.95740e+00_jprb, 4.94931e+00_jprb, 4.94124e+00_jprb, &
     & 4.93318e+00_jprb, 4.92513e+00_jprb, 4.91709e+00_jprb, 4.90907e+00_jprb, 4.90106e+00_jprb, &
     & 4.89307e+00_jprb, 4.88509e+00_jprb, 4.87712e+00_jprb, 4.86916e+00_jprb, 4.86122e+00_jprb, &
     & 4.85329e+00_jprb, 4.84537e+00_jprb, 4.83747e+00_jprb, 4.82958e+00_jprb/)
      kao_mco( 9, :, 7) = (/ &
     & 1.76121e+01_jprb, 1.75887e+01_jprb, 1.75653e+01_jprb, 1.75420e+01_jprb, 1.75187e+01_jprb, &
     & 1.74955e+01_jprb, 1.74722e+01_jprb, 1.74490e+01_jprb, 1.74259e+01_jprb, 1.74027e+01_jprb, &
     & 1.73796e+01_jprb, 1.73566e+01_jprb, 1.73335e+01_jprb, 1.73105e+01_jprb, 1.72875e+01_jprb, &
     & 1.72646e+01_jprb, 1.72416e+01_jprb, 1.72188e+01_jprb, 1.71959e+01_jprb/)
      kao_mco( 1, :, 8) = (/ &
     & 5.99126e+00_jprb, 6.08386e+00_jprb, 6.17790e+00_jprb, 6.27339e+00_jprb, 6.37035e+00_jprb, &
     & 6.46881e+00_jprb, 6.56880e+00_jprb, 6.67033e+00_jprb, 6.77343e+00_jprb, 6.87812e+00_jprb, &
     & 6.98443e+00_jprb, 7.09238e+00_jprb, 7.20201e+00_jprb, 7.31332e+00_jprb, 7.42636e+00_jprb, &
     & 7.54115e+00_jprb, 7.65771e+00_jprb, 7.77607e+00_jprb, 7.89626e+00_jprb/)
      kao_mco( 2, :, 8) = (/ &
     & 4.71621e+00_jprb, 4.78830e+00_jprb, 4.86149e+00_jprb, 4.93580e+00_jprb, 5.01124e+00_jprb, &
     & 5.08784e+00_jprb, 5.16561e+00_jprb, 5.24456e+00_jprb, 5.32473e+00_jprb, 5.40612e+00_jprb, &
     & 5.48875e+00_jprb, 5.57264e+00_jprb, 5.65782e+00_jprb, 5.74430e+00_jprb, 5.83211e+00_jprb, &
     & 5.92125e+00_jprb, 6.01176e+00_jprb, 6.10365e+00_jprb, 6.19694e+00_jprb/)
      kao_mco( 3, :, 8) = (/ &
     & 2.77067e+00_jprb, 2.81437e+00_jprb, 2.85876e+00_jprb, 2.90385e+00_jprb, 2.94965e+00_jprb, &
     & 2.99617e+00_jprb, 3.04343e+00_jprb, 3.09143e+00_jprb, 3.14019e+00_jprb, 3.18972e+00_jprb, &
     & 3.24003e+00_jprb, 3.29114e+00_jprb, 3.34305e+00_jprb, 3.39578e+00_jprb, 3.44934e+00_jprb, &
     & 3.50374e+00_jprb, 3.55901e+00_jprb, 3.61514e+00_jprb, 3.67216e+00_jprb/)
      kao_mco( 4, :, 8) = (/ &
     & 1.22388e+00_jprb, 1.24248e+00_jprb, 1.26136e+00_jprb, 1.28053e+00_jprb, 1.29999e+00_jprb, &
     & 1.31974e+00_jprb, 1.33979e+00_jprb, 1.36015e+00_jprb, 1.38082e+00_jprb, 1.40181e+00_jprb, &
     & 1.42311e+00_jprb, 1.44473e+00_jprb, 1.46669e+00_jprb, 1.48897e+00_jprb, 1.51160e+00_jprb, &
     & 1.53457e+00_jprb, 1.55789e+00_jprb, 1.58156e+00_jprb, 1.60559e+00_jprb/)
      kao_mco( 5, :, 8) = (/ &
     & 1.41479e+00_jprb, 1.43540e+00_jprb, 1.45631e+00_jprb, 1.47752e+00_jprb, 1.49904e+00_jprb, &
     & 1.52088e+00_jprb, 1.54303e+00_jprb, 1.56550e+00_jprb, 1.58831e+00_jprb, 1.61144e+00_jprb, &
     & 1.63491e+00_jprb, 1.65872e+00_jprb, 1.68288e+00_jprb, 1.70740e+00_jprb, 1.73227e+00_jprb, &
     & 1.75750e+00_jprb, 1.78310e+00_jprb, 1.80907e+00_jprb, 1.83542e+00_jprb/)
      kao_mco( 6, :, 8) = (/ &
     & 1.43154e+00_jprb, 1.46074e+00_jprb, 1.49053e+00_jprb, 1.52093e+00_jprb, 1.55196e+00_jprb, &
     & 1.58361e+00_jprb, 1.61591e+00_jprb, 1.64887e+00_jprb, 1.68250e+00_jprb, 1.71682e+00_jprb, &
     & 1.75184e+00_jprb, 1.78757e+00_jprb, 1.82403e+00_jprb, 1.86123e+00_jprb, 1.89920e+00_jprb, &
     & 1.93793e+00_jprb, 1.97746e+00_jprb, 2.01779e+00_jprb, 2.05895e+00_jprb/)
      kao_mco( 7, :, 8) = (/ &
     & 2.49358e+00_jprb, 2.56028e+00_jprb, 2.62875e+00_jprb, 2.69906e+00_jprb, 2.77124e+00_jprb, &
     & 2.84536e+00_jprb, 2.92146e+00_jprb, 2.99960e+00_jprb, 3.07982e+00_jprb, 3.16219e+00_jprb, &
     & 3.24677e+00_jprb, 3.33360e+00_jprb, 3.42276e+00_jprb, 3.51430e+00_jprb, 3.60829e+00_jprb, &
     & 3.70480e+00_jprb, 3.80388e+00_jprb, 3.90562e+00_jprb, 4.01007e+00_jprb/)
      kao_mco( 8, :, 8) = (/ &
     & 4.32513e+00_jprb, 4.39903e+00_jprb, 4.47420e+00_jprb, 4.55065e+00_jprb, 4.62841e+00_jprb, &
     & 4.70750e+00_jprb, 4.78794e+00_jprb, 4.86975e+00_jprb, 4.95296e+00_jprb, 5.03759e+00_jprb, &
     & 5.12367e+00_jprb, 5.21122e+00_jprb, 5.30027e+00_jprb, 5.39084e+00_jprb, 5.48295e+00_jprb, &
     & 5.57664e+00_jprb, 5.67193e+00_jprb, 5.76885e+00_jprb, 5.86743e+00_jprb/)
      kao_mco( 9, :, 8) = (/ &
     & 3.35160e+01_jprb, 3.36789e+01_jprb, 3.38425e+01_jprb, 3.40069e+01_jprb, 3.41722e+01_jprb, &
     & 3.43382e+01_jprb, 3.45050e+01_jprb, 3.46727e+01_jprb, 3.48412e+01_jprb, 3.50105e+01_jprb, &
     & 3.51806e+01_jprb, 3.53515e+01_jprb, 3.55233e+01_jprb, 3.56959e+01_jprb, 3.58693e+01_jprb, &
     & 3.60436e+01_jprb, 3.62187e+01_jprb, 3.63947e+01_jprb, 3.65715e+01_jprb/)
      kao_mco( 1, :, 9) = (/ &
     & 8.68159e-01_jprb, 9.13680e-01_jprb, 9.61587e-01_jprb, 1.01201e+00_jprb, 1.06507e+00_jprb, &
     & 1.12091e+00_jprb, 1.17969e+00_jprb, 1.24154e+00_jprb, 1.30664e+00_jprb, 1.37515e+00_jprb, &
     & 1.44726e+00_jprb, 1.52314e+00_jprb, 1.60300e+00_jprb, 1.68705e+00_jprb, 1.77551e+00_jprb, &
     & 1.86861e+00_jprb, 1.96658e+00_jprb, 2.06970e+00_jprb, 2.17822e+00_jprb/)
      kao_mco( 2, :, 9) = (/ &
     & 9.04391e-01_jprb, 9.49669e-01_jprb, 9.97214e-01_jprb, 1.04714e+00_jprb, 1.09956e+00_jprb, &
     & 1.15461e+00_jprb, 1.21242e+00_jprb, 1.27312e+00_jprb, 1.33685e+00_jprb, 1.40378e+00_jprb, &
     & 1.47406e+00_jprb, 1.54786e+00_jprb, 1.62535e+00_jprb, 1.70673e+00_jprb, 1.79217e+00_jprb, &
     & 1.88190e+00_jprb, 1.97611e+00_jprb, 2.07505e+00_jprb, 2.17893e+00_jprb/)
      kao_mco( 3, :, 9) = (/ &
     & 9.67479e-01_jprb, 1.01312e+00_jprb, 1.06092e+00_jprb, 1.11098e+00_jprb, 1.16339e+00_jprb, &
     & 1.21828e+00_jprb, 1.27576e+00_jprb, 1.33595e+00_jprb, 1.39898e+00_jprb, 1.46499e+00_jprb, &
     & 1.53411e+00_jprb, 1.60649e+00_jprb, 1.68228e+00_jprb, 1.76165e+00_jprb, 1.84476e+00_jprb, &
     & 1.93180e+00_jprb, 2.02294e+00_jprb, 2.11839e+00_jprb, 2.21833e+00_jprb/)
      kao_mco( 4, :, 9) = (/ &
     & 1.05240e+00_jprb, 1.09817e+00_jprb, 1.14592e+00_jprb, 1.19576e+00_jprb, 1.24776e+00_jprb, &
     & 1.30202e+00_jprb, 1.35864e+00_jprb, 1.41772e+00_jprb, 1.47937e+00_jprb, 1.54371e+00_jprb, &
     & 1.61084e+00_jprb, 1.68089e+00_jprb, 1.75398e+00_jprb, 1.83026e+00_jprb, 1.90985e+00_jprb, &
     & 1.99290e+00_jprb, 2.07957e+00_jprb, 2.17000e+00_jprb, 2.26437e+00_jprb/)
      kao_mco( 5, :, 9) = (/ &
     & 1.25800e+00_jprb, 1.30557e+00_jprb, 1.35494e+00_jprb, 1.40618e+00_jprb, 1.45935e+00_jprb, &
     & 1.51454e+00_jprb, 1.57181e+00_jprb, 1.63125e+00_jprb, 1.69293e+00_jprb, 1.75695e+00_jprb, &
     & 1.82339e+00_jprb, 1.89234e+00_jprb, 1.96389e+00_jprb, 2.03816e+00_jprb, 2.11523e+00_jprb, &
     & 2.19522e+00_jprb, 2.27823e+00_jprb, 2.36438e+00_jprb, 2.45378e+00_jprb/)
      kao_mco( 6, :, 9) = (/ &
     & 1.76509e+00_jprb, 1.80550e+00_jprb, 1.84683e+00_jprb, 1.88911e+00_jprb, 1.93235e+00_jprb, &
     & 1.97658e+00_jprb, 2.02183e+00_jprb, 2.06811e+00_jprb, 2.11546e+00_jprb, 2.16388e+00_jprb, &
     & 2.21342e+00_jprb, 2.26408e+00_jprb, 2.31591e+00_jprb, 2.36893e+00_jprb, 2.42315e+00_jprb, &
     & 2.47862e+00_jprb, 2.53536e+00_jprb, 2.59340e+00_jprb, 2.65277e+00_jprb/)
      kao_mco( 7, :, 9) = (/ &
     & 2.03543e+00_jprb, 2.05285e+00_jprb, 2.07042e+00_jprb, 2.08815e+00_jprb, 2.10602e+00_jprb, &
     & 2.12405e+00_jprb, 2.14223e+00_jprb, 2.16057e+00_jprb, 2.17907e+00_jprb, 2.19772e+00_jprb, &
     & 2.21654e+00_jprb, 2.23551e+00_jprb, 2.25465e+00_jprb, 2.27395e+00_jprb, 2.29342e+00_jprb, &
     & 2.31305e+00_jprb, 2.33285e+00_jprb, 2.35282e+00_jprb, 2.37296e+00_jprb/)
      kao_mco( 8, :, 9) = (/ &
     & 3.18883e+00_jprb, 3.20538e+00_jprb, 3.22200e+00_jprb, 3.23872e+00_jprb, 3.25552e+00_jprb, &
     & 3.27241e+00_jprb, 3.28939e+00_jprb, 3.30645e+00_jprb, 3.32360e+00_jprb, 3.34085e+00_jprb, &
     & 3.35818e+00_jprb, 3.37560e+00_jprb, 3.39311e+00_jprb, 3.41071e+00_jprb, 3.42841e+00_jprb, &
     & 3.44619e+00_jprb, 3.46407e+00_jprb, 3.48204e+00_jprb, 3.50011e+00_jprb/)
      kao_mco( 9, :, 9) = (/ &
     & 3.97585e+00_jprb, 3.96333e+00_jprb, 3.95084e+00_jprb, 3.93839e+00_jprb, 3.92598e+00_jprb, &
     & 3.91361e+00_jprb, 3.90128e+00_jprb, 3.88899e+00_jprb, 3.87674e+00_jprb, 3.86452e+00_jprb, &
     & 3.85235e+00_jprb, 3.84021e+00_jprb, 3.82811e+00_jprb, 3.81605e+00_jprb, 3.80402e+00_jprb, &
     & 3.79204e+00_jprb, 3.78009e+00_jprb, 3.76818e+00_jprb, 3.75631e+00_jprb/)
      kao_mco( 1, :,10) = (/ &
     & 8.62646e-01_jprb, 9.35164e-01_jprb, 1.01378e+00_jprb, 1.09900e+00_jprb, 1.19139e+00_jprb, &
     & 1.29154e+00_jprb, 1.40011e+00_jprb, 1.51781e+00_jprb, 1.64541e+00_jprb, 1.78373e+00_jprb, &
     & 1.93367e+00_jprb, 2.09623e+00_jprb, 2.27245e+00_jprb, 2.46348e+00_jprb, 2.67057e+00_jprb, &
     & 2.89507e+00_jprb, 3.13844e+00_jprb, 3.40227e+00_jprb, 3.68828e+00_jprb/)
      kao_mco( 2, :,10) = (/ &
     & 8.04693e-01_jprb, 8.72167e-01_jprb, 9.45298e-01_jprb, 1.02456e+00_jprb, 1.11047e+00_jprb, &
     & 1.20358e+00_jprb, 1.30450e+00_jprb, 1.41389e+00_jprb, 1.53244e+00_jprb, 1.66094e+00_jprb, &
     & 1.80021e+00_jprb, 1.95115e+00_jprb, 2.11476e+00_jprb, 2.29208e+00_jprb, 2.48427e+00_jprb, &
     & 2.69258e+00_jprb, 2.91835e+00_jprb, 3.16305e+00_jprb, 3.42827e+00_jprb/)
      kao_mco( 3, :,10) = (/ &
     & 7.66566e-01_jprb, 8.24651e-01_jprb, 8.87137e-01_jprb, 9.54358e-01_jprb, 1.02667e+00_jprb, &
     & 1.10447e+00_jprb, 1.18815e+00_jprb, 1.27818e+00_jprb, 1.37503e+00_jprb, 1.47922e+00_jprb, &
     & 1.59131e+00_jprb, 1.71189e+00_jprb, 1.84160e+00_jprb, 1.98114e+00_jprb, 2.13126e+00_jprb, &
     & 2.29275e+00_jprb, 2.46648e+00_jprb, 2.65337e+00_jprb, 2.85442e+00_jprb/)
      kao_mco( 4, :,10) = (/ &
     & 5.40305e-01_jprb, 5.77106e-01_jprb, 6.16414e-01_jprb, 6.58400e-01_jprb, 7.03245e-01_jprb, &
     & 7.51145e-01_jprb, 8.02307e-01_jprb, 8.56954e-01_jprb, 9.15323e-01_jprb, 9.77668e-01_jprb, &
     & 1.04426e+00_jprb, 1.11539e+00_jprb, 1.19136e+00_jprb, 1.27250e+00_jprb, 1.35918e+00_jprb, &
     & 1.45175e+00_jprb, 1.55064e+00_jprb, 1.65625e+00_jprb, 1.76906e+00_jprb/)
      kao_mco( 5, :,10) = (/ &
     & 8.22474e-01_jprb, 8.02911e-01_jprb, 7.83814e-01_jprb, 7.65171e-01_jprb, 7.46971e-01_jprb, &
     & 7.29204e-01_jprb, 7.11860e-01_jprb, 6.94928e-01_jprb, 6.78399e-01_jprb, 6.62263e-01_jprb, &
     & 6.46511e-01_jprb, 6.31133e-01_jprb, 6.16122e-01_jprb, 6.01467e-01_jprb, 5.87161e-01_jprb, &
     & 5.73195e-01_jprb, 5.59562e-01_jprb, 5.46252e-01_jprb, 5.33260e-01_jprb/)
      kao_mco( 6, :,10) = (/ &
     & 1.28162e+00_jprb, 1.25110e+00_jprb, 1.22131e+00_jprb, 1.19223e+00_jprb, 1.16384e+00_jprb, &
     & 1.13613e+00_jprb, 1.10908e+00_jprb, 1.08267e+00_jprb, 1.05689e+00_jprb, 1.03173e+00_jprb, &
     & 1.00716e+00_jprb, 9.83184e-01_jprb, 9.59774e-01_jprb, 9.36921e-01_jprb, 9.14613e-01_jprb, &
     & 8.92836e-01_jprb, 8.71577e-01_jprb, 8.50825e-01_jprb, 8.30567e-01_jprb/)
      kao_mco( 7, :,10) = (/ &
     & 1.92679e+00_jprb, 1.90551e+00_jprb, 1.88446e+00_jprb, 1.86365e+00_jprb, 1.84307e+00_jprb, &
     & 1.82271e+00_jprb, 1.80258e+00_jprb, 1.78267e+00_jprb, 1.76298e+00_jprb, 1.74351e+00_jprb, &
     & 1.72425e+00_jprb, 1.70520e+00_jprb, 1.68637e+00_jprb, 1.66774e+00_jprb, 1.64932e+00_jprb, &
     & 1.63111e+00_jprb, 1.61309e+00_jprb, 1.59527e+00_jprb, 1.57765e+00_jprb/)
      kao_mco( 8, :,10) = (/ &
     & 4.66485e+00_jprb, 4.60869e+00_jprb, 4.55320e+00_jprb, 4.49838e+00_jprb, 4.44423e+00_jprb, &
     & 4.39072e+00_jprb, 4.33786e+00_jprb, 4.28563e+00_jprb, 4.23404e+00_jprb, 4.18306e+00_jprb, &
     & 4.13270e+00_jprb, 4.08295e+00_jprb, 4.03379e+00_jprb, 3.98523e+00_jprb, 3.93725e+00_jprb, &
     & 3.88985e+00_jprb, 3.84301e+00_jprb, 3.79675e+00_jprb, 3.75104e+00_jprb/)
      kao_mco( 9, :,10) = (/ &
     & 1.41505e+00_jprb, 1.37820e+00_jprb, 1.34232e+00_jprb, 1.30736e+00_jprb, 1.27332e+00_jprb, &
     & 1.24016e+00_jprb, 1.20786e+00_jprb, 1.17641e+00_jprb, 1.14578e+00_jprb, 1.11594e+00_jprb, &
     & 1.08688e+00_jprb, 1.05858e+00_jprb, 1.03101e+00_jprb, 1.00416e+00_jprb, 9.78015e-01_jprb, &
     & 9.52547e-01_jprb, 9.27742e-01_jprb, 9.03583e-01_jprb, 8.80053e-01_jprb/)
      kao_mco( 1, :,11) = (/ &
     & 3.40468e-03_jprb, 4.05994e-03_jprb, 4.84130e-03_jprb, 5.77305e-03_jprb, 6.88412e-03_jprb, &
     & 8.20902e-03_jprb, 9.78890e-03_jprb, 1.16728e-02_jprb, 1.39194e-02_jprb, 1.65983e-02_jprb, &
     & 1.97927e-02_jprb, 2.36020e-02_jprb, 2.81444e-02_jprb, 3.35610e-02_jprb, 4.00200e-02_jprb, &
     & 4.77222e-02_jprb, 5.69067e-02_jprb, 6.78588e-02_jprb, 8.09187e-02_jprb/)
      kao_mco( 2, :,11) = (/ &
     & 3.85021e-02_jprb, 4.02208e-02_jprb, 4.20162e-02_jprb, 4.38918e-02_jprb, 4.58512e-02_jprb, &
     & 4.78980e-02_jprb, 5.00361e-02_jprb, 5.22697e-02_jprb, 5.46031e-02_jprb, 5.70405e-02_jprb, &
     & 5.95868e-02_jprb, 6.22468e-02_jprb, 6.50254e-02_jprb, 6.79282e-02_jprb, 7.09605e-02_jprb, &
     & 7.41282e-02_jprb, 7.74372e-02_jprb, 8.08940e-02_jprb, 8.45051e-02_jprb/)
      kao_mco( 3, :,11) = (/ &
     & 5.24852e-01_jprb, 5.10480e-01_jprb, 4.96501e-01_jprb, 4.82905e-01_jprb, 4.69681e-01_jprb, &
     & 4.56820e-01_jprb, 4.44310e-01_jprb, 4.32143e-01_jprb, 4.20310e-01_jprb, 4.08800e-01_jprb, &
     & 3.97606e-01_jprb, 3.86718e-01_jprb, 3.76128e-01_jprb, 3.65828e-01_jprb, 3.55810e-01_jprb, &
     & 3.46067e-01_jprb, 3.36590e-01_jprb, 3.27373e-01_jprb, 3.18409e-01_jprb/)
      kao_mco( 4, :,11) = (/ &
     & 5.86290e-01_jprb, 5.70241e-01_jprb, 5.54632e-01_jprb, 5.39450e-01_jprb, 5.24683e-01_jprb, &
     & 5.10321e-01_jprb, 4.96352e-01_jprb, 4.82765e-01_jprb, 4.69550e-01_jprb, 4.56697e-01_jprb, &
     & 4.44196e-01_jprb, 4.32036e-01_jprb, 4.20210e-01_jprb, 4.08708e-01_jprb, 3.97520e-01_jprb, &
     & 3.86638e-01_jprb, 3.76055e-01_jprb, 3.65761e-01_jprb, 3.55749e-01_jprb/)
      kao_mco( 5, :,11) = (/ &
     & 1.66977e+00_jprb, 1.61807e+00_jprb, 1.56798e+00_jprb, 1.51943e+00_jprb, 1.47239e+00_jprb, &
     & 1.42681e+00_jprb, 1.38264e+00_jprb, 1.33983e+00_jprb, 1.29835e+00_jprb, 1.25815e+00_jprb, &
     & 1.21920e+00_jprb, 1.18146e+00_jprb, 1.14488e+00_jprb, 1.10943e+00_jprb, 1.07509e+00_jprb, &
     & 1.04180e+00_jprb, 1.00955e+00_jprb, 9.78295e-01_jprb, 9.48008e-01_jprb/)
      kao_mco( 6, :,11) = (/ &
     & 1.96627e+00_jprb, 1.90948e+00_jprb, 1.85432e+00_jprb, 1.80076e+00_jprb, 1.74875e+00_jprb, &
     & 1.69823e+00_jprb, 1.64918e+00_jprb, 1.60155e+00_jprb, 1.55529e+00_jprb, 1.51036e+00_jprb, &
     & 1.46674e+00_jprb, 1.42437e+00_jprb, 1.38323e+00_jprb, 1.34328e+00_jprb, 1.30448e+00_jprb, &
     & 1.26680e+00_jprb, 1.23021e+00_jprb, 1.19467e+00_jprb, 1.16016e+00_jprb/)
      kao_mco( 7, :,11) = (/ &
     & 1.67574e+00_jprb, 1.63510e+00_jprb, 1.59544e+00_jprb, 1.55674e+00_jprb, 1.51898e+00_jprb, &
     & 1.48213e+00_jprb, 1.44618e+00_jprb, 1.41111e+00_jprb, 1.37688e+00_jprb, 1.34348e+00_jprb, &
     & 1.31090e+00_jprb, 1.27910e+00_jprb, 1.24808e+00_jprb, 1.21780e+00_jprb, 1.18826e+00_jprb, &
     & 1.15944e+00_jprb, 1.13132e+00_jprb, 1.10388e+00_jprb, 1.07710e+00_jprb/)
      kao_mco( 8, :,11) = (/ &
     & 2.00764e+00_jprb, 1.96233e+00_jprb, 1.91803e+00_jprb, 1.87474e+00_jprb, 1.83242e+00_jprb, &
     & 1.79106e+00_jprb, 1.75063e+00_jprb, 1.71111e+00_jprb, 1.67249e+00_jprb, 1.63474e+00_jprb, &
     & 1.59784e+00_jprb, 1.56177e+00_jprb, 1.52652e+00_jprb, 1.49206e+00_jprb, 1.45838e+00_jprb, &
     & 1.42546e+00_jprb, 1.39329e+00_jprb, 1.36184e+00_jprb, 1.33110e+00_jprb/)
      kao_mco( 9, :,11) = (/ &
     & 1.83026e+00_jprb, 1.77349e+00_jprb, 1.71849e+00_jprb, 1.66519e+00_jprb, 1.61355e+00_jprb, &
     & 1.56350e+00_jprb, 1.51501e+00_jprb, 1.46803e+00_jprb, 1.42250e+00_jprb, 1.37838e+00_jprb, &
     & 1.33563e+00_jprb, 1.29421e+00_jprb, 1.25407e+00_jprb, 1.21517e+00_jprb, 1.17749e+00_jprb, &
     & 1.14097e+00_jprb, 1.10558e+00_jprb, 1.07129e+00_jprb, 1.03807e+00_jprb/)
      kao_mco( 1, :,12) = (/ &
     & 3.90309e-04_jprb, 4.81310e-04_jprb, 5.93528e-04_jprb, 7.31909e-04_jprb, 9.02554e-04_jprb, &
     & 1.11299e-03_jprb, 1.37248e-03_jprb, 1.69247e-03_jprb, 2.08708e-03_jprb, 2.57368e-03_jprb, &
     & 3.17374e-03_jprb, 3.91370e-03_jprb, 4.82617e-03_jprb, 5.95140e-03_jprb, 7.33897e-03_jprb, &
     & 9.05007e-03_jprb, 1.11601e-02_jprb, 1.37621e-02_jprb, 1.69707e-02_jprb/)
      kao_mco( 2, :,12) = (/ &
     & 9.80585e-02_jprb, 9.77457e-02_jprb, 9.74339e-02_jprb, 9.71231e-02_jprb, 9.68132e-02_jprb, &
     & 9.65044e-02_jprb, 9.61965e-02_jprb, 9.58897e-02_jprb, 9.55838e-02_jprb, 9.52789e-02_jprb, &
     & 9.49749e-02_jprb, 9.46719e-02_jprb, 9.43699e-02_jprb, 9.40689e-02_jprb, 9.37688e-02_jprb, &
     & 9.34697e-02_jprb, 9.31715e-02_jprb, 9.28743e-02_jprb, 9.25780e-02_jprb/)
      kao_mco( 3, :,12) = (/ &
     & 3.15258e-01_jprb, 3.09936e-01_jprb, 3.04704e-01_jprb, 2.99560e-01_jprb, 2.94503e-01_jprb, &
     & 2.89532e-01_jprb, 2.84645e-01_jprb, 2.79840e-01_jprb, 2.75116e-01_jprb, 2.70472e-01_jprb, &
     & 2.65906e-01_jprb, 2.61417e-01_jprb, 2.57004e-01_jprb, 2.52666e-01_jprb, 2.48401e-01_jprb, &
     & 2.44207e-01_jprb, 2.40085e-01_jprb, 2.36032e-01_jprb, 2.32048e-01_jprb/)
      kao_mco( 4, :,12) = (/ &
     & 9.74407e-01_jprb, 9.46900e-01_jprb, 9.20170e-01_jprb, 8.94195e-01_jprb, 8.68952e-01_jprb, &
     & 8.44422e-01_jprb, 8.20585e-01_jprb, 7.97421e-01_jprb, 7.74910e-01_jprb, 7.53035e-01_jprb, &
     & 7.31777e-01_jprb, 7.11120e-01_jprb, 6.91046e-01_jprb, 6.71538e-01_jprb, 6.52581e-01_jprb, &
     & 6.34159e-01_jprb, 6.16257e-01_jprb, 5.98861e-01_jprb, 5.81956e-01_jprb/)
      kao_mco( 5, :,12) = (/ &
     & 1.04234e+00_jprb, 1.01364e+00_jprb, 9.85726e-01_jprb, 9.58581e-01_jprb, 9.32184e-01_jprb, &
     & 9.06514e-01_jprb, 8.81551e-01_jprb, 8.57275e-01_jprb, 8.33668e-01_jprb, 8.10710e-01_jprb, &
     & 7.88385e-01_jprb, 7.66675e-01_jprb, 7.45563e-01_jprb, 7.25032e-01_jprb, 7.05066e-01_jprb, &
     & 6.85650e-01_jprb, 6.66769e-01_jprb, 6.48408e-01_jprb, 6.30552e-01_jprb/)
      kao_mco( 6, :,12) = (/ &
     & 1.79052e+00_jprb, 1.73725e+00_jprb, 1.68557e+00_jprb, 1.63543e+00_jprb, 1.58678e+00_jprb, &
     & 1.53957e+00_jprb, 1.49377e+00_jprb, 1.44933e+00_jprb, 1.40622e+00_jprb, 1.36439e+00_jprb, &
     & 1.32380e+00_jprb, 1.28442e+00_jprb, 1.24621e+00_jprb, 1.20913e+00_jprb, 1.17316e+00_jprb, &
     & 1.13826e+00_jprb, 1.10440e+00_jprb, 1.07155e+00_jprb, 1.03967e+00_jprb/)
      kao_mco( 7, :,12) = (/ &
     & 2.99551e+00_jprb, 2.90366e+00_jprb, 2.81462e+00_jprb, 2.72831e+00_jprb, 2.64464e+00_jprb, &
     & 2.56355e+00_jprb, 2.48494e+00_jprb, 2.40874e+00_jprb, 2.33487e+00_jprb, 2.26328e+00_jprb, &
     & 2.19387e+00_jprb, 2.12660e+00_jprb, 2.06139e+00_jprb, 1.99818e+00_jprb, 1.93690e+00_jprb, &
     & 1.87751e+00_jprb, 1.81993e+00_jprb, 1.76413e+00_jprb, 1.71003e+00_jprb/)
      kao_mco( 8, :,12) = (/ &
     & 2.89665e+00_jprb, 2.81184e+00_jprb, 2.72951e+00_jprb, 2.64960e+00_jprb, 2.57202e+00_jprb, &
     & 2.49672e+00_jprb, 2.42362e+00_jprb, 2.35266e+00_jprb, 2.28378e+00_jprb, 2.21692e+00_jprb, &
     & 2.15201e+00_jprb, 2.08900e+00_jprb, 2.02784e+00_jprb, 1.96847e+00_jprb, 1.91084e+00_jprb, &
     & 1.85489e+00_jprb, 1.80059e+00_jprb, 1.74787e+00_jprb, 1.69669e+00_jprb/)
      kao_mco( 9, :,12) = (/ &
     & 1.03145e+00_jprb, 1.00335e+00_jprb, 9.76014e-01_jprb, 9.49426e-01_jprb, 9.23561e-01_jprb, &
     & 8.98402e-01_jprb, 8.73927e-01_jprb, 8.50120e-01_jprb, 8.26961e-01_jprb, 8.04433e-01_jprb, &
     & 7.82518e-01_jprb, 7.61201e-01_jprb, 7.40464e-01_jprb, 7.20293e-01_jprb, 7.00670e-01_jprb, &
     & 6.81583e-01_jprb, 6.63015e-01_jprb, 6.44953e-01_jprb, 6.27383e-01_jprb/)
      kao_mco( 1, :,13) = (/ &
     & 5.27769e-04_jprb, 6.65449e-04_jprb, 8.39047e-04_jprb, 1.05793e-03_jprb, 1.33392e-03_jprb, &
     & 1.68190e-03_jprb, 2.12066e-03_jprb, 2.67388e-03_jprb, 3.37143e-03_jprb, 4.25094e-03_jprb, &
     & 5.35990e-03_jprb, 6.75816e-03_jprb, 8.52117e-03_jprb, 1.07441e-02_jprb, 1.35470e-02_jprb, &
     & 1.70810e-02_jprb, 2.15370e-02_jprb, 2.71554e-02_jprb, 3.42395e-02_jprb/)
      kao_mco( 2, :,13) = (/ &
     & 1.08329e-01_jprb, 1.10179e-01_jprb, 1.12060e-01_jprb, 1.13974e-01_jprb, 1.15920e-01_jprb, &
     & 1.17899e-01_jprb, 1.19913e-01_jprb, 1.21960e-01_jprb, 1.24043e-01_jprb, 1.26161e-01_jprb, &
     & 1.28316e-01_jprb, 1.30507e-01_jprb, 1.32735e-01_jprb, 1.35002e-01_jprb, 1.37307e-01_jprb, &
     & 1.39652e-01_jprb, 1.42037e-01_jprb, 1.44462e-01_jprb, 1.46929e-01_jprb/)
      kao_mco( 3, :,13) = (/ &
     & 1.95992e-01_jprb, 1.94515e-01_jprb, 1.93049e-01_jprb, 1.91594e-01_jprb, 1.90150e-01_jprb, &
     & 1.88717e-01_jprb, 1.87294e-01_jprb, 1.85882e-01_jprb, 1.84481e-01_jprb, 1.83091e-01_jprb, &
     & 1.81711e-01_jprb, 1.80341e-01_jprb, 1.78982e-01_jprb, 1.77633e-01_jprb, 1.76294e-01_jprb, &
     & 1.74965e-01_jprb, 1.73646e-01_jprb, 1.72337e-01_jprb, 1.71038e-01_jprb/)
      kao_mco( 4, :,13) = (/ &
     & 4.49766e-01_jprb, 4.42749e-01_jprb, 4.35841e-01_jprb, 4.29042e-01_jprb, 4.22348e-01_jprb, &
     & 4.15759e-01_jprb, 4.09272e-01_jprb, 4.02887e-01_jprb, 3.96601e-01_jprb, 3.90414e-01_jprb, &
     & 3.84323e-01_jprb, 3.78327e-01_jprb, 3.72424e-01_jprb, 3.66614e-01_jprb, 3.60894e-01_jprb, &
     & 3.55264e-01_jprb, 3.49721e-01_jprb, 3.44265e-01_jprb, 3.38894e-01_jprb/)
      kao_mco( 5, :,13) = (/ &
     & 1.07498e+00_jprb, 1.04736e+00_jprb, 1.02045e+00_jprb, 9.94232e-01_jprb, 9.68686e-01_jprb, &
     & 9.43797e-01_jprb, 9.19547e-01_jprb, 8.95920e-01_jprb, 8.72900e-01_jprb, 8.50471e-01_jprb, &
     & 8.28619e-01_jprb, 8.07329e-01_jprb, 7.86585e-01_jprb, 7.66374e-01_jprb, 7.46683e-01_jprb, &
     & 7.27498e-01_jprb, 7.08805e-01_jprb, 6.90593e-01_jprb, 6.72849e-01_jprb/)
      kao_mco( 6, :,13) = (/ &
     & 1.66569e+00_jprb, 1.61740e+00_jprb, 1.57052e+00_jprb, 1.52500e+00_jprb, 1.48080e+00_jprb, &
     & 1.43787e+00_jprb, 1.39620e+00_jprb, 1.35573e+00_jprb, 1.31643e+00_jprb, 1.27827e+00_jprb, &
     & 1.24122e+00_jprb, 1.20524e+00_jprb, 1.17031e+00_jprb, 1.13638e+00_jprb, 1.10344e+00_jprb, &
     & 1.07146e+00_jprb, 1.04040e+00_jprb, 1.01025e+00_jprb, 9.80963e-01_jprb/)
      kao_mco( 7, :,13) = (/ &
     & 1.52948e+00_jprb, 1.48763e+00_jprb, 1.44693e+00_jprb, 1.40735e+00_jprb, 1.36884e+00_jprb, &
     & 1.33139e+00_jprb, 1.29497e+00_jprb, 1.25954e+00_jprb, 1.22508e+00_jprb, 1.19156e+00_jprb, &
     & 1.15896e+00_jprb, 1.12725e+00_jprb, 1.09641e+00_jprb, 1.06641e+00_jprb, 1.03724e+00_jprb, &
     & 1.00886e+00_jprb, 9.81259e-01_jprb, 9.54412e-01_jprb, 9.28301e-01_jprb/)
      kao_mco( 8, :,13) = (/ &
     & 3.81027e+00_jprb, 3.69422e+00_jprb, 3.58170e+00_jprb, 3.47261e+00_jprb, 3.36684e+00_jprb, &
     & 3.26429e+00_jprb, 3.16486e+00_jprb, 3.06847e+00_jprb, 2.97501e+00_jprb, 2.88439e+00_jprb, &
     & 2.79654e+00_jprb, 2.71136e+00_jprb, 2.62878e+00_jprb, 2.54871e+00_jprb, 2.47108e+00_jprb, &
     & 2.39581e+00_jprb, 2.32284e+00_jprb, 2.25209e+00_jprb, 2.18350e+00_jprb/)
      kao_mco( 9, :,13) = (/ &
     & 1.12516e+00_jprb, 1.09531e+00_jprb, 1.06625e+00_jprb, 1.03796e+00_jprb, 1.01042e+00_jprb, &
     & 9.83616e-01_jprb, 9.57520e-01_jprb, 9.32117e-01_jprb, 9.07387e-01_jprb, 8.83314e-01_jprb, &
     & 8.59880e-01_jprb, 8.37067e-01_jprb, 8.14859e-01_jprb, 7.93241e-01_jprb, 7.72196e-01_jprb, &
     & 7.51709e-01_jprb, 7.31766e-01_jprb, 7.12352e-01_jprb, 6.93454e-01_jprb/)
      kao_mco( 1, :,14) = (/ &
     & 4.79165e-04_jprb, 6.26966e-04_jprb, 8.20356e-04_jprb, 1.07340e-03_jprb, 1.40450e-03_jprb, &
     & 1.83772e-03_jprb, 2.40457e-03_jprb, 3.14628e-03_jprb, 4.11676e-03_jprb, 5.38660e-03_jprb, &
     & 7.04813e-03_jprb, 9.22216e-03_jprb, 1.20668e-02_jprb, 1.57889e-02_jprb, 2.06590e-02_jprb, &
     & 2.70314e-02_jprb, 3.53694e-02_jprb, 4.62792e-02_jprb, 6.05542e-02_jprb/)
      kao_mco( 2, :,14) = (/ &
     & 4.75367e-04_jprb, 6.21987e-04_jprb, 8.13830e-04_jprb, 1.06484e-03_jprb, 1.39328e-03_jprb, &
     & 1.82302e-03_jprb, 2.38530e-03_jprb, 3.12101e-03_jprb, 4.08365e-03_jprb, 5.34318e-03_jprb, &
     & 6.99121e-03_jprb, 9.14756e-03_jprb, 1.19690e-02_jprb, 1.56607e-02_jprb, 2.04909e-02_jprb, &
     & 2.68111e-02_jprb, 3.50806e-02_jprb, 4.59007e-02_jprb, 6.00580e-02_jprb/)
      kao_mco( 3, :,14) = (/ &
     & 1.39594e-01_jprb, 1.42407e-01_jprb, 1.45276e-01_jprb, 1.48202e-01_jprb, 1.51188e-01_jprb, &
     & 1.54234e-01_jprb, 1.57342e-01_jprb, 1.60512e-01_jprb, 1.63745e-01_jprb, 1.67044e-01_jprb, &
     & 1.70410e-01_jprb, 1.73843e-01_jprb, 1.77345e-01_jprb, 1.80918e-01_jprb, 1.84563e-01_jprb, &
     & 1.88282e-01_jprb, 1.92075e-01_jprb, 1.95945e-01_jprb, 1.99892e-01_jprb/)
      kao_mco( 4, :,14) = (/ &
     & 2.73418e-01_jprb, 2.71322e-01_jprb, 2.69242e-01_jprb, 2.67177e-01_jprb, 2.65129e-01_jprb, &
     & 2.63096e-01_jprb, 2.61079e-01_jprb, 2.59077e-01_jprb, 2.57091e-01_jprb, 2.55120e-01_jprb, &
     & 2.53164e-01_jprb, 2.51223e-01_jprb, 2.49297e-01_jprb, 2.47386e-01_jprb, 2.45489e-01_jprb, &
     & 2.43607e-01_jprb, 2.41739e-01_jprb, 2.39886e-01_jprb, 2.38047e-01_jprb/)
      kao_mco( 5, :,14) = (/ &
     & 5.01880e-01_jprb, 4.95066e-01_jprb, 4.88344e-01_jprb, 4.81713e-01_jprb, 4.75173e-01_jprb, &
     & 4.68721e-01_jprb, 4.62357e-01_jprb, 4.56079e-01_jprb, 4.49887e-01_jprb, 4.43778e-01_jprb, &
     & 4.37753e-01_jprb, 4.31809e-01_jprb, 4.25946e-01_jprb, 4.20163e-01_jprb, 4.14458e-01_jprb, &
     & 4.08830e-01_jprb, 4.03279e-01_jprb, 3.97804e-01_jprb, 3.92403e-01_jprb/)
      kao_mco( 6, :,14) = (/ &
     & 9.46125e-01_jprb, 9.29642e-01_jprb, 9.13447e-01_jprb, 8.97533e-01_jprb, 8.81897e-01_jprb, &
     & 8.66533e-01_jprb, 8.51437e-01_jprb, 8.36603e-01_jprb, 8.22029e-01_jprb, 8.07708e-01_jprb, &
     & 7.93636e-01_jprb, 7.79810e-01_jprb, 7.66225e-01_jprb, 7.52876e-01_jprb, 7.39760e-01_jprb, &
     & 7.26872e-01_jprb, 7.14209e-01_jprb, 7.01766e-01_jprb, 6.89540e-01_jprb/)
      kao_mco( 7, :,14) = (/ &
     & 2.47697e+00_jprb, 2.41183e+00_jprb, 2.34840e+00_jprb, 2.28664e+00_jprb, 2.22650e+00_jprb, &
     & 2.16795e+00_jprb, 2.11093e+00_jprb, 2.05541e+00_jprb, 2.00136e+00_jprb, 1.94872e+00_jprb, &
     & 1.89747e+00_jprb, 1.84757e+00_jprb, 1.79898e+00_jprb, 1.75167e+00_jprb, 1.70560e+00_jprb, &
     & 1.66074e+00_jprb, 1.61707e+00_jprb, 1.57454e+00_jprb, 1.53313e+00_jprb/)
      kao_mco( 8, :,14) = (/ &
     & 2.19323e+00_jprb, 2.13926e+00_jprb, 2.08662e+00_jprb, 2.03528e+00_jprb, 1.98519e+00_jprb, &
     & 1.93635e+00_jprb, 1.88870e+00_jprb, 1.84222e+00_jprb, 1.79689e+00_jprb, 1.75268e+00_jprb, &
     & 1.70955e+00_jprb, 1.66748e+00_jprb, 1.62645e+00_jprb, 1.58643e+00_jprb, 1.54739e+00_jprb, &
     & 1.50932e+00_jprb, 1.47218e+00_jprb, 1.43595e+00_jprb, 1.40062e+00_jprb/)
      kao_mco( 9, :,14) = (/ &
     & 5.68351e-01_jprb, 5.60360e-01_jprb, 5.52481e-01_jprb, 5.44712e-01_jprb, 5.37053e-01_jprb, &
     & 5.29502e-01_jprb, 5.22057e-01_jprb, 5.14716e-01_jprb, 5.07479e-01_jprb, 5.00343e-01_jprb, &
     & 4.93308e-01_jprb, 4.86372e-01_jprb, 4.79533e-01_jprb, 4.72790e-01_jprb, 4.66142e-01_jprb, &
     & 4.59588e-01_jprb, 4.53126e-01_jprb, 4.46754e-01_jprb, 4.40473e-01_jprb/)
      kao_mco( 1, :,15) = (/ &
     & 7.87937e-04_jprb, 1.01733e-03_jprb, 1.31351e-03_jprb, 1.69591e-03_jprb, 2.18964e-03_jprb, &
     & 2.82712e-03_jprb, 3.65018e-03_jprb, 4.71286e-03_jprb, 6.08493e-03_jprb, 7.85644e-03_jprb, &
     & 1.01437e-02_jprb, 1.30969e-02_jprb, 1.69098e-02_jprb, 2.18327e-02_jprb, 2.81889e-02_jprb, &
     & 3.63956e-02_jprb, 4.69915e-02_jprb, 6.06722e-02_jprb, 7.83357e-02_jprb/)
      kao_mco( 2, :,15) = (/ &
     & 7.87937e-04_jprb, 1.01733e-03_jprb, 1.31351e-03_jprb, 1.69591e-03_jprb, 2.18964e-03_jprb, &
     & 2.82712e-03_jprb, 3.65018e-03_jprb, 4.71286e-03_jprb, 6.08493e-03_jprb, 7.85644e-03_jprb, &
     & 1.01437e-02_jprb, 1.30969e-02_jprb, 1.69098e-02_jprb, 2.18327e-02_jprb, 2.81889e-02_jprb, &
     & 3.63956e-02_jprb, 4.69915e-02_jprb, 6.06722e-02_jprb, 7.83357e-02_jprb/)
      kao_mco( 3, :,15) = (/ &
     & 1.97281e-01_jprb, 2.00714e-01_jprb, 2.04206e-01_jprb, 2.07760e-01_jprb, 2.11375e-01_jprb, &
     & 2.15053e-01_jprb, 2.18795e-01_jprb, 2.22603e-01_jprb, 2.26476e-01_jprb, 2.30417e-01_jprb, &
     & 2.34426e-01_jprb, 2.38506e-01_jprb, 2.42656e-01_jprb, 2.46878e-01_jprb, 2.51174e-01_jprb, &
     & 2.55545e-01_jprb, 2.59992e-01_jprb, 2.64516e-01_jprb, 2.69119e-01_jprb/)
      kao_mco( 4, :,15) = (/ &
     & 4.47509e-01_jprb, 4.52222e-01_jprb, 4.56985e-01_jprb, 4.61799e-01_jprb, 4.66663e-01_jprb, &
     & 4.71578e-01_jprb, 4.76545e-01_jprb, 4.81565e-01_jprb, 4.86637e-01_jprb, 4.91763e-01_jprb, &
     & 4.96942e-01_jprb, 5.02177e-01_jprb, 5.07466e-01_jprb, 5.12811e-01_jprb, 5.18213e-01_jprb, &
     & 5.23671e-01_jprb, 5.29187e-01_jprb, 5.34761e-01_jprb, 5.40393e-01_jprb/)
      kao_mco( 5, :,15) = (/ &
     & 1.02732e+00_jprb, 1.02091e+00_jprb, 1.01453e+00_jprb, 1.00820e+00_jprb, 1.00191e+00_jprb, &
     & 9.95660e-01_jprb, 9.89447e-01_jprb, 9.83273e-01_jprb, 9.77137e-01_jprb, 9.71039e-01_jprb, &
     & 9.64980e-01_jprb, 9.58958e-01_jprb, 9.52974e-01_jprb, 9.47027e-01_jprb, 9.41118e-01_jprb, &
     & 9.35245e-01_jprb, 9.29409e-01_jprb, 9.23609e-01_jprb, 9.17846e-01_jprb/)
      kao_mco( 6, :,15) = (/ &
     & 1.13766e+00_jprb, 1.12944e+00_jprb, 1.12128e+00_jprb, 1.11318e+00_jprb, 1.10514e+00_jprb, &
     & 1.09715e+00_jprb, 1.08923e+00_jprb, 1.08136e+00_jprb, 1.07355e+00_jprb, 1.06579e+00_jprb, &
     & 1.05809e+00_jprb, 1.05045e+00_jprb, 1.04286e+00_jprb, 1.03532e+00_jprb, 1.02784e+00_jprb, &
     & 1.02042e+00_jprb, 1.01305e+00_jprb, 1.00573e+00_jprb, 9.98462e-01_jprb/)
      kao_mco( 7, :,15) = (/ &
     & 1.13268e+00_jprb, 1.12496e+00_jprb, 1.11730e+00_jprb, 1.10969e+00_jprb, 1.10213e+00_jprb, &
     & 1.09462e+00_jprb, 1.08717e+00_jprb, 1.07976e+00_jprb, 1.07241e+00_jprb, 1.06510e+00_jprb, &
     & 1.05785e+00_jprb, 1.05064e+00_jprb, 1.04349e+00_jprb, 1.03638e+00_jprb, 1.02932e+00_jprb, &
     & 1.02231e+00_jprb, 1.01535e+00_jprb, 1.00843e+00_jprb, 1.00156e+00_jprb/)
      kao_mco( 8, :,15) = (/ &
     & 1.11982e+00_jprb, 1.11343e+00_jprb, 1.10707e+00_jprb, 1.10075e+00_jprb, 1.09447e+00_jprb, &
     & 1.08822e+00_jprb, 1.08201e+00_jprb, 1.07583e+00_jprb, 1.06969e+00_jprb, 1.06358e+00_jprb, &
     & 1.05751e+00_jprb, 1.05147e+00_jprb, 1.04547e+00_jprb, 1.03950e+00_jprb, 1.03356e+00_jprb, &
     & 1.02766e+00_jprb, 1.02180e+00_jprb, 1.01596e+00_jprb, 1.01016e+00_jprb/)
      kao_mco( 9, :,15) = (/ &
     & 1.03561e+00_jprb, 1.02902e+00_jprb, 1.02246e+00_jprb, 1.01595e+00_jprb, 1.00948e+00_jprb, &
     & 1.00305e+00_jprb, 9.96667e-01_jprb, 9.90320e-01_jprb, 9.84013e-01_jprb, 9.77747e-01_jprb, &
     & 9.71520e-01_jprb, 9.65334e-01_jprb, 9.59186e-01_jprb, 9.53078e-01_jprb, 9.47008e-01_jprb, &
     & 9.40978e-01_jprb, 9.34985e-01_jprb, 9.29031e-01_jprb, 9.23115e-01_jprb/)
      kao_mco( 1, :,16) = (/ &
     & 1.22217e-03_jprb, 1.56836e-03_jprb, 2.01261e-03_jprb, 2.58271e-03_jprb, 3.31429e-03_jprb, &
     & 4.25310e-03_jprb, 5.45784e-03_jprb, 7.00383e-03_jprb, 8.98775e-03_jprb, 1.15336e-02_jprb, &
     & 1.48007e-02_jprb, 1.89931e-02_jprb, 2.43731e-02_jprb, 3.12771e-02_jprb, 4.01367e-02_jprb, &
     & 5.15059e-02_jprb, 6.60956e-02_jprb, 8.48179e-02_jprb, 1.08843e-01_jprb/)
      kao_mco( 2, :,16) = (/ &
     & 1.22217e-03_jprb, 1.56836e-03_jprb, 2.01261e-03_jprb, 2.58271e-03_jprb, 3.31429e-03_jprb, &
     & 4.25310e-03_jprb, 5.45784e-03_jprb, 7.00383e-03_jprb, 8.98775e-03_jprb, 1.15336e-02_jprb, &
     & 1.48007e-02_jprb, 1.89931e-02_jprb, 2.43731e-02_jprb, 3.12771e-02_jprb, 4.01367e-02_jprb, &
     & 5.15059e-02_jprb, 6.60956e-02_jprb, 8.48179e-02_jprb, 1.08843e-01_jprb/)
      kao_mco( 3, :,16) = (/ &
     & 1.22217e-03_jprb, 1.56836e-03_jprb, 2.01261e-03_jprb, 2.58271e-03_jprb, 3.31429e-03_jprb, &
     & 4.25310e-03_jprb, 5.45784e-03_jprb, 7.00383e-03_jprb, 8.98775e-03_jprb, 1.15336e-02_jprb, &
     & 1.48007e-02_jprb, 1.89931e-02_jprb, 2.43731e-02_jprb, 3.12771e-02_jprb, 4.01367e-02_jprb, &
     & 5.15059e-02_jprb, 6.60956e-02_jprb, 8.48179e-02_jprb, 1.08843e-01_jprb/)
      kao_mco( 4, :,16) = (/ &
     & 1.01221e+00_jprb, 1.01660e+00_jprb, 1.02101e+00_jprb, 1.02544e+00_jprb, 1.02989e+00_jprb, &
     & 1.03436e+00_jprb, 1.03884e+00_jprb, 1.04335e+00_jprb, 1.04788e+00_jprb, 1.05243e+00_jprb, &
     & 1.05699e+00_jprb, 1.06158e+00_jprb, 1.06619e+00_jprb, 1.07081e+00_jprb, 1.07546e+00_jprb, &
     & 1.08012e+00_jprb, 1.08481e+00_jprb, 1.08952e+00_jprb, 1.09425e+00_jprb/)
      kao_mco( 5, :,16) = (/ &
     & 1.01221e+00_jprb, 1.01660e+00_jprb, 1.02101e+00_jprb, 1.02544e+00_jprb, 1.02989e+00_jprb, &
     & 1.03436e+00_jprb, 1.03884e+00_jprb, 1.04335e+00_jprb, 1.04788e+00_jprb, 1.05243e+00_jprb, &
     & 1.05699e+00_jprb, 1.06158e+00_jprb, 1.06619e+00_jprb, 1.07081e+00_jprb, 1.07546e+00_jprb, &
     & 1.08012e+00_jprb, 1.08481e+00_jprb, 1.08952e+00_jprb, 1.09425e+00_jprb/)
      kao_mco( 6, :,16) = (/ &
     & 1.01221e+00_jprb, 1.01660e+00_jprb, 1.02101e+00_jprb, 1.02544e+00_jprb, 1.02989e+00_jprb, &
     & 1.03436e+00_jprb, 1.03884e+00_jprb, 1.04335e+00_jprb, 1.04788e+00_jprb, 1.05243e+00_jprb, &
     & 1.05699e+00_jprb, 1.06158e+00_jprb, 1.06619e+00_jprb, 1.07081e+00_jprb, 1.07546e+00_jprb, &
     & 1.08012e+00_jprb, 1.08481e+00_jprb, 1.08952e+00_jprb, 1.09425e+00_jprb/)
      kao_mco( 7, :,16) = (/ &
     & 1.01221e+00_jprb, 1.01660e+00_jprb, 1.02101e+00_jprb, 1.02544e+00_jprb, 1.02989e+00_jprb, &
     & 1.03436e+00_jprb, 1.03884e+00_jprb, 1.04335e+00_jprb, 1.04788e+00_jprb, 1.05243e+00_jprb, &
     & 1.05699e+00_jprb, 1.06158e+00_jprb, 1.06619e+00_jprb, 1.07081e+00_jprb, 1.07546e+00_jprb, &
     & 1.08012e+00_jprb, 1.08481e+00_jprb, 1.08952e+00_jprb, 1.09425e+00_jprb/)
      kao_mco( 8, :,16) = (/ &
     & 1.01221e+00_jprb, 1.01660e+00_jprb, 1.02101e+00_jprb, 1.02544e+00_jprb, 1.02989e+00_jprb, &
     & 1.03436e+00_jprb, 1.03884e+00_jprb, 1.04335e+00_jprb, 1.04788e+00_jprb, 1.05243e+00_jprb, &
     & 1.05699e+00_jprb, 1.06158e+00_jprb, 1.06619e+00_jprb, 1.07081e+00_jprb, 1.07546e+00_jprb, &
     & 1.08012e+00_jprb, 1.08481e+00_jprb, 1.08952e+00_jprb, 1.09425e+00_jprb/)
      kao_mco( 9, :,16) = (/ &
     & 1.01221e+00_jprb, 1.01660e+00_jprb, 1.02101e+00_jprb, 1.02544e+00_jprb, 1.02989e+00_jprb, &
     & 1.03436e+00_jprb, 1.03884e+00_jprb, 1.04335e+00_jprb, 1.04788e+00_jprb, 1.05243e+00_jprb, &
     & 1.05699e+00_jprb, 1.06158e+00_jprb, 1.06619e+00_jprb, 1.07081e+00_jprb, 1.07546e+00_jprb, &
     & 1.08012e+00_jprb, 1.08481e+00_jprb, 1.08952e+00_jprb, 1.09425e+00_jprb/)

!     the array kbo_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kbo_mo3(:, 1) = (/ &
     & 1.07596e-02_jprb, 1.12146e-02_jprb, 1.16887e-02_jprb, 1.21830e-02_jprb, 1.26981e-02_jprb, &
     & 1.32350e-02_jprb, 1.37946e-02_jprb, 1.43779e-02_jprb, 1.49858e-02_jprb, 1.56194e-02_jprb, &
     & 1.62799e-02_jprb, 1.69682e-02_jprb, 1.76857e-02_jprb, 1.84334e-02_jprb, 1.92129e-02_jprb, &
     & 2.00252e-02_jprb, 2.08719e-02_jprb, 2.17544e-02_jprb, 2.26743e-02_jprb/)
      kbo_mo3(:, 2) = (/ &
     & 9.48276e-02_jprb, 9.66591e-02_jprb, 9.85260e-02_jprb, 1.00429e-01_jprb, 1.02369e-01_jprb, &
     & 1.04346e-01_jprb, 1.06361e-01_jprb, 1.08416e-01_jprb, 1.10510e-01_jprb, 1.12644e-01_jprb, &
     & 1.14820e-01_jprb, 1.17037e-01_jprb, 1.19298e-01_jprb, 1.21602e-01_jprb, 1.23951e-01_jprb, &
     & 1.26345e-01_jprb, 1.28785e-01_jprb, 1.31273e-01_jprb, 1.33808e-01_jprb/)
      kbo_mo3(:, 3) = (/ &
     & 3.54721e-01_jprb, 3.55779e-01_jprb, 3.56841e-01_jprb, 3.57906e-01_jprb, 3.58973e-01_jprb, &
     & 3.60044e-01_jprb, 3.61119e-01_jprb, 3.62196e-01_jprb, 3.63277e-01_jprb, 3.64360e-01_jprb, &
     & 3.65448e-01_jprb, 3.66538e-01_jprb, 3.67631e-01_jprb, 3.68728e-01_jprb, 3.69828e-01_jprb, &
     & 3.70932e-01_jprb, 3.72038e-01_jprb, 3.73148e-01_jprb, 3.74262e-01_jprb/)
      kbo_mo3(:, 4) = (/ &
     & 6.46454e-01_jprb, 6.43823e-01_jprb, 6.41202e-01_jprb, 6.38593e-01_jprb, 6.35994e-01_jprb, &
     & 6.33405e-01_jprb, 6.30827e-01_jprb, 6.28260e-01_jprb, 6.25703e-01_jprb, 6.23156e-01_jprb, &
     & 6.20620e-01_jprb, 6.18094e-01_jprb, 6.15578e-01_jprb, 6.13073e-01_jprb, 6.10578e-01_jprb, &
     & 6.08093e-01_jprb, 6.05618e-01_jprb, 6.03153e-01_jprb, 6.00698e-01_jprb/)
      kbo_mo3(:, 5) = (/ &
     & 9.29832e-01_jprb, 9.22877e-01_jprb, 9.15975e-01_jprb, 9.09124e-01_jprb, 9.02324e-01_jprb, &
     & 8.95576e-01_jprb, 8.88877e-01_jprb, 8.82229e-01_jprb, 8.75631e-01_jprb, 8.69082e-01_jprb, &
     & 8.62582e-01_jprb, 8.56130e-01_jprb, 8.49727e-01_jprb, 8.43372e-01_jprb, 8.37064e-01_jprb, &
     & 8.30803e-01_jprb, 8.24589e-01_jprb, 8.18422e-01_jprb, 8.12301e-01_jprb/)
      kbo_mo3(:, 6) = (/ &
     & 1.43531e+00_jprb, 1.42616e+00_jprb, 1.41706e+00_jprb, 1.40802e+00_jprb, 1.39903e+00_jprb, &
     & 1.39010e+00_jprb, 1.38124e+00_jprb, 1.37242e+00_jprb, 1.36367e+00_jprb, 1.35496e+00_jprb, &
     & 1.34632e+00_jprb, 1.33773e+00_jprb, 1.32919e+00_jprb, 1.32071e+00_jprb, 1.31229e+00_jprb, &
     & 1.30391e+00_jprb, 1.29559e+00_jprb, 1.28733e+00_jprb, 1.27911e+00_jprb/)
      kbo_mo3(:, 7) = (/ &
     & 2.68664e+00_jprb, 2.67196e+00_jprb, 2.65736e+00_jprb, 2.64284e+00_jprb, 2.62840e+00_jprb, &
     & 2.61404e+00_jprb, 2.59975e+00_jprb, 2.58555e+00_jprb, 2.57142e+00_jprb, 2.55737e+00_jprb, &
     & 2.54340e+00_jprb, 2.52950e+00_jprb, 2.51568e+00_jprb, 2.50193e+00_jprb, 2.48826e+00_jprb, &
     & 2.47466e+00_jprb, 2.46114e+00_jprb, 2.44769e+00_jprb, 2.43432e+00_jprb/)
      kbo_mo3(:, 8) = (/ &
     & 2.45343e+00_jprb, 2.43442e+00_jprb, 2.41556e+00_jprb, 2.39684e+00_jprb, 2.37827e+00_jprb, &
     & 2.35984e+00_jprb, 2.34156e+00_jprb, 2.32342e+00_jprb, 2.30541e+00_jprb, 2.28755e+00_jprb, &
     & 2.26983e+00_jprb, 2.25224e+00_jprb, 2.23479e+00_jprb, 2.21747e+00_jprb, 2.20029e+00_jprb, &
     & 2.18324e+00_jprb, 2.16633e+00_jprb, 2.14954e+00_jprb, 2.13289e+00_jprb/)
      kbo_mo3(:, 9) = (/ &
     & 1.55879e-01_jprb, 1.55998e-01_jprb, 1.56118e-01_jprb, 1.56238e-01_jprb, 1.56358e-01_jprb, &
     & 1.56478e-01_jprb, 1.56599e-01_jprb, 1.56719e-01_jprb, 1.56840e-01_jprb, 1.56960e-01_jprb, &
     & 1.57081e-01_jprb, 1.57201e-01_jprb, 1.57322e-01_jprb, 1.57443e-01_jprb, 1.57564e-01_jprb, &
     & 1.57685e-01_jprb, 1.57806e-01_jprb, 1.57928e-01_jprb, 1.58049e-01_jprb/)
      kbo_mo3(:,10) = (/ &
     & 8.75149e-03_jprb, 8.88794e-03_jprb, 9.02651e-03_jprb, 9.16725e-03_jprb, 9.31018e-03_jprb, &
     & 9.45534e-03_jprb, 9.60276e-03_jprb, 9.75248e-03_jprb, 9.90454e-03_jprb, 1.00590e-02_jprb, &
     & 1.02158e-02_jprb, 1.03751e-02_jprb, 1.05368e-02_jprb, 1.07011e-02_jprb, 1.08680e-02_jprb, &
     & 1.10374e-02_jprb, 1.12095e-02_jprb, 1.13843e-02_jprb, 1.15618e-02_jprb/)
      kbo_mo3(:,11) = (/ &
     & 8.83874e-03_jprb, 8.97926e-03_jprb, 9.12201e-03_jprb, 9.26703e-03_jprb, 9.41436e-03_jprb, &
     & 9.56403e-03_jprb, 9.71608e-03_jprb, 9.87055e-03_jprb, 1.00275e-02_jprb, 1.01869e-02_jprb, &
     & 1.03488e-02_jprb, 1.05134e-02_jprb, 1.06805e-02_jprb, 1.08503e-02_jprb, 1.10228e-02_jprb, &
     & 1.11980e-02_jprb, 1.13761e-02_jprb, 1.15569e-02_jprb, 1.17407e-02_jprb/)
      kbo_mo3(:,12) = (/ &
     & 9.59461e-03_jprb, 9.70417e-03_jprb, 9.81498e-03_jprb, 9.92705e-03_jprb, 1.00404e-02_jprb, &
     & 1.01550e-02_jprb, 1.02710e-02_jprb, 1.03883e-02_jprb, 1.05069e-02_jprb, 1.06269e-02_jprb, &
     & 1.07482e-02_jprb, 1.08709e-02_jprb, 1.09951e-02_jprb, 1.11206e-02_jprb, 1.12476e-02_jprb, &
     & 1.13760e-02_jprb, 1.15059e-02_jprb, 1.16373e-02_jprb, 1.17702e-02_jprb/)
      kbo_mo3(:,13) = (/ &
     & 1.13077e-02_jprb, 1.14079e-02_jprb, 1.15089e-02_jprb, 1.16109e-02_jprb, 1.17138e-02_jprb, &
     & 1.18176e-02_jprb, 1.19223e-02_jprb, 1.20279e-02_jprb, 1.21344e-02_jprb, 1.22419e-02_jprb, &
     & 1.23504e-02_jprb, 1.24598e-02_jprb, 1.25702e-02_jprb, 1.26816e-02_jprb, 1.27939e-02_jprb, &
     & 1.29073e-02_jprb, 1.30216e-02_jprb, 1.31370e-02_jprb, 1.32534e-02_jprb/)
      kbo_mo3(:,14) = (/ &
     & 6.74844e-03_jprb, 6.82637e-03_jprb, 6.90519e-03_jprb, 6.98493e-03_jprb, 7.06558e-03_jprb, &
     & 7.14717e-03_jprb, 7.22970e-03_jprb, 7.31318e-03_jprb, 7.39762e-03_jprb, 7.48304e-03_jprb, &
     & 7.56945e-03_jprb, 7.65686e-03_jprb, 7.74527e-03_jprb, 7.83470e-03_jprb, 7.92517e-03_jprb, &
     & 8.01668e-03_jprb, 8.10925e-03_jprb, 8.20289e-03_jprb, 8.29761e-03_jprb/)
      kbo_mo3(:,15) = (/ &
     & 7.94595e-03_jprb, 8.00015e-03_jprb, 8.05472e-03_jprb, 8.10966e-03_jprb, 8.16497e-03_jprb, &
     & 8.22067e-03_jprb, 8.27674e-03_jprb, 8.33320e-03_jprb, 8.39004e-03_jprb, 8.44727e-03_jprb, &
     & 8.50489e-03_jprb, 8.56290e-03_jprb, 8.62130e-03_jprb, 8.68011e-03_jprb, 8.73932e-03_jprb, &
     & 8.79893e-03_jprb, 8.85895e-03_jprb, 8.91937e-03_jprb, 8.98021e-03_jprb/)
      kbo_mo3(:,16) = (/ &
     & 1.85967e-03_jprb, 1.86082e-03_jprb, 1.86197e-03_jprb, 1.86312e-03_jprb, 1.86428e-03_jprb, &
     & 1.86543e-03_jprb, 1.86658e-03_jprb, 1.86774e-03_jprb, 1.86889e-03_jprb, 1.87005e-03_jprb, &
     & 1.87121e-03_jprb, 1.87236e-03_jprb, 1.87352e-03_jprb, 1.87468e-03_jprb, 1.87584e-03_jprb, &
     & 1.87700e-03_jprb, 1.87816e-03_jprb, 1.87932e-03_jprb, 1.88049e-03_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &1.6586e-05_jprb,1.9995e-05_jprb,1.8582e-05_jprb,1.3988e-05_jprb,1.3650e-05_jprb,1.1079e-05_jprb, &
     &9.5855e-06_jprb,8.4062e-06_jprb,1.3558e-05_jprb,1.8620e-05_jprb,2.2652e-05_jprb,1.7883e-05_jprb, &
     &2.6241e-05_jprb,3.1171e-05_jprb,3.9386e-05_jprb,4.4415e-05_jprb/)
      forrefo(2,:) = (/ &
     &2.0730e-05_jprb,2.3258e-05_jprb,2.1543e-05_jprb,1.5660e-05_jprb,9.7872e-06_jprb,8.1078e-06_jprb, &
     &7.0246e-06_jprb,6.0428e-06_jprb,4.8793e-06_jprb,4.4937e-06_jprb,4.7078e-06_jprb,4.6898e-06_jprb, &
     &6.9481e-06_jprb,8.6269e-06_jprb,3.1761e-06_jprb,3.1440e-06_jprb/)
      forrefo(3,:) = (/ &
     &1.5737e-05_jprb,2.2501e-05_jprb,2.3520e-05_jprb,2.0288e-05_jprb,1.2083e-05_jprb,6.8256e-06_jprb, &
     &6.0637e-06_jprb,5.5434e-06_jprb,4.3888e-06_jprb,3.8435e-06_jprb,3.8477e-06_jprb,3.8314e-06_jprb, &
     &3.8251e-06_jprb,3.3637e-06_jprb,3.1950e-06_jprb,3.1440e-06_jprb/)
      forrefo(4,:) = (/ &
     &1.1400e-05_jprb,7.9751e-06_jprb,8.8659e-06_jprb,1.5884e-05_jprb,1.9118e-05_jprb,1.9429e-05_jprb, &
     &2.0532e-05_jprb,2.2155e-05_jprb,2.3894e-05_jprb,2.2984e-05_jprb,2.3731e-05_jprb,2.4538e-05_jprb, &
     &2.6697e-05_jprb,1.9329e-05_jprb,3.3306e-06_jprb,3.2018e-06_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

      selfrefo(:, 1) = (/ &
     & 9.62275e-03_jprb, 8.29909e-03_jprb, 7.15750e-03_jprb, 6.17294e-03_jprb, 5.32382e-03_jprb, &
     & 4.59150e-03_jprb, 3.95991e-03_jprb, 3.41520e-03_jprb, 2.94542e-03_jprb, 2.54026e-03_jprb/)
      selfrefo(:, 2) = (/ &
     & 9.76664e-03_jprb, 8.47783e-03_jprb, 7.35910e-03_jprb, 6.38799e-03_jprb, 5.54504e-03_jprb, &
     & 4.81331e-03_jprb, 4.17815e-03_jprb, 3.62680e-03_jprb, 3.14821e-03_jprb, 2.73277e-03_jprb/)
      selfrefo(:, 3) = (/ &
     & 9.53856e-03_jprb, 8.23750e-03_jprb, 7.11390e-03_jprb, 6.14356e-03_jprb, 5.30558e-03_jprb, &
     & 4.58190e-03_jprb, 3.95693e-03_jprb, 3.41720e-03_jprb, 2.95109e-03_jprb, 2.54856e-03_jprb/)
      selfrefo(:, 4) = (/ &
     & 8.47621e-03_jprb, 7.29518e-03_jprb, 6.27870e-03_jprb, 5.40385e-03_jprb, 4.65091e-03_jprb, &
     & 4.00287e-03_jprb, 3.44513e-03_jprb, 2.96510e-03_jprb, 2.55196e-03_jprb, 2.19638e-03_jprb/)
      selfrefo(:, 5) = (/ &
     & 6.71258e-03_jprb, 5.95346e-03_jprb, 5.28020e-03_jprb, 4.68307e-03_jprb, 4.15348e-03_jprb, &
     & 3.68377e-03_jprb, 3.26718e-03_jprb, 2.89770e-03_jprb, 2.57000e-03_jprb, 2.27937e-03_jprb/)
      selfrefo(:, 6) = (/ &
     & 6.29140e-03_jprb, 5.55557e-03_jprb, 4.90580e-03_jprb, 4.33203e-03_jprb, 3.82536e-03_jprb, &
     & 3.37795e-03_jprb, 2.98287e-03_jprb, 2.63400e-03_jprb, 2.32593e-03_jprb, 2.05389e-03_jprb/)
      selfrefo(:, 7) = (/ &
     & 6.00229e-03_jprb, 5.28180e-03_jprb, 4.64780e-03_jprb, 4.08990e-03_jprb, 3.59897e-03_jprb, &
     & 3.16696e-03_jprb, 2.78682e-03_jprb, 2.45230e-03_jprb, 2.15794e-03_jprb, 1.89891e-03_jprb/)
      selfrefo(:, 8) = (/ &
     & 5.78892e-03_jprb, 5.07191e-03_jprb, 4.44370e-03_jprb, 3.89330e-03_jprb, 3.41108e-03_jprb, &
     & 2.98858e-03_jprb, 2.61842e-03_jprb, 2.29410e-03_jprb, 2.00995e-03_jprb, 1.76100e-03_jprb/)
      selfrefo(:, 9) = (/ &
     & 4.96186e-03_jprb, 4.56767e-03_jprb, 4.20480e-03_jprb, 3.87076e-03_jprb, 3.56325e-03_jprb, &
     & 3.28017e-03_jprb, 3.01959e-03_jprb, 2.77970e-03_jprb, 2.55887e-03_jprb, 2.35559e-03_jprb/)
      selfrefo(:,10) = (/ &
     & 4.56849e-03_jprb, 4.35527e-03_jprb, 4.15200e-03_jprb, 3.95822e-03_jprb, 3.77348e-03_jprb, &
     & 3.59736e-03_jprb, 3.42946e-03_jprb, 3.26940e-03_jprb, 3.11681e-03_jprb, 2.97134e-03_jprb/)
      selfrefo(:,11) = (/ &
     & 4.47310e-03_jprb, 4.32453e-03_jprb, 4.18090e-03_jprb, 4.04204e-03_jprb, 3.90779e-03_jprb, &
     & 3.77799e-03_jprb, 3.65251e-03_jprb, 3.53120e-03_jprb, 3.41392e-03_jprb, 3.30053e-03_jprb/)
      selfrefo(:,12) = (/ &
     & 4.46459e-03_jprb, 4.24031e-03_jprb, 4.02730e-03_jprb, 3.82499e-03_jprb, 3.63284e-03_jprb, &
     & 3.45035e-03_jprb, 3.27702e-03_jprb, 3.11240e-03_jprb, 2.95605e-03_jprb, 2.80755e-03_jprb/)
      selfrefo(:,13) = (/ &
     & 4.43961e-03_jprb, 4.35658e-03_jprb, 4.27510e-03_jprb, 4.19514e-03_jprb, 4.11669e-03_jprb, &
     & 4.03969e-03_jprb, 3.96414e-03_jprb, 3.89000e-03_jprb, 3.81725e-03_jprb, 3.74585e-03_jprb/)
      selfrefo(:,14) = (/ &
     & 4.40512e-03_jprb, 4.41515e-03_jprb, 4.42520e-03_jprb, 4.43527e-03_jprb, 4.44537e-03_jprb, &
     & 4.45549e-03_jprb, 4.46563e-03_jprb, 4.47580e-03_jprb, 4.48599e-03_jprb, 4.49620e-03_jprb/)
      selfrefo(:,15) = (/ &
     & 3.21965e-03_jprb, 3.42479e-03_jprb, 3.64300e-03_jprb, 3.87512e-03_jprb, 4.12202e-03_jprb, &
     & 4.38466e-03_jprb, 4.66403e-03_jprb, 4.96120e-03_jprb, 5.27731e-03_jprb, 5.61355e-03_jprb/)
      selfrefo(:,16) = (/ &
     & 3.11402e-03_jprb, 3.35870e-03_jprb, 3.62260e-03_jprb, 3.90724e-03_jprb, 4.21424e-03_jprb, &
     & 4.54536e-03_jprb, 4.90250e-03_jprb, 5.28770e-03_jprb, 5.70317e-03_jprb, 6.15128e-03_jprb/)



if (lhook) call dr_hook('rrtm_kgb13',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb13:error reading file radrrtm")

end subroutine rrtm_kgb13
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

