! # 1 "ifsrrtm/rrtm_kgb9.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_kgb9.f90"
subroutine rrtm_kgb9

!     originally by eli j. mlawer, atmospheric & environmental research.
!     band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
use yommp0    , only : nproc, myproc

use yoerrto9 , only : kao     ,kbo     ,selfrefo   ,forrefo, fracrefao  ,&
 & fracrefbo, kao_mn2o, kbo_mn2o, kao_d, kbo_d


!     ------------------------------------------------------------------

implicit none
real(kind=jphook) :: zhook_handle



! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 31 "ifsrrtm/rrtm_kgb9.f90" 2

if (lhook) call dr_hook('rrtm_kgb9',0,zhook_handle)

if( myproc==1 )then
  read(nulrad,err=1001) kao_d,kbo_d
  kao = real(kao_d,jprb)
  kbo = real(kbo_d,jprb)
endif
if( nproc>1 )then
  call mpl_broadcast (kao,mtagrad,1,cdstring='rrtm_kgb9:')
  call mpl_broadcast (kbo,mtagrad,1,cdstring='rrtm_kgb9:')
endif

! planck fractions mapping level : p=212.7250 mb, t = 223.06 k
      fracrefao(:, 1) = (/ &
     &  1.8129e-01_jprb,1.6119e-01_jprb,1.3308e-01_jprb,1.2342e-01_jprb,1.1259e-01_jprb,9.7580e-02_jprb, &
     &  7.9176e-02_jprb,5.8541e-02_jprb,3.9084e-02_jprb,4.2419e-03_jprb,3.4314e-03_jprb,2.6935e-03_jprb, &
     &  1.9404e-03_jprb,1.2218e-03_jprb,4.5263e-04_jprb,6.0909e-05_jprb/)
      fracrefao(:, 2) = (/ &
     &  1.9665e-01_jprb,1.5640e-01_jprb,1.3101e-01_jprb,1.2153e-01_jprb,1.1037e-01_jprb,9.6043e-02_jprb, &
     &  7.7856e-02_jprb,5.7547e-02_jprb,3.8670e-02_jprb,4.1955e-03_jprb,3.4104e-03_jprb,2.6781e-03_jprb, &
     &  1.9245e-03_jprb,1.2093e-03_jprb,4.4113e-04_jprb,6.0913e-05_jprb/)
      fracrefao(:, 3) = (/ &
     &  2.0273e-01_jprb,1.5506e-01_jprb,1.3044e-01_jprb,1.2043e-01_jprb,1.0952e-01_jprb,9.5384e-02_jprb, &
     &  7.7157e-02_jprb,5.7176e-02_jprb,3.8379e-02_jprb,4.1584e-03_jprb,3.3836e-03_jprb,2.6412e-03_jprb, &
     &  1.8865e-03_jprb,1.1791e-03_jprb,4.2094e-04_jprb,4.7410e-05_jprb/)
      fracrefao(:, 4) = (/ &
     &  2.0272e-01_jprb,1.5963e-01_jprb,1.2913e-01_jprb,1.2060e-01_jprb,1.0820e-01_jprb,9.4685e-02_jprb, &
     &  7.6544e-02_jprb,5.6851e-02_jprb,3.8155e-02_jprb,4.0913e-03_jprb,3.3442e-03_jprb,2.6054e-03_jprb, &
     &  1.8875e-03_jprb,1.1263e-03_jprb,3.7743e-04_jprb,4.7410e-05_jprb/)
      fracrefao(:, 5) = (/ &
     &  2.0280e-01_jprb,1.6353e-01_jprb,1.2910e-01_jprb,1.1968e-01_jprb,1.0725e-01_jprb,9.4112e-02_jprb, &
     &  7.5828e-02_jprb,5.6526e-02_jprb,3.7972e-02_jprb,4.0205e-03_jprb,3.3063e-03_jprb,2.5681e-03_jprb, &
     &  1.8386e-03_jprb,1.0757e-03_jprb,3.5301e-04_jprb,4.7410e-05_jprb/)
      fracrefao(:, 6) = (/ &
     &  2.0294e-01_jprb,1.6840e-01_jprb,1.2852e-01_jprb,1.1813e-01_jprb,1.0724e-01_jprb,9.2946e-02_jprb, &
     &  7.5029e-02_jprb,5.6158e-02_jprb,3.7744e-02_jprb,3.9632e-03_jprb,3.2434e-03_jprb,2.5275e-03_jprb, &
     &  1.7558e-03_jprb,1.0080e-03_jprb,3.5301e-04_jprb,4.7410e-05_jprb/)
      fracrefao(:, 7) = (/ &
     &  2.0313e-01_jprb,1.7390e-01_jprb,1.2864e-01_jprb,1.1689e-01_jprb,1.0601e-01_jprb,9.1791e-02_jprb, &
     &  7.4224e-02_jprb,5.5500e-02_jprb,3.7374e-02_jprb,3.9214e-03_jprb,3.1984e-03_jprb,2.4162e-03_jprb, &
     &  1.6394e-03_jprb,9.7275e-04_jprb,3.5299e-04_jprb,4.7410e-05_jprb/)
      fracrefao(:, 8) = (/ &
     &  2.0332e-01_jprb,1.7800e-01_jprb,1.3286e-01_jprb,1.1555e-01_jprb,1.0407e-01_jprb,9.0475e-02_jprb, &
     &  7.2452e-02_jprb,5.4566e-02_jprb,3.6677e-02_jprb,3.7889e-03_jprb,3.0351e-03_jprb,2.2587e-03_jprb, &
     &  1.5764e-03_jprb,9.7270e-04_jprb,3.5300e-04_jprb,4.7410e-05_jprb/)
      fracrefao(:, 9) = (/ &
     &  1.9624e-01_jprb,1.6519e-01_jprb,1.3663e-01_jprb,1.1535e-01_jprb,1.0719e-01_jprb,9.4156e-02_jprb, &
     &  7.6745e-02_jprb,5.6987e-02_jprb,3.8135e-02_jprb,4.1626e-03_jprb,3.4243e-03_jprb,2.7116e-03_jprb, &
     &  1.7095e-03_jprb,9.7271e-04_jprb,3.5299e-04_jprb,4.7410e-05_jprb/)

! planck fraction mapping level : p=3.20e-2 mb, t = 197.92 k
      fracrefbo(:) = (/ &
     &  2.0914e-01_jprb,1.5077e-01_jprb,1.2878e-01_jprb,1.1856e-01_jprb,1.0695e-01_jprb,9.3048e-02_jprb, &
     &  7.7645e-02_jprb,6.0785e-02_jprb,4.0642e-02_jprb,4.0499e-03_jprb,3.3931e-03_jprb,2.6363e-03_jprb, &
     &  1.9151e-03_jprb,1.1963e-03_jprb,4.3471e-04_jprb,5.1421e-05_jprb/)


!     ------------------------------------------------------------------

!     the array kao contains absorption coefs at the 16 chosen g-values 
!     for a range of pressure levels> ~100mb, temperatures, and binary
!     species parameters (see taumol.f for definition).  the first 
!     index in the array, js, runs from 1 to 11, and corresponds to 
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
!     level below 100~ mb.   the first index in the array, js, runs
!     from 1 to 10, and corresponds to different gas column amount ratios,
!     as expressed through the binary species parameter eta, defined as
!     eta = gas1/(gas1 + (rat) * gas2), where rat is the 
!     ratio of the reference mls column amount value of gas 1 
!     to that of gas2.  the second index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the third index 
!     runs over the g-channel (1 to 16).

      kao_mn2o( 1, :, 1) = (/ &
     & 5.41078e-02_jprb, 5.59051e-02_jprb, 5.77620e-02_jprb, 5.96805e-02_jprb, 6.16628e-02_jprb, &
     & 6.37110e-02_jprb, 6.58272e-02_jprb, 6.80137e-02_jprb, 7.02728e-02_jprb, 7.26069e-02_jprb, &
     & 7.50185e-02_jprb, 7.75103e-02_jprb, 8.00848e-02_jprb, 8.27449e-02_jprb, 8.54933e-02_jprb, &
     & 8.83330e-02_jprb, 9.12670e-02_jprb, 9.42984e-02_jprb, 9.74306e-02_jprb/)
      kao_mn2o( 2, :, 1) = (/ &
     & 1.19602e-01_jprb, 1.22963e-01_jprb, 1.26417e-01_jprb, 1.29969e-01_jprb, 1.33621e-01_jprb, &
     & 1.37375e-01_jprb, 1.41235e-01_jprb, 1.45203e-01_jprb, 1.49283e-01_jprb, 1.53477e-01_jprb, &
     & 1.57789e-01_jprb, 1.62223e-01_jprb, 1.66780e-01_jprb, 1.71466e-01_jprb, 1.76284e-01_jprb, &
     & 1.81237e-01_jprb, 1.86329e-01_jprb, 1.91564e-01_jprb, 1.96946e-01_jprb/)
      kao_mn2o( 3, :, 1) = (/ &
     & 1.49614e-01_jprb, 1.53427e-01_jprb, 1.57337e-01_jprb, 1.61346e-01_jprb, 1.65457e-01_jprb, &
     & 1.69674e-01_jprb, 1.73997e-01_jprb, 1.78431e-01_jprb, 1.82978e-01_jprb, 1.87641e-01_jprb, &
     & 1.92422e-01_jprb, 1.97325e-01_jprb, 2.02354e-01_jprb, 2.07510e-01_jprb, 2.12798e-01_jprb, &
     & 2.18221e-01_jprb, 2.23781e-01_jprb, 2.29484e-01_jprb, 2.35331e-01_jprb/)
      kao_mn2o( 4, :, 1) = (/ &
     & 1.80029e-01_jprb, 1.84202e-01_jprb, 1.88472e-01_jprb, 1.92841e-01_jprb, 1.97311e-01_jprb, &
     & 2.01884e-01_jprb, 2.06564e-01_jprb, 2.11352e-01_jprb, 2.16252e-01_jprb, 2.21264e-01_jprb, &
     & 2.26393e-01_jprb, 2.31641e-01_jprb, 2.37010e-01_jprb, 2.42504e-01_jprb, 2.48126e-01_jprb, &
     & 2.53877e-01_jprb, 2.59762e-01_jprb, 2.65783e-01_jprb, 2.71944e-01_jprb/)
      kao_mn2o( 5, :, 1) = (/ &
     & 2.08279e-01_jprb, 2.13029e-01_jprb, 2.17888e-01_jprb, 2.22858e-01_jprb, 2.27941e-01_jprb, &
     & 2.33140e-01_jprb, 2.38458e-01_jprb, 2.43897e-01_jprb, 2.49460e-01_jprb, 2.55150e-01_jprb, &
     & 2.60969e-01_jprb, 2.66922e-01_jprb, 2.73010e-01_jprb, 2.79237e-01_jprb, 2.85606e-01_jprb, &
     & 2.92120e-01_jprb, 2.98783e-01_jprb, 3.05598e-01_jprb, 3.12568e-01_jprb/)
      kao_mn2o( 6, :, 1) = (/ &
     & 2.17336e-01_jprb, 2.22571e-01_jprb, 2.27931e-01_jprb, 2.33421e-01_jprb, 2.39043e-01_jprb, &
     & 2.44801e-01_jprb, 2.50697e-01_jprb, 2.56735e-01_jprb, 2.62918e-01_jprb, 2.69251e-01_jprb, &
     & 2.75735e-01_jprb, 2.82377e-01_jprb, 2.89178e-01_jprb, 2.96142e-01_jprb, 3.03275e-01_jprb, &
     & 3.10579e-01_jprb, 3.18060e-01_jprb, 3.25720e-01_jprb, 3.33565e-01_jprb/)
      kao_mn2o( 7, :, 1) = (/ &
     & 2.23903e-01_jprb, 2.29349e-01_jprb, 2.34926e-01_jprb, 2.40640e-01_jprb, 2.46493e-01_jprb, &
     & 2.52488e-01_jprb, 2.58628e-01_jprb, 2.64918e-01_jprb, 2.71361e-01_jprb, 2.77961e-01_jprb, &
     & 2.84721e-01_jprb, 2.91646e-01_jprb, 2.98739e-01_jprb, 3.06005e-01_jprb, 3.13447e-01_jprb, &
     & 3.21070e-01_jprb, 3.28879e-01_jprb, 3.36877e-01_jprb, 3.45071e-01_jprb/)
      kao_mn2o( 8, :, 1) = (/ &
     & 2.23400e-01_jprb, 2.28843e-01_jprb, 2.34418e-01_jprb, 2.40130e-01_jprb, 2.45980e-01_jprb, &
     & 2.51973e-01_jprb, 2.58112e-01_jprb, 2.64401e-01_jprb, 2.70843e-01_jprb, 2.77442e-01_jprb, &
     & 2.84202e-01_jprb, 2.91126e-01_jprb, 2.98219e-01_jprb, 3.05485e-01_jprb, 3.12928e-01_jprb, &
     & 3.20552e-01_jprb, 3.28362e-01_jprb, 3.36362e-01_jprb, 3.44557e-01_jprb/)
      kao_mn2o( 9, :, 1) = (/ &
     & 1.89279e-01_jprb, 1.94423e-01_jprb, 1.99707e-01_jprb, 2.05135e-01_jprb, 2.10710e-01_jprb, &
     & 2.16437e-01_jprb, 2.22319e-01_jprb, 2.28361e-01_jprb, 2.34568e-01_jprb, 2.40943e-01_jprb, &
     & 2.47492e-01_jprb, 2.54218e-01_jprb, 2.61127e-01_jprb, 2.68224e-01_jprb, 2.75514e-01_jprb, &
     & 2.83002e-01_jprb, 2.90694e-01_jprb, 2.98594e-01_jprb, 3.06709e-01_jprb/)
      kao_mn2o( 1, :, 2) = (/ &
     & 9.46669e-02_jprb, 9.77137e-02_jprb, 1.00858e-01_jprb, 1.04104e-01_jprb, 1.07455e-01_jprb, &
     & 1.10913e-01_jprb, 1.14483e-01_jprb, 1.18167e-01_jprb, 1.21971e-01_jprb, 1.25896e-01_jprb, &
     & 1.29948e-01_jprb, 1.34130e-01_jprb, 1.38447e-01_jprb, 1.42903e-01_jprb, 1.47502e-01_jprb, &
     & 1.52249e-01_jprb, 1.57149e-01_jprb, 1.62207e-01_jprb, 1.67427e-01_jprb/)
      kao_mn2o( 2, :, 2) = (/ &
     & 5.11901e-01_jprb, 5.24950e-01_jprb, 5.38331e-01_jprb, 5.52053e-01_jprb, 5.66125e-01_jprb, &
     & 5.80556e-01_jprb, 5.95354e-01_jprb, 6.10530e-01_jprb, 6.26093e-01_jprb, 6.42052e-01_jprb, &
     & 6.58418e-01_jprb, 6.75202e-01_jprb, 6.92413e-01_jprb, 7.10062e-01_jprb, 7.28162e-01_jprb, &
     & 7.46723e-01_jprb, 7.65758e-01_jprb, 7.85277e-01_jprb, 8.05294e-01_jprb/)
      kao_mn2o( 3, :, 2) = (/ &
     & 8.32946e-01_jprb, 8.45780e-01_jprb, 8.58813e-01_jprb, 8.72046e-01_jprb, 8.85482e-01_jprb, &
     & 8.99126e-01_jprb, 9.12980e-01_jprb, 9.27048e-01_jprb, 9.41332e-01_jprb, 9.55837e-01_jprb, &
     & 9.70565e-01_jprb, 9.85520e-01_jprb, 1.00070e+00_jprb, 1.01612e+00_jprb, 1.03178e+00_jprb, &
     & 1.04768e+00_jprb, 1.06382e+00_jprb, 1.08021e+00_jprb, 1.09686e+00_jprb/)
      kao_mn2o( 4, :, 2) = (/ &
     & 1.04032e+00_jprb, 1.05475e+00_jprb, 1.06937e+00_jprb, 1.08419e+00_jprb, 1.09922e+00_jprb, &
     & 1.11446e+00_jprb, 1.12991e+00_jprb, 1.14557e+00_jprb, 1.16145e+00_jprb, 1.17755e+00_jprb, &
     & 1.19387e+00_jprb, 1.21042e+00_jprb, 1.22720e+00_jprb, 1.24421e+00_jprb, 1.26146e+00_jprb, &
     & 1.27895e+00_jprb, 1.29668e+00_jprb, 1.31465e+00_jprb, 1.33287e+00_jprb/)
      kao_mn2o( 5, :, 2) = (/ &
     & 1.22685e+00_jprb, 1.24267e+00_jprb, 1.25870e+00_jprb, 1.27493e+00_jprb, 1.29137e+00_jprb, &
     & 1.30803e+00_jprb, 1.32490e+00_jprb, 1.34199e+00_jprb, 1.35930e+00_jprb, 1.37683e+00_jprb, &
     & 1.39459e+00_jprb, 1.41257e+00_jprb, 1.43079e+00_jprb, 1.44925e+00_jprb, 1.46794e+00_jprb, &
     & 1.48687e+00_jprb, 1.50605e+00_jprb, 1.52547e+00_jprb, 1.54515e+00_jprb/)
      kao_mn2o( 6, :, 2) = (/ &
     & 1.53781e+00_jprb, 1.55206e+00_jprb, 1.56645e+00_jprb, 1.58097e+00_jprb, 1.59563e+00_jprb, &
     & 1.61042e+00_jprb, 1.62535e+00_jprb, 1.64042e+00_jprb, 1.65562e+00_jprb, 1.67097e+00_jprb, &
     & 1.68646e+00_jprb, 1.70210e+00_jprb, 1.71788e+00_jprb, 1.73380e+00_jprb, 1.74987e+00_jprb, &
     & 1.76610e+00_jprb, 1.78247e+00_jprb, 1.79899e+00_jprb, 1.81567e+00_jprb/)
      kao_mn2o( 7, :, 2) = (/ &
     & 1.90476e+00_jprb, 1.91858e+00_jprb, 1.93251e+00_jprb, 1.94653e+00_jprb, 1.96065e+00_jprb, &
     & 1.97488e+00_jprb, 1.98921e+00_jprb, 2.00365e+00_jprb, 2.01819e+00_jprb, 2.03283e+00_jprb, &
     & 2.04758e+00_jprb, 2.06244e+00_jprb, 2.07741e+00_jprb, 2.09248e+00_jprb, 2.10767e+00_jprb, &
     & 2.12296e+00_jprb, 2.13837e+00_jprb, 2.15388e+00_jprb, 2.16951e+00_jprb/)
      kao_mn2o( 8, :, 2) = (/ &
     & 2.38211e+00_jprb, 2.39819e+00_jprb, 2.41438e+00_jprb, 2.43068e+00_jprb, 2.44709e+00_jprb, &
     & 2.46361e+00_jprb, 2.48024e+00_jprb, 2.49699e+00_jprb, 2.51384e+00_jprb, 2.53082e+00_jprb, &
     & 2.54790e+00_jprb, 2.56510e+00_jprb, 2.58242e+00_jprb, 2.59985e+00_jprb, 2.61741e+00_jprb, &
     & 2.63508e+00_jprb, 2.65287e+00_jprb, 2.67078e+00_jprb, 2.68881e+00_jprb/)
      kao_mn2o( 9, :, 2) = (/ &
     & 1.26464e+00_jprb, 1.28107e+00_jprb, 1.29772e+00_jprb, 1.31458e+00_jprb, 1.33166e+00_jprb, &
     & 1.34896e+00_jprb, 1.36649e+00_jprb, 1.38424e+00_jprb, 1.40223e+00_jprb, 1.42044e+00_jprb, &
     & 1.43890e+00_jprb, 1.45759e+00_jprb, 1.47653e+00_jprb, 1.49572e+00_jprb, 1.51515e+00_jprb, &
     & 1.53483e+00_jprb, 1.55478e+00_jprb, 1.57498e+00_jprb, 1.59544e+00_jprb/)
      kao_mn2o( 1, :, 3) = (/ &
     & 1.50082e-01_jprb, 1.59095e-01_jprb, 1.68650e-01_jprb, 1.78779e-01_jprb, 1.89516e-01_jprb, &
     & 2.00898e-01_jprb, 2.12963e-01_jprb, 2.25753e-01_jprb, 2.39311e-01_jprb, 2.53683e-01_jprb, &
     & 2.68919e-01_jprb, 2.85069e-01_jprb, 3.02190e-01_jprb, 3.20339e-01_jprb, 3.39577e-01_jprb, &
     & 3.59971e-01_jprb, 3.81590e-01_jprb, 4.04507e-01_jprb, 4.28801e-01_jprb/)
      kao_mn2o( 2, :, 3) = (/ &
     & 3.09551e+00_jprb, 3.09780e+00_jprb, 3.10008e+00_jprb, 3.10237e+00_jprb, 3.10466e+00_jprb, &
     & 3.10695e+00_jprb, 3.10925e+00_jprb, 3.11154e+00_jprb, 3.11384e+00_jprb, 3.11614e+00_jprb, &
     & 3.11844e+00_jprb, 3.12074e+00_jprb, 3.12305e+00_jprb, 3.12535e+00_jprb, 3.12766e+00_jprb, &
     & 3.12997e+00_jprb, 3.13228e+00_jprb, 3.13459e+00_jprb, 3.13691e+00_jprb/)
      kao_mn2o( 3, :, 3) = (/ &
     & 4.42661e+00_jprb, 4.40858e+00_jprb, 4.39062e+00_jprb, 4.37274e+00_jprb, 4.35493e+00_jprb, &
     & 4.33719e+00_jprb, 4.31953e+00_jprb, 4.30193e+00_jprb, 4.28441e+00_jprb, 4.26696e+00_jprb, &
     & 4.24958e+00_jprb, 4.23227e+00_jprb, 4.21503e+00_jprb, 4.19787e+00_jprb, 4.18077e+00_jprb, &
     & 4.16374e+00_jprb, 4.14678e+00_jprb, 4.12989e+00_jprb, 4.11307e+00_jprb/)
      kao_mn2o( 4, :, 3) = (/ &
     & 5.77864e+00_jprb, 5.74085e+00_jprb, 5.70331e+00_jprb, 5.66602e+00_jprb, 5.62897e+00_jprb, &
     & 5.59216e+00_jprb, 5.55559e+00_jprb, 5.51926e+00_jprb, 5.48317e+00_jprb, 5.44731e+00_jprb, &
     & 5.41169e+00_jprb, 5.37631e+00_jprb, 5.34115e+00_jprb, 5.30622e+00_jprb, 5.27152e+00_jprb, &
     & 5.23705e+00_jprb, 5.20281e+00_jprb, 5.16879e+00_jprb, 5.13499e+00_jprb/)
      kao_mn2o( 5, :, 3) = (/ &
     & 7.17294e+00_jprb, 7.10890e+00_jprb, 7.04542e+00_jprb, 6.98251e+00_jprb, 6.92017e+00_jprb, &
     & 6.85838e+00_jprb, 6.79714e+00_jprb, 6.73645e+00_jprb, 6.67630e+00_jprb, 6.61669e+00_jprb, &
     & 6.55761e+00_jprb, 6.49905e+00_jprb, 6.44102e+00_jprb, 6.38351e+00_jprb, 6.32651e+00_jprb, &
     & 6.27003e+00_jprb, 6.21404e+00_jprb, 6.15856e+00_jprb, 6.10357e+00_jprb/)
      kao_mn2o( 6, :, 3) = (/ &
     & 9.05082e+00_jprb, 8.96132e+00_jprb, 8.87271e+00_jprb, 8.78498e+00_jprb, 8.69811e+00_jprb, &
     & 8.61210e+00_jprb, 8.52694e+00_jprb, 8.44262e+00_jprb, 8.35914e+00_jprb, 8.27648e+00_jprb, &
     & 8.19464e+00_jprb, 8.11361e+00_jprb, 8.03338e+00_jprb, 7.95394e+00_jprb, 7.87529e+00_jprb, &
     & 7.79742e+00_jprb, 7.72032e+00_jprb, 7.64398e+00_jprb, 7.56839e+00_jprb/)
      kao_mn2o( 7, :, 3) = (/ &
     & 1.18749e+01_jprb, 1.17648e+01_jprb, 1.16556e+01_jprb, 1.15475e+01_jprb, 1.14403e+01_jprb, &
     & 1.13342e+01_jprb, 1.12290e+01_jprb, 1.11248e+01_jprb, 1.10216e+01_jprb, 1.09194e+01_jprb, &
     & 1.08180e+01_jprb, 1.07177e+01_jprb, 1.06182e+01_jprb, 1.05197e+01_jprb, 1.04221e+01_jprb, &
     & 1.03254e+01_jprb, 1.02296e+01_jprb, 1.01347e+01_jprb, 1.00407e+01_jprb/)
      kao_mn2o( 8, :, 3) = (/ &
     & 1.41428e+01_jprb, 1.40323e+01_jprb, 1.39227e+01_jprb, 1.38139e+01_jprb, 1.37060e+01_jprb, &
     & 1.35989e+01_jprb, 1.34927e+01_jprb, 1.33873e+01_jprb, 1.32827e+01_jprb, 1.31790e+01_jprb, &
     & 1.30760e+01_jprb, 1.29738e+01_jprb, 1.28725e+01_jprb, 1.27719e+01_jprb, 1.26722e+01_jprb, &
     & 1.25732e+01_jprb, 1.24750e+01_jprb, 1.23775e+01_jprb, 1.22808e+01_jprb/)
      kao_mn2o( 9, :, 3) = (/ &
     & 7.34993e+00_jprb, 7.29335e+00_jprb, 7.23720e+00_jprb, 7.18149e+00_jprb, 7.12620e+00_jprb, &
     & 7.07134e+00_jprb, 7.01690e+00_jprb, 6.96288e+00_jprb, 6.90928e+00_jprb, 6.85609e+00_jprb, &
     & 6.80331e+00_jprb, 6.75094e+00_jprb, 6.69897e+00_jprb, 6.64739e+00_jprb, 6.59622e+00_jprb, &
     & 6.54544e+00_jprb, 6.49505e+00_jprb, 6.44505e+00_jprb, 6.39543e+00_jprb/)
      kao_mn2o( 1, :, 4) = (/ &
     & 6.11248e-01_jprb, 6.37225e-01_jprb, 6.64306e-01_jprb, 6.92538e-01_jprb, 7.21970e-01_jprb, &
     & 7.52653e-01_jprb, 7.84639e-01_jprb, 8.17985e-01_jprb, 8.52749e-01_jprb, 8.88989e-01_jprb, &
     & 9.26770e-01_jprb, 9.66157e-01_jprb, 1.00722e+00_jprb, 1.05002e+00_jprb, 1.09465e+00_jprb, &
     & 1.14117e+00_jprb, 1.18967e+00_jprb, 1.24022e+00_jprb, 1.29293e+00_jprb/)
      kao_mn2o( 2, :, 4) = (/ &
     & 5.07253e+00_jprb, 5.05299e+00_jprb, 5.03353e+00_jprb, 5.01414e+00_jprb, 4.99483e+00_jprb, &
     & 4.97559e+00_jprb, 4.95642e+00_jprb, 4.93733e+00_jprb, 4.91831e+00_jprb, 4.89937e+00_jprb, &
     & 4.88050e+00_jprb, 4.86170e+00_jprb, 4.84297e+00_jprb, 4.82432e+00_jprb, 4.80573e+00_jprb, &
     & 4.78722e+00_jprb, 4.76878e+00_jprb, 4.75042e+00_jprb, 4.73212e+00_jprb/)
      kao_mn2o( 3, :, 4) = (/ &
     & 7.45829e+00_jprb, 7.42266e+00_jprb, 7.38719e+00_jprb, 7.35190e+00_jprb, 7.31677e+00_jprb, &
     & 7.28181e+00_jprb, 7.24702e+00_jprb, 7.21240e+00_jprb, 7.17794e+00_jprb, 7.14364e+00_jprb, &
     & 7.10951e+00_jprb, 7.07554e+00_jprb, 7.04173e+00_jprb, 7.00809e+00_jprb, 6.97461e+00_jprb, &
     & 6.94128e+00_jprb, 6.90812e+00_jprb, 6.87511e+00_jprb, 6.84226e+00_jprb/)
      kao_mn2o( 4, :, 4) = (/ &
     & 9.58893e+00_jprb, 9.54796e+00_jprb, 9.50716e+00_jprb, 9.46654e+00_jprb, 9.42609e+00_jprb, &
     & 9.38581e+00_jprb, 9.34571e+00_jprb, 9.30578e+00_jprb, 9.26602e+00_jprb, 9.22642e+00_jprb, &
     & 9.18700e+00_jprb, 9.14775e+00_jprb, 9.10866e+00_jprb, 9.06974e+00_jprb, 9.03099e+00_jprb, &
     & 8.99240e+00_jprb, 8.95398e+00_jprb, 8.91572e+00_jprb, 8.87762e+00_jprb/)
      kao_mn2o( 5, :, 4) = (/ &
     & 1.16344e+01_jprb, 1.16012e+01_jprb, 1.15681e+01_jprb, 1.15351e+01_jprb, 1.15022e+01_jprb, &
     & 1.14694e+01_jprb, 1.14366e+01_jprb, 1.14040e+01_jprb, 1.13715e+01_jprb, 1.13390e+01_jprb, &
     & 1.13067e+01_jprb, 1.12744e+01_jprb, 1.12422e+01_jprb, 1.12102e+01_jprb, 1.11782e+01_jprb, &
     & 1.11463e+01_jprb, 1.11145e+01_jprb, 1.10828e+01_jprb, 1.10511e+01_jprb/)
      kao_mn2o( 6, :, 4) = (/ &
     & 1.12460e+01_jprb, 1.12402e+01_jprb, 1.12344e+01_jprb, 1.12286e+01_jprb, 1.12228e+01_jprb, &
     & 1.12170e+01_jprb, 1.12112e+01_jprb, 1.12055e+01_jprb, 1.11997e+01_jprb, 1.11939e+01_jprb, &
     & 1.11882e+01_jprb, 1.11824e+01_jprb, 1.11766e+01_jprb, 1.11709e+01_jprb, 1.11651e+01_jprb, &
     & 1.11594e+01_jprb, 1.11536e+01_jprb, 1.11479e+01_jprb, 1.11421e+01_jprb/)
      kao_mn2o( 7, :, 4) = (/ &
     & 8.89265e+00_jprb, 8.91419e+00_jprb, 8.93578e+00_jprb, 8.95743e+00_jprb, 8.97913e+00_jprb, &
     & 9.00088e+00_jprb, 9.02268e+00_jprb, 9.04454e+00_jprb, 9.06645e+00_jprb, 9.08841e+00_jprb, &
     & 9.11043e+00_jprb, 9.13250e+00_jprb, 9.15462e+00_jprb, 9.17680e+00_jprb, 9.19903e+00_jprb, &
     & 9.22131e+00_jprb, 9.24365e+00_jprb, 9.26604e+00_jprb, 9.28849e+00_jprb/)
      kao_mn2o( 8, :, 4) = (/ &
     & 6.83933e+00_jprb, 6.86688e+00_jprb, 6.89453e+00_jprb, 6.92230e+00_jprb, 6.95018e+00_jprb, &
     & 6.97817e+00_jprb, 7.00627e+00_jprb, 7.03449e+00_jprb, 7.06282e+00_jprb, 7.09126e+00_jprb, &
     & 7.11982e+00_jprb, 7.14850e+00_jprb, 7.17729e+00_jprb, 7.20619e+00_jprb, 7.23521e+00_jprb, &
     & 7.26435e+00_jprb, 7.29361e+00_jprb, 7.32298e+00_jprb, 7.35248e+00_jprb/)
      kao_mn2o( 9, :, 4) = (/ &
     & 1.10637e+01_jprb, 1.10232e+01_jprb, 1.09829e+01_jprb, 1.09427e+01_jprb, 1.09026e+01_jprb, &
     & 1.08627e+01_jprb, 1.08230e+01_jprb, 1.07833e+01_jprb, 1.07439e+01_jprb, 1.07046e+01_jprb, &
     & 1.06654e+01_jprb, 1.06263e+01_jprb, 1.05875e+01_jprb, 1.05487e+01_jprb, 1.05101e+01_jprb, &
     & 1.04716e+01_jprb, 1.04333e+01_jprb, 1.03951e+01_jprb, 1.03571e+01_jprb/)
      kao_mn2o( 1, :, 5) = (/ &
     & 2.53460e+00_jprb, 2.56050e+00_jprb, 2.58667e+00_jprb, 2.61310e+00_jprb, 2.63980e+00_jprb, &
     & 2.66678e+00_jprb, 2.69403e+00_jprb, 2.72156e+00_jprb, 2.74937e+00_jprb, 2.77746e+00_jprb, &
     & 2.80585e+00_jprb, 2.83452e+00_jprb, 2.86348e+00_jprb, 2.89275e+00_jprb, 2.92231e+00_jprb, &
     & 2.95217e+00_jprb, 2.98234e+00_jprb, 3.01281e+00_jprb, 3.04360e+00_jprb/)
      kao_mn2o( 2, :, 5) = (/ &
     & 7.45650e+00_jprb, 7.44283e+00_jprb, 7.42919e+00_jprb, 7.41557e+00_jprb, 7.40198e+00_jprb, &
     & 7.38841e+00_jprb, 7.37487e+00_jprb, 7.36135e+00_jprb, 7.34786e+00_jprb, 7.33439e+00_jprb, &
     & 7.32095e+00_jprb, 7.30753e+00_jprb, 7.29413e+00_jprb, 7.28076e+00_jprb, 7.26742e+00_jprb, &
     & 7.25410e+00_jprb, 7.24080e+00_jprb, 7.22753e+00_jprb, 7.21428e+00_jprb/)
      kao_mn2o( 3, :, 5) = (/ &
     & 1.06311e+01_jprb, 1.06110e+01_jprb, 1.05909e+01_jprb, 1.05709e+01_jprb, 1.05509e+01_jprb, &
     & 1.05310e+01_jprb, 1.05111e+01_jprb, 1.04912e+01_jprb, 1.04713e+01_jprb, 1.04516e+01_jprb, &
     & 1.04318e+01_jprb, 1.04121e+01_jprb, 1.03924e+01_jprb, 1.03727e+01_jprb, 1.03531e+01_jprb, &
     & 1.03336e+01_jprb, 1.03140e+01_jprb, 1.02945e+01_jprb, 1.02751e+01_jprb/)
      kao_mn2o( 4, :, 5) = (/ &
     & 1.03924e+01_jprb, 1.03895e+01_jprb, 1.03867e+01_jprb, 1.03838e+01_jprb, 1.03809e+01_jprb, &
     & 1.03780e+01_jprb, 1.03751e+01_jprb, 1.03722e+01_jprb, 1.03693e+01_jprb, 1.03665e+01_jprb, &
     & 1.03636e+01_jprb, 1.03607e+01_jprb, 1.03578e+01_jprb, 1.03549e+01_jprb, 1.03521e+01_jprb, &
     & 1.03492e+01_jprb, 1.03463e+01_jprb, 1.03434e+01_jprb, 1.03406e+01_jprb/)
      kao_mn2o( 5, :, 5) = (/ &
     & 7.82277e+00_jprb, 7.83872e+00_jprb, 7.85471e+00_jprb, 7.87073e+00_jprb, 7.88678e+00_jprb, &
     & 7.90287e+00_jprb, 7.91899e+00_jprb, 7.93514e+00_jprb, 7.95132e+00_jprb, 7.96754e+00_jprb, &
     & 7.98379e+00_jprb, 8.00008e+00_jprb, 8.01639e+00_jprb, 8.03274e+00_jprb, 8.04913e+00_jprb, &
     & 8.06555e+00_jprb, 8.08200e+00_jprb, 8.09848e+00_jprb, 8.11500e+00_jprb/)
      kao_mn2o( 6, :, 5) = (/ &
     & 6.05225e+00_jprb, 6.06883e+00_jprb, 6.08545e+00_jprb, 6.10212e+00_jprb, 6.11883e+00_jprb, &
     & 6.13559e+00_jprb, 6.15240e+00_jprb, 6.16925e+00_jprb, 6.18615e+00_jprb, 6.20309e+00_jprb, &
     & 6.22008e+00_jprb, 6.23712e+00_jprb, 6.25420e+00_jprb, 6.27133e+00_jprb, 6.28851e+00_jprb, &
     & 6.30574e+00_jprb, 6.32301e+00_jprb, 6.34033e+00_jprb, 6.35769e+00_jprb/)
      kao_mn2o( 7, :, 5) = (/ &
     & 5.24135e+00_jprb, 5.25696e+00_jprb, 5.27261e+00_jprb, 5.28831e+00_jprb, 5.30405e+00_jprb, &
     & 5.31984e+00_jprb, 5.33568e+00_jprb, 5.35157e+00_jprb, 5.36750e+00_jprb, 5.38348e+00_jprb, &
     & 5.39951e+00_jprb, 5.41558e+00_jprb, 5.43171e+00_jprb, 5.44788e+00_jprb, 5.46410e+00_jprb, &
     & 5.48037e+00_jprb, 5.49668e+00_jprb, 5.51305e+00_jprb, 5.52946e+00_jprb/)
      kao_mn2o( 8, :, 5) = (/ &
     & 4.40240e+00_jprb, 4.40915e+00_jprb, 4.41591e+00_jprb, 4.42268e+00_jprb, 4.42946e+00_jprb, &
     & 4.43625e+00_jprb, 4.44305e+00_jprb, 4.44986e+00_jprb, 4.45668e+00_jprb, 4.46351e+00_jprb, &
     & 4.47035e+00_jprb, 4.47720e+00_jprb, 4.48407e+00_jprb, 4.49094e+00_jprb, 4.49782e+00_jprb, &
     & 4.50472e+00_jprb, 4.51162e+00_jprb, 4.51854e+00_jprb, 4.52547e+00_jprb/)
      kao_mn2o( 9, :, 5) = (/ &
     & 8.56554e+00_jprb, 8.59185e+00_jprb, 8.61824e+00_jprb, 8.64470e+00_jprb, 8.67125e+00_jprb, &
     & 8.69788e+00_jprb, 8.72460e+00_jprb, 8.75139e+00_jprb, 8.77827e+00_jprb, 8.80523e+00_jprb, &
     & 8.83227e+00_jprb, 8.85939e+00_jprb, 8.88660e+00_jprb, 8.91389e+00_jprb, 8.94127e+00_jprb, &
     & 8.96873e+00_jprb, 8.99627e+00_jprb, 9.02390e+00_jprb, 9.05161e+00_jprb/)
      kao_mn2o( 1, :, 6) = (/ &
     & 5.78695e+00_jprb, 5.78939e+00_jprb, 5.79182e+00_jprb, 5.79426e+00_jprb, 5.79670e+00_jprb, &
     & 5.79914e+00_jprb, 5.80158e+00_jprb, 5.80403e+00_jprb, 5.80647e+00_jprb, 5.80892e+00_jprb, &
     & 5.81136e+00_jprb, 5.81381e+00_jprb, 5.81626e+00_jprb, 5.81870e+00_jprb, 5.82115e+00_jprb, &
     & 5.82361e+00_jprb, 5.82606e+00_jprb, 5.82851e+00_jprb, 5.83096e+00_jprb/)
      kao_mn2o( 2, :, 6) = (/ &
     & 1.22893e+01_jprb, 1.22556e+01_jprb, 1.22221e+01_jprb, 1.21886e+01_jprb, 1.21552e+01_jprb, &
     & 1.21220e+01_jprb, 1.20888e+01_jprb, 1.20557e+01_jprb, 1.20227e+01_jprb, 1.19898e+01_jprb, &
     & 1.19569e+01_jprb, 1.19242e+01_jprb, 1.18915e+01_jprb, 1.18590e+01_jprb, 1.18265e+01_jprb, &
     & 1.17941e+01_jprb, 1.17618e+01_jprb, 1.17296e+01_jprb, 1.16975e+01_jprb/)
      kao_mn2o( 3, :, 6) = (/ &
     & 7.93118e+00_jprb, 7.94590e+00_jprb, 7.96065e+00_jprb, 7.97542e+00_jprb, 7.99022e+00_jprb, &
     & 8.00505e+00_jprb, 8.01990e+00_jprb, 8.03478e+00_jprb, 8.04970e+00_jprb, 8.06463e+00_jprb, &
     & 8.07960e+00_jprb, 8.09459e+00_jprb, 8.10961e+00_jprb, 8.12466e+00_jprb, 8.13974e+00_jprb, &
     & 8.15485e+00_jprb, 8.16998e+00_jprb, 8.18514e+00_jprb, 8.20033e+00_jprb/)
      kao_mn2o( 4, :, 6) = (/ &
     & 4.08899e+00_jprb, 4.11435e+00_jprb, 4.13988e+00_jprb, 4.16556e+00_jprb, 4.19140e+00_jprb, &
     & 4.21740e+00_jprb, 4.24357e+00_jprb, 4.26989e+00_jprb, 4.29638e+00_jprb, 4.32304e+00_jprb, &
     & 4.34985e+00_jprb, 4.37684e+00_jprb, 4.40399e+00_jprb, 4.43131e+00_jprb, 4.45880e+00_jprb, &
     & 4.48646e+00_jprb, 4.51430e+00_jprb, 4.54230e+00_jprb, 4.57048e+00_jprb/)
      kao_mn2o( 5, :, 6) = (/ &
     & 2.61358e+00_jprb, 2.64029e+00_jprb, 2.66728e+00_jprb, 2.69454e+00_jprb, 2.72209e+00_jprb, &
     & 2.74991e+00_jprb, 2.77802e+00_jprb, 2.80641e+00_jprb, 2.83510e+00_jprb, 2.86408e+00_jprb, &
     & 2.89335e+00_jprb, 2.92293e+00_jprb, 2.95280e+00_jprb, 2.98299e+00_jprb, 3.01348e+00_jprb, &
     & 3.04428e+00_jprb, 3.07540e+00_jprb, 3.10683e+00_jprb, 3.13859e+00_jprb/)
      kao_mn2o( 6, :, 6) = (/ &
     & 2.40720e+00_jprb, 2.43430e+00_jprb, 2.46169e+00_jprb, 2.48940e+00_jprb, 2.51741e+00_jprb, &
     & 2.54575e+00_jprb, 2.57440e+00_jprb, 2.60337e+00_jprb, 2.63267e+00_jprb, 2.66230e+00_jprb, &
     & 2.69226e+00_jprb, 2.72256e+00_jprb, 2.75320e+00_jprb, 2.78419e+00_jprb, 2.81552e+00_jprb, &
     & 2.84721e+00_jprb, 2.87925e+00_jprb, 2.91166e+00_jprb, 2.94443e+00_jprb/)
      kao_mn2o( 7, :, 6) = (/ &
     & 1.99607e+00_jprb, 2.01725e+00_jprb, 2.03865e+00_jprb, 2.06028e+00_jprb, 2.08214e+00_jprb, &
     & 2.10423e+00_jprb, 2.12655e+00_jprb, 2.14912e+00_jprb, 2.17192e+00_jprb, 2.19496e+00_jprb, &
     & 2.21825e+00_jprb, 2.24179e+00_jprb, 2.26557e+00_jprb, 2.28961e+00_jprb, 2.31390e+00_jprb, &
     & 2.33845e+00_jprb, 2.36326e+00_jprb, 2.38834e+00_jprb, 2.41368e+00_jprb/)
      kao_mn2o( 8, :, 6) = (/ &
     & 1.94150e+00_jprb, 1.96398e+00_jprb, 1.98671e+00_jprb, 2.00971e+00_jprb, 2.03298e+00_jprb, &
     & 2.05651e+00_jprb, 2.08032e+00_jprb, 2.10440e+00_jprb, 2.12876e+00_jprb, 2.15341e+00_jprb, &
     & 2.17834e+00_jprb, 2.20355e+00_jprb, 2.22906e+00_jprb, 2.25487e+00_jprb, 2.28097e+00_jprb, &
     & 2.30738e+00_jprb, 2.33409e+00_jprb, 2.36111e+00_jprb, 2.38844e+00_jprb/)
      kao_mn2o( 9, :, 6) = (/ &
     & 2.47259e+00_jprb, 2.48950e+00_jprb, 2.50653e+00_jprb, 2.52367e+00_jprb, 2.54093e+00_jprb, &
     & 2.55831e+00_jprb, 2.57581e+00_jprb, 2.59343e+00_jprb, 2.61117e+00_jprb, 2.62903e+00_jprb, &
     & 2.64701e+00_jprb, 2.66511e+00_jprb, 2.68334e+00_jprb, 2.70169e+00_jprb, 2.72017e+00_jprb, &
     & 2.73878e+00_jprb, 2.75751e+00_jprb, 2.77637e+00_jprb, 2.79536e+00_jprb/)
      kao_mn2o( 1, :, 7) = (/ &
     & 1.23417e+01_jprb, 1.22618e+01_jprb, 1.21823e+01_jprb, 1.21034e+01_jprb, 1.20250e+01_jprb, &
     & 1.19471e+01_jprb, 1.18697e+01_jprb, 1.17928e+01_jprb, 1.17164e+01_jprb, 1.16405e+01_jprb, &
     & 1.15651e+01_jprb, 1.14901e+01_jprb, 1.14157e+01_jprb, 1.13417e+01_jprb, 1.12683e+01_jprb, &
     & 1.11952e+01_jprb, 1.11227e+01_jprb, 1.10507e+01_jprb, 1.09791e+01_jprb/)
      kao_mn2o( 2, :, 7) = (/ &
     & 9.30957e+00_jprb, 9.32775e+00_jprb, 9.34597e+00_jprb, 9.36421e+00_jprb, 9.38250e+00_jprb, &
     & 9.40082e+00_jprb, 9.41918e+00_jprb, 9.43757e+00_jprb, 9.45600e+00_jprb, 9.47446e+00_jprb, &
     & 9.49296e+00_jprb, 9.51150e+00_jprb, 9.53007e+00_jprb, 9.54868e+00_jprb, 9.56732e+00_jprb, &
     & 9.58601e+00_jprb, 9.60472e+00_jprb, 9.62348e+00_jprb, 9.64227e+00_jprb/)
      kao_mn2o( 3, :, 7) = (/ &
     & 4.15867e+00_jprb, 4.19254e+00_jprb, 4.22668e+00_jprb, 4.26110e+00_jprb, 4.29581e+00_jprb, &
     & 4.33079e+00_jprb, 4.36606e+00_jprb, 4.40162e+00_jprb, 4.43747e+00_jprb, 4.47360e+00_jprb, &
     & 4.51004e+00_jprb, 4.54677e+00_jprb, 4.58380e+00_jprb, 4.62113e+00_jprb, 4.65876e+00_jprb, &
     & 4.69670e+00_jprb, 4.73495e+00_jprb, 4.77351e+00_jprb, 4.81239e+00_jprb/)
      kao_mn2o( 4, :, 7) = (/ &
     & 3.55634e+00_jprb, 3.59382e+00_jprb, 3.63169e+00_jprb, 3.66996e+00_jprb, 3.70863e+00_jprb, &
     & 3.74771e+00_jprb, 3.78720e+00_jprb, 3.82711e+00_jprb, 3.86744e+00_jprb, 3.90820e+00_jprb, &
     & 3.94938e+00_jprb, 3.99100e+00_jprb, 4.03305e+00_jprb, 4.07555e+00_jprb, 4.11850e+00_jprb, &
     & 4.16190e+00_jprb, 4.20576e+00_jprb, 4.25008e+00_jprb, 4.29486e+00_jprb/)
      kao_mn2o( 5, :, 7) = (/ &
     & 3.09468e+00_jprb, 3.12655e+00_jprb, 3.15876e+00_jprb, 3.19129e+00_jprb, 3.22416e+00_jprb, &
     & 3.25737e+00_jprb, 3.29092e+00_jprb, 3.32482e+00_jprb, 3.35907e+00_jprb, 3.39366e+00_jprb, &
     & 3.42862e+00_jprb, 3.46393e+00_jprb, 3.49961e+00_jprb, 3.53566e+00_jprb, 3.57208e+00_jprb, &
     & 3.60887e+00_jprb, 3.64604e+00_jprb, 3.68360e+00_jprb, 3.72154e+00_jprb/)
      kao_mn2o( 6, :, 7) = (/ &
     & 2.75473e+00_jprb, 2.78356e+00_jprb, 2.81268e+00_jprb, 2.84211e+00_jprb, 2.87185e+00_jprb, &
     & 2.90190e+00_jprb, 2.93227e+00_jprb, 2.96295e+00_jprb, 2.99395e+00_jprb, 3.02528e+00_jprb, &
     & 3.05694e+00_jprb, 3.08892e+00_jprb, 3.12125e+00_jprb, 3.15391e+00_jprb, 3.18691e+00_jprb, &
     & 3.22025e+00_jprb, 3.25395e+00_jprb, 3.28800e+00_jprb, 3.32240e+00_jprb/)
      kao_mn2o( 7, :, 7) = (/ &
     & 2.68587e+00_jprb, 2.71431e+00_jprb, 2.74306e+00_jprb, 2.77211e+00_jprb, 2.80146e+00_jprb, &
     & 2.83113e+00_jprb, 2.86111e+00_jprb, 2.89141e+00_jprb, 2.92203e+00_jprb, 2.95298e+00_jprb, &
     & 2.98425e+00_jprb, 3.01585e+00_jprb, 3.04779e+00_jprb, 3.08006e+00_jprb, 3.11268e+00_jprb, &
     & 3.14564e+00_jprb, 3.17896e+00_jprb, 3.21262e+00_jprb, 3.24664e+00_jprb/)
      kao_mn2o( 8, :, 7) = (/ &
     & 2.54778e+00_jprb, 2.57461e+00_jprb, 2.60173e+00_jprb, 2.62914e+00_jprb, 2.65683e+00_jprb, &
     & 2.68482e+00_jprb, 2.71310e+00_jprb, 2.74168e+00_jprb, 2.77056e+00_jprb, 2.79974e+00_jprb, &
     & 2.82923e+00_jprb, 2.85903e+00_jprb, 2.88915e+00_jprb, 2.91958e+00_jprb, 2.95033e+00_jprb, &
     & 2.98141e+00_jprb, 3.01282e+00_jprb, 3.04455e+00_jprb, 3.07662e+00_jprb/)
      kao_mn2o( 9, :, 7) = (/ &
     & 2.78137e+00_jprb, 2.80957e+00_jprb, 2.83805e+00_jprb, 2.86682e+00_jprb, 2.89589e+00_jprb, &
     & 2.92525e+00_jprb, 2.95491e+00_jprb, 2.98486e+00_jprb, 3.01512e+00_jprb, 3.04569e+00_jprb, &
     & 3.07657e+00_jprb, 3.10776e+00_jprb, 3.13927e+00_jprb, 3.17110e+00_jprb, 3.20324e+00_jprb, &
     & 3.23572e+00_jprb, 3.26852e+00_jprb, 3.30166e+00_jprb, 3.33513e+00_jprb/)
      kao_mn2o( 1, :, 8) = (/ &
     & 2.28384e+01_jprb, 2.27450e+01_jprb, 2.26519e+01_jprb, 2.25593e+01_jprb, 2.24670e+01_jprb, &
     & 2.23751e+01_jprb, 2.22835e+01_jprb, 2.21924e+01_jprb, 2.21016e+01_jprb, 2.20112e+01_jprb, &
     & 2.19211e+01_jprb, 2.18314e+01_jprb, 2.17421e+01_jprb, 2.16532e+01_jprb, 2.15646e+01_jprb, &
     & 2.14764e+01_jprb, 2.13885e+01_jprb, 2.13010e+01_jprb, 2.12139e+01_jprb/)
      kao_mn2o( 2, :, 8) = (/ &
     & 4.48608e+00_jprb, 4.52259e+00_jprb, 4.55939e+00_jprb, 4.59649e+00_jprb, 4.63389e+00_jprb, &
     & 4.67160e+00_jprb, 4.70961e+00_jprb, 4.74793e+00_jprb, 4.78656e+00_jprb, 4.82551e+00_jprb, &
     & 4.86478e+00_jprb, 4.90436e+00_jprb, 4.94427e+00_jprb, 4.98450e+00_jprb, 5.02506e+00_jprb, &
     & 5.06595e+00_jprb, 5.10717e+00_jprb, 5.14873e+00_jprb, 5.19062e+00_jprb/)
      kao_mn2o( 3, :, 8) = (/ &
     & 3.69928e+00_jprb, 3.72584e+00_jprb, 3.75259e+00_jprb, 3.77953e+00_jprb, 3.80666e+00_jprb, &
     & 3.83399e+00_jprb, 3.86152e+00_jprb, 3.88924e+00_jprb, 3.91717e+00_jprb, 3.94529e+00_jprb, &
     & 3.97361e+00_jprb, 4.00214e+00_jprb, 4.03088e+00_jprb, 4.05982e+00_jprb, 4.08896e+00_jprb, &
     & 4.11832e+00_jprb, 4.14789e+00_jprb, 4.17767e+00_jprb, 4.20766e+00_jprb/)
      kao_mn2o( 4, :, 8) = (/ &
     & 3.17856e+00_jprb, 3.19596e+00_jprb, 3.21345e+00_jprb, 3.23104e+00_jprb, 3.24872e+00_jprb, &
     & 3.26650e+00_jprb, 3.28437e+00_jprb, 3.30235e+00_jprb, 3.32042e+00_jprb, 3.33859e+00_jprb, &
     & 3.35686e+00_jprb, 3.37523e+00_jprb, 3.39370e+00_jprb, 3.41227e+00_jprb, 3.43095e+00_jprb, &
     & 3.44972e+00_jprb, 3.46860e+00_jprb, 3.48758e+00_jprb, 3.50667e+00_jprb/)
      kao_mn2o( 5, :, 8) = (/ &
     & 3.16549e+00_jprb, 3.18288e+00_jprb, 3.20037e+00_jprb, 3.21795e+00_jprb, 3.23563e+00_jprb, &
     & 3.25340e+00_jprb, 3.27128e+00_jprb, 3.28925e+00_jprb, 3.30732e+00_jprb, 3.32549e+00_jprb, &
     & 3.34376e+00_jprb, 3.36213e+00_jprb, 3.38060e+00_jprb, 3.39917e+00_jprb, 3.41785e+00_jprb, &
     & 3.43662e+00_jprb, 3.45551e+00_jprb, 3.47449e+00_jprb, 3.49358e+00_jprb/)
      kao_mn2o( 6, :, 8) = (/ &
     & 3.16612e+00_jprb, 3.18355e+00_jprb, 3.20108e+00_jprb, 3.21870e+00_jprb, 3.23643e+00_jprb, &
     & 3.25425e+00_jprb, 3.27217e+00_jprb, 3.29018e+00_jprb, 3.30830e+00_jprb, 3.32652e+00_jprb, &
     & 3.34483e+00_jprb, 3.36325e+00_jprb, 3.38177e+00_jprb, 3.40039e+00_jprb, 3.41911e+00_jprb, &
     & 3.43794e+00_jprb, 3.45687e+00_jprb, 3.47591e+00_jprb, 3.49505e+00_jprb/)
      kao_mn2o( 7, :, 8) = (/ &
     & 3.19644e+00_jprb, 3.21419e+00_jprb, 3.23203e+00_jprb, 3.24996e+00_jprb, 3.26800e+00_jprb, &
     & 3.28614e+00_jprb, 3.30438e+00_jprb, 3.32272e+00_jprb, 3.34116e+00_jprb, 3.35970e+00_jprb, &
     & 3.37835e+00_jprb, 3.39710e+00_jprb, 3.41596e+00_jprb, 3.43492e+00_jprb, 3.45398e+00_jprb, &
     & 3.47315e+00_jprb, 3.49243e+00_jprb, 3.51181e+00_jprb, 3.53130e+00_jprb/)
      kao_mn2o( 8, :, 8) = (/ &
     & 3.35759e+00_jprb, 3.37775e+00_jprb, 3.39804e+00_jprb, 3.41845e+00_jprb, 3.43899e+00_jprb, &
     & 3.45964e+00_jprb, 3.48042e+00_jprb, 3.50133e+00_jprb, 3.52236e+00_jprb, 3.54351e+00_jprb, &
     & 3.56480e+00_jprb, 3.58621e+00_jprb, 3.60775e+00_jprb, 3.62942e+00_jprb, 3.65122e+00_jprb, &
     & 3.67315e+00_jprb, 3.69521e+00_jprb, 3.71741e+00_jprb, 3.73974e+00_jprb/)
      kao_mn2o( 9, :, 8) = (/ &
     & 3.17378e+00_jprb, 3.19135e+00_jprb, 3.20901e+00_jprb, 3.22677e+00_jprb, 3.24462e+00_jprb, &
     & 3.26258e+00_jprb, 3.28063e+00_jprb, 3.29879e+00_jprb, 3.31704e+00_jprb, 3.33540e+00_jprb, &
     & 3.35386e+00_jprb, 3.37242e+00_jprb, 3.39108e+00_jprb, 3.40984e+00_jprb, 3.42871e+00_jprb, &
     & 3.44769e+00_jprb, 3.46677e+00_jprb, 3.48595e+00_jprb, 3.50524e+00_jprb/)
      kao_mn2o( 1, :, 9) = (/ &
     & 2.09106e+01_jprb, 2.08779e+01_jprb, 2.08452e+01_jprb, 2.08126e+01_jprb, 2.07800e+01_jprb, &
     & 2.07475e+01_jprb, 2.07150e+01_jprb, 2.06826e+01_jprb, 2.06502e+01_jprb, 2.06179e+01_jprb, &
     & 2.05856e+01_jprb, 2.05534e+01_jprb, 2.05213e+01_jprb, 2.04891e+01_jprb, 2.04571e+01_jprb, &
     & 2.04251e+01_jprb, 2.03931e+01_jprb, 2.03612e+01_jprb, 2.03293e+01_jprb/)
      kao_mn2o( 2, :, 9) = (/ &
     & 2.60494e+00_jprb, 2.62757e+00_jprb, 2.65040e+00_jprb, 2.67343e+00_jprb, 2.69665e+00_jprb, &
     & 2.72008e+00_jprb, 2.74372e+00_jprb, 2.76756e+00_jprb, 2.79160e+00_jprb, 2.81586e+00_jprb, &
     & 2.84032e+00_jprb, 2.86500e+00_jprb, 2.88989e+00_jprb, 2.91500e+00_jprb, 2.94033e+00_jprb, &
     & 2.96588e+00_jprb, 2.99164e+00_jprb, 3.01764e+00_jprb, 3.04386e+00_jprb/)
      kao_mn2o( 3, :, 9) = (/ &
     & 2.42238e+00_jprb, 2.44514e+00_jprb, 2.46811e+00_jprb, 2.49130e+00_jprb, 2.51471e+00_jprb, &
     & 2.53834e+00_jprb, 2.56219e+00_jprb, 2.58626e+00_jprb, 2.61056e+00_jprb, 2.63509e+00_jprb, &
     & 2.65985e+00_jprb, 2.68484e+00_jprb, 2.71007e+00_jprb, 2.73554e+00_jprb, 2.76124e+00_jprb, &
     & 2.78718e+00_jprb, 2.81337e+00_jprb, 2.83981e+00_jprb, 2.86649e+00_jprb/)
      kao_mn2o( 4, :, 9) = (/ &
     & 2.33681e+00_jprb, 2.35961e+00_jprb, 2.38263e+00_jprb, 2.40588e+00_jprb, 2.42935e+00_jprb, &
     & 2.45305e+00_jprb, 2.47699e+00_jprb, 2.50115e+00_jprb, 2.52556e+00_jprb, 2.55020e+00_jprb, &
     & 2.57508e+00_jprb, 2.60021e+00_jprb, 2.62558e+00_jprb, 2.65119e+00_jprb, 2.67706e+00_jprb, &
     & 2.70318e+00_jprb, 2.72955e+00_jprb, 2.75619e+00_jprb, 2.78308e+00_jprb/)
      kao_mn2o( 5, :, 9) = (/ &
     & 2.26420e+00_jprb, 2.28696e+00_jprb, 2.30996e+00_jprb, 2.33319e+00_jprb, 2.35665e+00_jprb, &
     & 2.38035e+00_jprb, 2.40429e+00_jprb, 2.42846e+00_jprb, 2.45288e+00_jprb, 2.47755e+00_jprb, &
     & 2.50246e+00_jprb, 2.52763e+00_jprb, 2.55304e+00_jprb, 2.57871e+00_jprb, 2.60465e+00_jprb, &
     & 2.63084e+00_jprb, 2.65729e+00_jprb, 2.68401e+00_jprb, 2.71100e+00_jprb/)
      kao_mn2o( 6, :, 9) = (/ &
     & 2.19628e+00_jprb, 2.21902e+00_jprb, 2.24199e+00_jprb, 2.26520e+00_jprb, 2.28865e+00_jprb, &
     & 2.31234e+00_jprb, 2.33628e+00_jprb, 2.36047e+00_jprb, 2.38491e+00_jprb, 2.40959e+00_jprb, &
     & 2.43454e+00_jprb, 2.45974e+00_jprb, 2.48521e+00_jprb, 2.51094e+00_jprb, 2.53693e+00_jprb, &
     & 2.56319e+00_jprb, 2.58973e+00_jprb, 2.61654e+00_jprb, 2.64363e+00_jprb/)
      kao_mn2o( 7, :, 9) = (/ &
     & 2.07829e+00_jprb, 2.10090e+00_jprb, 2.12375e+00_jprb, 2.14685e+00_jprb, 2.17021e+00_jprb, &
     & 2.19381e+00_jprb, 2.21767e+00_jprb, 2.24180e+00_jprb, 2.26618e+00_jprb, 2.29083e+00_jprb, &
     & 2.31575e+00_jprb, 2.34094e+00_jprb, 2.36640e+00_jprb, 2.39214e+00_jprb, 2.41816e+00_jprb, &
     & 2.44446e+00_jprb, 2.47105e+00_jprb, 2.49793e+00_jprb, 2.52510e+00_jprb/)
      kao_mn2o( 8, :, 9) = (/ &
     & 1.68230e+00_jprb, 1.70305e+00_jprb, 1.72404e+00_jprb, 1.74530e+00_jprb, 1.76681e+00_jprb, &
     & 1.78860e+00_jprb, 1.81065e+00_jprb, 1.83297e+00_jprb, 1.85557e+00_jprb, 1.87845e+00_jprb, &
     & 1.90161e+00_jprb, 1.92505e+00_jprb, 1.94878e+00_jprb, 1.97281e+00_jprb, 1.99713e+00_jprb, &
     & 2.02176e+00_jprb, 2.04668e+00_jprb, 2.07191e+00_jprb, 2.09746e+00_jprb/)
      kao_mn2o( 9, :, 9) = (/ &
     & 2.23224e+00_jprb, 2.25486e+00_jprb, 2.27771e+00_jprb, 2.30079e+00_jprb, 2.32411e+00_jprb, &
     & 2.34766e+00_jprb, 2.37145e+00_jprb, 2.39548e+00_jprb, 2.41975e+00_jprb, 2.44427e+00_jprb, &
     & 2.46904e+00_jprb, 2.49406e+00_jprb, 2.51933e+00_jprb, 2.54486e+00_jprb, 2.57065e+00_jprb, &
     & 2.59670e+00_jprb, 2.62301e+00_jprb, 2.64959e+00_jprb, 2.67644e+00_jprb/)
      kao_mn2o( 1, :,10) = (/ &
     & 1.30711e+01_jprb, 1.31853e+01_jprb, 1.33004e+01_jprb, 1.34166e+01_jprb, 1.35339e+01_jprb, &
     & 1.36521e+01_jprb, 1.37714e+01_jprb, 1.38917e+01_jprb, 1.40130e+01_jprb, 1.41355e+01_jprb, &
     & 1.42590e+01_jprb, 1.43835e+01_jprb, 1.45092e+01_jprb, 1.46360e+01_jprb, 1.47638e+01_jprb, &
     & 1.48928e+01_jprb, 1.50229e+01_jprb, 1.51542e+01_jprb, 1.52866e+01_jprb/)
      kao_mn2o( 2, :,10) = (/ &
     & 2.71206e-01_jprb, 2.90148e-01_jprb, 3.10413e-01_jprb, 3.32093e-01_jprb, 3.55287e-01_jprb, &
     & 3.80102e-01_jprb, 4.06649e-01_jprb, 4.35051e-01_jprb, 4.65436e-01_jprb, 4.97943e-01_jprb, &
     & 5.32721e-01_jprb, 5.69928e-01_jprb, 6.09733e-01_jprb, 6.52319e-01_jprb, 6.97879e-01_jprb, &
     & 7.46621e-01_jprb, 7.98767e-01_jprb, 8.54555e-01_jprb, 9.14239e-01_jprb/)
      kao_mn2o( 3, :,10) = (/ &
     & 2.65609e-01_jprb, 2.84236e-01_jprb, 3.04170e-01_jprb, 3.25501e-01_jprb, 3.48329e-01_jprb, &
     & 3.72758e-01_jprb, 3.98900e-01_jprb, 4.26875e-01_jprb, 4.56812e-01_jprb, 4.88849e-01_jprb, &
     & 5.23132e-01_jprb, 5.59820e-01_jprb, 5.99080e-01_jprb, 6.41095e-01_jprb, 6.86055e-01_jprb, &
     & 7.34169e-01_jprb, 7.85657e-01_jprb, 8.40756e-01_jprb, 8.99718e-01_jprb/)
      kao_mn2o( 4, :,10) = (/ &
     & 2.55277e-01_jprb, 2.73370e-01_jprb, 2.92745e-01_jprb, 3.13494e-01_jprb, 3.35714e-01_jprb, &
     & 3.59508e-01_jprb, 3.84989e-01_jprb, 4.12275e-01_jprb, 4.41496e-01_jprb, 4.72788e-01_jprb, &
     & 5.06298e-01_jprb, 5.42183e-01_jprb, 5.80611e-01_jprb, 6.21763e-01_jprb, 6.65831e-01_jprb, &
     & 7.13023e-01_jprb, 7.63560e-01_jprb, 8.17678e-01_jprb, 8.75633e-01_jprb/)
      kao_mn2o( 5, :,10) = (/ &
     & 2.41481e-01_jprb, 2.58840e-01_jprb, 2.77446e-01_jprb, 2.97390e-01_jprb, 3.18768e-01_jprb, &
     & 3.41682e-01_jprb, 3.66244e-01_jprb, 3.92571e-01_jprb, 4.20791e-01_jprb, 4.51039e-01_jprb, &
     & 4.83461e-01_jprb, 5.18215e-01_jprb, 5.55466e-01_jprb, 5.95396e-01_jprb, 6.38195e-01_jprb, &
     & 6.84071e-01_jprb, 7.33245e-01_jprb, 7.85954e-01_jprb, 8.42452e-01_jprb/)
      kao_mn2o( 6, :,10) = (/ &
     & 2.37173e-01_jprb, 2.54360e-01_jprb, 2.72792e-01_jprb, 2.92559e-01_jprb, 3.13759e-01_jprb, &
     & 3.36495e-01_jprb, 3.60878e-01_jprb, 3.87029e-01_jprb, 4.15074e-01_jprb, 4.45152e-01_jprb, &
     & 4.77409e-01_jprb, 5.12004e-01_jprb, 5.49105e-01_jprb, 5.88895e-01_jprb, 6.31569e-01_jprb, &
     & 6.77334e-01_jprb, 7.26416e-01_jprb, 7.79055e-01_jprb, 8.35508e-01_jprb/)
      kao_mn2o( 7, :,10) = (/ &
     & 2.27414e-01_jprb, 2.44231e-01_jprb, 2.62292e-01_jprb, 2.81689e-01_jprb, 3.02520e-01_jprb, &
     & 3.24891e-01_jprb, 3.48917e-01_jprb, 3.74720e-01_jprb, 4.02430e-01_jprb, 4.32190e-01_jprb, &
     & 4.64151e-01_jprb, 4.98475e-01_jprb, 5.35337e-01_jprb, 5.74926e-01_jprb, 6.17442e-01_jprb, &
     & 6.63102e-01_jprb, 7.12138e-01_jprb, 7.64801e-01_jprb, 8.21358e-01_jprb/)
      kao_mn2o( 8, :,10) = (/ &
     & 1.77234e-01_jprb, 1.92029e-01_jprb, 2.08060e-01_jprb, 2.25429e-01_jprb, 2.44248e-01_jprb, &
     & 2.64638e-01_jprb, 2.86730e-01_jprb, 3.10667e-01_jprb, 3.36601e-01_jprb, 3.64701e-01_jprb, &
     & 3.95147e-01_jprb, 4.28134e-01_jprb, 4.63875e-01_jprb, 5.02600e-01_jprb, 5.44557e-01_jprb, &
     & 5.90017e-01_jprb, 6.39272e-01_jprb, 6.92639e-01_jprb, 7.50461e-01_jprb/)
      kao_mn2o( 9, :,10) = (/ &
     & 2.41727e-01_jprb, 2.59094e-01_jprb, 2.77710e-01_jprb, 2.97662e-01_jprb, 3.19049e-01_jprb, &
     & 3.41972e-01_jprb, 3.66541e-01_jprb, 3.92877e-01_jprb, 4.21104e-01_jprb, 4.51359e-01_jprb, &
     & 4.83788e-01_jprb, 5.18547e-01_jprb, 5.55804e-01_jprb, 5.95737e-01_jprb, 6.38539e-01_jprb, &
     & 6.84417e-01_jprb, 7.33590e-01_jprb, 7.86297e-01_jprb, 8.42790e-01_jprb/)
      kao_mn2o( 1, :,11) = (/ &
     & 6.65287e+00_jprb, 6.69137e+00_jprb, 6.73008e+00_jprb, 6.76903e+00_jprb, 6.80820e+00_jprb, &
     & 6.84759e+00_jprb, 6.88721e+00_jprb, 6.92707e+00_jprb, 6.96715e+00_jprb, 7.00746e+00_jprb, &
     & 7.04801e+00_jprb, 7.08879e+00_jprb, 7.12981e+00_jprb, 7.17107e+00_jprb, 7.21256e+00_jprb, &
     & 7.25430e+00_jprb, 7.29628e+00_jprb, 7.33850e+00_jprb, 7.38096e+00_jprb/)
      kao_mn2o( 2, :,11) = (/ &
     & 2.06252e-01_jprb, 2.27731e-01_jprb, 2.51447e-01_jprb, 2.77633e-01_jprb, 3.06546e-01_jprb, &
     & 3.38470e-01_jprb, 3.73719e-01_jprb, 4.12638e-01_jprb, 4.55611e-01_jprb, 5.03058e-01_jprb, &
     & 5.55447e-01_jprb, 6.13292e-01_jprb, 6.77160e-01_jprb, 7.47680e-01_jprb, 8.25544e-01_jprb, &
     & 9.11517e-01_jprb, 1.00644e+00_jprb, 1.11125e+00_jprb, 1.22698e+00_jprb/)
      kao_mn2o( 3, :,11) = (/ &
     & 2.05840e-01_jprb, 2.27279e-01_jprb, 2.50952e-01_jprb, 2.77090e-01_jprb, 3.05950e-01_jprb, &
     & 3.37816e-01_jprb, 3.73002e-01_jprb, 4.11852e-01_jprb, 4.54748e-01_jprb, 5.02113e-01_jprb, &
     & 5.54411e-01_jprb, 6.12155e-01_jprb, 6.75915e-01_jprb, 7.46315e-01_jprb, 8.24047e-01_jprb, &
     & 9.09876e-01_jprb, 1.00465e+00_jprb, 1.10928e+00_jprb, 1.22482e+00_jprb/)
      kao_mn2o( 4, :,11) = (/ &
     & 2.04702e-01_jprb, 2.26031e-01_jprb, 2.49582e-01_jprb, 2.75587e-01_jprb, 3.04301e-01_jprb, &
     & 3.36008e-01_jprb, 3.71018e-01_jprb, 4.09675e-01_jprb, 4.52361e-01_jprb, 4.99495e-01_jprb, &
     & 5.51539e-01_jprb, 6.09007e-01_jprb, 6.72461e-01_jprb, 7.42528e-01_jprb, 8.19895e-01_jprb, &
     & 9.05323e-01_jprb, 9.99653e-01_jprb, 1.10381e+00_jprb, 1.21882e+00_jprb/)
      kao_mn2o( 5, :,11) = (/ &
     & 2.03481e-01_jprb, 2.24689e-01_jprb, 2.48108e-01_jprb, 2.73967e-01_jprb, 3.02522e-01_jprb, &
     & 3.34053e-01_jprb, 3.68871e-01_jprb, 4.07317e-01_jprb, 4.49771e-01_jprb, 4.96649e-01_jprb, &
     & 5.48414e-01_jprb, 6.05574e-01_jprb, 6.68691e-01_jprb, 7.38387e-01_jprb, 8.15347e-01_jprb, &
     & 9.00329e-01_jprb, 9.94168e-01_jprb, 1.09779e+00_jprb, 1.21221e+00_jprb/)
      kao_mn2o( 6, :,11) = (/ &
     & 2.01513e-01_jprb, 2.22529e-01_jprb, 2.45738e-01_jprb, 2.71367e-01_jprb, 2.99670e-01_jprb, &
     & 3.30924e-01_jprb, 3.65437e-01_jprb, 4.03550e-01_jprb, 4.45639e-01_jprb, 4.92117e-01_jprb, &
     & 5.43442e-01_jprb, 6.00120e-01_jprb, 6.62710e-01_jprb, 7.31827e-01_jprb, 8.08153e-01_jprb, &
     & 8.92439e-01_jprb, 9.85516e-01_jprb, 1.08830e+00_jprb, 1.20180e+00_jprb/)
      kao_mn2o( 7, :,11) = (/ &
     & 1.97136e-01_jprb, 2.17723e-01_jprb, 2.40459e-01_jprb, 2.65570e-01_jprb, 2.93304e-01_jprb, &
     & 3.23933e-01_jprb, 3.57762e-01_jprb, 3.95122e-01_jprb, 4.36385e-01_jprb, 4.81956e-01_jprb, &
     & 5.32287e-01_jprb, 5.87873e-01_jprb, 6.49264e-01_jprb, 7.17067e-01_jprb, 7.91950e-01_jprb, &
     & 8.74653e-01_jprb, 9.65992e-01_jprb, 1.06687e+00_jprb, 1.17828e+00_jprb/)
      kao_mn2o( 8, :,11) = (/ &
     & 1.79518e-01_jprb, 1.98371e-01_jprb, 2.19204e-01_jprb, 2.42224e-01_jprb, 2.67662e-01_jprb, &
     & 2.95772e-01_jprb, 3.26833e-01_jprb, 3.61157e-01_jprb, 3.99085e-01_jprb, 4.40996e-01_jprb, &
     & 4.87309e-01_jprb, 5.38486e-01_jprb, 5.95037e-01_jprb, 6.57526e-01_jprb, 7.26579e-01_jprb, &
     & 8.02883e-01_jprb, 8.87201e-01_jprb, 9.80373e-01_jprb, 1.08333e+00_jprb/)
      kao_mn2o( 9, :,11) = (/ &
     & 2.03481e-01_jprb, 2.24689e-01_jprb, 2.48108e-01_jprb, 2.73967e-01_jprb, 3.02522e-01_jprb, &
     & 3.34053e-01_jprb, 3.68871e-01_jprb, 4.07317e-01_jprb, 4.49771e-01_jprb, 4.96649e-01_jprb, &
     & 5.48414e-01_jprb, 6.05574e-01_jprb, 6.68691e-01_jprb, 7.38387e-01_jprb, 8.15347e-01_jprb, &
     & 9.00329e-01_jprb, 9.94168e-01_jprb, 1.09779e+00_jprb, 1.21221e+00_jprb/)
      kao_mn2o( 1, :,12) = (/ &
     & 5.89636e+00_jprb, 5.95081e+00_jprb, 6.00576e+00_jprb, 6.06121e+00_jprb, 6.11718e+00_jprb, &
     & 6.17366e+00_jprb, 6.23067e+00_jprb, 6.28820e+00_jprb, 6.34627e+00_jprb, 6.40487e+00_jprb, &
     & 6.46401e+00_jprb, 6.52369e+00_jprb, 6.58393e+00_jprb, 6.64472e+00_jprb, 6.70608e+00_jprb, &
     & 6.76800e+00_jprb, 6.83050e+00_jprb, 6.89357e+00_jprb, 6.95722e+00_jprb/)
      kao_mn2o( 2, :,12) = (/ &
     & 7.18699e-05_jprb, 9.48140e-05_jprb, 1.25083e-04_jprb, 1.65015e-04_jprb, 2.17695e-04_jprb, &
     & 2.87193e-04_jprb, 3.78877e-04_jprb, 4.99831e-04_jprb, 6.59400e-04_jprb, 8.69909e-04_jprb, &
     & 1.14762e-03_jprb, 1.51400e-03_jprb, 1.99733e-03_jprb, 2.63497e-03_jprb, 3.47616e-03_jprb, &
     & 4.58591e-03_jprb, 6.04993e-03_jprb, 7.98133e-03_jprb, 1.05293e-02_jprb/)
      kao_mn2o( 3, :,12) = (/ &
     & 7.20868e-05_jprb, 9.50993e-05_jprb, 1.25458e-04_jprb, 1.65508e-04_jprb, 2.18344e-04_jprb, &
     & 2.88046e-04_jprb, 3.80000e-04_jprb, 5.01307e-04_jprb, 6.61341e-04_jprb, 8.72462e-04_jprb, &
     & 1.15098e-03_jprb, 1.51841e-03_jprb, 2.00313e-03_jprb, 2.64260e-03_jprb, 3.48620e-03_jprb, &
     & 4.59911e-03_jprb, 6.06729e-03_jprb, 8.00416e-03_jprb, 1.05593e-02_jprb/)
      kao_mn2o( 4, :,12) = (/ &
     & 7.21734e-05_jprb, 9.52161e-05_jprb, 1.25616e-04_jprb, 1.65721e-04_jprb, 2.18630e-04_jprb, &
     & 2.88432e-04_jprb, 3.80519e-04_jprb, 5.02007e-04_jprb, 6.62282e-04_jprb, 8.73727e-04_jprb, &
     & 1.15268e-03_jprb, 1.52070e-03_jprb, 2.00621e-03_jprb, 2.64673e-03_jprb, 3.49174e-03_jprb, &
     & 4.60655e-03_jprb, 6.07727e-03_jprb, 8.01755e-03_jprb, 1.05773e-02_jprb/)
      kao_mn2o( 5, :,12) = (/ &
     & 7.22599e-05_jprb, 9.53329e-05_jprb, 1.25773e-04_jprb, 1.65933e-04_jprb, 2.18916e-04_jprb, &
     & 2.88818e-04_jprb, 3.81038e-04_jprb, 5.02706e-04_jprb, 6.63223e-04_jprb, 8.74992e-04_jprb, &
     & 1.15438e-03_jprb, 1.52298e-03_jprb, 2.00928e-03_jprb, 2.65085e-03_jprb, 3.49728e-03_jprb, &
     & 4.61398e-03_jprb, 6.08725e-03_jprb, 8.03094e-03_jprb, 1.05952e-02_jprb/)
      kao_mn2o( 6, :,12) = (/ &
     & 7.29962e-05_jprb, 9.63091e-05_jprb, 1.27067e-04_jprb, 1.67649e-04_jprb, 2.21191e-04_jprb, &
     & 2.91833e-04_jprb, 3.85036e-04_jprb, 5.08005e-04_jprb, 6.70247e-04_jprb, 8.84304e-04_jprb, &
     & 1.16672e-03_jprb, 1.53934e-03_jprb, 2.03096e-03_jprb, 2.67959e-03_jprb, 3.53537e-03_jprb, &
     & 4.66447e-03_jprb, 6.15417e-03_jprb, 8.11962e-03_jprb, 1.07128e-02_jprb/)
      kao_mn2o( 7, :,12) = (/ &
     & 7.47398e-05_jprb, 9.86139e-05_jprb, 1.30114e-04_jprb, 1.71677e-04_jprb, 2.26516e-04_jprb, &
     & 2.98872e-04_jprb, 3.94340e-04_jprb, 5.20305e-04_jprb, 6.86506e-04_jprb, 9.05797e-04_jprb, &
     & 1.19514e-03_jprb, 1.57690e-03_jprb, 2.08061e-03_jprb, 2.74522e-03_jprb, 3.62213e-03_jprb, &
     & 4.77915e-03_jprb, 6.30576e-03_jprb, 8.32001e-03_jprb, 1.09777e-02_jprb/)
      kao_mn2o( 8, :,12) = (/ &
     & 7.57487e-05_jprb, 9.99802e-05_jprb, 1.31963e-04_jprb, 1.74177e-04_jprb, 2.29896e-04_jprb, &
     & 3.03438e-04_jprb, 4.00506e-04_jprb, 5.28625e-04_jprb, 6.97729e-04_jprb, 9.20927e-04_jprb, &
     & 1.21553e-03_jprb, 1.60437e-03_jprb, 2.11759e-03_jprb, 2.79499e-03_jprb, 3.68909e-03_jprb, &
     & 4.86921e-03_jprb, 6.42684e-03_jprb, 8.48274e-03_jprb, 1.11963e-02_jprb/)
      kao_mn2o( 9, :,12) = (/ &
     & 7.22467e-05_jprb, 9.53177e-05_jprb, 1.25756e-04_jprb, 1.65915e-04_jprb, 2.18898e-04_jprb, &
     & 2.88800e-04_jprb, 3.81024e-04_jprb, 5.02700e-04_jprb, 6.63231e-04_jprb, 8.75024e-04_jprb, &
     & 1.15445e-03_jprb, 1.52311e-03_jprb, 2.00950e-03_jprb, 2.65121e-03_jprb, 3.49784e-03_jprb, &
     & 4.61483e-03_jprb, 6.08852e-03_jprb, 8.03280e-03_jprb, 1.05980e-02_jprb/)
      kao_mn2o( 1, :,13) = (/ &
     & 1.14265e+01_jprb, 1.16380e+01_jprb, 1.18534e+01_jprb, 1.20728e+01_jprb, 1.22962e+01_jprb, &
     & 1.25238e+01_jprb, 1.27556e+01_jprb, 1.29917e+01_jprb, 1.32322e+01_jprb, 1.34771e+01_jprb, &
     & 1.37265e+01_jprb, 1.39806e+01_jprb, 1.42394e+01_jprb, 1.45029e+01_jprb, 1.47714e+01_jprb, &
     & 1.50448e+01_jprb, 1.53232e+01_jprb, 1.56068e+01_jprb, 1.58957e+01_jprb/)
      kao_mn2o( 2, :,13) = (/ &
     & 7.97796e-05_jprb, 1.05659e-04_jprb, 1.39932e-04_jprb, 1.85324e-04_jprb, 2.45439e-04_jprb, &
     & 3.25054e-04_jprb, 4.30496e-04_jprb, 5.70140e-04_jprb, 7.55082e-04_jprb, 1.00002e-03_jprb, &
     & 1.32440e-03_jprb, 1.75401e-03_jprb, 2.32298e-03_jprb, 3.07651e-03_jprb, 4.07447e-03_jprb, &
     & 5.39614e-03_jprb, 7.14655e-03_jprb, 9.46475e-03_jprb, 1.25349e-02_jprb/)
      kao_mn2o( 3, :,13) = (/ &
     & 7.95035e-05_jprb, 1.05293e-04_jprb, 1.39449e-04_jprb, 1.84684e-04_jprb, 2.44592e-04_jprb, &
     & 3.23934e-04_jprb, 4.29013e-04_jprb, 5.68178e-04_jprb, 7.52486e-04_jprb, 9.96580e-04_jprb, &
     & 1.31985e-03_jprb, 1.74800e-03_jprb, 2.31502e-03_jprb, 3.06597e-03_jprb, 4.06052e-03_jprb, &
     & 5.37770e-03_jprb, 7.12214e-03_jprb, 9.43244e-03_jprb, 1.24922e-02_jprb/)
      kao_mn2o( 4, :,13) = (/ &
     & 7.92339e-05_jprb, 1.04938e-04_jprb, 1.38980e-04_jprb, 1.84065e-04_jprb, 2.43776e-04_jprb, &
     & 3.22857e-04_jprb, 4.27593e-04_jprb, 5.66305e-04_jprb, 7.50016e-04_jprb, 9.93322e-04_jprb, &
     & 1.31556e-03_jprb, 1.74233e-03_jprb, 2.30754e-03_jprb, 3.05612e-03_jprb, 4.04752e-03_jprb, &
     & 5.36055e-03_jprb, 7.09953e-03_jprb, 9.40262e-03_jprb, 1.24528e-02_jprb/)
      kao_mn2o( 5, :,13) = (/ &
     & 7.90000e-05_jprb, 1.04627e-04_jprb, 1.38566e-04_jprb, 1.83516e-04_jprb, 2.43046e-04_jprb, &
     & 3.21887e-04_jprb, 4.26303e-04_jprb, 5.64591e-04_jprb, 7.47738e-04_jprb, 9.90295e-04_jprb, &
     & 1.31154e-03_jprb, 1.73698e-03_jprb, 2.30044e-03_jprb, 3.04667e-03_jprb, 4.03498e-03_jprb, &
     & 5.34388e-03_jprb, 7.07737e-03_jprb, 9.37318e-03_jprb, 1.24137e-02_jprb/)
      kao_mn2o( 6, :,13) = (/ &
     & 7.76004e-05_jprb, 1.02776e-04_jprb, 1.36118e-04_jprb, 1.80278e-04_jprb, 2.38764e-04_jprb, &
     & 3.16224e-04_jprb, 4.18814e-04_jprb, 5.54686e-04_jprb, 7.34638e-04_jprb, 9.72970e-04_jprb, &
     & 1.28862e-03_jprb, 1.70668e-03_jprb, 2.26036e-03_jprb, 2.99367e-03_jprb, 3.96488e-03_jprb, &
     & 5.25118e-03_jprb, 6.95477e-03_jprb, 9.21105e-03_jprb, 1.21993e-02_jprb/)
      kao_mn2o( 7, :,13) = (/ &
     & 7.52813e-05_jprb, 9.97094e-05_jprb, 1.32064e-04_jprb, 1.74918e-04_jprb, 2.31677e-04_jprb, &
     & 3.06854e-04_jprb, 4.06426e-04_jprb, 5.38308e-04_jprb, 7.12984e-04_jprb, 9.44341e-04_jprb, &
     & 1.25077e-03_jprb, 1.65664e-03_jprb, 2.19420e-03_jprb, 2.90620e-03_jprb, 3.84923e-03_jprb, &
     & 5.09828e-03_jprb, 6.75263e-03_jprb, 8.94379e-03_jprb, 1.18460e-02_jprb/)
      kao_mn2o( 8, :,13) = (/ &
     & 6.87436e-05_jprb, 9.10605e-05_jprb, 1.20622e-04_jprb, 1.59781e-04_jprb, 2.11653e-04_jprb, &
     & 2.80364e-04_jprb, 3.71381e-04_jprb, 4.91946e-04_jprb, 6.51651e-04_jprb, 8.63203e-04_jprb, &
     & 1.14343e-03_jprb, 1.51464e-03_jprb, 2.00635e-03_jprb, 2.65769e-03_jprb, 3.52048e-03_jprb, &
     & 4.66337e-03_jprb, 6.17729e-03_jprb, 8.18269e-03_jprb, 1.08391e-02_jprb/)
      kao_mn2o( 9, :,13) = (/ &
     & 7.90357e-05_jprb, 1.04671e-04_jprb, 1.38622e-04_jprb, 1.83585e-04_jprb, 2.43132e-04_jprb, &
     & 3.21994e-04_jprb, 4.26434e-04_jprb, 5.64750e-04_jprb, 7.47931e-04_jprb, 9.90526e-04_jprb, &
     & 1.31181e-03_jprb, 1.73730e-03_jprb, 2.30081e-03_jprb, 3.04709e-03_jprb, 4.03543e-03_jprb, &
     & 5.34435e-03_jprb, 7.07782e-03_jprb, 9.37355e-03_jprb, 1.24139e-02_jprb/)
      kao_mn2o( 1, :,14) = (/ &
     & 1.61373e+01_jprb, 1.64784e+01_jprb, 1.68266e+01_jprb, 1.71822e+01_jprb, 1.75454e+01_jprb, &
     & 1.79162e+01_jprb, 1.82948e+01_jprb, 1.86814e+01_jprb, 1.90762e+01_jprb, 1.94794e+01_jprb, &
     & 1.98911e+01_jprb, 2.03114e+01_jprb, 2.07407e+01_jprb, 2.11790e+01_jprb, 2.16266e+01_jprb, &
     & 2.20836e+01_jprb, 2.25504e+01_jprb, 2.30269e+01_jprb, 2.35136e+01_jprb/)
      kao_mn2o( 2, :,14) = (/ &
     & 6.92866e-10_jprb, 9.24655e-10_jprb, 1.23398e-09_jprb, 1.64680e-09_jprb, 2.19771e-09_jprb, &
     & 2.93292e-09_jprb, 3.91409e-09_jprb, 5.22349e-09_jprb, 6.97093e-09_jprb, 9.30295e-09_jprb, &
     & 1.24151e-08_jprb, 1.65684e-08_jprb, 2.21111e-08_jprb, 2.95081e-08_jprb, 3.93796e-08_jprb, &
     & 5.25535e-08_jprb, 7.01346e-08_jprb, 9.35970e-08_jprb, 1.24908e-07_jprb/)
      kao_mn2o( 3, :,14) = (/ &
     & 6.94564e-10_jprb, 9.26928e-10_jprb, 1.23703e-09_jprb, 1.65088e-09_jprb, 2.20317e-09_jprb, &
     & 2.94024e-09_jprb, 3.92389e-09_jprb, 5.23661e-09_jprb, 6.98851e-09_jprb, 9.32650e-09_jprb, &
     & 1.24467e-08_jprb, 1.66107e-08_jprb, 2.21677e-08_jprb, 2.95839e-08_jprb, 3.94811e-08_jprb, &
     & 5.26894e-08_jprb, 7.03165e-08_jprb, 9.38407e-08_jprb, 1.25235e-07_jprb/)
      kao_mn2o( 4, :,14) = (/ &
     & 6.98644e-10_jprb, 9.32310e-10_jprb, 1.24413e-09_jprb, 1.66023e-09_jprb, 2.21551e-09_jprb, &
     & 2.95649e-09_jprb, 3.94531e-09_jprb, 5.26484e-09_jprb, 7.02570e-09_jprb, 9.37548e-09_jprb, &
     & 1.25112e-08_jprb, 1.66956e-08_jprb, 2.22795e-08_jprb, 2.97311e-08_jprb, 3.96748e-08_jprb, &
     & 5.29443e-08_jprb, 7.06518e-08_jprb, 9.42817e-08_jprb, 1.25815e-07_jprb/)
      kao_mn2o( 5, :,14) = (/ &
     & 7.03261e-10_jprb, 9.38472e-10_jprb, 1.25235e-09_jprb, 1.67121e-09_jprb, 2.23016e-09_jprb, &
     & 2.97605e-09_jprb, 3.97141e-09_jprb, 5.29968e-09_jprb, 7.07220e-09_jprb, 9.43754e-09_jprb, &
     & 1.25940e-08_jprb, 1.68062e-08_jprb, 2.24271e-08_jprb, 2.99280e-08_jprb, 3.99376e-08_jprb, &
     & 5.32951e-08_jprb, 7.11200e-08_jprb, 9.49066e-08_jprb, 1.26649e-07_jprb/)
      kao_mn2o( 6, :,14) = (/ &
     & 7.12478e-10_jprb, 9.50674e-10_jprb, 1.26850e-09_jprb, 1.69259e-09_jprb, 2.25845e-09_jprb, &
     & 3.01350e-09_jprb, 4.02096e-09_jprb, 5.36525e-09_jprb, 7.15896e-09_jprb, 9.55233e-09_jprb, &
     & 1.27459e-08_jprb, 1.70071e-08_jprb, 2.26928e-08_jprb, 3.02795e-08_jprb, 4.04025e-08_jprb, &
     & 5.39099e-08_jprb, 7.19330e-08_jprb, 9.59815e-08_jprb, 1.28070e-07_jprb/)
      kao_mn2o( 7, :,14) = (/ &
     & 7.28994e-10_jprb, 9.72644e-10_jprb, 1.29773e-09_jprb, 1.73147e-09_jprb, 2.31017e-09_jprb, &
     & 3.08230e-09_jprb, 4.11249e-09_jprb, 5.48700e-09_jprb, 7.32092e-09_jprb, 9.76777e-09_jprb, &
     & 1.30324e-08_jprb, 1.73883e-08_jprb, 2.31999e-08_jprb, 3.09540e-08_jprb, 4.12996e-08_jprb, &
     & 5.51032e-08_jprb, 7.35203e-08_jprb, 9.80928e-08_jprb, 1.30878e-07_jprb/)
      kao_mn2o( 8, :,14) = (/ &
     & 7.87604e-10_jprb, 1.05043e-09_jprb, 1.40097e-09_jprb, 1.86848e-09_jprb, 2.49201e-09_jprb, &
     & 3.32360e-09_jprb, 4.43271e-09_jprb, 5.91194e-09_jprb, 7.88479e-09_jprb, 1.05160e-08_jprb, &
     & 1.40253e-08_jprb, 1.87056e-08_jprb, 2.49478e-08_jprb, 3.32730e-08_jprb, 4.43764e-08_jprb, &
     & 5.91851e-08_jprb, 7.89356e-08_jprb, 1.05277e-07_jprb, 1.40408e-07_jprb/)
      kao_mn2o( 9, :,14) = (/ &
     & 7.03261e-10_jprb, 9.38472e-10_jprb, 1.25235e-09_jprb, 1.67121e-09_jprb, 2.23016e-09_jprb, &
     & 2.97605e-09_jprb, 3.97141e-09_jprb, 5.29968e-09_jprb, 7.07220e-09_jprb, 9.43754e-09_jprb, &
     & 1.25940e-08_jprb, 1.68062e-08_jprb, 2.24271e-08_jprb, 2.99280e-08_jprb, 3.99376e-08_jprb, &
     & 5.32951e-08_jprb, 7.11200e-08_jprb, 9.49066e-08_jprb, 1.26649e-07_jprb/)
      kao_mn2o( 1, :,15) = (/ &
     & 2.14029e+01_jprb, 2.16782e+01_jprb, 2.19571e+01_jprb, 2.22396e+01_jprb, 2.25257e+01_jprb, &
     & 2.28155e+01_jprb, 2.31090e+01_jprb, 2.34063e+01_jprb, 2.37074e+01_jprb, 2.40124e+01_jprb, &
     & 2.43213e+01_jprb, 2.46342e+01_jprb, 2.49511e+01_jprb, 2.52721e+01_jprb, 2.55972e+01_jprb, &
     & 2.59265e+01_jprb, 2.62600e+01_jprb, 2.65979e+01_jprb, 2.69400e+01_jprb/)
      kao_mn2o( 2, :,15) = (/ &
     & 5.68659e-10_jprb, 7.55629e-10_jprb, 1.00407e-09_jprb, 1.33421e-09_jprb, 1.77288e-09_jprb, &
     & 2.35579e-09_jprb, 3.13036e-09_jprb, 4.15960e-09_jprb, 5.52724e-09_jprb, 7.34455e-09_jprb, &
     & 9.75939e-09_jprb, 1.29682e-08_jprb, 1.72320e-08_jprb, 2.28978e-08_jprb, 3.04264e-08_jprb, &
     & 4.04304e-08_jprb, 5.37236e-08_jprb, 7.13875e-08_jprb, 9.48591e-08_jprb/)
      kao_mn2o( 3, :,15) = (/ &
     & 5.59573e-10_jprb, 7.43558e-10_jprb, 9.88035e-10_jprb, 1.31290e-09_jprb, 1.74457e-09_jprb, &
     & 2.31817e-09_jprb, 3.08037e-09_jprb, 4.09318e-09_jprb, 5.43900e-09_jprb, 7.22730e-09_jprb, &
     & 9.60360e-09_jprb, 1.27612e-08_jprb, 1.69570e-08_jprb, 2.25324e-08_jprb, 2.99409e-08_jprb, &
     & 3.97853e-08_jprb, 5.28665e-08_jprb, 7.02486e-08_jprb, 9.33459e-08_jprb/)
      kao_mn2o( 4, :,15) = (/ &
     & 5.50488e-10_jprb, 7.31486e-10_jprb, 9.71996e-10_jprb, 1.29158e-09_jprb, 1.71625e-09_jprb, &
     & 2.28055e-09_jprb, 3.03039e-09_jprb, 4.02676e-09_jprb, 5.35075e-09_jprb, 7.11005e-09_jprb, &
     & 9.44781e-09_jprb, 1.25542e-08_jprb, 1.66820e-08_jprb, 2.21670e-08_jprb, 2.94554e-08_jprb, &
     & 3.91402e-08_jprb, 5.20093e-08_jprb, 6.91098e-08_jprb, 9.18327e-08_jprb/)
      kao_mn2o( 5, :,15) = (/ &
     & 5.34010e-10_jprb, 7.09574e-10_jprb, 9.42858e-10_jprb, 1.25284e-09_jprb, 1.66473e-09_jprb, &
     & 2.21203e-09_jprb, 2.93927e-09_jprb, 3.90560e-09_jprb, 5.18963e-09_jprb, 6.89580e-09_jprb, &
     & 9.16290e-09_jprb, 1.21754e-08_jprb, 1.61782e-08_jprb, 2.14970e-08_jprb, 2.85645e-08_jprb, &
     & 3.79555e-08_jprb, 5.04340e-08_jprb, 6.70149e-08_jprb, 8.90470e-08_jprb/)
      kao_mn2o( 6, :,15) = (/ &
     & 5.08144e-10_jprb, 6.75221e-10_jprb, 8.97231e-10_jprb, 1.19224e-09_jprb, 1.58424e-09_jprb, &
     & 2.10513e-09_jprb, 2.79729e-09_jprb, 3.71703e-09_jprb, 4.93919e-09_jprb, 6.56317e-09_jprb, &
     & 8.72112e-09_jprb, 1.15886e-08_jprb, 1.53989e-08_jprb, 2.04620e-08_jprb, 2.71898e-08_jprb, &
     & 3.61297e-08_jprb, 4.80091e-08_jprb, 6.37943e-08_jprb, 8.47696e-08_jprb/)
      kao_mn2o( 7, :,15) = (/ &
     & 4.56716e-10_jprb, 6.06884e-10_jprb, 8.06427e-10_jprb, 1.07158e-09_jprb, 1.42391e-09_jprb, &
     & 1.89210e-09_jprb, 2.51422e-09_jprb, 3.34089e-09_jprb, 4.43938e-09_jprb, 5.89904e-09_jprb, &
     & 7.83864e-09_jprb, 1.04160e-08_jprb, 1.38408e-08_jprb, 1.83916e-08_jprb, 2.44387e-08_jprb, &
     & 3.24742e-08_jprb, 4.31517e-08_jprb, 5.73399e-08_jprb, 7.61932e-08_jprb/)
      kao_mn2o( 8, :,15) = (/ &
     & 2.78366e-10_jprb, 3.69881e-10_jprb, 4.91482e-10_jprb, 6.53061e-10_jprb, 8.67760e-10_jprb, &
     & 1.15304e-09_jprb, 1.53211e-09_jprb, 2.03581e-09_jprb, 2.70510e-09_jprb, 3.59441e-09_jprb, &
     & 4.77611e-09_jprb, 6.34629e-09_jprb, 8.43268e-09_jprb, 1.12050e-08_jprb, 1.48887e-08_jprb, &
     & 1.97835e-08_jprb, 2.62875e-08_jprb, 3.49296e-08_jprb, 4.64130e-08_jprb/)
      kao_mn2o( 9, :,15) = (/ &
     & 5.34010e-10_jprb, 7.09574e-10_jprb, 9.42858e-10_jprb, 1.25284e-09_jprb, 1.66473e-09_jprb, &
     & 2.21203e-09_jprb, 2.93927e-09_jprb, 3.90560e-09_jprb, 5.18963e-09_jprb, 6.89580e-09_jprb, &
     & 9.16290e-09_jprb, 1.21754e-08_jprb, 1.61782e-08_jprb, 2.14970e-08_jprb, 2.85645e-08_jprb, &
     & 3.79555e-08_jprb, 5.04340e-08_jprb, 6.70149e-08_jprb, 8.90470e-08_jprb/)
      kao_mn2o( 1, :,16) = (/ &
     & 2.90784e+01_jprb, 2.93787e+01_jprb, 2.96820e+01_jprb, 2.99885e+01_jprb, 3.02982e+01_jprb, &
     & 3.06110e+01_jprb, 3.09271e+01_jprb, 3.12464e+01_jprb, 3.15690e+01_jprb, 3.18950e+01_jprb, &
     & 3.22243e+01_jprb, 3.25571e+01_jprb, 3.28932e+01_jprb, 3.32329e+01_jprb, 3.35760e+01_jprb, &
     & 3.39227e+01_jprb, 3.42730e+01_jprb, 3.46269e+01_jprb, 3.49844e+01_jprb/)
      kao_mn2o( 2, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 3, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 4, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 5, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 6, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 7, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 8, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)
      kao_mn2o( 9, :,16) = (/ &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, &
     & 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb, 0.00000e+00_jprb/)

!     the array kbo_mxx contains the absorption coefficient for 
!     a minor species at the 16 chosen g-values for a reference pressure
!     level above 100~ mb.   the first index refers to temperature 
!     in 7.2 degree increments.  for instance, jt = 1 refers to a 
!     temperature of 188.0, jt = 2 refers to 195.2, etc. the second index 
!     runs over the g-channel (1 to 16).

      kbo_mn2o(:, 1) = (/ &
     & 8.42688e-03_jprb, 8.96787e-03_jprb, 9.54358e-03_jprb, 1.01563e-02_jprb, 1.08083e-02_jprb, &
     & 1.15021e-02_jprb, 1.22405e-02_jprb, 1.30263e-02_jprb, 1.38626e-02_jprb, 1.47525e-02_jprb, &
     & 1.56996e-02_jprb, 1.67075e-02_jprb, 1.77800e-02_jprb, 1.89215e-02_jprb, 2.01362e-02_jprb, &
     & 2.14289e-02_jprb, 2.28045e-02_jprb, 2.42685e-02_jprb, 2.58265e-02_jprb/)
      kbo_mn2o(:, 2) = (/ &
     & 2.24976e-02_jprb, 2.38935e-02_jprb, 2.53762e-02_jprb, 2.69508e-02_jprb, 2.86231e-02_jprb, &
     & 3.03991e-02_jprb, 3.22854e-02_jprb, 3.42887e-02_jprb, 3.64163e-02_jprb, 3.86760e-02_jprb, &
     & 4.10759e-02_jprb, 4.36246e-02_jprb, 4.63315e-02_jprb, 4.92064e-02_jprb, 5.22597e-02_jprb, &
     & 5.55024e-02_jprb, 5.89464e-02_jprb, 6.26040e-02_jprb, 6.64886e-02_jprb/)
      kbo_mn2o(:, 3) = (/ &
     & 5.93542e-02_jprb, 6.37312e-02_jprb, 6.84310e-02_jprb, 7.34774e-02_jprb, 7.88960e-02_jprb, &
     & 8.47141e-02_jprb, 9.09613e-02_jprb, 9.76692e-02_jprb, 1.04872e-01_jprb, 1.12605e-01_jprb, &
     & 1.20910e-01_jprb, 1.29826e-01_jprb, 1.39400e-01_jprb, 1.49680e-01_jprb, 1.60718e-01_jprb, &
     & 1.72570e-01_jprb, 1.85296e-01_jprb, 1.98961e-01_jprb, 2.13633e-01_jprb/)
      kbo_mn2o(:, 4) = (/ &
     & 1.98022e-01_jprb, 2.05895e-01_jprb, 2.14082e-01_jprb, 2.22594e-01_jprb, 2.31445e-01_jprb, &
     & 2.40647e-01_jprb, 2.50216e-01_jprb, 2.60164e-01_jprb, 2.70509e-01_jprb, 2.81265e-01_jprb, &
     & 2.92448e-01_jprb, 3.04076e-01_jprb, 3.16167e-01_jprb, 3.28738e-01_jprb, 3.41809e-01_jprb, &
     & 3.55400e-01_jprb, 3.69531e-01_jprb, 3.84224e-01_jprb, 3.99501e-01_jprb/)
      kbo_mn2o(:, 5) = (/ &
     & 6.41413e-01_jprb, 6.46239e-01_jprb, 6.51101e-01_jprb, 6.56000e-01_jprb, 6.60936e-01_jprb, &
     & 6.65910e-01_jprb, 6.70920e-01_jprb, 6.75968e-01_jprb, 6.81054e-01_jprb, 6.86179e-01_jprb, &
     & 6.91342e-01_jprb, 6.96544e-01_jprb, 7.01785e-01_jprb, 7.07065e-01_jprb, 7.12385e-01_jprb, &
     & 7.17746e-01_jprb, 7.23146e-01_jprb, 7.28587e-01_jprb, 7.34070e-01_jprb/)
      kbo_mn2o(:, 6) = (/ &
     & 1.47906e+00_jprb, 1.48768e+00_jprb, 1.49635e+00_jprb, 1.50507e+00_jprb, 1.51384e+00_jprb, &
     & 1.52267e+00_jprb, 1.53154e+00_jprb, 1.54047e+00_jprb, 1.54944e+00_jprb, 1.55847e+00_jprb, &
     & 1.56755e+00_jprb, 1.57669e+00_jprb, 1.58588e+00_jprb, 1.59512e+00_jprb, 1.60442e+00_jprb, &
     & 1.61377e+00_jprb, 1.62317e+00_jprb, 1.63263e+00_jprb, 1.64215e+00_jprb/)
      kbo_mn2o(:, 7) = (/ &
     & 3.53152e+00_jprb, 3.55492e+00_jprb, 3.57848e+00_jprb, 3.60219e+00_jprb, 3.62606e+00_jprb, &
     & 3.65008e+00_jprb, 3.67427e+00_jprb, 3.69862e+00_jprb, 3.72313e+00_jprb, 3.74780e+00_jprb, &
     & 3.77263e+00_jprb, 3.79763e+00_jprb, 3.82279e+00_jprb, 3.84812e+00_jprb, 3.87362e+00_jprb, &
     & 3.89929e+00_jprb, 3.92513e+00_jprb, 3.95114e+00_jprb, 3.97732e+00_jprb/)
      kbo_mn2o(:, 8) = (/ &
     & 9.06783e+00_jprb, 9.04597e+00_jprb, 9.02415e+00_jprb, 9.00239e+00_jprb, 8.98069e+00_jprb, &
     & 8.95903e+00_jprb, 8.93743e+00_jprb, 8.91588e+00_jprb, 8.89438e+00_jprb, 8.87293e+00_jprb, &
     & 8.85154e+00_jprb, 8.83020e+00_jprb, 8.80890e+00_jprb, 8.78766e+00_jprb, 8.76647e+00_jprb, &
     & 8.74533e+00_jprb, 8.72425e+00_jprb, 8.70321e+00_jprb, 8.68223e+00_jprb/)
      kbo_mn2o(:, 9) = (/ &
     & 3.88220e+01_jprb, 3.85805e+01_jprb, 3.83405e+01_jprb, 3.81019e+01_jprb, 3.78649e+01_jprb, &
     & 3.76293e+01_jprb, 3.73952e+01_jprb, 3.71625e+01_jprb, 3.69313e+01_jprb, 3.67016e+01_jprb, &
     & 3.64732e+01_jprb, 3.62463e+01_jprb, 3.60208e+01_jprb, 3.57967e+01_jprb, 3.55740e+01_jprb, &
     & 3.53527e+01_jprb, 3.51327e+01_jprb, 3.49142e+01_jprb, 3.46970e+01_jprb/)
      kbo_mn2o(:, 10) = (/ &
     & 1.14211e+02_jprb, 1.13955e+02_jprb, 1.13700e+02_jprb, 1.13445e+02_jprb, 1.13191e+02_jprb, &
     & 1.12938e+02_jprb, 1.12685e+02_jprb, 1.12433e+02_jprb, 1.12181e+02_jprb, 1.11930e+02_jprb, &
     & 1.11679e+02_jprb, 1.11429e+02_jprb, 1.11180e+02_jprb, 1.10931e+02_jprb, 1.10682e+02_jprb, &
     & 1.10434e+02_jprb, 1.10187e+02_jprb, 1.09940e+02_jprb, 1.09694e+02_jprb/)
      kbo_mn2o(:, 11) = (/ &
     & 1.60513e+02_jprb, 1.60857e+02_jprb, 1.61201e+02_jprb, 1.61547e+02_jprb, 1.61893e+02_jprb, &
     & 1.62240e+02_jprb, 1.62587e+02_jprb, 1.62936e+02_jprb, 1.63285e+02_jprb, 1.63635e+02_jprb, &
     & 1.63985e+02_jprb, 1.64337e+02_jprb, 1.64689e+02_jprb, 1.65041e+02_jprb, 1.65395e+02_jprb, &
     & 1.65749e+02_jprb, 1.66105e+02_jprb, 1.66460e+02_jprb, 1.66817e+02_jprb/)
      kbo_mn2o(:, 12) = (/ &
     & 1.71473e+02_jprb, 1.72766e+02_jprb, 1.74068e+02_jprb, 1.75381e+02_jprb, 1.76703e+02_jprb, &
     & 1.78035e+02_jprb, 1.79377e+02_jprb, 1.80729e+02_jprb, 1.82091e+02_jprb, 1.83464e+02_jprb, &
     & 1.84847e+02_jprb, 1.86240e+02_jprb, 1.87644e+02_jprb, 1.89059e+02_jprb, 1.90484e+02_jprb, &
     & 1.91920e+02_jprb, 1.93367e+02_jprb, 1.94824e+02_jprb, 1.96293e+02_jprb/)
      kbo_mn2o(:, 13) = (/ &
     & 2.71287e+01_jprb, 2.75538e+01_jprb, 2.79856e+01_jprb, 2.84241e+01_jprb, 2.88695e+01_jprb, &
     & 2.93219e+01_jprb, 2.97814e+01_jprb, 3.02480e+01_jprb, 3.07220e+01_jprb, 3.12035e+01_jprb, &
     & 3.16924e+01_jprb, 3.21890e+01_jprb, 3.26934e+01_jprb, 3.32058e+01_jprb, 3.37261e+01_jprb, &
     & 3.42546e+01_jprb, 3.47914e+01_jprb, 3.53365e+01_jprb, 3.58903e+01_jprb/)
      kbo_mn2o(:, 14) = (/ &
     & 1.70389e+01_jprb, 1.70899e+01_jprb, 1.71411e+01_jprb, 1.71924e+01_jprb, 1.72439e+01_jprb, &
     & 1.72955e+01_jprb, 1.73473e+01_jprb, 1.73992e+01_jprb, 1.74513e+01_jprb, 1.75035e+01_jprb, &
     & 1.75559e+01_jprb, 1.76085e+01_jprb, 1.76612e+01_jprb, 1.77141e+01_jprb, 1.77671e+01_jprb, &
     & 1.78203e+01_jprb, 1.78736e+01_jprb, 1.79271e+01_jprb, 1.79808e+01_jprb/)
      kbo_mn2o(:, 15) = (/ &
     & 2.49725e+00_jprb, 2.66861e+00_jprb, 2.85174e+00_jprb, 3.04743e+00_jprb, 3.25655e+00_jprb, &
     & 3.48003e+00_jprb, 3.71883e+00_jprb, 3.97403e+00_jprb, 4.24673e+00_jprb, 4.53815e+00_jprb, &
     & 4.84957e+00_jprb, 5.18236e+00_jprb, 5.53798e+00_jprb, 5.91801e+00_jprb, 6.32412e+00_jprb, &
     & 6.75809e+00_jprb, 7.22185e+00_jprb, 7.71742e+00_jprb, 8.24701e+00_jprb/)
      kbo_mn2o(:, 16) = (/ &
     & 1.82935e-03_jprb, 2.58912e-03_jprb, 3.66444e-03_jprb, 5.18637e-03_jprb, 7.34039e-03_jprb, &
     & 1.03890e-02_jprb, 1.47038e-02_jprb, 2.08106e-02_jprb, 2.94538e-02_jprb, 4.16865e-02_jprb, &
     & 5.89999e-02_jprb, 8.35040e-02_jprb, 1.18185e-01_jprb, 1.67270e-01_jprb, 2.36741e-01_jprb, &
     & 3.35065e-01_jprb, 4.74225e-01_jprb, 6.71180e-01_jprb, 9.49936e-01_jprb/)

!     the array forrefo contains the coefficient of the water vapor
!     foreign-continuum (including the energy term).  the first 
!     index refers to reference temperature (296,260,224,260) and 
!     pressure (970,475,219,3 mbar) levels.  the second index 
!     runs over the g-channel (1 to 16).

      forrefo(1,:) = (/ &
     &7.5352e-06_jprb,2.9812e-05_jprb,1.4497e-04_jprb,4.4006e-04_jprb,1.0492e-03_jprb,1.9676e-03_jprb, &
     &1.9989e-03_jprb,1.9099e-03_jprb,2.2121e-03_jprb,2.4491e-03_jprb,2.9573e-03_jprb,2.6344e-03_jprb, &
     &3.0629e-03_jprb,3.3547e-03_jprb,5.0643e-03_jprb,5.0642e-03_jprb/)
      forrefo(2,:) = (/ &
     &6.6070e-06_jprb,4.8618e-05_jprb,3.1112e-04_jprb,8.4235e-04_jprb,1.4179e-03_jprb,1.4315e-03_jprb, &
     &1.4685e-03_jprb,1.6554e-03_jprb,2.1171e-03_jprb,2.3545e-03_jprb,2.5165e-03_jprb,2.7680e-03_jprb, &
     &2.6985e-03_jprb,3.5345e-03_jprb,4.2924e-03_jprb,5.0712e-03_jprb/)
      forrefo(3,:) = (/ &
     &6.5962e-06_jprb,7.2595e-04_jprb,1.3429e-03_jprb,1.1675e-03_jprb,9.8384e-04_jprb,8.8787e-04_jprb, &
     &8.7557e-04_jprb,8.0589e-04_jprb,7.7024e-04_jprb,8.7518e-04_jprb,9.5213e-04_jprb,9.0849e-04_jprb, &
     &1.2596e-03_jprb,2.5106e-03_jprb,3.9471e-03_jprb,5.0742e-03_jprb/)
      forrefo(4,:) = (/ &
     &3.6217e-04_jprb,1.0709e-03_jprb,1.0628e-03_jprb,8.5640e-04_jprb,8.9332e-04_jprb,8.3372e-04_jprb, &
     &7.8539e-04_jprb,8.2828e-04_jprb,8.3329e-04_jprb,8.5118e-04_jprb,8.2878e-04_jprb,6.8570e-04_jprb, &
     &6.3815e-04_jprb,8.0648e-04_jprb,2.3236e-03_jprb,4.0321e-03_jprb/)


!     the array selfrefo contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

     selfrefo(:, 1) = (/ &
     & 2.83453e-02_jprb, 2.51439e-02_jprb, 2.23040e-02_jprb, 1.97849e-02_jprb, 1.75503e-02_jprb, &
     & 1.55681e-02_jprb, 1.38097e-02_jprb, 1.22500e-02_jprb, 1.08664e-02_jprb, 9.63912e-03_jprb/)
      selfrefo(:, 2) = (/ &
     & 3.05185e-02_jprb, 2.72374e-02_jprb, 2.43090e-02_jprb, 2.16955e-02_jprb, 1.93629e-02_jprb, &
     & 1.72811e-02_jprb, 1.54232e-02_jprb, 1.37650e-02_jprb, 1.22851e-02_jprb, 1.09643e-02_jprb/)
      selfrefo(:, 3) = (/ &
     & 4.23833e-02_jprb, 3.76250e-02_jprb, 3.34010e-02_jprb, 2.96512e-02_jprb, 2.63223e-02_jprb, &
     & 2.33672e-02_jprb, 2.07439e-02_jprb, 1.84150e-02_jprb, 1.63476e-02_jprb, 1.45123e-02_jprb/)
      selfrefo(:, 4) = (/ &
     & 5.76481e-02_jprb, 5.13686e-02_jprb, 4.57730e-02_jprb, 4.07870e-02_jprb, 3.63441e-02_jprb, &
     & 3.23851e-02_jprb, 2.88574e-02_jprb, 2.57140e-02_jprb, 2.29130e-02_jprb, 2.04171e-02_jprb/)
      selfrefo(:, 5) = (/ &
     & 6.92255e-02_jprb, 6.33521e-02_jprb, 5.79770e-02_jprb, 5.30580e-02_jprb, 4.85563e-02_jprb, &
     & 4.44365e-02_jprb, 4.06663e-02_jprb, 3.72160e-02_jprb, 3.40584e-02_jprb, 3.11687e-02_jprb/)
      selfrefo(:, 6) = (/ &
     & 6.07694e-02_jprb, 5.94182e-02_jprb, 5.80970e-02_jprb, 5.68052e-02_jprb, 5.55422e-02_jprb, &
     & 5.43072e-02_jprb, 5.30997e-02_jprb, 5.19190e-02_jprb, 5.07646e-02_jprb, 4.96358e-02_jprb/)
      selfrefo(:, 7) = (/ &
     & 6.23749e-02_jprb, 6.07744e-02_jprb, 5.92150e-02_jprb, 5.76956e-02_jprb, 5.62152e-02_jprb, &
     & 5.47728e-02_jprb, 5.33674e-02_jprb, 5.19980e-02_jprb, 5.06638e-02_jprb, 4.93638e-02_jprb/)
      selfrefo(:, 8) = (/ &
     & 6.90744e-02_jprb, 6.61811e-02_jprb, 6.34090e-02_jprb, 6.07530e-02_jprb, 5.82083e-02_jprb, &
     & 5.57702e-02_jprb, 5.34342e-02_jprb, 5.11960e-02_jprb, 4.90516e-02_jprb, 4.69970e-02_jprb/)
      selfrefo(:, 9) = (/ &
     & 8.08992e-02_jprb, 7.68876e-02_jprb, 7.30750e-02_jprb, 6.94514e-02_jprb, 6.60075e-02_jprb, &
     & 6.27344e-02_jprb, 5.96236e-02_jprb, 5.66670e-02_jprb, 5.38570e-02_jprb, 5.11864e-02_jprb/)
      selfrefo(:,10) = (/ &
     & 8.70197e-02_jprb, 8.27485e-02_jprb, 7.86870e-02_jprb, 7.48248e-02_jprb, 7.11522e-02_jprb, &
     & 6.76599e-02_jprb, 6.43389e-02_jprb, 6.11810e-02_jprb, 5.81781e-02_jprb, 5.53225e-02_jprb/)
      selfrefo(:,11) = (/ &
     & 8.84776e-02_jprb, 8.54262e-02_jprb, 8.24800e-02_jprb, 7.96354e-02_jprb, 7.68890e-02_jprb, &
     & 7.42373e-02_jprb, 7.16770e-02_jprb, 6.92050e-02_jprb, 6.68183e-02_jprb, 6.45139e-02_jprb/)
      selfrefo(:,12) = (/ &
     & 9.82552e-02_jprb, 9.25696e-02_jprb, 8.72130e-02_jprb, 8.21664e-02_jprb, 7.74118e-02_jprb, &
     & 7.29323e-02_jprb, 6.87121e-02_jprb, 6.47360e-02_jprb, 6.09900e-02_jprb, 5.74608e-02_jprb/)
      selfrefo(:,13) = (/ &
     & 9.32447e-02_jprb, 8.96818e-02_jprb, 8.62550e-02_jprb, 8.29592e-02_jprb, 7.97893e-02_jprb, &
     & 7.67405e-02_jprb, 7.38082e-02_jprb, 7.09880e-02_jprb, 6.82755e-02_jprb, 6.56667e-02_jprb/)
      selfrefo(:,14) = (/ &
     & 1.15363e-01_jprb, 1.08593e-01_jprb, 1.02220e-01_jprb, 9.62210e-02_jprb, 9.05741e-02_jprb, &
     & 8.52585e-02_jprb, 8.02549e-02_jprb, 7.55450e-02_jprb, 7.11115e-02_jprb, 6.69382e-02_jprb/)
      selfrefo(:,15) = (/ &
     & 1.23179e-01_jprb, 1.19247e-01_jprb, 1.15440e-01_jprb, 1.11755e-01_jprb, 1.08187e-01_jprb, &
     & 1.04734e-01_jprb, 1.01391e-01_jprb, 9.81540e-02_jprb, 9.50207e-02_jprb, 9.19875e-02_jprb/)
      selfrefo(:,16) = (/ &
     & 1.44104e-01_jprb, 1.36412e-01_jprb, 1.29130e-01_jprb, 1.22237e-01_jprb, 1.15712e-01_jprb, &
     & 1.09535e-01_jprb, 1.03688e-01_jprb, 9.81530e-02_jprb, 9.29135e-02_jprb, 8.79537e-02_jprb/)

if (lhook) call dr_hook('rrtm_kgb9',1,zhook_handle)
return

1001 continue
call abor1("rrtm_kgb9:error reading file radrrtm")

end subroutine rrtm_kgb9
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

