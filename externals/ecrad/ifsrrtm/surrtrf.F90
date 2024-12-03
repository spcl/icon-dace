! # 1 "ifsrrtm/surrtrf.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/surrtrf.f90"
! this file has been modified for the use in icon

subroutine surrtrf

!     adapted from eli j. mlawer, atmospheric & environmental research.
!     by jjmorcrette, ecmwf
!     ------------------------------------------------------------------

use parkind1  ,only : jprb
use ecradhook   ,only : lhook, dr_hook, jphook

use yoerrtrf , only : pref      ,preflog   ,tref,chi_mls

!     ------------------------------------------------------------------

!     these pressures are chosen such that the ln of the first pressure
!     has only a few non-zero digits (i.e. ln(pref(1)) = 6.96000) and
!     each subsequent ln(pressure) differs from the previous one by 0.2.

implicit none
real(kind=jphook) :: zhook_handle
if (lhook) call dr_hook('surrtrf',0,zhook_handle)
pref( :) = (/&
 & 1.05363e+03_jprb,8.62642e+02_jprb,7.06272e+02_jprb,5.78246e+02_jprb,4.73428e+02_jprb,&
 & 3.87610e+02_jprb,3.17348e+02_jprb,2.59823e+02_jprb,2.12725e+02_jprb,1.74164e+02_jprb,&
 & 1.42594e+02_jprb,1.16746e+02_jprb,9.55835e+01_jprb,7.82571e+01_jprb,6.40715e+01_jprb,&
 & 5.24573e+01_jprb,4.29484e+01_jprb,3.51632e+01_jprb,2.87892e+01_jprb,2.35706e+01_jprb,&
 & 1.92980e+01_jprb,1.57998e+01_jprb,1.29358e+01_jprb,1.05910e+01_jprb,8.67114e+00_jprb,&
 & 7.09933e+00_jprb,5.81244e+00_jprb,4.75882e+00_jprb,3.89619e+00_jprb,3.18993e+00_jprb,&
 & 2.61170e+00_jprb,2.13828e+00_jprb,1.75067e+00_jprb,1.43333e+00_jprb,1.17351e+00_jprb,&
 & 9.60789e-01_jprb,7.86628e-01_jprb,6.44036e-01_jprb,5.27292e-01_jprb,4.31710e-01_jprb,&
 & 3.53455e-01_jprb,2.89384e-01_jprb,2.36928e-01_jprb,1.93980e-01_jprb,1.58817e-01_jprb,&
 & 1.30029e-01_jprb,1.06458e-01_jprb,8.71608e-02_jprb,7.13612e-02_jprb,5.84256e-02_jprb,&
 & 4.78349e-02_jprb,3.91639e-02_jprb,3.20647e-02_jprb,2.62523e-02_jprb,2.14936e-02_jprb,&
 & 1.75975e-02_jprb,1.44076e-02_jprb,1.17959e-02_jprb,9.65769e-03_jprb/)  

preflog( :) = (/&
 & 6.9600e+00_jprb, 6.7600e+00_jprb, 6.5600e+00_jprb, 6.3600e+00_jprb, 6.1600e+00_jprb,&
 & 5.9600e+00_jprb, 5.7600e+00_jprb, 5.5600e+00_jprb, 5.3600e+00_jprb, 5.1600e+00_jprb,&
 & 4.9600e+00_jprb, 4.7600e+00_jprb, 4.5600e+00_jprb, 4.3600e+00_jprb, 4.1600e+00_jprb,&
 & 3.9600e+00_jprb, 3.7600e+00_jprb, 3.5600e+00_jprb, 3.3600e+00_jprb, 3.1600e+00_jprb,&
 & 2.9600e+00_jprb, 2.7600e+00_jprb, 2.5600e+00_jprb, 2.3600e+00_jprb, 2.1600e+00_jprb,&
 & 1.9600e+00_jprb, 1.7600e+00_jprb, 1.5600e+00_jprb, 1.3600e+00_jprb, 1.1600e+00_jprb,&
 & 9.6000e-01_jprb, 7.6000e-01_jprb, 5.6000e-01_jprb, 3.6000e-01_jprb, 1.6000e-01_jprb,&
 & -4.0000e-02_jprb,-2.4000e-01_jprb,-4.4000e-01_jprb,-6.4000e-01_jprb,-8.4000e-01_jprb,&
 & -1.0400e+00_jprb,-1.2400e+00_jprb,-1.4400e+00_jprb,-1.6400e+00_jprb,-1.8400e+00_jprb,&
 & -2.0400e+00_jprb,-2.2400e+00_jprb,-2.4400e+00_jprb,-2.6400e+00_jprb,-2.8400e+00_jprb,&
 & -3.0400e+00_jprb,-3.2400e+00_jprb,-3.4400e+00_jprb,-3.6400e+00_jprb,-3.8400e+00_jprb,&
 & -4.0400e+00_jprb,-4.2400e+00_jprb,-4.4400e+00_jprb,-4.6400e+00_jprb/)  

!     these are the temperatures associated with the respective 
!     pressures for the mls standard atmosphere. 
tref( :) = (/&
 & 2.9420e+02_jprb, 2.8799e+02_jprb, 2.7894e+02_jprb, 2.6925e+02_jprb, 2.5983e+02_jprb,&
 & 2.5017e+02_jprb, 2.4077e+02_jprb, 2.3179e+02_jprb, 2.2306e+02_jprb, 2.1578e+02_jprb,&
 & 2.1570e+02_jprb, 2.1570e+02_jprb, 2.1570e+02_jprb, 2.1706e+02_jprb, 2.1858e+02_jprb,&
 & 2.2018e+02_jprb, 2.2174e+02_jprb, 2.2328e+02_jprb, 2.2479e+02_jprb, 2.2655e+02_jprb,&
 & 2.2834e+02_jprb, 2.3113e+02_jprb, 2.3401e+02_jprb, 2.3703e+02_jprb, 2.4022e+02_jprb,&
 & 2.4371e+02_jprb, 2.4726e+02_jprb, 2.5085e+02_jprb, 2.5457e+02_jprb, 2.5832e+02_jprb,&
 & 2.6216e+02_jprb, 2.6606e+02_jprb, 2.6999e+02_jprb, 2.7340e+02_jprb, 2.7536e+02_jprb,&
 & 2.7568e+02_jprb, 2.7372e+02_jprb, 2.7163e+02_jprb, 2.6955e+02_jprb, 2.6593e+02_jprb,&
 & 2.6211e+02_jprb, 2.5828e+02_jprb, 2.5360e+02_jprb, 2.4854e+02_jprb, 2.4348e+02_jprb,&
 & 2.3809e+02_jprb, 2.3206e+02_jprb, 2.2603e+02_jprb, 2.2000e+02_jprb, 2.1435e+02_jprb,&
 & 2.0887e+02_jprb, 2.0340e+02_jprb, 1.9792e+02_jprb, 1.9290e+02_jprb, 1.8809e+02_jprb,&
 & 1.8329e+02_jprb, 1.7849e+02_jprb, 1.7394e+02_jprb, 1.7212e+02_jprb/)  

  chi_mls(1,1:12) = (/ &
     &  1.8760e-02_jprb, 1.2223e-02_jprb, 5.8909e-03_jprb, 2.7675e-03_jprb, 1.4065e-03_jprb, &
     &  7.5970e-04_jprb, 3.8876e-04_jprb, 1.6542e-04_jprb, 3.7190e-05_jprb, 7.4765e-06_jprb, &
     &  4.3082e-06_jprb, 3.3319e-06_jprb/)
  chi_mls(1,13:59) = (/ &
     &  3.2039e-06_jprb,  3.1619e-06_jprb,  3.2524e-06_jprb,  3.4226e-06_jprb,  3.6288e-06_jprb, &
     &  3.9148e-06_jprb,  4.1488e-06_jprb,  4.3081e-06_jprb,  4.4420e-06_jprb,  4.5778e-06_jprb, &
     &  4.7087e-06_jprb,  4.7943e-06_jprb,  4.8697e-06_jprb,  4.9260e-06_jprb,  4.9669e-06_jprb, &
     &  4.9963e-06_jprb,  5.0527e-06_jprb,  5.1266e-06_jprb,  5.2503e-06_jprb,  5.3571e-06_jprb, &
     &  5.4509e-06_jprb,  5.4830e-06_jprb,  5.5000e-06_jprb,  5.5000e-06_jprb,  5.4536e-06_jprb, &
     &  5.4047e-06_jprb,  5.3558e-06_jprb,  5.2533e-06_jprb,  5.1436e-06_jprb,  5.0340e-06_jprb, &
     &  4.8766e-06_jprb,  4.6979e-06_jprb,  4.5191e-06_jprb,  4.3360e-06_jprb,  4.1442e-06_jprb, &
     &  3.9523e-06_jprb,  3.7605e-06_jprb,  3.5722e-06_jprb,  3.3855e-06_jprb,  3.1988e-06_jprb, &
     &  3.0121e-06_jprb,  2.8262e-06_jprb,  2.6407e-06_jprb,  2.4552e-06_jprb,  2.2696e-06_jprb, &
     &  4.3360e-06_jprb,  4.1442e-06_jprb/)
  chi_mls(2,1:12) = (/ &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb/)
  chi_mls(2,13:59) = (/ &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb,  3.5500e-04_jprb, &
     &  3.5500e-04_jprb,  3.5471e-04_jprb,  3.5427e-04_jprb,  3.5384e-04_jprb,  3.5340e-04_jprb, &
     &  3.5500e-04_jprb,  3.5500e-04_jprb/)
  chi_mls(3,1:12) = (/ &
     &  3.0170e-08_jprb,  3.4725e-08_jprb,  4.2477e-08_jprb,  5.2759e-08_jprb,  6.6944e-08_jprb, &
     &  8.7130e-08_jprb,  1.1391e-07_jprb,  1.5677e-07_jprb,  2.1788e-07_jprb,  3.2443e-07_jprb, &
     &  4.6594e-07_jprb,  5.6806e-07_jprb/)
  chi_mls(3,13:59) = (/ &
     &  6.9607e-07_jprb,  1.1186e-06_jprb,  1.7618e-06_jprb,  2.3269e-06_jprb,  2.9577e-06_jprb, &
     &  3.6593e-06_jprb,  4.5950e-06_jprb,  5.3189e-06_jprb,  5.9618e-06_jprb,  6.5113e-06_jprb, &
     &  7.0635e-06_jprb,  7.6917e-06_jprb,  8.2577e-06_jprb,  8.7082e-06_jprb,  8.8325e-06_jprb, &
     &  8.7149e-06_jprb,  8.0943e-06_jprb,  7.3307e-06_jprb,  6.3101e-06_jprb,  5.3672e-06_jprb, &
     &  4.4829e-06_jprb,  3.8391e-06_jprb,  3.2827e-06_jprb,  2.8235e-06_jprb,  2.4906e-06_jprb, &
     &  2.1645e-06_jprb,  1.8385e-06_jprb,  1.6618e-06_jprb,  1.5052e-06_jprb,  1.3485e-06_jprb, &
     &  1.1972e-06_jprb,  1.0482e-06_jprb,  8.9926e-07_jprb,  7.6343e-07_jprb,  6.5381e-07_jprb, &
     &  5.4419e-07_jprb,  4.3456e-07_jprb,  3.6421e-07_jprb,  3.1194e-07_jprb,  2.5967e-07_jprb, &
     &  2.0740e-07_jprb,  1.9146e-07_jprb,  1.9364e-07_jprb,  1.9582e-07_jprb,  1.9800e-07_jprb, &
     &  7.6343e-07_jprb,  6.5381e-07_jprb/)
  chi_mls(4,1:12) = (/ &
     &  3.2000e-07_jprb,  3.2000e-07_jprb,  3.2000e-07_jprb,  3.2000e-07_jprb,  3.2000e-07_jprb, &
     &  3.1965e-07_jprb,  3.1532e-07_jprb,  3.0383e-07_jprb,  2.9422e-07_jprb,  2.8495e-07_jprb, &
     &  2.7671e-07_jprb,  2.6471e-07_jprb/)
  chi_mls(4,13:59) = (/ &
     &  2.4285e-07_jprb,  2.0955e-07_jprb,  1.7195e-07_jprb,  1.3749e-07_jprb,  1.1332e-07_jprb, &
     &  1.0035e-07_jprb,  9.1281e-08_jprb,  8.5463e-08_jprb,  8.0363e-08_jprb,  7.3372e-08_jprb, &
     &  6.5975e-08_jprb,  5.6039e-08_jprb,  4.7090e-08_jprb,  3.9977e-08_jprb,  3.2979e-08_jprb, &
     &  2.6064e-08_jprb,  2.1066e-08_jprb,  1.6592e-08_jprb,  1.3017e-08_jprb,  1.0090e-08_jprb, &
     &  7.6249e-09_jprb,  6.1159e-09_jprb,  4.6672e-09_jprb,  3.2857e-09_jprb,  2.8484e-09_jprb, &
     &  2.4620e-09_jprb,  2.0756e-09_jprb,  1.8551e-09_jprb,  1.6568e-09_jprb,  1.4584e-09_jprb, &
     &  1.3195e-09_jprb,  1.2072e-09_jprb,  1.0948e-09_jprb,  9.9780e-10_jprb,  9.3126e-10_jprb, &
     &  8.6472e-10_jprb,  7.9818e-10_jprb,  7.5138e-10_jprb,  7.1367e-10_jprb,  6.7596e-10_jprb, &
     &  6.3825e-10_jprb,  6.0981e-10_jprb,  5.8600e-10_jprb,  5.6218e-10_jprb,  5.3837e-10_jprb, &
     &  9.9780e-10_jprb,  9.3126e-10_jprb/)
  chi_mls(5,1:12) = (/ &
     &  1.5000e-07_jprb,  1.4306e-07_jprb,  1.3474e-07_jprb,  1.3061e-07_jprb,  1.2793e-07_jprb, &
     &  1.2038e-07_jprb,  1.0798e-07_jprb,  9.4238e-08_jprb,  7.9488e-08_jprb,  6.1386e-08_jprb, &
     &  4.5563e-08_jprb,  3.3475e-08_jprb/)
  chi_mls(5,13:59) = (/ &
     &  2.5118e-08_jprb,  1.8671e-08_jprb,  1.4349e-08_jprb,  1.2501e-08_jprb,  1.2407e-08_jprb, &
     &  1.3472e-08_jprb,  1.4900e-08_jprb,  1.6079e-08_jprb,  1.7156e-08_jprb,  1.8616e-08_jprb, &
     &  2.0106e-08_jprb,  2.1654e-08_jprb,  2.3096e-08_jprb,  2.4340e-08_jprb,  2.5643e-08_jprb, &
     &  2.6990e-08_jprb,  2.8456e-08_jprb,  2.9854e-08_jprb,  3.0943e-08_jprb,  3.2023e-08_jprb, &
     &  3.3101e-08_jprb,  3.4260e-08_jprb,  3.5360e-08_jprb,  3.6397e-08_jprb,  3.7310e-08_jprb, &
     &  3.8217e-08_jprb,  3.9123e-08_jprb,  4.1303e-08_jprb,  4.3652e-08_jprb,  4.6002e-08_jprb, &
     &  5.0289e-08_jprb,  5.5446e-08_jprb,  6.0603e-08_jprb,  6.8946e-08_jprb,  8.3652e-08_jprb, &
     &  9.8357e-08_jprb,  1.1306e-07_jprb,  1.4766e-07_jprb,  1.9142e-07_jprb,  2.3518e-07_jprb, &
     &  2.7894e-07_jprb,  3.5001e-07_jprb,  4.3469e-07_jprb,  5.1938e-07_jprb,  6.0407e-07_jprb, &
     &  6.8946e-08_jprb,  8.3652e-08_jprb/)
  chi_mls(6,1:12) = (/ &
     &  1.7000e-06_jprb,  1.7000e-06_jprb,  1.6999e-06_jprb,  1.6904e-06_jprb,  1.6671e-06_jprb, &
     &  1.6351e-06_jprb,  1.6098e-06_jprb,  1.5590e-06_jprb,  1.5120e-06_jprb,  1.4741e-06_jprb, &
     &  1.4385e-06_jprb,  1.4002e-06_jprb/)
  chi_mls(6,13:59) = (/ &
     &  1.3573e-06_jprb,  1.3130e-06_jprb,  1.2512e-06_jprb,  1.1668e-06_jprb,  1.0553e-06_jprb, &
     &  9.3281e-07_jprb,  8.1217e-07_jprb,  7.5239e-07_jprb,  7.0728e-07_jprb,  6.6722e-07_jprb, &
     &  6.2733e-07_jprb,  5.8604e-07_jprb,  5.4769e-07_jprb,  5.1480e-07_jprb,  4.8206e-07_jprb, &
     &  4.4943e-07_jprb,  4.1702e-07_jprb,  3.8460e-07_jprb,  3.5200e-07_jprb,  3.1926e-07_jprb, &
     &  2.8646e-07_jprb,  2.5498e-07_jprb,  2.2474e-07_jprb,  1.9588e-07_jprb,  1.8295e-07_jprb, &
     &  1.7089e-07_jprb,  1.5882e-07_jprb,  1.5536e-07_jprb,  1.5304e-07_jprb,  1.5072e-07_jprb, &
     &  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb, &
     &  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb, &
     &  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb,  1.5000e-07_jprb, &
     &  1.5000e-07_jprb,  1.5000e-07_jprb/)
  chi_mls(7,1:12) = (/ &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb/)
  chi_mls(7,13:59) = (/ &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb,  0.2090_jprb, &
     &  0.2090_jprb,  0.2090_jprb/)

  !$acc update device(pref, preflog, tref, chi_mls)

!     -----------------------------------------------------------------
if (lhook) call dr_hook('surrtrf',1,zhook_handle)
end subroutine surrtrf
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

