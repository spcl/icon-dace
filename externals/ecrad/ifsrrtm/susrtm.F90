! # 1 "ifsrrtm/susrtm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/susrtm.f90"
! this file has been modified for the use in icon

subroutine susrtm

!     adapted from e.j. mlawer, j. delamere, atmospheric & environmental research.
!     by jjmorcrette, ecmwf
!     modified to add arrays relevant to mapping for g-point reduction,
!     m.j. iacono, atmospheric & environmental research, inc. 
!     jjmorcrette 20010610 flexible configuration for number of g-points
!     ------------------------------------------------------------------

use parkind1  ,only : jprb ,   jpim
use ecradhook   ,only : lhook, dr_hook, jphook

use yoesrtm  , only : jpgpt, ngbsw, ngn
use yoesrtwn , only : ng      , nspa, nspb   , nmpsrtm, &
 & pref    , preflog , tref   , &
 & ngm     , wt      , ngc    , ngs
! & ngm     , wt      , ngc    , ngs , ngn    , ngbsw
! & wavenum1, wavenum2, delwave, pref, preflog, tref   , &

!     ------------------------------------------------------------------

implicit none

integer(kind=jpim) :: igc56(14), igc112(14) , igc224(14)
integer(kind=jpim) :: igs56(14), igs112(14) , igs224(14)

integer(kind=jpim) :: igm56(224),igm112(224), igm224(224)

integer(kind=jpim) :: ign56(56), ign112(112), ign224(224)
integer(kind=jpim) :: igb56(56), igb112(112), igb224(224)

real(kind=jphook) :: zhook_handle
!-----------------------------------------------------------------------
if (lhook) call dr_hook('susrtm',0,zhook_handle)

ng(:)     =(/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16 /)
nspa(:)   =(/  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1 /)
nspb(:)   =(/  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1 /)
nmpsrtm(:)=(/  6, 6, 5, 5, 5, 5, 5, 4, 4, 3, 2, 2, 1, 6 /)

!wavenum1( :) = (/&
! & 2600._jprb, 3250._jprb, 4000._jprb, 4650._jprb, 5150._jprb, 6150._jprb, 7700._jprb &
! & , 8050._jprb,12850._jprb,16000._jprb,22650._jprb,29000._jprb,38000._jprb,  820._jprb /)  
!wavenum2( :) = (/&
! & 3250._jprb, 4000._jprb, 4650._jprb, 5150._jprb, 6150._jprb, 7700._jprb, 8050._jprb &
! & ,12850._jprb,16000._jprb,22650._jprb,29000._jprb,38000._jprb,50000._jprb, 2600._jprb /)  
!delwave( :) = (/&
! & 650._jprb,  750._jprb,  650._jprb,  500._jprb, 1000._jprb, 1550._jprb,  350._jprb &
! & , 4800._jprb, 3150._jprb, 6650._jprb, 6350._jprb, 9000._jprb,12000._jprb, 1780._jprb /)  

!=====================================================================
! set arrays needed for the g-point reduction from 224 to 
! - either 112 for the high-resolution forecast model configuration
! - or 56 for the eps-type configuration  
! in the 14 sw bands:

! nb: this mapping from 224 to 112 points has been carefully selected to
! minimize the effect on the resulting fluxes and cooling rates, and
! caution should be used if the mapping is modified.
!     the further reduction to 56 for eps configuration is considered 
! acceptable, only because of the random perturbations introduced on 
! the total heating rates produced by the physical parametrization package.
! while a reduction to 56 obviously speeds up the model, it as obviously 
! reduces the accuracy that could be expected from the radiation scheme.

! jpgpt   the total number of new g-points (ngpt)
! ngc     the number of new g-points in each band (14)
! ngs     the cumulative sum of new g-points for each band (14)
! ngm     the index of each new g-point relative to the original
!         16 g-points for each band.
! ngn     the number of original g-points that are combined to make
!         each new g-point in each band.
! ngb     the band index for each new g-point.
! wt      rrtm weights for 16 g-points. (16)

!-- ecmwf eps model rrtm_sw configuration with 56 g-points
igc56(:) = (/ 3, 6, 4, 4, 5, 5, 1, 5, 4, 3, 3, 4, 3, 6 /)
igs56(:) = (/ 3, 9,13,17,22,27,28,33,37,40,43,47, 50, 56 /)

igm56(:) = (/ 1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,3, &            ! band 16
            & 1,1,2,2,3,3,3,4,4,4,5,5,5,6,6,6, &            ! band 17
            & 1,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4, &            ! band 18
            & 1,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4, &            ! band 19
            & 1,1,2,2,3,3,4,4,5,5,5,5,5,5,5,5, &            ! band 20
            & 1,1,2,2,3,3,4,4,5,5,5,5,5,5,5,5, &            ! band 21
            & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &            ! band 22
            & 1,1,1,1,2,2,3,3,4,4,5,5,5,5,5,5, &            ! band 23
            & 1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4, &            ! band 24
            & 1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3, &            ! band 25
            & 1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3, &            ! band 26
            & 1,1,2,2,3,3,4,4,4,4,4,4,4,4,4,4, &            ! band 27
            & 1,1,2,2,2,2,3,3,3,3,3,3,3,3,3,3, &            ! band 28
            & 1,1,2,2,3,3,3,3,4,4,4,4,5,5,6,6 /)            ! band 29

ign56(:) = (/ 4,4,8, &                                     ! band 16
            & 2,2,3,3,3,3, &                               ! band 17
            & 2,2,4,8, &                                   ! band 18
            & 2,2,4,8, &                                   ! band 19
            & 2,2,2,2,8, &                                 ! band 20
            & 2,2,2,2,8, &                                 ! band 21
            & 16, &                                        ! band 22
            & 4,2,2,2,6, &                                 ! band 23
            & 4,4,4,4, &                                   ! band 24
            & 2,4,10, &                                    ! band 25
            & 2,4,10, &                                    ! band 26
            & 2,2,2,10, &                                  ! band 27
            & 2,4,10, &                                    ! band 28
            & 2,2,4,4,2,2 /)                               ! band 29

igb56(:) = (/ 16,16,16, &                                  ! band 16
            & 17,17,17,17,17,17, &                         ! band 17
            & 18,18,18,18, &                               ! band 18
            & 19,19,19,19, &                               ! band 19
            & 20,20,20,20,20, &                            ! band 20
            & 21,21,21,21,21, &                            ! band 21
            & 22, &                                        ! band 22
            & 23,23,23,23,23, &                            ! band 23
            & 24,24,24,24, &                               ! band 24
            & 25,25,25, &                                  ! band 25
            & 26,26,26, &                                  ! band 26
            & 27,27,27,27, &                               ! band 27
            & 28,28,28, &                                  ! band 28
            & 29,29,29,29,29,29 /)                         ! band 29

!-------------------------------------------------------------------------------
!-- ecmwf high-resolution model rrtm_sw configuration with 112 g-points
! use this ngc, ngs, ngm, and ngn for reduced (112) g-point set
! (a related code change is required in modules parsrtm.f90 and yoesrtwn.f90)

igc112(:) = (/ 6,12, 8, 8,10,10, 2,10, 8, 6, 6, 8, 6,12 /)
igs112(:) = (/ 6,18,26,34,44,54,56,66,74,80,86,94,100,112 /)

!ngm(:)
igm112(:) = (/ 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 16
             & 1,2,3,4,5,6,6,7,8,8,9,10,10,11,12,12, &      ! band 17
             & 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 18
             & 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 19
             & 1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! band 20
             & 1,2,3,4,5,6,7,8,9,9,10,10,10,10,10,10, &     ! band 21
             & 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 22
             & 1,1,2,2,3,4,5,6,7,8,9,9,10,10,10,10, &       ! band 23
             & 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 24
             & 1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 25
             & 1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 26
             & 1,2,3,4,5,6,7,7,7,7,8,8,8,8,8,8, &           ! band 27
             & 1,2,3,3,4,4,5,5,5,5,6,6,6,6,6,6, &           ! band 28
             & 1,2,3,4,5,5,6,6,7,7,8,8,9,10,11,12 /)        ! band 29
!ngn(:)
ign112(:) = (/ 2,2,2,2,4,4, &                               ! band 16
             & 1,1,1,1,1,2,1,2,1,2,1,2, &                   ! band 17
             & 1,1,1,1,2,2,4,4, &                           ! band 18
             & 1,1,1,1,2,2,4,4, &                           ! band 19
             & 1,1,1,1,1,1,1,1,2,6, &                       ! band 20
             & 1,1,1,1,1,1,1,1,2,6, &                       ! band 21
             & 8,8, &                                       ! band 22
             & 2,2,1,1,1,1,1,1,2,4, &                       ! band 23
             & 2,2,2,2,2,2,2,2, &                           ! band 24
             & 1,1,2,2,4,6, &                               ! band 25
             & 1,1,2,2,4,6, &                               ! band 26
             & 1,1,1,1,1,1,4,6, &                           ! band 27
             & 1,1,2,2,4,6, &                               ! band 28
             & 1,1,1,1,2,2,2,2,1,1,1,1 /)                   ! band 29
!ngbsw(:)
igb112(:) = (/ 16,16,16,16,16,16, &                         ! band 16
             & 17,17,17,17,17,17,17,17,17,17,17,17, &       ! band 17
             & 18,18,18,18,18,18,18,18, &                   ! band 18
             & 19,19,19,19,19,19,19,19, &                   ! band 19
             & 20,20,20,20,20,20,20,20,20,20, &             ! band 20
             & 21,21,21,21,21,21,21,21,21,21, &             ! band 21
             & 22,22, &                                     ! band 22
             & 23,23,23,23,23,23,23,23,23,23, &             ! band 23
             & 24,24,24,24,24,24,24,24, &                   ! band 24
             & 25,25,25,25,25,25, &                         ! band 25
             & 26,26,26,26,26,26, &                         ! band 26
             & 27,27,27,27,27,27,27,27, &                   ! band 27
             & 28,28,28,28,28,28, &                         ! band 28
             & 29,29,29,29,29,29,29,29,29,29,29,29 /)       ! band 29

!-------------------------------------------------------------------------------
!-- original rrtm_sw configuration with 224 (14*16 g-points)
! use this ngc, ngs, ngm, and ngn for full (224) g-point set
! (a related code change is required in modules parsrtm.f90 and yoesrtwn.f90)
igc224(:) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16 /)
igs224(:) = (/ 16,32,48,64,80,96,112,128,144,160,176,192,208,224 /)

igm224(:) = (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 16
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 17
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 18
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 19
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 20
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 21
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 22
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 23
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 24
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 25
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 26
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 27
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 28
             & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 /)    ! band 29

ign224(:) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 16
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 17
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 18
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 19
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 20
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 21
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 22
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 23
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 24
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 25
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 26
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 27
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 28
             & 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /)           ! band 29

igb224(:) = (/ 16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16, &   ! band 16
             & 17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17, &   ! band 17
             & 18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18, &   ! band 18
             & 19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19, &   ! band 19
             & 20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20, &   ! band 20
             & 21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21, &   ! band 21
             & 22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22, &   ! band 22
             & 23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23, &   ! band 23
             & 24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24, &   ! band 24
             & 25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25, &   ! band 25
             & 26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26, &   ! band 26
             & 27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27, &   ! band 27
             & 28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28, &   ! band 28
             & 29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,28 /)   ! band 29

!=============================================================================

wt(:) =  (/ 0.1527534276_jprb, 0.1491729617_jprb, 0.1420961469_jprb, &
          & 0.1316886544_jprb, 0.1181945205_jprb, 0.1019300893_jprb, &
          & 0.0832767040_jprb, 0.0626720116_jprb, 0.0424925000_jprb, &
          & 0.0046269894_jprb, 0.0038279891_jprb, 0.0030260086_jprb, &
          & 0.0022199750_jprb, 0.0014140010_jprb, 0.0005330000_jprb, &
          & 0.0000750000_jprb /)

!=============================================================================

! these pressures are chosen such that the ln of the first pressure
! has only a few non-zero digits (i.e. ln(pref(1)) = 6.96000) and
!  each subsequent ln(pressure) differs from the previous one by 0.2.
pref = (/ &
 & 1.05363e+03_jprb,8.62642e+02_jprb,7.06272e+02_jprb,5.78246e+02_jprb,4.73428e+02_jprb, &
 & 3.87610e+02_jprb,3.17348e+02_jprb,2.59823e+02_jprb,2.12725e+02_jprb,1.74164e+02_jprb, &
 & 1.42594e+02_jprb,1.16746e+02_jprb,9.55835e+01_jprb,7.82571e+01_jprb,6.40715e+01_jprb, &
 & 5.24573e+01_jprb,4.29484e+01_jprb,3.51632e+01_jprb,2.87892e+01_jprb,2.35706e+01_jprb, &
 & 1.92980e+01_jprb,1.57998e+01_jprb,1.29358e+01_jprb,1.05910e+01_jprb,8.67114e+00_jprb, &
 & 7.09933e+00_jprb,5.81244e+00_jprb,4.75882e+00_jprb,3.89619e+00_jprb,3.18993e+00_jprb, &
 & 2.61170e+00_jprb,2.13828e+00_jprb,1.75067e+00_jprb,1.43333e+00_jprb,1.17351e+00_jprb, &
 & 9.60789e-01_jprb,7.86628e-01_jprb,6.44036e-01_jprb,5.27292e-01_jprb,4.31710e-01_jprb, &
 & 3.53455e-01_jprb,2.89384e-01_jprb,2.36928e-01_jprb,1.93980e-01_jprb,1.58817e-01_jprb, &
 & 1.30029e-01_jprb,1.06458e-01_jprb,8.71608e-02_jprb,7.13612e-02_jprb,5.84256e-02_jprb, &
 & 4.78349e-02_jprb,3.91639e-02_jprb,3.20647e-02_jprb,2.62523e-02_jprb,2.14936e-02_jprb, &
 & 1.75975e-02_jprb,1.44076e-02_jprb,1.17959e-02_jprb,9.65769e-03_jprb /)  
preflog = (/ &
 & 6.9600e+00_jprb, 6.7600e+00_jprb, 6.5600e+00_jprb, 6.3600e+00_jprb, 6.1600e+00_jprb, &
 & 5.9600e+00_jprb, 5.7600e+00_jprb, 5.5600e+00_jprb, 5.3600e+00_jprb, 5.1600e+00_jprb, &
 & 4.9600e+00_jprb, 4.7600e+00_jprb, 4.5600e+00_jprb, 4.3600e+00_jprb, 4.1600e+00_jprb, &
 & 3.9600e+00_jprb, 3.7600e+00_jprb, 3.5600e+00_jprb, 3.3600e+00_jprb, 3.1600e+00_jprb, &
 & 2.9600e+00_jprb, 2.7600e+00_jprb, 2.5600e+00_jprb, 2.3600e+00_jprb, 2.1600e+00_jprb, &
 & 1.9600e+00_jprb, 1.7600e+00_jprb, 1.5600e+00_jprb, 1.3600e+00_jprb, 1.1600e+00_jprb, &
 & 9.6000e-01_jprb, 7.6000e-01_jprb, 5.6000e-01_jprb, 3.6000e-01_jprb, 1.6000e-01_jprb, &
 & -4.0000e-02_jprb,-2.4000e-01_jprb,-4.4000e-01_jprb,-6.4000e-01_jprb,-8.4000e-01_jprb, &
 & -1.0400e+00_jprb,-1.2400e+00_jprb,-1.4400e+00_jprb,-1.6400e+00_jprb,-1.8400e+00_jprb, &
 & -2.0400e+00_jprb,-2.2400e+00_jprb,-2.4400e+00_jprb,-2.6400e+00_jprb,-2.8400e+00_jprb, &
 & -3.0400e+00_jprb,-3.2400e+00_jprb,-3.4400e+00_jprb,-3.6400e+00_jprb,-3.8400e+00_jprb, &
 & -4.0400e+00_jprb,-4.2400e+00_jprb,-4.4400e+00_jprb,-4.6400e+00_jprb /)  
! these are the temperatures associated with the respective 
! pressures for the mls standard atmosphere. 
tref = (/ &
 & 2.9420e+02_jprb, 2.8799e+02_jprb, 2.7894e+02_jprb, 2.6925e+02_jprb, 2.5983e+02_jprb, &
 & 2.5017e+02_jprb, 2.4077e+02_jprb, 2.3179e+02_jprb, 2.2306e+02_jprb, 2.1578e+02_jprb, &
 & 2.1570e+02_jprb, 2.1570e+02_jprb, 2.1570e+02_jprb, 2.1706e+02_jprb, 2.1858e+02_jprb, &
 & 2.2018e+02_jprb, 2.2174e+02_jprb, 2.2328e+02_jprb, 2.2479e+02_jprb, 2.2655e+02_jprb, &
 & 2.2834e+02_jprb, 2.3113e+02_jprb, 2.3401e+02_jprb, 2.3703e+02_jprb, 2.4022e+02_jprb, &
 & 2.4371e+02_jprb, 2.4726e+02_jprb, 2.5085e+02_jprb, 2.5457e+02_jprb, 2.5832e+02_jprb, &
 & 2.6216e+02_jprb, 2.6606e+02_jprb, 2.6999e+02_jprb, 2.7340e+02_jprb, 2.7536e+02_jprb, &
 & 2.7568e+02_jprb, 2.7372e+02_jprb, 2.7163e+02_jprb, 2.6955e+02_jprb, 2.6593e+02_jprb, &
 & 2.6211e+02_jprb, 2.5828e+02_jprb, 2.5360e+02_jprb, 2.4854e+02_jprb, 2.4348e+02_jprb, &
 & 2.3809e+02_jprb, 2.3206e+02_jprb, 2.2603e+02_jprb, 2.2000e+02_jprb, 2.1435e+02_jprb, &
 & 2.0887e+02_jprb, 2.0340e+02_jprb, 1.9792e+02_jprb, 1.9290e+02_jprb, 1.8809e+02_jprb, &
 & 1.8329e+02_jprb, 1.7849e+02_jprb, 1.7394e+02_jprb, 1.7212e+02_jprb /)  
!     -----------------------------------------------------------------

if (jpgpt == 56) then

!- 14
  ngc(:)=igc56(:)
  ngs(:)=igs56(:)
!- 14*16=224
  ngm(:)=igm56(:)

  ngn(1:56)=ign56(1:56)
  ngbsw(1:56)=igb56(1:56)

elseif (jpgpt == 112) then
!- 14
  ngc(:)=igc112(:)
  ngs(:)=igs112(:)
!- 14*16=224
  ngm(:)=igm112(:)

  ngn(1:112)=ign112(1:112)
  ngbsw(1:112)=igb112(1:112)

elseif (jpgpt == 224) then
!- 14
  ngc(:)=igc224(:)
  ngs(:)=igs224(:)
!- 14*16=224
  ngm(:)=igm224(:)

  ngn(1:224)=ign224(1:224)
  ngbsw(1:224)=igb224(1:224)

endif

!$acc update device(nspa, nspb, preflog, tref, ngc)

!     -----------------------------------------------------------------
if (lhook) call dr_hook('susrtm',1,zhook_handle)
end subroutine susrtm

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

