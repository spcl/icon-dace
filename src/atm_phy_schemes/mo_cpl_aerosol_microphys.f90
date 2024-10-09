! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! This module contains routines needed for basic coupling between cloud microphysics and the aerosol climatology

MODULE mo_cpl_aerosol_microphys


  USE mo_kind,               ONLY: ireals=>wp, iintegers=>i4
  USE mo_exception,          ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: lookupcreate_segalkhain, specccn_segalkhain, specccn_segalkhain_simple, &
             ncn_from_tau_aerosol_speccnconst, ncn_from_tau_aerosol_speccnconst_dust, &
             ice_nucleation

  INTERFACE ice_nucleation
    MODULE PROCEDURE ice_nucleation_demott
    MODULE PROCEDURE ice_nucleation_cooper
  END INTERFACE ice_nucleation


! Type declaration for a general 2D equidistant lookup table:
TYPE lookupt_2D
  INTEGER(KIND=iintegers) :: n1  ! number of grid points in x1-direction
  INTEGER(KIND=iintegers) :: n2  ! number of grid points in x2-direction
  REAL(KIND=ireals), DIMENSION(:), ALLOCATABLE :: x1  ! grid vector in x1-direction
  REAL(KIND=ireals), DIMENSION(:), ALLOCATABLE :: x2  ! grid vector in x1-direction
  REAL(KIND=ireals)                     :: dx1        ! dx1   (grid distance w.r.t. x1)
  REAL(KIND=ireals)                     :: dx2        ! dx2   (grid distance w.r.t. x2)
  REAL(KIND=ireals)                     :: odx1       ! one over dx 1
  REAL(KIND=ireals)                     :: odx2       ! one over dx 2
  REAL(KIND=ireals), DIMENSION(:,:), ALLOCATABLE :: ltable
END TYPE lookupt_2D

!------------------------------------------------------------------------------

! Type declaration for a general 4D equidistant lookup table:
TYPE lookupt_4D
  INTEGER(KIND=iintegers) :: n1  ! number of grid points in x1-direction
  INTEGER(KIND=iintegers) :: n2  ! number of grid points in x2-direction
  INTEGER(KIND=iintegers) :: n3  ! number of grid points in x3-direction
  INTEGER(KIND=iintegers) :: n4  ! number of grid points in x4-direction
  REAL(KIND=ireals), DIMENSION(:), ALLOCATABLE :: x1  ! grid vector in x1-direction
  REAL(KIND=ireals), DIMENSION(:), ALLOCATABLE :: x2  ! grid vector in x1-direction
  REAL(KIND=ireals), DIMENSION(:), ALLOCATABLE :: x3  ! grid vector in x1-direction
  REAL(KIND=ireals), DIMENSION(:), ALLOCATABLE :: x4  ! grid vector in x1-direction
  REAL(KIND=ireals)                     :: dx1        ! dx1   (grid distance w.r.t. x1)
  REAL(KIND=ireals)                     :: dx2        ! dx2   (grid distance w.r.t. x2)
  REAL(KIND=ireals)                     :: dx3        ! dx3   (grid distance w.r.t. x3)
  REAL(KIND=ireals)                     :: dx4        ! dx4   (grid distance w.r.t. x4)
  REAL(KIND=ireals)                     :: odx1       ! one over dx 1
  REAL(KIND=ireals)                     :: odx2       ! one over dx 2
  REAL(KIND=ireals)                     :: odx3       ! one over dx 3
  REAL(KIND=ireals)                     :: odx4       ! one over dx 4
  REAL(KIND=ireals), DIMENSION(:,:,:,:), ALLOCATABLE :: ltable
END TYPE lookupt_4D

TYPE(lookupt_4D) :: ltab4D
TYPE(lookupt_2D) :: ltab2D

! Parameters for Segal-Khain parameterization (aerosol-microphysics coupling)
! ---------------------------------------------------------------------------

REAL(KIND=ireals), PARAMETER ::       &
  r2_fix    = 0.03_ireals,            & ! Parameters for simplified lookup table computation; 
  lsigs_fix = 0.3_ireals                ! relevant for r2_lsigs_are_fixed = .TRUE.

LOGICAL,       PARAMETER ::       &
  r2_lsigs_are_fixed = .TRUE. ,   & !
  lincloud           = .FALSE.      ! ignore in-cloud nucleation

CONTAINS


SUBROUTINE nccn_lookupcreate_segalkhain_4D ( ltab )

  IMPLICIT NONE

  TYPE(lookupt_4D), INTENT(out) :: ltab

  ltab%n1 = 3    ! for R2
  ltab%n2 = 5    ! for log(sigma_s)
  ltab%n3 = 8    ! for N_CN
  ltab%n4 = 11   ! for w_cloudbase
  
  ALLOCATE( ltab%x1(ltab%n1) )
  ALLOCATE( ltab%x2(ltab%n2) )
  ALLOCATE( ltab%x3(ltab%n3) )
  ALLOCATE( ltab%x4(ltab%n4) )
  ALLOCATE( ltab%ltable(ltab%n1,ltab%n2,ltab%n3,ltab%n4) )

  ltab%x1 = (/  0.02_ireals, 0.03_ireals, 0.04_ireals /)
  ltab%x2 = (/  0.1_ireals, 0.2_ireals, 0.3_ireals, 0.4_ireals, 0.5_ireals /)
  ltab%x3 = LOG( (/  5.000E+07_ireals, 1.000E+08_ireals, 2.000E+08_ireals, 4.000E+08_ireals, 8.000E+08_ireals, &
                     1.600E+09_ireals, 3.200E+09_ireals, 6.400E+09_ireals /) )
  ltab%x4 = (/  0.0_ireals, 0.5_ireals, 1.0_ireals, 1.5_ireals, 2.0_ireals, &
                2.5_ireals, 3.0_ireals, 3.5_ireals, 4.0_ireals, 4.5_ireals, 5.0_ireals /)

  ltab%dx1 = ltab%x1(2) - ltab%x1(1)
  ltab%dx2 = ltab%x2(2) - ltab%x2(1)
  ltab%dx3 = ltab%x3(2) - ltab%x3(1)
  ltab%dx4 = ltab%x4(2) - ltab%x4(1)

  ltab%odx1 = 1.0_ireals / ltab%dx1
  ltab%odx2 = 1.0_ireals / ltab%dx2
  ltab%odx3 = 1.0_ireals / ltab%dx3
  ltab%odx4 = 1.0_ireals / ltab%dx4

  ! Ncn                            50              100                  200                 400        
  !                               800             1600                 3200                6400
  ! interpolated (R2=0.02mum, wcb=0.0m/s)
  ltab%ltable(1,1,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(1,2,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(1,3,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(1,4,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(1,5,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)

  ! table4a (R2=0.02mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 extrapolated)
  ltab%ltable(1,1,:, 2) = (/ 4.2200E+07_ireals, 7.0200E+07_ireals, 1.1220E+08_ireals, 1.7310E+08_ireals,  &
                             2.6370E+08_ireals, 3.9750E+08_ireals, 5.9919E+08_ireals, 9.0322E+08_ireals /)
  ltab%ltable(1,2,:, 2) = (/ 3.5500E+07_ireals, 6.0100E+07_ireals, 1.0000E+08_ireals, 1.6390E+08_ireals,  &
                             2.6450E+08_ireals, 4.1840E+08_ireals, 6.6185E+08_ireals, 1.0469E+09_ireals /)
  ltab%ltable(1,3,:, 2) = (/ 3.2600E+07_ireals, 5.6300E+07_ireals, 9.6700E+07_ireals, 1.6390E+08_ireals,  &
                             2.7200E+08_ireals, 4.3850E+08_ireals, 7.0692E+08_ireals, 1.1396E+09_ireals /)
  ltab%ltable(1,4,:, 2) = (/ 3.0900E+07_ireals, 5.4400E+07_ireals, 9.4600E+07_ireals, 1.6240E+08_ireals,  &
                             2.7190E+08_ireals, 4.3350E+08_ireals, 6.9114E+08_ireals, 1.1019E+09_ireals /)
  ltab%ltable(1,5,:, 2) = (/ 2.9400E+07_ireals, 5.1900E+07_ireals, 8.9900E+07_ireals, 1.5060E+08_ireals,  &
                             2.3650E+08_ireals, 3.6440E+08_ireals, 5.6147E+08_ireals, 8.6511E+08_ireals /)

  ! table4b (R2=0.02mum, wcb=1.0m/s) (for Ncn=50 interpolted and Ncn=6400 extrapolated)
  ltab%ltable(1,1,:, 3) = (/ 5.0000E+07_ireals, 9.1500E+07_ireals, 1.5870E+08_ireals, 2.6440E+08_ireals,  &
                             4.2310E+08_ireals, 6.7250E+08_ireals, 1.0182E+09_ireals, 1.5416E+09_ireals /)
  ltab%ltable(1,2,:, 3) = (/ 4.4695E+07_ireals, 7.7100E+07_ireals, 1.3300E+08_ireals, 2.2490E+08_ireals,  &
                             3.7650E+08_ireals, 6.1570E+08_ireals, 9.7560E+08_ireals, 1.5459E+09_ireals /)
  ltab%ltable(1,3,:, 3) = (/ 4.0000E+07_ireals, 7.0000E+07_ireals, 1.2250E+08_ireals, 2.1200E+08_ireals,  &
                             3.6210E+08_ireals, 6.0530E+08_ireals, 9.6480E+08_ireals, 1.5378E+09_ireals /)
  ltab%ltable(1,4,:, 3) = (/ 3.7196E+07_ireals, 6.5800E+07_ireals, 1.1640E+08_ireals, 2.0400E+08_ireals,  &
                             3.5060E+08_ireals, 5.8440E+08_ireals, 9.1430E+08_ireals, 1.4304E+09_ireals /)
  ltab%ltable(1,5,:, 3) = (/ 3.5252E+07_ireals, 6.2300E+07_ireals, 1.1010E+08_ireals, 1.9130E+08_ireals,  &
                             3.2060E+08_ireals, 5.0130E+08_ireals, 7.6040E+08_ireals, 1.1534E+09_ireals /)

  ! interpolated (R2=0.02mum, wcb=1.5m/s)
  ltab%ltable(1,1,:, 4) = (/ 5.0000E+07_ireals, 9.4333E+07_ireals, 1.7247E+08_ireals, 3.0063E+08_ireals,  &
                             5.0363E+08_ireals, 8.2593E+08_ireals, 1.3044E+09_ireals, 2.0190E+09_ireals /)
  ltab%ltable(1,2,:, 4) = (/ 4.6463E+07_ireals, 8.4358E+07_ireals, 1.4740E+08_ireals, 2.5460E+08_ireals,  &
                             4.3330E+08_ireals, 7.2427E+08_ireals, 1.1768E+09_ireals, 1.8780E+09_ireals /)
  ltab%ltable(1,3,:, 4) = (/ 4.3318E+07_ireals, 7.6365E+07_ireals, 1.3463E+08_ireals, 2.3580E+08_ireals,  &
                             4.0770E+08_ireals, 6.9217E+08_ireals, 1.1307E+09_ireals, 1.8105E+09_ireals /)
  ltab%ltable(1,4,:, 4) = (/ 4.0232E+07_ireals, 7.1461E+07_ireals, 1.2693E+08_ireals, 2.2420E+08_ireals,  &
                             3.8983E+08_ireals, 6.6070E+08_ireals, 1.0666E+09_ireals, 1.6660E+09_ireals /)
  ltab%ltable(1,5,:, 4) = (/ 3.8124E+07_ireals, 6.7543E+07_ireals, 1.1967E+08_ireals, 2.0983E+08_ireals,  &
                             3.5803E+08_ireals, 5.8013E+08_ireals, 8.9917E+08_ireals, 1.3466E+09_ireals /)

  ! interpolated (R2=0.02mum, wcb=2.0m/s)
  ltab%ltable(1,1,:, 5) = (/ 5.0000E+07_ireals, 9.7167E+07_ireals, 1.8623E+08_ireals, 3.3687E+08_ireals,  &
                             5.8417E+08_ireals, 9.7937E+08_ireals, 1.5906E+09_ireals, 2.4963E+09_ireals /)
  ltab%ltable(1,2,:, 5) = (/ 4.8232E+07_ireals, 9.1616E+07_ireals, 1.6180E+08_ireals, 2.8430E+08_ireals,  &
                             4.9010E+08_ireals, 8.3283E+08_ireals, 1.3780E+09_ireals, 2.2101E+09_ireals /)
  ltab%ltable(1,3,:, 5) = (/ 4.6636E+07_ireals, 8.2729E+07_ireals, 1.4677E+08_ireals, 2.5960E+08_ireals,  &
                             4.5330E+08_ireals, 7.7903E+08_ireals, 1.2967E+09_ireals, 2.0831E+09_ireals /)
  ltab%ltable(1,4,:, 5) = (/ 4.3267E+07_ireals, 7.7121E+07_ireals, 1.3747E+08_ireals, 2.4440E+08_ireals,  &
                             4.2907E+08_ireals, 7.3700E+08_ireals, 1.2190E+09_ireals, 1.9016E+09_ireals /)
  ltab%ltable(1,5,:, 5) = (/ 4.0995E+07_ireals, 7.2786E+07_ireals, 1.2923E+08_ireals, 2.2837E+08_ireals,  &
                             3.9547E+08_ireals, 6.5897E+08_ireals, 1.0379E+09_ireals, 1.5398E+09_ireals /)

  ! table4c (R2=0.02mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 interpolated)
  ltab%ltable(1,1,:, 6) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.7310E+08_ireals,  &
                             6.6470E+08_ireals, 1.1328E+09_ireals, 1.8768E+09_ireals, 2.9737E+09_ireals /)
  ltab%ltable(1,2,:, 6) = (/ 5.0000E+07_ireals, 9.8874E+07_ireals, 1.7620E+08_ireals, 3.1400E+08_ireals,  &
                             5.4690E+08_ireals, 9.4140E+08_ireals, 1.5792E+09_ireals, 2.5422E+09_ireals /)
  ltab%ltable(1,3,:, 6) = (/ 4.9954E+07_ireals, 8.9094E+07_ireals, 1.5890E+08_ireals, 2.8340E+08_ireals,  &
                             4.9890E+08_ireals, 8.6590E+08_ireals, 1.4626E+09_ireals, 2.3558E+09_ireals /)
  ltab%ltable(1,4,:, 6) = (/ 4.6303E+07_ireals, 8.2782E+07_ireals, 1.4800E+08_ireals, 2.6460E+08_ireals,  &
                             4.6830E+08_ireals, 8.1330E+08_ireals, 1.3713E+09_ireals, 2.1372E+09_ireals /)
  ltab%ltable(1,5,:, 6) = (/ 4.3866E+07_ireals, 7.8029E+07_ireals, 1.3880E+08_ireals, 2.4690E+08_ireals,  &
                             4.3290E+08_ireals, 7.3780E+08_ireals, 1.1767E+09_ireals, 1.7330E+09_ireals /)

  ! interpolated (R2=0.02mum, wcb=3.0m/s)
  ltab%ltable(1,1,:, 7) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.7848E+08_ireals,  &
                             6.8938E+08_ireals, 1.1969E+09_ireals, 2.0185E+09_ireals, 3.2555E+09_ireals /)
  ltab%ltable(1,2,:, 7) = (/ 5.0000E+07_ireals, 9.9099E+07_ireals, 1.8096E+08_ireals, 3.2572E+08_ireals,  &
                             5.6896E+08_ireals, 9.9368E+08_ireals, 1.6830E+09_ireals, 2.7451E+09_ireals /)
  ltab%ltable(1,3,:, 7) = (/ 4.9963E+07_ireals, 9.1275E+07_ireals, 1.6469E+08_ireals, 2.9424E+08_ireals,  &
                             5.2046E+08_ireals, 9.0842E+08_ireals, 1.5479E+09_ireals, 2.5260E+09_ireals /)
  ltab%ltable(1,4,:, 7) = (/ 4.7042E+07_ireals, 8.5573E+07_ireals, 1.5319E+08_ireals, 2.7422E+08_ireals,  &
                             4.8708E+08_ireals, 8.5070E+08_ireals, 1.4453E+09_ireals, 2.2918E+09_ireals /)
  ltab%ltable(1,5,:, 7) = (/ 4.5093E+07_ireals, 8.0693E+07_ireals, 1.4368E+08_ireals, 2.5584E+08_ireals,  &
                             4.5052E+08_ireals, 7.8246E+08_ireals, 1.2516E+09_ireals, 1.8753E+09_ireals /)

  ! interpolated (R2=0.02mum, wcb=3.5m/s)
  ltab%ltable(1,1,:, 8) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.8386E+08_ireals,  &
                             7.1406E+08_ireals, 1.2609E+09_ireals, 2.1601E+09_ireals, 3.5372E+09_ireals /)
  ltab%ltable(1,2,:, 8) = (/ 5.0000E+07_ireals, 9.9324E+07_ireals, 1.8572E+08_ireals, 3.3744E+08_ireals,  &
                             5.9102E+08_ireals, 1.0460E+09_ireals, 1.7867E+09_ireals, 2.9481E+09_ireals /)
  ltab%ltable(1,3,:, 8) = (/ 4.9973E+07_ireals, 9.3456E+07_ireals, 1.7048E+08_ireals, 3.0508E+08_ireals,  &
                             5.4202E+08_ireals, 9.5094E+08_ireals, 1.6332E+09_ireals, 2.6962E+09_ireals /)
  ltab%ltable(1,4,:, 8) = (/ 4.7782E+07_ireals, 8.8365E+07_ireals, 1.5837E+08_ireals, 2.8384E+08_ireals,  &
                             5.0586E+08_ireals, 8.8810E+08_ireals, 1.5192E+09_ireals, 2.4464E+09_ireals /)
  ltab%ltable(1,5,:, 8) = (/ 4.6319E+07_ireals, 8.3356E+07_ireals, 1.4856E+08_ireals, 2.6478E+08_ireals,  &
                             4.6814E+08_ireals, 8.2712E+08_ireals, 1.3265E+09_ireals, 2.0176E+09_ireals /)

  ! interpolated (R2=0.02mum, wcb=4.0m/s)
  ltab%ltable(1,1,:, 9) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.8924E+08_ireals,  &
                             7.3874E+08_ireals, 1.3250E+09_ireals, 2.3018E+09_ireals, 3.8190E+09_ireals /)
  ltab%ltable(1,2,:, 9) = (/ 5.0000E+07_ireals, 9.9550E+07_ireals, 1.9048E+08_ireals, 3.4916E+08_ireals,  &
                             6.1308E+08_ireals, 1.0982E+09_ireals, 1.8905E+09_ireals, 3.1510E+09_ireals /)
  ltab%ltable(1,3,:, 9) = (/ 4.9982E+07_ireals, 9.5638E+07_ireals, 1.7628E+08_ireals, 3.1592E+08_ireals,  &
                             5.6358E+08_ireals, 9.9346E+08_ireals, 1.7184E+09_ireals, 2.8665E+09_ireals /)
  ltab%ltable(1,4,:, 9) = (/ 4.8521E+07_ireals, 9.1156E+07_ireals, 1.6356E+08_ireals, 2.9346E+08_ireals,  &
                             5.2464E+08_ireals, 9.2550E+08_ireals, 1.5932E+09_ireals, 2.6009E+09_ireals /)
  ltab%ltable(1,5,:, 9) = (/ 4.7546E+07_ireals, 8.6019E+07_ireals, 1.5344E+08_ireals, 2.7372E+08_ireals,  &
                             4.8576E+08_ireals, 8.7178E+08_ireals, 1.4013E+09_ireals, 2.1600E+09_ireals /)

  ! interpolated (R2=0.02mum, wcb=4.5m/s)
  ltab%ltable(1,1,:,10) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9462E+08_ireals,  &
                             7.6342E+08_ireals, 1.3890E+09_ireals, 2.4434E+09_ireals, 4.1007E+09_ireals /)
  ltab%ltable(1,2,:,10) = (/ 5.0000E+07_ireals, 9.9775E+07_ireals, 1.9524E+08_ireals, 3.6088E+08_ireals,  &
                             6.3514E+08_ireals, 1.1505E+09_ireals, 1.9942E+09_ireals, 3.3540E+09_ireals /)
  ltab%ltable(1,3,:,10) = (/ 4.9991E+07_ireals, 9.7819E+07_ireals, 1.8207E+08_ireals, 3.2676E+08_ireals,  &
                             5.8514E+08_ireals, 1.0360E+09_ireals, 1.8037E+09_ireals, 3.0367E+09_ireals /)
  ltab%ltable(1,4,:,10) = (/ 4.9261E+07_ireals, 9.3948E+07_ireals, 1.6874E+08_ireals, 3.0308E+08_ireals,  &
                             5.4342E+08_ireals, 9.6290E+08_ireals, 1.6671E+09_ireals, 2.7555E+09_ireals /)
  ltab%ltable(1,5,:,10) = (/ 4.8773E+07_ireals, 8.8682E+07_ireals, 1.5833E+08_ireals, 2.8266E+08_ireals,  &
                             5.0338E+08_ireals, 9.1644E+08_ireals, 1.4762E+09_ireals, 2.3023E+09_ireals /)

  ! table4d (R2=0.02mum, wcb=5.0m/s) (for Ncn=50,100,200 interpolated)
  ltab%ltable(1,1,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.8810E+08_ireals, 1.4531E+09_ireals, 2.5851E+09_ireals, 4.3825E+09_ireals /)
  ltab%ltable(1,2,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.7260E+08_ireals,  &
                             6.5720E+08_ireals, 1.2028E+09_ireals, 2.0980E+09_ireals, 3.5569E+09_ireals /)
  ltab%ltable(1,3,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.8786E+08_ireals, 3.3760E+08_ireals,  &
                             6.0670E+08_ireals, 1.0785E+09_ireals, 1.8890E+09_ireals, 3.2069E+09_ireals /)
  ltab%ltable(1,4,:,11) = (/ 5.0000E+07_ireals, 9.6739E+07_ireals, 1.7393E+08_ireals, 3.1270E+08_ireals,  &
                             5.6220E+08_ireals, 1.0003E+09_ireals, 1.7411E+09_ireals, 2.9101E+09_ireals /)
  ltab%ltable(1,5,:,11) = (/ 5.0000E+07_ireals, 9.1345E+07_ireals, 1.6321E+08_ireals, 2.9160E+08_ireals,  &
                             5.2100E+08_ireals, 9.6110E+08_ireals, 1.5511E+09_ireals, 2.4446E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=0.0m/s)
  ltab%ltable(2,1,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(2,2,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(2,3,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(2,4,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(2,5,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)

  ! table5a (R2=0.03mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 extrapolated)
  ltab%ltable(2,1,:, 2) = (/ 5.0000E+07_ireals, 9.5800E+07_ireals, 1.7620E+08_ireals, 3.2160E+08_ireals,  &
                             5.6230E+08_ireals, 8.3550E+08_ireals, 1.2414E+09_ireals, 1.8446E+09_ireals /)
  ltab%ltable(2,2,:, 2) = (/ 4.4700E+07_ireals, 8.1400E+07_ireals, 1.4450E+08_ireals, 2.5150E+08_ireals,  &
                             4.2270E+08_ireals, 6.7780E+08_ireals, 1.0869E+09_ireals, 1.7428E+09_ireals /)
  ltab%ltable(2,3,:, 2) = (/ 4.0200E+07_ireals, 7.2800E+07_ireals, 1.2930E+08_ireals, 2.2590E+08_ireals,  &
                             3.7990E+08_ireals, 6.0650E+08_ireals, 9.6826E+08_ireals, 1.5458E+09_ireals /)
  ltab%ltable(2,4,:, 2) = (/ 3.7200E+07_ireals, 6.7100E+07_ireals, 1.1950E+08_ireals, 2.0670E+08_ireals,  &
                             3.4050E+08_ireals, 5.4940E+08_ireals, 8.8646E+08_ireals, 1.4303E+09_ireals /)
  ltab%ltable(2,5,:, 2) = (/ 3.3600E+07_ireals, 5.9000E+07_ireals, 9.9400E+07_ireals, 1.5030E+08_ireals,  &
                             2.5180E+08_ireals, 4.6600E+08_ireals, 8.6241E+08_ireals, 1.5960E+09_ireals /)

  ! table5b (R2=0.03mum, wcb=1.0m/s) (for Ncn=50 interpolted and Ncn=6400 extrapolated)
  ltab%ltable(2,1,:, 3) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9760E+08_ireals, 3.5720E+08_ireals,  &
                             6.8660E+08_ireals, 1.1864E+09_ireals, 1.8922E+09_ireals, 3.0179E+09_ireals /)
  ltab%ltable(2,2,:, 3) = (/ 5.0000E+07_ireals, 9.3300E+07_ireals, 1.7220E+08_ireals, 3.1210E+08_ireals,  &
                             5.5070E+08_ireals, 9.3160E+08_ireals, 1.4766E+09_ireals, 2.3404E+09_ireals /)
  ltab%ltable(2,3,:, 3) = (/ 4.6256E+07_ireals, 8.4400E+07_ireals, 1.5400E+08_ireals, 2.7630E+08_ireals,  &
                             4.8560E+08_ireals, 8.1120E+08_ireals, 1.2717E+09_ireals, 1.9936E+09_ireals /)
  ltab%ltable(2,4,:, 3) = (/ 4.2977E+07_ireals, 7.7900E+07_ireals, 1.4120E+08_ireals, 2.5180E+08_ireals,  &
                             4.3670E+08_ireals, 7.0870E+08_ireals, 1.1177E+09_ireals, 1.7627E+09_ireals /)
  ltab%ltable(2,5,:, 3) = (/ 3.9661E+07_ireals, 7.0100E+07_ireals, 1.2390E+08_ireals, 2.1020E+08_ireals,  &
                             3.2990E+08_ireals, 5.1190E+08_ireals, 9.3340E+08_ireals, 1.7020E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=1.5m/s)
  ltab%ltable(2,1,:, 4) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9840E+08_ireals, 3.7147E+08_ireals,  &
                             7.2320E+08_ireals, 1.2989E+09_ireals, 2.1886E+09_ireals, 3.5484E+09_ireals /)
  ltab%ltable(2,2,:, 4) = (/ 5.0000E+07_ireals, 9.5533E+07_ireals, 1.8087E+08_ireals, 3.3323E+08_ireals,  &
                             5.9990E+08_ireals, 1.0424E+09_ireals, 1.7187E+09_ireals, 2.7282E+09_ireals /)
  ltab%ltable(2,3,:, 4) = (/ 4.7504E+07_ireals, 8.9219E+07_ireals, 1.6373E+08_ireals, 2.9737E+08_ireals,  &
                             5.3003E+08_ireals, 9.0920E+08_ireals, 1.4751E+09_ireals, 2.3060E+09_ireals /)
  ltab%ltable(2,4,:, 4) = (/ 4.5318E+07_ireals, 8.2723E+07_ireals, 1.5057E+08_ireals, 2.7130E+08_ireals,  &
                             4.7763E+08_ireals, 7.9970E+08_ireals, 1.2823E+09_ireals, 1.9937E+09_ireals /)
  ltab%ltable(2,5,:, 4) = (/ 4.3107E+07_ireals, 7.9530E+07_ireals, 1.3353E+08_ireals, 2.1923E+08_ireals,  &
                             3.7770E+08_ireals, 5.9897E+08_ireals, 1.0116E+09_ireals, 1.7796E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=2.0m/s)
  ltab%ltable(2,1,:, 5) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9920E+08_ireals, 3.8573E+08_ireals,  &
                             7.5980E+08_ireals, 1.4115E+09_ireals, 2.4850E+09_ireals, 4.0788E+09_ireals /)
  ltab%ltable(2,2,:, 5) = (/ 5.0000E+07_ireals, 9.7767E+07_ireals, 1.8953E+08_ireals, 3.5437E+08_ireals,  &
                             6.4910E+08_ireals, 1.1533E+09_ireals, 1.9607E+09_ireals, 3.1159E+09_ireals /)
  ltab%ltable(2,3,:, 5) = (/ 4.8752E+07_ireals, 9.4039E+07_ireals, 1.7347E+08_ireals, 3.1843E+08_ireals,  &
                             5.7447E+08_ireals, 1.0072E+09_ireals, 1.6784E+09_ireals, 2.6185E+09_ireals /)
  ltab%ltable(2,4,:, 5) = (/ 4.7659E+07_ireals, 8.7547E+07_ireals, 1.5993E+08_ireals, 2.9080E+08_ireals,  &
                             5.1857E+08_ireals, 8.9070E+08_ireals, 1.4470E+09_ireals, 2.2246E+09_ireals /)
  ltab%ltable(2,5,:, 5) = (/ 4.6554E+07_ireals, 8.8960E+07_ireals, 1.4317E+08_ireals, 2.2827E+08_ireals,  &
                             4.2550E+08_ireals, 6.8603E+08_ireals, 1.0897E+09_ireals, 1.8573E+09_ireals /)

  ! table5c (R2=0.03mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 interpolated)
  ltab%ltable(2,1,:, 6) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.9640E+08_ireals, 1.5240E+09_ireals, 2.7814E+09_ireals, 4.6093E+09_ireals /)
  ltab%ltable(2,2,:, 6) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9820E+08_ireals, 3.7550E+08_ireals,  &
                             6.9830E+08_ireals, 1.2641E+09_ireals, 2.2028E+09_ireals, 3.5036E+09_ireals /)
  ltab%ltable(2,3,:, 6) = (/ 5.0000E+07_ireals, 9.8858E+07_ireals, 1.8320E+08_ireals, 3.3950E+08_ireals,  &
                             6.1890E+08_ireals, 1.1052E+09_ireals, 1.8818E+09_ireals, 2.9309E+09_ireals /)
  ltab%ltable(2,4,:, 6) = (/ 5.0000E+07_ireals, 9.2370E+07_ireals, 1.6930E+08_ireals, 3.1030E+08_ireals,  &
                             5.5950E+08_ireals, 9.8170E+08_ireals, 1.6116E+09_ireals, 2.4556E+09_ireals /)
  ltab%ltable(2,5,:, 6) = (/ 5.0000E+07_ireals, 9.8390E+07_ireals, 1.5280E+08_ireals, 2.3730E+08_ireals,  &
                             4.7330E+08_ireals, 7.7310E+08_ireals, 1.1679E+09_ireals, 1.9350E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=3.0m/s)
  ltab%ltable(2,1,:, 7) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.9712E+08_ireals, 1.5392E+09_ireals, 2.8548E+09_ireals, 4.8450E+09_ireals /)
  ltab%ltable(2,2,:, 7) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9856E+08_ireals, 3.8040E+08_ireals,  &
                             7.1404E+08_ireals, 1.3040E+09_ireals, 2.2988E+09_ireals, 3.7395E+09_ireals /)
  ltab%ltable(2,3,:, 7) = (/ 5.0000E+07_ireals, 9.9086E+07_ireals, 1.8656E+08_ireals, 3.4748E+08_ireals,  &
                             6.3686E+08_ireals, 1.1444E+09_ireals, 1.9723E+09_ireals, 3.1351E+09_ireals /)
  ltab%ltable(2,4,:, 7) = (/ 5.0000E+07_ireals, 9.3896E+07_ireals, 1.7376E+08_ireals, 3.1870E+08_ireals,  &
                             5.7716E+08_ireals, 1.0200E+09_ireals, 1.6992E+09_ireals, 2.6276E+09_ireals /)
  ltab%ltable(2,5,:, 7) = (/ 5.0000E+07_ireals, 9.8467E+07_ireals, 1.5766E+08_ireals, 2.5336E+08_ireals,  &
                             4.9254E+08_ireals, 8.1618E+08_ireals, 1.2574E+09_ireals, 2.0341E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=3.5m/s)
  ltab%ltable(2,1,:, 8) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.9784E+08_ireals, 1.5544E+09_ireals, 2.9281E+09_ireals, 5.0807E+09_ireals /)
  ltab%ltable(2,2,:, 8) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9892E+08_ireals, 3.8530E+08_ireals,  &
                             7.2978E+08_ireals, 1.3440E+09_ireals, 2.3947E+09_ireals, 3.9754E+09_ireals /)
  ltab%ltable(2,3,:, 8) = (/ 5.0000E+07_ireals, 9.9315E+07_ireals, 1.8992E+08_ireals, 3.5546E+08_ireals,  &
                             6.5482E+08_ireals, 1.1836E+09_ireals, 2.0628E+09_ireals, 3.3393E+09_ireals /)
  ltab%ltable(2,4,:, 8) = (/ 5.0000E+07_ireals, 9.5422E+07_ireals, 1.7822E+08_ireals, 3.2710E+08_ireals,  &
                             5.9482E+08_ireals, 1.0582E+09_ireals, 1.7868E+09_ireals, 2.7996E+09_ireals /)
  ltab%ltable(2,5,:, 8) = (/ 5.0000E+07_ireals, 9.8544E+07_ireals, 1.6253E+08_ireals, 2.6942E+08_ireals,  &
                             5.1178E+08_ireals, 8.5926E+08_ireals, 1.3470E+09_ireals, 2.1331E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=4.0m/s)
  ltab%ltable(2,1,:, 9) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.9856E+08_ireals, 1.5696E+09_ireals, 3.0015E+09_ireals, 5.3165E+09_ireals /)
  ltab%ltable(2,2,:, 9) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9928E+08_ireals, 3.9020E+08_ireals,  &
                             7.4552E+08_ireals, 1.3839E+09_ireals, 2.4907E+09_ireals, 4.2112E+09_ireals /)
  ltab%ltable(2,3,:, 9) = (/ 5.0000E+07_ireals, 9.9543E+07_ireals, 1.9328E+08_ireals, 3.6344E+08_ireals,  &
                             6.7278E+08_ireals, 1.2229E+09_ireals, 2.1533E+09_ireals, 3.5434E+09_ireals /)
  ltab%ltable(2,4,:, 9) = (/ 5.0000E+07_ireals, 9.6948E+07_ireals, 1.8268E+08_ireals, 3.3550E+08_ireals,  &
                             6.1248E+08_ireals, 1.0965E+09_ireals, 1.8745E+09_ireals, 2.9716E+09_ireals /)
  ltab%ltable(2,5,:, 9) = (/ 5.0000E+07_ireals, 9.8622E+07_ireals, 1.6739E+08_ireals, 2.8548E+08_ireals,  &
                             5.3102E+08_ireals, 9.0234E+08_ireals, 1.4365E+09_ireals, 2.2322E+09_ireals /)

  ! interpolated (R2=0.03mum, wcb=4.5m/s)
  ltab%ltable(2,1,:,10) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.9928E+08_ireals, 1.5848E+09_ireals, 3.0748E+09_ireals, 5.5522E+09_ireals /)
  ltab%ltable(2,2,:,10) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9964E+08_ireals, 3.9510E+08_ireals,  &
                             7.6126E+08_ireals, 1.4239E+09_ireals, 2.5866E+09_ireals, 4.4471E+09_ireals /)
  ltab%ltable(2,3,:,10) = (/ 5.0000E+07_ireals, 9.9772E+07_ireals, 1.9664E+08_ireals, 3.7142E+08_ireals,  &
                             6.9074E+08_ireals, 1.2621E+09_ireals, 2.2438E+09_ireals, 3.7476E+09_ireals /)
  ltab%ltable(2,4,:,10) = (/ 5.0000E+07_ireals, 9.8474E+07_ireals, 1.8714E+08_ireals, 3.4390E+08_ireals,  &
                             6.3014E+08_ireals, 1.1347E+09_ireals, 1.9621E+09_ireals, 3.1436E+09_ireals /)
  ltab%ltable(2,5,:,10) = (/ 5.0000E+07_ireals, 9.8699E+07_ireals, 1.7226E+08_ireals, 3.0154E+08_ireals,  &
                             5.5026E+08_ireals, 9.4542E+08_ireals, 1.5261E+09_ireals, 2.3312E+09_ireals /)

  ! table5d (R2=0.03mum, wcb=5.0m/s) (for Ncn=50,100,200 interpolated)
  ltab%ltable(2,1,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.6000E+09_ireals, 3.1482E+09_ireals, 5.7879E+09_ireals /)
  ltab%ltable(2,2,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             7.7700E+08_ireals, 1.4638E+09_ireals, 2.6826E+09_ireals, 4.6830E+09_ireals /)
  ltab%ltable(2,3,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.7940E+08_ireals,  &
                             7.0870E+08_ireals, 1.3013E+09_ireals, 2.3343E+09_ireals, 3.9518E+09_ireals /)
  ltab%ltable(2,4,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9160E+08_ireals, 3.5230E+08_ireals,  &
                             6.4780E+08_ireals, 1.1730E+09_ireals, 2.0497E+09_ireals, 3.3156E+09_ireals /)
  ltab%ltable(2,5,:,11) = (/ 5.0000E+07_ireals, 9.8777E+07_ireals, 1.7712E+08_ireals, 3.1760E+08_ireals,  &
                             5.6950E+08_ireals, 9.8850E+08_ireals, 1.6156E+09_ireals, 2.4303E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=0.0m/s)
  ltab%ltable(3,1,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(3,2,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(3,3,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(3,4,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)
  ltab%ltable(3,5,:, 1) = (/ 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals,  &
                             0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals, 0.0000E+00_ireals /)

  ! table6a (R2=0.04mum, wcb=0.5m/s) (for Ncn=3200  and Ncn=6400 extrapolated)
  ltab%ltable(3,1,:, 2) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9650E+08_ireals, 3.7470E+08_ireals,  &
                             6.7730E+08_ireals, 1.1389E+09_ireals, 1.9151E+09_ireals, 3.2203E+09_ireals /)
  ltab%ltable(3,2,:, 2) = (/ 4.8400E+07_ireals, 9.1900E+07_ireals, 1.7060E+08_ireals, 3.0690E+08_ireals,  &
                             5.2920E+08_ireals, 8.6240E+08_ireals, 1.4054E+09_ireals, 2.2903E+09_ireals /)
  ltab%ltable(3,3,:, 2) = (/ 4.4400E+07_ireals, 8.2500E+07_ireals, 1.5030E+08_ireals, 2.6640E+08_ireals,  &
                             4.4800E+08_ireals, 7.4070E+08_ireals, 1.2246E+09_ireals, 2.0247E+09_ireals /)
  ltab%ltable(3,4,:, 2) = (/ 4.0900E+07_ireals, 7.5000E+07_ireals, 1.3470E+08_ireals, 2.3190E+08_ireals,  &
                             3.8210E+08_ireals, 6.5760E+08_ireals, 1.1317E+09_ireals, 1.9477E+09_ireals /)
  ltab%ltable(3,5,:, 2) = (/ 3.4700E+07_ireals, 5.9300E+07_ireals, 9.3500E+07_ireals, 1.5680E+08_ireals,  &
                             3.0190E+08_ireals, 6.0380E+08_ireals, 1.2076E+09_ireals, 2.4152E+09_ireals /)

  ! table6b (R2=0.04mum, wcb=1.0m/s) (for Ncn=50 interpolted and Ncn=6400 extrapolated)
  ltab%ltable(3,1,:, 3) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9880E+08_ireals,  &
                             7.7370E+08_ireals, 1.4208E+09_ireals, 2.4118E+09_ireals, 4.0940E+09_ireals /)
  ltab%ltable(3,2,:, 3) = (/ 5.0000E+07_ireals, 9.8900E+07_ireals, 1.8970E+08_ireals, 3.5620E+08_ireals,  &
                             6.4950E+08_ireals, 1.1179E+09_ireals, 1.8052E+09_ireals, 2.9151E+09_ireals /)
  ltab%ltable(3,3,:, 3) = (/ 4.9138E+07_ireals, 9.1800E+07_ireals, 1.7150E+08_ireals, 3.1490E+08_ireals,  &
                             5.5900E+08_ireals, 9.3280E+08_ireals, 1.5016E+09_ireals, 2.4172E+09_ireals /)
  ltab%ltable(3,4,:, 3) = (/ 4.6047E+07_ireals, 8.4700E+07_ireals, 1.5580E+08_ireals, 2.8050E+08_ireals,  &
                             4.8190E+08_ireals, 7.7900E+08_ireals, 1.3219E+09_ireals, 2.2432E+09_ireals /)
  ltab%ltable(3,5,:, 3) = (/ 4.1788E+07_ireals, 7.2100E+07_ireals, 1.2440E+08_ireals, 1.9840E+08_ireals,  &
                             3.1910E+08_ireals, 6.0380E+08_ireals, 1.2076E+09_ireals, 2.4152E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=1.5m/s)
  ltab%ltable(3,1,:, 4) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9920E+08_ireals,  &
                             7.8247E+08_ireals, 1.4797E+09_ireals, 2.6319E+09_ireals, 4.5273E+09_ireals /)
  ltab%ltable(3,2,:, 4) = (/ 5.0000E+07_ireals, 9.9267E+07_ireals, 1.9313E+08_ireals, 3.6947E+08_ireals,  &
                             6.8647E+08_ireals, 1.2193E+09_ireals, 2.0426E+09_ireals, 3.2976E+09_ireals /)
  ltab%ltable(3,3,:, 4) = (/ 4.9426E+07_ireals, 9.4533E+07_ireals, 1.7893E+08_ireals, 3.3237E+08_ireals,  &
                             6.0067E+08_ireals, 1.0346E+09_ireals, 1.6968E+09_ireals, 2.7072E+09_ireals /)
  ltab%ltable(3,4,:, 4) = (/ 4.7365E+07_ireals, 8.8925E+07_ireals, 1.6413E+08_ireals, 2.9890E+08_ireals,  &
                             5.2500E+08_ireals, 8.7477E+08_ireals, 1.4524E+09_ireals, 2.4222E+09_ireals /)
  ltab%ltable(3,5,:, 4) = (/ 4.4525E+07_ireals, 7.7515E+07_ireals, 1.3477E+08_ireals, 2.2350E+08_ireals,  &
                             3.6447E+08_ireals, 6.3660E+08_ireals, 1.2153E+09_ireals, 2.4280E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=2.0m/s)
  ltab%ltable(3,1,:, 5) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9960E+08_ireals,  &
                             7.9123E+08_ireals, 1.5386E+09_ireals, 2.8521E+09_ireals, 4.9606E+09_ireals /)
  ltab%ltable(3,2,:, 5) = (/ 5.0000E+07_ireals, 9.9633E+07_ireals, 1.9657E+08_ireals, 3.8273E+08_ireals,  &
                             7.2343E+08_ireals, 1.3207E+09_ireals, 2.2800E+09_ireals, 3.6802E+09_ireals /)
  ltab%ltable(3,3,:, 5) = (/ 4.9713E+07_ireals, 9.7267E+07_ireals, 1.8637E+08_ireals, 3.4983E+08_ireals,  &
                             6.4233E+08_ireals, 1.1365E+09_ireals, 1.8921E+09_ireals, 2.9971E+09_ireals /)
  ltab%ltable(3,4,:, 5) = (/ 4.8682E+07_ireals, 9.3150E+07_ireals, 1.7247E+08_ireals, 3.1730E+08_ireals,  &
                             5.6810E+08_ireals, 9.7053E+08_ireals, 1.5829E+09_ireals, 2.6013E+09_ireals /)
  ltab%ltable(3,5,:, 5) = (/ 4.7263E+07_ireals, 8.2931E+07_ireals, 1.4513E+08_ireals, 2.4860E+08_ireals,  &
                             4.0983E+08_ireals, 6.6940E+08_ireals, 1.2230E+09_ireals, 2.4409E+09_ireals /)

  ! table6c (R2=0.04mum, wcb=2.5m/s) (for Ncn=50 and Ncn=100 interpolated)
  ltab%ltable(3,1,:, 6) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.5975E+09_ireals, 3.0722E+09_ireals, 5.3939E+09_ireals /)
  ltab%ltable(3,2,:, 6) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9600E+08_ireals,  &
                             7.6040E+08_ireals, 1.4221E+09_ireals, 2.5174E+09_ireals, 4.0628E+09_ireals /)
  ltab%ltable(3,3,:, 6) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9380E+08_ireals, 3.6730E+08_ireals,  &
                             6.8400E+08_ireals, 1.2383E+09_ireals, 2.0873E+09_ireals, 3.2871E+09_ireals /)
  ltab%ltable(3,4,:, 6) = (/ 5.0000E+07_ireals, 9.7375E+07_ireals, 1.8080E+08_ireals, 3.3570E+08_ireals,  &
                             6.1120E+08_ireals, 1.0663E+09_ireals, 1.7134E+09_ireals, 2.7803E+09_ireals /)
  ltab%ltable(3,5,:, 6) = (/ 5.0000E+07_ireals, 8.8346E+07_ireals, 1.5550E+08_ireals, 2.7370E+08_ireals,  &
                             4.5520E+08_ireals, 7.0220E+08_ireals, 1.2307E+09_ireals, 2.4537E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=3.0m/s)
  ltab%ltable(3,1,:, 7) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.5980E+09_ireals, 3.0978E+09_ireals, 5.5639E+09_ireals /)
  ltab%ltable(3,2,:, 7) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9680E+08_ireals,  &
                             7.6832E+08_ireals, 1.4492E+09_ireals, 2.6020E+09_ireals, 4.2923E+09_ireals /)
  ltab%ltable(3,3,:, 7) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9504E+08_ireals, 3.7318E+08_ireals,  &
                             6.9830E+08_ireals, 1.2735E+09_ireals, 2.1829E+09_ireals, 3.4873E+09_ireals /)
  ltab%ltable(3,4,:, 7) = (/ 5.0000E+07_ireals, 9.7900E+07_ireals, 1.8456E+08_ireals, 3.4294E+08_ireals,  &
                             6.2754E+08_ireals, 1.1054E+09_ireals, 1.8084E+09_ireals, 2.9165E+09_ireals /)
  ltab%ltable(3,5,:, 7) = (/ 5.0000E+07_ireals, 9.0677E+07_ireals, 1.6072E+08_ireals, 2.8284E+08_ireals,  &
                             4.7650E+08_ireals, 7.5254E+08_ireals, 1.2833E+09_ireals, 2.4559E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=3.5m/s)
  ltab%ltable(3,1,:, 8) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.5985E+09_ireals, 3.1233E+09_ireals, 5.7339E+09_ireals /)
  ltab%ltable(3,2,:, 8) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9760E+08_ireals,  &
                             7.7624E+08_ireals, 1.4762E+09_ireals, 2.6866E+09_ireals, 4.5217E+09_ireals /)
  ltab%ltable(3,3,:, 8) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9628E+08_ireals, 3.7906E+08_ireals,  &
                             7.1260E+08_ireals, 1.3088E+09_ireals, 2.2785E+09_ireals, 3.6875E+09_ireals /)
  ltab%ltable(3,4,:, 8) = (/ 5.0000E+07_ireals, 9.8425E+07_ireals, 1.8832E+08_ireals, 3.5018E+08_ireals,  &
                             6.4388E+08_ireals, 1.1446E+09_ireals, 1.9034E+09_ireals, 3.0527E+09_ireals /)
  ltab%ltable(3,5,:, 8) = (/ 5.0000E+07_ireals, 9.3007E+07_ireals, 1.6595E+08_ireals, 2.9198E+08_ireals,  &
                             4.9780E+08_ireals, 8.0288E+08_ireals, 1.3360E+09_ireals, 2.4581E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=4.0m/s)
  ltab%ltable(3,1,:, 9) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.5990E+09_ireals, 3.1489E+09_ireals, 5.9039E+09_ireals /)
  ltab%ltable(3,2,:, 9) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9840E+08_ireals,  &
                             7.8416E+08_ireals, 1.5033E+09_ireals, 2.7712E+09_ireals, 4.7512E+09_ireals /)
  ltab%ltable(3,3,:, 9) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9752E+08_ireals, 3.8494E+08_ireals,  &
                             7.2690E+08_ireals, 1.3440E+09_ireals, 2.3741E+09_ireals, 3.8877E+09_ireals /)
  ltab%ltable(3,4,:, 9) = (/ 5.0000E+07_ireals, 9.8950E+07_ireals, 1.9209E+08_ireals, 3.5742E+08_ireals,  &
                             6.6022E+08_ireals, 1.1837E+09_ireals, 1.9983E+09_ireals, 3.1888E+09_ireals /)
  ltab%ltable(3,5,:, 9) = (/ 5.0000E+07_ireals, 9.5338E+07_ireals, 1.7117E+08_ireals, 3.0112E+08_ireals,  &
                             5.1910E+08_ireals, 8.5322E+08_ireals, 1.3886E+09_ireals, 2.4603E+09_ireals /)

  ! interpolated (R2=0.04mum, wcb=4.5m/s)
  ltab%ltable(3,1,:,10) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.5995E+09_ireals, 3.1744E+09_ireals, 6.0739E+09_ireals /)
  ltab%ltable(3,2,:,10) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9920E+08_ireals,  &
                             7.9208E+08_ireals, 1.5303E+09_ireals, 2.8558E+09_ireals, 4.9806E+09_ireals /)
  ltab%ltable(3,3,:,10) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9876E+08_ireals, 3.9082E+08_ireals,  &
                             7.4120E+08_ireals, 1.3793E+09_ireals, 2.4697E+09_ireals, 4.0879E+09_ireals /)
  ltab%ltable(3,4,:,10) = (/ 5.0000E+07_ireals, 9.9475E+07_ireals, 1.9585E+08_ireals, 3.6466E+08_ireals,  &
                             6.7656E+08_ireals, 1.2229E+09_ireals, 2.0933E+09_ireals, 3.3250E+09_ireals /)
  ltab%ltable(3,5,:,10) = (/ 5.0000E+07_ireals, 9.7669E+07_ireals, 1.7640E+08_ireals, 3.1026E+08_ireals,  &
                             5.4040E+08_ireals, 9.0356E+08_ireals, 1.4413E+09_ireals, 2.4625E+09_ireals /)

  ! table6d (R2=0.04mum, wcb=5.0m/s) (for Ncn=50,100,200 interpolated)
  ltab%ltable(3,1,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.6000E+09_ireals, 3.2000E+09_ireals, 6.2439E+09_ireals /)
  ltab%ltable(3,2,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 4.0000E+08_ireals,  &
                             8.0000E+08_ireals, 1.5574E+09_ireals, 2.9404E+09_ireals, 5.2101E+09_ireals /)
  ltab%ltable(3,3,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 2.0000E+08_ireals, 3.9670E+08_ireals,  &
                             7.5550E+08_ireals, 1.4145E+09_ireals, 2.5653E+09_ireals, 4.2881E+09_ireals /)
  ltab%ltable(3,4,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.9961E+08_ireals, 3.7190E+08_ireals,  &
                             6.9290E+08_ireals, 1.2620E+09_ireals, 2.1883E+09_ireals, 3.4612E+09_ireals /)
  ltab%ltable(3,5,:,11) = (/ 5.0000E+07_ireals, 1.0000E+08_ireals, 1.8162E+08_ireals, 3.1940E+08_ireals,  &
                             5.6170E+08_ireals, 9.5390E+08_ireals, 1.4939E+09_ireals, 2.4647E+09_ireals /)

END SUBROUTINE nccn_lookupcreate_segalkhain_4D

!===========================================================================================

SUBROUTINE nccn_lookupcreate_segalkhain_2D ( ltab4D, lsigs, r2, ltab2D )

  IMPLICIT NONE

  TYPE(lookupt_4D), INTENT(in  ) :: ltab4D
  REAL(kind=ireals), INTENT(in)  :: lsigs, r2

  TYPE(lookupt_2D), INTENT(out)  :: ltab2D

  INTEGER(KIND=iintegers) :: i, j
  REAL(kind=ireals)  :: x1loc(ltab4D%n3)


  ltab2D%n1 = ltab4D%n3    ! for log(ncn)
  ltab2D%n2 = ltab4D%n4    ! for wcb

  
  ALLOCATE( ltab2D%x1(ltab2D%n1) )
  ALLOCATE( ltab2D%x2(ltab2D%n2) )
  ALLOCATE( ltab2D%ltable(ltab2D%n1,ltab2D%n2) )

  ! N_CN will be linearily interpolated in log-log-space:
  ! ncn in 1/m^3 : ltab2D%x1 = LOG(/  5.000E+07, 1.000E+08, 2.000E+08, 4.000E+08, 8.000E+08, 1.600E+09, 3.200E+09, 6.400E+09 /)
  ! wcb in m/s   : ltab2D%x2 = (/  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 /)

  ltab2D%x1 = ltab4D%x3
  ltab2D%x2 = ltab4D%x4

  ltab2D%dx1 = ltab4D%dx3
  ltab2D%dx2 = ltab4D%dx4

  ltab2D%odx1 = ltab4D%odx3
  ltab2D%odx2 = ltab4D%odx4

  ! .. Generate 2D table by interpolating w.r.t. the other two dimensions (lsigs, r2)
  DO j=1, ltab2D%n2
    DO i=1, ltab2D%n1

      CALL interpol_nccn_segalkhain_4D(ltab4D, r2, lsigs, EXP(ltab2D%x1(i)), ltab2D%x2(j), ltab2D%ltable(i,j))

    END DO
  END DO


END SUBROUTINE nccn_lookupcreate_segalkhain_2D

!===========================================================================================

SUBROUTINE lookupcreate_segalkhain

  ! Wrapper for init calls to create lookup tables for Segal-Khain parameterization
  CALL nccn_lookupcreate_segalkhain_4D ( ltab4D )
  IF (r2_lsigs_are_fixed) CALL nccn_lookupcreate_segalkhain_2D ( ltab4D, lsigs_fix, r2_fix, ltab2D )

END SUBROUTINE lookupcreate_segalkhain

!===========================================================================================

SUBROUTINE interpol_nccn_segalkhain_4D(ltab, r2, lsigs, ncn, wcb, nccn)

  IMPLICIT NONE

  ! Inputs:
  TYPE(lookupt_4D),  INTENT(in)  :: ltab
  REAL(kind=ireals), INTENT(in)  :: ncn, wcb, lsigs, r2

  ! Outputs:
  REAL(kind=ireals), INTENT(out) :: nccn

  ! Local variables:
  INTEGER(KIND=iintegers) :: iu, ju, ku, lu
  REAL(kind=ireals)       :: r2_loc, lsigs_loc, lncn_loc, wcb_loc, &
                             hilf1(2,2,2,2), hilf2(2,2,2), hilf3(2,2), hilf4(2)

  REAL(KIND=ireals),    PARAMETER :: eps  = 1e-8_ireals

  ! For information:
  ! r2 im mu-m   : ltab%x1 = (/  0.02, 0.03, 0.04 /)
  ! lsigs        : ltab%x2 = (/  0.1, 0.2, 0.3, 0.4, 0.5 /)
  ! ncn in 1/m^3 : ltab%x3 = LOG(/  5.000E+07, 1.000E+08, 2.000E+08, 4.000E+08, 8.000E+08, 1.600E+09, 3.200E+09, 6.400E+09 /)
  ! wcb in m/s   : ltab%x4 = (/  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 /)

  ! Interpolation of the look-up tables with respect to all 4 parameters:
  ! (clip values outside range to the marginal values where appropriate, otherwise extrapolate linearily)
  !
  ! .. r2 is clipped:
  r2_loc    = MIN(MAX(r2,     ltab%x1(1)), ltab%x1(ltab%n1))
  iu = MIN(FLOOR((r2_loc -    ltab%x1(1)) * ltab%odx1 ) + 1, ltab%n1-1)
  ! .. log(sigma_s) is clipped:
  lsigs_loc = MIN(MAX(lsigs,  ltab%x2(1)), ltab%x2(ltab%n2))
  ju = MIN(FLOOR((lsigs_loc - ltab%x2(1)) * ltab%odx2 ) + 1, ltab%n2-1)
  ! .. log(ncn) is not clipped, so that logarithmic extrapolation results for values outside table range:
  lncn_loc   = MAX(LOG(ncn),   0.0_ireals)

  ku = MAX(MIN(FLOOR((lncn_loc -   ltab%x3(1)) * ltab%odx3 ) + 1, ltab%n3-1), 1)
  ! .. wcb is only clipped at the lower table end:
  wcb_loc   = MAX(wcb,    ltab%x4(1))
  lu = MIN(FLOOR((wcb_loc -   ltab%x4(1)) * ltab%odx4 ) + 1, ltab%n4-1)

  hilf1 = MAX(ltab%ltable( iu:iu+1, ju:ju+1, ku:ku+1, lu:lu+1), eps)
  hilf2 = EXP( LOG(hilf1(:,:,1,:)) + LOG(hilf1(:,:,2,:)/hilf1(:,:,1,:)) * ltab%odx3 * ( lncn_loc-ltab%x3(ku)) )
  hilf3 = hilf2(1,:,:)   + (hilf2(2,:,:)   - hilf2(1,:,:)  ) * ltab%odx1 * ( r2_loc    - ltab%x1(iu) )
  hilf4 = hilf3(1,:)     + (hilf3(2,:)     - hilf3(1,:)    ) * ltab%odx2 * ( lsigs_loc - ltab%x2(ju) )
  nccn  = hilf4(1)       + (hilf4(2)       - hilf4(1)      ) * ltab%odx4 * ( wcb_loc   - ltab%x4(lu) )

  ! For safety, clip nccn to ncn, because it could become larger by extrapolation
  !  to values of ncn or wcb that are larger than the upper table range:
  nccn = MIN(nccn, ncn)

  IF (nccn < 2.0*eps) nccn = 0.0

END SUBROUTINE interpol_nccn_segalkhain_4D

!===========================================================================================

SUBROUTINE interpol_nccn_segalkhain_2D(ltab, ncn, wcb, nccn)

  IMPLICIT NONE

  ! Inputs:
  TYPE(lookupt_2D),  INTENT(in)  :: ltab
  REAL(kind=ireals), INTENT(in)  :: ncn, wcb

  ! Outputs:
  REAL(kind=ireals), INTENT(out) :: nccn

  ! Local variables:
  INTEGER(KIND=iintegers) :: lu, ku
  REAL(kind=ireals)       :: lncn_loc, wcb_loc, &
                             hilf1(2,2), hilf2(2)

  REAL(KIND=ireals),    PARAMETER :: eps  = 1e-8_ireals

  ! For information:
  ! ncn in 1/m^3 : ltab%x1 = LOG(/  5.000E+07, 1.000E+08, 2.000E+08, 4.000E+08, 8.000E+08, 1.600E+09, 3.200E+09, 6.400E+09 /)
  ! wcb in m/s   : ltab%x2 = (/  0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 /)

  ! Interpolation of the look-up tables with respect to 2 parameters (ncn and wcb):
  ! (clip values outside range to the marginal values where appropriate, otherwise extrapolate linearily)
  !
  ! .. log(ncn) is not clipped, so that logarithmic extrapolation results for values outside table range:
  lncn_loc   = MAX(LOG(ncn),   0.0_ireals)
  ku = MAX(MIN(FLOOR((lncn_loc -   ltab%x1(1)) * ltab%odx1 ) + 1, ltab%n1-1), 1)
  ! .. wcb is only clipped at the lower table end:
  wcb_loc   = MAX(wcb,    ltab%x2(1))
  lu = MIN(FLOOR((wcb_loc -   ltab%x2(1)) * ltab%odx2 ) + 1, ltab%n2-1)

  hilf1 = MAX(ltab%ltable( ku:ku+1, lu:lu+1), eps)
  hilf2 = EXP( LOG(hilf1(1,:)) + LOG(hilf1(2,:)/hilf1(1,:)) * ltab%odx1 * ( lncn_loc-ltab%x1(ku)) )
  nccn  = hilf2(1)             + (hilf2(2) - hilf2(1)     ) * ltab%odx2 * ( wcb_loc -ltab%x2(lu))

  ! For safety, clip nccn to ncn, because it could become larger by extrapolation
  !  to values of ncn or wcb that are larger than the upper table range:
  nccn = MIN(nccn, ncn)

  IF (nccn < 2.0_ireals*eps) nccn = 0.0_ireals

END SUBROUTINE interpol_nccn_segalkhain_2D

!===========================================================================================


SUBROUTINE specccn_segalkhain ( ie, ke, istart, iend, kstart, kend, ncn, w, qc, rho, hhl, cloud_num, r2, lsigs)

  ! ncn in 1/m^3, as resuling from subroutine interpol_nccn_segalkhain()
  ! qc in kg/kg
  ! cloud_num in 1/kg

  ! ncn should already contain the effects of the vertical profile!
  ! hhl are the heights of the model half levels, i.e., the points where w is defined.
  !
  ! w_cb is diagnosed from w and rhoc
  ! then an N_CCN is interpolated at w_cb, r2, lsigs and ncn at the height of the cloud base
  ! then a further exponential decrease with height of the cloud_num is prescribed
  ! If there is qc but w < 0, we set w to 0.1 m/s to avoid very small values of cloud_num

  IMPLICIT NONE
  
  ! Inputs:
  INTEGER(kind=iintegers), INTENT(in)                 :: ie, ke, istart, iend, kstart, kend
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie,ke)    :: ncn, rho, qc, hhl
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie,ke+1)  :: w
  REAL(kind=ireals), INTENT(in) , OPTIONAL            :: r2, lsigs ! needed if r2_lsigs_are_fixed=.FALSE.
  ! Outputs:
  REAL(kind=ireals), INTENT(out), DIMENSION(ie,ke)  :: cloud_num

  ! Local variables:
  INTEGER(kind=iintegers) :: i, k, kp1, kcb(ie,ke)
  REAL(kind=ireals)       :: wcb(ie,ke), zcb(ie,ke), z(ie,ke), cloud_num_min_loc
  LOGICAL                 :: found


  INTEGER(KIND=iintegers), PARAMETER :: itype_incloud_profile = 2
  REAL(KIND=ireals),    PARAMETER :: eps  = 1e-15_ireals
  REAL(KIND=ireals),    PARAMETER :: wcb_min  = 0.05_ireals
  REAL(KIND=ireals),    PARAMETER :: cloud_num_min  = 10.0e6_ireals  ! 1/kg

  IF (.NOT. r2_lsigs_are_fixed) THEN
    IF (.NOT. PRESENT(r2) .OR. .NOT. PRESENT(lsigs)) THEN
      CALL finish('specccn_segalkhain','argument r2 or lsigs is missing')
    ENDIF
  ENDIF

  wcb = 0.0_ireals
  zcb = 0.0_ireals
  kcb = 0_ireals
  cloud_num = 0.0_ireals
  z = 0.5 * (hhl(:,2:ke+1) + hhl(:,1:ke))

  IF (itype_incloud_profile == 1) THEN

    IF (r2_lsigs_are_fixed) THEN

      DO k = kend, kstart, -1
        DO i = istart, iend
          
          kp1 = MIN(k+1,ke)
          
          ! Sehr einfache Alternative:
          IF (qc(i,k) >= eps) THEN
            ! other cloudy grid points, which are not within a vertically "active" cloud:
            CALL interpol_nccn_segalkhain_2D(ltab2D, ncn(i,k), MAX(w(i,k),wcb_min), cloud_num(i,k))
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min*rho(i,k))
          END IF
          
        END DO
      END DO

    ELSE

      DO k = kend, kstart, -1
        DO i = istart, iend
          
          kp1 = MIN(k+1,ke)
          
          ! Sehr einfache Alternative:
          IF (qc(i,k) >= eps) THEN
            ! other cloudy grid points, which are not within a vertically "active" cloud:
            CALL interpol_nccn_segalkhain_4D(ltab4D, ncn(i,k), MAX(w(i,k),wcb_min), lsigs, r2, cloud_num(i,k))
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min*rho(i,k))
          END IF
          
        END DO
      END DO

    END IF

  ELSE IF (itype_incloud_profile == 2) THEN 

    ! Determine the w at true cloud bases, wcb, by scanning from bottom to top,
    ! interpolate NCCN from the lookup table.
    ! Hereby, all continuous cloudy layers above a
    ! cloud base with an updraft stronger than wcb_min are assumed to be vertically active
    ! clouds, and here an  exponentially decreasing profile NCCN(z) is assumed, depending
    ! on the vertical distance to the cloud base. If at some heights the vertical velocity
    ! is larger than that at cloud base, incloud nucleation is parameterized in a way,
    ! that the actual w is assumed to have been present already at cloud base (new NCCN-interpolation)
    ! and the exponential decrease is imposed as usual as function of distance to cloud base.
    !
    ! In all other (non-continuous) cloud layers, interpolate NCN from the NCCN-table
    ! assuming a w at cloud base of the value wcb_min.

    IF (r2_lsigs_are_fixed) THEN

      DO k = kend, kstart, -1
        DO i = istart, iend
      
          kp1 = MIN(k+1,ke)
          cloud_num_min_loc = cloud_num_min * rho(i,k)
          
          IF (w(i,kp1) >= wcb_min .AND. qc(i,k) >= eps .AND. (k == ke .OR. qc(i,kp1) < eps ) ) THEN
            ! either saturated updraft at cloud base or in the lowest model level:
            wcb(i,k) = w(i,kp1)
            zcb(i,k) = hhl(i,k+1)
            kcb(i,k) = k
            CALL interpol_nccn_segalkhain_2D(ltab2D, ncn(i,k), wcb(i,k), cloud_num(i,k))
            cloud_num(i,k) = cloud_num(i,k) * height_profile_simple(zcb(i,k),z(i,k))
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min_loc)
          ELSE IF (qc(i,k) >= eps .AND. k < ke .AND. kcb(i,kp1) > 0 .AND. wcb(i,kp1) >= wcb_min) THEN
            ! saturated grid point within a continuous cloud layer above a well-defined cloud base and above the lowest model level:
            wcb(i,k) = wcb(i,kp1)
            zcb(i,k) = zcb(i,kp1)
            kcb(i,k) = kcb(i,kp1)
            IF (lincloud .AND. w(i,kp1) > wcb(i,k)) THEN
              ! In this case, we assume additional incloud nucleation in an assumed coherent updraft.
              !  We parameterize
              !  the incloud nucleation in a simple way by calculating the cloud_num which would have resulted
              !  if the higher updraft speed would have been already present at cloud base
              !  (ncn from cloud base air), and by imposing the relative vertical profile starting
              !  from the half height between cloud base and the actual height.
              ! As alternative, one could have been calculated the difference of NCCN(w(k),ncn_cb) - NCCN(w_cb,ncn_cb)
              !  and add it to the existing cloud_num(k), maybe after superimposing a vertical profile
              !  on the difference from hhl(k) to z(k).
              ! However, in real clouds incloud nucleation is very unlikely because there is no
              !  secondary supersaturation maximum larger than the cloud base maximum, even if the
              !  w does increase with height, because the condensation rate ~ N*r_quer compensates for it.
              CALL interpol_nccn_segalkhain_2D(ltab2D, ncn(i,kcb(i,k)), w(i,kp1), cloud_num(i,k))
              cloud_num(i,k) = cloud_num(i,k) * height_profile_simple(0.5*(zcb(i,k)+z(i,k)),z(i,k))
              wcb(i,k) = w(i,kp1)
            ELSE
              cloud_num(i,k) = cloud_num(i,kp1) * height_profile_simple(z(i,kp1),z(i,k))
            END IF
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min_loc)
          ELSE IF (qc(i,k) >= eps) THEN
            ! other cloudy grid points, which are not within a vertically "active" cloud:
            CALL interpol_nccn_segalkhain_2D(ltab2D, ncn(i,k), MAX(w(i,k),wcb_min), cloud_num(i,k))
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min_loc)
          END IF
          
        END DO
      END DO

    ELSE

      ! The same again, but with full 4D table interpolation:

      DO k = kend, kstart, -1
        DO i = istart, iend
      
          kp1 = MIN(k+1,ke)
          cloud_num_min_loc = cloud_num_min * rho(i,k)
          
          IF (w(i,kp1) >= wcb_min .AND. qc(i,k) >= eps .AND. (k == ke .OR. qc(i,kp1) < eps ) ) THEN
            ! either saturated updraft at cloud base or in the lowest model level:
            wcb(i,k) = w(i,kp1)
            zcb(i,k) = hhl(i,k+1)
            kcb(i,k) = k
            CALL interpol_nccn_segalkhain_4D(ltab4D, r2, lsigs, ncn(i,k), wcb(i,k), cloud_num(i,k))
            cloud_num(i,k) = cloud_num(i,k) * height_profile_simple(zcb(i,k),z(i,k))
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min_loc)
          ELSE IF (qc(i,k) >= eps .AND. k < ke .AND. kcb(i,kp1) > 0 .AND. wcb(i,kp1) >= wcb_min) THEN
            ! saturated grid point within a continuous cloud layer above a well-defined cloud base and above the lowest model level:
            !  (otherwise, there would have been voids in the vertical qc-structure)
            wcb(i,k) = wcb(i,kp1)
            zcb(i,k) = zcb(i,kp1)
            kcb(i,k) = kcb(i,kp1)
            IF (lincloud .AND. w(i,kp1) > wcb(i,k)) THEN
              ! In this case, we assume additional incloud nucleation in an assumed coherent updraft.
              !  We parameterize
              !  the incloud nucleation in a simple way by calculating the cloud_num which would have resulted
              !  if the higher updraft speed would have been already present at cloud base
              !  (ncn from cloud base air), and by imposing the relative vertical profile starting
              !  from the half height between cloud base and the actual height.
              ! As alternative, one could have been calculated the difference of NCCN(w(k),ncn_cb) - NCCN(w_cb,ncn_cb)
              !  and add it to the existing cloud_num(k), maybe after superimposing a vertical profile
              !  on the difference from hhl(k) to z(k).
              ! However, in real clouds incloud nucleation is very unlikely because there is no
              !  secondary supersaturation maximum larger than the cloud base maximum, even if the
              !  w does increase with height, because the condensation rate ~ N*r_quer compensates for it.
              CALL interpol_nccn_segalkhain_4D(ltab4D, r2, lsigs, ncn(i,kcb(i,k)), w(i,kp1), cloud_num(i,k))
              cloud_num(i,k) = cloud_num(i,k) * height_profile_simple(0.5*(zcb(i,k)+z(i,k)),z(i,k))
              wcb(i,k) = w(i,kp1)
            ELSE
              cloud_num(i,k) = cloud_num(i,kp1) * height_profile_simple(z(i,kp1),z(i,k))
            END IF
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min_loc)
          ELSE IF (qc(i,k) >= eps) THEN
            ! other cloudy grid points, which are not within a vertically "active" cloud:
            CALL interpol_nccn_segalkhain_4D(ltab4D, r2, lsigs, ncn(i,k), MAX(w(i,k),wcb_min), cloud_num(i,k))
            cloud_num(i,k) = MAX(cloud_num(i,k), cloud_num_min_loc)
          END IF
          
        END DO
      END DO

    END IF

  END IF

  ! up to now, cloud_num is in 1/m^3. Convert to 1/kg and impose global lower limit:
  cloud_num = MAX(cloud_num/rho, cloud_num_min)

  ! IDEA: Instead of this, limit Dmean ~ (qc/cloud_num)^(1/3) to the range 5 um - 100 um

CONTAINS

  FUNCTION height_profile_simple (z1, z2) RESULT (f)
    REAL(KIND=ireals), INTENT(in) :: z1, z2
    REAL(KIND=ireals)             :: f
    REAL(KIND=ireals), PARAMETER  :: z1oe0 = 2000.0_ireals ! m
    f = MIN(EXP((z1-z2)/z1oe0), 1.0_ireals)
  END FUNCTION height_profile_simple

  ! Alternative formulation where the height decrease of cloud_num
  !  depends on w. The larger w, the slower the decrease because
  !  autoconversion and cloud selfcollection have less time to act.
  FUNCTION height_profile (z1, z2, w) RESULT (f)
    REAL(KIND=ireals), INTENT(in) :: z1, z2, w
    REAL(KIND=ireals)             :: f
    REAL(KIND=ireals)             :: z1oe
    REAL(KIND=ireals), PARAMETER  :: z1oe0 = 2000.0_ireals ! m
    REAL(KIND=ireals), PARAMETER  :: w0    = 1.0_ireals    ! m/s
    REAL(KIND=ireals), PARAMETER  :: c0    = z1oe0 / w0
    ! The 1/e-height of the exponential decrease scales linearily with w.
    !  This results, if the cloud_num-change by autoconversion and cloud selfcollection is assumed to
    !  be proportional to -cloud_num. This assumption is not entirely correct,
    !  but serves as a simple first-order approximation.
    ! c0 is a reference value, i.e., if the updraft is 1 m/s, the 1/e-height
    !  is assumed to be 2000 m.
    z1oe = w * c0
    f = MIN(EXP((z1-z2)/z1oe), 1.0_ireals)
  END FUNCTION height_profile

END SUBROUTINE specccn_segalkhain
  

SUBROUTINE specccn_segalkhain_simple (ie, istart, iend, ncn, cloud_num)

  ! highly simplified version of above routine disregarding the dependency on height and vertical wind speed

  IMPLICIT NONE
  
  ! Inputs:
  INTEGER(kind=iintegers), INTENT(in)            :: ie, istart, iend
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie)  :: ncn

  ! Outputs:
  REAL(kind=ireals), INTENT(out), DIMENSION(ie)  :: cloud_num

  ! Local variables:
  INTEGER(kind=iintegers) :: i

  REAL(KIND=ireals),    PARAMETER :: wcb  = 0.25_ireals ! assume wind speed of 25 cm7s at cloud base
  REAL(KIND=ireals),    PARAMETER :: cloud_num_min  = 10.0e6_ireals  ! 1/kg


  DO i = istart, iend
          
    CALL interpol_nccn_segalkhain_2D(ltab2D, ncn(i), wcb, cloud_num(i))
    cloud_num(i) = MAX(cloud_num(i), cloud_num_min)
          
  END DO

END SUBROUTINE specccn_segalkhain_simple

!===========================================================================================

! returns ncn as a partial number density, unit 1/m^3

SUBROUTINE ncn_from_tau_aerosol ( ie, ke, hhl, aer_ss, aer_so4, aer_org, aer_dust, &
     ncn )

  IMPLICIT NONE

  !-------------------------------------------------------------------------------------------
  !.. Inputs:

  INTEGER(kind=iintegers), INTENT(in)                 :: ie, ke
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie,ke+1)  :: hhl
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie)       :: aer_ss, aer_so4, aer_org, aer_dust


  !-------------------------------------------------------------------------------------------
  !.. Outputs:

  REAL(kind=ireals), INTENT(out), DIMENSION(ie,ke)    :: ncn

  !-------------------------------------------------------------------------------------------
  !.. Local variables:

  INTEGER(kind=iintegers) :: i, k
  REAL(kind=ireals)       :: z, z0(ie), nscale(ie)

  !-------------------------------------------------------------------------------------------
  !.. Local parameters:

  REAL(KIND=ireals), PARAMETER :: pi = 3.141592653589793

  !   Assumptions on the specific extinction coefficient from Tegen et al. (1997), Table 1:
  REAL(KIND=ireals),    PARAMETER :: &
       beta_ext_dust = 1.5d3, &    ! m^2/kg
       beta_ext_org  = 8.0d3, &
       beta_ext_so4  = 8.0d3, &
       beta_ext_ss   = 0.2d3

  !   Estimated value of the aerosol number concentation,
  !    assuming a certain mean mass radius of the aerosols and a certain bulk density:
  REAL(KIND=ireals),    PARAMETER :: &
       rho_ss   = 3000.0_ireals, &    ! aerosol bulk density kg/m^3 for ss
       r_ss     = 500e-9_ireals, &    ! aerosol mean mass radius in m for ss  (200 - 12000 nm)
       rho_so4  = 2000.0_ireals, &    ! aerosol bulk density kg/m^3 for so4
       r_so4    = 80e-9_ireals, &     ! aerosol mean mass radius in m for so4 (50 - 100 nm)
       rho_org  = 2000.0_ireals, &    ! aerosol bulk density kg/m^3 for organics
       r_org    = 80e-9_ireals, &     ! aerosol mean mass radius in m for organics (50 - 100 nm)
       rho_dust = 3000.0_ireals, &    ! aerosol bulk density kg/m^3 for dust
       r_dust   = 1000e-9_ireals      ! aerosol mean mass radius in m for dust
  
  !   Soluble_fraction_dust: comparatively low number  ~0.0 - 0.1
  !   Soluble_fraction_organics: comparatively high number ~0.9 - 1.0
  REAL(KIND=ireals),    PARAMETER :: &
       soluble_fraction_dust     = 0.1_ireals, &
       soluble_fraction_organics = 0.9_ireals

  !   Parameters of assumed vertical exponential profile:
  !    ztrans = transition height in m, constant value below, expon. profile above
  !    z1oe = 1/e - height of expon. decrease above z0
  REAL(KIND=ireals),    PARAMETER :: &
       ztrans = 3000.0_ireals, &
       z1oe   = 2000.0_ireals, &
       zmax   = 12000.0_ireals

  !-------------------------------------------------------------------------------------------

  ! In the parameterization of Tegen et al.,
  ! aer_org (organics) contains aer_bc (black carbon), so that
  ! only aer_org is used in the following.

  ! Total vertically integral aerosol number per m^2 from the optical thicknesses,
  !  assuming that the mean mass radius of each species is a constant everywhere:

  nscale(:) = 3.0 / (4.0*pi) * ( &
       aer_dust(:)/(beta_ext_dust*rho_dust*r_dust**3) * soluble_fraction_dust     + &
       aer_org(:)/(beta_ext_org*rho_org*r_org**3)     * soluble_fraction_organics + &
       aer_so4(:)/(beta_ext_so4*rho_so4*r_so4**3)     + &
       aer_ss(:)/(beta_ext_ss*rho_ss*r_ss**3)           &
       )
           
  ! From that, compute the value of the aerosol mass density in the surface layer,
  !  which serves as a scaling parameter in the vertical distribution:
  !  This value is computed under the assumptions that:
  !  - the vertical profile of the aerosol mass density is constant in a
  !    surface layer up to 3000 m MSL and decreases exponentially above,
  !  - the aerosols in the Tegen climatology are confined within the height layer from the ground
  !    to a zmax, which we assume to be 12 km.

  z0(:) = MAX(ztrans, hhl(:,ke+1))
  nscale(:) = nscale(:) / (z0(:) - hhl(:,ke+1) + z1oe*(1.0_ireals-EXP((z0(:)-zmax)/z1oe)))

  ! Now compute ncn as function of height from nscale and the vertical profile:
  DO k=1, ke
    DO i=1,ie
      z = 0.5 * (hhl(i,k) + hhl(i,k+1))
      IF (z <= z0(i)) THEN
        ncn(i,k) = nscale(i)
      ELSE
        ncn(i,k) = nscale(i) * EXP((z0(i)-z)/z1oe)
      END IF
    END DO
  END DO

END SUBROUTINE ncn_from_tau_aerosol

!===========================================================================================

! Alternative version which assumes an exponential profile and constant surface layer value
!  for the specific aerosol mass, not the aerosol partial density as in ncn_from_tau_aerosol.
! Also returns ncn as a partial number density, unit 1/m^3.

SUBROUTINE ncn_from_tau_aerosol_speccnconst ( ie, ke, istart, iend, kstart, kend, hhl, &
     aer_ss, aer_so4, aer_org, aer_dust, ncn )

  IMPLICIT NONE

  !-------------------------------------------------------------------------------------------
  !.. Inputs:

  INTEGER(kind=iintegers), INTENT(in)                 :: ie, ke, istart, iend, kstart, kend
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie,ke+1)  :: hhl
  REAL(kind=ireals), INTENT(in) , DIMENSION(ie)       :: aer_ss, aer_so4, aer_org, aer_dust


  !-------------------------------------------------------------------------------------------
  !.. Outputs:

  REAL(kind=ireals), INTENT(out), DIMENSION(ie,ke)    :: ncn

  !-------------------------------------------------------------------------------------------
  !.. Local variables:

  INTEGER(kind=iintegers) :: i, k
  REAL(kind=ireals)       :: z, z0(ie), rho_air_0(ie), nscale(ie)

  !-------------------------------------------------------------------------------------------
  !.. Local parameters:

  REAL(KIND=ireals), PARAMETER :: pi = 3.141592653589793

  !   Assumptions on the specific extinction coefficient from Tegen et al. (1997), Table 1:
  REAL(KIND=ireals),    PARAMETER :: &
       beta_ext_dust = 1.5d3, &    ! m^2/kg
       beta_ext_org  = 8.0d3, &
       beta_ext_so4  = 8.0d3, &
       beta_ext_ss   = 0.2d3

  !   Estimated value of the aerosol number concentation,
  !    assuming a certain mean mass radius of the aerosols and a certain bulk density:
  REAL(KIND=ireals),    PARAMETER :: &
       rho_ss   = 3000.0_ireals, &    ! aerosol bulk density kg/m^3 for ss
       r_ss     = 500e-9_ireals, &    ! aerosol mean mass radius in m for ss  (200 - 12000 nm)
       rho_so4  = 2000.0_ireals, &    ! aerosol bulk density kg/m^3 for so4
       r_so4    = 80e-9_ireals, &     ! aerosol mean mass radius in m for so4 (50 - 100 nm)
       rho_org  = 2000.0_ireals, &    ! aerosol bulk density kg/m^3 for organics
       r_org    = 80e-9_ireals, &     ! aerosol mean mass radius in m for organics (50 - 100 nm)
       rho_dust = 3000.0_ireals, &    ! aerosol bulk density kg/m^3 for dust
       r_dust   = 1000e-9_ireals      ! aerosol mean mass radius in m for dust
  
  !   Soluble_fraction_dust: comparatively low number  ~0.0 - 0.1
  !   Soluble_fraction_organics: comparatively high number ~0.9 - 1.0
  REAL(KIND=ireals),    PARAMETER :: &
       soluble_fraction_dust     = 0.1_ireals, &
       soluble_fraction_organics = 0.9_ireals

  !   Parameters of assumed vertical exponential profile of qcn:
  !    ztrans = transition height in m, constant value below, expon. profile above
  !    z1oe = 1/e - height of expon. decrease above z0
  REAL(KIND=ireals),    PARAMETER :: &
       ztrans = 3000.0_ireals, &
       z1oe   = 2000.0_ireals, &
       zmax   = 12000.0_ireals

  !   Parameters of assumed vertical exponential profile of rho_air:
  !    z12 = 1/2 - height of expon. decrease of air density
  REAL(KIND=ireals),    PARAMETER :: &
       z12   = 6000.0_ireals, &
       z1oe_rho = z12 / 0.693147180559945286    ! = z12/LOG(2)

  !  Derived parameters for the vertical profile computations:
  REAL(KIND=ireals),    PARAMETER :: &
       z1oe_eff    = (z1oe*z1oe_rho) / (z1oe+z1oe_rho), &
       rho_air_msl = 1.225_ireals

  !-------------------------------------------------------------------------------------------

  ! In the parameterization of Tegen et al.,
  ! aer_org (organics) contains aer_bc (black carbon), so that
  ! only aer_org is used in the following.

  ! Total vertically integral aerosol number per m^2 from the optical thicknesses,
  !  assuming that the mean mass radius of each species is a constant everywhere:

  nscale(:) = 3.0 / (4.0*pi) * ( &
       aer_dust(:)/(beta_ext_dust*rho_dust*r_dust**3) * soluble_fraction_dust     + &
       aer_org(:)/(beta_ext_org*rho_org*r_org**3)     * soluble_fraction_organics + &
       aer_so4(:)/(beta_ext_so4*rho_so4*r_so4**3)     + &
       aer_ss(:)/(beta_ext_ss*rho_ss*r_ss**3)           &
       )
           
  ! From that, compute the value of the specific aerosol number (n/rho_air, unit 1/kg) in the surface layer,
  !  which serves as a scaling parameter in the vertical distribution:
  !  This value is computed under the assumptions that:
  !  - the vertical profile of the specific aerosol number is constant in a
  !    surface layer up to 3000 m MSL and decreases exponentially above,
  !  - the aerosols in the Tegen climatology are confined within the height layer from the ground
  !    to a zmax, which we assume to be 12 km,
  !  - the air density has the standard profile (bisection every 6000 m)

  z0(:) = MAX(ztrans, hhl(:,ke+1))
  rho_air_0(:) = rho_air_msl * EXP(-z0(:)/z1oe_rho)
  nscale(:) = nscale(:) / ( rho_air_0(:) * ( z1oe_rho*(EXP((z0(:)-hhl(:,ke+1))/z1oe_rho) - 1.0_ireals) + &
       z1oe_eff*(1.0_ireals - EXP((z0(:)-zmax)/z1oe_eff)) ) )

  ! Now compute ncn in 1/m^3 as function of height from nscale (1/kg) and the vertical profile of
  !  the specific aerosol number:
  DO k=kstart, kend
    DO i=istart,iend
      z = 0.5 * (hhl(i,k) + hhl(i,k+1))
      IF (z <= z0(i)) THEN
        ncn(i,k) = nscale(i) * rho_air_msl * EXP(-z/z1oe_rho)
      ELSE
        ncn(i,k) = nscale(i) * EXP((z0(i)-z)/z1oe-z/z1oe_rho) * rho_air_msl
      END IF
    END DO
  END DO

END SUBROUTINE ncn_from_tau_aerosol_speccnconst

!===========================================================================================

! Same as ncn_from_tau_aerosol_speccnconst but only for dust

SUBROUTINE ncn_from_tau_aerosol_speccnconst_dust ( hhl,hhl1, aer_dust, ncn_dust )
  !$ACC ROUTINE SEQ
  IMPLICIT NONE

  !-------------------------------------------------------------------------------------------
  !.. Inputs:

  REAL(KIND=ireals), INTENT(in)   :: hhl, hhl1
  REAL(KIND=ireals), INTENT(in)   :: aer_dust


  !-------------------------------------------------------------------------------------------
  !.. Outputs:

  REAL(KIND=ireals), INTENT(out)  :: ncn_dust

  !-------------------------------------------------------------------------------------------
  !.. Local variables:

  REAL(KIND=ireals)       :: z, z0, rho_air_0, nscale

  !-------------------------------------------------------------------------------------------
  !.. Local parameters:

  REAL(KIND=ireals), PARAMETER :: pi = 3.141592653589793

  !   Assumptions on the specific extinction coefficient:
  REAL(KIND=ireals),    PARAMETER :: &
       beta_ext_dust = 1.5e3_ireals      ! m^2/kg Tegen et al. (1997), Table 1

  !   Estimated value of the aerosol number concentation,
  !    assuming a certain mean mass radius of the aerosols and a certain bulk density:
  REAL(KIND=ireals),    PARAMETER :: &
       rho_dust = 2500.0_ireals, &    ! aerosol bulk density kg/m^3 for dust Linke et al. (2006)
       r_dust   = 1.0e-6_ireals       ! aerosol mean mass radius in m for dust Wong et al. (2021)
  
  !   Parameters of assumed vertical exponential profile of qcn:
  !    ztrans = transition height in m, constant value below, expon. profile above
  !    z1oe = 1/e - height of expon. decrease above z0
  REAL(KIND=ireals),    PARAMETER :: &
       ztrans = 3000.0_ireals, &
       z1oe   = 2000.0_ireals, &
       zmax   = 12000.0_ireals

  !   Parameters of assumed vertical exponential profile of rho_air:
  !    z12 = 1/2 - height of expon. decrease of air density
  REAL(KIND=ireals),    PARAMETER :: &
       z12   = 6000.0_ireals, &
       z1oe_rho = z12 / 0.693147180559945286    ! = z12/LOG(2)

  !  Derived parameters for the vertical profile computations:
  REAL(KIND=ireals),    PARAMETER :: &
       z1oe_eff    = (z1oe*z1oe_rho) / (z1oe+z1oe_rho), &
       rho_air_msl = 1.225_ireals

  !-------------------------------------------------------------------------------------------

  ! In the parameterization of Tegen et al.,
  ! aer_org (organics) contains aer_bc (black carbon), so that
  ! only aer_org is used in the following.

  ! Total vertically integral aerosol number per m^2 from the optical thicknesses,
  !  assuming that the mean mass radius of each species is a constant everywhere:

  nscale = 3.0 / (4.0*pi) * (aer_dust/(beta_ext_dust*rho_dust*r_dust**3) )
        
  ! From that, compute the value of the specific aerosol number (n/rho_air, unit 1/kg) in the surface layer,
  !  which serves as a scaling parameter in the vertical distribution:
  !  This value is computed under the assumptions that:
  !  - the vertical profile of the specific aerosol number is constant in a
  !    surface layer up to 3000 m MSL and decreases exponentially above,
  !  - the aerosols in the Tegen climatology are confined within the height layer from the ground
  !    to a zmax, which we assume to be 12 km,
  !  - the air density has the standard profile (bisection every 6000 m)

  z0 = MAX(ztrans, hhl)
  rho_air_0 = rho_air_msl * EXP(-z0/z1oe_rho)

  nscale = nscale / ( rho_air_0 * ( z1oe_rho*(EXP((z0-hhl)/z1oe_rho) - 1.0_ireals) + z1oe_eff*(1.0_ireals - EXP((z0-zmax)/z1oe_eff)) ) )
  ! Now compute ncn in 1/m^3 as function of height from nscale (1/kg) and the vertical profile of
  !  the specific aerosol number:

      z = 0.5 * (hhl + hhl1)
      IF (z <= z0) THEN
        ncn_dust = nscale * rho_air_msl * EXP(-z/z1oe_rho)
      ELSE
        ncn_dust = nscale * EXP((z0-z)/z1oe-z/z1oe_rho) * rho_air_msl
      END IF

      ! units are [1/m^3] BUT we want to convert to cm^-3 to use in demott formula so we multiply by 10^-6
      ncn_dust = MAX(MIN(1.0E-6_ireals*ncn_dust,1.0E5_ireals),1.0_ireals)

END SUBROUTINE ncn_from_tau_aerosol_speccnconst_dust

!===========================================================================================

! Calculate ice nucleation based on cooper (1987) parametrizations
! Input for cooper is only Temperature

SUBROUTINE ice_nucleation_cooper ( t, znin )
  !$ACC ROUTINE SEQ
  IMPLICIT NONE

  !-------------------------------------------------------------------------------------------
  !.. Inputs:

  REAL(KIND=ireals),       INTENT(in)           :: t             ! temperature (K)

  !-------------------------------------------------------------------------------------------
  !.. Outputs:

  REAL(KIND=ireals), INTENT(out)                :: znin          !< number of cloud ice crystals at nucleation

  !-------------------------------------------------------------------------------------------
  !.. Local parameters:

  REAL(KIND=ireals), PARAMETER :: t0 = 273.15_ireals
  REAL(KIND=ireals), PARAMETER :: a = 5.0E+0_ireals
  REAL(KIND=ireals), PARAMETER :: b = 0.304_ireals

  znin  = a * EXP(b * (t0 - t))

END SUBROUTINE ice_nucleation_cooper

!===========================================================================================

! Calculate ice nucleation based on demott (2015) parametrizations
! Input for demott is dust number concentration and Temperature

SUBROUTINE ice_nucleation_demott ( t, aerncn , znin )
  !$ACC ROUTINE SEQ
  IMPLICIT NONE

  !-------------------------------------------------------------------------------------------
  !.. Inputs:

  REAL(KIND=ireals),       INTENT(in)           :: t             ! temperature (K)
  REAL(KIND=ireals),       INTENT(in)           :: aerncn        ! dust number concentrations (cm^-3)

  !-------------------------------------------------------------------------------------------
  !.. Outputs:

  REAL(KIND=ireals), INTENT(out)                :: znin          !< number of cloud ice crystals at nucleation

  !-------------------------------------------------------------------------------------------
  !.. Local parameters:

  REAL(KIND=ireals), PARAMETER :: t0 = 273.15_ireals
  REAL(KIND=ireals), PARAMETER :: alpha = 1.0E3_ireals
  REAL(KIND=ireals), PARAMETER :: beta = 1.25_ireals
  REAL(KIND=ireals), PARAMETER :: gama = 0.46_ireals
  REAL(KIND=ireals), PARAMETER :: delta = -11.6_ireals

  znin = alpha * (aerncn**beta) * (EXP(gama * (t0 - t)+delta))  

  znin = MAX(MIN(znin,1.0E7_ireals),1.0E-7_ireals)

END SUBROUTINE ice_nucleation_demott

END MODULE mo_cpl_aerosol_microphys

