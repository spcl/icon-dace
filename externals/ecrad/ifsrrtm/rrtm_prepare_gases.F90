! # 1 "ifsrrtm/rrtm_prepare_gases.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_prepare_gases.f90"
! this file has been modified for the use in icon

subroutine rrtm_prepare_gases &
 &( kidia, kfdia, klon, klev, &
 &  paph , pap , &
 &  pth  , pt  , &
 &  pq   , pco2 , pch4, pn2o  , pno2, pc11, pc12, pc22, pcl4, pozn, &
 &  pcoldry, pwbrodl, pwkl, pwx , &
 &  pavel  , ptavel , pz  , ptz , kreflect)  

!----compiled for cray with -h nopattern----

!     prepare the units of the gas concentrations for the longwave
!     rrtm gas absorption model.  this file is adapted from
!     rrtm_ecrt_140gp_mcica.f90, written mainly by jean-jacques
!     morcrette.

!- original
!     2015-07-15  robin hogan

!- modifications

use parkind1 , only : jpim, jprb
use ecradhook  , only : lhook, dr_hook, jphook
use yomcst   , only : rg
use parrrtm  , only : jpxsec, jpinpx  
use yomdyncore,only : rplrg

!------------------------------arguments--------------------------------

implicit none

integer(kind=jpim),intent(in)    :: klon! number of atmospheres (longitudes) 
integer(kind=jpim),intent(in)    :: klev! number of atmospheric layers 
integer(kind=jpim),intent(in)    :: kidia, kfdia 

real(kind=jprb)   ,intent(in)    :: paph(klon,klev+1) ! interface pressures (pa)
real(kind=jprb)   ,intent(in)    :: pap(klon,klev) ! layer pressures (pa)
real(kind=jprb)   ,intent(in)    :: pth(klon,klev+1) ! interface temperatures (k)
real(kind=jprb)   ,intent(in)    :: pt(klon,klev) ! layer temperature (k)
real(kind=jprb)   ,intent(in)    :: pq(klon,klev) ! h2o specific humidity (mmr)
real(kind=jprb)   ,intent(in)    :: pco2(klon,klev) ! co2 mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pch4(klon,klev) ! ch4 mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pn2o(klon,klev) ! n2o mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pno2(klon,klev) ! no2 mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pc11(klon,klev) ! cfc11 mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pc12(klon,klev) ! cfc12 mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pc22(klon,klev) ! cfc22 mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pcl4(klon,klev) ! ccl4  mass mixing ratio
real(kind=jprb)   ,intent(in)    :: pozn(klon,klev) ! o3 mass mixing ratio

real(kind=jprb)   ,intent(out)   :: pcoldry(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: pwbrodl(kidia:kfdia,klev) ! broadening gas column density (mol/cm2)
real(kind=jprb)   ,intent(out)   :: pwkl(kidia:kfdia,jpinpx,klev) 
real(kind=jprb)   ,intent(out)   :: pwx(kidia:kfdia,jpxsec,klev) ! amount of trace gases
real(kind=jprb)   ,intent(out)   :: pavel(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: ptavel(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: pz(kidia:kfdia,0:klev) 
real(kind=jprb)   ,intent(out)   :: ptz(kidia:kfdia,0:klev) 
integer(kind=jpim),intent(out)   :: kreflect(kidia:kfdia) 

!      real rch4                       ! ch4 mass mixing ratio
!      real rn2o                       ! n2o mass mixing ratio
!      real rcfc11                     ! cfc11 mass mixing ratio
!      real rcfc12                     ! cfc12 mass mixing ratio
!      real rcfc22                     ! cfc22 mass mixing ratio
!      real rccl4                      ! ccl4  mass mixing ratio
!- from profile             
!- from surface             
real(kind=jprb) :: zamd                  ! effective molecular weight of dry air (g/mol)
real(kind=jprb) :: zamw                  ! molecular weight of water vapor (g/mol)
real(kind=jprb) :: zamco2                ! molecular weight of carbon dioxide (g/mol)
real(kind=jprb) :: zamo                  ! molecular weight of ozone (g/mol)
real(kind=jprb) :: zamch4                ! molecular weight of methane (g/mol)
real(kind=jprb) :: zamn2o                ! molecular weight of nitrous oxide (g/mol)
real(kind=jprb) :: zamc11                ! molecular weight of cfc11 (g/mol) - cfcl3
real(kind=jprb) :: zamc12                ! molecular weight of cfc12 (g/mol) - cf2cl2
real(kind=jprb) :: zamc22                ! molecular weight of cfc22 (g/mol) - chf2cl
real(kind=jprb) :: zamcl4                ! molecular weight of ccl4  (g/mol) - ccl4
real(kind=jprb) :: zavgdro               ! avogadro's number (molecules/mole)
real(kind=jprb) :: zgravit               ! gravitational acceleration (cm/s**2)

real(kind=jprb) :: zsummol

! atomic weights for conversion from mass to volume mixing ratios; these
!  are the same values used in ecrt to assure accurate conversion to vmr
data zamd   /  28.970_jprb    /
data zamw   /  18.0154_jprb   /
data zamco2 /  44.011_jprb    /
data zamo   /  47.9982_jprb   /
data zamch4 /  16.043_jprb    /
data zamn2o /  44.013_jprb    /
data zamc11 / 137.3686_jprb   /
data zamc12 / 120.9140_jprb   /
data zamc22 /  86.4690_jprb   /
data zamcl4 / 153.8230_jprb   /
data zavgdro/ 6.02214e23_jprb /

integer(kind=jpim) :: iatm, jmol, ixmax, j1, j2, jk, jl
integer(kind=jpim), parameter :: itmol = 7

real(kind=jprb) :: zamm

real(kind=jphook) :: zhook_handle

! ***

! *** mji
! initialize all molecular amounts to zero here, 
! then pass ecrt amounts into rrtm arrays below.

!      data zwkl /maxprdw*0.0/
!      data zwx  /maxprod*0.0/
!      data kreflect /0/

! activate cross section molecules:
!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
!      data kxmol  /2/
!      data kxindx /0,2,3,0,31*0/

!      ireflect=kreflect
!      nxmol=kxmol

if (lhook) call dr_hook('rrtm_prepare_gases',0,zhook_handle)

!$acc data present(paph, pap, &
!$acc              pth, pt, &
!$acc              pq, pco2, pch4, pn2o, pno2, pc11, pc12, pc22, pcl4, pozn, &
!$acc              pcoldry, pwbrodl, pwkl, pwx , &
!$acc              pavel, ptavel, pz, ptz , kreflect)

zgravit=(rg/rplrg)*1.e2_jprb

!$acc parallel default(none) async(1)
!$acc loop gang vector
do jl = kidia, kfdia
  kreflect(jl)=0
enddo
!$acc end parallel

!$acc parallel default(none) async(1)
!do j1=1,35
! ixindx(j1)=0
!$acc loop gang vector collapse(3)
do j2=1,klev
  do j1=1,35
    do jl = kidia, kfdia
      pwkl(jl,j1,j2)=0.0_jprb 
    enddo
  enddo
enddo
!ixindx(2)=2
!ixindx(3)=3
!$acc end parallel

!     set parameters needed for rrtm execution:
iatm    = 0
!      ixsect  = 1
!      numangs = 0
!      iout    = -1
ixmax   = 4

!$acc parallel default(none) async(1)
!$acc loop gang(static:1) vector
do jl = kidia, kfdia
!     install ecrt arrays into rrtm arrays for pressure, temperature,
!     and molecular amounts.  pressures are converted from pascals
!     (ecrt) to mb (rrtm).  h2o, co2, o3 and trace gas amounts are 
!     converted from mass mixing ratio to volume mixing ratio.  co2
!     converted with same dry air and co2 molecular weights used in 
!     ecrt to assure correct conversion back to the proper co2 vmr.
!     the dry air column coldry (in molec/cm2) is calculated from 
!     the level pressures pz (in mb) based on the hydrostatic equation
!     and includes a correction to account for h2o in the layer.  the
!     molecular weight of moist air (amm) is calculated for each layer.
!     note: rrtm levels count from bottom to top, while the ecrt input
!     variables count from the top down and must be reversed 
  pz(jl,0) = paph(jl,klev+1)/100._jprb
  ptz(jl,0) = pth(jl,klev+1)
enddo

  !$acc loop seq
  do jk = 1, klev
    !$acc loop gang(static:1) vector private(zamm)
    do jl = kidia, kfdia
    pavel(jl,jk) = pap(jl,klev-jk+1)/100._jprb
    ptavel(jl,jk) = pt(jl,klev-jk+1)
    pz(jl,jk) = paph(jl,klev-jk+1)/100._jprb
    ptz(jl,jk) = pth(jl,klev-jk+1)
    ! rrtmg cannot cope with zero or negative water vapour
    pwkl(jl,1,jk) = max(pq(jl,klev-jk+1),1.0e-15)*zamd/zamw
    pwkl(jl,2,jk) = pco2(jl,klev-jk+1)*zamd/zamco2
    pwkl(jl,3,jk) = pozn(jl,klev-jk+1)*zamd/zamo
    pwkl(jl,4,jk) = pn2o(jl,klev-jk+1)*zamd/zamn2o
    pwkl(jl,6,jk) = pch4(jl,klev-jk+1)*zamd/zamch4
    pwkl(jl,7,jk) = 0.209488_jprb
    zamm = (1.0_jprb-pwkl(jl,1,jk))*zamd + pwkl(jl,1,jk)*zamw
    pcoldry(jl,jk) = (pz(jl,jk-1)-pz(jl,jk))*1.e3_jprb*zavgdro/(zgravit*zamm*(1.0_jprb+pwkl(jl,1,jk)))
  enddo
  enddo
  !$acc end parallel

  !$acc parallel default(none) async(1)
  !$acc loop gang vector collapse(3)
  do j2=1,klev
    do j1=1,jpxsec
      do jl = kidia, kfdia
        pwx(jl,j1,j2)=0.0_jprb
      enddo
    enddo
  enddo
  !$acc end parallel


  !$acc parallel default(none) async(1)
  !$acc loop seq
  do jk = 1, klev
!$acc loop gang vector private (zsummol) 
do jl = kidia, kfdia
!- set cross section molecule amounts from ecrt; convert to vmr
    pwx(jl,1,jk) = pcl4(jl,klev-jk+1) * zamd/zamcl4
    pwx(jl,2,jk) = pc11(jl,klev-jk+1) * zamd/zamc11
    pwx(jl,3,jk) = pc12(jl,klev-jk+1) * zamd/zamc12
    pwx(jl,4,jk) = pc22(jl,klev-jk+1) * zamd/zamc22
    pwx(jl,1,jk) = pcoldry(jl,jk) * pwx(jl,1,jk) * 1.e-20_jprb
    pwx(jl,2,jk) = pcoldry(jl,jk) * pwx(jl,2,jk) * 1.e-20_jprb
    pwx(jl,3,jk) = pcoldry(jl,jk) * pwx(jl,3,jk) * 1.e-20_jprb
    pwx(jl,4,jk) = pcoldry(jl,jk) * pwx(jl,4,jk) * 1.e-20_jprb

!- here, all molecules in wkl and wx are in volume mixing ratio; convert to
!  molec/cm2 based on coldry for use in rrtm

!cdir unroll=6
zsummol = 0.0_jprb
!ab broadening gases
    !$acc loop seq
    do jmol = 2, itmol
      zsummol = zsummol + pwkl(jl,jmol,jk)
    enddo
    pwbrodl(jl,jk) = pcoldry(jl,jk) * (1._jprb - zsummol)
    !$acc loop seq
    do jmol = 1, itmol
      pwkl(jl,jmol,jk) = pcoldry(jl,jk) * pwkl(jl,jmol,jk)
    enddo    
  enddo
enddo
!$acc end parallel
!$acc wait
!$acc end data

!     ------------------------------------------------------------------
if (lhook) call dr_hook('rrtm_prepare_gases',1,zhook_handle)

end subroutine rrtm_prepare_gases
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

