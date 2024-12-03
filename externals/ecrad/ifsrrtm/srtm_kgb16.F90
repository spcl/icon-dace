! # 1 "ifsrrtm/srtm_kgb16.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_kgb16.f90"
! this file has been modified for the use in icon

subroutine srtm_kgb16(directory)

!     originally by j.delamere, atmospheric & environmental research.
!     revision: 2.4
!     band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!     reformatted for f90 by jjmorcrette, ecmwf
!     r. elkhatib 12-10-2005 split for faster and more robust compilation.
!     g.mozdzynski march 2011 read constants from files
!     t. wilhelmsson and k. yessad (oct 2013) geometry and setup refactoring.
!     ------------------------------------------------------------------

use parkind1  , only : jprb
use ecradhook   , only : lhook, dr_hook
use yomlun    , only : nulrad
use yommp0    , only : nproc, myproc
use mpl_module, only : mpl_broadcast
use yomtag    , only : mtagrad
use yoesrta16 , only : ka, kb, ka_d, kb_d, selfref, forref, sfluxref, rayl, strrat1, layreffr

!     ------------------------------------------------------------------

implicit none

character(len=*), intent(in) :: directory

! kurucz
!character(len = 80) :: clzzz
character(len = 255) :: clf1
real(kind=jprb) :: zhook_handle


! # 1 "./include/abor1.intfb.h" 1
interface
subroutine abor1(cdtext)
character(len=*), intent(in) :: cdtext
end subroutine abor1
end interface
! # 34 "ifsrrtm/srtm_kgb16.f90" 2

if (lhook) call dr_hook('srtm_kgb16',0,zhook_handle)

if( myproc==1 )then
  !call getenv("data",clzzz)
  !if(clzzz /= " ") then
  !  clf1=trim(clzzz)//"/radsrtm"
  clf1 = directory // "/radsrtm"

  open(nulrad,file=trim(clf1),convert="big_endian",form="unformatted",action="read",err=1000)




  !else
  !  open(nulrad,file='radsrtm',form="unformatted",action="read",err=1000)
  !endif
 read(nulrad,err=1001) ka_d,kb_d
  ka = real(ka_d,jprb)
  kb = real(kb_d,jprb) 
endif
if( nproc>1 )then
  call mpl_broadcast (ka,mtagrad,1,cdstring='srtm_kgb16:')
  call mpl_broadcast (kb,mtagrad,1,cdstring='srtm_kgb16:')
endif

sfluxref = (/ &
 & 1.92269_jprb    , 1.72844_jprb    , 1.64326_jprb    , 1.58451_jprb     &
 & , 1.44031_jprb    , 1.25108_jprb    , 1.02724_jprb    , 0.776759_jprb    &
 & , 0.534444_jprb   , 5.87755e-02_jprb, 4.86706e-02_jprb, 3.87989e-02_jprb &
 & , 2.84532e-02_jprb, 1.82431e-02_jprb, 6.92320e-03_jprb, 9.70770e-04_jprb /)  

!     rayleigh extinction coefficient at v = 2925 cm-1.
rayl = 2.91e-10_jprb

strrat1 = 252.131_jprb

layreffr = 18

!     ------------------------------------------------------------------

!     the array ka contains absorption coefs at the 16 chosen g-values 
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
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
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
!     -----------------------------------------------------------------

forref(:, 1) = (/ 0.525585e-05_jprb, 0.527618e-05_jprb, 0.746929e-04_jprb /)
forref(:, 2) = (/ 0.794660e-05_jprb, 0.136902e-04_jprb, 0.849878e-04_jprb /)
forref(:, 3) = (/ 0.197099e-04_jprb, 0.733094e-04_jprb, 0.121687e-03_jprb /)
forref(:, 4) = (/ 0.148274e-03_jprb, 0.169776e-03_jprb, 0.164848e-03_jprb /)
forref(:, 5) = (/ 0.230296e-03_jprb, 0.210384e-03_jprb, 0.182028e-03_jprb /)
forref(:, 6) = (/ 0.280575e-03_jprb, 0.259217e-03_jprb, 0.196080e-03_jprb /)
forref(:, 7) = (/ 0.329034e-03_jprb, 0.291575e-03_jprb, 0.207044e-03_jprb /)
forref(:, 8) = (/ 0.349989e-03_jprb, 0.323471e-03_jprb, 0.225712e-03_jprb /)
forref(:, 9) = (/ 0.366097e-03_jprb, 0.321519e-03_jprb, 0.253150e-03_jprb /)
forref(:,10) = (/ 0.383589e-03_jprb, 0.355314e-03_jprb, 0.262555e-03_jprb /)
forref(:,11) = (/ 0.375933e-03_jprb, 0.372443e-03_jprb, 0.261313e-03_jprb /)
forref(:,12) = (/ 0.370652e-03_jprb, 0.382366e-03_jprb, 0.250070e-03_jprb /)
forref(:,13) = (/ 0.375092e-03_jprb, 0.379542e-03_jprb, 0.265794e-03_jprb /)
forref(:,14) = (/ 0.389705e-03_jprb, 0.384274e-03_jprb, 0.322135e-03_jprb /)
forref(:,15) = (/ 0.372084e-03_jprb, 0.390422e-03_jprb, 0.370035e-03_jprb /)
forref(:,16) = (/ 0.437802e-03_jprb, 0.373406e-03_jprb, 0.373222e-03_jprb /)

!     -----------------------------------------------------------------
!     the array selfref contains the coefficient of the water vapor
!     self-continuum (including the energy term).  the first index
!     refers to temperature in 7.2 degree increments.  for instance,
!     jt = 1 refers to a temperature of 245.6, jt = 2 refers to 252.8,
!     etc.  the second index runs over the g-channel (1 to 16).

selfref(:, 1) = (/ &
 & 0.126758e-02_jprb, 0.105253e-02_jprb, 0.873963e-03_jprb, 0.725690e-03_jprb, 0.602573e-03_jprb, &
 & 0.500344e-03_jprb, 0.415458e-03_jprb, 0.344973e-03_jprb, 0.286447e-03_jprb, 0.237849e-03_jprb /)  
selfref(:, 2) = (/ &
 & 0.144006e-02_jprb, 0.118514e-02_jprb, 0.975351e-03_jprb, 0.802697e-03_jprb, 0.660606e-03_jprb, &
 & 0.543667e-03_jprb, 0.447429e-03_jprb, 0.368226e-03_jprb, 0.303044e-03_jprb, 0.249400e-03_jprb /)  
selfref(:, 3) = (/ &
 & 0.294018e-02_jprb, 0.227428e-02_jprb, 0.175920e-02_jprb, 0.136077e-02_jprb, 0.105258e-02_jprb, &
 & 0.814189e-03_jprb, 0.629789e-03_jprb, 0.487153e-03_jprb, 0.376821e-03_jprb, 0.291478e-03_jprb /)  
selfref(:, 4) = (/ &
 & 0.395290e-02_jprb, 0.348405e-02_jprb, 0.307081e-02_jprb, 0.270658e-02_jprb, 0.238556e-02_jprb, &
 & 0.210261e-02_jprb, 0.185322e-02_jprb, 0.163341e-02_jprb, 0.143967e-02_jprb, 0.126891e-02_jprb /)  
selfref(:, 5) = (/ &
 & 0.419122e-02_jprb, 0.385638e-02_jprb, 0.354829e-02_jprb, 0.326481e-02_jprb, 0.300398e-02_jprb, &
 & 0.276399e-02_jprb, 0.254317e-02_jprb, 0.234000e-02_jprb, 0.215305e-02_jprb, 0.198104e-02_jprb /)  
selfref(:, 6) = (/ &
 & 0.495659e-02_jprb, 0.456777e-02_jprb, 0.420945e-02_jprb, 0.387924e-02_jprb, 0.357494e-02_jprb, &
 & 0.329450e-02_jprb, 0.303606e-02_jprb, 0.279790e-02_jprb, 0.257842e-02_jprb, 0.237615e-02_jprb /)  
selfref(:, 7) = (/ &
 & 0.526981e-02_jprb, 0.490687e-02_jprb, 0.456893e-02_jprb, 0.425426e-02_jprb, 0.396126e-02_jprb, &
 & 0.368844e-02_jprb, 0.343441e-02_jprb, 0.319788e-02_jprb, 0.297764e-02_jprb, 0.277256e-02_jprb /)  
selfref(:, 8) = (/ &
 & 0.575426e-02_jprb, 0.531597e-02_jprb, 0.491106e-02_jprb, 0.453699e-02_jprb, 0.419141e-02_jprb, &
 & 0.387216e-02_jprb, 0.357722e-02_jprb, 0.330475e-02_jprb, 0.305303e-02_jprb, 0.282048e-02_jprb /)  
selfref(:, 9) = (/ &
 & 0.549881e-02_jprb, 0.514328e-02_jprb, 0.481074e-02_jprb, 0.449970e-02_jprb, 0.420877e-02_jprb, &
 & 0.393665e-02_jprb, 0.368213e-02_jprb, 0.344406e-02_jprb, 0.322138e-02_jprb, 0.301310e-02_jprb /)  
selfref(:,10) = (/ &
 & 0.605357e-02_jprb, 0.561246e-02_jprb, 0.520349e-02_jprb, 0.482432e-02_jprb, 0.447278e-02_jprb, &
 & 0.414686e-02_jprb, 0.384469e-02_jprb, 0.356453e-02_jprb, 0.330479e-02_jprb, 0.306398e-02_jprb /)  
selfref(:,11) = (/ &
 & 0.640504e-02_jprb, 0.587858e-02_jprb, 0.539540e-02_jprb, 0.495194e-02_jprb, 0.454492e-02_jprb, &
 & 0.417136e-02_jprb, 0.382850e-02_jprb, 0.351382e-02_jprb, 0.322501e-02_jprb, 0.295993e-02_jprb /)  
selfref(:,12) = (/ &
 & 0.677803e-02_jprb, 0.615625e-02_jprb, 0.559152e-02_jprb, 0.507859e-02_jprb, 0.461271e-02_jprb, &
 & 0.418957e-02_jprb, 0.380524e-02_jprb, 0.345617e-02_jprb, 0.313913e-02_jprb, 0.285116e-02_jprb /)  
selfref(:,13) = (/ &
 & 0.690347e-02_jprb, 0.627003e-02_jprb, 0.569472e-02_jprb, 0.517219e-02_jprb, 0.469761e-02_jprb, &
 & 0.426658e-02_jprb, 0.387509e-02_jprb, 0.351953e-02_jprb, 0.319659e-02_jprb, 0.290328e-02_jprb /)  
selfref(:,14) = (/ &
 & 0.692680e-02_jprb, 0.632795e-02_jprb, 0.578087e-02_jprb, 0.528109e-02_jprb, 0.482452e-02_jprb, &
 & 0.440742e-02_jprb, 0.402638e-02_jprb, 0.367828e-02_jprb, 0.336028e-02_jprb, 0.306977e-02_jprb /)  
selfref(:,15) = (/ &
 & 0.754894e-02_jprb, 0.681481e-02_jprb, 0.615207e-02_jprb, 0.555378e-02_jprb, 0.501367e-02_jprb, &
 & 0.452609e-02_jprb, 0.408593e-02_jprb, 0.368857e-02_jprb, 0.332986e-02_jprb, 0.300603e-02_jprb /)  
selfref(:,16) = (/ &
 & 0.760689e-02_jprb, 0.709755e-02_jprb, 0.662232e-02_jprb, 0.617891e-02_jprb, 0.576519e-02_jprb, &
 & 0.537917e-02_jprb, 0.501899e-02_jprb, 0.468293e-02_jprb, 0.436938e-02_jprb, 0.407682e-02_jprb /)  

if (lhook) call dr_hook('srtm_kgb16',1,zhook_handle)
return

1000 continue
call abor1("srtm_kgb16:error opening file radsrtm")
1001 continue
call abor1("srtm_kgb16:error reading file radsrtm")

end subroutine srtm_kgb16
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

