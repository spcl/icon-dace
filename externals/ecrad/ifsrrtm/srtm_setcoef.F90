! # 1 "ifsrrtm/srtm_setcoef.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_setcoef.f90"
! this file has been modified for the use in icon

subroutine srtm_setcoef &
 & ( kidia   , kfdia    , klev    ,&
 &   pavel   , ptavel   ,&
 &   pcoldry , pwkl     ,&
 &   klaytrop,&
 &   pcolch4  , pcolco2 , pcolh2o , pcolmol  , pcolo2 , pcolo3 ,&
 &   pforfac , pforfrac , kindfor , pselffac, pselffrac, kindself ,&
 &   pfac00  , pfac01   , pfac10  , pfac11  ,&
 &   kjp     , kjt      , kjt1    , prmu0    &
 & )  

!     j. delamere, aer, inc. (version 2.5, 02/04/01)

!     modifications:
!     jjmorcrette 030224   rewritten / adapted to ecmwf f90 system
!        m.hamrud      01-oct-2003 cy28 cleaning
!        d.salmond  31-oct-2007 vector version in the style of rrtm from meteo france & nec

!     purpose:  for a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.

use parkind1 , only : jpim, jprb
use ecradhook  , only : lhook, dr_hook
use yoesrtwn , only : preflog, tref
!!  use yoeswn  , only : ndbug

implicit none

!-- input arguments

integer(kind=jpim),intent(in)    :: kidia, kfdia
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(in)    :: pavel(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: ptavel(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcoldry(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pwkl(kidia:kfdia,35,klev) 
integer(kind=jpim),intent(inout) :: klaytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(inout) :: pcolch4(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pcolco2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pcolh2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pcolmol(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pcolo2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pcolo3(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pforfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pforfrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(inout) :: kindfor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pselffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pselffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(inout) :: kindself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pfac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pfac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pfac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: pfac11(kidia:kfdia,klev) 
integer(kind=jpim),intent(inout) :: kjp(kidia:kfdia,klev) 
integer(kind=jpim),intent(inout) :: kjt(kidia:kfdia,klev) 
integer(kind=jpim),intent(inout) :: kjt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: prmu0(kidia:kfdia) 
!-- output arguments

!-- local integers

integer(kind=jpim) :: i_nlayers, jk, jl, jp1

!-- local reals

real(kind=jprb) :: z_stpfac, z_plog
real(kind=jprb) :: z_fp, z_ft, z_ft1, z_water, z_scalefac
real(kind=jprb) :: z_factor, z_co2reg, z_compfp
!real(kind=jprb) :: z_tbndfrac, z_t0frac
real(kind=jprb) :: zhook_handle


associate(nflevg=>klev)
if (lhook) call dr_hook('srtm_setcoef',0,zhook_handle)

z_stpfac = 296._jprb/1013._jprb
i_nlayers = klev

!$acc parallel default(none) present(pavel, ptavel, pcoldry, pwkl, klaytrop, pcolch4, pcolco2, pcolh2o, pcolmol, &
!$acc   pcolo2,pcolo3, pforfac, pforfrac, kindfor, pselffac, pselffrac, kindself, pfac00, pfac01, pfac10, pfac11, &
!$acc   kjp, kjt,kjt1, prmu0) async(1)
!$acc loop seq
do jk = 1, klev
  !$acc loop gang(static:1) vector
  do jl = kidia, kfdia
    pcolmol(jl,jk) = 0.0_jprb
  enddo
enddo
!$acc loop gang(static:1) vector
do jl = kidia, kfdia
  if (prmu0(jl) > 0.0_jprb) then
    klaytrop(jl)  = 0
  endif
enddo

!$acc loop seq
do jk = 1, i_nlayers
  !$acc loop gang(static:1) vector &
  !$acc   private(z_plog, z_fp, z_ft, z_ft1, z_water, z_scalefac, z_factor, z_co2reg, z_compfp, jp1)
  do jl = kidia, kfdia
    if (prmu0(jl) > 0.0_jprb) then
      !        find the two reference pressures on either side of the
      !        layer pressure.  store them in jp and jp1.  store in fp the
      !        fraction of the difference (in ln(pressure)) between these
      !        two values that the layer pressure lies.

      z_plog = log(pavel(jl,jk))
      kjp(jl,jk) = int(36._jprb - 5._jprb*(z_plog+0.04_jprb))
      if (kjp(jl,jk) < 1) then
        kjp(jl,jk) = 1
      elseif (kjp(jl,jk) > 58) then
        kjp(jl,jk) = 58
      endif
      jp1 = kjp(jl,jk) + 1
      z_fp = 5. * (preflog(kjp(jl,jk)) - z_plog)

      !        determine, for each reference pressure (jp and jp1), which
      !        reference temperature (these are different for each  
      !        reference pressure) is nearest the layer temperature but does
      !        not exceed it.  store these indices in jt and jt1, resp.
      !        store in ft (resp. ft1) the fraction of the way between jt
      !        (jt1) and the next highest reference temperature that the 
      !        layer temperature falls.

      kjt(jl,jk) = int(3. + (ptavel(jl,jk)-tref(kjp(jl,jk)))/15.)
      if (kjt(jl,jk) < 1) then
        kjt(jl,jk) = 1
      elseif (kjt(jl,jk) > 4) then
        kjt(jl,jk) = 4
      endif
      z_ft = ((ptavel(jl,jk)-tref(kjp(jl,jk)))/15.) - real(kjt(jl,jk)-3)
      kjt1(jl,jk) = int(3. + (ptavel(jl,jk)-tref(jp1))/15.)
      if (kjt1(jl,jk) < 1) then
        kjt1(jl,jk) = 1
      elseif (kjt1(jl,jk) > 4) then
        kjt1(jl,jk) = 4
      endif
      z_ft1 = ((ptavel(jl,jk)-tref(jp1))/15.) - real(kjt1(jl,jk)-3)

      z_water = pwkl(jl,1,jk)/pcoldry(jl,jk)
      z_scalefac = pavel(jl,jk) * z_stpfac / ptavel(jl,jk)

      !        if the pressure is less than ~100mb, perform a different
      !        set of species interpolations.

      if (z_plog <= 4.56_jprb) go to 5300
      klaytrop(jl) =  klaytrop(jl) + 1

      !        set up factors needed to separately include the water vapor
      !        foreign-continuum in the calculation of absorption coefficient.

      pforfac(jl,jk) = z_scalefac / (1.+z_water)
      z_factor = (332.0-ptavel(jl,jk))/36.0
      kindfor(jl,jk) = min(2, max(1, int(z_factor)))
      pforfrac(jl,jk) = z_factor - real(kindfor(jl,jk))

      !        set up factors needed to separately include the water vapor
      !        self-continuum in the calculation of absorption coefficient.

      pselffac(jl,jk) = z_water * pforfac(jl,jk)
      z_factor = (ptavel(jl,jk)-188.0)/7.2
      kindself(jl,jk) = min(9, max(1, int(z_factor)-7))
      pselffrac(jl,jk) = z_factor - real(kindself(jl,jk) + 7)

      !        calculate needed column amounts.

      pcolh2o(jl,jk) = 1.e-20 * pwkl(jl,1,jk)
      pcolco2(jl,jk) = 1.e-20 * pwkl(jl,2,jk)
      pcolo3(jl,jk) = 1.e-20 * pwkl(jl,3,jk)
      !         colo3(lay) = 0.
      !         colo3(lay) = colo3(lay)/1.16
      pcolch4(jl,jk) = 1.e-20 * pwkl(jl,6,jk)
      pcolo2(jl,jk) = 1.e-20 * pwkl(jl,7,jk)
      pcolmol(jl,jk) = 1.e-20 * pcoldry(jl,jk) + pcolh2o(jl,jk)
      !         colco2(lay) = 0.
      !         colo3(lay) = 0.
      !         colch4(lay) = 0.
      !         colo2(lay) = 0.
      !         colmol(lay) = 0.
      if (pcolco2(jl,jk) == 0.) pcolco2(jl,jk) = 1.e-32 * pcoldry(jl,jk)
      if (pcolch4(jl,jk) == 0.) pcolch4(jl,jk) = 1.e-32 * pcoldry(jl,jk)
      if (pcolo2(jl,jk) == 0.) pcolo2(jl,jk) = 1.e-32 * pcoldry(jl,jk)
      !        using e = 1334.2 cm-1.
      z_co2reg = 3.55e-24 * pcoldry(jl,jk)
      go to 5400

      !        above laytrop.
5300  continue

      !        set up factors needed to separately include the water vapor
      !        foreign-continuum in the calculation of absorption coefficient.

      pforfac(jl,jk) = z_scalefac / (1.+z_water)
      z_factor = (ptavel(jl,jk)-188.0)/36.0
      kindfor(jl,jk) = 3
      pforfrac(jl,jk) = z_factor - 1.0

      !        calculate needed column amounts.

      pcolh2o(jl,jk) = 1.e-20 * pwkl(jl,1,jk)
      pcolco2(jl,jk) = 1.e-20 * pwkl(jl,2,jk)
      pcolo3(jl,jk)  = 1.e-20 * pwkl(jl,3,jk)
      pcolch4(jl,jk) = 1.e-20 * pwkl(jl,6,jk)
      pcolo2(jl,jk)  = 1.e-20 * pwkl(jl,7,jk)
      pcolmol(jl,jk) = 1.e-20 * pcoldry(jl,jk) + pcolh2o(jl,jk)
      if (pcolco2(jl,jk) == 0.) pcolco2(jl,jk) = 1.e-32 * pcoldry(jl,jk)
      if (pcolch4(jl,jk) == 0.) pcolch4(jl,jk) = 1.e-32 * pcoldry(jl,jk)
      if (pcolo2(jl,jk) == 0.) pcolo2(jl,jk)  = 1.e-32 * pcoldry(jl,jk)
      z_co2reg = 3.55e-24 * pcoldry(jl,jk)

      pselffac(jl,jk) =0.0_jprb
      pselffrac(jl,jk)=0.0_jprb
      kindself(jl,jk) = 0

5400  continue

      !        we have now isolated the layer ln pressure and temperature,
      !        between two reference pressures and two reference temperatures 
      !        (for each reference pressure).  we multiply the pressure 
      !        fraction fp with the appropriate temperature fractions to get 
      !        the factors that will be needed for the interpolation that yields
      !        the optical depths (performed in routines taugbn for band n).

      z_compfp = 1. - z_fp
      pfac10(jl,jk) = z_compfp * z_ft
      pfac00(jl,jk) = z_compfp * (1. - z_ft)
      pfac11(jl,jk) = z_fp * z_ft1
      pfac01(jl,jk) = z_fp * (1. - z_ft1)

      !  if (ndbug.le.3) then
      !    print 9000,lay,laytrop,jp(lay),jt(lay),jt1(lay),tavel(lay) &
      !      &,fac00(lay),fac01(lay),fac10(lay),fac11(lay) &
      !      &,colmol(lay),colch4(lay),colco2(lay),colh2o(lay) &
      !      &,colo2(lay),colo3(lay),selffac(lay),selffrac(lay) &
      !      &,forfac(lay),forfrac(lay),indself(lay),indfor(lay)
9000  format(1x,2i3,3i4,f6.1,4f7.2,12e9.2,2i5)
      !  endif

    endif
  enddo
enddo
!$acc end parallel

!----------------------------------------------------------------------- 
if (lhook) call dr_hook('srtm_setcoef',1,zhook_handle)
end associate
end subroutine srtm_setcoef
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

