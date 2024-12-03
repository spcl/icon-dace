! # 1 "ifsrrtm/srtm_gas_optical_depth.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_gas_optical_depth.f90"
! this file has been modified for the use in icon




subroutine srtm_gas_optical_depth &
 & ( kidia   , kfdia   , klev    , poneminus, &
 &   prmu0, &
 &   klaytrop,&
 &   pcolch4  , pcolco2 , pcolh2o , pcolmol  , pcolo2 , pcolo3 ,&
 &   pforfac , pforfrac , kindfor , pselffac, pselffrac, kindself ,&
 &   pfac00  , pfac01   , pfac10  , pfac11 ,&
 &   kjp     , kjt      , kjt1 ,&
 !-- output arrays 
 &   pod, pssa, pincsol)


!**** *srtm_gas_optical_depth* - spectral loop to compute the shortwave radiation fluxes.

!     purpose.
!     --------

!          compute the gas optical depth at each shortwave g point

!**   interface.
!     ----------

!          *srtm_gas_optical_depth* is called from the new radiation scheme

!        implicit arguments :
!        --------------------

!     ==== inputs ===
!     ==== outputs ===

!     method.
!     -------

!     externals.
!     ----------

!     reference.
!     ----------

!        see radiation's part of the ecmwf research department
!        documentation
!     author.
!     -------
!        adapted from srtm_spcvrt_mcica (by jean-jacques morcrette) by
!        robin hogan
!
!     modifications.
!     --------------
!        original : 2015-07-16

!     ------------------------------------------------------------------

use parkind1 , only : jpim, jprb
use ecradhook  , only : lhook, dr_hook, jphook
use parsrtm  , only : jpb1, jpb2
use yoesrtm  , only : jpgpt
use yoesrtwn , only : ngc

implicit none

!     ------------------------------------------------------------------

!*       0.1   arguments
!              ---------

integer(kind=jpim),intent(in)    :: kidia, kfdia
integer(kind=jpim),intent(in)    :: klev
real(kind=jprb)   ,intent(in)    :: poneminus(kidia:kfdia)
real(kind=jprb)   ,intent(in)    :: prmu0(kidia:kfdia)
integer(kind=jpim),intent(in)    :: klaytrop(kidia:kfdia)
real(kind=jprb)   ,intent(in)    :: pcolch4(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pcolco2(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pcolh2o(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pcolmol(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pcolo2(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pcolo3(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pforfac(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pforfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in)    :: kindfor(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pselffac(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pselffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in)    :: kindself(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pfac00(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pfac01(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pfac10(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: pfac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in)    :: kjp(kidia:kfdia,klev)
integer(kind=jpim),intent(in)    :: kjt(kidia:kfdia,klev)
integer(kind=jpim),intent(in)    :: kjt1(kidia:kfdia,klev)

real(kind=jprb)   ,intent(inout) :: pod(kidia:kfdia,klev,jpgpt) ! optical depth
real(kind=jprb)   ,intent(inout) :: pssa(kidia:kfdia,klev,jpgpt) ! single scattering albedo
real(kind=jprb)   ,intent(inout) :: pincsol(kidia:kfdia,jpgpt) ! incoming solar flux


!     ------------------------------------------------------------------

integer(kind=jpim) :: ib1, ib2, ibm, igt, iw(kidia:kfdia), jb, jg, jk, jl, icount

!-- output of srtm_taumoln routines
real(kind=jprb) :: ztaug(kidia:kfdia,klev,16) ! absorption optical depth
real(kind=jprb) :: ztaur(kidia:kfdia,klev,16) ! rayleigh optical depth
real(kind=jprb) :: zsflxzen(kidia:kfdia,16) ! incoming solar flux

real(kind=jphook) :: zhook_handle



! # 1 "./include/srtm_taumol16.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol16&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colch4 , p_colmol,&
 & k_laytrop , p_selffac , p_selffrac, k_indself , p_forfac , p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol16
end interface
! # 114 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol17.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol17&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colco2 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol17
end interface
! # 115 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol18.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol18&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colch4 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol18
end interface
! # 116 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol19.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol19&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colco2 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol19
end interface
! # 117 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol20.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol20&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1,&
 & p_colh2o , p_colch4 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol20
end interface
! # 118 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol21.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol21&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colco2 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol21
end interface
! # 119 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol22.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol22&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colmol , p_colo2,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo2(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol22
end interface
! # 120 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol23.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol23&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1,&
 & p_colh2o , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol23
end interface
! # 121 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol24.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol24&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colh2o , p_colmol , p_colo2 , p_colo3,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol24
end interface
! # 122 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol25.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol25&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1,&
 & p_colh2o , p_colmol , p_colo3,&
 & k_laytrop,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol25
end interface
! # 123 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol26.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol26&
 & ( kidia , kfdia , klev,&
 & p_colmol ,k_laytrop,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol26
end interface
! # 124 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol27.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol27&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1,&
 & p_colmol , p_colo3,&
 & k_laytrop,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol27
end interface
! # 125 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol28.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol28&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1 , p_oneminus,&
 & p_colmol , p_colo2 , p_colo3,&
 & k_laytrop,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol28
end interface
! # 126 "ifsrrtm/srtm_gas_optical_depth.f90" 2

! # 1 "./include/srtm_taumol29.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_taumol29&
 & ( kidia , kfdia , klev,&
 & p_fac00 , p_fac01 , p_fac10 , p_fac11,&
 & k_jp , k_jt , k_jt1,&
 & p_colh2o , p_colco2 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug , p_taur , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
use parsrtm , only : jpg
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colmol(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg)
real(kind=jprb) ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_taumol29
end interface
! # 127 "ifsrrtm/srtm_gas_optical_depth.f90" 2

!     ------------------------------------------------------------------

if (lhook) call dr_hook('srtm_gas_optical_depth',0,zhook_handle)

!$acc data create(iw, ztaug, ztaur, zsflxzen) &
!$acc     present(poneminus, prmu0, klaytrop, pcolch4, pcolco2, pcolh2o, &
!$acc             pcolmol, pcolo2, pcolo3, pforfac, pforfrac, kindfor, pselffac, &
!$acc             pselffrac, kindself, pfac00, pfac01, pfac10, pfac11, kjp, &
!$acc             kjt, kjt1, pod, pssa, pincsol)

ib1=jpb1
ib2=jpb2

icount=0
!$acc wait
!$acc parallel loop gang vector default(none) reduction(+:icount)
do jl = kidia, kfdia
  if (prmu0(jl) > 0.0_jprb) then
    icount=icount+1
    iw(jl)=0
  endif
enddo
!$acc end parallel loop

if (icount/=0) then

  do jb = ib1, ib2
    ibm = jb-15
    igt = ngc(ibm)

    !-- for each band, computes the gaseous and rayleigh optical thickness 
    !  for all g-points within the band

    if (jb == 16) then
      call srtm_taumol16 &
      & ( kidia   , kfdia    , klev    ,&
      &   pfac00  , pfac01   , pfac10   , pfac11   ,&
      &   kjp     , kjt      , kjt1     , poneminus,&
      &   pcolh2o , pcolch4  , pcolmol  ,&
      &   klaytrop, pselffac , pselffrac, kindself, pforfac  , pforfrac, kindfor ,&
      &   zsflxzen, ztaug    , ztaur    , prmu0     &
      & )  

    elseif (jb == 17) then
      call srtm_taumol17 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     , poneminus ,&
      &   pcolh2o , pcolco2 , pcolmol  ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )

    elseif (jb == 18) then
      call srtm_taumol18 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     , poneminus ,&
      &   pcolh2o , pcolch4 , pcolmol  ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 19) then
      call srtm_taumol19 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     , poneminus ,&
      &   pcolh2o , pcolco2 , pcolmol  ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 20) then
      call srtm_taumol20 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     ,&
      &   pcolh2o , pcolch4 , pcolmol  ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 21) then
      call srtm_taumol21 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     , poneminus ,&
      &   pcolh2o , pcolco2 , pcolmol  ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 22) then
      call srtm_taumol22 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     , poneminus ,&
      &   pcolh2o , pcolmol , pcolo2   ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 23) then
      call srtm_taumol23 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     ,&
      &   pcolh2o , pcolmol ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 24) then
      call srtm_taumol24 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     , poneminus ,&
      &   pcolh2o , pcolmol , pcolo2   , pcolo3 ,&
      &   klaytrop, pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 25) then
      !--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
      call srtm_taumol25 &
      & ( kidia    , kfdia   , klev     ,&
      &   pfac00   , pfac01  , pfac10 , pfac11 ,&
      &   kjp      , kjt     , kjt1   ,&
      &   pcolh2o  , pcolmol , pcolo3 ,&
      &   klaytrop ,&
      &   zsflxzen, ztaug   , ztaur   , prmu0     &
      & )  

    elseif (jb == 26) then
      !--- uv-a 22650-29000 cm-1   0.3448 - 0.4415 um
      call srtm_taumol26 &
      & ( kidia   , kfdia   , klev    ,&
      &   pcolmol ,klaytrop,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 27) then
      !--- uv-b 29000-38000 cm-1   0.2632 - 0.3448 um
      call srtm_taumol27 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10   , pfac11 ,&
      &   kjp     , kjt     , kjt1     ,&
      &   pcolmol , pcolo3 ,&
      &   klaytrop ,&
      &   zsflxzen, ztaug   , ztaur    , prmu0     &
      & )  

    elseif (jb == 28) then
      !--- uv-c 38000-50000 cm-1   0.2000 - 0.2632 um
      call srtm_taumol28 &
      & ( kidia   , kfdia   , klev    ,&
      &   pfac00  , pfac01  , pfac10 , pfac11 ,&
      &   kjp     , kjt     , kjt1   , poneminus ,&
      &   pcolmol , pcolo2  , pcolo3 ,&
      &   klaytrop ,&
      &   zsflxzen, ztaug   , ztaur  , prmu0     &
      & )  

    elseif (jb == 29) then
      call srtm_taumol29 &
      & ( kidia    , kfdia   , klev     ,&
      &   pfac00   , pfac01  , pfac10   , pfac11 ,&
      &   kjp      , kjt     , kjt1     ,&
      &   pcolh2o  , pcolco2 , pcolmol  ,&
      &   klaytrop , pselffac, pselffrac, kindself  , pforfac, pforfrac, kindfor ,&
      &   zsflxzen , ztaug   , ztaur    , prmu0     &
      & )  

    endif
    
    !$acc parallel default(none) async(1)
    !$acc loop seq
    do jg=1,igt
! added for dwd (2020)
!nec$ ivdep
      !$acc loop gang(static:1) vector
      do jl = kidia, kfdia
        if (prmu0(jl) > 0.0_jprb) then
          iw(jl)=iw(jl)+1
          ! incoming solar flux into plane perp to incoming radiation
          pincsol(jl,iw(jl)) = zsflxzen(jl,jg)
        endif
      enddo

      !$acc loop seq
      do jk=1,klev
        !$acc loop gang(static:1) vector private(jl)
        do jl = kidia, kfdia
          if (prmu0(jl) > 0.0_jprb) then
            pod (jl,jk,iw(jl)) = ztaur(jl,jk,jg) + ztaug(jl,jk,jg)
            pssa(jl,jk,iw(jl)) = ztaur(jl,jk,jg) / pod(jl,jk,iw(jl))
          endif
        enddo
      enddo

    enddo   !-- end loop on jg (g point)
    !$acc end parallel

  enddo     !-- end loop on jb (band)

endif

!$acc wait
!$acc end data

!     ------------------------------------------------------------------

if (lhook) call dr_hook('srtm_gas_optical_depth',1,zhook_handle)

end subroutine srtm_gas_optical_depth
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

