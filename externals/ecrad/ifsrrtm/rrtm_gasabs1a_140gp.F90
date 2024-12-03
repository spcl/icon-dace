! # 1 "ifsrrtm/rrtm_gasabs1a_140gp.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_gasabs1a_140gp.f90"
!option! -pvctl no_on_adb
!option! -pvctl nocollapse
subroutine rrtm_gasabs1a_140gp (kidia,kfdia,klev,patr1,pod,ptf1,pavel, pcoldry,pcolbrd,pwx,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,pcolo3,pcoln2o,pcolch4,pcolo2,p_co2mult,&
 & klaytrop,klayswtch,klaylow,pselffac,pselffrac,kindself,pfrac, &
 & kindminor,pscaleminor,pscaleminorn2,pminorfrac,&
 &  prat_h2oco2, prat_h2oco2_1, prat_h2oo3, prat_h2oo3_1, &
 &  prat_h2on2o, prat_h2on2o_1, prat_h2och4, prat_h2och4_1, &
 &  prat_n2oco2, prat_n2oco2_1, prat_o3co2, prat_o3co2_1)  

!        nec/fc        05-oct-2009 optimisation 
!     reformatted for f90 by jjmorcrette, ecmwf, 980714
!        nec           25-oct-2007 optimisations
!        d. salmond    11-dec-2007 optimizations
!     jjmorcrette 20110613 flexible number of g-points
!     abozzo  201306 update to rrtmg-lw v4.85

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband   ,jpxsec
use yoerrtm  , only : jpgpt
use yoerrtab , only : trans    ,bpade


implicit none

integer(kind=jpim),intent(in)    :: kidia
integer(kind=jpim),intent(in)    :: kfdia
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(out)   :: patr1(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(out)   :: pod(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(in)    :: pavel(kidia:kfdia,klev) ! layer pressures (pa)
real(kind=jprb)   ,intent(out)   :: ptf1(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(in)    :: pcoldry(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pwx(kidia:kfdia,jpxsec,klev) ! amount of trace gases
real(kind=jprb)   ,intent(in)    :: ptauaerl(kidia:kfdia,klev,jpband) 
real(kind=jprb)   ,intent(in)    :: pfac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pfac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pfac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pfac11(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: kjp(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: kjt(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: kjt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: poneminus

real(kind=jprb)   ,intent(in)    :: pcolh2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcolco2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcolo3(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcoln2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcolch4(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcolo2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_co2mult(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: klaytrop(kidia:kfdia) 
integer(kind=jpim),intent(in)    :: klayswtch(kidia:kfdia) 
integer(kind=jpim),intent(in)    :: klaylow(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: pselffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pselffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: kindself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: pfrac(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(in)    :: pforfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pforfrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: kindfor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pminorfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pscaleminor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pscaleminorn2(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: kindminor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: pcolbrd(kidia:kfdia,klev)  
real(kind=jprb)  , intent(in) :: &                  !
                    &   prat_h2oco2(kidia:kfdia,klev),prat_h2oco2_1(kidia:kfdia,klev), &
                    &   prat_h2oo3(kidia:kfdia,klev),prat_h2oo3_1(kidia:kfdia,klev), & !    dimensions: (nlayers)
                    &   prat_h2on2o(kidia:kfdia,klev),prat_h2on2o_1(kidia:kfdia,klev), &
                    &   prat_h2och4(kidia:kfdia,klev),prat_h2och4_1(kidia:kfdia,klev), &
                    &   prat_n2oco2(kidia:kfdia,klev),prat_n2oco2_1(kidia:kfdia,klev), &
                    &   prat_o3co2(kidia:kfdia,klev),prat_o3co2_1(kidia:kfdia,klev)
!- from aer
!- from intfac      
!- from intind
!- from precise             
!- from profdata             
!- from self             
!- from sp             
real(kind=jprb) :: ztau   (kidia:kfdia,jpgpt,klev)

integer(kind=jpim) :: ji, itr, jlev
integer(kind=jpim) :: jlon

real(kind=jprb) :: zodepth, zsecang, ztf
real(kind=jprb) :: zhook_handle


! # 1 "./include/rrtm_taumol1.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol1 (kidia,kfdia,klev,p_tau,pavel,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,p_minorfrac,k_indminor,pscaleminorn2,pcolbrd) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng1
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: pavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_minorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pscaleminorn2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolbrd(kidia:kfdia,klev)
end subroutine rrtm_taumol1
end interface
! # 93 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol10.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol10 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng10 ,ngs9
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
end subroutine rrtm_taumol10
end interface
! # 94 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol11.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol11 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,p_colo2,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,p_minorfrac,kindminor,pscaleminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng11 ,ngs10
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo2(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_minorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pscaleminor(kidia:kfdia,klev)
end subroutine rrtm_taumol11
end interface
! # 95 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol12.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol12 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colco2,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & prat_h2oco2, prat_h2oco2_1) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng12 ,ngs11
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: prat_h2oco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: prat_h2oco2_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
end subroutine rrtm_taumol12
end interface
! # 96 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol13.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol13 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_coln2o,p_colco2,p_colo3,p_coldry,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & prat_h2on2o, prat_h2on2o_1,pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng13 ,ngs12
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coln2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: prat_h2on2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: prat_h2on2o_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol13
end interface
! # 97 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol14.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol14 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colco2,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ngs13 ,ng14
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
end subroutine rrtm_taumol14
end interface
! # 98 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol15.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol15 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colco2,p_coln2o,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & prat_n2oco2, prat_n2oco2_1,pminorfrac,kindminor,pscaleminor,pcolbrd) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ngs14 ,ng15
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coln2o(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: prat_n2oco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: prat_n2oco2_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pscaleminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolbrd(kidia:kfdia,klev)
end subroutine rrtm_taumol15
end interface
! # 99 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol16.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol16 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colch4,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & p_rat_h2och4,p_rat_h2och4_1) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ngs15 ,ng16
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colch4(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2och4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2och4_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
end subroutine rrtm_taumol16
end interface
! # 100 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol2.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol2 (kidia,kfdia,klev,p_tau,pavel,p_coldry,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng2 ,ngs1
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: pavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
end subroutine rrtm_taumol2
end interface
! # 101 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol3.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol3 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colco2,p_coln2o,p_coldry,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & prat_h2oco2, prat_h2oco2_1,pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng3 ,ngs2
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coln2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: prat_h2oco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: prat_h2oco2_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol3
end interface
! # 102 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol4.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol4 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colco2,p_colo3,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & p_rat_h2oco2, p_rat_h2oco2_1, p_rat_o3co2, p_rat_o3co2_1) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng4 ,ngs3
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2oco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2oco2_1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_o3co2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_o3co2_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
end subroutine rrtm_taumol4
end interface
! # 103 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol5.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol5 (kidia,kfdia,klev,p_tau,p_wx,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colco2, p_colo3,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & p_rat_h2oco2, p_rat_h2oco2_1, p_rat_o3co2, p_rat_o3co2_1,pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband ,jpxsec
use yoerrtm , only : jpgpt ,ng5 ,ngs4
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_wx(kidia:kfdia,jpxsec,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2oco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2oco2_1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_o3co2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_o3co2_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol5
end interface
! # 104 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol6.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol6 (kidia,kfdia,klev,p_tau,p_wx,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,p_colco2,p_coldry,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband ,jpxsec
use yoerrtm , only : jpgpt ,ng6 ,ngs5
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_wx(kidia:kfdia,jpxsec,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol6
end interface
! # 105 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol7.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol7 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_colo3,p_colco2,p_coldry,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & p_rat_h2oo3, p_rat_h2oo3_1,pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng7 ,ngs6
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2oo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_rat_h2oo3_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol7
end interface
! # 106 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol8.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol8 (kidia,kfdia,klev,p_tau,p_wx,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,p_colo3,p_coln2o,p_colco2,p_coldry,k_laytrop,p_selffac,p_selffrac,k_indself,pfrac,&
 & pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband ,jpxsec
use yoerrtm , only : jpgpt ,ng8 ,ngs7
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_wx(kidia:kfdia,jpxsec,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coln2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol8
end interface
! # 107 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

! # 1 "./include/rrtm_taumol9.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine rrtm_taumol9 (kidia,kfdia,klev,p_tau,&
 & p_tauaerl,p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,p_oneminus,&
 & p_colh2o,p_coln2o,p_colch4,p_coldry,k_laytrop,k_layswtch,k_laylow,p_selffac,p_selffrac,k_indself,pfrac,&
 & prat_h2och4,prat_h2och4_1,pminorfrac,kindminor) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband
use yoerrtm , only : jpgpt ,ng9 ,ngs8
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(inout) :: p_tau(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: p_tauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_fac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_oneminus
real(kind=jprb) ,intent(in) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coln2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_colch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_laytrop(kidia:kfdia)
integer(kind=jpim),intent(in) :: k_layswtch(kidia:kfdia)
integer(kind=jpim),intent(in) :: k_laylow(kidia:kfdia)
real(kind=jprb) ,intent(in) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: prat_h2och4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: prat_h2och4_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: k_indfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_forfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
end subroutine rrtm_taumol9
end interface
! # 108 "ifsrrtm/rrtm_gasabs1a_140gp.f90" 2

!cdir duplicate(trans,256)

!- secang is equal to the secant of the diffusivity angle.
associate(nflevg=>klev)
if (lhook) call dr_hook('rrtm_gasabs1a_140gp',0,zhook_handle)
zsecang = 1.66_jprb

call rrtm_taumol1  (kidia,kfdia,klev,ztau,pavel,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolh2o,klaytrop,pselffac,pselffrac,kindself,pfrac, pminorfrac, &
 & kindminor,pscaleminorn2,pcolbrd)  
call rrtm_taumol2  (kidia,kfdia,klev,ztau,pavel,pcoldry,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolh2o,klaytrop,pselffac,pselffrac,kindself,pfrac)  
call rrtm_taumol3  (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,pcoln2o,pcoldry,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2oco2, prat_h2oco2_1,pminorfrac,kindminor)  
call rrtm_taumol4  (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,pcolo3,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2oco2, prat_h2oco2_1, prat_o3co2, prat_o3co2_1)  
call rrtm_taumol5  (kidia,kfdia,klev,ztau,pwx,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,pcolo3,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2oco2, prat_h2oco2_1, prat_o3co2, prat_o3co2_1,pminorfrac,kindminor)   
call rrtm_taumol6  (kidia,kfdia,klev,ztau,pwx,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolh2o,pcolco2,pcoldry,klaytrop,pselffac,pselffrac,kindself,pfrac,pminorfrac,kindminor)  
call rrtm_taumol7  (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolo3,pcolco2,pcoldry,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2oo3, prat_h2oo3_1,pminorfrac,kindminor)  
call rrtm_taumol8  (kidia,kfdia,klev,ztau,pwx,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolh2o,pcolo3,pcoln2o,pcolco2,pcoldry,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & pminorfrac,kindminor)  
call rrtm_taumol9  (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcoln2o,pcolch4,pcoldry,klaytrop,klayswtch,klaylow,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2och4,prat_h2och4_1,pminorfrac,kindminor)  
call rrtm_taumol10 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolh2o,klaytrop,pselffac,pselffrac,kindself,pfrac)  
call rrtm_taumol11 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolh2o,pcolo2,klaytrop,pselffac,pselffrac,kindself,pfrac,pminorfrac,kindminor,pscaleminor)  
call rrtm_taumol12 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2oco2, prat_h2oco2_1)  
call rrtm_taumol13 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcoln2o,pcolco2,pcolo3,pcoldry,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2on2o, prat_h2on2o_1,pminorfrac,kindminor)  
call rrtm_taumol14 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,&
 & pcolco2,klaytrop,pselffac,pselffrac,kindself,pfrac)  
call rrtm_taumol15 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,pcoln2o,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_n2oco2, prat_n2oco2_1,pminorfrac,kindminor,pscaleminor,pcolbrd)  
call rrtm_taumol16 (kidia,kfdia,klev,ztau,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolch4,klaytrop,pselffac,pselffrac,kindself,pfrac, &
 & prat_h2och4,prat_h2och4_1)   

!to check total od for each band
    ! print*,'ztau2= ',sum(ztau(:,11:22,:),2)
    ! print*,'ztau3= ',sum(ztau(:,23:38,:),2)
    ! print*,'ztau4= ',sum(ztau(:,39:52,:),2)
    ! print*,'ztau5= ',sum(ztau(:,53:68,:),2)
    ! print*,'ztau6= ',sum(ztau(:,69:76,:),2)
    ! print*,'ztau7= ',sum(ztau(:,77:88,:),2)
    ! print*,'ztau8= ',sum(ztau(:,89:96,:),2)
    ! print*,'ztau9= ',sum(ztau(:,97:108,:),2)
    ! print*,'ztau10= ',sum(ztau(:,109:114,:),2)
    ! print*,'ztau11= ',sum(ztau(:,115:122,:),2)
    ! print*,'ztau12= ',sum(ztau(:,123:130,:),2)
    ! print*,'ztau13= ',sum(ztau(:,131:134,:),2)
    ! print*,'ztau14= ',sum(ztau(:,135:136,:),2)
    ! print*,'ztau15= ',sum(ztau(:,137:138,:),2)
    ! print*,'ztau16= ',sum(ztau(:,139:140,:),2)


do jlev = 1, klev
!cdir unroll=4
  do ji = 1, jpgpt
    do jlon = kidia, kfdia
      if (ztau(jlon,ji,jlev) < 0._jprb) then
9101    format(1x,'gasabs jlev,ji,jlon=',i3,i5,i9,' secang=',f9.6,' ztau=',e12.6)
      endif
    enddo
  enddo
enddo


!- loop over g-channels.
do jlev = 1, klev
!cdir unroll=4
  do ji = 1, jpgpt
    do jlon = kidia, kfdia
      zodepth = zsecang * ztau(jlon,ji,jlev)
      pod(jlon,ji,jlev) = zodepth
      zodepth=0.5d0*(abs(zodepth)+zodepth)

!-- revised code to get the pre-computed transmission
!          if (odepth.le.0.) print*, 'odepth = ',odepth
!!  if (odepth <= _zero_)then
!!    atr1(ji,lay) = _one_ - trans(0)
!!    tf1(ji,lay) = _zero_
!!  else

      ztf = zodepth/(bpade+zodepth)

      itr=int(5.e+03_jprb*ztf+0.5_jprb)
      patr1(jlon,ji,jlev) = 1.0_jprb - trans(itr)
      ptf1(jlon,ji,jlev) = ztf
!!  endif
    enddo
  enddo
enddo
!     -----------------------------------------------------------------

if (lhook) call dr_hook('rrtm_gasabs1a_140gp',1,zhook_handle)
end associate
end subroutine rrtm_gasabs1a_140gp
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

