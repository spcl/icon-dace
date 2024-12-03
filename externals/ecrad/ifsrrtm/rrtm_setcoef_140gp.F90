! # 1 "ifsrrtm/rrtm_setcoef_140gp.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_setcoef_140gp.f90"
! this file has been modified for the use in icon

subroutine rrtm_setcoef_140gp (kidia,kfdia,klev,p_coldry,p_wbroad,p_wkl,&
 & p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,p_colco2,p_colo3,p_coln2o,p_colch4, p_colo2,p_co2mult, p_colbrd, &
 & k_laytrop,k_layswtch,k_laylow,pavel,p_tavel,p_selffac,p_selffrac,k_indself,&
 & k_indminor,p_scaleminor,p_scaleminorn2,p_minorfrac,&
 & prat_h2oco2, prat_h2oco2_1, prat_h2oo3, prat_h2oo3_1, &
 & prat_h2on2o, prat_h2on2o_1, prat_h2och4, prat_h2och4_1, &
 & prat_n2oco2, prat_n2oco2_1, prat_o3co2, prat_o3co2_1)

!     reformatted for f90 by jjmorcrette, ecmwf, 980714
!        nec           25-oct-2007 optimisations
!     201305 abozzo updated to rrtmg_lw_v4.85
!     201507 rhogan bug fix: swapped p_colo2 & p_co2mult in argument list


!     purpose:  for a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.
!     also calculate the values of the integrated planck functions 
!     for each band at the level and layer temperatures.

use parkind1 , only : jpim, jprb
use ecradhook  , only : lhook, dr_hook, jphook
use parrrtm  , only : jpinpx
use yoerrtrf , only : preflog   ,tref, chi_mls

implicit none

integer(kind=jpim),intent(in)    :: kidia
integer(kind=jpim),intent(in)    :: kfdia
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(in)    :: p_coldry(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_wbroad(kidia:kfdia,klev)
real(kind=jprb)   ,intent(out)   :: p_colbrd(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: p_wkl(kidia:kfdia,jpinpx,klev) 
real(kind=jprb)   ,intent(out)   :: p_fac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_fac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_fac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_fac11(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(out)   :: k_jp(kidia:kfdia,klev) 
integer(kind=jpim),intent(out)   :: k_jt(kidia:kfdia,klev) 
integer(kind=jpim),intent(out)   :: k_jt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_colh2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_colco2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_colo3(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_coln2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_colch4(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_colo2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_co2mult(kidia:kfdia,klev) 
integer(kind=jpim),intent(out)   :: k_laytrop(kidia:kfdia) 
integer(kind=jpim),intent(out)   :: k_layswtch(kidia:kfdia) 
integer(kind=jpim),intent(out)   :: k_laylow(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: pavel(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_tavel(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(out)   :: k_indself(kidia:kfdia,klev) 
integer(kind=jpim),intent(out)   :: k_indfor(kidia:kfdia,klev) 
integer(kind=jpim),intent(out)   :: k_indminor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_scaleminor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_scaleminorn2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: p_minorfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(out)   :: &                 !
                     & prat_h2oco2(kidia:kfdia,klev),prat_h2oco2_1(kidia:kfdia,klev), &
                     & prat_h2oo3(kidia:kfdia,klev) ,prat_h2oo3_1(kidia:kfdia,klev), & 
                     & prat_h2on2o(kidia:kfdia,klev),prat_h2on2o_1(kidia:kfdia,klev), &
                     & prat_h2och4(kidia:kfdia,klev),prat_h2och4_1(kidia:kfdia,klev), &
                     & prat_n2oco2(kidia:kfdia,klev),prat_n2oco2_1(kidia:kfdia,klev), &
                     & prat_o3co2(kidia:kfdia,klev) ,prat_o3co2_1(kidia:kfdia,klev)
!- from intfac      
!- from intind
!- from profdata             
!- from profile             
!- from self             
integer(kind=jpim) :: jp1, jlay
integer(kind=jpim) :: jlon

real(kind=jprb) :: z_co2reg, z_compfp, z_factor, z_fp, z_ft, z_ft1, z_plog, z_scalefac, z_stpfac, z_water
real(kind=jphook) :: zhook_handle

if (lhook) call dr_hook('rrtm_setcoef_140gp',0,zhook_handle)

!$acc parallel default(none) present(p_coldry, p_wbroad, p_colbrd, p_wkl, p_fac00, p_fac01, p_fac10, p_fac11, &
!$acc   p_forfac, p_forfrac, k_jp, k_jt, k_jt1, p_colh2o, p_colco2, p_colo3, p_coln2o, p_colch4, p_colo2, p_co2mult, &
!$acc   k_laytrop, k_layswtch, k_laylow, pavel, p_tavel, p_selffac, p_selffrac, k_indself, k_indfor, k_indminor, &
!$acc   p_scaleminor, p_scaleminorn2, p_minorfrac, prat_h2oco2, prat_h2oco2_1, prat_h2oo3, prat_h2oo3_1, prat_h2on2o, &
!$acc   prat_h2on2o_1, prat_h2och4, prat_h2och4_1, prat_n2oco2, prat_n2oco2_1, prat_o3co2, prat_o3co2_1) async(1)
!$acc loop gang vector
do jlon = kidia, kfdia
  z_stpfac = 296._jprb/1013._jprb

  k_laytrop(jlon)  = 0
  k_layswtch(jlon) = 0
  k_laylow(jlon)   = 0
  !$acc loop seq
  do jlay = 1, klev
!        find the two reference pressures on either side of the
!        layer pressure.  store them in jp and jp1.  store in fp the
!        fraction of the difference (in ln(pressure)) between these
!        two values that the layer pressure lies.
    z_plog = log(pavel(jlon,jlay))
    k_jp(jlon,jlay) = int(36._jprb - 5*(z_plog+0.04_jprb))
    if (k_jp(jlon,jlay)  <  1) then
      k_jp(jlon,jlay) = 1
    elseif (k_jp(jlon,jlay)  >  58) then
      k_jp(jlon,jlay) = 58
    endif
    jp1 = k_jp(jlon,jlay) + 1
    z_fp = 5._jprb * (preflog(k_jp(jlon,jlay)) - z_plog)
! bound z_fp in case z_plog is outside range of ref. pressure preflog
! (in lvertfe, pressure at last full level is known, but not in finite diff (nh)
    z_fp = max(-1.0_jprb, min(1.0_jprb, z_fp))
!        determine, for each reference pressure (jp and jp1), which
!        reference temperature (these are different for each  
!        reference pressure) is nearest the layer temperature but does
!        not exceed it.  store these indices in jt and jt1, resp.
!        store in ft (resp. ft1) the fraction of the way between jt
!        (jt1) and the next highest reference temperature that the 
!        layer temperature falls.

    k_jt(jlon,jlay) = int(3._jprb + (p_tavel(jlon,jlay)-tref(k_jp(jlon,jlay)))/15._jprb)
    if (k_jt(jlon,jlay)  <  1) then
      k_jt(jlon,jlay) = 1
    elseif (k_jt(jlon,jlay)  >  4) then
      k_jt(jlon,jlay) = 4
    endif
    z_ft = ((p_tavel(jlon,jlay)-tref(k_jp(jlon,jlay)))/15._jprb) - real(k_jt(jlon,jlay)-3)
    k_jt1(jlon,jlay) = int(3._jprb + (p_tavel(jlon,jlay)-tref(jp1))/15._jprb)
    if (k_jt1(jlon,jlay)  <  1) then
      k_jt1(jlon,jlay) = 1
    elseif (k_jt1(jlon,jlay)  >  4) then
      k_jt1(jlon,jlay) = 4
    endif
    z_ft1 = ((p_tavel(jlon,jlay)-tref(jp1))/15._jprb) - real(k_jt1(jlon,jlay)-3)

    z_water = p_wkl(jlon,1,jlay)/p_coldry(jlon,jlay)
    z_scalefac = pavel(jlon,jlay) * z_stpfac / p_tavel(jlon,jlay)

!        if the pressure is less than ~100mb, perform a different
!        set of species interpolations.
!         if (plog .le. 4.56) go to 5300
!--------------------------------------         
    if (z_plog  >  4.56_jprb) then
      k_laytrop(jlon) =  k_laytrop(jlon) + 1
!        for one band, the "switch" occurs at ~300 mb. 
!      if (z_plog  >=  5.76_jprb) k_layswtch(jlon) = k_layswtch(jlon) + 1
!      if (z_plog  >=  6.62_jprb) k_laylow(jlon) = k_laylow(jlon) + 1


!        water vapor foreign continuum
      p_forfac(jlon,jlay) = z_scalefac / (1.0_jprb+z_water)
      z_factor = (332.0_jprb-p_tavel(jlon,jlay))/36.0_jprb
      k_indfor(jlon,jlay) = min(2, max(1, int(z_factor)))
      p_forfrac(jlon,jlay) = z_factor - real(k_indfor(jlon,jlay))

!        set up factors needed to separately include the water vapor
!        self-continuum in the calculation of absorption coefficient.
!c           selffac(lay) = water * scalefac / (1.+water)
      p_selffac(jlon,jlay) = z_water * p_forfac(jlon,jlay)
      z_factor = (p_tavel(jlon,jlay)-188.0_jprb)/7.2_jprb
      k_indself(jlon,jlay) = min(9, max(1, int(z_factor)-7))
      p_selffrac(jlon,jlay) = z_factor - real(k_indself(jlon,jlay) + 7)

!  set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
         p_scaleminor(jlon,jlay) = pavel(jlon,jlay)/p_tavel(jlon,jlay)
         p_scaleminorn2(jlon,jlay) = (pavel(jlon,jlay)/p_tavel(jlon,jlay)) &
           &  *(p_wbroad(jlon,jlay)/(p_coldry(jlon,jlay)+p_wkl(jlon,1,jlay)))
         z_factor = (p_tavel(jlon,jlay)-180.8_jprb)/7.2_jprb
         k_indminor(jlon,jlay) = min(18, max(1, int(z_factor)))
         p_minorfrac(jlon,jlay) = z_factor - real(k_indminor(jlon,jlay))

!  setup reference ratio to be used in calculation of binary
!  species parameter in lower atmosphere.
         prat_h2oco2(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay))/chi_mls(2,k_jp(jlon,jlay))
         prat_h2oco2_1(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay)+1)/chi_mls(2,k_jp(jlon,jlay)+1)

         prat_h2oo3(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay))/chi_mls(3,k_jp(jlon,jlay))
         prat_h2oo3_1(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay)+1)/chi_mls(3,k_jp(jlon,jlay)+1)

         prat_h2on2o(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay))/chi_mls(4,k_jp(jlon,jlay))
         prat_h2on2o_1(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay)+1)/chi_mls(4,k_jp(jlon,jlay)+1)

         prat_h2och4(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay))/chi_mls(6,k_jp(jlon,jlay))
         prat_h2och4_1(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay)+1)/chi_mls(6,k_jp(jlon,jlay)+1)

         prat_n2oco2(jlon,jlay)=chi_mls(4,k_jp(jlon,jlay))/chi_mls(2,k_jp(jlon,jlay))
         prat_n2oco2_1(jlon,jlay)=chi_mls(4,k_jp(jlon,jlay)+1)/chi_mls(2,k_jp(jlon,jlay)+1)



!        calculate needed column amounts.
      p_colh2o(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,1,jlay)
      p_colco2(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,2,jlay)
      p_colo3(jlon,jlay)  = 1.e-20_jprb * p_wkl(jlon,3,jlay)
      p_coln2o(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,4,jlay)
      p_colch4(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,6,jlay)
      p_colo2(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,7,jlay)
      p_colbrd(jlon,jlay) = 1.e-20_jprb * p_wbroad(jlon,jlay)
      if (p_colco2(jlon,jlay)  ==  0.0_jprb) p_colco2(jlon,jlay) = 1.e-32_jprb * p_coldry(jlon,jlay)
      if (p_coln2o(jlon,jlay)  ==  0.0_jprb) p_coln2o(jlon,jlay) = 1.e-32_jprb * p_coldry(jlon,jlay)
      if (p_colch4(jlon,jlay)  ==  0.0_jprb) p_colch4(jlon,jlay) = 1.e-32_jprb * p_coldry(jlon,jlay)
!        using e = 1334.2 cm-1.
      z_co2reg = 3.55e-24_jprb * p_coldry(jlon,jlay)
      p_co2mult(jlon,jlay)= (p_colco2(jlon,jlay) - z_co2reg) *&
       & 272.63_jprb*exp(-1919.4_jprb/p_tavel(jlon,jlay))/(8.7604e-4_jprb*p_tavel(jlon,jlay))  
!         go to 5400
!------------------
    else
!        above laytrop.
! 5300    continue

!        calculate needed column amounts.
      p_forfac(jlon,jlay) = z_scalefac / (1.0_jprb+z_water)
      z_factor = (p_tavel(jlon,jlay)-188.0_jprb)/36.0_jprb
      k_indfor(jlon,jlay) = 3
      p_forfrac(jlon,jlay) = z_factor - 1.0_jprb

!  set up factors needed to separately include the water vapor
!  self-continuum in the calculation of absorption coefficient.
      p_selffac(jlon,jlay) = z_water * p_forfac(jlon,jlay)

!  set up factors needed to separately include the minor gases
!  in the calculation of absorption coefficient
      p_scaleminor(jlon,jlay) = pavel(jlon,jlay)/p_tavel(jlon,jlay)         
      p_scaleminorn2(jlon,jlay) = (pavel(jlon,jlay)/p_tavel(jlon,jlay)) &
        &    * (p_wbroad(jlon,jlay)/(p_coldry(jlon,jlay)+p_wkl(jlon,1,jlay)))
      z_factor = (p_tavel(jlon,jlay)-180.8_jprb)/7.2_jprb
      k_indminor(jlon,jlay) = min(18, max(1, int(z_factor)))
      p_minorfrac(jlon,jlay) = z_factor - real(k_indminor(jlon,jlay))

!  setup reference ratio to be used in calculation of binary
!  species parameter in upper atmosphere.
      prat_h2oco2(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay))/chi_mls(2,k_jp(jlon,jlay))
      prat_h2oco2_1(jlon,jlay)=chi_mls(1,k_jp(jlon,jlay)+1)/chi_mls(2,k_jp(jlon,jlay)+1)         

      prat_o3co2(jlon,jlay)=chi_mls(3,k_jp(jlon,jlay))/chi_mls(2,k_jp(jlon,jlay))
      prat_o3co2_1(jlon,jlay)=chi_mls(3,k_jp(jlon,jlay)+1)/chi_mls(2,k_jp(jlon,jlay)+1)         


!  calculate needed column amounts.
      p_colh2o(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,1,jlay)
      p_colco2(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,2,jlay)
      p_colo3(jlon,jlay)  = 1.e-20_jprb * p_wkl(jlon,3,jlay)
      p_coln2o(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,4,jlay)
      p_colch4(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,6,jlay)
      p_colo2(jlon,jlay) = 1.e-20_jprb * p_wkl(jlon,7,jlay)
      p_colbrd(jlon,jlay) = 1.e-20_jprb * p_wbroad(jlon,jlay)
      if (p_colco2(jlon,jlay)  ==  0.0_jprb) p_colco2(jlon,jlay) = 1.e-32_jprb * p_coldry(jlon,jlay)
      if (p_coln2o(jlon,jlay)  ==  0.0_jprb) p_coln2o(jlon,jlay) = 1.e-32_jprb * p_coldry(jlon,jlay)
      if (p_colch4(jlon,jlay)  ==  0.0_jprb) p_colch4(jlon,jlay) = 1.e-32_jprb * p_coldry(jlon,jlay)
      z_co2reg = 3.55e-24_jprb * p_coldry(jlon,jlay)
      p_co2mult(jlon,jlay)= (p_colco2(jlon,jlay) - z_co2reg) *&
       & 272.63_jprb*exp(-1919.4_jprb/p_tavel(jlon,jlay))/(8.7604e-4_jprb*p_tavel(jlon,jlay))  
!----------------     
    endif
! 5400    continue

!        we have now isolated the layer ln pressure and temperature,
!        between two reference pressures and two reference temperatures 
!        (for each reference pressure).  we multiply the pressure 
!        fraction fp with the appropriate temperature fractions to get 
!        the factors that will be needed for the interpolation that yields
!        the optical depths (performed in routines taugbn for band n).

    z_compfp = 1.0_jprb - z_fp
    p_fac10(jlon,jlay) = z_compfp * z_ft
    p_fac00(jlon,jlay) = z_compfp * (1.0_jprb - z_ft)
    p_fac11(jlon,jlay) = z_fp * z_ft1
    p_fac01(jlon,jlay) = z_fp * (1.0_jprb - z_ft1)

!  rescale selffac and forfac for use in taumol
    p_selffac(jlon,jlay) = p_colh2o(jlon,jlay)*p_selffac(jlon,jlay)
    p_forfac(jlon,jlay) = p_colh2o(jlon,jlay)*p_forfac(jlon,jlay)


  enddo

! mt 981104 
!-- set laylow for profiles with surface pressure less than 750 hpa. 
  if (k_laylow(jlon) == 0) k_laylow(jlon)=1
enddo
!$acc end parallel

if (lhook) call dr_hook('rrtm_setcoef_140gp',1,zhook_handle)

end subroutine rrtm_setcoef_140gp
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

