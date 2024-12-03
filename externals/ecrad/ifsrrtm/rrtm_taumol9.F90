! # 1 "ifsrrtm/rrtm_taumol9.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_taumol9.f90"
! this file has been modified for the use in icon

!----------------------------------------------------------------------------
subroutine rrtm_taumol9 (kidia,kfdia,klev,taug,&
 & p_tauaerl,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,coln2o,colch4,coldry,laytrop,k_layswtch,k_laylow,selffac,selffrac,indself,fracs, &
 & rat_h2och4,rat_h2och4_1,minorfrac,indminor)  

!     band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)

!     author.
!     -------
!      jjmorcrette, ecmwf

!     modifications.
!     --------------
!      m.hamrud      01-oct-2003 cy28 cleaning
!      nec           25-oct-2007 optimisations
!      jjmorcrette 20110613 flexible number of g-points
!      201306 abozzo updated to rrtmg v4.85
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)
! ---------------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband
use yoerrtm  , only : jpgpt  ,ng9   ,ngs8
use yoerrtwn , only :      nspa   ,nspb
use yoerrta9 , only : absa   ,absb   ,fracrefa, fracrefb,selfref,forref,ka_mn2o, kb_mn2o   
use yoerrtrf, only : chi_mls

implicit none

integer(kind=jpim),intent(in)    :: kidia
integer(kind=jpim),intent(in)    :: kfdia
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(inout) :: taug(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(in)    :: p_tauaerl(kidia:kfdia,klev,jpband) 
real(kind=jprb)   ,intent(in)    :: fac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac11(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jp(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jt(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: oneminus
real(kind=jprb)   ,intent(in)    :: colh2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: coln2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: colch4(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: coldry(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: laytrop(kidia:kfdia) 
integer(kind=jpim),intent(in)    :: k_layswtch(kidia:kfdia) 
integer(kind=jpim),intent(in)    :: k_laylow(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: fracs(kidia:kfdia,jpgpt,klev) 

real(kind=jprb)   ,intent(in)   :: rat_h2och4(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: rat_h2och4_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indfor(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: forfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: minorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indminor(kidia:kfdia,klev)



! ---------------------------------------------------------------------------

integer(kind=jpim) :: ind0,ind1,inds,indf,indm
integer(kind=jpim) :: ig, js, lay, js1, jmn2o, jpl

real(kind=jprb) :: speccomb,speccomb1,speccomb_mn2o, &
& speccomb_planck
real(kind=jprb) :: refrat_planck_a, refrat_m_a

real(kind=jprb) :: fs, specmult, specparm,  &
 & fs1, specmult1, specparm1,   &
 & fmn2o, specmult_mn2o, specparm_mn2o,   &
 & fpl, specmult_planck, specparm_planck

real(kind=jprb) :: adjfac,adjcoln2o,ratn2o,chi_n2o
real(kind=jprb) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
real(kind=jprb) :: p, p4, fk0, fk1, fk2
real(kind=jprb) :: taufor,tauself,n2om1,n2om2,absn2o,tau_major(ng9),tau_major1(ng9)

    !     local integer arrays
    integer(kind=jpim) :: laytrop_min, laytrop_max
    integer(kind=jpim) :: ixc(klev), ixlow(kfdia,klev), ixhigh(kfdia,klev)
    integer(kind=jpim) :: ich, icl, ixc0, ixp, jc, jl



    !$acc data present(taug, p_tauaerl, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$acc             colh2o, coln2o, colch4, coldry, laytrop, k_layswtch, &
    !$acc             k_laylow, selffac, selffrac, indself, fracs, rat_h2och4, &
    !$acc             rat_h2och4_1, indfor, forfac, forfrac, minorfrac, indminor)


    laytrop_min = minval(laytrop)
    laytrop_max = maxval(laytrop)

    ixlow  = 0
    ixhigh = 0
    ixc    = 0

    ! create index lists for mixed layers
    do lay = laytrop_min+1, laytrop_max
      icl = 0
      ich = 0
      do jc = kidia, kfdia
        if ( lay <= laytrop(jc) ) then
          icl = icl + 1
          ixlow(icl,lay) = jc
        else
          ich = ich + 1
          ixhigh(ich,lay) = jc
        endif
      enddo
      ixc(lay) = icl
    enddo
! # 139 "ifsrrtm/rrtm_taumol9.f90"


      ! p = 212 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

      ! p = 706.272 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)

      ! compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! lower atmosphere loop
      !$acc wait
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, inds, indf, indm, js, js1, jmn2o, jpl, speccomb, & 
      !$acc   speccomb1, speccomb_mn2o, speccomb_planck, fs, specmult, specparm, fs1, specmult1, specparm1, fmn2o, &
      !$acc   specmult_mn2o, specparm_mn2o, fpl, specmult_planck, specparm_planck, adjfac, adjcoln2o, ratn2o, &
      !$acc   chi_n2o, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, &
      !$acc   fac211, p, p4, fk0, fk1, fk2, tau_major, tau_major1)
      do lay = 1, laytrop_min
        do jl = kidia, kfdia

          speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
          specparm = min(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._jprb*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - aint((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
          specparm1 = min(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._jprb*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - aint((specmult1)))

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colch4(jl,lay)
          specparm_mn2o = min(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._jprb*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - aint((specmult_mn2o)))

          !  in atmospheres where the amount of n2o is too great to be considered
          !  a minor species, adjust the column amount of n2o by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_jprb*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_jprb) then
            adjfac = 0.5_jprb+(ratn2o-0.5_jprb)**0.65_jprb
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = min(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._jprb*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - aint((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(9) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(9) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_jprb) then
            p = fs - 1._jprb
            p4 = p**4
            fk0 = p4
            fk1 = 1._jprb - p - 2.0_jprb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_jprb) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._jprb - p - 2.0_jprb*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._jprb - fs) * fac00(jl,lay)
            fac010 = (1._jprb - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._jprb
            fac210 = 0._jprb
          endif

          if (specparm1 .lt. 0.125_jprb) then
            p = fs1 - 1._jprb
            p4 = p**4
            fk0 = p4
            fk1 = 1._jprb - p - 2.0_jprb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_jprb) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._jprb - p - 2.0_jprb*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._jprb - fs1) * fac01(jl,lay)
            fac011 = (1._jprb - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._jprb
            fac211 = 0._jprb
          endif

          if (specparm .lt. 0.125_jprb) then
!$nec unroll(ng9)
            tau_major(1:ng9) = speccomb *    &
             (fac000 * absa(ind0,1:ng9)    + &
              fac100 * absa(ind0+1,1:ng9)  + &
              fac200 * absa(ind0+2,1:ng9)  + &
              fac010 * absa(ind0+9,1:ng9)  + &
              fac110 * absa(ind0+10,1:ng9) + &
              fac210 * absa(ind0+11,1:ng9))
          else if (specparm .gt. 0.875_jprb) then
!$nec unroll(ng9)
            tau_major(1:ng9) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng9) + &
              fac100 * absa(ind0,1:ng9)   + &
              fac000 * absa(ind0+1,1:ng9) + &
              fac210 * absa(ind0+8,1:ng9) + &
              fac110 * absa(ind0+9,1:ng9) + &
              fac010 * absa(ind0+10,1:ng9))
          else
!$nec unroll(ng9)
            tau_major(1:ng9) = speccomb *   &
             (fac000 * absa(ind0,1:ng9)   + &
              fac100 * absa(ind0+1,1:ng9) + &
              fac010 * absa(ind0+9,1:ng9) + &
              fac110 * absa(ind0+10,1:ng9))
          endif

          if (specparm1 .lt. 0.125_jprb) then
!$nec unroll(ng9)
            tau_major1(1:ng9) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng9)    + &
              fac101 * absa(ind1+1,1:ng9)  + &
              fac201 * absa(ind1+2,1:ng9)  + &
              fac011 * absa(ind1+9,1:ng9)  + &
              fac111 * absa(ind1+10,1:ng9) + &
              fac211 * absa(ind1+11,1:ng9))
          else if (specparm1 .gt. 0.875_jprb) then
!$nec unroll(ng9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng9) + &
              fac101 * absa(ind1,1:ng9)   + &
              fac001 * absa(ind1+1,1:ng9) + &
              fac211 * absa(ind1+8,1:ng9) + &
              fac111 * absa(ind1+9,1:ng9) + &
              fac011 * absa(ind1+10,1:ng9))
          else
!$nec unroll(ng9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac001 * absa(ind1,1:ng9)   + &
              fac101 * absa(ind1+1,1:ng9) + &
              fac011 * absa(ind1+9,1:ng9) + &
              fac111 * absa(ind1+10,1:ng9))
          endif

          !$acc loop seq private(taufor, tauself, n2om1, n2om2, absn2o)
!$nec unroll(ng9)
          do ig = 1, ng9
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

            taug(jl,ngs8+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,ngs8+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$acc end parallel

      ! upper atmosphere loop
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, indm, adjfac, adjcoln2o, ratn2o, chi_n2o) 
      do lay = laytrop_max+1, klev
        do jl = kidia, kfdia

          !  in atmospheres where the amount of n2o is too great to be considered
          !  a minor species, adjust the column amount of n2o by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_jprb*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_jprb) then
            adjfac = 0.5_jprb+(ratn2o-0.5_jprb)**0.65_jprb
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(9) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(9) + 1
          indm = indminor(jl,lay)
          !$acc loop seq private(absn2o)
!$nec unroll(ng9)
          do ig = 1, ng9
            absn2o = kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
            taug(jl,ngs8+ig,lay) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracs(jl,ngs8+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo
      !$acc end parallel

      if (laytrop_max /= laytrop_min) then
        ! mixed loop
        ! lower atmosphere part
        !$acc parallel default(none) async(1)
        !$acc loop gang vector collapse(2) private(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$acc   specmult1, js1, fs1, speccomb_mn2o, specparm_mn2o, specmult_mn2o, jmn2o, fmn2o, chi_n2o, ratn2o, &
        !$acc   adjfac, adjcoln2o, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$acc   indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, &
        !$acc   fac201, fac011, fac111, fac211, tau_major, tau_major1)
        do lay = laytrop_min+1, laytrop_max




          ixc0 = ixc(lay)

!$nec ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)


            speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
            specparm = min(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._jprb*(specparm)
            js = 1 + int(specmult)
            fs = ((specmult) - aint((specmult)))

            speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
            specparm1 = min(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._jprb*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = ((specmult1) - aint((specmult1)))

            speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colch4(jl,lay)
            specparm_mn2o = min(colh2o(jl,lay)/speccomb_mn2o,oneminus)
            specmult_mn2o = 8._jprb*specparm_mn2o
            jmn2o = 1 + int(specmult_mn2o)
            fmn2o = ((specmult_mn2o) - aint((specmult_mn2o)))

            !  in atmospheres where the amount of n2o is too great to be considered
            !  a minor species, adjust the column amount of n2o by an empirical factor
            !  to obtain the proper contribution.
            chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
            ratn2o = 1.e20_jprb*chi_n2o/chi_mls(4,jp(jl,lay)+1)
            if (ratn2o .gt. 1.5_jprb) then
              adjfac = 0.5_jprb+(ratn2o-0.5_jprb)**0.65_jprb
              adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
            else
              adjcoln2o = coln2o(jl,lay)
            endif

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
            specparm_planck = min(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._jprb*specparm_planck
            jpl= 1 + int(specmult_planck)
            fpl = ((specmult_planck) - aint((specmult_planck)))

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(9) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(9) + js1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)

            if (specparm .lt. 0.125_jprb) then
              p = fs - 1._jprb
              p4 = p**4
              fk0 = p4
              fk1 = 1._jprb - p - 2.0_jprb*p4
              fk2 = p + p4
              fac000 = fk0*fac00(jl,lay)
              fac100 = fk1*fac00(jl,lay)
              fac200 = fk2*fac00(jl,lay)
              fac010 = fk0*fac10(jl,lay)
              fac110 = fk1*fac10(jl,lay)
              fac210 = fk2*fac10(jl,lay)
            else if (specparm .gt. 0.875_jprb) then
              p = -fs
              p4 = p**4
              fk0 = p4
              fk1 = 1._jprb - p - 2.0_jprb*p4
              fk2 = p + p4
              fac000 = fk0*fac00(jl,lay)
              fac100 = fk1*fac00(jl,lay)
              fac200 = fk2*fac00(jl,lay)
              fac010 = fk0*fac10(jl,lay)
              fac110 = fk1*fac10(jl,lay)
              fac210 = fk2*fac10(jl,lay)
            else
              fac000 = (1._jprb - fs) * fac00(jl,lay)
              fac010 = (1._jprb - fs) * fac10(jl,lay)
              fac100 = fs * fac00(jl,lay)
              fac110 = fs * fac10(jl,lay)
              fac200 = 0._jprb
              fac210 = 0._jprb
            endif

            if (specparm1 .lt. 0.125_jprb) then
              p = fs1 - 1._jprb
              p4 = p**4
              fk0 = p4
              fk1 = 1._jprb - p - 2.0_jprb*p4
              fk2 = p + p4
              fac001 = fk0*fac01(jl,lay)
              fac101 = fk1*fac01(jl,lay)
              fac201 = fk2*fac01(jl,lay)
              fac011 = fk0*fac11(jl,lay)
              fac111 = fk1*fac11(jl,lay)
              fac211 = fk2*fac11(jl,lay)
            else if (specparm1 .gt. 0.875_jprb) then
              p = -fs1
              p4 = p**4
              fk0 = p4
              fk1 = 1._jprb - p - 2.0_jprb*p4
              fk2 = p + p4
              fac001 = fk0*fac01(jl,lay)
              fac101 = fk1*fac01(jl,lay)
              fac201 = fk2*fac01(jl,lay)
              fac011 = fk0*fac11(jl,lay)
              fac111 = fk1*fac11(jl,lay)
              fac211 = fk2*fac11(jl,lay)
            else
              fac001 = (1._jprb - fs1) * fac01(jl,lay)
              fac011 = (1._jprb - fs1) * fac11(jl,lay)
              fac101 = fs1 * fac01(jl,lay)
              fac111 = fs1 * fac11(jl,lay)
              fac201 = 0._jprb
              fac211 = 0._jprb
            endif

            if (specparm .lt. 0.125_jprb) then
!$nec unroll(ng9)
              tau_major(1:ng9) = speccomb *    &
              (fac000 * absa(ind0,1:ng9)    + &
                fac100 * absa(ind0+1,1:ng9)  + &
                fac200 * absa(ind0+2,1:ng9)  + &
                fac010 * absa(ind0+9,1:ng9)  + &
                fac110 * absa(ind0+10,1:ng9) + &
                fac210 * absa(ind0+11,1:ng9))
            else if (specparm .gt. 0.875_jprb) then
!$nec unroll(ng9)
              tau_major(1:ng9) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng9) + &
                fac100 * absa(ind0,1:ng9)   + &
                fac000 * absa(ind0+1,1:ng9) + &
                fac210 * absa(ind0+8,1:ng9) + &
                fac110 * absa(ind0+9,1:ng9) + &
                fac010 * absa(ind0+10,1:ng9))
            else
!$nec unroll(ng9)
              tau_major(1:ng9) = speccomb *   &
              (fac000 * absa(ind0,1:ng9)   + &
                fac100 * absa(ind0+1,1:ng9) + &
                fac010 * absa(ind0+9,1:ng9) + &
                fac110 * absa(ind0+10,1:ng9))
            endif

            if (specparm1 .lt. 0.125_jprb) then
!$nec unroll(ng9)
              tau_major1(1:ng9) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng9)    + &
                fac101 * absa(ind1+1,1:ng9)  + &
                fac201 * absa(ind1+2,1:ng9)  + &
                fac011 * absa(ind1+9,1:ng9)  + &
                fac111 * absa(ind1+10,1:ng9) + &
                fac211 * absa(ind1+11,1:ng9))
            else if (specparm1 .gt. 0.875_jprb) then
!$nec unroll(ng9)
              tau_major1(1:ng9) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng9) + &
                fac101 * absa(ind1,1:ng9)   + &
                fac001 * absa(ind1+1,1:ng9) + &
                fac211 * absa(ind1+8,1:ng9) + &
                fac111 * absa(ind1+9,1:ng9) + &
                fac011 * absa(ind1+10,1:ng9))
            else
!$nec unroll(ng9)
              tau_major1(1:ng9) = speccomb1 * &
              (fac001 * absa(ind1,1:ng9)   + &
                fac101 * absa(ind1+1,1:ng9) + &
                fac011 * absa(ind1+9,1:ng9) + &
                fac111 * absa(ind1+10,1:ng9))
            endif

!$nec unroll(ng9)
            !$acc loop seq private(tauself, taufor, n2om1, n2om2, absn2o)
            do ig = 1, ng9
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                  (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
              n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                  (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
              absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

              taug(jl,ngs8+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + adjcoln2o*absn2o
              fracs(jl,ngs8+ig,lay) = fracrefa(ig,jpl) + fpl * &
                  (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
            enddo



          enddo

          ! upper atmosphere part
          ixc0 = kfdia - kidia + 1 - ixc0
!$nec ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)


            !  in atmospheres where the amount of n2o is too great to be considered
            !  a minor species, adjust the column amount of n2o by an empirical factor
            !  to obtain the proper contribution.
            chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
            ratn2o = 1.e20_jprb*chi_n2o/chi_mls(4,jp(jl,lay)+1)
            if (ratn2o .gt. 1.5_jprb) then
              adjfac = 0.5_jprb+(ratn2o-0.5_jprb)**0.65_jprb
              adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
            else
              adjcoln2o = coln2o(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(9) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(9) + 1
            indm = indminor(jl,lay)
!$nec unroll(ng9)
            !$acc loop seq private(absn2o)
            do ig = 1, ng9
              absn2o = kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
              taug(jl,ngs8+ig,lay) = colch4(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) +  &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + adjcoln2o*absn2o
              fracs(jl,ngs8+ig,lay) = fracrefb(ig)
            enddo



          enddo

        enddo
        !$acc end parallel

      endif

      !$acc end data
      
end subroutine rrtm_taumol9
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define mod1(x) ((x) - aint((x)))
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

