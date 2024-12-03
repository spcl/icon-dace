! # 1 "ifsrrtm/rrtm_taumol13.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_taumol13.f90"
! this file has been modified for the use in icon

!----------------------------------------------------------------------------
subroutine rrtm_taumol13 (kidia,kfdia,klev,taug,&
 & p_tauaerl,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,coln2o,colco2,colo3,coldry,laytrop,selffac,selffrac,indself,fracs, &
 & rat_h2on2o, rat_h2on2o_1,minorfrac,indminor)  

!     band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)

!     author.
!     -------
!      jjmorcrette, ecmwf

  !     modifications.
!     --------------
!      m.hamrud      01-oct-2003 cy28 cleaning
!      nec           25-oct-2007 optimisations
!      jjmorcrette 20110613 flexible number of g-points
!      abozzo 201306 updated to rrtmg v4.85
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
! ---------------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband
use yoerrtm  , only : jpgpt  ,ng13  ,ngs12
use yoerrtwn , only : nspa   
use yoerrta13, only : absa   ,fracrefa,fracrefb,selfref,forref,ka_mco2, ka_mco, kb_mo3
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
real(kind=jprb)   ,intent(in)    :: colco2(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: colo3(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)    :: coldry(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: laytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: fracs(kidia:kfdia,jpgpt,klev) 

real(kind=jprb)   ,intent(in)   :: rat_h2on2o(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: rat_h2on2o_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indfor(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: forfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: minorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indminor(kidia:kfdia,klev)

! ---------------------------------------------------------------------------


real(kind=jprb) :: speccomb,speccomb1, speccomb_planck, &
                   & speccomb_mco2, speccomb_mco
integer(kind=jpim) :: ind0,ind1,inds,indf,indm
integer(kind=jpim) :: ig, js, lay, js1, jpl, jmco2, jmco

real(kind=jprb) :: refrat_planck_a, refrat_m_a, refrat_m_a3
real(kind=jprb) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
real(kind=jprb) :: p, p4, fk0, fk1, fk2

real(kind=jprb) :: taufor,tauself,tau_major(ng13),tau_major1(ng13), co2m1, co2m2, absco2
real(kind=jprb) :: com1, com2, absco, abso3
real(kind=jprb) :: chi_co2, ratco2, adjfac, adjcolco2

real(kind=jprb) :: fs, specmult, specparm,  &
& fs1, specmult1, specparm1, &
& fmco2, specmult_mco2, specparm_mco2, &
& fmco , specmult_mco , specparm_mco , &
& fpl, specmult_planck, specparm_planck

real(kind=jprb)   :: colco(kidia:kfdia,klev) !left =0 for now,not passed from rrtm_gasbas1a

    !     local integer arrays
    integer(kind=jpim) :: laytrop_min, laytrop_max
    integer(kind=jpim) :: ixc(klev), ixlow(kfdia,klev), ixhigh(kfdia,klev)
    integer(kind=jpim) :: ich, icl, ixc0, ixp, jc, jl




    !$acc data present(taug, p_tauaerl, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$acc             colh2o, coln2o, colco2, colo3, coldry, laytrop, selffac, &
    !$acc             selffrac, indself, fracs, rat_h2on2o, rat_h2on2o_1, &
    !$acc             indfor, forfac, forfrac, minorfrac, indminor) &
    !$acc       create(colco)

    !$acc parallel default(none) async(1)
    !$acc loop gang vector collapse(2)
    do lay = 1,klev
      do jc = kidia, kfdia
        colco(jc,lay) = 0._jprb
      end do
    end do
    !$acc end parallel 


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
! # 151 "ifsrrtm/rrtm_taumol13.f90"
 
      ! p = 473.420 mb (level 5)
      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)

      ! p = 1053. (level 1)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)

      ! p = 706. (level 3)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)

      ! compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! lower atmosphere loop
      !$acc wait
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(speccomb, speccomb1, speccomb_planck, speccomb_mco2, speccomb_mco, &
      !$acc   ind0, ind1, &
      !$acc   inds, indf, indm, js, js1, jpl, jmco2, jmco, fac000, fac100, fac200, fac010, &
      !$acc   fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, p4, fk0, fk1, fk2, chi_co2, ratco2, &
      !$acc   adjfac, adjcolco2, fs, specmult, specparm, fs1, specmult1, specparm1, fmco2, specmult_mco2, &
      !$acc   specparm_mco2, fmco, specmult_mco, specparm_mco, fpl, specmult_planck, specparm_planck, tau_major, &
      !$acc   tau_major1)
      do lay = 1, laytrop_min
        do jl = kidia, kfdia

          speccomb = colh2o(jl,lay) + rat_h2on2o(jl,lay)*coln2o(jl,lay)
          specparm = min(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._jprb*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - aint((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2on2o_1(jl,lay)*coln2o(jl,lay)
          specparm1 = min(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._jprb*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - aint((specmult1)))

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*coln2o(jl,lay)
          specparm_mco2 = min(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._jprb*specparm_mco2
          jmco2 = 1 + int(specmult_mco2)
          fmco2 = ((specmult_mco2) - aint((specmult_mco2)))

          !  in atmospheres where the amount of co2 is too great to be considered
          !  a minor species, adjust the column amount of co2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_jprb*chi_co2/3.55e-4_jprb
          if (ratco2 .gt. 3.0_jprb) then
            adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.68_jprb
            adjcolco2 = adjfac*3.55e-4_jprb*coldry(jl,lay)*1.e-20_jprb
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_mco = colh2o(jl,lay) + refrat_m_a3*coln2o(jl,lay)
          specparm_mco = min(colh2o(jl,lay)/speccomb_mco,oneminus)
          specmult_mco = 8._jprb*specparm_mco
          jmco = 1 + int(specmult_mco)
          fmco = ((specmult_mco) - aint((specmult_mco)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*coln2o(jl,lay)
          specparm_planck = min(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._jprb*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - aint((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(13) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(13) + js1
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
!$nec unroll(ng13)
            tau_major(1:ng13) = speccomb *    &
             (fac000 * absa(ind0,1:ng13)    + &
              fac100 * absa(ind0+1,1:ng13)  + &
              fac200 * absa(ind0+2,1:ng13)  + &
              fac010 * absa(ind0+9,1:ng13)  + &
              fac110 * absa(ind0+10,1:ng13) + &
              fac210 * absa(ind0+11,1:ng13))
          else if (specparm .gt. 0.875_jprb) then
!$nec unroll(ng13)
            tau_major(1:ng13) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng13) + &
              fac100 * absa(ind0,1:ng13)   + &
              fac000 * absa(ind0+1,1:ng13) + &
              fac210 * absa(ind0+8,1:ng13) + &
              fac110 * absa(ind0+9,1:ng13) + &
              fac010 * absa(ind0+10,1:ng13))
          else
!$nec unroll(ng13)
            tau_major(1:ng13) = speccomb *   &
             (fac000 * absa(ind0,1:ng13)   + &
              fac100 * absa(ind0+1,1:ng13) + &
              fac010 * absa(ind0+9,1:ng13) + &
              fac110 * absa(ind0+10,1:ng13))
          endif

          if (specparm1 .lt. 0.125_jprb) then
!$nec unroll(ng13)
            tau_major1(1:ng13) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng13)    + &
              fac101 * absa(ind1+1,1:ng13)  + &
              fac201 * absa(ind1+2,1:ng13)  + &
              fac011 * absa(ind1+9,1:ng13)  + &
              fac111 * absa(ind1+10,1:ng13) + &
              fac211 * absa(ind1+11,1:ng13))
          else if (specparm1 .gt. 0.875_jprb) then
!$nec unroll(ng13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng13) + &
              fac101 * absa(ind1,1:ng13)   + &
              fac001 * absa(ind1+1,1:ng13) + &
              fac211 * absa(ind1+8,1:ng13) + &
              fac111 * absa(ind1+9,1:ng13) + &
              fac011 * absa(ind1+10,1:ng13))
          else
!$nec unroll(ng13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac001 * absa(ind1,1:ng13)   + &
              fac101 * absa(ind1+1,1:ng13) + &
              fac011 * absa(ind1+9,1:ng13) + &
              fac111 * absa(ind1+10,1:ng13))
          endif

          !$acc loop seq private(taufor, tauself, co2m1, co2m2, absco2, &
          !$acc   com1, com2, absco)
!$nec unroll(ng13)
          do ig = 1, ng13
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)
            com1 = ka_mco(jmco,indm,ig) + fmco * &
                 (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
            com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                 (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
            absco = com1 + minorfrac(jl,lay) * (com2 - com1)

            taug(jl,ngs12+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(jl,lay)*absco
            fracs(jl,ngs12+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$acc end parallel

      ! upper atmosphere loop
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(indm)
      do lay = laytrop_max+1, klev
        do jl = kidia, kfdia

          indm = indminor(jl,lay)
          !$acc loop seq private(abso3)
!$nec unroll(ng13)
          do ig = 1, ng13
            abso3 = kb_mo3(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
            taug(jl,ngs12+ig,lay) = colo3(jl,lay)*abso3
            fracs(jl,ngs12+ig,lay) =  fracrefb(ig)
          enddo
        enddo

      enddo
      !$acc end parallel

      if (laytrop_max /= laytrop_min) then
        ! mixed loop
        ! lower atmosphere part
        !$acc parallel default(none) async(1)
        !$acc loop gang vector collapse(2) private(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$acc   specmult1, js1, fs1, speccomb_mco2, specparm_mco2, specmult_mco2, jmco2, fmco2, chi_co2, ratco2, adjfac, &
        !$acc   adjcolco2, speccomb_mco, specparm_mco, specmult_mco, jmco, fmco, speccomb_planck, specparm_planck, &
        !$acc   specmult_planck, jpl, fpl, ind0, ind1, inds, indf, indm, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, &
        !$acc   fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, tau_major, tau_major1)
        do lay = laytrop_min+1, laytrop_max




          ixc0 = ixc(lay)

!$nec ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)


            speccomb = colh2o(jl,lay) + rat_h2on2o(jl,lay)*coln2o(jl,lay)
            specparm = min(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._jprb*(specparm)
            js = 1 + int(specmult)
            fs = ((specmult) - aint((specmult)))

            speccomb1 = colh2o(jl,lay) + rat_h2on2o_1(jl,lay)*coln2o(jl,lay)
            specparm1 = min(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._jprb*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = ((specmult1) - aint((specmult1)))

            speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*coln2o(jl,lay)
            specparm_mco2 = min(colh2o(jl,lay)/speccomb_mco2,oneminus)
            specmult_mco2 = 8._jprb*specparm_mco2
            jmco2 = 1 + int(specmult_mco2)
            fmco2 = ((specmult_mco2) - aint((specmult_mco2)))

            !  in atmospheres where the amount of co2 is too great to be considered
            !  a minor species, adjust the column amount of co2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_jprb*chi_co2/3.55e-4_jprb
            if (ratco2 .gt. 3.0_jprb) then
              adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.68_jprb
              adjcolco2 = adjfac*3.55e-4_jprb*coldry(jl,lay)*1.e-20_jprb
            else
              adjcolco2 = colco2(jl,lay)
            endif

            speccomb_mco = colh2o(jl,lay) + refrat_m_a3*coln2o(jl,lay)
            specparm_mco = min(colh2o(jl,lay)/speccomb_mco,oneminus)
            specmult_mco = 8._jprb*specparm_mco
            jmco = 1 + int(specmult_mco)
            fmco = ((specmult_mco) - aint((specmult_mco)))

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*coln2o(jl,lay)
            specparm_planck = min(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._jprb*specparm_planck
            jpl = 1 + int(specmult_planck)
            fpl = ((specmult_planck) - aint((specmult_planck)))

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(13) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(13) + js1
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
!$nec unroll(ng13)
              tau_major(1:ng13) = speccomb *    &
              (fac000 * absa(ind0,1:ng13)    + &
                fac100 * absa(ind0+1,1:ng13)  + &
                fac200 * absa(ind0+2,1:ng13)  + &
                fac010 * absa(ind0+9,1:ng13)  + &
                fac110 * absa(ind0+10,1:ng13) + &
                fac210 * absa(ind0+11,1:ng13))
            else if (specparm .gt. 0.875_jprb) then
!$nec unroll(ng13)
              tau_major(1:ng13) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng13) + &
                fac100 * absa(ind0,1:ng13)   + &
                fac000 * absa(ind0+1,1:ng13) + &
                fac210 * absa(ind0+8,1:ng13) + &
                fac110 * absa(ind0+9,1:ng13) + &
                fac010 * absa(ind0+10,1:ng13))
            else
!$nec unroll(ng13)
              tau_major(1:ng13) = speccomb *   &
              (fac000 * absa(ind0,1:ng13)   + &
                fac100 * absa(ind0+1,1:ng13) + &
                fac010 * absa(ind0+9,1:ng13) + &
                fac110 * absa(ind0+10,1:ng13))
            endif

            if (specparm1 .lt. 0.125_jprb) then
!$nec unroll(ng13)
              tau_major1(1:ng13) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng13)    + &
                fac101 * absa(ind1+1,1:ng13)  + &
                fac201 * absa(ind1+2,1:ng13)  + &
                fac011 * absa(ind1+9,1:ng13)  + &
                fac111 * absa(ind1+10,1:ng13) + &
                fac211 * absa(ind1+11,1:ng13))
            else if (specparm1 .gt. 0.875_jprb) then
!$nec unroll(ng13)
              tau_major1(1:ng13) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng13) + &
                fac101 * absa(ind1,1:ng13)   + &
                fac001 * absa(ind1+1,1:ng13) + &
                fac211 * absa(ind1+8,1:ng13) + &
                fac111 * absa(ind1+9,1:ng13) + &
                fac011 * absa(ind1+10,1:ng13))
            else
!$nec unroll(ng13)
              tau_major1(1:ng13) = speccomb1 * &
              (fac001 * absa(ind1,1:ng13)   + &
                fac101 * absa(ind1+1,1:ng13) + &
                fac011 * absa(ind1+9,1:ng13) + &
                fac111 * absa(ind1+10,1:ng13))
            endif

!$nec unroll(ng13)
            !$acc loop seq private(tauself, taufor, co2m1, absco2, com1, com2, absco)
            do ig = 1, ng13
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                  (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
              co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                  (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
              absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)
              com1 = ka_mco(jmco,indm,ig) + fmco * &
                  (ka_mco(jmco+1,indm,ig) - ka_mco(jmco,indm,ig))
              com2 = ka_mco(jmco,indm+1,ig) + fmco * &
                  (ka_mco(jmco+1,indm+1,ig) - ka_mco(jmco,indm+1,ig))
              absco = com1 + minorfrac(jl,lay) * (com2 - com1)

              taug(jl,ngs12+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor &
                  + adjcolco2*absco2 &
                  + colco(jl,lay)*absco
              fracs(jl,ngs12+ig,lay) = fracrefa(ig,jpl) + fpl * &
                  (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
            enddo



          enddo

          ! upper atmosphere part
          ixc0 = kfdia - kidia + 1 - ixc0
!$nec ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)


            indm = indminor(jl,lay)
!$nec unroll(ng13)
            !$acc loop seq private(abso3)
            do ig = 1, ng13
              abso3 = kb_mo3(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
              taug(jl,ngs12+ig,lay) = colo3(jl,lay)*abso3
              fracs(jl,ngs12+ig,lay) =  fracrefb(ig)
            enddo



          enddo

        enddo
        !$acc end parallel

      endif

      !$acc wait
      !$acc end data

end subroutine rrtm_taumol13
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

