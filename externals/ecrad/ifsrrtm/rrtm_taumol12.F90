! # 1 "ifsrrtm/rrtm_taumol12.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_taumol12.f90"
! this file has been modified for the use in icon

!----------------------------------------------------------------------------
subroutine rrtm_taumol12 (kidia,kfdia,klev,taug,&
 & p_tauaerl,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,oneminus,&
 & colh2o,colco2,laytrop,selffac,selffrac,indself,fracs, &  
 & rat_h2oco2, rat_h2oco2_1)

!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)

!     author.
!     -------
!      jjmorcrette, ecmwf

!     modifications.
!     --------------
!      m.hamrud      01-oct-2003 cy28 cleaning
!      nec           25-oct-2007 optimisations
!      jjmorcrette 20110613 flexible number of g-points
!      abozzo updated to rrtmg v4.85
! ---------------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband
use yoerrtm  , only : jpgpt  ,ng12 ,ngs11
use yoerrtwn , only :      nspa    
use yoerrta12, only : absa   ,fracrefa,selfref,forref
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
real(kind=jprb)   ,intent(in)    :: colco2(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: laytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: fracs(kidia:kfdia,jpgpt,klev) 

real(kind=jprb)   ,intent(in)   :: rat_h2oco2(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: rat_h2oco2_1(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indfor(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: forfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: forfac(kidia:kfdia,klev) 


! ---------------------------------------------------------------------------

real(kind=jprb) :: speccomb,speccomb1,speccomb_planck
integer(kind=jpim) :: ind0,ind1,inds,indf

integer(kind=jpim) :: ig, js, lay,js1, jpl

! real(kind=jprb) :: fac000, fac001, fac010, fac011, fac100, fac101,&
!  & fac110, fac111
real(kind=jprb) :: fs, specmult, specparm,  &
& fs1, specmult1, specparm1, &
& fpl, specmult_planck, specparm_planck

real(kind=jprb) ::  fac000, fac100, fac200,&
 & fac010, fac110, fac210, &
 & fac001, fac101, fac201, &
 & fac011, fac111, fac211
real(kind=jprb) :: p, p4, fk0, fk1, fk2
real(kind=jprb) :: taufor,tauself,tau_major(ng12),tau_major1(ng12)
real(kind=jprb) :: refrat_planck_a

    !     local integer arrays
    integer(kind=jpim) :: laytrop_min, laytrop_max
    integer(kind=jpim) :: ixc(klev), ixlow(kfdia,klev), ixhigh(kfdia,klev)
    integer(kind=jpim) :: ich, icl, ixc0, ixp, jc, jl



    !$acc data present(taug, p_tauaerl, fac00, fac01, fac10, fac11, jp, jt, jt1, &
    !$acc             colh2o, colco2, laytrop, selffac, selffrac, indself, fracs, &
    !$acc             rat_h2oco2, rat_h2oco2_1, indfor, forfrac, forfac)


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
! # 128 "ifsrrtm/rrtm_taumol12.f90"

!  ----------------------------------------------------------

      ! p =   174.164 mb
      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)

      ! compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  below laytrop, the water
      ! vapor self-continuum adn foreign continuum is interpolated
      ! (in temperature) separately.

      ! lower atmosphere loop
      !$acc wait
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(speccomb, speccomb1, speccomb_planck, ind0,ind1,inds,indf, js, js1, &
      !$acc   jpl, fs, specmult, specparm, fs1, specmult1, specparm1, fpl, specmult_planck, specparm_planck, fac000, &
      !$acc   fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, fac011, fac111, fac211, p, p4, fk0, fk1, &
      !$acc   fk2, tau_major, tau_major1)
      do lay = 1, laytrop_min
        do jl = kidia, kfdia

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = min(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._jprb*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - aint((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = min(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._jprb*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - aint((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = min(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._jprb*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - aint((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(12) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(12) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

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
!$nec unroll(ng12)
            tau_major(1:ng12) = speccomb *    &
             (fac000 * absa(ind0,1:ng12)    + &
              fac100 * absa(ind0+1,1:ng12)  + &
              fac200 * absa(ind0+2,1:ng12)  + &
              fac010 * absa(ind0+9,1:ng12)  + &
              fac110 * absa(ind0+10,1:ng12) + &
              fac210 * absa(ind0+11,1:ng12))
          else if (specparm .gt. 0.875_jprb) then
!$nec unroll(ng12)
            tau_major(1:ng12) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng12) + &
              fac100 * absa(ind0,1:ng12)   + &
              fac000 * absa(ind0+1,1:ng12) + &
              fac210 * absa(ind0+8,1:ng12) + &
              fac110 * absa(ind0+9,1:ng12) + &
              fac010 * absa(ind0+10,1:ng12))
          else
!$nec unroll(ng12)
            tau_major(1:ng12) = speccomb *   &
             (fac000 * absa(ind0,1:ng12)   + &
              fac100 * absa(ind0+1,1:ng12) + &
              fac010 * absa(ind0+9,1:ng12) + &
              fac110 * absa(ind0+10,1:ng12))
          endif

          if (specparm1 .lt. 0.125_jprb) then
!$nec unroll(ng12)
            tau_major1(1:ng12) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng12)    + &
              fac101 * absa(ind1+1,1:ng12)  + &
              fac201 * absa(ind1+2,1:ng12)  + &
              fac011 * absa(ind1+9,1:ng12)  + &
              fac111 * absa(ind1+10,1:ng12) + &
              fac211 * absa(ind1+11,1:ng12))
          else if (specparm1 .gt. 0.875_jprb) then
!$nec unroll(ng12)
            tau_major1(1:ng12) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng12) + &
              fac101 * absa(ind1,1:ng12)   + &
              fac001 * absa(ind1+1,1:ng12) + &
              fac211 * absa(ind1+8,1:ng12) + &
              fac111 * absa(ind1+9,1:ng12) + &
              fac011 * absa(ind1+10,1:ng12))
          else
!$nec unroll(ng12)
            tau_major1(1:ng12) = speccomb1 * &
             (fac001 * absa(ind1,1:ng12)   + &
              fac101 * absa(ind1+1,1:ng12) + &
              fac011 * absa(ind1+9,1:ng12) + &
              fac111 * absa(ind1+10,1:ng12))
          endif

          !$acc loop seq private(taufor, tauself)
!$nec unroll(ng12)
          do ig = 1, ng12
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,ngs11+ig,lay) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,ngs11+ig,lay) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo
      !$acc end parallel

      ! upper atmosphere loop
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(3)
      do ig = 1, ng12
        do lay = laytrop_max+1, klev
          do jl = kidia, kfdia
            taug(jl,ngs11+ig,lay) = 0.0_jprb
            fracs(jl,ngs11+ig,lay) = 0.0_jprb
          enddo
        enddo

      enddo
      !$acc end parallel

      if (laytrop_max /= laytrop_min) then
        ! mixed loop
        ! lower atmosphere part
        !$acc parallel default(none) async(1)
        !$acc loop gang vector collapse(2) private(speccomb, specparm, specmult, js, fs, speccomb1, specparm1, &
        !$acc   specmult1, js1, fs1, speccomb_planck, specparm_planck, specmult_planck, jpl, fpl, ind0, ind1, inds, &
        !$acc   indf, p, p4, fk0, fk1, fk2, fac000, fac100, fac200, fac010, fac110, fac210, fac001, fac101, fac201, &
        !$acc   fac011, fac111, fac211, tau_major, tau_major1)
        do lay = laytrop_min+1, laytrop_max




          ixc0 = ixc(lay)

!$nec ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)


            speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
            specparm = min(colh2o(jl,lay)/speccomb,oneminus)
            specmult = 8._jprb*(specparm)
            js = 1 + int(specmult)
            fs = ((specmult) - aint((specmult)))

            speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
            specparm1 = min(colh2o(jl,lay)/speccomb1,oneminus)
            specmult1 = 8._jprb*(specparm1)
            js1 = 1 + int(specmult1)
            fs1 = ((specmult1) - aint((specmult1)))

            speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
            specparm_planck = min(colh2o(jl,lay)/speccomb_planck,oneminus)
            specmult_planck = 8._jprb*specparm_planck
            jpl = 1 + int(specmult_planck)
            fpl = ((specmult_planck) - aint((specmult_planck)))

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(12) + js
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(12) + js1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)

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
!$nec unroll(ng12)
              tau_major(1:ng12) = speccomb *    &
              (fac000 * absa(ind0,1:ng12)    + &
                fac100 * absa(ind0+1,1:ng12)  + &
                fac200 * absa(ind0+2,1:ng12)  + &
                fac010 * absa(ind0+9,1:ng12)  + &
                fac110 * absa(ind0+10,1:ng12) + &
                fac210 * absa(ind0+11,1:ng12))
            else if (specparm .gt. 0.875_jprb) then
!$nec unroll(ng12)
              tau_major(1:ng12) = speccomb *   &
              (fac200 * absa(ind0-1,1:ng12) + &
                fac100 * absa(ind0,1:ng12)   + &
                fac000 * absa(ind0+1,1:ng12) + &
                fac210 * absa(ind0+8,1:ng12) + &
                fac110 * absa(ind0+9,1:ng12) + &
                fac010 * absa(ind0+10,1:ng12))
            else
!$nec unroll(ng12)
              tau_major(1:ng12) = speccomb *   &
              (fac000 * absa(ind0,1:ng12)   + &
                fac100 * absa(ind0+1,1:ng12) + &
                fac010 * absa(ind0+9,1:ng12) + &
                fac110 * absa(ind0+10,1:ng12))
            endif

            if (specparm1 .lt. 0.125_jprb) then
!$nec unroll(ng12)
              tau_major1(1:ng12) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng12)    + &
                fac101 * absa(ind1+1,1:ng12)  + &
                fac201 * absa(ind1+2,1:ng12)  + &
                fac011 * absa(ind1+9,1:ng12)  + &
                fac111 * absa(ind1+10,1:ng12) + &
                fac211 * absa(ind1+11,1:ng12))
            else if (specparm1 .gt. 0.875_jprb) then
!$nec unroll(ng12)
              tau_major1(1:ng12) = speccomb1 * &
              (fac201 * absa(ind1-1,1:ng12) + &
                fac101 * absa(ind1,1:ng12)   + &
                fac001 * absa(ind1+1,1:ng12) + &
                fac211 * absa(ind1+8,1:ng12) + &
                fac111 * absa(ind1+9,1:ng12) + &
                fac011 * absa(ind1+10,1:ng12))
            else
!$nec unroll(ng12)
              tau_major1(1:ng12) = speccomb1 * &
              (fac001 * absa(ind1,1:ng12)   + &
                fac101 * absa(ind1+1,1:ng12) + &
                fac011 * absa(ind1+9,1:ng12) + &
                fac111 * absa(ind1+10,1:ng12))
            endif

!$nec unroll(ng12)
            !$acc loop seq private(tauself, taufor)
            do ig = 1, ng12
              tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))

              taug(jl,ngs11+ig,lay) = tau_major(ig) + tau_major1(ig) &
                  + tauself + taufor
              fracs(jl,ngs11+ig,lay) = fracrefa(ig,jpl) + fpl * &
                  (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
            enddo



          enddo

          ! upper atmosphere part
          ixc0 = kfdia - kidia + 1 - ixc0


          !$acc loop seq
          do ig = 1, ng12

!$nec ivdep
            do ixp = 1, ixc0
              jl = ixhigh(ixp,lay)


              taug(jl,ngs11+ig,lay) = 0.0_jprb
              fracs(jl,ngs11+ig,lay) = 0.0_jprb
            enddo



          enddo

        enddo
        !$acc end parallel

      endif

      !$acc end data

end subroutine rrtm_taumol12
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

