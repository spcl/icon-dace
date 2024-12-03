! # 1 "ifsrrtm/rrtm_taumol8.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_taumol8.f90"
! this file has been modified for the use in icon

!*******************************************************************************
subroutine rrtm_taumol8 (kidia,kfdia,klev,taug,wx,&
 & p_tauaerl,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,colo3,coln2o,colco2,coldry,laytrop,selffac,selffrac,indself,fracs, &
 & minorfrac,indminor)  

!     band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)

!     author.
!     -------
!      jjmorcrette, ecmwf

!     modifications.
!     --------------
!      m.hamrud      01-oct-2003 cy28 cleaning
!      nec           25-oct-2007 optimisations
!      jjmorcrette 20110613 flexible number of g-points
!      abozzo 201306 updated to rrtmg v4.85
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
! ---------------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband ,jpxsec
use yoerrtm  , only : jpgpt  ,ng8   ,ngs7
use yoerrtwn , only : nspa   ,nspb
use yoerrta8 , only : absa   ,absb   ,fracrefa, fracrefb,selfref,ka_mco2 ,kb_mco2  ,&
 & ka_mn2o , kb_mn2o,ka_mo3,cfc12  ,cfc22adj,forref  
use yoerrtrf, only : chi_mls

implicit none

integer(kind=jpim),intent(in)    :: kidia
integer(kind=jpim),intent(in)    :: kfdia
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(inout) :: taug(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(in)    :: wx(kidia:kfdia,jpxsec,klev) ! amount of trace gases
real(kind=jprb)   ,intent(in)    :: p_tauaerl(kidia:kfdia,klev,jpband) 
real(kind=jprb)   ,intent(in)    :: fac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac11(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jp(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jt(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: colh2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: colo3(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: coln2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: colco2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: coldry(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: laytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: fracs(kidia:kfdia,jpgpt,klev) 

integer(kind=jpim),intent(in)   :: indfor(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: forfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: minorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indminor(kidia:kfdia,klev)

! ---------------------------------------------------------------------------

integer(kind=jpim) :: ind0,ind1,inds,indf,indm

integer(kind=jpim) :: ig, lay

real(kind=jprb) :: chi_co2, ratco2, adjfac, adjcolco2
real(kind=jprb) :: taufor,tauself, abso3, absco2, absn2o
    !     local integer arrays
    integer(kind=jpim) :: laytrop_min, laytrop_max
    integer(kind=jpim) :: ixc(klev), ixlow(kfdia,klev), ixhigh(kfdia,klev)
    integer(kind=jpim) :: ich, icl, ixc0, ixp, jc, jl

    !$acc data present(taug, wx, p_tauaerl, fac00, fac01, fac10, fac11, jp, jt, &
    !$acc             jt1, colh2o, colo3, coln2o, colco2, coldry, laytrop, &
    !$acc             selffac, selffrac, indself, fracs, indfor, forfrac, forfac, &
    !$acc             minorfrac, indminor)


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
! # 119 "ifsrrtm/rrtm_taumol8.f90"

! minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

! compute the optical depth by interpolating in ln(pressure) and 
! temperature, and appropriate species.  below laytrop, the water vapor 
! self-continuum and foreign continuum is interpolated (in temperature) 
! separately.

      ! lower atmosphere loop
      !$acc wait
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, inds, indf, indm, chi_co2, ratco2, adjfac, adjcolco2)
      do lay = 1, laytrop_min
        do jl = kidia, kfdia

          !  in atmospheres where the amount of co2 is too great to be considered
          !  a minor species, adjust the column amount of co2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_jprb*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_jprb) then
            adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.65_jprb
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(8) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(8) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          !$acc loop seq private(taufor, tauself, abso3, absco2, absn2o)
!$nec unroll(ng8)
          do ig = 1, ng8
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            abso3 =  (ka_mo3(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
            absn2o =  (ka_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
            taug(jl,ngs7+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colo3(jl,lay) * abso3 &
                 + coln2o(jl,lay) * absn2o &
                 + wx(jl,3,lay) * cfc12(ig) &
                 + wx(jl,4,lay) * cfc22adj(ig)
            fracs(jl,ngs7+ig,lay) = fracrefa(ig)
          enddo
        enddo

      enddo
      !$acc end parallel

      ! upper atmosphere loop
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, indm, chi_co2, ratco2, adjfac, adjcolco2)
      do lay = laytrop_max+1, klev
        do jl = kidia, kfdia

          !  in atmospheres where the amount of co2 is too great to be considered
          !  a minor species, adjust the column amount of co2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/coldry(jl,lay)
          ratco2 = 1.e20_jprb*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_jprb) then
            adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.65_jprb
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1) * coldry(jl,lay)*1.e-20_jprb
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(8) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(8) + 1
          indm = indminor(jl,lay)
          !$acc loop seq private(absco2, absn2o)
!$nec unroll(ng8)
          do ig = 1, ng8
            absco2 =  (kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
            absn2o =  (kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
            taug(jl,ngs7+ig,lay) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(jl,lay)*absn2o &
                 + wx(jl,3,lay) * cfc12(ig) &
                 + wx(jl,4,lay) * cfc22adj(ig)
            fracs(jl,ngs7+ig,lay) = fracrefb(ig)
          enddo
        enddo

      enddo
      !$acc end parallel

      if (laytrop_max /= laytrop_min) then
        ! mixed loop
        ! lower atmosphere part
        !$acc parallel default(none) async(1)
        !$acc loop gang vector collapse(2) private(chi_co2, ratco2, adjfac, adjcolco2, ind0, ind1, inds, indf, indm)
        do lay = laytrop_min+1, laytrop_max




          ixc0 = ixc(lay)
!$nec ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)


            !  in atmospheres where the amount of co2 is too great to be considered
            !  a minor species, adjust the column amount of co2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_jprb*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_jprb) then
              adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.65_jprb
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
            else
              adjcolco2 = colco2(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(8) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(8) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
!$nec unroll(ng8)
            !$acc loop seq private(tauself, taufor, absco2, abso3, absn2o)
            do ig = 1, ng8
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) - forref(indf,ig)))
              absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
              abso3 =  (ka_mo3(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mo3(indm+1,ig) - ka_mo3(indm,ig)))
              absn2o =  (ka_mn2o(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mn2o(indm+1,ig) - ka_mn2o(indm,ig)))
              taug(jl,ngs7+ig,lay) = colh2o(jl,lay) * &
                  (fac00(jl,lay) * absa(ind0,ig) + &
                  fac10(jl,lay) * absa(ind0+1,ig) + &
                  fac01(jl,lay) * absa(ind1,ig) +  &
                  fac11(jl,lay) * absa(ind1+1,ig)) &
                  + tauself + taufor &
                  + adjcolco2*absco2 &
                  + colo3(jl,lay) * abso3 &
                  + coln2o(jl,lay) * absn2o &
                  + wx(jl,3,lay) * cfc12(ig) &
                  + wx(jl,4,lay) * cfc22adj(ig)
              fracs(jl,ngs7+ig,lay) = fracrefa(ig)
            enddo



          enddo

          ! upper atmosphere loop
          ixc0 = kfdia - kidia + 1 - ixc0
!$nec ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)


            !  in atmospheres where the amount of co2 is too great to be considered
            !  a minor species, adjust the column amount of co2 by an empirical factor
            !  to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/coldry(jl,lay)
            ratco2 = 1.e20_jprb*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_jprb) then
              adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.65_jprb
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1) * coldry(jl,lay)*1.e-20_jprb
            else
              adjcolco2 = colco2(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(8) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(8) + 1
            indm = indminor(jl,lay)
!$nec unroll(ng8)
            !$acc loop seq private(absco2, absn2o)
            do ig = 1, ng8
              absco2 =  (kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
              absn2o =  (kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                  (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
              taug(jl,ngs7+ig,lay) = colo3(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) + &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + adjcolco2*absco2 &
                  + coln2o(jl,lay)*absn2o &
                  + wx(jl,3,lay) * cfc12(ig) &
                  + wx(jl,4,lay) * cfc22adj(ig)
              fracs(jl,ngs7+ig,lay) = fracrefb(ig)
            enddo



          enddo

        enddo
        !$acc end parallel

      endif

      !$acc end data

end subroutine rrtm_taumol8
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

