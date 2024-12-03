! # 1 "ifsrrtm/rrtm_taumol6.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_taumol6.f90"
! this file has been modified for the use in icon

!----------------------------------------------------------------------------
subroutine rrtm_taumol6 (kidia,kfdia,klev,taug,wx,&
 & p_tauaerl,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,colco2,coldry,laytrop,selffac,selffrac,indself,fracs,minorfrac,indminor)  

!     band 6:  820-980 cm-1 (low - h2o; high - nothing)

!     author.
!     -------
!      jjmorcrette, ecmwf

!     modifications.
!     --------------
!      m.hamrud      01-oct-2003 cy28 cleaning
!      nec           25-oct-2007 optimisations
!      jjmorcrette 20110613 flexible number of g-points
!      abozzo 201306 updated to rrtmg v4.85
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
! ---------------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband ,jpxsec
use yoerrtm  , only : jpgpt  ,ng6   ,ngs5
use yoerrtwn , only : nspa   
use yoerrta6 , only : absa   ,ka_mco2 ,cfc11adj , cfc12  ,&
 & fracrefa,selfref,forref  
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
real(kind=jprb)   ,intent(in)    :: colco2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: coldry(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: laytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: fracs(kidia:kfdia,jpgpt,klev) 

integer(kind=jpim),intent(in)   :: indfor(kidia:kfdia,klev)
real(kind=jprb)   ,intent(in)   :: forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: forfrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)   :: minorfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in)   :: indminor(kidia:kfdia,klev)

! ---------------------------------------------------------------------------

integer(kind=jpim) :: ind0,ind1,inds,indf,indm

integer(kind=jpim) :: ig, lay

real(kind=jprb) :: adjfac,adjcolco2,ratco2,chi_co2
real(kind=jprb) :: taufor,tauself,absco2
    !     local integer arrays
    integer(kind=jpim) :: laytrop_min, laytrop_max
    integer(kind=jpim) :: ixc(klev), ixlow(kfdia,klev), ixhigh(kfdia,klev)
    integer(kind=jpim) :: ich, icl, ixc0, ixp, jc, jl

    !$acc data present(taug, wx, p_tauaerl, fac00, fac01, fac10, fac11, jp, jt,  &
    !$acc             jt1, colh2o, colco2, coldry, laytrop, selffac, selffrac , &
    !$acc             indself, fracs, indfor, forfac, forfrac, minorfrac, &
    !$acc             indminor)


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
! # 116 "ifsrrtm/rrtm_taumol6.f90"

! minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12


!     compute the optical depth by interpolating in ln(pressure) and
!     temperature. the water vapor self- and foreign- continuum is interpolated
!     (in temperature) separately.  

      ! lower atmosphere loop
      !$acc wait
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, inds, indf, indm, adjfac, adjcolco2, ratco2,chi_co2)
      do lay = 1, laytrop_min
        do jl = kidia, kfdia

          ! in atmospheres where the amount of co2 is too great to be considered
          ! a minor species, adjust the column amount of co2 by an empirical factor
          ! to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_jprb*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_jprb) then
            adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.77_jprb
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(6) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(6) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          !$acc loop seq private(taufor, tauself, absco2)
!$nec unroll(ng6)
          do ig = 1, ng6
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(jl,ngs5+ig,lay) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,ngs5+ig,lay) = fracrefa(ig)
          enddo
        enddo
      enddo
      !$acc end parallel

      ! upper atmosphere loop
      ! nothing important goes on above laytrop in this band.
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(3)
      do ig = 1, ng6
        do lay = laytrop_max+1, klev
          do jl = kidia, kfdia
            taug(jl,ngs5+ig,lay) = 0.0_jprb &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,ngs5+ig,lay) = fracrefa(ig)
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


            ! in atmospheres where the amount of co2 is too great to be considered
            ! a minor species, adjust the column amount of co2 by an empirical factor
            ! to obtain the proper contribution.
            chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
            ratco2 = 1.e20_jprb*chi_co2/chi_mls(2,jp(jl,lay)+1)
            if (ratco2 .gt. 3.0_jprb) then
              adjfac = 2.0_jprb+(ratco2-2.0_jprb)**0.77_jprb
              adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_jprb
            else
              adjcolco2 = colco2(jl,lay)
            endif

            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(6) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(6) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
!$nec unroll(ng6)
            !$acc loop seq private(tauself, taufor, absco2)
            do ig = 1, ng6
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                (forref(indf+1,ig) - forref(indf,ig)))
              absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                  (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
              taug(jl,ngs5+ig,lay) = colh2o(jl,lay) * &
                  (fac00(jl,lay) * absa(ind0,ig) + &
                  fac10(jl,lay) * absa(ind0+1,ig) + &
                  fac01(jl,lay) * absa(ind1,ig) +  &
                  fac11(jl,lay) * absa(ind1+1,ig))  &
                  + tauself + taufor &
                  + adjcolco2 * absco2 &
                  + wx(jl,2,lay) * cfc11adj(ig) &
                  + wx(jl,3,lay) * cfc12(ig)
              fracs(jl,ngs5+ig,lay) = fracrefa(ig)
            enddo



          enddo

          ! upper atmosphere part
          ! nothing important goes on above laytrop in this band.
          ixc0 = kfdia - kidia + 1 - ixc0


          !$acc loop seq
          do ig = 1, ng6

!$nec ivdep
            do ixp = 1, ixc0
              jl = ixhigh(ixp,lay)

              taug(jl,ngs5+ig,lay) = 0.0_jprb &
                  + wx(jl,2,lay) * cfc11adj(ig) &
                  + wx(jl,3,lay) * cfc12(ig)
              fracs(jl,ngs5+ig,lay) = fracrefa(ig)
            enddo



          enddo

        enddo
        !$acc end parallel

      endif

      !$acc end data

end subroutine rrtm_taumol6
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

