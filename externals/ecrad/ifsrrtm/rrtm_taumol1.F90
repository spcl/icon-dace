! # 1 "ifsrrtm/rrtm_taumol1.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/rrtm_taumol1.f90"
! this file has been modified for the use in icon

subroutine rrtm_taumol1 (kidia,kfdia,klev,taug,pavel,&
 & p_tauaerl,fac00,fac01,fac10,fac11,forfac,forfrac,indfor,jp,jt,jt1,&
 & colh2o,laytrop,selffac,selffrac,indself,fracs,minorfrac,indminor,&
 & scaleminorn2,colbrd)  

!******************************************************************************
!                                                                             *
!                  optical depths developed for the                           *
!                                                                             *
!                rapid radiative transfer model (rrtm)                        *
!                                                                             *
!            atmospheric and environmental research, inc.                     *
!                        840 memorial drive                                   *
!                        cambridge, ma 02139                                  *
!                                                                             *
!                           eli j. mlawer                                     *
!                         steven j. taubman                                   *
!                         shepard a. clough                                   *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        the authors wish to acknowledge the contributions of the             *
!        following people:  patrick d. brown, michael j. iacono,              *
!        ronald e. farren, luke chen, robert bergstrom.                       *
!                                                                             *
!******************************************************************************
!     taumol                                                                  *
!                                                                             *
!     this file contains the subroutines taugbn (where n goes from            *
!     1 to 16).  taugbn calculates the optical depths and planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  output:  optical depths (unitless)                                         *
!           fractions needed to compute planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     common /taugcom/  taug(mxlay,mg)                                        *
!     common /plankg/   fracs(mxlay,mg)                                       *
!                                                                             *
!  input                                                                      *
!                                                                             *
!     common /features/ ng(nbands),nspa(nbands),nspb(nbands)                  *
!     common /precise/  oneminus                                              *
!     common /profile/  klev,pavel(mxlay),tavel(mxlay),                    *
!    &                  pz(0:mxlay),tz(0:mxlay),tbound                        *
!     common /profdata/ laytrop,layswtch,laylow,                              *
!    &                  colh2o(mxlay),colco2(mxlay),                          *
!    &                  colo3(mxlay),coln2o(mxlay),colch4(mxlay),             *
!    &                  colo2(mxlay),co2mult(mxlay)                           *
!     common /intfac/   fac00(mxlay),fac01(mxlay),                            *
!    &                  fac10(mxlay),fac11(mxlay)                             *
!     common /intind/   jp(mxlay),jt(kidia:kfdia,mxlay),jt1(kidia:kfdia,mxlay)                        *
!     common /self/     selffac(mxlay), selffrac(mxlay), indself(kidia:kfdia,mxlay)       *
!                                                                             *
!     description:                                                            *
!     ng(iband) - number of g-values in band iband                            *
!     nspa(iband) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band iband per            *
!                   pressure level and temperature.  each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     nspb(iband) - same for upper atmosphere                                 *
!     oneminus - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     pavel - layer pressures (mb)                                            *
!     tavel - layer temperatures (degrees k)                                  *
!     pz - level pressures (mb)                                               *
!     tz - level temperatures (degrees k)                                     *
!     laytrop - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     colh2o, colco2, colo3, coln2o, colch4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     co2mult - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average co2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     facij(lay) - for layer lay, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  a value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     jp - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     jt, jt1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels jp and jp+1, respectively)                             *
!     selffac - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296k and        *
!               1013 mb)                                                      *
!     selffrac - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     indself - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  data input                                                                 *
!     common /kn/ ka(nspa(n),5,13,mg), kb(nspb(n),5,13:59,mg), selfref(10,mg) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     description:                                                            *
!     ka - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     kb - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     selfref - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below laytrop)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     dimension absa(65*nspa(n),mg), absb(235*nspb(n),mg)                     *
!     equivalence (ka,absa),(kb,absb)                                         *
!                                                                             *
!******************************************************************************

!     band 1:  10-250 cm-1 (low - h2o; high - h2o)
 
!     author.
!     -------
!      jjmorcrette, ecmwf, from
!      eli j. mlawer, atmospheric & environmental research.
!      (revised by michael j. iacono, atmospheric & environmental research.)

!     modifications.
!     --------------
!      d salmond   2000-05-15 speed-up
!      jjmorcrette 2000-05-17 speed-up
!      m.hamrud      01-oct-2003 cy28 cleaning
!      nec           25-oct-2007 optimisations
!      jjmorcrette 20110613 flexible number of g-points
!      abozzo 200130517 updated to rrtmg_lw_v4.85:
!*********
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1: 
!           10-250 cm-1 (low - h2o; high - h2o)
! ---------------------------------------------------------------------------

use parkind1  ,only : jpim     ,jprb
use ecradhook   ,only : lhook,   dr_hook

use parrrtm  , only : jpband
use yoerrtm  , only : jpgpt  ,ng1
use yoerrtwn , only : nspa   ,nspb
use yoerrta1 , only : absa   ,absb   ,fracrefa, fracrefb,&
 & forref   ,selfref,  ka_mn2, kb_mn2   

implicit none

integer(kind=jpim),intent(in)    :: kidia
integer(kind=jpim),intent(in)    :: kfdia
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(in)    :: pavel(kidia:kfdia,klev) ! layer pressures (hpa)
real(kind=jprb)   ,intent(inout) :: taug(kidia:kfdia,jpgpt,klev) 
real(kind=jprb)   ,intent(in)    :: p_tauaerl(kidia:kfdia,klev,jpband) 
real(kind=jprb)   ,intent(in)    :: fac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: fac11(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: forfrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jp(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jt(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: jt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: colh2o(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: laytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: selffrac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: minorfrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(inout) :: fracs(kidia:kfdia,jpgpt,klev) 

integer(kind=jpim),intent(in)    :: indfor(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: indminor(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: scaleminorn2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: colbrd(kidia:kfdia,klev)         
! ---------------------------------------------------------------------------

integer(kind=jpim) :: ind0,ind1,inds
integer(kind=jpim) :: indf,indm

integer(kind=jpim) :: ig, lay
real(kind=jprb) :: taufor,tauself,corradj,pp,scalen2, taun2
    !     local integer arrays
    integer(kind=jpim) :: laytrop_min, laytrop_max
    integer(kind=jpim) :: ixc(klev), ixlow(kfdia,klev), ixhigh(kfdia,klev)
    integer(kind=jpim) :: ich, icl, ixc0, ixp, jc, jl

    !$acc data present(pavel, taug, p_tauaerl, fac00, fac01, fac10, fac11, &
    !$acc             forfac, forfrac, jp, jt, jt1, colh2o, laytrop, selffac, &
    !$acc             selffrac, minorfrac, indself, fracs, &
    !$acc             indfor, indminor, scaleminorn2, colbrd)


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
! # 232 "ifsrrtm/rrtm_taumol1.f90"


! minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

!     compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  below laytrop, the water vapor self-continuum and
!     foreign continuum is interpolated (in temperature) separately.

      ! lower atmosphere loop
      !$acc wait
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, inds, indf, indm, pp, corradj)
      do lay = 1, laytrop_min
        do jl = kidia, kfdia

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(1) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(1) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1.0_jprb
          if (pp .lt. 250._jprb) then
            corradj = 1._jprb - 0.15_jprb * (250._jprb-pp) / 154.4_jprb
          endif

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
          !$acc loop seq private(tauself, taufor, taun2)
!$nec unroll(ng1)
          do ig = 1, ng1
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) -  forref(indf,ig)))
            taun2 = scalen2*(ka_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
            taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor + taun2)
            fracs(jl,ig,lay) = fracrefa(ig)
          enddo
        enddo
      enddo
      !$acc end parallel

      ! upper atmosphere loop
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2) private(ind0, ind1, indf, indm, pp, corradj, scalen2)
      do lay = laytrop_max+1, klev
        do jl = kidia, kfdia

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(1) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(1) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1._jprb - 0.15_jprb * (pp / 95.6_jprb)

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
          !$acc loop seq private(taufor, taun2)
!$nec unroll(ng1)
          do ig = 1, ng1
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taun2 = scalen2*(kb_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
            taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor + taun2)
            fracs(jl,ig,lay) = fracrefb(ig)
          enddo
        enddo
      enddo
      !$acc end parallel

      if (laytrop_max /= laytrop_min) then
        ! mixed loop
        ! lower atmosphere part
        !$acc parallel default(none) async(1)
        !$acc loop gang vector collapse(2) private(ind0, ind1, inds, indf, indm, pp, corradj, scalen2)
        do lay = laytrop_min+1, laytrop_max




          ixc0 = ixc(lay)
!$nec ivdep
          do ixp = 1, ixc0
            jl = ixlow(ixp,lay)


            ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(1) + 1
            ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(1) + 1
            inds = indself(jl,lay)
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
            pp = pavel(jl,lay)
            corradj =  1.0_jprb
            if (pp .lt. 250._jprb) then
              corradj = 1._jprb - 0.15_jprb * (250._jprb-pp) / 154.4_jprb
            endif

            scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$nec unroll(ng1)
            !$acc loop seq private(tauself, taufor, taun2)
            do ig = 1, ng1
              tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                  (selfref(inds+1,ig) - selfref(inds,ig)))
              taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                  (forref(indf+1,ig) -  forref(indf,ig)))
              taun2 = scalen2*(ka_mn2(indm,ig) + &
                  minorfrac(jl,lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
              taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                  (fac00(jl,lay) * absa(ind0,ig) + &
                  fac10(jl,lay) * absa(ind0+1,ig) + &
                  fac01(jl,lay) * absa(ind1,ig) + &
                  fac11(jl,lay) * absa(ind1+1,ig)) &
                  + tauself + taufor + taun2)
              fracs(jl,ig,lay) = fracrefa(ig)
            enddo



          enddo



          ! upper atmosphere part
          ixc0 = kfdia - kidia + 1 - ixc0
!$nec ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)


            ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(1) + 1
            ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(1) + 1
            indf = indfor(jl,lay)
            indm = indminor(jl,lay)
            pp = pavel(jl,lay)
            corradj =  1._jprb - 0.15_jprb * (pp / 95.6_jprb)

            scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$nec unroll(ng1)
            !$acc loop seq private(taufor, taun2)
            do ig = 1, ng1
              taufor = forfac(jl,lay) * (forref(indf,ig) + &
                  forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
              taun2 = scalen2*(kb_mn2(indm,ig) + &
                  minorfrac(jl,lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
              taug(jl,ig,lay) = corradj * (colh2o(jl,lay) * &
                  (fac00(jl,lay) * absb(ind0,ig) + &
                  fac10(jl,lay) * absb(ind0+1,ig) + &
                  fac01(jl,lay) * absb(ind1,ig) + &
                  fac11(jl,lay) * absb(ind1+1,ig)) &
                  + taufor + taun2)
              fracs(jl,ig,lay) = fracrefb(ig)
            enddo



          enddo

        enddo
        !$acc end parallel

      endif

      !$acc end data

end subroutine rrtm_taumol1
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

