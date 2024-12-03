! # 1 "ifsrrtm/srtm_taumol21.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "ifsrrtm/srtm_taumol21.f90"
! this file has been modified for the use in icon

subroutine srtm_taumol21 &
 & ( kidia   , kfdia    , klev,&
 & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
 & k_jp      , k_jt     , k_jt1     , p_oneminus,&
 & p_colh2o  , p_colco2 , p_colmol,&
 & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, k_indfor,&
 & p_sfluxzen, p_taug   , p_taur    , prmu0   &
 & )  

!     written by eli j. mlawer, atmospheric & environmental research.

!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)

! modifications
!        m.hamrud      01-oct-2003 cy28 cleaning

!     jjmorcrette 2003-02-24 adapted to ecmwf environment
!        d.salmond  31-oct-2007 vector version in the style of rrtm from meteo france & nec
!     jjmorcrette 20110610 flexible configuration for number of g-points

use parkind1 , only : jpim, jprb
use ecradhook  , only : lhook, dr_hook
use parsrtm  , only : jpg
use yoesrtm  , only : ng21
use yoesrta21, only : absa, absb, forrefc, selfrefc, sfluxrefc, rayl, layreffr, strrat  
use yoesrtwn , only : nspa, nspb

implicit none

!-- output
integer(kind=jpim),intent(in)    :: kidia, kfdia 
integer(kind=jpim),intent(in)    :: klev 
real(kind=jprb)   ,intent(in)    :: p_fac00(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_fac01(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_fac10(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_fac11(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_jp(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_jt(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_jt1(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_oneminus(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: p_colh2o(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_colco2(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_colmol(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_laytrop(kidia:kfdia) 
real(kind=jprb)   ,intent(in)    :: p_selffac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_selffrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_indself(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_forfac(kidia:kfdia,klev) 
real(kind=jprb)   ,intent(in)    :: p_forfrac(kidia:kfdia,klev) 
integer(kind=jpim),intent(in)    :: k_indfor(kidia:kfdia,klev) 

real(kind=jprb)   ,intent(inout) :: p_sfluxzen(kidia:kfdia,jpg) 
real(kind=jprb)   ,intent(inout) :: p_taug(kidia:kfdia,klev,jpg) 
real(kind=jprb)   ,intent(inout) :: p_taur(kidia:kfdia,klev,jpg) 
real(kind=jprb)   ,intent(in)    :: prmu0(kidia:kfdia)
!- from intfac      
!- from intind
!- from precise             
!- from profdata             
!- from self             
integer(kind=jpim) :: ig, ind0, ind1, inds, indf, js, i_lay, i_laysolfr(kidia:kfdia), i_nlayers, iplon
integer(kind=jpim) :: laytrop_min, laytrop_max
integer(kind=jpim) :: i_lay_next

real(kind=jprb) :: z_fac000, z_fac001, z_fac010, z_fac011, z_fac100, z_fac101,&
 & z_fac110, z_fac111, z_fs, z_speccomb, z_specmult, z_specparm, &
 & z_tauray  
real(kind=jprb) :: zhook_handle

    !$acc data create(i_laysolfr) &
    !$acc     present(p_fac00, p_fac01, p_fac10, p_fac11, k_jp, k_jt, k_jt1, &
    !$acc             p_oneminus, p_colh2o, p_colco2, p_colmol, k_laytrop, &
    !$acc             p_selffac, p_selffrac, k_indself, p_forfac, p_forfrac, &
    !$acc             k_indfor, p_sfluxzen, p_taug, p_taur, prmu0)


    laytrop_min = minval(k_laytrop(kidia:kfdia))
    laytrop_max = maxval(k_laytrop(kidia:kfdia))
! # 92 "ifsrrtm/srtm_taumol21.f90"

    i_nlayers = klev
    !$acc wait
    !$acc parallel default(none) async(1)
    !$acc loop gang(static:1) vector
    do iplon = kidia, kfdia
      i_laysolfr(iplon) = k_laytrop(iplon)
    enddo

    !$acc loop seq
    do i_lay = 1, laytrop_min
      !$acc loop gang(static:1) vector private(ind0, ind1, inds, indf, js, z_fs, &
      !$acc   z_speccomb, z_specmult, z_specparm, z_tauray)
       do iplon = kidia, kfdia
         if (k_jp(iplon,i_lay) < layreffr                           &
              &    .and. k_jp(iplon,i_lay+1) >= layreffr)           &
              &    i_laysolfr(iplon) = min(i_lay+1,k_laytrop(iplon))
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = min(p_oneminus(iplon),z_specparm)
         z_specmult = 8._jprb*(z_specparm)
         js = 1 + int(z_specmult)
         z_fs = z_specmult - aint(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(21) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(21) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

        !$acc loop seq
!$nec unroll(ng21)
         do ig = 1 , ng21
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._jprb- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                & ) +                                                       &
                & p_colh2o(iplon,i_lay) *                                   &
                & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(iplon,i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(iplon,i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           if (i_lay == i_laysolfr(iplon))                        &
                &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)      &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         enddo
       enddo
    enddo

    !$acc loop seq
    do i_lay = laytrop_min+1, laytrop_max
      !$acc loop gang(static:1) vector private(ind0, ind1, inds, indf, js, z_fs, &
      !$acc   z_speccomb, z_specmult, z_specparm, z_tauray)
       do iplon = kidia, kfdia
          if (i_lay <= k_laytrop(iplon)) then
            if (k_jp(iplon,i_lay) < layreffr                           &
                 &    .and. k_jp(iplon,i_lay+1) >= layreffr)           &
                 &    i_laysolfr(iplon) = min(i_lay+1,k_laytrop(iplon))
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = min(p_oneminus(iplon),z_specparm)
            z_specmult = 8._jprb*(z_specparm)
            js = 1 + int(z_specmult)
            z_fs = z_specmult - aint(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(21) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(21) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$nec unroll(ng21)
            !$acc loop seq
            do ig = 1 , ng21
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._jprb- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                   & ) +                                                       &
                   & p_colh2o(iplon,i_lay) *                                   &
                   & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(iplon,i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(iplon,i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              if (i_lay == i_laysolfr(iplon))                        &
                   &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)      &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            enddo
          else
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = min(p_oneminus(iplon),z_specparm)
            z_specmult = 4._jprb*(z_specparm)
            js = 1 + int(z_specmult)
            z_fs = z_specmult - aint(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(21) +js
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(21)+js
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$nec unroll(ng21)
            !$acc loop seq
            do ig = 1 , ng21
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._jprb - z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                   & ) +                                                       &
                   & p_colh2o(iplon,i_lay) *                                   &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(iplon,i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))
              p_taur(iplon,i_lay,ig) = z_tauray
            enddo
          endif
       enddo
    enddo

    !$acc loop seq
    do i_lay = laytrop_max+1, i_nlayers
      !$acc loop gang(static:1) vector private(ind0, ind1, indf, js, z_fs, z_speccomb, z_specmult, z_specparm, z_tauray)
       do iplon = kidia, kfdia
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = min(p_oneminus(iplon),z_specparm)
         z_specmult = 4._jprb*(z_specparm)
         js = 1 + int(z_specmult)
         z_fs = z_specmult - aint(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(21) +js
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(21)+js
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
         !$acc loop seq
!$nec unroll(ng21)
         do ig = 1 , ng21
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._jprb - z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                & ) +                                                       &
                & p_colh2o(iplon,i_lay) *                                   &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(iplon,i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))
           p_taur(iplon,i_lay,ig) = z_tauray
         enddo
       enddo
    enddo
    !$acc end parallel

    !$acc wait
    !$acc end data

end subroutine srtm_taumol21
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

