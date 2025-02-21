!----------------------------------------------------------------------------
! Copyright (c) 2002-2020, Atmospheric & Environmental Research, Inc. (AER)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC.,
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.
!                        (http://www.rtweb.aer.com/)
!----------------------------------------------------------------------------

! SPDX-License-Identifier: BSD-3-Clause
! Code has been modified for the use in ICON (and formerly in ECHAM)

!
! Optical depth calculations for RRTM Shortwave Radiation
!
! Description
!     Compute the optical depth by interpolating in ln(pressure),
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.
!
! Eli J. Mlawer, Atmospheric & Environmental Research.
!
!





MODULE mo_srtm_taumol

  USE mo_kind,  ONLY : wp,i4
  USE mo_srtm_config,  ONLY : &
       &  jpg , ng16, ng17, ng18, ng19, ng20, ng21, ng22, &
       &  ng23, ng24, ng25, ng26, ng27, ng28, ng29, nspa, nspb

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: srtm_taumol16, srtm_taumol17, srtm_taumol18, srtm_taumol19, &
       &    srtm_taumol20, srtm_taumol21, srtm_taumol22, srtm_taumol23, &
       &    srtm_taumol24, srtm_taumol25, srtm_taumol26, srtm_taumol27, &
       &    srtm_taumol28, srtm_taumol29


CONTAINS

  !-----------------------------------------------------------------------------
  !>
  !! BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
  !!
  SUBROUTINE srtm_taumol16 &
       & ( icount, kbdim     , klev,&
       & p_fac00   , p_fac01   , p_fac10   , p_fac11,&
       & k_jp      , k_jt      , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colch4  , p_colmol,&
       & k_laytrop , p_selffac , p_selffrac, k_indself , p_forfac  , &
       & p_forfrac , k_indfor  , p_sfluxzen, p_taug    , p_taur &
       & )

    USE mo_yoesrta16, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl, layreffr, strrat1

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colch4(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max


!     write(0,*) "icount=",icount
    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))
!     write(0,*) "laytrop_min, laytrop_max=",laytrop_min, laytrop_max

    i_nlayers = klev
    i_laysolfr(:) = i_nlayers

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         z_speccomb = p_colh2o(iplon,i_lay) + strrat1*p_colch4(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(16) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(16) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG16)
         DO ig = 1, ng16
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            z_speccomb = p_colh2o(iplon,i_lay) + strrat1*p_colch4(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(16) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(16) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG16)
            DO ig = 1, ng16
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr &
                 &  .AND. k_jp(iplon,i_lay) >= layreffr)  i_laysolfr(iplon) = i_lay
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(16) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(16)+ 1
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG16)
            DO ig = 1, ng16
              p_taug(iplon,i_lay,ig) = p_colch4(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0  ,ig) + &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig)  + &
                   & p_fac01(iplon,i_lay) * absb(ind1  ,ig)  + &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
              IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

!     write(0,*) "icount=",icount
!     write(0,*) "laytrop_max+1,i_nlayers=",laytrop_max+1,i_nlayers
    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay-1) < layreffr &
              &  .AND. k_jp(iplon,i_lay) >= layreffr)  i_laysolfr(iplon) = i_lay
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(16) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(16)+ 1
!          write(0,*) "ind0,ind1=",ind0,ind1
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG16)
         DO ig = 1, ng16
           p_taug(iplon,i_lay,ig) = p_colch4(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0  ,ig) + &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig)  + &
                & p_fac01(iplon,i_lay) * absb(ind1  ,ig)  + &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
           IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

!     write(0,*) "srtm_taumol16 ends!"
  END SUBROUTINE srtm_taumol16
  !-----------------------------------------------------------------------------
  !>
  !! BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !
  SUBROUTINE srtm_taumol17 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, &
       & p_forfrac , k_indfor , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta17, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl , layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(:) = i_nlayers

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(17) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(17) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG17)
         DO ig = 1, ng17
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                & )  +                                                      &
                & p_colh2o(iplon,i_lay) *                                   &
                & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(iplon,i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(iplon,i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))

           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(17) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(17) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG17)
            DO ig = 1, ng17
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                   & )  +                                                      &
                   & p_colh2o(iplon,i_lay) *                                   &
                   & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(iplon,i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(iplon,i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))

              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr &
                 & .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 4._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(17)+ js
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(17)+js
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG17)
            DO ig = 1, ng17
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay))+  &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                   & ) +                                                       &
                   & p_colh2o(iplon,i_lay) *                                   &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(iplon,i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))
              IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig,js) &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay-1) < layreffr &
              & .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 4._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(17)+ js
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(17)+js
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG17)
         DO ig = 1, ng17
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay))+  &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                & ) +                                                       &
                & p_colh2o(iplon,i_lay) *                                   &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(iplon,i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))
           IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig,js) &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol17
  !-----------------------------------------------------------------------------
  !>
  !! BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
  !
  SUBROUTINE srtm_taumol18 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colch4 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta18, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl , layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colch4(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                            &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colch4(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(18) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(18) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG18)
         DO ig = 1, ng18
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
           IF (i_lay == i_laysolfr(iplon))                           &
                &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                            &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colch4(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(18) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(18) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG18)
            DO ig = 1, ng18
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
              IF (i_lay == i_laysolfr(iplon))                           &
                   &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                   &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(18) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(18)+ 1
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG18)
            DO ig = 1, ng18
              p_taug(iplon,i_lay,ig) = p_colch4(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig)  +  &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig)   +  &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(18) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(18)+ 1
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG18)
         DO ig = 1, ng18
           p_taug(iplon,i_lay,ig) = p_colch4(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig)  +  &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(iplon,i_lay) * absb(ind1,ig)   +  &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol18
  !-----------------------------------------------------------------------------
  !>
  !! BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
  !
  SUBROUTINE srtm_taumol19 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta19, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl, layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) ::  z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                         &
              & .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
              & i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(19) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(19) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG19)
         DO ig = 1 , ng19
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
           IF (i_lay == i_laysolfr(iplon))                           &
                &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                         &
                 & .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
                 & i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(19) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(19) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG19)
            DO ig = 1 , ng19
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
              IF (i_lay == i_laysolfr(iplon))                           &
                   &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                   &    + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(19) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(19)+ 1
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG19)
            DO ig = 1 , ng19
              p_taug(iplon,i_lay,ig) = p_colco2(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig) +   &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +    &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(19) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(19)+ 1
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG19)
         DO ig = 1 , ng19
           p_taug(iplon,i_lay,ig) = p_colco2(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig) +   &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +    &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol19
  !-----------------------------------------------------------------------------
  !>
  !! BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
  !
  SUBROUTINE srtm_taumol20 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , &
       & p_colh2o  , p_colch4 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta20, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, absch4c, rayl, layreffr

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colch4(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                           &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(20) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(20) + 1
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
         DO ig = 1 , ng20
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *      &
                & ((p_fac00(iplon,i_lay) * absa(ind0,ig) +       &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +       &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +         &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +      &
                & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +  &
                & p_selffrac(iplon,i_lay) *                      &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +   &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +    &
                & p_forfrac(iplon,i_lay) *                       &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))      &
                & + p_colch4(iplon,i_lay) * absch4c(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                           &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(20) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(20) + 1
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
            DO ig = 1 , ng20
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *      &
                   & ((p_fac00(iplon,i_lay) * absa(ind0,ig) +       &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +       &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +         &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +      &
                   & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +  &
                   & p_selffrac(iplon,i_lay) *                      &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +   &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +    &
                   & p_forfrac(iplon,i_lay) *                       &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))      &
                   & + p_colch4(iplon,i_lay) * absch4c(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
            ENDDO
          ELSE
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(20) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(20) +1
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
            DO ig = 1 , ng20
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *   &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig) +     &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +    &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +      &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig) +    &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) + &
                   & p_forfrac(iplon,i_lay) *                    &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig)))) + &
                   & p_colch4(iplon,i_lay) * absch4c(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(20) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(20) +1
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG20)
         DO ig = 1 , ng20
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *   &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig) +     &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +    &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +      &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig) +    &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) + &
                & p_forfrac(iplon,i_lay) *                    &
                & (forrefc(indf+1,ig) - forrefc(indf,ig)))) + &
                & p_colch4(iplon,i_lay) * absch4c(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol20
  !-----------------------------------------------------------------------------
  !>
  !! BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !
  SUBROUTINE srtm_taumol21 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta21, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl, layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm,  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                           &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(21) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(21) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG21)
         DO ig = 1 , ng21
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
           IF (i_lay == i_laysolfr(iplon))                        &
                &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)      &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                           &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(21) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(21) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG21)
            DO ig = 1 , ng21
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
              IF (i_lay == i_laysolfr(iplon))                        &
                   &    p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)      &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 4._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(21) +js
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(21)+js
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG21)
            DO ig = 1 , ng21
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp - z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colco2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 4._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(21) +js
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(21)+js
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG21)
         DO ig = 1 , ng21
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp - z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol21
  !-----------------------------------------------------------------------------
  !>
  !! BAND 29:  7700-8050 cm-1 (low - H2O,O2; high - O2)
  !
  SUBROUTINE srtm_taumol22 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colmol , p_colo2,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta22, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, rayl, layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray, z_o2adj , &
         &      z_o2cont
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev

    !     The following factor is the ratio of total O2 band intensity (lines
    !     and Mate continuum) to O2 band intensity (line only).  It is needed
    !     to adjust the optical depths since the k's include only lines.
    z_o2adj = 1.6_wp
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                           &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         z_o2cont = 4.35e-4_wp*p_colo2(iplon,i_lay)/(350.0_wp*2.0_wp)
         z_speccomb = p_colh2o(iplon,i_lay) +          &
              &     z_o2adj*strrat*p_colo2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(22) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(22) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG22)
         DO ig = 1 , ng22
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))                 &
                & + z_o2cont
           IF (i_lay == i_laysolfr(iplon))                        &
                & p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                           &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)           &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            z_o2cont = 4.35e-4_wp*p_colo2(iplon,i_lay)/(350.0_wp*2.0_wp)
            z_speccomb = p_colh2o(iplon,i_lay) +          &
                 &     z_o2adj*strrat*p_colo2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(22) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(22) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG22)
            DO ig = 1 , ng22
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
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
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))                 &
                   & + z_o2cont
              IF (i_lay == i_laysolfr(iplon))                        &
                   & p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            z_o2cont = 4.35e-4_wp*p_colo2(iplon,i_lay)/(350.0_wp*2.0_wp)
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(22) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(22)+ 1
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG22)
            DO ig = 1 , ng22
              p_taug(iplon,i_lay,ig) = p_colo2(iplon,i_lay) * z_o2adj * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig) +            &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +           &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +             &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig)) +          &
                   & z_o2cont
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         z_o2cont = 4.35e-4_wp*p_colo2(iplon,i_lay)/(350.0_wp*2.0_wp)
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(22) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(22)+ 1
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG22)
         DO ig = 1 , ng22
           p_taug(iplon,i_lay,ig) = p_colo2(iplon,i_lay) * z_o2adj * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig) +            &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +           &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +             &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig)) +          &
                & z_o2cont
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol22
  !-----------------------------------------------------------------------------
  !>
  !! BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
  !
  SUBROUTINE srtm_taumol23 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , &
       & p_colh2o  , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta23, ONLY : absa, forrefc, selfrefc &
         & , sfluxrefc, raylc, layreffr, givfac

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) ::  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                            &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(23) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(23) + 1
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)

!$NEC unroll(NG23)
         DO ig = 1 , ng23
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *         &
                & (givfac * (p_fac00(iplon,i_lay) * absa(ind0,ig) + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +          &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +            &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +         &
                & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +     &
                & p_selffrac(iplon,i_lay) *                         &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +      &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +       &
                & p_forfrac(iplon,i_lay) *                          &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr(iplon)) &
                p_sfluxzen(iplon,ig) = sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                            &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(23) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(23) + 1
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)

!$NEC unroll(NG23)
            DO ig = 1 , ng23
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *         &
                   & (givfac * (p_fac00(iplon,i_lay) * absa(ind0,ig) + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +          &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +            &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +         &
                   & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +     &
                   & p_selffrac(iplon,i_lay) *                         &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +      &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +       &
                   & p_forfrac(iplon,i_lay) *                          &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr(iplon)) &
                   p_sfluxzen(iplon,ig) = sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
!$NEC unroll(NG23)
            DO ig = 1 , ng23
              p_taug(iplon,i_lay,ig) = 0.0_wp
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO ig = 1 , ng23
      DO i_lay = laytrop_max+1, i_nlayers
        DO iplon = 1, icount
          p_taug(iplon,i_lay,ig) = 0.0_wp
          p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol23
  !-----------------------------------------------------------------------------
  !>
  !! BAND 24:  12850-16000 cm-1 (low - H2O,O2; high - O2)
  !
  SUBROUTINE srtm_taumol24 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colh2o  , p_colmol , p_colo2   , p_colo3,&
       & k_laytrop , p_selffac, p_selffrac, k_indself  , p_forfac, p_forfrac, &
       & k_indfor  , p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta24, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, abso3ac, abso3bc, raylac, raylbc, layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) ::  z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr                            &
              &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
         z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(24) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(24) + js
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)

!$NEC unroll(NG24)
         DO ig = 1 , ng24
           z_tauray = p_colmol(iplon,i_lay) * (raylac(ig,js) + &
                & z_fs * (raylac(ig,js+1) - raylac(ig,js)))
           p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                & (                                                         &
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                & ) +                                                       &
                & p_colo3(iplon,i_lay) * abso3ac(ig) +                      &
                & p_colh2o(iplon,i_lay) *                                   &
                & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                & p_selffrac(iplon,i_lay) *                                 &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                & p_forfrac(iplon,i_lay) *                                  &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))
           IF (i_lay == i_laysolfr(iplon))                        &
                & p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) =  z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr                            &
                 &    .AND. k_jp(iplon,i_lay+1) >= layreffr)            &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            z_speccomb = p_colh2o(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
            z_specparm = p_colh2o(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(24) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(24) + js
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)

!$NEC unroll(NG24)
            DO ig = 1 , ng24
              z_tauray = p_colmol(iplon,i_lay) * (raylac(ig,js) + &
                   & z_fs * (raylac(ig,js+1) - raylac(ig,js)))
              p_taug(iplon,i_lay,ig) = z_speccomb *                            &
                   & (                                                         &
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                   & ) +                                                       &
                   & p_colo3(iplon,i_lay) * abso3ac(ig) +                      &
                   & p_colh2o(iplon,i_lay) *                                   &
                   & (p_selffac(iplon,i_lay) * (selfrefc(inds,ig) +            &
                   & p_selffrac(iplon,i_lay) *                                 &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +              &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +               &
                   & p_forfrac(iplon,i_lay) *                                  &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))
              IF (i_lay == i_laysolfr(iplon))                        &
                   & p_sfluxzen(iplon,ig) = sfluxrefc(ig,js)         &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) =  z_tauray
            ENDDO
          ELSE
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(24) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(24)+ 1

!$NEC unroll(NG24)
            DO ig = 1 , ng24
              z_tauray = p_colmol(iplon,i_lay) * raylbc(ig)
              p_taug(iplon,i_lay,ig) = p_colo2(iplon,i_lay) *  &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig) +   &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +    &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig)) + &
                   & p_colo3(iplon,i_lay) * abso3bc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(24) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(24)+ 1

!$NEC unroll(NG24)
         DO ig = 1 , ng24
           z_tauray = p_colmol(iplon,i_lay) * raylbc(ig)
           p_taug(iplon,i_lay,ig) = p_colo2(iplon,i_lay) *  &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig) +   &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +    &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig)) + &
                & p_colo3(iplon,i_lay) * abso3bc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol24
  !-----------------------------------------------------------------------------
  !>
  !! BAND 25:  16000-22650 cm-1 (low - H2O; high - nothing)
  !
  SUBROUTINE srtm_taumol25 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , &
       & p_colh2o  , p_colmol , p_colo3,&
       & k_laytrop,&
       & p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta25, ONLY : absa, sfluxrefc, abso3ac, abso3bc, raylc, layreffr

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) ::  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay) < layreffr .AND.   &
              &    k_jp(iplon,i_lay+1) >= layreffr) &
              &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
!$NEC unroll(NG25)
         DO ig = 1 , ng25
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                & p_colo3(iplon,i_lay) * abso3ac(ig)
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            IF (k_jp(iplon,i_lay) < layreffr .AND.   &
                 &    k_jp(iplon,i_lay+1) >= layreffr) &
                 &    i_laysolfr(iplon) = MIN(i_lay+1,k_laytrop(iplon))
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(25) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(25) + 1
!$NEC unroll(NG25)
            DO ig = 1 , ng25
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absa(ind0,ig)   + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig)  + &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig)    + &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) + &
                   & p_colo3(iplon,i_lay) * abso3ac(ig)
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
!$NEC unroll(NG25)
            DO ig = 1 , ng25
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * abso3bc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO ig = 1 , ng25
      DO i_lay = laytrop_max+1, i_nlayers
        DO iplon = 1, icount
          z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
          p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * abso3bc(ig)
          p_taur(iplon,i_lay,ig) = z_tauray
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol25
  !-----------------------------------------------------------------------------
  !>
  !! BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
  !
  SUBROUTINE srtm_taumol26 &
       & ( icount  , kbdim    , klev,&
       & p_colmol,&
       & k_laytrop , &
       & p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta26, ONLY : sfluxrefc, raylc

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, i_lay, i_nlayers, i_laysolfr(kbdim)
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(1:icount) = k_laytrop(1:icount)

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
!$NEC unroll(NG26)
         DO ig = 1 , ng26
           IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
           p_taug(iplon,i_lay,ig) = 0.0_wp
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
!$NEC unroll(NG26)
            DO ig = 1 , ng26
              IF(i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig)=sfluxrefc(ig)
              p_taug(iplon,i_lay,ig) = 0.0_wp
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ELSE
!$NEC unroll(NG26)
            DO ig = 1 , ng26
              p_taug(iplon,i_lay,ig) = 0.0_wp
              p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO ig = 1 , ng26
       DO i_lay = laytrop_max+1, i_nlayers
         DO iplon = 1, icount
           p_taug(iplon,i_lay,ig) = 0.0_wp
           p_taur(iplon,i_lay,ig) = p_colmol(iplon,i_lay) * raylc(ig)
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol26
  !-----------------------------------------------------------------------------
  !>
  !! BAND 27:  29000-38000 cm-1 (low - O3; high - O3)
  !
  SUBROUTINE srtm_taumol27 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_colmol  , p_colo3,&
       & k_laytrop, p_sfluxzen, p_taug    , p_taur &
       & )

    USE mo_yoesrta27, ONLY : absa, absb, sfluxrefc, raylc, layreffr, scalekur

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(:) = i_nlayers

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(27) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(27) + 1
!$NEC unroll(NG27)
         DO ig = 1 , ng27
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absa(ind0,ig)  + &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) + &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +   &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(27) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(27) + 1
!$NEC unroll(NG27)
            DO ig = 1 , ng27
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absa(ind0,ig)  + &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) + &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +   &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr &
                 &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(27) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(27)+ 1
!$NEC unroll(NG27)
            DO ig = 1 , ng27
              z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
              p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig)  + &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) + &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +   &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
              IF (i_lay == i_laysolfr(iplon)) &
                   & p_sfluxzen(iplon,ig) = scalekur * sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO

    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay-1) < layreffr &
              &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(27) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(27)+ 1
!$NEC unroll(NG27)
         DO ig = 1 , ng27
           z_tauray = p_colmol(iplon,i_lay) * raylc(ig)
           p_taug(iplon,i_lay,ig) = p_colo3(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig)  + &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) + &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +   &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))
           IF (i_lay == i_laysolfr(iplon)) &
                & p_sfluxzen(iplon,ig) = scalekur * sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol27
  !-----------------------------------------------------------------------------
  !>
  !! BAND 28:  38000-50000 cm-1 (low - O3,O2; high - O3,O2)
  !
  SUBROUTINE srtm_taumol28 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     , p_oneminus,&
       & p_colmol  , p_colo2  , p_colo3,&
       & k_laytrop,&
       & p_sfluxzen, p_taug   , p_taur &
       & )

    USE mo_yoesrta28, ONLY : absa, absb, sfluxrefc, rayl, layreffr, strrat

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_oneminus(kbdim)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colo3(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, js, i_lay, i_nlayers, i_laysolfr(kbdim)
    REAL(wp) :: z_fs, z_speccomb, z_specmult, z_specparm, z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(:) = i_nlayers

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
         z_specparm = p_colo3(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 8._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(28) + js
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(28) + js
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG28)
         DO ig = 1 , ng28
           p_taug(iplon,i_lay,ig) = z_speccomb * &
                & (&
                & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                & )
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
            z_specparm = p_colo3(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 8._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(28) + js
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(28) + js
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG28)
            DO ig = 1 , ng28
              p_taug(iplon,i_lay,ig) = z_speccomb * &
                   & (&
                   & (1._wp- z_fs) * ( absa(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absa(ind0+9,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absa(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absa(ind1+9,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absa(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absa(ind0+10,ig) * p_fac10(iplon,i_lay) + &
                   &                 absa(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absa(ind1+10,ig) * p_fac11(iplon,i_lay) ) &
                   & )
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr &
                 &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
            z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
            z_specparm = p_colo3(iplon,i_lay)/z_speccomb
            z_specparm = MIN(p_oneminus(iplon),z_specparm)
            z_specmult = 4._wp*(z_specparm)
            js = 1 + INT(z_specmult)
            z_fs = z_specmult - AINT(z_specmult)
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(28)+ js
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(28)+js
            z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG28)
            DO ig = 1 , ng28
              p_taug(iplon,i_lay,ig) = z_speccomb * &
                   & (&
                   & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                   &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                   &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay) )+ &
                   & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                   &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                   &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                   &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                   & )
              IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig,js) &
                   & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO


    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay-1) < layreffr &
              &    .AND. k_jp(iplon,i_lay) >= layreffr) i_laysolfr(iplon) = i_lay
         z_speccomb = p_colo3(iplon,i_lay) + strrat*p_colo2(iplon,i_lay)
         z_specparm = p_colo3(iplon,i_lay)/z_speccomb
         z_specparm = MIN(p_oneminus(iplon),z_specparm)
         z_specmult = 4._wp*(z_specparm)
         js = 1 + INT(z_specmult)
         z_fs = z_specmult - AINT(z_specmult)
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(28)+ js
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(28)+js
         z_tauray = p_colmol(iplon,i_lay) * rayl

!$NEC unroll(NG28)
         DO ig = 1 , ng28
           p_taug(iplon,i_lay,ig) = z_speccomb * &
                & (&
                & (1._wp- z_fs) * ( absb(ind0,ig) * p_fac00(iplon,i_lay) +    &
                &                 absb(ind0+5,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1,ig) * p_fac01(iplon,i_lay) +    &
                &                 absb(ind1+5,ig) * p_fac11(iplon,i_lay) )+ &
                & z_fs        * ( absb(ind0+1,ig) * p_fac00(iplon,i_lay) +  &
                &                 absb(ind0+6,ig) * p_fac10(iplon,i_lay) +  &
                &                 absb(ind1+1,ig) * p_fac01(iplon,i_lay) +  &
                &                 absb(ind1+6,ig) * p_fac11(iplon,i_lay) )  &
                & )
           IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig,js) &
                & + z_fs * (sfluxrefc(ig,js+1) - sfluxrefc(ig,js))
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE srtm_taumol28
  !-----------------------------------------------------------------------------
  !>
  !! BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
  !
  SUBROUTINE srtm_taumol29 &
       & ( icount, kbdim    , klev,&
       & p_fac00   , p_fac01  , p_fac10   , p_fac11,&
       & k_jp      , k_jt     , k_jt1     ,&
       & p_colh2o  , p_colco2 , p_colmol,&
       & k_laytrop , p_selffac, p_selffrac, k_indself , p_forfac, p_forfrac, &
       & k_indfor, p_sfluxzen, p_taug   , p_taur  &
       & )

    USE mo_yoesrta29, ONLY : absa, absb, forrefc, selfrefc &
         & , sfluxrefc, absh2oc, absco2c, rayl, layreffr

    INTEGER(i4),INTENT(in) :: icount
    INTEGER,INTENT(in)    :: kbdim
    INTEGER,INTENT(in)    :: klev
    INTEGER,INTENT(in)    :: k_jp(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt(kbdim,klev)
    INTEGER,INTENT(in)    :: k_jt1(kbdim,klev)
    INTEGER,INTENT(in)    :: k_laytrop(kbdim)
    INTEGER,INTENT(in)    :: k_indself(kbdim,klev)
    INTEGER,INTENT(in)    :: k_indfor(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac00(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac01(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac10(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_fac11(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colh2o(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colco2(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_colmol(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_selffrac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfac(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: p_forfrac(kbdim,klev)

    REAL(wp)   ,INTENT(out)   :: p_sfluxzen(kbdim,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taug(kbdim,klev,jpg)
    REAL(wp)   ,INTENT(out)   :: p_taur(kbdim,klev,jpg)

    INTEGER(i4) :: iplon
    INTEGER :: ig, ind0, ind1, inds, indf, i_lay, i_nlayers
    INTEGER :: i_laysolfr(kbdim)
    REAL(wp) ::  z_tauray
    INTEGER :: laytrop_min, laytrop_max

    laytrop_min = MINVAL(k_laytrop(1:icount))
    laytrop_max = MAXVAL(k_laytrop(1:icount))

    i_nlayers = klev
    i_laysolfr(:) = i_nlayers

    DO i_lay = 1, laytrop_min
       DO iplon = 1, icount
         ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(29) + 1
         ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(29) + 1
         inds = k_indself(iplon,i_lay)
         indf = k_indfor(iplon,i_lay)
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG29)
         DO ig = 1, ng29
           p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *     &
                & ((p_fac00(iplon,i_lay) * absa(ind0,ig) +      &
                & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +      &
                & p_fac01(iplon,i_lay) * absa(ind1,ig) +        &
                & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +     &
                & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) + &
                & p_selffrac(iplon,i_lay) *                     &
                & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +  &
                & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +   &
                & p_forfrac(iplon,i_lay) *                      &
                & (forrefc(indf+1,ig) - forrefc(indf,ig))))     &
                & + p_colco2(iplon,i_lay) * absco2c(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    DO i_lay = laytrop_min+1, laytrop_max
       DO iplon = 1, icount
          IF (i_lay <= k_laytrop(iplon)) THEN
            ind0 = ((k_jp(iplon,i_lay)-1)*5+(k_jt(iplon,i_lay)-1))*nspa(29) + 1
            ind1 = (k_jp(iplon,i_lay)*5+(k_jt1(iplon,i_lay)-1))*nspa(29) + 1
            inds = k_indself(iplon,i_lay)
            indf = k_indfor(iplon,i_lay)
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG29)
            DO ig = 1, ng29
              p_taug(iplon,i_lay,ig) = p_colh2o(iplon,i_lay) *     &
                   & ((p_fac00(iplon,i_lay) * absa(ind0,ig) +      &
                   & p_fac10(iplon,i_lay) * absa(ind0+1,ig) +      &
                   & p_fac01(iplon,i_lay) * absa(ind1,ig) +        &
                   & p_fac11(iplon,i_lay) * absa(ind1+1,ig)) +     &
                   & p_selffac(iplon,i_lay) * (selfrefc(inds,ig) + &
                   & p_selffrac(iplon,i_lay) *                     &
                   & (selfrefc(inds+1,ig) - selfrefc(inds,ig))) +  &
                   & p_forfac(iplon,i_lay) * (forrefc(indf,ig) +   &
                   & p_forfrac(iplon,i_lay) *                      &
                   & (forrefc(indf+1,ig) - forrefc(indf,ig))))     &
                   & + p_colco2(iplon,i_lay) * absco2c(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ELSE
            IF (k_jp(iplon,i_lay-1) < layreffr                                &
                 &   .AND. k_jp(iplon,i_lay) >= layreffr)  i_laysolfr(iplon) = i_lay
            ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(29) + 1
            ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(29)+ 1
            z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG29)
            DO ig = 1 , ng29
              p_taug(iplon,i_lay,ig) = p_colco2(iplon,i_lay) * &
                   & (p_fac00(iplon,i_lay) * absb(ind0,ig) +   &
                   & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                   & p_fac01(iplon,i_lay) * absb(ind1,ig) +    &
                   & p_fac11(iplon,i_lay) * absb(ind1+1,ig))   &
                   & + p_colh2o(iplon,i_lay) * absh2oc(ig)
              IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig)
              p_taur(iplon,i_lay,ig) = z_tauray
            ENDDO
          ENDIF
       ENDDO
    ENDDO



    DO i_lay = laytrop_max+1, i_nlayers
       DO iplon = 1, icount
         IF (k_jp(iplon,i_lay-1) < layreffr                                &
              &   .AND. k_jp(iplon,i_lay) >= layreffr)  i_laysolfr(iplon) = i_lay
         ind0 = ((k_jp(iplon,i_lay)-13)*5+(k_jt(iplon,i_lay)-1))*nspb(29) + 1
         ind1 = ((k_jp(iplon,i_lay)-12)*5+(k_jt1(iplon,i_lay)-1))*nspb(29)+ 1
         z_tauray = p_colmol(iplon,i_lay) * rayl
!$NEC unroll(NG29)
         DO ig = 1 , ng29
           p_taug(iplon,i_lay,ig) = p_colco2(iplon,i_lay) * &
                & (p_fac00(iplon,i_lay) * absb(ind0,ig) +   &
                & p_fac10(iplon,i_lay) * absb(ind0+1,ig) +  &
                & p_fac01(iplon,i_lay) * absb(ind1,ig) +    &
                & p_fac11(iplon,i_lay) * absb(ind1+1,ig))   &
                & + p_colh2o(iplon,i_lay) * absh2oc(ig)
           IF (i_lay == i_laysolfr(iplon)) p_sfluxzen(iplon,ig) = sfluxrefc(ig)
           p_taur(iplon,i_lay,ig) = z_tauray
         ENDDO
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
  END SUBROUTINE srtm_taumol29

END MODULE mo_srtm_taumol
