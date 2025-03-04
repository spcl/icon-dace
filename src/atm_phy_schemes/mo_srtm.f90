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

! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! For runs that check result consistency we fix the different
! contractions that the Intel compiler performs on some loops (at
! least in version 16.0) for the vectorized part and the
! non-vectorized parts
MODULE mo_srtm

  USE mo_exception,   ONLY: finish

  USE mo_kind, ONLY : wp, i4
  USE mo_physical_constants, ONLY : rd, rgrav

  USE mo_srtm_config, ONLY : delwave, wavenum2, wavenum1,            &
    &    ngc, jpinpx, jpb1, jpb2, preflog, tref, repclc, replog,     &
    &    ssi_default
  USE mo_srtm_taumol, ONLY :                                         &
    &    srtm_taumol16, srtm_taumol17, srtm_taumol18, srtm_taumol19, &
    &    srtm_taumol20, srtm_taumol21, srtm_taumol22, srtm_taumol23, &
    &    srtm_taumol24, srtm_taumol25, srtm_taumol26, srtm_taumol27, &
    &    srtm_taumol28, srtm_taumol29
  USE mo_radiation_config, ONLY: icld_overlap
  USE mo_nwp_tuning_config, ONLY: tune_difrad_3dcont

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: srtm_srtm_224gp


  REAL(wp), PARAMETER :: nir_vis_boundary   = 14500._wp
  REAL(wp), PARAMETER :: par_lower_boundary = 14285.7143_wp ! equivalent to
  !                                                         ! 700nm wavelength
  REAL(wp), PARAMETER :: par_upper_boundary = 25000._wp     ! equivalent to
  !                                                         ! 400nm wavelength

  ! replaced with namelist parameter icld_overlap
! INTEGER,  PARAMETER :: i_overlap = 1       ! 1: maximum-random overlap
                                             ! 2: generalized overlap (Hogan, Illingworth, 2000)
                                             ! 3: maximum overlap
                                             ! 4: random overlap
  REAL(wp), PARAMETER :: zdecorr = 2000.0_wp ! decorrelation length scale del(z0)

CONTAINS

  SUBROUTINE srtm_srtm_224gp                                                   &
                                !  input
    & (kproma          , kbdim           , klev            , ksw             , &
    &  alb_vis_dir     , alb_nir_dir     , alb_vis_dif     , alb_nir_dif     , &
    &  pm_fl_vr        , tk_fl_vr        , prmu0                             , &
    &  col_dry_vr      , wkl_vr                                              , &
    &  cld_frc_vr      , cld_tau_sw_vr   , cld_cg_sw_vr    , cld_piz_sw_vr   , &
    &  aer_tau_sw_vr   , aer_cg_sw_vr    , aer_piz_sw_vr                     , &
    &  ssi                                                                   , &
                                !  output
    &  flxd_sw         , flxu_sw         , flxd_sw_clr     , flxu_sw_clr     , &
                                ! optional output
    &  flxd_dff_sfc    , flxd_par_sfc    , vis_frc_sfc                       , &
    &  nir_dff_frc_sfc , vis_dff_frc_sfc , par_dff_frc_sfc                     )


    !-- Interface to RRTM_SW
    !     JJMorcrette 030225
    !     JJMorcrette 20071015 3D fields of CO2, CH4, N2O and NO2
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !     define bands 24-28 as visible
    !     TJ Raddatz  20091223 diagnosis of photosynthetically active radiation (PAR)
    !     TJ Raddatz  20100111 passing visible and NIR surface albedo to
    !        radiative transfer calculations
    !

    !-- Input arguments

    INTEGER,  INTENT(in)    :: kbdim
    INTEGER,  INTENT(in)    :: klev
    INTEGER,  INTENT(in)    :: ksw
    INTEGER,  INTENT(in)    :: kproma
    REAL(wp), INTENT(in)    :: alb_vis_dir(kbdim)
    REAL(wp), INTENT(in)    :: alb_nir_dir(kbdim)
    REAL(wp), INTENT(in)    :: alb_vis_dif(kbdim)
    REAL(wp), INTENT(in)    :: alb_nir_dif(kbdim)
    REAL(wp), INTENT(in)    :: prmu0(kbdim)
    REAL(wp), INTENT(in)    :: pm_fl_vr(kbdim,klev)
    REAL(wp), INTENT(in)    :: tk_fl_vr(kbdim,klev)
    REAL(wp), INTENT(in)    :: col_dry_vr(kbdim,klev)
    REAL(wp), INTENT(in)    :: wkl_vr(kbdim,jpinpx,klev)
    REAL(wp), INTENT(in)    :: cld_frc_vr(kbdim,klev)
    REAL(wp), INTENT(in)    :: cld_tau_sw_vr(kbdim,ksw,klev)
    REAL(wp), INTENT(in)    :: cld_cg_sw_vr(kbdim,ksw,klev)
    REAL(wp), INTENT(in)    :: cld_piz_sw_vr(kbdim,ksw,klev)
    ! hs: The order of indices in the aerosol optical properties was wrong. However,
    ! it should be checked if a fix in the interface would be more appropriate.
    !   REAL(wp), INTENT(in)    :: aer_tau_sw_vr(kbdim,ksw,klev)
    !   REAL(wp), INTENT(in)    :: aer_cg_sw_vr(kbdim,ksw,klev)
    !   REAL(wp), INTENT(in)    :: aer_piz_sw_vr(kbdim,ksw,klev)
    REAL(wp), INTENT(in)    :: aer_tau_sw_vr(kbdim,klev,ksw)
    REAL(wp), INTENT(in)    :: aer_cg_sw_vr(kbdim,klev,ksw)
    REAL(wp), INTENT(in)    :: aer_piz_sw_vr(kbdim,klev,ksw)
    REAL(wp), INTENT(in)    :: ssi(ksw)

    !-- Output arguments

    REAL(wp), INTENT(out)   :: flxd_sw(kbdim,klev+1)     !< downward flux total sky
    REAL(wp), INTENT(out)   :: flxd_sw_clr(kbdim,klev+1) !< downward flux clear sky
    REAL(wp), INTENT(out)   :: flxu_sw(kbdim,klev+1)     !< upward flux total sky
    REAL(wp), INTENT(out)   :: flxu_sw_clr(kbdim,klev+1) !< upward flux clear sky

    REAL(wp), INTENT(out), OPTIONAL :: &
      & flxd_dff_sfc(kbdim),     & !< surface downward diffuse rad
      & flxd_par_sfc(kbdim),     & !< surface downward photosynthetically active rad
      & vis_frc_sfc(kbdim),      & !< Visible fraction of net surface radiation
      & nir_dff_frc_sfc(kbdim),  & !< Diffuse fraction of downward surface near-infrared radiation
      & vis_dff_frc_sfc(kbdim),  & !< Diffuse fraction of downward surface visible radiation
      & par_dff_frc_sfc(kbdim)     !< Diffuse fraction of downward surface PAR

    !-- Local variables

    REAL(wp) :: z_colmol(kbdim,klev) ,  z_co2mult(kbdim,klev)
    REAL(wp) :: z_colch4(kbdim,klev) , z_colco2(kbdim,klev)
    REAL(wp) :: z_colh2o(kbdim,klev) , z_colo3(kbdim,klev)
    REAL(wp) :: z_coln2o(kbdim,klev) , z_colo2(kbdim,klev)
    REAL(wp) :: z_forfac(kbdim,klev) , z_forfrac(kbdim,klev)
    REAL(wp) :: z_selffrac(kbdim,klev), z_selffac(kbdim,klev)
    REAL(wp) :: z_fac00(kbdim,klev)  , z_fac01(kbdim,klev)
    REAL(wp) :: z_fac11(kbdim,klev)  , z_fac10(kbdim,klev)
    REAL(wp) :: zfrcl(kbdim,klev)
    REAL(wp) :: z_oneminus(kbdim)    , bnd_wght(ksw)
    REAL(wp) :: ztauc(kbdim,klev,ksw), ztaua(kbdim,klev,ksw)
    REAL(wp) :: zasyc(kbdim,klev,ksw), zasya(kbdim,klev,ksw)
    REAL(wp) :: zomgc(kbdim,klev,ksw), zomga(kbdim,klev,ksw)

    REAL(wp) :: zbbcd(kbdim,klev+1,ksw), zbbcu(kbdim,klev+1,ksw)
    REAL(wp) :: zbbfd(kbdim,klev+1,ksw), zbbfu(kbdim,klev+1,ksw)
    REAL(wp) :: zsudu(kbdim,ksw), zsuduc(kbdim,ksw)

    REAL(wp) :: zpm_fl_vr(kbdim,klev)
    REAL(wp) :: ztk_fl_vr(kbdim,klev)
    REAL(wp) :: zcol_dry_vr(kbdim,klev)
    REAL(wp) :: zwkl_vr(kbdim,jpinpx,klev)
    REAL(wp) :: zflxd_sw(kbdim,klev+1)
    REAL(wp) :: zflxd_sw_clr(kbdim,klev+1)
    REAL(wp) :: zflxd_sw_cld(kbdim,klev+1)
    REAL(wp) :: zflxu_sw(kbdim,klev+1)
    REAL(wp) :: zflxu_sw_clr(kbdim,klev+1)
    REAL(wp) :: zflxu_sw_cld(kbdim,klev+1)

    INTEGER  :: i_laytrop(kbdim), i_layswtch(kbdim), i_laylow(kbdim)
    INTEGER  :: indfor(kbdim,klev), indself(kbdim,klev)
    INTEGER  :: jp(kbdim,klev), jt(kbdim,klev), jt1(kbdim,klev)

    REAL(wp) :: zclear(kbdim), zcloud(kbdim), zeps, zfrcl_above(kbdim), zcloud_ll(kbdim), lcfac
    REAL(wp) :: zalbd(kbdim,ksw) , zalbp(kbdim,ksw)

    REAL(wp) :: frc_vis(ksw), frc_nir(ksw), frc_par
    REAL(wp) :: zfvis, zfnir, zfpar, total
    REAL(wp) :: zflxn_vis(kbdim), zflxn(kbdim), zflxd_vis(kbdim), zflxd_nir(kbdim), zflxd_par(kbdim), &
                zflxd_diff(kbdim), zflxd_vis_diff(kbdim), zflxd_nir_diff(kbdim), zflxd_par_diff(kbdim), &
                zrat_swdn(kbdim)

    INTEGER  :: icldatm, inflag, iceflag, i_liqflag, i_nstr
    INTEGER(i4) :: idx(kbdim)
    INTEGER(i4) :: icount, ic, jl, jk, jsw, jb, jk1, jkp1
    REAL(wp) :: zrmu0(kbdim)
    REAL(wp) :: ccmax, ccran, alpha, deltaz, ccrat
    LOGICAL  :: lcomp_fractions



    !-----------------------------------------------------------------------
    !-- calculate information needed ny the radiative transfer routine

    zeps       = 1.e-06_wp
    z_oneminus = 1.0_wp - zeps
    !

    !++hs
    ! --- weight radiation within a band for the solar cycle ---
    ! ssi contains the solar irradiation at 1 AU distance from the sun
    ! in each band. The sum over all bands of ssi(:) is the TSI.
    bnd_wght(:) = ssi(:) / ssi_default(:)
    !--hs

    icldatm   = 1
    inflag    = 2
    iceflag   = 3
    i_liqflag = 1
    i_nstr    = 2

    !-------------------------------
    ! scatter-gather in idx
    ic=0
    DO jl = 1,kproma
      ! The threshold value of 0.0 ensures that radiation is not calculated for night points
      IF (prmu0(jl) > 0.0_wp) THEN
        ic=ic+1
        idx(ic)=jl
        zrmu0(ic) = prmu0(jl)
      ENDIF
    ENDDO
    icount=ic

    DO jk=1,klev+1
      DO jl = 1, kproma
        flxu_sw(jl,jk)     = 0.0_wp
        flxd_sw(jl,jk)     = 0.0_wp
        flxu_sw_clr(jl,jk) = 0.0_wp
        flxd_sw_clr(jl,jk) = 0.0_wp
      END DO
    END DO

    IF (PRESENT(flxd_dff_sfc)) THEN
      DO jl = 1, kproma
        flxd_dff_sfc(jl) = 0.0_wp
      ENDDO
    ENDIF

    IF (PRESENT(flxd_par_sfc)) THEN
      DO jl = 1, kproma
        flxd_par_sfc(jl) = 0.0_wp
      ENDDO
    ENDIF

    IF (PRESENT(vis_frc_sfc) .AND. PRESENT(nir_dff_frc_sfc) .AND. &
        PRESENT(vis_dff_frc_sfc) .AND. PRESENT(par_dff_frc_sfc)) THEN
      lcomp_fractions = .TRUE.
      DO jl = 1, kproma
        vis_frc_sfc(jl)     = 0.0_wp
        nir_dff_frc_sfc(jl) = 0.0_wp
        vis_dff_frc_sfc(jl) = 0.0_wp
        par_dff_frc_sfc(jl) = 0.0_wp
      END DO
    ELSE
      lcomp_fractions = .FALSE.
    END IF

    IF (icount == 0) RETURN

    !-------------------------------

    DO jk = 1, klev
!$NEC ivdep
      DO ic = 1, icount
        jl = idx(ic)
        zfrcl(ic,jk)       = cld_frc_vr(jl,jk)
        zpm_fl_vr(ic,jk)   = pm_fl_vr(jl,jk)
        ztk_fl_vr(ic,jk)   = tk_fl_vr(jl,jk)
        zcol_dry_vr(ic,jk) = col_dry_vr(jl,jk)
!$NEC unroll(jpinpx)
        zwkl_vr(ic,1:jpinpx,jk) = wkl_vr(jl,1:jpinpx,jk)
      ENDDO

    END DO


    CALL srtm_setcoef                                                   &
                                !   input
      & ( icount, idx,   kbdim    , klev,                               &
      &   zpm_fl_vr  , ztk_fl_vr ,                                      &
      &   zcol_dry_vr, zwkl_vr   ,                                      &
                                !   output
      &   i_laytrop , i_layswtch  , i_laylow ,                          &
      &   z_co2mult , z_colch4 , z_colco2    , z_colh2o    , z_colmol , &
      &   z_coln2o  , z_colo2  , z_colo3     ,                          &
      &   z_forfac  , z_forfrac, indfor      ,                          &
      &   z_selffac , z_selffrac,indself     ,                          &
      &   z_fac00   , z_fac01  , z_fac10     , z_fac11     ,            &
      &   jp        , jt       , jt1                                    )


    !- call the radiation transfer routine

    ! surface albedo, direct/parallel beam (p) and diffuse (d)
    DO jsw=1,ksw
      frc_vis(jsw) = MAX(0.0_wp, MIN(1.0_wp, &
        (wavenum2(jsw+jpb1-1) - nir_vis_boundary) / delwave(jsw+jpb1-1)))
      frc_nir(jsw) = 1.0_wp - frc_vis(jsw)
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      DO ic = 1, icount
        jl = idx(ic)
        IF (nir_vis_boundary - wavenum1(jsw+jpb1-1) < 0.0_wp) THEN
          zalbd(ic,jsw) = alb_vis_dif(jl)
          zalbp(ic,jsw) = alb_vis_dir(jl)
        ELSEIF (nir_vis_boundary - wavenum2(jsw+jpb1-1) > 0.0_wp) THEN
          zalbd(ic,jsw) = alb_nir_dif(jl)
          zalbp(ic,jsw) = alb_nir_dir(jl)
        ELSE
          zalbd(ic,jsw) = frc_vis(jsw) * alb_vis_dif(jl) &
            &           + frc_nir(jsw) * alb_nir_dif(jl)
          zalbp(ic,jsw) = frc_vis(jsw) * alb_vis_dir(jl) &
            &           + frc_nir(jsw) * alb_nir_dir(jl)
        ENDIF
      ENDDO

      ! optical properties of clouds and aerosols
      DO jk=1,klev
!IBM* ASSERT(NODEPS)
        DO ic = 1, icount
          jl = idx(ic)
          ztauc(ic,jk,jsw) = cld_tau_sw_vr(jl,jsw,jk)
          zasyc(ic,jk,jsw) = cld_cg_sw_vr(jl,jsw,jk)
          zomgc(ic,jk,jsw) = cld_piz_sw_vr(jl,jsw,jk)
          ! hs: The order of indices in the aerosol optical properties was wrong. However,
          ! it should be checked if a fix in the interface would be more appropriate.
          !           ztaua(ic,jk,jsw) = aer_tau_sw_vr(jl,jsw,jk)
          !           zasya(ic,jk,jsw) = aer_cg_sw_vr(jl,jsw,jk)
          !           zomga(ic,jk,jsw) = aer_piz_sw_vr(jl,jsw,jk)
          ztaua(ic,jk,jsw) = aer_tau_sw_vr(jl,jk,jsw)
          zasya(ic,jk,jsw) = aer_cg_sw_vr(jl,jk,jsw)
          zomga(ic,jk,jsw) = aer_piz_sw_vr(jl,jk,jsw)
        ENDDO
      ENDDO
    ENDDO

    DO jsw=1,ksw
      DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
        DO ic = 1, icount
          zbbcu(ic,jk,jsw)=0.0_wp
          zbbcd(ic,jk,jsw)=0.0_wp
          zbbfu(ic,jk,jsw)=0.0_wp
          zbbfd(ic,jk,jsw)=0.0_wp
        ENDDO
      ENDDO
    ENDDO


    zsudu=0.0_wp
    zsuduc=0.0_wp

    CALL srtm_spcvrt                                                       &
                                !   input
      & ( icount,    kbdim       , klev       , ksw,                       &
      &  z_oneminus,                                                       &
      &   zalbd     , zalbp       , zfrcl      ,                           &
      &   ztauc     , zasyc       , zomgc      ,                           &
      &   ztaua     , zasya       , zomga      ,                           &
      &   zrmu0     , i_laytrop   ,                                        &
      &   z_colch4  , z_colco2    , z_colh2o   ,                           &
      &   z_colmol  , z_colo2     , z_colo3    ,                           &
      &   z_forfac  , z_forfrac   , indfor     ,                           &
      &   z_selffac , z_selffrac  , indself    ,                           &
      &   z_fac00   , z_fac01     , z_fac10    , z_fac11    ,              &
      &   jp        , jt          , jt1        ,                           &
      &   zbbfd     , zbbfu       , zbbcd      , zbbcu      ,              &
      &   zsudu     , zsuduc      )

    DO jk=1,klev+1
      DO ic = 1, icount
        zflxu_sw_cld(ic,jk) = 0.0_wp
        zflxd_sw_cld(ic,jk) = 0.0_wp
        zflxu_sw_clr(ic,jk) = 0.0_wp
        zflxd_sw_clr(ic,jk) = 0.0_wp
      END DO
    END DO


    IF (PRESENT(flxd_dff_sfc)) THEN ! compute diffuse parts of surface radiation
      zflxd_diff(:) = 0._wp
    ENDIF

    IF (PRESENT(flxd_par_sfc) .OR. lcomp_fractions) THEN ! compute photosynthetically active parts of surface radiation
      zflxd_par(:)  = 0._wp
    ENDIF

    IF (lcomp_fractions) THEN
      zflxn_vis(:)      = 0._wp
      zflxn(:)          = 0._wp
      zflxd_vis(:)      = 0._wp
      zflxd_nir(:)      = 0._wp
      zflxd_vis_diff(:) = 0._wp
      zflxd_nir_diff(:) = 0._wp
      zflxd_par_diff(:) = 0._wp
    END IF

    DO jb = 1,ksw

      ! sum up fluxes over all bands
      DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
        DO ic = 1, icount
          zflxu_sw_clr(ic,jk)=zflxu_sw_clr(ic,jk)+bnd_wght(jb)*zbbcu(ic,jk,jb)
          zflxd_sw_clr(ic,jk)=zflxd_sw_clr(ic,jk)+bnd_wght(jb)*zbbcd(ic,jk,jb)
          zflxu_sw_cld(ic,jk)=zflxu_sw_cld(ic,jk)+bnd_wght(jb)*zbbfu(ic,jk,jb)
          zflxd_sw_cld(ic,jk)=zflxd_sw_cld(ic,jk)+bnd_wght(jb)*zbbfd(ic,jk,jb)

        ENDDO
      ENDDO

    ENDDO


    !
    ! --- overlap computation
    !
!IBM* ASSERT(NODEPS)
    DO ic = 1, icount
      zclear(ic)     = 1.0_wp
      zcloud(ic)     = 0.0_wp
      zcloud_ll(ic)  = 0.0_wp
      zfrcl_above(ic)= 0.0_wp
    ENDDO

    SELECT CASE ( icld_overlap )

    CASE ( 1 )   ! maximum-random overlap
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
      DO jk = 1, klev
        DO ic = 1, icount
          zclear(ic) = zclear(ic)*(1.0_wp-MAX(zfrcl(ic,jk),zfrcl_above(ic))) &
            & /(1.0_wp-MIN(zfrcl_above(ic),1.0_wp-zeps))
          zfrcl_above(ic) = zfrcl(ic,jk)
          zcloud(ic) = 1.0_wp-zclear(ic)
        ENDDO
      ENDDO

    CASE ( 2 )   ! generalized maximum-random overlap (Hogan, Illingworth, 2000)
      DO ic = 1, icount
        zrat_swdn(ic) = zflxd_sw_cld(ic,klev+1)/zflxd_sw_clr(ic,klev+1)
      ENDDO

      DO jk = klev, 1, -1
        jkp1 = MIN(jk+1,klev)
        jk1 = klev+2-jk

        DO ic = 1, icount
          ! reduction factor for thin cirrus clouds lying above optically thicker water clouds
          ccrat = MIN(1._wp,8._wp*MAX(5.e-5_wp,1._wp-zflxd_sw_cld(ic,jk1)/zflxd_sw_clr(ic,jk1))/&
                                   MAX(1.e-3_wp,1._wp-zrat_swdn(ic)) )
          ccrat = MAX(ccrat,MIN(1._wp,0.1_wp*(ztk_fl_vr(ic,jkp1)-233.15_wp))) ! no reduction for T > -30C
          ccmax = MAX( ccrat*zfrcl(ic,jk),  zcloud(ic) )
          ccran =      ccrat*zfrcl(ic,jk) + zcloud(ic) - zfrcl(ic,jk) * zcloud(ic)

          ! layer thickness [m] between level jk and next upper (!) level jk+1
          deltaz = (zpm_fl_vr(ic,jk)-zpm_fl_vr(ic,jkp1))/(zpm_fl_vr(ic,jkp1)+zpm_fl_vr(ic,jk)) * &
                   (ztk_fl_vr(ic,jkp1)+ztk_fl_vr(ic,jk))*rd*rgrav

          alpha  = MIN(EXP(-deltaz/zdecorr), zfrcl(ic,jkp1)/MAX(zeps,zfrcl(ic,jk)) )

          zcloud(ic) = alpha * ccmax + (1-alpha) * ccran
          zclear(ic) = 1.0_wp-zcloud(ic)

          ! additional calculation of low-level cloud cover for postprocessing of diffuse fraction
          lcfac = MERGE(1.0_wp, 0.0_wp, zpm_fl_vr(ic,jk)/zpm_fl_vr(ic,1) > 0.75_wp)
          ccmax = MAX( zfrcl(ic,jk),  zcloud_ll(ic) )
          ccran =      zfrcl(ic,jk) + zcloud_ll(ic) - zfrcl(ic,jk) * zcloud_ll(ic)
          zcloud_ll(ic) = lcfac * (alpha * ccmax + (1-alpha) * ccran)

        ENDDO
      ENDDO

    CASE ( 3 )   ! maximum overlap
      DO jk = 1, klev
!IBM* ASSERT(NODEPS)
        DO ic = 1, icount
          zcloud(ic) = MAX(zcloud(ic),zfrcl(ic,jk))
          zclear(ic) = 1.0_wp-zcloud(ic)
        ENDDO
      ENDDO

    CASE ( 4 )   ! random overlap
      DO jk = 1, klev
!IBM* ASSERT(NODEPS)
        DO ic = 1, icount
          zclear(ic) = zclear(ic)*(1.0_wp-zfrcl(ic,jk))
          zcloud(ic) = 1.0_wp-zclear(ic)
        ENDDO
      ENDDO

    END SELECT

    DO jk = 1, klev+1
      DO ic = 1, icount
        zflxu_sw(ic,jk) = zcloud(ic)*zflxu_sw_cld(ic,jk) + zclear(ic)*zflxu_sw_clr(ic,jk)
        zflxd_sw(ic,jk) = zcloud(ic)*zflxd_sw_cld(ic,jk) + zclear(ic)*zflxd_sw_clr(ic,jk)
      ENDDO
    ENDDO

    DO jb = 1,ksw

      IF (jb == 9) THEN
        frc_par = 0.533725_wp
      ELSE IF (jb == 10) THEN
        frc_par = 1.0_wp
      ELSE IF (jb == 11) THEN
        frc_par = 0.550164_wp
      ELSE
        frc_par = 0._wp
      ENDIF

      IF (PRESENT(flxd_dff_sfc)) THEN ! compute diffuse parts of surface radiation
        DO ic = 1, icount
          zflxd_diff(ic) = zflxd_diff(ic) + bnd_wght(jb)*( zcloud(ic)*(zbbfd(ic,klev+1,jb)-zsudu(ic,jb))  &
            &                                            + zclear(ic)*(zbbcd(ic,klev+1,jb)-zsuduc(ic,jb)) )
        ENDDO
      ENDIF

      IF (PRESENT(flxd_par_sfc) .OR. lcomp_fractions) THEN ! compute photosynthetically active parts of surface radiation
        DO ic = 1, icount
          zflxd_par(ic)  = zflxd_par(ic)  + frc_par*bnd_wght(jb)*( zcloud(ic)*zbbfd(ic,klev+1,jb) &
            &                                                    + zclear(ic)*zbbcd(ic,klev+1,jb) )
        ENDDO
      ENDIF

      IF (lcomp_fractions) THEN
        ! VIS, NIR and PAR fractions of bands
        zfvis = bnd_wght(jb)*frc_vis(jb)
        zfnir = bnd_wght(jb)*frc_nir(jb)
        zfpar = bnd_wght(jb)*frc_par

        DO ic=1,icount
          zflxn_vis(ic)  = zflxn_vis(ic) + zfvis*( &
            &                 zcloud(ic)*zbbfd(ic,klev+1,jb)     &
            &               + zclear(ic)*zbbcd(ic,klev+1,jb)     &
            &               - zcloud(ic)*zbbfu(ic,klev+1,jb)     &
            &               - zclear(ic)*zbbcu(ic,klev+1,jb)     )
          zflxn(ic) = zflxd_sw(ic,klev+1) - zflxu_sw(ic,klev+1)
          zflxd_vis(ic)  = zflxd_vis(ic) + zfvis*(               &
            &                 zcloud(ic)*zbbfd(ic,klev+1,jb)     &
            &               + zclear(ic)*zbbcd(ic,klev+1,jb)     )
          zflxd_nir(ic)  = zflxd_nir(ic) + zfnir*(               &
            &                 zcloud(ic)*zbbfd(ic,klev+1,jb)     &
            &               + zclear(ic)*zbbcd(ic,klev+1,jb)     )
          zflxd_vis_diff(ic) = zflxd_vis_diff(ic) + zfvis*(      &
            &                 zcloud(ic)*zbbfd(ic,klev+1,jb)     &
            &               + zclear(ic)*zbbcd(ic,klev+1,jb)     &
            &               - zcloud(ic)*zsudu(ic,jb)            &
            &               - zclear(ic)*zsuduc(ic,jb))
          zflxd_nir_diff(ic) = zflxd_nir_diff(ic) + zfnir*(      &
            &                 zcloud(ic)*zbbfd(ic,klev+1,jb)     &
            &               + zclear(ic)*zbbcd(ic,klev+1,jb)     &
            &               - zcloud(ic)*zsudu(ic,jb)            &
            &               - zclear(ic)*zsuduc(ic,jb))
          zflxd_par_diff(ic) = zflxd_par_diff(ic) + zfpar*(      &
            &                 zcloud(ic)*zbbfd(ic,klev+1,jb)     &
            &               + zclear(ic)*zbbcd(ic,klev+1,jb)     &
            &               - zcloud(ic)*zsudu(ic,jb)            &
            &               - zclear(ic)*zsuduc(ic,jb))
        END DO
      END IF

    END DO ! jb

      DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
!$NEC ivdep
        DO ic = 1, icount
          jl = idx(ic)
          flxu_sw(jl,jk)     = zflxu_sw(ic,jk)
          flxd_sw(jl,jk)     = zflxd_sw(ic,jk)
          flxu_sw_clr(jl,jk) = zflxu_sw_clr(ic,jk)
          flxd_sw_clr(jl,jk) = zflxd_sw_clr(ic,jk)
        ENDDO
      ENDDO


      IF (PRESENT(flxd_dff_sfc)) THEN ! compute diffuse parts of surface radiation, including simple empirical parameterization
        DO ic = 1, icount             ! for 3D reflected radiation on low-level clouds
          jl = idx(ic)
          flxd_dff_sfc(jl) = MIN(zflxd_sw(ic,klev+1), &
                             zflxd_diff(ic) + tune_difrad_3dcont*zflxd_sw(ic,klev+1)*zcloud_ll(ic)*(1._wp-zcloud_ll(ic))**2)
   !       flxd_dff_sfc(jl) = zflxd_diff(ic)
        ENDDO
      ENDIF

      IF (PRESENT(flxd_par_sfc)) THEN ! compute photosynthetically active parts of surface radiation
        DO ic = 1, icount
          jl = idx(ic)
          flxd_par_sfc(jl) = zflxd_par(ic)
        ENDDO
      ENDIF

      IF (lcomp_fractions) THEN
        DO ic=1,icount
          jl = idx(ic)
          total = zflxn(ic) + zeps
          vis_frc_sfc(jl) = zflxn_vis(ic) / total
          total = zflxd_nir(ic) + zeps
          nir_dff_frc_sfc(jl) = zflxd_nir_diff(ic) / total
          total = zflxd_vis(ic) + zeps
          vis_dff_frc_sfc(jl) = zflxd_vis_diff(ic) / total
          total = zflxd_par(ic) + zeps
          par_dff_frc_sfc(jl) = zflxd_par_diff(ic) / total
        END DO
      END IF


  END SUBROUTINE srtm_srtm_224gp
  !-----------------------------------------------------------


  !-----------------------------------------------------------
  SUBROUTINE srtm_setcoef                                   &
    !   input
    & ( icount, idx,  kbdim    , klev    ,         &
    &   pavel   , ptavel   ,                                &
    &   pcoldry , pwkl     ,                                &
    !   output
    &   klaytrop, klayswtch, klaylow ,                      &
    &   pco2mult, pcolch4  , pcolco2 , pcolh2o , pcolmol  , &
    &   pcoln2o , pcolo2   , pcolo3  ,                      &
    &   pforfac , pforfrac , kindfor ,                      &
    &   pselffac, pselffrac, kindself,                      &
    &   pfac00  , pfac01   , pfac10  , pfac11  ,            &
    &   kjp     , kjt      , kjt1                          &
    !   input
    !&   zrmu0 )
    & )

    !     J. Delamere, AER, Inc. (version 2.5, 02/04/01)

    !     Modifications:
    !     JJMorcrette 030224   rewritten / adapted to ECMWF F90 system
    !        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC

    !     Purpose:  For a given atmosphere, calculate the indices and
    !     fractions related to the pressure and temperature interpolations.

    !-- Input arguments

    INTEGER(i4),INTENT(in)  :: icount
    INTEGER,  INTENT(in)    :: kbdim
    INTEGER(i4),INTENT(in)  :: idx(kbdim)
    INTEGER,  INTENT(in)    :: klev
    REAL(wp), INTENT(in)    :: pavel(kbdim,klev)
    REAL(wp), INTENT(in)    :: ptavel(kbdim,klev)
    REAL(wp), INTENT(in)    :: pcoldry(kbdim,klev)
    REAL(wp), INTENT(in)    :: pwkl(kbdim,jpinpx,klev)
    ! REAL(wp), INTENT(in)    :: zrmu0(kbdim)

    !-- Output arguments

    INTEGER,  INTENT(out)   :: klaytrop(kbdim)
    INTEGER,  INTENT(out)   :: klayswtch(kbdim)
    INTEGER,  INTENT(out)   :: klaylow(kbdim)
    REAL(wp), INTENT(out)   :: pco2mult(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcolch4(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcolco2(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcolh2o(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcolmol(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcoln2o(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcolo2(kbdim,klev)
    REAL(wp), INTENT(out)   :: pcolo3(kbdim,klev)
    REAL(wp), INTENT(out)   :: pforfac(kbdim,klev)
    REAL(wp), INTENT(out)   :: pforfrac(kbdim,klev)
    INTEGER,  INTENT(out)   :: kindfor(kbdim,klev)
    REAL(wp), INTENT(out)   :: pselffac(kbdim,klev)
    REAL(wp), INTENT(out)   :: pselffrac(kbdim,klev)
    INTEGER,  INTENT(out)   :: kindself(kbdim,klev)
    REAL(wp), INTENT(out)   :: pfac00(kbdim,klev)
    REAL(wp), INTENT(out)   :: pfac01(kbdim,klev)
    REAL(wp), INTENT(out)   :: pfac10(kbdim,klev)
    REAL(wp), INTENT(out)   :: pfac11(kbdim,klev)
    INTEGER,  INTENT(out)   :: kjp(kbdim,klev)
    INTEGER,  INTENT(out)   :: kjt(kbdim,klev)
    INTEGER,  INTENT(out)   :: kjt1(kbdim,klev)

    !-- local integers

    INTEGER :: i_nlayers, jk
    INTEGER :: jp1

    !-- local reals
    REAL(wp) :: z_stpfac
    REAL(wp) :: z_fp, z_ft, z_ft1, z_water, z_scalefac
    REAL(wp) :: z_factor, z_co2reg, z_compfp

    REAL(wp) :: rkindfor, rkindself
    REAL(wp) :: z_plog(kbdim),z_exptavel(kbdim),z_ptaveli(kbdim)
    INTEGER(i4) :: ic

    z_stpfac = 296._wp/1013._wp
    i_nlayers = klev

    DO ic = 1,icount
       klayswtch(ic) = 0
       klaytrop(ic)  = 0
       klaylow(ic)   = 0
    ENDDO

    DO jk = 1, i_nlayers

      DO ic = 1,icount
         z_plog(ic) = LOG(pavel(ic,jk))
         z_ptaveli(ic) = 1.0_wp/ptavel(ic,jk)
         z_exptavel(ic) = z_ptaveli(ic)*EXP(-1919.4_wp*z_ptaveli(ic))/8.7604e-4_wp
      END DO
      DO ic = 1,icount

        ! Find the two reference pressures on either side of the
        ! layer pressure.  Store them in JP and JP1.  Store in FP the
        ! fraction of the difference (in ln(pressure)) between these
        ! two values that the layer pressure lies.
        kjp(ic,jk) = INT(36._wp - 5._wp*(z_plog(ic)+0.04_wp))
        IF (kjp(ic,jk) < 1) THEN
          kjp(ic,jk) = 1
        ELSEIF (kjp(ic,jk) > 58) THEN
          kjp(ic,jk) = 58
        ENDIF
        jp1 = kjp(ic,jk) + 1
        z_fp = 5._wp * (preflog(kjp(ic,jk)) - z_plog(ic))

        !   Determine, for each reference pressure (JP and JP1), which
        !   reference temperature (these are different for each
        !   reference pressure) is nearest the layer temperature but does
        !   not exceed it.  Store these indices in JT and JT1, resp.
        !   Store in FT (resp. FT1) the fraction of the way between JT
        !   (JT1) and the next highest reference temperature that the
        !   layer temperature falls.
        kjt(ic,jk) = INT(3._wp + (ptavel(ic,jk)-tref(kjp(ic,jk)))/15._wp)
        IF (kjt(ic,jk) < 1) THEN
          kjt(ic,jk) = 1
        ELSEIF (kjt(ic,jk) > 4) THEN
          kjt(ic,jk) = 4
        ENDIF
        z_ft = ((ptavel(ic,jk)-tref(kjp(ic,jk)))/15._wp) - REAL(kjt(ic,jk)-3,wp)
        kjt1(ic,jk) = INT(3._wp + (ptavel(ic,jk)-tref(jp1))/15._wp)
        IF (kjt1(ic,jk) < 1) THEN
          kjt1(ic,jk) = 1
        ELSEIF (kjt1(ic,jk) > 4) THEN
          kjt1(ic,jk) = 4
        ENDIF
        z_ft1 = ((ptavel(ic,jk)-tref(jp1))/15._wp) - REAL(kjt1(ic,jk)-3,wp)

        z_water = pwkl(ic,1,jk)/pcoldry(ic,jk)
        z_scalefac = pavel(ic,jk) * z_stpfac * z_ptaveli(ic)


        !        If the pressure is less than ~100mb, perform a different
        !        set of species interpolations.

        IF (z_plog(ic) <= 4.56_wp) go to 5300

        klaytrop(ic) =  klaytrop(ic) + 1
        IF (z_plog(ic) >= 6.62_wp) klaylow(ic) = klaylow(ic) + 1

        !        Set up factors needed to separately include the water vapor
        !        foreign-continuum in the calculation of absorption coefficient.

        pforfac(ic,jk) = z_scalefac / (1._wp+z_water)
        z_factor = (332.0_wp-ptavel(ic,jk))/36.0_wp
        rkindfor = MIN(2.0_wp, MAX(1.0_wp, AINT(z_factor)))
        kindfor(ic,jk) = INT(rkindfor)
        pforfrac(ic,jk) = z_factor - rkindfor

        !        Set up factors needed to separately include the water vapor
        !        self-continuum in the calculation of absorption coefficient.

        pselffac(ic,jk) = z_water * pforfac(ic,jk)
        z_factor = (ptavel(ic,jk)-188.0_wp)/7.2_wp

        rkindself = MIN(9.0_wp, MAX(1.0_wp, AINT(z_factor)-7.0_wp))
        kindself(ic,jk) = INT(rkindself)
        pselffrac(ic,jk) = z_factor - rkindself - 7.0_wp

        !        Calculate needed column amounts.
        pcolh2o(ic,jk) = 1.e-20_wp * pwkl(ic,1,jk)
        pcolco2(ic,jk) = 1.e-20_wp * pwkl(ic,2,jk)
        pcolo3(ic,jk)  = 1.e-20_wp * pwkl(ic,3,jk)
        pcoln2o(ic,jk) = 1.e-20_wp * pwkl(ic,4,jk)
        pcolch4(ic,jk) = 1.e-20_wp * pwkl(ic,6,jk)
        pcolo2(ic,jk)  = 1.e-20_wp * pwkl(ic,7,jk)
        pcolmol(ic,jk) = 1.e-20_wp * pcoldry(ic,jk) + pcolh2o(ic,jk)
        IF (pcolco2(ic,jk) == 0._wp) pcolco2(ic,jk) = 1.e-32_wp * pcoldry(ic,jk)
        IF (pcoln2o(ic,jk) == 0._wp) pcoln2o(ic,jk) = 1.e-32_wp * pcoldry(ic,jk)
        IF (pcolch4(ic,jk) == 0._wp) pcolch4(ic,jk) = 1.e-32_wp * pcoldry(ic,jk)
        IF (pcolo2(ic,jk)  == 0._wp) pcolo2(ic,jk)  = 1.e-32_wp * pcoldry(ic,jk)
        !        Using E = 1334.2 cm-1.
        z_co2reg = 3.55e-24_wp * pcoldry(ic,jk)
        pco2mult(ic,jk)= (pcolco2(ic,jk) - z_co2reg) * &
          & 272.63_wp * z_exptavel(ic)

        go to 5400

        !        Above LAYTROP.
5300    CONTINUE

          !        Set up factors needed to separately include the water vapor
          !        foreign-continuum in the calculation of absorption coefficient.

          pforfac(ic,jk) = z_scalefac / (1._wp+z_water)
          z_factor = (ptavel(ic,jk)-188.0_wp)/36.0_wp
          kindfor(ic,jk) = 3
          pforfrac(ic,jk) = z_factor - 1.0_wp

          !        Calculate needed column amounts.

          pcolh2o(ic,jk) = 1.e-20_wp * pwkl(ic,1,jk)
          pcolco2(ic,jk) = 1.e-20_wp * pwkl(ic,2,jk)
          pcolo3(ic,jk)  = 1.e-20_wp * pwkl(ic,3,jk)
          pcoln2o(ic,jk) = 1.e-20_wp * pwkl(ic,4,jk)
          pcolch4(ic,jk) = 1.e-20_wp * pwkl(ic,6,jk)
          pcolo2(ic,jk)  = 1.e-20_wp * pwkl(ic,7,jk)
          pcolmol(ic,jk) = 1.e-20_wp * pcoldry(ic,jk) + pcolh2o(ic,jk)
          IF (pcolco2(ic,jk) == 0._wp) pcolco2(ic,jk) = 1.e-32_wp * pcoldry(ic,jk)
          IF (pcoln2o(ic,jk) == 0._wp) pcoln2o(ic,jk) = 1.e-32_wp * pcoldry(ic,jk)
          IF (pcolch4(ic,jk) == 0._wp) pcolch4(ic,jk) = 1.e-32_wp * pcoldry(ic,jk)
          IF (pcolo2(ic,jk)  == 0._wp) pcolo2(ic,jk)  = 1.e-32_wp * pcoldry(ic,jk)

          z_co2reg = 3.55e-24_wp * pcoldry(ic,jk)
          pco2mult(ic,jk)= (pcolco2(ic,jk) - z_co2reg) * &
            & 272.63_wp * z_exptavel(ic)

          pselffac(ic,jk) =0.0_wp
          pselffrac(ic,jk)=0.0_wp
          kindself(ic,jk) = 0

5400      CONTINUE

          !        We have now isolated the layer ln pressure and temperature,
          !        between two reference pressures and two reference temperatures
          !        (for each reference pressure).  We multiply the pressure
          !        fraction FP with the appropriate temperature fractions to get
          !        the factors that will be needed for the interpolation that yields
          !        the optical depths (performed in routines TAUGBn for band n).

          z_compfp = 1._wp - z_fp
          pfac10(ic,jk) = z_compfp * z_ft
          pfac00(ic,jk) = z_compfp * (1._wp - z_ft)
          pfac11(ic,jk) = z_fp * z_ft1
          pfac01(ic,jk) = z_fp * (1._wp - z_ft1)

      ENDDO
    ENDDO

    IF ( MAXVAL(klaytrop(1:icount)) == klev ) THEN
      CALL finish( 'mo_srtm:srtm_setcoef', 'Uppermost model layer too low for RRTM.' )
    ENDIF

    !-----------------------------------------------------------------------
  END SUBROUTINE srtm_setcoef
  !-----------------------------------------------------------------------

  SUBROUTINE srtm_spcvrt                                               &
    !   input
    & ( icount,   kbdim    , klev      , ksw,                          &
    &   poneminus,                                                     &
    &   palbd   , palbp    , zfrcl     ,                               &
    &   ptauc   , pasyc    , pomgc     ,                               &
    &   ptaua   , pasya    , pomga     ,                               &
    &   zrmu0   , klaytrop ,                                           &
    &   pcolch4 , pcolco2  , pcolh2o   , pcolmol  ,                    &
    &   pcolo2  , pcolo3   ,                                           &
    &   pforfac , pforfrac , kindfor   ,                               &
    &   pselffac, pselffrac, kindself  ,                               &
    &   pfac00  , pfac01   , pfac10    , pfac11   ,                    &
    &   kjp     , kjt      , kjt1      ,                               &
    !   output
    &   pbbfd   , pbbfu    , pbbcd     , pbbcu    ,                    &
    &   psudu   , psuduc )

    !**** *SRTM_SPCVRT* - SPECTRAL LOOP TO COMPUTE THE SHORTWAVE RADIATION FLUXES.

    !     PURPOSE.
    !     --------

    !          THIS ROUTINE COMPUTES THE TWO-STREAM METHOD OF BARKER

    !**   INTERFACE.
    !     ----------

    !          *SRTM_SPCVRT* IS CALLED FROM *SRTM_SRTM_224GP*

    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !     EXTERNALS.
    !     ----------

    !          *SWVRTQDR*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION
    !     AUTHOR.
    !     -------
    !        from Howard Barker
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 03-02-27
    !        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !        JJMorcrette 20070614 bug-fix for solar duration
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !        M.Puetz    20-Apr-2010 manual gather/scatter for prmu0 > 0, index passed to subroutines
    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    !-- Input arguments

    INTEGER(i4), INTENT(in) :: icount
    INTEGER,  INTENT(in)    :: ksw
    INTEGER,  INTENT(in)    :: kbdim
    INTEGER,  INTENT(in)    :: klev
    REAL(wp), INTENT(in)    :: poneminus(kbdim)
    REAL(wp), INTENT(in)    :: palbd(kbdim,ksw)
    REAL(wp), INTENT(in)    :: palbp(kbdim,ksw)
    REAL(wp), INTENT(in)    :: zfrcl(kbdim,klev)     ! bottom to top
    REAL(wp), INTENT(in)    :: ptauc(kbdim,klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pasyc(kbdim,klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pomgc(kbdim,klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: ptaua(kbdim,klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pasya(kbdim,klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: pomga(kbdim,klev,ksw) ! bottom to top
    REAL(wp), INTENT(in)    :: zrmu0(kbdim)
    INTEGER,  INTENT(in)    :: klaytrop(kbdim)
    REAL(wp), INTENT(in)    :: pcolch4(kbdim,klev)
    REAL(wp), INTENT(in)    :: pcolco2(kbdim,klev)
    REAL(wp), INTENT(in)    :: pcolh2o(kbdim,klev)
    REAL(wp), INTENT(in)    :: pcolmol(kbdim,klev)
    REAL(wp), INTENT(in)    :: pcolo2(kbdim,klev)
    REAL(wp), INTENT(in)    :: pcolo3(kbdim,klev)
    REAL(wp), INTENT(in)    :: pforfac(kbdim,klev)
    REAL(wp), INTENT(in)    :: pforfrac(kbdim,klev)
    INTEGER,  INTENT(in)    :: kindfor(kbdim,klev)
    REAL(wp), INTENT(in)    :: pselffac(kbdim,klev)
    REAL(wp), INTENT(in)    :: pselffrac(kbdim,klev)
    INTEGER,  INTENT(in)    :: kindself(kbdim,klev)
    REAL(wp), INTENT(in)    :: pfac00(kbdim,klev)
    REAL(wp), INTENT(in)    :: pfac01(kbdim,klev)
    REAL(wp), INTENT(in)    :: pfac10(kbdim,klev)
    REAL(wp), INTENT(in)    :: pfac11(kbdim,klev)
    INTEGER,  INTENT(in)    :: kjp(kbdim,klev)
    INTEGER,  INTENT(in)    :: kjt(kbdim,klev)
    INTEGER,  INTENT(in)    :: kjt1(kbdim,klev)

    !-- Input and output arguments

    REAL(wp) ,INTENT(inout) :: pbbfd(kbdim,klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: pbbfu(kbdim,klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: pbbcd(kbdim,klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: pbbcu(kbdim,klev+1,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: psudu(kbdim,jpb1:jpb2)
    REAL(wp) ,INTENT(inout) :: psuduc(kbdim,jpb1:jpb2)

    !-- Local variables

    LOGICAL :: llrtchk(kbdim,klev)

    REAL(wp) ::  zclear  , zcloud
    REAL(wp) ::                                                           &
      &   zdbt(kbdim,klev+1)                                              &
      & , zgcc(kbdim,klev)   , zgco(kbdim,klev)                           &
      & , zomcc(kbdim,klev)  , zomco(kbdim,klev)                          &
      & , zrdnd(kbdim,klev+1), zrdndc(kbdim,klev+1)                       &
      & , zref(kbdim,klev+1) , zrefc(kbdim,klev+1) , zrefo(kbdim,klev+1)  &
      & , zrefd(kbdim,klev+1), zrefdc(kbdim,klev+1), zrefdo(kbdim,klev+1) &
      & , zrup(kbdim,klev+1) , zrupd(kbdim,klev+1)                        &
      & , zrupc(kbdim,klev+1), zrupdc(kbdim,klev+1)                       &
      & , ztauc(kbdim,klev)  , ztauo(kbdim,klev)                          &
      & , ztdbt(kbdim,klev+1)                                             &
      & , ztra(kbdim,klev+1) , ztrac(kbdim,klev+1) , ztrao(kbdim,klev+1)  &
      & , ztrad(kbdim,klev+1), ztradc(kbdim,klev+1), ztrado(kbdim,klev+1)
    REAL(wp) ::                                                           &
      & zdbtc(kbdim,klev+1), ztdbtc(kbdim,klev+1)  , zincflx(kbdim)       &
      & , zincf14(kbdim,14)  , zinctot(kbdim)

    INTEGER :: ib1, ib2, ibm, igt, ikl, jb, jg, jk, i_kmodts
    REAL(wp) :: zdbtmc, zdbtmo, zf, zincflux(kbdim), zwf

    !-- Output of SRTM_TAUMOLn routines

    REAL(wp) :: ztaug(kbdim,klev,16), ztaur(kbdim,klev,16), zsflxzen(kbdim,16)

    !-- Output of SRTM_VRTQDR routine
    REAL(wp) ::                                            &
      &   zcd(kbdim,klev+1), zcu(kbdim,klev+1) &
      & , zfd(kbdim,klev+1), zfu(kbdim,klev+1)

    REAL(wp) :: zrmu0i(kbdim)
    INTEGER(i4) :: ic

    !     ------------------------------------------------------------------

    !-- Two-stream model 1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
    ! KMODTS is set in SWREFTRA
    !NDBUG=4
    ibm = -1
    igt = -1
    ib1=jpb1
    ib2=jpb2

!IBM* ASSERT(NODEPS)
    DO ic = 1,icount
       zincflux(ic)=0.0_wp
       zinctot(ic) =0.0_wp
       zrmu0i(ic)  = 1._wp / zrmu0(ic)
    ENDDO

    jb=ib1-1
    DO jb = ib1, ib2
!IBM* ASSERT(NODEPS)
      DO ic = 1,icount
        ibm = jb-15
        igt = ngc(ibm)
        zincf14(ic,ibm)=0.0_wp
      ENDDO

      !-- for each band, computes the gaseous and Rayleigh optical thickness
      !  for all g-points within the band
      IF (jb == 16) THEN
        CALL srtm_taumol16                                                          &
          & ( icount, kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolch4 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                      &
          & )

      ELSEIF (jb == 17) THEN
        CALL srtm_taumol17                                                          &
          & ( icount, kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 18) THEN
        CALL srtm_taumol18                                                          &
          & ( icount,   kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolch4 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 19) THEN
        CALL srtm_taumol19                                                          &
          & ( icount,   kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                  &
          & )

      ELSEIF (jb == 20) THEN
        CALL srtm_taumol20                                                          &
          & ( icount,  kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     ,                                        &
          &   pcolh2o , pcolch4 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                      &
          & )

      ELSEIF (jb == 21) THEN
        CALL srtm_taumol21                                                          &
          & ( icount,  kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                      &
          & )

      ELSEIF (jb == 22) THEN
        CALL srtm_taumol22                                                          &
          & ( icount,  kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolmol , pcolo2   ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 23) THEN
        CALL srtm_taumol23                                                          &
          & ( icount,   kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     ,                                        &
          &   pcolh2o , pcolmol ,                                                   &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 24) THEN
        CALL srtm_taumol24                                                          &
          & ( icount,   kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     , poneminus,                             &
          &   pcolh2o , pcolmol , pcolo2   , pcolo3   ,                             &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ELSEIF (jb == 25) THEN
        !--- visible 16000-22650 cm-1   0.4415 - 0.6250 um
        CALL srtm_taumol25                              &
          & ( icount,   kbdim   , klev     ,            &
          &   pfac00  , pfac01  , pfac10   , pfac11   , &
          &   kjp     , kjt     , kjt1     ,            &
          &   pcolh2o , pcolmol , pcolo3   ,            &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur          &
          & )

      ELSEIF (jb == 26) THEN
        !--- UV-A 22650-29000 cm-1   0.3448 - 0.4415 um
        CALL srtm_taumol26                              &
          & ( icount,  kbdim   , klev  ,       &
          &   pcolmol ,                                 &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur          &
          & )

      ELSEIF (jb == 27) THEN
        !--- UV-B 29000-38000 cm-1   0.2632 - 0.3448 um
        CALL srtm_taumol27                              &
          & ( icount,  kbdim   , klev     ,    &
          &   pfac00  , pfac01  , pfac10   , pfac11   , &
          &   kjp     , kjt     , kjt1     ,            &
          &   pcolmol , pcolo3  ,                       &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur                 &
          & )

      ELSEIF (jb == 28) THEN
        !--- UV-C 38000-50000 cm-1   0.2000 - 0.2632 um
        CALL srtm_taumol28 &
          & ( icount,   kbdim   , klev     ,   &
          &   pfac00  , pfac01  , pfac10   , pfac11   , &
          &   kjp     , kjt     , kjt1     , poneminus, &
          &   pcolmol , pcolo2  , pcolo3   ,            &
          &   klaytrop,                                 &
          &   zsflxzen, ztaug   , ztaur                 &
          & )

      ELSEIF (jb == 29) THEN
        CALL srtm_taumol29                                                          &
          & ( icount,  kbdim   , klev     ,                                        &
          &   pfac00  , pfac01  , pfac10   , pfac11   ,                             &
          &   kjp     , kjt     , kjt1     ,                                        &
          &   pcolh2o , pcolco2 , pcolmol  ,                                        &
          &   klaytrop, pselffac, pselffrac, kindself , pforfac, pforfrac, kindfor, &
          &   zsflxzen, ztaug   , ztaur                                     &
          & )

      ENDIF

      DO jg=1,igt

!IBM* ASSERT(NODEPS)
        DO ic = 1,icount

          zincflx(ic)     =zsflxzen(ic,jg)*zrmu0(ic)
          zincflux(ic)    =zincflux(ic)+zsflxzen(ic,jg)*zrmu0(ic)
          zinctot(ic)     =zinctot(ic)+zsflxzen(ic,jg)
          zincf14(ic,ibm) =zincf14(ic,ibm)+zsflxzen(ic,jg)
            !-- CALL to compute layer reflectances and transmittances for direct
            !  and diffuse sources, first clear then cloudy.
            !   Use direct/parallel albedo for direct radiation and diffuse albedo
            !   otherwise.

            ! ZREFC  direct albedo for clear
            ! ZREFO  direct albedo for cloud
            ! ZREFDC diffuse albedo for clear
            ! ZREFDO diffuse albedo for cloud
            ! ZTRAC  direct transmittance for clear
            ! ZTRAO  direct transmittance for cloudy
            ! ZTRADC diffuse transmittance for clear
            ! ZTRADO diffuse transmittance for cloudy

            ! ZREF   direct reflectance
            ! ZREFD  diffuse reflectance
            ! ZTRA   direct transmittance
            ! ZTRAD  diffuse transmittance

            ! ZDBTC  clear direct beam transmittance
            ! ZDBTO  cloudy direct beam transmittance
            ! ZDBT   layer mean direct beam transmittance
            ! ZTDBT  total direct beam transmittance at levels

            !-- clear-sky
            !----- TOA direct beam
            ztdbtc(ic,1)=1._wp
            !----- surface values
            zdbtc(ic,klev+1) =0.0_wp
            ztrac(ic,klev+1) =0.0_wp
            ztradc(ic,klev+1)=0.0_wp
            zrefc(ic,klev+1) =palbp(ic,ibm)
            zrefdc(ic,klev+1)=palbd(ic,ibm)
            zrupc(ic,klev+1) =palbp(ic,ibm)
            zrupdc(ic,klev+1)=palbd(ic,ibm)

            !-- total sky
            !----- TOA direct beam
            ztdbt(ic,1)=1._wp
            !----- surface values
            zdbt(ic,klev+1) =0.0_wp
            ztra(ic,klev+1) =0.0_wp
            ztrad(ic,klev+1)=0.0_wp
            zref(ic,klev+1) =palbp(ic,ibm)
            zrefd(ic,klev+1)=palbd(ic,ibm)
            zrup(ic,klev+1) =palbp(ic,ibm)
            zrupd(ic,klev+1)=palbd(ic,ibm)
          ENDDO

!IBM* NOUNROLL
          DO jk=1,klev
            ikl=klev+1-jk
            !-- NB: a two-stream calculations from top to bottom, but RRTM_SW quantities
            !       are given bottom to top (argh!)
            !       Inputs for clouds and aerosols are bottom to top as inputs
            DO ic = 1,icount

              !-- clear-sky optical parameters
              ! llrtchk(jl,jk)=.TRUE. here always true

              !-- clear-sky optical parameters including aerosols
              ztauc(ic,jk) = ztaur(ic,ikl,jg) + ztaug(ic,ikl,jg) + ptaua(ic,ikl,ibm)
              ztauo(ic,jk) = ztauc(ic,jk) + ptauc(ic,ikl,ibm)

              zomcc(ic,jk) = ztaur(ic,ikl,jg) + ptaua(ic,ikl,ibm)*pomga(ic,ikl,ibm)
              zomco(ic,jk) = zomcc(ic,jk) + ptauc(ic,ikl,ibm)*pomgc(ic,ikl,ibm)

              zgcc(ic,jk)  = pasya(ic,ikl,ibm)*pomga(ic,ikl,ibm)*ptaua(ic,ikl,ibm)
              zgco(ic,jk)  = zgcc(ic,jk) + ptauc(ic,ikl,ibm)*pomgc(ic,ikl,ibm)*pasyc(ic,ikl,ibm)

            ENDDO
            zgcc(1:icount,jk)  = zgcc(1:icount,jk)  / zomcc(1:icount,jk)
            zgco(1:icount,jk)  = zgco(1:icount,jk)  / zomco(1:icount,jk)
            zomcc(1:icount,jk) = zomcc(1:icount,jk) / ztauc(1:icount,jk)
            zomco(1:icount,jk) = zomco(1:icount,jk) / ztauo(1:icount,jk)

!IBM* NOVECTOR
            DO ic = 1,icount

              !-- Delta scaling for clear-sky / aerosol optical quantities
              zf=zgcc(ic,jk)*zgcc(ic,jk)
              zwf=zomcc(ic,jk)*zf
              ztauc(ic,jk)=MAX(0._wp,(1._wp-zwf)*ztauc(ic,jk))
              zomcc(ic,jk)=(zomcc(ic,jk)-zwf)/(1.0_wp-zwf)
              zgcc(ic,jk)=(zgcc(ic,jk)-zf)/(1.0_wp-zf)

              !-- Delta scaling for cloudy quantities
              zf=zgco(ic,jk)*zgco(ic,jk)
              zwf=zomco(ic,jk)*zf
              ztauo(ic,jk)=MAX(0._wp,(1._wp-zwf)*ztauo(ic,jk))
              zomco(ic,jk)=(zomco(ic,jk)-zwf)/(1._wp-zwf)
              zgco(ic,jk)=(zgco(ic,jk)-zf)/(1._wp-zf)
            ENDDO
        ENDDO

        ! mpuetz: llrtchk = TRUE always in this call to reftra()
        CALL srtm_reftra ( icount, kbdim, klev, i_kmodts ,&
          &   zgcc  , zrmu0, zrmu0i, ztauc , zomcc ,&
          &   zrefc  , zrefdc, ztrac, ztradc )

        DO jk=1,klev
          ikl=klev+1-jk
          DO ic = 1,icount
            llrtchk(ic,jk)=(zfrcl(ic,ikl) > repclc)
          END DO
        END DO

        CALL srtm_reftra ( icount, kbdim, klev, i_kmodts ,&
          &   zgco  , zrmu0, zrmu0i, ztauo , zomco ,&
          &   zrefo , zrefdo, ztrao, ztrado , llrtchk )

        DO jk=1,klev
          ikl=klev+1-jk
!IBM* ASSERT(NODEPS)
          DO ic = 1,icount

            !-- combine clear and cloudy contributions for total sky
            zclear   = 1.0_wp - zfrcl(ic,ikl)
            zcloud   = zfrcl(ic,ikl)

            zref(ic,jk) = zclear*zrefc(ic,jk) + zcloud*zrefo(ic,jk)
            zrefd(ic,jk)= zclear*zrefdc(ic,jk)+ zcloud*zrefdo(ic,jk)
            ztra(ic,jk) = zclear*ztrac(ic,jk) + zcloud*ztrao(ic,jk)
            ztrad(ic,jk)= zclear*ztradc(ic,jk)+ zcloud*ztrado(ic,jk)

            !-- direct beam transmittance

            zdbtmc     = EXP(-ztauc(ic,jk)*zrmu0i(ic))
            zdbtmo     = EXP(-ztauo(ic,jk)*zrmu0i(ic))
            zdbt(ic,jk)   = zclear*zdbtmc+zcloud*zdbtmo
            ztdbt(ic,jk+1)= zdbt(ic,jk)*ztdbt(ic,jk)

            !-- clear-sky
            zdbtc(ic,jk)   =zdbtmc
            ztdbtc(ic,jk+1)=zdbtc(ic,jk)*ztdbtc(ic,jk)
          ENDDO
        ENDDO


        !-- vertical quadrature producing clear-sky fluxes
        CALL srtm_vrtqdr ( icount, kbdim, klev, &
          &   zrefc, zrefdc, ztrac , ztradc ,&
          &   zdbtc, zrdndc, zrupc , zrupdc, ztdbtc ,&
          &   zcd  , zcu )

        !-- vertical quadrature producing cloudy fluxes

        CALL srtm_vrtqdr ( icount, kbdim, klev, &
          &   zref , zrefd , ztra , ztrad ,&
          &   zdbt , zrdnd , zrup , zrupd , ztdbt ,&
          &   zfd  , zfu )


        !-- up and down-welling fluxes at levels
        IF (icount > 0) THEN
!IBM* NOUNROLL
          DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
            DO ic = 1,icount
              !-- accumulation of spectral fluxes
              pbbfu(ic,jk,jb) = pbbfu(ic,jk,jb) + zincflx(ic)*zfu(ic,jk)
              pbbfd(ic,jk,jb) = pbbfd(ic,jk,jb) + zincflx(ic)*zfd(ic,jk)
              pbbcu(ic,jk,jb) = pbbcu(ic,jk,jb) + zincflx(ic)*zcu(ic,jk)
              pbbcd(ic,jk,jb) = pbbcd(ic,jk,jb) + zincflx(ic)*zcd(ic,jk)
            ENDDO
          ENDDO
        END IF
!IBM* ASSERT(NODEPS)
        DO ic = 1,icount
           psudu(ic,jb) = psudu(ic,jb)  +zincflx(ic)*ztdbt (ic,klev+1)
           psuduc(ic,jb)= psuduc(ic,jb) +zincflx(ic)*ztdbtc(ic,klev+1)
        ENDDO
        !
        ! -- direct flux
        !

      ENDDO
      !-- end loop on JG

    ENDDO
    !-- end loop on JB
    !------------------------------------------------------------------

  END SUBROUTINE srtm_spcvrt

!!!#ifdef RS6K
!!!@PROCESS HOT NOSTRICT
!!!#endif
  SUBROUTINE srtm_vrtqdr &
    &(icount, kbdim, klev , &
    & pref , prefd, ptra , ptrad,&
    & pdbt , prdnd, prup , prupd , ptdbt,&
    & pfd  , pfu  &
    & )

    !**** *SRTM_VRTQDR* - VERTICAL QUADRATURE

    !     PURPOSE.
    !     --------
    !
    !          THIS ROUTINE PERFORMS THE VERTICAL INTEGRATION
    !
    !     AUTHOR.
    !     -------
    !        from Howard Barker
    !        JEAN-JACQUES MORCRETTE  *ECMWF*
    !
    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 02-10-04
    !        M.Hamrud      01-Oct-2003 CY28 Cleaning
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !        M.Puetz    20-Apr-2010 Added IBM optimization directives, made indices 64 bit
    !              ---------

    INTEGER(i4),INTENT(in)    :: icount
    INTEGER,INTENT(in)        :: kbdim
    INTEGER,INTENT(in)        :: klev
    REAL(wp)   ,INTENT(in)    :: pref(kbdim,klev+1)
    REAL(wp)   ,INTENT(in)    :: prefd(kbdim,klev+1)
    REAL(wp)   ,INTENT(in)    :: ptra(kbdim,klev+1)
    REAL(wp)   ,INTENT(in)    :: ptrad(kbdim,klev+1)
    REAL(wp)   ,INTENT(in)    :: pdbt(kbdim,klev+1)
    REAL(wp)   ,INTENT(out)   :: prdnd(kbdim,klev+1)
    REAL(wp)   ,INTENT(inout) :: prup(kbdim,klev+1)
    REAL(wp)   ,INTENT(inout) :: prupd(kbdim,klev+1)
    REAL(wp)   ,INTENT(in)    :: ptdbt(kbdim,klev+1)
    REAL(wp)   ,INTENT(inout) :: pfd(kbdim,klev+1)
    REAL(wp)   ,INTENT(inout) :: pfu(kbdim,klev+1)
    !              ------------
    REAL(wp) :: ztdn(kbdim,klev+1)

    INTEGER(i4) :: ikp, ikx, jk, ic

    REAL(wp) :: zreflect

    !------------------------------------------------------------------
    ! PREF(JK)   direct reflectance
    ! PREFD(JK)  diffuse reflectance
    ! PTRA(JK)   direct transmittance
    ! PTRAD(JK)  diffuse transmittance
    ! PDBT(JK)   layer mean direct beam transmittance
    ! PTDBT(JK)  total direct beam transmittance at levels
    !
    !-- link lowest layer with surface


!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
    DO ic=1,icount
      zreflect=1.0_wp / (1.0_wp -prefd(ic,klev+1)*prefd(ic,klev))
      prup(ic,klev)=pref(ic,klev)+(ptrad(ic,klev)* &
           & ((ptra(ic,klev)-pdbt(ic,klev))*prefd(ic,klev+1)+ &
           & pdbt(ic,klev)*pref(ic,klev+1)))*zreflect
      prupd(ic,klev)=prefd(ic,klev)+ptrad(ic,klev)* &
           & ptrad(ic,klev)*prefd(ic,klev+1)*zreflect
    ENDDO

    !-- pass from bottom to top

!$NEC outerloop_unroll(8)
    DO jk=1,klev-1
      ikp=klev+1-jk
      ikx=ikp-1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO ic=1,icount
         zreflect=1.0_wp / (1.0_wp -prupd(ic,ikp)*prefd(ic,ikx))
         prup(ic,ikx)=pref(ic,ikx)+(ptrad(ic,ikx)* &
              & ((ptra(ic,ikx)-pdbt(ic,ikx))*prupd(ic,ikp)+ &
              & pdbt(ic,ikx)*prup(ic,ikp)))*zreflect
         prupd(ic,ikx)=prefd(ic,ikx)+ptrad(ic,ikx)* &
              & ptrad(ic,ikx)*prupd(ic,ikp)*zreflect
      ENDDO
    ENDDO

    !-- upper boundary conditions

!IBM* ASSERT(NODEPS)
    DO ic=1,icount
       ztdn(ic,1)=1.0_wp
       prdnd(ic,1)=0.0_wp
       ztdn(ic,2)=ptra(ic,1)
       prdnd(ic,2)=prefd(ic,1)
    ENDDO

    !-- pass from top to bottom
!$NEC outerloop_unroll(8)
    DO jk=2,klev
      ikp=jk+1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO ic=1,icount
         zreflect=1.0_wp / (1.0_wp -prefd(ic,jk)*prdnd(ic,jk))
         ztdn(ic,ikp)=ptdbt(ic,jk)*ptra(ic,jk)+ &
              & (ptrad(ic,jk)*((ztdn(ic,jk)-ptdbt(ic,jk))+ &
              & ptdbt(ic,jk)*pref(ic,jk)*prdnd(ic,jk))) * zreflect
         prdnd(ic,ikp)=prefd(ic,jk)+ptrad(ic,jk)*ptrad(ic,jk) &
              & *prdnd(ic,jk)*zreflect
      ENDDO
    ENDDO

    !-- up and down-welling fluxes at levels
    DO jk=1,klev+1
!IBM* ASSERT(NODEPS)
!IBM* NOVECTOR
      DO ic=1,icount
        zreflect=1.0_wp / (1.0_wp - prdnd(ic,jk)*prupd(ic,jk))
        pfu(ic,jk)=(ptdbt(ic,jk)*prup(ic,jk) + &
             & (ztdn(ic,jk)-ptdbt(ic,jk))*prupd(ic,jk))*zreflect
        pfd(ic,jk)=ptdbt(ic,jk) + (ztdn(ic,jk)-ptdbt(ic,jk)+ &
             & ptdbt(ic,jk)*prup(ic,jk)*prdnd(ic,jk))*zreflect
      ENDDO
    ENDDO

  END SUBROUTINE srtm_vrtqdr
  !------------------------------------------------------------------


  !------------------------------------------------------------------
!!!#ifdef RS6K
!!!@PROCESS HOT NOSTRICT
!!!#endif
  SUBROUTINE srtm_reftra &
       & ( icount, kbdim, klev  , kmodts, &
       &   pgg   , prmuz, prmuzi, ptau , pw, &
       &   pref  , prefd, ptra , ptrad , &
       &   ldrtchk &
       & )

    !**** *SRTM_REFTRA* - REFLECTIVITY AND TRANSMISSIVITY

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLEAR OR
    !     CLOUDY LAYER USING A CHOICE OF VARIOUS APPROXIMATIONS.

    !**   INTERFACE.
    !     ----------
    !          *SRTM_REFTRA* IS CALLED BY *SRTM_SPCVRT*

    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! INPUTS
    ! ------
    !      KMODTS  = 1 EDDINGTON (JOSEPH ET AL., 1976)
    !              = 2 PIFM (ZDUNKOWSKI ET AL., 1980)
    !              = 3 DISCRETE ORDINATES (LIOU, 1973)
    !      LDRTCHK = .T. IF CLOUDY
    !              = .F. IF CLEAR-SKY
    !      PGG     = ASSYMETRY FACTOR
    !      PRMUZ   = COSINE SOLAR ZENITH ANGLE
    !      PRMUZI  = INVERSE COSINE SOLAR ZENITH ANGLE
    !      PTAU    = OPTICAL THICKNESS
    !      PW      = SINGLE SCATTERING ALBEDO

    ! OUTPUTS
    ! -------
    !      PREF    : COLLIMATED BEAM REFLECTIVITY
    !      PREFD   : DIFFUSE BEAM REFLECTIVITY
    !      PTRA    : COLLIMATED BEAM TRANSMISSIVITY
    !      PTRAD   : DIFFUSE BEAM TRANSMISSIVITY

    !     METHOD.
    !     -------
    !          STANDARD DELTA-EDDINGTON, P.I.F.M., OR D.O.M. LAYER CALCULATIONS.

    !     EXTERNALS.
    !     ----------
    !          NONE

    !     REFERENCE.
    !     ----------

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 03-02-27
    !        M.Hamrud   01-Oct-2003      CY28 Cleaning
    !        Mike Iacono, AER, Mar 2004: bug fix
    !        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
    !        M. Puetz   20-Apr-2010 Gather/Scatter for better pipelining on scalar cpus
    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------
    INTEGER(i4),INTENT(in)  :: icount
    INTEGER,INTENT(in)        :: kbdim
    INTEGER,INTENT(in)        :: klev
    INTEGER,INTENT(out)       :: kmodts
    REAL(wp)   ,INTENT(in)    :: pgg(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: prmuz(kbdim)
    REAL(wp)   ,INTENT(in)    :: prmuzi(kbdim)
    REAL(wp)   ,INTENT(in)    :: ptau(kbdim,klev)
    REAL(wp)   ,INTENT(in)    :: pw(kbdim,klev)
    REAL(wp)   ,INTENT(out)   :: pref(kbdim,klev)
    REAL(wp)   ,INTENT(inout) :: prefd(kbdim,klev)
    REAL(wp)   ,INTENT(inout) :: ptra(kbdim,klev)
    REAL(wp)   ,INTENT(inout) :: ptrad(kbdim,klev)
    LOGICAL,INTENT(in),optional :: ldrtchk(kbdim,klev)
    !     ------------------------------------------------------------------
    INTEGER(i4) :: idxt(kbdim), idxf(kbdim), idxc(kbdim), idxn(kbdim)
    INTEGER(i4) :: jk, ic, jc, ict, icf, icc, icn

    REAL(wp) :: zgamma1(kbdim),zgamma2(kbdim),zgamma3(kbdim),zcrit(kbdim)
    REAL(wp) :: zrk(kbdim), zem1(kbdim), zem2(kbdim), zep1(kbdim), zep2(kbdim)
    REAL(wp) :: za, za1, za2, zemm
    REAL(wp) :: zbeta, zdend, zdenr, zdent
    REAL(wp) :: zg, zg3, zgamma4, zgt
    REAL(wp) :: zr1, zr2, zr3, zr4, zr5, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
    REAL(wp) :: zsr3, zt1, zt2, zt3 ! , zt4, zt5
    REAL(wp) :: zw, zwcrit
    REAL(wp) :: ztemp, zzz, zeps, zaux
    !------------------------------------------------------------------



    zeps = 1.e-20_wp

    zsr3=SQRT(3._wp)
    zwcrit=0.9995_wp
    kmodts=2
    ict = -1

    IF (.NOT.PRESENT(ldrtchk)) THEN
      ict = icount
      icf = 0
      DO ic=1,icount
        idxt(ic) = ic
      END DO
    END IF

!
    DO jk=1,klev
      IF (PRESENT(ldrtchk)) THEN
        ict = 0
        icf = 0
        IF (kmodts == 2) THEN ! optimization: check first if index lists need to be generated
          DO ic=1,icount
            IF (ldrtchk(ic,jk)) THEN
              ict = ict + 1
            END IF
          END DO
          IF (ict == 0) icf = icount
        END IF
        IF (kmodts /= 2 .OR. (ict < icount .AND. icf < icount)) THEN
          ict = 0
          DO ic=1,icount
            IF (ldrtchk(ic,jk)) THEN
              ict = ict + 1
              idxt(ict) = ic
            ELSE
              icf = icf + 1
              idxf(icf) = ic
            END IF
          END DO
        END IF

      END IF
       !-- GENERAL TWO-STREAM EXPRESSIONS

      IF (kmodts == 1) THEN
!IBM* ASSERT(NODEPS)
!$NEC ivdep
        DO jc = 1,ict
          ic = idxt(jc)

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)

          zg3= 3._wp * zg
          zgamma1(ic)= (7._wp - zw * (4._wp + zg3)) * 0.25_wp
          zgamma2(ic)=-(1._wp - zw * (4._wp - zg3)) * 0.25_wp
          zgamma3(ic)= (2._wp - zg3 * prmuz(ic) ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit(jc) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 2 .AND. ict == icount) THEN ! use direct addressing
        DO ic = 1,ict

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)

          zg3= 3._wp * zg
          zgamma1(ic)= (8._wp - zw * (5._wp + zg3)) * 0.25_wp
          zgamma2(ic)=  3._wp *(zw * (1._wp - zg )) * 0.25_wp
          zgamma3(ic)= (2._wp - zg3 * prmuz(ic) ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit(ic) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 2) THEN
!IBM* ASSERT(NODEPS)
!$NEC ivdep
!DIR$ IVDEP
        DO jc = 1,ict
          ic = idxt(jc)

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)

          zg3= 3._wp * zg
          zgamma1(ic)= (8._wp - zw * (5._wp + zg3)) * 0.25_wp
          zgamma2(ic)=  3._wp *(zw * (1._wp - zg )) * 0.25_wp
          zgamma3(ic)= (2._wp - zg3 * prmuz(ic) ) * 0.25_wp
          zzz=(1._wp - zg)**2
          zcrit(jc) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ELSEIF (kmodts == 3) THEN
!IBM* ASSERT(NODEPS)
!$NEC ivdep
        DO jc = 1,ict
          ic = idxt(jc)

          zw  =pw(ic,jk)
          zg  =pgg(ic,jk)

          zg3= 3._wp * zg
          zgamma1(ic)= zsr3 * (2._wp - zw * (1._wp + zg)) * 0.5_wp
          zgamma2(ic)= zsr3 * zw * (1._wp - zg ) * 0.5_wp
          zgamma3(ic)= (1._wp - zsr3 * zg * prmuz(ic) ) * 0.5_wp
          zzz=(1._wp - zg)**2
          zcrit(jc) = zw*zzz - zwcrit*(zzz - (1._wp - zw)*(zg **2))
        END DO
      ENDIF

      icc = 0
      icn = 0

      DO jc = 1,ict
        IF (zcrit(jc) >= 0._wp) icc = icc + 1
      END DO
      IF (icc == 0 .AND. ict == icount) icn = icount
      IF(icc < icount .AND. icn < icount) THEN
        icc = 0

        IF(ict == icount) THEN ! this case uses direct adressing so far
          DO jc = 1,ict
            IF (zcrit(jc) >= 0._wp) THEN
               icc = icc + 1
               idxc(icc) = jc
            ELSE
               icn = icn + 1
               idxn(icn) = jc
            END IF
          END DO
        ELSE ! further split existing index list
          DO jc = 1,ict
            ic = idxt(jc)
            IF (zcrit(jc) >= 0._wp) THEN
               icc = icc + 1
               idxc(icc) = ic
            ELSE
               icn = icn + 1
               idxn(icn) = ic
            END IF
          END DO
        ENDIF
      ENDIF

        !-- RECOMPUTE ORIGINAL S.S.A. TO TEST FOR CONSERVATIVE SOLUTION
        !   ZTEMP=(1._wp - ZG)**2
        !   ZWO= ZW*ZTEMP/ (ZTEMP - (1._wp - ZW)*(ZG **2))

        !       ZWO= ZW / (1._wp - (1._wp - ZW) * (ZG / (1._wp - ZG))**2)
        !       IF (ZWO >= ZWCRIT) THEN

       !!-- conservative scattering

        !-- Homogeneous reflectance and transmittance

      IF (icc == icount) THEN ! use direct addressing
        DO ic = 1,icc
          zem2(ic) = -MIN(ptau(ic,jk) * prmuzi(ic),500._wp)
        END DO
      ELSE
!IBM* ASSERT(NODEPS)
!$NEC ivdep
        DO jc = 1,icc
          ic = idxc(jc)
          zem2(jc) = -MIN(ptau(ic,jk) * prmuzi(ic),500._wp)
        END DO
      ENDIF

      zem2(1:icc) = EXP(zem2(1:icc))

      IF (icc == icount) THEN ! use direct addressing
!IBM* NOVECTOR
        DO ic = 1,icc
          za  = zgamma1(ic) * prmuz(ic)
          za1 = za - zgamma3(ic)
          zgt = zgamma1(ic) * ptau(ic,jk)

          ! collimated beam

          ztemp=1.0_wp/(1._wp + zgt)
          pref(ic,jk) = (zgt - za1 * (1._wp - zem2(ic))) *ztemp
          ptra(ic,jk) = 1._wp - pref(ic,jk)

          ! isotropic incidence

          prefd(ic,jk) = zgt *ztemp
          ptrad(ic,jk) = 1._wp - prefd(ic,jk)
        END DO
      ELSE
!IBM* NOVECTOR
!IBM* ASSERT(NODEPS)
!$NEC ivdep
!DIR$ IVDEP
        DO jc = 1,icc
          ic = idxc(jc)

          za  = zgamma1(ic) * prmuz(ic)
          za1 = za - zgamma3(ic)
          zgt = zgamma1(ic) * ptau(ic,jk)

          ! collimated beam

          ztemp=1.0_wp/(1._wp + zgt)
          pref(ic,jk) = (zgt - za1 * (1._wp - zem2(jc))) *ztemp
          ptra(ic,jk) = 1._wp - pref(ic,jk)

          ! isotropic incidence

          prefd(ic,jk) = zgt *ztemp
          ptrad(ic,jk) = 1._wp - prefd(ic,jk)
        END DO
      ENDIF

     !-- non-conservative scattering

     !-- Homogeneous reflectance and transmittance
      IF (icn == icount) THEN ! use direct addressing
!DIR$ IVDEP
        DO ic = 1,icn
          zzz = zgamma1(ic)**2 - zgamma2(ic)**2
          zrk(ic) = SQRT(MAX(zzz,replog))

          zep1(ic) = MIN(zrk(ic) * ptau(ic,jk), 500._wp)
          zep2(ic) = MIN(ptau(ic,jk) * prmuzi(ic),500._wp)
        END DO
      ELSE
!IBM* ASSERT(NODEPS)
!$NEC ivdep
!DIR$ IVDEP
        DO jc = 1,icn
          ic = idxn(jc)
          zzz = zgamma1(ic)**2 - zgamma2(ic)**2
          zrk(jc) = SQRT(MAX(zzz,replog))

          zep1(jc) = MIN(zrk(jc) * ptau(ic,jk), 500._wp)
          zep2(jc) = MIN(ptau(ic,jk) * prmuzi(ic),500._wp)
        END DO
      ENDIF

      zep1(1:icn) = EXP(zep1(1:icn))
      zep2(1:icn) = EXP(zep2(1:icn))
      zem1(1:icn) = 1.0_wp/zep1(1:icn)
      zem2(1:icn) = 1.0_wp/zep2(1:icn)

      IF (icn == icount) THEN ! use direct addressing
!$OMP SIMD
        DO ic = 1,icn
          zw  =pw(ic,jk)

          zgamma4 = 1._wp - zgamma3(ic)

          za1 = zgamma1(ic) * zgamma4 + zgamma2(ic) * zgamma3(ic)
          za2 = zgamma1(ic) * zgamma3(ic) + zgamma2(ic) * zgamma4

          zrp = zrk(ic) * prmuz(ic)
          zrp1 = 1._wp + zrp
          zrm1 = 1._wp - zrp
          zrk2 = 2._wp * zrk(ic)
          zrpp = 1._wp - zrp*zrp
          zrkg = zrk(ic) + zgamma1(ic)
          zr1  = zrm1 * (za2 + zrk(ic) * zgamma3(ic))
          zr2  = zrp1 * (za2 - zrk(ic) * zgamma3(ic))
          zr3  = zrk2 * (zgamma3(ic) - za2 * prmuz(ic) )
          zr4  = zrpp * zrkg
          zr5  = zrpp * (zrk(ic) - zgamma1(ic))
          zt1  = zrp1 * (za1 + zrk(ic) * zgamma4)
          zt2  = zrm1 * (za1 - zrk(ic) * zgamma4)
          zt3  = zrk2 * (zgamma4 + za1 * prmuz(ic) )
          ! zt4  = zr4
          ! zt5  = zr5
          ! GZ, 2015-06-08: another fix for potential division by zero
          zbeta = - zr5 / SIGN(MAX(zeps,ABS(zr4)),zr4)

          ! collimated beam

          ! GZ, 2014-07-03: provisional fix for potential division by zero
          zaux = zr4*zep1(ic) + zr5*zem1(ic)
          zdenr = SIGN(MAX(zeps,ABS(zaux)),zaux)
          pref(ic,jk) = zw  * (zr1*zep1(ic) - zr2*zem1(ic) - zr3*zem2(ic)) / zdenr

          ! zdent = zt4*zep1(ic) + zt5*zem1(ic)
          zdent = zdenr
          ptra(ic,jk) = zem2(ic) * &
            & (1._wp - zw  * (zt1*zep1(ic) - zt2*zem1(ic) - zt3*zep2(ic)) / zdent)

          ! diffuse beam

          zemm = zem1(ic)*zem1(ic)
          zdend = 1._wp / ( (1._wp - zbeta*zemm ) * zrkg)
          prefd(ic,jk) =  zgamma2(ic) * (1._wp - zemm) * zdend
          ptrad(ic,jk) =  zrk2*zem1(ic)*zdend

        END DO
      ELSE
!IBM* ASSERT(NODEPS)
!$NEC ivdep
!DIR$ IVDEP
        DO jc = 1,icn
          ic = idxn(jc)
          zw  =pw(ic,jk)

          zgamma4 = 1._wp - zgamma3(ic)

          za1 = zgamma1(ic) * zgamma4 + zgamma2(ic) * zgamma3(ic)
          za2 = zgamma1(ic) * zgamma3(ic) + zgamma2(ic) * zgamma4

          zrp = zrk(jc) * prmuz(ic)
          zrp1 = 1._wp + zrp
          zrm1 = 1._wp - zrp
          zrk2 = 2._wp * zrk(jc)
          zrpp = 1._wp - zrp*zrp
          zrkg = zrk(jc) + zgamma1(ic)
          zr1  = zrm1 * (za2 + zrk(jc) * zgamma3(ic))
          zr2  = zrp1 * (za2 - zrk(jc) * zgamma3(ic))
          zr3  = zrk2 * (zgamma3(ic) - za2 * prmuz(ic) )
          zr4  = zrpp * zrkg
          zr5  = zrpp * (zrk(jc) - zgamma1(ic))
          zt1  = zrp1 * (za1 + zrk(jc) * zgamma4)
          zt2  = zrm1 * (za1 - zrk(jc) * zgamma4)
          zt3  = zrk2 * (zgamma4 + za1 * prmuz(ic) )
          ! zt4  = zr4
          ! zt5  = zr5
          ! GZ, 2015-06-08: another fix for potential division by zero
          zbeta = - zr5 / SIGN(MAX(zeps,ABS(zr4)),zr4)

          ! collimated beam

          ! GZ, 2014-07-03: provisional fix for potential division by zero
          zaux  = zr4*zep1(jc) + zr5*zem1(jc)
          zdenr = SIGN(MAX(zeps,ABS(zaux)),zaux)
          pref(ic,jk) = zw  * (zr1*zep1(jc) - zr2*zem1(jc) - zr3*zem2(jc)) / zdenr

          ! zdent  = zt4*zep1(jc) + zt5*zem1(jc)
          zdent = zdenr
          ptra(ic,jk) = zem2(jc) * &
            & (1._wp - zw  * (zt1*zep1(jc) - zt2*zem1(jc) - zt3*zep2(jc)) / zdent)

          ! diffuse beam

          zemm = zem1(jc)*zem1(jc)
          zdend = 1._wp / ( (1._wp - zbeta*zemm ) * zrkg)
          prefd(ic,jk) =  zgamma2(ic) * (1._wp - zemm) * zdend
          ptrad(ic,jk) =  zrk2*zem1(jc)*zdend

        END DO
      ENDIF

      IF (icf == icount) THEN ! use direct addressing
!$OMP SIMD
        DO ic = 1,icf
          pref(ic,jk) =0.0_wp
          ptra(ic,jk) =1.0_wp
          prefd(ic,jk)=0.0_wp
          ptrad(ic,jk)=1.0_wp
        END DO
      ELSE IF (icf > 0) THEN
!IBM* ASSERT(NODEPS)
!$NEC ivdep
        DO jc = 1,icf
          ic=idxf(jc)
          pref(ic,jk) =0.0_wp
          ptra(ic,jk) =1.0_wp
          prefd(ic,jk)=0.0_wp
          ptrad(ic,jk)=1.0_wp
        END DO
      ENDIF
    ENDDO

  END SUBROUTINE srtm_reftra
  ! ------------------------------------------------------------------

END MODULE mo_srtm

