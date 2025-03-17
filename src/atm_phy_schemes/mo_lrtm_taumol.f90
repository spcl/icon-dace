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


!----------------------------
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

!----------------------------

module mo_lrtm_taumol

  ! ------- Modules -------

  use mo_kind,     only : wp
  use mo_lrtm_par, only : nspa, nspb

  implicit none

  private


  PUBLIC  :: chi_mls, lrtm_taumol

  real(wp), parameter :: oneminus = 1.0_wp - 1.0e-06_wp
  real(wp) :: chi_mls(7,59)

contains

  !----------------------------------------------------------------------------
  subroutine lrtm_taumol(kproma, nlayers, pavel, wx, coldry, &
       laytrop, jp, jt, jt1, &
       colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
       colbrd, fac00, fac01, fac10, fac11, &
       rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
       rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
       rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
       selffac, selffrac, indself, forfac, forfrac, indfor, &
       minorfrac, scaleminor, scaleminorn2, indminor, &
       fracs, taug)
    !----------------------------------------------------------------------------

    ! *******************************************************************************
    ! *                                                                             *
    ! *                  Optical depths developed for the                           *
    ! *                                                                             *
    ! *                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
    ! *                                                                             *
    ! *                                                                             *
    ! *            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
    ! *                        131 HARTWELL AVENUE                                  *
    ! *                        LEXINGTON, MA 02421                                  *
    ! *                                                                             *
    ! *                                                                             *
    ! *                           ELI J. MLAWER                                     *
    ! *                         JENNIFER DELAMERE                                   *
    ! *                         STEVEN J. TAUBMAN                                   *
    ! *                         SHEPARD A. CLOUGH                                   *
    ! *                                                                             *
    ! *                                                                             *
    ! *                                                                             *
    ! *                                                                             *
    ! *                       email:  mlawer@aer.com                                *
    ! *                       email:  jdelamer@aer.com                              *
    ! *                                                                             *
    ! *        The authors wish to acknowledge the contributions of the             *
    ! *        following people:  Karen Cady-Pereira, Patrick D. Brown,             *
    ! *        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
    ! *                                                                             *
    ! *******************************************************************************
    ! *                                                                             *
    ! *  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
    ! *                                                                             *
    ! *******************************************************************************
    ! *     TAUMOL                                                                  *
    ! *                                                                             *
    ! *     This file contains the subroutines TAUGBn (where n goes from            *
    ! *     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
    ! *     per g-value and layer for band n.                                       *
    ! *                                                                             *
    ! *  Output:  optical depths (unitless)                                         *
    ! *           fractions needed to compute Planck functions at every layer       *
    ! *               and g-value                                                   *
    ! *                                                                             *
    ! *     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
    ! *     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
    ! *                                                                             *
    ! *  Input                                                                      *
    ! *                                                                             *
    ! *     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
    ! *     COMMON /PRECISE/  ONEMINUS                                              *
    ! *     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
    ! *     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
    ! *     COMMON /PROFDATA/ LAYTROP,                                              *
    ! *    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
    ! *    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
    ! *    &                  COLO2(MXLAY)
    ! *     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
    ! *    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
    ! *     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
    ! *     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
    ! *                                                                             *
    ! *     Description:                                                            *
    ! *     NG(IBAND) - number of g-values in band IBAND                            *
    ! *     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
    ! *                   atmospheres that are stored for band IBAND per            *
    ! *                   pressure level and temperature.  Each of these            *
    ! *                   atmospheres has different relative amounts of the         *
    ! *                   key species for the band (i.e. different binary           *
    ! *                   species parameters).                                      *
    ! *     NSPB(IBAND) - same for upper atmosphere                                 *
    ! *     ONEMINUS - since problems are caused in some cases by interpolation     *
    ! *                parameters equal to or greater than 1, for these cases       *
    ! *                these parameters are set to this value, slightly < 1.        *
    ! *     PAVEL - layer pressures (mb)                                            *
    ! *     TAVEL - layer temperatures (degrees K)                                  *
    ! *     PZ - level pressures (mb)                                               *
    ! *     TZ - level temperatures (degrees K)                                     *
    ! *     LAYTROP - layer at which switch is made from one combination of         *
    ! *               key species to another                                        *
    ! *     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
    ! *               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
    ! *               respectively (molecules/cm**2)                                *
    ! *     FACij(LAY) - for layer LAY, these are factors that are needed to        *
    ! *                  compute the interpolation factors that multiply the        *
    ! *                  appropriate reference k-values.  A value of 0 (1) for      *
    ! *                  i,j indicates that the corresponding factor multiplies     *
    ! *                  reference k-value for the lower (higher) of the two        *
    ! *                  appropriate temperatures, and altitudes, respectively.     *
    ! *     JP - the index of the lower (in altitude) of the two appropriate        *
    ! *          reference pressure levels needed for interpolation                 *
    ! *     JT, JT1 - the indices of the lower of the two appropriate reference     *
    ! *               temperatures needed for interpolation (for pressure           *
    ! *               levels JP and JP+1, respectively)                             *
    ! *     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
    ! *               (water vapor density)/(atmospheric density at 296K and        *
    ! *               1013 mb)                                                      *
    ! *     SELFFRAC - factor needed for temperature interpolation of reference     *
    ! *                water vapor self-continuum data                              *
    ! *     INDSELF - index of the lower of the two appropriate reference           *
    ! *               temperatures needed for the self-continuum interpolation      *
    ! *     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
    ! *     FORFRAC - factor needed for temperature interpolation of reference      *
    ! *                water vapor foreign-continuum data                           *
    ! *     INDFOR  - index of the lower of the two appropriate reference           *
    ! *               temperatures needed for the foreign-continuum interpolation   *
    ! *                                                                             *
    ! *  Data input                                                                 *
    ! *     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
    ! *                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
    ! *        (note:  n is the band number,'MGAS' is the species name of the minor *
    ! *         gas)                                                                *
    ! *                                                                             *
    ! *     Description:                                                            *
    ! *     KA - k-values for low reference atmospheres (key-species only)          *
    ! *          (units: cm**2/molecule)                                            *
    ! *     KB - k-values for high reference atmospheres (key-species only)         *
    ! *          (units: cm**2/molecule)                                            *
    ! *     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
    ! *          (units: cm**2/molecule)                                            *
    ! *     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
    ! *          (units: cm**2/molecule)                                            *
    ! *     SELFREF - k-values for water vapor self-continuum for reference         *
    ! *               atmospheres (used below LAYTROP)                              *
    ! *               (units: cm**2/molecule)                                       *
    ! *     FORREF  - k-values for water vapor foreign-continuum for reference      *
    ! *               atmospheres (used below/above LAYTROP)                        *
    ! *               (units: cm**2/molecule)                                       *
    ! *                                                                             *
    ! *     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
    ! *     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
    ! *                                                                             *
    !*******************************************************************************

    ! ------- Declarations -------

    ! ----- Input -----
    integer, intent(in) :: kproma          ! number of columns
    integer, intent(in) :: nlayers         ! total number of layers
    real(wp), intent(in) :: pavel(:,:)     ! layer pressures (mb)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: wx(:,:,:)      ! cross-section amounts (mol/cm2)
    !    Dimensions: (maxxsec,nlayers)
    real(wp), intent(in) :: coldry(:,:)    ! column amount (dry air)
    !    Dimensions: (nlayers)

    integer, intent(in) :: laytrop(:)        ! tropopause layer index
    integer, intent(in) :: jp(:,:)           !
    !    Dimensions: (nlayers)
    integer, intent(in) :: jt(:,:)           !
    !    Dimensions: (nlayers)
    integer, intent(in) :: jt1(:,:)          !
    !    Dimensions: (nlayers)

    real(wp), intent(in) :: colh2o(:,:)          ! column amount (h2o)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: colco2(:,:)          ! column amount (co2)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: colo3(:,:)           ! column amount (o3)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: coln2o(:,:)          ! column amount (n2o)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: colco(:,:)           ! column amount (co)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: colch4(:,:)          ! column amount (ch4)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: colo2(:,:)           ! column amount (o2)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: colbrd(:,:)          ! column amount (broadening gases)
    !    Dimensions: (nlayers)

    integer, intent(in) :: indself(:,:)
    !    Dimensions: (nlayers)
    integer, intent(in) :: indfor(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: selffac(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: selffrac(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: forfac(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: forfrac(:,:)
    !    Dimensions: (nlayers)

    integer, intent(in) :: indminor(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: minorfrac(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: scaleminor(:,:)
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: scaleminorn2(:,:)
    !    Dimensions: (nlayers)

    real(wp), intent(in) :: &                  !
         fac00(:,:), fac01(:,:), &             !    Dimensions: (kproma,nlayers)
         fac10(:,:), fac11(:,:)
    real(wp), intent(in) :: &                  !
         rat_h2oco2(:,:),rat_h2oco2_1(:,:), &
         rat_h2oo3(:,:),rat_h2oo3_1(:,:), & !    Dimensions: (nlayers)
         rat_h2on2o(:,:),rat_h2on2o_1(:,:), &
         rat_h2och4(:,:),rat_h2och4_1(:,:), &
         rat_n2oco2(:,:),rat_n2oco2_1(:,:), &
         rat_o3co2(:,:),rat_o3co2_1(:,:)

    ! ----- Output -----
    real(wp), intent(out) :: fracs(:,:,:)        ! planck fractions
    !    Dimensions: (kproma,nlayers,ngptlw)
    real(wp), intent(out) :: taug(:,:,:)         ! gaseous optical depth
    !    Dimensions: (kproma,nlayers,ngptlw)

    integer :: icl, ich, jc, lay

    !     local integer arrays
    INTEGER :: laytrop_min, laytrop_max
    integer :: ixc(nlayers), ixlow(kproma,nlayers), ixhigh(kproma,nlayers)

    laytrop_min = MINVAL(laytrop)
    laytrop_max = MAXVAL(laytrop)

    ixlow  = 0
    ixhigh = 0
    ixc    = 0

    ! create index lists for mixed layers
    do lay = laytrop_min+1, laytrop_max
      icl = 0
      ich = 0
      do jc = 1, kproma
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


    ! Calculate gaseous optical depth and planck fractions for each spectral band.

    call taugb1
    call taugb2
    call taugb3
    call taugb4
    call taugb5
    call taugb6
    call taugb7
    call taugb8
    call taugb9
    call taugb10
    call taugb11
    call taugb12
    call taugb13
    call taugb14
    call taugb15
    call taugb16

  contains

    !----------------------------------------------------------------------------
    subroutine taugb1
      !----------------------------------------------------------------------------

      ! ------- Modifications -------
      !  Written by Eli J. Mlawer, Atmospheric & Environmental Research.
      !  Revised by Michael J. Iacono, Atmospheric & Environmental Research.
      !
      !     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
      !                          (high key - h2o; high minor - n2)
      !
      !     note: previous versions of rrtm band 1:
      !           10-250 cm-1 (low - h2o; high - h2o)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng1
      use mo_rrlw_kg01,   only : fracrefa, fracrefb, absa, absb, &
           ka_mn2, kb_mn2, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      real(wp) :: pp, corradj, scalen2, tauself, taufor, taun2


      ! Minor gas mapping levels:
      !     lower - n2, p = 142.5490 mbar, t = 215.70 k
      !     upper - n2, p = 142.5490 mbar, t = 215.70 k

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature.  Below laytrop, the water vapor self-continuum and
      ! foreign continuum is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(1) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(1) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1.0_wp
          if (pp .lt. 250._wp) then
            corradj = 1._wp - 0.15_wp * (250._wp-pp) / 154.4_wp
          endif

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$NEC unroll(NG1)
          do ig = 1, ng1
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) -  forref(indf,ig)))
            taun2 = scalen2*(ka_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
            taug(jl,lay,ig) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor + taun2)
            fracs(jl,lay,ig) = fracrefa(ig)
          enddo
        enddo
      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(1) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(1) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1._wp - 0.15_wp * (pp / 95.6_wp)

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$NEC unroll(NG1)
          do ig = 1, ng1
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taun2 = scalen2*(kb_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
            taug(jl,lay,ig) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor + taun2)
            fracs(jl,lay,ig) = fracrefb(ig)
          enddo
        enddo
      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(1) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(1) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1.0_wp
          if (pp .lt. 250._wp) then
            corradj = 1._wp - 0.15_wp * (250._wp-pp) / 154.4_wp
          endif

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$NEC unroll(NG1)
          do ig = 1, ng1
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) -  forref(indf,ig)))
            taun2 = scalen2*(ka_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (ka_mn2(indm+1,ig) - ka_mn2(indm,ig)))
            taug(jl,lay,ig) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor + taun2)
            fracs(jl,lay,ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(1) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(1) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          pp = pavel(jl,lay)
          corradj =  1._wp - 0.15_wp * (pp / 95.6_wp)

          scalen2 = colbrd(jl,lay) * scaleminorn2(jl,lay)
!$NEC unroll(NG1)
          do ig = 1, ng1
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taun2 = scalen2*(kb_mn2(indm,ig) + &
                 minorfrac(jl,lay) * (kb_mn2(indm+1,ig) - kb_mn2(indm,ig)))
            taug(jl,lay,ig) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor + taun2)
            fracs(jl,lay,ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb1

    !----------------------------------------------------------------------------
    subroutine taugb2
      !----------------------------------------------------------------------------
      !
      !     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
      !
      !     note: previous version of rrtm band 2:
      !           250 - 500 cm-1 (low - h2o; high - h2o)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng2, ngs1
      use mo_rrlw_kg02,   only : fracrefa, fracrefb, absa, absb, &
           selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf
      real(wp) :: pp, corradj, tauself, taufor


      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature.  Below laytrop, the water vapor self-continuum and
      ! foreign continuum is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(2) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(2) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          pp = pavel(jl,lay)
          corradj = 1._wp - .05_wp * (pp - 100._wp) / 900._wp
!$NEC unroll(NG2)
          do ig = 1, ng2
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs1+ig) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor)
            fracs(jl,lay,ngs1+ig) = fracrefa(ig)
          enddo
        enddo
      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(2) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(2) + 1
          indf = indfor(jl,lay)
!$NEC unroll(NG2)
          do ig = 1, ng2
            taufor =  forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs1+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(jl,lay,ngs1+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(2) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(2) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          pp = pavel(jl,lay)
          corradj = 1._wp - .05_wp * (pp - 100._wp) / 900._wp
!$NEC unroll(NG2)
          do ig = 1, ng2
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs1+ig) = corradj * (colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor)
            fracs(jl,lay,ngs1+ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(2) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(2) + 1
          indf = indfor(jl,lay)
!$NEC unroll(NG2)
          do ig = 1, ng2
            taufor =  forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs1+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(jl,lay,ngs1+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb2

    !----------------------------------------------------------------------------
    subroutine taugb3
      !----------------------------------------------------------------------------
      !
      !     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
      !                           (high key - h2o,co2; high minor - n2o)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng3, ngs2
      use mo_rrlw_kg03,   only : fracrefa, fracrefb, absa, absb, &
           ka_mn2o, kb_mn2o, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      integer :: js, js1, jmn2o, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, &
           fmn2o, chi_n2o, ratn2o, adjfac, adjcoln2o
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor, n2om1, n2om2, absn2o
      real(wp) :: refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b
      real(wp) :: tau_major(ng3), tau_major1(ng3)


      ! Minor gas mapping levels:
      !     lower - n2o, p = 706.272 mbar, t = 278.94 k
      !     upper - n2o, p = 95.58 mbar, t = 215.7 k

      !  P = 212.725 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)

      !  P = 95.58 mb
      refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)

      !  P = 706.270mb
      refrat_m_a = chi_mls(1,3)/chi_mls(2,3)

      !  P = 95.58 mb
      refrat_m_b = chi_mls(1,13)/chi_mls(2,13)

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature, and appropriate species.  Below laytrop, the water vapor
      ! self-continuum and foreign continuum is interpolated (in temperature)
      ! separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._wp*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - AINT((specmult_mn2o)))
          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(3) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(3) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *    &
             (fac000 * absa(ind0,1:ng3)    + &
              fac100 * absa(ind0+1,1:ng3)  + &
              fac200 * absa(ind0+2,1:ng3)  + &
              fac010 * absa(ind0+9,1:ng3)  + &
              fac110 * absa(ind0+10,1:ng3) + &
              fac210 * absa(ind0+11,1:ng3))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng3) + &
              fac100 * absa(ind0,1:ng3)   + &
              fac000 * absa(ind0+1,1:ng3) + &
              fac210 * absa(ind0+8,1:ng3) + &
              fac110 * absa(ind0+9,1:ng3) + &
              fac010 * absa(ind0+10,1:ng3))
          else
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *   &
             (fac000 * absa(ind0,1:ng3)   + &
              fac100 * absa(ind0+1,1:ng3) + &
              fac010 * absa(ind0+9,1:ng3) + &
              fac110 * absa(ind0+10,1:ng3))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng3)    + &
              fac101 * absa(ind1+1,1:ng3)  + &
              fac201 * absa(ind1+2,1:ng3)  + &
              fac011 * absa(ind1+9,1:ng3)  + &
              fac111 * absa(ind1+10,1:ng3) + &
              fac211 * absa(ind1+11,1:ng3))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng3) + &
              fac101 * absa(ind1,1:ng3)   + &
              fac001 * absa(ind1+1,1:ng3) + &
              fac211 * absa(ind1+8,1:ng3) + &
              fac111 * absa(ind1+9,1:ng3) + &
              fac011 * absa(ind1+10,1:ng3))
          else
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 * &
             (fac001 * absa(ind1,1:ng3)   + &
              fac101 * absa(ind1+1,1:ng3) + &
              fac011 * absa(ind1+9,1:ng3) + &
              fac111 * absa(ind1+10,1:ng3))
          endif

!$NEC unroll(NG3)
          do ig = 1, ng3
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

            taug(jl,lay,ngs2+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs2+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 4._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          fac000 = (1._wp - fs) * fac00(jl,lay)
          fac010 = (1._wp - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._wp - fs1) * fac01(jl,lay)
          fac011 = (1._wp - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_b*colco2(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 4._wp*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - AINT((specmult_mn2o)))
          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(3) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(3) + js1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG3)
          do ig = 1, ng3
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
            n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)
            taug(jl,lay,ngs2+ig) = speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) +  &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))  &
                 + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs2+ig) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max

        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._wp*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - AINT((specmult_mn2o)))
          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(3) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(3) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *    &
             (fac000 * absa(ind0,1:ng3)    + &
              fac100 * absa(ind0+1,1:ng3)  + &
              fac200 * absa(ind0+2,1:ng3)  + &
              fac010 * absa(ind0+9,1:ng3)  + &
              fac110 * absa(ind0+10,1:ng3) + &
              fac210 * absa(ind0+11,1:ng3))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng3) + &
              fac100 * absa(ind0,1:ng3)   + &
              fac000 * absa(ind0+1,1:ng3) + &
              fac210 * absa(ind0+8,1:ng3) + &
              fac110 * absa(ind0+9,1:ng3) + &
              fac010 * absa(ind0+10,1:ng3))
          else
!$NEC unroll(NG3)
            tau_major(1:ng3) = speccomb *   &
             (fac000 * absa(ind0,1:ng3)   + &
              fac100 * absa(ind0+1,1:ng3) + &
              fac010 * absa(ind0+9,1:ng3) + &
              fac110 * absa(ind0+10,1:ng3))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng3)    + &
              fac101 * absa(ind1+1,1:ng3)  + &
              fac201 * absa(ind1+2,1:ng3)  + &
              fac011 * absa(ind1+9,1:ng3)  + &
              fac111 * absa(ind1+10,1:ng3) + &
              fac211 * absa(ind1+11,1:ng3))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng3) + &
              fac101 * absa(ind1,1:ng3)   + &
              fac001 * absa(ind1+1,1:ng3) + &
              fac211 * absa(ind1+8,1:ng3) + &
              fac111 * absa(ind1+9,1:ng3) + &
              fac011 * absa(ind1+10,1:ng3))
          else
!$NEC unroll(NG3)
            tau_major1(1:ng3) = speccomb1 * &
             (fac001 * absa(ind1,1:ng3)   + &
              fac101 * absa(ind1+1,1:ng3) + &
              fac011 * absa(ind1+9,1:ng3) + &
              fac111 * absa(ind1+10,1:ng3))
          endif

!$NEC unroll(NG3)
          do ig = 1, ng3
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = ka_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm,ig) - ka_mn2o(jmn2o,indm,ig))
            n2om2 = ka_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (ka_mn2o(jmn2o+1,indm+1,ig) - ka_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)

            taug(jl,lay,ngs2+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs2+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 4._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          fac000 = (1._wp - fs) * fac00(jl,lay)
          fac010 = (1._wp - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._wp - fs1) * fac01(jl,lay)
          fac011 = (1._wp - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_b*colco2(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 4._wp*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - AINT((specmult_mn2o)))
          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/coldry(jl,lay)
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(3) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(3) + js1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG3)
          do ig = 1, ng3
            taufor = forfac(jl,lay) * (forref(indf,ig) + &
                 forfrac(jl,lay) * (forref(indf+1,ig) - forref(indf,ig)))
            n2om1 = kb_mn2o(jmn2o,indm,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm,ig)-kb_mn2o(jmn2o,indm,ig))
            n2om2 = kb_mn2o(jmn2o,indm+1,ig) + fmn2o * &
                 (kb_mn2o(jmn2o+1,indm+1,ig)-kb_mn2o(jmn2o,indm+1,ig))
            absn2o = n2om1 + minorfrac(jl,lay) * (n2om2 - n2om1)
            taug(jl,lay,ngs2+ig) = speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) +  &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))  &
                 + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs2+ig) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo

      enddo

    end subroutine taugb3

    !----------------------------------------------------------------------------
    subroutine taugb4
      !----------------------------------------------------------------------------
      !
      !     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng4, ngs3
      use mo_rrlw_kg04,   only : fracrefa, fracrefb, absa, absb, &
           selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf
      integer :: js, js1, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor
      real(wp) :: refrat_planck_a, refrat_planck_b
      real(wp) :: tau_major(ng4), tau_major1(ng4)


      ! P =   142.5940 mb
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)

      ! P = 95.58350 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated (in temperature)
      ! separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(4) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(4) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
         else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG4)
            tau_major(1:ng4) = speccomb *    &
             (fac000 * absa(ind0,1:ng4)    + &
              fac100 * absa(ind0+1,1:ng4)  + &
              fac200 * absa(ind0+2,1:ng4)  + &
              fac010 * absa(ind0+9,1:ng4)  + &
              fac110 * absa(ind0+10,1:ng4) + &
              fac210 * absa(ind0+11,1:ng4))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG4)
            tau_major(1:ng4) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng4) + &
              fac100 * absa(ind0,1:ng4)   + &
              fac000 * absa(ind0+1,1:ng4) + &
              fac210 * absa(ind0+8,1:ng4) + &
              fac110 * absa(ind0+9,1:ng4) + &
              fac010 * absa(ind0+10,1:ng4))
          else
!$NEC unroll(NG4)
             tau_major(1:ng4) = speccomb *   &
              (fac000 * absa(ind0,1:ng4)   + &
               fac100 * absa(ind0+1,1:ng4) + &
               fac010 * absa(ind0+9,1:ng4) + &
               fac110 * absa(ind0+10,1:ng4))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng4)    + &
              fac101 * absa(ind1+1,1:ng4)  + &
              fac201 * absa(ind1+2,1:ng4)  + &
              fac011 * absa(ind1+9,1:ng4)  + &
              fac111 * absa(ind1+10,1:ng4) + &
              fac211 * absa(ind1+11,1:ng4))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng4) + &
              fac101 * absa(ind1,1:ng4)   + &
              fac001 * absa(ind1+1,1:ng4) + &
              fac211 * absa(ind1+8,1:ng4) + &
              fac111 * absa(ind1+9,1:ng4) + &
              fac011 * absa(ind1+10,1:ng4))
          else
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 * &
             (fac001 * absa(ind1,1:ng4)   + &
              fac101 * absa(ind1+1,1:ng4) + &
              fac011 * absa(ind1+9,1:ng4) + &
              fac111 * absa(ind1+10,1:ng4))
          endif

!$NEC unroll(NG4)
          do ig = 1, ng4
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,lay,ngs3+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,lay,ngs3+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          speccomb = colo3(jl,lay) + rat_o3co2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colo3(jl,lay)/speccomb,oneminus)
          specmult = 4._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colo3(jl,lay) + rat_o3co2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colo3(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          fac000 = (1._wp - fs) * fac00(jl,lay)
          fac010 = (1._wp - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._wp - fs1) * fac01(jl,lay)
          fac011 = (1._wp - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_planck = colo3(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colo3(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(4) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(4) + js1
!$NEC unroll(NG4)
          do ig = 1, ng4
            taug(jl,lay,ngs3+ig) =  speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) +  &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))
            fracs(jl,lay,ngs3+ig) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo
      enddo

      ! Empirical modification to code to improve stratospheric cooling rates
      ! for co2.  Revised to apply weighting for g-point reduction in this band.
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma
          taug(jl,lay,ngs3+8)=taug(jl,lay,ngs3+8)*0.92_wp
          taug(jl,lay,ngs3+9)=taug(jl,lay,ngs3+9)*0.88_wp
          taug(jl,lay,ngs3+10)=taug(jl,lay,ngs3+10)*1.07_wp
          taug(jl,lay,ngs3+11)=taug(jl,lay,ngs3+11)*1.1_wp
          taug(jl,lay,ngs3+12)=taug(jl,lay,ngs3+12)*0.99_wp
          taug(jl,lay,ngs3+13)=taug(jl,lay,ngs3+13)*0.88_wp
          taug(jl,lay,ngs3+14)=taug(jl,lay,ngs3+14)*0.943_wp
        enddo
      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      DO lay = laytrop_min+1, laytrop_max

        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(4) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(4) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
         else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG4)
            tau_major(1:ng4) = speccomb *    &
             (fac000 * absa(ind0,1:ng4)    + &
              fac100 * absa(ind0+1,1:ng4)  + &
              fac200 * absa(ind0+2,1:ng4)  + &
              fac010 * absa(ind0+9,1:ng4)  + &
              fac110 * absa(ind0+10,1:ng4) + &
              fac210 * absa(ind0+11,1:ng4))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG4)
            tau_major(1:ng4) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng4) + &
              fac100 * absa(ind0,1:ng4)   + &
              fac000 * absa(ind0+1,1:ng4) + &
              fac210 * absa(ind0+8,1:ng4) + &
              fac110 * absa(ind0+9,1:ng4) + &
              fac010 * absa(ind0+10,1:ng4))
          else
!$NEC unroll(NG4)
             tau_major(1:ng4) = speccomb *   &
              (fac000 * absa(ind0,1:ng4)   + &
               fac100 * absa(ind0+1,1:ng4) + &
               fac010 * absa(ind0+9,1:ng4) + &
               fac110 * absa(ind0+10,1:ng4))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng4)    + &
              fac101 * absa(ind1+1,1:ng4)  + &
              fac201 * absa(ind1+2,1:ng4)  + &
              fac011 * absa(ind1+9,1:ng4)  + &
              fac111 * absa(ind1+10,1:ng4) + &
              fac211 * absa(ind1+11,1:ng4))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng4) + &
              fac101 * absa(ind1,1:ng4)   + &
              fac001 * absa(ind1+1,1:ng4) + &
              fac211 * absa(ind1+8,1:ng4) + &
              fac111 * absa(ind1+9,1:ng4) + &
              fac011 * absa(ind1+10,1:ng4))
          else
!$NEC unroll(NG4)
            tau_major1(1:ng4) = speccomb1 * &
             (fac001 * absa(ind1,1:ng4)   + &
              fac101 * absa(ind1+1,1:ng4) + &
              fac011 * absa(ind1+9,1:ng4) + &
              fac111 * absa(ind1+10,1:ng4))
          endif

!$NEC unroll(NG4)
          do ig = 1, ng4
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,lay,ngs3+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,lay,ngs3+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          speccomb = colo3(jl,lay) + rat_o3co2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colo3(jl,lay)/speccomb,oneminus)
          specmult = 4._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colo3(jl,lay) + rat_o3co2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colo3(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          fac000 = (1._wp - fs) * fac00(jl,lay)
          fac010 = (1._wp - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._wp - fs1) * fac01(jl,lay)
          fac011 = (1._wp - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_planck = colo3(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colo3(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(4) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(4) + js1
!$NEC unroll(NG4)
          do ig = 1, ng4
            taug(jl,lay,ngs3+ig) =  speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) +  &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))
            fracs(jl,lay,ngs3+ig) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo

        ! Empirical modification to code to improve stratospheric cooling rates
        ! for co2.  Revised to apply weighting for g-point reduction in this band.
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)
          taug(jl,lay,ngs3+8)=taug(jl,lay,ngs3+8)*0.92_wp
          taug(jl,lay,ngs3+9)=taug(jl,lay,ngs3+9)*0.88_wp
          taug(jl,lay,ngs3+10)=taug(jl,lay,ngs3+10)*1.07_wp
          taug(jl,lay,ngs3+11)=taug(jl,lay,ngs3+11)*1.1_wp
          taug(jl,lay,ngs3+12)=taug(jl,lay,ngs3+12)*0.99_wp
          taug(jl,lay,ngs3+13)=taug(jl,lay,ngs3+13)*0.88_wp
          taug(jl,lay,ngs3+14)=taug(jl,lay,ngs3+14)*0.943_wp
        enddo

      ENDDO

    end subroutine taugb4



    !----------------------------------------------------------------------------
    subroutine taugb5
      !----------------------------------------------------------------------------
      !
      !     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
      !                           (high key - o3,co2)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng5, ngs4
      use mo_rrlw_kg05,   only : fracrefa, fracrefb, absa, absb, &
           ka_mo3, selfref, forref, ccl4

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      integer :: js, js1, jmo3, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_mo3, specparm_mo3, specmult_mo3, fmo3
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor, o3m1, o3m2, abso3
      real(wp) :: refrat_planck_a, refrat_planck_b, refrat_m_a
      real(wp) :: tau_major(ng5), tau_major1(ng5)


      ! Minor gas mapping level :
      !     lower - o3, p = 317.34 mbar, t = 240.77 k
      !     lower - ccl4

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P = 473.420 mb
      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)

      ! P = 0.2369 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)

      ! P = 317.3480
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature, and appropriate species.  Below laytrop, the
      ! water vapor self-continuum and foreign continuum is
      ! interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mo3 = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mo3 = MIN(colh2o(jl,lay)/speccomb_mo3,oneminus)
          specmult_mo3 = 8._wp*specparm_mo3
          jmo3 = 1 + int(specmult_mo3)
          fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(5) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(5) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *    &
             (fac000 * absa(ind0,1:ng5)    + &
              fac100 * absa(ind0+1,1:ng5)  + &
              fac200 * absa(ind0+2,1:ng5)  + &
              fac010 * absa(ind0+9,1:ng5)  + &
              fac110 * absa(ind0+10,1:ng5) + &
              fac210 * absa(ind0+11,1:ng5))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng5) + &
              fac100 * absa(ind0,1:ng5)   + &
              fac000 * absa(ind0+1,1:ng5) + &
              fac210 * absa(ind0+8,1:ng5) + &
              fac110 * absa(ind0+9,1:ng5) + &
              fac010 * absa(ind0+10,1:ng5))
          else
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *   &
             (fac000 * absa(ind0,1:ng5)   + &
              fac100 * absa(ind0+1,1:ng5) + &
              fac010 * absa(ind0+9,1:ng5) + &
              fac110 * absa(ind0+10,1:ng5))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng5)    + &
              fac101 * absa(ind1+1,1:ng5)  + &
              fac201 * absa(ind1+2,1:ng5)  + &
              fac011 * absa(ind1+9,1:ng5)  + &
              fac111 * absa(ind1+10,1:ng5) + &
              fac211 * absa(ind1+11,1:ng5))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng5) + &
              fac101 * absa(ind1,1:ng5)   + &
              fac001 * absa(ind1+1,1:ng5) + &
              fac211 * absa(ind1+8,1:ng5) + &
              fac111 * absa(ind1+9,1:ng5) + &
              fac011 * absa(ind1+10,1:ng5))
          else
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 * &
             (fac001 * absa(ind1,1:ng5)   + &
              fac101 * absa(ind1+1,1:ng5) + &
              fac011 * absa(ind1+9,1:ng5) + &
              fac111 * absa(ind1+10,1:ng5))
          endif

!$NEC unroll(NG5)
          do ig = 1, ng5
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
            o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(jl,lay)*(o3m2-o3m1)

            taug(jl,lay,ngs4+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + abso3*colo3(jl,lay) &
                 + wx(jl,1,lay) * ccl4(ig)
            fracs(jl,lay,ngs4+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          speccomb = colo3(jl,lay) + rat_o3co2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colo3(jl,lay)/speccomb,oneminus)
          specmult = 4._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colo3(jl,lay) + rat_o3co2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colo3(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          fac000 = (1._wp - fs) * fac00(jl,lay)
          fac010 = (1._wp - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._wp - fs1) * fac01(jl,lay)
          fac011 = (1._wp - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_planck = colo3(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colo3(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(5) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(5) + js1
!$NEC unroll(NG5)
          do ig = 1, ng5
            taug(jl,lay,ngs4+ig) = speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) + &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))  &
                 + wx(jl,1,lay) * ccl4(ig)
            fracs(jl,lay,ngs4+ig) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max

        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mo3 = colh2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mo3 = MIN(colh2o(jl,lay)/speccomb_mo3,oneminus)
          specmult_mo3 = 8._wp*specparm_mo3
          jmo3 = 1 + int(specmult_mo3)
          fmo3 = ((specmult_mo3) - AINT((specmult_mo3)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(5) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(5) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *    &
             (fac000 * absa(ind0,1:ng5)    + &
              fac100 * absa(ind0+1,1:ng5)  + &
              fac200 * absa(ind0+2,1:ng5)  + &
              fac010 * absa(ind0+9,1:ng5)  + &
              fac110 * absa(ind0+10,1:ng5) + &
              fac210 * absa(ind0+11,1:ng5))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng5) + &
              fac100 * absa(ind0,1:ng5)   + &
              fac000 * absa(ind0+1,1:ng5) + &
              fac210 * absa(ind0+8,1:ng5) + &
              fac110 * absa(ind0+9,1:ng5) + &
              fac010 * absa(ind0+10,1:ng5))
          else
!$NEC unroll(NG5)
            tau_major(1:ng5) = speccomb *   &
             (fac000 * absa(ind0,1:ng5)   + &
              fac100 * absa(ind0+1,1:ng5) + &
              fac010 * absa(ind0+9,1:ng5) + &
              fac110 * absa(ind0+10,1:ng5))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng5)    + &
              fac101 * absa(ind1+1,1:ng5)  + &
              fac201 * absa(ind1+2,1:ng5)  + &
              fac011 * absa(ind1+9,1:ng5)  + &
              fac111 * absa(ind1+10,1:ng5) + &
              fac211 * absa(ind1+11,1:ng5))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng5) + &
              fac101 * absa(ind1,1:ng5)   + &
              fac001 * absa(ind1+1,1:ng5) + &
              fac211 * absa(ind1+8,1:ng5) + &
              fac111 * absa(ind1+9,1:ng5) + &
              fac011 * absa(ind1+10,1:ng5))
          else
!$NEC unroll(NG5)
            tau_major1(1:ng5) = speccomb1 * &
             (fac001 * absa(ind1,1:ng5)   + &
              fac101 * absa(ind1+1,1:ng5) + &
              fac011 * absa(ind1+9,1:ng5) + &
              fac111 * absa(ind1+10,1:ng5))
          endif

!$NEC unroll(NG5)
          do ig = 1, ng5
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            o3m1 = ka_mo3(jmo3,indm,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm,ig)-ka_mo3(jmo3,indm,ig))
            o3m2 = ka_mo3(jmo3,indm+1,ig) + fmo3 * &
                 (ka_mo3(jmo3+1,indm+1,ig)-ka_mo3(jmo3,indm+1,ig))
            abso3 = o3m1 + minorfrac(jl,lay)*(o3m2-o3m1)

            taug(jl,lay,ngs4+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + abso3*colo3(jl,lay) &
                 + wx(jl,1,lay) * ccl4(ig)
            fracs(jl,lay,ngs4+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          speccomb = colo3(jl,lay) + rat_o3co2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colo3(jl,lay)/speccomb,oneminus)
          specmult = 4._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colo3(jl,lay) + rat_o3co2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colo3(jl,lay)/speccomb1,oneminus)
          specmult1 = 4._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          fac000 = (1._wp - fs) * fac00(jl,lay)
          fac010 = (1._wp - fs) * fac10(jl,lay)
          fac100 = fs * fac00(jl,lay)
          fac110 = fs * fac10(jl,lay)
          fac001 = (1._wp - fs1) * fac01(jl,lay)
          fac011 = (1._wp - fs1) * fac11(jl,lay)
          fac101 = fs1 * fac01(jl,lay)
          fac111 = fs1 * fac11(jl,lay)

          speccomb_planck = colo3(jl,lay)+refrat_planck_b*colco2(jl,lay)
          specparm_planck = MIN(colo3(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 4._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(5) + js
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(5) + js1
!$NEC unroll(NG5)
          do ig = 1, ng5
            taug(jl,lay,ngs4+ig) = speccomb * &
                 (fac000 * absb(ind0,ig) + &
                 fac100 * absb(ind0+1,ig) + &
                 fac010 * absb(ind0+5,ig) + &
                 fac110 * absb(ind0+6,ig)) &
                 + speccomb1 * &
                 (fac001 * absb(ind1,ig) + &
                 fac101 * absb(ind1+1,ig) + &
                 fac011 * absb(ind1+5,ig) + &
                 fac111 * absb(ind1+6,ig))  &
                 + wx(jl,1,lay) * ccl4(ig)
            fracs(jl,lay,ngs4+ig) = fracrefb(ig,jpl) + fpl * &
                 (fracrefb(ig,jpl+1)-fracrefb(ig,jpl))
          enddo
        enddo

      enddo

    end subroutine taugb5

    !----------------------------------------------------------------------------
    subroutine taugb6
      !----------------------------------------------------------------------------
      !
      !     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
      !                           (high key - nothing; high minor - cfc11, cfc12)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng6, ngs5
      use mo_rrlw_kg06,   only : fracrefa, absa, ka_mco2, &
           selfref, forref, cfc11adj, cfc12

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      real(wp) :: chi_co2, ratco2, adjfac, adjcolco2
      real(wp) :: tauself, taufor, absco2


      ! Minor gas mapping level:
      !     lower - co2, p = 706.2720 mb, t = 294.2 k
      !     upper - cfc11, cfc12

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature. The water vapor self-continuum and foreign continuum
      ! is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          ! In atmospheres where the amount of CO2 is too great to be considered
          ! a minor species, adjust the column amount of CO2 by an empirical factor
          ! to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.77_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(6) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(6) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG6)
          do ig = 1, ng6
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(jl,lay,ngs5+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,lay,ngs5+ig) = fracrefa(ig)
          enddo
        enddo
      enddo

      ! Upper atmosphere loop
      ! Nothing important goes on above laytrop in this band.
      do ig = 1, ng6
        do lay = laytrop_max+1, nlayers
          do jl = 1, kproma
            taug(jl,lay,ngs5+ig) = 0.0_wp &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,lay,ngs5+ig) = fracrefa(ig)
          enddo
        enddo
      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          ! In atmospheres where the amount of CO2 is too great to be considered
          ! a minor species, adjust the column amount of CO2 by an empirical factor
          ! to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.77_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(6) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(6) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG6)
          do ig = 1, ng6
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
               (forref(indf+1,ig) - forref(indf,ig)))
            absco2 =  (ka_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mco2(indm+1,ig) - ka_mco2(indm,ig)))
            taug(jl,lay,ngs5+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) +  &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor &
                 + adjcolco2 * absco2 &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,lay,ngs5+ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ! Nothing important goes on above laytrop in this band.
        ixc0 = kproma - ixc0

        do ig = 1, ng6
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)
            taug(jl,lay,ngs5+ig) = 0.0_wp &
                 + wx(jl,2,lay) * cfc11adj(ig) &
                 + wx(jl,3,lay) * cfc12(ig)
            fracs(jl,lay,ngs5+ig) = fracrefa(ig)
          enddo
        enddo

      enddo

    end subroutine taugb6

    !----------------------------------------------------------------------------
    subroutine taugb7
      !----------------------------------------------------------------------------
      !
      !     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
      !                            (high key - o3; high minor - co2)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng7, ngs6
      use mo_rrlw_kg07,   only : fracrefa, fracrefb, absa, absb, &
           ka_mco2, kb_mco2, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      integer :: js, js1, jmco2, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor, co2m1, co2m2, absco2
      real(wp) :: chi_co2, ratco2, adjfac, adjcolco2
      real(wp) :: refrat_planck_a, refrat_m_a
      real(wp) :: tau_major(ng7), tau_major1(ng7)


      ! Minor gas mapping level :
      !     lower - co2, p = 706.2620 mbar, t= 278.94 k
      !     upper - co2, p = 12.9350 mbar, t = 234.01 k

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 706.2620 mb
      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)

      ! P = 706.2720 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2oo3(jl,lay)*colo3(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oo3_1(jl,lay)*colo3(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*colo3(jl,lay)
          specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._wp*specparm_mco2

          jmco2 = 1 + int(specmult_mco2)
          fmco2 = ((specmult_mco2) - AINT((specmult_mco2)))

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 3.0_wp+(ratco2-3.0_wp)**0.79_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colo3(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(7) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(7) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *    &
             (fac000 * absa(ind0,1:ng7)    + &
              fac100 * absa(ind0+1,1:ng7)  + &
              fac200 * absa(ind0+2,1:ng7)  + &
              fac010 * absa(ind0+9,1:ng7)  + &
              fac110 * absa(ind0+10,1:ng7) + &
              fac210 * absa(ind0+11,1:ng7))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng7) + &
              fac100 * absa(ind0,1:ng7)   + &
              fac000 * absa(ind0+1,1:ng7) + &
              fac210 * absa(ind0+8,1:ng7) + &
              fac110 * absa(ind0+9,1:ng7) + &
              fac010 * absa(ind0+10,1:ng7))
          else
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *   &
             (fac000 * absa(ind0,1:ng7)   + &
              fac100 * absa(ind0+1,1:ng7) + &
              fac010 * absa(ind0+9,1:ng7) + &
              fac110 * absa(ind0+10,1:ng7))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng7)    + &
              fac101 * absa(ind1+1,1:ng7)  + &
              fac201 * absa(ind1+2,1:ng7)  + &
              fac011 * absa(ind1+9,1:ng7)  + &
              fac111 * absa(ind1+10,1:ng7) + &
              fac211 * absa(ind1+11,1:ng7))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng7) + &
              fac101 * absa(ind1,1:ng7)   + &
              fac001 * absa(ind1+1,1:ng7) + &
              fac211 * absa(ind1+8,1:ng7) + &
              fac111 * absa(ind1+9,1:ng7) + &
              fac011 * absa(ind1+10,1:ng7))
          else
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 * &
             (fac001 * absa(ind1,1:ng7)   + &
              fac101 * absa(ind1+1,1:ng7) + &
              fac011 * absa(ind1+9,1:ng7) + &
              fac111 * absa(ind1+10,1:ng7))
          endif

!$NEC unroll(NG7)
          do ig = 1, ng7
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)

            taug(jl,lay,ngs6+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2
            fracs(jl,lay,ngs6+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.79_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(7) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(7) + 1
          indm = indminor(jl,lay)
!$NEC unroll(NG7)
          do ig = 1, ng7
            absco2 = kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
            taug(jl,lay,ngs6+ig) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2 * absco2
            fracs(jl,lay,ngs6+ig) = fracrefb(ig)
          enddo
        enddo
      enddo

      ! Empirical modification to code to improve stratospheric cooling rates
      ! for o3.  Revised to apply weighting for g-point reduction in this band.
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma
          taug(jl,lay,ngs6+6)=taug(jl,lay,ngs6+6)*0.92_wp
          taug(jl,lay,ngs6+7)=taug(jl,lay,ngs6+7)*0.88_wp
          taug(jl,lay,ngs6+8)=taug(jl,lay,ngs6+8)*1.07_wp
          taug(jl,lay,ngs6+9)=taug(jl,lay,ngs6+9)*1.1_wp
          taug(jl,lay,ngs6+10)=taug(jl,lay,ngs6+10)*0.99_wp
          taug(jl,lay,ngs6+11)=taug(jl,lay,ngs6+11)*0.855_wp
        enddo
      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2oo3(jl,lay)*colo3(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oo3_1(jl,lay)*colo3(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*colo3(jl,lay)
          specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._wp*specparm_mco2

          jmco2 = 1 + int(specmult_mco2)
          fmco2 = ((specmult_mco2) - AINT((specmult_mco2)))

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 3.0_wp+(ratco2-3.0_wp)**0.79_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colo3(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(7) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(7) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *    &
             (fac000 * absa(ind0,1:ng7)    + &
              fac100 * absa(ind0+1,1:ng7)  + &
              fac200 * absa(ind0+2,1:ng7)  + &
              fac010 * absa(ind0+9,1:ng7)  + &
              fac110 * absa(ind0+10,1:ng7) + &
              fac210 * absa(ind0+11,1:ng7))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng7) + &
              fac100 * absa(ind0,1:ng7)   + &
              fac000 * absa(ind0+1,1:ng7) + &
              fac210 * absa(ind0+8,1:ng7) + &
              fac110 * absa(ind0+9,1:ng7) + &
              fac010 * absa(ind0+10,1:ng7))
          else
!$NEC unroll(NG7)
            tau_major(1:ng7) = speccomb *   &
             (fac000 * absa(ind0,1:ng7)   + &
              fac100 * absa(ind0+1,1:ng7) + &
              fac010 * absa(ind0+9,1:ng7) + &
              fac110 * absa(ind0+10,1:ng7))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng7)    + &
              fac101 * absa(ind1+1,1:ng7)  + &
              fac201 * absa(ind1+2,1:ng7)  + &
              fac011 * absa(ind1+9,1:ng7)  + &
              fac111 * absa(ind1+10,1:ng7) + &
              fac211 * absa(ind1+11,1:ng7))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng7) + &
              fac101 * absa(ind1,1:ng7)   + &
              fac001 * absa(ind1+1,1:ng7) + &
              fac211 * absa(ind1+8,1:ng7) + &
              fac111 * absa(ind1+9,1:ng7) + &
              fac011 * absa(ind1+10,1:ng7))
          else
!$NEC unroll(NG7)
            tau_major1(1:ng7) = speccomb1 * &
             (fac001 * absa(ind1,1:ng7)   + &
              fac101 * absa(ind1+1,1:ng7) + &
              fac011 * absa(ind1+9,1:ng7) + &
              fac111 * absa(ind1+10,1:ng7))
          endif

!$NEC unroll(NG7)
          do ig = 1, ng7
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            co2m1 = ka_mco2(jmco2,indm,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm,ig) - ka_mco2(jmco2,indm,ig))
            co2m2 = ka_mco2(jmco2,indm+1,ig) + fmco2 * &
                 (ka_mco2(jmco2+1,indm+1,ig) - ka_mco2(jmco2,indm+1,ig))
            absco2 = co2m1 + minorfrac(jl,lay) * (co2m2 - co2m1)

            taug(jl,lay,ngs6+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2
            fracs(jl,lay,ngs6+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.79_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(7) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(7) + 1
          indm = indminor(jl,lay)
!$NEC unroll(NG7)
          do ig = 1, ng7
            absco2 = kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig))
            taug(jl,lay,ngs6+ig) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2 * absco2
            fracs(jl,lay,ngs6+ig) = fracrefb(ig)
          enddo
        enddo

        ! Empirical modification to code to improve stratospheric cooling rates
        ! for o3.  Revised to apply weighting for g-point reduction in this band.

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)
          taug(jl,lay,ngs6+6)=taug(jl,lay,ngs6+6)*0.92_wp
          taug(jl,lay,ngs6+7)=taug(jl,lay,ngs6+7)*0.88_wp
          taug(jl,lay,ngs6+8)=taug(jl,lay,ngs6+8)*1.07_wp
          taug(jl,lay,ngs6+9)=taug(jl,lay,ngs6+9)*1.1_wp
          taug(jl,lay,ngs6+10)=taug(jl,lay,ngs6+10)*0.99_wp
          taug(jl,lay,ngs6+11)=taug(jl,lay,ngs6+11)*0.855_wp
        enddo

      enddo

    end subroutine taugb7

    !----------------------------------------------------------------------------
    subroutine taugb8
      !----------------------------------------------------------------------------
      !
      !     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
      !                             (high key - o3; high minor - co2, n2o)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng8, ngs7
      use mo_rrlw_kg08,   only : fracrefa, fracrefb, absa, absb, &
           ka_mco2, ka_mn2o, ka_mo3, kb_mco2, kb_mn2o, &
           selfref, forref, cfc12, cfc22adj

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      real(wp) :: tauself, taufor, absco2, abso3, absn2o
      real(wp) :: chi_co2, ratco2, adjfac, adjcolco2


      ! Minor gas mapping level:
      !     lower - co2, p = 1053.63 mb, t = 294.2 k
      !     lower - o3,  p = 317.348 mb, t = 240.77 k
      !     lower - n2o, p = 706.2720 mb, t= 278.94 k
      !     lower - cfc12,cfc11
      !     upper - co2, p = 35.1632 mb, t = 223.28 k
      !     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature, and appropriate species.  Below laytrop, the water vapor
      ! self-continuum and foreign continuum is interpolated (in temperature)
      ! separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(8) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(8) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG8)
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
            taug(jl,lay,ngs7+ig) = colh2o(jl,lay) * &
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
            fracs(jl,lay,ngs7+ig) = fracrefa(ig)
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/coldry(jl,lay)
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1) * coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(8) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(8) + 1
          indm = indminor(jl,lay)
!$NEC unroll(NG8)
          do ig = 1, ng8
            absco2 =  (kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
            absn2o =  (kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
            taug(jl,lay,ngs7+ig) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(jl,lay)*absn2o &
                 + wx(jl,3,lay) * cfc12(ig) &
                 + wx(jl,4,lay) * cfc22adj(ig)
            fracs(jl,lay,ngs7+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(8) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(8) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
!$NEC unroll(NG8)
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
            taug(jl,lay,ngs7+ig) = colh2o(jl,lay) * &
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
            fracs(jl,lay,ngs7+ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere loop
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/coldry(jl,lay)
          ratco2 = 1.e20_wp*chi_co2/chi_mls(2,jp(jl,lay)+1)
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.65_wp
            adjcolco2 = adjfac*chi_mls(2,jp(jl,lay)+1) * coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(8) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(8) + 1
          indm = indminor(jl,lay)
!$NEC unroll(NG8)
          do ig = 1, ng8
            absco2 =  (kb_mco2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mco2(indm+1,ig) - kb_mco2(indm,ig)))
            absn2o =  (kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig)))
            taug(jl,lay,ngs7+ig) = colo3(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcolco2*absco2 &
                 + coln2o(jl,lay)*absn2o &
                 + wx(jl,3,lay) * cfc12(ig) &
                 + wx(jl,4,lay) * cfc22adj(ig)
            fracs(jl,lay,ngs7+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb8

    !----------------------------------------------------------------------------
    subroutine taugb9
      !----------------------------------------------------------------------------
      !
      !     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
      !                             (high key - ch4; high minor - n2o)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng9, ngs8
      use mo_rrlw_kg09,   only : fracrefa, fracrefb, absa, absb, &
           ka_mn2o, kb_mn2o, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      integer :: js, js1, jmn2o, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_mn2o, specparm_mn2o, specmult_mn2o, fmn2o
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor, n2om1, n2om2, absn2o
      real(wp) :: chi_n2o, ratn2o, adjfac, adjcoln2o
      real(wp) :: refrat_planck_a, refrat_m_a
      real(wp) :: tau_major(ng9), tau_major1(ng9)


      ! Minor gas mapping level :
      !     lower - n2o, p = 706.272 mbar, t = 278.94 k
      !     upper - n2o, p = 95.58 mbar, t = 215.7 k

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P = 212 mb
      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)

      ! P = 706.272 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colch4(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._wp*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - AINT((specmult_mn2o)))

          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(9) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(9) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *    &
             (fac000 * absa(ind0,1:ng9)    + &
              fac100 * absa(ind0+1,1:ng9)  + &
              fac200 * absa(ind0+2,1:ng9)  + &
              fac010 * absa(ind0+9,1:ng9)  + &
              fac110 * absa(ind0+10,1:ng9) + &
              fac210 * absa(ind0+11,1:ng9))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng9) + &
              fac100 * absa(ind0,1:ng9)   + &
              fac000 * absa(ind0+1,1:ng9) + &
              fac210 * absa(ind0+8,1:ng9) + &
              fac110 * absa(ind0+9,1:ng9) + &
              fac010 * absa(ind0+10,1:ng9))
          else
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *   &
             (fac000 * absa(ind0,1:ng9)   + &
              fac100 * absa(ind0+1,1:ng9) + &
              fac010 * absa(ind0+9,1:ng9) + &
              fac110 * absa(ind0+10,1:ng9))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng9)    + &
              fac101 * absa(ind1+1,1:ng9)  + &
              fac201 * absa(ind1+2,1:ng9)  + &
              fac011 * absa(ind1+9,1:ng9)  + &
              fac111 * absa(ind1+10,1:ng9) + &
              fac211 * absa(ind1+11,1:ng9))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng9) + &
              fac101 * absa(ind1,1:ng9)   + &
              fac001 * absa(ind1+1,1:ng9) + &
              fac211 * absa(ind1+8,1:ng9) + &
              fac111 * absa(ind1+9,1:ng9) + &
              fac011 * absa(ind1+10,1:ng9))
          else
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac001 * absa(ind1,1:ng9)   + &
              fac101 * absa(ind1+1,1:ng9) + &
              fac011 * absa(ind1+9,1:ng9) + &
              fac111 * absa(ind1+10,1:ng9))
          endif

!$NEC unroll(NG9)
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

            taug(jl,lay,ngs8+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs8+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(9) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(9) + 1
          indm = indminor(jl,lay)
!$NEC unroll(NG9)
          do ig = 1, ng9
            absn2o = kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
            taug(jl,lay,ngs8+ig) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs8+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mn2o = colh2o(jl,lay) + refrat_m_a*colch4(jl,lay)
          specparm_mn2o = MIN(colh2o(jl,lay)/speccomb_mn2o,oneminus)
          specmult_mn2o = 8._wp*specparm_mn2o
          jmn2o = 1 + int(specmult_mn2o)
          fmn2o = ((specmult_mn2o) - AINT((specmult_mn2o)))

          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl= 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(9) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(9) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *    &
             (fac000 * absa(ind0,1:ng9)    + &
              fac100 * absa(ind0+1,1:ng9)  + &
              fac200 * absa(ind0+2,1:ng9)  + &
              fac010 * absa(ind0+9,1:ng9)  + &
              fac110 * absa(ind0+10,1:ng9) + &
              fac210 * absa(ind0+11,1:ng9))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng9) + &
              fac100 * absa(ind0,1:ng9)   + &
              fac000 * absa(ind0+1,1:ng9) + &
              fac210 * absa(ind0+8,1:ng9) + &
              fac110 * absa(ind0+9,1:ng9) + &
              fac010 * absa(ind0+10,1:ng9))
          else
!$NEC unroll(NG9)
            tau_major(1:ng9) = speccomb *   &
             (fac000 * absa(ind0,1:ng9)   + &
              fac100 * absa(ind0+1,1:ng9) + &
              fac010 * absa(ind0+9,1:ng9) + &
              fac110 * absa(ind0+10,1:ng9))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng9)    + &
              fac101 * absa(ind1+1,1:ng9)  + &
              fac201 * absa(ind1+2,1:ng9)  + &
              fac011 * absa(ind1+9,1:ng9)  + &
              fac111 * absa(ind1+10,1:ng9) + &
              fac211 * absa(ind1+11,1:ng9))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng9) + &
              fac101 * absa(ind1,1:ng9)   + &
              fac001 * absa(ind1+1,1:ng9) + &
              fac211 * absa(ind1+8,1:ng9) + &
              fac111 * absa(ind1+9,1:ng9) + &
              fac011 * absa(ind1+10,1:ng9))
          else
!$NEC unroll(NG9)
            tau_major1(1:ng9) = speccomb1 * &
             (fac001 * absa(ind1,1:ng9)   + &
              fac101 * absa(ind1+1,1:ng9) + &
              fac011 * absa(ind1+9,1:ng9) + &
              fac111 * absa(ind1+10,1:ng9))
          endif

!$NEC unroll(NG9)
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

            taug(jl,lay,ngs8+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs8+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          !  In atmospheres where the amount of N2O is too great to be considered
          !  a minor species, adjust the column amount of N2O by an empirical factor
          !  to obtain the proper contribution.
          chi_n2o = coln2o(jl,lay)/(coldry(jl,lay))
          ratn2o = 1.e20_wp*chi_n2o/chi_mls(4,jp(jl,lay)+1)
          if (ratn2o .gt. 1.5_wp) then
            adjfac = 0.5_wp+(ratn2o-0.5_wp)**0.65_wp
            adjcoln2o = adjfac*chi_mls(4,jp(jl,lay)+1)*coldry(jl,lay)*1.e-20_wp
          else
            adjcoln2o = coln2o(jl,lay)
          endif

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(9) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(9) + 1
          indm = indminor(jl,lay)
!$NEC unroll(NG9)
          do ig = 1, ng9
            absn2o = kb_mn2o(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mn2o(indm+1,ig) - kb_mn2o(indm,ig))
            taug(jl,lay,ngs8+ig) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + adjcoln2o*absn2o
            fracs(jl,lay,ngs8+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb9

    !----------------------------------------------------------------------------
    subroutine taugb10
      !----------------------------------------------------------------------------
      !
      !     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng10, ngs9
      use mo_rrlw_kg10,   only : fracrefa, fracrefb, absa, absb, &
           selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf
      real(wp) :: tauself, taufor


      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature.  Below laytrop, the water vapor self-continuum and
      ! foreign continuum is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(10) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(10) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
!$NEC unroll(NG10)
          do ig = 1, ng10
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs9+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor
            fracs(jl,lay,ngs9+ig) = fracrefa(ig)
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(10) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(10) + 1
          indf = indfor(jl,lay)
!$NEC unroll(NG10)
          do ig = 1, ng10
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs9+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(jl,lay,ngs9+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(10) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(10) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
!$NEC unroll(NG10)
          do ig = 1, ng10
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs9+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig))  &
                 + tauself + taufor
            fracs(jl,lay,ngs9+ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(10) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(10) + 1
          indf = indfor(jl,lay)
!$NEC unroll(NG10)
          do ig = 1, ng10
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs9+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) +  &
                 fac11(jl,lay) * absb(ind1+1,ig)) &
                 + taufor
            fracs(jl,lay,ngs9+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb10

    !----------------------------------------------------------------------------
    subroutine taugb11
      !----------------------------------------------------------------------------
      !
      !     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
      !                              (high key - h2o; high minor - o2)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng11, ngs10
      use mo_rrlw_kg11,   only : fracrefa, fracrefb, absa, absb, &
           ka_mo2, kb_mo2, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      real(wp) :: scaleo2, tauself, taufor, tauo2


      ! Minor gas mapping level :
      !     lower - o2, p = 706.2720 mbar, t = 278.94 k
      !     upper - o2, p = 4.758820 mbarm t = 250.85 k

      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature.  Below laytrop, the water vapor self-continuum and
      ! foreign continuum is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(11) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(11) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
            taug(jl,lay,ngs10+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracs(jl,lay,ngs10+ig) = fracrefa(ig)
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(11) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(11) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
            taug(jl,lay,ngs10+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracs(jl,lay,ngs10+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(11) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(11) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (ka_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (ka_mo2(indm+1,ig) - ka_mo2(indm,ig)))
            taug(jl,lay,ngs10+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor &
                 + tauo2
            fracs(jl,lay,ngs10+ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(11) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(11) + 1
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)
          scaleo2 = colo2(jl,lay)*scaleminor(jl,lay)
!$NEC unroll(NG11)
          do ig = 1, ng11
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            tauo2 =  scaleo2 * (kb_mo2(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo2(indm+1,ig) - kb_mo2(indm,ig)))
            taug(jl,lay,ngs10+ig) = colh2o(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))  &
                 + taufor &
                 + tauo2
            fracs(jl,lay,ngs10+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb11

    !----------------------------------------------------------------------------
    subroutine taugb12
      !----------------------------------------------------------------------------
      !
      !     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng12, ngs11
      use mo_rrlw_kg12,   only : fracrefa, absa, &
           selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf
      integer :: js, js1, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor
      real(wp) :: refrat_planck_a
      real(wp) :: tau_major(ng12), tau_major1(ng12)


      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P =   174.164 mb
      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum adn foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      DO lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(12) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(12) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG12)
            tau_major(1:ng12) = speccomb *    &
             (fac000 * absa(ind0,1:ng12)    + &
              fac100 * absa(ind0+1,1:ng12)  + &
              fac200 * absa(ind0+2,1:ng12)  + &
              fac010 * absa(ind0+9,1:ng12)  + &
              fac110 * absa(ind0+10,1:ng12) + &
              fac210 * absa(ind0+11,1:ng12))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG12)
            tau_major(1:ng12) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng12) + &
              fac100 * absa(ind0,1:ng12)   + &
              fac000 * absa(ind0+1,1:ng12) + &
              fac210 * absa(ind0+8,1:ng12) + &
              fac110 * absa(ind0+9,1:ng12) + &
              fac010 * absa(ind0+10,1:ng12))
          else
!$NEC unroll(NG12)
            tau_major(1:ng12) = speccomb *   &
             (fac000 * absa(ind0,1:ng12)   + &
              fac100 * absa(ind0+1,1:ng12) + &
              fac010 * absa(ind0+9,1:ng12) + &
              fac110 * absa(ind0+10,1:ng12))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG12)
            tau_major1(1:ng12) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng12)    + &
              fac101 * absa(ind1+1,1:ng12)  + &
              fac201 * absa(ind1+2,1:ng12)  + &
              fac011 * absa(ind1+9,1:ng12)  + &
              fac111 * absa(ind1+10,1:ng12) + &
              fac211 * absa(ind1+11,1:ng12))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG12)
            tau_major1(1:ng12) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng12) + &
              fac101 * absa(ind1,1:ng12)   + &
              fac001 * absa(ind1+1,1:ng12) + &
              fac211 * absa(ind1+8,1:ng12) + &
              fac111 * absa(ind1+9,1:ng12) + &
              fac011 * absa(ind1+10,1:ng12))
          else
!$NEC unroll(NG12)
            tau_major1(1:ng12) = speccomb1 * &
             (fac001 * absa(ind1,1:ng12)   + &
              fac101 * absa(ind1+1,1:ng12) + &
              fac011 * absa(ind1+9,1:ng12) + &
              fac111 * absa(ind1+10,1:ng12))
          endif

!$NEC unroll(NG12)
          do ig = 1, ng12
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,lay,ngs11+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,lay,ngs11+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      ENDDO

      ! Upper atmosphere loop
      do ig = 1, ng12
        do lay = laytrop_max+1, nlayers
          do jl = 1, kproma
            taug(jl,lay,ngs11+ig) = 0.0_wp
            fracs(jl,lay,ngs11+ig) = 0.0_wp
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(12) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(12) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG12)
            tau_major(1:ng12) = speccomb *    &
             (fac000 * absa(ind0,1:ng12)    + &
              fac100 * absa(ind0+1,1:ng12)  + &
              fac200 * absa(ind0+2,1:ng12)  + &
              fac010 * absa(ind0+9,1:ng12)  + &
              fac110 * absa(ind0+10,1:ng12) + &
              fac210 * absa(ind0+11,1:ng12))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG12)
            tau_major(1:ng12) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng12) + &
              fac100 * absa(ind0,1:ng12)   + &
              fac000 * absa(ind0+1,1:ng12) + &
              fac210 * absa(ind0+8,1:ng12) + &
              fac110 * absa(ind0+9,1:ng12) + &
              fac010 * absa(ind0+10,1:ng12))
          else
!$NEC unroll(NG12)
            tau_major(1:ng12) = speccomb *   &
             (fac000 * absa(ind0,1:ng12)   + &
              fac100 * absa(ind0+1,1:ng12) + &
              fac010 * absa(ind0+9,1:ng12) + &
              fac110 * absa(ind0+10,1:ng12))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG12)
            tau_major1(1:ng12) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng12)    + &
              fac101 * absa(ind1+1,1:ng12)  + &
              fac201 * absa(ind1+2,1:ng12)  + &
              fac011 * absa(ind1+9,1:ng12)  + &
              fac111 * absa(ind1+10,1:ng12) + &
              fac211 * absa(ind1+11,1:ng12))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG12)
            tau_major1(1:ng12) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng12) + &
              fac101 * absa(ind1,1:ng12)   + &
              fac001 * absa(ind1+1,1:ng12) + &
              fac211 * absa(ind1+8,1:ng12) + &
              fac111 * absa(ind1+9,1:ng12) + &
              fac011 * absa(ind1+10,1:ng12))
          else
!$NEC unroll(NG12)
            tau_major1(1:ng12) = speccomb1 * &
             (fac001 * absa(ind1,1:ng12)   + &
              fac101 * absa(ind1+1,1:ng12) + &
              fac011 * absa(ind1+9,1:ng12) + &
              fac111 * absa(ind1+10,1:ng12))
          endif

!$NEC unroll(NG12)
          do ig = 1, ng12
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor = forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,lay,ngs11+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,lay,ngs11+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0

        do ig = 1, ng12
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)

            taug(jl,lay,ngs11+ig) = 0.0_wp
            fracs(jl,lay,ngs11+ig) = 0.0_wp
          enddo
        enddo

      enddo

    end subroutine taugb12

    !----------------------------------------------------------------------------
    subroutine taugb13
      !----------------------------------------------------------------------------
      !
      !     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng13, ngs12
      use mo_rrlw_kg13,   only : fracrefa, fracrefb, absa, &
           ka_mco2, ka_mco, kb_mo3, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      integer :: js, js1, jmco2, jmco, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_mco2, specparm_mco2, specmult_mco2, fmco2
      real(wp) :: speccomb_mco, specparm_mco, specmult_mco, fmco
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor, co2m1, co2m2, absco2
      real(wp) :: com1, com2, absco, abso3
      real(wp) :: chi_co2, ratco2, adjfac, adjcolco2
      real(wp) :: refrat_planck_a, refrat_m_a, refrat_m_a3
      real(wp) :: tau_major(ng13), tau_major1(ng13)


      ! Minor gas mapping levels :
      !     lower - co2, p = 1053.63 mb, t = 294.2 k
      !     lower - co, p = 706 mb, t = 278.94 k
      !     upper - o3, p = 95.5835 mb, t = 215.7 k

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower/upper atmosphere.

      ! P = 473.420 mb (Level 5)
      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)

      ! P = 1053. (Level 1)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)

      ! P = 706. (Level 3)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2on2o(jl,lay)*coln2o(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2on2o_1(jl,lay)*coln2o(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*coln2o(jl,lay)
          specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._wp*specparm_mco2
          jmco2 = 1 + int(specmult_mco2)
          fmco2 = ((specmult_mco2) - AINT((specmult_mco2)))

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/3.55e-4_wp
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.68_wp
            adjcolco2 = adjfac*3.55e-4_wp*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_mco = colh2o(jl,lay) + refrat_m_a3*coln2o(jl,lay)
          specparm_mco = MIN(colh2o(jl,lay)/speccomb_mco,oneminus)
          specmult_mco = 8._wp*specparm_mco
          jmco = 1 + int(specmult_mco)
          fmco = ((specmult_mco) - AINT((specmult_mco)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*coln2o(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(13) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(13) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *    &
             (fac000 * absa(ind0,1:ng13)    + &
              fac100 * absa(ind0+1,1:ng13)  + &
              fac200 * absa(ind0+2,1:ng13)  + &
              fac010 * absa(ind0+9,1:ng13)  + &
              fac110 * absa(ind0+10,1:ng13) + &
              fac210 * absa(ind0+11,1:ng13))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng13) + &
              fac100 * absa(ind0,1:ng13)   + &
              fac000 * absa(ind0+1,1:ng13) + &
              fac210 * absa(ind0+8,1:ng13) + &
              fac110 * absa(ind0+9,1:ng13) + &
              fac010 * absa(ind0+10,1:ng13))
          else
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *   &
             (fac000 * absa(ind0,1:ng13)   + &
              fac100 * absa(ind0+1,1:ng13) + &
              fac010 * absa(ind0+9,1:ng13) + &
              fac110 * absa(ind0+10,1:ng13))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng13)    + &
              fac101 * absa(ind1+1,1:ng13)  + &
              fac201 * absa(ind1+2,1:ng13)  + &
              fac011 * absa(ind1+9,1:ng13)  + &
              fac111 * absa(ind1+10,1:ng13) + &
              fac211 * absa(ind1+11,1:ng13))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng13) + &
              fac101 * absa(ind1,1:ng13)   + &
              fac001 * absa(ind1+1,1:ng13) + &
              fac211 * absa(ind1+8,1:ng13) + &
              fac111 * absa(ind1+9,1:ng13) + &
              fac011 * absa(ind1+10,1:ng13))
          else
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac001 * absa(ind1,1:ng13)   + &
              fac101 * absa(ind1+1,1:ng13) + &
              fac011 * absa(ind1+9,1:ng13) + &
              fac111 * absa(ind1+10,1:ng13))
          endif

!$NEC unroll(NG13)
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

            taug(jl,lay,ngs12+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(jl,lay)*absco
            fracs(jl,lay,ngs12+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          indm = indminor(jl,lay)
!$NEC unroll(NG13)
          do ig = 1, ng13
            abso3 = kb_mo3(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
            taug(jl,lay,ngs12+ig) = colo3(jl,lay)*abso3
            fracs(jl,lay,ngs12+ig) =  fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2on2o(jl,lay)*coln2o(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2on2o_1(jl,lay)*coln2o(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mco2 = colh2o(jl,lay) + refrat_m_a*coln2o(jl,lay)
          specparm_mco2 = MIN(colh2o(jl,lay)/speccomb_mco2,oneminus)
          specmult_mco2 = 8._wp*specparm_mco2
          jmco2 = 1 + int(specmult_mco2)
          fmco2 = ((specmult_mco2) - AINT((specmult_mco2)))

          !  In atmospheres where the amount of CO2 is too great to be considered
          !  a minor species, adjust the column amount of CO2 by an empirical factor
          !  to obtain the proper contribution.
          chi_co2 = colco2(jl,lay)/(coldry(jl,lay))
          ratco2 = 1.e20_wp*chi_co2/3.55e-4_wp
          if (ratco2 .gt. 3.0_wp) then
            adjfac = 2.0_wp+(ratco2-2.0_wp)**0.68_wp
            adjcolco2 = adjfac*3.55e-4_wp*coldry(jl,lay)*1.e-20_wp
          else
            adjcolco2 = colco2(jl,lay)
          endif

          speccomb_mco = colh2o(jl,lay) + refrat_m_a3*coln2o(jl,lay)
          specparm_mco = MIN(colh2o(jl,lay)/speccomb_mco,oneminus)
          specmult_mco = 8._wp*specparm_mco
          jmco = 1 + int(specmult_mco)
          fmco = ((specmult_mco) - AINT((specmult_mco)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*coln2o(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(13) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(13) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *    &
             (fac000 * absa(ind0,1:ng13)    + &
              fac100 * absa(ind0+1,1:ng13)  + &
              fac200 * absa(ind0+2,1:ng13)  + &
              fac010 * absa(ind0+9,1:ng13)  + &
              fac110 * absa(ind0+10,1:ng13) + &
              fac210 * absa(ind0+11,1:ng13))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng13) + &
              fac100 * absa(ind0,1:ng13)   + &
              fac000 * absa(ind0+1,1:ng13) + &
              fac210 * absa(ind0+8,1:ng13) + &
              fac110 * absa(ind0+9,1:ng13) + &
              fac010 * absa(ind0+10,1:ng13))
          else
!$NEC unroll(NG13)
            tau_major(1:ng13) = speccomb *   &
             (fac000 * absa(ind0,1:ng13)   + &
              fac100 * absa(ind0+1,1:ng13) + &
              fac010 * absa(ind0+9,1:ng13) + &
              fac110 * absa(ind0+10,1:ng13))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng13)    + &
              fac101 * absa(ind1+1,1:ng13)  + &
              fac201 * absa(ind1+2,1:ng13)  + &
              fac011 * absa(ind1+9,1:ng13)  + &
              fac111 * absa(ind1+10,1:ng13) + &
              fac211 * absa(ind1+11,1:ng13))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng13) + &
              fac101 * absa(ind1,1:ng13)   + &
              fac001 * absa(ind1+1,1:ng13) + &
              fac211 * absa(ind1+8,1:ng13) + &
              fac111 * absa(ind1+9,1:ng13) + &
              fac011 * absa(ind1+10,1:ng13))
          else
!$NEC unroll(NG13)
            tau_major1(1:ng13) = speccomb1 * &
             (fac001 * absa(ind1,1:ng13)   + &
              fac101 * absa(ind1+1,1:ng13) + &
              fac011 * absa(ind1+9,1:ng13) + &
              fac111 * absa(ind1+10,1:ng13))
          endif

!$NEC unroll(NG13)
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

            taug(jl,lay,ngs12+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + adjcolco2*absco2 &
                 + colco(jl,lay)*absco
            fracs(jl,lay,ngs12+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          indm = indminor(jl,lay)
!$NEC unroll(NG13)
          do ig = 1, ng13
            abso3 = kb_mo3(indm,ig) + minorfrac(jl,lay) * &
                 (kb_mo3(indm+1,ig) - kb_mo3(indm,ig))
            taug(jl,lay,ngs12+ig) = colo3(jl,lay)*abso3
            fracs(jl,lay,ngs12+ig) =  fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb13

    !----------------------------------------------------------------------------
    subroutine taugb14
      !----------------------------------------------------------------------------
      !
      !     band 14:  2250-2380 cm-1 (low - co2; high - co2)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng14, ngs13
      use mo_rrlw_kg14,   only : fracrefa, fracrefb, absa, absb, &
           selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf
      real(wp) :: tauself, taufor


      ! Compute the optical depth by interpolating in ln(pressure) and
      ! temperature.  Below laytrop, the water vapor self-continuum
      ! and foreign continuum is interpolated (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(14) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(14) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
!$NEC unroll(NG14)
          do ig = 1, ng14
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs13+ig) = colco2(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor
            fracs(jl,lay,ngs13+ig) = fracrefa(ig)
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(14) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(14) + 1
!$NEC unroll(NG14)
          do ig = 1, ng14
            taug(jl,lay,ngs13+ig) = colco2(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))
            fracs(jl,lay,ngs13+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(14) + 1
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(14) + 1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
!$NEC unroll(NG14)
          do ig = 1, ng14
            tauself = selffac(jl,lay) * (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            taug(jl,lay,ngs13+ig) = colco2(jl,lay) * &
                 (fac00(jl,lay) * absa(ind0,ig) + &
                 fac10(jl,lay) * absa(ind0+1,ig) + &
                 fac01(jl,lay) * absa(ind1,ig) + &
                 fac11(jl,lay) * absa(ind1+1,ig)) &
                 + tauself + taufor
            fracs(jl,lay,ngs13+ig) = fracrefa(ig)
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(14) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(14) + 1
!$NEC unroll(NG14)
          do ig = 1, ng14
            taug(jl,lay,ngs13+ig) = colco2(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))
            fracs(jl,lay,ngs13+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb14

    !----------------------------------------------------------------------------
    subroutine taugb15
      !----------------------------------------------------------------------------
      !
      !     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
      !                              (high - nothing)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng15, ngs14
      use mo_rrlw_kg15,   only : fracrefa, absa, &
           ka_mn2, selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf, indm
      integer :: js, js1, jmn2, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_mn2, specparm_mn2, specmult_mn2, fmn2
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: scalen2, tauself, taufor, n2m1, n2m2, taun2
      real(wp) :: refrat_planck_a, refrat_m_a
      real(wp) :: tau_major(ng15), tau_major1(ng15)


      ! Minor gas mapping level :
      !     Lower - Nitrogen Continuum, P = 1053., T = 294.

      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.
      ! P = 1053. mb (Level 1)
      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)

      ! P = 1053.
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature, and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = coln2o(jl,lay) + rat_n2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(coln2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = coln2o(jl,lay) + rat_n2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(coln2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mn2 = coln2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mn2 = MIN(coln2o(jl,lay)/speccomb_mn2,oneminus)
          specmult_mn2 = 8._wp*specparm_mn2
          jmn2 = 1 + int(specmult_mn2)
          fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))

          speccomb_planck = coln2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(coln2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(15) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(15) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          scalen2 = colbrd(jl,lay)*scaleminor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *    &
             (fac000 * absa(ind0,1:ng15)    + &
              fac100 * absa(ind0+1,1:ng15)  + &
              fac200 * absa(ind0+2,1:ng15)  + &
              fac010 * absa(ind0+9,1:ng15)  + &
              fac110 * absa(ind0+10,1:ng15) + &
              fac210 * absa(ind0+11,1:ng15))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng15) + &
              fac100 * absa(ind0,1:ng15)   + &
              fac000 * absa(ind0+1,1:ng15) + &
              fac210 * absa(ind0+8,1:ng15) + &
              fac110 * absa(ind0+9,1:ng15) + &
              fac010 * absa(ind0+10,1:ng15))
          else
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *   &
             (fac000 * absa(ind0,1:ng15)   + &
              fac100 * absa(ind0+1,1:ng15) + &
              fac010 * absa(ind0+9,1:ng15) + &
              fac110 * absa(ind0+10,1:ng15))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng15)    + &
              fac101 * absa(ind1+1,1:ng15)  + &
              fac201 * absa(ind1+2,1:ng15)  + &
              fac011 * absa(ind1+9,1:ng15)  + &
              fac111 * absa(ind1+10,1:ng15) + &
              fac211 * absa(ind1+11,1:ng15))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng15) + &
              fac101 * absa(ind1,1:ng15)   + &
              fac001 * absa(ind1+1,1:ng15) + &
              fac211 * absa(ind1+8,1:ng15) + &
              fac111 * absa(ind1+9,1:ng15) + &
              fac011 * absa(ind1+10,1:ng15))
          else
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 * &
             (fac001 * absa(ind1,1:ng15)   + &
              fac101 * absa(ind1+1,1:ng15) + &
              fac011 * absa(ind1+9,1:ng15) + &
              fac111 * absa(ind1+10,1:ng15))
          endif

!$NEC unroll(NG15)
          do ig = 1, ng15
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
            n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(jl,lay) * (n2m2 - n2m1))

            taug(jl,lay,ngs14+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + taun2
            fracs(jl,lay,ngs14+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo
      enddo

      ! Upper atmosphere loop
      do ig = 1, ng15
        do lay = laytrop_max+1, nlayers
          do jl = 1, kproma
            taug(jl,lay,ngs14+ig) = 0.0_wp
            fracs(jl,lay,ngs14+ig) = 0.0_wp
          enddo
        enddo
      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = coln2o(jl,lay) + rat_n2oco2(jl,lay)*colco2(jl,lay)
          specparm = MIN(coln2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = coln2o(jl,lay) + rat_n2oco2_1(jl,lay)*colco2(jl,lay)
          specparm1 = MIN(coln2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_mn2 = coln2o(jl,lay) + refrat_m_a*colco2(jl,lay)
          specparm_mn2 = MIN(coln2o(jl,lay)/speccomb_mn2,oneminus)
          specmult_mn2 = 8._wp*specparm_mn2
          jmn2 = 1 + int(specmult_mn2)
          fmn2 = ((specmult_mn2) - AINT((specmult_mn2)))

          speccomb_planck = coln2o(jl,lay)+refrat_planck_a*colco2(jl,lay)
          specparm_planck = MIN(coln2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(15) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(15) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)
          indm = indminor(jl,lay)

          scalen2 = colbrd(jl,lay)*scaleminor(jl,lay)

           if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *    &
             (fac000 * absa(ind0,1:ng15)    + &
              fac100 * absa(ind0+1,1:ng15)  + &
              fac200 * absa(ind0+2,1:ng15)  + &
              fac010 * absa(ind0+9,1:ng15)  + &
              fac110 * absa(ind0+10,1:ng15) + &
              fac210 * absa(ind0+11,1:ng15))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng15) + &
              fac100 * absa(ind0,1:ng15)   + &
              fac000 * absa(ind0+1,1:ng15) + &
              fac210 * absa(ind0+8,1:ng15) + &
              fac110 * absa(ind0+9,1:ng15) + &
              fac010 * absa(ind0+10,1:ng15))
          else
!$NEC unroll(NG15)
            tau_major(1:ng15) = speccomb *   &
             (fac000 * absa(ind0,1:ng15)   + &
              fac100 * absa(ind0+1,1:ng15) + &
              fac010 * absa(ind0+9,1:ng15) + &
              fac110 * absa(ind0+10,1:ng15))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 *  &
             (fac001 * absa(ind1,1:ng15)    + &
              fac101 * absa(ind1+1,1:ng15)  + &
              fac201 * absa(ind1+2,1:ng15)  + &
              fac011 * absa(ind1+9,1:ng15)  + &
              fac111 * absa(ind1+10,1:ng15) + &
              fac211 * absa(ind1+11,1:ng15))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng15) + &
              fac101 * absa(ind1,1:ng15)   + &
              fac001 * absa(ind1+1,1:ng15) + &
              fac211 * absa(ind1+8,1:ng15) + &
              fac111 * absa(ind1+9,1:ng15) + &
              fac011 * absa(ind1+10,1:ng15))
          else
!$NEC unroll(NG15)
            tau_major1(1:ng15) = speccomb1 * &
             (fac001 * absa(ind1,1:ng15)   + &
              fac101 * absa(ind1+1,1:ng15) + &
              fac011 * absa(ind1+9,1:ng15) + &
              fac111 * absa(ind1+10,1:ng15))
          endif

!$NEC unroll(NG15)
          do ig = 1, ng15
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))
            n2m1 = ka_mn2(jmn2,indm,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm,ig) - ka_mn2(jmn2,indm,ig))
            n2m2 = ka_mn2(jmn2,indm+1,ig) + fmn2 * &
                 (ka_mn2(jmn2+1,indm+1,ig) - ka_mn2(jmn2,indm+1,ig))
            taun2 = scalen2 * (n2m1 + minorfrac(jl,lay) * (n2m2 - n2m1))

            taug(jl,lay,ngs14+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor &
                 + taun2
            fracs(jl,lay,ngs14+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0

        do ig = 1, ng15
!$NEC ivdep
          do ixp = 1, ixc0
            jl = ixhigh(ixp,lay)

            taug(jl,lay,ngs14+ig) = 0.0_wp
            fracs(jl,lay,ngs14+ig) = 0.0_wp
          enddo
        enddo

      enddo

    end subroutine taugb15

    !----------------------------------------------------------------------------
    subroutine taugb16
      !----------------------------------------------------------------------------
      !
      !     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
      !----------------------------------------------------------------------------

      ! ------- Modules -------

      use mo_lrtm_par,    only : ng16, ngs15
      use mo_rrlw_kg16,   only : fracrefa, fracrefb, absa, absb, &
           selfref, forref

      ! ------- Declarations -------

      ! Local
      integer :: lay, ig, ixc0, ixp, jl
      integer :: ind0, ind1, inds, indf
      integer :: js, js1, jpl
      real(wp) :: speccomb, specparm, specmult, fs
      real(wp) :: speccomb1, specparm1, specmult1, fs1
      real(wp) :: speccomb_planck, specparm_planck, specmult_planck, fpl
      real(wp) :: p, p4, fk0, fk1, fk2
      real(wp) :: fac000, fac100, fac200
      real(wp) :: fac010, fac110, fac210
      real(wp) :: fac001, fac101, fac201
      real(wp) :: fac011, fac111, fac211
      real(wp) :: tauself, taufor
      real(wp) :: refrat_planck_a
      real(wp) :: tau_major(ng16), tau_major1(ng16)


      ! Calculate reference ratio to be used in calculation of Planck
      ! fraction in lower atmosphere.

      ! P = 387. mb (Level 6)
      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)

      ! Compute the optical depth by interpolating in ln(pressure),
      ! temperature,and appropriate species.  Below laytrop, the water
      ! vapor self-continuum and foreign continuum is interpolated
      ! (in temperature) separately.

      ! Lower atmosphere loop
      do lay = 1, laytrop_min
        do jl = 1, kproma

          speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(16) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(16) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *    &
             (fac000 * absa(ind0,1:ng16)    + &
              fac100 * absa(ind0+1,1:ng16)  + &
              fac200 * absa(ind0+2,1:ng16)  + &
              fac010 * absa(ind0+9,1:ng16)  + &
              fac110 * absa(ind0+10,1:ng16) + &
              fac210 * absa(ind0+11,1:ng16))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng16) + &
              fac100 * absa(ind0,1:ng16)   + &
              fac000 * absa(ind0+1,1:ng16) + &
              fac210 * absa(ind0+8,1:ng16) + &
              fac110 * absa(ind0+9,1:ng16) + &
              fac010 * absa(ind0+10,1:ng16))
          else
!$NEC unroll(NG16)
             tau_major(1:ng16) = speccomb *   &
              (fac000 * absa(ind0,1:ng16)   + &
               fac100 * absa(ind0+1,1:ng16) + &
               fac010 * absa(ind0+9,1:ng16) + &
               fac110 * absa(ind0+10,1:ng16))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG16)
             tau_major1(1:ng16) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng16)    + &
               fac101 * absa(ind1+1,1:ng16)  + &
               fac201 * absa(ind1+2,1:ng16)  + &
               fac011 * absa(ind1+9,1:ng16)  + &
               fac111 * absa(ind1+10,1:ng16) + &
               fac211 * absa(ind1+11,1:ng16))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng16) + &
              fac101 * absa(ind1,1:ng16)   + &
              fac001 * absa(ind1+1,1:ng16) + &
              fac211 * absa(ind1+8,1:ng16) + &
              fac111 * absa(ind1+9,1:ng16) + &
              fac011 * absa(ind1+10,1:ng16))
          else
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac001 * absa(ind1,1:ng16)   + &
              fac101 * absa(ind1+1,1:ng16) + &
              fac011 * absa(ind1+9,1:ng16) + &
              fac111 * absa(ind1+10,1:ng16))
          endif

!$NEC unroll(NG16)
          do ig = 1, ng16
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,lay,ngs15+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,lay,ngs15+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

      enddo

      ! Upper atmosphere loop
      do lay = laytrop_max+1, nlayers
        do jl = 1, kproma

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(16) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(16) + 1
!$NEC unroll(NG16)
          do ig = 1, ng16
            taug(jl,lay,ngs15+ig) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))
            fracs(jl,lay,ngs15+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

      IF (laytrop_max == laytrop_min) RETURN
      ! Mixed loop
      ! Lower atmosphere part
      do lay = laytrop_min+1, laytrop_max
        ixc0 = ixc(lay)

!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixlow(ixp,lay)

          speccomb = colh2o(jl,lay) + rat_h2och4(jl,lay)*colch4(jl,lay)
          specparm = MIN(colh2o(jl,lay)/speccomb,oneminus)
          specmult = 8._wp*(specparm)
          js = 1 + int(specmult)
          fs = ((specmult) - AINT((specmult)))

          speccomb1 = colh2o(jl,lay) + rat_h2och4_1(jl,lay)*colch4(jl,lay)
          specparm1 = MIN(colh2o(jl,lay)/speccomb1,oneminus)
          specmult1 = 8._wp*(specparm1)
          js1 = 1 + int(specmult1)
          fs1 = ((specmult1) - AINT((specmult1)))

          speccomb_planck = colh2o(jl,lay)+refrat_planck_a*colch4(jl,lay)
          specparm_planck = MIN(colh2o(jl,lay)/speccomb_planck,oneminus)
          specmult_planck = 8._wp*specparm_planck
          jpl = 1 + int(specmult_planck)
          fpl = ((specmult_planck) - AINT((specmult_planck)))

          ind0 = ((jp(jl,lay)-1)*5+(jt(jl,lay)-1))*nspa(16) + js
          ind1 = (jp(jl,lay)*5+(jt1(jl,lay)-1))*nspa(16) + js1
          inds = indself(jl,lay)
          indf = indfor(jl,lay)

          if (specparm .lt. 0.125_wp) then
            p = fs - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else if (specparm .gt. 0.875_wp) then
            p = -fs
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac000 = fk0*fac00(jl,lay)
            fac100 = fk1*fac00(jl,lay)
            fac200 = fk2*fac00(jl,lay)
            fac010 = fk0*fac10(jl,lay)
            fac110 = fk1*fac10(jl,lay)
            fac210 = fk2*fac10(jl,lay)
          else
            fac000 = (1._wp - fs) * fac00(jl,lay)
            fac010 = (1._wp - fs) * fac10(jl,lay)
            fac100 = fs * fac00(jl,lay)
            fac110 = fs * fac10(jl,lay)
            fac200 = 0._wp
            fac210 = 0._wp
          endif

          if (specparm1 .lt. 0.125_wp) then
            p = fs1 - 1._wp
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else if (specparm1 .gt. 0.875_wp) then
            p = -fs1
            p4 = p**4
            fk0 = p4
            fk1 = 1._wp - p - 2.0_wp*p4
            fk2 = p + p4
            fac001 = fk0*fac01(jl,lay)
            fac101 = fk1*fac01(jl,lay)
            fac201 = fk2*fac01(jl,lay)
            fac011 = fk0*fac11(jl,lay)
            fac111 = fk1*fac11(jl,lay)
            fac211 = fk2*fac11(jl,lay)
          else
            fac001 = (1._wp - fs1) * fac01(jl,lay)
            fac011 = (1._wp - fs1) * fac11(jl,lay)
            fac101 = fs1 * fac01(jl,lay)
            fac111 = fs1 * fac11(jl,lay)
            fac201 = 0._wp
            fac211 = 0._wp
          endif

          if (specparm .lt. 0.125_wp) then
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *    &
             (fac000 * absa(ind0,1:ng16)    + &
              fac100 * absa(ind0+1,1:ng16)  + &
              fac200 * absa(ind0+2,1:ng16)  + &
              fac010 * absa(ind0+9,1:ng16)  + &
              fac110 * absa(ind0+10,1:ng16) + &
              fac210 * absa(ind0+11,1:ng16))
          else if (specparm .gt. 0.875_wp) then
!$NEC unroll(NG16)
            tau_major(1:ng16) = speccomb *   &
             (fac200 * absa(ind0-1,1:ng16) + &
              fac100 * absa(ind0,1:ng16)   + &
              fac000 * absa(ind0+1,1:ng16) + &
              fac210 * absa(ind0+8,1:ng16) + &
              fac110 * absa(ind0+9,1:ng16) + &
              fac010 * absa(ind0+10,1:ng16))
          else
!$NEC unroll(NG16)
             tau_major(1:ng16) = speccomb *   &
              (fac000 * absa(ind0,1:ng16)   + &
               fac100 * absa(ind0+1,1:ng16) + &
               fac010 * absa(ind0+9,1:ng16) + &
               fac110 * absa(ind0+10,1:ng16))
          endif

          if (specparm1 .lt. 0.125_wp) then
!$NEC unroll(NG16)
             tau_major1(1:ng16) = speccomb1 *  &
              (fac001 * absa(ind1,1:ng16)    + &
               fac101 * absa(ind1+1,1:ng16)  + &
               fac201 * absa(ind1+2,1:ng16)  + &
               fac011 * absa(ind1+9,1:ng16)  + &
               fac111 * absa(ind1+10,1:ng16) + &
               fac211 * absa(ind1+11,1:ng16))
          else if (specparm1 .gt. 0.875_wp) then
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac201 * absa(ind1-1,1:ng16) + &
              fac101 * absa(ind1,1:ng16)   + &
              fac001 * absa(ind1+1,1:ng16) + &
              fac211 * absa(ind1+8,1:ng16) + &
              fac111 * absa(ind1+9,1:ng16) + &
              fac011 * absa(ind1+10,1:ng16))
          else
!$NEC unroll(NG16)
            tau_major1(1:ng16) = speccomb1 * &
             (fac001 * absa(ind1,1:ng16)   + &
              fac101 * absa(ind1+1,1:ng16) + &
              fac011 * absa(ind1+9,1:ng16) + &
              fac111 * absa(ind1+10,1:ng16))
          endif

!$NEC unroll(NG16)
          do ig = 1, ng16
            tauself = selffac(jl,lay)* (selfref(inds,ig) + selffrac(jl,lay) * &
                 (selfref(inds+1,ig) - selfref(inds,ig)))
            taufor =  forfac(jl,lay) * (forref(indf,ig) + forfrac(jl,lay) * &
                 (forref(indf+1,ig) - forref(indf,ig)))

            taug(jl,lay,ngs15+ig) = tau_major(ig) + tau_major1(ig) &
                 + tauself + taufor
            fracs(jl,lay,ngs15+ig) = fracrefa(ig,jpl) + fpl * &
                 (fracrefa(ig,jpl+1)-fracrefa(ig,jpl))
          enddo
        enddo

        ! Upper atmosphere part
        ixc0 = kproma - ixc0
!$NEC ivdep
        do ixp = 1, ixc0
          jl = ixhigh(ixp,lay)

          ind0 = ((jp(jl,lay)-13)*5+(jt(jl,lay)-1))*nspb(16) + 1
          ind1 = ((jp(jl,lay)-12)*5+(jt1(jl,lay)-1))*nspb(16) + 1
!$NEC unroll(NG16)
          do ig = 1, ng16
            taug(jl,lay,ngs15+ig) = colch4(jl,lay) * &
                 (fac00(jl,lay) * absb(ind0,ig) + &
                 fac10(jl,lay) * absb(ind0+1,ig) + &
                 fac01(jl,lay) * absb(ind1,ig) + &
                 fac11(jl,lay) * absb(ind1+1,ig))
            fracs(jl,lay,ngs15+ig) = fracrefb(ig)
          enddo
        enddo

      enddo

    end subroutine taugb16

  end subroutine lrtm_taumol

end module mo_lrtm_taumol

