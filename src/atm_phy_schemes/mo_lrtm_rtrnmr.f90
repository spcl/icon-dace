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


MODULE mo_lrtm_rtrnmr

  ! ------- Modules -------

  USE mo_kind,             ONLY : wp

  USE mo_math_constants,   ONLY : pi

  USE mo_lrtm_par,         ONLY : nbndlw, delwave, ngs
  USE mo_lrtm_setup,       ONLY : ntbl, bpade, exp_tbl, tfn_tbl
  USE mo_lrtm_setup,       ONLY : tau_tbl

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: lrtm_rtrnmr

  REAL(wp), PARAMETER :: fluxfac = 2.0e+04_wp * pi

  REAL(wp), PARAMETER :: tblint = REAL(ntbl, wp)

CONTAINS

  !-----------------------------------------------------------------------------
  SUBROUTINE lrtm_rtrnmr(&
    & kproma, nlayers, istart, iend, iout, semiss, ncbands, &
    & cldfrac, taucloud, planklay, planklev, plankbnd, &
    & pwvcm, fracs, taut, &
    & totuflux, totdflux, fnet, &
    & totuclfl, totdclfl, fnetc, &
    & idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
    !-----------------------------------------------------------------------------
    !
    !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
    !  Revision for GCMs:  Michael J. Iacono; October, 2002
    !  Revision for F90:  Michael J. Iacono; June, 2006
    !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
    !  Version for In-Order architectures: Thomas Jahns, DKRZ 2016
    !
    !  This program calculates the upward fluxes, downward fluxes, and
    !  heating rates for an arbitrary clear or cloudy atmosphere.  The input
    !  to this program is the atmospheric profile, all Planck function
    !  information, and the cloud fraction by layer.  A variable diffusivity
    !  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9
    !  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of
    !  the column water vapor, and other bands use a value of 1.66.  The Gaussian
    !  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that
    !  use of the emissivity angle for the flux integration can cause errors of
    !  1 to 4 W/m2 within cloudy layers.
    !  Clouds are treated with a maximum-random cloud overlap method.
    !  This subroutine also provides the optional capability to calculate
    !  the derivative of upward flux respect to surface temperature using
    !  the pre-tabulated derivative of the Planck function with respect to
    !  temperature integrated over each spectral band.
    !
    !  Setting the preprocessor switch LRTM_FULL_VECTORIZATION selects
    !  an equivalent (bit-reproducible results with NAG Fortran)
    !  algorithm that vectorizes much more thoroughly and gives
    !  parallelization-invariant results on platforms with FMA where
    !  the divergent branches of the default implementation lead to
    !  different contractions per branch. But it's only faster on
    !  architectures where branches and loop overhead are relatively
    !  expensive compared to computation and memory accesses. The
    !  alternative, incompletely-vectorized version was found to run
    !  faster on Intel Sandy Bridge and Ivy Bridge Xeons and about
    !  equally fast on later AVX2 CPUs (Haswell/Broadwell).

    !***************************************************************************

    ! ------- Declarations -------

    ! ----- Input -----
    INTEGER, INTENT(in) :: kproma          ! number of columns
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, INTENT(in) :: istart          ! beginning band of calculation
    INTEGER, INTENT(in) :: iend            ! ending band of calculation
    INTEGER, INTENT(in) :: iout            ! output option flag

    ! Atmosphere
    REAL(wp), INTENT(in) :: pwvcm(kproma)           ! precipitable water vapor (cm)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(in), CONTIGUOUS :: semiss(:,:)        ! lw surface emissivity
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in), CONTIGUOUS :: planklay(:,:,:)      !
    !    Dimensions: (nlayers,nbndlw)
    REAL(wp), INTENT(in), CONTIGUOUS :: planklev(:,0:,:)     !
    !    Dimensions: (0:nlayers,nbndlw)
    REAL(wp), INTENT(in), CONTIGUOUS :: plankbnd(:,:)        !
    !    Dimensions: (nbndlw)
    REAL(wp), INTENT(in), CONTIGUOUS :: fracs(:,:,:)         !
    !    Dimensions: (kproma,nlayers,ngptw)
    REAL(wp), INTENT(in), CONTIGUOUS :: taut(:,:,:)          ! gaseous + aerosol optical depths
    !    Dimensions: (kproma,nlayers,ngptlw)

    ! Clouds
    INTEGER, INTENT(in) :: ncbands(kproma)          ! number of cloud spectral bands
    REAL(wp), INTENT(in), CONTIGUOUS :: cldfrac(:,:)       ! layer cloud fraction
    !    Dimensions: (nlayers)
    REAL(wp), INTENT(in), CONTIGUOUS :: taucloud(:,:,:)      ! layer cloud optical depth
    !    Dimensions: (nlayers,nbndlw)
    INTEGER, INTENT(in) :: idrv            ! flag for calculation of dF/dt from
    ! Planck derivative [0=off, 1=on]
    REAL(wp), INTENT(in), CONTIGUOUS :: dplankbnd_dt(:,:)    ! derivative of Planck function wrt temp
    !    Dimensions: (nbndlw)

    ! ----- Output -----
    REAL(wp), INTENT(out), CONTIGUOUS :: totuflux(:,0:)      ! upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out), CONTIGUOUS :: totdflux(:,0:)      ! downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out), CONTIGUOUS :: fnet(:,0:)          ! net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out), CONTIGUOUS :: totuclfl(:,0:)      ! clear sky upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out) :: totdclfl(:,0:)      ! clear sky downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out), CONTIGUOUS :: fnetc(:,0:)         ! clear sky net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out), CONTIGUOUS :: dtotuflux_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)
    REAL(wp), INTENT(out), CONTIGUOUS :: dtotuclfl_dt(:,0:)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)

    ! Vectorized version implemented by Guenther Zaengl, DWD
    ! See above for variable descriptions
    ! ----- Local -----
    ! Declarations for radiative transfer
    REAL(wp) :: atot(kproma,nlayers)
    REAL(wp) :: atrans(kproma,nlayers)
    REAL(wp) :: bbugas(kproma,nlayers)
    REAL(wp) :: bbutot(kproma,nlayers)
    REAL(wp) :: clrurad(kproma,0:nlayers)
    REAL(wp) :: clrdrad(kproma,0:nlayers)
    REAL(wp) :: uflux
    REAL(wp) :: dflux
    REAL(wp) :: urad(kproma,0:nlayers)
    REAL(wp) :: drad(kproma,0:nlayers)
    REAL(wp) :: uclfl
    REAL(wp) :: dclfl

    REAL(wp) :: secdiff(kproma,nbndlw)          ! secant of diffusivity angle
    REAL(wp) :: dplankup(kproma,nlayers), dplankdn(kproma,nlayers)
    REAL(wp), PARAMETER :: rec_6 = 0.166667_wp, &
      ! diffusivity angle adjustment coefficients
      a0(nbndlw) = (/ 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
      & 1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
      & 1.668_wp,  1.66_wp,  1.66_wp,  1.66_wp, &
      & 1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /),&
      a1(nbndlw) = (/ 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
      & 0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
      & -0.006_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /),&
      a2(nbndlw) = (/ 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
      & -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
      & 0.414_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
      & 0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /), &
    ! This secant and weight corresponds to the standard diffusivity
    ! angle.  This initial value is redefined below for some bands.
      wtdiff = 0.5_wp
    REAL(wp) :: radld(kproma), radclrd(kproma), plfrac
    REAL(wp) :: odepth, odtot, gassrc, ttot
    REAL(wp) :: bbd(kproma)
    REAL(wp) :: rad0, reflect, radlu(kproma), radclru(kproma)

    REAL(wp) :: duflux_dt
    REAL(wp) :: duclfl_dt
    REAL(wp) :: d_urad_dt(kproma,0:nlayers)
    REAL(wp) :: d_clrurad_dt(kproma,0:nlayers)
    REAL(wp) :: d_rad0_dt, d_radlu_dt(kproma), d_radclru_dt(kproma)

    LOGICAL :: lcldlyr(kproma,0:nlayers+1)             ! flag for cloud in layer
    INTEGER :: ibnd, ib, iband, lay, lev    ! loop indices
    INTEGER :: igc                          ! g-point interval counter
    LOGICAL :: iclddn(kproma)               ! flag for cloud in down path
    INTEGER :: ittot, itgas            ! lookup table indices
    INTEGER, PARAMETER :: ipat(16,0:2) = RESHAPE((/&
      ! These arrays indicate the spectral 'region' (used in the
      ! calculation of ice cloud optical depths) corresponding
      ! to each spectral band.  See cldprop.f for more details.
      & 1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
      & 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
      & 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/), (/16, 3/))

    INTEGER :: jl

    INTEGER :: ibv(kproma)

    ! Declarations for cloud overlap adjustment
    REAL(wp) :: faccld1(kproma,nlayers+1),faccld2(kproma,nlayers+1)
    REAL(wp) :: facclr1(kproma,nlayers+1),facclr2(kproma,nlayers+1)
    REAL(wp) :: faccmb1(kproma,nlayers+1),faccmb2(kproma,nlayers+1)
    REAL(wp) :: faccld1d(kproma,0:nlayers),faccld2d(kproma,0:nlayers)
    REAL(wp) :: facclr1d(kproma,0:nlayers),facclr2d(kproma,0:nlayers)
    REAL(wp) :: faccmb1d(kproma,0:nlayers),faccmb2d(kproma,0:nlayers)
    !--------------------------------------------------------------------------
    ! Maximum/Random cloud overlap variables
    ! for upward radiative transfer
    !  facclr2  fraction of clear radiance from previous layer that needs to
    !           be switched to cloudy stream
    !  facclr1  fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  faccld2  fraction of cloudy radiance from previous layer that needs to
    !           be switched to clear stream
    !  faccld1  fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    ! for downward radiative transfer
    !  facclr2d fraction of clear radiance from previous layer that needs to
    !           be switched to cloudy stream
    !  facclr1d fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  faccld2d fraction of cloudy radiance from previous layer that needs to
    !           be switched to clear stream
    !  faccld1d fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    !--------------------------------------------------------------------------

    REAL(wp) :: cldsrc, oldclr, oldcld, radmod
    REAL(wp) :: odepth_rec, odtot_rec, tblind, tfactot, tfacgas, bbdtot, &
      transc, tausfac
    DIMENSION :: gassrc(kproma), cldsrc(kproma), bbdtot(kproma), &
      oldcld(kproma), oldclr(kproma), odepth(kproma), odtot(kproma), &
      radmod(kproma)
    INTEGER :: itr, icld1, npoints1, npoints2, npoints3, &
      ilist1(kproma), ilist2(kproma), ilist3(kproma)
    REAL(wp), DIMENSION(kproma) :: clrradd, cldradd, clrradu, cldradu, &
      & rad





    ! ------- Definitions -------
    ! input
    !    nlayers                      ! number of model layers
    !    ngptlw                       ! total number of g-point subintervals
    !    nbndlw                       ! number of longwave spectral bands
    !    ncbands                      ! number of spectral bands for clouds
    !    secdiff                      ! diffusivity angle
    !    wtdiff                       ! weight for radiance to flux conversion
    !    pavel                        ! layer pressures (mb)
    !    tavel                        ! layer temperatures (k)
    !    tz                           ! level (interface) temperatures(mb)
    !    tbound                       ! surface temperature (k)
    !    cldfrac                      ! layer cloud fraction
    !    taucloud                     ! layer cloud optical depth
    !    lcldlyr                      ! flag for cloudy layers
    !    iclddn                       ! flag for cloud in column at any layer
    !    semiss                       ! surface emissivities for each band
    !    reflect                      ! surface reflectance
    !    bpade                        ! 1/(pade constant)
    !    exp_tbl                      ! exponential look-up table for transmittance
    !    tfn_tbl                      ! tau transition function look-up table

    ! local
    !    atrans                       ! gaseous absorptivity
    !    atot                         ! combined gaseous and cloud absorptivity
    !    odclr                        ! clear sky (gaseous) optical depth
    !    odcld                        ! cloud optical depth
    !    odtot                        ! optical depth of gas and cloud
    !    tfacgas                      ! gas-only pade factor, used for planck fn
    !    tfactot                      ! gas and cloud pade factor, used for planck fn
    !    bbdgas                       ! gas-only planck function for downward rt
    !    bbugas                       ! gas-only planck function for upward rt
    !    bbdtot                       ! gas and cloud planck function for downward rt
    !    bbutot                       ! gas and cloud planck function for upward calc.
    !    gassrc                       ! source radiance due to gas only
    !    radlu                        ! spectrally summed upward radiance
    !    radclru                      ! spectrally summed clear sky upward radiance
    !    urad                         ! upward radiance by layer
    !    clrurad                      ! clear sky upward radiance by layer
    !    radld                        ! spectrally summed downward radiance
    !    radclrd                      ! spectrally summed clear sky downward radiance
    !    drad                         ! downward radiance by layer
    !    clrdrad                      ! clear sky downward radiance by layer
    !    d_radlu_dt                   ! spectrally summed upward radiance
    !    d_radclru_dt                 ! spectrally summed clear sky upward radiance
    !    d_urad_dt                    ! upward radiance by layer
    !    d_clrurad_dt                 ! clear sky upward radiance by layer

    ! output
    !    totuflux                     ! upward longwave flux (w/m2)
    !    totdflux                     ! downward longwave flux (w/m2)
    !    fnet                         ! net longwave flux (w/m2)
    !    totuclfl                     ! clear sky upward longwave flux (w/m2)
    !    totdclfl                     ! clear sky downward longwave flux (w/m2)
    !    fnetc                        ! clear sky net longwave flux (w/m2)
    !    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
    !                                 ! with respect to surface temperature
    !    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
    !

    ! Local variables for cloud / no cloud index lists
    !    icld_ind(1:n_cloudpoints(lev),lev) stores indices of cloudy points
    INTEGER :: icld_ind(kproma,nlayers)
    INTEGER :: icld, iclear, n_cloudpoints(nlayers)


    ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.

    DO ibnd = 1,nbndlw
      IF (ibnd.EQ.1 .OR. ibnd.EQ.4 .OR. ibnd.GE.10) THEN
!$OMP SIMD
        DO jl = 1, kproma  ! loop over columns
          secdiff(jl,ibnd) = 1.66_wp
        ENDDO
      ELSE
        DO jl = 1, kproma  ! loop over columns
          secdiff(jl,ibnd) = a0(ibnd) + a1(ibnd)*EXP(a2(ibnd)*pwvcm(jl))
          IF (secdiff(jl,ibnd) .GT. 1.80_wp) secdiff(jl,ibnd) = 1.80_wp
          IF (secdiff(jl,ibnd) .LT. 1.50_wp) secdiff(jl,ibnd) = 1.50_wp
        ENDDO
      ENDIF
    ENDDO

    urad(:,:)     = 0.0_wp
    drad(:,:)     = 0.0_wp
    totuflux(:,:) = 0.0_wp
    totdflux(:,:) = 0.0_wp
    clrurad(:,:)  = 0.0_wp
    clrdrad(:,:)  = 0.0_wp
    totuclfl(:,:) = 0.0_wp
    totdclfl(:,:) = 0.0_wp

    IF (idrv .EQ. 1) THEN
      d_urad_dt(:,:)    = 0.0_wp
      d_clrurad_dt(:,:) = 0.0_wp
      dtotuflux_dt(:,:) = 0.0_wp
      dtotuclfl_dt(:,:) = 0.0_wp
    ENDIF

    lcldlyr(:, 0) = .FALSE.
    DO lay = 1, nlayers

      icld   = 0
      iclear = 0

      DO jl = 1, kproma  ! loop over columns
        IF (cldfrac(jl,lay) .GE. 1.e-6_wp) THEN
          lcldlyr(jl,lay) = .TRUE.
          icld = icld + 1
          icld_ind(icld,lay) = jl
        ELSE
          lcldlyr(jl,lay) = .FALSE.
          icld_ind(kproma - iclear,lay) = jl
          iclear = iclear + 1
        ENDIF
      ENDDO

      n_cloudpoints(lay) = icld

    ENDDO
    lcldlyr(:, nlayers+1) = .FALSE.

    ! Maximum/Random cloud overlap parameter

    CALL cloud_overlap(kproma, nlayers, 1, cldfrac, &
      lcldlyr, n_cloudpoints, icld_ind, &
      1, nlayers, 1, &
      faccld1, faccld2, facclr1, facclr2, faccmb1, faccmb2)

    ! istcld(i,j) = .true. if j == 1 or .not. lcldlyr(i, j - 1)
    ! istcldd(i,j) = .true. if j == nlayers or .not. lcldlyr(i, j + 1)

    CALL cloud_overlap(kproma, nlayers, 0, cldfrac, &
      lcldlyr, n_cloudpoints, icld_ind, &
      nlayers, 1, -1, &
      faccld1d, faccld2d, facclr1d, facclr2d, faccmb1d, faccmb2d)

    igc = 1
    ! Loop over frequency bands.
    DO iband = istart, iend

      DO jl = 1, kproma
        IF (ncbands(jl) .EQ. 1) THEN
          ibv(jl) = ipat(iband,0)
        ELSEIF (ncbands(jl) .EQ. 16) THEN
          ibv(jl) = ipat(iband,2)
        ELSEIF (ncbands(jl) .EQ.  5) THEN
          ibv(jl) = ipat(iband,1)
        ENDIF
      ENDDO

      DO lev = 1, nlayers
!DIR$ SIMD
        DO jl = 1, kproma
          dplankup(jl,lev) = planklev(jl,lev,iband) - planklay(jl,lev,iband)
          dplankdn(jl,lev) = planklev(jl,lev-1,iband) - planklay(jl,lev,iband)
        ENDDO
      ENDDO

      ! Reinitialize g-point counter for each band if output for each band is requested.
      IF (iout.GT.0.AND.iband.GE.2) igc = ngs(iband-1)+1

      ! Loop over g-channels.
      1000    CONTINUE

      ! Radiative transfer starts here.
      radld(:) = 0._wp
      radclrd(:) = 0._wp
      iclddn(:) = .FALSE.

      ! Downward radiative transfer loop.
      DO lev = nlayers, 1, -1
        IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy

          DO jl = 1, kproma ! Thus, direct addressing can be used

            ib = ibv(jl)
            plfrac = fracs(jl,lev,igc)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            iclddn(jl) = .TRUE.
            odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)
            IF (odtot(jl) .LT. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                         + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

              atot(jl,lev) = odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)
              odtot_rec = rec_6*odtot(jl)
              bbdtot(jl) =  plfrac * (planklay(jl,lev,iband)+dplankdn(jl,lev)*odtot_rec)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)

              bbugas(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
              bbutot(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odtot_rec)
            ELSEIF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                         + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

              tblind = odtot(jl)/(bpade+odtot(jl))
              ittot = INT(tblint*tblind + 0.5_wp)
              tfactot = tfn_tbl(ittot)
              bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
              atot(jl,lev) = 1._wp - exp_tbl(ittot)

              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + dplankup(jl,lev)*odepth_rec)
              bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itgas = INT(tblint*tblind+0.5_wp)
              odepth(jl) = tau_tbl(itgas)
              atrans(jl,lev) = 1._wp - exp_tbl(itgas)
              tfacgas = tfn_tbl(itgas)
              gassrc(jl) = plfrac * (planklay(jl,lev,iband) &
                + tfacgas*dplankdn(jl,lev)) * atrans(jl,lev)

              tblind = odtot(jl)/(bpade+odtot(jl))
              ittot = INT(tblint*tblind + 0.5_wp)
              tfactot = tfn_tbl(ittot)
              bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+tfacgas*dplankdn(jl,lev))
              atot(jl,lev) = 1._wp - exp_tbl(ittot)

              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfacgas * dplankup(jl,lev))
              bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
            ENDIF

            IF (.NOT. lcldlyr(jl,lev+1)) THEN
              cldradd(jl) = cldfrac(jl,lev) * radld(jl)
              clrradd(jl) = radld(jl) - cldradd(jl)
              oldcld(jl) = cldradd(jl)
              oldclr(jl) = clrradd(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbdtot(jl) * atot(jl,lev)
            cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + &
              & (1._wp-cldfrac(jl,lev))*gassrc(jl)
            radld(jl) = cldradd(jl) + clrradd(jl)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

            radmod(jl) = rad(jl) * &
              & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
              & faccld1d(jl,lev-1) *  ttot) - &
              & faccmb1d(jl,lev-1) * gassrc(jl) + &
              & faccmb2d(jl,lev-1) * cldsrc(jl)

            oldcld(jl) = cldradd(jl) - radmod(jl)
            oldclr(jl) = clrradd(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2d(jl,lev-1)*oldclr(jl) -&
              &  faccld2d(jl,lev-1)*oldcld(jl)
            cldradd(jl) = cldradd(jl) + rad(jl)
            clrradd(jl) = clrradd(jl) - rad(jl)
          ENDDO

        ELSE IF (n_cloudpoints(lev) == 0) THEN ! all points are clear

!DIR$ IVDEP
          DO jl = 1, kproma ! Thus, direct addressing can be used

            plfrac = fracs(jl,lev,igc)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            IF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
              bbugas(jl,lev) = plfrac*(planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itr = INT(tblint*tblind+0.5_wp)
              transc = exp_tbl(itr)
              atrans(jl,lev) = 1._wp-transc
              tausfac = tfn_tbl(itr)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+tausfac*dplankdn(jl,lev))
              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tausfac * dplankup(jl,lev))
            ENDIF
            radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ELSE ! both cloudy and clear points are in the vector

          npoints1 = 0
          npoints2 = 0
          npoints3 = 0

          ! Cloudy layer
!$NEC ivdep
!DIR$ IVDEP
          DO icld = 1, n_cloudpoints(lev)
            jl = icld_ind(icld,lev)

            ib = ibv(jl)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            iclddn(jl) = .TRUE.
            odtot(jl) = odepth(jl) + secdiff(jl,ib) * taucloud(jl,lev,ib)

            IF (odtot(jl) .LT. 0.06_wp) THEN
              npoints1 = npoints1 + 1
              ilist1(npoints1) = jl
            ELSEIF (odepth(jl) .LE. 0.06_wp) THEN
              npoints2 = npoints2 + 1
              ilist2(npoints2) = jl
            ELSE
              npoints3 = npoints3 + 1
              ilist3(npoints3) = jl
            ENDIF
          ENDDO

!$NEC ivdep
          DO icld1 = 1, npoints1
            jl = ilist1(icld1)

            plfrac = fracs(jl,lev,igc)
            atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
            odepth_rec = rec_6*odepth(jl)
            gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                       + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

            atot(jl,lev) = odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)
            odtot_rec = rec_6*odtot(jl)
            bbdtot(jl) =  plfrac * (planklay(jl,lev,iband)+dplankdn(jl,lev)*odtot_rec)
            bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)

            bbugas(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
            bbutot(jl,lev) =  plfrac * (planklay(jl,lev,iband)+dplankup(jl,lev)*odtot_rec)
          ENDDO

!$NEC ivdep
!DIR$ IVDEP
          DO icld1 = 1, npoints2
            jl = ilist2(icld1)

            plfrac = fracs(jl,lev,igc)
            atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
            odepth_rec = rec_6*odepth(jl)
            gassrc(jl) = plfrac*(planklay(jl,lev,iband) &
                       + dplankdn(jl,lev)*odepth_rec)*atrans(jl,lev)

            tblind = odtot(jl)/(bpade+odtot(jl))
            ittot = INT(tblint*tblind + 0.5_wp)
            tfactot = tfn_tbl(ittot)
            bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
            bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
            atot(jl,lev) = 1._wp - exp_tbl(ittot)

            bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + dplankup(jl,lev)*odepth_rec)
            bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
          ENDDO

!$NEC ivdep
!DIR$ IVDEP
          DO icld1 = 1, npoints3
            jl = ilist3(icld1)

            plfrac = fracs(jl,lev,igc)
            tblind = odepth(jl)/(bpade+odepth(jl))
            itgas = INT(tblint*tblind+0.5_wp)
            odepth(jl) = tau_tbl(itgas)
            atrans(jl,lev) = 1._wp - exp_tbl(itgas)
            tfacgas = tfn_tbl(itgas)
            gassrc(jl) = plfrac * (planklay(jl,lev,iband) &
              + tfacgas*dplankdn(jl,lev)) * atrans(jl,lev)

            tblind = odtot(jl)/(bpade+odtot(jl))
            ittot = INT(tblint*tblind + 0.5_wp)
            tfactot = tfn_tbl(ittot)
            bbdtot(jl) = plfrac * (planklay(jl,lev,iband) + tfactot*dplankdn(jl,lev))
            bbd(jl) = plfrac*(planklay(jl,lev,iband)+tfacgas*dplankdn(jl,lev))
            atot(jl,lev) = 1._wp - exp_tbl(ittot)

            bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfacgas * dplankup(jl,lev))
            bbutot(jl,lev) = plfrac * (planklay(jl,lev,iband) + tfactot * dplankup(jl,lev))
          ENDDO

!$NEC ivdep
!DIR$ IVDEP
          DO icld = 1, n_cloudpoints(lev)
            jl = icld_ind(icld,lev)

            IF (.NOT. lcldlyr(jl, lev + 1)) THEN
              cldradd(jl) = cldfrac(jl,lev) * radld(jl)
              clrradd(jl) = radld(jl) - cldradd(jl)
              oldcld(jl) = cldradd(jl)
              oldclr(jl) = clrradd(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbdtot(jl) * atot(jl,lev)
            cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + &
              & (1._wp-cldfrac(jl,lev))*gassrc(jl)
            radld(jl) = cldradd(jl) + clrradd(jl)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

            radmod(jl) = rad(jl) * &
              & (facclr1d(jl,lev-1) * (1._wp-atrans(jl,lev)) + &
              & faccld1d(jl,lev-1) *  ttot) - &
              & faccmb1d(jl,lev-1) * gassrc(jl) + &
              & faccmb2d(jl,lev-1) * cldsrc(jl)

            oldcld(jl) = cldradd(jl) - radmod(jl)
            oldclr(jl) = clrradd(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2d(jl,lev-1)*oldclr(jl) -&
              &  faccld2d(jl,lev-1)*oldcld(jl)
            cldradd(jl) = cldradd(jl) + rad(jl)
            clrradd(jl) = clrradd(jl) - rad(jl)
          ENDDO

          ! Clear layer
!$NEC ivdep
!DIR$ IVDEP
          DO iclear = 0, kproma - n_cloudpoints(lev) - 1
            jl = icld_ind(kproma - iclear,lev)

            plfrac = fracs(jl,lev,igc)
            odepth(jl) = MAX(0.0_wp, secdiff(jl,iband) * taut(jl,lev,igc))

            IF (odepth(jl) .LE. 0.06_wp) THEN
              atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)
              odepth_rec = rec_6*odepth(jl)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+dplankdn(jl,lev)*odepth_rec)
              bbugas(jl,lev) = plfrac*(planklay(jl,lev,iband)+dplankup(jl,lev)*odepth_rec)
            ELSE
              tblind = odepth(jl)/(bpade+odepth(jl))
              itr = INT(tblint*tblind+0.5_wp)
              transc = exp_tbl(itr)
              atrans(jl,lev) = 1._wp-transc
              tausfac = tfn_tbl(itr)
              bbd(jl) = plfrac*(planklay(jl,lev,iband)+tausfac*dplankdn(jl,lev))
              bbugas(jl,lev) = plfrac * (planklay(jl,lev,iband) + tausfac * dplankup(jl,lev))
            ENDIF
            radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
            drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
          ENDDO

        ENDIF

        !  Set clear sky stream to total sky stream as long as layers
        !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
        !  and clear sky stream must be computed separately from that point.
!DIR$ IVDEP
        DO jl = 1, kproma
          IF (iclddn(jl)) THEN
            radclrd(jl) = radclrd(jl) + (bbd(jl)-radclrd(jl)) * atrans(jl,lev)
            clrdrad(jl,lev-1) = clrdrad(jl,lev-1) + radclrd(jl)
          ELSE
            radclrd(jl) = radld(jl)
            clrdrad(jl,lev-1) = drad(jl,lev-1)
          ENDIF
        ENDDO
      ENDDO

      DO jl = 1, kproma
        ! Spectral emissivity & reflectance
        !  Include the contribution of spectrally varying longwave emissivity
        !  and reflection from the surface to the upward radiative transfer.
        !  Note: Spectral and Lambertian reflection are identical for the
        !  diffusivity angle flux integration used here.
        !  Note: The emissivity is applied to plankbnd and dplankbnd_dt when
        !  they are defined in subroutine setcoef.

        rad0 = fracs(jl,1,igc) * plankbnd(jl,iband)

        !  Add in reflection of surface downward radiance.
        reflect = 1._wp - semiss(jl,iband)
        radlu(jl) = rad0 + reflect * radld(jl)
        radclru(jl) = rad0 + reflect * radclrd(jl)

        ! Upward radiative transfer loop.

        urad(jl,0) = urad(jl,0) + radlu(jl)
        clrurad(jl,0) = clrurad(jl,0) + radclru(jl)
      ENDDO
      IF (idrv .EQ. 1) THEN
!DIR$ SIMD
        DO jl = 1, kproma
          d_rad0_dt = fracs(jl,1,igc) * dplankbnd_dt(jl,iband)
          d_radlu_dt(jl) = d_rad0_dt
          d_urad_dt(jl,0) = d_urad_dt(jl,0) + d_radlu_dt(jl)
          d_radclru_dt(jl) = d_rad0_dt
          d_clrurad_dt(jl,0) = d_clrurad_dt(jl,0) + d_radclru_dt(jl)
        ENDDO
      ENDIF

      DO lev = 1, nlayers
        IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy

          DO jl = 1, kproma ! Thus, direct addressing can be used

            gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
            IF (.NOT. lcldlyr(jl,lev-1)) THEN
              cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
              clrradu(jl) = radlu(jl) - cldradu(jl)
              oldcld(jl) = cldradu(jl)
              oldclr(jl) = clrradu(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbutot(jl,lev) * atot(jl,lev)
            cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev))+(1._wp-cldfrac(jl,lev))*gassrc(jl)
            ! Total sky radiance
            radlu(jl) = cldradu(jl) + clrradu(jl)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            radmod(jl) = rad(jl) * &
              & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
              & faccld1(jl,lev+1) *  ttot) - &
              & faccmb1(jl,lev+1) * gassrc(jl) + &
              & faccmb2(jl,lev+1) * cldsrc(jl)
            oldcld(jl) = cldradu(jl) - radmod(jl)
            oldclr(jl) = clrradu(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2(jl,lev+1)*oldclr(jl) - faccld2(jl,lev+1)*oldcld(jl)
            cldradu(jl) = cldradu(jl) + rad(jl)
            clrradu(jl) = clrradu(jl) - rad(jl)
          ENDDO
          IF (idrv .EQ. 1) THEN
            DO jl = 1, kproma
              d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                & d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
              d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            ENDDO
          ENDIF

        ELSE IF (n_cloudpoints(lev) == 0) THEN ! all points are clear

!DIR$ IVDEP
          DO jl = 1, kproma ! thus, direct addressing can be used
            radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
          ENDDO
          IF (idrv .EQ. 1) THEN
            DO jl = 1, kproma
              d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
              d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            ENDDO
          ENDIF

        ELSE ! both cloudy and clear points are in the vector

          ! Cloudy layer
!$NEC ivdep
!DIR$ IVDEP
          DO icld = 1, n_cloudpoints(lev)
            jl = icld_ind(icld,lev)

            gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
            IF (.NOT. lcldlyr(jl, lev - 1)) THEN
              cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
              clrradu(jl) = radlu(jl) - cldradu(jl)
              oldcld(jl) = cldradu(jl)
              oldclr(jl) = clrradu(jl)
              rad(jl) = 0._wp
            ENDIF
            ttot = 1._wp - atot(jl,lev)
            cldsrc(jl) = bbutot(jl,lev) * atot(jl,lev)
            cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc(jl)
            clrradu(jl) = clrradu(jl) * (1.0_wp-atrans(jl,lev))+(1._wp-cldfrac(jl,lev))*gassrc(jl)
            ! Total sky radiance
            radlu(jl) = cldradu(jl) + clrradu(jl)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
            radmod(jl) = rad(jl) * &
              & (facclr1(jl,lev+1)*(1.0_wp-atrans(jl,lev))+ &
              & faccld1(jl,lev+1) *  ttot) - &
              & faccmb1(jl,lev+1) * gassrc(jl) + &
              & faccmb2(jl,lev+1) * cldsrc(jl)
            oldcld(jl) = cldradu(jl) - radmod(jl)
            oldclr(jl) = clrradu(jl) + radmod(jl)
            rad(jl) = -radmod(jl) + facclr2(jl,lev+1)*oldclr(jl) - faccld2(jl,lev+1)*oldcld(jl)
            cldradu(jl) = cldradu(jl) + rad(jl)
            clrradu(jl) = clrradu(jl) - rad(jl)
          ENDDO

          ! Clear layer
!$NEC ivdep
!DIR$ IVDEP
          DO iclear = 0, kproma - n_cloudpoints(lev) - 1
            jl = icld_ind(kproma - iclear,lev)

            radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
            urad(jl,lev) = urad(jl,lev) + radlu(jl)
          ENDDO

          IF (idrv .EQ. 1) THEN
            DO jl = 1, kproma
              IF (lcldlyr(jl,lev)) THEN
                d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                  & d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
                d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              ELSE
                d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
                d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              ENDIF
            ENDDO
          ENDIF

        ENDIF

        !  Set clear sky stream to total sky stream as long as all layers
        !  are clear (iclddn=true).  Streams must be calculated separately at
        !  all layers when a cloud is present (iclddn=false), because surface
        !  reflectance is different for each stream.
!DIR$ SIMD
        DO jl = 1, kproma
          radclru(jl) = MERGE(radclru(jl) + (bbugas(jl,lev)-radclru(jl))*atrans(jl,lev), &
               radlu(jl), iclddn(jl))
          clrurad(jl,lev) = MERGE(clrurad(jl,lev) + radclru(jl), urad(jl,lev), &
               iclddn(jl))
        ENDDO

        IF (idrv .EQ. 1) THEN
!DIR$ SIMD
          DO jl = 1, kproma
            d_radlu_dt(jl) = d_radlu_dt(jl) * MERGE(cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                 & (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev)), &
                 & (1.0_wp - atrans(jl,lev)), lcldlyr(jl,lev))
            d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
            d_radclru_dt(jl) = MERGE(d_radclru_dt(jl) * (1.0_wp - atrans(jl,lev)), &
                 d_radlu_dt(jl), iclddn(jl))
            d_clrurad_dt(jl,lev) = MERGE(d_clrurad_dt(jl,lev) + d_radclru_dt(jl), &
                 d_urad_dt(jl,lev), iclddn(jl))
          ENDDO
        ENDIF
      ENDDO

      ! Increment g-point counter
      igc = igc + 1
      ! Return to continue radiative transfer for all g-channels in present band
      IF (igc .LE. ngs(iband)) GO TO 1000

      ! Process longwave output from band.
      ! Calculate upward, downward, and net flux.
      DO lev = nlayers, 0, -1
!DIR$ SIMD
        DO jl = 1, kproma
          uflux = urad(jl,lev)*wtdiff
          dflux = drad(jl,lev)*wtdiff
          urad(jl,lev) = 0.0_wp
          drad(jl,lev) = 0.0_wp
          totuflux(jl,lev) = totuflux(jl,lev) + uflux * delwave(iband)
          totdflux(jl,lev) = totdflux(jl,lev) + dflux * delwave(iband)
          uclfl = clrurad(jl,lev)*wtdiff
          dclfl = clrdrad(jl,lev)*wtdiff
          clrurad(jl,lev) = 0.0_wp
          clrdrad(jl,lev) = 0.0_wp
          totuclfl(jl,lev) = totuclfl(jl,lev) + uclfl * delwave(iband)
          totdclfl(jl,lev) = totdclfl(jl,lev) + dclfl * delwave(iband)
        ENDDO
      ENDDO

      ! Calculate total change in upward flux wrt surface temperature
      IF (idrv .EQ. 1) THEN
        DO lev = nlayers, 0, -1
!DIR$ SIMD
          DO jl = 1, kproma
            duflux_dt = d_urad_dt(jl,lev) * wtdiff
            d_urad_dt(jl,lev) = 0.0_wp
            dtotuflux_dt(jl,lev) = dtotuflux_dt(jl,lev) + duflux_dt * delwave(iband) * fluxfac

            duclfl_dt = d_clrurad_dt(jl,lev) * wtdiff
            d_clrurad_dt(jl,lev) = 0.0_wp
            dtotuclfl_dt(jl,lev) = dtotuclfl_dt(jl,lev) + duclfl_dt * delwave(iband) * fluxfac
          ENDDO
        ENDDO
      ENDIF

      ! End spectral band loop
    ENDDO

    ! Calculate fluxes at surface (lev==0) and model levels
!
    DO lev = 0, nlayers
!DIR$ SIMD
      DO jl = 1, kproma  ! loop over columns
        totuflux(jl,lev) = totuflux(jl,lev) * fluxfac
        totdflux(jl,lev) = totdflux(jl,lev) * fluxfac
        fnet(jl,lev) = totuflux(jl,lev) - totdflux(jl,lev)
        totuclfl(jl,lev) = totuclfl(jl,lev) * fluxfac
        totdclfl(jl,lev) = totdclfl(jl,lev) * fluxfac
        fnetc(jl,lev) = totuclfl(jl,lev) - totdclfl(jl,lev)
      ENDDO
    ENDDO

    ! end vectorized version
  END SUBROUTINE lrtm_rtrnmr

  SUBROUTINE cloud_overlap(kproma, nlayers, ofs, cldfrac, &
       lcldlyr, n_cloudpoints, icld_ind, &
       start_lev, end_lev, lev_incr, &
       faccld1, faccld2, facclr1, facclr2, faccmb1, faccmb2)
    INTEGER, INTENT(in) :: kproma          ! number of columns
    INTEGER, INTENT(in) :: nlayers         ! total number of layers
    INTEGER, INTENT(in) :: ofs             ! layer offset to use for
                                           ! references to facc*
    INTEGER, INTENT(in) :: n_cloudpoints(nlayers)
    INTEGER, INTENT(in) :: icld_ind(kproma,nlayers)
    INTEGER, INTENT(in) :: start_lev, end_lev, lev_incr
    LOGICAL, INTENT(in) :: lcldlyr(kproma,0:nlayers+1)

    REAL(wp), INTENT(inout) :: faccld1(kproma,nlayers+1), &
         faccld2(kproma,nlayers+1), &
         facclr1(kproma,nlayers+1), &
         facclr2(kproma,nlayers+1), &
         faccmb1(kproma,nlayers+1), &
         faccmb2(kproma,nlayers+1)
    REAL(wp), INTENT(in), CONTIGUOUS :: cldfrac(:,:)       ! layer cloud fraction

    REAL(wp) :: fmax, fmin, rat1(kproma), rat2(kproma)
    INTEGER :: lev, jl, blev, clev, olev, icld, iclr

    faccld1(:,start_lev-lev_incr+ofs) = 0.0_wp
    faccld2(:,start_lev-lev_incr+ofs) = 0.0_wp
    facclr1(:,start_lev-lev_incr+ofs) = 0.0_wp
    facclr2(:,start_lev-lev_incr+ofs) = 0.0_wp
    faccmb1(:,start_lev-lev_incr+ofs) = 0.0_wp
    faccmb2(:,start_lev-lev_incr+ofs) = 0.0_wp

    DO lev = start_lev, end_lev - lev_incr, lev_incr
      olev = lev + ofs
      blev = lev + 1 - ofs
      clev = lev + lev_incr
      IF (n_cloudpoints(lev) == kproma) THEN ! all points are cloudy
        DO jl = 1, kproma ! Thus, direct addressing can be used
          ! Maximum/random cloud overlap
          IF (cldfrac(jl,clev) .GE. cldfrac(jl,lev)) THEN
            faccld1(jl,olev) = 0._wp
            faccld2(jl,olev) = 0._wp
            IF (.NOT. lcldlyr(jl,lev-lev_incr)) THEN
              facclr1(jl,olev) = 0._wp
              facclr2(jl,olev) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2(jl,olev) = &
                & (cldfrac(jl,clev)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2(jl,blev) = 0._wp
              faccld2(jl,blev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,clev) .GT. fmax) THEN
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = (cldfrac(jl,clev)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,clev) .LT. fmax) THEN
                facclr1(jl,olev) = (cldfrac(jl,clev)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev-lev_incr)-cldfrac(jl,lev))
                facclr2(jl,olev) = 0._wp
              ELSE
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat1(jl) = MERGE(1._wp, 0._wp, &
                 facclr1(jl,olev) > 0._wp .OR. facclr2(jl,olev) > 0._wp)
            rat2(jl) = 0._wp
          ELSE
            facclr1(jl,olev) = 0._wp
            facclr2(jl,olev) = 0._wp
            IF (.NOT. lcldlyr(jl,lev-lev_incr)) THEN
              faccld1(jl,olev) = 0._wp
              faccld2(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,clev))/cldfrac(jl,lev)

              facclr2(jl,blev) = 0._wp
              faccld2(jl,blev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,clev) .LE. fmin) THEN
                faccld1(jl,olev) = rat1(jl)
                faccld2(jl,olev) = (fmin-cldfrac(jl,clev))/fmin
              ELSE
                faccld1(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,clev))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat2(jl) = MERGE(1._wp, 0._wp, &
                 faccld1(jl,olev) > 0._wp .OR. faccld2(jl,olev) > 0._wp)
            rat1(jl) = 0._wp
          ENDIF
          IF (lev == start_lev) THEN
            faccmb1(jl,olev) = 0._wp
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,blev)
          ELSE
            faccmb1(jl,olev) = facclr1(jl,olev) * faccld2(jl,blev) * cldfrac(jl,lev-lev_incr)
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,blev) * (1._wp - cldfrac(jl,lev-lev_incr))
          ENDIF
        ENDDO
      ELSE IF (n_cloudpoints(lev) /= 0) THEN
        ! use index list for the case that not all points are cloudy
!$NEC ivdep
        DO icld = 1, n_cloudpoints(lev)
          jl = icld_ind(icld,lev)
          ! Maximum/random cloud overlap
          IF (cldfrac(jl,clev) .GE. cldfrac(jl,lev)) THEN
            faccld1(jl,olev) = 0._wp
            faccld2(jl,olev) = 0._wp
            IF (.NOT. lcldlyr(jl,lev-lev_incr)) THEN
              facclr1(jl,olev) = 0._wp
              facclr2(jl,olev) = 0._wp
              IF (cldfrac(jl,lev) .LT. 1._wp) facclr2(jl,olev) = &
                & (cldfrac(jl,clev)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
              facclr2(jl,blev) = 0._wp
              faccld2(jl,blev) = 0._wp
            ELSE
              fmax = MAX(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,clev) .GT. fmax) THEN
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = (cldfrac(jl,clev)-fmax)/(1._wp-fmax)
              ELSEIF (cldfrac(jl,clev) .LT. fmax) THEN
                facclr1(jl,olev) = (cldfrac(jl,clev)-cldfrac(jl,lev))/ &
                  & (cldfrac(jl,lev-lev_incr)-cldfrac(jl,lev))
                facclr2(jl,olev) = 0._wp
              ELSE
                facclr1(jl,olev) = rat2(jl)
                facclr2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat1(jl) = MERGE(1._wp, 0._wp, &
              facclr1(jl,olev) > 0._wp .OR. facclr2(jl,olev) > 0._wp)
            rat2(jl) = 0._wp
          ELSE
            facclr1(jl,olev) = 0._wp
            facclr2(jl,olev) = 0._wp
            IF (.NOT. lcldlyr(jl,lev-lev_incr)) THEN
              faccld1(jl,olev) = 0._wp
              faccld2(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,clev))/cldfrac(jl,lev)

              facclr2(jl,blev) = 0._wp
              faccld2(jl,blev) = 0._wp
            ELSE
              fmin = MIN(cldfrac(jl,lev),cldfrac(jl,lev-lev_incr))
              IF (cldfrac(jl,clev) .LE. fmin) THEN
                faccld1(jl,olev) = rat1(jl)
                faccld2(jl,olev) = (fmin-cldfrac(jl,clev))/fmin
              ELSE
                faccld1(jl,olev) = (cldfrac(jl,lev)-cldfrac(jl,clev))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,olev) = 0._wp
              ENDIF
            ENDIF
            rat2(jl) = MERGE(1._wp, 0._wp, &
              faccld1(jl,olev) > 0._wp .OR. faccld2(jl,olev) > 0._wp)
            rat1(jl) = 0._wp
          ENDIF
          IF (lev == start_lev) THEN
            faccmb1(jl,olev) = 0._wp
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,blev)
          ELSE
            faccmb1(jl,olev) = facclr1(jl,olev) * faccld2(jl,blev) * cldfrac(jl,lev-lev_incr)
            faccmb2(jl,olev) = faccld1(jl,olev) * facclr2(jl,blev) * (1._wp - cldfrac(jl,lev-lev_incr))
          ENDIF
        ENDDO
        DO iclr = n_cloudpoints(lev) + 1, kproma
          jl = icld_ind(iclr, lev)
          faccld1(jl,olev) = 0.0
          faccld2(jl,olev) = 0.0
          facclr1(jl,olev) = 0.0
          facclr2(jl,olev) = 0.0
          faccmb1(jl,olev) = 0.0
          faccmb2(jl,olev) = 0.0
        END DO
      ELSE
        faccld1(:,olev) = 0.0
        faccld2(:,olev) = 0.0
        facclr1(:,olev) = 0.0
        facclr2(:,olev) = 0.0
        faccmb1(:,olev) = 0.0
        faccmb2(:,olev) = 0.0
      ENDIF
    ENDDO
    olev = end_lev + ofs
    faccld1(:,olev) = 0._wp
    faccld2(:,olev) = 0._wp
    facclr1(:,olev) = 0._wp
    facclr2(:,olev) = 0._wp
    faccmb1(:,olev) = 0._wp
    faccmb2(:,olev) = 0._wp

  END SUBROUTINE cloud_overlap

END MODULE mo_lrtm_rtrnmr
!
! Local Variables:
! f90-continuation-indent: 2
! End:
!
