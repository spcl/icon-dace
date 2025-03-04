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


module mo_lrtm_setup

  ! ------- Modules -------
  use mo_kind,        only : wp
  use mo_lrtm_par
  use mo_lrtm_coeffs, only : lwatmref, lwavplank, lwavplankderiv

  implicit none

  private


  public  :: lrtm_setup, ntbl, bpade, tau_tbl, exp_tbl, tfn_tbl

  INTEGER,  PARAMETER :: ntbl = 10000
  REAL(wp), DIMENSION(0:ntbl) :: tau_tbl
  REAL(wp), DIMENSION(0:ntbl) :: exp_tbl
  REAL(wp), DIMENSION(0:ntbl) :: tfn_tbl
  REAL(wp), PARAMETER :: pade   = 0.278_wp   ! Smallest value for exponential table
  REAL(wp), PARAMETER :: bpade = 1.0_wp / pade

contains

  ! **************************************************************************
  subroutine lrtm_setup(data_filename)
    ! **************************************************************************
    !
    !  Original version:       Michael J. Iacono; July, 1998
    !  First revision for GCMs:   September, 1998
    !  Second revision for RRTM_V3.0:  September, 2002
    !
    !  This subroutine performs calculations necessary for the initialization
    !  of the longwave model.  Lookup tables are computed for use in the LW
    !  radiative transfer, and input absorption coefficient data for each
    !  spectral band are reduced from 256 g-point intervals to 140.
    ! **************************************************************************

    use mo_lrtm_netcdf, only: lrtm_read

    !> NetCDF file containing longwave absorption coefficients and other data
    !> for RRTMG_LW k-distribution model
    CHARACTER (LEN=*), INTENT(IN) :: data_filename

    ! ------- Local -------

    integer :: itr, ibnd, igc, ig, ind, ipr
    integer :: igcsm, iprsm

    real(wp) :: wtsum, wtsm(mg)        !
    real(wp) :: tfn                    !

    real(wp), parameter :: expeps = 1.e-20_wp   ! Smallest value for exponential table
    real(wp), parameter :: maxexp = ABS(LOG(expeps))

    ! GZ, 2013-12-04: Turn off vectorization and inlining for the Cray compiler. It generates incorrect code otherwise.




    ! ------- Definitions -------
    !     Arrays for 10000-point look-up tables:
    !     TAU_TBL Clear-sky optical depth (used in cloudy radiative transfer)
    !     EXP_TBL Exponential lookup table for ransmittance
    !     TFN_TBL Tau transition function; i.e. the transition of the Planck
    !             function from that for the mean layer temperature to that for
    !             the layer boundary temperature as a function of optical depth.
    !             The "linear in tau" method is used to make the table.
    !     PADE    Pade approximation constant (= 0.278)
    !     BPADE   Inverse of the Pade approximation constant
    !

    ! Initialize model data
    call lwdatinit
    call lwcmbdat                  ! g-point interval reduction data
    call lwatmref                  ! reference MLS profile
    call lwavplank                 ! Planck function
    call lwavplankderiv            ! Planck function derivative wrt temp
    call lrtm_read(data_filename)  ! molecular absorption coefficients

    ! Compute lookup tables for transmittance, tau transition function,
    ! and clear sky tau (for the cloudy sky radiative transfer).  Tau is
    ! computed as a function of the tau transition function, transmittance
    ! is calculated as a function of tau, and the tau transition function
    ! is calculated using the linear in tau formulation at values of tau
    ! above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables
    ! are computed at intervals of 0.001.  The inverse of the constant used
    ! in the Pade approximation to the tau transition function is set to b.

    tau_tbl(0) = 0.0_wp
    tau_tbl(ntbl) = 1.e10_wp
    exp_tbl(0) = 1.0_wp
    exp_tbl(ntbl) = expeps
    tfn_tbl(0) = 0.0_wp
    tfn_tbl(ntbl) = 1.0_wp

    do itr = 1, ntbl-1
      tfn = REAL(itr,wp) / REAL(ntbl,wp)
      tau_tbl(itr) = bpade * tfn / (1._wp - tfn)
      if (tau_tbl(itr) > maxexp) then
        exp_tbl(itr) = expeps
      else
        exp_tbl(itr) = exp(-tau_tbl(itr))
      end if
      if (tau_tbl(itr) .lt. 0.06_wp) then
        tfn_tbl(itr) = tau_tbl(itr)/6._wp
      else
        tfn_tbl(itr) = 1._wp-2._wp*((1._wp/tau_tbl(itr))-(exp_tbl(itr)/(1._wp-exp_tbl(itr))))
      endif
    enddo

    ! Perform g-point reduction from 16 per band (256 total points) to
    ! a band dependant number (140 total points) for all absorption
    ! coefficient input data and Planck fraction input data.
    ! Compute relative weighting for new g-point combinations.

    igcsm = 0
    do ibnd = 1,nbndlw
      iprsm = 0
      if (ngc(ibnd).lt.mg) then
        do igc = 1,ngc(ibnd)
          igcsm = igcsm + 1
          wtsum = 0._wp
          do ipr = 1, ngn(igcsm)
            iprsm = iprsm + 1
            wtsum = wtsum + wt(iprsm)
          enddo
          wtsm(igc) = wtsum
        enddo
        do ig = 1, ng(ibnd)
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = wt(ig)/wtsm(ngm(ind))
        enddo
      else
        do ig = 1, ng(ibnd)
          igcsm = igcsm + 1
          ind = (ibnd-1)*mg + ig
          rwgt(ind) = 1.0_wp
        enddo
      endif
    enddo

    ! Reduce g-points for absorption coefficient data in each LW spectral band.

    call cmbgb1
    call cmbgb2
    call cmbgb3
    call cmbgb4
    call cmbgb5
    call cmbgb6
    call cmbgb7
    call cmbgb8
    call cmbgb9
    call cmbgb10
    call cmbgb11
    call cmbgb12
    call cmbgb13
    call cmbgb14
    call cmbgb15
    call cmbgb16

  end subroutine lrtm_setup

  !***************************************************************************
  subroutine lwdatinit
    !***************************************************************************

    ! --------- Modules ----------


    save

    ! Longwave spectral band limits (wavenumbers)
    wavenum1(:) = (/ 10._wp, 350._wp, 500._wp, 630._wp, 700._wp, 820._wp, &
         980._wp,1080._wp,1180._wp,1390._wp,1480._wp,1800._wp, &
         2080._wp,2250._wp,2380._wp,2600._wp/)
    wavenum2(:) = (/350._wp, 500._wp, 630._wp, 700._wp, 820._wp, 980._wp, &
         1080._wp,1180._wp,1390._wp,1480._wp,1800._wp,2080._wp, &
         2250._wp,2380._wp,2600._wp,3250._wp/)
    delwave(:) =  (/340._wp, 150._wp, 130._wp,  70._wp, 120._wp, 160._wp, &
         100._wp, 100._wp, 210._wp,  90._wp, 320._wp, 280._wp, &
         170._wp, 130._wp, 220._wp, 650._wp/)

    ! Spectral band information
    ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
    nspa(:) = (/1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)
    nspb(:) = (/1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)

    !     nxmol     - number of cross-sections input by user
    !     ixindx(i) - index of cross-section molecule corresponding to Ith
    !                 cross-section specified by user
    !
    nxmol = 4
    ixindx(1) = 1 ! ccl4
    ixindx(2) = 2 ! cfc11
    ixindx(3) = 3 ! cfc12
    ixindx(4) = 4 ! cfc22
    ixindx(5:maxinpx) = 0

  end subroutine lwdatinit

  !***************************************************************************
  subroutine lwcmbdat
    !***************************************************************************

    save

    ! ------- Definitions -------
    !     Arrays for the g-point reduction from 256 to 140 for the 16 LW bands:
    !     This mapping from 256 to 140 points has been carefully selected to
    !     minimize the effect on the resulting fluxes and cooling rates, and
    !     caution should be used if the mapping is modified.  The full 256
    !     g-point set can be restored with ngptlw=256, ngc=16*16, ngn=256*1., etc.
    !     ngptlw  The total number of new g-points
    !     ngc     The number of new g-points in each band
    !     ngs     The cumulative sum of new g-points for each band
    !     ngm     The index of each new g-point relative to the original
    !             16 g-points for each band.
    !     ngn     The number of original g-points that are combined to make
    !             each new g-point in each band.
    !     ngb     The band index for each new g-point.
    !     wt      RRTM weights for 16 g-points.

    ! ------- Data statements -------
    ngc(:) = (/10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)
    ngs(:) = (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)
    ngm(:) = (/1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
         1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
         1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
         1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
         1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
         1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
         1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
         1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
         1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
         1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
         1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
         1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
         1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
         1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
         1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
         1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/)            ! band 16
    ngn(:) = (/1,1,2,2,2,2,2,2,1,1, &                       ! band 1
         1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
         1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
         2,2,2,2,2,2,2,2, &                           ! band 6
         2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
         2,2,2,2,2,2,2,2, &                           ! band 8
         1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
         2,2,2,2,4,4, &                               ! band 10
         1,1,2,2,2,2,3,3, &                           ! band 11
         1,1,1,1,2,2,4,4, &                           ! band 12
         3,3,4,6, &                                   ! band 13
         8,8, &                                       ! band 14
         8,8, &                                       ! band 15
         4,12/)                                       ! band 16
    ngb(:) = (/1,1,1,1,1,1,1,1,1,1, &                       ! band 1
         2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
         3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
         4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
         5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
         6,6,6,6,6,6,6,6, &                           ! band 6
         7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
         8,8,8,8,8,8,8,8, &                           ! band 8
         9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
         10,10,10,10,10,10, &                         ! band 10
         11,11,11,11,11,11,11,11, &                   ! band 11
         12,12,12,12,12,12,12,12, &                   ! band 12
         13,13,13,13, &                               ! band 13
         14,14, &                                     ! band 14
         15,15, &                                     ! band 15
         16,16/)                                      ! band 16
    wt(:) = (/ 0.1527534276_wp, 0.1491729617_wp, 0.1420961469_wp, &
         0.1316886544_wp, 0.1181945205_wp, 0.1019300893_wp, &
         0.0832767040_wp, 0.0626720116_wp, 0.0424925000_wp, &
         0.0046269894_wp, 0.0038279891_wp, 0.0030260086_wp, &
         0.0022199750_wp, 0.0014140010_wp, 0.0005330000_wp, &
         0.0000750000_wp/)

  end subroutine lwcmbdat

  !***************************************************************************
  subroutine cmbgb1
    !***************************************************************************
    !
    !  Original version:    MJIacono; July 1998
    !  Revision for GCMs:   MJIacono; September 1998
    !  Revision for RRTMG:  MJIacono, September 2002
    !  Revision for F90 reformatting:  MJIacono, June 2006
    !
    !  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
    !  data for each band, which are defined for 16 g-points and 16 spectral
    !  bands. The data are combined with appropriate weighting following the
    !  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
    !  in arrays FRACREFA and FRACREFB are combined without weighting.  All
    !  g-point reduced data are put into new arrays for use in RRTM.
    !
    !  band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
    !                       (high key - h2o; high minor - n2)
    !  note: previous versions of rrtm band 1:
    !        10-250 cm-1 (low - h2o; high - h2o)
    !***************************************************************************

    use mo_rrlw_kg01, only: fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, &
         selfref, forref

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumk1, sumk2, sumf1, sumf2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(1)
          sumk = 0.0_wp
          do ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(1)
          sumk = 0.0_wp
          do ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(1)
        sumk = 0.0_wp
        do ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(1)
        sumk = 0.0_wp
        do ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(1)
        sumk1 = 0.0_wp
        sumk2 = 0.0_wp
        do ipr = 1, ngn(igc)
          iprsm = iprsm + 1
          sumk1 = sumk1 + kao_mn2(jt,iprsm)*rwgt(iprsm)
          sumk2 = sumk2 + kbo_mn2(jt,iprsm)*rwgt(iprsm)
        enddo
        ka_mn2(jt,igc) = sumk1
        kb_mn2(jt,igc) = sumk2
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(1)
      sumf1 = 0.0_wp
      sumf2 = 0.0_wp
      do ipr = 1, ngn(igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      enddo
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    enddo

  end subroutine cmbgb1

  !***************************************************************************
  subroutine cmbgb2
    !***************************************************************************
    !
    !     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
    !
    !     note: previous version of rrtm band 2:
    !           250 - 500 cm-1 (low - h2o; high - h2o)
    !***************************************************************************

    use mo_rrlw_kg02, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, selfref, forref

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf1, sumf2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(2)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+16)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(2)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+16)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(2)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+16)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(2)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(1)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(2)
      sumf1 = 0.0_wp
      sumf2 = 0.0_wp
      do ipr = 1, ngn(ngs(1)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      enddo
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    enddo

  end subroutine cmbgb2

  !***************************************************************************
  subroutine cmbgb3
    !***************************************************************************
    !
    !     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
    !                           (high key - h2o,co2; high minor - n2o)
    !
    ! old band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
    !***************************************************************************

    use mo_rrlw_kg03, only: fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, &
         selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(3)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo
    do jn = 1,5
      do jt = 1,5
        do jp = 13,59
          iprsm = 0
          do igc = 1,ngc(3)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(2)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+32)
            enddo
            kb(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jn = 1,9
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(3)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
          enddo
          ka_mn2o(jn,jt,igc) = sumk
        enddo
      enddo
    enddo

    do jn = 1,5
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(3)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(2)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
          enddo
          kb_mn2o(jn,jt,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(3)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(3)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+32)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(3)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

    do jp = 1,5
      iprsm = 0
      do igc = 1,ngc(3)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(2)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefbo(iprsm,jp)
        enddo
        fracrefb(igc,jp) = sumf
      enddo
    enddo

  end subroutine cmbgb3

  !***************************************************************************
  subroutine cmbgb4
    !***************************************************************************
    !
    !     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
    !
    ! old band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
    !***************************************************************************

    use mo_rrlw_kg04, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(4)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo
    do jn = 1,5
      do jt = 1,5
        do jp = 13,59
          iprsm = 0
          do igc = 1,ngc(4)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(3)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+48)
            enddo
            kb(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(4)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(4)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+48)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(4)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

    do jp = 1,5
      iprsm = 0
      do igc = 1,ngc(4)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(3)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefbo(iprsm,jp)
        enddo
        fracrefb(igc,jp) = sumf
      enddo
    enddo

  end subroutine cmbgb4

  !***************************************************************************
  subroutine cmbgb5
    !***************************************************************************
    !
    !     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
    !                           (high key - o3,co2)
    !
    ! old band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
    !***************************************************************************

    use mo_rrlw_kg05, only: fracrefao, fracrefbo, kao, kbo, kao_mo3, ccl4o, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mo3, ccl4, &
         selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(5)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(4)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+64)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo
    do jn = 1,5
      do jt = 1,5
        do jp = 13,59
          iprsm = 0
          do igc = 1,ngc(5)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(4)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+64)
            enddo
            kb(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jn = 1,9
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(5)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mo3(jn,jt,iprsm)*rwgt(iprsm+64)
          enddo
          ka_mo3(jn,jt,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(5)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+64)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(5)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(5)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

    do jp = 1,5
      iprsm = 0
      do igc = 1,ngc(5)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(4)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefbo(iprsm,jp)
        enddo
        fracrefb(igc,jp) = sumf
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(5)
      sumk = 0.0_wp
      do ipr = 1, ngn(ngs(4)+igc)
        iprsm = iprsm + 1
        sumk = sumk + ccl4o(iprsm)*rwgt(iprsm+64)
      enddo
      ccl4(igc) = sumk
    enddo

  end subroutine cmbgb5

  !***************************************************************************
  subroutine cmbgb6
    !***************************************************************************
    !
    !     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
    !                           (high key - nothing; high minor - cfc11, cfc12)
    !
    ! old band 6:  820-980 cm-1 (low - h2o; high - nothing)
    !***************************************************************************

    use mo_rrlw_kg06, only: fracrefao, kao, kao_mco2, cfc11adjo, cfc12o, &
         selfrefo, forrefo, &
         fracrefa, ka, ka_mco2, cfc11adj, cfc12, &
         selfref, forref

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf, sumk1, sumk2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(6)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+80)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(6)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kao_mco2(jt,iprsm)*rwgt(iprsm+80)
        enddo
        ka_mco2(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(6)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+80)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(6)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(5)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(6)
      sumf = 0.0_wp
      sumk1= 0.0_wp
      sumk2= 0.0_wp
      do ipr = 1, ngn(ngs(5)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefao(iprsm)
        sumk1= sumk1+ cfc11adjo(iprsm)*rwgt(iprsm+80)
        sumk2= sumk2+ cfc12o(iprsm)*rwgt(iprsm+80)
      enddo
      fracrefa(igc) = sumf
      cfc11adj(igc) = sumk1
      cfc12(igc) = sumk2
    enddo

  end subroutine cmbgb6

  !***************************************************************************
  subroutine cmbgb7
    !***************************************************************************
    !
    !     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
    !                            (high key - o3; high minor - co2)
    !
    ! old band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
    !***************************************************************************

    use mo_rrlw_kg07, only: fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, &
         selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(7)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(6)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+96)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo
    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(7)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jn = 1,9
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(7)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+96)
          enddo
          ka_mco2(jn,jt,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(7)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kbo_mco2(jt,iprsm)*rwgt(iprsm+96)
        enddo
        kb_mco2(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(7)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(7)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+96)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(7)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(6)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(7)
      sumf = 0.0_wp
      do ipr = 1, ngn(ngs(6)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      enddo
      fracrefb(igc) = sumf
    enddo

  end subroutine cmbgb7

  !***************************************************************************
  subroutine cmbgb8
    !***************************************************************************
    !
    !     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
    !                             (high key - o3; high minor - co2, n2o)
    !
    ! old band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
    !***************************************************************************

    use mo_rrlw_kg08, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
         kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
         cfc12o, cfc22adjo, &
         fracrefa, fracrefb, ka, ka_mco2, ka_mn2o, &
         ka_mo3, kb, kb_mco2, kb_mn2o, selfref, forref, &
         cfc12, cfc22adj

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumk1, sumk2, sumk3, sumk4, sumk5, sumf1, sumf2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(8)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
    enddo
    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(8)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+112)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(8)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(8)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(8)
        sumk1 = 0.0_wp
        sumk2 = 0.0_wp
        sumk3 = 0.0_wp
        sumk4 = 0.0_wp
        sumk5 = 0.0_wp
        do ipr = 1, ngn(ngs(7)+igc)
          iprsm = iprsm + 1
          sumk1 = sumk1 + kao_mco2(jt,iprsm)*rwgt(iprsm+112)
          sumk2 = sumk2 + kbo_mco2(jt,iprsm)*rwgt(iprsm+112)
          sumk3 = sumk3 + kao_mo3(jt,iprsm)*rwgt(iprsm+112)
          sumk4 = sumk4 + kao_mn2o(jt,iprsm)*rwgt(iprsm+112)
          sumk5 = sumk5 + kbo_mn2o(jt,iprsm)*rwgt(iprsm+112)
        enddo
        ka_mco2(jt,igc) = sumk1
        kb_mco2(jt,igc) = sumk2
        ka_mo3(jt,igc) = sumk3
        ka_mn2o(jt,igc) = sumk4
        kb_mn2o(jt,igc) = sumk5
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(8)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      sumk1= 0.0_wp
      sumk2= 0.0_wp
      do ipr = 1, ngn(ngs(7)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
        sumk1= sumk1+ cfc12o(iprsm)*rwgt(iprsm+112)
        sumk2= sumk2+ cfc22adjo(iprsm)*rwgt(iprsm+112)
      enddo
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
      cfc12(igc) = sumk1
      cfc22adj(igc) = sumk2
    enddo

  end subroutine cmbgb8

  !***************************************************************************
  subroutine cmbgb9
    !***************************************************************************
    !
    !     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
    !                             (high key - ch4; high minor - n2o)!

    ! old band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
    !***************************************************************************

    use mo_rrlw_kg09, only: fracrefao, fracrefbo, kao, kao_mn2o, &
         kbo, kbo_mn2o, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, ka_mn2o, &
         kb, kb_mn2o, selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(9)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(8)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+128)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(9)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jn = 1,9
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(9)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+128)
          enddo
          ka_mn2o(jn,jt,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(9)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kbo_mn2o(jt,iprsm)*rwgt(iprsm+128)
        enddo
        kb_mn2o(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(9)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(9)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(9)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(8)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(9)
      sumf = 0.0_wp
      do ipr = 1, ngn(ngs(8)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      enddo
      fracrefb(igc) = sumf
    enddo

  end subroutine cmbgb9

  !***************************************************************************
  subroutine cmbgb10
    !***************************************************************************
    !
    !     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
    !
    ! old band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
    !***************************************************************************

    use mo_rrlw_kg10, only: fracrefao, fracrefbo, kao, kbo, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, &
         selfref, forref

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf1, sumf2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(10)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(10)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+144)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(10)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(9)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+144)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(10)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(9)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+144)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(10)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      do ipr = 1, ngn(ngs(9)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      enddo
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    enddo

  end subroutine cmbgb10

  !***************************************************************************
  subroutine cmbgb11
    !***************************************************************************
    !
    !     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
    !                              (high key - h2o; high minor - o2)
    !
    ! old band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
    !                              (high key - h2o; high minor - o2)
    !***************************************************************************

    use mo_rrlw_kg11, only: fracrefao, fracrefbo, kao, kao_mo2, &
         kbo, kbo_mo2, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, ka_mo2, &
         kb, kb_mo2, selfref, forref

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumk1, sumk2, sumf1, sumf2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(11)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+160)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
    enddo
    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(11)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+160)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(11)
        sumk1 = 0.0_wp
        sumk2 = 0.0_wp
        do ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          sumk1 = sumk1 + kao_mo2(jt,iprsm)*rwgt(iprsm+160)
          sumk2 = sumk2 + kbo_mo2(jt,iprsm)*rwgt(iprsm+160)
        enddo
        ka_mo2(jt,igc) = sumk1
        kb_mo2(jt,igc) = sumk2
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(11)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+160)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(11)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(10)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+160)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(11)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      do ipr = 1, ngn(ngs(10)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      enddo
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    enddo

  end subroutine cmbgb11

  !***************************************************************************
  subroutine cmbgb12
    !***************************************************************************
    !
    !     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    !
    ! old band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
    !***************************************************************************

    use mo_rrlw_kg12, only: fracrefao, kao, selfrefo, forrefo, &
         fracrefa, ka, selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(12)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(11)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+176)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(12)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+176)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(12)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+176)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(12)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(11)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

  end subroutine cmbgb12

  !***************************************************************************
  subroutine cmbgb13
    !***************************************************************************
    !
    !     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
    !
    ! old band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
    !***************************************************************************

    use mo_rrlw_kg13, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
         kbo_mo3, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, ka_mco2, ka_mco, &
         kb_mo3, selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumk1, sumk2, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(13)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(12)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jn = 1,9
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(13)
          sumk1 = 0.0_wp
          sumk2 = 0.0_wp
          do ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            sumk1 = sumk1 + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+192)
            sumk2 = sumk2 + kao_mco(jn,jt,iprsm)*rwgt(iprsm+192)
          enddo
          ka_mco2(jn,jt,igc) = sumk1
          ka_mco(jn,jt,igc) = sumk2
        enddo
      enddo
    enddo

    do jt = 1,19
      iprsm = 0
      do igc = 1,ngc(13)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumk = sumk + kbo_mo3(jt,iprsm)*rwgt(iprsm+192)
        enddo
        kb_mo3(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(13)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+192)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(13)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+192)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(13)
      sumf = 0.0_wp
      do ipr = 1, ngn(ngs(12)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      enddo
      fracrefb(igc) = sumf
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(13)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(12)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

  end subroutine cmbgb13

  !***************************************************************************
  subroutine cmbgb14
    !***************************************************************************
    !
    !     band 14:  2250-2380 cm-1 (low - co2; high - co2)
    !
    ! old band 14:  2250-2380 cm-1 (low - co2; high - co2)
    !***************************************************************************

    use mo_rrlw_kg14, only: fracrefao, fracrefbo, kao, kbo, &
         selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, &
         selfref, forref

    ! ------- Local -------
    integer :: jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf1, sumf2


    do jt = 1,5
      do jp = 1,13
        iprsm = 0
        do igc = 1,ngc(14)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
          enddo
          ka(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(14)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+208)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(14)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+208)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(14)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(13)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(14)
      sumf1= 0.0_wp
      sumf2= 0.0_wp
      do ipr = 1, ngn(ngs(13)+igc)
        iprsm = iprsm + 1
        sumf1= sumf1+ fracrefao(iprsm)
        sumf2= sumf2+ fracrefbo(iprsm)
      enddo
      fracrefa(igc) = sumf1
      fracrefb(igc) = sumf2
    enddo

  end subroutine cmbgb14

  !***************************************************************************
  subroutine cmbgb15
    !***************************************************************************
    !
    !     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
    !                              (high - nothing)
    !
    ! old band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
    !***************************************************************************

    use mo_rrlw_kg15, only: fracrefao, kao, kao_mn2, selfrefo, forrefo, &
         fracrefa, ka, ka_mn2, selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(15)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(14)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+224)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jn = 1,9
      do jt = 1,19
        iprsm = 0
        do igc = 1,ngc(15)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(14)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kao_mn2(jn,jt,iprsm)*rwgt(iprsm+224)
          enddo
          ka_mn2(jn,jt,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(15)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(14)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+224)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(15)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(14)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+224)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(15)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(14)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

  end subroutine cmbgb15

  !***************************************************************************
  subroutine cmbgb16
    !***************************************************************************
    !
    !     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
    !
    ! old band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
    !***************************************************************************

    use mo_rrlw_kg16, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
         fracrefa, fracrefb, ka, kb, selfref, forref

    ! ------- Local -------
    integer :: jn, jt, jp, igc, ipr, iprsm
    real(wp) :: sumk, sumf


    do jn = 1,9
      do jt = 1,5
        do jp = 1,13
          iprsm = 0
          do igc = 1,ngc(16)
            sumk = 0.0_wp
            do ipr = 1, ngn(ngs(15)+igc)
              iprsm = iprsm + 1
              sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+240)
            enddo
            ka(jn,jt,jp,igc) = sumk
          enddo
        enddo
      enddo
    enddo

    do jt = 1,5
      do jp = 13,59
        iprsm = 0
        do igc = 1,ngc(16)
          sumk = 0.0_wp
          do ipr = 1, ngn(ngs(15)+igc)
            iprsm = iprsm + 1
            sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+240)
          enddo
          kb(jt,jp,igc) = sumk
        enddo
      enddo
    enddo

    do jt = 1,10
      iprsm = 0
      do igc = 1,ngc(16)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(15)+igc)
          iprsm = iprsm + 1
          sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+240)
        enddo
        selfref(jt,igc) = sumk
      enddo
    enddo

    do jt = 1,4
      iprsm = 0
      do igc = 1,ngc(16)
        sumk = 0.0_wp
        do ipr = 1, ngn(ngs(15)+igc)
          iprsm = iprsm + 1
          sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+240)
        enddo
        forref(jt,igc) = sumk
      enddo
    enddo

    iprsm = 0
    do igc = 1,ngc(16)
      sumf = 0.0_wp
      do ipr = 1, ngn(ngs(15)+igc)
        iprsm = iprsm + 1
        sumf = sumf + fracrefbo(iprsm)
      enddo
      fracrefb(igc) = sumf
    enddo

    do jp = 1,9
      iprsm = 0
      do igc = 1,ngc(16)
        sumf = 0.0_wp
        do ipr = 1, ngn(ngs(15)+igc)
          iprsm = iprsm + 1
          sumf = sumf + fracrefao(iprsm,jp)
        enddo
        fracrefa(igc,jp) = sumf
      enddo
    enddo

  end subroutine cmbgb16

end module mo_lrtm_setup

