!<<<<<<<<<<<<<<<<fastJX codes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<only access to external variable thru cmn_FJX.f and calls<<<<<<
!<<<<<<<<<<<<<<<<<<version 7.4+  (8/2015, mjp)<<<<<<<<<<<<<<<<<<<<<<<<

!
! !DESCRIPTION: FJX_INIT contains variables and routines to input fast-JX data
!
!
! !INTERFACE:
!
! !USES:
!
!
! !USES:
!
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

      MODULE FJX_INIT_MOD
!
! !USES:
!
      USE FJX_CMN_MOD

      USE FJX_SUB_MOD, ONLY : EXITC
      USE mo_kind,                 ONLY: wp
      

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: INIT_FJX

      CONTAINS

!-----------------------------------------------------------------------
      subroutine INIT_FJX (TITLEJXX,NJXU,NJXX)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  ::NJXU
      integer, intent(out) ::NJXX
      character*6, intent(out), dimension(NJXU) :: TITLEJXX

      integer  JXUNIT,J, RANSEED

      if (W_ .ne. 18) then
        call EXITC(' INIT_JX: invalid no. wavelengths')
      endif
!
!      if(LRRTMG)LRRTMG_INI=.true.
!      if(LRRTMG_INI) call rrtmg_sw_ini(cpdir)

! Use channel 8 to read fastJX data files:
      JXUNIT  = 8

! Read in fast-J X-sections (spectral data)
      !call RD_XXX(JXUNIT,'FJX_spec.dat')
      call RD_XXX(JXUNIT,'FJX_spec_extended_lyman.dat')

! Read in cloud scattering data
      call RD_CLD(JXUNIT,'FJX_scat-cld.dat')

! Read in strat sulf aerosols scattering data
      call RD_SSA(JXUNIT,'FJX_scat-ssa.dat')

! Read in aerosols scattering data
      call RD_MIE(JXUNIT,'FJX_scat-aer.dat')

! Read in UMich aerosol scattering data
      call RD_UM (JXUNIT,'FJX_scat-UMa.dat')

! Read in T & O3 & CH4 & H2O climatology used to fill e.g. upper layers or if O3 not calc.
      call RD_PROF(JXUNIT,'atmos_std.dat')

      call RD_TRPROF(JXUNIT,'atmos_h2och4.dat')

        NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      enddo

! Read in photolysis rates used in chemistry code and mapping onto FJX J's
!---CTM call:  read in J-values names and link to fast-JX names
      !call RD_JS_JX(JXUNIT,'FJX_j2j.dat', TITLEJXX,NJXX)
      call RD_JS_JX(JXUNIT,'FJX_j2j_extended.dat', TITLEJXX,NJXX)

!---setup the random number sequence RAN4
        RANSEED = 66
      call RANSET (NRAN_,RAN4,RANSEED)


      END SUBROUTINE INIT_FJX


!-----------------------------------------------------------------------
      subroutine RD_XXX(NUN,NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-7.4  added Solar-J bins, W_=18 is only option, use NWBIN to zero strat fluxes
!>>>>NEW v-7.3  expanded input, full names & notes
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL
      integer  I, J, JJ, IW, LQ, NWWW, NSSS
      character*120  TIT_SPEC, TIT_J1N
      character*16 TIT_J1L
      character*6  TIT_J1S,TIT_J2S

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!   note that X_ = max # Xsects read in
!           NJX = # fast-JX J-values derived from this (.le. X_)
      if (W_ .ne. 18) then
       call EXITC(' no. wavelengths wrong: W_ .ne. 18')
      endif

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a120)',err=4) TIT_SPEC
      read (NUN,*,err=4)
      read (NUN,'(i5,5x,i5)',err=4) NWWW, NSSS
!print
       if (NWWW.gt.W_ .or. NSSS.gt.S_) then
        call EXITC(' W_ or S_ not large enough')
       endif

       NW1 = 1
       NW2 = NWWW
       NS1=1
       NS2=NSSS

!----w-params:  1=w-eff  2=w-bins, 3=solar(photons), 4=solar(W/m2), 4=Y-PAR,  5=Rayleigh, 
!               6=SJ sub-bins
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      !read (NUN,'(5x,6e10.3)',err=4)    (WL(IW),IW=1,NSSS)

      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)


      WL = (/   187. ,     191. ,     193.,      196. ,     202. ,     208., &
          211. ,     214. ,     261.,      267. ,     277. ,     295., &
          303. ,     310. ,     316.,      333. ,     383. ,     599., &
          973. ,    1267. ,    1448.,     1767. ,    2039. ,    2309., &
         2748. ,    3404. ,    5362./)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      !read (NUN,'(5x,6e10.3)',err=4)    (WBIN(IW),IW=1,NSSS)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)

      WBIN = (/       0.,        0. ,        0. ,        0.  , 0.,0., &
            0.,         0.,        0. ,        0. ,        0.  ,  0.29100,&
       0.29825,   0.30745,   0.31245 ,   0.32030 ,   0.34500  ,  0.41245,&
         0.778,     1.242,     1.299 ,     1.626 ,     1.942  ,    2.151,&
         2.500,     3.077,     3.846 ,     12.00 /)
        
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      !read (NUN,'(5x,6e10.3)',err=4)    (FL(IW),IW=1,NSSS)
            

      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)

      FL = (/    &
     1.391E+12_wp, 1.627E+12_wp, 1.664E+12_wp, 9.278E+11_wp,&
     7.842E+12_wp, 4.680E+12_wp, 9.918E+12_wp, 1.219E+13_wp,&
     6.364E+14_wp, 4.049E+14_wp, 3.150E+14_wp, 5.893E+14_wp,& 
     7.670E+14_wp, 5.041E+14_wp, 8.895E+14_wp, 3.852E+15_wp,&
     1.547E+16_wp, 1.805E+17_wp,       1.0_wp,       1.0_wp,&
           1.0_wp,       1.0_wp,       1.0_wp,       1.0_wp,& 
           1.0_wp,       1.0_wp,       1.0_wp /)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      !read (NUN,'(5x,6e10.3)',err=4)    (FW(IW),IW=1,NSSS)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)

      FW = (/ &
       0.01473,   0.01684,   0.01702,   0.00937,   0.07664,   0.04450,& 
       0.09299,   0.11278,     4.818,     2.962,     2.218,     3.703,& 
         4.670,     3.063,     5.414,     22.28,     77.17,    608.68,& 
        349.96,     25.59,    102.96,     56.01,     22.40,     23.50,& 
         20.20,     12.25,     12.58 /)


      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      !read (NUN,'(5x,6e10.3)',err=4)    (FP(IW),IW=1,NSSS)

      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)

      FP = (/ &
     0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp,&
     0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp,&
     0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp,&
     0.000E+00_wp, 0.000E+00_wp, 0.000E+00_wp, 5.138E-02_wp,&
     4.855E-01_wp, 6.760E-01_wp,       0.0_wp,       0.0_wp,&
           0.0_wp,       0.0_wp,       0.0_wp,       0.0_wp,& 
           0.0_wp,       0.0_wp,       0.0_wp /)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      !read (NUN,'(5x,6e10.3)',err=4)    (QRAYL(IW),IW=1,NSSS)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)
      read(NUN, *)

      QRAYL= (/&
     5.073E-25_wp, 4.479E-25_wp, 4.196E-25_wp, 3.906E-25_wp,&
     3.355E-25_wp, 2.929E-25_wp, 2.736E-25_wp, 2.581E-25_wp,&
     1.049E-25_wp, 9.492E-26_wp, 8.103E-26_wp, 6.135E-26_wp,& 
     5.424E-26_wp, 4.925E-26_wp, 4.516E-26_wp, 3.644E-26_wp,&
     2.082E-26_wp, 4.427E-27_wp, 5.380E-28_wp, 1.559E-28_wp,&
     9.578E-29_wp, 4.241E-29_wp, 2.347E-29_wp, 1.441E-29_wp,& 
     7.290E-30_wp, 3.117E-30_wp, 8.053E-31_wp /)

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      do I = 1,NSSS
        SJSUB(I, 1)=   1.0d0
        SJSUB(I, 2:15)=0.0d0
      enddo

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,1), (QO2(IW,3),IW=1,NWWW)
        TITLEJX(1) = TIT_J1S
        TITLEJL(1) = TIT_J1L
        LQQ(1) = 3

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,2), (QO3(IW,3),IW=1,NWWW)
        TITLEJX(2) = TIT_J1S
        TITLEJL(2) = TIT_J1L
        LQQ(2) = 3

      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,'(a6)',err=4) TIT_J2S
           if (TIT_J2s .ne. TIT_J1S) go to 4
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
                 TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)
        TITLEJX(3) = TIT_J1S
        TITLEJL(3) = TIT_J1L
        LQQ(3) = 3

!---Read remaining species:  X-sections at 1-2-3 T_s
!---read in 1 to 3 X-sects per J-value (JJ)
        JJ = 3
!-- read new Xsection block
    3 continue
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
        if (TIT_J1S .eq. 'endofJ') goto 1
!---try to add a new Xsect
    2 continue
       JJ = JJ+1
       LQ = 1
         if (JJ .gt. X_) call EXITC(' RD_XXX: X_ not large enough')
       TITLEJX(JJ) = TIT_J1S
       TITLEJL(JJ) = TIT_J1L
      read (NUN,'(a1,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
          SQQ(JJ),TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
!try to read a 2nd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N

        if (TIT_J1S .eq. 'endofJ') goto 1
      if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 2
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
!try to read a 3rd Temperature or Pressure
      read (NUN,'(a6,1x,a16,1x,a120)',err=4) TIT_J1S,TIT_J1L,TIT_J1N
         if (TIT_J1S .eq. 'endofJ') goto 1
       if (TIT_J1S .eq. TITLEJX(JJ)) then
        LQ = 3
      read (NUN,'(1x,f3.0,1x,6e10.3/5x,6e10.3/5x,6e10.3)',err=4)    &
        TQQ(LQ,JJ),(QQQ(IW,LQ,JJ),IW=1,NWWW)
        LQQ(JJ) = LQ
       else
        goto 2
       endif
      else
        goto 2
      endif
      goto 3
    4 continue
        call EXITC(' RD_XXX: error in read')
    1 continue
       NJX = JJ

!---read in complete, process Xsects for reduced wavelengths (Trop-Only)
!---    possibly also for WACCM >200nm-only version.
!---TROP-ONLY (W_ = 12 or 8) then drop the strat Xsects (labeled 'x')

      if (NWBIN .eq. 12 .or. NWBIN .eq. 8) then
        JJ = 3
        do J = 4,NJX
         if (SQQ(J) .ne. 'x') then
!---remove Xsects with 'x' label
          JJ = JJ+1
          if (JJ .lt. J) then
             TITLEJX(JJ) = TITLEJX(J)
             LQQ(JJ) = LQQ(J)
             SQQ(JJ) = SQQ(J)
           do LQ = 1,LQQ(J)
             TQQ(LQ,JJ) = TQQ(LQ,J)
            do IW = 1,NWWW
             QQQ(IW,LQ,JJ) = QQQ(IW,LQ,J)
            enddo
           enddo
          endif
         endif
        enddo
         NJX = JJ
      endif

!---need to check that TQQ (= T(K) or p(hPa)) is monotonically increasing:
      do J = 1,NJX
        if ((LQQ(J).eq.3) .and. (TQQ(2,J).ge.TQQ(3,J))) then
            call EXITC ('TQQ out of order')
        endif
        if ((LQQ(J).eq.2) .and. (TQQ(1,J).ge.TQQ(2,J))) then
            call EXITC ('TQQ out of order')
        endif
      enddo

!---zero strat SFluxes for TROP-ONLY (NWBIN = 12 or 8)
      if (NWBIN .eq. 12) then
        do IW = 1,4
          FL(IW) = 0.d0
        enddo
        do IW = 9,10
          FL(IW) = 0.d0
        enddo
      endif
      if (NWBIN .eq. 8) then
        do IW = 1,4
          FL(IW) = 0.d0
        enddo
          FL(5) = FL(5) * 2.d0
        do IW = 6,11
          FL(IW) = 0.d0
        enddo
      endif
      if (NWBIN.eq.12 .or. NWBIN.eq.8) then
      endif
      if (NSBIN .eq. 18) then
        do IW = 19,27
          FL(IW) = 0.d0
        enddo
      endif

      close(NUN)

      END SUBROUTINE RD_XXX


!-----------------------------------------------------------------------
      subroutine RD_CLD(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.4
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-cld.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, JCC
      character*120 TITLE0
      real(wp)     RCCJ,GCCJ, XNDR,XNDI

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

        read (NUN,'(a120)',err=4) TITLE0
        read (NUN,'(i4)')  NCC
        read (NUN,*)

      do K = 1,NCC
        read(NUN,'(4x,2f6.4,1x,a12)') DCC(K),RCC(K),TITLCC(K)
      enddo

      do J = 12,S_
        read (NUN,*)
       do K = 1,NCC
        read (NUN,'(i2,1x,3f5.3,f6.3,e8.1,f6.3,f8.5,7f6.3)',err=4) &
         JCC,WCC(J,K),RCCJ,GCCJ,XNDR,XNDI,QCC(J,K),SCC(J,K),(PCC(I,J,K),I=2,8)
         PCC(1,J,K) = 1.d0
       enddo
      enddo

! reproduce all cloud data for w < 295 nm from J=12 (= 295 nm), OK since trop clouds.
      do K = 1,NCC
        do J = 1,11
          WCC(J,K) = WCC(12,K)
          QCC(J,K) = QCC(12,K)
          SCC(J,K) = SCC(12,K)
         do I = 1,8
          PCC(I,J,K) = PCC(I,12,K)
         enddo
        enddo
      enddo

          goto 2

    4 continue
        call EXITC(' RD_CLD: error in read')

    2 continue
        close(NUN)


      END SUBROUTINE RD_CLD


!-----------------------------------------------------------------------
      subroutine RD_SSA(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX ver 7.4
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat-ssa.dat)
!     NUN      Channel number for reading data file
!     NSS      Number of categories for cloud scattering phase functions
!     QSS      Cloud scattering phase functions
!     WSS      5 Wavelengths for supplied phase functions
!     PSS      Phase function: first 8 terms of expansion
!     RSS      Effective radius associated with cloud type:  Integ(r^3 dr)/Integ(r^2 dr)
!     GSS      Effective geometric cross section:  Integ(pi r^2 dr)
!     SSS      Single scattering albedo
!     DSS      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, JSS, K
      character*120 TITLE0
      real(wp)     WJSS,XNDR,XNDI

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

        read (NUN,'(a120)',err=4) TITLE0
        read (NUN,*)
        read (NUN,'(i4)')  NSS
        read (NUN,*)

      do K = 1,NSS
        read (NUN,'(a12,5f8.2)')    &
             TITLSS(K),RSS(K),GSS(K),DSS(K),TSS(K),WSS(K)

       do J = 5,S_
         read(NUN,'(i2,2f8.4,e8.1,2f8.5,7f6.3)')      &
            JSS,WJSS,XNDR,XNDI,QSS(J,K),SSS(J,K),(PSS(I,J,K), I=2,8)
         PSS(1,J,K) = 1.d0
       enddo
      enddo

! reproduce all SSA data for J=1:4 with J=5
      do K = 1,NSS
        do J = 1,4
          QSS(J,K) = QSS(5,K)
          SSS(J,K) = SSS(5,K)
         do I = 1,8
          PSS(I,J,K) = PSS(I,5,K)
         enddo
        enddo
      enddo

          goto 2

    4 continue
        call EXITC(' RD_SSA: error in read')

    2 continue
        close(NUN)

      END SUBROUTINE RD_SSA


!-----------------------------------------------------------------------
      subroutine RD_MIE(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX ver 7.3+
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K , JAA
      character*120 TITLE0
!      character*12 TITLAA(A_)   ! TITLAA: Title for scattering data    NEEDS to be in COMMON
      Character*12 TITLAAJ
      real(wp)   RAAJ, DAAJ

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

      read (NUN,'(a120)',err=4) TITLE0
!print----
        read (NUN,*)
        read (NUN,*)
      do J = 1,A_
        read (NUN,'(i4,1x,a12,1x,2f6.3,1x,a120)',err=4) &
         JAA,TITLAAJ,RAAJ,DAAJ,TITLE0
       if (JAA.gt.0) then
         TITLAA(J) = TITLAAJ
         RAA(J) = RAAJ
         DAA(J) = DAAJ
        do K = 1,5
         read (NUN,'(f4.0,f7.4,f7.4,7f6.3)',err=4) &
          WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
         NAA = J
!print----
          ! J,TITLAAJ,RAAJ,DAAJ,(QAA(K,J),K=1,5),TITLE0
       else
          goto 2
       endif
      enddo
          goto 2

    4 continue
        call EXITC(' RD_MIE: error in read')
    2 continue
        close(NUN)

      END SUBROUTINE RD_MIE


!-----------------------------------------------------------------------
      subroutine RD_UM(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, L
      character*120 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      open (NUN,FILE=NAMFIL,status='old',form='formatted',err=4)

      read (NUN,'(a)',err=4) TITLE0
!print----
      read(NUN,'(5x,10f5.0)',err=4) WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
          read(NUN,'(a4)',err=4) TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
        do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NUN,'(18f9.5)',err=4)  ((UMAER(I,J,K,L),I=1,3),J=1,6)
        enddo
      enddo

      close(NUN)

!  collapse UM wavelengths, drop 550 nm
          WMM(4) = WMM(5)
          WMM(5) = WMM(6)
       do L=1,33
       do K=1,21
       do I=1,3
          UMAER(I,4,K,L) = UMAER(I,5,K,L)
          UMAER(I,5,K,L) = UMAER(I,6,K,L)
       enddo
       enddo
       enddo

!print----
      goto 2
    4 continue
        call EXITC(' RD_UM: error in read')
    2 continue

      END SUBROUTINE RD_UM


!-----------------------------------------------------------------------
      subroutine RD_PROF(NJ2,NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real(wp)  OFAC, OFAK

      character*78 TITLE0
!
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (T_REF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (O_REF(I,L,M), I=1,31)
      enddo
      close (NJ2)

!  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,LREF
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          O_REF(I,L,M) = O_REF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,LREF
        T_REF(I,L,M) = T_REF(41,L,M)
      enddo
      enddo
      enddo

 1000 format(1x,'std atmos profiles: T,O3',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_PROF


!-----------------------------------------------------------------------
      subroutine RD_TRPROF(NJ2,NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input H2O and CH4 reference profiles 'atmos_h2och4.dat'
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216

      character*78 TITLE0
!
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11E9.2)') (H2O_REF(I,L,M), I=1,31)
        read (NJ2,'(3X,11F9.2)') (CH4_REF(I,L,M), I=1,31)
      enddo
      close (NJ2)

!  Extend climatology to 100 km
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 32,LREF
        H2O_REF(I,L,M) = H2O_REF(31,L,M)
        CH4_REF(I,L,M) = CH4_REF(31,L,M)
      enddo
      enddo
      enddo

 1000 format(1x,'std atmos profiles: H2O,CH4',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_TRPROF


!-----------------------------------------------------------------------
      subroutine RD_JS_JX(NUNIT,NAMFIL,TITLEJX,NJX)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's JVMAP
!    including scaling factor JFACTA
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JVMAP
!           real(wp)  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JVMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JVMAP     label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
      implicit none
!
      integer, intent(in)                    ::  NUNIT, NJX
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) :: TITLEJX
      integer   J,JJ,K
      character*120 CLINE
      character*50 T_REACT
      character*6  T_FJX
      real(wp) F_FJX

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JVMAP(:) = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ .gt. JVN_) exit
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JVMAP(JJ) = T_FJX
        NRATJ = JJ
      enddo

      close(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
        if (JVMAP(K) .eq. TITLEJX(J)) then
         JIND(K) = J
        endif
       enddo
      enddo

      do K=1,NRATJ
       if (JVMAP(K) .ne. '------' ) then
        J = JIND(K)
       ! if (J.eq.0) then
       !        ' no mapping onto fast-JX',JVMAP(K)
       ! else
       !        ' mapped to FJX:',J,TITLEJX(J)
       ! endif
       endif
      enddo

      END SUBROUTINE RD_JS_JX


!-----------------------------------------------------------------------
      SUBROUTINE RANSET (ND,RAN4L,ISTART)
!-----------------------------------------------------------------------
!  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND)
!     program RAN3 from Press, based on Knuth
      implicit none
      integer, parameter ::  MBIG=1000000000
      integer, parameter ::  MSEED=161803398
      integer, parameter ::  MZ=0
      real(wp) , parameter ::  FAC=1.e-9_wp
      integer,intent(in)    :: ND
      real(wp), intent(out)   :: RAN4L(ND)
      integer,intent(inout) :: ISTART
      integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP
!---initialization and/or fix of ISEED < 0
        MJ = MSEED - abs(ISTART)
        MJ = mod(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        do I=1,54
          II = mod(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          if (MK.lt.MZ) then
            MK=MK+MBIG
          endif
          MJ = MA(II)
        enddo
        do K=1,4
         do I=1,55
           MA(I)=MA(I)-MA(1+MOD(I+30,55))
           if (MA(I) .lt. MZ) then
             MA(I) = MA(I)+MBIG
           endif
         enddo
        enddo
        INEXT = 0
        INEXTP = 31
        ISTART = 1
!---generate next ND pseudo-random numbers
      do J=1,ND
         INEXT = mod(INEXT,55) +1
         INEXTP = mod(INEXTP,55) +1
         MJ = MA(INEXT) - MA(INEXTP)
        if (MJ .lt. MZ) then
          MJ=MJ+MBIG
        endif
         MA(INEXT) = MJ
         RAN4L(J) = MJ*FAC
      enddo

      END SUBROUTINE RANSET

      END MODULE FJX_INIT_MOD

