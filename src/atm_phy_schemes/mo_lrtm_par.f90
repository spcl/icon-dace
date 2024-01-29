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


MODULE mo_lrtm_par

  USE mo_kind, ONLY : wp

  IMPLICIT NONE


  PUBLIC
  SAVE


  !------------------------------------------------------------------
  ! rrtmg_lw main parameters
  !
  ! Initial version:  JJMorcrette, ECMWF, Jul 1998
  ! Revised: MJIacono, AER, Jun 2006
  ! Revised: MJIacono, AER, Aug 2007
  ! Revised: MJIacono, AER, Aug 2008
  !------------------------------------------------------------------

  !  name     type     purpose
  ! -----  :  ----   : ----------------------------------------------
  ! mxlay  :  integer: maximum number of layers
  ! mg     :  integer: number of original g-intervals per spectral band
  ! nbndlw :  integer: number of spectral bands
  ! maxxsec:  integer: maximum number of cross-section molecules
  !                    (e.g. cfcs)
  ! maxinpx:  integer:
  ! ngptlw :  integer: total number of reduced g-intervals for rrtmg_lw
  ! ngNN   :  integer: number of reduced g-intervals per spectral band
  ! ngsNN  :  integer: cumulative number of g-intervals per band
  !------------------------------------------------------------------

  INTEGER, PARAMETER :: mxlay  = 203
  INTEGER, PARAMETER :: mg     = 16
  INTEGER, PARAMETER :: nbndlw = 16
  INTEGER, PARAMETER :: maxxsec= 4
  INTEGER, PARAMETER :: mxmol  = 38
  INTEGER, PARAMETER :: maxinpx= 38
  INTEGER, PARAMETER :: nmol   = 7

  INTEGER, PARAMETER :: ngptlw = 140 ! set to 256 for 256 gpt model
  INTEGER, PARAMETER :: ng1  = 10
  INTEGER, PARAMETER :: ng2  = 12
  INTEGER, PARAMETER :: ng3  = 16
  INTEGER, PARAMETER :: ng4  = 14
  INTEGER, PARAMETER :: ng5  = 16
  INTEGER, PARAMETER :: ng6  = 8
  INTEGER, PARAMETER :: ng7  = 12
  INTEGER, PARAMETER :: ng8  = 8
  INTEGER, PARAMETER :: ng9  = 12
  INTEGER, PARAMETER :: ng10 = 6
  INTEGER, PARAMETER :: ng11 = 8
  INTEGER, PARAMETER :: ng12 = 8
  INTEGER, PARAMETER :: ng13 = 4
  INTEGER, PARAMETER :: ng14 = 2
  INTEGER, PARAMETER :: ng15 = 2
  INTEGER, PARAMETER :: ng16 = 2

  INTEGER, PARAMETER :: ngs1  = 10
  INTEGER, PARAMETER :: ngs2  = 22
  INTEGER, PARAMETER :: ngs3  = 38
  INTEGER, PARAMETER :: ngs4  = 52
  INTEGER, PARAMETER :: ngs5  = 68
  INTEGER, PARAMETER :: ngs6  = 76
  INTEGER, PARAMETER :: ngs7  = 88
  INTEGER, PARAMETER :: ngs8  = 96
  INTEGER, PARAMETER :: ngs9  = 108
  INTEGER, PARAMETER :: ngs10 = 114
  INTEGER, PARAMETER :: ngs11 = 122
  INTEGER, PARAMETER :: ngs12 = 130
  INTEGER, PARAMETER :: ngs13 = 134
  INTEGER, PARAMETER :: ngs14 = 136
  INTEGER, PARAMETER :: ngs15 = 138

  !------------------------------------------------------------------
  ! rrtmg_lw spectral information

  ! Initial version:  JJMorcrette, ECMWF, jul1998
  ! Revised: MJIacono, AER, jun2006
  ! Revised: MJIacono, AER, aug2008
  !------------------------------------------------------------------

  !  name     type     purpose
  ! -----  :  ----   : ----------------------------------------------
  ! ng     :  integer: Number of original g-intervals in each spectral band
  ! nspa   :  integer: For the lower atmosphere, the number of reference
  !                    atmospheres that are stored for each spectral band
  !                    per pressure level and temperature.  Each of these
  !                    atmospheres has different relative amounts of the
  !                    key species for the band (i.e. different binary
  !                    species parameters).
  ! nspb   :  integer: Same as nspa for the upper atmosphere
  !wavenum1:  real   : Spectral band lower boundary in wavenumbers
  !wavenum2:  real   : Spectral band upper boundary in wavenumbers
  ! delwave:  real   : Spectral band width in wavenumbers
  ! totplnk:  real   : Integrated Planck value for each band; (band 16
  !                    includes total from 2600 cm-1 to infinity)
  !                    Used for calculation across total spectrum
  !totplk16:  real   : Integrated Planck value for band 16 (2600-3250 cm-1)
  !                    Used for calculation in band 16 only if
  !                    individual band output requested
  !totplnkderiv: real: Integrated Planck function derivative with respect
  !                    to temperature for each band; (band 16
  !                    includes total from 2600 cm-1 to infinity)
  !                    Used for calculation across total spectrum
  !totplk16deriv:real: Integrated Planck function derivative with respect
  !                    to temperature for band 16 (2600-3250 cm-1)
  !                    Used for calculation in band 16 only if
  !                    individual band output requested
  !
  ! ngc    :  integer: The number of new g-intervals in each band
  ! ngs    :  integer: The cumulative sum of new g-intervals for each band
  ! ngm    :  integer: The index of each new g-interval relative to the
  !                    original 16 g-intervals in each band
  ! ngn    :  integer: The number of original g-intervals that are
  !                    combined to make each new g-intervals in each band
  ! ngb    :  integer: The band index for each new g-interval
  ! wt     :  real   : RRTM weights for the original 16 g-intervals
  ! rwgt   :  real   : Weights for combining original 16 g-intervals
  !                    (256 total) into reduced set of g-intervals
  !                    (140 total)
  ! nxmol  :  integer: Number of cross-section molecules
  ! ixindx :  integer: Flag for active cross-sections in calculation
  !------------------------------------------------------------------

  INTEGER :: ng(nbndlw)
  INTEGER :: nspa(nbndlw)
  INTEGER :: nspb(nbndlw)

  REAL(wp) :: wavenum1(nbndlw)
  REAL(wp) :: wavenum2(nbndlw)
  REAL(wp) :: delwave(nbndlw)

  REAL(wp) :: totplnk(181,nbndlw)
  REAL(wp) :: totplk16(181)

  REAL(wp) :: totplnkderiv(181,nbndlw)
  REAL(wp) :: totplk16deriv(181)

  INTEGER :: ngc(nbndlw)
  INTEGER :: ngs(nbndlw)
  INTEGER :: ngn(ngptlw)
  INTEGER :: ngb(ngptlw)
  INTEGER :: ngm(nbndlw*mg)

  REAL(wp) :: wt(mg)
  REAL(wp) :: rwgt(nbndlw*mg)

  INTEGER :: nxmol
  INTEGER :: ixindx(maxinpx)

END MODULE mo_lrtm_par

