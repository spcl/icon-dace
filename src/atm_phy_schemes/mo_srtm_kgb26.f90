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

!option! -Nv -NO
#ifdef VPP
!OCL SCALAR
#endif
#ifdef __xlC__
@PROCESS NOOPTIMIZE
#endif
MODULE mo_srtm_kgb26
PUBLIC :: srtm_kgb26
CONTAINS

SUBROUTINE srtm_kgb26

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF

!     ------------------------------------------------------------------

USE mo_kind, ONLY : wp

USE mo_yoesrta26, ONLY : sfluxref, rayl

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
sfluxref = (/ &
 !  &     129.462_wp, 15*_ZERO_ /)
 & 29.0079_wp,  28.4088_wp,     20.3099_wp,  13.0283_wp &
 & ,  11.8619_wp,  9.95840_wp,     6.68696_wp,  5.38987_wp &
 & ,  3.49829_wp, 0.407693_wp,    0.299027_wp, 0.236827_wp &
 & , 0.188502_wp, 0.163489_wp, 4.64335e-02_wp, 2.72662e-03_wp /)

!     Rayleigh extinction coefficient at all v
rayl = (/ &
 & 1.21263e-06_wp,1.43428e-06_wp,1.67677e-06_wp,1.93255e-06_wp &
 & , 2.19177e-06_wp,2.44195e-06_wp,2.66926e-06_wp,2.85990e-06_wp &
 & , 3.00380e-06_wp,3.06996e-06_wp,3.08184e-06_wp,3.09172e-06_wp &
 & , 3.09938e-06_wp,3.10456e-06_wp,3.10727e-06_wp,3.10818e-06_wp /)

!     ------------------------------------------------------------------
END SUBROUTINE srtm_kgb26
END MODULE mo_srtm_kgb26
