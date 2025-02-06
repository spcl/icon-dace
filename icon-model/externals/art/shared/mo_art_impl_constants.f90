!
! mo_art_impl_constants
! This module provides parameters that are required for the ART implementation
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

MODULE mo_art_impl_constants
! ICON
  USE mo_kind,                          ONLY: wp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: UNDEF_REAL_ART, UNDEF_INT_ART
  PUBLIC :: IART_VARNAMELEN, IART_XMLTAGLEN
  PUBLIC :: IART_AERO_TR, IART_CHEM_TR
  PUBLIC :: IART_MODE_1MOM, IART_MODE_RADIO, IART_MODE_POLL, IART_MODE_VOLC
  PUBLIC :: IART_MODE_2MOM, IART_EMISS2TRACER
  PUBLIC :: IART_MAX_RADIOACT, std_radioact_names
  PUBLIC :: IART_MAX_DIAGCONT, IART_ACC_DRYDEPO, IART_ACC_SEDIM, IART_ACC_EMISS, IART_EMISS
  PUBLIC :: IART_ACC_WETDEPO_GSCP, IART_ACC_WETDEPO_CON, IART_ACC_WETDEPO_RRSFC
  PUBLIC :: IART_LINOZ_ANA, IART_CHEM_NO
  PUBLIC :: IART_POLARCHEM, IART_LINOZ_LT
  PUBLIC :: IART_SIMNOY_SEDI
  PUBLIC :: IART_SIMNOY_PRES, IART_SIMNOY_WMO, IART_SIMNOY_EXTP
  PUBLIC :: IART_QV, IART_QC, IART_QI

! Determine undefined values
  REAL(wp), PARAMETER   :: &
    &  UNDEF_REAL_ART = -999._wp
  INTEGER,  PARAMETER   :: &
    &  UNDEF_INT_ART  = -999

! Character length definitions
  INTEGER,  PARAMETER    ::  &
    &  IART_VARNAMELEN = 50, &
    &  IART_XMLTAGLEN  = 400

! Tracer types
  INTEGER,  PARAMETER   :: &
    &  IART_AERO_TR   = 1, & !< Aerosol tracer
    &  IART_CHEM_TR   = 2    !< Chemical tracer

! Mode types
  INTEGER,  PARAMETER     :: &
    &  IART_MODE_1MOM  = 10, &
    &  IART_MODE_VOLC  = 11, &
    &  IART_MODE_RADIO = 12, &
    &  IART_MODE_POLL  = 13, &
    &  IART_MODE_2MOM  = 20
    
! Map types
  INTEGER, PARAMETER     :: &
    &  IART_EMISS2TRACER = 31

! Radioactive tracer names
  INTEGER,  PARAMETER     :: &
    &  IART_MAX_RADIOACT = 20
  CHARACTER(LEN=6), PARAMETER :: &
    &  std_radioact_names(9) =   &
    &     (/'Cs_137','I_131a',   &
    &       'Te_132','Zr_95 ',   &
    &       'Xe_133','I_131g',   &
    &       'I_131o','Ba_140',   &
    &       'Ru_103'/)

! Diagnostic container types
  INTEGER, PARAMETER      :: &
    &  IART_ACC_DRYDEPO       = 1, &
    &  IART_ACC_SEDIM         = 2, &
    &  IART_ACC_WETDEPO_GSCP  = 3, &
    &  IART_ACC_WETDEPO_CON   = 4, &
    &  IART_ACC_WETDEPO_RRSFC = 5, &
    &  IART_ACC_EMISS         = 6, &
    &  IART_EMISS             = 7, &
    &  IART_MAX_DIAGCONT      = 7

! Parameter for chemistry
  INTEGER, PARAMETER ::          &
    &  IART_CHEM_NO  = 0,        &
    &  IART_POLARCHEM = 1
! Parameters for water tracers
  INTEGER, PARAMETER :: &
    &  IART_QV = 1,     &
    &  IART_QC = 2,     &
    &  IART_QI = 3

! Parameters for linoz chemistry
  INTEGER, PARAMETER ::          &
    &  IART_LINOZ_ANA = 2,       &
    &  IART_LINOZ_LT = 3

! Parameters for Simnoy chemistry
  INTEGER, PARAMETER ::          &
    &  IART_SIMNOY_SEDI = 2,     &
    &  IART_SIMNOY_PRES = 3,     &
    &  IART_SIMNOY_WMO  = 4,     &
    &  IART_SIMNOY_EXTP = 5
END MODULE mo_art_impl_constants
