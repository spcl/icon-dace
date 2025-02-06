! Definition of different levels of accuracy!
!
! Definition of different levels of accuracy!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: GPL-3.0-only  
! ---------------------------------------------------------------

MODULE messy_mecca_kpp_precision 
  USE mo_kind,                 ONLY: dp

!
! Definition of different levels of accuracy
! for REAL variables using KIND parameterization 
 USE mo_kind,                 ONLY: wp, dp, sp
!
! KPP SP - Single precision kind
 !INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
! KPP DP - Double precision kind
 !INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
! KPP QP - Quadruple precision kind
 !INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(18,400)

 END MODULE messy_mecca_kpp_precision 


