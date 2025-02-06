!
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

MODULE messy_cmn_photol_mem
  IMPLICIT NONE

  ! ip_* = index of photolysis
  INTEGER, PUBLIC  :: &
    ! Photolysis reactions currently implemented in CloudJ
    ip_NO       =  1,  &
    ip_O2       =  2,  &
    ip_O3P      =  3,  &
    ip_O1D      =  4,  &
    ip_CHOH     =  5,  &
    ip_COH2     =  6,  &
    ip_H2O2     =  7,  &
    ip_CH3OOH   =  8,  &
    ip_NO2      =  9,  &
    ip_NOO2     =  11, &
    ip_NO2O     =  12, &
    ip_N2O5     =  13, &
    ip_HONO     =  14, &
    ip_HNO3     =  15, &
    ip_HNO4     =  16, &
    ip_ClNO3    =  17, &
    ip_ClONO2   =  18, &
    ip_Cl2      =  19, &
    ip_HOCl     =  20, &
    ip_OClO     =  21, &
    ip_Cl2O2    =  22, &
    ip_BrO      =  24, &
    ip_BrNO3    =  25, & 
    ip_HOBr     =  26, &
    ip_BrCl     =  27, &
    ip_OCS      =  28, &
    ip_SO2      =  29, &
    ip_N2O      =  30, &
    ip_CFCl3    =  31, &
    ip_CF2Cl2   =  32, &
    ip_CCl4     =  36, &
    ip_CH3Cl    =  37, &
    ip_CH3CCl3  =  38, &
    ip_CH3Br    =  44, &
    ip_CF2ClBr  =  45, &
    ip_CH2Br2   =  48, &
    ip_CHBr3    =  49, &
    ip_CH3I     =  50, &
    ip_PAN      =  52, &
    ip_CH3NO3   =  53, &
    ip_CH3CHO   =  54, &
    ip_MVK      =  55, &
    ip_CF3Br    =  56, &
    ip_MACR     =  57, &
    ip_HOCH2CHO =  58, &
    ip_KET      =  61, &
    ip_MGLYOX   =  64, &
    ip_GLYOX    =  65, &
    ip_GLYOXb   =  65, &
    ip_GLYOXa   =  67, &
    ip_CH3COCH3 =  68, &
    ip_H2O      =  71, &
    ip_CO2      =  72, &

    ! indices not in CloudJ but derived from other photolysis rates in
    ! mo_mecicion.f90
    ip_Br2      =  73, &
    ip_ClNO2    =  74, &

    ! indices not implemented in mo_mecicon.f90 (or not available in CloudJ)
    ip_CH3CO3H  = 0, &
    ip_BrNO2    = 0, &
    ip_C3H7I    = 0, &
    ip_CH2ClI   = 0, &
    ip_CH2I2    = 0, &
    ip_IO       = 0, &
    ip_HOI      = 0, &
    ip_I2       = 0, &
    ip_ICl      = 0, &
    ip_IBr      = 0, &
    ip_INO2     = 0, &
    ip_INO3     = 0, &
    ip_SO3      = 0, &
    ip_CS2      = 0, &
    ip_HCl      = 0, &
    ip_CHCl2Br  = 0, &
    ip_CHClBr2  = 0, &
    ip_CH2ClBr  = 0, &
    ip_SF6      = 0, &
    ip_NO3NOO   = 0, & 
    ip_CH4      = 0, &
    ip_O2_b1b2  = 0, &
    ip_O2_b1    = 0, &
    ip_O2_b2    = 0, &
    ip_O3PO1D   = 0, &
    ip_O3Pp     = 0, &
    ip_H2O1D    = 0, &
    ip_N2       = 0, &
    ip_N2_b1    = 0, &
    ip_N2_b2    = 0, &
    ip_N2_b3    = 0, &
    ip_NN2D     = 0, &
    ip_NOp      = 0, &
    ip_Op_em    = 0, &
    ip_O2p_em   = 0, &
    ip_Op_O_em  = 0, &
    ip_N2p_em   = 0, &
    ip_Np_N_em  = 0, &
    ip_Np_N2D_em= 0, &
    ip_N_N2D_em = 0, &
    ip_Op_em_b  = 0, &
    ip_se_O2_b1 = 0, &
    ip_se_O2_b2 = 0, &
    ip_se_N2_b1 = 0, &
    ip_se_N2_b2 = 0, &
    ip_se_N2_b3 = 0, &
    ip_se_N2_b4 = 0, &
    ip_se_Op_em = 0, &
    ip_O2_aurq  = 0, &
    ip_N2_aurq  = 0, &
    ip_H2SO4    = 0, &
    ip_C3O2     = 0, &
    ip_CH3ONO   = 0, &
    ip_CH3O2NO2 = 0, &
    ip_CH3O2    = 0, &
    ip_HCOOH    = 0 
! IP_MAX must be set to the highest ip_* value from the definitions above:
  INTEGER, PUBLIC, PARAMETER :: IP_MAX = 111

  CHARACTER(LEN=9), PUBLIC, PARAMETER, DIMENSION(IP_MAX) :: jname = (/ &
    'O2       ', 'O3P      ', 'O1D      ', 'H2O2     ', &
    'NO2      ', 'NO2O     ', 'NOO2     ', 'N2O5     ', &
    'HNO3     ', 'HNO4     ', 'PAN      ', 'HONO     ', &
    'CH3OOH   ', 'COH2     ', 'CHOH     ', 'CH3CO3H  ', &
    'CH3CHO   ', 'CH3COCH3 ', 'MGLYOX   ', 'HOCl     ', &
    'OClO     ', 'Cl2O2    ', 'ClNO3    ', 'ClNO2    ', &
    'Cl2      ', 'BrO      ', 'HOBr     ', 'BrCl     ', &
    'BrNO3    ', 'BrNO2    ', 'Br2      ', 'CCl4     ', &
    'CH3Cl    ', 'CH3CCl3  ', 'CFCl3    ', 'CF2Cl2   ', &
    'CH3Br    ', 'CF2ClBr  ', 'CF3Br    ', 'CH3I     ', &
    'C3H7I    ', 'CH2ClI   ', 'CH2I2    ', 'IO       ', &
    'HOI      ', 'I2       ', 'ICl      ', 'IBr      ', &
    'INO2     ', 'INO3     ', 'SO2      ', 'SO3      ', &
    'OCS      ', 'CS2      ', 'H2O      ', 'N2O      ', &
    'NO       ', 'CO2      ', 'HCl      ', 'CHCl2Br  ', &
    'CHClBr2  ', 'CH2ClBr  ', 'CH2Br2   ', 'CHBr3    ', &
    'SF6      ', 'NO3NOO   ', 'ClONO2   ', 'MACR     ', &
    'MVK      ', 'GLYOX    ', 'HOCH2CHO ', 'CH4      ', &
    'O2_b1b2  ', 'O2_b1    ', 'O2_b2    ', 'O3PO1D   ', &
    'O3Pp     ', 'H2O1D    ', 'N2       ', 'N2_b1    ', &
    'N2_b2    ', 'N2_b3    ', 'NN2D     ', 'NOp      ', &
    'Op_em    ', 'O2p_em   ', 'Op_O_em  ', 'N2p_em   ', &
    'Np_N_em  ', 'Np_N2D_em', 'N_N2D_em ', 'Op_em_b  ', &
    'se_O2_b1 ', 'se_O2_b2 ', 'se_N2_b1 ', 'se_N2_b2 ', &
    'se_N2_b3 ', 'se_N2_b4 ', 'se_Op_em ', 'O2_aurq  ', &
    'N2_aurq  ', 'H2SO4    ', 'C3O2     ', 'KET      ', &
    'GLYOXa   ', 'GLYOXb   ', 'CH3NO3   ', 'CH3ONO   ', &
    'CH3O2NO2 ', 'CH3O2    ', 'HCOOH    '/)
END MODULE messy_cmn_photol_mem

!*****************************************************************************
