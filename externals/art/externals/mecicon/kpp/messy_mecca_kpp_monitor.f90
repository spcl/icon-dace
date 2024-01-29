! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Utility Data Module File
!
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

MODULE messy_mecca_kpp_monitor 
  USE mo_kind,                 ONLY: dp 
 PUBLIC


  CHARACTER(LEN=32), PARAMETER, DIMENSION(14) :: SPC_NAMES = (/ &
     'N2O                             ','N2O5                            ','HO2                             ', &
     'H2O                             ','NO                              ','NO3                             ', &
     'HNO3                            ','O3P                             ','NO2                             ', &
     'OH                              ','O3                              ','O1D                             ', &
     'O2                              ','N2                              ' /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=32), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(22) :: EQN_NAMES = (/ &
     ' O1D + O2 --> O3P + O2                                                                              ', &
     ' O3P + O2 --> O3                                                                                    ', &
     ' O3P + O3 --> 2 O2                                                                                  ', &
     '  OH + O3 --> HO2 + O2                                                                              ', &
     ' HO2 + O3 --> OH + 2 O2                                                                             ', &
     ' HO2 + OH --> H2O + O2                                                                              ', &
     'H2O + O1D --> 2 OH                                                                                  ', &
     ' O1D + N2 --> O3P + N2                                                                              ', &
     'N2O + O1D --> 2 NO                                                                                  ', &
     '  NO + O3 --> NO2 + O2                                                                              ', &
     'O3P + NO2 --> NO + O2                                                                               ', &
     ' NO2 + O3 --> NO3 + O2                                                                              ', &
     'NO3 + NO2 --> N2O5                                                                                  ', &
     ' NO2 + OH --> HNO3                                                                                  ', &
     'HNO3 + OH --> H2O + NO3                                                                             ', &
     '       O2 --> 2 O3P                                                                                 ', &
     '       O3 --> O1D + O2                                                                              ', &
     '       O3 --> O3P + O2                                                                              ', &
     '      NO2 --> NO + O3P                                                                              ', &
     '      NO3 --> O3P + NO2                                                                             ', &
     '     N2O5 --> NO3 + NO2                                                                             ', &
     '     HNO3 --> NO2 + OH                                                                              ' /)

  CHARACTER(LEN=32), PARAMETER, DIMENSION(22) :: EQN_TAGS = (/ &
     'G1000                           ','G1001                           ','G1003                           ', &
     'G2104                           ','G2107                           ','G2109                           ', &
     'G2111                           ','G3101                           ','G3102a                          ', &
     'G3103                           ','G3105                           ','G3106                           ', &
     'G3109                           ','G3202                           ','G3206                           ', &
     'J1000a                          ','J1001a                          ','J1001b                          ', &
     'J3101                           ','J3103a                          ','J3104                           ', &
     'J3201                           ' /)

! INLINED global variables

! End INLINED global variables


 END MODULE messy_mecca_kpp_monitor 
