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

MODULE mo_radiation_general

  USE mo_kind, ONLY: wp

  PUBLIC

  INTEGER, PARAMETER :: nbndlw = 16, nbndsw = 14, seed_size = 4

  ! Shortwave spectral band limits (wavenumbers)
  REAL(wp), PARAMETER :: wavenum1(nbndsw) = (/ &
       2600._wp, 3250._wp, 4000._wp, 4650._wp, 5150._wp, 6150._wp, 7700._wp, &
       8050._wp,12850._wp,16000._wp,22650._wp,29000._wp,38000._wp,  820._wp/)
  REAL(wp), PARAMETER :: wavenum2(nbndsw) = (/ &
       3250., 4000., 4650., 5150., 6150., 7700., 8050., &
       12850.,16000.,22650.,29000.,38000.,50000., 2600./), &
    delwave(nbndsw)  = (/ &
       650.,  750.,  650.,  500., 1000., 1550.,  350., &
       4800., 3150., 6650., 6350., 9000.,12000., 1780./)

END MODULE mo_radiation_general
