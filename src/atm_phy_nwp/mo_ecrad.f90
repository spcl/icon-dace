!
! This module is the interface between ICON and ECMWFs radiation code
! ecRad which is provided as a library.
! - Modules, variables and data types from ecRad are used and provided
!   to ICON. Possibly ambiguous names get an ecrad_ prefix.
! - By this approach, no other ICON module has to make a direct use
!   statement on an ecRad module.
! - This module also holds a few ecRad configuration objects:
!   ecrad_conf, nweight_par_ecrad, iband_par_ecrad, weight_par_ecrad
!
! Literature references:
! Coddington et al. (2016) - Coddington, O., Lean, J. L., Pilewskie, P., Snow, M., & Lindholm, D. (2016).
!                            A solar irradiance climate data record.
!                            Bulletin of the American Meteorological Society, 97(7), 1265-1282.
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

MODULE mo_ecrad

  USE mo_kind,                    ONLY: wp

  IMPLICIT NONE

  PRIVATE



END MODULE mo_ecrad
