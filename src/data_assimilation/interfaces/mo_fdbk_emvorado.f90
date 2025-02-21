! 3DVAR/COSMO source module for feedback file interface to COSMO
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

MODULE mo_fdbk_emvorado

!-------------------------------------------------------------------------------
!
! Description:
!   COSMO interface to write NetCDF feedobs (or feedback) file (FOF, 
!   Common format for FOF in 3DVAR and COSMO), especially for
!   the radar forward OPERATOR emvorado.
!
!   This module contains the following module procedures:
!     - write_report_radar_1
!     - write_report_radar_2
!     - add_data : interface for: add_inte_vala, add_real_vala_1D,
!                                 add_text_vala, add_real_vala_2D
!     - fill_bodybuf_int, fill_bodybuf_real
!
!   This module has been derived from the original "mo_fdbk_cosmo.f90"
!   and contains somewhat optimized clones of the original procedure
!   "write_report", namely "write_report_radar_1" and "write_report_radar_2".
!   These differ from the original in the sense of a perfomance-optimized
!   interface and vectorizable loop structures for radar data.
!
!-------------------------------------------------------------------------------


IMPLICIT  NONE


END MODULE mo_fdbk_emvorado
