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

! Parameters for the VDIFF turbulence scheme.
!
! References:
!     Angevine, W. M., Jiang, H., & Mauritsen T. (2010).
!           Performance of an eddy diffusivity mass flux scheme for shallow cumulus boundary layers.
!           Monthly Weather Review, 138(7), 2895-2912. https://doi.org/10.1175/2010MWR3142.1
!     Mauritsen, T., & Svensson, G. (2007).
!           Observations of stably stratified shear-driven atmospheric turbulence at low and high Richardson numbers.
!           Journal of the Atmospheric Sciences, 64(2), 645-655. https://doi.org/10.1175/JAS3856.1

MODULE mo_turb_vdiff_params
  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: VDIFF_TURB_TTE = 1 !< TTE scheme, use in t_vdiff_config%turb
  INTEGER, PARAMETER :: VDIFF_TURB_3DSMAGORINSKY = 2 !< 3D Smagorinsky scheme, use in t_vdiff_config%turb

  ! Parameters
  REAL(wp),PARAMETER :: ckap    = 0.4_wp       !< karman constant.
  REAL(wp),PARAMETER :: cchar   = 0.018_wp     !< charnock constant.
  REAL(wp),PARAMETER :: cb      = 5._wp        !< stability parameter near neutrality.
  REAL(wp),PARAMETER :: cc      = 5._wp        !< stability parameter for unstable cases.

  REAL(wp),PARAMETER :: eps_shear = 1.e-5_wp   !< zepshr in sbr. vdiff of ECHAM6
  REAL(wp),PARAMETER :: eps_corio = 5.e-5_wp   !< zepcor in sbr. vdiff of ECHAM6
  REAL(wp),PARAMETER :: totte_min = 1.e-10_wp  !< minimum total turbulent energy

  REAL(wp),PARAMETER :: chneu  = 0.3_wp

  REAL(wp),PARAMETER :: cons5  = 3._wp*cb*cc

  REAL(wp),PARAMETER :: shn = 2.22_wp*0.22_wp*SQRT(2._wp)
  REAL(wp),PARAMETER :: smn = shn*1.24_wp*2.37_wp/3.69_wp
  REAL(wp),PARAMETER :: da1 = 1._wp/smn**3

  ! Parameters related to time step weighting in *rhs* of *vdiff* and *scv*

  REAL(wp),PARAMETER :: cvdifts = 1.5_wp
  REAL(wp),PARAMETER :: tpfac1  = cvdifts
  REAL(wp),PARAMETER :: tpfac2  = 1._wp / tpfac1
  REAL(wp),PARAMETER :: tpfac3  = 1._wp - tpfac2
  REAL(wp),PARAMETER :: tpfac4  = 1._wp + tpfac3

END MODULE mo_turb_vdiff_params
