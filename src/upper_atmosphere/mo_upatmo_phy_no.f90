!
! Computation of NO heating using the parameterization of Kockarts [GRL, 1980]
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

MODULE mo_upatmo_phy_no

  USE mo_kind,               ONLY: wp
  USE mo_impl_constants,     ONLY: SUCCESS
  USE mo_physical_constants, ONLY: ak, avo, argas 

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: no_heating

CONTAINS
  
  !>
  !! Compute heating rate due to radiative cooling at 5.3 x 10^-6 m by nitric oxide (NO)
  !!
  !! Literature:
  !! - Kockarts, G. (1980) Nitric oxide cooling in the terrestrial thermosphere. 
  !!   Geophys. Res. Lett., 7, 137-140.
  !!
  SUBROUTINE no_heating(jcs, jce, kbdim, klev, zo, zno, cp, tm1, apm1, amu, ptte, &
    &                   opt_istartlev, opt_iendlev)

    ! in/out variables
    INTEGER,  INTENT(IN)  :: jcs, jce, kbdim, klev
    REAL(wp), INTENT(IN)  :: zo(kbdim,klev)         ! volume mixing ratio (vmr) of O (in m3/m3)
    REAL(wp), INTENT(IN)  :: zno(kbdim,klev)        ! vmr of NO
    REAL(wp), INTENT(IN)  :: cp(kbdim,klev)         ! cp in J/K/kg
    REAL(wp), INTENT(IN)  :: tm1(kbdim,klev)        ! temperature in K
    REAL(wp), INTENT(IN)  :: apm1(kbdim,klev)       ! hydrostatic full level pressure in Pa
    REAL(wp), INTENT(IN)  :: amu(kbdim,klev)        ! mol mass of dry air in g/mol

    REAL(wp), INTENT(OUT) :: ptte(kbdim,klev)       ! temperature tendency K/s

    INTEGER, OPTIONAL, INTENT(IN)  :: opt_istartlev, opt_iendlev  ! optional vertical start and end indices

    ! local variables 
    REAL(wp)              :: nd_o                   ! number density (1/m3)
    REAL(wp)              :: nd_no                  ! number density (1/m3)
    REAL(wp)              :: en_rate                ! energy rate (W/m3)
    REAL(wp)              :: rho_air                ! air density (kg/m3)
    REAL(wp)              :: amc                    ! molar concentration of (dry) air (mol/m3)
    REAL(wp)              :: avo_amc
    REAL(wp)              :: inv_tm1                ! inverse temperature (1/K)

    INTEGER :: jk, jl, istartlev, iendlev

    REAL(wp), PARAMETER   :: hv  = 3.726e-20_wp     ! in Joule
    REAL(wp), PARAMETER   :: k10 = 6.5e-17_wp       ! O quenching in m^3 s^(-1)
    REAL(wp), PARAMETER   :: a10 = 13.3_wp          ! in s^(-1)
    !
    REAL(wp), PARAMETER   :: hv_k10_a10 = hv * k10 * a10
    REAL(wp), PARAMETER   :: n_hv_o_ak  = -hv / ak
    REAL(wp), PARAMETER   :: inv_argas  = 1._wp / argas

    !---------------------------------------------------------

    ! please do not limit range of assignment 
    ! (e.g., ptte(jcs:jce,istartlev:iendlev) = 0._wp)), 
    ! since tendencies have attribute INTENT(OUT)
    ptte(:, :) = 0._wp

    ! determine start and end indices of vertical grid layers,
    ! for which tendencies should be computed
    IF (PRESENT(opt_istartlev)) THEN
      istartlev = MIN(MAX(1, opt_istartlev), klev)
    ELSE
      istartlev = 1
    ENDIF

    IF (PRESENT(opt_iendlev)) THEN
      iendlev = MIN(MAX(1, opt_iendlev), klev)
    ELSE
      iendlev = klev
    ENDIF

    IF (istartlev > iendlev) RETURN 

    DO jk = istartlev, iendlev
      DO jl = jcs, jce

        inv_tm1 = 1._wp / tm1(jl,jk)
        
        ! (hydrostatic) molar concentration of air = pres / ( R_ideal * temp ), 
        ! with the ideal gas constant R_ideal -> argas
        amc = inv_argas * inv_tm1 * apm1(jl,jk)

        ! calculation of number densities of NO and O
        avo_amc = avo * amc
        nd_o    = avo_amc * zo(jl,jk)
        nd_no   = avo_amc * zno(jl,jk)

        ! calculation of air density
        ! (factor 10^(-3) is for: [amu] = g/mol -> kg/mol)
        rho_air = 1.E-3_wp * amc * amu(jl,jk)
        
        ! calculation of energy rate
        en_rate = ( ( hv_k10_a10 * nd_no * nd_o * EXP( n_hv_o_ak * inv_tm1 ) ) &
          &     / ( k10 * nd_o + a10 ) )
        
        ! calculation of NO heating rate
        ptte(jl,jk) = -en_rate / ( cp(jl,jk) * rho_air )
        
      ENDDO !jl
    ENDDO !jk
    
  END SUBROUTINE no_heating

END MODULE mo_upatmo_phy_no
