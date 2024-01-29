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
MODULE mo_snow_ice_reff

  USE mo_kind              , ONLY: wp
  USE mo_physical_constants, ONLY: rhoi, tmelt
  USE mo_aes_graupel       , ONLY: snow_number, snow_lambda, ice_number

  IMPLICIT NONE
  PRIVATE
  PUBLIC                         :: init_reff_funeedles, snow_reff_funeedles, &
                                  & ice_reff_funeedles, snow_x, ice_x, ice_reff_moss

  REAL(wp), PARAMETER            :: a_geo_snow=0.069_wp,& ! Formfactor in the mass-size relation of snow particles
                                  & a_geo_ice=130._wp     ! and ice particles
                                                          ! for graupel scheme or cloud ice scheme
  REAL(wp), PARAMETER            :: b_geo_snow=2.000_wp,& ! Exponent in the mass-size relation of snow particles
                                  & b_geo_ice=3.000_wp    ! and ice particles
  REAL(wp), PARAMETER            :: nu   = 0.0_wp,  & ! Marshall Palmer distribution (exponential)
                                    mu   = 1.0_wp
  REAL(wp)                       :: reff_coeff_fu_ice(1:5), reff_coeff_fu_snow(1:5), x_min, x_max
  !$ACC DECLARE CREATE(reff_coeff_fu_ice, reff_coeff_fu_snow)

CONTAINS

  SUBROUTINE init_reff_funeedles
    REAL(wp)                     :: bf1, bf3

    ! Fu Random Hexagonal needles:  Dge = 1/(c1 * x**[c2] + c3 * x**[c4])
    ! Parameterization based on Fu, 1996; Fu et al., 1998; Fu ,2007

    ! ice particles (monodisperse)

    reff_coeff_fu_ice(5)=0.5_wp

    reff_coeff_fu_ice(1) = SQRT( 3.0_wp *SQRT(3.0_wp) * rhoi / 8.0_wp ) *a_geo_ice**(-1.0_wp/2.0_wp/b_geo_ice)
    reff_coeff_fu_ice(2) = (1.0_wp-b_geo_ice)/2.0_wp/b_geo_ice
    reff_coeff_fu_ice(3) = SQRT(3.0_wp)/4.0_wp*a_geo_ice**(1.0_wp/b_geo_ice)
    reff_coeff_fu_ice(4) = -1.0_wp/b_geo_ice

    ! snow particles (not monodisperse)

    reff_coeff_fu_snow(5)=0.5_wp

    reff_coeff_fu_snow(1) = SQRT( 3.0_wp *SQRT(3.0_wp) * rhoi / 8.0_wp ) *a_geo_snow**(-1.0_wp/2.0_wp/b_geo_snow)
    reff_coeff_fu_snow(2) = (1.0_wp-b_geo_snow)/2.0_wp/b_geo_snow
    reff_coeff_fu_snow(3) = SQRT(3.0_wp)/4.0_wp*a_geo_snow**(1.0_wp/b_geo_snow)
    reff_coeff_fu_snow(4) = -1.0_wp/b_geo_snow

    ! Broadening for not monodisperse. Generalized gamma distribution
    bf1  =  GAMMA( ( b_geo_snow + 2.0_wp * nu + 3.0_wp)/mu/2.0_wp ) / GAMMA( (b_geo_snow + nu + 1.0_wp)/ mu) * &
         & ( GAMMA( (nu + 1.0_wp)/mu) / GAMMA( (b_geo_snow + nu + 1.0_wp)/mu) )** &
         & ( (1.0_wp-b_geo_snow)/2.0_wp/b_geo_snow)

    bf3 =  GAMMA( (b_geo_snow + nu )/mu ) / GAMMA( (b_geo_snow + nu + 1.0_wp)/mu) * &
         & ( GAMMA( (nu + 1.0_wp)/mu ) / GAMMA( (b_geo_snow + nu + 1.0_wp)/mu) )**( -1.0_wp/b_geo_snow)

    reff_coeff_fu_snow(1) = reff_coeff_fu_snow(1)*bf1
    reff_coeff_fu_snow(3) = reff_coeff_fu_snow(3)*bf3
    !$ACC UPDATE DEVICE(reff_coeff_fu_ice, reff_coeff_fu_snow)
  END SUBROUTINE init_reff_funeedles

  ELEMENTAL FUNCTION snow_reff_funeedles(temp, rho, qs)
  !$ACC ROUTINE SEQ
    REAL(wp)                          :: snow_reff_funeedles
    REAL(wp), INTENT(IN)              :: temp, rho, qs
    REAL(wp)                          :: x

    x=snow_x(temp, rho, qs)
    snow_reff_funeedles=reff_funeedles_snow(x)
  END FUNCTION snow_reff_funeedles

  ELEMENTAL FUNCTION snow_x(temp, rho, qs)
  !$ACC ROUTINE SEQ
    ! average mass of snow flake
    REAL(wp)                          :: snow_x
    REAL(wp), INTENT(IN)              :: temp, rho, qs
    REAL(wp)                          :: ncn
    REAL(wp), PARAMETER               :: eps = 1.0e-8    ! Epsilon constant    

    ! number of snow flakes per volume air given by N_0/lambda
    ncn=snow_number(temp,rho,qs) 
    ncn=ncn/snow_lambda(rho,qs,ncn) 
    snow_x=rho*qs/(ncn+eps) ! average mass of snow flake in [kg]
  END FUNCTION snow_x

  ELEMENTAL FUNCTION ice_reff_funeedles(temp, rho, qi)
  !$ACC ROUTINE SEQ
    ! effective radius of ice crystals for optical properties according to
    ! Fu needles in [m]
    REAL(wp)                          :: ice_reff_funeedles
    REAL(wp), INTENT(IN)              :: temp, rho, qi
    REAL(wp)                          :: x

    x=ice_x(temp, rho, qi)
    ice_reff_funeedles=reff_funeedles_ice(x)
  END FUNCTION ice_reff_funeedles
  ELEMENTAL FUNCTION ice_x(temp, rho, qi)
  !$ACC ROUTINE SEQ
    ! average mass of ice crystal in kg
    REAL(wp)                          :: ice_x
    REAL(wp), INTENT(IN)              :: temp, rho, qi
    REAL(wp)                          :: ncn
    REAL(wp), PARAMETER               :: eps = 1.0e-8    ! Epsilon constant

    ice_x=0._wp
    IF (temp < tmelt) THEN
      ncn=ice_number(temp,rho) ! number of ice crystals per mass air
      ice_x=qi/(ncn+eps)   ! [qi]: specific ice mass (mass of ice per mass air)
    END IF
  END FUNCTION ice_x
  ELEMENTAL FUNCTION reff_funeedles_snow(x)
  !$ACC ROUTINE SEQ
    ! effective radius of solid hydrometeors according to the Fu-formula for snow crystals in
    ! form of randomly oriented needles.
    ! reff_funeedles= c5/(c1*x**c2 + c3*x**c4)
    ! x average mass of crystal in kg
    REAL(wp)                          :: reff_funeedles_snow ! effective radius in [m]
    REAL(wp), INTENT(IN)              :: x              ! average mass of crystal in [kg]

    IF (x > 0._wp) THEN
       reff_funeedles_snow=reff_coeff_fu_snow(5)/( reff_coeff_fu_snow(1) * x**reff_coeff_fu_snow(2) + &
            & reff_coeff_fu_snow(3) * x**reff_coeff_fu_snow(4) )
    ELSE
       reff_funeedles_snow=1.0e-8
    END IF
  END FUNCTION reff_funeedles_snow
  ELEMENTAL FUNCTION reff_funeedles_ice(x)
  !$ACC ROUTINE SEQ
    ! effective radius of solid hydrometeors according to the Fu-formula for ice crystals in
    ! form of randomly oriented needles.
    ! reff_funeedles= c5/(c1*x**c2 + c3*x**c4)
    ! x average mass of crystal in kg
    REAL(wp)                          :: reff_funeedles_ice ! effective radius in [m]
    REAL(wp), INTENT(IN)              :: x              ! average mass of crystal in [kg]

    IF (x > 0._wp) THEN
       reff_funeedles_ice=reff_coeff_fu_ice(5)/( reff_coeff_fu_ice(1) * x**reff_coeff_fu_ice(2) + &
            & reff_coeff_fu_ice(3) * x**reff_coeff_fu_ice(4) )
    ELSE
       reff_funeedles_ice=1.0e-8
    END IF
  END FUNCTION reff_funeedles_ice

  ELEMENTAL FUNCTION ice_reff_moss(temp, rho, qi)
  !$ACC ROUTINE SEQ
    ! Effective radius of ice crystals (in [m]).
    !
    ! Reference:
    ! Moss, S. J., Francis, P. N., & Johnson, D. G. (1996).
    ! Calculation and parameterization of the effective radius of ice particles using aircraft data.
    ! In Proc. 12th Int. Conf. on Clouds and Precipitation (pp. 1255-1258).

    REAL(wp)                          :: ice_reff_moss
    REAL(wp), INTENT(IN)              :: temp, rho, qi  ! temp is unused
    REAL(wp)                          :: iwc

    iwc = 1000.0_wp * qi * rho  ! [g / m^3]
    IF (iwc > 0._wp) THEN
        ice_reff_moss = 83.8_wp * iwc**0.216_wp * 1.0e-6_wp  ! [m]
    ELSE
        ice_reff_moss = 1.0e-8
    END IF
  END FUNCTION ice_reff_moss

END MODULE mo_snow_ice_reff
