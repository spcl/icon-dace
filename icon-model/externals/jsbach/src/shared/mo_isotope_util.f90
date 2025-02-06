!> Routines for isotope discrimination (QUINCY)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Contains functions and subroutines to calculate isotope discrimination
!>

MODULE mo_isotope_util
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp
  USE mo_jsb_math_constants,   ONLY: eps8

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_delta_C13, calc_mixing_ratio_C13C12,   &
            calc_delta_C14, calc_mixing_ratio_C14C,     &
            calc_delta_N15, calc_mixing_ratio_N15N14,   &
            calc_fractionation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_isotope_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Calculates the isotopic enrichment of C13 in delta notation
  !!
  !! Input: molar masses of element and its isotope
  !!
  !! Output: delta (per mill)
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_delta_C13(element, isotope) RESULT(delta)

    USE mo_jsb_physical_constants, ONLY: PDB_STANDARD_C13

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: element             !< base element (including the isotope)
    REAL(wp), INTENT(in) :: isotope             !< isotope of element
    REAL(wp)             :: delta               !< isotope ratio against standard (per mill)
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_C13'


    IF(ABS(element - isotope) > eps8 .AND. (element * isotope) > 0._wp) THEN
       delta = ( isotope / (element-isotope) / PDB_STANDARD_C13 - 1.0_wp ) * 1000._wp
    ELSE
       delta = 0.0_wp
    ENDIF

  END FUNCTION calc_delta_C13


  !-----------------------------------------------------------------------------------------------------
  !> Calculates the isotopic mixing ratio of C13 for a given delta
  !!
  !! Input: delta C13 relative to PDB standard
  !!
  !! Output: molar mixing ratio of 13C to 12C
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_mixing_ratio_C13C12(delta) RESULT(mixing_ratio)

    USE mo_jsb_physical_constants, ONLY: PDB_STANDARD_C13

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: delta               !< isotope ratio in delta notation (per mill)
    REAL(wp)             :: mixing_ratio        !< isotopic mixing ratio 13C / 12C
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mixing_ratio_C13C12'


    mixing_ratio = ( delta / 1000._wp + 1.0_wp ) * PDB_STANDARD_C13

  END FUNCTION calc_mixing_ratio_C13C12


  !-----------------------------------------------------------------------------------------------------
  !> Calculates the isotopic enrichment of C14 in delta notation
  !!
  !! Input: molar masses of element and its isotope, as well as the corresponding delta C13
  !!
  !! Output: Delta (per mill)
  !!
  !! Based on equation 1-3 in Levin et al. (2010), Tellus 62B, 26-46
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_delta_C14(element, isotope, delta_C13) RESULT(delta)

    USE mo_jsb_physical_constants,    ONLY: AvogadroC_times_lambda_C14, molar_mass_C, AAB_STANDARD_C14, scale_correction_C14

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: element             !< base element (including the isotope) (mol)
    REAL(wp), INTENT(in) :: isotope             !< isotope of element (micro-mol)
    REAL(wp), INTENT(in) :: delta_C13           !< delta C13 of the pool/flux (per mill)
    REAL(wp)             :: Delta               !< C14 isotope ratio against standard, corrected for delta C13 (per mill)
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_C14'


    IF(element > eps8 .AND. (element * isotope) > 0._wp) THEN
       Delta = (AvogadroC_times_lambda_C14 / ( AAB_STANDARD_C14 * molar_mass_C ) * isotope / scale_correction_C14 / element * &
               (1._wp - (2._wp* (25._wp + delta_C13)) /1000._wp) - 1._wp) * 1000._wp
    ELSE
       Delta = 0.0_wp
    ENDIF

  END FUNCTION calc_delta_C14


  !-----------------------------------------------------------------------------------------------------
  !> Calculates the isotopic mixing ratio of C14
  !!
  !! Input: delta C13, and Delta C14
  !!
  !! Output: C14 mixing ratio (C14 against all C)
  !!
  !! Based on equation 1-3 in Levin et al. (2010), Tellus 62B, 26-46
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_mixing_ratio_C14C(delta_C13,Delta_C14) RESULT(mixing_ratio)

    USE mo_jsb_physical_constants,        ONLY: AvogadroC_times_lambda_C14, molar_mass_C, AAB_STANDARD_C14, scale_correction_C14

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: delta_C13           !< delta C13 of the pool/flux
    REAL(wp), INTENT(in) :: Delta_C14           !< Delta C14 of the pool/flux
    REAL(wp)             :: mixing_ratio        !< isotope molar mixing ratio 14C / C
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mixing_ratio_C14C'


    mixing_ratio = ( Delta_C14 / 1000._wp + 1._wp ) / &
                   ( ( AvogadroC_times_lambda_C14 / ( AAB_STANDARD_C14 * molar_mass_C ) ) * &
                     (1._wp - (2._wp* (25._wp + delta_C13)) /1000._wp) ) * &
                   scale_correction_C14

  END FUNCTION calc_mixing_ratio_C14C


  !-----------------------------------------------------------------------------------------------------
  !> Calculates the isotopic enrichment of N15 in delta notation
  !!
  !! Input: molar masses of element and its isotope
  !!
  !! Output: delta (per mill)
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_delta_N15(element, isotope) RESULT(delta)

    USE mo_jsb_physical_constants, ONLY: ATM_STANDARD_N15

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: element             !< base element (including the isotope)
    REAL(wp), INTENT(in) :: isotope             !< isotope of element
    REAL(wp)             :: delta               !< isotope ratio against standard (per mill)
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_N15'


    IF(ABS(element - isotope) > eps8 .AND. (element * isotope) > 0._wp) THEN
       delta = ( isotope / (element-isotope) / ATM_STANDARD_N15 - 1.0_wp ) * 1000._wp
    ELSE
       delta = 0.0_wp
    ENDIF

  END FUNCTION calc_delta_N15


  !-----------------------------------------------------------------------------------------------------
  !> Calculates the isotopic mixing ratio of N15 for a given delta
  !!
  !! Input: delta N15 relative to ATM standard
  !!
  !! Output: molar mixing ratio of 15N to 14N
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_mixing_ratio_N15N14(delta) RESULT(mixing_ratio)

    !$ACC ROUTINE SEQ

    USE mo_jsb_physical_constants, ONLY: ATM_STANDARD_N15

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: delta               !< isotope ratio in delta notation (per mill)
    REAL(wp)             :: mixing_ratio        !< isotopic mixing ratio 13C / 12C
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mixing_ratio_N15N14'


    mixing_ratio = ( delta / 1000._wp + 1.0_wp ) * ATM_STANDARD_N15

  END FUNCTION calc_mixing_ratio_N15N14


  !-----------------------------------------------------------------------------------------------------
  !> Calculates fractionation of a processes given the process' discrimination value
  !!
  !! Input: molar masses of the sum of all isotopes of an element and the isotope of the source pool,
  !!   as well as the process discrimination
  !!
  !! Output: fraction of the transfer flux that will be the isotope
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_fractionation(source_element,source_isotope,discrimination,source_isotope2) RESULT(fractionation)

     IMPLICIT NONE
     ! ---------------------------
     ! 0.1 InOut
     REAL(wp), INTENT(in)           :: source_element     !< mol of element (sum of all isotopes) in source material
     REAL(wp), INTENT(in)           :: source_isotope     !< mol of the discriminated isotope in source material
     REAL(wp), INTENT(in)           :: discrimination     !< discrimination of transfer process (per mill)
     REAL(wp), INTENT(in), OPTIONAL :: source_isotope2    !< mol of the co-occurring discriminated isotope in source material (13C)
     REAL(wp)                       :: fractionation      !< fraction of the total flux that is isotope
     ! ---------------------------
     ! 0.2 Local
     REAL(wp)                    :: Rsource,Rsink,hlp1,hlp2
     CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fractionation'


     !> 1.0 normal mixing model
     !!
     IF(.NOT.PRESENT(source_isotope2)) THEN
        ! molar mixing ratio of the source pool
        IF((source_element - source_isotope) > eps8) THEN
           Rsource = source_isotope / ( source_element - source_isotope )
        ELSE
           Rsource = 0.0_wp
        ENDIF

        ! molar mixing ratio of the isotope with the main element in the product of the reaction
        ! from Discrimination = ( Rsource/Rsink - 1 ) * 1000
        Rsink = Rsource / (discrimination/1000._wp + 1._wp)

        ! Since here the processes fluxes F are represented as the sum of all isotopes of an element
        ! (main: E and its isotope: I)
        ! F = E + I, with I/E = Rsink follows that
        !   = I/Rsink + I = ( 1/ Rsink + 1 ) * I => I = F * 1 / ( 1 + 1 / Rsink )
        IF( Rsink > eps8) THEN
          fractionation = 1._wp / (1._wp + 1._wp/Rsink)
        ELSE
          fractionation = 0.0_wp
        ENDIF
     !> 2.0  ..
     !!
     ELSE
        ! source dC13
        hlp1 = calc_delta_C13(source_element,source_isotope2)
        ! source DC14
        hlp2 = calc_delta_C14(source_element,source_isotope,hlp1)
        ! mixing ratio dependent on input C13 and C14 - discrimination
        fractionation = calc_mixing_ratio_C14C(hlp1,hlp2-discrimination)
     ENDIF
  END FUNCTION calc_fractionation

#endif
END MODULE mo_isotope_util
