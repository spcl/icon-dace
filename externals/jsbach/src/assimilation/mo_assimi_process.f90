!> Contains the routines for the assimilation processes
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
MODULE mo_assimi_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_assimilation_scaling_factors, calc_assimilation_waterunlimited, &
    &       calc_assimilation_waterlimited, calc_NPP_pot_rate, calc_NPPbuf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_process'

CONTAINS

  SUBROUTINE calc_assimilation_scaling_factors( &
    & declination,                                        & ! in
    & canopy_bound_lai,                                   & ! in
    & lai,                                                & ! in
    & lat,                                                & ! in
    & scaling_fact_cl                                     & ! out
    & )

    ! Declarations
    USE mo_rad_constants,      ONLY: LaiLimit
    USE mo_jsb_math_constants, ONLY: deg2rad

    ! Arguments
    REAL(wp), INTENT(in) :: &
      & declination,        & ! Solar Declination angle
      & canopy_bound_lai
    REAL(wp), INTENT(in), DIMENSION(:) :: &
      & lai,                & ! Leaf area index [-]
      & lat

    REAL(wp), INTENT(out), DIMENSION(:) :: &
      & scaling_fact_cl

    ! Local variables
    REAL(wp) ::         &
      & k12,            &
      & cos_zenith_noon

    INTEGER :: nc, ic

    nc = SIZE(lai)

    !$ACC DATA PRESENT(lai, lat, scaling_fact_cl)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(cos_zenith_noon, k12)
    DO ic=1,nc

      ! initialize local arguments and variables
      scaling_fact_cl(ic) = 1._wp

      ! Cosine of zenith angle at local noon:
      cos_zenith_noon = COS(declination) * COS(deg2rad * lat(ic)) +  SIN(declination) * SIN(deg2rad * lat(ic))

      ! Extinction factor
      k12 = 0.5_wp / MAX(cos_zenith_noon, 1.E-3_wp) ! K = 1 / 2mue , Knorr (119c)

      ! R: Todo when offline and online versions exist...
      ! vg: in the current echam version inquire_declination at radiation time steps returns the declination at radiation time
      !     and not the current declination. This leads to slightly different result for offline and online jsbach simulations
      !     (https://code.zmaw.de/issues/4803). For that reason we use a k12 independent of the declination for the interface test.
      ! IF (interface_test) k12(:)=1._wp ! interface_test=.FALSE./.TRUE.

      ! Condition: LAI>LaiLimit
      IF (lai(ic) >= LaiLimit) THEN
        ! R: Note this could be packed into an ELEMENTAL procedure, but I think it would make things only complicated without
        ! improving security
        scaling_fact_cl(ic) = MAX(1.e-10_wp,EXP( -k12 * canopy_bound_lai * lai(ic))) !! see Eqs. (107), (108) in Knorr
      END IF

    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE calc_assimilation_scaling_factors


  !> ! Coomputes canopy conductance per canopy layer for water-unlimited case
  REAL(wp) FUNCTION calc_assimilation_waterunlimited( &
    & C4Flag,                                            & ! in
    & CarboxRate,                                        & ! in
    & ETransport,                                        & ! in
    & t_air,                                             & ! in
    & press_srf,                                         & ! in
    & par_down_mol,                                      & ! in
    & apar_per_lai_cl,                                   & ! in
    & CO2_mol_mixing_ratio,                              & ! in
    & N_scaling_factors                                  & ! in
    & ) RESULT(canopy_cond_cl)

    !$ACC ROUTINE SEQ

    !-----------------------------------------------------------------------
    USE mo_assimi_constants                                     ! we need all of them, except of f_aut_leaf and cCost.
    USE mo_jsb_physical_constants, ONLY: tmelt, R_gas => argas  ! tmelt =273.15_wp;  argas =Universal gas constant [J/(mole*K)]

    LOGICAL,  INTENT(in) ::   &
      & C4Flag
    REAL(wp), INTENT(in) ::   &
      & CarboxRate,           &
      & ETransport
    REAL(wp), INTENT(in) ::   &
      & t_air,                &  ! Atmosphere temperature (lowest layer) in Kelvin
      & press_srf,            &  ! Surface pressure
      & par_down_mol,         &  ! Downward PAR flux in mol (photons)/(m^2 s)
      & apar_per_lai_cl,      &  ! Absorbed PAR of canopy layer.
      & CO2_mol_mixing_ratio, &  ! CO2 substance mixing ratio [molCO2/molDryAir]
      & N_scaling_factors        ! Enzymatic scaling with canopy height

    ! Local variables

    REAL(wp) ::           &
      & CO2_conc_leaf_cl, &
      & T0,               & ! Temperature relative to 25 degree Celcius, i.e. T - 25
      & T1,               & ! 25 degress Celsius expressed in Kelvin
      & KC,               & ! Michaelis-Menten constant for CO2
      & KO,               & ! Michaelis-Menten constant for O2
      & GAM,              & ! COMPENSATION POINT WITHOUT DARK RESPIRATION [MOL(CO2) / MOL(AIR)]
      & HITINHIB,         & ! Hight temperature inhibition
                            ! Value to inhibit Assimilation and Respiration at temperatures above
                            ! 55 Celsius from Collatz et al., Physiological and environmental regulation
                            ! of stomatal conductance, photosynthesis and transpiration: a model that
                            ! includes a laminar boundary layer,
                            ! Agricultural and Forest Meteorology, 54, pp. 107-136, 1991
      & DARKINHIB           ! Dark inhibition
                            ! Value to inhibit Dark-Respiration in light
                            ! after Brooks and Farquhar, Effect of temperature on the CO2/O2 specifity on RBisCO
                            ! and the rate of respiration in the light, Planta 165, 397-406, 1985
                            ! It is fitted to inhibit the dark-respiration to 50% of it's uninhibited value
                            ! up from 50 umol/m^2s total solar downward radiation, i.e. 25 umol/m^2s downward PAR
                            ! (idea for improvement: darkinhib should be closer to 1 at very high latitudes)
                            ! The FUNCTION is an own creation.

    REAL(wp) ::           & ! variables with canopy layer
      & VCMAX,            & ! MAXIMUM CARBOXYLATION RATE [MOL(CO2)/M^2 S]
      & JMAX,             & ! MAXIMUM RATE OF ELECTRON TRANSPORT [MOL(CO2)/M^2 S]
      & K,                & ! CO2 SPECIFICITY OF PECASE
      & TC,               & ! Canopy Temperature (leaf) in Celsius.
      & GASS,             & ! gross assimilation (formerly called A)
      & DARK_RESP,        & ! dark respiration (formerly called RD).
      & JE,               & ! Electron transport rate (Knorr, Eq. 102c)
      & JC,               & ! Carboxylation rate (Knorr, Eq. 102b)
      & J0,               &
      & J1

    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_assimilation_waterunlimited'

    ! Prepare CO2_conc_leaf_cl
    ! R: siehe mo_bethy.f90, Zeile 709
    IF (C4Flag) THEN ! landcover type is C4 plants
      CO2_conc_leaf_cl = FCI1C4 * CO2_mol_mixing_ratio
    ELSE       ! landcover type is C3 plants
      CO2_conc_leaf_cl = FCI1C3 * CO2_mol_mixing_ratio
    END IF

    T1 =  25._wp + tmelt   ! 25 degress Celsius expressed in Kelvin

    TC = t_air - tmelt     ! change temperature dimension from Kelvin to Celcius

    T0 = t_air - T1        ! Atmosphere temperature (lowest layer) in Kelvin minus
                           ! 25 degress Celsius expressed in Kelvin
                           ! => relative T to 25 degree Celcius

    !---------------------------------------------------------------------------------
    !                  FIRST ENTRY, TC=TA, CI0 => GC0, AC0
    !---------------------------------------------------------------------------------

    IF (.NOT. C4Flag) THEN

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! CASE: No water limitation and C3 !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---------------------------------------------------------------------------------
      !   C3     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT
      !
      ! Rate (with Aktivationenergy) vegetation temperature dependance is
      !  k = k(25C) * EXP((Veg Temp - 25) * aktivation energy
      !                    / 298.16 * R * (Veg Temp + 273.16))
      ! => k = k0 * EXP( T0 * E / T1 / R / t_air ),  WHERE R is the gas constant (8.314)
      ! This holds for OX-Oxygen partial pressure, KC-Michaelis-Menten constant for CO2,
      !    KO-Michaelis-Menten constant for O2, PVM-carboxylation capacity,
      !    RD-Dark respiration, K-PEPcase CO2 specivity
      !    Knorr (106)
      !---------------------------------------------------------------------------------
      KC = KC0 * EXP(EC / R_gas * T0 * T1**(-1) * t_air**(-1))  ! Division "/T1 /t_air" are expressed strange here (and in the
                                                                ! following) because if this subroutine is not used as elemental
                                                                ! (e.g. for tests) it makes accuracy of calculation much better!
                                                                ! R: Abweichungen des Ergebnisses im Vergl. zu JSBACH3 ab 14. 
                                                                !    Stelle nach dem Komma, obwohl Inputs (R_gas, T0, T1, t_air)
                                                                !    genau gleich. Ursache: epsilon ist 2.22E-16 (JSBACH3
                                                                !    und JSBACH4) und EXP verstärkt Fehler.
      KO = KO0 * EXP(EO / R_gas * T0 * T1**(-1)  * t_air**(-1)) ! see above

      !---------------------------------------------------------------------------------!
      ! CO2 compensation point without leaf respiration, Gamma* is assumed to be linearly
      ! dependant on vegetation temperature, Gamma* = 1.7 * TC (IF Gamma* in microMol/Mol)
      ! Here, Gam in Mol/Mol,       Knorr (105)
      !---------------------------------------------------------------------------------
      GAM = MAX (1.7E-6_wp * TC, 0._wp)

      !---------------------------------------------------------------------------------
      ! PVM and PJM are not only temperature dependant but also differ inside the canopy.
      ! This is due to the fact that the plant distributes its Nitrogen content and
      ! therefore Rubisco content so that, the place with the most incoming light got the
      ! most Rubisco. Therefore, it is assumed that the Rubisco content falls
      ! exponentially inside the canopy. This is reflected directly in the values of PVM
      ! and PJM at 25 Celsius (PVM * nscl),  Knorr (107/108)
      !---------------------------------------------------------------------------------
      VCMAX =   CarboxRate * N_scaling_factors &
            & * EXP(EV / R_gas * T0 * T1**(-1) * t_air**(-1))

      !---------------------------------------------------------------------------------
      ! The temperature dependance of the electron transport capacity follows
      ! Farqhuar(1988) with a linear temperature dependance according to the vegetation
      ! temperature
      !  J = J(25C) * TC / 25 WHERE J(25) = J0 * N_scaling_factors
      ! ?????????????????????/ minOfMaxCarboxrate=1E-12 vielleicht etwas gering
      !---------------------------------------------------------------------------------
      JMAX = ETransport * N_scaling_factors * TC/25._wp
      JMAX = MAX(JMAX,minOfMaxCarboxrate) ! R: Abweichungen des Ergebnisses im Vergl. zu JSBACH3 ab 21. Stelle nach
                                          !    Komma. Ursache: epsilon ist 2.22E-16 (JSBACH3 und JSBACH4).

      !---------------------------------------------------------------------------------
      !
      !  C3       GROSS PHOTOSYNTHESIS AT GIVEN CI
      !
      !---------------------------------------------------------------------------------
      !c  The assimilation follows the Farqhuar (1980) formulation for C3 plants
      !  A = min{JC, JE} - RD
      !  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
      !  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)      with
      !   J = alpha * I * PJM / sqrt(PJM^2 + alpha^2 * I^2) with I=PAR in Mol(Photons)
      !        Knorr (102a-c, 103)
      !  Here J = J1 and A is the gross photosynthesis (GASS), i.e. still including the
      !          respiratory part RD
      !---------------------------------------------------------------------------------
      IF (JMAX .GT. minOfMaxCarboxrate) THEN
        J1 =   ALPHA * apar_per_lai_cl * JMAX &
          &  / SQRT(JMAX**2 + (ALPHA * apar_per_lai_cl)**2)  ! R: Abweichungen des Ergebnisses im Verg. zu JSBACH3 ab 15./16. Stelle
      ELSE                                                   !    nach Komma obwohl ALPHA, apar_per_lai_cl, JMAX genau gleich.
        J1 = 0._wp
      END IF

      JE = J1 * (CO2_conc_leaf_cl - GAM) / 4._wp / (CO2_conc_leaf_cl + 2._wp * GAM) ! R: keine Abweichung wenn Input gleich

      JC =   VCMAX * (CO2_conc_leaf_cl - GAM) &
        &  / ( CO2_conc_leaf_cl + KC * (1._wp + OX * KO**(-1)) ) ! R: obwohl CO2_conc_leaf_cl, GAM, KC, OX, KO genau gleich
                                                                 !    Abweichung ab 15. Stelle nach Komma. Wenn man OX-Parameter
                                                                 !    mit 0.21_wp ersetzt, dann ist die Abweichung gleich Null.
                                                                 !    Grund: epsilon Abweichung entsteht durch das usen von OX
                                                                 !    auch wenn dasselbe OX verwendet wird. Wahrscheinlich nimmt
                                                                 !    Kompiler einen Kopiervorgang vor.

      HITINHIB = 1._wp / ( 1._wp + EXP( 1.3_wp * ( TC - 55._wp ) ) ) ! R: bei gleichem Input keine Abweichung (22 Stell nach Komma)
      GASS = MIN (JE, JC) * HITINHIB                                 ! R: Bei gleichem Input keine Abweichung (22 Stell nach Komma)

      !---------------------------------------------------------------------------------
      !
      !  C3      COMPUTE 'DARK' RESPIRATION, PHOTORESPIRATION, STOMATAL CONDUCTANCE
      !
      !---------------------------------------------------------------------------------
      !---------------------------------------------------------------------------------
      ! Following Farqhuar et al. (1980), the dark respiration at 25C is proportional to
      ! PVM at 25C, therefore RD = const * PVM, but the temperature dependance goes with
      ! ER (for respiration) and not with EV (for PVM)
      !---------------------------------------------------------------------------------
      DARKINHIB = 0.5_wp + 0.5_wp * EXP(-2.e5_wp * MAX(par_down_mol,0.0_wp)) ! R: bei gleichem Input keine Abweichung (22 Stellen
                                                                             !    nach Komma)

      DARK_RESP =   FRDC3 * CarboxRate * N_scaling_factors         &
                & * EXP(ER * T0 / R_gas  * T1**(-1) * t_air**(-1)) &
                & * HITINHIB * DARKINHIB                          ! R: auch wenn FRDC3, CarboxRate, N_scaling_factors, ER, T0,
                                                                  !    T1, R_gas, t_air, HITINHIB, DARKINHIB gleich sind
                                                                  !    entstehen ab der 16. Stelle nach Komma Abweichungen wg. EXP.

      !---------------------------------------------------------------------------------
      ! Diffusion equation Flux = (CA - CI) / resistence, rs
      !   CI=CO2_conc_leaf
      !   conductance canopy_conductance = 1 / rs  =>  Flux = (CA-CI) * canopy_conductance
      !   (CA ... CO2mixingRatio)
      !   Flux of CO2 is A * amount, Assimilation rate * amount
      !   A is here, Gross Assimilation, though A-RD = (net) Assimilation rate
      !   the amount comes from the ideal gas equation pV=nRT => n/V = p / RT
      !   the stomatal conductance for CO2 is less THEN the conductance of H2O by
      !   the factor of 1.6: canopy_conductance(CO2) = canopy_conductance(H2O) / 1.6, due to its lower mobiblity due
      !   to its higher mass
      !   => A (net) = canopy_conductance/1.6 * (CA-CI) * p/RT
      !   => canopy_conductance = A(net)*1.6*RT/p/(CA-CI)
      !---------------------------------------------------------------------------------
      canopy_cond_cl = MAX ( &
        & 1.6_wp * (GASS - DARK_RESP) / (CO2_mol_mixing_ratio - CO2_conc_leaf_cl) * R_gas * t_air * press_srf**(-1), &
        & minStomaConductance)          ! R: bei gleichem Input keine Abweichung (22 Stellen nach Komma).
                                        !    Mit Abweichung die von DARK_RESP kommt ist die Maximal-Abweichung dann schon
                                        !    ab 12. Stelle nach dem Komma.

    ELSE IF (C4Flag) THEN

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! CASE: No water limitation and C4 !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !---------------------------------------------------------------------------------
      !
      !   C4     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT
      !           AND 'DARK' RESPIRATION
      !
      !---------------------------------------------------------------------------------
      ! For C4 plants the Farquhar equations are replaced by the set of equations of
      !  Collatz et al. 1992:
      !  A = min{JC, JE} - RD
      !  JC = k * CI
      !  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]      with
      !  Ji = alphai * Ipar / Epar with Ji=PAR in Mol(Photons)
      !        Knorr (114a-d)
      !  alphai is the integrated quantum efficiency for C4 plants (ALC4 = 0.04,
      !    compared to the efficiency of C3 plants, ALPHA = 0.28)
      !  Theta is the curve PARAMETER (0.83) which makes the change between
      !   PVM and K limitation smooth
      !  K is the PECase CO2 specifity instead of the electron transport capacity
      !   within C3 plants
      !  Ci is the stomatal CO2 concentration = Cimin + (Ci0 - Cimin)* GC/GC0 with
      !    Cimin = 0.3 and 0.15 CA respectivly (CA is the CO2 mixing ratio)
      !
      ! The factor 1E3 comes that PJM for C3 is in microMol and K is in milliMol,
      !   which is not considered in INITVEGDATA
      ! K scales of course with EK
      !---------------------------------------------------------------------------------
      K =   ETransport * 1.E3_wp * N_scaling_factors &
        & * EXP(EK * T0 / R_gas * T1**(-1) * t_air**(-1)) ! R: Abweichung ab der 16. Stelle nach dem Komma also im Bereich EPSILON

      !---------------------------------------------------------------------------------
      ! same as C3
      !---------------------------------------------------------------------------------
      VCMAX =   CarboxRate * N_scaling_factors &
            & * EXP(EV / R_gas * T0 * T1**(-1) * t_air**(-1))  ! R: Abw. ab der 19. Stelle nach dem Komma also im Bereich EPSILON

      !---------------------------------------------------------------------------------
      !  same as C3, just the 25 degree Celsius proportional factor is different
      !    0.011 for C3,  0.0042 for C4
      !---------------------------------------------------------------------------------
      DARKINHIB = 0.5_wp + 0.5_wp * EXP(-2.e5_wp * MAX(par_down_mol,0.0_wp)) ! R: keine Abweichung
      HITINHIB  = 1._wp / ( 1._wp + EXP( 1.3_wp * ( TC - 55._wp ) ) )        ! R: keine Abweichung
      DARK_RESP =    FRDC4 * CarboxRate * N_scaling_factors &
                & *  EXP(ER * T0 / R_gas * T1**(-1) * t_air**(-1)) * HITINHIB * DARKINHIB ! R: Abweichung ab der 15. Stelle nach
                                                                                          ! dem Komma also im Bereich EPSILON

      !---------------------------------------------------------------------------------
      !
      !  C4       GROSS PHOTOSYNTHESIS AT GIVEN CI                                     C
      !
      !---------------------------------------------------------------------------------
      !  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]
      !    Ji = ALC4 * PAR
      !  J0 is the sum of the first two terms in JE
      !---------------------------------------------------------------------------------
      J0 = (ALC4 * apar_per_lai_cl + VCMAX) /  2._wp / THETA
      !---------------------------------------------------------------------------------
      !  last 2 terms:  with J0^2 = 1/4/Theta^2*(PVM+Ji)^2
      !       sqrt(1/4/Theta^2)*sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji))
      !   = sqrt (J0^2 - PVM*Ji/Theta)
      !---------------------------------------------------------------------------------
      JE =   J0 - SQRT (J0**2 - VCMAX &
        &  * ALC4 * apar_per_lai_cl / THETA) ! R: Abweichung ab der 13. Stelle nach dem Komma
      !---------------------------------------------------------------------------------
      !         see above
      !---------------------------------------------------------------------------------
      JC = K * CO2_conc_leaf_cl ! R: Abweichung ab der 20. Stelle nach dem Komma also im Bereich EPSILON
      !---------------------------------------------------------------------------------
      ! same as C3, Farquhar Assimilation
      !---------------------------------------------------------------------------------
      GASS = MIN(JE, JC) * HITINHIB ! R: Abweichung ab der 14. Stelle nach dem Komma also im Bereich EPSILON
      !---------------------------------------------------------------------------------
      !
      !   C4     COMPUTE PHOTORESPIRATION, STOMATAL CONDUCTANCE
      !
      !---------------------------------------------------------------------------------
      ! same as C3 (diffusion equation)
      !---------------------------------------------------------------------------------
      canopy_cond_cl =   1.6_wp * (GASS - DARK_RESP) / (CO2_mol_mixing_ratio &
        &              - CO2_conc_leaf_cl) * R_gas * t_air * press_srf**(-1)
      canopy_cond_cl = MAX(canopy_cond_cl, minStomaConductance) ! R: Abweichung ab der 10. Stelle nach dem Komma

    END IF  ! C3/C4

  END FUNCTION calc_assimilation_waterunlimited

  SUBROUTINE calc_assimilation_waterlimited( &
    & C4Flag,                                          & ! in
    & CarboxRate,                                      & ! in
    & ETransport,                                      & ! in
    & t_air,                                           & ! in
    & press_srf,                                       & ! in
    & par_down_mol,                                    & ! in
    & apar_per_lai_cl,                                 & ! in
    & CO2_mol_mixing_ratio,                            & ! in
    & N_scaling_factors,                               & ! in
    & CO2_conc_leaf_cl,                                & ! out
    & canopy_cond_cl,                                  & ! in
    & gross_assimilation_cl,                           & ! out
    & dark_respiration_cl,                             & ! out
    & carbox_rate_max_cl,                              & ! out
    & e_transport_rate_max_cl,                         & ! out
    & carbox_rate_cl,                                  & ! out
    & e_transport_rate_cl                              & ! out
    & )


    !-----------------------------------------------------------------------
    !  DECLARATIONS
    USE mo_assimi_constants                                      ! we need all of them, except of f_aut_leaf and cCost.
    USE mo_physical_constants, ONLY: tmelt, R_gas => argas       ! tmelt =273.15_wp;  argas =Universal gas constant [J/(mole*K)]

    !-----------------------------------------------------------------------
    !  ARGUMENTS
    LOGICAL,  INTENT(in) ::      &
      & C4Flag
    REAL(wp), INTENT(in) ::      &
      & CarboxRate,              &
      & ETransport

    REAL(wp), INTENT(in), DIMENSION(:) ::      &
      & t_air,                   & ! Atmosphere temperature (lowest layer) in Kelvin
      & press_srf,               & ! Surface pressure
      & par_down_mol,            & ! Downward PAR flux in mol (photons)/(m^2 s)
      & apar_per_lai_cl,         & ! Absorbed PAR of canopy layer.
      & CO2_mol_mixing_ratio,    & ! CO2 substance mixing ratio [molCO2/molDryAir]
      & N_scaling_factors,       & ! Enzymatic scaling with canopy height
      & canopy_cond_cl             ! Input only for the water limited case here

    REAL(wp), INTENT(out), DIMENSION(:) ::     &
      & CO2_conc_leaf_cl,        & ! CO2 concentration inside leaf [mol(CO2)/mol(Air)] (identical for all layers):
                                   ! Output only for the water limited case here
      & gross_assimilation_cl,   & ! Gross assimilation of C per leaf area [ mol (CO2) / (m^2(leaf area) * s) ] and canopy layer
      & dark_respiration_cl,     & ! Dark respiration per leaf area and per canopy layer
      & carbox_rate_max_cl,      & ! Maximum carboxylation rate of canopy layer (=VCmax)
                                   ! [mol(CO2)/m^2(leaf area)/s]. carbox_rate_max_leaf adapted to canopy layer
                                   ! and actual temperature.
      & e_transport_rate_max_cl, & !
      & carbox_rate_cl,          & !
      & e_transport_rate_cl        !

    !-----------------------------------------------------------------------
    ! Local variables
! !     INTEGER                           :: icanopy
! !     INTEGER                           :: i

    REAL(wp) ::    &
      & T0,        & ! Temperature relative to 25 degree Celcius, i.e. T - 25
      & T1,        & ! 25 degress Celsius expressed in Kelvin
      & KC,        & ! Michaelis-Menten constant for CO2
      & KO,        & ! Michaelis-Menten constant for O2
      & GAM,       & ! COMPENSATION POINT WITHOUT DARK RESPIRATION [MOL(CO2) / MOL(AIR)]
      & HITINHIB,  & ! Hight temperature inhibition
                     ! Value to inhibit Assimilation and Respiration at temperatures above
                     ! 55 Celsius from Collatz et al., Physiological and environmental regulation
                     ! of stomatal conductance, photosynthesis and transpiration: a model that
                     ! includes a laminar boundary layer,
                     ! Agricultural and Forest Meteorology, 54, pp. 107-136, 1991
      & DARKINHIB, & ! Dark inhibition
                     ! Value to inhibit Dark-Respiration in light
                     ! after Brooks and Farquhar, Effect of temperature on the CO2/O2 specifity on RBisCO
                     ! and the rate of respiration in the light, Planta 165, 397-406, 1985
                     ! It is fitted to inhibit the dark-respiration to 50% of it's uninhibited value
                     ! up from 50 umol/m^2s total solar downward radiation, i.e. 25 umol/m^2s downward PAR
                     ! (idea for improvement: darkinhib should be closer to 1 at very high latitudes)
                     ! The FUNCTION is an own creation.
      & K1,        & ! helper
      & K2           ! helper

    REAL(wp) ::    & ! variables with canopy layer
      & VCMAX,     & ! MAXIMUM CARBOXYLATION RATE [MOL(CO2)/M^2 S]
      & JMAX,      & ! MAXIMUM RATE OF ELECTRON TRANSPORT [MOL(CO2)/M^2 S]
      & K,         & ! CO2 SPECIFICITY OF PECASE
      & TC,        & ! Canopy Temperature (leaf) in Celsius.
      & GASS,      & ! gross assimilation (formerly called A)
      & DARK_RESP, & ! dark respiration (formerly called RD).
      & JE,        & ! Electron transport rate (Knorr, Eq. 102c)
      & JC,        & ! Carboxylation rate (Knorr, Eq. 102b)
      & J0,        &
      & J1,        &
      & G0,        &
      & W1,        & ! helper
      & W2,        & ! helper
      & B,         & ! helper
      & C            ! helper

    INTEGER :: nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_assimilation_waterlimited'

    !-----------------------------------------------------------------------
    ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
    ! IF (one_of('calc_assimilation', tile%components) < 1) RETURN

    nc = SIZE(apar_per_lai_cl)

    !$ACC DATA &
    !$ACC   PRESENT(t_air, press_srf, par_down_mol, apar_per_lai_cl, co2_mol_mixing_ratio, n_scaling_factors) &
    !$ACC   PRESENT(canopy_cond_cl, co2_conc_leaf_cl, gross_assimilation_cl, dark_respiration_cl, carbox_rate_max_cl) &
    !$ACC   PRESENT(e_transport_rate_max_cl, carbox_rate_cl, e_transport_rate_cl)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR

    DO ic=1,nc

      T1 =  25._wp + tmelt   ! 25 degress Celsius expressed in Kelvin

      TC = t_air(ic) - tmelt     ! change temperature dimension from Kelvin to Celcius

      T0 = t_air(ic) - T1        ! Atmosphere temperature (lowest layer) in Kelvin minus
                            ! 25 degress Celsius expressed in Kelvin
                            ! => relative T to 25 degree Celcius

      !---------------------------------------------------------------------------------
      !                  FIRST ENTRY, TC=TA, CI0 => GC0, AC0
      !---------------------------------------------------------------------------------

      IF (.NOT. C4Flag) THEN
        !WHERE (.not. C4flag) ! R: only for test issues when this subroutine is not elemental

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! CASE: water limitation and C3 !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !---------------------------------------------------------------------------------
        !   C3     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT
        !           And 'DARK' RESPIRATION
        !---------------------------------------------------------------------------------
        KC = KC0 * EXP(EC / R_gas * T0 * T1**(-1) * t_air(ic)**(-1) )  ! R: Abweichung ab der 13. Stelle nach dem Komma,
                                                                  !    obwohl T0 und T1 und t_air genau gleich. Ursache:
                                                                  !    epsilon ist 2.22E-16 (JSBACH3 und JSBACH4) und EXP verstärkt
                                                                  !    Fehler.
        KO = KO0 * EXP(EO / R_gas * T0 * T1**(-1) * t_air(ic)**(-1) )  ! R: Abweichung ab der 9. Stelle nach dem Komma, obwohl T0 und T1
                                                                  !    und t_air genau gleich. (Grund warum hier Abweichung größer
                                                                  !    als oben: KO0 ist um Faktor 3 größer als KC0.)
        GAM = MAX (1.7E-6_wp * TC, 0._wp)                          ! R: Keine Abweichungen
        VCMAX = CarboxRate * N_scaling_factors(ic) &
                  * EXP(EV * T0 / R_gas * T1**(-1) * t_air(ic)**(-1) ) ! PVM=PVM0*N_scaling_factors*EXP(Energy scaled relative 25C)
                                                                  ! R: Abweichung ab der 17. Stelle nach dem Komma

        DARKINHIB = 0.5_wp + 0.5_wp * EXP(-2.e5_wp * MAX(par_down_mol(ic),0.0_wp)) ! R: keine Abweichungen
        HITINHIB  = 1._wp / ( 1._wp + EXP( 1.3_wp * ( TC - 55._wp ) ) )        ! R: keine Abweichungen
        DARK_RESP = FRDC3 * CarboxRate * N_scaling_factors(ic)          &     ! RD=RD0*N_scaling_factors*EXP(Energy scaled relative 25C)
                  * EXP(ER * T0  / R_gas * T1**(-1) * t_air(ic)**(-1) ) &            ! with RD0 proportional PVM0 at 25C
                  * HITINHIB * DARKINHIB                                         ! R: Abweichung ab der 15. Stelle nach dem Komma
        JMAX = ETransport * N_scaling_factors(ic)       &  ! PJM=PJM0*N_scaling_factors*EXP(Energy scaled relative 25C)
                  * TC / 25._wp                        ! R: bis hierher keinerlei Abweichung. N_scaling_factors und ETransport genau
                                                      !    gleich.
        JMAX = MAX(JMAX, minOfMaxCarboxrate) ! R: Abweichung nur für Fall minOfmaxCarboxrate und dann ab der 13. Stelle nach dem
                                            !    Komma weil minOfMaxCarboxrate = 1.0e-12 ist

        !---------------------------------------------------------------------------------
        !
        !  C3       GROSS PHOTOSYNTHESIS AT GIVEN TC
        !
        !---------------------------------------------------------------------------------
        !  Remember:
        !  A = min{JC, JE} - RD
        !  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
        !  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)      with
        !   J = alpha * I * PJM / sqrt(PJM^2 + alpha^2 * I^2) with I=PAR in Mol(Photons)
        !        Knorr (102a-c, 103)
        ! J = J1
        !---------------------------------------------------------------------------------
        IF (JMAX .GT. minOfMaxCarboxrate) THEN
          J1 = ALPHA * apar_per_lai_cl(ic) * JMAX &
                    / SQRT(JMAX**2 + (ALPHA * apar_per_lai_cl(ic))**2) ! R: Abweichung ab der 15. Stelle nach dem Komma.
        ELSE
          J1 = 0._wp
        END IF
        !---------------------------------------------------------------------------------
        !         Helping friends K1, W1, W2, K2
        !---------------------------------------------------------------------------------
        K1 = 2._wp * GAM
        W1 = J1 / 4._wp
        W2 = VCMAX
        K2 = KC * (1._wp + OX / KO)
        !---------------------------------------------------------------------------------
        ! A = canopy_conductance / 1.6 * (CA - Ci) * p / R / t_air
        ! <=> Ci = CA - 1.6 * R * t_air / p / canopy_conductance * A = CA - A / G0
        !---------------------------------------------------------------------------------
        !if (R_gas < 0.1 .OR. t_air < 0.1 .OR. press_srf < 0.1) &
        !print*,iGP,R_gas,t_air,press_srf
        G0 = canopy_cond_cl(ic) / 1.6_wp / R_gas * t_air(ic)**(-1) * press_srf(ic) ! R: Die folgenden Berechnungen führen zu keinen Abweichungen
                                                                      !    im Ergebnis CO2_conc_leaf_cl und bei GASS zu einer
                                                                      !    Abweichung ab der 16. Stelle nach Komma
        !---------------------------------------------------------------------------------
        ! A = min{JC, JE} - RD
        ! => A = JC - RD ^ A = JE - RD
        ! Set this (A =) in Ci formula above
        ! Set Ci in
        !  JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)
        ! => quadratic formula in JE
        ! 0 = JE^2 -(RD+J/4+G0*(CA+2*GAMMA))*JE +J/4*G0*(CA-GAMMA)+J/4*RD
        !---------------------------------------------------------------------------------
        B = DARK_RESP + W1 + G0 * (CO2_mol_mixing_ratio(ic) + K1)
        C = W1 * G0 * (CO2_mol_mixing_ratio(ic) - GAM) + W1 * DARK_RESP
        !---------------------------------------------------------------------------------
        ! with 0 = x^2 - bx + c
        !      x1/2 = b/2 +/- sqrt(b^2/4 - c)
        !  take '-' as minimum value of formula
        !---------------------------------------------------------------------------------
        IF (JMAX .GT. minOfMaxCarboxrate) THEN
          JE = B / 2._wp - SQRT ( MAX (B**2 / 4._wp - C, 0._wp))
        ELSE
          JE = 0._wp
        END IF
        !---------------------------------------------------------------------------------
        ! Set Ci in
        !  JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
        ! WRITE JC = PVM * (Ci - Gam) / (Ci + K2)
        ! => quadratic formula in JC
        ! 0 = JC^2 -(RD+PVM+G0*(CA+K2))*JC +PVM*G0*(CA-GAMMA)+RD*PVM
        !---------------------------------------------------------------------------------
        B = DARK_RESP + W2 + G0 * (CO2_mol_mixing_ratio(ic) + K2)
        C = W2 * G0 * (CO2_mol_mixing_ratio(ic) - GAM) + W2 * DARK_RESP
        JC = B / 2._wp - SQRT ( MAX(B**2 / 4._wp - C, 0._wp))
        !---------------------------------------------------------------------------------
        ! A = min{JC, JE} - RD
        !  but here A = Gross Photosynthesis = iGPP
        !  => A = min{JC, JE}
        !---------------------------------------------------------------------------------
        GASS = MIN(JE, JC) * HITINHIB ! R: Abweichung ab der 16. Stelle nach Komma
        !---------------------------------------------------------------------------------
        !
        !   C3     COMPUTE PHOTORESPIRATION AND INTER STOMATE CO2 CONCENTRATION
        !
        !---------------------------------------------------------------------------------
        ! A = canopy_conductance / 1.6 * (CA - Ci) * p / R / t_air
        ! <=> Ci = CA - 1.6 * R * t_air / p / canopy_conductance * A = CA - A / G0
        !   with A = Net assimilation = NPP
        ! (CA is the CO2 mixing ratio)
        !---------------------------------------------------------------------------------
        ! R: Note: CO2_conc_leaf acts as output only for water-limited case. Then only the
        !   result for the actual canopy layer goes back to update_bethy as bethy%CO2_conc_leaf has no canopy layers.
        !   Therefore bethy%CO2_conc_leaf_acc is calculated inside the canopy loop!
        CO2_conc_leaf_cl(ic) = CO2_mol_mixing_ratio(ic) - &
          & MAX((GASS - DARK_RESP) / MAX(G0, 1.E-6_wp), 0._wp) ! R: keine Abweichung auch nicht aus der Summe aus Abweichungen oben
        !---------------------------------------------------------------------------------
        ! Photorespiration
        ! Carboxylilation controlled assimilation:
        ! JC = PVM * (Ci - Gam) / (Ci + KC * (1 + OX/KO))
        ! Photorespiration = PVM * Gam / (Ci + KC * (1 + OX/KO))
        ! Light limited Assimilation:
        ! JE = J * (Ci - Gam) / 4 / (Ci + 2 * Gam)
        ! Photorespiration = J * Gam / 4 / (Ci + 2 * Gam)
  !!              PHOTO_RESP = MIN(VCMAX * GAM / ( CO2_conc_leaf_cl + KC * ( 1._wp + OX / KO )), &
  !!                   J1 * GAM / 4._wp / (CO2_conc_leaf_cl + 2._wp * GAM))                      &
  !!                * HITINHIB

      ELSE

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! CASE: water limitation and C4 !!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !---------------------------------------------------------------------------------
        !
        !                               C 4
        !
        !
        !   C4     DETERMINE TEMPERATURE-DEPENDENT RATES AND COMPENSATION POINT
        !           AND 'DARK' RESPIRATION
        !
        !---------------------------------------------------------------------------------
        ! Remind:
        !  Collatz et al. 1992:
        !  A = min{JC, JE} - RD
        !  JC = k * Ci
        !  JE = 1/2/Theta *[PVM + Ji - sqrt((PVM+Ji)^2 - 4*Theta*PVM*Ji)]      with
        !   Ji = alphai * Ipar / Epar with Ji=PAR in Mol(Photons)           and
        !   PAR = Ipar / Epar;  ALC4=alphai; J0=1/2/Theta *(PVM + Ji);
        !   Ci = CA - 1.6 * R * t_air / p / canopy_conductance * A = CA - A / G0                and
        !   => A = JC - RD ^ A = JE - RD
        !---------------------------------------------------------------------------------

        K = ETransport * 1.E3_wp * N_scaling_factors(ic)          &      ! K=K0*N_scaling_factors*EXP(Energie scaled relative 25C)
                * EXP(EK * T0 / R_gas * T1**(-1) * t_air(ic)**(-1) )     ! in mmol CO2 Pecase for C4.
        VCMAX = CarboxRate * N_scaling_factors(ic) &                     ! PVM=PVM0*N_scaling_factors*EXP(Energie
                  * EXP(EV * T0 / R_gas * T1**(-1) * t_air(ic)**(-1) )  ! R: Abweichung an der 22. Stelle.

        DARKINHIB = 0.5_wp + 0.5_wp * EXP(-2.e5_wp * MAX(par_down_mol(ic),0.0_wp)) ! R: keine Abweichungen
        HITINHIB  = 1._wp / ( 1._wp + EXP( 1.3_wp * ( TC - 55._wp ) ) )        ! R: keine Abweichungen
        DARK_RESP = FRDC4 * CarboxRate * N_scaling_factors(ic)  &           ! DARK_RESP=RD0*N_scaling_factors*EXP(Energie scaled relative
                  * EXP(ER * T0 / R_gas * T1**(-1) * t_air(ic)**(-1) )   & !  25C) with RD0 proportional PVM0 at 25C
                  * HITINHIB * DARKINHIB                               ! Abweichung ab der 16. Stelle nach dem Komma
        !---------------------------------------------------------------------------------
        !
        !  C4       GROSS PHOTOSYNTHESIS AT GIVEN CI (=CO2_conc_leaf)
        !
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !  Ci = CA - 1.6 * R * t_air / p / canopy_conductance * A = CA - A / G0
        !---------------------------------------------------------------------------------
        G0 = canopy_cond_cl(ic) / 1.6_wp / R_gas * t_air(ic)**(-1) * press_srf(ic)
        !---------------------------------------------------------------------------------
        !  J0=1/2/Theta *(PVM + Ji) = (alphai * PAR + PVM) / 2 / Theta
        !---------------------------------------------------------------------------------
        J0 = (ALC4 * apar_per_lai_cl(ic) + VCMAX) /  2._wp / THETA
        !---------------------------------------------------------------------------------
        !  JE = J0 - sqrt( J0^2 - PVM*alphai*PAR/Theta)
        !---------------------------------------------------------------------------------
        JE = J0 - SQRT (J0**2 - VCMAX * ALC4 * apar_per_lai_cl(ic) / THETA) ! R: Abweichung ab der 13. Stelle.
        !---------------------------------------------------------------------------------
        !  JC = (G0*CA + RD) / (1 + G0/K)
        !---------------------------------------------------------------------------------
        JC = (G0 * CO2_mol_mixing_ratio(ic) + DARK_RESP) / (1._wp + G0 / K) ! R: Abweichung ab der 16. Stelle.
        !---------------------------------------------------------------------------------
        !  A = min{JC, JE} - RD, but here: A = iGPP
        !---------------------------------------------------------------------------------
        GASS = MIN(JE, JC) * HITINHIB ! R: Abweichung ab der 14. Stelle (?sollte aber 13. sein wg. JE. ?)
        !---------------------------------------------------------------------------------
        !   C4     COMPUTE PHOTORESPIRATION AND INTER STOMATE CO2 CONCENTRATION
        !---------------------------------------------------------------------------------
        !  Ci = CA - 1.6 * R * t_air / p / canopy_conductance * A = CA - A / G0, with A = NPP
        !---------------------------------------------------------------------------------
        ! R: Note: CO2_conc_leaf acts as output only for water-limited case. Then only the
        !    result for the actual canopy layer goes back to update_bethy as bethy%CO2_conc_leaf has no canopy layers.
        !    Therefore bethy%CO2_conc_leaf_acc is calculated inside the canopy loop!
        CO2_conc_leaf_cl(ic) = CO2_mol_mixing_ratio(ic) - &
                MAX((GASS-DARK_RESP) / MAX(G0, 1.E-6_wp), 0._wp)  ! R: Keine Abweichungen
        !---------------------------------------------------------------------------------
        ! No Photorespiration with C4 plants
        !---------------------------------------------------------------------------------
        ! PHOTO_RESP = 0._wp
        JMAX = 0._wp
      END IF !  C3/C4

      !------------------------------------------------------------------------------------------------------------------
      ! Diagnostic output for water-limited case
      !------------------------------------------------------------------------------------------------------------------
      gross_assimilation_cl(ic)   = GASS             ! gross assimilation
      dark_respiration_cl(ic)     = DARK_RESP        ! dark respiration
      ! photo_respiration(ic)      = PHOTO_RESP        ! photo respiration
      carbox_rate_max_cl(ic)      = VCMAX            ! max carboxylation rate
      e_transport_rate_max_cl(ic) = JMAX             ! max e-transport rate
      carbox_rate_cl(ic)          = JC               ! carboxylation rate
      e_transport_rate_cl(ic)     = JE               ! e-transport rate

    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE calc_assimilation_waterlimited


  SUBROUTINE calc_NPP_pot_rate( &
    & gross_assimilation_ca,              &
    & dark_respiration_ca,                &
    & NPP_pot_rate_ca )

    ! Use Declarations
    USE mo_assimi_constants,    ONLY: f_aut_leaf, cCost

    ! Arguments
    REAL(wp),  INTENT(in), DIMENSION(:) ::   &
      & gross_assimilation_ca, & ! Gross primary product (where photorespiration has already  ..
                                 ! .. been accounted for) [mol(CO2)/(m^2 (canopy ground) s)]
      & dark_respiration_ca      ! Dark respiration of leaves [mol(CO2)/(m^2 (canopy ground) s)]

    REAL(wp),  INTENT(out), DIMENSION(:) :: &
      & NPP_pot_rate_ca            ! Net primary production rate [mol(CO2)/(m^2 (canopy ground) s)]

    ! Local variables
    CHARACTER(len=*), PARAMETER :: routine = modname//':calc_NPP_pot_rate'

    REAL(wp),PARAMETER   :: hlp =  (cCost -1._wp) / cCost
    REAL(wp)             :: maintenanceRespiration ! R_m
    REAL(wp)             :: growthRespiration      ! R_g
    INTEGER :: nc, ic

    !-----------------------------------------------------------------------
    ! DOES THIS PROCESS EXIST FOR THE CURRENT TILE?
    ! IF (one_of('calc_assimilation', tile%components) < 1) RETURN

    ! ---------------------------
    ! Go

    nc = SIZE(gross_assimilation_ca)

    !$ACC DATA PRESENT(gross_assimilation_ca, dark_respiration_ca, npp_pot_rate_ca)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR &
    !$ACC   PRIVATE(growthRespiration, maintenanceRespiration)
    DO ic=1,nc
      maintenanceRespiration = dark_respiration_ca(ic) / f_aut_leaf
      growthRespiration = max(0._wp,hlp * (gross_assimilation_ca(ic) - maintenanceRespiration))
      NPP_pot_rate_ca(ic) = gross_assimilation_ca(ic) - maintenanceRespiration - growthRespiration
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE calc_NPP_pot_rate

  !
  ! ===============================================================================================================================
  !>
  !>   Compute running average of NPP for NLCC process
  !!

  PURE ELEMENTAL SUBROUTINE calc_NPPbuf (lstart, init_running_means, new_year, delta_time, NPP, &
                                         seconds_year, NPP_sum_year, NPP_mean_5year)

    LOGICAL,  INTENT(in)    :: lstart             ! first time step in a simulation
    LOGICAL,  INTENT(in)    :: init_running_means ! initialize the min/max temperature climatology
    LOGICAL,  INTENT(in)    :: new_year           ! first time step in a year
    REAL(wp), INTENT(in)    :: delta_time         ! time step length [s]
    REAL(wp), INTENT(in)    :: NPP                ! potential NPP [mol(C)m-2(canopy)s-1]
    REAL(wp), INTENT(inout) :: seconds_year       ! sum of seconds till new_year
    REAL(wp), INTENT(inout) :: NPP_sum_year       ! 
    REAL(wp), INTENT(inout) :: NPP_mean_5year     ! NPP

!   local variables
    REAL(wp) :: prev_year_mean_NPP ! NPP of previous year

   ! update NPP_sum_year and NPP_mean_5year
   IF (.NOT. new_year .OR. lstart) THEN
      NPP_sum_year = NPP_sum_year + NPP * delta_time
   ELSE
      prev_year_mean_NPP = NPP_sum_year / seconds_year
      NPP_sum_year = NPP * delta_time
      IF (init_running_means) THEN
         NPP_mean_5year = prev_year_mean_NPP
      ELSE
         NPP_mean_5year = (NPP_mean_5year * 4._wp + prev_year_mean_NPP) / 5._wp
      END IF
   ENDIF

   ! count time in seconds that has passed till summing of NPP started
   seconds_year = seconds_year + delta_time
   IF (new_year) seconds_year = delta_time

  END SUBROUTINE calc_NPPbuf

#endif
END MODULE mo_assimi_process
