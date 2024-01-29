!
! Source module for computing diffusion coefficients and implicit vertical diffusion:
!
! Description of *turb_diffusion*:
!   This  module calculates the tendencies for turbulent
!   vertical transport of momentum and heat and the coefficients
!   for turbulent diffusion as well.
!
!   The clousure is made on lever 2.5 (Mellor/Yamada) using a prognostic
!   TKE-equation and includes the formulation of a flow through a porous
!   medium (roughness layer)
!
!   The turbulence model (with some Prandtl-layer approximations is used
!   for the calculation of turbulent transfer between atmosphere and the
!   lower boundary too.
!
! The module contains the public subroutines :
!
!   turbdiff
!
! called from the turbulence interface routine of the model.
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
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE turb_diffusion

!-------------------------------------------------------------------------------
!
! Documentation History:
!
!  The history of all these modifications is as follows, where those belonging to the fomal
!   reorganization of the whole package (atmospheric turbulence and surface-to-atmpsphere transfer)
!   are now in the header of MODULE 'turb_utilities', containing various common SUBs for 'turbdiff'
!   and 'turtran' (related moist thermodynamicds and the treatment of turbulent budget equations)
!   and also the blocked code for semi-implicit vertical diffusion. The new blocked version of SUB 'turbtran'
!   is now in MODULE 'turb_transfer':
!
!              2010/12/17 Matthias Raschendorfer
!  Introduction of a TKE-source term due to scale interaction with sub grid scale convection
!   using 'ltkecon' and the convective buoyant heat flux density in 'tket_conv'.
!              2011/02/18 Matthias Raschendorfer
!  Introduction of some minor formal modif. of some parts of the code mainly to achiev better vectorization
!   (results will be modyfied only because of numerical effects).
!  Introduction of OPTIONAL extra output fields 'tket_sso' and 'tket_hshr'.
!              2011/03/23 Matthias Raschendorfer
!  Substitution of run time allocations because of bad performance on some computers.
!              2011/08/26 Matthias Raschendorfer
!  Discriminating the variable for standard deviation of oversaturation (input) and cloud fraction (output)
!   in the CALL of 'turb_cloud'.
!  Changing the the definition of LOGICAL 'lcircterm'.
!  Introduction of a preconditioning of the tridiagonal system for "imode_turb=3".
!              2011/09/23 Matthias Raschendorfer
!  Changing 'k' to 'ke' in the initialization section of 'turbdiff' in case of "itype_tran.EQ.3".
!  Formulating temp-grad for 'ltmpcor' for "k.EQ.ke1" similar to the other levels.
!  Removing a wrong multiplication with Exner-factor in case of "imode_turb.EQ.4".
!  Adding a missing multiplication by 'dicke()' in case of ".NOT.limpltkediff".
!  Dirscriminating between two effective values of effictive Prandtl-layer depth,
!   by introducing the arrays 'dz0(:,mom|sca)' and 'vh0(:,:)
!   removing wrong surface conditions in case "imode_turb.GE.3 .AND. itype_tran.EQ.2".
!              2011/12/08 Matthias Raschendorfer
!   Diffusion tendency (rather than flux-density) of pot. temp. is converted in that of ordinary temp.
!              2012/01/10 Matthias Raschendorfer
!  Correction of some bugs:
!   In case of "limpltkediff" 'ko' needs to be added by "1".
!   Coding mistake in n-loop bounds for surface layer gradients in 'turbdiff'.
!  Introduction of array 'rhoh' (density on main levels) to avoid a back interpolation from half levels.
!  Using a mass weighted density interpolation onto half levels as for all other variables.
!              2012/01/11 Matthias Raschendorfer
!  Some rearrangement and simplification:
!   Removing explicit vertical diffusion and expressing explicit (eg. moist) corrections
!    by corrected vertical profiles still using implicit formulation of vertical diffusion.
!   Now "imode_turb=2" is for semi-implicit diffusion within 'turbdiff' using a concentration condition
!    at the surface and "imode_turb=3" is for the same with a flux condition at the surface.
!    "imode_turb_turb=4" is obsolet now.
!   One loop for diffusion of all first order model variables.
!              2012/01/20 Matthias Raschendorfer
!  Reformulation of roughness layer drag by an implicit equation (keeping sign of wind components)
!  Expressing explicit circulation term (as far as it has divergence form) by a respective
!   q**2-profile correction within the implicit formulation of diffusion similar to the treatment
!   of the former explicit (moist) corrections
!  Introduction of effective flux profile calculation into SUB 'calc_impl_vert_diff'.
!  Setup of a positiv definit solution avoiding some of the former limitations.
!              2012/03/20 Matthias Raschendorfer
!  Rearrangement, modularization and revision of numerical treatment, such as:
!  Intoduction of the SUBs 'solve_turb_budgets' and 'adjust_satur_equil' in order to use the same code
!  for the same purpose in 'turbtran' and 'turbdiff.
!  Removing the explicit diffusion option.
!  Introducing the new driving SUB 'vert_grad_diff' organizing vertical diffusion when called in a loop of
!  several full level variables and managing options how to treat non-gradient-fluxes and additional
!  explicit time tendencies numerically, including smoothing options and preconditioning the linear system.
!  Parameters 'imode_tran' and 'imode_turb' define whether the TKE-equation is solved in a dianostic (1)
!  or a prognostic (2) form.
!  New parameter 'lsflcnd' is a flag for using a lower flux condition for vertical diffusion.
!  Introducing the flags 'lturatm', 'ltursrf', 'lmomdif', 'lscadif' arranging what tasks of 
!  'organize_turbdiff' are required.
!  Introduction of 'ltkeshs' and removing the case "itype_sher.EQ.3".
!              2014/07/28 Matthias Raschendorfer
!  Introduction of precalculated 'hdef2', 'hdiv', 'dwdx', 'dwdy' used for calculation of horizontal shear
!   including the scale separated non-turbulent part controlled by 'itype_sher'.
!  Simpler (physically identical) formulation of moist flux conversion
!   -> only numerical differences in the case of 'lexplcor'
!  Numerically more efficient formulation of Blackadar 'len_scale'
!   -> only numerical differences
!  Eliminating array 'wind', as 'u' and 'v' are already defined at mass positions.
!              2015/08/25 Matthias Raschendorfer
! Adopting other development within the ICON-version (namely by Guenther Zaengl) as switchable options
!  related to the following new selectors and switches:
!   imode_pat_len, imode_frcsmot, imode_shshear, imode_tkvmini, imode_charpar.
! Rearranging the development by Matthias Raschendorfer that had not yet been transferred to COSMO as 
!  switchable options related to the following switches:
!   lsflcnd, ldynimp, lprecnd, ltkeshs, loutshs
!  and selectors:
!   imode_tkediff, imode_adshear, imode_tkemini
!  and a partly new (more consistent) interpretation of:
!   imode_turb, icldm_turb and  itype_sher
! Controlling numerical restrictions gradually namely by the parameters:
!  tndsmot, frcsmot
! Correction of some bugs (of the latest ICON-version) in particular related to the optional lower 
!  concentration condition.
! Using the arrays 'tvm', 'tvh' and 'tkm', allowing an easier formulation of transfer-resistances.
!              2016-05-10 Ulrich Schaettler 
! Splitting this module from the original module 'organize_turbdiff' as it was used by ICON before.
! Moving declarations, allocation and deallocations of ausxilary arrays into MODULE 'turb_data'.
!-------------------------------------------------------------------------------------------------------


! Modules used:

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif

!-------------------------------------------------------------------------------
! Parameter for precision
!-------------------------------------------------------------------------------

USE mo_kind,         ONLY : wp, vp

!-------------------------------------------------------------------------------
! Mathematical and physical constants
!-------------------------------------------------------------------------------

USE mo_mpi,                ONLY : get_my_global_mpi_id
USE mo_exception,          ONLY : finish
USE mo_physical_constants, ONLY : &
!
! Physical constants and related variables:
! -------------------------------------------
!
    r_d      => rd,       & ! gas constant for dry air
    rvd_m_o  => vtmpc1,   & ! r_v/r_d - 1
    cp_d     => cpd,      & ! specific heat for dry air
    lh_v     => alv,      & ! evaporation heat
    lhocp    => alvdcp,   & ! lh_v / cp_d
    con_m,                & ! kinematic vsicosity of dry air (m2/s)
    con_h,                & ! scalar conductivity of dry air (m2/s)
    grav                    ! acceleration due to gravity

!-------------------------------------------------------------------------------
! From Flake model
!-------------------------------------------------------------------------------

USE sfc_flake_data, ONLY: &
    h_Ice_min_flk      ! Minimum ice thickness [m]

!-------------------------------------------------------------------------------
! Turbulence data (should be the same in ICON and COSMO)
!-------------------------------------------------------------------------------

USE turb_data, ONLY : &

! Numerical constants and parameters:
! -----------------------------------

    tkhmin,       & ! minimal diffusion coefficients for heat
    tkmmin,       & ! minimal diffusion coefficients for momentum
    tkhmin_strat, & ! additional minimal diffusion coefficients for heat for stratosphere
    tkmmin_strat, & ! additional minimal diffusion coefficients for momentum for stratosphere
    tndsmot,      & ! vertical smoothing factor for diffusion tendencies
    frcsmot,      & ! vertical smoothing factor for TKE forcing
    imode_frcsmot,& ! mode for TKE forcing
    epsi,         & ! relative limit of accuracy for comparison of numbers
    it_end,       & ! number of initialization iterations (>=0)

! Parameters describing physical properties of the lower boundary 
! of the atmosphere:
!---------------------------------------------------------------

    pat_len,      & ! lenth scale of subscale patterns over land [m]
    len_min,      & ! minimal turbulent length scale [m]
    vel_min,      & ! minimal velocity scale [m/s]
    akt,          & ! von Karman-constant
!
    d_h=>d_heat,  & ! factor for turbulent heat dissipation
    d_m=>d_mom,   & ! factor for turbulent momentum dissipation
!
    c_diff,       & ! factor for turbulent diffusion of TKE
    a_hshr,       & ! factor for horizontal shear production of TKE
    c_scld,       & ! factor for liquid water flux density in sub grid scale clouds

    ! derived parameters calculated in 'turb_setup'
    tet_g, rim, b_m, b_h, sm_0, sh_0,   &
    d_1, d_2, d_3, d_4, d_5, d_6,       &
    a_3, a_5 ,a_6,                      &
    tur_rcpv, tur_rcpl,                 &

    ! used derived types
    modvar, turvar, varprf, & !

! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

    ltkesso,      & ! consider SSO-wake turbulence production of TKE
    ltkecon,      & ! consider convective buoyancy production of TKE
    ltkeshs,      & ! consider separ. horiz. shear production of TKE 
    loutshs,      & ! consider separ. horiz. shear production of TKE for output
    lnonloc,      & ! nonlocal calculation of vertical gradients used for turb. diff.
    lprfcor,      & ! using the profile values of the lowest main level instead of
    ltmpcor,      & ! consideration of thermal TKE-sources in the enthalpy budget
    lexpcor,      & ! explicit corrections of the implicit calculated turbul. diff.

!   for semi-implicit vertical diffusion:
    lsflcnd,      & ! lower flux condition
    ldynimp,      & ! dynamical calculation of implicit weights
    lprecnd         ! preconditioning of tridiagonal matrix

USE turb_data, ONLY : &

! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------
!
    imode_turb,   & ! mode of TKE-equation in turbulence scheme                  (compare 'imode_tran')
                    !  0: diagnostic equation
                    !  1: prognostic equation (default)
                    !  2: prognostic equation (implicitly positive definit)
    icldm_turb,   & ! mode of water cloud representation in turbulence parametr. (compare 'icldm_tran')
                    ! -1: ignoring cloud water completely (pure dry scheme)
                    !  0: no clouds considered (all cloud water is evaporated)
                    !  1: only grid scale condensation possible
                    !  2: also sub grid (turbulent) condensation considered
    itype_sher,   & ! type of shear production for TKE
                    ! 0: only vertical shear of horizontal wind
                    ! 1: previous plus horizontal shear correction
                    ! 2: previous plus shear from vertical velocity
    ilow_def_cond,& !type of the default condition at the lower boundary
                    ! 1: zero surface flux density
                    ! 2: zero surface value
    imode_calcirc,& ! mode of treating the circulation term (related to 'pat_len', imode_pat_len')
                    ! 1: explicit calculation of the flux convergence
                    ! 2: quasi implicit treatment by calculation of effective TKE-gradients
    imode_shshear,& ! mode of calculat. the separated horizontal shear mode related to 'ltkeshs', 'a_hshr')
                    ! 1: with a constant lenght scale 
                    ! 2: with a Ri-dependent length sclale correction
    imode_tkesso,&  ! mode of calculat. the SSO source term for TKE production
                    ! 1: original implementation
                    ! 2: with a Ri-dependent reduction factor for Ri>1
                    ! 3: as 2, but additional reduction for mesh sizes < 2 km
    imode_tkvmini,& ! mode of calculating the minimal turbulent diff. coeffecients
                    ! 1: with a constant value
                    ! 2: with a stability dependent correction
    imode_syndiag,& ! mode of diagnostics at the synoptic near surface levels (related to 'itype_diag_t2m')
                    ! 1: direct interpolation of temperature and specific humidity
                    ! 2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                    !    allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
    imode_trancnf,& ! mode of configuring the transfer-scheme 
                    ! 1: old version: start. with lamin. diffus.; with a lamin. correct. for profile-funct.;
                    !    interpol. T_s rather then Tet_l onto zero-level; calcul. only approx. Tet_l-grads.;
                    !    using an upper bound for TKE-forcing; without transmit. skin-layer depth to turbul.
                    ! 2: 1-st ConSAT: start. with estim. Ustar, without a laminar correct. for prof.-funct.;
                    !    interpol. Tet_l onto zero-level; calcul. Tet_l-gradients directly; 
                    !    without an upper bound for TKE-forcing; with transmit. skin-layer depth to turbul.
                    ! 3: 2-nd ConSAT: as "2", but with a hyperbolic interpol. of profile function
                    !    for stable stratification
                    ! 4: 3-rd ConSAT: as "3", but without using an upper interpolation node
    imode_tkediff,& ! mode of implicit TKE-Diffusion (related to 'c_diff')
                    ! 1: in terms of q=SQRT(2*TKE))
                    ! 2; in terms of TKE=0.5*TKE**2
    imode_adshear,& ! mode of considering additional shear by scale interaction
                    ! 1: not considered for stability functions
                    ! 2:     considered for stability functions
    imode_tkemini   ! mode of fixing a lower limit of q=2TKE**2
                    ! 1: by using 'vel_min' only
                    ! 2: by adapting to minimal diffusion coefficients

USE turb_data, ONLY : &

    ! numbers and indices

    nvel    ,     & ! number of velocity components
    naux    ,     & ! number of auxilary variables
    ndim    ,     & !
    nmvar   ,     & ! number of included prognostic model variables
                  !Verrueckungen invarianten Groessen ('tet_l', 'h2o_g')

    ntyp    ,     & ! number of variable types (mom) und (sca)

    mom     ,     & ! index for a momentum variable
    sca     ,     & ! index for a scalar   variable
    u_m     ,     & ! index for mass centered zonal      velocity compont
    v_m     ,     & ! index for mass centered meridional  ,,         ,,
    tet_l   ,     & ! index for liquid water potential temperature
    tet     ,     & ! index for potential temperature
    tem     ,     & ! index for temperature
    h2o_g   ,     & ! index for toatal water
    vap     ,     & ! index for water vapor
    liq             ! index for liquid water


!-------------------------------------------------------------------------------
! Control parameters for the run
!-------------------------------------------------------------------------------

! ICON data have to be declared for these variables, which is done later on
!-------------------------------------------------------------------------------

USE turb_utilities,          ONLY:   &
    turb_setup,                      &
    adjust_satur_equil,              &
    solve_turb_budgets,              &
    vert_grad_diff,                  &
    prep_impl_vert_diff,             &
    calc_impl_vert_diff,             &
    vert_smooth,                     &
    bound_level_interp,              &
    zexner
    
!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
    lsclm, latmflu, i_cal, i_mod, imb, &
    SHF, LHF
#endif
!SCLM---------------------------------------------------------------------------

!===============================================================================

IMPLICIT NONE

PUBLIC  :: turbdiff

!===============================================================================

REAL (KIND=wp), PARAMETER :: &
    z0 = 0.0_wp,    &
    z1 = 1.0_wp,    &
    z2 = 2.0_wp,    &
    z3 = 3.0_wp,    &
    z4 = 4.0_wp,    &
    z5 = 5.0_wp,    &
    z6 = 6.0_wp,    &
    z7 = 7.0_wp,    &
    z8 = 8.0_wp,    &
    z9 = 9.0_wp,    &
    z10=10.0_wp,    &

    z1d2=z1/z2     ,&
    z1d3=z1/z3     ,&
    z2d3=z2/z3     ,&
    z3d2=z3/z2

!===============================================================================

CONTAINS

!===============================================================================

#  define err_args

SUBROUTINE turbdiff ( &
!
          iini, ltkeinp, lstfnct, l3dturb,                           &
!
          lrunsso, lruncnv, lrunscm, lsfluse,                        &
!
          dt_var,dt_tke, nprv, ntur, ntim,                           &
!
          nvec, ke, ke1, kcm, iblock, ivstart, ivend,                &
!
          l_hori, hhl, dp0, trop_mask, innertrop_mask,               &
!
          gz0, l_pat, c_big, c_sml, r_air,                           &
!
          t_g, qv_s, ps,                                             &
          u, v, w, t, qv, qc, prs, rhoh, rhon, epr,                  &
!
          impl_weight,                                               &
!
          tvm, tvh, tfm, tfh, tkred_sfc, tkred_sfc_h,                &
          tke, tkvm, tkvh, rcld, tkhm, tkhh,                         &
          hdef2, hdiv, dwdx, dwdy,                                   &
!
          edr, tket_sso, tket_conv, tket_hshr,                       &
          u_tens, v_tens, t_tens,                                    &
          qv_tens, qc_tens,                                          &
          tketens, tketadv,                                          &
          ut_sso, vt_sso,                                            &
!
          shfl_s, qvfl_s,                                            &
!
          zvari,                                                     &
!
          tet_flux, vap_flux, liq_flux,l_3d_turb_fluxes,             &
          tket_therm, tket_mech, tket_shear                          &
!
          err_args)

!-------------------------------------------------------------------------------
!
! 
! All tendency parameters are OPTIONAL (except 'tketens' in case of "lturatm=T". If they are missing,
!  calculated tendencies of SUB 'turbdiff' are automatically added to the related prognostic variables.
!
! Note:
! It is also possible to use only one time level for TKE using "ntim=1" and thus "nprv=1=ntur".
!
! Description:
!
!     Es werden die Diffusionskoeffizienten berechnet und ggf. Anteile
!     der zeitlichen Tendenzen der turbulenten Diffusion bestimmt
!     und zu den Tendenzfeldern hinzuaddiert.
!     Optional wird eine explizite oder (teil-)implizite Berechnung der
!     Diffusionstendenzen oder aber nur eine Berechnung der Diffusions-
!     koeffizienten durchgefuehrt. Im letzten Fall wird dann ein
!     implizit zu berechnender Anteil der Diffusionstendenzen an
!     anderer Stelle (slow_tendencies) bestimmt.
!     Allerdings koennen dann zusaetzliche explizite Korrekturtendenzen
!     hier in tubdiff bestimmt werden.
!
! Method:
!
!     Die Berechnung basiert auf einer Schliessung 2-ter Ordnung auf
!     dem level 2.5 (nach Mellor/Yamada). Demnach wird also eine
!     prognostische Gleichung fuer die TKE geloest.
!     Ausser der TKE-Advektion, die zusammen mit den Advektionstendenzen
!     der anderen prognostischen Variablen an anderer Stelle berechnet
!     wird, geschieht die gesamte TKE-Prognose in diesem Unterprogramm.

!     Die Formulierung des Schemas erfolgt mit thermodynamischen
!     Variablen, die bei feuchtadiabatischen vertikalen Verrueckungen
!     erhalten bleiben (pot. Fluessigw.temp. und Gesamtwassergehalt),
!     so dass der Kondesationseffekt auf subskalige Vertikalbewegungen
!     beruecksichtigt wird.
!     Die turbulenten Flussdichten der Erhaltungsgroessen werden in
!     solche der Modellvariablen konvertiert, so dass die thermodyn.
!     Kopplung der Flussdichten richtig erhalten bleibt.

!     Angeschlossen ist auch ein optionales statistisches Wolkenschema
!     (nach Sommeria und Deardorff), sub turb_cloud, welches auch
!     subskalige Bewoelkung mit Hilfe der ueber das Feld rcld ausge-
!     gebenen Standardabweichung des Saettigungsdefizites berechnet.

!     Das Turbulenzschema wurde so verallgemeinert, dass es auch bei
!     einer vertikal aufgeloesten Bestandesschicht gueltig ist, indem
!     idealisierend von der Durchstroemung eines poroesen Mediums
!     ausgegangen wird. Die Bilanzgleichungen 1-ter und 2-ter Ordnung
!     enthalten dann zusaetzliche Terme, welche die Wechselwirkungen
!     mit der Bestandes-Matrix beschreiben. Dies wirkt sich zum einen
!     auf die Stabilitaetsfunktionen und zum anderen vor allem auf die
!     TKE-Gleichung aus, welche einen auf den Formwiderstand der
!     Bestandeselemente zurueckzufuehrenden zusaetzlichen Quellterm
!     (Nachlaufturbulenz) enthaelt. Ausserdem werden die turbulenten
!     Flussdichtedivergenzen noch um einen Zusatzterm, welcher der
!     Reduktion des lufterfuellten Volumens im Gitterelement Rechnung
!     traegt, erweitert. Der Effekt des Formwiderstandes in der
!     Impulsgleichung ist ebenfalls beruecksichtigt. Die zusaetzlichen
!     Tendenzterme, die auf die Flussdichten zwischen Bestandes-Matrix
!     und umgebender Luft zurueckzufuehren sind (Bestandesquellen),
!     muessen noch in einem separaten Bestandesmodell parametrisiert
!     werden und sind nicht Gegenstand des Turbulenzmodells.
!     Schliesslich wird auch der Effekt der Transformation von Turbulenz
!     auf der dominierenden Skala in kleinskalige dissipative Turbulenz
!     durch Wirbelbrechen an Koerpern mit
!          Laengenskalen der Abmessungen << turbulente Laengenskala
!     beruecksichtigt; was sich durch eine (von der Laengenskala und
!     Volumendichte jener sehr kleinen Bestandeselemente aubhaengige)
!     Modifikation der Modellkonstanten ausdruecken laesst.

!     Es wird auch versucht den Effekt thermisch induzierter
!     Zirkulationen auf die TKE-Produktion zu beruecksichtigen
!     Hierdurch wird (vor allem) der Austauch in der naechtlichen
!     Grenzschicht erhoeht, was der Tendenz des alten Schemas,
!     in Bodennaehe zu kalte und nicht schnell genug anwachsende
!     Inversionen zu produzieren, entgegenwirkt.

!     Optional kann die Berechnung der vertikalen Gradienten durch eine
!     nicht-lokale Variante erfolgen. Hierbei werden die Gradienten
!     mit Profilen gebildet, die mit einem ueber die stabilitaets-
!     abhaengige Laengenskala gebildeten gleitenden Mittel behandelt
!     wurden.

!     Die Bildung der Anteile der Diffusionstendenzen, die numerisch
!     durch die Multiplikation mit einer Tridiagonalmatrix ausdrueckbar
!     sind, kann (neben der expliziten Variante) auch implizit erfolgen
!     (im Falle der Berechnung von nich-lokalen Gradienten und fuer die
!     TKE allerdings nur explizit).
!     Bei expliziter Rechnung ist, um auch bei Zeitschritten von
!     mehreren Minuten numerisch stabil zu bleiben, eine Limitierung
!     der Groesse der Diffusionskoeffezienten und eine im vertikalen Integral
!     quellenfreie numerische Glaettung der Vertikalprofile der Diffusions-
!     tendenzen erforderlich, sowie eine teilimplizite Behandlung der
!     Diffusionstendenzen in der untersten Modellschicht, erforderlich.

!     Die unteren Randwerte der turbulenten Flussdichten werden ueber
!     die Transferkoeffizienten zwischen Erdboden und unterster
!     Modellschicht (tcm und tch) bestimmt.
!     Optional koennen die Transferkoeffizienten auch mit diesem
!     Unterprogramm bestimmt werden, indem das Turbulenzmodell auch auf
!     das Niveau d+z0 angewandt wird, wobei vertikale Gradienten in
!     diesem Niveau mit Hilfe der Prandtl-Schicht-Hypothese berechnet
!     werden.
!     In diesem Zusammenhang wird auch die Wirkung der laminaren
!     Grenzschicht behandelt.
!
!     Turbulente Horizontaldiffusion (um ein 3-d-Schema zu erhalten)
!     ist noch nicht enthalten, kann aber integriert werden.
!     Uebergabevariablen:
!
!-------------------------------------------------------------------------------

! Declarations
!-------------------------------------------------------------------------------

!Formal Parameters:
!-------------------------------------------------------------------------------

! Parameters controlling the call of 'turbtran':
! ----------------------------------------------

LOGICAL, INTENT(IN) :: &
  lstfnct,      & !calculation of stability function required

  l3dturb,      & !a model run with 3D-(turbulent)-diffusion

  ltkeinp,      & !TKE present as input (at level k=ke1 for current time level 'ntur')

  lruncnv,      & !convection scheme is active
  lrunsso,      & !SSO-Scheme is active
  lrunscm,      & !a Single Column run (default: FALSE)
  lsfluse         !use explicit heat flux densities at the surface

REAL (KIND=wp), INTENT(IN) :: &
  dt_var,       & !time step for ordinary prognostic variables
  dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER,        INTENT(IN) :: &
  iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
  ntur,         & !current  time level of 'tke' valid after  prognostic incrementation
  nprv,         & !previous time level of 'tke valid before prognostic incrementation
  ntim            !number of 'tke' time levels


! Horizontal and vertical sizes of the fields and related variables:
! ------------------------------------------------------------------

INTEGER,        INTENT(IN) :: &
  nvec,         & ! number of horizontal grid points in the nproma-vector
  ke,           & ! index of the lowest main model level
  ke1,          & ! index of the lowest model half level (=ke+1)
  kcm,          & ! level index of the upper canopy bound
  iblock


! Start- and end-indices for the computations in the horizontal layers:
! ---------------------------------------------------------------------

INTEGER,        INTENT(IN) :: &
  ivstart,      & ! start index in the nproma vector
  ivend           ! end index in the nproma vector


! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! -----------------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
  hhl             ! height of model half levels                   ( m )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
  dp0             ! pressure thickness of layer                   (pa )

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
  l_pat  ,      & ! effektive Laengenskala der therm. Inhomogenitaeten 
                  ! der Erdbodenoberflaeche
  l_hori          ! horizontal grid spacing (m)
 
REAL (KIND=wp), DIMENSION(:,kcm:), TARGET, OPTIONAL, INTENT(IN) :: &
  c_big,        & ! effective drag coefficient of canopy elements
                  ! larger than or equal to the turbulent length scale (1/m)
  c_sml           ! effective drag coefficient of canopy elements
                  ! smaller than the turbulent length scale            (1/m)

REAL (KIND=wp), DIMENSION(:,kcm-1:), TARGET, OPTIONAL, INTENT(IN) :: &
  r_air           ! log of air containing fraction of a gridbox inside
                  ! the canopy                                          (1)


! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
  ps,           & ! surface pressure                              ( pa  )
  qv_s,         & ! specific water vapor content on the surface   (kg/kg)
  t_g             ! weighted surface temperature                  (  k  )



! Atmospheric model variables:
! ---------------------------------

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
  u,            & ! zonal wind speed       (at mass positions)    ( m/s )
  v,            & ! meridional wind speed  (at mass positions)    ( m/s )
  t,            & ! temperature                                   (  k  )
  qv,           & ! specific water vapor content                  (kg/kg)
  qc              ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
  prs             ! atmospheric pressure                          ( pa  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
  rhoh,         & ! total density of air                          (kg/m3)
  epr             ! exner pressure                                 (1)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(OUT) :: &
  rhon            ! total density of air                          (kg/m3)

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
  w               ! vertical wind speed (defined on half levels)  ( m/s )

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
  impl_weight     ! profile of precalculated implicit weights 



! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
  gz0,           & ! roughness length * g of the vertically not
                   ! resolved canopy                               (m2/s2)
  !Achtung: Der g-Faktor ist ueberfluessig!

  ! turbulent (transfer) velocity scales at the surface
  tvm,           & ! for momentum                                  ( m/s)
  tvh              ! for heat and moisture                         ( m/s)

  !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
  !vecolities  'tvm' and 'tvh' makes live much easier!!               

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
  ! turbulent transfer factors for laminar- and roughness-layer transfer
  tfm,           & ! of momentum                                     --
  tfh              ! of scalars                                      --

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(IN) :: &
  tkred_sfc, tkred_sfc_h   ! reduction factors for minimum diffusion coefficients near the surface


! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=wp), DIMENSION(nvec,ke1,ntim), TARGET, INTENT(INOUT) :: &
  tke              ! q:=SQRT(2*TKE); TKE='turbul. kin. energy'     ( m/s )
                   ! (defined on half levels)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
  tkvm,          & ! turbulent diffusion coefficient for momentum  (m2/s )
  tkvh             ! turbulent diffusion coefficient for heat      (m2/s )
                   ! (and other scalars)

REAL (KIND=wp), DIMENSION(nvec,ke1,ndim), TARGET, INTENT(OUT) :: &
  zvari            ! to give values to vertical diffusion

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
  rcld             ! standard deviation of the local oversaturation
                   ! (as input and output)
                   ! fractional cloud cover (in turbdiff)            --

REAL (KIND=vp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
  hdef2,         & ! horizontal deformation square at half levels  ( 1/s2 )
  hdiv,          & ! horizontal divergence                   ,,    ( 1/s )

  dwdx,          & ! zonal      derivative of vertical wind  ,,    ( 1/s )
  dwdy             ! meridional derivative of vertical wind  ,,    ( 1/s )


! Tendency fields for the prognostic variables:
! -----------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
  u_tens,        & ! u-tendency                                    ( m/s2)
  v_tens,        & ! v-tendency                                    ( m/s2)
  t_tens,        & ! t-tendency                                    ( K/s )
  qv_tens,       & ! qv-tendency                                   ( 1/s )
  qc_tens          ! qc-tendency                                   ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
  tketens,       & ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)
  tketadv          ! advection tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:),   TARGET, OPTIONAL, INTENT(IN)    :: &
  trop_mask,     & ! mask-factor (1: within tropics; 0: within extra-tropics)
                   ! used for vertical smoothing of TKE forcing terms
  innertrop_mask

REAL (KIND=wp), DIMENSION(:,:),         OPTIONAL, INTENT(IN)    :: &
  ut_sso,        & ! u-tendency due to the SSO-Scheme              ( 1/s )
  vt_sso           ! v-tendency due to the SSO-Scheme              ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(OUT)    :: &
  edr,           & ! eddy dissipation rate of TKE (EDR)            (m2/s3)
  tket_sso,      & ! TKE-tendency due to SSO wake production       (m2/s3)
  tket_hshr,     & ! TKE-tendency due to separ. horiz. shear       (m2/s3)
  tkhm,          & ! horizontal diffusion coefficient for momentum ( m2/s )
  tkhh             ! horizontal diffusion coefficient for scalars  ( m2/s )

REAL (KIND=wp), DIMENSION(:,:),         OPTIONAL, INTENT(IN)    :: &
  tket_conv        ! TKE-tendency due to convective buoyancy       (m2/s3)

REAL (KIND=wp), DIMENSION(:),   TARGET, OPTIONAL, INTENT(INOUT) :: &
  shfl_s,        & ! sensible heat flux at the surface             (W/m2)    (positive downward)
  qvfl_s           ! water vapor   flux at the surface             (kg/m2/s) (positive downward)

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(OUT)    :: &
  tet_flux,      & ! vertical flux of theta (K/m2s)  (positive upward)   
  vap_flux,      & ! vertical flux of qv    (kg/m2s) (positive upward)   
  liq_flux,      & ! vertical flux of qc    (kg/m2s) (positive upward)   
  tket_therm,    & ! tke tendency due to buoyancy (m2/s3)
  tket_mech,     & ! tke tendency due to mechanical forcing (m2/s3)
  tket_shear       ! tke tendency due to vertical shear (m2/s3)

LOGICAL, OPTIONAL, INTENT(IN) :: l_3d_turb_fluxes


!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

INTEGER ::        &
  i, k,           & !horizontaler und vertikaler Laufindex
  kem,            & !ke oder ke1
  nvor,           & !laufende Zeittstufe des bisherigen TKE-Feldes
  it_start,       & !Startindex der Iterationen
  it_durch,       & !Durchgangsindex der Iterationen
  ndiff             !number of 1-st order variables

LOGICAL ::        &
  lini,           & !initialization required
  lssintact         !trenne Skalen-Interaktionsterme vom mech. Forcing ab

REAL (KIND=wp) :: &
  fr_tke              ! z1/dt_tke

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
! pointer for density and eddy dissipation rate:
  prhon(:,:), prhoh(:,:), ediss(:,:)


! Lokale logical Variablen:

LOGICAL ::          &
  ldotkedif,        & !berechne (teil-)implizite Vert.diff von TKE
  lcircterm,        & !Zirkulationsterm muss berechnet werden
  lcircdiff           !Zirkulationsterm wird zusammen mit TKE-Diffusion bestimmt

! Lokale Integer-Hilfsvariablen:

INTEGER ::          &
  ii, n, kk,        & !Indices fuer diverse Schleifen
  ku,k1,k2,         & !Schicht-Indices
  itndcon             !Index fuer Modus der  Tendenzberuecksichtigung

! Lokale real Variablen:

REAL (KIND=wp) ::   &

! Hilfsvariablen:
  wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
  fakt, fakt2,      & !  ,,         ,,     ,,     Faktoren

! Platzh. fuer thermodynamische Hilfsgreossen
  flw_h2o_g,        & !                 rc/(1+(lh_v/cp_d)*d_qsat/d_T)
  flw_tet_l           !epr*d_qsat/d_T*rc/(1+(lh_v/cp_d)*d_qsat/d_T)

REAL (KIND=wp) ::   &

! Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2
  fh2,fm2,          &

! Platzh. fuer horiz. Geschw.-Komponenten und bel. Geschw.:
  vel1,vel2,velo,   &

! Platzh. fuer die Hoehe ueber Grund, Hoehendifferenzen, obere und
! untere Hoehe, turbulente Laengaenskalen, Kohaerenzlaenge,
! Dicke der laminaren Grenzschicht,sowie eine Laenge allgemein:
  h,hu,l_turb,lh,   &
  lm,kohae_len,     &
  com_len,          &
  edh,              & ! Kehrwert von Schichtdicken

! Zwischenspeicher fuer
  thermik,          & !(negative) Auftriebsproduktion an TKE
  phasdif,          & !Temperaturtendenz durch Phasendiffusion

! Tuning
  x4, x4i


! Local arrays:

INTEGER ::          &
  ivtp(nmvar)         ! index of variable type

! Time increment and inverse time increment of ordinary prognostic variables:
REAL (KIND=wp), TARGET :: &
     tinc(nmvar)      !

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
! Pointer fuer Tendenzfelder:
  utens(:,:), vtens(:,:), ttens(:,:), qvtens(:,:), qctens(:,:)

! Note:
! The following buffers wouldn't be necessary, if the related pointers above
! were allowed to be allocated at run time:

REAL (KIND=wp), DIMENSION(:,:), POINTER, CONTIGUOUS :: &
  cur_prof, upd_prof, sav_prof, &
  expl_mom, impl_mom, invs_mom, &
  eff_flux

LOGICAL ::          &
  ltend(nmvar),     &  !calculation of tendencies required
  lsfli(nmvar),     &  !surface value input is a flux density instead of a concentration
  lcalc_frcsmot        !local control variable if smoothing of TKE forcing terms needs to be calculated

TYPE (varprf) :: pvar(naux+1) !vertical variable profiles

! these fields are still taken as local arrays, because the CRAY compiler cannot do the
! same optimizations with OMP threadprivate variables

REAL (KIND=wp), TARGET ::  &
  ! targets of used pointers
  diss_tar   (nvec,ke1)      ! target for eddy dissipation rate (m2/s3)

REAL (KIND=wp), TARGET ::  &
  ! internal atmospheric variables
  len_scale(nvec,ke1),     & ! turbulent length-scale (m)
  hor_scale(nvec,ke),      & ! effective hoprizontal length-scale used for sep. horiz. shear calc. (m)

  l_scal   (nvec),         & ! reduced maximal turbulent length scale due to horizontal grid spacing (m)

  fc_min   (nvec),         & ! minimal value for TKE-forcing (1/s2)

  shv      (nvec,ke1),     & ! velocity scale of the separated horiz. shear mode (m/s)
  frh      (nvec,ke1),     & ! thermal forcing (1/s2) or thermal acceleration (m/s2)
  frm      (nvec,ke1),     & ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)
  ftm      (nvec,ke1),     & ! mechan. forcing (1/s2) by pure turbulent shear 

  prss     (nvec,ke1:ke1), & ! surface pressure (Pa)
  tmps     (nvec,ke1:ke1), & ! surface temperature (K)
  vaps     (nvec,ke1:ke1), & ! surface specific humidity
  liqs     (nvec,ke1:ke1), & ! liquid water content at the surface

  dicke    (nvec,ke1),     & ! any (effective) depth of model layers (m) or other auxilary variables
  hlp      (nvec,ke1),     & ! any 'help' variable

  zaux     (nvec,ke1,ndim),& ! auxilary array containing thermodynamical properties
                             ! (dQs/dT,ex_fakt,cp_fakt,g_tet,g_vap) or various
                             ! auxilary variables for calculation of implicit vertical diffusion

  can      (nvec,kcm:ke1), & ! auxilary array valid for the vertically resolved canopy
  lay      (nvec),         & ! any variable at a specific layer
  lays     (nvec,2),       & ! any (2-D) vector of variables at a specific layer

  src      (nvec),         & ! effective depth of Prandtl-layer applied to scalars  (m)

  dzsm     (nvec),         & ! effective depth of Prandtl-layer applied to momentum (m)
  dzsh     (nvec)            ! effective depth of Prandtl-layer applied to scalars  (m)

REAL (KIND=wp)         ::  &
  grad     (nvec,nmvar),   &  ! any vertical gradient
  xri      (nvec,ke),      &  ! tunning
  hig      (nvec,2)           ! obere und untere Referenzhoehe bei der Bildung
                                ! nicht-lokaler Gradienten


! Eingrenzende Hoehenvieaus bei der Berechnung der
! gemittelten Profile bei nicht-lokaler Gradient-Berechnung:
INTEGER                ::  &
  lev(nvec,2)

LOGICAL, PARAMETER :: ldebug=.FALSE.

INTEGER :: my_cart_id, my_thrd_id

LOGICAL :: lzacc

!---- End of header ------------------------------------------------------------

!===============================================================================

!All variables and their tendencies are defined at horizontal mass positions.

 lzacc = iini == 0

 lssintact=((ltkesso.OR.ltkeshs.OR.ltkecon) .AND. imode_adshear.EQ.1)

 ndiff=nmvar      !number of 1-st order variables used in the turbulence model
                  !without additional tracer: these are treated in vertdiff

    ! from turb_data:
    !   nmvar = nscal+nvel
    !   nvel  = 2    active horizontal wind components:  'u_m', 'v_m'
    !                u_m = 1
    !                v_m = 2
    !   nscal = 3    active scalar variables 1st order:  'tem', 'vap', 'liq'
    !                tem   = 3:   temperature
    !                vap   = 4:   water vapor mixing ration
    !                liq   = 5:   liquid water
    !
    !      but also: tem_l = 3:   liquid water temperature
    !                tet   = 3:   potential temperature
    !                tet_l = 3:   moist (liquid water?) potential temperature
    !                h2o_g = 4:   total water content

 kem=ke

 IF (PRESENT(edr)) THEN
    ediss => edr
 ELSE
    ediss => diss_tar
 END IF

!-------------------------------------------------------------------------------

!===============================================================================

!     Fuer die Turb.par. benutzter Variablensatz auf Hauptflaechen:
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)

!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     und hat einen der Werte u_m,v_m,tet_l,h2o_g,liq;
!                        bzw. tem,vap,liq.
!     Der letzte Index bezeichnet die physik. Bedeutung der Variablen
!     Bei k=ke1 stehen die unteren Randwerte der Prandtlschicht
!     (skin-layer-Werte)
!     Am Ende des Programms wird zvari mit den (Co-Varianzen) der
!     Geschwindigkeitskomponenten ueberschrieben, die fuer die
!     Berechung der Horizontaldiffusion benoetigt werden.
!     zvari() enthaelt spaeter auch die nmvar (nichtlokalen) vertikalen
!     Gradienten und auch die durch Wirkung der subskaligen Kondensation
!     veraenderten (effektiven) vertikalen Gradienten.
!     Zum Schluss enthaelt zvari() fuer die turbulente Horizontaldiff.
!     benoetigte Komponenten des turbulenten Spannungstensors.

!########################################################################

  ltend(u_m)=PRESENT(u_tens)
  IF (ltend(u_m)) THEN    ! calculation of tendencies required
    utens => u_tens       ! 'utens' points to the tendency
  ELSE                    ! update of ordinary prognostic variables required
    utens => u            ! 'utens' points to the prognostic variables
  END IF
  ltend(v_m)=PRESENT(v_tens)
  IF (ltend(v_m)) THEN
    vtens => v_tens
  ELSE
    vtens => v
  END IF
  ltend(tem)=PRESENT(t_tens)
  IF (ltend(tem)) THEN
    ttens => t_tens
  ELSE
    ttens => t
  END IF
  ltend(vap)=PRESENT(qv_tens)
  IF (ltend(vap)) THEN
    qvtens => qv_tens
  ELSE
    qvtens => qv
  END IF
  ltend(liq)=PRESENT(qc_tens)
  IF (ltend(liq)) THEN
    qctens => qc_tens
  ELSE
    qctens => qc
  END IF

  ! check if vertical smoothing of TKE forcing terms is needed
  IF (frcsmot > z0) THEN
    IF (.NOT. PRESENT(trop_mask) .OR. imode_frcsmot == 1) THEN
      lcalc_frcsmot = .TRUE.
    ELSE IF (ANY(trop_mask(ivstart:ivend) > z0)) THEN
      lcalc_frcsmot = .TRUE.
    ELSE
      lcalc_frcsmot = .FALSE.
    ENDIF
  ELSE
    lcalc_frcsmot = .FALSE.
  ENDIF

  lsfli(:)=.FALSE. !surface values are concentrations by default

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
  IF (lsclm) THEN
    IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
      !measured SHF has to be used for forcing:
      lsfli(tem)=.TRUE.
    END IF
    IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN
      !measured LHF has to be used for forcing:
      lsfli(vap)=.TRUE.
    END IF
  END IF
  !Note: the measured SHF and LHF have to be present by shfl_s and qvfl_s!
#endif
!SCLM --------------------------------------------------------------------------------

  IF (lsfluse) THEN !use explicit heat flux densities at the surface
    lsfli(tem)=.TRUE.; lsfli(vap)=.TRUE.
  END IF

  !Begin of GPU data region
  !Input
  !$ACC DATA &
  !Working arrays                                                        !
  !$ACC   CREATE(ivtp, tinc, hig, ltend, lsfli) &
  !$ACC   CREATE(diss_tar) PRESENT(c_big, c_sml, r_air) &
  !$ACC   CREATE(len_scale, hor_scale, xri, l_scal, fc_min, ediss) &
  !$ACC   CREATE(shv, frh, frm, ftm, prss, tmps, vaps, liqs, dicke) &
  !$ACC   CREATE(hlp, zaux, can, lay, lays, src, dzsm, dzsh, grad, hig, lev) &
  !$ACC   IF(lzacc)

  !Note ACC : optional hdef2,hdiv,dwdx,dwdy,tketens,tketadv,trop_mask,ut_sso,vt_sso,edr,
  ! tket_sso,tket_hshr,tkhm,tkhh,tket_conv,shfl_s,qvfl_s have separate data region 

  !Note ACC: pointer not assigned cur_prof, upd_prof, sav_prof have separate data region       


  !Note:
  !If a tendency field of an ordinary prognostic variable is not present,
  !the related time step increment due to turbulent diffusion will be
  !added to the prognostic variable directly.

  fakt=z1/dt_var

  DO n=1,ndiff
    IF (ltend(n)) THEN  ! calculation of tendencies required
      tinc(n)=z1        ! no time increment multiplication for tendencies
    ELSE                ! update of prognostic variables required
      tinc(n)=dt_var    ! time increment multiplication for tendencies
    END IF
    IF (n.LE.nvel) THEN
      ivtp(n)=mom
    ELSE
      ivtp(n)=sca
    END IF
  END DO
  !$ACC UPDATE DEVICE(tinc, ivtp) ASYNC(1) IF(lzacc)

  IF (l3dturb .AND..NOT. (PRESENT(tkhm) .AND. PRESENT(tkhh))) THEN
    CALL finish("", 'ERROR *** 3D-diffusion with not present horiz. diff.coeffs. ***')
  END IF

!-------------------------------------------------------------------------------
  CALL turb_setup (ivstart=ivstart, ivend=ivend, ke1=ke1, &
                   iini=iini, dt_tke=dt_tke, nprv=nprv, l_hori=l_hori, qc_a=qc(:,ke), &
                   lini=lini, it_start=it_start, nvor=nvor, fr_tke=fr_tke, &
                   l_scal=l_scal, fc_min=fc_min, liqs=liqs(:,ke1), rcld=rcld, tfm=tfm, tfh=tfh, &
                   lacc=lzacc)
  !$ACC WAIT(1)
!-------------------------------------------------------------------------------

my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
my_thrd_id = omp_get_thread_num()
#endif

! Just do some check printouts:
  IF (ldebug) THEN
    DO i = ivstart, ivend
      IF (i ==  1 .AND. iblock == 1 .AND. my_cart_id == 0) THEN
        WRITE(*,'(A       )') '  '
        WRITE(*,'(A,2I5   )') 'TURB-DIAGNOSIS diffusion: iblock = ', iblock, i

        WRITE(*,'(A       )') ' Control Switches and Variables'
        WRITE(*,'(A,I28   )') '   iini             :  ', iini
        WRITE(*,'(A,L28   )') '   ltkeinp          :  ', ltkeinp
        WRITE(*,'(A,L28   )') '   lstfnct          :  ', lstfnct
        WRITE(*,'(A,L28   )') '   l3dturb          :  ', l3dturb
        WRITE(*,'(A,L28   )') '   lrunsso          :  ', lrunsso
        WRITE(*,'(A,L28   )') '   lruncnv          :  ', lruncnv
        WRITE(*,'(A,L28   )') '   lsfluse          :  ', lsfluse
        WRITE(*,'(A,F28.16)') '   dt_tke           :  ', dt_tke
        WRITE(*,'(A,F28.16)') '   dt_var           :  ', dt_var
        WRITE(*,'(A,I28   )') '   nprv             :  ', nprv
        WRITE(*,'(A,I28   )') '   ntur             :  ', ntur
        WRITE(*,'(A,I28   )') '   ntim             :  ', ntim
        WRITE(*,'(A,I28   )') '   nvec             :  ', nvec
        WRITE(*,'(A,I28   )') '   kcm              :  ', kcm
        WRITE(*,'(A       )') ' Other input parameters:'
        WRITE(*,'(A,F28.16)') '   l_hori           :  ', l_hori      (i)
        WRITE(*,'(A,F28.16)') '   hhl    ke        :  ', hhl         (i,ke)
        WRITE(*,'(A,F28.16)') '   hhl    ke1       :  ', hhl         (i,ke1)
        WRITE(*,'(A,F28.16)') '   gz0              :  ', gz0         (i)
        WRITE(*,'(A,F28.16)') '   lpat             :  ', l_pat       (i)
        WRITE(*,'(A,F28.16)') '   t_g              :  ', t_g         (i)
        WRITE(*,'(A,F28.16)') '   qv_s             :  ', qv_s        (i)
        WRITE(*,'(A,F28.16)') '   ps               :  ', ps          (i)
        WRITE(*,'(A,F28.16)') '   u     ke         :  ', u           (i,ke)
        WRITE(*,'(A,F28.16)') '   v     ke         :  ', v           (i,ke)
!       WRITE(*,'(A,F28.16)') '   w     ke         :  ', w           (i,ke)
        WRITE(*,'(A,F28.16)') '   t     ke         :  ', t           (i,ke)
        WRITE(*,'(A,F28.16)') '   qv    ke         :  ', qv          (i,ke)
        WRITE(*,'(A,F28.16)') '   qc    ke         :  ', qc          (i,ke)
        WRITE(*,'(A,F28.16)') '   prs   ke         :  ', prs         (i,ke)
        WRITE(*,'(A,F28.16)') '   tvm              :  ', tvm         (i)
        WRITE(*,'(A,F28.16)') '   tvh              :  ', tvh         (i)
        WRITE(*,'(A,F28.16)') '   tfm              :  ', tfm         (i)
        WRITE(*,'(A,F28.16)') '   tfh              :  ', tfh         (i)
!       WRITE(*,'(A,F28.16)') '   tkred_sfc        :  ', tkred_sfc   (i)
        WRITE(*,'(A,F28.16)') '   l_scal           :  ', l_scal      (i)
        WRITE(*,'(A,F28.16)') '   fc_min           :  ', fc_min      (i)
        WRITE(*,'(A,F28.16)') '   liqs             :  ', liqs        (i,ke1)
        WRITE(*,'(A,F28.16)') '   tke   ke         :  ', tke         (i,ke,ntur)
        WRITE(*,'(A,F28.16)') '   tke   ke1        :  ', tke         (i,ke1,ntur)
        WRITE(*,'(A,F28.16)') '   tkvm  ke         :  ', tkvm        (i,ke)
        WRITE(*,'(A,F28.16)') '   tkvm  ke1        :  ', tkvm        (i,ke1)
        WRITE(*,'(A,F28.16)') '   tkvh  ke         :  ', tkvh        (i,ke)
        WRITE(*,'(A,F28.16)') '   tkvh  ke1        :  ', tkvh        (i,ke1)
        WRITE(*,'(A       )') '  '
      ENDIF
    ENDDO
  ENDIF

!------------------------------------------------------------------------------------
! 0)  Berechnung der Erhaltungsvariablen (auf 'zvari') samt des Bedeckungsgrades
!     und thermodynamischer Hilfgroessen, sowie der turbulenten Laengenskalen:
!------------------------------------------------------------------------------------

  ! Achtung: Bei T-Gleichung in cv-Form gesamte Thermodynamik ueberpruefen auch gegenueber satad

  lcircterm=(pat_len.GT.z0)
  ldotkedif=(c_diff .GT.z0)
  lcircdiff=(lcircterm .AND. imode_calcirc.EQ.2)

  ! Thermodynamische Hilfsvariablen auf Hauptflaechen:
  CALL adjust_satur_equil ( khi=1, ktp=1,                           &
 
           i_st=ivstart, i_en=ivend, k_st=1, k_en=ke, i1dim=nvec,   &
 
           lcalrho=.FALSE., lcalepr=.FALSE.,                        &
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        &
 
           icldmod=icldm_turb,                                      &
 
           zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          &
 
           prs=prs, t=t,     qv=qv,    qc=qc,                       &
 
           psf=ps,                                                  &
 
           exner=epr, rcld=rcld, dens=rhoh,                         &
 
           r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3), g_tet=zaux(:,:,4), &
                                                 g_h2o=zaux(:,:,5), &
 
           tet_l=zvari(:,:,tet_l), q_h2o=zvari(:,:,h2o_g),          &
                                   q_liq=zvari(:,:,liq),            &
           lacc=lzacc )

  ! Thermodynamische Hilfsvariablen auf Unterrand der Prandtl-Schicht:

!DIR$ IVDEP
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
    prss(i,ke1)=ps(i)
    tmps(i,ke1)=t_g(i)
    vaps(i,ke1)=qv_s(i)
  END DO
  !$ACC END PARALLEL

  !$ACC WAIT

  CALL adjust_satur_equil ( khi=ke1, ktp=ke,                        &
 
           i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1, i1dim=nvec,&
 
           lcalrho=.TRUE., lcalepr=.TRUE.,                          &
           lcaltdv=.TRUE., lpotinp=.FALSE., ladjout=.FALSE.,        &
 
           icldmod=icldm_turb,                                      &
 
           zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          &
 
           !Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
           !         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
           prs=prss, t=tmps, qv=vaps, qc=liqs,                      &
 
           fip=tfh,                                                 &
 
           exner=zaux(:,ke1:ke1,1), rcld=rcld(:,ke1:ke1),           &
           dens=rhon(:,ke1:ke1),                                    &
 
           r_cpd=zaux(:,ke1:ke1,2), qst_t=zaux(:,ke1:ke1,3),        &
           g_tet=zaux(:,ke1:ke1,4), g_h2o=zaux(:,ke1:ke1,5),        &
 
           tet_l=zvari(:,ke:ke1,tet_l), q_h2o=zvari(:,ke:ke1,h2o_g),&
                                        q_liq=zvari(:,ke:ke1,liq),  &
           lacc=lzacc )

  ! Note: 
  !     After a proper rearrangement, it should no longer be necessary that surface layer
  !      calculations are done in 'turbdiff'. 

  ! Beachte:
  !     'zvari(:,ke1,tet_l)' und 'zvari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen
  !      am Unterrand der Prandtl-Schicht (zero-level). Die Werte im Niveau 'ke' wurden dabei
  !      zur Interpolation benutzt.
  !     'zaux(:,ke1,1)' enthaelt den Exner-Faktor im zero-level. 
  !     Das Feld 'zaux(:,:,1) wird im Folgenden mit dem Exner-Faktor auf NF belegt.

  ! Kommentar: 
  !     Der 2-te Aufruf von 'adjust_satur_equil' stellt die unteren Randwerte der thermodyn. Variablen
  !      zur Verfuegung. Dies koennte in den 1-ten Aufruf integriert werden, wenn alle thermodyn.
  !      Modell-Variablen bis "k=ke1" allociert waeren. Dies wuere Vieles vereinfachen!

  IF (imode_trancnf.EQ.1) THEN !old version of zero-level-values requested
    !Transformation of Tet_l at zero-level into the value following from the old
    !treatment of interpolation in terms of T_l (rather than Tet_l):

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
       zvari(i,ke1,tet_l) = zvari(i,ke1,tet_l)  &
              + ( (epr(i,ke)-zaux(i,ke1,1))*zvari(i,ke,tet_l)*(z1-tfh(i)) ) &
                    / zaux(i,ke1,1)
    END DO
    !$ACC END PARALLEL
  END IF

  ! Berechnung der horizontalen Windgeschwindigkeiten im Massenzentrum der Gitterbox:
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=1,ke
!DIR$ IVDEP
    DO i=ivstart, ivend
      zvari(i,k,u_m)=u(i,k)
      zvari(i,k,v_m)=v(i,k)
    END DO
  END DO
  !$ACC END PARALLEL

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
         zvari(i,ke1,u_m)=zvari(i,ke,u_m)*(z1-tfm(i))
         zvari(i,ke1,v_m)=zvari(i,ke,v_m)*(z1-tfm(i))
  END DO
  !$ACC END PARALLEL

  ! Berechnung der Schichtdicken und der Dichte auf Nebenflaechen:
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP SEQ
  DO k=1,ke
!DIR$ IVDEP
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      dicke(i,k)=hhl(i,k)-hhl(i,k+1)
    END DO
  END DO
  !$ACC END PARALLEL

  pvar(1)%bl => rcld       ; pvar(1)%ml => rcld  !NF-Werte wieder nach 'rcld'
  pvar(2)%bl => zaux(:,:,1); pvar(2)%ml => epr   !NF-Werte nach 'zaux(:,:,1)', weil
                                                     !'epr' auf HF noch benoetigt wird.
  DO n=2,naux
    pvar(1+n)%bl => zaux(:,:,n) ; pvar(1+n)%ml => zaux(:,:,n)
  END DO

  ! Interpolation der thermodyn. Hilfsgroessen im Feld zaux(),
  ! der Wolkendichte und der Luftdichte auf Nebenflaechen:

  prhon => rhon
  prhoh => rhoh

  CALL bound_level_interp( ivstart, ivend, 2,ke, &
    !-----------------------------------------------------
    !test: mass weighted interpolation
    !  nvars=1, pvar=(/varprf(prhon,prhoh)/), depth=dicke)
    !Achtung: Macht minimale Unterschiede
       nvars=1, pvar=(/varprf(prhon,prhoh)/), depth=dp0, lacc=lzacc)
    !-----------------------------------------------------

  !Note: 
  !Internal order of level looping in 'bound_level_interp' allows to store the 
  !'bl'-values (output) at the same place as the 'ml'-values (input).
      
  CALL bound_level_interp( ivstart, ivend, 2, ke,               &
                           nvars=naux+1, pvar=pvar, depth=dp0, rpdep=hlp, lacc=lzacc)

  !Spezifische effektive Dicke der Prandtlschicht:

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend

    !Achtung: < zum Vergleich mit alter Variante:
    ! velo=MAX( vel_min, SQRT(zvari(i,ke,u_m)**2+zvari(i,ke,v_m)**2) )
    ! tvm(i)=tcm(i)*velo
    ! tvh(i)=tch(i)*velo
    !> zum Vergleich: in ICON werden zwischen turbtran und turbdiff 'u' und 'v' incrementiert!
    !Achtung: Modifikation tcm -> tvm; tch -> tvh: macht Unterschiede
    dzsm(i)=tkvm(i,ke1)*tfm(i)/tvm(i)
    dzsh(i)=tkvh(i,ke1)*tfh(i)/tvh(i)
  END DO
  !$ACC END PARALLEL

  ! Beachte: Auch wenn 'turbtran' nicht als Transfer-Schema benutzt wird, muessen 'tkvm', tkvh'
  !         und 'tke' und (falls es PRESENT ist) auch 'edr' fuer "k=ke1" belegt sein!

  ! Berechnung der turbulenten Laengenscalen:

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
    len_scale(i,ke1)=gz0(i)/grav
  END DO
  !$ACC END PARALLEL

  IF (PRESENT(c_big) .AND. PRESENT(r_air)) THEN

    !US up to now it is kcm = ke+1 and the next vertical loop will not be executed!!
    !   if a canopy layer is implemented, kcm will be <= ke

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=ke,kcm,-1 !Innerhalb des Bestandesmodells
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR PRIVATE(l_turb)
      DO i=ivstart, ivend
        IF (c_big(i,k).GT.z0) THEN
          ! Die turbulente Laengenskala wird durch die Laengenskala 
          ! der lufterfuellten Zwischenraeume limitiert:

          l_turb=z1/(c_big(i,k)*SQRT(z1/EXP(r_air(i,k))-z1))
          len_scale(i,k)=MIN( dicke(i,k)+len_scale(i,k+1), l_turb )
        ELSE
          len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
  ENDIF

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP SEQ
  DO k=kcm-1,1,-1
!DIR$ IVDEP
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      len_scale(i,k)=dicke(i,k)+len_scale(i,k+1)
    END DO
  END DO
  !$ACC END PARALLEL


  ! Uebergang von der maximalen turbulenten Laengenskala zur
  ! effektiven turbulenten Laengenskala:

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP SEQ
  DO k=ke1,1,-1
!DIR$ IVDEP
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      lay(i)=l_scal(i)
      len_scale(i,k)=akt*MAX( len_min, &
                            ! len_scale(i,k)/(z1+len_scale(i,k)/lay(i)) )
                       lay(i)*len_scale(i,k)/(lay(i)+len_scale(i,k)) )
    END DO
  END DO
  !$ACC END PARALLEL

  ! Initialisierung der Felder fuer tke,tkvh,tkvm:

  IF (lini) THEN  !nur beim allerersten Durchgang

    ! Erste Schaetzwerte aus vereinfachtem TKE-Gleichgewicht:
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=2,kem
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR PRIVATE(com_len, edh, fh2, fm2, fakt, lm, lh, val1, val2)
      DO i=ivstart, ivend

        !      der Einfachheit halber nur lokale Berechnung der
        !      vertikalen Gradienten:

        com_len=len_scale(i,k)
        edh=z2/(hhl(i,k+1)-hhl(i,k-1))

        grad(i,u_m  )=(zvari(i,k,u_m  )-zvari(i,k-1,u_m  ))*edh
        grad(i,v_m  )=(zvari(i,k,v_m  )-zvari(i,k-1,v_m  ))*edh
        grad(i,tet_l)=(zvari(i,k,tet_l)-zvari(i,k-1,tet_l))*edh
        grad(i,h2o_g)=(zvari(i,k,h2o_g)-zvari(i,k-1,h2o_g))*edh

        fh2=zaux(i,k,4)*grad(i,tet_l)+zaux(i,k,5)*grad(i,h2o_g)
        fm2=MAX( grad(i,u_m)**2+grad(i,v_m)**2, fc_min(i) )

        ! Vereinfachte Loesung mit Rf=Ri:
        IF (fh2.GE.(z1-rim)*fm2) THEN
          ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
          ! werden durch lm bei der krit. Ri-Zahl angenaehert:
          fakt=z1/rim-z1
          lm=com_len*(sm_0-(a_6+a_3)*fakt)
          lh=lm
        ELSE
          fakt=fh2/(fm2-fh2)
          lm=com_len*(sm_0-(a_6+a_3)*fakt)
          lh=com_len*(sh_0-a_5*fakt)
        END IF

        val1=lm*fm2
        val2=lh*fh2
        hlp(i,k)=MAX( val1-val2, rim*val1 )

        IF (ltkeinp) THEN
          tke(i,k,1)=tke(i,k,ntur)
        ELSE
          tke(i,k,1)=MAX( vel_min, SQRT(d_m*com_len*hlp(i,k)) ) !Initialwert fuer SQRT(2TKE)
        END IF

        val1=MAX ( con_m, tkmmin ); tkvm(i,k)=lm*tke(i,k,1)
        val2=MAX ( con_h, tkhmin ); tkvh(i,k)=lh*tke(i,k,1)
        IF (imode_tkemini.EQ.2) THEN
          tke(i,k,1)=tke(i,k,1)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                         val2/tkvh(i,k) )
        END IF

        ! Achtung: Bislang fehlte die Beschraenkung: Macht Unterschiede geg. ICON
        tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv
        tkvh(i,k)=MAX( val2, tkvh(i,k) )

        ! Am Anfang konnte noch keine Diffusion von q=SQRT(2*TKE) berechnet werden:
        tketens(i,k)=z0

      END DO
    END DO
    !$ACC END PARALLEL

!DIR$ IVDEP
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      tke(i,1,1)=tke(i,2,1)
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(3)
    DO n=2,ntim
      DO k=1,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          tke(i,k,n)=tke(i,k,1)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT

  END IF   ! IF (lini)

!------------------------------------------------------------------------------------
! 1)  Berechnung der benoetigten vertikalen Gradienten und Abspeichern auf 'zvari':
!------------------------------------------------------------------------------------

  ! Am unteren Modellrand:

!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
    lays(i,mom)=z1/dzsm(i)
    lays(i,sca)=z1/dzsh(i)
  END DO
  !$ACC END PARALLEL

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO n=1,nmvar
!DIR$ IVDEP
    DO i=ivstart, ivend
      zvari(i,ke1,n)=(zvari(i,ke,n)-zvari(i,ke1,n))*lays(i,ivtp(n))
    END DO
  END DO
  !$ACC END PARALLEL

  ! An den darueberliegenden Nebenflaechen:
  IF (lnonloc) THEN   ! nonlocal calculation of vertical gradients used for turb. diff.

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      hlp(i,ke1)=z0
    END DO
    !$ACC END PARALLEL

    ! Berechnung nicht-lokaler Gradienten:
    DO n=1,nmvar  ! nmvar = 5

      ! Berechnung der vertikalen Integralfunktionen in hlp():

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=ke,2,-1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          hlp(i,k)=hlp(i,k+1)+zvari(i,k,n)*(hhl(i,k)-hhl(i,k+1))
        END DO
      END DO
      !$ACC END PARALLEL
 
!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(h, hu, kk, k1, k2, ku, wert)
      DO i=ivstart, ivend
        k1=1
        k2=2
        lays(i,k1)=hhl(i,1)-hhl(i,2)

        !$ACC LOOP SEQ
        DO k=2,ke
          ! Berechnung der nicht-lokalen Gradienten als mittlere
          ! Differenzenquotienten ueber die stabilitaetsabhaengige
          ! Laengenskala in tkvh() bzw tkvm():

          ! Speichern der stab.abh. Laengenskala unter lay():

          IF (n.LE.nvel) THEN
            lay(i)=tkvm(i,k)/tke(i,k,nvor)
          ELSE
            lay(i)=tkvh(i,k)/tke(i,k,nvor)
          END IF

          ! Bestimmung der nicht-lokalen Gradienten und
          ! Zwischenspeichern derselben auf dem Feld dicke():
          lays(i,k2)=hhl(i,k)-hhl(i,k+1)

          IF ( lay(i) <= z1d2 * MIN (lays(i,k1), lays(i,k2)) ) THEN

            ! Die vertikalen Diffusionswege schneiden weder
            ! eine untere noch eine obere Hauptflaeche. Es
            ! koennen somit die lokalen Gradienten genommen
            ! werden. Bei sehr kleinen Diffusionslaengen, muss
            ! aus num. Gruenden sogar die lokale Berechnung
            ! gewaehlt werden:

            dicke(i,k)=z2*(zvari(i,k-1,n)-zvari(i,k,n)) / (hhl(i,k-1)-hhl(i,k+1))

          ELSE

            ! Berechn. der benoetigten Referenzhoehen und -level:
            h=hhl(i,k)
            hu=MAX( h-lay(i), hhl(i,ke1) )
            hig(i,1)=hu+lay(i)
            hig(i,2)=h+lay(i)

            kk=k
            DO WHILE (hhl(i,kk).GT.hu)
              kk=kk+1
            END DO
            ku=kk
            DO ii=1,2
              IF (kk.GT.1) THEN
                DO WHILE (hhl(i,kk).LE.hig(i,ii))
                  kk=kk-1
                END DO
              END IF
              lev(i,ii)=kk+1
            END DO

            ! Berechnung der gemittelten Differenzenquotienten
            ! als Ausdruck fuer die nicht-lokalen Gradienten:
            wert=hlp(i,ku)-hlp(i,k) &
                +hlp(i,lev(i,2))-hlp(i,lev(i,1)) &
                +zvari(i,ku-1,n)*(hu-hhl(i,ku)) &
                -zvari(i,lev(i,1)-1,n)*(hig(i,1)-hhl(i,lev(i,1)))&
                +zvari(i,lev(i,2)-1,n)*(hig(i,2)-hhl(i,lev(i,2)))

            dicke(i,k)=wert/(lay(i)*(h-hu))
          END IF

        END DO   ! vertical loop

        kk=k1
        k1=k2
        k2=kk

      END DO     ! horizontal loop
      !$ACC END PARALLEL

      ! Sichern der nicht-lokalen Gradienten im Feld zvari():
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      DO k=2,ke
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          zvari(i,k,n)=dicke(i,k)
        END DO
      END DO
      !$ACC END PARALLEL

    END DO   ! loop over n=1,nmvar

    ! Belegung von dicke() mit den Schichtdicken*rhon/dt_tke
    ! bzgl. Nebenflaechen:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=2,ke
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        dicke(i,k)=rhon(i,k)*z1d2*(hhl(i,k-1)-hhl(i,k+1))*fr_tke
      END DO
    END DO
    !$ACC END PARALLEL

  ELSE    ! lnonloc

    ! Berechnung lokaler Gradienten:
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=ke,2,-1
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR PRIVATE(com_len)
      DO i=ivstart, ivend
        com_len=(hhl(i,k-1)-hhl(i,k+1))*z1d2
        hlp(i,k)=z1/com_len
        dicke(i,k)=rhon(i,k)*com_len*fr_tke
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    DO n=1,nmvar

#ifdef __INTEL_COMPILER
      FORALL(k=2:ke,i=ivstart:ivend)                        &
               zvari(i,k,n)=(zvari(i,k-1,n)-zvari(i,k,n))*hlp(i,k)
#else
      !$ACC LOOP SEQ
      DO k=ke,2,-1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          zvari(i,k,n)=(zvari(i,k-1,n)-zvari(i,k,n))*hlp(i,k)
        END DO
      END DO
#endif

    END DO
    !$ACC END PARALLEL

    !$ACC WAIT

  END IF    ! lnonloc

!------------------------------------------------------------------------------------
! 2)  Berechnung der verallgemeinerten Antriebsfunktionen einschliesslich der
!     Korrekturen innerhalb der Rauhigkeitsschicht (samt der Windtendenz durch Formreibung)
!     und der Scherung durch nicht-turbulente subskalige Stroemungen:
!------------------------------------------------------------------------------------

  ! Thermal forcing:
  !-------------------------

  ! Achtung:
  !'frh'(ke1) wird fuer Zirkulationsterm und Temperaturkorrektur benoetigt

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=2,ke1 
!DIR$ IVDEP
    DO i=ivstart, ivend
      frh(i,k)=zaux(i,k,4)*zvari(i,k,tet_l) + zaux(i,k,5)*zvari(i,k,h2o_g)
    END DO
  END DO
  !$ACC END PARALLEL

  !Note: zaux(:,:,4) and zaux(:,:,5) are free now.

  ! Total mechanical forcing:  
  !--------------------------

  !hdef2 = (d1v2+d2v1)^2 + (d1v1-d2v2)^2 !horizontal deformation square     (at half levels)
  !hdiv  = (d1v1+d2v2)                   !horizontal wind-divergence            ,,

  !dwdx !zonal      derivation of vertical wind                                 ,,
  !dwdy !meridional derivation of vertical wind                                 ,,

  !vel_div=hdiv+dzdw=0 !Incomressibility

  !itype_sher = 0 : only single column vertical shear
  !             1 : previous and additional 3D horiz. shear correction
  !             2 : previous and additional 3D vertc. shear correction

  !ltkeshs: consider separated non-turbulent horizontal shear mode for TKE forcing
  !loutshs: consider separated non-turbulent horizontal shear mode for output

  ! Mechanical forcing by vertical shear:

  IF (itype_sher.EQ.2 .AND. (PRESENT(dwdx) .AND. PRESENT(dwdy) .AND. PRESENT(hdiv))) THEN

    !Include 3D-shear correction by the vertical wind (employing incomressibility):

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=2,kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        frm(i,k)=MAX( (zvari(i,k,u_m)+dwdx(i,k))**2+(zvari(i,k,v_m)+dwdy(i,k))**2 &
                                                         +z3*hdiv(i,k)**2, fc_min(i) )
      END DO
    END DO
    !$ACC END PARALLEL

  ELSE   

    !Load pure single column shear:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    DO k=2,kem 
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        frm(i,k)=MAX( zvari(i,k,u_m)**2+zvari(i,k,v_m)**2, fc_min(i))
      END DO
    END DO
    !$ACC END PARALLEL

  END IF    

  ! Mechanical forcing by horizontal shear:

  IF (PRESENT(hdef2)) THEN 
    IF (itype_sher.GE.1) THEN   !Apply horizontal 3D-shear correction:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,kem 
!DIR$ IVDEP
        DO i=ivstart, ivend
          frm(i,k)=frm(i,k)+hdef2(i,k) !extended shear
        END DO
      END DO
      !$ACC END PARALLEL
    END IF   
  END IF

  IF (lssintact) THEN !save pure turbulent shear
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k=2,kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        ftm(i,k)=frm(i,k)
      END DO
    END DO
    !$ACC END PARALLEL
  END IF

  !For_Tuning>
  !     Preparation for Richardson-number-dependent factor used for correcting 
  !     the minimum diffusion coefficient and the horizontal shear production term:

  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=2,ke
!DIR$ IVDEP
    DO i=ivstart, ivend
      !Achtung: Mit Hilfe von 'frm' und 'frh' auszudruecken: <
      !x1 = z1d2*(hhl(i,k-1)-hhl(i,k+1))
      !x2 = MAX( 1.e-6_wp, ((u(i,k-1)-u(i,k))**2+(v(i,k-1)-v(i,k))**2)/x1**2 )          ! |du/dz|**2
      !x3 = MAX( 1.e-5_wp, grav/(z1d2*(t(i,k-1)+t(i,k)))*((t(i,k-1)-t(i,k))/x1+tet_g) ) !       N**2

      !xri(i,k) = EXP(z2d3*LOG(x2/x3))  ! 1/Ri**(2/3)

      xri(i,k)=EXP( z2d3*LOG( MAX( 1.e-6_wp, frm(i,k) ) / & !1/Ri**(2/3)
                              MAX( 1.e-5_wp, frh(i,k) ) ) )
      !>
    END DO
  END DO
  !$ACC END PARALLEL

  !>For_Tuning
  IF (PRESENT(hdef2)) THEN

    !Additional impact by separated horizontal shear:
    IF ((ltkeshs .OR. (loutshs .AND. PRESENT(tket_hshr))) .AND. PRESENT(hdiv)) THEN 
      !Include separated horizontal shear mode:

      fakt=z1/(z2*sm_0)**2; wert=a_hshr*akt*z1d2

!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        lay(i)=wert*l_hori(i) !uncorrected effective horizontal length scale
      END DO
      !$ACC END PARALLEL

      IF (imode_shshear.EQ.2) THEN
        !>Tuning

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(x4, x4i)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            ! Factor for variable 3D horizontal-vertical length scale proportional to 1/SQRT(Ri),
            ! decreasing to zero in the lowest two kilometer above ground

            x4i = MIN( 1._wp, 0.5e-3_wp*(hhl(i,k)-hhl(i,ke1)) )
            x4 = 3._wp*x4i**2 - 2._wp*x4i**3                   ! low-level reduction factor
            hor_scale(i,k) = lay(i)*MIN( 5.0_wp, MAX( 0.01_wp, x4*xri(i,k) ) )/MAX(1._wp,0.2_wp*tke(i,k,1))
          END DO
        END DO
        !$ACC END PARALLEL

        !>Tuning: This kind of correction can be substituded by a less ad-hoc approach.

      ELSE

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            hor_scale(i,k) = lay(i)
          END DO
        END DO
        !$ACC END PARALLEL

      ENDIF

      !strain velocity (shv) of the separated horizontal shear mode:
      IF (imode_shshear.EQ.0) THEN !former variant based on 3D-shear and incompressibility
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            shv(i,k)=hor_scale(i,k)*SQRT(hdef2(i,k)+hdiv(i,k)**2) !not equal to trace of 2D-strain tensor
          END DO
        END DO
        !$ACC END PARALLEL
      ELSE !new variant in accordance with the trace constraint for the separated horizontal strain tensor
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(wert)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            wert=fakt*hdiv(i,k)
            shv(i,k)=hor_scale(i,k)*(SQRT(wert**2+hdef2(i,k))-wert) !equal to trace of 2D-strain tensor
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,kem
!DIR$ IVDEP
        DO i=ivstart, ivend
          hlp(i,k)=(shv(i,k))**3/hor_scale(i,k) !additional TKE-source by related shear
        END DO
      END DO
      !$ACC END PARALLEL

      IF (loutshs .AND. PRESENT(tket_hshr)) THEN
        !Load output variable for the TKE-source by separated horiz. shear:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem 
!DIR$ IVDEP
          DO i=ivstart, ivend
            tket_hshr(i,k)=hlp(i,k)
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      IF (ltkeshs) THEN 
        !Consider separated horizontal shear mode in mechanical forcing:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem 
!DIR$ IVDEP
          DO i=ivstart, ivend
            frm(i,k)=frm(i,k)+hlp(i,k)/tkvm(i,k) !extended shear
          END DO
        END DO
        !$ACC END PARALLEL

        IF (l3dturb) THEN
          ! Load related horizontal diffusion coefficients:
          fakt=sh_0/sm_0
          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO k=2,kem
!DIR$ IVDEP
            DO i=ivstart, ivend
              tkhm(i,k)=hor_scale(i,k)*shv(i,k) !for momentum related to the sep. shear mode
              tkhh(i,k)=fakt*tkhm(i,k)          !for scalars    ,,       ,,            ,,
            END DO
          END DO
          !$ACC END PARALLEL
        END IF  ! l3dturb

      END IF    ! ltkeshs

    END IF      ! IF ((ltkeshs .OR. (loutshs .AND. PRESENT(tket_hshr))) .AND. PRESENT(hdiv))
 
  END IF        ! IF (PRESENT (hdef2))

  !Addition verallgemeinerter Scherterme durch nicht-turbulente subskalige Stroemung:

  IF (iini.NE.1) THEN !nicht bei der separaten Initialisierung

    !Special data regions for optional variables
        
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=2,kem
      !$ACC LOOP GANG VECTOR PRIVATE(vel1, vel2)
!DIR$ IVDEP
      DO i=ivstart, ivend

        IF (lrunsso .AND. PRESENT(ut_sso) .AND. PRESENT(vt_sso)) THEN
          !SSO-Schema ist aktiv und SSO-Tendenzen des Windes sind vorhanden:

          !Berechnung der TKE-Tendenz durch Nachlaufwirbelproduktion aus SSO-Tendenzen:

          !Achtung: horizontale Pos. beachten!
          vel1=-(ut_sso(i,k)  *u(i,k  )+vt_sso(i,k)  *v(i,k  ))*dp0(i,k-1)
          vel2=-(ut_sso(i,k-1)*u(i,k-1)+vt_sso(i,k-1)*v(i,k-1))*dp0(i,k)

          src(i)=MAX( z0, (vel1+vel2)/(dp0(i,k-1)+dp0(i,k)) )

          !Beachte:
          !Die SSO-Tendenzen beziehen sich tatsaechlich auf Massenpunkte, sie werden
          !erst spaeter in SUB 'organize_physics' auf die Zwischenpositionen interpoliert!
          !Obwohl vel1 und vel2 immer positiv sein muessten, wird zur Sicherheit MAX benutzt!

          IF (PRESENT(tket_sso)) THEN
            tket_sso(i,k)=src(i)
          END IF

          !Addition des Scherterms durch Nachlaufwirbel:

          IF (imode_tkesso == 1) THEN !Nachlaufwirbeltendenzen sollen beruecksichtigt werden
            frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)
          ELSE IF (imode_tkesso == 2) THEN ! Reduce TKE production in the presence of large Richardson numbers
            frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)*MIN(1.0_wp,MAX(0.01_wp,xri(i,k)))
          ELSE IF (imode_tkesso == 3) THEN ! Reduce TKE production in the presence of large Richardson numbers
            frm(i,k)=frm(i,k) + src(i)/tkvm(i,k)*MIN(1.0_wp,MAX(0.01_wp,xri(i,k)))*MIN(1.0_wp,l_hori(i)/2000._wp)
          END IF

        END IF

        IF (lruncnv .AND. ltkecon .AND. PRESENT(tket_conv)) THEN
          !Konvektionsschema ist aktiv, soll mit Turbulenz interagieren und conv. TKE-Tend. ist vorhanden:

          !Addition des Scherterms durch konvektive Zirkulation:
          frm(i,k)=frm(i,k) + MAX( z0, tket_conv(i,k)/tkvm(i,k) )
          !Beachte:  Obwohl tket_conv immer positiv sein muesste, wird zur Sicherheit MAX benutzt!
        END IF
      END DO
    END DO
    !$ACC END PARALLEL

  END IF  ! IF (iini /= 1)

  ! Berechnung von Korrekturtermen innerhalb der Rauhigkeitsschicht
  ! (ausser Volumenterme, die zur Diffusion gehoeren):

  IF (PRESENT(c_big) .AND. PRESENT(c_sml)) THEN

    !US: at the moment kcm = kem+1, so this vertical loop is never executed
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=kcm,kem !von oben nach unten durch Rauhiggkeitsschicht
      !$ACC LOOP GANG VECTOR PRIVATE(velo, wert)
!DIR$ IVDEP
      DO i=ivstart, ivend
        ! Achtung: Neue Behandlung der Rauhigkeitsschicht einfuehren
        ! Vertikalwind und Formreibungskoeff. auf Hauptflaechen:
        lay(i)=z1d2*(w(i,k)+w(i,k+1))
        src(i)=z1d2*(c_big(i,k) + c_sml(i,k) + c_big(i,k+1) + c_sml(i,k+1))

        ! Windbetrag auf der angrenzenden Hauptflaeche:
        IF (k.EQ.kcm) THEN
          velo=z1d2*(w(i,k-1)+w(i,k))
          lays(i,1)=SQRT(u(i,k-1)**2+v(i,k-1)**2+velo**2)
        END IF

        ! Reduzierte Konstanten und implizite Horizontalwind-Tendenzen durch Formreibung:

        ! Windbetrag auf der aktuellen Hauptflaeche:
        lays(i,2)=SQRT(u(i,k)**2+v(i,k)**2+lay(i)**2)

        ! Aufaddieren der Windtendenzen durch Fromreibung:
        wert=src(i)*lays(i,2)

        utens(i,k)=utens(i,k)-tinc(u_m)*wert*u(i,k)/(z1+dt_var*wert)
        vtens(i,k)=vtens(i,k)-tinc(v_m)*wert*v(i,k)/(z1+dt_var*wert)

        ! Windbetrag auf Nebenflaechen:
        can(i,k)=(lays(i,2)*dp0(i,k-1)+lays(i,1)*dp0(i,k))/(dp0(i,k-1)+dp0(i,k))

        ! Windbetrag unter Wirkung der Formreibung:
        velo=can(i,k)/(z1+can(i,k)*(c_big(i,k)+c_sml(i,k))*dt_var)

        ! Addition des Scherterms durch Nachlaufwirbel an grossen Rauhigkeitselementen:
        frm(i,k)=frm(i,k)+c_big(i,k)*velo**3/tkvm(i,k)

        ! Frequenz der kleinskaligen Rauhigkeitselemente:
        can(i,k)=c_sml(i,k)*can(i,k)
      END DO
    END DO  !k=kcm,kem !von oben nach unten druch Rauhigkeitsschicht
    !$ACC END PARALLEL

  ENDIF   ! IF (PRESENT(c_big) .AND. PRESENT(c_sml))

  ! optional output of source terms
  IF(PRESENT(tket_therm).AND.PRESENT(tket_mech).AND.PRESENT(tket_shear)) THEN
     !$ACC PARALLEL ASYNC(1) IF(lzacc)
     DO k=2, ke1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend !istartpar,iendpar
           tket_therm(i,k) = frh(i,k) * tkvh(i,k)
           tket_mech(i,k)  = frm(i,k) * tkvm(i,k)
           tket_shear(i,k) = ftm(i,k) * tkvm(i,k)
        END DO
     END DO
  !$ACC END PARALLEL
  END IF

  !$ACC WAIT

  ! Optionale vertikale Glaettung des mechanischen Antriebs:
  IF (lcalc_frcsmot) THEN
    CALL vert_smooth (i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
                      disc_mom=dicke, cur_tend=frm, vertsmot=frcsmot, smotfac=trop_mask, lacc=lzacc )
  END IF

  ! Optionale vertikale Glaettung des thermischen Antriebs:
  IF (lcalc_frcsmot) THEN
    CALL vert_smooth (i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
                      disc_mom=dicke, cur_tend=frh, vertsmot=frcsmot, smotfac=trop_mask, lacc=lzacc )
  END IF

  ! Belegung von tkvh und tkvm mit den stabilitaetsabhaengigen Laengenmassen:
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR COLLAPSE(2)
  DO k=2,kem
!DIR$ IVDEP
    DO i=ivstart, ivend
      tkvh(i,k)=tkvh(i,k)/tke(i,k,nvor)
      tkvm(i,k)=tkvm(i,k)/tke(i,k,nvor)
    END DO
  END DO
  !$ACC END PARALLEL

!------------------------------------------------------------------------------------
! 3)  Loesung der turbulenten Bilanzgleichungen (Bestimmung von TKE und der Stabilitaetsfunktionen)
!     und Berechnung der turbulenten Diffusionskoeffizienten:
!------------------------------------------------------------------------------------

  DO it_durch=it_start, it_end

    !Die Schleife wird nur bei der Initialisierung (d.h. beim ersten Aufruf) wiederholt,
    !um TKE-Gleichgewicht anzunaehern. Die resultierenden TKE-Werte der Zeitstufe 'ntur'
    !gehoeren in diesem Fall dann zur Zeitstufe der uebrigen prognostischen Variablen
    !('nold' bei "leap-frog" oder 'nnow' bei 2-Zeitebenen).
    !Fuer die folgenden Aufrufe wird die Schleife nur einmal durchlaufen und liefert TKE-Werte
    !die gegenueber den Vorgaengerwerten um einen Zeitschritt weiter in der Zukunft liegen,
    !also wieder zur Zeitstufe der uebrigen prognostischen Variablen gehoeren.

    !$ACC WAIT
    
    CALL solve_turb_budgets (khi=1, it_s=it_durch, it_start=it_start,                 &
                             i_st=ivstart, i_en=ivend, k_st=2, k_en=kem,              &
                             i1dim=nvec, kcm=kcm, ntur=ntur, nvor=nvor,               &
                             lssintact=lssintact, lupfrclim=.FALSE.,                  &
                             lpresedr=PRESENT(edr), lstfnct=lstfnct, ltkeinp=ltkeinp, &
                             imode_stke=imode_turb, imode_vel_min=1,                  &
                             dt_tke=dt_tke, fr_tke=fr_tke,                            &
                             tke=tke, ediss=ediss,                                    &
                             fm2=frm, fh2=frh, ft2=ftm, lsm=tkvm, lsh=tkvh,           &
#ifdef SCLM
                             grd=zvari,                                               &
#endif
                             fcd=can, tls=len_scale, tvt=tketens, avt=tketadv,        &
                             lacc=lzacc)

    IF (it_durch.LT.it_end .AND. .NOT.ltkeinp) THEN
      nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
    END IF

  END DO

  ! Kein TKE-Gradient am Oberrand:
!DIR$ IVDEP
  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
  !$ACC LOOP GANG VECTOR
  DO i=ivstart, ivend
    tke(i,1,ntur)=tke(i,2,ntur)
  END DO
  !$ACC END PARALLEL

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  IF (iini.EQ.1) THEN !only for separate initialization before the time loop

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(val1, val2)
    DO k=2, kem
!DIR$ IVDEP
      DO i=ivstart, ivend
        !Achtung: Bislang fehtlte die laminare Beschraenkung
        val1=MAX ( con_m, tkmmin ); tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
        val2=MAX ( con_h, tkhmin ); tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
        IF (imode_tkemini.EQ.2) THEN
                  tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                                       val2/tkvh(i,k) )
        END IF
        tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv
        tkvh(i,k)=MAX( val2, tkvh(i,k) )
      END DO
    END DO
    !$ACC END PARALLEL

  ELSE

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!US the rest of the program is NOT executed for the initialization phase !!

!------------------------------------------------------------------------------------
!  4) Berechnung der effektiven turbulenten Vertikalgradienten,
!     Standardabweichnung des Saettigungsdefizites und der Diffusionskoeffizienten:
!------------------------------------------------------------------------------------

    IF (ltmpcor) THEN

      !  Berechnung des vert. Temp.grad. fuer den Phasendiffusionsterm:

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          frm(i,k)=zaux(i,k,1)*zvari(i,k,tet_l)-tet_g  !vertical temperature gradient
        END DO
      END DO

      !$ACC END PARALLEL
      IF (icldm_turb.NE.-1) THEN !water phase changes are possible
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2, ke1
!DIR$ IVDEP
          DO i=ivstart, ivend
            frm(i,k)=frm(i,k)+lhocp*zvari(i,k,liq) !liquid water correction
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      !  Dies geschieht schon hier, weil im naechsten Schritt das Feld zvari()
      !  durch die effiktiven Gradienten ueberschrieben wird.
    END IF

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG WORKER VECTOR COLLAPSE(2)
    DO k=2, ke1
!DIR$ IVDEP
      DO i=ivstart, ivend
        zaux(i,k,5)=zaux(i,k,1)*zaux(i,k,3)                     !epr*d_qsat/d_T
        zaux(i,k,4)=c_scld*rcld(i,k)/(z1+rcld(i,k)*(c_scld-z1)) !effective cloud cover
        rcld(i,k)=SQRT(len_scale(i,k)*tkvh(i,k)*d_h)* &         !standard deviation
                       ABS(zaux(i,k,5)*zvari(i,k,tet_l) &       !of local
                          -zvari(i,k,h2o_g))                    !oversaturation
      END DO   
    END DO   
    !$ACC END PARALLEL

    IF (icldm_turb.NE.-1) THEN !consideration of water phase changes

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG WORKER VECTOR COLLAPSE(2) PRIVATE(flw_h2o_g, flw_tet_l)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          ! Effective vertical gradient of liquid water content:
   
          flw_h2o_g=zaux(i,k,4)/(z1+lhocp*zaux(i,k,3))  !weight of h2o_g-flux
          flw_tet_l=-flw_h2o_g*zaux(i,k,5)              !weight of tet_l-flux
   
          zvari(i,k,liq)=  flw_tet_l*zvari(i,k,tet_l) & ! eff_grad(liq)
                        +  flw_h2o_g*zvari(i,k,h2o_g)   
   
          ! Effective vertical gradient of water vapor content and pot. temper.:
   
          zvari(i,k,vap)= zvari(i,k,h2o_g)-zvari(i,k,liq)                      ! eff_grad(vap)
          zvari(i,k,tet)=(zvari(i,k,tet_l)+zvari(i,k,liq)*lhocp/zaux(i,k,1)) & ! eff_grad(tet) 
                           *zaux(i,k,2)                                        !*(Cp/Cpd)

            ! Note: 
            ! -flux(h2o_g)/(rho*K)=grad(h2o_g)=eff_grad(vap)+eff_grad(liq)    
            ! -flux(tet_l)/(rho*K)=grad(tet_l)=eff_grad(tet)-eff_grad(liq)*lh_v/(cp_d*epr)
            ! Treating "lh_v/(cp_d*epr)" like a constnat, besides numerical effects, 
            !  vertical gradient diffusion of non conserved variables without a moist gradient correction
            !  does not change the resulting diffusion of conserved variables. Thus the redistribution of
            !  water constituents can (or should) be left to a final (sub grid scale) saturation adjustment.
        END DO
      END DO
      !$ACC END PARALLEL

    ELSE !no water phase change possible

      !'zvari(:,:,liq)' bleibt unveraendert, genauso wie auch
      !'zvari(:,:,tet)'='zvari(:,:,tet_l)' und 'zvari(:,:,vap)'='zvari(:,:,h2o_g)'.
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          zvari(i,k,tet)=zvari(i,k,tet_l)*zaux(i,k,2) ! grad(tet)*(Cp/Cpd)
        END DO
      END DO
      !$ACC END PARALLEL

    END IF   ! IF (icldm_turb /= -1)

    !Beachte:
    !Die eff_grad ergeben multipliziert mit 'rho*tkvh' die vertikalen Flussdichten
    ! von 'tet', 'vap' und 'liq' unter Beruecksichtigung der turbulenten Phasenuebergaenge.
    !'zvari(:,:,tet)' ist der in der Teta-Gleichung benoetigte effective Gradient 
    ! und entspricht der 'tet'-Flussdichte*(Cp/Cpd).
    !Die Indices 'tet', 'tem' und 'tet_l' haben den gleichen Wert
    ! so wie auch 'vap' und 'h2o_g' den gleichen Wert haben.
    !Die Unterscheidungen sollen nur den jeweiligen Inhalt verdeutlichen.
    !Im Falle "icldm_turb.EQ.-1" verschwindet der eff. Bedeckungsgrad "zaux(:,:,4)=0",
    ! troltzdem steht auch dann die Standardabweichung der Uebersaettigung in 'rcld'.

    ! Beschraenkung der Diffusionskoeffizienten:

    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(fakt, fakt2, val1, val2, x4, x4i)
    DO k=2, ke
!DIR$ IVDEP
      DO i=ivstart, ivend
        ! Berechn. der Diffusionskoeffizienten:
        val1=tkmmin; val2=tkhmin !default minimum diffusion values

        IF (imode_tkvmini.EQ.2) THEN
        !>Tuning
          ! Factor for variable minimum diffusion coefficient proportional to 1/SQRT(Ri);
          ! the namelist parameters tkhmin/tkmmin specify the value for Ri=1:
          IF (gz0(i) < 0.01_wp .AND. l_pat(i) > 0._wp) THEN ! glaciers
            fakt =MIN( z1, tkred_sfc(i)*4.e-3_wp*(hhl(i,k)-hhl(i,ke1)) ) !low-level red.-fact.
            fakt2=MIN( z1, tkred_sfc_h(i)*tkred_sfc(i)*4.e-3_wp*(hhl(i,k)-hhl(i,ke1)) ) !low-level red.-fact.
          ELSE
            fakt =MIN( z1, tkred_sfc(i)*(0.25_wp+7.5e-3_wp*(hhl(i,k)-hhl(i,ke1))) ) !low-level red.-fact.
            fakt2=MIN( z1, tkred_sfc_h(i)*tkred_sfc(i)*(0.25_wp+7.5e-3_wp*(hhl(i,k)-hhl(i,ke1))) ) !low-level red.-fact.
          ENDIF
          fakt =MIN( 2.5_wp, MAX( 0.01_wp, fakt*xri(i,k) ) )
          fakt2=MIN( 2.5_wp, MAX( 0.01_wp, fakt2*xri(i,k) ) )

          val1=tkmmin*fakt; val2=tkhmin*fakt2

          IF (tkhmin_strat.GT.z0 .OR. tkmmin_strat.GT.z0) THEN
            IF (PRESENT(innertrop_mask)) THEN
              ! from ICON Version 180206
              ! Enhanced diffusion in the stratosphere - needed primarily for momentum because 
              ! there is otherwise too little dynamic coupling between adjacent model levels
              fakt = MIN( z1, 2.e-4_wp*MAX( z0, hhl(i,k) - 12500._wp ) ) ! transition zone between 12.5 and 17.5 km

              ! Wider transition zone in the tropics in order to avoid too strong diffusion in the tropopause region
              x4  = z1-z1d3*trop_mask(i)*MIN(z1, 2.e-4_wp*MAX(z0, 22500._wp-hhl(i,k)) )
              x4i = z1-z2d3*innertrop_mask(i)*MIN(z1, 2.e-4_wp*MAX(z0, 27500._wp-hhl(i,k)) )
              fakt = fakt*MIN( x4*1.5_wp, MAX( 0.25_wp, SQRT(xri(i,k)) ) )
              val1=MAX( val1, tkmmin_strat*MIN(x4,x4i)*fakt ) ; val2=MAX( val2, tkhmin_strat*x4*fakt )
            ELSE
              ! necessary for COSMO
              ! Enhanced diffusion in the stratosphere - very important for the data assimilation cycle,
              ! but can also be used in forecasting mode because there is no detectable detrimental
              ! impact on gravity waves:
              fakt = MIN( z1, 2.e-4_wp*MAX( z0, hhl(i,k) - 25000._wp ) ) !lin. incr. betw. 25 and 30 km
              fakt = fakt*MIN( 7.5_wp, MAX( 0.125_wp, xri(i,k) ) )

              val1=MAX( val1, tkmmin_strat*fakt ) ; val2=MAX( val2, tkhmin_strat*fakt )
            ENDIF
          END IF
        !>Tuning: This kind of correction can be substituded by a less ad-hoc approach.
        ! Remark (GZ): The enhanced stratospheric diffusion seems to parameterize a missing 
        !              process outside the turbulence scheme, maybe momentum transports due 
        !              to non-stationary gravity waves. This may also explain why we need a 
        !              much larger minimum diffusion coefficient for momentum than for heat.

        END IF   !  IF (imode_tkvmini.EQ.2)

        !Achtung: Beschraenkung mit lam. diff.coef. fehlte bislang auch in ICON; macht ev. Unterschiede
        val1=MAX ( con_m, val1 ); tkvm(i,k)=tkvm(i,k)*tke(i,k,ntur)
        val2=MAX ( con_h, val2 ); tkvh(i,k)=tkvh(i,k)*tke(i,k,ntur)
        IF (imode_tkemini.EQ.2) THEN
          tke(i,k,ntur)=tke(i,k,ntur)*MAX( z1, val1/tkvm(i,k), & !adapted tke
                                               val2/tkvh(i,k) )
        END IF
        tkvm(i,k)=MAX( val1, tkvm(i,k) ) !corrected tkv's
        tkvh(i,k)=MAX( val2, tkvh(i,k) )

        !test: ohne tk?min (wird aber bei Diffusion benutzt!)
        !test: ohne tk?min und ohen lam. Diffkof. (wird aber bei Diffusion benutzt!)

        ! tkvh und tkvm enthalten jetzt nicht mehr Diffusionslaengen
        ! sondern wieder Diffusionskoeffizienten in m^2/s!

      END DO  ! i
    END DO    ! k
    !$ACC END PARALLEL

    IF (l3dturb) THEN
      !Consider horizontal diffusion coefficients:
      IF (PRESENT(hdef2) .AND. PRESENT(hdiv) .AND. ltkeshs) THEN
        !Add isotropic turbulent part to that part due to the sep. horiz. shear mode:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            tkhm(i,k)=tkhm(i,k)+tkvm(i,k)
            tkhh(i,k)=tkhh(i,k)+tkvh(i,k)
          END DO
        END DO
        !$ACC END PARALLEL
      ELSE !no treatment of sep. horiz. shear mode has taken place
        !Load only the isotropic turbulent part:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,kem
!DIR$ IVDEP
          DO i=ivstart, ivend
            tkhm(i,k)=tkvm(i,k); tkhh(i,k)=tkvh(i,k)
          END DO
        END DO
        !$ACC END PARALLEL
      END IF
    END IF


!------------------------------------------------------------------------------------
! 5) Untere Randwerte der Turbulenzgroessen:
!------------------------------------------------------------------------------------

    !Beachte: Auch wenn 'turbtran' nicht als Transfer-Schema benutzt wird,
    !         muessen 'tkvm', tkvh' und 'tke' und (falls es PRESENT ist)
    !         auch 'edr' fuer "k=ke1" belegt sein!

    IF (lsfli(tem)) THEN !use explicit shfl_s
!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        zvari(i,ke1,tet)=shfl_s(i)/(cp_d*rhon(i,ke1)*tkvh(i,ke1)*zaux(i,ke1,1))
      END DO
      !$ACC END PARALLEL

      !Note:'zvari(tet)' belongs to potential temperature, and
      !        'zaux(1)' contains Exner factor on half levels.
    END IF

    IF (lsfli(vap)) THEN !use explicit qvfl_s
!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        zvari(i,ke1,vap)=qvfl_s(i)/(rhon(i,ke1)*tkvh(i,ke1))
      END DO
      !$ACC END PARALLEL
    END IF
    !Note: "tkvh(ke1) > 0.0" is required!

    !OPTIONAL: diagnostic output of fluxes
    IF ( (PRESENT(tet_flux).AND.PRESENT(vap_flux).AND.PRESENT(liq_flux)) &
        .AND. PRESENT(l_3d_turb_fluxes)                           )  THEN
       IF (l_3d_turb_fluxes) THEN	 
          !$ACC PARALLEL ASYNC(1) IF(lzacc)
!DIR$ IVDEP
          DO k=2, ke1
          !$ACC LOOP GANG VECTOR
             DO i=ivstart, ivend !    gradient     rho     coefficient
  	          tet_flux(i,k) = zvari(i,k,tet)*rhon(i,k)*tkvh(i,k) 
  	          vap_flux(i,k) = zvari(i,k,vap)*rhon(i,k)*tkvh(i,k)  
  	          liq_flux(i,k) = zvari(i,k,liq)*rhon(i,k)*tkvh(i,k)  
             END DO
          END DO
          !$ACC END PARALLEL
       END IF !3d flux output?
    END IF ! all variables present?
    
!------------------------------------------------------------------------------------
! 6)  Berechnung der zu TKE-Quellen gehoerigen Temperaturtendenzen
!     (ausser der Divergenz des Zirkulationsflusses):
!------------------------------------------------------------------------------------

    IF (ltmpcor) THEN
      IF (.NOT.PRESENT(edr)) THEN
!DIR$ IVDEP
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          ediss(i,ke1)=tke(i,ke1,ntur)**3/(d_m*len_scale(i,ke1))
        END DO
        !$ACC END PARALLEL
      END IF

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(thermik, phasdif)
      DO k=2, ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          thermik=tkvh(i,k)*frh(i,k)
          phasdif=tkvh(i,k)*frm(i,k) *(tur_rcpv*zvari(i,k,vap)+tur_rcpl*zvari(i,k,liq))

          tketens(i,k)=len_scale(i,k)/zaux(i,k,2) *((ediss(i,k)+thermik)/cp_d+phasdif)
        END DO

        ! Beachte:
        ! In 'frm' steht an dieser Stelle der Temp.-Gradient!
        ! Wegen der nachfolgenden Interpolation auf Hauptflaechen,
        ! wird 'tketens' mit der turbulenten Laengenskala multipliziert.
      END DO
      !$ACC END PARALLEL

!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        ttens(i,1)=ttens(i,1)+tinc(tem)*tketens(i,2) /(len_scale(i,1)+len_scale(i,2))

      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=2,ke
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          ttens(i,k)=ttens(i,k)+tinc(tem)*(tketens(i,k)+tketens(i,k+1)) &
                            /(len_scale(i,k)+len_scale(i,k+1))
        END DO
      END DO
      !$ACC END PARALLEL
    END IF   ! ltmpcor

!------------------------------------------------------------------------------------
! 7)  Bestimmung des Zirkulationstermes als zusaetzliche TKE-Flussdichte:
!------------------------------------------------------------------------------------
    !Achtung: Zirkulationsterm revidieren:

    IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=2,ke1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR PRIVATE(fakt, com_len, kohae_len)
        DO i=ivstart, ivend

          IF (k.LT.ke1) THEN
            ! Interpolation des Druckes auf die naechst hoehere
            ! Nebenflaeche:
            lay(i)=(prs(i,k)*dp0(i,k-1)+prs(i,k-1)*dp0(i,k)) &
                      /(dp0(i,k)+dp0(i,k-1))
          ELSE
            lay(i)=ps(i)
          END IF

          ! Achtung: Variation durch Wolkenbedeckung
          fakt=z1-z2*ABS(zaux(i,k,4)-z1d2)
          com_len=MAX( l_pat(i), SQRT(fakt*len_scale(i,k)*l_hori(i)) )

          ! Berechnung der lokalen Kohaerenzlaenge fuer die Parameterisierung
          ! des Zirkulationstermes (kohae_len=l_pat*exnr*grad(tet_v)*r/g):

          fakt=frh(i,k)*lay(i)/(rhon(i,k)*grav**2)
          kohae_len=com_len*SIGN(z1,fakt)*MIN( ABS(fakt), z1 )

          ! Belegung von 'frh' mit der TKE-Flussdichte des Zirkulationstermes
          ! (Zirkulationsfluss oder Drucktransport):

          frh(i,k)=tkvh(i,k)*kohae_len*frh(i,k)

          ! Die Divergenz dieser Flussdichte ist gleichzeitig
          ! eine weitere Quelle fuer thermische Energie.

          ! Addition des Teta-Gradienten, welcher zur Teta-Flussdichte
          ! durch den Zirkulationsterm gehoert:

          IF (ltmpcor) THEN
            zvari(i,k,tet)=zvari(i,k,tet) &
                                -frh(i,k)/(tkvh(i,k)*zaux(i,k,1)*zaux(i,k,2)*cp_d)
          END IF

          !Die kinetische Energie der Zirkulation wird der thermischen Energie
          !entzogen!
        END DO
      END DO
      !$ACC END PARALLEL

    END IF     ! lcircterm

!------------------------------------------------------------------------------------
! 8) Berechnung der Diffusionstendenz von q=SQRT(2*TKE) einschliesslich der
!    q-Tendenz durch den Zirkulationsterm:
!------------------------------------------------------------------------------------

    ! Vorbereitung zur Bestimmung der zugehoerigen Incremente von TKE=(q**2)/2:

    cur_prof => hlp
    upd_prof => zaux(:,:,1)
    sav_prof => zaux(:,:,5)

    IF (ldotkedif .OR. lcircdiff) THEN
              ! ldotkedif: partly implicit vertical diffusion for TKE:  c_diff > 0.0
              ! lcircdiff: circulation term, computed together with TKE diffusion: 
              !                     = lcircterm .and. imode_calcirc==2

      expl_mom => zaux(:,:,2)
      impl_mom => zaux(:,:,3)
      invs_mom => zaux(:,:,4)

      ! Diffusions-Koeffizienten auf NF:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      DO k=2, ke1
        !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
        DO i=ivstart, ivend
          IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
            !------------------------------------------------------
            !test: TKE-Diffusion mit Stab.fnkt. fuer Skalare:
            ! sav_prof(i,k)=c_diff*tkvh(i,k)
            sav_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)
            !------------------------------------------------------
          ELSE !Diffusion in terms of q=SQRT(2*TKE)
            sav_prof(i,k)=c_diff*len_scale(i,k)*tke(i,k,ntur)**2
          END IF
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP SEQ
      DO k=3, ke1
!DIR$ IVDEP
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          expl_mom(i,k)=rhoh(i,k-1)*z1d2*(sav_prof(i,k-1)+sav_prof(i,k)) &
                                             /(hhl(i,k-1)-hhl(i,k))
        END DO
        !Beachte: 
        !'expl_mom' bezieht sich auf HF, also die Fluss-Niveaus fuer die TKE (bzw. q-)-Diffusion.
        !Wegen der spaeteren Nutzung der SUBs 'prep_impl_vert_diff' und 'calc_impl_vert_diff'
        ! muss ein Fluss-Niveau (hier HF) ueber dem Variabl.-Niveau (hier NF) mit gleichem Index liegen.

        !Achtung: In der COSMO-Version werden 'rhon', 'len_scale' und 'tke' einzeln auf HF interpoliert,
        !         was geringe Unterschiede verursacht!
      END DO
      !$ACC END PARALLEL

      ! end data section for expl_mom,impl_mom,invs_mom
    END IF  ! ldotkedif .OR. lcircdiff


    IF (ldotkedif .OR. lcircterm) THEN
      ! lcircterm: circulation term has to be computed: =(pat_len > 0.0)

      IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2, ke1
!DIR$ IVDEP
          DO i=ivstart, ivend
            sav_prof(i,k)=z1d2*tke(i,k,ntur)**2 !TKE
          END DO
        END DO
        !$ACC END PARALLEL
      ELSE !Diffusion in terms of q=SQRT(2*TKE)
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2, ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            sav_prof(i,k)=tke(i,k,ntur)         !q=SQRT(2*TKE)
            dicke(i,k)=dicke(i,k)*tke(i,k,ntur) !related effective discretization momentum
          END DO
        END DO
        !$ACC END PARALLEL

        !Das Feld 'dicke' wird bei "k=ke1" nicht benoetigt und war zuvor dort auch nicht belegt!
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
!DIR$ IVDEP
        DO i=ivstart, ivend
          sav_prof(i,ke1)=tke(i,ke1,ntur) !q at the surface layer 
        END DO
        !$ACC END PARALLEL
      END IF

      !Beachte: 
      !Das Feld 'tke' enthaelt nicht TKE sondern q!

    END IF  ! ldotkedif .or. lcircterm


    ! Aufnahme des Zirkulationstermes:
    IF (lcircterm) THEN !Der Zirkulationsterm muss berechnet werden

      !  Interpolation der Zirulationsflussdichte auf Hauptflaechen:

      !Beachte:
      !'frm' ist frei. 
      !'frh' enthaelt bislang eine TKE-Flussdichte (bis auf den Faktor 'rhon'),
      ! deren Vertikalprofil im wesentlichen durch d_z(tet_v)**2 bestimmt ist,
      ! was zumindest in der Prandtl-Schicht prop. zu 1/len_scale ist.
      !Die Interpolation von "rhon*frh" auf Hauptflaechen erfolgt daher mit
      ! 'rhon*frh*len_scale'. Anderenfalls ist mit grossen Interpolationsfehlern zu rechnen.
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          frh(i,k)=rhon(i,k)*frh(i,k)*len_scale(i,k) !skalierte Flussdichte auf NF

            !Achtung: In COSMO-Version ist "imode_tkediff=1".
            !         Es wird in dieser Version aber 'frh' an dieser Stelle durch 'tke' dividiert,
            !         was geringe Unterschiede verursacht:"
            !frh(i,k)=rhon(i,k)*frh(i,k)/tke(i,k,ntur)*len_scale(i,k) !skalierte Flussdichte auf NF

            !if (k.lt.ke1) then
            !   frh(i,k)=rhon(i,k)*grav*tkvh(i,k)*t(i,k)/t_g(i)/tke(i,k,ntur)*len_scale(i,k)
            !else
            !   frh(i,k)=rhon(i,k)*grav*tkvh(i,k)/tke(i,k,ntur)*len_scale(i,k)
            !endif
        END DO
      END DO
      !$ACC END PARALLEL

      ! Korrektur der TKE-Profile durch die Zirkulations-Tendenz:
      IF (imode_calcirc.EQ.1) THEN !expliziten Berechnung der Zirkulationstendenz

        cur_prof => sav_prof
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP SEQ
        DO k=3,ke1
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR
          DO i=ivstart, ivend
            frm(i,k)=(frh(i,k)+frh(i,k-1))/(len_scale(i,k)+len_scale(i,k-1)) !interpolierte Flussdichte auf HF
          END DO
        END DO
        !$ACC END PARALLEL

        k=2
!DIR$ IVDEP
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          upd_prof(i,k)=sav_prof(i,k)+frm(i,k+1)/dicke(i,k)

            !Achtung: In COSMO-Version ist "imode_tkediff=1". 
            !Da 'frh' bereits durch 'tke' dividiert wurde, muss fuer die exakte COSMO-Version der 'tke'-Faktor in 'dicke'
            !wieder beseitigt werden:
            !upd_prof(i,k)=sav_prof(i,k)+frm(i,k+1)*tke(i,k,ntur)/dicke(i,k)

        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP SEQ
        DO k=3,ke
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR
          DO i=ivstart, ivend
            upd_prof(i,k)=sav_prof(i,k)-(frm(i,k)-frm(i,k+1))/dicke(i,k)

              !upd_prof(i,k)=sav_prof(i,k)-(frm(i,k)-frm(i,k+1))*tke(i,k,ntur)/dicke(i,k)

              !Achtung: In der COSMO-Version wird die explizite q-Tendenz durch den Zirkulationsterm
              !erst nach der TKE-Diffusion hinzuaddiert (verursacht geringe Unterschiede):
              !upd_prof(i,k)=sav_prof(i,k)
              !tketens(i,k)=-(frm(i,k)-frm(i,k+1))*tke(i,k,ntur)/dicke(i,k)
          END DO
        END DO
        !$ACC END PARALLEL

        !Beachte:
        !'dicke' enthaelt bereits den Faktor "1/dt_tke" (sowie den Faktor 'tke' bei "imode_tkediff=1").
        !'upd_prof' enthaelt das um die Zirkulations-Tendenz aufdatierte TKE-Profil
        ! (oder q-Profile bei "imode_tkediff=1")

        IF (PRESENT(r_air)) THEN

          !Zuschlag durch Volumenterm aus der Divergenzbildung:
#ifdef __INTEL_COMPILER
          FORALL(k=kcm:ke, i=ivstart:ivend) &
            upd_prof(i,k)=upd_prof(i,k)+frh(i,k)*z1d2*(r_air(i,k-1)-r_air(i,k+1)) &
                                                        /(len_scale(i,k)*dicke(i,k))
#else
          !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
          !$ACC LOOP SEQ
          DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
            !$ACC LOOP GANG VECTOR
            DO i=ivstart, ivend
!Achtung: Korrektur
            ! upd_prof(i,k)=upd_prof(i,k)+dt_tke*frh(i,k)*z1d2*(r_air(i,k-1)-r_air(i,k+1)) &
              upd_prof(i,k)=upd_prof(i,k)+frh(i,k)*z1d2*(r_air(i,k-1)-r_air(i,k+1)) &
                                                             /(len_scale(i,k)*dicke(i,k))
              !'frh' enthaelt die Zirkultions-Flussdichte der TKE auf NF (skaliert mit 'len_scale').
            END DO
          END DO
          !$ACC END PARALLEL
#endif
        ENDIF ! PRESENT(r_air)

        !Bereucksichtige Zirkulations-Tendenz:
        itndcon=1 !indem 'upd_prof' auf rechter Seite der impliz. Diff.-Gl. benutzt wird.

            !Achtung: Um die COSMO-Version exakt nachzubilden, darf die Zirkulationstendenz nicht bei der
            !impliziten Diffusions-Gleichung benutzt werden:
            !itndcon=0
            !Fuer die expliziten Diff.-Anteile wird aber 'cur_prof' benutzt.

      ELSE ! now: imode_calcirc /= 1
           ! quasi implizite Berechnung der Zirkulatinstendenz (entspricht "lcircdiff=T")

!DIR$ IVDEP
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
          cur_prof(i,2)=sav_prof(i,2)
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP SEQ
        DO k=3,ke1
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR PRIVATE(wert)
          DO i=ivstart, ivend
            wert=(frh(i,k)+frh(i,k-1))/((len_scale(i,k)+len_scale(i,k-1))*expl_mom(i,k))
            cur_prof(i,k)=(cur_prof(i,k-1)-sav_prof(i,k-1)+wert)+sav_prof(i,k)
          END DO
        END DO
        !$ACC END PARALLEL

        !Beachte:
        !'cur_prof' enthaelt ein virtuelles TKE-Profil (oder q-Profile bei "imode_tkediff=1"), 
        ! dessen Diffusions-Tendenz die Zirkulations-Tendenz einschliesst.

        !Bereucksichtige Zirkulations-Tendenz:
        itndcon=0 !indem 'cur_prof' auf rechter Seite der impliz. Diff.-Gl. benutzt wird.

        !Fuer die expliziten Diff.-Anteile wird ebenfalls 'cur_prof' benutzt.

      END IF    ! imode_calcirc

    ELSEIF (ldotkedif) THEN
 
      cur_prof => sav_prof

      itndcon=0 !'cur_prof' wird auf rechter Seite der impliz. Diff.-Gl.
                ! und fuer explizite Diff.-Anteile benutzt.
    END IF

    ! Aufdatieren des TKE-Profils durch die (erweiterte) Diffusions-Tendenz 

    !$ACC WAIT
    IF (ldotkedif .OR. lcircdiff) THEN

      !'frm', 'frh' und 'len_scale' sind frei.
      !In den Diffusionroutinen wird vorausgesetzt, dass ein Flussniveau mit gleichem
      !Vertikalindex wie ein Konzentrationsniveau gerade ueber letzterem liegt.
      !Die bisherige Hauptflaechenindizierung musste daher fuer die Uebergabefelder
      !der Routinen 'prep_impl_vert_diff' und 'calc_impl_vert_diff' um "1" verschoben werden.

      CALL prep_impl_vert_diff( lsflucond=.FALSE., ldynimpwt=ldynimp, lprecondi=lprecnd,   &
                                i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1,                &
                                disc_mom=dicke,    expl_mom=expl_mom,                      &
                                impl_mom=impl_mom, invs_mom=invs_mom,                      &
                                invs_fac=frh, scal_fac=frm, impl_weight=impl_weight,       &
                                lacc=lzacc )

           !Achtung: q_Diff:    disc_mom=sav_prof*dicke !!!

      ! Berechnung der vertikalen Diffusionstendenzen von TKE=z1d2*q**2:

      eff_flux => len_scale
         
      CALL calc_impl_vert_diff ( lsflucond=.FALSE.,lprecondi=lprecnd,                      &
                                 leff_flux=(kcm.LE.ke), itndcon=-itndcon,                  &
                                 i_st=ivstart, i_en=ivend,k_tp=1, k_sf=ke1,                &
                                 disc_mom=dicke,    expl_mom=expl_mom,                     &
                                 impl_mom=impl_mom, invs_mom=invs_mom,                     &
                                 invs_fac=frh, scal_fac=frm, cur_prof=cur_prof,            &
                                 upd_prof=upd_prof, eff_flux=eff_flux, lacc=lzacc )

           !Achtung: q_Diff:    disc_mom=sav_prof*dicke !!!

      !Beachte:
      !Bei "imode_diff=1" erfolgt q-Diffusion, so dass 'disc_mom' und 'expl_mom' den zusaetzlichen 
      ! Faktor 'sav_prof'='tke'=q enthalten!
      !'upd_prof' enthaelt jetzt die mit der Diffusionstendenz aufdatierten (modifizierte) TKE-Werte.
      !Weil "itndcon<=0", bleiben auf 'cur_prof' die Eingangsprofile erhalten.
      !'eff_flux' enthaelt die effektiven Flussdichten (positiv abwaerts) der (semi-)impliziten
      ! Vertikaldiffusion.

      IF (lcircdiff) THEN !es wurden virtuelle Effektiv-Profile benutzt
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            upd_prof(i,k)=sav_prof(i,k)+upd_prof(i,k)-cur_prof(i,k) !aufdatierte echte Profile
          END DO
        END DO
        !$ACC END PARALLEL
      END IF

      IF (PRESENT(r_air)) THEN
        ! Zuschlag durch Volumenterm innerhalb der Rauhigkeitsschicht:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP SEQ
        DO k=ke,kcm,-1 !innerhalb der Rauhigkeitsschicht
!DIR$ IVDEP
          !$ACC LOOP GANG VECTOR PRIVATE(wert)
          DO i=ivstart, ivend
            wert=(eff_flux(i,k)*dp0(i,k)+eff_flux(i,k+1)*dp0(i,k-1))/(dp0(i,k)+dp0(i,k-1))
              !effektive TKE-Flussdichte interpoliert auf die k-te Nebenflaeche,
              !wobei sich 'eff_flux(:,k)' auf die (k-1)-te Hauptflaeche bezieht!
          !Achtung: Korrektur
          ! upd_prof(i,k)=upd_prof(i,k)+dt_tke*wert*z1d2*(r_air(i,k-1)-r_air(i,k+1))/dicke(i,k)
            upd_prof(i,k)=upd_prof(i,k)+wert*z1d2*(r_air(i,k-1)-r_air(i,k+1))/dicke(i,k)
          END DO
        END DO
        !$ACC END PARALLEL
        !'upd_prof' enthaelt das mit dem Volunmenterm beaufschlagte aufdatierte Profil.
      ENDIF

    END IF   ! IF (ldotkedif .OR. lcircdiff)

!------------------------------------------------------------------------------------
! 9)  Speichern der zugehoerigen q-Tendenzen:
!------------------------------------------------------------------------------------

    IF (ldotkedif .OR. lcircterm) THEN   

      IF (imode_tkediff.EQ.2) THEN !Diffusion in terms of TKE
        !'upd_prof' ist ein TKE-Profil:
        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,ke 
!DIR$ IVDEP
          DO i=ivstart, ivend
            !----------------------------------------------------------------------------
            !test:
          ! tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - sav_prof(i,k) )*fr_tke/tke(i,k,ntur)
            tketens(i,k)=( SQRT( 2*MAX( upd_prof(i,k), z0 ) ) - tke(i,k,ntur) )*fr_tke
            !----------------------------------------------------------------------------
          END DO
        END DO
        !$ACC END PARALLEL

      ELSE !Diffusion in terms of q=SQRT(2*TKE)
        !'upd_prof' ist ein q-Profil:

        !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k=2,ke
!DIR$ IVDEP
          DO i=ivstart, ivend
            !Achtung:
            tketens(i,k)=( MAX( upd_prof(i,k), z0 ) - tke(i,k,ntur) )*fr_tke

            !Achtung: Bei der COSMO-Version wird erst hier die explizite Zirkulations-Tendenz hinzuaddiert:
            !tketens(i,k)=MAX( -tke(i,k,ntur), tketens(i,k)+( upd_prof(i,k) - tke(i,k,ntur) ) )*fr_tke
          END DO
        END DO
        !$ACC END PARALLEL
      END IF
      !'tketens' enthaelt jetzt immer eine q-Tendenz!

      !Am Unterrand gibt es keine q-Tendenz durch Diffusionsterme:
!DIR$ IVDEP
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        tketens(i,ke1)=z0
      END DO
      !$ACC END PARALLEL

      ! Optionale vertikale Glaettung der erweiterten Diffusionstendenz von q=SQRT(2*TKE):
      IF (tndsmot.GT.z0) THEN
        !$ACC WAIT
        CALL vert_smooth ( i_st=ivstart, i_en=ivend, k_tp=1, k_sf=ke1, &
                           disc_mom=dicke, cur_tend=tketens, vertsmot=tndsmot, lacc=lzacc )
      END IF

    ELSE !keine q-Tendenzen, weder durch TKE-Diffusion noch durch den Zirkulationsterm

      ! Zuruecksetzen der q-Tendenzen:
      !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k=2,ke1
!DIR$ IVDEP
        DO i=ivstart, ivend
          tketens(i,k)=z0
        END DO
      END DO
      !$ACC END PARALLEL

    END IF   ! ldotkedif .OR. lcircterm


!------------------------------------------------------------------------------------
! 10) Interpolationen auf Hauptflaechen fuer die Standardabweichnung
!     des Saettigungsdefizites:
!------------------------------------------------------------------------------------

!DIR$ IVDEP
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO i=ivstart, ivend
      rcld(i,1)=rcld(i,2)
    END DO
    !$ACC END PARALLEL

#ifdef __INTEL_COMPILER
    FORALL(k=2:kem-1, i=ivstart:ivend) &
        rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
#else
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) IF(lzacc)
    !$ACC LOOP SEQ
    DO k=2,kem-1
!DIR$ IVDEP
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
        rcld(i,k)=(rcld(i,k)+rcld(i,k+1))*z1d2
      END DO
    END DO
    !$ACC END PARALLEL
#endif

    ! Fuer die unterste Hauptflaeche (k=ke) wird bei kem=ke
    ! der Wert auf der entspr. Nebenflaeche beibehalten.

    ! Just do some check printouts:
    IF (ldebug) THEN
      DO i = ivstart, ivend
        IF (i ==  1 .AND. iblock == 1 .AND. my_cart_id == 0) THEN
          WRITE(*,'(A       )') '  '
          WRITE(*,'(A,2I5   )') 'TURB-DIAGNOSIS diffusion: iblock = ', iblock, i

          WRITE(*,'(A       )') ' Some output parameters:'
          WRITE(*,'(A,F28.16)') '   gz0              :  ', gz0         (i)
          WRITE(*,'(A,F28.16)') '   tvm              :  ', tvm         (i)
          WRITE(*,'(A,F28.16)') '   tvh              :  ', tvh         (i)
          WRITE(*,'(A,F28.16)') '   tfm              :  ', tfm         (i)
          WRITE(*,'(A,F28.16)') '   tfh              :  ', tfh         (i)
          WRITE(*,'(A,F28.16)') '   tke   ke         :  ', tke         (i,ke,ntur)
          WRITE(*,'(A,F28.16)') '   tke   ke1        :  ', tke         (i,ke1,ntur)
          do k = 1, ke
            WRITE(*,'(A,I4,A,F28.16)') '   tkvh  ',k,'       :  ', tkvh        (i,k)
          enddo
          do k = 1, ke
            WRITE(*,'(A,I4,A,F28.16)') '   tkvm  ',k,'       :  ', tkvm        (i,k)
          enddo
          WRITE(*,'(A,F28.16)') '   rcld  ke         :  ', rcld        (i,ke)
          WRITE(*,'(A,F28.16)') '   rcld  ke1        :  ', rcld        (i,ke1)
          WRITE(*,'(A,F28.16)') '   utens ke         :  ', utens       (i,ke)
          WRITE(*,'(A,F28.16)') '   utens ke1        :  ', utens       (i,ke1)
          WRITE(*,'(A,F28.16)') '   vtens ke         :  ', vtens       (i,ke)
          WRITE(*,'(A,F28.16)') '   vtens ke1        :  ', vtens       (i,ke1)
          WRITE(*,'(A,F28.16)') '   ttens ke         :  ', ttens       (i,ke)
          WRITE(*,'(A,F28.16)') '   ttens ke1        :  ', ttens       (i,ke1)
          WRITE(*,'(A,F28.16)') '  qvtens ke         :  ', qvtens      (i,ke)
          WRITE(*,'(A,F28.16)') '  qvtens ke1        :  ', qvtens      (i,ke1)
          WRITE(*,'(A,F28.16)') '  qctens ke         :  ', qctens      (i,ke)
          WRITE(*,'(A,F28.16)') '  qctens ke1        :  ', qctens      (i,ke1)
          WRITE(*,'(A,F28.16)') '   zvari ke   1     :  ', zvari       (i,ke ,1)
          WRITE(*,'(A,F28.16)') '   zvari ke1  1     :  ', zvari       (i,ke1,1)
          WRITE(*,'(A,F28.16)') '   zvari ke   2     :  ', zvari       (i,ke ,2)
          WRITE(*,'(A,F28.16)') '   zvari ke1  2     :  ', zvari       (i,ke1,2)
          WRITE(*,'(A,F28.16)') '   zvari ke   3     :  ', zvari       (i,ke ,3)
          WRITE(*,'(A,F28.16)') '   zvari ke1  3     :  ', zvari       (i,ke1,3)
          WRITE(*,'(A,F28.16)') '   zvari ke   4     :  ', zvari       (i,ke ,4)
          WRITE(*,'(A,F28.16)') '   zvari ke1  4     :  ', zvari       (i,ke1,4)
          WRITE(*,'(A,F28.16)') '   zvari ke   5     :  ', zvari       (i,ke ,5)
          WRITE(*,'(A,F28.16)') '   zvari ke1  5     :  ', zvari       (i,ke1,5)
          WRITE(*,'(A       )') '  '
        ENDIF
      ENDDO
    ENDIF

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  END IF   ! IF iini == 1

  !$ACC END DATA
  !$ACC WAIT

END SUBROUTINE turbdiff

!==============================================================================

END MODULE turb_diffusion
