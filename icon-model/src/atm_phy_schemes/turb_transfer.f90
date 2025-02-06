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

! Source module for computing the coefficients for turbulent transfer
!
! Description of *turb_transfer*:
!   This  module calculates the coefficients for turbulent transfer.
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
!   turbtran
!
! called from the turbulence interface routine of the model.
!
!-------------------------------------------------------------------------------

MODULE turb_transfer

!-------------------------------------------------------------------------------
!
! Documentation History:
!
!  The history of all these modifications is as follows, where those belonging to the formal
!   reorganization of the whole package (atmospheric turbulence and surface-to-atmosphere transfer)
!   are now in the header of MODULE 'turb_utilities', containing various common SUBs for 'turbdiff'
!   and 'turtran' (related moist thermodynamicds and the treatment of turbulent budget equations)
!   and also the blocked code for semi-implicit vertical diffusion. The new blocked version of SUB 'turbtran'
!   is now in MODULE 'turb_transfer':
!
!              2010/09/30 Matthias Raschendorfer
!  Substitution of 'itype_diag_t2m' by 'itype_synd' being already used for that purpose
!   "1": alternative SYNOP-digansostics related to previous Lewis scheme (not included here)
!   "2": SYNOP-diagnostics according to current transfer scheme using SYNOP-z0 'z0d'
!   "3": like "2" but using 'z0d' only for 10m-wind and a specific roughness layer profile for
!        2m-temperature and -humidity.
!  Including the adiabatic lapse rate correction for 't_2m' for "itype_synd.EQ.3" as well.
!              2011/03/23 Matthias Raschendorfer
!  Correction of two bugs introduced together with the last modifications in SUB 'turbtran'
!   (related to SUB 'diag_level' and the 'tet_l'-gradient used for the flux output).
!  Substitution of run time allocations because of bad performance on some computers.
!              2014/07/28 Matthias Raschendorfer
!  Removing a bug in formular for 'rcld(:,ke1)' in SUB 'turbtran'
!   -> influence on near-surface temperature and - humidity
!              2015/08/25 Matthias Raschendorfer
! Adopting other development within the ICON-version (namely by Guenther Zaengl) as switchable options
!  related to the following new selectors and switches:
!   imode_rat_sea, imode_vel_min, imode_charpar and lfreeslip.
! Rearranging the development by Matthias Raschendorfer that had not yet been transferred to COSMO as switchable
!  options related to the following switches:
!  and selectors:
!   itype_diag_t2m, imode_syndiag, imode_trancnf, imode_tkemini, imode_lamdiff
!  and a partly new (more consistent) interpretation of:
!   imode_tran, icldm_tran, itype_sher, and itype_diag_t2m
! Controlling numerical restrictions gradually namely by the parameter:
!  ditsmot
! Using the arrays 'tvm', 'tvh' and 'tkm', allowing an easier formulation of transfer-resistances.
!              2016-05-10 Ulrich Schaettler 
! Splitting this module from the original module 'organize_turbdiff' as it was used by ICON before.
! Moving declarations, allocation and deallocations of ausxilary arrays into MODULE 'turb_data'.
!
!-------------------------------------------------------------------------------

! Modules used:

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif
USE mo_exception,         ONLY: message_text, message
!-------------------------------------------------------------------------------
! Parameter for precision
!-------------------------------------------------------------------------------

USE mo_kind,         ONLY :   &
    wp              ! KIND-type parameter for real variables

!-------------------------------------------------------------------------------
! Mathematical and physical constants
!-------------------------------------------------------------------------------

USE mo_mpi,                ONLY : get_my_global_mpi_id

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
    t0_melt  => tmelt,    & ! absolute zero for temperature (K)
    b3       => tmelt,    & !          -- " --

    rdv, con_m, con_h, grav      

USE mo_lookup_tables_constants, ONLY : &
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    b1       => c1es,     & ! variables for computing the saturation steam pressure
    b2w      => c3les,    & ! over water (w) and ice (e)
    b4w      => c4les       !               -- " --

!-------------------------------------------------------------------------------
! Turbulence data (should be the same in ICON and COSMO)
!-------------------------------------------------------------------------------

USE turb_data, ONLY : &

! Numerical constants and parameters:
! -----------------------------------

    tkhmin,       & ! minimal diffusion coefficients for heat
    tkmmin,       & ! minimal diffusion coefficients for momentum
    ditsmot,      & ! smoothing factor for direct time-step iteration
    epsi,         & ! relative limit of accuracy for comparison of numbers
    it_end,       & ! number of initialization iterations (>=0)

! Parameters describing physical properties of the lower boundary 
! of the atmosphere:
!---------------------------------------------------------------
!
    rlam_mom,     & ! scaling factor of the laminar boudary layer for momentum
    rlam_heat,    & ! scaling factor of the laminar boudary layer for heat

    rat_lam,      & ! vapour/heat ratio of laminar scaling factors (over land)
    rat_sea,      & ! sea/land ratio of laminar scaling factors for heat (and vapor)
    rat_glac,     & ! glacier/land ratio of laminar scaling factors for heat (and vapor)

    rat_can,      & ! factor for the canopy height

    z0m_dia,      & ! roughness length of a typical synoptic station [m]

    alpha0,       & ! lower bound for Charnock-parameter
    alpha1,       & ! parameter scaling the molek. roughness of water waves

    z0_ice,       & ! roughness length of sea ice

    len_min,      & ! minimal turbulent length scale [m]
    tur_len,      & ! maximal turbulent length scale [m]
    vel_min,      & ! minimal velocity scale [m/s]
    prfsecu,      & ! relat. secur. fact. for prof. funct.  (out of ]0; 1[)

    akt,          & ! von Karman-constant
    d_h=>d_heat,  & ! factor for turbulent heat dissipation
    d_m=>d_mom,   & ! factor for turbulent momentum dissipation

    rsur_sher,    & ! scaling factor for additional shear-forcing by Non-Turbulent subgrid Circulations (NTCs)
                    ! or via Lower Limits of Diffusion-Coefficients (LLDCs) in the surface layer

    ! derived parameters calculated in 'turb_setup'
    tet_g, rim, b_m, b_h, sm_0, sh_0,   &
    a_3, a_5 ,a_6,                      &
    tur_rcpv, tur_rcpl,                 &

    ! used derived types
    modvar,       & !

! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

    lprfcor,      & ! using the profile values of the lowest main level instead of
    ltst2ml, &      ! test required, whether  2m-level is  above the lowest main-level
                    ! F:  2m-level is assumed to be always below the lowest main-level
    ltst10ml, &     ! test required, whether 10m-level is  above the lowest half-level   
                    ! F: 10m-level is assumed to be always below the lowest half-level
    lfreeslip       ! free-slip lower boundary condition (enforeced zero-flux condition for
                    ! for all diffused variables, only for idealized test cases)

USE turb_data, ONLY : &

! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------
!
    imode_tran,   & ! mode of TKE-equation in transfer scheme                    (compare 'imode_turb')
                    !  0: diagnostic equation
                    !  1: prognostic equation (default)
                    !  2: prognostic equation (implicitly positive definit)
    icldm_tran,   & ! mode of water cloud representation in transfer parametr.   (compare 'icldm_turb)
                    ! -1: ignoring cloud water completely (pure dry scheme)
                    !  0: no clouds considered (all cloud water is evaporated)
                    !  1: only grid scale condensation possible
                    !  2: also sub grid (turbulent) condensation considered
    itype_sher,   & ! type of shear production for TKE
                    ! 0: only vertical shear of horizontal wind
                    ! 1: previous plus horizontal shear correction
                    ! 2: previous plus shear from vertical velocity
    ilow_def_cond,& !type of the default condition at the lower boundary
                    ! 1: zero surface gradient
                    ! 2: zero surface value
    imode_rat_sea,& ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                    ! 1: constant ratio compared to land surface
                    ! 2: with a correction for a strongly overheated SST
    imode_vel_min,& ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                    ! 1: with the constant value "tkesecu*vel_min"
                    ! 2: with a stability dependent correction
    imode_charpar,& ! mode of estimating the Charnock-Parameter
                    ! 1: use a constant value 
                    ! 2: use a wind-dependent value with a constant upper bound
    imode_2m_diag,& ! mode of 2m-diagnostics of temperature and dew-point (related to 'itype_2m_diag')
                    ! (-)1: direct interpolation of temperature and specific humidity
                    ! (-)2: interpol. of conserved quantities and subsequent statistical saturation adjustm.,
                    !       allowing particularly for the diagnostic of cloud water at the 2m-level (fog)
                    !  > 0: extra pressure-calculat. at 2m-level
                    !  < 0: surface pressure applied at 2m-level
   imode_nsf_wind,& ! mode of local wind-definition at near-surface levels (applied for 10m wind
                    !  diagnostics as well as for calculation of sea-surface roughness)
                    ! 1: ordinary wind speed (magnitude of grid-scale averaged wind-vector)
                    ! 2: including relative wind-speed amplification by NTCs  (at "rsur_shear>0")
                    !                                        or even by LLDCs (at "imode_suradap=3")
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

    imode_tkemini,& ! mode of adapting q=2TKE**2 and the TMod. to Lower Limits for Diff. Coeffs. (LLDCs)
                    ! 1: LLDC treated as corrections of stability length without any further adaptation
                    ! 2: TKE adapted to that part of LLDC representing so far missing shear forcing, while the
                    !     assumed part of LLDC representing missing drag-forces has no feedback to the TMod.
    imode_suradap,& ! mode of adapting surface-layer profile-functions to Lower Limits of Diff. Coeffs. (LLDCs)
                    ! 0: no adaptations at all
                    ! 1: removing the artific. drag contrib. by the LLDC for momentum at level "k=ke"
                    ! 2: "1" and also removing  shear contrib. by LLDCs at level "k=ke" 
                    ! 3: "1" and employing shear contrib. by LLDCs at surface-level "k=ke1"
                    !Notice:
                    !Any shear contrib. by NTCs or LLDCs is only considered at surf.-lev, if "rsur_sher.GT.0".
    imode_lamdiff   ! mode of considering the laminar limit for diff.coeffs. at surface layer
                    ! 1: not applied for the profile functions in case of "imode_trancnf.GE.2"
                    ! 2: always applied

USE turb_data, ONLY : &

    ! numbers and indices

    ninv    ,     & ! number of invariant scalar variables being conserved during 'vap'<->'liq'-transistions ('tet_l', 'h2o_g')
    nvel    ,     & ! number of velocity-components active for turbulece ('u_m', 'v_m')
    nmvar   ,     & ! number of progn. variables beding active for turbulence
    ndim    ,     & ! (positive) limit of last dimension used for 'zaux' and 'zvari'
    mom, sca,     & ! indices for momentum- and scalar-variables
    u_m     ,     & ! zonal velocity-component at the mass center
    v_m     ,     & ! meridional ,,      ,,    ,, ,,   ,,    ,,
    tet_l   ,     & ! liquid-water potential temperature
    h2o_g   ,     & ! total water ('qv+qc')
    tet     ,     & ! potential temperature
    vap     ,     & ! water vapor (specific humidity 'qv')
    liq             ! liquid water (mass fraction 'qc') 

    !Note: It always holds: tem=tet=tet_l=tem_l and vap=h2o_g (respective usage of equal indices)!


!-------------------------------------------------------------------------------
! Control parameters for the run
!-------------------------------------------------------------------------------

! ICON data have to be declared for these variables, which is done later on
!-------------------------------------------------------------------------------

! Switches controlling other physical parameterizations:
USE mo_lnd_nwp_config,       ONLY: lseaice, llake, lterra_urb, itype_kbmo
USE mo_atm_phy_nwp_config,   ONLY: lcuda_graph_turb_tran
!   
USE turb_data,         ONLY:   &
    itype_2m_diag      ! type of 2m-diagnostics for temperature and -dewpoint
                       ! 1: Considering a fictive surface roughness of a SYNOP lawn
                       ! 2: Considering the mean surface roughness of a grid box
                       !    and using an exponential roughness layer profile
!
USE turb_utilities,          ONLY:   &
    turb_setup,                      &
    adjust_satur_equil,              &
    solve_turb_budgets,              &
    zexner, zpsat_w,                 &
    alpha0_char

!-------------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
    lsclm, lsurflu, i_cal, i_upd, i_mod, imb, &
    SHF, LHF
#endif
!SCLM---------------------------------------------------------------------------

USE mo_fortran_tools, ONLY: set_acc_host_or_device
!===============================================================================

IMPLICIT NONE

PUBLIC  :: turbtran

!===============================================================================

!-------------------------------------------------------------------------------

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

SUBROUTINE turbtran (                                                         &
!
          iini, ltkeinp, lgz0inp, lstfnct, lsrflux, lnsfdia, lrunscm,         &
          ladsshr,                                                            &
!
          dt_tke, nprv, ntur, ntim,                                           &
! 
          nvec, ke, ke1, kcm, iblock, ivstart, ivend,                         &
!
          l_pat, l_hori, hhl, fr_land, l_lake, l_sice, gz0,                   &
          rlamh_fac, sai, urb_isa,                                            &  
          t_g, qv_s, ps, u, v, t, qv, qc, epr,                                &
!
          tcm, tch, tvm, tvh, tfm, tfh, tfv, tkr,                             &

          tke, tkvm, tkvh, rcld,                                              &
          hdef2, dwdx, dwdy,           & ! optional for itype_sher=2
!
          edr, tketens,                                                       &
!
          t_2m, qv_2m, td_2m, rh_2m, u_10m, v_10m,                            &
          shfl_s, qvfl_s, umfl_s, vmfl_s,                                     &
!
          lacc, opt_acc_async_queue)
!-------------------------------------------------------------------------------
!
! Note:
! It is also possible to use only one time level for TKE using "ntim=1" and thus "nprv=1=ntur".
!
! Description:
!
!     Es werden die Transfer-Geschwindigkeiten fuer den Austausch von Impuls,
!     sowie fuehlbarer und latenter Waerme bestimmt und die Modellwerte
!     fuer die bodennahen Messwerte (in 2m und 10m) berechnet.
!
! Method:
!
!     Hierzu wird der gesamte Bereich von den festen Oberflachen am
!     Unterrand des Modells bis hin zur untersten Hauptflaeche in
!     die drei Teilbereiche:
!
!     - laminare Grenzschicht (L-Schicht)
!     - turbulente Bestandesschicht (B-Schicht)
!     - turbulente Prandtl-Schicht (P-Schicht)
!
!     aufgeteilt. Fuer jeden dieser Teilbereiche wird (getrennt nach
!     skalaren Eigenschaften und Impuls) ein zugehoeriger Transport-
!     widerstand berechnet, der gleich einer effektiven Widerstands-
!     laenge ( dz_(sg, g0, 0a)_(h,m) ) dividiert durch den Diffusions-
!     koeffizienten am Unterrand der P-Schicht (Niveau '0') ist.
!     Die Konzentrationen am Unterrand der B-Schicht, also im
!     Abstand der L-Schicht-Dicke entlang der festen Oberflaechen,
!     haben den Index 'g' (ground) und die Oberflaechenkonzentrationen
!     den Index 's'. Groessen fuer den Imoulst haben den Index 'm' (momentum)
!     und solche fuer skalare Eigenschaften 'h' (heat).
!     Der Widerstand der P-Schicht vom Niveau '0' bis zum Niveau 'a'
!     (atmospheric) der untersten Hauptflaeche wird durch vertikale
!     Integration der Modellgleichungen in P-Schicht-Approximation
!     (vertikal konstante Flussdichten, turbulente Laengenskala lin. Funkt.
!      von der Hoehe) gewonnen.
!     Dabei wird das atmosphaerische Turbulenzschema aus der Subroutine
!     'turbdiff' benutzt, so dass alo keine empirischen Profilfunktionen
!     benutzt werden. Zur Vereinfachung der Integration  wird das Produkt
!     aus turbulenter Geschwindigkeitsskala 'q' und der Stabilitaetsfunktion
!     s(h,m), alo die stabilitaetsabhaengige turb. Geschwindigkeitsskala
!     'v' innerhalb der P-Schicht als linear angesehen.
!     Die turb. Laengenskala im Niveau '0' wird mit der Rauhigkeitslaenge 'z0'
!     (multipliziert mit der v.Kaman-Konstanten) gleichgesetzt. Formal werden
!     dann fuer das Nieveau '0' Vertikalgradienten und auch Diffusions-
!     koeffizienten abgeleitet.
!     Unter der Annahme, dass 'v' innerhalb der B-Schicht konstant bleibt,
!     ergibt sich die laminare Widerstandslaenge dz_sg als prop. zu 'z0'
!     und die Widerstandsstrecke durch die B-Schicht als prop. zu
!     'z0*ln(delta/z0)', wobei 'delta' die Dicke der L-Schicht ist, die der
!     Abstand von einer ebenen Wand sein soll in dem der turbulente
!     Diffusionskoeffizient fuer Impuls gleich dem molekularen ist.
!     Ferner wird angenommen, dass die Widerstaende durch die L- und
!     B-Schicht prop. zur effektiven Quellflaech der Bestandeselemente
!     zuzueglich der Grundflaeche des Erdbodens sind. Die Bestandesoberflaechen
!     werden durch den Wert 'sai' (surface area index) ausgedrueckt und setzt
!     sich aus dem Flaechenindex der transpirierenden Oberflaechen 'lai'
!     (leaf area index) und dem fuer die nicht transpirierenden Flaechen
!     zusammen. Im Falle nicht benetzter Oberlfaechen hat die latente Waerme
!     i.a. eine kleinere Quellflaeche als die fuehlbare Waerme, so dass die
!     Wiederstaende fuer beide Groessen unterschieden werden muessten.
!     Um dies zu vermeiden, wird nur der Widerstand fuer die fuehlbare Waerme
!     berechnet. Dafuer wird aber bei der Berechnung der effektiven
!     Oberflaechenkonzentration 'qv_s' der spez. Feuchtigkeit in Subroutine
!     'terra1' dieser Effekt beruecksichtigt.
!     Beim vertikalen Impulstransport ist aber noch die zusaetzliche
!     Impulssenke innerhalb der B-Schicht durch die Wirkung der Formreibungs-
!     kraft zu beruecksichtigen, was durch einen zusaetzlichen Flaechenindex
!     'dai' (drag area index) bewerkstelligt wird.
!
!     Die Vertikalprofile aller Eigenschaften innerhalb der P-Schicht ergeben
!     sich aus dem vertikal integrierten Turbulenzmodell in P-Schicht-
!     Approximation zu logarithmischen Funktionen, welche durch die
!     thermische Schichtung modifiziert sind. Wie bereits erwaehnt, ist die
!     Stabilitaetsfunktion nur noch von Konstanten des atmosphaerischen
!     Turbulenzmodells abhaengig. Das Transferschema ist somit auch automatisch
!     konsistent zum oben anschliessenden Turbulenzmodell formuliert.
!
!     Die Profilfunktionen innerhalb der B-Schicht ergeben sich aus der
!     Annahme eines Gleichgewichtes zwischen vertikalen Flussdichtedivergenzen
!     und Quellstaerken durch die laminaren Grenzschichten der Rauhigkeits-
!     elemente bei vertikal konstanten Bestandeseigenschaften zu exponentiellen
!     Funktionen. Durch die Bedingung eines glatten Ueberganges zwischen beiden
!     Profiltypen im Niveau '0' und der Bedingung, dass im Abstand einer
!     effektiven Bestandesdicke 'Hb' unterhalb des Nieveaus '0' die Bestandes-
!     profile in die Konzentration am Unterrand der B-Schicht (Niveau mit
!     Index 'g') uebergehen, ist das gesamte Transferschema geschlossen und
!     es kann auch der "drag area index" 'dai', sowie die Bestandeshoehe
!     'Hb' selbst eliminiert werden.
!
!     Zur Charakterisierung des Oberflaechentransfers werden dann nur die
!     externen Parameter 'z0', 'sai', 'lai' und je ein globaler Parameter
!     fuer den laminaren Grenzschichtwiderstand des skalaren - und des
!     Impulstransportes benoetigt '(lam_(h,m)'. Hieraus koennte auch eine
!     aequivalente Rauhigkeitslaenge fuer Skalare 'z0h' berechnet werden.
!     Die Oberfalaechenkonzentrationen (Niveau mit Index 's') fuer die skalaren
!     Groessen werden im Modul 'terra' berechnet. Fuer den Impuls gilt die
!     Haftbedingung. Im Grundniveau des atmosphaerischen Modells
!     (Niveau '0') verschwindet also der Wind i.a. nicht; dies ist erst
!     entlang der festen Oberfalechen der Fall. Die bodennahen synoptischen
!     Niveaus werden nun vom Niveau 'z=-Hb', also von der effektiven Bestandes-
!     grundflaeche (Umsatzniveau) aus gezaehlt. Ist z.B. 'Hb>2m', werden
!     die 2m-Werte entlang der exponentiellen Bestandesprofile ausgeweret.
!     Ist 'Hb<2m', wird das logarithmische Profil in der Hoehe '2m-Hb' entlang
!     dinnerhalb der P-Schicht ausgewertet.
!     Die resultierenden Transfer-Geschwindigkeiten 'tv(h,m)' sind die Kehrwerte
!     des Gesamtwiderstandes von den festen Oberflaechen (Neviau 's') bis
!     zur untersten Modellhauptflaeche (Niveau 'a').
!     Die turbulenten Diffusionskoeffizienten 'tkv(h,m)' fuer den vertikalen
!     Index 'ke1', beziehen sich aber auf den Unterrand des atmosphaerischen
!     Modells (Niveau '0').
!     Mit Hilfe der Felder 'tf(mh)' werden noch Reduktionsfaktoren der
!     Transfer-Geschwindigkeiten durch die Wirkung der L-Schicht uebertragen.
!     Diese koennen im Modul 'terra' benutzt werden, um ev. das effektive
!     'qv_s' so zu bestimmen, als gaebe es fuer fuehlbare und latente Waerme
!     unterschiedliche Parameter fuer den laminaren Transportwiderstand.
!     Zu beachten ist, dass im Falle eines vertikal vom atmosphaerischen Modell
!     aufgeloesten 'Makrobestandes' (z.B. Bebauung, Wald oder subskalige
!     Orographie) das Transferschema genauso wie im Falle eines nicht
!     aufgeloesten Bestandes angewendet wird. Allerdings beziehen sich die
!     den Bestand des Transferschemas charakterisierenden externen Parameter
!     dann auf den nicht vertikal aufgeloesten verbleibenden 'Mikrobestand',
!     der ev. allein durch niedrigen Bewuchs gebildet wird.
!     Im Transferschema eingearbeitet ist auch dei iterative Bestimmmung der
!     Rauhigkeitslaenge der Meeresoberflaeche gemaess einer modifizierten
!     Charnock-Formel, bei der die Wellenerzeugung bei verschwindenden
!     mittleren Wind mit hilfe der zur TKE ausgedrueckt wird.
!
!-------------------------------------------------------------------------------

! Declarations
!-------------------------------------------------------------------------------

!Formal Parameters:
!-------------------------------------------------------------------------------

! 0. Parameters controlling the call of 'organize_turbdiff':

LOGICAL, INTENT(IN) :: &

   lnsfdia,      & !calculation of (synoptical) near-surface variables required
   lsrflux,      & !calculation of surface flux densities in 'trubtran'

   lstfnct,      & !calculation of stability function required

   ltkeinp,      & !TKE present as input (at level k=ke1 for current time level 'ntur')
   lgz0inp,      & !gz0 present as input

   lrunscm,      & !a Single Column run (default: FALSE)
   ladsshr         !treatment of additional shear by NTCs or LLDCs active

REAL (KIND=wp), INTENT(IN) :: &

   dt_tke          !time step for the 2-nd order porgnostic variable 'tke'

INTEGER,        INTENT(IN) :: &

   iini,         & !type of initialization (0: no, 1: separate before the time loop
                   !                             , 2: within the first time step)
   ntur,         & !current  time level of 'tke' valid after  prognostic incrementation
   nprv,         & !previous time level of 'tke valid before prognostic incrementation
   ntim            !number of 'tke' time levels

INTEGER,        INTENT(IN) :: &

! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------

    nvec,         & ! number of grid points in the nproma-vector
    ke,           & ! index of the lowest main model level
    ke1,          & ! index of the lowest model half level (=ke+1)
    kcm,          & ! level index of the upper vertically-resolved canopy bound
    iblock

INTEGER,        INTENT(IN) :: &

! Start- and end-indices for the computations in the horizontal layers:
! -----------------------------------------------------------------------

    ivstart,      & ! start index in the nproma vector
    ivend           ! end index in the nproma vector

! Constants related to the earth, the coordinate system
! and the reference atmosphere:
! --------------------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: &
!
    hhl             ! height of model half levels                   ( m )

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
    l_pat,        & ! effective length scale of near-surface circulation patterns [m]
                    !  (scaling the near-surface circulation acceleration)
    l_hori,       & ! horizontal grid spacing (m)
    rlamh_fac,    & ! scaling factor for rlam_heat
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area             ( 1 )
    sai,          & ! surface area index                            ( 1 )
    urb_isa         ! urban impervious surface area                 ( 1 )

! Fields for surface values and soil|canopy model variables:
! ------------------------------------------------------------

LOGICAL, DIMENSION(:), INTENT(IN) :: &
!
    l_lake,       & ! a lake surface
    l_sice          ! an ice surface

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
!
    ps,           & ! surface pressure                              ( pa  )
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
! Atmospheric model variables:
! ---------------------------------
                    ! main-level values of:
     u,           & ! zonal wind speed       (at mass positions)    ( m/s )
     v,           & ! meridional wind speed  (at mass positions)    ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
     epr            ! exner pressure (at main levels)                (1)


REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------
!
     gz0             ! roughness length * g of the vertically not
                     ! resolved canopy                               (m2/s2)
!Achtung: Der g-Faktor ist ueberfluessig!

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(OUT) :: &
     !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
     ! velocities 'tvm' and 'tvh' makes live much easier!!               

     !turbulent transfer coefficients at the surface
     tcm,          & ! OUT: ... for momentum                         ( -- )
                     ! AUX: specific length-scale fraction           ( -- )
     tch             ! OUT: ... for scalars (heat and moisture)      ( -- )
                     ! AUX: specific length-scale fraction           ( -- )

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
     !turbulent (transfer) velocity scales at the surface 
     tvm,          & ! ... for momentum                              ( m/s)
     tvh,          & ! ... for heat and moisture                     ( m/s)

     !turbulent reduction- or amplification factors: 
          ! INP: {pure 'tkvm|h(:,ke)' from TMod}/{'tkvm|h(:,ke)' incl. LLDCs}
          ! OUT:            {Prtl-layer resist.}/{total transfer-layer resist.}   
     tfm,          & ! ... for momentum                              ( -- )
     tfh,          & ! ... for scalars                               ( -- )

     tfv  ! INP: shear-factor due to NTCs at half-level "k=ke"       ( -- )
          ! OUT: transfer-reduct.-fact. for water-vap. comp. to heat ( -- )

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(INOUT) :: &
     !reciprocal dimensionless diffusion coefficient at top of RL:
     tkr             ! Ustar/(q*Sm)_0                                ( -- )
                     ! INOUT: only, if "imode_trancnf.GE.4"
                     ! AUX: also at initialization, if "imode_trancnf.EQ.2 .OR. imode_trancnf.EQ.3"

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=wp), DIMENSION(nvec,ke1,ntim), TARGET, INTENT(INOUT) :: &
                     ! half-level values of:
     tke             ! q:=SQRT(2*TKE) with TKE='turbul. kin. energy' ( m/s )
                     ! (defined on half levels)
     !Note:
     !'tke' is the "turbulent velocity" (in m/s) and NOT the (mass-density) of turb. kin. energy,
     ! which has the dimension m2/s2!
     !In case of "ntim=1", the actual parameter for 'tke' may be a 2-dim. array for a fix time level.

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
                     ! half-level values of:
     tkvm,         & ! turbulent diffusion coefficient for momentum  (m2/s )
     tkvh            ! turbulent diffusion coefficient for heat      (m2/s )
                     ! (and other scalars)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
     rcld            ! standard deviation of local super-saturation (SDSS)
                     !  at MAIN levels including the lower boundary  (---)
                     ! AUX: cloud-cover at zero-level (as output of SUB 'adjust_satur_equil'  
                     !                                 and input of SUB 'solve_turb_budgets')       

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
                     ! half-level values of:
     hdef2,        & ! horizontal deformation square at half levels  ( 1/s2 )
     dwdx,         & ! zonal      derivative of vertical wind  ,,    ( 1/s )
     dwdy            ! meridional derivative of vertical wind  ,,    ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
                     ! half-level values of:
     tketens         ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(OUT) :: &
                     ! half-level values of:
     edr             ! eddy dissipation rate of TKE (EDR)            (m2/s3)

REAL (KIND=wp), DIMENSION(:), OPTIONAL, INTENT(OUT) :: &
!
! Diagnostic near surface variables:
! -----------------------------------------------
!
     t_2m,         & ! temperature in 2m                             (  K  )
     qv_2m,        & ! specific water vapor content in 2m            (kg/kg)
     td_2m,        & ! dew-point in 2m                               (  K  )
     rh_2m,        & ! relative humidity in 2m                       (  %  )
     u_10m,        & ! zonal wind in 10m                             ( m/s )
     v_10m           ! meridional wind in 10m                        ( m/s )

REAL (KIND=wp), DIMENSION(:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
     shfl_s,       & ! sensible heat flux at the surface             (W/m2)    (positive downward)
     qvfl_s,       & ! water vapor   flux at the surface             (kg/m2/s) (positive downward)
     umfl_s,       & ! u-momentum flux at the surface                (N/m2)    (positive downward)
     vmfl_s          ! v-momentum flux at the surface                (N/m2)    (positive downward)

LOGICAL, OPTIONAL, INTENT(IN) :: lacc
INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

LOGICAL :: lzacc
INTEGER :: acc_async_queue

INTEGER            :: my_cart_id, my_thrd_id

!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

INTEGER ::      &
    i, k,       & !horizontaler und vertikaler Laufindex
    k1,k2,ks,   & !spezifische Level-Indices
    n,          & !Index fuer diverse Schleifen
!
    nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes (wird zwischen Iterationen zu 'ntur')
    it_durch,   & !Durchgangsindex der Iterationen
    it_start,   & !Startindex der Iterationen
!
    imode_syndiag != ABS(imode_2m_diag)

REAL (KIND=wp) :: &
    fr_tke,             & ! z1/dt_tke
    wert, val1, val2,   & ! Platzhalter fuer beliebige Zwischenergebnisse
    fakt, fac_m, fac_h, & !  ,,         ,,     ,,      Faktoren

!   Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 (fh2,fm2):
    fh2,fm2, &

!   Platzh. fuer horiz. Geschw.-Komponenten, bel. Geschw. und Druck:
    vel1,vel2,velo, patm, &

!   Platzh. fuer Hoehendifferenzen und Laengaenskalen:
    dh,l_turb,lh,lm,z_surf,len1,len2, &
    fr_sd_h, &     !dimensionless resistance for scalars between the surface and the synoptic zero-level
    h_2m, h_10m, & !level heights (equal 2m and 10m)
    a_2m, a_10m, & !turbulent distance of 2m- and 10m-level with respect to diag. roughness
    a_atm, &       !turbulent distance of the atmosp. level with respect to land-use roughness
    h_atm, &       !mid-level heigth of transfer layer (atm. level)
    edgrav,  &     ! 1/grav

!   for TERRA_URB:
    zkbmo_dia, zkbmo_urb, zustar, & ! inverse Stanton number

!   Sonstiges:
    ren_m,ren_h, & !specific Re-numbers at top of land-use roughness
    g_z0_ice, g_len_min, g_alpha1_con_m, & !used for calculating sea-surface roughness
    edprfsecu  , & !"1/prfsecu" (out of ]0; 1[)
    xf             !scaling factor representing the volume-height ratio of the lowest atmospheric 
                   ! full and half level

LOGICAL            ::  &
    lini,      & !initialization required
    lgz0ini,   & ! initialization of roughness lenght over water and ice
    lnswindia, & !diagnostics of near-surface wind active in in this SUB 
    lnswinamp    !amplification of near-surface wind due to NTCs or LLDCs active

! Local arrays:

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!   pointer for variable layers:
    vel1_2d  (:),      &
    vel2_2d  (:),      &
    ta_2d    (:),      &
    qda_2d   (:),      &
!
    g_tet    (:),      &
    g_vap    (:),      &
    qsat_dT  (:),      &
!
    epr_2d   (:),      &
!
    l_tur_z0 (:),      &
!
    z_mom_tot(:),      & ! total roughness length
    a_atm_tot(:),      & ! total turb. dist. (of the adapted virtual atm. level)
    a_atm_mod(:),      & ! modified total turb. dist.
    h_atm_mod(:),      & ! modified heigth of the atm. level
    prf_ren_m(:),      & ! log. profile Re-number for momentum
    prf_ren_h(:),      & ! log. profile Re-number for scalars
!
    prss     (:,:),    & ! near-surface pressure (Pa)
    tmps     (:,:),    & ! near-surface temperature-varible (K)
    vaps     (:,:),    & ! near-surface humidity-variable
    liqs     (:,:),    & ! near-surface liquid water content
!
    tvt      (:,:),    & ! surface tendency of turbulent velocity scale
    ediss    (:,:)       ! surface eddy-dissipation rate

REAL (KIND=wp), TARGET    :: &
  ! targets of used pointers
  diss_tar    (nvec,ke1:ke1), & ! eddy dissipation rate (m2/s3)
  tketens_tar (nvec,ke1:ke1), & ! turbulent transport of SQRT(TKE) (m/s2)

  ! internal atmospheric variables
  len_scale   (nvec,ke1:ke1), & ! turbulent length-scale (m)
  l_scal      (nvec),         & ! reduced maximal turbulent length scale due to 
                                !   horizontal grid spacing (m)

  fc_min      (nvec),         & ! minimal value for TKE-forcing (1/s2)

  rhon        (nvec,ke1:ke1), & ! boundary level air density (at surface level) (Kg/m3)
  frh         (nvec,ke1:ke1), & ! thermal forcing (1/s2) or thermal acceleration (m/s2)
  frm         (nvec,ke1:ke1), & ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)

  zaux        (nvec,ke1:ke1,ndim), & 
                                ! auxilary array containing thermodynamical properties on boundary levels:
                                ! (1:ex_fakt, 2:cp_fakt, 3:dQs/dT, 4:g_tet l, 5:g_vap)
  zvari       (nvec,ke-1:ke1,0:ndim), & 
                                ! set of variables used in the turbulent 2-nd order equations
                                ! and later their effective vertical gradients
                                ! 'zvari(:,:,0)' is reserved for half-level pressure or "circulation acceleration"
                                ! (finally, 'zvari(:,ke1,tet|vap)' also contains T2m and qv_2m)

  rcls        (nvec,ke1),     & ! double for standard deviation of the local super-saturation (SDSS)
                                ! ('rcls(:,ke1)' provides 2m-SDSS to, and receives 2m-cl_cov from, SUB 'adjust_satur_equil')

  tl_s_2d     (nvec),         & ! surface level value of conserved temperature 
                                !   (liquid water temperature) (K)
  qt_s_2d     (nvec),         & ! surface level value of conserved humidity    
                                !   (total water)
  vel_2d      (nvec,ke-1:ke), & ! wind speed (m/s) at the two lowest full model levels

  velmin      (nvec),         & ! modified 'vel_min' used for tuning corrections (m/s)
                                ! (hyper-parameterizations)
  ratsea      (nvec),         & ! modified 'rat_sea' used for tuning corrections 
                                ! (hyper-parameterizations)

  ! internal variables for the resistance model
  hk_2d       (nvec),         & ! mid level height above ground belonging to 'k_2d' (m)
  hk1_2d      (nvec),         & ! mid level height above ground of the previous layer (below) (m)
  hk2_2d      (nvec),         & ! 'hk1_2d, hk1_2d and hk2_2d' are also used a shelf for any level height

  a_atm_2d    (nvec),         & ! turbulent distance of the lowermost full model level

  h_top_2d    (nvec),         & ! boundary level height of transfer layer (top  level)
  h_atm_2d    (nvec),         & ! mid      level heigth of transfer layer (atm. level)
  h_can_2d    (nvec),         & ! effective canopy height (m) used as depth of the R-layer, applied
                                !  for calculation of the roughness layer resistance of momentum or 
                                !  for determination of the synoptic 2m-level in case of a diagnostic
                                !  exponential R-layer profile (at "itype_2m_diag.EQ.2")

  edh         (nvec),         & ! reciprocal of any layer depth
  z0m_2d      (nvec),         & ! mean  roughness length
  z0d_2d      (nvec),         & ! diag. roughness length (for the SYNOP lawn)
  z2m_2d      (nvec),         & ! height of 2m  level (above the surface) or total roughn. length for momentum
  z10m_2d     (nvec),         & ! height of 10m level (above the surface)

  rat_m_2d    (nvec),         & ! any surface layer ratio for momentum 
                                !   (like Re-number or Stability factor)
  rat_h_2d    (nvec),         & ! any surface layer ratio for scalars  
                                !   (like Re-number or Stability factor)
  fac_h_2d    (nvec),         & ! surface layer profile factor for scalars
  fac_m_2d    (nvec),         & ! surface layer profile factor for momentum

  frc_2d      (nvec),         & ! saved amplfification factor for wind-shear forcing

  val_m       (nvec),         & ! any specific value for momentum (Re-numer or resistance length)
  val_h       (nvec),         & ! any specific value for scalars  (Re-numer or resistance length) 

  dz_sg_m           ,         & ! laminar resistance lenght for momentum (m)
  dz_sg_h     (nvec),         & ! laminar resistance length for scalars (m)
  dz_g0_h     (nvec),         & ! turbulent roughness layer resistance length for scalars (m)
  dz_0a_m     (nvec),         & ! turbulent Prandtl-layer resistance length for momentum (m)
  dz_0a_h     (nvec),         & ! turbulent Prandtl-layer resistance length for scalars (m)
  dz_sa_m           ,         & ! otal transfer resistance length for momentum (m)
  dz_sa_h     (nvec),         & ! total transfer resistance length for scalars (m)
  dz_s0_m     (nvec),         & ! total roughness layer resistance lenght for momentum (m)
  dz_s0_h     (nvec)            ! total roughness layer resistance lenght for scalars (m)

REAL (KIND=wp) ::             &
  grad        (nvec,nmvar)      ! any vertical gradient

INTEGER        ::             &
  k_2d        (nvec)            ! index field of the upper level index to be used
                                !   for near surface diagn.
LOGICAL        ::   ldebug = .FALSE.

!---- End of header -----------------------------------------------------------

  CALL set_acc_host_or_device(lzacc, lacc)

  IF(PRESENT(opt_acc_async_queue)) THEN
    acc_async_queue = opt_acc_async_queue
  ELSE
    acc_async_queue = 1
  ENDIF

!==============================================================================
! Begin subroutine turbtran
!------------------------------------------------------------------------------

! 1)  Vorbereitungen:

  !GPU data region of all local variables except pointers which are set later on
  !$ACC DATA &
! local variables
  !$ACC   CREATE(diss_tar, tketens_tar, len_scale, l_scal, fc_min) &
  !$ACC   CREATE(rhon, frh, frm, zaux, zvari, rcls) &
  !$ACC   CREATE(tl_s_2d, qt_s_2d, vel_2d, velmin, ratsea) &
  !$ACC   CREATE(hk_2d, hk1_2d, hk2_2d, h_top_2d, h_atm_2d, a_atm_2d, h_can_2d) &
  !$ACC   CREATE(edh, val_m, val_h, z0m_2d, z0d_2d, z2m_2d) &
  !$ACC   CREATE(z10m_2d, rat_m_2d, rat_h_2d, fac_h_2d, fac_m_2d) &
  !$ACC   CREATE(frc_2d, dz_s0_m, dz_sg_h, dz_g0_h) &
  !$ACC   CREATE(dz_0a_m, dz_0a_h, dz_sa_h, dz_s0_h, grad) &
  !$ACC   CREATE(k_2d) &
  !$ACC   COPYIN(ivend) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)


  ! take care that all pointers have a target
  IF (PRESENT(edr)) THEN
     ediss => edr
  ELSE
     ediss => diss_tar
  END IF

  IF (PRESENT(tketens)) THEN
     tvt => tketens
  ELSE
     tvt => tketens_tar
     ! This loop can't be OpenACC-collapsed yet because ivend is on the GPU
     !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(tvt) ASYNC(acc_async_queue) IF(lzacc)
     !$ACC LOOP SEQ
     DO k=ke1, ke1
        !$ACC LOOP GANG VECTOR
        DO i=ivstart, ivend
           tvt(i,k) = z0
        END DO
     END DO
     !$ACC END PARALLEL
  END IF

  prss(1:,ke-1:) => zvari(:,:,0)    ! near-surface pressure (Pa)
  tmps(1:,ke-1:) => zvari(:,:,tet)  ! near-surface temperature-variable (K)
  vaps(1:,ke-1:) => zvari(:,:,vap)  ! near-surface humidity-variable
  liqs(1:,ke-1:) => zvari(:,:,liq)  ! near-surface liquid water content

  vel1_2d => zvari(:,ke ,u_m)
  vel2_2d => zvari(:,ke ,v_m)
  
  ta_2d   => zvari(:,ke ,tet_l)
  qda_2d  => zvari(:,ke ,h2o_g)

  epr_2d  => zaux(:,ke1,1)
  qsat_dT => zaux(:,ke1,3)
  g_tet   => zaux(:,ke1,4)
  g_vap   => zaux(:,ke1,5)

  l_tur_z0 => len_scale(:,ke1)

  prf_ren_h => val_h
  z_mom_tot => z0m_2d   !total roughness length                     equals local roughness length
  a_atm_tot => a_atm_2d !total turbulent distance                   equals local turbulent dist.
  a_atm_mod => a_atm_2d !modified total turbulent distance          equals local turbulent dist.
  h_atm_mod => h_atm_2d !modified heigth of the atmospheric level   equals local one
  prf_ren_m => val_h    !logarithmic profile Re-number for momentum equals that one for scalars

  !SUB-arguments and pointers explicitly operating in this code: 
  !$ACC DATA PRESENT(fr_land, t_g, l_lake, l_sice) &
  !$ACC   PRESENT(prss, tmps, vaps, liqs) &
  !$ACC   PRESENT(z_mom_tot, a_atm_tot, a_atm_mod, h_atm_mod, prf_ren_m, prf_ren_h) &
  !$ACC   PRESENT(hhl, epr_2d, u, v, t) &
  !$ACC   PRESENT(ediss, l_pat) &
  !$ACC   PRESENT(epr, tke, tkvm, tkvh, gz0, tkr) &
  !$ACC   PRESENT(l_tur_z0, ps, sai, urb_isa, rlamh_fac) &
  !$ACC   PRESENT(tfm, tfh, tfv, qv_s, qv, qc) &
  !$ACC   PRESENT(dwdx, dwdy, hdef2, g_tet) &
  !$ACC   PRESENT(g_vap, tcm, tch, z0d_2d, shfl_s) &
  !$ACC   PRESENT(rcld, qsat_dT, qvfl_s, umfl_s, vmfl_s) &
  !$ACC   PRESENT(edr, tvh, tvm, t_2m) &
  !$ACC   PRESENT(qv_2m, qda_2d, ta_2d) &
  !$ACC   PRESENT(v_10m, u_10m, vel1_2d, vel2_2d, rh_2m, td_2m) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)

!-------------------------------------------------------------------------------
  CALL turb_setup (i_st=ivstart, i_en=ivend, k_st=ke, k_en=ke1, &
                   iini=iini, dt_tke=dt_tke, nprv=nprv, l_hori=l_hori, &
                   ps=ps, t_g=t_g, qv_s=qv_s, qc_a=qc(:,ke), &
                   lini=lini, it_start=it_start, nvor=nvor, fr_tke=fr_tke,  &
                   l_scal=l_scal, fc_min=fc_min, &
                   prss=prss(:,ke1), tmps=tmps(:,ke1), vaps=vaps(:,ke1), liqs=liqs(:,ke1), rcld=rcld, &
                   lacc=lzacc, opt_acc_async_queue=acc_async_queue)

#ifdef ICON_USE_CUDA_GRAPH
   IF (lzacc .AND. lini .AND. lcuda_graph_turb_tran ) THEN
      CALL finish ('turbtran', 'initialization is not supported when capturing a graph with OpenACC')
   END IF
#endif
!-------------------------------------------------------------------------------

  my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
  my_thrd_id = omp_get_thread_num()
#endif

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren:

      ! Berechnung einiger Hilfsgroessen und Initialisierung der Diffusionskoeff.:

      ! Unterste Hauptflaeche halbiert den Abstand zur untersten Nebenflaeche:
      xf=z2
      xf=z1/xf

      ! Hoehe des 2m- und 10m-Niveaus:
      h_2m  = z2
      h_10m = z10
 
      edprfsecu = z1/prfsecu !used for limitation of the profile funtion
      edgrav = z1/grav

      ! Fixed parameter used for calculating sea-surface roughness:
      g_z0_ice=grav*z0_ice; g_len_min=grav*len_min; g_alpha1_con_m=grav*alpha1*con_m

      ks=ke

      lnswindia=(lnsfdia .AND. .NOT.lfreeslip) !diagnostics of near-surface wind active in in this SUB
      lnswinamp=(ladsshr .AND. imode_nsf_wind.EQ.2) !amplific. of near-surface wind by NTCs or LLDCs active

      imode_syndiag=ABS(imode_2m_diag)
 
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ! Dicke der Modell-Prandtl-Schicht
         h_top_2d(i) = hhl(i,ke)-hhl(i,ke1)
         h_atm_2d(i) = h_top_2d(i)*xf

         ! Surface-Exner-pressure:
         epr_2d(i) = zexner(ps(i))
 
      END DO

      !$ACC LOOP SEQ
      DO k=ks, ke
!DIR$ IVDEP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            vel_2d(i,k) = MAX( vel_min, SQRT( u(i,k)**2+v(i,k)**2 ) ) !wind speed
         END DO
      END DO
  
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ratsea(i) = rat_sea
      END DO

      IF (imode_rat_sea.EQ.2) THEN
!<Tuning
!---------------------------------------------------------------------------------------
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            IF (t_g(i) - t(i,ke) > 8._wp) THEN
               ! Increase rat_sea for very large temperature differences between water and adjacent air
               ! in order to reduce excessive peaks in latent heat flux

               ratsea(i) = rat_sea*(1._wp + 0.05_wp*(t_g(i) - t(i,ke) - 8._wp))
            END IF
         END DO
      END IF

      IF (imode_vel_min.EQ.2) THEN
!<Tuning    
!---------------------------------------------------------------------------------------
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            ! stability-dependent minimum velocity serving as lower limit on surface TKE
            ! (parameterizes small-scale circulations developing over a strongly heated surface;
            ! tuned to get 0.75 m/s when the land surface is at least 7.5 K warmer than the air in the
            ! lowest model level; nothing is set over water because this turned out to induce
            ! detrimental effects in NH winter)

            velmin(i) = MAX( vel_min, MIN(0.75_wp, fr_land(i)*(t_g(i)/epr_2d(i) - t(i,ke)/epr(i,ke))/ &
                        LOG(2.e3_wp*h_atm_2d(i))) )
         END DO
      END IF
!---------------------------------------------------------------------------------------
!>Tuning: This kind of correction should be substituded by a less ad-hoc approach.

      IF (lini) THEN !only for initialization

         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(lgz0ini, l_turb, dh, vel1, vel2, fm2, fh2, fakt, lm, lh, wert, val1, val2)
         DO i=ivstart, ivend

            lgz0ini=(.NOT.lgz0inp .AND. fr_land(i) <= z1d2)
            !Note: This definition of a non-land surface is now in line with the ICON-definition 
            !       using "frlnd_thrhld=z1d2"

            IF (imode_trancnf.GE.2 .OR. (lgz0ini .AND. .NOT.l_sice(i))) THEN !initial 'tkvm|h(:ke)' required
!              Einfachste Schaetzung der Schubspannung als Impusls-
!              flussdichte durch die Nebenflaeche ke mit Hilfe
!              einer diagnostischen TKE ohne Beruecksichtigung von
!              Feuchte-Effekten und mit neuchtralen Stabilitaets-
!              funktion:

               l_turb=h_top_2d(i) !approx. turb. length scale at level ke

               l_turb=akt*MAX( len_min, l_turb/( z1+l_turb/l_scal(i) ) )

               dh=z1d2*(hhl(i,ke-1)-hhl(i,ke1))

               vel1=u(i,ke-1)
               vel2=u(i,ke  )
               grad(i,u_m)=(vel1-vel2)/dh
   
               vel1=v(i,ke-1)
               vel2=v(i,ke  )
               grad(i,v_m)=(vel1-vel2)/dh

               grad(i,tet_l)=(t(i,ke-1)-t(i,ke))/dh + tet_g

               fm2=MAX( grad(i,u_m)**2+grad(i,v_m)**2, fc_min(i) )
               fh2=grav*grad(i,tet_l)/t(i,ke)

               ! Vereinfachte Loesung mit Rf=Ri:
               IF (fh2.GE.(z1-rim)*fm2) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  lm=l_turb*(sm_0-(a_6+a_3)*fakt)
                  lh=lm
               ELSE
                  fakt=fh2/(fm2-fh2)
                  lm=l_turb*(sm_0-(a_6+a_3)*fakt)
                  lh=l_turb*(sh_0-a_5*fakt)
               END IF

               val1=lm*fm2; val2=lh*fh2
               wert=MAX( val1-val2, rim*val1 )

               IF (ltkeinp) THEN
                  tke(i,ke,nvor)=tke(i,ke,ntur)
               ELSE
                  tke(i,ke,nvor)=SQRT(d_m*l_turb*wert)
               END IF

               val1=MAX ( con_m, tkmmin ); tkvm(i,ke)=lm*tke(i,ke,nvor)
               val2=MAX ( con_h, tkhmin ); tkvh(i,ke)=lh*tke(i,ke,nvor)

               IF (imode_tkemini.EQ.2 .OR. rsur_sher.GT.z0) THEN !adaptation of TKE and TMod. to lower limits
                  tke(i,ke,nvor)=tke(i,ke,nvor)*MAX( z1, val2/tkvh(i,ke) ) !adapted 'tke'
               END IF

               IF (imode_suradap.GE.1) THEN  !using 'tkvm(:,ke)' relieved from artific. drag impact by LLDCs
                  wert=tkvm(i,ke)/tkvh(i,ke) !turbulent Prandtl-number before application of lower limits
                  tkvh(i,ke)=MAX(val2, tkvh(i,ke)) !'tkvh(:,ke)' with lower limit
                  tkvm(i,ke)=wert*tkvh(i,ke)       !'tkvm(:,ke)' with adapted lower limit
               ELSE !employing lower limit of 'tkvm(:,ke)' without any correction
                  tkvh(i,ke)=MAX(val2, tkvh(i,ke)) !'tkvh(:,ke)' with lower limit
                  tkvm(i,ke)=MAX(val1, tkvm(i,ke)) !'tkvm(:,ke)' with full lower limit 
               END IF 
                     
               val2=MAX( epsi, tkvm(i,ke)*SQRT(fm2) ) !estimate of Ustar**2
               val1=SQRT(val2) !Ustar

            ELSE !needed for vectorization
               l_turb = z0
               val1   = z0
               val2   = z0
            END IF !initial 'tkvm|h(:ke)' required

            tfm(i)=z1 !'tkvm(:,ke)' is already adapted
            tfh(i)=z1 !'tkhm(:,ke)' is already adapted
            tfv(i)=z1 !no amplif.-factor of wind-shear due to NTCs at initialization

            IF (lgz0ini) THEN !initialization of roughness length for water- or ice-covered surface:
               IF ( l_sice(i) ) THEN !ice-covered surface
                  gz0(i)=g_z0_ice
               ELSE !water-covered surface
                  fakt=MERGE( alpha0_char(   & !extended value at "imode_charpar>1" apart from lakes
                                           vel_2d(i,ke) ), & !as a function of vel_ke (substiting vel_10m at initializtion)
                              MERGE( alpha0, & !standard value at "imode_charpar=1"
                                     0.1_wp, & !enhanced value at "imode_charpar>1" over lakes (with non-equlib. wave spectrum)
                                     imode_charpar.EQ.1 ), &
                              imode_charpar.GT.1 .AND. .NOT.l_lake(i) )
                  gz0(i)=MAX( g_len_min, fakt*val2+g_alpha1_con_m/val1 )
               END IF
            END IF

            IF (imode_trancnf.GE.2) THEN !new version of init. using estimated Ustar
               tkr(i)=l_turb*val1                          !l_0*Ustar

               rat_m_2d(i)= tkr(i)/tkvm(i,ke)              !Ustar/(q*Sm)_p
               rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke)*sm_0) !Ustar/(q*Sh)_p*Sh(0)/Sm(0)
            END IF

         END DO

      END IF !only for initialization

!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         z0m_2d(i)  = gz0(i)*edgrav !mean roughness-length
         l_tur_z0(i)= akt*z0m_2d(i) !turbulent length scale
         tcm(i) = z0m_2d(i)/(h_top_2d(i)+z0m_2d(i)) !unspecific length scale fraction
         tch(i) = tcm(i)       !default setting with unspecific length scale fraction

         a_atm_2d(i)=h_atm_2d(i)+z0m_2d(i) !turbulent distance of the lowermost full model level

         !Note: The additional 'grav'-factor in 'gz0' is obsolete and should be removed!
      END DO

      IF (lini) THEN !only for initialization

         IF (imode_trancnf.GE.2) THEN !new version of init. using estimated Ustar
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tkvm(i,ke1)=tkvm(i,ke)*tcm(i)*(tcm(i)+(z1-tcm(i))*rat_m_2d(i))
               tkvh(i,ke1)=tkvh(i,ke)*tch(i)*(tch(i)+(z1-tch(i))*rat_h_2d(i))

               tkr(i)=tkr(i)/tkvm(i,ke1) !Ustar/(q*Sm)_0
            END DO
         ELSE !old version of init. using laminar diff. coeffs. 
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tkvm(i,ke) =con_m; tkvh(i,ke) =con_h
               tkvm(i,ke1)=con_m; tkvh(i,ke1)=con_h
            END DO
         END IF    

      ELSEIF (imode_trancnf.LT.4 .AND. imode_suradap.GT.0) THEN !not for initialization and profile-factors
                                                                ! calculated by means of an upper node
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            IF (imode_suradap.EQ.1 .OR. imode_suradap.EQ.3) THEN !using 'tkvm(:,ke)' relieved from 
                                                                 ! artific. drag impact by LLDCs
               tcm(i)=tcm(i)*(tfm(i)/tfh(i)) !length-scale fraction including related reduction factor
            ELSEIF (imode_suradap.EQ.2) THEN !using 'tkvm|h(:,ke)' relieved from any LLDC-impact
               tcm(i)=tcm(i)*tfm(i) !length-scale fraction including full reduction factor for 'tkvm(:,ke)'
               tch(i)=tch(i)*tfh(i) !length-scale fraction including full reduction factor for 'tkvh(:,ke)'
            END IF
         END DO
      END IF

      IF (imode_trancnf.GE.4) THEN
         !Calculation the profile-factors without an upper node of diffusion coefficients (above the RL)
         !but based on previous values of Ustar and diffusion-coefficients (at the top of the RL):
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            rat_m_2d(i)=tkr(i)                                       !Ustar/(q*Sm)_0
            rat_h_2d(i)=tkr(i)*(tkvm(i,ke1)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(Sh(0)/Sm(0))
         END DO
      ELSEIF (imode_trancnf.GE.2) THEN
         !Profile-factors by using the previous diffusion coefficients
         !without a laminar correction, but still based on the upper node:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            rat_m_2d(i)=tcm(i)*tkvm(i,ke)/tkvm(i,ke1) !(q*Sm)_p/(q*Sm)_0
            rat_h_2d(i)=tch(i)*tkvh(i,ke)/tkvh(i,ke1) !(q*Sh)_p/(q*Sh)_0
         END DO 
      END IF

      IF (ladsshr) THEN !treatment of additional shear by NTCs or LLDCS active
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            frc_2d(i)=tfv(i) !saving the amplification-factor for wind-shear forcing
            IF (imode_suradap.EQ.3) THEN !employing shear contrib. by LLDCs at surface-level "k=ke1"
               frc_2d(i)=frc_2d(i)/tfh(i)**2 !total amplific.-factor for wind-shear forcing
            END IF 
!           frc_2d(i)=rsur_sher*frc_2d(i) !scaled overall amplification factor for wind shear forcing
            frc_2d(i)=rsur_sher*(frc_2d(i)-z1)+z1 !scaled overall amplification factor for wind shear forcing
         END DO
      END IF

      !$ACC END PARALLEL

! 4)  Berechnung der Transfer-Geschwindigkeiten:

      IF (lprfcor) THEN   
         ks=ke-1
      ELSE
         ks=ke
      END IF

!----------------------------------------------
      DO it_durch=it_start, it_end !Iterationen
!----------------------------------------------

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z_surf, fakt, ren_m, ren_h, dz_sg_m)
         DO i=ivstart, ivend
            z_surf= z0m_2d(i)/sai(i) !effektive Rauhigkeitslaenge

            ! Laminare Korrektur der Diffusionskoeffizienten:
            tkvm(i,ke1)=MAX( con_m, tkvm(i,ke1) )
            tkvh(i,ke1)=MAX( con_h, tkvh(i,ke1) )

            fakt=z1+(z1-REAL(NINT(fr_land(i)),wp))*(ratsea(i)-z1)

            ren_m=tkvm(i,ke1)/con_m
            ren_h=tkvh(i,ke1)/con_h

            ! Effective resistance length of local R-layer for scalars:

            dz_sg_h(i)=fakt*rlam_heat*rlamh_fac(i)*z_surf*(ren_h/ren_m) & !through laminar layer
                           *MERGE( rat_glac, 1._wp, & !particul. scaling by 'rat_glac'
                                   gz0(i)<0.01_wp .AND. fr_land(i)>z1d2 ) !over glaciers
            !Note: This definition of a land surface is now in line with the ICON-definition 
            !       using "frlnd_thrhld=z1d2"

            dz_g0_h(i)=z_surf*LOG(ren_m) !through turbulent R-layer

            dz_s0_h(i)=dz_sg_h(i)+dz_g0_h(i) !through full R-layer

            ! Effective height of the R-layer

            IF (rlam_mom.GT.z0 .OR. itype_2m_diag.EQ.2 ) THEN !R-height required
               h_can_2d(i)=rat_can*MERGE( sai(i)*z0m_2d(i), &
                                          dz_s0_h(i)*LOG(dz_s0_h(i)/dz_sg_h(i)), &
                                          dz_sg_h(i).EQ.z0 )
            END IF

            ! Effective resistance length of local R-layer for momentum:

            IF (rlam_mom.GT.z0) THEN !including R-layer resistance for momentum
               dz_sg_m=rlam_mom*z_surf !through laminar layer
               wert=z1d2*dz_sg_m
               dz_s0_m(i)=wert+SQRT(wert**2+h_can_2d(i)*dz_sg_m) !through full R-layer
            ELSE
               dz_s0_m(i) = z0 !no R-layer resistance for momentum
            END IF
         END DO 

!--------------------------------------------------------------------------
         IF (lterra_urb .AND. (.NOT. itype_kbmo == 1)) THEN
!DIR$ IVDEP
           !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(zkbmo_dia, zustar, zkbmo_urb)
           DO i=ivstart, ivend
              IF (urb_isa(i) > 0.0_wp) THEN
                zkbmo_dia    = dz_s0_h(i)/z0m_2d(i)
                zustar       = SQRT(tvm(i)*vel_2d(i,ke))
                IF (itype_kbmo == 2) THEN
                  ! Brutsaert Kanda parameterisation for bluff-body elements
                  ! for some reason, it doesn't withstand zero values of ustar
                  zkbmo_urb = MAX(0.1_wp, 1.29_wp * (z0m_2d(i)*zustar/con_m)**0.25_wp - z2)
                ELSE !(itype_kbmo == 3)
                  ! Zilitinkevich
                  zkbmo_urb = MAX(0.1_wp, 0.13_wp * (z0m_2d(i)*zustar/con_m)**0.45_wp)
                END IF

                dz_s0_h(i) = (urb_isa(i)*zkbmo_urb + (1.0_wp - urb_isa(i))*zkbmo_dia) * z0m_2d(i)
              ENDIF
           END DO
         END IF
!--------------------------------------------------------------------------

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
         DO i=ivstart, ivend

!           Profilfakoren der turbulenten Prandtl-Schicht:

            fakt=z0m_2d(i)/h_top_2d(i)

            IF (imode_trancnf.LT.4) THEN
               !Profile-factors by employing previous values of the diffusion-coefficients
               !at the top fo the roughness-layer (0) and also at the upper bound of the
               !lowest atm. model layer (p) as an upper node:


               IF (imode_trancnf.EQ.1) THEN !first version
                  !Profile factors by using the previous diffusion coefficients
                  !including a laminar correction:
                  rat_m_2d(i)=tcm(i)*tkvm(i,ke)/tkvm(i,ke1)
                  rat_h_2d(i)=tch(i)*tkvh(i,ke)/tkvh(i,ke1)
               END IF

               rat_m_2d(i)=MIN( edprfsecu, MAX( prfsecu, rat_m_2d(i) ) ) !limitted (q*Sm)_p/(q*Sm)_0
               rat_h_2d(i)=MIN( edprfsecu, MAX( prfsecu, rat_h_2d(i) ) ) !limitted (q*Sh)_p/(q*Sh)_0

               fac_m_2d(i)=(rat_m_2d(i)-z1)*fakt !non-stab. profile-factor for momentum
               fac_h_2d(i)=(rat_h_2d(i)-z1)*fakt !non-stab. profile-factor for scalars
             
            ELSE !Profile-factors without using the upper node

               fac_m_2d(i)=z1-rat_m_2d(i) !profile-factor for momentum
               fac_h_2d(i)=z1-rat_h_2d(i) !profile-factor for scalars
            END IF

         END DO

!        Effektive Widerstandslaengen der turb. Prandtl-Schicht:

        ! Preparations with regard to the hyperbolic interpolation of profile-function:

         IF (imode_trancnf.GE.3) THEN !hyperbolic interpolation of profile-function for stable strat.
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               prf_ren_h(i)=LOG(a_atm_2d(i)/z0m_2d(i)) !log. profile Re-number for scalars
            END DO
         END IF

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(a_atm, z_surf, h_atm, fac_h, fac_m)
         DO i=ivstart, ivend
            a_atm=a_atm_2d(i)
            z_surf=z0m_2d(i)
            h_atm=h_atm_2d(i)

            fac_m=fac_m_2d(i)

            IF (fac_m.GE.z0 .OR. imode_trancnf.LT.3) THEN
               !non-stable stratfication or based on linear interpolation of profile-function
               !for the velocity scale (q*Sm):
               dz_0a_m(i)=z_surf*MERGE( (z_surf*h_atm)/(a_atm_tot(i)*z_mom_tot(i)), &
                                        LOG( a_atm_mod(i)/(z_surf+fac_m*h_atm_mod(i)) )/(z1-fac_m), &
                                        fac_m.EQ.z1 )
            ELSEIF (imode_trancnf.EQ.3) THEN 
               !based on hyperbolic interpolation of profile-function (q*Sm) for stable stratification,
               ! at which the upper node for diffusion-coefficients is used:
               fac_m_2d(i)=-fac_m/rat_m_2d(i) !transformed profile-factor for stable stratification
               dz_0a_m(i)=z_surf*(z1-fac_m)*prf_ren_m(i)+fac_m*h_atm
            ELSE !(imode_trancnf.EQ.4): without using the upper node for diff.-coefs.
               dz_0a_m(i)=(z_surf*prf_ren_m(i)-fac_m*h_atm)/(z1-fac_m)
            END IF

            fac_h=fac_h_2d(i)

            IF (fac_h.GE.z0 .OR. imode_trancnf.LT.3) THEN
               dz_0a_h(i)=z_surf*MERGE( h_atm/a_atm, &
                                        LOG( a_atm/(z_surf+fac_h*h_atm) )/(z1-fac_h), &
                                        fac_h.EQ.z1 )
            ELSEIF (imode_trancnf.EQ.3) THEN !the upper node for diffusion-coefficients is used
               fac_h_2d(i)=-fac_h/rat_h_2d(i) !transformed profile-factor for stable stratification
               dz_0a_h(i)=z_surf*(z1-fac_h)*prf_ren_h(i)+fac_h*h_atm
            ELSE !(imode_trancnf.EQ.4): without using the upper node for diff.-coefs.
               dz_0a_h(i)=(z_surf*prf_ren_h(i)-fac_h*h_atm)/(z1-fac_h)
            END IF
         END DO

         ! Kombination of resistance-length values for momentum and scalars:

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(dz_sa_m, wert)
         DO i=ivstart, ivend

!           Effektive Widerstandslaengen von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):

            ! total resistance-path of the transfer-layer for momentum:
            dz_sa_m    = dz_0a_m(i) + dz_s0_m(i)

            ! total resistance-path of the transfer-layer for scalars:
            dz_sa_h(i) = dz_s0_h(i) + dz_0a_h(i) 

!           Reduktionsfaktoren fuer die Bestandesschicht incl. lam. Grenzschicht:

            tfm(i)=dz_0a_m(i)/dz_sa_m    !for momentum
            tfh(i)=dz_0a_h(i)/dz_sa_h(i) !for scalars

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den Faktor 'rat_lam'
!           gegenueber fuehlbarer Waerme vergroesserten laminaren Transportwiderstandes:

            tfv(i)=z1/(z1+(rat_lam-z1)*dz_sg_h(i)/dz_sa_h(i))

            !Note: 
            !So far, this reduction factor for evaporation is only applied in 'terra',
            ! hence, it is not affecting evaporation of the open sea or sea-ice!
         END DO

!        Berechnung der Erhaltungsgroessen in der Prandtl-Schicht:

         IF (icldm_tran.EQ.-1 .OR. ilow_def_cond.EQ.2) THEN
            !conserved values at the rigid surface are temperature and humidity
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tl_s_2d(i)=t_g(i); qt_s_2d(i)=qv_s(i)
            END DO
         ELSE !conserved variables at the rigid surface depend on liquid water
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tl_s_2d(i)= t_g(i) - lhocp*liqs(i,ke1)
               qt_s_2d(i)=qv_s(i) +       liqs(i,ke1)
            END DO
         END IF

         !$ACC LOOP SEQ
         DO k=ks, ke
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,k,u_m)=u(i,k)
               zvari(i,k,v_m)=v(i,k)
            END DO      
         END DO 

         IF (icldm_tran.EQ.-1) THEN !no water phase change possible
            !$ACC LOOP SEQ
            DO k=ks, ke
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  zvari(i,k,tet_l)=t(i,k)/epr(i,k)
                  zvari(i,k,h2o_g)=qv(i,k)
               END DO
            END DO
         ELSE !water phase changes are possible
            !$ACC LOOP SEQ
            DO k=ks, ke
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  zvari(i,k,tet_l)=(t(i,k) - lhocp*qc(i,k))/epr(i,k)
                  zvari(i,k,h2o_g)=qv(i,k) +       qc(i,k)
               END DO
            END DO
         END IF

         IF (lprfcor) THEN
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(len1, len2, lm, lh)
            DO i=ivstart, ivend
               len1=z2*h_top_2d(i)
               len2=(h_top_2d(i)-h_atm_2d(i))**2 &
                   /((hhl(i,ke-1)+hhl(i,ke))*z1d2-hhl(i,ke1)-h_atm_2d(i))
               lm=len1-tfm(i)*h_atm_2d(i)-len2
               lh=len1-tfh(i)*h_atm_2d(i)-len2

               zvari(i,ke,u_m  )=(len1*zvari(i,ke  ,u_m  ) &
                                 -len2*zvari(i,ke-1,u_m  ))/lm
               zvari(i,ke,v_m  )=(len1*zvari(i,ke  ,v_m  ) &
                                 -len2*zvari(i,ke-1,v_m  ))/lm
               zvari(i,ke,tet_l)=(len1*zvari(i,ke  ,tet_l)-h_atm_2d(i)*tfh(i)* t_g(i)/epr_2d(i) &
                                 -len2*zvari(i,ke-1,tet_l))/lh
               zvari(i,ke,h2o_g)=(len1*zvari(i,ke  ,h2o_g)-h_atm_2d(i)*tfh(i)*qv_s(i) &
                                 -len2*zvari(i,ke-1,h2o_g))/lh
            END DO
         END IF

!        Thermodynamische Hilfsvariablen auf dem Unterrand der Prandtl-Schicht:

         !$ACC END PARALLEL
         
         CALL adjust_satur_equil ( i1dim=nvec, khi=ke1, ktp=ke-1,      & !in
!
              i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1,            & !in
!
              lcalrho=.TRUE.,  lcalepr=.FALSE., lcaltdv=.TRUE.,        & !in
              lpotinp=.FALSE., ladjout=.FALSE.,                        & !in
!
              icldmod=icldm_tran,                                      & !in
!
              zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          & !in
!
              prs=prss, t=tmps, qv=vaps, qc=liqs,                      & !in (surface values at level 'ke1')
!
              psf=ps, fip=tfh,                                         & !in
!
              rcld=rcld,  & !inp: std. deviat. of local super-saturat.
                            !out: saturation fraction (cloud-cover)
!
              dens=rhon,         exner=zaux(:,:,1),                    & !out
              r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3),                    & !out
              g_tet=zaux(:,:,4), g_h2o=zaux(:,:,5),                    & !out                                   
!
              tet_liq=zvari(:,:,tet_l), q_h2o=zvari(:,:,h2o_g),        & !inout (inp as target of 'tmps, vaps')
                                        q_liq=zvari(:,:,liq),          & !out
!
              lacc=lzacc, opt_acc_async_queue=acc_async_queue )

         !Beachte:
         !'zvari(:,ke1,tet_l)' und 'zvari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen am Unterrand der 
         ! Prandtl-Schicht, waehrend  ta_2d' => 'zvari(:,ke,tet_l) und 'qda_2d  => 'zvari(:,ke,h2o_g)' auf diese 
         ! Groessen bzgl. der untersten Hauptflaeche zeigen, welche zur Interpolation der Groessen an der Oberflaeche 
         ! auf jenen Unterrand der Prandtl-Schicht (Oberrand der Rauhigkeitsschicht) benutzt werden.

         ! Berechnung der benoetigten Vertikalgradienten und der TKE-Antriebe:

         !Vertikalgradienten des Horizontalwindes:

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            val_m(i)=tfm(i)/dz_0a_m(i) !reciprocal resistance lenght of roughness- and laminar-layer
         END DO

         !$ACC LOOP SEQ
         DO n=1, nvel
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,ke1,n)=zvari(i,ke,n)*val_m(i) !vertical gradient of wind-component
            END DO
            !Beachte: Dies ist die Darstellung ohne Nutzung der unteren Randwerte der Prandtl-Schicht
         END DO

         !Scherungs-Antrieb der TKE:
         IF (itype_sher.EQ.2 .AND. PRESENT(dwdx) .AND. PRESENT(dwdy)) THEN
            !Einschliesslich der 3D-Korrektur durch den Vertikalwind bzgl. der mittleren Hangneigung:
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=MAX( (zvari(i,ke1,u_m)+dwdx(i,ke1)*val_m(i))**2 &
                              +(zvari(i,ke1,v_m)+dwdy(i,ke1)*val_m(i))**2 &
                              +hdef2(i,ke1)*val_m(i)**2, fc_min(i) )
            END DO
            !Beachte: 
            !'dwdx(ke1)', 'dwdy(ke1)' und 'hdef2(ke1)' beziehen sich auf die vorlaeufige Schichtdicke "1m".
            !Diese Felder sind in ICON fuer den Level 'ke1' nicht vorhanden!
         ELSE
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=MAX( zvari(i,ke1,u_m)**2+zvari(i,ke1,v_m)**2, fc_min(i) )
            END DO
         END IF    
         IF (ladsshr) THEN !treatment of additional shear by NTCs or LLDCs active
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=frc_2d(i)*frm(i,ke1)
            END DO
         END IF   

         ! Vertikalgradienten der dynamisch wirksamen Skalare:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            edh(i)=z1/dz_0a_h(i)   !reciprocal resistance lenght for momentum of turb. Prandt.-layer
            val_h(i)=tfh(i)*edh(i) !reciprocal resistance lenght of roughness- and laminar-layer
         END DO

         IF (imode_trancnf.EQ.1) THEN !old version of zero-level-gradients requested
            !Transformation of Tet_l-gradient into the old form following from interpolation
            !onto the zero-level in terms of T_l (rather than Tet_l) and correcting the
            !calculated T_l-Gradient by the adiabatic laps rate:

!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert, val1, val2)
            DO i=ivstart, ivend

               !Estimated zero-level values according to the old definition:
               wert=(epr(i,ke)-epr_2d(i))*zvari(i,ke,tet_l)                          !temp.-deviation due to interpol.
               val1=epr_2d(i)*zvari(i,ke1,tet_l)+wert*(z1-tfh(i))                    !liqu. water temp. 
               val2=zvari(i,ke1,h2o_g)-zvari(i,ke1,liq)                              !specfic humidity 
               val2=(val1+lhocp*zvari(i,ke1,liq))*(z1+rvd_m_o*val2-zvari(i,ke1,liq)) !virt. temperature

               !Tet_l-gradient according to the old definition:
               wert=(wert*val_h(i)+tet_g*val1/val2)/epr_2d(i)                        !deviation of Tet_l-grad.
               zvari(i,ke1,tet_l)=(zvari(i,ke,tet_l)-zvari(i,ke1,tet_l))*edh(i)+wert !old repr. of Tet_l-grad.

               !H2O_g-gradient:
               zvari(i,ke1,h2o_g)=(zvari(i,ke,h2o_g)-zvari(i,ke1,h2o_g))*edh(i)
            END DO
         ELSE
           !$ACC LOOP SEQ
           DO n=tet_l, h2o_g
              !$ACC LOOP GANG(STATIC: 1) VECTOR
              DO i=ivstart, ivend
                 zvari(i,ke1,n)=(zvari(i,ke,n)-zvari(i,ke1,n))*edh(i)
              END DO
           END DO
         END IF
         !'zvari(:,ke1,n)' enthaelt jetzt die Vertikalgradienten der Erhaltungsvariablen

         !Auftriebs-Antrieb der TKE:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            frh(i,ke1)=g_tet(i)*zvari(i,ke1,tet_l)+g_vap(i)*zvari(i,ke1,h2o_g)
         END DO
!US !$acc end parallel
!US !$acc parallel async(1) default(none) if(lzacc)

         ! Berechnung der Stabilitaetslaengen:

         IF (it_durch.EQ.it_start .AND. lini) THEN !Startinitialisierung

            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt, val1, val2, wert)
            DO i=ivstart, ivend
               IF (frh(i,ke1).GE.(z1-rim)*frm(i,ke1)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und 'lm', sowie 'lh'
                  ! werden durch 'lm' bei der krit. Ri-Zahl angenaehert:
                  fakt=z1/rim-z1
                  tkvm(i,ke1)=l_tur_z0(i)*(sm_0-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=tkvm(i,ke1)
               ELSE
                  fakt=frh(i,ke1)/(frm(i,ke1)-frh(i,ke1))
                  tkvm(i,ke1)=l_tur_z0(i)*(sm_0-(a_6+a_3)*fakt)
                  tkvh(i,ke1)=l_tur_z0(i)*(sh_0-a_5*fakt)
               END IF

               val1=tkvm(i,ke1)*frm(i,ke1)
               val2=tkvh(i,ke1)*frh(i,ke1)
               wert=MAX( val1-val2, rim*val1 )

               IF (.NOT.ltkeinp) THEN !TKE not present as input
                  tke(i,ke1,nvor)=MAX( SQRT(d_m*l_tur_z0(i)*wert), vel_min )
               END IF

!                 Retrieving this peace of out-commented code:
               !Note:
               !'tkvm|h' are stability-dependent length-scales here, which should not be dependent on LLDCs.
            END DO    

         ELSE ! mit Hilfe der vorhergehenden TKE-Werte

!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
            DO i=ivstart, ivend
               wert=z1/tke(i,ke1,nvor)
               tkvm(i,ke1)=tkvm(i,ke1)*wert
               tkvh(i,ke1)=tkvh(i,ke1)*wert
            END DO
         END IF

         !$ACC END PARALLEL

! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes, der Stabilitaetsfuntionen und des SDSS:

         CALL solve_turb_budgets( it_s=it_durch, it_start=it_start,                             &

                                  i1dim=nvec, i_st=ivstart, i_en=ivend,                         & !in 
 
                                  khi=ke1, ktp=ke-1, kcm=kcm, k_st=ke1, k_en=ke1, k_sf=ke1,     & !in

                                  ntur=ntur, nvor=nvor,                                         & !in

                                  lssintact=.FALSE.,      lupfrclim=(imode_trancnf.EQ.1),       & !in
                                  lpres_edr=PRESENT(edr),                                       & !in
                                  ltkeinp=ltkeinp,        lstfnct=lstfnct,                      & !in

                                  imode_stke=imode_tran,  imode_vel_min=imode_vel_min,          & !in

                                  dt_tke=dt_tke, fr_tke=fr_tke,                                 & !in

                                  fm2=frm, fh2=frh, ft2=frm,                                    & !in
                                  lsm=tkvm, lsh=tkvh, tls=len_scale,                            & !in(out)

                                  tvt=tvt, velmin=velmin,                                       & !in 
                                  tke=tke, ediss=ediss,                                         & !inout, out

                                  lactcnv=(icldm_tran.NE.-1 .AND. lsrflux),                     & !in (act. flux conversion)
                                  laddcnv=lsrflux,                                              & !in (add. flux-conversion)
                                  exner=zaux(:,:,1), r_cpd=zaux(:,:,2), qst_t=zaux(:,:,3),      & !in

                                  rcld=rcld, & !inp: effective saturation fraction (cloud-cover)
                                               !out: std. deviat. of local super-saturat.
                                               !     (only for last iteration step)      

                                  lcircterm=.FALSE.,                                            & !in
                                  dens=rhon, l_pat=l_pat, l_hori=l_hori,                        & !in

                                  grd=zvari,  & !inp: vert. grads. (incl. those of tet_l, h2o_g and liq)
                                                !out: vert. grads. (incl. those of tet,   vap   and liq
                                                !                   resulting from flux conversion)
                                  !'zvari'-output is only calculated at the last iteration step; and it's still equal 
                                  ! to the input, if these calculations are not executed.
                                  !As "lcircterm=F", 'prss=>zvari(:,:,0)' is still near-surface pressure.

                                  lacc=lzacc, opt_acc_async_queue=acc_async_queue               ) !in

         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(val1, val2, wert)
         DO i=ivstart, ivend

! 4h)       Bestimmung der durch Wirkung der L-Schicht korrigierten Diffusionskoeffizienten
!           und der zugehoerigen Transfer-Geschwindigkeiten:

!           Unkorrigierte Diffusionskoeffizienten:

            val1=con_m; tkvm(i,ke1)=tke(i,ke1,ntur)*tkvm(i,ke1)
            val2=con_h; tkvh(i,ke1)=tke(i,ke1,ntur)*tkvh(i,ke1)

            IF (imode_lamdiff.EQ.2) THEN !surface-layer diffusion coeff. always at least at laminar value
               IF (imode_tkemini.EQ.2) THEN !adaptation of TKE and TMod. to lower limits
                  tke(i,ke1,ntur)=tke(i,ke1,ntur)*MAX( z1, val2/tkvh(i,ke1) ) !adapted 'tke'
                  wert=tkvm(i,ke1)/tkvh(i,ke1) !turbulent Prandtl-number
                  tkvh(i,ke1)=MAX( val2, tkvh(i,ke1) ) !'tkvh' with lower limit
                  tkvm(i,ke1)=wert*tkvh(i,ke1)         !'tkvm' with adapted lower limit
               ELSE
                  tkvh(i,ke1)=MAX( val2, tkvh(i,ke1) ) !'tkvh' with lower limit
                  tkvm(i,ke1)=MAX( val1, tkvm(i,ke1) ) !'tkvm' with lower limit
               END IF   
            END IF

         END DO
         IF (imode_trancnf.GE.4 .OR. (imode_trancnf.GE.2 .AND. it_durch.LT.it_end)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
            DO i=ivstart, ivend
               wert=l_tur_z0(i)*SQRT(SQRT(frm(i,ke1))/tkvm(i,ke1)) !updated tkr=Ustar/(q*Sm)_0
               tkr(i)=MERGE( ditsmot*tkr(i) + (z1-ditsmot)*wert, wert, ditsmot.GT.z0 )
            END DO
         END IF

         IF (it_durch.LT.it_end) THEN !at least one additional iteration will take place

            IF (imode_trancnf.EQ.2 .OR. imode_trancnf.EQ.3) THEN
               !new version of initializing the profile-factors using Ustar,
               !but still epressing this factor in terms of "(q*Sx)_p/(q*Sx)_0":
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
               DO i=ivstart, ivend
                  fakt=h_top_2d(i)/z0m_2d(i) !(l_p-l_0)/l_0; l_0=akt*z0m
                  rat_m_2d(i)=z1+fakt*(z1-tkr(i))                                       !(q*Sm)_p/(q*Sm)_0
                  rat_h_2d(i)=z1+fakt*(z1-tkr(i)*(tkvm(i,ke1)*sh_0)/(tkvh(i,ke1)*sm_0)) !(q*Sh)_p/(q*Sh)_0
               END DO

            ELSEIF (imode_trancnf.GE.4) THEN
               !new version of initializing the profile-factors and already expressing
               !them in terms of "Ustar/(q*Sh)_0*(Sh(0)/Sm(0))":

               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  rat_m_2d(i)=tkr(i)                                       !Ustar/(q*Sm)_0
                  rat_h_2d(i)=tkr(i)*(tkvm(i,ke1)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(sh(0)/sm(0))
               END DO
            END IF   

         END IF

         !$ACC END PARALLEL

         ! This should happen outside the OpenACC paralle region
         IF ( it_durch.LT.it_end .AND. .NOT.ltkeinp) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF

!----------------------------------------------
      END DO !Iterationen
!----------------------------------------------

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

! 4i) Belegung der Felder fuer die Transfer-Geschwindigkeiten:

      IF (.NOT.lini .AND. ditsmot.GT.0) THEN
         !previous values of 'tvm|h' need to be saved
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            tcm(i)=tvm(i) !previous tvm
            tch(i)=tvh(i) !previous tvh
         END DO
      END IF
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ! Transfer-Velocities:
         tvm(i)=tkvm(i,ke1)*val_m(i) !to be used
         tvh(i)=tkvh(i,ke1)*val_h(i) !to be used
      END DO
      IF (.NOT.lini .AND. ditsmot.GT.0) THEN !smoothing of transfer velocity required   
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            tvm(i)=ditsmot*tcm(i) + (z1-ditsmot)*tvm(i) !smoothed new tvm
            tvh(i)=ditsmot*tch(i) + (z1-ditsmot)*tvh(i) !smoothed new tvh
         END DO
      END IF
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
      DO i=ivstart, ivend
         ! Transfer-Coefficients:
         fakt=z1/vel_2d(i,ke)
         tcm(i)=tvm(i)*fakt
         tch(i)=tvh(i)*fakt
      END DO

!-----------------------------------------------
!-----------------------------------------------

! 4j) Berechnung der Enthalpie- und Impulsflussdichten sowie der EDR am Unterrand:

      IF ((lsrflux.OR.lrunscm) .AND. (PRESENT(shfl_s).OR.PRESENT(qvfl_s))) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
         DO i=ivstart, ivend
            wert=rhon(i,ke1)*tkvh(i,ke1)

            IF (PRESENT(shfl_s)) shfl_s(i)=cp_d*wert*zvari(i,ke1,tet)*epr_2d(i)
            IF (PRESENT(qvfl_s)) qvfl_s(i)=wert*zvari(i,ke1,vap)
            !Note: 'shfl_s' and 'qvfl_s' are positive downward and 'shfl_s' belogns to the T-equation!
        END DO
      END IF
  
      IF ((lsrflux.OR.lrunscm) .AND. (PRESENT(umfl_s).OR.PRESENT(vmfl_s))) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
         DO i=ivstart, ivend
            wert=rhon(i,ke1)*tkvm(i,ke1)
            IF (PRESENT(umfl_s)) umfl_s(i)=wert*zvari(i,ke1,u_m)
            IF (PRESENT(vmfl_s)) vmfl_s(i)=wert*zvari(i,ke1,v_m)
            !Note: 'umfl_s' and 'vmfl_s' are positive downward!
         END DO
      END IF

      !$ACC END PARALLEL

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN
         IF (PRESENT(shfl_s)) THEN
            IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
               !measured SHF has to be used for forcing:
               shfl_s(imb)=SHF%mod(0)%val
            ELSEIF (lsurflu) THEN !SHF defined by explicit surface flux density
               SHF%mod(0)%val=shfl_s(imb)
               SHF%mod(0)%vst=MAX(i_upd, SHF%mod(0)%vst) !SHF is at least updated
            END IF
         END IF
         IF (PRESENT(qvfl_s)) THEN
            IF (LHF%mod(0)%vst.GT.i_cal .AND. LHF%mod(0)%ist.EQ.i_mod) THEN
               !measured LHF has to be used for forcing:
               qvfl_s(imb)=LHF%mod(0)%val / lh_v
            ELSEIF (lsurflu) THEN !LHF defined by explicit surface flux density
               LHF%mod(0)%val=qvfl_s(imb) * lh_v
               LHF%mod(0)%vst=MAX(i_upd, LHF%mod(0)%vst) !LHF is at least updated
            END IF
            !Note: LHF always is the latent heat flux connected with evaporation by definition,
            !      independent whether the surface is frozen or not!
         END IF
      END IF
#endif
!SCLM --------------------------------------------------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      IF (lnsfdia) THEN !diagnostics at near surface level required at this place

!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend

         !Einschraenkung von z0m_dia ueber Land:

         IF (fr_land(i) <= z1d2) THEN
            !Ueber See gibt es keinen synoptischen Garten
            z0d_2d(i)=z0m_2d(i)
         ELSE
            !Die Rauhigkeitslaenge einer SYNOP Station soll immer
            !kleiner als 10m bleiben:
            z0d_2d(i)=MIN( h_10m, z0m_dia )
         END IF

         !Festlegung der synoptischen Niveaus: 

         IF (itype_2m_diag.EQ.2) THEN !using an exponetial rougness layer profile
           z2m_2d (i) = h_2m -h_can_2d(i) !2m ueber dem Bodenniveau des Bestandes
           z10m_2d(i) = h_10m-z0m_2d(i)   !Hoehe, in der die turbulente Distanz 10m betraegt
         ELSE !using only a logarithmic profile above a SYNOP lawn
           z2m_2d (i) = h_2m
           z10m_2d(i) = h_10m
         END IF

         !Erste Belegung zweier benachbarter Modellniveaus:

         hk_2d(i)=h_atm_2d(i)
         hk1_2d(i)=z0
         k_2d(i)=ke

      END DO
      !$ACC END PARALLEL

!     Diagnose der 2m-Groessen:

      IF (ltst2ml) THEN !test required, whether 2m-level is above the lowest main-level
#ifdef _OPENACC
        CALL diag_level_gpu(ivstart, ivend, ke1, z2m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                            lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#else
        CALL diag_level(ivstart, ivend, ke1, z2m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                        lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#endif
      END IF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

      IF (itype_2m_diag.EQ.2) THEN !using an exponential rougness layer profile

         val2=z1/epsi
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(val1, fakt, wert, z_surf, fac_h)
         DO i=ivstart, ivend
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               IF (z2m_2d(i).LT.z0) THEN
!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

                  val1=z2m_2d(i)/dz_s0_h(i)
                  IF (-val1.LE.val2) THEN
                    fakt=dz_s0_h(i)/dz_sa_h(i)
                    fakt=MIN( z1, MAX( z0, fakt*EXP(val1) ) )
                  ELSE
                    fakt=z0
                  ENDIF
               ELSE
!                 2m-Niveau liegt innerhalb der Modell_Prandtl-Schicht
!                 mit logarithmischen Vertikalprofilen:

                  z_surf=z0m_2d(i)
                  fac_h=fac_h_2d(i)

                  IF (ABS(z1-fac_h) < epsi ) THEN
                     wert=z_surf*z2m_2d(i)/(z2m_2d(i)+z_surf)
                  ELSEIF (fac_h.GE.z0 .OR. imode_trancnf.LT.3) THEN
                     !non-stable strat. or using only linear interpolation of profile-function
                     !for the velocity scale (q*Sh):
                     wert=z_surf*LOG((z2m_2d(i)+z_surf)/(z_surf+fac_h*z2m_2d(i)))/(z1-fac_h)
                  ELSE !hyperbolic interpolation of profile-function (q*Sh) for stable stratification
                     wert=z2m_2d(i)/z_surf
                     IF (imode_trancnf.EQ.3) THEN !only if the upper node for diffusion coefficients is used
                        wert=z_surf*((z1-fac_h)*LOG(wert+z1)+fac_h*wert)/(z1-fac_h)
                     ELSE !(imode_trancnf.GE.4): without using the upper node for diff.-coefs.
                        wert=z_surf*(LOG(wert+z1)-fac_h*wert)/(z1-fac_h)
                     END IF
                  END IF
                  fakt=(dz_s0_h(i)+wert)/dz_sa_h(i)
               END IF

               IF (imode_syndiag.EQ.1) THEN !direkte interpol. von temperatur und spezifischer Feuchte
                  tmps(i,ke1) = t_g(i) + (t(i,ke)-t_g(i))*fakt &
                              + tet_g*( (h_atm_2d(i)+h_can_2d(i) )*fakt-h_2m ) !Achtung: mit 'tet_g'-Korrektur
                  vaps(i,ke1)= qv_s(i) + (qv(i,ke)-qv_s(i))*fakt
               ELSE !interpolation von Erhaltungsvariablen
                  tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*tl_s_2d(i)/epr_2d(i)
                  vaps(i,ke1) = qt_s_2d(i) + fakt*(qda_2d(i)-qt_s_2d(i))
               END IF 

            END IF
         END DO

      ELSE !using only a logarithmic profile above a SYNOP lawn

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(a_atm, a_2m, fr_sd_h, val1, val2, fakt, wert, z_surf, fac_h)
         DO i=ivstart, ivend
            IF (k_2d(i).EQ.ke .OR. .NOT.ltst2ml) THEN
               !2m-Niveau unterhalb der untersten Modell-Hauptflaeche
               !in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z_surf=z0d_2d(i)

               a_atm=h_atm_2d(i)+z_surf
               a_2m=h_2m+z_surf

               !Dimensionsloser Widerstand des Rauhigkeits-Schicht der SYNOP-Wiese (zwischen dem Unterrand
               !an der Erdoberflaeche und Oberrand bzgl. der Rauheigkeitslaenge 'z0d' einer SYNOP-Wiese
               !(wobei die turbulente Geschwindigkeitsskala und der Oberflaechenindex die fuer das gesamte
               !Gitterelemnt gueltigen Werte behalten):
               fr_sd_h=MAX( z0, dz_s0_h(i)/z0m_2d(i)+LOG(z_surf/z0m_2d(i)) )

               IF (imode_trancnf.LT.4) THEN !only if the  upper node for diffusion coefficients is used
                  fac_h_2d(i)=(rat_h_2d(i)-z1)*z_surf/h_top_2d(i) !re-defined profile-factor employing 'z0d'
                  !Attention:
                  !Possibly, the original profile-factor 'fac_h_2d' should be considered as a constant
                  ! of the vertical profile, rather than the re-defined one.
               END IF

               !Verhaeltnis der dimensionslosen Widerstaende:

               fac_h=fac_h_2d(i)

               IF (fac_h.GE.z0 .OR. imode_trancnf.LT.3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sh):
                  IF (fac_h.EQ.z1) THEN
                     val1=fr_sd_h+h_2m/a_2m
                     val2=fr_sd_h+h_atm_2d(i)/a_atm
                  ELSE
                     fakt=z1/(z1-fac_h)
                     val1=fr_sd_h+LOG(a_2m /(z_surf+fac_h*h_2m       ))*fakt
                     val2=fr_sd_h+LOG(a_atm/(z_surf+fac_h*h_atm_2d(i)))*fakt
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sh) for stable stratification
                  wert=z1/z_surf
                  IF (imode_trancnf.EQ.3) THEN !only if the upper node for diffusion coefficients is used
                     fac_h_2d(i)=-fac_h/rat_h_2d(i) !transformed profile-factor for stable strat.
                     val1=fr_sd_h+(z1-fac_h)*LOG(a_2m *wert)+fac_h*h_2m       *wert
                     val2=fr_sd_h+(z1-fac_h)*LOG(a_atm*wert)+fac_h*h_atm_2d(i)*wert
                  ELSE !(imode_trancnf.GE.4): without using the upper node for diff.-coefs.
                     fakt=z1/(z1-fac_h)
                     val1=fr_sd_h+(LOG(a_2m *wert)-fac_h*h_2m       *wert)*fakt
                     val2=fr_sd_h+(LOG(a_atm*wert)-fac_h*h_atm_2d(i)*wert)*fakt
                  END IF 
               END IF

               fakt=val1/val2

               !Interpolationswerte fuer das synoptische 2m-Niveau:

               IF (imode_syndiag.EQ.1) THEN
                  tmps(i,ke1) = t_g(i) + (t(i,ke)-t_g(i))*fakt &
                              + tet_g*(h_atm_2d(i)*fakt-h_2m)
                  vaps(i,ke1) = qv_s(i) + (qv(i,ke)-qv_s(i))*fakt
               ELSE
                  tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*tl_s_2d(i)/epr_2d(i)
                  vaps(i,ke1) = qt_s_2d(i) + fakt*(qda_2d(i)-qt_s_2d(i))
                  IF (icldm_tran.GT.-1) THEN !water phase change is possible
                     fakt=h_2m/h_atm_2d(i)
                     rcls(i,ke1)=rcld(i,ke1)+fakt*(rcld(i,ke)-rcld(i,ke1))
                  END IF
               END IF
            END IF
         END DO
      END IF !using only a logarithmic profile above a SYNOP lawn

      !$ACC END PARALLEL

      IF (ltst2ml) THEN !test required, whether 2m-level is above the lowest main-level
      !$ACC KERNELS ASYNC(acc_async_queue) DEFAULT(PRESENT) IF(lzacc)
      k=MINVAL(k_2d(ivstart:ivend))
      !$ACC END KERNELS

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      IF (k.LT.ke) THEN !2m-level is above the lowest main-level at least for one grid point
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(k2, k1, fakt, wert, val1, val2)
      DO i=ivstart, ivend
         IF (k_2d(i).LT.ke) THEN
!           2m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!           trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommenen
!           linearen Profile der progn. Modellvariablen eine logarith. Interpolation:

            k2=k_2d(i); k1=k2+1

            fakt=z1/(hk1_2d(i)+z0d_2d(i))
            wert=(h_2m    +z0d_2d(i))*fakt
            fakt=(hk_2d(i)+z0d_2d(i))*fakt
            fakt=LOG(wert)/LOG(fakt)

            IF (imode_syndiag.EQ.1) THEN
               tmps(i,ke1)= t(i,k1)+fakt*( t(i,k2)- t(i,k1))
               vaps(i,ke1)=qv(i,k1)+fakt*(qv(i,k2)-qv(i,k1))
            ELSEIF (icldm_tran.EQ.-1) THEN !no water phase change possible
               val2=qv(i,k2)          ; val1= qv(i,k1)        ; vaps(i,ke1)=val1+fakt*(val2-val1)
               val2= t(i,k2)/epr(i,k2); val1=t(i,k1)/epr(i,k1); tmps(i,ke1)=val1+fakt*(val2-val1)
            ELSE !water phase changes are possible
               val2=qv(i,k2)+qc(i,k2) ; val1=qv(i,k1)+qc(i,k1); vaps(i,ke1)=val1+fakt*(val2-val1)
               val2=(t(i,k2)-lhocp*qc(i,k2))/epr(i,k2); val1=(t(i,k1)-lhocp*qc(i,k1))/epr(i,k1)
                                                               tmps(i,ke1)=val1+fakt*(val2-val1)
               rcls(i,ke1)=rcld(i,k1)+fakt*(rcld(i,k2)-rcld(i,k1))
            END IF
         END IF
         !Note: In the case "k_2d(i).EQ.ke" 'tmps' and 'vaps' have already been
         !      calculated above.
      END DO
      END IF !At least for one grid-point, the 2m-level is above the lowest main-level.
      !$ACC END PARALLEL
      END IF !test required, whether 2m-level is above the lowest main-level

      !Druck im 2m-Niveau:
      
      IF (imode_2m_diag.GT.0) THEN !pressure correction of 2m-level desired
!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
      DO i=ivstart, ivend
!test<
         wert=tmps(i,ke1)*(z1+rvd_m_o*vaps(i,ke1)) !angenaeherte virt. Temp.
         prss(i,ke1)=prss(i,ke1)                &  !Druck
                    *EXP(-(z2m_2d(i)-hk1_2d(i))*grav/(r_d*wert))
!prss(i,ke1)=prss(i,ke1)-(z2m_2d(i)-hk1_2d(i))*grav*rhon(i,ke1)
!test>
      END DO
      !$ACC END PARALLEL
      END IF

      IF (imode_syndiag.EQ.1 .OR. icldm_tran.EQ.-1) THEN
!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
             t_2m(i)=tmps(i,ke1)
            qv_2m(i)=vaps(i,ke1)
         END DO
         !$ACC END PARALLEL
      ELSE
!        Berechnung der zugehoerigen Modell- und Feuchtevariablen im 2m-Niveau
!        aus den Erhalturngsvariablen.

         CALL adjust_satur_equil ( i1dim=nvec, khi=ke1, ktp=ke-1,              & !in
!
              i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1,                    & !in
!
              lcalrho=.FALSE., lcalepr=.TRUE., lcaltdv=.FALSE.,                & !in
              lpotinp=.TRUE. , ladjout=.TRUE.,                                 & !in
!
              icldmod=icldm_tran,                                              & !in
!
              zrcpv=tur_rcpv, zrcpl=tur_rcpl,                                  & !in
!
              prs=prss, t=tmps, qv=vaps,                                       & !in (pres. and conserved variabs. at 2m-level)
!
              psf=ps,                                                          & !in
!
              ! Note: The follwing REAL variables will be overwritten by values valid for the 2m-level:
!
              rcld=rcls,  & !inp: std. deviat. of local super-saturat. at 2m-level
                            !out: saturation fraction (cloud-cover)    at 2m-level (not yet used)
!
              exner=zaux(:,:,1),                                               & !out (not usd)
                                                                                 !aux (internally) 
!
              tet_liq=zvari(:,:,tet), q_h2o=zvari(:,:,vap),                    & !inp: as target of 'tmps, vaps'
                                                                                 !out: adjusted  variables at 2m-level
                                      q_liq=zvari(:,:,liq),                    & !out: cloud-water of 2m-fog (not yet used)
!
              lacc=lzacc, opt_acc_async_queue=acc_async_queue )

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
             t_2m(i)=zvari(i,ke1,tet)
            qv_2m(i)=zvari(i,ke1,vap)
         END DO
         !$ACC END PARALLEL
      END IF

      IF (lfreeslip) THEN ! only for idealized dry runs with free-slip condition
!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
            qv_2m(i)=z0
            rh_2m(i)=z0
            td_2m(i)=z0
            u_10m(i)=z0
            v_10m(i)=z0
         END DO
         !$ACC END PARALLEL

      ELSE !continuation of near-surface diagnostics required

!        Finale 2m-Diagnose:

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR PRIVATE(patm, fakt, wert)
         DO i=ivstart, ivend
            patm=prss(i,ke1)*qv_2m(i) &
                /(rdv+(z1-rdv)*qv_2m(i))          !Wasserdampfdruck

            fakt=patm/zpsat_w( t_2m(i) )
            rh_2m(i)=100.0_wp*MIN( fakt, z1 )     !relative Feuchte
   
            !UB: old formulation           
            wert=LOG(patm/b1)
            !UB: For dry atmosphere, the Teten's formula is not defined at vapor 
            !    pressure patm=0 because of log(patm/b1). However, it converges 
            !    to a dew point of the value of parameter b4w in case we impose 
            !    a very small positive lower bound on  patm:
            ! erst mal wieder weggenommen, da man den Absturz will:  
            !wert=LOG(MAX(patm,1.0E-16_wp)/b1)

            td_2m(i)=MIN( (b2w*b3-b4w*wert) &
                         /(b2w-wert), t_2m(i) )   !Taupunktstemperatur
         END DO
         !$ACC END PARALLEL

!        Diagnose der 10m-Groessen:

#ifdef _OPENACC
         CALL diag_level_gpu(ivstart, ivend, ke1, z10m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                             lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#else
         CALL diag_level(ivstart, ivend, ke1, z10m_2d, hhl, k_2d, hk_2d, hk1_2d, &
                         lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#endif


         !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(vel2_2d) ASYNC(acc_async_queue) IF(lzacc)

!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR PRIVATE(a_atm, a_10m, val1, val2, fakt, wert, z_surf, fac_m)
         DO i=ivstart, ivend

            IF (k_2d(i).EQ.ke) THEN

!              10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z_surf=z0d_2d(i)

               a_atm=h_atm_2d(i)+z_surf
               a_10m=h_10m+z_surf

               IF (imode_trancnf.LT.4) THEN !further corrected profile-factor using an upper node
                  fac_m_2d(i)=(rat_m_2d(i)-z1)*z_surf/h_top_2d(i)
                  !Attention:
                  !Possibly, the original profile-factor 'fac_m_2d' should be considered as a constant
                  ! of the vertical profile, rather than the re-defined one.
               END IF

               !Verhaeltnis der dimensionslosen Widerstaende:

               fac_m=fac_m_2d(i)

               IF (fac_m.GE.z0 .OR. imode_trancnf.LT.3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sm):
                  IF (fac_m.EQ.z1) THEN
                     val1=h_10m/a_10m
                     val2=h_atm_2d(i)/a_atm
                  ELSE
                     val1=LOG(a_10m/(z_surf+fac_m*h_10m))
                     val2=LOG(a_atm/(z_surf+fac_m*h_atm_2d(i)))
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sm) for stable stratification
                  wert=z1/z_surf
                  IF (imode_trancnf.EQ.3) THEN !only if the upper node for diffusion coefficients is used
                     fac_m_2d(i)=-fac_m/rat_m_2d(i) !transformed profile-factor for stable stratification
                     val1=(z1-fac_m)*LOG(a_10m*wert)+fac_m*h_10m      *wert                
                     val2=(z1-fac_m)*LOG(a_atm*wert)+fac_m*h_atm_2d(i)*wert
                  ELSE !(imode_trancnf.GE.4): without using the upper node for diff.-coefs.
                     val1=LOG(a_10m*wert)-fac_m*h_10m      *wert 
                     val2=LOG(a_atm*wert)-fac_m*h_atm_2d(i)*wert
                  END IF
               END IF
   
               fakt=val1/val2

               !Interpolationswerte fuer das synoptische 10m-Niveau:
   
               u_10m(i)=vel1_2d(i)*fakt; v_10m(i)=vel2_2d(i)*fakt
 
            END IF   
         END DO
  
!DIR$ IVDEP
         !$ACC LOOP GANG VECTOR PRIVATE(fakt, wert, k1, k2)
         DO i=ivstart, ivend
            IF (k_2d(i).LT.ke) THEN
!              10m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!              trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommen
!              linearen Profile der progn. Modellvariablen eine logarithm. Interpolation:

               IF (ltst10ml) THEN
                  k2=k_2d(i); k1=k2+1
               ELSE
                  k2=ke-1; k1=ke
               END IF

               fakt=z1/(hk1_2d(i)+z0d_2d(i))
               wert=(h_10m   +z0d_2d(i))*fakt
               fakt=(hk_2d(i)+z0d_2d(i))*fakt
               fakt=LOG(wert)/LOG(fakt)
               u_10m(i)=u(i,k1)+fakt*(u(i,k2)-u(i,k1))
               v_10m(i)=v(i,k1)+fakt*(v(i,k2)-v(i,k1))
            END IF
         END DO
 
         !$ACC END PARALLEL

      END IF !continuation of near-surface diagnostics required

      !Notes:
      !'u|v_10m' always belong to mass points!
      !In case of ".NOT.lnsfdia" this kind of diagnostics is done at another place.

      END IF !diagnostics at near surface level required at this place

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

      IF (lnswinamp) THEN !amplification of near-surface wind by NTCs or LLDCs active
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
         DO i=ivstart, ivend
            fakt=SQRT(frc_2d(i)) !related amplification factor of wind speed
            vel_2d(i,ke)=fakt*vel_2d(i,ke)
            IF (lnswindia) THEN !diagnostics of near-surface wind active in in this SUB
               u_10m(i)=fakt*u_10m(i); v_10m(i)=fakt*v_10m(i)
            END IF
            frc_2d(i)=fakt
         END DO
      END IF

      IF (.NOT.lgz0inp .OR. lini) THEN
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(velo, wert, fakt)
      DO i=ivstart, ivend

!        Diagnose von 'gz0' (fuer den naechsten Zeitschritt)
!        ueber Wasserflaechen mit der (angepassten) Charnock-Formel

         IF (fr_land(i) <= z1d2) THEN !a non-land surface
            !Note: This definition of a non-land surface is now in line with the ICON-definition 
            !       using "frlnd_thrhld=z1d2"

            IF ( l_sice(i) ) THEN !ice-covered surface
               gz0(i)=g_z0_ice
            ELSE !water-covered surface
               !Bug_emlation: Previous treatment, which may get numerically unstable:
               velo=tke(i,ke,nvor)*z1d2
               !Note:
               !This refers to linear interpolation of turbulent velocity at the lowest atmospheric boundary-level 
               ! onto the lowest main-level, with vanishing magnitude at zero-level
               !The length-scales contributing to 'tke(i,ke1,ntur)' are assumed to be too small as to generate
               ! surface waves that can consturct the roughness layer below. Hence, 'velo' is not affected by them.
               !However, through 'tke(i,ke,nvor)', Ustar may also develop without vertical shear of mean wind,
               ! (e.g. at the convective limit).
               wert=MAX( epsi, tvm(i)*SQRT(vel_2d(i,ke)**2+velo**2) ) !effective Ustar**2
               !Note: 'vel_2d(i,ke)' may include velocity amplification by NTCs or LLDCs.

               ! Charnock-parameter:
               fakt=MERGE( alpha0_char(   & !extended value at "imode_charpar>1" apart from lakes, as a function of vel_10m
                                        MERGE( vel_2d(i,ke), & !substitute vel_ke by vel_10m at "lini .AND. .NOT.lnswindia'
                                               SQRT( u_10m(i)**2+v_10m(i)**2 ), lini .AND. .NOT.lnswindia ) ), &
                           MERGE( alpha0, & !standard value at "imode_charpar=1"
                                  0.1_wp, & !enhanced value at "imode_charpar>1" over lakes (with non-equlibrium wave spectrum)
                                  imode_charpar.EQ.1 ),                                                      &
                           imode_charpar.GT.1 .AND. .NOT.l_lake(i) )

               wert=MAX( g_len_min, fakt*wert+g_alpha1_con_m/SQRT(wert) )
               gz0(i)=MERGE( ditsmot*gz0(i)+(z1-ditsmot)*wert, wert, ditsmot.GT.z0 )
            END IF
         END IF
      END DO
      ENDIF  !lgz0inp

      !$ACC END PARALLEL

      !$ACC END DATA ! from acc data present
      !$ACC END DATA ! from acc data create

END SUBROUTINE turbtran

!==============================================================================

!+ Module procedure 'diag_level' for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, ke1, zdia_2d, hhl, k_2d, hk_2d, hk1_2d, lacc, opt_acc_async_queue)
   INTEGER, INTENT(IN) :: &
!
      ke1, &
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=wp), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER, INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                  !to be used for near surface diagnostics
 
   REAL (KIND=wp), INTENT(IN) :: &
!  
     hhl(:,:)

   REAL (KIND=wp), INTENT(INOUT) :: &
!
     hk_2d(:), & !mid level height above ground belonging to 'k_2d'
     hk1_2d(:)    !mid level height above ground of the previous layer (below)

   LOGICAL, OPTIONAL, INTENT(IN) :: lacc ! If true, use openacc
   INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

   LOGICAL :: lzacc ! non-optional version of lacc
   INTEGER :: acc_async_queue

   INTEGER :: i

   LOGICAL :: lcheck

   CALL set_acc_host_or_device(lzacc, lacc)

   IF(PRESENT(opt_acc_async_queue)) THEN
       acc_async_queue = opt_acc_async_queue
   ELSE
       acc_async_queue = 1
   ENDIF

   lcheck=.TRUE. !check whether a diagnostic level is above the current layer

   !XL_ACCTMP:this could be implemented with an explcit K loop, may be faster on GPU (no atomic)

   DO WHILE (lcheck) !loop while previous layer had to be checked
      lcheck=.FALSE. !check next layer ony, if diagnostic level is at least once
                     !above the current layer

      !$ACC PARALLEL ASYNC(acc_async_queue) PRESENT(hhl, zdia_2d, k_2d, hk_2d, hk1_2d) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=i_st,i_en
         IF (hk_2d(i)<zdia_2d(i) .AND. k_2d(i)>1) THEN !diagnostic level is above current layer
            !$ACC ATOMIC WRITE
            lcheck=.TRUE. !for this point or any previous one, the diagnostic level is above the current layer
            !$ACC END ATOMIC
            k_2d(i)=k_2d(i)-1
            hk1_2d(i)=hk_2d(i)
            hk_2d(i)=(hhl(i,k_2d(i))+hhl(i,k_2d(i)+1))*z1d2-hhl(i,ke1)
          END IF
       END DO
       !$ACC END PARALLEL
   END DO

END SUBROUTINE diag_level

!==============================================================================

!+ Module procedure 'diag_level_gpu' for computing the upper level index
!+ used for near surface diganostics (GPU-version)

SUBROUTINE diag_level_gpu (i_st, i_en, ke1, zdia_2d, hhl, k_2d, hk_2d, hk1_2d, lacc, opt_acc_async_queue)

   INTEGER, INTENT(IN) :: &
!
      ke1, &
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=wp), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER, INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                  !to be used for near surface diagnostics

   REAL (KIND=wp), INTENT(IN) :: &
!
     hhl(:,:)
     
   REAL (KIND=wp), INTENT(INOUT) :: &
!
     hk_2d(:), & !mid level height above ground belonging to 'k_2d'
     hk1_2d(:)    !mid level height above ground of the previous layer (below)

   LOGICAL, OPTIONAL, INTENT(IN) :: lacc
   INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

   LOGICAL :: lzacc
   INTEGER :: acc_async_queue

   INTEGER :: i, k, ke1_war

   CALL set_acc_host_or_device(lzacc, lacc)

   IF(PRESENT(opt_acc_async_queue)) THEN
       acc_async_queue = opt_acc_async_queue
   ELSE
       acc_async_queue = 1
   ENDIF

   ! NV HPC 23.3 workaround
   ke1_war = ke1
   ! Need to keep data section separate for now
   !$ACC DATA PRESENT(i_en) IF(lzacc)
   !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
   DO i=i_st,i_en
      IF (i>i_en) CYCLE ! NVHPC compiler WAR
      !$ACC LOOP SEQ
      DO k=k_2d(i)-1, 0, -1   
         IF (hk_2d(i)<zdia_2d(i) .AND. k_2d(i)>1) THEN !diagnostic level is above current layer
            k_2d(i)=k 
            hk1_2d(i)=hk_2d(i)
            hk_2d(i)=(hhl(i,k_2d(i))+hhl(i,k_2d(i)+1))*z1d2-hhl(i,ke1_war)
         ELSE
            EXIT
         END IF
      END DO
   END DO
   !$ACC END DATA
   
END SUBROUTINE diag_level_gpu
   
!==============================================================================

END MODULE turb_transfer
