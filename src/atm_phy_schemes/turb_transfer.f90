!
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
                tf_salt,  & ! salt water freezing point

    rdv, con_m, con_h, grav      

USE mo_convect_tables, ONLY : &
!
! Parameters for auxilary parametrizations:
! ------------------------------------------
!
    b1       => c1es,     & ! variables for computing the saturation steam pressure
    b2w      => c3les,    & ! over water (w) and ice (e)
    b4w      => c4les       !               -- " --

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
    ditsmot,      & ! smoothing factor for direct time-step iteration
    epsi,         & ! relative limit of accuracy for comparison of numbers
    it_end,       & ! number of initialization iterations (>=0)

! Parameters describing physical properties of the lower boundary 
! of the atmosphere:
!---------------------------------------------------------------
!
    rlam_mom,     & ! scaling factor of the laminar boudary layer for momentum
    rlam_heat,    & ! scaling factor of the laminar boudary layer for heat

    rat_can,      & ! factor for the canopy height
    rat_sea,      & ! ratio of laminar scaling factors for heat over sea and land
    rat_lam,      & ! ratio of laminar scaling factors for vapour and heat
    rat_glac,     & ! ratio of laminar scaling factors for heat over glaciers

    z0m_dia,      & ! roughness length of a typical synoptic station

    alpha0,       & ! lower bound for Charnock-parameter
    alpha1,       & ! parameter scaling the molek. roughness of water waves

    z0_ice,       & ! roughness length of sea ice

    len_min,      & ! minimal turbulent length scale [m]
    tur_len,      & ! maximal turbulent length scale [m]
    vel_min,      & ! minimal velocity scale [m/s]

    akt,          & ! von Karman-constant
    d_h=>d_heat,  & ! factor for turbulent heat dissipation
    d_m=>d_mom,   & ! factor for turbulent momentum dissipation

    ! derived parameters calculated in 'turb_setup'
    tet_g, rim, b_m, b_h, sm_0, sh_0,   &
    d_1, d_2, d_3, d_4, d_5, d_6,       &
    a_3, a_5 ,a_6,                      &
    tur_rcpv, tur_rcpl,                 &

    ! used derived types
    modvar,       & !

! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

    lprfcor,      & ! using the profile values of the lowest main level instead of
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
                    ! 1: zero surface flux density
                    ! 2: zero surface value
    imode_rat_sea,& ! mode of scaling the laminar resistance for heat over sea (related to 'rat_sea')
                    ! 1: constant ratio compared to land surface
                    ! 2: with a correction for a strongly overheated SST
    imode_vel_min,& ! mode of calculating the minimal turbulent velocity scale (in the surface layer only)
                    ! 1: with a constant value
                    ! 2: with a stability dependent correction
    imode_charpar,& ! mode of estimating the Charnock-Parameter
                    ! 1: use a constant value 
                    ! 2: use a wind-dependent value with a constant lower bound
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
    imode_tkemini,& ! mode of fixing a lower limit of q=2TKE**2
                    ! 1: by using 'vel_min' only
                    ! 2: by adapting to minimal diffusion coefficients
    imode_lamdiff   ! mode of considering laminar diffusion at surface layer
                    ! 1: only when calculating the profile functions
                    ! 2: surface-layer diffusion coeff. always at least at laminar value

USE turb_data, ONLY : &

    ! numbers and indices

    ninv    ,     & ! Anzahl der invarianten skalren Groessen (Erhaltungsvariablen)
    nvel    ,     & ! Anzahl der Geschwindigkeitskomponenten
    nred    ,     & ! = ninv + nvel
    ndim    ,     & !
    nmvar   ,     & !
    u_m     ,     & ! zonale Geschw.komp. im Massenzentrum
    v_m     ,     & ! meridionale  ,,      ,,     ,,
    tet_l   ,     & ! feucht-potentielle Temperatur
    h2o_g   ,     & ! Gesamtwasseergehalt
    liq             ! Fluessigwasser  ,,

!-------------------------------------------------------------------------------
! Control parameters for the run
!-------------------------------------------------------------------------------

! ICON data have to be declared for these variables, which is done later on
!-------------------------------------------------------------------------------

! Switches controlling other physical parameterizations:
USE mo_lnd_nwp_config,       ONLY: lseaice, llake, lterra_urb, itype_kbmo
!   
USE turb_data,         ONLY:   &
    itype_diag_t2m    !

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

#ifdef ICON_USE_CUDA_GRAPH
  LOGICAL, PARAMETER :: using_cuda_graph = .TRUE.
#else
  LOGICAL, PARAMETER :: using_cuda_graph = .FALSE.
#endif

!===============================================================================

CONTAINS

!===============================================================================

SUBROUTINE turbtran (                                                         &
!
          iini, ltkeinp, lgz0inp, lstfnct, lsrflux, lnsfdia, lrunscm,         &
!
          dt_tke, nprv, ntur, ntim,                                           &
!
          nvec, ke, ke1, kcm, iblock, ivstart, ivend,                         &
!
          l_hori, hhl, fr_land, depth_lk, h_ice, rlamh_fac, sai, urb_isa, gz0,&
!
          t_g, qv_s, ps, u, v, w, t, qv, qc, prs, epr,                        &
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
!     Es werden die Transferkoeffizienten fuer den Austausch von Impuls,
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
!     Die resultierenden Transferkoeffizienten 'tc(h,m)' sind die Kehrwerte
!     des Gesamtwiderstandes von den festen Oberflaechen (Neviau 's') bis
!     zur untersten Modellhauptflaeche (Niveau 'a').
!     Die turbulenten Diffusionskoeffizienten 'tkv(h,m)' fuer den vertikalen
!     Index 'ke1', beziehen sich aber auf den Unterrand des atmosphaerischen
!     Modells (Niveau '0').
!     Mit Hilfe der Felder 'tf(mh)' werden noch Reduktionsfaktoren der
!     Transferkoeffizienten durch die Wirkung der L-Schicht uebertragen.
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

   lrunscm         !a Single Column run (default: FALSE)

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
    kcm,          & ! level index of the upper canopy bound
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
    l_hori,       & ! horizontal grid spacing (m)
    rlamh_fac,    & ! scaling factor for rlam_heat
!
! External parameter fields:
! ----------------------------
    fr_land,      & ! land portion of a grid point area                    ( 1 )
    depth_lk,     & ! lake depth                                           ( m )
    sai,          & ! surface area index                                   ( 1 )
    urb_isa         ! impervious surface area fraction of the urban canopy ( 1 )

! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
    h_ice           ! ice thickness                                 (  m  )

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
!
    ps,           & ! surface pressure                              ( pa  )
    qv_s,         & ! specific water vapor content on the surface   (kg/kg)
    t_g             ! weighted surface temperature                  (  k  )

!XL: why not only INTENT IN ?
!REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
! Atmospheric model variables:
! ---------------------------------
     u,           & ! zonal wind speed       (at mass positions)    ( m/s )
     v,           & ! meridional wind speed  (at mass positions)    ( m/s )
     t,           & ! temperature                                   (  k  )
     qv,          & ! specific water vapor content                  (kg/kg)
     qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
     prs            ! atmospheric pressure                          ( pa  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
     epr            ! exner pressure                                 (1)

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
     w              ! vertical wind speed (defined on half levels)  ( m/s )

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------
!
     gz0             ! roughness length * g of the vertically not
                     ! resolved canopy                               (m2/s2)
!Achtung: Der g-Faktor ist ueberfluessig!

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(OUT) :: &
!    turbulent (transfer) velocity scales at the surface
     tvm,          & ! for momentum                                  ( m/s)
     tvh,          & ! for heat and moisture                         ( m/s)

     !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
     !velocities  'tvm' and 'tvh' makes live much easier!!               

!    turbulent transfer factors for laminar- and roughness-layer transfer
     tfm,          & ! of momentum                                     --
     tfh,          & ! of scalars                                      --
     tfv             ! of water vapor compared to heat                 --

!    turbulent transfer coefficients at the surface
REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
     tcm             ! for momentum                                  ( -- )
REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(OUT) :: &
     tch             ! for scalars (heat and moisture)               ( -- )

REAL (KIND=wp), DIMENSION(:), TARGET, OPTIONAL, INTENT(INOUT) :: &

!    reference surface diffusion coefficient
!    (only if "imode_trancnf.GE.2"):
     tkr             ! l*Ustar                                       (m2/s)

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=wp), DIMENSION(nvec,ke1,ntim), TARGET, INTENT(INOUT) :: &
     tke             ! q:=SQRT(2*TKE); TKE='turbul. kin. energy'     ( m/s )
                     ! (defined on half levels)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
     tkvm,         & ! turbulent diffusion coefficient for momentum  (m2/s )
     tkvh            ! turbulent diffusion coefficient for heat      (m2/s )
                     ! (and other scalars)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
     rcld            ! standard deviation of the local oversaturation
                     ! (as input and output)
                     ! fractional cloud cover (in turbdiff)            --

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(IN) :: &
     hdef2,        & ! horizontal deformation square at half levels  ( 1/s2 )
     dwdx,         & ! zonal      derivative of vertical wind  ,,    ( 1/s )
     dwdy            ! meridional derivative of vertical wind  ,,    ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(INOUT) :: &
     tketens         ! diffusion tendency of q=SQRT(2*TKE)           ( m/s2)

REAL (KIND=wp), DIMENSION(:,:), TARGET, OPTIONAL, INTENT(OUT) :: &
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
    n, nvar,    & !Index fuer diverse Schleifen oder Variablen
!
    nvor,       & !laufende Zeittstufe des bisherigen TKE-Feldes (wird zwischen Iterationen zu 'ntur')
    it_durch,   & !Durchgangsindex der Iterationen
    it_start      !Startindex der Iterationen

REAL (KIND=wp) :: &
    fr_tke,           & ! z1/dt_tke
    wert, val1, val2, & ! Platzhalter fuer beliebige Zwischenergebnisse
    fakt,             & !  ,,         ,,     ,,      Faktoren

!   Platzh. fuer therm. und mech. Antrieb der Turbulenz in (1/s)**2 (fh2,fm2):
    fh2,fm2, &

!   Platzh. fuer horiz. Geschw.-Komponenten, bel. Geschw. und Druck:
    vel1,vel2,velo, patm, &

!   Platzh. fuer Hoehendifferenzen und Laengaenskalen:
    dh,l_turb,lh,lm,z0d,z_surf,len1,len2, &
    dz_s0_m, dz_sa_m, &
    h_2m, h_10m, & !level heights (equal 2m and 10m)
    a_2m, a_10m, & !turbulent distance of 2m- and 10m-level (with respect to diag. roughness)
    a_atm, &       !turbulent distance of the atmosp. level (with respect to mean  roughness)

!   for TERRA_URB
    zkbmo_dia, zkbmo_urb, zustar, & ! inverse Stanton number
!   zdz_s0_h , &

!   sonstiges:
    rin_m,rin_h, fr_sd_h, xf

LOGICAL            ::  &
    lini,      & !initialization required
    lgz0ini      ! initialization of roughness lenght over water and ice

! Local arrays:

     REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!
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
        rcl_2d   (:),      &
!
        l_tur_z0 (:),      &
!
        tvt      (:,:)

REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!
!    pointer for eddy dissipation rate:
     ediss(:,:)

! these fields are still taken as local arrays, because the CRAY compiler cannot
! do the same optimizations with OMP threadprivate variables

REAL (KIND=wp), TARGET    :: &
  ! targets of used pointers
  diss_tar    (nvec,ke1),     & ! eddy dissipation rate (m2/s3)
  tketens_tar (nvec,ke1:ke1), & ! turbulent transport of SQRT(TKE) (m/s2)

  ! internal atmospheric variables
  len_scale   (nvec,ke1:ke1), & ! turbulent length-scale (m)
  l_scal      (nvec),         & ! reduced maximal turbulent length scale due to 
                                !   horizontal grid spacing (m)

  fc_min      (nvec),         & ! minimal value for TKE-forcing (1/s2)

  rclc        (nvec,ke1:ke1), & ! cloud cover
  rhon        (nvec,ke1:ke1), & ! boundary level air density (including surface level) (Kg/m3)
  frh         (nvec,ke1:ke1), & ! thermal forcing (1/s2) or thermal acceleration (m/s2)
  frm         (nvec,ke1:ke1), & ! mechan. forcing (1/s2) or mechan. accelaration (m/s2)

  zaux        (nvec,ke1,ndim),& ! auxilary array containing thermodynamical properties
                                ! (dQs/dT,ex_fakt,cp_fakt,g_tet,g_vap) or various auxiliary
                                ! variables for calculation of implicit vertical diffusion
  zvari       (nvec,ke1,ndim),& ! set of variables used in the turbulent 2-nd order equations
                                ! and later their effective vertical gradients

  prss        (nvec,ke1:ke1), & ! surface pressure (Pa)
  eprs        (nvec,ke1:ke1), & ! surface Exner-factor
  tmps        (nvec,ke1:ke1), & ! surface temperature (K)
  vaps        (nvec,ke1:ke1), & ! surface specific humidity
  liqs        (nvec,ke1:ke1), & ! liquid water content at the surface

  tl_s_2d     (nvec),         & ! surface level value of conserved temperature 
                                !   (liquid water temperature) (K)
  qt_s_2d     (nvec),         & ! surface level value of conserved humidity    
                                !   (total water)
  vel_2d      (nvec),         & ! wind speed (m/s) at the top of the transfer layer 
                                !   (lowest mid level of atm. model)

  velmin      (nvec),         & ! modified 'vel_min' used for tuning corrections 
                                !   (hyper-parameterizations) (m/s)
  ratsea      (nvec),         & ! modified 'rat_sea' used for tuning corrections 
                                !   (hyper-parameterizations)

  ! internal variables for the resistance model
  hk_2d       (nvec),         & ! mid level height above ground belonging to 'k_2d' (m)
  hk1_2d      (nvec),         & ! mid level height above ground of the previous layer
                                !   (below) (m)
  h_top_2d    (nvec),         & ! boundary level height of transfer layer (top  level)
  h_atm_2d    (nvec),         & ! mid      level heigth of transfer layer (atm. level)
  h_can_2d    (nvec),         & ! effective canopy height (m) used for calculation of 
                                !   roughness layer resistance for momentum
  edh         (nvec),         & ! reciprocal of a layer depth

  z0m_2d      (nvec),         & ! mean  roughness length
  z0d_2d      (nvec),         & ! diag. roughness length
  z2m_2d      (nvec),         & ! height of 2m  level (above the surface)
  z10m_2d     (nvec),         & ! height of 10m level (above the surface)

  rat_m_2d    (nvec),         & ! any surface layer ratio for momentum 
                                !   (like Re-number or Stability factor)
  rat_h_2d    (nvec),         & ! any surface layer ratio for scalars  
                                !   (like Re-number or Stability factor)
  fac_h_2d    (nvec),         & ! surface layer profile factor for scalars
  fac_m_2d    (nvec),         & ! surface layer profile factor for momentum

  frc_2d      (nvec),         & ! length scale fraction

  dz_sg_m     (nvec),         & ! laminar resistance lenght for momentum (m)
  dz_sg_h     (nvec),         & ! laminar resistance length for scalars (m)
  dz_g0_h     (nvec),         & ! turbulent roughness layer resistance length for scalars (m)
  dz_0a_m     (nvec),         & ! turbulent Prandtl-layer resistance length for momentum (m)
  dz_0a_h     (nvec),         & ! turbulent Prandtl-layer resistance length for scalars (m)
  dz_sa_h     (nvec),         & ! total transfer resistance length for scalars (m)
  dz_s0_h     (nvec)            ! total roughness layer resistance lenght for scalars (m)

REAL (KIND=wp) ::             &
  grad        (nvec,nmvar)           ! any vertical gradient

INTEGER        ::             &
  k_2d        (nvec)            ! index field of the upper level index to be used
                                !   for near surface diagn.
LOGICAL        ::             &
  lo_ice      (nvec)            ! logical sea ice indicator

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

  !GPU data region of all variables except pointers which are set later on
  !$ACC DATA CREATE(diss_tar, tketens_tar, len_scale, l_scal) &
  !$ACC   CREATE(fc_min, rclc, rhon, frh, frm, zaux, zvari, prss) &
  !$ACC   CREATE(eprs, tmps, vaps, liqs, tl_s_2d, qt_s_2d, vel_2d) &
  !$ACC   CREATE(velmin, ratsea, hk_2d, hk1_2d, h_top_2d) &
  !$ACC   CREATE(h_atm_2d, h_can_2d, edh, z0m_2d, z0d_2d, z2m_2d) &
  !$ACC   CREATE(z10m_2d, rat_m_2d, rat_h_2d, fac_h_2d, fac_m_2d) &
  !$ACC   CREATE(frc_2d, dz_sg_m, dz_sg_h, dz_g0_h) &
  !$ACC   CREATE(dz_0a_m, dz_0a_h, dz_sa_h, dz_s0_h, grad) &
  !$ACC   CREATE(k_2d, lo_ice) &
  !$ACC   COPYIN(ivend) &
  !$ACC   ASYNC(acc_async_queue) IF(lzacc)
!-------------------------------------------------------------------------------
  CALL turb_setup (ivstart=ivstart, ivend=ivend, ke1=ke1, &
                   iini=iini, dt_tke=dt_tke, nprv=nprv, l_hori=l_hori, qc_a=qc(:,ke), &
                   lini=lini, it_start=it_start, nvor=nvor, fr_tke=fr_tke,  &
                   l_scal=l_scal, fc_min=fc_min, liqs=liqs(:,ke1), rcld=rcld, tfm=tfm, tfh=tfh, &
                   lacc=lzacc, opt_acc_async_queue=acc_async_queue)

#ifdef ICON_USE_CUDA_GRAPH
   IF (lzacc .AND. lini) THEN
      CALL finish ('turbtran', 'initialization is not supported when capturing a graph with OpenACC')
   END IF
#endif
!-------------------------------------------------------------------------------

my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
my_thrd_id = omp_get_thread_num()
#endif

! Just do some check printouts:
  IF (ldebug) THEN
    DO i = ivstart, ivend
      IF (i ==  1 .AND. iblock == 1 .AND. my_cart_id == 0) THEN
        WRITE(*,'(A       )') '             '
        WRITE(*,'(A,2I5   )') 'TURB-DIAGNOSIS transfer:  iblock = ', iblock, i

        WRITE(*,'(A       )') ' Control Switches and Variables'
        WRITE(*,'(A,I28   )') '   iini             :  ', iini
        WRITE(*,'(A,L28   )') '   ltkeinp          :  ', ltkeinp
        WRITE(*,'(A,L28   )') '   lgz0inp          :  ', lgz0inp
        WRITE(*,'(A,L28   )') '   lstfnct          :  ', lstfnct
        WRITE(*,'(A,L28   )') '   lsrflux          :  ', lsrflux
        WRITE(*,'(A,L28   )') '   lnsfdia          :  ', lnsfdia
        WRITE(*,'(A,F28.16)') '   dt_tke           :  ', dt_tke
        WRITE(*,'(A,I28   )') '   nprv             :  ', nprv
        WRITE(*,'(A,I28   )') '   ntur             :  ', ntur
        WRITE(*,'(A,I28   )') '   ntim             :  ', ntim
        WRITE(*,'(A,I28   )') '   nvec             :  ', nvec
        WRITE(*,'(A,I28   )') '   kcm              :  ', kcm
        WRITE(*,'(A       )') ' Other input parameters:'
        WRITE(*,'(A,F28.16)') '   l_hori           :  ', l_hori      (i)
        WRITE(*,'(A,F28.16)') '   hhl    ke        :  ', hhl         (i,ke)
        WRITE(*,'(A,F28.16)') '   hhl    ke1       :  ', hhl         (i,ke1)
        WRITE(*,'(A,F28.16)') '   fr_land          :  ', fr_land     (i)
        WRITE(*,'(A,F28.16)') '   depth_lk         :  ', depth_lk    (i)
        WRITE(*,'(A,F28.16)') '   h_ice            :  ', h_ice       (i)
        WRITE(*,'(A,F28.16)') '   sai              :  ', sai         (i)
        WRITE(*,'(A,F28.16)') '   urb_isa          :  ', urb_isa     (i)
        WRITE(*,'(A,F28.16)') '   gz0              :  ', gz0         (i)
        WRITE(*,'(A,F28.16)') '   t_g              :  ', t_g         (i)
        WRITE(*,'(A,F28.16)') '   qv_s             :  ', qv_s        (i)
        WRITE(*,'(A,F28.16)') '   ps               :  ', ps          (i)
        WRITE(*,'(A,F28.16)') '   u     ke         :  ', u           (i,ke)
        WRITE(*,'(A,F28.16)') '   v     ke         :  ', v           (i,ke)
        WRITE(*,'(A,F28.16)') '   t     ke         :  ', t           (i,ke)
        WRITE(*,'(A,F28.16)') '   qv    ke         :  ', qv          (i,ke)
        WRITE(*,'(A,F28.16)') '   qc    ke         :  ', qc          (i,ke)
        WRITE(*,'(A,F28.16)') '   prs   ke         :  ', prs         (i,ke)
        WRITE(*,'(A,F28.16)') '   l_scal           :  ', l_scal      (i)
        WRITE(*,'(A,F28.16)') '   fc_min           :  ', fc_min      (i)
        WRITE(*,'(A,F28.16)') '   liqs             :  ', liqs        (i,ke1)
        WRITE(*,'(A,F28.16)') '   tke   ke         :  ', tke         (i,ke,ntur)
        WRITE(*,'(A,F28.16)') '   tke   ke1        :  ', tke         (i,ke1,ntur)
        WRITE(*,'(A,F28.16)') '   tkvm  ke         :  ', tkvm        (i,ke)
        WRITE(*,'(A,F28.16)') '   tkvm  ke1        :  ', tkvm        (i,ke1)
        WRITE(*,'(A,F28.16)') '   tkvh  ke         :  ', tkvh        (i,ke)
        WRITE(*,'(A,F28.16)') '   tkvh  ke1        :  ', tkvh        (i,ke1)
        WRITE(*,'(A       )') '             '
      ENDIF
    ENDDO
  ENDIF


  ! take care that all pointers have a target
  IF (PRESENT(edr)) THEN
     ediss => edr
  ELSE
     ediss => diss_tar
  END IF

! Unterste Hauptflaeche halbiert den Abstand zur
! untersten Nebenflaeche:

      xf=z2
      xf=z1/xf

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

      vel1_2d => zvari(:,ke ,u_m)
      vel2_2d => zvari(:,ke ,v_m)
      ta_2d   => zvari(:,ke ,tet_l)
      qda_2d  => zvari(:,ke ,h2o_g)

      qsat_dT => zaux(:,ke1,3)
      g_tet   => zaux(:,ke1,4)
      g_vap   => zaux(:,ke1,5)

      epr_2d  => eprs(:,ke1)
      rcl_2d  => rclc(:,ke1)

      l_tur_z0 => len_scale(:,ke1)

! 2)  Initialisierung der z0-Werte ueber Meer
!     und der laminaren Transferfaktoren:

      ! Set the logical mask lo_ice to distinguish between ice covered
      ! and open water sea or lake grid points.
      
      !$ACC DATA PRESENT(fr_land, t_g, h_ice, depth_lk) &
      !$ACC   PRESENT(hhl, epr_2d, u, v, t) &
      !$ACC   PRESENT(epr, tke, tkvm, tkvh, gz0, tkr) &
      !$ACC   PRESENT(l_tur_z0, ps, sai, urb_isa, rlamh_fac) &
      !$ACC   PRESENT(tfm, tfh, tfv, qv_s, qv, qc) &
      !$ACC   PRESENT(dwdx, dwdy, hdef2, g_tet) &
      !$ACC   PRESENT(g_vap, tcm, tch, z0d_2d, shfl_s) &
      !$ACC   PRESENT(rcld, qsat_dT, qvfl_s, umfl_s, vmfl_s) &
      !$ACC   PRESENT(edr, tvh, tvm, qsat_dt, t_2m) &
      !$ACC   PRESENT(qv_2m, rcl_2d, prs, qda_2d, ta_2d) &
      !$ACC   PRESENT(v_10m, u_10m, vel1_2d, vel2_2d, rh_2m, td_2m) &
      !$ACC   ASYNC(acc_async_queue) IF(lzacc)
      
! Berechnung einiger Hilfsgroessen und Initialisierung der Diffusionskoeff.:

      ! Hoehe des 2m- und 10m-Niveaus:
      h_2m  = z2
      h_10m = z10


      IF (lprfcor) THEN
         ks=ke-1
      ELSE
         ks=ke
      END IF

!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend

         IF (fr_land(i) < z1d2) THEN
            ! Water point.
            IF (.NOT. lseaice) THEN
              ! Sea ice model is not used.
              ! Ice surface if SST is less than the salt water freezing temperature.
              lo_ice(i) = t_g(i) < tf_salt
            ELSE
              ! Sea ice model is used.
              ! Ice surface if ice is present.
              lo_ice(i) = h_ice(i) > z0
            END IF
            IF (llake) THEN
              ! Lake model is used.
              ! Ice surface if this is a lake point AND ice is present.
              IF ((depth_lk(i) > z0) .AND. (h_ice(i) >= h_Ice_min_flk)) &
              lo_ice(i) = .TRUE.
            END IF
         END IF

      END DO

!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ! Dicke der Modell-Prandtl-Schicht
         h_top_2d(i) = hhl(i,ke)-hhl(i,ke1)
         h_atm_2d(i) = h_top_2d(i)*xf

         ! Surface-Exner-pressure:
         epr_2d(i) = zexner(ps(i))

         ! Scaling velocity = Wind at lowest full level:
         vel_2d(i) = MAX( vel_min, SQRT(u(i,ke)**2+v(i,ke)**2) )
      END DO

!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         ratsea(i) = rat_sea
      END DO

      IF (imode_rat_sea.EQ.2) THEN
!<Tuning
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            IF (t_g(i) - t(i,ke) > 8._wp) THEN
               ! Increase rat_sea for very large temperature differences between water and adjacent air
               ! in order to reduce excessive peaks in latent heat flux

               ratsea(i) = rat_sea*(1._wp + 0.05_wp*(t_g(i) - t(i,ke) - 8._wp))
            END IF
         END DO
!>Tuning: This kind of correction can be substituded by a less ad-hoc approach.
      END IF

      IF (imode_vel_min.EQ.2) THEN
!<Tuning    
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
!>Tuning: his kind of correction can be substituded by a less ad-hoc approach.
      END IF

      IF (lini) THEN !only for initialization
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(lgz0ini, l_turb, dh, vel1, vel2, fm2, fh2, fakt, lm, lh, wert, val1, val2)
         DO i=ivstart, ivend

            lgz0ini=(.NOT.lgz0inp .AND. fr_land(i) < z1d2)

            IF (imode_trancnf.GE.2 .OR. (lgz0ini .AND. .NOT.lo_ice(i))) THEN
!              Einfachste Schaetzung der Schubspannung als Impusls-
!              flussdichte durch die Nebenflaeche ke mit Hilfe
!              einer diagnostischen TKE ohne Beruecksichtigung von
!              Feuchte-Effekten und mit neuchtralen Stabilitaets-
!              funktion:

               l_turb=h_top_2d(i) !approx. turb. length scale at level ke

!test: different turb. length scale
               l_turb=akt*MAX( len_min, l_turb/( z1+l_turb/l_scal(i) ) )
!  l_turb=akt*MIN( l_scal(i), hhl(i,ke)-hhl(i,ke1) )
!test

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

               IF (imode_tkemini.EQ.2) THEN
                  tke(i,ke,1)=tke(i,ke,1)*MAX( z1, val1/tkvm(i,ke), & !adapted tke
                                                   val2/tkvh(i,ke) )
               END IF
               tkvm(i,ke)=MAX(val1, tkvm(i,ke)) !corrected tkv
               tkvh(i,ke)=MAX(val2, tkvh(i,ke))
   
!Achtung: Die 'epsi'-Beschraenkung ist recht willkuerlich und fehlt in COSMO-Version!
               val2=MAX( epsi, tkvm(i,ke)*SQRT(fm2) ) !estimate of Ustar**2
               val1=SQRT(val2) !Ustar
            ELSE ! needed for vectorization
               l_turb = 0._wp
               val1   = 0._wp
               val2   = 0._wp
            END IF   

            IF (lgz0ini) THEN 
              ! Iinitialization of roughness length for water or ice covered surface:

              ! Use ice surface roughness or open-water surface roughness
              ! according to lo_ice:

               IF ( lo_ice(i) ) THEN ! ice covered surface
                  gz0(i)=grav*z0_ice
               ELSE !water covered surface
!                 Schubspannung abhaengiger Wellenhoehe:
                  IF (imode_charpar.EQ.1) THEN !constant Charnock-Parameter
                     fakt=alpha0
!US from ICON version 044780ed>
                  ELSE IF (depth_lk(i) > z0) THEN
                     ! enhanced Charnock parameter over lakes, parameterizing a non-equlibrium wave spectrum
                     fakt = 0.1_wp
!US<
                  ELSE
                     fakt=alpha0_char(vel_2d(i))
                     !Note: The argument of alpha0_char should be 'vel_10m', which might not yet be
                     !      present during initialization!
                  END IF
                  gz0(i)=MAX( grav*len_min, fakt*val2+alpha1*grav*con_m/val1 )
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
         z0m_2d(i)  = gz0(i)/grav   !mean roughness length
         l_tur_z0(i)= akt*z0m_2d(i) !turbulent length scale
         frc_2d(i)  = z0m_2d(i)/(h_top_2d(i)+z0m_2d(i)) !length scale fraction
      END DO

      IF (lini) THEN

         IF (imode_trancnf.GE.2) THEN !new version of init. using estimated Ustar
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tkvm(i,ke1)=tkvm(i,ke)*frc_2d(i)*(frc_2d(i)+(z1-frc_2d(i))*rat_m_2d(i))
               tkvh(i,ke1)=tkvh(i,ke)*frc_2d(i)*(frc_2d(i)+(z1-frc_2d(i))*rat_h_2d(i))
            END DO
         ELSE    
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               tkvm(i,ke) =con_m; tkvh(i,ke) =con_h
               tkvm(i,ke1)=con_m; tkvh(i,ke1)=con_h
            END DO
         END IF    

      ELSEIF (imode_trancnf.GE.4) THEN
         !Not for initialization, but for calculation the profile-factors
         !without an upper node based on previous values of transfer-velocity
         !and diffusion-coefficients (at the top of the roughness layer):
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            rat_m_2d(i)= tkr(i)/tkvm(i,ke1)              !Ustar/(q*Sm)_0
            rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(Sh(0)/Sm(0))
            tkr(i)=rat_m_2d(i) !saved profile factor for momentum
                               !to be optionally smoothed during direct iteration
         END DO
      END IF

      IF (imode_trancnf.EQ.2 .OR. imode_trancnf.EQ.3) THEN
         !Profile-factors by using the previous diffusion coefficients
         !without a laminar correction, but still based on the upper node:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            rat_m_2d(i)=frc_2d(i)*tkvm(i,ke)/tkvm(i,ke1) !(q*Sm)_p/(q*Sm)_0
            rat_h_2d(i)=frc_2d(i)*tkvh(i,ke)/tkvh(i,ke1) !(q*Sh)_p/(q*Sh)_0
         END DO 
      END IF   

      !$ACC END PARALLEL

! 4)  Berechnung der Transferkoeffizienten:

      DO it_durch=it_start, it_end !Iterationen
!print *,"it_durch=",it_durch

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z_surf, fakt, rin_m, rin_h, zkbmo_dia, zkbmo_urb, zustar)
         DO i=ivstart, ivend

            z_surf= z0m_2d(i)/sai(i) !effektive Rauhigkeitslaenge

!           Laminare Korrektur der Diffusionskoeffizienten:
            tkvm(i,ke1)=MAX( con_m, tkvm(i,ke1) )
            tkvh(i,ke1)=MAX( con_h, tkvh(i,ke1) )

            fakt=z1+(z1-REAL(NINT(fr_land(i)),wp))*(ratsea(i)-z1)

            rin_m=tkvm(i,ke1)/con_m
            rin_h=tkvh(i,ke1)/con_h
!rin_m=tkvm(i,ke1)/con_m+z1
!rin_m=tkvh(i,ke1)/con_h+z1

!           Effektiven Widerstandslaengen der Rauhigkeits-Schicht:

            dz_sg_m(i)=rlam_mom*z_surf
            dz_sg_h(i)=fakt*rlam_heat*rlamh_fac(i)*z_surf*(rin_h/rin_m) * &
              MERGE(rat_glac, 1._wp, gz0(i)<0.01_wp .AND. fr_land(i)>=0.5_wp)  ! enhanced by a factor of 'rat_glac' over glaciers

          ! ohne lam. Grenzschicht fuer Skalare:
            dz_g0_h(i)=z_surf*LOG(rin_m)

          ! inclusive lam. Grenzschicht fuer Skalare:
            dz_s0_h(i)=dz_sg_h(i)+dz_g0_h(i)

            IF (lterra_urb .AND. (.NOT. itype_kbmo == 1)) THEN

              IF (urb_isa(i) > 0.0_wp) THEN
                zkbmo_dia    = dz_s0_h(i)/z0m_2d(i)
                zustar       = SQRT(tcm(i))*vel_2d(i)

                IF (itype_kbmo == 2) THEN
                  ! Brutsaert Kanda parameterisation for bluff-body elements
                  ! for some reason, it doesn't withstand zero values of ustar
                  zkbmo_urb = MAX(0.1_wp, 1.29_wp * (z0m_2d(i)*zustar/con_m)**0.25_wp - z2)
                ELSEIF (itype_kbmo == 3) THEN
                  ! Zilitinkevich
                  zkbmo_urb = MAX(0.1_wp, 0.13_wp * (z0m_2d(i)*zustar/con_m)**0.45_wp)
                END IF

                dz_s0_h(i) = (urb_isa(i)*zkbmo_urb + (1.0_wp - urb_isa(i))*zkbmo_dia) * z0m_2d(i)

              ENDIF

            ENDIF

         END DO    


!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
!           Effektiven Bestandeshoehe:
            IF (dz_sg_h(i).EQ.z0) THEN
               h_can_2d(i)=rat_can*z0m_2d(i)
            ELSE
               h_can_2d(i)=rat_can*dz_s0_h(i)*LOG(dz_s0_h(i)/dz_sg_h(i))
            END IF
         END DO

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert, dz_s0_m, fakt)
         DO i=ivstart, ivend

          ! inclusive lam. Grenzschicht fuer Impuls:
            wert=z1d2*dz_sg_m(i)
            dz_s0_m=wert+SQRT(wert**2+h_can_2d(i)*dz_sg_m(i))

          ! ohne lam. Grenzschicht fuer Impuls:
          ! this variable is never used later on
          ! dz_g0_m(i)=dz_s0_m-dz_sg_m(i)

!------------------------------------------------------------------------------

!           Profilfakoren der turbulenten Prandtl-Schicht:

            fakt=z0m_2d(i)/h_top_2d(i)

            IF (imode_trancnf.LT.4) THEN
               !Profile-factors by employing previous values of the diffusion-coefficients
               !at the top fo the roughness-layer (0) and also at the upper bound of the
               !lowest atm. model layer (p) as an upper node:


               IF (imode_trancnf.EQ.1) THEN !first version
                  !Profile factors by using the previous diffusion coefficients
                  !including a laminar correction:
                  rat_m_2d(i)=frc_2d(i)*tkvm(i,ke)/tkvm(i,ke1)
                  rat_h_2d(i)=frc_2d(i)*tkvh(i,ke)/tkvh(i,ke1)
               END IF

!Achtung: Die Beschraenkung ist recht willkuerlich
               rat_m_2d(i)=MIN( z2, MAX( z1d2, rat_m_2d(i) ) ) !limitted (q*Sm)_p/(q*Sm)_0
               rat_h_2d(i)=MIN( z2, MAX( z1d2, rat_h_2d(i) ) ) !limitted (q*Sh)_p/(q*Sh)_0

               fac_m_2d(i)=(rat_m_2d(i)-z1)*fakt !non-stab. profile-factor for momentum
               fac_h_2d(i)=(rat_h_2d(i)-z1)*fakt !non-stab. profile-factor for scalars
             
            ELSE !Profile-factors without using the upper node
               fac_m_2d(i)=z1-rat_m_2d(i) !profile-factor for momentum
               fac_h_2d(i)=z1-rat_h_2d(i) !profile-factor for scalars
!Achtung: Analoge Beschraenkung wie bei "imode_trancnf.LT.4":
!fac_m_2d(i)=Min( fakt, MAX( -z1d2*fakt, fac_m_2d(i) ) )
!fac_h_2d(i)=Min( fakt, MAX( -z1d2*fakt, fac_h_2d(i) ) )
!Achtung: Neutrale Profile:
!fac_m_2d(i)=z0
!fac_h_2d(i)=z0
            END IF

         END DO

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(a_atm)
         DO i=ivstart, ivend

            a_atm = h_atm_2d(i)+z0m_2d(i) !turbulente Distanz auf der untersten Hauptflaeche

!           Effektive Widerstandslaengen der turb. Prandtl-Schicht:

            IF (fac_m_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
               !non-stable strat. or based on linear interpolation of profile-function
               !for the velocity scale (q*Sm):
               IF (fac_m_2d(i).EQ.z1) THEN
                  dz_0a_m(i)=z0m_2d(i)*h_atm_2d(i)/a_atm
               ELSE
                  dz_0a_m(i)=z0m_2d(i)*LOG(a_atm/(z0m_2d(i)+fac_m_2d(i)*h_atm_2d(i))) &
                                      /(z1-fac_m_2d(i))
               END IF
            ELSE !based on hyperbolic interpolation of (q*Sm) for stable stratification
               IF (imode_trancnf.EQ.3) THEN !correction if using the upper node
                  fac_m_2d(i)=fac_m_2d(i)/(fac_m_2d(i)+rat_m_2d(i))
               END IF
               dz_0a_m(i)=z0m_2d(i)*(LOG(a_atm/z0m_2d(i))-fac_m_2d(i)*h_atm_2d(i)/z0m_2d(i)) &
                                   /(z1-fac_m_2d(i))
            END IF
            IF (fac_h_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
               !non-stable strat. or using only linear interpolation of profile-function
               !for the velocity scale (q*Sh):
               IF (fac_h_2d(i).EQ.z1) THEN
                  dz_0a_h(i)=z0m_2d(i)*h_atm_2d(i)/a_atm
               ELSE
                  dz_0a_h(i)=z0m_2d(i)*LOG(a_atm/(z0m_2d(i)+fac_h_2d(i)*h_atm_2d(i))) &
                                      /(z1-fac_h_2d(i))
               END IF
            ELSE !hyperbolic interpolation of (q*Sh) for stable stratification
               IF (imode_trancnf.EQ.3) THEN !correction if using the upper node
                  fac_h_2d(i)=fac_h_2d(i)/(fac_h_2d(i)+rat_h_2d(i))
               END IF
               dz_0a_h(i)=z0m_2d(i)*(LOG(a_atm/z0m_2d(i))-fac_h_2d(i)*h_atm_2d(i)/z0m_2d(i)) &
                                   /(z1-fac_h_2d(i))
            END IF

         END DO    
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(dz_sa_m, wert, dz_s0_m)
         DO i=ivstart, ivend

!           Effektive Widerstandslaengen von den Oberflaechen bis zum Oberrand der Prandtl-Schicht
!           (unterste Modell-Hauptflaeche):

            ! inclusive lam. Grenzschicht fuer Impuls:
            wert=z1d2*dz_sg_m(i)
            dz_s0_m=wert+SQRT(wert**2+h_can_2d(i)*dz_sg_m(i))

            dz_sa_m    = dz_s0_m    + dz_0a_m(i)
            dz_sa_h(i) = dz_s0_h(i) + dz_0a_h(i)

!           Reduktionsfaktoren fuer die Bestandesschicht incl. lam. Grenzschicht:

            tfm(i)=dz_0a_m(i)/dz_sa_m
            tfh(i)=dz_0a_h(i)/dz_sa_h(i)

!           Reduktionsfaktor fuer die Verdunstung aufgrund eines um den Faktor 'rat_lam'
!           gegenueber fuehlbarer Waerme vergroesserten laminaren Transportwiderstandes:

            tfv(i)=z1/(z1+(rat_lam-z1)*dz_sg_h(i)/dz_sa_h(i))
         END DO

!        Berechnung der Erhaltungsgroessen in der Prandtl-Schicht:

!Achtung: <Korrketur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
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
!Achtung: Fehler
             ! tl_s_2d(i)=(t_g(i) - lhocp*liqs(i,ke1))
               tl_s_2d(i)=(t_g(i) - lhocp*liqs(i,ke1))/eprs(i,ke1)
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
!> Korrektur

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

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            prss(i,ke1)=ps(i)
            tmps(i,ke1)=t_g(i)
            vaps(i,ke1)=qv_s(i)
         END DO
         !$ACC END PARALLEL
         
         CALL adjust_satur_equil ( khi=ke1, ktp=ke, &
!
              i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1, i1dim=nvec, &
!
              lcalrho=.TRUE.,  lcalepr=.FALSE., lcaltdv=.TRUE.,        &
              lpotinp=.FALSE., ladjout=.FALSE.,                        &
!
              icldmod=icldm_tran,                                      &
!
              zrcpv=tur_rcpv, zrcpl=tur_rcpl,                          &
!
!Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
!         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
              prs=prss, t=tmps, qv=vaps, qc=liqs,                      &
!
              fip=tfh,                                                 &
!
              exner=eprs, rcld=rcld(:,ke1:ke1), dens=rhon(:,ke1:ke1),  &
!
              qst_t=zaux(:,ke1:ke1,3), g_tet=zaux(:,ke1:ke1,4),        &
                                       g_h2o=zaux(:,ke1:ke1,5),        &
!
              tet_l=zvari(:,ke:ke1,tet_l), q_h2o=zvari(:,ke:ke1,h2o_g),&
                                           q_liq=zvari(:,ke:ke1,liq),  &
              lacc=lzacc, opt_acc_async_queue=acc_async_queue )

!        Beachte: 
!        'zvari(:,ke1,tet_l)' und 'zvari(:,ke1,h2o_g) sind jetzt die Erhaltungsvariablen
!        am Unterrand der Prandtl-Schicht, waehrend  ta_2d' => 'zvari(:,ke,tet_l)
!        und 'qda_2d  => 'zvari(:,ke,h2o_g)' auf diese Groessen bzgl. der untersten
!        Hauptflaeche zeigen, welche zur Interpolation der Groessen an der Oberflaeche
!        auf jenen Unterrand der Prandtl-Schicht (Oberrand der Rauhigkeitsschicht)
!        benutzt werden.

!        Berechnung der benoetigten Vertikalgradienten und der TKE-Antriebe:

         IF (imode_trancnf.EQ.1) THEN 
            nvar=tet_l
         ELSE
            nvar=nvel
         END IF

         !Vertikalgradienten des Horizontalwindes:
!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            edh(i)=tfm(i)/dz_0a_m(i)
         END DO
  !LST_LOOP
         !$ACC LOOP SEQ
         DO n=1, nvel
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,ke1,n)=zvari(i,ke,n)*edh(i)
            END DO
            !Beachte: Dies ist die Darstellung ohne Nutzung der unteren Randwerte der Prandtl-Schicht
         END DO

         !Scherungs-Antrieb der TKE:
         IF (itype_sher.EQ.2 .AND. PRESENT(dwdx) .AND. PRESENT(dwdy)) THEN
            !Einschliesslich der 3D-Korrektur durch den Vertikalwind bzgl. der mittleren Hangneigung:
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=MAX( (zvari(i,ke1,u_m)+dwdx(i,ke1)*edh(i))**2 &
                              +(zvari(i,ke1,v_m)+dwdy(i,ke1)*edh(i))**2 &
                              +hdef2(i,ke1)*edh(i)**2, fc_min(i) )
            END DO
            !Beachte: dwdx(ke1), dwdy(ke1) und hdef2(ke1) beziehen sich auf die vorlaeufige Schichtdicke 1m.
         ELSE
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               frm(i,ke1)=MAX( zvari(i,ke1,u_m)**2+zvari(i,ke1,v_m)**2, fc_min(i) )
            END DO
         END IF    

         !Vertikalgradienten der dynamisch wirksamen Skalare:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            edh(i)=z1/dz_0a_h(i)
         END DO

!MR: 06.07.2018: Correction of representing the former tet_l-gradient.<
         IF (imode_trancnf.EQ.1) THEN !old version of zero-level-gradients requested
            !Transformation of Tet_l-gradient into the old form following from interpolation
            !onto the zero-level in terms of T_l (rather than Tet_l) and correcting the
            !calculated T_l-Gradient by the adiabatic laps rate:

!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               wert=zvari(i,ke1,liq)                                      !liquid water at zero-level
               val1=zvari(i,ke1,h2o_g)-zvari(i,ke1,liq)                   !spec. humid. at zero-level
               val2=eprs(i,ke1)*zvari(i,ke1,tet_l)+lhocp*zvari(i,ke1,liq) !temperature  at zero-level
               wert=val2*(z1+rvd_m_o*val1-wert)                           !virt. temp.  at zero-level

               wert=zvari(i,ke1,tet_l)*(tet_g/wert+(epr(i,ke)/eprs(i,ke1)-z1)*edh(i))
               zvari(i,ke1,tet_l)=(zvari(i,ke,tet_l)-zvari(i,ke1,tet_l))*edh(i)+wert
            END DO

         END IF

  !LST_LOOP
         !$ACC LOOP SEQ
         DO n=nvar+1,nred
!MR:06.07.2018:Correction of representing the former tet_l-gradient.>
!DIR$ IVDEP
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO i=ivstart, ivend
               zvari(i,ke1,n)=(zvari(i,ke,n)-zvari(i,ke1,n))*edh(i)
            END DO
         END DO
         !'zvari(:,ke1,n)' enthaelt jetzt die Vertikalgradienten der Erhaltungsvariablen
     
         !Auftriebs-Antrieb der TKE:
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            frh(i,ke1)=g_tet(i)*zvari(i,ke1,tet_l)+g_vap(i)*zvari(i,ke1,h2o_g)
         END DO

!        Berechnung der Stabilitaetslaengen:

         IF (it_durch.EQ.it_start .AND. lini) THEN !Startinitialisierung
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt, val1, val2, wert)
            DO i=ivstart, ivend
               IF (frh(i,ke1).GE.(z1-rim)*frm(i,ke1)) THEN
                  ! Die krit. Ri-Zahl wird ueberschritten und lm, sowie lh
                  ! werden durch lm bei der krit. Ri-Zahl angenaehert:
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

               IF (imode_tkemini.EQ.2) THEN
                  val1=con_m; val2=con_h
                  tke(i,ke1,nvor)=MAX( tke(i,ke1,nvor), val1/tkvm(i,ke1), &
                                                        val2/tkvh(i,ke1) )
               END IF          

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

! 4f)    Bestimmung des neuen SQRT(2*TKE)-Wertes:
         CALL solve_turb_budgets( khi=ke1, it_s=it_durch, it_start=it_start,                   &
                                  i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1,                &
                                  i1dim=nvec, kcm=kcm, ntur=ntur, nvor=nvor,                   &
                                  lssintact=.FALSE., lupfrclim=(imode_trancnf.EQ.1),           &
                                  lpresedr=PRESENT(edr), lstfnct=lstfnct, ltkeinp=ltkeinp,     &
                                  imode_stke=imode_tran, imode_vel_min=imode_vel_min,          &
                                  dt_tke=dt_tke, fr_tke=fr_tke,                                &
                                  tke=tke, ediss=ediss,                                        &
                                  fm2=frm(:,ke1:ke1), fh2=frh(:,ke1:ke1), ft2=frm(:,ke1:ke1),  &
                                  lsm=tkvm(:,ke1:ke1), lsh=tkvh(:,ke1:ke1),                    &
#ifdef SCLM
                                  grd=zvari(:,ke1:ke1,:),                                      &
#endif
                                  tls=len_scale(:,ke1:ke1), tvt=tvt(:,ke1:ke1),                &
                                  velmin=velmin(:),                                            &
                                  lacc=lzacc, opt_acc_async_queue=acc_async_queue              )

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(val1, val2)
         DO i=ivstart, ivend
! 4h)       Bestimmung der durch Wirkung der L-Schicht
!           korrigierten Diffusionskoeffizienten
!           und der zugehoerigen Transferkoeffizienten:

!           Unkorrigierte Diffusionskoeffizienten:

            val1=con_m; tkvm(i,ke1)=tke(i,ke1,ntur)*tkvm(i,ke1)
            val2=con_h; tkvh(i,ke1)=tke(i,ke1,ntur)*tkvh(i,ke1)

            IF (imode_lamdiff.EQ.2) THEN !surface-layer diffusion coeff. always at least at laminar value
               IF (imode_tkemini.EQ.2) THEN
                   tke(i,ke1,ntur)=tke(i,ke1,ntur)*MAX( z1, val1/tkvm(i,ke1), &
                                                            val2/tkvh(i,ke1) )
               END IF
               tkvm(i,ke1)=MAX( val1, tkvm(i,ke1) )
               tkvh(i,ke1)=MAX( val2, tkvh(i,ke1) )
            END IF

!           Belegung der Felder fuer die Transferkoeffizienten:

!Achtung: <zum Vergleich mit alter Variante:
!tvm(i)=tkvm(i,ke1)*tfm(i)/(dz_0a_m(i)*vel_2d(i))
!tvh(i)=tkvh(i,ke1)*tfh(i)/(dz_0a_h(i)*vel_2d(i))
!Achtung: Modifikation tcm -> tvm; tch -> tvh: macht Unterschiede
            tvm(i)=tkvm(i,ke1)*tfm(i)/dz_0a_m(i)
            tvh(i)=tkvh(i,ke1)*tfh(i)/dz_0a_h(i)
         END DO

         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
!tcm(i)=tvm(i); tvm(i)=tcm(i)*vel_2d(i)
!tch(i)=tvh(i); tvh(i)=tch(i)*vel_2d(i)
            tcm(i)=tvm(i)/vel_2d(i)
            tch(i)=tvh(i)/vel_2d(i)
!>zum Vergleich
         END DO

         IF (imode_trancnf.GE.4 .OR. (imode_trancnf.GE.2 .AND. it_durch.LT.it_end)) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
            DO i=ivstart, ivend
               wert=l_tur_z0(i)*SQRT(tkvm(i,ke1)*SQRT(frm(i,ke1))) !updated l_0*Ustar
               IF (ditsmot.GT.z0) THEN
                  tkr(i)=ditsmot*tkr(i)*tkvm(i,ke1) + (z1-ditsmot)*wert
               ELSE
                  tkr(i)=wert
               END IF
            END DO
         END IF

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
! 4i)       Einschraenkung von z0m_dia ueber Land:

            IF (fr_land(i) <= z1d2) THEN
               !Ueber See gibt es keinen synoptischen Garten
               z0d_2d(i)=z0m_2d(i)
            ELSE
               !Die Rauhigkeitslaenge einer SYNOP Station soll immer
               !kleiner als 10m bleiben:
               z0d_2d(i)=MIN( h_10m, z0m_dia )
            END IF
         END DO

         IF (it_durch.LT.it_end) THEN !at least one additional iteration will take place

            IF (imode_trancnf.EQ.2 .OR. imode_trancnf.EQ.3) THEN
               !new version of initializing the profile-factors using Ustar,
               !but still epressing this factor in terms of "(q*Sx)_p/(q*Sx)_0":
!DIR$ IVDEP
               !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(fakt)
               DO i=ivstart, ivend
                  fakt=h_top_2d(i)/z0m_2d(i) !(l_p-l_0)/l_0; l_0=akt*z0m
                  rat_m_2d(i)=z1+fakt*(z1-(tkr(i)      / tkvm(i,ke1)     )) !(q*Sm)_p/(q*Sm)_0
                  rat_h_2d(i)=z1+fakt*(z1-(tkr(i)*sh_0)/(tkvh(i,ke1)*sm_0)) !(q*Sh)_p/(q*Sh)_0
               END DO

            ELSEIF (imode_trancnf.GE.4) THEN
               !new version of initializing the profile-factors and already expressing
               !them in terms of "Ustar/(q*Sh)_0*(Sh(0)/Sm(0))":
               !$ACC LOOP GANG(STATIC: 1) VECTOR
               DO i=ivstart, ivend
                  rat_m_2d(i)= tkr(i)/tkvm(i,ke1)              !Ustar/(q*Sm)_0
                  rat_h_2d(i)=(tkr(i)*sh_0)/(tkvh(i,ke1)*sm_0) !Ustar/(q*Sh)_0*(sh(0)/sm(0))
               END DO
            END IF   

         END IF
         !$ACC END PARALLEL

         ! This should happen outside the OpenACC paralle region
         IF ( it_durch.LT.it_end .AND. .NOT.ltkeinp) THEN
            nvor=ntur !benutze nun aktuelle TKE-Werte als Vorgaengerwerte
         END IF

      END DO !Iteration


! 4j) Berechnung der Standardabweichnung des Saettigungsdefizites:

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)

!k->ke1
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO i=ivstart, ivend
         rcld(i,ke1)=SQRT(l_tur_z0(i)*tkvh(i,ke1)*d_h)* &
!modif: previously zvari' contained NOT gradient values
                       ABS(epr_2d(i)*qsat_dT(i)*zvari(i,ke1,tet_l)-zvari(i,ke1,h2o_g))
!modif
      END DO
! 4h) Berechnung der Enthalpie- und Impulsflussdichten sowie der EDR am Unterrand:

      IF ((lsrflux .AND. PRESENT(shfl_s)) .OR. lrunscm) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            shfl_s(i)=cp_d*rhon(i,ke1)*tkvh(i,ke1)*zvari(i,ke1,tet_l)*epr_2d(i)
            !Note: shfl_s is positive downward and belogns to the T-equation!
         END DO
!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
         IF (lsclm) THEN
            IF (SHF%mod(0)%vst.GT.i_cal .AND. SHF%mod(0)%ist.EQ.i_mod) THEN
               !measured SHF has to be used for forcing:
               shfl_s(imb)=SHF%mod(0)%val
            ELSEIF (lsurflu) THEN !SHF defined by explicit surface flux density
               SHF%mod(0)%val=shfl_s(imb)
               SHF%mod(0)%vst=MAX(i_upd, SHF%mod(0)%vst) !SHF is at least updated
            END IF
         END IF
#endif
!SCLM --------------------------------------------------------------------------------
      END IF

      IF ((lsrflux .AND. PRESENT(qvfl_s)) .OR. lrunscm) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            qvfl_s(i)=rhon(i,ke1)*tkvh(i,ke1)*zvari(i,ke1,h2o_g)
            !Note: qvfl_s is positive downward!
         END DO
      END IF   

!SCLM --------------------------------------------------------------------------------
#ifdef SCLM
      IF (lsclm) THEN
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
#endif
!SCLM --------------------------------------------------------------------------------

      !Note:
      !Both shfl_s and qvfl_s contain the flux densities of the related conserved variables now,
      ! which are equal to the usual definition only, if the cloud-water flux at the surface vanishes.

      IF (lsrflux .AND. PRESENT(umfl_s)) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            umfl_s(i)=rhon(i,ke1)*tkvm(i,ke1)*zvari(i,ke1,u_m)
         END DO
      END IF
      IF (lsrflux .AND. PRESENT(vmfl_s)) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            vmfl_s(i)=rhon(i,ke1)*tkvm(i,ke1)*zvari(i,ke1,v_m)
         END DO
      END IF

      IF (PRESENT(edr)) THEN
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
            edr(i,ke1)=tke(i,ke1,ntur)**3/(d_m*l_tur_z0(i))
         END DO
      END IF
      !$ACC END PARALLEL !!

!----------------------------------------

! 5)  Diagnose der meteorologischen Groessen im 2m- und 10m-Niveau:

      IF (lnsfdia) THEN !diagnostics at near surface level required at this place

!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR
      DO i=ivstart, ivend
         IF (itype_diag_t2m.EQ.2) THEN !using an exponetial rougness layer profile
           z2m_2d (i) = h_2m -h_can_2d(i) !2m ueber dem Bodenniveau des Bestandes
           z10m_2d(i) = h_10m-z0m_2d(i)   !Hoehe, in der turbulente Distanz 10m betraegt
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

#ifdef _OPENACC
      CALL diag_level_gpu(ivstart, ivend, z2m_2d, k_2d, hk_2d, hk1_2d, &
         lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#else
      CALL diag_level(ivstart, ivend, z2m_2d, k_2d, hk_2d, hk1_2d, &
         lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#endif

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      IF (itype_diag_t2m.EQ.2) THEN !using an exponential rougness layer profile

         val2=z1/epsi
!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(val1, fakt, wert)
         DO i=ivstart, ivend
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               IF (z2m_2d(i).LT.z0) THEN
!                 2m-Niveau liegt innerhalb der Bestandesschicht
!                 mit exponentiellen Vertikalprofilen:

!modif: Vereinfachung der Sicherheitsabfrage:
                  val1=z2m_2d(i)/dz_s0_h(i)
                  IF (-val1.LE.val2) THEN
                    fakt=dz_s0_h(i)/dz_sa_h(i)
                    fakt=MIN( z1, MAX( z0, fakt*EXP(val1) ) )
                  ELSE
                    fakt=z0
                  ENDIF
!modif
               ELSE
!                 2m-Niveau liegt innerhalb der Modell_Prandtl-Schicht
!                 mit logarithmischen Vertikalprofilen:

                  IF (ABS(z1-fac_h_2d(i)) < epsi ) THEN
                     wert=z0m_2d(i)*z2m_2d(i)/(z2m_2d(i)+z0m_2d(i))
                  ELSEIF (fac_h_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
                     !non-stable strat. or using only linear interpolation of profile-function
                     !for the velocity scale (q*Sh):
                     wert=z0m_2d(i)*LOG((z2m_2d(i)+z0m_2d(i))/            &
                         (z0m_2d(i)+fac_h_2d(i)*z2m_2d(i)))/(z1-fac_h_2d(i))
                  ELSE !hyperbolic interpolation of (q*Sh) for stable stratification
                     wert=z0m_2d(i)*(LOG(z2m_2d(i)/z0m_2d(i)+z1)-fac_h_2d(i)*z2m_2d(i)/z0m_2d(i)) &
                                   /(z1-fac_h_2d(i))
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
               prss(i,ke1) = ps(i)

            END IF
         END DO

      ELSE !using only a logarithmic profile above a SYNOP lawn

!DIR$ IVDEP
         !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(z0d, a_atm, a_2m, fr_sd_h, val1, val2, fakt)
         DO i=ivstart, ivend
            IF (k_2d(i).EQ.ke) THEN
!              2m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z0d=z0d_2d(i)
               a_atm=h_atm_2d(i)+z0d
               a_2m=h_2m+z0d

!              Dimensionsloser Widerstand des Rauhigkeitsbestandes der SYNOP-Wiese,
!              wobei die turbulente Geschwindigkeitsskala und der Oberflaechenindex
!              gegenueber dem mittleren Rauhigkeitsbestand gleich bleiben:

               fr_sd_h=MAX( z0, dz_s0_h(i)/z0m_2d(i)+LOG(z0d/z0m_2d(i)) )

               IF (imode_trancnf.LT.4) THEN !further corrected profile-factor using an upper node
                  fac_h_2d(i)=(rat_h_2d(i)-z1)*z0d/h_top_2d(i)
!Achtung:
!'fac_h_2d' sollte wohl eher auch eine Profil-Konstante sein
! und nicht mit "z0d/z0m_2d" skaliert werden (also hier nicht mehr veraendert werden)
               END IF

               IF (fac_h_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sh):
                  IF (fac_h_2d(i).EQ.z1) THEN
                     val1=fr_sd_h+h_2m/a_2m
                     val2=fr_sd_h+h_atm_2d(i)/a_atm
                  ELSE
                     fakt=z1/(z1-fac_h_2d(i))
                     val1=fr_sd_h+LOG(a_2m /(z0d+fac_h_2d(i)*h_2m       ))*fakt
                     val2=fr_sd_h+LOG(a_atm/(z0d+fac_h_2d(i)*h_atm_2d(i)))*fakt
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sh) for stable stratification
                  fakt=z1/(z1-fac_h_2d(i))
                  val1=fr_sd_h+(LOG(a_2m /z0d)-fac_h_2d(i)*h_2m       /z0d)*fakt
                  val2=fr_sd_h+(LOG(a_atm/z0d)-fac_h_2d(i)*h_atm_2d(i)/z0d)*fakt
               END IF

               fakt=val1/val2

               IF (imode_syndiag.EQ.1) THEN
                  tmps(i,ke1) = t_g(i) + (t(i,ke)-t_g(i))*fakt &
                              + tet_g*(h_atm_2d(i)*fakt-h_2m)
                  vaps(i,ke1)= qv_s(i) + (qv(i,ke)-qv_s(i))*fakt
               ELSE
                  tmps(i,ke1) = fakt*ta_2d(i) + (z1-fakt)*tl_s_2d(i)/epr_2d(i)
                  vaps(i,ke1) = qt_s_2d(i) + fakt*(qda_2d(i)-qt_s_2d(i))
               END IF
               prss(i,ke1) = ps(i)

            END IF
         END DO
      END IF
      
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
!test: different interpol. weight
!  fakt=(z2m_2d(i)-hk1_2d(i))/(hk_2d(i)-hk1_2d(i))
!test

            IF (imode_syndiag.EQ.1) THEN
               tmps(i,ke1)= t(i,k1)+fakt*( t(i,k2)- t(i,k1))
               vaps(i,ke1)=qv(i,k1)+fakt*(qv(i,k2)-qv(i,k1))
            ELSE
               val2=qv(i,k2)+qc(i,k2)
               val1=qv(i,k1)+qc(i,k1)

               vaps(i,ke1)=val1+fakt*(val2-val1)

               IF (icldm_tran.EQ.-1) THEN !no water phase change possible
                  val2=t(i,k2)/epr(i,k2)
                  val1=t(i,k1)/epr(i,k1)
               ELSE !water phase changes are possible
                  val2=(t(i,k2)-lhocp*qc(i,k2))/epr(i,k2)
                  val1=(t(i,k1)-lhocp*qc(i,k1))/epr(i,k1)
               END IF
               tmps(i,ke1)=val1+fakt*(val2-val1)

               rcl_2d(i)=rcld(i,k1)+fakt*(rcld(i,k2)-rcld(i,k1))
            END IF
            prss(i,ke1)=prs(i,k1)
         ELSEIF (z2m_2d(i).LE.z0) THEN
            rcl_2d(i)=rcld(i,ke1)
         ELSE
            fakt=z2m_2d(i)/h_atm_2d(i)
            rcl_2d(i)=rcld(i,ke1)+fakt*(rcld(i,ke)-rcld(i,ke1))
         END IF
         !Note: In the case "k_2d(i).EQ.ke" 'tmps' and 'prss' have already been
         !      caclculated above.
      END DO

      !Druck im 2m-Niveau:
!DIR$ IVDEP
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(wert)
      DO i=ivstart, ivend
         wert=tmps(i,ke1)*(z1+rvd_m_o*vaps(i,ke1)) !angenaeherte virt. Temp.
         prss(i,ke1)=prss(i,ke1)                &  !Druck
                    *EXP(-(z2m_2d(i)-hk1_2d(i)) &
                    *grav/(r_d*wert))
      END DO
      !$ACC END PARALLEL
      IF (imode_syndiag.EQ.1) THEN
!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG(STATIC: 1) VECTOR
         DO i=ivstart, ivend
             t_2m(i)=tmps(i,ke1)
            qv_2m(i)=vaps(i,ke1)
         END DO
         !$ACC END PARALLEL
      ELSE

!        Berechnung der zugehoerigen Modell- und Feuchtevariablen im 2m-Niveau
!        aus den Erhalturngsvariablen.

         CALL adjust_satur_equil ( khi=ke1, ktp=ke, &
!
              i_st=ivstart, i_en=ivend, k_st=ke1, k_en=ke1, i1dim=nvec,        &
!
              lcalrho=.FALSE., lcalepr=.TRUE., lcaltdv=.FALSE.,                &
              lpotinp=.TRUE. , ladjout=.TRUE.,                                 &
!
              icldmod=icldm_tran,                                              &
!
              zrcpv=tur_rcpv, zrcpl=tur_rcpl,                                  &
!
!Achtung: Korrektur: Konsistente Behandlung der unteren Null-Fluss-Randbedingung fuer qc
!         und an COSMO-Version angepasste Interpretation von "icldmod=-1":
              prs=prss, t=tmps, qv=vaps, qc=liqs,                              &
!
              psf=ps,                                                          &
!
              exner=eprs, rcld=rclc,                                           &
!
              qst_t=zaux(:,ke1:ke1,3), g_tet=zaux(:,ke1:ke1,4),                &
                                       g_h2o=zaux(:,ke1:ke1,5),                &
!
              tet_l=zvari(:,ke:ke1,tet_l), q_h2o=zvari(:,ke:ke1,h2o_g),        &
                                           q_liq=zvari(:,ke:ke1,liq),          &
              lacc=lzacc, opt_acc_async_queue=acc_async_queue )
!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR
         DO i=ivstart, ivend
             t_2m(i)=zvari(i,ke1,tet_l)
            qv_2m(i)=zvari(i,ke1,h2o_g)
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
      ELSE

!        Finale 2m-Diagnose:

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR PRIVATE(patm, fakt, wert)
         DO i=ivstart, ivend
!Achtung: Macht minimale Unterschiede
            patm=prss(i,ke1)*qv_2m(i) &
                /(rdv+(z1-rdv)*qv_2m(i))          !Wasserdampfdruck
!patm=prss(i,ke1) &
!    /(rdv/MAX(qv_2m(i),epsi)+(z1-rdv))           !Wasserdampfdruck (alt)

            fakt=patm/zpsat_w( t_2m(i) )
            rh_2m(i)=100.0_wp*MIN( fakt, z1 ) !relative Feuchte
   
            !UB: old formulation           
            wert=LOG(patm/b1)
            !UB: For dry atmosphere, the Teten's formula is not defined at vapor 
            !    pressure patm=0 because of log(patm/b1). However, it converges 
            !    to a dew point of the value of parameter b4w in case we impose 
            !    a very small positive lower bound on  patm:
            ! erst mal wieder weggenommen, da man den Absturz will:  
            !wert=LOG(MAX(patm,1.0E-16_wp)/b1)

            td_2m(i)=MIN( (b2w*b3-b4w*wert) &
                      /(b2w-wert), t_2m(i) )      !Taupunktstemperatur
         END DO
         !$ACC END PARALLEL

!        Diagnose der 10m-Groessen:

#ifdef _OPENACC
         CALL diag_level_gpu(ivstart, ivend, z10m_2d, k_2d, hk_2d, hk1_2d, &
            lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#else
         CALL diag_level(ivstart, ivend, z10m_2d, k_2d, hk_2d, hk1_2d, &
            lacc=lzacc, opt_acc_async_queue=acc_async_queue)
#endif

!DIR$ IVDEP
         !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(vel2_2d) ASYNC(acc_async_queue) IF(lzacc)
         !$ACC LOOP GANG VECTOR PRIVATE(z0d, a_atm, a_10m, val1, val2, fakt, k1, k2, wert)
         DO i=ivstart, ivend

            IF (k_2d(i).EQ.ke) THEN

!              10m-Niveau unterhalb der untersten Modell-Hauptflaeche
!              in der lokalen Prandtl-Schicht mit Rauhigkeitslaenge z0d:

               z0d=z0d_2d(i)
               a_atm=h_atm_2d(i)+z0d
               a_10m=h_10m+z0d

               IF (imode_trancnf.LT.4) THEN !further corrected profile-factor using an upper node
                  fac_m_2d(i)=(rat_m_2d(i)-z1)*z0d/h_top_2d(i)
!Achtung:
!'fac_m_2d' sollte wohl eher auch eine Profil-Konstante sein
! und nicht mit "z0d/z0m_2d" skaliert werden (also hier nicht mehr veraendert werden)
               END IF

               IF (fac_m_2d(i).GE.z0 .OR. imode_trancnf.LT.3) THEN
                  !non-stable strat. or based on linear interpolation of profile-function
                  !for the velocity scale (q*Sm):
                  IF (fac_m_2d(i).EQ.z1) THEN
                     val1=h_10m/a_10m
                     val2=h_atm_2d(i)/a_atm
                  ELSE
                     val1=LOG(a_10m/(z0d+fac_m_2d(i)*h_10m))
                     val2=LOG(a_atm/(z0d+fac_m_2d(i)*h_atm_2d(i)))
                  END IF
               ELSE !based on hyperbolic interpolation of (q*Sm) for stable stratification
                  val1=(LOG(a_10m/z0d)-fac_m_2d(i)*h_10m      /z0d)
                  val2=(LOG(a_atm/z0d)-fac_m_2d(i)*h_atm_2d(i)/z0d)
               END IF
   
               fakt=val1/val2
   
               u_10m(i)=vel1_2d(i)*fakt
               v_10m(i)=vel2_2d(i)*fakt

            ELSE
!              10m-Niveau liegt oberhalb der untersten Hauptflaeche und wir nutzen
!              trotz der allgemein zwischen atm. Modellneveaus als gueltig angenommen
!              lenearen Profile der progn. Modellvariablen eine logarithm. Interpolation:

               k2=k_2d(i); k1=k2+1

               fakt=z1/(hk1_2d(i)+z0d_2d(i))
               wert=(h_10m   +z0d_2d(i))*fakt
               fakt=(hk_2d(i)+z0d_2d(i))*fakt
               fakt=LOG(wert)/LOG(fakt)
!test: different interpol. weight
!  fakt=(z10m_2d(i)-hk1_2d(i))/(hk_2d(i)-hk1_2d(i))
!test
               val2=u(i,k2)
               val1=u(i,k1)
               u_10m(i)=val1+fakt*(val2-val1)
   
               val2=v(i,k2)
               val1=v(i,k1)
               v_10m(i)=val1+fakt*(val2-val1)
            END IF

         END DO
         !$ACC END PARALLEL

      END IF
      !Note: '(u, v)_10m' always belong to mass points!

      END IF !in case of ".NOT.lnsfdia" this kind of diagnostics is done at another place
      
!DIR$ IVDEP
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(acc_async_queue) IF(lzacc)
      !$ACC LOOP GANG VECTOR PRIVATE(velo, wert, fakt)
      DO i=ivstart, ivend

!        Diagnose von gz0 (fuer den naechsten Zeitschritt)
!        ueber Wasserflaechen mit der (angepassten) Charnock-Formel

         IF (fr_land(i) < z1d2) THEN

           ! Use ice surface roughness or open-water surface roughness
           ! according to lo_ice
            IF ( lo_ice(i) ) THEN
               ! Ice-covered grid box
               gz0(i)=grav*z0_ice
            ELSE !water covered surface
               velo=(tke(i,ke1,ntur)+tke(i,ke,nvor))*z1d2
!Achtung: Die 'epsi'-Beschraenkung ist recht willkuerlich und fehlt in COSMO-Version!
!Achtung: Modifikation tcm -> tvm: macht Unterschiede
               wert=MAX( epsi, tvm(i)*SQRT(vel_2d(i)**2+velo**2) ) !effective Ustar**2
               IF (imode_charpar.EQ.1) THEN !constant Charnock-Parameter
                  fakt=alpha0
               ELSE
                  IF (lini .AND. .NOT.lnsfdia) THEN
                     velo=vel_2d(i)
                  ELSE
                     velo=SQRT(u_10m(i)**2+v_10m(i)**2)
                  END IF
!US from ICON version 044780ed>
                  IF (depth_lk(i) > z0) THEN
                    ! enhanced Charnock parameter over lakes, parameterizing a non-equlibrium wave spectrum
                    fakt = 0.1_wp
                  ELSE
                    fakt=alpha0_char(velo)
                  ENDIF
!US<
               END IF
               wert=MAX( grav*len_min, fakt*wert+grav*alpha1*con_m/SQRT(wert) )
               IF (ditsmot.GT.z0) THEN
                  gz0(i)=ditsmot*gz0(i)+(z1-ditsmot)*wert
               ELSE
                  gz0(i)=wert
               END IF
            END IF
         END IF
      END DO
      !$ACC END PARALLEL

      !$ACC END DATA ! from acc data present
      !$ACC END DATA ! from acc data create
! Just do some checkout prints:
  IF (ldebug) THEN
    DO i = ivstart, ivend
      IF (i ==  1 .AND. iblock == 1 .AND. my_cart_id == 0) THEN
        WRITE(*,'(A       )') '  '
        WRITE(*,'(A,3I5   )') 'TURB-DIAGNOSIS transfer:  iblock = ', iblock, i

        WRITE(*,'(A       )') ' Some output parameters:'
        WRITE(*,'(A,F28.16)') '   gz0              :  ', gz0         (i)
        WRITE(*,'(A,F28.16)') '   tcm              :  ', tcm         (i)
        WRITE(*,'(A,F28.16)') '   tch              :  ', tch         (i)
        WRITE(*,'(A,F28.16)') '   tvm              :  ', tvm         (i)
        WRITE(*,'(A,F28.16)') '   tvh              :  ', tvh         (i)
        WRITE(*,'(A,F28.16)') '   tfm              :  ', tfm         (i)
        WRITE(*,'(A,F28.16)') '   tfh              :  ', tfh         (i)
        WRITE(*,'(A,F28.16)') '   tfv              :  ', tfv         (i)
        WRITE(*,'(A,F28.16)') '   tkr              :  ', tkr         (i)
        WRITE(*,'(A,F28.16)') '   tke   ke         :  ', tke         (i,ke,ntur)
        WRITE(*,'(A,F28.16)') '   tke   ke1        :  ', tke         (i,ke1,ntur)
        WRITE(*,'(A,F28.16)') '   tkvm  ke         :  ', tkvm        (i,ke)
        WRITE(*,'(A,F28.16)') '   tkvm  ke1        :  ', tkvm        (i,ke1)
        WRITE(*,'(A,F28.16)') '   tkvh  ke         :  ', tkvh        (i,ke)
        WRITE(*,'(A,F28.16)') '   tkvh  ke1        :  ', tkvh        (i,ke1)
!       WRITE(*,'(A,F28.16)') '   edr   ke         :  ', edr         (i,ke)
!       WRITE(*,'(A,F28.16)') '   edr   ke1        :  ', edr         (i,ke1)
        WRITE(*,'(A,F28.16)') '   t_2m             :  ', t_2m        (i)
        WRITE(*,'(A,F28.16)') '   qv_2m            :  ', qv_2m       (i)
        WRITE(*,'(A,F28.16)') '   td_2m            :  ', td_2m       (i)
        WRITE(*,'(A,F28.16)') '   rh_2m            :  ', rh_2m       (i)
        WRITE(*,'(A,F28.16)') '   u_10m            :  ', u_10m       (i)
        WRITE(*,'(A,F28.16)') '   v_10m            :  ', v_10m       (i)
        WRITE(*,'(A       )') '  '
      ENDIF
    ENDDO
  ENDIF




!==============================================================================

CONTAINS

!********************************************************************************

!+ Module procedure diag_level for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level (i_st, i_en, zdia_2d, k_2d, hk_2d, hk1_2d, lacc, opt_acc_async_queue)
   INTEGER, INTENT(IN) :: &
!
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=wp), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER, INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                  !to be used for near surface diagnostics

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
            !            lcheck=lcheck .OR. .TRUE. !for this point or any previous one the
            !$ACC ATOMIC WRITE
            lcheck=.TRUE.             !for this point or any previous one the
                                      !diagnostic level is above the current layer
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

!********************************************************************************

!+ Module procedure diag_level for computing the upper level index
!+ used for near surface diganostics

SUBROUTINE diag_level_gpu (i_st, i_en, zdia_2d, k_2d, hk_2d, hk1_2d, lacc, opt_acc_async_queue)
   INTEGER, INTENT(IN) :: &
!
      i_st, i_en  !start end end indices of horizontal domain

   REAL (KIND=wp), INTENT(IN) :: &
!
      zdia_2d(:)  !diagnostic height

   INTEGER, INTENT(INOUT) :: &
!
      k_2d(:)     !index field of the upper level index
                  !to be used for near surface diagnostics

   REAL (KIND=wp), INTENT(INOUT) :: &
!
     hk_2d(:), & !mid level height above ground belonging to 'k_2d'
     hk1_2d(:)    !mid level height above ground of the previous layer (below)

   LOGICAL, OPTIONAL, INTENT(IN) :: lacc
   INTEGER, OPTIONAL, INTENT(IN) :: opt_acc_async_queue

   LOGICAL :: lzacc
   INTEGER :: acc_async_queue

   INTEGER :: i, k, ke1_war

   LOGICAL :: lcheck

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


END SUBROUTINE turbtran

!==============================================================================

END MODULE turb_transfer
