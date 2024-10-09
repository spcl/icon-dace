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

! Source module for computing implicit vertical diffusion

MODULE turb_vertdiff

!-------------------------------------------------------------------------------

! Modules used:

#ifdef _OPENMP
  USE omp_lib,            ONLY: omp_get_thread_num
#endif

!-------------------------------------------------------------------------------
! Parameter for precision
!-------------------------------------------------------------------------------

USE mo_kind,         ONLY :   &
    wp              ! KIND-type parameter for real variables

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
    lh_v     => alv         ! evaporation heat

!-------------------------------------------------------------------------------
! Turbulence data (should be the same in ICON and COSMO)
!-------------------------------------------------------------------------------

USE turb_data, ONLY : &

    ! used derived types
    modvar, turvar, varprf, & !

! Switches controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------

    lnonloc,      & ! nonlocal calculation of vertical gradients used for turb. diff.

!   for semi-implicit vertical diffusion:
    lsflcnd,      & ! lower flux condition
    ldynimp,      & ! dynamical calculation of implicit weights
    lprecnd,      & ! preconditioning of tridiagonal matrix

! Selectors controlling the turbulence model, turbulent transfer and diffusion:
! ----------------------------------------------------------------------------
!
    ilow_def_cond,& !type of the default condition at the lower boundary
                    ! 1: zero surface gradient
                    ! 2: zero surface value

    ! numbers and indices

    nvel    ,     & ! number of velocity components
    naux    ,     & ! number of auxilary variables
    nmvar   ,     & ! number of included prognostic model-variables
    ntyp    ,     & ! number of variable types (mom) und (sca)
    ndim    ,     & !
    mom     ,     & ! index for a momentum variable
    sca     ,     & ! index for a scalar   variable
    u_m     ,     & ! index for mass centered zonal      velocity-compont
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
    vert_grad_diff,                  &
    bound_level_interp,              &
    zexner

!SCLM---------------------------------------------------------------------------
#ifdef SCLM
USE data_1d_global, ONLY : &
    lsclm, latmflu, i_cal, i_mod, imb, &
    SHF, LHF
#endif
!SCLM---------------------------------------------------------------------------

!===============================================================================

IMPLICIT NONE

PUBLIC  :: vertdiff

!===============================================================================

REAL (KIND=wp), PARAMETER :: &
    z0 = 0.0_wp,    &
    z1 = 1.0_wp      

!===============================================================================

CONTAINS

!===============================================================================

#  define err_args

SUBROUTINE vertdiff ( &
!
          itnd, lscadif, lum_dif, lvm_dif,   &
                lsfluse, lqvcrst, lrunscm,   &
          ldoexpcor, ldocirflx,              &
          l3dflxout,                         &
!
          dt_var, nvec, ke, ke1,             &
!
          kcm, kstart_cloud,                 &
          iblock, ivstart, ivend,            &
!
          hhl, dp0, r_air, zvari,            &
!
          t_g, qv_s, ps,                     &
          u, v, t, qv, qc, prs,              &
          rhoh, rhon, epr,                   &
!
          impl_weight,                       &
          ptr, ndtr,                         &
!
          tvm, tvh, tkvm, tkvh,              &
          u_tens, v_tens, t_tens,            &
          qv_tens, qc_tens,                  &
          qv_conv,                           &
!
          umfl_s, vmfl_s, shfl_s, qvfl_s     &
          err_args)

!-------------------------------------------------------------------------------
!
! 
! Description:
!
! Method:
!
!-------------------------------------------------------------------------------

! Declarations:
!--------------

! Formal Parameters:
!-------------------

! 0. Parameters controlling the call of 'organize_turbdiff':

LOGICAL, INTENT(IN) :: &

  lum_dif,      & !running vertical gradient diffusion of horizontal u-momenum
  lvm_dif,      & !running vertical gradient diffusion of horizontal v-momenum
  lscadif,      & !running vertical gradient diffusion of scalar properties

  lsfluse,      & !use explicit heat flux densities at the surface
  lqvcrst,      & !qv-flux-divergence reset requested (only if 'qv_conv' is present)

  lrunscm,      & !a Single Column run (default: FALSE)
  ldoexpcor,    & !consider explicit warm-cloud correct. for turb. scalar fluxes
  ldocirflx,    & !consider circulation heat-flux
  l3dflxout       !3D-output of effective vertical fluxes of dynamically active scalars requested

REAL (KIND=wp), INTENT(IN) :: &

  dt_var          !time step for ordinary prognostic variables

INTEGER,        INTENT(IN) :: &

  itnd            !type of tendency cons. (0: no, 1: in implicit vertical diffusion equation
                  !                               2: by adding to current profile before vertical diffusion
                  !                               3: by using corrected virtual vertical profiles

INTEGER,        INTENT(IN) :: &

! Horizontal and vertical sizes of the fields and related variables:
! --------------------------------------------------------------------

  nvec,         & ! number of grid points in the nproma-vector
  ke,           & ! index of the lowest main level
  ke1,          & ! index of the lowest half level (=ke+1)
  kcm,          & ! level index of the upper canopy bound
  iblock,       & ! index of the current block

  kstart_cloud    ! start level index for vertical diffusion of cloud water

  !Note: Through 'kstart_cloud' vertical diffusion of the respective properties
  !       is artificially restricted above this very level

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

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: &
!
  dp0             ! pressure thickness of layer                   (pa )

 REAL (KIND=wp), DIMENSION(:,:,0:), TARGET, INTENT(INOUT) :: & 
  zvari           ! input: effective (possibly non-local) vertical gradients of regular model variables
                  !        (including the effect of turbulent saturation adjustment);
                  !        used to calculate the vertical divergence of non-gradient vertical fluxes, 
                  !        if "ldoexpcor=T" or "lnonloc" or "ldocirflx=T" is valid!
                  ! outpt: effective vertical gradients as resulting from the applied semi-implicit procedure
                  !        for vertical diffusion.

REAL (KIND=wp), DIMENSION(:,kcm-1:), TARGET, OPTIONAL, INTENT(IN) :: &
  r_air           ! log of air containing fraction of a gridbox inside
                  ! the canopy                                          (1)

! Fields for surface values and soil/canopy model variables:
! ------------------------------------------------------------

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(IN) :: &
!
  ps,           & ! surface pressure                              ( pa  )
  qv_s,         & ! specific water vapor content on the surface   (kg/kg)
  t_g             ! weighted surface temperature                  (  k  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
!
! Atmospheric model variables:
! ---------------------------------
!
  u,           & ! zonal wind speed       (at mass positions)    ( m/s )
  v,           & ! meridional wind speed  (at mass positions)    ( m/s )
  t,           & ! temperature                                   (  k  )
  qv,          & ! specific water vapor content                  (kg/kg)
  qc             ! specific cloud water content                  (kg/kg)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
  prs            ! atmospheric pressure                          ( pa  )

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(IN) :: &
!
  rhoh,        & ! total density of air                          (kg/m3)
  epr            ! exner pressure                                 (1)

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &
  rhon           ! total density of air (at half levels)         (kg/m3)

REAL (KIND=wp), DIMENSION(:), INTENT(IN) :: &
!
  impl_weight    ! profile of precalculated implicit weights 

TYPE (modvar), OPTIONAL, INTENT(INOUT) :: ptr(:) ! passive tracers
INTEGER                , INTENT(IN)    :: ndtr   ! number of tracers to be diffused

REAL (KIND=wp), DIMENSION(:), TARGET, INTENT(INOUT) :: &
!
! Diagnostic surface variable of the turbulence model:
! -----------------------------------------------------

! turbulent (transfer) velocity scales at the surface
  tvm,          & ! for momentum                                  ( m/s)
  tvh             ! for heat and moisture                         ( m/s)

  !Notice that 'tcm' and 'tch' are dispensable. The common use of the related
  !vecolities  'tvm' and 'tvh' makes live much easier!!               

! Atmospheric variables of the turbulence model:
! ------------------------------------------------

REAL (KIND=wp), DIMENSION(:,:), TARGET, INTENT(INOUT) :: &

  tkvm,         & ! turbulent diffusion coefficient for momentum  (m2/s )
  tkvh            ! turbulent diffusion coefficient for heat      (m2/s )
                     ! (and other scalars)

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
! Tendency fields for the prognostic variables:
! -----------------------------------------------
!
  u_tens,       & ! u-tendency                                    ( m/s2)
  v_tens,       & ! v-tendency                                    ( m/s2)
  t_tens,       & ! t-tendency                                    ( K/s )
  qv_tens,      & ! qv-tendency                                   ( 1/s )
  qc_tens         ! qc-tendency                                   ( 1/s )

REAL (KIND=wp), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: &
!
  qv_conv         ! qv-flux-convergence                            ( 1/s )

REAL (KIND=wp), DIMENSION(:), OPTIONAL, TARGET, INTENT(INOUT) :: &
!
  shfl_s,       & ! sensible heat flux at the surface             (W/m2)    (positive downward)
  qvfl_s,       & ! water vapor   flux at the surface             (kg/m2/s) (positive downward)
  umfl_s,       & ! u-momentum flux at the surface                (N/m2)    (positive downward)
  vmfl_s          ! v-momentum flux at the surface                (N/m2)    (positive downward)

!-------------------------------------------------------------------------------
!Local Parameters:
!-------------------------------------------------------------------------------

! Local logicals:

  LOGICAL ::   &
    ldovardif, & !berechne (teil-)implizite Vert.diff von Mod.var 1-ter Ordnung
    ldogrdcor, & !mache Gradientkorrektur bei Berechnung der vertikalen Diffusion
    linisetup, & !initiales setup bei impliziter Vertikaldiffusion
    lnewvtype, & !neuer Variablentyp muss vorbereitet werden
    lsflucond, & !untere Flussrandbedingung
    leff_flux    !calculation of effective flux density required

  LOGICAL ::  &
    ltend(nmvar+ndtr), &  !calculation of tendencies required
    lsfli(nmvar+ndtr)     !surface value input is a flux density instead of a concentration

! Local integers:

  INTEGER ::   &
    n,m,       & !Indices fuer diverse Schleifen
    i,k,       & !horizontaler und vertikaler Laufindex

    igrdcon,   & !Index fuer Modus der Gradientberuecksichtigung
    itndcon,   & !Index fuer Modus der  Tendenzberuecksichtigung
    ivtype,    & !Index fuer Variablentyp

    ncorr,     & !Start-Index der Variablen mit Gradientkorrektur
    mcorr,     & !End-Index   der Variablen mit Gradientkorrektur
    ndiff,     & !number of 1-st order variables
    ntrac,     & !number of included passive tracers in tracer vector 'ptr'
    itrac,     &  !index of a current tracer provided for vertical turbulent diffusion

    k_st_up,   & !First used level-index of an individual vertical variable-profile (bottom-up integration)
    k_st_pp      !First level-index for a present individual vertical variable-profile

! Local reals:

  REAL (KIND=wp) ::  &
    fakt, wert,      & !for any factors
    virt               !z1+(Rv/Rd-1)*qv-qc

! Time increment and inverse time increment of ordinary prognostic variables:
  REAL (KIND=wp), TARGET :: &
    tinc(nmvar+ndtr)


  REAL (KIND=wp), TARGET ::   &
    ! Already allocated variables of the turbulence model used for various purposes :

    eprs     (nvec,ke1:ke1),  & ! surface Exner-factor

    len_scale(nvec,ke1),      & ! storage for    diffusion        momentum (diff_mom)

    zaux     (nvec,ke1,ndim), & ! storage for 1: discretization   momentum (disc_mom)
                                !             2: explcit  part of diff_mom (expl_mom)
                                !             3: implicit part of diff_mom (impl_mom)
                                !             4: inverted momentum         (invs_mom)
                                !             5: diffusion depth           (diff_dep)

    frh      (nvec,ke1),      & ! storage for an inversion factor          (invs_fac)
    frm      (nvec,ke1),      & ! storage for a  scaling   factor          (scal_fac)

    dicke    (nvec,ke1),      & ! storage mainly for diffusion tendencies  (dif_tend)
    hlp      (nvec,ke1)         ! storage for current variable profiles    (cur_prof)

  REAL (KIND=wp), POINTER, CONTIGUOUS :: &
!
!   Pointer fuer Tendenzfelder:
    utens(:,:), vtens(:,:), &
    ttens(:,:), qvtens(:,:), qctens(:,:)

! Note:
! The following buffers wouldn't be necessary, if the related pointers above
! were allowed to be allocated at run time:

  REAL (KIND=wp), DIMENSION(:,:), POINTER, CONTIGUOUS :: &
!
    cur_prof, dvar_av, dvar_at, vtyp_tkv

  REAL (KIND=wp), POINTER, CONTIGUOUS :: &
  dvar_sv(:)

! Local parameters of derived types:

  TYPE (modvar) :: dvar(nmvar+ndtr) !model variables to be diffused

  TYPE (turvar) :: vtyp(ntyp)       !variable types (momentum and scalars)

! TYPE (varprf) :: pvar(0) !vertical variable profile at main- and boundary levels
  TYPE (varprf) :: pvar(0:0) !vertical variable profile at main- and boundary levels

! technical parameters:

LOGICAL :: ldebug=.FALSE.

INTEGER :: my_cart_id, my_thrd_id

!---- End of header ------------------------------------------------------------

!===============================================================================

!All variables and their tendencies are defined at horizontal mass positions.

  ldogrdcor=(ldoexpcor.OR.lnonloc.OR.ldocirflx)  !gradient correction has to be done
  ldovardif=(lum_dif .OR. lvm_dif .OR. lscadif)  !some variable has to be diffused

  IF (PRESENT(ptr)) THEN !passive tracers are present
    ntrac = UBOUND(ptr,1)
  ELSE
    ntrac=0
  END IF
  IF (ndtr.GT.ntrac) THEN
    CALL finish('', 'ERROR *** Number of tracers larger than dimension of tracer vector ''prt'' ***')
  END IF

  ndiff=nmvar+ndtr !number of 1-st order variables used in the turbulence model
                   !note that cloud ice is treated like a passive trace here

  !According to the setting in 'turb_data' it holds:
  !     nmvar = nscal+nvel: number of model variables being dynamically active for turbulence
  !             nvel  = 2    active horizontal wind components:  'u_m', 'v_m'
  !                          u_m = 1:     zonal      wind
  !                          v_m = 2:     meridional wind
  !             nscal = 3    active 1-st order scalar variables: 'tem', 'vap', 'liq'
  !                          tem   = 3:   temperature
  !                          vap   = 4:   water vapor mixing ration
  !                          liq   = 5:   liquid water
  !
  !                but also: tem_l = 3:   liquid water temperature
  !                          tet   = 3:   potential temperature
  !                          tet_l = 3:   moist (liquid water?) potential temperature
  !                          h2o_g = 4:   total water content
  !Further, according to the INPUT list, it holds:
  !     ndtr: number of (in this sense) passive tracers to be additionally diffused

  !Begin of GPU data region
  !$ACC DATA &
  !$ACC   CREATE(len_scale, frh, frm, eprs, dicke, hlp, zaux, tinc)

  ltend(u_m)=PRESENT(u_tens)
  IF (ltend(u_m)) THEN !calculation of tendencies required
    utens => u_tens    !'utens' points to the tendency
  ELSE                 !update of ordinary prognostic variables required
    utens => u         !'utens' points to the prognostic variables
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

  lsfli(:)=.FALSE. !surface values are concentrations by default

  dvar(u_m)%av  => u  ; dvar(u_m)%at => utens  ; dvar(u_m)%sv => NULL() ; dvar(u_m)%kstart = 1
  dvar(v_m)%av  => v  ; dvar(v_m)%at => vtens  ; dvar(v_m)%sv => NULL() ; dvar(v_m)%kstart = 1

!Note: Use                                       dvar(u_m)%sv => u(:,ke)
!      and                                       dvar(v_m)%sv => v(:,ke)
!      in order to force a "free-slip condition"!

  dvar(tem)%av  => t  ; dvar(tem)%at => ttens  ; dvar(tem)%sv => t_g    ; dvar(tem)%kstart = 1
  dvar(vap)%av  => qv ; dvar(vap)%at => qvtens ; dvar(vap)%sv => qv_s   ; dvar(vap)%kstart = 1
  dvar(liq)%av  => qc ; dvar(liq)%at => qctens ; dvar(liq)%sv => NULL() ; dvar(liq)%kstart = kstart_cloud

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

  IF ((lsfli(tem) .AND. .NOT.PRESENT(shfl_s)) .OR. &
      (lsfli(vap) .AND. .NOT.PRESENT(qvfl_s))) THEN
    CALL finish('', 'ERROR *** forcing with not present surface heat flux densities  ***')
  ENDIF

  IF (lsfli(tem)) dvar(tem)%sv => shfl_s
  IF (lsfli(vap)) dvar(vap)%sv => qvfl_s

  IF (PRESENT(ptr) .AND. ndtr .GE. 1) THEN !passive tracers are present
    DO m=1, ndtr
      n=liq+m
      dvar(n)%av => ptr(m)%av
      ltend(n)=ASSOCIATED(ptr(m)%at)
      IF (ltend(n)) THEN
        dvar(n)%at => ptr(m)%at
      ELSE
        dvar(n)%at => ptr(m)%av
      END IF
      IF (ASSOCIATED(ptr(m)%sv)) THEN
        dvar(n)%sv => ptr(m)%sv; lsfli(n)=ptr(m)%fc
      ELSE
        dvar(n)%sv => NULL()   ; lsfli(n)=.FALSE.
      END IF
      dvar(n)%kstart = ptr(m)%kstart
    END DO
  END IF

  vtyp(mom)%tkv => tkvm ; vtyp(mom)%tsv => tvm
  vtyp(sca)%tkv => tkvh ; vtyp(sca)%tsv => tvh

  !Note:
  !If a tendency field of an ordinary prognostic variable is not present,
  ! the related time step increment due to turbulent diffusion will be
  ! added to the prognostic variable directly.
  !It always holds: "lsfli(liq)=F"!

  fakt=z1/dt_var

  DO n=1,ndiff
    IF (ltend(n)) THEN  !calculation of tendencies required
      tinc(n)=z1        !no time increment multiplication for tendencies
    ELSE                !update of prognostic variables required
      tinc(n)=dt_var    !time increment multiplication for tendencies
    END IF
  END DO
  !$ACC UPDATE DEVICE(tinc) ASYNC(1)

!--------------------------------------------------
  IF (ldovardif .OR. ldogrdcor) THEN !Vertikaldiffusion wird hier berechnet
!--------------------------------------------------

my_cart_id = get_my_global_mpi_id()
#ifdef _OPENMP
my_thrd_id = omp_get_thread_num()
#endif

!########################################################################

         !Note: 
         !If ".NOT.ldovardif .AND. ldogrdcor", only a correction of pure vertical gradient diffusion
         ! due to sub grid scale condensation or non-local gradients is performed.
      
!        Berechnung der Luftdichte und des Exner-Faktors am Unterrand:
!DIR$ IVDEP
!$NEC ivdep
         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
         !$ACC LOOP GANG VECTOR PRIVATE(virt)
         DO i=ivstart, ivend
            virt=z1+rvd_m_o*qv_s(i) !virtueller Faktor
            rhon(i,ke1)=ps(i)/(r_d*virt*t_g(i))
            eprs(i,ke1)=zexner(ps(i))
         END DO
         !$ACC END PARALLEL
         !Note:
         !In the turbulence model 'rhon(:,ke1)' belongs to the lower boundary of the
         !Prandtl-layer, rather than to the surface level.
         !However, for the calculation in 'vert_grad_diff' only real surface-level values are used!

!        Setzen von Steuerparametern:

         mcorr=nmvar !end index for gradient correction belongs to the last dynamically active scalar by default
         IF (ldogrdcor) THEN !gradient correction for dynamically active prognostic variables:
            IF (lnonloc) THEN !in case of non-local gradient calculations:
               ncorr=1 !starting with the first variable in the list
            ELSE ! in case of local gradient calculations:
               ncorr=nvel+1 !only for scalar variables
               IF (.NOT.ldoexpcor) THEN !if just and only 'ldocirflx' is true:
                  mcorr=ncorr !only for the first scalar, which is (potential) temperature
               END IF
            END IF
         ELSE !no gradient correction at all
            ncorr=ndiff+1 !gradient correction must not start at any variable in the list
         END IF

         ivtype=0

!-----------------------------------------------------------------
!        Berechnung der Vertikaldiffusion von Modellvariablen auf Hauptflaechen:
!-----------------------------------------------------------------

!        DO n=nprim, nlast !loop over all variables to be diffused
         DO n=1, ndiff !loop over all variables to be diffused potentially

         ! define start index for vertical diffusion
         k_st_pp = dvar(n)%kstart

         IF ( (lum_dif .AND. n.EQ.u_m)   .OR. &                   !u_m-diffusion or
              (lvm_dif .AND. n.EQ.v_m)   .OR. &                   !v_m-diffusion or
              (lscadif .AND. n.GT.nvel)  .OR. &                   !sca-diffusion or
            (ldogrdcor .AND. n.GE.ncorr .AND. n.LE.mcorr) ) THEN  !gradient correction

            m=MIN(n,nmvar)

            IF (ivtype.EQ.0) THEN
               linisetup=.TRUE.; lnewvtype=.TRUE.
            ELSE
               linisetup=.FALSE.;lnewvtype=.FALSE.
            END IF

            IF (n.LE.nvel) THEN !a wind component
               IF (ivtype.EQ.sca) lnewvtype=.TRUE.

               !Attention:
               lsflucond=.FALSE. !no lower flux condition for momentum!

               ivtype=mom
            ELSE !a scalar property
               IF (ivtype.EQ.mom) lnewvtype=.TRUE.

               lsflucond=lsflcnd !use chosen type of lower boundary condition

               ivtype=sca
            END IF

            IF (n.LT.ncorr .OR. n.GT.mcorr) THEN !no gradient correction for this variable
               !Notice that this is particularly the case, if "ldogrdcor=F".
               igrdcon=0 !keine Gradientkorrektur der Profile
            ELSEIF ( (.NOT.lscadif .AND. ivtype.EQ.sca) .OR. &
                     (.NOT.lum_dif .AND.      n.EQ.u_m) .OR. &
                     (.NOT.lvm_dif .AND.      n.EQ.v_m) ) THEN !only a diffusion correction required
               igrdcon=1 !verwende nur Profil aus Gradientkorrektur
            ELSEIF (ldoexpcor .OR. lnonloc) THEN !full vertical diffusion of given non-gradient fluxes
               igrdcon=2 !verwende korrigiertes Profil aus effektiven Gradienten
            ELSE !full vertical diffusion including effective gradients of an extra non-local flux-contribution
               igrdcon=3 !addiere Gradientkorrektur zum vorhandenen Profil
            END IF

            IF (.NOT.ltend(n)) THEN !tendency array not present
               itndcon=0 !no explicit tendency consideration
            ELSE
               itndcon=itnd !use chosen mode of tendency consideration
            END IF
            IF (igrdcon.EQ.2) THEN !full vertical diffusion of given non-gradient fluxes
               k_st_up=ke !only level "k=ke" needs to be provided for bottom-up integration
            ELSE !vertical profiles needs to be provided 
               k_st_up=k_st_pp !up to the first present level
            END IF


            IF (lsfli(n)) THEN !surface value input is a flux density instead of a concentration
               !Load effective surface layer gradients due to given flux values:

               vtyp_tkv => vtyp(ivtype)%tkv
               dvar_sv  => dvar(n)%sv
!DIR$ IVDEP
!$NEC ivdep
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  wert=rhon(i,ke1)*vtyp_tkv(i,ke1)
                  IF (n.EQ.tem) THEN !flux density is that of sensible heat
                     zvari(i,ke1,m)=dvar_sv(i)/(wert*cp_d*eprs(i,ke1))
                  ELSE
                     zvari(i,ke1,m)=dvar_sv(i)/wert
                  END IF
               END DO
               !$ACC END PARALLEL

               !Attention:
               !In this case, not the current surface concentration, but the current flux-density at the surface,
               ! is being used in 'vert_grad_diff'!
               !Nevertheless, 'zvari' contains vertical gradients at this place; and for ".NOT.lsflucond", 
               ! a related surface concentration will be recalculated in 'vert_grad_diff' so that, in this case,
               ! an implicit deviation of the surface flux may still develop!
               !For the above calculation, 'tkv(ke1)' needs to be ">0", which is always the case, if it is 
               ! calculated by 'turbtran'; thus "tkvh(ke1)=0.0" should never be forced, if "lsfli=.TRUE"!
               !For tracers, it is always "m=nmvar"!
              
            END IF !surface value input is a flux density instead of a concentration

            ! Providing the profiles of concentrations and their current tendencies:

            cur_prof => hlp
            dvar_av  => dvar(n)%av    ! OpenACC issue with derived type

            !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR COLLAPSE(2)
            DO k=k_st_up,ke !necessary profile-levels (for bottom-up integration at "igrdcon.EQ.2")
!DIR$ IVDEP
!$NEC ivdep
               DO i=ivstart, ivend
                  IF (n.EQ.tem) THEN !temperature needs to be transformed
                     cur_prof(i,k)=dvar_av(i,k)/epr(i,k) !potential temperature
                     !Note: Here, '=dvar_av' points to ordinary temperature 't'.
                  ELSE
                     cur_prof(i,k)=dvar_av(i,k)
                  END IF
               END DO
            END DO
            !$ACC END PARALLEL

            !Surface concentrations:

            IF (ASSOCIATED(dvar(n)%sv)) THEN !surface variable is present
                 dvar_sv=>dvar(n)%sv !XL_CHANGE: OpenACC issue with derived type
!DIR$ IVDEP
!$NEC ivdep
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=dvar_sv(i)
               END DO
               !$ACC END PARALLEL
            ELSEIF (n.LE.nvel .OR. ilow_def_cond.EQ.2) THEN
               !No-slip-condition for momentum or zero-concentr. condition as a default:
!DIR$ IVDEP
!$NEC ivdep
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=z0
               END DO
               !$ACC END PARALLEL
            ELSE !(ilow_def_cond.EQ.1) 
                !Enforce a zero-flux condition as a default:
!DIR$ IVDEP
!$NEC ivdep
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=dvar_av(i,ke)
               END DO
               !$ACC END PARALLEL
            END IF
            IF (n.EQ.tem) THEN !temperature needs to be transformed
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  cur_prof(i,ke1)=cur_prof(i,ke1)/eprs(i,ke1) !potential surface temperature
                  !Note: Here, '=dvar_av' points to ordinary temperature 't'.
               END DO
               !$ACC END PARALLEL
            END IF
            !Current tendencies:

            IF (itndcon.GT.0) THEN !explicit tendencies have to be considered
               dvar_at => dvar(n)%at
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR COLLAPSE(2)
               DO k=k_st_pp,ke
!DIR$ IVDEP
!$NEC ivdep
                  DO i=ivstart, ivend
                     IF (n.EQ.tem) THEN !temperature needs to be transformed
                        dicke(i,k)=dvar_at(i,k)/eprs(i,ke1) !related to potential surface temperature
                     ELSE
                        dicke(i,k)=dvar_at(i,k)
                     END IF
                  END DO
               END DO
               !$ACC END PARALLEL
            END IF


!Achtung: Belegung von 'cur_prof' ev. verzichtbar bei "igrdcon=2"!

            leff_flux= .FALSE. !so far no extra vertical calculation of final fluxes necessary
            IF (.NOT.(lsfluse .AND. lsflcnd)) THEN ! calculation of effective flux density required
               IF ( ( n.EQ.tet .AND. PRESENT(shfl_s) ) .OR. ( n.EQ.vap .AND. PRESENT(qvfl_s) ) ) THEN
                  leff_flux = .TRUE.
               END IF
            END IF

            IF ( ( n.EQ.u_m .AND. PRESENT(umfl_s) ) .OR. ( n.EQ.v_m .AND. PRESENT(vmfl_s) ) ) THEN
               leff_flux = .TRUE.
               !Note:
               !Surface fluxes for momentum are always implicitly calculated.
            ENDIF

            IF (n.LE.nmvar) THEN !for any dynamically active variable
               leff_flux = l3dflxout !extra flux calculation, if requested
            END IF

!           Berechnung der vertikalen Diffusionstendenzen:
!XL_COMMENTS : this print seems to occurs for any debug level, on purpose ?
!            print*, ivtype, associated(vtyp(ivtype)%tkv)

            !$ACC WAIT(1)

            CALL vert_grad_diff( kcm,                                &
!
                 i_st=ivstart, i_en=ivend, k_tp=k_st_pp-1, k_sf=ke1, &
!
                 dt_var=dt_var, ivtype=ivtype, igrdcon=igrdcon, itndcon=itndcon, &
!
                 linisetup=linisetup, lnewvtype=lnewvtype,            &
                 lsflucond=lsflucond, lsfgrduse=lsfli(n),             &
                 ldynimpwt=ldynimp  , lprecondi=lprecnd,              &
                 leff_flux=leff_flux,                                 &
!
                 rho=rhoh, rho_n=rhon, hhl=hhl, r_air=r_air,          &
!
                 tkv=vtyp(ivtype)%tkv, tsv=vtyp(ivtype)%tsv,          &
!
                 impl_weight=impl_weight,                             &
!
                 disc_mom=zaux(:,:,1), expl_mom=zaux(:,:,2),          &
                 impl_mom=zaux(:,:,3), invs_mom=zaux(:,:,4),          &
                 diff_dep=zaux(:,:,5), diff_mom=len_scale,            &
                 invs_fac=frh, scal_fac=frm,                          &
!
                 dif_tend=dicke, cur_prof=cur_prof, eff_flux=zvari(:,:,m), &
                 lacc=.TRUE. )

            !Beachte:
            !'frh', 'frm' und 'len_scale' sind genauso wie 'zaux(:,:,1:5)' Hilfsspeicher in 'vert_grad_diff'.
            !Weil Fluesse ab "n>=liq=nmvar" nicht mehr benoetigt werden, bleibt 'zvari' nur bis 
            ! 'nmvar' dimensioniert und zvari(nmvar) wird auch fuer "n>nmvar" benutzt.

!           Sichern der Tendenzen:

            IF (n.EQ.tem) THEN
               dvar_at => dvar(n)%at
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR COLLAPSE(2)
               DO k=k_st_pp,ke
!DIR$ IVDEP
!$NEC ivdep
                  DO i=ivstart, ivend
                     dvar_at(i,k)=dvar_at(i,k)+epr(i,k)*dicke(i,k)*tinc(n)
                  END DO
               END DO
               !$ACC END PARALLEL
            ELSE
               dvar_at => dvar(n)%at
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR COLLAPSE(2)
               DO k=k_st_pp,ke
!DIR$ IVDEP
!$NEC ivdep
                  DO i=ivstart, ivend
                     dvar_at(i,k)=dvar_at(i,k)+dicke(i,k)*tinc(n)
                  END DO
               END DO
               !$ACC END PARALLEL
            END IF

            IF (n.EQ.vap .AND. PRESENT(qv_conv)) THEN
               !qv-flux-convergence (always a tendency) needs to be adapted:
               
               IF (lqvcrst) THEN !by initializing 'qv_conv' with vertical qv-diffusion
                  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
                  !$ACC LOOP GANG VECTOR COLLAPSE(2)
                  DO k=k_st_pp,ke
!DIR$ IVDEP
!$NEC ivdep
                     DO i=ivstart, ivend
                        qv_conv(i,k)=dicke(i,k)
                     END DO 
                  END DO
                  !$ACC END PARALLEL
               ELSE !by adding vertical qv-diffusion to 'qv_conv'
                  !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
                  !$ACC LOOP GANG VECTOR COLLAPSE(2)
                  DO k=k_st_pp,ke
!DIR$ IVDEP
!$NEC ivdep
                     DO i=ivstart, ivend
                        qv_conv(i,k)=qv_conv(i,k)+dicke(i,k)
                     END DO
                  END DO
                  !$ACC END PARALLEL
               END IF
            END IF   
                    
         END IF !diffusion calculation requested
         END DO !1, ndiff 

!-----------------------------------------------------------------

!Achtung:
!Ist cp-Fluss tatsaechlich der thermische Erdbodenantrieb?
!Was gilt im Falle der T-Gleichung in cv-Form?

!        Update of surface fluxes, if the vertical diffusion has determined them implicitly:

!Achtung: "lscadif" ergaenzt
         IF (.NOT.(lsfluse .AND. lsflcnd) .AND. lscadif) THEN 
            !effektive Oberfl.flussdichten wurden neu bestimmt

            IF (PRESENT(shfl_s) .OR. lrunscm) THEN
!DIR$ IVDEP
!$NEC ivdep
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  shfl_s(i)=eprs(i,ke1)*cp_d*zvari(i,ke1,tet)
               END DO
               !$ACC END PARALLEL
            END IF
            IF (PRESENT(qvfl_s) .OR. lrunscm) THEN
!DIR$ IVDEP
!$NEC ivdep
               !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
               !$ACC LOOP GANG VECTOR
               DO i=ivstart, ivend
                  qvfl_s(i)=zvari(i,ke1,vap)
               END DO
               !$ACC END PARALLEL
            END IF

!---------------------------------------------------------------------------------------
#ifdef SCLM
            IF (lsclm .AND. latmflu) THEN
               !Berechnung der Enthalpieflussdichten:

               SHF%mod(0)%val=shfl_s(imb)     ; SHF%mod(0)%vst=i_cal
               LHF%mod(0)%val=qvfl_s(imb)*lh_v; LHF%mod(0)%vst=i_cal

               !Note:
               !IF ".NOT.latmflu", SHF and LHF either are loaded by the fluxes used for
               ! the soil budget (lertflu) or they have been loaded above by the explicit 
               ! SHF and LHF at the surface (lsurflu).
               !SHF and LHF are positive downward and they may have been corrected with
               ! vertical integrated correction tendencies.
               !Thus they always refer to the used flux densities, which are only then equal
               ! to the explicit surface flux density, if a lower flux condition is used "lsflcnd=.TRUE.".
            END IF
#endif
!SCLM-----------------------------------------------------------------------------------

            !Bem: shfl_s und qvfl_s, sowie SHF und LHF sind positiv abwaerts!

         END IF

         IF (lum_dif .AND. PRESENT(umfl_s)) THEN
!DIR$ IVDEP
!$NEC ivdep
            !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO i=ivstart, ivend
               umfl_s(i)=zvari(i,ke1,u_m)
            END DO
            !$ACC END PARALLEL
         END IF
         IF (lvm_dif .AND. PRESENT(vmfl_s)) THEN
!DIR$ IVDEP
!$NEC ivdep
            !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT)
            !$ACC LOOP GANG VECTOR
            DO i=ivstart, ivend
               vmfl_s(i)=zvari(i,ke1,v_m)
            END DO
            !$ACC END PARALLEL
         END IF

         !Note:
         !The fluxes updated here are those effectively used for the atmospheric budgets and may slightly
         ! differ from the (aggregated) explicit surface fluxes used in the surface schemes!
         !The latent heat flux is not included here, since the required vaporization heat depends on the
         ! the surface state of of each tile.

!--------------------------------------------------
  END IF !Vertikaldiffusion wird hier berechnet
!--------------------------------------------------

  !$ACC WAIT(1)
  !$ACC END DATA

END SUBROUTINE vertdiff

!==============================================================================

END MODULE turb_vertdiff
