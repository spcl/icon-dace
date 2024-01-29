!
! mo_art_fplume_bpt
! This module solves the 1-D radially-averaged equations for a plume.
! Original code by A. Folch, G. Macedonio, A. Costa (2016)
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

MODULE mo_art_fplume_bpt 
  !*****************************************************************************
  !*
  !*    Module PlumeBPT : solves the 1-D radially-averaged equations for a plume
  !*
  !*****************************************************************************
  USE mo_kind,                          ONLY: wp
  USE mo_physical_constants,            ONLY: argas,amw,amd,t3,grav,ak,alv,alf,rhoh2o,rhoice
  USE mo_math_constants,                ONLY: pi
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH  
  USE mo_exception,                     ONLY: finish, message, message_text 
  USE mo_fplume_dlsode,                 ONLY: dlsode
  USE mo_art_fplume_types,              ONLY: t_fplume_init,t_fplume_phases

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initialize_plume_BPT, initialize_plume_wind, solve_plume_BPT
  SAVE
  INTEGER, PARAMETER :: s_work = 32  
  !
  !*** plume structure
  !
  TYPE volcanic_plume
    !
    LOGICAL               :: aggregation    
    LOGICAL               :: moist_air
    LOGICAL               :: wind_coupling 
    LOGICAL               :: reentrainment 
    LOGICAL               :: latent_heat   
    !
    CHARACTER(LEN=s_work) :: type_aggr
    CHARACTER(LEN=s_work) :: solve_for  ! HEIGHT or MFR
    CHARACTER(LEN=s_work) :: type_as    ! CONSTANT (value jet, value plume) / KAMINSKI-R / KAMINSKI-C / OLD
    CHARACTER(LEN=s_work) :: type_av    ! CONSTANT (value) / TATE
    !
    INTEGER               :: stat     ! Status code
    INTEGER               :: neq        ! number of equations for lsode
    INTEGER               :: ns         ! number of plume source points (plume+umbrella)
    INTEGER               :: np         ! number of plume source points (plume)
    INTEGER               :: nc         ! ++ number of particle classes
    INTEGER               :: modv       ! terminal velocity model
    !
    REAL(wp)              :: zv         ! vent altitude (m)
    REAL(wp)              :: n_MFR(2)   ! MER search range
    REAL(wp)              :: MER        ! MER (in kg/s) during this phase
    REAL(wp)              :: mass       ! Total erupted mass (in kg)
    REAL(wp)              :: xi         ! Factor (Bursik 2001).
    REAL(wp)              :: zmin_wind  ! Ignore wind entrainment below this zvalue (low jet region)
    REAL(wp)              :: c_umbrella ! Thickness of umbrella relative to Hb (>1)
    REAL(wp)              :: a_s_jet    ! Default (constant) value in jet   region
    REAL(wp)              :: a_s_plume  ! Default (constant) value in plume region
    REAL(wp)              :: a_v        ! Default (constant) value
    !
    REAL(wp), ALLOCATABLE :: x   (:)    ! lon(ns)
    REAL(wp), ALLOCATABLE :: y   (:)    ! lat(ns)
    REAL(wp), ALLOCATABLE :: z   (:)    ! z(ns) (elevation in m above terrain)
    REAL(wp), ALLOCATABLE :: Q   (:)    ! Bulk mass flow rate (ns)
    REAL(wp), ALLOCATABLE :: En  (:)    ! Total energy flow rate (ns)
    REAL(wp), ALLOCATABLE :: M   (:,:)  ! Particle mass flow rate (nc,ns)
    REAL(wp), ALLOCATABLE :: Mf  (:,:)  ! Mass flow rate of particles that fall from the eruption column (nc,ns)
    REAL(wp), ALLOCATABLE :: Magr(:,:)  ! Mass flow rate of particles that aggregate (nc,ns)
    REAL(wp), ALLOCATABLE :: Mair(:)    ! Mass flow rate of air (ns)
    REAL(wp), ALLOCATABLE :: Mw  (:)    ! Mass flow rate of volatiles (vapor + liquid + ice) (ns)
    REAL(wp), ALLOCATABLE :: L   (:)    ! Coordinate s (along the plume centerline)
    REAL(wp), ALLOCATABLE :: H   (:)    ! Theta angle
    REAL(wp), ALLOCATABLE :: U   (:)    ! Bulk velocity
    REAL(wp), ALLOCATABLE :: T   (:)    ! Bulk temperature
    REAL(wp), ALLOCATABLE :: D   (:)    ! Bulk density
    REAL(wp), ALLOCATABLE :: D_rhoa(:)  ! Bulk density/air density
    REAL(wp), ALLOCATABLE :: R   (:)    ! Plume radius
    REAL(wp), ALLOCATABLE :: xv  (:)    ! water vapor  mass fraction
    REAL(wp), ALLOCATABLE :: xl  (:)    ! liquid water mass fraction
    REAL(wp), ALLOCATABLE :: xs  (:)    ! ice (solid)  mass fraction
    REAL(wp), ALLOCATABLE :: as  (:)    ! a_shear
    REAL(wp), ALLOCATABLE :: av  (:)    ! a_vortex
    !
  END TYPE volcanic_plume
  TYPE(volcanic_plume) :: plume
  !
  !*** Numeric aspects
  !
  REAL(wp), ALLOCATABLE :: A_p(:),A_m(:)  ! aggregation coefficients
  REAL(wp), ALLOCATABLE :: fo (:),f  (:)  ! work arrays for BC at inlet and solutions
  REAL(wp)              :: s_old,r_old
  !
  REAL(wp) :: R0
  REAL(wp) :: u0
  REAL(wp) :: w0
  REAL(wp) :: P0
  REAL(wp) :: rho0
  REAL(wp) :: rhopart0
  REAL(wp) :: Enthalp0
  !
  !*** Plume status 
  !
  INTEGER,     PARAMETER :: STATUS_OK       = 0
  INTEGER,     PARAMETER :: STATUS_COLLAPSE = 1  ! Column collapse
  INTEGER,     PARAMETER :: STATUS_ERROR    = 2  ! Generic error
  !
  !*** Particles structure
  !
  TYPE particle_properties
    INTEGER               :: iaggr = 0 
    REAL(wp)              :: vset_aggr
    REAL(wp)              :: Dfo
    REAL(wp)              :: diam_aggr
    !
    REAL(wp), ALLOCATABLE :: fc    (:)                ! fc  (nc)
    REAL(wp), ALLOCATABLE :: rhop  (:)                ! rhop(nc)
    REAL(wp), ALLOCATABLE :: diam  (:)                ! diam(nc)
    REAL(wp), ALLOCATABLE :: sphe  (:)                ! sphe(nc)
    REAL(wp), ALLOCATABLE :: psi   (:)                ! psi (nc)
    REAL(wp), ALLOCATABLE :: vlimit(:)                ! vlimit(nc)
     !
  END TYPE particle_properties
  TYPE(particle_properties) :: part 
  !
  !*** Profile structure  
  !
  TYPE meteo_data
    INTEGER               :: nz
    REAL(wp), ALLOCATABLE :: z(:)       ! Height (a.s.l.)
    REAL(wp), ALLOCATABLE :: rho(:)     ! Air density
    REAL(wp), ALLOCATABLE :: p(:)       ! Pressure
    REAL(wp), ALLOCATABLE :: T(:)       ! Temperature
    REAL(wp), ALLOCATABLE :: sh(:)      ! Specific humidity
    REAL(wp), ALLOCATABLE :: u(:)       ! Wind speed
    REAL(wp), ALLOCATABLE :: PSIa(:)    ! Wind direction
    !
  END TYPE meteo_data
  TYPE(meteo_data)         :: profile
  !
  !*** Constants
  !
  REAL(wp) :: Cw   = 2000.0_wp      ! water (generic) used if latent_heat=.false.
  REAL(wp) :: Cp   = 1600.0_wp      ! solid particles (pyroclasts)
  REAL(wp) :: Ca   = 1000.0_wp      ! air
  !
  REAL(wp), PRIVATE, PARAMETER :: rvapour = argas/(amw*1.e-3_wp)    ! Gas constant of water vapour
  REAL(wp), PRIVATE, PARAMETER :: rair    = argas/(amd*1.e-3_wp)    ! Gas constant of air
  REAL(wp), PRIVATE, PARAMETER :: Pref    = 611.22_wp      ! Reference pressure (Triple point)
  REAL(wp), PRIVATE, PARAMETER :: Cv0     = 1996.0_wp      ! Specific heat of water vapour
  REAL(wp), PRIVATE, PARAMETER :: Cl0     = 4187.0_wp      ! Specific heat of water liquid
  REAL(wp), PRIVATE, PARAMETER :: Cs0     = 2108.0_wp      ! Specific heat of water solid (ice)
  !
  REAL(wp) :: Cv = Cv0  ! Set the default values
  REAL(wp) :: Cl = Cl0
  REAL(wp) :: Cs = Cs0
  !
  REAL(wp), PRIVATE, PARAMETER :: hs0  = t3*Cs0       ! Enthalpy of ice at T=Tref (it is set such that enthalpy_ice(0)=0)
  REAL(wp), PRIVATE            :: hl0  = hs0+alf       ! Enthalpy of liquid at T=Tref. Default value for latent_heat=.true.
  REAL(wp), PRIVATE            :: hv0  = hs0+alf+alv  ! Enthalpy of vapour at T=Tref. Default value for latent_heat=.true.
  !
  !
  !*** Internal function (keep private)
  !
  PRIVATE :: get_zone
  !
  !*** Interface for rhoa(z) and rhoa(z,T)
  !
  INTERFACE rhoa
    MODULE PROCEDURE rhoa_z
    MODULE PROCEDURE rhoa_z_T
  END INTERFACE rhoa
  !
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE initialize_plume_BPT(plume_ns,plume_np,fplume_init)
  !************************************************************************
  !*
  !*  Initialize the structures of the module and allocates memory
  !*
  !************************************************************************
  IMPLICIT NONE

  INTEGER , INTENT(IN) :: plume_ns               !< Number of integration steps up to top
  INTEGER , INTENT(IN) :: plume_np               !< Number of integration steps up to NBL
  TYPE(t_fplume_phases),INTENT(IN)  :: fplume_init
  INTEGER,     SAVE                 :: ipass = 0
  INTEGER                           :: ic
  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:initialize_plume_BPT:"  
  !
  !*** If necessary deallocates memory. This is done in order to call fplumeBPT
  !*** with different granulometries or solving strategies
  !
  IF (ipass == 0) THEN
    ipass = 1
  ELSE
    DEALLOCATE(part%fc)
    DEALLOCATE(part%diam)
    DEALLOCATE(part%rhop)
    DEALLOCATE(part%sphe)
    DEALLOCATE(part%psi)
    DEALLOCATE(part%vlimit)
    DEALLOCATE(A_p)
    DEALLOCATE(A_m)
    !
    DEALLOCATE(plume%x)
    DEALLOCATE(plume%y)
    DEALLOCATE(plume%z)
    DEALLOCATE(plume%L)
    DEALLOCATE(plume%H)
    DEALLOCATE(plume%U)
    DEALLOCATE(plume%T)
    DEALLOCATE(plume%D)
    DEALLOCATE(plume%D_rhoa)
    DEALLOCATE(plume%R)
    DEALLOCATE(plume%Q)
    DEALLOCATE(plume%En)
    DEALLOCATE(plume%Mw)
    DEALLOCATE(plume%Mair)
    DEALLOCATE(plume%xv)
    DEALLOCATE(plume%xl)
    DEALLOCATE(plume%xs)
    DEALLOCATE(plume%as)
    DEALLOCATE(plume%av)
    DEALLOCATE(plume%M)
    DEALLOCATE(plume%Mf)
    DEALLOCATE(plume%Magr)
    !
    DEALLOCATE(fo)
    DEALLOCATE(f )
  END IF
  !
  !*** Initializations and memory allocation 
  plume%nc = fplume_init%nc
  plume%ns = plume_ns          ! number of plume sources (plume+umbrella)
  plume%np = plume_np          ! number of plume sources (up to NBL)
  !
  !*** Number of equations for lsode. Note that for the COSTA aggragation model plume%nc additional
  !*** equations are necessary to determine the partition between fallen mass and aggregating mass
  IF (TRIM(fplume_init%plume_type_aggr)=='COSTA') THEN
    plume%neq = 9+2*plume%nc
  ELSE
    plume%neq = 9+plume%nc
  END IF
  !
  ALLOCATE(part%fc    (plume%nc))
  ALLOCATE(part%diam  (plume%nc))
  ALLOCATE(part%rhop  (plume%nc))
  ALLOCATE(part%sphe  (plume%nc))
  ALLOCATE(part%psi   (plume%nc))
  ALLOCATE(part%vlimit(plume%nc))
  ALLOCATE(A_p        (plume%nc))
  ALLOCATE(A_m        (plume%nc))
  A_p(:) = 0.0_wp
  A_m(:) = 0.0_wp
  !
  ALLOCATE(plume%x     (plume%ns))
  ALLOCATE(plume%y     (plume%ns))
  ALLOCATE(plume%z     (plume%ns))
  ALLOCATE(plume%L     (plume%ns))
  ALLOCATE(plume%H     (plume%ns))
  ALLOCATE(plume%U     (plume%ns))
  ALLOCATE(plume%T     (plume%ns))
  ALLOCATE(plume%D     (plume%ns))
  ALLOCATE(plume%D_rhoa(plume%ns))
  ALLOCATE(plume%R     (plume%ns))
  ALLOCATE(plume%Q     (plume%ns))
  ALLOCATE(plume%En    (plume%ns))
  ALLOCATE(plume%Mw    (plume%ns))
  ALLOCATE(plume%Mair  (plume%ns))
  ALLOCATE(plume%xv    (plume%ns))
  ALLOCATE(plume%xl    (plume%ns))
  ALLOCATE(plume%xs    (plume%ns))
  ALLOCATE(plume%as    (plume%ns))
  ALLOCATE(plume%av    (plume%ns))
  ALLOCATE(plume%M     (plume%nc,plume%ns))
  ALLOCATE(plume%Mf    (plume%nc,plume%ns))
  ALLOCATE(plume%Magr  (plume%nc,plume%ns))
  !
  ALLOCATE(fo(plume%neq))
  ALLOCATE(f (plume%neq))
  !
  !*** Initializations (fill the plume structure)
  plume%aggregation   = fplume_init%plume_aggregation
  plume%moist_air     = fplume_init%plume_moist_air 
  plume%wind_coupling = fplume_init%plume_wind_coupling
  plume%reentrainment = fplume_init%plume_reentrainment 
  plume%latent_heat   = fplume_init%plume_latent_heat 
  IF (.NOT. plume%latent_heat) CALL set_latent_heat_flag(.FALSE.)
  !
  plume%modv         = fplume_init%plume_modv
  plume%zv           = fplume_init%plume_zv
  plume%solve_for    = fplume_init%plume_solve_for
  plume%n_MFR(1:2)   = fplume_init%plume_n_MFR(1:2)
  plume%xi           = fplume_init%plume_xi
  plume%zmin_wind    = fplume_init%plume_zmin_wind
  plume%c_umbrella   = fplume_init%plume_c_umbrella
  plume%a_s_jet      = fplume_init%plume_a_s_jet
  plume%a_s_plume    = fplume_init%plume_a_s_plume
  plume%a_v          = fplume_init%plume_a_v
  !
  plume%type_aggr    = fplume_init%plume_type_aggr
  plume%type_as      = fplume_init%plume_type_as
  plume%type_av      = fplume_init%plume_type_av
  !
  !*** Initializations (fill the part structure)
  part%fc    (1:plume%nc) = fplume_init%fc  (1:plume%nc)
  part%diam  (1:plume%nc) = fplume_init%diam(1:plume%nc)
  part%rhop  (1:plume%nc) = fplume_init%rhop(1:plume%nc)
  part%sphe  (1:plume%nc) = fplume_init%sphe(1:plume%nc)
  part%psi   (1:plume%nc) = fplume_init%psi (1:plume%nc)
  part%vlimit(1:plume%nc) = 0.0_wp
  !
  part%vset_aggr     = fplume_init%vset_aggr
  part%diam_aggr     = fplume_init%diam_aggr
  part%Dfo           = fplume_init%Dfo
  IF (TRIM(plume%type_aggr)=='COSTA') THEN
    part%iaggr = 0
    DO ic = 1,plume%nc-1
      IF (part%diam_aggr>part%diam(ic)) THEN
        IF (part%iaggr==0) part%iaggr = ic    ! index of the first aggregating class
      END IF
    END DO
      IF (part%iaggr == 0) part%iaggr = plume%nc-1  ! ensure at least 1 aggregating class
  END IF
  !
  !*** Initializations (constants)
  Cw = fplume_init%Cw
  Cp = fplume_init%Cp
  Ca = fplume_init%Ca
  !
  RETURN
END SUBROUTINE initialize_plume_BPT
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE initialize_plume_wind(profile_nz,profile_z,profile_rho,profile_p,profile_T,  &
                                 profile_sh,profile_u,profile_PSIa) 
  !************************************************************************
  !*
  !*  Initialize the wind profile structure and allocates memory
  !*
  !************************************************************************
  IMPLICIT NONE
  INTEGER    , INTENT(IN)    :: profile_nz               !< Number of z-layers 
  REAL   (wp), INTENT(INOUT) :: profile_z   (profile_nz-1) !< Height of layers
  REAL   (wp), INTENT(IN)    :: profile_rho (profile_nz-1) !< Density
  REAL   (wp), INTENT(IN)    :: profile_p   (profile_nz-1) !< Pressure
  REAL   (wp), INTENT(IN)    :: profile_T   (profile_nz-1) !< Temperature
  REAL   (wp), INTENT(IN)    :: profile_sh  (profile_nz-1) !< Specific humidity
  REAL   (wp), INTENT(IN)    :: profile_u   (profile_nz-1) !< Wind speed
  REAL   (wp), INTENT(IN)    :: profile_PSIa(profile_nz-1) !< Wind direction 
  !
  INTEGER    , SAVE       :: ipass = 0
  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:initialize_plume_wind:" 
  !
  !*** If necessary deallocates memory. This is done in order to call fplumeBPT
  !*** with different wind profiles
  
  IF (ipass==0) THEN
    ipass = 1
  ELSE
    DEALLOCATE(profile%z   )
    DEALLOCATE(profile%rho )
    DEALLOCATE(profile%p   )
    DEALLOCATE(profile%T   )
    DEALLOCATE(profile%sh  )
    DEALLOCATE(profile%u   )
    DEALLOCATE(profile%PSIa)
  END IF
  !
  !*** Allocates and fills the structure
  profile%nz = profile_nz-1
  ALLOCATE(profile%z   (profile%nz))
  ALLOCATE(profile%rho (profile%nz))
  ALLOCATE(profile%p   (profile%nz))
  ALLOCATE(profile%T   (profile%nz))
  ALLOCATE(profile%sh  (profile%nz))
  ALLOCATE(profile%u   (profile%nz))  
  ALLOCATE(profile%PSIa(profile%nz))  
  !
  profile%z   (1:profile%nz) = profile_z   (1:profile%nz)
  profile%rho (1:profile%nz) = profile_rho (1:profile%nz)
  profile%p   (1:profile%nz) = profile_p   (1:profile%nz)
  profile%T   (1:profile%nz) = profile_T   (1:profile%nz)
  profile%sh  (1:profile%nz) = profile_sh  (1:profile%nz)
  profile%u   (1:profile%nz) = profile_u   (1:profile%nz)
  profile%PSIa(1:profile%nz) = profile_PSIa(1:profile%nz)
  !
  RETURN
END SUBROUTINE initialize_plume_wind
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE solve_plume_BPT(fplume_phase_value,plume_status,plume_radius,plume_MER,plume_x,plume_y,&
                         & plume_z,plume_Q,plume_En,plume_M,plume_Mf,plume_Magr,plume_Mair,       &
                         & plume_Mw,plume_L,plume_H,plume_U,plume_T,plume_D,plume_D_rhoa,plume_R, &
                         & plume_xv,plume_xl,plume_xs,plume_as,plume_av)
  !************************************************************************
  !*
  !*    Solves the 1-D radially-averaged equations for a plume.
  !*
  !*      np        Number of plume points for output. Points are equally spaced
  !*                along the s direction (i.e. not in the vertical z due
  !*                to plume bent-over by wind). The points start at s=0
  !*                and finish at the neutral bouyancy level.
  !*      ns        Total number of points for output (Plume+Umbrella). Mass
  !*                in the umbrella region is distributed following a Gaussian.
  !*      nc        Number of classes (including aggregates).
  !*                Plume equations are solved for the particles only.
  !*
  !*************************************************************************
  IMPLICIT NONE
  !
  TYPE(t_fplume_init), INTENT(IN) :: fplume_phase_value
  REAL(wp)    :: M0             !< MER at vent
  !
  INTEGER , INTENT(INOUT) :: plume_status                     !< Exit status
  REAL(wp), INTENT(INOUT) ::  &
      &  plume_radius,                   &  !< Vend radius
      &  plume_MER,                      &  !< MER
      &  plume_x     (plume%ns),         &  !< X-coordinates
      &  plume_y     (plume%ns),         &  !< Y-coordinates
      &  plume_z     (plume%ns),         &  !< Z-coordinates (m above terrain)
      &  plume_Q     (plume%ns),         &  !< Bulk mass flow rate
      &  plume_En    (plume%ns),         &  !< Total energy flow rate
      &  plume_M     (plume%nc,plume%ns),&  !< Particle mass flow rate
      &  plume_Mf    (plume%nc,plume%ns),&  !< MFR of particles that fall from the eruption column
      &  plume_Magr  (plume%nc,plume%ns),&  !< MFR of particles that aggregate
      &  plume_Mair  (plume%ns),         &  !< Mass flow rate of air
      &  plume_Mw    (plume%ns),         &  !< Mass flow rate of volatiles (vapor+liquid+ice)
      &  plume_L     (plume%ns),         &  !< Coordinate s (along the plume centerline)
      &  plume_H     (plume%ns),         &  !< Theta angle
      &  plume_U     (plume%ns),         &  !< Bulk velocity
      &  plume_T     (plume%ns),         &  !< Bulk temperature
      &  plume_D     (plume%ns),         &  !< Bulk density
      &  plume_D_rhoa(plume%ns),         &  !< Bulk density/air density
      &  plume_R     (plume%ns),         &  !< Plume radius
      &  plume_xv    (plume%ns),         &  !< water vapor mass fraction
      &  plume_xl    (plume%ns),         &  !< liquid water mass fraction
      &  plume_xs    (plume%ns),         &  !< ice (solid) mass fraction
      &  plume_as    (plume%ns),         &  !< a_shear
      &  plume_av    (plume%ns)             !< a_vortex
  !
  LOGICAL     :: lexit,go_on_s,go_on_MFR,jet_zone
  INTEGER     :: ind,indc
  INTEGER     :: MFR_iter,ic,istate,is
  REAL(wp)    :: n_MFR_min,n_MFR_max,s,so,ds,sb,sc,rhorhoa,rhorhoamin
  REAL(wp)    :: x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xa,xw,xv,xl,xs,xp
  REAL(wp)    :: rhoair,rhoaT,rhopart,Vaair,a_shear,a_vortex
  REAL(wp)    :: Hb,Hc,dz
  REAL(wp)    :: xnbl,ynbl,Tnbl,Qnbl,Ennbl,znbl,zpp,dxdz,dydz,dxnbl,dynbl,dznbl,time_u
  REAL(wp)    :: dMp,dEn,enthalpy,Pair,Pvap
  REAL(wp)    :: err
  !
  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:solve_plume_BPT:" 
  !
  !*** Initializations
  plume%Magr(:,:) = 0.0_wp
  plume%Mf  (:,:) = 0.0_wp
  plume%M   (:,:) = 0.0_wp
  !
  !*** Loop over MFR. It depends on the solving strategy (solve for MFR or height)
  go_on_MFR = .TRUE.
  MFR_iter  = 1
  n_MFR_min = plume%n_MFR(1)   ! min MFR (e.g. 10**2 . From input file)
  n_MFR_max = plume%n_MFR(2)   ! max MFR (e.g. 10**10. From input file)
  !
  DO WHILE(go_on_MFR)
    !
    plume%stat = STATUS_OK     ! Clear status flag
    !
    SELECT CASE(plume%solve_for)
      CASE('HEIGHT')
        M0 = fplume_phase_value%plume_Mdt
        CONTINUE                      ! MFR given, a single iteration with MFR=M0 is sufficient
      CASE('MFR')
        ! Method of bisection
        SELECT CASE(MFR_iter)
          CASE(1)
            M0 = 10.0_wp**n_MFR_min
          CASE(2)
            M0 = 10.0_wp**n_MFR_max
          CASE DEFAULT
            M0 = 10.0_wp**(0.5_wp*(n_MFR_min+n_MFR_max))
        END SELECT
      END SELECT
    !
    ! *** Set the initial conditions at the vent
    w0 = fplume_phase_value%plume_wvdt  &
     & + fplume_phase_value%plume_wldt  &
     & + fplume_phase_value%plume_wsdt                 ! Water mass fraction
    u0 = fplume_phase_value%plume_udt            ! Velocity at the vent (store for later use) 
    P0 = Pa(plume%zv)        ! Pressure at the vent
    Enthalp0 = (1._wp-w0)*enthalpy_particles(fplume_phase_value%plume_Tdt)  &  ! Enthalpy (per unit mass) at the vent
          &  + fplume_phase_value%plume_wvdt*enthalpy_vapour(fplume_phase_value%plume_Tvdt)         &  
          &  + fplume_phase_value%plume_wldt*enthalpy_liquid(fplume_phase_value%plume_Tldt)         &
          &  + fplume_phase_value%plume_wsdt*enthalpy_ice(fplume_phase_value%plume_Tsdt)
    rhopart0 = 1._wp/SUM(part%fc(1:plume%nc)/part%rhop(1:plume%nc))      
                                                  ! Average density of the particles at the vent
    rho0     = 1._wp/((1._wp-w0)/rhopart0           & ! Density of mixture at the vent
          &  + fplume_phase_value%plume_wvdt/density_vapour(fplume_phase_value%plume_Tvdt,P0)     & !(particles + vapour + liquid + ice)
          &  + fplume_phase_value%plume_wldt/rhoh2o+fplume_phase_value%plume_wsdt/rhoice)  
    R0 = sqrt(M0/(pi*rho0*U0))    ! Vent radius
    !
    !*** Load boundary conditions for this MFR iteration
    fo(1)  = M0         !  Eq(1).  TOTAL Mass flow rate   Q = pi*r*r*rho*u
    fo(2)  = M0*u0      !  Eq(2).  Axial momentum         P = pi*r*r*rho*u*u
    fo(3)  = 0.5_wp*pi  !  Eq(3).  Radial momentum        Theta
    fo(4)  = M0*(Enthalp0+grav*plume%zv+0.5_wp*u0*u0)     !  Eq(4).  Energy   J = pi*r*r*rho*u*E
    fo(5)  = 0._wp      !  Eq(5).  Mass flow rate of air         Ma
    fo(6)  = M0*w0      !  Eq(6).  Mass flow rate of volatiles   Mw = Q*xw
    fo(7)  = 0._wp      !  Eq(7).  Trajectory X relative to vent location
    fo(8)  = 0._wp      !  Eq(8).  Trajectory y relative to vent location
    fo(9)  = plume%zv   !  Eq(9).  Z  in m a.s.l.
    DO ic = 1,plume%nc
      fo(9+ic) = M0*(1._wp-w0)*part%fc(ic)   !  Eq(9+ic). Mass flow rate of particles Mi
    END DO
    !
    IF (plume%type_aggr=='COSTA') THEN
      DO ic = 1,plume%nc
        fo(9+plume%nc+ic) = 0.0_wp                !  Eq(9+plume%nc+ic). Mass flow rate of particles that aggregate
      END DO
    END IF
    !
    !***  Loop over s to find Sb, the Neutral Buoyancy Level in s coordinate
    go_on_s = .TRUE.
    s       =  0.0_wp
    so      =  0.0_wp
    ds      = 10.0_wp   ! In meters: arbitrary choice (change with caution)
    f(:)    = fo(:)     ! Set initial conditions
    !
    rhorhoamin = 1.0e9_wp
    jet_zone   = .TRUE.
    !
    s_old   = so
    r_old   = 0.0_wp
    !
    ind = 0
    DO WHILE(go_on_s)
      s  = s + ds
      CALL splum(so,f,s,plume%neq,istate)      ! integration from so to s
      ind = ind + 1                            ! Number of iteration
      !
      IF (istate == 2) THEN
        IF (plume%stat /= STATUS_OK) THEN
        ! Plume collapse/errors detected
          go_on_s = .FALSE.
          lexit   = .FALSE.
         ELSE
         ! No errors detected
           CALL getplume(f,x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xv,xl,xs,xp)  !s
           rhoair     = rhoa(z)
           rhorhoa    = rho/rhoair
           rhorhoamin = min(rhorhoamin,rhorhoa)
           !
           sb = s   ! Save reached plume length
           ! Check jet/convective transition
           IF (rhorhoa <= 1.0_wp .AND. jet_zone) jet_zone = .FALSE.
           ! Check if we have reached the NBL
           IF (rhorhoa > 1.0_wp .AND.(.NOT.jet_zone)) THEN
             go_on_s = .FALSE.
             lexit   = .TRUE.
           END IF
           !
           !*** Update iteration values
           s_old   = s
           r_old   = r
         END IF
         !
      ELSE
      !*** Wrong termination (istate /= 2)
        go_on_s = .FALSE.
        lexit   = .FALSE.
      END IF
      ! 
      !*** Handle collapse (for any return status)
      IF (is_collapse(plume)) THEN
        go_on_s = .FALSE.
        lexit   = .FALSE.
        indc = ind-1  ! Save index before collapse
        sc = s - ds   ! Save plume length before collapse
      END IF
      !
    END DO ! go_on_s
    !
    SELECT CASE(plume%solve_for)
      CASE('HEIGHT')
        go_on_MFR = .FALSE.   ! Do not need other iterations
        IF (.NOT.lexit) THEN
          IF (is_collapse(plume)) THEN
            WRITE(message_text,*) 'plume collapse detected'
            CALL message(thisroutine,message_text)
          ELSE
            CALL finish(thisroutine,'wrong termination')
          END IF
        END IF
        !
      CASE('MFR')
      !
        IF (lexit) THEN
        !
        !***  Checks if the plume height is close to the required value.
        !***  Note that z is the NBL in m a.s.l. and HPlume is the total
        !***  required height (including umbrella) in m a.v.
          err = fplume_phase_value%plume_Hdt - plume%c_umbrella*(z - plume%zv + 8d0*R0)
          !
          IF (ABS(err) <= ds) THEN
            go_on_MFR = .FALSE.   ! Done
          ELSE IF (err < 0.0_wp .OR. is_collapse(plume)) THEN
            n_MFR_max = LOG10(M0)
          ELSE
            n_MFR_min = LOG10(M0)
          END IF
             !
        ELSE
        ! Go here if column collapse or wrong termination of splum
        ! n_MFR_max = 0.5_wp*(n_MFR_min+n_MFR_max)
          n_MFR_max = LOG10(M0)
        END IF
        !
        MFR_iter = MFR_iter + 1
        !
      IF (MFR_iter == 100) THEN
        IF (ABS(err) <= 1000.0_wp) THEN
          go_on_MFR = .FALSE.
        ELSE
          WRITE(message_text,*) 'MFR iterations = Itermax. Convergence not achieved:', err 
          CALL finish(thisroutine,message_text)          
        ENDIF
      ENDIF
        !
    END SELECT
       !
  END DO   
  !
  !*** Integrates again up to the NBL (sb) and store the variables
  !*** in the plume%np points
  s  = 0.0_wp
  so = 0.0_wp
  ds = sb/plume%np                          ! np points
  f(:) = fo(:)                              ! Set initial conditions
  plume%MER = M0                            ! Store mass eruption rate
  !
  s_old   = so
  r_old   = 0.0_wp
  !
  ! COLLAPSING COLUMN: Set maximum path along the column
  IF (is_collapse(plume)) THEN
    plume%np = MIN(indc,plume%ns) ! Use maximum memory space
    ds = sc/dble(plume%np+1)      ! Add 1 to compensate numerical errors
  END IF
  !
  DO is = 1,plume%np
    s  = s + ds
    CALL splum(so,f,s,plume%neq,istate)
    IF (istate /= 2) CALL finish(thisroutine,'Error in the second loop?')

    CALL getplume(f,x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xv,xl,xs,xp)
    !
    IF (plume%wind_coupling) THEN
      Vaair  = Va(z)
      rhoair = rhoa(z)
    ELSE
      Vaair  = 0.0_wp
      rhoair = rhoa(z)
    END IF
    !
    !*** call get_entrainment_coef for storing a_shear and a_vortex
    CALL get_entrainment_coef(z-plume%zv,R0,U0,rhoair,rho,r,u,Vaair,th,a_shear,a_vortex)
    !
    !*** Update iteration values
    s_old   = s
    r_old   = r
    !
    !*** Store values to the plume structure
    CALL setplume(is,f,x,y,z,s,th,u,T,Q,En,r,rho,Ma,Mw,xv,xl,xs,a_shear,a_vortex)
    !
  END DO
  !
  ! COLLAPSING COLUMN: return
  !
  IF (is_collapse(plume)) THEN  
    CALL finish(thisroutine,'COLLAPSE of volcanic plume')
  ELSE
    !
    !*** UMBRELLA REGION (empirical model)
    !
    !*** Variables stored from M(ic,np) to M(ic,ns), i.e. on ns-np+1 points.
    !*** The total column height Hc is assumed to be Hc=c_umbrella*(Hb+8R0)
    !*** Gaussian profile with radius e^(-2) ar z=Hc

    Hb = plume%z(plume%np)-plume%zv       ! NBL   height (above vent)
    Hc = plume%c_umbrella*(Hb+8.0_wp*R0)     ! Total height (above vent)
    dz = (Hc-Hb)/(plume%ns-plume%np)      ! ns-np+1 points
    !
    xnbl = plume%x(plume%np)
    ynbl = plume%y(plume%np)
    znbl = plume%z(plume%np)
    Qnbl = plume%Q(plume%np)
    Tnbl = plume%T(plume%np)
    Ennbl= plume%En(plume%np)                     ! Energy flow rate at NBL
    dxnbl=(plume%x(plume%np)-plume%x(plume%np-1)) ! DX at NBL
    dynbl=(plume%y(plume%np)-plume%y(plume%np-1)) ! DY at NBL
    dznbl=(plume%z(plume%np)-plume%z(plume%np-1)) ! DZ at NBL
    dxdz= dxnbl/dznbl                             ! dx/dz at NBL
    dydz= dynbl/dznbl                             ! dy/dz at NBL
    !
    DO is = 1,plume%ns-plume%np
      !
      !*** Coordinates: linear extrapolation
      plume%x(plume%np+is) = xnbl + dxdz*is*dz
      plume%y(plume%np+is) = ynbl + dydz*is*dz
      plume%z(plume%np+is) = plume%z(plume%np)+(is)*dz
      ! Normalized height above NBL (zpp=0 at Hb, zpp=1 at Hc)
      zpp = dble(is)/dble(plume%ns-plume%np)
      !
      !*** Length. L = L + dL where dL^2=dx^2+dy^2+dz^2
      plume%L(plume%np+is) = plume%L(plume%np+is-1) + dz*sqrt( dxdz**2   &
                         & + dydz**2 + 1._wp)
      !
      !*** Radius. Decreases with zpp as a Gaussian with sigma=1/2
      plume%R(plume%np+is) = plume%R(plume%np)*exp(-2._wp*zpp**2)
      !
      !*** Theta = theta(NBL)
      plume%H(plume%np+is) = plume%H(plume%np)
      !
      !*** Velocity: U^2 = Unbl^2 - k*z
      plume%u(plume%np+is) = plume%u(plume%np)*sqrt(1._wp-zpp)
      !
      !*** Mair and volatiles (constant)
      plume%Mair(plume%np+is) = plume%Mair(plume%np)
      plume%Mw  (plume%np+is) = plume%Mw  (plume%np)
      !
      !*** Particle vertical mass flow rate (decrease exponentially)
      DO ic = 1,plume%nc
        plume%M(ic,plume%np+is) = plume%M(ic,plume%np)*exp(-2._wp*zpp**2)
      END DO
      !
      Mp = SUM(plume%M(1:plume%nc,plume%np+is))
      Ma = plume%Mair(plume%np+is)
      Mw = plume%Mw(plume%np+is)
      plume%Q(plume%np+is) = Mp+Mw+Ma                   ! Total mass flow rate
      xp = Mp/(Mp+Mw+Ma)                                ! particle  mass fraction
      xw = Mw/(Mp+Mw+Ma)                                ! volatiles mass fraction
      xa = Ma/(Mp+Mw+Ma)                                ! air mass fraction
      dMp = Mp - SUM(plume%M(1:plume%nc,plume%np+is-1)) ! Loss of particles
      dEn = dMp*(enthalpy_particles(plume%T(plume%np+is-1))  &
        & + 0.5_wp*plume%u(plume%np+is-1)**2 + grav*plume%z(plume%np+is-1))
      plume%En(plume%np+is) = plume%En(plume%np+is-1) + dEn
      ! Enthalpy per unit mass
      enthalpy = plume%En(plume%np+is)/(Mp+Mw+Ma)-grav*plume%z(plume%np+is) &
             & - 0.5_wp*plume%u(plume%np+is)**2
      !
      CALL temperature_mixture(xa,xp,xw,plume%T(plume%np+is), &
         & Pa(plume%u(plume%np+is)),plume%xv(plume%np+is), &
         & plume%xl(plume%np+is),plume%xs(plume%np+is),Pair,Pvap,enthalpy)
      !
      !*** Density
      ! Density of air at (z,T_plume)
      rhoaT   = rhoa(plume%z(plume%np+is),plume%T(plume%np+is))
      !
      ! Average particle density
      rhopart = Mp/SUM(plume%M(1:plume%nc,plume%np+is)/part%rhop(1:plume%nc))
      !
      ! Bulk density
      plume%D(plume%np+is)=1._wp/(xp/rhopart+xl/rhoh2o+xs/rhoice+(1._wp-xp-xl-xs)/rhoaT)
      plume%D_rhoa(plume%np+is) = plume%D(plume%np+is)/rhoa(plume%z(plume%np+is))
      !
      !*** Entrainment
      plume%as(plume%np+is) = plume%as(plume%np)
      plume%av(plume%np+is) = plume%av(plume%np)
    !
    ENDDO
    !
    !*** If necessary, modify mass of particles to account for aggregation in the umbrella region
    IF (TRIM(plume%type_aggr)=='COSTA') THEN
      !
      time_u = 2._wp  ! time factor for the umbrella region (up and down)
      !
      DO is = 1,plume%ns-plume%np
        f(1:plume%nc) = plume%M(1:plume%nc,plume%np+is)                                             ! Particle MFR
        Q = SUM(plume%M(1:plume%nc,plume%np+is)) + plume%Mair(plume%np+is) + plume%Mw(plume%np+is)  ! Total MFR
        CALL costa(f, plume%z (plume%np+is),plume%u (plume%np+is),plume%R(plume%np+is), &
                &  plume%T (plume%np+is),plume%D (plume%np+is),Q, &
                &  plume%xv(plume%np+is),plume%xl(plume%np+is),plume%xs(plume%np+is))
        !
        DO ic = part%iaggr,plume%nc
          IF (ic<plume%nc) THEN
            Mp = time_u*(A_p(ic)-A_m(ic))*dz
            Mp = MIN(Mp,plume%M(ic,plume%np+is))  ! limit the amount of aggregates
          ELSE
            Mp = SUM(plume%Magr(1:plume%nc-1,plume%np+is))
            Mp = -Mp
          ENDIF
          !
          plume%Magr(ic,plume%np+is) = Mp
          plume%M   (ic,plume%np+is) = plume%M(ic,plume%np+is) + Mp
        ENDDO
      ENDDO
    ENDIF
    !
    !*** Finally, computes the mass that FALLS from the column as:
    !*** Mf = DM/Ds - A(+) + A(-)
    !
    DO ic = 1,plume%nc
      plume%Mf(ic,plume%ns) = plume%M(ic,plume%ns-1)
    ENDDO
    !
    DO ic = 1,plume%nc
      !
      ! Umbrella region (Magr already substracted)
      DO is = plume%ns-1,plume%np+1,-1
        Mp =  plume%M(ic,is-1) - plume%M(ic,is)
        plume%Mf(ic,is) = MAX(Mp,0.0_wp)
      ENDDO
      !
      ! plume region. Note that plume%Magr(ic,is-1)-plume%Magr(ic,is) approaches the aggragated
      ! mass in the slab only. Small mass imbalance (<1%) can ocurr in case of aggregation because of this...
      !
      DO is = plume%np,2,-1
        Mp =  plume%M(ic,is-1) - plume%M(ic,is) - (plume%Magr(ic,is-1)-plume%Magr(ic,is))
        plume%Mf(ic,is) = MAX(Mp,0.0_wp)
      ENDDO
    ENDDO
    !
    DO ic = 1,plume%nc
      Mp = M0*(1.0_wp-w0)*part%fc(ic) - plume%M(ic,1) + plume%Magr(ic,1)
      plume%Mf(ic,1)= MAX(Mp,0.0_wp)
    ENDDO
  !
  !*** Finally, load values to return
  !
    plume_status    = plume%stat
    plume_radius    = R0
    plume_MER       = plume%MER
    !
    plume_x(:)      = plume%x(:)
    plume_y(:)      = plume%y(:)
    plume_z(:)      = plume%z(:)
    plume_Q(:)      = plume%Q(:)
    plume_En(:)     = plume%En(:)
    plume_M(:,:)    = plume%M(:,:)
    plume_Mf(:,:)   = plume%Mf(:,:)
    plume_Magr(:,:) = plume%Magr(:,:)
    plume_Mair(:)   = plume%Mair(:)  
    plume_Mw(:)     = plume%Mw(:)
    plume_L(:)      = plume%L(:)
    plume_H(:)      = plume%H(:)
    plume_U(:)      = plume%U(:)
    plume_T(:)      = plume%T(:)
    plume_D(:)      = plume%D(:)
    plume_D_rhoa(:) = plume%D_rhoa(:)
    plume_R(:)      = plume%R(:)  
    plume_xv(:)     = plume%xv(:)
    plume_xl(:)     = plume%xl(:)
    plume_xs(:)     = plume%xs(:)
    plume_as(:)     = plume%as(:)
    plume_av(:)     = plume%av(:)
    !
    RETURN
  ENDIF
END SUBROUTINE solve_plume_BPT
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE getplume(f,x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xv,xl,xs,xp)
  !********************************************************************
  !*
  !*    Extracts physical variables from values in the ODE's at
  !*    a certain value of the arc parameter s.
  !*
  !*    OUTPUTS : f,x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,xv,xl,xs,xp
  !*
  !********************************************************************
  IMPLICIT NONE
  !
  REAL(wp), INTENT(INOUT) :: f(plume%neq)
  REAL(wp), INTENT(INOUT) :: x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xv,xl,xs,xp
  !
  INTEGER     :: ic
  REAL(wp)    :: P,Pv,Pair,xw,xa,rhopart,rhoaT,enthalpy
  !
  !*** Limit mass of particles (in case of aggregation a sink term appears and
  !*** mass of aggregating classes could become negative if not limited)
  !
  IF (TRIM(plume%type_aggr)=='COSTA') THEN
    DO ic = part%iaggr,plume%nc-1
      f(9+ic) = MAX(f(9+ic),0.0_wp)
    END DO
  END IF
  !
  !***  Extract values
  !
  Q  = f(1)                     ! Total mass flow rate
  u  = f(2)/f(1)                ! Bulk velocity
  th = f(3)                     ! Plume bent-over angle (in rad).
  En = f(4)                     ! Energy flow rate
  Ma = f(5)                     ! Mass flux of entrained air
  Mw = f(6)                     ! Mass flux of volatiles (vapor+liquid+ice)
  x  = f(7)                     ! Coordinate-X
  y  = f(8)                     ! Coordinate-Y
  z  = f(9)                     ! Height above sea level
  Mp = SUM(f(9+1:9+plume%nc)) ! Mass flux of particles
  !
  xp = Mp/(Mp+Mw+Ma)            ! particle  mass fraction
  xw = Mw/(Mp+Mw+Ma)            ! volatiles mass fraction
  xa = Ma/(Mp+Mw+Ma)            ! air mass fraction
  !
  enthalpy = f(4)/f(1)-grav*z-0.5_wp*u*u ! Enthalpy per unit mass
  !
  ! Total pressure (P of external air)
  P  = Pa(z)
  ! Evaluate bulk temperature and volatile fractions
  CALL temperature_mixture(xa,xp,xw,T,P,xv,xl,xs,Pair,Pv,enthalpy)
  !
  !*** Air density at P(z),T
  !
  rhoaT = rhoa(z,T)
  !
  !*** Bulk density and plume radius
  !
  ! Average particle density
  rhopart = Mp/SUM(f(9+1:9+plume%nc)/part%rhop(1:plume%nc))
  !
  ! Bulk density
  rho = 1._wp/(xp/rhopart + xl/rhoh2o + xs/rhoice + xa/rhoaT + xv/density_vapour(T,P))
  !
  r  = sqrt(Q/(pi*rho*u))   ! Plume radius
  !
  RETURN
END SUBROUTINE getplume
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE setplume(is,f,x,y,z,s,th,u,T,Q,En,r,rho,Ma,Mw,xv,xl,xs,a_shear,a_vortex)
  !
  !*** Store values in the plume structure
  !
  IMPLICIT NONE
  REAL(wp),    INTENT(IN) :: x,y,z,s,th,u,T,Q,En,r,rho,Ma,Mw,xv,xl,xs,a_shear,a_vortex
  REAL(wp),    INTENT(IN) :: f(plume%neq)
  INTEGER,     INTENT(IN) :: is
  !
  plume%x     (is) = x
  plume%y     (is) = y
  plume%z     (is) = z
  plume%L     (is) = s
  plume%H     (is) = th
  plume%u     (is) = u
  plume%T     (is) = T
  plume%Q     (is) = Q
  plume%En    (is) = En
  plume%R     (is) = r
  plume%D     (is) = rho
  plume%D_rhoa(is) = rho/rhoa(z)
  plume%Mair  (is) = Ma
  plume%Mw    (is) = Mw
  plume%xv    (is) = xv
  plume%xl    (is) = xl
  plume%xs    (is) = xs
  plume%as    (is) = a_shear
  plume%av    (is) = a_vortex
  plume%M(1:plume%nc,is) = f(9+1:9+plume%nc)
  IF (plume%type_aggr=='COSTA') plume%Magr(1:plume%nc,is) = f(9+plume%nc+1:9+plume%nc+plume%nc)
  !
  RETURN
END SUBROUTINE setplume
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE splum(so,f,s,neq,istate)
  !************************************************************************
  !*
  !*    Returns the value of f(neq,s ) (at point s ) given the
  !*    initial conditions   f(neq,so) (at point so)
  !*
  !************************************************************************
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)    ::  neq
  INTEGER, INTENT(INOUT) ::  istate
  REAL(wp),INTENT(INOUT) ::  f(neq)
  REAL(wp),INTENT(INOUT) ::  so,s
  !
  INTEGER  ::  iwork(20)
  REAL(wp), ALLOCATABLE :: rwork(:)
  INTEGER  ::  nneq(1),itol,itask,iopt,lrw,liw,mf
  REAL(wp) ::  rtol(1),atol(1)
  !
  !***  Defines lsode parameters
  !
  nneq(1) = neq
  itol   = 1
  itask  = 1              ! Normal computation (overshooting)
  istate = 1
  iopt   = 0
  lrw    = 20+16*neq      ! Dimension of rwork array
  liw    = 20
  mf     = 10             ! Jacobian is not supplied
  rtol(1) = 1.e-6_wp          ! Relative error
  atol(1) = 1.e-6_wp          ! Absolute error
  !
  ! Allocate work array (Allocation/deallocation should be moved out of splum)
  !
  ALLOCATE(rwork(lrw))
  IF (itask==2) rwork(1)=s  ! Set TCRIT = TOUT (no overshooting)
  CALL dlsode(dfds,nneq,f,so,s,itol,rtol,atol,itask, &
     & istate,iopt,rwork,lrw,iwork,liw,jac,mf)
  !
  ! Deallocate work array
  !
  DEALLOCATE(rwork)
  !
  RETURN
END SUBROUTINE splum
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE dfds(neq,s,f,E)
  !*************************************************************************
  !*
  !*    Defines the system of equations for lsode loading E(i)
  !*
  !*    df
  !*    -- = E(s,f)    i = 1:neq
  !*    ds
  !*
  !*   f(1)  = Total mass flow rate
  !*   f(2)  = Axial momentum flow rate
  !*   f(3)  = Radial momentum flow rate (theta angle)
  !*   f(4)  = Total energy (thermal + potential + kinetic)
  !*   f(5)  = Mass flow rate of entrained air
  !*   f(6)  = Mass flow rate of volatiles
  !*   f(7)  = X
  !*   f(8)  = Y
  !*   f(9)  = Z
  !*   f(9+1:9+nc)       = Mass flow rate for each particle class
  !*   f(9+nc+1:9+nc+nc) = Mass flow rate of aggregating particles
  !
  !*************************************************************************
  IMPLICIT NONE
  !
  INTEGER     :: neq
  REAL(wp)    :: s,f(neq),E(neq)
  !
  LOGICAL     :: conve
  INTEGER     :: ic
  REAL(wp)    :: x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xv,xl,xs,xp    ! Derivated variables
  REAL(wp)    :: Taair,rhoair,muair,rhoaT,wair,Vaair,PSIair,hair,ue,dvaadz
  REAL(wp)    :: a_shear,a_vortex
  REAL(wp)    :: dsdr,Po,Fo,fre,RHS,xw

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:dfds" 
  !
  !*** Current values of derived variables
  !
  CALL getplume(f,x,y,z,th,u,T,r,rho,Ma,Mw,Mp,Q,En,xv,xl,xs,xp)
  !
  ! Check for collapsing column
  IF (f(2) < 0._wp .OR. f(3) < 0._wp) THEN
    plume%stat = STATUS_COLLAPSE
    E = 0._wp   ! Set all derivatives to 0 and return
    RETURN
  END IF
  ! Check for errors
  IF (T < 0._wp) THEN
    plume%stat = STATUS_ERROR
    E = 0._wp   ! Set derivatives to 0 and return
    RETURN
  END IF
  !
  Taair  = Ta   (z)                   !  Temperature of the entrained air at z
  rhoair = rhoa (z)                   !  Density     of the entrained air at z
  muair  = mua  (T)                   !  Air viscosity (at bulk temperature)
  rhoaT  = rhoa (z,T)                 !  Air density   (at bulk temperature)
  !
  !*** Compute the settling velocity (using air density and viscosity at bulk temperature)
  !
  DO ic = 1,plume%nc

    CALL vsettl(part%diam(ic),part%rhop(ic),rhoaT,muair,part%vlimit(ic),plume%modv,part%psi(ic),conve)

    IF (.NOT.conve) THEN
      CALL finish(thisroutine,'No convergence in vsettl routine')
    END IF
  END DO
  IF (plume%aggregation) part%vlimit(plume%nc) = part%vset_aggr*part%vlimit(plume%nc)  ! modify settling velocity of aggregating class
  !
  !***  Relative humidity of entrained air at z (kg/kg)
  !
  IF (plume%moist_air) THEN
    wair = sh_a (z)
  ELSE
    wair = 0.0_wp
  END IF
  !
  IF (plume%wind_coupling) THEN
    Vaair  = Va   (z)
    PSIair = Da   (z)
    dVaadz = dVadz(z)
  ELSE
    Vaair  = 0.0_wp
    PSIair = 0.0_wp
    dVaadz = 0.0_wp
  END IF
  !
  !*** Computes aggregation coefficients A_p(nc) and A_m(nc)
  !*** Note that this is necessary only for the Costa model because in the
  !*** other aggregation models aggregates are "formed" at the vent
  !
  IF (TRIM(plume%type_aggr).EQ.'COSTA') CALL costa(f(9+1),z,u,r,T,rho,Q,xv,xl,xs)
  !
  !*** Computes the entrainment velocity
  !
  CALL get_entrainment_coef(z-plume%zv,R0,U0,rhoair,rho,r,u,Vaair,th,a_shear,a_vortex)
  !
  IF ((z-plume%zv)<plume%zmin_wind) THEN
    ue = 0.0_wp                      ! Ignore entrainment in the low jet zone
  ELSE
    ue = a_shear*ABS(u-Vaair*COS(th)) + a_vortex*ABS(Vaair*SIN(th))
  END IF
  !
  !*** Defines the equations (note that the order in the array E needs to
  !*** be changed to account for dependencies)
  !
  !  E(9+1   :9+nc   ): Mass flux for each particle class
  !  E(9+nc+1:9+nc+nc): Mass flux of aggregating particles
  !
  ! These parameters are used with reentrainment

  Po  = pi*R0*R0*u0*u0
  Fo  = pi*R0*R0*u0*Enthalp0   ! Enthalp0 = Cbo*T0
  !
  DO ic = 1,plume%nc
    !
    IF (plume%reentrainment) THEN
      IF ( ABS((r-r_old))>1.e-6_wp) THEN
        dsdr = (s-s_old)/(r-r_old)
      ELSE
        dsdr = 0._wp
      END IF
        fre = 0.43_wp/(1._wp+(0.78_wp*part%vlimit(ic)*(Po**0.25_wp)/sqrt(Fo))**6._wp)
    ELSE
      fre  = 0._wp
    END IF
    !
    RHS = 1._wp/(1._wp+fre*ue*dsdr/part%vlimit(ic)) ! Modified reentraimment factor
    E(9+ic) = -plume%xi*f(9+ic)*part%vlimit(ic)*RHS/(u*r) + A_p(ic) - A_m(ic)
    ! Used only for fallen mass later on

    IF (TRIM(plume%type_aggr)=='COSTA') E(9+plume%nc+ic) = A_p(ic) - A_m(ic)
  END DO
  !
  ! E(1): Total mass flow rate
  !
  E(1) = 2.0_wp*pi*r*rhoair*ue + SUM(E(9+1:9+plume%nc))
  !
  ! E(2): Axial momentum flow rate
  !
  E(2) = pi*r*r*(rhoair-rho)*grav*SIN(th)+Vaair*COS(th)*2.0_wp*pi*r*rhoair*ue &
     & + u*SUM(E(9+1:9+plume%nc))
  !
  ! E(3): Radial momentum flow rate (theta)
  !
  E(3) = (pi*r*r*(rhoair-rho)*grav*COS(th)-Vaair*SIN(th)*2.0_wp*pi*r*rhoair*ue)/(pi*r*r*rho*u*u)
  !
  ! E(5): Mass flow rate of entrained (dry) air
  !
  E(5) = 2.0_wp*pi*r*rhoair*ue*(1.0_wp-wair)
  !
  ! E(6): Mass flow rate of volatiles
  !
  E(6) = 2.0_wp*pi*r*rhoair*ue*wair
  !
  ! E(7) : X
  !
  E(7) = COS(th)*COS(PSIair)
  !
  ! E(8) : Y
  !
  E(8) = COS(th)*SIN(PSIair)
  !
  ! E(9) : Z
  !
  E(9) = SIN(th)
  !
  xw = Mw/(Mp+Mw+Ma)                  ! volatiles mass fraction
  !
  !  E(4) : Total energy (thermal + potential + kinetic)
  !
  !  Assumes that water in the atmosphere is all vapour
  hair = (1._wp-wair)*enthalpy_air(Taair) + wair*enthalpy_vapour(Taair)
  E(4) = 2._wp*pi*r*rhoair*ue*(hair+grav*z+0.5_wp*ue*ue)+Cp*T*SUM(E(9+1:9+plume%nc))
  !
  RETURN
END SUBROUTINE dfds
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE get_entrainment_coef(z,R0,U0,rhoair,rho,r,u,Vaair,theta,a_shear,a_vortex)
  !****************************************************************
  !*
  !*   Computes the entrainment coefficients
  !*
  !****************************************************************
  IMPLICIT NONE
  !
  REAL(wp) :: z,R0,U0,rhoair,rho,r,u,Vaair,theta,a_shear,a_vortex
  REAL(wp) :: zs,A,dlnAdz,Ri,c0,c1,c2,c3,c4,h
  !
  zs      = MAX(z/(2.0_wp*R0),0.0_wp)                          ! zstar = z/Dvent
  !
  Ri      = MAX(grav*(rhoair-rho)*r/(rhoair*u*u+1e-6_wp),1e-6_wp      ) ! Ri
  !
  !*** a_vortex
  !
  SELECT CASE(plume%type_av)
    CASE('CONSTANT')
      a_vortex = plume%a_v
    CASE('TATE')
      a_vortex= 0.34_wp*( Vaair*(SQRT(2._wp*ABS(Ri))/U0) )**(-0.125_wp)
      a_vortex = MIN(a_vortex,1.0_wp)
      a_vortex = MAX(a_vortex,0.0_wp)
  END SELECT
  !
  !*** a_shear
  !
  SELECT CASE(plume%type_as)
    CASE('CONSTANT')
      IF (rhoair<rho) THEN  ! jet
        a_shear = plume%a_s_jet
      ELSE                    ! plume
        a_shear = plume%a_s_plume
      END IF
      RETURN  ! we are done
      !
    CASE('KAMINSKI-R')
      IF (rhoair<rho) THEN  ! jet
        c0 = 1.92003_wp
        c1 = 3737.26_wp
        c2 = 4825.98_wp
        c3 = 2d0*(c2-c1)
        c4 = 0.00235_wp
      ELSE                    ! plume
        c0 = 1.61717_wp
        c1 = 478.374_wp
        c2 = 738.348_wp
        c3 = 2._wp*(c2-c1)
        c4 = -0.00145_wp
      END IF
    !
    CASE('KAMINSKI-C')
    !
      IF (rhoair<rho) THEN  ! jet
        c0 = 1.92003_wp
        c1 = 3737.26_wp
        c2 = 4825.98_wp
        c3 = 2._wp*(c2-c1)
        c4 = 0.00235_wp
      ELSE                    ! plume
        c0 = 1.55_wp
        c1 = 329.0_wp
        c2 = 504.5_wp
        c3 = 2._wp*(c2-c1)
        c4 = -0.00145_wp
      END IF
      !
    CASE('OLD')
      IF (rhoair<rho) THEN  ! jet
        c0 = 1.6_wp
        c1 = 1657.0_wp
        c2 = 2411.0_wp
        c3 = 2._wp*(c2-c1)
        c4 = 0.0_wp
      ELSE                    ! plume
        c0 = 1.6_wp
        c1 = 1657.0_wp
        c2 = 2411.0_wp
        c3 = 2._wp*(c2-c1)
        c4 = 0.0_wp
      END IF
       !
  END SELECT
  !
  A      = c0*(zs*zs+c1)/(zs*zs+c2)
  dlnAdz = c3*zs/((zs*zs+c1)*(zs*zs+c2))
  h      = 1._wp/(1._wp+c4*EXP(-5._wp*(zs/10._wp-1._wp)))
  A      = A*h
  !
  a_shear = 0.0675_wp + (1._wp-(1._wp/A))*Ri*SIN(theta) + (0.5_wp*r*dlnAdz)/(2._wp*R0)
  a_shear  = MIN(a_shear,0.17_wp)
  a_shear  = MAX(a_shear,0.0_wp)
  !
  RETURN
END SUBROUTINE get_entrainment_coef
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE costa(f,z,u,r,T,rho,Q,xv,xl,xs)
  !************************************************************
  !*
  !*  Computes aggregation coefficients according to Costa et al. (2010)
  !*
  !************************************************************
  IMPLICIT NONE
  !
  REAL(wp), INTENT(IN) :: f(plume%nc)               ! Class mass flow rate
  REAL(wp), INTENT(IN) :: z,u,r,T,rho,Q,xv,xl,xs
  !
  INTEGER    , SAVE              :: ipass = 0
  REAL   (wp), SAVE, ALLOCATABLE :: Ni(:)          ! Number of particles of size i in an aggregate
  REAL   (wp), SAVE, ALLOCATABLE :: Df(:)          ! Fractal exponent  (size dependent) 
  REAL   (wp), SAVE, ALLOCATABLE :: ka(:)          ! Fractal prefactor (size dependent) 
  REAL   (wp), SAVE, ALLOCATABLE :: psi3j(:)       ! Diameter to volume fractal relationship **3
  REAL   (wp), SAVE, ALLOCATABLE :: psi4j(:)       ! Diameter to volume fractal relationship **4
  !
  REAL   (wp), SAVE              :: sumNi,psi3,psi4
  !
  CHARACTER(LEN=6)  :: wphase
  LOGICAL           :: conve
  INTEGER           :: ic,jc
  REAL   (wp)       :: muair,rhoaT,rhomean
  REAL   (wp)       :: alfa_mean,sumfc,sumff,w,Stij,viscl,alfa,Stcr,qq,conc
  REAL   (wp)       :: Ab,As,Ad,At,gammas,epsilon,dudr,ntot,dntot,vij,fi
  REAL   (wp)       :: expb,exps,expd,factor_t,factor_c,Dmin,dstar
  REAL   (wp)       :: deltaphi,dja,djb,work
  !
  !***  Allocates memory and stores some variables (first time only)
  !
  IF (ipass == 0) THEN
    ipass = 1
    !
    !*** Fractal exponent D_f(d) and fractal prefactor ka(d)
    !
    ALLOCATE(Df(plume%nc))
    ALLOCATE(ka(plume%nc))
    Df(:) = 0.0_wp
    ka(:) = 0.0_wp
    !
    Dmin  = 1.6_wp
    dstar = 2.e-6_wp
    DO ic = part%iaggr,plume%nc-1
      Df(ic) = part%Dfo - (1.36788_wp*(part%Dfo-Dmin))  &
           & / (1.0_wp+EXP((part%diam(ic)-dstar)/dstar))
      !
      ka(ic) = (SQRT(1.56_wp-((1.728_wp-Df(ic)/2.0_wp)**2)) - 0.228_wp)**Df(ic)
      ka(ic) = ka(ic)*( ((2.0_wp+Df(ic))/Df(ic))**(Df(ic)/2.0_wp) )  
    END DO
    !
    !*** Diameter to volume fractal relationship (note that this is actually size dependent)
    !
    ALLOCATE(psi3j(plume%nc))
    ALLOCATE(psi4j(plume%nc))
    psi3j(:) = 0.0_wp
    psi4j(:) = 0.0_wp
    DO ic = part%iaggr,plume%nc-1
      work = part%diam(ic)*((pi*part%diam(ic)*part%diam(ic)   &
         & * part%diam(ic)/6.0_wp)**(-1._wp/Df(ic)))
      psi3j(ic) = work*work*work
      psi4j(ic) = psi3j(ic)*work
    END DO
    !
    !*** However, because kernels are integrated we use by now for simplicity constant values
    psi3 = 6.0_wp/pi
    psi4 =  psi3*((6.0_wp/pi)**(1.0_wp/3.0_wp))   
    !
    !*** Number of particles in an aggregate
    !
    ALLOCATE(Ni(plume%nc))
    Ni(:) = 0.0_wp
    !
    sumNi = 0.0_wp
    DO ic = part%iaggr,plume%nc-1
      Ni(ic) = ka(ic)*( (part%diam(plume%nc)/part%diam(ic))**Df(ic) )
      sumNi  = sumNi + Ni(ic)
    END DO
    !
  END IF
  !
  !*** Initializations (necessary to ensure zero aggregation in vapour)
  !
  A_p(:) = 0.0_wp
  A_m(:) = 0.0_wp
  !
  !*** Select water phase (and do nothing for pure vapor phase)
  !*** Note that we assume that liquid and ice do not coexist
  !
  IF ((xl==0.0_wp).AND.(xs==0.0_wp)) THEN
    wphase = 'vapour'
    RETURN
  ELSE If(xs==0.0_wp) THEN
    wphase = 'liquid'
  ELSE
    wphase = 'ice'
  END IF
  !
  muair  = mua (T)      !  Air viscosity (at bulk temperature)
  rhoaT  = rhoa(z,T)    !  Air density   (at bulk temperature)
  !
  !*** Calculates the mean density of the aggregating particles
  !
  rhomean = 0.0_wp
  sumfc   = 0.0_wp
  DO ic = part%iaggr,plume%nc-1
    rhomean = rhomean + part%fc(ic)*part%rhop(ic)
    sumfc   = sumfc   + part%fc(ic)
  END DO
  rhomean = rhomean/sumfc
  !
  !*** Collision frequency kernel for Brownian motion at each point
  !
  Ab = -4.0_wp*ak*T/(3.0_wp*muair)
  !
  !*** Collision frequency kernel for laminar and turbulent fluid shear
  !
  dudr    = SQRT(2.0_wp*pi)*u/r                 ! laminar
  epsilon = 0.0724_wp*(u*u*u)/r                 ! Smagorinsky-Lilly, factor = 2*sqrt(2)*cs*cs with cs = 0.16
  gammas  = SQRT(epsilon/(muair/rhoaT))    ! turbulent
  gammas  = MAX(gammas,dudr)
  !
  As = -2.0_wp*gammas*psi3/3.0_wp
  !
  !*** Collision frequency kernel for differential settling velocity at each point
  !
  Ad   = -pi*(rhomean-rho)*grav*psi4/(48.0_wp*muair)
  !
  !*** Collision frequency kernel for turbulent inertial kernel at each point
  !
  At = 1.82_wp*(epsilon**(0.75_wp))/(grav*(muair/rhoaT)**(0.25_wp))*Ad
  !
  !*** Compute the settling velocity (using air density and viscosity at bulk temperature)
  !
  DO ic = 1,plume%nc
    CALL vsettl(part%diam(ic),part%rhop(ic),rhoaT,muair,part%vlimit(ic),plume%modv,part%psi(ic),conve)
    IF (.NOT.conve) part%vlimit(ic) = 0.0_wp
  END DO
  IF (plume%aggregation) part%vlimit(plume%nc) = part%vset_aggr*part%vlimit(plume%nc)
  !
  !*** Computes the number of available particles per unit volume and
  !*** the solid volume fraction fi
  !
  ntot  = 0.0_wp
  fi    = 0.0_wp
  DO ic = part%iaggr,plume%nc-1
    IF (ic == part%iaggr) THEN
      dja = part%diam(ic)
      djb = part%diam(ic)*(part%diam(ic)/part%diam(ic+1))         ! Assume same interval
      deltaphi = LOG(djb/dja)/LOG(2._wp)
    ELSE
      dja  = part%diam(ic)
      djb  = part%diam(ic-1)
      deltaphi = LOG(djb/dja)/LOG(2._wp)
    END IF
      conc = rho*f(ic)/Q
      ntot = ntot + 6.0_wp*conc*(1._wp/dja**3._wp-1._wp/djb**3._wp)   & 
         & / (pi*part%rhop(ic)*deltaphi)
      fi = fi + conc/part%rhop(ic)
    END DO
    ntot = ntot/(3._wp*LOG(2._wp))
  !
  !*** Computes the class-averaged sticking efficiency
  !
  SELECT CASE(wphase)
    CASE('ice')
      alfa_mean = 0.09_wp
    CASE('liquid')
      sumff = 0.0_wp
      DO ic = part%iaggr,plume%nc-1
        DO jc = part%iaggr,plume%nc-1
          sumff = sumff + part%fc(ic)*part%fc(jc)
        END DO
      END DO
      !
      viscl = 2.414_wp*(10.0_wp**(247.7_wp/(T-140.0_wp)))  ! water viscosity at bulk temperature
      Stcr  = 1.3_wp                                       ! Critical Stokes number
      qq    = 0.8_wp
      !
      alfa_mean = 0.0_wp
      DO ic = part%iaggr,plume%nc-1
        DO jc = part%iaggr,plume%nc-1
          vij = ABS(part%vlimit(ic)-part%vlimit(jc)) + (8.0_wp*ak*T)   &
            & / (3.0_wp*pi*muair*part%diam(ic)*part%diam(jc))        &
            & +  2.0_wp*gammas*(part%diam(ic)+part%diam(jc))/(3.0_wp*pi)
          w = part%fc(ic)*part%fc(jc)/sumff
          Stij = 8.0_wp*rho*vij*part%diam(ic)*part%diam(jc)            &
             & / (9.0_wp*viscl*(part%diam(ic)+part%diam(jc)))
          alfa = 1.0_wp+((Stij/Stcr)**qq)
          alfa = 1.0_wp/alfa
          alfa_mean = alfa_mean + w*alfa
        END DO
      END DO
  END SELECT
  !
  !*** Computes total particle decay per unit volume and time
  !
  expb = 2.0_wp                      ! Brownian     ntot exponent
  exps = 2.0_wp-(3.0_wp/part%Dfo)    ! Shear        ntot exponent
  expd = 2.0_wp-(4.0_wp/part%Dfo)    ! Differential ntot exponent
  !
  dntot = Ab*(ntot**expb) + As*(ntot**exps)*(fi**(3.0_wp/part%Dfo)) &
      & + (Ad+At)*(ntot**expd)*(fi**(4.0_wp/part%Dfo))
  dntot = ABS(alfa_mean*dntot)  ! positive sign in coefficent A
  !
  !** factor_t for estimating time that single particles spend in a control volume (lagrangian vs eulerian time)
  !
  factor_c = 2.0_wp  ! correction from local formulation to hop-hat  for ^2 terms
  factor_t = 2.0_wp  ! 2-7  also from top-hat to Gaussian
  !
  dntot = factor_c*factor_t*dntot
  !
  !*** Computes aggregation coefficients
  !
  DO ic = part%iaggr,plume%nc-1
    IF (f(ic)<=0.0_wp) THEN
      A_m(ic) = 0.0_wp
    ELSE
      A_m(ic) = (dntot*Ni(ic)/sumNi)*(part%rhop(ic)*pi*part%diam(ic)**3.0_wp/6.0_wp)*pi*r*r
    END IF
  END DO
  A_p(plume%nc) = SUM(A_m(:))
  !
  RETURN
END SUBROUTINE costa
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE vsettl(diam,rhop,rhoa,visc,vset,model,psi,conve)
  !*****************************************************************************
  !*
  !*    Set fall velocity, as a function of particles diameter and
  !*    air density and viscosity.
  !*
  !*    diam     - Particle diameter in meters (Input)
  !*    rhop     - Particle density (Input)
  !*    rhoa     - Air density (Input vector)
  !*    visc     - Air viscosity (Input)
  !*    vset     - Particle settling velocity (Output)
  !*    model    - Settling velocity model (Input)
  !*    psi      - Model factor
  !*
  !*****************************************************************************
  IMPLICIT NONE
  !
  INTEGER :: model, it
  LOGICAL :: conve
  !
  INTEGER :: maxit
  REAL(wp), INTENT(IN) ::  diam,rhop,rhoa,visc,psi
  REAL(wp), INTENT(OUT) :: vset
  REAL(wp) ::  gi,eps,rey,rk1,rk2,a,b,vold,cd100   ! work variables
  REAL(wp), SAVE ::  cd = 1.0_wp
  !
  conve = .true.
  gi=9.81_wp                ! Gravity acceleration
  eps=1e-3_wp               ! Tolerance
  maxit=100                 ! Maximum number of iterations
  vold =1e5_wp              ! A large number
  !
  !***  Begin iteration to compute vset
  SELECT CASE(model)
    CASE(1) ! Model of Arastoopour et al. 1982
      DO it=1,maxit
        vset=SQRT(4.0_wp*gi*diam*rhop/(3.0_wp*cd*rhoa))
        rey=rhoa*vset*diam/visc
        IF (rey <= 988.947_wp) THEN  ! This is the actual transition point
          cd=24.0_wp/rey*(1.0_wp+0.15_wp*rey**0.687_wp)
        ELSE
          cd=0.44_wp
        ENDIF
        IF (it > 1 .AND. ABS(vset-vold) <= eps) THEN
          vset=SQRT(4.0_wp*gi*diam*rhop/(3.0_wp*cd*rhoa))  ! vset with new cd
          EXIT
        ELSE
          vold=vset
        ENDIF
      ENDDO
      IF (vold/=vset) RETURN
    CASE(2) ! Model of Ganser 1993
      DO it=1,maxit
        vset=SQRT(4.0_wp*gi*diam*rhop/(3.0_wp*cd*rhoa))
        rey=rhoa*vset*diam/visc
        rk1=3.0_wp/(1.0_wp+2.0_wp/SQRT(psi))
        rk2=10.0_wp**(1.8148_wp*(-LOG10(psi))**0.5743_wp)
        cd=24.0_wp/(rey*rk1)*(1.0_wp+0.1118_wp           &
        & * (rey*rk1*rk2)**0.6567_wp)+0.4305_wp*rk2      &
        & / (1.0_wp+3305.0_wp/(rey*rk1*rk2))
        IF (it > 1 .AND. ABS(vset-vold) <= eps) THEN
          vset=SQRT(4.0_wp*gi*diam*rhop/(3.0_wp*cd*rhoa))  ! vset with new cd
          EXIT
        ELSE
          vold=vset
        ENDIF
      ENDDO
      IF (vold/=vset) RETURN
    CASE(3) ! Model of Wilson & Huang 1979
      DO it=1,maxit
        vset=SQRT(4.0_wp*gi*diam*rhop/(3.0_wp*cd*rhoa))
        rey=rhoa*vset*diam/visc
        IF (rey <= 100.0_wp) THEN
          cd=24.0_wp/rey*psi**(-0.828_wp)+2.0_wp*SQRT(1.07_wp-psi)
        ELSE IF (rey > 100.0_wp .AND. rey < 1000.0_wp) THEN
          cd100=0.24_wp*psi**(-0.828_wp)+2.0_wp*SQRT(1.07_wp-psi)
          a=(1.0_wp-cd100)/900.0_wp
          b=1.0_wp-1000.0_wp*a
          cd=a*rey+b
        ELSE
          cd=1.0_wp
        ENDIF
        !
        IF (it > 1 .AND. ABS(vset-vold) <= eps) THEN
          vset=SQRT(4.0_wp*gi*diam*rhop/(3.0_wp*cd*rhoa))  ! vset with new cd
          EXIT
        ELSE
          vold=vset
        ENDIF
      ENDDO
      IF (vold/=vset) RETURN
    CASE(4) ! Model of Dellino et al.
      vset=((diam*diam*diam*gi*(rhop-rhoa)*rhoa*(psi**1.6_wp))/(visc*visc))**0.5206_wp
      vset=(1.2065_wp*visc*vset)/(diam*rhoa)
      RETURN
  END SELECT
  !
  conve = .FALSE.
RETURN
  !
END SUBROUTINE vsettl
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE jac
  ! This is a dummy routine needed by lsode
  IMPLICIT NONE
  RETURN
END SUBROUTINE jac
!!
!!-------------------------------------------------------------------------
!!
LOGICAL FUNCTION is_collapse(plume)
  ! Check if column is collapsed (actually checks status flag)
  IMPLICIT NONE
  TYPE(volcanic_plume), INTENT(IN) :: plume
  is_collapse = .FALSE.
  IF (iand(plume%stat,STATUS_COLLAPSE) /= 0) is_collapse = .TRUE.
  RETURN
END FUNCTION is_collapse
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION Ta(z)
  !****************************************************************
  !*
  !*   Gets air temperature (K) at z (in m a.s.l.)
  !*
  !****************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAL     :: go_on
  INTEGER     :: iz1,iz2
  REAL(wp)    :: f1,f2

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:Ta(z):"  
  !
  !*** Values interpolated from a profile
  !
  IF (z<profile%z(1)) THEN
    Ta = profile%T(1)
    RETURN
  END IF
  !
  IF (z>profile%z(profile%nz)) THEN
    Ta = profile%T(profile%nz)
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2  and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1) .AND. z<profile%z(iz2)) THEN
      go_on = .FALSE.
      f1 = 1.0_wp-((z-profile%z(iz1))/(profile%z(iz2)-profile%z(iz1)))
      f2 = 1.0_wp-f1
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function Ta')
    END IF
  END DO
  !
  Ta = f1*profile%T(iz1)+f2*profile%T(iz2)
  !
  RETURN
END FUNCTION Ta
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION Pa(z)
  !**************************************************************
  !*
  !*   Gets air pressure (in Pa) at z (in m a.s.l.)
  !*
  !**************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAL     :: go_on
  INTEGER     :: iz1,iz2
  REAL(wp)    :: f1,f2

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:Pa(z):" 
  !
  !*** Values interpolated from a profile
  !
  IF (z<profile%z(1)) THEN
    Pa = profile%P(1)
    RETURN
  END IF
  !
  IF (z>profile%z(profile%nz)) THEN
    Pa = profile%P(profile%nz)
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2  and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1).AND.z<profile%z(iz2)) THEN
      go_on = .FALSE.
      f1 = 1.0_wp-((z-profile%z(iz1))/(profile%z(iz2)-profile%z(iz1)))
      f2 = 1.0_wp-f1
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function Pa')
    END IF
  END DO
  !
  Pa = f1*profile%P(iz1)+f2*profile%P(iz2)
  !
  RETURN
END FUNCTION Pa
!!
!!-------------------------------------------------------------------------
!!
REAl(wp) FUNCTION rhoa_z(z)
  !**************************************************************
  !*
  !*    Interpolates air density at z (in m a.s.l.)
  !*
  !**************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAL     :: go_on
  INTEGER     :: iz1,iz2
  REAL(wp)    :: f1,f2

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:rhoa_z(z):" 
  !
  !*** Values interpolated
  !
  IF (z<profile%z(1)) THEN
    rhoa_z = profile%rho(1)
    RETURN
  END IF
  !
  IF (z>profile%z(profile%nz)) THEN
    rhoa_z = profile%rho(profile%nz)
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2 and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1).AND.z<profile%z(iz2)) THEN
      go_on = .FALSE.
      f1 = 1.0_wp-((z-profile%z(iz1))/(profile%z(iz2)-profile%z(iz1)))
      f2 = 1.0_wp-f1
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function rhoa')
    END IF
  END DO
  !
  rhoa_z = f1*profile%rho(iz1)+f2*profile%rho(iz2)
  !
  RETURN
END FUNCTION rhoa_z
  !
  ! Air density (function of z and T)
  !
REAL(wp) FUNCTION rhoa_z_T(z,T)
  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: z,T
  rhoa_z_T = rhoa(z)*Ta(z)/T
  !
  RETURN
END FUNCTION rhoa_z_T
  !
  !
  !
REAL(wp) FUNCTION sh_a(z)
  !**************************************************************
  !*
  !*    Interpolates air specific humidity at z (in m a.s.l.)
  !*
  !**************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAL     :: go_on
  INTEGER     :: iz1,iz2
  REAL(wp)    :: f1,f2

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:sh_a(z):" 
  !
  !*** Values interpolated
  !
  IF (z<profile%z(1)) THEN
    sh_a = profile%sh(1)
    RETURN
  END IF
  !
  IF (z>profile%z(profile%nz)) THEN
    sh_a = profile%sh(profile%nz)
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2 and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1).AND.z<profile%z(iz2)) THEN
      go_on = .FALSE.
      f1 = 1.0_wp-((z-profile%z(iz1))/(profile%z(iz2)-profile%z(iz1)))
      f2 = 1.0_wp-f1
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function sh_a')
    END IF
  END DO
  !
  sh_a = f1*profile%sh(iz1)+f2*profile%sh(iz2)
  !
  RETURN
END FUNCTION sh_a
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION Va(z)
  !**************************************************************
  !*
  !*    Interpolates air velocity (m/s) at z (in m a.s.l.)
  !*    from Vair
  !*
  !**************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAl     :: go_on
  INTEGER     :: iz1,iz2
  REAL(wp)    :: f1,f2

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:Va(z):" 
  !
  !*** Values interpolated from a profile
  !
  IF (z<profile%z(1)) THEN
    Va = profile%u(1)
    RETURN
  END IF
  !
  IF (z>profile%z(profile%nz)) THEN
    Va = profile%u(profile%nz)
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2  and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1).AND.z<profile%z(iz2)) THEN
      go_on = .FALSE.
      f1 = 1.0_wp-((z-profile%z(iz1))/(profile%z(iz2)-profile%z(iz1)))
      f2 = 1.0_wp-f1
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function Va')
    END IF
  END DO
  !
  Va = f1*profile%u(iz1)+f2*profile%u(iz2)
  !
  RETURN
END FUNCTION Va
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION Da(z)
  !**************************************************************
  !*
  !*    Interpolates air direction
  !*
  !**************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAL     :: go_on
  INTEGER     :: iz1,iz2
  REAL(wp)    :: f1,f2

  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:Da:" 
  !
  !*** Values interpolated from a profile
  !
  IF (z<profile%z(1)) THEN
    Da = profile%PSIa(1)
    RETURN
  END IF
  !
  IF (z>profile%z(profile%nz)) THEN
    Da = profile%PSIa(profile%nz)
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2  and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1).AND.z<profile%z(iz2)) THEN
      go_on = .FALSE.
      f1 = 1.0_wp-((z-profile%z(iz1))/(profile%z(iz2)-profile%z(iz1)))
      f2 = 1.0_wp-f1
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function Da')
    END IF
  END DO
  !
  Da = f1*profile%PSIa(iz1)+f2*profile%PSIa(iz2)
  !
  RETURN
END FUNCTION Da
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION dVadz(z)
  !**************************************************************
  !*
  !*    Interpolates air velocity gradient
  !*    from Vair
  !*
  !**************************************************************
  IMPLICIT NONE
  REAL(wp) :: z
  !
  LOGICAL     :: go_on
  INTEGER     :: iz1,iz2
  
  CHARACTER(LEN=MAX_CHAR_LENGTH)    :: thisroutine="mo_art_fplume_bpt:dVadz:" 
  !
  !*** Values interpolated
  !
  IF (z<profile%z(1)) THEN
    dVadz = profile%u(1)/profile%z(1)
    RETURN
  END IF
  !
  IF (z>=profile%z(profile%nz)) THEN
    dVadz = (profile%u(profile%nz)- profile%u(profile%nz-1)) &
        & / (profile%z(profile%nz)- profile%z(profile%nz-1))
    RETURN
  END IF
  !
  !*** Gets the position indexes iz1,iz2  and the weights
  !
  iz1 = 0
  go_on = .TRUE.
  DO WHILE (go_on)
    iz1 = iz1 + 1
    iz2 = iz1 + 1
    IF (z>=profile%z(iz1).AND.z<profile%z(iz2)) THEN
      go_on = .FALSE.
      dVadz = (profile%u(iz2)- profile%u(iz1)) &
          & / (profile%z(iz2)- profile%z(iz1))
    ELSE IF (iz2==profile%nz) THEN
      CALL finish(thisroutine,'Source position not found in function dVadz')
    END IF
  END DO
  !
  RETURN
END FUNCTION dVadz
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION mua(T)
  !***********************************************************
  !*
  !*    Computes air viscosity at temperature T
  !*
  !************************************************************
  IMPLICIT NONE
  REAL(wp) :: T,mua0,Tair0
  !
  mua0  = 1.827e-5_wp  ! reference viscosity
  Tair0 = 291.15_wp    ! reference temperature
  !
  !*** Sutherland's law
  !
  mua = mua0*((Tair0+120.0_wp)/(T+120.0_wp))*((T/Tair0)**1.5_wp)
  !
  RETURN
END FUNCTION mua
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE set_latent_heat_flag(lh)
  ! This subroutine sets the logical value of the variable latent_heat
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: lh
  plume%latent_heat = lh
  IF (plume%latent_heat .EQV. .TRUE.) THEN
    hl0=hs0+alf
    hv0=hl0+alv
    Cv = Cv0
    Cl = Cl0
    Cs = Cs0
  ELSE
  ! No latent heat
    hl0= Cw*t3
    hv0= Cw*t3
    Cv = Cw
    Cl = Cw
    Cs = Cw
  END IF
  !
  RETURN
END SUBROUTINE set_latent_heat_flag
!!
!!-------------------------------------------------------------------------
!!
REAL(wp) FUNCTION enthalpy_particles(temp)
  ! Enthalpy per unit mass of solid particles
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: temp  ! Temperature
  enthalpy_particles = Cp*temp
  RETURN
END FUNCTION enthalpy_particles

REAL(wp) FUNCTION enthalpy_air(temp)
  ! Enthalpy per unit mass of air
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: temp  ! Temperature
  enthalpy_air = Ca*temp
  RETURN
END FUNCTION enthalpy_air

REAL(wp) FUNCTION enthalpy_vapour(temp)
  ! Enthalpy per unit mass of water vapor
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: temp  ! Temperature
  enthalpy_vapour = hv0 + Cv*(temp-t3)
  RETURN
END FUNCTION enthalpy_vapour

REAL(wp) FUNCTION enthalpy_liquid(temp)
  ! Enthalpy per unit mass of liquid water
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: temp  ! Temperature
  enthalpy_liquid = hl0 + Cl*(temp-t3)
  RETURN
END FUNCTION enthalpy_liquid

REAL(wp) FUNCTION enthalpy_ice(temp)
  ! Enthalpy per unit mass of ice
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: temp  ! Temperature
  enthalpy_ice = hs0 + Cs*(temp-t3)
  RETURN
END FUNCTION enthalpy_ice

REAL(wp) FUNCTION density_air(temp,pres)
  ! Density of air
  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: temp,pres
  density_air = pres/(rair*temp)
  RETURN
END FUNCTION density_air

REAL(wp) FUNCTION density_vapour(temp,pres)
  ! Density of vapour
  IMPLICIT NONE
  REAL(wp), INTENT(IN) :: temp,pres
  density_vapour = pres/(rvapour*temp)
  RETURN
END FUNCTION density_vapour

SUBROUTINE get_gas_molar_fractions(xa,xv,na,nv)
  ! Get molar fractions of vapour and air in the gas phase (na+nv=1)
  ! Safe result for xa=xv=0
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: xa,xv
  REAL(wp), INTENT(OUT) :: na,nv
  REAL(wp) :: ntot
  ntot = xa/(amd*1.e-3_wp)+xv/(amw*1.e-3_wp)
  IF (ntot /= 0._wp) THEN
    na = xa/(amd*1.e-3_wp*ntot)
    nv = xv/(amw*1.e-3_wp*ntot)
  ELSE
    na = 0._wp
    nv = 0._wp
  END IF
  RETURN
END SUBROUTINE get_gas_molar_fractions

SUBROUTINE saturation_pressure_over_liquid(T,psat)
  ! Vapor saturation pressure over liquid in Pa (T > Tref)
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: T
  REAL(wp), INTENT(OUT) :: psat
  psat = 611.22_wp*EXP(17.67_wp*(T-273.16_wp)/(T-29.65_wp))
  RETURN
END SUBROUTINE saturation_pressure_over_liquid

SUBROUTINE saturation_pressure_over_ice(T,psat)
  ! Vapor saturation pressure over solid in Pa (T < Tref)
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: T
  REAL(wp), INTENT(OUT) :: psat
  psat = -9.09718_wp*(273.16_wp/T-1._wp)-3.56654_wp*Log10(273.16_wp/T)   &
     & + 0.87679_wp*(1._wp-T/273.16_wp)
  ! psat = 610.71d0*(10.0d0**psat)  ! Original
  psat = 611.22_wp*(10.0_wp**psat)    ! Modified by G.Macedonio
  RETURN
END SUBROUTINE saturation_pressure_over_ice

SUBROUTINE saturation_pressure(temp,psat)
  ! Vapor saturation pressure (over liquid/solid)
  IMPLICIT NONE
  REAL(wp), INTENT(IN)  :: temp
  REAL(wp), INTENT(OUT) :: psat
  IF (temp <= t3) THEN
    CALL saturation_pressure_over_ice(temp,psat)
  ELSE
    CALL saturation_pressure_over_liquid(temp,psat)
  END IF
  RETURN
END SUBROUTINE saturation_pressure

SUBROUTINE temperature_at_saturation_over_liquid(T,P)
  ! This is T(Psat): the inverse of Psat(T) for T>Tref
  ! Defined for P>0
  IMPLICIT NONE
  REAL(wp), INTENT(OUT) :: T
  REAL(wp), INTENT(IN)  :: P
  REAL(wp) :: a
  a = 17.67_wp/LOG(P/611.22_wp)
  T = (29.65_wp-a*273.16_wp)/(1._wp-a)
  RETURN
END SUBROUTINE temperature_at_saturation_over_liquid

SUBROUTINE temperature_at_saturation_over_ice(T,P)
  ! This is T(Psat): the inverse of Psat(T) for T<Tref and P < Pref=611.2 Pa
  IMPLICIT NONE
  REAL(wp), INTENT(OUT) :: T
  REAL(wp), INTENT(IN)  :: P
  REAL(wp) :: pleft,pright,pLeftMin,ttry,ptry
  REAL(wp) :: tleft,tright,tleftMin
  REAL(wp) :: tol = 1.e-6_wp   ! Tolerance
  INTEGER  :: iiter
  !
  IF (P > Pref) THEN
    T = t3 
    RETURN
  END IF
  tLeftMin = 100.0_wp   ! Set minimum search Temperature
  CALL saturation_pressure_over_ice(tLeftMIn,pLeftMin)
  IF (P < pLeftMin) THEN
    T = TLeftMin
    RETURN
  END IF
  tleft  = tLeftMin         ! Search between tleft and tright
  tright = 274.0_wp
  CALL saturation_pressure_over_ice(tleft,pleft)
  IF (pright == P) THEN
    T = tleft
    RETURN
  END IF
  ! Try right value
  CALL saturation_pressure_over_ice(tright,pright)
  IF (pright == P) THEN
    T = tright
    RETURN
  END IF
  IF (P < pleft .or. P > pright) THEN
    ! Out of range (should not occur)
    WRITE(*,'(''Error in temperature_at_saturation_over_ice'')')
    WRITE(*,'(''P:'',3(1x,e12.5))') P,pleft,pright
    WRITE(*,'(''T:'',3(1x,e12.5))') T,tleft,tright
    CALL finish('mo_art_fplume_bpt:temperature_at_saturation_over_ice',      &
             &  'Error in temperature_at_saturation_over_ice')
  END IF
  ! Bisection
  iiter = 0
  DO 
    iiter = iiter + 1
    ttry = 0.5_wp*(tleft+tright)
    CALL saturation_pressure_over_ice(ttry,ptry)
    IF (ptry == P) THEN
      T = ttry
      RETURN
    END IF
    IF (ptry > P) THEN
      tright = ttry
    ELSE
      tleft = ttry
    END IF
    IF (ABS(tleft-tright) < tol) EXIT
      IF (iiter==100) EXIT
  END DO
  !
  T = 0.5_wp*(tleft+tright)
  !
  RETURN
END SUBROUTINE temperature_at_saturation_over_ice

SUBROUTINE temperature_mixture(xa,xp,xw,temp,pres,xv,xl,xs,pa,pv,enthalpy)
  ! Returns the temperature of air+particles+water mixture with known
  ! total enthalpy. It also evaluates xv,xl,xs
  ! The sum xw+xl+xp must be 1 (no check)
  IMPLICIT NONE
  REAL(wp), INTENT(IN)    :: xa,xp,xw  ! Air, particles and water fractions
  REAL(wp), INTENT(INOUT) :: temp      ! Temperature
  REAL(wp), INTENT(IN)    :: pres      ! Total pressure
  REAL(wp), INTENT(INOUT) :: xv,xl,xs  ! Mass fractions of vapour,liquid,ice
  REAL(wp), INTENT(INOUT) :: pa,pv     ! Partial pressures of air and vapour
  REAL(wp), INTENT(IN)    :: enthalpy  ! Enthalpy of the mixture
  REAL(wp) :: na,nv                  ! Molar fractions of air and vapour
  REAL(wp) :: tsi,tsl
  REAL(wp) :: h1,h2,h3,h4,hpi,hpl,hai,hal,hli,hlv,hs,hv
  INTEGER  :: i,inzone,outzone
  !
  IF (xa+xw == 0._wp) THEN
    ! Only particles
    xl = 0._wp
    xs = 0._wp
    xv = 0._wp
    pv = 0._wp
    pa = 0._wp
    temp = enthalpy/Cp
    RETURN
  END IF
  IF (xw == 0._wp) THEN
    ! Air + particles (no water)
    xl = 0._wp
    xs = 0._wp
    xv = 0._wp
    pv = 0._wp
    pa = pres
    temp = enthalpy/(xp*Cp+xa*Ca)   ! Assume constant specific heats
    RETURN
  END IF
  !
  DO i=1,2  ! Almost two tentatives
    CALL get_gas_molar_fractions(xa,xw,na,nv)
    pv = nv*pres ! Guess pressure (assume all water is vapour)
    IF (i==1) inzone = get_zone(pv,Pref)   ! Tentative zone at first step
    IF (inzone == 1) THEN       
      ! Only ice + vapor
      CALL temperature_at_saturation_over_ice(tsi,pv)
      hs = enthalpy_ice(tsi)      ! Enthalpy of ice
      hv = enthalpy_vapour(tsi)   ! Enthalpy of vapour
      hpi = enthalpy_particles(tsi)
      hai = enthalpy_air(tsi)
      h1 = xp*hpi + xa*hai + xw*hs
      h2 = xp*hpi + xa*hai + xw*hv
      IF (enthalpy <= h1) THEN
        ! Only ice + particles + air
        xv = 0._wp
        xl = 0._wp
        xs = xw
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres   ! Ensures that pa=0 when xa=0
        temp = enthalpy/(xp*Cp+xa*Ca+xs*Cs)
        RETURN
      ELSE IF (enthalpy <= h2) THEN
        ! Vapour + ice + particles + air at T=tsat_ice (h1 < enthalpy <= h2)
        xv = (enthalpy-xp*hpi-xa*hai-xw*hs)/(hv-hs)
        xl = 0._wp
        xs = xw - xv
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = tsi
      ELSE
        ! Only vapour + particles + air
        xv = xw
        xl = 0._wp
        xs = 0._wp
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = (enthalpy-xv*(hv0-Cv*t3))/(xp*Cp+xa*Ca+xv*Cv)
      END IF
    ELSE     ! pv > Pref
      CALL temperature_at_saturation_over_ice(tsi,pv)
      CALL temperature_at_saturation_over_liquid(tsl,pv)
      hs  = enthalpy_ice(tsi)     ! Enthalpy of ice at T=tsi
      hli = enthalpy_liquid(tsi)  ! Enthalpy of liquid at T=tsi
      hlv = enthalpy_liquid(tsl)  ! Enthalpy of liquid at T=tsl
      hv  = enthalpy_vapour(tsl)  ! Enthalpy of vapour at T=tsl
      hpi = enthalpy_particles(tsi)
      hpl = enthalpy_particles(tsl)
      hai = enthalpy_air(tsi)
      hal = enthalpy_air(tsl)
      h1 = xp*hpi + xa*hai + xw*hs
      h2 = xp*hpi + xa*hai + xw*hli
      h3 = xp*hpl + xa*hal + xw*hlv
      h4 = xp*hpl + xa*hal + xw*hv
      IF (enthalpy <= h1) THEN
        ! Only ice + particles + air
        xv = 0._wp
        xl = 0._wp
        xs = xw
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = enthalpy/(xp*Cp+xa*Ca+xs*Cs)
        RETURN
      ELSE IF (enthalpy <= h2) THEN  ! (h1 < enthalpy <= h2)
        ! Ice + liquid + particles + air at T=tsi
        xv = 0._wp
        xs = (xw*hli+xp*Cp*tsi+xa*Ca*tsi-enthalpy)/(hli-hs)
        xl = xw-xs
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = tsi
        RETURN
      ELSE IF (enthalpy <= h3) THEN  ! (h2 < enthalpy <= h3)
        ! Liquid + particles + air
        xv = 0._wp
        xl = xw
        xs = 0._wp
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = (enthalpy-xw*(hl0-Cl*t3))/(xp*Cp+xa*Ca+xw*Cl)
        RETURN
      ELSE IF (enthalpy <= h4) THEN  ! (h3 < enthalpy <= h4)
        ! Liquid + vapour + particles + air at T=tsl
        xv = (enthalpy-xp*Cp*tsl-xa*Ca*tsl-xw*hlv)/(hv-hlv)
        xl = xw-xv
        xs = 0._wp
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = tsl
      ELSE ! (enthalpy > h4)
        ! Vapour + particles + air
        xv = xw
        xl = 0._wp
        xs = 0._wp
        CALL get_gas_molar_fractions(xa,xv,na,nv)
        pv = nv*pres
        pa = na*pres
        temp = (enthalpy-xw*(hv0-Cv*t3))/(xp*Cp+xa*Ca+xw*Cv)
      END IF
    END IF
    outzone = get_zone(pv,pref)
    IF (outzone == inzone) EXIT
    inzone = outzone   ! Try with the new zone
  END DO
  RETURN
END SUBROUTINE temperature_mixture

FUNCTION get_zone(pv,pref)
  ! Internal function
  IMPLICIT NONE
  REAL(wp)   :: pv,Pref
  INTEGER    :: get_zone
  IF (pv <= Pref) THEN
    get_zone = 1
  ELSE
    get_zone = 2
  END IF
  RETURN
END FUNCTION get_zone
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_fplume_bpt
