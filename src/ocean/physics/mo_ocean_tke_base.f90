! This module contains the main computations of diffusivities based on
! TKE (following Gaspar'90)  with the calculation of the mixing length following (Blanke, B., P. Delecluse)
!
! @see  Gaspar, P., Y. Gregoris, and J.-M. Lefevre
!       J. Geophys. Res., 95(C9), 16179-16193, doi:10.1029/JC095iC09p16179.
!
! @see  Blanke, B., P. Delecluse
!       J. Phys. Oceanogr., 23, 1363-1388. doi: http://dx.doi.org/10.1175/1520-0485(1993)023<1363:VOTTAO>2.0.CO;2
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
MODULE mo_ocean_tke_base
  USE mo_kind,               ONLY: wp
  USE mo_ocean_math_utils,   ONLY: solve_tridiag, solve_tridiag_block
  USE mo_exception,          ONLY: finish
  USE mo_fortran_tools,      ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE
  SAVE

  !public member functions
  PUBLIC :: init_tke
  PUBLIC :: coeffs_tke
  PUBLIC :: put_tke


!=================================================================================
!---------------------------------------------------------------------------------
! Interface to call the TKE parameterization!
!---------------------------------------------------------------------------------

interface coeffs_tke
  module procedure integrate_tke              ! calculation if prognostic TKE equation
  module procedure integrate_tke_gpu
  module procedure integrate_tke_block
end interface coeffs_tke

!---------------------------------------------------------------------------------
! Interface to put values to TKE variables
!---------------------------------------------------------------------------------

interface put_tke
  module procedure tke_put_tke_int
  module procedure tke_put_tke_real
  module procedure tke_put_tke_logical
end interface put_tke

!=================================================================================

! types for TKE
type, public :: tke_type
private

real(wp)             ::  &
  c_k                   ,& !
  c_eps                 ,& ! dissipation parameter
  cd                    ,& !
  alpha_tke             ,& !
  clc                   ,& ! factor for Langmuir turbulence
  mxl_min               ,& ! minimum value for mixing length
  KappaM_min            ,& ! minimum value for Kappa momentum
  KappaH_min            ,& ! minimum value for Kappa tracer
  KappaM_max            ,& ! maximum value for Kappa momentum
  tke_surf_min          ,& ! minimum value for surface TKE
  tke_min                  ! minimum value for TKE, necessary to set this value when
                           ! run without IDEMIX, since there are no sources for TKE in the deep ocean otherwise

integer               :: &
  tke_mxl_choice        ,& ! choice of calculation of mixing length;
  handle_old_vals          ! Flag for what to do with old values of Vmix_vars%TKE,tke_diss,KappaM_iface,KappaH_iface

logical                :: &
  only_tke               ,&
  l_lc                   ,&
  use_Kappa_min          ,&
  use_ubound_dirichlet   ,&
  use_lbound_dirichlet

end type tke_type

type(tke_type), target :: tke_constants_saved

CHARACTER(LEN=*), PARAMETER :: module_name = 'tke'


 contains

!=================================================================================

subroutine init_tke(c_k, c_eps, cd, alpha_tke, mxl_min, &
                    use_Kappa_min, KappaM_min, KappaH_min, KappaM_max, &
                    tke_mxl_choice, use_ubound_dirichlet, use_lbound_dirichlet, &
                    handle_old_vals, only_tke, l_lc, clc, tke_min, tke_surf_min, &
                    tke_userdef_constants)

! This subroutine sets user or default values for TKE parameters

real(wp),optional, intent(in)                 ::  &
  c_k                                            ,&
  c_eps                                          ,&
  cd                                             ,&
  alpha_tke                                      ,&
  mxl_min                                        ,&
  KappaM_min                                     ,&
  KappaH_min                                     ,&
  KappaM_max                                     ,&
  tke_surf_min                                   ,&
  clc                                            ,&
  tke_min

integer, intent(in),optional                   :: &
  tke_mxl_choice                                 ,&
  handle_old_vals

logical, intent(in), optional                  :: &
  only_tke                                       ,&
  l_lc                                           ,&
  use_Kappa_min                                  ,&
  use_ubound_dirichlet                           ,&
  use_lbound_dirichlet

type(tke_type), intent(inout),target, optional :: &
  tke_userdef_constants

CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':init_tke'


! FIXME: not sure about the allowed ranges for TKE parameters
if (present(c_k)) then
  if(c_k.lt. 0.0d0 .or. c_k .gt. 1.5d0) then
!    print*, "ERROR:c_k can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:c_k can only be allowed_range')
  end if
  call put_tke('c_k', c_k, tke_userdef_constants)
else
  call put_tke('c_k',0.1d0 , tke_userdef_constants)
end if

if (present(c_eps)) then
  if(c_eps.lt. 0.d0 .or. c_eps .gt. 10.d0) then
!    print*, "ERROR:c_eps can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:c_eps can only be allowed_range')
  end if
  call put_tke('c_eps', c_eps, tke_userdef_constants)
else
  call put_tke('c_eps', 0.7d0, tke_userdef_constants)
end if

if (present(cd)) then
  if(cd.lt. 0.1d0 .or. cd .gt. 30.d0) then
!    print*, "ERROR:cd can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:cd can only be allowed_range')
  end if
  call put_tke('cd', cd, tke_userdef_constants)
else
  call put_tke('cd', 3.75d0, tke_userdef_constants)
end if

if (present(alpha_tke)) then
  if(alpha_tke.lt. 1.d0 .or. alpha_tke .gt. 90.d0) then
!    print*, "ERROR:alpha_tke can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:alpha_tke can only be allowed_range')
 end if
  call put_tke('alpha_tke', alpha_tke, tke_userdef_constants)
else
  call put_tke('alpha_tke', 30.d0, tke_userdef_constants)
end if

if (present(mxl_min)) then
  if(mxl_min.lt. 1.d-12 .or. mxl_min .gt. 0.4d0) then
!    print*, "ERROR:mxl_min can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:mxl_min can only be allowed_range')
  end if
  call put_tke('mxl_min', mxl_min, tke_userdef_constants)
else
  call put_tke('mxl_min', 1.d-8, tke_userdef_constants)
end if

if (present(use_Kappa_min)) then
  call put_tke('use_Kappa_min', use_Kappa_min, tke_userdef_constants)
else
  call put_tke('use_Kappa_min', .false., tke_userdef_constants)
end if

if (present(KappaM_min)) then
!  if(KappaM_min.lt. 0.d0 .or. KappaM_min .gt. 1.d0) then
!!    print*, "ERROR:KappaM_min can only be allowed_range"
!!    stop 1
!    CALL finish(method_name,'ERROR:KappaM_min can only be allowed_range')
!  end if
  call put_tke('KappaM_min', KappaM_min, tke_userdef_constants)
else
  call put_tke('KappaM_min', 1.d-4, tke_userdef_constants)
end if

if (present(KappaH_min)) then
!  if(KappaH_min.lt. 0.d0 .or. KappaH_min .gt. 1.d0) then
!!    print*, "ERROR:KappaH_min can only be allowed_range"
!!    stop 1
!    CALL finish(method_name,'ERROR:KappaH_min can only be allowed_range')
!  end if
  call put_tke('KappaH_min', KappaH_min, tke_userdef_constants)
else
  call put_tke('KappaH_min', 1.d-5, tke_userdef_constants)
end if

if (present(KappaM_max)) then
  if(KappaM_max.lt. 10.d0 .or. KappaM_max .gt. 1000.d0) then
!    print*, "ERROR:KappaM_max can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:KappaM_max can only be allowed_range')
  end if
  call put_tke('KappaM_max', KappaM_max, tke_userdef_constants)
else
  call put_tke('KappaM_max', 100.d0, tke_userdef_constants)
end if

if (present(tke_mxl_choice)) then
  if(tke_mxl_choice.lt. 1 .or. tke_mxl_choice .gt. 2 ) then
!    print*, "ERROR:tke_mxl_choice can only be 1 or 2"
!    stop 1
    CALL finish(method_name,'ERROR:tke_mxl_choice can only be 1 or 2')
  end if
  call put_tke('tke_mxl_choice', tke_mxl_choice, tke_userdef_constants)
else
  call put_tke('tke_mxl_choice', 2, tke_userdef_constants)
end if

if (present(handle_old_vals)) then
  if(handle_old_vals.lt. 1 .or. handle_old_vals.gt. 3 ) then
!    print*, "ERROR:handle_old_vals can only be 1 to 3"
!    stop 1
    CALL finish(method_name,'ERROR:handle_old_vals can only be 1 to 3')
  end if
  call put_tke('handle_old_vals', handle_old_vals, tke_userdef_constants)
else
  call put_tke('handle_old_vals', 1, tke_userdef_constants)
end if

if (present(clc)) then
  if(clc.lt. 0.0 .or. clc .gt. 30.0) then
    print*, "ERROR:clc can only be allowed_range"
    stop 1
  end if
  call put_tke('clc', clc, tke_userdef_constants)
else
  call put_tke('clc',0.15d0 , tke_userdef_constants)
end if


if (present(tke_min)) then
  if(tke_min.lt. 1.d-9 .or. tke_min.gt. 1.d-2 ) then
!    print*, "ERROR:tke_min can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:tke_min can only be allowed_range')
  end if
  call put_tke('tke_min', tke_min, tke_userdef_constants)
else
  call put_tke('tke_min', 1.d-6, tke_userdef_constants)
end if

if (present(tke_surf_min)) then
  if(tke_surf_min.lt. 1.d-7 .or. tke_surf_min.gt. 1.d-2 ) then
!    print*, "ERROR:tke_surf_min can only be allowed_range"
!    stop 1
    CALL finish(method_name,'ERROR:tke_surf_min can only be allowed_range')
  end if
  call put_tke('tke_surf_min', tke_surf_min, tke_userdef_constants)
else
  call put_tke('tke_surf_min', 1.d-4, tke_userdef_constants)
end if

if (present(use_ubound_dirichlet)) then
  call put_tke('use_ubound_dirichlet', use_ubound_dirichlet, tke_userdef_constants)
else
  call put_tke('use_ubound_dirichlet', .false., tke_userdef_constants)
end if

if (present(use_lbound_dirichlet)) then
  call put_tke('use_lbound_dirichlet', use_lbound_dirichlet, tke_userdef_constants)
else
  call put_tke('use_lbound_dirichlet', .false., tke_userdef_constants)
end if

if (present(only_tke)) then

  call put_tke('only_tke', only_tke, tke_userdef_constants)
else
  call put_tke('only_tke', .true., tke_userdef_constants)
end if

if (present(l_lc)) then

  call put_tke('l_lc', l_lc, tke_userdef_constants)
else
  call put_tke('l_lc', .false., tke_userdef_constants)
end if

end subroutine init_tke

!=================================================================================

subroutine integrate_tke_gpu(                      &
                             start_index,          &
                             end_index,            &
                             nproma,               &
                             levels,               &
                             tke_old,              &
                             tke_new,              &
                             KappaM_out,           &
                             KappaH_out,           &
                             vmix_int_1,           &
                             vmix_int_2,           &
                             vmix_int_3,           &
                             dzw,                  &
                             dzt,                  &
                             max_nlev,             &
                             Ssqr,                 &
                             Nsqr,                 &
                             tke_Tbpr,             & ! diagnostic
                             tke_Tspr,             & ! diagnostic
                             tke_Tdif,             & ! diagnostic
                             tke_Tdis,             & ! diagnostic
                             tke_Twin,             & ! diagnostic
                             tke_Tiwf,             & ! diagnostic
                             tke_Tbck,             & ! diagnostic
                             tke_Ttot,             & ! diagnostic
                             tke_Lmix,             & ! diagnostic
                             tke_Pr,               & ! diagnostic
                             tke_plc,              & ! langmuir turbulence
                             forc_tke_surf,        &
                             E_iw,                 &
                             dtime,                &
                             iw_diss,              &
                             forc_rho_surf,        &
                             rho_ref,              &
                             grav,                 &
                             alpha_c,              &
                             lacc)

  integer, intent(in)                  :: start_index
  integer, intent(in)                  :: end_index
  integer, intent(in)                  :: nproma
  integer, intent(in)                  :: max_nlev
  integer, intent(in)                  :: levels(nproma)
  real(wp), intent(in)                 :: tke_old(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_new(nproma,max_nlev+1)
  real(wp), intent(out)                :: KappaM_out(nproma,max_nlev+1)
  real(wp), intent(out)                :: KappaH_out(nproma,max_nlev+1)
  real(wp), intent(out)                :: vmix_int_1(nproma,max_nlev+1)
  real(wp), intent(out)                :: vmix_int_2(nproma,max_nlev+1)
  real(wp), intent(out)                :: vmix_int_3(nproma,max_nlev+1)
  real(wp), intent(in)                 :: dzw(nproma,max_nlev)
  real(wp), intent(in)                 :: dzt(nproma,max_nlev+1)
  real(wp), intent(in)                 :: Ssqr(nproma,max_nlev+1)
  real(wp), intent(in)                 :: Nsqr(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Tbpr(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Tspr(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Tdif(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Tdis(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Twin(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Tiwf(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Tbck(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Ttot(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Lmix(nproma,max_nlev+1)
  real(wp), intent(out)                :: tke_Pr(nproma,max_nlev+1)
  real(wp), intent(in), optional       :: tke_plc(nproma,max_nlev+1)
  real(wp), intent(in)                 :: forc_tke_surf(nproma)
  real(wp), intent(in), optional       :: E_iw(nproma,max_nlev+1)
  real(wp), intent(in)                 :: dtime
  real(wp), intent(in), optional       :: iw_diss(nproma,max_nlev+1)
  real(wp), intent(in)                 :: forc_rho_surf(nproma)
  real(wp), intent(in)                 :: rho_ref
  real(wp), intent(in)                 :: grav
  real(wp), intent(in), optional       :: alpha_c(nproma,max_nlev+1)
  logical, intent(in), optional        :: lacc

  ! local variables
  type(tke_type), pointer              :: tke_constants_in
  logical :: lzacc

  real(wp), dimension(nproma,max_nlev+1)       :: &
    tke_unrest                                  , & ! copy of tke before restorring to background value
    tke_upd                                     , & ! copy of tke before in which surface/bottom values given by Dirichlet boundary conditions
    mxl                                         , & ! mixing length scale (m)
    sqrttke                                     , & ! square root of TKE (m/s)
    prandtl                                     , & ! Prandtl number
    Rinum                                       , & ! Richardson number
    K_diss_v                                    , & ! shear production of TKE (m^2/s^3)
    P_diss_v                                    , & ! buoyancy production of TKE (m^2/s^3)
    forc                                            ! combined forcing for TKE (m^2/s^3)

  real(wp) :: tke_surf, tke_bott

  real(wp)                                                     :: &
    alpha_tke                                                   , &
    c_eps                                                       , &
    cd                                                          , &
    KappaM_max                                                  , &
    KappaM_min                                                  , &
    KappaH_min                                                  , &
    mxl_min                                                     , &
    c_k                                                         , &
    clc                                                         , &
    tke_surf_min                                                , &
    tke_min

  integer :: tke_mxl_choice

  logical :: only_tke, use_ubound_dirichlet, use_lbound_dirichlet, l_lc, use_Kappa_min

  real(wp)                                                     :: &
    zzw                                                         , & ! depth of interface k
    depth                                                       , & ! total water depth
    diff_surf_forc                                              , &
    diff_bott_forc

  ! input to tri-diagonal solver
  real(wp), dimension(nproma,max_nlev+1)                       :: &
    a_dif                                                       , & !
    b_dif                                                       , & !
    c_dif                                                       , & !
    a_tri                                                       , & !
    b_tri                                                       , & !
    c_tri                                                       , & !
    d_tri                                                       , & !
    ke                                                              !  diffusivity for tke

  real(wp), dimension(nproma,max_nlev+1)                       :: &
    cp                                                          , & !
    dp

  integer :: k, kk, kp1, jc, nlev, i
  real(wp) :: m, fxa

  CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':integrate_tke_gpu'

  CALL set_acc_host_or_device(lzacc, lacc)

  tke_constants_in => tke_constants_saved

  !$ACC DATA COPYIN(tke_constants_in) &
  !$ACC   CREATE(tke_unrest, tke_upd, mxl, sqrttke, prandtl, Rinum, K_diss_v, P_diss_v, forc) &
  !$ACC   CREATE(a_dif, b_dif, c_dif, a_tri, b_tri, c_tri, d_tri, ke, cp, dp) &
  !$ACC   IF(lzacc)

  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  DO jc = start_index, end_index
    IF (levels(jc) > 0) THEN

      ! initialize diagnostics
      tke_Tbpr(jc,:) = 0.0
      tke_Tspr(jc,:) = 0.0
      tke_Tdif(jc,:) = 0.0
      tke_Tdis(jc,:) = 0.0
      tke_Twin(jc,:) = 0.0
      tke_Tiwf(jc,:) = 0.0
      tke_Tbck(jc,:) = 0.0
      tke_Ttot(jc,:) = 0.0

      vmix_int_1(jc,:) = 0.0
      vmix_int_2(jc,:) = 0.0
      vmix_int_3(jc,:) = 0.0

      tke_new(jc,:) = 0.0
      tke_upd(jc,:) = 0.0

      a_dif(jc,:) = 0.0
      b_dif(jc,:) = 0.0
      c_dif(jc,:) = 0.0
      a_tri(jc,:) = 0.0
      b_tri(jc,:) = 0.0
      c_tri(jc,:) = 0.0

      alpha_tke  = tke_constants_in%alpha_tke
      c_eps      = tke_constants_in%c_eps
      cd         = tke_constants_in%cd
      use_Kappa_min = tke_constants_in%use_Kappa_min
      KappaM_min = tke_constants_in%KappaM_min
      KappaH_min = tke_constants_in%KappaH_min
      KappaM_max = tke_constants_in%KappaM_max
      mxl_min    = tke_constants_in%mxl_min
      c_k        = tke_constants_in%c_k
      tke_min    = tke_constants_in%tke_min
      tke_surf_min   = tke_constants_in%tke_surf_min
      tke_mxl_choice = tke_constants_in%tke_mxl_choice
      only_tke = tke_constants_in%only_tke
      l_lc     = tke_constants_in%l_lc
      clc      = tke_constants_in%clc
      use_ubound_dirichlet = tke_constants_in%use_ubound_dirichlet
      use_lbound_dirichlet = tke_constants_in%use_lbound_dirichlet

      nlev = levels(jc)

      !---------------------------------------------------------------------------------
      ! Part 1: calculate mixing length scale
      !---------------------------------------------------------------------------------
      sqrttke(jc,:) = sqrt(max(0d0,tke_old(jc,:)))

      ! turbulent mixing length
      mxl(jc,:) = sqrt(2D0)*sqrttke(jc,:)/sqrt(max(1d-12,Nsqr(jc,:)))

      ! constrain mixing length scale as in MITgcm
      if (tke_mxl_choice==2) then
        mxl(jc,1) = 0.d0
        mxl(jc,nlev+1) = 0.d0
        do k=2,nlev
          mxl(jc,k) = min(mxl(jc,k), mxl(jc,k-1) + dzw(jc,k-1))
        end do
        mxl(jc,nlev) = min(mxl(jc,nlev), mxl_min + dzw(jc,nlev))
        do k=nlev-1,2,-1
          mxl(jc,k) = min(mxl(jc,k), mxl(jc,k+1) + dzw(jc,k))
        end do
        mxl(jc,:) = max(mxl(jc,:), mxl_min)
      else if (tke_mxl_choice==3) then
        depth = sum(dzw(jc,1:nlev))
        do k=2,nlev+1
          zzw = sum(dzw(jc,1:k-1))
          mxl(jc,k) = min(zzw, mxl(jc,k), depth-zzw)
        end do
        mxl(jc,1) = mxl(jc,2)
        mxl(jc,:) = max(mxl(jc,:), mxl_min)
      else
!        CALL finish(method_name,'Wrong choice of tke_mxl_choice. Aborting...')
      end if

      !---------------------------------------------------------------------------------
      ! Part 2: calculate diffusivities
      !---------------------------------------------------------------------------------
      KappaM_out(jc,:) = min(KappaM_max, c_k*mxl(jc,:)*sqrttke(jc,:))
      Rinum(jc,:) = Nsqr(jc,:) / max(Ssqr(jc,:), 1d-12)

      if (.not.only_tke) &
        Rinum(jc,:) = min(Rinum(jc,:), KappaM_out(jc,:) * Nsqr(jc,:) / max(1d-12, alpha_c(jc,:)*E_iw(jc,:)**2))

      prandtl(jc,:) = max(1d0, min(10d0, 6.6*Rinum(jc,:)))
      KappaH_out(jc,:) = KappaM_out(jc,:) / prandtl(jc,:)

      ! restrict to minimum values
      if (use_Kappa_min) then
        KappaM_out(jc,:) = max(KappaM_min, KappaM_out(jc,:))
        KappaH_out(jc,:) = max(KappaH_min, KappaH_out(jc,:))
      end if

      !---------------------------------------------------------------------------------
      ! Part 3: tke forcing
      !---------------------------------------------------------------------------------
      ! initialize forcing
      forc(jc,:) = 0.0

      ! --- forcing by shear and buoycancy production
      K_diss_v(jc,:) = Ssqr(jc,:) * KappaM_out(jc,:)
      P_diss_v(jc,:) = Nsqr(jc,:) * KappaH_out(jc,:)
      P_diss_v(jc,1) = -forc_rho_surf(jc) * grav / rho_ref
      forc(jc,:) = forc(jc,:) + K_diss_v(jc,:) - P_diss_v(jc,:)

      ! --- additional langmuir turbulence term
      if (l_lc) forc(jc,:) = forc(jc,:) + tke_plc(jc,:)

      ! --- forcing by internal wave dissipation
      if (.not.only_tke) forc(jc,:) = forc(jc,:) + iw_diss(jc,:)

      !---------------------------------------------------------------------------------
      ! Part 4: vertical diffusion and dissipation is solved implicitely
      !---------------------------------------------------------------------------------
      ke(jc,:) = 0.d0
      do k = 1, nlev
        kp1      = min(k+1,nlev)
        kk       = max(k,2)
        ke(jc,k) = alpha_tke * 0.5 * (KappaM_out(jc,kp1) + KappaM_out(jc,kk))
      end do

      !--- c is lower diagonal of matrix
      do k=1,nlev
        c_dif(jc,k) = ke(jc,k)/( dzt(jc,k)*dzw(jc,k) )
      end do
      c_dif(jc,nlev+1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

      !--- b is main diagonal of matrix
      do k=2,nlev
        b_dif(jc,k) = ke(jc,k-1)/( dzt(jc,k)*dzw(jc,k-1) ) + ke(jc,k)/( dzt(jc,k)*dzw(jc,k) )
      end do

      !--- a is upper diagonal of matrix
      do k=2,nlev+1
        a_dif(jc,k) = ke(jc,k-1)/( dzt(jc,k)*dzw(jc,k-1) )
      end do
      a_dif(jc,1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

      ! copy tke_old
      tke_upd(jc,1:nlev+1) = tke_old(jc,1:nlev+1)

      ! upper boundary condition
      if (use_ubound_dirichlet) then
        sqrttke(jc,1)   = 0.d0 ! to suppres dissipation for k=1
        forc(jc,1)      = 0.d0 ! to suppres forcing for k=1
        tke_surf        = max(tke_surf_min, cd*forc_tke_surf(jc))
        tke_upd(jc,1)   = tke_surf
        diff_surf_forc  = a_dif(jc,2) * tke_surf
        forc(jc,2)      = forc(jc,2) + diff_surf_forc
        a_dif(jc,2)     = 0.d0 ! and set matrix element to zero
        b_dif(jc,1)     = 0.d0 ! 0 line in matrix for k=1
        c_dif(jc,1)     = 0.d0 ! 0 line in matrix for k=1
      else
        ! add wind forcing
        forc(jc,1)      = forc(jc,1) + (cd*forc_tke_surf(jc)**(3./2.))/(dzt(jc,1))
        b_dif(jc,1)     = ke(jc,1)/( dzt(jc,1)*dzw(jc,1) )
        diff_surf_forc  = 0.0
      endif

      ! lower boundary condition
      if (use_lbound_dirichlet) then
        sqrttke(jc,nlev+1) = 0.d0 ! to suppres dissipation for k=nlev+1
        forc(jc,nlev+1)    = 0.d0 ! to suppres forcing for k=nlev+1
        tke_bott           = tke_min
        tke_upd(jc,nlev+1) = tke_bott
        diff_bott_forc     = c_dif(jc,nlev)*tke_bott
        forc(jc,nlev)      = forc(jc,nlev)+diff_bott_forc
        c_dif(jc,nlev)     = 0.d0 ! and set matrix element to zero
        b_dif(jc,nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
        a_dif(jc,nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
      else
        b_dif(jc,nlev+1)   = ke(jc,nlev)/( dzt(jc,nlev+1)*dzw(jc,nlev) )
        diff_bott_forc     = 0.0
      end if

      !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
      a_tri(jc,:)      = -dtime * a_dif(jc,:)
      b_tri(jc,:)      = 1+dtime * b_dif(jc,:)
      b_tri(jc,2:nlev) = b_tri(jc,2:nlev) + dtime * c_eps * sqrttke(jc,2:nlev) / mxl(jc,2:nlev)
      c_tri(jc,:)      = -dtime * c_dif(jc,:)

      !--- d is r.h.s. of implicite equation (d: new tke with only explicite tendencies included)
      d_tri(jc,1:nlev+1)  = tke_upd(jc,1:nlev+1) + dtime*forc(jc,1:nlev+1)

      ! solve the tri-diag matrix
      cp(jc,1) = c_tri(jc,1) / b_tri(jc,1)
      dp(jc,1) = d_tri(jc,1) / b_tri(jc,1)

      do i = 2,nlev+1
        m = b_tri(jc,i) - cp(jc,i-1) * a_tri(jc,i)
        fxa = 1D0/m
        cp(jc,i) = c_tri(jc,i) * fxa
        dp(jc,i) = (d_tri(jc,i) - dp(jc,i-1) * a_tri(jc,i)) * fxa
      end do

      tke_new(jc,nlev+1) = dp(jc,nlev+1)

      do i = nlev,1,-1
        tke_new(jc,i) = dp(jc,i) - cp(jc,i)* tke_new(jc,i+1)
      end do

      ! --- diagnose implicite tendencies (only for diagnostics)
      ! vertical diffusion of TKE
      do k=2,nlev
        tke_Tdif(jc,k)    = a_dif(jc,k) * tke_new(jc,k-1) - b_dif(jc,k) * tke_new(jc,k) + c_dif(jc,k) * tke_new(jc,k+1)
      end do
      tke_Tdif(jc,1)      = - b_dif(jc,1) * tke_new(jc,1) + c_dif(jc,1) * tke_new(jc,2)
      tke_Tdif(jc,nlev+1) = a_dif(jc,nlev+1) * tke_new(jc,nlev) - b_dif(jc,nlev+1) * tke_new(jc,nlev+1)
      tke_Tdif(jc,2)      = tke_Tdif(jc,2) + diff_surf_forc
      tke_Tdif(jc,nlev)   = tke_Tdif(jc,nlev) + diff_bott_forc

      ! flux out of first box due to diffusion with Dirichlet boundary value of TKE
      ! (tke_surf=tke_upd(1)) and TKE of box below (tke_new(2))
      if (use_ubound_dirichlet) &
        tke_Tdif(jc,1) = - ke(jc,1) / dzw(jc,1) / dzt(jc,1) &
                        * (tke_surf - tke_new(jc,2))

      if (use_lbound_dirichlet) then
        k = nlev+1
        tke_Tdif(jc,k) = ke(jc,k-1)/dzw(jc,k-1)/dzt(jc,k) &
                        * (tke_new(jc,k-1)-tke_bott)
      end if

      ! dissipation of TKE
      tke_Tdis(jc,:) = 0.d0
      tke_Tdis(jc,2:nlev) = -c_eps / mxl(jc,2:nlev) * sqrttke(jc,2:nlev) * tke_new(jc,2:nlev)

      !---------------------------------------------------------------------------------
      ! Part 5: reset tke to bounding values
      !---------------------------------------------------------------------------------
      ! copy of unrestored tke to diagnose energy input by restoring
      tke_unrest(jc,:) = tke_new(jc,:)

      ! restrict values of TKE to tke_min, if IDEMIX is not used
      if (only_tke) &
        tke_new(jc,1:nlev+1) = MAX(tke_new(jc,1:nlev+1), tke_min)

      !---------------------------------------------------------------------------------
      ! Part 6: Assign diagnostic variables
      !---------------------------------------------------------------------------------
      tke_Tbpr(jc,1:nlev+1) = -P_diss_v(jc,1:nlev+1)
      tke_Tspr(jc,1:nlev+1) = K_diss_v(jc,1:nlev+1)
      tke_Tbck(jc,:) = (tke_new(jc,:) - tke_unrest(jc,:)) / dtime
      if (use_ubound_dirichlet) then
        tke_Twin(jc,1) = (tke_new(jc,1) - tke_old(jc,1)) / dtime - tke_Tdif(jc,1)
        tke_Tbck(jc,1) = 0.0
      else
        tke_Twin(jc,1) = (cd * forc_tke_surf(jc)**(3./2.)) / (dzt(jc,1))
      end if

      if (use_lbound_dirichlet) then
        tke_Twin(jc,nlev+1) = (tke_new(jc,nlev+1) - tke_old(jc,nlev+1)) /  dtime - tke_Tdif(jc,nlev+1)
        tke_Tbck(jc,nlev+1) = 0.0
      else
        tke_Twin(jc,nlev+1) = 0.0
      end if

      tke_Tiwf(jc,1:nlev+1) = iw_diss(jc,1:nlev+1)
      tke_Ttot(jc,:)        = (tke_new(jc,:) - tke_old(jc,:)) / dtime
      tke_Lmix(jc,nlev+1:)  = 0.0
      tke_Lmix(jc,1:nlev+1) = mxl(jc,1:nlev+1)
      tke_Pr(jc,nlev+1:)    = 0.0
      tke_Pr(jc,1:nlev+1)   = prandtl(jc,1:nlev+1)

      ! -----------------------------------------------
      ! the rest is for debugging
      ! -----------------------------------------------
      vmix_int_1(jc,:) = KappaH_out(jc,:)
      vmix_int_2(jc,:) = KappaM_out(jc,:)
      vmix_int_3(jc,:) = Nsqr(jc,:)

    END IF
  END DO
  !$ACC END PARALLEL LOOP
  !$ACC WAIT(1)

  !$ACC END DATA

end subroutine integrate_tke_gpu


subroutine integrate_tke_block( &
                         si, &
                         ei, &
                         nproma, &
                         j,                    & ! FIXME: nils: for debuging
                         tstep_count ,         & ! FIXME: nils: for debuging
                         tke_old,              & ! FIXME: nils: today: rename?
                         tke_new,              & ! FIXME: nils: today: rename?
                         KappaM_out,           &
                         KappaH_out,           &
                         vmix_int_1,           & ! FIXME: nils: for debuging
                         vmix_int_2,           & ! FIXME: nils: for debuging
                         vmix_int_3,           & ! FIXME: nils: for debuging
                         dzw,                  &
                         dzt,                  &
                         nlev,                 &
                         max_nlev,             &
                         !old_tke_diss,         & ! FIXME: nils: today: delete?
                         Ssqr,                 &
                         Nsqr,                 &
                         tke_Tbpr,             & ! diagnostic
                         tke_Tspr,             & ! diagnostic
                         tke_Tdif,             & ! diagnostic
                         tke_Tdis,             & ! diagnostic
                         tke_Twin,             & ! diagnostic
                         tke_Tiwf,             & ! diagnostic
                         tke_Tbck,             & ! diagnostic
                         tke_Ttot,             & ! diagnostic
                         !tke,                  &
                         tke_Lmix,             & ! diagnostic
                         tke_Pr,               & ! diagnostic
                         tke_plc,              & ! langmuir turbulence
                         forc_tke_surf,        &
                         E_iw,                 &
                         dtime,                &
                         bottom_fric,          &
                         old_KappaM,           &
                         old_KappaH,           &
                         iw_diss,              & ! FIXME: nils: rename?
                         forc_rho_surf,        &
                         rho_ref,              & ! FIXME: today: put to initialize
                         grav,                 & ! FIXME: today: put to initialize
                         alpha_c,              & ! FIXME: today: put to initialize
                         tke_userdef_constants,&
                         lacc)
!subroutine integrate_tke(jc, blockNo, tstep_count)
!NEC$ always_inline

  type(tke_type), intent(in), optional, target                 :: &
    tke_userdef_constants

  integer, dimension(:), contiguous, intent(in) :: &
    nlev
  integer,intent(in)                                           :: &
    max_nlev                                                      !
  integer,intent(in)                                           :: &
    si, ei, j, nproma, tstep_count

  ! OLD values
  real(wp), contiguous, dimension(:,:), intent(in)                :: &
    tke_old                                                      ,& !
    !old_tke_diss                                                 ,& !
    old_KappaM                                                   ,& !
    old_KappaH                                                   ,& !
    dzt                                                             !
  real(wp), contiguous, dimension(:,:), intent(in)             :: &
    Ssqr                                                         ,& !
    Nsqr

  real(wp), contiguous, dimension(:,:), intent(in)                  :: &
    dzw                                                             !

  ! Langmuir turbulence
  real(wp), contiguous, dimension(:,:), intent(in), optional  :: &
    tke_plc

  ! IDEMIX variables, if run coupled iw_diss is added as forcing to TKE
  real(wp), contiguous, dimension(:,:), intent(in), optional  :: &
    E_iw                                                         ,& !
    alpha_c                                                      ,& !
    iw_diss                                                         !

  real(wp), intent(in)                                   :: &
    rho_ref                                                      ,& !
    dtime                                                        ,& ! time step
    grav                                                            ! gravity constant

  ! FIXME: nils: why inout? why not only in?
  real(wp), contiguous, dimension(:), intent(in)                                   :: &
    bottom_fric                                                  ,& !
    forc_rho_surf                                                ,& ! FIXME: nils: what is this?
    forc_tke_surf

  !real(wp),dimension(max_nlev+1), intent(in), optional       :: &
  !  Kappa_GM                                                        !

  ! NEW values
  real(wp), contiguous, dimension(:,:), intent(inout)             :: &
    tke_new                                                      ,& !
    KappaM_out                                                   ,& !
    KappaH_out

  ! diagnostics
  real(wp), contiguous, dimension(:,:), intent(inout) ::                &
     tke_Tbpr                                                     ,&
     tke_Tspr                                                     ,&
     tke_Tdif                                                     ,&
     tke_Tdis                                                     ,&
     tke_Twin                                                     ,&
     tke_Tiwf                                                     ,&
     tke_Tbck                                                     ,&
     tke_Ttot                                                     ,&
     !tke                                                          ,&
     tke_Lmix                                                     ,&
     tke_Pr                                                       !,&

  real(wp), contiguous, dimension(:,:), intent(inout) ::           &
     vmix_int_1                                                   ,&
     vmix_int_2                                                   ,&
     vmix_int_3

  LOGICAL, OPTIONAL, INTENT(IN) :: lacc

  ! local variables
  real(wp), dimension(nproma,max_nlev+1)                            :: &
    tke_unrest                                                   ,& ! copy of tke before restorring to background value
    tke_upd                                                      ,& ! copy of tke before in which surface/bottom values given by Dirichlet boundary conditions
    mxl                                                          ,& ! mixing length scale (m)
    sqrttke                                                      ,& ! square root of TKE (m/s)
    prandtl                                                      ,& ! Prandtl number
    Rinum                                                        ,& ! Richardson number
    K_diss_v                                                     ,& ! shear production of TKE (m^2/s^3)
    P_diss_v                                                     ,& ! buoyancy production of TKE (m^2/s^3)
    forc                                                            ! combined forcing for TKE (m^2/s^3)

  integer :: &
    nlev_p1(si:ei)

  real(wp), dimension(nproma) :: tke_surf, tke_bott

  real(wp)                                               :: &
    ! FIXME: nils: today: delete tke_surf_corr?
    !tke_surf_corr                                                ,& ! correction of surface density according to surface buoyancy flux
    alpha_tke                                                    ,& ! {30}
    c_eps                                                        ,& ! {0.7}
    cd                                                           ,& ! (3.75}
    KappaM_max                                                   ,& !
    KappaM_min                                                   ,& !
    KappaH_min                                                   ,& !
    mxl_min                                                      ,& ! {1e-8}
    c_k                                                          ,& ! {0.1}
    clc                                                          ,& ! {0.15}
    tke_surf_min                                                 ,& ! {1e-4}
    tke_min                                                         ! {1e-6}
  integer :: tke_mxl_choice

  integer :: max_n

  logical :: only_tke, use_ubound_dirichlet, use_lbound_dirichlet, l_lc, use_Kappa_min

  real(wp), dimension(nproma) :: &
    depth, &
    zzw, &
    diff_surf_forc                                               ,&
    diff_bott_forc

  ! input to tri-diagonal solver
  real(wp), dimension(nproma,max_nlev+1)                            :: &
    a_dif                                                        ,& !
    b_dif                                                        ,& !
    c_dif                                                        ,& !
    a_tri                                                        ,& !
    b_tri                                                        ,& !
    c_tri                                                        ,& !
    d_tri                                                        ,& !
    ke                                                              !  diffusivity for tke

  integer :: k, kk, kp1, jc
  LOGICAL :: lzacc

  type(tke_type), pointer :: tke_constants_in

  CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':integrate_tke'

  CALL set_acc_host_or_device(lzacc, lacc)

  max_n = MAXVAL(nlev(si:ei))

  ! FIXME: nils: What happens here?
  tke_constants_in => tke_constants_saved
  if (present(tke_userdef_constants)) then
   tke_constants_in => tke_userdef_constants
  end if

  ! FIXME: nils: What should we do with height of last grid box dzt(max_nlev+1)?
  !              This should not be as thick as the distance to the next tracer
  !              point (which is a dry point).
  !              Be careful if you divide by 0.5 here. Maybe later we use ddpo
  !              or something like this that might include already this factor.

  !$ACC DATA CREATE(tke_unrest, tke_upd, mxl, sqrttke, prandtl, Rinum, K_diss_v, P_diss_v, forc) &
  !$ACC   CREATE(tke_surf, tke_bott, depth, zzw, diff_surf_forc, diff_bott_forc) &
  !$ACC   CREATE(a_dif, b_dif, c_dif, a_tri, b_tri, c_tri, d_tri, ke, nlev_p1) IF(lzacc)

  ! initialize diagnostics
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  tke_Tbpr = 0.0
  tke_Tspr = 0.0
  tke_Tdif = 0.0
  tke_Tdis = 0.0
  tke_Twin = 0.0
  tke_Tiwf = 0.0
  tke_Tbck = 0.0
  tke_Ttot = 0.0
  vmix_int_1 = 0.0
  vmix_int_2 = 0.0
  vmix_int_3 = 0.0

  tke_new = 0.0
  tke_upd = 0.0

  a_dif = 0.0
  b_dif = 0.0
  c_dif = 0.0
  a_tri = 0.0
  b_tri = 0.0
  c_tri = 0.0

  nlev_p1(si:ei) = nlev(si:ei) + 1

  !$ACC END KERNELS

  !$ACC UPDATE HOST(nlev_p1) WAIT(1)

  !---------------------------------------------------------------------------------
  ! set tke_constants locally
  !---------------------------------------------------------------------------------

  alpha_tke  = tke_constants_in%alpha_tke
  c_eps      = tke_constants_in%c_eps
  cd         = tke_constants_in%cd
  use_Kappa_min = tke_constants_in%use_Kappa_min
  KappaM_min = tke_constants_in%KappaM_min
  KappaH_min = tke_constants_in%KappaH_min
  KappaM_max = tke_constants_in%KappaM_max
  mxl_min    = tke_constants_in%mxl_min
  c_k        = tke_constants_in%c_k
  tke_min    = tke_constants_in%tke_min
  tke_surf_min   = tke_constants_in%tke_surf_min
  tke_mxl_choice = tke_constants_in%tke_mxl_choice
  only_tke = tke_constants_in%only_tke
  l_lc     = tke_constants_in%l_lc
  clc      = tke_constants_in%clc
  use_ubound_dirichlet = tke_constants_in%use_ubound_dirichlet
  use_lbound_dirichlet = tke_constants_in%use_lbound_dirichlet

  !c_k        = 0.1
  !c_eps      = 0.7
  !alpha_tke  = 30.0
  !mxl_min    = 1.d-8
  !KappaM_min = 0.0
  !KappaM_max = 100.0
  !cd         = 3.75
  !tke_min    = 1.d-6
  !tke_mxl_choice = 2
  !tke_surf_min = 1.d-4
  !only_tke = .true.
  !use_ubound_dirichlet = .false.
  !use_lbound_dirichlet = .false.

  !---------------------------------------------------------------------------------
  ! Part 1: calculate mixing length scale
  !---------------------------------------------------------------------------------
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  sqrttke(si:ei,:) = sqrt(max(0d0,tke_old(si:ei,:)))

  ! turbulent mixing length
  mxl(si:ei,:) = sqrt(2D0)*sqrttke(si:ei,:)/sqrt(max(1d-12,Nsqr(si:ei,:)))
  !$ACC END KERNELS

  ! constrain mixing length scale as in MITgcm
  if (tke_mxl_choice==2) then
    !FIXME: What should we do at the surface and bottom?
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    do jc = si, ei
      if (nlev(jc) > 0) then
        mxl(jc,1) = 0.d0
        mxl(jc,nlev(jc)+1) = 0.d0
      endif
    enddo

    !$ACC LOOP SEQ
    do k = 2, max_n
      !$ACC LOOP GANG VECTOR
      do jc = si, ei
        if (k <= nlev(jc)) then
          mxl(jc,k) = min(mxl(jc,k), mxl(jc,k-1)+dzw(jc,k-1))
        endif
      enddo
    enddo

    !$ACC LOOP GANG VECTOR
    do jc = si, ei
      if (nlev(jc) > 0) then
        mxl(jc,nlev(jc)) = min(mxl(jc,nlev(jc)), mxl_min+dzw(jc,nlev(jc)))
      endif
    enddo

    !$ACC LOOP SEQ
    do k = max_n-1, 2, -1
      !$ACC LOOP GANG VECTOR
      do jc = si, ei
        if (k <= nlev(jc) - 1) then
          mxl(jc,k) = min(mxl(jc,k), mxl(jc,k+1)+dzw(jc,k))
        endif
      enddo
    enddo

    !$ACC LOOP SEQ
    do k = 1, max_n + 1
      !$ACC LOOP GANG VECTOR
      do jc = si, ei
        if (k <= nlev(jc) + 1) then
          mxl(jc,k) = max(mxl(jc,k),mxl_min)
        endif
      enddo
    enddo
    !$ACC END PARALLEL

  ! bounded by the distance to surface/bottom
  elseif (tke_mxl_choice==3) then
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    depth(si:ei) = 0.
    zzw(si:ei) = 0.
    !$ACC END KERNELS

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    do k = 1, max_n
      do jc = si, ei
        if (k <= nlev(jc)) then
          depth(jc) = depth(jc) + dzw(jc,k)
        endif
      enddo
    enddo

    !$ACC LOOP SEQ
    do k = 2, max_n + 1
      !$ACC LOOP GANG VECTOR
      do jc = si, ei
        if (k <= nlev(jc) + 1) then
          zzw(jc) = zzw(jc) + dzw(jc,k-1)
          mxl(jc,k) = min(zzw(jc),mxl(jc,k),depth(jc)-zzw(jc))
        endif
      enddo
    enddo
    !$ACC END PARALLEL

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    mxl(si:ei,1) = mxl(si:ei,2)
    mxl(si:ei,:) = max(mxl(si:ei,:), mxl_min)
    !$ACC END KERNELS
  else
!    write(*,*) 'Wrong choice of tke_mxl_choice. Aborting...'
!    stop
    CALL finish(method_name,'Wrong choice of tke_mxl_choice. Aborting...')
  endif

  !---------------------------------------------------------------------------------
  ! Part 2: calculate diffusivities
  !---------------------------------------------------------------------------------
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 1, max_n + 1
    do jc = si, ei
      if (k <= nlev(jc) + 1) then
        KappaM_out(jc,k) = min(KappaM_max,c_k*mxl(jc,k)*sqrttke(jc,k))
        Rinum(jc,k) = Nsqr(jc,k)/max(Ssqr(jc,k),1d-12)
        ! FIXME: nils: Check this later if IDEMIX is coupled.
        ! FIXME: nils: Why E_iw**2 and not dissipation with mixed time level?
        ! FIXME: nils: Why not passing Rinum as Rinum_idemix to tke scheme?
        if (.not.only_tke) then  !IDEMIX is on
          Rinum(jc,k) = min(Rinum(jc,k),KappaM_out(jc,k)*Nsqr(jc,k)/max(1d-12,alpha_c(jc,k)*E_iw(jc,k)**2))
        end if

        prandtl(jc,k)=max(1d0,min(10d0,6.6*Rinum(jc,k)))
        KappaH_out(jc,k)=KappaM_out(jc,k)/prandtl(jc,k)

        ! restrict to minimum values
        if (use_Kappa_min) then
          KappaM_out(jc,k) = max(KappaM_min, KappaM_out(jc,k))
          KappaH_out(jc,k) = max(KappaH_min, KappaH_out(jc,k))
        end if
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP

  !---------------------------------------------------------------------------------
  ! Part 3: tke forcing
  !---------------------------------------------------------------------------------
  ! initialize forcing
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) NO_CREATE(tke_plc) ASYNC(1) IF(lzacc)
  do k = 1, max_n + 1
    do jc = si, ei
      if (k <= nlev(jc) + 1) then

        forc(jc,k) = 0.

        ! --- forcing by shear and buoycancy production
        K_diss_v(jc,k)   = Ssqr(jc,k)*KappaM_out(jc,k)
        if (k > 1) then
          P_diss_v(jc,k)   = Nsqr(jc,k)*KappaH_out(jc,k)
        else
          ! FIXME: nils: Is forc_rho_surf set somewhere?
          ! FIXME: nils: What does forc_rho_surf mean?
          P_diss_v(jc,1) = -forc_rho_surf(jc)*grav/rho_ref
        endif

        forc(jc,k) = K_diss_v(jc,k) - P_diss_v(jc,k)

        ! --- additional langmuir turbulence term
        if (l_lc) then
          forc(jc,k) = forc(jc,k) + tke_plc(jc,k)
        endif

        ! --- forcing by internal wave dissipation
        if (.not.only_tke) then
          forc(jc,k) = forc(jc,k) + iw_diss(jc,k)
        endif
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP


  !---------------------------------------------------------------------------------
  ! Part 4: vertical diffusion and dissipation is solved implicitely
  !---------------------------------------------------------------------------------
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  ke = 0.d0
  !$ACC END KERNELS
  
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP SEQ
  do k = 1, max_n
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    do jc = si, ei
      if (k <= nlev(jc)) then
        kp1 = min(k+1,nlev(jc))
        kk  = max(k,2)
        ke(jc,k) = alpha_tke*0.5*(KappaM_out(jc,kp1)+KappaM_out(jc,kk))
        c_dif(jc,k) = ke(jc,k)/( dzt(jc,k)*dzw(jc,k) )
      endif
    enddo
  enddo

! !--- c is lower diagonal of matrix
! !$ACC LOOP GANG(STATIC:1) VECTOR
! do k = 1, max_n
!   do jc = si, ei
!     if (k <= nlev(jc)) then
!       c_dif(jc,k) = ke(jc,k)/( dzt(jc,k)*dzw(jc,k) )
!     endif
!   enddo
! enddo

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  do jc = si, ei
    if (nlev(jc) > 0) then
      c_dif(jc,nlev(jc)+1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary
    endif
  enddo

  !--- b is main diagonal of matrix
  !$ACC LOOP SEQ
  do k = 2, max_n
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    do jc = si, ei
      if (k <= nlev(jc)) then
        b_dif(jc,k) = ke(jc,k-1)/( dzt(jc,k)*dzw(jc,k-1) ) + ke(jc,k)/( dzt(jc,k)*dzw(jc,k) )
      endif
    enddo
  enddo

  !--- a is upper diagonal of matrix
  !$ACC LOOP SEQ
  do k= 2, max_n + 1
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    do jc = si, ei
      if (k <= nlev(jc) + 1) then
        a_dif(jc,k) = ke(jc,k-1)/( dzt(jc,k)*dzw(jc,k-1) )
      endif
    enddo
  enddo

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  do jc = si, ei
    if (nlev(jc) > 0) then
      a_dif(jc,1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary
    endif
  enddo
  !$ACC END PARALLEL

  ! copy tke_old
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  tke_upd(si:ei,1:max_nlev+1) = tke_old(si:ei,1:max_nlev+1)
  !$ACC END KERNELS

  ! upper boundary condition
  if (use_ubound_dirichlet) then
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        sqrttke(jc,1)      = 0.d0 ! to suppres dissipation for k=1
        forc(jc,1)         = 0.d0 ! to suppres forcing for k=1
        tke_surf(jc)       = max(tke_surf_min, cd*forc_tke_surf(jc))
        tke_upd(jc,1)      = tke_surf(jc)
        ! add diffusive part that depends on tke_surf to forcing
        diff_surf_forc(jc) = a_dif(jc,2)*tke_surf(jc)
        forc(jc,2)         = forc(jc,2)+diff_surf_forc(jc)
        a_dif(jc,2)        = 0.d0 ! and set matrix element to zero
        b_dif(jc,1)        = 0.d0 ! 0 line in matrix for k=1
        c_dif(jc,1)        = 0.d0 ! 0 line in matrix for k=1
      endif
    enddo
    !$ACC END PARALLEL LOOP
  else
    ! add wind forcing
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        forc(jc,1) = forc(jc,1) + (cd*forc_tke_surf(jc)**(3./2.))/(dzt(jc,1))
        b_dif(jc,1)        = ke(jc,1)/( dzt(jc,1)*dzw(jc,1) )
        diff_surf_forc(jc)  = 0.0
      endif
    enddo
    !$ACC END PARALLEL LOOP
  endif

  ! lower boundary condition
  if (use_lbound_dirichlet) then
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        sqrttke(jc,nlev(jc)+1) = 0.d0 ! to suppres dissipation for k=nlev+1
        forc(jc,nlev(jc)+1)    = 0.d0 ! to suppres forcing for k=nlev+1
        ! FIXME: make tke_bott dependend on bottom friction?
        tke_bott(jc)        = tke_min
        !tke_old(nlev+1) = tke_bott
        tke_upd(jc,nlev(jc)+1) = tke_bott(jc)
        ! add diffusive part that depends on tke_bott to forcing
        diff_bott_forc(jc)  = c_dif(jc,nlev(jc))*tke_bott(jc)
        forc(jc,nlev(jc))      = forc(jc,nlev(jc))+diff_bott_forc(jc)
        c_dif(jc,nlev(jc))     = 0.d0 ! and set matrix element to zero
        b_dif(jc,nlev(jc)+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
        a_dif(jc,nlev(jc)+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
      endif
    enddo
    !$ACC END PARALLEL LOOP
  else
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        b_dif(jc,nlev(jc)+1)   = ke(jc,nlev(jc))/( dzt(jc,nlev(jc)+1)*dzw(jc,nlev(jc)) )
        diff_bott_forc(jc)  = 0.0
      endif
    enddo
    !$ACC END PARALLEL LOOP
  endif

  !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  a_tri(si:ei,:) = -dtime*a_dif(si:ei,:)
  b_tri(si:ei,:) = 1+dtime*b_dif(si:ei,:)
  c_tri(si:ei,:) = -dtime*c_dif(si:ei,:)
  !$ACC END KERNELS

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 2, max_n
    do jc = si, ei
      if (k <= nlev(jc)) then
        b_tri(jc,k) = b_tri(jc,k) + dtime*c_eps*sqrttke(jc,k)/mxl(jc,k)
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP

  !--- d is r.h.s. of implicite equation (d: new tke with only explicite tendencies included)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 1, max_n + 1
    do jc = si, ei
      if (k <= nlev(jc) + 1) then
        d_tri(jc,k)  = tke_upd(jc,k) + dtime*forc(jc,k)
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP

  ! solve the tri-diag matrix
  call solve_tridiag_block(a_tri(si:ei,:), b_tri(si:ei,:), c_tri(si:ei,:), d_tri(si:ei,:), &
    tke_new(si:ei,:), nlev_p1, eliminate_upper=.TRUE., lacc=lzacc)

  ! --- diagnose implicite tendencies (only for diagnostics)
  ! vertical diffusion of TKE
  !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  !$ACC LOOP SEQ
  do k=2,max_n
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    do jc = si, ei
      if (k <= nlev(jc)) then
        tke_Tdif(jc,k) = a_dif(jc,k)*tke_new(jc,k-1) - b_dif(jc,k)*tke_new(jc,k) + c_dif(jc,k)*tke_new(jc,k+1)
      endif
    enddo
  enddo

  !$ACC LOOP GANG(STATIC: 1) VECTOR
  do jc = si, ei
    if (nlev(jc) > 0) then
      tke_Tdif(jc,1) = - b_dif(jc,1)*tke_new(jc,1) + c_dif(jc,1)*tke_new(jc,2)
      tke_Tdif(jc,2) = tke_Tdif(jc,2) + diff_surf_forc(jc)
      tke_Tdif(jc,nlev(jc)) = tke_Tdif(jc,nlev(jc)) + diff_bott_forc(jc)
      tke_Tdif(jc,nlev(jc)+1) = a_dif(jc,nlev(jc)+1)*tke_new(jc,nlev(jc)) - b_dif(jc,nlev(jc)+1)*tke_new(jc,nlev(jc)+1)
    endif
  enddo
  !$ACC END PARALLEL

  ! flux out of first box due to diffusion with Dirichlet boundary value of TKE
  ! (tke_surf=tke_upd(1)) and TKE of box below (tke_new(2))
  if (use_ubound_dirichlet) then
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        tke_Tdif(jc,1) = - ke(jc,1)/dzw(jc,1)/dzt(jc,1) &
                      * (tke_surf(jc)-tke_new(jc,2))
      endif
    enddo
    !$ACC END PARALLEL LOOP
  endif
  if (use_lbound_dirichlet) then
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        k = nlev(jc) + 1
        tke_Tdif(jc,k) = ke(jc,k-1)/dzw(jc,k-1)/dzt(jc,k) &
                        * (tke_new(jc,k-1)-tke_bott(jc))
      endif
    enddo
    !$ACC END PARALLEL LOOP
  endif

  ! dissipation of TKE
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  tke_Tdis(si:ei,:) = 0.d0
  !$ACC END KERNELS

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 2, max_n
    do jc = si, ei
      if (k <= nlev(jc)) then
        tke_Tdis(jc,k) = -c_eps/mxl(jc,k)*sqrttke(jc,k)*tke_new(jc,k)
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP
  !tke_diss_out(1:nlev+1) = c_eps/mxl(1:nlev+1)*sqrttke(1:nlev+1)*tke_new(1:nlev)

  !---------------------------------------------------------------------------------
  ! Part 5: reset tke to bounding values
  !---------------------------------------------------------------------------------
  ! copy of unrestored tke to diagnose energy input by restoring
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  tke_unrest(si:ei,:) = tke_new(si:ei,:)
  !$ACC END KERNELS

  ! FIXME: nils: today: delete
  ! add TKE if surface density flux drains TKE in uppermost box
  !tke_surf_corr = 0.0

  ! restrict values of TKE to tke_min, if IDEMIX is not used
  if (only_tke) then
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do k = 1, max_n + 1
      do jc = si, ei
        if (k <= nlev(jc) + 1) then
          tke_new(jc,k) = MAX(tke_new(jc,k), tke_min)
        endif
      enddo
    enddo
    !$ACC END PARALLEL LOOP
  end if

  !---------------------------------------------------------------------------------
  ! Part 6: Assign diagnostic variables
  !---------------------------------------------------------------------------------
  ! tke_Ttot =   tke_Tbpr + tke_Tspr + tke_Tdif + tke_Tdis
  !            + tke_Twin + tke_Tiwf
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 1, max_n + 1
    do jc = si, ei
      if (k <= nlev(jc) + 1) then
        tke_Tbpr(jc,k) = -P_diss_v(jc,k)
        tke_Tspr(jc,k) = K_diss_v(jc,k)
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP

  !tke_Tdif is set above
  !tke_Tdis = -tke_diss_out
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  tke_Tbck(si:ei,:) = (tke_new(si:ei,:)-tke_unrest(si:ei,:))/dtime
  !$ACC END KERNELS
  if (use_ubound_dirichlet) then
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    tke_Twin(si:ei,1) = (tke_new(si:ei,1)-tke_old(si:ei,1))/dtime - tke_Tdif(si:ei,1)
    tke_Tbck(si:ei,1) = 0.0
    !$ACC END KERNELS
  else
    !tke_Twin(1) = forc_tke_surf/(dzt(1))
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    tke_Twin(si:ei,1) = (cd*forc_tke_surf(si:ei)**(3./2.))/(dzt(si:ei,1))
    !$ACC END KERNELS
  endif
  ! FIXME: Find better name for tke_Twin either tke_Tbou or use tke_Tsur
  ! tke_Tbot
  if (use_lbound_dirichlet) then
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        tke_Twin(jc,nlev(jc)+1) = (tke_new(jc,nlev(jc)+1)-tke_old(jc,nlev(jc)+1))/dtime - tke_Tdif(jc,nlev(jc)+1)
        tke_Tbck(jc,nlev(jc)+1) = 0.0
      endif
    enddo
    !$ACC END PARALLEL LOOP
  else
    !FIXME: no flux condition so far, add bottom friction later
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do jc = si, ei
      if (nlev(jc) > 0) then
        tke_Twin(jc,nlev(jc)+1) = 0.0
      endif
    enddo
    !$ACC END PARALLEL LOOP
  endif

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 1, max_n + 1
    do jc = si, ei
      if (k <= nlev(jc) + 1) then
        tke_Tiwf(jc,k) = iw_diss(jc,k)
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP

  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  tke_Ttot(si:ei,:) = (tke_new(si:ei,:)-tke_old(si:ei,:))/dtime
  !$ACC END KERNELS
  !tke = tke_new

  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  do k = 1, max_n + 1
    do jc = si, ei
      if (k <= nlev(jc) + 1) then
        tke_Lmix(jc,k) = mxl(jc,k)
        tke_Pr(jc,k) = prandtl(jc,k)
      else
        tke_Pr(jc,k) = 0.0
        tke_Lmix(jc,k) = 0.0
      endif
    enddo
  enddo
  !$ACC END PARALLEL LOOP

  ! -----------------------------------------------
  ! the rest is for debugging
  ! -----------------------------------------------
  !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  vmix_int_1(si:ei,:) = KappaH_out(si:ei,:)
  vmix_int_2(si:ei,:) = KappaM_out(si:ei,:)
  vmix_int_3(si:ei,:) = Nsqr(si:ei,:)
  !$ACC END KERNELS
  !vmix_int_1 = forc
  !vmix_int_2 = Nsqr
  !vmix_int_3 = Ssqr

!  if (.false.) then
!    write(*,*) 'i = ', i, 'j = ', j, 'tstep_count = ', tstep_count
!  if (i==8 .and. j==10) then
  !if (i==45 .and. j==10 .and. tstep_count==10) then
  ! -----------------------------------------------

!    write(*,*) '================================================================================'
!    write(*,*) 'i = ', i, 'j = ', j, 'tstep_count = ', tstep_count
!!!    write(*,*) 'nlev = ', nlev
!!!    write(*,*) 'dtime = ', dtime
!!!    write(*,*) 'dzt = ', dzt
!!!    write(*,*) 'dzw = ', dzw
!!!    write(*,*) 'Nsqr = ', Nsqr
!!!    write(*,*) 'Ssqr = ', Ssqr
!!!    !write(*,*) 'tho = ', tho(i,j,1:nlev)
!!!    !write(*,*) 'sao = ', sao(i,j,1:nlev)
!!!    !write(*,*) 'bottom_fric = ', bottom_fric
!!!    !write(*,*) 'forc_tke_surf = ', forc_tke_surf
!    write(*,*) 'sqrttke = ', sqrttke
!!!    write(*,*) 'mxl = ', mxl
!!!    write(*,*) 'KappaM_out = ', KappaM_out
!!!    write(*,*) 'KappaH_out = ', KappaH_out
!!!    write(*,*) 'forc = ', forc
!!!    !write(*,*) 'Rinum = ', Rinum
!!!    write(*,*) 'prandtl = ', prandtl
!!!    !write(*,*) 'checkpoint d_tri'
!!!    !write(*,*) 'K_diss_v = ', K_diss_v
!!!    !write(*,*) 'P_diss_v = ', P_diss_v
!!!    !write(*,*) 'delta = ', delta
!!!    write(*,*) 'ke = ', ke
!!!    write(*,*) 'a_tri = ', a_tri
!!!    write(*,*) 'b_tri = ', b_tri
!!!    write(*,*) 'c_tri = ', c_tri
!!!    write(*,*) 'd_tri = ', d_tri
!!!    !write(*,*) 'tke_old = ', tke_old
!!!    write(*,*) 'tke_new = ', tke_new
!!!    write(*,*) 'tke_Tbpr = ', tke_Tbpr
!!!    write(*,*) 'tke_Tspr = ', tke_Tspr
!!!    write(*,*) 'tke_Tdif = ', tke_Tdif
!!!    write(*,*) 'tke_Tdis = ', tke_Tdis
!!!    write(*,*) 'tke_Twin = ', tke_Twin
!!!    write(*,*) 'tke_Tiwf = ', tke_Tiwf
!!!    write(*,*) 'tke_Ttot = ', tke_Ttot
!!!    write(*,*) 'tke_Ttot - tke_Tsum = ', &
!!!      tke_Ttot-(tke_Tbpr+tke_Tspr+tke_Tdif+tke_Tdis+tke_Twin+tke_Tiwf)
!!!    !write(*,*) 'dzw = ', dzw
!!!    !write(*,*) 'dzt = ', dzt
!!!    ! FIXME: partial bottom cells!!
!!!    ! namelist parameters
!!!    write(*,*) 'c_k = ', c_k
!!!    write(*,*) 'c_eps = ', c_eps
!!!    write(*,*) 'alpha_tke = ', alpha_tke
!!!    write(*,*) 'mxl_min = ', mxl_min
!!!    write(*,*) 'KappaM_min = ', KappaM_min
!!!    write(*,*) 'KappaM_max = ', KappaM_max
!!!    ! FIXME: Make tke_mxl_choice available!
!!!    !write(*,*) 'tke_mxl_choice = ', tke_mxl_choice
!!!    !write(*,*) 'cd = ', cd
!!!    write(*,*) 'tke_min = ', tke_min
!!!    write(*,*) 'tke_surf_min = ', tke_surf_min
!!!    write(*,*) 'only_tke = ', only_tke
!!!    write(*,*) 'use_ubound_dirichlet = ', use_ubound_dirichlet
!!!    write(*,*) 'use_lbound_dirichlet = ', use_lbound_dirichlet
!!!    !write(*,*) 'tke(nlev) = ', tke(nlev), 'tke(nlev+1) = ', tke(nlev+1)
!!!    !write(*,*) 'tke(nlev+2) = ', tke(nlev+2)
!    write(*,*) '================================================================================'
  !end if
  !if (i==45 .and. j==10 .and. tstep_count==10) then
!    stop
!  end if ! if (i==, j==, tstep==)
!  end if ! if (.true./.false.)

  !$ACC WAIT(1)
  !$ACC END DATA
end subroutine integrate_tke_block


!=================================================================================

subroutine integrate_tke( &
                         i,                    & ! FIXME: nils: for debuging
                         j,                    & ! FIXME: nils: for debuging
                         tstep_count ,         & ! FIXME: nils: for debuging
                         !tke_diss_out,         & ! FIXME: nils: today: delete?
                         tke_old,              & ! FIXME: nils: today: rename?
                         tke_new,              & ! FIXME: nils: today: rename?
                         KappaM_out,           &
                         KappaH_out,           &
                         vmix_int_1,           & ! FIXME: nils: for debuging
                         vmix_int_2,           & ! FIXME: nils: for debuging
                         vmix_int_3,           & ! FIXME: nils: for debuging
                         dzw,                  &
                         dzt,                  &
                         nlev,                 &
                         max_nlev,             &
                         !old_tke_diss,         & ! FIXME: nils: today: delete?
                         Ssqr,                 &
                         Nsqr,                 &
                         tke_Tbpr,             & ! diagnostic
                         tke_Tspr,             & ! diagnostic
                         tke_Tdif,             & ! diagnostic
                         tke_Tdis,             & ! diagnostic
                         tke_Twin,             & ! diagnostic
                         tke_Tiwf,             & ! diagnostic
                         tke_Tbck,             & ! diagnostic
                         tke_Ttot,             & ! diagnostic
                         !tke,                  &
                         tke_Lmix,             & ! diagnostic
                         tke_Pr,               & ! diagnostic
                         tke_plc,              & ! langmuir turbulence
                         forc_tke_surf,        &
                         E_iw,                 &
                         dtime,                &
                         bottom_fric,          &
                         iw_diss,              & ! FIXME: nils: rename?
                         forc_rho_surf,        &
                         !Kappa_GM,             & ! FIXME: nils: today: delete?
                         rho_ref,              & ! FIXME: today: put to initialize
                         grav,                 & ! FIXME: today: put to initialize
                         alpha_c,              & ! FIXME: today: put to initialize
                         tke_userdef_constants)
!subroutine integrate_tke(jc, blockNo, tstep_count)
!NEC$ always_inline

  type(tke_type), intent(in), optional, target                 :: &
    tke_userdef_constants

  integer,intent(in)                                           :: &
    nlev                                                         ,& !
    max_nlev                                                      !
  integer,intent(in)                                           :: &
    i, j, tstep_count

  ! OLD values
  real(wp), dimension(max_nlev+1), intent(in)                  :: &
    tke_old                                                      ,& !
    !old_tke_diss                                                 ,& !
    dzt                                                             !
  real(wp), dimension(max_nlev+1), intent(in)                  :: &
    Ssqr                                                         ,& !
    Nsqr

  real(wp), dimension(max_nlev), intent(in)                    :: &
    dzw                                                             !

  ! Langmuir turbulence
  real(wp), dimension(max_nlev+1), intent(in), optional        :: &
    tke_plc

  ! IDEMIX variables, if run coupled iw_diss is added as forcing to TKE
  real(wp), dimension(max_nlev+1), intent(in), optional        :: &
    E_iw                                                         ,& !
    alpha_c                                                      ,& !
    iw_diss                                                         !

  real(wp), intent(in)                                         :: &
    bottom_fric                                                  ,& !
    forc_rho_surf                                                ,& ! FIXME: nils: what is this?
    rho_ref                                                      ,& !
    dtime                                                        ,& ! time step
    grav                                                            ! gravity constant

  ! FIXME: nils: why inout? why not only in?
  real(wp), intent(in)                                         :: &
    forc_tke_surf

  !real(wp),dimension(max_nlev+1), intent(in), optional           :: &
  !  Kappa_GM                                                        !

  ! NEW values
  real(wp), dimension(max_nlev+1), intent(out)                 :: &
    tke_new                                                      ,& !
    !tke_diss_out                                                 ,& !
    KappaM_out                                                   ,& !
    KappaH_out

  ! diagnostics
  real(wp), dimension(max_nlev+1), intent(out) ::                  &
     tke_Tbpr                                                     ,&
     tke_Tspr                                                     ,&
     tke_Tdif                                                     ,&
     tke_Tdis                                                     ,&
     tke_Twin                                                     ,&
     tke_Tiwf                                                     ,&
     tke_Tbck                                                     ,&
     tke_Ttot                                                     ,&
     !tke                                                          ,&
     tke_Lmix                                                     ,&
     tke_Pr                                                       !,&

  real(wp), dimension(max_nlev+1), intent(out) ::                  &
     vmix_int_1                                                   ,&
     vmix_int_2                                                   ,&
     vmix_int_3

  ! local variables
  real(wp), dimension(max_nlev+1)                              :: &
    tke_unrest                                                   ,& ! copy of tke before restorring to background value
    tke_upd                                                      ,& ! copy of tke before in which surface/bottom values given by Dirichlet boundary conditions
    mxl                                                          ,& ! mixing length scale (m)
    sqrttke                                                      ,& ! square root of TKE (m/s)
    prandtl                                                      ,& ! Prandtl number
    Rinum                                                        ,& ! Richardson number
    K_diss_v                                                     ,& ! shear production of TKE (m^2/s^3)
    P_diss_v                                                     ,& ! buoyancy production of TKE (m^2/s^3)
    forc                                                            ! combined forcing for TKE (m^2/s^3)

  real(wp) :: tke_surf, tke_bott

  real(wp)                                                     :: &
    ! FIXME: nils: today: delete tke_surf_corr?
    !tke_surf_corr                                                ,& ! correction of surface density according to surface buoyancy flux
    alpha_tke                                                    ,& ! {30}
    c_eps                                                        ,& ! {0.7}
    cd                                                           ,& ! (3.75}
    KappaM_max                                                   ,& !
    KappaM_min                                                   ,& !
    KappaH_min                                                   ,& !
    mxl_min                                                      ,& ! {1e-8}
    c_k                                                          ,& ! {0.1}
    clc                                                          ,& ! {0.15}
    tke_surf_min                                                 ,& ! {1e-4}
    tke_min                                                         ! {1e-6}
  integer :: tke_mxl_choice

  logical :: only_tke, use_ubound_dirichlet, use_lbound_dirichlet, l_lc, use_Kappa_min

  real(wp)                                                     :: &
    zzw                                                          ,& ! depth of interface k
    depth                                                        ,&! total water depth
    diff_surf_forc                                               ,&
    diff_bott_forc

  ! input to tri-diagonal solver
  real(wp), dimension(max_nlev+1)                              :: &
    a_dif                                                        ,& !
    b_dif                                                        ,& !
    c_dif                                                        ,& !
    a_tri                                                        ,& !
    b_tri                                                        ,& !
    c_tri                                                        ,& !
    d_tri                                                        ,& !
    ke                                                              !  diffusivity for tke

  integer :: k, kk, kp1

  type(tke_type), pointer :: tke_constants_in

  CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':integrate_tke'

  ! FIXME: nils: What happens here?
  tke_constants_in => tke_constants_saved
  if (present(tke_userdef_constants)) then
   tke_constants_in => tke_userdef_constants
  end if

  ! FIXME: nils: What should we do with height of last grid box dzt(max_nlev+1)?
  !              This should not be as thick as the distance to the next tracer
  !              point (which is a dry point).
  !              Be careful if you divide by 0.5 here. Maybe later we use ddpo
  !              or something like this that might include already this factor.

  ! initialize diagnostics
  tke_Tbpr = 0.0
  tke_Tspr = 0.0
  tke_Tdif = 0.0
  tke_Tdis = 0.0
  tke_Twin = 0.0
  tke_Tiwf = 0.0
  tke_Tbck = 0.0
  tke_Ttot = 0.0

  vmix_int_1 = 0.0
  vmix_int_2 = 0.0
  vmix_int_3 = 0.0

  tke_new = 0.0
  tke_upd = 0.0

  a_dif = 0.0
  b_dif = 0.0
  c_dif = 0.0
  a_tri = 0.0
  b_tri = 0.0
  c_tri = 0.0

  !---------------------------------------------------------------------------------
  ! set tke_constants locally
  !---------------------------------------------------------------------------------

  alpha_tke  = tke_constants_in%alpha_tke
  c_eps      = tke_constants_in%c_eps
  cd         = tke_constants_in%cd
  use_Kappa_min = tke_constants_in%use_Kappa_min
  KappaM_min = tke_constants_in%KappaM_min
  KappaH_min = tke_constants_in%KappaH_min
  KappaM_max = tke_constants_in%KappaM_max
  mxl_min    = tke_constants_in%mxl_min
  c_k        = tke_constants_in%c_k
  tke_min    = tke_constants_in%tke_min
  tke_surf_min   = tke_constants_in%tke_surf_min
  tke_mxl_choice = tke_constants_in%tke_mxl_choice
  only_tke = tke_constants_in%only_tke
  l_lc     = tke_constants_in%l_lc
  clc      = tke_constants_in%clc
  use_ubound_dirichlet = tke_constants_in%use_ubound_dirichlet
  use_lbound_dirichlet = tke_constants_in%use_lbound_dirichlet

  !c_k        = 0.1
  !c_eps      = 0.7
  !alpha_tke  = 30.0
  !mxl_min    = 1.d-8
  !KappaM_min = 0.0
  !KappaM_max = 100.0
  !cd         = 3.75
  !tke_min    = 1.d-6
  !tke_mxl_choice = 2
  !tke_surf_min = 1.d-4
  !only_tke = .true.
  !use_ubound_dirichlet = .false.
  !use_lbound_dirichlet = .false.

  !---------------------------------------------------------------------------------
  ! Part 1: calculate mixing length scale
  !---------------------------------------------------------------------------------
  sqrttke = sqrt(max(0d0,tke_old))

  ! turbulent mixing length
  mxl = sqrt(2D0)*sqrttke/sqrt(max(1d-12,Nsqr))

  ! constrain mixing length scale as in MITgcm
  if (tke_mxl_choice==2) then
    !FIXME: What should we do at the surface and bottom?
    mxl(1) = 0.d0
    mxl(nlev+1) = 0.d0
    do k=2,nlev
      mxl(k) = min(mxl(k), mxl(k-1)+dzw(k-1))
    enddo
    mxl(nlev) = min(mxl(nlev), mxl_min+dzw(nlev))
    do k=nlev-1,2,-1
      mxl(k) = min(mxl(k), mxl(k+1)+dzw(k))
    enddo
    mxl= max(mxl,mxl_min)
  ! bounded by the distance to surface/bottom
  elseif (tke_mxl_choice==3) then
    depth = sum(dzw(1:nlev))
    do k=2,nlev+1
     zzw = sum(dzw(1:k-1))
     mxl(k) = min(zzw,mxl(k),depth-zzw)
    enddo
    mxl(1) = mxl(2)
    mxl= max(mxl,mxl_min)
  else
!    write(*,*) 'Wrong choice of tke_mxl_choice. Aborting...'
!    stop
    CALL finish(method_name,'Wrong choice of tke_mxl_choice. Aborting...')
  endif

  !---------------------------------------------------------------------------------
  ! Part 2: calculate diffusivities
  !---------------------------------------------------------------------------------
  KappaM_out = min(KappaM_max,c_k*mxl*sqrttke)
  Rinum = Nsqr/max(Ssqr,1d-12)

  ! FIXME: nils: Check this later if IDEMIX is coupled.
  ! FIXME: nils: Why E_iw**2 and not dissipation with mixed time level?
  ! FIXME: nils: Why not passing Rinum as Rinum_idemix to tke scheme?
  if (.not.only_tke) then  !IDEMIX is on
    Rinum = min(Rinum,KappaM_out*Nsqr/max(1d-12,alpha_c*E_iw**2))
  end if

  prandtl=max(1d0,min(10d0,6.6*Rinum))
  KappaH_out=KappaM_out/prandtl

  ! restrict to minimum values
  if (use_Kappa_min) then
    KappaM_out = max(KappaM_min, KappaM_out)
    KappaH_out = max(KappaH_min, KappaH_out)
  end if

  !---------------------------------------------------------------------------------
  ! Part 3: tke forcing
  !---------------------------------------------------------------------------------
  ! initialize forcing
  forc = 0.0

  ! --- forcing by shear and buoycancy production
  K_diss_v   = Ssqr*KappaM_out
  P_diss_v   = Nsqr*KappaH_out
  ! FIXME: nils: Is forc_rho_surf set somewhere?
  ! FIXME: nils: What does forc_rho_surf mean?
  P_diss_v(1) = -forc_rho_surf*grav/rho_ref
  forc = forc + K_diss_v - P_diss_v

  ! --- additional langmuir turbulence term
  if (l_lc) then
    forc = forc + tke_plc
  endif

  ! --- forcing by internal wave dissipation
  if (.not.only_tke) then
    forc = forc + iw_diss
  endif

  !---------------------------------------------------------------------------------
  ! Part 4: vertical diffusion and dissipation is solved implicitely
  !---------------------------------------------------------------------------------
  ke = 0.d0
  do k = 1, nlev
    kp1 = min(k+1,nlev)
    kk  = max(k,2)
    ke(k) = alpha_tke*0.5*(KappaM_out(kp1)+KappaM_out(kk))
  enddo

  !--- c is lower diagonal of matrix
  do k=1,nlev
    !c_dif(k) = delta(k+1)/dzt(k)
    c_dif(k) = ke(k)/( dzt(k)*dzw(k) )
  enddo
  c_dif(nlev+1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

  !--- b is main diagonal of matrix
  do k=2,nlev
    !b_dif(k) = delta(k)/dzt(k)+delta(k+1)/dzt(k)
    b_dif(k) = ke(k-1)/( dzt(k)*dzw(k-1) ) + ke(k)/( dzt(k)*dzw(k) )
  enddo

  !--- a is upper diagonal of matrix
  do k=2,nlev+1
    !a_dif(k) = delta(k)/dzt(k)
    a_dif(k) = ke(k-1)/( dzt(k)*dzw(k-1) )
  enddo
  a_dif(1) = 0.d0 ! not part of the diffusion matrix, thus value is arbitrary

  ! copy tke_old
  tke_upd(1:nlev+1) = tke_old(1:nlev+1)

  ! upper boundary condition
  if (use_ubound_dirichlet) then
    sqrttke(1)      = 0.d0 ! to suppres dissipation for k=1
    forc(1)         = 0.d0 ! to suppres forcing for k=1
    tke_surf        = max(tke_surf_min, cd*forc_tke_surf)
    tke_upd(1)      = tke_surf
    ! add diffusive part that depends on tke_surf to forcing
    diff_surf_forc  = a_dif(2)*tke_surf
    forc(2)         = forc(2)+diff_surf_forc
    a_dif(2)        = 0.d0 ! and set matrix element to zero
    b_dif(1)        = 0.d0 ! 0 line in matrix for k=1
    c_dif(1)        = 0.d0 ! 0 line in matrix for k=1
  else
    ! add wind forcing
    forc(1) = forc(1) + (cd*forc_tke_surf**(3./2.))/(dzt(1))
    b_dif(1)        = ke(1)/( dzt(1)*dzw(1) )
    diff_surf_forc  = 0.0
  endif

  ! lower boundary condition
  if (use_lbound_dirichlet) then
    sqrttke(nlev+1) = 0.d0 ! to suppres dissipation for k=nlev+1
    forc(nlev+1)    = 0.d0 ! to suppres forcing for k=nlev+1
    ! FIXME: make tke_bott dependend on bottom friction?
    tke_bott        = tke_min
    !tke_old(nlev+1) = tke_bott
    tke_upd(nlev+1) = tke_bott
    ! add diffusive part that depends on tke_bott to forcing
    diff_bott_forc  = c_dif(nlev)*tke_bott
    forc(nlev)      = forc(nlev)+diff_bott_forc
    c_dif(nlev)     = 0.d0 ! and set matrix element to zero
    b_dif(nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
    a_dif(nlev+1)   = 0.d0 ! 0 line in matrix for k=nlev+1
  else
    b_dif(nlev+1)   = ke(nlev)/( dzt(nlev+1)*dzw(nlev) )
    diff_bott_forc  = 0.0
  endif

  !--- construct tridiagonal matrix to solve diffusion and dissipation implicitely
  a_tri = -dtime*a_dif
  b_tri = 1+dtime*b_dif
  b_tri(2:nlev) = b_tri(2:nlev) + dtime*c_eps*sqrttke(2:nlev)/mxl(2:nlev)
  c_tri = -dtime*c_dif

  !--- d is r.h.s. of implicite equation (d: new tke with only explicite tendencies included)
  d_tri(1:nlev+1)  = tke_upd(1:nlev+1) + dtime*forc(1:nlev+1)

  ! solve the tri-diag matrix
  call solve_tridiag(a_tri, b_tri, c_tri, d_tri, tke_new, nlev+1)

  ! --- diagnose implicite tendencies (only for diagnostics)
  ! vertical diffusion of TKE
  do k=2,nlev
    tke_Tdif(k) = a_dif(k)*tke_new(k-1) - b_dif(k)*tke_new(k) + c_dif(k)*tke_new(k+1)
  enddo
  tke_Tdif(1) = - b_dif(1)*tke_new(1) + c_dif(1)*tke_new(2)
  tke_Tdif(nlev+1) = a_dif(nlev+1)*tke_new(nlev) - b_dif(nlev+1)*tke_new(nlev+1)
  tke_Tdif(2) = tke_Tdif(2) + diff_surf_forc
  tke_Tdif(nlev) = tke_Tdif(nlev) + diff_bott_forc

  ! flux out of first box due to diffusion with Dirichlet boundary value of TKE
  ! (tke_surf=tke_upd(1)) and TKE of box below (tke_new(2))
  if (use_ubound_dirichlet) then
    tke_Tdif(1) = - ke(1)/dzw(1)/dzt(1) &
                    * (tke_surf-tke_new(2))
  endif
  if (use_lbound_dirichlet) then
    k = nlev+1
    tke_Tdif(k) = ke(k-1)/dzw(k-1)/dzt(k) &
                    * (tke_new(k-1)-tke_bott)
  endif

  ! dissipation of TKE
  tke_Tdis = 0.d0
  tke_Tdis(2:nlev) = -c_eps/mxl(2:nlev)*sqrttke(2:nlev)*tke_new(2:nlev)
  !tke_diss_out(1:nlev+1) = c_eps/mxl(1:nlev+1)*sqrttke(1:nlev+1)*tke_new(1:nlev)

  !---------------------------------------------------------------------------------
  ! Part 5: reset tke to bounding values
  !---------------------------------------------------------------------------------
  ! copy of unrestored tke to diagnose energy input by restoring
  tke_unrest = tke_new

  ! FIXME: nils: today: delete
  ! add TKE if surface density flux drains TKE in uppermost box
  !tke_surf_corr = 0.0

  ! restrict values of TKE to tke_min, if IDEMIX is not used
  if (only_tke) then
    tke_new(1:nlev+1) = MAX(tke_new(1:nlev+1), tke_min)
  end if

  !---------------------------------------------------------------------------------
  ! Part 6: Assign diagnostic variables
  !---------------------------------------------------------------------------------
  ! tke_Ttot =   tke_Tbpr + tke_Tspr + tke_Tdif + tke_Tdis
  !            + tke_Twin + tke_Tiwf
  tke_Tbpr(1:nlev+1) = -P_diss_v(1:nlev+1)
  tke_Tspr(1:nlev+1) = K_diss_v(1:nlev+1)
  !tke_Tdif is set above
  !tke_Tdis = -tke_diss_out
  tke_Tbck = (tke_new-tke_unrest)/dtime
  if (use_ubound_dirichlet) then
    tke_Twin(1) = (tke_new(1)-tke_old(1))/dtime - tke_Tdif(1)
    tke_Tbck(1) = 0.0
  else
    !tke_Twin(1) = forc_tke_surf/(dzt(1))
    tke_Twin(1) = (cd*forc_tke_surf**(3./2.))/(dzt(1))
  endif
  ! FIXME: Find better name for tke_Twin either tke_Tbou or use tke_Tsur
  ! tke_Tbot
  if (use_lbound_dirichlet) then
    tke_Twin(nlev+1) = (tke_new(nlev+1)-tke_old(nlev+1))/dtime - tke_Tdif(nlev+1)
    tke_Tbck(nlev+1) = 0.0
  else
    !FIXME: no flux condition so far, add bottom friction later
    tke_Twin(nlev+1) = 0.0
  endif

  tke_Tiwf(1:nlev+1) = iw_diss(1:nlev+1)
  tke_Ttot = (tke_new-tke_old)/dtime
  !tke = tke_new
  tke_Lmix(nlev+1:) = 0.0
  tke_Lmix(1:nlev+1) = mxl(1:nlev+1)
  tke_Pr(nlev+1:) = 0.0
  tke_Pr(1:nlev+1) = prandtl(1:nlev+1)

  ! -----------------------------------------------
  ! the rest is for debugging
  ! -----------------------------------------------
  vmix_int_1 = KappaH_out
  vmix_int_2 = KappaM_out
  vmix_int_3 = Nsqr
  !vmix_int_1 = forc
  !vmix_int_2 = Nsqr
  !vmix_int_3 = Ssqr

!  if (.false.) then
!    write(*,*) 'i = ', i, 'j = ', j, 'tstep_count = ', tstep_count
!  if (i==8 .and. j==10) then
  !if (i==45 .and. j==10 .and. tstep_count==10) then
  ! -----------------------------------------------

!    write(*,*) '================================================================================'
!    write(*,*) 'i = ', i, 'j = ', j, 'tstep_count = ', tstep_count
!!!    write(*,*) 'nlev = ', nlev
!!!    write(*,*) 'dtime = ', dtime
!!!    write(*,*) 'dzt = ', dzt
!!!    write(*,*) 'dzw = ', dzw
!!!    write(*,*) 'Nsqr = ', Nsqr
!!!    write(*,*) 'Ssqr = ', Ssqr
!!!    !write(*,*) 'tho = ', tho(i,j,1:nlev)
!!!    !write(*,*) 'sao = ', sao(i,j,1:nlev)
!!!    !write(*,*) 'bottom_fric = ', bottom_fric
!!!    !write(*,*) 'forc_tke_surf = ', forc_tke_surf
!    write(*,*) 'sqrttke = ', sqrttke
!!!    write(*,*) 'mxl = ', mxl
!!!    write(*,*) 'KappaM_out = ', KappaM_out
!!!    write(*,*) 'KappaH_out = ', KappaH_out
!!!    write(*,*) 'forc = ', forc
!!!    !write(*,*) 'Rinum = ', Rinum
!!!    write(*,*) 'prandtl = ', prandtl
!!!    !write(*,*) 'checkpoint d_tri'
!!!    !write(*,*) 'K_diss_v = ', K_diss_v
!!!    !write(*,*) 'P_diss_v = ', P_diss_v
!!!    !write(*,*) 'delta = ', delta
!!!    write(*,*) 'ke = ', ke
!!!    write(*,*) 'a_tri = ', a_tri
!!!    write(*,*) 'b_tri = ', b_tri
!!!    write(*,*) 'c_tri = ', c_tri
!!!    write(*,*) 'd_tri = ', d_tri
!!!    !write(*,*) 'tke_old = ', tke_old
!!!    write(*,*) 'tke_new = ', tke_new
!!!    write(*,*) 'tke_Tbpr = ', tke_Tbpr
!!!    write(*,*) 'tke_Tspr = ', tke_Tspr
!!!    write(*,*) 'tke_Tdif = ', tke_Tdif
!!!    write(*,*) 'tke_Tdis = ', tke_Tdis
!!!    write(*,*) 'tke_Twin = ', tke_Twin
!!!    write(*,*) 'tke_Tiwf = ', tke_Tiwf
!!!    write(*,*) 'tke_Ttot = ', tke_Ttot
!!!    write(*,*) 'tke_Ttot - tke_Tsum = ', &
!!!      tke_Ttot-(tke_Tbpr+tke_Tspr+tke_Tdif+tke_Tdis+tke_Twin+tke_Tiwf)
!!!    !write(*,*) 'dzw = ', dzw
!!!    !write(*,*) 'dzt = ', dzt
!!!    ! FIXME: partial bottom cells!!
!!!    ! namelist parameters
!!!    write(*,*) 'c_k = ', c_k
!!!    write(*,*) 'c_eps = ', c_eps
!!!    write(*,*) 'alpha_tke = ', alpha_tke
!!!    write(*,*) 'mxl_min = ', mxl_min
!!!    write(*,*) 'KappaM_min = ', KappaM_min
!!!    write(*,*) 'KappaM_max = ', KappaM_max
!!!    ! FIXME: Make tke_mxl_choice available!
!!!    !write(*,*) 'tke_mxl_choice = ', tke_mxl_choice
!!!    !write(*,*) 'cd = ', cd
!!!    write(*,*) 'tke_min = ', tke_min
!!!    write(*,*) 'tke_surf_min = ', tke_surf_min
!!!    write(*,*) 'only_tke = ', only_tke
!!!    write(*,*) 'use_ubound_dirichlet = ', use_ubound_dirichlet
!!!    write(*,*) 'use_lbound_dirichlet = ', use_lbound_dirichlet
!!!    !write(*,*) 'tke(nlev) = ', tke(nlev), 'tke(nlev+1) = ', tke(nlev+1)
!!!    !write(*,*) 'tke(nlev+2) = ', tke(nlev+2)
!    write(*,*) '================================================================================'
  !end if
  !if (i==45 .and. j==10 .and. tstep_count==10) then
!    stop
!  end if ! if (i==, j==, tstep==)
!  end if ! if (.true./.false.)
end subroutine integrate_tke

!=================================================================================

subroutine tke_put_tke_int(varname,val,tke_userdef_constants)
!This subroutine puts integer values to TKE variables
!IN
    character(len=*),           intent(in) :: varname
    integer,                    intent(in) :: val
!OUT
    type(tke_type), intent(inout), target, optional:: tke_userdef_constants
    type(tke_type), pointer :: tke_constants_out

  tke_constants_out=>tke_constants_saved
  if (present(tke_userdef_constants)) then
  tke_constants_out=> tke_userdef_constants
 end if

  select case(trim(varname))

    case('handle_old_vals')
    tke_constants_out%handle_old_vals=val
    case ('tke_mxl_choice')
    tke_constants_out%tke_mxl_choice=val
  end select

end subroutine tke_put_tke_int

!=================================================================================

subroutine tke_put_tke_logical(varname,val,tke_userdef_constants)
!This subroutine puts logicals to TKE variables
!IN
    character(len=*),           intent(in) :: varname
    logical,                    intent(in) :: val
!OUT
    type(tke_type), intent(inout), target, optional:: tke_userdef_constants
    type(tke_type), pointer :: tke_constants_out

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tke_put_tke_logical'

  tke_constants_out=>tke_constants_saved
  if (present(tke_userdef_constants)) then
  tke_constants_out=> tke_userdef_constants
 end if

  select case(trim(varname))

    case('only_tke')
      tke_constants_out%only_tke=val
    case('use_Kappa_min')
      tke_constants_out%use_Kappa_min = val
    case('use_ubound_dirichlet')
      tke_constants_out%use_ubound_dirichlet=val
    case('use_lbound_dirichlet')
      tke_constants_out%use_lbound_dirichlet=val
    case('l_lc')
      tke_constants_out%l_lc=val
    case DEFAULT
!      print*, "ERROR:", trim(varname), " not a valid choice"
!      stop 1
      CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice')
  end select
    !!enable_GM etc can go in here
end subroutine tke_put_tke_logical

!=================================================================================

subroutine tke_put_tke_real(varname,val,tke_userdef_constants)
!This subroutine puts real values to TKE variables
!IN
    character(len=*),           intent(in) :: varname
    real(wp),                   intent(in) :: val
!OUT
    type(tke_type), intent(inout), target, optional:: tke_userdef_constants
    type(tke_type), pointer :: tke_constants_out

    CHARACTER(LEN=*), PARAMETER :: method_name = module_name//':tke_put_tke_real'

  tke_constants_out=>tke_constants_saved
  if (present(tke_userdef_constants)) then
  tke_constants_out=> tke_userdef_constants
  end if

  select case(trim(varname))

    case('c_k')
      tke_constants_out%c_k= val
    case('c_eps')
      tke_constants_out%c_eps= val
    case('cd')
      tke_constants_out%cd= val
    case('alpha_tke')
      tke_constants_out%alpha_tke = val
    case('mxl_min')
      tke_constants_out%mxl_min = val
    case('KappaM_min')
      tke_constants_out%KappaM_min = val
    case('KappaH_min')
      tke_constants_out%KappaH_min = val
    case('KappaM_max')
      tke_constants_out%KappaM_max = val
    case('tke_min')
      tke_constants_out%tke_min = val
    case('clc')
      tke_constants_out%clc = val
    case('tke_surf_min')
      tke_constants_out%tke_surf_min = val
    case DEFAULT
!      print*, "ERROR:", trim(varname), " not a valid choice"
!      stop 1
      CALL finish(method_name,'ERROR: '//TRIM(varname)//' not a valid choice')

  end select

end subroutine tke_put_tke_real

!=================================================================================

END MODULE mo_ocean_tke_base
