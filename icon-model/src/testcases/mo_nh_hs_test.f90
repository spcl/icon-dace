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

! Held-Suarez test for the NH-Core
!
! Literature
! - Held, I. M. and Suarez, M. J. (1994): A Proposal for the Intercomparison
!   of the Dynamical Cores of Atmospheric General Circulation Models

MODULE mo_nh_hs_test

!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2008
!
!-------------------------------------------------------------------------
!
!
!

  USE mo_kind,                ONLY: wp
  USE mo_physical_constants,  ONLY: rd, rd_o_cpd, p0ref, grav, cpd, cvd
  USE mo_time_base,           ONLY: rdaylen
  USE mo_model_domain,        ONLY: t_patch
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_parallel_config,     ONLY: nproma


  IMPLICIT NONE

  PRIVATE

  REAL(wp), PARAMETER :: zp_hs      = 100000._wp            !< surface pressure
  REAL(wp), PARAMETER :: zt_hs      = 300._wp               !< atmospheric temperature
!REAL(wp), PARAMETER :: zt_hs      = 250._wp               !< atmospheric temperature

  ! Module parameters

  REAL(wp), PARAMETER :: HSsigb    = 0.7_wp
  REAL(wp), PARAMETER :: HSvcoeff1 = 1._wp/(1._wp-HSsigb)
  REAL(wp), PARAMETER :: HSvcoeff2 = -HSsigb/(1._wp-HSsigb)

  REAL(wp), PARAMETER :: HSkf    = 1._wp/ rdaylen
  REAL(wp), PARAMETER :: HSka    = 1._wp/(rdaylen*40._wp)
  REAL(wp), PARAMETER :: HSks    = 1._wp/(rdaylen* 4._wp )

  REAL(wp), PARAMETER :: HSdty   = 60._wp
  REAL(wp), PARAMETER :: HSdthz  = 10._wp

  REAL(wp), PARAMETER :: HScappa = rd/cpd
  REAL(wp), PARAMETER :: HSp0    = 1.E5_wp
  REAL(wp), PARAMETER :: HSt0    = 200._wp
  REAL(wp), PARAMETER :: HSt1    = 315._wp

  PUBLIC :: init_nh_state_prog_held_suarez
  PUBLIC :: held_suarez_forcing_vn
  PUBLIC :: held_suarez_forcing_temp

!--------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!
  !>
  !!               Initialization of prognostic state vector.
  !!
  !!
  SUBROUTINE init_nh_state_prog_held_suarez( ptr_patch, ptr_nh_prog, ptr_nh_diag, &
    &                                   ptr_ext_data, p_metrics )

    TYPE(t_patch), TARGET,INTENT(INOUT) :: &  !< patch on which computation is performed
      &  ptr_patch

    TYPE(t_nh_prog), INTENT(INOUT)      :: &  !< prognostic state vector
      &  ptr_nh_prog

    TYPE(t_nh_diag), INTENT(INOUT)      :: &  !< diagnostic state vector
      &  ptr_nh_diag

    TYPE(t_external_data), INTENT(INOUT):: &  !< external data
      &  ptr_ext_data

    TYPE(t_nh_metrics), INTENT(IN)      :: p_metrics !< NH metrics state


    INTEGER               :: nblks_c, npromz_c, nlen
    INTEGER               :: nlev                 !< number of full levels
    INTEGER               :: jk, jb, jc  ! loop variables

    REAL(wp)              :: zscale_h             !< initialized variables
    REAL(wp), ALLOCATABLE :: z_sfc(:,:,:)



!--------------------------------------------------------------------
!
    zscale_h = rd*zt_hs/grav

    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    nlev = ptr_patch%nlev

    ALLOCATE (z_sfc(nproma, 1, ptr_patch%nblks_c))

    ! topography
    ptr_ext_data%atm%topography_c(:,:) = 0.0_wp

    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = zp_hs

    ! init temperature
    ptr_nh_diag%temp(:,:,:)   = zt_hs


!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,nlen)
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      ENDIF

      DO jk = 1, nlev
        DO jc = 1, nlen

          ! compute full level pressure
          z_sfc(jc,1,jb) = ptr_ext_data%atm%topography_c(jc,jb)
          ptr_nh_diag%pres(jc,jk,jb) = zp_hs                                      &
            &                  * exp(-(p_metrics%z_mc(jc,jk,jb)-z_sfc(jc,1,jb))/zscale_h)

          ! init virtual potential temperature
          ptr_nh_prog%theta_v(jc,jk,jb) = zt_hs                                   &
            & * (p0ref/ptr_nh_diag%pres(jc,jk,jb))**rd_o_cpd

          ! init density field rho
          ptr_nh_prog%rho(jc,jk,jb) = ptr_nh_diag%pres(jc,jk,jb)            &
            &                           / (rd * zt_hs)

          ! init exner pressure
          ptr_nh_prog%exner(jc,jk,jb) = (ptr_nh_diag%pres(jc,jk,jb)/p0ref)**rd_o_cpd

        ENDDO !jc
      ENDDO !jk
    ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !
    ! initialize horizontal velocity field
    !

    ptr_nh_prog%vn(:,:,:) =     0.0_wp

    !
    ! initialize vertical velocity field
    !

    ptr_nh_prog%w(:,:,:) =      0.0_wp


  END SUBROUTINE init_nh_state_prog_held_suarez


  !>
  !! Linear damping of velocity in the 'boundary layer'.
  !!
  SUBROUTINE held_suarez_forcing_vn( pvn, psigma,           & !in
                                   & nlev, nproma, is, ie,  & !in
                                   & fvn_hs )                 !out

    INTEGER, INTENT(IN) :: nlev            !< number of vertical layers
    INTEGER, INTENT(IN) :: nproma, is, ie  !<

    REAL(wp),INTENT(IN) :: pvn      ( nproma, nlev ) !< normal velocity
    REAL(wp),INTENT(IN) :: psigma   ( nproma, nlev ) !< sigma = pres/pres_sfc

    REAL(wp),INTENT(INOUT) :: fvn_hs( nproma, nlev ) !< forcing on velocity

    REAL(wp) :: ztmp
    INTEGER  :: i, jk
    !---

    !$ACC DATA PRESENT(pvn, psigma, fvn_hs)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk=1,nlev
      DO i = is, ie
        ztmp = psigma(i,jk)*HSvcoeff1 + HSvcoeff2
        fvn_hs(i,jk) = -HSkf*MAX( 0._wp,ztmp )*pvn(i,jk)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC END DATA
   

  END SUBROUTINE held_suarez_forcing_vn
  !-------------
  !>
  !! Newtonian cooling.
  !!
  SUBROUTINE held_suarez_forcing_temp( ptemp_mc, ppres_mc,   &! in
                                     & psigma, plat,         &! in
                                     & nlev, nproma, is, ie, &! in
                                     & fT_hs,                &! out
                                     & opt_ekinh, opt_ldissip_heat)

    INTEGER, INTENT(IN) :: nlev                   !< number of vertical layers
    INTEGER, INTENT(IN) :: nproma, is, ie
    REAL(wp),INTENT(IN) :: ptemp_mc (nproma,nlev) !< temperature in Kelvin
    REAL(wp),INTENT(IN) :: ppres_mc (nproma,nlev) !< pressure in Pa
    REAL(wp),INTENT(IN) :: psigma   (nproma,nlev) !< sigma = pres/pres_sfc
    REAL(wp),INTENT(IN) :: plat     (nproma)      !< latitide in radians

    REAL(wp),INTENT(IN), OPTIONAL :: opt_ekinh(nproma,nlev) !< kinetic energy
    LOGICAL, OPTIONAL  :: opt_ldissip_heat        !< dissipative heating or not

    REAL(wp),INTENT(OUT):: fT_hs (nproma,nlev) !< forcing on temperature

    INTEGER :: jk  !vertical layer index
    INTEGER :: i   ! nproma index
    REAL(wp) :: zsinlat2(nproma), zcoslat2(nproma), zcoslat4(nproma)
    REAL(wp) :: zsigma0, zTempEq, ztmp, kT_hs
    LOGICAL  :: l_friheat

    !------------------------------

    ! check, whether dissipative heating should be performed
    IF ( PRESENT(opt_ldissip_heat) ) THEN
      l_friheat=opt_ldissip_heat
    ELSE
      l_friheat=.FALSE.
    ENDIF

    ! latitude related parameters

    !$ACC DATA CREATE(zsinlat2, zcoslat2, zcoslat4) PRESENT(ptemp_mc, ppres_mc, psigma, plat, fT_hs)
    !$ACC DATA PRESENT(opt_ekinh) IF(PRESENT( opt_ekinh ))
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR
    DO i=is, ie
      zsinlat2(i) = SIN(plat(i))**2
      zcoslat2(i) = 1._wp - zsinlat2(i)
      zcoslat4(i) = zcoslat2(i)**2
    ENDDO
    !$ACC END PARALLEL
   
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(kT_hs, zsigma0, zTempEq, ztmp)
    DO jk = 1,nlev

      DO i=is,ie
        ! equilibrium temperature

        zsigma0 = ppres_mc(i,jk)/HSp0
        zTempEq = HSt1 - HSdty*zsinlat2(i)               &
                       & - HSdthz*LOG(zsigma0)*zcoslat2(i)
        zTempEq = zTempEq * zsigma0**HScappa
        zTempEq = MAX( 200._wp, zTempEq )

        ! pressure-dependent coefficient

        ztmp   = psigma(i,jk)*HSvcoeff1 + HSvcoeff2
        kT_hs  = HSka + (HSks-HSka)*zcoslat4(i)*MAX(0._wp,ztmp)

        ! Newtonian cooling

        fT_hs(i,jk) = -kT_hs *( ptemp_mc(i,jk)-zTempEq )
      ENDDO

    ENDDO !vertical layer loop
    !$ACC END PARALLEL

    IF (l_friheat) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(ztmp)
      DO jk = 1,nlev
        DO i = is,ie
          ztmp = psigma(i,jk)*HSvcoeff1 + HSvcoeff2
          fT_hs(i,jk) = fT_hs(i,jk) + HSkf*MAX( 0._wp,ztmp)*2.0_wp*opt_ekinh(i,jk)/cvd
        ENDDO
      ENDDO
      !$ACC END PARALLEL
    ENDIF
    !$ACC WAIT(1)

    !$ACC END DATA
    !$ACC END DATA

  END SUBROUTINE held_suarez_forcing_temp

!--------------------------------------------------------------------
  END MODULE mo_nh_hs_test
