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

! Module containing subroutines for the initialization of a bubble on a torus
!
! Literature  cr2021_08_03_jsr for further documentation

MODULE mo_aes_bubble
  USE mo_kind,                ONLY: wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_physical_constants,  ONLY: rd, grav, p0ref,rd_o_cpd, o_m_rdv
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics, t_nh_ref
  USE mo_parallel_config,     ONLY: nproma
  USE mo_aes_bubble_config,   ONLY: aes_bubble_config
  USE mo_run_config,          ONLY: iqv
  USE mo_aes_thermo,          ONLY: specific_humidity, sat_pres_water
  USE mo_hydro_adjust,        ONLY: hydro_adjust_iterative

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_aes_bubble

  CONTAINS
  
  SUBROUTINE init_aes_bubble( ptr_patch, ptr_nh_prog, ptr_nh_ref, ptr_nh_diag, ptr_metrics)

    ! INPUT PARAMETERS:
    TYPE(t_patch),TARGET,  INTENT(IN)   :: &  !< patch on which computation is performed
      &  ptr_patch
    TYPE(t_nh_prog),       INTENT(INOUT):: &  !< prognostic state vector
      &  ptr_nh_prog
    TYPE(t_nh_diag),       INTENT(INOUT):: &  !< diagnostic state vector
      &  ptr_nh_diag
    TYPE(t_nh_metrics),    INTENT(IN)   :: &  !< NH metrics state
      &  ptr_metrics
    TYPE(t_nh_ref),        INTENT(INOUT):: &  !< reference state vector
      &  ptr_nh_ref

    INTEGER  :: jb,jk,jl,nblks_c,npromz_c,nlen,nlev  !<loop indices and control

    REAL(wp), ALLOCATABLE :: temp(:,:,:), rh(:,:,:), zeta_xy(:,:)
    REAL(wp)              :: pres(nproma), tt(nproma)
    REAL(wp)              :: qv, x, y, z, dz, sat_pres, wat_pres, tv, sigma_x, sigma_z, zeta

    REAL(wp), POINTER     :: psfc,        t_am,      t0,        gamma0,      &
                           & gamma1,      z0,        t_perturb, relhum_bg,   &
                           & relhum_mx,   hw_x,      hw_z,      x_center
    LOGICAL, POINTER      :: lgaussxy

    psfc        => aes_bubble_config%psfc
    t_am        => aes_bubble_config%t_am
    t0          => aes_bubble_config%t0
    gamma0      => aes_bubble_config%gamma0
    gamma1      => aes_bubble_config%gamma1
    z0          => aes_bubble_config%z0
    t_perturb   => aes_bubble_config%t_perturb
    relhum_bg   => aes_bubble_config%relhum_bg
    relhum_mx   => aes_bubble_config%relhum_mx
    hw_x        => aes_bubble_config%hw_x
    hw_z        => aes_bubble_config%hw_z
    x_center    => aes_bubble_config%x_center
    lgaussxy    => aes_bubble_config%lgaussxy
    
    ! values for the blocking
    nblks_c  = ptr_patch%nblks_c
    npromz_c = ptr_patch%npromz_c

    ! number of vertical levels
    nlev   = ptr_patch%nlev
    
    ALLOCATE(temp(nproma,nlev,nblks_c))
    ALLOCATE(rh(nproma,nlev,nblks_c))
    ALLOCATE(zeta_xy(nproma,nblks_c))
    
    ! init surface pressure
    ptr_nh_diag%pres_sfc(:,:) = psfc
  
    ! Tracers: all zero by default
    ptr_nh_prog%tracer(:,:,:,:) = 0._wp

    ! Parameters for horizontal Gaussian profile
    sigma_x = hw_x / (2._wp * SQRT(2._wp * LOG(2._wp)))
    sigma_z = hw_z / (2._wp * SQRT(2._wp * LOG(2._wp)))
    DO jb = 1, nblks_c
      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      END IF
      ! Calculate gaussians in horizontal direction first
      ! Either in x direction only (lgaussxy=.FALSE.) or
      ! in x and y direction with the same sigma_x and x_center
      IF (lgaussxy) THEN
         DO jl = 1, nlen
            x = ptr_patch%cells%cartesian_center(jl, jb)%x(1)
            y = ptr_patch%cells%cartesian_center(jl, jb)%x(2)
            zeta_xy(jl,jb) = gaussian(x_center, sigma_x, 1._wp, x) * &
                           & gaussian(x_center, sigma_x, 1._wp, y)
         END DO
      ELSE
         DO jl = 1, nlen
            x = ptr_patch%cells%cartesian_center(jl, jb)%x(1)
            zeta_xy(jl,jb) = gaussian(x_center, sigma_x, 1._wp, x)
         END DO
      END IF
    END DO
    
    DO jb = 1, nblks_c

      IF (jb /= nblks_c) THEN
         nlen = nproma
      ELSE
         nlen = npromz_c
      END IF

      tt(1:nlen)   = t0 + t_perturb * zeta_xy(1:nlen,jb)
      pres(1:nlen) = psfc
      
      DO jk = nlev, 1, -1
        ! Calculation of temperature profile, pressure, Exner pressure, and potential virtual temperature 
        ! We perform a simple integration of the hydrostatic equation from bottom to top of atm.
          DO jl = 1, nlen
            z = ptr_metrics%z_mc(jl, jk, jb)
            IF ( jk == nlev ) THEN
              dz = z
            ELSE
              dz = z-ptr_metrics%z_mc(jl, jk+1, jb)
            END IF
            ! Set a 2D patch of moisture at height 'bub_ver_width' [m] above the
            ! centre of the experiment domain. The max moisture is set to
            ! saturation pressure and decreases exponentially towards 0.
            ! Calculate the temperature profile, (instable in the PBL)
            zeta=zeta_xy(jl,jb)*gaussian(0._wp, sigma_z, 1._wp, z)
            IF ( z <= z0 ) THEN ! we are within the PBL
              temp(jl, jk, jb) = t0 - gamma0 * z
            ELSE
              temp(jl, jk, jb) = t0 - gamma0 * z0 - gamma1 * (z-z0)
            END IF
            temp(jl, jk, jb) = MAX(temp(jl, jk, jb), t_am) + t_perturb * zeta
            sat_pres = sat_pres_water(temp(jl, jk, jb))
            wat_pres = sat_pres * ((relhum_mx-relhum_bg) * zeta + relhum_bg) 
              ! The gaussian is defined such that the third argument being =1 means gaussian(0,sigma_x,1,x)=1
            tt(jl)   = 0.5_wp * (tt(jl) + temp(jl, jk, jb))
            pres(jl) = pres(jl) * EXP(-grav/rd/tt(jl) * dz)
            tt(jl)   = temp(jl,jk,jb)
            tv                                  = tt(jl) / (1._wp - (wat_pres / pres(jl)) * o_m_rdv)
            qv                                  = specific_humidity(wat_pres, pres(jl))
            rh(jl, jk, jb)                      = wat_pres/sat_pres
            ptr_nh_prog%exner(jl, jk, jb)       = (pres(jl) / p0ref)**rd_o_cpd
            ptr_nh_diag%pres(jl, jk, jb)        = pres(jl)
            ptr_nh_prog%rho(jl, jk, jb)         = pres(jl) / rd / tv
            ptr_nh_prog%theta_v(jl, jk, jb)     = tv / ptr_nh_prog%exner(jl, jk, jb)
            ptr_nh_prog%tracer(jl, jk, jb, iqv) = qv
          END DO !jl

      END DO !jk
    END DO ! jb       

    ! Adjust preliminary profiles to numerics of ICON dynamical core
    ! Relative humidity is constant throughout domain and levels, but needs to be stored in a 3d var.
    CALL hydro_adjust_iterative(                                                                  &
         & ptr_patch,                             ptr_metrics,                                    &
         & temp,                                  rh, ptr_nh_prog%exner,                          &
         & ptr_nh_prog%theta_v,                   ptr_nh_prog%rho,                                &
         & ptr_nh_prog%tracer(:,:,:,iqv),         luse_exner_fg=.TRUE.,                           &
         & opt_exner_lbc=ptr_nh_prog%exner(:,nlev,:)                                              )

    
  !meridional and zonal wind
  ptr_nh_prog%vn = 0._wp
  ptr_nh_ref%vn_ref = ptr_nh_prog%vn

  !vertical wind
  ptr_nh_prog%w = 0._wp
  ptr_nh_ref%w_ref = ptr_nh_prog%w

  END SUBROUTINE init_aes_bubble

!!!=============================================================================================

  FUNCTION gaussian(mu, sigma, delta, x)
    ! calculate value of gaussian at x, for gaussian centered at mu with
    ! standard deviation sigma, and amplitude delta
    !
    REAL(wp), INTENT(in) :: mu    !< expectation of the Gaussian
    REAL(wp), INTENT(in) :: sigma !< std-dev. of Gaussian
    REAL(wp), INTENT(in) :: delta !< max. diviation from 0
    REAL(wp), INTENT(in) :: x     !< x value
    REAL(wp) :: gaussian, xx

    xx = ((x - mu) / sigma)
    gaussian = EXP(-.5_wp * xx * xx) * delta
  END FUNCTION gaussian
  
END MODULE mo_aes_bubble
