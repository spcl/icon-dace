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
!
! Classes and functions for the turbulent mixing package (tmx)
!


!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_tmx_smagorinsky

  USE mo_kind,                ONLY: wp, vp
  USE mo_tmx_field_class,     ONLY: t_domain
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: min_rlcell_int
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_variable_list,       ONLY: t_variable_set
  USE mo_sync,                ONLY: SYNC_C, sync_patch_array
  USE mo_physical_constants,  ONLY: grav,rgrav

#ifdef _OPENACC
  use openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Smagorinsky_init, Smagorinsky_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tmx_smagorinsky'

  ! PROCEDURE(stability_interface), POINTER :: compute_stability_term => NULL()
  ! ABSTRACT INTERFACE 
  !  FUNCTION stability_interface(mech_prod,bruvais,rturb_prandtl,jb,jc,jk) result(stability_term)
  !    IMPORT    :: wp
  !    REAL(wp), INTENT(in), POINTER :: mech_prod(:,:,:), bruvais(:,:,:)
  !    REAL(wp), INTENT(in), POINTER :: rturb_prandtl
  !    INTEGER,  INTENT(in)          :: jb,jc,jk
  !    REAL(wp) :: stability_term 
  !  END FUNCTION stability_interface
  ! END INTERFACE
  
  CONTAINS
    !============================================================================
    SUBROUTINE Smagorinsky_init(domain,config,inputs,diagnostics)

      TYPE(t_domain),        INTENT(in) :: domain
      CLASS(t_variable_set), INTENT(in), TARGET :: config
      CLASS(t_variable_set), INTENT(in), TARGET :: inputs

      CLASS(t_variable_set), INTENT(inout), TARGET :: diagnostics

      REAL(wp), POINTER, DIMENSION(:,:,:) :: mixing_length_sq, gepot_agl_ic
      REAL(vp), POINTER, DIMENSION(:,:,:) :: dzh
      REAL(wp), POINTER, DIMENSION(:,:)   :: scaling_factor_louis
      REAL(wp), POINTER :: smag_constant, max_turb_scale
      LOGICAL,  POINTER :: use_louis
     
      mixing_length_sq  => diagnostics%list%get_ptr_r3d('square of mixing length for Smagorinsky model')

      use_louis         => config%list%Get_ptr_l0d('switch to activate Louis formula')
      smag_constant     => config%list%Get_ptr_r0d('Smagorinsky constant')
      max_turb_scale    => config%list%Get_ptr_r0d('maximum turbulence length scale')
      
#ifdef __MIXED_PRECISION
      dzh               => inputs%list%get_ptr_s3d('layer thickness half') 
#else
      dzh               => inputs%list%get_ptr_r3d('layer thickness half') 
#endif

      gepot_agl_ic      => inputs%list%get_ptr_r3d('geopotential above groundlevel at interface and cell center')

      ! Compute mixing length for Smagorinsky model
      CALL compute_mixing_length(domain, dzh, gepot_agl_ic, smag_constant, max_turb_scale, mixing_length_sq)

      ! Computation of variables for stability correction function 
      IF (use_louis) THEN
        ! -> set procedure pointer for computation of stability term
        ! PROCEDURE POINTER ARE NOT SUPPORTED YET BY NVIDIA COMPILER
        !compute_stability_term => compute_stability_term_louis
       
        scaling_factor_louis  => diagnostics%list%Get_ptr_r2d('scaling factor for Louis constant b')
        __acc_attach(scaling_factor_louis)

        ! compute scaling_factor_louis in init!
        CALL compute_scaling_factor_louis(domain,scaling_factor_louis)

      !ELSE 
        !compute_stability_term => compute_stability_term_classic
      END IF

    END SUBROUTINE Smagorinsky_init
    !============================================================================
    !
    ! Computes the square of the SGS mixing length for the Smagorinsky model. 
    !
    ! lambda^2 = (Cs * Delta)^2 *(kappa*x_3)^2 / ((Cs * Delta)^2 + (kappa*x_3)^2) 
    !          = (Cs * Delta * x_3)^2 / ( (Cs * Delta / kappa)^2 + x_3^2 )
    !
    ! with   Cs    : Smagorinsky constant
    !        Delta : filter/grid width
    !        x_3   : vertical coordinate
    !        kappa : von Karman constant (kappa=0.40)
    !
    ! reference see Dipankar et al. (2015)
    !
    SUBROUTINE compute_mixing_length(domain,dzh, gepot_agl_ic,smag_constant,max_turb_scale,mixing_length_sq)

      TYPE(t_domain), INTENT(in)    :: domain
      REAL(wp),       INTENT(in)    :: smag_constant, max_turb_scale
      REAL(wp), DIMENSION(:,:,:), INTENT(in)    :: gepot_agl_ic
      REAL(vp), DIMENSION(:,:,:), INTENT(in)    :: dzh
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) :: mixing_length_sq

      REAL(wp)  :: kappa, les_filter, z_mc
      INTEGER :: jg, jk, jb, jc
      INTEGER :: nlevp1

      ! von Karman constant 
      kappa = 0.4_wp

      jg     = domain%patch%id             
      nlevp1 = domain%nlev + 1

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,z_mc,les_filter) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c, domain%i_endblk_c
        DO jk = 1 , nlevp1
          DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
            z_mc  = gepot_agl_ic(jc,jk,jb) * rgrav

            les_filter =  smag_constant           & 
                          * MIN( max_turb_scale,  &
                                (dzh(jc,jk,jb) * domain%area(jc,jb))**0.33333_wp &
                               )
            ! 
            mixing_length_sq(jc,jk,jb) = (les_filter*z_mc)**2    &
                                         / ((les_filter/kappa)**2+z_mc**2)

          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    !$ACC UPDATE DEVICE(mixing_length_sq)

    END SUBROUTINE compute_mixing_length
    !============================================================================
    !
    ! Computes the scaling factor for Louis constant b. 
    ! Scaling factor for Louis constant b is designed to be 1 with R2B8 setup.
    !
    SUBROUTINE compute_scaling_factor_louis(domain,scaling_factor_louis)   
      TYPE(t_domain),           INTENT(in)    :: domain
      REAL(wp), DIMENSION(:,:), INTENT(inout) :: scaling_factor_louis

      INTEGER                               :: jb,jc

      ! Global mean of cell area for R2B8 [m]
      REAL(wp), PARAMETER :: mean_area_R2B8 = 97294071.23714285_wp
      
!$OMP PARALLEL DO PRIVATE(jb,jc) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = domain%i_startblk_c, domain%i_endblk_c
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jc = domain%i_startidx_c(jb), domain%i_endidx_c(jb)
          scaling_factor_louis(jc,jb) = mean_area_R2B8 / domain%area(jc,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    END SUBROUTINE compute_scaling_factor_louis
    !============================================================================
    !
    ! This subroutine calls the models to compute the eddy viscosity and 
    ! diffusivity based on the Smagorinsky-Lilly eddy viscosity model. 
    ! Depending on the configuration the classical version (Lilly 1962) or the 
    ! Louis forumlation (Louis 1979) for the stability correction function is used.  
    !
    SUBROUTINE Smagorinsky_model( &
      domain,                     &
      mech_prod,                  &
      bruvais,                    &
      rho_ic,                     &
      mixing_length_sq,           &
      rturb_prandtl,              &
      use_louis,                  &
      louis_constant_b,           &
      scaling_factor_louis,       &
      patch,                      &
      km_ic,                      &
      kh_ic,                      &
      stability_function          &
      )

      TYPE(t_domain), INTENT(in)    :: domain
      REAL(wp), INTENT(in), DIMENSION(:,:,:)  :: mech_prod, bruvais, rho_ic, mixing_length_sq
      REAL(wp), INTENT(in), DIMENSION(:,:)    :: scaling_factor_louis
      REAL(wp), INTENT(in) :: rturb_prandtl, louis_constant_b
      LOGICAL,  INTENT(in) :: use_louis
      TYPE(t_patch), INTENT(in) :: patch

      REAL(wp), POINTER, INTENT(inout), DIMENSION(:,:,:) :: km_ic, kh_ic, stability_function

      IF (use_louis) THEN
        CALL Smagorinsky_model_louis(domain,mech_prod,bruvais,rho_ic,mixing_length_sq, &
                                     rturb_prandtl,louis_constant_b,scaling_factor_louis, &
                                     patch,km_ic,kh_ic,stability_function)
      ELSE
        CALL Smagorinsky_model_classic(domain,mech_prod,bruvais,rho_ic,mixing_length_sq, &
                                       rturb_prandtl,patch,km_ic,kh_ic,stability_function)
      END IF

      !$ACC WAIT

      CALL sync_patch_array(SYNC_C, patch, kh_ic)
      CALL sync_patch_array(SYNC_C, patch, km_ic)

    END SUBROUTINE Smagorinsky_model
    !============================================================================
    !
    ! This subroutine computes the eddy viscosity and diffusivity based on the 
    ! classical formulation of the stability correction term:
    !
    !    Km = rho * lambda^2 * stability_term
    !
    ! where the the stability term includes the strain rate into the stability
    ! correction function: 
    !
    !   stability term = sqrt(|S|^2 - N^2 / Pr_t )
    !   
    !   Km     : eddy viscosity
    !   |S|    : magnitude of strain rate  (-> |S|^2 = 0.5*mech_prod)
    !   N^2    : bruvais 
    !   Pr_t   : turbulent Prandtl number
    !   lambda : mixing length
    !
    SUBROUTINE Smagorinsky_model_classic( &
      domain,                     &
      mech_prod,                  &
      bruvais,                    &
      rho_ic,                     &
      mixing_length_sq,           &
      rturb_prandtl,              &
      patch,                      &
      km_ic,                      &
      kh_ic,                      &
      stability_function          &
      )

      TYPE(t_domain), INTENT(in)    :: domain
      REAL(wp), INTENT(in), DIMENSION(:,:,:) :: mech_prod, bruvais, rho_ic, mixing_length_sq
      REAL(wp), INTENT(in) :: rturb_prandtl
      TYPE(t_patch), INTENT(in) :: patch

      REAL(wp), INTENT(inout), DIMENSION(:,:,:) :: km_ic, kh_ic, stability_function

      INTEGER :: jb,jc,jk,nlev,nlevp1
      INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
      INTEGER :: rl_start, rl_end
      REAL(wp) :: stability_term     

      nlev = domain%nlev 
      nlevp1 = nlev+1

      rl_start   = 3
      rl_end     = min_rlcell_int
      i_startblk = patch%cells%start_block(rl_start)
      i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, stability_term) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                              i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
          DO jc = i_startidx, i_endidx
#endif            
            stability_term = SQRT(MAX( 0._wp, 0.5_wp * mech_prod(jc,jk,jb) - rturb_prandtl * bruvais(jc,jk,jb) ))

            km_ic(jc,jk,jb) = rho_ic(jc,jk,jb)               &
                              * mixing_length_sq(jc,jk,jb)   &
                              * stability_term

            kh_ic(jc,jk,jb) = km_ic(jc,jk,jb) * rturb_prandtl
            
          END DO
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jc = i_startidx, i_endidx
          kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
          kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
          km_ic(jc,1,jb)      = km_ic(jc,2,jb)
          km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    END SUBROUTINE Smagorinsky_model_classic
    !============================================================================
    !
    ! This subroutine computes the eddy viscosity and diffusivity based on the 
    ! Louis formulation of the stability correction term:
    !
    !    Km = rho * lambda^2 * stability_term
    !
    ! where the the stability term includes the strain rate into the stability
    ! correction function: 
    !
    !     -> stability_term = sqrt(|S|^2 * stability_factor_louis )
    !
    !     -> stability_factor_louis = 1 / (1 + b * Ri )**n
    !
    ! Km        : eddy viscosity
    ! |S|^2     : square of strain rate  (-> |S|^2 = 0.5*mech_prod)
    ! N^2       : bruvais 
    ! Pr_t      : turbulent Prandtl number
    ! lambda    : mixing length
    ! Ri        : Richardson number 
    ! b         : Louis constant 
    !
    SUBROUTINE Smagorinsky_model_louis( &
      domain,                     &
      mech_prod,                  &
      bruvais,                    &
      rho_ic,                     &
      mixing_length_sq,           &
      rturb_prandtl,              &
      louis_constant_b,           &
      scaling_factor_louis,       &
      patch,                      &
      km_ic,                      &
      kh_ic,                      &
      stability_function          &
      )
      
      TYPE(t_domain), INTENT(in)    :: domain
      REAL(wp), INTENT(in), DIMENSION(:,:,:) :: mech_prod, bruvais, rho_ic, mixing_length_sq
      REAL(wp), INTENT(in), DIMENSION(:,:)    :: scaling_factor_louis
      REAL(wp), INTENT(in) :: rturb_prandtl, louis_constant_b
      TYPE(t_patch), INTENT(in) :: patch

      REAL(wp), INTENT(inout), DIMENSION(:,:,:) :: km_ic, kh_ic, stability_function

      INTEGER :: jb,jc,jk,nlev,nlevp1
      INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
      INTEGER :: rl_start, rl_end

      REAL(wp) :: Ri, stability_term     

      nlev   = domain%nlev
      nlevp1 = domain%nlev + 1

      rl_start   = 3
      rl_end     = min_rlcell_int
      i_startblk = patch%cells%start_block(rl_start)
      i_endblk   = patch%cells%end_block(rl_end)

!$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, Ri, stability_term) ICON_OMP_DEFAULT_SCHEDULE
      DO jb = i_startblk,i_endblk
        CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
                              i_startidx, i_endidx, rl_start, rl_end)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
#ifdef __LOOP_EXCHANGE
        DO jc = i_startidx, i_endidx
          DO jk = 2 , nlev
#else
        DO jk = 2 , nlev
          DO jc = i_startidx, i_endidx
#endif
            Ri  = 2._wp * bruvais(jc,jk,jb) / mech_prod(jc,jk,jb) 

            stability_function(jc,jk,jb) =  MAX(  1.0_wp - Ri*rturb_prandtl,                &
                                                  MIN(1._wp,                                & 
                                                      1._wp/(1._wp+louis_constant_b         &
                                                            *scaling_factor_louis(jc,jb)    &
                                                            *ABS(Ri))**4                    &
                                                     ))
      
            stability_term = SQRT( 0.5_wp * mech_prod(jc,jk,jb) * stability_function(jc,jk,jb) )

            km_ic(jc,jk,jb) = rho_ic(jc,jk,jb)              &
                              * mixing_length_sq(jc,jk,jb)  &
                              * stability_term

            kh_ic(jc,jk,jb) = km_ic(jc,jk,jb) * rturb_prandtl

          END DO
        END DO
        !$ACC END PARALLEL LOOP

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO jc = i_startidx, i_endidx
          kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
          kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
          km_ic(jc,1,jb)      = km_ic(jc,2,jb)
          km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!$OMP END PARALLEL DO

    END SUBROUTINE Smagorinsky_model_louis
    !============================================================================
    !
    ! This function computes the stability correction term for the eddy viscosity. 
    ! The stability term includes the strain rate into the stability
    ! correction function, e.g. 
    !
    !   stability term = sqrt(|S|^2 - N^2 / Pr_t )
    !   
    !   in  Km = rho * lambda^2 * stability_term
    !
    ! Km     : eddy viscosity
    ! |S|    : magnitude of strain rate  (-> |S|^2 = 0.5*mech_prod)
    ! N^2    : bruvais 
    ! Pr_t   : turbulent Prandtl number
    ! lambda : mixing length
    !
    !============================================================================
    !
    ! FUNCTION compute_stability_term_classic(mech_prod,bruvais,rturb_prandtl,jb,jc,jk) result(stability_term)
    ! !$ACC ROUTINE SEQ

    !   REAL(wp), INTENT(in), POINTER :: mech_prod(:,:,:), bruvais(:,:,:)
    !   !REAL(wp), INTENT(in), POINTER :: rturb_prandtl,  mech_prod, bruvais
    !   REAL(wp), INTENT(in), POINTER :: rturb_prandtl
    !   INTEGER,  INTENT(in)          :: jb,jc,jk
    !   REAL(wp) :: stability_term 

    !   stability_term = SQRT(MAX( 0._wp, 0.5_wp * mech_prod(jc,jk,jb) - rturb_prandtl * bruvais(jc,jk,jb) ))

    ! END FUNCTION compute_stability_term_classic
    !
    !============================================================================
    !
    ! This function computes the stability term for the eddy viscosity based on 
    ! the stability correction function of Louis (1979). 
    ! The stability term includes the strain rate into the stability
    ! correction function, e.g. 
    !     Km = rho * lambda^2 * stability_term
    !
    !     -> stability_term = sqrt(|S|^2 * stability_factor_louis )
    !
    !     -> stability_factor_louis = 1 / (1 + b * Ri )**n
    !
    ! Km        : eddy viscosity
    ! |S|^2     : square of strain rate  (-> |S|^2 = 0.5*mech_prod)
    ! N^2       : bruvais 
    ! Pr_t      : turbulent Prandtl number
    ! lambda    : mixing length
    ! Ri        : Richardson number 
    ! b         : Louis constant 
    !
    ! FUNCTION compute_stability_term_louis(mech_prod,bruvais,rturb_prandtl,jb,jc,jk) result(stability_term)
    ! !$ACC ROUTINE SEQ

    !   REAL(wp), INTENT(in), POINTER :: mech_prod(:,:,:), bruvais(:,:,:)
    !   !REAL(wp), INTENT(in), POINTER :: mech_prod, bruvais
    !   REAL(wp), INTENT(in), POINTER :: rturb_prandtl
    !   !REAL(wp), INTENT(in), POINTER :: scaling_factor_louis
    !   INTEGER,  INTENT(in)          :: jb,jc,jk
    !   REAL(wp) :: stability_term, Ri, stability_factor_louis 

    !   ! Ri  = 2._wp * bruvais(jc,jk,jb)/ mech_prod(jc,jk,jb)

    !   ! stability_factor_louis = MAX(1.0_wp - Ri*rturb_prandtl,                      &
    !   !                             MIN(1._wp,                                       & 
    !   !                                 1._wp/(1._wp+louis_constant_b                &
    !   !                                        *scaling_factor_louis(jc,jb)          &
    !   !                                        *ABS(Ri))**4))
      
    !   ! stability_term = SQRT( 0.5_wp * mech_prod(jc,jk,jb) * stability_factor_louis )


    !   Ri  = 2._wp * bruvais(jc,jk,jb) / mech_prod(jc,jk,jb) 

    !   stability_function(jc,jk,jb) =  MAX(  1.0_wp - Ri*rturb_prandtl,                &
    !                                         MIN(1._wp,                                & 
    !                                             1._wp/(1._wp+louis_constant_b         &
    !                                                   *scaling_factor_louis(jc,jb)    &
    !                                                   *ABS(Ri))**4                    &
    !                                            ))
      
    !   stability_term = SQRT( 0.5_wp * mech_prod(jc,jk,jb) * stability_function(jc,jk,jb) )

    ! END FUNCTION compute_stability_term_louis
    !============================================================================
    !
    ! This subroutine computes the SGS eddy viscosity (Km) and diffusivity (Kh) 
    ! at interface cell centers.
    ! Note: at this point mech_prod is twice the actual mechanical production term.
    !--------------------------------------------------------------------------
    !   Km = rho * lambda^2 * |S| * sqrt(1 - Ri / Pr_t)
    !   Km = rho * lambda^2 * sqrt(|S|^2 - N^2 / Pr_t )
    !
    !   Kh = Km / Pr_t
    !
    !   with  lambda : mixing length 
    !         N      : Brunt-Vaisala frequency
    !         Pr_t   : turbulent Prandtl number
    !         Ri     : Richardson number (Ri = N^2/|S|^2)
    !
    !
!     SUBROUTINE Smagorinsky_model_not_working( &
!       mech_prod,                  &
!       bruvais,                    &
!       rho_ic,                     &
!       rturb_prandtl,              &
!       patch,                      &
!       p_nh_metrics,               &
!       km_ic,                      &
!       kh_ic                       &
!       )

!       REAL(wp), INTENT(in), POINTER :: mech_prod(:,:,:), bruvais(:,:,:), rho_ic(:,:,:)
!       REAL(wp), INTENT(in), POINTER :: km_ic(:,:,:), kh_ic(:,:,:)
!       REAL(wp), INTENT(in), POINTER :: rturb_prandtl
!       TYPE(t_patch), INTENT(in), POINTER :: patch

!       INTEGER :: jb,jc,jk
!       INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
!       INTEGER :: rl_start, rl_end

!       REAL(wp) :: stability_term     

!       rl_start   = 3
!       rl_end     = min_rlcell_int
!       i_startblk = patch%cells%start_block(rl_start)
!       i_endblk   = patch%cells%end_block(rl_end)

! !$OMP PARALLEL DO PRIVATE(jb, jk, jc, i_startidx, i_endidx, stability_term) ICON_OMP_DEFAULT_SCHEDULE
!       DO jb = i_startblk,i_endblk
!         CALL get_indices_c(patch, jb, i_startblk, i_endblk, &
!                               i_startidx, i_endidx, rl_start, rl_end)

!       !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
!       !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
! #ifdef __LOOP_EXCHANGE
!         DO jc = i_startidx, i_endidx
!           DO jk = 2 , nlev
! #else
!         DO jk = 2 , nlev
!           DO jc = i_startidx, i_endidx
! #endif

!             stability_term = compute_stability_term( mech_prod,bruvais,       &
!                                                      rturb_prandtl,jb,jc,jk)

!             km_ic(jc,jk,jb) = rho_ic(jc,jk,jb) *                          &
!                               mixing_length_sq(jc,jk,jb) *   &
!                               stability_term

!             kh_ic(jc,jk,jb) = km_ic(jc,jk,jb) * rturb_prandtl
!           END DO
!         END DO
!         !$ACC END PARALLEL LOOP

!         !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
!         !$ACC LOOP GANG VECTOR
!         DO jc = i_startidx, i_endidx
!           kh_ic(jc,1,jb)      = kh_ic(jc,2,jb)
!           kh_ic(jc,nlevp1,jb) = kh_ic(jc,nlev,jb)
!           km_ic(jc,1,jb)      = km_ic(jc,2,jb)
!           km_ic(jc,nlevp1,jb) = km_ic(jc,nlev,jb)
!         END DO
!         !$ACC END PARALLEL
!       END DO
! !$OMP END PARALLEL DO

!       !$ACC WAIT

!       CALL sync_patch_array(SYNC_C, patch, kh_ic)
!       CALL sync_patch_array(SYNC_C, patch, km_ic)

!     END SUBROUTINE Smagorinsky_model_not_working
!     !============================================================================
END MODULE mo_tmx_smagorinsky
