!
! Module to compute effective radius consistent with microphysics, cloud scheme
! and convection scheme choice (not yet!).
! The effective radius calculated here can be used by the radiation module ECRAD,
! and by some satellite forward operators (like VISOP)
!
! The idea is to keep a consistent effective radius for the whole model and forward
! operators. For this reason the coefficients and concentration can be claculated
! in the microphysics (the module only provides an interface).
!
! Description:
! The module also contains adapted versions of the  routines developed
! by Simon Gruber and Uli Blahak for the optical properties in RRTM
! (only the effective radius, not the optical porperties),
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

MODULE mo_reff_main

  USE mo_kind              ,   ONLY: wp, i4
  USE mo_math_constants    ,   ONLY: pi
  USE mo_physical_constants,   ONLY: rhoh2o, t0 => tmelt, rhoi
  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_reff_types,           ONLY: t_reff_calc
  USE mo_2mom_mcrph_driver,    ONLY: two_mom_reff_coefficients 
  USE mo_parallel_config,      ONLY: nproma
  USE mo_radiation_config,     ONLY: irad_aero, iRadAeroTegen, iRadAeroCAMSclim, iRadAeroCAMStd
  USE mo_cpl_aerosol_microphys,ONLY: ncn_from_tau_aerosol_speccnconst_dust, ice_nucleation
  USE microphysics_1mom_schemes, ONLY: get_params_for_ncn_calculation, get_params_for_reff_coefficients, get_cloud_number
  USE mo_index_list,           ONLY: generate_index_list_batched


  IMPLICIT NONE
  PRIVATE
  
  PUBLIC:: init_reff_calc, mapping_indices, calculate_ncn, calculate_reff, set_max_reff, combine_reff


  CONTAINS

! ------------------------------------------------------------------------------------------

  ! Subroutine that provides coefficients for the effective radius calculations
  ! consistent with microphysics
  SUBROUTINE one_mom_reff_coefficients( reff_calc ,return_fct)
    TYPE(t_reff_calc), INTENT(INOUT)  :: reff_calc           ! Structure with options and coefficiencts
    LOGICAL,           INTENT(INOUT)  :: return_fct          ! Return code of the subroutine

    ! Parameters used in the paramaterization of reff (the same for all)
    REAL(wp)                          :: a_geo, b_geo        ! Geometrical factors x =a_geo D**[b_geo]
    REAL(wp)                          :: mu,nu,N0            ! Parameters if the gamma distribution      
    REAL(wp)                          :: bf, bf2             ! Broadening factors of reff
    REAL(wp)                          :: x_min,x_max         ! Maximum and minimum mass of particles (kg)
    LOGICAL                           :: monodisperse        ! .true. for monodisperse DSD assumption

    REAL(wp) :: zami, zmi0, zmimax, zn0r, mu_rain, ageo_snow, zbms, zmsmin

    CHARACTER(len=*), PARAMETER :: routine = 'one_mom_reff_coefficients'
    
    N0 = -1.0_wp

    ! Check input return_fct
    IF (.NOT. return_fct) THEN
      WRITE (message_text,*) 'Reff: Function one_mom_provide_reff_coefficients (1mom) entered with previous error'
      CALL message('',message_text)
      RETURN
    END IF

    CALL get_params_for_reff_coefficients(zami_arg=zami, &
                                          zmi0_arg=zmi0, &
                                          zmimax_arg=zmimax, &
                                          zn0r_arg=zn0r, &
                                          mu_rain_arg=mu_rain, &
                                          ageo_snow_arg=ageo_snow, &
                                          zbms_arg=zbms, &
                                          zmsmin_arg=zmsmin)

    ! Default values
    x_min           = reff_calc%x_min
    x_max           = reff_calc%x_max

    ! Extract parameterization parameters
    SELECT CASE ( reff_calc%hydrometeor ) ! Select Hydrometeor
    CASE (0)                              ! Cloud water
      a_geo         = pi/6.0_wp * rhoh2o
      b_geo         = 3.0_wp
      monodisperse  = .true.              ! Monodisperse assumption by default
    CASE (1)                              ! Ice
      a_geo         = zami
      b_geo         = 3.0_wp              ! According to COSMO Documentation
      monodisperse  = .true.              ! Monodisperse assumption by default
      x_min         = zmi0                ! Limits to crystal mass set by the scheme
      x_max         = zmimax 
    CASE (2)                              ! Rain
      a_geo         = pi/6.0_wp * rhoh2o  ! Assume spherical rain
      b_geo         = 3.0_wp
      monodisperse  = .false.  
      N0            = zn0r 
      nu            = mu_rain             ! This is right, there are different mu/nu notations
      mu            = 1.0   
    CASE (3)                              ! Snow
      a_geo         = ageo_snow 
      b_geo         = zbms
      monodisperse  = .false.  
      N0            = 1.0_wp              ! Complex dependency for N0 (set in calculate_ncn)
      nu            = 0.0_wp              ! Marshall Palmer distribution (exponential)
      mu            = 1.0_wp
      x_min         = zmsmin  
    CASE (4)                              ! Graupel: values from Documentation (not in micro. code)
      a_geo         = 169.6_wp   
      b_geo         = 3.1_wp
      monodisperse  = .false.  
      N0            = 4.0E6_wp 
      nu            = 0.0_wp              ! Marshall Palmer distribution (exponential)
      mu            = 1.0_wp
    CASE DEFAULT
      CALL finish(TRIM(routine),'wrong value for reff_calc%hydrometeor')      
    END SELECT

    ! Set values if changed
    reff_calc%x_min = x_min
    reff_calc%x_max = x_max

    ! Overwrite dsd options if they are provided
    IF ( reff_calc%dsd_type /= 0) THEN
      ! Overwrite monodisperse/polydisperse according to options
      SELECT CASE (reff_calc%dsd_type)
      CASE (1)
        monodisperse  = .true.
      CASE (2)
        monodisperse  = .false.
      CASE DEFAULT
        CALL finish(TRIM(routine),'wrong value for reff_calc%dsd_type')      
      END SELECT
     
      IF ( reff_calc%dsd_type == 2) THEN    ! Overwrite mu and nu coefficients
        mu            = reff_calc%mu
        nu            = reff_calc%nu
      END IF
    END IF

    ! Calculate parameters to calculate effective radius
    SELECT CASE ( reff_calc%reff_param )  ! Select Parameterization
    CASE(0)                               ! Spheroids  Dge = c1 * x**[c2], which x = mean mass
      ! First calculate monodisperse
      reff_calc%reff_coeff(1) = a_geo**(-1.0_wp/b_geo)
      reff_calc%reff_coeff(2) = 1.0_wp/b_geo

      ! Broadening for not monodisperse
      IF ( .NOT. monodisperse ) THEN 
        bf =  GAMMA( (nu + 4.0_wp)/ mu) / GAMMA( (nu + 3.0_wp)/ mu) * &
          & ( GAMMA( (nu + 1.0_wp)/ mu) / GAMMA( (b_geo + nu + 1.0_wp)/ mu) )**(1.0_wp/b_geo)
        reff_calc%reff_coeff(1) = reff_calc%reff_coeff(1)*bf
      END IF       
     
    CASE (1) ! Fu Random Hexagonal needles:  Dge = 1/(c1 * x**[c2] + c3 * x**[c4])
             ! Parameterization based on Fu, 1996; Fu et al., 1998; Fu ,2007

      ! First calculate monodisperse
      reff_calc%reff_coeff(1) = SQRT( 3.0_wp *SQRT(3.0_wp) * rhoi / 8.0_wp ) *a_geo**(-1.0_wp/2.0_wp/b_geo)
      reff_calc%reff_coeff(2) = (1.0_wp-b_geo)/2.0_wp/b_geo
      reff_calc%reff_coeff(3) = SQRT(3.0_wp)/4.0_wp*a_geo**(1.0_wp/b_geo)
      reff_calc%reff_coeff(4) = -1.0_wp/b_geo

      ! Broadening for not monodisperse. Generalized gamma distribution
      IF ( .NOT. monodisperse ) THEN 
        bf  =  GAMMA( ( b_geo + 2.0_wp * nu + 3.0_wp)/ mu/2.0_wp ) / GAMMA( (b_geo + nu + 1.0_wp)/ mu) * &
           & ( GAMMA( (nu + 1.0_wp)/ mu) / GAMMA( (b_geo + nu + 1.0_wp)/ mu) )**( (1.0_wp-b_geo)/2.0_wp/b_geo)

        bf2 =  GAMMA( (b_geo + nu )/ mu ) / GAMMA( (b_geo + nu + 1.0_wp)/ mu) * &
           & ( GAMMA( (nu + 1.0_wp)/ mu ) / GAMMA( (b_geo + nu + 1.0_wp)/ mu) )**( -1.0_wp/b_geo)

        reff_calc%reff_coeff(1) = reff_calc%reff_coeff(1)*bf
        reff_calc%reff_coeff(3) = reff_calc%reff_coeff(3)*bf2
      END IF

    CASE DEFAULT
      CALL finish(TRIM(routine),'wrong value for reff_calc%reff_param')      
    END SELECT

    ! Calculate coefficients to calculate n from qn, in case N0 is available
    IF ( N0 > 0.0_wp) THEN 
      reff_calc%ncn_coeff(1) =  GAMMA ( ( nu + 1.0_wp)/mu) / a_geo / GAMMA (( b_geo + nu + 1.0_wp)/mu) * & 
           &( a_geo * N0 / mu * GAMMA ( (b_geo + nu + 1.0_wp)/mu) ) ** ( b_geo / (nu + b_geo +1.0_wp)) 
      reff_calc%ncn_coeff(2) = (nu +1.0_wp)/(nu + b_geo + 1.0_wp)  
      ! Scaling coefficent of N0 (only for snow)
      reff_calc%ncn_coeff(3) = b_geo/(nu + b_geo + 1.0_wp)         
    END IF
    
  END SUBROUTINE one_mom_reff_coefficients

  ! This function provides the number concentration of hydrometoers consistent with the
  ! one moment scheme. 
  ! It contains copied code from the 1 moment scheme, because many functions are hard-coded.
  SUBROUTINE one_mom_calculate_ncn( ncn, return_fct, reff_calc, k_start,             &
       &                           k_end, indices, n_ind,                            &
       &                           icpl_aero_ice, cams5, cams6, z_ifc, aer_dust,     &
       &                           q, t, rho, surf_cloud_num)

    REAL(wp)         , INTENT(INOUT)     ::  ncn(:,:)           ! Number concentration
    LOGICAL          , INTENT(INOUT)     ::  return_fct         ! Return code of the subroutine
    TYPE(t_reff_calc), INTENT(IN)        ::  reff_calc          ! Structure with options and coefficiencts
    INTEGER          , INTENT(IN)        ::  k_start, k_end     ! Start, end total indices    
    INTEGER (KIND=i4), INTENT(IN)        ::  indices(:,:)       ! Mapping for going through array
    INTEGER (KIND=i4), INTENT(IN)        ::  n_ind(:)

    INTEGER, INTENT(IN)                          :: icpl_aero_ice   ! aerosols ice nucleation scheme
    REAL(wp),INTENT(IN), POINTER, DIMENSION(:,:) :: cams5, cams6    ! CAMS dust mixing ratios
    REAL(wp),INTENT(IN), DIMENSION(:,:)          :: z_ifc           ! height at interface levels
    REAL(wp),INTENT(IN), POINTER, DIMENSION(:)   :: aer_dust        ! Tegen dust total column mass
    REAL(wp),INTENT(IN), POINTER, DIMENSION(:,:) :: q               ! Mixing ratio of hydrometeor
    
    REAL(wp), OPTIONAL, INTENT(IN)            ::  t(:,:)             ! Temperature
    REAL(wp), OPTIONAL, INTENT(IN)            ::  rho(:,:)           ! Mass density of air


    REAL(wp), OPTIONAL, INTENT(IN)       ::  surf_cloud_num(:)  ! Number concentration at surface
                                                                !CALL WITH prm_diag%cloud_num(is:ie,:) 
    ! --- End of input/output variables.

    INTEGER                              ::  jc, k, ic          ! Running indices
    ! Indices array vectorization 
    LOGICAL                              ::  well_posed         ! Logical that indicates if enough data for calculations

    ! Variables for Ice parameterization 
    REAL(wp)                             ::  znimax, znimix     ! Maximum and minimum of ice concentration
    REAL(wp)                             ::  aerncn             ! CAMS dust aerosols number concentration

    ! This is constant in both cloudice and graupel
    LOGICAL                              ::  lsuper_coolw = .true.   

    ! Variables for snow parameterization
    REAL(wp)                             ::  ztc, zn0s, nnr, hlp, alf, bet, m2s, m3s, zlog_10

    REAL(wp)            :: ageo_snow, zn0s1, zn0s2, &
                             znimax_Thom, mma(10),  &
                             mmb(10), dummy
    INTEGER             :: isnow_n0temp
    REAL(wp)            :: cloud_num

    zlog_10 = LOG(10._wp) ! logarithm of 10

    ! Check input return_fct

    IF (.NOT. return_fct) THEN
      WRITE (message_text,*) 'Reff: Function one_mom_calculate_ncn in mo_reff_main entered with previous error'
      CALL message('',message_text)
      RETURN
    END IF

    CALL get_params_for_ncn_calculation(isnow_n0temp_arg=isnow_n0temp, &
                                        ageo_snow_arg=ageo_snow,&
                                        zn0s1_arg=zn0s1, &
                                        zn0s2_arg=zn0s2, &
                                        znimax_Thom_arg=znimax_Thom, &
                                        mma_arg=mma, &
                                        mmb_arg=mmb)
      
    SELECT CASE ( reff_calc%hydrometeor )   ! Select Hydrometeor
    CASE (0)   ! Cloud water from surface field cloud_num field or fixed

#ifdef _OPENACC
      CALL finish('one_moment_calculate_ncn:','CASE hydrometeor=0 not available on GPU')
#endif

      IF (PRESENT(surf_cloud_num)) THEN
        
        DO k = k_start,k_end          
          DO ic  = 1,n_ind(k)
            jc =  indices(ic,k)
            ncn(jc,k) = surf_cloud_num(jc) ! Notice no vertical dependency
          END DO
        END DO
        
      ELSE
        
        CALL get_cloud_number(cloud_num)
        DO k = k_start,k_end
          DO ic  = 1,n_ind(k)
            jc = indices(ic,k)
            ncn(jc,k) = cloud_num ! Set constant value
          END DO
        END DO
      
      ENDIF

    CASE (1)   ! Ice
      IF (icpl_aero_ice == 1) THEN
        SELECT CASE(irad_aero)
          CASE (iRadAeroCAMSclim, iRadAeroCAMStd)! use DeMott with CAMS dust aerosols
            !$ACC DATA PRESENT(n_ind, indices, ncn, rho, t, cams5, cams6)
            !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end)
            !$ACC LOOP SEQ
            DO k = k_start,k_end
              !$ACC LOOP GANG VECTOR PRIVATE(jc, aerncn, dummy)
              DO ic  = 1,n_ind(k)
                jc        = indices(ic,k)
                aerncn = 1.0E-6_wp*rho(jc,k)*( cams5(jc,k)/4.72911E-16_wp + cams6(jc,k)/1.55698E-15_wp )
                CALL ice_nucleation ( t(jc,k), aerncn=aerncn , znin=dummy )
                ncn(jc,k) = dummy
              ENDDO
            ENDDO
            !$ACC END PARALLEL
            !$ACC END DATA
          CASE (iRadAeroTegen) ! use Tegen dust with DeMott formula
            !$ACC DATA PRESENT(n_ind, indices, ncn, t, z_ifc, aer_dust)
            !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end)
            !$ACC LOOP SEQ
            DO k = k_start,k_end
              !$ACC LOOP GANG VECTOR PRIVATE(jc, aerncn, dummy)
              DO ic  = 1,n_ind(k)
                jc        = indices(ic,k)
                CALL ncn_from_tau_aerosol_speccnconst_dust (z_ifc(jc,k), z_ifc(jc,k+1), aer_dust(jc), aerncn)
                CALL ice_nucleation ( t(jc,k), aerncn=aerncn , znin=dummy )
                ncn(jc,k) = dummy
              ENDDO
            ENDDO
            !$ACC END PARALLEL
            !$ACC END DATA
          CASE DEFAULT
            CALL finish('mo_reff_main', 'icpl_aero_ice = 1 only available for irad_aero = 6,7,8.')
        END SELECT
      ELSE ! FR: Cooper (1986) used by Greg Thompson(2008)
        ! Some constant coefficients
        IF( lsuper_coolw) THEN
          znimax = znimax_Thom         !znimax_Thom = 250.E+3_wp,
        ELSE
          znimax = 150.E+3_wp     ! from previous ICON code
        END IF

        !$ACC DATA PRESENT(indices, ncn, t, n_ind)
        !$ACC PARALLEL ASYNC(1)
        !$ACC LOOP SEQ
        DO k = k_start,k_end
          !$ACC LOOP GANG VECTOR PRIVATE(jc, dummy)
          DO ic  = 1,n_ind(k)
            jc =  indices(ic,k)
             CALL ice_nucleation ( t(jc,k), znin=dummy )
             ncn(jc,k) = MIN(dummy,znimax)
           END DO
        END DO
        !$ACC END PARALLEL
        !$ACC END DATA

      ENDIF

    CASE (2,4) ! Rain, Graupel done with 1 mom param. w. fixed N0

#ifdef _OPENACC
      CALL finish('one_moment_calculate_ncn:','CASE hydrometeor=2,4 not available on GPU')
#endif

      well_posed = PRESENT(rho) .AND. ASSOCIATED(q)
      IF (.NOT. well_posed) THEN
        WRITE (message_text,*) 'Reff: Rho ans q  needs to be provided to one_mom_calculate_ncn'
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF

      DO k = k_start,k_end
        DO ic  = 1,n_ind(k)
          jc = indices(ic,k)
          ncn(jc,k) = reff_calc%ncn_coeff(1) * EXP ( reff_calc%ncn_coeff(2) * LOG( rho(jc,k)*q(jc,k)  ) )
        END DO
      END DO

    CASE (3) ! Snow, complex parameterization of N0 (copy paste and adapt from graupel)

#ifdef _OPENACC
      CALL finish('one_moment_calculate_ncn:','CASE hydrometeor=3 not available on GPU')
#endif

      well_posed = PRESENT(rho) .AND. ASSOCIATED(q) .AND. PRESENT(t)
      IF (.NOT. well_posed) THEN
        WRITE (message_text,*) 'Reff: Rho ans q  needs to be provided to one_mom_calculate_ncn->snow'
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF

      DO k = k_start,k_end
        DO ic  = 1,n_ind(k)
          jc = indices(ic,k)

          IF (isnow_n0temp == 1) THEN
            ! Calculate n0s using the temperature-dependent
            ! formula of Field et al. (2005)
            ztc = t(jc,k) - t0
            ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
            zn0s = zn0s1*EXP(zn0s2*ztc)
            zn0s = MIN(zn0s,1e9_wp)
            zn0s = MAX(zn0s,1e6_wp)
          ELSEIF (isnow_n0temp == 2) THEN
            ! Calculate n0s using the temperature-dependent moment
            ! relations of Field et al. (2005)
            ztc = t(jc,k) - t0
            ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)

            nnr  = 3._wp
            hlp = mma(1) + mma(2)*ztc + mma(3)*nnr + mma(4)*ztc*nnr &
                 & + mma(5)*ztc**2 + mma(6)*nnr**2 + mma(7)*ztc**2*nnr &
                 & + mma(8)*ztc*nnr**2 + mma(9)*ztc**3 + mma(10)*nnr**3
            alf = EXP(hlp*zlog_10) ! 10.0_wp**hlp
            bet = mmb(1) + mmb(2)*ztc + mmb(3)*nnr + mmb(4)*ztc*nnr &
                 & + mmb(5)*ztc**2 + mmb(6)*nnr**2 + mmb(7)*ztc**2*nnr &
                 & + mmb(8)*ztc*nnr**2 + mmb(9)*ztc**3 + mmb(10)*nnr**3

            ! assumes bms=2.0
            m2s = q(jc,k) * rho(jc,k) / ageo_snow
            m3s = alf*EXP(bet*LOG(m2s))

            hlp  = zn0s1*EXP(zn0s2*ztc)
            zn0s = 13.50_wp * m2s * (m2s / m3s) **3
            zn0s = MAX(zn0s,0.5_wp*hlp)
            zn0s = MIN(zn0s,1e2_wp*hlp)
            zn0s = MIN(zn0s,1e9_wp)
            zn0s = MAX(zn0s,1e6_wp)
          ELSE
            ! Old constant n0s
            zn0s = 8.0e5_wp
          ENDIF

          ncn(jc,k) =  reff_calc%ncn_coeff(1) * EXP( reff_calc%ncn_coeff(3)*LOG(zn0s)) &
               & * EXP (  reff_calc%ncn_coeff(2) * LOG( rho(jc,k)*q(jc,k)  ) )
        END DO
      END DO

    END SELECT

  END SUBROUTINE one_mom_calculate_ncn


! Init parameters for one effective radius calculation
  SUBROUTINE init_reff_calc (  reff_calc, hydrometeor, grid_scope, microph_param, &
                      &        p_q,p_reff,                                        &
                      &        return_fct,                                        &
                      &        p_qtot, p_ncn3D, p_ncn2D,                          &
                      &        ncn_param, dsd_type,    reff_param,                &
                      &        x_min, x_max, mu, nu, r_max, r_min  )

    ! Output
    TYPE(t_reff_calc), INTENT(INOUT)            :: reff_calc     ! Reff calculation parameters and pointers     
    LOGICAL,           INTENT(INOUT)            :: return_fct    ! Return code if some param was found (.true.)

    ! Obligatory parameters
    INTEGER,           INTENT(IN)               :: hydrometeor   ! Hydrometeor type
    INTEGER,  INTENT(IN)                        :: grid_scope    ! Total/Grid/Subgrid
    INTEGER,  INTENT(IN)                        :: microph_param ! Microphysics Parameterization

    ! Obligatory fields
    REAL(wp), DIMENSION(:,:,:),TARGET           :: p_q           ! Pointer to mixing ratio of hydrometeor
    REAL(wp), DIMENSION(:,:,:),TARGET           :: p_reff        ! Pointer to effective radius output   

    ! Extra fields
    REAL(wp), DIMENSION(:,:,:),TARGET, OPTIONAL :: p_qtot        ! Pointer to total (grid+subgrid) mixing ratio
    REAL(wp), DIMENSION(:,:,:),TARGET, OPTIONAL :: p_ncn3D       ! Pointer to 3D hydro. condensation nuclei
    REAL(wp), DIMENSION(:,:),  TARGET, OPTIONAL :: p_ncn2D       ! Pointer to 2D surface hydro. cond. nuc. 

    ! These parameters are needed by some param
    INTEGER, OPTIONAL,  INTENT(IN)              :: ncn_param     ! Parameterization for the number density
    INTEGER, OPTIONAL,  INTENT(IN)              :: dsd_type      ! Assumed Droplet Size Distribution. 
    INTEGER, OPTIONAL,  INTENT(IN)              :: reff_param    ! Parameterization type of reff

    ! These parameters should be set by micro, but the can also be overwritten by the user
    REAL(wp), OPTIONAL, INTENT(IN)              :: x_min         ! Min particle mass admited by the parameterization
    REAL(wp), OPTIONAL, INTENT(IN)              :: x_max         ! Max particle mass admited by the parameterization
    REAL(wp), OPTIONAL, INTENT(IN)              :: r_min         ! Min mass/radius admited by the parameterization
    REAL(wp), OPTIONAL, INTENT(IN)              :: r_max         ! Max mass/radius admited by the parameterization
    REAL(wp), OPTIONAL, INTENT(IN)              :: mu            ! Given Gamma parameter in DSD (only for dsd_type=2)
    REAL(wp), OPTIONAL, INTENT(IN)              :: nu            ! Given Nu parameter in DSD    (only for dsd_type=2)

    ! End of subroutine variable declaration 
    


    REAL(wp)                                    :: bf            ! Increase in reff due to DSD broadening

    ! Check input return_fct
    IF (.NOT. return_fct) THEN
      WRITE (message_text,*) 'Reff: Function init_reff_calc entered with previous error'
      CALL message('',message_text)
      RETURN
    END IF

    ! Allocate memory and nullify pointers
    CALL reff_calc%construct()

    ! Fill the type with initation
    reff_calc%hydrometeor   = hydrometeor   
    reff_calc%microph_param = microph_param
    reff_calc%grid_scope    = grid_scope

    ! Init standard coeffients
    reff_calc%mu            = -999.0
    reff_calc%nu            = -999.0
    reff_calc%r_min         = 1.e-6_wp                                            ! Minimum radius (1 mum)
    reff_calc%r_max         = 1.e-1_wp                                            ! Maximum radius (10cm)
    reff_calc%x_min         = 4.0_wp/3.0_wp*pi*rhoh2o* (reff_calc%r_min)**3.0_wp
    reff_calc%x_max         = 4.0_wp/3.0_wp*pi*rhoh2o* (reff_calc%r_max)**3.0_wp
    reff_calc%reff_param    = 0                                                   ! Spheroids
    reff_calc%dsd_type      = 0                                                   ! Consistent with param
    reff_calc%ncn_param     = 0                                                   ! Constant incloud-number

    ! Set pointers
    IF(PRESENT(p_qtot))     reff_calc%p_qtot=>p_qtot
    IF(PRESENT(p_ncn3D))    reff_calc%p_ncn3D=>p_ncn3D
    IF(PRESENT(p_ncn2D))    reff_calc%p_ncn2D=>p_ncn2D



    reff_calc%p_q=>p_q
    reff_calc%p_reff=>p_reff


    
    IF (PRESENT(dsd_type) )  reff_calc%dsd_type    = dsd_type
    IF (PRESENT(mu))         reff_calc%mu          = mu 
    IF (PRESENT(nu))         reff_calc%nu          = nu 
    IF (PRESENT(reff_param)) reff_calc%reff_param  = reff_param 
    IF (PRESENT(ncn_param))  reff_calc%ncn_param   = ncn_param 


    ! Consistency checks
    IF      ( (( reff_calc%ncn_param >= 4 .AND. reff_calc%ncn_param <= 7) .OR. reff_calc%ncn_param == 9 ) &
     &  .AND. ( .NOT. PRESENT(p_ncn3D) ) ) THEN
      ! Error: ncn pointer needs to be provided for these parameterizations
      WRITE (message_text,*) 'Error in reff: pointer to hydrometor number concentration needs to be provided'
      CALL message('',message_text)
      return_fct = .false.
      RETURN
    END IF
 
    IF (        (  (reff_calc%ncn_param >= 1    ) .AND. (reff_calc%ncn_param <= 3    ) ) .AND. &
     &    .NOT. (  (reff_calc%microph_param >= 1) .AND. (reff_calc%microph_param <= 3) ) .AND. &
                    reff_calc%hydrometeor >= 2  )  THEN
      ! Error: 1-mom ncn is only allowed with conistent param  for rain, graupel, snow   
      WRITE (message_text,*)       'Error in reff: the ncn 1 moment parameterization only runs with &
                                    &1 moment microphysical scheme'
      CALL message('',message_text)
      return_fct = .false.
      RETURN
    END IF



    IF (  reff_calc%dsd_type == 2 .AND. ( reff_calc%mu < -900.0_wp .OR.  reff_calc%nu < -900.0_wp )  ) THEN
      ! Error: parameters needed for predefined DSD      
      WRITE (message_text,*) 'Error in reff: Insufficent parameters to initiate reff calculations for the choosen DSD'
      CALL message('',message_text)
      return_fct = .false.
      RETURN
    END IF

    ! Grid scale does not neccesarily needs total quantities, but it is then the same as total
    IF ( (reff_calc%grid_scope == 1) .AND. (.NOT. ASSOCIATED(reff_calc%p_qtot )) ) THEN 
      reff_calc%grid_scope = 0
    END IF

    ! Grid scale does not neccesarily needs total quantities, but it is then the same as total
    IF ( (reff_calc%grid_scope == 2) .AND. (.NOT. ASSOCIATED(reff_calc%p_qtot )) ) THEN 
      ! Error: total fields needed for subgrid-scale radius
      WRITE (message_text,*)       'Error in reff: a total field is needed to caculate subgrid effective radius'
      CALL message('',message_text)
      return_fct = .false.
      RETURN
    END IF
    

! -----------------------------
! Calculate coefficients
! -----------------------------

    SELECT CASE ( microph_param ) ! Choose which microphys scheme

    CASE (1,2,3)        ! One-Moment schemes
      CALL  one_mom_reff_coefficients( reff_calc,return_fct )  
      IF (.NOT. return_fct) THEN
          WRITE (message_text,*) 'Error in init reff: the 1 mom scheme could not initiate coefficients. Check options'
          CALL message('',message_text)
          return_fct = .false.
          RETURN
      END IF

    CASE (4,5,6,7,9)      ! Two-Moment schemes
      CALL  two_mom_reff_coefficients( reff_calc,return_fct )  
      IF (.NOT. return_fct) THEN
        WRITE (message_text,*) 'Error in init reff: the 2 mom scheme could not initiate coefficients. Check options'
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF

    CASE (8)      ! sbm scheme
      CALL  two_mom_reff_coefficients( reff_calc,return_fct )
      IF (.NOT. return_fct) THEN
        WRITE (message_text,*) 'Error in init reff: the SBM-2M Piggybacking scheme could not initiate coefficients. Check options'
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF

    CASE (101)             ! RRTM Param.
      SELECT CASE ( hydrometeor ) 
      CASE(0)  ! Cloud water
        ! Base is monodisperse. Broadening factor in calculations because depends on coeff.
        CALL  reff_coeff_monodisperse_spherical (reff_calc )   ! RRTM Parameterization
        reff_calc%r_min         = 2.e-6_wp  ! Minimum radius
        reff_calc%r_max         = 32.e-6_wp ! Maximum radius

      CASE(1)  !Ice, see ECHAM5 documentation (Roeckner et al, MPI report 349)
        reff_calc%reff_coeff(1) = 83.8e-6_wp
        reff_calc%reff_coeff(2) = 0.216_wp
        reff_calc%ncn_param     = -1 ! No ncn parameteriyation is needed

        ! Extra limits added for eccrad
        reff_calc%r_min         = 4.e-6_wp  ! Minimum radius 
        reff_calc%r_max         = 99.e-6_wp ! Maximum radius
       CASE DEFAULT
         WRITE (message_text,*) 'Error in init reff: RRTM is only defined for cloud and ice (no rain, graupel...)'
         CALL message('',message_text)
         return_fct = .false.
         RETURN
       END SELECT

     CASE (100)    ! Spherical liquid particles
      ! Base is monodisperse
       CALL  reff_coeff_monodisperse_spherical (reff_calc )                            

       IF ( reff_calc%dsd_type == 2) THEN  ! Polydisperse
         ! Broadening due to choosing a radial gamma distribution with fixed gamma, nu 
         bf = GAMMA ( (nu + 4.0_wp)/mu ) / GAMMA ( (nu + 3.0_wp)/mu ) * &
              & ( GAMMA ( (nu + 1.0_wp)/mu ) / GAMMA ( (nu + 4.0_wp)/mu ) )**(1.0_wp/3.0_wp)
         reff_calc%reff_coeff(1) = reff_calc%reff_coeff(1) * bf
        
       END IF

     END SELECT

     ! Overwrite xmin and xmax if present
     IF ( PRESENT (x_min) )  reff_calc%x_min = x_min
     IF ( PRESENT (x_max) )  reff_calc%x_max = x_max
     IF ( PRESENT (r_min) )  reff_calc%r_min = r_min
     IF ( PRESENT (r_max) )  reff_calc%r_max = r_max


     ! Select if ncn parameterization is incloud or grid scale
     IF (PRESENT(ncn_param)) THEN
       SELECT CASE ( reff_calc%ncn_param )   ! Select NCN parameterization
        CASE (1,2,3) ! 1 momment microphysics
          SELECT CASE ( hydrometeor ) 
          CASE (0,1)   ! Cloud water or ice
            reff_calc%ncn_param_incloud = 1  ! All NCN param provides incloud values
          CASE (2,3,4)
            reff_calc%ncn_param_incloud = 0  ! Grid scale values for graupel, snow, rain      
          END SELECT
        CASE (4,5,6,7,8,9)  ! 2 Moment/SBM microphysics
          reff_calc%ncn_param_incloud = 0  ! Grid scale values for all param.        

        CASE DEFAULT
          reff_calc%ncn_param_incloud = 1  ! Default params. are incloud
        END SELECT

      END IF

      !$ACC UPDATE DEVICE(reff_calc%reff_coeff) ASYNC(1)
    
  END SUBROUTINE init_reff_calc




!------------------------------------------------------------------------------------------------------------

! Coefficients for monodisperse spheres. 
    SUBROUTINE reff_coeff_monodisperse_spherical (reff_calc ) 
      
      TYPE(t_reff_calc) , INTENT(INOUT)   :: reff_calc       ! Reff calculation parameters and pointers     

      REAL(wp)                            :: a,b             ! Geometric factors

      ! Geometric factors x = a D**[b]
      a                       = pi/6.0_wp * rhoh2o
      b                       = 3.0_wp
      reff_calc%reff_coeff(1) = a**(-1.0_wp/b)
      reff_calc%reff_coeff(2) = 1.0_wp/b

    END SUBROUTINE reff_coeff_monodisperse_spherical

!------------------------------------------------------------------------------------------------------------



!! Calculte running indices for reff and n of a parameterization differentiating between grid and subgrid
    SUBROUTINE mapping_indices ( indices, n_ind, reff_calc, k_start, k_end, is, ie, jb, return_fct )


      INTEGER (KIND=i4), INTENT(INOUT)   ::     indices(:,:) ! Mapping for going through array
      INTEGER (KIND=i4), INTENT(INOUT)   ::     n_ind(:)     ! Number of indices for each k level
      TYPE(t_reff_calc), INTENT(IN)      ::     reff_calc    ! Reff calculation parameters and pointers

      INTEGER, INTENT(IN)                ::     k_start, k_end, is, ie    ! Start, end total indices
      INTEGER, INTENT(IN)                ::     jb            ! Domain index
      LOGICAL, INTENT(INOUT)             ::     return_fct    ! Function return. .true. for right

      ! End of subroutine variable declaration

      REAL(wp) ,POINTER, DIMENSION(:,:)  ::     q            ! Mixing ratio of hydrometeor
      REAL(wp) ,POINTER, DIMENSION(:,:)  ::     q_tot        ! Mixing ratio of hydrometeor
                                                             ! From ICON Scientific Documentaion (Cloud scheme)
      REAL(wp), PARAMETER                ::     qmin = 1E-6  ! Difference between cloud/nocloud in kg/kg
      REAL(wp), PARAMETER                ::     qsub = 1E-6  ! Difference between grid/subgrid in kg/kg

      INTEGER                            ::     k, jc        ! Counters
      INTEGER, DIMENSION(ie,k_end)       ::     llq          ! logical conditions for cloud/no cloud and grid/subgrid


    ! Check input return_fct
      IF (.NOT. return_fct) THEN
        WRITE (message_text,*) 'Reff: Function init_reff_calc entered with previous error'
        CALL message('',message_text)
        RETURN
      END IF

      !$ACC ENTER DATA CREATE(llq)

      ! Initialize inidices
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(k_end, ie)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO k = 1,k_end
        DO jc = 1,ie
          indices(jc,k) = 0
          llq(jc,k) = 0
        END DO
      END DO
      !$ACC END PARALLEL

      SELECT CASE ( reff_calc%grid_scope )

      CASE (0) ! Total parameterization ( no differentation grid/subgrid)

        ! Use total if available
        IF ( ASSOCIATED(reff_calc%p_qtot)) THEN 
          q=>reff_calc%p_qtot(:,:,jb)
        ELSE  ! In case grid scale only or no total available
          q=>reff_calc%p_q(:,:,jb)
        END IF

        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(k_start, k_end, is, ie)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO k = k_start,k_end
          DO jc = is, ie
            IF (q(jc,k) > qmin) THEN
              llq(jc,k) = 1
            ENDIF
          END DO
        END DO
        !$ACC END PARALLEL

        CALL generate_index_list_batched(llq, indices, 1, ie, n_ind, 1)

      CASE (1) ! Only grid scale (with same subgrid/grid criteria as subgrid)
        
        IF ( ASSOCIATED(reff_calc%p_qtot ) ) THEN
          q_tot=>reff_calc%p_qtot(:,:,jb)
          q=>reff_calc%p_q(:,:,jb)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(k_start, k_end, is, ie)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO k = k_start,k_end
            DO jc = is, ie
              IF ((q(jc,k) > qsub) .AND. (q_tot(jc,k) > qmin)) THEN
                llq(jc,k) = 1
              ENDIF
            END DO
          END DO
          !$ACC END PARALLEL

          CALL generate_index_list_batched(llq, indices, 1, ie, n_ind, 1)

        ELSE
          WRITE (message_text,*) 'Warning: Reff does not have information for generate inidices for subgrid ncn'
          CALL message('',message_text)
          return_fct = .false.              
          RETURN
        END IF

      CASE (2) ! Only subgrid scale
        
        IF ( ASSOCIATED(reff_calc%p_qtot ) ) THEN
          q_tot=>reff_calc%p_qtot(:,:,jb)
          q=>reff_calc%p_q(:,:,jb)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) FIRSTPRIVATE(k_start, k_end, is, ie)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO k = k_start,k_end
            DO jc = is, ie
              IF ((q_tot(jc,k) > qmin) .AND. (q(jc,k) < qsub)) THEN
                llq(jc,k) = 1
              ENDIF
            END DO
          END DO
          !$ACC END PARALLEL

          CALL generate_index_list_batched(llq, indices, 1, ie, n_ind, 1)

        ELSE
          WRITE (message_text,*) 'Warning: Reff does not have information for generate inidices for subgrid ncn'
          CALL message('',message_text)
          return_fct = .false.              
          RETURN
        END IF

    END SELECT

    !$ACC WAIT
    !$ACC EXIT DATA DELETE(llq)

  END SUBROUTINE mapping_indices


! -----------------------------------------------------------------------------------------------------------
  
  !! Calculate reff based on the parameters and ncn and indices previosly calculated
  SUBROUTINE calculate_reff ( reff_calc, indices, n_ind, rho, k_start,      &
                            & k_end, jb, return_fct, ncn, clc, fr_gl, fr_land )

    TYPE(t_reff_calc)  ,INTENT(INOUT)    ::    reff_calc        ! Reff calculation parameters and pointers
    INTEGER (KIND=i4)  ,INTENT(IN)       ::    indices(:,:)     ! Mapping for going through array
    INTEGER (KIND=i4)  ,INTENT(IN)       ::    n_ind(:)         ! Number of indices for each k level
    REAL(wp) ,INTENT(IN)                 ::    rho(:,:)         ! Densityof air

    INTEGER ,INTENT(IN)                  ::    k_start, k_end   ! Start, end total indices
    INTEGER ,INTENT(IN)                  ::    jb               ! Domain
    LOGICAL ,INTENT(INOUT)               ::    return_fct       ! Function return. .true. for right

    REAL(wp) ,INTENT(INOUT), OPTIONAL    ::    ncn(:,:)         ! Number concentration
    REAL(wp) ,INTENT(IN),OPTIONAL        ::    clc(:,:)         ! Cloud fraction
    REAL(wp) ,INTENT(IN),OPTIONAL        ::    fr_gl(:)         ! Fraction of glaciers (for RRTM)
    REAL(wp) ,INTENT(IN),OPTIONAL        ::    fr_land(:)       ! Fraction of land (for RRTM)


    ! End of subroutine variable declarations

    REAL(wp) ,POINTER, DIMENSION(:,:)    ::     reff            ! Pointer to effective radius
    REAL(wp) ,POINTER, DIMENSION(:,:)    ::     q               ! Pointer to mixing ratio
    INTEGER                              ::     k,ic,jc         ! Counters
    REAL(wp)                             ::     x, x_max,x_min  ! Mean mass of particle, maximum,minimum
    REAL(wp)                             ::     r_min, r_max    ! Minimum and maximum radius 
    REAL(wp)                             ::     bf              ! Broadening factor of DSD (for RRTM)
    REAL(wp) ,PARAMETER                  ::     eps = 1.0e-8    ! Epsilon constant

    ! RRTM Parameters (Author of the original code: Bjorn Stevens, MPI-M, Hamburg)
    REAL (wp), PARAMETER                 ::    &
                                &  zkap_cont = 1.143_wp, &      ! continental (Martin et al.JAS 1994? ) breadth param
                                &  zkap_mrtm = 1.077_wp         !  maritime (Martin et al.) breadth parameter

    ! Check input return_fct
    IF (.NOT. return_fct) THEN
      WRITE (message_text,*) 'Reff: Function calculate_reff entered with previous error'
      CALL message('',message_text)
      RETURN
    END IF

    reff=>reff_calc%p_reff(:,:,jb)

    SELECT CASE ( reff_calc%grid_scope ) ! Select subgrid or grid field
    CASE (0) ! Grid and subgrid
      IF ( ASSOCIATED(reff_calc%p_qtot)) THEN
        q=>reff_calc%p_qtot(:,:,jb)
      ELSE
        q=>reff_calc%p_q(:,:,jb)
      END IF
    CASE (1)
      q=>reff_calc%p_q(:,:,jb)
    CASE(2)
      q=>reff_calc%p_qtot(:,:,jb)
    END SELECT

    ! Translate ncn values from incloud to grid-scale values
    IF ( reff_calc%ncn_param_incloud == 1 .AND. PRESENT(ncn) .AND. reff_calc%ncn_param >= 0 ) THEN 
      !$ACC DATA PRESENT(indices, ncn, clc, n_ind)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end)
      !$ACC LOOP SEQ
      DO k = k_start,k_end
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic  = 1,n_ind(k)
          jc        =  indices(ic,k)
          ncn(jc,k) =  ncn(jc,k)*clc(jc,k) 
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA
    END IF


    SELECT CASE ( reff_calc%microph_param ) ! Choose which microphys param

    CASE (0,1,2,3,4,5,6,7,8,9,100)      ! Currently all cases except for RRTM follow same scheme as function of mean mass
      x_max = reff_calc%x_max
      x_min = reff_calc%x_min

      SELECT CASE (reff_calc%reff_param )        
      CASE (0)    !Spheroid  : reff = 0.5*c1*x**c2 (x= mean mass)
        !$ACC DATA PRESENT(indices, ncn, n_ind, rho, q, reff_calc, reff, reff_calc%reff_coeff(2))
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end, x_max, x_min)
        !$ACC LOOP GANG VECTOR PRIVATE(jc, x)
        DO k = k_start,k_end
          DO ic  = 1,n_ind(k)
            jc         =  indices(ic,k)
            x          =  rho(jc,k)* q(jc,k) / ( ncn(jc,k) + eps )
            x          =  MAX( MIN( x,x_max),x_min)
            reff(jc,k) =  0.5_wp* reff_calc%reff_coeff(1) * EXP( reff_calc%reff_coeff(2) * LOG( x ) )
          END DO
        END DO
        !$ACC END PARALLEL
        !$ACC END DATA

      !Fu Needles: reff= c5/(c1*x**c2 + c3*x**c4) 
      ! Here c5=0.5 fixed (different from libRadtran documentation c5=3*sqrt(3)/8=0.65)
      CASE (1)  
        !$ACC DATA PRESENT(indices, ncn, n_ind, rho, q, reff_calc, reff, reff_calc%reff_coeff(4))
        !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end, x_max, x_min)
        !$ACC LOOP SEQ
        DO k = k_start,k_end
          !$ACC LOOP GANG VECTOR PRIVATE(jc, x)
          DO ic  = 1,n_ind(k)
            jc         =  indices(ic,k)
            x          =  rho(jc,k)* q(jc,k) / ( ncn(jc,k) + eps )                
            x          =  MAX( MIN( x,x_max),x_min)
            reff(jc,k) =  0.5_wp/( reff_calc%reff_coeff(1) * EXP( reff_calc%reff_coeff(2) * LOG( x ) ) + &
                        & reff_calc%reff_coeff(3) * EXP( reff_calc%reff_coeff(4) * LOG( x ) ) )
          END DO
        END DO
        !$ACC END PARALLEL
        !$ACC END DATA
      END SELECT

    CASE(101)  ! RRTM

#ifdef _OPENACC
      CALL finish('calculate_reff:','CASE microph_param=101 not available on GPU')
#endif

      r_max = reff_calc%r_max
      r_min = reff_calc%r_min
      
      SELECT CASE (reff_calc%hydrometeor )
 
      CASE (0)  ! Cloud water
        DO k = k_start,k_end
          DO ic  = 1,n_ind(k)
            jc         =  indices(ic,k)
            x          =  rho(jc,k)* q(jc,k) / ( ncn(jc,k) + eps )                
            ! Broadening factor depending on sea-land
            bf         =  zkap_cont*(fr_land(jc)-fr_gl(jc)) + zkap_mrtm*(1.0_wp-fr_land(jc)+fr_gl(jc))
            reff(jc,k) =  0.5_wp* bf*reff_calc%reff_coeff(1) * EXP( reff_calc%reff_coeff(2) * LOG( x ) )
            reff(jc,k) =  MAX( MIN(reff(jc,k) ,r_max),r_min)
          END DO
        END DO
        
      CASE (1)  !Ice,  see ECHAM5 documentation (Roeckner et al, MPI report 349)
        DO k = k_start,k_end
          DO ic  = 1,n_ind(k)
            jc         =  indices(ic,k)
            x          =  rho(jc,k)* q(jc,k)*1000.0_wp / ( clc(jc,k) + eps ) ! There is no N dependency
            reff(jc,k) =  reff_calc%reff_coeff(1) * EXP ( reff_calc%reff_coeff(2)* LOG( x ))
            reff(jc,k) =  MAX( MIN( reff(jc,k) ,r_max),r_min)             
          END DO
        END DO
        
      END SELECT

    END SELECT


  END SUBROUTINE calculate_reff

! --------------------------------------------------------------------------------------------------------


  !! Calculte number concentraion of a hydrometeor
  SUBROUTINE calculate_ncn( ncn, reff_calc, indices, n_ind , k_start, k_end ,jb, rho, t,  &
       &                    icpl_aero_ice, z_ifc, cams5, cams6, aer_dust,                 &
       &                    return_fct )

    REAL(wp)          , INTENT(INOUT), DIMENSION(:,:) :: ncn             ! Number concentration
    TYPE(t_reff_calc) , INTENT(IN)                    :: reff_calc       ! Reff calculation parameters and pointers
    INTEGER (KIND=i4) , INTENT(IN)   , DIMENSION(:,:) :: indices         ! Mapping for going through array
    INTEGER (KIND=i4) , INTENT(IN)   , DIMENSION(:)   :: n_ind           ! Number of indices for each k level

    INTEGER           , INTENT(IN)                    :: k_start, k_end  ! Start, end total indices
    INTEGER           , INTENT(IN)                    :: jb              ! Domain index
    INTEGER           , INTENT(IN)                    :: icpl_aero_ice   ! aerosols ice nucleation scheme

    REAL(wp), INTENT(IN)        , DIMENSION(:,:)       :: rho             ! Density of air
    REAL(wp), INTENT(IN)        , DIMENSION(:,:)       :: t               ! Temperature
    REAL(wp), INTENT(IN), POINTER, DIMENSION(:,:)      :: cams5, cams6    ! CAMS dust mixing ratios
    REAL(wp), INTENT(IN)        , DIMENSION(:,:)       :: z_ifc           ! height at interface levels
    REAL(wp), INTENT(IN), POINTER, DIMENSION(:)        :: aer_dust        ! Tegen dust total column mass
    
    LOGICAL           , INTENT(INOUT)                 :: return_fct      ! Function return. .true. for right
    
    
    ! End of subroutine variable declarations

    REAL(wp), POINTER                , DIMENSION(:)   :: surf_cloud_num  ! Number concentration at surface (cloud_num)
    REAL(wp), POINTER                , DIMENSION(:,:) :: space_cloud_num ! Number concentration
    REAL(wp), POINTER                , DIMENSION(:,:) :: q               ! Mixing ratio of hydrometeor
    INTEGER                                           :: k, ic, jc       ! Counters
    LOGICAL                                           :: well_posed      ! Logical check
!    REAL(wp)                                          :: aerncn          ! CAMS dust aerosols number concentration
    REAL(wp)                                          :: cloud_num

    ! Check input return_fct
    IF (.NOT. return_fct) THEN
      WRITE (message_text,*) 'Reff: Function calculate_ncn entered with previous error'
      CALL message('',message_text)
      RETURN
    END IF

    ! In case tot is avalaible, use it by default
    IF ( ASSOCIATED(reff_calc%p_qtot) .AND. reff_calc%grid_scope .NE. 1 ) THEN
      q=>reff_calc%p_qtot(:,:,jb)
    ELSE  ! In case grid scale only or no total available
      q=>reff_calc%p_q(:,:,jb)
    END IF

    IF ( ASSOCIATED(reff_calc%p_ncn3D ) ) THEN
      space_cloud_num=>reff_calc%p_ncn3D(:,:,jb)
    END IF

    IF ( ASSOCIATED(reff_calc%p_ncn2D ) ) THEN
      surf_cloud_num=>reff_calc%p_ncn2D(:,jb)
    END IF

    ncn = 0.0

    SELECT CASE ( reff_calc%ncn_param ) ! Choose which microphys param

    CASE (0)      ! Constant number. Use cloud_num variable.
#ifdef _OPENACC
      CALL finish('calculate_ncn:','CASE ncn_param=0 not available on GPU')
#endif
      CALL get_cloud_number(cloud_num)
      DO k = k_start,k_end
        DO ic  = 1,n_ind(k)
          jc        = indices(ic,k)
          ncn(jc,k) = cloud_num
        END DO
      END DO

    CASE (1,2,3)  ! 1 mom microphysics


      SELECT CASE ( reff_calc%hydrometeor)
      CASE (0) ! Cloud water. It currently uses cloud_num 2D for the cloud water.
        IF (ASSOCIATED(reff_calc%p_ncn2D)) THEN
          ! Cloud_num field
          CALL one_mom_calculate_ncn( ncn, return_fct, reff_calc, k_start,             &
               &                      k_end, indices, n_ind,                           &
               &                      icpl_aero_ice, cams5, cams6, z_ifc, aer_dust, q, &
                                      surf_cloud_num = surf_cloud_num )
        ELSE
          ! Constant cloud_num
          CALL one_mom_calculate_ncn( ncn, return_fct, reff_calc, k_start,              &
               &                      k_end, indices, n_ind,                            &
               &                      icpl_aero_ice, cams5, cams6, z_ifc, aer_dust, q)
        END IF

      CASE DEFAULT
        well_posed = ASSOCIATED(reff_calc%p_q)
        IF (.NOT. well_posed) THEN
          WRITE (message_text,*) 'Reff: Insufficient arguments to call calculate ncn from 1 moment scheme'
          CALL message('',message_text)
          return_fct = .false.
          RETURN
        END IF
        CALL one_mom_calculate_ncn( ncn, return_fct, reff_calc, k_start,             &
             &                      k_end, indices, n_ind,                           &
             &                      icpl_aero_ice, cams5, cams6, z_ifc, aer_dust, q, &
             &                      t = t, rho =rho) 
      END SELECT


      IF (.NOT. return_fct) THEN
        WRITE (message_text,*) 'Error in init reff: the 1 mom scheme can not calculate the ncn. Check options'
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF


    CASE (4,5,6,7,8,9,101) ! Use acdnc or other field (from radiation)
      well_posed = ASSOCIATED(reff_calc%p_ncn3D)
      IF (.NOT. well_posed) THEN
        WRITE (message_text,*) 'Reff: A 3D clound number field (cdnc/qn) needs to be provided '
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF

      !$ACC DATA PRESENT(n_ind, indices, ncn, space_cloud_num)
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end)
      !$ACC LOOP SEQ
      DO k = k_start,k_end
        !$ACC LOOP GANG VECTOR PRIVATE(jc)
        DO ic  = 1,n_ind(k)
          jc        = indices(ic,k)
          ncn(jc,k) = space_cloud_num(jc,k)
        END DO
      END DO
      !$ACC END PARALLEL
      !$ACC END DATA


    CASE (102) ! Use cloud_num (from radiation). This is current default for 1 mom microphysics.
#ifdef _OPENACC
      CALL finish('calculate_ncn:','CASE ncn_param=102 not available on GPU')
#endif
      well_posed = ASSOCIATED(reff_calc%p_ncn2D)
      IF (.NOT. well_posed) THEN
        WRITE (message_text,*) 'Reff: surf_cloud needs to be provided to calculate cloud number calculate_ncn'
        CALL message('',message_text)
        return_fct = .false.
        RETURN
      END IF

      DO k = k_start,k_end
        DO ic  = 1,n_ind(k)
          jc        = indices(ic,k)
          ncn(jc,k) = surf_cloud_num(jc)
        END DO
      END DO

    CASE(-1) ! No calculation of ncn, but also not error (because not neccesary)

    CASE DEFAULT

      WRITE (message_text,*) 'Reff: Insufficient arguments to call calculate ncn'
      CALL message('',message_text)
      return_fct = .false.
      RETURN

    END SELECT



  END SUBROUTINE calculate_ncn


! Combine two hydrometeors fields into one, keeping qtot/rtot = q1/r1 + q2/r2
  SUBROUTINE  combine_reff( q1, reff_1, q2, reff_2, clc, k_start, k_end, is, ie )

    REAL(wp)          , INTENT(INOUT)         :: q1(:,:)                ! Mass concentration of smaller hydromet. (also store results)
    REAL(wp)          , INTENT(INOUT)         :: reff_1(:,:)            ! Effective radius of smaller hydromet. (also store results)
    REAL(wp)          , INTENT(IN)            :: q2(:,:)                ! Mass concentration of larger hydromet (ususally not in cloud cover).
    REAL(wp)          , INTENT(IN)            :: reff_2(:,:)            ! Effective radius of larger hydromet.
    REAL(wp)          , INTENT(INOUT)         :: clc(:,:)               ! Modified cloud cover: It is set to 1 if reff_2 > 1e-5  and q2>qcrit_reff
    INTEGER           , INTENT(IN)            :: k_start, k_end, is, ie ! Start, end total indices

    REAL(wp)                                  :: q_ov_reff              ! Local cross section

    INTEGER                                   :: k, jc                  ! Local counters

    REAL(wp)          , PARAMETER             :: qcrit_reff = 5e-5

    !$ACC DATA PRESENT(q1, reff_1, q2, reff_2, clc)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end, is, ie)
    !$ACC LOOP GANG VECTOR COLLAPSE(2) PRIVATE(q_ov_reff)
    DO k = k_start,k_end
      DO jc = is,ie
        IF ( reff_2(jc,k) > 1e-5_wp .AND. q2(jc,k) > qcrit_reff) THEN  ! Combine only when there is something in second phase
          clc(jc,k) = 1.0_wp        ! Set cloud cover to 1.0 if there is somethin in the larger phase
          IF ( reff_1(jc,k) > 1e-6_wp)  THEN ! Also something in first phase
            q_ov_reff = q1(jc,k)/reff_1(jc,k) + q2(jc,k)/reff_2(jc,k)
            q1(jc,k)  = q1(jc,k) + q2(jc,k)
            IF ( q_ov_reff > 1E-6) THEN
              reff_1(jc,k) = q1(jc,k)/q_ov_reff
            ELSE
              q1(jc,k)     = 0.0_wp
              reff_1(jc,k) = 0.0_wp  ! Set to 0 micro, nominally for negligible extinction
            END IF
          ELSE  ! Something in second phase, but not in first
            q1(jc,k) = q2(jc,k)
            reff_1(jc,k) = reff_2(jc,k)
          END IF
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE combine_reff

! Set a maximum effective radius, keeping qend/rmax = qini/rini
  SUBROUTINE set_max_reff( q, reff, reff_max, k_start, k_end, is, ie )
    REAL(wp)          , INTENT(INOUT)         :: q(:,:)                 ! Mass concentration of hydromet. (also store results)
    REAL(wp)          , INTENT(INOUT)         :: reff(:,:)              ! Effective radius of hydromet. (also store results)
    REAL(wp)          , INTENT(IN)            :: reff_max               ! Maximum effective radius
    INTEGER           , INTENT(IN)            :: k_start, k_end, is, ie ! Start, end total indices    

    REAL(wp)                                  :: q_ov_reff     ! Local cross section
    INTEGER                                   :: k,jc           ! Local counters 
    
    !$ACC DATA PRESENT(q, reff)
    !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(k_start, k_end, is, ie, reff_max)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO k = k_start,k_end
      DO jc  = is,ie
        IF ( reff(jc,k) > reff_max) THEN
          q(jc,k)    = q(jc,k)*reff_max/reff(jc,k)
          reff(jc,k) = reff_max
        END IF
      END DO
    END DO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE set_max_reff

END MODULE mo_reff_main
