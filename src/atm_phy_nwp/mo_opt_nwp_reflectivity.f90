!NEC$ options "-finline-max-depth=3 -finline-max-function-size=1000"
!
! Routines for optional diagnostic radar reflectivity on the model grid in NWP
! (formerly located in mo_opt_nwp_diagnostics.f90)
!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_opt_nwp_reflectivity

  USE mo_kind,                  ONLY: wp
  USE mo_math_constants,        ONLY: pi
  USE gscp_data,                ONLY: isnow_n0temp, zami, mu_rain, zams_ci, zams_gr, zbms, &
    &                                 znimax_Thom, zthn, mma, mmb, zcnue
  USE mo_2mom_mcrph_main,       ONLY: init_2mom_scheme,      &
    &                                 rain_coeffs  ! contains the parameters for the mue-Dm-relation
  USE mo_2mom_mcrph_types,      ONLY: particle, particle_frozen
  USE mo_2mom_mcrph_processes,  ONLY: moment_gamma, rain_mue_dm_relation
  USE mo_exception,             ONLY: finish, message
  USE mo_fortran_tools,         ONLY: set_acc_host_or_device
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: compute_field_dbz_1mom
  PUBLIC :: compute_field_dbz_2mom

  !> module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_opt_nwp_reflectivity'

CONTAINS

  !>
  !! Calculate radar reflectivity for the 1-moment microphysics scheme in linear units mm^6/m^3
  !!
  SUBROUTINE compute_field_dbz_1mom( npr, nlev, nblks, startblk, endblk, jk_start,       &
                                     startidx1, endidx2,                                 &
                                     lmessage_light, lmessage_full, my_id_for_message,   &
                                     rho_w, rho_ice,                                     &
                                     K_w, K_ice, T_melt, igscp, q_crit_radar,            &
                                     T, rho, q_cloud, q_rain, q_ice, q_snow, z_radar,    &
                                     q_graupel, n_cloud_s, lacc )
 
   !------------------------------------------------------------------------------
    !
    ! Description:  Calculation of grid point values for effective radar
    !               reflectivity factor Z in dBZ.
    !
    ! Method:       Rayleigh-Approximation for the Back-Scattering
    !               (no attenuation!), Debye-Approximation for the
    !               effective refractive index of two-component ice-air-
    !               mixture particles (dry ice, snow, and graupel).
    !               Melting particles by substitution of the ice substance
    !               with water when T_air > 273.16 K. For the applied
    !               Rayleigh scattering, this is equal to assuming instantaneous
    !               melting, because only the square of the total water mass of the particle
    !               counts ( = the total square of its dipole moment).
    !               
    ! Inputs:       npr, nlev, nblks : field dimensions
    !               startblk, endblk, jk_start, startidx1, endidx2   : loop start and end indices
    !               rho_w        : bulk density of pure water [kg/m**3]
    !               rho_ice      : bulk density of pure ice   [kg/m**3]
    !               K_w          : dielectric constant of water
    !               K_ice        : dielectric constant of ice
    !               T_melt       : melting temperature of ice
    !               igscp        : itype of grid scale precip (1,2)
    !               q_crit_radar : threshold for the q's to compute reflectivity  [kg/m**3]
    !               T            : temperature field          [K]
    !               rho          : air density                [kg/m**3]
    !               q_cloud      : cloud water mixing ratio   [kg/kg] 
    !               q_rain       : rain water mixing ratio    [kg/kg] 
    !               q_ice        : cloud ice mixing ratio     [kg/kg]
    !               q_snow       : snow mixing ratio          [kg/kg]  
    !   OPTIONAL:   q_graupel    : graupel mixing ratio       [kg/kg] 
    !   OPTIONAL:   n_cloud_s    : surface value of cloud droplet number concentration [1/kg] 
    !
    ! Output:       z_radar      : 3D field of Z              [mm^6/m^3]
    ! 
    !
    !------------------------------------------------------------------------------
    
    ! Input/Output parameters:
    !-------------------------

    INTEGER , INTENT(IN)       :: npr, nlev, nblks
    INTEGER,  INTENT(IN)       :: startblk, endblk, jk_start, startidx1, endidx2, my_id_for_message
    LOGICAL,  INTENT(in)       :: lmessage_light, lmessage_full
    REAL(wp), INTENT(IN)       :: rho_w, rho_ice, K_w, K_ice, T_melt, q_crit_radar
    INTEGER , INTENT(IN)       :: igscp
    REAL(wp), INTENT(IN)       :: T(:,:,:),          &
                                  rho(:,:,:),        &
                                  q_cloud(:,:,:),    &
                                  q_rain(:,:,:),     &
                                  q_ice(:,:,:),      &
                                  q_snow(:,:,:)
    REAL(wp), INTENT(IN), OPTIONAL :: q_graupel(:,:,:), n_cloud_s(:,:)
    REAL(wp), INTENT(OUT)      :: z_radar(:,:,:)
    LOGICAL, OPTIONAL, INTENT(IN)  :: lacc              !< initialization flag

    ! Local Variables
    !----------------

    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_dbz_1mom'
    
    REAL(wp) ::  rho_c, rho_i, rho_r, rho_s, rho_g, z_cloud, z_rain, z_snow, z_ice, z_grau

    INTEGER        :: jc, jk, jb, i_startidx, i_endidx
    LOGICAL, SAVE  :: firstcall = .TRUE.
    logical        :: lqnc_input

    REAL(wp)       :: z_fac_ice_dry, z_fac_ice_wet
    REAL(wp)       :: ztc, m2s, m3s, alf, bet, hlp, zn0s, x_c, x_i_mono, n_i
    INTEGER        :: nn, pe_center


    REAL(wp), SAVE :: z_r, z_g, z_s, p_r, p_s, p_g, mom_fac, z_fac_c
    REAL(wp), SAVE :: nor, nog, amg, bmg, nos, ams, bms, &
                      ami, bmi, D_c_fix, x_c_fix, mue_rain_c

    REAL(wp), PARAMETER :: zn0s1       = 13.5_wp * 5.65E5_wp, & ! parameter in N0S(T)
                           zn0s2       = -0.107_wp,           & ! parameter in N0S(T), Field et al
                           eps         = 1.00E-15_wp,         &
                           n_cloud_min = 1.00E7_wp,           & ! min. allowed cloud number conc. in externally provided n_cloud_s
                           convfac     = 1.E18_wp,            &
                           z10olog10   = 10.0_wp / LOG(10.0_wp)

    TYPE(particle) :: cloud_tmp

    ! Variables for cloud ice parameterization 
    REAL(wp)                             ::  fxna               ! statement function for ice crystal number, Cooper(1986)
    REAL(wp)                             ::  fxna_cooper        ! statement function for ice crystal number, Cooper(1986)
    REAL(wp)                             ::  ztx, ztmelt        ! dummy arguments for statement functions
    REAL(wp), SAVE                       ::  znimax             ! Maximum ice concentration
    ! This is constant in both cloudice and graupel
    LOGICAL, PARAMETER                   ::  lsuper_coolw = .TRUE.

    REAL(wp) :: zdebug

    CHARACTER(len=250) :: message_text

    LOGICAL :: lzacc             ! OpenACC flag

!------------------------------------------------------------------------------

    ! Number of activate ice crystals;  ztx is temperature. Statement functions
    fxna       (ztx,ztmelt) = 1.0E2_wp  * EXP(0.2_wp   * (ztmelt - ztx))  ! 1/m^3
    fxna_cooper(ztx,ztmelt) = 5.0E+0_wp * EXP(0.304_wp * (ztmelt - ztx))  ! 1/m^3

!------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (lmessage_light .OR. lmessage_full) THEN
      message_text(:) = ' '
      WRITE(message_text,'(a,i0)') 'Computing dbz3d_lin for inw_gscp = ', igscp
      CALL message(TRIM(routine), TRIM(message_text), all_print=.TRUE.)
    ENDIF
      
    lqnc_input = PRESENT(n_cloud_s)

    z_fac_ice_dry = (rho_w/rho_ice)**2 * K_ice/K_w
    z_fac_ice_wet = 1.0_wp
    
    IF (firstcall) THEN

      mom_fac = (6.0_wp / (pi * rho_w))**2
      
      ! Parameters for cloud droplets:
      !   PSD consistent to autoconversion parameterization of SB2001:
      !   gamma distribution w.r.t. mass with gamma shape parameter zcnue from gscp_data.f90
      D_c_fix = 2.0E-5_wp                         ! constant mean mass diameter of cloud droplets (if n_cloud_s(:,:)
      x_c_fix = pi * rho_w / 6.0_wp * D_c_fix**3  !   is not present in argument list)
      cloud_tmp%nu = zcnue     ! gamma shape parameter mu for cloud droplets, in Axel's notation this is called nu instead!
      cloud_tmp%mu = 1.0_wp    ! 2nd shape parameter of generalized gamma distribution 
      z_fac_c = moment_gamma(cloud_tmp,2) * mom_fac

      ! Parameters for cloud ice:
      ami = zami              ! mass-size-relation prefactor
      bmi = 3.0_wp            ! mass-size-relation exponent
      IF( lsuper_coolw) THEN
        znimax = znimax_Thom        ! znimax_Thom = 250.E+3_wp 1/m^3
      ELSE
        znimax = fxna(zthn,T_melt)  ! Maximum number of cloud ice crystals in 1/m^3
      END IF

      ! Parameters for rain:
      mue_rain_c = mu_rain
      nor = 8.0e6_wp * EXP(3.2_wp*mue_rain_c) * (0.01_wp)**(-mue_rain_c)
      p_r = (7.0_wp+mue_rain_c) / (4.0_wp+mue_rain_c)
      z_r = nor*GAMMA(7.0_wp+mue_rain_c) * (pi*rho_w*nor*GAMMA(4.0_wp+mue_rain_c)/6.0_wp)**(-p_r)
      
      ! Parameters for snow and graupel:
      IF (igscp == 1 .or. igscp == 3) THEN
        
        ams = zams_ci
        bms = zbms
        p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
        z_s = mom_fac*ams**2 * GAMMA(2.0_wp*bms+1.0_wp) * (ams*GAMMA(bms+1.0_wp))**(-p_s)

        IF (lmessage_light) THEN
          WRITE (*, *) TRIM(routine)//": cloud ice scheme (using rain and snow)"
          WRITE (*,'(A,F10.3)') '     p_r = ', p_r
          WRITE (*,'(A,F10.3)') '     z_r = ', z_r
          WRITE (*,'(A,F10.3)') '     p_s = ', p_s
          WRITE (*,'(A,F10.3)') '     z_s = ', z_s
        ENDIF
      
        firstcall = .FALSE.

      ELSEIF (igscp == 2) THEN

        IF ( .NOT. PRESENT(q_graupel) ) THEN
          message_text(:) = ' '
          WRITE(message_text, '(a,i0,a)') 'optional argument q_graupel needed for inwp_gscp = ',igscp,' but is not present!'
          CALL finish(TRIM(routine), TRIM(message_text))
        END IF

        ams = zams_gr         
        bms = zbms        
        amg = 169.6_wp
        bmg = 3.1_wp
        nog = 4.E6_wp
        p_s = (2.0_wp*bms+1.0_wp)/(bms+1.0_wp)
        p_g = (2.0_wp*bmg+1.0_wp)/(bmg+1.0_wp)
        z_s = mom_fac*ams**2 * GAMMA(2.0_wp*bms+1.0_wp) *       &
                         (ams*GAMMA(bms+1.0_wp))**(-p_s)
        z_g = mom_fac*amg**2 * nog*GAMMA(2.0_wp*bmg+1.0_wp) *       &
                         (amg*nog*GAMMA(bmg+1.0_wp))**(-p_g)
        
        IF (lmessage_light) THEN
          WRITE (*, *) TRIM(routine)//": graupel scheme (using rain, snow, graupel)"
          WRITE (*,'(A,F10.3)') '     p_r = ', p_r
          WRITE (*,'(A,F10.3)') '     z_r = ', z_r
          WRITE (*,'(A,F10.3)') '     p_s = ', p_s
          WRITE (*,'(A,F10.3)') '     z_s = ', z_s
          WRITE (*,'(A,F10.3)') '     p_g = ', p_g
          WRITE (*,'(A,F10.3)') '     z_g = ', z_g
        ENDIF
        firstcall = .FALSE.
        
      ELSE
        
        message_text(:) = ' '
        WRITE(message_text, '(a,i0,a)') 'inwp_gscp = ',igscp,' not implemented for DBZ calculation!'
        CALL finish(TRIM(routine), trim(message_text))

      ENDIF
      
    ENDIF

    !$ACC DATA COPYIN(n_cloud_s) PRESENT(q_cloud, q_graupel, q_ice, q_rain, q_snow, rho, T, z_radar) IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,ztc,zn0s,nn,hlp,alf,bet,m2s,m3s, &
!$OMP            rho_c,rho_i,rho_r,rho_s,rho_g,z_cloud,z_rain,z_snow,z_ice,z_grau,x_c,n_i,x_i_mono), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = startblk, endblk

      IF (jb == startblk) THEN
        i_startidx = startidx1
      ELSE
        i_startidx = 1
      ENDIF

      IF (jb == endblk) THEN
        i_endidx = endidx2
      ELSE
        i_endidx = npr
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2) &
      !$ACC   PRIVATE(alf, bet, hlp, m2s, m3s, nn, n_i, rho_c, rho_i, rho_r, rho_s, x_c, x_i_mono, zn0s, ztc) &
      !$ACC   PRIVATE(z_cloud, z_ice, z_rain, z_snow)
      DO jk = jk_start, nlev
        DO jc = i_startidx, i_endidx

          ! convert mixing ratios into densities
          rho_c = q_cloud(jc,jk,jb) * rho(jc,jk,jb)
          rho_i = q_ice(jc,jk,jb)   * rho(jc,jk,jb)
          rho_r = q_rain(jc,jk,jb)  * rho(jc,jk,jb)
          rho_s = q_snow(jc,jk,jb)  * rho(jc,jk,jb)

          ! .. cloud droplets (gamma distribution w.r.t. mass x and mu_x=zcnue from gscp_data.f90):
          IF (rho_c >= q_crit_radar) THEN
            IF (lqnc_input) THEN
              x_c = q_cloud(jc,jk,jb) / MAX(n_cloud_s(jc,jb), n_cloud_min)
            ELSE
              x_c = x_c_fix
            END IF
            z_cloud = z_fac_c * rho_c * x_c
            z_radar(jc,jk,jb) = z_cloud * convfac
          ELSE
            z_radar(jc,jk,jb) = 0._wp
          END IF
          
          ! .. rain:
          IF (rho_r >= q_crit_radar) THEN
            z_rain  = z_r * EXP(p_r * LOG(rho_r))
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_rain * convfac
          ENDIF

          ! .. cloud ice (monodisperse size distribution):
          IF (rho_i >= q_crit_radar) THEN
            IF( lsuper_coolw) THEN
              n_i   = MIN( fxna_cooper(T(jc,jk,jb),T_melt), znimax )  ! 1/m^3
            ELSE
              n_i   = MIN( fxna(T(jc,jk,jb),T_melt), znimax )         ! 1/m^3
            END IF
            x_i_mono = rho_i / n_i
            z_ice = MERGE(z_fac_ice_dry, z_fac_ice_wet, T(jc,jk,jb) < T_melt) * mom_fac * rho_i * x_i_mono
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_ice * convfac
          END IF

          ! .. snow:
          IF (rho_s >= q_crit_radar) THEN

            IF (isnow_n0temp == 1) THEN
              ! Calculate n0s using the temperature-dependent
              ! formula of Field et al. (2005)
              ztc = T(jc,jk,jb) - T_melt
              ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
              zn0s = zn0s1*EXP(zn0s2*ztc)
              zn0s = MIN(zn0s,1e9_wp)
              zn0s = MAX(zn0s,1e6_wp)
            ELSEIF (isnow_n0temp == 2) THEN
              ! Calculate n0s using the temperature-dependent moment
              ! relations of Field et al. (2005)
              ztc = T(jc,jk,jb) - T_melt
              ztc = MAX(MIN(ztc,0.0_wp),-40.0_wp)
              nn  = 3
              hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
                  + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
              alf = EXP(hlp * LOG(10.0_wp))
              bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
                  + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
!!$ UB: caution! Here is the exponent bms=2.0 hardwired! not ideal!
              m2s = rho_s / ams
              m3s = alf*EXP(bet * LOG(m2s))
              hlp  = zn0s1*EXP(zn0s2*ztc)
!!$ UB: the 13.5 is actually 3^3 / gamma(3) ...
              zn0s = 13.50_wp * m2s*(m2s/m3s)**3
              zn0s = MAX(zn0s,0.5_wp*hlp)
              zn0s = MIN(zn0s,1e2_wp*hlp)
              zn0s = MIN(zn0s,1e9_wp)
              zn0s = MAX(zn0s,8e5_wp)
            ELSE
              ! Old constant n0s
              zn0s = 8.0e5_wp
            ENDIF
            hlp = z_s * EXP((1.0_wp-p_s) * LOG(zn0s))

            z_snow = hlp * EXP(p_s*LOG(rho_s)) * MERGE(z_fac_ice_dry, z_fac_ice_wet, T(jc,jk,jb) < T_melt)
            
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_snow * convfac

          END IF

        ENDDO
      ENDDO

      IF (igscp == 2) THEN
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2) PRIVATE(rho_g, z_grau)
        DO jk = jk_start, nlev
          DO jc = i_startidx, i_endidx
            rho_g = q_graupel(jc,jk,jb) * rho(jc,jk,jb)
            IF (rho_g >= q_crit_radar) THEN
              z_grau = z_g * EXP(p_g * LOG(rho_g)) * MERGE(z_fac_ice_dry, z_fac_ice_wet, T(jc,jk,jb) < T_melt)
              z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_grau * convfac
            ENDIF
          ENDDO
        ENDDO

      END IF
      !$ACC END PARALLEL

    ENDDO
!$OMP END DO
!$OMP END PARALLEL

    IF ( lmessage_light .OR. lmessage_full ) THEN
      !$ACC UPDATE HOST(z_radar) ASYNC(1) IF(lzacc)
      !$ACC WAIT(1)
      zdebug = MAXVAL(z10olog10 * LOG(z_radar + eps))
      message_text(:) = ' '
      WRITE (message_text, '(a,i4,a,1x,f6.1)') 'reflectivity statistics on proc ',my_id_for_message, &
           ': MAX dBZ tot : ', zdebug
      CALL message (TRIM(routine), TRIM(message_text), all_print=.TRUE.)
    ENDIF

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE compute_field_dbz_1mom


  !>
  !! Calculate radar reflectivity for the 2-moment microphysics scheme in linear units mm^6/m^3
  !!
  SUBROUTINE compute_field_dbz_2mom( npr, nlev, nblks, startblk, endblk, jk_start,          &
                                     startidx1, endidx2,                                    &
                                     lmessage_light, lmessage_full, my_id_for_message,      &
                                     rho_w, rho_ice,                                        &
                                     K_w, K_ice, T_melt, q_crit_radar,                      &
                                     luse_mu_Dm_rain, T, rho,                               &
                                     q_cloud, q_rain, q_ice, q_snow, q_graupel, q_hail,     &
                                     n_cloud, n_rain, n_ice, n_snow, n_graupel, n_hail,     &
                                     ql_graupel, ql_hail, z_radar, lacc )

    !------------------------------------------------------------------------------
    !
    ! Description:  Calculation of grid point values for effective radar
    !               reflectivity factor Z in dBZ.
    !
    ! Method:       Rayleigh-Approximation for the Back-Scattering
    !               (no attenuation!), Debye-Approximation for the
    !               effective refractive index of two-component ice-air-
    !               mixture particles (dry ice, snow, and graupel).
    !               Melting particles by substitution of the ice substance
    !               with water when T_air > 273.16 K. For the applied
    !               Rayleigh scattering, this is equal to assuming instantaneous
    !               melting, because only the square of the total water mass of the particle
    !               counts ( = the total square of its dipole moment).
    !               
    ! Inputs:       npr, nlev, nblks : field dimensions
    !               startblk, endblk, jk_start, startidx1, endidx2   : loop start and end indices
    !               rho_w        : bulk density of pure water [kg/m**3]
    !               rho_ice      : bulk density of pure ice   [kg/m**3]
    !               K_w          : dielectric constant of water
    !               K_ice        : dielectric constant of ice
    !               T_melt       : melting temperature of ice
    !               q_crit_radar : threshold for the q's to compute reflectivity  [kg/m**3]
    !               luse_mu_Dm_rain : switch to enable usage of mu-D-relation for rain outside cloud cores (Seifert, 2008)
    !               T            : temperature field          [K]
    !               rho          : air density                [kg/m**3]
    !               q_cloud      : cloud water mixing ratio   [kg/kg] 
    !               q_rain       : rain water mixing ratio    [kg/kg] 
    !               q_ice        : cloud ice mixing ratio     [kg/kg]
    !               q_snow       : snow mixing ratio          [kg/kg]  
    !               q_graupel    : graupel mixing ratio       [kg/kg]
    !               q_hail       : hail mixing ratio          [kg/kg]
    !               n_cloud      : cloud water number density [1/kg] 
    !               n_rain       : rain water number density  [1/kg] 
    !               n_ice        : cloud ice number density   [1/kg] 
    !               n_snow       : snow number density        [1/kg] 
    !               n_graupel    : graupel number density     [1/kg]
    !               n_hail       : hail number density        [1/kg]
    !   OPTIONAL:   ql_graupel   : liquid water on graupel mixing ratio  [kg/kg]
    !   OPTIONAL:   ql_hail      : liquid water on hail mixing ratio     [kg/kg]
    !
    ! Output:       z_radar      : 3D field of Z              [mm^6/m^3]
    ! 
    !
    !------------------------------------------------------------------------------

    ! Input/Output parameters:
    !-------------------------

    INTEGER,  INTENT(IN) :: npr, nlev, nblks
    INTEGER,  INTENT(IN) :: startblk, endblk, jk_start, startidx1, endidx2, my_id_for_message
    LOGICAL,  INTENT(in) :: lmessage_light, lmessage_full
    REAL(wp), INTENT(in) :: K_w, K_ice, T_melt, rho_w, rho_ice, q_crit_radar
    LOGICAL,  INTENT(in) :: luse_mu_Dm_rain

    REAL(wp), INTENT(IN) :: T(:,:,:),          &
                            rho(:,:,:),        &
                            q_cloud(:,:,:),    &
                            q_rain(:,:,:),     &
                            q_ice(:,:,:),      &
                            q_snow(:,:,:),     &
                            q_graupel(:,:,:),  &
                            q_hail(:,:,:),     &
                            n_cloud(:,:,:),    &
                            n_rain(:,:,:),     &
                            n_ice(:,:,:),      &
                            n_snow(:,:,:),     &
                            n_graupel(:,:,:),  &
                            n_hail(:,:,:)

    REAL(wp), INTENT(IN), OPTIONAL ::          &
                            ql_graupel(:,:,:), &
                            ql_hail(:,:,:)
    
    REAL(wp), INTENT(OUT) :: z_radar(:,:,:)

    LOGICAL, OPTIONAL, INTENT(IN)  :: lacc              !< initialization flag

    ! Local Variables
    !----------------

    CHARACTER(len=*), PARAMETER :: routine = modname//': compute_field_dbz_2mom'

    CHARACTER(len=250) :: message_text

    REAL(wp), PARAMETER :: eps  = 1.00E-15_wp
    REAL(wp), PARAMETER :: eps2 = 1.00E-20_wp
    REAL(wp), PARAMETER :: convfac = 1.E18_wp

    INTEGER       :: jc,jk, jb, pe_center, i_startidx, i_endidx
    LOGICAL, SAVE :: firstcall = .TRUE.
    LOGICAL       :: llwf_scheme

    REAL(wp)      :: T_a,                            &
                     q_c, n_c, x_c,                  &
                     q_r, n_r, x_r,                  &
                     q_g, n_g, x_g,                  &
                     q_h, n_h, x_h,                  &
                     q_s, n_s, x_s,                  &
                     q_i, n_i, x_i,                  &
                     z_fac_ice_dry, z_fac_ice_wet,   &
                     d_r, muD, z_fac_r_muD

    REAL(wp), SAVE :: z_fac_c, z_fac_r, z_fac_i, z_fac_s, z_fac_g, z_fac_h, mom_fac

    TYPE(particle)        :: cloud, rain
    TYPE(particle_frozen) :: ice, snow, graupel, hail

    LOGICAL :: lzacc             ! OpenACC flag
    CALL set_acc_host_or_device(lzacc, lacc)

!!$ LWF scheme not yet implemented:    CLASS(particle_lwf)    :: graupel_lwf, hail_lwf

    IF (lmessage_light .OR. lmessage_full) THEN
      message_text(:) = ' '
      WRITE(message_text,'(a,i0)') 'Computing dbz3d_lin for inw_gscp = 4, 5, 6 or 7'
      CALL message(TRIM(routine), TRIM(message_text))
    ENDIF

    IF (PRESENT(ql_graupel) .AND. PRESENT(ql_hail)) THEN
      llwf_scheme = .TRUE.
    ELSE
      llwf_scheme = .FALSE.
    END IF

    CALL init_2mom_scheme(cloud, rain, ice, snow, graupel, hail)

    ! .. K_ice and K_w  might change from call to call, so not in firstcall only:
    z_fac_ice_dry = (rho_w/rho_ice)**2 * K_ice/K_w
    z_fac_ice_wet = 1.0_wp

    IF (firstcall) THEN
      mom_fac = (6.0_wp / (pi * rho_w))**2
      z_fac_c = moment_gamma(cloud,2)   * mom_fac
      z_fac_r = moment_gamma(rain,2)    * mom_fac
      z_fac_i = moment_gamma(ice,2)     * mom_fac
      z_fac_s = moment_gamma(snow,2)    * mom_fac
      z_fac_g = moment_gamma(graupel,2) * mom_fac
      z_fac_h = moment_gamma(hail,2)    * mom_fac
      IF (lmessage_light) THEN
        WRITE (*, *) TRIM(routine)//":"
        WRITE (*,'(A,D10.3)') "     z_fac_c = ",z_fac_c
        WRITE (*,'(A,D10.3)') "     z_fac_r = ",z_fac_r
        WRITE (*,'(A,D10.3)') "     z_fac_i = ",z_fac_i
        WRITE (*,'(A,D10.3)') "     z_fac_s = ",z_fac_s
        WRITE (*,'(A,D10.3)') "     z_fac_g = ",z_fac_g
        WRITE (*,'(A,D10.3)') "     z_fac_h = ",z_fac_h
        WRITE (*,'(A,D10.3)') "     fak_ice = ",z_fac_ice_dry
      ENDIF
      firstcall = .FALSE.
    ENDIF

    !$ACC DATA COPYIN(cloud, graupel, hail, ice, rain, snow) &
    !$ACC   PRESENT(n_cloud, n_graupel, n_hail, n_ice, n_rain, n_snow, ql_graupel, ql_hail) &
    !$ACC   PRESENT(q_cloud, q_graupel, q_hail, q_ice, q_rain, q_snow, rho, T, z_radar) &
    !$ACC   IF(lzacc)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,i_startidx,i_endidx,T_a,q_c,n_c,q_r,n_r,q_i,n_i,q_s,n_s,q_g,n_g,q_h,n_h, &
!$OMP&           x_c,x_r,x_i,x_s,x_g,x_h,d_r,muD,z_fac_r_muD), ICON_OMP_RUNTIME_SCHEDULE
    DO jb = startblk, endblk

      IF (jb == startblk) THEN
        i_startidx = startidx1
      ELSE
        i_startidx = 1
      ENDIF

      IF (jb == endblk) THEN
        i_endidx = endidx2
      ELSE
        i_endidx = npr
      ENDIF

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      !$ACC LOOP GANG VECTOR COLLAPSE(2) &
      !$ACC   PRIVATE(d_r, muD, n_c, n_g, n_h, n_i, n_r, n_s, q_c, q_g, q_h, q_i, q_r, q_s, T_a) &
      !$ACC   PRIVATE(x_c, x_g, x_h, x_i, x_r, x_s, z_fac_r_muD)
      DO jk = jk_start, nlev
        DO jc = i_startidx, i_endidx

          z_radar(jc,jk,jb) = 0.0_wp

          T_a = T(jc,jk,jb)
          q_c = MAX(q_cloud(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
          n_c = MAX(n_cloud(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
          q_r = MAX(q_rain(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          n_r = MAX(n_rain(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          q_i = MAX(q_ice(jc,jk,jb)  , 0.0_wp) * rho(jc,jk,jb)
          n_i = MAX(n_ice(jc,jk,jb)  , 0.0_wp) * rho(jc,jk,jb)
          q_s = MAX(q_snow(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          n_s = MAX(n_snow(jc,jk,jb) , 0.0_wp) * rho(jc,jk,jb)
          IF (llwf_scheme) THEN
            q_g = MAX(q_graupel(jc,jk,jb) + ql_graupel(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
            q_h = MAX(q_hail(jc,jk,jb)    + ql_hail(jc,jk,jb)   , 0.0_wp) * rho(jc,jk,jb)
          ELSE
            q_g = MAX(q_graupel(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
            q_h = MAX(q_hail(jc,jk,jb)   , 0.0_wp) * rho(jc,jk,jb)
          END IF
          n_g = MAX(n_graupel(jc,jk,jb), 0.0_wp) * rho(jc,jk,jb)
          n_h = MAX(n_hail(jc,jk,jb)   , 0.0_wp) * rho(jc,jk,jb)

          x_c = MIN( MAX(q_c/(n_c+eps2),cloud%x_min),cloud%x_max )
          x_r = MIN( MAX(q_r/(n_r+eps2),rain%x_min),rain%x_max )
          x_i = MIN( MAX(q_i/(n_i+eps2),ice%x_min),ice%x_max )
          x_s = MIN( MAX(q_s/(n_s+eps2),snow%x_min),snow%x_max )
          x_g = MIN( MAX(q_g/(n_g+eps2),graupel%x_min),graupel%x_max )
          x_h = MIN( MAX(q_h/(n_h+eps2),hail%x_min),hail%x_max )

          ! .. Cloud water reflectivity:
          IF (q_c < q_crit_radar) x_c = 0.0_wp
          z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_fac_c * q_c * x_c

          ! .. Rain water reflectivity:
          IF (q_r >= q_crit_radar) THEN
            IF (q_c > q_crit_radar .OR. .NOT.luse_mu_Dm_rain) THEN
              ! Inside of cloud cores assume generalized gamma DSD:
              z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_fac_r * q_r * x_r
            ELSE
              ! Outside of cloud cores assume mu-D-relation Seifert (2008):
              d_r = rain%a_geo * EXP(rain%b_geo * LOG(x_r))
              muD = rain_mue_dm_relation(rain_coeffs,d_r)
              z_fac_r_muD = mom_fac*(muD+6.0)*(muD+5.0)*(muD+4.0)/((muD+3.0)*(muD+2.0)*(muD+1.0))
              z_radar(jc,jk,jb) = z_radar(jc,jk,jb) + z_fac_r_muD * q_r * x_r
            END IF
          END IF

          ! .. Ice species reflectivity:
          IF (q_i < q_crit_radar) x_i = 0.0_wp
          IF (q_s < q_crit_radar) x_s = 0.0_wp
          IF (q_g < q_crit_radar) x_g = 0.0_wp
          IF (q_h < q_crit_radar) x_h = 0.0_wp

          IF (T_a < T_melt) THEN
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_i * q_i * x_i * z_fac_ice_dry
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_s * q_s * x_s * z_fac_ice_dry
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_g * q_g * x_g * z_fac_ice_dry
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_h * q_h * x_h * z_fac_ice_dry
          ELSE
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_i * q_i * x_i * z_fac_ice_wet
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_s * q_s * x_s * z_fac_ice_wet
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_g * q_g * x_g * z_fac_ice_wet
            z_radar(jc,jk,jb) = z_radar(jc,jk,jb)                            &
                 &       + z_fac_h * q_h * x_h * z_fac_ice_wet
          ENDIF

          ! conversion of output unit
          z_radar(jc,jk,jb) = z_radar(jc,jk,jb) * convfac

        END DO
      END DO
      !$ACC END PARALLEL

    END DO
!$OMP END DO
!$OMP END PARALLEL

    IF (lmessage_light .OR. lmessage_full ) THEN
      !$ACC UPDATE HOST(z_radar) ASYNC(1) IF(lzacc)
      !$ACC WAIT(1)
      message_text(:) = ' '
      WRITE (message_text, '(A,i4,2(A,F10.1))') 'on proc ',my_id_for_message,': '// &
           'MAX dBZ = ', &
           MAXVAL(10.0_wp / LOG(10.0_wp) * LOG(z_radar + eps)), &
           '  MIN dBZ = ', &
           MINVAL(10.0_wp / LOG(10.0_wp) * LOG(z_radar + eps))
      CALL message(TRIM(routine), TRIM(message_text), all_print=.TRUE.)
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA
  
  END SUBROUTINE compute_field_dbz_2mom

END MODULE mo_opt_nwp_reflectivity

