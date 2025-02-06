!
! mo_art_nucleation_warm
! This routine calculates the number of activated cloud droplets
! Original code by A. Nenes et al.
! Rewritten for the use in ICON-ART and COSMO-ART by D. Rieger (KIT)
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

MODULE mo_art_nucleation_warm
! ICON
  USE mo_kind,                          ONLY: wp
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: activation_master
  
  REAL(wp)                    :: &
    &    ama    = 29e-3_wp,      &       !< Air molecular weight
    &    r_gas  = 8.31_wp,       &       !< Universal gas constant
    &    rho_w  = 1000._wp,      &       !< Water density
    &    dhv    = 2.25e6_wp,     &       !< Water enthalpy of vaporization
    &    cp_air = 1.0061e3_wp,   &       !< Specific heat capacity of air
    &    accom  = 0.06_wp,       &       !< Default accommodation coef
    &    amw    = 18.e-3_wp,     &       !< Water molecular weight
    &    pi     = 3.1415927e0_wp,&       !< PI
    &    dw     = 2.75e-10_wp,   &       !< Water Molecule Diameter
    &    a_fhh  = 2.50_wp,       &       !< FHH parameter
    &    b_fhh  = 1.20_wp,       &       !< FHH parameter
    &    grav   = 9.81_wp
    
    
  REAL(wp)                    :: &
    &    th_cond_air,            &       !< Air thermal conductivity
    &    rho_air,                &       !< Air density
    &    dv,                     &       !< Water vapor diffusivity in air
    &    sfc_tens,               &       !< Surface Tension for water (J m-2)
    &    p_sat,                  &       !< Saturated water pressure
    &    akoh,                   &       !< Curvature param
    &    alfa,                   &       !< 
    &    bet1,                   &       !< 
    &    bet2                            !< 
    
  REAL(wp),ALLOCATABLE  ::       &
    &    dpc(:),                 &    !< (dim: nmodes) 
    &    sg(:),                  &    !< Critical Super Saturation (dim: nmodes)
    &    wgs(:),                 &    !< gaussian weights (dim: Npgauss)
    &    xgs(:)                       !< gaussian weights (dim: Npgauss)
    
  INTEGER               ::          &
    &    npgauss=10                      !<
    
  LOGICAL               ::          &
    &    l_do_scav                       !< if m_scav,n_scav and mas_conc are passed
                                         !< those are calculated and passed as optional
                                         !< arguments
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE activation_master(nmodes,l_koehler,dpg,amfs,vhf,ams,   &  !< INTENT(IN)
          &                  rho_sol,rho_insol,temp,pres,w_parc,  &  !< INTENT(IN) 
          &                  sig_w,nmb_conc,sig_g,                &  !< INTENT(IN)
          &                  n_act,n_act_fhh,s_max,               &  !< INTENT(OUT)
          &                  mas_conc,n_scav,m_scav)                 !< INTENT(IN,OUT,OUT) 
                                                                     !  FOR WET-PHASE-CHEM
  
  !=================================================================
  ! THIS SUBROUTINE ACTS AS INTERFACE FOR THE ACTIVATION.
  ! ALL THE SINGLE ACTIVATION ROUTINES ARE CALLED WITHIN THIS ROUTINE
  ! ---
  ! Written by DANIEL RIEGER (KIT)
  !=================================================================
  
  INTEGER, INTENT(IN)     ::    &
     &      nmodes                   !< Number of modes
           
  LOGICAL, INTENT(IN)     ::    &
     &      l_koehler(:)             !< If true, use koehler activation, else use fhh (dim: nmodes)
  
  REAL(wp), INTENT(IN)    ::    &
     &      dpg(:),             &    !< Median diameter of mode (dim: nmodes)
     &      amfs(:),            &    !< Soluble mass fraction (dim: nmodes)
     &      vhf(:),             &    !< Van't Hoff factor for soluble frac.(ions molec-1) 
                                     !    (dim: nmodes)
     &      ams(:),             &    !< Molar mass of Soluble fraction (kg mol-1) (dim: nmodes)
     &      rho_sol(:),         &    !< Density of Soluble fraction (kg m-3) (dim: nmodes)
     &      rho_insol(:),       &    !< Density of Insoluble fraction (kg m-3) (dim: nmodes)
     &      temp,               &    !< Parcel Temperature
     &      pres,               &    !< Parcel Pressure
     &      w_parc,             &    !< parcel updraft velocity 
     &      sig_w,              &    !< sigma
     &      nmb_conc(:),        &    !< total number concentration of particles
     &      sig_g(:)                 !< geometric dispersion
     
  REAL(wp), INTENT(OUT)   ::    &
     &      n_act,              &    !< Total number of activated particles
     &      n_act_fhh,          &    !< Number of activated particles (FHH)
     &      s_max                    !< maximum supersaturation
     
  REAL(wp), INTENT(OUT),OPTIONAL :: &
     &      n_scav(:), m_scav(:)     !< number and mass conc. of scavenged particles
     
  REAL(wp), INTENT(IN), OPTIONAL :: &
     &      mas_conc(:)              !< total mass concentration of particles
  
     
  !-------------------------------------
  ! Section 1.0
  !   Allocations and initializations
  !-------------------------------------
  
  CALL util_memory_activation('alloc',nmodes)
  
  IF(PRESENT(n_scav) .AND. PRESENT(m_scav) .AND. PRESENT(mas_conc)) THEN
    l_do_scav = .TRUE.
  ELSE 
    l_do_scav = .FALSE.
  ENDIF
  
  !-------------------------------------
  ! Section 2.0
  !   Points and weights for n-point Gauss-Integration
  !-------------------------------------
  
  CALL gauleg
  
  !-------------------------------------
  ! Section 3.0
  !   Calculate the CCN spectrum of the aerosol
  !-------------------------------------
  
  CALL ccn_spectrum(dpg,amfs,vhf,ams,rho_sol,rho_insol,l_koehler, &
    &               temp,pres,nmodes)
  
  !-------------------------------------
  ! Section 4.0
  !   Calculate number of activated particles and maximum supersaturation
  !-------------------------------------
  
  IF(l_do_scav) THEN
    CALL pdf_activation(w_parc,sig_w,n_act,n_act_fhh,s_max,temp,pres, &
      &                 l_koehler,nmodes,nmb_conc,sig_g,mas_conc,n_scav,m_scav)
  ELSE 
    CALL pdf_activation(w_parc,sig_w,n_act,n_act_fhh,s_max,temp,pres, &
      &                 l_koehler,nmodes,nmb_conc,sig_g)
  ENDIF      
  !-------------------------------------
  ! Section 5.0
  !   Deallocation
  !-------------------------------------
             
  CALL util_memory_activation('dealloc')
  
END SUBROUTINE activation_master
!
! ----------------------------------------------------------------------
!
SUBROUTINE ccn_spectrum(dpg,amfs,vhf,ams,rho_sol,rho_insol,l_koehler,temp,pres,nmodes)

  !=================================================================
  ! THIS SUBROUTINE CALCULATES THE CCN SPECTRUM OF THE AEROSOL USING
  ! THE APPROPRIATE FORM OF KOHLER THEORY
  ! ORIGINALLY WRITTEN BY ATHANASIOS NENES FOR ONLY KOHLER PARTICLES
  ! MODIFIED BY PRASHANT KUMAR AND ATHANSIOS NENES TO INCLUDE 
  ! ACTIVATION BY FHH PARTICLES
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================
  
  INTEGER, INTENT(IN) ::      &
     &      nmodes               !< Number of modes
           
  LOGICAL, INTENT(IN)     ::  &
     &      l_koehler(:)         !< If true, use koehler activation, else use fhh (dim: nmodes)
                                 
  REAL(wp), INTENT(IN)    ::  &  
     &      dpg(:),           &  !< Median diameter of mode (dim: nmodes)
     &      amfs(:),          &  !< Soluble mass fraction (dim: nmodes)
     &      vhf(:),           &  !< Van't Hoff factor for soluble frac.(ions molec-1) (dim: nmodes)
     &      ams(:),           &  !< Molar mass of Soluble fraction (kg mol-1) (dim: nmodes)
     &      rho_sol(:),       &  !< Density of Soluble fraction (kg m-3) (dim: nmodes)
     &      rho_insol(:),     &  !< Density of Insoluble fraction (kg m-3) (dim: nmodes)
     &      temp,             &  !< Parcel Temperature
     &      pres                 !< Parcel Pressure
                                 
  REAL(wp)                ::  &  
     &      amfi,             &  !< Insoluble mass fraction
     &      rho_p,            &  !< Particle density
     &      par1,             &  !< Parameter
     &      par2,             &  !< Parameter
     &      dc,               &  !< local dpc(i)
     &      dg,               &  
     &      vol_fr_salt          !< Vol.Fr.Salt
                                 
  INTEGER                 ::  &  
     &      i                    !< counting variable
  
  !-------------------------------------
  ! Start Routine
  !-------------------------------------
  
  !-------------------------------------
  ! Calculate thermophysical properties
  !-------------------------------------
  CALL thermophysical_props(pres,temp)

  !-------------------------------------
  ! Calculate critical properties
  !-------------------------------------

  DO i=1,nmodes
    IF (l_koehler(i)) THEN ! Koehler Modes
      IF(amfs(i) < 1.0_wp) THEN
        amfi  = MAX( (1._wp-amfs(i)) , 0._wp )
        rho_p = amfs(i)*rho_sol(i) + amfi*rho_insol(i)
        vol_fr_salt  = amfs(i)/rho_sol(i)/(amfs(i)/rho_sol(i)+amfi/rho_insol(i))
      ELSE
        amfi  = 0.0_wp
        rho_p = rho_sol(i) 
        vol_fr_salt  = 1.0_wp
      ENDIF
      par1  = 4._wp*rho_w*ams(i)/27._wp/vhf(i)/rho_p/amw/dpg(i)**3
      par1  = par1/vol_fr_salt
      par2  = SQRT(par1*akoh**3)
      sg(i) = EXP(par2) - 1._wp
    ELSE                   ! FHH Modes  
      dg=dpg(i)    
      CALL dpcfhh(dg,temp,pres,dc)
      dpc(i)=dc
      ! Calculating Critical Super Saturation by Taylor Series Expansion 
      sg(i) = (akoh/dpc(i))+(-a_fhh*(((dpc(i)-dpg(i))/(2._wp*dw))**(-b_fhh)))
    ENDIF
  ENDDO   
  
  !-------------------------------------
  ! End Routine
  !-------------------------------------
  
END SUBROUTINE ccn_spectrum
!
! ----------------------------------------------------------------------
!
SUBROUTINE dpcfhh(Ddry,temp,pres,dc)
  
  !=================================================================
  ! SUBROUTINE DpcFHH
  ! THIS SUBROUTINE CALCULATES THE CRITICAL PARTICLE DIAMETER
  ! ACCORDING TO THE FHH ADSOSPRTION ISOTHERM THEORY.
  !
  ! WRITTEN BY PRASHANT KUMAR AND ATHANASIOS NENES
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================
  
  REAL(wp),INTENT(in)     ::    &
     &  ddry,                   &
     &  temp,                   &
     &  pres
     
  REAL(wp),INTENT(out)    ::    &
     &  dc                       
     
  LOGICAL                 ::    &
     &  l_iterate
     
  REAL(wp)                ::    &
     &  mu1,mu2,                &         !< 
     &  x1,x2l,x3l,             &         !<
     &  dpcm,dpcl,dpcu,         &         !<
     &  fdpcl,fdpcm,            &         !<
     &  x2m,x3m

  CALL thermophysical_props(pres,temp)
              
  mu1= (akoh*2._wp*dw) / ((a_fhh*b_fhh) * ( (2._wp*dw)**(b_fhh+1._wp) ) )
  mu2= 1._wp / mu1

  dpcl = 0._wp         !Lower Limit
  dpcu = 10.e-4_wp     !Upper Limit
  l_iterate = .TRUE.
  DO WHILE (l_iterate)

    x1    = mu2**(1._wp/(1._wp+b_fhh))
    x2l   = dpcl**(2._wp/(1._wp+b_fhh))
    x3l   = x1*x2l
    fdpcl = ((dpcl-x3l)/ddry)-1._wp

    dpcm = (dpcu+dpcl)/2._wp

    x2m= dpcm**(2._wp/(1._wp+b_fhh))
    x3m= x1*x2m
    fdpcm=((dpcm-x3m)/ddry)-1
    IF ((fdpcl*fdpcm) <= 0._wp) THEN
      IF (ABS(fdpcm) <= 10.e-8_wp) THEN
        l_iterate = .FALSE.
      ELSE
        dpcl = dpcl
        dpcu = dpcm
      ENDIF
    ELSE IF ((fdpcl*fdpcm) >= 0._wp) THEN
      IF (ABS(fdpcm) <= 10.e-8_wp) THEN
        l_iterate = .FALSE.
      ELSE
        dpcl = dpcm
        dpcu = dpcu
      ENDIF
    ENDIF
  ENDDO
  
  dc = dpcm

END SUBROUTINE dpcfhh
!
! ----------------------------------------------------------------------
!
SUBROUTINE thermophysical_props(pres,temp)
  
  !=================================================================
  ! THIS SUBROUTINE CALCULATES THE THERMOPHYSICAL PROPERTIES
  ! WRITTEN BY ATHANASIOS NENES
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================
  
  REAL(wp), INTENT(IN)    ::    &
     &      temp,               &      !< Parcel Temperature
     &      pres
     
  REAL(wp)                ::    &
     &      presa,              &      !< 
     &      dlow,               &      !< 
     &      dbig,               &      !<
     &      coef,               &      !< coefficient
     &      a_data(7)                  !< data block for calculation of saturated water pressure
     
  INTEGER                 ::    &
     &      i                          !< counting variable
  
  !-------------------------------------
  ! Start Routine
  !-------------------------------------
  
  a_data   = (/6.107799610e+0_wp, 4.436518521e-1_wp, 1.428945805e-2_wp,   &
     &       2.650648471e-4_wp, 3.031240396e-6_wp, 2.034080948e-8_wp,     &
     &       6.136820929e-11_wp/)
     
  presa    = pres/1.013e5_wp
  
  !-------------------------------------
  ! Air density
  !-------------------------------------
     
  rho_air  = pres*ama/r_gas/temp
  
  !-------------------------------------
  ! Air thermal conductivity
  !-------------------------------------

  th_cond_air   = (4.39_wp+0.071_wp*temp)*1.e-3_wp

  !-------------------------------------
  ! Water vapor diffusivity in air dv
  !-------------------------------------
  
  dv    = (0.211_wp/(presa))*(temp/273._wp)**1.94_wp
  dv    = dv*1.e-4_wp
  dbig  = 5.0e-6_wp
  dlow  = 0.207683_wp*((accom)**(-0.33048_wp))
  dlow  = dlow*1e-6_wp 
  coef  = SQRT(2._wp*pi*amw/(r_gas*temp))

  ! Average
  dv    = (dv/(dbig-dlow))*((dbig-dlow)-(2._wp*dv/accom)*coef         &    ! Non-continuum effects
   &    * (DLOG((dbig+(2._wp*dv/accom)*coef)/(dlow+(2._wp*dv/accom)   &
   &    * coef))))

  !-------------------------------------
  ! Saturated water pressure as fct(T)
  ! Valid from -50 C to 50 C
  !-------------------------------------

  ! Calculate polynomial (without exponentiation).
  p_sat = a_data(7)*(temp-273._wp)
  DO i=6,2,-1
    p_sat = (p_sat + a_data(i))*(temp-273._wp)
  ENDDO
  p_sat = (p_sat + a_data(1))*(1.e5_wp/1.e3_wp)

  !-------------------------------------
  ! Surface Tension for water
  !-------------------------------------
  
  sfc_tens   = 0.0761_wp - (1.55e-4_wp * (temp-273._wp) )
  
  !-------------------------------------
  ! Curvature param
  !-------------------------------------
  
  akoh = 4._wp*amw*sfc_tens/(r_gas*temp*rho_w)
  
  !-------------------------------------
  ! End Routine
  !-------------------------------------
  
END SUBROUTINE thermophysical_props
!
! ----------------------------------------------------------------------
!
SUBROUTINE util_memory_activation(yaction,nmodes)
  !=================================================================
  ! THIS SUBROUTINE HANDLES ALLOCATION AND DEALLOCATION ACCORDING
  ! TO THE NUMBER OF MODES. 
  ! ---
  ! Written by DANIEL RIEGER (KIT)
  !=================================================================
  
  INTEGER, INTENT(IN),OPTIONAL ::    &
     &      nmodes                        !< Number of modes
  
  CHARACTER(LEN=*),INTENT(IN)  ::    &
     &      yaction                       !< allocation or deallocation?
     
  IF (TRIM(yaction) == 'alloc') THEN
    IF(PRESENT(nmodes)) THEN
  !-------------------------------------
  ! ALLOCATION
  !-------------------------------------
      ALLOCATE(dpc(nmodes),sg(nmodes))
      ALLOCATE(xgs(npgauss),wgs(npgauss))
    ELSE
      print *,'ERROR: Action is: ',yaction,', but no number of modes is passed!'
    ENDIF
  ELSE IF (TRIM(yaction) == 'dealloc') THEN
  !-------------------------------------
  ! DEALLOCATION
  !-------------------------------------
      IF(ALLOCATED(sg))  DEALLOCATE(sg)
      IF(ALLOCATED(dpc)) DEALLOCATE(dpc)
      IF(ALLOCATED(wgs)) DEALLOCATE(wgs)
      IF(ALLOCATED(xgs)) DEALLOCATE(xgs)
  ELSE
    print *,'ERROR: Wrong action: ',yaction
  ENDIF
  
  END SUBROUTINE util_memory_activation
!
! ----------------------------------------------------------------------
!
SUBROUTINE pdf_activation(w_parc,sig_w,n_act,n_act_fhh,s_max,temp,pres, &
                   &        l_koehler,nmodes,nmb_conc,sig_g,mas_conc,n_scav,m_scav)
  
  !=================================================================
  ! THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
  ! TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
  ! MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
  ! and Nenes (2004). THIS ROUTINE CALCULATES FOR A PDF OF
  ! UPDRAFT VELOCITIES.
  !
  ! WRITTEN BY ATHANASIOS NENES
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================
  
  REAL(wp), INTENT(IN)    ::    &
     &      w_parc,             &    !< parcel updraft velocity 
     &      sig_w,              &    !< sigma
     &      temp,               &    !< Parcel Temperature
     &      pres,               &    !< Parcel Pressure
     &      nmb_conc(:),        &    !< total number concentration of particles
     &      sig_g(:)                 !< geometric dispersion
     
  REAL(wp), INTENT(OUT)   ::    &
     &      n_act,              &    !< 
     &      n_act_fhh,          &    !< 
     &      s_max                    !< 
  
  INTEGER, INTENT(IN) ::    &
     &      nmodes                   !< Number of modes
           
  LOGICAL, INTENT(IN)     ::    &
     &      l_koehler(:)             !< If true, use koehler activation, else use fhh (dim: nmodes)
     
  REAL(wp), INTENT(OUT),OPTIONAL :: &
     &      n_scav(:), m_scav(:)     !< number and mass conc. of scavenged particles
     
  REAL(wp), INTENT(IN), OPTIONAL :: &
     &      mas_conc(:)              !< total mass concentration of particles
     
  REAL(wp)                ::    &
     &      pdf,                &
     &      w_pdf,              &    !< pdf updraft velocity
     &      p_updr_lmt,         &    !< probability of high updraft limit
     &      p_updr,             &    !< probability of high updraft
     &      w_lmt_hi,           &    !< Upper updraft limit
     &      w_lmt_lo,           &    !< Low updraft limit
     &      n_act_new,          &    !< 
     &      n_act_fhh_new,      &    !< 
     &      s_max_new,          &    !<
     &      scl_fctr                 !< Scaling factor
    
  REAL(wp),allocatable    ::    &
     &      n_scav_new(:),      &    !<
     &      m_scav_new(:)
 
  REAL(wp)                ::    &
     &      sq2pi=2.5066282746_wp
     
  INTEGER                 ::    &
     &      i                        !< counting variable
    
  !-------------------------------------
  ! Case where updraft is very small
  !-------------------------------------

  IF (w_parc <= 1.e-6_wp) THEN
    s_max     = 0._wp
    n_act     = 0._wp
    n_act_fhh = 0._wp
  ENDIF
  
  !-------------------------------------
  ! Single updraft case
  !-------------------------------------
  
  IF (sig_w < 1.e-14_wp) THEN
    IF(l_do_scav) THEN
      CALL activation(w_parc,n_act,n_act_fhh,s_max,temp,pres, &
               &      l_koehler,nmodes,nmb_conc,sig_g,mas_conc,n_scav,m_scav)
    ELSE
      CALL activation(w_parc,n_act,n_act_fhh,s_max,temp,pres, &
               &      l_koehler,nmodes,nmb_conc,sig_g)
    ENDIF
    
  !-------------------------------------
  ! PDF of updrafts
  !-------------------------------------
  ELSE
  
    n_act     = 0.0_wp
    n_act_fhh = 0.0_wp
    s_max     = 0.0_wp
    IF(l_do_scav) THEN
      ALLOCATE(n_scav_new(nmodes),m_scav_new(nmodes))
      n_scav(:) = 0.0_wp
      m_scav(:) = 0.0_wp
    ENDIF
    p_updr_lmt = 1.e-3_wp
    p_updr     = SQRT(-2.0_wp*LOG(p_updr_lmt*sig_w*sq2pi))
    w_lmt_hi   = w_parc + sig_w*p_updr
    w_lmt_lo   = MAX(0.05_wp, (w_parc - sig_w*p_updr) )
    
    IF(w_lmt_lo>=w_lmt_hi) THEN
      PRINT*,'Warning: w_lmt_lo>=w_lmt_hi',w_lmt_lo,'>=',w_lmt_hi
      w_lmt_hi=w_lmt_lo+0.1_wp
    ENDIF
    
    scl_fctr  = 0.5_wp*(w_lmt_hi-w_lmt_lo)                             ! Scaling for updrafts
    
    DO i=1,npgauss
      w_pdf     = w_lmt_lo + scl_fctr*(1.0-xgs(i))                       ! Updraft
      
      IF(l_do_scav) THEN
        CALL activation(w_pdf,n_act_new,n_act_fhh_new,s_max_new,temp,pres,l_koehler, &
                 &      nmodes,nmb_conc,sig_g,mas_conc,n_scav_new,m_scav_new)         ! # of drops
      ELSE
        CALL activation(w_pdf,n_act_new,n_act_fhh_new,s_max_new,temp,pres, &
                 &      l_koehler,nmodes,nmb_conc,sig_g) ! # of drops
      ENDIF
      pdf       = (1.0/sq2pi/sig_w)*EXP(-0.5*((w_pdf-w_parc)/sig_w)**2)  ! Prob. of updrafts
      n_act     = n_act + wgs(i)*(pdf*n_act_new)                         ! Integral for drops
      n_act_fhh = n_act_fhh + wgs(i)*(pdf*n_act_fhh_new)                 ! Integral for fhh drops
      s_max     = s_max + wgs(i)*(pdf*s_max_new)                         ! Integral for Smax
      IF(l_do_scav) THEN
        n_scav(:) = n_scav(:) + wgs(i)*(pdf*n_scav_new(:))
        m_scav(:) = m_scav(:) + wgs(i)*(pdf*m_scav_new(:))
      ENDIF
      IF (pdf < p_updr_lmt) EXIT
    ENDDO
    
    n_act = n_act*scl_fctr                                             ! Scale Integrals
    n_act_fhh = n_act_fhh*scl_fctr
    s_max = s_max*scl_fctr
    IF (l_do_scav) THEN
      n_scav(:) = n_scav(:)*scl_fctr
      m_scav(:) = m_scav(:)*scl_fctr
    ENDIF
  ENDIF

  !-------------------------------------
  ! Deallocations
  !-------------------------------------

  IF(l_do_scav) DEALLOCATE(n_scav_new,m_scav_new)


END SUBROUTINE pdf_activation
!
! ----------------------------------------------------------------------
!
SUBROUTINE activation (w_parc,n_drpl,sinteg3,s_max,temp,pres, &
                 &       l_koehler,nmodes,nmb_conc,sig_g,mas_conc,n_scav,m_scav)
  
  !=================================================================
  ! THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
  ! TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
  ! MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
  ! and Nenes (in preparation).
  !
  ! WRITTEN BY ATHANASIOS NENES FOR KOHLER PARTICLES
  ! MODIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH 
  ! PARTICLES 
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================

  REAL(wp), INTENT(IN)    ::    &
     &      w_parc,             &    !< parcel updraft velocity
     &      temp,               &    !< Parcel Temperature
     &      pres,               &    !< Parcel Pressure
     &      nmb_conc(:),        &    !< total number concentration of particles
     &      sig_g(:)                 !< geometric dispersion
     
  REAL(wp), INTENT(OUT)   ::    &
     &      n_drpl,             &    !< number of droplets
     &      sinteg3,            &    !< 
     &      s_max                    !< 
     
  INTEGER, INTENT(IN) ::        &
     &      nmodes                   !< Number of modes
           
  LOGICAL, INTENT(IN)     ::    &
     &      l_koehler(:)             !< If true, use koehler activation, else use fhh (dim: nmodes)
     
  REAL(wp), INTENT(OUT),OPTIONAL :: &
     &      n_scav(:), m_scav(:)     !< number and mass conc. of scavenged particles
     
  REAL(wp), INTENT(IN), OPTIONAL :: &
     &      mas_conc(:)              !< total mass concentration of particles
     
  REAL(wp)                ::    &
     &      beta,               &    !< 
     &      cf1,                &    !< 
     &      cf2,                &    !< 
     &      x1,x2,x3,           &    !<
     &      y1,y3,              &    !<
     &      sinteg1,            &    !<
     &      sinteg2
     
  REAL(wp)                ::    &
     &      eps=1.e-6_wp             !< Convergence criterion
     
  INTEGER                 ::    &
     &      i                        !< counting variable
     
  INTEGER                 ::    &
     &      max_iterations=100

  !-------------------------------------
  ! Setup constants
  !-------------------------------------
  alfa = grav*amw*dhv/cp_air/r_gas/temp/temp - grav*ama/r_gas/temp
  bet1 = pres*ama/p_sat/amw + amw*dhv*dhv/cp_air/r_gas/temp/temp
  bet2 = r_gas*temp*rho_w/p_sat/dv/amw/4._wp                              &
    &    + dhv*rho_w/4._wp/th_cond_air/temp*(dhv*amw/r_gas/temp - 1._wp)
  beta = 0.5_wp*pi*bet1*rho_w/bet2/alfa/w_parc/rho_air
  
  !DR: DEBUG: On of these operations is susceptible to floating invalid errors
  cf1  = 1._wp/bet2

  cf1  = cf1/alfa

  cf1  = cf1/w_parc

  cf1  = SQRT(cf1)

  cf1  = 0.5_wp*cf1

  !DR: END DEBUG
! cf1  = 0.5_wp*(((1._wp/bet2)/(alfa*w_parc))**0.5_wp)
  cf2  = akoh/3._wp
  
  !-------------------------------------
  ! Initial values for besection
  !-------------------------------------
   
  x1   = 1.0e-5_wp   ! Min cloud supersaturation -> 0
  CALL sintegral (x1,n_drpl,sinteg1,sinteg2,sinteg3,w_parc, &
          &       l_koehler,nmodes,nmb_conc,sig_g)
  y1   = (sinteg1*cf1+sinteg2*cf2+sinteg3*cf1)*beta*x1 - 1._wp

  x2   = 1._wp       ! MAX cloud supersaturation = 100%
  CALL sintegral (x2,n_drpl,sinteg1,sinteg2,sinteg3,w_parc, &
          &       l_koehler,nmodes,nmb_conc,sig_g)
  !-------------------------------------
  ! Perform besection
  !-------------------------------------
  DO i=1,max_iterations
    x3   = 0.5*(x1+x2)

    CALL sintegral (x3,n_drpl,sinteg1,sinteg2,sinteg3,w_parc, &
            &       l_koehler,nmodes,nmb_conc,sig_g)
    y3 = (sinteg1*cf1+sinteg2*cf2+sinteg3*cf1)*beta*x3 - 1._wp

    IF (SIGN(1._wp,y1)*SIGN(1._wp,y3) <= 0._wp) THEN  ! (y1*y3 .LE. ZERO)
      x2 = x3
    ELSE
      y1 = y3
      x1 = x3
    ENDIF

    IF (ABS(x2-x1) <= eps*x1) EXIT
  ENDDO

  x3   = 0.5_wp*(x1+x2)
  
  IF(l_do_scav) THEN
    CALL sintegral (x3,n_drpl,sinteg1,sinteg2,sinteg3,w_parc, &
            &       l_koehler,nmodes,nmb_conc,sig_g,mas_conc,n_scav,m_scav)
  ELSE
    CALL sintegral (x3,n_drpl,sinteg1,sinteg2,sinteg3,w_parc, &
            &       l_koehler,nmodes,nmb_conc,sig_g)
  ENDIF
  y3   = (sinteg1*cf1+sinteg2*cf2+sinteg3*cf1)*beta*x3 - 1._wp
      
  s_max = x3
  
END SUBROUTINE activation
!
! ----------------------------------------------------------------------
!
SUBROUTINE sintegral (s_parc, sum_drpl, sinteg1, sinteg2, sinteg3,w_parc, &
              &       l_koehler,nmodes,nmb_conc,sig_g,mas_conc,n_scav,m_scav)

  !=================================================================
  ! THIS SUBROUTINE CALCULATES THE CONDENSATION INTEGRALS, ACCORDING
  ! TO THE POPULATION SPLITTING ALGORITHM OF Nenes and Seinfeld (2003)
  ! Modal formulation according to Fountoukis and Nenes (2004)
  !
  ! WRITTEN BY ATHANASIOS NENES for Kohler Particles
  ! MODFIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH
  ! PARTICLES
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================
 
  REAL(wp), INTENT(IN)    ::    &
     &      w_parc,             &    !< parcel updraft velocity
     &      s_parc,             &    !< supersaturation s_max
     &      nmb_conc(:),        &    !< total number concentration of particles
     &      sig_g(:)                 !< geometric dispersion
     
  REAL(wp), INTENT(OUT)   ::    &
     &      sinteg1,            &    !< Contribution of integral 1 for Kohler 
     &      sinteg2,            &    !< Contribution of integral 2 for kohler
     &      sinteg3,            &    !< Contribution of FHH integral
     &      sum_drpl                 !< Variable that stores all droplets
     
  INTEGER, INTENT(IN) ::    &
     &      nmodes                   !< Number of modes
           
  LOGICAL, INTENT(IN)     ::    &
     &      l_koehler(:)             !< If true, use koehler activation, else use fhh (dim: nmodes)
     
  REAL(wp), INTENT(OUT),OPTIONAL :: &
     &      n_scav(:), m_scav(:)     !< number and mass conc. of scavenged particles
     
  REAL(wp), INTENT(IN), OPTIONAL :: &
     &      mas_conc(:)              !< total mass concentration of particles
  
  INTEGER                 ::    &
     &      i                        !< counting variable
  
  REAL(wp)                ::    &
    &    c1,c2,c3,c4,           &    !< factors to calculate fhh exponent
    &    xfhh,                  &    !< fhh exponent
    &    descr,                 &    !< descriminant
    &    s_split,               &    !< 
    &    sqtwo,                 &    !< sqrt(2)
    &    dlgsg,                 &    !< ln(sigma_i) for Koehler and FHH
    &    dlgsp,                 &    !< ln(sg/smax) for Koehler and FHH
    &    orism1,orism2,         &    !< part of the integral for Koehler and FHH
    &    orism3,orism4,         &    !< part of the integral for Koehler and FHH
    &    orism5,orism6,         &    !< part of the integral for Koehler and FHH
    &    orism7,orism8,         &    !< part of the integral for FHH
    &    ekth,dw3,deq
    
  REAL(wp), ALLOCATABLE   ::    &
    &    integ1(:),integ2(:),   &    !< part of sinteg1 and sinteg2 for a single mode
    &    nmb_drpl(:)                 !< number of droplets for a single mode
  
  REAL(wp)         ::   &
    &    d11=-0.1907,   &
    &    d12=-1.6929,   &
    &    d13= 1.4963,   &
    &    d14=-0.5644,   &
    &    d15= 0.0711,   &
    &    d21=-3.9310,   &
    &    d22= 7.0906,   &
    &    d23=-5.3436,   &
    &    d24= 1.8025,   &
    &    d25=-0.2131,   &
    &    d31= 8.4825,   &
    &    d32=-14.9297,  &
    &    d33= 11.4552,  &
    &    d34=-3.9115,   &
    &    d35= 0.4647,   &
    &    d41=-5.1774,   &
    &    d42= 8.8725,   &
    &    d43=-6.8527,   &
    &    d44= 2.3514,   &
    &    d45=-0.2799    

  sinteg1   = 0._wp
  sinteg2   = 0._wp
  sum_drpl  = 0._wp
  sinteg3   = 0._wp
  sqtwo     = SQRT(2._wp)
  
  ALLOCATE ( integ1(nmodes),integ2(nmodes),nmb_drpl(nmodes))
    
  !-------------------------------------
  ! Determination of exponent (edit this in future: does not need to be done everytime)
  !-------------------------------------

  c1 = (d11)+(d12/a_fhh)+(d13/(a_fhh**2))+(d14/(a_fhh**3))+(d15/(a_fhh**4))
  c2 = (d21)+(d22/a_fhh)+(d23/(a_fhh**2))+(d24/(a_fhh**3))+(d25/(a_fhh**4))
  c3 = (d31)+(d32/a_fhh)+(d33/(a_fhh**2))+(d34/(a_fhh**3))+(d35/(a_fhh**4))
  c4 = (d41)+(d42/a_fhh)+(d43/(a_fhh**2))+(d44/(a_fhh**3))+(d45/(a_fhh**4))

  xfhh = (c1) + (c2/b_fhh) + (c3/(b_fhh**2)) + (c4/(b_fhh**3))
  
  !-------------------------------------
  ! Here is where the criterion with the descriminant is put.
  !-------------------------------------
  
  descr  = 1._wp - (16._wp/9._wp)*alfa*w_parc*bet2*(akoh/s_parc**2)**2

  IF (descr <= 0._wp) THEN
    s_split  = (2.e7_wp/3._wp)*akoh*s_parc**(-0.3824_wp)
    IF (s_split > 1._wp) s_split = 1._wp
    s_split = s_parc*s_split
  ELSE ! DR: As only the former SSPLIT2 was used, only this is left in the code
    s_split = 0.5_wp*(1._wp+SQRT(descr)) ! max root of both
    s_split = SQRT(s_split)*s_parc
  ENDIF
  
  !-------------------------------------
  ! Calculate integrals
  !-------------------------------------
  DO i = 1, nmodes
    IF (l_koehler(i)) THEN ! Koehler Modes
      dlgsg     = DLOG(sig_g(i))
      dlgsp     = DLOG(sg(i)/s_parc)
      orism1    = 2._wp*DLOG(sg(i)/s_split)/(3._wp*sqtwo*dlgsg) !upart
      orism2    = orism1 - 3._wp*dlgsg/(2._wp*sqtwo) 
      orism3    = 2._wp*dlgsp/(3._wp*sqtwo*dlgsg)-3._wp*dlgsg/(2._wp*sqtwo)
      orism4    = orism1 + 3._wp*dlgsg/sqtwo
      orism5    = 2._wp*dlgsp/(3._wp*sqtwo*dlgsg)
      ekth      = EXP(9._wp/2._wp*dlgsg*dlgsg)
      integ1(i) = nmb_conc(i)*s_parc*((1._wp-ERF(orism1)) - &
     &            0.5_wp*((sg(i)/s_parc)**2._wp)*ekth*(1._wp-ERF(orism4)))
      integ2(i) = (EXP(9._wp/8._wp*dlgsg*dlgsg)*nmb_conc(i)/sg(i)) &
     &          * (ERF(orism2) - ERF(orism3))
     
  !-------------------------------------
  ! Correction for kinetically limited large particles (D. Barahona et al. ACPD 2009)
  !-------------------------------------
  
      dlgsp     = DLOG(sg(i)/s_split)
      orism6    = (sqtwo*dlgsp/3._wp/dlgsg)-(1.5_wp*dlgsg/sqtwo)
      deq       = akoh*2._wp/sg(i)/3._wp/SQRT(3._wp)
      dw3       = nmb_conc(i)*deq*EXP(9._wp/8._wp*dlgsg*dlgsg)  &
     &          * (1._wp-ERF(orism6))*SQRT(bet2*alfa*w_parc)
     
  !-------------------------------------
  ! Calculate number of Drops
  !-------------------------------------

      nmb_drpl(i) = (nmb_conc(i)/2._wp)*(1._wp-ERF(orism5))


      sinteg1   = sinteg1  + integ1(i) + dw3      !SUM OF INTEGRAL 1 FOR KOHLER
      sinteg2   = sinteg2  + integ2(i)            !SUM OF INTEGRAL 2 FOR KOHLER
      sum_drpl  = sum_drpl + nmb_drpl(i)          !SUM OF ACTIVATED KOHLER PARTICLES
      
      ! l_do_scav does not work here
      IF(PRESENT(n_scav) .AND. PRESENT(m_scav) .AND. PRESENT(mas_conc)) THEN
        m_scav(i) = (mas_conc(i)/2._wp)*(1._wp-ERF(orism5-(3._wp*dlgsg/sqtwo)))
        n_scav(i) = nmb_drpl(i)
      ENDIF
  
    ELSE ! FHH
    
      dlgsg  = DLOG(sig_g(i))                  !ln(sigma,i)
      dlgsp  = DLOG(sg(i)/s_parc)              !ln(sg/smax)
      orism1 = (sg(i)*sg(i))/(s_parc*s_parc)     !(sg/smax)^2
      orism2 = EXP(2._wp*xfhh*xfhh*dlgsg*dlgsg) !exp(term)
      orism3 = sqtwo*xfhh*dlgsg             !sqrt(2).x.ln(sigma,i)
      orism4 = dlgsp/(-1._wp*orism3)           !Umax
      orism5 = orism3 - orism4
      orism6 = ERF(orism5)
      orism7 = 0.5_wp*orism1*orism2*(orism6+1._wp)
      orism8 = orism7 + ERF(orism4) - 1._wp

      integ1(i) =-1._wp*nmb_conc(i)*s_parc*orism8
      
  !-------------------------------------
  ! Correction for kinetically limited large particles (D. Barahona et al. ACPD 2009)
  !-------------------------------------
!      ORISM66  = (-dlgsp/xfhh/sqtwo/dlgsg)-(-xfhh*dlgsg/sqtwo)
!      dw3=nmb_conc(i)*Dpc(i)*EXP(xfhh**2*dlgsg*dlgsg)*
!     &    (1d0-ERF(ORISM66))*((bet2*alfa*w_parc)**0.5d0)   
  !-------------------------------------
  ! Calculate number of drops activated by FHH theory
  !-------------------------------------

      nmb_drpl(i) = (nmb_conc(i)/2._wp)*(1._wp-ERF(orism4))
      sinteg3     = sinteg3 + integ1(i) !+ dw3  !Sum of Integral 1 for FHH
      sum_drpl    = sum_drpl + nmb_drpl(i)         !Sum of ACTIVATED Kohler + FHH particles
      
      IF(PRESENT(n_scav) .AND. PRESENT(m_scav) .AND. PRESENT(mas_conc)) THEN
        m_scav(i) = (mas_conc(i)/2._wp)*(1._wp-erf(orism4-(3._wp*dlgsg/sqtwo)))
        n_scav(i) = nmb_drpl(i)
      ENDIF
    ENDIF
  ENDDO
  
  DEALLOCATE(integ1,integ2,nmb_drpl)
  
END SUBROUTINE sintegral
!
! ----------------------------------------------------------------------
!
SUBROUTINE gauleg

  !=================================================================
  ! Calculation of points and weights for N point GAUSS integration
  ! ---
  ! Rewritten in FORTRAN90 for the use in COSMO-ART by DANIEL RIEGER (KIT)
  !=================================================================

  INTEGER              ::   &
     &    j,i,              &   !< counting variable
     &    m                     !< limit
  
  REAL(wp)             ::   &
     &    x1=-1._wp,        &   !<
     &    x2= 1._wp,        &
     &    epsgaul=1.e-6_wp
     
  REAL(wp)             ::   &
     &    xm,xl,            &
     &    z,z1,             &
     &    p1,p2,p3,pp

  !-------------------------------------
  ! Calculation
  !-------------------------------------

  m  = (npgauss+1)/2
  xm = (x2+x1)/2._wp
  xl = (x2-x1)/2._wp
  DO i=1,m
    z  = COS(pi*(i-.25_wp)/(npgauss+.5_wp))
    z1 = 0.9_wp*z-epsgaul
    DO WHILE(ABS(z-z1) > epsgaul)
      p1=1._wp
      p2=0._wp
      DO j=1,npgauss
        p3=p2
        p2=p1
        p1=(((2._wp*REAL(j))-1._wp)*z*p2-(REAL(j)-1._wp)*p3)/REAL(j)
      ENDDO
      pp=npgauss*(z*p1-p2)/(z*z-1._wp)
      z1=z
      z=z1-p1/pp
    ENDDO
    xgs(i)=xm-xl*z
    xgs(npgauss+1-i)=xm+xl*z
    wgs(i)=2._wp*xl/((1._wp-z*z)*pp*pp)
    wgs(npgauss+1-i)=wgs(i)
  ENDDO
END SUBROUTINE gauleg
!
! ----------------------------------------------------------------------
!
END MODULE mo_art_nucleation_warm
