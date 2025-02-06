!
! mo_art_nucleation_interface
! This module serves as an interface between the two-moment scheme
! (mo_art_two_mom_main:art_clouds_twomoment) and the nucleation parameterizations
! for warm phase (mo_art_nucleation_warm:activation_master) and cold phase
! (mo_art_nucleation_cold:IceParam)
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

MODULE mo_art_nucleation_interface
! ICON
  USE mo_kind,                     ONLY: wp
  USE mo_exception,                ONLY: finish
  USE mo_satad,                    ONLY: sat_pres_water,  & ! saturation pressure over liquid water
    &                                    sat_pres_ice       ! saturation pressure over ice
  USE mo_physical_constants,       ONLY: rv
  USE mo_run_config,               ONLY: iqv,iqc,iqnc,iqi,iqni,iqns
! ART
  USE mo_art_data,                 ONLY: p_art_data
  USE mo_art_diag_types,           ONLY: t_art_diag
  USE mo_art_nucleation_warm,      ONLY: activation_master
  USE mo_art_nucleation_cold,      ONLY: IceParam
  USE mo_art_modes,                ONLY: t_fields_2mom, t_fields_1mom
  USE mo_art_modes_linked_list,    ONLY: p_mode_state,t_mode
  
  IMPLICIT NONE
    
  PRIVATE
  
  PUBLIC :: art_nuc_warm_interface
  PUBLIC :: art_nuc_cold_interface
  PUBLIC :: t_nuc_mode_cold
  
  
  TYPE t_nuc_mode_warm
    LOGICAL          :: do_nuc
    LOGICAL          :: l_koehler
    INTEGER          :: &
      &  i_numb,        & !< Index of tracer field with number density of this mode
      &  jsp(50),       & !< Array corresponding mass density species indizes are saved at
      &  njsp             !< Number of corresponding mass density species
    REAL(wp),POINTER :: &
      &  diam_pointer(:,:) !< Pointer to diameter field
    REAL(wp)         :: &
      &  sigma,         & !< Standard deviation
      &  number,        & !< Number concentration
      &  diameter,      & !< Median diameter
      &  dissfac_mean,  & !< Dissocitation factor
      &  molweight_mean,& !< Mean molercular weight
      &  solmassfr,     & !< Massfraction of soluble
      &  rhosol_mean,   & !< Density of the soluble fraction
      &  rhoinsol_mean    !< Density of the insoluble fraction
  END TYPE t_nuc_mode_warm
  
  TYPE t_nuc_mode_cold
    LOGICAL          :: &
      &  do_homfreez,   & !< Do homogeneous freezing
      &  is_dust,       & !< Do heterogeneous nucleation (Dust)
      &  is_soot,       & !< Do heterogeneous nucleation (Soot)
      &  is_org           !< Do heterogeneous nucleation (Organic material)
    INTEGER          :: &
      &  i_numb,        & !< Index of tracer field with number density of this mode
      &  jsp(50),       & !< Array corresponding mass density species indizes are saved at
      &  njsp             !< Number of corresponding mass density species
    REAL(wp),POINTER     :: &
      &  diam_pointer(:,:), & !< Pointer to diameter field
      &  mom3_pointer(:,:)    !< Pointer to third moment
    REAL(wp)         :: &
      &  sigma,         & !< Standard deviation
      &  number,        & !< Number concentration
      &  diameter         !< Median diameter
  END TYPE t_nuc_mode_cold
  
  TYPE(t_nuc_mode_warm),ALLOCATABLE :: &
    &  nuc_modes_warm(:)
  TYPE(t_nuc_mode_cold),ALLOCATABLE :: &
    &  nuc_modes_cold(:)
  
  REAL(KIND=wp), PARAMETER :: &
    &  wcb_min    = 0.1_wp,    & !< Minimum vertical velocity at cloud base 
                                 !    TODO: Check, should be lower in global model
    &  scb_min    = 0.0_wp,    & !< Minimum supersaturation at cloud base
    &  scb_max    =-0.0_wp,    & !< Maximum supersaturation at cloud base
    &  T_nuc      = 268.15_wp, & !< lower temperature threshold for ice nucleation, -5 C
    &  ni_het_max = 50.0e3_wp    !< max number of IN between 1-10 per liter, i.e. 1d3-10d3
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_nuc_warm_interface(w,t,pres,rho,tke,tkvh,dz,p_trac,jg,jb,istart,iend, &
  &                               kstart,kend,dtime,cloud_xmin,lextended_out)
!<
! SUBROUTINE art_nuc_warm_interface
! This SR sets name and the index of the number concentrations of
! the according mode and sets undefined values as initial values for
! other variables
! Based on: -
! Part of Module: mo_art_nucleation_interface
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  REAL(wp), INTENT(in), TARGET    :: &
    &  w(:,:),                       & !< Vertical velocity
    &  t(:,:),                       & !< Temperature
    &  pres(:,:),                    & !< Pressure
    &  rho(:,:),                     & !< Air density
    &  tke(:,:),                     & !< Turbulent kinetic energy
    &  tkvh(:,:),                    & !< turbulent diffusion coefficients for heat     (m/s2 )
    &  dz(:,:)                         !< height of layer
  REAL(wp), INTENT(inout), TARGET :: &
    &  p_trac(:,:,:)                   !< Tracer fields
  INTEGER, INTENT(in)             :: &
    &  jg,                           & !< Domain index
    &  jb,                           & !< block index
    &  istart, iend,                 & !< start/end indices
    &  kstart, kend                    !< start/end indices
  REAL(wp), INTENT(in)            :: &
    &  dtime,                        & !< Time step
    &  cloud_xmin                      !< Minimal mass of cloud droplet
  LOGICAL, INTENT(in)             :: &
    &  lextended_out                   !< Extended output variables
! NOTE: is this #if 0 really needed??????
#if 0
  ! Local variables
  TYPE(t_mode), POINTER           :: &
    &  this_mode
  TYPE(t_art_diag),POINTER        :: &
    &  art_diag                        !< Pointer to ART diagnostic fields
  REAL(wp)           :: &
    &  e_v,             & !< Partial pressure of water vapor
    &  e_s,             & !< Saturation pressure at ambient temperature
    &  s_sw,            & !< Supersaturation over water
    &  s_sw_kp1,        & !< Supersaturation over water one level below (i.e. at k+1)
    &  w_cb,            & !< Vertical velocity at cloud base
    &  solmass,         & !< Soluble mass of a mode
    &  wgt_mass,        & !< Mass weighted by molar mass
    &  totmass,         & !< Total mass of a mode (soluble+insoluble)
    &  volsol,          & !< Volume of the soluble part
    &  volinsol,        & !< Volume of the insoluble part
    &  ncn_tot,         & !< Total number of condensation nuclei
    &  n_act,           & !< Number of activated particles
    &  n_act_fhh,       & !< Number of activated particlesdue to FHH absorption theory
    &  s_max,           & !< Maximum supersaturation
    &  sig_w,           & !< Standard deviation of the vertical velocity PDF distribution
    &  adv_factor,      & !< Advection factor, depends on cloudtype
    &  npact,           & !< Number of activated particles
    &  nuc_n,nuc_q        !< Nucleated number and mass
  LOGICAL             :: &
    &  l_noaerosol
  LOGICAL,ALLOCATABLE :: &
    &  l_koehler(:)       !< True: Activation via Koehler theory, False: FHH absorption theory
  REAL(wp),ALLOCATABLE :: &
    &  dpg(:),            & !< Median diameter of the mode passed to activation
    &  amfs(:),           & !< Soluble mass fraction of the mode passed to activation
    &  vhf(:),            & !< Dissociation factor of the mode passed to activation
    &  ams(:),            & !< Mean molweight of the mode passed to activation
    &  rho_sol(:),        & !< Mean density of soluble fraction of the mode passed to activation
    &  rho_insol(:),      & !< Mean density of insoluble fraction of the mode passed to activation
    &  nmb_conc(:),       & !< Number concentration of the mode passed to activation
    &  sig_g(:)             !< Standard deviation of the mode passed to activation
  INTEGER            :: &
    &  ii,kk,k,         & !< Loop indices
    &  nmodes,          & !< Number of modes contained in mode list
    &  icloud_mode,     & !< Integer to have access to cloud droplets in nuc_modes_warm
                          !   for in-cloud activation
    &  imodes,itracer,  & !< Loop indice for modes and tracer
    &  cloudtype,       & !< Type of activation, 1: New cloud, 2: Cloud base, 3: In-cloud
    &  jsp,             & !< counter for species
    &  nact_modes,      & !< number of modes that contribute to activation
    &  iact               !< counter

  ! Associate Pointer for short references
  art_diag => p_art_data(jg)%diag
  
  nmodes = p_mode_state(jg)%p_mode_list%p%nmodes
  icloud_mode = nmodes+1 ! last spot is reserved for cloud mode
  ! Allocate object of structure t_nuc_mode_warm with nmodes + 1 (for cloud+rain droplets)
  ALLOCATE(nuc_modes_warm(nmodes+1))
  ! Initialization: Modes do not count as CN by default
  nuc_modes_warm(:)%do_nuc = .FALSE.
  nuc_modes_warm(:)%number = 0._wp
  
  imodes = 0
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(this_mode))
    imodes = imodes + 1
    ! Select type of mode // This may be not optimal to do a select type within a jb loop.
    ! But the object nuc_modes_warm cant be allocated outside the two-mom scheme and we would
    ! then suffer the problem that it needs the dimension jg
    ! ->Maybe a better solution can be found
    SELECT TYPE (fields=>this_mode%fields)
      CLASS is (t_fields_2mom)
        ! Calculate modal parameters (i.e. diameter, 3rd moment)
        CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),                &
          &                     istart, iend, kstart, kend, jb, p_trac(:,:,:))
        ! Some initializations
        nuc_modes_warm(imodes)%do_nuc       = .TRUE.
        nuc_modes_warm(imodes)%dissfac_mean = 1._wp
        nuc_modes_warm(imodes)%l_koehler    = .TRUE.
        ! The following values should be set by Metadata and may also not reflect the mean
        ! as in the future e.g. sulphate is allowed to condense on sea salt aerosol
        IF (TRIM(fields%name) == 'seasa') nuc_modes_warm(imodes)%dissfac_mean = 2._wp
        IF (TRIM(fields%name) == 'seasb') nuc_modes_warm(imodes)%dissfac_mean = 2._wp
        IF (TRIM(fields%name) == 'seasc') nuc_modes_warm(imodes)%dissfac_mean = 2._wp
        IF (TRIM(fields%name) == 'dusta') nuc_modes_warm(imodes)%l_koehler = .FALSE.
        IF (TRIM(fields%name) == 'dustb') nuc_modes_warm(imodes)%l_koehler = .FALSE.
        IF (TRIM(fields%name) == 'dustc') nuc_modes_warm(imodes)%l_koehler = .FALSE.
        nuc_modes_warm(imodes)%i_numb       =  fields%info%i_number_conc
        nuc_modes_warm(imodes)%sigma        =  fields%info%sg_ini
        nuc_modes_warm(imodes)%diam_pointer => fields%diameter(:,:,jb)
        nuc_modes_warm(imodes)%jsp(:)       =  fields%info%jsp(:)
        nuc_modes_warm(imodes)%njsp         =  fields%info%njsp
      CLASS is (t_fields_1mom) ! Modes with only mass concentration considered can not serve as CN
        nuc_modes_warm(imodes)%do_nuc   =.FALSE.
      CLASS DEFAULT
        ! Should not happen...
        CALL finish('mo_art_nucleation_interface:art_nuc_warm_interface', &
             &      'ART: Unknown class')
    END SELECT
    this_mode => this_mode%next_mode
  END DO
  
  
  ! One giant loop around the complete activation interface
  DO kk = kstart, kend
    DO ii = istart, iend
      
      ! Ideal gas law (calculate partial pressure of water vapor)
      e_v  = p_trac(ii,kk,iqv) * rv * t(ii,kk)
      ! Calculate saturation pressure over water
      e_s  = sat_pres_water(t(ii,kk))
      ! Calculate supersaturation
      s_sw = e_v / e_s - 1.0_wp
      ! And the same for the level below
      e_v  = p_trac(ii,MIN(kk+1,kend),iqv) * rv * t(ii,MIN(kk+1,kend))
      e_s  = sat_pres_water(t(ii,MIN(kk+1,kend)))
      s_sw_kp1 = e_v / e_s - 1.0_wp
    
      ! ----------------------------------
      ! --- Check for ambient conditions. (supersaturations and updraft velocity)
      ! ----------------------------------
      IF ( w(ii,kk) > wcb_min .AND. s_sw >= scb_min .AND.  s_sw > s_sw_kp1)  THEN
        w_cb = w(ii,kk)
      ELSE
        w_cb = 0.0_wp
      ENDIF
     
      ! Do nucleation only if ambient conditions fit
      IF (p_trac(ii,kk,iqc) > 0.0_wp .AND. w_cb > 0.0_wp &
         & .AND. t(ii,kk) > 223._wp .AND. t(ii,kk) < 323._wp ) THEN
        
        ! ----------------------------------
        ! --- Determine activation case:
        ! --- 1: New cloud
        ! --- 2: Cloud base
        ! --- 3: In-cloud
        ! ----------------------------------
        
        IF(s_sw_kp1 < scb_max) THEN
          cloudtype=2    ! cloud base
        ELSE
          cloudtype=3    ! incloud activation
        ENDIF
        ! If n_cloud 'too low' do activation as if new cloud
        IF( p_trac(ii,kk,iqnc) <= 10.0e6_wp )          THEN ! new cloud
          cloudtype=1
        ENDIF
        IF(cloudtype==2 .AND. kk<kend) THEN
          k = kk+1 ! below cloud aerosol activation one level below
        ELSE
          k = kk   ! In-Cloud or new cloud at this level
        ENDIF
        
        ! ----------------------------------
        ! --- Get the values needed by the nucleation parameterization and 
        ! --- store it inside nuc_modes_warm
        ! ----------------------------------
        
        nact_modes = 0
        
        DO imodes = 1, nmodes
          IF (nuc_modes_warm(imodes)%do_nuc) THEN
            nuc_modes_warm(imodes)%diameter = nuc_modes_warm(imodes)%diam_pointer(ii,kk)
            solmass  = 0._wp
            totmass  = 0._wp
            wgt_mass = 0._wp
            volsol   = 0._wp
            volinsol = 0._wp
            ! Sum up all species in mode
            DO itracer = 1, nuc_modes_warm(imodes)%njsp
              jsp = nuc_modes_warm(imodes)%jsp(itracer)
              IF(ASSOCIATED(p_art_data(jg)%p_meta(jsp)%p_aero)) THEN
                solmass  = solmass  + p_art_data(jg)%p_meta(jsp)%p_aero%solubility            &
                  &                 * p_trac(ii,kk,jsp)*rho(ii,kk)
                totmass  = totmass  + p_trac(ii,kk,jsp)*rho(ii,kk)
                wgt_mass = wgt_mass + p_trac(ii,kk,jsp)*rho(ii,kk)                            &
                  &                 / p_art_data(jg)%p_meta(jsp)%p_aero%mol_weight
                volsol   = volsol   + (p_art_data(jg)%p_meta(jsp)%p_aero%solubility           &
                  &                 * p_trac(ii,kk,jsp)*rho(ii,kk)                            &
                  &                 / p_art_data(jg)%p_meta(jsp)%p_aero%rho)
                volinsol = volinsol + ((1._wp - p_art_data(jg)%p_meta(jsp)%p_aero%solubility) &
                  &                 * p_trac(ii,kk,jsp)*rho(ii,kk)                            &
                  &                 / p_art_data(jg)%p_meta(jsp)%p_aero%rho)
              ELSE
                ! This should not happen, as the modes should have stored only indices of aerosols
                CALL finish('mo_art_nucleation_interface:art_nuc_warm_interface', &
                  &      'p_aero is not associated')
              ENDIF
            ENDDO
            nuc_modes_warm(imodes)%number  = p_trac(ii,kk,nuc_modes_warm(imodes)%i_numb)*rho(ii,kk)
            nuc_modes_warm(imodes)%molweight_mean = solmass / wgt_mass
            nuc_modes_warm(imodes)%solmassfr      = solmass / totmass
            nuc_modes_warm(imodes)%rhosol_mean    = solmass           / MAX(volsol  ,1.0e-20_wp)
            nuc_modes_warm(imodes)%rhoinsol_mean  = (totmass-solmass) / MAX(volinsol,1.0e-20_wp)
            IF(nuc_modes_warm(imodes)%number > 1.e6) THEN ! TUNING PARAMETER
              nact_modes = nact_modes + 1
            ELSE ! Save some computational time
              nuc_modes_warm(imodes)%do_nuc = .FALSE.
            ENDIF
          ENDIF ! %do_nuc
        ENDDO
        
        
        ! ----------------------------------
        ! --- Check if there is enough aerosol to perform the activation
        ! --- otherwise use prescribed continental scenario
        ! ----------------------------------
        ! TUNING PARAMETER: Maybe this determination is too hard
        IF (nact_modes == 0) THEN
          l_noaerosol = .TRUE.
          nact_modes = 1
        ELSE
          l_noaerosol = .FALSE.
        ENDIF
        
        ! ----------------------------------
        ! --- Consider cloud droplets in in-cloud activation calculation
        ! ----------------------------------
        
        IF(cloudtype==3) THEN
          nuc_modes_warm(icloud_mode)%do_nuc         = .TRUE.
          ! preexisting cloud droplets
          nuc_modes_warm(icloud_mode)%l_koehler      = .TRUE.
          nuc_modes_warm(icloud_mode)%number         = p_trac(ii,kk,iqnc)
          nuc_modes_warm(icloud_mode)%sigma          = 1.5_wp
          nuc_modes_warm(icloud_mode)%solmassfr      = 1.0_wp
          nuc_modes_warm(icloud_mode)%diameter       = 2.0e-06_wp
          nuc_modes_warm(icloud_mode)%dissfac_mean   = 1.0
          nuc_modes_warm(icloud_mode)%molweight_mean = 100.0e-03_wp
          nuc_modes_warm(icloud_mode)%rhosol_mean    = 1.8e+03_wp
          nuc_modes_warm(icloud_mode)%rhoinsol_mean  = 1.0e+03
          nact_modes = nact_modes + 1
        ELSE
          nuc_modes_warm(icloud_mode)%do_nuc = .FALSE.
          nuc_modes_warm(icloud_mode)%number = 0._wp
        ENDIF
        
        ncn_tot = 0._wp
        DO imodes = 1, (nmodes+1)
          IF (nuc_modes_warm(imodes)%do_nuc) THEN
            ncn_tot = ncn_tot + nuc_modes_warm(imodes)%number
          ENDIF
        ENDDO
        
        ! ----------------------------------
        ! --- At this place, certain predefined scenarios may be added in the future
        ! ----------------------------------
        
        ALLOCATE(l_koehler(nact_modes),dpg(nact_modes),amfs(nact_modes),vhf(nact_modes),          &
          &      ams(nact_modes),rho_sol(nact_modes),rho_insol(nact_modes),nmb_conc(nact_modes),  &
          &      sig_g(nact_modes))
        iact = 0
        DO imodes = 1, nmodes + 1
          IF (nuc_modes_warm(imodes)%do_nuc) THEN
            iact = iact + 1
            ! Debug check - If facing problems uncomment it!
            !IF (iact > nact_modes) THEN
            !  CALL finish('mo_art_nucleation_interface:art_nuc_warm_interface', &
            !   &      'ART: iact > nact_modes')
            !ENDIF
            l_koehler(iact)= nuc_modes_warm(imodes)%l_koehler
            dpg(iact)      = nuc_modes_warm(imodes)%diameter
            amfs(iact)     = nuc_modes_warm(imodes)%solmassfr
            vhf(iact)      = nuc_modes_warm(imodes)%dissfac_mean
            ams(iact)      = nuc_modes_warm(imodes)%molweight_mean
            rho_sol(iact)  = nuc_modes_warm(imodes)%rhosol_mean
            rho_insol(iact)= nuc_modes_warm(imodes)%rhoinsol_mean
            nmb_conc(iact) = nuc_modes_warm(imodes)%number
            sig_g(iact)    = nuc_modes_warm(imodes)%sigma
          ENDIF ! %do_nuc
        ENDDO
        
        IF (l_noaerosol) THEN 
          ! No prognostic aerosol available, so use a prescribed scenario (extreme maritime)
          ! according to Segal and Kain
          ! TUNING PARAMETER: What scenario is to be used?
          iact = iact + 1
          ! Debug check - If facing problems uncomment it!
          !IF (iact /= nact_modes) THEN
          !  CALL finish('mo_art_nucleation_interface:art_nuc_warm_interface', &
          !   &      'ART: iact /= nact_modes')
          !ENDIF
          l_koehler(iact)= .TRUE.
          dpg(iact)      = 0.08E-06_wp  !< Modal diameter (m)
          amfs(iact)     = 1.0_wp       !< Soluble mass fraction
          vhf(iact)      = 1.0_wp       !< Van't Hoff factor for soluble frac.(ions molec-1)
          ams(iact)      = 58.4E-03_wp  !< Molar mass of Soluble fraction (kg mol-1)
          rho_sol(iact)  = 2.2E+03_wp   !< Density of Soluble fraction (kg m-3)
          rho_insol(iact)= 2.2E+03_wp   !< Density of Insoluble fraction (kg m-3)
          nmb_conc(iact) = 100.0E+06_wp !< Total concentration (# m-3) extr mar
          sig_g(iact)    = 2.5E+00_wp   !< Geometric dispersion (sigma_g) extr mar
        ENDIF
        
        ! ----------------------------------
        ! --- Call the activation
        ! ----------------------------------
        
        sig_w = MAX(SQRT(2._wp*tke(ii,kk)) ,0.1_wp)
        
        
        CALL activation_master(nact_modes,l_koehler,dpg,amfs,vhf,ams,  &  !< INTENT(IN)
                &              rho_sol,rho_insol,                      &  !< INTENT(IN)
                &              t(ii,kk),pres(ii,kk),                   &  !< INTENT(IN)
                &              MAX(w(ii,kk),0.01_wp),sig_w,            &  !< INTENT(IN)
                &              nmb_conc,sig_g,                         &  !< INTENT(IN)
                &              n_act,n_act_fhh,s_max)                     !< INTENT(OUT)
        
        DEALLOCATE(l_koehler,dpg,amfs,vhf,ams,  &
                &  rho_sol,rho_insol,nmb_conc,sig_g)
        
        ! In case of in-cloud activation substract the cloud droplets again
        IF(cloudtype==3) THEN
          n_act=MAX(0.0_wp, (n_act-nuc_modes_warm(icloud_mode)%number) )
        ENDIF
        
        IF(cloudtype==1) adv_factor = 1.0_wp / dtime
        IF(cloudtype==3) adv_factor = 1.0_wp / dtime
        IF(cloudtype==2) THEN   ! should be changed to a consistent formulation with updraft PDF,
                                ! w -> w_eff?! without turb
          adv_factor = min( (w(ii,kk)   / (dz(ii,kk-1))                &
            &        + tkvh(ii,kk)                                     &
            &        /  (0.5_wp* (dz(ii,kk-1)) -  0.5_wp* (dz(ii,kk))) &
            &        /  (dz(ii,kk-1))) , 1.0_wp/ dtime)
        ENDIF !cloudtype==2
        
        
        npact = n_act * adv_factor
        
        nuc_n = dtime * npact ! calculate nucleation rate
        
        IF(nuc_n+p_trac(ii,kk,iqnc) > ncn_tot) THEN
          nuc_n=ncn_tot-p_trac(ii,kk,iqnc)
        ENDIF
        
        nuc_n = MAX(nuc_n,0.0_wp)
        nuc_q = MIN(nuc_n * cloud_xmin, p_trac(ii,kk,iqv))
        nuc_n = nuc_q / cloud_xmin
        
        p_trac(ii,kk,iqnc) = p_trac(ii,kk,iqnc) + nuc_n
        p_trac(ii,kk,iqc)  = p_trac(ii,kk,iqc)  + nuc_q
        p_trac(ii,kk,iqv)  = p_trac(ii,kk,iqv)  - nuc_q
        
        IF (lextended_out) THEN
          art_diag%ncalls_warm(ii,kk,jb)      = art_diag%ncalls_warm(ii,kk,jb)       &
             &                                + 1._wp
          art_diag%aci_nnuctot_warm(ii,kk,jb) = art_diag%aci_nnuctot_warm(ii,kk,jb)  &
             &                                + n_act
          art_diag%aci_nnucfhh_warm(ii,kk,jb) = art_diag%aci_nnucfhh_warm(ii,kk,jb)  &
             &                                + n_act_fhh
          art_diag%smax_water(ii,kk,jb)       = s_max
        ENDIF
        
      ELSE
        ! Constant cloud droplet number as with inwp_gscp=4
        p_trac(ii,kk,iqnc) = 200.0e6_wp
      ENDIF ! qc and w > 0. t between 223 and 323
    ENDDO !ii
  ENDDO !kk
  
  ! ----------------------------------
  ! --- Some clean up
  ! ----------------------------------
  DEALLOCATE(nuc_modes_warm)
#endif
END SUBROUTINE art_nuc_warm_interface
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_nuc_cold_interface(w,t,pres,rho,tke,p_trac,ice_xmin,istart,iend,kstart,kend,jg,jb, &
  &                               iaci_cold,lextended_out,n_inact)
!<
! SUBROUTINE art_nuc_cold_interface
! This SR prepares the aerosol for the usage in the nucleation routine
! Based on: -
! Part of Module: mo_art_nucleation_interface
! Author: Daniel Rieger, KIT
! Initial Release: 2014-11-10
! Modifications:
!! YYYY-MM-DD: <name>, <institution>
!! - ...
!>
  REAL(wp), INTENT(in), TARGET    :: &
    &  w(:,:),                       & !< Vertical velocity
    &  t(:,:),                       & !< Temperature
    &  pres(:,:),                    & !< Pressure
    &  rho(:,:),                     & !< Air density
    &  tke(:,:),                     & !< Turbulent kinetic energy
    &  ice_xmin                        !< Minimal mass of ice particle
  INTEGER, INTENT(in)     :: &
    &  istart,iend,          & !< Start and end of loop horizontal
    &  kstart,kend,          & !< Start and end of loop vertical
    &  jg,jb                   !< Domain index, block index
  REAL(wp), INTENT(inout), TARGET :: &
    &  p_trac(:,:,:)           !< Tracer fields
  INTEGER, INTENT(in)             :: &
    &  iaci_cold               !< Heterogeneous ice nucleation scheme
  LOGICAL, INTENT(in)             :: &
    &  lextended_out           !< Extended output variables
  REAL(wp),INTENT(inout),OPTIONAL :: &
    &  n_inact(:,:)            !< Number of already nucleated ice crystals (tracking variable to 
                               !    account for in-cloud scavenging), unit here: [m-3]
! Local Variables
  TYPE(t_mode), POINTER           :: this_mode
  TYPE(t_art_diag),POINTER        :: &
    &  art_diag                        !< Pointer to ART diagnostic fields
  INTEGER                 :: &
    &  ii, kk,               & !< Loop indizes
    &  nmodes,imodes,        & !< Number and counter of modes
    &  idx                     !< Indexing variable
  REAL(wp)                :: &
    &  inp_t, inp_p,        & !< Temperature and pressure prepared for call to nucleation
    &  nhom,dhom,           & !< Number and diameter of CN for homogeneous freezing
    &  mom3,                & !< Third moment of CN for homogeneous freezing
    &  esa36,dgmin,one3,    & !< Constants used for diameter calculation (details see definition
                              !   of variables)
    &  ndust(3),ddust(3),   & !< Number and diameter of mineral dust (dim=3 -> 3 dust modes)
    &  nsoot, dsoot,        & !< Number and diameter of soot
    &  norg, dorg,          & !< Number and diameter of organic particles
    &  sigdust(3),          & !< Standard deviations of mineral dust distributions
    &  sigsoot,sigorg,      & !< Standard deviations of soot and organic particle distributions
    &  inp_sigma,inp_miuv,  & !< Standard deviation and ? prepared for call to nucleation
    &  e_v,sigrid,          & !< Saturation vapor pressure/Supersaturation over ice (grid-scale)
    &  nuc_n,nuc_q            !< Nucleated number and mass
  LOGICAL                :: &
    &  lhet, lhom,          & !< Switches for pure(!) heterogeneous/homogeneous nucleation
    &  luse_prog_liquidaero   !< Use prognostic aerosol from ART
  REAL(wp)               :: &
    &  nice, nhet,          & !< Output of nucleation: Number of ice crystals, 
                              !     number of heterogeneously formed ice crystals
    &  smaxice, nlim          !< Output of nucleation: Maximum supersaturation over ice /
                              !     limiting number for homogeneous freezing
  
  ! ----------------------------------
  ! --- Initial Values / Constants
  ! ----------------------------------
  
  esa36   = EXP( 0.125_wp * LOG( 2.0_wp)**2 ) ** 36
  dgmin   = 1.0E-09_wp    !< Minimum diameter [m]
  one3    = 1._wp / 3._wp !< One third
  imodes  = 0
  
  luse_prog_liquidaero = .FALSE. ! currently, there is no sulfate aerosol in ART, so use prescribed scenario
  
  ! Associate Pointer for short references
  art_diag => p_art_data(jg)%diag
  ! ----------------------------------
  ! --- Start routine
  ! ----------------------------------
  
  ! Drieg: Include something like aerostart?
  ! Drieg: Predefined scenarios?
  
  ! Loop through modes
  nmodes = p_mode_state(jg)%p_mode_list%p%nmodes
  ! Allocate object of structure t_nuc_mode_warm with nmodes + 1 (for cloud+rain droplets)
  ALLOCATE(nuc_modes_cold(nmodes))
  
  this_mode => p_mode_state(jg)%p_mode_list%p%first_mode
  DO WHILE(ASSOCIATED(this_mode))
    imodes = imodes + 1
    ! Select type of mode // This may be not optimal to do a select type within a jb loop.
    ! But the object nuc_modes_cold cant be allocated outside the two-mom scheme and we would
    ! then suffer the problem that it needs the dimension jg
    ! ->Maybe a better solution can be found
    SELECT TYPE (fields=>this_mode%fields)
      CLASS is (t_fields_2mom)
        ! Calculate modal parameters (i.e. diameter, 3rd moment)
        CALL fields%modal_param(p_art_data(jg)%air_prop%art_free_path(:,:,jb),                &
          &                     istart, iend, kstart, kend, jb, p_trac(:,:,:))
        ! Drieg: do_homfreez/do_hetnuc needs to be set via Metadata in the future!!!
        nuc_modes_cold(imodes)%do_homfreez = .TRUE.
        nuc_modes_cold(imodes)%is_dust     = .FALSE.
        nuc_modes_cold(imodes)%is_soot     = .FALSE.
        nuc_modes_cold(imodes)%is_org      = .FALSE.
        ! Exceptions (Not necessary if using Metadata in the future)
        IF (TRIM(fields%name) == 'dusta') THEN 
          nuc_modes_cold(imodes)%do_homfreez = .FALSE.
          nuc_modes_cold(imodes)%is_dust     = .TRUE.
        ENDIF
        IF (TRIM(fields%name) == 'dustb') THEN
          nuc_modes_cold(imodes)%do_homfreez = .FALSE.
          nuc_modes_cold(imodes)%is_dust     = .TRUE.
        ENDIF
        IF (TRIM(fields%name) == 'dustc') THEN
          nuc_modes_cold(imodes)%do_homfreez = .FALSE.
          nuc_modes_cold(imodes)%is_dust     = .TRUE.
        ENDIF
        nuc_modes_cold(imodes)%i_numb   = fields%itr0
        nuc_modes_cold(imodes)%sigma    = fields%info%sg_ini
        nuc_modes_cold(imodes)%diam_pointer => fields%diameter(:,:,jb)
        nuc_modes_cold(imodes)%mom3_pointer => fields%third_moment(:,:,jb)
        nuc_modes_cold(imodes)%jsp(:)   = fields%itr3(:)
        nuc_modes_cold(imodes)%njsp     = fields%ntr-1
      CLASS is (t_fields_1mom)
        nuc_modes_cold(imodes)%do_homfreez = .FALSE.
        nuc_modes_cold(imodes)%is_dust     = .FALSE.
        nuc_modes_cold(imodes)%is_soot     = .FALSE.
        nuc_modes_cold(imodes)%is_org      = .FALSE.
      CLASS DEFAULT
        ! Should not happen...
        CALL finish('mo_art_nucleation_interface:art_nuc_cold_interface', &
             &      'ART: Unknown class')
    END SELECT
    this_mode => this_mode%next_mode
  ENDDO
  
  
  
  ! One giant loop around the complete activation interface
  DO kk = kstart, kend
    DO ii = istart, iend
      
      ! Calculate sigrid
      e_v    = p_trac(ii,kk,iqv)*t(ii,kk)*rv
      sigrid = e_v/sat_pres_ice(t(ii,kk)) - 1.0_wp
      
      IF (t(ii,kk) < T_nuc .AND. sigrid > 0.0_wp &
              .AND. ( (p_trac(ii,kk,iqni)+p_trac(ii,kk,iqns)) < ni_het_max ) )THEN
        
        ! ----------------------------------
        ! --- Get aerosols number/diameter which is 
        ! --- available for homogeneous freezing
        ! ----------------------------------
        ! Note: At this point a rough approximation is performed:
        !       The Barahona and Nenes ice nucleation routine can not
        !       handle the competition of multimodal distributions for
        !       homogeneous freezing. Therefore, all modes are combined
        !       into one log-normal distribution (characterized by 
        !       number and diameter) to be passed to the Barahona and 
        !       Nenes routines.
        ! ----------------------------------
        
        IF (luse_prog_liquidaero) THEN ! Use prognostic liquid aerosol from ART
          nhom = 0.0_wp
          dhom = 0.0_wp
          mom3 = 0.0_wp
          
          DO imodes = 1, nmodes
            IF (nuc_modes_cold(imodes)%do_homfreez) THEN
              nhom = nhom + p_trac(ii,kk,nuc_modes_cold(imodes)%i_numb) * rho(ii,kk)
              mom3 = mom3 + nuc_modes_cold(imodes)%mom3_pointer(ii,kk) * rho(ii,kk)
            ENDIF
          ENDDO
          
          IF(nhom>10.0_wp) THEN
            dhom =  MAX(dgmin,( mom3 / ( nhom * esa36 ) ) ** one3)
          ELSE ! Use minim aerosol concentration 10 m-3
            nhom = 10.0_wp
            dhom = 200.0e-09_wp
          ENDIF
        ELSE
          ! Regime that is not aerosol-limited ("typical number density according to 
          !                                      Koehler and Seifert, 2015)
          nhom = 1.0e+09_wp ! 1000 cm -3 
          dhom = 200.0e-09_wp
        ENDIF
        ! ----------------------------------
        ! --- Get aerosols number/diameter which is 
        ! --- available for heterogeneous nucleation:
        ! --- Dust, Soot, Organic Particles
        ! ----------------------------------
        
        ndust(:)   = 0._wp
        ddust(:)   = 0._wp
        sigdust(:) = 0._wp
        
        idx = 1
        DO imodes = 1, nmodes
          IF (nuc_modes_cold(imodes)%is_dust) THEN
            IF (idx > 3) THEN
              CALL finish('mo_art_nucleation_interface:art_nuc_cold_interface', &
                   &      'ART: More dust species than 3 found')
            ENDIF
            ndust(idx)   = p_trac(ii,kk,nuc_modes_cold(imodes)%i_numb)* rho(ii,kk)
            ddust(idx)   = nuc_modes_cold(imodes)%diam_pointer(ii,kk)
            sigdust(idx) = nuc_modes_cold(imodes)%sigma
            idx = idx+1
          ENDIF
        ENDDO
        
        IF (idx == 1) THEN ! No dust mode has been found -> Predefined scenario
          ndust(1)   = 1.0e+06_wp
          ndust(2)   = 1.0e+06_wp
          ndust(3)   = 1.0e+06_wp
          ddust(1)   = 200.0e-09_wp
          ddust(2)   = 400.0e-09_wp
          ddust(3)   = 600.0e-09_wp
          sigdust(1) = 1.7_wp
          sigdust(2) = 1.6_wp
          sigdust(3) = 1.5_wp
        ENDIF
        
        IF (PRESENT(n_inact)) THEN
          ! If there are more activated IN than IN available, 
          !set the activated to the value of the available IN
          n_inact(ii,kk) = MIN(n_inact(ii,kk),(ndust(1)+ndust(2)+ndust(3))) ! Units all: [m-3]
          ! At this point we assume, that the dust mode indices are sorted according
          ! to their diameter. The mode with the highest diameter has the highest index.
          ndust(3) = ndust(3) - n_inact(ii,kk)
          IF (ndust(3) < 0._wp) THEN
            ndust(2) = ndust(2) + ndust(3) ! ndust(3) is negative
            ndust(3) = 0._wp
            IF (ndust(2) < 0._wp) THEN
              ndust(1) = ndust(1) + ndust(2) ! ndust(2) is negative
              ndust(2) = 0._wp
              IF (ndust(1) < 0._wp) THEN
                ! Should not happen, but to be safe...
                ndust(1) = 0._wp
              ENDIF ! ndust(1) < 0
            ENDIF ! ndust(2) < 0
          ENDIF ! ndust(3) < 0
        ENDIF ! PRESENT(n_inact)
        
        ! Drieg: At this point, nsoot, dsoot, norg, dorg have to be calculated.
        !        As they are not part of ICON-ART yet, this is skipped
        !        We need some kind of profile for the time being
        
        nsoot   = 0.0_wp
        dsoot   = 10.0e-09_wp
        norg    = 0.0_wp
        dorg    = 1.0e-09_wp
        sigorg  = 1.5_wp
        sigsoot = 1.4_wp
        
        ! ----------------------------------
        ! --- Setup atmospheric conditions
        ! ----------------------------------
        
        inp_t     =   t(ii,kk)
        inp_p     =   pres(ii,kk)
        inp_sigma =   0.3_wp * SQRT(tke(ii,kk)) ! see Dissertation Rieger
        inp_miuv  =   w(ii,kk)
        
        ! ----------------------------------
        ! --- Freezing mechanisms
        ! ----------------------------------
        
        IF(p_trac(ii,kk,iqc) > 0.0_wp) THEN
          lhet = .TRUE.  ! pure het. freezing  maybe not neccessary
          lhom = .FALSE.
        ELSE ! mixed freezing
          lhet = .FALSE. 
          lhom = .FALSE.
        ENDIF
        
        ! Reset output variables
        nhet   = 0.0_wp
        nice   = 0.0_wp
        smaxice= 0.0_wp
        nlim   = 0.0_wp
        
        
        CALL IceParam (inp_miuv,inp_sigma, inp_T, inp_p,sigrid,            & !<  input
                       nhom,dhom, ndust, ddust, nsoot, dsoot, norg, dorg,  & !<  input
                       lhet, lhom,sigdust,sigsoot,sigorg,                  & !<  input
                       nhet, nice, smaxice, nlim,                          & !< output
                       iaci_cold)           !< Type of heterogeneous nucleation scheme
        
        ! ----------------------------------
        ! --- Update mass and number concentrations of ice and vapor
        ! ----------------------------------
        nuc_n = MAX(nice - (p_trac(ii,kk,iqni)+p_trac(ii,kk,iqns)),0._wp)
        nuc_q = MIN(nuc_n * ice_xmin, p_trac(ii,kk,iqv))
        p_trac(ii,kk,iqni) = p_trac(ii,kk,iqni) + nuc_n
        p_trac(ii,kk,iqi)  = p_trac(ii,kk,iqi)  + nuc_q
        p_trac(ii,kk,iqv)  = p_trac(ii,kk,iqv)  - nuc_q
        
        IF (PRESENT(n_inact)) THEN
          ! Add number of heterogeneously nucleated ice crystals to tracking variable n_inact
          n_inact(ii,kk) = n_inact(ii,kk) + nhet ! Units: [m-3]
        ENDIF
        
        IF (lextended_out) THEN
          art_diag%ncalls_cold(ii,kk,jb)      = art_diag%ncalls_cold(ii,kk,jb)       &
             &                                + 1._wp
          art_diag%aci_nnuctot_cold(ii,kk,jb) = art_diag%aci_nnuctot_cold(ii,kk,jb)  &
             &                                + nice
          art_diag%aci_nnuchet_cold(ii,kk,jb) = art_diag%aci_nnuchet_cold(ii,kk,jb)  &
             &                                + nhet
          art_diag%smax_ice(ii,kk,jb)         = smaxice
        ENDIF
        
      ENDIF ! t, sigrid, ni_het_max
    ENDDO ! ii
  ENDDO ! kk
  
  DEALLOCATE(nuc_modes_cold)
END SUBROUTINE art_nuc_cold_interface
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_nucleation_interface
