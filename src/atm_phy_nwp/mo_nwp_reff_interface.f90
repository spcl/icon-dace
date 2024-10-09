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

! Interface to compute effective radius consistent with microphysics, cloud scheme
! and convection scheme choice (not yet!).
!
! Description:
! The module provides the interface between the ICON routine and the effective
! radius module

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_nwp_reff_interface

  USE mo_kind,                 ONLY: wp, i4
  USE mo_exception,            ONLY: finish, message, message_text

  USE mo_run_config,           ONLY: msg_level, iqc, iqi, iqr, iqs,       &
                                       iqni, iqg, iqh, iqnr, iqns,     &
                                       iqng, iqnh ,iqnc

  USE mo_nonhydro_types,       ONLY: t_nh_prog,t_nh_diag,t_nh_metrics
  USE mo_nwp_phy_types,        ONLY: t_nwp_phy_diag
  USE mo_parallel_config,      ONLY: nproma
  USE mo_model_domain,         ONLY: t_patch
  USE mo_ext_data_types,       ONLY: t_external_data
  USE mo_impl_constants,       ONLY: min_rlcell_int, idu
  USE mo_impl_constants_grf,   ONLY: grf_bdywidth_c
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_nonhydrostatic_config,ONLY: kstart_moist
  USE mo_atm_phy_nwp_config,   ONLY: atm_phy_nwp_config, iprog_aero, icpl_aero_ice
  USE mo_radiation_config,     ONLY: irad_aero, iRadAeroTegen, iRadAeroCAMSclim, iRadAeroCAMStd
  USE mo_nwp_tuning_config,    ONLY: tune_zceff_min, tune_v0snow, tune_zvz0i, tune_icesedi_exp

  USE mo_reff_types,           ONLY: t_reff_calc_dom,  nreff_max_calc
  USE mo_reff_main,            ONLY: init_reff_calc, mapping_indices, calculate_ncn, calculate_reff, combine_reff, set_max_reff
  USE mo_impl_constants,       ONLY: max_dom  
  USE microphysics_1mom_schemes, ONLY: microphysics_1mom_init

  IMPLICIT NONE
  PRIVATE
  
  PUBLIC:: set_reff, init_reff, reff_calc_dom, combine_phases_radiation_reff

  TYPE(t_reff_calc_dom) ::  reff_calc_dom(max_dom)        ! Calculations array 

  CONTAINS
! ------------------------------------------------------------------------------------------

! Initialize parameters for reff calculations. 
! It should be calculated only once after the microphysics
  SUBROUTINE init_reff (prm_diag, p_patch, p_prog, return_reff) 
    
    TYPE(t_nwp_phy_diag),INTENT(in)    :: prm_diag         ! Diagnostic statistics from physics (inlcuding reff)
    TYPE(t_patch)       ,INTENT(in)    :: p_patch          ! Grid/patch info.
    TYPE(t_nh_prog)     ,INTENT(in)    :: p_prog           ! The dyn prog vars
    LOGICAL, OPTIONAL   ,INTENT(out)   :: return_reff      ! Return call from function .true. if right

    ! Local scalars:
    INTEGER :: icalc_reff_loc   ! Predefined set of reff parame (usually equal to igscp) (NAMELIST)
    INTEGER :: igscp            ! Microphysical scheme called by microphysics
    INTEGER :: jg               ! patch indices
    INTEGER :: ii               ! counter
    LOGICAL :: available_acdnc   ! Available cloud water cloud droplet concentration from radiation
    LOGICAL :: return_fct       ! Return values

    INTEGER :: nreff_calc       ! Number of effective radius calculations 
    
    return_fct = .true.

    ! domain ID
    jg = p_patch%id

    !$ACC ENTER DATA CREATE(reff_calc_dom(jg:jg))


! The initialization is needed only once and could be set in a different routine

! Initialize 1 moment scheme coefficients in case, they were not inititated
    SELECT CASE (  atm_phy_nwp_config(jg)%inwp_gscp )
    CASE (1,2)
      IF (msg_level >= 15) THEN 
        WRITE (message_text,*) "Reff: one-moment scheme already initialized"
        CALL message('',message_text)
      END IF

    CASE DEFAULT  ! Initialize 1 moment with graupel schme (igscp=2) for subgrid clouds
      CALL microphysics_1mom_init(            igscp = 2,                             & 
           &                        tune_zceff_min = tune_zceff_min,                &
           &                        tune_v0snow    = tune_v0snow,                   &
           &                        tune_zvz0i     = tune_zvz0i,                    &
           &                      tune_icesedi_exp = tune_icesedi_exp,              &
           &                        tune_mu_rain   = atm_phy_nwp_config(jg)%mu_rain,&
           &                   tune_rain_n0_factor = atm_phy_nwp_config(jg)%rain_n0_factor)
    END SELECT

    IF ( ANY ( atm_phy_nwp_config(jg)%icpl_aero_gscp == (/1, 3/) ) ) THEN  ! Only defined if aerosol coupling is on
      IF (iprog_aero == 0) THEN  ! Take CCN from cloud_num or acdnc
        available_acdnc = .true.
      ELSE  ! Not yet developed fucntion
        WRITE (message_text,*) 'Warnning Reff: 1 mom cannot generate cloud droplet &
                            &numbers for current options and icpl_aero_gscp=1.     &
                            &Using constant number'

        ! NOT YET IMPLEMENTED MAYBE EXTRACT TO INDEPENDENT FUNCTION (Called by radiation, microphysics and reff)
        !     ELSE   
        
        !         CALL ncn_from_tau_aerosol_speccnconst (nproma, nlev, i_startidx, i_endidx, nlev, nlev, &
        !              p_metrics%z_ifc(:,:,jb), prm_diag%aerosol(:,iss,jb), prm_diag%aerosol(:,iso4,jb),    &
        !              prm_diag%aerosol(:,iorg,jb), prm_diag%aerosol(:,idu,jb), zncn)
        
        !         CALL specccn_segalkhain_simple (nproma, i_startidx, i_endidx, zncn(:,nlev), prm_diag%cloud_num(:,jb))
        !         qnc_s(i_startidx:i_endidx) = prm_diag%cloud_num(i_startidx:i_endidx,jb)


        CALL message('',message_text)
        available_acdnc = .false.
      END IF
    ELSE 
      IF ( ASSOCIATED ( prm_diag%acdnc ) ) THEN
        available_acdnc = .true.
      ELSE ! Use constant number cloud_num if acdnc is not allocated        
        available_acdnc = .false.
        WRITE (message_text,*) 'Warnning Reff: 1 mom cannot generate cloud droplet &
                            &numbers for current options (acdnc is not allocated). &
                            &Using constant number'
      END IF
    END IF

    ! Current Microphysical scheme
    igscp =  atm_phy_nwp_config(jg)%inwp_gscp

    ! Refff paramaterization consistent with microphysics
    IF (  atm_phy_nwp_config(jg)%icalc_reff == 100) THEN
      icalc_reff_loc =  igscp
    ELSE
      icalc_reff_loc =   atm_phy_nwp_config(jg)%icalc_reff
    END IF


    ! Maximum 10 parameterizations. If you need more, change nreff_max_calc in mo_reff_types.f90
    nreff_calc = 0

    SELECT CASE ( icalc_reff_loc )


    ! RRTM. 
    CASE(101)
      ! It uses acdns for number concentration of cloud water. No number concentration of ice needed.
      ! No distinction between grid and subgrid
      nreff_calc = nreff_calc +  1

      IF ( available_acdnc ) THEN 
        CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                      &     hydrometeor   = 0,                             & ! Cloud Water DSD and geom. properties
                      &     grid_scope    = 0,                             & ! Grid and Subgrid
                      &     microph_param = 101,                           & ! RRTM Param
                      &     return_fct    = return_fct,                    & ! Return parameter
                      &     p_q           = p_prog%tracer(:,:,:,iqc),      & ! Grid Cloud water
                      &     p_qtot        = prm_diag%tot_cld(:,:,:,iqc),   & ! Total Cloud water
                      &     p_ncn3D       = prm_diag%acdnc(:,:,:),         & ! Number concentration acdnc
                      &     reff_param    = 0 ,                            & ! Spheroid
                      &     ncn_param     = 101,                           & ! External acdnc field for ncn (default RRTM)
                      &     p_reff        = prm_diag%reff_qc(:,:,:) )        ! Output
      ELSE 
        WRITE (message_text,*) 'WANING: Reff not defined for RRTM when acdnc is not available.'
        CALL message('',message_text)
        IF ( PRESENT (return_reff) ) return_reff = .false.  
      END IF


      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 1,                             & ! Ice DSD and geom. properties
                    &     grid_scope    = 0 ,                            & ! Grid and Subgrid
                    &     microph_param = 101,                           & ! RRTM Param
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqi),      & ! Grid Ice 
                    &     p_qtot        = prm_diag%tot_cld(:,:,:,iqi),   & ! Total Ice
                    &     p_reff        = prm_diag%reff_qi(:,:,:) )        ! Output



    ! 1 Moment Scheme
    CASE (1,2)
      ! Cloud water using mono-modal spherical particles and ncnd from acdnc (different from standard with cloud_num)
      ! No distinction between grid and subgrid
      nreff_calc = nreff_calc + 1
      IF ( available_acdnc ) THEN 
        CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                      &     hydrometeor   = 0,                             & ! Cloud Water (here overruled by microph_param)
                      &     grid_scope    = 0,                             & ! Grid and Subgrid
                      &     microph_param = 100,                           & ! Spherical particles
                      &     dsd_type      = 2,                             & ! Generalized gamma distribution
                      &     nu            = 5.0_wp,                        & ! Consisitent with nu_mass = 1.0 in 2 moment 
                      &     mu            = 3.0_wp,                        & ! Consisitent with mu_mass = 1.0 in 2 moment 
                      &     return_fct    = return_fct,                    & ! Return parameter
                      &     p_q           = p_prog%tracer(:,:,:,iqc),      & ! Total Cloud water
                      &     p_qtot        = prm_diag%tot_cld(:,:,:,iqc),   & ! Grid Cloud water
                      &     p_ncn3D       = prm_diag%acdnc(:,:,:),         & ! Number concentration acdnc
                      &     reff_param    = 0,                             & ! Spheroid
                      &     ncn_param     = 101,                           & ! External acdnc field for ncn
                      &     p_reff        = prm_diag%reff_qc(:,:,:) )        ! Output
     ! Constant cloud number
      ELSE 
        CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                      &     hydrometeor   = 0,                             & ! Cloud Water (here overruled by microph_param)
                      &     grid_scope    = 0,                             & ! Grid and Subgrid
                      &     microph_param = 100,                           & ! Spherical particles
                      &     dsd_type      = 2,                             & ! Generalized gamma distribution
                      &     nu            = 5.0_wp,                        & ! Consisitent with nu_mass = 1.0 in 2 moment 
                      &     mu            = 3.0_wp,                        & ! Consisitent with mu_mass = 1.0 in 2 moment 
                      &     return_fct    = return_fct,                    & ! Return parameter
                      &     p_q           = p_prog%tracer(:,:,:,iqc),      & ! Total Cloud water
                      &     p_qtot        = prm_diag%tot_cld(:,:,:,iqc),   & ! Grid Cloud water
                      &     reff_param    = 0,                             & ! Spheroid
                      &     ncn_param     = 0,                             & ! Constant ccn (given by cloud_num)
                      &     p_reff        = prm_diag%reff_qc(:,:,:) )        ! Output

      ENDIF

      ! Ice using standard 1mom (no distinction between grid and subgrid
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 1,                             & ! Ice DSD and geom. properties
                    &     grid_scope    = 0 ,                            & ! Grid and Subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 1 moment scheme
                    &     return_fct    = return_fct,                    & ! Return function
                    &     p_q           = p_prog%tracer(:,:,:,iqi),      & ! Total Ice
                    &     p_qtot        = prm_diag%tot_cld(:,:,:,iqi),   & ! Total Cloud water
                    &     reff_param    = 1 ,                            & ! Fu
                    &     ncn_param     = 1,                             & ! 1 moment ncn
                    &     p_reff        = prm_diag%reff_qi(:,:,:) )        ! Output
    
      ! Rain using 1 mom 
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 2,                             & ! Rain DSD and geom. properties
                    &     grid_scope    = 0 ,                            & ! Grid and Subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 1 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqr),      & ! Rain mixing ratio
                    &     reff_param    = 0 ,                            & ! Spheroid
                    &     ncn_param     = 1,                             & ! 1 moment ncn
                    &     p_reff        = prm_diag%reff_qr(:,:,:) )        ! Output

      ! Snow using 1 mom
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 3,                             & ! Snow DSD and geom. properties
                    &     grid_scope    = 0 ,                            & ! Grid and Subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 1 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqs),      & ! Snow mixing ratio
                    &     reff_param    = 1 ,                            & ! Fu
                    &     ncn_param     = 1,                             & ! 1 moment ncn
                    &     p_reff        = prm_diag%reff_qs(:,:,:) )        ! Output

      ! One moment scheme with graupel
      IF ( icalc_reff_loc == 2 ) THEN
        ! Graupel using 1 mom
        nreff_calc = nreff_calc + 1
        CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                      &     hydrometeor   = 4,                             & ! Graupel DSD and geom. properties
                      &     grid_scope    = 0 ,                            & ! Grid and Subgrid
                      &     microph_param = icalc_reff_loc,                & ! Chosen 1 moment scheme
                      &     return_fct    = return_fct,                    & ! Return parameter
                      &     p_q           = p_prog%tracer(:,:,:,iqg),      & ! Graupel mixing ratio
                      &     reff_param    = 1 ,                            & ! Fu
                      &     ncn_param     = 1,                             & ! 1 moment ncn
                      &     p_reff        = prm_diag%reff_qg(:,:,:) )        ! Output
      END IF

      ! 2 Moment Scheme
    CASE ( 4,5,6,7,8) 

      ! Grid cloud water from 2 moment scheme
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 0,                             & ! Clour Water DSD and geom. properties
                    &     grid_scope    = 1,                             & ! Grid clouds only
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqc),      & ! Grid Cloud water
                    &     p_qtot        = prm_diag%tot_cld(:,:,:,iqc),   & ! Total Cloud water
                    &     p_ncn3D       = p_prog%tracer(:,:,:,iqnc),     & ! Number concentration from 2 mom
                    &     reff_param    = 0,                             & ! Spheroid
                    &     ncn_param     = 4,                             & ! 2 moment ncn
                    &     p_reff        = prm_diag%reff_qc(:,:,:) )        ! Output

      ! Subgrid cloud water (same geometry and DSD, ncn from acdnc)
      nreff_calc = nreff_calc + 1
       IF ( available_acdnc ) THEN 
         CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                      &     hydrometeor   = 0,                              & ! Cloud Water DSD and geom. properties
                      &     grid_scope    = 2 ,                             & ! Subgrid clouds only
                      &     microph_param = icalc_reff_loc,                 & ! Chosen 2 moment scheme
                      &     return_fct    = return_fct,                     & ! Return parameter
                      &     p_q           = p_prog%tracer(:,:,:,iqc),       & ! Grid Cloud water
                      &     p_qtot        = prm_diag%tot_cld(:,:,:,iqc),    & ! Total Cloud water
                      &     p_ncn3D       = prm_diag%acdnc(:,:,:),          & ! Number concentration from acdnc
                      &     reff_param    = 0,                              & ! Spheroid
                      &     ncn_param     = 101,                            & ! External acdnc field for ncn
                      &     p_reff        = prm_diag%reff_qc(:,:,:) )         ! Output
       ! Constant ncn
       ELSE
         CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                      &     hydrometeor   = 0,                              & ! Cloud Water DSD and geom. properties
                      &     grid_scope    = 2 ,                             & ! Subgrid clouds only
                      &     microph_param = icalc_reff_loc,                 & ! Chosen 2 moment scheme
                      &     return_fct    = return_fct,                     & ! Return parameter
                      &     p_q           = p_prog%tracer(:,:,:,iqc),       & ! Grid Cloud water
                      &     p_qtot        = prm_diag%tot_cld(:,:,:,iqc),    & ! Total Cloud water
                      &     reff_param    = 0,                              & ! Spheroid
                      &     ncn_param     = 0,                              & ! Constant ccn (given by cloud_num)
                      &     p_reff        = prm_diag%reff_qc(:,:,:) )         ! Output
       END IF

      ! Grid ice from 2 moment scheme
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 1,                             & ! Ice DSD and geom. properties
                    &     grid_scope    = 1 ,                            & ! Grid clouds only
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqi),      & ! Total ice
                    &     p_qtot        = prm_diag%tot_cld(:,:,:,iqi),   & ! Grid ice
                    &     p_ncn3D       = p_prog%tracer(:,:,:,iqni),     & ! Number concentration from 2 mom
                    &     reff_param    = 1 ,                            & !  Fu param
                    &     ncn_param     = 4,                             & ! 2 mom ncn
                    &     p_reff = prm_diag%reff_qi(:,:,:) )               ! Output

      ! Sub Grid ice with 2 moment microphysics and DSD and ncn Cooper formula from 1D scheme
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 1,                             & ! Ice DSD and geom. properties
                    &     grid_scope    = 2 ,                            & ! Subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqi),      & ! Grid ice
                    &     p_qtot        = prm_diag%tot_cld(:,:,:,iqi),   & ! Total ice
                    &     reff_param    = 1 ,                            & ! Fu param
                    &     ncn_param     = 1,                             & ! ncn from 1 moment scheme
                    &     p_reff        = prm_diag%reff_qi(:,:,:) )        ! Output

      ! Rain
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 2,                             & ! Rain DSD and geom. properties
                    &     grid_scope    = 0,                             & ! Grid and subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqr),      & ! Rain mixing ratio
                    &     p_ncn3D       =  p_prog%tracer(:,:,:,iqnr),    & ! Rain Number concentration
                    &     reff_param    = 0,                             & ! Spheroid
                    &     ncn_param     = 4,                             & ! 2 mom ncn
                    &     p_reff        = prm_diag%reff_qr(:,:,:) )      ! Output

      ! Snow
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 3,                             & ! Snow DSD and geom. properties
                    &     grid_scope    = 0,                             & ! Grid and subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqs),      & ! Snow mixing ratio
                    &     p_ncn3D       = p_prog%tracer(:,:,:,iqns),     & ! Snow Number concentration
                    &     reff_param    = 1,                             & ! Fu
                    &     ncn_param     = 4,                             & ! 2 mom ncn
                    &     p_reff        = prm_diag%reff_qs(:,:,:) )        ! Output

      ! Graupel
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 4,                             & ! Graupel DSD and geom. properties
                    &     grid_scope    = 0 ,                            & ! Grid and subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqg),      & ! Graupel mixing ratio
                    &     p_ncn3D       = p_prog%tracer(:,:,:,iqng),     & ! Graupel Number concentration
                    &     reff_param    = 1 ,                            & ! Fu
                    &     ncn_param     = 4,                             & ! 2 mom ncn
                    &     p_reff        = prm_diag%reff_qg(:,:,:) )        ! Output
   

      ! Hail
      nreff_calc = nreff_calc + 1
      CALL  init_reff_calc ( reff_calc_dom(jg)%reff_calc_arr(nreff_calc),&
                    &     hydrometeor   = 5,                             & ! Hail DSD and geom. properties
                    &     grid_scope    = 0 ,                            & ! Grid and subgrid
                    &     microph_param = icalc_reff_loc,                & ! Chosen 2 moment scheme
                    &     return_fct    = return_fct,                    & ! Return parameter
                    &     p_q           = p_prog%tracer(:,:,:,iqh),      & ! Hail mixing ratio
                    &     p_ncn3D       = p_prog%tracer(:,:,:,iqnh),     & ! Hail Number concentration
                    &     reff_param    = 1 ,                            & ! Fu
                    &     ncn_param     = 4,                             & ! 2 mom ncn
                    &     p_reff        = prm_diag%reff_qh(:,:,:) )        ! Output


    END SELECT

    reff_calc_dom(jg)%nreff_calc = nreff_calc

    IF ( return_fct ) THEN
      IF ( PRESENT (return_reff) ) return_reff = .true.
      IF (msg_level >= 15) THEN 
        WRITE (message_text,*) "Reff: init_reff finished"
        CALL message('',message_text)
      END IF
    ELSE
      IF ( PRESENT (return_reff) ) return_reff = .false.
      WRITE (message_text,*) 'WARNING Reff: Something went wrong in the initialization.'
      CALL message('',message_text)
    END IF
  
  END SUBROUTINE init_reff
  ! End initiation


! Main routine to calculate effective radius
! Calculate effective radius according to different parameterizations
  SUBROUTINE set_reff( prm_diag, p_patch, p_prog, p_nh_diag, ext_data, p_metrics, return_reff )

    TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag         !Diagnostic statistics from physics (inlcuding reff)
    TYPE(t_patch)          , INTENT(in)   :: p_patch          !<grid/patch info.
    TYPE(t_nh_prog)        , INTENT(in)   :: p_prog           !<the dyn prog vars
    TYPE(t_nh_diag)        , INTENT(in)   :: p_nh_diag        !<the dyn prog vars
    TYPE(t_external_data)  , INTENT(in)   :: ext_data         ! External data (like land fraction for RRTM)
    TYPE(t_nh_metrics)     , INTENT(in), OPTIONAL    :: p_metrics
    LOGICAL, OPTIONAL      , INTENT(out)  :: return_reff      ! Return call from function .true. if right

    ! End of subroutine variables

    ! Auxiliary arrays
    REAL(wp)              :: ncn(nproma,p_patch%nlev)     ! Number concentration
    INTEGER (KIND=i4)     :: indices(nproma,p_patch%nlev) ! Mapping for going through array
    INTEGER (KIND=i4)     :: n_ind(p_patch%nlev)          ! Number of indices for each k level

    ! Local scalars:
    INTEGER               :: jb, jg, ireff                !<block indices
    INTEGER               :: i_startblk, i_endblk         !< blocks
    INTEGER               :: is, ie                       !< slices
    INTEGER               :: i_rlstart, i_rlend           ! blocks limits 
    INTEGER               :: nlev                         ! Number of grid levels in vertical
    INTEGER               :: jc, jk                       ! loop indices
    LOGICAL               :: return_fct(nreff_max_calc)   ! Return values
    INTEGER               :: nreff_calc                   ! Number of effective radius calculations 

    ! Pointers for IN parameterization
    REAL(wp), POINTER, DIMENSION(:,:) :: ptr_cams5=>NULL(), ptr_cams6=>NULL()
    REAL(wp), POINTER, DIMENSION(:)   :: ptr_aer_dust=>NULL()

    ! Initiate proper return
    IF ( PRESENT (return_reff) ) return_reff = .true.


    ! domain ID
    jg         = p_patch%id
    ! Number of levels
    nlev       = p_patch%nlev

  ! Number of reff calculations
    nreff_calc = reff_calc_dom(jg)%nreff_calc

! Prepare the loop over grid points
! exclude boundary interpolation zone of nested domains
    i_rlstart  = grf_bdywidth_c+1
    i_rlend    = min_rlcell_int

    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

    return_fct = .true.   ! Start with right return

    ! Main loop HERE OMP

    !$ACC ENTER DATA CREATE(indices, n_ind, ncn)
    !$ACC DATA PRESENT(reff_calc_dom(jg:jg), reff_calc_dom(jg)%reff_calc_arr(nreff_calc)) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(1)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(2)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(3)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(4)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(5)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(6)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(7)%p_reff) &
    !$ACC   PRESENT(reff_calc_dom(jg)%reff_calc_arr(8)%p_reff)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,is,ie,ncn,indices,n_ind,ireff,ptr_cams5,ptr_cams6, ptr_aer_dust) FIRSTPRIVATE(return_fct) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk     

      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &                is, ie, i_rlstart, i_rlend)

! Set to zero all
      !$ACC PARALLEL DEFAULT(NONE) ASYNC(1) FIRSTPRIVATE(nreff_calc, jg, jb)
      !$ACC LOOP GANG VECTOR
      DO ireff = 1, nreff_calc
        reff_calc_dom(jg)%reff_calc_arr(ireff)%p_reff(:,:,jb) = 0.0_wp  ! Clean values
      END DO
      !$ACC END PARALLEL

      DO ireff = 1, nreff_calc

      ! Obtain indices using qc_dia
        CALL mapping_indices (  indices, n_ind,                                             &
             &                  reff_calc=reff_calc_dom(jg)%reff_calc_arr(ireff),           &
             &                  k_start=kstart_moist(jg), k_end=nlev, is=is, ie=ie, jb=jb,  &
             &                  return_fct=return_fct(ireff) )

        IF ( .NOT. return_fct(ireff) ) THEN
          WRITE(*,*) "WARNING: Something went wrong with mapping_indices for calculation ", &
             & ireff," in domain ", jg
          IF ( PRESENT (return_reff) ) return_reff = .false.
        END IF

        IF ( icpl_aero_ice == 1) THEN
          SELECT CASE(irad_aero)
            CASE (iRadAeroCAMStd, iRadAeroCAMSclim)
              ptr_cams5=>p_nh_diag%camsaermr(:,:,jb,5)
              ptr_cams6=>p_nh_diag%camsaermr(:,:,jb,6)
            CASE (iRadAeroTegen)
              ptr_aer_dust=>prm_diag%aerosol(:,idu,jb)
            CASE DEFAULT
              CALL finish('mo_nwp_reff_interface', 'icpl_aero_ice = 1 only available for irad_aero = 6,7,8.')
            END SELECT
        ENDIF

        CALL calculate_ncn (    ncn, reff_calc_dom(jg)%reff_calc_arr(ireff) ,indices,       &
             &                  n_ind, kstart_moist(jg), nlev, jb,                          &
             &                  p_prog%rho(:,:,jb),  p_nh_diag%temp(:,:,jb),                &
             &                  icpl_aero_ice,                                              &
             &                  p_metrics%z_ifc(:,:,jb),                                    &
             &                  ptr_cams5, ptr_cams6,                                       &
             &                  ptr_aer_dust,                                               &
             &                  return_fct=return_fct(ireff))

        IF ( .NOT. return_fct(ireff) )  THEN
          WRITE(*,*) "WARNING: Something went wrong with calculate_ncn", ireff,             &
             &       " in domain ", jg      
          IF ( PRESENT (return_reff) ) return_reff = .false.
        END IF

        CALL calculate_reff (   reff_calc_dom(jg)%reff_calc_arr(ireff), indices, n_ind,     &
             &                  rho=p_prog%rho(:,:,jb),                                     &
             &                  ncn=ncn, clc= prm_diag%clc(:,:,jb),                         &
             &                  k_start=kstart_moist(jg), k_end=nlev, jb=jb,                &
             &                  fr_gl = ext_data%atm%fr_glac(:,jb),                         &
             &                  fr_land = ext_data%atm%fr_land(:,jb),                       &
             &                  return_fct=return_fct(ireff) )

        IF ( .NOT. return_fct(ireff) )  THEN 
          WRITE(*,*) "WARNING: Something went wrong with calculate_reff for calculation ",  &
             &      ireff ," in domain ", jg      
          IF ( PRESENT (return_reff) ) return_reff = .false.        
        END IF

      END DO

    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  
  IF (msg_level >= 15) THEN 
    WRITE (message_text,*) 'Reff calculated for domain ', jg
    CALL message('',message_text)
  END IF
  
  IF ( PRESENT (return_reff) ) THEN
    IF ( .NOT. return_reff ) THEN
      WRITE (message_text,*) 'WARNING: Something wrong in Reff calculations for domain ', jg
      CALL message('',message_text)
    END IF
  END IF
  
  !$ACC WAIT
  !$ACC END DATA ! p_reff
  !$ACC EXIT DATA DELETE(indices, n_ind, ncn)

  END SUBROUTINE set_reff

! Main routine to combine all hidropmeteors into a frozen phase and a liquid phase for radiation
! It changes the effective radius fields, as well as the tot_cld fields
  SUBROUTINE combine_phases_radiation_reff( prm_diag, p_patch, p_prog, return_reff )
      TYPE(t_nwp_phy_diag)   , INTENT(inout):: prm_diag         !Diagnostic statistics from physics (inlcuding reff)
      TYPE(t_patch)          , INTENT(in)   :: p_patch          !<grid/patch info.
      TYPE(t_nh_prog)        , INTENT(in)   :: p_prog           !<the dyn prog vars
      LOGICAL, OPTIONAL      , INTENT(out)  :: return_reff      ! Return call from function .true. if right

      ! Limitors of effective radius.
      REAL(wp), PARAMETER   :: reff_max_qc = 30.0e-6_wp      ! Maximum reff cloud water in tables from RRTM
      REAL(wp), PARAMETER   :: reff_max_qi = 99.0e-6_wp      ! Maximum reff ice in ECCRAD. It Crashes for larger values
!      REAL(wp), PARAMETER   :: reff_max_qi = 120.0e-6       ! Maximum reff ice in tables from RRTM

      ! Local scalars:
      INTEGER               :: jc,jk,jb,jg, ireff           !<block indices
      INTEGER               :: i_startblk, i_endblk         !< blocks
      INTEGER               :: is, ie                       !< slices
      INTEGER               :: i_rlstart, i_rlend           ! blocks limits
      INTEGER               :: nlev                         ! Number of grid levels in vertical

      !

      ! Initiare proper return
      IF ( PRESENT (return_reff) ) return_reff = .true.


      ! domain ID
      jg         = p_patch%id
      ! Number of levels
      nlev       = p_patch%nlev

      ! Prepare the loop over grid points
      ! exclude boundary interpolation zone of nested domains
      i_rlstart  = grf_bdywidth_c+1
      i_rlend    = min_rlcell_int

      i_startblk = p_patch%cells%start_block(i_rlstart)
      i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jk,jc,is,ie) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
           &                is, ie, i_rlstart, i_rlend)

      ! .. Branch clc_rad from the normal clc and modify it for input to the radiation scheme:
      
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO jk = kstart_moist(jg), nlev
        DO jc = is, ie
          prm_diag%clc_rad(jc,jk,jb) = prm_diag%clc(jc,jk,jb)
        END DO
      END DO
      !$ACC END PARALLEL
      
! Combine rain into the liquid phase and modify clc_rad accordingly:
      IF ( ASSOCIATED( prm_diag%reff_qr ) ) THEN
        CALL combine_reff( prm_diag%tot_cld(:,:,jb,iqc), prm_diag%reff_qc(:,:,jb), &
          &                p_prog%tracer(:,:,jb,iqr), prm_diag%reff_qr(:,:,jb),    &
          &                prm_diag%clc_rad(:,:,jb), k_start=kstart_moist(jg),     &
          &                k_end=nlev, is=is, ie=ie                                )
      END IF

! Combine snow into the ice phase and modify clc_rad accordingly:
      IF ( ASSOCIATED( prm_diag%reff_qs ) ) THEN
        CALL combine_reff( prm_diag%tot_cld(:,:,jb,iqi), prm_diag%reff_qi(:,:,jb), &
          &                p_prog%tracer(:,:,jb,iqs), prm_diag%reff_qs(:,:,jb),    &
          &                prm_diag%clc_rad(:,:,jb), k_start=kstart_moist(jg),     &
          &                k_end=nlev, is=is, ie=ie                                )
      END IF
! Combine graupel into the ice phase and modify clc_rad accordingly:
      IF ( (iqg > 0 ) .AND. ASSOCIATED( prm_diag%reff_qg ) ) THEN
        CALL combine_reff( prm_diag%tot_cld(:,:,jb,iqi), prm_diag%reff_qi(:,:,jb), &
          &                p_prog%tracer(:,:,jb,iqg), prm_diag%reff_qg(:,:,jb),    &
          &                prm_diag%clc_rad(:,:,jb), k_start=kstart_moist(jg),     &
          &                k_end=nlev, is=is, ie=ie                                )
      END IF
! Set maximum radius in the liquid phase keeping the ratio q/r constant
       CALL set_max_reff ( prm_diag%tot_cld(:,:,jb,iqc), prm_diag%reff_qc(:,:,jb), reff_max_qc, &
           &               k_start=kstart_moist(jg), k_end=nlev, is=is, ie=ie )

! Set maximum radius in the frozen phase keeping the ratio q/r constant
       CALL set_max_reff ( prm_diag%tot_cld(:,:,jb,iqi), prm_diag%reff_qi(:,:,jb), reff_max_qi, &
           &               k_start=kstart_moist(jg), k_end=nlev, is=is, ie=ie )


    END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  END SUBROUTINE combine_phases_radiation_reff


END MODULE mo_nwp_reff_interface

