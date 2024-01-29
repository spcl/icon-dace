!  Namelist for configuration of 2-moment cloud microphysics scheme
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

MODULE mo_2mom_mcrph_nml

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_io_units,            ONLY: nnml, nnml_output, filename_max
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,    &
    &                               open_and_restore_namelist, close_tmpfile

  USE mo_atm_phy_nwp_config,  ONLY: atm_phy_nwp_config

  USE mo_2mom_mcrph_config_default,ONLY: cfg_2mom_default

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_2mom_mcrph_namelist

  ! module name
  CHARACTER(*), PARAMETER :: modname = "mo_2mom_mcrph_nml"
  
CONTAINS

  !-------------------------------------------------------------------------
  !
  !! Read Namelist for the 2-moment cloud microphysics. 
  !!
  !! This subroutine 
  !! - reads the Namelist for NWP physics
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - performs sanity checks
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)    
  !!
  SUBROUTINE read_2mom_mcrph_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit, jg
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER ::  &
         &  routine = modname//':read_2mom_mcrph_namelist'

    !-------------------------------------------------------------------------
    ! Namelist variables
    !-------------------------------------------------------------------------

    ! .. The namelist parameters for the config parameters in the container:
    ! (Not yet domain dependent; same config for all domains)
    INTEGER  :: i2mom_solver ! 0) explicit (1) semi-implicit solver
    INTEGER  :: ccn_type     ! if not set by namelist, the ccn_type_gscp4 or ccn_type_gscp5 will win
    REAL(wp) :: alpha_spacefilling  !..factor involved in the conversion of ice/snow to graupel by riming
    REAL(wp) :: D_conv_ii    ! D-threshold for conversion to snow ice_selfcollection 
    REAL(wp) :: D_rainfrz_ig ! rain --> ice oder graupel
    REAL(wp) :: D_rainfrz_gh ! rain --> graupel oder hail
    LOGICAL  :: luse_mu_Dm_rain ! Use mu-Dm-Relation of Seifert (2008). If false use the constant in rain type
    REAL(wp) :: nu_r         ! nu for rain, N(x) = N0 * D^nu * exp(-lambda*x^mu)
    REAL(wp) :: rain_cmu0    ! asymptotic mue-value for small D_m in the mu-Dm-Relation of Seifert (2008)
    REAL(wp) :: rain_cmu1    ! asymptotic mue-value for large D_m in the mu-Dm-Relation of Seifert (2008)
    REAL(wp) :: rain_cmu3    ! D_br: equilibrium diameter for breakup and selfcollection
    REAL(wp) :: rain_cmu4    ! mue-value at D_br in the mu-Dm-Relation of Seifert (2008) 
    REAL(wp) :: melt_h_tune_fak ! Factor to increase/decrease hail melting rate of hail
    REAL(wp) :: Tmax_gr_rime    ! Allow formation of graupel by riming ice/snow only at T < this threshold [K]
    LOGICAL  :: lturb_enhc   ! Enhancesment of collisons by turbulence (only warm microphysics)
    REAL(wp) :: ecoll_gg        ! Collision efficiency for graupel autoconversion (dry graupel)
    REAL(wp) :: ecoll_gg_wet    ! Collision efficiency for graupel autoconversion (wet graupel)
    REAL(wp) :: Tcoll_gg_wet ! Temperature threshold for switching to wet graupel autoconversion
    REAL(wp) :: melt_g_tune_fak     ! Factor multiplying melting of graupel
    INTEGER  :: iice_stick     ! Formulation for sticking efficiency of cloud ice
    INTEGER  :: isnow_stick     ! Formulation for sticking efficiency of snow/graupel
    INTEGER  :: iparti_stick     ! Formulation for sticking efficiency of frozen inter-categorical collisions
    REAL(wp) :: nu_h      ! nu for hail, N(x) = N0 * D^nu * exp(-lambda*x^mu)
    REAL(wp) :: mu_h      ! mu for hail, N(x) = N0 * D^nu * exp(-lambda*x^mu)
    REAL(wp) :: ageo_h    ! ageo for hail, D = ageo*x^bgeo
    REAL(wp) :: bgeo_h    ! bgeo for hail, D = ageo*x^bgeo
    REAL(wp) :: avel_h    ! avel for hail, v = avel*x^bvel
    REAL(wp) :: bvel_h    ! bvel for hail, v = avel*x^bvel
    REAL(wp) :: cap_snow  ! capacitance for snow deposition/sublimation
    REAL(wp) :: vsedi_max_s ! max fallspeed limit for snow
    REAL(wp) :: nu_g      ! nu for graupel, N(x) = N0 * D^nu * exp(-lambda*x^mu)
    REAL(wp) :: mu_g      ! mu for graupel, N(x) = N0 * D^nu * exp(-lambda*x^mu)
    REAL(wp) :: ageo_g    ! ageo for graupel, D = ageo*x^bgeo
    REAL(wp) :: bgeo_g    ! bgeo for graupel, D = ageo*x^bgeo
    REAL(wp) :: avel_g    ! avel for graupel, v = avel*x^bvel
    REAL(wp) :: bvel_g    ! bvel for graupel, v = avel*x^bvel
    INTEGER  :: iicephase    ! (0) warm-phase 2M, (1) mixed-phase 2M

    NAMELIST /twomom_mcrph_nml/ i2mom_solver, ccn_type, alpha_spacefilling,             &
         &                      D_conv_ii, D_rainfrz_ig, D_rainfrz_gh,                  &
         &                      luse_mu_Dm_rain, nu_r, rain_cmu0, rain_cmu1, rain_cmu3, &
         &                      rain_cmu4, melt_h_tune_fak, Tmax_gr_rime, lturb_enhc,   &
         &                      ecoll_gg, ecoll_gg_wet, Tcoll_gg_wet,                   &
         &                      melt_g_tune_fak, isnow_stick, iice_stick, iparti_stick, &
         &                      nu_h, mu_h, ageo_h, bgeo_h, avel_h, bvel_h, cap_snow,   &
         &                      vsedi_max_s,                &
         &                      nu_g, mu_g, ageo_g, bgeo_g, avel_g, bvel_g,             &
         &                      iicephase 

    !----------------------------------------------------------
    ! 1. default settings from module mo_2mom_mcrph_processes:
    !----------------------------------------------------------

    ! .. Initialize the container with the defaults from mo_2mom_mcrph_config_default:
    DO jg=1,max_dom
      atm_phy_nwp_config(jg) % cfg_2mom = cfg_2mom_default
    END DO

    ! .. Initialize the namelist parameters which are later put back into the container:
    i2mom_solver       = cfg_2mom_default % i2mom_solver
    ccn_type           = cfg_2mom_default % ccn_type
    alpha_spacefilling = cfg_2mom_default % alpha_spacefilling
    D_conv_ii          = cfg_2mom_default % D_conv_ii
    D_rainfrz_ig       = cfg_2mom_default % D_rainfrz_ig
    D_rainfrz_gh       = cfg_2mom_default % D_rainfrz_gh
    luse_mu_Dm_rain    = cfg_2mom_default % luse_mu_Dm_rain
    nu_r               = cfg_2mom_default % nu_r
    rain_cmu0          = cfg_2mom_default % rain_cmu0     
    rain_cmu1          = cfg_2mom_default % rain_cmu1
    rain_cmu3          = cfg_2mom_default % rain_cmu3
    rain_cmu4          = cfg_2mom_default % rain_cmu4
    melt_g_tune_fak    = cfg_2mom_default % melt_g_tune_fak     
    melt_h_tune_fak    = cfg_2mom_default % melt_h_tune_fak     
    Tmax_gr_rime       = cfg_2mom_default % Tmax_gr_rime        
    lturb_enhc         = cfg_2mom_default % lturb_enhc
    ecoll_gg           = cfg_2mom_default % ecoll_gg
    ecoll_gg_wet       = cfg_2mom_default % ecoll_gg_wet
    Tcoll_gg_wet       = cfg_2mom_default % Tcoll_gg_wet
    iice_stick         = cfg_2mom_default % iice_stick 
    isnow_stick        = cfg_2mom_default % isnow_stick 
    iparti_stick       = cfg_2mom_default % iparti_stick 
    nu_h               = cfg_2mom_default % nu_h
    mu_h               = cfg_2mom_default % mu_h
    ageo_h             = cfg_2mom_default % ageo_h
    bgeo_h             = cfg_2mom_default % bgeo_h
    avel_h             = cfg_2mom_default % avel_h
    bvel_h             = cfg_2mom_default % bvel_h
    cap_snow           = cfg_2mom_default % cap_snow
    vsedi_max_s        = cfg_2mom_default % vsedi_max_s
    nu_g               = cfg_2mom_default % nu_g
    mu_g               = cfg_2mom_default % mu_g
    ageo_g             = cfg_2mom_default % ageo_g
    bgeo_g             = cfg_2mom_default % bgeo_g
    avel_g             = cfg_2mom_default % avel_g
    bvel_g             = cfg_2mom_default % bvel_g
    iicephase          = cfg_2mom_default % iicephase

    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, twomom_mcrph_nml)   ! write defaults to temporary text file
    END IF

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above 
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('twomom_mcrph_nml')
      READ(funit,NML=twomom_mcrph_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('twomom_mcrph_nml', status=istat)

    SELECT CASE (istat)
    CASE (POSITIONED)


      READ (nnml, twomom_mcrph_nml)   ! overwrite default settings

      ! SHOULD THERE BE DOMAIN DEPENDENCE IN THE FUTURE:
      ! Restore default values for global domain WHERE nothing at all has been specified
     
      ! Copy values of parent domain (in case of linear nesting) to nested domains where nothing has been specified

      ! Is currently not needed, because cfg_2mom should be the same for all domains, but here is a blueprint:
!!$      DO jg = 2, max_dom
!!$
!!$        ! Physics packages
!!$        IF (ccn_type(jg)       < 0) ccn_type(jg)       = ccn_type(jg-1)
!!$
!!$      ENDDO


      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, twomom_mcrph_nml)   ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !----------------------------------------------------
    ! 4. Sanity check
    !----------------------------------------------------
    
    ! check for valid parameters in namelists:

    IF (ALL(i2mom_solver /= (/0, 1/)) ) THEN
      CALL finish( TRIM(routine), 'Incorrect setting for cfg_2mom%i2mom_solver. Must be 0, or 1.')
    END IF

    IF (ALL(ccn_type /= (/-1, 6, 7, 8, 9/)) ) THEN
      CALL finish( TRIM(routine), 'Incorrect setting for cfg_2mom%ccn_type. Must be -1, 6, 7, 8, or 9.')
    END IF

    IF (ALL(iicephase /= (/0, 1/)) ) THEN
      CALL finish( TRIM(routine), 'Incorrect setting for cfg_2mom%iicephase. Must be 0, or 1.')
    END IF

    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg=1,max_dom

      atm_phy_nwp_config(jg) % cfg_2mom % i2mom_solver        = i2mom_solver
      atm_phy_nwp_config(jg) % cfg_2mom % ccn_type            = ccn_type
      atm_phy_nwp_config(jg) % cfg_2mom % alpha_spacefilling  = alpha_spacefilling
      atm_phy_nwp_config(jg) % cfg_2mom % D_conv_ii           = D_conv_ii
      atm_phy_nwp_config(jg) % cfg_2mom % D_rainfrz_ig        = D_rainfrz_ig
      atm_phy_nwp_config(jg) % cfg_2mom % D_rainfrz_gh        = D_rainfrz_gh
      atm_phy_nwp_config(jg) % cfg_2mom % luse_mu_Dm_rain     = luse_mu_Dm_rain
      atm_phy_nwp_config(jg) % cfg_2mom % nu_r                = nu_r
      atm_phy_nwp_config(jg) % cfg_2mom % rain_cmu0           = rain_cmu0
      atm_phy_nwp_config(jg) % cfg_2mom % rain_cmu1           = rain_cmu1
      atm_phy_nwp_config(jg) % cfg_2mom % rain_cmu3           = rain_cmu3
      atm_phy_nwp_config(jg) % cfg_2mom % rain_cmu4           = rain_cmu4
      atm_phy_nwp_config(jg) % cfg_2mom % melt_h_tune_fak     = melt_h_tune_fak     
      atm_phy_nwp_config(jg) % cfg_2mom % Tmax_gr_rime        = Tmax_gr_rime        
      atm_phy_nwp_config(jg) % cfg_2mom % lturb_enhc          = lturb_enhc
      atm_phy_nwp_config(jg) % cfg_2mom % ecoll_gg            = ecoll_gg
      atm_phy_nwp_config(jg) % cfg_2mom % ecoll_gg_wet        = ecoll_gg_wet
      atm_phy_nwp_config(jg) % cfg_2mom % Tcoll_gg_wet        = Tcoll_gg_wet
      atm_phy_nwp_config(jg) % cfg_2mom % melt_g_tune_fak     = melt_g_tune_fak
      atm_phy_nwp_config(jg) % cfg_2mom % iice_stick          = iice_stick    
      atm_phy_nwp_config(jg) % cfg_2mom % isnow_stick         = isnow_stick    
      atm_phy_nwp_config(jg) % cfg_2mom % iparti_stick        = iparti_stick    
      atm_phy_nwp_config(jg) % cfg_2mom % nu_h                = nu_h
      atm_phy_nwp_config(jg) % cfg_2mom % mu_h                = mu_h
      atm_phy_nwp_config(jg) % cfg_2mom % ageo_h              = ageo_h
      atm_phy_nwp_config(jg) % cfg_2mom % bgeo_h              = bgeo_h
      atm_phy_nwp_config(jg) % cfg_2mom % avel_h              = avel_h
      atm_phy_nwp_config(jg) % cfg_2mom % bvel_h              = bvel_h
      atm_phy_nwp_config(jg) % cfg_2mom % cap_snow            = cap_snow
      atm_phy_nwp_config(jg) % cfg_2mom % vsedi_max_s         = vsedi_max_s
      atm_phy_nwp_config(jg) % cfg_2mom % nu_g                = nu_g
      atm_phy_nwp_config(jg) % cfg_2mom % mu_g                = mu_g
      atm_phy_nwp_config(jg) % cfg_2mom % ageo_g              = ageo_g
      atm_phy_nwp_config(jg) % cfg_2mom % bgeo_g              = bgeo_g
      atm_phy_nwp_config(jg) % cfg_2mom % avel_g              = avel_g
      atm_phy_nwp_config(jg) % cfg_2mom % bvel_g              = bvel_g
      atm_phy_nwp_config(jg) % cfg_2mom % iicephase           = iicephase

    ENDDO

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=twomom_mcrph_nml)                    
      CALL store_and_close_namelist(funit, 'twomom_mcrph_nml') 
    ENDIF
    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=twomom_mcrph_nml)

  END SUBROUTINE read_2mom_mcrph_namelist

END MODULE mo_2mom_mcrph_nml

