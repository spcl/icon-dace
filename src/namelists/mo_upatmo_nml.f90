! This module controls the namelist input for:
! - Upper-atmosphere physics
! - Upper-atmosphere extrapolation
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

MODULE mo_upatmo_nml

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: finish
  USE mo_upatmo_config,           ONLY: upatmo_exp_config, upatmo_phy_config
  USE mo_namelist,                ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_io_units,                ONLY: nnml, nnml_output, filename_max
  USE mo_mpi,                     ONLY: my_process_is_stdio
  USE mo_master_control,          ONLY: use_restart_namelists
  USE mo_restart_nml_and_att,     ONLY: open_tmpfile, store_and_close_namelist, &
    &                                   open_and_restore_namelist, close_tmpfile
  USE mo_nml_annotate,            ONLY: temp_defaults, temp_settings
  USE mo_impl_constants,          ONLY: SUCCESS, MAX_CHAR_LENGTH, max_dom
  USE mo_util_string,             ONLY: int2string
  USE mo_upatmo_impl_const,       ONLY: isolvar, isolvardat, iorbit, icycle,              &
    &                                   iUpatmoPrcMode, iUpatmoPrcId, iUpatmoGasMode,     & 
    &                                   iUpatmoGrpId, iUpatmoGasId, iUpatmoExtdatId,      &
    &                                   iThermdynCoupling
  USE mo_physical_constants,      ONLY: amd, amco2, amo2, amo3, amo, amno
  USE mtime,                      ONLY: MAX_DATETIME_STR_LEN
  USE mo_upatmo_utils,            ONLY: isInInterval

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_upatmo_namelist

  CHARACTER(LEN = *), PARAMETER :: modname = 'mo_upatmo_nml'

  !------------------------------------------------------------
  !                    Namelist parameters
  !------------------------------------------------------------

  !---------------
  ! Extrapolation
  !---------------

  ! Variables for extrapolation methods
  !
  REAL(wp) :: expol_start_height     ! [m] Height above which extrapolation (blending) of initial data starts
  REAL(wp) :: expol_blending_scale   ! [m] Blending scale height
  REAL(wp) :: expol_vn_decay_scale   ! [m] Scale height for (exponential) decay of extrapolated 
                                     ! horizontal wind component (for stability reasons)
  REAL(wp) :: expol_temp_infty       ! [K] Climatological temperature of exosphere (z -> infinity)  
  LOGICAL  :: lexpol_sanitycheck     ! .TRUE. -> Apply sanity check to extrapolated fields

  !---------------
  !   Physics
  !---------------

  INTEGER  :: orbit_type   ! Orbit model: 
                           ! * 1: vsop87 -> standard and accurate model (vsop87)
                           ! * 2: kepler -> simple model, appropriate for idealized work 
  INTEGER  :: solvar_type  ! Solar activity: 
                           ! * 1: normal activity
                           ! * 2: low --,,--
                           ! * 3: high --,,--
  INTEGER  :: solvar_data  ! Solar activity data type:
                           ! * 1: G. Rottman data
                           ! * 2: J. Lean data
  INTEGER  :: solcyc_type  ! Solar cycle:
                           ! * 1: standard cycle
                           ! * 2: 27day cycle
  ! The following parameters 
  ! are not controllable by namelist
  ! for the time being
  REAL(wp) :: cecc         ! Eccentricity of orbit 
  REAL(wp) :: cobld        ! Obliquity of Earth axis
  REAL(wp) :: clonp        ! Longitude of perihelion 
  LOGICAL  :: lyr_perp     ! .TRUE.: Earth orbit of year 'yr_perp' 
                           ! of the VSOP87 orbit model is perpetuated
                           ! .FALSE.: transient Earth orbit following VSOP87
  INTEGER  :: yr_perp      ! Year used for 'lyr_perp = .TRUE.' 
  LOGICAL  :: lsanitycheck ! Switch for applying sanity checks

  !--------------------------------

  ! The following types are conceptual copies 
  ! from the process control of ECHAM and NWP

  ! Start heights, above which physical processes compute tendencies in case of ECHAM-physics

  TYPE t_aes_start_height_nml
    REAL(wp) :: all       ! For all processes
    ! 
    REAL(wp) :: rad       ! For processes of the group RAD (-> SRBC, EUV, NO, CHEMHEAT)
    REAL(wp) :: imf       ! For processes of the group IMF (-> IONDRAG, VDFMOL, FRIC, JOULE-HEATING)
    !
    REAL(wp) :: srbc      ! For heating due to Schumann-Runge bands and continuum of O2
    REAL(wp) :: euv       ! For extreme-ultraviolet heating
    ! NLTE is not modifiable in that regard
    REAL(wp) :: no        ! For near-infrared heating by NO
    REAL(wp) :: chemheat  ! For chemical heating
    REAL(wp) :: iondrag   ! For ion drag / Joule heating
    REAL(wp) :: vdfmol    ! For molecular diffusion
    REAL(wp) :: fric      ! For frictional heating    
  END TYPE t_aes_start_height_nml

  TYPE(t_aes_start_height_nml) :: aes_start_height

  ! Auxiliary parameters for evaluation of start heights
  INTEGER, PARAMETER :: iNonNegVal = 1

  !--------------------------------

  ! Control type for the groups in which the single upper-atmosphere physics processes are clustered, 
  ! following the example set by 'src/configre_model/mo_aes_phy_config' (only required for NWP-physics)  
  ! (Unfortunately, the compilers seem to be unable to cope with string arrays in derived types 
  ! dedicated for the namelist read-in, so for the time being start and end time apply to all domains.)

  TYPE t_nwp_prc_nml
    INTEGER                             :: imode(max_dom)  ! Select mode of process group
    REAL(wp)                            :: dt(max_dom)     ! Time step for process group
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: t_start         ! Start time of process group
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: t_end           ! End time of process group
    REAL(wp)                            :: start_height    ! Start height of process group
  END TYPE t_nwp_prc_nml

  TYPE(t_nwp_prc_nml) :: nwp_grp_imf  ! Ion drag (I), molecular diffusion (M) and frictional heating (F)
  TYPE(t_nwp_prc_nml) :: nwp_grp_rad  ! Radiation and chemical heating 

  !--------------------------------

  ! Namelist input required for the radiatively active gases following 
  ! * src/configure_model/mo_radiation_config
  ! * src/configure_model/mo_aes_rad_config
  ! Please note that the upper-atmosphere settings here are independent from the settings there!
  ! NWP physics and upper-atmosphere physics work with different gas climatologies.  
  ! A harmonization is not straightforward and postponed to the future. 
  ! (No domain-dependence (-> 'max_dom') envisaged for the time being.)

  TYPE t_nwp_gas_nml
    INTEGER  :: imode    ! Gas mode (corresponds to 'irad' in 'mo_radiation_config')
    REAL(wp) :: vmr      ! Volume mixing ratio ((m3/m3), should be equal to mole fraction (mol/mol))
    REAL(wp) :: fscale   ! scaling factor for mixing ratio ('frad_<gas>' in 'mo_aes_rad_config')
  END TYPE t_nwp_gas_nml

  TYPE(t_nwp_gas_nml) :: nwp_gas_o3  ! Ozone
  TYPE(t_nwp_gas_nml) :: nwp_gas_o2  ! Dioxygen
  TYPE(t_nwp_gas_nml) :: nwp_gas_o   ! Atomic oxygen
  TYPE(t_nwp_gas_nml) :: nwp_gas_co2 ! Carbon dioxide
  TYPE(t_nwp_gas_nml) :: nwp_gas_no  ! Nitric oxide

  ! Dinitrogen (n2) is currently designated for the diagnostic mode 
  ! => no namelist entry required

  !--------------------------------

  ! Type for external data.
  ! (No domain-dependence envisaged for the time being.)

  TYPE t_nwp_extdat_nml
    REAL(wp)                    :: dt        ! Update period for time interpolation
    CHARACTER(LEN=filename_max) :: filename  ! Name of the file containing the external data
  END TYPE t_nwp_extdat_nml

  TYPE(t_nwp_extdat_nml) :: nwp_extdat_gases
  TYPE(t_nwp_extdat_nml) :: nwp_extdat_chemheat

  !--------------------------------

  ! Miscellaneous switches. 
  ! (These are "unofficial" namelist switches, 
  ! i.e. they do not and shall not appear in Namelist_overview.pdf.)
  INTEGER :: nwp_thermdyn_cpl         ! Type of thermodynamic coupling of physics & dynamics:
                                      ! * 1: isobaric coupling
                                      ! * 2: isochoric coupling
                                      ! * 3: entropic coupling
  LOGICAL :: nwp_ldiss_from_heatdiff  ! Switch for considering heat source from heat diffusion

  !------------------------------------------------------------

  NAMELIST /upatmo_nml/ expol_start_height,     & 
    &                   expol_blending_scale,   &
    &                   expol_vn_decay_scale,   &
    &                   expol_temp_infty,       &
    &                   lexpol_sanitycheck,     &
    &                   orbit_type,             &
    &                   solvar_type,            &
    &                   solvar_data,            &
    &                   solcyc_type,            &
    &                   aes_start_height,       &
    &                   nwp_grp_imf,            &
    &                   nwp_grp_rad,            &
    &                   nwp_gas_o3,             &
    &                   nwp_gas_o2,             &
    &                   nwp_gas_o,              &
    &                   nwp_gas_co2,            &
    &                   nwp_gas_no,             &
    &                   nwp_extdat_gases,       &
    &                   nwp_extdat_chemheat,    &
    ! "unofficial" switches
    &                   nwp_thermdyn_cpl,       &
    &                   nwp_ldiss_from_heatdiff

  !------------------------------------------------------------

CONTAINS !..................................................................................

  !------------------------------------------------------------
  !                       Subroutines
  !------------------------------------------------------------

  !>
  !! Read Namelist for upper atmosphere. 
  !!
  !! This subroutine: 
  !! - Reads the Namelist for upper atmosphere
  !! - Sets default values
  !! - Potentially overwrites the defaults by values used in a 
  !!   previous integration (if this is a resumed run)
  !! - Reads the user's (new) specifications
  !! - Does a consistency check
  !! - Stores the Namelist for restart
  !! - Fills the configuration state (partly) 
  !!
  SUBROUTINE read_upatmo_namelist( filename )

    ! In/out variables
    CHARACTER(LEN=*), INTENT(IN) :: filename

    ! Local variables
    REAL(wp) :: startheightall, startheightimf, startheightrad
    REAL(wp) :: sum_mmr
    INTEGER  :: istat, funit
    INTEGER  :: iunit
    INTEGER  :: jg, igrp, igas, iext, istartitem, ienditem
    LOGICAL  :: lvalid
    REAL(wp), PARAMETER :: ramd = 1._wp / amd
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':read_upatmo_namelist'

    !------------------------------------------------

    !------------------------------------------------------------
    !                    Set default values
    !------------------------------------------------------------

    ! PLEASE, remember to update doc/icon_nml_tables.tex 
    ! and run/upatmo_nml_elements, if you would change the defaults, 
    ! add new variables etc. Thank you! 

    ! (Where not otherwise stated, the settings apply to all domains)

    !---------------
    ! Extrapolation
    !---------------

    ! (Note: to switch on the upper-atmosphere extrapolation set 'initicon_nml: itype_vert_expol = 2'.)

    expol_start_height   = 70000.0_wp  ! (m) Height above which extrapolation (blending) of initial data starts
    expol_blending_scale = 10000.0_wp  ! (m) Blending scale height
    expol_vn_decay_scale = 10000.0_wp  ! (m) Horizontal wind decay scale height
    expol_temp_infty     = 400.0_wp    ! (K) Exospheric mean reference temperature 
                                       ! (e.g. climatological value from Hedin (1983) would be 1035 K)
    lexpol_sanitycheck   = .FALSE.     ! No sanity check of extrapolated fields 
                                       ! (bacause it is computationally expensive)

    !---------------
    !    Physics
    !---------------  

    orbit_type   = iorbit%vsop87    ! 1 -> Standard orbit model
    solvar_type  = isolvar%norm     ! 1 -> Normal solar activity
    solvar_data  = isolvardat%lean  ! 2 -> J. Lean data
    solcyc_type  = icycle%day27     ! 2 -> 27-day cycle
    ! NO namelist input:
    cecc         = 0.016715_wp      ! Default from 'src/configure_model/mo_aes_rad_config: init_aes_rad_config'
    cobld        = 23.44100_wp      ! --,,--
    clonp        = 282.7000_wp      ! --,,--
    lyr_perp     = .FALSE.          ! --,,--
    yr_perp      = -99999           ! --,,--
    lsanitycheck = .FALSE.          ! .TRUE.: apply sanity checks (e.g.: Is mass of diagnostic gas N2 >=0 ?). 
                                    ! Switch on with care, these checks come along with significant additional 
                                    ! computational costs!

    ! ECHAM-specific settings:

    ! Start heights
    ! -------------
    ! Negative values mean that the default start heights, set in
    ! 'src/upper_atmosphere/mo_upatmo_phy_config: init_upatmo_phy_config'
    ! will be used
    aes_start_height%all      = -999._wp  ! (m) For all processes
    !
    aes_start_height%rad      = -999._wp  ! (m) For processes of the group RAD 
                                            ! (-> SRBC, EUV, NO, CHEMHEAT)
    aes_start_height%imf      = -999._wp  ! (m) For processes of the group IMF 
                                            ! (-> IONDRAG, VDFMOL, FRIC, JOULE-HEATING)    
    ! For each process separately:
    aes_start_height%srbc     = -999._wp  ! (m) Heating due to Schumann-Runge bands and continuum of O2
    aes_start_height%euv      = -999._wp  ! (m) Extreme-ultraviolet heating
    aes_start_height%no       = -999._wp  ! (m) Near-infrared heating by NO
    aes_start_height%chemheat = -999._wp  ! (m) Chemical heating
    aes_start_height%iondrag  = -999._wp  ! (m) Ion drag / Joule heating
    aes_start_height%vdfmol   = -999._wp  ! (m) Molecular diffusion
    aes_start_height%fric     = -999._wp  ! (m) Frictional heating    

    ! NWP-specific settings:

    ! Please note: 
    ! * the corresponding settings for the ECHAM mode 
    !   will be added in a second step
    ! * to switch on the NWP + upatmo physics set 
    !   'nwp_phy_nml: lupatmo_phy = .TRUE.'
    ! * if 'lupatmo_phy(1) = .TRUE.', the below settings for 'imode(1)'
    !   will result in ALL processes being switched on by default 
    !   on the primary domain, and on all secondary domains as well, 
    !   if they are not explicitly switched off by 'lupatmo_phy(>1) = .FALSE.'
    !   or 'imode(>1) = iUpatmoPrcMode%off = 0'!
    ! * be careful, if you touch the defaults below, 
    !   some of them are necessary to determine the settings 
    !   for the nests following the example of several other namelists, 
    !   e.g., 'src/namelists/mo_nwp_phy_nml'

    ! Process group: ion drag (I), molecular diffusion (M) and frictional heating (F):
    ! --------------------------------------------------------------------------------
    nwp_grp_imf%imode(:)     = iUpatmoPrcMode%unassigned  ! -1 -> This means that the group will be switched on 
                                                          ! on nests, if the following conditions are met: 
                                                          ! * lupatmo_phy(>1) = .TRUE., 
                                                          !   which follows automatically, if the entry
                                                          !   in the namelist file would be: 
                                                          !      &nwp_phy_nml
                                                          !        lupatmo_phy = .TRUE.
                                                          !      /
                                                          ! * imode(1) = iUpatmoPrcMode%on = 1, 
                                                          !   which is the subsequent default setting
    nwp_grp_imf%imode(1)     = iUpatmoPrcMode%on          ! 1 -> Switched on by default, if 'lupatmo_phy(1) = .TRUE.'
    !
    nwp_grp_imf%dt(:)        = -999._wp
    nwp_grp_imf%dt(1)        =  300._wp  ! (s) Default value for tendency update period.
                                         ! The value for 'dt(1)' sets the upper bound for 'dt(>1)'!
    !
    nwp_grp_imf%t_start      = ' '       ! Start time is simulation start time
    nwp_grp_imf%t_end        = ' '       ! End time is simulation end time
    ! 
    nwp_grp_imf%start_height = -999._wp  ! (m) Negative value: take default start heights of single processes

    ! Process group: radiation and chemical heating:
    ! ----------------------------------------------
    nwp_grp_rad%imode(:)     = iUpatmoPrcMode%unassigned
    nwp_grp_rad%imode(1)     = iUpatmoPrcMode%on
    !
    nwp_grp_rad%dt(:)        = -999._wp
    nwp_grp_rad%dt(1)        =  600._wp 
    !
    nwp_grp_rad%t_start      = ' '
    nwp_grp_rad%t_end        = ' '
    ! 
    nwp_grp_rad%start_height = -999._wp
    
    ! Radiatively active gas: ozone 
    ! -----------------------------
    nwp_gas_o3%imode   = iUpatmoGasMode%extdat  ! 2 -> Gas concentration is read from file (=> vmr = 0)
    nwp_gas_o3%vmr     = 0._wp
    nwp_gas_o3%fscale  = 1._wp                  ! No rescaling of the entire gas concentration

    ! Radiatively active gas: dioxygen 
    ! --------------------------------
    nwp_gas_o2%imode   = iUpatmoGasMode%extdat
    nwp_gas_o2%vmr     = 0._wp
    nwp_gas_o2%fscale  = 1._wp   

    ! Radiatively active gas: atomic oxygen
    ! -------------------------------------
    nwp_gas_o%imode    = iUpatmoGasMode%extdat
    nwp_gas_o%vmr      = 0._wp
    nwp_gas_o%fscale   = 1._wp

    ! Radiatively active gas: carbon dioxide
    ! --------------------------------------
    nwp_gas_co2%imode  = iUpatmoGasMode%extdat
    nwp_gas_co2%vmr    = 0._wp
    nwp_gas_co2%fscale = 1._wp

    ! Radiatively active gas: nitric oxide
    ! ------------------------------------
    nwp_gas_no%imode   = iUpatmoGasMode%extdat
    nwp_gas_no%vmr     = 0._wp
    nwp_gas_no%fscale  = 1._wp

    ! External data: radiatively active gases
    ! ---------------------------------------
    nwp_extdat_gases%dt       = 86400._wp                   ! Update period for time interpolation 
                                                            ! of gas concentrations from external data: every day
    nwp_extdat_gases%filename = "upatmo_gases_chemheat.nc"  ! Name of file containing external data

    ! External data: chemical heating tendency
    ! ----------------------------------------
    nwp_extdat_chemheat%dt       = 86400._wp
    nwp_extdat_chemheat%filename = "upatmo_gases_chemheat.nc"

    ! "Unofficial" switches (to stay away from doc)
    ! ---------------------------------------------
    nwp_thermdyn_cpl        = iThermdynCoupling%isochoric
    nwp_ldiss_from_heatdiff = .FALSE.

    !------------------------------------------------------------
    !  If this is a resumed integration, overwrite the defaults 
    !     above by values used in the previous integration
    !------------------------------------------------------------
    
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('upatmo_nml')
      READ(funit,NML=upatmo_nml)
      CALL close_tmpfile(funit)
    END IF

    !------------------------------------------------------------
    !             Read user's (new) specifications 
    !            (Done so far by all MPI processes)
    !------------------------------------------------------------

    CALL open_nml(TRIM(filename))
    CALL position_nml ('upatmo_nml', status=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, upatmo_nml)    ! Write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, upatmo_nml)     ! Overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, upatmo_nml)  ! Write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml

    !------------------------------------------------------------
    !                     Consistency checks
    !------------------------------------------------------------

    !---------------
    ! Extrapolation
    !---------------

    IF( expol_start_height < 0._wp ) THEN 
      ! The height above which the extrapolation starts
      ! has to be equal to or greater than zero
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_start_height, it has to be >= 0.' )
    ELSEIF( expol_blending_scale < 0._wp ) THEN 
      ! The blending scale height has to be equal to 
      ! or greater than zero
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_blending_scale, it has to be >= 0.' )
    ELSEIF( expol_vn_decay_scale < 0._wp ) THEN 
      ! Likewise the decay scale height for the horizontal wind
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_vn_decay_scale, it has to be >= 0.' )
    ELSEIF( expol_temp_infty < 0._wp ) THEN 
      ! The exospheric mean temperature has to be equal to 
      ! or greater than zero
      CALL finish( TRIM(routine), & 
        & 'Invalid value for expol_temp_infty, it has to be >= 0.' )
    ENDIF

    !---------------
    !    Physics
    !---------------

    ! Check 'orbit_type'
    ! (Valid entries range from 1 to 'iorbit%nitem')
    lvalid = isInInterval(number=orbit_type, opt_clbnd=1, opt_cubnd=iorbit%nitem)
    IF (.NOT. lvalid) CALL finish(TRIM(routine), 'Invalid orbit_type: ' &
      & //TRIM(int2string(orbit_type)))

    ! Check 'solvar_type'
    lvalid = isInInterval(number=solvar_type, opt_clbnd=1, opt_cubnd=isolvar%nitem)
    IF (.NOT. lvalid) CALL finish(TRIM(routine), 'Invalid solvar_type: ' &
      & //TRIM(int2string(solvar_type)))

    ! Check 'solvar_data'
    lvalid = isInInterval(number=solvar_data, opt_clbnd=1, opt_cubnd=isolvardat%nitem)
    IF (.NOT. lvalid) CALL finish(TRIM(routine), 'Invalid solvar_data: ' &
      & //TRIM(int2string(solvar_data)))

    ! Check 'solcyc_type'
    lvalid = isInInterval(number=solcyc_type, opt_clbnd=1, opt_cubnd=icycle%nitem)
    IF (.NOT. lvalid) CALL finish(TRIM(routine), 'Invalid solcyc_type: ' &
      & //TRIM(int2string(solcyc_type)))

    ! Check mode of physics groups
    istartitem = iUpatmoPrcMode%startitem
    ienditem   = iUpatmoPrcMode%enditem
    DO jg = 1, max_dom
      ! IMF
      ! ---
      lvalid = isInInterval(number=nwp_grp_imf%imode(jg), opt_clbnd=istartitem, opt_cubnd=ienditem)
      ! 'unassigned' is no valid option for the primary domain
      IF (jg == 1) lvalid = lvalid .AND. .NOT. (nwp_grp_imf%imode(jg) == iUpatmoPrcMode%unassigned)
      IF (.NOT. lvalid) CALL finish( TRIM(routine), & 
        & 'Invalid value for nwp_grp_imf%imode('//TRIM(int2string(jg))//').' )
      ! RAD
      ! ---
      lvalid = isInInterval(number=nwp_grp_rad%imode(jg), opt_clbnd=istartitem, opt_cubnd=ienditem)
      IF (jg == 1) lvalid = lvalid .AND. .NOT. (nwp_grp_rad%imode(jg) == iUpatmoPrcMode%unassigned)
      IF (.NOT. lvalid) CALL finish( TRIM(routine), & 
        & 'Invalid value for nwp_grp_rad%imode('//TRIM(int2string(jg))//').' )      
    ENDDO  !jg

    ! Check tendency update period on primary domain
    ! IMF
    ! ---
    IF (.NOT. (nwp_grp_imf%dt(1) > 0._wp)) THEN
      CALL finish( TRIM(routine), 'Invalid value for nwp_grp_imf%dt(1) (dt > 0 required).' )
    ENDIF
    ! RAD
    ! ---
    IF (.NOT. (nwp_grp_rad%dt(1) > 0._wp)) THEN
      CALL finish( TRIM(routine), 'Invalid value for nwp_grp_rad%dt(1) (dt > 0 required).' )
    ENDIF

    ! Check mode of radiatively active gases
    istartitem = iUpatmoGasMode%startitem
    ienditem   = iUpatmoGasMode%enditem
    ! O3
    ! --
    lvalid = isInInterval(number=nwp_gas_o3%imode, opt_clbnd=istartitem, opt_cubnd=ienditem)
    IF (.NOT. lvalid)  CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o3%imode.' )
    ! O2
    ! --
    lvalid = isInInterval(number=nwp_gas_o2%imode, opt_clbnd=istartitem, opt_cubnd=ienditem)
    IF (.NOT. lvalid)  CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o2%imode.' )    
    ! O
    ! -
    lvalid = isInInterval(number=nwp_gas_o%imode, opt_clbnd=istartitem, opt_cubnd=ienditem)
    IF (.NOT. lvalid)  CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o%imode.' )
    ! CO2
    ! ---
    lvalid = isInInterval(number=nwp_gas_co2%imode, opt_clbnd=istartitem, opt_cubnd=ienditem)
    IF (.NOT. lvalid)  CALL finish( TRIM(routine), 'Invalid value for nwp_gas_co2%imode.' )
    ! NO
    ! --
    lvalid = isInInterval(number=nwp_gas_no%imode, opt_clbnd=istartitem, opt_cubnd=ienditem)
    IF (.NOT. lvalid)  CALL finish( TRIM(routine), 'Invalid value for nwp_gas_no%imode.' )

    ! Check volume mixing ratios (vmr >= 0 required, applies to 'const' mode only)
    ! O3
    ! --
    IF (nwp_gas_o3%imode == iUpatmoGasMode%const .AND. nwp_gas_o3%vmr < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o3%vmr (vmr >= 0 required).' )
    ENDIF
    ! O2
    ! --
    IF (nwp_gas_o2%imode == iUpatmoGasMode%const .AND. nwp_gas_o2%vmr < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o2%vmr (vmr >= 0 required).' )
    ENDIF
    ! O
    ! -
    IF (nwp_gas_o%imode == iUpatmoGasMode%const .AND. nwp_gas_o%vmr < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o%vmr (vmr >= 0 required).' )
    ENDIF
    ! CO2
    ! ---
    IF (nwp_gas_co2%imode == iUpatmoGasMode%const .AND. nwp_gas_co2%vmr < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_co2%vmr (vmr >= 0 required).' )
    ENDIF
    ! NO
    ! --
    IF (nwp_gas_no%imode == iUpatmoGasMode%const .AND. nwp_gas_no%vmr < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_no%vmr (vmr >= 0 required).' )
    ENDIF

    ! Check scaling factor for gas concentrations (fscale >= 0 required)
    ! O3
    ! --
    IF (nwp_gas_o3%fscale < 0._wp) THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o3%fscale (fscale >= 0 required).' )
    ENDIF
    ! O2
    ! --
    IF (nwp_gas_o2%fscale < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o2%fscale (fscale >= 0 required).' )
    ENDIF
    ! O
    ! -
    IF (nwp_gas_o%fscale < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_o%fscale (fscale >= 0 required).' )
    ENDIF
    ! CO2
    ! ---
    IF (nwp_gas_co2%fscale < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_co2%fscale (fscale >= 0 required).' )
    ENDIF
    ! NO
    ! --
    IF (nwp_gas_no%fscale < 0._wp)  THEN 
      CALL finish( TRIM(routine), 'Invalid value for nwp_gas_no%fscale (fscale >= 0 required).' )
    ENDIF    

    ! Check time interpolation update period for external data (dt > 0 required)
    ! gases
    ! -----
    IF (.NOT. (nwp_extdat_gases%dt > 0._wp)) THEN
      CALL finish( TRIM(routine), 'Invalid value for nwp_extdat_gases%dt (dt > 0 required).' )
    ENDIF
    ! chemheat
    ! --------
    IF (.NOT. (nwp_extdat_chemheat%dt > 0._wp)) THEN
      CALL finish( TRIM(routine), 'Invalid value for nwp_extdat_chemheat%dt (dt > 0 required).' )
    ENDIF

    ! Check "unofficial" switches
    lvalid = isInInterval(number=nwp_thermdyn_cpl, opt_clbnd=1, opt_cubnd=iThermdynCoupling%nitem)
    IF (.NOT. lvalid) CALL finish(TRIM(routine), 'Invalid nwp_thermdyn_cpl: ' &
      & //TRIM(int2string(nwp_thermdyn_cpl)))    

    !------------------------------------------------------------
    !              Fill the configuration state
    !------------------------------------------------------------

    DO jg = 0, max_dom

      !---------------
      ! Extrapolation
      !---------------

      upatmo_exp_config(jg)%expol_start_height   = expol_start_height  
      upatmo_exp_config(jg)%expol_blending_scale = expol_blending_scale
      upatmo_exp_config(jg)%expol_vn_decay_scale = expol_vn_decay_scale
      upatmo_exp_config(jg)%expol_temp_infty     = expol_temp_infty  
      upatmo_exp_config(jg)%lexpol_sanitycheck   = lexpol_sanitycheck
      ! Change status
      upatmo_exp_config(jg)%lset                 = .TRUE. 

      !---------------
      !    Physics
      !---------------

      upatmo_phy_config(jg)%orbit_type   = orbit_type
      upatmo_phy_config(jg)%solvar_type  = solvar_type
      upatmo_phy_config(jg)%solvar_data  = solvar_data
      upatmo_phy_config(jg)%solcyc_type  = solcyc_type
      upatmo_phy_config(jg)%cecc         = cecc
      upatmo_phy_config(jg)%cobld        = cobld
      upatmo_phy_config(jg)%clonp        = clonp
      upatmo_phy_config(jg)%lyr_perp     = lyr_perp
      upatmo_phy_config(jg)%yr_perp      = yr_perp
      upatmo_phy_config(jg)%lsanitycheck = lsanitycheck
      !
      ! ECHAM-specific:
      !
      startheightall = aes_start_height%all  ! }
      startheightimf = aes_start_height%imf  ! } Just for convenience
      startheightrad = aes_start_height%rad  ! }
      ! -----------------------
      ! Processes of group: IMF
      ! -----------------------
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%vdfmol )    =       & 
        & eval_start_height( start_height_process = aes_start_height%vdfmol,   &
        &                    start_height_group   = startheightimf,            &
        &                    start_height_all     = startheightall,            &
        &                    detect_entry_by      = iNonNegVal                 )
      !
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%fric )    =         &
        & eval_start_height(aes_start_height%fric, startheightimf, startheightall, iNonNegVal)
      !
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%iondrag ) =         &
        & eval_start_height(aes_start_height%iondrag, startheightimf, startheightall, iNonNegVal)
      ! Joule heating gets same start height as ion drag
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%joule )   =         &
        & upatmo_phy_config(jg)%aes_start_height(iUpatmoPrcId%iondrag)
      ! -----------------------
      ! Processes of group: RAD
      ! -----------------------
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%srbc )     =  & 
        & eval_start_height(aes_start_height%srbc, startheightrad, startheightall, iNonNegVal)
      !
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%euv )      =  & 
        & eval_start_height(aes_start_height%euv, startheightrad, startheightall, iNonNegVal)
      !
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%no )       =  & 
        & eval_start_height(aes_start_height%no, startheightrad, startheightall, iNonNegVal)
      !
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%chemheat ) =  & 
        & eval_start_height(aes_start_height%chemheat, startheightrad, startheightall, iNonNegVal)
      ! NLTE is not modifiable with regard to the start height
      upatmo_phy_config(jg)%aes_start_height( iUpatmoPrcId%nlte )     = -999._wp
      !
      ! NWP-specific:
      !
      IF (jg == 1) THEN  ! Primary domain

        ! ------------------
        ! Process group: IMF
        ! ------------------
        igrp = iUpatmoGrpId%imf
        upatmo_phy_config(jg)%nwp_grp(igrp)%imode        = nwp_grp_imf%imode(jg)
        upatmo_phy_config(jg)%nwp_grp(igrp)%dt           = nwp_grp_imf%dt(jg)
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_start      = nwp_grp_imf%t_start
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_end        = nwp_grp_imf%t_end
        upatmo_phy_config(jg)%nwp_grp(igrp)%start_height = nwp_grp_imf%start_height
        ! ------------------
        ! Process group: RAD
        ! ------------------
        igrp = iUpatmoGrpId%rad
        upatmo_phy_config(jg)%nwp_grp(igrp)%imode        = nwp_grp_rad%imode(jg)
        upatmo_phy_config(jg)%nwp_grp(igrp)%dt           = nwp_grp_rad%dt(jg)
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_start      = nwp_grp_rad%t_start
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_end        = nwp_grp_rad%t_end
        upatmo_phy_config(jg)%nwp_grp(igrp)%start_height = nwp_grp_rad%start_height

      ELSEIF (jg > 1) THEN  ! Nests

        ! ------------------
        ! Process group: IMF
        ! ------------------
        igrp = iUpatmoGrpId%imf
        IF (nwp_grp_imf%imode(jg) == iUpatmoPrcMode%unassigned) THEN
          ! No explicit namelist entry. Take value from previous domain
          upatmo_phy_config(jg)%nwp_grp(igrp)%imode = upatmo_phy_config(jg-1)%nwp_grp(igrp)%imode
        ELSE
          ! Explicit namelist entry
          upatmo_phy_config(jg)%nwp_grp(igrp)%imode = nwp_grp_imf%imode(jg)
        ENDIF
        IF (isInInterval(nwp_grp_imf%dt(jg), opt_olbnd=0._wp, opt_cubnd=nwp_grp_imf%dt(jg-1))) THEN
          ! Explicit (and valid) namelist entry
          upatmo_phy_config(jg)%nwp_grp(igrp)%dt = nwp_grp_imf%dt(jg)
        ELSE
          ! No explicit or valid namelist entry. Take value from prev. dom. 
          upatmo_phy_config(jg)%nwp_grp(igrp)%dt = upatmo_phy_config(jg-1)%nwp_grp(igrp)%dt
        ENDIF
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_start      = nwp_grp_imf%t_start
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_end        = nwp_grp_imf%t_end
        upatmo_phy_config(jg)%nwp_grp(igrp)%start_height = nwp_grp_imf%start_height
        ! ------------------
        ! Process group: RAD
        ! ------------------
        igrp = iUpatmoGrpId%rad
        IF (nwp_grp_rad%imode(jg) == iUpatmoPrcMode%unassigned) THEN
          upatmo_phy_config(jg)%nwp_grp(igrp)%imode = upatmo_phy_config(jg-1)%nwp_grp(igrp)%imode
        ELSE
          upatmo_phy_config(jg)%nwp_grp(igrp)%imode = nwp_grp_rad%imode(jg)
        ENDIF
        IF (isInInterval(nwp_grp_rad%dt(jg), opt_olbnd=0._wp, opt_cubnd=nwp_grp_rad%dt(jg-1))) THEN
          upatmo_phy_config(jg)%nwp_grp(igrp)%dt = nwp_grp_rad%dt(jg)
        ELSE
          upatmo_phy_config(jg)%nwp_grp(igrp)%dt = upatmo_phy_config(jg-1)%nwp_grp(igrp)%dt
        ENDIF
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_start      = nwp_grp_rad%t_start
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_end        = nwp_grp_rad%t_end
        upatmo_phy_config(jg)%nwp_grp(igrp)%start_height = nwp_grp_rad%start_height

      ELSEIF (jg == 0) THEN  ! Radiation domain

        ! Upper-atmosphere physics do not apply 
        ! ------------------
        ! Process group: IMF
        ! ------------------
        igrp = iUpatmoGrpId%imf
        upatmo_phy_config(jg)%nwp_grp(igrp)%imode        = iUpatmoPrcMode%off
        upatmo_phy_config(jg)%nwp_grp(igrp)%dt           = -999._wp
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_start      = ' '
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_end        = ' '
        upatmo_phy_config(jg)%nwp_grp(igrp)%start_height = -999._wp
        ! ------------------
        ! Process group: RAD
        ! ------------------
        igrp = iUpatmoGrpId%rad
        upatmo_phy_config(jg)%nwp_grp(igrp)%imode        = iUpatmoPrcMode%off
        upatmo_phy_config(jg)%nwp_grp(igrp)%dt           = -999._wp
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_start      = ' '
        upatmo_phy_config(jg)%nwp_grp(igrp)%t_end        = ' '
        upatmo_phy_config(jg)%nwp_grp(igrp)%start_height = -999._wp

      ELSE
        CALL finish(TRIM(routine), "Domain(s) jg < 0 not covered.")
      ENDIF

      sum_mmr = 0._wp
      ! --------------------------
      ! Radiatively active gas: O3
      ! --------------------------
      igas = iUpatmoGasId%o3
      upatmo_phy_config(jg)%nwp_gas(igas)%imode  = nwp_gas_o3%imode
      upatmo_phy_config(jg)%nwp_gas(igas)%vmr    = nwp_gas_o3%vmr * nwp_gas_o3%fscale ! Scaled volume mixing ratio [mol mol-1]
      upatmo_phy_config(jg)%nwp_gas(igas)%mmr    = ramd * amo3 * nwp_gas_o3%vmr       ! Mass mixing ratio [kg/kg]
      upatmo_phy_config(jg)%nwp_gas(igas)%fscale = nwp_gas_o3%fscale
      sum_mmr = sum_mmr + upatmo_phy_config(jg)%nwp_gas(igas)%mmr
      ! --------------------------
      ! Radiatively active gas: O2
      ! --------------------------
      igas = iUpatmoGasId%o2
      upatmo_phy_config(jg)%nwp_gas(igas)%imode  = nwp_gas_o2%imode
      upatmo_phy_config(jg)%nwp_gas(igas)%vmr    = nwp_gas_o2%vmr * nwp_gas_o2%fscale
      upatmo_phy_config(jg)%nwp_gas(igas)%mmr    = ramd * amo2 * nwp_gas_o2%vmr
      upatmo_phy_config(jg)%nwp_gas(igas)%fscale = nwp_gas_o2%fscale
      sum_mmr = sum_mmr + upatmo_phy_config(jg)%nwp_gas(igas)%mmr
      ! -------------------------
      ! Radiatively active gas: O
      ! -------------------------
      igas = iUpatmoGasId%o
      upatmo_phy_config(jg)%nwp_gas(igas)%imode  = nwp_gas_o%imode
      upatmo_phy_config(jg)%nwp_gas(igas)%vmr    = nwp_gas_o%vmr * nwp_gas_o%fscale
      upatmo_phy_config(jg)%nwp_gas(igas)%mmr    = ramd * amo * nwp_gas_o%vmr
      upatmo_phy_config(jg)%nwp_gas(igas)%fscale = nwp_gas_o%fscale
      sum_mmr = sum_mmr + upatmo_phy_config(jg)%nwp_gas(igas)%mmr
      ! ---------------------------
      ! Radiatively active gas: CO2
      ! ---------------------------
      igas = iUpatmoGasId%co2
      upatmo_phy_config(jg)%nwp_gas(igas)%imode  = nwp_gas_co2%imode
      upatmo_phy_config(jg)%nwp_gas(igas)%vmr    = nwp_gas_co2%vmr * nwp_gas_co2%fscale
      upatmo_phy_config(jg)%nwp_gas(igas)%mmr    = ramd * amco2 * nwp_gas_co2%vmr
      upatmo_phy_config(jg)%nwp_gas(igas)%fscale = nwp_gas_co2%fscale
      sum_mmr = sum_mmr + upatmo_phy_config(jg)%nwp_gas(igas)%mmr
      ! --------------------------
      ! Radiatively active gas: NO
      ! --------------------------
      igas = iUpatmoGasId%no
      upatmo_phy_config(jg)%nwp_gas(igas)%imode  = nwp_gas_no%imode
      upatmo_phy_config(jg)%nwp_gas(igas)%vmr    = nwp_gas_no%vmr * nwp_gas_no%fscale
      upatmo_phy_config(jg)%nwp_gas(igas)%mmr    = ramd * amno * nwp_gas_no%vmr
      upatmo_phy_config(jg)%nwp_gas(igas)%fscale = nwp_gas_no%fscale
      sum_mmr = sum_mmr + upatmo_phy_config(jg)%nwp_gas(igas)%mmr
      ! ---------------------------------------------
      ! Radiatively active gas: N2 => diagnostic mode
      ! ---------------------------------------------
      ! One of the radiatively active gases, currently N2, 
      ! is determined diagnostically:
      !
      ! mass(dry air) = sum[ mass(gas) ] = sum\[ mass(gas) ] + mass(N2) <=>
      !
      ! mass(N2) = mass(dry air) - sum\[ mass(gas) ] 
      !          = { 1 - sum\[ mmr(gas) ] } * mass(dry air), 
      !
      ! where sum\ denotes the sum over all gases except for the diagnostic gas N2, 
      ! and mmr is the mass mixing ratio of a gas. 
      ! From mass(N2), mass(dry air) >= 0, the following constraint follows:
      !
      ! sum\[ mmr(gas) ] <= 1.
      !
      ! If already those gases, for which the user made an entry for vmr in the namelist, 
      ! violate this constraint, we can stop here
      IF (sum_mmr > 1._wp) THEN
        CALL finish(TRIM(routine), 'Mass mixing ratios of gases sum up to >1, ' &
          & //'please reconsider your settings for vmr and fscale in upatmo_nml.')
      ENDIF
      igas = iUpatmoGasId%n2
      upatmo_phy_config(jg)%nwp_gas(igas)%imode  = iUpatmoGasMode%diag
      upatmo_phy_config(jg)%nwp_gas(igas)%vmr    = 0._wp
      upatmo_phy_config(jg)%nwp_gas(igas)%mmr    = 0._wp
      upatmo_phy_config(jg)%nwp_gas(igas)%fscale = 1._wp

      ! --------------------
      ! External data: gases
      ! --------------------
      iext = iUpatmoExtdatId%gases
      upatmo_phy_config(jg)%nwp_extdat(iext)%dt       = nwp_extdat_gases%dt
      upatmo_phy_config(jg)%nwp_extdat(iext)%filename = nwp_extdat_gases%filename
      ! -----------------------
      ! External data: chemheat
      ! -----------------------
      iext = iUpatmoExtdatId%chemheat
      upatmo_phy_config(jg)%nwp_extdat(iext)%dt       = nwp_extdat_chemheat%dt
      upatmo_phy_config(jg)%nwp_extdat(iext)%filename = nwp_extdat_chemheat%filename

      ! ----------------------
      ! Miscellaneous switches
      ! ----------------------
      upatmo_phy_config(jg)%nwp_thermdyn_cpl        = nwp_thermdyn_cpl
      upatmo_phy_config(jg)%nwp_ldiss_from_heatdiff = nwp_ldiss_from_heatdiff

      ! Change status
      upatmo_phy_config(jg)%lset = .TRUE.

    ENDDO  !jg

    !------------------------------------------------------------
    !              Store the namelist for restart
    !------------------------------------------------------------

    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=upatmo_nml)
      CALL store_and_close_namelist(funit, 'upatmo_nml') 
    ENDIF

    !------------------------------------------------------------
    !     Write the contents of the namelist to an ASCII file
    !------------------------------------------------------------

    IF(my_process_is_stdio()) WRITE(nnml_output,nml=upatmo_nml)        

  END SUBROUTINE read_upatmo_namelist

  !------------------------------------------------------------

  !>
  !! Function to evaluate the namelist input for the start heights 
  !! above which the upper-atmosphere physics processes compute tendencies 
  !! in case of ECHAM-forcing.
  !!
  FUNCTION eval_start_height( start_height_process, &
    &                         start_height_group,   &
    &                         start_height_all,     &
    &                         detect_entry_by       ) RESULT(start_height)
    
    ! In/out variables
    REAL(wp), INTENT(IN) :: start_height_process  ! Start height for single process
    REAL(wp), INTENT(IN) :: start_height_group    ! Start height for process group
    REAL(wp), INTENT(IN) :: start_height_all      ! Start height for all processes
    INTEGER,  INTENT(IN) :: detect_entry_by       ! How to detect a (valid) namelist entry
    ! 
    REAL(wp)             :: start_height
    
    ! Local variables
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
      &  routine = modname//':eval_start_height'

    !------------------------------------------------

    ! Which of the three start heights becomes the return value?
    !
    !                                   Entry for start_height_process?
    !                            Yes                   |
    !            --------------------------------------|
    !           |                                      | No
    !   start_height =                                 |
    !   start_height_process             Entry for start_height_group?
    !                            Yes                   |
    !            --------------------------------------|
    !           |                                      | No
    !   start_height =                                 |
    !   start_height_group                Entry for start_height_all? 
    !                            Yes                   |
    !            --------------------------------------|
    !           |                                      | No
    !   start_height =                                 |
    !   start_height_all                      start_height =   
    !                                         start_height_process
    !

    SELECT CASE(detect_entry_by)
    CASE(iNonNegVal)

      ! (Valid) namelist entries for the three start heights 
      ! are detected by non-negative values
      IF (.NOT. (start_height_process < 0._wp)) THEN
        start_height = start_height_process
      ELSEIF (.NOT. (start_height_group < 0._wp)) THEN
        start_height = start_height_group
      ELSEIF (.NOT. (start_height_all < 0._wp)) THEN
        start_height = start_height_all
      ELSE
        start_height = start_height_process
      ENDIF

    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid entry for detect_entry_by.')
    END SELECT

  END FUNCTION eval_start_height

END MODULE mo_upatmo_nml
