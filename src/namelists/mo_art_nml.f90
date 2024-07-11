! Namelist for ART-package
!
!
! Subroutine is called by read_atmo_namelists for setting up the ART-package
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

MODULE mo_art_nml
 
  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_run_config,          ONLY: lart
  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist,     &
    &                               open_and_restore_namelist, close_tmpfile
  USE mo_art_config,          ONLY: art_config, IART_PATH_LEN
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_art_namelist

  !----------------------------------!
  ! art_nml namelist variables       !
  !----------------------------------!

  ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_input_folder             !< Absolute Path to ART source code
  INTEGER :: iart_init_aero(1:max_dom)          !< Initialization of aerosol species
  INTEGER :: iart_init_gas(1:max_dom)           !< Initialization of gaseous species
  INTEGER :: iart_fplume             !< run FPlume model (Volcanic Plumes)
  INTEGER :: iart_volc_numb          !< number of volcanoes
  CHARACTER(LEN=IART_PATH_LEN)  :: cart_fplume_inp          
                                     !< path to FPlume input files (use without file extension)
  LOGICAL :: lart_diag_out           !< Enable output of diagnostic fields
  LOGICAL :: lart_pntSrc             !< Enables point sources
  LOGICAL :: lart_emiss_turbdiff     !< Switch if emissions should be included as surface flux condition
  CHARACTER(LEN=20) :: & 
   &  cart_io_suffix(1:max_dom)      !< user given suffix instead of automatically generated grid number 
                                     !  in ICON-ART input filename convention: 
                                     !  ART_iconR<n>B<kk>-grid-<yyyy-mm-dd-hh>_<grid_suffix>.nc

  ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
  LOGICAL :: lart_chem               !< Main switch to enable chemistry
  LOGICAL :: lart_chemtracer         !< switch for parametrised chemtracers
  LOGICAL :: lart_mecca              !< switch for MECCA chemistry
  LOGICAL :: lart_psc                !< switch for computation of PSCs 
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_vortex_init_date         !< Date of vortex initialization
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_cheminit_file(max_dom)   !< Path to chemical initialization file
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_cheminit_coord           !< Path to chemical initialization coordinate file
  CHARACTER(LEN=IART_PATH_LEN)  :: &
      &  cart_cheminit_type          !< type of chemical initialization coordinate file
  ! Paths and filenames of XML configuration
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_chemtracer_xml           !< Path to XML file for parametrised chemtracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_mecca_xml                !< Path to XML file for MECCA tracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_aerosol_xml              !< Path to XML file for aerosol tracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_modes_xml                !< Path to XML file for modes
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_pntSrc_xml               !< Path to XML file for point sources
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_coag_xml                 !< Path to XML file for coagulation(-matrix)
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_diagnostics_xml          !< Path to XML file for aerosol diagnostics (GRIB2 meta data)
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_emiss_xml_file           !< path and file name of the xml files for emission metadata
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_ext_data_xml             !< Path to XML file for metadata of datasets 
                                     !  that can prescribe tracers
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_aero_emiss_xml           !< Path to XML file for aerosol emission routines

  ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
  LOGICAL :: lart_aerosol            !< Main switch for the treatment of atmospheric aerosol
  INTEGER :: iart_seasalt            !< Treatment of sea salt aerosol
  INTEGER :: iart_dust               !< Treatment of mineral dust aerosol
  INTEGER :: iart_anthro             !< Treatment of anthropogenic aerosol
  INTEGER :: iart_fire               !< Treatment of wildfire aerosol
  INTEGER :: iart_volcano            !< Treatment of volcanic ash aerosol
  INTEGER :: iart_nonsph             !< Treatment of nonspherical particles
  INTEGER :: iart_isorropia          !< Treatment of aerosol gas partioning
  INTEGER :: iart_seas_water         !< Calculation of seasalt water content
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_volcano_file             !< Absolute path + filename of input file for volcanoes
  INTEGER :: iart_radioact           !< Treatment of radioactive particles
  CHARACTER(LEN=IART_PATH_LEN)  :: &
    &  cart_radioact_file            !< Absolute path + filename of input file for radioactive emissions
  INTEGER :: iart_pollen             !< Treatment of pollen
  INTEGER :: iart_modeshift          !< Doing mode shift (only temporary switch for debug)
    
  ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
  INTEGER :: iart_aci_warm           !< Nucleation of aerosol to cloud droplets
  INTEGER :: iart_aci_cold           !< Nucleation of aerosol to cloud ice
  INTEGER :: iart_ari                !< Direct interaction of aerosol with radiation

  LOGICAL :: lart_dusty_cirrus       !< Dusty cirrus parameterization in cloud cover scheme
  REAL(wp):: rart_dustyci_crit            !< Dust threshold for dusty cirrus  [mug/kg]
  REAL(wp):: rart_dustyci_rhi             !< RHi  threshold for dusty cirrus  [-]

  ! Treatment of grid scale and convective precipitation in dust washout
  INTEGER :: iart_aero_washout       !< 0:gscp+con; 1:gscp,con; 2:gscp,rcucov*con

  ! Number of substeps for sedimentation
  INTEGER :: nart_substeps_sedi(1:max_dom)
  CHARACTER(LEN=4) :: cart_type_sedim  !< type of sedimentation scheme: "expl": explicit, "impl": implicit

  ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
  LOGICAL :: lart_conv               !< Convection of aerosol (TRUE/FALSE)
  LOGICAL :: lart_turb               !< Turbulent diffusion of aerosol (TRUE/FALSE)

  ! Restart-DEBUG: Write DEBUG-Restartfile
  LOGICAL :: lart_debugRestart

  ! Time interval over which maximum of air concentration of radionuclides is taken
  REAL(wp):: radioact_maxtint(1:max_dom)


  NAMELIST/art_nml/ cart_input_folder, lart_chem, lart_chemtracer, lart_mecca,          &
   &                cart_io_suffix, lart_pntSrc, lart_aerosol, iart_seasalt, iart_dust, &
   &                iart_anthro, iart_fire, iart_volcano, cart_volcano_file,            &
   &                iart_fplume, iart_volc_numb, cart_fplume_inp, iart_radioact,        &
   &                cart_radioact_file, iart_pollen, iart_nonsph, iart_isorropia,       &
   &                iart_seas_water, lart_dusty_cirrus, rart_dustyci_crit,              &
   &                rart_dustyci_rhi, &
   &                iart_modeshift, iart_aci_warm, iart_aci_cold, iart_ari,             &
   &                iart_aero_washout, lart_conv, lart_turb, iart_init_aero,            &
   &                iart_init_gas, lart_diag_out, cart_emiss_xml_file,                  &
   &                cart_ext_data_xml, cart_vortex_init_date , cart_cheminit_file,      &
   &                cart_cheminit_coord, cart_cheminit_type,                            &
   &                lart_emiss_turbdiff, nart_substeps_sedi,                            &
   &                cart_chemtracer_xml, cart_mecca_xml, cart_aerosol_xml,              &
   &                cart_modes_xml, cart_pntSrc_xml, cart_diagnostics_xml,              &
   &                lart_psc, cart_coag_xml, cart_aero_emiss_xml, cart_type_sedim,      &
   &                lart_debugRestart, radioact_maxtint

CONTAINS
  !-------------------------------------------------------------------------
  !>
  !! Read Namelist for ART-package.
  !!
  !! This subroutine
  !! - reads the Namelist for the ART-package
  !! - sets default values
  !! - potentially overwrites the defaults by values used in a
  !!   previous integration (if this is a resumed run)
  !! - reads the user's (new) specifications
  !! - stores the Namelist for restart
  !! - fills the configuration state (partly)
  !!
  SUBROUTINE read_art_namelist( filename )

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat, funit
    INTEGER :: jg            !< patch loop index
    LOGICAL :: l_exist       !< variable for inquiring if the xml file 
                                        !   and emission base path exist.
    CHARACTER(len=*), PARAMETER ::  &
      &  routine = 'mo_art_nml: read_art_nml'
    INTEGER :: iunit

    !-----------------------
    ! 1. default settings
    !-----------------------

    ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
    cart_input_folder          = ''
    iart_init_aero(:)          = 0
    iart_init_gas(:)           = 0
    lart_diag_out              = .FALSE.
    lart_pntSrc                = .FALSE.
    lart_emiss_turbdiff        = .FALSE.
    cart_io_suffix(1:max_dom)  = 'grid-number'
    iart_fplume                = 0
    iart_volc_numb             = 0
    cart_fplume_inp            = ''
 
    ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
    lart_chem             = .FALSE.
    lart_chemtracer       = .FALSE.
    lart_mecca            = .FALSE.
    lart_psc              = .FALSE.
    cart_vortex_init_date = ''
    cart_cheminit_file(:) = ''
    cart_cheminit_coord   = ''
    cart_cheminit_type    = ''

    ! Paths and filenames of XML configuration
    cart_chemtracer_xml   = ''
    cart_mecca_xml        = ''
    cart_aerosol_xml      = ''
    cart_modes_xml        = ''
    cart_pntSrc_xml       = ''
    cart_coag_xml         = ''
    cart_diagnostics_xml  = ''
    cart_emiss_xml_file   = ''
    cart_ext_data_xml     = ''
    cart_aero_emiss_xml   = ''

    ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
    lart_aerosol        = .FALSE.
    iart_seasalt        = 0
    iart_dust           = 0
    iart_anthro         = 0
    iart_fire           = 0
    iart_volcano        = 0
    cart_volcano_file   = ''
    iart_radioact       = 0
    cart_radioact_file  = ''
    iart_pollen         = 0
    iart_modeshift      = 0
    iart_nonsph         = 0
    iart_isorropia      = 0
    iart_seas_water     = 0

    ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
    iart_aci_warm       = 0
    iart_aci_cold       = 0
    iart_ari            = 0
    ! Dusty cirrus
    lart_dusty_cirrus   = .FALSE.
    rart_dustyci_crit   = 70.0_wp
    rart_dustyci_rhi    = 0.90_wp

    ! Treatment of grid scale and convective precipitation in dust washout
    iart_aero_washout   = 0

    ! Number of substeps for sedimentation
    nart_substeps_sedi(:) = 2

    cart_type_sedim       = "expl"

    ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
    lart_conv           = .TRUE.
    lart_turb           = .TRUE.

    ! Write DEBUG-Restartfile
    lart_debugRestart   = .FALSE.

    ! Time interval over which maximum of air concentration of radionuclides is taken
    radioact_maxtint(:) = 3600._wp

    !------------------------------------------------------------------
    ! 2. If this is a resumed integration, overwrite the defaults above
    !    by values used in the previous integration.
    !------------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('art_nml')
      READ(funit,NML=art_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! 3. Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('art_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, art_nml)    ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)

      ! Set array parameters to dummy values to determine which ones are actively set in the namelist
      iart_init_aero(:)     = -1
      iart_init_gas(:)      = -1
      nart_substeps_sedi(:) = -1
      radioact_maxtint(:)   = -1._wp

      READ (nnml, art_nml)                                        ! overwrite default settings

      ! Restore default values for global domain where nothing at all has been specified
      IF (iart_init_aero(1) < 0)     iart_init_aero(1)     = 0
      IF (iart_init_gas(1) < 0)      iart_init_gas(1)      = 0
      IF (nart_substeps_sedi(1) < 0) nart_substeps_sedi(1) = 2    ! number of substeps for sedimentation
      IF (radioact_maxtint(1) < 0._wp) radioact_maxtint(1) = 3600._wp

      ! Copy values of parent domain (in case of linear nesting) to nested domains where nothing has been specified
      DO jg = 2, max_dom
        IF (iart_init_aero(jg) < 0)     iart_init_aero(jg)     = iart_init_aero(jg-1)
        IF (iart_init_gas(jg) < 0)      iart_init_gas(jg)      = iart_init_gas(jg-1)
        IF (nart_substeps_sedi(jg) < 0) nart_substeps_sedi(jg) = nart_substeps_sedi(jg-1)
        IF (radioact_maxtint(jg) < 0._wp) radioact_maxtint(jg) = radioact_maxtint(jg-1)
      ENDDO

      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, art_nml)    ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml



    !----------------------------------------------------
    ! 4. Sanity check (only if lart is true)
    !----------------------------------------------------

    IF (lart) THEN
    
      IF (iart_aci_cold == 6 .AND. iart_dust == 0) THEN
        CALL finish('mo_art_nml:read_art_namelist',  &
          &         'Invalid combination: iart_aci_cold = 6 and iart_dust = 0')
      ENDIF
      IF (iart_aci_cold == 7 .AND. iart_dust == 0) THEN
        CALL finish('mo_art_nml:read_art_namelist',  &
          &         'Invalid combination: iart_aci_cold = 7 and iart_dust = 0')
      ENDIF
      IF (lart_dusty_cirrus .AND. iart_dust == 0) THEN
        CALL finish('mo_art_nml:read_art_namelist',  &
          &         'Invalid combination: lart_dusty_cirrus = .TRUE. and iart_dust = 0')
      ENDIF
  
      ! Emission paths and file
      IF (TRIM(cart_emiss_xml_file) /= '') THEN
        INQUIRE(file = TRIM(cart_emiss_xml_file), EXIST = l_exist)
      
        IF (.NOT. l_exist) THEN
          CALL finish('mo_art_nml:read_art_namelist',  &
                      TRIM(cart_emiss_xml_file)//  &
                      & ' could not be found. Check cart_emiss_xml_file.')
        END IF
      END IF
  

      ! Diagnostics paths and file
      IF (TRIM(cart_diagnostics_xml) /= '') THEN
        INQUIRE(file = TRIM(cart_diagnostics_xml), EXIST = l_exist)

        IF (.NOT. l_exist) THEN
          CALL finish('mo_art_nml:read_art_namelist',  &
                      TRIM(cart_diagnostics_xml)//  &
                      & ' could not be found. Check cart_diagnostics_xml.')
        END IF
      END IF

      ! FPLUME input path
      IF (iart_fplume>=1) THEN
        IF(TRIM(cart_fplume_inp) == '') THEN
          CALL finish('mo_art_nml:read_art_namelist','namelist parameter cart_fplume_inp' &       
                    //' has to be given for iart_fplume>=1')
        END IF
        IF (iart_volc_numb==0) iart_volc_numb = 1
      END IF

    END IF  ! lart


    !----------------------------------------------------
    ! 5. Fill the configuration state
    !----------------------------------------------------

    DO jg= 1,max_dom !< Do not take into account reduced radiation grid
      ! General variables (Details: cf. Tab. 2.1 ICON-ART User Guide)
      art_config(jg)%cart_input_folder   = TRIM(cart_input_folder)
      art_config(jg)%iart_init_aero      = iart_init_aero(jg)
      art_config(jg)%iart_init_gas       = iart_init_gas(jg)
      art_config(jg)%lart_diag_out       = lart_diag_out
      art_config(jg)%lart_pntSrc         = lart_pntSrc
      art_config(jg)%lart_emiss_turbdiff = lart_emiss_turbdiff
      art_config(jg)%cart_io_suffix      = TRIM(cart_io_suffix(jg))
      art_config(jg)%iart_fplume         = iart_fplume
      art_config(jg)%iart_volc_numb      = iart_volc_numb
      art_config(jg)%cart_fplume_inp     = TRIM(cart_fplume_inp)
      
      ! Atmospheric Chemistry (Details: cf. Tab. 2.2 ICON-ART User Guide)
      art_config(jg)%lart_chem             = lart_chem
      art_config(jg)%lart_chemtracer       = lart_chemtracer
      art_config(jg)%lart_mecca            = lart_mecca
      art_config(jg)%lart_psc              = lart_psc
      art_config(jg)%cart_vortex_init_date = TRIM(cart_vortex_init_date)
      art_config(jg)%cart_cheminit_file    = TRIM(cart_cheminit_file(jg))
      art_config(jg)%cart_cheminit_coord   = TRIM(cart_cheminit_coord)
      art_config(jg)%cart_cheminit_type    = TRIM(cart_cheminit_type)

      ! Paths and filenames of XML configuration
      art_config(jg)%cart_chemtracer_xml   = TRIM(cart_chemtracer_xml)
      art_config(jg)%cart_mecca_xml        = TRIM(cart_mecca_xml)
      art_config(jg)%cart_aerosol_xml      = TRIM(cart_aerosol_xml)
      art_config(jg)%cart_modes_xml        = TRIM(cart_modes_xml)
      art_config(jg)%cart_pntSrc_xml       = TRIM(cart_pntSrc_xml)
      art_config(jg)%cart_coag_xml         = TRIM(cart_coag_xml)
      art_config(jg)%cart_diagnostics_xml  = TRIM(cart_diagnostics_xml)
      art_config(jg)%cart_emiss_xml_file   = TRIM(cart_emiss_xml_file)
      art_config(jg)%cart_ext_data_xml     = TRIM(cart_ext_data_xml)
      art_config(jg)%cart_aero_emiss_xml   = TRIM(cart_aero_emiss_xml)

      ! Atmospheric Aerosol (Details: cf. Tab. 2.3 ICON-ART User Guide)
      art_config(jg)%lart_aerosol        = lart_aerosol
      art_config(jg)%iart_seasalt        = iart_seasalt
      art_config(jg)%iart_dust           = iart_dust
      art_config(jg)%iart_anthro         = iart_anthro
      art_config(jg)%iart_fire           = iart_fire
      art_config(jg)%iart_volcano        = iart_volcano
      art_config(jg)%iart_nonsph         = iart_nonsph
      art_config(jg)%iart_isorropia      = iart_isorropia
      art_config(jg)%iart_seas_water     = iart_seas_water
      art_config(jg)%cart_volcano_file   = TRIM(cart_volcano_file)
      art_config(jg)%iart_radioact       = iart_radioact
      art_config(jg)%cart_radioact_file  = TRIM(cart_radioact_file)
      art_config(jg)%iart_pollen         = iart_pollen
      art_config(jg)%iart_modeshift      = iart_modeshift
      
      ! Feedback processes (Details: cf. Tab. 2.4 ICON-ART User Guide)
      art_config(jg)%iart_aci_warm       = iart_aci_warm
      art_config(jg)%iart_aci_cold       = iart_aci_cold
      art_config(jg)%iart_ari            = iart_ari
      ! Dusty cirrus
      art_config(jg)%lart_dusty_cirrus   = lart_dusty_cirrus
      art_config(jg)%rart_dustyci_crit   = rart_dustyci_crit
      art_config(jg)%rart_dustyci_rhi    = rart_dustyci_rhi

      ! Treatment of grid scale and convective precipitation in dust washout
      art_config(jg)%iart_aero_washout   = iart_aero_washout

      ! Number of substeps for sedimentation
      art_config(jg)%nart_substeps_sedi  = nart_substeps_sedi(jg)

      art_config(jg)%cart_type_sedim     = cart_type_sedim

      ! Fast Physics Processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
      art_config(jg)%lart_conv           = lart_conv
      art_config(jg)%lart_turb           = lart_turb

      ! Write DEBUG-Restartfile
      art_config(jg)%lart_debugRestart   = lart_debugRestart

      ! Time interval over which maximum of air concentration of radionuclides is taken
      art_config(jg)%radioact_maxtint    = radioact_maxtint(jg)
    ENDDO !jg

    !-----------------------------------------------------
    ! 6. Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=art_nml)
      CALL store_and_close_namelist(funit, 'art_nml')
    ENDIF

    ! 7. write the contents of the namelist to an ASCII file
    !
    IF(my_process_is_stdio()) WRITE(nnml_output,nml=art_nml)


  END SUBROUTINE read_art_namelist

END MODULE mo_art_nml
