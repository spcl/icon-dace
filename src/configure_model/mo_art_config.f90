! @brief configuration setup for ART-package
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

MODULE mo_art_config

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: max_dom
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: nart_tendphy
  PUBLIC :: npreslay
  PUBLIC :: t_art_config
  PUBLIC :: art_config
  PUBLIC :: configure_art
  PUBLIC :: IART_PATH_LEN

  PUBLIC :: ctracer_art                                                  !generic tracer list, contains name of tracer

  !!--------------------------------------------------------------------------
  !! Basic configuration setup for ICON-ART
  !!--------------------------------------------------------------------------
  INTEGER, PARAMETER  :: npreslay      = 7   !Number of pressure layers for diagnostic output of maximum concentration
  INTEGER, PARAMETER  :: IART_PATH_LEN = 200 !Maximum length of file- and pathnames

  INTEGER             :: nart_tendphy  = 0  !Maximum number of tracers that are effected by deep convective transport 
  
  
  TYPE t_art_config ! Namelist variables for ART
    
    ! Namelist variables
    
    ! General control variables (Details: cf. Tab. 2.2 ICON-ART User Guide)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_input_folder               !< Absolute Path to ART input files 
    INTEGER :: iart_ntracer              !< number of transported ART tracers
    INTEGER :: iart_init_aero            !< Initialization of aerosol species
    INTEGER :: iart_init_gas             !< Initialization of gaseous species
    INTEGER :: iart_fplume               !< run FPlume model (Volcanic plumes)
    INTEGER :: iart_volc_numb            !< number of volcanoes in FPlume
    CHARACTER(LEN=IART_PATH_LEN)  ::   &
      &  cart_fplume_inp                 !< Path to FPlume input files (insert without extension)
    LOGICAL :: lart_diag_out             !< Enable output of diagnostic fields
    LOGICAL :: lart_pntSrc               !< Enables point sources
    LOGICAL :: lart_emiss_turbdiff       !< Switch if emission should be included as surface flux condition
    CHARACTER(LEN=20) :: cart_io_suffix  !< user given suffix instead of automatically generated grid number 
                                         !  in ICON-ART input filename convention: 
                                         !  ART_iconR<n>B<kk>-grid-<yyyy-mm-dd-hh>_<grid_suffix>.nc
    
    ! Atmospheric Chemistry (Details: cf. Tab. 2.3 ICON-ART User Guide)
    LOGICAL :: lart_chem               !< Main switch to enable chemistry
    LOGICAL :: lart_chemtracer         !< switch for parametrised tracers
    LOGICAL :: lart_mecca              !< switch for MECCA chemistry
    LOGICAL :: lart_psc                !< switch for computation of PSCs
    INTEGER :: O3_feedback             !< O3 radiation feedback from ART? (0 = false, 1 = true)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_vortex_init_date        !< Date of vortex initialization
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_cheminit_file           !< Path to chemical initialization file
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_cheminit_coord          !< Path to chemical initialization coordinate file
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_cheminit_type           !< Type of  chemical initialization coordinate file
    ! Paths and filenames of XML configuration
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_chemtracer_xml           !< Path to XML file for parametrised chemtracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_mecca_xml               !< Path to XML file for MECCA chemistry tracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_aerosol_xml             !< Path to XML file for aerosol tracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_modes_xml               !< Path to XML file for modes
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_pntSrc_xml              !< Path to XML file for point sources
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_coag_xml                !< Path to XML file for coagulation(-matrix)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_diagnostics_xml         !< Path to XML file for aerosol diagnostics (GRIB2 meta data)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_emiss_xml_file          !< Path and file name of the xml file containing meta information of the emissions.
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_ext_data_xml            !< Path to XML file with metadata of
                                      !  datasets which can prescribe tracers
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_aero_emiss_xml          !< Path to XML file for aerosol emission routines

    ! Atmospheric Aerosol (Details: cf. Tab. 2.4 ICON-ART User Guide)
    LOGICAL :: lart_aerosol            !< Main switch for the treatment of atmospheric aerosol
    INTEGER :: iart_seasalt            !< Treatment of sea salt aerosol
    INTEGER :: iart_dust               !< Treatment of mineral dust aerosol
    INTEGER :: iart_anthro             !< Treatment of anthropogenic aerosol
    INTEGER :: iart_fire               !< Treatment of wildfire aerosol
    INTEGER :: iart_volcano            !< Treatment of volcanic ash aerosol
    INTEGER :: iart_nonsph             !< Treatment of nonspherical particles
    INTEGER :: iart_isorropia          !< Treatment of gas aerosol partioning
    INTEGER :: iart_seas_water         !< Calculation of seasalt water content
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_volcano_file             !< Absolute path + filename of input file for volcanoes
    INTEGER :: iart_radioact           !< Treatment of radioactive particles
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  cart_radioact_file            !< Absolute path + filename of input file for radioactive emissions
    INTEGER :: iart_pollen             !< Treatment of pollen
    INTEGER :: iart_modeshift          !< Doing mode shift (only temporary switch for debug)
    
    ! Feedback processes (Details: cf. Tab. 2.5 ICON-ART User Guide)
    INTEGER :: iart_aci_warm           !< Nucleation of aerosol to cloud droplets
    INTEGER :: iart_aci_cold           !< Nucleation of aerosol to cloud ice
    INTEGER :: iart_ari                !< Direct interaction of aerosol with radiation

    LOGICAL :: lart_dusty_cirrus       !< Dusty cirrus parameterization in cloud cover scheme
    REAL(wp):: rart_dustyci_crit       !< Dust threshold for dusty cirrus  [mug/kg]
    REAL(wp):: rart_dustyci_rhi        !< RHi  threshold for dusty cirrus  [-]

    ! Treatment of grid scale and convective precipitation in dust washout
    INTEGER :: iart_aero_washout       !< 0:gscp+con; 1:gscp,con; 2:gscp,rcucov*con

    ! Number of substeps for sedimentation
    INTEGER :: nart_substeps_sedi

    CHARACTER(LEN=4) :: cart_type_sedim  !< type of sedimentation scheme: "expl": explicit, "impl": implicit

    ! Fast Physics Processes (Details: cf. Tab. 2.6 ICON-ART User Guide)
    LOGICAL :: lart_conv               !< Convection of aerosol (TRUE/FALSE)
    INTEGER :: nconv_tracer            !< number of tracers in convection 
    LOGICAL :: lart_turb               !< Turbulent diffusion of aerosol (TRUE/FALSE)
    INTEGER :: nturb_tracer            !< number of tracers in turbulence

    ! Write DEBUG-Restartfile
    LOGICAL :: lart_debugRestart
    
    ! Time interval over which maximum of air concentration of radionuclides is taken
    REAL(wp):: radioact_maxtint
  END TYPE t_art_config

  TYPE(t_art_config), TARGET :: art_config(0:max_dom)


    !------------- generic tracer list for NWP and ECHAM physics
  CHARACTER(len=MAX_CHAR_LENGTH), ALLOCATABLE :: ctracer_art(:) !< tracer acronyms

CONTAINS

  !>
  !! setup components of ICON-ART depending on this namelist
  !!
  !! Setup of additional ICON-ART control variables depending on the 
  !! art-NAMELIST and potentially other namelists. This routine is 
  !! called, after all namelists have been read and a synoptic consistency 
  !! check has been done.
  !!
  SUBROUTINE configure_art(jg)
    INTEGER, INTENT(IN) :: jg          !< patch 

    art_config(jg)%nconv_tracer=0
    art_config(jg)%nturb_tracer=0
 
  END SUBROUTINE configure_art


END MODULE mo_art_config
