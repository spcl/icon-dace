!
! mo_art_io_constants
! This module provides constants for the I/O implementations
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

MODULE mo_art_io_constants
  USE mo_exception,        ONLY: finish

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: IART_FILENAMELEN

  PUBLIC :: IEMISS_ANTHRO, IEMISS_BIO, IEMISS_BBE
  PUBLIC :: IEXT_BIO, IEXT_SOIL, IEXT_POL, IEXT_PAL, IEXT_BCF
  PUBLIC :: IEXT_POCON, IEXT_POVAR, IEXT_POAMB
  PUBLIC :: IEXT_CHEM_SPEC
  PUBLIC :: IINIT_AERO, IINIT_CHEM
  PUBLIC :: IEXT_DMS
   
  PUBLIC :: IRES_LEN

  PUBLIC :: STORAGE_ATTR_SEP               !< separator for attributes in the XML file 
                                           !  stored as t_storage element

  PUBLIC :: art_get_abbr_string            !< returns the three character 
                                           !  abbreviation of the emission
  PUBLIC :: iemiss                         !< indices of prescribed emissions
  PUBLIC :: prescr_emiss_elements          !< possible element names for prescribed emissions
  PUBLIC :: art_get_all_emission_elements  !< function that writes prescr_emiss_elements into a 
                                           !  comma separated list with "or" at the end

  INTEGER, PARAMETER   ::  &
    & IART_FILENAMELEN = 240 !< Length of file names

  ! Integers to differentiate between i/o dataset types used to determine file name
  INTEGER, PARAMETER   ::  &
    & IEMISS_ANTHRO  = 10, & !< Anthropogenic emission dataset
    & IEMISS_BIO     = 11, & !< Biogenic emission dataset
    & IEMISS_BBE     = 12, & !< biomass burning (= fire) emission dataset
    & IEXT_BIO       = 20, & !< External parameters for biogenic emissions
    & IEXT_SOIL      = 21, & !< External parameters for soil (i.e. mineral dust emissions)
    & IEXT_POL       = 22, & !< External parameters for pollen emission
    & IEXT_PAL       = 23, & !< External parameters for altitude correction of pollen emission
    & IEXT_CHEM_SPEC = 24, & !< External data sets for chemical species
    & IEXT_BCF       = 25, & !< External data for flux of black carbon, GFAS-product
    & IINIT_AERO     = 30, & !< Initial data for aerosol
    & IINIT_CHEM     = 31, & !< Initial data for chemtracers
    & IEXT_DMS       = 32, & !< External dataset for onlie DMS calculation
    & IEXT_POCON     = 34, & !< External parameters for Pollen
    & IEXT_POVAR     = 35, & !< Variable parameters for Pollen emission
    & IEXT_POAMB     = 36    !< Daily sdes for ambrosia

  CHARACTER(LEN=1), PARAMETER ::   &
    &  STORAGE_ATTR_SEP = '@'

  INTEGER, PARAMETER   ::  &
    & IRES_LEN       = 6     !< length of resolution string

  INTEGER, PARAMETER  :: &
    &  iemiss(3) = (/ IEMISS_ANTHRO, IEMISS_BIO, IEMISS_BBE /)
  CHARACTER(LEN=9), PARAMETER   ::  &
    &  prescr_emiss_elements(3)  = (/ 'emiss_ANT',   &
                                      'emiss_BIO',   &
                                      'emiss_BBE'    /)

CONTAINS

CHARACTER(LEN=3) FUNCTION art_get_abbr_string(io_type)
!<
! FUNCTION art_get_io_abbr
! This function returns the three character abbreviated string of the integers
!
! Part of Module: mo_art_create_filenames
! Author: Michael Weimer, KIT
! Initial Release: 2016-10-27
! Modifications:
! 2016-12-21: Jonas Straub, KIT
! - included POL, PAL
!>
  IMPLICIT NONE
  INTEGER, INTENT(in)  :: io_type
  CHARACTER(LEN=3) :: typestring

  SELECT CASE(io_type)
    CASE(IEMISS_ANTHRO)
      typestring = 'ANT'
    CASE(IEMISS_BIO)
      typestring = 'BIO'
    CASE(IEMISS_BBE)
      typestring = 'BBE'
    CASE(IEXT_BIO)
      typestring = 'PFT'
    CASE(IEXT_SOIL)
      typestring = 'STY'
    CASE(IEXT_CHEM_SPEC)
      typestring = 'ECS'
    CASE(IEXT_POL)
      typestring = 'POL'
    CASE(IEXT_PAL)
      typestring = 'PAL'
    CASE(IEXT_BCF)
      typestring = 'BCF'  
    CASE(IINIT_AERO)
      typestring = 'IAE'
    CASE(IINIT_CHEM)
      typestring = 'ICE'
    CASE(IEXT_DMS)
      typestring = 'BIO'
    CASE(IEXT_POCON)
      typestring = 'POC'
    CASE(IEXT_POVAR)
      typestring = 'POV'
    CASE(IEXT_POAMB)
      typestring = 'AMB'
    CASE DEFAULT
      CALL finish('mo_io_constants:art_get_abbr_string', &
         &      'Unknown io_type')
  END SELECT

  art_get_abbr_string = typestring
  
END FUNCTION art_get_abbr_string
!
!--------------------------------------------------------------------------------------
!
CHARACTER(LEN = 100) FUNCTION art_get_all_emission_elements()
!<
! FUNCTION art_get_all_emission_elements
! This function writes the variable prescr_emiss_elements of mo_art_io_constants
! into a comma separated list with 'or' before the last array element
! Part of Module: mo_art_io_constants
! Author: Michael Weimer, KIT
! Initial Release: 2016-11-07
! Modifications:
! yyyy-mm-dd:
! - 
!>
  IMPLICIT NONE
  CHARACTER(LEN = 10) :: format_all_emissions
  INTEGER :: idx

  WRITE(format_all_emissions,'(A1,I5,A2)') '(',SIZE(prescr_emiss_elements),'A)'
  WRITE(art_get_all_emission_elements,format_all_emissions)  TRIM(prescr_emiss_elements(1)),  &
          & ((', '//TRIM(prescr_emiss_elements(idx))), idx=2,SIZE(prescr_emiss_elements)-1),  &
          & ' or '//TRIM(prescr_emiss_elements(SIZE(prescr_emiss_elements)))
END FUNCTION art_get_all_emission_elements
END MODULE mo_art_io_constants
