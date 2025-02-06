!
! mo_art_prescribed_types
! This module provides data structures for general prescribed data
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

MODULE mo_art_prescribed_types
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_art_config,                    ONLY: IART_PATH_LEN
  USE mtime,                            ONLY: datetime
! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_chem_init_types,           ONLY: t_chem_init_state

  IMPLICIT NONE

  PUBLIC :: t_art_emiss_prescribed
  PUBLIC :: t_art_chem_prescribed
  PUBLIC :: t_art_prescribed
  PUBLIC :: t_art_prescr_list
  PUBLIC :: t_art_prescr_list_element
  PUBLIC :: t_art_prescr_var_dep

  PRIVATE

  TYPE dt_container
    TYPE(datetime), POINTER :: &
      &  ptr => NULL()   !< datetimes before, at and after simulation time
  END TYPE dt_container
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_prescr_var_dep  ! variable dependent values (tracer index etc.)
    INTEGER ::        &
      &  tracer_idx,  &  !< index of the tracer
      &  bot_idx,     &  !< bottom and upper indices where prescribing has to be performed
      &  upp_idx
    REAL(wp) ::       &
      &  mol_weight      !< molar weight in kg / mol
  END TYPE t_art_prescr_var_dep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_prescribed
    TYPE(dt_container) ::  &
      & datetime(3)       !< datetimes before, at and after simulation time
    LOGICAL  ::  &
      &  type_is_init = .FALSE.,  &  !< flag if this is initialised
      &  first_date_is_leap_year, &  !< flag if first date is a leap year
      &  last_date_is_leap_year      !< flag if last  date is a leap year

    TYPE(datetime), POINTER  :: &
      &  first_date => NULL(),  &   !< start date of the dataset
      &  last_date  => NULL()       !< end date of the dataset

    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  path                 !< full path to the emission dataset
    INTEGER :: &
      &  num_vars,          & !< number of variables in the netCDF file
      &  iType_data           !< integer of emission type (see io/mo_art_io_constants)
    CHARACTER(:), ALLOCATABLE :: &
      &  id                   !< id in the emission xml file
    CHARACTER(LEN=IART_VARNAMELEN), ALLOCATABLE  ::  &
      &  vname(:)             !< names of the variables
  END TYPE t_art_prescribed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE, EXTENDS(t_art_prescribed) :: t_art_chem_prescribed
    REAL(wp), ALLOCATABLE :: &
      &  vinterp_3d(:,:,:,:,:)     !< prescribed data (:,:,:,num_vars,1): emission before sim. time
                                   !                  (:,:,:,num_vars,2): emission at     sim. time
                                   !                  (:,:,:,num_vars,3): emission after  sim. time
    CHARACTER(LEN=IART_VARNAMELEN) :: &
      &  unit_vertical        !< unit of the vertical coordinate in the netCDF file 
                              !  (currently one of 'mbar','hPa','Pa')
    TYPE(t_chem_init_state) :: &
      &  prescribed_state     !< containing all information about the vertical interpolation
                              !  (especially weights wrt geopotential)
    TYPE(t_art_prescr_var_dep), ALLOCATABLE  ::  &
      &  var_dep(:)           !< variable dependent values (tracer index etc.)
    CHARACTER(LEN=IART_PATH_LEN) :: &
      &  vert_coord_filename  !< (only) name of the file containing information 
                              !  about the vertical hybrid coordinate
  END TYPE t_art_chem_prescribed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE, EXTENDS(t_art_prescribed) :: t_art_emiss_prescribed
    REAL(wp), ALLOCATABLE :: &
      &  emiss_2d_read(:,:)   !< just temporal storage of read emissions

    REAL(wp), ALLOCATABLE :: &
      &  num_dims(:),        &!< this is only temporarily. If it is 3 it should be included 
                              !  to the list of t_art_chem_prescribed
      &  emiss_2d(:,:,:)      !< emission in kg m-2 s-1 (:,:,1): emission before simulation time
                              !                         (:,:,2): emission at     simulation time
                              !                         (:,:,3): emission after  simulation time
    INTEGER  :: &
      &  num_emiss_lev        !< number of lowest model levels into which emission shall be included
    REAL(wp) :: &
      &  scaling_factor       !< scaling factor with which the emission is multipled (usually 1)
  END TYPE t_art_emiss_prescribed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! linked list for prescribed data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_prescr_list_element
    TYPE(t_art_chem_prescribed) :: &
      &  prescr                        !< part of the linked list element that contains 
                                       !  the meta information about the data set
    TYPE(t_art_prescr_list_element), POINTER :: &
      &  next_prescr_element => NULL() !< next element
  END TYPE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_prescr_list
    TYPE(t_art_prescr_list_element), POINTER :: &
      &  first_prescr_element => NULL() !< first element
    INTEGER :: &
      &  num_elements = 0               !< number of elements in the list
  END TYPE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE mo_art_prescribed_types
