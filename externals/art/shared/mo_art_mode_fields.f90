!
! mo_art_mode_fields
! This module provides the base type t_mode_fields, its parameters and initialization structures
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

MODULE mo_art_mode_fields
! ICON
  USE mo_exception,                     ONLY: message, finish, message_text
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
! ART aerosol dynamics processes routines

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_mode_fields'


!-------------------------------------------------------------------------------
!-------       Hierarchy tree for t_mode_fields polymorphic objects      -------
!-------------------------------------------------------------------------------
!                                  t_mode_fields
!                t_fields_2mom <|                  |> t_fields_1mom
!                                   t_fields_pollen | t_fields_volc | t_fields_radio
!-------------------------------------------------------------------------------
  
  TYPE :: t_mode_fields
    CHARACTER(LEN=IART_VARNAMELEN) :: name    !< Mode name
    LOGICAL                        :: linit   !< initial values have been set
    !
    CONTAINS
      PROCEDURE :: create   => create_t_mode_fields
      PROCEDURE :: set_meta => set_meta_t_mode_fields
      PROCEDURE :: print    => print_meta_t_mode_fields
  END TYPE t_mode_fields
  !
  ! -------------------------------------------------------------------------
  !
  TYPE, extends (t_mode_fields) :: t_art_map2tracer
    INTEGER               :: &
      &  ntr     = 0,        & !< Number of tracers for which emissions are provided
      &  nmodes  = 0,        & !< Number of modes required by this routine (i.e. routine outputs)
      &  nsub    = 0           !< Number of substances in XML-File (most of the time equal to ntr)
    CHARACTER(LEN=IART_VARNAMELEN), ALLOCATABLE :: &
      &  substance(:)          !< substances to be emitted
    CHARACTER(LEN=IART_VARNAMELEN), ALLOCATABLE :: &
      &  modenames(:)          !< Names of modes which are emitted into(nmodes)
    INTEGER, ALLOCATABLE  :: &
      &  itr0(:),            & !< Index of number tracers in container (nmodes)
      &  itr3(:,:),          & !< Index of mass tracers in container (nmodes,ntr)
      &  ntrpermode(:)         !< support-structure for determining ntr
  END TYPE t_art_map2tracer


  PUBLIC :: t_mode_fields, t_art_map2tracer
  
  CONTAINS

  SUBROUTINE create_t_mode_fields(this_fields, modename, idims)
  !<
  ! SUBROUTINE create_t_mode_fields
  ! This subroutine performs a setup for the fields of type t_mode_fields
  ! Part of Module: mo_art_mode_fields (moved from mo_art_modes)
  ! Author: Daniel Rieger, KIT
  ! Initial Release: 2016-08-22
  ! Modifications:
  ! 2020-02-14: Sven Werchner, KIT
  ! - moved to mo_art_mode_fields.f90
  !>
    CLASS(t_mode_fields),INTENT(inout) :: &
      &  this_fields                        !< Container with fields
    CHARACTER(LEN=*),INTENT(in)        :: &
      &  modename                           !< Name of current mode
    INTEGER,INTENT(in)                 :: &
      &  idims(3)                           !< Dimensions to allocate fields
  
    this_fields%linit = .FALSE.
    IF(LEN_TRIM(modename) > IART_VARNAMELEN) THEN
      CALL finish(TRIM(routine)//':create_t_mode_fields',          &
        &         'Length of name for '//TRIM(modename)//' to long.')
    ENDIF
    this_fields%name  = TRIM(modename)
  
  END SUBROUTINE create_t_mode_fields
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE set_meta_t_mode_fields(this_fields, meta_key_value_store)
  !<
  ! SUBROUTINE set_meta_t_mode_fields
  ! This subroutine sets the metadata associated to type t_mode_fields
  ! Part of Module: mo_art_mode_fields (moved from mo_art_modes)
  ! Author: Daniel Rieger, KIT
  ! Initial Release: 2016-08-22
  ! Modifications:
  ! 2020-02-14: Sven Werchner, KIT
  ! - moved to mo_art_mode_fields.f90
  !>
    CLASS(t_mode_fields),INTENT(inout) :: &
      &  this_fields                        !< Container with fields
    TYPE(t_key_value_store),INTENT(in) :: &
      &  meta_key_value_store               !< Metadata container
  
    this_fields%linit = .TRUE.
  
  END SUBROUTINE set_meta_t_mode_fields
  !!
  !!-------------------------------------------------------------------------
  !!
  SUBROUTINE print_meta_t_mode_fields(this_fields)
  !<
  ! SUBROUTINE print_meta_t_mode_fields
  ! This subroutine prints the metadata contained in t_mode_fields
  ! Part of Module: mo_art_mode_fields (moved from mo_art_modes)
  ! Author: Daniel Rieger, KIT
  ! Initial Release: 2016-09-15
  ! Modifications:
  ! 2020-02-14: Sven Werchner, KIT
  ! - moved to mo_art_mode_fields.f90
  !>
    CLASS(t_mode_fields),INTENT(in)  :: &
      &  this_fields                      !< Fields of current mode
  
    IF(this_fields%linit) THEN
      WRITE (message_text,*) '==========ART: PRINTOUT MODEL METADATA FOR MODE '// &
        &                    TRIM(this_fields%name)//'=========='
      CALL message ('', message_text)
      WRITE (message_text,*) 'NAME            DATA'
      CALL message ('', message_text)
    ELSE
      CALL finish(TRIM(routine)//':print_meta_t_mode_fields',          &
        &         'Print not possible: Mode '//TRIM(this_fields%name)//' was not initialized.')
    ENDIF
  
  END SUBROUTINE print_meta_t_mode_fields
  !!
  !!-------------------------------------------------------------------------
  !!
END MODULE mo_art_mode_fields
