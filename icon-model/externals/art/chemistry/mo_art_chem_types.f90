!
! mo_art_chem_types
! This module provides data storage structures and constants for chem
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

MODULE mo_art_chem_types
! ICON 
  USE mo_kind,                          ONLY: wp
  USE mo_impl_constants,                ONLY: SUCCESS
  USE mo_exception,                     ONLY: finish
  USE mo_tracer_metadata_types,         ONLY: t_chem_meta
  USE mo_physical_constants,            ONLY: amd

! ART
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string
  
 

  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_chem_meta_param ! type of parametrised tracers
  PUBLIC :: t_chem_meta_passive
  PUBLIC :: t_chem_meta_mecca
  PUBLIC :: t_prod_list


  !#################################
  ! linked list for production rates
  !#################################

  TYPE t_prod_list
    REAL(wp) ::   &
      &  factor               !< stocheometric coefficient or branching coefficient
                              !  educt + other_educts => factor * this_tracer + others
    CLASS(t_chem_meta_param), POINTER :: &
      &  educt => NULL()      !< educt tracer, see above comment
    TYPE(t_prod_list), POINTER ::   &
      &  next_prod => NULL()  !< linked list of educts for this tracer
  END TYPE t_prod_list

  !###########################################################################
  ! type of chemical tracers in general
  !###########################################################################

  TYPE, extends(t_chem_meta) :: t_chem_meta_param
    REAL(wp) :: &
      &  mol_weight              !< molecular weight of the species (kg / mol)
    TYPE(t_prod_list), POINTER ::  &
      &  first_prod => NULL()    !< linked list of products
    REAL(wp), POINTER ::  &
      &  tracer(:,:,:)           !< tracer mass mixing ratio (kg/kg)
    REAL(wp), ALLOCATABLE ::  &
      &  prod(:,:,:)             !< 3-D production (sum of all educt)  ([tracer]/s)
    REAL(wp), ALLOCATABLE ::   &
      &  des_3d(:,:,:)  !< calculated 3-D destruction rate (s-1)

    CONTAINS

      PROCEDURE :: set_tracer => set_p_tracer_now
      PROCEDURE :: init_param => init_param_meta
      PROCEDURE :: append_prod => param_append_prod
      PROCEDURE :: get_prod => param_get_prod
      PROCEDURE :: add_prod => param_add_prod
      PROCEDURE :: convert_mmr_Nconc => param_convert_mmr_Nconc
      PROCEDURE :: convert_Nconc_mmr => param_convert_Nconc_mmr
      PROCEDURE :: convert_vmr_mmr => param_convert_vmr_mmr
      PROCEDURE :: convert_mmr_vmr => param_convert_mmr_vmr
  END TYPE t_chem_meta_param

  !###########################################################################
  ! type of passive chemical tracers
  !###########################################################################

  TYPE, EXTENDS(t_chem_meta) :: t_chem_meta_passive
    REAL(wp), POINTER ::  &
      &  tracer(:,:,:)           !< tracer mass mixing ratio (kg/kg)
    CONTAINS

      PROCEDURE :: set_tracer => set_p_tracer_now_passive
  END TYPE

  !###########################################################################
  ! type of MECCA tracers
  !###########################################################################

  TYPE, EXTENDS(t_chem_meta) :: t_chem_meta_mecca
    REAL(wp) :: &
      &  mol_weight
    INTEGER :: &
      &  number,  &
      &  init_mode

    CONTAINS

      PROCEDURE :: init => init_mecca_meta
      PROCEDURE :: convert_mmr_Nconc => mecca_convert_mmr_Nconc
      PROCEDURE :: convert_Nconc_mmr => mecca_convert_Nconc_mmr
      PROCEDURE :: convert_vmr_mmr => mecca_convert_vmr_mmr
      PROCEDURE :: convert_mmr_vmr => mecca_convert_mmr_vmr
  END TYPE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for t_chem_meta_param
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE set_p_tracer_now(this, p_tracer_now)
!<
! SUBROUTINE set_p_tracer_now
! This subroutine sets the current tracer value to the internal structure.
! This is needed because the tracer array depends on nnow and nnew
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter and Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_param),INTENT(inout) :: &
    &  this                        !< Container with fields
  REAL(wp), INTENT(IN), TARGET ::  &
    &  p_tracer_now(:,:,:)   !< tracer mass mixing ratio (kg/kg)
  
  this%tracer => p_tracer_now

END SUBROUTINE set_p_tracer_now


SUBROUTINE init_param_meta(this,nproma,nlev,nblks)
!<
! SUBROUTINE init_param_meta
! Initialisation of the arrays corresponding to t_chem_meta_param
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter and Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(inout) :: &
    &  this                   !< container with fields
  INTEGER, INTENT(in) :: &
    &  nproma,nlev,nblks  !< dimensions
  ! local variables
  INTEGER :: &
    &  ierror
  CHARACTER(LEN=IART_VARNAMELEN) :: &
    &  tracer_name
  CHARACTER(:), ALLOCATABLE      :: &
    &  c_tmp

  CALL this%opt_meta%get('mol_weight',this%mol_weight, ierror)

  IF (ierror /= SUCCESS) THEN
    CALL key_value_storage_as_string(this%opt_meta,'name',c_tmp)
    WRITE(tracer_name,'(A)') c_tmp
    CALL finish('mo_art_chem_types:init_param_meta', &
           &    'Element ''mol_weight'' missing for '//TRIM(tracer_name))
  END IF

  this%tracer => NULL()
  this%first_prod => NULL()
  IF (.NOT. ALLOCATED(this%prod)) ALLOCATE(this%prod(nproma,nlev,nblks))

  this%prod(:,:,:) = 0._wp
END SUBROUTINE init_param_meta



SUBROUTINE param_append_prod(this,factor,educt)
!<
! SUBROUTINE param_append_prod
! appending an element for the production list
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(inout) :: &
    &  this          !< container with fields
  REAL(wp), INTENT(in) :: &
    &  factor        !< factor as described above
  CLASS(t_chem_meta_param), INTENT(in), TARGET :: &
    &  educt         !< educt as described above
  ! local variables
  TYPE(t_prod_list), POINTER :: &
    &  current_prod  !< element of linked list

  IF (ASSOCIATED(this%first_prod)) THEN
    current_prod => this%first_prod

    DO WHILE (ASSOCIATED(current_prod%next_prod))
      current_prod => current_prod%next_prod
    END DO

    ALLOCATE(current_prod%next_prod)
    current_prod => current_prod%next_prod
  ELSE
    ALLOCATE(this%first_prod)
    current_prod => this%first_prod
  END IF

  current_prod%next_prod => NULL()

  current_prod%factor = factor
  current_prod%educt => educt

END SUBROUTINE param_append_prod


SUBROUTINE param_get_prod(this)
!<
! SUBROUTINE param_get_prod
! Calculates the three-dimensional production of this tracer
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(inout) :: &
    &  this            !< container with fields
  ! local variables
  TYPE(t_prod_list), POINTER :: &
    &  current_prod    !< element in prod list

  IF (ASSOCIATED(this%first_prod)) THEN
    this%prod = 0._wp

    current_prod => this%first_prod

    DO WHILE(ASSOCIATED(current_prod))
      this%prod = this%prod + current_prod%factor * current_prod%educt%des_3d  &
             &                 * current_prod%educt%tracer

      current_prod => current_prod%next_prod
    END DO
  END IF

END SUBROUTINE param_get_prod

SUBROUTINE param_add_prod(this, p_dtime)
!<
! SUBROUTINE param_add_prod
! This adds the production in case that it is not further parametrised
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_param),INTENT(inout) :: &
    &  this         !< Container with fields
  REAL(wp), INTENT(IN) :: &
    &  p_dtime      !< model time step

  this%tracer   = this%tracer + this%prod * p_dtime 

END SUBROUTINE param_add_prod

SUBROUTINE param_convert_mmr_Nconc(this, p_tracer_now, vmr2Nconc)
!<
! SUBROUTINE param_convert_mmr_Nconc
! This converts the tracer values from mass mixing ratio (kg/kg) 
! to number concentration (cm-3)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(inout) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)

  p_tracer_now = p_tracer_now * (amd * 1.e-3_wp) / this%mol_weight * vmr2Nconc

END SUBROUTINE param_convert_mmr_Nconc

SUBROUTINE param_convert_Nconc_mmr(this, p_tracer_now, vmr2Nconc)
!<
! SUBROUTINE param_convert_Nconc_mmr
! This converts the tracer values from number concentration (cm-3)
! to mass mixing ratio (kg/kg)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(inout) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)

  p_tracer_now = p_tracer_now / (amd * 1.e-3_wp) * this%mol_weight / vmr2Nconc

END SUBROUTINE param_convert_Nconc_mmr

SUBROUTINE param_convert_vmr_mmr(this, p_tracer_now)
!<
! SUBROUTINE param_convert_vmr_mmr
! This converts the tracer values from volume mixing ratio (mol/mol)
! to mass mixing ratio (kg/kg)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(in) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)

  p_tracer_now = p_tracer_now * this%mol_weight / (amd * 1.e-3_wp)

END SUBROUTINE param_convert_vmr_mmr

SUBROUTINE param_convert_mmr_vmr(this, p_tracer_now)
!<
! SUBROUTINE param_convert_mmr_vmr
! This converts the tracer values from mass mixing ratio (kg/kg)
! to volume mixing ratio (mol/mol)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_param), INTENT(in) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)

  p_tracer_now = p_tracer_now / this%mol_weight * (amd * 1.e-3_wp)

END SUBROUTINE param_convert_mmr_vmr

!#############################################################################################
! Passive tracer routines
!#############################################################################################


SUBROUTINE set_p_tracer_now_passive(this, p_tracer_now)
!<
! SUBROUTINE set_p_tracer_now_passive
! This subroutine sets the current tracer value to the internal structure.
! This is needed because the tracer array depends on nnow and nnew
! Part of Module: mo_art_chem_types
! Author: Jennifer Schroeter and Michael Weimer, KIT
! Initial Release: around 2018-10 
! Modifications:
!>
  CLASS(t_chem_meta_passive),INTENT(inout) :: &
    &  this                  !< Container with fields
  REAL(wp), INTENT(IN), TARGET ::  &
    &  p_tracer_now(:,:,:)   !< tracer value
  
  this%tracer => p_tracer_now

END SUBROUTINE set_p_tracer_now_passive
!#############################################################################################
! MECCA tracers routines
!#############################################################################################

SUBROUTINE init_mecca_meta(this)
!<
! SUBROUTINE init_mecca_meta
! This initialises the MECCA tracer meta information
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-16
! Modifications:
!>
  CLASS(t_chem_meta_mecca), INTENT(inout) :: &
    &  this
  INTEGER :: &
    &  ierror
  CHARACTER(LEN = IART_VARNAMELEN) :: &
    &  tracer_name
  CHARACTER(:), ALLOCATABLE  :: &
    &  c_tmp

  CALL this%opt_meta%get('number',this%number, ierror)

  IF (ierror /= SUCCESS) THEN
    CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
    WRITE(tracer_name,'(A)') c_tmp
    CALL finish('mo_art_chem_types:init_mecca_meta',  &
      &         'Element number missing which is required for the MECCA tracer ' &
      &       //TRIM(tracer_name)//'.')
  END IF

  CALL this%opt_meta%get('init_mode',this%init_mode, ierror)

  IF (ierror /= SUCCESS) THEN
    CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
    WRITE(tracer_name,'(A)') c_tmp
    CALL finish('mo_art_chem_types:init_mecca_meta',  &
      &         'Element init_mode missing which is required for the MECCA tracer ' &
      &       //TRIM(tracer_name)//'.')
  END IF

  CALL this%opt_meta%get('mol_weight',this%mol_weight, ierror)

  IF (ierror /= SUCCESS) THEN
    CALL key_value_storage_as_string(this%opt_meta,'name', c_tmp)
    WRITE(tracer_name,'(A)') c_tmp
    CALL finish('mo_art_chem_types:init_mecca_meta',  &
      &         'Element mol_weight missing which is required for the MECCA tracer ' &
      &       //TRIM(tracer_name)//'.')
  END IF

END SUBROUTINE init_mecca_meta

SUBROUTINE mecca_convert_mmr_Nconc(this, p_tracer_now, vmr2Nconc)
!<
! SUBROUTINE mecca_convert_mmr_Nconc
! This converts the tracer values from mass mixing ratio (kg/kg) 
! to number concentration (cm-3)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_mecca), INTENT(inout) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)

  p_tracer_now = p_tracer_now * (amd * 1.e-3_wp) / this%mol_weight * vmr2Nconc

END SUBROUTINE mecca_convert_mmr_Nconc

SUBROUTINE mecca_convert_Nconc_mmr(this, p_tracer_now, vmr2Nconc)
!<
! SUBROUTINE mecca_convert_Nconc_mmr
! This converts the tracer values from number concentration (cm-3)
! to mass mixing ratio (kg/kg)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_mecca), INTENT(inout) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)
  REAL(wp), INTENT(in) :: &
    &  vmr2Nconc(:,:,:)

  p_tracer_now = p_tracer_now / (amd * 1.e-3_wp) * this%mol_weight / vmr2Nconc

END SUBROUTINE mecca_convert_Nconc_mmr

SUBROUTINE mecca_convert_vmr_mmr(this, p_tracer_now)
!<
! SUBROUTINE mecca_convert_vmr_mmr
! This converts the tracer values from volume mixing ratio (mol/mol)
! to mass mixing ratio (kg/kg)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_mecca), INTENT(in) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)

  p_tracer_now = p_tracer_now * this%mol_weight / (amd * 1.e-3_wp)

END SUBROUTINE mecca_convert_vmr_mmr


SUBROUTINE mecca_convert_mmr_vmr(this, p_tracer_now)
!<
! SUBROUTINE mecca_convert_mmr_vmr
! This converts the tracer values from mass mixing ratio (kg/kg)
! to volume mixing ratio (mol/mol)
! Part of Module: mo_art_chem_types
! Author: Michael Weimer, KIT
! Initial Release: 2020-04-20
! Modifications:
!>
  CLASS(t_chem_meta_mecca), INTENT(in) :: &
    &  this
  REAL(wp), INTENT(inout) :: &
    &  p_tracer_now(:,:,:)

  p_tracer_now = p_tracer_now / this%mol_weight * (amd * 1.e-3_wp)

END SUBROUTINE mecca_convert_mmr_vmr


END MODULE mo_art_chem_types
  
