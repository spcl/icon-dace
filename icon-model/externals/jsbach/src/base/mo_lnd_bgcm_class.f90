!> define indices and helper functions for bgc material
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define indices (enumerator) and helper functions for bgc material and their components
!>
MODULE mo_lnd_bgcm_class
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_util,                ONLY: int2string
  USE mo_jsb_pool_class,      ONLY: t_jsb_pool
  USE mo_exception,           ONLY: message, message_text, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, FIRST_ELEM_ID, LAST_ELEM_ID

  PUBLIC :: get_name_for_element_id, get_max_number_of_elements, get_element_index, get_short_name_for_element_id
  PUBLIC :: assert_no_compartment_bgc_materials, get_dimension_of_element_variables

  !> IDs of elements (since elements are dynamic these can only be used as idx in arrays with all possible elements)
  ENUM, BIND(C)
    ENUMERATOR ::           &
      & ELEM_C_ID = 1,      &
      & ELEM_N_ID,          &
      & ELEM_P_ID,          &
      & ELEM_C13_ID,        &
      & ELEM_C14_ID,        &
      & ELEM_N15_ID,        &
      & ONEPLUS_LAST_ELEM_ID  ! needs to be the last -- it is only used to determine the max number of elements
  END ENUM

  !> ID of first and last bgcm elements
  INTEGER, PARAMETER :: FIRST_ELEM_ID = ELEM_C_ID
  INTEGER, PARAMETER :: LAST_ELEM_ID  = ONEPLUS_LAST_ELEM_ID - 1

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_bgcm_class'

CONTAINS

  ! ====================================================================================================== !
  !>
  !> Returns the number of elements defined in the enumerator
  !>
  FUNCTION get_max_number_of_elements() RESULT(last_type_id)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: last_type_id
    ! -------------------------------------------------------------------------------------------------- !
    ! The elem id enum  starts with 1 to be able to use the ids as indices in arrays with all possible elements
    last_type_id = LAST_ELEM_ID

  END FUNCTION get_max_number_of_elements


  ! ======================================================================================================= !
  !>
  !> Helper: assert that bgc material with elements does not contain further bgc materials
  !>
  SUBROUTINE assert_no_compartment_bgc_materials(this, tile_name, routine)
    USE mo_exception,           ONLY: finish
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool), POINTER, INTENT(IN) :: this      !< bgc material to test assertion for
    CHARACTER(len=*),           INTENT(IN) :: tile_name !< name of the tile with this bgc material
    CHARACTER(len=*),           INTENT(IN) :: routine   !< name of the calling routine
    ! -------------------------------------------------------------------------------------------------- !
    IF(ASSOCIATED(this%pool_list)) THEN
      IF(.NOT. SIZE(this%pool_list) .EQ. 0) THEN
        CALL finish(TRIM(routine), 'Found bgcm '//TRIM(this%name)//' in tile '//TRIM(tile_name)// &
          & ' with elements and compartments, but no behaviour is implemented for this so far, please check!')
      END IF
    END IF

  END SUBROUTINE assert_no_compartment_bgc_materials

  ! ======================================================================================================= !
  !>
  !> returns the dimension of the element variables contained in the given bgc material
  !> asserts that it is a dimension which can be handled by the current implementation
  !> and asserts that all elements have the same dimension
  !>
  FUNCTION get_dimension_of_element_variables(this, tile_name, routine) RESULT(this_dim)
    USE mo_exception,           ONLY: finish
    USE mo_jsb_var_class,       ONLY: REAL2D, REAL3D
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool), POINTER, INTENT(IN) :: this      !< bgc material to test assertion for
    CHARACTER(len=*),           INTENT(IN) :: tile_name !< name of the tile with this bgc material
    CHARACTER(len=*),           INTENT(IN) :: routine   !< name of the calling routine
    INTEGER                                :: this_dim  !< dimension
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: nr_of_elements, i_element
    ! -------------------------------------------------------------------------------------------------- !
    nr_of_elements = SIZE(this%element_list)
    this_dim = -1

    ! assert that all variables have the same dimensionality
    DO i_element = 1, nr_of_elements
      IF(this_dim == -1) THEN
        this_dim = this%element_list(i_element)%p%type
      ELSEIF(this_dim /= this%element_list(i_element)%p%type) THEN
        CALL finish(TRIM(routine), 'Found elements with different dimensions for bgcm '//TRIM(this%name)// &
          & ' in tile '//TRIM(tile_name)//', but no behaviour is implemented for this so far, please check!')
      ENDIF
    ENDDO

    IF(this_dim /= REAL2D .AND. this_dim /= REAL3D) THEN
      CALL finish(TRIM(routine), 'Found elements with other than 2D or 3D type for bgcm '//TRIM(this%name)// &
        & ' in tile '//TRIM(tile_name)//' for this so far no behaviour is implemented, please check!')
    ENDIF

  END FUNCTION get_dimension_of_element_variables


  ! ======================================================================================================= !
  !>
  !> get the short name for this element ID
  !>
  FUNCTION get_short_name_for_element_id(id_element) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in)           :: id_element
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id_element)
    CASE (ELEM_C_ID)
      return_value = 'C'
    CASE (ELEM_N_ID)
      return_value = 'N'
    CASE (ELEM_P_ID)
      return_value = 'P'
    CASE (ELEM_C13_ID)
      return_value = 'C13'
    CASE (ELEM_C14_ID)
      return_value = 'C14'
    CASE (ELEM_N15_ID)
      return_value = 'N15'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION get_short_name_for_element_id


  ! ======================================================================================================= !
  !>
  !> get the name for this element ID
  !>
  FUNCTION get_name_for_element_id(id_element) RESULT(return_value)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: id_element
    CHARACTER(len=:), ALLOCATABLE :: return_value
    ! -------------------------------------------------------------------------------------------------- !

    SELECT CASE(id_element)
    CASE (ELEM_C_ID)
      return_value = 'carbon'
    CASE (ELEM_N_ID)
      return_value = 'nitrogen'
    CASE (ELEM_P_ID)
      return_value = 'phosphorus'
    CASE (ELEM_C13_ID)
      return_value = 'carbon13'
    CASE (ELEM_C14_ID)
      return_value = 'carbon14'
    CASE (ELEM_N15_ID)
      return_value = 'nitrogen15'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION get_name_for_element_id


  ! ======================================================================================================= !
  !>
  !> Returns (and checks!) the element index in the index map for the queried id
  !>
  FUNCTION get_element_index(element_id, elements_index_map, caller) RESULT(element_idx)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,          INTENT(in) :: element_id            !< id of the queried element
    INTEGER,          INTENT(in) :: elements_index_map(:) !< id-index map
    CHARACTER(len=*), INTENT(in) :: caller                !< calling routine
    INTEGER                      :: element_idx           !< index of the queried element
    ! -------------------------------------------------------------------------------------------------- !

    element_idx = elements_index_map(element_id)

    ! Assert that the model is running with this element
    IF ( element_idx == 0 ) THEN
      CALL finish(TRIM(caller), 'Unknown or unused element id: ' // int2string(element_id) &
        & // ' - element name: "' // get_name_for_element_id(element_id) // '". Please check!' )
    ENDIF
  END FUNCTION get_element_index

#endif
END MODULE mo_lnd_bgcm_class
