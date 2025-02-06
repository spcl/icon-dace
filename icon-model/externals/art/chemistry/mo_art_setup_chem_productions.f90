!
! mo_art_string_tools
! This module provides routines for parsing arrays and production elements from
! XML strings
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

MODULE mo_art_setup_chem_productions
  ! ICON
  USE mo_kind,                       ONLY: wp
  USE mo_exception,                  ONLY: finish
  USE mo_var_list,                   ONLY: t_var_list_ptr
  USE mo_var,                        ONLY: t_var
  USE mo_impl_constants,             ONLY: SUCCESS
  USE mo_key_value_store,            ONLY: t_key_value_store
  ! ART
  USE mo_art_string_tools,           ONLY: art_parse_production, &
                                       &   key_value_storage_as_string
  USE mo_art_chem_types,             ONLY: t_chem_meta_param
  USE mo_art_impl_constants,         ONLY: IART_VARNAMELEN

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_setup_chem_productions
CONTAINS
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_setup_chem_productions(p_prog_list, dict_tracer)
!<
! SUBROUTINE art_setup_chem_productions
! setup the production lists of each tracer based on the XML element products
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  TYPE(t_var_list_ptr), INTENT(inout) :: &
    &  p_prog_list
  TYPE(t_key_value_store), INTENT(in) :: &
    &  dict_tracer

  ! local variables
  INTEGER ::              &
    &  ierror, iv
  REAL(wp), ALLOCATABLE :: &
    &  factors(:)
  CHARACTER(:), ALLOCATABLE :: &
    &  products_str
  CHARACTER(LEN = IART_VARNAMELEN), ALLOCATABLE :: &
    &  tracer_names(:)

  DO iv = 1, p_prog_list%p%nvars
    SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      CLASS IS (t_chem_meta_param)
        CALL key_value_storage_as_string(meta%opt_meta,'products',products_str, ierror)
 
        IF (ierror == SUCCESS) THEN
          CALL art_parse_production(factors,tracer_names,TRIM(products_str))
          CALL art_setup_product_list(meta,p_prog_list,factors,tracer_names,dict_tracer)
          DEALLOCATE(factors)
          DEALLOCATE(tracer_names)
        END IF
    END SELECT
  END DO

END SUBROUTINE art_setup_chem_productions
!!
!!-------------------------------------------------------------------
!!
SUBROUTINE art_setup_product_list(educt,p_prog_list,factors,tracer_names,dict_tracer)
!<
! SUBROUTINE art_setup_product_list
! appends all the elements from XML products to the respective tracers
! Based on: -
! Part of Module: mo_art_string_tools
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-18
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE
  CLASS(t_chem_meta_param), INTENT(in), TARGET :: &
    &  educt
  TYPE(t_var_list_ptr), INTENT(inout) :: &
    &  p_prog_list
  REAL(wp), INTENT(in) :: &
    &  factors(:)
  CHARACTER(LEN=*), INTENT(in) :: &
    &  tracer_names(:)
  TYPE(t_key_value_store), INTENT(in) :: &
    &  dict_tracer
  ! local variables
  INTEGER ::        &
    &  num_elem, i, &
    &  ierror, iTR, iv
  CHARACTER(:), ALLOCATABLE :: &
    &  educt_name,             &
    &  product_name

  num_elem = SIZE(factors)

  DO i = 1,num_elem
    CALL dict_tracer%get(tracer_names(i),iTR,ierror)
    IF (ierror /= SUCCESS) THEN
      CALL key_value_storage_as_string(educt%opt_meta,'name',educt_name)

      CALL finish('mo_art_setup_chem_productions:art_setup_product_list',  &
              &   'Tracer '//TRIM(tracer_names(i))                         &
              & //' not found which is in products of '//TRIM(educt_name))
    END IF


    DO iv = 1, p_prog_list%p%nvars
      SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
        CLASS IS (t_chem_meta_param)
          CALL key_value_storage_as_string(meta%opt_meta,'name',product_name)

          IF (TRIM(product_name) == TRIM(ADJUSTL(tracer_names(i)))) THEN
            CALL meta%append_prod(factors(i),educt)
            EXIT
          END IF
      END SELECT
    END DO
  END DO

END SUBROUTINE art_setup_product_list
!!
!!-------------------------------------------------------------------
!!
END MODULE mo_art_setup_chem_productions
