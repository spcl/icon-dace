!
! mo_art_chem_init_meta
! This module provides subroutines for initialising
! and cleaning up the tracer meta data
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

MODULE mo_art_chem_init_meta

  ! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: finish
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_fortran_tools,                 ONLY: init
  ! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_types_param,          ONLY: t_chem_meta_OH,     &
                                          &   t_chem_meta_lt,     &
                                          &   t_chem_meta_cold,   &
                                          &   t_chem_meta_linoz,  &
                                          &   t_chem_meta_simnoy, &
                                          &   simnoy_fill_init,   &
                                          &   linoz_fill_init
  USE mo_art_chem_types,                ONLY: t_chem_meta_param,   &
                                          &   t_chem_meta_passive, &
                                          &   t_chem_meta_mecca
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_OH_chem_types,             ONLY: t_art_OH_chem_meta
  USE mo_art_read_linoz,                ONLY: art_linoz_read
#ifdef __ART_GPL
  USE mo_art_psc_init,                  ONLY: art_psc_set_number_bins
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC ::                       &
    &  art_chem_init_meta,        &
    &  art_delete_OH_chem_meta


CONTAINS

!
!------------------------------------------------------------------------------------
!

SUBROUTINE art_chem_init_meta(jg,p_prog_list,OH_chem_meta,  &
                    &         lart_chemtracer, lart_mecca,  &
                    &         lart_psc)
!<
! SUBROUTINE art_chem_init_meta
! This subroutines goes through the list of tracers and checks if it takes part
! in the OH chemisty or is a lifetime tracer. Then, it creates some meta information
!
! Based on: 
! Part of Module: mo_art_chem_init_meta
! Author: Michael Weimer, KIT
! Initial Release: 2018-10-16
! Modifications:
!>
  IMPLICIT NONE
  INTEGER, INTENT(in)  ::   &
    &  jg                   !< patch on which computation is performed
  TYPE(t_var_list_ptr),INTENT(in) :: &
    &  p_prog_list          !< current list: prognostic
  TYPE(t_art_OH_chem_meta), INTENT(inout) ::  &
    &  OH_chem_meta         !< general meta data of OH chemistry
  LOGICAL, INTENT(in) :: &
    &  lart_chemtracer,  &  !< switch for chemtracers
    &  lart_psc,         &  !< switch for PSCs
    &  lart_mecca           !< switch for MECCA tracers

  !local variables
  CHARACTER(:), ALLOCATABLE :: &
    &  tracer_name         !< name of the tracer
  LOGICAL ::             &
    &  CH4_is_OH_tracer, & !< flag if CH4 is part of OH chemistry
    &  CO_is_OH_tracer     !< flag if CO  is part of OH chemistry
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo    !< ART atmo fields
  INTEGER :: &
    &  number_param_tracers,  & !< number of parametrised chemical tracers
    &  iv                       !< loop index 

  art_atmo => p_art_data(jg)%atmo

  IF (lart_psc) THEN
#ifdef __ART_GPL
    CALL art_psc_set_number_bins(p_art_data(jg)%chem%PSC_meta, jg)
#endif
  END IF

  CH4_is_OH_tracer = .FALSE.
  CO_is_OH_tracer = .FALSE.
  number_param_tracers = 0

  DO iv = 1, p_prog_list%p%nvars

    ! consistency check and setting the number of parametrised tracers
    SELECT TYPE(meta =>  p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      CLASS IS (t_chem_meta_param)
        IF (.NOT. lart_chemtracer) THEN
          CALL finish('mo_art_chem_init_meta:art_chem_init_meta',        &
                 &    'Found chemtracer '//TRIM(meta%name)//' although ' &
                 &  //'lart_chemtracer is set to .FALSE.')
        END IF

        number_param_tracers = number_param_tracers + 1
      CLASS IS (t_chem_meta_passive)
        IF (.NOT. lart_chemtracer) THEN
          CALL finish('mo_art_chem_init_meta:art_chem_init_meta',      &
                 &    'Found passive '//TRIM(meta%name)//' although '  &
                 &  //'lart_chemtracer is set to .FALSE.')
        END IF
      CLASS IS (t_chem_meta_mecca)
        IF (.NOT. lart_mecca) THEN
          CALL finish('mo_art_chem_init_meta:art_chem_init_meta',          &
                 &    'Found MECCA tracer '//TRIM(meta%name)//' although ' &
                 &  //'lart_mecca is set to .FALSE.')
        END IF
    END SELECT


    SELECT TYPE(meta => p_prog_list%p%vl(iv)%p%info_dyn%tracer)
      TYPE IS (t_chem_meta_OH)
        CALL meta%get_tracer_name(tracer_name)
        CALL meta%init(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks)

        IF (.NOT. (OH_chem_meta%is_init)) THEN
          CALL art_init_OH_chem(jg,OH_chem_meta)
        END IF

        ! this part below just checks if the tracers CH4 and CO are included in OH
        ! chemistry. Otherwise an error is raised
        IF ('TRCH4' == TRIM(tracer_name)) THEN
          CH4_is_OH_tracer = .TRUE.
        ELSE IF ('TRCO' == TRIM(tracer_name)) THEN
          CO_is_OH_tracer = .TRUE.
        END IF


      TYPE IS (t_chem_meta_lt)
        CALL meta%init(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks)
        CALL meta%get_destruct

      TYPE IS (t_chem_meta_cold)
        CALL meta%init(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks)

      TYPE IS (t_chem_meta_linoz)
        IF (.NOT. p_art_data(jg)%chem%param%linoz%is_init) THEN
          CALL art_linoz_read(p_art_data(jg)%chem%param%linoz,  &
                     &        art_atmo%nproma,                  &
                     &        art_atmo%nlev)
        END IF

        CALL meta%init(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks)
        CALL linoz_fill_init(meta,jg,p_prog_list%p%vl(iv)%p%info%ncontained)

      TYPE IS (t_chem_meta_simnoy)
        CALL meta%init(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks)
        CALL simnoy_fill_init(meta,jg,p_prog_list)

      TYPE IS (t_chem_meta_param)
        CALL meta%init_param(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks)

      TYPE IS (t_chem_meta_mecca)
        CALL meta%init()

    END SELECT

  END DO

  IF (OH_chem_meta%is_init) THEN
    IF (.NOT. ((CH4_is_OH_tracer) .AND. (CO_is_OH_tracer) )) THEN
      CALL finish('mo_art_chem_init_meta:art_chem_init_meta',                       &
       &          'Please include <c_solve>OH</c_solve> element for TRCO and TRCH4' &
       &        //' in the XML file if you want to use the simplified OH '          &
       &        //'chemistry with other tracers.')
    END IF
  END IF

  p_art_data(jg)%chem%param%number_param_tracers = number_param_tracers

END SUBROUTINE art_chem_init_meta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE art_init_OH_chem(jg,OH_chem_meta)
!<
! SUBROUTINE art_init_OH_chem
! This subroutine appends a list element to the previously created sublist of tracers for
! predictor-corrector method, and includes some persistent metadata.
! Based on: 
! Part of Module: mo_art_chem_init_meta
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
! 2016-11-16: Michael Weimer, KIT
! - moved reaction rates to mo_art_predcorr_method
!>
  IMPLICIT NONE
  INTEGER, INTENT(in) ::  &
    &  jg           !< patch id
  TYPE(t_art_OH_chem_meta), INTENT(inout) ::   &
    &  OH_chem_meta !< general meta data of OH chemistry

  !local variables
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo     !< ART atmo fields

  art_atmo => p_art_data(jg)%atmo

  ! oh_nconc already handled by add_var. 

  OH_chem_meta%num_elements = 0

  ALLOCATE(OH_chem_meta%ozone_Nconc(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ALLOCATE(OH_chem_meta%k_CH4_OH(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ALLOCATE(OH_chem_meta%k_CO_OH(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ALLOCATE(OH_chem_meta%level_CH4_gt_1ppm(art_atmo%nproma,art_atmo%nblks))


!  CALL init(OH_chem_meta%OH_Nconc)
  CALL init(OH_chem_meta%ozone_Nconc, lacc=.FALSE.)
  CALL init(OH_chem_meta%k_CH4_OH, lacc=.FALSE.)
  CALL init(OH_chem_meta%k_CO_OH, lacc=.FALSE.)
  OH_chem_meta%level_CH4_gt_1ppm(:,:) = art_atmo%nlev

  OH_chem_meta%is_init = .TRUE.

END SUBROUTINE art_init_OH_chem

!
!----------------------------------------------------------------------------------------
!

SUBROUTINE art_delete_OH_chem_meta(OH_chem_meta)
!<
! SUBROUTINE art_delete_OH_chem_meta
! This subroutine deletes the meta data of the OH chemistry
! Based on: 
! Part of Module: mo_art_chem_init_meta
! Author: Michael Weimer, KIT
! Initial Release: 2016-04-05
! Modifications:
!>
  IMPLICIT NONE
  TYPE(t_art_OH_chem_meta), INTENT(inout) ::   &
    &  OH_chem_meta           !< general meta data of OH chemistry

  IF (ALLOCATED(OH_chem_meta%ozone_Nconc))  DEALLOCATE(OH_chem_meta%ozone_Nconc)
  IF (ALLOCATED(OH_chem_meta%k_CH4_OH))     DEALLOCATE(OH_chem_meta%k_CH4_OH)
  IF (ALLOCATED(OH_chem_meta%k_CO_OH))      DEALLOCATE(OH_chem_meta%k_CO_OH)
  IF (ALLOCATED(OH_chem_meta%level_CH4_gt_1ppm))  &
           &    DEALLOCATE(OH_chem_meta%level_CH4_gt_1ppm)

END SUBROUTINE art_delete_OH_chem_meta

!!
!!-----------------------------------------------------------------------------
!!
END MODULE mo_art_chem_init_meta
