!
! mo_art_emission_full_chemistry
! This module provides the emission subroutines for chemical tracer processes.
! This Module is only envoked if lart_mecca == .TRUE.
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

MODULE mo_art_emission_full_chemistry
    USE mo_kind,                     ONLY: wp
    USE mo_var_list,                 ONLY: t_var_list_ptr
    USE mo_var,                      ONLY: t_var
    USE mo_var_metadata_types,       ONLY: t_var_metadata_dynamic, t_var_metadata

    USE mtime,                       ONLY: datetime

    ! ART
    USE mo_art_cell_loop,            ONLY: art_loop_cell_tracer_td
    USE mo_art_chem_types_param,     ONLY: t_chem_meta_lt
    USE mo_art_external_types,       ONLY: t_art_online_dms
    USE mo_art_data,                 ONLY: p_art_data
    USE mo_art_atmo_data,            ONLY: t_art_atmo
    USE mo_art_emission_online_dms,  ONLY: art_add_dms_emission_to_tracer
    USE messy_mecca_kpp_parameters,  ONLY: ind_DMS
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC   ::                        &
      &   art_emiss_full_chemistry
         
CONTAINS

SUBROUTINE art_emiss_full_chemistry(current_date,dtime,p_tracer_now,jg,p_prog_list)

! SUBROUTINE art_emiss_full_chemistry
! This subroutine actually does nothing at the current stage
! Part of Module: mo_art_emission_full_chemistry
! Author: Jennifer Schroeter, KIT
! Initial Release: ????
! Modifications:
! 2018-06-27: Michael Weimer, KIT
! - Added online DMS emissions
!>
  IMPLICIT NONE

  INTEGER, INTENT(IN)               :: jg

  REAL(wp), INTENT(INOUT) ::  &  !< tracer mixing ratios (specific concentrations)
    &  p_tracer_now(:,:,:,:)     !< at current time level n (before transport)
                                 !< [kg/kg]
                                 !< dim: (nproma,nlev,nblks_c,ntracer)

  REAL(wp), INTENT(IN) :: dtime  ! time step    

  TYPE(datetime), POINTER, INTENT(IN) :: current_date
  TYPE(t_var_list_ptr), INTENT(INOUT)   ::  &
    &  p_prog_list         !< current prognostic state list

  !local variables
  TYPE(t_var_metadata_dynamic), POINTER ::  &
    &  info_dyn         !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata), POINTER ::  &
    &  info             !< returns reference to tracer  metadata of current element
  INTEGER, POINTER ::   &
    &  jsp               !< returns index of element
  INTEGER          ::   &
    &  iv

  TYPE(t_art_online_dms), POINTER :: &
    &  p_online_dms
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  INTEGER, POINTER :: &
    &  mapping_indices_kpp(:)

  !!------------------------------------------------------------------------------
  !! --- Allocation
  !!------------------------------------------------------------------------------

  p_online_dms => p_art_data(jg)%ext%online_dms
  art_atmo => p_art_data(jg)%atmo
  mapping_indices_kpp => p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp

  DO iv = 1, p_prog_list%p%nvars

    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    IF (info_dyn%tracer%lis_tracer) THEN

      jsp=>info%ncontained

      SELECT TYPE(tracer => info_dyn%tracer)
        CLASS IS (t_chem_meta_lt)
          ! Online DMS emissions
          IF (ind_DMS /= 0) THEN
            IF ((jsp == mapping_indices_kpp(ind_DMS)) .AND. (p_online_dms%lcalc_onl)) THEN
              CALL tracer%set_tracer(p_tracer_now(:,:,:,jsp))
              CALL art_loop_cell_tracer_td(jg, tracer, current_date,   &
                            &              dtime, art_add_dms_emission_to_tracer)
            END IF
          END IF
      END SELECT
    END IF

  END DO

  NULLIFY(art_atmo)
END SUBROUTINE art_emiss_full_chemistry
!
! -----------------------------------------------------------
!
END MODULE mo_art_emission_full_chemistry


