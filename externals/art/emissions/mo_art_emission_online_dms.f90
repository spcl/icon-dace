!
! mo_art_emission_online_dms
! This module provides the allocation and calculation of the emission flux
! density for online dms.
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

MODULE mo_art_emission_online_dms
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_tracer_metadata_types,         ONLY: t_chem_meta
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic
  USE mo_impl_constants,                ONLY: SUCCESS

  USE mtime,                            ONLY: datetime

! ART
  USE mo_art_chem_types_param,          ONLY: t_chem_meta_lt
  USE mo_art_impl_constants,            ONLY: IART_VARNAMELEN
  USE mo_art_external_types,            ONLY: t_art_online_dms
  USE mo_art_read_emissions,            ONLY: art_convert_emission_to_mmr
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_string_tools,              ONLY: key_value_storage_as_string
  IMPLICIT NONE
    
  PRIVATE
  
  PUBLIC :: art_extinit_dms
  PUBLIC :: art_add_dms_emission_to_tracer

  CHARACTER(LEN=*), PARAMETER :: routine = 'mo_art_emission_online_dms'

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_dms(p_online_dms, jg, tracer_name_DMS, p_prog_list)
!<
! SUBROUTINE art_extinit_dms
! This SR allocates an array for online dms.
! Part of Module: mo_art_emission_online_dms
! Author: Carmen Ullwer, KIT
! Initial Release: 2018-01-09
! Modifications:
! YYYY-MM-DD: <author>, <institution>
! - <Description>
!>
  TYPE(t_art_online_dms), INTENT(INOUT) :: &
    &  p_online_dms           !< online dms
  INTEGER, INTENT(IN) ::    &
    &  jg              !< patch id
  CHARACTER(LEN=*), INTENT(IN) :: &
    &  tracer_name_DMS
  TYPE(t_var_list_ptr), INTENT(IN)       :: &
    &  p_prog_list            !< current prognostic state list
  ! local variables
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn 
  INTEGER :: ierror, iv
  CHARACTER(LEN=IART_VARNAMELEN)       :: &
    &  calc_onl,        &     !< switch for online calculation 
    &  tracer_name            !< name of the tracer
  CHARACTER(:), ALLOCATABLE :: &
    &  c_tmp
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo


  art_atmo  => p_art_data(jg)%atmo

  ! ----------------------------------
  ! --- Start DO-loop over elements in list
  ! ----------------------------------

  DO iv = 1, p_prog_list%p%nvars
    
    ! ----------------------------------
    ! --- Get meta data of current element and assure that current element is
    ! tracer
    ! ----------------------------------
  
    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
  
  
    IF (info_dyn%tracer%lis_tracer) THEN
  
      SELECT TYPE(meta => info_dyn%tracer)
        CLASS IS (t_chem_meta)
  
          CALL key_value_storage_as_string(meta%opt_meta,'name', c_tmp)
          WRITE(tracer_name,'(A)') c_tmp

          IF (TRIM(ADJUSTL(tracer_name)) == TRIM(ADJUSTL(tracer_name_DMS))) THEN
            ! ----------------------------------
            ! --- Check for XML parameter 'online' (DMS only)
            ! ----------------------------------
  
            CALL key_value_storage_as_string(meta%opt_meta,'online', c_tmp, ierror)
            IF (ierror== SUCCESS) THEN
              WRITE(calc_onl,'(A)') c_tmp
            ELSE 
              calc_onl = 'OFF'
            END IF

            IF (TRIM(ADJUSTL(calc_onl)) == 'ON') THEN 
              p_online_dms%lcalc_onl = .TRUE.

              ! it could happen that this is called for chemtracer and MECCA
              IF (.NOT. ALLOCATED(p_online_dms%dms_month)) THEN
                ALLOCATE( p_online_dms%dms_month(art_atmo%nproma,    &
                                       &         art_atmo%nblks,     &
                                       &         p_online_dms%ndms_months) )
              END IF

              IF (.NOT. ALLOCATED(p_online_dms%mmr_onl_dms)) THEN
                ALLOCATE( p_online_dms%mmr_onl_dms(art_atmo%nlev))
              END IF
            ENDIF
          END IF
      END SELECT !type meta
    ENDIF
  ENDDO

END SUBROUTINE art_extinit_dms
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_add_dms_emission_to_tracer(jg,jb,jcs,jce ,&
            &             tracer,current_date, p_dtime)
!<
! SUBROUTINE art_add_dms_emission_to_tracer
! This SR calculates the emission flux of online DMS and adds it to the
! corresponding DMS tracer
! Part of Module: mo_art_emission_online_dms
! Author: Michael Weimer, KIT
! Initial Release: 2018-06-26
! Modifications:
! YYYY-MM-DD: <author>, <institution>
! - <Description>
!>
  IMPLICIT NONE
  INTEGER        ,INTENT(IN) :: jg,jb,jcs,jce
  TYPE(t_chem_meta_lt), INTENT(INOUT) :: tracer
  TYPE(datetime), POINTER, INTENT(IN)  ::  current_date       !< current date
  REAL(wp), INTENT(IN) :: p_dtime
  ! local variables
  INTEGER :: &
    &  i, nlev               !< loop index
  REAL(wp) ::         &
    &  emiss_onl_dms, &      !< emission flux of DMS (kg m-2 s-1)
    &  vb                    !< absolute value of wind speed (m/s)
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  TYPE(t_art_online_dms), POINTER :: &
    &  p_online_dms           !< online dms

  art_atmo     => p_art_data(jg)%atmo
  p_online_dms => p_art_data(jg)%ext%online_dms

  nlev = art_atmo%nlev
  
  DO i = jcs,jce
    IF (.NOT. art_atmo%llsm(i,jb)) THEN
      vb = SQRT(art_atmo%u_10m(i,jb) * art_atmo%u_10m(i,jb)   &
            & + art_atmo%v_10m(i,jb) * art_atmo%v_10m(i,jb))
    
      emiss_onl_dms = (tracer%mol_weight                                        &
           &           * p_online_dms%dms_month(i,jb,current_date%date%month)   &
           &           * (0.222_wp * vb*vb + 0.3334_wp * vb  )    )             &
           &           / (3.6e2_wp * 1e9_wp)

      CALL art_convert_emission_to_mmr(p_online_dms%mmr_onl_dms(:), &
           &        emiss_onl_dms,                                  &
           &        art_atmo%temp(i,:,jb),                          &
           &        art_atmo%pres(i,:,jb),                          &
           &        art_atmo%dz(i,:,jb),                            &
           &        p_dtime,1,nlev)

      tracer%tracer(i,nlev,jb) = tracer%tracer(i,nlev,jb) + p_online_dms%mmr_onl_dms(nlev)
    ENDIF
  END DO

END SUBROUTINE art_add_dms_emission_to_tracer
!!
!!-------------------------------------------------------------------------
!!
  
END MODULE mo_art_emission_online_dms
