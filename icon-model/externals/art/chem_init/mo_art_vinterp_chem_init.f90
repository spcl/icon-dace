!
! mo_art_vinterp_chem_init
! This module provides interpolation of external input data onto ICON grid
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

#include "omp_definitions.inc"

MODULE mo_art_vinterp_chem_init
  ! ICON
  USE mo_kind,                         ONLY: wp
  USE mo_run_config,                   ONLY: ntracer
  USE mo_model_domain,                 ONLY: t_patch
  USE mo_exception,                    ONLY: message,finish
  USE mo_io_units,                     ONLY: filename_max
  USE mo_impl_constants,               ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_var_list,                     ONLY: t_var_list_ptr
  USE mo_var,                          ONLY: t_var
  USE mo_var_metadata_types,           ONLY: t_var_metadata, t_var_metadata_dynamic
  USE mo_nh_vert_interp,               ONLY: prepare_lin_intp,prepare_extrap, lin_intp
  USE mo_read_interface,               ONLY: t_stream_id,  openInputFile, closeFile, &
    &                                        read_3d_1time, on_cells
  USE mo_io_config,                    ONLY: default_read_method
  USE mo_tracer_metadata_types,        ONLY: t_chem_meta
  ! ART
  USE mo_art_data,                     ONLY: p_art_data
  USE mo_art_atmo_data,                ONLY: t_art_atmo
  USE mo_art_chem_init_types,          ONLY: t_chem_init_state
  USE mo_art_string_tools,             ONLY: key_value_storage_as_string
 

  IMPLICIT NONE  
  

  INCLUDE 'netcdf.inc'


  PRIVATE


  PUBLIC  ::  art_vinterp_chem_init 
  PUBLIC  ::  art_read_in_chem_init

 

CONTAINS





SUBROUTINE art_read_in_chem_init(chem_init,p_patch, p_prog_list, &
                    &            chem_init_path)
!<
! SUBROUTINE art_read_in_chem_init                   
! This subroutine reads in external data set and stores ist
! Part of Module: mo_art_vinterp_chem_init
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-11-15                
! Modifications:
! 2017-06-09: Michael Weimer, KIT
! - included chem_init as parameter
!>

  TYPE(t_patch), INTENT(in)               ::  &    !< patch on which computation
    &  p_patch                                     !< is performed
  CHARACTER(LEN=filename_max), INTENT(IN) ::  &
    &  chem_init_path                                !< path of the external file
  TYPE(t_var_list_ptr), INTENT(in)       :: &
    &  p_prog_list      !< current prognostic state list
  TYPE(t_chem_init_state), INTENT(inout)    :: &
    &  chem_init                    !< Pointer to ART diagnostic fields

  ! local variables
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER :: &
    &  routine = 'art_read_in_chem_init'
  INTEGER ::      &
    &  init_mode, & !< mode for initialisation of the tracer (from xml)
    &  ierror,    & !< flag if meta storage element exists
    &  iv,        & !< loop variable
    &  jg           !< patch id
  CHARACTER(:), ALLOCATABLE          ::  &
    &  init_name    !< variable name in external intialisation file
  TYPE(t_var_metadata), POINTER      :: & 
    &  info                     !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata_dynamic), POINTER :: & 
    &  info_dyn                 !< returns reference to tracer metadata of current element 
  INTEGER, POINTER ::   &
    &  jsp                      !< returns index of element
  TYPE(t_stream_id) ::  &
    &  stream_id                !< id of the netCDF file
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  jg = p_patch%id

  art_atmo => p_art_data(jg)%atmo

  chem_init%chem_init_chem%spec(:,:,:,:) = 0.0_wp
  chem_init%chem_init_chem%spec_interp(:,:,:,:) = 0.0_wp

  CALL openInputFile(stream_id, TRIM(chem_init_path), p_patch, &
    &                default_read_method)
  ! first get species names that should be initialized:

  DO iv = 1,p_prog_list%p%nvars

    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
    info=>p_prog_list%p%vl(iv)%p%info

    IF (info_dyn%tracer%lis_tracer) THEN

      jsp=>info%ncontained
      SELECT TYPE(meta=>info_dyn%tracer)
        CLASS IS (t_chem_meta)
          CALL meta%opt_meta%get('init_mode',init_mode,ierror)

          IF (ierror/= SUCCESS) CALL finish(TRIM(routine)//':art_read_in_chem_init',       &
                                    &       'Metadata init_mode not available for tracer ' &
                                    &     //TRIM(meta%name)//'.')
          IF (init_mode == 1) THEN

            CALL key_value_storage_as_string(meta%opt_meta,'init_name',init_name,ierror)
                IF (ierror/= SUCCESS) CALL finish(TRIM(routine)//':art_read_in_chem_init',   &
                                    &         'Metadata init_name not available for tracer ' &
                                    &        //TRIM(meta%name)//'.')

            CALL message('mo_art_vinterp ', 'reading in ' // TRIM(init_name) // ' ...' )
            CALL read_3d_1time(stream_id,                                            &
                       &       on_cells,                                             &
                       &       TRIM(init_name),                                      &
                       &       fill_array= chem_init%chem_init_chem%spec(:,:,:,jsp))

         ENDIF

     END SELECT
    ENDIF

  ENDDO 

  CALL closeFile(stream_id)
  
END SUBROUTINE art_read_in_chem_init

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

SUBROUTINE  art_vinterp_chem_init(chem_init,num_vars, jg, p_prog_list)
!<
! SUBROUTINE art_vinterp_chem_init                  
! This subroutine interpolates MOZART dataset
! Part of Module: mo_art_vinterp_chem_init
! Author: Jennifer Schroeter, KIT
! Initial Release: 2015-11-15                
! Modifications:
! 2017-06-09: Michael Weimer, KIT
! - included chem_init as parameter
!>

  INTEGER, INTENT(in)   ::  &
    &  jg                      !< patch id
  TYPE(t_chem_init_state),INTENT(inout)     ::  &
    & chem_init                !< contains all meta data for vertical interpolation of external data
  INTEGER, INTENT(in) :: &
    &  num_vars                !< number of variables in the file
  TYPE(t_var_list_ptr), INTENT(in),optional       :: &
    &  p_prog_list             !< current prognostic state list

  ! local variables
  CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
     &    routine = 'art_vinterp_chem_init'

  INTEGER  ::     &
    &  no_levels, &    !< number of levels in the external data
    &  i,iv            !< loop indices

  REAL(wp), ALLOCATABLE :: &
    &  wfac_lin(:,:,:)    !< weighting factor of upper level
  INTEGER , ALLOCATABLE ::   &
    &  idx0_lin(:,:,:),      &  !< index of upper level
    &  bot_idx_lin(:,:),     &  !< index of lowest level for which interpolation is possible
    &  kpbl1(:,:), kpbl2(:,:)   !< Indices of model levels lying immediately above
  REAL(wp), ALLOCATABLE ::  &
    &  wfacpbl1(:,:),wfacpbl2(:,:) !< Corresponding interpolation coefficients

  INTEGER ::   &
    &  init_mode, ierror   !< identifier how to initialise the tracer
  CHARACTER(:), ALLOCATABLE          ::  &
    &  init_name           !< name of the variable in the external dataset
  TYPE(t_var_metadata), POINTER      :: & 
    &  info                !< returns reference to tracer metadata of current element
  TYPE(t_var_metadata_dynamic), POINTER :: & 
    &  info_dyn            !< returns reference to tracer metadata of current element 
  INTEGER, POINTER ::  &
    &  jsp                 !< returns index of element 
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo            !< ART atmo fields

  art_atmo => p_art_data(jg)%atmo


  ALLOCATE(wfac_lin(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ALLOCATE(idx0_lin(art_atmo%nproma,art_atmo%nlev,art_atmo%nblks))
  ALLOCATE(bot_idx_lin(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(kpbl1(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(kpbl2(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(wfacpbl1(art_atmo%nproma,art_atmo%nblks))
  ALLOCATE(wfacpbl2(art_atmo%nproma,art_atmo%nblks))


  no_levels = chem_init%chem_init_in%nlev_chem_init

  ! l_extrapol = .FALSE. means that of the uppermost level in the external data
  ! is used above top of external data
  CALL prepare_lin_intp(chem_init%chem_init_in%z3d, art_atmo%z_mc,     &
                      & art_atmo%nblks, art_atmo%npromz, no_levels,    &
                      & art_atmo%nlev,                                 &
                      & wfac_lin, idx0_lin, bot_idx_lin,               &
                      & lextrap = .FALSE.)

  CALL prepare_extrap(chem_init%chem_init_in%z3d,                     &
                       & art_atmo%nblks, art_atmo%npromz, no_levels,  &
                       & kpbl1, wfacpbl1, kpbl2, wfacpbl2 )


  ! if not for initialisation
  IF (.NOT. PRESENT(p_prog_list)) THEN

    ! go through all variables and interpolate vertically
    ! ATTENTION: This always includes an extrapolation downwards. With
    !            "l_extrapol" you can only choose how it is done: .false. means
    !            that the values of the bottom external level are copied to
    !            the ICON levels that are below this level
    DO i=1,num_vars

      CALL lin_intp(chem_init%chem_init_chem%spec(:,:,:,i),                &
                &   chem_init%chem_init_chem%spec_interp(:,:,:,i),         &
                &   art_atmo%nblks, art_atmo%npromz, no_levels,            &
                &   art_atmo%nlev,    &
                &   wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,      &
                &   wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.False., &
                &   l_pd_limit=.False.)
      ! DISCLAIMER:
      ! This line is PRELIMINARY done to prevent RESTART-issues.
      ! It kind of resets the interpolation result, which is currently operating
      ! differently when in a RESTART-environment due to a changed reference
      ! field. THIS IS NOT THE INTENDED SOLUTION FOR THE PROBLEM.
      ! ---------------------------------------------------------------------
      chem_init%chem_init_chem%spec_interp(:,:,:,i) = chem_init%chem_init_chem%spec(:,:,:,i) !temp RESTART
    END DO
  ELSE
    DO iv = 1,p_prog_list%p%nvars

      info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn
      info=>p_prog_list%p%vl(iv)%p%info

      IF (info_dyn%tracer%lis_tracer) THEN

        jsp=>info%ncontained
        SELECT TYPE(meta=>info_dyn%tracer)
          CLASS IS (t_chem_meta)
            CALL meta%opt_meta%get('init_mode',init_mode,ierror)
            IF (ierror /= SUCCESS) CALL finish(TRIM(routine)//':art_vinterp_chem_init',   &
                                 &         'Metadata init_mode not available for tracer ' &
                                 &        //TRIM(meta%name)//'.')

            IF (init_mode == 1) THEN

              CALL key_value_storage_as_string(meta%opt_meta,'init_name',init_name,ierror)
              IF (ierror/= SUCCESS) CALL finish(TRIM(routine)//':art_vinterp_chem_init',    &
                                  &         'Metadata init_name not available for tracer '  &
                                  &        //TRIM(meta%name)//'.')

              CALL message('mo_art_vinterp ', 'interpolate ' // TRIM(init_name) // ' ...' )

              CALL lin_intp(chem_init%chem_init_chem%spec(:,:,:,jsp),                   &
                      &     chem_init%chem_init_chem%spec_interp(:,:,:,jsp),            &
                      &     art_atmo%nblks, art_atmo%npromz, no_levels, art_atmo%nlev,  &
                      &     wfac_lin, idx0_lin, bot_idx_lin, wfacpbl1, kpbl1,           &
                      &     wfacpbl2, kpbl2, l_loglin=.FALSE., l_extrapol=.False.,      &
                      &     l_pd_limit=.False.)

           ENDIF

        END SELECT
      ENDIF
    ENDDO 
  END IF

END SUBROUTINE art_vinterp_chem_init

END MODULE mo_art_vinterp_chem_init


