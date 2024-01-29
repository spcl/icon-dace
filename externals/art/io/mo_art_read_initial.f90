!
! mo_art_read_initial
! This module sets initial values for aerosol and chemical tracers.
! Note, that the check for namelist settings (i.e. lart_aerosol, lart_chem)
! is performed before the according subroutines are called
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

MODULE mo_art_read_initial
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_nonhydrostatic_config,         ONLY: kstart_tracer
  USE mo_exception,                     ONLY: message, message_text
  USE mo_util_cdi,                      ONLY: t_inputparameters
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: t_var_metadata,t_var_metadata_dynamic
  USE mo_tracer_metadata_types,         ONLY: t_aero_meta, t_chem_meta
  USE mo_impl_constants,                ONLY: VARNAME_LEN => vname_len
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_create_filenames,          ONLY: art_create_filenames
  USE mo_art_config,                    ONLY: art_config
  USE mo_art_read,                      ONLY: art_read, art_open_cdi, art_close_cdi
  USE mo_art_io_constants,              ONLY: IINIT_AERO, IINIT_CHEM

  USE mo_bc_sulfate_aerosol,            ONLY: read_echam_bc_sulfate_aerosol, &
       &                                      destruct_echam_bc_sulfate_aerosol, &
       &                                      ext_sulfate_aerosol

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_read_initial'

  PUBLIC :: art_init_aero
  PUBLIC :: art_init_chem

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init_aero(jg, p_prog_list, tracer, use_echam_climatology)
!<
! SUBROUTINE art_init_aero
! This subroutine performs the aerosol tracer initialization from files.
! Part of Module: mo_art_read_initial
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-10
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  INTEGER, INTENT(in)                   :: &
    &  jg                                    !< patch on which computation is performed
  TYPE(t_var_list_ptr), INTENT(IN)          :: &
    &  p_prog_list                           !< current prognostic state list
  REAL(wp), INTENT(inout)               :: &
    &  tracer(:,:,:,:)                       !< Tracer mixing ratios [kg kg-1]
  LOGICAL, INTENT(in), OPTIONAL         :: &
    &  use_echam_climatology
! Local variables
  TYPE(t_inputParameters)               :: &
    &  cdi_param                             !< Parameters for read_cdi call
  TYPE(t_var_metadata), POINTER         :: &
    &  info                                  !< returns reference to current element metadata
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn                              !< returns reference to current element tracer metadata
  INTEGER                :: &
    &  nlen                   !< nproma or npromz_c
  INTEGER                :: &
    &  cdi_artdataset_id,   & !< CDI stream ID (for each domain)
    &  cdi_filetype           !< One of CDI's FILETYPE_XXX constants
  INTEGER, POINTER       :: &
    &  jsp                    !< returns index of element
  CHARACTER(LEN=200)     :: &
    &  aero_dataset           !< Path+filename of dataset for initialization of aerosol species
  INTEGER                :: &
    &  jb, jk, jc, iv         !< loop indices
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo
  CHARACTER(LEN=VARNAME_LEN) :: &
    &  varname

  !TODO: add namelist support later
  LOGICAL :: read_aerosol_background_from_echam

  IF (PRESENT(use_echam_climatology)) THEN
    read_aerosol_background_from_echam = use_echam_climatology
  ELSE
    read_aerosol_background_from_echam = .FALSE.
  ENDIF

  IF (read_aerosol_background_from_echam) THEN

    CALL read_echam_bc_sulfate_aerosol

    DO iv = 1, p_prog_list%p%nvars

      info=>p_prog_list%p%vl(iv)%p%info
      info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn

      IF (info_dyn%tracer%lis_tracer) THEN

        jsp => info%ncontained

        SELECT TYPE(meta=>info_dyn%tracer)

        CLASS IS (t_aero_meta)
          ! With 1:(LEN(TRIM(meta%name))-4) cut .TL1 of meta%name
          varname = meta%name(1:(LEN(TRIM(meta%name))-4))
          WRITE (message_text,*) 'ART: Reading aes climatology as initial data for '//TRIM(varname)
          CALL message(TRIM(routine)//':art_init_aero', message_text)

          SELECT CASE (TRIM(varname))

          CASE('nmb_insol_acc')
            tracer(:,:,:,jsp) = 0.0_wp
          CASE('so4_insol_acc')
            tracer(:,:,:,jsp) = 0.0_wp
          CASE('nmb_insol_coa')
            tracer(:,:,:,jsp) = 0.0_wp
          CASE('so4_insol_coa')
            tracer(:,:,:,jsp) = 0.0_wp
          CASE('nmb_sol_acc')
            tracer(:,:,:,jsp) = ext_sulfate_aerosol(jg)%num_as(:,:,:)
          CASE('so4_sol_acc')
            tracer(:,:,:,jsp) = ext_sulfate_aerosol(jg)%so4_as(:,:,:)
          CASE('nmb_sol_ait')
            tracer(:,:,:,jsp) = ext_sulfate_aerosol(jg)%num_ks(:,:,:)
          CASE('so4_sol_ait')
            tracer(:,:,:,jsp) = ext_sulfate_aerosol(jg)%so4_ks(:,:,:)
          CASE('nmb_mixed_coa')
            tracer(:,:,:,jsp) = ext_sulfate_aerosol(jg)%num_cs(:,:,:)
          CASE('so4_mixed_coa')
            tracer(:,:,:,jsp) = ext_sulfate_aerosol(jg)%so4_cs(:,:,:)
          CASE DEFAULT
            tracer(:,:,:,jsp) = 0.0_wp
          END SELECT

          ! tracer mapping can only be done be predefined association
          ! of modes, so we can hard code and finish if target mode is
          ! not available

          ! ART tracers are supposed to be zero in the region where tracer-related processes
          ! are turned off (above htop_art_proc, or rather above htop specified in tracer.xml)
          ! --compare copy_initicon2prog_atm for kstart_moist--
          ! really???
          ! DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          !   IF (jb /= art_atmo%nblks) THEN
          !     nlen = art_atmo%nproma
          !   ELSE
          !     nlen = art_atmo%npromz
          !   ENDIF
          !   DO jk = 1, kstart_tracer(jg,jsp)-1
          !     DO jc = 1, nlen
          !       tracer(jc,jk,jb,jsp) = 0.0_wp
          !     ENDDO
          !   ENDDO
          ! ENDDO

        CLASS DEFAULT
          ! Nothing to do here
        END SELECT
      ENDIF
    ENDDO


    CALL destruct_echam_bc_sulfate_aerosol

  ELSE

    art_atmo => p_art_data(jg)%atmo

    CALL art_create_filenames(jg,TRIM(art_config(jg)%cart_input_folder), &
         &                       IINIT_AERO,aero_dataset,.FALSE.)

    CALL art_open_cdi(jg,aero_dataset,cdi_artdataset_id,cdi_param,cdi_filetype)

    ! Loop through tracers and initialize
    DO iv = 1, p_prog_list%p%nvars

      info=>p_prog_list%p%vl(iv)%p%info
      info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn

      IF (info_dyn%tracer%lis_tracer) THEN

        jsp=>info%ncontained

        SELECT TYPE(meta=>info_dyn%tracer)

          CLASS IS (t_aero_meta)
            ! With 1:(LEN(TRIM(meta%name))-4) cut .TL1 of meta%name
            varname = meta%name(1:(LEN(TRIM(meta%name))-4))
            WRITE (message_text,*) 'ART: Reading initial data for '//TRIM(varname)
            CALL message (TRIM(routine)//':art_init_aero', message_text)
            CALL art_read(jg, cdi_artdataset_id, cdi_param, varname,       &
                 &           tracer(:,:,:,jsp), art_atmo%nlev, cdi_filetype)


            ! ART tracers are supposed to be zero in the region where tracer-related processes
            ! are turned off (above htop_art_proc, or rather above htop specified in tracer.xml)
            ! --compare copy_initicon2prog_atm for kstart_moist--
            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
              IF (jb /= art_atmo%nblks) THEN
                nlen = art_atmo%nproma
              ELSE
                nlen = art_atmo%npromz
              ENDIF
              DO jk = 1, kstart_tracer(jg,jsp)-1
                DO jc = 1, nlen
                  tracer(jc,jk,jb,jsp) = 0.0_wp
                ENDDO
              ENDDO
            ENDDO

          CLASS DEFAULT
          ! Nothing to do here
        END SELECT
      ENDIF
    ENDDO

    ! Cleanup
    CALL art_close_cdi(cdi_artdataset_id, cdi_param)

  ENDIF

END SUBROUTINE art_init_aero
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_init_chem(jg, p_prog_list, tracer)
!<
! SUBROUTINE art_init_chem
! This subroutine performs the chemical tracer initialization from files.
! Part of Module: mo_art_read_initial
! Author: Daniel Rieger, KIT
! Initial Release: 2016-08-10
! Modifications:
! YYYY-MM-DD: <name>,<institution>
! - <description>
!>
  INTEGER, INTENT(in)                   :: &
    &  jg                                    !< patch on which computation is performed
  TYPE(t_var_list_ptr), INTENT(IN)          :: &
    &  p_prog_list                           !< current prognostic state list
  REAL(wp), INTENT(inout)               :: &
    &  tracer(:,:,:,:)                       !< Tracer mixing ratios [kg kg-1]
! Local variables
  TYPE(t_inputParameters)               :: &
    &  cdi_param                             !< Parameters for read_cdi call
  TYPE(t_var_metadata), POINTER         :: &
    &  info                                  !< returns reference to current element metadata
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn                              !< returns reference to current element tracer metadata
  INTEGER                :: &
    &  cdi_artdataset_id,   & !< CDI stream ID (for each domain)
    &  cdi_filetype           !< One of CDI's FILETYPE_XXX constants
  INTEGER, POINTER       :: &
    &  jsp                    !< returns index of element
  CHARACTER(LEN=200)     :: &
    &  chem_dataset           !< Path+filename of dataset for initialization of chemical species
  INTEGER                :: &
    &  jb, jk, jc, iv         !< loop indices
  INTEGER                :: &
    &  nlen                   !< nproma or npromz_c
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo

  art_atmo => p_art_data(jg)%atmo

  CALL art_create_filenames(jg,TRIM(art_config(jg)%cart_input_folder), &
    &                       IINIT_CHEM, chem_dataset,.FALSE.)

  CALL art_open_cdi(jg,chem_dataset,cdi_artdataset_id,cdi_param,cdi_filetype)

! Loop through tracers and initialize

  DO iv = 1, p_prog_list%p%nvars

    info=>p_prog_list%p%vl(iv)%p%info
    info_dyn=>p_prog_list%p%vl(iv)%p%info_dyn

    IF (info_dyn%tracer%lis_tracer) THEN

      jsp=>info%ncontained

      SELECT TYPE(meta=>info_dyn%tracer)

        CLASS IS (t_chem_meta)
          WRITE (message_text,*) 'ART: Reading initial data for '//TRIM(meta%name)
          CALL message (TRIM(routine)//':art_init_chem', message_text)
          CALL art_read(jg, cdi_artdataset_id, cdi_param, meta%name,  &
             &           tracer(:,:,:,jsp), art_atmo%nlev, cdi_filetype)

          ! ART tracers are supposed to be zero in the region where tracer-related processes
          ! are turned off (above htop_art_proc, or rather above htop specified in tracer.xml)
          ! --compare copy_initicon2prog_atm for kstart_moist--
          DO jb = art_atmo%i_startblk, art_atmo%i_endblk
            IF (jb /= art_atmo%nblks) THEN
              nlen = art_atmo%nproma
            ELSE
              nlen = art_atmo%npromz
            ENDIF

            DO jk = 1, kstart_tracer(jg,jsp)-1
              DO jc = 1, nlen
                tracer(jc,jk,jb,jsp) = 0.0_wp
              ENDDO
            ENDDO
          ENDDO

        CLASS DEFAULT
          ! Nothing to do here
      END SELECT
    ENDIF
  ENDDO

! Cleanup
  CALL art_close_cdi(cdi_artdataset_id, cdi_param)

END SUBROUTINE art_init_chem
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_read_initial
