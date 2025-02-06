!
! mo_art_init_full_chemistry
! This module provides initialization routines for full gas phase routine
! Using lart_mecca
!
! Species are either initialized using fixed values or 3d Init data
!  (e.g. MOZART, EMAC, ...)
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

MODULE mo_art_init_full_chemistry
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: message,finish, message_text
  USE mo_var_list,                      ONLY: t_var_list_ptr
  USE mo_var,                           ONLY: t_var
  USE mo_var_metadata_types,            ONLY: t_var_metadata_dynamic,    &
                                          &   t_var_metadata
  USE mo_tracer_metadata_types,         ONLY: t_tracer_meta
  USE mo_run_config,                    ONLY: ntracer
  USE mo_art_config,                    ONLY: art_config
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_chem_data,                 ONLY: t_art_chem
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_read_xml,                  ONLY: t_xml_file, art_open_xml_file,  &
                                          &   art_close_xml_file,             &
                                          &   art_get_text_attribute_xml
  USE mo_art_chem_types,                ONLY: t_chem_meta_mecca
  USE mo_art_chem_init_types,           ONLY: t_chem_init_state
  USE mo_mecicon_init,                  ONLY: mecicon_init, initialize_kpp_variables
  USE messy_mecca_kpp_parameters
  USE messy_mecca_kpp_global,           ONLY: c

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_init_full_chemistry'

  PUBLIC :: art_init_full_chemistry, art_init_mapping_kpp_indices

CONTAINS


SUBROUTINE art_init_full_chemistry(jg, p_prog_list, p_tracer_now)
!<
! SUBROUTINE art_init_full_chemistry
! This subroutine sets initial values for MECCA tracers
! Part of Module: mo_art_init_full_chemistry
! Author: Jennifer Schroeter, KIT
! Initial Release: 2016-08-01
! Modifications:
! 2016-11-08: Michael Weimer, KIT
! - removed emission initialisation
!>
  IMPLICIT NONE

  INTEGER, INTENT(in)  ::  &
    &  jg                              !< patch id on which computation is performed
  TYPE(t_var_list_ptr), INTENT(IN)  ::  &
    &  p_prog_list                     !< current prognostic state list
  REAL(wp), INTENT(INOUT)  ::  &
    &  p_tracer_now(:,:,:,:)           !< tracer volume mxing ratios (mol/mol)

  ! local variables
  TYPE(t_var_metadata), POINTER  ::  &
    &  info                            !< returns reference to tracer
  TYPE(t_chem_init_state),POINTER       :: &
    &  chem_init                       !< reference to chem_init variables
  TYPE(t_art_atmo), POINTER             :: &
    &  art_atmo                        !< ART atmo fields
  TYPE(t_art_chem), POINTER             :: &
    &  art_chem                        !< ART chem fields
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn                        !< dynamic tracer metadata
  INTEGER ::                      &
    &    jc, jk, jb, iv,          &    !< loop indizes
    &    i_startidx, i_endidx
  INTEGER ::    &
    &  kpp_ind                         !< MECCA tracer index
  INTEGER, POINTER  ::  &
    &  jsp                             !< returns index of element
  REAL(wp)  ::  &
    &  factor                          !< factor when initialising a tracer with another one
  REAL(wp), ALLOCATABLE ::   &
    &  vmr_kpp(:)                      !< volume mixing ratio for initialisation
  INTEGER, POINTER :: &
    &  mapping_indices_kpp(:)


  ! ######################################
  ! Init
  ! ######################################

  chem_init           =>  p_art_data(jg)%chem%init
  art_atmo            =>  p_art_data(jg)%atmo
  art_chem            =>  p_art_data(jg)%chem
  mapping_indices_kpp => art_chem%mecicon%utils%mapping_indices_kpp

  CALL message(routine, 'init gasphase')

  !#######################################
  ! set initial fields for those tracers with init_mode = 1 by external dataset
  ! set fixed values from MECCA for init_mode = 0
  !#######################################


  IF (.NOT. ALLOCATED(vmr_kpp)) ALLOCATE(vmr_kpp(ntracer))
           ! allocation of c

  c(:) = 0.0e-0_wp

  DO iv = 1, p_prog_list%p%nvars

    info_dyn => p_prog_list%p%vl(iv)%p%info_dyn
    info     => p_prog_list%p%vl(iv)%p%info
    jsp      => info%ncontained

    SELECT TYPE(meta => info_dyn%tracer)
      CLASS IS (t_chem_meta_mecca)

        IF ((meta%init_mode == 0) .OR. (art_config(jg)%iart_init_gas == 0)) THEN

          kpp_ind = meta%number

          CALL message(routine, 'initializing ' // TRIM(meta%name)  &
                   & //' by MECCA value ...')

          DO jb = art_atmo%i_startblk, art_atmo%i_endblk
            CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
            DO jk=1,art_atmo%nlev
              DO jc=i_startidx,i_endidx

              CALL mecicon_init( art_atmo%pres(jc,jk,jb),     &
                                   art_atmo%temp(jc,jk,jb),     &
                                   vmr_kpp(kpp_ind), kpp_ind)

                p_tracer_now(jc,jk,jb,jsp) = vmr_kpp(kpp_ind)

              ENDDO
            ENDDO
          ENDDO

        ENDIF
     END SELECT

     ENDDO

  ! ######################################################
  ! special treatment for some species 
  ! ######################################################

  DO iv = 1, p_prog_list%p%nvars

    info_dyn  =>  p_prog_list%p%vl(iv)%p%info_dyn
    info      =>  p_prog_list%p%vl(iv)%p%info
    jsp       =>  info%ncontained
    SELECT TYPE(meta => info_dyn%tracer)
      CLASS IS (t_chem_meta_mecca)

        SELECT CASE (ADJUSTL(TRIM(meta%name)))

          CASE("H2O_full")

          !  CALL message(routine, 'Setting H2O to qv values')
          !  DO jb = art_atmo%i_startblk, art_atmo%i_endblk
          !    CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
          !    DO jk=1,art_atmo%nlev
          !      DO jc=i_startidx,i_endidx

          !         p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb, iqv) * (amd / amw)

          !       ENDDO
          !     ENDDO
          !    ENDDO

          CASE("O3_full")
            ! In general, iTRO3 is a parametrised tracer, but it is set here
            ! for MECCA chemistry for radiation feedback
            art_chem%indices%iTRO3 = jsp
!             CALL message(routine, 'init o3 with gems')
!             
!
!             !  CALL calc_o3_gems_linoz(p_patch,current_date,p_diag,p_tracer_now(:,:,:,jsp))
!             DO jb = art_atmo%i_startblk, art_atmo%i_endblk
!               CALL art_get_indices_c(jg, jb,  i_startidx, i_endidx)
!
!               DO jk=1,art_atmo%nlev
!                 DO jc=i_startidx,i_endidx
!
!                   p_tracer_now(jc,jk,jb,jsp)  &
!                      &  = chem_init%chem_init_chem%spec_interp(jc,jk,jb,jsp)/( ppmv2gg * 1.e6)
!
!                 ENDDO
!               ENDDO
!             ENDDO 
!
!         CASE("CFCl3_full")
!
!           CALL message(routine, 'Scaling CFCL3 by CH4 profile')
!
!            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
!              CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
!                DO jc=i_startidx,i_endidx
!
!                  factor = p_tracer_now(jc,art_atmo%nlev-1,jb,jsp)  &
!                   &       / p_tracer_now(jc,art_atmo%nlev-1,jb,mapping_indices_kpp(ind_CH4))
!   
!                  DO jk=1,art_atmo%nlev
!
!                     p_tracer_now(jc,jk,jb,jsp) = & 
!                         factor*p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_CH4))
!
!                  ENDDO
!                ENDDO
!              ENDDO


          CASE("CF3Br_full")

            CALL message(routine, 'Scaling CF3Br by CH4 profile')
            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
              CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

              DO jc = i_startidx,i_endidx
                factor = p_tracer_now(jc,art_atmo%nlev-1,jb,jsp)   &
                  &   /  p_tracer_now(jc,art_atmo%nlev-1,jb,mapping_indices_kpp(ind_CH4))
                
                DO jk = 1,art_atmo%nlev
                  p_tracer_now(jc,jk,jb,jsp) = factor  & 
                       &    * p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_CH4))
                ENDDO
              ENDDO
            ENDDO


          CASE("CF2ClBr_full")

            CALL message(routine, 'Scaling CF2ClBr by CH4 profile')
            DO jb = art_atmo%i_startblk, art_atmo%i_endblk
              CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
              DO jc = i_startidx,i_endidx
                factor = p_tracer_now(jc,art_atmo%nlev-1,jb,jsp)  &
                    &  / p_tracer_now(jc,art_atmo%nlev-1,jb,mapping_indices_kpp(ind_CH4))

                DO jk = 1,art_atmo%nlev

                  p_tracer_now(jc,jk,jb,jsp) = factor & 
                       &    * p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_CH4))

                ENDDO
              ENDDO
            ENDDO

          CASE("OCS_full")
            IF (ind_N2O /= 0) THEN
              CALL message(routine,'Scaling OCS by N2O profile')
              DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

                DO jc = i_startidx,i_endidx
                  p_tracer_now(jc,art_atmo%nlev,jb,jsp) = 4.e-10_wp
  
                  factor = p_tracer_now(jc,art_atmo%nlev,jb,jsp)  &
                    &    / p_tracer_now(jc,art_atmo%nlev,jb,mapping_indices_kpp(ind_N2O))
  
                  DO jk = 1,art_atmo%nlev

                     p_tracer_now(jc,jk,jb,jsp) = factor &
                         &    * p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_N2O))

                  ENDDO
                ENDDO
              ENDDO
             
            END IF

           CASE("CF2Cl2_full")
             CALL message(routine, 'Scaling CF2Cl2 by CH4 profile')
               DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                 CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

                 DO jc = i_startidx,i_endidx
                  factor = p_tracer_now(jc,art_atmo%nlev-1,jb,jsp) &
                     &   / p_tracer_now(jc,art_atmo%nlev-1,jb,mapping_indices_kpp(ind_CH4))

                  DO jk = 1,art_atmo%nlev

                     p_tracer_now(jc,jk,jb,jsp) = factor &
                          &    * p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_CH4))

                   ENDDO
                 ENDDO
               ENDDO
!           CASE("N2O_full")
!             CALL message(routine, 'FOUND N2O in tracers')
!             DO jb = art_atmo%i_startblk, art_atmo%i_endblk
!               CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)
!                 DO jc=i_startidx,i_endidx
!
!                   factor = p_tracer_now(jc,art_atmo%nlev-1,jb,jsp) &
!                     &       / p_tracer_now(jc,art_atmo%nlev-1,jb,mapping_indices_kpp(ind_CH4))
!                 
!                   DO jk=1,art_atmo%nlev
!
!                      p_tracer_now(jc,jk,jb,jsp) = & 
!                         factor*p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_CH4))
!
!                   ENDDO
!                 ENDDO
!               ENDDO

               END SELECT
            ! additionally initialise some passive tracers if they are there
      ! additionally initialise some passive tracers if they are there
      CLASS IS (t_tracer_meta)
        SELECT CASE (ADJUSTL(TRIM(meta%name)))
          CASE('TRO3_passive')
            IF (ind_O3 /= 0) THEN
              CALL message(routine,    &
                    &      'initialise passive ozone by gas phase chemistry ozone')

              DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

                DO jk = 1,art_atmo%nlev
                  DO jc = i_startidx,i_endidx
                    p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,  &
                                               &  mapping_indices_kpp(ind_O3))
                  END DO
                END DO
              END DO
            END IF
          CASE('TRNOy_passive')
            IF (ANY( (/ ind_HNO3, ind_NO, ind_NO2, ind_ClNO2, ind_N,  &
                 &      ind_HNO4, ind_N2O5, ind_BrNO3, ind_ClNO3, ind_NO3 /) /= 0 )) THEN
              CALL message(routine,    &
                    &      'initialise passive NOy by other substances')

              DO jb = art_atmo%i_startblk, art_atmo%i_endblk
                CALL art_get_indices_c(jg, jb, i_startidx, i_endidx)

                DO jk = 1,art_atmo%nlev
                  DO jc = i_startidx,i_endidx
                    p_tracer_now(jc,jk,jb,jsp) = 0._wp
                    IF (ind_HNO3 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_HNO3))
                    END IF
            
                    IF (ind_NO /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_NO))
                    END IF
            
                    IF (ind_NO2 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_NO2))
                    END IF

                    IF (ind_ClNO2 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_ClNO2))
                    END IF
            
                    IF (ind_N /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_N))
                    END IF
            
                    IF (ind_HNO4 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_HNO4))
                    END IF
            
                    IF (ind_N2O5 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + 2._wp * p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_N2O5))
                    END IF
            
                    IF (ind_BrNO3 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_BrNO3))
                    END IF

                    IF (ind_ClNO3 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_ClNO3))
                    END IF
            
                    IF (ind_NO3 /= 0) THEN
                      p_tracer_now(jc,jk,jb,jsp) = p_tracer_now(jc,jk,jb,jsp) &
                                &   + p_tracer_now(jc,jk,jb,mapping_indices_kpp(ind_NO3))
                    END IF
            
                  END DO
                END DO
              END DO
            END IF

        END SELECT
              
    END SELECT

  ENDDO

  NULLIFY(chem_init)

END SUBROUTINE art_init_full_chemistry
!
!###################################
!

SUBROUTINE art_init_mapping_kpp_indices(jg, p_prog_list, cart_mecca_xml)
!<
! SUBROUTINE art_init_mapping_kpp_indices
! This subroutine sets the mapping array between ICON and MECCA tracers
! Part of Module: mo_art_init_full_chemistry
! Author: Michael Weimer, KIT
! Initial Release: 2018-03-29
! Modifications:
!>
  IMPLICIT NONE
  INTEGER :: &
    &  jg                !< patch id
  TYPE(t_var_list_ptr), INTENT(in) :: &
    &  p_prog_list       !< list of prognostic variables
  CHARACTER(LEN=*), INTENT(in) :: &
    &  cart_mecca_xml    !< XML file
  !local variables
  TYPE(t_var_metadata), POINTER         :: &
    &   info            !< returns reference to tracer
  TYPE(t_var_metadata_dynamic), POINTER :: &
    &  info_dyn         !< returns reference to dynamic meta data
  INTEGER, POINTER :: &
    &  jsp              !< tracer index
  INTEGER ::   &
    &  kpp_ind, iTR,  & !< kpp tracer index and ICON-ART tracer index
    &  iv               !< loop index
  TYPE(t_xml_file) :: &
    &  tixi_file        !< XML file handle
  CHARACTER(LEN=50) :: &
    &  checksum_xml     !< checksum of the XML file
  INTEGER ::   &
    &  nvars_chem       !< number of MECCA tracers in the XML file

  ! ######################################
  ! consistency check mechanisms
  ! ######################################

  CALL art_open_xml_file(TRIM(cart_mecca_xml),tixi_file)

  CALL art_get_text_attribute_xml(tixi_file,'/tracers','checksum',checksum_xml)

  CALL art_close_xml_file(tixi_file)

  IF (TRIM(checksum_xml) /= TRIM(KPP_XML_CHECKSUM)) THEN
    CALL finish(TRIM(routine)//':art_init_mapping_kpp_indices',                                  &
           &    'Your KPP chemistry mechanism and your chemistry XML file do not fit together.'  &
           &  //' Please consider updating them: xml:'//TRIM(checksum_xml)//' vs. kpp:'          &
           &  //TRIM(TRIM(KPP_XML_CHECKSUM)))
  END IF

  ! ######################################
  ! create the mapping array
  ! ######################################

  IF (.NOT. ALLOCATED(p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp)) THEN
    ALLOCATE( p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp(ntracer))

    DO iTR = 1,ntracer
      p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp(iTR) = -1
    END DO

    DO iv = 1, p_prog_list%p%nvars

      info_dyn  =>  p_prog_list%p%vl(iv)%p%info_dyn
      info      =>  p_prog_list%p%vl(iv)%p%info
      jsp       =>  info%ncontained

      SELECT TYPE(meta => info_dyn%tracer)
        CLASS IS (t_chem_meta_mecca)

          kpp_ind = meta%number

          IF ((kpp_ind <= 0) .OR. (kpp_ind > ntracer)   &
            &  .OR. (p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp(kpp_ind) /= -1)) THEN
            CALL finish(routine//':art_init_mapping_kpp_indices', &
                     &  'Meta data number for '//TRIM(meta%name)//' inconsistent in tracer XML. ' &
                     &//'Please consider creating the XML file again. '  &
                     &//'Reason is either that the number is not unique or not in [1,ntracer].')
          END IF

          p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp(kpp_ind) = jsp
      END SELECT


    ENDDO
  END IF ! mapping_indices_kpp allocated

  CALL initialize_kpp_variables()

  ! consistency check: equals the number of chemical tracers that of the
  ! mechanism?
  nvars_chem = COUNT(p_art_data(jg)%chem%mecicon%utils%mapping_indices_kpp > 0)

  IF (nvars_chem /= NSPEC) THEN
    WRITE(message_text,*) 'Number of MECCA tracers in ', TRIM(cart_mecca_xml),   &
             &            ' does not equal that of MECCA mechanism: xml:', nvars_chem,  &
             &            ' vs. KPP:', NSPEC
    CALL finish(routine//':art_init_mapping_kpp_indices', message_text)
  END IF
END SUBROUTINE art_init_mapping_kpp_indices

END MODULE mo_art_init_full_chemistry

