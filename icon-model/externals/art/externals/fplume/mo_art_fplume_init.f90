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

MODULE mo_art_fplume_init
  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: message,finish,message_text
  USE mtime,                            ONLY: timedelta, datetime
  USE mo_art_config,                    ONLY: IART_PATH_LEN
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH
  USE mo_model_domain,                  ONLY: t_patch

  USE mo_art_fplume_types,              ONLY: t_art_all_volc_fplume
  USE mo_art_fplume_read_inp,           ONLY: reainp,reagrn
  USE mo_art_fplume_utilities,          ONLY: get_input_npar,get_input_rea

  IMPLICIT NONE

  PRIVATE

  PUBLIC:: art_fplume_init

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_init'

  CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_fplume_init(p_patch, fplume_init_all, volc_numb, inp_path, tc_exp_refdate,        &
                         & tc_dt_model)
  TYPE(t_patch), INTENT(in)      :: &
    &  p_patch
  TYPE(t_art_all_volc_fplume), POINTER, INTENT(in)   :: &
    & fplume_init_all
  INTEGER, INTENT(in)                                :: &
    & volc_numb
  CHARACTER(LEN=IART_PATH_LEN), INTENT(in)                              :: &
    & inp_path
  TYPE(timedelta), POINTER, INTENT(in)               :: &
    &  tc_dt_model                 !< Model timestep
  TYPE(datetime), POINTER, INTENT(in)                :: &
    &  tc_exp_refdate              !< Experiment reference date
  CHARACTER(LEN=IART_PATH_LEN)                       :: &
    &  input_file1, input_file2
  CHARACTER(LEN=MAX_CHAR_LENGTH)  :: mess
  INTEGER :: iphase,ivolcs,istat
  CHARACTER(1)  :: cvolc
  LOGICAL :: l_exist1,l_exist2
  REAL(wp)                       :: &
    &  rvoid(2),                    &
    &  lon_all,                     &
    &  lat_all

    DO ivolcs = 1, volc_numb
      ! Construct filenames
      IF (volc_numb==1) THEN
        input_file1=TRIM(inp_path)//'.inp'
        input_file2=TRIM(inp_path)//'.tgsd'
        INQUIRE(file = input_file1, EXIST = l_exist1)
        INQUIRE(file = input_file2, EXIST = l_exist2)
      ELSE
        WRITE(cvolc,'(I1)') ivolcs
        input_file1=TRIM(inp_path)//'0'//TRIM(cvolc)//'.inp'
        input_file2=TRIM(inp_path)//'0'//TRIM(cvolc)//'.tgsd'
        INQUIRE(file = input_file1, EXIST = l_exist1)
        INQUIRE(file = input_file2, EXIST = l_exist2)
      ENDIF

      IF (l_exist1 .AND. l_exist2) THEN
        WRITE(message_text,*) ': read ',input_file1, ' and ',input_file2
        CALL message(routine,message_text)
      ELSE
        CALL finish(routine,': FPlume input files missing')
      ENDIF

      ! read all locations
      CALL get_input_rea(input_file1,'SOURCE','LON_VENT',rvoid,1,istat,mess)
      IF (istat>0) THEN
        WRITE(message_text,*) mess
        CALL message('',message_text)
      ENDIF
      IF(istat<0) CALL finish(routine,mess)
      lon_all=rvoid(1)
      !
      CALL get_input_rea(input_file1,'SOURCE','LAT_VENT',rvoid,1,istat,mess)
      IF (istat>0) THEN
        WRITE(message_text,*) mess
        CALL message('',message_text)
      ENDIF
      IF (istat<0) CALL finish(routine,mess)
      lat_all=rvoid(1)

      ! initialize volcano
      IF (ivolcs==1) ALLOCATE(fplume_init_all%p(volc_numb))
      CALL fplume_init_all%p(ivolcs)%init(p_patch,lon_all,        &
              &                           lat_all)
      ! read phase independent variables
      CALL get_input_npar(input_file1,'SOURCE','EXIT_VELOCITY_(MS)',                             &
                        & fplume_init_all%p(ivolcs)%nphases,istat,message_text)
      IF (istat>0) THEN
        CALL message(routine,message_text)
      ENDIF
      IF (istat<0) CALL finish(routine,message_text)

      CALL reagrn(input_file1,input_file2,fplume_init_all%p(ivolcs))

      ! initialize eruption phase
      DO iphase = 1, fplume_init_all%p(ivolcs)%nphases
        CALL reainp(input_file1,iphase,tc_exp_refdate,tc_dt_model,fplume_init_all%p(ivolcs))
      ENDDO

    ENDDO
END SUBROUTINE
!!
!!-------------------------------------------------------------------------
!!
END MODULE
