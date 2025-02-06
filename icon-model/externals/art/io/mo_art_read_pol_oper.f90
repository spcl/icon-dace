!
! mo_art_read_pol_oper
! This module provides is a wrapper to read pollen data which are used in the parametrization
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

MODULE mo_art_read_pol_oper
! ICON
  USE mo_exception,                     ONLY: message
  USE mo_util_cdi,                      ONLY: t_inputParameters
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_config,                    ONLY: art_config !! TEMPORARY!!!
  USE mo_art_read,                      ONLY: art_read, art_open_cdi, art_close_cdi
  USE mo_art_external_types,            ONLY: t_art_pollen_table
  USE mo_art_create_filenames,          ONLY: art_create_filenames
  USE mo_art_io_constants,              ONLY: IEXT_POVAR, IART_FILENAMELEN

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_read_pol_oper'

  PUBLIC :: art_read_pol_oper

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_pol_oper(jg)
!<
! SUBROUTINE art_read_pol_oper
! Read distribution of pollen plants and parameter for altitude correction 
! Part of Module: mo_art_read_extdata
! Author: Jonas Straub, KIT
! Initial Release: 2016-07-21
! Modifications:
! 2016-12-21: Jonas Straub, KIT
! - Adaption to new structure
!>
  INTEGER, INTENT(in)                       :: &
    &  jg                                        !< domain index
                        
  ! Local variables:
  TYPE(t_inputParameters)                   :: &
    &  cdi_param                                 !< Parameters for read_cdi call
  TYPE(t_art_pollen_table), POINTER         :: &
    &  this_pollen_table                         !< pointer to pollen table
  INTEGER                                   :: &
    &  cdi_artdataset_id                         !< CDI stream ID (for each domain)
  INTEGER                                   :: &
    &  n, iv                                     !< counter
  CHARACTER(LEN=IART_FILENAMELEN)           :: &
    &  pollendataset                             !< Path and filename of pollen coverage dataset and altitude correction
  CHARACTER(LEN=30)                         :: &
    &  vpref,                                  & !< Name of coverage variable and altitude correction to be read
    &  vname                                  
  CHARACTER(LEN=10), ALLOCATABLE            :: &
    &  vsuff(:)                                !< Name of coverage variable and altitude correction to be read

  ! Loop through pollen types and read plant coverage
  DO n=1, p_art_data(jg)%ext%pollen_prop%npollen_types
    this_pollen_table => p_art_data(jg)%ext%pollen_prop%pollen_type(n)
    vname = this_pollen_table%sname
    vpref = ''
    ! Last pollentype was reached and all data reading is finished
    !_jf: on NEC vector host 'this_pollen_table%fr_cov' is associated?!
    !_jf: as workaround test if variable shortname is set.
    IF (.NOT. ASSOCIATED(this_pollen_table%fr_cov) .OR. TRIM(this_pollen_table%shortname) == '') THEN
      CALL message(TRIM(routine)//': art_read_pol_oper', 'Left reading loop')
      EXIT
    END IF 
  
    IF (TRIM(vname) /= '') THEN

      SELECT CASE(TRIM(vname))
        CASE('pollbetu')
          vpref      = 'BETU'
        CASE('pollalnu')
          vpref      = 'ALNU'
        CASE('pollcory')
          vpref      = 'CORY'
        CASE('pollpoac')
          vpref      = 'POAC'
        CASE('pollambr')
          vpref      = 'AMBR'
        CASE DEFAULT
          CALL message(TRIM(routine)//': art_read_pol_oper', &
            &          'pollentype is not yet implemented'//TRIM(vname))
          CONTINUE
      END SELECT

      ALLOCATE(vsuff(11))
      vsuff      = (/'rprec','reso ','ress ','sdes ','ctsum','saisn', 'saisa', 'tthrs', 'tthre', 'saisl','tune '/) 

      CALL message(TRIM(routine)//':art_read_pol_oper',      &
        &          'Reading extern data for pollentype: '//TRIM(vname))

      CALL art_create_filenames(jg, TRIM(art_config(jg)%cart_input_folder), IEXT_POVAR, pollendataset, &
        &                       .FALSE.)
      CALL art_open_cdi(jg,pollendataset,cdi_artdataset_id,cdi_param)

      DO iv=1, SIZE(vsuff)

        SELECT CASE (TRIM(vsuff(iv)))
          CASE ('rprec')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%r_precip)
          CASE ('reso')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%res_old)
          CASE ('ress')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%res_new_sum)
          CASE ('sdes')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%f_q_seas)
          CASE ('ctsum')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%ctsum)
          CASE ('saisn')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%saisn)
          CASE ('saisa')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%saisa)
          CASE ('tthrs')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%tthrs)
          CASE ('tthre')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%tthre)
          CASE ('saisl')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%saisl)
          CASE ('tune')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%tune)
          CASE DEFAULT
            CALL message(TRIM(routine)//': art_read_pol_oper', &
              &          'input field is not defined'//TRIM(vpref)//TRIM(vsuff(iv)))
        END SELECT

      END DO

      CALL art_close_cdi(cdi_artdataset_id, cdi_param)
      DEALLOCATE( vsuff )
 
    END IF 
  ENDDO

  NULLIFY (this_pollen_table)

END SUBROUTINE art_read_pol_oper
!!
!!
END MODULE mo_art_read_pol_oper
