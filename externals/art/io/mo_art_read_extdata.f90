!
! mo_art_read_extdata
! This module provides is a wrapper to read soildata
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

MODULE mo_art_read_extdata
! ICON
  USE mo_kind,                          ONLY: wp
  USE mtime,                            ONLY: datetime, datetimeToPosixString, &
     &                                        MAX_DATETIME_STR_LEN
  USE mo_exception,                     ONLY: message
  USE mo_util_cdi,                      ONLY: t_inputParameters
  USE mo_model_domain,                  ONLY: p_patch
  USE mo_read_interface,                ONLY: t_stream_id, openInputFile, closeFile,  &
    &                                     read_2d_1time, read_2d, on_cells
  USE mo_io_config,                     ONLY: default_read_method
! ART
  USE mo_art_data,                      ONLY: p_art_data
  USE mo_art_atmo_data,                 ONLY: t_art_atmo
  USE mo_art_wrapper_routines,          ONLY: art_get_indices_c
  USE mo_art_config,                    ONLY: art_config !! TEMPORARY!!!
  USE mo_art_read,                      ONLY: art_read, art_open_cdi, art_close_cdi
  USE mo_art_external_types,            ONLY: t_art_soil_properties, t_art_soil_table, &
    &                                         t_art_pollen_table, t_art_pollen_properties, &
    &                                         t_art_online_dms, t_art_biomBurn_properties
  USE mo_art_create_filenames,          ONLY: art_create_filenames, art_check_filename
  USE mo_art_bvoc,                      ONLY: pftn
  USE mo_art_io_constants,              ONLY: IEXT_BIO, IEXT_SOIL, IEXT_DMS,  &
    &                                         IEXT_POCON, IEXT_POAMB,         &
    &                                         IART_FILENAMELEN, IEXT_BCF
  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER :: routine = 'mo_art_read_extdata'

  PUBLIC :: art_read_ext_soil, art_read_PFTdata, art_read_pollendata, &
          & art_read_sdes_ambrosia
  PUBLIC :: art_read_ext_dms, art_read_biomBurndata

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_ext_soil(jg,soil_prop,input_folder)
!<
! SUBROUTINE art_read_ext_soil
! Read soil particle size distribution
! Part of Module: mo_art_read_extdata
! Author: Daniel Rieger, KIT
! Initial Release: 2014-05-14
! Modifications:
! 2016-11-16: Daniel Rieger, KIT
! - Adaption to new mo_art_external_state module, reduced dependencies
!>
  INTEGER, INTENT(in)                    :: &
    &  jg                                     !< patch on which computation is performed
  TYPE(t_art_soil_properties),INTENT(inout) :: &
    &  soil_prop                              !< pointer to soil properties (stored in p_art_data)
  CHARACTER(LEN=*),INTENT(in)            :: &
    &  input_folder                           !< Folder containing the soil data
  ! Local variables:
  TYPE(t_inputParameters)             :: &
    &  cdi_param                           !< Parameters for read_cdi call
  TYPE(t_art_soil_table),POINTER      :: &
    &  this_soil_table                     !< pointer to soil table
  INTEGER                             :: &
    &  cdi_artdataset_id                   !< CDI stream ID (for each domain)
  INTEGER                             :: &
    &  n                                   !< counter
  CHARACTER(LEN=IART_FILENAMELEN)     :: &
    &  soildataset                         !< Path and filename of soil dataset
  CHARACTER(LEN=7)                    :: &
    &  vname                               !< Name of current variable to be read

  CALL art_create_filenames(jg, TRIM(input_folder), IEXT_SOIL, soildataset, .FALSE.)
! JF:   CALL art_check_filename(soildataset)

  CALL art_open_cdi(jg,soildataset,cdi_artdataset_id,cdi_param)

  ! Loop through soil types and read external data
  DO n=1, soil_prop%nsoil_types
    this_soil_table => soil_prop%soil_type(n)
    vname = 'fr_'//TRIM(this_soil_table%sname)
    CALL message(TRIM(routine)//':art_read_ext_soil','Reading soil data for soiltype: '//TRIM(vname))
    CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vname), &
      &           this_soil_table%fr_soil)
  ENDDO

  CALL art_close_cdi(cdi_artdataset_id, cdi_param)
  NULLIFY (this_soil_table)

END SUBROUTINE art_read_ext_soil
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_pollendata(jg, p_pollen_prop, input_folder)
!<
! SUBROUTINE art_read_pollendata
! Read distribution of pollen plants and parameter for altitude correction 
! Part of Module: mo_art_read_extdata
! Author: Jonas Straub, KIT
! Initial Release: 2016-07-21
! Modifications:
! 2016-12-21: Jonas Straub, KIT
! - Adaption to new structure
!>
  INTEGER, INTENT(in)                       :: &
    &  jg                                        !< patch on which computation is performed
  TYPE(t_art_pollen_properties), INTENT(inout) :: &
    &  p_pollen_prop                             !< pointer to pollen properties (stored in p_art_data)
  CHARACTER(LEN=*), INTENT(in)              :: &
    &  input_folder                              !< Folder containing the pollen data
  ! Local variables:
  TYPE(t_inputParameters)                   :: &
    &  cdi_param                                 !< Parameters for read_cdi call
  TYPE(t_art_pollen_table), POINTER         :: &
    &  this_pollen_table                         !< pointer to pollen table
  INTEGER                                   :: &
    &  cdi_artdataset_id                         !< CDI stream ID (for each domain)
  INTEGER                                   :: &
    &  n, iv                                      !< counter
  CHARACTER(LEN=IART_FILENAMELEN)           :: &
    &  pollendataset                             !< Path and filename of pollen coverage dataset and altitude correction
  CHARACTER(LEN=30)                         :: &
    &  vpref,                                  & !< Name of coverage variable and altitude correction to be read
    &  vname                                  
  CHARACTER(LEN=10), ALLOCATABLE            :: &
    &  vsuff(:)                                !< Name of coverage variable and altitude correction to be read


  ! Loop through pollen types and read plant coverage
  DO n=1, p_pollen_prop%npollen_types
    this_pollen_table => p_pollen_prop%pollen_type(n)
    vname = this_pollen_table%sname
    vpref = ''
    ! Last pollentype was reached and all data reading is finished
    !_jf: on NEC vector host 'this_pollen_table%fr_cov' is associated?!
    !_jf: as workaround test if variable shortname is set.
    IF (.NOT. ASSOCIATED(this_pollen_table%fr_cov) .OR.  &
      & TRIM(this_pollen_table%shortname) == '') THEN 
      CALL message(TRIM(routine)//': art_read_pollendata', 'Left reading loop')
      EXIT
    END IF 

    IF ((this_pollen_table%linit .EQV. .TRUE.) .AND. (TRIM(vname) /= '' )) THEN
      SELECT CASE(TRIM(vname))
        CASE('pollbetu')
          vpref      = 'BETU' 
          ALLOCATE(vsuff(2))
          vsuff      = (/'fr   ','hcem '/)  ! be aware of the whitespaces        
        CASE('pollalnu')
          vpref      = 'ALNU'
          ALLOCATE(vsuff(2))
          vsuff      = (/'fr   ','hcem '/)  ! be aware of the whitespaces        
        CASE('pollcory')
          vpref      = 'CORY'
          ALLOCATE(vsuff(2))
          vsuff      = (/'fr   ','hcem '/)  ! be aware of the whitespaces        
        CASE('pollpoac')
          vpref      = 'POAC'
          ALLOCATE(vsuff(2))
          vsuff      = (/'fr   ','hcem '/)  ! be aware of the whitespaces        
        CASE('pollambr')
          vpref      = 'AMBR'
          ALLOCATE(vsuff(2))
          vsuff      = (/'fr   ','hcem '/)  ! be aware of the whitespaces        
        CASE DEFAULT
          CALL message(TRIM(routine)//': art_read_pollendata', &
            &          'pollentype is not yet implemented'//TRIM(vname))
          CONTINUE
      END SELECT
    
      CALL message(TRIM(routine)//':art_read_pollendata',      &
        &          'Reading extern data for pollentype: '//TRIM(vname))
      
      CALL art_create_filenames(jg, TRIM(input_folder), IEXT_POCON, pollendataset, &
        &                       .FALSE.)
! JF:       CALL art_check_filename(pollendataset)

      CALL art_open_cdi(jg,pollendataset,cdi_artdataset_id,cdi_param)
    
      DO iv=1, SIZE(vsuff)

        SELECT CASE (TRIM(vsuff(iv)))
          CASE ('fr')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%fr_cov)
          CASE ('hcem')
            CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vpref)//TRIM(vsuff(iv)),  &
              &           this_pollen_table%f_q_alt)
          CASE DEFAULT 
            CALL message(TRIM(routine)//': art_read_pollendata', &
              &          'input field is not defined'//TRIM(vpref)//TRIM(vsuff(iv)))
        END SELECT

      END DO

      CALL art_close_cdi(cdi_artdataset_id, cdi_param)
      DEALLOCATE( vsuff )

    END IF 
  ENDDO

  NULLIFY (this_pollen_table)

END SUBROUTINE art_read_pollendata

SUBROUTINE art_read_sdes_ambrosia(jg, p_pollen_prop, input_folder, current_date)
!<
! SUBROUTINE art_read_sdes_ambrosia
! Read daily sdes files for ambrosia
! Part of Module: mo_art_read_extdata
! Author: Stefan Muthers, DWD
! Initial Release: 2019-12-03
! Modifications:
!>
  INTEGER, INTENT(in)                       :: &
    &  jg                                        !< patch on which computation is performed
  TYPE(t_art_pollen_properties), INTENT(inout) :: &
    &  p_pollen_prop                             !< pointer to pollen properties (stored in p_art_data)
  CHARACTER(LEN=*), INTENT(in)              :: &
    &  input_folder                              !< Folder containing the pollen data
  TYPE(datetime), POINTER, INTENT(in)        :: &
   &   current_date
  ! Local variables:
  TYPE(t_inputParameters)                   :: &
    &  cdi_param                                 !< Parameters for read_cdi call
  TYPE(t_art_pollen_table), POINTER         :: &
    &  this_pollen_table                         !< pointer to pollen table
  INTEGER                                   :: &
    &  cdi_artdataset_id                         !< CDI stream ID (for each domain)
  INTEGER                                   :: &
    &  n                                         !< counter
  CHARACTER(LEN=IART_FILENAMELEN)           :: &
    &  pollendataset                             !< Path and filename of pollen coverage dataset and altitude correction
  CHARACTER(LEN=30)                         :: &
    &  vname, varname                            !< Name of pollen type and variable to be read
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)       :: &
    &  sdatetime                                 !< Date String for the filename


  CALL datetimeToPosixString(current_date, sdatetime, "%Y-%m-%d")

  ! Loop through pollen types and sdes when AMBR is found
  DO n=1, p_pollen_prop%npollen_types
    this_pollen_table => p_pollen_prop%pollen_type(n)
    vname = this_pollen_table%sname
  
    IF ((this_pollen_table%linit .EQV. .TRUE.) .AND. (TRIM(vname) == 'pollambr' )) THEN
          varname      = 'AMBRsdes'
          CALL message(TRIM(routine)//':art_read_pollendata',      &
            &          'Reading extern data for pollentype: '//TRIM(varname))
              
          CALL art_create_filenames(jg, TRIM(input_folder), IEXT_POAMB, pollendataset, &
                        &                       .TRUE., TRIM(sdatetime(1:10)))  ! FIXME: why is 1:10 needed? 
! JF:           CALL art_check_filename(pollendataset)

          CALL art_open_cdi(jg,pollendataset,cdi_artdataset_id,cdi_param)
            
          CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(varname),  &
            &           this_pollen_table%f_q_seas)
          CALL art_close_cdi(cdi_artdataset_id, cdi_param)

          !$ACC UPDATE DEVICE(this_pollen_table%f_q_seas) ASYNC(1)

    END IF  

  ENDDO

  NULLIFY (this_pollen_table)

END SUBROUTINE art_read_sdes_ambrosia
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_PFTdata(pft,jg)
!<
! SUBROUTINE art_read_PFTdata
! Read plant functional type distribution
! Part of Module: mo_art_read_extdata
! Author: Michael Weimer, KIT
! Initial Release: YYYY-MM-DD
! Modifications:
! YYYY-MM-DD: <author>, <institution>
! - <Description>
!>
  IMPLICIT NONE
  REAL(wp), INTENT(out)               :: &
    &  pft(:,:,:)                          !< Plant functional type distribution
  INTEGER, INTENT(in)                 :: &
    &  jg                                  !< domain index
    ! Local variables:
  INTEGER                             :: &
    &  pft_idx                             !< counter
  CHARACTER(LEN=200)                  :: &
    &  pftdataset                          !< Path and filename of pft dataset
  CHARACTER(LEN=2)                    :: &
    &  pft_idx_str                         !< Name of current variable to be read
  TYPE(t_stream_id) ::  &
    &  stream_id                           !< stream id of the netCDF file


  CALL art_create_filenames(jg,TRIM(art_config(jg)%cart_input_folder)//'/PFT', &
    &                       IEXT_BIO,pftdataset,.FALSE.)
  CALL art_check_filename(pftdataset)

  CALL openInputFile( stream_id, TRIM(pftdataset), p_patch(jg), &
    &                 default_read_method)

  DO pft_idx = 1 , pftn
    WRITE(pft_idx_str,'(I2)') pft_idx
    pft_idx_str =  ADJUSTL(pft_idx_str)
    CALL read_2d(stream_id, on_cells,'PFT'//TRIM(pft_idx_str), fill_array= pft(:,:,pft_idx))
  ENDDO

  CALL closeFile(stream_id)
END SUBROUTINE art_read_PFTdata
!!
!!
SUBROUTINE art_read_ext_dms(jg,p_online_dms,input_folder)
!<
! SUBROUTINE art_read_ext_dms
! Read distribution of DMS from the ocean
! Part of Module: mo_art_read_extdata
! Author: Carmen Ullwer, KIT
! Initial Release: YYYY-MM-DD
! Modifications:
! YYYY-MM-DD: <author>, <institution>
! - <Description>
!>
  INTEGER, INTENT(in)                    :: &
    &  jg                                     !< patch on which computation is performed
  TYPE(t_art_online_dms),INTENT(inout)   :: &
    &  p_online_dms                           !< pointer to online dms (stored in p_art_data)
  CHARACTER(LEN=*),INTENT(in)            :: &
    &  input_folder                           !< Folder containing the dms data
  ! Local variables:
  INTEGER                             :: &
    &  n, jc, jb,                        & !< counter
    &  i_startidx, i_endidx
  CHARACTER(LEN=IART_FILENAMELEN)     :: &
    &  dmsdataset                         !< Path and filename of dms dataset
  CHARACTER(LEN=2)                    :: &
    &  dms_month_str
  CHARACTER(LEN=7)                    :: &
    &  vname                               !< Name of current variable to be read
  TYPE(t_stream_id) ::  &
    &  stream_id                           !< stream id of the netCDF file
  TYPE(t_art_atmo), POINTER :: &
    &  art_atmo                            !< POINTER to ART atmo structure

  art_atmo => p_art_data(jg)%atmo

  CALL art_create_filenames(jg, TRIM(input_folder)//'/DMS', IEXT_DMS,   &
             &              dmsdataset, .TRUE., '2010-01-01-00')
  
  CALL art_check_filename(dmsdataset)

  CALL openInputFile( stream_id, TRIM(dmsdataset), p_patch(jg), &
    &                 default_read_method)

  DO n=1, p_online_dms%ndms_months
    WRITE(dms_month_str,'(I2)') n
    dms_month_str = ADJUSTL(dms_month_str)
    vname = 'DMSe_'//TRIM(dms_month_str)
    CALL message(TRIM(routine)//':art_read_ext_dms',                            &
           &     'Reading DMS data for DMS online calculation: '//TRIM(vname))

    CALL read_2d_1time(stream_id, on_cells,TRIM(vname), fill_array= p_online_dms%dms_month(:,:,n))
    ! There are negative values on land points of the dataset as missing values.
    ! As they are used as a factor in the calculation of online DMS emissions,
    ! set it to 0
    DO jb = art_atmo%i_startblk,art_atmo%i_endblk
      CALL art_get_indices_c(jg, jb,i_startidx, i_endidx)

      DO jc = i_startidx,i_endidx
        IF (p_online_dms%dms_month(jc,jb,n) < 0._wp) THEN
          p_online_dms%dms_month(jc,jb,n) = 0._wp
        END IF
      END DO ! j
    END DO ! i
  ENDDO  ! n

  CALL closeFile(stream_id)

END SUBROUTINE art_read_ext_dms
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_read_biomBurndata(jg,p_biomBurn_prop, input_folder)
!<
! SUBROUTINE art_read_biomBurndata
! Read GFAS data as an proxy for wildfire locations 
! Part of Module: mo_art_read_extdata
! Author: Jonas Straub, KIT
! Initial Release: 2018-02-22
! Modifications:
! YYYY-MM-DD: <author>, <institution>
! - <Description>
!>
  INTEGER, INTENT(IN)                       :: &
    &  jg                                        !< domain index
  TYPE(t_art_biomBurn_properties), INTENT(inout) :: &
    &  p_biomBurn_prop                           !< pointer to biomass burning properties 
                                                 ! (stored in p_art_data)
  CHARACTER(LEN=*),INTENT(in)               :: &
    &  input_folder                              !< Folder containing the biomass burning data
  ! Local variables:
  TYPE(t_inputParameters)                   :: &
    &  cdi_param                                 !< Parameters for read_cdi call
  INTEGER                                   :: &
    &  cdi_artdataset_id                         !< CDI stream ID (for each domain)
  CHARACTER(LEN=IART_FILENAMELEN)           :: &
    &  biomBurndataset                           !< Path and filename of GFAS-data
  CHARACTER(LEN=13)                         :: &
    &  vname                                     !< Name of fire variable

  CALL art_create_filenames(jg, TRIM(input_folder), IEXT_BCF, biomBurndataset, .FALSE.)
  CALL art_check_filename(biomBurndataset)

  CALL art_open_cdi(jg,biomBurndataset,cdi_artdataset_id,cdi_param)

  vname = 'bcfire'
  CALL message(TRIM(routine)//':art_read_biomBurndata',                               &
    &          'Reading flux of black carbon as wildfire locations: '//TRIM(vname))
  CALL art_read(jg, cdi_artdataset_id, cdi_param, TRIM(vname),         &
    &           p_biomBurn_prop%flux_bc)

  CALL art_close_cdi(cdi_artdataset_id, cdi_param)

END SUBROUTINE art_read_biomBurndata
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_read_extdata
