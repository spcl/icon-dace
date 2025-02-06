!
! mo_art_external_init_volc
! This module initializes the volcanoes
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

MODULE mo_art_external_init_volc
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_math_constants,                ONLY: pi
  USE mo_model_domain,                  ONLY: t_patch
  USE mo_parallel_config,               ONLY: nproma,p_test_run
  USE mo_gnat_gridsearch,               ONLY: gnat_init_grid, gnat_destroy, t_gnat_tree, &
    &                                         gnat_query_containing_triangles,           &
    &                                         gnat_merge_distributed_queries, gk
  USE mo_grid_config,                   ONLY: grid_sphere_radius
  USE mo_io_units,                      ONLY: find_next_free_unit
  USE mo_exception,                     ONLY: message_text, finish
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_external_types,            ONLY: t_art_volcdata
  USE mo_art_emiss_types,               ONLY: t_art_emiss2tracer
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: art_extinit_volcanoes
  
CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_extinit_volcanoes(p_patch,volclist_file,dict_tracer, &
  &                              emiss2tracer,volc_data,ash_scheme)
!<
! SUBROUTINE art_extinit_volcanoes
! The subroutine initializes the volcanoes
! Part of Module: mo_art_external_init_volc
! Based on: -
! Author: Daniel Rieger, KIT
! Initial Release: 2015-02-06
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
  TYPE(t_patch),INTENT(in)           :: &
    &  p_patch
  CHARACTER(LEN=*),INTENT(in)        :: &
    &  volclist_file       !< Path and name of file including the volcano list
  TYPE(t_key_value_store),INTENT(in) :: &
    &  dict_tracer                        !< Dictionary for tracer indizes
  TYPE(t_art_emiss2tracer),INTENT(IN)      :: &
    &  emiss2tracer    !< Aerosol tracer to emission dictionary
  TYPE(t_art_volcdata),INTENT(inout) :: &
    &  volc_data           !< Volcano data storage type
  INTEGER, INTENT(in)                :: &
    &  ash_scheme          !< 1 for bins, 2 for modes
  ! Local Variables
  TYPE(t_gnat_tree)   :: &
    &  gnat                !< Necessary for getting volcano locations via gnat grid search
  REAL(wp)            :: &
    &  pi_180              !< PI / 180
  INTEGER             :: &
    &  nlev,             & !< Number of vertical levels in current domain
    &  nvolc_curr,       & !< Total number of volcanoes in input file
    &  jvolcano,         & !< Loop index for nvolc_curr
    &  iunit,iostat,     & !< i/o unit number of input file / error variable for i/o
    &  idx_slash,        & !< Index of last slash in volclist_file, needed to get path
    &  pos(16),          & !< Position of information in line_volcdata_trim (starting character)
    &  ilength(16),      & !< pos(:) + ilength(16) = ending character of information 
                           !    in line_volcdata_trim
    &  nwords,           & !< Number of words (i.e. information) in line_volcdata_trim
    &  gnat_nblks,       & !< Number of blocks
    &  gnat_npromz,      & !< Number of proma
    &  gnat_jb, gnat_jc, & !< Loop index for blocks / proma
    &  iasha, iashb,     & !< ash index mode a/b
    &  iashc, iash1,     & !< ash index mode c/bin 1
    &  iash2, iash3,     & !< ash index bin 2/3
    &  iash4, iash5,     & !< ash index bin 4/5
    &  iash6, ierror,    & !< ash index bin 6/ error return value of mo_key_value_store
    &  min_loc_tracer,   & !< Lowest index of ash in tracer container
    &  max_loc_tracer      !< Highest index of ash in tracer container
  CHARACTER(LEN=MAX_CHAR_LENGTH) :: &
    &  thisroutine = "mo_art_external_init_volc:art_extinit_volcanoes"
  CHARACTER(LEN=255)  :: &
    &  line_volcdata,    & !< One line of volcanodata is read into this string
    &  line_volcdata_trim  !< TRIM(line_volcdata)
  CHARACTER(LEN=10)   :: &
    &  latitude,         & !< latitude of volcano
    &  longitude,        & !< longitude of volcano
    &  elev                !< elevation of volcano
  REAL(gk),ALLOCATABLE :: &
    &  in_points(:,:,:), & !< geographical locations
    &  min_dist(:,:)       !< minimal distance
  INTEGER,ALLOCATABLE :: &
    &  tri_idx(:,:,:)      !< tri_idx(2,nproma,(maxvolcs/nproma+1))
  
  ! ----------------------------------
  ! Initializations
  ! ----------------------------------
  
  pi_180     = pi/180._wp
  nvolc_curr = 0
  nlev       = p_patch%nlev
  ALLOCATE(volc_data%info(volc_data%maxvolcs))
  
  SELECT CASE(ash_scheme)
    CASE(1)
      CALL dict_tracer%get('ash1',iash1,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash1 not found in dictionary.')
      CALL dict_tracer%get('ash2',iash2,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash2 not found in dictionary.')
      CALL dict_tracer%get('ash3',iash3,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash3 not found in dictionary.')
      CALL dict_tracer%get('ash4',iash4,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash4 not found in dictionary.')
      CALL dict_tracer%get('ash5',iash5,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash5 not found in dictionary.')
      CALL dict_tracer%get('ash6',iash6,ierror)
        IF(ierror /= SUCCESS) CALL finish (thisroutine, 'iash6 not found in dictionary.')
    CASE(2)
      IF (emiss2tracer%lcalcemiss) THEN
        iasha = emiss2tracer%itr3(1,1) 
        iashb = emiss2tracer%itr3(2,1) 
        iashc = emiss2tracer%itr3(3,1) 
      ELSE
        CALL finish(thisroutine, &
          &         'For iart_volcano=2, all emissions need to be specified in tracer XML.')
      ENDIF
    CASE DEFAULT
      CALL finish (thisroutine, 'Unknown volcanic ash scheme.')
  END SELECT
  
  ! ----------------------------------
  ! Open volcano file, read and store data
  ! ----------------------------------
  
  iunit = find_next_free_unit(100,1000)
  IF (iunit < 0) THEN
    WRITE (message_text,'(a,a)') TRIM(volclist_file), " -- Failed call to find_next_free_unit."
    CALL finish (thisroutine, message_text)
  END IF
  
  OPEN(iunit, file= TRIM(volclist_file), FORM='FORMATTED',  &
    &  STATUS='OLD', ACTION='READ', IOSTAT=iostat)
    IF (iostat /= 0) THEN
      WRITE (message_text,'(a,a,a)') 'ERROR *** opening file ', TRIM(volclist_file), ' failed ***'
      CALL finish (thisroutine, message_text)
    ENDIF
  
  idx_slash = SCAN( TRIM(volclist_file), "/", .TRUE. )
  
  IF ( idx_slash == 0 ) THEN
    volc_data%volcanofile_path = "./"
  ELSE
    volc_data%volcanofile_path = volclist_file(1:idx_slash)
  END IF

  ! READ VOLCANO DATA AND SET DEFAULT ERUPTION PARAMETERS
  volclist: DO
    ! read line of volcano list
    READ(iunit, '(A255)', IOSTAT=iostat) line_volcdata
    
    line_volcdata_trim = TRIM(line_volcdata)
    IF (iostat /= 0) THEN
      EXIT volclist  ! end of file is reached
    ELSE
      CALL split_volcstring(line_volcdata_trim, nwords, pos, ilength)
      IF ( nwords <= 15 ) CYCLE volclist
      ! complete volcano specification line found -- increase number of volcanoes by one
      nvolc_curr = nvolc_curr + 1
      jvolcano = nvolc_curr
      ! Get all data available in the volcano database.
      ! By default only the latitude, longitude and eruption type will be used
      ! following 'Preliminary Spreadsheet of Eruption Source 
      !                       Parameters for Volcanoes of the World'
      volc_data%info(jvolcano)%id      = line_volcdata_trim(pos(1):pos(1)+ilength(1))
      volc_data%info(jvolcano)%rn      = line_volcdata_trim(pos(2):pos(2)+ilength(2))
      volc_data%info(jvolcano)%sn      = line_volcdata_trim(pos(3):pos(3)+ilength(3))
      volc_data%info(jvolcano)%vn      = line_volcdata_trim(pos(4):pos(4)+ilength(4))
      volc_data%info(jvolcano)%name    = line_volcdata_trim(pos(5):pos(5)+ilength(5))
      volc_data%info(jvolcano)%location= line_volcdata_trim(pos(6):pos(6)+ilength(6))
      volc_data%info(jvolcano)%status  = line_volcdata_trim(pos(7):pos(7)+ilength(7))
      latitude                         = line_volcdata_trim(pos(8):pos(8)+ilength(8))
      READ(latitude,*) volc_data%info(jvolcano)%latitude
      volc_data%info(jvolcano)%ns      = line_volcdata_trim(pos(9):pos(9)+ilength(9))
      IF (volc_data%info(jvolcano)%ns == 'S') THEN
        volc_data%info(jvolcano)%latitude = - volc_data%info(jvolcano)%latitude
      ENDIF
      volc_data%info(jvolcano)%vf      = line_volcdata_trim(pos(10):pos(10)+ilength(10))
      longitude                        = line_volcdata_trim(pos(11):pos(11)+ilength(11))
      READ(longitude,*) volc_data%info(jvolcano)%longitude
      volc_data%info(jvolcano)%ew      = line_volcdata_trim(pos(12):pos(12)+ilength(12))
      IF (volc_data%info(jvolcano)%ew == 'W') THEN
        volc_data%info(jvolcano)%longitude = - volc_data%info(jvolcano)%longitude
      ENDIF
      elev                             = line_volcdata_trim(pos(13):pos(13)+ilength(13))
      IF (elev == 'X' ) THEN
        volc_data%info(jvolcano)%elev = -1000.0_wp ! elevation not available DO NOT USE
      ELSE
        READ(elev,*) volc_data%info(jvolcano)%elev
      END IF
      volc_data%info(jvolcano)%type          = line_volcdata_trim(pos(14):pos(14)+ilength(14))
      volc_data%info(jvolcano)%timeframe     = line_volcdata_trim(pos(15):pos(15)+ilength(15))
      volc_data%info(jvolcano)%eruption_type = line_volcdata_trim(pos(16):pos(16)+ilength(16))
      
      ! ----------------------------------
      ! Default vaules: Distribute ash emission to the individual size classes
      !                 based on aircraft measurements
      ! ----------------------------------
      IF (ash_scheme == 1) THEN
        min_loc_tracer = MINVAL((/iash1,iash2,iash3,iash4,iash5,iash6/))
        max_loc_tracer = MAXVAL((/iash1,iash2,iash3,iash4,iash5,iash6/))
        ALLOCATE( volc_data%info(jvolcano)%cfactor(min_loc_tracer:max_loc_tracer))
        volc_data%info(jvolcano)%cfactor(iash1) = 0.014884_wp
        volc_data%info(jvolcano)%cfactor(iash2) = 0.080372_wp
        volc_data%info(jvolcano)%cfactor(iash3) = 0.186047_wp
        volc_data%info(jvolcano)%cfactor(iash4) = 0.372093_wp
        volc_data%info(jvolcano)%cfactor(iash5) = 0.226047_wp
        volc_data%info(jvolcano)%cfactor(iash6) = 0.120558_wp
      END IF
      IF (ash_scheme == 2) THEN
        min_loc_tracer = MINVAL((/iasha,iashb,iashc/))
        max_loc_tracer = MAXVAL((/iasha,iashb,iashc/))
        ALLOCATE( volc_data%info(jvolcano)%cfactor(min_loc_tracer:max_loc_tracer))
        volc_data%info(jvolcano)%cfactor(iasha) = 0.0585_wp
        volc_data%info(jvolcano)%cfactor(iashb) = 0.1661_wp
        volc_data%info(jvolcano)%cfactor(iashc) = 0.7754_wp
      END IF
      
      
      ! ----------------------------------
      !  Default vaules: Mastin et al. 2009 eruption rate / height above vent relation
      ! ----------------------------------
      volc_data%info(jvolcano)%mastin_fac = 3.295_wp
      volc_data%info(jvolcano)%mastin_exp =  4.15_wp

      ! --------- ASSIGN DEFAULT ERUPTION PARAMETERS BASED ON THE ERUPTION TYPE -------------------

      SELECT CASE(TRIM(volc_data%info(jvolcano)%eruption_type))
        CASE ('M0')
          ! Mafic, Standard
          volc_data%info(jvolcano)%height_above_vent = 7000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =   60.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =  1.0e5_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =   0.01_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =   0.05_wp ! [] mass fraction of ash smaller
                                                                 !   than 63 mu m
        CASE ('M1')
          ! Mafic, Small
          volc_data%info(jvolcano)%height_above_vent = 2000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =  100.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =  5.0e3_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =  0.001_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =   0.02_wp ! [] mass fraction of ash smaller
                                                                 !   than 63 mu m
        CASE ('M2')
          ! Mafic, Medium
          volc_data%info(jvolcano)%height_above_vent = 7000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =   60.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =  1.0e5_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =   0.01_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =   0.05_wp ! [] mass fraction of ash smaller
                                                                 !   than 63 mu m
        CASE ('M3')
          ! Mafic, Large
          volc_data%info(jvolcano)%height_above_vent = 10000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =     5.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =   1.0e6_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =    0.17_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =     0.1_wp ! [] mass fraction of ash smaller
                                                                  !   than 63 mu m
        CASE ('S0')
          ! Silicic, Standard
          volc_data%info(jvolcano)%height_above_vent = 11000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =     3.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =   4.0e6_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =   0.015_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =     0.4_wp ! [] mass fraction of ash smaller
                                                                  !   than 63 mu m
        CASE ('S1')
          ! Silicic, Small
          volc_data%info(jvolcano)%height_above_vent = 5000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =   12.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =  2.0e5_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =  0.003_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =    0.1_wp ! [] mass fraction of ash smaller
                                                                 !   than 63 mu m
        CASE ('S2')
          ! Silicic, Medium
          volc_data%info(jvolcano)%height_above_vent = 11000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =     3.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =   4.0e6_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =   0.015_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =     0.4_wp ! [] mass fraction of ash smaller
                                                                  !   than 63 mu m
        CASE ('S3')
          ! Silicic, Large
          volc_data%info(jvolcano)%height_above_vent = 15000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =     8.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =   1.0e7_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =    0.15_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =     0.5_wp ! [] mass fraction of ash smaller
                                                                  !   than 63 mu m
        CASE ('S8')
          !  Co-ignimbrite cloud
          volc_data%info(jvolcano)%height_above_vent = 25000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =     0.5_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =   1.0e8_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =    0.05_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =     0.5_wp ! [] mass fraction of ash smaller
                                                                  !   than 63 mu m
        CASE ('S9')
          ! Brief
          volc_data%info(jvolcano)%height_above_vent = 10000.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     =    0.01_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     =   3.0e6_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   =  0.0003_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      =     0.6_wp ! [] mass fraction of ash smaller
                                                                  !   than 63 mu m
        CASE ('U0')
          ! Submarine
          volc_data%info(jvolcano)%height_above_vent = 0.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     = 0.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     = 0.0_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   = 0.0_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      = 0.0_wp ! [] mass fraction of ash smaller
                                                              !   than 63 mu m
        CASE default
          ! ERROR?
          volc_data%info(jvolcano)%height_above_vent = 0.0_wp ! [m]
          volc_data%info(jvolcano)%erup_duration     = 0.0_wp ! [h]
          volc_data%info(jvolcano)%eruption_rate     = 0.0_wp ! [kg/s]
          volc_data%info(jvolcano)%eruption_volume   = 0.0_wp ! [km^3]
          volc_data%info(jvolcano)%ash_fraction      = 0.0_wp ! [] mass fraction of ash smaller
                                                              !   than 63 mu m
      END SELECT
      
      volc_data%info(jvolcano)%eruption_rate = volc_data%info(jvolcano)%eruption_rate  &
        &                                  * volc_data%info(jvolcano)%ash_fraction      
                                            ! only ash emission (particles < 63 microns)
    END IF ! iostat /= 0
  ENDDO volclist
  
  CLOSE(iunit)
  
  ! ----------------------------------
  ! Locate volcanoes on ICON grid
  ! ----------------------------------
  
  gnat_nblks  = nvolc_curr/nproma + 1
  gnat_npromz = nvolc_curr - nproma*(gnat_nblks-1)

  ALLOCATE(tri_idx(2,nproma,gnat_nblks))
  ALLOCATE(volc_data%owner(nvolc_curr))
  ALLOCATE(in_points(nproma,gnat_nblks,2))
  ALLOCATE(min_dist(nproma,gnat_nblks))
  ALLOCATE(volc_data%lat_idx(nproma,gnat_nblks) )
  ALLOCATE(volc_data%lon_idx(nproma,gnat_nblks) )
  ALLOCATE(volc_data%height_factor(nlev,nvolc_curr))
  volc_data%height_factor = 0.0_wp

  !----------------
  ! Fill in_points array for finding grid indices
  !----------------

  gnat_jc = 0
  gnat_jb = 1
  DO jvolcano=1, nvolc_curr
    gnat_jc=gnat_jc+1
    IF (gnat_jc>nproma) THEN
      gnat_jc=1
      gnat_jb=gnat_jb+1
    ENDIF

    volc_data%lat_idx(gnat_jc,gnat_jb) = volc_data%info(jvolcano)%latitude
    volc_data%lon_idx(gnat_jc,gnat_jb) = volc_data%info(jvolcano)%longitude

    in_points(gnat_jc,gnat_jb,1) = volc_data%lon_idx(gnat_jc,gnat_jb) * pi_180
    in_points(gnat_jc,gnat_jb,2) = volc_data%lat_idx(gnat_jc,gnat_jb) * pi_180
  END DO
  
  !----------------
  ! Calculating grid indeces -> fill tri_idx
  !----------------
  
  ! build GNAT data structure
  CALL gnat_init_grid(gnat, p_patch)

  ! perform proximity query
  CALL gnat_query_containing_triangles(gnat, p_patch, in_points(:,:,:),  &
    &                                  nproma, gnat_nblks, gnat_npromz,  &
    &                                  grid_sphere_radius,p_test_run,    &
    &                                  tri_idx(:,:,:), min_dist(:,:))

  CALL gnat_merge_distributed_queries(p_patch, nvolc_curr, nproma, gnat_nblks, min_dist,  &
    &                                 tri_idx(:,:,:), in_points(:,:,:),                   &
    &                                 volc_data%owner(:), volc_data%ithis_nlocal_pts)

  ! clean up
  CALL gnat_destroy(gnat)


  gnat_jc = 0
  gnat_jb = 1
  DO jvolcano = 1, nvolc_curr
    gnat_jc=gnat_jc+1
    IF (gnat_jc>nproma) THEN
      gnat_jc=1
      gnat_jb=gnat_jb+1
    ENDIF
    volc_data%info(jvolcano)%tri_iidx_loc=tri_idx(1,gnat_jc,gnat_jb)
    volc_data%info(jvolcano)%tri_iblk_loc=tri_idx(2,gnat_jc,gnat_jb)
  END DO
  
END SUBROUTINE art_extinit_volcanoes
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE split_volcstring(zline, n, pos, ilength)
!<
! SUBROUTINE split_volcstring
! The subroutine splits the string from the input file
! Part of Module: mo_art_init_volc
! Based on: -
! Author: Kristina Lundgren, KIT
! Initial Release: 2012-01-30
! Modifications:
! YYYY-MM-DD: <name>, <instituttion>
! - ...
!>
IMPLICIT NONE

  CHARACTER(len=*), INTENT(in)      :: zline              ! string containing list
  INTEGER,          INTENT(out)     :: n                  ! number of parts
  INTEGER,          INTENT(inout)   :: pos(:), ilength(:) ! position, lengths of parts
  ! local variables
  INTEGER       :: i           ! index position
  LOGICAL       :: l_word_open ! flag. if true, index "i" is part of a word
  INTEGER       :: istart

  l_word_open = .FALSE.
  n           = 0
  istart      = 1
  DO i = 1, LEN(zline)
    !                 space                         tab
    IF (.NOT. (IACHAR(zline(i:i)) == 32 .OR. IACHAR(zline(i:i)) == 9) ) THEN
      l_word_open = .TRUE.
    ELSE
      IF (l_word_open) THEN
        n = n + 1
        pos(n)  = istart
        ilength(n) = LEN( TRIM( zline(istart:(i-1)) ) )
      END IF
      istart = i+1
      l_word_open = .FALSE.
    END IF
  END DO

END SUBROUTINE split_volcstring
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_external_init_volc
