!
! mo_art_emission_volc_2mom
! This module provides the routine for calculating the emission of volcanic ash
! Based on: Rieger et al. (2015): ICON-ART 1.0 - a new online-coupled model system
!                                 from the global to regional scale
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

MODULE mo_art_emission_volc_2mom
! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_io_units,                      ONLY: find_next_free_unit
  USE mo_exception,                     ONLY: message, message_text, finish
  USE mo_impl_constants,                ONLY: MAX_CHAR_LENGTH, SUCCESS
  USE mtime,                            ONLY: datetime, newdatetime, OPERATOR(<=),      &
    &                                         deallocatedatetime, MAX_DATETIME_STR_LEN
  USE mo_key_value_store,               ONLY: t_key_value_store
! ART
  USE mo_art_external_types,            ONLY: t_art_volcdata
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC   :: art_prepare_emission_volc
  PUBLIC   :: art_calculate_emission_volc

CONTAINS
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE art_prepare_emission_volc(current_date, jb, nlev, hlev, dict_tracer, volc_data, &
  &                                  iasha, iashb, iashc)
!<
! SUBROUTINE art_prepare_emission_volc
! Main routine for volcanic ash emissions
! Based on: Rieger et al. (2015): ICON-ART 1.0 - a new online-coupled model system 
!                                 from the global to regional scale
! Part of Module: mo_art_emission_volc_2mom
! Author: Kristina Lundgren, KIT
! Initial Release: 2012-01-30
! Modifications:
! 2015-05-13: Carolin Walter, KIT
! - Implementation of a 2MOM version
! 2015-07-16: Carolin Walter, KIT
! - bugfix: oberserved height = top_of_plume = height_above_vent + elevation
! 2017-04-19: Carolin Walter, KIT
! - bugfix: height_above_vent needs to be set to the value from volcspec file even when it is 0.0
!>
  TYPE(datetime), INTENT(IN), POINTER :: &
    &  current_date                         !< Current date in mtime format
  INTEGER, INTENT(IN)                 :: &  !<
    &  jb,                               &  !< block index
    &  nlev                                 !< number of full levels
  REAL(wp), INTENT(IN)                :: &  !<
    &  hlev(:,:)                            !< height of level (=z_ifc)
  TYPE(t_key_value_store),INTENT(IN)          :: &
    &  dict_tracer                          !< Dictionary for tracer indizes
  TYPE(t_art_volcdata), INTENT(inout) :: &
    &  volc_data              !< Volcano data storage type
  INTEGER, INTENT(IN)   :: &
    &  iasha, iashb, iashc
!Local variables
  INTEGER               :: &
    &  jk,                 & !< loop index over vertical levels
    &  cfactor_check,      & !<
    &  iunit,iostat,       & !< unit and status of ASCII-input file
    &  jvolcano,           & !< loop index over number of volcanoes
    &  jlocal,             & !< loop index over local points
    &  pos(16),            & !<
    &  ilength(16),        & !<
    &  idx_word, nwords,   & !<
    &  idx_hash
  INTEGER, POINTER      :: &
    &  idx,blk
  LOGICAL               :: & 
    &  calc_erup_rate,     & !< calculate emiss rate or use database
    &  lerror_volcspec,    & !<
    &  lnew_volcspec         !<
  CHARACTER(LEN=255)    :: &
    &  line_volcdata,      & !< One line of volcanodata is read into this string
    &  line_volcdata_trim    !< TRIM(line_volcdata)
  CHARACTER(LEN=MAX_CHAR_LENGTH) :: &
    &  volcspec_file,      & !<
    &  thisroutine           !<
  TYPE(datetime), POINTER :: &
    &  mtime_inp_date        !< Input date in mtime format
  CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: &
    &  mtime_inp_datetime_str, & !<
    &  mtime_inp_datetime_str_save
  CHARACTER(LEN=5)      :: &
    &  volcspec_keyword      !<
  CHARACTER(LEN=20)     :: &
    &  volcspec_value        !<
  CHARACTER(LEN=256)    :: &
    &  volcspec_text         !<
  REAL(wp)              :: &
    &  cfactor_sum,        & !<
    &  frac_save,          & !< 
    &  mfac_save,          & !< 
    &  mexp_save,          & !< 
    &  top_of_plume          !< 
  REAL(wp), PARAMETER   :: &
    &  cfactor_eps = 0.000002_wp !<
  REAL(wp), ALLOCATABLE :: &
    &  cfac_save(:)
  
  ! ----------------------------------
  ! ---  Start Routine
  ! ----------------------------------
  thisroutine = "mo_art_emission_volc_2mom:art_prepare_emission_volc"
  calc_erup_rate = .FALSE.
  ALLOCATE(cfac_save(iasha:iashc))

  !-----------------------------------
  ! -- Read specific input file for each volcano
  ! -- Calculate vertical emission distribution for each volcano
  !-----------------------------------
  DO jlocal = 1, volc_data%ithis_nlocal_pts
    jvolcano = volc_data%owner(jlocal)
    idx => volc_data%info(jlocal)%tri_iidx_loc
    blk => volc_data%info(jlocal)%tri_iblk_loc
    IF (blk == jb) THEN 
      volcspec_file = "volcspec_2mom_"//TRIM(volc_data%info(jvolcano)%NAME)//".txt"
      iunit = find_next_free_unit(100,1000)
      IF (iunit < 0) THEN
        WRITE (message_text,'(a,a,a)') TRIM(volc_data%volcanofile_path), TRIM(volcspec_file), &
          &  " -- Failed call to find_next_free_unit."
        CALL finish (thisroutine, message_text)
      ENDIF
      
      OPEN(iunit, file=TRIM(volc_data%volcanofile_path)//TRIM(volcspec_file), form='FORMATTED',  &
        &  status='OLD', action='READ', iostat=iostat )
      
      ! if file is present...
      IF (iostat == 0) THEN
        lerror_volcspec = .FALSE.
        lnew_volcspec   = .FALSE.

        volcspec: DO
          ! read line of specific input data for volcano
          READ(iunit, '(A255)', IOSTAT=iostat) line_volcdata
          IF (iostat /= 0) THEN
            EXIT volcspec  ! end of file is reached
          ELSE
            line_volcdata_trim = TRIM(ADJUSTL(line_volcdata))
            idx_hash = SCAN( line_volcdata_trim, "#" )
            IF ( idx_hash /= 1 ) THEN
              CALL split_volcstring(line_volcdata_trim, nwords, pos, ilength)
              IF ( nwords <= 1 ) CYCLE volcspec
            ELSE
              CYCLE volcspec
            ENDIF

            mtime_inp_datetime_str =  line_volcdata_trim(pos(1):pos(1)+ilength(1))
            mtime_inp_date         => newdatetime(TRIM(mtime_inp_datetime_str))

            IF ( mtime_inp_date <= current_date ) THEN  ! date check
              lerror_volcspec = .FALSE.
              lnew_volcspec   = .TRUE.
              mtime_inp_datetime_str_save = mtime_inp_datetime_str
              calc_erup_rate = .TRUE.
      
              ! set new height above vent
              volcspec_value = line_volcdata_trim(pos(2):pos(2)+ilength(2))
              READ(volcspec_value,*) top_of_plume
              IF (top_of_plume > volc_data%info(jvolcano)%elev) THEN
                volc_data%info(jvolcano)%height_above_vent =  &
                  &     top_of_plume - volc_data%info(jvolcano)%elev
              ELSE
                volc_data%info(jvolcano)%height_above_vent = 0.0
              END IF
      
              WRITE(volcspec_text,'(a,f7.1)') TRIM(volcspec_file)//": "//  &
                &   TRIM(mtime_inp_datetime_str)//"  ", top_of_plume
      
              ! optional values are present
              IF (nwords >= 3 ) THEN
      
                ! save previous values
                frac_save = volc_data%info(jvolcano)%ash_fraction
                cfac_save(iasha:iashc) = volc_data%info(jvolcano)%cfactor(iasha:iashc)
                mfac_save = volc_data%info(jvolcano)%mastin_fac
                mexp_save = volc_data%info(jvolcano)%mastin_exp
                cfactor_check = 0
                cfactor_sum   = 0.0_wp
                ! read optional values
                DO idx_word = 3, nwords
                  volcspec_keyword = line_volcdata_trim(pos(idx_word):pos(idx_word)+4)
                  volcspec_value   = line_volcdata_trim(pos(idx_word)+  &
                    &                                   5:pos(idx_word)+ilength(idx_word))
                  SELECT CASE( volcspec_keyword )
                    CASE('frac=')
                      READ(volcspec_value,*) volc_data%info(jvolcano)%ash_fraction
                    CASE('cfa1=')
                      READ(volcspec_value,*) volc_data%info(jvolcano)%cfactor(iasha)
                      cfactor_check = cfactor_check + 1
                      cfactor_sum   = cfactor_sum + volc_data%info(jvolcano)%cfactor(iasha)
                    CASE('cfa2=')
                      READ(volcspec_value,*) volc_data%info(jvolcano)%cfactor(iashb)
                      cfactor_check = cfactor_check + 10
                      cfactor_sum   = cfactor_sum + volc_data%info(jvolcano)%cfactor(iashb)
                    CASE('cfa3=')
                      READ(volcspec_value,*) volc_data%info(jvolcano)%cfactor(iashc)
                      cfactor_check = cfactor_check + 100
                      cfactor_sum   = cfactor_sum + volc_data%info(jvolcano)%cfactor(iashc)
                    CASE('mfac=')
                      READ(volcspec_value,*) volc_data%info(jvolcano)%mastin_fac
                    CASE('mexp=')
                      READ(volcspec_value,*) volc_data%info(jvolcano)%mastin_exp
                    CASE default
                      WRITE(message_text,'(a)')                                                   &
                        &  "WARNING: undefined keword '"//TRIM(volcspec_keyword)//"' found in "// &
                        &  TRIM(volcspec_file)
                      CALL message(thisroutine, message_text, all_print=.TRUE.)
                      lerror_volcspec = .TRUE.
                  END SELECT
                  IF ( .NOT. lerror_volcspec ) THEN
                    WRITE(volcspec_text,'(a)')                                                    &
                      & TRIM(volcspec_text)//"  "//line_volcdata_trim(pos(idx_word):pos(idx_word) &
                      & + ilength(idx_word))
                  ENDIF
                END DO
                ! check if all cfactor-values were given and if sum equals 1.0
                IF (( cfactor_check /= 0 .AND. cfactor_check /= 111 ) .OR.  &
                  &  ( cfactor_check == 111 .AND. ABS( cfactor_sum - 1.0_wp ) > cfactor_eps )) THEN
                  volc_data%info(jvolcano)%cfactor(iasha:iashc) = cfac_save(iasha:iashc)
                  lerror_volcspec = .TRUE.
                  WRITE(message_text,'(a,f8.6,a)')  &
                    &  "WARNING: uncomplete or wrong (cfactor_sum=", cfactor_sum,      &
                    &  ") specification of cfactor(:) found in "//TRIM(volcspec_file)
                  CALL message(thisroutine, message_text, all_print=.TRUE.)
                ENDIF
      
              ENDIF ! optional values
      
              IF (volc_data%lprint_once) THEN   ! ... do ONLY ONCE
                WRITE(volcspec_text,'(a,f7.1,a,f7.1,4(a,f8.6),2(a,f5.3))')  &
                  & TRIM(volcspec_file)//": "//                             &
                  & TRIM(mtime_inp_datetime_str_save)//"  ",                &
                  & volc_data%info(jvolcano)%height_above_vent,             &
                  & "  elevation=", volc_data%info(jvolcano)%elev,          &
                  & "  frac=", volc_data%info(jvolcano)%ash_fraction,       &
                  & "  cfa1=", volc_data%info(jvolcano)%cfactor(iasha),     &
                  & "  cfa2=", volc_data%info(jvolcano)%cfactor(iashb),     &
                  & "  cfa3=", volc_data%info(jvolcano)%cfactor(iashc),     &
                  & "  mfac=", volc_data%info(jvolcano)%mastin_fac,         &
                  & "  mexp=", volc_data%info(jvolcano)%mastin_exp   
                volc_data%lprint_once = .FALSE.
              ENDIF
      
            ELSE

              CALL deallocatedatetime(mtime_inp_date)
              EXIT volcspec  ! input date > current date found
      
            ENDIF  ! date check

            CALL deallocatedatetime(mtime_inp_date)
      
          ENDIF ! iostat /= 0
        END DO volcspec
        CLOSE(iunit)
      
        ! specific input for volcano was read
        IF (lnew_volcspec) THEN
          ! an error has occured --> reset values to latest correct ones
          IF (lerror_volcspec) THEN
            volc_data%info(jvolcano)%ash_fraction = frac_save
            volc_data%info(jvolcano)%cfactor(iasha:iashc) = cfac_save(iasha:iashc)
            volc_data%info(jvolcano)%mastin_fac   = mfac_save
            volc_data%info(jvolcano)%mastin_exp   = mexp_save
            WRITE(message_text,'(a,f7.1,a,f7.1,4(a,f8.6),2(a,f5.3))')  &
              & TRIM(volcspec_file)//": "//                            &
              & TRIM(mtime_inp_datetime_str_save)//"  ",               &
              & volc_data%info(jvolcano)%height_above_vent,            &
              & "  elevation=", volc_data%info(jvolcano)%elev,         &
              & "  frac=", frac_save,                                  &
              & "  cfa1=", cfac_save(iasha),                           &
              & "  cfa2=", cfac_save(iashb),                           &
              & "  cfa3=", cfac_save(iashc),                           &
              & "  mfac=", mfac_save,                                  &
              & "  mexp=", mexp_save 
          ELSE
            message_text = volcspec_text
          ENDIF
          CALL message(thisroutine, message_text, all_print=.TRUE.)
        ENDIF ! lnew_volcspec

      ELSE
        calc_erup_rate = .FALSE.
      
        IF (volc_data%lprint_once) THEN   ! ... do ONLY ONCE
          WRITE(message_text,'(a,f10.1,a,f7.1,a,f7.1,4(a,f8.6),2(a,f5.3))')  &
            & "REMARK: NO '"//TRIM(volcspec_file)//"' file -- ",             &
            & volc_data%info(jvolcano)%eruption_rate, "  ",                  &
            & volc_data%info(jvolcano)%height_above_vent,                    &
            & "  frac=", volc_data%info(jvolcano)%ash_fraction,              &
            & "  elevation=", volc_data%info(jvolcano)%elev,                 &
            & "  cfa1=", volc_data%info(jvolcano)%cfactor(iasha),            &
            & "  cfa2=", volc_data%info(jvolcano)%cfactor(iashb),            &
            & "  cfa3=", volc_data%info(jvolcano)%cfactor(iashc),            &
            & "  mfac=", volc_data%info(jvolcano)%mastin_fac,                &
            & "  mexp=", volc_data%info(jvolcano)%mastin_exp
          CALL message(thisroutine, message_text, all_print=.TRUE.)
          volc_data%lprint_once = .FALSE.
          CALL message(thisroutine, "Complete path: "//TRIM(volc_data%volcanofile_path)// &
            &                       TRIM(volcspec_file),all_print=.TRUE.)
        ENDIF ! volc_data%lprint_once
      
      ENDIF
      
      !----------------------------------------------------------------------------------
      ! Calculate emission profile on ICON grid based on height_above_vent,
      !                                                  eruption_rate and ash_fraction
      !----------------------------------------------------------------------------------
      IF (volc_data%info(jvolcano)%height_above_vent > 0.0_wp) THEN ! calculate emission profile
        IF( calc_erup_rate ) THEN  ! calculate eruption rate following Mastin (2009)
          volc_data%info(jvolcano)%eruption_rate =                         &
            &     volc_data%info(jvolcano)%ash_fraction                    &
            &  *  ( volc_data%info(jvolcano)%mastin_fac                    &
            &  *  volc_data%info(jvolcano)%height_above_vent / 1000.0_wp ) &
            &  ** volc_data%info(jvolcano)%mastin_exp
        ENDIF

        IF (.NOT.iostat == 0) THEN
          top_of_plume = volc_data%info(jvolcano)%height_above_vent + volc_data%info(jvolcano)%elev
        ENDIF

        DO jk=1, nlev
! JF:     IF (top_of_plume >= hlev(idx, jk+1) .AND. hlev(idx, jk+1) >= &
! JF:       &                                                  volc_data%info(jvolcano)%elev) THEN 
! JF:     ! use this if emission profil should end at height of volcanic vent
          IF (top_of_plume >= hlev(idx, jk+1)) THEN
            ! replace  hlev(idx, nlev+1) by volc_data%info(jvolcano)%elev if emission profil 
            !   should end at height of volcanic vent
            CALL emiss_vert_normaldistri(top_of_plume, hlev(idx, nlev+1), &
              &  hlev(idx, jk+1),                                         &   ! level bottom
              &  hlev(idx, jk  ), volc_data%height_factor(jk,jvolcano))
          ENDIF
        ENDDO
      ENDIF
    ENDIF 
  ENDDO

END SUBROUTINE art_prepare_emission_volc 
!!
!!--------------------------------------------------------------------------------  
!!
SUBROUTINE art_calculate_emission_volc(jb, dz, cellarea, volc_data, nlev, nmodes, imodes, & 
  &                                    emiss_rate)
!<
! SUBROUTINE art_calculate_emission_volc
! Calculates the emission flux from the emission rate
! Part of Module: mo_art_emission_volc_2mom
! Author: Carolin Walter, KIT
! Initial Release: 2015-05-13
! Modifications:
! - ...
!>
  INTEGER, INTENT(IN)                 :: &  !<
    &  jb,                               &  !< block index
    &  nlev,                             &  !< number of full levels
    &  imodes(:),                        &  !< indices for mass mixing ratios
    &  nmodes                               !< number of modes
  REAL(wp), INTENT(IN)                :: &  !<
    &  dz(:,:)                              !< layer thickness
  REAL(wp), INTENT(IN)                :: &  !<
    &  cellarea(:)                          !< cell area  
  REAL(wp), INTENT(inout)             :: &  !<
    &  emiss_rate(:,:,:)                    !< output emission rate [ug m-3 s-1]
  TYPE(t_art_volcdata), INTENT(inout) :: &
    &  volc_data             !< Volcano data storage type
!Local variables
  INTEGER               :: & !<
    &  jk,                 & !< loop index over vertical levels
    &  imod,               & !< loop index over modes
    &  jvolcano,           & !< loop index over number of volcanoes
    &  jlocal                !< loop index over local points
  INTEGER, POINTER      :: &
    &  idx, blk
  REAL(wp)              :: &
    &  gridvol               !< volume of grid cell

    DO jlocal = 1, volc_data%ithis_nlocal_pts
      jvolcano = volc_data%owner(jlocal)
      idx => volc_data%info(jlocal)%tri_iidx_loc
      blk => volc_data%info(jlocal)%tri_iblk_loc
      
      IF (blk == jb) THEN
        DO jk = 1, nlev
          gridvol = cellarea(idx)*dz(idx,jk)
          DO imod = 1, nmodes
            emiss_rate(idx,jk,imod) = volc_data%info(jvolcano)%eruption_rate         &
              &                     * volc_data%height_factor(jk,jvolcano)           &
              &                     * volc_data%info(jvolcano)%cfactor(imodes(imod)) &
              &                     / gridvol*1.0e9_wp          ! converts output rate [kg m-3 s-1]
          ENDDO ! imod
        ENDDO ! jk
      ENDIF ! blk == jb
    ENDDO ! jlocal
  
END SUBROUTINE art_calculate_emission_volc
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE split_volcstring(zline, n, pos, ilength)
!<
! SUBROUTINE split_volcstring
! Splitting of the string from the input dataset
! Part of Module: mo_art_emission_volc_2mom
! Author: Kristina Lundgren, KIT
! Initial Release: 2012-01-30
! Modifications:
! YYYY-MM-DD: <name>, <institution>
! - ...
!>
  IMPLICIT NONE

  CHARACTER(len=*), INTENT(IN)      :: zline              ! string containing list
  INTEGER,          INTENT(OUT)     :: n                  ! number of parts
  INTEGER,          INTENT(INOUT)   :: pos(:), ilength(:) ! position, lengths of parts
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
      ENDIF
      istart = i+1
      l_word_open = .FALSE.
    ENDIF
  END DO

END SUBROUTINE split_volcstring
!!
!!-------------------------------------------------------------------------
!!
SUBROUTINE emiss_vert_normaldistri(emiss_height,oro,lev_bottom, lev_top,vert_fac)
!<
! SUBROUTINE emiss_vert_normaldistri
! Calculate vertical distribution of volcanic ash plume
! Based on: Rieger et al. (2015): ICON-ART 1.0 - a new online-coupled model system 
!                                 from the global to regional scale
! Part of Module: mo_art_emission_volc_2mom
! Author: Kristina Lundgren, KIT
! Initial Release: 2012-01-30
! Modifications:
! 2015-07-16: Carolin Walter, KIT
! - bugfix: oberserved height = top_of_plume = height_above_vent + elevation
!>
IMPLICIT NONE

  REAL(wp), INTENT(IN) ::  emiss_height  ! in
  REAL(wp), INTENT(IN) ::  oro
  REAL(wp), INTENT(IN) ::  lev_bottom
  REAL(wp), INTENT(IN) ::  lev_top

  REAL(wp), INTENT(INOUT) ::  vert_fac  ! out

  REAL(wp) lev_bottom_star ! normalized heights
  REAL(wp) lev_top_star

  REAL(wp) a1, a2, a3, a4, normfac, sqrtpi

  sqrtpi = 1.772453850905516_wp

  a1 = 0.0076_wp
  a2 = 0.9724_wp
  a3 = 0.4481_wp
  a4 = 0.3078_wp
  
  lev_bottom_star = (lev_bottom - oro) / (emiss_height-oro)
  IF(lev_top > emiss_height) THEN
    lev_top_star    = 1.0_wp
  ELSE
    lev_top_star = (lev_top - oro) / (emiss_height-oro)
  ENDIF

  ! wolframalpha integral solutions
  ! integral_c^d (a1+a2 exp(-((x-a3)/a4)^2)) dx = -a1 c+a1 d+1/2 sqrtpi a2 a4 
  !                                               erf((a3-c)/a4)-erf((a3-d)/a4))
  ! integral_0^1 (a1+a2 exp(-((x-a3)/a4)^2)) dx = a1+1/2 sqrtpi a2 a4
  !                                               (erf((1-a3)/a4)+erf(a3/a4))

  normfac= a1 + 0.5_wp *  sqrtpi *  a2 * a4 * (ERF((1.0_wp-a3)/a4)+ERF(a3/a4))


  vert_fac = ( -a1 *lev_bottom_star + a1 *lev_top_star   &
    &      + 0.5_wp *sqrtpi *a2 *a4 *(ERF((a3-lev_bottom_star)/a4)-ERF((a3-lev_top_star)/a4)))/normfac

END SUBROUTINE emiss_vert_normaldistri
!!
!!-------------------------------------------------------------------------
!!
END MODULE mo_art_emission_volc_2mom
