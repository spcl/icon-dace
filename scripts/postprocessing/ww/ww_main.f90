!+ Calculate WW for ICON
!
PROGRAM ww_icon
!
! Description:
! <Say what this module is for>
!
! Current Code Owner: DWD, <Name of person responsible for this code>
!    <smail, phone, fax and email>
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! @VERSION@    @DATE@     <Your name>
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!
! Modules used:
!
  USE grib_api

#ifdef ONLYWW
  USE mo_wwonly, ONLY: wp, vct_a_from_hhl, kstart_moist
  USE mo_nwp_ww, ONLY: configure_ww, ww_diagnostics, rkgrenz1, rkgrenz2,            &
                           rf_fog, clc_fog,                                         &
                           rain_l_m,  rain_m_s,  snow_l_m, snow_m_s,                &
                           rash_lm_s, rash_s_vs, snsh_l_ms,                         &
                           driz_l_m, driz_m_s, drif_l_ms, raif_l_ms,                &
                           rgdiff_th1, rgdiff_th2                            

#else
  USE mo_kind,   ONLY: wp
  USE mo_nwp_ww, ONLY: configure_ww, ww_diagnostics
#endif


  IMPLICIT NONE

  REAL(wp), ALLOCATABLE:: t  (:,:),  &  ! temperature (K)
                          qv (:,:),  &  ! specific humidity (kg/kg)
                          qc (:,:),  &  ! specific cloud liquid water (kg/kg)
                          u  (:,:),  &  ! zonal comp. of wind velocity      (m/s)
                          v  (:,:),  &  ! meridional comp. of wind velocity (m/s)
                          clc(:,:),  &  ! meridional comp. of wind velocity (m/s)
                          pf (:,:),  &  ! pressure on full levels
                          ph (:,:),  &  ! pressure on half levels (Pa)
                          hhl(:,:)      ! height of half levels (m)

  REAL(wp), ALLOCATABLE:: t_2m (:),  &  ! temperature at 2m above the ground (K)
                          td_2m(:),  &  ! dewpoint temperature at 2m above ground (K)
                          t_g  (:),  &  ! surface temperature (weighted from t_s and t_snow) (K)
                          ps   (:),  &  ! surface pressure
                          clct (:),  &  ! total cloud cover
                          clcm (:),  &  ! medium cloud cover
                          u_10m(:),  &  ! zonal wind component at 10m above the ground (m/s)
                          v_10m(:),  &  ! meridional wind component at 10m above the ground (m/s)
                          rain_gsp (:), & ! grid scale rain accum. until ivv
                          rain_con (:), & ! convective rain accum. until ivv
                          snow_gsp (:), & ! grid scale snow accum. until ivv
                          snow_con (:)    ! convective snow accum. until ivv

  REAL(wp), ALLOCATABLE:: hbas_con(:),  & ! height of base of main convective cloud
                          htop_con(:)     ! height of top of main convective cloud
  INTEGER , ALLOCATABLE:: ibas_con(:),  & ! index of base of main convective cloud
                          itop_con(:),  & ! index of top of main convective cloud
                          iww     (:)     ! WW values

  REAL(wp), ALLOCATABLE:: rain_gsp0(:), & ! grid scale rain accum. until ivv-hinc
                          rain_con0(:), & ! convective rain accum. until ivv-hinc
                          snow_gsp0(:), & ! grid scale snow accum. until ivv-hinc
                          snow_con0(:)    ! convective snow accum. until ivv-hinc
!
  CHARACTER(LEN=4)  :: ymodel  ! NWP model
  CHARACTER(LEN=240):: ycat    ! directory with input GRIB files
  CHARACTER(LEN=10) :: ydate   ! initial date/time
  INTEGER           :: ivv     ! forecast step in hours

  INTEGER           :: ke      ! number of layers

  INTEGER           :: ie      ! number of values per GRIB record. For the
                                           ! first call ie<=0 to indicate that the
                                           ! arrays will be allocated.
  INTEGER           :: igrib   ! GRIB handle
  INTEGER           :: ierr    ! error indicator
  INTEGER           :: verbosity  ! debug level, controls amount of output

  INTEGER           :: hstart, hstop, hinc ! forecast range in hours for which
                                           ! WW is calculated

  REAL(wp)          :: dhour 
  CHARACTER(LEN=256):: ww_output
  INTEGER           :: iunit
  CHARACTER(LEN=11) :: y_namelist = 'NAMELIST_WW'
  INTEGER           :: idate_time(8)

  NAMELIST / ww_namelist/  ymodel, ycat, ydate, hstart, hstop, hinc, ke, verbosity, &
                           ww_output
  NAMELIST / ww_tune_nml/  rkgrenz1, rkgrenz2,                                      &
                           rf_fog, clc_fog,                                         &
                           rain_l_m,  rain_m_s,  snow_l_m, snow_m_s,                &
                           rash_lm_s, rash_s_vs, snsh_l_ms,                         &
                           driz_l_m, driz_m_s, drif_l_ms, raif_l_ms,                &
                           rgdiff_th1, rgdiff_th2                            


  ymodel = 'ICON'
  CALL DATE_AND_TIME( values=idate_time)
  WRITE( ydate, '(I4.4,3I2.2)') idate_time(1:3), 0
  ycat   = '.'
  hstop  = 12
  hinc   = 3
  hstart = -99
  ke     = -1
  verbosity = 0
  ww_output = ''


! Read the namelist
  OPEN( 35, file=y_namelist, form='formatted')

  READ( 35, ww_namelist, IOSTAT=ierr)
  IF ( ierr /= 0) THEN
    PRINT *, 'Error ', ierr,' reading namelist ww_namelist in ww_icon!'
    STOP 'Error in ww_namelist'
  END IF

  WRITE(*,'(/A)') 'Namelist ww_namelist:'
  WRITE(*,'(A13,A)')     ' ymodel    = ', ymodel
  WRITE(*,'(A13,A)')     ' ycat      = ', TRIM(ycat)
  WRITE(*,'(A13,A)')     ' ydate     = ', ydate
  WRITE(*,'(A13,I4)')    ' hstart    = ', hstart
  WRITE(*,'(A13,I4)')    ' hstop     = ', hstop
  WRITE(*,'(A13,I4)')    ' hinc      = ', hinc
  WRITE(*,'(A13,I4)')    ' ke        = ', ke
  WRITE(*,'(A13,I4/)')   ' verbosity = ', verbosity


! Check namelist parameters
  IF ( hstop <= 0) THEN
    PRINT *, 'Error in ww_icon: hstop =', hstop, 'h < 0'
    STOP 'Incorrect hstop'
  END IF
  IF ( hstop < hstart) THEN
    PRINT *, 'Error in ww_icon: hstop =', hstop, 'h < hstart =', hstart
    STOP 'Incorrect hstop'
  END IF
  IF ( hinc <= 0) THEN
    PRINT *, 'Error in ww_icon: hinc =', hinc, ' <= 0'
    STOP 'Incorrect hinc'
  END IF
  IF ( hstart < 0 ) hstart = hinc
  dhour = hinc

  IF ( mod(hstart,hinc) /= 0 .AND. hstart < hinc) THEN
    PRINT *, 'Error in ww_icon: Time intervall hinc = ', hinc,                 &
             'h not possible with first value hstart = ', hstart,' h < hinc'
    STOP 'Incorrect hstart and hinc'
  END IF

  IF ( ke <= 0) THEN
    PRINT *, 'Error in ww_icon: Number of levels ke =', ke, ' <= 0'
    STOP 'Incorrect ke'
  END IF

  IF ( LEN_TRIM(ww_output) == 0) THEN
      ww_output = 'ww_' // TRIM(ymodel) // '.grb'
  END IF
  CALL grib_open_file( iunit, ww_output, 'w', ierr)
  IF ( ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'grib_open_file error ', ierr,' file ', TRIM(ww_output),  &
             ' on unit ', iunit
    STOP 'Error opening WW output file'
  END IF
  IF ( verbosity > 9) PRINT *, 'Openend output GRIB file ', TRIM(ww_output)

! Read HHL. In addition read RAIN_GSP, RAIN_CON, SNOW_GSP, SNOW_CON if hstart > hinc
  ie = 0
  CALL readgribs4ww( 'INIT', ycat, ydate, hstart-hinc, ke, verbosity, ie, &
                     igrib, ierr)

! Allocate rain_gsp0, ..
  IF ( verbosity > 0) PRINT *, 'Allocate rain_gsp0, rain_con0,... with dimension ', ie
  ALLOCATE( ph(ie,ke+1), rain_gsp0(ie), rain_con0(ie), snow_gsp0(ie), snow_con0(ie), &
            ibas_con(ie), itop_con(ie), iww(ie), STAT = ierr)
  IF ( ierr /= 0) THEN
    PRINT *, 'Error ', ierr, ' allocating ph, rain_gsp0, ...'
    STOP  'Error allocating ph, rain_gsp0, ...'
  END IF
  ibas_con (:) = 0
  itop_con (:) = 0

  IF ( hstart > hinc) THEN
!   Check if precipitation values hinc hours before hstart are known
    IF ( rain_gsp(1) < -100._wp) THEN
      PRINT *, 'Error! RAIN_GSP at step ', hstart-hinc, ' was not read!'
      STOP 'Missing rain_gsp0'
    END IF
    IF ( rain_con(1) < -100._wp) THEN
      PRINT *, 'Error! RAIN_CON at step ', hstart-hinc, ' was not read!'
      STOP 'Missing rain_con0'
    END IF
    IF ( snow_gsp(1) < -100._wp) THEN
      PRINT *, 'Error! SNOW_GSP at step ', hstart-hinc, ' was not read!'
      STOP 'Missing snow_gsp0'
    END IF
    IF ( snow_con(1) < -100._wp) THEN
      PRINT *, 'Error! SNOW_CON at step ', hstart-hinc, ' was not read!'
      STOP 'Missing snow_con0'
    END IF
  ELSE
    rain_gsp0(:) = 0._wp
    rain_con0(:) = 0._wp
    snow_gsp0(:) = 0._wp
    snow_con0(:) = 0._wp
  END IF

! Find a reference height profile over sea level
  CALL vct_a_from_hhl( ie, ke+1, hhl, verbosity)
  CALL configure_ww( -1, 1, ke, 0, ymodel)

! Read tuning variables for WW
  READ( 35, ww_tune_nml, IOSTAT=ierr)
  IF ( ierr /= 0) THEN
    PRINT *, 'Error ', ierr,' reading namelist ww_tune_nml in ww_icon!'
    STOP 'Error in ww_namelist'
  END IF
  CLOSE( 35)

  WRITE(*,'(/A)') 'Namelist WW_TUNE_NML:'
  WRITE(*,'(A13,F5.2)') ' rkgrenz1  = ', rkgrenz1
  WRITE(*,'(A13,F5.2)') ' rkgrenz2  = ', rkgrenz2
  WRITE(*,'(A13,F6.3)') ' rf_fog    = ', rf_fog
  WRITE(*,'(A13,F6.3)') ' clc_fog   = ', clc_fog
  WRITE(*,'(A13,F6.3)') ' rain_l_m  = ', rain_l_m
  WRITE(*,'(A13,F6.3)') ' rain_m_s  = ', rain_m_s
  WRITE(*,'(A13,F6.3)') ' snow_l_m  = ', snow_l_m
  WRITE(*,'(A13,F6.3)') ' snow_m_s  = ', snow_m_s
  WRITE(*,'(A13,F6.3)') ' rash_lm_s = ', rash_lm_s
  WRITE(*,'(A13,F6.3)') ' rash_s_vs = ', rash_s_vs
  WRITE(*,'(A13,F6.3)') ' snsh_lm_s = ', rash_lm_s
  WRITE(*,'(A13,F6.3)') ' driz_l_m  = ', driz_l_m
  WRITE(*,'(A13,F6.3)') ' driz_m_s  = ', driz_m_s
  WRITE(*,'(A13,F6.3)') ' drif_l_ms = ', drif_l_ms
  WRITE(*,'(A13,F6.3)') ' raif_l_ms = ', raif_l_ms
  WRITE(*,'(A13,F6.3)') ' rgdiff_th1= ', rgdiff_th1
  WRITE(*,'(A13,F6.3)') ' rgdiff_th2= ', rgdiff_th2

! Read data and calculate WW
  DO ivv = hstart, hstop, hinc

    WRITE(*,'(/3A,i4,A)')  'Model ', ymodel, ': step = ', ivv, 'h'
    CALL readgribs4ww( ymodel, ycat, ydate, ivv, ke, verbosity, ie, &
                       igrib, ierr)
    IF ( ierr /= 0) THEN
      PRINT *, 'Error calling readgribs4ww in ww_icon: ierr =', ierr
      STOP  'Error in readgribs4ww'
    END IF

!   Find indices ibas_con and itop_con from hbas_con and htop_con
    CALL prep4ww( ie, ke, hhl, pf, ph, ps,                   &
                  hbas_con, htop_con, ibas_con, itop_con)

!     Calculate WW
    IF ( verbosity > 2) PRINT *, 'Call ww_diagnostics'
    CALL ww_diagnostics( ie, ke, ke+1, 1, ie, 1,                             &
                         t   , qv   , qc , u   , v   , clc  , pf   , ph,     &
                         t_2m, td_2m, t_g, clct, clcm, u_10m, v_10m,         &
                         rain_gsp0, rain_gsp, rain_con0, rain_con,           &
                         snow_gsp0, snow_gsp, snow_con0, snow_con,           &
                         ibas_con, itop_con, dhour, iww)

    CALL write_ww( ie, iww, iunit, igrib)
    CALL print_stat( 6, ivv, iww, ie, 'WW')

!   Save old values of rain_gsp, ... for next call of ww_diagnostics
    rain_gsp0(:) = rain_gsp(:)
    rain_con0(:) = rain_con(:)
    snow_gsp0(:) = snow_gsp(:)
    snow_con0(:) = snow_con(:)
    
  END DO

  IF ( verbosity > 9) PRINT *, 'Close output GRIB file'
  CALL grib_close_file( iunit)

CONTAINS

  SUBROUTINE write_ww( ie, iww, iunit, igrib)
!
!   Write GRIB records of WW to file
!
    INTEGER, INTENT(IN) :: ie, iunit, igrib
    INTEGER, INTENT(IN) :: iww(ie)

    INTEGER :: iedi

    CALL grib_set( igrib, 'shortName', 'WW')
    CALL grib_set( igrib, 'typeOfLevel', 'surface')

!   set level to surface and height to 0
    CALL grib_set( igrib, 'values', iww)
    CALL grib_get( igrib, 'edition', iedi)
    CALL DATE_AND_TIME( values=idate_time)
    IF ( iedi == 1) THEN
      CALL grib_set( igrib, 'topLevel',    0)
      CALL grib_set( igrib, 'bottomLevel', 0)
      CALL grib_set( igrib,'localDecodeDateYear'  , idate_time(1)-1900)
      CALL grib_set( igrib,'localDecodeDateMonth' , idate_time(2) )
      CALL grib_set( igrib,'localDecodeDateDay'   , idate_time(3) )
      CALL grib_set( igrib,'localDecodeDateHour'  , idate_time(5) )
      CALL grib_set( igrib,'localDecodeDateMinute', idate_time(6) )
    ELSE
      CALL grib_set_missing( igrib, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing( igrib, 'scaledValueOfFirstFixedSurface')
      CALL grib_set_missing( igrib, 'scaleFactorOfFirstFixedSurface')
      CALL grib_set_missing( igrib, 'scaledValueOfFirstFixedSurface')
      CALL grib_set( igrib,'localCreationDateYear' , idate_time(1) )
      CALL grib_set( igrib,'localCreationDateMonth', idate_time(2) )
      CALL grib_set( igrib,'localCreationDateDay'  , idate_time(3) )
      CALL grib_set( igrib,'localCreationDateHour',  idate_time(5) )
      CALL grib_set( igrib,'localCreationDateMinute',idate_time(6) )
      CALL grib_set( igrib,'localCreationDateSecond',idate_time(7) )
    END IF

    IF ( verbosity > 0) THEN
      PRINT '(A,i2,a,I4.4,2(A1,I2.2),I3,2(A1,I2.2)/)', 'Write WW: edition',     &
            iedi, ' at ', idate_time(1),'-', idate_time(2),'-', idate_time(3), &
            idate_time(5),':', idate_time(6), ':', idate_time(7)
    END IF
    CALL grib_write( igrib, iunit)
    CALL grib_release( igrib)
  END SUBROUTINE write_ww


  SUBROUTINE prep4ww ( ie, ke, hhl, pf, ph, ps,                 &
                       hbas_con, htop_con, bas_con , top_con)
!
!  Prepare for ww_diagnostics, i.e. calculate pf from ph, bas_con and
!  top_con from hbas_con, htop_con and hhl.
!
    INTEGER,  INTENT(IN) :: ie, ke
    REAL(wp), INTENT(IN) :: hhl(ie,ke+1), ps(ie), hbas_con(ie), htop_con(ie)

    REAL(wp), INTENT(INOUT):: pf(ie,ke), ph(ie,ke+1)
    INTEGER , INTENT(OUT):: bas_con(ie), top_con(ie)

    INTEGER :: i, k
    REAL(wp):: hh
    REAL(wp), PARAMETER :: delta_h = 0.5_wp 

  
    IF ( verbosity > 2) PRINT *, 'Calculate ibas_con, itop_con'
    bas_con(:) = 0
    top_con(:) = 0
    DO k = 1, ke
      DO i = 1, ie
        hh = 0.5*( hhl(i,k) + hhl(i,k+1) )
        IF ( ABS( hbas_con(i) - hhl(i,k)) < 0.5_wp) bas_con(i) = k
        IF ( ABS( htop_con(i) - hhl(i,k)) < 0.5_wp) top_con(i) = k
      END DO
    END DO
    IF ( verbosity > 2) THEN
      PRINT *, 'Statistics of BAS_CON'
      CALL print_stat( 6, ivv, bas_con, ie, 'BAS_CON')
      PRINT *, 'Statistics of TOP_CON'
      CALL print_stat( 6, ivv, top_con, ie, 'TOP_CON')
    END IF

    IF ( ymodel(1:2) == 'IC' .OR. ymodel(1:2) == 'LM' .OR. ymodel == 'TEST' ) THEN
!     For test calculate ph from pf only approximately
      IF ( verbosity > 2) PRINT *, 'Calculate ph from pf and ps'
      ph(:,1) = 0.5_wp*pf(:,1)   ! top level not important for WW
      DO k = 2, ke
        DO i = 2, ie
          ph(i,k) = SQRT( pf(i,k)*pf(i,k-1) )
        END DO
      END DO
      ph(:,ke+1) = ps(:)
    ELSE
      IF ( verbosity > 2) PRINT *, 'Calculate pf from ph'
      DO k = 1, ke
        DO i = 1, ie
          pf(i,k) = 0.5_wp*( ph(i,k) + ph(i,k+1) )
        END DO
      END DO
    END IF
  END SUBROUTINE prep4ww


  SUBROUTINE readgribs4ww( ymodel, ycat, ydate, ivv, ke, iverb, ie, &
                           igrib_t_2m, ierr)
!
!  Read fields needed to calculate WW.
!  For atmospheric fields (2d) only levels lev >= kstart_moist
!  are required.
!

  CHARACTER(LEN=*),  INTENT(IN) :: ymodel  ! NWP model
                                           ! ymodel = 'HHL' to read only hhl
                                           ! rain_gsp, rain_con, snow_gsp, snow_con
  CHARACTER(LEN=*),  INTENT(IN) :: ycat    ! directory with input GRIB files
  CHARACTER(LEN=10), INTENT(IN) :: ydate   ! initial date/time
  INTEGER, INTENT(IN)           :: ivv     ! forecast step in hours

  INTEGER, INTENT(IN)           :: ke      ! number of layers
  INTEGER, INTENT(IN)           :: iverb   ! verbosity level

  INTEGER, INTENT(INOUT)        :: ie      ! number of values per GRIB record. For the
                                           ! first call ie<=0 to indicate that the
                                           ! arrays will be allocated.
  INTEGER, INTENT(OUT)          :: igrib_t_2m ! GRIB handle for T_2M. This is used to
                                              ! clone a GRIB handle for WW
  INTEGER, INTENT(OUT)          :: ierr       ! error indicator

! Local arrays:
  REAL(wp), ALLOCATABLE :: values(:)

  INTEGER, PARAMETER    :: t_idx        =  1, &
                           qv_idx       =  2, &
                           qc_idx       =  3, &
                           u_idx        =  4, &
                           v_idx        =  5, &
                           c_idx        =  6, &
                           p_idx        =  7
!                          hhl_idx      =  8
  INTEGER, PARAMETER ::    ps_idx       =  1, &
                           t_2m_idx     =  2, &
                           td_2m_idx    =  3, &
                           t_g_idx      =  4, &
                           clct_idx     =  5, &
                           clcm_idx     =  6, &
                           u_10m_idx    =  7, &
                           v_10m_idx    =  8, &
                           rain_gsp_idx =  9, &
                           rain_con_idx = 10, &
                           snow_gsp_idx = 11, &
                           snow_con_idx = 12, &
                           hbas_con_idx = 13, &
                           htop_con_idx = 14   

   CHARACTER(LEN=8), PARAMETER :: varname2d(p_idx) = (/   &
  &                  'T       ', 'QV      ', 'QC      ', 'U       ', 'V       ', \
  &                  'CLC     ', 'P       ' /)
   CHARACTER(LEN=8), PARAMETER :: varname1d(1:htop_con_idx) = (/   &
  &                  'PS      ', 'T_2M    ', 'TD_2M   ', 'T_G     ', \
  &                  'CLCT    ', 'CLCM    ', 'U_10M   ', 'V_10M   ', 'RAIN_GSP', \
  &                  'RAIN_CON', 'SNOW_GSP', 'SNOW_CON', 'BAS_CON ', 'TOP_CON ' /)

  LOGICAL, ALLOCATABLE  :: lvar2d(:,:), lvar1d(:), l_hhl(:)

  INTEGER            :: iunit, igrib, iedi, mdim, lev, istep
  CHARACTER(LEN=13)  :: ww_input      ! file name for input data
  CHARACTER(LEN=256) :: in_file
  CHARACTER(LEN=20)  :: shortname

  INTEGER            :: idate, itime, i
  CHARACTER(LEN=10)  :: ydatetime, yform

  IF ( iverb > 1) PRINT *, 'Start readgribs4ww: ', ymodel, ydate, ivv, ie
  ierr = 0

  IF ( ymodel == 'INIT') THEN
    ie = 0
    in_file =  TRIM(ycat) // '/ww_input_hhl'
    ALLOCATE( l_hhl(ke+1), STAT=ierr)
    IF ( ierr /= 0) THEN
      PRINT *, 'Error ', ierr,  &
               ' allocating arrays lvar1d, l_hhl in readgribs4ww. ', ke+1
      RETURN
    END IF
    l_hhl ( :)   = .false.
  ELSE
    WRITE( ww_input, '(A9,I4.4)') 'ww_input_', ivv
    in_file =  TRIM(ycat) // '/' // ww_input
  END IF
  IF ( iverb > 8) PRINT *, 'Open input file ', TRIM(in_file)
  CALL grib_open_file( iunit, TRIM(in_file), 'r', ierr)
  IF ( ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'grib_open_file error ', ierr,' file ', TRIM(in_file),  &
             ' on unit ', iunit
    ierr = 1
    RETURN
  ELSE IF ( iverb > 8) THEN
    PRINT *, 'Opened input file ', TRIM(in_file), ' iunit=',iunit
  END IF

  IF ( ie > 0) THEN
    ALLOCATE( values(ie), lvar2d( 1:ke, p_idx), lvar1d( ps_idx:htop_con_idx), &
              STAT=ierr)
    IF ( ierr /= 0) THEN
      PRINT *, 'Error ', ierr,  &
               ' allocating arrays values, lvar1d, lvar2d in readgribs4ww. ie, ke+1', ie, ke+1
      ierr = 8
      RETURN
    END IF
    lvar1d(:)   = .false.
    lvar2d(:,:) = .false.
  END IF

  DO
    CALL grib_new_from_file( iunit, igrib, ierr)
    IF ( ierr == GRIB_END_OF_FILE) EXIT
    IF ( ierr /= GRIB_SUCCESS) THEN
      PRINT *, 'grib_new_from_file error ', ierr
      CYCLE
    END IF

!   Check date, time
    CALL grib_get( igrib, 'dataDate', idate)
    CALL grib_get( igrib, 'dataTime', itime)
    WRITE( ydatetime,'(i8.8,i2.2)') idate, itime/100
    IF ( ydatetime /= ydate) THEN
      PRINT *, 'Wrong date/time ', ydatetime, ' /= ', ydate,'! Read next field.'
      CYCLE
    END IF

    CALL grib_get( igrib, 'shortName', shortname )

!   Check forecast hour
    CALL grib_set( igrib, 'stepUnits', 'h')
    CALL grib_get( igrib, 'step', istep, ierr)
    IF ( ierr /= GRIB_SUCCESS) THEN
      PRINT *, 'Error ', ierr, ' reading step'
      CYCLE
    END IF
    IF ( istep /= ivv .AND. shortName /= 'HHL' ) THEN
      PRINT *, shortname, ': Wrong step ', istep, 'h /= ', ivv,'h! Read next field.'
      CYCLE
    END IF

    CALL grib_get_size( igrib, 'values', mdim)

    IF ( ie <= 0) THEN
!     First call of readgribs4ww:
!     Allocate the arrays, read HHL
      ie = mdim
      PRINT *, 'Horizontal dimension ie =', ie
      ALLOCATE( values(ie), STAT=ierr )
      IF ( ierr /= 0) THEN
        PRINT *, 'Error ', ierr, ' allocating array values in readgribs4ww. ie=', ie
        ierr = 9
        RETURN
      ELSE IF ( iverb > 2) THEN
        PRINT *, 'Allocated array values with dimension ', ie
      END IF

      ALLOCATE( hhl(ie,ke+1), STAT=ierr)
      IF ( ierr /= 0) THEN
        PRINT *, 'Error ', ierr, ' allocating HHL in readgribs4ww'
        PRINT *, 'Dimensions: ', ie, ke
        ierr = 10
        RETURN
      ELSE IF ( iverb > 2) THEN
        PRINT *, 'Allocated hhl with dimension ', ie, ke+1
      END IF

      ALLOCATE( t(ie,ke), qv(ie,ke), qc(ie,ke), u(ie,ke), v(ie,ke), clc(ie,ke),     &
                pf(ie,ke), STAT=ierr)
      IF ( ierr /= 0) THEN
        PRINT *, 'Error ', ierr, ' allocating t, qv, qc, u, v, clc, pf in readgribs4ww'
        PRINT *, 'Dimensions: ', ie, ke
        ierr = 10
        RETURN
      END IF

      ALLOCATE( t_2m(ie), td_2m(ie), t_g(ie), ps(ie),                   &
                clct(ie), clcm(ie), u_10m(ie), v_10m(ie),               &
                rain_gsp(ie), rain_con(ie), snow_gsp(ie), snow_con(ie), &
                hbas_con(ie), htop_con(ie), STAT=ierr)
      IF ( ierr /= 0) THEN
        PRINT *, 'Error ', ierr, ' allocating t_2m, ..., htop_con in readgribs4ww'
        PRINT *, 'Dimensions: ', ie, ke
        ierr = 10
        RETURN
      END IF

      ALLOCATE( lvar2d( 1:ke,htop_con_idx), lvar1d( ps_idx:htop_con_idx), STAT=ierr)
      IF ( ierr /= 0) THEN
        PRINT *, 'Error ', ierr, ' allocating array lvar2d, lvar1d in readgribs4ww'
        ierr = 11
        RETURN
      END IF
      lvar2d(:,:) = .false.
      lvar1d(:)   = .false.

    ELSE IF ( mdim /= ie) THEN

      PRINT *, 'ERROR in readgribs4ww! ', TRIM(shortname), ' with ', mdim, ' values!'
      PRINT *, 'Expected ', ie, ' values.'
      CYCLE

    END IF

    CALL grib_get( igrib, 'values', values)
    CALL grib_get( igrib, 'edition', iedi)
    IF ( iedi == 1) THEN
        CALL grib_get( igrib, 'topLevel', lev)
    ELSE
        CALL grib_get( igrib, 'scaledValueOfFirstFixedSurface', lev)
    END IF
    IF ( iverb > 0) THEN
      WRITE(*,'(3A,I8.8,A,I4.4,A,i2,A,I5)') 'Read ', shortname, ' at ', &
               idate, ' ', itime, ', edition', iedi,', level ', lev
    END IF

!   Ignore the level of CLCM
    IF ( shortname == 'CLCM') lev = 0

    IF ( lev > ke+1 .OR. &
   &     (shortname /= 'HHL' .AND. shortname /= 'P' .AND. lev > ke) ) THEN
      IF ( shortname == 'HHL' .OR. shortname == 'P') THEN
        PRINT *, 'Warning: level ', lev, ' > ke+1 = ', ke+1
      ELSE
        PRINT *, 'Warning: level ', lev, ' > ke = ', ke
      END IF
      CYCLE
    END IF

    SELECT CASE (TRIM(shortname))
      CASE ('T')
        t  (:,lev) = values(:)
        lvar2d(lev, t_idx) = .true.
      CASE ('QV')
        qv (:,lev) = values(:)
        lvar2d(lev,qv_idx) = .true.
      CASE ('QC')
        qc (:,lev) = values(:)
        lvar2d(lev,qc_idx) = .true.
      CASE ('U')
        u  (:,lev) = values(:)
        lvar2d(lev, u_idx) = .true.
      CASE ('V')
        v  (:,lev) = values(:)
        lvar2d(lev, v_idx) = .true.
      CASE ('CLC')
        clc(:,lev) = values(:)
        lvar2d(lev, c_idx) = .true.
      CASE ('P')
        pf (:,lev) = values(:)
        lvar2d(lev, p_idx) = .true.
      CASE ('HHL')
        hhl(:,lev)= values(:)
        l_hhl(lev) = .true.
      CASE ('T_2M')
        t_2m(:) = values(:)
        lvar1d( t_2m_idx) = .true.
        CALL grib_clone( igrib, igrib_t_2m)  ! used to clone GRIB handle for WW.
      CASE ('TD_2M')
        td_2m (:) = values(:)
        lvar1d( td_2m_idx) = .true.
      CASE ('T_G')
        t_g (:) = values(:)
        lvar1d( t_g_idx) = .true.
      CASE ('PS')
        ps (:) = values(:)
        lvar1d( ps_idx) = .true.
!     CLCT, and CLCM must be rescaled to values from 0 to 1.
      CASE ('CLCT')
        clct (:) = values(:)*0.01
        lvar1d( clct_idx) = .true.
      CASE ('CLCM')
        clcm (:) = values(:)*0.01
        lvar1d( clcm_idx) = .true.
      CASE ('U_10M')
        u_10m (:) = values(:)
        lvar1d( u_10m_idx) = .true.
      CASE ('V_10M')
        v_10m (:) = values(:)
        lvar1d( v_10m_idx) = .true.
      CASE ('RAIN_GSP')
        rain_gsp (:) = values(:)
        lvar1d( rain_gsp_idx) = .true.
      CASE ('RAIN_CON')
        rain_con (:) = values(:)
        lvar1d( rain_con_idx) = .true.
      CASE ('SNOW_GSP')
        snow_gsp (:) = values(:)
        lvar1d( snow_gsp_idx) = .true.
      CASE ('SNOW_CON')
        snow_con (:) = values(:)
        lvar1d( snow_con_idx) = .true.
      CASE ('HBAS_CON')
        hbas_con (:) = values(:)
        lvar1d( hbas_con_idx) = .true.
      CASE ('HTOP_CON')
        htop_con (:) = values(:)
        lvar1d( htop_con_idx) = .true.
      CASE DEFAULT
        PRINT *, 'Field ', TRIM(shortname), ' not needed!'
    END SELECT

    CALL grib_release(igrib)
  END DO


  IF ( iverb > 8) PRINT *, 'Close input file ', TRIM(in_file), ' Unit ', iunit
  CALL grib_close_file( iunit, ierr)
  IF ( ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'grib_close_file error ', ierr,' on unit ', iunit
    ierr = 2
    RETURN
  END IF

!  Check if all fiels were read
    IF ( ymodel == 'INIT') THEN
      IF ( .NOT. ALL( l_hhl( MAX(ke-30,1):ke+1)) ) THEN
        PRINT *, 'Error! Cannot determine kstart_moist. Only the following HHL levels were read:'
        DO i = 1, ke+1
          IF ( l_hhl(i) ) WRITE(*,'(i4)',ADVANCE='NO') i
        END DO
        PRINT *
        PRINT *, 'Expected at least levels ', MAX(ke-30,1), ' to ', ke+1
        STOP 'Read not enough levels of HHL!'
      END IF
      IF ( .NOT. lvar1d( rain_gsp_idx)) rain_gsp(:) = -999._wp
      IF ( .NOT. lvar1d( rain_con_idx)) rain_con(:) = -999._wp
      IF ( .NOT. lvar1d( snow_gsp_idx)) snow_gsp(:) = -999._wp
      IF ( .NOT. lvar1d( snow_con_idx)) snow_con(:) = -999._wp

    ELSE
!!!   IF ( .NOT. ALL(lvar2d(kstart_moist(1):ke+1,:)) .AND. &
      IF ( .NOT. ALL(lvar2d(kstart_moist(1):ke,:)) .AND. &
           .NOT. ALL(lvar1d(:)) ) THEN
        PRINT *, 'Error in readgribs4ww! Not all necessary data found'
        WRITE( yform,'(A4,i3,a3)') '(A8,', ke+2-kstart_moist(1),'i3)'
!!      WRITE(*,yform) 'Var\lev ', ( i, i = kstart_moist(1), ke+1)
!!      WRITE( yform,'(A4,i3,a3)') '(A8,', ke+2-kstart_moist(1),'L3)'
        WRITE(*,yform) 'Var\lev ', ( i, i = kstart_moist(1), ke)
        WRITE( yform,'(A4,i3,a3)') '(A8,', ke+1-kstart_moist(1),'L3)'
        DO i = 1, p_idx
!!        WRITE(*,yform) varname2d(i), lvar2d(kstart_moist(1):ke+1,i)
          WRITE(*,yform) varname2d(i), lvar2d(kstart_moist(1):ke,i)
        END DO
        DO i = ps_idx, htop_con_idx
          WRITE(*,yform) varname1d(i), lvar1d(i)
        END DO
        ierr = 20
        RETURN
      END IF
    END IF

    DEALLOCATE( values, STAT=ierr )
    IF ( ierr /= 0) THEN
      PRINT *, 'Error ', ierr,' deallocating array values!'
    END IF
    DEALLOCATE( lvar2d, lvar1d, STAT=ierr )
    IF ( ierr /= 0) THEN
      PRINT *, 'Error ', ierr,' deallocating array lvar2d, lvar1d!'
    END IF

  END SUBROUTINE readgribs4ww

  SUBROUTINE print_stat( iun, ivv, ivar, ldim, yname)
!
! doing some statistics of calculated ww 
!
! Current Code Owner: DWD, Ulrich Pflueger
! phone: 069 8062 2753
! email: ulrich.pflueger@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V1_0         2014/02/27 Ulrich Pflueger
!  Initial release
!
! Code Description:
! Language: Fortran 90.
!=======================================================================
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iun,ivv,ldim
  INTEGER, INTENT(IN) :: ivar(ldim)
  CHARACTER(LEN=*), INTENT(IN) :: yname

  INTEGER :: i
  INTEGER :: limits1(12) = (/ -9, 0,  1,  2,  3,  4, 40, 50, 60, 70, 80, 90 /)
! INTEGER :: limits2(24) = (/ 0,  1,  2,  3, 45, 48, 50, 56, 60, 63, 65,  &
!                        &   66, 67, 70, 73, 75, 80, 81, 82, 85, 86, 95,  &
!                        &   96, 99 /)
  INTEGER :: limits2(32) = (/ -9, 0,  1,  2,  3, 45, 48, 50, 51, 53, 55, 56, 57, &
                         &   60, 61, 63, 65, 66, 67, 70, 71, 73, 75, 77, 80, &
                         &   81, 82, 85, 86, 95, 96, 99 /)
  REAL(wp) :: rldim

  rldim = 1._wp/REAL(ldim,KIND=wp)
! Control output
  WRITE(iun,'(a)') '*************************************'
  WRITE(iun,'(a)') 'Control Output for a quick check'
  WRITE(iun,*) '  '
  WRITE(iun,'(a,i6)'  ) 'statistics for forecast hour: ',ivv
  WRITE(iun,'(2a,f9.2)') 'mean over all ',yname        , sum(ivar(:))*rldim
  DO i = 1, 12
    WRITE(iun,'(3a,i2,a,i9)') 'number of all ',yname, '(:) >=',limits1(i),  &
   &                            '  :', count(ivar(:) >= limits1(i))
  END DO
  WRITE(iun,'(a)') '-------------------------------------'
  WRITE(iun,'(3a,i9)')'   max of all ', yname, ' :' , maxval(ivar(:))
  WRITE(iun,'(a,i9)') 'location of maximum        :', maxloc(ivar(:))
  IF ( yname == 'TOP_CON' .OR. yname == 'BAS_CON') THEN
  WRITE(iun,'(3a,i9)')'   min of all ', yname, ' > 0:', MINVAL( ivar(:), MASK=( ivar(:)>0))
  END IF
  WRITE(iun,'(a)') '*************************************'

  WRITE(iun,*) ' '
  rldim = 100._wp/REAL(ldim,KIND=wp)
  IF (yname =='WW') THEN
    DO i = 1, 32
      WRITE(iun,'(3a,i2,a,i9,F8.2,a2)') 'number of all ',yname,'(:) = ',limits2(i),      &
   &    '  :', count( ivar(:)==limits2(i)), count( ivar(:)==limits2(i))*rldim,' %'
    END DO
  ELSE
    DO i = 1, 100
      IF (count( ivar(:)==i) > 0) THEN
        WRITE(iun,'(3a,i2,a,i9,F8.2,a2)') 'number of all ',yname,'(:) = ',i,      &
   &      '  :', count( ivar(:)==i), count( ivar(:)==i)*rldim,' %'
      ENDIF  
    END DO
  ENDIF

  END SUBROUTINE print_stat

END PROGRAM ww_icon
