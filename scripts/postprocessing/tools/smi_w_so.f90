PROGRAM smi_w_so
!-------------------------------------------------------------------
!
!  Calculate soil moisture index, SMI, from W_SO for the soil types
!  of GME or COSMO. Or, calculate W_SO from SMI. It is assumed that
!  SMI has the same shortName 'W_SO'.
!  Set W_SO, W_SO_ICE, SMI for ice, rock, sea water, sea ice to 0.
!
!  SOILTYP must be read before W_SO.
!  Other GRIB records than W_SO, W_SO_ICE are just copied to the
!  output file.
!
!    Helmut Frank, DWD, 25.10.2013
!
!-------------------------------------------------------------------
  USE grib_api
  IMPLICIT NONE
  INTEGER, PARAMETER :: ireals = KIND(1.D0)

  CHARACTER(255)    :: gribfile_in, gribfile_out
  CHARACTER(20)     :: which_way
  CHARACTER(32)     :: gridTypeST
  CHARACTER(80)     :: uuidST

  INTEGER, PARAMETER :: nsoil = 10  ! number of soil types
!
!    element    soil type
!     1         ice
!     2         rock
!     3         sand
!     4         sandy loam
!     5         loam
!     6         clay loam 
!     7         clay 
!     8         peat       
!     9         sea water  
!    10         sea ice    
!
  CHARACTER :: ysoiltyp(nsoil)*10 =                                                    &
             (/ 'ice       ', 'rock      ', 'sand      ', 'sandyLoam ', 'loam      ',  &
                'clay loam ', 'clay      ', 'peat      ', 'sea water ', 'sea ice   '/)

  REAL(KIND=ireals) :: Cfcap(nsoil), Cporv(nsoil), Cpwp(nsoil), Cadp(nsoil), &
                       w_so_mean(nsoil)
  INTEGER           :: n_mean(nsoil), n_miss(nsoil)

  INTEGER, ALLOCATABLE           :: soiltyp(:)
  REAL(KIND=ireals), ALLOCATABLE :: zgrib(:)

  INTEGER           :: igrib, ierr
  REAL(KIND=ireals) :: d, wmiss              ! thickness of soil layer, missing value
  INTEGER           :: z1, z2, ifac1, ifac2, oc
  INTEGER           :: nn, np, iednr, nrec, ibitmap, i, isoil
  CHARACTER(20)     :: yname
  INTEGER           :: nu_in = 10, nu_out=20
  LOGICAL           :: ldebug, l_porous, l_soiltyp = .false.

! Soil properties
       Cporv = (/  & ! pore volume    (fraction of volume)
       1.E-10_ireals,1.E-10_ireals, .364_ireals, .445_ireals, .455_ireals,  &
        .475_ireals, .507_ireals, .863_ireals,1.E-10_ireals,1.E-10_ireals/)


       Cfcap(1:nsoil) = (/  & ! field capacity (fraction of volume)
       1.E-10_ireals, 1.E-10_ireals, .196_ireals, 0.260_ireals , 0.340_ireals ,  &
       0.370_ireals ,  0.463_ireals, .763_ireals, 1.E-10_ireals, 1.E-10_ireals /)

       Cpwp =  (/  & ! plant wilting point  (fraction of volume)
        0.0_ireals  , 0.0_ireals  , .042_ireals, .100_ireals, .110_ireals,  &
        .185_ireals, .257_ireals, .265_ireals, 0.0_ireals  , 0.0_ireals  /)

       Cadp =  (/  & ! air dryness point  (fraction of volume)
        0.0_ireals  , 0.0_ireals  , .012_ireals, .030_ireals, .035_ireals,  &
       .060_ireals, .065_ireals, .098_ireals, 0.0_ireals  , 0.0_ireals  /)

  PRINT *, '   Convert W_SO to SMI or SMI to W_SO'
  PRINT *, '   =================================='
  PRINT *, 'Enter input GRIB file. It must contain SOILTYP and W_SO or SMI'
  READ (*,'(A)') gribfile_in
  ierr = INDEX( gribfile_in, ' -D')
  IF ( ierr > 0) THEN
    ldebug = .true.
    gribfile_in(ierr:ierr+2) = '   '
  END IF
  WRITE(*,'(2A)') 'Input file: ', TRIM(gribfile_in)
  PRINT *, 'Which way to convert? Enter one of the following words.'
  PRINT *, ' w_so2smi: From W_SO to SMI'
  PRINT *, ' smi2w_so: From SMI to W_SO'
  PRINT *, ' smi2w_solimit: SMI to W_SO. Limit W_SO by air dryness point and pore volume.'
  PRINT *, ' smi2w_sodef: SMI to W_SO. Limit W_SO by air dryness point and pore volume,'
  PRINT *, '              set missing values over land to mean values of the same soil type.'
  READ (*,'(A)') which_way
  IF ( INDEX(which_way,'P') + INDEX(which_way,'p') > 0 ) THEN
    l_porous = .true.
  ELSE
    l_porous = .false.
  END IF
  IF ( INDEX(which_way, '2S') + INDEX(which_way,'2s') > 0) THEN
    which_way = 'w_so2smi'
  ELSE IF ( INDEX(which_way, 'L') + INDEX(which_way,'l') > 0) THEN
    which_way = 'smi2w_solimit'
  ELSE IF ( INDEX(which_way, 'D') + INDEX(which_way,'d') > 0) THEN
    which_way = 'smi2w_sodef'
  ELSE
    which_way = 'smi2w_so'
  END IF
  PRINT *, 'Conversion ', which_way
  PRINT *, 'Enter output GRIB file.'
  READ (*,'(A)') gribfile_out
  WRITE(*,'(2A)') 'Output file: ', TRIM(gribfile_out)
!
!     Open GRIB files
!
  CALL grib_open_file( nu_in, gribfile_in, 'r', ierr)
  IF (ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'ERROR ',ierr,' opening GRIB file ', TRIM(gribfile_in),' for reading'
    STOP 'ERROR opening GRIB file for reading'
  END IF
  CALL grib_open_file( nu_out,gribfile_out, 'w', ierr)
  IF (ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'ERROR ',ierr,' opening GRIB file ', TRIM(gribfile_out),' for writing'
    STOP 'ERROR opening GRIB file for writing'
  END IF

  nrec = 0
  DO
    nrec = nrec + 1
!
!   Read the GRIB data
!
    CALL grib_new_from_file( nu_in, igrib, ierr)
    IF ( ierr == GRIB_END_OF_FILE) EXIT

    CALL grib_get( igrib, 'numberOfDataPoints', np)
    PRINT *, 'numberOfDataPoints =', np

    IF ( nrec == 1 ) THEN
!
!     Get dimensions of arrays
!
      nn = np
      ALLOCATE( zgrib(nn), soiltyp(nn), STAT=ierr)
      IF ( ierr/=0) THEN
        PRINT *, 'Error allocating zgrib or soiltyp with dimension', nn
        STOP 'Allocation error for zgribin'
      END IF

    ELSE

        IF ( np /= nn) THEN
          PRINT *, 'Error! ', np, ' numberOfDataPoints different from allocated array length ', nn
          STOP 'Different numberOfDataPoints'
        END IF

    ENDIF
!
!   Find the quantities
!
    CALL grib_get ( igrib, 'shortName', yname)

    WRITE(*,'(I3,A10,A12)',ADVANCE='YES') nrec,'. record: ', TRIM(yname)
    CALL grib_get( igrib, 'values', zgrib)

    SELECT CASE ( yname)
      CASE ( 'SOILTYP')
        soiltyp (:) = NINT( zgrib(:) )
        l_soiltyp = .TRUE.
        CALL grib_get( igrib, 'gridType', gridTypeST)
        IF ( gridTypeST == 'unstructured_grid') THEN
          CALL grib_get( igrib, 'uuidOfHGrid', uuidST)
        END IF
        IF( l_porous ) WHERE ( zgrib(:) < 2.5_ireals ) zgrib(:) = zgrib(:) + 10._ireals

      CASE ( 'W_SO', 'SMI')
        IF ( .NOT. l_soiltyp) THEN
          PRINT *, 'Error in smi_w_so: SOILTYP not known, yet! Cannot convert W_SO!'
          STOP 'Error! SOILTYP not known!'
        END IF
        CALL grib_get( igrib, 'edition', iednr)
        IF ( ldebug) PRINT *, '   GRIB edition ', iednr
        IF ( iednr == 2) THEN
!         depth below land in m
          CALL grib_get (igrib, 'scaledValueOfFirstFixedSurface',  z1)
          CALL grib_get (igrib, 'scaledValueOfSecondFixedSurface', z2)
          CALL grib_get (igrib, 'scaleFactorOfFirstFixedSurface',  ifac1)
          CALL grib_get (igrib, 'scaleFactorOfSecondFixedSurface', ifac2)
          IF ( ldebug) PRINT *, '   scaledValueOfFirstFixedSurface  ', z1,    &
                                '   scaledValueOfSecondFixedSurface ', z2,    &
                                '   scaleFactorOfFirstFixedSurface  ', ifac1, &
                                '   scaleFactorOfSecondFixedSurface ', ifac2
!         Thickness of the layer in mm
          d = ( z2*0.1_ireals**ifac2 - z1*0.1_ireals**ifac1) * 1000._ireals
          CALL grib_get (igrib, 'centre',  oc)
        ELSE  ! GRIB1
!         depth below land in cm
          CALL grib_get( igrib, 'level', z1)
          IF ( ldebug) PRINT *, '   level ', z1
          d = z1*10._ireals
          oc = 0
        END IF
        PRINT *, '     layer thickness ', d, ' mm'
        IF ( d <= 0._ireals) STOP 'Error! Layer thickness is 0.'

!!      CALL check_grid( igrib, gridTypeST, uuidST)
        CALL check_grid( igrib)

!       Replace missing values by 0.
        CALL grib_get ( igrib, 'bitmapPresent', ibitmap)
        IF ( ibitmap > 0) THEN
          CALL grib_get ( igrib, 'missingValue', wmiss)
          IF ( ldebug) PRINT *, '   missingValue ', wmiss
          IF ( which_way /= 'smi2w_sodef') WHERE( zgrib(:) == wmiss) zgrib(:) = 0._ireals
          CALL grib_set ( igrib, 'bitmapPresent', 0)
        ELSE
          wmiss = -1.e9   ! set wmiss in case which_way is 'smi2w_sodef'
        END IF

        IF ( yname == 'W_SO' .AND. which_way == 'w_so2smi') THEN
          WHERE ( soiltyp(:) > 2 .AND. soiltyp(:) < 9 )
            zgrib(:) = ( zgrib(:)/d - Cpwp(soiltyp(:)) )/( Cfcap(soiltyp(:))-Cpwp(soiltyp(:)) )
          ELSE WHERE
            zgrib(:) = 0._ireals
          END WHERE
          IF ( oc == 78) CALL grib_set( igrib, 'shortName', 'SMI')

        ELSE IF ( yname /= 'SMI') THEN
          PRINT *, 'shortName = ', yname, ' should not get here!'
          STOP 'Error! Unclear if shortName is SMI or not SMI'

        ELSE  ! yname = 'SMI'

          SELECT CASE ( which_way)
          CASE ( 'smi2w_so')
!           From SMI to W_SO
            WHERE ( soiltyp(:) > 2 .AND. soiltyp(:) < 9 ) 
              zgrib(:) = d*( zgrib(:)*( Cfcap(soiltyp(:))-Cpwp(soiltyp(:)) ) + Cpwp(soiltyp(:)) )
            ELSE WHERE
              zgrib(:) = 0._ireals
            END WHERE
            IF ( oc == 78) CALL grib_set( igrib, 'shortName', 'W_SO')

          CASE ( 'smi2w_solimit')
!           From SMI to W_SO and limit W_SO by Cadp, Cporv
            WHERE ( soiltyp(:) > 2 .AND. soiltyp(:) < 9 ) 
              zgrib(:) = MAX( d*Cadp(soiltyp(:)), MIN( d*Cporv(soiltyp(:)),  &
                              d*( zgrib(:)*( Cfcap(soiltyp(:))-Cpwp(soiltyp(:)) ) + Cpwp(soiltyp(:)) ) ))
            ELSE WHERE
              zgrib(:) = 0._ireals
            END WHERE
            IF ( oc == 78) CALL grib_set( igrib, 'shortName', 'W_SO')

          CASE ( 'smi2w_sodef')
!           From SMI to W_SO, limit W_SO by Cadp, Cporv, and set missing values over land
!           to mean values for that soil type
            w_so_mean(:) = 0._ireals
            n_mean(:)    = 0
            n_miss(:)    = 0
            DO i = 1, nn
              IF ( zgrib(i) /= wmiss) THEN
                isoil = soiltyp(i)
                IF ( isoil > 2 .AND. isoil < 9) THEN
                  zgrib(i) = MAX( d*Cadp(isoil), MIN( d*Cporv(isoil),  &
                             d*( zgrib(i)*( Cfcap(isoil)-Cpwp(isoil) ) + Cpwp(isoil) ) ))
                  n_mean(isoil)    = n_mean(isoil) + 1
                  w_so_mean(isoil) = w_so_mean(isoil) + zgrib(i)
                ELSE
                  zgrib(i) = 0._ireals
                END IF
              END IF
            END DO

            DO isoil = 3, 8
              IF ( n_mean(isoil) > 0) THEN
                w_so_mean(isoil) = w_so_mean(isoil)/REAL(n_mean(isoil),KIND=ireals)
              ELSE
!!            w_so_mean(isoil) = d*0.5_ireals*( Cadp(isoil) + Cporv(isoil) )
                w_so_mean(isoil) = d*0.5_ireals*( Cfcap(isoil) + Cpwp(isoil) )
              END IF
            END DO

            DO i = 1, nn
              IF ( zgrib(i) == wmiss) THEN
                zgrib(i) = w_so_mean(soiltyp(i))
                n_miss(soiltyp(i)) = n_miss(soiltyp(i))+1
              END IF
            END DO
            PRINT *, ' Missing values replaced by means:'
            DO i = 3, 8
              WRITE(*,'(i4,tr2,A10,i6,A,F10.1)') i, ysoiltyp(i), n_miss(i), ' misses set to ', w_so_mean(i)
            END DO
            IF ( oc == 78) CALL grib_set( igrib, 'shortName', 'W_SO')
          END SELECT
        ENDIF

      CASE ( 'W_SO_ICE')
        IF ( .NOT. l_soiltyp) THEN
          PRINT *, 'Error in smi_w_so: SOILTYP not known, yet! Cannot check W_SO_ICE!'
          STOP 'Error! SOILTYP not known!'
        END IF
!!      CALL check_grid( igrib, gridTypeST, uuidST)
        CALL check_grid( igrib)
!       Enforce value 0 for water, ice, or rock points
        WHERE ( soiltyp(:) <= 2 .OR. soiltyp(:) >= 9 ) zgrib(:) = 0._ireals

      END SELECT

      CALL grib_set( igrib, 'values', zgrib)
      CALL grib_write(igrib,nu_out)
      CALL grib_release( igrib)

  ENDDO
!
!   End of file
!
  CALL grib_release( igrib)
  CALL grib_close_file( nu_in, ierr)
  IF (ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'Error ', ierr, ' closing GRIB file ', TRIM(gribfile_in)
  END IF
! IF ( ldebug) PRINT *, 'Closed file: ', TRIM(gribfile_in)
  CALL grib_close_file( nu_out, ierr)
  IF (ierr /= GRIB_SUCCESS) THEN
    PRINT *, 'Error ', ierr, ' closing GRIB file ', TRIM(gribfile_out)
  END IF

  CONTAINS

!!SUBROUTINE check_grid( igrib, gridTypeST, uuidST)
  SUBROUTINE check_grid( igrib)
!
!  Check if gridType and UUID are the same as the input values
!
    INTEGER, INTENT(IN)           :: igrib
!!  CHARACTER(LEN=32), INTENT(IN) :: gridTypeST, uuidST

    CHARACTER(LEN=32)             :: gridType
    CHARACTER(LEN=80)             :: uuid

    CALL grib_get( igrib, 'gridType', gridType)
    IF ( gridType /= gridTypeST) THEN
      PRINT *, 'Error! This GRIB record has gridType ', TRIM(gridType), &
                     ' SOILTYP has gridType ', TRIM(gridTypeST)
      PRINT *, 'No transformation between W_SO and SMI possible'
      STOP 'Different grids for SOILTYP and W_SO or W_SO_ICE'
    END IF

    IF ( gridType == 'unstructured_grid') THEN
      CALL grib_get( igrib, 'uuidOfHGrid', uuid)
      IF ( uuid /= uuidST) THEN
        PRINT *, 'Error! Different uuidOfHGrid for SOILTYP and W_SO, W_SO_ICE, or SMI! '
        PRINT *, 'Cannot convert between W_SO and SMI.'
        STOP 'Different UUID for SOILTYP and W_SO or W_SO_ICE'
      END IF
    END IF
  END SUBROUTINE check_grid

END PROGRAM smi_w_so
