!> Program to scale the LUH2 transitions
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Scaling of LUH2 transitions
!>

! -- to be executed after scaling of the LUH2 states (scale_with_veg-ratio-max.bash)
!    [which is to be executed after processing the LUH2 data in the "CMIP5-way"
!   The scaling of states is the "CMIP6-way":
!   In order to conserve as much as possible of the LUH2 crop + pasture area
!   the states S_t are scaled to S'_t using the veg_ratio_max,
!   i.e. crop and pasture are scaled up (divided by veg_ratio_max)
!   The additional area is then subtracted from the natural vegetation
!   (the desert fraction used in the jsbach simulation is assumed to only be part of LUH2 natural vegetation
!   -- this scaling is only done up to a remaining natural vegetation fraction of 0.01
!   => for years and cells where nat < 0.01 no scaling is done, even if veg_ratio_max < 1)
!
!   Scaling of states makes the previously absolute area fractions of crop and pasture 'relative' to JSBACHs desert fraction.
!   The aim of this program is to also scale the transitions accordingly, thus to also make them 'relative' to the deserts
!   (Note, after the 'CMIP5' like processing the transitions are relative to the unscaled states,
!   thus first the absolute unscaled transitions are calculated, then they are scaled with the desert fraction
!   and then they are made relative to the area fraction of the source vegetation type of the transition again)
!
!   Scaling of the transitions is done in several steps
!   - unscaled transtions:
!       1. get the cyclic part (a. identify potential bilateral ones, b. identify potential remaining full cycles)
!       2. calculate net from absolute and cyclic transitions
!       3. calculate 'direct' transitions from the difference in the states
!       4. determine the 'error' by comparing net and direct transitions -> 'correction factor'
!   - scaled transitions
!       1. calculate direct transitions from the difference in the scaled states
!       2. scale the correction factor
!       3. assemble the scaled transitions: direct transitions + correction + cyclic
!
!   + all parts of the scaled transitions are subject to boundary conditions related with the remaining amount of natural vegetation
!   + transitions T_t need to be adapted to T'_t such that states and transitions are still consistent, i.e. S'_t+1 = S'_t + T'_t
!   + Additional aim: where possible the total area taking part in the transitions should also be conserved
!
!
! NOTE: uses a namelist called "scaling_of_transition"
!
!---- inspired by old fortran programs of TR/SW ++ see discussions with TR, e.g. TRs emails from 02.12.2016 and 27.01.2017
!-- JN 09.12.16 -- julia.nabel@mpimet.mpg.de
!-- JN 31.07.17: changed check on remaining natveg vs delta transitions (during testing rangelands as pastures for Trendy v6)
!
! compile locally:
! nagfor -colour -w=uep Adapting_transitions_after_scaling_of_states.f90 -o Adapting_transitions_after_scaling_of_states.x -I/sw/jessie-x64/netcdf_fortran-4.4.2-static-nag60/include -L/sw/jessie-x64/netcdf-4.3.3.1-static-gccsys/lib -L/sw/jessie-x64/netcdf_fortran-4.4.2-static-nag60/lib -lnetcdff -lnetcdf -L/sw/jessie-x64/hdf5-1.8.16-static-gccsys/lib -lhdf5_hl -lhdf5 -lz -L/sw/jessie-x64/szip-2.1-static-gccsys/lib/ -lsz -ldl
PROGRAM Adapting_transitions_after_scaling_of_states
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  LOGICAL  ::  INFO = .FALSE.

  ! dp & miss: cf. TRs scripts
  INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12,307)

  ! parameters
  INTEGER, PARAMETER  :: NTRANS = 6           ! number of transitions
  INTEGER, PARAMETER  :: NTRANSSTATES = 3     ! number of states in transitions
  INTEGER, PARAMETER  :: NSTATES = 4          ! number of states in state files (here primary and secondary are split)
  INTEGER, PARAMETER  :: NDIMSOFVARS = 3      ! expected number of dimensions of the variables (expected lat,lon and time...)

  REAL(dp), PARAMETER :: FRACT_SMALL = 1.e-10_dp ! from jsbach
  REAL(dp), PARAMETER :: MISS_VALUE = 2.0000000400817547E+20 ! missing value in state files

  ! Minimum fraction of natural vegetation that was conserved in scaled cells
  REAL(dp), PARAMETER :: MIN_NAT_FRACT_IN_SCALING = 0.01

  ! NOTE: Unfortunately the accuracy of the LUH data is only about ~1.e-6,
  ! Louise Parsons Chini - Email 05.01.17: "these differences are on the order of 1e-5 so we were not too concerned about these
  !                                         differences since the accuracy/precision of our datasets is about that order anyway."
  ! Thus, e.g.
  ! - changes in states do not equalise from one time step to another but differ of up to +-1.e-6
  ! - there are changes in the states in the magnitude of 1.e-6 although transitions are zero
  REAL(dp), PARAMETER :: EXPECTED_ACCURACY = 1.e-6 ! Note: needed to be slightly increased for GCB2021 (3.e-5)
  ! - NOTE: depending on the accepted error, imbalances between transitions and states might add up
  REAL(dp), PARAMETER :: RESULTING_ACCURACY_ERRORS = 3.e-5 !5.e-6
  REAL(dp), PARAMETER :: SMALL_RELATIVE_ERROR = 1.e-3

  ! Threshold for detours to be considered (i.e. do not consider unscaled detours below this threshold)
  REAL(dp), PARAMETER :: DETOUR_THRESHOLD = EXPECTED_ACCURACY

  ! resolution dependant
  INTEGER :: nLon           ! number of longitudes
  INTEGER :: nLat           ! number of latitudes

  ! error status return
  INTEGER :: errorStatus

  ! netCDF id
  INTEGER :: ncFileID

  ! dimensions
  INTEGER :: nVar
  INTEGER, ALLOCATABLE, DIMENSION(:) :: dimLength
  CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: dimName
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NDIMINID

  ! variables
  INTEGER :: writeVarIDs(NTRANS+NDIMSOFVARS), writeDimIDs(NDIMSOFVARS)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: varTyp
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: varDimIDs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nVarAttr
  CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:) :: varName

  ! attributes
  INTEGER :: maxNrOfReadAttrs
  CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:,:) :: attrName
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: attrType
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: attrLen
  CHARACTER(LEN=128), ALLOCATABLE, DIMENSION(:,:) :: attrContentChar
  INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: attrContentInt
  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: attrContentReal

  ! data
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: readLons, readLats, readTime
  REAL(dp),ALLOCATABLE, DIMENSION(:,:) :: veg_ratio_max
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: readData, writeData
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: readTrans, scaledTrans
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: readScaledSThisY, readScaledSNextY
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: readUnscaledSThisY, readUnscaledSNextY
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: sortedScaledSThisY, sortedScaledSNextY, deltaScaledS
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: sortedUnscaledSThisY, sortedUnscaledSNextY, deltaUnscaledS
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:,:) :: sortedTrans
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:,:) :: absTrans_UnscaledS, absTrans_deltaUnscaledS, absTrans_netUnscaledS
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:,:) :: absTrans_circ_bilateral, absTrans_circ_fullCycles, absTrans_circ_via_reduction
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:,:) :: trans_bounded_circ
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:,:) :: absTrans_diffNetDelta_UnscaledS, trans_scaled_diffNetDelta
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:,:) :: trans_deltaScaledS, trans_netScaledS, trans_ScaledS

  ! indices, counter, other ints
  INTEGER :: lon, lat, i, j, yearCounter, indexTransition
  INTEGER :: indexCrop2Pasture, indexPasture2Crop, indexNatural2Crop
  INTEGER :: indexNatural2Pasture, indexCrop2Natural, indexPasture2Natural
  INTEGER :: indexOfPrimary, indexOfSecondary, indexOfCrop, indexOfPasture
  INTEGER :: year_start, year_end, nyear, resolution, nrOfVarsVegMaxFile, nrOfDimsVegMaxFile
  INTEGER :: countTransWhereStateIsNAN, countScaledDetoursLeadingToNegTrans

  ! work path
  CHARACTER(LEN=512) :: workingPath

  ! file names and other chars
  CHARACTER(LEN=512) :: vegRatioMaxFileNameWithPath
  CHARACTER(LEN=128) :: folderStateFiles      ! directory with scaled states files
  CHARACTER(LEN=128) :: folderUnscaledStateFiles
  CHARACTER(LEN=128) :: folderTransitionFiles ! directory with unscaled transition files
  CHARACTER(LEN=128) :: outputFolder          ! directory for scaled transition files
  CHARACTER(LEN=128) :: stateTagScaled, stateTagUnscaled, transitionsTag
  CHARACTER(LEN=512) :: fileName, outFileName
  CHARACTER(LEN=4) :: char_year, thisYYYY, nextYYYY
  CHARACTER(LEN=3) :: char_res
  CHARACTER(LEN=15) :: vrmVarName
  CHARACTER(LEN=8) :: date
  CHARACTER(LEN=10) :: time
  CHARACTER(LEN=256) :: calledFrom
  CHARACTER(LEN=128) :: luhRelease
  CHARACTER(LEN=512) :: history

  ! logicals
  LOGICAL :: checkCrop2Pasture, checkPasture2Crop, checkNatural2Crop
  LOGICAL :: checkNatural2Pasture, checkCrop2Natural, checkPasture2Natural
  LOGICAL :: wasModified

  !--------------------------------------------------------------------------------------------------------------------------
  NAMELIST /scaling_of_transition/ year_start, year_end, resolution, workingPath, vegRatioMaxFileNameWithPath, vrmVarName, &
      & nrOfVarsVegMaxFile, nrOfDimsVegMaxFile, folderStateFiles, folderUnscaledStateFiles, folderTransitionFiles, outputFolder, &
      & stateTagScaled, stateTagUnscaled, transitionsTag, calledFrom, luhRelease

  ! Read time info from namelist
  resolution = 63
  year_start = -1
  year_end   = -2
  workingPath = "."
  vegRatioMaxFileNameWithPath = ""
  vrmVarName = ""
  nrOfVarsVegMaxFile = -1
  nrOfDimsVegMaxFile = -1
  folderStateFiles = "/"
  folderUnscaledStateFiles = "/"
  folderTransitionFiles = "/"
  outputFolder = "/"
  stateTagScaled = ""
  stateTagUnscaled = ""
  transitionsTag = ""
  calledFrom = ""
  luhRelease = ""
  open(11,file='scaling_of_transition.nml',status='old',position='rewind')
  read(11,scaling_of_transition,iostat=i)
  if (i /= 0) then
    write(*,*) '"scaling_of_transition" namelist not found or corrupted', i
    stop
  endif
  close(11)
  if (year_end < year_start) then
    write(*,*) 'Bad combination of namelist parameters'
    stop
  endif

  year_start = year_start - 1            ! To be adjusted in the time loop
  nyear      = year_end - year_start

  write(*,*) 'Work path: ',trim(workingPath)
  write(*,*) 'Processing years: ',year_start+1,year_end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Gauss resolution to expected LAT/LON !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (resolution == 21) THEN; nLat = 32
  ELSEIF (resolution == 31) THEN; nLat = 48
  ELSEIF (resolution == 42) THEN; nLat = 64
  ELSEIF (resolution == 63) THEN; nLat = 96
  ELSEIF (resolution == 85) THEN; nLat = 128
  ELSEIF (resolution == 106) THEN; nLat = 160
  ELSEIF (resolution == 127) THEN; nLat = 192
  ELSEIF (resolution == 159) THEN; nLat = 240
  ELSEIF (resolution == 255) THEN; nLat = 384
  ELSEIF (resolution == 319) THEN; nLat = 480
  ELSE
    STOP 'Resolution not prepared!'
  END IF
  nLon = nLat * 2

  IF (resolution <= 99) THEN
    WRITE(char_res,'(i2)')  resolution
  ELSE IF (resolution >=100 .AND. resolution <= 999) THEN
    WRITE(char_res,'(i3)')  resolution
  ELSE
    STOP 'CHAR for resolution not prepared!'
  END IF

  write(*,*) 'Expected resolution: T',resolution,'LONS',nLon,'LATS',nLat

  !----- allocate data arrays
  ALLOCATE(veg_ratio_max(nLon,nLat))
  ALLOCATE(readData(nLon,nLat,1))

  ALLOCATE(readScaledSNextY(nLon,nLat,NSTATES))
  ALLOCATE(readScaledSThisY(nLon,nLat,NSTATES))
  ALLOCATE(readUnscaledSNextY(nLon,nLat,NSTATES))
  ALLOCATE(readUnscaledSThisY(nLon,nLat,NSTATES))
  ALLOCATE(readTrans(nLon,nLat,NTRANS))
  ALLOCATE(scaledTrans(nLon,nLat,NTRANS))
  ALLOCATE(sortedScaledSThisY(nLon,nLat,NTRANSSTATES))
  ALLOCATE(sortedScaledSNextY(nLon,nLat,NTRANSSTATES))
  ALLOCATE(sortedUnscaledSThisY(nLon,nLat,NTRANSSTATES))
  ALLOCATE(sortedUnscaledSNextY(nLon,nLat,NTRANSSTATES))
  ALLOCATE(deltaScaledS(nLon,nLat,NTRANSSTATES))
  ALLOCATE(deltaUnscaledS(nLon,nLat,NTRANSSTATES))
  ALLOCATE(sortedTrans(nLon,nLat,NTRANSSTATES,NTRANSSTATES))

  ALLOCATE(absTrans_UnscaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(absTrans_deltaUnscaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(absTrans_netUnscaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))

  ALLOCATE(absTrans_diffNetDelta_UnscaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(trans_scaled_diffNetDelta(nLon,nLat,NTRANSSTATES,NTRANSSTATES))

  ALLOCATE(absTrans_circ_bilateral(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(absTrans_circ_fullCycles(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(absTrans_circ_via_reduction(nLon,nLat,NTRANSSTATES,NTRANSSTATES))

  ALLOCATE(trans_bounded_circ(nLon,nLat,NTRANSSTATES,NTRANSSTATES))

  ALLOCATE(trans_deltaScaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(trans_netScaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))
  ALLOCATE(trans_ScaledS(nLon,nLat,NTRANSSTATES,NTRANSSTATES))

  veg_ratio_max(:,:) = 0.
  readData(:,:,:) = 0.
  readScaledSNextY(:,:,:) = 0.
  readScaledSThisY(:,:,:) = 0.
  readUnscaledSNextY(:,:,:) = 0.
  readUnscaledSThisY(:,:,:) = 0.
  readTrans(:,:,:) = 0.

  sortedUnscaledSThisY(:,:,:) = 0.
  sortedScaledSThisY(:,:,:) = 0.

  !----- get veg_ratio_max
  ALLOCATE (varName(nrOfVarsVegMaxFile))
  ALLOCATE (varTyp(nrOfVarsVegMaxFile))
  ALLOCATE (nVarAttr(nrOfVarsVegMaxFile))
  ALLOCATE (dimLength(nrOfDimsVegMaxFile))
  ALLOCATE (dimName(nrOfDimsVegMaxFile))
  ALLOCATE (NDIMINID(nrOfDimsVegMaxFile))
  ALLOCATE (varDimIDs(nrOfVarsVegMaxFile,nrOfDimsVegMaxFile))
  CALL prepare_input_file(trim(vegRatioMaxFileNameWithPath), nrOfVarsVegMaxFile, &
      & nrOfDimsVegMaxFile, ncFileID, varName, varTyp, nVarAttr, dimLength, varDimIDs, dimName)
  DO i = 1,nrOfVarsVegMaxFile
    IF (TRIM(varName(i)) == vrmVarName ) THEN
      errorStatus = nf_get_var_double(ncFileID,i,readData)
      CALL check_err(errorStatus)
      call check_dim_order(nrOfVarsVegMaxFile, varDimIDs(i,:), dimName, NDIMINID)
      if(ALL(NDIMINID .eq. (/ 1,2,3 /))) THEN
        veg_ratio_max = readData(:,:,1)
      else
        write(*,*) TRIM(varName(i)), NDIMINID
        write(*,*) nrOfVarsVegMaxFile
        write(*,*) varDimIDs
        write(*,*) dimName
        STOP "Wrong order of dimensions! Expected: lon, lat, time"
      endif
    ENDIF
  ENDDO

  DEALLOCATE(varName, varTyp, nVarAttr, dimLength, dimName, varDimIDs)

  !----- time loop
  DO yearCounter = 1,nyear

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! preparation: init/reset arrays !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sortedUnscaledSNextY(:,:,:) = 0.
    sortedScaledSNextY(:,:,:) = 0.
    sortedTrans(:,:,:,:) = 0.

    deltaScaledS(:,:,:) = 0.
    deltaUnscaledS(:,:,:) = 0.

    absTrans_UnscaledS(:,:,:,:) = 0.
    absTrans_deltaUnscaledS(:,:,:,:) = 0.
    absTrans_netUnscaledS(:,:,:,:) = 0.

    absTrans_circ_bilateral(:,:,:,:) = 0.
    absTrans_circ_fullCycles(:,:,:,:) = 0.
    absTrans_circ_via_reduction(:,:,:,:) = 0.
    trans_bounded_circ(:,:,:,:) = 0.

    absTrans_diffNetDelta_UnscaledS(:,:,:,:) = 0.
    trans_scaled_diffNetDelta(:,:,:,:) = 0.

    trans_deltaScaledS(:,:,:,:) = 0.
    trans_netScaledS(:,:,:,:) = 0.
    trans_ScaledS(:,:,:,:) = 0.

    scaledTrans(:,:,:) = 0.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! preparation: get annual data !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !----- format year (jsbach files have years with 4 digits)
    IF (year_start + yearCounter .LT. 999) THEN
      WRITE(char_year,'(i3)') year_start + yearCounter
      thisYYYY(1:1) = '0'
      thisYYYY(2:4) = TRIM(char_year)
      WRITE(char_year,'(i3)') year_start + yearCounter + 1
      nextYYYY(1:1) = '0'
      nextYYYY(2:4) = TRIM(char_year)
    ELSEIF (year_start + yearCounter .EQ. 999) THEN
      WRITE(char_year,'(i3)') year_start + yearCounter
      thisYYYY(1:1) = '0'
      thisYYYY(2:4) = TRIM(char_year)
      WRITE(char_year,'(i4)') year_start + yearCounter + 1
      nextYYYY(1:4) = TRIM(char_year)
    ELSEIF (year_start + yearCounter .GT. 9999) THEN
      STOP 'Only input files with years <= 9999 expected!'
    ELSEIF (year_start + yearCounter .LE. 99) THEN
      STOP 'Only input files with years >= 99 expected!'
    ELSE
      WRITE(char_year,'(i4)') year_start + yearCounter
      thisYYYY(1:4) = TRIM(char_year)
      WRITE(char_year,'(i4)') year_start + yearCounter + 1
      nextYYYY(1:4) = TRIM(char_year)
    END IF


    !----- get land cover states of next year
    fileName = trim(workingPath) // trim(folderStateFiles) // 'LUH_states_T' // &
        & TRIM(char_res) // TRIM(stateTagScaled) // '_' // TRIM(nextYYYY) // '.nc'
    CALL read_states_file(fileName, nLon, nLat, readScaledSNextY)

    fileName = trim(workingPath) // trim(folderUnscaledStateFiles) // 'LUH_states_T' // &
        & TRIM(char_res) // TRIM(stateTagUnscaled) // '_' // TRIM(nextYYYY) // '.nc'
    CALL read_states_file(fileName, nLon, nLat, readUnscaledSNextY)

    ! in the first year we have to additionally get the state data of the current year
    IF (yearCounter .EQ. 1) THEN
      fileName = trim(workingPath) // trim(folderStateFiles) // 'LUH_states_T' // &
          & TRIM(char_res) // TRIM(stateTagScaled) // '_' // TRIM(thisYYYY) // '.nc'
      CALL read_states_file(fileName, nLon, nLat, readScaledSThisY)

      fileName = trim(workingPath) // trim(folderUnscaledStateFiles) // 'LUH_states_T' // &
          & TRIM(char_res) // TRIM(stateTagUnscaled) // '_' // TRIM(thisYYYY) // '.nc'
      CALL read_states_file(fileName, nLon, nLat, readUnscaledSThisY)

    ENDIF


    !----- get transitions (a lot of the file information is used later for the new file with scaled transitions)
    fileName = trim(workingPath) // trim(folderTransitionFiles) &
        & // 'LUH_transitions_T' // TRIM(char_res) // TRIM(transitionsTag) // '_' // TRIM(thisYYYY) // '.nc'

    nVar = NTRANS + NDIMSOFVARS
    ALLOCATE (varName(nVar))
    ALLOCATE (varTyp(nVar))
    ALLOCATE (nVarAttr(nVar))
    ALLOCATE (dimLength(NDIMSOFVARS))
    ALLOCATE (dimName(NDIMSOFVARS))
    ALLOCATE (varDimIDs(nVar,NDIMSOFVARS))

    CALL prepare_input_file(TRIM(fileName), nVar, NDIMSOFVARS, ncFileID, varName, varTyp, nVarAttr, dimLength, varDimIDs, dimName)

    ! get attribute names and then attribute types and lengths
    maxNrOfReadAttrs=0
    DO i=1,nVar
      IF (nVarAttr(i).GT.maxNrOfReadAttrs) maxNrOfReadAttrs=nVarAttr(i)
    END DO
    IF (INFO) WRITE (*,*) 'maxNrOfReadAttrs', maxNrOfReadAttrs
    ALLOCATE (attrName(nVar,maxNrOfReadAttrs))
    ALLOCATE (attrType(nVar,maxNrOfReadAttrs))
    ALLOCATE (attrLen(nVar,maxNrOfReadAttrs))
    DO i=1,nVar
      DO j=1,nVarAttr(i)
        errorStatus = nf_inq_attname(ncFileID,i,j,attrName(i,j))
        CALL check_err(errorStatus)
        errorStatus = nf_inq_att(ncFileID,i,attrName(i,j),attrType(i,j),attrLen(i,j))
        CALL check_err(errorStatus)
        IF (INFO) WRITE (*,*) 'i,j,attrNameE,ATTTYP,attrLen',i,j,TRIM(attrName(i,j)),attrType(i,j),attrLen(i,j)
      END DO
    END DO
    ! get attributes
    ALLOCATE (attrContentChar(nVar,maxNrOfReadAttrs))
    ALLOCATE (attrContentInt(nVar,maxNrOfReadAttrs))
    ALLOCATE (attrContentReal(nVar,maxNrOfReadAttrs))
    attrContentChar(:,:) = ''
    DO i=1,nVar
      DO j=1,nVarAttr(i)
        IF (attrType(i,j).EQ.2) THEN
          errorStatus = nf_get_att_text(ncFileID,i,attrName(i,j),attrContentChar(i,j))
          CALL check_err(errorStatus)
          IF (INFO) WRITE (*,*) 'i,attrName,attrContentChar',i,TRIM(attrName(i,j)),TRIM(attrContentChar(i,j))
        ELSE IF (attrType(i,j).EQ.4) THEN
          errorStatus = nf_get_att_int(ncFileID,i,attrName(i,j),attrContentInt(i,j))
          CALL check_err(errorStatus)
          IF (INFO) WRITE (*,*) 'i,attrName,attrContentInt',i,TRIM(attrName(i,j)),attrContentInt(i,j)
        ELSE IF ((attrType(i,j).EQ.5) .OR. (attrType(i,j).EQ.6)) THEN
          errorStatus = nf_get_att_double(ncFileID,i,attrName(i,j),attrContentReal(i,j))
          CALL check_err(errorStatus)
          IF (INFO) WRITE (*,*) 'i,attrName,attrContentReal',i,TRIM(attrName(i,j)),attrContentReal(i,j)
        ELSE
          WRITE (*,*) 'Type of attribute not prepared in the code, failed to read grid',i,j
          STOP
        END IF
      END DO
    END DO
    ! get data and check for expected transitions
    checkCrop2Pasture = .false.
    checkPasture2Crop = .false.
    checkNatural2Crop = .false.
    checkNatural2Pasture = .false.
    checkCrop2Natural = .false.
    checkPasture2Natural = .false.
    indexTransition = 0
    DO i = 1,nVar
      IF (varName(i) == 'lon' ) THEN
        ALLOCATE(readLons(dimLength(varDimIDs(i,1))))
        errorStatus = nf_get_var_double(ncFileID,i,readLons)
        CALL check_err(errorStatus)
      ELSE IF (varName(i) == 'lat') THEN
        ALLOCATE(readLats(dimLength(varDimIDs(i,1))))
        errorStatus = nf_get_var_double(ncFileID,i,readLats)
        CALL check_err(errorStatus)
      ELSE IF (varName(i) == 'time') THEN
        ALLOCATE(readTime(dimLength(varDimIDs(i,1))))
        errorStatus = nf_get_var_double(ncFileID,i,readTime)
        CALL check_err(errorStatus)
      ELSE
        indexTransition = indexTransition + 1
        IF (varName(i) == 'crop2pasture') THEN
          checkCrop2Pasture = .true.
          indexCrop2Pasture = indexTransition
        ELSE IF (varName(i) == 'pasture2crop') THEN
          checkPasture2Crop = .true.
          indexPasture2Crop = indexTransition
        ELSE IF (varName(i) == 'natural2crop') THEN
          checkNatural2Crop = .true.
          indexNatural2Crop = indexTransition
        ELSE IF (varName(i) == 'natural2pasture') THEN
          checkNatural2Pasture = .true.
          indexNatural2Pasture = indexTransition
        ELSE IF (varName(i) == 'crop2natural') THEN
          checkCrop2Natural = .true.
          indexCrop2Natural = indexTransition
        ELSE IF (varName(i) == 'pasture2natural') THEN
          checkPasture2Natural = .true.
          indexPasture2Natural = indexTransition
        END IF

        errorStatus = nf_get_var_double(ncFileID,i,readData)
        CALL check_err(errorStatus)
        readTrans(:,:,indexTransition) = readData(:,:,1)
      END IF
    END DO

    IF (indexTransition /= NTRANS) STOP 'Unexpected number of transitions!'
    IF (.NOT. checkCrop2Pasture) STOP 'Variable crop2pasture not found'
    IF (.NOT. checkPasture2Crop) STOP 'Variable pasture2crop not found'
    IF (.NOT. checkNatural2Crop) STOP 'Variable natural2pasture not found'
    IF (.NOT. checkNatural2Pasture) STOP 'Variable natural2pasture not found'
    IF (.NOT. checkCrop2Natural) STOP 'Variable crop2natural not found'
    IF (.NOT. checkPasture2Natural) STOP 'Variable pasture2natural not found'

    errorStatus = nf_close(ncFileID)
    CALL check_err(errorStatus)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! "Main part:" scale transitions !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !--- sort the read 4 states in the 3 trans states {n,c,p}
    sortedScaledSNextY(:,:,1) = readScaledSNextY(:,:,indexOfPrimary) + readScaledSNextY(:,:,indexOfSecondary)  ! natural
    sortedScaledSNextY(:,:,2) = readScaledSNextY(:,:,indexOfCrop)                                        ! crop
    sortedScaledSNextY(:,:,3) = readScaledSNextY(:,:,indexOfPasture)                                     ! pasture

    sortedUnscaledSNextY(:,:,1) = readUnscaledSNextY(:,:,indexOfPrimary) + readUnscaledSNextY(:,:,indexOfSecondary)  ! natural
    sortedUnscaledSNextY(:,:,2) = readUnscaledSNextY(:,:,indexOfCrop)                                        ! crop
    sortedUnscaledSNextY(:,:,3) = readUnscaledSNextY(:,:,indexOfPasture)

    ! TODO: revisit: numerical issues can lead to natural vegetation with cf > 1
    CALL check_natural_veg_for_numerical_issues(sortedUnscaledSNextY, 'unscaled states next year')
    CALL check_natural_veg_for_numerical_issues(sortedScaledSNextY, 'scaled states next year')

    IF (yearCounter .EQ. 1) THEN
      sortedScaledSThisY(:,:,1) = readScaledSThisY(:,:,indexOfPrimary) + readScaledSThisY(:,:,indexOfSecondary)  ! natural
      sortedScaledSThisY(:,:,2) = readScaledSThisY(:,:,indexOfCrop)                                        ! crop
      sortedScaledSThisY(:,:,3) = readScaledSThisY(:,:,indexOfPasture)                                     ! pasture

      sortedUnscaledSThisY(:,:,1) = readUnscaledSThisY(:,:,indexOfPrimary) + readUnscaledSThisY(:,:,indexOfSecondary)  ! natural
      sortedUnscaledSThisY(:,:,2) = readUnscaledSThisY(:,:,indexOfCrop)                                        ! crop
      sortedUnscaledSThisY(:,:,3) = readUnscaledSThisY(:,:,indexOfPasture)                                     ! pasture

      ! TODO: revisit: numerical issues can lead to natural vegetation with cf > 1
      CALL check_natural_veg_for_numerical_issues(sortedUnscaledSThisY, 'unscaled states this year')
      CALL check_natural_veg_for_numerical_issues(sortedScaledSThisY, 'scaled states this year')

    ENDIF


    !--- sort the read 6 3D transitions into a 4D array with 3 entries {{n->n,n->c,n->p},{c->n,c->c,c->p},{p->n,p->c,p->p}}
    sortedTrans = 0.
    !x sortedTrans(:,:,1,1) = 0.0 nat-nat
    sortedTrans(:,:,1,2) = readTrans(:,:,indexNatural2Crop)
    sortedTrans(:,:,1,3) = readTrans(:,:,indexNatural2Pasture)
    sortedTrans(:,:,2,1) = readTrans(:,:,indexCrop2Natural)
    !x sortedTrans(:,:,2,2) = 0.0 crop-crop
    sortedTrans(:,:,2,3) = readTrans(:,:,indexCrop2Pasture)
    sortedTrans(:,:,3,1) = readTrans(:,:,indexPasture2Natural)
    sortedTrans(:,:,3,2) = readTrans(:,:,indexPasture2Crop)
    !x sortedTrans(:,:,3,3) = 0.0 pastr-pastr
    !ASSERTION: sum of outgoing transition of any state may not exceed one
    CALL check_sum_of_transitions(sortedTrans, sortedUnscaledSThisY, sortedUnscaledSNextY, 'org relative trans')

    !--- get absolute transitions i -> j (but only were states not NAN)
    countTransWhereStateIsNAN = COUNT(ANY(sortedUnscaledSThisY(:,:,:) .EQ. MISS_VALUE, DIM=3) &
        & .AND. ANY(ANY(sortedTrans(:,:,:,:) .NE. 0., DIM=4), DIM=3))
    IF(countTransWhereStateIsNAN .GT. 0) THEN
      WRITE(*,*) '### Counted non zero transitions at NAN states: ', countTransWhereStateIsNAN
    ENDIF
    FORALL(lon=1:nLon, lat=1:nLat, i=1:NTRANSSTATES, ALL(sortedUnscaledSThisY(lon,lat,:) .NE. MISS_VALUE))
      absTrans_UnscaledS(lon,lat,i,:) = sortedTrans(lon,lat,i,:) * sortedUnscaledSThisY(lon,lat,i)
    END FORALL
    !ASSERTION: sum of outgoing transition of any state may not exceed one
    CALL check_sum_of_transitions(absTrans_UnscaledS, sortedUnscaledSThisY, sortedUnscaledSNextY, 'org absolute trans')

    !--- calculate delta for scaled and unscaled states
    deltaUnscaledS = sortedUnscaledSNextY - sortedUnscaledSThisY
    deltaScaledS = sortedScaledSNextY - sortedScaledSThisY
    ! ASSERTION: changes in states need to equal out
    CALL check_delta_of_states(deltaUnscaledS, 'unscaled')
    CALL check_delta_of_states(deltaScaledS, 'scaled')
    ! calculate direct transitions from delta
    CALL calculate_transitions_from_deltaS(absTrans_UnscaledS, deltaUnscaledS, absTrans_deltaUnscaledS)
    ! ASSERTION: delta and transitions need to be consistent
    ! - TODO: reconsider: under certain circumstances transitions are replaced by delta transitions in case of inconsistency
    wasModified = .FALSE.
    CALL check_consistency_delta_and_transitions(deltaUnscaledS, absTrans_UnscaledS, absTrans_deltaUnscaledS, wasModified)

    !--- disentangle unscaled transitions
    !- get circular transitions from abs transitions via reduction
    CALL calculate_bilateral_and_remaining_cycling(nLon, nLat, &
        & absTrans_UnscaledS, absTrans_circ_bilateral, absTrans_circ_fullCycles)
    absTrans_circ_via_reduction = absTrans_circ_fullCycles + absTrans_circ_bilateral
    CALL check_circularity_of_circular_transitions(absTrans_circ_via_reduction, 'unconstraint circular transitions')

    !- derive net transitions as remaining transitions
    absTrans_netUnscaledS = absTrans_UnscaledS - absTrans_circ_via_reduction
    ! clean net trans (and abs trans) where possible
    CALL remove_unscaled_where_zero(deltaUnscaledS, absTrans_netUnscaledS, absTrans_UnscaledS)

    !- get a correction required for direct transitions (detours)
    ! the difference between the direct transitions from delta and the net transitions remaining when eliminating the circular
    ! transitions from the actual transitions are the detours contained in the actual transitions
    absTrans_diffNetDelta_UnscaledS = absTrans_netUnscaledS - absTrans_deltaUnscaledS
    !TODO: reconsider: detours are only considered if they are larger than a threshold
    FORALL(lon=1:nLon,lat=1:nLat, &
        & ALL(ABS(absTrans_netUnscaledS(lon,lat,:,:) - absTrans_deltaUnscaledS(lon,lat,:,:)) .LE. DETOUR_THRESHOLD))
      absTrans_diffNetDelta_UnscaledS(lon,lat,:,:) = 0.
    END FORALL


    !--- assemble scaled transitions
    !- calculate delta transitions
    ! NOTE: in the scaled case there are no abs trans (since this is the target variable)
    !       therefore, absTrans_UnscaledS is used for plausibility checks (e.g. transitions should not be zero when delta is not)
    CALL calculate_transitions_from_deltaS(absTrans_UnscaledS, deltaScaledS, trans_deltaScaledS)

    !- derive net transitions
    ! Note: as for the unscaled case, delta does not contain the detours from the original transitions
    ! therefore: scale the detours found in the unscaled case and use them as detours for the scaled case
    ! TR: t / veg_ratio_max!
    FORALL(lon=1:nLon,lat=1:nLat, veg_ratio_max(lon,lat) .GT. FRACT_SMALL)
      trans_scaled_diffNetDelta(lon,lat,:,:) = absTrans_diffNetDelta_UnscaledS(lon,lat,:,:) / veg_ratio_max(lon,lat)
    END FORALL
    ! boundary conditions: ignore detours in cells which were scaled but are at the limit after scaling
    !                      + detours can only be fulfilled, if enough natural vegetation is available
    CALL constrain_detours_with_natural_vegetation(sortedScaledSThisY, sortedUnscaledSThisY, &
        & trans_deltaScaledS, trans_scaled_diffNetDelta)

    ! now use the scaled correction to derive the scaled net transitions
    FORALL(lon=1:nLon,lat=1:nLat, ALL(trans_deltaScaledS(lon,lat,:,:) + trans_scaled_diffNetDelta(lon,lat,:,:) .GE. 0.))
      trans_netScaledS(lon,lat,:,:) = trans_deltaScaledS(lon,lat,:,:) + trans_scaled_diffNetDelta(lon,lat,:,:)
    END FORALL
    countScaledDetoursLeadingToNegTrans = COUNT(trans_deltaScaledS + trans_scaled_diffNetDelta .LT. 0.)
    IF ( countScaledDetoursLeadingToNegTrans .GT. 0 ) THEN
      WRITE(*,*) '### Counted detours leading to negative transitions - these were ignored: ', countScaledDetoursLeadingToNegTrans
      FORALL(lon=1:nLon,lat=1:nLat, ANY(trans_deltaScaledS(lon,lat,:,:) + trans_scaled_diffNetDelta(lon,lat,:,:) .LT. 0.))
        trans_netScaledS(lon,lat,:,:) = trans_deltaScaledS(lon,lat,:,:)
      END FORALL
    ENDIF
    ! ASSERTION: transitions are not allowed to be (more) negative (than expected accuracy) and need to be smaller than 1
    CALL check_transitions_for_valid_range(trans_netScaledS, sortedScaledSNextY, 'net scaled transitions')

    !- add net and circular transitions
    trans_ScaledS = trans_netScaledS + absTrans_circ_via_reduction
    ! boundary condition: transitions involving natural vegetation are constrained by the available area
    ! (the veg_ratio_max scaling might reduce natural vegetation strongly)
    ! -> i.e. it is not possible to move more natural vegetation in a year than available in the scaled state
    CALL constrain_with_scaled_natural_vegetation(sortedScaledSThisY, &
        & trans_netScaledS, trans_ScaledS, absTrans_circ_fullCycles, absTrans_circ_bilateral, trans_bounded_circ)
    ! reassemble transitions
    trans_ScaledS = trans_netScaledS + trans_bounded_circ
    ! ASSERTION: transitions are not allowed to be negative and need to be smaller than 1
    CALL check_transitions_for_valid_range(trans_ScaledS, sortedScaledSNextY, 'scaled transitions')
    !ASSERTION: sum of outgoing transition of any state may not exceed one
    CALL check_sum_of_transitions(trans_ScaledS, sortedScaledSThisY, sortedScaledSNextY, 'scaled absolute trans')
    ! check that sum trans_ScaledS still equals delta
    wasModified = .TRUE.
    CALL check_consistency_delta_and_transitions(deltaScaledS, trans_ScaledS, trans_deltaScaledS, wasModified)


    ! Calculate transitions relative to source area fraction
    sortedTrans = 0.
    FORALL(lon=1:nLon,lat=1:nLat,i = 1:NTRANSSTATES, sortedScaledSThisY(lon,lat,i) .GT. FRACT_SMALL )
      sortedTrans(lon,lat,i,:) = trans_ScaledS(lon,lat,i,:) / sortedScaledSThisY(lon,lat,i)
    END FORALL
    ! ASSERTION: transitions are not allowed to be negative and need to be smaller than 1
    CALL check_transitions_for_valid_range(sortedTrans, sortedScaledSNextY, 'final relative transitions')
    !ASSERTION: sum of outgoing transition of any state may not exceed one
    CALL check_sum_of_transitions(sortedTrans, sortedScaledSThisY, sortedScaledSNextY, 'scaled relative trans')

    !-- sort the 4D array back into a 3D array for writing
    scaledTrans(:,:,indexNatural2Crop) = sortedTrans(:,:,1,2)
    scaledTrans(:,:,indexNatural2Pasture) = sortedTrans(:,:,1,3)
    scaledTrans(:,:,indexCrop2Natural) = sortedTrans(:,:,2,1)
    scaledTrans(:,:,indexCrop2Pasture) = sortedTrans(:,:,2,3)
    scaledTrans(:,:,indexPasture2Natural) = sortedTrans(:,:,3,1)
    scaledTrans(:,:,indexPasture2Crop) = sortedTrans(:,:,3,2)

    ! store this next year arrays for the next year in this year arrays ;)
    sortedScaledSThisY = sortedScaledSNextY
    sortedUnscaledSThisY = sortedUnscaledSNextY

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Output file: scaled transitions!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    outFileName = trim(workingPath) // trim(outputFolder) // 'LUH_scaled_transitions_T' &
        & // TRIM(char_res) // '_' // TRIM(thisYYYY) // '.nc'

    WRITE (*,*) 'scaled transitions filename: ',TRIM(outFileName)
    CALL check_err(nf_create(outFileName,NF_CLOBBER,ncFileID))
    ! define dimensions
    DO i = 1,nVar
      IF ((varName(i) == 'lon') .OR. (varName(i) == 'lat')) THEN
        CALL check_err(nf_def_dim(ncFileID,dimName(i),dimLength(i),writeDimIDs(i)))
        CALL check_err(nf_def_var(ncFileID,varName(i),varTyp(i),1,writeDimIDs(i),writeVarIDs(i)))
      ELSEIF (varName(i) == 'time') THEN
        CALL check_err(nf_def_dim(ncFileID,dimName(i),NF_UNLIMITED,writeDimIDs(i)))
        CALL check_err(nf_def_var(ncFileID,varName(i),varTyp(i),1,writeDimIDs(i),writeVarIDs(i)))
      ELSE
        ! define non-dimension variables
        CALL check_err(nf_def_var(ncFileID,varName(i),varTyp(i),NDIMSOFVARS,varDimIDs(i,:),writeVarIDs(i)))
      ENDIF
      CALL PUTATT(nVar,i,nVarAttr,attrType,attrName,attrLen,attrContentChar,ncFileID, &
            & writeVarIDs(i),attrContentInt,attrContentReal,INFO)
    END DO

    !Add global attributes
    CALL check_err(nf_put_att_text(ncFileID, NF_GLOBAL, 'title', 27, 'transitions file for JSBACH'))
    CALL DATE_AND_TIME(date,time)
    history = date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '//time(1:2)//':'//time(3:4)//CHAR(10) &
        & //' created with: '//TRIM(calledFrom)//CHAR(10)//' from: ' // TRIM(luhRelease)
    CALL check_err(nf_put_att_text(ncFileID, NF_GLOBAL, 'history', LEN_TRIM(history), TRIM(history)))

    CALL check_err(nf_enddef(ncFileID))

    ! put longitudes, latitudes and time
    DO i=1, 3
      IF (varName(i) == 'lon') then
        CALL check_err(nf_put_var_double(ncFileID,writeVarIDs(i),readLons))
      elseif (varName(i) == 'lat') THEN
        CALL check_err(nf_put_var_double(ncFileID,writeVarIDs(i),readLats))
      else
        CALL check_err(nf_put_vara_double(ncFileID,writeVarIDs(1),(/1/),(/1/),readTime))
      end if
    end do

    ! put data
    ALLOCATE(writeData(dimLength(2),dimLength(3),1))

    DO i = 1,NTRANS
      writeData(:,:,1) = scaledTrans(:,:,i)
      CALL check_err(nf_put_vara_double(ncFileID,writeVarIDs(i+NDIMSOFVARS),&
          & (/1,1,1/),(/dimLength(2),dimLength(3),1/),writeData(:,:,1)))
    END DO

    ! close the output-file
    CALL check_err(nf_close(ncFileID))

    DEALLOCATE(nVarAttr, attrName, attrType, attrLen, attrContentChar, attrContentInt, attrContentReal)
    DEALLOCATE(readLons, readLats, readTime, writeData)
    DEALLOCATE(dimLength, dimName, varDimIDs, varName, varTyp)

    WRITE (*,*) '__________________________'

  END DO  ! close time loop

  ! deallocate remaining arrays
  DEALLOCATE(readData)
  DEALLOCATE(readScaledSThisY, readScaledSNextY, readUnscaledSThisY, readUnscaledSNextY)
  DEALLOCATE(readTrans, sortedScaledSThisY, sortedScaledSNextY, sortedTrans)
  DEALLOCATE(sortedUnscaledSThisY, sortedUnscaledSNextY, deltaScaledS, deltaUnscaledS)
  DEALLOCATE(absTrans_UnscaledS, absTrans_deltaUnscaledS, absTrans_netUnscaledS)
  DEALLOCATE(absTrans_diffNetDelta_UnscaledS, trans_scaled_diffNetDelta)
  DEALLOCATE(absTrans_circ_bilateral, absTrans_circ_fullCycles, absTrans_circ_via_reduction, trans_bounded_circ)
  DEALLOCATE(trans_deltaScaledS, trans_netScaledS, trans_ScaledS, scaledTrans)


CONTAINS

  ! --- read states file --------------------------------------------------------------------------------------------
  ! reads the state file and detects the indices of the expected states
  SUBROUTINE read_states_file(fileName, nLon, nLat, readStates)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

    !-- args
    CHARACTER(LEN=*),INTENT(IN) :: fileName
    INTEGER,INTENT(IN) :: nLon, nLat
    REAL(dp), INTENT(OUT), DIMENSION(nLon, nLat, NSTATES) :: readStates

    !-- vars
    INTEGER :: errorStatus, i, indexState, ncFileID
    CHARACTER(LEN=128), DIMENSION(NSTATES + NDIMSOFVARS) :: varName
    INTEGER , DIMENSION(NSTATES + NDIMSOFVARS) :: nVarAttr, varTyp
    REAL(dp), DIMENSION(nLon, nLat, 1) :: readData
    INTEGER, DIMENSION(NDIMSOFVARS) :: dimLength, NDIMINID
    CHARACTER(LEN=4), DIMENSION(NDIMSOFVARS) :: dimName
    INTEGER, DIMENSION(NSTATES + NDIMSOFVARS, NDIMSOFVARS) :: varDimIDs
    LOGICAL :: checkPrimary, checkSecondary, checkCrop, checkPastrure

    CALL prepare_input_file(TRIM(fileName), NSTATES + NDIMSOFVARS, NDIMSOFVARS, &
        & ncFileID, varName, varTyp, nVarAttr, dimLength, varDimIDs, dimName)

    ! get data and check for expected states
    checkPrimary = .false.
    checkSecondary = .false.
    checkCrop = .false.
    checkPastrure = .false.
    indexState = 0
    DO i = 1, NSTATES + NDIMSOFVARS
      IF (.NOT. varName(i) == 'lon' .AND. .NOT. varName(i) == 'lat' .AND. .NOT. varName(i) == 'time') THEN
        indexState = indexState + 1
        errorStatus = nf_get_var_double(ncFileID,i,readData)
        CALL check_err(errorStatus)

        call check_dim_order(NDIMSOFVARS, varDimIDs(i,:), dimName, NDIMINID)

        if (ALL(NDIMINID .eq. (/ 1,2,3 /))) then
          readStates(:,:,indexState) = readData(:,:,1)
        else
          write(*,*) "Error ", trim(varName(i))
          STOP "Wrong order of dimensions! Expected: lon, lat, time"
        end if

        ! expected states
        IF (varName(i) == 'gothr') THEN         ! 'gothr' is primary vegetation
          checkPrimary = .true.
          indexOfPrimary = indexState
        ELSE IF (varName(i) == 'gsecd') THEN    ! 'gsecd' is secondary vegetation
          checkSecondary = .true.
          indexOfSecondary = indexState
        ELSE IF (varName(i) == 'gcrop') THEN
          checkCrop = .true.
          indexOfCrop = indexState
        ELSE IF (varName(i) == 'gpast') THEN
          checkPastrure = .true.
          indexOfPasture = indexState
        END IF

      END IF
    END DO

    IF (indexState /= NSTATES) STOP 'Unexpected number of states!'
    IF (.NOT. checkPrimary) STOP 'Variable gothr (primary) not found'
    IF (.NOT. checkSecondary) STOP 'Variable gsecd (secondary) not found'
    IF (.NOT. checkCrop) STOP 'Variable gcrop not found'
    IF (.NOT. checkPastrure) STOP 'Variable gpast not found'

    ! close the input-file states
    errorStatus = nf_close(ncFileID)
    CALL check_err(errorStatus)

  END SUBROUTINE read_states_file

  ! --- calculate_transitions_from_deltaS -------------------------------------------------------
  ! calculates direct transitions from the delta calculated from the states:
  ! TR: "Only two cases are possible regarding deltaS (if we don't consider the easy cases with one or all of the three entries of deltaS equals zero):
  !      Either 2 vegetation groups increase in area and the other group is shrinking (one transition from the shrinking to each of the two extending,
  !      the other transitions are zero) or only one vegetation group is increasing and the other two are shrinking (one transition from each of the
  !      two shrinking to the extending, the other transitions are zero). So, it's obvious how to construct Tabs,n(t -> t+1) from deltaS(t -> t+1)."
  ! --> problem as always: the numerical accuracy of the data
  SUBROUTINE calculate_transitions_from_deltaS(absTrans, deltaS, netAbsTrans)
  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: absTrans(:,:,:,:)
  REAL(dp), INTENT(INOUT) :: deltaS(:,:,:)
  REAL(dp), INTENT(OUT) :: netAbsTrans(:,:,:,:)

  INTEGER :: countTransStateIssues
  INTEGER :: countStateStateIssues
  INTEGER :: lon, lat, i, j

  LOGICAL :: setDeltaSToZero, setTransToZero
  LOGICAL :: isIncreasing(NTRANSSTATES)
  LOGICAL :: isDecreasing(NTRANSSTATES)

  countTransStateIssues = 0
  countStateStateIssues = 0
  DO lon = 1,nLon
    DO lat = 1,nLat

      ! Identify increasing and decreasing states
      isIncreasing = .FALSE.
      isDecreasing = .FALSE.
      DO i = 1,NTRANSSTATES
        IF (deltaS(lon,lat,i) .GT. 0.0) THEN
          isIncreasing(i) = .TRUE.
        ELSEIF (deltaS(lon,lat,i) .LT. 0.0) THEN
          isDecreasing(i) = .TRUE.
        ENDIF
      ENDDO

      ! Check and correct issues (here: numerical issues)
      setDeltaSToZero = .FALSE.
      setTransToZero = .FALSE.
      IF (COUNT(isIncreasing) .GT. 2) THEN
        CALL check_expected_accuracy_problems(lon, lat, absTrans(lon,lat,:,:), deltaS(lon,lat,:), &
            & 'ALL states are increasing', setDeltaSToZero, setTransToZero)
        isIncreasing = .FALSE.
        countTransStateIssues = countTransStateIssues + 1
      ELSEIF (COUNT(isDecreasing) .GT. 2) THEN
        CALL check_expected_accuracy_problems(lon, lat, absTrans(lon,lat,:,:), deltaS(lon,lat,:), &
            & 'ALL states are decreasing', setDeltaSToZero, setTransToZero)
        isDecreasing = .FALSE.
        countTransStateIssues = countTransStateIssues + 1
      ELSEIF ((COUNT(isDecreasing) .GE. 1) .AND. (COUNT(isIncreasing) .EQ. 0)) THEN
        CALL check_expected_accuracy_problems(lon, lat, absTrans(lon,lat,:,:), deltaS(lon,lat,:), &
            & 'ONE/TWO states decreasing, but NO state increases', setDeltaSToZero, setTransToZero)
        isDecreasing = .FALSE.
        countTransStateIssues = countTransStateIssues + 1
      ELSEIF ((COUNT(isIncreasing) .GE. 1) .AND. (COUNT(isDecreasing) .EQ. 0)) THEN
        CALL check_expected_accuracy_problems(lon, lat, absTrans(lon,lat,:,:), deltaS(lon,lat,:), &
            & 'ONE/TWO states increasing, but NO state decreases', setDeltaSToZero, setTransToZero)
        isIncreasing = .FALSE.
        countTransStateIssues = countTransStateIssues + 1
      ENDIF

      IF (setDeltaSToZero) deltaS(lon,lat,:) = 0.
      IF (setTransToZero) absTrans(lon,lat,:,:) = 0.

      ! Construct net abs transitions from deltaS
      IF ((COUNT(isDecreasing) .EQ. 2) .AND. (COUNT(isIncreasing) .EQ. 1)) THEN
        ! Two decreasing: thus using the change in the decreasing states
        DO i = 1,NTRANSSTATES
          DO j = 1,NTRANSSTATES
            IF (isDecreasing(i) .AND. isIncreasing(j)) THEN
              netAbsTrans(lon,lat,i,j) = ABS(deltaS(lon,lat,i))
            ENDIF
          ENDDO
        ENDDO
      ELSEIF ((COUNT(isIncreasing) .EQ. 2) .AND. (COUNT(isDecreasing) .EQ. 1)) THEN
        ! Two increasing: thus using the change in the increasing states
        DO i = 1,NTRANSSTATES
          DO j = 1,NTRANSSTATES
            IF (isDecreasing(i) .AND.  isIncreasing(j)) THEN
              netAbsTrans(lon,lat,i,j) = ABS(deltaS(lon,lat,j))
            ENDIF
          ENDDO
        ENDDO
      ELSEIF ((COUNT(isIncreasing) .EQ. 1) .AND. (COUNT(isDecreasing) .EQ. 1)) THEN
        ! Trivial case: transition from decreasing to increasing set to delta of one of them
        DO i = 1,NTRANSSTATES
          DO j = 1,NTRANSSTATES
            IF (isDecreasing(i) .AND.  isIncreasing(j)) THEN
              netAbsTrans(lon,lat,i,j) = ABS(deltaS(lon,lat,j))
              IF( ABS(deltaS(lon,lat,j)) .NE. ABS(deltaS(lon,lat,i))) THEN
                countStateStateIssues = countStateStateIssues + 1
                IF( ABS(deltaS(lon,lat,j) + deltaS(lon,lat,i)) .GT. EXPECTED_ACCURACY) THEN
                  WRITE(*,*) 'lon, lat, delta', lon, lat, deltaS(lon,lat,:)
                  STOP 'Difference in delta states higher than EXPECTED_ACCURACY!'
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ELSEIF ((COUNT(isIncreasing) .EQ. 0) .AND. (COUNT(isDecreasing) .EQ. 0)) THEN
        ! Trivial case: Nothing to do, netAbsTrans = 0.
      ELSE
        STOP 'ERROR: all cases should have been treated above: unexpected case, please check!'
      ENDIF

    ENDDO !lat
  ENDDO !lon

  WRITE(*,*) '### Counted transition vs state issues while calculating net transitions: ', countTransStateIssues
  WRITE(*,*) '### Counted state vs state issues while calculating net transitions: ', countStateStateIssues

  END SUBROUTINE calculate_transitions_from_deltaS


  ! --- check_expected_accuracy_problems ---------------------------------------------------
  ! According to Louis the accuracy of the LUH data is only about ~1.e-5, thus some accuracy problems are expected
  SUBROUTINE check_expected_accuracy_problems(lon, lat, &
      & absTrans, deltaS, checkedCase, setDeltaSToZero, setTransToZero)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lon, lat
  REAL(dp), INTENT(IN) :: absTrans(:,:)
  REAL(dp), INTENT(IN) :: deltaS(:)
  LOGICAL, INTENT(OUT) :: setDeltaSToZero, setTransToZero

  CHARACTER(LEN=*),INTENT(IN) :: checkedCase

  setDeltaSToZero = .FALSE.
  setTransToZero = .FALSE.

  ! check known cases:
  !- all transitions zero but delta smaller than EXPECTED_ACCURACY
  IF (ALL(absTrans(:,:) .EQ. 0.) .AND. ALL(ABS(deltaS(:)) .LT. EXPECTED_ACCURACY)) THEN
    setDeltaSToZero = .TRUE.
    IF (INFO) THEN
      WRITE(*,*) 'WARNING: Set deltaS to zero '&
          & // '-- all transitions were zero but '// trim(checkedCase) // ' [less than EXPECTED_ACCURACY] for (lon,lat) ', lon, lat
    ENDIF
  !- all transitions and deltaS smaller than EXPECTED_ACCURACY
  ELSEIF (ALL(ABS(absTrans(:,:)) .LT. EXPECTED_ACCURACY) &
      & .AND. ALL(deltaS(:) .LT. EXPECTED_ACCURACY)) THEN
    setDeltaSToZero = .TRUE.
    setTransToZero = .TRUE.
    IF (INFO) THEN
      WRITE(*,*) 'WARNING: Set deltaS and transitions to zero '&
          & // '-- all transitions and detlaS less than EXPECTED_ACCURACY and '// trim(checkedCase) // ' for (lon,lat) ', lon, lat
    ENDIF
  !- deltaS smaller than EXPECTED_ACCURACY
  ELSEIF (ALL(deltaS(:) .LT. EXPECTED_ACCURACY)) THEN
    setDeltaSToZero = .TRUE.
    IF (INFO) THEN
      WRITE(*,*) 'WARNING: Set deltaS to zero '&
          & // '-- detlaS less than EXPECTED_ACCURACY and '// trim(checkedCase) // ' for (lon,lat) ', lon, lat
    ENDIF
  ELSE
    WRITE(*,*) 'ERROR -'// trim(checkedCase) // '- (lon,lat)' , lon, lat
    DO i = 1,NTRANSSTATES
      WRITE(*,*) '*** delta: ', i, deltaS(i)
      WRITE(*,*) '*** transitions to this state: ', absTrans(:,i)
      WRITE(*,*) '*** transitions from this state: ', absTrans(i,:)
      WRITE(*,*) deltaS(i) - sum(absTrans(:,i)) + sum(absTrans(i,:))
    ENDDO
    STOP
  ENDIF

  END SUBROUTINE check_expected_accuracy_problems


  ! --- check_consistency_delta_and_transitions ---------------------------------------------------
  ! Required: delta = sum of transitions = incoming - outgoing => delta - incoming + outgoing = 0
  !
  !TODO maybe adapt transitions here and remove in the other function?
  SUBROUTINE check_consistency_delta_and_transitions(deltaS, trans, deltaTrans, wasModified)
  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: trans(:,:,:,:)
  REAL(dp), INTENT(IN) :: deltaS(:,:,:), deltaTrans(:,:,:,:)
  LOGICAL, INTENT(IN) :: wasModified

  REAL(dp) :: difference
  INTEGER :: lon, lat, countSeriousIssues, countSmallRelativeErrors, countSmallDiffs
  INTEGER :: countReplacedInconsistentUnmodifiedData
  LOGICAL :: hasSeriousIssue, hasSmallAbsDiff, hasSmallRelativeError, unmodifedDataInconsistent
  CHARACTER(LEN=128) :: checkedCase

  IF (wasModified) THEN
    checkedCase = 'final scaled transitions and scaled states'
  ELSE
    checkedCase = 'unscaled not modified states and transitions'
  ENDIF

  countSeriousIssues = 0
  countSmallRelativeErrors = 0
  countSmallDiffs = 0
  countReplacedInconsistentUnmodifiedData = 0
  DO lon = 1,nLon
    DO lat = 1,nLat
      hasSeriousIssue = .FALSE.
      hasSmallRelativeError = .FALSE.
      hasSmallAbsDiff = .FALSE.
      unmodifedDataInconsistent = .FALSE.
      DO i = 1,NTRANSSTATES
        difference = ABS(deltaS(lon,lat,i) + SUM(trans(lon,lat,i,:)) - SUM(trans(lon,lat,:,i)))
        !Two tests: one for the absolute error
        IF(difference .GT. RESULTING_ACCURACY_ERRORS) THEN
          ! and the second for the relative error
          IF((deltaS(lon,lat,i) .GT. 0.) .AND. (difference / deltaS(lon,lat,i) .LT. SMALL_RELATIVE_ERROR)) THEN
            hasSmallRelativeError = .TRUE.
          ELSE
            IF(.NOT. wasModified) THEN
              ! i.e. inconsistency in unmodified transitions and delta from unscaled states => inconsistency already in the data!
              unmodifedDataInconsistent = .TRUE.
              EXIT
            ELSE
              hasSeriousIssue = .TRUE.
              EXIT
            ENDIF
          ENDIF
        ELSEIF(difference .GT. 0.) THEN
          hasSmallAbsDiff = .TRUE.
        ENDIF
      ENDDO

      IF(hasSeriousIssue) THEN
        countSeriousIssues = countSeriousIssues + 1
        WRITE(*,*) '######################################################################################################'
        WRITE(*,*) 'WARNING: Found serious imbalance between delta and final transitions - Replaced with trans from delta!'
        CALL show_values_for_lon_lat(lon,lat)
        ! Replace by transitions derived from delta!
        trans(lon,lat,:,:) = deltaTrans(lon,lat,:,:)
      ELSEIF(unmodifedDataInconsistent) THEN
        countReplacedInconsistentUnmodifiedData = countReplacedInconsistentUnmodifiedData + 1
        IF(INFO)THEN
          WRITE(*,*) '######################################################################'
          WRITE(*,*) 'WARNING: inconsistent unmodified data:'
          CALL show_values_for_lon_lat(lon,lat)
        ENDIF
        ! Replace by transitions derived from delta!
        trans(lon,lat,:,:) = deltaTrans(lon,lat,:,:)
      ELSEIF(hasSmallRelativeError) THEN
        countSmallRelativeErrors = countSmallRelativeErrors + 1
      ELSEIF(hasSmallAbsDiff) THEN
        countSmallDiffs = countSmallDiffs + 1
      ENDIF
    ENDDO
  ENDDO

  IF(countReplacedInconsistentUnmodifiedData .GT. 0) THEN
    WRITE(*,*) '### WARNING: Serious imbalances between unmodified delta and transitions! Replaced with trans from delta ', &
        & countReplacedInconsistentUnmodifiedData
  ENDIF
  IF(countSeriousIssues .GT. 0) THEN
    WRITE(*,*) '### WARNING: Serious imbalances between delta and final transitions! Replaced with trans from delta ', &
        & countSeriousIssues
  ENDIF
  IF(countSmallRelativeErrors .GT. 0) THEN
    WRITE(*,*) '### Counted imbalances between delta and transitions > expected accuracy; but with a small relative error: ', &
        & countSmallRelativeErrors, trim(checkedCase)
  ENDIF
  IF(countSmallDiffs .GT. 0) THEN
    WRITE(*,*) '### Counted imbalances between delta and transitions < expected accuracy: ', countSmallDiffs, trim(checkedCase)
  ENDIF

  END SUBROUTINE check_consistency_delta_and_transitions


  ! --- check_delta_of_states --------------------------------------------------------------------------------------------
  ! ASSERTION: changes in states need to equal out (apart from accuracy)
  SUBROUTINE check_delta_of_states(deltaS, checkedCase)
  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: deltaS(:,:,:)

  CHARACTER(LEN=*),INTENT(IN) :: checkedCase
  INTEGER :: lon, lat, countSmallDiffs

  IF(ANY(ABS(SUM(deltaS, DIM=3)) .GT. EXPECTED_ACCURACY)) THEN
    DO lon = 1,nLon
      DO lat = 1,nLat
        IF(SUM(deltaS(lon, lat,:)) .GT. EXPECTED_ACCURACY) THEN
          WRITE(*,*) 'Found imbalanced changes in states larger than EXPECTED_ACCURACY for lon,lat,checkedCase ', &
              & lon, lat, checkedCase
          WRITE(*,*) 'deltaS(lon,lat,:)', lon, lat, deltaS(lon,lat,:)
          STOP
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  countSmallDiffs = COUNT(ABS(SUM(deltaS, DIM=3)) .GT. 0.)
  IF(countSmallDiffs .GT. 0) THEN
    WRITE(*,*) '### Counted delta state sums which were larger than zero but smaller than expected accuracy for case ', &
        & countSmallDiffs, checkedCase
  ENDIF

  END SUBROUTINE check_delta_of_states


  ! --- check_sum_of_transitions ----------------------------------------------------------------------------------
  ! ASSERTION: the sum of outgoing transitions may not exceed 1
  SUBROUTINE check_sum_of_transitions(transitions, stateThisYear, stateNextYear, checkedCase)
  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: transitions(:,:,:,:)
  REAL(dp), INTENT(IN) ::  stateThisYear(:,:,:), stateNextYear(:,:,:)
  CHARACTER(LEN=*),INTENT(IN) :: checkedCase

  INTEGER :: lon, lat, i, j, countSmallExceedance, countExceedanceSmallStates
  REAL(dp) :: thisSum

  countSmallExceedance = 0
  countExceedanceSmallStates = 0
  IF(ANY(SUM(transitions(:,:,:,:), DIM=4) .GT. 1.0)) THEN
    DO lon = 1,nLon
      DO lat = 1,nLat
        DO i = 1,NTRANSSTATES
          IF(SUM(transitions(lon,lat,i,:)) .GT. 1.0) THEN
            thisSum = SUM(transitions(lon,lat,i,:))
            IF(thisSum .LT. 1.0 + SMALL_RELATIVE_ERROR) THEN
              countSmallExceedance = countSmallExceedance + 1
              DO j = 1,NTRANSSTATES
                IF(i .NE. j) THEN
                  transitions(lon,lat,i,j) = transitions(lon,lat,i,j) - ( (thisSum-1.0) * transitions(lon,lat,i,j) / thisSum)
                ENDIF
              ENDDO
            ELSEIF((stateThisYear(lon,lat,i) .LT. EXPECTED_ACCURACY) .OR. (stateNextYear(lon,lat,i) .LT. EXPECTED_ACCURACY)) THEN
              countExceedanceSmallStates = countExceedanceSmallStates + 1
              DO j = 1,NTRANSSTATES
                IF(i .NE. j) THEN
                  transitions(lon,lat,i,j) = transitions(lon,lat,i,j) - ( (thisSum-1.0) * transitions(lon,lat,i,j) / thisSum)
                ENDIF
              ENDDO
            ELSE
              WRITE(*,*) 'ERROR: Found sum of outgoing transitions exceeding 1. + SMALL_RELATIVE_ERROR - for case ', checkedCase
              WRITE(*,*) 'i, sum', i, thisSum
              CALL show_values_for_lon_lat(lon,lat)
              STOP
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF(countSmallExceedance .GT. 0) THEN
    WRITE(*,*) '### Corrected sum of outgoing transitions exceeding one (but < SMALL_RELATIVE_ERROR) -- reduced proportionally ', &
        & countSmallExceedance, checkedCase
  ENDIF

  IF(countExceedanceSmallStates .GT. 0) THEN
    WRITE(*,*) '### Corrected sum of outgoing transitions exceeding one (state < expected accuracy) -- reduced proportionally ', &
        & countExceedanceSmallStates, checkedCase
  ENDIF

  END SUBROUTINE check_sum_of_transitions

  ! --- check_transitions_for_valid_range ---------------------------------------------------------------------------
  ! ASSERTION: transitions should be larger than zero and smaller than 1 (appart from accuracy)
  SUBROUTINE check_transitions_for_valid_range(trans, state, checkedCase)
  IMPLICIT NONE

  REAL(dp), INTENT(IN) :: state(:,:,:)
  CHARACTER(LEN=*),INTENT(IN) :: checkedCase

  REAL(dp), INTENT(INOUT) :: trans(:,:,:,:)

  INTEGER :: lon, lat, i, j, countExceedingButSmallStates, countSmallNegTrans, countSmallDiffTo1

  IF(ANY(trans .LT. -EXPECTED_ACCURACY)) THEN
    DO lon = 1,nLon
      DO lat = 1,nLat
        IF(ANY(trans(lon,lat,:,:) .LT. -EXPECTED_ACCURACY)) THEN
          WRITE(*,*) 'ERROR: Found transitions more negative than EXPECTED_ACCURACY ', checkedCase
          CALL show_values_for_lon_lat(lon,lat)
          STOP
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  countExceedingButSmallStates = 0
  IF(ANY(trans .GT. 1._dp + SMALL_RELATIVE_ERROR)) THEN
    DO lon = 1,nLon
      DO lat = 1,nLat
        DO i = 1,NTRANSSTATES
          DO j = 1,NTRANSSTATES
            IF(trans(lon,lat,i,j) .GT. 1._dp + SMALL_RELATIVE_ERROR) THEN
              ! test for the case that the associated state is smaller than EXPECTED_ACCURACY
              IF(state(lon,lat,i) .LT. EXPECTED_ACCURACY)THEN
                countExceedingButSmallStates = countExceedingButSmallStates + 1
                trans(lon,lat,i,j) = 0.0 !TODO: setting to one would require additional test for other transitions from this state!
              ELSE
                WRITE(*,*) 'ERROR: Found transitions larger than 1.0 (+ SMALL_RELATIVE_ERROR) ', checkedCase
                CALL show_values_for_lon_lat(lon,lat)
                STOP
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF(countExceedingButSmallStates .GT. 0) THEN
    WRITE(*,*) '### Corrected transition larger than 1, with states smaller than EXPECTED_ACCURACY to 0.0: ', &
        & countExceedingButSmallStates
  ENDIF

  countSmallDiffTo1 = COUNT(trans .GT. 1.)
  IF(countSmallDiffTo1 .GT. 0) THEN
    FORALL(lon=1:nLon,lat=1:nLat,i=1:NTRANSSTATES,j=1:NTRANSSTATES, trans(lon,lat,i,j) .GT. 1.)
      trans(lon,lat,i,j) = 1.
    END FORALL
    WRITE(*,*) '### Corrected transition larger than 1, but smaller than 1 + SMALL_RELATIVE_ERROR to 1.0: ', &
      & countSmallDiffTo1, checkedCase
  ENDIF

  countSmallNegTrans = COUNT(trans .LT. 0.)
  IF(countSmallNegTrans .GT. 0) THEN
    FORALL(lon=1:nLon,lat=1:nLat,i=1:NTRANSSTATES,j=1:NTRANSSTATES, trans(lon,lat,i,j) .LT. 0.)
      trans(lon,lat,i,j) = 0.
    END FORALL
    WRITE(*,*) '### Corrected negative transition which were smaller than expected accuracy to zero: ', &
        & countSmallNegTrans, checkedCase
  ENDIF

  END SUBROUTINE check_transitions_for_valid_range


  ! --- check_natural_veg_for_numerical_issues ------------------------------------------------
  SUBROUTINE check_natural_veg_for_numerical_issues(states, checkedCase)
  IMPLICIT NONE

  CHARACTER(LEN=*),INTENT(IN) :: checkedCase
  REAL(dp), DIMENSION(:,:,:), INTENT(INOUT) :: states

  INTEGER :: lon, lat, countNumericalIssues

  IF(ANY(states(:,:,1) .GT. 1._dp + EXPECTED_ACCURACY .AND. states(:,:,1) .NE. MISS_VALUE)) THEN
    DO lon = 1,nLon
      DO lat = 1,nLat
        IF(states(lon,lat,1) .GT. 1._dp + EXPECTED_ACCURACY .AND. states(lon,lat,1) .NE. MISS_VALUE) THEN
          WRITE(*,*) 'ERROR: Found natural area fraction larger than 1.0 + EXPECTED_ACCURACY'
          CALL show_values_for_lon_lat(lon,lat)
          STOP
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  countNumericalIssues = COUNT(states(:,:,1) .GT. 1._dp .AND. states(:,:,1) .NE. MISS_VALUE)
  FORALL(lon=1:nLon, lat=1:nLat, states(lon,lat,1) .GT. 1._dp .AND. states(lon,lat,1) .NE. MISS_VALUE)
    states(lon,lat,1) = 1._dp
  END FORALL

  IF(countNumericalIssues .GT. 0) THEN
    WRITE(*,*) '### Corrected nat area fracts larger than 1.0 (but < 1.0 + expected accuracy) to 1.0: ', &
      countNumericalIssues, checkedCase
  ENDIF

  END SUBROUTINE check_natural_veg_for_numerical_issues


  ! --- constrain_with_scaled_natural_vegetation ------------------------------------------------
  ! boundary condition: transitions involving natural vegetation are constrained by the available area
  ! (the veg_ratio_max scaling might reduce natural vegetation strongly)
  ! -> i.e. it is not possible to move more natural vegetation in a year than available in the scaled state
  ! Here: net transitions need to be conserved, here: second full cycles are conserved and only then bilateral cycles
  SUBROUTINE constrain_with_scaled_natural_vegetation(sortedScaledSThisY, &
        & trans_netScaledS, trans_ScaledS, absTrans_circ_fullCycles, absTrans_circ_bilateral, trans_bounded_circ)
  IMPLICIT NONE

  REAL(dp), DIMENSION(:,:,:), INTENT(IN) :: sortedScaledSThisY
  REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT) :: trans_netScaledS, trans_ScaledS, absTrans_circ_fullCycles, &
    absTrans_circ_bilateral, trans_bounded_circ

  INTEGER :: lon, lat
  INTEGER :: countReducedNetSmallTrans, countReducedNetSmallDiff
  INTEGER :: countRemovedNatBilateralAndReducedFullCycle, countReducedNatBilateral
  REAL(dp) :: remainingNatVeg, sumOfBilateralTransitions, netTransSum

  countReducedNetSmallTrans = 0.
  countReducedNetSmallDiff = 0.
  countRemovedNatBilateralAndReducedFullCycle = 0.
  countReducedNatBilateral = 0.

  DO lon = 1,nLon
    DO lat = 1,nLat

      remainingNatVeg = sortedScaledSThisY(lon,lat,1)
      IF(SUM(trans_ScaledS(lon,lat,1,:)) .GT. remainingNatVeg) THEN
        ! TODO: reconsider: constrain to available natural vegetation while trying to maintain net transitions,
        !       second full cycles and finally bilateral transitions involving nat vegetation

        !ASSERTION: since state delta was used to construct the net transitions,
        !           net transitions should not exceed the state by more than EXPECTED_ACCURACY
        IF(SUM(trans_netScaledS(lon,lat,1,:)) .GT. remainingNatVeg) THEN
          IF(SUM(trans_netScaledS(lon,lat,1,:)) .LT. EXPECTED_ACCURACY) THEN
            countReducedNetSmallTrans = countReducedNetSmallTrans + 1
            trans_netScaledS(lon,lat,1,:) = 0.
          ELSEIF(ABS(SUM(trans_netScaledS(lon,lat,1,:)) - remainingNatVeg) .LT. EXPECTED_ACCURACY) THEN
            !try to maintain net transitions
            !-> anticipated problems: 2 increasing 1 decreasing or 1 increasing and 1 decreasing
            !   -> delta decreasing does not match
            !   (since the transitions were constructed from the increasing deltas)
            !=> thus outgoing transitions from nat have to be constrained proportionally
            IF(trans_netScaledS(lon,lat,1,2) .GT. 0. .OR. trans_netScaledS(lon,lat,1,3) .GT. 0.) THEN
              countReducedNetSmallDiff = countReducedNetSmallDiff + 1
              trans_netScaledS(lon,lat,1,2) = remainingNatVeg * trans_netScaledS(lon,lat,1,2)/remainingNatVeg
              trans_netScaledS(lon,lat,1,3) = remainingNatVeg * trans_netScaledS(lon,lat,1,3)/remainingNatVeg
            ELSE
              WRITE(*,*) 'ERROR: Not anticipated problem for net transitions exceeding delta of scaled nat veg, please check!'
              CALL show_values_for_lon_lat(lon,lat)
              STOP
            ENDIF
          ELSE
            WRITE(*,*) 'ERROR: Found net transitions exceeding delta of scaled natural vegetation'
            CALL show_values_for_lon_lat(lon,lat)
            STOP
          ENDIF
        ENDIF
        remainingNatVeg = remainingNatVeg - SUM(trans_netScaledS(lon,lat,1,:))

        ! Reduce circular transitions until the boundary condition is met
        ! First reduce bilateral cycles, then full cycles
        ! Note: by construction there is only one full cycle per cell and year
        IF(SUM(absTrans_circ_fullCycles(lon,lat,1,:)) .GE. remainingNatVeg) THEN
          countRemovedNatBilateralAndReducedFullCycle = countRemovedNatBilateralAndReducedFullCycle + 1
          ! In this case none of the bilateral transitions involving natural vegetation can be realised
          absTrans_circ_bilateral(lon,lat,1,2) = 0.
          absTrans_circ_bilateral(lon,lat,2,1) = 0.
          absTrans_circ_bilateral(lon,lat,1,3) = 0.
          absTrans_circ_bilateral(lon,lat,3,1) = 0.
          ! and the full cycle needs to be reduced
          WHERE(absTrans_circ_fullCycles(lon,lat,:,:) .NE. 0.)
            absTrans_circ_fullCycles(lon,lat,:,:) = &
                & MIN(SUM(absTrans_circ_fullCycles(lon,lat,:,:)), remainingNatVeg)
          ENDWHERE
        ELSEIF(absTrans_circ_bilateral(lon,lat,1,2) + absTrans_circ_bilateral(lon,lat,1,3) .GT. 0.) THEN
          countReducedNatBilateral = countReducedNatBilateral + 1
          ! reduce the bilateral transitions involving natural vegetation (i.e. changes between crop and pasture are not constrained)
          sumOfBilateralTransitions = absTrans_circ_bilateral(lon,lat,1,2) + absTrans_circ_bilateral(lon,lat,1,3)
          remainingNatVeg = remainingNatVeg - SUM(absTrans_circ_fullCycles(lon,lat,1,:)) - sumOfBilateralTransitions

          absTrans_circ_bilateral(lon,lat,1,2) = absTrans_circ_bilateral(lon,lat,1,2) + &
            & (remainingNatVeg * (absTrans_circ_bilateral(lon,lat,1,2) / sumOfBilateralTransitions))
          absTrans_circ_bilateral(lon,lat,2,1) = absTrans_circ_bilateral(lon,lat,2,1) + &
            & (remainingNatVeg * (absTrans_circ_bilateral(lon,lat,2,1) / sumOfBilateralTransitions))
          absTrans_circ_bilateral(lon,lat,1,3) = absTrans_circ_bilateral(lon,lat,1,3) + &
            & (remainingNatVeg * (absTrans_circ_bilateral(lon,lat,1,3) / sumOfBilateralTransitions))
          absTrans_circ_bilateral(lon,lat,3,1) = absTrans_circ_bilateral(lon,lat,3,1) + &
            & (remainingNatVeg * (absTrans_circ_bilateral(lon,lat,3,1) / sumOfBilateralTransitions))
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  IF(countReducedNatBilateral .GT. 0) THEN
    WRITE(*,*) '### Counted reduction of nat bilateral transitions due to limited nat veg area: ', countReducedNatBilateral
  ENDIF
  IF(countRemovedNatBilateralAndReducedFullCycle .GT. 0) THEN
    WRITE(*,*) '### Counted removed nat bilateral and reduced full cylces due to limited nat veg area: ', &
        & countRemovedNatBilateralAndReducedFullCycle
  ENDIF
  IF(countReducedNetSmallTrans .GT. 0) THEN
    WRITE(*,*) '### Counted reduction of net transitions (< EXPECTED_ACCURACY) due to limited nat veg area: ', &
        & countReducedNetSmallTrans
  ENDIF
  IF(countReducedNetSmallDiff .GT. 0) THEN
    WRITE(*,*) '### Counted reduction of net transitions due to limited nat veg area (diff < EXPECTED_ACCURACY): ', &
        & countReducedNetSmallDiff
  ENDIF


  trans_bounded_circ = absTrans_circ_fullCycles + absTrans_circ_bilateral
  CALL check_circularity_of_circular_transitions(trans_bounded_circ, 'constraint circular transitions')

  END SUBROUTINE constrain_with_scaled_natural_vegetation


  ! --- constrain_detours_with_natural_vegetation --------------------------------------------------
  ! boundary conditions: ignore detours in cells which were scaled but are at the limit after scaling
  !                      + detours can only be fulfilled, if enough natural vegetation is available
  SUBROUTINE constrain_detours_with_natural_vegetation(sortedScaledSThisY, sortedUnscaledSThisY, &
      & trans_deltaScaledS, trans_scaled_diffNetDelta)
  IMPLICIT NONE

  REAL(dp), DIMENSION(:,:,:), INTENT(IN) :: sortedScaledSThisY, sortedUnscaledSThisY
  REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT) :: trans_deltaScaledS, trans_scaled_diffNetDelta

  INTEGER :: lon, lat, countViolatingDetours, countDetoursInCellsWithRestrictedScaling

  countViolatingDetours = 0
  countDetoursInCellsWithRestrictedScaling = 0

  DO lon = 1,nLon
    DO lat = 1,nLat
      IF ((sortedScaledSThisY(lon,lat,1) .LE. MIN_NAT_FRACT_IN_SCALING + EXPECTED_ACCURACY &
            & .AND. sortedUnscaledSThisY(lon,lat,1) .GT. MIN_NAT_FRACT_IN_SCALING + EXPECTED_ACCURACY) &
          & .OR. (sortedScaledSNextY(lon,lat,1) .LE. MIN_NAT_FRACT_IN_SCALING + EXPECTED_ACCURACY &
            & .AND. sortedUnscaledSNextY(lon,lat,1) .GT. MIN_NAT_FRACT_IN_SCALING + EXPECTED_ACCURACY) &
          & .AND. ANY(trans_scaled_diffNetDelta(lon,lat,:,:) .GT. 0.)) THEN
        ! TODO: reconsider: ignore detours when the scaled natural vegetation in this or in the next year
        !       is smaller than MIN_NAT_FRACT_IN_SCALING + EXPECTED_ACCURACY, but was larger in the unscaled case
        !       i.e. do not use detours when the scaling was affected by the 'remaining nat fraction' boundary condition rule
        countDetoursInCellsWithRestrictedScaling = countDetoursInCellsWithRestrictedScaling + 1
        trans_scaled_diffNetDelta(lon,lat,:,:) = 0.
      ELSEIF(SUM(trans_deltaScaledS(lon,lat,1,:) + trans_scaled_diffNetDelta(lon,lat,1,:)) &
          & .GT. sortedScaledSThisY(lon,lat,1)) THEN
        ! boundary condition: a detour is only possible if it does not violate the amount of nat vegetation available after scaling
        countViolatingDetours = countViolatingDetours + 1
        !TODO: instead of setting them to zero they should maybe set to the minimum possible,
        !      -> this needs some thought -- i.e. are several detours possible?
        trans_scaled_diffNetDelta(lon,lat,:,:) = 0.
      ENDIF
    ENDDO
  ENDDO

  IF(countDetoursInCellsWithRestrictedScaling .GT. 0) THEN
    WRITE(*,*) '### Counted detours not accounted for because cells were scaled down to minimum nat. fraction: ', &
        & countDetoursInCellsWithRestrictedScaling
  ENDIF
  IF(countViolatingDetours .GT. 0) THEN
    WRITE(*,*) '### Counted detours not accounted for due to violation of available nat. fraction: ', countViolatingDetours
  ENDIF

  END SUBROUTINE constrain_detours_with_natural_vegetation


  ! --- remove_unscaled_where_zero --------------------------------------------------
  ! when delta is zero, there should not be any transitions (check for EXPECTED_ACCURACY)
  SUBROUTINE remove_unscaled_where_zero(deltaUnscaledS, absTrans_netUnscaledS, absTrans_UnscaledS)
  IMPLICIT NONE

  REAL(dp), DIMENSION(:,:,:), INTENT(IN) :: deltaUnscaledS
  REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT) :: absTrans_netUnscaledS, absTrans_UnscaledS

  INTEGER :: lon, lat, i, countRemovedtransitions

  countRemovedtransitions = 0.
  DO lon = 1,nLon
    DO lat = 1,nLat
      DO i = 1,NTRANSSTATES
        IF(deltaUnscaledS(lon,lat,i) .EQ. 0. &
            & .AND. (SUM(absTrans_netUnscaledS(lon,lat,i,:)) - SUM(absTrans_netUnscaledS(lon,lat,:,i)) .GT. EXPECTED_ACCURACY)) THEN
          IF(ALL(absTrans_netUnscaledS(lon,lat,i,:) .LT. EXPECTED_ACCURACY) &
                & .AND. ALL(absTrans_netUnscaledS(lon,lat,:,i) .LT. EXPECTED_ACCURACY)) THEN
            countRemovedtransitions = countRemovedtransitions + 1
            absTrans_UnscaledS(lon,lat,i,:) = absTrans_UnscaledS(lon,lat,i,:) - absTrans_netUnscaledS(lon,lat,i,:)
            absTrans_netUnscaledS(lon,lat,i,:) = 0.
            absTrans_UnscaledS(lon,lat,:,i) = absTrans_UnscaledS(lon,lat,:,i) - absTrans_netUnscaledS(lon,lat,:,i)
            absTrans_netUnscaledS(lon,lat,:,i) = 0.
          ELSE
            WRITE(*,*) 'ERROR: delta zero, but unscaled transitions larger than EXPECTED_ACCURACY'
            CALL show_values_for_lon_lat(lon,lat)
            STOP
          ENDIF
        ENDIF
      ENDDO
    ENDDO !lat
  ENDDO !lon

  IF(countRemovedtransitions .GT. 0) THEN
    WRITE(*,*) '### Counted removed net transitions where delta states were zero: ', countRemovedtransitions
  ENDIF

  END SUBROUTINE remove_unscaled_where_zero



  ! --- check_circularity_of_circular_transitions --------------------------------------------------
  SUBROUTINE check_circularity_of_circular_transitions(absTrans_circ, checkedCase)
  IMPLICIT NONE

  REAL(dp), INTENT(IN) :: absTrans_circ(:,:,:,:)
  CHARACTER(LEN=*),INTENT(IN) :: checkedCase

  INTEGER :: lon, lat, i

  DO lon = 1,nLon
    DO lat = 1,nLat
      ! For each state the incoming transitions must equal the outgoing transitions
      DO i = 1,NTRANSSTATES
        IF(SUM(absTrans_circ(lon,lat,i,:)) .NE. SUM(absTrans_circ(lon,lat,:,i))) THEN
          IF(SUM(absTrans_circ(lon,lat,i,:)) - SUM(absTrans_circ(lon,lat,:,i)) .LT. EXPECTED_ACCURACY) THEN
            IF (INFO) THEN
              WRITE(*,*) 'WARNING: '// trim(checkedCase) // ' not really circular, '&
                  & // 'but difference less than EXPECTED_ACCURACY -- (i, lon, lat) ', i, lon, lat
            ENDIF
          ELSE
            WRITE(*,*) 'ERROR - '// trim(checkedCase) // ' not circular -- (i, lon, lat)' , i, lon, lat
            WRITE(*,*) '*** incoming: ', absTrans_circ(lon,lat,i,:)
            WRITE(*,*) '*** outgoing: ', absTrans_circ(lon,lat,:,j)
            CALL show_values_for_lon_lat(lon,lat)
            STOP
          ENDIF
        ENDIF
      ENDDO
    ENDDO !lat
  ENDDO !lon

  END SUBROUTINE check_circularity_of_circular_transitions


  ! --- calculate_bilateral_and_remaining_cycling --------------------------------------------
  ! By construction there is only one full cycle per cell and year
  ! Here: first get bilateral ones, than check, if a full cylce is left
  SUBROUTINE calculate_bilateral_and_remaining_cycling(nLon, nLat, &
      & absTrans_UnscaledS, absTrans_circ_bilateral, absTrans_circ_fullCycles)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nLon, nLat
  REAL(dp), INTENT(IN) :: absTrans_UnscaledS(:,:,:,:)
  REAL(dp), INTENT(OUT) :: absTrans_circ_bilateral(:,:,:,:), absTrans_circ_fullCycles(:,:,:,:)

  REAL(dp) :: absTrans_netUnscaledS_via_reduction(nLon,nLat,NTRANSSTATES,NTRANSSTATES)
  REAL(dp) :: minTrans
  INTEGER :: lon, lat, i, j
  LOGICAL :: aStateHasMoreThanTwoOutgoingTransitions

  absTrans_netUnscaledS_via_reduction = absTrans_UnscaledS
  DO lon = 1,nLon
    DO lat = 1,nLat
      !get bilateral ones
      DO i = 1,NTRANSSTATES
        DO j = i+1,NTRANSSTATES
          minTrans = MIN(absTrans_netUnscaledS_via_reduction(lon,lat,i,j), absTrans_netUnscaledS_via_reduction(lon,lat,j,i))
          absTrans_circ_bilateral(lon,lat,i,j) = minTrans
          absTrans_circ_bilateral(lon,lat,j,i) = minTrans

          absTrans_netUnscaledS_via_reduction(lon,lat,i,j) = absTrans_netUnscaledS_via_reduction(lon,lat,i,j) - minTrans
          absTrans_netUnscaledS_via_reduction(lon,lat,j,i) = absTrans_netUnscaledS_via_reduction(lon,lat,j,i) - minTrans
        ENDDO
      ENDDO
      ! get full cycles
      aStateHasMoreThanTwoOutgoingTransitions = .FALSE.
      minTrans = 1.
      DO i = 1,NTRANSSTATES
        IF (COUNT(absTrans_netUnscaledS_via_reduction(lon,lat,i,:) .GT. 0.) .GT. 1) THEN
          aStateHasMoreThanTwoOutgoingTransitions = .TRUE.
          EXIT
        ELSE
          DO j = 1,NTRANSSTATES
            IF (absTrans_netUnscaledS_via_reduction(lon,lat,i,j) .GT. 0. &
                & .AND. absTrans_netUnscaledS_via_reduction(lon,lat,i,j) .LT. minTrans) THEN
              minTrans = absTrans_netUnscaledS_via_reduction(lon,lat,i,j)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      ! we are done if any state has more than two outgoing transitions after reducing by the bilateral ones
      ! or if less than three transitions are left -- so only if none has two and three are left in total:
      IF ((.NOT. aStateHasMoreThanTwoOutgoingTransitions) &
            & .AND. (COUNT(absTrans_netUnscaledS_via_reduction(lon,lat,:,:) .GT. 0.) .EQ. 3)) THEN
        ! the smallest transition left is the amount which is fully cylced
        DO i = 1,NTRANSSTATES
          DO j = 1,NTRANSSTATES
            IF(absTrans_netUnscaledS_via_reduction(lon,lat,i,j) .GT. 0.) THEN
              absTrans_circ_fullCycles(lon,lat,i,j) = minTrans
              absTrans_netUnscaledS_via_reduction(lon,lat,i,j) = absTrans_netUnscaledS_via_reduction(lon,lat,i,j) - minTrans
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE calculate_bilateral_and_remaining_cycling


  ! --- show_values_for_lon_lat --------------------------------------------
  SUBROUTINE show_values_for_lon_lat(lon, lat)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: lon, lat

  WRITE(*,*) '############################################################################################'
  WRITE(*,*) 'lon, lat', lon, lat
  WRITE(*,*) '--- relative to either scaled or unscaled states, depending on when called...'
  WRITE(*,*) sortedTrans(lon,lat,1,:)
  WRITE(*,*) sortedTrans(lon,lat,2,:)
  WRITE(*,*) sortedTrans(lon,lat,3,:)
  WRITE(*,*) '--- scaled'
  WRITE(*,*) trans_ScaledS(lon,lat,1,:)
  WRITE(*,*) trans_ScaledS(lon,lat,2,:)
  WRITE(*,*) trans_ScaledS(lon,lat,3,:)
  WRITE(*,*) '--- scaled net'
  WRITE(*,*) trans_netScaledS(lon,lat,1,:)
  WRITE(*,*) trans_netScaledS(lon,lat,2,:)
  WRITE(*,*) trans_netScaledS(lon,lat,3,:)
  WRITE(*,*) '--- scaled delta'
  WRITE(*,*) trans_deltaScaledS(lon,lat,1,:)
  WRITE(*,*) trans_deltaScaledS(lon,lat,2,:)
  WRITE(*,*) trans_deltaScaledS(lon,lat,3,:)
  WRITE(*,*) 'delta scaled', deltaScaledS(lon,lat,:)
  WRITE(*,*) 'scaled next y', sortedScaledSNextY(lon,lat,:)
  WRITE(*,*) 'scaled this y', sortedScaledSThisY(lon,lat,:)
  WRITE(*,*) '--- scaled diff'
  WRITE(*,*) trans_scaled_diffNetDelta(lon,lat,1,:)
  WRITE(*,*) trans_scaled_diffNetDelta(lon,lat,2,:)
  WRITE(*,*) trans_scaled_diffNetDelta(lon,lat,3,:)
  WRITE(*,*) '--- net un'
  WRITE(*,*) absTrans_netUnscaledS(lon,lat,1,:)
  WRITE(*,*) absTrans_netUnscaledS(lon,lat,2,:)
  WRITE(*,*) absTrans_netUnscaledS(lon,lat,3,:)
  WRITE(*,*) '--- delta un'
  WRITE(*,*) absTrans_deltaUnscaledS(lon,lat,1,:)
  WRITE(*,*) absTrans_deltaUnscaledS(lon,lat,2,:)
  WRITE(*,*) absTrans_deltaUnscaledS(lon,lat,3,:)
  WRITE(*,*) 'delta unscaled', deltaUnscaledS(lon,lat,:)
  WRITE(*,*) 'unscaled next y', sortedUnscaledSNextY(lon,lat,:)
  WRITE(*,*) 'unscaled this y', sortedUnscaledSThisY(lon,lat,:)
  WRITE(*,*) 'vrm', veg_ratio_max(lon,lat)
  WRITE(*,*) '--- circ'
  WRITE(*,*) absTrans_circ_via_reduction(lon,lat,1,:)
  WRITE(*,*) absTrans_circ_via_reduction(lon,lat,2,:)
  WRITE(*,*) absTrans_circ_via_reduction(lon,lat,3,:)
  WRITE(*,*) '--- org'
  WRITE(*,*) absTrans_UnscaledS(lon,lat,1,:)
  WRITE(*,*) absTrans_UnscaledS(lon,lat,2,:)
  WRITE(*,*) absTrans_UnscaledS(lon,lat,3,:)

  END SUBROUTINE show_values_for_lon_lat

  subroutine check_dim_order(expectedNrOfDims, varDimIDs, nameofDim, orderOfInputDim)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: expectedNrOfDims
    INTEGER, INTENT(IN), Dimension(expectedNrOfDims) :: varDimIDs
    CHARACTER(LEN=4), INTENT(IN), Dimension(expectedNrOfDims)  ::  nameofDim
    INTEGER, INTENT(OUT), Dimension(expectedNrOfDims) :: orderOfInputDim
    ! local variabel
    INTEGER :: I

    DO I = 1,expectedNrOfDims
      IF (nameofDim(varDimIDs(I)) == 'lon') THEN
        orderOfInputDim(1) = I
      ELSE IF (nameofDim(varDimIDs(I)) == 'lat') THEN
        orderOfInputDim(2) = I
      ELSE IF (nameofDim(varDimIDs(I)) == 'time') THEN
        orderOfInputDim(3) = I
      END IF
    END DO

  end SUBROUTINE check_dim_order

  ! --- prepare_input_file --------------------------------------------------------------------------------------------
  ! TODO: this is very clumsy -- should be refactored
  SUBROUTINE prepare_input_file(fileName, expectedNrVars, expectedNrOfDims, &
      & ncFileID, varName, varTyp, nVarAttr, dimLength, varDimIDs, dimName)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'

  CHARACTER(LEN=*), INTENT(IN) ::  fileName
  INTEGER, INTENT(IN) :: expectedNrVars, expectedNrOfDims
  INTEGER, INTENT(OUT) :: ncFileID
  CHARACTER(LEN=128), INTENT(OUT) :: varName(:)
  INTEGER, INTENT(OUT) :: nVarAttr(:), varTyp(:)
  INTEGER, INTENT(OUT) :: dimLength(:)
  CHARACTER(LEN=4), INTENT(OUT) :: dimName(:)
  INTEGER, INTENT(OUT) :: varDimIDs(:,:)

  INTEGER :: i, j, nDim, nVar, nAttr, unLimDimID, errorStatus

  INTEGER :: readDimIDs(expectedNrOfDims)
  INTEGER :: nVarDims(expectedNrVars)

  LOGICAL :: checkLonDim, checkLatDim

  ! open the file
  WRITE (*,*) 'open file: ',TRIM(fileName)
  errorStatus = nf_open(fileName,nf_nowrite,ncFileID)
  CALL check_err(errorStatus)
  ! check what is in the input-file
  errorStatus = nf_inq(ncFileID,nDim,nVar,nAttr,unLimDimID)
  CALL check_err(errorStatus)
  IF (INFO) WRITE (*,*) 'nDim,nVar,nAttr,unLimDimID',nDim,nVar,nAttr,unLimDimID

  IF (nDim .NE. expectedNrOfDims) STOP 'Unexpected number of dimensions'
  IF (nVar .NE. expectedNrVars) STOP 'Unexpected number of variables'

  ! get the dimension name and length
  DO i=1,nDim
    errorStatus = nf_inq_dim(ncFileID,i,dimName(i),dimLength(i))
    CALL check_err(errorStatus)
    IF (INFO) WRITE (*,*) 'i,dimName,dimLength',i,TRIM(dimName(i)),dimLength(i)
  END DO

  ! get variable names, types and shapes
  DO i=1,nVar
    errorStatus = nf_inq_var(ncFileID, i, varName(i), varTyp(i), nVarDims(i), readDimIDs, nVarAttr(i))
    CALL check_err(errorStatus)
    DO j=1,expectedNrOfDims
      varDimIDs(i,j)=readDimIDs(j)
    END DO

    IF (INFO) WRITE (*,*) 'i,varName,varTyp,nVarDims,varDimIDs(1),nVarAttr', &
          i,TRIM(varName(i)),varTyp(i),nVarDims(i),varDimIDs(i,1),nVarAttr(i)
  END DO

  ! check number of longitudes and latitudes
  checkLonDim = .false.
  checkLatDim = .false.
  DO i = 1,nDim
      IF (TRIM(dimName(i)) == 'lon') THEN
        IF (dimLength(i) == nLon) checkLonDim = .true.
      ELSEIF (TRIM(dimName(i)) == 'lat') THEN
        IF (dimLength(i) == nLat) checkLatDim = .true.
      END IF
  END DO

  IF (.NOT. checkLonDim) STOP 'Dimension lon not found or number of longitudes not correct'
  IF (.NOT. checkLatDim) STOP 'Dimension lat not found or number of latitudes not correct'

  END SUBROUTINE prepare_input_file


  ! --- check_err --------------------------------------------------------------------------------------------
  SUBROUTINE check_err(iret,text)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER,INTENT(IN)          :: iret
  CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: text
  IF (iret /= NF_NOERR) THEN
    IF (present(text)) WRITE(*,*) 'check_err(): '//TRIM(text)
    WRITE(*,*) 'check_err(): netcdf says: '//TRIM(nf_strerror(iret))
    STOP
  END IF

  END SUBROUTINE check_err


  !COPY from TRs prog
   SUBROUTINE PUTATT(NVAR,I,NVARATT,NATTTYP,FATTNAM,NATTDIM,FATT,NCIDOUT,NVARID,NATT,RATT,INFO)
!  This routine adds the attributes to the output-file
   IMPLICIT NONE
   INCLUDE 'netcdf.inc'
   INTEGER, INTENT(IN)  ::  NVAR
   INTEGER, INTENT(IN)  ::  I
   INTEGER, INTENT(IN), DIMENSION(NVAR)  ::  NVARATT
   INTEGER, INTENT(IN), DIMENSION(NVAR,NVARATT(I))  :: NATTTYP
   CHARACTER(LEN=*), INTENT(IN), DIMENSION(NVAR,NVARATT(I))  ::  FATTNAM
   INTEGER, INTENT(IN), DIMENSION(NVAR,NVARATT(I)) ::  NATTDIM
   CHARACTER(LEN=*), INTENT(IN), DIMENSION(NVAR,NVARATT(I))  ::  FATT
   INTEGER, INTENT(IN)  ::  NCIDOUT, NVARID
   INTEGER, INTENT(IN), DIMENSION(NVAR,NVARATT(I))  ::  NATT
   REAL(dp), INTENT(IN), DIMENSION(NVAR,NVARATT(I)) ::  RATT
   LOGICAL, INTENT(IN)  ::  INFO
   INTEGER  ::  III,IRET

   DO III=1,NVARATT(I)
     IF (NATTTYP(I,III).EQ.2) THEN
       IF (INFO) WRITE (*,*) 'put_att:',TRIM(FATTNAM(I,III)),NATTDIM(I,III),TRIM(FATT(I,III))
       IRET = nf_put_att_text(NCIDOUT,NVARID,FATTNAM(I,III),NATTDIM(I,III),FATT(I,III))
       CALL check_err(IRET)
     ELSE IF (NATTTYP(I,III).EQ.4) THEN
        IF (INFO) WRITE (*,*) 'put_att:',TRIM(FATTNAM(I,III)),NATTTYP(I,III),NATTDIM(I,III),NATT(I,III)
        IRET = nf_put_att_int(NCIDOUT,NVARID,FATTNAM(I,III),NATTTYP(I,III),NATTDIM(I,III),NATT(I,III))
        CALL check_err(IRET)
     ELSE IF ((NATTTYP(I,III).EQ.5) .OR. (NATTTYP(I,III).EQ.6)) THEN
        IF (INFO) WRITE (*,*) 'put_att:',TRIM(FATTNAM(I,III)),NATTTYP(I,III),NATTDIM(I,III),RATT(I,III)
        IRET = nf_put_att_double(NCIDOUT,NVARID,FATTNAM(I,III),NATTTYP(I,III),NATTDIM(I,III),RATT(I,III))
        CALL check_err(IRET)
      ELSE
        WRITE (*,*) 'Type of attribute not prepared in the code, failed to write',I,III,NATTTYP(I,III)
      END IF
   END DO

   END SUBROUTINE PUTATT


END PROGRAM Adapting_transitions_after_scaling_of_states
