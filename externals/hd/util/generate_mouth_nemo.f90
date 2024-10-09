! generate_mouth_nemo.f90 - Generates mask of coastal ocean points from a NEMO sea-land mask
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

      PROGRAM GEN_MOUTH
!
!     ******** Programm that generates potential mouths as coastal ocean points from NEMO sea mask ( = 1)
!     *** that neighbor land points (currently that neighbor missing values) (= 0) at the NSEW side.
!
!     ******** Version 1.0 - Jan. 2017
!              Programmierung und Entwicklung: Stefan Hagemann 
!              It writes file: river_mouths_NEMO.nc
!
!     ******** Version 1.1 - Oct. 2018 - Generalization
!              Input file: ocean_mask.nc - Variablename: mask
!              It writes file: river_mouths_oceangrid.nc
!
!     ******** Variablenliste:
!     ***  
!     ***   DNAM = Name of Mask file: Set to ocean_mask.nc
!     ***   IQUE = Kommentarvariable ( 0 = Kein Kommentar )
!     ***     LU = Logical Unit der zu lesenden Datei
!     *** MASK_NEMO = Muendungsmaske on NEMO Grid
!     ***   FMAS = Sea mask NEMO
!     ***  
      use netcdf
!
      PARAMETER(XMISS=-9999.) 
      INTEGER, PARAMETER  :: IQUE =0     ! Print out some debug statements for IQUE=1
      INTEGER NX, NY, IHDDIM, I
      INTEGER, DIMENSION(2) :: ICOORD
!     Ocean grid arrays
      REAL, DIMENSION(:,:), ALLOCATABLE :: XNEMO    ! Longitudes
      REAL, DIMENSION(:,:), ALLOCATABLE :: YNEMO    ! Latitudes
      REAL, DIMENSION(:,:), ALLOCATABLE :: FMAS
      REAL, DIMENSION(:), ALLOCATABLE :: XDUM   ! Longitudes in ocean source file if 1-D
      REAL, DIMENSION(:), ALLOCATABLE :: YDUM   ! Latitudes in ocean source file if 1-D

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_NEMO
      REAL XMAX, XMIN, YMAX, YMIN
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, dimid, varid, dimids(2), numdims
!  
      CHARACTER (len=128) :: DNAM, ORDER, DNOUT
!
      DNAM="ocean_mask.nc"
      DNOUT="river_mouths_oceangrid.nc"
      LU=50
!	  
!     ******** READ oceangrid Mask points and coordinates
      ierr = nf90_open(DNAM, NF90_NOWRITE, ncid)
      ierr = nf90_inq_dimid(ncid,'x',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NX)
      ierr = nf90_inq_dimid(ncid,'y',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NY)
      IF (NX.EQ.0 .OR. NY.EQ.0) THEN
        ierr = nf90_inq_dimid(ncid,'lon',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NX)
        ierr = nf90_inq_dimid(ncid,'lat',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NY)
      ENDIF
      ! Testing --> NX = 619, NY=523
      WRITE(*,*) "Oceangrid, e.g. NEMO:  NX =" , NX, "  NY = ", NY
!
!     *** Feld Dimensionierung
      ALLOCATE(XNEMO(NX,NY))
      ALLOCATE(YNEMO(NX,NY))
      ALLOCATE(FMAS(NX,NY))
      ALLOCATE(MASK_NEMO(NX,NY))
!
!     Reading arrays from NETCDF file
      ierr = nf90_inq_varid(ncid,'mask',varid)
      ierr = nf90_get_var(ncid,varid,FMAS)

      ierr = nf90_inq_varid(ncid,'lon',varid)
      ierr = nf90_inquire_variable(ncid, varid, ndims = numdims)
      IF (numdims.EQ.2) THEN
        ierr = nf90_get_var(ncid,varid,XNEMO)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,YNEMO)
      ELSE IF (numdims.EQ.1) THEN
        ALLOCATE(XDUM(NX))
        ALLOCATE(YDUM(NY))
        ierr = nf90_get_var(ncid,varid,XDUM)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,YDUM)
        DO I=1, NY
          XNEMO(:,I) = XDUM(:)
        ENDDO
        DO I=1, NX
          YNEMO(I,:) = YDUM(:)
        ENDDO
        DEALLOCATE(XDUM)
        DEALLOCATE(YDUM)
      ELSE
        STOP 'No. of coordinate dimension neither 1 nor 2 --> STOP!'
      ENDIF
      ierr = nf90_close(ncid)
      XMIN = MINVAL(XNEMO)
      XMAX = MAXVAL(XNEMO)
      YMIN = MINVAL(YNEMO)
      YMAX = MAXVAL(YNEMO)
      WRITE(*,*) "Ocean-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "              Lat ", YMIN,' - ', YMAX 
      WRITE(*,*) "Mask ranges from ", MINVAL(FMAS),' - ', MAXVAL(FMAS)

!
!     *** Note that the oceangrid mouth masks defines sea points as 1, land and missing values as Zero.
!
!     ******** Defining coastal ocean points as Zero.
      MASK_NEMO(:,:)=1
      NMOU=0
      DO JB = 1, NY
      DO JL = 1, NX
      IF ( ABS(FMAS(JL,JB)-1).LE.0.1 ) THEN
        ICOAST=0
        IF (JB.GT.1 .AND. ABS(FMAS(JL,JB-1)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (JB.LT.NY .AND. ABS(FMAS(JL,JB+1)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (JL.GT.1 .AND. ABS(FMAS(JL-1,JB)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (JL.LT.NX .AND. ABS(FMAS(JL+1,JB)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (ICOAST.GT.0.5) THEN
          MASK_NEMO(JL,JB) = 0
          NMOU = NMOU + 1
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      WRITE(*,*) NMOU, ' coastal ocean points in oceangrid region found'
      IF (NMOU.EQ.0) STOP 'No coastal ocean point found/generated --> ERROR'
!
!     ******** Schreiben des potential mouth array
!
!     *** WRITE NETCDF output
      ierr = nf90_create(DNOUT, NF90_CLOBBER, ncid)
      ierr = nf90_def_dim(ncid, 'x', NX, dimids(1))
      ierr = nf90_def_dim(ncid, 'y', NY, dimids(2))

      ierr = nf90_def_var(ncid, 'nmhd_msk', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 731)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','Coastal ocean points on ocean grid')
      ierr = nf90_put_att(ncid,varid,'standard_name','MASK_NEMO')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_put_att(ncid,varid,'missing_value',-9999.)
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, MASK_NEMO)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_FLOAT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_east')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, XNEMO)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_FLOAT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_north')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, YNEMO)
      ierr = nf90_close(ncid)
      WRITE(*,*) "Coastal ocean points were written: ", TRIM(DNOUT)
!
      DEALLOCATE(MASK_NEMO)
      DEALLOCATE(FMAS)
      DEALLOCATE(XNEMO)
      DEALLOCATE(YNEMO)
!
!     ******** Programmende
  999 END
!
