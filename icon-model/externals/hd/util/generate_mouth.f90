! generate_mouth.f90 - Generate coastal ocean points from a given ocean model land sea mask
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
!     ******** Programm that generates potential mouths as coastal ocean points from an 
!              ocean model (e.g. NEMO) sea mask (sea = 1) that neighbor land points and
!              (currently that neighbor missing values) (= 0) at the NSEW side.
!              Note that the oceangrid mouth masks usually define sea points as 1, 
!              land and missing values as Zero. Hence, missing values must be set to Zero.
!              The variable name of the mask must be 'mask'.
!
!     ******** Version 2.0 - June 2020
!              Programmierung und Entwicklung: Stefan Hagemann, HZG
!              Further development of generate_mouth_nemo, but now mouth points are
!              indicate by 1, with 0 elsewhere
!
!              Input file: ocean_mask.nc - Variablename: mask
!              It writes file: coast_oceangrid.nc
!
!     ******** Version 2.1 - Feb. 2022 - Stefan Hagemann, Hereon
!              Extended classification:
!                 1  Coastal ocean points
!                 0  other ocean points
!                 2  Land points
!
!     ******** Variablenliste:
!     ***  
!     ***       dnam = Name of Mask file: Set to ocean_mask.nc
!     ***       IQUE = Kommentarvariable ( 0 = Kein Kommentar )
!     ***         LU = Logical Unit der zu lesenden Datei
!     *** mask_mouth = Muendungsmaske on ocean grid
!     ***      mask_src%value = Sea mask NEMO
!     ***  
      use netcdf
      use mo_grid,         ONLY: domain, model, read_grid_info, read_coordinates
!
      DOUBLE PRECISION, PARAMETER :: xmiss=-9999. 
      REAL, PARAMETER :: zeps = 1.E-6 
      INTEGER, PARAMETER  :: IQUE =0     ! Print out some debug statements for IQUE=1
      INTEGER I
!     Ocean grid arrays

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_mouth
      REAL, DIMENSION(:,:), ALLOCATABLE :: fdat
      REAL XMAX, XMIN, YMAX, YMIN

      TYPE(domain)        :: gr_src      ! Source grid
      TYPE(model)         :: mask_src    ! Mask with inflow points and source coordinates
      CHARACTER (len=20) :: clon, clat, cdimlon, cdimlat
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, varid, dimids(2)
!  
      CHARACTER (len=128) :: dnam, ORDER, DNOUT
!
      dnam="ocean_mask.nc"
      DNOUT="coast_oceangrid.nc"
      LU=50
      clon = 'lon' ; clat = 'lat'
      cdimlon = 'x' ; cdimlat = 'y'
!	  
!     ******** READ oceangrid Mask points and coordinates
!!      CALL read_grid_info(dnam, gr_src, xmiss, clon, clat, cdimlon,cdimlat, idim)
      CALL read_grid_info(dnam, gr_src, xmiss, clon, clat, cdimlon,cdimlat)
      mask_src%nlon = gr_src%nlon
      mask_src%nlat = gr_src%nlat
      WRITE(*,*) "Oceangrid, e.g. NEMO:  gr_src%nlon =" , gr_src%nlon, "  gr_src%nlat = ", gr_src%nlat
!
!     *** Feld Dimensionierung
      ALLOCATE(mask_src%value(gr_src%nlon,gr_src%nlat))
      ALLOCATE(mask_src%xlon(gr_src%nlon,gr_src%nlat))
      ALLOCATE(mask_src%xlat(gr_src%nlon,gr_src%nlat))
      ALLOCATE(mask_mouth(gr_src%nlon,gr_src%nlat))
      ALLOCATE(fdat(gr_src%nlon,gr_src%nlat))
!
!     ******** READ arrays on source grid, e.g. HD: mask_src%value
      CALL read_coordinates(dnam, mask_src, clon, clat)
      XMIN = MINVAL(mask_src%xlon)
      XMAX = MAXVAL(mask_src%xlon)
      YMIN = MINVAL(mask_src%xlat)
      YMAX = MAXVAL(mask_src%xlat)
      WRITE(*,*) "Ocean-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "              Lat ", YMIN,' - ', YMAX 
!
      ierr = nf90_open(dnam, NF90_NOWRITE, ncid)
      ierr = nf90_inq_varid(ncid,'mask',varid)
      ierr = nf90_get_var(ncid,varid,fdat)
      mask_src%value = NINT(fdat)
      ierr = nf90_close(ncid)

      WRITE(*,*) 'Mask array on source grid: ', &
                MINVAL(mask_src%value, ABS(mask_src%value-xmiss).GT.zeps), ' - ', MAXVAL(mask_src%value)
      WRITE(*,*) '         SUM: ', SUM(mask_src%value, ABS(mask_src%value-xmiss).GT.zeps)
!
!     *** Note that the oceangrid mouth masks defines sea points as 1, land and missing values as Zero.
!
!     ******** Defining coastal ocean points as ONE (NOT Zero as previously done) in main N-S and W-E directions.
      mask_mouth(:,:)=0
      NMOU=0
      DO JB = 1, gr_src%nlat
      DO JL = 1, gr_src%nlon
      IF ( ABS(mask_src%value(JL,JB)-1).LE.0.1 ) THEN
        ICOAST=0
        IF (JB.GT.1 .AND. ABS(mask_src%value(JL,JB-1)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (JB.LT.gr_src%nlat .AND. ABS(mask_src%value(JL,JB+1)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (JL.GT.1 .AND. ABS(mask_src%value(JL-1,JB)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (JL.LT.gr_src%nlon .AND. ABS(mask_src%value(JL+1,JB)).LE.0.1 ) ICOAST = ICOAST + 1
        IF (ICOAST.GT.0.5) THEN
          mask_mouth(JL,JB) = 1
          NMOU = NMOU + 1
        ENDIF
      ELSE    ! Land box
        mask_mouth(JL,JB) = 2
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
      ierr = nf90_def_dim(ncid, cdimlon, gr_src%nlon, dimids(1))
      ierr = nf90_def_dim(ncid, cdimlat, gr_src%nlat, dimids(2))

      ierr = nf90_def_var(ncid, 'mask_coast_ocean', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 731)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','Coastal ocean points on ocean grid')
      ierr = nf90_put_att(ncid,varid,'standard_name','mask_coast_ocean')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_put_att(ncid,varid,'missing_value',-9999.)

      ierr = nf90_put_att(ncid, NF90_GLOBAL, 'Comment', "1 = Coastal Ocean Point, 0 = Other ocean points, 2 = land point")
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_mouth)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_east')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_src%xlon)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_north')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_src%xlat)
      ierr = nf90_close(ncid)
      WRITE(*,*) "Coastal ocean points were written: ", TRIM(DNOUT)
!
      DEALLOCATE(mask_mouth)
      DEALLOCATE(mask_src%value)
      DEALLOCATE(mask_src%xlon)
      DEALLOCATE(mask_src%xlat)
      DEALLOCATE(fdat)
!
!     ******** Programmende
  999 END
!
