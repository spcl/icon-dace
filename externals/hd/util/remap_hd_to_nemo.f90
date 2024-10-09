! remap_hd_to_nemo.f90 - Program that generates a remapping (coupling) file from HD river mouths to NEMO ocean grid 
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

      PROGRAM remap_hd_to_nemo
!
!     ******** Programm that remaps a HD mask to the NEMO grid using nearest neighbor remapping.
!              As OASIS and cdos fail to produce a general and exact 1 to 1 NN remapping, the NN remapping
!              is taken from the self-programmed module mo_interpol
!
!     ******** Version 1.0 - May 2018
!              Programmierung und Entwicklung: Stefan Hagemann 
!              Input files: mask_hd.nc, nemo_grid.nc, Output file: mask_on_nemo.nc
!
!     ******** Version 1.1 - Oct. 2018 - Generalization
!              Input files: mask_hd.nc, ocean_grid.nc
!              Output file: mask_on_oceangrid.nc
!
!     ******** Variablenliste:
!     ***  
!     ***   DNHD = Name of HD Mask file = mask_hd.nc
!     *** DNNEMO = Name of file with NEMO grid information and coordinates = nemo_grid.nc
!     ***  DNOUT = Name of Output on NEMO grid = mask_on_nemo.nc
!     ***   IQUE = Kommentarvariable ( 0 = Kein Kommentar )
!     ***     LU = Logical Unit der zu lesenden Datei
!     *** MASK_NEMO = Muendungsmaske on NEMO Grid
!     ***   FMAS = Sea mask NEMO
!     ***  
      use netcdf
      use mo_interpol,              ONLY: remapnn
!
      REAL, PARAMETER  :: xmiss=-9999.
      INTEGER, PARAMETER  :: IQUE =0     ! Print out some debug statements for IQUE=1
      INTEGER nx, ny, IHDDIM, nlon, nlat
      INTEGER, DIMENSION(2) :: ICOORD
!
!     HD arrays
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xlon
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xlat
      REAL, DIMENSION(:,:), ALLOCATABLE :: fdata
!!      INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_hd
!     NEMO arrayS
      REAL, DIMENSION(:,:), ALLOCATABLE :: xnemo
      REAL, DIMENSION(:,:), ALLOCATABLE :: ynemo
      REAL, DIMENSION(:,:), ALLOCATABLE :: fnemo
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: mask_nemo
      REAL, DIMENSION(:), ALLOCATABLE :: xdum   ! Longitudes in ocean source file if 1-D
      REAL, DIMENSION(:), ALLOCATABLE :: ydum   ! Latitudes in ocean source file if 1-D
      REAL XMAX, XMIN, YMAX, YMIN
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, dimid, varid, dimids(2)
      CHARACTER (len=128) :: DNHD = "mask_hd.nc"
      CHARACTER (len=128) :: DNNEMO = "ocean_grid.nc"
      CHARACTER (len=128) :: DNOUT = "mask_on_oceangrid.nc"
      CHARACTER (len=128) :: ORDER
!
      LU=50
!	  
!     ******** READ HD Mask and coordinates
      ierr = nf90_open(DNHD, NF90_NOWRITE, ncid)
      ierr = nf90_inq_dimid(ncid,'lon',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=nlon)
      ierr = nf90_inq_dimid(ncid,'lat',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=nlat)

      WRITE(*,*) "HD:  nlon =" , nlon, "  nlat = ", nlat
!
!     *** Feld Dimensionierung
      ALLOCATE(xlon(nlon))
      ALLOCATE(xlat(nlat))
      ALLOCATE(fdata(nlon,nlat))
!!      ALLOCATE(mask_hd(nlon,nlat))
!
!     Reading HD arrays from NETCDF file
      ierr = nf90_inq_varid(ncid,'FDIR',varid)
      ierr = nf90_get_var(ncid,varid,fdata)
      ierr = nf90_inq_varid(ncid,'lon',varid)
      ierr = nf90_get_var(ncid,varid,xlon)
      ierr = nf90_inq_varid(ncid,'lat',varid)
      ierr = nf90_get_var(ncid,varid,xlat)
      ierr = nf90_close(ncid)
      XMIN = MINVAL(xlon)
      XMAX = MAXVAL(xlon)
      YMIN = MINVAL(xlat)
      YMAX = MAXVAL(xlat)
      WRITE(*,*) "HD-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "           Lat ", YMIN,' - ', YMAX 
!      mask_hd(:,:) = FLOOR(FDAT(:,:) + 0.001)
!	  
!     ******** READ NEMO grid info and coordinates
      ierr = nf90_open(DNNEMO, NF90_NOWRITE, ncid)
      ierr = nf90_inq_dimid(ncid,'x',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=nx)
      ierr = nf90_inq_dimid(ncid,'y',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=ny)
      IF (nx.EQ.0 .OR. ny.EQ.0) THEN
        ierr = nf90_inq_dimid(ncid,'lon',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=nx)
        ierr = nf90_inq_dimid(ncid,'lat',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=ny)
      ENDIF

      ! Testing --> nx = 619, ny=523
      WRITE(*,*) "NEMO:  nx =" , nx, "  ny = ", ny
!
!     *** Feld Dimensionierung
      ALLOCATE(xnemo(nx,ny))
      ALLOCATE(ynemo(nx,ny))
      ALLOCATE(fnemo(nx,ny))
      ALLOCATE(mask_nemo(nx,ny))
!
!     Reading arrays from NETCDF file
      ierr = nf90_inq_varid(ncid,'lon',varid)
      ierr = nf90_inquire_variable(ncid, varid, ndims = numdims)
      IF (numdims.EQ.2) THEN
        ierr = nf90_get_var(ncid,varid,xnemo)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,ynemo)
      ELSE IF (numdims.EQ.1) THEN
        ALLOCATE(xdum(nx))
        ALLOCATE(ydum(ny))
        ierr = nf90_get_var(ncid,varid,xdum)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,ydum)
        DO I=1, ny
          xnemo(:,I) = xdum(:)
        ENDDO
        DO I=1, nx
          ynemo(I,:) = ydum(:)
        ENDDO
        DEALLOCATE(xdum)
        DEALLOCATE(ydum)
      ELSE
        STOP 'No. of coordinate dimension neither 1 nor 2 --> STOP!'
      ENDIF

      ierr = nf90_close(ncid)
      XMIN = MINVAL(xnemo)
      XMAX = MAXVAL(xnemo)
      YMIN = MINVAL(ynemo)
      YMAX = MAXVAL(ynemo)
      WRITE(*,*) "NEMO-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "             Lat ", YMIN,' - ', YMAX 
!
!     ******** NN interpolation of HD mask to NEMO grid
      CALL remapnn(nlon, nlat, xlon, xlat, fdata, xmiss, nx, ny, xnemo, ynemo, fnemo)
      mask_nemo(:,:) = FLOOR(fnemo(:,:) + 0.001)
!
!     ******** Writing of NN interpolated mask on NEMO grid
!
!     *** WRITE NETCDF output
      ierr = nf90_create(DNOUT, NF90_CLOBBER, ncid)
      ierr = nf90_def_dim(ncid, 'x', nx, dimids(1))
      ierr = nf90_def_dim(ncid, 'y', ny, dimids(2))

      ierr = nf90_def_var(ncid, 'mask_hd_mouth', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 731)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','HD mouth points on NEMO grid')
      ierr = nf90_put_att(ncid,varid,'standard_name','MASK_HD_MOUTH')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_put_att(ncid,varid,'missing_value',-9999.)
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, mask_nemo)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_FLOAT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_east')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, xnemo)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_FLOAT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_north')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, ynemo)
      ierr = nf90_close(ncid)
      WRITE(*,*) "HD mouths on NEMO grid are written to: ", TRIM(DNOUT)
!
      DEALLOCATE(mask_nemo)
      DEALLOCATE(fnemo)
      DEALLOCATE(xnemo)
      DEALLOCATE(ynemo)
      DEALLOCATE(fdata)
      DEALLOCATE(xlon)
      DEALLOCATE(xlat)
!
!     ******** Programmende
  999 END
!
