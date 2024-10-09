! generate_coast.f90 - Generates mask of coastal areas over land within a radius from the coast
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

      PROGRAM GENERATE_COAST
!
!     ******** Programm that generates coastal area over land within distance of 
!              dmax = ?? km (e.g. 25 km) from the next ocean grid box
!
!     ******** Version 1.0 - June 2019
!              Programmierung und Entwicklung: Stefan Hagemann 
!              It writes file: coastal.nc
!
!
!     ******** Variablenliste:
!     ***  
!     ***  dninp = Name of input land sea mask: Set to lsm.nc
!     ***  dnout = Name of output coastal mask: Set to coast.nc
!     ***   IQUE = Kommentarvariable ( 0 = Kein Kommentar )
!     ***     LU = Logical Unit der zu lesenden Datei
!     *** fcoast = Muendungsmaske on NEMO Grid
!     ***   fmas = Sea mask NEMO
!     ***  
      use netcdf
      USE mo_interpol,              ONLY: distance
!
      PARAMETER(xmiss=-9999.) 
!      DOUBLE PRECISION, PARAMETER  :: dmax = 25000.    ! Maximum distance to coast in [m]
!      DOUBLE PRECISION, PARAMETER  :: dmax = 100000.    ! Maximum distance to coast in [m]
      DOUBLE PRECISION  :: dmax
      INTEGER NX, NY, IHDDIM, I
      INTEGER jb, jl, jl1, jl2, jb1, jb2
      INTEGER, DIMENSION(2) :: ICOORD
!     Grid arrays
      REAL, DIMENSION(:,:), ALLOCATABLE :: flon    ! Longitudes
      REAL, DIMENSION(:,:), ALLOCATABLE :: flat    ! Latitudes
      REAL, DIMENSION(:,:), ALLOCATABLE :: fmas
      REAL, DIMENSION(:,:), ALLOCATABLE :: fcoast
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  fdist 
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: lland
      REAL, DIMENSION(:), ALLOCATABLE :: vlon   ! Longitudes in ocean source file 
      REAL, DIMENSION(:), ALLOCATABLE :: vlat   ! Latitudes in ocean source file
      DOUBLE PRECISION :: xlon, xlat
      CHARACTER (len=4) :: cmax
      REAL XMAX, XMIN, YMAX, YMIN
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, dimid, varid, dimids(2), numdims
!  
      CHARACTER (len=128) :: dninp, ORDER, dnout
!
      dninp="lsm.nc"
      dnout="coast.nc"
      LU=50
      ! getarg  for dmax in [km]
      CALL GETARG(1,cmax)    ! Exp. No. for simulated/forcing precip. 
      READ(cmax,'(I4)') I
      dmax = DBLE(I) * 1000.
      WRITE(*,*) 'dmax = ', dmax
!	  
!     ******** READ lsm Mask points and coordinates
      ierr = nf90_open(dninp, NF90_NOWRITE, ncid)
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
      IF (NX.EQ.0 .OR. NY.EQ.0) THEN
        ierr = nf90_inq_dimid(ncid,'rlon',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NX)
        ierr = nf90_inq_dimid(ncid,'rlat',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NY)
      ENDIF
      ! Testing --> NX = 619, NY=523
      WRITE(*,*) "Sourcegrid of land sea mask:  NX =" , NX, "  NY = ", NY
!
!     *** Feld Dimensionierung
      ALLOCATE(flon(NX,NY))
      ALLOCATE(flat(NX,NY))
      ALLOCATE(fmas(NX,NY))
      ALLOCATE(fcoast(NX,NY))
      ALLOCATE(fdist(NX,NY))
      ALLOCATE(lland(NX,NY))
!
!     Reading arrays from NETCDF file
      ierr = nf90_inq_varid(ncid,'mask',varid)
      ierr = nf90_get_var(ncid,varid,fmas)

      ierr = nf90_inq_varid(ncid,'lon',varid)
      ierr = nf90_inquire_variable(ncid, varid, ndims = numdims)
      IF (numdims.EQ.2) THEN
        ierr = nf90_get_var(ncid,varid,flon)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,flat)
      ELSE IF (numdims.EQ.1) THEN
        ALLOCATE(vlon(NX))
        ALLOCATE(vlat(NY))
        ierr = nf90_get_var(ncid,varid,vlon)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,vlat)
        DO I=1, NY
          flon(:,I) = vlon(:)
        ENDDO
        DO I=1, NX
          flat(I,:) = vlat(:)
        ENDDO
        DEALLOCATE(vlon)
        DEALLOCATE(vlat)
      ELSE
        STOP 'No. of coordinate dimension neither 1 nor 2 --> STOP!'
      ENDIF
      ierr = nf90_close(ncid)
      XMIN = MINVAL(flon)
      XMAX = MAXVAL(flon)
      YMIN = MINVAL(flat)
      YMAX = MAXVAL(flat)
      WRITE(*,*) "  LSM-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "              Lat ", YMIN,' - ', YMAX 
      WRITE(*,*) "             Mask ", MINVAL(fmas),' - ', MAXVAL(fmas)
!
!     *** Land MASK
      lland(:,:) = .FALSE.
      WHERE (ABS(fmas(:,:)-1).LE.0.1)
        lland = .TRUE.
      END WHERE 
      WRITE(*,*) ' Land boxes: ', COUNT(lland) 
!
!     ******** 
      fcoast(:,:) = xmiss
      icoast=0
      DO jb = 1, NY
      DO jl = 1, NX
      IF (.NOT. lland(jl,jb)) THEN           ! Ocean point ?
        jb1 = jb-1 ; jb2 = jb + 1
        jl1 = jl-1 ; jl2 = jl + 1
        IF (jb.EQ.1) jb1=1
        IF (jb.EQ.NY) jb2=NY
        IF (jl.EQ.1) jl1=1
        IF (jl.EQ.NX) jl2=NX
 
!!        WRITE(*,*) jl,jb, lland(jl1:jl2,jb1:jb2)

        IF (ANY(lland(jl1:jl2,jb1:jb2))) THEN
          icoast = icoast + 1
          fcoast(jl,jb) = 0
          xlon = DBLE(flon(jl,jb))
          xlat = DBLE(flat(jl,jb))
          CALL distance(xlon, xlat, NX, NY, lland, flon, flat, fdist)
          WHERE (lland .AND. fdist.LE.dmax)
            fcoast(:,:) = 1
          END WHERE
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      ncoast = icoast
      WRITE(*,*) ncoast, ' coastal ocean points in lsm grid found'
      IF (ncoast.EQ.0) STOP 'No coastal ocean point found/generated --> ERROR'
!
!     ******** Schreiben des coastal array
!
!     *** WRITE NETCDF output
      ierr = nf90_create(DNOUT, NF90_CLOBBER, ncid)
      ierr = nf90_def_dim(ncid, 'x', NX, dimids(1))
      ierr = nf90_def_dim(ncid, 'y', NY, dimids(2))

      ierr = nf90_def_var(ncid, 'coast', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 731)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','Coastal area over land')
      ierr = nf90_put_att(ncid,varid,'standard_name','fcoast')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_put_att(ncid,varid,'missing_value',-9999.)
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, fcoast)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_FLOAT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_east')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, flon)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_FLOAT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_north')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, flat)
      ierr = nf90_close(ncid)
      WRITE(*,*) "Coastal land points were written: ", TRIM(dnout)
!
      DEALLOCATE(fcoast)
      DEALLOCATE(fmas)
      DEALLOCATE(flon)
      DEALLOCATE(flat)
!
!     ******** Programmende
  999 END
!
