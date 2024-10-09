! hdtonemo.f90 - Converts HD discharge to NEMO inflows
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

      PROGRAM HDTONEMO
!
!     ******** Programm that takes the HD river mouth points that have been created from the
!              River-Direction-File (HD para Array FDIR) beforehand, either on the HD grid
!              or already transformed to the NEMO grid. Then it finds the 
!              the nearest river mouth points on the NEMO grid. 
!              The latter is read in from file into MASK_NEMO. The program 
!              creates an HD array with the NEMO lon and lat indices that
!              allocate these nearest mouth points to each HD mouth point.
!
!     ******** Output
!       hd_to_nemo_mouth.nc  --> hd_to_ocean_mouth.nc
!           Variable: FMOU_HD_TO_NEMO   HD Mouth point with an associated ocean model, e.g. NEMO, mouth point
!                                       on the source grid (HD or ocean model grid)
!                     INDEXX            Ocean model longitude index of associated mouth for the HD mouth point
!                     INDEXY            Ocean model latitude index of associated mouth for the HD mouth point
!       nemo_hdmouth.nc  --> hdmouth_on_oceangrid.nc
!           Variable: FMOU_HD_ON_NEMO   Mask of associated NEMO mouth points on NEMO grid
!
!     ******** Version 1.0 - Dec. 2016
!              Programmierung und Entwicklung: Stefan Hagemann 
!              Input file:    rivmouth_source.nc, river_mouths_NEMO.nc
!              Output files:  hd_nemo_mouth.nc, nemo_hdmouth,nc
!
!     ******** Version 1.1 - Jan. 2017
!              Output files:  hd_to_nemo_mouth.nc, hdmouth_on_nemo.nc
!
!     ******** Version 1.2 - Nov. 2017
!              Coupling adaptations
!
!     ******** Version 1.3 - April 2018
!              Based on command line parameter INEMOU, two NEMO masks may be used (INEMOU=3):
!              MASK_NEMO_PRESET (primary mask) and MASK_NEMO (secondary mask), e.g.:
!                MASK_NEMO_PRESET (predifined mask) and MASK_NEMO (coastal points based on sea mask)
!
!     ******** Version 1.4 - Oct. 2018 - Generalization
!              Input file:   rivmouth_source.nc
!                            river_mouths_oceangrid.nc     0: Mouth, 1: others
!              Output file:  hd_to_ocean_mouth.nc, hdmouth_on_oceangrid.nc
!
!     ******** Version 1.5 - Apr. 2020 - Include utilization of ICON ocean grid
!
!     ******** HD Direction format
!
!                               7  8  9
!                                \ | /
!                                 \|/
!                               4--5--6
!                                 /|\
!                                / | \
!                               1  2  3   
!                                  
!     ***            
!     ***        Anmerkung: Richtung 5 = Discharge-Trap
!     ***                   Richtung 0 = Muendungspunkt
!     ***                           -1 = Seepunkt, aber keine Muendung
!
!
!     ******** Variablenliste:
!     ***  
!     ***  DNDIR = Name des RDF
!     ***   IQUE = Kommentarvariable ( 0 = Kein Kommentar )
!     *** FMOUTH = Muendungsmaske HD
!     ***   MASK_NEMO = Muendungsmaske NEMO
!     ***   FDIR = HD Direction array (not used)
!     *** INDEXX = x-Indices of nearest NEMO mouth points
!     *** INDEXY = y-Indices of nearest NEMO mouth points
!     ***  
!     *** DISMAX = Distance threshold
!     *** DEGMAX = Maximum Distance in degree from outside the NEMO region
!     ***
      use netcdf
!
      PARAMETER(XMISS=-9999.) 
      PARAMETER(PI=3.14159265) 
      INTEGER, PARAMETER  :: IQUE =0     ! Print out some debug statements for IQUE=1
      INTEGER :: NL, NB, NX, NY
      INTEGER :: ITYP=0                  ! Unique directions 0/1 = NO/YES (Ha wants 1)
      INTEGER, DIMENSION(2) :: ICOORD
!     HD arrays
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_SRC
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_MOU   ! Mask with HD mouths with a valid nearest NEMO mouth
      REAL, DIMENSION(:,:), ALLOCATABLE :: FMOUTH
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FLON
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FLAT
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VLON
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VLAT
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDEXX
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: INDEXY
!     Ocean model (e.g. NEMO) arrays
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_NEMO        ! Note: this is an inverted mask, 0=mouth, 1 otherwise
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_NEMO_PRESET ! Note: this is an inverted mask, 0=mouth, 1 otherwise
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: MASK_ON_OCEAN
      REAL, DIMENSION(:,:), ALLOCATABLE :: XNEMO
      REAL, DIMENSION(:,:), ALLOCATABLE :: YNEMO
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, dimid, varid, dimids(2), numdims
!  
      CHARACTER (len=128) :: DNDIR, DNNEMO, ORDER, DNOUT
      CHARACTER (len=2) :: CNEMOU  ! Command line input parameter for program   
      INTEGER :: INEMOU            !  --> Method of using potential masks
      INTEGER :: NR_MASK = 1       ! No. of masks (1 for INEMOU =1,2 or 2 for INEMOU = 3)
      CHARACTER (len=2) :: COCEAN  ! Command line input parameter for ID IOCEAN   
      INTEGER :: IOCEAN = 1        ! Ocean model setup ID: 
                                   !   1    NEMO Vs. 3.3 - North & Baltic Seas
                                   !   2    ECOSMO Vs. 3 - North & Baltic Seas
                                   !   3    SCHISM
                                   !   4    ECOSMO Vs. II - North & Baltic Seas
                                   !   5    ICON
      CHARACTER (len=10), DIMENSION(5) :: CNAME_OC   ! Ocean model names Array
      CHARACTER (len=10)               :: CONAME   ! Ocean model name
!
      DATA CNAME_OC / 'NEMO', 'ECOSMO3', 'SCHISM', ' ECOSMO2', 'ICON' / 
!
      CALL GETARG(1,CNEMOU)
      READ(CNEMOU, '(I2)') INEMOU
!
      IF (INEMOU.LT.1 .OR. INEMOU.GT.3) THEN
         WRITE(*,*) ' INEMOU out of range [1,3]: ', INEMOU
         STOP 'ERROR --> TERMINATED!'
      ENDIF
      CALL GETARG(2, COCEAN)
      READ(COCEAN, '(I2)') IOCEAN
      IF (IOCEAN.LT.1 .OR. IOCEAN.GT.5) THEN
         WRITE(*,*) ' IOCEAN out of range [1,5]: ', IOCEAN
         STOP 'ERROR --> TERMINATED!'
      ELSE IF (IOCEAN.EQ.1) THEN
!!        ITYP=1
!!        WRITE(*,*) ' --> ONLY unique directions are allowed'
      ENDIF
      CONAME = CNAME_OC(IOCEAN)
!
      DNDIR="rivmouth_source.nc"
      DNNEMO="river_mouths_oceangrid.nc"
!
!     ******** READ HD array FMOUTH
      ierr = nf90_open(DNDIR, NF90_NOWRITE, ncid)

      ierr = nf90_inq_dimid(ncid,'lon',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NL)
      ierr = nf90_inq_dimid(ncid,'lat',dimid)
      ierr = nf90_inquire_dimension(ncid,dimid,len=NB)
      WRITE(*,*) "HD Model: NL =" , NL, "  NB = ", NB
!
!     *** Feld Dimensionierung
      ALLOCATE(FMOUTH(NL,NB))
      ALLOCATE(MASK_SRC(NL,NB))
      ALLOCATE(MASK_MOU(NL,NB))
      ALLOCATE(INDEXX(NL,NB))
      ALLOCATE(INDEXY(NL,NB))
      ALLOCATE(FLON(NL,NB))
      ALLOCATE(FLAT(NL,NB))

      ierr = nf90_inq_varid(ncid,'FMOUTH',varid)
      ierr = nf90_get_var(ncid,varid,FMOUTH)
!
!     *** Dimension of Lon/Lat arrays of HD model mouth array: currently 1
      WRITE(*,*) "Read logitudinal coordinates, e.g. VLON"
      ierr = nf90_inq_varid(ncid,'lon',varid)
      ierr = nf90_inquire_variable(ncid, varid, ndims = numdims)
      IF (numdims.EQ.2) THEN
        ierr = nf90_get_var(ncid,varid,FLON)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,FLAT)
      ELSE IF (numdims.EQ.1) THEN
        ALLOCATE(VLON(NL))
        ALLOCATE(VLAT(NB))
        ierr = nf90_get_var(ncid,varid,VLON)
        WRITE(*,*) "Read VLAT"
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,VLAT)

        WRITE(*,*) "VLON:", VLON(1), VLON(NL)
        WRITE(*,*) "VLAT:", VLAT(1), VLAT(NB)
        DO JL=1, NL
          FLAT(JL,:) = VLAT(:)
        ENDDO
        DO JB=1, NB
          FLON(:,JB) = VLON(:)
        ENDDO
      ELSE
        STOP 'No. of coordinate dimension neither 1 nor 2 --> STOP!'
      ENDIF

      ierr = nf90_close(ncid)
      WRITE(*,*) TRIM(DNDIR)," is closed."

      WRITE(*,*) "FMOUTH ausgelesen --> Origin: ", FLON(1,1), FLAT(1,1)
      WRITE(*,*) "                 Lower right: ", FLON(NL,NB), FLAT(NL,NB)
!
!!        FLON(JL,JB) = FLOAT(JL)/ 2. - 180.25 
!!        FLAT(JL,JB) = 90.25 - FLOAT(JB) / 2.
!	  
!     ******** READ Ocean grid Mouth points and coordinates
      ierr = nf90_open(DNNEMO, NF90_NOWRITE, ncid)
!
!     *** 2D grid or ICON grid
      IF (IOCEAN.EQ.5) THEN    ! ICON grid --> NX=ncells, NY=1
        ierr = nf90_inq_dimid(ncid,'ncells',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NX)
        NY = 1
      ELSE    ! 2D grid
        ierr = nf90_inq_dimid(ncid,'x',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NX)
        ierr = nf90_inq_dimid(ncid,'y',dimid)
        ierr = nf90_inquire_dimension(ncid,dimid,len=NY)
      ENDIF
      ! Testing --> NX = 619, NY=523
      WRITE(*,*) "Ocean model:  NX =" , NX, "  NY = ", NY
!
!     *** Feld Dimensionierung
      ALLOCATE(MASK_NEMO(NX,NY))
      IF (INEMOU.EQ.3) ALLOCATE(MASK_NEMO_PRESET(NX,NY))
      ALLOCATE(MASK_ON_OCEAN(NX,NY))
      ALLOCATE(XNEMO(NX,NY))
      ALLOCATE(YNEMO(NX,NY))
!
!     Reading arrays from NETCDF file
      ierr = nf90_inq_varid(ncid,'nmhd_msk',varid)
      ierr = nf90_get_var(ncid,varid,MASK_NEMO)

      WRITE(*,*) 'Mouth array on ocean grid - MASK_NEMO: ', &
                  MINVAL(MASK_NEMO), ' - ', MAXVAL(MASK_NEMO)
      WRITE(*,*) '         SUM: ', SUM(MASK_NEMO)

      IF (INEMOU.EQ.3) THEN
        NR_MASK=2
        ierr = nf90_inq_varid(ncid,'nemo_mask_preset',varid)
        ierr = nf90_get_var(ncid,varid,MASK_NEMO_PRESET)
      ENDIF

      IF (IOCEAN.EQ.5) THEN    ! ICON grid --> NX=ncells, NY=1
        ierr = nf90_inq_varid(ncid,'clon',varid)
        ierr = nf90_get_var(ncid,varid,XNEMO)
        ierr = nf90_inq_varid(ncid,'clat',varid)
        ierr = nf90_get_var(ncid,varid,YNEMO)
        XNEMO = XNEMO / PI * 180.
        WHERE (XNEMO.GT.180) XNEMO = XNEMO - 360.
        WHERE (XNEMO.LT.-180) XNEMO = XNEMO + 360.
        YNEMO = YNEMO / PI * 180.
      ELSE
        ierr = nf90_inq_varid(ncid,'lon',varid)
        ierr = nf90_get_var(ncid,varid,XNEMO)
        ierr = nf90_inq_varid(ncid,'lat',varid)
        ierr = nf90_get_var(ncid,varid,YNEMO)
      ENDIF
      ierr = nf90_close(ncid)

      WRITE(*,*) "Ocean Longitudes:", MINVAL(XNEMO), ' - ', MAXVAL(XNEMO)
      WRITE(*,*) "Ocean Latitudes: ", MINVAL(YNEMO), ' - ', MAXVAL(YNEMO)
!
!     ******** Source and target grid are the same
      IF (NL.EQ.NX .AND. NB.EQ.NY) THEN
!
!       *** Prepare HD mouth mask as in Ocean model: 0 for mouth, 1 otherwise
        WHERE( ABS(FMOUTH(:,:)-1).LE.0.1 ) 
           MASK_SRC(:,:) = 0
        ELSEWHERE
           MASK_SRC(:,:) = 1
        END WHERE
        WRITE(*,*) 'HD Mouth points on Ocean model grid: ', NX*NY - SUM(MASK_SRC)
!
!!        print*," XNEMO=" , XNEMO
!!        pause
!!        print*," YNEMO=" , YNEMO
!!       pause
!       ******** Search for nearest NEMO mouth point and generate Index arrays
        CALL FINDNEMOMOUTH(NX, NY, ITYP, XNEMO, YNEMO, MASK_SRC, &
             NR_MASK, MASK_NEMO_PRESET, MASK_NEMO, XMISS, IQUE,  &
             MASK_MOU, INDEXX, INDEXY, MASK_ON_OCEAN)
      ELSE     ! Source grid differs from target grid
!
!       ******** Search for nearest NEMO mouth point and generate Index arrays
        CALL FINDNEMO_FROMHD(NL, NB, ITYP, IOCEAN, FLON, FLAT, FMOUTH, &
             NX, NY, XNEMO, YNEMO, &
             NR_MASK, MASK_NEMO_PRESET, MASK_NEMO, XMISS, IQUE,  &
             MASK_MOU, INDEXX, INDEXY, MASK_ON_OCEAN)
        WRITE(*,*) 'HD Mouth points with unique ocean model gridbox target: ', SUM(MASK_MOU, MASK_MOU.EQ.1)
        WRITE(*,*) 'HD Mouth points with shared ocean model gridbox targets: ', COUNT(MASK_MOU.gt.1.5)
      ENDIF
!
!     ********* Correction for Specific grid boxes
      CALL CORROCEAN(IOCEAN, NL, NB, MASK_MOU, INDEXX, INDEXY, NX, NY, MASK_ON_OCEAN)

      WRITE(*,*) 'Maximum sources an ocean model gridbox is sharing: ', MAXVAL(MASK_MOU)
!
!     ******** Schreiben der Index Arrays with ocean model River-Mouth targets
!
!     *** WRITE NETCDF output
      DNOUT="hd_to_ocean_mouth.nc"
      ierr = nf90_create(DNOUT, NF90_CLOBBER, ncid)
      ierr = nf90_def_dim(ncid, 'lon', NL, dimids(1))
      ierr = nf90_def_dim(ncid, 'lat', NB, dimids(2))

      ierr = nf90_def_var(ncid, 'FMOU_HD_TO_NEMO', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 731)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','HD river mouths pointing towards ocean grid')
      ierr = nf90_put_att(ncid,varid,'standard_name','FMOU_HD_TO_NEMO')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_put_att(ncid,varid,'missing_value',-9999.)
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, MASK_MOU)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'INDEXX', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 732)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','Ocean model longitude index of HD river mouth')
      ierr = nf90_put_att(ncid,varid,'standard_name','INDEXX')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, INDEXX)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'INDEXY', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 733)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','Ocean model latitude index of HD river mouth')
      ierr = nf90_put_att(ncid,varid,'standard_name','INDEXY')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, INDEXY)
      IF (numdims.EQ.1) THEN
        ierr = nf90_redef(ncid)
        ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimids(1), varid)
        ierr = nf90_put_att(ncid,varid,'units','degrees_east')
        ierr = nf90_put_att(ncid,varid,'long_name','longitude')
        ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
        ierr = nf90_put_att(ncid,varid,'axis','X')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, VLON)

        ierr = nf90_redef(ncid)
        ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimids(2), varid)
        ierr = nf90_put_att(ncid,varid,'units','degrees_north')
        ierr = nf90_put_att(ncid,varid,'long_name','latitude')
        ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
        ierr = nf90_put_att(ncid,varid,'axis','Y')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, VLAT)
        DEALLOCATE(VLON)
        DEALLOCATE(VLAT)
      ELSE
        ierr = nf90_redef(ncid)
        ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimids, varid)
        ierr = nf90_put_att(ncid,varid,'units','degrees_east')
        ierr = nf90_put_att(ncid,varid,'long_name','longitude')
        ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, FLON)

        ierr = nf90_redef(ncid)
        ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimids, varid)
        ierr = nf90_put_att(ncid,varid,'units','degrees_north')
        ierr = nf90_put_att(ncid,varid,'long_name','latitude')
        ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
        ierr = nf90_enddef(ncid)
        ierr = nf90_put_var(ncid, varid, FLAT)
       ENDIF
       ierr = nf90_close(ncid)
      WRITE(*,*) "River-Mouths pointing towards ocean grid were written: ", TRIM(DNOUT)
!
      DNOUT="hdmouth_on_oceangrid.nc"
      ierr = nf90_create(DNOUT, NF90_CLOBBER, ncid)
      ierr = nf90_def_dim(ncid, 'x', NX, dimids(1))
      ierr = nf90_def_dim(ncid, 'y', NY, dimids(2))

      ierr = nf90_def_var(ncid, 'FMOU_HD_ON_NEMO', NF90_INT, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'code', 741)
      ierr = nf90_put_att(ncid,varid,'units','[]')
      ierr = nf90_put_att(ncid,varid,'long_name','HD river mouth points on ocean model grid')
      ierr = nf90_put_att(ncid,varid,'standard_name','FMOU_HD_ON_NEMO')
      ierr = nf90_put_att(ncid,varid,'coordinates','lon lat')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, MASK_ON_OCEAN)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lon', NF90_DOUBLE, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_E')
      ierr = nf90_put_att(ncid,varid,'long_name','longitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','longitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, XNEMO)

      ierr = nf90_redef(ncid)
      ierr = nf90_def_var(ncid, 'lat', NF90_DOUBLE, dimids, varid)
      ierr = nf90_put_att(ncid,varid,'units','degrees_N')
      ierr = nf90_put_att(ncid,varid,'long_name','latitude')
      ierr = nf90_put_att(ncid,varid,'standard_name','latitude')
      ierr = nf90_enddef(ncid)
      ierr = nf90_put_var(ncid, varid, YNEMO)

      ierr = nf90_close(ncid)
      WRITE(*,*) "River-Mouths on ocean model grid were written: ", TRIM(DNOUT)

      DEALLOCATE(MASK_MOU)
      DEALLOCATE(FMOUTH)
      DEALLOCATE(FLON)
      DEALLOCATE(FLAT)
      DEALLOCATE(INDEXX)
      DEALLOCATE(INDEXY)
      DEALLOCATE(MASK_NEMO)
      IF (INEMOU.EQ.3) DEALLOCATE(MASK_NEMO_PRESET)
      DEALLOCATE(XNEMO)
      DEALLOCATE(YNEMO)
!
!     ******** Programmende
      CONTAINS
!
!     *********************************************************************
      SUBROUTINE DISTANCE(XLON, XLAT, NX, NY, LMASK,    &
                 XNEMO, YNEMO, FDIS)
!     *********************************************************************
!
!     Distance of XNEMO,YNEMO to XLON, XLAT for Points in MASK_NEMO

      DOUBLE PRECISION, INTENT(in) :: XLON, XLAT
      INTEGER, INTENT(in) :: NX, NY
      LOGICAL, DIMENSION(NX,NY), INTENT(in) :: LMASK
      REAL, DIMENSION(NX,NY), INTENT(in) :: XNEMO
      REAL, DIMENSION(NX,NY), INTENT(in) :: YNEMO
      DOUBLE PRECISION, DIMENSION(NX,NY), INTENT(out) :: FDIS
      DOUBLE PRECISION, PARAMETER  ::   RERDE = 6371000.
      DOUBLE PRECISION, PARAMETER  ::   PI = 3.141592653589793
      DOUBLE PRECISION, PARAMETER  ::   PIFAK = PI/180.

      DOUBLE PRECISION :: DLON, DLAT
      INTEGER :: IX, IY
!
      FDIS(:,:) = -9999.
      DO IX=1, NX 
      DO IY=1, NY 
         IF (LMASK(IX, IY)) THEN
	    ! Pay regard to sphere
           DLON = ABS(XNEMO(IX,IY)-XLON)
           IF (DLON.GT.180.) DLON = 360. - ABS(XNEMO(IX,IY)) - ABS(XLON)
           DLON = DCOS( (YNEMO(IX,IY)+XLAT) / 2. *PIFAK ) * RERDE * DLON * PIFAK
           DLAT = ABS(YNEMO(IX,IY)-XLAT) * PIFAK * RERDE
            
	   FDIS(IX, IY) = SQRT(DLON*DLON + DLAT*DLAT)
	 ENDIF
      ENDDO	  
      ENDDO	  
!
      END SUBROUTINE
!
!     *********************************************************************
      SUBROUTINE FINDNEMOMOUTH(NX, NY, ITYP, XNEMO, YNEMO, MASK_SRC, &
             NR_MASK, MASK_NEMO1, MASK_NEMO2, XMISS, IQUE,  &
             MASK_MOU, INDEXX, INDEXY, MASK_ON_OCEAN)
!     *********************************************************************
!
      IMPLICIT NONE
!
      INTEGER, INTENT(in) :: NX, NY
      INTEGER, INTENT(in) :: ITYP                          ! Unique directions 0/1 = NO/YES
      REAL, DIMENSION(NX,NY), INTENT(in) :: XNEMO          ! NEMO Longitudes
      REAL, DIMENSION(NX,NY), INTENT(in) :: YNEMO          ! NEMO Latitudes
      INTEGER, DIMENSION(NX,NY), INTENT(in) :: MASK_SRC    ! HD mouths 0, otherwise 1
      INTEGER, INTENT(in) :: NR_MASK                       ! No. of masks: 1 or 2
                                                           ! 1 mask -> only secondary mask is used
      INTEGER, DIMENSION(NX,NY), INTENT(in) :: MASK_NEMO1  ! Primary mask: Nemo mouths 0, otherwise 1
      INTEGER, DIMENSION(NX,NY), INTENT(in) :: MASK_NEMO2  ! Secondary mask: Nemo mouths 0, otherwise 1
      REAL, INTENT(in)    :: XMISS    ! Missing Value, -9999. is suitable, depends of definition.
      INTEGER, INTENT(in) :: IQUE                       ! Write some debug statements: 0/1 = No/yes

      INTEGER, DIMENSION(NX,NY), INTENT(out) :: MASK_MOU   ! Mask with HD mouths with a valid nearest NEMO mouth
      INTEGER, DIMENSION(NX,NY), INTENT(out) :: INDEXX     ! x-Indices of nearest NEMO mouth points
      INTEGER, DIMENSION(NX,NY), INTENT(out) :: INDEXY     ! y-Indices of nearest NEMO mouth points
      INTEGER, DIMENSION(NX,NY), INTENT(out) :: MASK_ON_OCEAN   ! Mask with associated NEMO mouth points
!
!
      DOUBLE PRECISION :: XLON, XLAT
!     *** Maximum allowed distances to next mask point
      DOUBLE PRECISION, PARAMETER :: DISMAX1=25000.    ! 25 km for primary mask
      DOUBLE PRECISION, PARAMETER :: DISMAX2=100000.   ! 100 km for secondary mask (default for only 1 mask)
      DOUBLE PRECISION, PARAMETER :: DEGMAX=0.5       ! 0.5 degree 
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LMASK1  ! Logical array of primary mask
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LMASK2  ! Logical array of secondary mask
      LOGICAl  :: LFOUND
      INTEGER :: JL, JB, NMES, NMOU
      REAL XMAX, XMIN, YMAX, YMIN
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FDIS
      DOUBLE PRECISION :: DMIN
      INTEGER, DIMENSION(2) :: ICOORD
!     *** Variables for achieving unique pointers
      INTEGER  NMIX                                   ! NMIX = Maximum number of HD mouths pointing to a single NEMO mouth
      INTEGER JX, JY, IBOX, I, NACT, IL,IB
      REAL, DIMENSION(:,:), ALLOCATABLE :: DMIX
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ICMIX
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LFREE
!
      ALLOCATE(FDIS(NX,NY))
!
!     *** Domain Boundaries
      XMIN = MINVAL(XNEMO)
      XMAX = MAXVAL(XNEMO)
      YMIN = MINVAL(YNEMO)
      YMAX = MAXVAL(YNEMO)
      WRITE(*,*) "Ocean model-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "                    Lat ", YMIN,' - ', YMAX 
!
!     *** Note that the NEMO mouth masks defines mouth points as Zero.
      ALLOCATE(LMASK2(NX,NY))
      LMASK2(:,:) = .FALSE.
      WHERE(MASK_NEMO2.EQ.0)        
        LMASK2=.true.
      END WHERE
      IF (NR_MASK.EQ.2) THEN
        ALLOCATE(LMASK1(NX,NY))
        LMASK1(:,:) = .FALSE.
        WHERE(MASK_NEMO1.EQ.0)        
          LMASK1=.true.
        END WHERE
      ENDIF
!
!     ******** Searching for nearest NEMO mouth point
      INDEXX(:,:) = 0
      INDEXY(:,:) = 0
      MASK_MOU(:, :) = NINT(XMISS)
      MASK_ON_OCEAN(:,:)=0
      NMES=0
      NMOU=0
      DO JB = 1, NY
      DO JL = 1, NX
      IF ( MASK_SRC(JL,JB).EQ.0 ) THEN
        IF (XNEMO(JL,JB).GE.XMIN-DEGMAX .AND. XNEMO(JL,JB).LE.XMAX+DEGMAX .AND. &
            YNEMO(JL,JB).GE.YMIN-DEGMAX .AND. YNEMO(JL,JB).LE.YMAX+DEGMAX) THEN
          NMOU = NMOU+1
          MASK_MOU(JL, JB) = 0
          XLON = XNEMO(JL,JB)
          XLAT = YNEMO(JL,JB)

          LFOUND = .FALSE.
          IF (NR_MASK.EQ.2) THEN
  	    CALL DISTANCE(XLON, XLAT, NX, NY, LMASK1,  &
                 XNEMO, YNEMO, FDIS)
            ICOORD = MINLOC(FDIS, MASK=LMASK1)
            DMIN = MINVAL(FDIS, MASK=LMASK1)
            IF (DMIN.LE.DISMAX1) LFOUND = .TRUE.
            IF (JB.GT.730 .AND. JL.GT.700) THEN
              WRITE(*, '(A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                  'M1: - Lon/Lat: ', JL, JB, ' = ', XLON,XLAT, ' Distance HD-NEMO mouth = ', DMIN/1000., ' km at ',  ICOORD
              WRITE(*,*) '     -> Lon = ', XNEMO(ICOORD(1), ICOORD(2)),  ' Lat = ', YNEMO(ICOORD(1),ICOORD(2))
            ENDIF
          ENDIF
          IF (.NOT. LFOUND) THEN
            CALL DISTANCE(XLON, XLAT, NX, NY, LMASK2,  &
                 XNEMO, YNEMO, FDIS)
            ICOORD = MINLOC(FDIS, MASK=LMASK2)
            DMIN = MINVAL(FDIS, MASK=LMASK2)
            IF (DMIN.LE.DISMAX2) LFOUND = .TRUE.
            IF (JB.GT.730 .AND. JL.GT.700) THEN
              WRITE(*, '(A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                  'M2 - Lon/Lat: ', JL, JB, ' = ', XLON,XLAT,  ' Distance HD-NEMO mouth = ', DMIN/1000., ' km at ',  ICOORD
              WRITE(*,*) '     -> Lon = ', XNEMO(ICOORD(1), ICOORD(2)),  ' Lat = ', YNEMO(ICOORD(1),ICOORD(2))
            ENDIF
          ENDIF

          IF (LFOUND) THEN
            MASK_MOU(JL, JB) = 1
            NMES = NMES+1
            INDEXX(JL,JB) = ICOORD(1)
            INDEXY(JL,JB) = ICOORD(2)
            MASK_ON_OCEAN(ICOORD(1), ICOORD(2)) = MASK_ON_OCEAN(ICOORD(1), ICOORD(2)) + 1
            IF (IQUE.EQ.1) THEN
              WRITE(*,*) NMES, '. Distance HD-NEMO mouth = ', DMIN/1000., ' km at ',  ICOORD
              WRITE(*,*) '     -> Lon = ', XNEMO(ICOORD(1), ICOORD(2)),  ' Lat = ', YNEMO(ICOORD(1),ICOORD(2))
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      WRITE(*,*) NMOU, '. HD Points on ', TRIM(CONAME), ' grid converted to ', NMES, ' ', &
                 TRIM(CONAME), ' mouth points'
      NMIX = MAXVAL(MASK_ON_OCEAN)
      WRITE(*,*) 'Maximum number of HD mouths pointing towards a single ', TRIM(CONAME),' mouth: ', NMIX
!
!     *** Indicate HD mouths that point towards a Ocean model mouth with more than one input  
      DO JB = 1, NY
      DO JL = 1, NX
      IF ( MASK_MOU(JL,JB).EQ.1 ) THEN
        IF (MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB)).GT.1) MASK_MOU(JL, JB) = 2
      ENDIF
      ENDDO
      ENDDO
!
      IF (ITYP.EQ.1) THEN
        IF (NMIX.GT.2) THEN
          WRITE(*,*) "****** ERROR - Cases with more more than 2 HD mouths: up to ", NMIX
          WRITE(*,*) "       pointing to the ocean model box are NOT implemented"
          STOP       " --> Program stops right here!!!!"
        ENDIF
        ALLOCATE( DMIX(NMIX, NMIX) )
        ALLOCATE( ICMIX(NMIX, NMIX, 2) )
        ALLOCATE( LFREE(NX,NY) )
        WHERE (MASK_ON_OCEAN(:,:).EQ.1) LMASK2(:,:) = .FALSE.
        DO JB = 1, NY
        DO JL = 1, NX
        IF ( MASK_MOU(JL,JB).EQ.2 ) THEN
          JX = INDEXX(JL,JB)
          JY = INDEXY(JL,JB)
          NACT = MASK_ON_OCEAN(JX,JY)
          XLON = XNEMO(JL,JB)
          XLAT = YNEMO(JL,JB)
          LFREE(:,:) = LMASK2(:,:)
          I=0

          WRITE(*, '(A,2(X,I4),A,2(X,F7.2) ,A,2(X,F7.2),A, 2(X,I4) )') &
                   'Mixpoint - Lon/Lat: ', JL,JB, ' = ', XLON, XLAT, ' pointing towards ', &
                    XNEMO(JX,JY), YNEMO(JX,JY), ' --> JX, JY: ', JX, JY
!
!         *** Distances 1, NACT
          DO WHILE(I.LT.NACT)
            I=I+1
	    CALL DISTANCE(XLON, XLAT, NX, NY, LFREE,  &
               XNEMO, YNEMO, FDIS)
            ICMIX(1,I,:) = MINLOC(FDIS, MASK=LFREE)
            DMIX(1,I) = MINVAL(FDIS, MASK=LFREE)
            LFREE(ICMIX(1,I,1),ICMIX(1,I,2)) = .FALSE.
          ENDDO
         
          IF (IQUE.EQ.1) WRITE(*,*) 'Dist1 ', DMIX(1,1:NACT) 
!
!         *** HD mouths poiting to the same Ocean model mouth
          IBOX=2
          DO IB = 1, NY
          DO IL = 1, NX
          IF (MASK_MOU(IL,IB).EQ.2 .AND. (IL.NE.JL .OR. IB.NE.JB) .AND.  &
             INDEXX(IL,IB).EQ.JX .AND. INDEXY(IL,IB).EQ.JY ) THEN

            IF (IQUE.EQ.1) WRITE(*,*) 'IL, IB ', IL, IB,  ' -> ', INDEXX(IL,IB), INDEXY(IL,IB)
            XLON = XNEMO(IL,IB)
            XLAT = YNEMO(IL,IB)
            LFREE(:,:) = LMASK2(:,:)
            I=0
!
!           *** Distances IBOX, NACT
            DO WHILE(I.LT.NACT)
              I=I+1
	      CALL DISTANCE(XLON, XLAT, NX, NY, LFREE,  &
                 XNEMO, YNEMO, FDIS)
              ICMIX(IBOX,I,:) = MINLOC(FDIS, MASK=LFREE)
              DMIX(IBOX,I) = MINVAL(FDIS, MASK=LFREE)
              LFREE(ICMIX(IBOX,I,1),ICMIX(IBOX,I,2)) = .FALSE.
            ENDDO

            IF (IQUE.EQ.1) WRITE(*,*) 'Dist2 ', DMIX(2,1:NACT) 
!
!           *** Only NMIX = NACT = 2 can be dealt with
            IF (DMIX(1,1).LE.DMIX(2,1)) THEN           ! Box 1 closer than 2
              MASK_MOU(JL,JB) = 1           
              IF (DMIX(2,2).LE.DISMAX2) THEN
                INDEXX(IL,IB) = ICMIX(2,2,1)
                INDEXY(IL,IB) = ICMIX(2,2,2)
                MASK_ON_OCEAN(ICMIX(2,2,1), ICMIX(2,2,2)) = 1
                LMASK2(ICMIX(2,2,1), ICMIX(2,2,2)) = .FALSE.
!!                IF (IQUE.EQ.1) THEN
                WRITE(*,*) 'New Distance HD-NEMO mouth = ', DMIX(2,2)/1000., ' km at ',  ICMIX(2,2,:)
                IF (IQUE.EQ.1) WRITE(*,*) '     -> Lon = ', XNEMO(ICMIX(2,2,1), ICMIX(2,2,2)),  ' Lat = ', YNEMO(ICMIX(2,2,1), ICMIX(2,2,2))
!!                ENDIF
                MASK_MOU(IL,IB) = 1
              ELSE
                MASK_MOU(IL,IB) = 0
                INDEXX(IL,IB) = 0
                INDEXY(IL,IB) = 0
                WRITE(*,*) 'Second closest NEMO mouth is out of range: ', DMIX(2,2)/1000., ' km at ',  ICMIX(2,2,:)
              ENDIF
            ELSE                                       ! Box 2 closer than 1
              MASK_MOU(IL,IB) = 1           
              IF (DMIX(1,2).LE.DISMAX2) THEN
                INDEXX(JL,JB) = ICMIX(1,2,1)
                INDEXY(JL,JB) = ICMIX(1,2,2)
                MASK_ON_OCEAN(ICMIX(1,2,1), ICMIX(1,2,2)) = 1
                LMASK2(ICMIX(1,2,1), ICMIX(1,2,2)) = .FALSE.
!!                IF (IQUE.EQ.1) THEN
                WRITE(*,*) 'New Distance HD-NEMO mouth = ', DMIX(1,2)/1000., ' km at ',  ICMIX(1,2,:)
                IF (IQUE.EQ.1) WRITE(*,*) '     -> Lon = ', XNEMO(ICMIX(1,2,1), ICMIX(1,2,2)),  ' Lat = ', YNEMO(ICMIX(1,2,1), ICMIX(1,2,2))
!!                ENDIF
                MASK_MOU(JL,JB) = 1
              ELSE
                MASK_MOU(JL,JB) = 0
                INDEXX(JL,JB) = 0
                INDEXY(JL,JB) = 0
                WRITE(*,*) 'Second closest NEMO mouth is out of range: ', DMIX(1,2)/1000., ' km at ',  ICMIX(1,2,:)
              ENDIF
            ENDIF                           
!
!           *** Mixed problem solved for JX,JY
            MASK_ON_OCEAN(JX, JY) = 1
            LMASK2(JX, JY) = .FALSE.
            EXIT
          ENDIF
          ENDDO
          ENDDO
        ENDIF
        ENDDO
        ENDDO
        DEALLOCATE(DMIX)
        DEALLOCATE(ICMIX)
        DEALLOCATE(LFREE)
      ENDIF
!
      DEALLOCATE(FDIS)
      IF (NR_MASK.EQ.2) DEALLOCATE(LMASK1)
      DEALLOCATE(LMASK2)
      END SUBROUTINE FINDNEMOMOUTH
!
!     *********************************************************************
      SUBROUTINE FINDNEMO_FROMHD(NL, NB, ITYP, IOCEAN, FLON, FLAT, FMOUTH, &
             NX, NY, XNEMO, YNEMO, &
             NR_MASK, MASK_NEMO1, MASK_NEMO2, XMISS, IQUE,  &
             MASK_MOU, INDEXX, INDEXY, MASK_ON_OCEAN)
!     *********************************************************************
!
      INTEGER, INTENT(in) :: NL, NB
      INTEGER, INTENT(in) :: ITYP                             ! Unique directions 0/1 = NO/YES
      INTEGER, INTENT(in) :: IOCEAN  ! Ocean model setup ID: 
                                     !   1    NEMO Vs. 3.3 - North & Baltic Seas
                                     !   2    ECOSMO Vs. 3 - North & Baltic Seas
                                     !   3    SCHISM
                                     !   4    ECOSMO Vs. II - North & Baltic Seas
                                     !   5    ICON
      DOUBLE PRECISION, DIMENSION(NL,NB), INTENT(in) :: FLON  ! HD longitudes
      DOUBLE PRECISION, DIMENSION(NL,NB), INTENT(in) :: FLAT  ! HD latitudes 
      REAL, DIMENSION(NL,NB), INTENT(in) :: FMOUTH            ! HD mouths 1, otherwise 0

      INTEGER, INTENT(in) :: NX, NY
      REAL, DIMENSION(NX,NY), INTENT(in) :: XNEMO          ! NEMO Longitudes
      REAL, DIMENSION(NX,NY), INTENT(in) :: YNEMO          ! NEMO Latitudes
      INTEGER, INTENT(in) :: NR_MASK                       ! No. of masks: 1 or 2
                                                           ! 1 mask -> only secondary mask is used
      INTEGER, DIMENSION(NX,NY), INTENT(in) :: MASK_NEMO1  ! Primary mask: Nemo mouths 0, otherwise 1
      INTEGER, DIMENSION(NX,NY), INTENT(in) :: MASK_NEMO2  ! Secondary mask: Nemo mouths 0, otherwise 1
      REAL, INTENT(in)    :: XMISS    ! Missing Value, -9999. is suitable, depends of definition.
      INTEGER, INTENT(in) :: IQUE                       ! Write some debug statements: 0/1 = No/yes

      INTEGER, DIMENSION(NL,NB), INTENT(out) :: MASK_MOU   ! Mask with HD mouths with a valid nearest NEMO mouth
      INTEGER, DIMENSION(NL,NB), INTENT(out) :: INDEXX     ! x-Indices of nearest NEMO mouth points
      INTEGER, DIMENSION(NL,NB), INTENT(out) :: INDEXY     ! y-Indices of nearest NEMO mouth points
      INTEGER, DIMENSION(NX,NY), INTENT(out) :: MASK_ON_OCEAN   ! Mask with associated NEMO mouth points
!
      DOUBLE PRECISION :: XLON, XLAT
!
!     *** Maximum allowed distances to next mask point
      DOUBLE PRECISION :: DISMAX1=25000.    ! 25 km for primary mask
      DOUBLE PRECISION :: DISMAX2=100000.   ! 100 km for secondary mask (default for only 1 mask)
      DOUBLE PRECISION :: DEGMAX=0.5        ! 0.5 degree for outside ocean domain
      REAL :: HDRES                         ! Resolution HD model
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LMASK1  ! Logical array of primary mask
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LMASK2  ! Logical array of secondary mask
      LOGICAl  :: LFOUND
      INTEGER :: JL, JB, NMES, NMOU
      INTEGER :: NOMO                       ! Number of HJD points with no ocean model mouth found
      REAL XMAX, XMIN, YMAX, YMIN
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FDIS
      DOUBLE PRECISION :: DMIN
      INTEGER, DIMENSION(2) :: ICOORD
!     *** Variables for achieving unique pointers
      INTEGER  NMIX            ! NMIX = Maximum number of HD mouths pointing to a single Oceanmodel mouth
      INTEGER JX, JY, IBOX, I, NACT, IL,IB
      REAL, DIMENSION(:,:), ALLOCATABLE :: DMIX
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ICMIX
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LFREE
      REAL, PARAMETER :: zeps = 1.E-4
!
      ALLOCATE(FDIS(NX,NY))
!
!     *** Domain Boundaries
      XMIN = MINVAL(XNEMO)
      XMAX = MAXVAL(XNEMO)
      YMIN = MINVAL(YNEMO)
      YMAX = MAXVAL(YNEMO)
      WRITE(*,*) "Ocean Model-Region: Lon ", XMIN,' - ', XMAX 
      WRITE(*,*) "                    Lat ", YMIN,' - ', YMAX 
      WRITE(*,*)        "First - Last Lat ", YNEMO(1,1),' - ', YNEMO(NX,NY) 
!
!     *** HD model resolution
      HDRES = FLON(2,1) - FLON(1,1)
      WRITE(*,*) '         Resolution HD model: ', HDRES, ' degree'
      DEGMAX=HDRES      ! 1 Gridbox outside ocean domain
      IF (ABS(HDRES-0.5) .LT.zeps) THEN     ! 0.5 degree
        DISMAX1=25000.    ! 25 km for primary mask
        DISMAX2=100000.   ! 100 km for secondary mask (default for only 1 mask)
        IF (IOCEAN.EQ.5) DISMAX2=400000.    ! ICON ocean 60km-coast --> 200 km 
                                            ! For all bgc inflows, 1000 km necessary
      ELSE IF (ABS(HDRES-0.5/6.) .LT.zeps) THEN    ! 5 Min.
        DISMAX1=4000.     ! 4 km for primary mask
        DISMAX2=16000.    ! 16 km for secondary mask (default for only 1 mask)
        IF (IOCEAN.EQ.1) DISMAX2=200000.    ! NEMO ocean coast too smooth --> 200 km 
        IF (IOCEAN.EQ.4) DISMAX2=100000.    ! ECOSMO-10 km ocean coast very smooth --> 200 km 
      ELSE IF (ABS(HDRES-0.0625) .LT.zeps) THEN          ! mHm
        DISMAX1=4000.     ! 4 km for primary mask
        DISMAX2=16000.    ! 16 km for secondary mask (default for only 1 mask)
        DEGMAX=0.5
        IF (IOCEAN.EQ.1) DISMAX2=200000.    ! NEMO ocean coast too smooth --> 200 km 
        IF (IOCEAN.EQ.2) DISMAX2=100000.     ! UFZ data contain sinks that should be 
                                            !     connected to the coast 
        IF (IOCEAN.EQ.4) DISMAX2=100000.    ! ECOSMO-10 km ocean coast very smooth --> 200 km 
      ENDIF
      WRITE(*,*) "Maximum distances: Primary = ", DISMAX1 
      WRITE(*,*) "                 Secondary = ", DISMAX2 
      WRITE(*,*) "                  Boundary = ", DEGMAX 
!
!     *** Note that the NEMO mouth masks defines mouth points as Zero.
      ALLOCATE(LMASK2(NX,NY))
      LMASK2(:,:) = .FALSE.
      WHERE(MASK_NEMO2.EQ.0)        
         LMASK2=.true.
      END WHERE
      IF (NR_MASK.EQ.2) THEN
        ALLOCATE(LMASK1(NX,NY))
        LMASK1(:,:) = .FALSE.
        WHERE(MASK_NEMO1.EQ.0)        
          LMASK1=.true.
        END WHERE
      ENDIF
!
!     ******** Searching for nearest NEMO mouth point
      INDEXX(:,:) = 0
      INDEXY(:,:) = 0
      MASK_MOU(:, :) = NINT(XMISS)
      MASK_ON_OCEAN(:,:)=0
      NMES=0
      NMOU=0
      NOMO = 0
      DO JB = 1, NB
      DO JL = 1, NL
      IF ( ABS(FMOUTH(JL,JB)-1).LE.0.1 ) THEN
        IF (FLON(JL,JB).GE.XMIN-DEGMAX .AND. FLON(JL,JB).LE.XMAX+DEGMAX .AND. &
            FLAT(JL,JB).GE.YMIN-DEGMAX .AND. FLAT(JL,JB).LE.YMAX+DEGMAX) THEN
          NMOU = NMOU+1
          MASK_MOU(JL, JB) = 0
          LFOUND = .FALSE.
          IF (NR_MASK.EQ.2) THEN
  	    CALL DISTANCE(FLON(JL,JB), FLAT(JL,JB), NX, NY, LMASK1,  &
                 XNEMO, YNEMO, FDIS)
            ICOORD = MINLOC(FDIS, MASK=LMASK1)
            DMIN = MINVAL(FDIS, MASK=LMASK1)
            IF (DMIN.LE.DISMAX1) LFOUND = .TRUE.
            IF (ICOORD(2).GT.730 .AND. ICOORD(1).GT.700) THEN
              WRITE(*, '(A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                  'M1: - Lon/Lat: ', JL, JB, ' = ', FLON(JL,JB), FLAT(JL,JB), &
                  ' Distance HD-Oceanmodel mouth = ', DMIN/1000., ' km at ', ICOORD
              WRITE(*,*) '     -> Lon = ', XNEMO(ICOORD(1), ICOORD(2)),  ' Lat = ', YNEMO(ICOORD(1),ICOORD(2))
            ENDIF
          ENDIF
          IF (.NOT. LFOUND) THEN

  	    CALL DISTANCE(FLON(JL,JB), FLAT(JL,JB), NX, NY, LMASK2,  &
                 XNEMO, YNEMO, FDIS)
            ICOORD = MINLOC(FDIS, MASK=LMASK2)
            DMIN = MINVAL(FDIS, MASK=LMASK2)
            IF (DMIN.LE.DISMAX2) LFOUND = .TRUE.
            IF (ICOORD(2).GT.730 .AND. ICOORD(1).GT.700) THEN
              WRITE(*, '(A,2(X,I4),A, 2(X,F7.2), A,X,F6.2,A,2(X,I4) )') &
                  'M2 - Lon/Lat: ', JL, JB, ' = ', FLON(JL,JB), FLAT(JL,JB),  &
                  ' Distance HD-Oceanmodel mouth = ', DMIN/1000., ' km at ', ICOORD
              WRITE(*,*) '     -> Lon = ', XNEMO(ICOORD(1), ICOORD(2)),  ' Lat = ', YNEMO(ICOORD(1),ICOORD(2))
            ENDIF
          ENDIF

          IF (LFOUND) THEN
            MASK_MOU(JL, JB) = 1
            NMES = NMES+1
            INDEXX(JL,JB) = ICOORD(1)
            INDEXY(JL,JB) = ICOORD(2)
            MASK_ON_OCEAN(ICOORD(1), ICOORD(2)) = MASK_ON_OCEAN(ICOORD(1), ICOORD(2)) + 1
            IF (IQUE.EQ.1) THEN
              WRITE(*,*) NMES, '. Distance HD-Oceanmodel mouth = ', DMIN/1000., ' km at ',  ICOORD
              WRITE(*,*) '     -> Lon = ', XNEMO(ICOORD(1), ICOORD(2)),  ' Lat = ', YNEMO(ICOORD(1),ICOORD(2))
            ENDIF
          ELSE
            NOMO = NOMO + 1
          ENDIF
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      WRITE(*,*) NMOU, ' HD Points in ',TRIM(CONAME),' region --> Converted to ', NMES,  &
                 ' ', TRIM(CONAME),' mouth points'
      NMIX = MAXVAL(MASK_ON_OCEAN)
      WRITE(*,*) 'Maximum number of HD mouths pointing towards a single ', TRIM(CONAME),' mouth: ', NMIX
      WRITE(*,*) 'No ', TRIM(CONAME), ' mouth was found for ', NOMO, ' HD mouths'
!
!     *** Indicate HD mouths that point towards a Ocean model mouth with more than one input  
      DO JB = 1, NB
      DO JL = 1, NL
      IF ( MASK_MOU(JL,JB).EQ.1 ) THEN
        IF (MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB)).GT.1)   &
            MASK_MOU(JL, JB) = MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB))
      ENDIF
      ENDDO
      ENDDO
!
      IF (ITYP.EQ.1) THEN
        IF (NMIX.GT.2) THEN
          WRITE(*,*) "****** ERROR - Cases with more more than 2 HD mouths: up to ", NMIX
          WRITE(*,*) "       pointing to the ",TRIM(CONAME)," NEMo box are NOT implemented"
          STOP       " --> Program stops right here!!!!"
        ENDIF
        ALLOCATE( DMIX(NMIX, NMIX) )
        ALLOCATE( ICMIX(NMIX, NMIX, 2) )
        ALLOCATE( LFREE(NX,NY) )
        WHERE (MASK_ON_OCEAN(:,:).EQ.1) LMASK2(:,:) = .FALSE.
        DO JB = 1, NB
        DO JL = 1, NL
        IF ( MASK_MOU(JL,JB).EQ.2 ) THEN
          JX = INDEXX(JL,JB)
          JY = INDEXY(JL,JB)
          NACT = MASK_ON_OCEAN(JX,JY)
          XLON = FLON(JL,JB)
          XLAT = FLAT(JL,JB)
          LFREE(:,:) = LMASK2(:,:)
          I=0

          WRITE(*, '(A,2(X,I4),A,2(X,F7.2) ,A,2(X,F7.2),A, 2(X,I4) )') &
                   'Mixpoint - Lon/Lat: ', JL,JB, ' = ', XLON, XLAT, ' pointing towards ', &
                    XNEMO(JX,JY), YNEMO(JX,JY), ' --> JX, JY: ', JX, JY
!
!         *** Distances 1, NACT
          DO WHILE(I.LT.NACT)
            I=I+1
	    CALL DISTANCE(XLON, XLAT, NX, NY, LFREE,  &
               XNEMO, YNEMO, FDIS)
            ICMIX(1,I,:) = MINLOC(FDIS, MASK=LFREE)
            DMIX(1,I) = MINVAL(FDIS, MASK=LFREE)
            LFREE(ICMIX(1,I,1),ICMIX(1,I,2)) = .FALSE.
          ENDDO
         
          IF (IQUE.EQ.1) WRITE(*,*) 'Dist1 ', DMIX(1,1:NACT) 
!
!         *** HD mouths poiting to the same Ocean model mouth
          IBOX=2
          DO IB = 1, NB
          DO IL = 1, NL
          IF (MASK_MOU(IL,IB).EQ.2 .AND. (IL.NE.JL .OR. IB.NE.JB) .AND.  &
             INDEXX(IL,IB).EQ.JX .AND. INDEXY(IL,IB).EQ.JY ) THEN

            IF (IQUE.EQ.1) WRITE(*,*) 'IL, IB ', IL, IB,  ' -> ', INDEXX(IL,IB), INDEXY(IL,IB)
            XLON = FLON(IL,IB)
            XLAT = FLAT(IL,IB)
            LFREE(:,:) = LMASK2(:,:)
            I=0
!
!           *** Distances IBOX, NACT
            DO WHILE(I.LT.NACT)
              I=I+1
	      CALL DISTANCE(XLON, XLAT, NX, NY, LFREE,  &
                 XNEMO, YNEMO, FDIS)
              ICMIX(IBOX,I,:) = MINLOC(FDIS, MASK=LFREE)
              DMIX(IBOX,I) = MINVAL(FDIS, MASK=LFREE)
              LFREE(ICMIX(IBOX,I,1),ICMIX(IBOX,I,2)) = .FALSE.
            ENDDO

            IF (IQUE.EQ.1) WRITE(*,*) 'Dist2 ', DMIX(2,1:NACT) 
!
!           *** Only NMIX = NACT = 2 can be dealt with
            IF (DMIX(1,1).LE.DMIX(2,1)) THEN           ! Box 1 closer than 2
              MASK_MOU(JL,JB) = 1           
              IF (DMIX(2,2).LE.DISMAX2) THEN
                INDEXX(IL,IB) = ICMIX(2,2,1)
                INDEXY(IL,IB) = ICMIX(2,2,2)
                MASK_ON_OCEAN(ICMIX(2,2,1), ICMIX(2,2,2)) = 1
                LMASK2(ICMIX(2,2,1), ICMIX(2,2,2)) = .FALSE.
!!                IF (IQUE.EQ.1) THEN
                WRITE(*,*) 'New Distance HD-NEMO mouth = ', DMIX(2,2)/1000., ' km at ',  ICMIX(2,2,:)
                IF (IQUE.EQ.1) WRITE(*,*) '     -> Lon = ', XNEMO(ICMIX(2,2,1), ICMIX(2,2,2)),  ' Lat = ', YNEMO(ICMIX(2,2,1), ICMIX(2,2,2))
!!                ENDIF
                MASK_MOU(IL,IB) = 1
              ELSE
                MASK_MOU(IL,IB) = 0
                INDEXX(IL,IB) = 0
                INDEXY(IL,IB) = 0
                WRITE(*,*) 'Second closest NEMO mouth is out of range: ', DMIX(2,2)/1000., ' km at ',  ICMIX(2,2,:)
              ENDIF
            ELSE                                       ! Box 2 closer than 1
              MASK_MOU(IL,IB) = 1           
              IF (DMIX(1,2).LE.DISMAX2) THEN
                INDEXX(JL,JB) = ICMIX(1,2,1)
                INDEXY(JL,JB) = ICMIX(1,2,2)
                MASK_ON_OCEAN(ICMIX(1,2,1), ICMIX(1,2,2)) = 1
                LMASK2(ICMIX(1,2,1), ICMIX(1,2,2)) = .FALSE.
!!                IF (IQUE.EQ.1) THEN
                WRITE(*,*) 'New Distance HD-NEMO mouth = ', DMIX(1,2)/1000., ' km at ',  ICMIX(1,2,:)
                IF (IQUE.EQ.1) WRITE(*,*) '     -> Lon = ', XNEMO(ICMIX(1,2,1), ICMIX(1,2,2)),  ' Lat = ', YNEMO(ICMIX(1,2,1), ICMIX(1,2,2))
!!                ENDIF
                MASK_MOU(JL,JB) = 1
              ELSE
                MASK_MOU(JL,JB) = 0
                INDEXX(JL,JB) = 0
                INDEXY(JL,JB) = 0
                WRITE(*,*) 'Second closest NEMO mouth is out of range: ', DMIX(1,2)/1000., ' km at ',  ICMIX(1,2,:)
              ENDIF
            ENDIF                           
!
!           *** Mixed problem solved for JX,JY
            MASK_ON_OCEAN(JX, JY) = 1
            LMASK2(JX, JY) = .FALSE.
            EXIT
          ENDIF
          ENDDO
          ENDDO
        ENDIF
        ENDDO
        ENDDO
        DEALLOCATE(DMIX)
        DEALLOCATE(ICMIX)
        DEALLOCATE(LFREE)
      ENDIF

      DEALLOCATE(FDIS)
      IF (NR_MASK.EQ.2) DEALLOCATE(LMASK1)
      DEALLOCATE(LMASK2)
      END SUBROUTINE FINDNEMO_FROMHD


      SUBROUTINE CORROCEAN(IOCEAN, NL, NB, MASK_MOU, INDEXX, INDEXY, NX, NY, MASK_ON_OCEAN)
!
!     ***  Correction of specific grid boxes depending on the ocean model grid.
!     ***  1) Removal of inflow points, e.g. near the ocean model domain boundary 
!     ***  2) Choosing a dedicated inflow point, e.g. for the Elbe
! 
      INTEGER, INTENT(in) :: IOCEAN    ! Ocean model setup ID
      INTEGER, INTENT(in) :: NL, NB
      INTEGER, DIMENSION(NL,NB), INTENT(inout) :: MASK_MOU  ! Mask with HD mouths with a valid 
                                                            ! nearest mouth on ocean grid
      INTEGER, DIMENSION(NL,NB), INTENT(inout) :: INDEXX    ! x-Indices of nearest ocean grid mouth points
      INTEGER, DIMENSION(NL,NB), INTENT(inout) :: INDEXY    ! y-Indices of nearest ocean grid mouth points

      INTEGER, INTENT(in) :: NX, NY
      INTEGER, DIMENSION(NX,NY), INTENT(inout) :: MASK_ON_OCEAN  ! Mask with associated ocean mouth points

      INTEGER :: JL, JB, I, NDUM
      INTEGER :: NREM                              ! Number of mouth points to be removed
      INTEGER, DIMENSION(:), ALLOCATABLE :: IREMX  ! X/Lon ocean grid coordinate to be removed
      INTEGER, DIMENSION(:), ALLOCATABLE :: IREMY  ! Y/Lat ocean grid coordinate to be removed
      INTEGER :: NCORR                             ! Number of mouth points to be corrected
      INTEGER, DIMENSION(:), ALLOCATABLE :: IX     ! Longitude on discharge grid with changing target
      INTEGER, DIMENSION(:), ALLOCATABLE :: IY     ! Latitude on discharge grid with changing target
      INTEGER, DIMENSION(:), ALLOCATABLE :: INEWX  ! New target longitude on ocean grid
      INTEGER, DIMENSION(:), ALLOCATABLE :: INEWY  ! New target latitude on ocean grid
!
!     ***

      SELECT CASE (IOCEAN)
        CASE (2)
          NREM=2
          ALLOCATE(IREMX(NREM)) ; ALLOCATE(IREMY(NREM))
          IREMX(1) = 281  ; IREMY(1) = 334
          IREMX(2) = 294  ; IREMY(2) = 341
        CASE (4)
          NCORR=1
          ALLOCATE(IX(NCORR)) ; ALLOCATE(IY(NCORR))
          ALLOCATE(INEWX(NCORR)) ; ALLOCATE(INEWY(NCORR))
          IF (NL.EQ.832 .AND. NB.EQ.592) THEN        ! mRm source
            IX(1) = 334  ; IY(1) = 296         ! Elbe
          ELSE IF (NL.EQ.960 .AND. NB.EQ.540) THEN   ! HD Euro 5 Min. source
!
!           *** Removal HD
!            NREM=8
!            ALLOCATE(IREMX(NREM)) ; ALLOCATE(IREMY(NREM))
!            IREMX(1) = 18  ; IREMY(1) = 66   ! Orkney
!            IREMX(2) = 57  ; IREMY(2) = 66
!            IREMX(3) = 59  ; IREMY(3) = 69
!            IREMX(4) = 9   ; IREMY(4) = 73
!            IREMX(5) = 8   ; IREMY(5) = 75
!            IREMX(6) = 16  ; IREMY(6) = 156
!            IREMX(7) = 19  ; IREMY(7) = 166
!            IREMX(8) = 20  ; IREMY(8) = 167
!
!           *** Correction
            IX(1) = 251  ; IY(1) = 222         ! Elbe
          ENDIF
!!          INEWX(1) = 77  ; INEWY(1) = 122    ! Elbe
          INEWX(1) = 78  ; INEWY(1) = 120    ! Elbe
        CASE DEFAULT
           RETURN
      END SELECT
!
!     *** Removal of ocean mouths
      IF (NREM.GE.1) THEN
      DO I = 1, NREM
        NDUM = MASK_ON_OCEAN(IREMX(I), IREMY(I))
        DO JB=1, NB
        DO JL=1, NL
        IF (INDEXX(JL,JB).EQ.IREMX(I) .AND. INDEXY(JL,JB).EQ.IREMY(I)) THEN
          INDEXX(JL,JB) = 0
          INDEXY(JL,JB) = 0
          MASK_MOU(JL,JB) = 0
          MASK_ON_OCEAN(IREMX(I), IREMY(I)) = MASK_ON_OCEAN(IREMX(I), IREMY(I)) - 1
        ENDIF
        ENDDO
        ENDDO
        WRITE(*,*) I, '. point reduced from ', NDUM, ' to ', MASK_ON_OCEAN(IREMX(I), IREMY(I))
      ENDDO
      DEALLOCATE(IREMX)
      DEALLOCATE(IREMY)
      ENDIF
!
!     *** Correction of target points
      IF (NCORR.GE.1) THEN
        DO I = 1, NCORR
          JL = IX(I) ; JB = IY(I)
          MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB)) = MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB)) - 1
          INDEXX(JL,JB) = INEWX(I)
          INDEXY(JL,JB) = INEWY(I)
          MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB)) = MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB)) + 1
          MASK_MOU(JL,JB) = MASK_ON_OCEAN(INDEXX(JL,JB), INDEXY(JL,JB))
          WRITE(*,*) I, '. point at ', JL,JB, ' corrected towards ', INDEXX(JL,JB), INDEXY(JL,JB)
        ENDDO
        DEALLOCATE(IX) ; DEALLOCATE(IY)
        DEALLOCATE(INEWX) ; DEALLOCATE(INEWY)
      ENDIF

      END SUBROUTINE CORROCEAN
!	  
    END PROGRAM HDTONEMO
