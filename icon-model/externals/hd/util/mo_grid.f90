! mo_grid.f90 - Utilities for reading and handling gridded data stored in Netcdf files
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

MODULE mo_grid

  !
  ! Authors:
  !
  ! S. Hagemann, HZG-IfK, November 2017- November 2018, original source

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: domain, indexcalc, areacalc, read_grid_info, index_nn, distance, &
            calc_area_reg, calc_area_irreg, model, read_coordinates, generate_grid_coordinates

  !-----------------------------------------------------------------------------
  TYPE domain
    INTEGER :: nlon                         ! Longitude dimension
    INTEGER :: nlat                         ! Latitude dimension
    DOUBLE PRECISION  :: resolution         ! = resolution_lon = resolution_lat 
    DOUBLE PRECISION  :: resolution_lon     ! Longitudinal resolution
    DOUBLE PRECISION  :: resolution_lat     ! Latitudinal resolution
    DOUBLE PRECISION  :: origin_lon, origin_lat   ! NW coordinates of grid (1,1), not centre.
    DOUBLE PRECISION  :: end_lon, end_lat         ! SE corner coordinates of grid not centre.
    INTEGER :: kshift_n, kshift_s                 ! Box shift of grid from N and S pole
    INTEGER :: kshift_w, kshift_e                 ! Box shift of grid W and E date line

  END TYPE domain

  TYPE model
    INTEGER :: nlon                                        ! Longitude dimension
    INTEGER :: nlat                                        ! Latitude dimension
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: value          ! Mask value, usually 0 or 1.
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: xlon ! Longitude coordinate
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: xlat ! Latitude coordinate
    CHARACTER (LEN=20) :: name          ! Mask name or name of respective model
    INTEGER :: idim                     ! No. of coordinate dimensions (1 or 2)
  END TYPE model


  INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  REAL(dp), PARAMETER :: PI = 3.141592653589793_dp
  REAL(dp), PARAMETER :: rerde = 6371229._dp          ! Referenz-Mittelwert des Erdradius, as in ICON
  REAL(dp), PARAMETER :: zeps=1.E-10_dp 

CONTAINS

!   *******************************************************************************
    SUBROUTINE indexcalc(grid, iglob, xlon, xlat, ilon, ilat)
!   *******************************************************************************

      TYPE(domain),INTENT(in)      :: grid

      INTEGER, INTENT(in)          :: iglob                    ! provide indices in global grid (0/1=no/yes)
      REAL(dp), INTENT(in)         :: xlon, xlat
      INTEGER, INTENT(out)         :: ilon, ilat

      REAL(dp) :: orglon, orglat

      IF (iglob.eq.1) THEN
        orglon = -180._dp
        orglat = 90._dp
      ELSE
        orglon = grid%origin_lon
        orglat = grid%origin_lat
      ENDIF

      ilon = FLOOR( (xlon - orglon)/grid%resolution_lon + 1.0001)
      IF (ilon.LE.0) ilon = ilon+grid%nlon
      IF (ilon.GT.grid%nlon) ilon = ilon-grid%nlon
      ilat = FLOOR(1.0001 + (orglat-xlat)/grid%resolution_lat )

!!      WRITE(*,*) 'indexcalc: ',  xlon, xlat , ' --> ', ilon, ilat

    END SUBROUTINE indexcalc

!   *******************************************************************************
    SUBROUTINE areacalc(grid, vlat, area, dlat, dlon)
!   *******************************************************************************

!      *** Routine calculates area of each grid cell in a regular grid 
!          with a constant resolution in degree

      TYPE(domain),INTENT(in)      :: grid

      DOUBLE PRECISION, DIMENSION(grid%nlat), INTENT(in)  :: vlat  ! latitude coordinates (centre)
      DOUBLE PRECISION, DIMENSION(grid%nlon, grid%nlat), INTENT(out) :: area  ! Grid cell area in m2
      DOUBLE PRECISION, INTENT(out)                       :: dlat  ! Latitudinal distance in m
      DOUBLE PRECISION, DIMENSION(grid%nlat), INTENT(out) :: dlon  ! Longitudinal distance in m

      DOUBLE PRECISION :: pifak

      INTEGER :: jb

      pifak = PI/180._dp
      dlat = grid%resolution_lat * pifak * rerde

      DO jb=1, grid%nlat

        dlon(jb) = grid%resolution_lon * pifak * cos(vlat(jb)*pifak) * rerde
        area(:, jb) = dlat * rerde *  &
            ABS( SIN((vlat(jb)+grid%resolution_lat/2._dp) * pifak) - SIN((vlat(jb)-grid%resolution_lat/2._dp)*pifak) )

        IF (jb.LE.5 .OR. jb.GE.grid%nlat-4 .OR. jb.EQ.grid%nlat/2  &
      	     .OR. jb.EQ.INT(grid%nlat/2) +1) THEN
           WRITE(*, '(I5,X,F9.4,X,F11.6,X, F9.4,A4,F11.2,A4)')  &
                jb, vlat(jb), COS( vlat(jb)*PI/180._dp ), dlon(jb)/1000.,' km ', area(1,jb)*1.E-6, ' km2'
        ENDIF		  
      ENDDO

    END SUBROUTINE areacalc

!   *******************************************************************************
    SUBROUTINE calc_area_reg(nlon,nlat, vlon, vlat, farea)
!   *******************************************************************************

!      *** Routine calculates area of each grid cell in a regular grid where lon and lat are
!          stored in 1D arrays
!          Area of gridbox with borders lat1, lat2, lon1, lon2 (unit degree) is:
!          A = (pi/180) * R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|

      INTEGER, INTENT(in) :: nlon
      INTEGER, INTENT(in) :: nlat
      DOUBLE PRECISION, DIMENSION(nlon), INTENT(in)  :: vlon  ! longitude coordinates (centre)
      DOUBLE PRECISION, DIMENSION(nlat), INTENT(in)  :: vlat  ! latitude coordinates (centre)
      DOUBLE PRECISION, DIMENSION(nlon, nlat), INTENT(out) :: farea  ! Grid cell area in m2

      DOUBLE PRECISION :: pifak

      DOUBLE PRECISION :: dlat, dlon, vlat1, vlat2
      INTEGER :: jl, jb

      pifak = PI/180._dp

      DO jb=1, nlat
        IF (jb.eq.1) THEN
          vlat2 = (vlat(jb)+vlat(jb+1))/2._dp
          vlat1 = vlat(jb) - (vlat2 - vlat(jb))
        ELSE IF (jb.eq.nlat) THEN
          vlat1 = (vlat(jb)+vlat(jb-1))/2._dp
          vlat2 = vlat(jb) - (vlat1 - vlat(jb)) 
        ELSE 
          vlat1 = (vlat(jb)+vlat(jb-1))/2._dp
          vlat2 = (vlat(jb)+vlat(jb+1))/2._dp
        ENDIF
!!         dlat = ABS(vlat1 - vlat2) * pifak * rerde
        DO jl=1, nlon
          IF (jl.eq.1) THEN
            dlon = ABS((vlon(jl)+vlon(jl+1))/2._dp - vlon(jl)) * 2._dp * pifak * rerde
          ELSE IF (jl.eq.nlon) THEN
            dlon = ABS((vlon(jl)+vlon(jl-1))/2._dp - vlon(jl)) * 2._dp * pifak * rerde
          ELSE 
            dlon = ABS((vlon(jl)+vlon(jl-1))/2._dp - (vlon(jl)+vlon(jl+1))/2._dp) * pifak * rerde
          ENDIF

          farea(jl, jb) = rerde * dlon *  &
             ABS( SIN(vlat1 * pifak) - SIN(vlat2*pifak) )
        ENDDO
      ENDDO

    END SUBROUTINE calc_area_reg

!   *******************************************************************************
    SUBROUTINE calc_area_irreg(nlon,nlat, flon, flat, farea)
!   *******************************************************************************

!      *** Routine calculates area of each grid cell in a regular grid where lon and lat are
!          stored in 2D arrays
!          Area of gridbox with borders lat1, lat2, lon1, lon2 (unit degree) is:
!          A = (pi/180) * R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|

      INTEGER, INTENT(in) :: nlon
      INTEGER, INTENT(in) :: nlat
      DOUBLE PRECISION, DIMENSION(nlon, nlat), INTENT(in)  :: flon   ! longitude coordinates (centre)
      DOUBLE PRECISION, DIMENSION(nlon, nlat), INTENT(in)  :: flat   ! latitude coordinates (centre)
      DOUBLE PRECISION, DIMENSION(nlon, nlat), INTENT(out) :: farea  ! Grid cell area in m2

      DOUBLE PRECISION :: pifak

      DOUBLE PRECISION :: dlat, dlon, vlat1, vlat2
      INTEGER :: jl, jb

      pifak = PI/180._dp

      DO jb=1, nlat
      DO jl=1, nlon
        IF (jb.eq.1) THEN
          vlat2 = (flat(jl,jb)+flat(jl,jb+1))/2._dp
          vlat1 = flat(jl,jb) - (vlat2 - flat(jl,jb)) / 2._dp
        ELSE IF (jb.eq.nlat) THEN
          vlat1 = (flat(jl,jb)+flat(jl,jb-1))/2._dp
          vlat2 = flat(jl,jb) - (vlat1 - flat(jl,jb)) / 2._dp
        ELSE 
          vlat1 = (flat(jl,jb)+flat(jl,jb-1))/2._dp
          vlat2 = (flat(jl,jb)+flat(jl,jb+1))/2._dp
        ENDIF
!!         dlat = ABS(vlat1 - vlat2) * pifak * rerde
        IF (jl.eq.1) THEN
          dlon = ABS((flon(jl,jb)+flon(jl+1,jb))/2._dp - flon(jl,jb)) * 2._dp * pifak * rerde
        ELSE IF (jl.eq.nlon) THEN
          dlon = ABS((flon(jl,jb)+flon(jl-1,jb))/2._dp - flon(jl,jb)) * 2._dp * pifak * rerde
        ELSE 
          dlon = ABS((flon(jl,jb)+flon(jl-1,jb))/2._dp - (flon(jl,jb)+flon(jl+1,jb))/2._dp) * pifak * rerde
        ENDIF

        farea(jl, jb) = rerde * dlon *  &
             ABS( SIN(vlat1 * pifak) - SIN(vlat2*pifak) )
      ENDDO
      ENDDO

    END SUBROUTINE calc_area_irreg

!   *******************************************************************************
    SUBROUTINE read_grid_info(dnam, grid, xmiss, clon,clat, cdimlon,cdimlat, idim)
!   *******************************************************************************
      use netcdf
!
!     *** For 1-D mask arrays (such as ICON), cnb shoul be: cnb='-'  
! 	  
      CHARACTER (LEN=*), INTENT(in)  :: dnam    ! File name with grid info
      TYPE(domain),INTENT(out)       :: grid    ! grid info structure
      REAL(dp), INTENT(in)           :: xmiss   ! Missing value 
!     *** NetCDF names
      CHARACTER (LEN=*), INTENT(in)  :: clon    ! Array name of longitudinal coordinates
      CHARACTER (LEN=*), INTENT(in)  :: clat    ! Array name of latitudinal coordinates
      CHARACTER (LEN=*), INTENT(inout), OPTIONAL :: cdimlon ! Name of longitudinal dimension 
      CHARACTER (LEN=*), INTENT(inout), OPTIONAL :: cdimlat ! Name of latitudinal dimension
      INTEGER, INTENT(inout), OPTIONAL  :: idim    ! Number of coordinate dimensions: 1 or 2 (default)
      LOGICAL                        :: LOG2    ! True for 2 dimensions
      LOGICAL                        :: log_cdim ! Are dimension names provided?

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vlon
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: vlat
      CHARACTER (LEN=20)  :: cnl
      CHARACTER (LEN=20)  :: cnb
      CHARACTER (LEN=20)  :: cunit
!
      INTEGER, PARAMETER  :: nltest = 5
      CHARACTER (LEN=20), DIMENSION(nltest) :: clontest
      CHARACTER (LEN=20), DIMENSION(nltest) :: clattest
      INTEGER :: i, idimtest, numdims
!
      DATA clontest / 'lon', 'x', 'nx', 'nlon', 'ncells' /
      DATA clattest / 'lat', 'y', 'ny', 'nlat', '-' /
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, lat_dimid, lon_dimid, varid(2), dimid, ierr2
!
!     **********************************************
      log_cdim = .FALSE.
      IF (PRESENT(cdimlon)) THEN
        cnl = TRIM(cdimlon)
        log_cdim = .TRUE.
      ELSE
        cnl = 'lon'
      ENDIF
      IF (PRESENT(cdimlat)) THEN
        cnb = cdimlat
      ELSE
        cnb = 'lat'
      ENDIF
      idimtest = 2
!				 
!     ******** Open and read dimensions
      ierr = nf90_open(dnam, NF90_NOWRITE, ncid)
      ierr = nf90_inq_dimid(ncid,TRIM(cnl),lon_dimid)
      ierr2 = 1
      DO i=1, nltest
        IF (ierr.NE.NF90_NOERR) THEN
          cnl = clontest(i)
          cnb = clattest(i)
          ierr = nf90_inq_dimid(ncid,TRIM(cnl),lon_dimid)
          IF (i.EQ.nltest) idimtest = 1
        ELSE
          exit
        ENDIF
      ENDDO

      IF (ierr.NE.NF90_NOERR) THEN
        WRITE(*,*) 'First coordinate dimension not found. Tests for ', clontest
        STOP 'ABBRUCH'
      ENDIF 
      ierr = nf90_inquire_dimension(ncid,lon_dimid,len=grid%nlon)
      IF (log_cdim) cdimlon = cnl

      IF (idimtest.NE.1 .AND. TRIM(cnb).NE.'-') THEN
        ierr = nf90_inq_dimid(ncid,TRIM(cnb),lat_dimid)
        IF (ierr.NE.NF90_NOERR) THEN
          WRITE(*,*) 'Second coordinate dimension not found for cnb = ', TRIM(cnb)
          STOP 'ABBRUCH'
        ENDIF 
        ierr = nf90_inquire_dimension(ncid,lat_dimid,len=grid%nlat)
        WRITE(*,*) "Dimensions: NL =" , grid%nlon, "  NB = ", grid%nlat
        IF (log_cdim) cdimlat = cnb
      ELSE
        WRITE(*,*) "1-D Array dimension size =" , grid%nlon
        grid%nlat = 1
      ENDIF

      LOG2 = .TRUE.
      IF (PRESENT(idim)) THEN
        idimtest = idim
      ENDIF

!
!     *** Dimension of Lon/Lat arrays
      ierr = nf90_inq_varid(ncid,TRIM(clon),varid(1))
      IF (ierr.NE.NF90_NOERR) THEN
        WRITE(*,*) 'Variable ', TRIM(clon), ' not included in ', TRIM(dnam) 
        STOP 'ABBRUCH'
      ENDIF      
      ierr = nf90_inquire_variable(ncid, varid(1), ndims = numdims)
      IF (numdims.EQ.1) LOG2 = .FALSE.

      IF (LOG2) THEN
        ALLOCATE(vlon(grid%nlon, grid%nlat))
        ALLOCATE(vlat(grid%nlon, grid%nlat))
        WRITE(*,*) 'Grid has 2-dimensional coordinates'
      ELSE
        ALLOCATE(vlon(grid%nlon, 1))
        IF (grid%nlat.EQ.1) THEN
          ALLOCATE(vlat(grid%nlon, 1))
          WRITE(*,*) 'Grid has 1-D coordinates with', grid%nlon, ' longitudes'
        ELSE
          ALLOCATE(vlat(grid%nlat, 1))
          WRITE(*,*) 'Grid has 1-D coordinates with ', grid%nlat, ' latitudes'
        ENDIF 
      ENDIF
      ierr = nf90_get_var(ncid,varid(1),vlon)

      ierr = nf90_inq_varid(ncid,TRIM(clat),varid(2))
      IF (ierr.NE.NF90_NOERR) THEN
        WRITE(*,*) 'Variable ', TRIM(clat), ' not included in ', TRIM(dnam) 
        STOP 'ABBRUCH'
      ENDIF      
      ierr = nf90_get_var(ncid,varid(2),vlat)
!
!     Calculate coordinates in degree if original unit is radian
      ierr = nf90_inquire_attribute(ncid, varid(1), 'units')
      IF (ierr.EQ.NF90_NOERR) THEN
        ierr = nf90_get_att(ncid, varid(1), 'units', cunit)
        IF (TRIM(cunit).EQ.'radian') THEN
          vlon = vlon / PI * 180._dp
          WHERE (vlon.GT. 180._dp) vlon = vlon - 360._dp
          WHERE (vlon.LT.-180._dp) vlon = vlon + 360._dp
        ENDIF
      ENDIF
      ierr = nf90_inquire_attribute(ncid, varid(2), 'units')
      IF (ierr.EQ.NF90_NOERR) THEN
        ierr = nf90_get_att(ncid, varid(2), 'units', cunit)
        IF (TRIM(cunit).EQ.'radian') THEN
          vlat = vlat / PI * 180._dp
        ENDIF
      ENDIF

      ierr = nf90_close(ncid)
      WRITE(*,*) TRIM(dnam)," is closed."
!
!     Calculate grid info
      WRITE(*,*) "Longitude range: ", TRIM(clon), MINVAL(vlon, ABS(vlon-xmiss).GT.zeps), &
                  MAXVAL(vlon, ABS(vlon-xmiss).GT.zeps)
      WRITE(*,*) "Latitude range vlat:", TRIM(clat), MINVAL(vlat, ABS(vlat-xmiss).GT.zeps), &
                  MAXVAL(vlat, ABS(vlat-xmiss).GT.zeps)
      IF (LOG2) THEN
        WRITE(*,*) 'Grid centre: ', vlon(grid%nlon/2, grid%nlat/2), vlat(grid%nlon/2, grid%nlat/2)
      ENDIF
      grid%resolution_lon = (MAXVAL(vlon) - MINVAL(vlon)) / (grid%nlon-1)
      grid%resolution_lat = (MAXVAL(vlat) - MINVAL(vlat)) / (grid%nlat-1)
      IF (ABS(grid%resolution_lon-grid%resolution_lat).LT.zeps) THEN
        grid%resolution = grid%resolution_lon
        WRITE(*,*) "Res. ", grid%resolution, '°'
      ELSE
        grid%resolution = -1._dp
        WRITE(*,*) "Lon Res. ", grid%resolution_lon, "°  Lat Res. ", grid%resolution_lat, "°"
      ENDIF

      grid%origin_lon = MINVAL(vlon) - grid%resolution_lon/2._dp
      grid%origin_lat = MAXVAL(vlat) + grid%resolution_lat/2._dp
      WRITE(*,*) 'Origin ', grid%origin_lon, '°E, ', grid%origin_lat, '°N'
      grid%end_lon = MAXVAL(vlon) + grid%resolution_lon/2._dp
      grid%end_lat = MINVAL(vlat) - grid%resolution_lat/2._dp
      WRITE(*,*) "SE corner: ", grid%end_lon, '°E, ', grid%end_lat, '°S'

      grid%kshift_n = NINT((90._dp - grid%origin_lat) / grid%resolution_lat)
      grid%kshift_s = NINT((90._dp + MINVAL(vlat) - grid%resolution_lat/2._dp) / grid%resolution_lat)
      WRITE(*,*) "N and S borders are shifted from the poles by ", grid%kshift_n, ' and ', grid%kshift_s
      grid%kshift_w = NINT((grid%origin_lon + 180._dp ) / grid%resolution_lon)
      grid%kshift_e = NINT((180._dp - (MAXVAL(vlon) + grid%resolution_lon/2._dp)) / grid%resolution_lon )
      WRITE(*,*) "W and E borders are shifted from the dateline by ", grid%kshift_w, ' and ', grid%kshift_e
	  
      DEALLOCATE(vlon)
      DEALLOCATE(vlat)

    END SUBROUTINE read_grid_info

!   *******************************************************************************
    SUBROUTINE read_coordinates(dnam, msk, clon, clat)
!   *******************************************************************************
      use netcdf
	  
      CHARACTER (LEN=*), INTENT(in)  :: dnam    ! File name with grid info
      TYPE(model), INTENT(inout)      :: msk     ! grid info structure
      CHARACTER (LEN=*), INTENT(in)  :: clon    ! Array name of longitudinal coordinates
      CHARACTER (LEN=*), INTENT(in)  :: clat    ! Array name of latitudinal coordinates

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vlon
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vlat
      INTEGER :: jl, jb
      CHARACTER (LEN=20)  :: cunit
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, numdims, varid(2)
!				 
!     ******** Open file 
      ierr = nf90_open(dnam, NF90_NOWRITE, ncid)
!
!     *** Dimension of Lon/Lat arrays of HD model mouth array: currently 1
      ierr = nf90_inq_varid(ncid,TRIM(clon),varid(1))
      ierr = nf90_inquire_variable(ncid, varid(1), ndims = numdims)
      IF (numdims.EQ.2 .AND. msk%nlat.NE.1) THEN
        ierr = nf90_get_var(ncid,varid(1),msk%xlon)
        ierr = nf90_inq_varid(ncid,TRIM(clat),varid(2))
        ierr = nf90_get_var(ncid,varid(2),msk%xlat)
!
!       Calculate coordinates in degree if original unit is radian
        ierr = nf90_inquire_attribute(ncid, varid(1), 'units')
        IF (ierr.EQ.NF90_NOERR) THEN
          ierr = nf90_get_att(ncid, varid(1), 'units', cunit)
          IF (TRIM(cunit).EQ.'radian') THEN
            msk%xlon = msk%xlon / PI * 180._dp
            WHERE (msk%xlon.GT. 180._dp) msk%xlon = msk%xlon - 360._dp
            WHERE (msk%xlon.LT.-180._dp) msk%xlon = msk%xlon + 360._dp
          ENDIF
        ENDIF
        ierr = nf90_inquire_attribute(ncid, varid(2), 'units')
        IF (ierr.EQ.NF90_NOERR) THEN
          ierr = nf90_get_att(ncid, varid(2), 'units', cunit)
          IF (TRIM(cunit).EQ.'radian') THEN
            msk%xlat = msk%xlat / PI * 180._dp
          ENDIF
        ENDIF

      ELSE IF (numdims.EQ.1 .OR. msk%nlat.EQ.1) THEN
        ALLOCATE(vlon(msk%nlon))
        IF (msk%nlat.EQ.1) THEN
          ALLOCATE(vlat(msk%nlon))
        ELSE
          ALLOCATE(vlat(msk%nlat))
        ENDIF
        ierr = nf90_get_var(ncid,varid(1),vlon)
        ierr = nf90_inq_varid(ncid,TRIM(clat),varid(2))
        ierr = nf90_get_var(ncid,varid(2),vlat)
!
!       Calculate coordinates in degree if original unit is radian
        ierr = nf90_inquire_attribute(ncid, varid(1), 'units')
        IF (ierr.EQ.NF90_NOERR) THEN
          ierr = nf90_get_att(ncid, varid(1), 'units', cunit)
          IF (TRIM(cunit).EQ.'radian') THEN
            vlon = vlon / PI * 180._dp
            WHERE (vlon.GT. 180._dp) vlon = vlon - 360._dp
            WHERE (vlon.LT.-180._dp) vlon = vlon + 360._dp
          ENDIF
        ENDIF
        ierr = nf90_inquire_attribute(ncid, varid(2), 'units')
        IF (ierr.EQ.NF90_NOERR) THEN
          ierr = nf90_get_att(ncid, varid(2), 'units', cunit)
          IF (TRIM(cunit).EQ.'radian') THEN
            vlat = vlat / PI * 180._dp
          ENDIF
        ENDIF

        WRITE(*,'(A,A,X,F10.4,A,F10.4)') "READ COORDINATES - vlon: ", TRIM(clon), MINVAL(vlon),' - ',MAXVAL(vlon)
        IF (msk%nlat.EQ.1) THEN  ! e.g. ICON
          WRITE(*,'(A,A,X,F10.4,A,F10.4)') "READ COORDINATES - vlat: ", TRIM(clat), MINVAL(vlat),' - ',MAXVAL(vlat)
          DO jl=1, msk%nlon
            msk%xlon(:,1) = vlon(:)
            msk%xlat(:,1) = vlat(:)
          ENDDO
        ELSE
          WRITE(*,'(A,A,X,F10.4,A,F10.4)') "READ COORDINATES - vlat: ", TRIM(clat), MINVAL(vlat),' - ',MAXVAL(vlat)
          DO jl=1, msk%nlon
            msk%xlat(jl,:) = vlat(:)
          ENDDO
          DO jb=1, msk%nlat
            msk%xlon(:,jb) = vlon(:)
          ENDDO
        ENDIF
        DEALLOCATE(vlon)
        DEALLOCATE(vlat)
      ELSE
        STOP 'No. of coordinate dimension neither 1 nor 2 --> STOP!'
      ENDIF
      ierr = nf90_close(ncid)
      msk%idim = numdims
      WRITE(*,*) TRIM(dnam)," is closed in read_coordinates."

    END SUBROUTINE read_coordinates
!
!  *********************************************************************
   SUBROUTINE distance(xcoord, ycoord, nx, ny, lmask, xlon, xlat, fdist)
!     *********************************************************************
!
!     Calculates distance of xcoord,ycoord to xlon, xlat for all points 
!     that are true LMASK.
!

      DOUBLE PRECISION, INTENT(in) :: xcoord, ycoord
      INTEGER, INTENT(in) :: nx, ny
      LOGICAL, DIMENSION(nx,ny), INTENT(in) :: lmask
      DOUBLE PRECISION, DIMENSION(nx, ny), INTENT(in) :: xlon
      DOUBLE PRECISION, DIMENSION(nx, ny), INTENT(in) :: xlat
      DOUBLE PRECISION, DIMENSION(nx, ny), INTENT(out) :: fdist
      DOUBLE PRECISION, PARAMETER  ::   PIFAK = PI/180.

      DOUBLE PRECISION :: DLON, DLAT
      INTEGER :: IX, IY
!
      fdist(:,:) = -9999.
      DO IX=1, nx 
      DO IY=1, ny 
         IF (lmask(IX, IY)) THEN
	    ! Pay regard to sphere
           DLON = ABS(xlon(IX,IY)-xcoord)
           IF (DLON.GT.180.) DLON = 360. - ABS(xlon(IX,IY)) - ABS(xcoord)
           DLON = DCOS( (xlat(IX,IY)+ycoord) / 2. *PIFAK ) * RERDE * DLON * PIFAK
           DLAT = ABS(xlat(IX,IY)-ycoord) * PIFAK * RERDE
            
	   fdist(IX, IY) = SQRT(DLON*DLON + DLAT*DLAT)
	 ENDIF
      ENDDO	  
      ENDDO	  
!
   END SUBROUTINE distance

!   *******************************************************************************
    SUBROUTINE index_nn(nlon, nlat, xlon, xlat, lmask, nco, xcoord, ycoord, icox, icoy)
!   *******************************************************************************

      INTEGER, INTENT(in)             :: nlon     ! No. of longitudes in grid
      INTEGER, INTENT(in)             :: nlat     ! No. of latitudes in grid
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)  :: xlon  ! Longitudes of grid
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)  :: xlat  ! Latitudes of grid
      LOGICAL, DIMENSION(:, :), INTENT(in)      :: lmask     ! Mask of valid points where index should be searched.
      INTEGER, INTENT(in)             :: nco                 ! No. of given coordinates 
      DOUBLE PRECISION, DIMENSION(:), INTENT(in)  :: xcoord  ! Array of longitude coordinates
      DOUBLE PRECISION, DIMENSION(:), INTENT(in)  :: ycoord  ! Array of latitude coordinates
      INTEGER, DIMENSION(:), INTENT(out)  :: icox          ! Array of longitude indices nearest to the coordinates
      INTEGER, DIMENSION(:), INTENT(out)  :: icoy          ! Array of latitude indices nearest to the coordinates
!
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: fdist
      DOUBLE PRECISION                     :: dmin
      REAL :: xmin, xmax, ymin, ymax
      INTEGER  :: jl, jb, i
      INTEGER, DIMENSION(2) :: icoord
!
      xmin = MINVAL(xlon, lmask)
      xmax = MAXVAL(xlon, lmask)
      ymin = MINVAL(xlat, lmask)
      ymax = MAXVAL(xlat, lmask)
      ALLOCATE(fdist(nlon, nlat))
      icox(:) = 0
      icoy(:) = 0
!
      DO i=1, nco
      IF (xcoord(i).GE.xmin .AND. xcoord(i).LE.xmax .AND. &
          ycoord(i).GE.ymin .AND. ycoord(i).LE.ymax ) THEN

         CALL distance(xcoord(i), ycoord(i), nlon, nlat, lmask,    &
              xlon, xlat, fdist)
         dmin = MINVAL(fdist, lmask)
         icoord = MINLOC(fdist, lmask)
         icox(i) = icoord(1)
         icoy(i) = icoord(2)
         IF (icoord(1).LE.0 .OR. icoord(1).GT.nlon .OR. icoord(2).LE.0 .OR. icoord(2).GT. nlat) THEN
           WRITE(*,'(A,I4, 2F7.2, A, F8.2)') 'Index_nn: ', i, xcoord(i), ycoord(i), ' Dist:', fdist
           STOP
         ENDIF
      ENDIF
      ENDDO

      DEALLOCATE(fdist)

    END SUBROUTINE index_nn

!   *******************************************************************************
    SUBROUTINE generate_grid_coordinates(grid, msk)
!   *******************************************************************************

      TYPE(domain),INTENT(in)      :: grid
      TYPE(model), INTENT(out)     :: msk     ! grid info structure
!
      INTEGER :: jl,jb

      msk%nlon = grid%nlon
      msk%nlat = grid%nlat
!
      ALLOCATE(msk%xlon(msk%nlon,msk%nlat))
      ALLOCATE(msk%xlat(msk%nlon,msk%nlat))

      DO jb=1, msk%nlat
        msk%xlat(:,jb) = grid%origin_lat - (jb-0.5_dp) * grid%resolution_lat
      ENDDO
      DO jl=1, msk%nlon
        msk%xlon(jl,:) = grid%origin_lon + (jl-0.5_dp) * grid%resolution_lon 
      ENDDO

    END SUBROUTINE generate_grid_coordinates

END MODULE mo_grid

