! mo_grid.f90 - Access and handle grid information of HD model domain
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
  ! Modifications:
  !
  ! S. Hagemann, HZG-IfK, November 2017, original source
  !                       November 2018 - adapted for HD model        
  !

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: domain, indexcalc, areacalc, read_grid_info

  !-----------------------------------------------------------------------------
  TYPE domain
    INTEGER :: nlon                      ! Longitude dimension
    INTEGER :: nlat                      ! Latitude dimension
    REAL(dp) :: resolution
    REAL(dp) :: origin_lon, origin_lat   ! NW coordinates of grid (1,1), not centre.
    INTEGER :: kshift_n, kshift_s        ! Box shift of grid from N and S pole
    INTEGER :: kshift_w, kshift_e        ! Box shift of grid W and E date line

  END TYPE domain


CONTAINS

!   *******************************************************************************
    SUBROUTINE indexcalc(grid, iglob, xlon, xlat, ilon, ilat)
!   *******************************************************************************

      TYPE(domain),INTENT(in)      :: grid

      INTEGER, INTENT(in)          :: iglob                    ! provide indices in global grid (0/1=no/yes)
      REAL(dp), INTENT(in) :: xlon, xlat
      INTEGER, INTENT(out)         :: ilon, ilat

      REAL(dp) :: fskal, orglon, orglat

      IF (iglob.eq.1) THEN
        orglon = -180.
        orglat = 90.
      ELSE
        orglon = grid%origin_lon
        orglat = grid%origin_lat
      ENDIF
      fskal = 1./grid%resolution

      ilon = FLOOR( fskal*(xlon - orglon) + 1.0001)
      IF (ilon.LE.0) ilon = ilon+grid%nlon
      IF (ilon.GT.grid%nlon) ilon = ilon-grid%nlon
      ilat = FLOOR(1.0001 + fskal*(orglat-xlat) )

!      WRITE(*,*) xlon, xlat , ' --> ', ilon, ilat

    END SUBROUTINE indexcalc

!   *******************************************************************************
    SUBROUTINE areacalc(grid, vlat, area, dlat, dlon)
!   *******************************************************************************

      USE mo_constants,     ONLY: a, api

      TYPE(domain),INTENT(in)      :: grid

      REAL(dp), DIMENSION(grid%nlat), INTENT(in)  :: vlat  ! latitude coordinates (centre)
      REAL(dp), DIMENSION(grid%nlat), INTENT(out) :: area  ! Grid cell area in m2
      REAL(dp), INTENT(out), OPTIONAL                        :: dlat  ! Latitudinal distance in m
      REAL(dp), DIMENSION(grid%nlat), INTENT(out), OPTIONAL  :: dlon  ! Longitudinal distance in m

      REAL(dp) :: pifak
      REAL(dp) :: dist_lat
      REAL(dp), ALLOCATABLE :: dist_lon(:)

      INTEGER :: jb

      pifak = api/180._dp
      dist_lat = grid%resolution * pifak * a
      ALLOCATE(dist_lon(grid%nlat))

      DO jb=1, grid%nlat

        dist_lon(jb) = grid%resolution * pifak * cos(vlat(jb)*pifak) * a
        area(jb) = dist_lat * a *  &
            ABS( SIN((vlat(jb)+grid%resolution/2._dp) * pifak) - SIN((vlat(jb)-grid%resolution/2._dp)*pifak) )

        IF (jb.LE.5 .OR. jb.GE.grid%nlat-4 .OR. jb.EQ.grid%nlat/2  &
      	     .OR. jb.EQ.INT(grid%nlat/2) +1) THEN
           WRITE(*, '(I5,1X,F9.4,1X,F11.6,1X, F9.4,A4,F11.2,A4)')  &
                jb, vlat(jb), COS( vlat(jb)*pifak ), dist_lon(jb)/1000.,' km ', area(jb)*1.E-6, ' km2'
        ENDIF		  
      ENDDO

      IF (PRESENT(dlat)) dlat = dist_lat
      IF (PRESENT(dlon)) dlon(1:grid%nlat) = dist_lon(:)
      DEALLOCATE(dist_lon)

    END SUBROUTINE areacalc

!   *******************************************************************************
    SUBROUTINE read_grid_info(dnam, grid, clon, clat, iname, idim)
!   *******************************************************************************
      use netcdf
	  
      CHARACTER (LEN=*), INTENT(in)  :: dnam    ! File name with grid info
      TYPE(domain),INTENT(out)       :: grid    ! grid info structure
      CHARACTER (LEN=*), INTENT(in)  :: clon    ! Array name of longitudinal coordinates
      CHARACTER (LEN=*), INTENT(in)  :: clat    ! Array name of latitudinal coordinates
      INTEGER, INTENT(in), OPTIONAL  :: iname   ! Switch for of grid dimension names in Netcdf
	                                        !  1 = lon, lat (default)
	                                        !  2 = x, y
      INTEGER, INTENT(in), OPTIONAL  :: idim    ! Number of coordinate dimensions: 1 or 2 (default)
      LOGICAL                        :: LOG2    ! True for 2 dimensions

      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vlon
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vlat
!
!     *** NETCDF variables
      INTEGER :: ierr, ncid, lat_dimid, lon_dimid, varid(2), dimid
      CHARACTER (LEN=3)              :: cnl     ! Name of longitudinal dimension in netcdf
      CHARACTER (LEN=3)              :: cnb     ! Name of latitudinal dimension in netcdf
!
!     **********************************************
      cnl = 'lon'
      cnb = 'lat'
	  IF (PRESENT(iname)) THEN
	    IF (iname.EQ.2) THEN
		  cnl = 'x'
		  cnb = 'y'
        ENDIF
      ENDIF
!				 
!     ******** Open and read dimensions
      ierr = nf90_open(dnam, NF90_NOWRITE, ncid)

      ierr = nf90_inq_dimid(ncid,TRIM(cnl),lon_dimid)
      ierr = nf90_inquire_dimension(ncid,lon_dimid,len=grid%nlon)
      ierr = nf90_inq_dimid(ncid,TRIM(cnb),lat_dimid)
      ierr = nf90_inquire_dimension(ncid,lat_dimid,len=grid%nlat)
      WRITE(*,*) "Dimensions: NL =" , grid%nlon, "  NB = ", grid%nlat

      LOG2 = .TRUE.
      IF (PRESENT(idim)) THEN
        IF (idim.EQ.1) LOG2 = .FALSE.
      ENDIF
      IF (LOG2) THEN
        ALLOCATE(VLON(grid%nlon, grid%nlat))
        ALLOCATE(VLAT(grid%nlon, grid%nlat))
        WRITE(*,*) 'Grid has 2-dimensional coordinates'
      ELSE
        ALLOCATE(VLON(grid%nlon, 1))
        ALLOCATE(VLAT(grid%nlat, 1))
        WRITE(*,*) 'Grid has 1-dimensional coordinates'
      ENDIF
      ierr = nf90_inq_varid(ncid,TRIM(CLON),varid(1))
      ierr = nf90_get_var(ncid,varid(1),vlon)
      ierr = nf90_inq_varid(ncid,TRIM(CLAT),varid(2))
      ierr = nf90_get_var(ncid,varid(2),vlat)

      ierr = nf90_close(ncid)
      WRITE(*,*) TRIM(dnam)," is closed."
!
!     Calculate grid info
      WRITE(*,*) "Longitude range vlon:", MINVAL(vlon), MAXVAL(vlon)
      WRITE(*,*) "Latitude range vlat:", MINVAL(vlat), MAXVAL(vlat)
      IF (LOG2) THEN
        WRITE(*,*) 'Grid centre: ', vlon(grid%nlon/2, grid%nlat/2), vlat(grid%nlon/2, grid%nlat/2)
      ENDIF
      grid%resolution = (MAXVAL(vlon) - MINVAL(vlon)) / (grid%nlon-1)
      grid%origin_lon = MINVAL(vlon) - grid%resolution/2.
      grid%origin_lat = MAXVAL(vlat) + grid%resolution/2.
      WRITE(*,*) "Res. ", grid%resolution, '°  Origin ', grid%origin_lon, '°E, ', grid%origin_lat, '°N'
      grid%kshift_n = NINT((90. - grid%origin_lat) / grid%resolution)
      grid%kshift_s = NINT((90. + MINVAL(vlat) - grid%resolution/2.) / grid%resolution)
      WRITE(*,*) "N and S borders are shifted from the poles by ", grid%kshift_n, ' and ', grid%kshift_s

      grid%kshift_w = NINT((grid%origin_lon + 180. ) / grid%resolution)
      grid%kshift_e = NINT((180. - (MAXVAL(vlon) + grid%resolution/2.)) / grid%resolution )
      WRITE(*,*) "W and E borders are shifted from the dateline by ", grid%kshift_w, ' and ', grid%kshift_e
	  
      DEALLOCATE(VLON)
      DEALLOCATE(VLAT)

    END SUBROUTINE read_grid_info

END MODULE mo_grid

